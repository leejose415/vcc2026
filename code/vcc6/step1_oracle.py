"""Step 1, oracle diagnosis on the RPE1 dev 300 (eval RPE1dev, scored by the frozen harness).

  arms     A real cells (the reference's own perturbed cells), B zero (lfc 0, nothing called), C oracle (the
           replogle_rpe1 table's lfc and de_sig) -- B and C through stage1.prepare_generator_inputs and the frozen
           generator with knockdown_seed(name) (default salt), n_cells 400, library_size 20000 (vcc_2026_6's settings).
           Each prediction = the knockdown cells + the eval's context controls; scored with harness score.score_one.
  analyze  (separate module step1_analyze.py)

  ~/venv/bin/python code/vcc6/step1_oracle.py arms [--only A B C] [--workers 4]
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
import rpe1_eval as E                                                                      # noqa: E402  (patches the harness paths)
import gen_cells as GC                                                                     # noqa: E402

OUT = E.OUT / 'step1'
V5 = V / 'lib/effects_v5'
QUANT = V5 / 'replogle_rpe1_control_quantiles.npz'
BASE = V5 / 'replogle_rpe1_baseline.csv'


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


def targets():
    return [t for t in E.LIST.read_text().split() if t]


def table_rows(tg, genes):
    """The replogle_rpe1 table's lfc, de_sig, se for the knockdowns, on the panel in harness order."""
    import anndata as ad
    t = ad.read_h5ad(V5 / 'replogle_rpe1.h5ad')
    assert list(t.var_names.astype(str)) == list(genes), 'table var is not in panel order'
    row = {s: i for i, s in enumerate(t.obs['target_gene'].astype(str))}
    idx = np.array([row[s] for s in tg])
    lfc = np.asarray(t.layers['lfc'][idx], np.float64)
    se = np.asarray(t.layers['se'][idx], np.float64)
    sig = np.asarray(t.layers['de_sig'][idx].toarray(), bool)
    return lfc, sig, se, t.obs.iloc[idx].copy()


def control_mean(genes):
    base = pd.read_csv(BASE, index_col=0)
    assert list(base.index.astype(str)) == list(genes), 'baseline csv is not in panel order'
    return base['control_mean_cpm'].to_numpy(float)


def context_controls():
    import anndata as ad
    ctx = ad.read_h5ad(E.H.eval_dir(E.NAME) / 'context.h5ad')
    return ctx.X.tocsr()


def arm_real(path):
    import anndata as ad
    from common import write_h5ad
    real = ad.read_h5ad(E.H.eval_dir(E.NAME) / 'real.h5ad')
    lab = real.obs['target_gene'].astype(str).to_numpy()
    pert = np.flatnonzero(lab != E.H.CONTROL)
    ctx = context_controls()
    import scipy.sparse as sp
    X = sp.vstack([real.X[pert], ctx]).tocsr()
    labels = list(lab[pert]) + [E.H.CONTROL] * ctx.shape[0]
    write_h5ad(path, X, pd.DataFrame(dict(target=labels, target_gene=labels), index=[f'c{i}' for i in range(len(labels))]),
               real.var_names.astype(str).tolist())


def arms(only, workers):
    OUT.mkdir(parents=True, exist_ok=True)
    genes = E.H.panel_genes(); tg = targets()
    lfc, sig, _, _ = table_rows(tg, genes)
    measured = np.isfinite(lfc)
    m0 = control_mean(genes)
    for arm in only:
        path = OUT / f'arm{arm}_pred.h5ad'
        log(f'arm {arm}: writing {path}')
        if arm == 'A':
            arm_real(path)
        else:
            L = np.zeros_like(lfc) if arm == 'B' else np.where(measured, lfc, 0.0)
            S = np.zeros_like(sig) if arm == 'B' else sig
            GC.write_arm(path, tg, genes, L, S, measured.astype(np.int16), QUANT, m0, context_controls(), salt='',
                         workers=workers, log=log)
        s = E.score(f'step1_arm{arm}', path, threads=4, keep=True)
        log(f'arm {arm}: scaled {s["scaled"]}  raw {s["raw"]}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['arms'])
    p.add_argument('--only', nargs='+', default=['C', 'B', 'A'])
    p.add_argument('--workers', type=int, default=4)
    a = p.parse_args()
    arms(a.only, a.workers)

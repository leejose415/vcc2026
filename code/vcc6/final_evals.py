"""Final-target end-to-end evals (Jose 2026-10-27): the frozen harness (scripts/harness, not edited) on the knockdowns of
the D/E/F target list that a source line holds as QC'd cells, so the sources that vote are the ones that will vote on
D/E/F. As rpe1_eval.py: the harness's own build_cells / finish / score_one, its folders pointed at
results/stage1/final_evals/harness/.

  HCT116final  lib/qc_cells/orion_hct116_qc.h5ad   QC'd cells; knockdown = target_gene, kept by the Orion table QC
                                                    (results/qc/orion_hct116_table/knockdown_qc.csv); genes by symbol
  K562final    lib/qc_cells/replogle_k562_gw_qc.h5ad  kept cells of kept constructs of kept knockdowns (build.k562_kept,
                                                    as Eval B); genes by Ensembl ID through the panel table
Targets: results/target_lists/final_<line>.txt (the final list x the line's table, written by screen.py --set final),
restricted to knockdowns with kept cells. Per knockdown a seeded 400 of its cells (all when fewer); controls split into
context / scoring halves; thinning to the 20k regime -- all the harness's own rules.

  ~/venv/bin/python code/vcc6/final_evals.py build HCT116final|K562final|HEK293Tfinal [--threads 4] [--cells-only]
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

V = Path.home() / 'vcc'
OUT = V / 'results/stage1/final_evals'
sys.path.insert(0, str(V / 'scripts/harness'))
import common as H                                                                       # noqa: E402
H.LIB, H.RES = OUT / 'harness/lib', OUT / 'harness/results'
import build as HB                                                                       # noqa: E402
import score as HS                                                                       # noqa: E402
HS.RES = H.RES

EVALS = {'HCT116final': dict(cells=V / 'lib/qc_cells/orion_hct116_qc.h5ad', line='HCT116', kind='symbol',
                             qc=V / 'results/qc/orion_hct116_table/knockdown_qc.csv',
                             targets=V / 'results/target_lists/final_HCT116.txt'),
         'K562final': dict(cells=HB.K562, line='K562', kind='ensembl', qc=None,
                           targets=V / 'results/target_lists/final_K562.txt'),
         'HEK293Tfinal': dict(cells=V / 'lib/qc_cells/orion_hek293t_qc.h5ad', line='HEK293T', kind='symbol',
                              qc=V / 'results/qc/orion_hek293t_table/knockdown_qc.csv',
                              targets=V / 'results/target_lists/final_HEK293T.txt')}


def labels(name, obs):
    e = EVALS[name]
    ctrl = obs['is_control'].to_numpy(bool)
    if name == 'K562final':
        kept, label, ctrl = HB.k562_kept(obs)
        return np.where(ctrl, H.CONTROL, label)
    q = pd.read_csv(e['qc']); kept = set(q.loc[q.knockdown_keep.astype(bool), 'target_gene'].astype(str))
    tg = obs['target_gene'].astype(str).to_numpy()
    return np.where(ctrl, H.CONTROL, np.where(np.isin(tg, list(kept)), tg, ''))


def build(name, threads=4, cells_only=False, finish_only=False):
    import shutil
    e = EVALS[name]
    if finish_only:                     # scorer-space part 2 (2026-10-28): full DE bundle on the existing cells-only build
        d = H.eval_dir(name); keep = d / 'msepds_2026-10-28'
        if not keep.exists():           # a restart leaves the saved MSE/PDS set alone; finish clears its partial outputs
            keep.mkdir()
            for s in ('real_cache', 'baseline', 'bundle'):
                if (d / s).exists(): shutil.move(str(d / s), str(keep / s))
        HB.finish(name, H.read_json(d / 'split.json'), threads)
        return
    H.LIB.mkdir(parents=True, exist_ok=True)
    for f in ('panel.csv', 'generator_settings.json'):            # the harness's own files, read from its LIB
        shutil.copy(V / 'lib/harness' / f, H.LIB / f)
    with h5py.File(e['cells'], 'r') as fh:
        obs = read_elem(fh['obs']); var = read_elem(fh['var'])
    to_panel, genes = HB.panel_map(var, e['kind'])
    label = labels(name, obs)
    have = set(label[(label != '') & (label != H.CONTROL)])
    want = [t for t in e['targets'].read_text().split() if t]
    targets = sorted(t for t in want if t in have)
    print(f'{name}: {len(want)} listed final targets, {len(targets)} with kept cells; '
          f'dropped {sorted(set(want) - set(targets))[:20]}', flush=True)
    extra = dict(batch=obs['batch'].astype(str).to_numpy())
    if 'guide_id' in obs: extra['guide_id'] = obs['guide_id'].astype(str).to_numpy()
    split = HB.build_cells(name, e['cells'], label, targets, to_panel, genes, extra, name)
    split.update(line=e['line'], rule=f'final D/E/F targets held as QC\'d cells ({e["cells"].name})',
                 knockdown_list=str(e['targets']), listed=len(want))
    if cells_only:                      # scorer-space study (2026-10-28): real + context cells only, no DE bundle
        HB.save_json(H.eval_dir(name) / 'split.json', split)
        return
    HB.finish(name, split, threads)


def score(name, run, cells, threads=4, keep=False):
    """The harness's score_one on a final eval; waits for its build to finish."""
    import time
    while not (H.eval_dir(name) / 'build_invocations.json').exists():
        time.sleep(30)
    return HS.score_one(name, run, cells=cells, threads=threads, keep=keep)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['build'])
    p.add_argument('name', choices=list(EVALS))
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--cells-only', action='store_true', help='real.h5ad + context.h5ad + split.json only (no cache, '
                   'ceiling or bundle): the scorer-space study builds its own pseudobulk real side from them')
    p.add_argument('--finish-only', action='store_true', help='the harness finish (real cache, ceiling, DE bundle) on an '
                   'existing --cells-only build; its MSE/PDS-only cache, baseline and bundle move to msepds_2026-10-28/')
    a = p.parse_args()
    build(a.name, a.threads, a.cells_only, a.finish_only)

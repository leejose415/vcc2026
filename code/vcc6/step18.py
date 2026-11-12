"""Step 18 (Jose 2026-10-30): the D/E/F file for config 1' + K0. No upload.
Outputs in results/stage1/2026-10-30_step18/; write-up in the scorer-space report.md.

Config 1' + K0 with all 117 source tables (nothing held out): called genes (C+G's calls) 0.5 x C+G LFC; non-called
0.6 x c x Delta_pool, c = z_a^2 / (z_a^2 + 1), single-source c1 = 1, no source -> control (the generator input step as
in 15c, unchanged). Per context D, E, F: its own 18,400 controls, the 300 final targets, salt "".

  model    step15c.model (vcc_2026_7's data-loader + model code, config 1') for D, E, F in one process (one table load), with
           background None. Checks against the uploaded C+G file's inputs (results/stage1/2026-10-28_DEF/<x>/input.npz):
           identical calls and n_sources, called LFC = exactly 0.5 x C+G's.
  stats    before generation: per context, requested calls per knockdown (median, mean) and median |LFC| of called
           genes, next to the uploaded C+G file's -> tables/pre_generation.json
  (generate: run_line.py predict --lfc-npz, unchanged)
  final    per context: the submission file (controls stripped, context set, run_line.submission_file); merged on disk
           (anndata.experimental.concat_on_disk); shape / context checks; requested vs realised calls from the generator
           report -> tables/final.json. vcc prep runs after this, in step18.sh.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
for p in (V / 'code/vcc6', V / 'scripts', V / 'scripts/harness', V / 'code/cellgen'):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import step15c as C                                                                        # noqa: E402

OUT = V / 'results/stage1/2026-10-30_step18'
CFG = 'config1p_K0'
DEF = V / 'data_final/extracted'
CG = V / 'results/stage1/2026-10-28_DEF'
LINES = {x: dict(controls=DEF / f'context_{x}.h5ad', knockdowns=V / 'results/target_lists/final_DEF_300.txt',
                 background=None) for x in 'DEF'}


def log(msg):
    print(f'[{time.strftime("%H:%M:%S", time.gmtime())}] {msg}', flush=True)


def cg_checks(l, targets, B, A, n_src):
    z = np.load(CG / l / 'input.npz', allow_pickle=True)
    assert list(map(str, z['knockdowns'])) == targets
    sig = B['K0']
    return dict(cg_file_calls_identical=bool((z['de_sig'] == sig).all()),
                cg_file_n_sources_identical=bool((z['n_sources'] == n_src).all()),
                cg_file_lfc_max_abs_diff=float(np.abs(z['lfc'] - A['cg_lfc']).max()),
                called_lfc_vs_half_cg_max_abs_diff=float(np.abs(A['lfc'][sig] - 0.5 * z['lfc'][sig]).max()))


def model(a):
    a.config = CFG
    C.model(a, lines=LINES, out=OUT, extra_checks=cg_checks)


def req(z):
    """What the generator sees (stage1.prepare_generator_inputs): no-source genes unchanged and uncalled."""
    un = z['n_sources'] <= 0
    return np.where(un, 0.0, z['lfc'].astype(float)), z['de_sig'] & ~un


def stats(a):
    (OUT / 'tables').mkdir(parents=True, exist_ok=True); res = {}
    for x in 'DEF':
        row = {}
        for name, f in ((CFG, OUT / CFG / x / 'input.npz'), ('CG_uploaded', CG / x / 'input.npz')):
            lfc, sig = req(np.load(f, allow_pickle=True)); n = sig.sum(1)
            row[name] = dict(calls_median=float(np.median(n)), calls_mean=float(n.mean()), calls_min=int(n.min()),
                             calls_max=int(n.max()), abs_lfc_called_median=float(np.median(np.abs(lfc[sig]))),
                             abs_lfc_uncalled_median=float(np.median(np.abs(lfc[~sig & (lfc != 0)]))))
        res[x] = row
        log(f'{x}: {json.dumps(row)}')
    (OUT / 'tables/pre_generation.json').write_text(json.dumps(res, indent=1) + '\n')


def final(a):
    import anndata as ad
    import anndata.experimental as AE
    import run_line as RL
    fin = OUT / 'final'; fin.mkdir(parents=True, exist_ok=True); res = {}
    for x in 'DEF':
        sub = fin / f'submission_{x}.h5ad'
        if not sub.exists():
            kept, dropped = RL.submission_file(OUT / CFG / x / 'e2e/salt_none/prediction.h5ad', sub, context=x)
            log(f'{x}: {kept} knockdown cells kept, {dropped} control cells removed')
        g = pd.read_csv(OUT / CFG / x / 'e2e/salt_none/prediction_generator.csv')
        res[x] = dict(requested_median=float(g.sig_requested.median()), requested_mean=float(g.sig_requested.mean()),
                      realised_median=float(g.sig_realized.median()), realised_mean=float(g.sig_realized.mean()),
                      recall_median=float(g.sig_recall.median()), precision_median=float(g.sig_precision.median()))
    merged = fin / 'prediction_DEF.h5ad'
    if not merged.exists():
        AE.concat_on_disk([str(fin / f'submission_{x}.h5ad') for x in 'DEF'], str(merged), join='inner', merge='same',
                          keys=list('DEF'), index_unique='-')
    m = ad.read_h5ad(merged, backed='r')
    genes = pd.read_csv(DEF / 'gene_names.csv').iloc[:, 0].astype(str).tolist()
    tg = [t for t in (V / 'results/target_lists/final_DEF_300.txt').read_text().split() if t]
    obs = m.obs
    out = dict(shape=list(m.shape), genes_match_final_list=list(m.var_names.astype(str)) == genes,
               obs_names_unique=bool(m.obs_names.is_unique), contexts=obs['context'].value_counts().to_dict(),
               per_context={x: dict(cells=int((obs.context == x).sum()),
                                    knockdowns=int(obs.loc[obs.context == x, 'target_gene'].nunique()),
                                    cells_per_knockdown=sorted(set(obs.loc[obs.context == x, 'target_gene'].value_counts().tolist())),
                                    targets_match_list=sorted(obs.loc[obs.context == x, 'target_gene'].astype(str).unique()) == sorted(tg),
                                    controls=int((obs.loc[obs.context == x, 'target_gene'] == 'non-targeting').sum()))
                            for x in 'DEF'},
               calls=res, size_bytes=merged.stat().st_size)
    m.file.close()
    assert out['shape'] == [360000, len(genes)] and out['genes_match_final_list'] and out['obs_names_unique']
    assert all(v['cells'] == 120000 and v['knockdowns'] == 300 and v['targets_match_list'] and v['controls'] == 0
               for v in out['per_context'].values()), out['per_context']
    (OUT / 'tables/final.json').write_text(json.dumps(out, indent=1, default=str) + '\n')
    log(json.dumps(out, default=str))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['model', 'stats', 'final'])
    a = p.parse_args()
    {'model': model, 'stats': stats, 'final': final}[a.cmd](a)

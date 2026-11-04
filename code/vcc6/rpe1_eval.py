"""RPE1 dev eval for stage 1, built and scored by the frozen harness (scripts/harness), never edited.

The harness hard-codes its folders (common.LIB / score.RES). This module points them at
results/stage1/2026-10-27_rpe1/harness/ and then calls the harness's own build_cells, finish and score_one,
so the reference is built exactly as Eval B's is (per knockdown a seeded 400 of its kept cells, all when fewer;
controls split into context / scoring halves; thinning to the 20k regime; panel reindex).

  ~/venv/bin/python code/vcc6/rpe1_eval.py sample          # results/target_lists/rpe1_dev_300.txt
  ~/venv/bin/python code/vcc6/rpe1_eval.py build           # eval RPE1dev (reference, cache, ceiling, bundle)
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

V = Path.home() / 'vcc'
OUT = V / 'results/stage1/2026-10-27_rpe1'
sys.path.insert(0, str(V / 'scripts/harness'))
import common as H                                                                       # noqa: E402
H.LIB, H.RES = OUT / 'harness/lib', OUT / 'harness/results'
import build as HB                                                                       # noqa: E402
import score as HS                                                                       # noqa: E402
HS.RES = H.RES

NAME = 'RPE1dev'
CELLS = V / 'lib/qc_cells/replogle_rpe1_qc.h5ad'
QC = V / 'results/qc/replogle_rpe1'
TABLE = V / 'lib/effects_v5/replogle_rpe1.h5ad'
LIST = V / 'results/target_lists/rpe1_dev_300.txt'
N_DEV, SEED = 300, 2026


def rpe1_kept(obs):
    """build.k562_kept with RPE1's QC: cell -> kept guide -> kept knockdown symbol ('' = not kept)."""
    guides = pd.read_csv(QC / 'guide_qc.csv'); kds = pd.read_csv(QC / 'knockdown_qc.csv')
    kept = kds[kds.knockdown_keep.astype(bool)]
    sym_of = dict(zip(kept.target_gene_id.astype(str), kept.target_gene.astype(str)))
    g = guides[guides.keep.astype(bool)]
    guide_kd = {gid: sym_of.get(t, '') for gid, t in zip(g.guide_id.astype(str), g.target_gene_id.astype(str))}
    ctrl = obs['is_control'].to_numpy(bool)
    label = np.array(['' if c else guide_kd.get(x, '') for c, x in zip(ctrl, obs['guide_id'].astype(str))], dtype=object)
    return set(kept.target_gene.astype(str)), label, ctrl


def cells_obs():
    with h5py.File(CELLS, 'r') as fh:
        return read_elem(fh['obs']), read_elem(fh['var'])


def sample():
    obs, _ = cells_obs()
    _, label, _ = rpe1_kept(obs)
    counts = pd.Series(label[label != '']).value_counts()
    with h5py.File(TABLE, 'r') as fh:
        tab = read_elem(fh['obs'])
    in_table = set(tab['target_gene'].astype(str))
    pool = sorted(in_table & set(counts.index))
    picked = sorted(np.random.default_rng(SEED).choice(np.array(pool, dtype=object), N_DEV, replace=False).tolist())
    LIST.parent.mkdir(parents=True, exist_ok=True)
    LIST.write_text('\n'.join(picked) + '\n')
    print(f'{len(in_table)} table knockdowns, {len(counts)} with kept cells, {len(pool)} in both; '
          f'{len(picked)} drawn (default_rng({SEED}).choice of the sorted pool) -> {LIST}')
    return picked, counts


def build(threads=8):
    import shutil
    H.LIB.mkdir(parents=True, exist_ok=True)
    for f in ('panel.csv', 'generator_settings.json'):            # the harness's own files, read from its LIB
        shutil.copy(V / 'lib/harness' / f, H.LIB / f)
    obs, var = cells_obs()
    to_panel, genes = HB.panel_map(var, 'ensembl')
    _, label, ctrl = rpe1_kept(obs)
    label = np.where(ctrl, H.CONTROL, label)
    targets = [t for t in LIST.read_text().split() if t]
    split = HB.build_cells(NAME, CELLS, label, targets, to_panel, genes,
                           dict(guide_id=obs['guide_id'].astype(str).to_numpy(), batch=obs['batch'].astype(str).to_numpy()), NAME)
    split.update(line='RPE1', rule="QC'd cells of kept constructs of kept knockdowns (guide_qc.csv -> knockdown_qc.csv)",
                 knockdown_list=str(LIST))
    HB.finish(NAME, split, threads)


def score(run, cells, threads=8, keep=False):
    """The harness's score_one on eval RPE1dev (a cell prediction: knockdown cells + control cells). Waits for the
    eval's build to finish (finish() writes build_invocations.json last)."""
    import time
    while not (H.eval_dir(NAME) / 'build_invocations.json').exists():
        time.sleep(30)
    return HS.score_one(NAME, run, cells=cells, threads=threads, keep=keep)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['sample', 'build'])
    p.add_argument('--threads', type=int, default=8)
    a = p.parse_args()
    sample() if a.cmd == 'sample' else build(a.threads)

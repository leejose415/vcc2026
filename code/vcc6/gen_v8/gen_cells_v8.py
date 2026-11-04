"""Step 22 Part E1: a COPY of code/vcc6/gen_cells.py (not edited) that drives generator v8 (gen_v8/cell_generator_v8.py)
with a per-(knockdown, gene) rank input; everything else as gen_cells.

A whole arm through the frozen generator (code/cellgen/cell_generator.py, not edited), in parallel.

Per knockdown: stage1.prepare_generator_inputs(lfc, de_sig, n_sources, p_sign) -> generate_cells(..., seed =
stage1.knockdown_seed(name, salt)). The prediction file is every knockdown's cells plus the given control cells
(label 'non-targeting'), the layout the harness's trivial arms submit.
"""
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

V = Path.home() / 'vcc'
for p in (V / 'code/vcc6/gen_v8', V / 'code/vcc6'):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import cell_generator_v8 as C                                                              # noqa: E402
from stage1 import knockdown_seed, prepare_generator_inputs                                # noqa: E402

CONTROL = 'non-targeting'
_W = {}                                    # shared with forked workers (copy-on-write)


def _one(i):
    w = _W
    name = w['targets'][i]
    packed = prepare_generator_inputs(w['lfc'][i], w['de_sig'][i], w['n_sources'][i], w['p_sign'][i])
    seed = knockdown_seed(name, w['salt'])
    a = C.generate_cells(packed['lfc'], packed['de_sig'], w['cq'], w['control_mean'], n_cells=w['n_cells'],
                         library_size=w['library_size'], seed=seed, target=name, genes=w['genes'], rank=w['rank'][i])
    r = a.uns['cellgen']
    rep = dict(knockdown=name, seed=seed, n_called_in=int(packed['n_called']), n_unmeasured=int(packed['n_unmeasured']),
               sig_requested=r['sig_requested'], sig_realized=r['sig_realized'], sig_recall=r['sig_recall'],
               sig_precision=r['sig_precision'], composition=r['composition'], library_median=r['library_median'],
               mean_count_max_abs_error=r['mean_count_max_abs_error'])
    return i, sp.csr_matrix(a.X), rep


def write_arm(path, targets, genes, lfc, de_sig, n_sources, quantiles, control_mean, controls_X, salt='',
              p_sign=None, n_cells=400, library_size=20000, workers=4, log=print, rank=None):
    """Generate every knockdown and write the prediction h5ad (knockdowns, then controls). Returns the per-knockdown
    generator report as a DataFrame (also written next to `path` as <stem>_generator.csv)."""
    sys.path.insert(0, str(V / 'scripts/harness'))
    from common import write_h5ad
    K, G = len(targets), len(genes)
    _W.update(targets=list(map(str, targets)), genes=list(map(str, genes)), lfc=np.asarray(lfc, np.float32),
              de_sig=np.asarray(de_sig, bool), n_sources=np.asarray(n_sources),
              p_sign=np.full((K, G), 0.5, np.float32) if p_sign is None else np.asarray(p_sign, np.float32),
              cq=quantiles if isinstance(quantiles, C.ControlQuantiles) else C.ControlQuantiles(quantiles),
              control_mean=np.asarray(control_mean, float), salt=str(salt), n_cells=int(n_cells),
              library_size=float(library_size), rank=np.asarray(rank, np.float32))
    blocks, reps = [None] * K, [None] * K
    t0 = time.time()
    with mp.get_context('fork').Pool(workers) as pool:
        for done, (i, X, rep) in enumerate(pool.imap_unordered(_one, range(K)), 1):
            blocks[i], reps[i] = X, rep
            if done % 25 == 0 or done == K:
                log(f'  generated {done}/{K} knockdowns ({time.time() - t0:.0f}s)')
    _W.clear()
    labels = [t for t in targets for _ in range(n_cells)] + [CONTROL] * controls_X.shape[0]
    X = sp.vstack(blocks + [sp.csr_matrix(controls_X)]).tocsr()
    obs = pd.DataFrame(dict(target=labels, target_gene=labels), index=[f'c{i}' for i in range(len(labels))])
    write_h5ad(path, X, obs, list(map(str, genes)))
    rep = pd.DataFrame(reps)
    rep.to_csv(Path(path).with_name(Path(path).stem + '_generator.csv'), index=False)
    return rep

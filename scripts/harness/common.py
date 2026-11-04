"""Shared pieces of the local eval harness (see lib/harness/README.md).

Evals (lib/harness/eval<name>/):
  A         line 2 = H1, the whole line held out; reference from the RAW unfiltered lib/h1/H1_full_2025.h5ad
  A_shared  Eval A restricted to the H1 knockdowns that K562 also keeps (same cells, its own cache and bundle)
  B         line 1 = K562, 100 kept knockdowns held out (from those shared with H1; filled at random if fewer)
Regime (both): per knockdown a seeded 400 of its real cells (all when fewer); controls capped at 2 x 18,400 and split
into two seeded halves, the context (what a model / the generator conditions on) and the scoring controls (in the
reference); every cell binomially thinned by one factor f = min(1, 20,000 / median control panel library), seeded per
raw row; genes reindexed to the 18,533-gene panel (data/gene_names.csv order; unmeasured = 0 counts).
Scoring: official cell-eval2 (vcc2026 preset, pdex backend) through vcc3_generate.OfficialEvaluator; the cache, ceiling
and bundle (generic baseline = 0, replicate anchor = 1) are built once per eval. Effect tables go through the
existing generator (vcc3_generate: calibrate / pick_cells / write_arm) with fixed settings (generator_settings.json).
"""
import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

V = Path.home() / 'vcc'
LIB, RES = V / 'lib/harness', V / 'results/harness'
sys.path.insert(0, str(V / 'code/vcc_2026_3c'))
import vcc3_generate as G                                                                   # noqa: E402
from vcc2_common import seed_for, expression_gate, save_json                                 # noqa: E402

CONTROL, SEED, DEPTH, N_CELLS, N_CTRL = 'non-targeting', 2026, 20000, 400, 18400
EVALS = ('A', 'A_shared', 'B')
CHUNK = 50000


def log(msg):
    import time
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


def eval_dir(name):
    return LIB / f'eval{name}'


def panel_genes():
    return pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).tolist()


def evaluator(name, out=None, threads=8):
    """OfficialEvaluator on one eval's reference; `out` is where a scored arm goes (cache and bundle stay in the eval)."""
    exe = Path(sys.executable).parent / 'cell-eval2'                # the venv's, whatever PATH the job runs with
    a = argparse.Namespace(cell_eval=str(exe) if exe.exists() else 'cell-eval2', eval_profile='vcc2026',
                           num_threads=threads, anchor_splits=5)
    d = eval_dir(name)
    return G.OfficialEvaluator(a, d / 'real.h5ad', out or d, reuse=None if out is None else d)


def thin(counts, row, f):
    """The harness's one depth draw for a cell: binomial thinning by f, seeded by the cell's row in its source file."""
    if f >= 1.0:
        return counts
    return np.random.default_rng([SEED, int(row)]).binomial(counts, f).astype(counts.dtype)


def read_rows(path, rows, to_panel, n_panel):
    """Rows (sorted) of an h5ad's CSR X, columns mapped onto the panel (to_panel[src] = panel col or -1), raw counts.
    One streaming pass in CHUNK-row blocks into buffers preallocated at the rows' stored size (never two copies).
    `rows` must be sorted and unique. Returns a CSR (int32) with the rows in that order."""
    rows = np.asarray(rows)
    with h5py.File(path, 'r') as fh:
        X = fh['X']; ip = X['indptr'][:]; n_obs = len(ip) - 1
        want = np.zeros(n_obs, bool); want[rows] = True
        ub = int(np.diff(ip)[rows].sum())
        data = np.empty(ub, np.int32); idx = np.empty(ub, np.int32); ptr = np.zeros(len(rows) + 1, np.int64)
        w = 0; k_row = 0
        for s in range(0, n_obs, CHUNK):
            e = min(s + CHUNK, n_obs); sel = np.flatnonzero(want[s:e])
            if not len(sel): continue
            a, b = int(ip[s]), int(ip[e]); D = X['data'][a:b]; I = X['indices'][a:b]
            for r in sel + s:
                lo, hi = int(ip[r] - a), int(ip[r + 1] - a)
                c = to_panel[I[lo:hi]]; k = c >= 0; m = int(k.sum())
                data[w:w + m] = D[lo:hi][k]; idx[w:w + m] = c[k]; w += m
                k_row += 1; ptr[k_row] = w
            if (s // CHUNK) % 2 == 0: log(f'  read {e:,} / {n_obs:,} rows')
    M = sp.csr_matrix((data[:w], idx[:w], ptr), shape=(len(rows), n_panel))
    M.sort_indices()
    return M


def write_h5ad(path, X, obs, genes):
    """gzip-compressed h5ad: integer counts CSR, obs with both `target` and `target_gene`, var = the panel symbols."""
    import anndata as ad
    obs = obs.copy()
    for c in obs.columns:
        if pd.api.types.is_string_dtype(obs[c]) or obs[c].dtype == object: obs[c] = obs[c].astype(str).astype('category')
    a = ad.AnnData(X=X.astype(np.int32), obs=obs, var=pd.DataFrame(index=pd.Index(genes, name='gene_name')))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(path, compression='gzip')


def targets_of(name):
    """The eval's held-out knockdowns and their real cell counts (from its reference)."""
    return G.reference_targets(eval_dir(name) / 'real.h5ad')


def generator_line(name):
    """The generator's input for one eval: its context controls as a vcc3_generate.Line, with the gate and the reference
    size taken from the scoring controls, exactly as vcc3_generate.h1_controls sets them."""
    import anndata as ad
    from vcc2_common import moments, alignment
    genes = panel_genes(); d = eval_dir(name)
    ctx = ad.read_h5ad(d / 'context.h5ad')
    guide = ctx.obs['guide_id'].astype(str).values if 'guide_id' in ctx.obs else None
    line = G.Line(G.aligned_X(ctx, genes), guide, genes); del ctx
    real = ad.read_h5ad(d / 'real.h5ad', backed='r')
    try:
        score_rows = np.flatnonzero(real.obs['target_gene'].astype(str).values == CONTROL)
        line.gate = expression_gate(moments(real, score_rows, alignment(real.var_names.astype(str), genes))['mean_cpm'])
        line.n_reference = len(score_rows)
    finally:
        real.file.close()
    return line


def read_json(path):
    return json.loads(Path(path).read_text())

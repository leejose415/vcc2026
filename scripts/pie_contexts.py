"""Control-only PIE contexts and queries for the PIE-clean lines (K562, RPE1, HepG2, Jurkat) and challenge lines A-C.

Per line: lib/pie/contexts/<line>.h5ad with the line's control cells only, X = log1p(CP10k) (per-cell library = the cell's
total UMI where the file records it, else its row sum), obs context / target_gene = 'non-targeting', and var = the panel
genes the line measured that PIE's gene_text has under the same symbol (PIE looks genes up by symbol, without aliases).
Contexts are named as PIE's own Replogle contexts (k562, rpe1, hepg2, jurkat) so PIE's context_text applies; A-C are
lineA / lineB / lineC (no context_text: masked). Query lib/pie/queries/<line>.json: {"vcc_<line>.<context>": [targets]}
with every effects_v4 knockdown of the line (K562 ... Jurkat) or our 300 targets (A-C).
PIE is leaked for HCT116, HEK293T and H1: they are refused here.
Usage: python scripts/pie_contexts.py --lines A K562 ...  [--max-controls 20000]
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

V = Path.home() / 'vcc'
OUT = V / 'lib/pie'
LINES = {'K562': ('replogle_k562_gw', 'k562'), 'RPE1': ('replogle_rpe1', 'rpe1'), 'HepG2': ('nadig_hepg2', 'hepg2'),
         'Jurkat': ('nadig_jurkat', 'jurkat'), 'A': (None, 'lineA'), 'B': (None, 'lineB'), 'C': (None, 'lineC')}
LEAKED = {'HCT116', 'HEK293T', 'H1'}
CHUNK = 20_000


def gene_text_symbols():
    cov = pd.read_csv(V / 'lib/pie_features/coverage.csv', index_col=0)
    return set(cov.index[cov['gene_text'] == 'symbol'])


def build(line, max_controls, seed=2026):
    if line in LEAKED: raise SystemExit(f'{line}: PIE is leaked for this line; no PIE context is made')
    source, context = LINES[line]
    panel = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy()
    gt = gene_text_symbols()
    if source is None:
        path = V / f'data/context_{line}.h5ad'
        with h5py.File(path, 'r') as f:
            obs = read_elem(f['obs']); var_sym = read_elem(f['var']).index.astype(str).to_numpy()
        ctrl = np.flatnonzero(obs['target_gene'].astype(str).to_numpy() == 'non-targeting'); lib_col = None
        measured = set(panel)
        targets = sorted(pd.read_csv(V / 'data/pert_counts.csv')['target_gene'].astype(str))
    else:
        path = V / f'lib/qc_cells/{source}_qc.h5ad'
        with h5py.File(path, 'r') as f:
            obs = read_elem(f['obs']); var = read_elem(f['var'])
        with h5py.File(V / f'lib/effects_v4/{source}.h5ad', 'r') as f:
            pvar = read_elem(f['var']); kobs = read_elem(f['obs'])
        id2sym = {g: s for s, g in zip(pvar.index.astype(str), pvar['gene_id'].astype(str)) if g}
        src_ids = var.index.astype(str).to_numpy() if 'gene_name' in var else var['gene_id'].astype(str).to_numpy()
        var_sym = np.array([id2sym.get(g, '') for g in src_ids])
        ctrl = np.flatnonzero(obs['is_control'].to_numpy(bool))
        lib_col = next((c for c in ('n_umi', 'UMI_count') if c in obs), None)
        base = pd.read_csv(V / f'lib/effects_v4/{source}_baseline.csv', index_col=0)['control_mean_cpm']
        measured = set(base.index[base.notna()])
        targets = sorted(set(kobs['target_gene'].astype(str)))
    if len(ctrl) > max_controls:
        ctrl = np.sort(np.random.default_rng(seed).choice(ctrl, max_controls, replace=False))
    keep_cols = np.array([i for i, s in enumerate(var_sym) if s and s in gt and s in measured])
    _, first = np.unique(var_sym[keep_cols], return_index=True); keep_cols = np.sort(keep_cols[first])
    lib = obs[lib_col].to_numpy(float)[ctrl] if lib_col else None
    blocks = []
    with h5py.File(path, 'r') as f:
        X = f['X']; ip = X['indptr']; n_vars = int(X.attrs['shape'][1])
        for s in range(0, len(ctrl), CHUNK):
            rows = ctrl[s:s + CHUNK]
            parts = []
            for r in rows:                       # control rows are scattered; read each row's slice of the CSR arrays
                a, b = int(ip[r]), int(ip[r + 1])
                parts.append(sp.csr_matrix((X['data'][a:b], X['indices'][a:b], [0, b - a]), shape=(1, n_vars)))
            blocks.append(sp.vstack(parts).tocsr())
    M = sp.vstack(blocks).tocsr().astype(np.float64)
    total = lib if lib is not None else np.asarray(M.sum(1)).ravel()
    M = sp.diags(1e4 / np.maximum(total, 1)) @ M[:, keep_cols]
    M.data = np.log1p(M.data)
    a = ad.AnnData(X=M.astype(np.float32).tocsr(),
                   obs=pd.DataFrame(dict(context=context, target_gene='non-targeting'), index=[f'c{i}' for i in range(len(ctrl))]),
                   var=pd.DataFrame(index=pd.Index(var_sym[keep_cols], name='gene')))
    (OUT / 'contexts').mkdir(parents=True, exist_ok=True); (OUT / 'queries').mkdir(parents=True, exist_ok=True)
    a.write_h5ad(OUT / f'contexts/{line}.h5ad')
    q = {f'vcc_{line.lower()}.{context}': targets}
    (OUT / f'queries/{line}.json').write_text(json.dumps(q))
    info = dict(line=line, context=context, controls=len(ctrl), genes=len(keep_cols), library='total UMI' if lib_col else 'row sum',
                queries=len(targets))
    print(json.dumps(info), flush=True)
    return info


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--lines', nargs='+', required=True)
    ap.add_argument('--max-controls', type=int, default=20000); a = ap.parse_args()
    for line in a.lines: build(line, a.max_controls)


if __name__ == '__main__':
    main()

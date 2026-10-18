"""Downstream-change sanity check for the on-target threshold: do knockdowns that pass it change other genes more?

Per knockdown (tables) or construct (cells), against the controls, every measured readout gene except the target:
  log2FC = log2((mean CPM + 0.5) / (control mean CPM + 0.5))
  SE     = the effect tables' SE: delta method, sqrt(relvar_kd + relvar_ctrl) / ln 2, relvar = max(var of the mean,
           (mean + u) * u) / (mean + 0.5)^2, u = the CPM one count adds to the group mean; floored at 0.05
  counts: all = |log2FC| > 0.5;  noise_aware = |log2FC| > 0.5 and |log2FC| / SE > 3
Groups (measurable knockdowns / constructs only): dropped (enough cells, failed the on-target rule), band_-1_to_-0.5
(kept, -1 < on-target lfc <= -0.5), le_-1 (kept, on-target lfc <= -1); dropped_few_cells reported apart.

Modes:
  --table HCT116|HEK293T   Orion per-knockdown sums (lib/orion/<line>_pseudobulk.h5ad) + results/qc/orion_<line>_table/knockdown_qc.csv
  --cells <dataset>        a cell-QC'd file + results/qc/<source>/guide_qc.csv (constructs; every measurable one)
Writes downstream_counts.csv and downstream_summary.json next to the QC table.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

sys.path.insert(0, str(Path(__file__).parent))
from qc_replogle_umi import SOURCES                                                        # noqa: E402

V = Path.home() / 'vcc'
PSEUDO, SE_FLOOR, LFC_MIN, Z_MIN, MAX_ONTARGET_LFC = 0.5, 0.05, 0.5, 3.0, -0.5


def relvar(mean, var_of_mean, unit):
    return np.maximum(var_of_mean, (np.maximum(mean, 0) + unit) * unit) / (mean + PSEUDO) ** 2


def counts(m, vm, u, m0, vm0, u0, exclude, use_genes):
    """Per row: genes (in use_genes, not the row's own target) with |lfc| > 0.5, and with |lfc| > 0.5 and |lfc| / SE > 3."""
    lfc = np.log2((m + PSEUDO) / (m0 + PSEUDO))
    se = np.maximum(np.sqrt(relvar(m, vm, u) + relvar(m0, vm0, u0)) / np.log(2), SE_FLOOR)
    big = (np.abs(lfc) > LFC_MIN) & use_genes[None, :]
    rows = np.flatnonzero(exclude >= 0); big[rows, exclude[rows]] = False
    return big.sum(1), (big & (np.abs(lfc) / se > Z_MIN)).sum(1)


def summarise(df, out, label):
    df.to_csv(out / 'downstream_counts.csv', index=False)
    res = {}
    for g, d in df.groupby('group'):
        res[g] = dict(n=len(d), median_all=float(d['n_all'].median()), median_noise_aware=float(d['n_noise_aware'].median()),
                      mean_noise_aware=float(d['n_noise_aware'].mean()))
    summary = dict(source=label, rule=f'|log2FC| > {LFC_MIN} and |log2FC| / SE > {Z_MIN} (SE as the effect tables, floor {SE_FLOOR})', groups=res)
    (out / 'downstream_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    order = ['dropped', 'band_-1_to_-0.5', 'le_-1', 'dropped_few_cells']
    print(f'\n{label}: median downstream genes per group (all |lfc| > 0.5 / noise-aware)')
    print(pd.DataFrame(res).T.reindex([o for o in order if o in res]).to_string())


def group_of(measurable, keep, lfc, enough):
    if not measurable: return None
    if not enough: return 'dropped_few_cells'
    if not keep: return 'dropped'
    return 'le_-1' if lfc <= -1 else 'band_-1_to_-0.5'


def table_mode(line, chunk=256):
    out = V / f'results/qc/orion_{line.lower()}_table'
    qc = pd.read_csv(out / 'knockdown_qc.csv')
    qc['group'] = [group_of(m, k, l, e) for m, k, l, e in zip(qc.measurable, qc.knockdown_keep, qc.ontarget_lfc, qc.enough_cells)]
    qc = qc[qc.group.notna()]
    with h5py.File(V / f'lib/orion/{line}_pseudobulk.h5ad', 'r') as f:
        obs = read_elem(f['obs']); var = read_elem(f['var'])
        names = obs.index.astype(str).to_numpy(); pos = {t: i for i, t in enumerate(names)}; n = obs['n_cells'].to_numpy(float)
        genes = var.index.astype(str).to_numpy(); gpos = {g: i for i, g in enumerate(genes)}; use = var['measured'].to_numpy(bool)
        L = {k: f['layers'][k] for k in ('sum_cpm', 'sum_cpm2', 'sum_counts')}
        def moments(rows):
            rows = np.sort(rows); nn = n[rows][:, None]
            m = L['sum_cpm'][rows].astype(float) / nn
            vm = np.maximum(L['sum_cpm2'][rows].astype(float) / nn - m ** 2, 0) / nn
            u = 1e6 / np.maximum(L['sum_counts'][rows].astype(float).sum(1, keepdims=True), 1)
            return rows, m, vm, u
        _, m0, vm0, u0 = moments(np.array([pos['non-targeting']]))
        rows_all = np.array([pos[t] for t in qc.target_gene]); res = {}
        t0 = time.time()
        for s in range(0, len(rows_all), chunk):
            rows, m, vm, u = moments(rows_all[s:s + chunk])
            excl = np.array([gpos.get(names[r], -1) for r in rows])
            a, b = counts(m, vm, u, m0, vm0, u0, excl, use)
            for r, x, y in zip(rows, a, b): res[names[r]] = (int(x), int(y))
            print(f'  {min(s + chunk, len(rows_all)):,} / {len(rows_all):,} knockdowns ({time.time() - t0:.0f}s)', flush=True)
    qc['n_all'] = qc.target_gene.map(lambda t: res[t][0]); qc['n_noise_aware'] = qc.target_gene.map(lambda t: res[t][1])
    summarise(qc[['target_gene', 'group', 'n_cells', 'ontarget_lfc', 'n_all', 'n_noise_aware']], out, f'X-Atlas/Orion {line} (all knockdowns, from sums)')


def cells_mode(dataset, chunk=50000):
    _, out, path, title = SOURCES[dataset]
    guides = pd.read_csv(out / 'guide_qc.csv')
    guides['group'] = [group_of(m, k, l, True) for m, k, l in zip(guides.measurable, guides.keep, guides.ontarget_lfc)]
    guides = guides[guides.group.notna()].reset_index(drop=True)
    gidx = {g: i for i, g in enumerate(guides.guide_id.astype(str))}
    t0 = time.time()
    with h5py.File(path, 'r') as f:
        obs = read_elem(f['obs']); var = read_elem(f['var']); X = f['X']; n_obs = len(obs); n_vars = int(X.attrs['shape'][1])
        ids = var.index.astype(str).to_numpy() if 'gene_name' in var else var['gene_id'].astype(str).to_numpy()
        col = {g: i for i, g in enumerate(ids)}
        ctrl = obs['is_control'].to_numpy(bool); lib = obs['n_umi'].to_numpy(float)
        grp = np.where(ctrl, len(guides), [gidx.get(g, -1) for g in obs['guide_id'].astype(str)])
        k = len(guides) + 1
        s1 = np.zeros((k, n_vars)); s2 = np.zeros((k, n_vars)); nn = np.zeros(k); depth = np.zeros(k)
        ip = X['indptr']
        for s in range(0, n_obs, chunk):
            e = min(s + chunk, n_obs); g = grp[s:e]; r = np.flatnonzero(g >= 0)
            if len(r):
                p = ip[s:e + 1]; a, b = int(p[0]), int(p[-1])
                M = sp.csr_matrix((X['data'][a:b], X['indices'][a:b], p - a), shape=(e - s, n_vars))[r]
                M = sp.diags(1e6 / lib[s:e][r]) @ M.astype(np.float64)
                G = sp.csr_matrix((np.ones(len(r)), (g[r], np.arange(len(r)))), shape=(k, len(r)))
                s1 += (G @ M).toarray(); s2 += (G @ M.multiply(M)).toarray()
                nn += np.bincount(g[r], minlength=k); depth += np.bincount(g[r], weights=lib[s:e][r], minlength=k)
            print(f'  {e:,} / {n_obs:,} cells ({time.time() - t0:.0f}s)', flush=True)
    m = s1 / np.maximum(nn, 1)[:, None]
    vm = np.maximum(s2 / np.maximum(nn, 1)[:, None] - m ** 2, 0) * nn[:, None] / np.maximum(nn - 1, 1)[:, None] / np.maximum(nn, 1)[:, None]
    u = 1e6 / np.maximum(depth, 1)[:, None]
    excl = np.array([col.get(t, -1) for t in guides.target_gene_id.astype(str)])
    a, b = counts(m[:-1], vm[:-1], u[:-1], m[-1:], vm[-1:], u[-1:], excl, np.ones(n_vars, bool))
    guides['n_all'] = a; guides['n_noise_aware'] = b
    summarise(guides[['guide_id', 'target_gene', 'group', 'n_cells', 'ontarget_lfc', 'n_all', 'n_noise_aware']], out, f'{title} (every measurable construct, from cells)')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--table', choices=['HCT116', 'HEK293T']); g.add_argument('--cells', choices=list(SOURCES))
    a = p.parse_args()
    table_mode(a.table) if a.table else cells_mode(a.cells)

"""Perturbation-level QC on a cell-QC'd file (lib/qc_cells/<source>_qc.h5ad). Read only; writes tables and a plot.

Parameters: MIN_CELLS = 30, MIN_CTRL_CPM = 5, MAX_ONTARGET_LFC = -0.5, P_MAX = 0.05. Every readout gene is kept (no CPM mask).

Per cell: CPM = counts / n_umi * 1e6, computed on the fly. The control cells' counts are loaded once; every other cell is
streamed in row chunks and keeps only the CPM of its own target gene.

1. Per construct (guide_id):
     target = target_gene_id; not in .var -> measurable = False
     ctrl_cpm = mean CPM of the target in control cells; measurable = ctrl_cpm >= MIN_CTRL_CPM
       (decides only whether the knockdown can be verified; removes no genes)
     ontarget_lfc = log2((mean CPM in the construct's cells + 0.5) / (ctrl_cpm + 0.5))
     p = one-sided Mann-Whitney U on the target's CPM, construct cells < control cells
     keep = (not measurable) or (ontarget_lfc <= MAX_ONTARGET_LFC and p < P_MAX)
2. Per knockdown (target_gene_id): cells pooled from kept constructs only.
     n_cells_kept; knockdown_keep = n_cells_kept >= MIN_CELLS; efficiency = ontarget_lfc of the pooled kept cells
     (NaN if not measurable)
3. Sanity check, 50 random measurable constructs from each of three groups: dropped; kept with -1 < ontarget_lfc <= -0.5;
   kept with ontarget_lfc <= -1. Per construct: genes other than the target with |log2FC| > 0.5
   (log2((construct mean CPM + 0.5) / (control mean CPM + 0.5))) over every readout gene, and, for interpretation only,
   over the genes >= 5 CPM in controls (lowly expressed genes cross 0.5 by noise alone).

Writes results/qc/<source>/: guide_qc.csv, knockdown_qc.csv, ontarget_hist.png, perturbation_qc.json.
Usage: ~/venv/bin/python scripts/qc_perturbations.py --dataset K562|RPE1|HepG2|Jurkat|H1
"""
import argparse
import json
import sys
import time
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).parent))
from qc_replogle_umi import CONTROL, REFERENCE, SOURCES                                    # noqa: E402

MIN_CELLS, MIN_CTRL_CPM, MAX_ONTARGET_LFC, P_MAX, PSEUDO = 30, 5.0, -0.5, 0.05, 0.5
CHUNK = 50000
N_SANITY, SANITY_LFC, SEED = 50, 0.5, 2026


def gene_ids(var):
    """Ensembl ID per readout gene: the var index for Replogle / Nadig files, var['gene_id'] for H1."""
    if 'gene_name' in var: return var.index.astype(str).to_numpy(), var['gene_name'].astype(str).to_numpy()
    return var['gene_id'].astype(str).to_numpy(), var.index.astype(str).to_numpy()


def chunks(X, n_obs, n_vars):
    ip = X['indptr']
    for s in range(0, n_obs, CHUNK):
        e = min(s + CHUNK, n_obs); p = ip[s:e + 1]; a, b = int(p[0]), int(p[-1])
        yield s, e, sp.csr_matrix((X['data'][a:b], X['indices'][a:b], p - a), shape=(e - s, n_vars))


def main(dataset):
    _, out, path, title = SOURCES[dataset]
    t0 = time.time()
    with h5py.File(path, 'r') as f:
        obs = read_elem(f['obs']); var = read_elem(f['var']); X = f['X']; n_vars = int(X.attrs['shape'][1]); n_obs = len(obs)
        ids, symbols = gene_ids(var); col_of = {g: i for i, g in enumerate(ids)}
        ctrl = obs['is_control'].to_numpy(bool); lib = obs['n_umi'].to_numpy(float)
        tid = obs['target_gene_id'].astype(str).to_numpy(); guide = obs['guide_id'].astype(str).to_numpy()
        cell_col = np.array([col_of.get(t, -1) if not c else -1 for t, c in zip(tid, ctrl)])
        own_cpm = np.full(n_obs, np.nan, np.float32); ctrl_parts = []
        # REFERENCE 'rest' (Feng genome-wide, ~no non-targeting cells): each target is tested against every other cell of the
        # screen, so pass 1 keeps every cell's CPM of every target gene (cells x targets, sparse) and the all-cell mean CPM
        rest = REFERENCE.get(dataset) == 'rest'
        tcols = np.unique(cell_col[cell_col >= 0]); tpos = {c: i for i, c in enumerate(tcols)}; rest_parts = []; all_sum = np.zeros(n_vars)
        for s, e, M in chunks(X, n_obs, n_vars):                    # pass 1: controls' counts; each cell's own-target CPM
            scale = 1e6 / lib[s:e]; cm = ctrl[s:e]
            if rest:
                P = sp.diags(scale) @ M; rest_parts.append(P[:, tcols].astype(np.float32)); all_sum += np.asarray(P.sum(0)).ravel()
            elif cm.any(): ctrl_parts.append(sp.diags(scale[cm]) @ M[np.flatnonzero(cm)])
            rows = np.flatnonzero(cell_col[s:e] >= 0)
            if len(rows): own_cpm[s + rows] = np.asarray(M[rows, cell_col[s:e][rows]]).ravel() * scale[rows]
            print(f'  pass 1: {e:,} / {n_obs:,} cells ({time.time() - t0:.0f}s)', flush=True)
    if rest:
        R = sp.vstack(rest_parts).tocsc(); n_ctrl = n_obs; ctrl_mean = all_sum / n_obs    # ctrl_mean: sanity check and plots only
    else:
        C = sp.vstack(ctrl_parts).tocsc().astype(np.float32); n_ctrl = C.shape[0]
        ctrl_mean = np.asarray(C.mean(0)).ravel()
    ref_cache = {}

    def reference(j, t):
        """CPM of gene j in the reference cells for target t: the controls, or (rest) every cell not carrying target t."""
        if rest:
            if t not in ref_cache: ref_cache[t] = R[:, tpos[j]].toarray().ravel()[tid != t]
            return ref_cache[t]
        if j not in ctrl_col_cache: ctrl_col_cache[j] = C[:, j].toarray().ravel()
        return ctrl_col_cache[j]

    # 1. constructs
    cells = pd.DataFrame(dict(guide=guide, tid=tid, target=obs['target_gene'].astype(str).to_numpy(), col=cell_col, cpm=own_cpm))[~ctrl]
    cells['row'] = np.flatnonzero(~ctrl)
    per_target = cells.groupby('guide')['tid'].nunique()
    if (per_target > 1).any(): raise ValueError(f'constructs with more than one target: {per_target[per_target > 1].index[:5].tolist()}')
    ctrl_col_cache = {}
    rows = []
    for g, d in cells.groupby('guide', sort=False):
        j = int(d['col'].iloc[0]); in_var = j >= 0
        ref_vals = reference(j, d['tid'].iloc[0]) if in_var else None
        c_cpm = (float(ref_vals.mean()) if rest else float(ctrl_mean[j])) if in_var else np.nan
        measurable = bool(in_var and c_cpm >= MIN_CTRL_CPM)
        k_mean = float(d['cpm'].mean()) if in_var else np.nan
        lfc = np.log2((k_mean + PSEUDO) / (c_cpm + PSEUDO)) if in_var else np.nan
        p = np.nan
        if in_var:
            p = float(mannwhitneyu(d['cpm'].to_numpy(), ref_vals, alternative='less', method='asymptotic').pvalue)
        keep = (not measurable) or (lfc <= MAX_ONTARGET_LFC and p < P_MAX)
        rows.append(dict(guide_id=g, target_gene=d['target'].iloc[0], target_gene_id=d['tid'].iloc[0], n_cells=len(d), ctrl_cpm=c_cpm,
                         measurable=measurable, ontarget_lfc=lfc, p=p, keep=bool(keep)))
    guides = pd.DataFrame(rows)
    guides.to_csv(out / 'guide_qc.csv', index=False)

    # 2. knockdowns: kept constructs pooled
    kept_cells = cells[cells['guide'].isin(set(guides.loc[guides.keep, 'guide_id']))]
    first = cells.groupby('tid', sort=False)[['target', 'col']].first()                      # grouped, not one scan per knockdown
    per_guide = guides.groupby('target_gene_id')
    kd = pd.DataFrame(dict(target_gene=first['target'], target_gene_id=first.index,
                           n_constructs_kept=per_guide['keep'].sum().reindex(first.index).fillna(0).astype(int),
                           n_cells_kept=kept_cells.groupby('tid').size().reindex(first.index).fillna(0).astype(int),
                           measurable=per_guide['measurable'].any().reindex(first.index).fillna(False).astype(bool)))
    kept_mean = kept_cells.groupby('tid')['cpm'].mean().reindex(first.index)
    j = first['col'].to_numpy(int)
    ref_mean = (np.array([float(ref_cache[t].mean()) if t in ref_cache else np.nan for t in first.index]) if rest
                else np.where(j >= 0, ctrl_mean[np.maximum(j, 0)], np.nan))
    kd['efficiency'] = np.where(kd.measurable & (kd.n_cells_kept > 0),
                                np.log2((kept_mean.to_numpy() + PSEUDO) / (ref_mean + PSEUDO)), np.nan)
    kd['knockdown_keep'] = kd.n_cells_kept >= MIN_CELLS
    kd = kd[['target_gene', 'target_gene_id', 'n_constructs_kept', 'n_cells_kept', 'efficiency', 'measurable', 'knockdown_keep']]
    kd = kd.reset_index(drop=True).sort_values(['knockdown_keep', 'efficiency'], ascending=[False, True])
    kd.to_csv(out / 'knockdown_qc.csv', index=False)

    # 3. plot
    m = guides[guides.measurable]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    lo = np.nanpercentile(m.ontarget_lfc, 0.5) if len(m) else -5
    ax.hist(m.ontarget_lfc.clip(lower=lo), bins=80, color='#2a78d6', histtype='stepfilled', alpha=0.8)
    ax.axvline(MAX_ONTARGET_LFC, ls='--', color='grey', lw=1); ax.axvline(-1, ls=':', color='grey', lw=1)
    ax.text(MAX_ONTARGET_LFC, ax.get_ylim()[1] * 0.95, f' {MAX_ONTARGET_LFC}', color='grey', fontsize=8)
    ax.set_xlabel('on-target log2 fold change of mean CPM (construct vs controls)'); ax.set_ylabel('constructs')
    ax.set_title(f'{title}: {len(m):,} measurable constructs (dotted: -1, a 50% drop)')
    fig.tight_layout(); fig.savefig(out / 'ontarget_hist.png', dpi=150); plt.close(fig)

    # 4. sanity check: downstream changes in three groups of measurable constructs
    groups = {'dropped': m[~m.keep], 'band_-1_to_-0.5': m[m.keep & (m.ontarget_lfc > -1)], 'le_-1': m[m.keep & (m.ontarget_lfc <= -1)]}
    pick = {k: (v.sample(min(N_SANITY, len(v)), random_state=SEED)['guide_id'].tolist() if len(v) else []) for k, v in groups.items()}
    chosen = [g for v in pick.values() for g in v]; slot = {g: i for i, g in enumerate(chosen)}
    sums = np.zeros((len(chosen), n_vars)); counts = np.zeros(len(chosen))
    if chosen:
        cell_slot = np.full(n_obs, -1); sel = cells[cells['guide'].isin(slot)]
        cell_slot[sel['row'].to_numpy()] = sel['guide'].map(slot).to_numpy()
        with h5py.File(path, 'r') as f:
            for s, e, M in chunks(f['X'], n_obs, n_vars):         # pass 2: full profiles of the chosen constructs only
                r = np.flatnonzero(cell_slot[s:e] >= 0)
                if not len(r): continue
                P = sp.diags(1e6 / lib[s:e][r]) @ M[r]
                k = cell_slot[s:e][r]
                for i in np.unique(k):
                    sums[i] += np.asarray(P[k == i].sum(0)).ravel(); counts[i] += int((k == i).sum())
    expressed = ctrl_mean >= MIN_CTRL_CPM
    n_all, n_expr = {}, {}
    for g in chosen:
        i = slot[g]; lfc = np.abs(np.log2((sums[i] / counts[i] + PSEUDO) / (ctrl_mean + PSEUDO))) > SANITY_LFC
        j = int(guides.loc[guides.guide_id == g, 'target_gene_id'].map(col_of).iloc[0]); lfc[j] = False
        n_all[g] = int(lfc.sum()); n_expr[g] = int((lfc & expressed).sum())
    med = lambda d, gs: float(np.median([d[g] for g in gs])) if gs else None
    sanity = {k: dict(n=len(gs), median_all_genes=med(n_all, gs), median_expressed_genes=med(n_expr, gs)) for k, gs in pick.items()}

    loose = m[m.keep & (m.ontarget_lfc > -1)]
    near = kd[(kd.n_cells_kept >= 20) & (kd.n_cells_kept < MIN_CELLS)]
    summary = dict(source=out.name, controls=int(n_ctrl), readout_genes=int(n_vars),
                   reference=('rest of screen: every cell not carrying the target (almost no non-targeting cells); controls = all cells'
                              if rest else 'non-targeting controls'),
                   constructs_total=int(len(guides)), constructs_measurable=int(guides.measurable.sum()), constructs_kept=int(guides.keep.sum()),
                   constructs_kept_unmeasurable=int((guides.keep & ~guides.measurable).sum()),
                   constructs_kept_between_minus1_and_minus05=int(len(loose)),
                   knockdowns_total=int(len(kd)), knockdowns_kept=int(kd.knockdown_keep.sum()),
                   knockdowns_20_29_cells=int(len(near)),
                   median_efficiency_kept=float(kd.loc[kd.knockdown_keep, 'efficiency'].median()),
                   sanity_downstream_genes=dict(groups=sanity, threshold=SANITY_LFC, expressed_genes=int(expressed.sum())),
                   parameters=dict(MIN_CELLS=MIN_CELLS, MIN_CTRL_CPM=MIN_CTRL_CPM, MAX_ONTARGET_LFC=MAX_ONTARGET_LFC, P_MAX=P_MAX),
                   minutes=round((time.time() - t0) / 60, 1))
    (out / 'perturbation_qc.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--dataset', required=True, choices=list(SOURCES))
    main(p.parse_args().dataset)

"""Knockdown strength vs downstream effect, per training dataset (descriptive, read-only).

  cd ~/vcc && ~/venv/bin/python analysis/knockdown_strength/ks.py

Normalisation = the transformer's effect tables (scripts/build_effects_v4.py, lib/effects_v4): per-cell CPM (counts /
n_umi x 1e6), mean per-cell CPM per group, delta = log2((perturbed mean CPM + 0.5) / (control mean CPM + 0.5)).
Perturbed and control means come from results/qc/<source>/construct_totals.npz (per construct: n, sum CPM, sum CPM^2,
depth; last row = the source's non-targeting controls; cells after cell QC = lib/qc_cells). Every construct in
guide_qc.csv is pooled by target gene, so every perturbation is included, also those the knockdown QC later drops.
knockdown_pct = 100 (1 - perturbed mean CPM / control mean CPM) of the target gene. Unmeasurable: target not among the
source's genes or detected in < 10% of control cells.
n_DE_z: the pipeline's DE count for every perturbation (BH 5%, two-sided normal p from delta / SE, delta-method SE,
target excluded). n_DE_wilcoxon: Wilcoxon rank-sum on per-cell CPM vs 1,000 random control cells, BH 5% over the
source's genes, target excluded, for a stratified sample (<= 30 measurable perturbations per bin, n_cells >= 20).
"""
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
from scipy.stats import false_discovery_control, mannwhitneyu, norm, pearsonr, spearmanr

V = Path.home() / 'vcc'
OUT = V / 'analysis/knockdown_strength'
DATASETS = {'H1 (Arc)': 'arc_h1', 'K562 (Replogle)': 'replogle_k562_gw', 'RPE1 (Replogle)': 'replogle_rpe1',
            'HepG2 (Nadig)': 'nadig_hepg2', 'Jurkat (Nadig)': 'nadig_jurkat'}
PSEUDO, SE_FLOOR, FDR = 0.5, 0.05, 0.05
BINS = ['Strong', 'Medium', 'Weak', 'Minimal']
COLORS = {'Strong': '#1b7837', 'Medium': '#5aae61', 'Weak': '#f1a340', 'Minimal': '#d73027'}
N_WILCOX_PER_BIN, N_CTRL, MIN_CELLS_WILCOX = 30, 1000, 20


def bin_of(p):
    return 'Strong' if p >= 80 else 'Medium' if p >= 50 else 'Weak' if p >= 30 else 'Minimal'


def bh_n(p):
    if not len(p): return 0
    return int((false_discovery_control(p, method='bh') < FDR).sum())


def per_perturbation(source):
    qc = V / 'results/qc' / source
    guides = pd.read_csv(qc / 'guide_qc.csv'); guides['target_gene_id'] = guides['target_gene_id'].astype(str)
    kq = pd.read_csv(qc / 'knockdown_qc.csv'); kq['target_gene_id'] = kq['target_gene_id'].astype(str)
    tot = np.load(qc / 'construct_totals.npz', allow_pickle=False)
    assert list(tot['guides'].astype(str)) == guides.guide_id.astype(str).tolist()
    with h5py.File(V / f'lib/qc_cells/{source}_qc.h5ad', 'r') as f: var = read_elem(f['var'])
    src_ids = var.index.astype(str).to_numpy() if 'gene_name' in var else var['gene_id'].astype(str).to_numpy()
    col = {g: i for i, g in enumerate(src_ids)}
    targets = guides.drop_duplicates('target_gene_id')[['target_gene', 'target_gene_id']].reset_index(drop=True)
    row = {t: i for i, t in enumerate(targets.target_gene_id)}
    P = sp.csr_matrix((np.ones(len(guides)), ([row[t] for t in guides.target_gene_id], np.arange(len(guides)))),
                      shape=(len(targets), len(guides)))
    s1, s2, n, depth = P @ tot['s1'][:-1], P @ tot['s2'][:-1], P @ tot['n'][:-1], P @ tot['depth'][:-1]
    n0 = tot['n'][-1]; m0 = tot['s1'][-1] / n0
    nn = np.maximum(n, 1)[:, None]; m1 = s1 / nn

    def relvar(s1_, s2_, n_, depth_, mean):                     # build_effects_v4.stats: delta-method, Poisson floor
        nn_ = np.maximum(n_, 1)[:, None]
        var_mean = np.maximum(s2_ / nn_ - mean ** 2, 0) * nn_ / np.maximum(nn_ - 1, 1) / nn_
        u = 1e6 / np.maximum(depth_, 1)[:, None]
        return np.maximum(var_mean, (np.maximum(mean, 0) + u) * u) / (mean + PSEUDO) ** 2
    rv1 = relvar(s1, s2, n, depth, m1); rv0 = relvar(tot['s1'][-1:], tot['s2'][-1:], tot['n'][-1:], tot['depth'][-1:], m0[None])
    delta = np.log2((m1 + PSEUDO) / (m0 + PSEUDO)); se = np.maximum(np.sqrt(rv1 + rv0) / np.log(2), SE_FLOOR)
    det = tot['detect'] / n0
    recs = []
    for i, t in enumerate(targets.itertuples()):
        j = col.get(t.target_gene_id, -1)
        keep = np.ones(len(src_ids), bool)
        if j >= 0: keep[j] = False
        d = delta[i, keep]; l2 = float(np.sqrt((d ** 2).sum()))
        l2_nc = float(np.sqrt(max((d ** 2).sum() - (se[i, keep] ** 2).sum(), 0)))
        p = 2 * norm.sf(np.abs(d / se[i, keep]))
        r = dict(target_gene=t.target_gene, target_gene_id=t.target_gene_id, n_cells=int(n[i]), target_in_genes=j >= 0,
                 control_mean_cpm=float(m0[j]) if j >= 0 else np.nan, control_detect_frac=float(det[j]) if j >= 0 else np.nan,
                 pert_mean_cpm=float(m1[i, j]) if j >= 0 else np.nan, downstream_l2=l2, downstream_l2_noise_corrected=l2_nc,
                 n_DE_z=bh_n(p), target_delta=float(delta[i, j]) if j >= 0 else np.nan,
                 target_share=float(delta[i, j] ** 2 / (delta[i] ** 2).sum()) if j >= 0 else np.nan)
        r['measurable'] = bool(j >= 0 and det[j] >= 0.10 and m0[j] > 0)
        r['knockdown_pct'] = 100 * (1 - r['pert_mean_cpm'] / r['control_mean_cpm']) if j >= 0 and m0[j] > 0 else np.nan
        r['bin'] = bin_of(r['knockdown_pct']) if r['measurable'] else 'Unmeasurable'
        recs.append(r)
    df = pd.DataFrame(recs)
    df['pipeline_keeps'] = df.target_gene_id.map(kq.set_index('target_gene_id').knockdown_keep.astype(bool).to_dict()).fillna(False)
    return df, guides, src_ids, col


def wilcoxon_sample(source, df, guides, src_ids, col, rng):
    path = V / f'lib/qc_cells/{source}_qc.h5ad'
    with h5py.File(path, 'r') as f:
        obs = read_elem(f['obs']); X = f['X']; ip = X['indptr'][:]; G = int(X.attrs['shape'][1])
        lib = obs['n_umi'].to_numpy(float); gid = obs['guide_id'].astype(str).to_numpy()
        by_target = guides.groupby('target_gene_id').guide_id.apply(lambda s: set(s.astype(str)))
        pick = []
        for b in BINS:
            c = df[(df.bin == b) & (df.n_cells >= MIN_CELLS_WILCOX)]
            pick += list(c.sample(min(N_WILCOX_PER_BIN, len(c)), random_state=0).index)
        ctrl_rows = rng.choice(np.flatnonzero(obs['is_control'].to_numpy(bool)), N_CTRL, replace=False)
        kd_rows = {i: np.flatnonzero(np.isin(gid, list(by_target[df.at[i, 'target_gene_id']]))) for i in pick}
        need = np.unique(np.concatenate([ctrl_rows] + list(kd_rows.values())))
        # one streaming pass in large row blocks (row-by-row reads decompress a chunk each time)
        blocks, t0 = [], time.time()
        for s0 in range(0, len(obs), 50000):
            e0 = min(s0 + 50000, len(obs)); r = need[(need >= s0) & (need < e0)]
            if not len(r): continue
            a, b = int(ip[s0]), int(ip[e0])
            M = sp.csr_matrix((X['data'][a:b], X['indices'][a:b], ip[s0:e0 + 1] - a), shape=(e0 - s0, G))[r - s0]
            blocks.append(M.astype(np.float32))
        C = sp.vstack(blocks).tocsr(); pos = {k: i for i, k in enumerate(need)}
        print(f'  read {len(need):,} cells in {time.time() - t0:.0f}s', flush=True)

        def rows(r):
            idx = np.array([pos[k] for k in r])
            return (sp.diags((1e6 / lib[r]).astype(np.float32)) @ C[idx]).toarray()
        ctrl = rows(ctrl_rows)
        out = {}; t0 = time.time()
        for i in pick:
            t = df.at[i, 'target_gene_id']
            A = rows(kd_rows[i]); keep = np.ones(G, bool); j = col.get(t, -1)
            if j >= 0: keep[j] = False
            tested = keep & ((A > 0).any(0) | (ctrl > 0).any(0))          # all-zero genes: p = 1, still counted in BH
            p = np.ones(G); p[tested] = mannwhitneyu(A[:, tested], ctrl[:, tested], axis=0).pvalue
            out[i] = bh_n(np.nan_to_num(p[keep], nan=1.0))
        print(f'  wilcoxon: {len(pick)} perturbations in {time.time() - t0:.0f}s', flush=True)
    return pd.Series(out)


def partial_spearman(x, y, z):
    rx, ry, rz = (pd.Series(v).rank().to_numpy() for v in (x, y, z))
    res = lambda a: a - np.polyval(np.polyfit(rz, a, 1), rz)
    return pearsonr(res(rx), res(ry))


def iqr_str(s):
    s = s.dropna()
    return f'{s.median():.3g} [{s.quantile(.25):.3g}, {s.quantile(.75):.3g}]' if len(s) else ''


def main():
    OUT.mkdir(parents=True, exist_ok=True); rng = np.random.default_rng(2026)
    counts, summ, cors = [], [], []
    for name, source in DATASETS.items():
        t0 = time.time(); print(name, flush=True)
        df, guides, src_ids, col = per_perturbation(source); print(f'  per-perturbation table {len(df)} rows ({time.time() - t0:.0f}s)', flush=True)
        df['n_DE_wilcoxon'] = wilcoxon_sample(source, df, guides, src_ids, col, rng)
        df.insert(0, 'dataset', name); df.to_csv(OUT / f'per_perturbation_{source}.csv', index=False)
        m = df[df.measurable]; c = {'dataset': name, 'perturbations': len(df), 'measurable': len(m)}
        for b in BINS: c[b] = f'{(m.bin == b).sum()} ({100 * (m.bin == b).mean():.1f}%)'
        c['Unmeasurable'] = int((~df.measurable).sum()); counts.append(c)
        for b in BINS + ['Unmeasurable']:
            g = df[df.bin == b]
            w = g.n_DE_wilcoxon.dropna()
            summ.append(dict(dataset=name, bin=b, n=len(g), kept_by_pipeline=int(g.pipeline_keeps.sum()), median_n_cells=g.n_cells.median(),
                             downstream_l2=iqr_str(g.downstream_l2), downstream_l2_noise_corrected=iqr_str(g.downstream_l2_noise_corrected),
                             n_DE_z=iqr_str(g.n_DE_z), target_share=iqr_str(g.target_share),
                             null_like_z=f'{100 * (g.n_DE_z <= 5).mean():.0f}%' if len(g) else '',
                             n_wilcoxon=len(w), n_DE_wilcoxon=iqr_str(w), null_like_wilcoxon=f'{100 * (w <= 5).mean():.0f}%' if len(w) else ''))
        x, z = m.knockdown_pct.to_numpy(), np.log(m.n_cells.to_numpy())
        for yname in ('downstream_l2', 'downstream_l2_noise_corrected', 'n_DE_z'):
            y = m[yname].to_numpy(); r, p = spearmanr(x, y); pr, pp = partial_spearman(x, y, z)
            cors.append(dict(dataset=name, y=yname, n=len(m), spearman=r, p=p, partial_spearman_given_log_n_cells=pr, partial_p=pp))
        mw = m.dropna(subset=['n_DE_wilcoxon'])
        if len(mw) > 10:
            r, p = spearmanr(mw.knockdown_pct, mw.n_DE_wilcoxon); pr, pp = partial_spearman(mw.knockdown_pct, mw.n_DE_wilcoxon, np.log(mw.n_cells))
            cors.append(dict(dataset=name, y='n_DE_wilcoxon (sample)', n=len(mw), spearman=r, p=p, partial_spearman_given_log_n_cells=pr, partial_p=pp))
            cors.append(dict(dataset=name, y='agreement: n_DE_z vs n_DE_wilcoxon', n=len(mw), spearman=spearmanr(mw.n_DE_z, mw.n_DE_wilcoxon)[0]))
        # scatter
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        for b in BINS:
            g = m[m.bin == b]; ax.scatter(g.knockdown_pct.clip(-100, 100), g.downstream_l2, s=6, alpha=0.55, c=COLORS[b], label=f'{b} (n={len(g)})', linewidths=0)
        ax.set_yscale('log'); ax.set_xlabel('knockdown_pct (clipped to -100..100)'); ax.set_ylabel('downstream_l2 (log scale, target excluded)')
        ax.set_title(f'{name}: knockdown vs downstream effect'); ax.axvline(30, c='grey', lw=0.6, ls=':'); ax.axvline(50, c='grey', lw=0.6, ls=':'); ax.axvline(80, c='grey', lw=0.6, ls=':')
        ax.legend(fontsize=8, markerscale=2, frameon=False); fig.tight_layout(); fig.savefig(OUT / f'scatter_{source}.png', dpi=130); plt.close(fig)
        print(f'  done in {time.time() - t0:.0f}s', flush=True)
    pd.DataFrame(counts).to_csv(OUT / 'bin_counts.csv', index=False)
    pd.DataFrame(summ).to_csv(OUT / 'per_bin_summary.csv', index=False)
    pd.DataFrame(cors).to_csv(OUT / 'correlations.csv', index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
    print(pd.DataFrame(counts).to_string(index=False)); print(pd.DataFrame(summ).to_string(index=False)); print(pd.DataFrame(cors).round(3).to_string(index=False))


if __name__ == '__main__':
    main()

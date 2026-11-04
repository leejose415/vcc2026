"""Model 7, Step 0: per-line effect tables in the scorer's space, rebuilt from raw counts where cells exist.

  cd ~/vcc && PYTHONPATH=scripts ~/venv/bin/python analysis/model7/m7_tables.py <line> [<line> ...]
  lines: K562 RPE1 HepG2 Jurkat H1 HCT116 HEK293T CD4_Rest CD4_Stim8hr CD4_Stim48hr

Scorer space: S = log1p(5e4 x P), P = summed raw counts / summed total counts over a group's cells (cell-eval2
bulk_lognorm), totals over the genes the line measures on the challenge panel. Delta = S_k - S_ctrl.
Groups: every construct with guide_qc keep, pooled by target gene; controls = the line's non-targeting cells.
Variance of S per group and gene: the scorer's delete-1 cell jackknife, linearised (first order in one cell's share):
  var_g = (n-1)/n * (TS / (1 + TS P_g))^2 / S_tot^2 * sum_i (y_ig - P_g l_i)^2,
from per-cell sufficient statistics (sum y, sum y^2, sum y l, sum l, sum l^2) in one streaming pass.
Table lines (no per-cell counts here):
  HCT116 / HEK293T (X-Atlas/Orion pseudobulk sums): P from summed raw counts (exact); variance from the per-cell CPM
    second moments (delta method, mean-of-ratios) - an approximation.
  CD4 Rest / Stim8hr / Stim48hr (Marson tables: log2 ratio of pseudobulk CPM, pseudocount 0.5, and control CPM):
    P_k CPM = 2^lfc (C + 0.5) - 0.5 mapped into S; variance carried by the delta method from the table se.
Knockdown % = 100 (1 - P_target / C_target) (CD4: 100 (1 - 2^target_lfc)); unmeasurable = target not measured or
detected in < 10% of control cells (CD4: control CPM < 1); weights: < 30% dropped (0), 30-50% and unmeasurable 0.5.
Output: lib/vcc_2026_5/m7/<line>.npz (float16 Delta and variance on the 18,533-gene panel, NaN = not measured).
"""
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

V = Path.home() / 'vcc'
OUT = V / 'lib/vcc_2026_5/m7'
TS = 5e4
CELL = {'K562': 'replogle_k562_gw', 'RPE1': 'replogle_rpe1', 'HepG2': 'nadig_hepg2', 'Jurkat': 'nadig_jurkat', 'H1': 'arc_h1'}
TABLE = {'HCT116': 'orion_hct116', 'HEK293T': 'orion_hek293t'}
CD4 = {'CD4_Rest': 'Rest', 'CD4_Stim8hr': 'Stim8hr', 'CD4_Stim48hr': 'Stim48hr'}


def log(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)


def weights(kd, measurable):
    return np.where(~measurable | ~np.isfinite(kd), 0.5, np.where(kd < 30, 0.0, np.where(kd < 50, 0.5, 1.0)))


def save(name, tg, S, varS, S0, var0, n, n0, measured, kd, measurable, extra):
    G = S.shape[1]
    S[:, ~measured] = np.nan; varS[:, ~measured] = np.nan
    D = (S - S0[None, :]).astype(np.float32); vD = (varS + var0[None, :]).astype(np.float32)
    w = weights(kd, measurable)
    R = np.nanmean(np.where(w[:, None] > 0, D, np.nan), 0)
    keep = w > 0
    noise = float(np.nansum(vD[keep]) / np.nansum(D[keep] ** 2))
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(OUT / f'{name}.npz', tg=np.asarray(tg, dtype=str), D=D.astype(np.float16), vD=vD.astype(np.float16), S0=S0, var0=var0,
             n=np.asarray(n, float), n0=float(n0), measured=measured, kd_pct=kd, measurable=measurable, w=w, R=R, **extra)
    rep = dict(line=name, knockdowns=len(tg), kept_w_gt0=int(keep.sum()), dropped_lt30=int((w == 0).sum()), half_weight=int((w == 0.5).sum()),
               measured_genes=int(measured.sum()), median_cells=float(np.median(n)), controls=float(n0), noise_share_of_sum_Delta2=noise)
    (OUT / f'{name}_report.json').write_text(json.dumps(rep, indent=1)); log(json.dumps(rep))


def cell_line(name, src, genes, var_panel):
    qc = V / 'results/qc' / src
    guides = pd.read_csv(qc / 'guide_qc.csv'); guides = guides[guides.keep.astype(bool)].copy()
    guides['target_gene_id'] = guides['target_gene_id'].astype(str); guides['guide_id'] = guides['guide_id'].astype(str)
    targets = guides.drop_duplicates('target_gene_id')[['target_gene', 'target_gene_id']].reset_index(drop=True)
    tpos = {t: i for i, t in enumerate(targets.target_gene_id)}; g2t = dict(zip(guides.guide_id, guides.target_gene_id))
    path = V / f'lib/qc_cells/{src}_qc.h5ad'
    with h5py.File(path, 'r') as f:
        obs = read_elem(f['obs']); var = read_elem(f['var']); X = f['X']; ip = X['indptr'][:]; Gs = int(X.attrs['shape'][1])
        src_ids = var.index.astype(str).to_numpy() if 'gene_name' in var else var['gene_id'].astype(str).to_numpy()
        pid = {g: i for i, g in enumerate(var_panel.gene_id) if g}; to_panel = np.array([pid.get(g, -1) for g in src_ids])
        ok = to_panel >= 0; keepcols = np.flatnonzero(ok)
        ctrl = obs['is_control'].to_numpy(bool); gid = obs['guide_id'].astype(str).to_numpy()
        K = len(targets); grp = np.array([K if c else tpos.get(g2t.get(g, ''), -1) for c, g in zip(ctrl, gid)])
        Gk = len(keepcols); s1 = np.zeros((K + 1, Gk)); s2 = np.zeros((K + 1, Gk)); syl = np.zeros((K + 1, Gk))
        nl, nl2, n = np.zeros(K + 1), np.zeros(K + 1), np.zeros(K + 1); det = np.zeros(Gk); t0 = time.time()
        for s0 in range(0, len(obs), 50000):
            e0 = min(s0 + 50000, len(obs)); g = grp[s0:e0]; r = np.flatnonzero(g >= 0)
            if not len(r): continue
            a, b = int(ip[s0]), int(ip[e0])
            M = sp.csr_matrix((X['data'][a:b], X['indices'][a:b], ip[s0:e0 + 1] - a), shape=(e0 - s0, Gs))[r][:, keepcols].astype(np.float64)
            l = np.asarray(M.sum(1)).ravel(); gg = g[r]
            Gm = sp.csr_matrix((np.ones(len(r)), (gg, np.arange(len(r)))), shape=(K + 1, len(r)))
            s1 += (Gm @ M).toarray(); s2 += (Gm @ M.multiply(M)).toarray(); syl += (Gm @ sp.diags(l) @ M).toarray()
            nl += np.bincount(gg, weights=l, minlength=K + 1); nl2 += np.bincount(gg, weights=l * l, minlength=K + 1); n += np.bincount(gg, minlength=K + 1)
            c = gg == K
            if c.any(): det += np.asarray((M[np.flatnonzero(c)] > 0).sum(0)).ravel()
            log(f'  {name}: {e0:,}/{len(obs):,} cells ({time.time() - t0:.0f}s)')
    Ptot = np.maximum(nl, 1)[:, None]; P = s1 / Ptot
    resid = np.maximum(s2 - 2 * P * syl + P ** 2 * nl2[:, None], 0)
    varS = ((n - 1) / np.maximum(n, 1))[:, None] * (TS / (1 + TS * P)) ** 2 / Ptot ** 2 * resid
    S = np.log1p(TS * P)
    G = len(genes)
    def panelize(A):
        out = np.full((A.shape[0], G), np.nan, np.float32); out[:, to_panel[keepcols]] = A; return out
    Sp, Vp = panelize(S), panelize(varS); measured = np.zeros(G, bool); measured[to_panel[keepcols]] = True
    P0 = panelize(P[-1:])[0]; detp = np.full(G, np.nan); detp[to_panel[keepcols]] = det / max(n[-1], 1)
    gp = {gg: i for i, gg in enumerate(genes)}; own = np.array([gp.get(t, -1) for t in targets.target_gene])
    kcol = np.full(G, -1); kcol[to_panel[keepcols]] = np.arange(len(keepcols))
    kd = np.array([100 * (1 - P[i, kcol[j]] / P0[j]) if j >= 0 and measured[j] and P0[j] > 0 else np.nan for i, j in enumerate(own)])
    measurable = np.array([j >= 0 and measured[j] and detp[j] >= 0.10 for j in own])
    save(name, targets.target_gene.astype(str).to_numpy(), Sp[:-1], Vp[:-1], Sp[-1], Vp[-1], n[:-1], n[-1], measured, kd, measurable,
         dict(ctrl_detect=detp, P0=P0))


def table_line(name, src, genes):
    with h5py.File(V / f'lib/orion/{name}_pseudobulk.h5ad', 'r') as f:
        obs = read_elem(f['obs']); var = read_elem(f['var'])
        names = obs.index.astype(str).to_numpy(); ci = int(np.flatnonzero(names == 'non-targeting')[0])
        rows = np.array([i for i in range(len(names)) if i != ci])
        ncell = obs['n_cells'].to_numpy(float); measured = var['measured'].to_numpy(bool)
        sc, s1c, s2c, nex = f['layers/sum_counts'], f['layers/sum_cpm'], f['layers/sum_cpm2'], f['layers/n_expr']
        def block(ds, r): return ds[r].astype(np.float64)
        S, Vv = np.full((len(rows), len(genes)), np.nan, np.float32), np.full((len(rows), len(genes)), np.nan, np.float32)
        def stats(cnt, cpm, cpm2, n):
            cnt = np.where(measured[None, :], cnt, 0); tot = cnt.sum(1, keepdims=True); P = cnt / np.maximum(tot, 1)
            mr = cpm / 1e6; m2 = cpm2 / 1e12; nn = np.maximum(n, 2)[:, None]
            varP = np.maximum(m2 - mr ** 2 / nn, 0) / (nn * (nn - 1))
            return np.log1p(TS * P), (TS / (1 + TS * P)) ** 2 * varP, P
        S0, V0, P0 = (x[0] for x in stats(block(sc, [ci]), block(s1c, [ci]), block(s2c, [ci]), ncell[[ci]]))
        det = block(nex, [ci])[0] / ncell[ci]
        for a in range(0, len(rows), 1000):
            r = rows[a:a + 1000]
            s_, v_, _ = stats(block(sc, r), block(s1c, r), block(s2c, r), ncell[r]); S[a:a + 1000], Vv[a:a + 1000] = s_, v_
            log(f'  {name}: {a + len(r)}/{len(rows)}')
    assert list(var.index.astype(str)) == list(genes)
    tg = names[rows]; gp = {g: i for i, g in enumerate(genes)}; own = np.array([gp.get(t, -1) for t in tg])
    kd = np.array([100 * (1 - np.expm1(float(S[i, j])) / TS / P0[j]) if j >= 0 and measured[j] and P0[j] > 0 else np.nan for i, j in enumerate(own)])
    measurable = np.array([j >= 0 and measured[j] and det[j] >= 0.10 for j in own])
    save(name, tg, S, Vv, S0, V0, ncell[rows], ncell[ci], measured, kd, measurable, dict(ctrl_detect=det, P0=P0))


def cd4_line(name, cond, genes):
    with h5py.File(V / f'lib/tcell/TCELL_{cond}_effects.h5ad', 'r') as f:
        obs = read_elem(f['obs']); var = read_elem(f['var']); L = f['X'][:].astype(np.float32); SE = f['layers/se'][:].astype(np.float32)
    assert list(var.index.astype(str)) == list(genes)
    measured = var['measured'].to_numpy(bool); C = var['ctrl_cpm'].to_numpy(float)
    Pc = np.maximum(np.exp2(L) * (C + 0.5) - 0.5, 0)                      # perturbed CPM implied by the table's lfc
    S = np.log1p(TS * Pc / 1e6); S0 = np.log1p(TS * C / 1e6)
    dS = TS / 1e6 * np.log(2) * np.exp2(L) * (C + 0.5) / (1 + TS * Pc / 1e6)  # dS / dlfc
    Vv = (dS * SE) ** 2
    tg = obs.index.astype(str).to_numpy()
    kd = 100 * (1 - np.exp2(obs['target_lfc'].to_numpy(float)))
    gp = {g: i for i, g in enumerate(genes)}; own = np.array([gp.get(t, -1) for t in tg])
    measurable = np.array([j >= 0 and measured[j] and C[j] >= 1.0 for j in own]) & np.isfinite(kd)
    save(name, tg, S, Vv, S0, np.zeros(len(genes)), obs['n_cells'].to_numpy(float), np.nan, measured, kd, measurable,
         dict(ctrl_detect=np.full(len(genes), np.nan), P0=C / 1e6))


if __name__ == '__main__':
    sys.path.insert(0, str(V / 'scripts'))
    from build_effects_v4 import panel
    genes = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy()
    var_panel, _ = panel()
    for name in sys.argv[1:]:
        log('building', name)
        if name in CELL: cell_line(name, CELL[name], genes, var_panel)
        elif name in TABLE: table_line(name, TABLE[name], genes)
        else: cd4_line(name, CD4[name], genes)

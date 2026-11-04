"""v6-lin: pooled PCA programs + tercile-calibrated linear transfer from the training lines to H1 (pseudobulk scoring).

  cd ~/vcc && PYTHONPATH=scripts/v5 ~/venv/bin/python analysis/v6lin/v6lin.py

Lines (training and source): the seven Model 2 used (K562 genome-wide, RPE1, HepG2, Jurkat, HCT116, HEK293T,
CD4_Stim48hr as the one T-cell state), effects_v4 tables (lfc = log2 mean-CPM ratio, delta-method se). H1 enters no
fitting step: its table and the 75-target reference are read only for scoring and the diagnostics' truth side.

Definitions: Delta_{L,k} = row lfc; R_L = mean over L's kept knockdowns; d = Delta - R_L. Strength = ||d||^2 - sum se^2
over measured genes (noise-corrected), as a within-line percentile rank; a knockdown's predicted strength is the mean of
its rank over the source lines that measured it; terciles cut at 1/3 and 2/3 of that rank.
Row weights (PCA and every regression): knockdown < 30% dropped, 30-50% and unmeasurable 0.5, else 1 (knockdown % =
100 (1 - 2^efficiency), efficiency = the target's own lfc). PCA only: K562 rows need strength > 0 and n_de_bh > 5, and
K562's total weight is capped at the largest other line's; rows are also precision-weighted (1 / mean se^2 on the
panel, mean 1 within a line).
R_hat_{-L} = mean over L's knockdowns of Model 2's prediction for them (random-effects mean over the other lines,
weight 1 / (se^2 + tau2_g), tau2 DerSimonian-Laird over the seven lines, own gene excluded).
Similarity weights = the v3 transfer model's: exp((corr - max corr) / 0.5), corr = centred correlation of the lines'
baseline log1p CPM (effects_v4 *_baseline_v2.csv).
"""
import json
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy.stats import spearmanr

from pb import H1, row_r

V = Path.home() / 'vcc'
OUT = V / 'analysis/v6lin'
E = V / 'lib/effects_v4'
LINES = {'K562': 'replogle_k562_gw', 'RPE1': 'replogle_rpe1', 'HepG2': 'nadig_hepg2', 'Jurkat': 'nadig_jurkat',
         'HCT116': 'orion_hct116', 'HEK293T': 'orion_hek293t', 'CD4_Stim48hr': 'marson_cd4_stim48hr'}
BASAL = dict(LINES, H1='arc_h1', A='challenge_A', B='challenge_B', C='challenge_C')
KMAX, K_MAIN, CHUNK, TEMP = 100, 30, 128, 0.5
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)


def log(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)


class Line:
    def __init__(self, name, src, genes):
        with h5py.File(E / f'{src}.h5ad', 'r') as f:
            obs = read_elem(f['obs']); self.lfc = f['layers/lfc'][:].astype(np.float16); self.se = f['layers/se'][:].astype(np.float16)
        self.name = name; self.tg = obs['target_gene'].astype(str).to_numpy(); self.n = len(self.tg)
        self.row = {}
        for i, t in enumerate(self.tg): self.row.setdefault(t, i)
        eff = np.where(obs['efficiency_known'].to_numpy(bool), obs['efficiency'].to_numpy(float), np.nan)
        self.kd_pct = 100 * (1 - np.exp2(eff))
        self.w = np.where(~np.isfinite(self.kd_pct), 0.5, np.where(self.kd_pct < 30, 0.0, np.where(self.kd_pct < 50, 0.5, 1.0)))
        self.n_de = obs['n_de_bh'].to_numpy(float)
        gp = {g: i for i, g in enumerate(genes)}; self.own = np.array([gp.get(t, -1) for t in self.tg])
        self.measured = np.zeros(len(genes), bool); s1 = np.zeros(len(genes)); c = np.zeros(len(genes))
        for a in range(0, self.n, 2000):
            L = self.lfc[a:a + 2000].astype(np.float64); ok = np.isfinite(L); self.measured |= ok.any(0)
            k = self.w[a:a + 2000] > 0
            s1 += np.where(ok[k], L[k], 0).sum(0); c += ok[k].sum(0)
        self.R = np.where(c > 0, s1 / np.maximum(c, 1), np.nan)
        st = np.zeros(self.n)
        for a in range(0, self.n, 2000):
            d = self.lfc[a:a + 2000].astype(np.float64) - self.R; s = self.se[a:a + 2000].astype(np.float64)
            ok = np.isfinite(d) & np.isfinite(s)
            st[a:a + 2000] = np.where(ok, d * d - s * s, 0).sum(1)
        self.strength = st; self.pct = pd.Series(st).rank(pct=True).to_numpy()

    def d(self, rows):
        out = np.full((len(rows), self.lfc.shape[1]), np.nan, np.float32); ok = rows >= 0
        out[ok] = self.lfc[rows[ok]].astype(np.float32) - self.R.astype(np.float32); return out

    def raw(self, rows, which):
        arr = self.lfc if which == 'lfc' else self.se
        out = np.full((len(rows), arr.shape[1]), np.nan, np.float32); ok = rows >= 0
        out[ok] = arr[rows[ok]].astype(np.float32); return out


def basal_profile(stem, genes):
    b = pd.read_csv(E / f'{stem}_baseline_v2.csv', index_col=0).reindex(genes)
    return b['log1p_cpm_S'].to_numpy(float)


def centred_similarity(a, b):
    ok = np.isfinite(a) & np.isfinite(b); x, y = a[ok] - a[ok].mean(), b[ok] - b[ok].mean()
    return float(np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y)))


def tau2(lines, genes):
    """DerSimonian-Laird per gene over the seven lines, pooled over targets with >= 2 lines (own gene excluded)."""
    tgs = sorted(set().union(*[set(l.row) for l in lines.values()]))
    gp = {g: i for i, g in enumerate(genes)}; G = len(genes)
    Q, df, C = np.zeros(G), np.zeros(G), np.zeros(G)
    multi = [t for t in tgs if sum(t in l.row for l in lines.values()) >= 2]
    for a in range(0, len(multi), CHUNK):
        ts = multi[a:a + CHUNK]
        y = np.stack([l.raw(np.array([l.row.get(t, -1) for t in ts]), 'lfc') for l in lines.values()]).astype(np.float64)
        v = np.stack([l.raw(np.array([l.row.get(t, -1) for t in ts]), 'se') for l in lines.values()]).astype(np.float64) ** 2
        for i, t in enumerate(ts):
            if t in gp: y[:, i, gp[t]] = np.nan
        ok = np.isfinite(y) & np.isfinite(v) & (v > 0); k = ok.sum(0); use = k >= 2
        w = np.where(ok, 1 / np.where(ok, v, 1), 0); yy = np.where(ok, y, 0); sw = w.sum(0)
        mu = np.divide((w * yy).sum(0), sw, out=np.zeros_like(sw), where=sw > 0)
        Q += np.where(use, (w * (yy - mu) ** 2).sum(0), 0).sum(0); df += np.where(use, k - 1, 0).sum(0)
        C += np.where(use, sw - np.divide((w ** 2).sum(0), sw, out=np.zeros_like(sw), where=sw > 0), 0).sum(0)
    t2 = np.where(C > 0, np.maximum(0, (Q - df) / np.where(C > 0, C, 1)), np.nan)
    log(f'tau2 from {len(multi)} targets measured in >= 2 lines; median tau2 {np.nanmedian(t2):.4f}')
    return t2


def gather(targets, sources, lines, simw, panel, Vk, t2, genes):
    """For each target: similarity-weighted z_bar (Vk' d on the panel), e_bar (off-panel d), mean source strength rank,
    Model 2's random-effects prediction, and whether any source measured it."""
    gp = {g: i for i, g in enumerate(genes)}; G = len(genes); off = ~panel
    n = len(targets); zbar = np.zeros((n, Vk.shape[1])); ebar = np.full((n, G), np.nan, np.float32)
    pct = np.full(n, np.nan); re = np.full((n, G), np.nan, np.float32); has = np.zeros(n, bool)
    for a in range(0, n, CHUNK):
        ts = targets[a:a + CHUNK]; m = len(ts)
        zs, ws, ds, ys, vs, ps = [], [], [], [], [], []
        for s in sources:
            l = lines[s]; rows = np.array([l.row.get(t, -1) for t in ts]); ok = rows >= 0
            d = l.d(rows); zs.append(np.nan_to_num(d[:, panel]) @ Vk); ws.append(np.where(ok, simw[s], 0.0)); ds.append(d)
            ys.append(l.raw(rows, 'lfc')); vs.append(l.raw(rows, 'se')); ps.append(np.where(ok, l.pct[np.maximum(rows, 0)], np.nan))
        W = np.stack(ws); sw = W.sum(0); has[a:a + m] = sw > 0
        zbar[a:a + m] = np.einsum('sn,snj->nj', W, np.stack(zs)) / np.maximum(sw, 1e-12)[:, None]
        D = np.stack(ds); okD = np.isfinite(D); Wg = W[:, :, None] * okD          # all genes: e_bar off-panel, d_bar (1d) everywhere
        ebar[a:a + m] = np.where(Wg.sum(0) > 0, (Wg * np.nan_to_num(D)).sum(0) / np.maximum(Wg.sum(0), 1e-12), np.nan)
        with np.errstate(invalid='ignore'): pct[a:a + m] = np.nanmean(np.stack(ps), 0)
        Y = np.stack(ys).astype(np.float64); Vv = np.stack(vs).astype(np.float64) ** 2 + t2
        for i, t in enumerate(ts):
            if t in gp: Y[:, i, gp[t]] = np.nan
        ok = np.isfinite(Y) & np.isfinite(Vv); w = np.where(ok, 1 / np.where(ok, Vv, 1), 0); sw2 = w.sum(0)
        re[a:a + m] = np.where(sw2 > 0, (w * np.nan_to_num(Y)).sum(0) / np.maximum(sw2, 1e-30), np.nan)
    return dict(zbar=zbar, ebar=ebar, pct=pct, re=re, has=has)


def tercile(p):
    return np.where(np.isfinite(p), np.minimum((np.nan_to_num(p) * 3).astype(int), 2), -1)


def main():
    OUT.mkdir(parents=True, exist_ok=True); res = {}
    genes = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy(); G = len(genes)
    lines = {}
    for k, s in LINES.items(): lines[k] = Line(k, s, genes); log(f'loaded {k}: {lines[k].n} rows')
    h1 = Line('H1', 'arc_h1', genes)                                     # truth only (scoring / diagnostics)
    # ---------------- 1b panel
    panel = np.ones(G, bool)
    for l in lines.values(): panel &= l.measured
    names = list(lines)
    pair = pd.DataFrame([[int((lines[a].measured & lines[b].measured).sum()) for b in names] for a in names], index=names, columns=names)
    pair.to_csv(OUT / 'panel_pairwise.csv')
    res['1b'] = dict(program_panel_size=int(panel.sum()), per_line_measured={k: int(l.measured.sum()) for k, l in lines.items()},
                     model3_panel='7,283 challenge genes with control CPM >= 1.0 in each of HCT116, HEK293T and K562 genome-wide '
                                  '(code/vcc_2026_3c/vcc3_ridge.py program_genes, lines 90-97, GENE_SET_SOURCES line 67, MIN_CPM line 73)')
    log('panel', res['1b']['program_panel_size']); print(pair.to_string())
    # ---------------- strength + 1a
    simcache = {k: basal_profile(BASAL[k], genes) for k in BASAL}
    def simw(target, sources):
        c = np.array([centred_similarity(simcache[target], simcache[s]) for s in sources]); w = np.exp((c - c.max()) / TEMP)
        return dict(zip(sources, w)), dict(zip(sources, c))
    a1 = {}
    for k, l in lines.items():
        others = [s for s in names if s != k]
        mp = np.array([np.nanmean([lines[s].pct[lines[s].row[t]] for s in others if t in lines[s].row]) if any(t in lines[s].row for s in others) else np.nan for t in l.tg])
        ok = np.isfinite(mp); a1[k] = dict(n=int(ok.sum()), spearman=float(spearmanr(l.strength[ok], mp[ok])[0]))
    mpH = np.array([np.nanmean([lines[s].pct[lines[s].row[t]] for s in names if t in lines[s].row]) if any(t in lines[s].row for s in names) else np.nan for t in h1.tg])
    okH = np.isfinite(mpH)
    a1['H1'] = dict(n=int(okH.sum()), spearman=float(spearmanr(h1.strength[okH], mpH[okH])[0]))
    conf = pd.crosstab(pd.Series(tercile(h1.pct), name='real tercile (H1 own rank)'), pd.Series(tercile(mpH), name='predicted tercile'))
    conf.to_csv(OUT / 'h1_tercile_confusion.csv')
    res['1a'] = dict(spearman=a1, confusion=conf.to_dict()); print(pd.DataFrame(a1).T.round(3).to_string()); print(conf.to_string())
    # ---------------- 1c PCA
    rows_used = {}; C = np.zeros((panel.sum(), panel.sum())); tw = {}
    for k, l in lines.items():
        keep = l.w > 0
        if k == 'K562': keep &= (l.strength > 0) & (l.n_de > 5)
        idx = np.flatnonzero(keep)
        se_p = np.array([np.nanmean(l.se[i, panel].astype(np.float64) ** 2) for i in idx])
        prec = 1 / se_p; prec /= np.nanmean(prec); w = l.w[idx] * np.nan_to_num(prec)
        rows_used[k] = (idx, w); tw[k] = float(w.sum())
    cap = max(v for k, v in tw.items() if k != 'K562')
    if tw['K562'] > cap:
        idx, w = rows_used['K562']; rows_used['K562'] = (idx, w * cap / tw['K562'])
    for k, (idx, w) in rows_used.items():
        l = lines[k]
        for a in range(0, len(idx), 1000):
            X = np.nan_to_num(l.d(idx[a:a + 1000])[:, panel]).astype(np.float64) * np.sqrt(w[a:a + 1000])[:, None]
            C += X.T @ X
    ev, U = np.linalg.eigh(C); order = np.argsort(-ev); ev, U = ev[order], U[:, order]; Vfull = U[:, :KMAX]
    res['1c_rows'] = {k: dict(rows_in_table=lines[k].n, rows_used=int(len(idx)), total_weight=float(w.sum()),
                              dropped_kd_lt_30=int((lines[k].w == 0).sum()), half_weight=int((lines[k].w == 0.5).sum()))
                      for k, (idx, w) in rows_used.items()}
    res['1c_rows']['K562']['note'] = f'strength > 0 and n_de_bh > 5; total weight capped at {cap:.1f} (was {tw["K562"]:.1f})'
    print(pd.DataFrame(res['1c_rows']).T.to_string())
    # H1 capture by top K (panel genes): raw and noise-corrected
    dH = h1.d(np.arange(h1.n)).astype(np.float64); seH = h1.se.astype(np.float64)
    okH2 = np.isfinite(dH) & np.isfinite(seH)
    dP = np.where(okH2[:, panel], dH[:, panel], 0); s2P = np.where(okH2[:, panel], seH[:, panel] ** 2, 0)
    tot, tot_nc = (dP ** 2).sum(), (dP ** 2).sum() - s2P.sum()
    cap_tab = {}
    for K in (5, 10, 30, 50, 100):
        Z = dP @ Vfull[:, :K]; noise = (s2P @ (Vfull[:, :K] ** 2)).sum()
        cap_tab[K] = dict(raw_share=float((Z ** 2).sum() / tot), noise_corrected_share=float(((Z ** 2).sum() - noise) / tot_nc))
    res['1c_capture'] = cap_tab; res['1c_pooled_var_share_top'] = {K: float(ev[:K].sum() / ev.sum()) for K in (5, 10, 30, 50, 100)}
    print(pd.DataFrame(cap_tab).T.round(4).to_string())
    # off-panel share of H1's real sum ||d_k||^2 (lfc space, all 300; and the scorer's pseudobulk space, 75 targets)
    full2 = np.where(okH2, dH ** 2, 0); full_nc = np.where(okH2, dH ** 2 - seH ** 2, 0)
    h = H1(); run_set = h.run_set
    dPB = (h.R - h.R.mean(0)) ** 2
    res['1b_offpanel'] = dict(lfc_raw_off_program_panel=float(full2[:, ~panel].sum() / full2.sum()),
                              lfc_noise_corrected_off_program_panel=float(full_nc[:, ~panel].sum() / full_nc.sum()),
                              pseudobulk_off_program_panel=float(dPB[:, ~panel].sum() / dPB.sum()),
                              lfc_raw_off_transformer_5000=float(full2[:, ~run_set].sum() / full2.sum()),
                              pseudobulk_off_transformer_5000=float(dPB[:, ~run_set].sum() / dPB.sum()))
    print(json.dumps(res['1b_offpanel'], indent=1))
    # ---------------- Model 2 tau2, similarity, gather
    t2 = np.nan_to_num(tau2(lines, genes), nan=0.0).astype(np.float64)
    fits = {}
    acc = {K: dict(bzz=np.zeros((4, K)), bz2=np.zeros((4, K)), gde=np.zeros(4), ge2=np.zeros(4), ade=np.zeros(4), ae2=np.zeros(4)) for K in (10, 30, 50)}
    sims = {}
    for L in names:
        srcs = [s for s in names if s != L]; w_s, c_s = simw(L, srcs); sims[L] = c_s
        l = lines[L]; idx = np.flatnonzero(l.w > 0); tg = l.tg[idx]
        g = gather(list(tg), srcs, lines, w_s, panel, Vfull[:, :50], t2, genes)
        t = tercile(g['pct']); ok = g['has']; ww = l.w[idx]
        Rhat = np.nanmean(np.where(ok[:, None], g['re'], np.nan), 0); Rhat = np.nan_to_num(Rhat)
        for a in range(0, len(idx), 1000):
            sl = slice(a, a + 1000); dL = l.d(idx[sl]).astype(np.float64); lfcL = l.raw(idx[sl], 'lfc').astype(np.float64)
            zL = np.nan_to_num(dL[:, panel]) @ Vfull[:, :50]
            e = g['ebar'][sl].astype(np.float64); okg = np.isfinite(dL) & np.isfinite(e) & ~panel[None]
            de = np.where(okg, dL * e, 0).sum(1); e2 = np.where(okg, e * e, 0).sum(1)
            okR = np.isfinite(lfcL); dr = np.where(okR, lfcL * Rhat, 0).sum(1); r2 = np.where(okR, Rhat ** 2, 0).sum(1)
            for K in (10, 30, 50):
                A = acc[K]; zb = g['zbar'][sl][:, :K]; zl = zL[:, :K]
                for tt in (0, 1, 2, 3):
                    s = ok[sl] & ((t[sl] == tt) if tt < 3 else (t[sl] >= 0))
                    w_ = ww[sl][s]
                    A['bzz'][tt] += (w_[:, None] * zl[s] * zb[s]).sum(0); A['bz2'][tt] += (w_[:, None] * zb[s] ** 2).sum(0)
                    A['gde'][tt] += (w_ * de[s]).sum(); A['ge2'][tt] += (w_ * e2[s]).sum()
                    A['ade'][tt] += (w_ * dr[s]).sum(); A['ae2'][tt] += (w_ * r2[s]).sum()
        log(f'fit pass {L}: {len(idx)} rows, {int(ok.sum())} with a source; tercile counts {np.bincount(t[t >= 0], minlength=3).tolist()}')
    for K, A in acc.items():
        fits[K] = dict(beta=A['bzz'] / np.maximum(A['bz2'], 1e-30), gamma=A['gde'] / np.maximum(A['ge2'], 1e-30), a=A['ade'] / np.maximum(A['ae2'], 1e-30))
    labels = ['weak', 'mid', 'strong', 'pooled']
    bt = pd.DataFrame(fits[30]['beta'][:, :10], index=labels, columns=[f'PC{j + 1}' for j in range(10)])
    ga = pd.DataFrame({K: dict(zip(labels, fits[K]['gamma'])) for K in fits}).T; aa = pd.DataFrame({K: dict(zip(labels, fits[K]['a'])) for K in fits}).T
    bt.to_csv(OUT / 'beta_K30_first10.csv'); ga.to_csv(OUT / 'gamma.csv'); aa.to_csv(OUT / 'a.csv')
    res['2'] = dict(beta_K30_first10=bt.to_dict('index'), gamma=ga.to_dict('index'), a=aa.to_dict('index'), similarity_corr=sims)
    print('beta (K=30, first 10 PCs)'); print(bt.round(3).to_string()); print('gamma'); print(ga.round(3).to_string()); print('a'); print(aa.round(3).to_string())
    # ---------------- 3 predict H1
    w_h, c_h = simw('H1', names); sims['H1'] = c_h
    g = gather(list(h1.tg), names, lines, w_h, panel, Vfull[:, :50], t2, genes)
    Rhat = np.nan_to_num(np.nanmean(np.where(g['has'][:, None], g['re'], np.nan), 0))
    t_pred = tercile(g['pct']); t_real = tercile(h1.pct)
    # 1d evidence ceiling: H1's real d_k vs the similarity-weighted mean of the sources' d_k, full gene space
    r1d = np.full(h1.n, np.nan)
    for i in range(h1.n):
        a_, b_ = dH[i], g['ebar'][i].astype(np.float64); ok_ = np.isfinite(a_) & np.isfinite(b_)
        if ok_.sum() > 10 and np.std(b_[ok_]) > 0: r1d[i] = np.corrcoef(a_[ok_], b_[ok_])[0, 1]
    res['1d'] = {lab: dict(n=int(np.isfinite(r1d[t_real == q]).sum()), median_r=float(np.nanmedian(r1d[t_real == q])))
                 for q, lab in ((0, 'weak'), (1, 'mid'), (2, 'strong'))}
    res['1d']['all'] = dict(n=int(np.isfinite(r1d).sum()), median_r=float(np.nanmedian(r1d)))
    print('1d', json.dumps(res['1d']))
    def predict(K, tsel, pooled=False):
        f = fits[K]; D = np.zeros((h1.n, G), np.float32)
        for i in range(h1.n):
            if not g['has'][i] or tsel[i] < 0:
                D[i] = f['a'][3] * Rhat; continue
            tt = 3 if pooled else tsel[i]
            D[i] = f['a'][tt] * Rhat
            D[i, panel] += (g['zbar'][i, :K] * f['beta'][tt]) @ Vfull[:, :K].T
            D[i, ~panel] += f['gamma'][tt] * np.nan_to_num(g['ebar'][i, ~panel])
        return D
    preds = {'v6-lin K=30 (main)': predict(30, t_pred), 'v6-lin K=10': predict(10, t_pred), 'v6-lin K=50': predict(50, t_pred),
             'v6-lin pooled beta/gamma/a': predict(30, t_pred, pooled=True), 'v6-lin real H1 terciles (oracle)': predict(30, t_real),
             'Model 2 on effects_v4 (7 lines)': np.nan_to_num(g['re'])}
    tg_ch = pd.read_csv(V / 'data/pert_counts.csv')['target_gene'].astype(str).tolist()
    res['3_counts'] = dict(H1_with_source=int(g['has'].sum()), H1_total=h1.n,
                           challenge_with_source_7_lines=int(sum(any(t in lines[s].row for s in names) for t in tg_ch)),
                           challenge_with_source_7_lines_plus_H1=int(sum(any(t in lines[s].row for s in names) or t in h1.row for t in tg_ch)),
                           challenge_total=len(tg_ch))
    print(json.dumps(res['3_counts']))
    # ---------------- 4 score (pseudobulk)
    pos = {t: i for i, t in enumerate(h1.tg)}; order = np.array([pos[t] for t in h.t300])
    rows = []
    def score(name, D300):
        P = h.bulk(D300); r = {'arm': name}
        for q, lab in [(None, 'all'), (0, 'weak'), (1, 'mid'), (2, 'strong')]:
            sel = None if q is None else np.flatnonzero(h.terc == q); s = h.score(P, sel)
            for m in ('mse_raw', 'mse_scaled', 'mse_unclamped', 'pds_raw', 'pds_scaled'): r[f'{lab}:{m}'] = s[m]
        d_real = h.R - h.R.mean(0); r['median_r_demeaned'] = float(np.median(row_r(P - P.mean(0), d_real, h.pan)))
        rows.append(r); log(f"{name}: mse {r['all:mse_raw']:.4f} ({r['all:mse_unclamped']:+.4f}) pds {r['all:pds_raw']:.4f}")
    for k, D in preds.items(): score(k, D[order])
    score('dev transformer (dev_nopie)', h.full(h.delta))
    df = pd.DataFrame(rows).set_index('arm'); df.to_csv(OUT / 'scores_pseudobulk.csv')
    res['4'] = df.to_dict('index')
    print(df[[c for c in df.columns if c.startswith('all:')] + ['median_r_demeaned']].round(4).to_string())
    print(df[[c for c in df.columns if c.endswith('mse_raw') or c.endswith('pds_raw')]].round(4).to_string())
    # save the main prediction as a generator folder (H1 order of the dev predictions)
    pdir = OUT / 'h1_pred'; pdir.mkdir(exist_ok=True)
    src = V / 'lib/vcc_2026_4/runs/dev_nopie/h1_pred_g3'
    for f in ('genes.npy', 'targets.npy', 'ontarget_lfc.npy'): (pdir / f).write_bytes((src / f).read_bytes())
    np.save(pdir / 'v6lin.npy', preds['v6-lin K=30 (main)'][order].astype(np.float32))
    (pdir / 'prediction_report.json').write_text(json.dumps(dict(schema='vcc_2026_3_ridge_v1', heldout_background='H1', source_backgrounds=names,
                                                                  models=['v6lin'], note='v6-lin main (K=30), analysis/v6lin/v6lin.py'), indent=1))
    (OUT / 'results.json').write_text(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()

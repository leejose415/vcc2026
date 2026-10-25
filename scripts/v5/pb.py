"""Pseudobulk H1 harness (the scorer's mse / pds formulas on noise-free pseudobulks; same real side as the cell harness:
lib/vcc_2026_3/val_final, 75 targets) and the v5 parts that run on it: A1 (pb view), A2, B, C, D.

  PYTHONPATH=~/v5 ~/venv/bin/python ~/v5/pb.py a1|a2|b|c|d [--best JSON]

Predicted side = delta -> the generator's calibrate (alpha 1, gate, own gene 0, +/-8 clip; floor only where stated) ->
project_factors on the generator-half H1 controls with the measured knockdown residual. The line mean dbar is over all
300 H1 knockdowns (as for A-C); the 75 scored targets are a subset. Nothing here reads H1 truth to choose genes,
strength groups or weights: the only H1-truth reads are the real pseudobulks used for scoring and the tercile labels
used for reporting.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy.stats import spearmanr

from v5lib import (EFFECTS, H1_DIR, H1_SRC_DIR, RUN, SOURCE_FILES, V, arm, full_array, panel_genes,
                   source_line_means)

sys.path.insert(0, str(V / 'code/vcc_2026_3c'))
from vcc3_generate import (Cfg, calibrate, h1_controls, kd_residual, prediction_ontarget, project_factors,  # noqa: E402
                           reference_targets)

REF = V / 'lib/vcc_2026_3/val_final'
CTRL = 'non-targeting'
BTS = 50000.0
CACHE = Path.home() / 'v5/cache'
ARMS_CSV = V / 'results/tuning/v5_arms.csv'
H1_SOURCES = ['K562', 'RPE1', 'HepG2', 'Jurkat', 'HCT116', 'HEK293T', 'CD4_Rest', 'CD4_Stim8hr', 'CD4_Stim48hr']
COLS = ['part', 'arm', 'level', 's', 'beta', 'ext', 'detail', 'mse_raw', 'mse_scaled', 'mse_unclamped', 'pds_raw',
        'pds_scaled', 'slope_shared_all', 'slope_shared_run', 'nmae', 'fid', 'reach', 'jac', 'overall', 'note']


def log_arm(row):
    ARMS_CSV.parent.mkdir(parents=True, exist_ok=True); new = not ARMS_CSV.exists()
    with open(ARMS_CSV, 'a', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction='ignore')
        if new: w.writeheader()
        w.writerow({k: (round(v, 5) if isinstance(v, float) else v) for k, v in row.items()})


class H1:
    def __init__(self):
        from cell_eval2.score import _replicate_entries
        from cell_eval2.scoring import score_one
        self._score_one = score_one
        genes = panel_genes(); self.genes = genes; G = self.G = len(genes); gp = {g: i for i, g in enumerate(genes)}
        mz = np.load(next((REF / 'real_cache').glob('*.moments.npz')))
        self.perts, self.means, jk = mz['perts'].astype(str), mz['means'].astype(np.float64), mz['jk'].astype(np.float64)
        if (REF / 'real.h5ad').exists(): targets, _ = reference_targets(REF / 'real.h5ad')
        else: targets = sorted(set(self.perts) - {CTRL})          # vcc2026 copy: the reference's targets from its cache
        self.targets = targets
        ridx = {p: i for i, p in enumerate(self.perts)}
        self.R = self.means[[ridx[t] for t in targets]]; self.Creal = jk[[ridx[t] for t in targets]]
        self.ctrl = self.means[ridx[CTRL]]; Cctrl = jk[ridx[CTRL]]
        self.own = np.array([gp.get(t, -1) for t in targets]); self.nrow = np.where(self.own >= 0, G - 1, G).astype(float)
        T = self.T = len(targets); self.rows = np.arange(T)
        dr = self.R - self.ctrl[None]
        self.den = (np.einsum('ij,ij->i', dr, dr) - self.Creal - Cctrl - self.own_sq(self.R, np.broadcast_to(self.ctrl, self.R.shape))) / self.nrow
        anchor = pd.read_parquet(REF / 'bundle/anchor_agg.parquet')
        base = pd.read_csv(REF / 'bundle/baseline_agg.csv').set_index('statistic').loc['mean'].to_dict()
        self.entries = _replicate_entries(base, anchor); rep = anchor.set_index('metric')['replicate'].to_dict()
        self.MSE, self.PDS = 'expr_mse_unbiased_capped_norm', 'pds_cosine'
        self.bm, self.rm = float(base[self.MSE]), float(rep[self.MSE])
        h = ad.read_h5ad(V / 'lib/h1/H1_full_2025.h5ad', backed='r')
        try: self.line = h1_controls(h, h.obs['target_gene'].astype(str).values, list(genes), 2026, 'guide_id')[0]
        finally: h.file.close()
        sel = json.loads((RUN / 'h1_bench_g3/selected_generator.json').read_text())['settings']
        self.cfg = Cfg(argparse.Namespace(**dict(sel, detect_floor=False, alpha=1.0, call='z', z=0.0, min_effect=0.0)))
        self.cfg_floor = Cfg(argparse.Namespace(**dict(sel, detect_floor=True, alpha=1.0, call='z', z=0.0, min_effect=0.0)))
        self.kd = kd_residual(prediction_ontarget(H1_SRC_DIR, targets), self.cfg)
        z = np.load(H1_SRC_DIR / 'raw_predictions.npz')
        self.delta, self.p, self.gcols = z['delta'].astype(np.float64), z['p'].astype(np.float64), z['gene_cols']
        self.t300 = np.load(H1_SRC_DIR / 'targets.npy').astype(str); pos = {t: i for i, t in enumerate(self.t300)}
        self.idx = np.array([pos[t] for t in targets])
        self.run_set = np.zeros(G, bool); self.run_set[self.gcols] = True
        self.pan = np.ones(G, bool); self.pan[[gp[t] for t in targets if t in gp]] = False     # pds feature space
        d = self.R - self.R.mean(0)
        self.size = np.sqrt(np.maximum((d ** 2).sum(1) - self.Creal, 0))                    # real ||d_k||, noise-corrected
        self.terc = np.argsort(np.argsort(self.size)) * 3 // T                                 # reporting only

    def own_sq(self, A, B):
        return np.array([(A[i, self.own[i]] - B[i, self.own[i]]) ** 2 if self.own[i] >= 0 else 0.0 for i in self.rows])

    def bulk(self, D300, floor=False):
        """(300, G) full-panel predictions -> (75, G) noise-free predicted pseudobulks."""
        Dc = calibrate(np.asarray(D300, np.float32)[self.idx], self.line, self.targets, self.cfg_floor if floor else self.cfg)
        P = np.empty((self.T, self.G))
        for i in self.rows:
            _, proj, _ = project_factors(self.line.share, Dc[i], self.own[i] if self.own[i] >= 0 else None, self.kd[i])
            P[i] = np.log1p(BTS * proj)
        return P

    def scale(self, metric, u):
        e = self.entries[metric]; return float(self._score_one(float(u), e.base, e.scoring))

    def pds_rows(self, P):
        from cell_eval2.metrics.discrimination import discrimination_score
        out = discrimination_score(pred_bulk=(np.array(self.targets), P), real_bulk=(self.perts, self.means), distance='cosine',
                                   rank_denominator='n-1', tie_policy='midrank', exclude_target_gene=True,
                                   exclusion_scope='panel', control_source='real', genes=self.genes, control=CTRL)
        return np.array([out[t] for t in self.targets])

    def score(self, P, sel=None):
        sel = self.rows if sel is None else sel
        d = P - self.R
        num = ((np.einsum('ij,ij->i', d, d) - self.Creal - self.own_sq(P, self.R)) / self.nrow)[sel]
        u = num.sum() / self.den[sel].sum(); pv = float(self.pds_rows(P)[sel].mean())
        mh = P.mean(0); m = self.R.mean(0); c = self.ctrl
        sl = {}
        for nm, msk in (('all', np.ones(self.G, bool)), ('run', self.run_set)):
            a, b = (mh - c)[msk], (m - c)[msk]; sl[nm] = float(a @ b / (b @ b))
        return dict(mse_raw=float(u), mse_scaled=self.scale(self.MSE, u), mse_unclamped=(self.bm - u) / (self.bm - self.rm),
                    pds_raw=pv, pds_scaled=self.scale(self.PDS, pv), slope_shared_all=sl['all'], slope_shared_run=sl['run'])

    def full(self, run_vals, ext=None):
        return full_array(run_vals, self.gcols, self.G, ext)


def ext_vector(sources, beta, tag):
    """beta x the per-gene mean over source lines of each line's mean knockdown lfc (>= 3 lines, else 0). Cached."""
    f = CACHE / f'line_means_{tag}.npz'
    if not f.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        m, P, cnt = source_line_means(sources, panel_genes()); np.savez(f, mean=m, per=P, cnt=cnt, sources=np.array(sources))
    return beta * np.load(f)['mean']


def report(df, cols=('mse_raw', 'mse_scaled', 'mse_unclamped', 'pds_raw', 'pds_scaled', 'slope_shared_all', 'slope_shared_run')):
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
    print(df[[c for c in cols if c in df]].round(4).to_string(), flush=True)


def cmd_a1(h):
    rows = {}
    for s in (1, 1.5, 2):
        D = h.full(arm(h.delta, s=s))
        for floor in ((True, False) if s == 1 else (False,)):
            r = h.score(h.bulk(D, floor)); lab = f's{s:g}' + ('_floor' if floor else ''); rows[lab] = r
            log_arm(dict(part='A1', arm=lab, level='pseudobulk', s=s, beta=1, ext='no', detail='floor' if floor else '', **r))
    report(pd.DataFrame(rows).T)


def cmd_a2(h, s, floor=False):
    ext1 = ext_vector(H1_SOURCES, 1.0, 'h1sources'); rows = {}
    for ext in (False, True):
        for beta in (1, 1.5, 2, 2.5, 3, 4):
            D = h.full(arm(h.delta, s=s, beta=beta), beta * ext1 if ext else None)
            r = h.score(h.bulk(D, floor)); lab = f's{s:g}_b{beta:g}' + ('_ext' if ext else '') + ('_floor' if floor else ''); rows[lab] = r
            log_arm(dict(part='A2', arm=lab, level='pseudobulk', s=s, beta=beta, ext='yes' if ext else 'no', detail='floor' if floor else '', **r))
    report(pd.DataFrame(rows).T)
    print('genes outside the run set with an extension value:', int(((ext1 != 0) & ~h.run_set).sum()))


# ------------------------------------------------------------------ source effects for the H1 targets (B, C)
def source_rows(targets, sources):
    """{source: (lfc (K, G) NaN where unmeasured or no row, se, efficiency (K,), has_row (K,))} for the given target genes."""
    out = {}; G = len(panel_genes())
    for src in sources:
        f = CACHE / f'rows_{src}_{len(targets)}.npz'
        if not f.exists():
            with h5py.File(EFFECTS / f'{SOURCE_FILES[src]}.h5ad', 'r') as fh:
                obs = read_elem(fh['obs']); tg = obs['target_gene'].astype(str).to_numpy()
                first = {}
                for i, t in enumerate(tg): first.setdefault(t, i)       # one row per target gene (first if repeated)
                ri = np.array([first.get(t, -1) for t in targets]); has = ri >= 0
                L = np.full((len(targets), G), np.nan, np.float32); S = np.full_like(L, np.nan)
                if has.any():
                    order = np.argsort(ri[has]); rr = ri[has][order]
                    Lh = fh['layers/lfc'][rr]; Sh = fh['layers/se'][rr]
                    tmpL = np.empty_like(Lh); tmpS = np.empty_like(Sh); tmpL[order] = Lh; tmpS[order] = Sh
                    L[has] = tmpL; S[has] = tmpS
                eff = np.full(len(targets), np.nan)
                ek = obs['efficiency_known'].to_numpy().astype(bool); ev = obs['efficiency'].to_numpy().astype(float)
                eff[has] = np.where(ek[ri[has]], ev[ri[has]], np.nan)
            CACHE.mkdir(parents=True, exist_ok=True); np.savez(f, L=L, S=S, eff=eff, has=has)
        z = np.load(f); out[src] = (z['L'], z['S'], z['eff'], z['has'])
        print(f'  {src}: {int(out[src][3].sum())}/{len(targets)} targets', flush=True)
    return out


def row_r(A, B, msk):
    out = np.full(len(A), np.nan)
    for i in range(len(A)):
        a, b = A[i, msk], B[i, msk]
        if np.std(a) > 0 and np.std(b) > 0: out[i] = np.corrcoef(a, b)[0, 1]
    return out


def cmd_b(h):
    src = source_rows(list(h.targets), H1_SOURCES)
    ext_vector(H1_SOURCES, 1.0, 'h1sources')                            # builds the line-mean cache if missing
    lm = np.load(CACHE / 'line_means_h1sources.npz'); per = {s: lm['per'][i] for i, s in enumerate(lm['sources'])}
    T, G = h.T, h.G
    dem = {}; se = {}
    for s, (L, S, _, has) in src.items():
        dem[s] = L - np.nan_to_num(per[s])[None, :]                    # de-meaned by the source line's own mean
        se[s] = S
    # (a) precision-weighted mean of the de-meaned source effects
    num = np.zeros((T, G)); wsum = np.zeros((T, G))
    for s in src:
        ok = np.isfinite(dem[s]) & np.isfinite(se[s]) & (se[s] > 0); w = np.where(ok, 1 / np.where(ok, se[s], 1) ** 2, 0)
        num += w * np.nan_to_num(dem[s]); wsum += w
    a = np.where(wsum > 0, num / np.maximum(wsum, 1e-300), 0.0)
    # (b) per knockdown the source with the lowest median SE (precision), not chosen on H1
    medse = np.array([[np.nanmedian(se[s][i]) if src[s][3][i] else np.inf for s in src] for i in range(T)])
    names = list(src); b = np.zeros((T, G)); bsrc = []
    for i in range(T):
        j = int(np.argmin(medse[i])); bsrc.append(names[j] if np.isfinite(medse[i, j]) else 'none')
        if np.isfinite(medse[i, j]): b[i] = np.nan_to_num(dem[names[j]][i])
    # real d_k and the model
    d = h.R - h.R.mean(0)
    dh_lfc = h.full(arm(h.delta, s=1) - h.delta.mean(0))[h.idx]                    # delta_k - dbar (dbar over 300)
    P = h.bulk(h.full(arm(h.delta, s=1))); dh_pb = P - P.mean(0)
    spaces = {'run5000_expressed_pds': h.run_set & h.line.gate & h.pan, 'all_pds_genes': h.pan}
    res = {}; rows = []
    for sp, msk in spaces.items():
        r = dict(a=row_r(a, d, msk), b=row_r(b, d, msk), d_lfc=row_r(dh_lfc, d, msk), d_pseudobulk=row_r(dh_pb, d, msk))
        cands = np.stack([row_r(np.nan_to_num(dem[s]), d, msk) for s in names])          # (n_src, T)
        cands = np.where(np.stack([src[s][3] for s in names]), cands, np.nan)
        r['c_oracle'] = np.nanmax(np.where(np.isfinite(cands), cands, -np.inf), 0); r['c_oracle'][~np.isfinite(r['c_oracle'])] = np.nan
        for k, v in r.items():
            for q, name in [(None, 'all'), (0, 'weak'), (1, 'mid'), (2, 'strong')]:
                x = v if q is None else v[h.terc == q]; x = x[np.isfinite(x)]
                rows.append(dict(space=sp, predictor=k, tercile=name, n=len(x), median=float(np.median(x)) if len(x) else np.nan,
                                 q25=float(np.percentile(x, 25)) if len(x) else np.nan, q75=float(np.percentile(x, 75)) if len(x) else np.nan))
    df = pd.DataFrame(rows); out = V / 'results/tuning/v5_partB.csv'; df.to_csv(out, index=False)
    piv = df.pivot_table(index=['space', 'predictor'], columns='tercile', values='median')[['all', 'weak', 'mid', 'strong']]
    print(piv.round(3).to_string())
    print(df[df.tercile == 'all'].round(3).to_string(index=False))
    print('targets with >= 1 source row:', int(sum(any(src[s][3][i] for s in names) for i in range(T))), 'of', T,
          '| (b) chosen sources:', pd.Series(bsrc).value_counts().to_dict())


# ------------------------------------------------------------------ C: strength
def strengths(h):
    """Predicted strength for all 300 H1 knockdowns, without H1 data."""
    src = source_rows(list(h.t300), H1_SOURCES)
    st = {}
    st['i_sum_pde'] = h.p.sum(1)
    st['ii_norm_specific'] = np.linalg.norm(h.delta - h.delta.mean(0), axis=1)
    effs = np.stack([src[s][2] for s in src]); st['iii_median_ontarget'] = -np.nanmedian(effs, 0)   # stronger = larger
    rms = []; wts = []
    for s, (L, S, _, has) in src.items():
        with np.errstate(invalid='ignore'):
            r = np.sqrt(np.maximum(np.nanmean(L.astype(np.float64) ** 2 - S.astype(np.float64) ** 2, 1), 0))
            w = 1 / np.nanmedian(S.astype(np.float64) ** 2, 1)
        rms.append(np.where(has, r, np.nan)); wts.append(np.where(has, w, 0))
    rms, wts = np.stack(rms), np.stack(wts)
    with np.errstate(invalid='ignore'):
        st['iv_pw_source_effect'] = np.nansum(np.nan_to_num(rms) * wts, 0) / np.where(wts.sum(0) > 0, wts.sum(0), np.nan)
    return pd.DataFrame(st, index=h.t300)


def cmd_c(h, best):
    st = strengths(h); st.to_csv(V / 'results/tuning/v5_partC_strength.csv')
    cor = {}
    for k in st:
        x = st[k].to_numpy()[h.idx]; ok = np.isfinite(x)
        cor[k] = dict(spearman=float(spearmanr(x[ok], h.size[ok])[0]), n=int(ok.sum()), missing_of_300=int((~np.isfinite(st[k])).sum()))
    print(pd.DataFrame(cor).T.round(3).to_string())
    pick = max(cor, key=lambda k: cor[k]['spearman'])
    json.dump(dict(correlations=cor, chosen=pick), open(V / 'results/tuning/v5_partC_choice.json', 'w'), indent=1)
    if cor[pick]['spearman'] < 0.3:
        print(f'all below 0.3 (best {pick} {cor[pick]["spearman"]:.3f}): C2 skipped'); return
    x = st[pick].to_numpy().copy(); x[~np.isfinite(x)] = np.nanmedian(x)            # missing -> middle third
    third = np.argsort(np.argsort(x)) * 3 // len(x)                                     # within-line thirds, 300 knockdowns
    beta, ext = best['beta'], best['ext']; ext1 = ext_vector(H1_SOURCES, 1.0, 'h1sources') if ext else None
    grid = [0.25, 0.5, 1, 1.5, 2, 3]; rows = {}
    for lo in grid:
        for mid in grid:
            for hi in grid:
                s = np.array([lo, mid, hi])[third]
                r = h.score(h.bulk(h.full(arm(h.delta, s=s, beta=beta), beta * ext1 if ext else None), best.get('floor', False)))
                rows[(lo, mid, hi)] = r
    df = pd.DataFrame(rows).T; df.index.names = ['s_lo', 's_mid', 's_hi']; df['sum_scaled'] = df.mse_scaled + df.pds_scaled
    df.to_csv(V / 'results/tuning/v5_partC_grid.csv')
    picks = {'best_mse_unclamped': df.mse_unclamped.idxmax(), 'best_pds': df.pds_raw.idxmax(), 'best_sum_scaled': df.sum_scaled.idxmax(),
             'uniform_s': (best['s'],) * 3}
    for nm, key in picks.items():
        r = df.loc[key].to_dict()
        log_arm(dict(part='C2', arm=f'{nm} lo/mid/hi={key}', level='pseudobulk', s=str(key), beta=beta, ext='yes' if ext else 'no',
                     detail=f'thirds by {pick}' + (' floor' if best.get('floor') else ''), **{k: v for k, v in r.items() if k != 'sum_scaled'}))
    report(df.loc[list(picks.values())].assign(pick=list(picks)).set_index('pick', append=True),
           cols=('mse_raw', 'mse_scaled', 'mse_unclamped', 'pds_raw', 'pds_scaled', 'sum_scaled'))


def cmd_d(h, best):
    sv = np.load(Path.home() / 'h1_decomp/train_signal_var.npy')      # training lines only (six target lines)
    v = sv[h.gcols]; rk = pd.Series(v).rank(method='average').to_numpy().copy(); rk[~np.isfinite(v)] = 1
    r = (rk - np.nanmin(rk)) / (np.nanmax(rk) - np.nanmin(rk))
    beta, ext = best['beta'], best['ext']; ext1 = ext_vector(H1_SOURCES, 1.0, 'h1sources') if ext else None
    rows = {}
    for s0 in (1, 1.5, 2):
        for al in (0, 0.5, 1, 2):
            res = h.score(h.bulk(h.full(arm(h.delta, s=s0, beta=beta, per_gene=1 + al * r), beta * ext1 if ext else None), best.get('floor', False)))
            rows[(s0, al)] = res
            log_arm(dict(part='D', arm=f's0={s0:g} alpha={al:g}', level='pseudobulk', s=s0, beta=beta, ext='yes' if ext else 'no',
                         detail='per-gene (1 + alpha r_g)' + (' floor' if best.get('floor') else ''), **res))
    df = pd.DataFrame(rows).T; df.index.names = ['s0', 'alpha']; df['sum_scaled'] = df.mse_scaled + df.pds_scaled
    df.to_csv(V / 'results/tuning/v5_partD_grid.csv')
    report(df, cols=('mse_raw', 'mse_scaled', 'mse_unclamped', 'pds_raw', 'pds_scaled', 'sum_scaled'))
    print('best by sum_scaled:', df.sum_scaled.idxmax())


def cmd_e1(h, tag):
    """Per tercile of real ||d_k|| (H1, 75): per-knockdown slope of d_k on d_hat_k, Pearson r, ||d_hat_k|| / ||d_k||
    (pds gene space, noise-corrected ||d_k||), plus the pooled slope. d_hat = the s = 1 pseudobulk minus its mean."""
    P = h.bulk(h.full(arm(h.delta, s=1))); dh = (P - P.mean(0))[:, h.pan]; d = (h.R - h.R.mean(0))[:, h.pan]
    slope = (dh * d).sum(1) / (dh * dh).sum(1); r = row_r(dh, d, np.ones(dh.shape[1], bool))
    ratio = np.linalg.norm(dh, axis=1) / np.maximum(h.size, 1e-12)
    rows = []
    for q, name in [(None, 'all'), (0, 'weak'), (1, 'mid'), (2, 'strong')]:
        sel = h.rows if q is None else np.flatnonzero(h.terc == q)
        rows.append(dict(run=tag, tercile=name, n=len(sel), median_slope=float(np.median(slope[sel])),
                         pooled_slope=float((dh[sel] * d[sel]).sum() / (dh[sel] * dh[sel]).sum()),
                         median_r=float(np.median(r[sel])), r_q25=float(np.percentile(r[sel], 25)), r_q75=float(np.percentile(r[sel], 75)),
                         median_norm_ratio=float(np.median(ratio[sel]))))
    df = pd.DataFrame(rows); f = V / 'results/tuning/v5_E1_h1_terciles.csv'
    df.to_csv(f, mode='a', header=not f.exists(), index=False); print(df.round(4).to_string(index=False))


def cmd_ens(h):
    """Item 3: TTA (mean over no-mask + each source masked), checkpoint average (best + last), Model 2 blends.
    Each variant -> pseudobulk; arm D s grid in pseudobulk space (m_hat + s (P - m_hat)); per-knockdown median r of
    d_hat with real d_k (pds gene space)."""
    from vcc3_generate import load_deltas
    z = np.load(Path.home() / 'v5/ens_preds.npz'); assert np.array_equal(z['gene_cols'], h.gcols)
    masks = [k for k in z.files if k.startswith('best_mask_')]
    run = {'transformer (best.pt)': z['best_none'], 'TTA mean (none + %d masked)' % len(masks): np.mean([z['best_none']] + [z[k] for k in masks], 0),
           'last.pt': z['last_none'], 'checkpoint mean (best + last)': (z['best_none'] + z['last_none']) / 2}
    for k in masks: run[f'mask {k[10:]} only'] = z[k]
    variants = {k: h.full(v.astype(np.float64)) for k, v in run.items()}
    m2 = load_deltas(V / 'lib/vcc_2026_3/validation/pred/H1', list(h.t300), list(h.genes), 're')
    for w in (0, 0.25, 0.5, 0.75):
        variants[f'blend w={w:g} (w transformer + (1-w) Model 2)'] = w * variants['transformer (best.pt)'] + (1 - w) * m2
    d = (h.R - h.R.mean(0)); rows = []
    for name, D in variants.items():
        P = h.bulk(D); mh = P.mean(0); dh = P - mh
        r = row_r(dh, d, h.pan)
        for s in (0.5, 1, 1.5, 2, 3):
            res = h.score(mh + s * dh)
            rows.append(dict(variant=name, s=s, median_r=float(np.median(r)), r_q25=float(np.percentile(r, 25)),
                             r_q75=float(np.percentile(r, 75)), **res))
            log_arm(dict(part='3 ensembling', arm=name, level='pseudobulk', s=s, beta=1, ext='no', detail='arm D: m_hat + s d_hat',
                         note=f'median r {np.median(r):.4f}', **res))
    df = pd.DataFrame(rows); df.to_csv(V / 'results/tuning/v5_item3_ensembling.csv', index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
    print(df[['variant', 's', 'median_r', 'r_q25', 'r_q75', 'mse_raw', 'mse_scaled', 'mse_unclamped', 'pds_raw', 'pds_scaled']].round(4).to_string(index=False))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('cmd', choices=['a1', 'a2', 'b', 'c', 'd', 'e1', 'ens'])
    ap.add_argument('--s', type=float, nargs='+', default=[1.0], help='a2: the A1 winner s (several: one grid each)')
    ap.add_argument('--pred', help='raw_predictions.npz of another validation run (E3): replaces delta and p')
    ap.add_argument('--tag', default='dev_nopie')
    ap.add_argument('--floor', action='store_true', help='a2: detection floor (the A1 winner has it)')
    ap.add_argument('--best', help='c / d: JSON {"s":..,"beta":..,"ext":true|false} of the Part A best')
    a = ap.parse_args(); h = H1()
    if a.pred:
        z = np.load(a.pred); assert np.array_equal(z['gene_cols'], h.gcols)
        h.delta, h.p = z['delta'].astype(np.float64), z['p'].astype(np.float64)
    if a.cmd == 'a1': cmd_a1(h)
    elif a.cmd == 'a2':
        for s in a.s: cmd_a2(h, s, a.floor)
    elif a.cmd == 'b': cmd_b(h)
    elif a.cmd == 'e1': cmd_e1(h, a.tag)
    elif a.cmd == 'ens': cmd_ens(h)
    else:
        best = json.loads(a.best); {'c': cmd_c, 'd': cmd_d}[a.cmd](h, best)

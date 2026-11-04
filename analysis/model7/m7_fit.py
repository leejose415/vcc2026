"""Model 7, Steps 1-3: diagnostics, ridge fit in the scorer's space, predictions for H1 (validation) and A-C (final).

  cd ~/vcc && PYTHONPATH=scripts/v5:scripts ~/venv/bin/python analysis/model7/m7_fit.py diag|validate|final

Tables: lib/vcc_2026_5/m7/<line>.npz (m7_tables.py). Model for destination D, knockdown k, gene g in group q:
  Delta_hat = a_q R_hat_{-D,g} + sum_{L != D} beta_{L,q} d_{L,k,g}   (no similarity weighting anywhere)
R_hat_{-D} = SE-weighted mean of the other lines' R (weight 1 / var(R_L,g), var(R_L,g) = sum_k var / n^2).
Gene groups: deciles of D's control expression (scorer space) over D's measured genes; 1-4 pooled, then 5..10 (7).
Fit: pooled over destinations and their kept knockdowns (knockdown weights 1 / 0.5; no precision weights), own gene
left out; per group a ridge with the penalty chosen by GCV. Variants: one beta per source (groups pooled), and the
program denoiser (Denoiser) applied to every source's d before refitting.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
T7 = V / 'lib/vcc_2026_5/m7'
OUT = V / 'analysis/model7'
TRAIN = ['K562', 'RPE1', 'HepG2', 'Jurkat', 'HCT116', 'HEK293T', 'CD4_Rest', 'CD4_Stim8hr', 'CD4_Stim48hr']
BASAL = {'K562': 'replogle_k562_gw', 'RPE1': 'replogle_rpe1', 'HepG2': 'nadig_hepg2', 'Jurkat': 'nadig_jurkat', 'HCT116': 'orion_hct116',
         'HEK293T': 'orion_hek293t', 'CD4_Rest': 'marson_cd4_rest', 'CD4_Stim8hr': 'marson_cd4_stim8hr', 'CD4_Stim48hr': 'marson_cd4_stim48hr',
         'H1': 'arc_h1', 'A': 'challenge_A', 'B': 'challenge_B', 'C': 'challenge_C'}
LAMBDAS = [-2, 0, 1, 2, 4, 8, 16]
ALPHAS = np.logspace(-4, 4, 33)
CHUNK = 96
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)


def log(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)


class Tab:
    def __init__(self, name):
        z = np.load(T7 / f'{name}.npz', allow_pickle=False)
        self.name = name; self.tg = z['tg'].astype(str); self.D = z['D']; self.vD = z['vD']; self.var0 = z['var0'].astype(np.float32)
        self.S0 = z['S0']; self.w = z['w']; self.R = z['R']; self.measured = z['measured']; self.n = self.D.shape[0]
        self.row = {}
        for i, t in enumerate(self.tg):
            if self.w[i] > 0: self.row.setdefault(t, i)
        k = self.w > 0
        vk = np.where(np.isfinite(self.vD[k]), self.vD[k].astype(np.float64), 0)
        self.varR = vk.sum(0) / max(k.sum(), 1) ** 2

    def d(self, rows):
        out = np.full((len(rows), self.D.shape[1]), np.nan, np.float32); ok = rows >= 0
        out[ok] = self.D[rows[ok]].astype(np.float32) - self.R.astype(np.float32); return out


def genes_list():
    return pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy()


def basal(stem, genes):
    return pd.read_csv(V / f'lib/effects_v4/{stem}_baseline_v2.csv', index_col=0).reindex(genes)['log1p_cpm_S'].to_numpy(float)


def centred_similarity(a, b):
    ok = np.isfinite(a) & np.isfinite(b); x, y = a[ok] - a[ok].mean(), b[ok] - b[ok].mean()
    return float(np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y)))


def r_hat(tabs, names):
    num = np.zeros(len(tabs[names[0]].R)); den = np.zeros_like(num)
    for n in names:
        t = tabs[n]; ok = np.isfinite(t.R) & (t.varR > 0)
        w = np.where(ok, 1 / np.where(ok, t.varR, 1), 0); num += w * np.nan_to_num(t.R); den += w
    return np.where(den > 0, num / np.maximum(den, 1e-300), np.nan)


def groups_of(S0, measured):
    q = np.zeros(len(S0), int)
    ok = measured & np.isfinite(S0)
    rk = pd.Series(S0[ok]).rank(pct=True, method='first').to_numpy()
    dec = np.minimum((rk * 10).astype(int), 9)                 # 0..9
    q[ok] = np.maximum(dec - 3, 0)                              # deciles 1-4 -> 0, 5..10 -> 1..6
    return q


def rows_r(A, B, M):
    out = np.full(len(A), np.nan)
    for i in range(len(A)):
        m = M[i]
        if m.sum() > 10:
            a, b = A[i, m], B[i, m]
            if a.std() > 0 and b.std() > 0: out[i] = np.corrcoef(a, b)[0, 1]
    return out


# ------------------------------------------------------------------ Step 1 diagnostics
def pair(S, Dt, genes):
    shared = [t for t in S.row if t in Dt.row]
    if len(shared) < 10: return None
    iS = np.array([S.row[t] for t in shared]); iD = np.array([Dt.row[t] for t in shared])
    gp = {g: i for i, g in enumerate(genes)}
    m = S.measured & Dt.measured
    DS = S.D[iS].astype(np.float32); DD = Dt.D[iD].astype(np.float32); vD = Dt.vD[iD].astype(np.float32)
    M = np.isfinite(DS) & np.isfinite(DD) & m[None]
    for r, t in enumerate(shared):
        if t in gp: M[r, gp[t]] = False
    vk = np.where(M, vD - Dt.var0[None], 0); vall = np.where(M, vD, 0)
    err = np.where(M, (DS - DD) ** 2, 0).sum(1)
    num = err - vk.sum(1); den = np.where(M, DD ** 2, 0).sum(1) - vall.sum(1)
    dS = DS - S.R[None].astype(np.float32); dD = DD - Dt.R[None].astype(np.float32)
    return dict(source=S.name, dest=Dt.name, n_shared=len(shared), median_r_Delta=float(np.nanmedian(rows_r(DS, DD, M))),
                median_r_d=float(np.nanmedian(rows_r(dS, dD, M))), raw_mse=float(num.sum() / den.sum()))


def shared_response(tabs, Dn, others, genes):
    Dt = tabs[Dn]; k = np.flatnonzero(Dt.w > 0); gp = {g: i for i, g in enumerate(genes)}
    Rh = r_hat(tabs, others); out = {}
    for lab, pred in (('own R_D (oracle)', Dt.R), ('R_hat_-D (SE-weighted)', Rh)):
        num = den = 0.0; dot = rr = 0.0
        for a in range(0, len(k), 500):
            rr_ = k[a:a + 500]; DD = Dt.D[rr_].astype(np.float64); vD = Dt.vD[rr_].astype(np.float64)
            M = np.isfinite(DD) & Dt.measured[None]
            for r, i in enumerate(rr_):
                j = gp.get(Dt.tg[i], -1)
                if j >= 0: M[r, j] = False
            P = np.where(np.isfinite(pred), pred, 0)[None]
            num += (np.where(M, (P - DD) ** 2 - (vD - Dt.var0[None]), 0)).sum(); den += (np.where(M, DD ** 2 - vD, 0)).sum()
            dot += np.where(M, P * DD, 0).sum(); rr += np.where(M, P * P, 0).sum()
        out[lab] = num / den
        if lab.startswith('R_hat'):
            a_ = dot / rr; num2 = den2 = 0.0
            for a in range(0, len(k), 500):
                rr_ = k[a:a + 500]; DD = Dt.D[rr_].astype(np.float64); vD = Dt.vD[rr_].astype(np.float64)
                M = np.isfinite(DD) & Dt.measured[None]
                for r, i in enumerate(rr_):
                    j = gp.get(Dt.tg[i], -1)
                    if j >= 0: M[r, j] = False
                P = a_ * np.where(np.isfinite(pred), pred, 0)[None]
                num2 += (np.where(M, (P - DD) ** 2 - (vD - Dt.var0[None]), 0)).sum(); den2 += (np.where(M, DD ** 2 - vD, 0)).sum()
            out['a * R_hat_-D'] = num2 / den2; out['a'] = a_
    return out


def diag():
    genes = genes_list(); names = TRAIN + ['H1']; tabs = {n: Tab(n) for n in names}; res = {}
    res['step0'] = {n: json.loads((T7 / f'{n}_report.json').read_text()) for n in names}
    print(pd.DataFrame(res['step0']).T.to_string())
    rows = []
    for Dn in names:
        if Dn != 'CD4_Stim48hr':
            r = pair(tabs['CD4_Stim48hr'], tabs[Dn], genes)
            if r: rows.append(r)
    for a, b in (('CD4_Rest', 'CD4_Stim8hr'), ('CD4_Stim8hr', 'CD4_Stim48hr'), ('CD4_Rest', 'CD4_Stim48hr'),
                 ('CD4_Stim8hr', 'CD4_Rest'), ('CD4_Stim48hr', 'CD4_Stim8hr'), ('CD4_Stim48hr', 'CD4_Rest')):
        r = pair(tabs[a], tabs[b], genes)
        if r: rows.append(r)
    df = pd.DataFrame(rows); df.to_csv(OUT / 'step1a_cd4.csv', index=False); res['1a'] = df.to_dict('records'); print(df.round(3).to_string(index=False))
    sr = {}
    for Dn in names:
        others = [n for n in TRAIN if n != Dn]; sr[Dn] = shared_response(tabs, Dn, others, genes); log(Dn, sr[Dn])
    res['1b'] = sr; print(pd.DataFrame(sr).T.round(4).to_string())
    (OUT / 'step1.json').write_text(json.dumps(res, indent=1, default=float))


# ------------------------------------------------------------------ Step 2 fit
class Denoiser:
    """Programs: PCA of precision-weighted d over RPE1, HepG2 and Jurkat (kept knockdowns, knockdown weights x 1 / mean
    variance on the program genes, mean 1 within a line), on the genes all three measure; top K components V.
    d' = d + shrink (V V' d - d) on those genes, unchanged elsewhere."""
    def __init__(self, tabs, lines=('RPE1', 'HepG2', 'Jurkat'), K=30, shrink=0.5):
        self.genes = np.ones(len(tabs[lines[0]].R), bool)
        for l in lines: self.genes &= tabs[l].measured
        gi = np.flatnonzero(self.genes); C = np.zeros((len(gi), len(gi)))
        for l in lines:
            t = tabs[l]; k = np.flatnonzero(t.w > 0)
            prec = 1 / np.maximum(np.nanmean(t.vD[k][:, gi].astype(np.float64), 1), 1e-12); prec /= prec.mean()   # mean 1 within the line
            for a in range(0, len(k), 500):
                rr = k[a:a + 500]; d = np.nan_to_num(t.d(rr)[:, gi]).astype(np.float64)
                ww = t.w[rr] * prec[a:a + 500]
                C += (d * ww[:, None]).T @ d
        self.lines = list(lines)
        ev, U = np.linalg.eigh(C); self.V = U[:, np.argsort(-ev)[:K]].astype(np.float32); self.gi = gi; self.shrink = shrink

    def apply(self, d):
        x = np.nan_to_num(d[:, self.gi]); proj = (x @ self.V) @ self.V.T
        out = d.copy(); fin = np.isfinite(d[:, self.gi])
        out[:, self.gi] = np.where(fin, d[:, self.gi] + self.shrink * (proj - x), d[:, self.gi]); return out

    def share(self, t):
        k = np.flatnonzero(t.w > 0); num = den = 0.0
        for a in range(0, len(k), 500):
            x = np.nan_to_num(t.d(k[a:a + 500])[:, self.gi]).astype(np.float64)
            num += ((x @ self.V.astype(np.float64)) ** 2).sum(); den += (x ** 2).sum()
        return num / den


def accumulate(tabs, dests, sources_all, genes, den=None):
    """Per gene group: Gram of [R_hat_-D, d_L for L in sources_all] (0 where L = D or L lacks the knockdown), X'y,
    y'y, n; rows = each destination's kept knockdowns x measured genes (own gene out), knockdown weights."""
    nf = 1 + len(sources_all); Q = 7
    Gm = np.zeros((Q, nf, nf)); Xy = np.zeros((Q, nf)); yy = np.zeros(Q); nn = np.zeros(Q)
    gp = {g: i for i, g in enumerate(genes)}; G = len(genes)
    for Dn in dests:
        Dt = tabs[Dn]; srcs = [s for s in sources_all if s != Dn and s in tabs]
        Rh = np.nan_to_num(r_hat(tabs, srcs)).astype(np.float32)
        q = groups_of(Dt.S0, Dt.measured)
        k = np.flatnonzero(Dt.w > 0); t0 = time.time()
        for a in range(0, len(k), CHUNK):
            rr = k[a:a + CHUNK]; m = len(rr); y = Dt.D[rr].astype(np.float32)
            M = np.isfinite(y) & Dt.measured[None]
            for r, i in enumerate(rr):
                j = gp.get(Dt.tg[i], -1)
                if j >= 0: M[r, j] = False
            sw = np.sqrt(Dt.w[rr]).astype(np.float32)[:, None]
            F = np.zeros((nf, m, G), np.float32); F[0] = Rh[None]
            for si, s in enumerate(sources_all):
                if s not in srcs: continue
                rows = np.array([tabs[s].row.get(t, -1) for t in Dt.tg[rr]])
                if (rows < 0).all(): continue
                dl = tabs[s].d(rows)
                if den is not None: dl = den.apply(dl)
                F[1 + si] = np.nan_to_num(dl)
            y = np.where(M, y, 0) * sw; F = F * sw[None] * M[None]
            for qq in range(Q):
                cols = q == qq
                Fq = F[:, :, cols].reshape(nf, -1).astype(np.float64); yq = y[:, cols].ravel().astype(np.float64)
                Gm[qq] += Fq @ Fq.T; Xy[qq] += Fq @ yq; yy[qq] += yq @ yq; nn[qq] += M[:, cols].sum()
        log(f'  accumulated {Dn}: {len(k)} knockdowns ({time.time() - t0:.0f}s)')
    return dict(G=Gm, Xy=Xy, yy=yy, n=nn)


def solve(acc, pooled=False):
    """Ridge per gene group (or one pooled fit), penalty by GCV."""
    nf = acc['G'].shape[1]
    groups = [list(range(7))] if pooled else [[q] for q in range(7)]
    coef, info = np.zeros((7, nf)), []
    for gs in groups:
        G = sum(acc['G'][q] for q in gs); Xy = sum(acc['Xy'][q] for q in gs); yy = sum(acc['yy'][q] for q in gs); n = sum(acc['n'][q] for q in gs)
        ev, U = np.linalg.eigh(G); best = None
        for al in ALPHAS * np.trace(G) / nf:
            b = U @ ((U.T @ Xy) / (ev + al)); rss = yy - 2 * b @ Xy + b @ G @ b; df = (ev / (ev + al)).sum()
            gcv = n * rss / (n - df) ** 2
            if best is None or gcv < best[0]: best = (gcv, al, b, rss)
        for q in gs: coef[q] = best[2]
        info.append(dict(groups=gs, alpha=float(best[1]), n=float(n), r2=float(1 - best[3] / yy)))
    return coef, info


def predict(tabs, Dn, tg, sources, coef, genes, S0, measured_D, den=None):
    """Delta_hat (len(tg) x G) for destination Dn's knockdowns tg, every panel gene; genes no source measures get a_q R_hat."""
    G = len(genes); srcs = [s for s in sources if s in tabs and s != Dn]
    Rh = np.nan_to_num(r_hat(tabs, srcs)); q = groups_of(S0, measured_D)
    out = (coef[q, 0] * Rh)[None].repeat(len(tg), 0).astype(np.float32)
    for si, s in enumerate(sources):
        if s not in srcs: continue
        rows = np.array([tabs[s].row.get(t, -1) for t in tg])
        if (rows < 0).all(): continue
        dl = tabs[s].d(rows)
        if den is not None: dl = den.apply(dl)
        out += coef[q, 1 + si][None] * np.nan_to_num(dl)
    return out


QLAB = ['dec1-4', 'dec5', 'dec6', 'dec7', 'dec8', 'dec9', 'dec10']


def fit_all(tabs, dests, sources, genes, tag, den=None):
    acc = accumulate(tabs, dests, sources, genes, den); np.savez(OUT / f'acc_{tag}.npz', **acc)
    fits = {'main': solve(acc), 'single_beta': solve(acc, pooled=True)}
    tabc = pd.DataFrame(fits['main'][0], index=QLAB, columns=['a'] + [f'beta_{s}' for s in sources])
    log(tag); print(tabc.round(4).to_string()); print(pd.DataFrame(fits['main'][1]).to_string())
    print('single beta:', dict(zip(['a'] + list(sources), np.round(fits['single_beta'][0][0], 4))))
    return fits, tabc


def validate():
    genes = genes_list(); names = TRAIN + ['H1']; tabs = {n: Tab(n) for n in names}; tr = {n: tabs[n] for n in TRAIN}
    fits, tabc = fit_all(tr, TRAIN, TRAIN, genes, 'validation'); tabc.to_csv(OUT / 'step2_coefficients_validation.csv')
    den = Denoiser(tr)
    shares = {n: dict(program_share_of_sum_d2=float(den.share(tr[n])),
                      noise_share_step0=json.loads((T7 / f'{n}_report.json').read_text())['noise_share_of_sum_Delta2']) for n in TRAIN}
    print(pd.DataFrame(shares).T.round(4).to_string())
    fits_dn, tabc_dn = fit_all(tr, TRAIN, TRAIN, genes, 'validation_denoised', den); tabc_dn.to_csv(OUT / 'step2_coefficients_validation_denoised.csv')
    h1 = tabs['H1']; tg = h1.tg
    arrays = {'model7': predict(tr, 'H1', tg, TRAIN, fits['main'][0], genes, h1.S0, h1.measured),
              'model7_single_beta': predict(tr, 'H1', tg, TRAIN, fits['single_beta'][0], genes, h1.S0, h1.measured),
              'model7_denoised': predict(tr, 'H1', tg, TRAIN, fits_dn['main'][0], genes, h1.S0, h1.measured, den)}
    np.savez(OUT / 'h1_predictions_scorer_space.npz', tg=tg, S0=h1.S0, **arrays)
    (OUT / 'step2_validation.json').write_text(json.dumps(dict(
        fits={k: dict(coef=v[0].tolist(), info=v[1]) for k, v in fits.items()},
        fits_denoised={k: dict(coef=v[0].tolist(), info=v[1]) for k, v in fits_dn.items()},
        denoiser=dict(genes=int(den.genes.sum()), K=30, shrink=0.5, lines=den.lines), program_shares=shares), indent=1, default=float))


def final():
    """Refit with H1 as source and destination; predict A-C (300 targets, every panel gene)."""
    import anndata as ad
    genes = genes_list(); names = TRAIN + ['H1']; tabs = {n: Tab(n) for n in names}
    fits, tabc = fit_all(tabs, names, names, genes, 'final'); tabc.to_csv(OUT / 'step2_coefficients_final.csv')
    den = Denoiser(tabs); fits_dn, tabc_dn = fit_all(tabs, names, names, genes, 'final_denoised', den)
    tabc_dn.to_csv(OUT / 'step2_coefficients_final_denoised.csv')
    tg = pd.read_csv(V / 'data/pert_counts.csv')['target_gene'].astype(str).to_numpy(); out = {}
    for c in 'ABC':
        a = ad.read_h5ad(V / f'data/context_{c}.h5ad'); assert list(a.var_names.astype(str)) == list(genes)
        tot = np.asarray(a.X.sum(0)).ravel().astype(np.float64); S0 = np.log1p(5e4 * tot / tot.sum()); allg = np.ones(len(genes), bool)
        out[f'{c}_model7'] = predict(tabs, c, tg, names, fits['main'][0], genes, S0, allg)
        out[f'{c}_model7_single_beta'] = predict(tabs, c, tg, names, fits['single_beta'][0], genes, S0, allg)
        out[f'{c}_model7_denoised'] = predict(tabs, c, tg, names, fits_dn['main'][0], genes, S0, allg, den)
        out[f'{c}_S0'] = S0; log(f'predicted {c}')
    np.savez(OUT / 'abc_predictions_scorer_space.npz', tg=tg, **out)
    (OUT / 'step2_final.json').write_text(json.dumps(dict(
        fits={k: dict(coef=v[0].tolist(), info=v[1]) for k, v in fits.items()},
        fits_denoised={k: dict(coef=v[0].tolist(), info=v[1]) for k, v in fits_dn.items()}), indent=1, default=float))


if __name__ == '__main__':
    {'diag': diag, 'validate': validate, 'final': final}[sys.argv[1]]()

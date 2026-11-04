"""Part 2 (addendum) on the pseudobulk H1 harness: partial projection onto training-line programs, the amplitude hybrid,
and the H1R oracle. All in LFC space on the run's 5,000 genes (genes outside it unchanged), then the generator path.

  cd ~/v5 && PYTHONPATH=~/v5 ~/venv/bin/python part2.py

W (programs): the v3 loadings (lib/vcc_2026_3/models/ridge_validation, H1 held out) cover 7,283 program genes and stop
at K = 20, so W is rebuilt here from the six training lines only, on the run genes: per line, lfc minus the line mean
(NaN -> 0), covariance X'X / n_line summed over lines (each line weighted equally), top 50 eigenvectors.
Partial projection: d' = d + lam (P_K d - d), d = delta - dbar; delta' = dbar + s d'.
Amplitude hybrid: R = dbar (control = 0 in lfc space), d_perp = d - (d.R / R.R) R,
pred_k = dbar + (c_hat_k - mean c_hat) R + lam P_K d_perp_k. c_k = (t_k . R_line) / (R_line . R_line) in a line;
c_hat from a ridge on [mean c_k over the other source lines, DepMap mean gene effect, median on-target lfc over the other
source lines, mean noise-corrected effect RMS over the other source lines], trained on the six training lines with each
line's features built leaving that line out; H1 features use the nine non-H1 sources. H1 c_k (truth) only for reporting.
H1R (oracle, upper bound only): dbar := H1's real mean lfc response.
"""
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from pb import CACHE, H1, H1_SOURCES, V, log_arm, row_r
from v5lib import EFFECTS, SOURCE_FILES

TRAIN = ['K562', 'RPE1', 'HepG2', 'Jurkat', 'HCT116', 'HEK293T']
LAMS = [0, 0.25, 0.5, 0.75, 1]
KS = [5, 10, 20, 50]


def load_line(src, gcols):
    """Run-gene lfc / se of every knockdown of a source line, its target genes and on-target efficiency."""
    f = CACHE / f'line_{src}_run.npz'
    if not f.exists():
        with h5py.File(EFFECTS / f'{SOURCE_FILES[src]}.h5ad', 'r') as fh:
            obs = read_elem(fh['obs']); L = fh['layers/lfc'][:, :][:, gcols]; S = fh['layers/se'][:, :][:, gcols]
        eff = np.where(obs['efficiency_known'].to_numpy().astype(bool), obs['efficiency'].to_numpy().astype(float), np.nan)
        np.savez(f, L=L.astype(np.float32), S=S.astype(np.float32), tg=np.array(obs['target_gene'].astype(str).tolist(), dtype=str), eff=eff)
    z = np.load(f, allow_pickle=False); return z['L'].astype(np.float64), z['S'].astype(np.float64), z['tg'].astype(str), z['eff']


def amplitude(L):
    """c_k = (t_k . R) / (R . R) over the row's measured genes, R = the line mean."""
    M = np.isfinite(L); R = np.where(M, L, 0).sum(0) / np.maximum(M.sum(0), 1)
    T = np.where(M, L, 0); RM = np.where(M, R[None, :], 0)
    return (T * RM).sum(1) / np.maximum((RM * RM).sum(1), 1e-12), R


def effect_rms(L, S):
    with np.errstate(invalid='ignore'):
        return np.sqrt(np.maximum(np.nanmean(L ** 2 - S ** 2, 1), 0))


def main():
    h = H1(); gc = h.gcols; res = {}
    lines = {s: load_line(s, gc) for s in H1_SOURCES}
    # ---- W from the training lines
    C = np.zeros((len(gc), len(gc)))
    for s in TRAIN:
        L = lines[s][0]; M = np.isfinite(L); mean = np.where(M, L, 0).sum(0) / np.maximum(M.sum(0), 1)
        X = np.where(M, L - mean, 0); C += X.T @ X / len(X)
    w, U = np.linalg.eigh(C); U = U[:, ::-1][:, :max(KS)]; w = w[::-1]
    res['W'] = dict(built='training lines only (six), run genes', var_share={K: float(w[:K].sum() / w.sum()) for K in KS})
    print('W variance share', res['W']['var_share'], flush=True)

    def proj(X, K): Uk = U[:, :K]; return X @ Uk @ Uk.T

    # ---- per-source amplitude, efficiency, effect size per target gene
    per = {}
    for s, (L, S, tg, eff) in lines.items():
        c, _ = amplitude(L); rms = effect_rms(L, S)
        d = pd.DataFrame(dict(tg=tg, c=c, eff=eff, rms=rms)).groupby('tg').first(); per[s] = d
    dep = np.load(V / 'lib/pie_features/depmap_gene_effect.npy').astype(np.float64)
    cov = np.load(V / 'lib/pie_features/coverage_mask.npy')
    genes = list(h.genes); gp = {g: i for i, g in enumerate(genes)}
    dcol = json.loads((V / 'lib/pie_features/sources.json').read_text())['columns'].index('depmap_gene_effect')
    dep_has = cov[:, dcol].astype(bool)
    dep_mean = np.where(dep_has, dep.mean(1), np.nan)

    def feats(targets, exclude):
        rows = []
        for t in targets:
            others = [s for s in H1_SOURCES if s not in exclude and t in per[s].index]
            cs = [per[s].loc[t, 'c'] for s in others]; es = [per[s].loc[t, 'eff'] for s in others]; rs = [per[s].loc[t, 'rms'] for s in others]
            rows.append([np.nanmean(cs) if cs else np.nan, dep_mean[gp[t]] if t in gp else np.nan,
                         np.nanmedian(es) if np.isfinite(es).any() else np.nan, np.nanmean(rs) if rs else np.nan])
        return np.array(rows, float)

    Xtr, ytr = [], []
    for s in TRAIN:
        d = per[s]; X = feats(list(d.index), {s}); Xtr.append(X); ytr.append(d['c'].to_numpy())
    Xtr, ytr = np.vstack(Xtr), np.concatenate(ytr)
    fill = np.nanmean(Xtr, 0)
    def prep(X): return np.hstack([np.where(np.isfinite(X), X, fill), np.isfinite(X).astype(float)])
    ok = np.isfinite(ytr)
    ridge = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-3, 3, 13))).fit(prep(Xtr[ok]), ytr[ok])
    Xh = feats(list(h.t300), set()); chat = ridge.predict(prep(Xh))
    with h5py.File(EFFECTS / 'arc_h1.h5ad', 'r') as fh:
        obs = read_elem(fh['obs']); LH = fh['layers/lfc'][:]
    ch1, Rh1_full = amplitude(LH); h1c = pd.Series(ch1, index=obs['target_gene'].astype(str)).groupby(level=0).first()
    real_c = h1c.reindex(h.t300).to_numpy()
    rho = spearmanr(chat, real_c, nan_policy='omit')[0]
    rho75 = spearmanr(chat[h.idx], real_c[h.idx], nan_policy='omit')[0]
    res['c_hat'] = dict(train_rows=int(ok.sum()), ridge_alpha=float(ridge[-1].alpha_), coef=ridge[-1].coef_.round(4).tolist(),
                        features='[mean c other lines, DepMap mean gene effect, median on-target lfc other lines, mean effect RMS other lines] + 4 has-flags',
                        spearman_chat_vs_H1_c_300=float(rho), spearman_chat_vs_H1_c_75=float(rho75),
                        chat_sd=float(np.std(chat)), H1_c_sd=float(np.nanstd(real_c)))
    print(json.dumps(res['c_hat'], indent=1), flush=True)

    dbar = h.delta.mean(0); d = h.delta - dbar
    real_d = h.R - h.R.mean(0); out = []

    def run(label, part, run_vals, ext=None, **kw):
        P = h.bulk(h.full(run_vals, ext)); dh = P - P.mean(0)
        r = float(np.median(row_r(dh, real_d, h.pan))); sc = h.score(P)
        rec = dict(part=part, arm=label, median_r=r, **sc, **kw); out.append(rec)
        log_arm(dict(part=f'2 {part}', arm=label, level='pseudobulk', s=kw.get('s', ''), beta=1, ext='H1 real mean' if 'H1R' in part else 'no',
                     detail=json.dumps(kw), note=f'median r {r:.4f}', **sc))
        print(f"{part:12s} {label:34s} r {r:.4f} mse {sc['mse_raw']:.4f} ({sc['mse_unclamped']:+.4f}) pds {sc['pds_raw']:.4f} ({sc['pds_scaled']:.4f})", flush=True)

    # ---- partial projection
    for s in (1, 1.5, 2):
        run(f's={s:g} lam=0', 'projection', dbar + s * d, s=s, lam=0, K='-')
        for K in KS:
            Pd = proj(d, K)
            for lam in LAMS[1:]:
                run(f's={s:g} K={K} lam={lam:g}', 'projection', dbar + s * (d + lam * (Pd - d)), s=s, lam=lam, K=K)
    # ---- amplitude hybrid
    R = dbar; dperp = d - np.outer(d @ R / (R @ R), R); cc = (chat - chat.mean())[:, None] * R[None, :]
    for K in KS + ['all']:
        Pd = dperp if K == 'all' else proj(dperp, K)
        for lam in LAMS:
            if lam == 0 and K != KS[0]: continue
            run(f'K={K} lam={lam:g}', 'hybrid', dbar + cc + lam * Pd, lam=lam, K=K)
    # ---- H1R oracle: dbar := H1's real mean lfc
    Rh1 = np.nan_to_num(Rh1_full)
    for s in (1, 1.5, 2):
        run(f'H1R s={s:g}', 'H1R oracle', Rh1[gc] + s * d, ext=Rh1, s=s)
    Rr = Rh1[gc]; dperp_r = d - np.outer(d @ Rr / (Rr @ Rr), Rr); ccr = (chat - chat.mean())[:, None] * Rr[None, :]
    for lam in (0, 0.5, 1):
        run(f'H1R hybrid K=all lam={lam:g}', 'H1R oracle', Rr + ccr + lam * dperp_r, ext=Rh1, lam=lam, K='all')
    df = pd.DataFrame(out); df.to_csv(V / 'results/tuning/v5_part2.csv', index=False)
    (V / 'results/tuning/v5_part2.json').write_text(json.dumps(res, indent=1, default=float))
    for part, g in df.groupby('part', sort=False):
        g = g.assign(sum_scaled=g.mse_scaled + g.pds_scaled)
        print(f'\n== {part}: best by mse_unclamped / pds_raw / sum_scaled')
        print(g.loc[[g.mse_unclamped.idxmax(), g.pds_raw.idxmax(), g.sum_scaled.idxmax()],
                    ['arm', 'median_r', 'mse_raw', 'mse_unclamped', 'pds_raw', 'pds_scaled', 'sum_scaled']].round(4).to_string(index=False))


if __name__ == '__main__':
    main()

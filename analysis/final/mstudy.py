"""M-study (M1-M3) and U-arm (U1-U3), scorer space (lib/vcc_2026_5/m7 tables: summed raw counts, jackknife SEs; CD4
via Marson limma SEs mapped, approximate).

  cd ~/vcc && PYTHONPATH=scripts/v5:analysis/model7:analysis/final ~/venv/bin/python analysis/final/mstudy.py

Sources: K562, RPE1, HepG2, Jurkat, HCT116, HEK293T, CD4_Stim48hr (one T-cell vote). Destinations: each source line
(the other six as sources) and H1 (all seven). QC weights from the tables (rows with weight 0 dropped).
M1 signature S_{D,k}: genes where >= n_min sources (not D) have |Delta / SE| >= z_min with the same sign.
M2 on S: sign agreement of D's real Delta with the consensus sign (null: the signature of a random other knockdown of
D, same genes' consensus signs); no-intercept slope of D's real Delta on the sources' SE-weighted mean Delta; share of
D's sum ||Delta_k||^2 inside S. By tercile of D's real strength (noise-corrected ||Delta_k||^2).
M3 (H1): Delta_hat = k_q x (SE-weighted mean) on S, 0 elsewhere; k_q on half 1, scored on half 2, swapped (seed 0
halves, as calib.py); alone and added to F1's calibrated prediction (scorer space).
U1: per line U_L = mean Delta over knockdowns of DepMap-essential targets (gene effect < -0.5 in >= 50% of lines)
minus mean over the other targets; pairwise Pearson, and against H1's own U. U2: U = SE-weighted mean of the training
lines' U_L; e_k = share of DepMap lines with gene effect < -0.5. U3 (H1): k (e_k - mean e) U, k on half 1, swapped.
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from calib import groups
from m7_fit import Tab
from m7_score import to_log2_factor  # noqa: F401  (A-C conversion is done in the build scripts)
from pb import H1, row_r

V = Path.home() / 'vcc'
F = V / 'analysis/final'
SRC = ['K562', 'RPE1', 'HepG2', 'Jurkat', 'HCT116', 'HEK293T', 'CD4_Stim48hr']
GRID = [(3, 2), (3, 3), (5, 2), (5, 3)]
CH = 64


def log(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)


def gather(tabs, srcs, tg, gp):
    D = np.stack([tabs[s].D[np.maximum(np.array([tabs[s].row.get(t, -1) for t in tg]), 0)].astype(np.float32) for s in srcs])
    Vv = np.stack([tabs[s].vD[np.maximum(np.array([tabs[s].row.get(t, -1) for t in tg]), 0)].astype(np.float32) for s in srcs])
    has = np.stack([np.array([t in tabs[s].row for t in tg]) for s in srcs])
    ok = has[:, :, None] & np.isfinite(D) & np.isfinite(Vv) & (Vv > 0)
    for i, t in enumerate(tg):
        if t in gp: ok[:, i, gp[t]] = False
    return np.where(ok, D, 0), np.where(ok, Vv, 1), ok


def signature(D, Vv, ok, z, n):
    zz = np.where(ok, D / np.sqrt(Vv), 0); pos = (zz >= z).sum(0); neg = (zz <= -z).sum(0)
    sgn = np.where(pos >= n, 1, np.where(neg >= n, -1, 0)); sgn[(pos >= n) & (neg >= n)] = 0
    w = np.where(ok, 1 / Vv, 0); mean = (w * D).sum(0) / np.maximum(w.sum(0), 1e-30)
    return sgn, mean


def m12(tabs, dest, srcs, gp, rng):
    Dt = tabs[dest]; k = np.flatnonzero(Dt.w > 0); tg = Dt.tg[k]
    keep = np.array([any(t in tabs[s].row for s in srcs) for t in tg]); k, tg = k[keep], tg[keep]
    real_all = Dt.D[k].astype(np.float32); var_all = Dt.vD[k].astype(np.float32)
    strength = np.nansum(real_all.astype(np.float64) ** 2 - var_all, 1); terc = np.argsort(np.argsort(strength)) * 3 // max(len(k), 1)
    out = {g: dict(S=[], agree=[], agree_null=[], num=[], den=[], inS=[], tot=[], terc=[]) for g in GRID}
    perm = rng.permutation(len(k))
    for a in range(0, len(k), CH):
        sl = slice(a, a + CH); D, Vv, ok = gather(tabs, srcs, tg[sl], gp)
        real = real_all[sl]; fin = np.isfinite(real) & Dt.measured[None]; r0 = np.where(fin, real, 0)
        for g in GRID:
            sgn, mean = signature(D, Vv, ok, *g); S = (sgn != 0) & fin
            o = out[g]; o['S'] += S.sum(1).tolist(); o['terc'] += terc[sl].tolist()
            ag = np.where(S, np.sign(r0) == sgn, False).sum(1); o['agree'] += [ag[i] / max(S[i].sum(), 1) if S[i].sum() else np.nan for i in range(len(ag))]
            o['num'] += (np.where(S, r0 * mean, 0)).sum(1).tolist(); o['den'] += (np.where(S, mean * mean, 0)).sum(1).tolist()
            o['inS'] += np.where(S, r0 ** 2, 0).sum(1).tolist(); o['tot'] += np.where(fin, r0 ** 2, 0).sum(1).tolist()
        # null: each knockdown scored against another knockdown's signature (permuted within the destination)
        pk = perm[sl]; Dn, Vn, okn = gather(tabs, srcs, tg[pk], gp)
        for g in GRID:
            sgn, _ = signature(Dn, Vn, okn, *g); S = (sgn != 0) & fin
            ag = np.where(S, np.sign(r0) == sgn, False).sum(1)
            out[g]['agree_null'] += [ag[i] / max(S[i].sum(), 1) if S[i].sum() else np.nan for i in range(len(ag))]
    rows = []
    for g, o in out.items():
        df = pd.DataFrame({k_: np.asarray(v, float) for k_, v in o.items()})
        for tq, lab in [(None, 'all'), (0, 'weak'), (1, 'mid'), (2, 'strong')]:
            s = df if tq is None else df[df.terc == tq]
            rows.append(dict(dest=dest, z_min=g[0], n_min=g[1], tercile=lab, knockdowns=len(s), median_S=float(s.S.median()),
                             with_S=int((s.S >= 1).sum()), sign_agree=float(np.nanmean(s.agree)), sign_agree_null=float(np.nanmean(s.agree_null)),
                             slope=float(s.num.sum() / max(s.den.sum(), 1e-12)), share_in_S=float(s.inS.sum() / max(s.tot.sum(), 1e-12))))
    log(f'M2 {dest}: {len(k)} knockdowns')
    return rows


def main():
    genes = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy(); gp = {g: i for i, g in enumerate(genes)}
    tabs = {n: Tab(n) for n in SRC + ['H1']}; rng = np.random.default_rng(0); res = {}
    rows = []
    for dest in SRC + ['H1']:
        srcs = [s for s in SRC if s != dest]; rows += m12(tabs, dest, srcs, gp, rng)
    m2 = pd.DataFrame(rows); m2.to_csv(F / 'm12.csv', index=False)
    print(m2[m2.tercile == 'all'].round(3).to_string(index=False))
    # ---------------- M3 (H1)
    h = H1(); h1 = tabs['H1']; tg = list(h1.tg); pos = {t: i for i, t in enumerate(tg)}
    n = len(tg); half = np.zeros(n, int); half[np.random.default_rng(0).permutation(n)[: n // 2]] = 1
    # calib.py's halves are over Model 2 v4's target order (arc_h1 effects table); map them by target
    z = np.load(F / 'm2v4_h1.npz', allow_pickle=False); ztg = z['targets'].astype(str)
    hz = np.zeros(len(ztg), int); hz[np.random.default_rng(0).permutation(len(ztg))[: len(ztg) // 2]] = 1
    half = np.array([hz[list(ztg).index(t)] if t in set(ztg) else 0 for t in tg])
    real = h1.D.astype(np.float32); fin = np.isfinite(real) & h1.measured[None]
    q = groups(h1.S0, h1.measured)
    D, Vv, ok = gather(tabs, SRC, h1.tg, gp)
    idx75 = np.array([pos[t] for t in h.targets]); cal = F / 'h1_pred/f1_cal.npy'
    f1cal = None
    if cal.exists():                                   # F1 calibrated (log2 factors on the generator share) back to scorer space
        C = h.line.share[None].astype(np.float64); fz = np.load(cal).astype(np.float64)
        fs = np.log1p(5e4 * C * np.exp2(fz)) - np.log1p(5e4 * C)
        zpos = {t: i for i, t in enumerate(ztg)}
        f1cal = np.stack([fs[zpos[t]] if t in zpos else np.zeros(len(genes)) for t in tg]).astype(np.float32)
    m3 = {}
    for g in GRID:
        sgn, mean = signature(D, Vv, ok, *g); S = (sgn != 0)
        pred = np.where(S, mean, 0)
        kq = {}
        for hh in (0, 1):
            m = (half == hh)[:, None] & S & fin
            kq[hh] = np.array([float((np.where(m & (q == qq)[None], real * pred, 0)).sum() / max((np.where(m & (q == qq)[None], pred ** 2, 0)).sum(), 1e-12)) for qq in range(7)])
        K = np.stack([kq[1 - half[i]] for i in range(n)])
        sparse = K[np.arange(n)[:, None], q[None, :]] * pred
        for lab, A in (('alone', sparse), ('+ F1 calibrated', None if f1cal is None else f1cal + sparse)):
            if A is None: continue
            P = np.maximum(h.ctrl[None] + A[idx75], 0); s = h.score(P)
            m3[f'z{g[0]} n{g[1]} {lab}'] = dict(mse_raw=s['mse_raw'], mse_scaled=s['mse_scaled'], mse_unclamped=s['mse_unclamped'], pds=s['pds_raw'],
                                               median_S=float(np.median(S.sum(1))), k_q_half0=kq[0].round(3).tolist())
        np.save(F / f'm3_sparse_z{g[0]}_n{g[1]}_h1.npy', sparse.astype(np.float32))
    if f1cal is not None:
        P = np.maximum(h.ctrl[None] + f1cal[idx75], 0); s = h.score(P)
        m3['F1 calibrated alone'] = dict(mse_raw=s['mse_raw'], mse_scaled=s['mse_scaled'], mse_unclamped=s['mse_unclamped'], pds=s['pds_raw'])
    res['M3'] = m3; print(pd.DataFrame(m3).T.round(4).to_string())
    # ---------------- U-arm
    dep = np.load(V / 'lib/pie_features/depmap_gene_effect.npy').astype(np.float32); cov = np.load(V / 'lib/pie_features/coverage_mask.npy')
    dcol = json.loads((V / 'lib/pie_features/sources.json').read_text())['columns'].index('depmap_gene_effect')
    frac = np.where(cov[:, dcol], (dep < -0.5).mean(1), np.nan); ess = frac >= 0.5
    def U_of(t):
        k = np.flatnonzero(t.w > 0); e = np.array([ess[gp[x]] if x in gp and np.isfinite(frac[gp[x]]) else False for x in t.tg[k]])
        ne = np.array([(x in gp) and np.isfinite(frac[gp[x]]) and not ess[gp[x]] for x in t.tg[k]])
        A = t.D[k].astype(np.float32); Vv_ = t.vD[k].astype(np.float32)
        u = np.nanmean(A[e], 0) - np.nanmean(A[ne], 0)
        var = np.nansum(Vv_[e], 0) / max(e.sum(), 1) ** 2 + np.nansum(Vv_[ne], 0) / max(ne.sum(), 1) ** 2
        return np.where(t.measured, u, np.nan), var, int(e.sum()), int(ne.sum())
    U = {}
    for nm in SRC + ['H1']:
        U[nm] = U_of(tabs[nm])
    names = SRC + ['H1']
    corr = pd.DataFrame([[float(pd.Series(U[a][0]).corr(pd.Series(U[b][0]))) for b in names] for a in names], index=names, columns=names)
    corr.to_csv(F / 'u1_corr.csv'); print(corr.round(3).to_string())
    res['U1'] = dict(n_essential={k: U[k][2] for k in names}, n_other={k: U[k][3] for k in names}, corr=corr.round(4).to_dict())
    num = np.zeros(len(genes)); den = np.zeros(len(genes))
    for nm in SRC:
        u, var, _, _ = U[nm]; okk = np.isfinite(u) & (var > 0); w = np.where(okk, 1 / np.where(okk, var, 1), 0); num += w * np.nan_to_num(u); den += w
    Ug = np.where(den > 0, num / np.maximum(den, 1e-30), 0)
    res['U1']['pooled_vs_H1'] = float(pd.Series(Ug).corr(pd.Series(U['H1'][0])))
    e = np.array([frac[gp[t]] if t in gp and np.isfinite(frac[gp[t]]) else np.nan for t in tg]); e = np.nan_to_num(e, nan=np.nanmean(e))
    ec = e - e.mean(); base = ec[:, None] * Ug[None, :]
    ks = {}
    for hh in (0, 1):
        m = (half == hh)[:, None] & fin
        ks[hh] = float((np.where(m, real * base, 0)).sum() / max((np.where(m, base ** 2, 0)).sum(), 1e-12))
    Kv = np.array([ks[1 - half[i]] for i in range(n)])
    uterm = Kv[:, None] * base; np.save(F / 'u3_term_h1.npy', uterm.astype(np.float32))
    u3 = {}
    for lab, A in (('U alone', uterm), ('F1 calibrated + U', None if f1cal is None else f1cal + uterm)):
        if A is None: continue
        P = np.maximum(h.ctrl[None] + A[idx75], 0); s = h.score(P)
        u3[lab] = dict(mse_raw=s['mse_raw'], mse_scaled=s['mse_scaled'], mse_unclamped=s['mse_unclamped'], pds=s['pds_raw'])
    res['U3'] = dict(k=ks, scores=u3, e_mean=float(e.mean()), essential_targets_H1=int((e >= 0.5).sum()))
    np.savez(F / 'u_signature.npz', Ug=Ug, frac=frac)
    print(json.dumps(res['U3'], indent=1), 'pooled U vs H1 U r =', res['U1']['pooled_vs_H1'])
    (F / 'mstudy.json').write_text(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()

"""C (calibration, scorer space) + F1 arrays.

  cd ~/vcc && PYTHONPATH=scripts/v5:analysis/model7 ~/venv/bin/python analysis/final/calib.py

Truth: H1's scorer-space table (Model 7 Step 0: lib/vcc_2026_5/m7/H1.npz, 300 knockdowns, summed raw counts of QC'd
cells). Model 2 v4's lfc prediction is put into scorer space on H1's control share: S(P0 2^lfc) - S(P0).
H1's knockdowns are split in half (seed 0); every fit is made on one half and applied to the other, then swapped.
C1 sign confidence: logistic regression of [sign(pred) == sign(real)] on sign agreement, number of sources, log1p re_z
  and the control-expression group (7 one-hot), rows = (knockdown, gene) with a nonzero prediction, gene measured in
  H1, own gene out; pred' = pred x (2p - 1).
C2 magnitude: per group q (control deciles 1-4, 5, ..., 10) the no-intercept slope k_q of real Delta on pred'
  (and, for the k-only arm, on pred).
Scoring: the 75-target reference at pseudobulk level (reference control + prediction). Full-H1 fits saved for A-C.
F1 arrays (log2 factors on the generator control share, via m7_score.to_log2_factor) + the called-gene masks
(top 2,000 by re_z of the unshrunk prediction among genes > 5 CPM in the destination's controls, own gene out).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from m7_score import to_log2_factor
from pb import H1, row_r
from v5lib import write_dir

V = Path.home() / 'vcc'
F = V / 'analysis/final'
TS = 5e4


def groups(S0, measured):
    q = np.zeros(len(S0), int); ok = measured & np.isfinite(S0)
    rk = pd.Series(S0[ok]).rank(pct=True, method='first').to_numpy(); q[ok] = np.maximum(np.minimum((rk * 10).astype(int), 9) - 3, 0)
    return q


def to_scorer(lfc, P0):
    P0 = np.nan_to_num(P0)[None]
    return (np.log1p(TS * P0 * np.exp2(lfc)) - np.log1p(TS * P0)).astype(np.float32)


def features(sa, ns, z, q):
    oh = np.eye(7)[q]
    return np.column_stack([sa, ns, np.log1p(z), oh])


def fit_sign(rows_X, y, seed=0):
    rng = np.random.default_rng(seed); idx = rng.choice(len(y), min(len(y), 400000), replace=False)
    return LogisticRegression(max_iter=500, C=1.0).fit(rows_X[idx], y[idx])


def main():
    res = {}
    z = np.load(F / 'm2v4_h1.npz', allow_pickle=False); tgt = z['targets'].astype(str)
    t7 = np.load(V / 'lib/vcc_2026_5/m7/H1.npz', allow_pickle=False); pos7 = {t: i for i, t in enumerate(t7['tg'].astype(str))}
    rows = np.array([pos7.get(t, -1) for t in tgt]); assert (rows >= 0).all()
    real = t7['D'][rows].astype(np.float32); measured = t7['measured']; P0 = t7['P0']; S0 = t7['S0']
    genes = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy(); gp = {g: i for i, g in enumerate(genes)}
    q = groups(S0, measured)
    pred = to_scorer(z['re'], P0); sa, ns, zz = z['sign_agree'], z['n_src'], z['re_z']
    n = len(tgt); half = np.zeros(n, int); half[np.random.default_rng(0).permutation(n)[: n // 2]] = 1
    mask = (pred != 0) & np.isfinite(real) & measured[None]
    for i, t in enumerate(tgt):
        if t in gp: mask[i, gp[t]] = False
    pprob = np.full(pred.shape, 0.5, np.float32); k_cal = {}; k_only = {}
    curve = []
    for h in (0, 1):
        tr = (half == h)[:, None] & mask; te = (half != h)
        Xtr = features(sa[tr], ns[tr], zz[tr], np.broadcast_to(q, pred.shape)[tr]); ytr = (np.sign(pred[tr]) == np.sign(real[tr])).astype(int)
        lr = fit_sign(Xtr, ytr)
        for i in np.flatnonzero(te):
            Xi = features(sa[i], ns[i], zz[i], q); pprob[i] = np.where(pred[i] != 0, lr.predict_proba(Xi)[:, 1], 0.5)
        sh = pred * (2 * pprob - 1)
        for nm, src, store in (('cal', sh, k_cal), ('konly', pred, k_only)):
            ks = np.zeros(7)
            for qq in range(7):
                m = tr & (q == qq)[None]
                ks[qq] = float((real[m] * src[m]).sum() / max((src[m] ** 2).sum(), 1e-12))
            store[h] = ks
        tem = te[:, None] & mask
        bins = np.linspace(0.5, 1.0, 11)
        pv = pprob[tem]; ov = (np.sign(pred[tem]) == np.sign(real[tem]))
        for lo, hi in zip(bins[:-1], bins[1:]):
            s = (pv >= lo) & (pv < hi)
            if s.sum(): curve.append(dict(fit_half=h, p_lo=lo, p_hi=hi, n=int(s.sum()), mean_p=float(pv[s].mean()), observed=float(ov[s].mean())))
        res[f'logistic_half{h}'] = dict(coef=lr.coef_.ravel().round(4).tolist(), intercept=float(lr.intercept_[0]),
                                        features=['sign_agree', 'n_src', 'log1p_re_z'] + [f'q{i}' for i in range(7)])
    # held-out k application: knockdown in half h' uses the fit from the other half
    kq_cal = np.stack([k_cal[1 - half[i]] for i in range(n)]); kq_only = np.stack([k_only[1 - half[i]] for i in range(n)])
    arms = {'pred': pred, "pred' = pred (2p-1)": pred * (2 * pprob - 1),
            "k_q pred' (calibrated)": kq_cal[np.arange(n)[:, None], q[None, :]] * pred * (2 * pprob - 1),
            'k_q pred (k-only)': kq_only[np.arange(n)[:, None], q[None, :]] * pred}
    res['k_q'] = dict(cal={h: v.round(4).tolist() for h, v in k_cal.items()}, konly={h: v.round(4).tolist() for h, v in k_only.items()})
    res['calibration_curve'] = curve
    # score on the 75-target reference
    h = H1(); idx = np.array([list(tgt).index(t) for t in h.targets])
    d_real = h.R - h.R.mean(0); sc = {}
    for nm, A in arms.items():
        P = np.maximum(h.ctrl[None] + A[idx], 0); s = h.score(P)
        sc[nm] = {k: s[k] for k in ('mse_raw', 'mse_scaled', 'mse_unclamped', 'pds_raw', 'pds_scaled')}
        sc[nm]['median_r_demeaned'] = float(np.nanmedian(row_r(P - P.mean(0), d_real, h.pan)))
    res['scores_pseudobulk'] = sc
    print(pd.DataFrame(curve).round(3).to_string(index=False)); print(json.dumps(res['k_q'], indent=1)); print(pd.DataFrame(sc).T.round(4).to_string())
    # full-H1 fits for A-C
    Xall = features(sa[mask], ns[mask], zz[mask], np.broadcast_to(q, pred.shape)[mask]); yall = (np.sign(pred[mask]) == np.sign(real[mask])).astype(int)
    lr_full = fit_sign(Xall, yall); p_full = np.full(pred.shape, 0.5, np.float32)
    for i in range(n): p_full[i] = np.where(pred[i] != 0, lr_full.predict_proba(features(sa[i], ns[i], zz[i], q))[:, 1], 0.5)
    shf = pred * (2 * p_full - 1)
    kfc = np.array([float((real[mask & (q == qq)[None]] * shf[mask & (q == qq)[None]]).sum() / max((shf[mask & (q == qq)[None]] ** 2).sum(), 1e-12)) for qq in range(7)])
    kfo = np.array([float((real[mask & (q == qq)[None]] * pred[mask & (q == qq)[None]]).sum() / max((pred[mask & (q == qq)[None]] ** 2).sum(), 1e-12)) for qq in range(7)])
    res['full_fit'] = dict(coef=lr_full.coef_.ravel().tolist(), intercept=float(lr_full.intercept_[0]), k_cal=kfc.tolist(), k_only=kfo.tolist())
    # ---- F1 arrays: H1 (held-out fits) and A-C (full fits)
    called = np.zeros(pred.shape, bool); gate = h.line.gate
    zr = np.where(gate[None], zz, -1)
    for i, t in enumerate(tgt):
        if t in gp: zr[i, gp[t]] = -1
    top = np.argpartition(-zr, 1999, axis=1)[:, :2000]; called[np.arange(n)[:, None], top] = True; called &= zr > 0
    share = h.line.share
    f1 = {'f1_cal': to_log2_factor(arms["k_q pred' (calibrated)"], share), 'f1_konly': to_log2_factor(arms['k_q pred (k-only)'], share),
          'called_top2000': called}
    write_dir(F / 'h1_pred', V / 'lib/vcc_2026_4/runs/dev_nopie/h1_pred_g3', f1, 'H1',
              ['K562', 'RPE1', 'HepG2', 'Jurkat', 'HCT116', 'HEK293T', 'CD4_Stim48hr'])
    base = dict(call='z', z=0.0, min_effect=0.0, alpha=1.0, detect_floor=False, called_from='called_top2000')
    v = [dict(base, method=a, narrow_m=m, label=f'{a}_m{m:g}') for a in ('f1_cal', 'f1_konly') for m in (0.5, 0.25)]
    (F / 'variants_f1.json').write_text(json.dumps(v, indent=1))
    np.savez(F / 'calib_full.npz', coef=lr_full.coef_.ravel(), intercept=lr_full.intercept_, k_cal=kfc, k_only=kfo)
    (F / 'calib.json').write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res['full_fit'], indent=1))


if __name__ == '__main__':
    main()

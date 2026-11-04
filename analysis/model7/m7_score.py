"""Model 7, Steps 4-5 (H1): convert scorer-space predictions to generator factors, check the round trip, score at
pseudobulk level, and write the prediction folder for the cell-based H1 benchmark.

  cd ~/vcc && PYTHONPATH=scripts/v5 ~/venv/bin/python analysis/model7/m7_score.py

Conversion (Step 4): P_hat = (exp(log1p(5e4 C) + Delta_hat) - 1) / 5e4, factor = P_hat / C on the generator's control
pool share C; a gene with C = 0 and P_hat > 0 takes the generator's zero-baseline induction (share 0.5e-6 (f - 1)).
The generator re-projects the composition, so the realised pseudobulk can differ from S0 + Delta_hat; the check
generates 5 knockdowns (floor off, 400 cells) and reports the error.
Pseudobulk scoring: (a) direct: the reference's control pseudobulk + Delta_hat; (b) the generator's noise-free path
(calibrate: gate, own gene, +/-8 clip, optional floor; then composition projection), floor off and on.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from pb import H1, row_r
from v5lib import arm

V = Path.home() / 'vcc'
OUT = V / 'analysis/model7'
sys.path.insert(0, str(V / 'code/vcc_2026_3c'))
from vcc3_generate import generate_one, pick_cells  # noqa: E402
from vcc2_common import seed_for  # noqa: E402

TS = 5e4
QLAB = ['dec1-4', 'dec5', 'dec6', 'dec7', 'dec8', 'dec9', 'dec10']


def to_log2_factor(Dhat, share):
    """Scorer-space Delta_hat (rows x G) -> log2 factors on the control share."""
    C = share[None, :].astype(np.float64)
    Ph = (np.exp(np.log1p(TS * C) + Dhat.astype(np.float64)) - 1) / TS
    Ph = np.maximum(Ph, 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        f = np.where(C > 0, Ph / np.where(C > 0, C, 1), 1 + Ph / 0.5e-6)
    return np.clip(np.log2(np.maximum(f, 2 ** -8)), -8, 8).astype(np.float32)


def groups_of(S0):
    rk = pd.Series(S0).rank(pct=True, method='first').to_numpy(); dec = np.minimum((rk * 10).astype(int), 9)
    return np.maximum(dec - 3, 0)


def main():
    h = H1(); res = {}
    z = np.load(OUT / 'h1_predictions_scorer_space.npz', allow_pickle=True); tg = z['tg'].astype(str)
    pos = {t: i for i, t in enumerate(tg)}; idx300 = np.array([pos.get(t, -1) for t in h.t300])
    missing = int((idx300 < 0).sum())
    arms_scorer = {k: z[k] for k in ('model7', 'model7_single_beta', 'model7_denoised')}
    def to300(A):
        out = np.zeros((len(h.t300), A.shape[1]), np.float32); ok = idx300 >= 0; out[ok] = A[idx300[ok]]; return out
    arms_scorer = {k: to300(v) for k, v in arms_scorer.items()}
    share = h.line.share
    arms_log2 = {k: to_log2_factor(v, share) for k, v in arms_scorer.items()}
    m2 = np.load(OUT / 'h1_model2_effects_v4.npz', allow_pickle=True); p2 = {t: i for i, t in enumerate(m2['tg'].astype(str))}
    arms_log2['model2_effects_v4'] = np.stack([m2['re'][p2[t]] if t in p2 else np.zeros(h.G, np.float32) for t in h.t300])
    arms_log2['file1_transformer_s2'] = h.full(arm(h.delta, s=2)).astype(np.float32)
    # ---- Step 4 round trip on 5 knockdowns (floor off)
    rng = np.random.default_rng(0); pick5 = rng.choice(h.T, 5, replace=False)
    common = pick_cells(h.line, np.zeros(h.G), h.cfg, np.random.default_rng(seed_for(2026, 'H1', 'base_controls')))
    from vcc3_generate import calibrate
    Dc = calibrate(arms_log2['model7'][h.idx], h.line, h.targets, h.cfg)
    errs = []
    for i in pick5:
        X = generate_one(h.line, Dc[i], h.own[i] if h.own[i] >= 0 else None, h.kd[i], h.cfg, pick=common[:400], line_key='H1', target_key=h.targets[i])[0]
        x = np.asarray(X.sum(0)).ravel().astype(np.float64); S_gen = np.log1p(TS * x / x.sum())
        S_int = np.log1p(TS * share) + arms_scorer['model7'][h.idx[i]]
        keep = np.ones(h.G, bool)
        if h.own[i] >= 0: keep[h.own[i]] = False
        expr = keep & h.line.gate
        S_pool = np.log1p(TS * np.asarray(h.line.X[common[:400]].sum(0)).ravel() / h.line.X[common[:400]].sum())
        errs.append(dict(target=h.targets[i], max_abs_err_all=float(np.abs(S_gen - S_int)[keep].max()),
                         max_abs_err_expressed=float(np.abs(S_gen - S_int)[expr].max()), median_abs_err_expressed=float(np.median(np.abs(S_gen - S_int)[expr])),
                         pool_sampling_median_abs_expressed=float(np.median(np.abs(S_pool - np.log1p(TS * share))[expr])),
                         intended_shift_median_abs_expressed=float(np.median(np.abs(arms_scorer['model7'][h.idx[i]])[expr]))))
    res['step4_roundtrip'] = errs; print(pd.DataFrame(errs).round(4).to_string(index=False))
    # ---- Step 5 pseudobulk scoring
    d_real = h.R - h.R.mean(0); rows = []
    qH = groups_of(h.ctrl)
    def record(name, path, P):
        r = dict(arm=name, path=path); s = h.score(P)
        r.update({k: s[k] for k in ('mse_raw', 'mse_scaled', 'mse_unclamped', 'pds_raw', 'pds_scaled')})
        r['median_r_demeaned'] = float(np.nanmedian(row_r(P - P.mean(0), d_real, h.pan)))
        for q, lab in [(0, 'weak'), (1, 'mid'), (2, 'strong')]:
            r[f'mse_raw_{lab}'] = h.score(P, np.flatnonzero(h.terc == q))['mse_raw']
        E = (P - h.R) ** 2
        for i in h.rows:
            if h.own[i] >= 0: E[i, h.own[i]] = 0
        sden = h.den.sum()
        for qq, lab in enumerate(QLAB): r[f'plugin_{lab}'] = float((E[:, qH == qq].sum(1) / h.nrow).sum() / sden)
        rows.append(r); print(f"{name:28s} {path:22s} mse {r['mse_raw']:.4f} ({r['mse_unclamped']:+.4f}) pds {r['pds_raw']:.4f} r {r['median_r_demeaned']:.3f}", flush=True)
    for k, v in arms_scorer.items(): record(k, 'direct (ctrl + Delta)', h.ctrl[None] + v[h.idx])
    for k, v in arms_log2.items():
        record(k, 'generator path, no floor', h.bulk(v)); record(k, 'generator path, floor', h.bulk(v, floor=True))
    df = pd.DataFrame(rows); df.to_csv(OUT / 'step5_pseudobulk.csv', index=False); res['step5_pseudobulk'] = df.to_dict('records')
    res['h1_targets_missing_from_model7'] = missing
    # ---- prediction folder for the cell benchmark
    pdir = OUT / 'h1_pred'; pdir.mkdir(exist_ok=True)
    src = V / 'lib/vcc_2026_4/runs/dev_nopie/h1_pred_g3'
    for f in ('genes.npy', 'targets.npy', 'ontarget_lfc.npy'): (pdir / f).write_bytes((src / f).read_bytes())
    models = []
    for k in ('model7', 'model7_single_beta', 'model7_denoised', 'model2_effects_v4'):
        np.save(pdir / f'{k}.npy', arms_log2[k]); models.append(k)
    (pdir / 'prediction_report.json').write_text(json.dumps(dict(schema='vcc_2026_3_ridge_v1', heldout_background='H1',
        source_backgrounds=['K562', 'RPE1', 'HepG2', 'Jurkat', 'HCT116', 'HEK293T', 'CD4_Rest', 'CD4_Stim8hr', 'CD4_Stim48hr'],
        models=models, note='Model 7 (analysis/model7), log2 factors on the generator control share'), indent=1))
    base = dict(call='z', z=0.0, min_effect=0.0, alpha=1.0)
    v1 = [dict(base, method=m, detect_floor=f, label=f'{m}_{"floor" if f else "nofloor"}') for m in ('model7', 'model7_single_beta') for f in (True, False)]
    v2 = [dict(base, method=m, detect_floor=f, label=f'{m}_{"floor" if f else "nofloor"}') for m in ('model7_denoised', 'model2_effects_v4') for f in (True, False)]
    (OUT / 'variants_p1.json').write_text(json.dumps(v1, indent=1)); (OUT / 'variants_p2.json').write_text(json.dumps(v2, indent=1))
    (OUT / 'step45.json').write_text(json.dumps(res, indent=1, default=float))
    print(df[['arm', 'path', 'mse_raw', 'mse_unclamped', 'pds_raw', 'median_r_demeaned', 'mse_raw_weak', 'mse_raw_mid', 'mse_raw_strong']].round(4).to_string(index=False))


if __name__ == '__main__':
    main()

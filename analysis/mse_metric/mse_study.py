"""vcc2026 mse metric study on H1's 75-target harness reference (lib/vcc_2026_3/val_final). Read-only; no upload.

  cd ~/vcc && PYTHONPATH=scripts/v5 ~/venv/bin/python analysis/mse_metric/mse_study.py [s1 s2 s3 s4]

Uses the v5 pseudobulk harness (scripts/v5/pb.py: the scorer's own per-knockdown formula on cached real moments, the
generator-half H1 controls with guide stratification, dev_nopie's zero-shot H1 predictions).
Pseudobulk space = cell-eval2 bulk_lognorm: log1p(5e4 x summed counts / total). "True" quantities in S2 are the
reference's own pseudobulks R (so they share the reference's sampling noise; the unbiased correction then makes a
perfect copy read slightly below 0).
"""
import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

from pb import H1, V, row_r

sys.path.insert(0, str(V / 'code/vcc_2026_3c'))
from vcc3_generate import Cfg, calibrate, generate_one, kd_residual, pick_cells, project_factors  # noqa: E402
from vcc2_common import seed_for  # noqa: E402

OUT = V / 'analysis/mse_metric'
BTS = 5e4
REF = V / 'lib/vcc_2026_3/val_final'
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200); pd.set_option('display.max_columns', 30)


def per_kd(h, P):
    """Per-knockdown numerator (noise-free prediction: C_pred = 0) and the reference's denominator."""
    d = P - h.R
    num = (np.einsum('ij,ij->i', d, d) - h.Creal - h.own_sq(P, h.R)) / h.nrow
    return num, h.den


def bulk_with(h, D75, kd, cfg):
    Dc = calibrate(np.asarray(D75, np.float32), h.line, h.targets, cfg)
    P = np.empty((h.T, h.G))
    for i in h.rows:
        _, proj, _ = project_factors(h.line.share, Dc[i], h.own[i] if h.own[i] >= 0 else None, kd[i])
        P[i] = np.log1p(BTS * proj)
    return P, Dc


# ------------------------------------------------------------------ S1
def s1(h, res):
    P_dev = h.bulk(h.full(h.delta)); num, den = per_kd(h, P_dev)
    bl = pd.read_csv(REF / 'baseline/baseline_results.csv'); bl.columns = ['perturbation', 'metric', 'value']
    bl = bl.pivot_table(index='perturbation', columns='metric', values='value').loc[h.targets]
    bnum, bden = bl['expr_mse_unbiased_capped'].to_numpy(), bl['expr_distance_unbiased'].to_numpy()
    out = dict(F_ratio_of_sums=float(num.sum() / den.sum()), F_harness=h.score(P_dev)['mse_raw'],
               baseline_ratio_of_sums=float(bnum.sum() / bden.sum()),
               den_matches_baseline_den=float(np.abs(bden - den).max()),
               F_mean_of_per_kd_ratios=float(np.mean(num / den)), baseline_mean_of_per_kd_ratios=float(np.mean(bnum / bden)),
               F_median_per_kd_ratio=float(np.median(num / den)), baseline_median_per_kd_ratio=float(np.median(bnum / bden)))
    rows = []
    for q, name in enumerate(['weak', 'mid', 'strong']):
        s = h.terc == q
        rows.append(dict(tercile=name, n=int(s.sum()), den_share=den[s].sum() / den.sum(),
                         ours_ratio_of_sums=num[s].sum() / den[s].sum(), base_ratio_of_sums=bnum[s].sum() / den[s].sum(),
                         ours_mean_ratio=np.mean(num[s] / den[s]), base_mean_ratio=np.mean(bnum[s] / den[s]),
                         ours_num_share=num[s].sum() / num.sum(), base_num_share=bnum[s].sum() / bnum.sum()))
    t = pd.DataFrame(rows)
    out['mean_of_tercile_ratios'] = dict(ours=float(t.ours_ratio_of_sums.mean()), base=float(t.base_ratio_of_sums.mean()))
    t.to_csv(OUT / 's1_terciles.csv', index=False)
    pd.DataFrame(dict(target=h.targets, num_dev=num, num_baseline=bnum, den=den, tercile=h.terc)).to_csv(OUT / 's1_per_knockdown.csv', index=False)
    res['S1'] = out; print(json.dumps(out, indent=1)); print(t.round(4).to_string(index=False))


# ------------------------------------------------------------------ S2
def s2(h, res):
    c = h.ctrl[None, :]; T = h.R - c; m = T.mean(0); rng = np.random.default_rng(0)
    P_dev = h.bulk(h.full(h.delta)); Ddev = P_dev - c
    d_real = h.R - h.R.mean(0)
    arms = {}
    def add(name, P, pearson=False):
        P = np.maximum(np.broadcast_to(P, h.R.shape), 0.0)                    # log1p pseudobulks are >= 0
        r = h.score(P); r = {k: r[k] for k in ('mse_raw', 'mse_scaled', 'mse_unclamped', 'pds_raw', 'pds_scaled')}
        if pearson:
            r['median_pearson_delta'] = float(np.median(row_r(P - c, T, h.pan)))
            r['median_pearson_demeaned'] = float(np.median(row_r(P - P.mean(0), d_real, h.pan)))
        arms[name] = r
    add('true R_k (reference itself)', h.R)
    add('dev_nopie (s=1, no floor)', P_dev, True)
    add('zero effect (control)', np.broadcast_to(c, h.R.shape).copy())
    add('true shared response m for every knockdown', c + m[None])
    for target in (0.2, 0.4, 0.6, 0.8):
        sig = np.sqrt(1 / target ** 2 - 1)
        eps = rng.standard_normal(T.shape) * (T.std(1, keepdims=True) * sig)
        add(f'true Delta + noise (target r {target})', c + T + eps, True)
        # the calibrated version: one global regression shrink a = r^2 (the mse-optimal scalar for isotropic noise)
        add(f'true Delta + noise (target r {target}), shrunk x r^2', c + target ** 2 * (T + eps), True)
    for target in (0.9, 0.95):
        sig = np.sqrt(1 / target ** 2 - 1)
        eps = rng.standard_normal(T.shape) * (T.std(1, keepdims=True) * sig)
        add(f'true Delta + noise (target r {target}), shrunk x r^2', c + target ** 2 * (T + eps), True)
    for a in (0.5, 0.8, 1.2, 1.5): add(f'true Delta x {a}', c + a * T)
    ck = (T @ m) / (m @ m); add('amplitude only: c_k m, c_k = Delta_k.m/|m|^2', c + ck[:, None] * m[None])
    A = np.abs(T).copy()
    for i in h.rows:
        if h.own[i] >= 0: A[i, h.own[i]] = -1
    for N in (100, 1000, 5000):
        thr = -np.sort(-A, 1)[:, N - 1:N]; add(f'true Delta on top {N} |Delta_k| genes', c + np.where(A >= thr, T, 0))
    add('true Delta on the 5,000-gene panel only', c + T * h.run_set[None])
    add('true Delta off the panel only', c + T * ~h.run_set[None])
    add('true sign x dev_nopie magnitude', c + np.sign(T) * np.abs(Ddev), True)
    add('dev_nopie sign x true magnitude', c + np.sign(Ddev) * np.abs(T), True)
    df = pd.DataFrame(arms).T; df.index.name = 'arm'; df.to_csv(OUT / 's2_arms.csv')
    b, rp = h.bm, h.rm
    lb = {f'rank {k}': dict(scaled=s, raw=rp + (1 - s) * (b - rp)) for k, s in (('1', 0.786), ('10', 0.481))}
    res['S2'] = dict(arms=df.to_dict('index'), anchors=dict(baseline_raw=b, replicate_raw=rp), leaderboard_to_raw=lb,
                     c_k=dict(mean=float(ck.mean()), sd=float(ck.std()), min=float(ck.min()), max=float(ck.max())))
    print(df.round(4).to_string()); print(json.dumps(res['S2']['leaderboard_to_raw'], indent=1), res['S2']['c_k'])


# ------------------------------------------------------------------ S3 (generator)
_G = None


def _gen(args):
    """One knockdown through the submission generator -> (pseudobulk, jackknife C_pred)."""
    from cell_eval2.moments import jackknife_correction
    arm, i = args
    h, Dc, kd, common, cfg, libsizes = _G[arm]
    local = Cfg(argparse.Namespace(**vars(cfg))); local.n_cells = 400
    local.library_size = libsizes[i] if libsizes is not None else None
    t = h.targets[i]
    X = generate_one(h.line, Dc[i], h.own[i] if h.own[i] >= 0 else None, kd[i], local, pick=common[:400], line_key='H1', target_key=t)[0]
    x = np.asarray(X.sum(0)).ravel().astype(np.float64)
    jk = float(np.asarray(jackknife_correction(X.astype(np.float64), np.zeros(X.shape[0], np.int64), 1, BTS)).ravel()[0])
    return arm, i, np.log1p(BTS * x / x.sum()), jk


def score_generated(h, P, jk_pred):
    """The capped numerator with a predicted-side correction: #247 cap per row, #348 budget scaling (whole panel)."""
    from cell_eval2.metrics.delta import _across_pert_budget, _row_weights
    excl = {i: int(h.own[i]) for i in h.rows if h.own[i] >= 0}
    term = np.minimum(jk_pred, h.Creal)
    w = _row_weights(h.T, h.G, excl); claim = float((term * w).sum())
    budget = _across_pert_budget(P, excl, w); ratio = min(1.0, budget / claim) if claim > 0 else 1.0
    d = P - h.R
    num = (np.einsum('ij,ij->i', d, d) - term * ratio - h.Creal - h.own_sq(P, h.R)) / h.nrow
    u = num.sum() / h.den.sum()
    pv = float(h.pds_rows(P).mean())
    return dict(mse_raw=float(u), mse_scaled=h.scale(h.MSE, u), mse_unclamped=(h.bm - u) / (h.bm - h.rm), pds_raw=pv,
                pds_scaled=h.scale(h.PDS, pv), pred_correction_capped_rows=int((jk_pred > h.Creal).sum()),
                budget_over_claim=float(budget / claim) if claim > 0 else np.nan,
                pred_correction_share_of_den=float((term * ratio / h.nrow).sum() / h.den.sum())), num


def s3(h, res, workers=4):
    global _G
    genes = list(h.genes)
    with h5py.File(V / 'lib/effects_v4/arc_h1.h5ad', 'r') as f:
        obs = read_elem(f['obs']); L = f['layers/lfc'][:]
    pos = {t: i for i, t in enumerate(obs['target_gene'].astype(str))}
    Dtrue = np.nan_to_num(L[[pos[t] for t in h.targets]]).astype(np.float32)
    own_lfc = np.array([Dtrue[i, h.own[i]] if h.own[i] >= 0 else np.nan for i in h.rows])
    kd_true = kd_residual(own_lfc, h.cfg)
    Ddev = h.full(h.delta)[h.idx].astype(np.float32)
    with h5py.File(V / 'lib/qc_cells/arc_h1_qc.h5ad', 'r') as f:
        qobs = read_elem(f['obs'])
    med_umi = qobs.groupby(qobs['target_gene'].astype(str))['n_umi'].median()
    libsizes = np.array([float(med_umi.get(t, np.nan)) for t in h.targets])
    mz = np.load(next((REF / 'real_cache').glob('*.moments.npz'))); cnt = dict(zip(mz['perts'].astype(str), mz['counts']))
    pool_cfg = Cfg(argparse.Namespace(**vars(h.cfg))); pool_cfg.n_cells = int(max([400] + [cnt[t] for t in h.targets]))
    common = pick_cells(h.line, np.zeros(h.G), pool_cfg, np.random.default_rng(seed_for(2026, 'H1', 'base_controls')))
    arms, rows, gen_nums = {}, [], {}
    for name, D, kd in (('true', Dtrue, kd_true), ('dev_nopie', Ddev, h.kd)):
        for floor in (False, True):
            cfg = h.cfg_floor if floor else h.cfg
            P, Dc = bulk_with(h, D, kd, cfg); r = h.score(P)
            rows.append(dict(delta=name, floor=floor, path='direct pseudobulk', libsize='-', **{k: r[k] for k in ('mse_raw', 'mse_scaled', 'mse_unclamped', 'pds_raw', 'pds_scaled')}))
            for lib in (False, True):
                arms[(name, floor, lib)] = (h, Dc, kd, common, cfg, libsizes if lib else None)
    _G = arms; t0 = time.time()
    jobs = [(a, i) for a in arms for i in h.rows]
    with mp.get_context('fork').Pool(workers) as pool:
        got = pool.map(_gen, jobs, chunksize=4)
    print(f'  generated {len(jobs)} knockdowns in {time.time() - t0:.0f}s', flush=True)
    Pg = {a: np.empty((h.T, h.G)) for a in arms}; jk = {a: np.empty(h.T) for a in arms}
    for a, i, p, j in got: Pg[a][i] = p; jk[a][i] = j
    for (name, floor, lib) in arms:
        r, num = score_generated(h, Pg[(name, floor, lib)], jk[(name, floor, lib)]); gen_nums[(name, floor, lib)] = num
        rows.append(dict(delta=name, floor=floor, path='generator (400 cells)', libsize='on' if lib else 'off', **r))
    df = pd.DataFrame(rows); df.to_csv(OUT / 's3_generator.csv', index=False)
    sub = json.loads((V / 'results/vcc_2026_5_v5_file1_s2_floor/generator_report.json').read_text())['settings']
    old = json.loads((V / 'results/vcc_2026_4_final_nopie_top5000/generator_report.json').read_text())['settings']
    res['S3'] = dict(table=df.to_dict('records'), submission_settings=dict(file1_detect_floor=sub['detect_floor'], file1_library_size=sub['library_size'],
                     v5_0638_detect_floor=old['detect_floor'], v5_0638_library_size=old['library_size'], z_req=sub['z_req']),
                     libsize_median_umi=dict(median=float(np.nanmedian(libsizes)), min=float(np.nanmin(libsizes)), max=float(np.nanmax(libsizes))),
                     generator_control_median_umi=float(np.median(h.line.lib[common[:400]])))
    np.savez(OUT / 's3_generated_dev.npz', P=Pg[('dev_nopie', True, False)], jk=jk[('dev_nopie', True, False)],
             P_nofloor=Pg[('dev_nopie', False, False)], jk_nofloor=jk[('dev_nopie', False, False)])
    print(df.round(4).to_string(index=False)); print(json.dumps({k: v for k, v in res['S3'].items() if k != 'table'}, indent=1))


# ------------------------------------------------------------------ S4
def decompose(h, P, label):
    c = h.ctrl; E = P - h.R
    mask = np.ones_like(E, bool)
    for i in h.rows:
        if h.own[i] >= 0: mask[i, h.own[i]] = False
    E0 = np.where(mask, E, 0.0); Ek = (E0 ** 2).sum(1)
    sden = h.den.sum(); unit = lambda x: float(np.sum(x / h.nrow) / sden)      # squared-error sums -> raw mse units
    num = (Ek - h.Creal) / h.nrow; total = float(num.sum() / sden)
    out = dict(model=label, raw_mse=total, plugin_error=unit(Ek), real_noise_correction=-unit(h.Creal))
    for q, name in enumerate(['weak', 'mid', 'strong']):
        s = h.terc == q; out[f'tercile_{name}'] = float(num[s].sum() / sden)
    worst = np.argsort(-num)[:max(1, int(round(0.1 * h.T)))]
    out['worst_10pct_share_of_numerator'] = float(num[worst].sum() / num.sum())
    out['worst_10pct_targets'] = [h.targets[i] for i in worst]
    G2 = (E0 ** 2) / h.nrow[:, None]
    out['on_panel_5000'] = float(G2[:, h.run_set].sum() / sden); out['off_panel'] = float(G2[:, ~h.run_set].sum() / sden)
    rk = pd.Series(c).rank(method='first').to_numpy(); dec = np.minimum((rk - 1) * 10 // len(c), 9).astype(int)
    for k in range(10): out[f'ctrl_decile_{k + 1}'] = float(G2[:, dec == k].sum() / sden)
    # components: shared (m_hat - m), amplitude along the true shared response, specific (rest); own gene left out
    mh, m = np.where(mask, P, 0).sum(0) / mask.sum(0), np.where(mask, h.R, 0).sum(0) / mask.sum(0)
    shared = np.where(mask, (mh - m)[None, :], 0.0)
    dd = E0 - shared
    u = (m - c) / np.linalg.norm(m - c)
    amp = (dd @ u)[:, None] * u[None, :] * mask
    spec = dd - amp
    for nm, X in (('shared', shared), ('amplitude', amp), ('specific', spec)):
        out[f'component_{nm}'] = unit((X ** 2).sum(1))
    out['component_cross_terms'] = out['plugin_error'] - out['component_shared'] - out['component_amplitude'] - out['component_specific']
    return out


def s4(h, res):
    P_dev = h.bulk(h.full(h.delta))
    from vcc3_generate import load_deltas
    m2 = load_deltas(V / 'lib/vcc_2026_3/validation/pred/H1', list(h.t300), list(h.genes), 're')
    z = np.load(V / 'lib/vcc_2026_4/runs/dev_nopie/h1_pred_g3/raw_predictions.npz')
    models = {'dev_nopie transformer': P_dev, 'Model 2 (vcc_2026_3 re, held-out H1)': h.bulk(m2),
              'effects_v4 evidence mean (pw_mean)': h.bulk(h.full(z['pw_mean'].astype(np.float64)))}
    rows = [decompose(h, P, k) for k, P in models.items()]
    g = OUT / 's3_generated_dev.npz'
    if g.exists():
        zz = np.load(g)
        for key, lab in (('P', 'dev_nopie through generator (submission: floor on, libsize off)'), ('P_nofloor', 'dev_nopie through generator (floor off)')):
            r = decompose(h, zz[key], lab)
            jk = zz['jk' if key == 'P' else 'jk_nofloor']
            r['raw_mse_scored_with_pred_correction'] = score_generated(h, zz[key], jk)[0]['mse_raw']
            rows.append(r)
    df = pd.DataFrame(rows).set_index('model'); df.to_csv(OUT / 's4_decomposition.csv')
    res['S4'] = df.drop(columns=['worst_10pct_targets']).to_dict('index')
    print(df.drop(columns=['worst_10pct_targets']).T.round(4).to_string())


def s4b(h, res):
    """H1 batch: each scored knockdown's batch (gem group) composition, how much of d_k the batch-mean control profiles
    explain, and how far apart batch-wise control means are relative to ||d_k||. Batch composition comes from the
    unfiltered H1 obs (the reference's own cells); batch control profiles from the QC'd H1 controls (lib/qc_cells, the
    same space: log1p(5e4 x summed counts / total)); noise of a batch mean from two random halves of its controls."""
    import scipy.sparse as sp
    with h5py.File(V / 'lib/h1/H1_full_2025.h5ad', 'r') as f:
        obs = read_elem(f['obs'])
    lab = obs['target_gene'].astype(str).to_numpy(); bat = obs['batch'].astype(str).to_numpy()
    batches = sorted(set(bat)); B = len(batches); bpos = {b: i for i, b in enumerate(batches)}
    ctrl_frac = pd.Series(bat[lab == 'non-targeting']).value_counts(normalize=True).reindex(batches).fillna(0).to_numpy()
    W = np.zeros((h.T, B))
    for i, t in enumerate(h.targets):
        W[i] = pd.Series(bat[lab == t]).value_counts(normalize=True).reindex(batches).fillna(0).to_numpy()
    tv = 0.5 * np.abs(W - ctrl_frac[None]).sum(1)
    comp = dict(n_batches=B, runs=sorted({b.rsplit('_', 1)[0] for b in batches}),
                cells_per_kd_median=float(np.median([np.sum(lab == t) for t in h.targets])),
                batches_per_kd_median=float(np.median((W > 0).sum(1))), max_batch_share_median=float(np.median(W.max(1))),
                max_batch_share_max=float(W.max()), tv_distance_to_control_composition=dict(median=float(np.median(tv)), max=float(tv.max())),
                uniform_share=1 / B)
    # batch-wise control pseudobulks (and two random halves for their sampling noise), QC'd controls
    path = V / 'lib/qc_cells/arc_h1_qc.h5ad'; rng = np.random.default_rng(0)
    with h5py.File(path, 'r') as f:
        qo = read_elem(f['obs']); X = f['X']; ip = X['indptr'][:]; Gs = int(X.attrs['shape'][1])
        qv = read_elem(f['var']); qgenes = qv.index.astype(str).to_numpy()
        cr = np.flatnonzero(qo['is_control'].to_numpy(bool)); half = rng.integers(0, 2, len(cr))
        code = np.full(len(qo), -1); code[cr] = np.array([bpos[b] for b in qo['batch'].astype(str).to_numpy()[cr]]) * 2 + half
        S = np.zeros((2 * B, Gs))
        for s0 in range(0, len(qo), 50000):
            e0 = min(s0 + 50000, len(qo)); r = np.flatnonzero(code[s0:e0] >= 0)
            if not len(r): continue
            a, b = int(ip[s0]), int(ip[e0])
            M = sp.csr_matrix((X['data'][a:b], X['indices'][a:b], ip[s0:e0 + 1] - a), shape=(e0 - s0, Gs))[r]
            Gm = sp.csr_matrix((np.ones(len(r)), (code[s0:e0][r], np.arange(len(r)))), shape=(2 * B, len(r)))
            S += (Gm @ M).toarray()
    gp = {g: i for i, g in enumerate(qgenes)}; cols = np.array([gp.get(g, -1) for g in h.genes]); ok = cols >= 0
    def to_panel(Y):
        out = np.zeros((Y.shape[0], h.G)); out[:, ok] = Y[:, cols[ok]]; return out
    S = to_panel(S); SB = S[0::2] + S[1::2]
    bulk = lambda Y: np.log1p(BTS * Y / Y.sum(1, keepdims=True))
    Cb, CA, CBh = bulk(SB), bulk(S[0::2]), bulk(S[1::2]); Cbar = bulk(SB.sum(0, keepdims=True))[0]
    pan = h.pan
    dist_raw = np.linalg.norm((Cb - Cbar)[:, pan], axis=1)
    noise_b = ((CA - CBh)[:, pan] ** 2).sum(1) / 4                                   # sampling variance of each batch mean
    dist_nc = np.sqrt(np.maximum(dist_raw ** 2 - noise_b, 0))
    dk_nc = h.size                                                                      # noise-corrected ||d_k|| of the reference
    runs = np.array([b.rsplit('_', 1)[0] for b in batches])
    run_prof = {r: bulk(SB[runs == r].sum(0, keepdims=True))[0] for r in sorted(set(runs))}
    run_dist = {r: float(np.linalg.norm((v - Cbar)[pan])) for r, v in run_prof.items()}
    ctrl = dict(batch_control_dist_to_pooled_raw=dict(median=float(np.median(dist_raw)), max=float(dist_raw.max())),
                batch_control_dist_noise_corrected=dict(median=float(np.median(dist_nc)), max=float(dist_nc.max())),
                batch_mean_sampling_noise_norm_median=float(np.median(np.sqrt(noise_b))),
                run_control_dist_to_pooled=run_dist,
                typical_dk_noise_corrected=dict(median=float(np.median(dk_nc)), q25=float(np.percentile(dk_nc, 25)), q75=float(np.percentile(dk_nc, 75))),
                ratio_median_batch_dist_nc_to_median_dk=float(np.median(dist_nc) / np.median(dk_nc)))
    # regress d_k on batch profiles
    d = (h.R - h.R.mean(0))[:, pan]
    E = (W @ (Cb - Cbar))[:, pan]; E = E - E.mean(0)                                   # composition-implied batch deviation
    dd = (d * d).sum(1); dnc = np.maximum(dd - h.Creal, 1e-12)
    expl_fixed = 1 - ((d - E) ** 2).sum(1) / dd                                         # beta = 1, no fitting
    beta = (d * E).sum(1) / np.maximum((E * E).sum(1), 1e-30); expl_beta = beta ** 2 * (E * E).sum(1) / dd
    Bmat = (Cb - Cbar)[:, pan].T                                                         # G x B, free fit (48 params / kd)
    Q, _ = np.linalg.qr(Bmat); proj = d @ Q; expl_free = (proj ** 2).sum(1) / dd
    noise_share = h.Creal / dd
    reg = dict(pooled_R2_composition_beta1=float(1 - ((d - E) ** 2).sum() / dd.sum()),
               median_R2_composition_beta1=float(np.median(expl_fixed)),
               median_R2_composition_fitted_beta=float(np.median(expl_beta)), median_beta=float(np.median(beta)),
               pooled_R2_free_48_batch_profiles=float((proj ** 2).sum() / dd.sum()), median_R2_free_48=float(np.median(expl_free)),
               chance_R2_free_48_random_directions=float(B / pan.sum()),
               median_real_noise_share_of_dk2=float(np.median(noise_share)),
               pooled_free_R2_of_noise_corrected_energy=float((proj ** 2).sum() / dnc.sum()))
    by_t = pd.DataFrame(dict(tercile=h.terc, R2_free=expl_free, R2_beta=expl_beta, tv=tv, maxshare=W.max(1))).groupby('tercile').median()
    res['S4b'] = dict(composition=comp, controls=ctrl, regression=reg, by_tercile=by_t.to_dict('index'),
                      lines_A_C_batch_labels='none: obs = target_gene, context, ntc_id (46 control guide ids); var and uns empty')
    pd.DataFrame(dict(target=h.targets, tv_to_control_comp=tv, max_batch_share=W.max(1), R2_composition=expl_fixed,
                      R2_fitted_beta=expl_beta, R2_free48=expl_free, real_noise_share=noise_share, dk_nc=dk_nc)).to_csv(OUT / 's4b_per_knockdown.csv', index=False)
    print(json.dumps(res['S4b'], indent=1, default=float))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('parts', nargs='*', default=['s1', 's2', 's3', 's4'])
    a = ap.parse_args(); OUT.mkdir(parents=True, exist_ok=True)
    h = H1(); res = {}
    f = OUT / 'results.json'
    if f.exists(): res = json.loads(f.read_text())
    for p in a.parts:
        print(f'===== {p.upper()}', flush=True); {'s1': s1, 's2': s2, 's3': s3, 's4': s4, 's4b': s4b}[p](h, res)
        f.write_text(json.dumps(res, indent=1, default=float))

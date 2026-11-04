"""Step 15a/15b (Jose 2026-10-30; continues results/stage1/2026-10-28_scorer_space/; outputs in
results/stage1/2026-10-29_step15/). Noise-free proxy only (scorer_proxy.Ref.proxy, the Step 9 proxy for all six
components), no cells, run on vcc2027. Inputs: the scorer-space fold caches (cache/fold_<line>.npz, s11_<line>.npz).

Config 1' (F0): called genes 0.5 x C+G LFC; non-called 0.6 x c x Delta_pool, c = z_a^2 / (z_a^2 + 1), single-source
c1 = 1, no source -> control. The magnitude split is always C+G's calls (K0); a call variant changes only the call set.
  K0       C+G's calls
  K4_1200  per knockdown the top 1,200 tested genes by C+G z = |LFC| / SE (scorer_step16.call_sets)
H1 halves: scorer_step16.halves (seed 2026, pairs by rank of true DE count). Halves are scored with the full panel's
scaling, as boot_score does.

  step15a                      config 1' (+ config 1, C+G for reference) x {K0, K4_1200}: H1 full + halves, K562,
                               HCT116, HEK293T (fold caches, held out as before)
  oracle --part P --nparts N   15b.2 on H1: O1 (s_k per true-strength tercile, {0, .25, .5, 1, 1.5}^3, K4_1200 calls),
                               O2 (N_k = clip(lambda x true DE count, 0, 4000), lambda {0.5, 1, 2}, config 1'
                               magnitudes), O1 + O2; rows split over N worker processes
  predictors                   15b.3: Spearman of predictors (a)-(f) with the true DE count on H1 (full, halves) and K562
  predrun --pred X --part P --nparts N   15b.4: O1 + O2 with the predictor X in place of true strength
  o3 [--pred X] --part P --nparts N   addendum O3: N_k = round(sum_g P(DE | C+G z)), the z -> true-DE-rate map
                               fitted on one H1 half, applied to the other; ranked by C+G z; O3 alone and O1 + O3. With
                               --pred: O1 + O3 on X's terciles, and O3s (the map fitted within each X tercile)
  fit                          fits on one half / scores on the other, both ways -> tables/step15b.json
"""
import argparse
import json
import re
import sys
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
import scorer_study as SS                                                                  # noqa: E402
import scorer_proxy as SP                                                                  # noqa: E402
import scorer_step11 as S11                                                                # noqa: E402
import scorer_step12 as S12                                                                # noqa: E402
import scorer_step16 as S16                                                                # noqa: E402

OUT = V / 'results/stage1/2026-10-29_step15'
TAB = OUT / 'tables'
COMPS = S11.COMPS
CONFIGS = {'config1p': dict(m=0.6, kappa=1.0, s_c=0.5, c1=1.0),        # config 1'
           'config1': dict(m=0.5, kappa=1.0, s_c=0.5, c1=1.0)}         # frozen config 1, reference only
S_GRID = [0.0, 0.25, 0.5, 1.0, 1.5]
LAMBDA = [0.5, 1.0, 2.0]
N_MAX = 4000
TERC = ('weak', 'mid', 'strong')


def row(line, s, **keys):
    r = dict(line=line, **keys, overall=s['overall'], mse_unclamped=s['mse_unclamped'])
    for m in SP.SIX:
        r[SP.SHORT[m]] = s['scaled'][m]; r[f'{SP.SHORT[m]}_raw'] = s['raw'][m]
    return r


def cgz_order(G):
    """Tested genes (gated, measured, LFC != 0, own gene out) and their order by C+G z, as scorer_step16.call_sets."""
    import stage1 as S
    fo, h = G.fo, G.fo.h
    ok = fo.gate[None, :] & ~fo.unmeas & (fo.lfc != 0)
    for i in h.rows:
        if h.own[i] >= 0: ok[i, h.own[i]] = False
    zc = np.where(ok, np.abs(fo.lfc) / np.maximum(fo.se, S.SE_MIN), -1.0)
    return ok, np.argsort(-zc, axis=1, kind='stable')


def top_n(G, ok, order, N):
    """Per knockdown the top N[i] tested genes by C+G z (N scalar or per knockdown)."""
    N = np.broadcast_to(np.asarray(N), (G.fo.h.T,))
    sig = np.zeros_like(G.fo.sig)
    for i in G.fo.h.rows:
        n = int(min(N[i], ok[i].sum())); sig[i, order[i, :n]] = True
    return sig


def mu_scaled(G, p, s_k):
    """Config-family means with every knockdown's called and non-called magnitudes x s_k."""
    s = np.asarray(s_k, float)[:, None]
    return G.mu(p['s_c'] * s, 1.0, p['m'] * s, G.c_kappa(p['kappa'], p['c1'], 0.0))


def terciles(x):
    """0 weak / 1 mid / 2 strong by rank of x over the panel."""
    r = pd.Series(x).rank(method='first').to_numpy(); n = len(r)
    return np.where(r <= n / 3, 0, np.where(r <= 2 * n / 3, 1, 2))


# ---------------------------------------------------------------- 15a
def step15a(a):
    TAB.mkdir(parents=True, exist_ok=True); rows = []
    for line in ('H1', 'K562', 'HCT116', 'HEK293T'):
        G = S16.Mag(line); T = S12.truth_calls(G.ref)
        ok, order = cgz_order(G)
        calls = {'K0': G.fo.sig, 'K4_1200': top_n(G, ok, order, 1200)}
        assert (calls['K4_1200'] == S16.call_sets(G, [1200])['K4_1200']).all()
        parts = {'full': np.arange(G.fo.h.T)}
        if line == 'H1':
            hs = S16.halves(T); parts.update(half0=hs[0], half1=hs[1])
        mus = {c: G.build('i', p) for c, p in CONFIGS.items()}
        mus['CG'] = G.fo.mu_lfc(G.fo.lfc)
        for cfg, mu in mus.items():
            for v, sig in calls.items():
                if cfg == 'CG' and v != 'K0':
                    continue
                out = G.ref.proxy(mu=mu, sig=sig, per_kd=True)
                for part, rr in parts.items():
                    s = G.ref.summarise(out['per_kd'], rows=rr)
                    rows.append(row(line, s, config=cfg, calls=v, part=part, calls_median=float(np.median(sig[rr].sum(1))),
                                    truth_calls_median=float(np.median(T[rr].sum(1)))))
                    print(json.dumps({k: rows[-1][k] for k in ['line', 'config', 'calls', 'part', 'overall'] + COMPS},
                                     default=float), flush=True)
        pd.DataFrame(rows).to_csv(TAB / 'step15a.csv', index=False)
    # config 1 must reproduce its scorer-space rows (Step 16 B, K0 / K4_1200 on the H1 halves)
    t = pd.DataFrame(rows); ref = pd.read_parquet(SS.OUT / 'tables/step16_calls_rows_config1.parquet')
    for v in ('K0', 'K4_1200'):
        for h in (0, 1):
            x = t[(t.line == 'H1') & (t.config == 'config1') & (t.calls == v) & (t.part == f'half{h}')].overall.iloc[0]
            y = ref[(ref.line == 'H1') & (ref.variant == v) & (ref.half == h)].overall.iloc[0]
            print(f'config 1 {v} half {h} reproduction |diff| {abs(x - y):.2e}', flush=True); assert abs(x - y) < 1e-9


# ---------------------------------------------------------------- 15b.2 / 15b.4 grids
def grid_rows(G, strength_terc, n_de, label, part, nparts):
    """O1, O2, O1+O2 rows on H1 for a strength tercile assignment and a DE-count estimate n_de (per knockdown)."""
    p = CONFIGS['config1p']; ok, order = cgz_order(G)
    base_sig = top_n(G, ok, order, 1200)
    hs = S16.halves(S12.truth_calls(G.ref))
    jobs = [('base', None, None)] + [('O2', None, l) for l in LAMBDA] + [('O1', s, None) for s in product(S_GRID, repeat=3)] \
        + [('O12', s, l) for s in product(S_GRID, repeat=3) for l in LAMBDA]
    jobs = jobs[part::nparts]
    f_out = TAB / f'grid_{label}_p{part}.parquet'
    rows = pd.read_parquet(f_out).to_dict('records') if f_out.exists() else []
    done = {r['key'] for r in rows}
    sigs = {l: top_n(G, ok, order, np.clip(np.round(l * n_de), 0, N_MAX)) for l in LAMBDA}
    for arm, s, l in jobs:
        key = f'{arm}|{s}|{l}'
        if key in done:
            continue
        s_k = np.ones(G.fo.h.T) if s is None else np.array(s)[strength_terc]
        mu = mu_scaled(G, p, s_k); sig = base_sig if l is None else sigs[l]
        out = G.ref.proxy(mu=mu, sig=sig, per_kd=True)
        for hi, rr in enumerate(hs):
            sm = G.ref.summarise(out['per_kd'], rows=rr)
            rows.append(row('H1', sm, label=label, key=key, arm=arm, half=hi, lam=np.nan if l is None else l,
                            **{f's_{t}': (np.nan if s is None else s[j]) for j, t in enumerate(TERC)},
                            calls_median=float(np.median(sig[rr].sum(1)))))
        print(json.dumps({k: rows[-1][k] for k in ['key', 'half', 'overall'] + COMPS}, default=float), flush=True)
        pd.DataFrame(rows).to_parquet(f_out)


def oracle(a):
    G = S16.Mag('H1'); T = S12.truth_calls(G.ref); n = T.sum(1).astype(float)
    grid_rows(G, terciles(n), n, 'oracle', a.part, a.nparts)


# ---------------------------------------------------------------- O3 (addendum): the model's own expected DE count
def calibrate(z, y, nbins=200):
    """Monotone map C+G z -> observed true-DE rate: quantile bins of z, then isotonic in the bins' median z, weighted by
    the bin size; linear between bin medians, flat past the ends."""
    from sklearn.isotonic import IsotonicRegression
    q = np.unique(np.quantile(z, np.linspace(0, 1, nbins + 1)))
    b = np.clip(np.searchsorted(q, z, side='right') - 1, 0, len(q) - 2)
    cnt = np.bincount(b, minlength=len(q) - 1); ok = cnt > 0
    xm = np.array([np.median(z[b == i]) for i in np.flatnonzero(ok)])
    rate = np.bincount(b, weights=y.astype(float), minlength=len(q) - 1)[ok] / cnt[ok]
    ir = IsotonicRegression(increasing=True, y_min=0.0, y_max=1.0).fit(xm, rate, sample_weight=cnt[ok])
    ky = ir.predict(xm)
    return (lambda v: np.interp(v, xm, ky)), dict(z=xm.tolist(), rate=ky.tolist())


def o3_counts(G, T, fit_rows, strata=None):
    """N_k = round(sum over tested genes of P(DE | C+G z)), the map fitted on the knockdowns fit_rows (pooled, or per
    stratum when strata gives each knockdown a stratum) and applied to every knockdown."""
    import stage1 as S
    ok, _ = cgz_order(G); fo = G.fo
    zc = np.abs(fo.lfc) / np.maximum(fo.se, S.SE_MIN)
    strata = np.zeros(fo.h.T, int) if strata is None else np.asarray(strata)
    N = np.zeros(fo.h.T); maps = {}
    for s_ in np.unique(strata):
        fr = np.intersect1d(fit_rows, np.flatnonzero(strata == s_)); m = np.zeros_like(ok); m[fr] = ok[fr]
        f_, knots = calibrate(zc[m], T[m]); maps[int(s_)] = knots
        rows = np.flatnonzero(strata == s_)
        N[rows] = np.round(np.where(ok[rows], f_(zc[rows]), 0.0).sum(1))
    return N, maps


def o3_rows(G, terc, N, label, part, nparts, single=True):
    """Rows for one fit direction: base (config 1' + K4 1,200), O3 (config 1' magnitudes, top N_k by C+G z) when single,
    and O1 + O3 (s_k per tercile of terc)."""
    p = CONFIGS['config1p']; ok, order = cgz_order(G)
    sigs = {'K4_1200': top_n(G, ok, order, 1200), 'O3': top_n(G, ok, order, np.clip(N, 0, N_MAX))}
    hs = S16.halves(S12.truth_calls(G.ref))
    jobs = [('base', None, 'K4_1200')] + ([('O3', None, 'O3')] if single else []) + \
        [('O13', s, 'O3') for s in product(S_GRID, repeat=3)]
    jobs = jobs[part::nparts]
    f_out = TAB / f'grid_{label}_p{part}.parquet'
    rows = pd.read_parquet(f_out).to_dict('records') if f_out.exists() else []
    done = {r['key'] for r in rows}
    for arm, s, cs in jobs:
        key = f'{arm}|{s}|{cs}'
        if key in done:
            continue
        s_k = np.ones(G.fo.h.T) if s is None else np.array(s)[terc]
        out = G.ref.proxy(mu=mu_scaled(G, p, s_k), sig=sigs[cs], per_kd=True)
        for hi, rr in enumerate(hs):
            sm = G.ref.summarise(out['per_kd'], rows=rr)
            rows.append(row('H1', sm, label=label, key=key, arm=arm, half=hi, lam=np.nan,
                            **{f's_{t}': (np.nan if s is None else s[j]) for j, t in enumerate(TERC)},
                            calls_median=float(np.median(sigs[cs][rr].sum(1)))))
        print(json.dumps({k: rows[-1][k] for k in ['key', 'half', 'overall'] + COMPS}, default=float), flush=True)
        pd.DataFrame(rows).to_parquet(f_out)


def o3(a):
    """O3 per fit half: oracle (true-strength terciles for O1 + O3); with --pred X, (i) O1 + O3 with X's terciles and
    (ii) O3s: the z map fitted within each X tercile (O3s alone and O1 + O3s with X's terciles)."""
    G = S16.Mag('H1'); T = S12.truth_calls(G.ref); hs = S16.halves(T); n = T.sum(1).astype(float)
    info = {}
    for fit_h in (0, 1):
        if not a.pred:
            N, maps = o3_counts(G, T, hs[fit_h])
            o3_rows(G, terciles(n), N, f'o3_fit{fit_h}', a.part, a.nparts)
            info[f'fit{fit_h}'] = dict(N=N.tolist(), maps=maps)
        else:
            h1 = pd.read_csv(TAB / 'step15b_predictors_H1.csv'); assert list(h1.target) == G.fo.h.targets
            x = h1[a.pred].to_numpy(float); x = np.where(np.isfinite(x), x, np.nanmedian(x)); tx = terciles(x)
            N, _ = o3_counts(G, T, hs[fit_h])
            o3_rows(G, tx, N, f'o3p_{a.pred}_fit{fit_h}', a.part, a.nparts, single=False)
            Ns, maps = o3_counts(G, T, hs[fit_h], strata=tx)
            o3_rows(G, tx, Ns, f'o3s_{a.pred}_fit{fit_h}', a.part, a.nparts)
            info[f'fit{fit_h}'] = dict(N_strat=Ns.tolist(), maps_strat=maps)
    if a.part == 0:
        lab = f'_{a.pred}' if a.pred else ''
        (TAB / f'step15b_o3{lab}_counts.json').write_text(json.dumps(dict(targets=G.fo.h.targets, true_de=n.tolist(),
                                                                          halves=[h.tolist() for h in hs], **info)) + '\n')


# ---------------------------------------------------------------- 15b.3 predictors
HGNC_RULES = {   # 15b.3 (e): HGNC gene-group names matched (case-insensitive)
    'chromatin_remodelling': r'chromatin-remodelling|chromatin accessibility complex|^BAF complex|^PBAF complex|'
                             r'^GBAF complex|^INO80 complex|^NuRD complex|^SRCAP complex',
    'transcription_factor': r'(?<!General )transcription factor|homeobox|^Basic helix-loop-helix proteins|bZIP|'
                            r'^Basic leucine zipper|^Forkhead boxes|^Nuclear receptor subfamily|^SRY-box',
    'mediator': r'^Mediator complex',
    'TFIID': r'^General transcription factor IID complex|^TATA-box binding protein',
    'm6A': r'^m6A methyltransferase complex',
}


def hgnc_flags(targets):
    grp = S11.hgnc_groups(); allg = sorted({g for gs in grp.values() for g in gs})
    matched = {c: [g for g in allg if re.search(rx, g, re.I)] for c, rx in HGNC_RULES.items()}
    flags = {c: np.array([any(g in set(ms) for g in grp.get(t, [])) for t in targets]) for c, ms in matched.items()}
    return flags, matched


def source_n_real(line, targets):
    """(b): per knockdown, the mean over source tables (the line's own background excluded) of each table's mean n_real
    over its rows for that knockdown; tables with de_labels_known False or no finite n_real are skipped."""
    import h5py
    import stage1 as S
    from anndata.io import read_elem
    key = SS.LINE_KEY[line]; per = {t: [] for t in targets}; used = []
    for p in S.discover_tables([str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')]):
        with h5py.File(p, 'r') as f:
            o = read_elem(f['obs'])
        if S._line_key(o['background'].iloc[0]) == key:
            continue
        o = o[o['de_labels_known'].astype(bool) & np.isfinite(o['n_real'].astype(float))]
        if not len(o):
            continue
        used.append(p.stem)
        m = o.groupby(o['target_gene'].astype(str))['n_real'].mean()
        for t in targets:
            if t in m.index: per[t].append(float(m[t]))
    return np.array([np.mean(per[t]) if per[t] else np.nan for t in targets]), used


def predictor_table(line):
    G = S16.Mag(line); fo = G.fo; T = S12.truth_calls(G.ref); targets = fo.h.targets
    w = np.where(fo.ok, fo.J2[None, :], 0.0)
    a_ = (w * fo.dpool ** 2).sum(1)                                                       # (a) as Step 11 S1
    b_, used = source_n_real(line, targets)                                               # (b)
    zok = fo.ok & G.has2                                                                  # (c)
    q90 = float(np.quantile(G.za[zok], 0.9))
    c_ = np.where(zok.sum(1) > 0, (zok & (np.nan_to_num(G.za) > q90)).sum(1) / np.maximum(zok.sum(1), 1), np.nan)
    gp = {g: j for j, g in enumerate(fo.h.genes)}                                         # (d) target CPM in the controls
    d_ = np.array([fo.m0[gp[t]] if t in gp else np.nan for t in targets])
    flags, matched = hgnc_flags(targets)                                                  # (e)
    e_ = np.any(np.stack(list(flags.values())), 0).astype(float)
    df = pd.DataFrame(dict(target=targets, true_de=T.sum(1), a=a_, b=b_, c=c_, d=d_, e=e_,
                           **{f'e_{k}': v.astype(int) for k, v in flags.items()}))
    meta = dict(line=line, z_agree_q90=q90, n_real_tables=used, hgnc_groups_matched=matched,
                hgnc_counts={k: int(v.sum()) for k, v in flags.items()})
    return df, meta


FEATS = ['a', 'b', 'c', 'd', 'e']


def design(df, mu=None, sd=None, fill=None):
    X = np.column_stack([np.log1p(df.a), np.log1p(df.b), df.c, np.log1p(df.d), df.e]).astype(float)
    fill = np.nanmedian(X, 0) if fill is None else fill
    X = np.where(np.isfinite(X), X, fill)
    mu = X.mean(0) if mu is None else mu; sd = X.std(0) if sd is None else sd
    return (X - mu) / np.where(sd > 0, sd, 1.0), mu, sd, fill


def ridge_fit(df, alpha=1.0):
    X, mu, sd, fill = design(df); y = np.log1p(df.true_de.to_numpy(float)); ym = y.mean()
    beta = np.linalg.solve(X.T @ X + alpha * np.eye(X.shape[1]), X.T @ (y - ym))
    return lambda d: design(d, mu, sd, fill)[0] @ beta + ym, dict(beta=dict(zip(FEATS, beta.tolist())), alpha=alpha)


def predictors(a):
    TAB.mkdir(parents=True, exist_ok=True)
    h1, m1 = predictor_table('H1'); k5, m5 = predictor_table('K562')
    G = S16.Mag('H1'); hs = S16.halves(S12.truth_calls(G.ref)); del G
    h1['half'] = -1
    for hi, rr in enumerate(hs): h1.loc[rr, 'half'] = hi
    res = dict(meta=dict(H1=m1, K562=m5), rho={}, ridge={})
    f_held = np.full(len(h1), np.nan)
    for hi in (0, 1):
        fit, info = ridge_fit(h1[h1.half == hi]); sc = h1.half == 1 - hi
        f_held[sc.to_numpy()] = fit(h1[sc]); res['ridge'][f'fit_half{hi}'] = info
    h1['f'] = f_held                                    # each knockdown's ridge score from the half it was not fitted on
    fit_all, info_all = ridge_fit(h1); k5['f'] = fit_all(k5); res['ridge']['fit_H1_all_for_K562'] = info_all
    for x in FEATS + ['f']:
        r = dict(H1_full=S11.spearman(h1[x].to_numpy(float), h1.true_de.to_numpy(float)),
                 **{f'H1_half{hi}': S11.spearman(h1[x][h1.half == hi].to_numpy(float),
                                                h1.true_de[h1.half == hi].to_numpy(float)) for hi in (0, 1)},
                 K562=S11.spearman(k5[x].to_numpy(float), k5.true_de.to_numpy(float)),
                 H1_n_finite=int(np.isfinite(h1[x]).sum()), K562_n_finite=int(np.isfinite(k5[x]).sum()))
        r['H1_heldout_mean'] = float(np.mean([r['H1_half0'], r['H1_half1']]))
        res['rho'][x] = r; print(x, json.dumps(r, default=float), flush=True)
    for k in [c for c in h1.columns if c.startswith('e_')]:
        res['rho'][k] = dict(H1_full=S11.spearman(h1[k].to_numpy(float), h1.true_de.to_numpy(float)),
                             K562=S11.spearman(k5[k].to_numpy(float), k5.true_de.to_numpy(float)), H1_n=int(h1[k].sum()))
    best = max(FEATS + ['f'], key=lambda x: res['rho'][x]['H1_heldout_mean']); res['best'] = best
    h1.to_csv(TAB / 'step15b_predictors_H1.csv', index=False); k5.to_csv(TAB / 'step15b_predictors_K562.csv', index=False)
    (TAB / 'step15b_predictors.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print('best', best, flush=True)


def predrun(a):
    """15b.4: the predictor in place of true strength. Terciles by the predictor's rank (no truth). DE count: the
    predictor's panel percentile mapped onto the fit half's true DE-count quantiles, so N_k depends on the half the
    parameters are fitted on: rows are built for each fit half and scored on the other."""
    h1 = pd.read_csv(TAB / 'step15b_predictors_H1.csv'); x = h1[a.pred].to_numpy(float)
    x = np.where(np.isfinite(x), x, np.nanmedian(x)); G = S16.Mag('H1')
    assert list(h1.target) == G.fo.h.targets
    pct = (pd.Series(x).rank(method='average').to_numpy() - 0.5) / len(x)
    for fit_h in (0, 1):
        tr = h1.true_de[h1.half == fit_h].to_numpy(float)
        n_de = np.quantile(tr, pct)
        grid_rows(G, terciles(x), n_de, f'pred_{a.pred}_fit{fit_h}', a.part, a.nparts)


def lamext(a):
    """Side check past the O2 grid edge (not part of the task's fits): O2 at lambda {3, 4, 6} with the true count
    (oracle) and with predictor X's mapped count (per fit half), config 1' magnitudes. Rows -> tables/side_lambda.parquet."""
    G = S16.Mag('H1'); T = S12.truth_calls(G.ref); hs = S16.halves(T); n = T.sum(1).astype(float)
    ok, order = cgz_order(G); p = CONFIGS['config1p']; mu = mu_scaled(G, p, np.ones(G.fo.h.T))
    h1 = pd.read_csv(TAB / 'step15b_predictors_H1.csv'); x = h1[a.pred].to_numpy(float)
    x = np.where(np.isfinite(x), x, np.nanmedian(x)); pct = (pd.Series(x).rank(method='average').to_numpy() - 0.5) / len(x)
    counts = {'oracle': n, **{f'pred_{a.pred}_fit{h}': np.quantile(n[hs[h]], pct) for h in (0, 1)}}
    rows = []
    for lab, c in counts.items():
        for l in (3.0, 4.0, 6.0):
            sig = top_n(G, ok, order, np.clip(np.round(l * c), 0, N_MAX))
            out = G.ref.proxy(mu=mu, sig=sig, per_kd=True)
            for hi, rr in enumerate(hs):
                rows.append(row('H1', G.ref.summarise(out['per_kd'], rows=rr), label=lab, arm='O2', lam=l, half=hi,
                                calls_median=float(np.median(sig[rr].sum(1)))))
                print(json.dumps({k: rows[-1][k] for k in ['label', 'lam', 'half', 'overall', 'calls_median']}, default=float), flush=True)
    pd.DataFrame(rows).to_parquet(TAB / 'side_lambda.parquet')


# ---------------------------------------------------------------- fits
def fit_grid(t, fit_h):
    """Per arm, the key with the best overall on fit half fit_h, and its row on the other half."""
    out = {}
    for arm in sorted(set(t.arm) - {'base'}):
        sub = t[(t.arm == arm) & (t.half == fit_h)]
        k = sub.sort_values('overall').iloc[-1].key
        out[arm] = t[(t.key == k) & (t.half == 1 - fit_h)].iloc[0].to_dict()
    out['base'] = t[(t.arm == 'base') & (t.half == 1 - fit_h)].iloc[0].to_dict()
    return out


def summarise_fit(label_of_fit):
    """label_of_fit(fit_h) -> grid label. Held-out rows per fit direction, means and gains over the base."""
    res = {}
    for fit_h in (0, 1):
        lab = label_of_fit(fit_h)
        t = pd.concat([pd.read_parquet(p) for p in sorted(TAB.glob(f'grid_{lab}_p*.parquet'))])
        res[f'fit_half{fit_h}'] = fit_grid(t, fit_h)
    arms = [k for k in res['fit_half0'] if k in res['fit_half1']]
    mean = {arm: {k: float(np.mean([res[f'fit_half{h}'][arm][k] for h in (0, 1)])) for k in ['overall'] + COMPS}
            for arm in arms}
    gain = {arm: {k: mean[arm][k] - mean['base'][k] for k in ['overall'] + COMPS} for arm in arms if arm != 'base'}
    return dict(per_fit=res, heldout_mean=mean, gain_vs_base=gain)


def fit(a):
    out = dict(oracle=summarise_fit(lambda h: 'oracle'))
    if list(TAB.glob('grid_o3_fit0_p*.parquet')):
        out['oracle_o3'] = summarise_fit(lambda h: f'o3_fit{h}')
    pj = TAB / 'step15b_predictors.json'
    if pj.exists():
        best = json.loads(pj.read_text())['best']
        og = {**out['oracle']['gain_vs_base'], **out.get('oracle_o3', {}).get('gain_vs_base', {})}
        for lab in (f'pred_{best}', f'o3p_{best}', f'o3s_{best}'):
            if list(TAB.glob(f'grid_{lab}_fit0_p*.parquet')):
                pr = summarise_fit(lambda h: f'{lab}_fit{h}')
                pr['share_of_oracle_gain'] = {arm: (pr['gain_vs_base'][arm]['overall'] / og[arm]['overall']
                                                    if og.get(arm) and og[arm]['overall'] else None)
                                              for arm in pr['gain_vs_base']}
                out[lab] = pr
    (TAB / 'step15b.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(json.dumps({k: dict(heldout_mean=v['heldout_mean'], gain=v['gain_vs_base']) for k, v in out.items()},
                     indent=1, default=float))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['step15a', 'oracle', 'o3', 'predictors', 'predrun', 'lamext', 'fit'])
    p.add_argument('--part', type=int, default=0)
    p.add_argument('--nparts', type=int, default=1)
    p.add_argument('--pred', default='')
    a = p.parse_args()
    {'step15a': step15a, 'oracle': oracle, 'predictors': predictors, 'predrun': predrun, 'o3': o3, 'lamext': lamext, 'fit': fit}[a.cmd](a)

"""Scorer-space study part 2, Step 16 (Jose 2026-10-29, results/stage1/2026-10-28_scorer_space/): noise-free proxy on
vcc2027, no cells. Base: config 1 (F0 / arm i / K0).

A, magnitude (all four lines; LOFO over the three folds as Steps 10-11):
  i   called genes s_c x LFC; non-called m x c x Delta_pool, c = z_a^2 / (z_a^2 + kappa) with the agreement z (>= 2
      sources), c1 for one source, 0 for none
  j   as i, and the called genes also x c (c1 for <= 1 source)
  k   as i, c replaced by a monotone curve of z_a fixed at 0: the per-decile beta of agreement z on the two training
      folds (each line normalised to its top decile, folds weighted equally), isotonic in z, 0 at z = 0, 1 past the
      top decile. Every line is scored with each held-out fold's curve, so a fold's fit never sees its own beta.
  kj  the curve on both sides, as j
  grid: m {0.5, 0.6, 0.75}, kappa {0.25, 0.5, 1, 2} (i, j), s_c {0.35, 0.5, 0.65, 0.8}, c1 {0.5, 1}
B, calls (on A's carried config, or with --magnitude config1 on config 1; the magnitudes stay, only the calls change):
  H1's knockdowns split in two halves (seed 2026, balanced by true DE count: pairs by rank, one of each pair per half);
  fit on one half, score the other, swap. Halves are scored with the full panel's scaling, as boot_score does.
  K0  C+G's calls
  K4  per knockdown the top N tested genes by C+G z = |LFC| / SE
  K5  the top N by agreement z among genes with >= 2 sources, then the rest by C+G z
  N in {400, 600, 800, 1200}; every variant also scored on all of HCT116 and HEK293T (guardrails).

  step16 --line L | step16fit | step16calls
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
import scorer_study as SS                                                                  # noqa: E402
import scorer_proxy as SP                                                                  # noqa: E402
import scorer_part2 as P2                                                                  # noqa: E402
import scorer_step11 as S11                                                                # noqa: E402
import scorer_step12 as S12                                                                # noqa: E402

OUT = SS.OUT
COMPS = S11.COMPS
M, KAPPA, SC, C1 = [0.5, 0.6, 0.75], [0.25, 0.5, 1.0, 2.0], [0.35, 0.5, 0.65, 0.8], [0.5, 1.0]
GRID = dict(m=M, kappa=KAPPA, s_c=SC, c1=C1)
KEYS = {'i': ['m', 'kappa', 's_c', 'c1'], 'j': ['m', 'kappa', 's_c', 'c1'], 'k': ['m', 's_c', 'c1'], 'kj': ['m', 's_c', 'c1']}
BUDGETS = [400, 600, 800, 1200]


def fold_of(line):
    return next(f_ for f_, fl in SS.FOLD_LINES.items() if line in fl)


class Mag:
    """Config-1 family magnitudes for one line (F0)."""

    def __init__(self, line):
        self.line = line; self.ref = P2.refs_for(line); self.fo = P2.Fold(line, 'F0', self.ref)
        s = np.load(OUT / 'cache' / f's11_{line}.npz', allow_pickle=True)
        self.n = s['n_bg']; self.za = s['z_agree'].astype(np.float64)
        self.has2 = (self.n >= 2) & np.isfinite(self.za)

    def c_kappa(self, kappa, c1, none):
        z2 = np.nan_to_num(self.za) ** 2
        return np.where(self.has2, z2 / (z2 + kappa), np.where(self.n == 1, c1, none))

    def c_curve(self, curve, c1, none):
        return np.where(self.has2, curve(np.nan_to_num(self.za)), np.where(self.n == 1, c1, none))

    def mu(self, s_c, c_called, m, c_non):
        fo = self.fo
        return np.where(fo.sig, fo.mu_lfc(s_c * c_called * fo.lfc), fo.mu_delta(np.where(fo.dn > 0, m * c_non * fo.dpool, 0.0)))

    def build(self, arm, p, curve=None):
        if arm in ('i', 'j'):
            c_non = self.c_kappa(p['kappa'], p['c1'], 0.0)
            c_cal = self.c_kappa(p['kappa'], p['c1'], p['c1']) if arm == 'j' else 1.0
        else:
            c_non = self.c_curve(curve, p['c1'], 0.0)
            c_cal = self.c_curve(curve, p['c1'], p['c1']) if arm == 'kj' else 1.0
        return self.mu(p['s_c'], c_cal, p['m'], c_non)

    def zbeta(self):
        """Per-decile (median agreement z, beta) of true on pooled Delta, gated non-called genes (Step 11 S4)."""
        fo = self.fo; ok = fo.ok & ~fo.sig & (fo.dn > 0) & self.has2
        x, y, z = fo.dpool[ok], fo.true_delta[ok], self.za[ok]
        q = np.quantile(z, np.linspace(0, 1, 11)); d = np.clip(np.searchsorted(q, z, side='right') - 1, 0, 9)
        return [dict(decile=k, z_med=float(np.median(z[d == k])), beta=float((x[d == k] * y[d == k]).sum() / (x[d == k] ** 2).sum()),
                     n=int((d == k).sum())) for k in range(10)]


def zbeta_file(line):
    return OUT / 'tables' / f'step16_zbeta_{line}.json'


def curve_for(held_fold, wait=True):
    """The monotone c(z) fitted on the folds other than held_fold, plus its knots."""
    from sklearn.isotonic import IsotonicRegression
    train = [f_ for f_ in SS.FOLD_LINES if f_ != held_fold]
    xs, ys, ws = [], [], []
    for f_ in train:
        for l in SS.FOLD_LINES[f_]:
            while wait and not zbeta_file(l).exists():
                time.sleep(20)
            zb = json.loads(zbeta_file(l).read_text()); top = zb[-1]['beta']
            for r in zb:
                xs.append(r['z_med']); ys.append(r['beta'] / top); ws.append(1.0 / (len(train) * len(SS.FOLD_LINES[f_]) * len(zb)))
    o = np.argsort(xs); xs, ys, ws = np.array(xs)[o], np.array(ys)[o], np.array(ws)[o]
    ir = IsotonicRegression(increasing=True, y_min=0.0).fit(xs, ys, sample_weight=ws)
    kx = np.unique(xs); ky = ir.predict(kx); ky = ky / ky[-1]
    kx, ky = np.concatenate([[0.0], kx]), np.concatenate([[0.0], ky])
    return (lambda z: np.interp(z, kx, ky, right=1.0)), dict(z=kx.tolist(), c=ky.tolist())


def step16(a):
    line = a.line; G = Mag(line); fo = G.fo
    zbeta_file(line).write_text(json.dumps(G.zbeta()) + '\n')
    f_out = OUT / 'tables' / f'step16_rows_{line}.parquet'
    rows = pd.read_parquet(f_out).to_dict('records') if f_out.exists() else []
    done = {r['key'] for r in rows}

    def add(arm, key, mu, **p):
        if key in done: return
        rows.append(dict(line=line, arm=arm, key=key, **p, **fo.score(mu, fo.sig)))
        print(json.dumps({k: rows[-1][k] for k in ['arm', 'key', 'overall'] + COMPS}, default=float), flush=True)
        pd.DataFrame(rows).to_parquet(f_out)

    for arm in ('i', 'j'):
        for m in M:
            for kappa in KAPPA:
                for s_c in SC:
                    for c1 in C1:
                        p = dict(m=m, kappa=kappa, s_c=s_c, c1=c1)
                        add(arm, f'{arm}_{m}_{kappa}_{s_c}_{c1}', G.build(arm, p), curve='', **p)
    # Step 11's config 1 (m 0.5, kappa 1, called 0.5, c1 1) must reproduce its Step 11 row
    s11 = pd.read_parquet(OUT / 'tables' / f'step11_rows_{line}.parquet').set_index('key')
    t = pd.DataFrame(rows).set_index('key')
    d = abs(t.loc['i_0.5_1.0_0.5_1.0', 'overall'] - s11.loc['i0.5_1.0_1.0', 'overall'])
    print(f'config-1 reproduction |diff| {d:.2e}', flush=True); assert d < 1e-9, d
    for held in SS.FOLD_LINES:
        curve, knots = curve_for(held)
        (OUT / 'tables' / f'step16_curve_{held}.json').write_text(json.dumps(knots) + '\n')
        for arm in ('k', 'kj'):
            for m in M:
                for s_c in SC:
                    for c1 in C1:
                        p = dict(m=m, s_c=s_c, c1=c1)
                        add(arm, f'{arm}_{held}_{m}_{s_c}_{c1}', G.build(arm, p, curve), curve=held, kappa=np.nan, **p)


def fit_arm(t, arm):
    """Per held-out fold, the keys maximising the training folds' mean expected overall (k/kj: rows on that fold's curve)."""
    keys = KEYS[arm]; out = {}
    for f_ in SS.FOLD_LINES:
        sub = t[(t.arm == arm) & ((t.curve == f_) if arm in ('k', 'kj') else True)]
        train = [g for g in SS.FOLD_LINES if g != f_]
        best = max(sub.groupby(keys), key=lambda kv: P2.fold_mean(kv[1], train))[0]
        out[f_] = dict(zip(keys, map(float, best)))
    return out


def held_row(t, arm, f_, p, line):
    m = (t.arm == arm) & (t.line == line) & ((t.curve == f_) if arm in ('k', 'kj') else True)
    for k, v in p.items(): m &= np.isclose(t[k], v)
    return t[m].iloc[0].to_dict()


def spreads(fits, arm):
    out = {}
    for k in KEYS[arm]:
        v = [fits[f_][k] for f_ in SS.FOLD_LINES]
        out[k] = dict(folds=dict(zip(SS.FOLD_LINES, v)), max_over_min=(max(v) / min(v)) if min(v) > 0 else None,
                      edge=sorted({x for x in v if x in (min(GRID[k]), max(GRID[k]))}))
    return out


def step16fit(a):
    t = pd.concat([pd.read_parquet(OUT / 'tables' / f'step16_rows_{l}.parquet') for l in P2.LINES])
    held, fits = [], {}
    for arm in KEYS:
        fits[arm] = fit_arm(t, arm)
        for f_, p in fits[arm].items():
            for l in SS.FOLD_LINES[f_]:
                r = held_row(t, arm, f_, p, l); r['arm'] = f'{arm}_lofo'; held.append(r)
    H = pd.DataFrame(held); H.to_csv(OUT / 'tables/step16_heldout.csv', index=False)
    ev = H[H.line == 'H1'].copy(); ev['qualifies'] = (ev.mse > 0) & (ev[COMPS] > 0).all(1)
    q = ev[ev.qualifies]; pick = (q if len(q) else ev).sort_values('overall', ascending=False).iloc[0]
    arm = pick['arm'].replace('_lofo', '')
    insample = {a_: t[(t.line == 'H1') & (t.arm == a_)].sort_values('overall').iloc[-1][['key', 'overall'] + COMPS].to_dict()
                for a_ in KEYS}
    res = dict(fits=fits, spreads={a_: spreads(fits[a_], a_) for a_ in KEYS},
               carried=dict(arm=arm, params=fits[arm], evalA_overall=pick['overall'], qualified=bool(len(q) > 0)),
               evalA=ev[['arm', 'key', 'overall'] + COMPS + ['mse_unclamped', 'r_J2_delta', 'qualifies']].to_dict('records'),
               h1_in_sample_best=insample,
               curves={f_: json.loads((OUT / 'tables' / f'step16_curve_{f_}.json').read_text()) for f_ in SS.FOLD_LINES})
    (OUT / 'tables/step16.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    pd.set_option('display.width', 250)
    print(H[['arm', 'line', 'key', 'overall'] + COMPS + ['r_J2_delta']].round(4).to_string())
    print(json.dumps(res['carried'], default=float)); print(json.dumps(res['spreads'], default=float))


def halves(T, seed=2026):
    """Two halves of the knockdowns, balanced by true DE count (rank pairs, one of each pair per half)."""
    order = np.argsort(T.sum(1), kind='stable'); rng = np.random.default_rng(seed); h = np.zeros(len(order), int)
    for i in range(0, len(order) - 1, 2):
        h[order[i:i + 2]] = rng.permutation(2)
    if len(order) % 2: h[order[-1]] = 0
    return [np.flatnonzero(h == 0), np.flatnonzero(h == 1)]


def call_sets(G, budgets=BUDGETS):
    import stage1 as S
    fo, h = G.fo, G.fo.h
    ok = fo.gate[None, :] & ~fo.unmeas & (fo.lfc != 0)
    for i in h.rows:
        if h.own[i] >= 0: ok[i, h.own[i]] = False
    zc = np.where(ok, np.abs(fo.lfc) / np.maximum(fo.se, S.SE_MIN), -1.0)
    o4 = np.argsort(-zc, axis=1, kind='stable')
    k5 = np.where(ok & G.has2, 1e9 + np.nan_to_num(G.za), zc)                             # >= 2 sources first, by z_a
    o5 = np.argsort(-k5, axis=1, kind='stable')
    out = {'K0': fo.sig}
    for N in budgets:
        for v, o in (('K4', o4), ('K5', o5)):
            sig = np.zeros_like(fo.sig)
            for i in h.rows:
                n = min(N, int(ok[i].sum())); sig[i, o[i, :n]] = True
            out[f'{v}_{N}'] = sig
    return out


CONFIG1 = dict(arm='i', params={f_: dict(m=0.5, kappa=1.0, s_c=0.5, c1=1.0) for f_ in SS.FOLD_LINES})


def step16calls(a):
    car = CONFIG1 if a.magnitude == 'config1' else json.loads((OUT / 'tables/step16.json').read_text())['carried']
    arm = car['arm']; suf = ('_config1' if a.magnitude == 'config1' else '') + ('_ext' if a.budgets else '')
    budgets = [int(x) for x in a.budgets.split(',')] if a.budgets else BUDGETS
    rows = []
    for line in ('H1', 'HCT116', 'HEK293T'):
        G = Mag(line); f_ = fold_of(line); p = car['params'][f_]
        mu = G.build(arm, p, curve_for(f_, wait=False)[0] if arm in ('k', 'kj') else None)
        sets = call_sets(G, budgets); T = S12.truth_calls(G.ref)
        parts = halves(T) if line == 'H1' else [np.arange(G.fo.h.T)]
        if line == 'H1':
            hs_truth = [float(np.median(T[x].sum(1))) for x in parts]
        for v, sig in sets.items():
            out = G.ref.proxy(mu=mu, sig=sig, per_kd=True)
            for hi, rr in enumerate(parts):
                s = G.ref.summarise(out['per_kd'], rows=rr)
                r = dict(line=line, variant=v, half=hi if line == 'H1' else -1, overall=s['overall'],
                         calls_median=float(np.median(sig[rr].sum(1))), truth_calls_median=float(np.median(T[rr].sum(1))))
                for m in SP.SIX:
                    r[SP.SHORT[m]] = s['scaled'][m]; r[f'{SP.SHORT[m]}_raw'] = s['raw'][m]
                rows.append(r)
                print(json.dumps({k: r[k] for k in ['line', 'variant', 'half', 'overall', 'calls_median', 'fidelity', 'reach', 'jaccard']},
                                 default=float), flush=True)
    t = pd.DataFrame(rows); t.to_parquet(OUT / f'tables/step16_calls_rows{suf}.parquet')
    h1 = t[t.line == 'H1']
    fam = {'K4': [f'K4_{N}' for N in budgets], 'K5': [f'K5_{N}' for N in budgets],
           'any': ['K0'] + [f'K{k}_{N}' for k in (4, 5) for N in budgets]}
    fits = {}
    for fname, vs in fam.items():
        res = []
        for fit_h in (0, 1):
            sub = h1[(h1.half == fit_h) & h1.variant.isin(vs)]; best = sub.sort_values('overall').iloc[-1].variant
            sc = h1[(h1.half == 1 - fit_h) & (h1.variant == best)].iloc[0]
            res.append(dict(fit_half=fit_h, chosen=best, **{k: float(sc[k]) for k in ['overall', 'fidelity', 'reach', 'jaccard', 'calls_median']}))
        fits[fname] = dict(swaps=res, heldout_mean={k: float(np.mean([r[k] for r in res])) for k in ['overall', 'fidelity', 'reach', 'jaccard']})
    k0 = h1[h1.variant == 'K0']
    out = dict(carried_magnitude=car, halves_truth_median=hs_truth, fits=fits,
               K0_heldout_mean={k: float(k0[k].mean()) for k in ['overall', 'fidelity', 'reach', 'jaccard']},
               guardrails=t[t.line != 'H1'][['line', 'variant', 'overall', 'calls_median'] + COMPS].to_dict('records'))
    (OUT / f'tables/step16_calls{suf}.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    pd.set_option('display.width', 250)
    print(t[['line', 'variant', 'half', 'calls_median', 'truth_calls_median', 'overall'] + COMPS].round(4).to_string())
    print(json.dumps(fits, default=float))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['step16', 'step16fit', 'step16calls'])
    p.add_argument('--line', choices=P2.LINES)
    p.add_argument('--budgets', default='', help='step16calls: comma-separated N past the grid (side check, writes *_ext)')
    p.add_argument('--magnitude', choices=['carried', 'config1'], default='carried', help='step16calls: Step 16 carry or config 1')
    a = p.parse_args()
    {'step16': step16, 'step16fit': step16fit, 'step16calls': step16calls}[a.cmd](a)

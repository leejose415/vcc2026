"""Step 17 (Jose 2026-10-30): volume-matched call allocation by predictor (b), mean source DE count. Noise-free proxy
only (scorer_proxy.Ref.proxy), no cells, one VM (vcc2027, after its 15c job). Inputs: the scorer-space fold caches and
Step 15's predictor table; helpers from step15.py (unchanged). Outputs in results/stage1/2026-10-30_step17/ (tables,
logs); write-up in the scorer-space report.md.

Base: config 1' magnitudes (step15.CONFIGS['config1p']) for every row; calls are the top N_k tested genes by C+G z.
  b|h|lam        N_k = clip(round(lam x q_h(pct_k)), 0, 4000): pct_k = (b)'s panel percentile, q_h = the true DE-count
                 quantile function of H1 half h (the fit half); lam in {4, 6, 8, 10, 14}. On H1 scored on the other half
                 (lam is fitted on half h, scored on 1 - h); on K562 / HCT116 / HEK293T the guardrails, q_h from H1.
  oracle|lam     N_k = clip(round(lam x true count_k), 0, 4000), lam in {1, 2, 3, 4, 6, 8, 10, 14}; H1 halves.
  bfull|lam      volume-free variant: q = H1's full true DE-count distribution (what D/E/F would borrow); lam {4, 6, 8}.
  flat|N         top N by C+G z for every knockdown (K4 at N): N = the median calls of an allocation row on the part it
                 is scored on, so each allocation row has a volume-matched flat row.

  predictor                  run first: 1c (Spearman of (b) per line, own tables out; H1 cross-fit calibration by
                             decile) and 4 ((b) for the 300 D/E/F targets with all sources); also writes
                             tables/h1_truth.json and (b) for HCT116 / HEK293T, which `rows` reads
  rows --line L --part P --nparts N   proxy rows for one line (all six components and overall per part, per-knockdown
                             values on H1) -> tables/rows_<line>_p<P>of<N>.parquet
  summary                    held-out tables, flat matches, oracle shares, item 3 terciles -> tables/step17.json
  meanmatch                  side check: flat K4 at each allocation row's MEAN volume -> tables/meanmatch.csv
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
import step15 as S15                                                                       # noqa: E402
import scorer_proxy as SP                                                                  # noqa: E402
import scorer_step11 as S11                                                                # noqa: E402
import scorer_step12 as S12                                                                # noqa: E402
import scorer_step16 as S16                                                                # noqa: E402

OUT = V / 'results/stage1/2026-10-30_step17'
TAB = OUT / 'tables'
LAM_B = [4.0, 6.0, 8.0, 10.0, 14.0]
LAM_O = [1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 14.0]
LAM_F = [4.0, 6.0, 8.0]
GUARD = ('K562', 'HCT116', 'HEK293T')
COMPS = S15.COMPS


def pct_of(x):
    x = np.where(np.isfinite(x), x, np.nanmedian(x))
    return (pd.Series(x).rank(method='average').to_numpy() - 0.5) / len(x)


def b_values(line, targets):
    """(b) for a line's knockdowns, the line's own tables out (Step 15's table for H1 / K562; computed otherwise)."""
    f = S15.TAB / f'step15b_predictors_{line}.csv'
    if f.exists():
        t = pd.read_csv(f); assert list(t.target) == list(targets); return t.b.to_numpy(float)
    c = TAB / f'b_{line}.csv'
    if not c.exists():
        b, used = S15.source_n_real(line, targets)
        pd.DataFrame(dict(target=targets, b=b)).to_csv(c, index=False)
    t = pd.read_csv(c); assert list(t.target) == list(targets); return t.b.to_numpy(float)


class Line:
    def __init__(self, line):
        self.line = line; self.G = G = S16.Mag(line); self.T = S12.truth_calls(G.ref); self.n = self.T.sum(1).astype(float)
        self.ok, self.order = S15.cgz_order(G)
        self.mu = S15.mu_scaled(G, S15.CONFIGS['config1p'], np.ones(G.fo.h.T))
        self.targets = G.fo.h.targets; self.b = b_values(line, self.targets); self.pct = pct_of(self.b)
        self.parts = {'full': np.arange(G.fo.h.T)}

    def sig(self, N):
        return S15.top_n(self.G, self.ok, self.order, np.clip(np.round(N), 0, S15.N_MAX))


def h1_halves(L):
    hs = S16.halves(L.T); return hs


def h1_truth():
    """H1's true DE counts and halves (written by `predictor`), so a guardrail process need not load H1."""
    z = json.loads((TAB / 'h1_truth.json').read_text())
    return np.array(z['n'], float), [np.array(h) for h in z['halves']]


def jobs_for(L):
    """(key, N per knockdown, parts to score) for one line; q_h from H1's truth."""
    n1, hs = h1_truth(); J = []
    q = {f'fit{h}': n1[hs[h]] for h in (0, 1)}; q['full'] = n1
    for h in (0, 1):
        for l in LAM_B:
            J.append((f'b|{h}|{l}', l * np.quantile(q[f'fit{h}'], L.pct)))
    for l in LAM_F:
        J.append((f'bfull|{l}', l * np.quantile(q['full'], L.pct)))
    if L.line == 'H1':
        for l in LAM_O:
            J.append((f'oracle|{l}', l * L.n))
    J.append(('flat|1200', np.full(len(L.targets), 1200.0)))
    if L.line == 'H1':
        assert np.array_equal(L.n, n1)
    parts = L.parts if L.line != 'H1' else dict(full=np.arange(len(L.targets)), half0=hs[0], half1=hs[1])
    flats = set()
    for _, N in J:
        s = np.clip(np.round(N), 0, S15.N_MAX)
        cap = L.ok.sum(1)
        for rr in parts.values():
            flats.add(int(np.median(np.minimum(s, cap)[rr])))
    J += [(f'flat|{n}', np.full(len(L.targets), float(n))) for n in sorted(flats) if f'flat|{n}' not in {k for k, _ in J}]
    return J, parts


def rows(a):
    TAB.mkdir(parents=True, exist_ok=True)
    ln = a.line; L = Line(ln); J, parts = jobs_for(L)
    mine = J[a.part::a.nparts]
    tag = f'{ln}_p{a.part}of{a.nparts}'
    f_out = TAB / f'rows_{tag}.parquet'; f_pk = TAB / f'perkd_{tag}.npz'
    rec = pd.read_parquet(f_out).to_dict('records') if f_out.exists() else []
    done = {(r['line'], r['key']) for r in rec}
    pk = dict(np.load(f_pk, allow_pickle=True)) if f_pk.exists() else {}
    for key, N in mine:
        if (ln, key) in done:
            continue
        sig = L.sig(N)
        out = L.G.ref.proxy(mu=L.mu, sig=sig, per_kd=True)
        for part, rr in parts.items():
            s = L.G.ref.summarise(out['per_kd'], rows=rr)
            rec.append(S15.row(ln, s, key=key, part=part, calls_median=float(np.median(sig[rr].sum(1))),
                               calls_mean=float(sig[rr].sum(1).mean()), truth_median=float(np.median(L.n[rr]))))
        if ln == 'H1':
            for m, v in out['per_kd'].items():
                pk[f'{key}::{m}'] = np.asarray(v, float)
        print(json.dumps({k: rec[-1][k] for k in ['line', 'key', 'part', 'overall', 'calls_median']}, default=float), flush=True)
        pd.DataFrame(rec).to_parquet(f_out)
        if pk:
            np.savez(f_pk, **pk)


# ---------------------------------------------------------------- 1c + 4
def source_n_real_all(targets):
    """(b) with every source table (none held out), as step15.source_n_real."""
    import h5py
    import stage1 as S
    from anndata.io import read_elem
    per = {t: [] for t in targets}; used = []
    for p in S.discover_tables([str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')]):
        with h5py.File(p, 'r') as f:
            o = read_elem(f['obs'])
        o = o[o['de_labels_known'].astype(bool) & np.isfinite(o['n_real'].astype(float))]
        if not len(o):
            continue
        used.append(p.stem); m = o.groupby(o['target_gene'].astype(str))['n_real'].mean()
        for t in targets:
            if t in m.index: per[t].append(float(m[t]))
    return np.array([np.mean(per[t]) if per[t] else np.nan for t in targets]), np.array([len(per[t]) for t in targets]), used


def qs(x):
    x = x[np.isfinite(x)]
    return {f'q{int(q * 100):02d}': float(np.quantile(x, q)) for q in (0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)} | \
        dict(n=int(len(x)), mean=float(x.mean()), min=float(x.min()), max=float(x.max()))


def predictor(a):
    TAB.mkdir(parents=True, exist_ok=True); res = dict(rho={}, calibration={})
    h1 = Line('H1'); hs = h1_halves(h1)
    (TAB / 'h1_truth.json').write_text(json.dumps(dict(targets=list(h1.targets), n=h1.n.tolist(), halves=[h.tolist() for h in hs])))
    res['rho']['H1_full'] = S11.spearman(h1.b, h1.n)
    for h in (0, 1):
        res['rho'][f'H1_half{h}'] = S11.spearman(h1.b[hs[h]], h1.n[hs[h]])
    res['n_finite'] = {'H1': int(np.isfinite(h1.b).sum())}
    for g in GUARD:
        L = Line(g); res['rho'][g] = S11.spearman(L.b, L.n); res['n_finite'][g] = int(np.isfinite(L.b).sum())
        res.setdefault('b_dist', {})[g] = qs(L.b)
        res.setdefault('truth', {})[g] = qs(L.n)
        del L
    res.setdefault('truth', {})['H1'] = qs(h1.n)
    # cross-fit calibration: predicted counts for half 1 - h from half h's quantiles; deciles of the predicted count
    cal = []
    for h in (0, 1):
        s = hs[1 - h]; pred = np.quantile(h1.n[hs[h]], h1.pct[s]); tr = h1.n[s]
        r = pd.Series(pred).rank(method='first').to_numpy(); dec = np.minimum((10 * (r - 1) / len(r)).astype(int), 9)
        for d in range(10):
            m = dec == d
            cal.append(dict(fit_half=h, scored_half=1 - h, decile=d, knockdowns=int(m.sum()),
                            pred_median=float(np.median(pred[m])), true_median=float(np.median(tr[m])),
                            true_q25=float(np.quantile(tr[m], 0.25)), true_q75=float(np.quantile(tr[m], 0.75))))
    cal = pd.DataFrame(cal); cal.to_csv(TAB / 'calibration_H1.csv', index=False)
    res['calibration'] = cal.groupby('decile')[['pred_median', 'true_median']].mean().to_dict('index')
    # 4: D/E/F
    tg = [t for t in (V / 'results/target_lists/final_DEF_300.txt').read_text().split() if t]
    bd, ntab, used = source_n_real_all(tg)
    pd.DataFrame(dict(target=tg, b=bd, tables_with_count=ntab)).to_csv(TAB / 'b_DEF.csv', index=False)
    res['DEF'] = dict(targets=len(tg), no_source_count=int((~np.isfinite(bd)).sum()),
                      no_source_targets=[t for t, v in zip(tg, bd) if not np.isfinite(v)],
                      tables_used=len(used), b=qs(bd), H1_b=qs(h1.b), tables_with_count=qs(ntab.astype(float)),
                      share_above_H1_max=float(np.mean(bd[np.isfinite(bd)] > np.nanmax(h1.b))),
                      share_below_H1_min=float(np.mean(bd[np.isfinite(bd)] < np.nanmin(h1.b))))
    for g in GUARD:
        res['DEF'][f'{g}_b'] = res['b_dist'][g]
    (TAB / 'step17_predictor.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(json.dumps(res, indent=1, default=float))


# ---------------------------------------------------------------- summary
def load_rows():
    t = pd.concat([pd.read_parquet(p) for p in sorted(TAB.glob('rows_*.parquet'))]).drop_duplicates(['line', 'key', 'part'])
    pk = {}
    for p in sorted(TAB.glob('perkd_*.npz')):
        pk.update(dict(np.load(p, allow_pickle=True)))
    return t, pk


def get(t, line, key, part):
    r = t[(t.line == line) & (t.key == key) & (t.part == part)]
    assert len(r) == 1, (line, key, part, len(r)); return r.iloc[0]


def summary(a):
    t, pk = load_rows()
    cols = ['overall'] + COMPS
    h1 = Line('H1'); hs = h1_halves(h1)

    def pair(line, key, part):
        r = get(t, line, key, part); f = get(t, line, f'flat|{int(r.calls_median)}', part)
        d = dict(line=line, key=key, part=part, calls_median=r.calls_median, flat_N=int(r.calls_median))
        for m in cols:
            d[m] = r[m]; d[f'flat_{m}'] = f[m]; d[f'd_{m}'] = r[m] - f[m]
        d['raw_mse'] = r['mse_raw']; return d

    res = {}
    # item 1: (b) per lam, fitted on half h, scored on 1 - h; guardrails (full)
    it1 = [pair('H1', f'b|{h}|{l}', f'half{1 - h}') | dict(fit_half=h, lam=l) for h in (0, 1) for l in LAM_B]
    it1g = [pair(g, f'b|{h}|{l}', 'full') | dict(fit_half=h, lam=l) for g in GUARD for h in (0, 1) for l in LAM_B]
    # lam fitted on the fit half (in-sample there), scored held out
    fitted = []
    for h in (0, 1):
        ins = [(l, get(t, 'H1', f'b|{h}|{l}', f'half{h}').overall) for l in LAM_B]
        lbest = max(ins, key=lambda x: x[1])[0]
        fitted.append(pair('H1', f'b|{h}|{lbest}', f'half{1 - h}') | dict(fit_half=h, lam=lbest,
                                                                         lam_grid_position=f'{LAM_B.index(lbest) + 1}/{len(LAM_B)}'))
    # lam fitted on the fit half by the gain over volume-matched flat K4 (the allocation question), scored held out
    fitted_gain = []
    for h in (0, 1):
        ins = [(l, pair('H1', f'b|{h}|{l}', f'half{h}')['d_overall']) for l in LAM_B]
        lbest = max(ins, key=lambda x: x[1])[0]
        fitted_gain.append(pair('H1', f'b|{h}|{lbest}', f'half{1 - h}') | dict(
            fit_half=h, lam=lbest, lam_grid_position=f'{LAM_B.index(lbest) + 1}/{len(LAM_B)}',
            insample_gain_by_lam={str(l): g for l, g in ins}))
    # 1b oracle per half, its flat match, and (b)'s share
    it1b = [pair('H1', f'oracle|{l}', f'half{h}') | dict(scored_half=h, lam=l) for h in (0, 1) for l in LAM_O]
    ob = pd.DataFrame(it1b)
    share = []
    for r in it1:
        s = 1 - r['fit_half']; o = ob[(ob.scored_half == s) & (ob.lam == r['lam'])]
        obest = ob[ob.scored_half == s].sort_values('d_overall').iloc[-1]
        og = float(o.d_overall.iloc[0]) if len(o) else np.nan
        share.append(dict(fit_half=r['fit_half'], scored_half=s, lam=r['lam'], b_gain_vs_flat=r['d_overall'],
                          oracle_gain_vs_flat_same_lam=og, share_same_lam=r['d_overall'] / og if og else np.nan,
                          oracle_best_lam=float(obest.lam), oracle_best_gain_vs_flat=float(obest.d_overall),
                          share_vs_oracle_best=r['d_overall'] / float(obest.d_overall)))
    # 2: volume-free
    it2 = [pair('H1', f'bfull|{l}', f'half{h}') | dict(scored_half=h, lam=l) for h in (0, 1) for l in LAM_F]
    it2g = [pair(g, f'bfull|{l}', 'full') | dict(lam=l) for g in GUARD for l in LAM_F]
    base = {ln: {m: get(t, ln, 'flat|1200', 'full')[m] for m in cols} for ln in ('H1',) + GUARD}
    base_h = {p: {m: get(t, 'H1', 'flat|1200', p)[m] for m in cols} for p in ('half0', 'half1')}
    # 3: per tercile of (b), held out: the fitted (b) rows, each knockdown from the direction it is scored in, vs flat
    terc = S15.terciles(np.where(np.isfinite(h1.b), h1.b, np.nanmedian(h1.b)))

    def kd_mix(keys):
        """per-knockdown dict assembled from key[h] on half 1 - h."""
        out = {}
        for m in [k.split('::', 1)[1] for k in pk if k.startswith(keys[0] + '::')]:
            v = np.full(len(h1.targets), np.nan)
            for h in (0, 1):
                v[hs[1 - h]] = pk[f'{keys[h]}::{m}'][hs[1 - h]]
            out[m] = v
        return out

    def item3(fit_rows, label):
        bk = [f'b|{h}|{fit_rows[h]["lam"]}' for h in (0, 1)]
        fk = [f'flat|{fit_rows[h]["flat_N"]}' for h in (0, 1)]
        kb, kf = kd_mix(bk), kd_mix(fk)
        rng = np.random.default_rng(2026); it3 = []
        for j, name in enumerate(S15.TERC):
            rr = np.flatnonzero(terc == j)
            sb, sf = h1.G.ref.summarise(kb, rows=rr), h1.G.ref.summarise(kf, rows=rr)
            boots = []
            for _ in range(1000):
                ii = rr[rng.integers(0, len(rr), len(rr))]
                xb, xf = h1.G.ref.summarise(kb, rows=ii), h1.G.ref.summarise(kf, rows=ii)
                boots.append({SP.SHORT[m]: xb['scaled'][m] - xf['scaled'][m] for m in SP.SIX} | {'overall': xb['overall'] - xf['overall']})
            B = pd.DataFrame(boots)
            for m in ['jaccard', 'fidelity', 'reach', 'overall', 'mse', 'pds', 'nmae']:
                mm = [k for k, v in SP.SHORT.items() if v == m]
                it3.append(dict(selection=label, tercile=name, knockdowns=len(rr), metric=m,
                                b=sb['overall'] if m == 'overall' else sb['scaled'][mm[0]],
                                flat=sf['overall'] if m == 'overall' else sf['scaled'][mm[0]],
                                b_raw=None if m == 'overall' else sb['raw'][mm[0]], flat_raw=None if m == 'overall' else sf['raw'][mm[0]],
                                diff_ci_lo=float(B[m].quantile(0.025)), diff_ci_hi=float(B[m].quantile(0.975)), true_de_median=float(np.median(h1.n[rr]))))
        for x in it3:
            x['diff'] = x['b'] - x['flat']
        return it3, bk, fk

    it3a, bk, fk = item3(fitted_gain, 'gain_fitted')
    it3b, bk2, fk2 = item3(fitted, 'overall_fitted')
    it3 = it3a + it3b
    res.update(item1=it1, item1_fitted=fitted, item1_fitted_gain=fitted_gain, item1_guard=it1g, item1b=it1b, item1b_share=share, item2=it2, item2_guard=it2g,
               item3=it3, flat1200=dict(full=base, halves=base_h), fitted_keys=dict(gain=dict(b=bk, flat=fk), overall=dict(b=bk2, flat=fk2)))
    for k in ('item1', 'item1_fitted', 'item1_fitted_gain', 'item1_guard', 'item1b', 'item1b_share', 'item2', 'item2_guard', 'item3'):
        pd.DataFrame(res[k]).to_csv(TAB / f'{k}.csv', index=False)
    (TAB / 'step17.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        for k in ('item1', 'item1_fitted', 'item1_fitted_gain', 'item1b_share', 'item2', 'item1_guard', 'item2_guard', 'item3'):
            print(k); print(pd.DataFrame(res[k]).round(4).to_string())


# ---------------------------------------------------------------- side check: flat K4 at the same MEAN volume
def meanmatch(a):
    """(b) concentrates calls on strong knockdowns, so at the same median its mean (total) volume is far above the flat
    row's. Side check (not in the task's grid): for every H1 allocation row, flat K4 at N = round(the row's mean calls on
    the part it is scored on). Guardrails: flat rows interpolated linearly in N between the existing flat rows."""
    t, pk = load_rows(); L = Line('H1'); n1, hs = h1_truth(); parts = dict(half0=hs[0], half1=hs[1])
    J, _ = jobs_for(L); JN = dict(J)
    want = {}
    for key in JN:
        if key.startswith(('b|', 'bfull|', 'oracle|')):
            sig = L.sig(JN[key])
            for p_, rr in parts.items():
                want[(key, p_)] = int(round(sig[rr].sum(1).mean()))
    f_out = TAB / 'rows_H1_meanmatch.parquet'
    rec = pd.read_parquet(f_out).to_dict('records') if f_out.exists() else []
    done = {r['key'] for r in rec}
    for N in sorted(set(want.values())):
        key = f'flat|{N}'
        if key in done or len(t[(t.line == 'H1') & (t.key == key)]):
            continue
        sig = L.sig(np.full(len(L.targets), float(N))); out = L.G.ref.proxy(mu=L.mu, sig=sig, per_kd=True)
        for p_, rr in parts.items():
            rec.append(S15.row('H1', L.G.ref.summarise(out['per_kd'], rows=rr), key=key, part=p_,
                               calls_median=float(np.median(sig[rr].sum(1))), calls_mean=float(sig[rr].sum(1).mean())))
        print(key, flush=True); pd.DataFrame(rec).to_parquet(f_out)
    t = pd.concat([t, pd.DataFrame(rec)]) if rec else t
    cols = ['overall'] + COMPS; out = []
    for (key, p_), N in want.items():
        if key.startswith('b|') and p_ != f'half{1 - int(key.split("|")[1])}':
            continue                                                       # (b) rows: held-out half only
        r = get(t, 'H1', key, p_); f = get(t, 'H1', f'flat|{N}', p_)
        out.append(dict(line='H1', key=key, part=p_, calls_mean=N, calls_median=r.calls_median,
                        **{m: r[m] for m in cols}, **{f'flatmean_{m}': f[m] for m in cols},
                        **{f'd_{m}': r[m] - f[m] for m in cols}))
    for g in GUARD:                                                        # interpolated
        fl = t[(t.line == g) & (t.part == 'full') & t.key.str.startswith('flat|')].copy()
        fl['N'] = fl.key.str.split('|').str[1].astype(float); fl = fl.sort_values('N')
        for r in t[(t.line == g) & (t.part == 'full') & t.key.str.startswith(('b|', 'bfull|'))].itertuples():
            fo = np.interp(r.calls_mean, fl.N, fl.overall)
            out.append(dict(line=g, key=r.key, part='full', calls_mean=r.calls_mean, calls_median=r.calls_median,
                            overall=r.overall, flatmean_overall=float(fo), d_overall=r.overall - float(fo), interpolated=True))
    o = pd.DataFrame(out); o.to_csv(TAB / 'meanmatch.csv', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        print(o[['line', 'key', 'part', 'calls_median', 'calls_mean', 'overall', 'flatmean_overall', 'd_overall']].round(4).to_string())


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['rows', 'predictor', 'summary', 'meanmatch'])
    p.add_argument('--part', type=int, default=0)
    p.add_argument('--nparts', type=int, default=1)
    p.add_argument('--line', default='H1')
    a = p.parse_args()
    {'rows': rows, 'predictor': predictor, 'summary': summary, 'meanmatch': meanmatch}[a.cmd](a)

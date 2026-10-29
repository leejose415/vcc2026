"""Scorer-space study part 2, Steps 10-12 (Jose 2026-10-28, results/stage1/2026-10-28_scorer_space/), on vcc2027.

Every arm is a noise-free requested prediction (means mu on the CPM scale + a call set) ranked by its six-component
proxy (scorer_proxy.Ref.proxy; "expected overall" = the mean of the six scaled proxies, validated in Step 9 to within
0.003 of the end-to-end overall). Per fold line the inputs come from scorer_study's fold cache (cache/fold_<line>
[_F1|_F2].npz: C+G and the Delta pool with the fold's lines removed, under source filter F0/F1/F2).

  step10 --line L     every (filter, arm, parameter) row for one fold line -> tables/step10_rows_<L>.parquet
  step10fit           leave-one-fold-out fits (arm e), per-filter summary, the carried (filter, arm) -> tables/step10.json

Fits: on the two training folds, the parameters maximising the mean expected overall over the training folds (each
fold weighted equally; the HCT116+HEK293T fold is the mean of its two lines), applied to the held-out fold.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
import scorer_study as SS                                                                  # noqa: E402
import scorer_arms as SA                                                                   # noqa: E402
import scorer_proxy as SP                                                                  # noqa: E402

OUT = SS.OUT
LINES = ['K562', 'HCT116', 'HEK293T', 'H1']
FILTERS = ['F0', 'F1', 'F2']
E_MC, E_M = [0.25, 0.5, 0.75, 1.0], [0.0, 0.1, 0.2]
E_M_EXT = [0.3, 0.5, 0.75, 1.0]    # side check past the Step 10 grid edge, F0 only


def suffix(filt):
    return '' if filt == 'F0' else f'_{filt}'


class Fold:
    """One fold line's requested inputs under one source filter."""

    def __init__(self, line, filt, ref):
        z = np.load(OUT / 'cache' / f'fold_{line}{suffix(filt)}.npz', allow_pickle=True)
        assert list(z['targets']) == ref.h.targets, (line, filt)
        self.ref, self.h = ref, ref.h
        self.unmeas = z['n_sources'] <= 0
        self.n_src = z['n_sources']
        self.lfc = np.where(self.unmeas, 0.0, z['lfc'].astype(np.float64))
        self.sig = z['de_sig'] & ~self.unmeas
        self.se = z['se_measurement'].astype(np.float64)
        self.m0 = z['m0'].astype(np.float64); self.gate = z['gate']; self.cprof = z['ctrl_profile'].astype(np.float64)
        self.dn = z['dpool_n']; self.dpool = np.where(self.dn > 0, z['dpool'].astype(np.float64), 0.0)
        v = self.m0 / 20.0; self.J2 = (np.log(2) * v / (1 + v)) ** 2
        ok = np.broadcast_to(self.gate[None, :], self.lfc.shape).copy()
        for i in self.h.rows:
            if self.h.own[i] >= 0: ok[i, self.h.own[i]] = False
        self.ok = ok
        self.true_delta = self.h.R - self.h.ctrl[None, :]
        self.ctrl_bulk = np.log1p(SS.BTS * self.m0 / self.m0.sum())

    def mu_lfc(self, lfc):
        return SA.mu_from_lfc(lfc, self.m0)

    def mu_delta(self, D):
        return SA.mu_from_profile(self.cprof[None, :] + D)

    def arm(self, name, m=1.0, mc=1.0):
        """(mu, calls) for an arm. a: m x LFC everywhere. c: called genes at C+G's LFC, non-called genes control profile
        + m x Delta_pool (control where no Delta source). e: called genes mc x LFC, non-called as c."""
        if name == 'CG':
            return self.mu_lfc(self.lfc), self.sig
        if name == 'a':
            return self.mu_lfc(m * self.lfc), self.sig
        if name in ('c', 'e'):
            mc_ = 1.0 if name == 'c' else mc
            return np.where(self.sig, self.mu_lfc(mc_ * self.lfc), self.mu_delta(m * self.dpool)), self.sig
        raise ValueError(name)

    def r_delta(self, mu):
        """J^2-weighted Pearson r between predicted and true scorer-space Delta over gated genes (own gene out)."""
        P = SA.bulk(mu)
        x = (P - self.ctrl_bulk[None, :])[self.ok]; y = self.true_delta[self.ok]
        w = np.broadcast_to(self.J2[None, :], P.shape)[self.ok]
        return SA.jr(x, y, w)

    def score(self, mu, sig):
        out = self.ref.proxy(mu=mu, sig=sig)
        row = dict(overall=out['overall'], mse_unclamped=out['mse_unclamped'], r_J2_delta=self.r_delta(mu),
                   calls_median=float(np.median(sig.sum(1))))
        for m in SP.SIX:
            row[f'{SP.SHORT[m]}_raw'] = out['raw'][m]; row[f'{SP.SHORT[m]}'] = out['scaled'][m]
        return row

    def coverage(self):
        g = self.ok
        return dict(share_kg_no_source=float((self.unmeas & g).sum() / g.sum()),
                    share_kg_no_delta_source=float(((self.dn <= 0) & g).sum() / g.sum()),
                    knockdowns_no_source=int((~(~self.unmeas & g).any(1)).sum()))


def refs_for(line):
    return SP.Ref(SS.EVAL_DIRS[line], SA.ratio_from_step4(), name=line)


def step10(a):
    line = a.line; ref = refs_for(line)
    f_out = OUT / 'tables' / f'step10_rows_{line}.parquet'
    rows = pd.read_parquet(f_out).to_dict('records') if f_out.exists() else []
    done = {(r['filter'], r['arm'], r['m'], r['mc']) for r in rows}
    for filt in FILTERS:
        fo = Fold(line, filt, ref); cov = fo.coverage()
        grid = [('CG', 1.0, 1.0), ('a', 0.5, 1.0), ('c', 0.1, 1.0)] + [('e', m, mc) for mc in E_MC for m in E_M]
        for arm, m, mc in grid:
            if (filt, arm, m, mc) in done:
                continue
            mu, sig = fo.arm(arm, m, mc)
            rows.append(dict(line=line, filter=filt, arm=arm, m=m, mc=mc, **fo.score(mu, sig), **cov))
            print(json.dumps({k: rows[-1][k] for k in ('line', 'filter', 'arm', 'm', 'mc', 'overall', 'mse', 'pds', 'nmae',
                                                       'fidelity', 'reach', 'jaccard', 'r_J2_delta')}, default=float), flush=True)
            pd.DataFrame(rows).to_parquet(f_out)


def fold_mean(t, folds):
    """Mean over folds of each fold's mean over its lines (t indexed by line)."""
    return np.mean([t[t.line.isin(SS.FOLD_LINES[f])].overall.mean() for f in folds])


def lofo(t, keys):
    """Per held-out fold: the keys maximising the training folds' mean expected overall."""
    out = {}
    for f_ in SS.FOLD_LINES:
        train = [g for g in SS.FOLD_LINES if g != f_]
        best = max(t.groupby(keys), key=lambda kv: fold_mean(kv[1], train))
        out[f_] = best[0] if isinstance(best[0], tuple) else (best[0],)
    return out


def step10fit(a):
    t = pd.concat([pd.read_parquet(OUT / 'tables' / f'step10_rows_{l}.parquet') for l in LINES])
    held, fits = [], {}
    for filt in FILTERS:
        tf = t[t['filter'] == filt]
        for arm in ('CG', 'a', 'c'):
            for f_, fl in SS.FOLD_LINES.items():
                for l in fl:
                    held.append(tf[(tf.arm == arm) & (tf.line == l)].iloc[0].to_dict())
        te = tf[tf.arm == 'e']; best = lofo(te, ['mc', 'm']); fits[filt] = {f_: dict(mc=b[0], m=b[1]) for f_, b in best.items()}
        for f_, fl in SS.FOLD_LINES.items():
            mc, m = best[f_]
            for l in fl:
                r = te[(te.line == l) & (te.mc == mc) & (te.m == m)].iloc[0].to_dict(); r['arm'] = 'e_lofo'; held.append(r)
    H = pd.DataFrame(held); H.to_csv(OUT / 'tables/step10_heldout.csv', index=False)
    comps = [SP.SHORT[m] for m in SP.SIX]
    ev = H[H.line == 'H1'].copy()
    ev['qualifies'] = (ev.mse > 0) & (ev[comps] > 0).all(1)
    q = ev[ev.qualifies]
    pick = (q if len(q) else ev).sort_values('overall', ascending=False).iloc[0]
    res = dict(fits_e=fits, carried=dict(filter=pick['filter'], arm=pick['arm'], m=pick['m'], mc=pick['mc'],
                                         evalA_overall=pick['overall'], qualified=bool(len(q) > 0)),
               evalA=ev[['filter', 'arm', 'm', 'mc', 'overall'] + comps + ['qualifies']].to_dict('records'))
    for filt in FILTERS:
        res[f'coverage_{filt}'] = json.loads((OUT / 'tables' / f'step10_coverage_{filt}.json').read_text())
    (OUT / 'tables/step10.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    pd.set_option('display.width', 250)
    print(H[['filter', 'arm', 'line', 'm', 'mc', 'overall'] + comps + ['mse_unclamped', 'r_J2_delta', 'share_kg_no_source']].to_string())
    print(json.dumps(res['carried'], default=float), json.dumps(fits, default=float))


def step10ext(a):
    """Arm (e) with non-called m past the Step 10 grid (E_M_EXT), F0 only. A side check: not part of the Step 10 carry rule."""
    line = a.line; ref = refs_for(line)
    f_out = OUT / 'tables' / f'step10ext_rows_{line}.parquet'
    rows = pd.read_parquet(f_out).to_dict('records') if f_out.exists() else []
    done = {(r['m'], r['mc']) for r in rows}
    fo = Fold(line, 'F0', ref); cov = fo.coverage()
    for mc in E_MC:
        for m in E_M_EXT:
            if (m, mc) in done:
                continue
            mu, sig = fo.arm('e', m, mc)
            rows.append(dict(line=line, filter='F0', arm='e', m=m, mc=mc, **fo.score(mu, sig), **cov))
            print(json.dumps({k: rows[-1][k] for k in ('line', 'm', 'mc', 'overall', 'mse', 'pds', 'nmae', 'fidelity', 'reach',
                                                       'jaccard', 'r_J2_delta')}, default=float), flush=True)
            pd.DataFrame(rows).to_parquet(f_out)


def step10extfit(a):
    """LOFO over the full arm (e) grid (Step 10 plus E_M_EXT), F0; held-out rows and Eval A qualification."""
    t = pd.concat([pd.read_parquet(OUT / 'tables' / f'step10_rows_{l}.parquet') for l in LINES])
    t = pd.concat([t[(t['filter'] == 'F0') & (t.arm == 'e')]] +
                  [pd.read_parquet(OUT / 'tables' / f'step10ext_rows_{l}.parquet') for l in LINES])
    best = lofo(t, ['mc', 'm']); held = []
    for f_, fl in SS.FOLD_LINES.items():
        mc, m = best[f_]
        for l in fl:
            r = t[(t.line == l) & (t.mc == mc) & (t.m == m)].iloc[0].to_dict(); r['arm'] = 'e_ext_lofo'; held.append(r)
    H = pd.DataFrame(held); comps = [SP.SHORT[m] for m in SP.SIX]
    ev = t[t.line == 'H1'].copy(); ev['qualifies'] = (ev.mse > 0) & (ev[comps] > 0).all(1)
    res = dict(fits={f_: dict(mc=b[0], m=b[1]) for f_, b in best.items()},
               heldout=H[['line', 'm', 'mc', 'overall'] + comps + ['r_J2_delta']].to_dict('records'),
               evalA_grid=ev[['m', 'mc', 'overall'] + comps + ['qualifies']].to_dict('records'))
    (OUT / 'tables/step10ext.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    pd.set_option('display.width', 250)
    print(H[['line', 'm', 'mc', 'overall'] + comps + ['r_J2_delta']].to_string())
    print(ev.pivot_table(index='mc', columns='m', values='overall').round(4).to_string())
    print(ev.pivot_table(index='mc', columns='m', values='nmae').round(4).to_string())
    print(json.dumps(res['fits'], default=float))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['step10', 'step10fit', 'step10ext', 'step10extfit'])
    p.add_argument('--line', choices=LINES)
    a = p.parse_args()
    {'step10': step10, 'step10fit': step10fit, 'step10ext': step10ext, 'step10extfit': step10extfit}[a.cmd](a)

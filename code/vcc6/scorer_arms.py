"""Scorer-space study, Steps 3, 5 and 6 (Jose 2026-10-28; results/stage1/2026-10-28_scorer_space/), on vcc2027, from
scorer_study.py's cache (cache/fold_<line>.npz: C+G per fold as run_line builds it, and Delta pooled with C+G's weights).

Predictions are noise-free pseudobulks in the scorer's space, P = log1p(5e4 x mu / sum mu), with
  mu from an LFC   (control_mean_cpm + 1e-9) x 2^lfc - 1e-9   (cell_generator.target_means: the generator's own inversion)
  mu from Delta    20 x (exp(ctrl_profile + Delta) - 1), floored at 0 (ctrl_profile = the target's context controls:
                   counts summed -> 50k -> log1p)
and every one is scored by analysis/mse_metric/mse_study.py:score_generated with C_pred = C0 x the Step 4 median ratio.
Genes no source measured are 0 (LFC) / control level (Delta), and are never called, as stage1.prepare_generator_inputs.

  ~/venv/bin/python code/vcc6/scorer_arms.py step3|step5|step6
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

OUT = SS.OUT
MS = [0, 0.1, 0.2, 0.3, 0.5, 0.75, 1]
AS = [0, 0.5, 1, 1.5]
LINES = ['K562', 'HCT116', 'HEK293T', 'H1']


def mu_from_lfc(lfc, m0):
    return np.maximum((m0[None, :] + 1e-9) * np.exp2(lfc) - 1e-9, 0.0)


def mu_from_profile(prof):
    return np.maximum(20.0 * np.expm1(prof), 0.0)


def bulk(mu):
    return np.log1p(SS.BTS * mu / mu.sum(1, keepdims=True))


class Line:
    def __init__(self, line, ratio, suffix=''):
        import mse_study as M
        self.M = M
        self.line = line
        self.h = SS.EvalPB(SS.EVAL_DIRS[line])
        z = np.load(OUT / 'cache' / f'fold_{line}{suffix}.npz', allow_pickle=True)
        assert list(z['targets']) == self.h.targets
        self.unmeas = z['n_sources'] <= 0
        self.lfc = np.where(self.unmeas, 0.0, z['lfc'].astype(np.float64))
        self.sig = z['de_sig'] & ~self.unmeas
        self.m0 = z['m0'].astype(np.float64); self.gate = z['gate']; self.cprof = z['ctrl_profile'].astype(np.float64)
        self.dpool = z['dpool'].astype(np.float64); self.dn = z['dpool_n']
        self.dpool = np.where(self.dn > 0, self.dpool, 0.0)
        self.sbar = self.dpool.mean(0)
        self.true_delta = self.h.R - self.h.ctrl[None, :]
        C0s, _ = SS.control_C0(SS.EVAL_DIRS[line])
        self.C0 = float(C0s.mean()); self.Cpred = np.full(self.h.T, self.C0 * ratio)
        self.strength = np.sqrt((np.where(self.gate[None, :], self.lfc, 0.0) ** 2).sum(1) / self.gate.sum())
        r = pd.Series(self.strength).rank(method='first', ascending=False).to_numpy(); n = len(r)
        self.terc = np.where(r <= n / 3, 'strong', np.where(r <= 2 * n / 3, 'mid', 'weak'))

    def score(self, P):
        s, num = self.M.score_generated(self.h, P, self.Cpred)
        return s, num

    # ------------------------------------------------------------------ arms
    def arm(self, name, m, a=None, ceiling=False):
        lfc, sig = self.lfc, self.sig
        if name == 'a':
            return bulk(mu_from_lfc(m * lfc, self.m0))
        if name == 'b':
            return bulk(mu_from_lfc(np.where(sig, lfc, m * lfc), self.m0))
        if name in ('c', 'd'):
            if name == 'c':
                D = m * self.dpool
            else:
                sb = self.true_delta.mean(0) if ceiling else self.sbar
                D = a * sb[None, :] + m * (self.dpool - sb[None, :])
                D = np.where(self.dn > 0, D, 0.0)
            mu_nc = mu_from_profile(self.cprof[None, :] + D)
            mu_c = mu_from_lfc(lfc, self.m0)
            return bulk(np.where(sig, mu_c, mu_nc))
        if name == 'control':
            return bulk(np.broadcast_to(self.m0, lfc.shape).copy())
        if name == 'sbar':
            return bulk(mu_from_profile(np.broadcast_to(self.cprof + self.sbar, lfc.shape).copy()))
        raise ValueError(name)


def ratio_from_step4():
    d = pd.read_parquet(OUT / 'tables/step4_per_knockdown.parquet')
    return float(np.median(d.ratio))


def load_lines():
    r = ratio_from_step4()
    return {l: Line(l, r) for l in LINES}, r


def truth_lfc(line, targets, genes):
    import stage1 as S
    import screen as SC
    by = {p.stem: p for p in S.discover_tables(SC.DIRS)}
    t = S.truth_signature(by[SS.LINE_TABLE[line]], targets, np.array(genes, dtype=object))
    return np.stack([t[k][0] if k in t else np.full(len(genes), np.nan) for k in targets])


# ---------------------------------------------------------------------------------------------------- step 3
def step3(a):
    L, ratio = load_lines()
    res, split_rows = {}, []
    for line, x in L.items():
        h = x.h
        rec = {}
        for nm, P in (('CG', x.arm('a', 1.0)), ('control', x.arm('control', 0)), ('sbar', x.arm('sbar', 0))):
            s, num = x.score(P)
            rec[nm] = dict(mse_raw=s['mse_raw'], mse_scaled=s['mse_scaled'], mse_unclamped=s['mse_unclamped'],
                           pds_scaled=s['pds_scaled'], budget_over_claim=s['budget_over_claim'])
        rec['anchors'] = dict(baseline_raw=h.bm, replicate_raw=h.rm, C0=x.C0, C_pred_used=float(x.Cpred[0]), ratio=ratio)
        # split of the plug-in term ||P - R||^2 (own gene out) for C+G
        P = x.arm('a', 1.0); d2 = (P - h.R) ** 2
        for i in h.rows:
            if h.own[i] >= 0: d2[i, h.own[i]] = 0.0
        tot = d2.sum()
        dec = np.full(h.G, -1); g = np.flatnonzero(x.gate)
        edges = np.quantile(x.m0[g], np.linspace(0.1, 0.9, 9)); dec[g] = np.searchsorted(edges, x.m0[g], side='right')
        parts = dict(called=float(d2[x.sig].sum() / tot), non_called=float(d2[~x.sig].sum() / tot),
                     ungated_genes=float(d2[:, dec < 0].sum() / tot),
                     by_ctrl_cpm_decile=[float(d2[:, dec == k].sum() / tot) for k in range(10)],
                     by_strength_tercile={t: float(d2[x.terc == t].sum() / tot) for t in ('strong', 'mid', 'weak')},
                     corrections_over_plugin=float(((x.Cpred.clip(max=h.Creal)) + h.Creal).sum() / tot))
        rec['plugin_split'] = parts
        # top 5% (k, g) by |predicted LFC|
        flat = np.abs(x.lfc).ravel(); cut = np.quantile(flat[flat > 0], 0.95)
        kk, gg = np.nonzero(np.abs(x.lfc) >= cut)
        rec['top5pct'] = dict(n=int(len(kk)), median_target_ctrl_cpm=float(np.median(x.m0[gg])),
                              share_in_gate=float(x.gate[gg].mean()))
        res[line] = rec
        print(line, json.dumps(rec, default=float), flush=True)
    (OUT / 'tables/step3.json').write_text(json.dumps(res, indent=1, default=float) + '\n')


def step3b(a):
    """(1) Top 5% (k, g) by |C+G LFC| per fold line: median control CPM in the TARGET's controls and in the SOURCES that
    hold k (each source table's own baseline control_mean_cpm, median over the sources voting for k).
    (2) From the harness DE tables of the generated C+G runs (HCT116 C+G 14:48; Eval A C+G 10-28): realised calls
    (p_adj < 0.05, own gene excluded) not in the requested call set, by target control-CPM decile, and their share."""
    import stage1 as S
    import screen as SC
    by = {p.stem: p for p in S.discover_tables(SC.DIRS)}
    genes = list(pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str))
    base_cpm, holds = {}, {}
    for name, p in by.items():
        b = pd.read_csv(p.with_name(name + '_baseline.csv'), index_col=0)['control_mean_cpm'].reindex(genes).to_numpy(float)
        base_cpm[name] = b
        z = OUT / 'delta' / f'{name}.npz'
    import h5py
    from anndata.io import read_elem
    for name, p in by.items():
        with h5py.File(p, 'r') as f:
            holds[name] = (set(read_elem(f['obs'])['target_gene'].astype(str)), S._line_key(str(read_elem(f['obs'])['background'].iloc[0])))
    res = {}
    for line in LINES:
        z = np.load(OUT / 'cache' / f'fold_{line}.npz', allow_pickle=True)
        tables = json.loads(str(z['tables'])); targets = list(z['targets'])
        lfc = np.where(z['n_sources'] > 0, z['lfc'], 0.0); flat = np.abs(lfc).ravel(); cut = np.quantile(flat[flat > 0], 0.95)
        kk, gg = np.nonzero(np.abs(lfc) >= cut)
        src_med = []
        for k, g in zip(kk, gg):
            vals = [base_cpm[t][g] for t in tables if targets[k] in holds[t][0] and np.isfinite(base_cpm[t][g])]
            src_med.append(np.median(vals) if vals else np.nan)
        res[line] = dict(top5_n=int(len(kk)), median_target_ctrl_cpm=float(np.median(z['m0'][gg])),
                         median_source_ctrl_cpm=float(np.nanmedian(src_med)))
    runs = {'HCT116': (V / 'results/stage1/e2e_final/CG/HCT116', 'evalHCT116final'),
            'H1': (V / 'results/stage1/2026-10-27_incumbent/evalA/CG', 'evalA')}
    for line, (d, ev) in runs.items():
        m = np.load(d / 'model' / 'predictions.npz', allow_pickle=True)
        req = {str(k): set(np.flatnonzero(m['de_sig'][i])) for i, k in enumerate(m['knockdowns'])}
        de = pd.read_parquet((d / 'salt_none' if (d / 'salt_none').exists() else d / 'e2e' / 'salt_none') / 'harness' / ev / 'de_pred.parquet')
        de = de[(de.p_adj < 0.05) & (de.feature != de.target)]
        gp = {g: j for j, g in enumerate(genes)}
        ctrl = pd.read_csv((d / 'controls' if (d / 'controls').exists() else d / 'e2e' / 'controls') / 'controls_baseline.csv', index_col=0)
        m0 = ctrl['control_mean_cpm'].reindex(genes).to_numpy(float); gate = ctrl['de_tested'].reindex(genes).fillna(False).to_numpy(bool)
        gi = np.flatnonzero(gate); edges = np.quantile(m0[gi], np.linspace(0.1, 0.9, 9))
        de['g'] = de.feature.map(gp); de = de[de.g.notna()]; de['g'] = de.g.astype(int)
        de['requested'] = [g in req.get(t, set()) for t, g in zip(de.target, de.g)]
        de['decile'] = np.where(gate[de.g], np.searchsorted(edges, m0[de.g], side='right'), -1)
        extra = de[~de.requested]
        per_kd = de.groupby('target').requested.agg(['size', 'sum'])
        res[f'{line}_realised_calls'] = dict(
            realised_total=int(len(de)), not_requested=int(len(extra)), share_not_requested=float(len(extra) / max(len(de), 1)),
            per_knockdown_not_requested_median=float((per_kd['size'] - per_kd['sum']).median()),
            per_knockdown_not_requested_mean=float((per_kd['size'] - per_kd['sum']).mean()),
            not_requested_by_ctrl_cpm_decile={int(k): int(v) for k, v in extra.decile.value_counts().sort_index().items()})
    print(json.dumps(res, indent=1, default=float))
    (OUT / 'tables/step3b.json').write_text(json.dumps(res, indent=1, default=float) + '\n')


# ---------------------------------------------------------------------------------------------------- step 5
def lofo(table, key):
    """m* (or (a*, m*)) per held-out line's fold: minimise the pooled numerator over the other two folds' lines."""
    out = {}
    for f_, fl in SS.FOLD_LINES.items():
        train = table[~table.line.isin(fl)]
        g = train.groupby(key)['num_sum'].sum()
        best = g.idxmin()
        out[f_] = best
    return out


def step5(a):
    L, ratio = load_lines()
    rows = []
    for line, x in L.items():
        for arm in ('a', 'b', 'c'):
            for m in MS:
                s, num = x.score(x.arm(arm, m))
                rows.append(dict(line=line, arm=arm, a=np.nan, m=m, num_sum=float(num.sum()), den_sum=float(x.h.den.sum()), **s))
        for ceiling in (False, True):
            for av in AS:
                for m in MS:
                    s, num = x.score(x.arm('d', m, av, ceiling))
                    rows.append(dict(line=line, arm='d_ceiling' if ceiling else 'd', a=av, m=m, num_sum=float(num.sum()),
                                     den_sum=float(x.h.den.sum()), **s))
        print(line, 'grid done', flush=True)
    t = pd.DataFrame(rows); t['budget_binds'] = t.budget_over_claim < 1
    t.to_parquet(OUT / 'tables/step5_grid.parquet')
    # leave-one-fold-out fits
    fits = {}
    for arm in ('a', 'b', 'c'):
        sub = t[t.arm == arm]; best = lofo(sub, 'm')
        fits[arm] = {f_: dict(m=float(m_), held_out=[dict(line=l, **sub[(sub.line == l) & (sub.m == m_)].iloc[0][
            ['mse_raw', 'mse_unclamped', 'pds_scaled']].to_dict()) for l in SS.FOLD_LINES[f_]]) for f_, m_ in best.items()}
    for arm in ('d', 'd_ceiling'):
        sub = t[t.arm == arm]; best = lofo(sub, ['a', 'm'])
        fits[arm] = {f_: dict(a=float(am[0]), m=float(am[1]), held_out=[dict(line=l, **sub[(sub.line == l) & (sub.a == am[0]) & (sub.m == am[1])].iloc[0][
            ['mse_raw', 'mse_unclamped', 'pds_scaled']].to_dict()) for l in SS.FOLD_LINES[f_]]) for f_, am in best.items()}
    # calibration: three betas and R0^2 per line, binned Delta calibration, Delta beta by tercile / source count
    cal = {}
    for line, x in L.items():
        h = x.h; G = h.G
        y = truth_lfc(line, h.targets, list(h.genes))
        ok = np.isfinite(y) & x.gate[None, :] & ~x.unmeas
        for i in h.rows:
            if h.own[i] >= 0: ok[i, h.own[i]] = False
        v = x.m0 / 20.0; J = np.log(2) * v / (1 + v)
        def fit(xv, yv, wv=None):
            wv = np.ones_like(xv) if wv is None else wv
            b = float((wv * xv * yv).sum() / (wv * xv * xv).sum())
            r0 = float(1 - (wv * (yv - b * xv) ** 2).sum() / (wv * yv * yv).sum())
            return dict(beta=b, R0_sq=r0, n=int(xv.size))
        X1, Y1 = x.lfc[ok], y[ok]; W = np.broadcast_to(J ** 2, y.shape)[ok]
        okd = (x.dn > 0) & x.gate[None, :] & ~x.sig
        for i in h.rows:
            if h.own[i] >= 0: okd[i, h.own[i]] = False
        Xd, Yd = x.dpool[okd], x.true_delta[okd]
        rec = dict(lfc_unweighted=fit(X1, Y1), lfc_J2=fit(X1, Y1, W), delta=fit(Xd, Yd))
        q = np.quantile(np.abs(Xd), np.linspace(0, 1, 11)); b_ = np.clip(np.searchsorted(q, np.abs(Xd), side='right') - 1, 0, 9)
        rec['delta_binned'] = [dict(decile=k, mean_pred=float(Xd[b_ == k].mean()), mean_true=float(Yd[b_ == k].mean()),
                                    beta=float((Xd[b_ == k] * Yd[b_ == k]).sum() / (Xd[b_ == k] ** 2).sum())) for k in range(10)]
        tk = np.broadcast_to(x.terc[:, None], okd.shape)[okd]
        rec['delta_beta_by_tercile'] = {t_: fit(Xd[tk == t_], Yd[tk == t_]) for t_ in ('strong', 'mid', 'weak')}
        nk = x.dn[okd]
        rec['delta_beta_by_sources'] = {lab: fit(Xd[(nk >= lo) & (nk <= hi)], Yd[(nk >= lo) & (nk <= hi)])
                                        for lab, lo, hi in (('1', 1, 1), ('2-3', 2, 3), ('4+', 4, 10 ** 6)) if ((nk >= lo) & (nk <= hi)).any()}
        ctrl_raw = float(t[(t.line == line) & (t.arm == 'a') & (t.m == 0)].mse_raw.iloc[0])
        rec['needed_R0_sq_to_beat_baseline'] = float(1 - h.bm / ctrl_raw)
        rec['control_raw'] = ctrl_raw
        cal[line] = rec
        print(line, json.dumps({k: rec[k] for k in ('lfc_unweighted', 'lfc_J2', 'delta', 'needed_R0_sq_to_beat_baseline')}), flush=True)
    (OUT / 'tables/step5.json').write_text(json.dumps(dict(ratio=ratio, fits=fits, calibration=cal), indent=1, default=float) + '\n')


def jr(x, y, w):
    """J^2-weighted Pearson r."""
    mx, my = np.average(x, weights=w), np.average(y, weights=w)
    c = np.average((x - mx) * (y - my), weights=w)
    return float(c / np.sqrt(np.average((x - mx) ** 2, weights=w) * np.average((y - my) ** 2, weights=w)))


def step7arc(a):
    """Step 7's Arc-first sub-arm: arm (c) at its m* (0.1 on every fold, Step 5) with Arc H1's weight x w_arc."""
    import h5py
    from anndata.io import read_elem
    import stage1 as S
    import screen as SC
    ratio = ratio_from_step4()
    m = json.loads((OUT / 'tables/step5.json').read_text())['fits']['c']
    by = {p.stem: p for p in S.discover_tables(SC.DIRS)}
    with h5py.File(by['arc_h1'], 'r') as f:
        arc = set(read_elem(f['obs'])['target_gene'].astype(str))
    defs = set(pd.read_csv(V / 'data_final/extracted/pert_counts.csv')['target_gene'].astype(str))
    res = dict(arc_h1_covers_DEF=len(arc & defs), rows=[])
    for f_ in ('K562', 'HCT116+HEK293T'):
        mstar = m[f_]['m']
        for line in SS.FOLD_LINES[f_]:
            y = None
            for w in (1, 2, 5, 20):
                x = Line(line, ratio, '' if w == 1 else f'_arcw{w}')
                if y is None:
                    y = truth_lfc(line, x.h.targets, list(x.h.genes))
                h = x.h; P = x.arm('c', mstar)
                v = x.m0 / 20.0; J2 = (np.log(2) * v / (1 + v)) ** 2
                for scope, rows in (('all', h.rows), ('arc_holds', np.array([i for i, t in enumerate(h.targets) if t in arc]))):
                    if len(rows) == 0: continue
                    sel = np.zeros(h.T, bool); sel[rows] = True
                    s_all, num = x.score(P)
                    u = float(num[sel].sum() / h.den[sel].sum())
                    pds = float(h.pds_rows(P)[sel].mean())
                    ok = np.isfinite(y) & x.gate[None, :] & ~x.unmeas & sel[:, None]
                    okd = (x.dn > 0) & x.gate[None, :] & ~x.sig & sel[:, None]
                    W = np.broadcast_to(J2, y.shape)
                    res['rows'].append(dict(fold=f_, line=line, w_arc=w, m=mstar, scope=scope, knockdowns=int(sel.sum()),
                                            mse_raw=u, mse_unclamped=(u - h.bm) / (h.rm - h.bm),
                                            pds_raw=pds, pds_scaled=h.scale(h.PDS, pds),
                                            r_J2_lfc=jr(x.lfc[ok], y[ok], W[ok]), r_J2_delta=jr(x.dpool[okd], x.true_delta[okd], W[okd])))
                    print(json.dumps(res['rows'][-1], default=float), flush=True)
    pd.DataFrame(res['rows']).to_parquet(OUT / 'tables/step7_arc.parquet')
    (OUT / 'tables/step7_arc.json').write_text(json.dumps(res, indent=1, default=float) + '\n')


# ---------------------------------------------------------------------------------------------------- step 6
def step6(a):
    L, ratio = load_lines()
    s5 = json.loads((OUT / 'tables/step5.json').read_text())
    t = pd.read_parquet(OUT / 'tables/step5_grid.parquet')
    # the arm + pooled m* with the lowest pooled unclamped MSE (over all lines, LOFO m* per fold)
    best_arm, best_val = None, np.inf
    for arm, fits in s5['fits'].items():
        if arm == 'd_ceiling': continue
        v = np.mean([h_['mse_unclamped'] for f_ in fits.values() for h_ in f_['held_out']])
        v = -v                                                     # unclamped scaled: higher is better
        if v < best_val: best_arm, best_val = arm, v
    res = dict(best_arm=best_arm, lines={})
    for f_, fl in SS.FOLD_LINES.items():
        fit = s5['fits'][best_arm][f_]
        for line in fl:
            x = L[line]; h = x.h
            P0 = x.arm(best_arm[0], fit['m'], fit.get('a'))
            nc = ~x.sig & x.gate[None, :]
            ctrl_b = np.log1p(SS.BTS * x.m0 / x.m0.sum())
            def med_delta(P):
                return np.array([np.median((P[i] - ctrl_b)[nc[i]]) for i in h.rows])
            med_pred, med_true = med_delta(P0), np.array([np.median((h.R[i] - h.ctrl)[nc[i]]) for i in h.rows])
            s0, num0 = x.score(P0)
            # offset share: n_nc x (median offset gap)^2 over the plug-in total
            d2 = ((P0 - h.R) ** 2); tot = d2.sum()
            share = float((nc.sum(1) * (med_pred - med_true) ** 2).sum() / tot)
            # kappa per knockdown on non-called genes, iterated through the renormalisation
            P = P0.copy()
            for _ in range(5):
                k = med_delta(P)
                mu = np.expm1(P); mu = np.where(nc, mu * np.exp(-k)[:, None], mu)
                P = bulk(mu)
            s1, _ = x.score(P)
            res['lines'][line] = dict(arm=best_arm, m=fit['m'], a=fit.get('a'), med_delta_pred=np.quantile(med_pred, [0.1, 0.5, 0.9]).tolist(),
                                      med_delta_true=np.quantile(med_true, [0.1, 0.5, 0.9]).tolist(), offset_share_of_plugin=share,
                                      before=dict(mse_raw=s0['mse_raw'], mse_unclamped=s0['mse_unclamped'], pds_scaled=s0['pds_scaled']),
                                      after_kappa=dict(mse_raw=s1['mse_raw'], mse_unclamped=s1['mse_unclamped'], pds_scaled=s1['pds_scaled']),
                                      residual_median_after=float(np.median(np.abs(med_delta(P)))))
            print(line, json.dumps(res['lines'][line], default=float), flush=True)
    (OUT / 'tables/step6.json').write_text(json.dumps(res, indent=1, default=float) + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['step3', 'step3b', 'step5', 'step6', 'step7arc'])
    a = p.parse_args()
    {'step3': step3, 'step3b': step3b, 'step5': step5, 'step6': step6, 'step7arc': step7arc}[a.cmd](a)

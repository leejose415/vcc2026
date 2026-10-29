"""Scorer-space study part 2, Step 11: adding signal on the carried (filter, arm) of Step 10 (Jose 2026-10-28).

The carried arm fixes the called genes (m_c x C+G's LFC; arm a scales non-called genes by the same m in LFC space).
Arms f-i replace the NON-CALLED genes with a Delta-space rule (control profile + Delta'; control where no Delta source):
  f  Delta' = m_t x Delta_pool, m_t per predicted-strength tercile, t in {strong, mid, weak}, grid {0.05, 0.1, 0.2}^3
  g  Delta' = m x (Delta_pool projected onto the fold sources' top-K Delta PCs on the line's gated genes; ungated genes
     unprojected), K in {5, 10, 20, 50}, m in {0.1, 0.2, 0.3, 0.5, 1}
  h  Delta' = m_h x Delta_pool with the n nearest complex partners' votes at weight w_nb, n in {3, 10},
     w_nb in {0.25, 0.5, 1}; m_h = the carried m (Delta arms), or {0.1, 0.2} when the carried arm has none
  i  Delta' = m x c x Delta_pool, c = z_a^2 / (z_a^2 + kappa) from the agreement z (>= 2 sources), c = c1 for a
     single source; m in {0.1, 0.2, 0.3, 0.5}, kappa in {1, 4, 16}, c1 in {0, 0.5, 1}
Fits as Step 10 (scorer_part2.lofo). Predicted strength = J^2-weighted ||Delta_pool,k||^2 over gated genes, own gene out;
true strength the same on the eval's R - control.

  step11 --line L   rows + per-line statistics -> tables/step11_rows_<L>.parquet, tables/step11_stats_<L>.json
  step11fit         fits, carried, top-30 H1 knockdowns with HGNC gene groups -> tables/step11.json
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
import scorer_proxy as SP                                                                  # noqa: E402
import scorer_part2 as P2                                                                  # noqa: E402

OUT = SS.OUT
COMPS = [SP.SHORT[m] for m in SP.SIX]
TERC = ('strong', 'mid', 'weak')


def carried10():
    return json.loads((OUT / 'tables/step10.json').read_text())['carried']


def spearman(x, y):
    from scipy.stats import spearmanr
    ok = np.isfinite(x) & np.isfinite(y)
    return float(spearmanr(x[ok], y[ok]).statistic) if ok.sum() > 2 else float('nan')


class Arms:
    def __init__(self, fo, s11, car):
        self.fo, self.s, self.car = fo, s11, car
        mc = 1.0 if car['arm'] in ('CG', 'c') else (car['m'] if car['arm'] == 'a' else car['mc'])
        self.mu_called = fo.mu_lfc(mc * fo.lfc)
        self.m_car = car['m'] if car['arm'] in ('c', 'e') else None
        D = fo.dpool; w = np.where(fo.ok, fo.J2[None, :], 0.0)
        self.pred_strength = (w * D ** 2).sum(1); self.true_strength = (w * fo.true_delta ** 2).sum(1)
        r = pd.Series(self.pred_strength).rank(method='first', ascending=False).to_numpy(); n = len(r)
        self.terc = np.where(r <= n / 3, 0, np.where(r <= 2 * n / 3, 1, 2))

    def mu(self, Dp):
        fo = self.fo
        return np.where(fo.sig, self.mu_called, fo.mu_delta(np.where(fo.dn > 0, Dp, 0.0)))

    def carried(self):
        fo, car = self.fo, self.car
        if car['arm'] in ('CG', 'a'):
            m = 1.0 if car['arm'] == 'CG' else car['m']
            return fo.mu_lfc(m * fo.lfc)
        return self.mu(self.m_car * fo.dpool)

    def project(self, K):
        fo, g = self.fo, self.s['gate']
        Vt = self.s['pc_Vt'][:K].astype(np.float64)
        D = fo.dpool.copy(); Dg = D[:, g]; D[:, g] = (Dg @ Vt.T) @ Vt
        return D

    def conf(self, kappa, c1):
        z2 = self.s['z_agree'].astype(np.float64) ** 2; n = self.s['n_bg']
        return np.where(n >= 2, z2 / (z2 + kappa), np.where(n == 1, c1, 0.0))


def stats(fo, s11, A, line):
    """S1-S4 statistics for one line."""
    out = dict(spearman_pred_vs_true_strength=spearman(A.pred_strength, A.true_strength))
    # S2: energy of the true Delta (gated, own out) captured by the source top-K PCs, own top-K (ceiling), random (floor)
    g = s11['gate']; T = np.where(fo.ok, fo.true_delta, 0.0)[:, g]; E = float((T ** 2).sum())
    _, sv, _ = np.linalg.svd(T, full_matrices=False)
    rng = np.random.default_rng(2026); en = {}
    for K in (5, 10, 20, 50):
        Vt = s11['pc_Vt'][:K].astype(np.float64)
        Q, _ = np.linalg.qr(rng.standard_normal((T.shape[1], K)))
        en[K] = dict(source_pcs=float(((T @ Vt.T) ** 2).sum() / E), own_top_k=float((sv[:K] ** 2).sum() / E),
                     random=float(((T @ Q) ** 2).sum() / E))
    out['pc_energy'] = en
    # S4: beta / R0^2 in Delta space on gated non-called genes, by agreement-z decile, C+G z decile and source tier
    ok = fo.ok & ~fo.sig & (fo.dn > 0)
    x, y = fo.dpool[ok], fo.true_delta[ok]
    za = s11['z_agree'][ok].astype(np.float64); zc = (np.abs(fo.lfc) / np.maximum(fo.se, 1e-12))[ok]; nb = s11['n_bg'][ok]

    def fit(m):
        if m.sum() < 50: return dict(beta=None, R0_sq=None, n=int(m.sum()))
        b = float((x[m] * y[m]).sum() / (x[m] ** 2).sum())
        return dict(beta=b, R0_sq=float(1 - ((y[m] - b * x[m]) ** 2).sum() / (y[m] ** 2).sum()), n=int(m.sum()))

    def by_decile(v, m0):
        q = np.nanquantile(v[m0], np.linspace(0, 1, 11)); d = np.clip(np.searchsorted(q, v, side='right') - 1, 0, 9)
        return [dict(decile=k, lo=float(q[k]), **fit(m0 & (d == k))) for k in range(10)]

    out['beta_by_agreement_z_decile'] = by_decile(za, np.isfinite(za))
    out['beta_by_CG_z_decile'] = by_decile(zc, np.isfinite(zc))
    out['beta_by_source_tier'] = {lab: fit((nb >= lo) & (nb <= hi)) for lab, lo, hi in (('1', 1, 1), ('2-3', 2, 3), ('4+', 4, 99))}
    out['sign_share_median'] = float(np.nanmedian(s11['sign_share'][ok]))
    return out


def step11(a):
    line, car = a.line, carried10(); filt = car['filter']; suf = '' if filt == 'F0' else f'_{filt}'
    ref = P2.refs_for(line); fo = P2.Fold(line, filt, ref)
    s11 = dict(np.load(OUT / 'cache' / f's11_{line}{suf}.npz', allow_pickle=True))
    A = Arms(fo, s11, car)
    st = stats(fo, s11, A, line)
    st.update(pool_check_maxdiff=float(s11['dpool_check_maxdiff']), true_strength=A.true_strength.tolist(),
              pred_strength=A.pred_strength.tolist(), targets=ref.h.targets)
    (OUT / 'tables' / f'step11_stats_{line}.json').write_text(json.dumps(st, default=float) + '\n')
    print(json.dumps({k: st[k] for k in ('spearman_pred_vs_true_strength', 'pc_energy', 'beta_by_source_tier', 'pool_check_maxdiff')},
                     default=float), flush=True)
    f_out = OUT / 'tables' / f'step11_rows_{line}.parquet'
    rows = pd.read_parquet(f_out).to_dict('records') if f_out.exists() else []
    done = {r['key'] for r in rows}

    def add(arm, key, mu, **p):
        if key in done: return
        rows.append(dict(line=line, arm=arm, key=key, **p, **fo.score(mu, fo.sig)))
        print(json.dumps({k: rows[-1][k] for k in ['arm', 'key', 'overall'] + COMPS}, default=float), flush=True)
        pd.DataFrame(rows).to_parquet(f_out)

    add('carried', 'carried', A.carried())
    grid = [0.05, 0.1, 0.2]
    for ms in [(x, y, z) for x in grid for y in grid for z in grid]:
        add('f', f'f{ms}', A.mu(np.array(ms)[A.terc][:, None] * fo.dpool), m_strong=ms[0], m_mid=ms[1], m_weak=ms[2])
    for K in (5, 10, 20, 50):
        D = A.project(K)
        for m in (0.1, 0.2, 0.3, 0.5, 1.0):
            add('g', f'g{K}_{m}', A.mu(m * D), K=K, m=m)
    for n in (3, 10):
        for w in (0.25, 0.5, 1.0):
            D = s11[f'dpool_nb{n}_w{w}'].astype(np.float64)
            for m in ([A.m_car] if A.m_car is not None else [0.1, 0.2]):
                add('h', f'h{n}_{w}_{m}', A.mu(m * D), n=n, w_nb=w, m=m)
    if A.m_car is None:
        for m in (0.1, 0.2):
            add('h', f'h0_0_{m}', A.mu(m * fo.dpool), n=0, w_nb=0.0, m=m)
    for m in (0.1, 0.2, 0.3, 0.5):
        for kappa in (1.0, 4.0, 16.0):
            for c1 in (0.0, 0.5, 1.0):
                c = A.conf(kappa, c1)
                add('i', f'i{m}_{kappa}_{c1}', A.mu(m * c * fo.dpool), m=m, kappa=kappa, c1=c1)


def build_mu(A, fo, s11, arm, p):
    """The means of a Step 11 arm with parameters p (as the grid in step11 builds them)."""
    if arm == 'carried':
        return A.carried()
    if arm == 'f':
        return A.mu(np.array([p['m_strong'], p['m_mid'], p['m_weak']])[A.terc][:, None] * fo.dpool)
    if arm == 'g':
        return A.mu(p['m'] * A.project(int(p['K'])))
    if arm == 'h':
        D = fo.dpool if int(p['n']) == 0 else s11[f"dpool_nb{int(p['n'])}_w{p['w_nb']}"].astype(np.float64)
        return A.mu(p['m'] * D)
    if arm == 'i':
        return A.mu(p['m'] * A.conf(p['kappa'], p['c1']) * fo.dpool)
    raise ValueError(arm)


def carried_mu(line):
    """(Fold, Arms, s11, mu) of Step 11's carried config for one line (its held-out fold's fitted parameters)."""
    car10 = carried10(); filt = car10['filter']; suf = '' if filt == 'F0' else f'_{filt}'
    s = json.loads((OUT / 'tables/step11.json').read_text()); arm = s['carried']['arm'].replace('_lofo', '')
    fold = next(f_ for f_, fl in SS.FOLD_LINES.items() if line in fl)
    p = s['fits'][arm][fold] if arm != 'carried' else {}
    ref = P2.refs_for(line); fo = P2.Fold(line, filt, ref)
    s11 = dict(np.load(OUT / 'cache' / f's11_{line}{suf}.npz', allow_pickle=True))
    A = Arms(fo, s11, car10)
    return fo, A, s11, build_mu(A, fo, s11, arm, p), dict(arm=arm, params=p, filter=filt)


def hgnc_groups():
    t = pd.read_csv(V / 'lib/downloads/hgnc_2026-10-28/hgnc_complete_set.txt', sep='\t', usecols=['symbol', 'gene_group'],
                    dtype=str)
    return {r.symbol: [g for g in str(r.gene_group).split('|') if g and g != 'nan'] for r in t.itertuples()}


def step11fit(a):
    t = pd.concat([pd.read_parquet(OUT / 'tables' / f'step11_rows_{l}.parquet') for l in P2.LINES])
    car10 = carried10(); held = []
    fam_keys = {'f': ['m_strong', 'm_mid', 'm_weak'], 'g': ['K', 'm'], 'h': ['n', 'w_nb', 'm'], 'i': ['m', 'kappa', 'c1']}
    fits = {}
    for l in P2.LINES:
        held.append(t[(t.arm == 'carried') & (t.line == l)].iloc[0].to_dict())
    for arm, keys in fam_keys.items():
        sub = t[t.arm == arm]
        if not len(sub): continue
        best = P2.lofo(sub, keys); fits[arm] = {f_: dict(zip(keys, map(float, b))) for f_, b in best.items()}
        for f_, b in best.items():
            for l in SS.FOLD_LINES[f_]:
                m = (sub.line == l)
                for k, v in zip(keys, b): m &= (sub[k] == v)
                r = sub[m].iloc[0].to_dict(); r['arm'] = f'{arm}_lofo'; held.append(r)
    H = pd.DataFrame(held); H.to_csv(OUT / 'tables/step11_heldout.csv', index=False)
    ev = H[H.line == 'H1'].copy(); ev['qualifies'] = (ev.mse > 0) & (ev[COMPS] > 0).all(1)
    q = ev[ev.qualifies]; pick = (q if len(q) else ev).sort_values('overall', ascending=False).iloc[0]
    st = {l: json.loads((OUT / 'tables' / f'step11_stats_{l}.json').read_text()) for l in P2.LINES}
    # strength between lines on shared knockdowns (true strengths)
    between = {}
    for i, l1 in enumerate(P2.LINES):
        for l2 in P2.LINES[i + 1:]:
            a1 = dict(zip(st[l1]['targets'], st[l1]['true_strength'])); a2 = dict(zip(st[l2]['targets'], st[l2]['true_strength']))
            sh = sorted(set(a1) & set(a2))
            between[f'{l1}~{l2}'] = dict(shared=len(sh), rho=spearman(np.array([a1[k] for k in sh]), np.array([a2[k] for k in sh])))
    # the 30 strongest H1 knockdowns and their HGNC gene groups
    grp = hgnc_groups(); h1 = pd.DataFrame(dict(target=st['H1']['targets'], strength=st['H1']['true_strength']))
    top = h1.sort_values('strength', ascending=False).head(30).copy()
    top['groups'] = top.target.map(lambda g: grp.get(g, []))
    cnt = pd.Series([g for gs in top.groups for g in gs]).value_counts()
    base = pd.Series([g for k in h1.target for g in grp.get(k, [])]).value_counts()
    shared = {g: dict(top30=int(n), panel=int(base.get(g, 0))) for g, n in cnt.items() if n >= 2}
    P = json.load(open(OUT / 'tables/step11_partners.json'))['replogle_k562_gw']
    top['k562_partners_cos_gt_0_5'] = top.target.map(lambda k: [p for p, c in P.get(k, []) if c > 0.5])
    top['partners_in_top30'] = top.k562_partners_cos_gt_0_5.map(lambda ps: [p for p in ps if p in set(top.target)])
    res = dict(carried_step10=car10, fits=fits, carried=dict(arm=pick['arm'], key=pick['key'], evalA_overall=pick['overall'],
                                                          qualified=bool(len(q) > 0)),
               evalA=ev[['arm', 'key', 'overall'] + COMPS + ['mse_unclamped', 'r_J2_delta', 'qualifies']].to_dict('records'),
               strength_between_lines=between,
               strength_pred_vs_true={l: st[l]['spearman_pred_vs_true_strength'] for l in P2.LINES},
               pc_energy={l: st[l]['pc_energy'] for l in P2.LINES},
               s4={l: {k: st[l][k] for k in ('beta_by_agreement_z_decile', 'beta_by_CG_z_decile', 'beta_by_source_tier', 'sign_share_median')}
                   for l in P2.LINES},
               partner_share=json.loads((OUT / 'tables/step11_partner_share.json').read_text()),
               h1_top30=top.to_dict('records'), h1_top30_shared_groups=shared)
    (OUT / 'tables/step11.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    pd.set_option('display.width', 250)
    print(H[['arm', 'line', 'key', 'overall'] + COMPS + ['mse_unclamped', 'r_J2_delta']].round(4).to_string())
    print(json.dumps(res['carried'], default=float)); print(json.dumps(fits, default=float))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['step11', 'step11fit'])
    p.add_argument('--line', choices=P2.LINES)
    a = p.parse_args()
    {'step11': step11, 'step11fit': step11fit}[a.cmd](a)

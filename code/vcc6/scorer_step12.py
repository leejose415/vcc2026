"""Scorer-space study part 2, Steps 12 and 13 (Jose 2026-10-28, results/stage1/2026-10-28_scorer_space/), on vcc2027.

Step 12, call rule on Step 11's carried magnitudes (the means stay; only the call set changes):
  K0  C+G's calls
  K1  BH at 0.05 on |m_call x LFC| / SE_target (stage1.predictive_calls on m_call x LFC; tested = the target's gated
      genes, own gene out), m_call in {1, 0.5, 0.25, 0.15}; K1 at 1 must reproduce K0 (checked)
  K2  K1, keeping calls with >= 2 (K2_2) or >= 4 (K2_4) source backgrounds whose Delta has the call's sign
  K3  K1, restricted to the target control-CPM deciles whose Eval A call precision (C+G's realised calls that are true
      DE, by decile) is at least the median over the deciles
  K4  per knockdown the top N_k genes by |LFC| / SE among tested genes, N_k = lambda x predicted strength (Step 11 S1),
      lambda set per line so the median N_k is in {25, 50, 100, 200, 400, 800}
  K2 and K3 use K1 at its fitted m_call (fit per fold, as Step 10); K4's median N is fitted the same way.
Step 13, generator check (measurement only): library sizes, per-gene zero fractions of generated vs real cells for the
same knockdown, and the realised calls that were not requested, on the existing Eval A C+G and HCT116 C+G cells.

  step12 --line L | step12fit | step13
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
import scorer_step11 as S11                                                                # noqa: E402

OUT = SS.OUT
COMPS = S11.COMPS
M_CALL = [1.0, 0.5, 0.25, 0.15]
N_MED = [25, 50, 100, 200, 400, 800]


def deciles(m0, gate):
    gi = np.flatnonzero(gate); edges = np.quantile(m0[gi], np.linspace(0.1, 0.9, 9))
    return np.where(gate, np.searchsorted(edges, m0, side='right'), -1)


def truth_calls(ref):
    d = ref.real_df.filter(ref.real_df['p_adj'] < 0.05)
    h = ref.h; gp = {g: j for j, g in enumerate(h.genes)}; tp = {t: i for i, t in enumerate(h.targets)}
    T = np.zeros((h.T, h.G), bool)
    T[[tp[t] for t in d['target'].to_list()], [gp[f] for f in d['feature'].to_list()]] = True
    for i in h.rows:
        if h.own[i] >= 0: T[i, h.own[i]] = False
    return T


def evalA_precision():
    """Eval A C+G: realised calls (harness de_pred, p_adj < 0.05, own gene out) that are true DE, by the target's
    context-control CPM decile."""
    f = OUT / 'tables/step12_evalA_precision.json'
    if f.exists():
        return json.loads(f.read_text())
    r = SP.RUNS9['EvalA_CG']; ref = SP.Ref(r['eval'], None, name='H1'); h = ref.h
    de = pd.read_parquet(r['harness'] / 'de_pred.parquet', columns=['target', 'feature', 'p_adj'])
    de = de[(de.p_adj < 0.05) & (de.feature != de.target)]
    gp = {g: j for j, g in enumerate(h.genes)}; tp = {t: i for i, t in enumerate(h.targets)}
    R = np.zeros((h.T, h.G), bool); R[de.target.map(tp).to_numpy(), de.feature.map(gp).to_numpy()] = True
    T = truth_calls(ref)
    cb = pd.read_csv(r['controls'], index_col=0)
    dec = deciles(cb['control_mean_cpm'].reindex(h.genes).to_numpy(float), cb['de_tested'].reindex(h.genes).fillna(False).to_numpy(bool))
    prec = {int(d): float((R & T)[:, dec == d].sum() / max(R[:, dec == d].sum(), 1)) for d in range(10)}
    med = float(np.median(list(prec.values())))
    res = dict(precision_by_decile=prec, median=med, kept_deciles=[d for d, p in prec.items() if p >= med],
               realised=int(R.sum()), true_positive=int((R & T).sum()))
    f.write_text(json.dumps(res, indent=1) + '\n')
    return res


def step12(a):
    import stage1 as S
    line = a.line
    fo, A, s11, mu, info = S11.carried_mu(line)
    h = fo.h; T = truth_calls(fo.ref)
    tested = fo.gate.copy()
    def k1(m):
        out = np.zeros_like(fo.sig)
        for i in h.rows:
            t = tested.copy()
            if h.own[i] >= 0: t[h.own[i]] = False
            out[i] = S.predictive_calls(m * fo.lfc[i], fo.se[i], t)[0]
        return out & ~fo.unmeas
    sgn = np.sign(fo.lfc); dsg = np.sign(s11['dpool_check']); share = np.nan_to_num(s11['sign_share']); nb = s11['n_bg']
    n_same = np.round(np.where(dsg == sgn, share * nb, (1 - share) * nb))
    prec = evalA_precision(); keep_dec = np.isin(deciles(fo.m0, fo.gate), prec['kept_deciles'])
    z = np.abs(fo.lfc) / np.maximum(fo.se, S.SE_MIN)
    rows = []

    def add(var, sig, **p):
        r = dict(line=line, variant=var, **p, **fo.score(mu, sig),
                 calls_mean=float(sig.sum(1).mean()), truth_calls_median=float(np.median(T.sum(1))),
                 truth_calls_mean=float(T.sum(1).mean()))
        rows.append(r)
        print(json.dumps({k: r[k] for k in ['variant', 'calls_median', 'calls_mean', 'overall'] + COMPS} | p, default=float), flush=True)

    add('K0', fo.sig, m_call=np.nan, support=np.nan, n_med=np.nan)
    for m in M_CALL:
        s1 = k1(m)
        if m == 1.0:
            chk = int((s1 != fo.sig).sum())
        add('K1', s1, m_call=m, support=np.nan, n_med=np.nan)
        for sup in (2, 4):
            add(f'K2_{sup}', s1 & (n_same >= sup), m_call=m, support=sup, n_med=np.nan)
        add('K3', s1 & keep_dec[None, :], m_call=m, support=np.nan, n_med=np.nan)
    ok = tested[None, :] & ~fo.unmeas & (fo.lfc != 0)
    for i in h.rows:
        if h.own[i] >= 0: ok[i, h.own[i]] = False
    zz = np.where(ok, z, -1.0); order = np.argsort(-zz, axis=1)
    for nm in N_MED:
        lam = nm / np.median(A.pred_strength)
        N = np.minimum(np.round(lam * A.pred_strength).astype(int), ok.sum(1))
        sig = np.zeros_like(fo.sig)
        for i in h.rows:
            sig[i, order[i, :N[i]]] = True
        add('K4', sig, m_call=np.nan, support=np.nan, n_med=nm)
    t = pd.DataFrame(rows); t['k1_at_1_mismatch_vs_K0'] = chk; t['magnitude'] = json.dumps(info, default=float)
    t.to_parquet(OUT / 'tables' / f'step12_rows_{line}.parquet')


def step12fit(a):
    t = pd.concat([pd.read_parquet(OUT / 'tables' / f'step12_rows_{l}.parquet') for l in P2.LINES])
    held = []
    for l in P2.LINES:
        held.append(t[(t.variant == 'K0') & (t.line == l)].iloc[0].to_dict())
    k1 = P2.lofo(t[t.variant == 'K1'], ['m_call']); k4 = P2.lofo(t[t.variant == 'K4'], ['n_med'])
    for f_, fl in SS.FOLD_LINES.items():
        m = k1[f_][0]
        for l in fl:
            for v in ('K1', 'K2_2', 'K2_4', 'K3'):
                r = t[(t.line == l) & (t.variant == v) & (t.m_call == m)].iloc[0].to_dict(); r['variant'] = f'{v}_lofo'; held.append(r)
            r = t[(t.line == l) & (t.variant == 'K4') & (t.n_med == k4[f_][0])].iloc[0].to_dict(); r['variant'] = 'K4_lofo'; held.append(r)
    H = pd.DataFrame(held); H.to_csv(OUT / 'tables/step12_heldout.csv', index=False)
    ev = H[H.line == 'H1'].copy(); ev['qualifies'] = (ev.mse > 0) & (ev[COMPS] > 0).all(1)
    q = ev[ev.qualifies]; pick = (q if len(q) else ev).sort_values('overall', ascending=False).iloc[0]
    res = dict(precision=evalA_precision(), fits=dict(K1_m_call={f: float(v[0]) for f, v in k1.items()},
                                                      K4_n_med={f: float(v[0]) for f, v in k4.items()}),
               carried=dict(variant=pick['variant'], m_call=pick['m_call'], support=pick['support'], n_med=pick['n_med'],
                            evalA_overall=pick['overall'], qualified=bool(len(q) > 0)),
               k1_at_1_mismatch={l: int(t[t.line == l].k1_at_1_mismatch_vs_K0.iloc[0]) for l in P2.LINES},
               evalA=ev[['variant', 'm_call', 'support', 'n_med', 'calls_median', 'calls_mean', 'truth_calls_median', 'overall'] + COMPS + ['qualifies']].to_dict('records'))
    (OUT / 'tables/step12.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    pd.set_option('display.width', 250)
    print(H[['variant', 'line', 'm_call', 'support', 'n_med', 'calls_median', 'calls_mean', 'truth_calls_median', 'overall'] + COMPS].round(4).to_string())
    print(json.dumps(res['carried'], default=float), json.dumps(res['fits'], default=float), res['k1_at_1_mismatch'])


# ---------------------------------------------------------------------------------------------------- step 13
def group_stats(X, labels, targets):
    """Per target: cell count, per-gene zero fraction, library sizes."""
    X = X.tocsr()
    lib = np.asarray(X.sum(1)).ravel()
    zf = np.zeros((len(targets), X.shape[1]), np.float32); libs = {}
    for i, t in enumerate(targets):
        idx = np.flatnonzero(labels == t); sub = X[idx]
        nnz = np.bincount(sub.indices, minlength=X.shape[1])
        zf[i] = 1.0 - nnz / max(len(idx), 1); libs[t] = lib[idx]
    return zf, libs, lib


def step13(a):
    import anndata as ad
    res = {}
    for run in ('EvalA_CG', 'HCT116_CG'):
        r = SP.RUNS9[run]; ref = SP.Ref(r['eval'], None, name=r['line']); h = ref.h
        kds, lfc, sig = SP.requested(r['inp'])
        cb = pd.read_csv(r['controls'], index_col=0)
        m0 = cb['control_mean_cpm'].reindex(h.genes).to_numpy(float); gate = cb['de_tested'].reindex(h.genes).fillna(False).to_numpy(bool)
        dec = deciles(m0, gate)
        pred_path = r['harness'].parent.parent / 'prediction.h5ad' if (r['harness'].parent.parent / 'prediction.h5ad').exists() else r['harness'].parent / 'prediction.h5ad'
        g = ad.read_h5ad(pred_path); gl = g.obs['target_gene' if 'target_gene' in g.obs else 'target'].astype(str).to_numpy()
        assert list(g.var_names.astype(str)) == list(h.genes)
        zg, libg, _ = group_stats(g.X, gl, h.targets)
        gctrl = np.asarray(g.X[gl == SS.CTRL].sum(1)).ravel() if (gl == SS.CTRL).any() else np.array([])
        del g
        rr = ad.read_h5ad(r['eval'] / 'real.h5ad'); rl = rr.obs['target_gene' if 'target_gene' in rr.obs else 'target'].astype(str).to_numpy()
        zr, libr, _ = group_stats(rr.X, rl, h.targets)
        rctrl = np.asarray(rr.X[rl == SS.CTRL].sum(1)).ravel()
        del rr
        q = [0.05, 0.25, 0.5, 0.75, 0.95]
        gk = np.concatenate(list(libg.values())); rk = np.concatenate(list(libr.values()))
        lib = dict(generated_knockdown_cells=np.quantile(gk, q).tolist(), generated_control_cells=np.quantile(gctrl, q).tolist() if len(gctrl) else None,
                   real_scoring_controls=np.quantile(rctrl, q).tolist(), real_knockdown_cells=np.quantile(rk, q).tolist(),
                   ratio_median_generated_kd_over_real_kd=float(np.median(gk) / np.median(rk)),
                   ratio_median_generated_kd_over_real_controls=float(np.median(gk) / np.median(rctrl)))
        dz = zg - zr
        by = [dict(decile=d, mean_zero_frac_gen=float(zg[:, dec == d].mean()), mean_zero_frac_real=float(zr[:, dec == d].mean()),
                   mean_diff=float(dz[:, dec == d].mean()), share_absdiff_gt_0_05=float((np.abs(dz[:, dec == d]) > 0.05).mean()))
              for d in range(-1, 10)]
        # realised calls not requested
        de = pd.read_parquet(r['harness'] / 'de_pred.parquet', columns=['target', 'feature', 'p_adj'])
        de = de[(de.p_adj < 0.05) & (de.feature != de.target)]
        gp = {x: j for j, x in enumerate(h.genes)}; tp = {t: i for i, t in enumerate(h.targets)}
        R = np.zeros((h.T, h.G), bool); R[de.target.map(tp).to_numpy(), de.feature.map(gp).to_numpy()] = True
        req = sig.copy()
        for i in h.rows:
            if h.own[i] >= 0: req[i, h.own[i]] = False
        T = truth_calls(ref)
        import scorer_arms as SA
        P = SA.bulk(SA.mu_from_lfc(lfc, m0)); ctrl_b = np.log1p(SS.BTS * m0 / m0.sum()); Dreq = np.abs(P - ctrl_b[None, :])
        leak = R & ~req; false_ = R & ~T; miss = req & ~R
        def describe(M):
            kk, gg = np.nonzero(M)
            return dict(n=int(M.sum()), per_kd_median=float(np.median(M.sum(1))),
                        by_decile={int(d): int((dec[gg] == d).sum()) for d in range(-1, 10)},
                        abs_requested_delta_q=np.quantile(Dreq[kk, gg], [0.1, 0.5, 0.9]).tolist() if len(kk) else None,
                        zero_frac_absdiff_gt_0_05=int((np.abs(dz[kk, gg]) > 0.05).sum()),
                        zero_frac_diff_median=float(np.median(dz[kk, gg])) if len(kk) else None)
        res[run] = dict(library_size=lib, zero_fraction_by_decile=by, leaked_not_requested=describe(leak),
                        false_vs_truth=describe(false_), leaked_and_false=describe(leak & false_),
                        requested_not_realised=describe(miss),
                        realised=int(R.sum()), requested=int(req.sum()), truth=int(T.sum()))
        print(run, json.dumps({k: res[run][k] for k in ('library_size', 'leaked_not_requested', 'leaked_and_false')}, default=float), flush=True)
    (OUT / 'tables/step13.json').write_text(json.dumps(res, indent=1, default=float) + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['step12', 'step12fit', 'step13'])
    p.add_argument('--line', choices=P2.LINES)
    a = p.parse_args()
    {'step12': step12, 'step12fit': step12fit, 'step13': step13}[a.cmd](a)

"""Step 19 (Jose 2026-11-01): studies before the next D/E/F file. vcc2026, CPU only; no cells generated, no harness
scoring of new cells, no D/E/F file, no upload. Write-up: section 11 of results/stage1/2026-10-28_scorer_space/report.md.
Outputs: results/stage1/2026-11-01_step19/{tables,logs,figures}. Frozen code (generator, harness, run_line, boot_score,
config 1' files, vcc_2026_7, Step 18 outputs) is imported, never edited.

  part0   checks: (1) combining the two 15c runs reproduces each run at the endpoints (all six components, 1e-6);
          (2) config 1' + K0 proxy on the generator's input (15c sec. 4, 0.2456) and on the fold cache (15a, 0.2470,
          halves 0.2679 / 0.2261)
  partA   hybrid calls from the existing 15c cells: knockdown k takes K4 1,200's outputs if b_k >= T, else K0's
  partB   power and depth: true DE count vs real cell count (partial on log1p(b)); control-cell depth and offset
  partC   fingerprint D/E/F controls: marker CPM / detection; Spearman to source and A/B/C control profiles
  partD   called-gene scale by expression tier (proxy): D1 diagnostics on the uploaded file's input; D2 375-row grid
  partE   program ceiling: PCA bases from source knockdowns, captured energy of H1's truth
  test    unit tests of the summary and fit functions on small slices

  ~/venv/bin/python code/vcc6/step19.py {test,part0,partA,partB,partC,partD,partE} [--line L]
"""
import argparse
import gc
import json
import os
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
for p in (V / 'code/vcc6', V / 'scripts', V / 'scripts/harness', V / 'code/cellgen', V / 'analysis/mse_metric'):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

OUT = V / 'results/stage1/2026-11-01_step19'
TAB, LOG, FIG = OUT / 'tables', OUT / 'logs', OUT / 'figures'
S15 = V / 'results/stage1/2026-10-29_step15'
C15 = S15 / 'c'
S17 = V / 'results/stage1/2026-10-30_step17'
S18 = V / 'results/stage1/2026-10-30_step18'
FE = V / 'results/stage1/final_evals/harness/lib'
EVAL = {'EvalA': V / 'lib/harness/evalA', 'HCT116final': FE / 'evalHCT116final',
        'K562final': FE / 'evalK562final', 'HEK293Tfinal': FE / 'evalHEK293Tfinal'}
PROXY_NAME = {'EvalA': 'H1', 'HCT116final': 'HCT116', 'K562final': 'K562', 'HEK293Tfinal': 'HEK293T'}
CFG = {'K0': 'config1p_K0', 'K4': 'config1p_K4_1200'}
LINES_A = ('EvalA', 'HCT116final')
T_GRID = [1, 2, 5, 10, 20, 35, 50, 100, 200, 500]
COMPS = ['mse', 'pds', 'nmae', 'fidelity', 'reach', 'jaccard']
TERC = ('weak', 'mid', 'strong')
N_BOOT, SEED = 2000, 2026


def log(msg):
    print(f'[{time.strftime("%H:%M:%S", time.gmtime())}] {msg}', flush=True)


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def avail_gb():
    for ln in open('/proc/meminfo'):
        if ln.startswith('MemAvailable'):
            return int(ln.split()[1]) / 1e6
    return float('nan')


def mem_before(label, est_gb):
    """Log the estimated peak of a load before it starts, with what is available."""
    log(f'MEM before {label}: estimated peak +{est_gb:.1f} GB; process peak so far {rss_gb():.1f} GB; available {avail_gb():.1f} GB')
    if est_gb > avail_gb() - 2:
        raise MemoryError(f'{label}: estimated {est_gb:.1f} GB exceeds available {avail_gb():.1f} GB minus 2 GB headroom')


def mem_after(label):
    log(f'MEM after {label}: process peak {rss_gb():.1f} GB; available {avail_gb():.1f} GB')


def save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=1, default=float) + '\n')


def grid_pos(values, v):
    """Position of v in an ordered grid: 'k/n (edge|interior)'."""
    values = list(values); i = values.index(v)
    return f'{i + 1}/{len(values)} ({"edge" if i in (0, len(values) - 1) else "interior"})'


# ================================================================ shared: refs, runs, scaling
_REFS = {}


def ref_for(line):
    """scorer_proxy.Ref (all six components) for an eval, cached per process."""
    import scorer_proxy as SP
    import scorer_arms as SA
    if line not in _REFS:
        mem_before(f'Ref {line}', 4.0)
        _REFS[line] = SP.Ref(EVAL[line], SA.ratio_from_step4(), name=PROXY_NAME[line])
        mem_after(f'Ref {line}')
    return _REFS[line]


def hdir(key, line):
    import step15c
    return step15c.harness_dir(CFG[key], line)


def results_wide(key, line):
    return pd.read_csv(hdir(key, line) / 'results.csv').pivot(index='perturbation', columns='metric', values='value')


class Fast:
    """boot_score's aggregate + scaling, vectorised over resamples. Per-knockdown arrays in h.targets order:
    MSE num / den (ratio of sums), PDS and the four DE components (nanmean). Scaling = EvalPB.scale (score_one on the
    bundle's replicate entries), the same map boot_score.Run.scaled reaches through score_metrics (unit-tested)."""
    import scorer_proxy as _SP
    KEYS = {'mse': _SP.MSE, 'pds': _SP.PDS, 'nmae': _SP.DE[0], 'fidelity': _SP.DE[1], 'reach': _SP.DE[2], 'jaccard': _SP.DE[3]}

    def __init__(self, h, num, den, comp):
        self.h, self.num, self.den = h, np.asarray(num, float), np.asarray(den, float)
        self.comp = {k: np.asarray(v, float) for k, v in comp.items()}            # pds, nmae, fidelity, reach, jaccard

    @classmethod
    def from_wide(cls, h, W):
        W = W.reindex(h.targets)
        return cls(h, W['expr_mse_unbiased_capped'], W['expr_distance_unbiased'],
                   {k: W[m] for k, m in cls.KEYS.items() if k != 'mse'})

    def raw(self, I):
        """I: (draws, n) integer positions -> dict of raw aggregates per draw."""
        I = np.atleast_2d(I)
        out = {'mse': self.num[I].sum(1) / self.den[I].sum(1)}
        for k, v in self.comp.items():
            x = v[I]; ok = np.isfinite(x)
            out[k] = np.where(ok, x, 0.0).sum(1) / ok.sum(1)
        return out

    def scaled(self, I):
        r = self.raw(I)
        s = {k: np.array([self.h.scale(self.KEYS[k], u) for u in r[k]]) for k in COMPS}
        s['overall'] = np.mean([s[k] for k in COMPS], 0)
        return s, r

    def point(self, rows=None):
        rows = np.arange(len(self.num)) if rows is None else np.asarray(rows)
        s, r = self.scaled(rows[None, :])
        return {k: float(v[0]) for k, v in s.items()}, {k: float(v[0]) for k, v in r.items()}


def boot_index(n, n_boot=N_BOOT, seed=SEED):
    """boot_score's resample sequence: one rng.integers(0, n, n) per draw, seed 2026."""
    rng = np.random.default_rng(seed)
    return np.stack([rng.integers(0, n, n) for _ in range(n_boot)])


# ================================================================ Part 0 / A: the 15c cells
def predside(key, line):
    """Predicted pseudobulk means P and jackknife C from the 15c cells (cell-eval2's own prep), cached."""
    import scorer_study as SS
    f = TAB / f'predside_{key}_{line}.npz'
    if f.exists():
        z = np.load(f, allow_pickle=True); return z['P'], z['C']
    import anndata as ad
    from cell_eval2.moments import correction_for
    from cell_eval2.prep import pseudobulk_bulk_lognorm_with_moments
    h = ref_for(line).h
    pred = hdir(key, line).parent.parent / 'prediction.h5ad'
    # The whole file (Eval A: ~850 M nnz incl. controls) would need ~20 GB; group means and their jackknife are
    # per-group, so the targets are read in chunks of CHUNK knockdowns (Part 0 check 1 validates against the harness).
    A = ad.read_h5ad(pred, backed='r')
    col = 'target_gene' if 'target_gene' in A.obs else 'target'
    lab = A.obs[col].astype(str).to_numpy()
    nnz_per_cell = A.file['X/data'].shape[0] / A.n_obs
    CHUNK = 40
    P = np.full((h.T, h.G), np.nan); C = np.full(h.T, np.nan)
    pos = {t: i for i, t in enumerate(h.targets)}
    for c0 in range(0, h.T, CHUNK):
        ts = h.targets[c0:c0 + CHUNK]
        rows = np.flatnonzero(np.isin(lab, ts))
        if c0 == 0:
            mem_before(f'cells {key} {line}, chunk of {len(ts)} knockdowns ({len(rows)} cells)', len(rows) * nnz_per_cell * 12 * 3 / 1e9 + 1)
        sub = A[rows].to_memory()
        perts, means, mom = pseudobulk_bulk_lognorm_with_moments(sub, col, bulk_target_sum=SS.BTS)
        Cc = np.asarray(correction_for(mom, perts, means, comparator='bulk_lognorm'), np.float64)
        for j, p_ in enumerate(map(str, perts)):
            if p_ in pos:
                P[pos[p_]] = np.asarray(means[j], np.float64); C[pos[p_]] = Cc[j]
        del sub, means, mom; gc.collect()
    A.file.close()
    assert np.isfinite(C).all() and np.isfinite(P).all(), f'{key} {line}: knockdowns missing from the prediction'
    mem_after(f'cells {key} {line}')
    np.savez_compressed(f, P=P, C=C, targets=np.array(h.targets))
    gc.collect()
    return P, C


def combined(line, onK4, PC=None):
    """Per-knockdown values of the combined panel: DE components from the run that supplies each knockdown; MSE and
    PDS recomputed through mse_study.score_generated / EvalPB.pds_rows on the combined P (each knockdown's C_pred from
    its run; the #348 budget formed on the combined panel)."""
    import mse_study as M
    h = ref_for(line).h
    onK4 = np.asarray(onK4, bool)
    W0, W4 = results_wide('K0', line).reindex(h.targets), results_wide('K4', line).reindex(h.targets)
    W = W0.copy(); W.loc[onK4] = W4.loc[onK4].to_numpy()
    (P0, C0), (P4, C4) = PC if PC is not None else (predside('K0', line), predside('K4', line))
    P = np.where(onK4[:, None], P4, P0); C = np.where(onK4, C4, C0)
    s, num = M.score_generated(h, P, C)
    W['expr_mse_unbiased_capped'] = num; W['expr_distance_unbiased'] = h.den; W['pds_cosine'] = h.pds_rows(P)
    return W, s


def part0(a):
    import scorer_proxy as SP
    import scorer_step16 as S16
    import step15
    res = dict(check1={}, check2={})
    ok_all = True
    for line in LINES_A:
        h = ref_for(line).h
        PC = (predside('K0', line), predside('K4', line))
        dP = float(np.abs(PC[0][0] - PC[1][0]).max())
        res.setdefault('pred_means_max_abs_diff_K0_K4', {})[line] = dP
        log(f'{line}: max |P_K0 - P_K4| = {dP:.3e}')
        for key, mask in (('K0', np.zeros(h.T, bool)), ('K4', np.ones(h.T, bool))):
            W, _ = combined(line, mask, PC)
            got, _ = Fast.from_wide(h, W).point()
            summ = json.loads((hdir(key, line) / 'summary.json').read_text())['scaled']
            want = {k: summ[m] for k, m in Fast.KEYS.items()}; want['overall'] = summ['overall']
            orig = results_wide(key, line).reindex(h.targets)
            perkd = {m: float(np.nanmax(np.abs(W[m].to_numpy(float) - orig[m].to_numpy(float))))
                     for m in ('expr_mse_unbiased_capped', 'expr_distance_unbiased', 'pds_cosine')}
            diff = {k: abs(got[k] - want[k]) for k in want}
            ok = all(v < 1e-6 for v in diff.values())
            ok_all &= ok
            res['check1'][f'{line}_{key}'] = dict(rebuilt=got, summary=want, abs_diff=diff, per_kd_max_abs_diff=perkd, pass_1e6=ok)
            log(f'check1 {line} all-{key}: overall {got["overall"]:.6f} vs {want["overall"]:.6f}; max comp diff '
                f'{max(diff.values()):.2e}; per-kd {perkd}; {"PASS" if ok else "FAIL"}')
    # check 2a: generator input (15c sec. 4: 0.2456); check 2b: fold cache (15a: 0.2470, halves 0.2679 / 0.2261)
    ref = ref_for('EvalA')
    kds, lfc, sig = SP.requested(C15 / 'config1p_K0/EvalA/input.npz'); assert kds == ref.h.targets
    m0 = pd.read_csv(C15 / 'config1p_K0/EvalA/e2e/controls/controls_baseline.csv', index_col=0)['control_mean_cpm'].reindex(ref.h.genes).to_numpy(float)
    g = ref.proxy(lfc, sig, m0)['overall']
    want_g = json.loads((C15 / 'config1p_K0/diag.json').read_text())['EvalA']['proxy_minus_actual']['overall']['proxy']
    res['check2']['generator_input'] = dict(proxy=g, diag_15c=want_g, rounded_report=0.2456, abs_diff=abs(g - want_g), pass_1e6=abs(g - want_g) < 1e-6)
    G = S16.Mag('H1'); T = __import__('scorer_step12').truth_calls(G.ref)
    out = G.ref.proxy(mu=G.build('i', step15.CONFIGS['config1p']), sig=G.fo.sig, per_kd=True)
    t15a = pd.read_csv(S15 / 'tables/step15a.csv')
    t15a = t15a[(t15a.line == 'H1') & (t15a.config == 'config1p') & (t15a.calls == 'K0')].set_index('part').overall
    for part, rows in (('full', None), ('half0', S16.halves(T)[0]), ('half1', S16.halves(T)[1])):
        v = G.ref.summarise(out['per_kd'], rows=rows)['overall']
        res['check2'][f'fold_cache_{part}'] = dict(proxy=v, step15a=float(t15a[part]), abs_diff=abs(v - t15a[part]), pass_1e6=abs(v - t15a[part]) < 1e-6)
    for k, v in res['check2'].items():
        ok_all &= v['pass_1e6']; log(f'check2 {k}: {v}')
    res['all_pass'] = bool(ok_all)
    save_json(TAB / 'part0.json', res)
    log(f'PART0 {"PASS" if ok_all else "FAIL"}')
    if not ok_all:
        sys.exit(2)


# ---------------------------------------------------------------- Part A
def b_for(line, targets):
    if line == 'EvalA':
        b = pd.read_csv(S15 / 'tables/step15b_predictors_H1.csv').set_index('target')['b']
    elif line == 'HCT116final':
        b = pd.read_csv(S17 / 'tables/b_HCT116.csv').set_index('target')['b']
    elif line == 'HEK293Tfinal':
        b = pd.read_csv(S17 / 'tables/b_HEK293T.csv').set_index('target')['b']
    elif line == 'K562final':
        b = pd.read_csv(S15 / 'tables/step15b_predictors_K562.csv').set_index('target')['b']
    return b.reindex(targets).to_numpy(float)


def assignments(b):
    """Ordered grid from K4-everywhere (lowest T) to K0-everywhere (highest T): label -> onK4 mask."""
    out = {'K4_all': np.ones(len(b), bool)}
    for T in T_GRID:
        out[f'T{T}'] = np.nan_to_num(b, nan=-np.inf) >= T
    out['K0_all'] = np.zeros(len(b), bool)
    return out


def requested_calls(key, line, targets):
    p = pd.read_parquet(C15 / CFG[key] / line / 'calls_per_kd.parquet').set_index('target')['requested']
    return p.reindex(targets).to_numpy(float)


def fit_pick(table, fit_part, order):
    """Highest overall on the fit part; ties go to the higher T (fewer K4 knockdowns). table: label -> {part: overall}."""
    best = max(order[::-1], key=lambda lab: table[lab][fit_part])     # max() keeps the first maximum in iteration order
    return best


def partA(a):
    import scorer_step12 as S12
    import scorer_step16 as S16
    rows, paired, res = [], [], dict(fits={}, terciles={}, distributions={})
    keep = {}
    for line in LINES_A:
        ref = ref_for(line); h = ref.h
        PC = (predside('K0', line), predside('K4', line))
        b = b_for(line, h.targets)
        A = assignments(b); order = list(A)
        req = {k: requested_calls(k, line, h.targets) for k in CFG}
        parts = {'full': np.arange(h.T)}
        Ttrue = S12.truth_calls(ref)
        if line == 'EvalA':
            hs = S16.halves(Ttrue); parts.update(half0=hs[0], half1=hs[1])
        I = boot_index(h.T)
        F = {}
        for lab, m in A.items():
            W, _ = combined(line, m, PC); F[lab] = Fast.from_wide(h, W)
            calls = np.where(m, req['K4'], req['K0'])
            for part, rr in parts.items():
                sc, raw = F[lab].point(rr)
                rows.append(dict(line=line, rule=lab, part=part, knockdowns=len(rr), share_K4=float(m[rr].mean()),
                                 calls_median=float(np.median(calls[rr])), calls_mean=float(np.mean(calls[rr])),
                                 **sc, **{f'{k}_raw': raw[k] for k in COMPS}))
            log(f'A {line} {lab}: share K4 {m.mean():.3f} full overall {rows[-len(parts)]["overall"]:.4f}')
        S0, _ = F['K0_all'].scaled(I)
        for lab in order:
            S1, _ = F[lab].scaled(I) if lab != 'K0_all' else (S0, None)
            p1 = F[lab].point()[0]; p0 = F['K0_all'].point()[0]
            for k in COMPS + ['overall']:
                d = S1[k] - S0[k]
                paired.append(dict(line=line, rule=lab, metric=k, diff=p1[k] - p0[k], ci_lo=float(np.quantile(d, 0.025)),
                                   ci_hi=float(np.quantile(d, 0.975)), p_better=float((d > 0).mean())))
            res.setdefault('ci_K0_overall', {})[line] = [float(np.quantile(S0['overall'], 0.025)), float(np.quantile(S0['overall'], 0.975))]
            log(f'A {line} {lab}: paired overall vs K0 {paired[-1]["diff"]:+.4f} [{paired[-1]["ci_lo"]:+.4f}, {paired[-1]["ci_hi"]:+.4f}]')
        keep[line] = dict(F=F, A=A, b=b, Ttrue=Ttrue, I=I, order=order, h=h)
    R = pd.DataFrame(rows); D = pd.DataFrame(paired)
    R.to_csv(TAB / 'partA_rows.csv', index=False); D.to_csv(TAB / 'partA_paired.csv', index=False)
    # fit on one H1 half, score the other
    E = keep['EvalA']; order = E['order']
    tab = {lab: R[(R.line == 'EvalA') & (R['rule'] == lab)].set_index('part').overall.to_dict() for lab in order}
    for fit_h in (0, 1):
        pick = fit_pick(tab, f'half{fit_h}', order)
        held = f'half{1 - fit_h}'
        r = R[(R.line == 'EvalA') & (R['rule'] == pick) & (R.part == held)].iloc[0]
        k0 = R[(R.line == 'EvalA') & (R['rule'] == 'K0_all') & (R.part == held)].iloc[0]
        hct = R[(R.line == 'HCT116final') & (R['rule'] == pick) & (R.part == 'full')].iloc[0]
        ci0 = res['ci_K0_overall']['HCT116final']
        hp = D[(D.line == 'HCT116final') & (D['rule'] == pick) & (D.metric == 'overall')].iloc[0]
        res['fits'][f'fit_half{fit_h}'] = dict(
            picked=pick, grid_position=grid_pos(order, pick), held_out_part=held, held_out=r[COMPS + ['overall']].to_dict(),
            K0_held_out=k0[COMPS + ['overall']].to_dict(), gain_vs_K0={k: float(r[k] - k0[k]) for k in COMPS + ['overall']},
            share_K4_held_out=float(r.share_K4), HCT116_overall=float(hct.overall), HCT116_vs_K0=float(hp['diff']),
            HCT116_vs_K0_ci=[float(hp.ci_lo), float(hp.ci_hi)], HCT116_K0_ci=ci0,
            HCT116_within_K0_ci=bool(ci0[0] <= hct.overall <= ci0[1]))
        log(f'A fit half{fit_h}: pick {pick} ({grid_pos(order, pick)}); held-out {held} {r.overall:.4f} vs K0 {k0.overall:.4f}')
    fits = res['fits']
    res['criterion'] = dict(heldout_gain_positive_both=bool(all(f['gain_vs_K0']['overall'] > 0 for f in fits.values())),
                            HCT116_within_K0_ci_both=bool(all(f['HCT116_within_K0_ci'] for f in fits.values())))
    res['criterion']['met'] = all(res['criterion'].values())
    # terciles at the picked T (15c sec. 1b terciles: rank thirds of the true DE count over the 292)
    nt = E['Ttrue'].sum(1); rk = pd.Series(nt).rank(method='first').to_numpy(); n = len(rk)
    tt = np.where(rk <= n / 3, 0, np.where(rk <= 2 * n / 3, 1, 2))
    for pick in sorted({f['picked'] for f in fits.values()}):
        F1, F0, m = E['F'][pick], E['F']['K0_all'], E['A'][pick]
        out = {}
        for j, t in enumerate(TERC):
            rr = np.flatnonzero(tt == j)
            rng = np.random.default_rng(SEED); Ij = np.stack([rr[rng.integers(0, len(rr), len(rr))] for _ in range(N_BOOT)])
            s1, _ = F1.scaled(Ij); s0, _ = F0.scaled(Ij)
            d = s1['overall'] - s0['overall']
            out[t] = dict(knockdowns=len(rr), share_K4=float(m[rr].mean()), true_de_median=float(np.median(nt[rr])),
                          overall_diff=float(F1.point(rr)[0]['overall'] - F0.point(rr)[0]['overall']),
                          ci=[float(np.quantile(d, 0.025)), float(np.quantile(d, 0.975))])
        res['terciles'][pick] = out
        log(f'A terciles at {pick}: {out}')
    # distributions: share of targets with b >= T
    bdef = pd.read_csv(S17 / 'tables/b_DEF.csv').set_index('target')['b']
    targets_def = pd.read_csv(V / 'data_final/extracted/pert_counts.csv')['target_gene'].astype(str).tolist()
    sets = {'DEF_300': bdef.reindex(targets_def).to_numpy(float), 'H1_panel': keep['EvalA']['b'],
            'HCT116_panel': keep['HCT116final']['b']}
    for nm, x in sets.items():
        res['distributions'][nm] = dict(n=len(x), n_with_b=int(np.isfinite(x).sum()),
                                        share_b_ge_T={T: float((np.nan_to_num(x, nan=-np.inf) >= T).mean()) for T in T_GRID})
    save_json(TAB / 'partA.json', res)
    with pd.option_context('display.width', 250, 'display.max_rows', 500):
        print(R.round(4).to_string()); print(D[D.metric == 'overall'].round(4).to_string())
    print(json.dumps(res, indent=1, default=float))


# ================================================================ Part B / C: real cells and controls
CTRL_SETS = {
    'D': (V / 'data_final/extracted/context_D.h5ad', None), 'E': (V / 'data_final/extracted/context_E.h5ad', None),
    'F': (V / 'data_final/extracted/context_F.h5ad', None),
    'EvalA_context': (EVAL['EvalA'] / 'context.h5ad', None), 'EvalA_scoring': (EVAL['EvalA'] / 'real.h5ad', 'non-targeting'),
    **{f'{l}_{k}': (EVAL[l] / f, lab) for l in ('HCT116final', 'K562final', 'HEK293Tfinal')
       for k, f, lab in (('context', 'context.h5ad', None), ('scoring', 'real.h5ad', 'non-targeting'))},
    'A': (V / 'data/context_A.h5ad', None), 'B': (V / 'data/context_B.h5ad', None), 'C': (V / 'data/context_C.h5ad', None),
}
STEP18_OFFSET = {'D': 0.060, 'E': 0.021, 'F': 0.025}


def control_stats(name):
    """Per control set: cells, per-cell UMI quantiles and CV, mean per-cell CPM (build_effects_v5.cpm_csr, as the model's
    control profile), detection fraction, pooled CPM, genes > 5 CPM, and the Step 18 baseline offset = median over
    gated genes of |log2(pooled CPM / mean per-cell CPM)|. Cached as tables/controls_<name>.npz."""
    import anndata as ad
    import scipy.sparse as sp
    import build_effects_v5 as B
    import run_line as RL
    f = TAB / f'controls_{name}.npz'
    if f.exists():
        z = np.load(f, allow_pickle=True); return {k: z[k] for k in z.files}
    path, lab = CTRL_SETS[name]
    a = ad.read_h5ad(path, backed='r')
    col = 'target' if 'target' in a.obs else 'target_gene'
    rows = np.arange(a.n_obs) if lab is None else np.flatnonzero(a.obs[col].astype(str).to_numpy() == lab)
    nnz_cell = a.file['X/data'].shape[0] / a.n_obs if 'X/data' in a.file else a.n_vars
    mem_before(f'controls {name} ({len(rows)} cells)', len(rows) * nnz_cell * 8 * 4 / 1e9 + 1)
    genes = RL.panel_genes()
    assert list(a.var_names.astype(str)) == genes, name
    X = a[rows].to_memory().X; a.file.close()
    X = sp.csr_matrix(X); X.sort_indices()
    lib = np.asarray(X.sum(1)).ravel().astype(float)
    Cm, _ = B.cpm_csr(X.data, X.indices, X.indptr, X.shape[1], np.ones(X.shape[1], bool))
    m0 = np.asarray(Cm.sum(0)).ravel() / X.shape[0]
    det = np.asarray((X > 0).sum(0)).ravel() / X.shape[0]
    raw = np.asarray(X.sum(0)).ravel().astype(float); pooled = 1e6 * raw / raw.sum()
    gate = m0 > B.MIN_CPM
    off = float(np.median(np.abs(np.log2(pooled[gate] / m0[gate]))))
    out = dict(n_cells=len(rows), lib_q=np.quantile(lib, [0.05, 0.25, 0.5, 0.75, 0.95]), lib_cv=lib.std() / lib.mean(),
               m0=m0, det=det, pooled=pooled, gate=gate, n_gated=int(gate.sum()), offset=off, min_cpm=B.MIN_CPM)
    np.savez_compressed(f, **out)
    mem_after(f'controls {name}')
    del X, Cm; gc.collect()
    return out


def de_counts(line):
    """True DE count per knockdown: real p_adj < 0.05, own gene out (as scorer_step12.truth_calls), from the eval's
    cached real DE table; and real cells per knockdown from real.h5ad's obs."""
    import glob
    import anndata as ad
    import polars as pl
    t = pl.read_parquet(glob.glob(str(EVAL[line] / 'real_cache' / 'de_wilcoxon_table-*.parquet'))[0], columns=['target', 'feature', 'p_adj'])
    targets = sorted(t['target'].unique().to_list())
    n = (t.filter((pl.col('p_adj') < 0.05) & (pl.col('feature') != pl.col('target'))).group_by('target').len()
         .to_pandas().set_index('target')['len'].reindex(targets).fillna(0).astype(int))
    a = ad.read_h5ad(EVAL[line] / 'real.h5ad', backed='r')
    col = 'target' if 'target' in a.obs else 'target_gene'
    cells = a.obs[col].astype(str).value_counts().reindex(targets).fillna(0).astype(int)
    a.file.close()
    return pd.DataFrame(dict(true_de=n, cells=cells, b=b_for(line, targets)), index=pd.Index(targets, name='target'))


def partial_spearman(y, x, z):
    """Spearman of y and x after removing rank(z) from both ranks (OLS residuals); rows with all three finite."""
    from scipy.stats import rankdata
    ok = np.isfinite(y) & np.isfinite(x) & np.isfinite(z)
    ry, rx, rz = (rankdata(v[ok]) for v in (y, x, z))
    Z = np.c_[np.ones(ok.sum()), rz]
    res = lambda v: v - Z @ np.linalg.lstsq(Z, v, rcond=None)[0]
    return float(np.corrcoef(res(ry), res(rx))[0, 1]), int(ok.sum())


def tercile_of(x):
    r = pd.Series(x).rank(method='first').to_numpy(); n = len(r)
    return np.where(r <= n / 3, 0, np.where(r <= 2 * n / 3, 1, 2))


def partB(a):
    from scipy.stats import spearmanr
    res = dict(knockdowns={}, controls={}, final_truth={})
    for line in ('EvalA', 'HCT116final', 'K562final', 'HEK293Tfinal'):
        d = de_counts(line); d.to_csv(TAB / f'partB_kd_{line}.csv')
        rho = spearmanr(d.true_de, d.cells).statistic
        pr, n_ok = partial_spearman(d.true_de.to_numpy(float), d.cells.to_numpy(float), np.log1p(d.b.to_numpy(float)))
        ok = d[np.isfinite(d.b)].copy()
        ok['b_terc'] = tercile_of(ok.b.to_numpy()); ok['cell_terc'] = -1
        for j in range(3):
            m = ok.b_terc == j; ok.loc[m, 'cell_terc'] = tercile_of(ok.loc[m, 'cells'].to_numpy())
        grid = ok.groupby(['b_terc', 'cell_terc']).true_de.median().unstack()
        res['knockdowns'][line] = dict(
            knockdowns=len(d), with_b=n_ok, cells_median=float(d.cells.median()), cells_q10=float(d.cells.quantile(0.1)),
            cells_q90=float(d.cells.quantile(0.9)), cells_min=int(d.cells.min()), cells_max=int(d.cells.max()),
            true_de_median=float(d.true_de.median()), spearman_true_cells=float(rho), partial_spearman_given_log1p_b=pr,
            spearman_true_b=float(spearmanr(ok.true_de, ok.b).statistic), spearman_cells_b=float(spearmanr(ok.cells, ok.b).statistic),
            median_true_by_b_tercile_x_cell_tercile={f'b{i}': {f'cells{j}': float(grid.loc[i, j]) for j in grid.columns} for i in grid.index})
        log(f'B1 {line}: {res["knockdowns"][line]}')
    for name in [n for n in CTRL_SETS if n not in ('A', 'B', 'C')]:
        s = control_stats(name)
        res['controls'][name] = dict(cells=int(s['n_cells']), umi_q05_q25_q50_q75_q95=[float(x) for x in s['lib_q']],
                                     umi_cv=float(s['lib_cv']), genes_above_5cpm=int(s['n_gated']), baseline_offset=float(s['offset']))
        log(f'B2 {name}: {res["controls"][name]}')
    chk = {x: dict(ours=res['controls'][x]['baseline_offset'], step18=v,
                   match_3dp=bool(round(res['controls'][x]['baseline_offset'], 3) == v)) for x, v in STEP18_OFFSET.items()}
    res['offset_check'] = chk; res['offset_check_pass'] = all(c['match_3dp'] for c in chk.values())
    log(f'B2 offset check vs Step 18: {chk}')
    man = json.loads((V / 'data_final/extracted/manifest.json').read_text())
    pc = pd.read_csv(V / 'data_final/extracted/pert_counts.csv')
    res['final_truth'] = dict(manifest={k: man[k] for k in ('cells_per_pert', 'n_constructs', 'per_context')},
                              pert_counts_columns=pc.columns.tolist(), pert_counts_rows=len(pc))
    log(f'B3 final truth: {res["final_truth"]}')
    save_json(TAB / 'partB.json', res)
    if not res['offset_check_pass']:
        log('B2 CHECK FAILED: baseline offsets do not reproduce Step 18'); sys.exit(3)


MARKERS = {
    'pluripotent': ['POU5F1', 'NANOG', 'LIN28A', 'DPPA4', 'TDGF1', 'DNMT3B', 'L1TD1', 'ZFP42'],
    'neural': ['SOX2', 'PAX6', 'NES', 'SOX1'],
    'blood': ['PTPRC', 'CD3E', 'CD3D', 'LCK', 'HBG1', 'HBG2', 'GATA1', 'TAL1', 'SPI1'],
    'epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1'],
    'mesenchymal': ['VIM', 'COL1A1', 'COL1A2', 'THY1'],
    'liver': ['ALB', 'APOA1', 'AFP'],
    'HEK293': ['LIN28B', 'NEFL', 'NEFM'],
    'sex': ['XIST', 'RPS4Y1', 'DDX3Y', 'KDM5D', 'UTY'],
    'proliferation_p53': ['MKI67', 'TOP2A', 'CDKN2A', 'TP53', 'MDM2', 'CDKN1A'],
}


def source_baselines():
    """(table, background) -> baseline DataFrame (control_mean_cpm, detection_rate) on the panel."""
    import stage1 as S
    import run_line as RL
    import h5py
    from anndata.io import read_elem
    genes = RL.panel_genes(); out = {}
    for p in S.discover_tables([str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')]):
        with h5py.File(p, 'r') as f:
            bg = str(read_elem(f['obs'])['background'].iloc[0])
        out[(p.stem, bg)] = pd.read_csv(p.with_name(p.stem + '_baseline.csv'), index_col=0).reindex(genes)
    return out


def partC(a):
    from scipy.stats import spearmanr
    import run_line as RL
    genes = RL.panel_genes(); gp = {g: j for j, g in enumerate(genes)}
    mk = [(grp, g) for grp, gs in MARKERS.items() for g in gs]
    missing = [g for _, g in mk if g not in gp]
    sets = ['D', 'E', 'F', 'EvalA_context', 'A', 'B', 'C']
    on_disk = {x: CTRL_SETS[x][0].exists() for x in ('A', 'B', 'C')}
    S = {x: control_stats(x) for x in sets if CTRL_SETS[x][0].exists()}
    rows = []
    for x, s in S.items():
        for grp, g in mk:
            j = gp.get(g)
            rows.append(dict(profile=x, kind='controls_h5ad', group=grp, gene=g, mean_cpm=float(s['m0'][j]) if j is not None else np.nan,
                             detection=float(s['det'][j]) if j is not None else np.nan))
    SB = source_baselines()
    for (tab, bg), d in SB.items():
        for grp, g in mk:
            rows.append(dict(profile=tab, kind=f'source_baseline:{bg}', group=grp, gene=g,
                             mean_cpm=float(d['control_mean_cpm'].get(g, np.nan)), detection=float(d['detection_rate'].get(g, np.nan))))
    M = pd.DataFrame(rows); M.to_csv(TAB / 'partC_markers.csv', index=False)
    # Spearman of log1p mean CPM over the D/E/F context's gated genes (> 5 CPM) with a finite value in the other profile
    cor = []
    for x in ('D', 'E', 'F'):
        g = S[x]['gate']; lx = np.log1p(S[x]['m0'])
        others = {y: (S[y]['m0'], 'controls_h5ad') for y in S if y != x}
        others.update({tab: (d['control_mean_cpm'].to_numpy(float), bg) for (tab, bg), d in SB.items()})
        for y, (v, kind) in others.items():
            ok = g & np.isfinite(v)
            if ok.sum() < 100:
                continue
            cor.append(dict(context=x, profile=y, kind=kind, genes=int(ok.sum()), spearman=float(spearmanr(lx[ok], np.log1p(v[ok])).statistic)))
    Cr = pd.DataFrame(cor); Cr.to_csv(TAB / 'partC_spearman.csv', index=False)
    wide = M[M.kind == 'controls_h5ad'].pivot_table(index=['group', 'gene'], columns='profile', values='mean_cpm', sort=False)
    det = M[M.kind == 'controls_h5ad'].pivot_table(index=['group', 'gene'], columns='profile', values='detection', sort=False)
    res = dict(abc_on_disk=on_disk, markers_missing_from_panel=missing,
               top_sources={x: Cr[Cr.context == x].sort_values('spearman', ascending=False).head(12).to_dict('records') for x in ('D', 'E', 'F')},
               among_contexts=Cr[Cr.kind == 'controls_h5ad'].to_dict('records'))
    save_json(TAB / 'partC.json', res)
    with pd.option_context('display.width', 250, 'display.max_rows', 500):
        print(wide.round(1).to_string()); print(det.round(3).to_string())
        print(Cr.sort_values(['context', 'spearman'], ascending=[True, False]).groupby('context').head(15).round(3).to_string())


# ================================================================ Part D: called-gene scale by expression tier (proxy)
S_TIER = [0.35, 0.5, 0.65, 0.8, 1.0]
CAPS = ['none', 2.0, 1.0]
D_LINES = ['H1', 'K562', 'HCT116', 'HEK293T']


def tiers(m0, gate):
    """-1 ungated; 0 / 1 / 2 = low / mid / high tercile of the line's own context-control mean CPM among gated genes."""
    t = np.full(len(m0), -1)
    gi = np.flatnonzero(gate); t[gi] = tercile_of(m0[gi])
    return t


def scaled_called(cg, sig, tier, s, cap):
    """Called gated genes: s[tier] x C+G LFC, capped at |cap| keeping the sign; returns (values, mask)."""
    sg = np.where(tier >= 0, np.asarray(s, float)[np.maximum(tier, 0)], np.nan)
    v = sg[None, :] * cg
    if cap != 'none':
        v = np.clip(v, -cap, cap)
    return v, sig & (tier >= 0)[None, :]


class DLine:
    """Part D inputs for one line. H1: the uploaded file's config (15c config 1' + K0 generator input, no-source genes
    zeroed). Guardrails: the fold pools (scorer-space fold caches, held out as before), config 1' + K0 as in 15a."""

    def __init__(self, line):
        import scorer_proxy as SP
        import scorer_step12 as S12
        import scorer_step16 as S16
        import step15
        self.line = line
        if line == 'H1':
            self.ref = ref_for('EvalA'); h = self.ref.h
            z = np.load(C15 / 'config1p_K0/EvalA/input.npz', allow_pickle=True)
            assert list(map(str, z['knockdowns'])) == h.targets
            un = z['n_sources'] <= 0
            self.lfc = np.where(un, 0.0, z['lfc'].astype(np.float64)); self.sig = z['de_sig'] & ~un
            self.cg = np.where(un, 0.0, z['cg_lfc'].astype(np.float64))
            cb = pd.read_csv(C15 / 'config1p_K0/EvalA/e2e/controls/controls_baseline.csv', index_col=0)
            self.m0 = cb['control_mean_cpm'].reindex(h.genes).to_numpy(float)
            self.gate = cb['de_tested'].reindex(h.genes).fillna(False).to_numpy(bool)
            self.base_mu = None
        else:
            G = S16.Mag(line); fo = G.fo; self.ref = G.ref; self.G = G
            p = step15.CONFIGS['config1p']
            self.sig, self.cg, self.m0, self.gate = fo.sig, fo.lfc, fo.m0, fo.gate
            self.mu_non = fo.mu_delta(np.where(fo.dn > 0, p['m'] * G.c_kappa(p['kappa'], p['c1'], 0.0) * fo.dpool, 0.0))
            self.base_mu = G.build('i', p)
        self.tier = tiers(self.m0, self.gate)
        T = S12.truth_calls(self.ref)
        self.parts = {'full': np.arange(self.ref.h.T)}
        if line == 'H1':
            hs = S16.halves(T); self.parts.update(half0=hs[0], half1=hs[1])

    def mu(self, s, cap):
        v, M = scaled_called(self.cg, self.sig, self.tier, s, cap)
        if self.line == 'H1':
            L = np.where(M, v, self.lfc)
            return np.maximum((self.m0[None, :] + 1e-9) * np.exp2(L) - 1e-9, 0.0)
        Lc = np.where(M, v, 0.5 * self.cg)
        return np.where(self.sig, self.G.fo.mu_lfc(Lc), self.mu_non)

    def score(self, s, cap):
        out = self.ref.proxy(mu=self.mu(s, cap), sig=self.sig, per_kd=True)
        rows = []
        for part, rr in self.parts.items():
            sm = self.ref.summarise(out['per_kd'], rows=rr)
            r = dict(line=self.line, s_low=s[0], s_mid=s[1], s_high=s[2], cap=str(cap), part=part, overall=sm['overall'])
            for m, k in Fast.KEYS.items():
                r[m] = sm['scaled'][k]; r[f'{m}_raw'] = sm['raw'][k]
            rows.append(r)
        return rows, out


def d_base_reference(line):
    """The value the base row must reproduce: H1 = the 15c generator-input proxy (Part 0); guardrails = 15a."""
    if line == 'H1':
        return json.loads((TAB / 'part0.json').read_text())['check2']['generator_input']['proxy']
    t = pd.read_csv(S15 / 'tables/step15a.csv')
    return float(t[(t.line == line) & (t.config == 'config1p') & (t.calls == 'K0') & (t.part == 'full')].overall.iloc[0])


def partD(a):
    """Grid worker: one line per worker (0 H1, 1 K562, 2 HCT116, 3 HEK293T); resumable."""
    from itertools import product
    line = D_LINES[a.worker]
    os.environ.setdefault('POLARS_MAX_THREADS', '1')
    f = TAB / f'partD_rows_{line}.parquet'
    rows = pd.read_parquet(f).to_dict('records') if f.exists() else []
    done = {(r['s_low'], r['s_mid'], r['s_high'], r['cap']) for r in rows}
    D = DLine(line)
    base = [r for r in rows if (r['s_low'], r['s_mid'], r['s_high'], r['cap']) == (0.5, 0.5, 0.5, 'none') and r['part'] == 'full']
    if not base:
        rr, _ = D.score((0.5, 0.5, 0.5), 'none'); rows += rr; done.add((0.5, 0.5, 0.5, 'none')); base = [rr[0]]
    want = d_base_reference(line); diff = abs(base[0]['overall'] - want)
    log(f'D {line}: base row {base[0]["overall"]:.6f} vs reference {want:.6f} |diff| {diff:.2e}')
    if diff > 1e-9:
        log(f'D {line}: BASE CHECK FAILED'); sys.exit(4)
    t0 = time.time(); n = 0
    for cap in CAPS:
        for s in product(S_TIER, S_TIER, S_TIER):
            if (s[0], s[1], s[2], str(cap)) in done:
                continue
            rr, _ = D.score(s, cap); rows += rr; n += 1
            if n % 10 == 0:
                pd.DataFrame(rows).to_parquet(f)
                log(f'D {line}: {len(rows) // len(D.parts)} rows, {(time.time() - t0) / n:.1f} s/row')
    pd.DataFrame(rows).to_parquet(f)
    log(f'D {line}: grid done, {len(rows) // len(D.parts)} rows')


def nmae_parts(ref, mu, sig, P, cats):
    """de_lfc_nmae's per-knockdown numerator sum|pred - real| and denominator sum|real| over its gate (real p_adj < 0.05,
    finite real LFC, own gene out; knockdowns with < 10 gate genes omitted), split by gene category.
    cats: name -> (T x G) bool. Returns per-category mean over knockdowns of (num_c / den_total) and (den_c / den_total),
    and the per-knockdown nmae (whose mean is the raw NMAE)."""
    import polars as pl
    pt = ref.pred_table(mu, sig, P)
    rl = ref.real_df['log2_fold_change'].to_numpy().astype(float); rp = ref.real_df['p_adj'].to_numpy().astype(float)
    pl_ = pt['log2_fold_change'].to_numpy().astype(float); pl_ = np.where(np.isfinite(pl_), pl_, 0.0)
    own = np.array([ref.h.own[i] for i in ref.ri]) == ref.rg
    g = (rp < 0.05) & np.isfinite(rl) & ~own
    T = ref.h.T
    cnt = np.bincount(ref.ri[g], minlength=T)
    den = np.bincount(ref.ri[g], weights=np.abs(rl[g]), minlength=T)
    num = np.bincount(ref.ri[g], weights=np.abs(pl_[g] - rl[g]), minlength=T)
    keep = (cnt >= 10) & (den > 0)
    out = {}
    for c, M in cats.items():
        mc = M[ref.ri, ref.rg] & g
        nc = np.bincount(ref.ri[mc], weights=np.abs(pl_[mc] - rl[mc]), minlength=T)
        dc = np.bincount(ref.ri[mc], weights=np.abs(rl[mc]), minlength=T)
        out[c] = dict(num_share_of_nmae=float((nc[keep] / den[keep]).mean()), den_share=float((dc[keep] / den[keep]).mean()),
                      gate_genes=int(mc[np.isin(ref.ri, np.flatnonzero(keep))].sum()))
    return out, num[keep] / den[keep], dict(pred_lfc=pl_, real_lfc=rl, real_p=rp, own=own)


def partD1(a):
    """D1 diagnostics on the base (the uploaded file's config on Eval A)."""
    import scorer_arms as SA
    import scorer_step12 as S12
    D = DLine('H1'); ref = D.ref; h = ref.h
    mu = D.mu((0.5, 0.5, 0.5), 'none'); P = SA.bulk(mu)
    prox = ref.proxy(mu=mu, sig=D.sig)
    gi = np.flatnonzero(D.gate); edges = np.quantile(D.m0[gi], np.linspace(0.1, 0.9, 9))
    dec = np.where(D.gate, np.searchsorted(edges, D.m0, side='right'), -1)                    # as scorer_proxy.call_diff
    DEC = np.broadcast_to(dec[None, :], (h.T, h.G))
    sig = D.sig
    cats = {f'decile_{k}': DEC == k for k in range(-1, 10)}
    cats.update(called=sig, uncalled=~sig)
    pt = ref.pred_table(mu, sig, P); plfc = np.zeros((h.T, h.G)); plfc[ref.ri, ref.rg] = pt['log2_fold_change'].to_numpy()
    ab = np.abs(plfc); bins = [(0, 0.1), (0.1, 0.25), (0.25, 0.5), (0.5, 1.0), (1.0, np.inf)]
    cats.update({f'pred_abs_lfc_{lo}-{hi}': (ab >= lo) & (ab < hi) for lo, hi in bins})
    cats.update({f'decile_{k}_{c}': (DEC == k) & (sig if c == 'called' else ~sig) for k in range(-1, 10) for c in ('called', 'uncalled')})
    parts, per_kd, arr = nmae_parts(ref, mu, sig, P, cats)
    chk = abs(per_kd.mean() - prox['raw'][Fast.KEYS['nmae']])
    log(f'D1 NMAE rebuild {per_kd.mean():.6f} vs proxy raw {prox["raw"][Fast.KEYS["nmae"]]:.6f} |diff| {chk:.2e}')
    if chk > 1e-9:
        log('D1 CHECK FAILED: NMAE split does not reproduce the proxy'); sys.exit(5)
    # sign accuracy on true DE genes (real p_adj < 0.05, finite non-zero real LFC, own gene out)
    td = (arr['real_p'] < 0.05) & np.isfinite(arr['real_lfc']) & (arr['real_lfc'] != 0) & ~arr['own']
    sg = sig[ref.ri, ref.rg]; dd = dec[ref.rg]
    right = np.sign(arr['pred_lfc']) == np.sign(arr['real_lfc'])
    sign = [dict(decile=k, called=c, genes=int((m := td & (dd == k) & (sg if c else ~sg)).sum()),
                 sign_accuracy=float(right[m].mean()) if m.any() else np.nan) for k in range(-1, 10) for c in (True, False)]
    # MSE plug-in error in scorer space, own gene out, by decile x called (raw-MSE units: sum_k sum_g d^2 / |G_k| / sum den)
    d2 = (P - h.R) ** 2
    for i in h.rows:
        if h.own[i] >= 0: d2[i, h.own[i]] = 0.0
    unit = d2 / h.nrow[:, None] / h.den.sum()
    mse = [dict(decile=k, called=c, plugin_raw_units=float(unit[(DEC == k) & (sig if c else ~sig)].sum()),
                share_of_plugin=float(unit[(DEC == k) & (sig if c else ~sig)].sum() / unit.sum())) for k in range(-1, 10) for c in (True, False)]
    res = dict(proxy_overall=prox['overall'], nmae_raw=float(per_kd.mean()), nmae_knockdowns=int(len(per_kd)), nmae=parts,
               sign=sign, mse_plugin=mse, plugin_total_raw_units=float(unit.sum()), mse_raw=prox['raw'][Fast.KEYS['mse']],
               decile_edges_cpm=edges.tolist(), calls_per_kd_median=float(np.median(sig.sum(1))))
    save_json(TAB / 'partD1.json', res)
    with pd.option_context('display.width', 250, 'display.max_rows', 500):
        print(pd.DataFrame(parts).T.round(4).to_string()); print(pd.DataFrame(sign).round(3).to_string()); print(pd.DataFrame(mse).round(4).to_string())


def partDsum(a):
    """Fit on one H1 half (rows with all six components > 0 on that half), score the other; surface; guardrails."""
    R = pd.concat([pd.read_parquet(TAB / f'partD_rows_{l}.parquet') for l in D_LINES])
    key = ['s_low', 's_mid', 's_high', 'cap']
    for l in D_LINES:
        n = R[(R.line == l) & (R.part == 'full')].drop_duplicates(key).shape[0]
        assert n == 375, (l, n)
    H = R[R.line == 'H1']
    isbase = lambda d: (d.s_low == 0.5) & (d.s_mid == 0.5) & (d.s_high == 0.5) & (d.cap == 'none')
    res = dict(fits={}, surface={}, guardrails={})
    for fit_h in (0, 1):
        fit = H[H.part == f'half{fit_h}']; held = f'half{1 - fit_h}'
        elig = fit[(fit[COMPS] > 0).all(1)]
        pick = elig.sort_values('overall', ascending=False).iloc[0]
        sel = lambda d: (d.s_low == pick.s_low) & (d.s_mid == pick.s_mid) & (d.s_high == pick.s_high) & (d.cap == pick.cap)
        ho = H[(H.part == held) & sel(H)].iloc[0]; b0 = H[(H.part == held) & isbase(H)].iloc[0]
        g = {}
        for l in D_LINES[1:]:
            L = R[(R.line == l) & (R.part == 'full')]
            x, y = L[sel(L)].iloc[0], L[isbase(L)].iloc[0]
            g[l] = dict(pick=x[COMPS + ['overall']].to_dict(), base=y[COMPS + ['overall']].to_dict(), diff=float(x.overall - y.overall))
        res['fits'][f'fit_half{fit_h}'] = dict(
            picked={k: pick[k] for k in key}, eligible_rows=int(len(elig)), fit_overall=float(pick.overall),
            grid_position={k: grid_pos(S_TIER, pick[k]) for k in key[:3]} | {'cap': grid_pos([str(c) for c in CAPS], pick.cap)},
            held_out=ho[COMPS + ['overall']].to_dict(), base_held_out=b0[COMPS + ['overall']].to_dict(),
            gain={k: float(ho[k] - b0[k]) for k in COMPS + ['overall']}, guardrails=g)
        log(f'D fit half{fit_h}: pick {res["fits"][f"fit_half{fit_h}"]["picked"]} held-out {ho.overall:.4f} vs base {b0.overall:.4f}')
    F = H[H.part == 'full']
    res['surface']['base_full'] = F[isbase(F)].iloc[0][COMPS + ['overall']].to_dict()
    res['surface']['best_per_cap'] = {c: F[F.cap == c].sort_values('overall').iloc[-1][key + COMPS + ['overall']].to_dict() for c in F.cap.unique()}
    res['surface']['best_per_cap_all_positive'] = {c: (lambda d: d.sort_values('overall').iloc[-1][key + COMPS + ['overall']].to_dict() if len(d) else None)(F[(F.cap == c) & (F[COMPS] > 0).all(1)]) for c in F.cap.unique()}
    res['surface']['marginal_mean_overall'] = {c: {k: F[F.cap == c].groupby(k).overall.mean().to_dict() for k in key[:3]} for c in F.cap.unique()}
    res['surface']['marginal_mean_nmae'] = {c: {k: F[F.cap == c].groupby(k).nmae.mean().to_dict() for k in key[:3]} for c in F.cap.unique()}
    fits = res['fits'].values()
    res['criterion'] = dict(heldout_gain_ge_0005_both=bool(all(f['gain']['overall'] >= 0.005 for f in fits)),
                            nmae_pos_both=bool(all(f['held_out']['nmae'] > 0 for f in fits)),
                            guardrails_within_0005=bool(all(g['diff'] >= -0.005 for f in fits for g in f['guardrails'].values())))
    res['criterion']['met'] = all(res['criterion'].values())
    save_json(TAB / 'partD.json', res)
    print(json.dumps(res, indent=1, default=float))


# ================================================================ Part E: program ceiling (PCA on source knockdowns)
E_K = [5, 10, 20, 50]
N_COMP = 100
HGNC_E = {   # Part E keyword rules on HGNC gene-group names (exact names found in hgnc_complete_set 2026-10-28)
    'histones': r'^H(1|2A|2B|3|4) histones$',
    'MCM_GINS_DNApol': r'^(MCM family|minichromosome maintenance 2-7 complex|GINS complex|DNA polymerases)$',
    'mitosis_kinetochore_APC_cyclins': r'^(Constitutive centromere associated network|MIS12 kinetochore complex|NDC80 kinetochore complex|Anaphase promoting complex|Cyclins)\s*$',
    'TP53_pathway': r'(p53|TP53)',
}


def e_tables():
    """Every source table except the Arc H1 ones: (path, stem, background key, is_tahoe)."""
    import h5py
    import stage1 as S
    from anndata.io import read_elem
    out = []
    for p in S.discover_tables([str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')]):
        with h5py.File(p, 'r') as f:
            o = read_elem(f['obs'])
        bg = str(o['background'].iloc[0])
        if bg == 'H1':
            continue
        tahoe = p.stem.startswith('tahoe_')
        out.append((p, p.stem, 'Tahoe' if tahoe else bg, tahoe))
    return out


def e_entries(lfc, se, gated, space, full=None, ctrl=None):
    """Raw basis entries for one block of rows (before centring), and the 'measured' mask.
    lfc space: LFC x z^2 / (z^2 + 1), z = LFC / SE. delta space: full-profile log1p CP50K minus the control profile."""
    if space == 'lfc':
        z2 = (lfc / se) ** 2
        x = lfc * z2 / (z2 + 1.0)
    else:
        x = full - ctrl[None, :]
    meas = np.isfinite(x) & gated[None, :]
    return np.where(meas, x, 0.0).astype(np.float32), meas


def e_build(variant):
    """The basis matrix (weighted rows) and its row metadata. variant: E1 / E2 (n_real >= 50) / E3 (20 CPM gate) /
    E4 (delta space). E5 is E1's rows minus Tahoe's (same centring and weights; built from E1)."""
    import h5py
    import run_line as RL
    from anndata.io import read_elem
    genes = RL.panel_genes(); G = len(genes)
    gate_cpm = 20.0 if variant == 'E3' else 5.0
    space = 'delta' if variant == 'E4' else 'lfc'
    plan = []
    for p, stem, bg, tahoe in e_tables():
        with h5py.File(p, 'r') as f:
            o = read_elem(f['obs'])
            if space == 'delta' and not np.isfinite(f['layers/log1p_bulk_50k_full'][:min(5, len(o))]).any():
                continue
        keep = np.ones(len(o), bool)
        if variant == 'E2':
            keep = np.isfinite(o['n_real'].astype(float).to_numpy()) & (o['n_real'].astype(float).to_numpy() >= 50)
        if keep.any():
            plan.append((p, stem, bg, tahoe, o, np.flatnonzero(keep)))
    n_rows = sum(len(k) for *_, k in plan)
    n_bg = pd.Series({}, dtype=float)
    for _, _, bg, _, _, k in plan:
        n_bg[bg] = n_bg.get(bg, 0) + len(k)
    mem_before(f'{variant} basis {n_rows} x {G}', n_rows * G * 4 / 1e9 + 3)
    M = np.zeros((n_rows, G), np.float32); meta = []; r0 = 0
    for p, stem, bg, tahoe, o, k in plan:
        b = pd.read_csv(p.with_name(stem + '_baseline.csv'), index_col=0).reindex(genes)
        gated = (b['control_mean_cpm'].to_numpy(float) > gate_cpm)
        ctrl = b['log1p_bulk_50k_full'].to_numpy(float) if space == 'delta' else None
        blocks, meass = [], []
        with h5py.File(p, 'r') as f:
            for c0 in range(0, len(k), 2000):
                kk = k[c0:c0 + 2000]
                sl = slice(kk[0], kk[-1] + 1)
                if space == 'lfc':
                    l_ = f['layers/lfc'][sl].astype(np.float64)[kk - kk[0]]; s_ = f['layers/se'][sl].astype(np.float64)[kk - kk[0]]
                    x, m = e_entries(l_, s_, gated, 'lfc')
                else:
                    fu = f['layers/log1p_bulk_50k_full'][sl].astype(np.float64)[kk - kk[0]]
                    x, m = e_entries(None, None, gated, 'delta', fu, ctrl)
                blocks.append(x); meass.append(m)
        X = np.concatenate(blocks); Mm = np.concatenate(meass)
        mean = np.where(Mm, X, 0.0).sum(0) / np.maximum(Mm.sum(0), 1)                      # measured entries only
        X = np.where(Mm, X - mean[None, :], 0.0) * (1.0 / np.sqrt(n_bg[bg]))
        M[r0:r0 + len(k)] = X; r0 += len(k)
        nr = o['n_real'].astype(float).to_numpy()[k]
        meta.append(pd.DataFrame(dict(table=stem, background=bg, tahoe=tahoe, modality=o['modality'].astype(str).to_numpy()[k],
                                      target=o['target_gene'].astype(str).to_numpy()[k], n_real=nr, weight=1.0 / np.sqrt(n_bg[bg]))))
        del X, Mm, blocks, meass
    mem_after(f'{variant} basis')
    return M, pd.concat(meta, ignore_index=True), n_bg


def e_svd(M, label):
    from sklearn.utils.extmath import randomized_svd
    t0 = time.time()
    U, s, Vt = randomized_svd(M, n_components=N_COMP, n_oversamples=20, n_iter=5, random_state=SEED)
    tot = float(np.einsum('ij,ij->', M, M, dtype=np.float64))
    cum = np.cumsum(s.astype(np.float64) ** 2) / tot
    k25 = int(np.searchsorted(cum, 0.25) + 1) if cum[-1] >= 0.25 else None
    log(f'E {label}: SVD {M.shape} in {time.time() - t0:.0f} s; cum share ' +
        ', '.join(f'K{K} {cum[K - 1]:.3f}' for K in E_K + [N_COMP] if K <= len(cum)) + f'; K at 25% {k25}')
    return (U * s[None, :]).astype(np.float32), s, Vt.astype(np.float32), cum, k25


def e_truth():
    """H1 (Eval A, 292 knockdowns): true LFC on the real DE table's tested genes, true Delta (harness real cache), the
    true-DE (NMAE gate) mask, H1's control CP50K and the 15c true-strength terciles."""
    import scorer_step12 as S12
    ref = ref_for('EvalA'); h = ref.h
    X = np.full((h.T, h.G), np.nan)
    rl = ref.real_df['log2_fold_change'].to_numpy().astype(float); rp = ref.real_df['p_adj'].to_numpy().astype(float)
    X[ref.ri, ref.rg] = rl
    DE = np.zeros((h.T, h.G), bool); DE[ref.ri, ref.rg] = (rp < 0.05) & np.isfinite(rl)
    for i in h.rows:
        if h.own[i] >= 0: X[i, h.own[i]] = np.nan; DE[i, h.own[i]] = False
    gated = np.isfinite(X)
    Dt = h.R - h.ctrl[None, :]
    nt = S12.truth_calls(ref).sum(1)
    return dict(X=X, gated=gated, DE=DE, Delta=Dt, v=np.expm1(h.ctrl), terc=tercile_of(nt), n_true=nt, h=h, ref=ref)


def e_ceiling(Vt, K, tr, space='lfc'):
    """Per knockdown least-squares fit of its truth on the top-K basis rows over its gated genes; returns the per-
    knockdown squared-error and energy sums for (a) all gated genes, (b) its true-DE genes (coefficients from (a)),
    (c) Delta space (LFC fit converted gene by gene with H1's control CP50K; Delta-space bases fit Delta directly)."""
    T = tr['X'].shape[0]
    out = {k: np.zeros((T, 2)) for k in 'abc'}
    B = Vt[:K].astype(np.float64)
    for i in range(T):
        g = tr['gated'][i]
        if space == 'lfc':
            x = tr['X'][i, g]; A = B[:, g].T
            c, *_ = np.linalg.lstsq(A, x, rcond=None); xh = A @ c
            out['a'][i] = [((x - xh) ** 2).sum(), (x ** 2).sum()]
            de = tr['DE'][i, g]; out['b'][i] = [((x - xh)[de] ** 2).sum(), (x[de] ** 2).sum()]
            v = tr['v'][g]; dh = np.log1p(v * np.exp2(xh)) - np.log1p(v); dt = tr['Delta'][i, g]
            out['c'][i] = [((dt - dh) ** 2).sum(), (dt ** 2).sum()]
        else:
            dt = tr['Delta'][i, g]; A = B[:, g].T
            c, *_ = np.linalg.lstsq(A, dt, rcond=None); dh = A @ c
            out['c'][i] = [((dt - dh) ** 2).sum(), (dt ** 2).sum()]
            de = tr['DE'][i, g]; out['b'][i] = [((dt - dh)[de] ** 2).sum(), (dt[de] ** 2).sum()]
            out['a'][i] = [np.nan, np.nan]
    return out


def captured(sums, terc):
    """1 - sum err / sum energy over all knockdowns and per true-strength tercile."""
    r = {}
    for nm, rr in (('all', slice(None)),) + tuple((t, terc == j) for j, t in enumerate(TERC)):
        e = sums[rr]; r[nm] = float(1 - e[:, 0].sum() / e[:, 1].sum()) if np.isfinite(e).all() else None
    return r


FOLD_DROP = {'E1-K562': ('K562',), 'E1-Orion': ('HCT116', 'HEK293T'), 'E5': ('Tahoe',)}   # row subsets of E1


def partE(a):
    """Bases E1 (main), E1-K562 / E1-Orion (E-reg fold bases), E5, E2, E3, E4; ceilings (E1-E5); programs; anchors.
    The subsets drop whole backgrounds, so per-table centring and the other backgrounds' 1/sqrt(n_b) weights are
    unchanged: each is exactly E1 built without those tables. For E1 and the fold bases, every E1 row's unweighted
    program scores (row / weight @ Vt^T) are saved for E-reg."""
    import scorer_step11 as S11
    import re
    tr = e_truth(); res = dict(bases={}, ceiling=[], programs={}, targets={})
    E1_store = None
    for variant in ('E1', 'E1-K562', 'E1-Orion', 'E5', 'E2', 'E3', 'E4'):
        if variant in FOLD_DROP:
            M1, meta1, nb1 = E1_store
            keep = ~meta1.background.isin(FOLD_DROP[variant]).to_numpy()
            mem_before(f'{variant}: copy of E1 rows ({keep.sum()} of {len(keep)})', keep.sum() * M1.shape[1] * 4 / 1e9 + 1)
            M = M1[keep]; meta = meta1[keep].reset_index(drop=True); n_bg = nb1.drop(list(FOLD_DROP[variant]))
        else:
            M, meta, n_bg = e_build(variant)
        US, s, Vt, cum, k25 = e_svd(M, variant)
        if variant == 'E1':
            E1_store = (M, meta, n_bg)
        else:
            del M; gc.collect()
        if variant in ('E1', 'E1-K562', 'E1-Orion'):
            M1, meta1, _ = E1_store
            RS = np.empty((len(meta1), N_COMP), np.float32); w = meta1.weight.to_numpy(np.float32)
            for r0 in range(0, len(meta1), 10000):
                RS[r0:r0 + 10000] = (M1[r0:r0 + 10000] / w[r0:r0 + 10000, None]) @ Vt.T
            np.save(TAB / f'partE_{variant}_E1rows_scores.npy', RS); del RS
        if variant == 'E5':
            E1_store = None; gc.collect(); mem_after('E1 freed')
        np.savez_compressed(TAB / f'partE_{variant}_basis.npz', Vt=Vt, s=s, cum=cum, US=US)
        meta.to_parquet(TAB / f'partE_{variant}_rows.parquet')
        res['bases'][variant] = dict(rows=int(len(meta)), rows_per_background={k: int(v) for k, v in n_bg.items()},
                                     tables=int(meta.table.nunique()), cum_share={K: float(cum[K - 1]) for K in E_K + [N_COMP]},
                                     K_at_25pct=k25)
        if variant == 'E1':
            res['bases']['E1']['tahoe_share_of_row_score_energy'] = {
                c + 1: float((US[meta.tahoe.to_numpy(), c].astype(np.float64) ** 2).sum() / (US[:, c].astype(np.float64) ** 2).sum()) for c in range(10)}
        if variant in ('E1-K562', 'E1-Orion'):
            continue                                                       # fold bases: no ceiling (E-reg only)
        space = 'delta' if variant == 'E4' else 'lfc'
        for K in sorted(set(E_K + [k25 or N_COMP])):
            sums = e_ceiling(Vt, K, tr, space)
            res['ceiling'].append(dict(basis=variant, K=K, K_is_25pct=K == (k25 or N_COMP), **{f'{m}_{g}': v for m in 'abc' for g, v in captured(sums[m], tr['terc']).items()}))
            log(f'E ceiling {variant} K {K}: {res["ceiling"][-1]}')
        if variant == 'E1':
            # programs: top 10 components
            import run_line as RL
            genes = np.array(RL.panel_genes())
            for c in range(10):
                o = np.argsort(Vt[c]); sc = US[:, c]
                top = np.argsort(-sc)[:12]; bot = np.argsort(sc)[:12]
                res['programs'][c + 1] = dict(
                    var_share=float(cum[c] - (cum[c - 1] if c else 0.0)),
                    pos_genes=genes[o[::-1][:20]].tolist(), neg_genes=genes[o[:20]].tolist(),
                    top_pos_rows=[f'{meta.table[i]}:{meta.target[i]}' for i in top], top_neg_rows=[f'{meta.table[i]}:{meta.target[i]}' for i in bot])
            # D/E/F targets: mean unweighted score (weighted row score / row weight) on comps 1-5 over genetic source tables
            gen = (meta.modality != 'drug').to_numpy()
            raw = US[:, :5] / meta.weight.to_numpy()[:, None]
            p90 = np.quantile(np.abs(raw[gen]), 0.9, axis=0)
            targets = pd.read_csv(V / 'data_final/extracted/pert_counts.csv')['target_gene'].astype(str).tolist()
            dg = pd.DataFrame(raw[gen], columns=[f'pc{c + 1}' for c in range(5)]); dg['target'] = meta.target.to_numpy()[gen]
            dg['table'] = meta.table.to_numpy()[gen]
            ms = dg[dg.target.isin(targets)].groupby(['target', 'table']).mean().groupby('target').mean()   # per table, then over tables
            hi = ms[(ms.abs() > p90[None, :]).any(axis=1)]
            bdef = pd.read_csv(S17 / 'tables/b_DEF.csv').set_index('target')['b']
            grp = S11.hgnc_groups()
            res['targets'] = dict(p90_abs_score=p90.tolist(), targets_with_score=int(len(ms)), n_above_p90_any=int(len(hi)),
                                  above=[dict(target=t, **{k: float(v) for k, v in r.items()}, b=float(bdef.get(t, np.nan)),
                                              hgnc='; '.join(grp.get(t, []))) for t, r in hi.iterrows()])
            ms.to_csv(TAB / 'partE_DEF_target_scores.csv')
            log(f'E targets: {len(hi)} of {len(ms)} D/E/F targets above p90 on any of PCs 1-5')
    # anchors: random orthonormal (mean of 5 draws), H1's own top-K PCs (in sample), config 1' (15c generator input)
    rng = np.random.default_rng(SEED); G = tr['X'].shape[1]
    for K in E_K:
        acc = []
        for d in range(5):
            Q, _ = np.linalg.qr(rng.standard_normal((G, K)))
            acc.append({m: captured(v, tr['terc']) for m, v in e_ceiling(Q.T, K, tr).items()})
        res['ceiling'].append(dict(basis='random', K=K, **{f'{m}_{g}': float(np.mean([x[m][g] for x in acc])) for m in 'abc' for g in ('all',) + TERC}))
        log(f'E ceiling random K {K}: {res["ceiling"][-1]}')
    X0 = np.nan_to_num(tr['X'])
    _, _, Vh = np.linalg.svd(X0, full_matrices=False)
    for K in E_K:
        sums = e_ceiling(Vh, K, tr)
        res['ceiling'].append(dict(basis='H1_own_PCs', K=K, **{f'{m}_{g}': v for m in 'abc' for g, v in captured(sums[m], tr['terc']).items()}))
    import scorer_arms as SA
    D = DLine('H1'); mu = D.mu((0.5, 0.5, 0.5), 'none'); P = SA.bulk(mu)
    pt = D.ref.pred_table(mu, D.sig, P); pl_ = np.full_like(tr['X'], np.nan); pl_[D.ref.ri, D.ref.rg] = pt['log2_fold_change'].to_numpy()
    sums = {k: np.zeros((tr['X'].shape[0], 2)) for k in 'abc'}
    for i in range(tr['X'].shape[0]):
        g = tr['gated'][i]; x = tr['X'][i, g]; xh = np.nan_to_num(pl_[i, g]); de = tr['DE'][i, g]
        dt = tr['Delta'][i, g]; dh = (P - tr['h'].ctrl[None, :])[i, g]
        sums['a'][i] = [((x - xh) ** 2).sum(), (x ** 2).sum()]; sums['b'][i] = [((x - xh)[de] ** 2).sum(), (x[de] ** 2).sum()]
        sums['c'][i] = [((dt - dh) ** 2).sum(), (dt ** 2).sum()]
    res['ceiling'].append(dict(basis="config1p_prediction", K=None, **{f'{m}_{g}': v for m in 'abc' for g, v in captured(sums[m], tr['terc']).items()}))
    res['S2_cited_H1_delta'] = {5: 0.043, 10: 0.087, 20: 0.125, 50: 0.213}
    # HGNC keyword counts among the 300 D/E/F targets
    grp = S11.hgnc_groups(); targets = pd.read_csv(V / 'data_final/extracted/pert_counts.csv')['target_gene'].astype(str).tolist()
    allg = sorted({x for v in grp.values() for x in v})
    res['hgnc_rules'] = {k: [g for g in allg if re.search(rx, g)] for k, rx in HGNC_E.items()}
    res['hgnc_counts'] = {k: sorted(t for t in targets if set(grp.get(t, [])) & set(gs)) for k, gs in res['hgnc_rules'].items()}
    C = pd.DataFrame(res['ceiling']); C.to_csv(TAB / 'partE_ceiling.csv', index=False)
    st = C[(C.basis == 'E1') & (C.K <= 20)]
    cp = C[C.basis == 'config1p_prediction'].iloc[0]
    res['criterion'] = dict(E1_strong_b_by_K={int(r.K): r.b_strong for r in st.itertuples()}, E1_strong_c_by_K={int(r.K): r.c_strong for r in st.itertuples()},
                            config1p_strong_b=cp.b_strong, config1p_strong_c=cp.c_strong, S2_range=[0.043, 0.213])
    save_json(TAB / 'partE.json', res)
    with pd.option_context('display.width', 250, 'display.max_rows', 500):
        print(C.round(4).to_string())
    print(json.dumps({k: v for k, v in res.items() if k != 'ceiling'}, indent=1, default=float))


# ================================================================ Part D addendum: source-support gate rows
GATE_ROWS = {'G1': (5.0, 0.5), 'G2': (5.0, 0.0), 'G3': (20.0, 0.5), 'G4': (20.0, 0.0)}
OWN_BG = {'H1': ('H1',), 'K562': ('K562',), 'HCT116': ('HCT116', 'HEK293T'), 'HEK293T': ('HCT116', 'HEK293T')}


def source_support_tables():
    """Per genetic background: control CPM on the panel (its tables' baseline control_mean_cpm averaged), and the
    knockdowns it holds (genetic rows only: Tahoe's drug rows carry a drug-target gene but are not knockdowns of it)."""
    import h5py
    import stage1 as S
    import run_line as RL
    from anndata.io import read_elem
    genes = RL.panel_genes(); cpm, holds = {}, {}
    for p in S.discover_tables([str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')]):
        with h5py.File(p, 'r') as f:
            o = read_elem(f['obs'])
        o = o[o['modality'].astype(str) != 'drug']
        if not len(o):
            continue
        bg = str(o['background'].iloc[0])
        b = pd.read_csv(p.with_name(p.stem + '_baseline.csv'), index_col=0)['control_mean_cpm'].reindex(genes).to_numpy(float)
        cpm.setdefault(bg, []).append(b); holds.setdefault(bg, set()).update(o['target_gene'].astype(str))
    with np.errstate(all='ignore'):
        cpm = {bg: np.nanmean(np.stack(v), 0) for bg, v in cpm.items()}
    return cpm, holds


def support_values(line, targets, cpm, holds):
    """(T x G) median over the backgrounds holding each knockdown (the fold's own lines excluded) of their control CPM;
    NaN where no holding background measures the gene."""
    bgs = [bg for bg in cpm if bg not in OWN_BG[line]]
    C = np.stack([cpm[bg] for bg in bgs])
    out = np.full((len(targets), C.shape[1]), np.nan)
    with np.errstate(all='ignore'):
        for i, t in enumerate(targets):
            rows = [j for j, bg in enumerate(bgs) if t in holds[bg]]
            if rows:
                out[i] = np.nanmedian(C[rows], 0)
    return out


def partDgate(a):
    """G1-G4: called LFC x factor where the source-support value is below the threshold; base = Part D's base."""
    import scorer_arms as SA
    cpm, holds = source_support_tables()
    rows, touch, res = [], {}, dict(lines={})
    D1 = {}
    for line in D_LINES:
        D = DLine(line); h = D.ref.h
        sup = support_values(line, h.targets, cpm, holds)
        called = D.sig
        mag = np.abs(0.5 * D.cg)                                                          # called |LFC| at the base
        info = dict(called_pairs=int(called.sum()), called_without_support_value=int((called & ~np.isfinite(sup)).sum()))
        for thr in (5.0, 20.0):
            below = called & ~(sup >= thr)                                                # NaN support counts as below
            info[f'share_called_below_{thr:g}'] = float(below.sum() / called.sum())
            info[f'share_called_abs_lfc_below_{thr:g}'] = float(mag[below].sum() / mag[called].sum())
            touch[(line, thr)] = below
        res['lines'][line] = info
        log(f'G {line}: {info}')
        for name, (thr, fac) in [('base', (None, 1.0))] + list(GATE_ROWS.items()):
            mult = np.where(touch[(line, thr)], fac, 1.0) if thr else np.ones_like(D.cg)
            if line == 'H1':
                mu = SA.mu_from_lfc(np.where(D.sig, D.lfc * mult, D.lfc), D.m0)
            else:
                mu = np.where(D.sig, D.G.fo.mu_lfc(0.5 * D.cg * mult), D.mu_non)
            out = D.ref.proxy(mu=mu, sig=D.sig, per_kd=True)
            for part, rr in D.parts.items():
                sm = D.ref.summarise(out['per_kd'], rows=rr)
                r = dict(line=line, row=name, part=part, overall=sm['overall'])
                for m, k in Fast.KEYS.items():
                    r[m] = sm['scaled'][k]
                rows.append(r)
            log(f'G {line} {name}: full overall {rows[-len(D.parts)]["overall"]:.4f}')
        if line == 'H1':
            D1 = dict(D=D, sup=sup)
        else:
            del D; gc.collect()
    R = pd.DataFrame(rows)
    base = R[R.row == 'base'].set_index(['line', 'part'])
    for c in COMPS + ['overall']:
        R[f'd_{c}'] = R[c].to_numpy() - base.loc[list(zip(R.line, R.part)), c].to_numpy()
    want = d_base_reference('H1'); bh = float(base.loc[('H1', 'full'), 'overall'])
    res['base_check_H1'] = dict(base=bh, reference=want, abs_diff=abs(bh - want))
    assert abs(bh - want) < 1e-9, res['base_check_H1']
    R.to_csv(TAB / 'partDgate_rows.csv', index=False)
    # D1 placement of the touched genes on Eval A
    D = D1['D']; ref = D.ref; h = ref.h
    mu = D.mu((0.5, 0.5, 0.5), 'none'); P = SA.bulk(mu)
    gi = np.flatnonzero(D.gate); edges = np.quantile(D.m0[gi], np.linspace(0.1, 0.9, 9))
    dec = np.where(D.gate, np.searchsorted(edges, D.m0, side='right'), -1)
    cats = {f'below_{thr:g}': touch[('H1', thr)] for thr in (5.0, 20.0)}
    cats['called'] = D.sig
    parts, _, arr = nmae_parts(ref, mu, D.sig, P, cats)
    td = (arr['real_p'] < 0.05) & np.isfinite(arr['real_lfc']) & (arr['real_lfc'] != 0) & ~arr['own']
    right = np.sign(arr['pred_lfc']) == np.sign(arr['real_lfc'])
    place = {}
    for nm, M in cats.items():
        m = M[ref.ri, ref.rg] & td
        place[nm] = dict(nmae=parts[nm], true_de_pairs=int(m.sum()), sign_accuracy=float(right[m].mean()) if m.any() else None,
                         decile_distribution={int(k): float(v) for k, v in pd.Series(np.broadcast_to(dec, M.shape)[M]).value_counts(normalize=True).sort_index().items()})
    res['D1_placement_EvalA'] = place
    crit = {}
    for nm in GATE_ROWS:
        h0 = R[(R.line == 'H1') & (R.row == nm) & (R.part == 'half0')].iloc[0]; h1_ = R[(R.line == 'H1') & (R.row == nm) & (R.part == 'half1')].iloc[0]
        g = R[(R.line != 'H1') & (R.row == nm)]
        crit[nm] = dict(gain_half0=float(h0.d_overall), gain_half1=float(h1_.d_overall), nmae_half0=float(h0.nmae), nmae_half1=float(h1_.nmae),
                        worst_guardrail_diff=float(g.d_overall.min()))
        crit[nm]['met'] = bool(min(crit[nm]['gain_half0'], crit[nm]['gain_half1']) >= 0.005 and min(crit[nm]['nmae_half0'], crit[nm]['nmae_half1']) > 0
                               and crit[nm]['worst_guardrail_diff'] >= -0.005)
    res['criterion'] = crit
    save_json(TAB / 'partDgate.json', res)
    with pd.option_context('display.width', 250, 'display.max_rows', 200):
        print(R.round(4).to_string())
    print(json.dumps({k: v for k, v in res.items()}, indent=1, default=float))


# ================================================================ Part E-reg: program regression, scored on the proxy
EREG_SCALES = [0.25, 0.5, 0.75, 1.0]
EREG_BASIS = {'H1': 'E1', 'K562': 'E1-K562', 'HCT116': 'E1-Orion', 'HEK293T': 'E1-Orion'}


def ereg_khat(RS, meta, targets, exclude):
    """s_hat (T x 100): per knockdown, genetic background votes (tables averaged first, then backgrounds), each the
    projection of the unweighted centred confidence-weighted row; 0 where no background holds it. Also the vote count."""
    ok = (meta.modality != 'drug').to_numpy() & ~meta.background.isin(exclude).to_numpy() & meta.target.isin(targets).to_numpy()
    d = pd.DataFrame(RS[ok].astype(np.float64)); d['table'] = meta.table.to_numpy()[ok]
    d['background'] = meta.background.to_numpy()[ok]; d['target'] = meta.target.to_numpy()[ok]
    per_t = d.groupby(['target', 'background', 'table']).mean(numeric_only=True)
    per_bg = per_t.groupby(level=['target', 'background']).mean()
    sh = per_bg.groupby(level='target').mean().reindex(targets)
    nb = per_bg.groupby(level='target').size().reindex(targets).fillna(0).astype(int)
    return np.nan_to_num(sh.to_numpy()), nb.to_numpy()


def e_coefs(Vt, K, X, gated):
    """Least-squares coefficients of each row of X (NaN = ungated) on the top-K basis rows over its gated genes."""
    B = Vt[:K].astype(np.float64); out = np.zeros((X.shape[0], K))
    for i in range(X.shape[0]):
        g = gated[i]; out[i] = np.linalg.lstsq(B[:, g].T, X[i, g], rcond=None)[0]
    return out


def fit_beta(s, sh, rows):
    num = (s[rows] * sh[rows]).sum(0); den = (sh[rows] ** 2).sum(0)
    return np.where(den > 0, num / np.where(den > 0, den, 1), 0.0)


def energy_sums(tr, Lhat):
    """(a), (b), (c) sums for a predicted LFC matrix (H1 truth's gated genes)."""
    T = tr['X'].shape[0]; out = {k: np.zeros((T, 2)) for k in 'abc'}
    for i in range(T):
        g = tr['gated'][i]; x = tr['X'][i, g]; xh = Lhat[i, g]; de = tr['DE'][i, g]
        v = tr['v'][g]; dh = np.log1p(v * np.exp2(xh)) - np.log1p(v); dt = tr['Delta'][i, g]
        out['a'][i] = [((x - xh) ** 2).sum(), (x ** 2).sum()]; out['b'][i] = [((x - xh)[de] ** 2).sum(), (x[de] ** 2).sum()]
        out['c'][i] = [((dt - dh) ** 2).sum(), (dt ** 2).sum()]
    return out


def ereg_mu(D, P, arm, scale):
    """R1: uncalled gated genes = scale x P_prog; R2: config 1' + scale x P_prog on all gated genes. LFC space, converted
    to the generator's means as config 1''s called values are (mu = (m0 + 1e-9) 2^LFC - 1e-9)."""
    import scorer_arms as SA
    Lb = D.lfc_base
    g = D.gate[None, :]
    if arm == 'R1':
        L = np.where(~D.sig & g, scale * P, Lb)
    else:
        L = np.where(g, Lb + scale * P, Lb)
    return SA.mu_from_lfc(L, D.m0)


def partEreg(a):
    import scorer_arms as SA
    E = json.loads((TAB / 'partE.json').read_text())
    tr = e_truth(); H = tr['h']
    terc = tr['terc']
    import scorer_step12 as S12
    import scorer_step16 as S16
    halves = S16.halves(S12.truth_calls(tr['ref']))
    res = dict(K={}, prediction=[], programs=[], arms={}, fits={}, guardrails={})
    bases = {}
    for b in ('E1', 'E1-K562', 'E1-Orion'):
        z = np.load(TAB / f'partE_{b}_basis.npz')
        k25 = E['bases'][b]['K_at_25pct'] or N_COMP
        bases[b] = dict(Vt=z['Vt'], RS=np.load(TAB / f'partE_{b}_E1rows_scores.npy'), Ks=sorted(set([5, 10, 20, int(k25)])), k25=int(k25))
        res['K'][b] = dict(K_grid=bases[b]['Ks'], K_at_25pct=E['bases'][b]['K_at_25pct'])
    meta = pd.read_parquet(TAB / 'partE_E1_rows.parquet')
    C = pd.read_csv(TAB / 'partE_ceiling.csv')
    # ---- H1: s, s_hat, beta per half, cross-fitted prediction report
    b1 = bases['E1']
    sh_all, nb = ereg_khat(b1['RS'], meta, H.targets, OWN_BG['H1'])
    res['H1_knockdowns_with_votes'] = int((nb > 0).sum())
    P_dir = {}                                                         # (K, fit half) -> P for all 292 with beta_fit
    for K in b1['Ks']:
        s = e_coefs(b1['Vt'], K, tr['X'], tr['gated']); sh = sh_all[:, :K]
        betas = {f: fit_beta(s, sh, halves[f]) for f in (0, 1)}
        for f in (0, 1):
            ho = halves[1 - f]; bta = betas[f]
            for j in range(K):
                x, y = sh[ho, j], s[ho, j]
                res['programs'].append(dict(K=K, fit_half=f, program=j + 1, beta=float(bta[j]),
                                            r=float(np.corrcoef(x, y)[0, 1]) if x.std() > 0 and y.std() > 0 else np.nan,
                                            r2_heldout=float(1 - ((y - bta[j] * x) ** 2).sum() / (y ** 2).sum())))
            P_dir[(K, f)] = (bta[None, :] * sh) @ b1['Vt'][:K].astype(np.float64)
        Pcf = np.empty_like(P_dir[(K, 0)])                                 # cross-fitted: each half from the other's beta
        Pcf[halves[1]] = P_dir[(K, 0)][halves[1]]; Pcf[halves[0]] = P_dir[(K, 1)][halves[0]]
        np.save(TAB / f'partEreg_P_H1_K{K}_crossfit.npy', Pcf.astype(np.float32))      # for Part F's combination rows
        es = energy_sums(tr, Pcf)
        cl = C[(C.basis == 'E1') & (C.K == K)]
        res['prediction'].append(dict(K=K, **{f'{m}_{g}': v for m in 'abc' for g, v in captured(es[m], terc).items()},
                                      ceiling={c: float(cl[c].iloc[0]) for c in cl.columns if c[:2] in ('a_', 'b_', 'c_')} if len(cl) else None))
        log(f'E-reg H1 K {K}: held-out energy {res["prediction"][-1]}')
    cp = C[C.basis == 'config1p_prediction'].iloc[0]
    res['config1p_energy'] = {c: float(cp[c]) for c in C.columns if c[:2] in ('a_', 'b_', 'c_')}
    pd.DataFrame(res['programs']).to_csv(TAB / 'partEreg_programs.csv', index=False)
    # ---- proxy arms
    rows = []

    def add_rows(D, P, K, tag, extra):
        for arm in ('R1', 'R2'):
            for sc in EREG_SCALES:
                out = D.ref.proxy(mu=ereg_mu(D, P, arm, sc), sig=D.sig, per_kd=True)
                for part, rr in D.parts.items():
                    sm = D.ref.summarise(out['per_kd'], rows=rr)
                    r = dict(line=D.line, arm=arm, scale=sc, K=K, beta_from=tag, part=part, overall=sm['overall'], **extra)
                    for m, k in Fast.KEYS.items():
                        r[m] = sm['scaled'][k]
                    rows.append(r)
        pd.DataFrame(rows).to_csv(TAB / 'partEreg_rows.csv', index=False)

    def prep(D):
        if D.line == 'H1':
            D.lfc_base = D.lfc
        else:
            D.lfc_base = np.log2((D.base_mu + 1e-9) / (D.m0[None, :] + 1e-9))
        mu0 = SA.mu_from_lfc(D.lfc_base, D.m0)
        out = D.ref.proxy(mu=mu0, sig=D.sig, per_kd=True)
        base = {part: D.ref.summarise(out['per_kd'], rows=rr) for part, rr in D.parts.items()}
        want = d_base_reference(D.line); got = base['full']['overall']
        log(f'E-reg {D.line}: base via LFC round trip {got:.6f} vs {want:.6f} |diff| {abs(got - want):.1e}')
        assert abs(got - want) < 1e-6
        return {p: dict(overall=v['overall'], **{m: v['scaled'][k] for m, k in Fast.KEYS.items()}) for p, v in base.items()}

    D = DLine('H1'); res['base'] = {'H1': prep(D)}
    for K in b1['Ks']:
        for f in (0, 1):
            add_rows(D, P_dir[(K, f)] * D.gate[None, :], K, f'half{f}', dict(fit_half=f))
        log(f'E-reg H1 K {K}: rows done')
    del D; gc.collect()
    # guardrails: beta on all 292 H1 under the fold's basis (fold lines excluded from the H1 votes too)
    for line in D_LINES[1:]:
        bn = EREG_BASIS[line]; bb = bases[bn]
        D = DLine(line); res['base'][line] = prep(D)
        shH, _ = ereg_khat(bb['RS'], meta, H.targets, OWN_BG['H1'] + OWN_BG[line])
        shL, nbL = ereg_khat(bb['RS'], meta, D.ref.h.targets, OWN_BG[line])
        res.setdefault('fold_votes', {})[line] = int((nbL > 0).sum())
        for K in bb['Ks']:
            s = e_coefs(bb['Vt'], K, tr['X'], tr['gated'])
            bta = fit_beta(s, shH[:, :K], np.arange(H.T))
            P = (bta[None, :] * shL[:, :K]) @ bb['Vt'][:K].astype(np.float64) * D.gate[None, :]
            np.save(TAB / f'partEreg_P_{line}_K{K}.npy', P.astype(np.float32))
            add_rows(D, P, K, 'all_H1', dict(fit_half=-1))
        log(f'E-reg {line}: rows done'); del D; gc.collect()
    R = pd.DataFrame(rows)
    # ---- fit on one H1 half (rows fitted with that half's beta), score the other
    Hh = R[R.line == 'H1']
    for f in (0, 1):
        fit = Hh[(Hh.fit_half == f) & (Hh.part == f'half{f}')]
        elig = fit[(fit[COMPS] > 0).all(1)]
        if not len(elig):
            res['fits'][f'fit_half{f}'] = dict(picked=None, note='no row with all six components > 0 on the fit half'); continue
        pk = elig.sort_values('overall').iloc[-1]
        ho = Hh[(Hh.fit_half == f) & (Hh.part == f'half{1 - f}') & (Hh.arm == pk.arm) & (Hh.scale == pk.scale) & (Hh.K == pk.K)].iloc[0]
        b0 = res['base']['H1'][f'half{1 - f}']
        g = {}
        for line in D_LINES[1:]:
            L = R[(R.line == line) & (R.arm == pk.arm) & (R.scale == pk.scale)]
            Kg = pk.K if pk.K in set(L.K) else (bases[EREG_BASIS[line]]['k25'] if pk.K == b1['k25'] else pk.K)
            x = L[L.K == Kg].iloc[0]
            g[line] = dict(K=int(Kg), overall=float(x.overall), base=res['base'][line]['full']['overall'],
                           diff=float(x.overall - res['base'][line]['full']['overall']), **{c: float(x[c]) for c in COMPS})
        res['fits'][f'fit_half{f}'] = dict(
            picked=dict(arm=pk.arm, scale=float(pk.scale), K=int(pk.K)), eligible=int(len(elig)),
            grid_position=dict(arm=grid_pos(['R1', 'R2'], pk.arm), scale=grid_pos(EREG_SCALES, pk.scale), K=grid_pos(b1['Ks'], pk.K)),
            held_out={c: float(ho[c]) for c in COMPS + ['overall']}, base_held_out=b0,
            gain={c: float(ho[c] - b0[c]) for c in COMPS + ['overall']}, guardrails=g)
        log(f'E-reg fit half{f}: {res["fits"][f"fit_half{f}"]["picked"]} held-out gain {res["fits"][f"fit_half{f}"]["gain"]["overall"]:+.4f}')
    fv = [v for v in res['fits'].values() if v.get('picked')]
    res['criterion'] = dict(both_directions_picked=len(fv) == 2,
                            heldout_gain_ge_0005_both=bool(len(fv) == 2 and all(v['gain']['overall'] >= 0.005 for v in fv)),
                            all_components_pos_both=bool(len(fv) == 2 and all(all(v['held_out'][c] > 0 for c in COMPS) for v in fv)),
                            guardrails_within_0005=bool(len(fv) == 2 and all(g['diff'] >= -0.005 for v in fv for g in v['guardrails'].values())))
    res['criterion']['met'] = all(res['criterion'].values())
    # ---- D/E/F: P_prog for the 300 targets per context (all sources incl. Arc H1; E1 basis; beta on all of H1)
    targets = pd.read_csv(V / 'data_final/extracted/pert_counts.csv')['target_gene'].astype(str).tolist()
    shD, nbD = ereg_khat(b1['RS'], meta, targets, ())
    arc = arc_h1_scores(b1['Vt'], targets)                                  # Arc H1 as one more background
    has = np.isfinite(arc).all(1)
    shD = np.where(has[:, None], (shD * nbD[:, None] + np.nan_to_num(arc)) / np.maximum(nbD + has, 1)[:, None], shD)
    nbD = nbD + has
    out = dict(targets=np.array(targets), n_backgrounds=nbD)
    for K in b1['Ks']:
        s = e_coefs(b1['Vt'], K, tr['X'], tr['gated'])
        bta = fit_beta(s, sh_all[:, :K], np.arange(H.T)); out[f'beta_K{K}'] = bta
        P = (bta[None, :] * shD[:, :K]) @ b1['Vt'][:K].astype(np.float64)
        for x in ('D', 'E', 'F'):
            out[f'P_{x}_K{K}'] = (P * control_stats(x)['gate'][None, :]).astype(np.float32)
    np.savez_compressed(TAB / 'partEreg_DEF_Pprog.npz', **out)
    res['DEF'] = dict(file=str(TAB / 'partEreg_DEF_Pprog.npz'), targets_with_votes=int((nbD > 0).sum()), median_backgrounds=float(np.median(nbD)))
    save_json(TAB / 'partEreg.json', res)
    with pd.option_context('display.width', 250, 'display.max_rows', 400):
        print(R.round(4).to_string())
    print(json.dumps({k: v for k, v in res.items() if k != 'programs'}, indent=1, default=float))


def arc_h1_scores(Vt, targets):
    """Arc H1's confidence-weighted centred rows (5 CPM gate on its own controls), projected on Vt; NaN rows where
    Arc H1 lacks the target. Duplicate rows of a target are averaged."""
    import h5py
    import run_line as RL
    from anndata.io import read_elem
    p = V / 'lib/effects_v5/arc_h1.h5ad'; genes = RL.panel_genes()
    with h5py.File(p, 'r') as f:
        o = read_elem(f['obs']); l_ = f['layers/lfc'][:].astype(np.float64); s_ = f['layers/se'][:].astype(np.float64)
    b = pd.read_csv(p.with_name('arc_h1_baseline.csv'), index_col=0).reindex(genes)
    X, Mm = e_entries(l_, s_, b['control_mean_cpm'].to_numpy(float) > 5.0, 'lfc')
    mean = np.where(Mm, X, 0.0).sum(0) / np.maximum(Mm.sum(0), 1)
    X = np.where(Mm, X - mean[None, :], 0.0)
    sc = pd.DataFrame(X @ Vt.T.astype(np.float64)); sc['target'] = o['target_gene'].astype(str).to_numpy()
    return sc.groupby('target').mean().reindex(targets).to_numpy()


# ================================================================ Part F: combination rows on the proxy
def combo_mu(D, spec, sup=None):
    """Config 1' + K0 base (Part D's), then in order: tier row (spec['tier'] = [s_low, s_mid, s_high, cap]), gate row
    (spec['gate'] = [threshold, factor]), the no-source fix (H1: gated no-source genes keep config 1''s pooled-profile
    LFC, i.e. n_sources >= 1), and an E-reg arm (spec['ereg'] = [arm, scale, K]). The magnitude split is K0's."""
    import scorer_arms as SA
    if D.line == 'H1':
        L = D.lfc.copy()
        if spec.get('nosource_fix'):
            z = np.load(C15 / 'config1p_K0/EvalA/input.npz', allow_pickle=True)
            un = (z['n_sources'] <= 0) & D.gate[None, :]
            L = np.where(un, z['lfc'].astype(np.float64), L)
    else:
        L = np.log2((D.base_mu + 1e-9) / (D.m0[None, :] + 1e-9))       # fold pools already carry config 1''s value
    called = np.where(D.sig, 0.5 * D.cg, 0.0)
    if spec.get('tier'):
        s_l, s_m, s_h, cap = spec['tier']
        v, M = scaled_called(D.cg, D.sig, D.tier, (s_l, s_m, s_h), cap if cap == 'none' else float(cap))
        called = np.where(M, v, called)
    if spec.get('gate'):
        thr, fac = spec['gate']
        called = np.where(D.sig & ~(sup >= thr), fac * called, called)
    L = np.where(D.sig, called, L)
    if spec.get('ereg'):
        arm, sc, K = spec['ereg']
        f = TAB / (f'partEreg_P_H1_K{K}_crossfit.npy' if D.line == 'H1' else f'partEreg_P_{D.line}_K{K}.npy')
        P = np.load(f).astype(np.float64) * D.gate[None, :]
        g = D.gate[None, :]
        L = np.where(~D.sig & g, sc * P, L) if arm == 'R1' else np.where(g, L + sc * P, L)
    return SA.mu_from_lfc(L, D.m0)


def partF(a):
    """Score one combination spec (JSON) on H1 (full, halves) and the three fold pools, next to the base."""
    import scorer_step16 as S16
    spec = json.loads(a.spec); name = spec.get('name', a.spec)
    f = TAB / 'partF_rows.csv'
    rows = pd.read_csv(f).to_dict('records') if f.exists() else []
    sup_t = source_support_tables() if spec.get('gate') else None
    for line in D_LINES:
        D = DLine(line)
        sup = support_values(line, D.ref.h.targets, *sup_t) if sup_t else None
        sig = D.sig
        if spec.get('calls') in ('K4', 'hybrid'):
            if line == 'H1':
                z = np.load(C15 / 'config1p_K0/EvalA/input.npz', allow_pickle=True); k4 = z['de_sig_K4'] & ~(z['n_sources'] <= 0)
            else:
                k4 = S16.call_sets(D.G, [1200])['K4_1200']
            if spec['calls'] == 'K4':
                sig = k4
            else:                                   # Part A's rule on the proxy: K4 1,200 where b >= T, else K0
                b = b_for({'H1': 'EvalA', 'K562': 'K562final', 'HCT116': 'HCT116final', 'HEK293T': 'HEK293Tfinal'}[line], D.ref.h.targets)
                on = np.nan_to_num(b, nan=-np.inf) >= spec['T']
                sig = np.where(on[:, None], k4, D.sig)
        for nm, sp in (('base', {}), (name, spec)):
            if nm == 'base' and any(r['line'] == line and r['row'] == 'base' for r in rows):
                continue
            out = D.ref.proxy(mu=combo_mu(D, sp, sup), sig=sig if nm != 'base' else D.sig, per_kd=True)
            for part, rr in D.parts.items():
                sm = D.ref.summarise(out['per_kd'], rows=rr)
                r = dict(line=line, row=nm, spec=json.dumps(sp), part=part, overall=sm['overall'])
                for m, k in Fast.KEYS.items():
                    r[m] = sm['scaled'][k]
                rows.append(r)
            log(f'F {line} {nm}: full {rows[-len(D.parts)]["overall"]:.4f}')
        del D; gc.collect()
        pd.DataFrame(rows).to_csv(f, index=False)
    if spec.get('nosource_fix') and spec.get('genin_check'):
        # the fix on HCT116final's generator input (15c config1p_K0), where the fold pool cannot show it
        import scorer_proxy as SP
        ref = ref_for('HCT116final'); h = ref.h
        z = np.load(C15 / 'config1p_K0/HCT116final/input.npz', allow_pickle=True)
        names = list(map(str, z['knockdowns'])); r = [names.index(t) for t in h.targets]
        cb = pd.read_csv(C15 / 'config1p_K0/HCT116final/e2e/controls/controls_baseline.csv', index_col=0)
        m0 = cb['control_mean_cpm'].reindex(h.genes).to_numpy(float); gate = cb['de_tested'].reindex(h.genes).fillna(False).to_numpy(bool)
        un = z['n_sources'][r] <= 0; lf = z['lfc'][r].astype(np.float64); sg = z['de_sig'][r] & ~un
        out = {}
        for nm, L in (('genin_base', np.where(un, 0.0, lf)), ('genin_fix', np.where(un & ~gate[None, :], 0.0, lf))):
            pr = ref.proxy(L, sg, m0)
            out[nm] = dict(overall=pr['overall'], **{m: pr['scaled'][k] for m, k in Fast.KEYS.items()})
        save_json(TAB / 'partF_nosource_HCT116final_genin.json', out)
        log(f'F HCT116final generator input: base {out["genin_base"]["overall"]:.4f} fix {out["genin_fix"]["overall"]:.4f}')


# ================================================================ Follow-up 2026-11-02: MSE by true-strength tercile, residual programs
def mse_terc_rows(h, num, den, terc):
    """Raw MSE (ratio of sums) per tercile and overall, scaled with the full panel's anchors (as boot_score), clamped
    and unclamped."""
    out = {}
    for nm, rr in (('all', np.arange(h.T)),) + tuple((t, np.flatnonzero(terc == j)) for j, t in enumerate(TERC)):
        u = float(num[rr].sum() / den[rr].sum())
        out[nm] = dict(raw=u, scaled=h.scale(Fast.KEYS['mse'], u), unclamped=float((h.bm - u) / (h.bm - h.rm)),
                       num_share=float(num[rr].sum() / num.sum()), den_share=float(den[rr].sum() / den.sum()))
    return out


def mseterc(a):
    import mse_study as M
    import scorer_arms as SA
    import scorer_step12 as S12
    import scorer_step16 as S16
    ref = ref_for('EvalA'); h = ref.h
    terc = tercile_of(S12.truth_calls(ref).sum(1)); halves = S16.halves(S12.truth_calls(ref))
    res = dict(mse={}, residual={}, programs={}, own_pcs={})
    # ---- arms: (P, C) in scorer space
    arms = {}
    z = np.load(V / 'results/stage1/2026-10-28_scorer_space/tables/pred_cells_CG_H1.npz', allow_pickle=True)
    assert list(map(str, z['targets'])) == h.targets
    arms['C+G (cells)'] = (z['P'], z['C'])
    z = np.load(TAB / 'predside_K0_EvalA.npz'); arms["config 1' + K0, uploaded (cells)"] = (z['P'], z['C'])
    D = DLine('H1'); D.lfc_base = D.lfc
    arms["config 1' + K0 (proxy, generator input)"] = (SA.bulk(D.mu((0.5, 0.5, 0.5), 'none')), ref.Cpred)
    for sc in (0.5, 0.75):
        P = np.load(TAB / 'partEreg_P_H1_K100_crossfit.npy').astype(np.float64) * D.gate[None, :]
        arms[f'E-reg R2 gamma {sc} K100, cross-fitted (proxy)'] = (SA.bulk(ereg_mu(D, P, 'R2', sc)), ref.Cpred)
    nums = {}
    for nm, (P, C) in arms.items():
        _, num = M.score_generated(h, P, C); nums[nm] = num
        res['mse'][nm] = mse_terc_rows(h, num, h.den, terc)
    b = pd.read_csv(EVAL['EvalA'] / 'baseline/baseline_results.csv').pivot(index='perturbation', columns='metric', values='value').reindex(h.targets)
    res['mse']['bundle baseline'] = mse_terc_rows(h, b['expr_mse_unbiased_capped'].to_numpy(float), b['expr_distance_unbiased'].to_numpy(float), terc)
    chk = abs(res['mse']["config 1' + K0, uploaded (cells)"]['all']['raw'] - json.loads((hdir('K0', 'EvalA') / 'summary.json').read_text())['raw'][Fast.KEYS['mse']])
    log(f'mseterc: uploaded cells full raw MSE reproduces summary.json |diff| {chk:.1e}'); assert chk < 1e-9
    for nm, r in res['mse'].items():
        log(f'MSE {nm}: ' + ', '.join(f'{t} {v["scaled"]:.3f} (unclamped {v["unclamped"]:+.3f}, raw {v["raw"]:.3f})' for t, v in r.items()))
    # ---- residual of the uploaded config's cells
    P, C = arms["config 1' + K0, uploaded (cells)"]
    Rr = P - h.R
    for i in h.rows:
        if h.own[i] >= 0: Rr[i, h.own[i]] = 0.0
    cb = pd.read_csv(C15 / 'config1p_K0/EvalA/e2e/controls/controls_baseline.csv', index_col=0)
    gate = cb['de_tested'].reindex(h.genes).fillna(False).to_numpy(bool)
    E2 = Rr ** 2
    for j, t in enumerate(TERC):
        rr = terc == j; tot = E2[rr].sum()
        m = Rr[rr].mean(0)
        res['residual'][t] = dict(
            plugin_raw_units=float((E2[rr].sum(1) / h.nrow[rr]).sum() / h.den[rr].sum()),
            share_called=float(E2[rr][D.sig[rr]].sum() / tot), share_uncalled_gated=float(E2[rr][~D.sig[rr] & gate[None, :]].sum() / tot),
            share_ungated=float(E2[rr][:, ~gate].sum() / tot),
            share_tercile_mean=float(rr.sum() * (m ** 2).sum() / tot),
            share_panel_mean=float(rr.sum() * (Rr.mean(0) ** 2).sum() / tot))
    # ---- program capture of the residual (perfect coefficients) and the MSE it would leave
    bases = {'E4 (delta)': np.load(TAB / 'partE_E4_basis.npz')['Vt'].astype(np.float64),
             'E1 (LFC)': np.load(TAB / 'partE_E1_basis.npz')['Vt'].astype(np.float64)}
    rng = np.random.default_rng(SEED); Q, _ = np.linalg.qr(rng.standard_normal((h.G, 50))); bases['random'] = Q.T
    cap = []
    for bn, Vt in bases.items():
        for K in [k for k in (5, 10, 20, 50, 96, 100) if k <= Vt.shape[0]]:
            Cf = Rr @ Vt[:K].T; Rh = Cf @ Vt[:K]
            for i in h.rows:
                if h.own[i] >= 0: Rh[i, h.own[i]] = 0.0
            row = dict(basis=bn, K=K)
            for j, t in enumerate(TERC):
                rr = terc == j; row[f'captured_{t}'] = float(1 - ((Rr[rr] - Rh[rr]) ** 2).sum() / E2[rr].sum())
            _, num2 = M.score_generated(h, P - Rh, C)
            o = mse_terc_rows(h, num2, h.den, terc)
            row.update({f'mse_{t}': o[t]['scaled'] for t in ('all',) + TERC}); row.update({f'mse_unclamped_{t}': o[t]['unclamped'] for t in ('all',) + TERC})
            cap.append(row); log(f'capture {bn} K{K}: {row}')
    res['capture'] = cap
    # ---- which E4 programs carry the mid + strong residual
    import run_line as RL
    genes = np.array(RL.panel_genes())
    Vt = bases['E4 (delta)']; Cf = Rr @ Vt.T
    ms = terc >= 1
    tot = E2[ms].sum(); share = (Cf[ms] ** 2).sum(0) / tot
    meta4 = pd.read_parquet(TAB / 'partE_E4_rows.parquet'); US4 = np.load(TAB / 'partE_E4_basis.npz')['US']
    blog = np.log1p(b_for('EvalA', h.targets))
    from scipy.stats import spearmanr
    top = np.argsort(-share)[:10]
    for c in top:
        o = np.argsort(Vt[c]); x = Cf[ms, c]
        okb = np.isfinite(blog[ms])
        res['programs'][int(c) + 1] = dict(
            residual_share_mid_strong=float(share[c]), share_strong=float((Cf[terc == 2, c] ** 2).sum() / E2[terc == 2].sum()),
            share_weak=float((Cf[terc == 0, c] ** 2).sum() / E2[terc == 0].sum()),
            mean_coef_mid_strong=float(x.mean()), share_positive=float((x > 0).mean()),
            mean_sq_over_sq_mean=float((x ** 2).mean() / max(x.mean() ** 2, 1e-30)),
            spearman_coef_vs_log1p_b=float(spearmanr(x[okb], blog[ms][okb]).statistic),
            pos_genes=genes[o[::-1][:15]].tolist(), neg_genes=genes[o[:15]].tolist(),
            top_source_rows=[f'{meta4.table[i]}:{meta4.target[i]}' for i in np.argsort(-np.abs(US4[:, c]))[:8]])
    # ---- H1's own residual programs (mid + strong): fit on one half, captured on the other; are they in the sources?
    E4span = bases['E4 (delta)']; E1span = bases['E1 (LFC)']
    for f in (0, 1):
        fit = np.intersect1d(halves[f], np.flatnonzero(ms)); ho = np.intersect1d(halves[1 - f], np.flatnonzero(ms))
        _, sv, W = np.linalg.svd(Rr[fit], full_matrices=False)
        rows_ = []
        for K in (1, 2, 5, 10, 20):
            Rh = (Rr[ho] @ W[:K].T) @ W[:K]
            rows_.append(dict(K=K, captured_heldout=float(1 - ((Rr[ho] - Rh) ** 2).sum() / (Rr[ho] ** 2).sum()),
                              captured_fit=float((sv[:K] ** 2).sum() / (Rr[fit] ** 2).sum())))
        pcs = []
        for c in range(5):
            w = W[c]; o = np.argsort(w)
            pcs.append(dict(pc=c + 1, in_E4_span=float(((E4span @ w) ** 2).sum()), in_E1_span=float(((E1span @ w) ** 2).sum()),
                            max_cos_E4=float(np.abs(E4span @ w).max()), heldout_mean_coef=float((Rr[ho] @ w).mean()),
                            heldout_share_positive=float(((Rr[ho] @ w) > 0).mean() if (Rr[ho] @ w).mean() > 0 else ((Rr[ho] @ w) < 0).mean()),
                            pos_genes=genes[o[::-1][:15]].tolist(), neg_genes=genes[o[:15]].tolist()))
        res['own_pcs'][f'fit_half{f}'] = dict(capture=rows_, pcs=pcs)
        log(f'own residual PCs fit half {f}: {rows_}')
    save_json(TAB / 'mse_terciles.json', res)
    pd.DataFrame(cap).to_csv(TAB / 'mse_terciles_capture.csv', index=False)
    print(json.dumps(res, indent=1, default=float))


def mseterc2(a):
    """Held-out MSE of corrections learned on one H1 half and applied to the other (uploaded config's cells):
    (i) the fit half's mean residual per tercile, assigned by true tercile (oracle) or by (b) tercile (source-only);
    (ii) the fit half's top-K residual PCs (mid + strong) with oracle per-knockdown coefficients (a ceiling);
    (iii) those PCs with each tercile's mean fit-half coefficient (a shift along them, by true or (b) tercile)."""
    import mse_study as M
    import scorer_step12 as S12
    import scorer_step16 as S16
    ref = ref_for('EvalA'); h = ref.h
    T = S12.truth_calls(ref); terc = tercile_of(T.sum(1)); halves = S16.halves(T)
    b = b_for('EvalA', h.targets); bt = np.full(h.T, 0)
    ok = np.isfinite(b); bt[ok] = tercile_of(b[ok])                              # no b -> weak
    z = np.load(TAB / 'predside_K0_EvalA.npz'); P, C = z['P'], z['C']
    Rr = P - h.R
    for i in h.rows:
        if h.own[i] >= 0: Rr[i, h.own[i]] = 0.0
    rows = []

    def score(name, f, Rh):
        _, num = M.score_generated(h, P - Rh, C)
        ho = halves[1 - f]
        for t in ('all',) + TERC:
            rr = ho if t == 'all' else np.intersect1d(ho, np.flatnonzero(terc == TERC.index(t)))
            u = num[rr].sum() / h.den[rr].sum()
            rows.append(dict(arm=name, fit_half=f, part=t, raw=float(u), scaled=h.scale(Fast.KEYS['mse'], u), unclamped=float((h.bm - u) / (h.bm - h.rm))))

    for f in (0, 1):
        fit, ho = halves[f], halves[1 - f]
        Z = np.zeros_like(Rr); score('uploaded (no correction)', f, Z)
        for lab, tt in (('true tercile', terc), ('(b) tercile', bt)):
            Rh = np.zeros_like(Rr)
            for j in (1, 2):
                mfit = Rr[np.intersect1d(fit, np.flatnonzero(terc == j))].mean(0)     # shift learned on true terciles
                Rh[np.intersect1d(ho, np.flatnonzero(tt == j))] = mfit
            score(f'tercile-mean shift (mid, strong), assigned by {lab}', f, Rh)
        ms_fit = np.intersect1d(fit, np.flatnonzero(terc >= 1))
        _, _, W = np.linalg.svd(Rr[ms_fit], full_matrices=False)
        for K in (1, 2, 5, 10, 20):
            Wk = W[:K]; Rh = np.zeros_like(Rr); hm = np.intersect1d(ho, np.flatnonzero(terc >= 1))
            Rh[hm] = (Rr[hm] @ Wk.T) @ Wk
            score(f'own residual PCs K{K}, oracle coefficients (mid, strong)', f, Rh)
            for lab, tt in (('true tercile', terc), ('(b) tercile', bt)):
                Rh = np.zeros_like(Rr)
                for j in (1, 2):
                    cbar = (Rr[np.intersect1d(fit, np.flatnonzero(terc == j))] @ Wk.T).mean(0)
                    Rh[np.intersect1d(ho, np.flatnonzero(tt == j))] = cbar @ Wk
                score(f'own residual PCs K{K}, fit-half mean coefficient per tercile, by {lab}', f, Rh)
        log(f'mseterc2 fit half {f} done')
    R = pd.DataFrame(rows); R.to_csv(TAB / 'mse_terciles_corrections.csv', index=False)
    W_ = R.pivot_table(index='arm', columns=['part', 'fit_half'], values='scaled', sort=False)
    with pd.option_context('display.width', 250, 'display.max_rows', 100, 'display.max_colwidth', 90):
        print(W_.round(3).to_string())


# ================================================================ unit tests
def test(a):
    """Fast vs boot_score.Run.scaled on random resamples and on halves; fit_pick on a toy table; combined endpoints on
    a 20-knockdown slice of the per-knockdown arrays."""
    import boot_score as BS
    line = 'EvalA'; h = ref_for(line).h
    run = BS.Run(hdir('K0', line))
    F = Fast.from_wide(h, run.per_kd)
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(15):
        idx = rng.integers(0, h.T, h.T)
        want = run.scaled(pd.Index(h.targets)[idx])
        got, _ = F.scaled(idx[None, :])
        worst = max(worst, max(abs(got[k][0] - want[Fast.KEYS[k]]) for k in COMPS), abs(got['overall'][0] - want['overall']))
    sub = np.sort(rng.choice(h.T, 20, replace=False))
    want = run.scaled(pd.Index(h.targets)[sub]); got = F.point(sub)[0]
    worst = max(worst, abs(got['overall'] - want['overall']))
    log(f'test Fast vs Run.scaled: max |diff| {worst:.2e}'); assert worst < 1e-12, worst
    I = boot_index(5, 3); rng = np.random.default_rng(SEED)
    assert all((I[i] == rng.integers(0, 5, 5)).all() for i in range(3)), 'boot_index differs from boot_score sequence'
    tab = {'K4_all': {'half0': .2}, 'T1': {'half0': .3}, 'T2': {'half0': .3}, 'K0_all': {'half0': .1}}
    assert fit_pick(tab, 'half0', list(tab)) == 'T2', 'ties must go to the higher T'
    assert grid_pos(list(tab), 'T2') == '3/4 (interior)' and grid_pos(list(tab), 'K0_all') == '4/4 (edge)'
    b = np.array([np.nan, 0, 1, 5, 600.0]); A = assignments(b)
    assert A['T1'].tolist() == [False, False, True, True, True] and A['K4_all'].all() and not A['K0_all'].any()
    # Part B / D helpers on small slices
    x = np.array([3., 1, 2, 5, 4, 6]); assert tercile_of(x).tolist() == [1, 0, 0, 2, 1, 2]
    t = tiers(np.array([10., 1, 20, 30, 6, 7, 8]), np.array([1, 0, 1, 1, 1, 1, 1], bool)); assert t.tolist() == [1, -1, 2, 2, 0, 0, 1], t
    cg = np.array([[1.0, -4.0, 2.0, 3.0, 0.5, 1, 1]]); sg = np.array([[1, 1, 1, 0, 1, 0, 0]], bool)
    v, M = scaled_called(cg, sg, t, (0.35, 0.5, 1.0), 2.0)
    assert M.tolist() == [[True, False, True, False, True, False, False]] and np.allclose(v[0, [0, 2, 4]], [0.5, 2.0, 0.175])
    rng = np.random.default_rng(0); z = rng.normal(size=200); y = z + rng.normal(size=200); xx = z + rng.normal(size=200)
    pr, n = partial_spearman(y, xx, z); assert n == 200 and abs(pr) < 0.2, pr          # common cause removed
    pr2, _ = partial_spearman(y, y + 0.01 * rng.normal(size=200), z); assert pr2 > 0.95
    D = DLine('H1'); import scorer_arms as SA
    assert np.array_equal(D.mu((0.5, 0.5, 0.5), 'none'), SA.mu_from_lfc(D.lfc, D.m0)), 'H1 base mu differs from the generator input'
    m0n = np.load(C15 / 'config1p_K0/EvalA/input.npz')['m0']; rel = float(np.nanmax(np.abs(D.m0 - m0n) / np.maximum(m0n, 1e-12)))
    log(f'test: controls_baseline.csv m0 vs input.npz m0, max rel diff {rel:.1e} (CSV rounding; the 15c proxy used the CSV)'); assert rel < 1e-9
    log(f'test: H1 tiers {np.bincount(D.tier + 1).tolist()} (ungated, low, mid, high); called ungated pairs {int((D.sig & (D.tier < 0)[None]).sum())}')
    # Part E helpers on synthetic slices
    l_ = np.array([[1.0, 2.0, np.nan], [0.5, -1.0, 3.0]]); se_ = np.array([[1.0, 1.0, 1.0], [0.5, np.inf, 1.0]])
    x, m = e_entries(l_, se_, np.array([True, True, False]), 'lfc')
    assert np.allclose(x, [[0.5, 1.6, 0], [0.25, 0.0, 0]]) and m.tolist() == [[True, True, False], [True, True, False]], (x, m)
    rng = np.random.default_rng(3); Q, _ = np.linalg.qr(rng.standard_normal((40, 3))); Bt = Q.T
    Xs = (rng.standard_normal((6, 3)) @ Bt); gated = np.ones((6, 40), bool); gated[0, :5] = False
    trs = dict(X=np.where(gated, Xs, np.nan), gated=gated, DE=np.abs(Xs) > 0.2, Delta=Xs * 0.3, v=np.full(40, 2.0))
    sm = e_ceiling(Bt, 3, trs); assert abs(captured(sm['a'], np.array([0, 0, 1, 1, 2, 2]))['all'] - 1) < 1e-12
    assert abs(captured(sm['b'], np.array([0, 0, 1, 1, 2, 2]))['strong'] - 1) < 1e-12
    sm1 = e_ceiling(Bt, 1, trs); assert captured(sm1['a'], np.zeros(6, int))['all'] < 1
    # gate / E-reg helpers on toy slices
    cpm = {'X': np.array([1.0, 30, np.nan]), 'Y': np.array([10.0, 2, 4]), 'H1': np.array([100.0, 100, 100])}
    holds = {'X': {'a', 'b'}, 'Y': {'a'}, 'H1': {'a', 'b'}}
    sv = support_values('H1', ['a', 'b', 'c'], cpm, holds)
    assert np.allclose(sv[0], [5.5, 16, 4]) and np.allclose(sv[1][:2], [1, 30]) and np.isnan(sv[1][2]) and np.isnan(sv[2]).all(), sv
    mt = pd.DataFrame(dict(modality=['sc', 'sc', 'sc', 'drug', 'sc'], background=['P', 'P', 'Q', 'T', 'H1'],
                           table=['p1', 'p2', 'q1', 't1', 'h'], target=['a', 'a', 'a', 'a', 'a']))
    RSt = np.array([[1.0, 0], [3.0, 0], [10.0, 2], [100.0, 100], [1000.0, 0]], np.float32)
    sh, nb = ereg_khat(RSt, mt, ['a', 'z'], ('H1',))
    assert np.allclose(sh[0], [(2 + 10) / 2, 1]) and nb.tolist() == [2, 0] and (sh[1] == 0).all(), (sh, nb)
    sv_ = rng.normal(size=(30, 3)); assert np.allclose(fit_beta(2 * sv_, sv_, np.arange(30)), 2)
    D.lfc_base = D.lfc
    assert np.array_equal(ereg_mu(D, np.zeros_like(D.lfc), 'R2', 0.5), SA.mu_from_lfc(D.lfc, D.m0)), 'R2 with P = 0 must be the base'
    log('test: all unit tests pass')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['test', 'part0', 'partA', 'partB', 'partC', 'partD1', 'partD', 'partDsum', 'partDgate', 'partE', 'partEreg', 'partF', 'mseterc', 'mseterc2'])
    p.add_argument('--line')
    p.add_argument('--worker', type=int, default=-1)
    p.add_argument('--spec', default='{}', help='partF: JSON combination spec')
    a = p.parse_args()
    globals()[a.cmd](a)

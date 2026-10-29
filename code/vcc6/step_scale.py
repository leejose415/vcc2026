"""Step 1, the magnitude of C+G: two arms, each scored against C+G on the three-line screen (effect space).

  E  (option a) per-knockdown scale by the knocked-down gene's expression. On held-out SOURCE lines, every knockdown's
     least-squares scale of true on predicted lfc, s_k = sum(x y) / sum(x^2) over its tested genes, is regressed on
     e = log1p(control mean CPM of the knocked-down gene in that line) with weighted isotonic regression (increasing,
     weights sum(x^2)). The target's lfc = s(e_target) x G, e_target from the target line's own controls (its baseline).
  M  the magnitude curve. On the same held-out source lines, isotonic regression (increasing) of true lfc on predicted
     lfc, gene by gene, pooled; the target's lfc = f(G) where a source measured the gene (0 elsewhere, as G).
Calibration lines: the single-cell CRISPRi lines with tables (K562, RPE1, H1, HepG2, Jurkat, KOLF2.1J, HCT116,
HEK293T), minus the target's own line. For a target T and a calibration line L, L is predicted with C+G from every
table except T's and L's backgrounds (T never informs its own curve). Up to N_CAL random knockdowns per L (seed 2026).
Calls follow the lfc (stage1.predictive_calls). R^2 and beta (calibration_stats) are reported in sample (the pooled
source-line pairs each curve was fitted on) and out of sample (the target line): a drop flags overfitting. Keep: NMAE and MSE better than C+G on >= 2 of the 3 lines, no tercile
worse (bootstrap CI of the paired median difference entirely on the worse side).

--set final (Jose 2026-10-27): the final-target screen (screen.FINAL, the 2026-10-27_final_baseline C+G predictions).
  beta  the reference: C+G x one beta = sum(x y) / sum(x^2) over the pooled held-out source-line pairs (no intercept).
        E and M must beat it; any extra parameter has to earn its place against this one.
  Every variant keeps G's call set (calls from the C+G lfc), so only the magnitude differs. Scorer-space pb_nmse
  (screen.extra_metrics) per variant. Keep: K562, HCT116, HEK293T count; H1 (28 targets) is reported only.
  One json per line: results/stage1/step1_scale/final/<line>.json.

  ~/venv/bin/python code/vcc6/step_scale.py [--lines RPE1 H1 K562] [--n-cal 200]
  ~/venv/bin/python code/vcc6/step_scale.py --set final [--lines K562 HCT116 HEK293T H1]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
import stage1 as S                                                                         # noqa: E402
import screen as SC                                                                        # noqa: E402

OUT = V / 'results/stage1/step1_scale'
BASE = V / 'results/stage1/screen/2026-10-27_baseline'
BASE_FINAL = V / 'results/stage1/screen/2026-10-27_final_baseline'
CAL = ['replogle_k562_gw', 'replogle_rpe1', 'arc_h1', 'nadig_hepg2', 'nadig_jurkat', 'kolf_kolf21j',
       'orion_hct116', 'orion_hek293t']
MAX_PAIRS = 3_000_000


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


def own_expression(basal_cpm, genes, knockdowns):
    col = {str(g): j for j, g in enumerate(genes)}
    return np.array([np.log1p(basal_cpm[col[k]]) if k in col and np.isfinite(basal_cpm[col[k]]) else 0.0
                     for k in knockdowns])


def calibration_pairs(sources, by_stem, target_bg, n_cal, out):
    """For every calibration line L (not the target's): C+G predictions of up to n_cal of L's knockdowns from every
    table except T's and L's backgrounds, against L's table. Returns per-knockdown records and pooled gene pairs."""
    pool = [s for s in sources if S._line_key(s.background) != target_bg]
    recs, xs, ys, es = [], [], [], []
    for L in CAL:
        src = next((s for s in sources if s.name == L), None)
        if src is None or S._line_key(src.background) == target_bg:
            continue
        line = S.line_for_target(pool, by_stem, L)
        kds = sorted(src.rows_of)
        rng = np.random.default_rng(2026)
        pick = sorted(rng.choice(np.array(kds, dtype=object), min(n_cal, len(kds)), replace=False).tolist())
        kf = out / f'cal_{L}.txt'; kf.write_text('\n'.join(pick) + '\n')
        d = out / f'cal_{L}'
        S.predict_line(line, str(d), knockdowns_file=str(kf), arm='C', magnitude='gaussian')
        z = np.load(d / 'predictions.npz', allow_pickle=True)
        truth = S.truth_signature(by_stem[L], pick, line['genes'])
        e = own_expression(line['basal_cpm'], line['genes'], pick)
        gate = line['tested']
        col = {str(g): j for j, g in enumerate(line['genes'])}
        for i, k in enumerate(pick):
            if k not in truth: continue
            x = z['lfc'][i].astype(float); y = truth[k][0]
            ok = gate & np.isfinite(y) & np.isfinite(x) & (z['n_sources'][i] > 0)
            if k in col: ok[col[k]] = False
            sxx = float((x[ok] ** 2).sum())
            if sxx <= 0: continue
            recs.append(dict(line=L, knockdown=k, expr=float(e[i]), scale=float((x[ok] * y[ok]).sum() / sxx), weight=sxx))
            xs.append(x[ok]); ys.append(y[ok]); es.append(np.full(int(ok.sum()), e[i]))
        log(f'  calibration {L}: {len(pick)} knockdowns')
    return pd.DataFrame(recs), np.concatenate(xs), np.concatenate(ys), np.concatenate(es)


def fit_curves(recs, x, y):
    iso_e = IsotonicRegression(increasing=True, out_of_bounds='clip', y_min=0.0)
    iso_e.fit(recs.expr.to_numpy(), recs.scale.to_numpy(), sample_weight=recs.weight.to_numpy())
    if x.size > MAX_PAIRS:
        take = np.random.default_rng(2026).choice(x.size, MAX_PAIRS, replace=False); x, y = x[take], y[take]
    y = np.clip(y, -S.LFC_CAP, S.LFC_CAP)                          # the epsilon-floor lfc of the truth, as the pool drops
    iso_m = IsotonicRegression(increasing=True, out_of_bounds='clip')
    iso_m.fit(x, y)
    return iso_e, iso_m


def score_variant(name, lfc, z, truth, gate, want, genes, calls_lfc=None):
    """score_against_truth for one variant; calls follow its lfc, or calls_lfc's (G's call set) when given."""
    rows = []
    col = {str(g): j for j, g in enumerate(genes)}
    for i, k in enumerate(want):
        if k not in truth: continue
        tested = gate.copy()
        if k in col: tested[col[k]] = False
        sig, _ = S.predictive_calls(lfc[i] if calls_lfc is None else calls_lfc[i], z['se_measurement'][i],
                                    tested & (z['n_sources'][i] > 0))
        pred = dict(lfc=lfc[i].astype(np.float32), lfc_unshrunk=z['lfc_unshrunk'][i], de_sig=sig,
                    v_post=np.zeros(len(genes), np.float32), se_measurement=z['se_measurement'][i])
        y, e, tsig = truth[k]
        rows.append(dict(knockdown=k, **S.score_against_truth(pred, y, e, tsig, gate)))
    fit = S.calibration_stats(lfc, {i: truth[k][0] for i, k in enumerate(want) if k in truth}, gate)
    return rows, fit


def main(lines, n_cal, line_set='v1'):
    final = line_set == 'final'
    root, base = (OUT / 'final', BASE_FINAL) if final else (OUT, BASE)
    root.mkdir(parents=True, exist_ok=True)
    tables = S.discover_tables(SC.DIRS); by_stem = {p.stem: p for p in tables}
    sources = S.load_tables(tables)
    res = {}
    for name in lines:
        if final:
            target, lst = SC.FINAL[name], V / f'results/target_lists/final_{name}.txt'
        else:
            target, lst = SC.LINES[name]
        want = [t for t in lst.read_text().split() if t]
        out = root / name; out.mkdir(parents=True, exist_ok=True)
        tsrc = next(s for s in sources if s.name == target); tbg = S._line_key(tsrc.background)
        log(f'{name}: calibration on held-out source lines (target background {tbg} excluded)')
        recs, x, y, ex = calibration_pairs(sources, by_stem, tbg, n_cal, out)
        recs.to_csv(out / 'calibration_knockdowns.csv', index=False)
        iso_e, iso_m = fit_curves(recs, x, y)
        beta = float((x * y).sum() / (x * x).sum())
        def pooled(pred, truth):
            return S.calibration_stats(pred[None, :], {0: truth}, np.ones(pred.size, bool))
        in_sample = dict(G=pooled(x, y), beta=pooled(beta * x, y), E=pooled(iso_e.predict(ex) * x, y),
                         M=pooled(iso_m.predict(x), y))
        grid = np.linspace(0, 8, 33)
        xq = np.quantile(x, [0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999])
        curves = dict(E=dict(expr_log1p_cpm=grid.tolist(), scale=iso_e.predict(grid).tolist()),
                      M=dict(pred_quantiles=xq.tolist(), curve=iso_m.predict(xq).tolist()),
                      calibration_lines=sorted(recs.line.unique().tolist()), calibration_knockdowns=int(len(recs)),
                      pooled_pairs=int(x.size), pooled_slope_true_on_pred=float((x * y).sum() / (x * x).sum()))
        # the target: the baseline screen's own C+G predictions (identical baseline), truth from its table
        z = np.load(base / name / 'C_gaussian' / 'predictions.npz', allow_pickle=True)
        assert [str(k) for k in z['knockdowns']] == want
        genes = [str(g) for g in z['genes']]
        basal, gate = S.load_baseline(by_stem[target].with_name(target + '_baseline.csv'), np.array(genes, dtype=object))
        truth = S.truth_signature(by_stem[target], want, np.array(genes, dtype=object))
        G = z['lfc'].astype(float)
        e_t = own_expression(basal, genes, want)
        meas = z['n_sources'] > 0
        variants = {'C+G': G, 'E': iso_e.predict(e_t)[:, None] * G, 'M': np.where(meas, iso_m.predict(G.ravel()).reshape(G.shape), 0.0)}
        if final:
            variants = {'C+G': G, 'beta': beta * G, 'E': variants['E'], 'M': variants['M']}
        terc = pd.read_csv(base / f'{name}_terciles.csv', index_col=0)['tercile']
        rows, summ = {}, {}
        for v, L in variants.items():
            rows[v], fit = score_variant(v, L, z, truth, gate, want, genes, calls_lfc=G if final else None)
            if final:
                extra = {r['knockdown']: r for r in SC.extra_metrics(base / name / 'C_gaussian' / 'predictions.npz',
                                                                     by_stem[target], want, genes, lfc_all=L)}
                for r in rows[v]:
                    r.update({k: w for k, w in extra.get(r['knockdown'], {}).items() if k != 'knockdown'})
            summ[v] = SC.summarise(rows[v], dict(shrunk=fit), terc)
        out_of_sample = {v: summ[v]['fit']['shrunk'] for v in variants}
        pairs = [('E', 'C+G'), ('M', 'C+G')]
        if final:
            pairs = [('beta', 'C+G'), ('E', 'beta'), ('M', 'beta'), ('E', 'C+G'), ('M', 'C+G')]
        res[name] = dict(curves=dict(curves, beta=beta), fit_in_sample_source_lines=in_sample,
                         fit_out_of_sample_target=out_of_sample,
                         scale_applied=dict(median=float(np.median(iso_e.predict(e_t))),
                                                           q10_q90=np.quantile(iso_e.predict(e_t), [0.1, 0.9]).tolist()),
                         summary=summ, compare={f'{c} vs {b}': SC.compare(rows[b], rows[c], terc) for c, b in pairs})
        if final:
            (root / f'{name}.json').write_text(json.dumps(res[name], indent=1, default=str) + '\n')
        else:
            (OUT / 'step1_scale.json').write_text(json.dumps(res, indent=1, default=str) + '\n')
        log(f'{name}: done (beta {beta:.3f})')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--set', default='v1', choices=['v1', 'final'])
    p.add_argument('--lines', nargs='+', default=None)
    p.add_argument('--n-cal', type=int, default=200)
    a = p.parse_args()
    main(a.lines or list(SC.FINAL if a.set == 'final' else SC.LINES), a.n_cal, a.set)

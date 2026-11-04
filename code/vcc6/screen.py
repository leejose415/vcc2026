"""The three-line effect-space screen (keep rule, Jose 2026-10-27): a change is kept if it improves effect-space NMAE,
MSE and sign on at least 2 of the 3 lines, and no tercile of the line it was tuned on gets worse by more than the noise.
Fidelity and call count never count toward a keep.

  lines   RPE1  replogle_rpe1     results/target_lists/rpe1_dev_300.txt
          H1    arc_h1            results/target_lists/h1_targets_300.txt   (the arc_h1 effect table; effect space only)
          K562  replogle_k562_gw  results/target_lists/k562_heldout_evalB_100.txt
  per line  stage1.line_for_target (every table of the line's background held out) -> predict_line(arm, magnitude);
            truth = the line's own effect table. Metrics per knockdown (score_against_truth), summarised as the
            median over knockdowns with a bootstrap 95% CI (2,000 knockdown resamples, seed 2026); r and both slopes
            pooled (calibration_stats); terciles by the table's pseudobulk distance from control
            (sum of squared log1p_bulk_50k_matched - control over the measured genes, own gene excluded).
  compare   paired bootstrap of the median difference (candidate - baseline) per line and tercile.

Line set "final" (Jose 2026-10-27, after the D/E/F target list): the same sources vote as will vote on D/E/F.
  K562     replogle_k562_gw          its final targets (267)   every K562 table held out
  H1       arc_h1                    its final targets (28)    reported; too few to count toward a keep (proposed)
  HCT116   orion_hct116              its final targets (278)   both Orion HCT116 tables held out (single-cell since
  HEK293T  orion_hek293t             its final targets (271)   the 10-27 re-stream: QC'd cells, DE labels)
  The lists are built at run time: data_final/extracted/pert_counts.csv x the line's table, written to
  results/target_lists/final_<line>.txt.
Every line also gets, per knockdown: sign_z3 (sign accuracy on the genes whose true |lfc / se| > 3; the sign metric
where a table has no DE labels) and scorer-space errors: the predicted lfc turned into the scorer's pseudobulk with the
generator's own inversion on the line's control mean (log1p(5e4 x m1 / sum m1)), against the table's pseudobulk
(log1p_bulk_50k_matched, or _full where matched is absent), own gene excluded. Line level: pb_nmse = sum of squared
errors / sum of the true distances from control (no noise correction), for the prediction and for the mean-response
baseline (each knockdown predicted as the mean true pseudobulk over the scored list: the analogue of the scorer's 0).

  ~/venv/bin/python code/vcc6/screen.py --tag NAME --variant C:gaussian [--variant C:slab ...] [--set final]
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
import stage1 as S                                                                         # noqa: E402

DIRS = [str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')]
LINES = {'RPE1': ('replogle_rpe1', V / 'results/target_lists/rpe1_dev_300.txt'),
         'H1': ('arc_h1', V / 'results/target_lists/h1_targets_300.txt'),
         'K562': ('replogle_k562_gw', V / 'results/target_lists/k562_heldout_evalB_100.txt')}
# Orion: since the 2026-10-27 re-stream the final targets are in the single-cell tables (QC'd cells, DE labels)
FINAL = {'K562': 'replogle_k562_gw', 'H1': 'arc_h1', 'HCT116': 'orion_hct116', 'HEK293T': 'orion_hek293t'}
FINAL_TARGETS = V / 'data_final/extracted/pert_counts.csv'
KEEP = ('nmae', 'mse', 'sign_acc_true_de', 'sign_z3')  # lower, lower, higher, higher
HIGHER = {'sign_acc_true_de', 'sign_z3', 'auroc'}
OUT = V / 'results/stage1/screen'
N_BOOT, SEED = 2000, 2026


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


def bulk_layer(a):
    """log1p_bulk_50k_matched, or _full where the table has no matched profile (the Orion pseudobulk tables)."""
    m = np.asarray(a.layers['log1p_bulk_50k_matched'][0], float)
    return 'log1p_bulk_50k_matched' if np.isfinite(m).any() else 'log1p_bulk_50k_full'


def final_list(name, target, by_stem):
    """The line's final targets: data_final pert_counts x its table, written to results/target_lists/final_<line>.txt."""
    import h5py
    from anndata.io import read_elem
    fin = sorted(set(pd.read_csv(FINAL_TARGETS)['target_gene'].astype(str)))
    with h5py.File(by_stem[target], 'r') as fh:
        have = set(map(str, read_elem(fh['obs'])['target_gene']))
    path = V / f'results/target_lists/final_{name}.txt'
    path.write_text('\n'.join(k for k in fin if k in have) + '\n')
    return path


def extra_metrics(pred_npz, table_path, want, genes, lfc_all=None):
    """Per knockdown: sign_z3, and scorer-space squared errors of the prediction and of the mean-response baseline.
    lfc_all: (K, G) lfc to use instead of the npz's own (step_scale's rescaled variants)."""
    import anndata as ad
    sys.path.insert(0, str(V / 'code/cellgen'))
    import cell_generator as C
    z = np.load(pred_npz, allow_pickle=True)
    a = ad.read_h5ad(table_path, backed='r')
    base = pd.read_csv(Path(table_path).with_name(Path(table_path).stem + '_baseline.csv'), index_col=0).reindex(genes)
    layer = bulk_layer(a)
    m0 = base['control_mean_cpm'].to_numpy(float)
    ctrl = base[layer].to_numpy(float)
    row = {s: i for i, s in enumerate(a.obs['target_gene'].astype(str))}
    col = {g: j for j, g in enumerate(genes)}
    meas = np.isfinite(m0) & np.isfinite(ctrl)
    truth_pb = np.vstack([np.asarray(a.layers[layer][row[k]], float) for k in want])
    mean_pb = np.nanmean(truth_pb, axis=0)
    out = []
    for i, k in enumerate(want):
        lfc = np.asarray(z['lfc'][i] if lfc_all is None else lfc_all[i], float)
        m1 = np.where(np.isfinite(m0), C.target_means(lfc, m0), 0.0)
        pb = np.log1p(5e4 * m1 / m1.sum())
        y = np.asarray(a.layers['lfc'][row[k]], float); se = np.asarray(a.layers['se'][row[k]], float)
        zt = np.divide(y, se, out=np.zeros_like(y), where=np.isfinite(se) & (se > 0))
        strong = np.isfinite(zt) & (np.abs(zt) > 3) & np.isfinite(y)
        if k in col: strong[col[k]] = False
        ok = meas & np.isfinite(truth_pb[i])
        if k in col: ok[col[k]] = False
        t = truth_pb[i]
        out.append(dict(knockdown=k,
                        sign_z3=float((np.sign(lfc[strong]) == np.sign(y[strong])).mean()) if strong.any() else None,
                        pb_err=float(((pb[ok] - t[ok]) ** 2).sum()), pb_den=float(((t[ok] - ctrl[ok]) ** 2).sum()),
                        pb_err_meanresp=float(((mean_pb[ok] - t[ok]) ** 2).sum())))
    return out


def pb_summary(rows, terc):
    """pb_nmse of the prediction and of the mean-response baseline, overall and per tercile, with a bootstrap CI of
    their difference (prediction - baseline; negative = the prediction beats the scorer's 0)."""
    df = pd.DataFrame(rows).set_index('knockdown')
    rng = np.random.default_rng(SEED)
    out = {}
    for name, sub in [('all', df)] + [(t, df[terc.reindex(df.index) == t]) for t in ('strong', 'mid', 'weak')]:
        e, b, d = sub.pb_err.to_numpy(), sub.pb_err_meanresp.to_numpy(), sub.pb_den.to_numpy()
        if d.sum() <= 0: continue
        ii = rng.integers(0, len(sub), size=(N_BOOT, len(sub)))
        diff = (e[ii].sum(1) - b[ii].sum(1)) / d[ii].sum(1)
        out[name] = dict(pb_nmse=float(e.sum() / d.sum()), pb_nmse_meanresp=float(b.sum() / d.sum()),
                         diff=float((e.sum() - b.sum()) / d.sum()), diff_ci=[float(np.quantile(diff, 0.025)), float(np.quantile(diff, 0.975))])
    return out


def terciles(path, want, genes):
    """strong / mid / weak by the table's pseudobulk distance from its control."""
    import anndata as ad
    a = ad.read_h5ad(path, backed='r')
    base = pd.read_csv(Path(path).with_name(Path(path).stem + '_baseline.csv'), index_col=0)
    layer = bulk_layer(a)
    ctrl = base.reindex(list(map(str, genes)))[layer].to_numpy(float)
    row = {s: i for i, s in enumerate(a.obs['target_gene'].astype(str))}
    col = {g: j for j, g in enumerate(map(str, genes))}
    d = {}
    for k in want:
        b = np.asarray(a.layers[layer][row[k]], float)
        ok = np.isfinite(b) & np.isfinite(ctrl)
        if k in col: ok[col[k]] = False
        d[k] = float(((b[ok] - ctrl[ok]) ** 2).sum())
    s = pd.Series(d).rank(method='first', ascending=False)
    n = len(s)
    return s.map(lambda r: 'strong' if r <= n / 3 else ('mid' if r <= 2 * n / 3 else 'weak')), pd.Series(d)


def boot_median(x, rng):
    x = np.asarray([v for v in x if v is not None and np.isfinite(v)], float)
    if x.size == 0: return None, None, None
    idx = rng.integers(0, x.size, size=(N_BOOT, x.size))
    m = np.median(x[idx], axis=1)
    return float(np.median(x)), float(np.quantile(m, 0.025)), float(np.quantile(m, 0.975))


def summarise(rows, fit, terc):
    df = pd.DataFrame(rows).set_index('knockdown')
    rng = np.random.default_rng(SEED)
    out = {}
    for name, sub in [('all', df)] + [(t, df[terc.reindex(df.index) == t]) for t in ('strong', 'mid', 'weak')]:
        rec = dict(knockdowns=len(sub))
        for m in ('nmae', 'mse', 'sign_acc_true_de', 'sign_z3', 'auroc', 'nmae_zero', 'mse_zero', 'n_called', 'n_true'):
            if m in sub:
                med, lo, hi = boot_median(sub[m].tolist(), rng)
                rec[m] = med; rec[f'{m}_ci'] = [lo, hi]
        out[name] = rec
    out['fit'] = fit
    if rows and 'pb_err' in rows[0]:
        out['scorer_space'] = pb_summary(rows, terc)
    es = S.summarize_scores(rows)
    out['call_count_pearson'] = es['call_count_pearson']
    return out


NON_CELL = ('pseudobulk_table', 'replicate_table')     # Orion pseudobulk x2, Marson CD4 x3: no cells, no DE labels


def run_variant(line, name, want, out_dir, arm, magnitude, list_path=None, cells_only=False):
    lst = list_path or LINES[name][1]
    if cells_only:                                          # step 5a: the five non-cell tables never vote
        line = dict(line, others=[s for s in line['others'] if s.modality not in NON_CELL])
    S.predict_line(line, str(out_dir), knockdowns_file=str(lst), arm=arm, magnitude=magnitude)
    kd = json.loads((out_dir / 'per_knockdown_scores.json').read_text())
    summ = json.loads((out_dir / 'summary.json').read_text())
    extra = {r['knockdown']: r for r in extra_metrics(out_dir / 'predictions.npz', line['by_stem'][line['target']], want,
                                                     [str(g) for g in line['genes']])}
    for r in kd:
        r.update({k: v for k, v in extra.get(r['knockdown'], {}).items() if k != 'knockdown'})
    return kd, summ['fit'], summ['prior']


def compare(base_rows, cand_rows, terc):
    """Paired bootstrap of median(candidate) - median(baseline) for the keep metrics, per tercile."""
    b = pd.DataFrame(base_rows).set_index('knockdown'); c = pd.DataFrame(cand_rows).set_index('knockdown')
    common = b.index.intersection(c.index)
    rng = np.random.default_rng(SEED)
    out = {}
    for name, idx in [('all', common)] + [(t, common[terc.reindex(common).to_numpy() == t]) for t in ('strong', 'mid', 'weak')]:
        rec = {}
        for m in KEEP + ('auroc',):
            if m not in b.columns or m not in c.columns:
                continue
            x, y = b.loc[idx, m].to_numpy(float), c.loc[idx, m].to_numpy(float)
            ok = np.isfinite(x) & np.isfinite(y); x, y = x[ok], y[ok]
            if x.size == 0: continue
            ii = rng.integers(0, x.size, size=(N_BOOT, x.size))
            d = np.median(y[ii], axis=1) - np.median(x[ii], axis=1)
            better = (d > 0) if m in HIGHER else (d < 0)
            rec[m] = dict(diff=float(np.median(y) - np.median(x)), ci=[float(np.quantile(d, 0.025)), float(np.quantile(d, 0.975))],
                          p_better=float(better.mean()))
        out[name] = rec
    return out


def main(tag, variants, lines, line_set='v1'):
    out = OUT / tag; out.mkdir(parents=True, exist_ok=True)
    tables = S.discover_tables(DIRS); by_stem = {p.stem: p for p in tables}
    sources = S.load_tables(tables)
    res = {}
    for name in lines:
        if line_set == 'final':
            target = FINAL[name]; lst = final_list(name, target, by_stem)
        else:
            target, lst = LINES[name]
        want = [t for t in lst.read_text().split() if t]
        log(f'{name}: {target}, {len(want)} knockdowns ({lst})')
        line = S.line_for_target(sources, by_stem, target)
        terc, dist = terciles(by_stem[target], want, line['genes'])
        pd.DataFrame(dict(tercile=terc, pseudobulk_distance=dist)).to_csv(out / f'{name}_terciles.csv')
        res[name] = {}
        rows_of = {}
        for v in variants:
            arm, mag, *opt = v.split(':')
            log(f'{name}: {v}')
            kd, fit, prior = run_variant(line, name, want, out / name / v.replace(':', '_'), arm, mag, lst,
                                         cells_only='cells' in opt)
            rows_of[v] = kd
            res[name][v] = dict(summary=summarise(kd, fit, terc), prior=prior)
        if len(variants) > 1:
            res[name]['compare'] = {f'{v} vs {variants[0]}': compare(rows_of[variants[0]], rows_of[v], terc) for v in variants[1:]}
        (out / 'screen.json').write_text(json.dumps(res, indent=1, default=str) + '\n')
        del line
    log(f'wrote {out / "screen.json"}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--tag', required=True)
    p.add_argument('--variant', action='append', required=True, help='ARM:MAGNITUDE[:cells], the first is the baseline (cells: no non-cell tables)')
    p.add_argument('--set', default='v1', choices=['v1', 'final'], help='line set (final: the D/E/F target list)')
    p.add_argument('--lines', nargs='+', default=None)
    a = p.parse_args()
    main(a.tag, a.variant, a.lines or list(FINAL if a.set == 'final' else LINES), a.set)

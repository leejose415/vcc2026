"""Screen H1 effect arrays without generating cells, using the generator's exact calibration.

Each target contributes equally. Differences and paired 95% bootstrap CIs are
arm minus re (alpha=1, min_effect=0); lower MSE and higher other metrics are better.
Only finite measured effects on H1 ctrl_cpm > 5 genes are scored; own genes are
excluded. Constant-vector Pearson is defined as zero.

  top50_overlap       share of the 50 strongest NONZERO predictions that are among the 50 strongest measured
                      genes (denominator 50): an all-zero prediction scores 0, never a chance fill by gene order
  direction_accuracy  sign agreement on significant genes (|lfc| / SE > 3) that the arm actually predicts (nonzero)
  direction_coverage  share of the significant genes the arm predicts at all; accuracy and coverage are kept apart
  MSE is reported, not used for selection: H1's measured effects are mostly noise, so squared error rewards
  shrinking toward zero. Pearson, top-50 overlap and direction accuracy do not depend on alpha, so the screen
  chooses the array and minimum size; alpha is left to the official benchmark.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from types import SimpleNamespace

import anndata as ad
import h5py
import numpy as np
import pandas as pd
try:
    from anndata.io import read_elem
except ImportError:  # older anndata
    from anndata.experimental import read_elem

from vcc2_common import expression_gate, read_names, save_json
from vcc3_generate import Cfg, calibrate, load_deltas

ALPHAS = (0.5, 0.75, 1.0, 1.5, 2.0)
MINIMUMS = (0.0, 0.1, 0.2)
METRICS = ('mse', 'pearson', 'top50_overlap', 'direction_accuracy', 'direction_coverage')


def unique_names(values, label):
    names = pd.Index(np.asarray(values, dtype=str))
    if not names.is_unique: raise ValueError(f'{label}: duplicate names')
    return names


def required_indices(source, requested, label):
    idx = source.get_indexer(requested)
    if (idx < 0).any():
        raise ValueError(f'{label}: missing {np.asarray(requested)[idx < 0].tolist()}')
    return idx


def load_inputs(predictions, truth, filters, targets_csv):
    p = Path(predictions)
    report = json.loads((p / 'prediction_report.json').read_text())
    if report.get('heldout_background') != 'H1' or 'H1' in report.get('source_backgrounds', []):
        raise ValueError('screen requires predictions fitted with H1 held out')
    genes = unique_names(np.load(p / 'genes.npy'), 'prediction genes')
    pred_targets = unique_names(np.load(p / 'targets.npy'), 'prediction targets')
    requested = unique_names(read_names(targets_csv, 'target_gene'), 'validation targets')
    required_indices(pred_targets, requested, 'prediction targets')
    panel = ad.read_h5ad(truth)
    ti = required_indices(unique_names(panel.obs_names, 'truth targets'), requested, 'truth targets')
    gi = required_indices(unique_names(panel.var_names, 'truth genes'), genes, 'truth genes')
    kept = panel.obs.iloc[ti]['status'].astype(str).to_numpy() == 'kept'
    targets = requested[kept].to_numpy(dtype=str); ti = ti[kept]
    if not len(targets): raise ValueError('no kept validation targets')
    measured = np.asarray(panel.layers['lfc'])[np.ix_(ti, gi)].astype(np.float64)
    se = np.asarray(panel.layers['se'])[np.ix_(ti, gi)].astype(np.float64)
    with h5py.File(filters, 'r') as f: var = read_elem(f['var'])
    unique_names(var.index, 'filter genes')
    cpm = var['ctrl_cpm'].reindex(genes).to_numpy(float)
    gate = expression_gate(cpm)
    eligible = np.isfinite(measured) & gate[None, :] & (targets[:, None] != genes.to_numpy()[None, :])
    if (eligible.sum(1) == 0).any(): raise ValueError('validation target has no eligible measured genes')
    significant = eligible & np.isfinite(se) & (se >= 0) & (np.abs(measured) > 3 * se)
    info = dict(requested_targets=len(requested), n_targets=len(targets),
                excluded_not_kept=requested[~kept].tolist(), genes=len(genes), gated_genes=int(gate.sum()),
                genes_absent_from_filter=int((var.index.get_indexer(genes) < 0).sum()),
                n_significant_targets=int(significant.any(1).sum()),
                measured_genes_min=int(eligible.sum(1).min()), measured_genes_max=int(eligible.sum(1).max()))
    return genes.to_numpy(dtype=str), targets, measured, se, eligible, significant, gate, info


def target_metrics(predicted, measured, eligible, significant):
    out = np.full((len(measured), len(METRICS)), np.nan)
    for i in range(len(measured)):
        x, y = predicted[i, eligible[i]].astype(float), measured[i, eligible[i]]
        out[i, 0] = np.mean((x - y) ** 2)
        xc, yc = x - x.mean(), y - y.mean()
        denom = np.linalg.norm(xc) * np.linalg.norm(yc)
        out[i, 1] = np.dot(xc, yc) / denom if denom > 0 else 0.0
        n = min(50, len(x))
        nonzero = np.flatnonzero(x != 0)
        top_x = nonzero[np.argsort(-np.abs(x[nonzero]), kind='stable')[:n]]       # only genes the arm predicts
        top_y = np.argsort(-np.abs(y), kind='stable')[:n]
        out[i, 2] = np.isin(top_x, top_y).sum() / n
        sig = significant[i]
        if sig.any():
            called = sig & (predicted[i] != 0)
            out[i, 4] = called.sum() / sig.sum()
            if called.any(): out[i, 3] = np.mean(np.sign(predicted[i, called]) == np.sign(measured[i, called]))
    return out


def paired_summary(values, reference, draws):
    diff = values - reference
    with np.errstate(invalid='ignore', divide='ignore'):
        sampled = diff[draws]
        valid = np.isfinite(sampled)
        boot = np.divide(np.where(valid, sampled, 0).sum(1), valid.sum(1),
                         out=np.full((len(draws), len(METRICS)), np.nan), where=valid.sum(1) > 0)
    out = {}
    for j, metric in enumerate(METRICS):
        finite = np.isfinite(values[:, j]); b = boot[np.isfinite(boot[:, j]), j]
        out[metric] = float(values[finite, j].mean()) if finite.any() else np.nan
        out[metric + '_n_targets'] = int(finite.sum())
        out[metric + '_diff'] = float(diff[finite, j].mean()) if finite.any() else np.nan
        lo, hi = np.percentile(b, [2.5, 97.5]) if len(b) else (np.nan, np.nan)
        out[metric + '_diff_ci_low'], out[metric + '_diff_ci_high'] = float(lo), float(hi)
    return out


def shortlist(frame):
    """Best array and minimum size per family by Pearson at alpha = 1 (Pearson does not depend on alpha), ties by
    top-50 overlap; plus any arm whose top-50 overlap is clearly above Model 2's."""
    chosen = []
    at_one = frame[frame.alpha == 1.0]
    for family in ('re', 'programs', 'corrected'):
        candidates = at_one[at_one.family == family].sort_values(['pearson', 'top50_overlap', 'arm'], ascending=[False, False, True])
        if candidates.empty: raise ValueError(f'no candidates for {family}')
        row = candidates.iloc[0]
        chosen.append(dict(arm=row.arm, delta_array=row.delta_array, alpha=float(row.alpha),
                           min_effect=float(row.min_effect), reason=f'best {family} Pearson'))
    winners = at_one[(at_one.family != 'zero') & (at_one.top50_overlap_diff_ci_low > 0)]
    if not winners.empty:
        row = winners.sort_values(['top50_overlap', 'mse', 'arm'], ascending=[False, True, True]).iloc[0]
        if row.arm not in {x['arm'] for x in chosen}:
            chosen.append(dict(arm=row.arm, delta_array=row.delta_array, alpha=float(row.alpha),
                               min_effect=float(row.min_effect), reason='best top-50 overlap with CI above re'))
    return chosen


def screen(predictions, truth, filters, targets_csv, output, seed=2026, resamples=1000):
    if resamples < 1: raise ValueError('resamples must be positive')
    genes, targets, measured, se, eligible, significant, gate, info = load_inputs(predictions, truth, filters, targets_csv)
    p = Path(predictions); out = Path(output); out.mkdir(parents=True, exist_ok=True)
    models = sorted(f.stem for f in p.glob('*.npy') if re.fullmatch(r're|zero|(?:programs|corrected)_k\d+', f.stem))
    if not {'re', 'zero'}.issubset(models): raise ValueError('screen requires re.npy and zero.npy')
    line = SimpleNamespace(genes=genes, gate=gate)
    values, settings = {}, []
    expected_shape = (len(np.load(p / 'targets.npy')), len(genes))
    for model in models:
        raw = np.load(p / f'{model}.npy', mmap_mode='r')
        if raw.shape != expected_shape: raise ValueError(f'{model}: wrong array shape {raw.shape}')
        if model == 'zero' and not np.all(raw == 0): raise ValueError('zero.npy is not all zero')
        D = load_deltas(p, targets, genes, model)  # generator convention: all non-finite predictions -> 0
        for alpha in ALPHAS:
            for minimum in MINIMUMS:
                arm = f'{model}_a{alpha:g}_m{minimum:g}'
                cfg = Cfg(SimpleNamespace(alpha=alpha, min_effect=minimum, z=0))
                transformed = calibrate(D, line, targets, cfg)
                values[arm] = target_metrics(transformed, measured, eligible, significant)
                settings.append(dict(arm=arm, delta_array=model, family=model.split('_k')[0], alpha=alpha, min_effect=minimum))
        print(f'screen: {model}, {len(targets)} targets, 15 settings', flush=True)
    draws = np.random.default_rng(seed).integers(0, len(targets), size=(resamples, len(targets)))
    np.save(out / 'bootstrap_indices.npy', draws)
    pd.DataFrame(dict(target_gene=targets, n_genes=eligible.sum(1), n_significant_genes=significant.sum(1))).to_csv(out / 'targets.csv', index=False)
    reference = values['re_a1_m0']
    rows, per_target = [], []
    zero = values['zero_a1_m0']
    for setting in settings:
        metrics = values[setting['arm']]
        rows.append(dict(**setting, **paired_summary(metrics, reference, draws)))
        t = pd.DataFrame(metrics, columns=METRICS)
        t.insert(0, 'target_gene', targets); t.insert(0, 'arm', setting['arm']); per_target.append(t)
    frame = pd.DataFrame(rows).sort_values(['mse', 'arm']).reset_index(drop=True)
    frame.to_csv(out / 'screen.csv', index=False)
    pd.concat(per_target, ignore_index=True).to_csv(out / 'per_target.csv', index=False)
    finalists = shortlist(frame)
    # Alignment gate: zero must be worst on the signal metrics (Pearson and top-50 overlap). Every shortlisted arm
    # must beat it with a paired 95% CI above zero on both; MSE is not a gate (it favours shrinking to zero).
    zero_comparisons = []
    for setting in settings:
        if setting['family'] == 'zero': continue
        comp = paired_summary(values[setting['arm']], zero, draws)
        zero_comparisons.append(dict(arm=setting['arm'], mse_diff=comp['mse_diff'],
                                     pearson_diff=comp['pearson_diff'], pearson_ci_low=comp['pearson_diff_ci_low'],
                                     top50_diff=comp['top50_overlap_diff'], top50_ci_low=comp['top50_overlap_diff_ci_low']))
    shortlisted = {x['arm'] for x in finalists}
    worse = [x for x in zero_comparisons if x['pearson_diff'] <= 0 or x['top50_diff'] <= 0]
    unclear = [x for x in zero_comparisons if x['arm'] in shortlisted and (x['pearson_ci_low'] <= 0 or x['top50_ci_low'] <= 0)]
    passed = not unclear
    save_json(out / 'zero_comparisons.json', zero_comparisons)
    save_json(out / 'shortlist.json', dict(arms=finalists, ready_for_benchmark=passed))
    info.update(seed=seed, bootstrap_resamples=resamples, bootstrap_reference='re_a1_m0',
                bootstrap_indices_sha256=hashlib.sha256(draws.tobytes()).hexdigest(),
                inputs=dict(predictions=str(p.resolve()), truth=str(Path(truth).resolve()),
                            filters=str(Path(filters).resolve()), targets=str(Path(targets_csv).resolve())),
                n_arms=len(frame), zero_mse=float(zero[:, 0].mean()), zero_check_passed=passed,
                arms_not_better_than_zero=len(worse), shortlisted_not_clearly_better_than_zero=len(unclear),
                zero_check='Every shortlisted arm must beat zero on Pearson and top-50 overlap with paired 95% CIs above 0; MSE is reported, not gated.',
                transform='vcc3_generate.load_deltas + calibrate: minimum size, alpha, +/-8, ctrl_cpm > 5, own gene zero; nonfinite predictions -> 0',
                scoring='Gate genes with finite truth, excluding own gene; equal target weight; no clipping of measured truth.',
                pearson_constant=0, top50='nonzero predictions only, ties in stable gene order; denominator min(50, eligible genes)',
                direction='accuracy on significant genes the arm predicts; coverage = share of significant genes predicted',
                selection='array and minimum size by Pearson at alpha = 1; alpha is chosen by the official benchmark')
    save_json(out / 'screen_report.json', info)
    cols = ['arm', *METRICS, 'pearson_diff_ci_low', 'pearson_diff_ci_high', 'top50_overlap_diff_ci_low']
    ordered = frame.sort_values(['pearson', 'arm'], ascending=[False, True])
    summary = ordered[cols].to_string(index=False, float_format=lambda x: f'{x:.6f}')
    summary += '\n\nShortlist:\n' + json.dumps(finalists, indent=2)
    summary += (f'\n\nZero gate: {"PASS" if passed else "FAIL"}; {len(worse)} arms are not above zero on Pearson and top-50, '
                f'{len(unclear)} shortlisted arms are not clearly above it.\n')
    (out / 'summary.txt').write_text(summary)
    print(summary, flush=True)
    return info


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('predictions', 'truth', 'filters', 'targets', 'out'): parser.add_argument('--' + arg, required=True)
    parser.add_argument('--seed', type=int, default=2026)
    parser.add_argument('--resamples', type=int, default=1000)
    args = parser.parse_args()
    report = screen(args.predictions, args.truth, args.filters, args.targets, args.out, args.seed, args.resamples)
    if not report['zero_check_passed']:
        raise SystemExit('STOP: a shortlisted arm does not clearly beat zero on Pearson and top-50. Inspect screen_report.json and alignment before benchmarking.')

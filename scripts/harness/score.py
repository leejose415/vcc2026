"""Score a prediction on the local harness evals with official cell-eval2 (vcc2026, pdex) and the eval's cached b / r.

  # 400-cell count prediction (h5ad): obs `target` or `target_gene` (controls 'non-targeting', at least one control
  # cell), every held-out knockdown of the eval present; genes by symbol (reindexed to the panel; missing = 0)
  ~/venv/bin/python scripts/harness/score.py --run my_model --eval A B --cells pred_A.h5ad pred_B.h5ad
  # effect table (the generator's prediction-folder format: targets.npy, genes.npy, <array>.npy = T x G log2 effects,
  # optional <array>_se.npy, ontarget_lfc.npy, and nochange_targets.txt = knockdowns emitted as pure control cells),
  # turned into 400 cells per knockdown by the existing generator with the fixed settings in generator_settings.json
  ~/venv/bin/python scripts/harness/score.py --run my_model --eval A A_shared B --effects predA/ predA/ predB/ --array re

One prediction per eval (the same path may repeat; extra knockdowns in it are dropped, e.g. Eval A's file for
A_shared). Writes results/harness/<run>/eval<name>/: agg_results.csv (raw, cell-eval2), scaled.csv (b / r from the
eval's bundle, 'from_replicate' = the leaderboard scale), per_target.csv, summary.json; and
results/harness/<run>/summary.csv (one row per eval: overall + the six metrics, scaled and raw). Generated / normalized
cell files are deleted after scoring unless --keep-cells.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import (LIB, RES, CONTROL, G, EVALS, eval_dir, panel_genes, evaluator, targets_of, generator_line,  # noqa: E402
                    write_h5ad, log, save_json, seed_for, read_json)

SETTINGS = LIB / 'generator_settings.json'


def normalize_cells(name, pred, out):
    """A conforming copy of a cell prediction for one eval (or the input itself when it already conforms)."""
    import anndata as ad
    genes = panel_genes(); targets, _ = targets_of(name)
    a = ad.read_h5ad(pred, backed='r')
    try:
        cols = [c for c in ('target_gene', 'target') if c in a.obs]
        if not cols: raise ValueError(f'{pred}: obs needs `target` or `target_gene`')
        lab = a.obs[cols[0]].astype(str).values
        if len(cols) == 2 and not (a.obs['target'].astype(str).values == lab).all():
            raise ValueError(f'{pred}: obs `target` and `target_gene` disagree')
        present = set(lab)
        missing = sorted(set(targets) - present)
        if missing: raise ValueError(f'{pred}: {len(missing)} held-out knockdowns of eval{name} missing, e.g. {missing[:5]}')
        if CONTROL not in present: raise ValueError(f'{pred}: no {CONTROL!r} cells (cell-eval2 needs the control label)')
        keep = np.isin(lab, list(targets) + [CONTROL])
        counts = pd.Series(lab[keep]).value_counts()
        cells = {t: int(counts[t]) for t in targets}
        conforms = (keep.all() and cols[0] == 'target_gene' and a.var_names.astype(str).tolist() == genes)
        if conforms:
            return Path(pred), cells, False
        sub = a[np.flatnonzero(keep)].to_memory()
    finally:
        a.file.close()
    X = G.aligned_X(sub, genes); G.validate_counts(X, str(pred))
    path = out / 'pred_cells.h5ad'
    labels = sub.obs[cols[0]].astype(str).values
    write_h5ad(path, X, pd.DataFrame(dict(target=labels, target_gene=labels), index=sub.obs_names.astype(str)), genes)
    return path, cells, True


def generate_cells(name, folder, array, out, settings):
    """An effect-table prediction -> 400 cells per held-out knockdown, through the existing generator, exactly as
    vcc3_generate.cmd_benchmark builds an arm (calibrate -> common control pick -> write_arm)."""
    genes = panel_genes(); targets, real_counts = targets_of(name)
    cfg = G.Cfg(argparse.Namespace(**settings))
    line = generator_line(name)
    D = G.load_deltas(folder, targets, genes, array); info = G.load_info(folder, targets, genes, array, cfg.gate_from)
    nc_file = Path(folder) / 'nochange_targets.txt'
    nochange = set(nc_file.read_text().split()) if nc_file.exists() else set()
    kd = G.kd_residual(G.prediction_ontarget(folder, targets), cfg)
    kd_list = [None if t in nochange else kd[i] for i, t in enumerate(targets)]
    for i, t in enumerate(targets):
        if t in nochange: D[i] = 0.0
    calibrated = G.calibrate(D, line, targets, cfg, info=info)
    pool = G.Cfg(argparse.Namespace(**settings)); pool.n_cells = max([cfg.n_cells] + list(real_counts.values()))
    common = G.pick_cells(line, np.zeros(len(genes)), pool, np.random.default_rng(seed_for(cfg.seed, 'H1', 'base_controls')))
    path = out / 'pred_cells.h5ad'
    G.write_arm(path, line, calibrated, kd_list, targets, {t: cfg.n_cells for t in targets}, common, cfg)
    return path, {t: cfg.n_cells for t in targets}, dict(nochange_targets=sorted(nochange & set(targets)),
                                                         called_genes_median=float(np.median((calibrated != 0).sum(1))))


def raw_metrics(arm_dir):
    d = pd.DataFrame({'metric': list(G.VCC2026_METRICS), 'direction': list(G.VCC2026_METRICS.values())})
    path = arm_dir / 'metric_directions.csv'; d.to_csv(path, index=False)
    frame, _ = G.metric_rank_table({'arm': str(arm_dir)}, path); path.unlink()
    return dict(zip(frame.metric, frame.value))


def score_one(name, run, cells=None, effects=None, array='re', settings=None, threads=8, keep=False):
    out = RES / run / f'eval{name}'; out.mkdir(parents=True, exist_ok=True)
    if cells is not None:
        pred, n_cells, temp = normalize_cells(name, cells, out); extra = dict(input='cells', source=str(cells))
    else:
        settings = settings or read_json(SETTINGS)
        pred, n_cells, gen = generate_cells(name, effects, array, out, settings); temp = True
        extra = dict(input='effects', source=str(effects), array=array, generator_settings=settings, **gen)
    log(f'eval{name}: scoring {pred}')
    ev = evaluator(name, out=out, threads=threads)
    ev.score(pred, out); scaled = ev.scale(out)
    try:
        G.per_target_table(out).to_csv(out / 'per_target.csv')
    except Exception as e:                                                    # diagnostic only
        log(f'  per-target table skipped: {e}')
    raw = raw_metrics(out)
    if temp and not keep: Path(pred).unlink()
    summary = dict(eval=name, run=run, scaled=scaled, raw=raw, n_knockdowns=len(n_cells),
                   cells_per_knockdown=dict(min=min(n_cells.values()), max=max(n_cells.values())),
                   bundle=str(eval_dir(name) / 'bundle'), **extra)
    save_json(out / 'summary.json', summary)
    log(f'eval{name}: scaled overall {scaled["overall"]:.4f}')
    return summary


def write_run_summary(run):
    rows = []
    for name in EVALS:
        p = RES / run / f'eval{name}' / 'summary.json'
        if not p.exists(): continue
        s = read_json(p)
        rows.append(dict(eval=name, overall=s['scaled']['overall'], **{f'scaled:{m}': s['scaled'].get(m) for m in G.VCC2026_METRICS},
                         **{f'raw:{m}': s['raw'].get(m) for m in G.VCC2026_METRICS}))
    pd.DataFrame(rows).to_csv(RES / run / 'summary.csv', index=False)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--run', required=True); p.add_argument('--eval', nargs='+', required=True, choices=list(EVALS))
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--cells', nargs='+'); g.add_argument('--effects', nargs='+')
    p.add_argument('--array', default='re'); p.add_argument('--settings', help=f'generator settings JSON (default {SETTINGS})')
    p.add_argument('--threads', type=int, default=8); p.add_argument('--keep-cells', action='store_true')
    a = p.parse_args()
    preds = a.cells or a.effects
    if len(preds) != len(a.eval): sys.exit('error: give one prediction per --eval (repeat a path to reuse it)')
    settings = read_json(a.settings) if a.settings else None
    for name, pred in zip(a.eval, preds):
        score_one(name, a.run, cells=pred if a.cells else None, effects=pred if a.effects else None,
                  array=a.array, settings=settings, threads=a.threads, keep=a.keep_cells)
    write_run_summary(a.run)

"""One run, no salt (Jose 2026-10-27): the noise that matters is which knockdowns are scored, not the generator's draw.
On RPE1dev C+G the three salts gave identical per-knockdown values on 8 of 9 metrics (only direction_reach moved,
r = 0.991), and the knockdown bootstrap SD was 3-100x the salt range. So an end-to-end run is scored once (salt "")
and its uncertainty is a bootstrap over knockdowns.

Per scored harness run (the folder with results.csv / agg_results.csv / run_meta.json / summary.json): the per-knockdown
values of cell-eval2's results.csv are resampled (2,000 knockdown resamples, seed 2026), the aggregate rebuilt as
cell-eval2 builds it (mean over the knockdowns with a value; expr_mse_unbiased_capped_norm = the ratio of sums
capped / distance), and scaled by cell-eval2's own score_metrics against the eval's bundle (from_replicate, the
leaderboard scale). The bundle's baseline and replicate anchors stay fixed. The full-set rebuild must reproduce the
harness's summary.json exactly (checked on every run). Two runs on the same eval are compared paired (the same
resampled knockdowns): candidate - reference per metric and overall, with a 95% CI and P(candidate better).

  ~/venv/bin/python code/vcc6/boot_score.py RUN_DIR [CANDIDATE_RUN_DIR ...] [--out FILE.json]
"""
import argparse
import json
import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

N_BOOT, SEED = 2000, 2026
RATIO = ('expr_mse_unbiased_capped_norm', 'expr_mse_unbiased_capped', 'expr_distance_unbiased')


class Run:
    def __init__(self, path):
        self.dir = Path(path)
        self.summary = json.loads((self.dir / 'summary.json').read_text())
        self.meta = json.loads((self.dir / 'run_meta.json').read_text())
        self.bundle = self.summary['bundle']
        self.template = pl.read_csv(self.dir / 'agg_results.csv').filter(pl.col('statistic') == 'mean')
        self.per_kd = (pd.read_csv(self.dir / 'results.csv')
                       .pivot(index='perturbation', columns='metric', values='value'))

    def aggregate(self, idx):
        """cell-eval2's mean row from the knockdowns idx (labels; repeats allowed)."""
        P = self.per_kd.loc[idx]
        row = {}
        for m in self.template.columns[1:]:
            if m == RATIO[0]:
                row[m] = float(P[RATIO[1]].sum() / P[RATIO[2]].sum())
            elif m in P:
                row[m] = float(P[m].mean())                      # pandas skips NaN, as cell-eval2's n_used does
            else:
                row[m] = float(self.template[m][0])
        return pl.DataFrame({'statistic': ['mean'], **{m: [row[m]] for m in self.template.columns[1:]}},
                            schema=self.template.schema)

    def scaled(self, idx):
        from cell_eval2.score import score_metrics
        logging.getLogger('cell_eval2').setLevel(logging.ERROR)     # 'metric not scored' notes, once per call
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            out = score_metrics(self.aggregate(idx), real_bundle=self.bundle, user_meta=self.meta)
        d = dict(zip(out['metric'].to_list(), out['from_replicate'].to_list()))
        d['overall'] = d.pop('avg_score')
        return d

    def check(self):
        """The full set must reproduce the harness's own scaled numbers."""
        full = self.scaled(self.per_kd.index)
        ref = self.summary['scaled']
        bad = {m: (full[m], ref[m]) for m in ref if not np.isclose(full[m], ref[m], rtol=0, atol=1e-9)}
        if bad:
            raise AssertionError(f'{self.dir}: rebuild does not reproduce summary.json: {bad}')
        return full


def boot(runs, n_boot=N_BOOT):
    """Point estimates, bootstrap CIs, and paired differences of runs[1:] against runs[0]."""
    common = runs[0].per_kd.index
    for r in runs[1:]:
        common = common.intersection(r.per_kd.index)
    rng = np.random.default_rng(SEED)
    draws = [[] for _ in runs]
    for _ in range(n_boot):
        idx = common[rng.integers(0, len(common), len(common))]
        for i, r in enumerate(runs):
            draws[i].append(r.scaled(idx))
    frames = [pd.DataFrame(d) for d in draws]
    res = dict(knockdowns=int(len(common)), n_boot=n_boot, runs={})
    for r, f in zip(runs, frames):
        full = r.check()
        res['runs'][str(r.dir)] = {m: dict(value=full[m], ci=[float(f[m].quantile(0.025)), float(f[m].quantile(0.975))],
                                           sd=float(f[m].std())) for m in f.columns}
    res['compare'] = {}
    for r, f in zip(runs[1:], frames[1:]):
        d = f - frames[0]
        res['compare'][f'{r.dir} vs {runs[0].dir}'] = {
            m: dict(diff=float(r.check()[m] - runs[0].check()[m]),
                    ci=[float(d[m].quantile(0.025)), float(d[m].quantile(0.975))], p_better=float((d[m] > 0).mean()))
            for m in d.columns}
    return res


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('runs', nargs='+', help='scored harness folders; the first is the reference')
    p.add_argument('--n-boot', type=int, default=N_BOOT)
    p.add_argument('--out')
    a = p.parse_args()
    res = boot([Run(x) for x in a.runs], a.n_boot)
    text = json.dumps(res, indent=1) + '\n'
    if a.out:
        Path(a.out).write_text(text)
    for name, r in res['runs'].items():
        print(f'{name}: overall {r["overall"]["value"]:.4f} [{r["overall"]["ci"][0]:.4f}, {r["overall"]["ci"][1]:.4f}]')
    for name, c in res['compare'].items():
        print(f'{name}: overall {c["overall"]["diff"]:+.4f} [{c["overall"]["ci"][0]:+.4f}, {c["overall"]["ci"][1]:+.4f}] '
              f'P(better) {c["overall"]["p_better"]:.3f}')
    sys.exit(0)

"""Tables for the RPE1 dev run (results/stage1/2026-10-27_rpe1/): harness summaries, salt spread, effect-size terciles.

Per tercile only raw values exist (the scaled scale is anchored on the whole panel): the normalised MSE as the scorer
forms it (sum of expr_mse_unbiased_capped / sum of expr_distance_unbiased over the tercile's knockdowns), every other
member as the mean of its per-knockdown values in results.csv (each member's own aggregation, mean, NaN-skipping).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
R = V / 'results/stage1/2026-10-27_rpe1'
EVAL = 'evalRPE1dev'
METRICS = ['pds_cosine', 'expr_mse_unbiased_capped_norm', 'de_wilcoxon_lfc_nmae', 'de_wilcoxon_direction_fidelity_yield_raw',
           'de_wilcoxon_direction_reach_raw', 'de_wilcoxon_sig_jaccard']
SHORT = dict(pds_cosine='pds', expr_mse_unbiased_capped_norm='mse_norm', de_wilcoxon_lfc_nmae='lfc_nmae',
             de_wilcoxon_direction_fidelity_yield_raw='dir_fidelity', de_wilcoxon_direction_reach_raw='dir_reach',
             de_wilcoxon_sig_jaccard='sig_jaccard')


def run_dirs():
    """{label: harness output dir} for every scored run that exists."""
    out = {}
    for arm in ('A', 'B', 'C'):
        d = R / 'harness/results' / f'step1_arm{arm}' / EVAL
        if (d / 'summary.json').exists(): out[f'step1 {arm}'] = d
    for arm in ('zero', 'oracle', 'armB', 'armC'):
        for s in ('none', 's1', 's2'):
            d = R / 'step3' / arm / f'salt_{s}' / 'harness' / EVAL
            if (d / 'summary.json').exists(): out[f'step3 {arm} {s}'] = d
    return out


def summary_table(dirs=None):
    rows = []
    for label, d in (dirs or run_dirs()).items():
        s = json.loads((d / 'summary.json').read_text())
        r = dict(run=label, scaled_overall=s['scaled'].get('overall'))
        r.update({f'scaled_{SHORT[m]}': s['scaled'].get(m) for m in METRICS})
        r.update({f'raw_{SHORT[m]}': s['raw'].get(m) for m in METRICS})
        rows.append(r)
    return pd.DataFrame(rows).set_index('run') if rows else pd.DataFrame()


def salt_spread(table):
    """Step 3: per arm, mean and spread (max - min, sd) across the three salts."""
    t = table[table.index.str.startswith('step3')].copy()
    if t.empty: return t
    t['arm'] = t.index.str.split().str[1]
    g = t.groupby('arm')
    num = t.select_dtypes('number').columns
    return pd.concat({'mean': g[num].mean(), 'range': g[num].max() - g[num].min(), 'sd': g[num].std(ddof=1)}, axis=1)


def terciles():
    p = R / 'step1/knockdowns_300.csv'
    return pd.read_csv(p, index_col=0)['tercile'] if p.exists() else None


def tercile_table(dirs=None):
    terc = terciles()
    if terc is None: return pd.DataFrame()
    rows = []
    for label, d in (dirs or run_dirs()).items():
        rc = pd.read_csv(d / 'results.csv').pivot_table(index='perturbation', columns='metric', values='value')
        rc['tercile'] = terc.reindex(rc.index)
        for k, g in rc.groupby('tercile'):
            r = dict(run=label, tercile=k, knockdowns=len(g),
                     mse_norm=float(g['expr_mse_unbiased_capped'].sum() / g['expr_distance_unbiased'].sum()))
            for m in METRICS:
                if m != 'expr_mse_unbiased_capped_norm' and m in g:
                    r[SHORT[m]] = float(g[m].mean())
            rows.append(r)
    return pd.DataFrame(rows).set_index(['run', 'tercile']) if rows else pd.DataFrame()

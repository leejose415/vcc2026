"""Perturbation-level QC for sources we hold only as per-knockdown tables, not cells: X-Atlas/Orion HCT116 / HEK293T (all
knockdowns, from their per-knockdown sums) and the Marson CD4 T-cell tables (Rest, Stim8hr, Stim48hr).

The on-target rule of scripts/qc_perturbations.py, applied to each knockdown's effect on its own target gene in the table
built by scripts/knockdown_filters.py (lib/filters/_cache/<name>.csv), plus MIN_CELLS where the table records cell counts:
  measurable = the target is measured and ctrl_cpm >= MIN_CTRL_CPM
  keep = n_cells >= MIN_CELLS and ((not measurable) or (ontarget_lfc <= MAX_ONTARGET_LFC and p_one_sided < P_MAX))
The tables carry a two-sided p from their own test (Orion: Welch t on per-cell CPM; T cells: moderated t on replicate
SEs), not the per-cell Mann-Whitney U used for cell files; it is converted to one-sided 'lower': p / 2 when the change is
negative, 1 - p / 2 otherwise. There are no constructs here: one row per knockdown (guides already pooled), and no cells
for a downstream sanity check.

Writes results/qc/<name>/knockdown_qc.csv and perturbation_qc.json.
Usage: ~/venv/bin/python scripts/qc_table_perturbations.py [--sources HCT116 HEK293T CD4_Rest CD4_Stim8hr CD4_Stim48hr]
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
CACHE = V / 'lib/filters/_cache'
OUT = {'HCT116': 'orion_hct116_table', 'HEK293T': 'orion_hek293t_table',
       'CD4_Rest': 'marson_cd4_rest', 'CD4_Stim8hr': 'marson_cd4_stim8hr', 'CD4_Stim48hr': 'marson_cd4_stim48hr'}
MIN_CELLS, MIN_CTRL_CPM, MAX_ONTARGET_LFC, P_MAX = 30, 5.0, -0.5, 0.05


def main(name):
    t = pd.read_csv(CACHE / f'{name}.csv')
    lfc = t['ontarget_lfc'].astype(float); p2 = t['pvalue'].astype(float)
    t['p_one_sided'] = np.where(lfc < 0, p2 / 2, 1 - p2 / 2)
    measured = t['target_measured'].astype(bool)
    t['measurable'] = measured & (t['ctrl_cpm'].astype(float) >= MIN_CTRL_CPM)
    has_cells = t['n_cells'].notna()
    t['enough_cells'] = ~has_cells | (t['n_cells'] >= MIN_CELLS)
    on_target = (lfc <= MAX_ONTARGET_LFC) & (t['p_one_sided'] < P_MAX)
    t['knockdown_keep'] = t['enough_cells'] & (~t['measurable'] | on_target)
    out = V / 'results/qc' / OUT[name]; out.mkdir(parents=True, exist_ok=True)
    cols = ['target', 'n_cells', 'ctrl_cpm', 'measurable', 'ontarget_lfc', 'pvalue', 'p_one_sided', 'enough_cells', 'knockdown_keep', 'test']
    t[cols].rename(columns={'target': 'target_gene', 'pvalue': 'p_two_sided_source_test'}).to_csv(out / 'knockdown_qc.csv', index=False)
    m = t[t.measurable]; kept = t[t.knockdown_keep]
    summary = dict(source=OUT[name], table=str(CACHE / f'{name}.csv'), test=str(t['test'].iloc[0]),
                   knockdowns_total=len(t), cell_counts_recorded=int(has_cells.sum()), knockdowns_measurable=int(t.measurable.sum()),
                   knockdowns_kept=int(t.knockdown_keep.sum()), kept_unmeasurable=int((t.knockdown_keep & ~t.measurable).sum()),
                   kept_between_minus1_and_minus05=int((t.knockdown_keep & t.measurable & (lfc > -1)).sum()),
                   kept_le_minus1=int((t.knockdown_keep & t.measurable & (lfc <= -1)).sum()),
                   dropped_measurable=int((t.measurable & ~t.knockdown_keep).sum()),
                   knockdowns_20_29_cells=int(((t.n_cells >= 20) & (t.n_cells < MIN_CELLS)).sum()),
                   median_efficiency_kept=float(kept.loc[kept.measurable, 'ontarget_lfc'].median()),
                   parameters=dict(MIN_CELLS=MIN_CELLS, MIN_CTRL_CPM=MIN_CTRL_CPM, MAX_ONTARGET_LFC=MAX_ONTARGET_LFC, P_MAX=P_MAX))
    (out / 'perturbation_qc.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sources', nargs='+', default=list(OUT), choices=list(OUT))
    rows = [main(n) for n in p.parse_args().sources]
    print(pd.DataFrame(rows).drop(columns=['parameters', 'table']).set_index('source').T.to_string())

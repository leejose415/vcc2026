"""Apply cell-level QC to the Replogle 2022 K562 genome-wide file; writes barcode lists, never an .h5ad.

Reads only obs (UMI_count, gene, gem_group) and the knockdowns flagged by qc_replogle_umi.py. On x = log10(UMI_count + 1),
with the scaled MAD (1.4826 x MAD, floor 0.05; 1 scaled MAD = 1 SD for a normal distribution, as scater's isOutlier):
  low = median - 3 MAD,  high = median + 4 MAD
Reference cells: per gem group, its controls ('non-targeting'), or all controls when the lane has fewer than 200;
for each flagged knockdown, its own cells (low cutoff only).
Rules: 1. every cell: removed below 200 UMIs ('floor') or above its lane's high cutoff ('high');
       2. controls and unflagged knockdowns: removed below the lane's low cutoff ('low_control');
       3. flagged knockdowns: removed below the knockdown's own low cutoff instead ('low_own').

Writes results/qc/replogle_k562_gw/: kept_cells.txt, removed_cells.csv, qc_method.md.
Usage: ~/venv/bin/python scripts/qc_replogle_apply.py   (after scripts/qc_replogle_umi.py)
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from qc_replogle_umi import (CONTROL, FILE, HIGH_K, LOW_K, MAD_FLOOR, MAD_SCALE, MIN_LANE_CONTROLS, MIN_KD_CELLS, OUT,  # noqa: E402
                             UMI_FLOOR, cutoffs)


def main():
    flagged_table = pd.read_csv(OUT / 'knockdowns_flagged.csv')
    flagged = set(flagged_table.loc[flagged_table.flagged, 'knockdown'].astype(str))
    a = ad.read_h5ad(FILE, backed='r'); obs = a.obs[['gene', 'gem_group', 'UMI_count']].copy(); a.file.close()
    obs['gene'] = obs['gene'].astype(str); obs['gem_group'] = obs['gem_group'].astype(str)
    umi = obs['UMI_count'].to_numpy(float); x = np.log10(umi + 1); ctrl = (obs['gene'] == CONTROL).to_numpy()

    whole = cutoffs(x[ctrl])
    lane = {}
    for g, idx in obs.groupby('gem_group').indices.items():
        c = idx[ctrl[idx]]
        lane[g] = cutoffs(x[c]) if len(c) >= MIN_LANE_CONTROLS else whole
    lane_low = obs['gem_group'].map({g: v['low_umi'] for g, v in lane.items()}).to_numpy(float)
    lane_high = obs['gem_group'].map({g: v['high_umi'] for g, v in lane.items()}).to_numpy(float)
    own = {k: cutoffs(x[idx])['low_umi'] for k, idx in obs.groupby('gene').indices.items() if k in flagged}
    is_flagged = obs['gene'].isin(flagged).to_numpy()
    own_low = obs['gene'].map(own).to_numpy(float)

    floor = umi < UMI_FLOOR
    high = ~floor & (umi > lane_high)
    low_control = ~floor & ~high & ~is_flagged & (umi < lane_low)
    low_own = ~floor & ~high & is_flagged & (umi < own_low)
    removed = floor | high | low_control | low_own
    reason = np.select([floor, high, low_control, low_own], ['floor', 'high', 'low_control', 'low_own'], '')
    before = floor | (umi > lane_high) | (umi < lane_low)                      # the same cells under rules 1-2 only

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'kept_cells.txt').write_text('\n'.join(obs.index[~removed].astype(str)) + '\n')
    out = obs.loc[removed, ['gene', 'gem_group', 'UMI_count']].copy(); out.insert(0, 'cell_barcode', out.index.astype(str))
    out['reason'] = reason[removed]; out.to_csv(OUT / 'removed_cells.csv', index=False)

    kd = ~ctrl
    stats = dict(cells=len(obs), kept=int((~removed).sum()), controls=int(ctrl.sum()), knockdown_cells=int(kd.sum()),
                 ctrl_removed=float(removed[ctrl].mean()), kd_removed=float(removed[kd].mean()),
                 by_reason={r: int((reason == r).sum()) for r in ('floor', 'high', 'low_control', 'low_own')},
                 lanes=len(lane), lanes_fallback=sum(v is whole for v in lane.values()), flagged=len(flagged))
    per = pd.DataFrame(dict(gene=obs['gene'], before=before, after=removed))[is_flagged].groupby('gene').mean()
    top = flagged_table.loc[flagged_table.flagged].head(10)['knockdown'].astype(str)
    comparison = per.loc[top].rename(columns=dict(before='removed_before_rule3', after='removed_after_rule3'))
    comparison.insert(0, 'cells', obs['gene'].value_counts().loc[top]); comparison['own_low_umi'] = [own[g] for g in top]

    print(f"cells kept: {stats['kept']:,} of {stats['cells']:,} ({stats['kept'] / stats['cells']:.2%})")
    print(f"controls removed: {stats['ctrl_removed']:.2%}; knockdown cells removed: {stats['kd_removed']:.2%}; by reason: {stats['by_reason']}")
    print('top 10 flagged knockdowns, share removed before vs after rule 3:'); print(comparison.round(3).to_string())

    rows = '\n'.join(f"| {g} | {int(r.cells)} | {r.removed_before_rule3:.1%} | {r.removed_after_rule3:.1%} | {r.own_low_umi:,.0f} |"
                     for g, r in comparison.iterrows())
    (OUT / 'qc_method.md').write_text(f"""# Cell QC: Replogle 2022 K562 genome-wide

Source: `{FILE}` (read only; obs columns `UMI_count`, `gene`, `gem_group`). Scripts: `scripts/qc_replogle_umi.py` (distribution,
cutoffs, flagged knockdowns), `scripts/qc_replogle_apply.py` (this filter). No .h5ad is written; apply `kept_cells.txt` when loading.

## Definitions
On x = log10(UMI_count + 1), over the reference cells:
- median = median(x)
- scaled MAD = {MAD_SCALE} × median(|x − median|), floor {MAD_FLOOR} (1 scaled MAD = 1 SD for a normal distribution, as scater's `isOutlier`); the raw MAD is never used
- low cutoff = median − {LOW_K} × scaled MAD; high cutoff = median + {HIGH_K} × scaled MAD

Reference cells: per gem group (10x lane), that lane's control cells (`{CONTROL}`); a lane with fewer than {MIN_LANE_CONTROLS} controls
uses all controls in the file ({stats['lanes_fallback']} of {stats['lanes']} lanes). For each of the {stats['flagged']} flagged knockdowns
(`knockdowns_flagged.csv`: knockdowns with ≥ {MIN_KD_CELLS} cells losing more than max(2 × the control rate, control rate + 5 points) to
the lane cutoffs), that knockdown's own cells, for the low cutoff only.

## Rules
1. Every cell: removed if UMI_count < {UMI_FLOOR} (`floor`, empty droplets) or above its lane's high cutoff (`high`, likely doublets).
2. Controls and unflagged knockdowns: also removed below their lane's low cutoff (`low_control`).
3. Flagged knockdowns: removed below the knockdown's own low cutoff instead (`low_own`). Their low UMI counts are biology
   (loss of core transcription machinery: Mediator, SUPT5H / SUPT6H, TAF1, GATA1), not poor quality.

## Whole-file cutoffs (all controls)
low {whole['low_umi']:,.0f} UMIs, high {whole['high_umi']:,.0f} UMIs (median {10 ** whole['median'] - 1:,.0f}, scaled MAD {whole['mad']:.3f} log10).
Lane cutoffs vary with depth: control medians span about 2,900–20,800 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: {stats['cells']:,}; kept {stats['kept']:,} ({stats['kept'] / stats['cells']:.2%}); removed {stats['cells'] - stats['kept']:,}
- removed by reason: floor {stats['by_reason']['floor']:,}, high {stats['by_reason']['high']:,}, low_control {stats['by_reason']['low_control']:,}, low_own {stats['by_reason']['low_own']:,}
- controls removed: {stats['ctrl_removed']:.2%} of {stats['controls']:,}; knockdown cells removed: {stats['kd_removed']:.2%} of {stats['knockdown_cells']:,}

Top 10 flagged knockdowns, share of cells removed with the lane low cutoff (before rule 3) and with their own (after):

| knockdown | cells | before | after | own low cutoff (UMIs) |
|---|---|---|---|---|
{rows}
""")
    print(f'wrote {OUT}')


if __name__ == '__main__':
    main()

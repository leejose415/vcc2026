"""Step 2 of cell QC: apply the cutoffs from qc_replogle_umi.py; writes barcode lists, never an .h5ad.

Reads only the cell metadata (gene, lane, UMIs: see qc_replogle_umi.load_cells) and step 1's outputs. On x = log10(UMI + 1),
with the scaled MAD (1.4826 x MAD, floor 0.05; 1 scaled MAD = 1 SD for a normal distribution, as scater's isOutlier):
  low = median - 3 MAD,  high = median + 4 MAD
Reference cells: per lane, its controls ('non-targeting'), or all controls when the lane has fewer than LANE_MIN;
for each flagged knockdown, its own cells (low cutoff only).
Rules: 1. every cell: removed below 200 UMIs ('floor') or above its lane's high cutoff ('high');
       2. controls and unflagged knockdowns: removed below the lane's low cutoff ('low_control');
       3. flagged knockdowns: removed below the knockdown's own low cutoff instead ('low_own').
A source with a sharp left edge (left_edge.json: pre-filtered by its authors) gets rule 1 only: no low cut.

Writes results/qc/<source>/: kept_cells.txt, removed_cells.csv, qc_method.md.
Usage: ~/venv/bin/python scripts/qc_replogle_apply.py --dataset K562|RPE1|HepG2|Jurkat|H1   (after scripts/qc_replogle_umi.py)
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from qc_replogle_umi import (CONTROL, HIGH_K, LANE_MIN, LOW_K, MAD_FLOOR, MAD_SCALE, MIN_LANE_CONTROLS, MIN_KD_CELLS, SOURCES,  # noqa: E402
                             UMI_FLOOR, cutoffs, load_cells)


def main(dataset='K562'):
    FILE, OUT, _, title = SOURCES[dataset]; min_lane = LANE_MIN.get(dataset, MIN_LANE_CONTROLS)
    edge = json.loads((OUT / 'left_edge.json').read_text()); skip_low = bool(edge['sharp'])
    flagged_table = pd.read_csv(OUT / 'knockdowns_flagged.csv')
    flagged = set() if skip_low else set(flagged_table.loc[flagged_table.flagged, 'knockdown'].astype(str))
    _, obs = load_cells(dataset)
    umi = obs['UMI_count'].to_numpy(float); x = np.log10(umi + 1); ctrl = (obs['gene'] == CONTROL).to_numpy()

    whole = cutoffs(x[ctrl])
    lane = {}
    for g, idx in obs.groupby('gem_group').indices.items():
        c = idx[ctrl[idx]]
        lane[g] = cutoffs(x[c]) if len(c) >= min_lane else whole
    lane_low = obs['gem_group'].map({g: v['low_umi'] for g, v in lane.items()}).to_numpy(float)
    lane_high = obs['gem_group'].map({g: v['high_umi'] for g, v in lane.items()}).to_numpy(float)
    own = {k: cutoffs(x[idx])['low_umi'] for k, idx in obs.groupby('gene').indices.items() if k in flagged}
    is_flagged = obs['gene'].isin(flagged).to_numpy()
    own_low = obs['gene'].map(own).to_numpy(float)

    floor = umi < UMI_FLOOR
    high = ~floor & (umi > lane_high)
    low_control = ~floor & ~high & ~is_flagged & (umi < lane_low) & (not skip_low)
    low_own = ~floor & ~high & is_flagged & (umi < own_low) & (not skip_low)
    removed = floor | high | low_control | low_own
    reason = np.select([floor, high, low_control, low_own], ['floor', 'high', 'low_control', 'low_own'], '')
    before = floor | (umi > lane_high) | (umi < lane_low)                      # rules 1-2 for every cell (the first pass)

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'kept_cells.txt').write_text('\n'.join(obs.index[~removed].astype(str)) + '\n')
    out = obs.loc[removed, ['gene', 'gem_group', 'UMI_count']].copy(); out.insert(0, 'cell_barcode', out.index.astype(str))
    out['reason'] = reason[removed]; out.to_csv(OUT / 'removed_cells.csv', index=False)

    kd = ~ctrl
    stats = dict(cells=len(obs), kept=int((~removed).sum()), controls=int(ctrl.sum()), knockdown_cells=int(kd.sum()),
                 ctrl_removed=float(removed[ctrl].mean()), kd_removed=float(removed[kd].mean()),
                 by_reason={r: int((reason == r).sum()) for r in ('floor', 'high', 'low_control', 'low_own')},
                 lanes=len(lane), lanes_fallback=sum(v is whole for v in lane.values()), flagged=len(flagged))
    # the knockdowns step 1 flagged, first-pass share vs share removed by the final rules
    top = flagged_table.loc[flagged_table.flagged].head(10)['knockdown'].astype(str).tolist()
    per = pd.DataFrame(dict(gene=obs['gene'], before=before, after=removed)).groupby('gene').mean()
    comparison = per.reindex(top).rename(columns=dict(before='removed_first_pass', after='removed_final'))
    comparison.insert(0, 'cells', obs['gene'].value_counts().reindex(top).to_numpy())
    comparison['own_low_umi'] = [own.get(g, np.nan) for g in top]

    print(f"cells kept: {stats['kept']:,} of {stats['cells']:,} ({stats['kept'] / stats['cells']:.2%})"
          f"{' | low cut SKIPPED: sharp left edge (pre-filtered by the authors)' if skip_low else ''}")
    print(f"controls removed: {stats['ctrl_removed']:.2%}; knockdown cells removed: {stats['kd_removed']:.2%}; by reason: {stats['by_reason']}")
    print('top 10 flagged knockdowns (step 1), share removed in the first pass vs by the final rules:'); print(comparison.round(3).to_string())

    medians = [v['median'] for v in lane.values() if v is not whole]
    lane_lo, lane_hi = (10 ** min(medians) - 1, 10 ** max(medians) - 1) if medians else (0, 0)
    rows = '\n'.join(f"| {g} | {r.cells:,.0f} | {r.removed_first_pass:.1%} | {r.removed_final:.1%} | "
                     f"{'—' if np.isnan(r.own_low_umi) else f'{r.own_low_umi:,.0f}'} |" for g, r in comparison.iterrows())
    prefiltered = f"""
## Pre-filtered source: no low cut
The UMI histogram has a sharp left edge (`left_edge.json`): the lowest cell has {edge['min_umi']:,.0f} UMIs and the lowest bin still holds
{edge['first_bin_share_of_peak']:.1%} of the peak bin, so the authors filtered at about {round(edge['min_umi'], -2):,.0f} UMIs before release. Only the
{UMI_FLOOR}-UMI floor and the per-lane high cut are applied (rules 2 and 3 skipped). Low-UMI cells of strong knockdowns were removed
before we received the file, so some strong effects (loss of core transcription machinery and similar) may look weaker here than in K562,
where those cells were kept by rule 3.
""" if skip_low else ''
    rules = (f"""1. Every cell: removed if UMI < {UMI_FLOOR} (`floor`, empty droplets) or above its lane's high cutoff (`high`, likely doublets).
2. (skipped: pre-filtered source) Controls and unflagged knockdowns below their lane's low cutoff.
3. (skipped: pre-filtered source) Flagged knockdowns below their own low cutoff.""" if skip_low else
             f"""1. Every cell: removed if UMI < {UMI_FLOOR} (`floor`, empty droplets) or above its lane's high cutoff (`high`, likely doublets).
2. Controls and unflagged knockdowns: also removed below their lane's low cutoff (`low_control`).
3. Flagged knockdowns: removed below the knockdown's own low cutoff instead (`low_own`). Their low UMI counts are biology
   (in K562, loss of core transcription machinery: Mediator, SUPT5H / SUPT6H, TAF1, GATA1), not poor quality.
   Flagged here: {', '.join(top)}{' ...' if len(flagged) > 10 else ''}.""")
    (OUT / 'qc_method.md').write_text(f"""# Cell QC: {title}

Source: `{FILE}` (read only; cell metadata only; UMIs from {obs.attrs['umi_source']}). Scripts: `scripts/qc_replogle_umi.py --dataset {dataset}`
(distribution, cutoffs, flagged knockdowns, left edge), `scripts/qc_replogle_apply.py --dataset {dataset}` (this filter),
`scripts/write_replogle_qc.py --dataset {dataset}` (the filtered file).

## Definitions
On x = log10(UMI + 1), over the reference cells:
- median = median(x)
- scaled MAD = {MAD_SCALE} × median(|x − median|), floor {MAD_FLOOR} (1 scaled MAD = 1 SD for a normal distribution, as scater's `isOutlier`); the raw MAD is never used
- low cutoff = median − {LOW_K} × scaled MAD; high cutoff = median + {HIGH_K} × scaled MAD

Reference cells: per lane (10x gem group / batch), that lane's control cells (`{CONTROL}`); a lane with fewer than **{min_lane}** controls
uses all controls in the file ({stats['lanes_fallback']} of {stats['lanes']} lanes). For each flagged knockdown
(`knockdowns_flagged.csv`: knockdowns with ≥ {MIN_KD_CELLS} cells losing more than max(2 × the control rate, control rate + 5 points) to
the lane cutoffs), that knockdown's own cells, for the low cutoff only.

## Rules
{rules}
{prefiltered}
## Whole-file cutoffs (all controls)
low {whole['low_umi']:,.0f} UMIs, high {whole['high_umi']:,.0f} UMIs (median {10 ** whole['median'] - 1:,.0f}, scaled MAD {whole['mad']:.3f} log10).
Lane cutoffs vary with depth: control medians span {lane_lo:,.0f}–{lane_hi:,.0f} UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: {stats['cells']:,}; kept {stats['kept']:,} ({stats['kept'] / stats['cells']:.2%}); removed {stats['cells'] - stats['kept']:,}
- removed by reason: floor {stats['by_reason']['floor']:,}, high {stats['by_reason']['high']:,}, low_control {stats['by_reason']['low_control']:,}, low_own {stats['by_reason']['low_own']:,}
- controls removed: {stats['ctrl_removed']:.2%} of {stats['controls']:,}; knockdown cells removed: {stats['kd_removed']:.2%} of {stats['knockdown_cells']:,}

Top 10 knockdowns flagged in step 1: share of cells removed in the first pass (lane low and high cuts) and by the final rules:

| knockdown | cells | first pass | final | own low cutoff (UMIs) |
|---|---|---|---|---|
{rows}
""")
    print(f'wrote {OUT}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--dataset', default='K562', choices=list(SOURCES))
    main(p.parse_args().dataset)

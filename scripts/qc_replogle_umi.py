"""Per-cell UMI distribution of the Replogle 2022 K562 genome-wide file, and a control-referenced UMI cutoff.

Reads only the cell metadata (obs); never loads the counts and never writes to the source file.

Cutoff, computed on control cells ('non-targeting') on x = log10(UMI_count + 1):
  low  = median - 3 MAD,  high = median + 4 MAD     MAD = median absolute deviation x 1.4826, so 3 MAD = 3 standard
                                                      deviations for a normal distribution; floored at 0.05
  low_umi = max(200, 10**low - 1)                    an absolute floor of 200 UMIs
(a) once for the whole file and (b) per gem_group (10x lane), falling back to the whole-file controls when a lane has
fewer than 200 control cells.

Outputs (results/qc/replogle_k562_gw/): umi_hist.png, gemgroup_medians.png, cutoffs_by_gemgroup.csv,
knockdowns_flagged.csv (knockdowns with >= 30 cells: share removed by the per-lane cutoffs, flagged when above
max(2 x the control rate, control rate + 5 points)).

Usage: ~/venv/bin/python scripts/qc_replogle_umi.py [--file PATH]
"""
import argparse
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
FILE = V / 'lib/replogle/K562_gwps_raw_singlecell_01.h5ad'
OUT = V / 'results/qc/replogle_k562_gw'
CONTROL = 'non-targeting'
LOW_K, HIGH_K, MAD_FLOOR, UMI_FLOOR, MIN_LANE_CONTROLS, MIN_KD_CELLS = 3, 4, 0.05, 200, 200, 30
MAD_SCALE = 1.4826                # MAD -> standard deviation for a normal distribution


def cutoffs(x, scale=MAD_SCALE):
    """Median, MAD (floored) and the low / high cutoffs in UMIs for log10 values x."""
    med = float(np.median(x)); mad = max(float(np.median(np.abs(x - med))) * scale, MAD_FLOOR)
    return dict(median=med, mad=mad, low_umi=max(UMI_FLOOR, 10 ** (med - LOW_K * mad) - 1), high_umi=10 ** (med + HIGH_K * mad) - 1)


def main(path, scale=MAD_SCALE):
    a = ad.read_h5ad(path, backed='r'); obs = a.obs[['UMI_count', 'gene', 'gem_group']].copy(); a.file.close()
    print('obs["gene"].value_counts().head():'); print(obs['gene'].value_counts().head().to_string())
    if CONTROL not in set(obs['gene'].astype(str)): raise ValueError(f'no {CONTROL!r} cells: check the control label')
    obs['gene'] = obs['gene'].astype(str); obs['gem_group'] = obs['gem_group'].astype(str)
    obs['x'] = np.log10(obs['UMI_count'].to_numpy(float) + 1); ctrl = obs['gene'] == CONTROL
    whole = cutoffs(obs.loc[ctrl, 'x'].to_numpy(), scale)

    rows = []
    for lane, g in obs[ctrl].groupby('gem_group', sort=False):
        c = cutoffs(g['x'].to_numpy(), scale) if len(g) >= MIN_LANE_CONTROLS else dict(whole)
        rows.append(dict(gem_group=lane, n_control_cells=len(g), **c, fallback=len(g) < MIN_LANE_CONTROLS))
    for lane in set(obs['gem_group']) - {r['gem_group'] for r in rows}:          # a lane with no controls at all
        rows.append(dict(gem_group=lane, n_control_cells=0, **whole, fallback=True))
    lanes = pd.DataFrame(rows).sort_values('gem_group', key=lambda s: pd.to_numeric(s, errors='coerce')).reset_index(drop=True)
    OUT.mkdir(parents=True, exist_ok=True)
    lanes[['gem_group', 'n_control_cells', 'median', 'mad', 'low_umi', 'high_umi']].to_csv(OUT / 'cutoffs_by_gemgroup.csv', index=False)

    lim = obs[['gem_group']].merge(lanes[['gem_group', 'low_umi', 'high_umi']], on='gem_group', how='left')
    umi = obs['UMI_count'].to_numpy(float)
    removed = (umi < lim['low_umi'].to_numpy()) | (umi > lim['high_umi'].to_numpy())
    obs['removed'] = removed
    ctrl_rate = float(obs.loc[ctrl, 'removed'].mean()); kd = ~ctrl
    per = obs[kd].groupby('gene').agg(n_cells=('removed', 'size'), n_removed=('removed', 'sum'))
    per = per[per.n_cells >= MIN_KD_CELLS].copy(); per['share_removed'] = per.n_removed / per.n_cells
    bar = max(2 * ctrl_rate, ctrl_rate + 0.05); per['flagged'] = per.share_removed > bar
    per = per.sort_values('share_removed', ascending=False)
    per.reset_index().rename(columns={'gene': 'knockdown'}).to_csv(OUT / 'knockdowns_flagged.csv', index=False)

    # 1. histogram
    low_x, high_x = np.log10(whole['low_umi'] + 1), np.log10(whole['high_umi'] + 1)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    bins = np.linspace(obs['x'].min(), obs['x'].max(), 81)
    ax.hist(obs.loc[ctrl, 'x'], bins=bins, density=True, histtype='step', color='#2a78d6', lw=1.6, label=f'controls (n={ctrl.sum():,})')
    ax.hist(obs.loc[kd, 'x'], bins=bins, density=True, histtype='step', color='#eb6834', lw=1.6, label=f'knockdown cells (n={kd.sum():,})')
    for v, name in ((low_x, 'low'), (high_x, 'high')):
        ax.axvline(v, ls='--', color='grey', lw=1); ax.text(v, ax.get_ylim()[1] * 0.97, f' {name}: {10 ** v - 1:,.0f}', color='grey', va='top', fontsize=8)
    ax.set_xlabel('log10 total UMIs per cell'); ax.set_ylabel('density'); ax.set_title('Replogle 2022 K562 genome-wide: UMIs per cell')
    ax.legend(frameon=False); fig.tight_layout(); fig.savefig(OUT / 'umi_hist.png', dpi=150); plt.close(fig)

    # 2. lane medians
    m = lanes[lanes.n_control_cells > 0].sort_values('median').reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(np.arange(len(m)), m['median'], 'o', ms=3, color='#2a78d6')
    ax.axhline(whole['median'], ls='--', color='grey', lw=1)
    ax.set_xlabel(f'gem_group (10x lane), sorted; {len(m)} lanes'); ax.set_ylabel('median control log10 UMIs')
    ax.set_title(f'Depth by lane: control medians span {10 ** m["median"].min():,.0f}–{10 ** m["median"].max():,.0f} UMIs')
    fig.tight_layout(); fig.savefig(OUT / 'gemgroup_medians.png', dpi=150); plt.close(fig)

    # summary
    x = obs['x'].to_numpy(); lo_edge = np.sort(x)[:1000]
    counts, edges = np.histogram(x, bins=bins)
    first = int(np.argmax(counts > 0)); peak = int(counts.argmax())
    print(f'\ncells: {len(obs):,}; controls: {ctrl.sum():,}; lanes: {obs.gem_group.nunique()} ({int(lanes.fallback.sum())} use the whole-file cutoffs)')
    print(f'MAD x{scale} (3 MAD = 3 SD); whole-file cutoffs: low {whole["low_umi"]:,.0f} UMIs, high {whole["high_umi"]:,.0f} UMIs '
          f'(median {10 ** whole["median"] - 1:,.0f}, MAD {whole["mad"]:.3f} log10)')
    print(f'removed by the per-lane cutoffs: controls {ctrl_rate:.2%}, knockdown cells {obs.loc[kd, "removed"].mean():.2%} '
          f'(below low {int((removed & (umi < lim.low_umi.to_numpy())).sum()):,}, above high {int((removed & (umi > lim.high_umi.to_numpy())).sum()):,})')
    print(f'flag bar: share removed > {bar:.2%}; flagged knockdowns: {int(per.flagged.sum())} of {len(per):,} with >= {MIN_KD_CELLS} cells')
    print(per.head(10)[['n_cells', 'n_removed', 'share_removed', 'flagged']].to_string())
    print(f'\nleft edge: min UMI {obs.UMI_count.min():,.0f}; the first non-empty bin holds {counts[first]:,} cells, '
          f'{counts[first] / counts[peak]:.1%} of the peak bin; lowest 1,000 cells span {10 ** lo_edge.min():,.0f}–{10 ** lo_edge.max():,.0f} UMIs')
    print(f'wrote {OUT}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--file', default=str(FILE))
    main(Path(p.parse_args().file))

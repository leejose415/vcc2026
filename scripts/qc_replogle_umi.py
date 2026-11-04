"""Step 1 of cell QC for a single-cell perturbation source: per-cell UMI distribution and control-referenced cutoffs.

Sources (SOURCES): Replogle 2022 K562 / RPE1, Nadig 2025 HepG2 / Jurkat (dense counts; obs gene, gem_group, UMI_count,
sgID_AB) and Arc 2025 H1 (sparse counts; obs target_gene, batch, guide_id; no UMI column, so depth = the matrix row sum).
Reads metadata only (plus row sums for H1); never modifies the source file.

Cutoff, computed on control cells ('non-targeting') on x = log10(UMI_count + 1):
  low  = median - 3 MAD,  high = median + 4 MAD     MAD = median absolute deviation x 1.4826, so 3 MAD = 3 standard
                                                      deviations for a normal distribution; floored at 0.05
  low_umi = max(200, 10**low - 1)                    an absolute floor of 200 UMIs
(a) once for the whole file and (b) per lane (gem_group / batch), falling back to the whole-file controls when a lane has
fewer than LANE_MIN controls: 200 for K562; 100 for RPE1 (22 of 56 lanes hold 119-199 controls); 50 for the rest
(HepG2: every lane holds < 200, median 91, lane depth varies 5x; the standing rule for later sources).

Left edge: if the lowest histogram bin holds >= 1% of the peak bin, the file was pre-filtered by its authors (a sharp
left edge); step 2 then skips the low cut (floor and high cut only). Saved to left_edge.json.

Outputs (results/qc/<source>/): umi_hist.png, gemgroup_medians.png, cutoffs_by_gemgroup.csv, left_edge.json,
knockdowns_flagged.csv (knockdowns with >= 30 cells: share removed by the per-lane cutoffs, flagged when above
max(2 x the control rate, control rate + 5 points)).

Usage: ~/venv/bin/python scripts/qc_replogle_umi.py --dataset K562|RPE1|HepG2|Jurkat|H1
"""
import argparse
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata.io import read_elem

V = Path.home() / 'vcc'
# dataset -> (raw file, QC output folder, filtered file, title)
SOURCES = {
    'K562':   (V / 'lib/replogle/K562_gwps_raw_singlecell_01.h5ad', V / 'results/qc/replogle_k562_gw',
               V / 'lib/qc_cells/replogle_k562_gw_qc.h5ad', 'Replogle 2022 K562 genome-wide'),
    'RPE1':   (V / 'lib/downloads/rpe1/rpe1_raw_singlecell_01.h5ad', V / 'results/qc/replogle_rpe1',
               V / 'lib/qc_cells/replogle_rpe1_qc.h5ad', 'Replogle 2022 RPE1 essential'),
    'HepG2':  (V / 'lib/downloads/nadig_hepg2/GSE264667_hepg2_raw_singlecell_01.h5ad', V / 'results/qc/nadig_hepg2',
               V / 'lib/qc_cells/nadig_hepg2_qc.h5ad', 'Nadig 2025 HepG2'),
    'Jurkat': (V / 'lib/downloads/nadig_jurkat/GSE264667_jurkat_raw_singlecell_01.h5ad', V / 'results/qc/nadig_jurkat',
               V / 'lib/qc_cells/nadig_jurkat_qc.h5ad', 'Nadig 2025 Jurkat'),
    'H1':     (V / 'lib/h1/H1_full_2025.h5ad', V / 'results/qc/arc_h1',
               V / 'lib/qc_cells/arc_h1_qc.h5ad', 'Arc 2025 H1'),
    # X-Atlas/Orion: only the subsets we hold (the 300 A-C targets + 20,000 subsampled controls, on the 18,533 challenge genes)
    'HCT116': (V / 'lib/orion/HCT116_cells_targets.h5ad', V / 'results/qc/orion_hct116',
               V / 'lib/qc_cells/orion_hct116_qc.h5ad', 'X-Atlas/Orion 2025 HCT116 (A-C targets subset)'),
    'HEK293T': (V / 'lib/orion/HEK293T_cells_targets.h5ad', V / 'results/qc/orion_hek293t',
                V / 'lib/qc_cells/orion_hek293t_qc.h5ad', 'X-Atlas/Orion 2025 HEK293T (A-C targets subset)'),
}
BARCODE_PREFIX = '<barcode prefix>'
# source obs column for each role; umi None = the matrix row sum; gene_id None = look the target up in var['gene_id']
# (or, without one, in the H1 / K562 gene tables); guide None = no guide column released
LAYOUTS = {
    'replogle': dict(target='gene', lane='gem_group', umi='UMI_count', guide='sgID_AB', gene_id='gene_id'),
    'h1':       dict(target='target_gene', lane='batch', umi=None, guide='guide_id', gene_id=None),
    # lane = the batch in the barcode ('HCT116_Batch1.parquet:ACGT...' -> 'HCT116_Batch1'); no guide column was released
    'orion':    dict(target='target_gene', lane=BARCODE_PREFIX, umi=None, guide=None, gene_id=None),
}
LAYOUT_OF = {'H1': 'h1', 'HCT116': 'orion', 'HEK293T': 'orion'}      # everything else: 'replogle'
FILE, OUT = SOURCES['K562'][:2]
CONTROL = 'non-targeting'
LOW_K, HIGH_K, MAD_FLOOR, UMI_FLOOR, MIN_LANE_CONTROLS, MIN_KD_CELLS = 3, 4, 0.05, 200, 200, 30
MAD_SCALE = 1.4826                # MAD -> standard deviation for a normal distribution
LANE_MIN = {'RPE1': 100, 'HepG2': 50, 'Jurkat': 50, 'H1': 50, 'HCT116': 50, 'HEK293T': 50}   # minimum controls for a lane's own cutoffs; MIN_LANE_CONTROLS otherwise
SHARP_EDGE = 0.01                 # lowest bin >= 1% of the peak bin -> pre-filtered by the authors

# New CRISPRi contexts (2026-10-26), converted by scripts/convert_new_contexts.py into lib/cells/<name>.h5ad (single-target
# cells, integer counts, obs target_gene / guide_id / batch, var gene_id): the 'h1' layout, the standing 50-control lane
# minimum, and the true cell line as the background. The entries above are unchanged.
FENG_TARGETED_LINES = ('eipl_1', 'eipl_3', 'fiaj_1', 'fiaj_3', 'iudw_1', 'iudw_4', 'jejf_2', 'jejf_3', 'kolf_2', 'kolf_3',
                       'oikd_2', 'oikd_5', 'paab_3', 'paab_4', 'pipw_4', 'pipw_5', 'tolg_4', 'tolg_6', 'zapk_3')
NEW_CONTEXTS = {   # key: (file name, title, background)
    'KOLF2.1J': ('kolf_kolf21j', 'KOLF 2025 KOLF2.1J hiPSC CRISPRi (author-selected strong perturbations)', 'KOLF2.1J'),
    'Adamson_10X001': ('adamson_k562_10x001', 'Adamson 2016 K562 CRISPRi (10X001)', 'K562'),
    'Adamson_10X005': ('adamson_k562_10x005', 'Adamson 2016 K562 CRISPRi (10X005, single-gene constructs)', 'K562'),
    'Adamson_10X010': ('adamson_k562_10x010', 'Adamson 2016 K562 CRISPRi (10X010)', 'K562'),
    **{f'FengT_{l}': (f'feng_targeted_{l}', f'Feng HipSci iPSC CRISPRi targeted screen ({l})', l) for l in FENG_TARGETED_LINES},
    'FengGW_fitness': ('feng_gw_fitness', 'Feng HipSci iPSC CRISPRi genome-wide fitness genes (donor lines pooled)', 'HipSci_pooled'),
    'FengGW_nonfitness': ('feng_gw_nonfitness', 'Feng HipSci iPSC CRISPRi genome-wide non-fitness genes (donor lines pooled)', 'HipSci_pooled'),
}
# PIE-relevant extra sources (2026-10-26; their effect tables go to lib/effects_v5_extra/, not the scored-task tables):
# Tahoe-100M per cell line (DMSO_TF = the controls, lane = plate; drugs at their highest dose, cell QC only) and Jiang
# pathway Perturb-seq per cell line x pathway (CRISPRko; NT = the controls, lane = Parse sample well; standard on-target QC)
TAHOE_LINES = ('A-172', 'A-427', 'A498', 'A549', 'AN3 CA', 'AsPC-1', 'BT-474', 'C-33 A', 'C32', 'CFPAC-1', 'CHP-212', 'COLO 205',
               'H4', 'HCT15', 'HEC-1-A', 'HOP62', 'HS-578T', 'HT-29', 'HepG2/C3A', 'Hs 766T', 'J82', 'KATO III', 'LOX-IMVI', 'LS 180',
               'LoVo', 'MIA PaCa-2', 'NCI-H1573', 'NCI-H1792', 'NCI-H2030', 'NCI-H2122', 'NCI-H23', 'NCI-H2347', 'NCI-H460', 'NCI-H596',
               'NCI-H661', 'PANC-1', 'Panc 03.27', 'RKO', 'RPMI-7951', 'SHP-77', 'SK-MEL-2', 'SNU-1', 'SNU-423', 'SW 1088', 'SW 1271',
               'SW 900', 'SW1417', 'SW48', 'SW480', 'hTERT-HPNE')
JIANG_LINES, JIANG_PATHWAYS = ('A549', 'BXPC3', 'HAP1', 'HT29', 'K562', 'MCF7'), ('IFNB', 'IFNG', 'INS', 'TGFB', 'TNFA')


def _slug(s):
    return ''.join(c.lower() if c.isalnum() else '_' for c in str(s)).strip('_')


NEW_CONTEXTS.update({f'Tahoe_{_slug(l)}': (f'tahoe_{_slug(l)}', f'Tahoe-100M {l} (DMSO_TF + single-target drugs, highest dose)', l)
                     for l in TAHOE_LINES})
NEW_CONTEXTS.update({f'Jiang_{p}_{l}': (f'jiang_{p.lower()}_{l.lower()}', f'Jiang 2025 pathway Perturb-seq {p} {l} (CRISPRko)', l)
                     for p in JIANG_PATHWAYS for l in JIANG_LINES})
for _key, (_name, _title, _bg) in NEW_CONTEXTS.items():
    SOURCES[_key] = (V / f'lib/cells/{_name}.h5ad', V / f'results/qc/{_name}', V / f'lib/qc_cells/{_name}_qc.h5ad', _title)
    LAYOUT_OF[_key] = 'h1'; LANE_MIN[_key] = 50
BACKGROUND = {k: bg for k, (_, _, bg) in NEW_CONTEXTS.items()}    # write_replogle_qc's background; the dataset key otherwise
# Reference cells for the UMI cutoffs, the flag rate and the on-target test: the non-targeting controls, except for the
# Feng genome-wide screens (~0-2 non-targeting cells per donor line), where every other assigned cell stands in for them
# (rest of the lane for the UMI cutoffs; rest of the screen for the on-target test). User decision, 2026-10-25/26.
REFERENCE = {'FengGW_fitness': 'rest', 'FengGW_nonfitness': 'rest'}


def reference_mask(dataset, gene):
    """The cells the UMI cutoffs and the flag rate are referenced to: controls, or every cell under REFERENCE 'rest'."""
    gene = pd.Series(gene)
    return np.ones(len(gene), bool) if REFERENCE.get(dataset) == 'rest' else (gene == CONTROL).to_numpy()


def layout(dataset):
    return LAYOUTS[LAYOUT_OF.get(dataset, 'replogle')]


def row_sums(path, cache=None):
    """Per-cell sum of the stored counts (CSR or dense), read in chunks; cached next to the QC outputs."""
    if cache is not None and Path(cache).exists():
        stored = np.load(cache, allow_pickle=False)
        if int(stored['mtime_ns']) == Path(path).stat().st_mtime_ns: return stored['sums']
    with h5py.File(path, 'r') as f:
        X = f['X']
        if isinstance(X, h5py.Group) and X.attrs.get('encoding-type') == 'csr_matrix':
            ip = X['indptr'][:]; n = len(ip) - 1; out = np.zeros(n)
            for s in range(0, n, 20000):
                e = min(s + 20000, n); a, b = int(ip[s]), int(ip[e])
                out[s:e] = np.bincount(np.repeat(np.arange(e - s), np.diff(ip[s:e + 1])), weights=X['data'][a:b].astype(float), minlength=e - s)
        else:
            n = X.shape[0]; out = np.zeros(n)
            for s in range(0, n, 20000): out[s:s + 20000] = X[s:s + 20000].sum(1)
    if cache is not None:
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        np.savez(cache, sums=out, mtime_ns=np.int64(Path(path).stat().st_mtime_ns))
    return out


def load_cells(dataset):
    """The source's obs plus the QC columns under fixed names: gene, gem_group, UMI_count (and umi_source)."""
    path, out, _, _ = SOURCES[dataset]; lay = layout(dataset)
    with h5py.File(path, 'r') as f: obs = read_elem(f['obs'])
    cells = pd.DataFrame(index=obs.index)
    cells['gene'] = obs[lay['target']].astype(str).to_numpy()
    if lay['lane'] == BARCODE_PREFIX:
        cells['gem_group'] = obs.index.astype(str).str.split(':').str[0].str.replace('.parquet', '', regex=False).to_numpy()
    else:
        cells['gem_group'] = obs[lay['lane']].astype(str).to_numpy()
    if lay['umi']: cells['UMI_count'] = obs[lay['umi']].to_numpy(float); cells.attrs['umi_source'] = f"obs['{lay['umi']}']"
    else: cells['UMI_count'] = row_sums(path, out / 'row_sums.npz'); cells.attrs['umi_source'] = 'matrix row sum (the source has no UMI column)'
    return obs, cells


def cutoffs(x, scale=MAD_SCALE):
    """Median, MAD (floored) and the low / high cutoffs in UMIs for log10 values x."""
    med = float(np.median(x)); mad = max(float(np.median(np.abs(x - med))) * scale, MAD_FLOOR)
    return dict(median=med, mad=mad, low_umi=max(UMI_FLOOR, 10 ** (med - LOW_K * mad) - 1), high_umi=10 ** (med + HIGH_K * mad) - 1)


def main(dataset='K562', scale=MAD_SCALE):
    path, OUT, _, title = SOURCES[dataset]; min_lane = LANE_MIN.get(dataset, MIN_LANE_CONTROLS)
    _, obs = load_cells(dataset)
    print('gene value counts (head):'); print(obs['gene'].value_counts().head().to_string())
    if CONTROL not in set(obs['gene']): raise ValueError(f'no {CONTROL!r} cells: check the control label')
    obs['x'] = np.log10(obs['UMI_count'].to_numpy(float) + 1); ctrl = obs['gene'] == CONTROL
    ref = pd.Series(reference_mask(dataset, obs['gene']), index=obs.index)          # == ctrl unless REFERENCE is 'rest'
    whole = cutoffs(obs.loc[ref, 'x'].to_numpy(), scale)

    rows = []
    for lane, g in obs[ref].groupby('gem_group', sort=False):
        c = cutoffs(g['x'].to_numpy(), scale) if len(g) >= min_lane else dict(whole)
        rows.append(dict(gem_group=lane, n_control_cells=len(g), **c, fallback=len(g) < min_lane))
    for lane in set(obs['gem_group']) - {r['gem_group'] for r in rows}:          # a lane with no controls at all
        rows.append(dict(gem_group=lane, n_control_cells=0, **whole, fallback=True))
    lanes = pd.DataFrame(rows).sort_values('gem_group', key=lambda s: pd.to_numeric(s, errors='coerce')).reset_index(drop=True)
    OUT.mkdir(parents=True, exist_ok=True)
    lanes[['gem_group', 'n_control_cells', 'median', 'mad', 'low_umi', 'high_umi']].to_csv(OUT / 'cutoffs_by_gemgroup.csv', index=False)

    lim = obs[['gem_group']].merge(lanes[['gem_group', 'low_umi', 'high_umi']], on='gem_group', how='left')
    umi = obs['UMI_count'].to_numpy(float)
    removed = (umi < lim['low_umi'].to_numpy()) | (umi > lim['high_umi'].to_numpy())
    obs['removed'] = removed
    ctrl_rate = float(obs.loc[ref, 'removed'].mean()); kd = ~ctrl
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
    ax.set_xlabel('log10 total UMIs per cell'); ax.set_ylabel('density'); ax.set_title(f'{title}: UMIs per cell')
    ax.legend(frameon=False); fig.tight_layout(); fig.savefig(OUT / 'umi_hist.png', dpi=150); plt.close(fig)

    # 2. lane medians
    m = lanes[lanes.n_control_cells > 0].sort_values('median').reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(np.arange(len(m)), m['median'], 'o', ms=3, color='#2a78d6')
    ax.axhline(whole['median'], ls='--', color='grey', lw=1)
    ax.set_xlabel(f'lane (gem group / batch), sorted; {len(m)} lanes'); ax.set_ylabel('median control log10 UMIs')
    ax.set_title(f'Depth by lane: control medians span {10 ** m["median"].min():,.0f}–{10 ** m["median"].max():,.0f} UMIs')
    fig.tight_layout(); fig.savefig(OUT / 'gemgroup_medians.png', dpi=150); plt.close(fig)

    # left edge: pre-filtered by the authors?
    x = obs['x'].to_numpy(); lo_edge = np.sort(x)[:1000]
    counts, _ = np.histogram(x, bins=bins)
    first = int(np.argmax(counts > 0)); peak = int(counts.argmax()); share = float(counts[first] / counts[peak])
    edge = dict(min_umi=float(obs.UMI_count.min()), first_bin_cells=int(counts[first]), first_bin_share_of_peak=share,
                sharp=share >= SHARP_EDGE, low_cut='skipped (pre-filtered by the authors)' if share >= SHARP_EDGE else 'applied',
                rule=f'sharp when the lowest of 80 bins holds >= {SHARP_EDGE:.0%} of the peak bin', umi_source=obs.attrs['umi_source'])
    (OUT / 'left_edge.json').write_text(json.dumps(edge, indent=1) + '\n')

    print(f'\ncells: {len(obs):,}; controls: {ctrl.sum():,}; lanes: {obs.gem_group.nunique()} ({int(lanes.fallback.sum())} with < {min_lane} controls use the whole-file cutoffs)')
    print(f'UMIs from {obs.attrs["umi_source"]}')
    print(f'MAD x{scale} (3 MAD = 3 SD); whole-file cutoffs: low {whole["low_umi"]:,.0f} UMIs, high {whole["high_umi"]:,.0f} UMIs '
          f'(median {10 ** whole["median"] - 1:,.0f}, MAD {whole["mad"]:.3f} log10)')
    print(f'first pass (per-lane low and high): controls {ctrl_rate:.2%}, knockdown cells {obs.loc[kd, "removed"].mean():.2%} '
          f'(below low {int((removed & (umi < lim.low_umi.to_numpy())).sum()):,}, above high {int((removed & (umi > lim.high_umi.to_numpy())).sum()):,})')
    print(f'flag bar: share removed > {bar:.2%}; flagged knockdowns: {int(per.flagged.sum())} of {len(per):,} with >= {MIN_KD_CELLS} cells')
    print(per.head(10)[['n_cells', 'n_removed', 'share_removed', 'flagged']].to_string())
    print(f'\nleft edge: min UMI {obs.UMI_count.min():,.0f}; the lowest bin holds {counts[first]:,} cells, {share:.1%} of the peak bin; '
          f'lowest 1,000 cells span {10 ** lo_edge.min():,.0f}–{10 ** lo_edge.max():,.0f} UMIs -> {"SHARP: pre-filtered, low cut skipped" if edge["sharp"] else "no sharp edge, low cut applied"}')
    print(f'wrote {OUT}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--dataset', default='K562', choices=list(SOURCES))
    main(p.parse_args().dataset)

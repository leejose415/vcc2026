"""Step 3 of cell QC: write the filtered file (Replogle 2022 K562 / RPE1, Nadig 2025 HepG2 / Jurkat, Arc 2025 H1): the cells
in kept_cells.txt, counts unchanged, with standardised obs columns. Never modifies the source file.

Paths per dataset: qc_replogle_umi.SOURCES (raw file -> lib/qc_cells/<source>_qc.h5ad; QC outputs in results/qc/<source>/).
Run after scripts/qc_replogle_umi.py and scripts/qc_replogle_apply.py.

0. Checks on the source: every count is a whole number (per cell); where the source has a UMI column, it is the raw depth:
   at least 99.9% of cells have row sum <= UMI and >= 0.9 x UMI (file-level; the released matrix omits lowly expressed
   genes, so it holds ~98-99% of a cell's UMIs). Cells outside are kept and listed in checks.json.
1. Output: kept cells in source order; X as int32 CSR (gzip); var unchanged; every original obs column kept, plus
     source, background, condition ('none'), batch (= the lane), target_gene, is_control, target_gene_id (Ensembl ID;
     'non-targeting' for controls), guide_id (the guide-pair column: sgID_AB, one-to-one with gene_transcript; H1:
     guide_id), n_umi (the UMIs used for QC), qc_low_ref ('control', 'own', or 'none' when the low cut was skipped)
   Columns the source already has with the same meaning (H1: target_gene, guide_id, batch) are kept as they are.
   uns['qc']: definitions, rules, lane minimum, left edge, per-lane and own cutoffs, counts, source checks, the source's
   size and modification time.
2. Checks after writing: cells == lines in kept_cells.txt and the same barcodes; 100 random kept cells' counts identical to
   the source; the source's size and modification time unchanged. Saved to results/qc/<source>/checks.json and qc_method.md.

--metadata-only: for a file already written, replace obs and uns with the standardised ones (counts are not rewritten),
then run the checks in 2.

Usage: ~/venv/bin/python scripts/write_replogle_qc.py --dataset K562|RPE1|HepG2|Jurkat|H1 [--metadata-only]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem, write_elem

sys.path.insert(0, str(Path(__file__).parent))
from build_raw_cells import CSRWriter                                                     # noqa: E402
from qc_replogle_umi import (CONTROL, HIGH_K, LANE_MIN, LOW_K, MAD_FLOOR, MAD_SCALE, MIN_KD_CELLS, MIN_LANE_CONTROLS,  # noqa: E402
                             SOURCES, UMI_FLOOR, cutoffs, layout, load_cells)

CHUNK = 20000
N_CHECK = 100
UMI_PASS_SHARE = 0.999


def source_stat(path):
    s = Path(path).stat(); return dict(size=s.st_size, mtime=time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(s.st_mtime)), mtime_ns=s.st_mtime_ns)


def read_rows(X, start, stop, n_vars):
    """Rows start:stop of a dense or CSR source matrix, as CSR."""
    if isinstance(X, h5py.Group):
        ip = X['indptr'][start:stop + 1]; a, b = int(ip[0]), int(ip[-1])
        return sp.csr_matrix((X['data'][a:b], X['indices'][a:b], ip - a), shape=(stop - start, n_vars))
    return sp.csr_matrix(X[start:stop])


def gene_id_lookup():
    """Symbol -> Ensembl ID from the H1 and K562 gene tables, for sources whose var has no gene_id (Orion subsets)."""
    lookup = {}
    for path, sym in ((SOURCES['K562'][0], 'gene_name'), (SOURCES['H1'][0], None)):
        if not Path(path).exists(): continue
        with h5py.File(path, 'r') as f: v = read_elem(f['var'])
        names = v[sym].astype(str) if sym else pd.Series(v.index.astype(str), index=v.index)
        ids = pd.Series(v.index.astype(str), index=v.index) if sym else v['gene_id'].astype(str)
        lookup.update({n: i for n, i in zip(names, ids) if n not in lookup})
    return lookup


def standard_obs(obs, cells, var, dataset, source, flagged, skip_low):
    """The original columns plus the standard ones (same names in every source)."""
    lay = layout(dataset); gene = cells['gene']; ctrl = gene == CONTROL
    if lay['gene_id']: gid = obs[lay['gene_id']].astype(str)
    else:
        lookup = dict(zip(var.index.astype(str), var['gene_id'].astype(str))) if 'gene_id' in var else gene_id_lookup()
        gid = gene.map(lookup).fillna('')
    std = pd.DataFrame(index=obs.index)
    std['source'] = source; std['background'] = dataset; std['condition'] = 'none'
    std['batch'] = cells['gem_group']; std['target_gene'] = gene; std['is_control'] = ctrl.to_numpy()
    std['target_gene_id'] = np.where(ctrl, CONTROL, gid)
    std['guide_id'] = obs[lay['guide']].astype(str) if lay['guide'] else 'not_available'; std['n_umi'] = cells['UMI_count'].to_numpy()
    std['qc_low_ref'] = 'none' if skip_low else np.where(gene.isin(flagged), 'own', 'control')
    for c in list(std.columns):
        if c in obs.columns:
            if (obs[c].astype(str).to_numpy() == std[c].astype(str).to_numpy()).all(): std = std.drop(columns=c)   # same meaning
            else: raise ValueError(f'obs already has {c!r} with different values')
    for c in ('source', 'background', 'condition', 'batch', 'target_gene', 'target_gene_id', 'guide_id', 'qc_low_ref'):
        if c in std: std[c] = std[c].astype('category')
    return pd.concat([obs, std], axis=1)


def qc_record(dataset, cells, keep, out, src_stat):
    """uns['qc']: the method, its parameters, the cutoffs used and the counts."""
    edge = json.loads((out / 'left_edge.json').read_text()); skip_low = bool(edge['sharp'])
    lanes = pd.read_csv(out / 'cutoffs_by_gemgroup.csv', dtype={'gem_group': str})
    flagged_table = pd.read_csv(out / 'knockdowns_flagged.csv')
    flagged = [] if skip_low else flagged_table.loc[flagged_table.flagged, 'knockdown'].astype(str).tolist()
    removed = pd.read_csv(out / 'removed_cells.csv')
    x = np.log10(cells['UMI_count'].to_numpy(float) + 1); gene = cells['gene'].to_numpy()
    own = {k: cutoffs(x[gene == k])['low_umi'] for k in flagged}
    ctrl = gene == CONTROL; whole = cutoffs(x[ctrl]); lane_min = int(LANE_MIN.get(dataset, MIN_LANE_CONTROLS))
    rules = [f'every cell: removed if UMI < {UMI_FLOOR} (floor) or above its lane high cutoff (high)']
    rules += (['low cut skipped: the source was pre-filtered by its authors (sharp left edge, see left_edge)'] if skip_low else
              ['controls and unflagged knockdowns: also removed below their lane low cutoff (low_control)',
               "flagged knockdowns: removed below the knockdown's own low cutoff instead (low_own)",
               f'flagged: knockdowns with >= {MIN_KD_CELLS} cells whose first-pass removal share exceeds max(2 x control rate, control rate + 0.05)'])
    rules.append(f'reference cells per lane: its controls; fewer than {lane_min} -> all controls in the file')
    return dict(
        definitions=dict(x='log10(UMI + 1)', umi_source=cells.attrs['umi_source'],
                         scaled_mad=f'max({MAD_SCALE} * median(|x - median(x)|), {MAD_FLOOR})',
                         low=f'median - {LOW_K} * scaled_mad', high=f'median + {HIGH_K} * scaled_mad',
                         note='scaled MAD everywhere (1 scaled MAD = 1 SD for a normal distribution, as scater isOutlier); never the raw MAD'),
        rules=rules, lane_min_controls=lane_min, left_edge=edge,
        whole_file_cutoffs=dict(low_umi=float(whole['low_umi']), high_umi=float(whole['high_umi']), median_umi=float(10 ** whole['median'] - 1),
                                scaled_mad=float(whole['mad'])),
        lane_cutoffs={c: lanes[c].tolist() for c in ('gem_group', 'n_control_cells', 'low_umi', 'high_umi')},
        own_low_cutoffs=dict(knockdown=list(own), low_umi=[float(v) for v in own.values()]),
        counts=dict(cells_source=int(len(cells)), cells_kept=int(keep.sum()), controls=int(ctrl.sum()),
                    controls_removed=float(1 - keep[ctrl].mean()), knockdown_cells_removed=float(1 - keep[~ctrl].mean()),
                    flagged_knockdowns=int(len(flagged)), knockdowns_with_30_cells=int(len(flagged_table)),
                    removed_by_reason={str(k): int(v) for k, v in removed['reason'].value_counts().items()}),
        source=dict(file=str(SOURCES[dataset][0]), **{k: v for k, v in src_stat.items() if k != 'mtime_ns'}),
        method_note=str(out / 'qc_method.md')), set(flagged), skip_low


def check_source(f, cells, n_vars, has_umi_column):
    """Every count a whole number; where a UMI column exists, >= 99.9% of cells have 0.9 x UMI <= row sum <= UMI."""
    X = f['X']; umi = cells['UMI_count'].to_numpy(float); bad_int = 0; sums = np.zeros(len(cells))
    for start in range(0, len(cells), CHUNK):
        stop = min(start + CHUNK, len(cells)); block = read_rows(X, start, stop, n_vars)
        bad_int += int(np.unique(block.nonzero()[0][block.data != np.round(block.data)]).size)
        sums[start:stop] = np.asarray(block.sum(1)).ravel()
    result = dict(cells=len(cells), non_integer_cells=bad_int)
    if has_umi_column:
        ratio = sums / umi; ok = (sums <= umi + 1) & (sums >= 0.9 * umi)
        bad = np.flatnonzero(~ok)
        result.update(umi_check='file-level: >= 99.9% of cells with 0.9 x UMI <= row sum <= UMI', cells_passing=float(ok.mean()),
                      row_sum_over_umi_median=float(np.median(ratio)), row_sum_over_umi_min=float(ratio.min()),
                      exceptions=[dict(cell=str(cells.index[i]), gene=str(cells['gene'].iloc[i]), umi=float(umi[i]), row_sum=float(sums[i]),
                                       ratio=round(float(ratio[i]), 4)) for i in bad[:200]], exceptions_total=int(len(bad)))
    else:
        result.update(umi_check='skipped: no UMI column; depth is the matrix row sum')
    print('source checks:', {k: v for k, v in result.items() if k != 'exceptions'}, flush=True)
    if bad_int: raise ValueError(f'source has {bad_int} cells with non-integer counts (normalised data?)')
    if has_umi_column and result['cells_passing'] < UMI_PASS_SHARE: raise ValueError(f'UMI column check failed: {result}')
    return result


def post_checks(dest, src_path, out, src_before, source_checks):
    kept = pd.read_csv(out / 'kept_cells.txt', header=None)[0].astype(str).to_numpy()
    with h5py.File(dest, 'r') as g, h5py.File(src_path, 'r') as f:
        names = read_elem(g['obs']).index.astype(str).to_numpy()
        src_names = pd.Index(read_elem(f['obs']).index.astype(str))
        ip = g['X']['indptr']; n_vars = int(g['X'].attrs['shape'][1])
        pick = np.sort(np.random.default_rng(2026).choice(len(names), min(N_CHECK, len(names)), replace=False))
        same = 0
        for i in pick:
            a, b = int(ip[i]), int(ip[i + 1])
            row = np.zeros(n_vars, np.int64); row[g['X']['indices'][a:b]] = g['X']['data'][a:b]
            j = src_names.get_loc(names[i])
            same += int(np.array_equal(row, read_rows(f['X'], j, j + 1, n_vars).toarray().ravel().astype(np.int64)))
    after = source_stat(src_path)
    checks = dict(cells=int(len(names)), kept_lines=int(len(kept)), count_matches=bool(len(names) == len(kept)),
                  barcodes_match=bool(np.array_equal(names, kept)), random_cells_checked=int(len(pick)), random_cells_identical=int(same),
                  source_size_unchanged=bool(after['size'] == src_before['size']), source_mtime_unchanged=bool(after['mtime_ns'] == src_before['mtime_ns']),
                  source_size=after['size'], source_mtime=after['mtime'])
    if source_checks is not None: checks['source_checks'] = source_checks
    (out / 'checks.json').write_text(json.dumps(checks, indent=1) + '\n')
    sc = source_checks or {}
    umi_line = (f"UMI column check (file-level): {sc['cells_passing']:.4%} of cells pass; {sc['exceptions_total']} exception(s), kept and listed in checks.json."
                if 'cells_passing' in sc else (f"UMI column check: {sc['umi_check']}." if sc else 'Source checks: run when the file was first written.'))
    md = (out / 'qc_method.md').read_text().split('\n## Filtered file')[0].rstrip()
    (out / 'qc_method.md').write_text(md + f"""

## Filtered file
`{dest}`: {checks['cells']:,} cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
the standard ones (`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id`, `n_umi`, `qc_low_ref`);
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: {checks['count_matches']}; same barcodes in order: {checks['barcodes_match']};
{checks['random_cells_identical']} of {checks['random_cells_checked']} random kept cells identical to the source;
source size unchanged: {checks['source_size_unchanged']}, modification time unchanged: {checks['source_mtime_unchanged']}
({checks['source_size']:,} bytes, {checks['source_mtime']} UTC). {umi_line}
""")
    print('checks:', {k: v for k, v in checks.items() if k != 'source_checks'}, flush=True)
    if not (checks['count_matches'] and checks['barcodes_match'] and same == len(pick)
            and checks['source_size_unchanged'] and checks['source_mtime_unchanged']):
        raise ValueError(f'post-write checks failed: {checks}')
    return checks


def main(dataset, metadata_only=False):
    src_path, out, dest, _ = SOURCES[dataset]
    source = out.name; lay = layout(dataset)
    src_before = source_stat(src_path)
    kept = pd.Index(pd.read_csv(out / 'kept_cells.txt', header=None)[0].astype(str))
    obs, cells = load_cells(dataset)
    t0 = time.time(); source_checks = None
    with h5py.File(src_path, 'r') as f:
        var = read_elem(f['var']); n_vars = len(var)
        keep = obs.index.astype(str).isin(kept)
        if keep.sum() != len(kept): raise ValueError(f'{len(kept) - keep.sum()} kept barcodes are not in the source file')
        record, flagged, skip_low = qc_record(dataset, cells, keep, out, src_before)
        new_obs = standard_obs(obs[keep], cells[keep], var, dataset, source, flagged, skip_low)
        if metadata_only:
            if not dest.exists(): raise FileNotFoundError(dest)
        else:
            source_checks = check_source(f, cells, n_vars, bool(lay['umi']))
            # uns stores the exceptions column-wise (h5ad cannot hold a list of records); checks.json keeps them as records
            exc = source_checks.get('exceptions', [])
            record['source_checks'] = {k: v for k, v in source_checks.items() if k != 'exceptions'}
            if exc: record['source_checks']['exceptions'] = {k: [e[k] for e in exc] for k in exc[0]}
            dest.parent.mkdir(parents=True, exist_ok=True)
            w = CSRWriter(dest, n_vars)
            for start in range(0, len(obs), CHUNK):
                stop = min(start + CHUNK, len(obs)); m = keep[start:stop]
                if m.any(): w.append(read_rows(f['X'], start, stop, n_vars)[m])
                if (start // CHUNK) % 10 == 0: print(f'  {stop:,} / {len(obs):,} source cells read, {w.rows:,} written ({time.time() - t0:.0f}s)', flush=True)
            w.close(new_obs, var, dict(qc=record))
    if metadata_only:
        with h5py.File(dest, 'r+') as g:
            old = read_elem(g['obs'])
            if not np.array_equal(old.index.astype(str), new_obs.index.astype(str)): raise ValueError('existing file holds different cells')
            for key in ('obs', 'uns'): del g[key]
            write_elem(g, 'obs', new_obs); write_elem(g, 'uns', dict(qc=record))
        print(f'replaced obs and uns of {dest}', flush=True)
    else:
        print(f'wrote {dest}: {int(keep.sum()):,} of {len(obs):,} cells, {dest.stat().st_size / 1e9:.1f} GB ({(time.time() - t0) / 60:.0f} min)', flush=True)
    post_checks(dest, src_path, out, src_before, source_checks)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--dataset', required=True, choices=list(SOURCES))
    p.add_argument('--metadata-only', action='store_true', help='replace obs / uns of an already written file; counts untouched')
    a = p.parse_args(); main(a.dataset, a.metadata_only)

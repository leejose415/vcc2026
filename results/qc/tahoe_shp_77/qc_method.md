# Cell QC: Tahoe-100M SHP-77 (DMSO_TF + single-target drugs, highest dose)

Source: `/home/mooneysad/vcc/lib/cells/tahoe_shp_77.h5ad` (read only; cell metadata only; UMIs from matrix row sum (the source has no UMI column)). Scripts: `scripts/qc_replogle_umi.py --dataset Tahoe_shp_77`
(distribution, cutoffs, flagged knockdowns, left edge), `scripts/qc_replogle_apply.py --dataset Tahoe_shp_77` (this filter),
`scripts/write_replogle_qc.py --dataset Tahoe_shp_77` (the filtered file).

## Definitions
On x = log10(UMI + 1), over the reference cells:
- median = median(x)
- scaled MAD = 1.4826 × median(|x − median|), floor 0.05 (1 scaled MAD = 1 SD for a normal distribution, as scater's `isOutlier`); the raw MAD is never used
- low cutoff = median − 3 × scaled MAD; high cutoff = median + 4 × scaled MAD

Reference cells: per lane (10x gem group / batch), that lane's control cells (`non-targeting`); a lane with fewer than **50** controls
uses all controls in the file (0 of 14 lanes). For each flagged knockdown
(`knockdowns_flagged.csv`: knockdowns with ≥ 30 cells losing more than max(2 × the control rate, control rate + 5 points) to
the lane cutoffs), that knockdown's own cells, for the low cutoff only.

## Rules
1. Every cell: removed if UMI < 200 (`floor`, empty droplets) or above its lane's high cutoff (`high`, likely doublets).
2. (skipped: pre-filtered source) Controls and unflagged knockdowns below their lane's low cutoff.
3. (skipped: pre-filtered source) Flagged knockdowns below their own low cutoff.

## Pre-filtered source: no low cut
The UMI histogram has a sharp left edge (`left_edge.json`): the lowest cell has 700 UMIs and the lowest bin still holds
45.5% of the peak bin, so the authors filtered at about 700 UMIs before release. Only the
200-UMI floor and the per-lane high cut are applied (rules 2 and 3 skipped). Low-UMI cells of strong knockdowns were removed
before we received the file, so some strong effects (loss of core transcription machinery and similar) may look weaker here than in K562,
where those cells were kept by rule 3.

## Whole-file cutoffs (all controls)
low 383 UMIs, high 8,159 UMIs (median 1,422, scaled MAD 0.190 log10).
Lane cutoffs vary with depth: control medians span 974–1,746 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: 73,166; kept 72,892 (99.63%); removed 274
- removed by reason: floor 0, high 274, low_control 0, low_own 0
- controls removed: 0.14% of 17,039; knockdown cells removed: 0.45% of 56,127

Top 10 knockdowns flagged in step 1: share of cells removed in the first pass (lane low and high cuts) and by the final rules:

| knockdown | cells | first pass | final | own low cutoff (UMIs) |
|---|---|---|---|---|

## Filtered file
`/home/mooneysad/vcc/lib/qc_cells/tahoe_shp_77_qc.h5ad`: 72,892 cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
the standard ones (`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id`, `n_umi`, `qc_low_ref`);
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: True; same barcodes in order: True;
100 of 100 random kept cells identical to the source;
source size unchanged: True, modification time unchanged: True
(206,198,901 bytes, 2026-10-26 04:47:49 UTC). UMI column check: skipped: no UMI column; depth is the matrix row sum.

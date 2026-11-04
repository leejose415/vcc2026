# Cell QC: Tahoe-100M C-33 A (DMSO_TF + single-target drugs, highest dose)

Source: `/home/mooneysad/vcc/lib/cells/tahoe_c_33_a.h5ad` (read only; cell metadata only; UMIs from matrix row sum (the source has no UMI column)). Scripts: `scripts/qc_replogle_umi.py --dataset Tahoe_c_33_a`
(distribution, cutoffs, flagged knockdowns, left edge), `scripts/qc_replogle_apply.py --dataset Tahoe_c_33_a` (this filter),
`scripts/write_replogle_qc.py --dataset Tahoe_c_33_a` (the filtered file).

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
49.7% of the peak bin, so the authors filtered at about 700 UMIs before release. Only the
200-UMI floor and the per-lane high cut are applied (rules 2 and 3 skipped). Low-UMI cells of strong knockdowns were removed
before we received the file, so some strong effects (loss of core transcription machinery and similar) may look weaker here than in K562,
where those cells were kept by rule 3.

## Whole-file cutoffs (all controls)
low 351 UMIs, high 9,605 UMIs (median 1,450, scaled MAD 0.205 log10).
Lane cutoffs vary with depth: control medians span 993–2,024 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: 73,030; kept 72,504 (99.28%); removed 526
- removed by reason: floor 0, high 526, low_control 0, low_own 0
- controls removed: 0.24% of 16,902; knockdown cells removed: 0.87% of 56,128

Top 10 knockdowns flagged in step 1: share of cells removed in the first pass (lane low and high cuts) and by the final rules:

| knockdown | cells | first pass | final | own low cutoff (UMIs) |
|---|---|---|---|---|

## Filtered file
`/home/mooneysad/vcc/lib/qc_cells/tahoe_c_33_a_qc.h5ad`: 72,504 cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
the standard ones (`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id`, `n_umi`, `qc_low_ref`);
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: True; same barcodes in order: True;
100 of 100 random kept cells identical to the source;
source size unchanged: True, modification time unchanged: True
(199,249,951 bytes, 2026-10-26 04:30:08 UTC). UMI column check: skipped: no UMI column; depth is the matrix row sum.

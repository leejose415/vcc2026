# Cell QC: Adamson 2016 K562 CRISPRi (10X001)

Source: `/home/mooneysad/vcc/lib/cells/adamson_k562_10x001.h5ad` (read only; cell metadata only; UMIs from matrix row sum (the source has no UMI column)). Scripts: `scripts/qc_replogle_umi.py --dataset Adamson_10X001`
(distribution, cutoffs, flagged knockdowns, left edge), `scripts/qc_replogle_apply.py --dataset Adamson_10X001` (this filter),
`scripts/write_replogle_qc.py --dataset Adamson_10X001` (the filtered file).

## Definitions
On x = log10(UMI + 1), over the reference cells:
- median = median(x)
- scaled MAD = 1.4826 × median(|x − median|), floor 0.05 (1 scaled MAD = 1 SD for a normal distribution, as scater's `isOutlier`); the raw MAD is never used
- low cutoff = median − 3 × scaled MAD; high cutoff = median + 4 × scaled MAD

Reference cells: per lane (10x gem group / batch), that lane's control cells (`non-targeting`); a lane with fewer than **50** controls
uses all controls in the file (0 of 1 lanes). For each flagged knockdown
(`knockdowns_flagged.csv`: knockdowns with ≥ 30 cells losing more than max(2 × the control rate, control rate + 5 points) to
the lane cutoffs), that knockdown's own cells, for the low cutoff only.

## Rules
1. Every cell: removed if UMI < 200 (`floor`, empty droplets) or above its lane's high cutoff (`high`, likely doublets).
2. (skipped: pre-filtered source) Controls and unflagged knockdowns below their lane's low cutoff.
3. (skipped: pre-filtered source) Flagged knockdowns below their own low cutoff.

## Pre-filtered source: no low cut
The UMI histogram has a sharp left edge (`left_edge.json`): the lowest cell has 1,975 UMIs and the lowest bin still holds
2.3% of the peak bin, so the authors filtered at about 2,000 UMIs before release. Only the
200-UMI floor and the per-lane high cut are applied (rules 2 and 3 skipped). Low-UMI cells of strong knockdowns were removed
before we received the file, so some strong effects (loss of core transcription machinery and similar) may look weaker here than in K562,
where those cells were kept by rule 3.

## Whole-file cutoffs (all controls)
low 4,307 UMIs, high 30,817 UMIs (median 10,010, scaled MAD 0.122 log10).
Lane cutoffs vary with depth: control medians span 10,010–10,010 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: 5,752; kept 5,752 (100.00%); removed 0
- removed by reason: floor 0, high 0, low_control 0, low_own 0
- controls removed: 0.00% of 1,769; knockdown cells removed: 0.00% of 3,983

Top 10 knockdowns flagged in step 1: share of cells removed in the first pass (lane low and high cuts) and by the final rules:

| knockdown | cells | first pass | final | own low cutoff (UMIs) |
|---|---|---|---|---|

## Filtered file
`/home/mooneysad/vcc/lib/qc_cells/adamson_k562_10x001_qc.h5ad`: 5,752 cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
the standard ones (`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id`, `n_umi`, `qc_low_ref`);
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: True; same barcodes in order: True;
100 of 100 random kept cells identical to the source;
source size unchanged: True, modification time unchanged: True
(36,638,131 bytes, 2026-10-26 01:24:41 UTC). UMI column check: skipped: no UMI column; depth is the matrix row sum.

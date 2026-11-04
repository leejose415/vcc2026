# Cell QC: Feng HipSci iPSC CRISPRi genome-wide fitness genes (donor lines pooled)

Source: `/home/mooneysad/vcc/lib/cells/feng_gw_fitness.h5ad` (read only; cell metadata only; UMIs from matrix row sum (the source has no UMI column)). Scripts: `scripts/qc_replogle_umi.py --dataset FengGW_fitness`
(distribution, cutoffs, flagged knockdowns, left edge), `scripts/qc_replogle_apply.py --dataset FengGW_fitness` (this filter),
`scripts/write_replogle_qc.py --dataset FengGW_fitness` (the filtered file).

## Definitions
On x = log10(UMI + 1), over the reference cells:
- median = median(x)
- scaled MAD = 1.4826 × median(|x − median|), floor 0.05 (1 scaled MAD = 1 SD for a normal distribution, as scater's `isOutlier`); the raw MAD is never used
- low cutoff = median − 3 × scaled MAD; high cutoff = median + 4 × scaled MAD

Reference cells: per lane (10x gem group / batch), that lane's control cells (`non-targeting`); a lane with fewer than **50** controls
uses all controls in the file (0 of 24 lanes). For each flagged knockdown
(`knockdowns_flagged.csv`: knockdowns with ≥ 30 cells losing more than max(2 × the control rate, control rate + 5 points) to
the lane cutoffs), that knockdown's own cells, for the low cutoff only.

**Rest-of-lane reference** (`REFERENCE` = rest): this source has almost no non-targeting cells, so every cell of the lane (all assigned cells) stands in for the controls, for the cutoffs and for the flag rate.

## Rules
1. Every cell: removed if UMI < 200 (`floor`, empty droplets) or above its lane's high cutoff (`high`, likely doublets).
2. Controls and unflagged knockdowns: also removed below their lane's low cutoff (`low_control`).
3. Flagged knockdowns: removed below the knockdown's own low cutoff instead (`low_own`). Their low UMI counts are biology
   (in K562, loss of core transcription machinery: Mediator, SUPT5H / SUPT6H, TAF1, GATA1), not poor quality.
   Flagged here: ALG2, SEC31B, MAP1S, RPP30, POU5F1B, MRPL52, HPS5.

## Whole-file cutoffs (all controls)
low 9,165 UMIs, high 91,918 UMIs (median 24,618, scaled MAD 0.143 log10).
Lane cutoffs vary with depth: control medians span 22,022–31,948 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: 104,130; kept 103,601 (99.49%); removed 529
- removed by reason: floor 0, high 81, low_control 441, low_own 7
- controls removed: 0.00% of 36; knockdown cells removed: 0.51% of 104,094

Top 10 knockdowns flagged in step 1: share of cells removed in the first pass (lane low and high cuts) and by the final rules:

| knockdown | cells | first pass | final | own low cutoff (UMIs) |
|---|---|---|---|---|
| ALG2 | 38 | 7.9% | 0.0% | 6,050 |
| SEC31B | 42 | 7.1% | 9.5% | 12,461 |
| MAP1S | 42 | 7.1% | 2.4% | 6,798 |
| RPP30 | 31 | 6.5% | 0.0% | 5,615 |
| POU5F1B | 34 | 5.9% | 11.8% | 14,005 |
| MRPL52 | 35 | 5.7% | 2.9% | 5,897 |
| HPS5 | 54 | 5.6% | 0.0% | 6,879 |

## Filtered file
`/home/mooneysad/vcc/lib/qc_cells/feng_gw_fitness_qc.h5ad`: 103,601 cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
the standard ones (`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id`, `n_umi`, `qc_low_ref`);
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: True; same barcodes in order: True;
100 of 100 random kept cells identical to the source;
source size unchanged: True, modification time unchanged: True
(1,332,372,220 bytes, 2026-10-26 02:42:54 UTC). UMI column check: skipped: no UMI column; depth is the matrix row sum.

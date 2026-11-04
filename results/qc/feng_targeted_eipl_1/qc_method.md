# Cell QC: Feng HipSci iPSC CRISPRi targeted screen (eipl_1)

Source: `/home/mooneysad/vcc/lib/cells/feng_targeted_eipl_1.h5ad` (read only; cell metadata only; UMIs from matrix row sum (the source has no UMI column)). Scripts: `scripts/qc_replogle_umi.py --dataset FengT_eipl_1`
(distribution, cutoffs, flagged knockdowns, left edge), `scripts/qc_replogle_apply.py --dataset FengT_eipl_1` (this filter),
`scripts/write_replogle_qc.py --dataset FengT_eipl_1` (the filtered file).

## Definitions
On x = log10(UMI + 1), over the reference cells:
- median = median(x)
- scaled MAD = 1.4826 × median(|x − median|), floor 0.05 (1 scaled MAD = 1 SD for a normal distribution, as scater's `isOutlier`); the raw MAD is never used
- low cutoff = median − 3 × scaled MAD; high cutoff = median + 4 × scaled MAD

Reference cells: per lane (10x gem group / batch), that lane's control cells (`non-targeting`); a lane with fewer than **50** controls
uses all controls in the file (72 of 72 lanes). For each flagged knockdown
(`knockdowns_flagged.csv`: knockdowns with ≥ 30 cells losing more than max(2 × the control rate, control rate + 5 points) to
the lane cutoffs), that knockdown's own cells, for the low cutoff only.

## Rules
1. Every cell: removed if UMI < 200 (`floor`, empty droplets) or above its lane's high cutoff (`high`, likely doublets).
2. Controls and unflagged knockdowns: also removed below their lane's low cutoff (`low_control`).
3. Flagged knockdowns: removed below the knockdown's own low cutoff instead (`low_own`). Their low UMI counts are biology
   (in K562, loss of core transcription machinery: Mediator, SUPT5H / SUPT6H, TAF1, GATA1), not poor quality.
   Flagged here: FA2H, CCND1, FBXO46, ZNF257, CHAF1A.

## Whole-file cutoffs (all controls)
low 6,565 UMIs, high 74,034 UMIs (median 18,543, scaled MAD 0.150 log10).
Lane cutoffs vary with depth: control medians span 0–0 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: 38,699; kept 37,954 (98.07%); removed 745
- removed by reason: floor 0, high 1, low_control 736, low_own 8
- controls removed: 1.61% of 496; knockdown cells removed: 1.93% of 38,203

Top 10 knockdowns flagged in step 1: share of cells removed in the first pass (lane low and high cuts) and by the final rules:

| knockdown | cells | first pass | final | own low cutoff (UMIs) |
|---|---|---|---|---|
| FA2H | 90 | 12.2% | 0.0% | 3,383 |
| CCND1 | 44 | 9.1% | 0.0% | 5,231 |
| FBXO46 | 39 | 7.7% | 2.6% | 5,938 |
| ZNF257 | 68 | 7.4% | 5.9% | 5,258 |
| CHAF1A | 43 | 7.0% | 7.0% | 6,552 |

## Filtered file
`/home/mooneysad/vcc/lib/qc_cells/feng_targeted_eipl_1_qc.h5ad`: 37,954 cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
the standard ones (`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id`, `n_umi`, `qc_low_ref`);
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: True; same barcodes in order: True;
100 of 100 random kept cells identical to the source;
source size unchanged: True, modification time unchanged: True
(429,475,814 bytes, 2026-10-26 02:42:56 UTC). UMI column check: skipped: no UMI column; depth is the matrix row sum.

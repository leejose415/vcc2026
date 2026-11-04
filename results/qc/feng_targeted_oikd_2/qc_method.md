# Cell QC: Feng HipSci iPSC CRISPRi targeted screen (oikd_2)

Source: `/home/mooneysad/vcc/lib/cells/feng_targeted_oikd_2.h5ad` (read only; cell metadata only; UMIs from matrix row sum (the source has no UMI column)). Scripts: `scripts/qc_replogle_umi.py --dataset FengT_oikd_2`
(distribution, cutoffs, flagged knockdowns, left edge), `scripts/qc_replogle_apply.py --dataset FengT_oikd_2` (this filter),
`scripts/write_replogle_qc.py --dataset FengT_oikd_2` (the filtered file).

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
   Flagged here: MRPL18, TARDBP.

## Whole-file cutoffs (all controls)
low 5,461 UMIs, high 92,953 UMIs (median 18,402, scaled MAD 0.176 log10).
Lane cutoffs vary with depth: control medians span 0–0 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: 14,127; kept 14,008 (99.16%); removed 119
- removed by reason: floor 0, high 0, low_control 116, low_own 3
- controls removed: 0.00% of 177; knockdown cells removed: 0.85% of 13,950

Top 10 knockdowns flagged in step 1: share of cells removed in the first pass (lane low and high cuts) and by the final rules:

| knockdown | cells | first pass | final | own low cutoff (UMIs) |
|---|---|---|---|---|
| MRPL18 | 38 | 5.3% | 0.0% | 4,134 |
| TARDBP | 57 | 5.3% | 5.3% | 5,419 |

## Filtered file
`/home/mooneysad/vcc/lib/qc_cells/feng_targeted_oikd_2_qc.h5ad`: 14,008 cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
the standard ones (`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id`, `n_umi`, `qc_low_ref`);
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: True; same barcodes in order: True;
100 of 100 random kept cells identical to the source;
source size unchanged: True, modification time unchanged: True
(160,470,818 bytes, 2026-10-26 02:43:03 UTC). UMI column check: skipped: no UMI column; depth is the matrix row sum.

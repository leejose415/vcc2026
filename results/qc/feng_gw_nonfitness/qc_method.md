# Cell QC: Feng HipSci iPSC CRISPRi genome-wide non-fitness genes (donor lines pooled)

Source: `/home/mooneysad/vcc/lib/cells/feng_gw_nonfitness.h5ad` (read only; cell metadata only; UMIs from matrix row sum (the source has no UMI column)). Scripts: `scripts/qc_replogle_umi.py --dataset FengGW_nonfitness`
(distribution, cutoffs, flagged knockdowns, left edge), `scripts/qc_replogle_apply.py --dataset FengGW_nonfitness` (this filter),
`scripts/write_replogle_qc.py --dataset FengGW_nonfitness` (the filtered file).

## Definitions
On x = log10(UMI + 1), over the reference cells:
- median = median(x)
- scaled MAD = 1.4826 × median(|x − median|), floor 0.05 (1 scaled MAD = 1 SD for a normal distribution, as scater's `isOutlier`); the raw MAD is never used
- low cutoff = median − 3 × scaled MAD; high cutoff = median + 4 × scaled MAD

Reference cells: per lane (10x gem group / batch), that lane's control cells (`non-targeting`); a lane with fewer than **50** controls
uses all controls in the file (0 of 34 lanes). For each flagged knockdown
(`knockdowns_flagged.csv`: knockdowns with ≥ 30 cells losing more than max(2 × the control rate, control rate + 5 points) to
the lane cutoffs), that knockdown's own cells, for the low cutoff only.

**Rest-of-lane reference** (`REFERENCE` = rest): this source has almost no non-targeting cells, so every cell of the lane (all assigned cells) stands in for the controls, for the cutoffs and for the flag rate.

## Rules
1. Every cell: removed if UMI < 200 (`floor`, empty droplets) or above its lane's high cutoff (`high`, likely doublets).
2. Controls and unflagged knockdowns: also removed below their lane's low cutoff (`low_control`).
3. Flagged knockdowns: removed below the knockdown's own low cutoff instead (`low_own`). Their low UMI counts are biology
   (in K562, loss of core transcription machinery: Mediator, SUPT5H / SUPT6H, TAF1, GATA1), not poor quality.
   Flagged here: TUBB6.

## Whole-file cutoffs (all controls)
low 3,930 UMIs, high 286,757 UMIs (median 24,712, scaled MAD 0.266 log10).
Lane cutoffs vary with depth: control medians span 12,551–59,596 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: 115,124; kept 114,989 (99.88%); removed 135
- removed by reason: floor 0, high 4, low_control 131, low_own 0
- controls removed: 0.00% of 12; knockdown cells removed: 0.12% of 115,112

Top 10 knockdowns flagged in step 1: share of cells removed in the first pass (lane low and high cuts) and by the final rules:

| knockdown | cells | first pass | final | own low cutoff (UMIs) |
|---|---|---|---|---|
| TUBB6 | 30 | 6.7% | 0.0% | 3,956 |

## Filtered file
`/home/mooneysad/vcc/lib/qc_cells/feng_gw_nonfitness_qc.h5ad`: 114,989 cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
the standard ones (`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id`, `n_umi`, `qc_low_ref`);
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: True; same barcodes in order: True;
100 of 100 random kept cells identical to the source;
source size unchanged: True, modification time unchanged: True
(1,522,280,302 bytes, 2026-10-26 02:42:55 UTC). UMI column check: skipped: no UMI column; depth is the matrix row sum.

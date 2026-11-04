# Cell QC: Feng HipSci iPSC CRISPRi targeted screen (iudw_1)

Source: `/home/mooneysad/vcc/lib/cells/feng_targeted_iudw_1.h5ad` (read only; cell metadata only; UMIs from matrix row sum (the source has no UMI column)). Scripts: `scripts/qc_replogle_umi.py --dataset FengT_iudw_1`
(distribution, cutoffs, flagged knockdowns, left edge), `scripts/qc_replogle_apply.py --dataset FengT_iudw_1` (this filter),
`scripts/write_replogle_qc.py --dataset FengT_iudw_1` (the filtered file).

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
   Flagged here: TAF8, CLN3, WDR1.

## Whole-file cutoffs (all controls)
low 6,242 UMIs, high 102,372 UMIs (median 20,701, scaled MAD 0.174 log10).
Lane cutoffs vary with depth: control medians span 0–0 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: 44,469; kept 43,693 (98.25%); removed 776
- removed by reason: floor 0, high 1, low_control 771, low_own 4
- controls removed: 1.22% of 575; knockdown cells removed: 1.75% of 43,894

Top 10 knockdowns flagged in step 1: share of cells removed in the first pass (lane low and high cuts) and by the final rules:

| knockdown | cells | first pass | final | own low cutoff (UMIs) |
|---|---|---|---|---|
| TAF8 | 68 | 10.3% | 0.0% | 4,164 |
| CLN3 | 42 | 7.1% | 2.4% | 5,515 |
| WDR1 | 87 | 6.9% | 3.4% | 5,594 |

## Filtered file
`/home/mooneysad/vcc/lib/qc_cells/feng_targeted_iudw_1_qc.h5ad`: 43,693 cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
the standard ones (`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id`, `n_umi`, `qc_low_ref`);
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: True; same barcodes in order: True;
100 of 100 random kept cells identical to the source;
source size unchanged: True, modification time unchanged: True
(496,649,613 bytes, 2026-10-26 02:42:59 UTC). UMI column check: skipped: no UMI column; depth is the matrix row sum.

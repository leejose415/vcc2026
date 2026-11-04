# Cell QC: Feng HipSci iPSC CRISPRi targeted screen (jejf_3)

Source: `/home/mooneysad/vcc/lib/cells/feng_targeted_jejf_3.h5ad` (read only; cell metadata only; UMIs from matrix row sum (the source has no UMI column)). Scripts: `scripts/qc_replogle_umi.py --dataset FengT_jejf_3`
(distribution, cutoffs, flagged knockdowns, left edge), `scripts/qc_replogle_apply.py --dataset FengT_jejf_3` (this filter),
`scripts/write_replogle_qc.py --dataset FengT_jejf_3` (the filtered file).

## Definitions
On x = log10(UMI + 1), over the reference cells:
- median = median(x)
- scaled MAD = 1.4826 × median(|x − median|), floor 0.05 (1 scaled MAD = 1 SD for a normal distribution, as scater's `isOutlier`); the raw MAD is never used
- low cutoff = median − 3 × scaled MAD; high cutoff = median + 4 × scaled MAD

Reference cells: per lane (10x gem group / batch), that lane's control cells (`non-targeting`); a lane with fewer than **50** controls
uses all controls in the file (54 of 54 lanes). For each flagged knockdown
(`knockdowns_flagged.csv`: knockdowns with ≥ 30 cells losing more than max(2 × the control rate, control rate + 5 points) to
the lane cutoffs), that knockdown's own cells, for the low cutoff only.

## Rules
1. Every cell: removed if UMI < 200 (`floor`, empty droplets) or above its lane's high cutoff (`high`, likely doublets).
2. Controls and unflagged knockdowns: also removed below their lane's low cutoff (`low_control`).
3. Flagged knockdowns: removed below the knockdown's own low cutoff instead (`low_own`). Their low UMI counts are biology
   (in K562, loss of core transcription machinery: Mediator, SUPT5H / SUPT6H, TAF1, GATA1), not poor quality.
   Flagged here: LAMTOR5, EIF3B.

## Whole-file cutoffs (all controls)
low 6,778 UMIs, high 106,718 UMIs (median 22,089, scaled MAD 0.171 log10).
Lane cutoffs vary with depth: control medians span 0–0 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: 34,851; kept 34,246 (98.26%); removed 605
- removed by reason: floor 0, high 0, low_control 603, low_own 2
- controls removed: 2.07% of 484; knockdown cells removed: 1.73% of 34,367

Top 10 knockdowns flagged in step 1: share of cells removed in the first pass (lane low and high cuts) and by the final rules:

| knockdown | cells | first pass | final | own low cutoff (UMIs) |
|---|---|---|---|---|
| LAMTOR5 | 73 | 9.6% | 1.4% | 4,920 |
| EIF3B | 63 | 7.9% | 1.6% | 4,525 |

## Filtered file
`/home/mooneysad/vcc/lib/qc_cells/feng_targeted_jejf_3_qc.h5ad`: 34,246 cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
the standard ones (`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id`, `n_umi`, `qc_low_ref`);
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: True; same barcodes in order: True;
100 of 100 random kept cells identical to the source;
source size unchanged: True, modification time unchanged: True
(413,441,108 bytes, 2026-10-26 02:43:01 UTC). UMI column check: skipped: no UMI column; depth is the matrix row sum.

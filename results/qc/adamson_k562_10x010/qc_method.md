# Cell QC: Adamson 2016 K562 CRISPRi (10X010)

Source: `/home/mooneysad/vcc/lib/cells/adamson_k562_10x010.h5ad` (read only; cell metadata only; UMIs from matrix row sum (the source has no UMI column)). Scripts: `scripts/qc_replogle_umi.py --dataset Adamson_10X010`
(distribution, cutoffs, flagged knockdowns, left edge), `scripts/qc_replogle_apply.py --dataset Adamson_10X010` (this filter),
`scripts/write_replogle_qc.py --dataset Adamson_10X010` (the filtered file).

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
2. Controls and unflagged knockdowns: also removed below their lane's low cutoff (`low_control`).
3. Flagged knockdowns: removed below the knockdown's own low cutoff instead (`low_own`). Their low UMI counts are biology
   (in K562, loss of core transcription machinery: Mediator, SUPT5H / SUPT6H, TAF1, GATA1), not poor quality.
   Flagged here: HARS, XBP1, CARS, QARS, EIF2S1, SARS.

## Whole-file cutoffs (all controls)
low 5,877 UMIs, high 64,543 UMIs (median 16,413, scaled MAD 0.149 log10).
Lane cutoffs vary with depth: control medians span 16,413–16,413 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: 62,623; kept 60,119 (96.00%); removed 2,504
- removed by reason: floor 0, high 1, low_control 2,501, low_own 2
- controls removed: 3.63% of 7,295; knockdown cells removed: 4.05% of 55,328

Top 10 knockdowns flagged in step 1: share of cells removed in the first pass (lane low and high cuts) and by the final rules:

| knockdown | cells | first pass | final | own low cutoff (UMIs) |
|---|---|---|---|---|
| HARS | 325 | 11.4% | 0.3% | 3,573 |
| XBP1 | 62 | 11.3% | 0.0% | 2,825 |
| CARS | 501 | 11.2% | 0.0% | 2,757 |
| QARS | 379 | 9.8% | 0.3% | 3,314 |
| EIF2S1 | 511 | 9.2% | 0.0% | 3,208 |
| SARS | 473 | 8.9% | 0.0% | 3,240 |

## Filtered file
`/home/mooneysad/vcc/lib/qc_cells/adamson_k562_10x010_qc.h5ad`: 60,119 cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
the standard ones (`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id`, `n_umi`, `qc_low_ref`);
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: True; same barcodes in order: True;
100 of 100 random kept cells identical to the source;
source size unchanged: True, modification time unchanged: True
(503,646,081 bytes, 2026-10-26 01:26:25 UTC). UMI column check: skipped: no UMI column; depth is the matrix row sum.

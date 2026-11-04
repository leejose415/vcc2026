# Cell QC: Replogle 2022 RPE1 essential

Source: `/home/mooneysad/vcc/lib/downloads/rpe1/rpe1_raw_singlecell_01.h5ad` (read only; obs columns `UMI_count`, `gene`, `gem_group`). Scripts: `scripts/qc_replogle_umi.py` (distribution,
cutoffs, flagged knockdowns), `scripts/qc_replogle_apply.py --dataset RPE1` (this filter). No .h5ad is written; apply `kept_cells.txt` when loading.

## Definitions
On x = log10(UMI_count + 1), over the reference cells:
- median = median(x)
- scaled MAD = 1.4826 × median(|x − median|), floor 0.05 (1 scaled MAD = 1 SD for a normal distribution, as scater's `isOutlier`); the raw MAD is never used
- low cutoff = median − 3 × scaled MAD; high cutoff = median + 4 × scaled MAD

Reference cells: per gem group (10x lane), that lane's control cells (`non-targeting`); a lane with fewer than 100 controls
uses all controls in the file (0 of 56 lanes). For each of the 128 flagged knockdowns
(`knockdowns_flagged.csv`: knockdowns with ≥ 30 cells losing more than max(2 × the control rate, control rate + 5 points) to
the lane cutoffs), that knockdown's own cells, for the low cutoff only.

## Rules
1. Every cell: removed if UMI_count < 200 (`floor`, empty droplets) or above its lane's high cutoff (`high`, likely doublets).
2. Controls and unflagged knockdowns: also removed below their lane's low cutoff (`low_control`).
3. Flagged knockdowns: removed below the knockdown's own low cutoff instead (`low_own`). Their low UMI counts are biology
   (in K562, loss of core transcription machinery: Mediator, SUPT5H / SUPT6H, TAF1, GATA1), not poor quality.
   Flagged here: EIF2S1, GPKOW, PABPC1, SART1, CTR9, SF3B1, FAM102B, MMS22L, PABPN1, BUD13 ....

## Whole-file cutoffs (all controls)
low 4,127 UMIs, high 48,731 UMIs (median 11,890, scaled MAD 0.153 log10).
Lane cutoffs vary with depth: control medians span 5,937–34,727 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: 247,914; kept 243,500 (98.22%); removed 4,414
- removed by reason: floor 0, high 375, low_control 3,964, low_own 75
- controls removed: 1.53% of 11,485; knockdown cells removed: 1.79% of 236,429

Top 10 flagged knockdowns, share of cells removed with the lane low cutoff (before rule 3) and with their own (after):

| knockdown | cells | before | after | own low cutoff (UMIs) |
|---|---|---|---|---|
| EIF2S1 | 83 | 24.1% | 0.0% | 2,184 |
| GPKOW | 42 | 23.8% | 0.0% | 1,219 |
| PABPC1 | 102 | 23.5% | 0.0% | 1,675 |
| SART1 | 105 | 22.9% | 0.0% | 861 |
| CTR9 | 41 | 22.0% | 0.0% | 1,293 |
| SF3B1 | 43 | 18.6% | 0.0% | 1,207 |
| FAM102B | 81 | 18.5% | 0.0% | 1,405 |
| MMS22L | 33 | 18.2% | 18.2% | 3,521 |
| PABPN1 | 41 | 17.1% | 0.0% | 2,066 |
| BUD13 | 36 | 16.7% | 0.0% | 2,123 |

## Filtered file
`/home/mooneysad/vcc/lib/qc_cells/replogle_rpe1_qc.h5ad`: 243,500 cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id` (= `sgID_AB`), `n_umi`, `qc_low_ref`;
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: True; same barcodes in order: True;
100 of 100 random kept cells identical to the source;
source size unchanged: True, modification time unchanged: True
(8,700,873,216 bytes, 2026-10-18 13:19:23 UTC).

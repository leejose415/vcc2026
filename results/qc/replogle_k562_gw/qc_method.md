# Cell QC: Replogle 2022 K562 genome-wide

Source: `/home/mooneysad/vcc/lib/replogle/K562_gwps_raw_singlecell_01.h5ad` (read only; obs columns `UMI_count`, `gene`, `gem_group`). Scripts: `scripts/qc_replogle_umi.py` (distribution,
cutoffs, flagged knockdowns), `scripts/qc_replogle_apply.py` (this filter). No .h5ad is written; apply `kept_cells.txt` when loading.

## Definitions
On x = log10(UMI_count + 1), over the reference cells:
- median = median(x)
- scaled MAD = 1.4826 × median(|x − median|), floor 0.05 (1 scaled MAD = 1 SD for a normal distribution, as scater's `isOutlier`); the raw MAD is never used
- low cutoff = median − 3 × scaled MAD; high cutoff = median + 4 × scaled MAD

Reference cells: per gem group (10x lane), that lane's control cells (`non-targeting`); a lane with fewer than 200 controls
uses all controls in the file (8 of 267 lanes). For each of the 358 flagged knockdowns
(`knockdowns_flagged.csv`: knockdowns with ≥ 30 cells losing more than max(2 × the control rate, control rate + 5 points) to
the lane cutoffs), that knockdown's own cells, for the low cutoff only.

## Rules
1. Every cell: removed if UMI_count < 200 (`floor`, empty droplets) or above its lane's high cutoff (`high`, likely doublets).
2. Controls and unflagged knockdowns: also removed below their lane's low cutoff (`low_control`).
3. Flagged knockdowns: removed below the knockdown's own low cutoff instead (`low_own`). Their low UMI counts are biology
   (loss of core transcription machinery: Mediator, SUPT5H / SUPT6H, TAF1, GATA1), not poor quality.

## Whole-file cutoffs (all controls)
low 3,221 UMIs, high 63,546 UMIs (median 11,563, scaled MAD 0.185 log10).
Lane cutoffs vary with depth: control medians span about 2,900–20,800 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: 1,989,578; kept 1,936,807 (97.35%); removed 52,771
- removed by reason: floor 0, high 930, low_control 51,237, low_own 604
- controls removed: 2.65% of 75,328; knockdown cells removed: 2.65% of 1,914,250

Top 10 flagged knockdowns, share of cells removed with the lane low cutoff (before rule 3) and with their own (after):

| knockdown | cells | before | after | own low cutoff (UMIs) |
|---|---|---|---|---|
| MED30 | 114 | 71.1% | 0.0% | 592 |
| MED12 | 186 | 64.5% | 0.5% | 919 |
| SUPT6H | 105 | 60.0% | 0.0% | 438 |
| GATA1 | 46 | 54.3% | 0.0% | 200 |
| MED19 | 299 | 50.8% | 0.3% | 769 |
| SUPT5H | 139 | 50.4% | 0.0% | 324 |
| MED9 | 39 | 48.7% | 0.0% | 918 |
| MED7 | 122 | 42.6% | 0.0% | 445 |
| TAF1 | 88 | 38.6% | 0.0% | 698 |
| MED17 | 80 | 36.2% | 0.0% | 431 |

## Filtered file
`/home/mooneysad/vcc/lib/qc_cells/replogle_k562_gw_qc.h5ad`: 1,936,807 cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id` (= `sgID_AB`), `n_umi`, `qc_low_ref`;
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: True; same barcodes in order: True;
100 of 100 random kept cells identical to the source;
source size unchanged: True, modification time unchanged: True
(65,830,941,948 bytes, 2026-10-18 11:30:50 UTC).

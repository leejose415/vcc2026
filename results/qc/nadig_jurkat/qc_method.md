# Cell QC: Nadig 2025 Jurkat

Source: `/home/mooneysad/vcc/lib/downloads/nadig_jurkat/GSE264667_jurkat_raw_singlecell_01.h5ad` (read only; cell metadata only; UMIs from obs['UMI_count']). Scripts: `scripts/qc_replogle_umi.py --dataset Jurkat`
(distribution, cutoffs, flagged knockdowns, left edge), `scripts/qc_replogle_apply.py --dataset Jurkat` (this filter),
`scripts/write_replogle_qc.py --dataset Jurkat` (the filtered file).

## Definitions
On x = log10(UMI + 1), over the reference cells:
- median = median(x)
- scaled MAD = 1.4826 × median(|x − median|), floor 0.05 (1 scaled MAD = 1 SD for a normal distribution, as scater's `isOutlier`); the raw MAD is never used
- low cutoff = median − 3 × scaled MAD; high cutoff = median + 4 × scaled MAD

Reference cells: per lane (10x gem group / batch), that lane's control cells (`non-targeting`); a lane with fewer than **50** controls
uses all controls in the file (1 of 55 lanes). For each flagged knockdown
(`knockdowns_flagged.csv`: knockdowns with ≥ 30 cells losing more than max(2 × the control rate, control rate + 5 points) to
the lane cutoffs), that knockdown's own cells, for the low cutoff only.

## Rules
1. Every cell: removed if UMI < 200 (`floor`, empty droplets) or above its lane's high cutoff (`high`, likely doublets).
2. (skipped: pre-filtered source) Controls and unflagged knockdowns below their lane's low cutoff.
3. (skipped: pre-filtered source) Flagged knockdowns below their own low cutoff.

## Pre-filtered source: no low cut
The UMI histogram has a sharp left edge (`left_edge.json`): the lowest cell has 1,751 UMIs and the lowest bin still holds
1.6% of the peak bin, so the authors filtered at about 1,800 UMIs before release. Only the
200-UMI floor and the per-lane high cut are applied (rules 2 and 3 skipped). Low-UMI cells of strong knockdowns were removed
before we received the file, so some strong effects (loss of core transcription machinery and similar) may look weaker here than in K562,
where those cells were kept by rule 3.

## Whole-file cutoffs (all controls)
low 3,013 UMIs, high 60,591 UMIs (median 10,906, scaled MAD 0.186 log10).
Lane cutoffs vary with depth: control medians span 7,753–32,884 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: 262,956; kept 262,215 (99.72%); removed 741
- removed by reason: floor 0, high 741, low_control 0, low_own 0
- controls removed: 0.15% of 12,013; knockdown cells removed: 0.29% of 250,943

Top 10 knockdowns flagged in step 1: share of cells removed in the first pass (lane low and high cuts) and by the final rules:

| knockdown | cells | first pass | final | own low cutoff (UMIs) |
|---|---|---|---|---|
| TAF12 | 91 | 54.9% | 0.0% | — |
| MYC | 120 | 50.8% | 0.0% | — |
| MED27 | 35 | 48.6% | 0.0% | — |
| TAF5 | 80 | 47.5% | 0.0% | — |
| MED30 | 97 | 47.4% | 0.0% | — |
| MED11 | 95 | 47.4% | 1.1% | — |
| SUPT5H | 56 | 46.4% | 0.0% | — |
| SMC1A | 59 | 45.8% | 0.0% | — |
| MED7 | 47 | 44.7% | 0.0% | — |
| TAF8 | 174 | 43.7% | 0.0% | — |

## Filtered file
`/home/mooneysad/vcc/lib/qc_cells/nadig_jurkat_qc.h5ad`: 262,215 cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
the standard ones (`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id`, `n_umi`, `qc_low_ref`);
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: True; same barcodes in order: True;
100 of 100 random kept cells identical to the source;
source size unchanged: True, modification time unchanged: True
(9,366,490,264 bytes, 2026-10-18 11:49:50 UTC). UMI column check (file-level): 99.9951% of cells pass; 13 exception(s), kept and listed in checks.json.

# Cell QC: Nadig 2025 HepG2

Source: `/home/mooneysad/vcc/lib/downloads/nadig_hepg2/GSE264667_hepg2_raw_singlecell_01.h5ad` (read only; cell metadata only; UMIs from obs['UMI_count']). Scripts: `scripts/qc_replogle_umi.py --dataset HepG2`
(distribution, cutoffs, flagged knockdowns, left edge), `scripts/qc_replogle_apply.py --dataset HepG2` (this filter),
`scripts/write_replogle_qc.py --dataset HepG2` (the filtered file).

## Definitions
On x = log10(UMI + 1), over the reference cells:
- median = median(x)
- scaled MAD = 1.4826 × median(|x − median|), floor 0.05 (1 scaled MAD = 1 SD for a normal distribution, as scater's `isOutlier`); the raw MAD is never used
- low cutoff = median − 3 × scaled MAD; high cutoff = median + 4 × scaled MAD

Reference cells: per lane (10x gem group / batch), that lane's control cells (`non-targeting`); a lane with fewer than **50** controls
uses all controls in the file (3 of 56 lanes). For each flagged knockdown
(`knockdowns_flagged.csv`: knockdowns with ≥ 30 cells losing more than max(2 × the control rate, control rate + 5 points) to
the lane cutoffs), that knockdown's own cells, for the low cutoff only.

## Rules
1. Every cell: removed if UMI < 200 (`floor`, empty droplets) or above its lane's high cutoff (`high`, likely doublets).
2. (skipped: pre-filtered source) Controls and unflagged knockdowns below their lane's low cutoff.
3. (skipped: pre-filtered source) Flagged knockdowns below their own low cutoff.

## Pre-filtered source: no low cut
The UMI histogram has a sharp left edge (`left_edge.json`): the lowest cell has 3,001 UMIs and the lowest bin still holds
2.1% of the peak bin, so the authors filtered at about 3,000 UMIs before release. Only the
200-UMI floor and the per-lane high cut are applied (rules 2 and 3 skipped). Low-UMI cells of strong knockdowns were removed
before we received the file, so some strong effects (loss of core transcription machinery and similar) may look weaker here than in K562,
where those cells were kept by rule 3.

## Whole-file cutoffs (all controls)
low 3,606 UMIs, high 116,965 UMIs (median 16,019, scaled MAD 0.216 log10).
Lane cutoffs vary with depth: control medians span 8,294–40,026 UMIs across lanes (`cutoffs_by_gemgroup.csv`).

## Result
- cells: 145,473; kept 145,469 (100.00%); removed 4
- removed by reason: floor 0, high 4, low_control 0, low_own 0
- controls removed: 0.00% of 4,976; knockdown cells removed: 0.00% of 140,497

Top 10 knockdowns flagged in step 1: share of cells removed in the first pass (lane low and high cuts) and by the final rules:

| knockdown | cells | first pass | final | own low cutoff (UMIs) |
|---|---|---|---|---|
| CWC22 | 34 | 14.7% | 0.0% | — |
| RPL38 | 38 | 13.2% | 0.0% | — |
| COPS3 | 40 | 12.5% | 0.0% | — |
| TVP23C | 33 | 12.1% | 0.0% | — |
| CTBP2 | 45 | 11.1% | 0.0% | — |
| RPAP2 | 55 | 10.9% | 0.0% | — |
| PABPC1 | 66 | 10.6% | 0.0% | — |
| PSMB4 | 70 | 10.0% | 0.0% | — |
| DHX8 | 61 | 9.8% | 0.0% | — |
| RUVBL1 | 53 | 9.4% | 0.0% | — |

## Filtered file
`/home/mooneysad/vcc/lib/qc_cells/nadig_hepg2_qc.h5ad`: 145,469 cells, counts as int32 CSR (values identical to the source). obs keeps every original column and adds
the standard ones (`source`, `background`, `condition`, `batch`, `target_gene`, `is_control`, `target_gene_id`, `guide_id`, `n_umi`, `qc_low_ref`);
`uns['qc']` holds the definitions, rules, cutoffs and counts.

Checks: cells = lines in kept_cells.txt: True; same barcodes in order: True;
100 of 100 random kept cells identical to the source;
source size unchanged: True, modification time unchanged: True
(5,614,460,941 bytes, 2026-10-18 11:48:51 UTC). UMI column check (file-level): 99.9993% of cells pass; 1 exception(s), kept and listed in checks.json.

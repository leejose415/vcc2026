# Effect tables (v5): 117 sources in the scorer's regime
Built by `scripts/build_effects_v5.py`; Orion and the new contexts by `build_effects_v5_orion.py`, `build_effects_v5_new.py` and `build_effects_extra.py`. One file per source: `lib/effects_v5/<source>.h5ad` (`lib/effects_v5_extra/{tahoe,jiang}/` for drugs and knockouts), plus `<source>_baseline.csv`, `<source>_control_quantiles.npz` and `<source>_build.json`.

## Definitions: the same quantities `cell-eval2` scores (vcc2026 profile)
- **CPM:** per cell, counts / (library / 10⁶). The library is the cell's sum over the source's genes that map to the 18,533-gene panel.
- **LFC:** log2((mean CPM of the knockdown's cells + 1e-9) / (mean CPM of controls + 1e-9)), arithmetic means.
- **Test:** two-sided Mann–Whitney U on log1p(CPM), normal approximation, tie- and continuity-corrected.
  - **Gate:** a gene is tested only when its control mean CPM > 5.
  - **BH** per knockdown over the gated genes; significant = p_adj < 0.05.
  - The target's own gene is in the BH universe, then removed from the significant set.
- **Profile:** log1p(5·10⁴ × summed counts / summed total) per group. That is the scorer's pseudobulk, used by MSE and PDS.

## The scorer's regime
The scored contexts have 18,400 controls at a median of about 20k UMI, and 400 cells per knockdown. The tables are computed in that regime so their significance matches what the scorer can detect:
1. **Depth:** every cell is binomially thinned by one factor, min(1, 20,000 / the median control UMI), with a seed per cell. Sources below 20k are left as they are (no upsampling).
2. **Cells:** a random 400 cells per kept knockdown and 18,400 controls (all when fewer), seed 2026.

LFC, SE, the gate, the test and BH are computed on these sampled, thinned cells.

## Layers and annotations
- **Rows:** kept knockdowns.
- **`obs`:** `source`, `background` (cell line), `condition`, `target_gene`, `target_gene_id`, `n_cells_kept`, `n_cells_used`, `efficiency` (the target's own LFC), `n_real` (true DE count), `modality` (`single_cell`, `replicate_table`, `drug`), `significance_unknown`.
- **`layers`:**
  - `lfc`;
  - `se` (delta method on per-cell CPM; +inf where a mean is 0; NaN where not measured);
  - `de_sig` (sparse bool);
  - `log1p_bulk_50k_matched` (the sampled, thinned cells);
  - `log1p_bulk_50k_full` (all kept cells, full depth).
- **`_baseline.csv`:** per panel gene: control mean CPM, detection rate, both control profiles, `de_tested` (the gate).
- **`_control_quantiles.npz`:** for each gated gene, the sorted unique non-zero log1p(CPM) values of the sampled controls, their multiplicities and the zero count. With these, any predicted cell values can be rank-tested exactly as the scorer would, without re-reading cells. The cell generator uses them.

## Special sources
- **Marson CD4 T:** replicate-level tables only, so no Wilcoxon labels (`de_sig` empty, `n_real` NaN, `significance_unknown`).
- **Tahoe-100M:** 50 lines × drugs at 5 µM. Rows carry the drug's annotated target, but a drug is not a knockdown, so `modality = drug` gives them low weight in the model. Drugs sharing a target stay separate observations (never averaged; an agonist and an antagonist would cancel).
- **Jiang:** CRISPR knockouts under 5 cytokine conditions in 6 lines (`modality` knockout, weight ×4).

## Validation
- **Every build:** the cells are re-read and re-thinned independently, and the layers are recomputed and compared (`_build.json`).
- **`tests/test_effects_v5.py`:** checks the format.
- **Provenance:** `scripts/check_provenance.py` checks that every table was built from the QC'd cells or tables.
- **Freeze:** the 117 sorted stems are pinned by hash in `code/vcc6/CONFIG1_FROZEN.md`.

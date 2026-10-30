# Knockdown strength vs downstream effect (training data, 21 Oct 2026)

Script: `ks.py` in this folder (read-only on lib/qc_cells and results/qc). Runtime: 21 min on vcc2026.

**Outputs**
- `per_perturbation_<source>.csv`
- `bin_counts.csv`
- `per_bin_summary.csv`
- `correlations.csv`
- `scatter_<source>.png`

## Method
**Normalisation.** The same as the transformer's effect tables (`scripts/build_effects_v4.py`): per-cell CPM (counts / n_umi × 1e6), then the mean per-cell CPM per group. The delta is the LFC, log2((perturbed mean CPM + 0.5) / (control mean CPM + 0.5)).

**Cells and grouping.** Cells are after cell QC (lib/qc_cells). Every construct in `guide_qc.csv` is pooled by target gene, so every perturbation is included, also those the knockdown QC drops. `pipeline_keeps` marks the ones the transformer trains on.

**knockdown_pct** = 100 (1 − perturbed / control mean CPM) of the target gene.

**Unmeasurable** means the target is not among the source's genes, or it is detected in fewer than 10% of control cells.

**downstream_l2** is the L2 norm of the delta over the source's genes, target excluded. `downstream_l2_noise_corrected` = sqrt(max(Σ delta² − Σ SE², 0)), using the pipeline's delta-method SE.

**n_DE has two versions:**
- **n_DE_z:** the pipeline's own DE count, for every perturbation. BH 5%, two-sided normal p from delta / SE, target excluded.
- **n_DE_wilcoxon:** a Wilcoxon rank-sum on per-cell CPM against 1,000 random control cells, BH 5%, target excluded. It covers a stratified sample of at most 30 perturbations per bin with n_cells ≥ 20 (H1: all 50 eligible non-Strong plus 30 Strong).

The two counts agree well (Spearman 0.88–0.99 on the sample).

**Partial correlation:** Spearman, controlling for log n_cells (residualised ranks).

## Results
See the CSV files; the summary is in the chat report of 21 Oct 2026.

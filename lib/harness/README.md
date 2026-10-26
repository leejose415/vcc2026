# Local eval harness

Scores a submission the way cell-eval2 scores the leaderboard (vcc2026 preset, rule_version 3), on two fixed
hold-outs from our data. No model is trained here.

## Evals

| eval | line | held out | reference cells |
|---|---|---|---|
| **A** | **H1** (line 2) | the whole line: every H1 knockdown with >= 100 raw cells (**292** of 300; 8 have fewer, the generator benchmark's eligibility rule) | **raw, unfiltered** `lib/h1/H1_full_2025.h5ad` |
| **A_shared** | H1 | Eval A restricted to the H1 knockdowns that K562 also keeps | the same cells as Eval A |
| **B** | **K562** (line 1) | 100 kept K562 knockdowns, drawn from those shared with H1 (filled at random from the other kept ones if fewer than 100 are shared) | `lib/qc_cells/replogle_k562_gw_qc.h5ad`, kept cells only (`guide_qc.csv` keep → `knockdown_qc.csv` keep) |

**H1 rule:** Eval A's reference comes from the raw file. It applies no cell QC and no knockdown QC. The QC'd
`lib/qc_cells/arc_h1_qc.h5ad` and `lib/effects_v5/arc_h1.h5ad` are training data only.

**Same regime in every eval**, chosen to match what the leaderboard scores (contexts A–C: 18,400 controls, median
~20k UMI):
- **Cells per knockdown:** a seeded 400 real cells (all when fewer).
- **Controls:** a seeded 2 × 18,400 (all when fewer), split into two seeded halves:
  - `context.h5ad`: what a model and the generator condition on;
  - the **scoring controls**, appended to `real.h5ad`.
- **Depth:** every cell is binomially thinned by one factor, f = min(1, 20,000 / median control panel library), seeded per source row.
  - H1 is thinned to a ~20k median.
  - K562 is already shallower (panel library median ~7.8k) and is not thinned. Eval B therefore has less DE power than the leaderboard.
- **Genes:** reindexed to the 18,533-gene panel (`data/gene_names.csv` order). Genes a source never measured are 0 counts.
  - K562 measures only 7,681 panel genes. The other ~10.8k genes are zero in Eval B's reference and context, so they fall outside the DE gate and score as unexpressed in pds / mse.
  - The panel leaves out the ribosomal (RPL/RPS), EEF1A1 and MALAT1 genes, about 32% of K562's UMIs, so K562's panel library is smaller than its whole-matrix library.

Every seed, cell count, thinning factor and source row is in `eval<name>/split.json` and `split_rows.npz`.

## Files

```
lib/harness/
  README.md, panel.csv (panel symbol + Ensembl ID, the v5 mapping), generator_settings.json
  eval{A,A_shared,B}/
    split.json         knockdowns, cells per knockdown, seeds, thinning factor, control split, generator file hashes
    split_rows.npz     source rows of every context / scoring-control / knockdown cell
    context.h5ad       context controls (obs target = target_gene = 'non-targeting', guide_id, batch)
    real.h5ad          held-out knockdowns' real cells, then the scoring controls (obs target, target_gene)
    real_cache/ ceiling/ bundle/   cell-eval2: real side cached once; bundle = generic baseline (0) + replicate anchor (1)
results/harness/<run>/eval<name>/  agg_results.csv (raw), scaled.csv, per_target.csv, summary.json
results/harness/<run>/summary.csv  one row per eval
results/harness/trivial_summary.{csv,md}
```

The large files (`real.h5ad`, `context.h5ad`, caches) were built on **vcc2027**, under the same paths. The small
ones (`split.json`, the README, settings, scores) are copied back to vcc2026.

## Scoring

- **Scorer.** Official `cell-eval2` 0.18.0, preset `vcc2026`, DE backend pdex, called through
  `vcc3_generate.OfficialEvaluator` (`code/vcc_2026_3c/vcc3_generate.py`). The steps are:
  1. `prep-cache --side real`;
  2. `run --ceiling`;
  3. `baseline` + `prep-real-bundle --anchor-splits 5` (built once per eval);
  4. per prediction, `run -ap … -ar … --cache-real …`, then `score --real-bundle …`.
- **Reported value.** The leaderboard value is `from_replicate`: per metric, and `avg_score` = overall.
- **Engine.** pdex runs the same Wilcoxon test as gpudge, the engine on the official competition path, and pdex is what every earlier H1 result used.
- **Effect tables.** These go through the **existing generator**: `vcc3_generate`'s `calibrate`, then a common control pick, then `write_arm`, exactly as `cmd_benchmark` builds an arm.
  - Settings are fixed in `generator_settings.json`: those of the vcc_2026_5 FILE 1 leaderboard file (`results/vcc_2026_5_v5_file1_s2_floor`), with detection floor on, z = 0, alpha 1, 400 cells and seed 2026.
  - The generator's controls are the eval's context, stratified by guide. Its gate comes from the scoring controls (as `h1_controls` sets it).
  - Generator file hashes are recorded in `split.json`.

## Re-scoring a new prediction

```bash
# 400-cell count h5ad per eval: obs `target` (or `target_gene`), controls 'non-targeting' (>= 1 cell), every held-out
# knockdown present; genes by symbol (reindexed to the panel). Extra knockdowns are dropped.
~/venv/bin/python scripts/harness/score.py --run <name> --eval A A_shared B --cells predA.h5ad predA.h5ad predB.h5ad

# effect table per eval, in the generator's prediction-folder format: targets.npy, genes.npy, <array>.npy (T x G log2
# effects), optional <array>_se.npy, ontarget_lfc.npy, and nochange_targets.txt (knockdowns emitted as pure controls)
~/venv/bin/python scripts/harness/score.py --run <name> --eval A A_shared B --effects dirA dirA dirB --array <array>
```

The held-out knockdowns are `split.json["knockdowns"]`. The cells to condition on are `eval<name>/context.h5ad`.
Never train on Eval A's H1 cells, or on Eval B's knockdowns in K562.

**Memory.** Scoring one prediction on Eval A or A_shared takes about 28 GB of RAM in `cell-eval2 run`. On vcc2026
(31 GB, no swap) it is OOM-killed, so score A / A_shared on vcc2027 (32 GB swap) or add swap first. Eval B fits
anywhere. Building a bundle takes about 3.5 h for A and A_shared and 35 min for B.

## Trivial submissions

`scripts/harness/trivial.py`; numbers in `results/harness/trivial_summary.md`.
- **control_paste:** every knockdown is a seeded 400 of the context controls.
- **mean_response:** cell-eval2's own generic baseline (`cell-eval2 baseline`, seed 0). It is the bundle's 0 end, so it scores ~0 by construction.
  - The baseline writes fractional counts, which the vcc2026 preset refuses as a counts submission. The harness applies unbiased stochastic rounding (floor + Bernoulli of the fraction, seeded), so the arm scores ~0 rather than exactly 0.
- **cross_copy:** the other line's v5 effect for the same knockdown, through the generator.
  - **What is copied:** the source's official LFC on its significant genes (`de_sig`), 0 elsewhere, plus its on-target LFC as the knockdown. Knockdowns the source lacks become pure control cells.
  - **Sources:** Eval A and A_shared copy **K562**. Eval B copies **H1**, because Eval B's knockdowns are K562's own, so copying K562 would copy the answer.
  - **Coverage:** Eval A reports all held-out knockdowns, and A_shared reports the shared subset. All of Eval B's knockdowns are shared with H1.

Results (2026-10-26, scaled `avg_score`; per-metric columns in `results/harness/trivial_summary.md`):

| arm | Eval A (H1, 292) | Eval A_shared (H1, 261) | Eval B (K562, 100) |
|---|---|---|---|
| control_paste | -0.0517 | -0.0307 | 0.0063 |
| mean_response | -0.0031 | 0.0008 | 0.0071 |
| cross_copy | 0.0673 | 0.0818 | 0.2376 |

- control_paste and mean_response sit at ~0, as intended.
- cross_copy is the bar a model must beat. On Eval A it copies 261 of 292 knockdowns from K562 (median 3 significant genes each); almost all of its score is pds_cosine (0.45). Copying H1 into K562 (Eval B) is worth far more (0.24), mostly from pds_cosine (0.68) and direction fidelity (0.86).
- Each Eval A / A_shared score takes 20-45 min on CPU (pdex Wilcoxon, 8 threads).

## Build

```bash
~/venv/bin/python scripts/harness/build.py --eval B          # then A, then A_shared (needs A)
~/venv/bin/python scripts/harness/trivial.py
```

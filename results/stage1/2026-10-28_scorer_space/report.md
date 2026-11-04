# VCC 2026 — research report (all studies, oldest first)
Combined on 2026-10-30 from the six separate reports, at Jose's request: **all research is written into this one report from now on**, as new dated sections at the end.
- Each study keeps its own text. Its headings sit one level lower under a dated title, and paths that were relative to its folder are now repo-relative.
- Each study's tables, logs and outputs stay in its own folder, named under its title. The old report files are now one-line pointers here; the untouched originals are in `results/archive/2026-10-30_reports_premerge/`.
- **Tidied on 2026-11-04:** compute logistics (VM starts, swap and disk budgets, out-of-memory incidents, waiters, thread limits, file staging, runtimes) were removed so the report reads as research. The untidied version is `results/archive/2026-11-04_report_pre_tidy/report.md`.
- Within a study the order is the original one. Study 2 opens with its "current state" summary and ends with a dated log that runs from 10-26 to 10-28 and points ahead to studies 4–6.

**Contents:**
1. 2026-10-21 — Knockdown strength vs downstream effect (training data)
2. 2026-10-26 to 10-28 — Stage 1: RPE1 dev run, validation setup and the dated log
3. 2026-10-27 — Shrinkage and source-dependence study, 3 folds (archived)
4. 2026-10-27 to 10-28 — D/E/F incumbent: Model 2 (effects_v5, no floor) vs C+G on Eval A
5. 2026-10-28 — D/E/F final file, C+G (submitted 10:07 UTC)
6. 2026-10-28 to 10-30 — Scorer-space study (pseudobulk MSE, proxies, configs 1 and 1′, knockdown strength)
7. 2026-10-30 — Illustrations: profiles, changes and DE calls (H1, HCT116, context D)
8. 2026-10-30 — Step 15c: config 1′ on cells (Eval A, HCT116final)
9. 2026-10-30 — Step 17: call allocation by predictor (b), volume-matched (proxy)
10. 2026-10-30 — Step 18: D/E/F file for config 1′ + K0 (entry RTXvB7VtfATH5NQLobww, resubmitted 14:54 UTC)
11. 2026-11-01 — Step 19: studies before the next D/E/F file (hybrid calls, power, fingerprints, tier and gate rows, programs and E-reg)
12. 2026-11-02 — Step 20: lab effects, the strong-knockdown shift and program regression v2 (lab vs biology, same-lab upweighting, S and S_H1, bases, gene features, co-expression prior, recommendation)

---

# 1. 2026-10-21 — Knockdown strength vs downstream effect (training data)
*Originally `analysis/knockdown_strength/report.md` (title: "Knockdown strength vs downstream effect (training data, 21 Oct 2026)"). Its tables, logs and outputs stay in `analysis/knockdown_strength/`.*

Script: `analysis/knockdown_strength/ks.py` in this folder (read-only on lib/qc_cells and results/qc).

**Outputs**
- `analysis/knockdown_strength/per_perturbation_<source>.csv`
- `analysis/knockdown_strength/bin_counts.csv`
- `analysis/knockdown_strength/per_bin_summary.csv`
- `analysis/knockdown_strength/correlations.csv`
- `analysis/knockdown_strength/scatter_<source>.png`

### Method
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

### Results
See the CSV files; the summary is in the chat report of 21 Oct 2026.

---

# 2. 2026-10-26 to 10-28 — Stage 1: RPE1 dev run, validation setup and the dated log
*Originally `results/stage1/2026-10-27_rpe1/REPORT.md` (title: "Stage 1 — report (2026-10-27; started as the RPE1 dev run)"). Its tables, logs and outputs stay in `results/stage1/2026-10-27_rpe1/`.*

### Current state (updated 2026-10-27 14:10 UTC)

The RPE1 dev run (Steps 0–3, below) is finished. Validation has since moved to the D/E/F final target list. The dated log at the end holds every observation and decision.

#### Validation setup (since 2026-10-27 11:00)
**Effect-space screen** (`screen.py --set final`). One line is held out at a time: every table of that cell line is dropped, and all other lines vote. It is scored on the final targets the line holds.

| Line | Truth table | Final targets | Role |
|---|---|---|---|
| K562 | `replogle_k562_gw` | 267 | counts toward the keep rule |
| HCT116 | `orion_hct116` (per-cell, rebuilt 10-27) | 278 | counts |
| HEK293T | `orion_hek293t` (per-cell, rebuilt 10-27) | 271 | counts |
| H1 | `arc_h1` | 28 | reported only |

**End to end** (frozen generator → frozen harness, salt "" only; uncertainty from `boot_score.py`, a paired bootstrap over knockdowns):
- HCT116final: 278 final targets, the weak-response line. This is the one end-to-end eval in use.
- H1 Eval A (292 knockdowns, strong response): set aside on 10-27 at 14:00 (Jose), to finish HCT116 sooner.
- K562final: cancelled. RPE1dev: retired.

**Keep rule:** NMAE, MSE and sign must improve on at least 2 of K562 / HCT116 / HEK293T, with no tercile worse beyond noise. Fidelity and call count never count.

**How a prediction is built (C+G):**
1. **Sources.** All 117 tables, minus every table from the target's cell line.
2. **One vote per cell line.** A line's tables are combined by a random-effects average weighted by 1 / (SE² × multiplier). Multipliers: single-cell ×1; rest-of-line or knockout ×4; Marson, drug and Orion pseudobulk ×16. A line's variance is never below its best single row's.
3. **Pool across lines.** The line votes are averaged with 1/SE² weights (DerSimonian–Laird) to give mu.
4. **Shrink.** lfc = mu × v0 / (v0 + var_fe + τ²).
5. **Calls.** BH on |lfc| / the target's 400-cell measurement SE.
6. **The reference model** is C+G × β, with one β fitted on held-out source lines (through the origin).

Note: τ² sits at its 1e-4 floor on every screen, so step 4 barely shrinks; β ≈ 0.2 does that work instead.

#### Baselines on the new setup

| | K562 | HCT116 | HEK293T | H1 (reported) |
|---|---|---|---|---|
| NMAE — C+G / C+G × β (predict nothing = 1) | 1.121 / 1.004 | 1.151 / 1.006 | 1.338 / 1.017 | 1.406 / 1.018 |
| MSE — C+G / C+G × β (predict nothing) | 0.0755 / 0.0591 (0.0585) | 0.0895 / 0.0737 (0.0735) | 0.0714 / 0.0414 (0.0402) | 0.0417 / 0.0248 (0.0249) |
| Sign on true DE (same for both) | 0.692 | 0.848 | 0.926 | 0.559 |
| pb_nmse — C+G / C+G × β (mean response) | 1.35 / 1.02 (0.99) | 1.39 / 1.04 (0.99) | 1.86 / 1.11 (0.99) | 1.60 / 1.00 (0.92) |
| β (held-out source lines) / scale the line wants | 0.195 / 0.071 | 0.225 / 0.055 | 0.223 / 0.028 | 0.205 / 0.127 |
| r², out of sample / noise ceiling | 0.0006 / 0.136 | 0.0002 / 0.123 | 0.0002 / 0.143 | 0.0098 / 0.606 |
| End-to-end scaled (C+G) | — | 0.190 [0.144, 0.235] (MSE 0, Jaccard −1.57) | — | set aside |

Step 1 arms: E is identical to β (a flat curve), and M loses sign on every line; neither is kept. Figures are in `results/archive/2026-10-27/step1_scale/final/figures/` (archived 10-27 15:07).

#### Running now
- **HCT116final end to end:** built (14:28); C+G scored 0.190 [0.144, 0.235].
- **Nothing is running** (15:05, Jose: all study, harness and sweep jobs stopped; study outputs archived in `results/archive/2026-10-27/`).

---

## RPE1 dev run (Steps 0–3)

**Bottom line.**
- **The scoring chain works.** The effect table's own LFC and DE calls, passed through the frozen generator and the frozen harness, score a scaled overall of +2.10 to +2.24, with raw MSE −0.54. Step 1's premise (an oracle raw MSE of 0.983–1) does not hold for this harness.
- **The model was the problem.** The precision-weighted model, as written in `vcc_2026_6`, scored −0.20 end to end, below the generic baseline. Its slab prior multiplies effects by about 0.003.
- **The magnitude fix works.** Dropping the slab and keeping only the Gaussian shrinkage (G) gives **+0.122**. That is the first stage-1 model above the baseline.
- **Arms A, B and C are within noise of each other.** Tahoe covers only 1 of the 300 knockdowns.
- **Decision (Jose):** the baseline is arm C with the Gaussian shrink.

Setup:
- **Scored list:** 300 fixed RPE1 knockdowns, `results/target_lists/rpe1_dev_300.txt` (numpy `default_rng(2026)`, drawn once).
- **Held out:** the whole RPE1 line.
- **Code:** `code/vcc6/`.
- **Runner:** `vcc_2026_6.ipynb`.
- **Scorer:** cell-eval2 0.18.0, vcc2026 profile, through the harness in `scripts/harness/` (not edited).
- **Generator:** the frozen `code/cellgen/cell_generator.py`, sha256 `422f8a3d…`, unchanged.

---

### Step 0 — one copy of the model, a fixed list

**(a) The move.** The notebook's model cells went into `code/vcc6/stage1.py` verbatim, and the notebook imports it.
- **Changes:** `knockdown_seed(name, salt="")` (the default reproduces the old seeds), and `calibration_stats` written to `summary["fit"]` for the shrunk and the unshrunk predictions.
- **`calibration_stats`:** this function did not exist anywhere on vcc2026: not in the deleted `stage1_model.py`, the notebook or the bytecode. It is new. It reports Pearson r, slope of true on predicted (calibration) and slope of predicted on true (transmission), pooled over each knockdown's tested genes.
- **Self-test:** passes, including new checks for the salt, `calibration_stats` and the Step 2 arms.
- **Identity:** `predictions.npz` on 20 RPE1 knockdowns is **byte-identical** before and after the move (`results/stage1/2026-10-27_rpe1/step0/identity.json`; numpy writes fixed zip timestamps). The "after" run used the later, refactored `stage1.py` (arms, votes, `load_line`), so arm A is unchanged by everything since.

**(b) Archive.** `stage1_model.py` had already been deleted from `~/vcc`.
- My snapshot of it, and its stale `.pyc`, are in `archive/2026-10-26/` with a README.
- `stage1_rpe1.ipynb` did not exist on vcc2026.
- The only file that mentioned either was `vcc_2026_6.ipynb`, in markdown, and it has since been rewritten.

**(c) The list.** All 1,968 knockdowns in the `replogle_rpe1` table have kept cells, so the 300 were drawn from all 1,968. Per knockdown (`results/stage1/2026-10-27_rpe1/step1/knockdowns_300.csv`):
- **Cells:** median 77.5 real cells (IQR 49–117). 295 of the 300 have fewer than 400.
- **Table cells = reference cells:** the cells the table was built from are the same number as the reference cells in all 300 cases. For the 295 under 400, they are the very same cells.
- **DE calls:** the table's `de_sig` count has median 407.5 (IQR 31–1,018).
- **Terciles:** by the scorer's MSE denominator (`expr_distance_unbiased`, each knockdown's noise-corrected distance from control):

| tercile | median real cells | median table de_sig | median denominator |
|---|---|---|---|
| strong | 52.5 | 1,380 | 0.0288 |
| mid | 72.5 | 408 | 0.0069 |
| weak | 114.5 | 11.5 | 0.0007 |

**Knockdowns with fewer than 400 cells.** The harness keeps all of their cells (`scripts/harness/build.py:75`, `sample_rows`: a seeded 400 when more, all of them otherwise).
- **Controls:** shuffled, capped at 2 × 18,400, and split into a context half and a scoring half (`build.py:78`). RPE1 has 11,309 controls, so 5,654 / 5,655.
- **Thinning:** the factor is `min(1, 20000 / median control panel library)` (`build.py:95`). RPE1's median panel library is 8,482, so the factor is 1 and the reference stays at RPE1's own depth, about 8.5k UMI.

**The eval** (`results/stage1/2026-10-27_rpe1/harness/lib/evalRPE1dev`) was built with the harness's own `build_cells` and `finish`. `code/vcc6/rpe1_eval.py` only points `common.LIB` and `score.RES` at this folder.

---

### Step 1 — oracle diagnosis

| arm | scaled overall | pds | mse_norm (raw) | lfc_nmae (raw) | dir fidelity (raw) | dir reach (raw) | sig Jaccard (raw) |
|---|---|---|---|---|---|---|---|
| A, real cells (the reference's own) | 2.53 | 1.34 | 1.0 (−1.14) | 6.64 (0.000) | 1.55 (1.00) | 1.08 (1.00) | 3.54 (1.00) |
| B, zero | −0.40 | 0.00 | 0.0 (1.011) | −1.40 (1.000) | −0.39 (0.10) | −0.50 (0.07) | −0.13 (0.00) |
| C, oracle through the generator | 2.24 | 1.34 | 1.0 (−0.54) | 6.60 (0.004) | 1.55 (1.00) | 1.08 (1.00) | 1.86 (0.54) |
| D, no cells: pseudobulk of control mean × 2^lfc | — | — | clamped 0 (−0.44) | — | — | — | — |
| D2, no cells: the table's stored `log1p_bulk_50k_matched` | — | — | clamped 0 (−0.57) | — | — | — | — |

Scaled values first, raw in brackets.
- **Above 1 scaled:** scaled metrics above 1 are possible because the oracle can beat the harness's replicate anchor. Here it partly does so by construction: the reference cells are the table's own cells.
- **Arm A** is the harness check: raw squared error 0, and every per-knockdown value reproduced.

**MSE decomposition** (summed over the 300 in the metric's units, using cell-eval2's own pseudobulk and helpers in `_numerator`'s order; it matches the harness's `results.csv` to 1e-9):

| arm | raw squared error | C_pred used | C_real | numerator | denominator | ρ (budget ratio) |
|---|---|---|---|---|---|---|
| A | 0.000 | 2.309 | 2.309 | −4.617 | 4.058 | 1.0 |
| B | 6.413 | 0 | 2.309 | 4.104 | 4.058 | 0 (no spread to budget) |
| C | 0.517 | 0.413 | 2.309 | −2.205 | 4.058 | 1.0 |
| D | 0.516 | 0 | 2.309 | −1.792 | 4.058 | — |
| D2 | 0.006 | 0 | 2.309 | −2.302 | 4.058 | — |

**Per knockdown, arm C** (`results/stage1/2026-10-27_rpe1/step1/armC_per_knockdown.csv`), in scorer space over the genes RPE1 measures, own gene excluded:
- **Generated Δ against real Δ:** median Pearson r 0.953 (strong 0.974, mid 0.954, weak 0.927).
- **Slope of generated Δ on input Δ (transmission):** 0.998.
- **Slope of real Δ on input Δ (definition match):** 0.910 (strong 0.936, weak 0.883).

**Where arm C's error sits:**
- **Inside vs outside the RPE1 table panel:** all of it inside. Genes RPE1 never measured are 0 on both sides.
- **Called vs uncalled genes:** 12% / 88% of the raw squared error.
- **By control-expression decile:** almost flat (10.7% in decile 1 down to 6.3% in decile 10).
- **Shared error** (each gene's mean residual across the 300): 8.5% of the squared error, against 24% for the zero arm. The top shared genes are proliferation genes (UBE2C, TOP2A, CDK1) plus IGFBP3 and PAPSS2.
- **By tercile:** the shared share is 8.8% (strong), 9.2% (mid) and 14.1% (weak).

**Checks:**
- **Gene order:** 0 mismatches between the table's lfc, `replogle_rpe1_baseline.csv`, the scorer's `real.h5ad` var_names, the prediction, and the panel (`data/gene_names.csv`), all 18,533 long.
- **Quantile sidecar:** all 8,259 gated columns are genes the table measures, and they match the baseline's `de_tested` gate exactly.
- **Units of the table's lfc:** `log2((m1 + 1e-9) / (m0 + 1e-9))` (`scripts/build_effects_v5.py:444`). m is the **mean of per-cell CPM**, each cell normalised by its own panel library (`cpm_csr`, `build_effects_v5.py:235`).
- **Units the generator assumes:** `generate_cells` inverts exactly that, mean CPM = (m0 + 1e-9)·2^lfc − 1e-9 (`code/cellgen/cell_generator.py:168`), then counts at `library_size` (`:185`).
- **Units the scorer uses:** log1p(5e4 × **group sum / group total**), a pooled ratio (`cell_eval2/moments.py:156`, TS at `:43`).
- **Seeds:** the harness's own effect-table path never calls `code/cellgen`. It uses the older vcc3 generator with **one common control pick** for every knockdown, `seed_for(cfg.seed, 'H1', 'base_controls')` (`scripts/harness/score.py:77`). Every arm here is generated outside the harness by `code/vcc6/gen_cells.py:31`, `seed = knockdown_seed(name, salt)`, and submitted as cells. That gives 300 distinct seeds for 300 knockdowns.
- **Depth:** Step 1 generated at the notebook's `library_size=20000` against a reference at about 8.5k. Generator calls realised were a median of 871.5 against 407.5 requested (400 generated cells against a median of 77.5 real). At RPE1's own depth (Step 3) realised matched requested: 402.5 against 403.5.
- **Unresolved own genes:** 26 of the 300 target labels do not resolve to a measured gene (MRPL10, MRPL16, …), so cell-eval2 scores them without own-gene exclusion.

**Diagnosis: no link loses the signal on this harness.**
- **The generator transmits its input exactly** (slope 0.998).
- **The oracle's remaining raw error (0.517) is the definition gap, not sampling.** Mean-CPM lfc applied to the control mean (D, 0.516) against the scorer's pooled pseudobulk (D2, 0.006) is about an 85× difference. Real Δ is about 9% smaller than the lfc implies.
- **The gap doesn't move today's score.** Both versions clamp at the best end, because at about 78 cells per knockdown the reference's own noise term (C_real = 2.31) is over half the denominator (4.06).
- **The link I'd fix first** is that definition gap: predict the pooled pseudobulk (`log1p_bulk_50k`) instead of mean-CPM lfc. It matters only once a model gets close to the clamp. The 0.983–1 that motivated Step 1 matches the earlier `trivial_mean_response` run on Eval B (raw 0.9830), not an oracle on this harness.

---

### Step 2 — how sources vote (effect space, the 300)

**Structure** (`results/stage1/2026-10-27_rpe1/step2/backgrounds.csv`, `tables.csv`, `structure.json`):
- **Backgrounds:** 83. All but RPE1 vote for RPE1, and all but H1 vote for H1.
- **Tahoe:** 50 lines, 94 distinct target genes in total. It covers 1 of the 300 (GSK3B). On that knockdown, Tahoe holds a median 10.8% of the precision weight across genes, and 15.7% in RPE1's bottom expression decile.
- **`fit_prior`'s 150 sampled knockdowns:** 0 have a Tahoe row.
- **Orion:** `orion_hct116_pseudobulk` (background HCT116) and `orion_hek293t_pseudobulk` (HEK293T) hold 245 and 237 of the 300. They are ×8 in arm A (the default for an unlisted modality) and ×16 in B and C. The Orion single-cell tables hold 0 of the 300.
- **Background names aren't harmonised:** `ht-29` (Tahoe) vs `ht29` (Jiang), `hepg2/c3a` vs `hepg2`. Fixing this is in the plan.

| arm | NMAE (zero 1.0) | NMAE unshrunk | MSE (zero 0.432) | sign acc. true DE | AUROC | calls / true | call-count r | r shrunk / unshrunk | calibration slope, shrunk / unshrunk | transmission slope, shrunk | τ² | v0 | π |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A every row votes | 0.9992 | 1.044 | 0.426 | 0.714 | 0.618 | 12.5 / 407.5 | 0.42 | 0.053 / 0.100 | 1.24 / 0.51 | 0.0023 | 1e-4 | 0.51 | 0.005 |
| B one vote per background, no drugs | 0.9988 | 1.055 | 0.426 | 0.714 | 0.612 | 21 / 407.5 | 0.46 | 0.065 / 0.097 | 1.21 / 0.48 | 0.0035 | 1e-4 | 0.52 | 0.020 |
| C, B + Tahoe | 0.9988 | 1.055 | 0.426 | 0.714 | 0.612 | 21 / 407.5 | 0.46 | 0.065 / 0.097 | 1.21 / 0.48 | 0.0035 | 1e-4 | 0.52 | 0.020 |

- **On the Tahoe-touched knockdown (n = 1):** NMAE A 1.0034, B 1.0012, C 1.0032. This is not interpretable.
- **What Tahoe adds:** nothing measurable on these lists. It covers 1 of the RPE1 300, 6 of H1's 300, 3 of the Eval B 100, and 1 of the competition list.

---

### Step 3 — end to end (`run_line.py`: controls only → model → frozen generator → frozen harness)

**Controls-only profile against the sidecar** (`results/stage1/2026-10-27_rpe1/step3/profile_check.json`):
- **All 11,309 RPE1 controls:** reproduces the sidecar exactly. Max |Δ| of the control mean is 0.0, the quantiles are identical, and the 8,259 gated genes are the same.
- **The context half `run_line` is given (5,654 controls):** max |Δ| of the control mean is 34.8 CPM (38% relative, on low-count genes). The per-gene quantile max |Δ| has a median of 0.17 log1p-CPM and a maximum of 4.5, in the tails. The gate is identical.

| arm | salt | scaled overall | pds | mse_norm | lfc_nmae | dir fidelity | dir reach | sig Jaccard |
|---|---|---|---|---|---|---|---|---|
| zero | "", s1, s2 | −0.405 / −0.407 / −0.406 | 0.003 | 0 | −1.405 | −0.605 | −0.51 to −0.52 | 0.087 |
| oracle | all three | 2.097 | 1.339 | 1 | 6.488 | 1.547 | 1.080 | 1.130 |
| B | "", s1, s2 | −0.197 / −0.200 / −0.199 | 0.584 | 0 | −1.194 | −0.22 | −0.38 to −0.39 | 0.02–0.03 |
| C | "", s1, s2 | −0.196 / −0.200 / −0.199 | 0.586 | 0 | −1.194 | −0.22 | −0.38 to −0.39 | 0.02–0.03 |

- **Salt spread, the noise floor:** the range of scaled overall is 0.0035 (B, C), 0.0018 (zero) and 0.0001 (oracle). B against C (0.0005) is noise.
- **By tercile, raw** (`report.tercile_table()`), for B:
  - strong: pds 0.80, normalised MSE 0.958
  - mid: pds 0.64, normalised MSE 0.953
  - weak: pds 0.71, normalised MSE 0.736

**`vcc prep` correction.** The first Step 3 runs called `vcc prep` on the harness file and it failed every time. I had not checked its exit code.
- **The format needs two things:** no control cells in the file, and a `context` column covering A, B and C in one file.
- **Fix:** `run_line.py` now writes `prediction.h5ad` (scored by the harness, controls kept) and, with `--prep` only, `submission.h5ad` (no controls; `--context A|B|C`, or an A/B/C split of a dev line's knockdowns as a format-only check).
- **Result:** G's file passes (120,000 cells, 18,533 genes, counts preserved).
- **Opt-in:** `vcc prep` does not affect the scores and is only for files to upload, as the plan says.

---

### Magnitude test (no refit; arrays from the Step 3 arm B model run)

| variant | r | calibration slope | transmission slope | NMAE | MSE | sign acc. | AUROC | calls / 407.5 |
|---|---|---|---|---|---|---|---|---|
| S, as run | 0.065 | 1.21 | 0.0035 | 0.999 | 0.426 | 0.714 | 0.612 | 20.5 |
| U, unshrunk | 0.097 | 0.48 | 0.020 | 1.055 | 0.439 | 0.714 | 0.663 | 663 |
| **G**, Gaussian shrink only | **0.099** | 0.67 | 0.015 | 1.028 | 0.424 | 0.714 | 0.652 | **403.5** |
| H, ½ × U | 0.097 | 0.96 | 0.010 | 1.006 | 0.417 | 0.714 | 0.663 | 100.5 |

G was chosen over H by the rule fixed beforehand (higher r, then calibration slope nearer 1). End to end, one salt:

| | zero | S | **G** |
|---|---|---|---|
| scaled overall | −0.405 | −0.197 | **+0.122** |
| pds | 0.00 | 0.58 | 0.88 |
| dir fidelity | −0.60 | −0.22 | 0.46 |
| sig Jaccard | 0.09 | 0.03 | 0.21 |
| lfc_nmae | −1.41 | −1.19 | −0.58 |
| dir reach | −0.51 | −0.38 | −0.24 |
| mse_norm | 0 | 0 | 0 (raw 0.952 → 0.867) |

- **By tercile:** G helps the strong tercile most (pds 0.80 → 0.95, fidelity 0.11 → 0.49, normalised MSE 0.96 → 0.84). It makes the weak tercile's normalised MSE worse (0.74 → 1.43).
- **On transmission:** a transmission slope near 0.5 is not reachable at r ≈ 0.1. It equals r × sd(pred) / sd(true). The figure that was about 0.5 is U's calibration slope.

---

### Recommendation

**Adopt arm C with the Gaussian shrink (decided).** It is the first stage-1 model above the generic baseline, and the gain is far outside the salt noise.
- **Weakest metrics:** lfc NMAE (−0.58) and direction reach (−0.24), together with the weak tercile's MSE.
- **What would move them:** a calibrated magnitude (G is still about 1.5× too large, calibration slope 0.67), and a gate on knockdowns that do nothing in the target line.

### What I would have done differently

- **Tested the premise first.** A ten-minute oracle-only run before building Step 1 would have shown the oracle scores 2.2, and redirected Step 1 before it was built.
- **Acted on the first sign of the shrinkage problem.** Step 0's 20-knockdown run already showed π at its floor and NMAE 0.9997. The magnitude test belonged straight after Step 0, not after Step 3.
- **Checked every exit code.** I only checked `vcc prep`'s after Step 3. Its failures sat in 12 run records unnoticed.
- **Changed one thing at a time between Steps 1 and 3.** Steps 1 and 3 differ in two ways at once: library size (20k against the line's own about 8.5k) and the controls used (the sidecar's 11,309 against the context half's 5,654). So their oracle difference (2.24 against 2.10, Jaccard 0.54 against 0.34) can't be attributed to either.
- **Harmonised background names before Step 2,** since one-vote-per-background and same-line exclusion both depend on them.

---

### Study plan (agreed 2026-10-27)

**Keep rule (Jose).** A change is kept only if:
- it improves effect-space **NMAE, MSE and sign** on at least **2 of 3 lines**: RPE1 (dev 300), H1 (300, against the `arc_h1` effect table) and K562 (Eval B 100);
- and no tercile of the line it was tuned on gets worse by more than the noise.

Fidelity and call count never count towards a keep. End-to-end runs (3 salts on the RPE1 300) happen only at the milestones after 0, 1 and 2.

| # | step | cost | gates |
|---|---|---|---|
| 0 | **Lock the baseline.** C+G as a model option (`magnitude="gaussian"`, calls from it). Background merges `ht-29`=`ht29` and `hepg2/c3a`=`hepg2`, with every background string printed first; `kolf_2`, `kolf_3` and `kolf2.1j` stay separate. The three-line screen. One end-to-end run with 3 salts. | about ½ day, mostly end to end | everything |
| 1 | **Shrinkage calibration.** τ² per gene class from cross-background pairs. Earlier models got τ² ≈ 0 here, so a diagnostic first; if τ² is still about 0, document it and fall back to calibrating G's scale on held-out source backgrounds. | minutes per try after one table load | runs before 2 but does not gate it |
| 2 | **Target-expression gate on G,** then a monotone multiplier on the log expression ratio (target line against source line), fitted on source-line pairs, with wider uncertainty at low expression. | minutes, effect space | — |
| 3 | **Low-rank test,** alongside 2. C+G's error on rotating held-out source lines, inside and outside the program subspace (the `vcc3_ridge` programs and `k80` parallel analysis, unchanged), in scorer-space Δ and in LFC. Expected negative. | minutes | 5–7 |
| 4 | **Uncovered targets:** related knockdowns, with complex membership (STRING) plus response similarity in the source lines. Coverage measured on the list being scored (RPE1 300 first, then H1 300 and K562 100). Scored separately, with a CI. | about 1 h | — |
| 5 | λ / w grid: choose on one line, confirm on another | effect space | only if 3 is positive |
| 6 | Program ridge (v3) | effect space | only if 3 is positive and 6 beats 5 |
| 7 | ReLU, then transformer | large | only if 3 finds nonlinear structure and 6 wins |

**Parked:**
- **The definition gap** (mean-CPM against pooled pseudobulk): it matters only near the MSE clamp.
- **Tahoe-specific work:** no test list has enough drug targets.
- **The timed dry run:** done, in Step 3.

---

### Study plan v2 (2026-10-27, after the external research brief)

The brief changes the frame more than the method. **The final deadline is Nov 5, 23:59 UTC (9 days).** The final test lines are D, E and F, released Oct 21–22. Two submissions are allowed per UTC day, only the last counts, and scores stay hidden until Nov 19.

Two gaps on this machine:
- **The CLI is too old.** The `vcc` CLI is **0.2.2**; the final phase needs **0.3.0 or newer**.
- **The final data isn't here.** `controls-final` has not been downloaded, so `data/` holds only the validation A/B/C contexts and `pert_counts.csv`.

Without these, nothing we build can be submitted. The plan therefore gains a submission track ahead of the model steps.

#### What the brief confirms in our setup
- **Scorer:** cell-eval2 vcc2026, with the six metrics as implemented, scaled from baseline 0 to replicate 1. MSE is a ratio of sums with the noise credit capped, and pds is cosine with all 300 targets removed. Our harness runs that package (0.18.0) unchanged.
- **The MSE baseline sits at raw 0.986–0.992, and strong entries reach raw 0.70–0.85.** C+G is at raw 0.867 on RPE1, though RPE1 isn't D/E/F.
- **Tahoe is low value** (CC0, drug data; useful mostly for basal states). That matches Step 2.
- **PIE and text embeddings add little.** That matches the decision to drop PIE.
- **Shrinkage and amplitude can't be tuned on D/E/F;** they must come from cross-validation or a rule learned on source lines. That is exactly Step 1, arms E and M.

#### Rules risk in our current stack (needs Jose's decision)
1. **The frozen generator shapes significance.** For every gene in the predicted `de_sig`, it pulls the control-shaped cells toward the mean until the Wilcoxon test rejects. The brief lists "DE-significance shaping: narrow or reshape cell spread at a fixed mean so chosen genes test significant" as the **highest rules risk** of all, unanswered by staff. It is also the source of our Fidelity and Jaccard gains.
2. **Building cells from the target line's control distribution** (resampled or quantile-based NTC cells): grey, no staff ruling. Peers say "almost everybody" does it.
3. **Precision-weighted pooling of published deltas, plus shrinkage.** The 9/17 rule allows training or fine-tuning on published results, but not pasting them in or revising predictions with them. A learned model (fitted scale and weights, cross-validated) is on the safer side. A plain weighted average of public deltas is closer to "pasting". Arms E and M, and any learned weighting, move us toward "learned".
- **Suggested action:** one written question to help@virtualcellchallenge.org covering 1–3. Keep a fallback generator mode that doesn't shape significance ready. That has to be a new module with its own oracle test, because `cell_generator.py` is frozen.

#### Ordered plan

| # | Step | Why now | Cost |
|---|---|---|---|
| **S1** | **Submission readiness.** Upgrade the CLI to 0.3.0 or newer, download `controls-final`, and get the D/E/F target lists and context labels. Coverage of the D/E/F targets in our sources gives the real "uncovered targets" count. | Nothing can be submitted without it; the deadline is fixed | < 1 h |
| **S2** | **Rules question and decision** (items 1–3 above) | Decides whether the generator can be used as it is | Jose |
| **S3** | **Safe incumbent:** C+G (or Step 1's winner if ready) on D/E/F through `run_line.py --context`, then `vcc prep --prep`, then **one submission (with Jose's approval)**. Afterwards, one candidate per day at most, and only if it wins locally by more than the noise (salt range 0.0035; the brief puts replicate seed noise at 0.5–4.8% of each metric's span). | The final is blind and the last entry counts, so keep a good one in place | ~2 h compute |
| 1 | **Step 1, magnitude (running):** E (scale by the knocked-down gene's expression) and M (isotonic curve), each against C+G, under the agreed keep rule. If neither wins, add **L, a per-line amplitude predicted from control features** (the brief's "learned amplitude"). The three lines need about 0.67, 0.18 and 0.22, and a curve on the knocked-down gene's expression may not see a line-level difference. | Biggest lever on NMAE and MSE; G is 1.5–5.5× too large | ~30 min per arm |
| 2 | **Target-expression gate on the Step 1 winner** (agreed) | Strongest piece of biology | minutes |
| 3 | **Lab-aware screen.** Add a fourth held-out line from another lab (Nadig Jurkat or HepG2). When a line is held out, also drop its same-lab siblings (Replogle K562 / RPE1; the Nadig pair), at least as a sensitivity check. Weight H1 most: D/E/F are fresh Arc lines and H1 is Arc's own Flex data. | The brief: "same line from different labs less similar than different lines from one lab"; local scorers overestimate | ~15 min per screen |
| 4 | **Scorer alignment: predict on the pooled profile.** Promoted from parked. Step 1 found the oracle's whole raw error is the mean-CPM vs pooled-pseudobulk gap (0.516 vs 0.006), and the brief says to build on the scorer's pooled control profile. Investigate passing pooled-space targets through the frozen generator's mean inversion. | Raises MSE once a model gets near the clamp | half a day |
| 5 | **Low-rank test, extended with a shared response plus per-gene responsiveness** (the HKU compendium's finding; the brief notes cross-context similarity is barely penalised while cross-perturbation similarity is penalised). Expected negative on knockdown-specific structure. A per-gene responsiveness prior may still help the magnitude. | Gates the ridge, ReLU and transformer | minutes |
| 6 | Uncovered targets (now counted on D/E/F) | Pure gain where we predict 0 today | ~1 h |
| 7 | Nuisance-axis removal of the target's control variation from the predicted deltas (brief: PDS +0.05 to +0.08) | PDS is our strongest metric; optional | half a day |
| — | λ / w grid, v3 ridge, ReLU, transformer | Gated as before. With 9 days left, the transformer is out of scope. | — |

**Not doing:**
- **More DE shaping** than the frozen generator already does, pending S2.
- **Gaussian-noise exploits** (unverified claim).
- **Direct pasting of public deltas.**
- **PIE.**

**Treated as unverified** (peer claims in the brief):
- the guessed line identities (A/B/C = Jurkat, HeLa, Cal-33);
- "gains over 0.003 transfer";
- the 0.25 noise score.

#### Update 2026-10-27 08:10 UTC

**S1 is done.**
- The `vcc` CLI is now 0.3.0.
- `controls-final` is downloaded and checksum-verified, in `data_final/extracted/`: context_D/E/F.h5ad (18,400 controls each), and one list of 300 targets for all three contexts (`pert_counts.csv`; 129 shared with validation A–C). `gene_names.csv` is identical.

**Final-target coverage:**
- 299 of 300 targets are covered by at least one table. **Only ANG is uncovered.** The median is 6 tables per target.
- The sources that hold them are not the dev screen's:
  - Marson CD4 replicate tables (×16, about 275 targets each);
  - Replogle K562 (267);
  - Orion pseudobulk HCT116 / HEK293T (160 / 152);
  - Orion single-cell (about 118);
  - KOLF (50);
  - H1 (28);
  - 48 Tahoe tables cover at least one target.
  - **RPE1, HepG2 and Jurkat cover none.**
- **Implication:** the dev screen should be rebuilt around the final targets. Hold out K562 scored on its 267 final targets, and H1 on its 28. The Orion pseudobulk lines are possible held-out lines for effect space, with no DE labels. That way the sources voting in development are the ones that will vote on D/E/F.

**MSE: the generic baseline on our RPE1 eval is raw 0.664.** It's 0.986–0.992 on the validation contexts and 0.954 on H1. RPE1's knockdowns share a strong common response. C+G's raw 0.867 is worse than that baseline, so its scaled MSE clamps to 0. To beat MSE, the model must predict the shared response of the target line as well as the knockdown-specific part. The baseline gets the shared response from held-out perturbed cells, which we never see. That makes the "shared response plus per-gene responsiveness" step an MSE lever, not optional.

**Every regression now reports R² and β** (`calibration_stats`: `beta` = slope of true on predicted, `r2_fit` = r², `r2_raw` = 1 − SS_res/SS_tot of the predictions as they are). In-sample is the fitting population; out-of-sample is the held-out target line.

---

### Log of observations and decisions (UTC)

Kept in time order. "Decision" means Jose decided it; "Observation" is a measured result; "Action" is something done.

#### 2026-10-26 — source tables
- **Decision:** Jiang genes with two symbols on one Ensembl ID keep the panel's symbol and drop the alias; never sum them.
- **Decision:** Jiang's blank Ensembl IDs are filled from the panel by symbol (617 genes per file).
- **Action:** `write_replogle_qc` now maps `''` to `unresolved:<symbol>` for new contexts. QC had pooled 2–4 blank-ID targets into one knockdown in all 30 Jiang sources. QC and the 30 tables were redone, all passing.
- **Observation:** Tahoe 50/50 and Jiang 30/30 tables built and passing. `qc_effects_extra.py`: 80 tables, 0 failing.
- **Decision:** a source never votes on a target of its own cell line (`obs.background`), whatever its stem.
- **Decision:** Tahoe drugs that share a target are separate observations (×16), never averaged.
- **Decision:** small Tahoe tables are kept.
- **Action:** `stage1_model.py` moved into `vcc_2026_6.ipynb`. The source list was checked: 117 of 117 tables. Orion pseudobulk falls to the ×8 default.

#### 2026-10-27 — stage-1 RPE1 dev run (Steps 0–3)
- **Action (Step 0):** model moved to `code/vcc6/stage1.py`; predictions byte-identical to the notebook's on 20 RPE1 knockdowns. Fixed list `results/target_lists/rpe1_dev_300.txt` (seed 2026). Eval RPE1dev built with the frozen harness.
- **Observation (Step 1):** the oracle scores +2.10 to +2.24 overall, raw MSE −0.54. The generator transmits its input (slope 0.998). The oracle's raw error is the mean-CPM vs pooled-pseudobulk gap (0.516 vs 0.006). **No link loses signal; the model does.**
- **Observation (Step 2):** arms A, B and C are within noise of each other. τ² sits at its 1e-4 floor; π at or near its floor. The slab shrinks effects about 300-fold. Tahoe covers 1 of the 300.
- **Observation (Step 3):** end to end, zero −0.41, B and C −0.20, oracle +2.10. Salt noise range 0.0035.
- **Observation:** `vcc prep` rejects files that contain controls, and needs one file covering contexts A, B and C. Fixed in `run_line.py`, and made opt-in (`--prep`).
- **Observation (magnitude test, no refit):** G (Gaussian shrink only) gives r 0.099, 403.5 calls against 407.5 true, and **+0.122** end to end.
- **Decision:** the stage-1 baseline is **arm C with the Gaussian shrink (C+G)**.
- **Decision (keep rule):** a change is kept if it improves effect-space NMAE, MSE and sign on ≥ 2 of the screen's lines, and no tercile of the line it was tuned on gets worse by more than the noise. Fidelity and call count never count. End to end only at milestones, with 3 salts.
- **Decision:** background merges `ht-29` = `ht29` and `hepg2/c3a` = `hepg2` only. `kolf_2`, `kolf_3` and `kolf2.1j` stay separate.
- **Decision:** H1's effect-space truth is the `arc_h1` table. Eval A end to end stays on raw H1.
- **Decision:** uncovered targets are counted on whichever list is being scored.
- **Action:** both merges applied in `stage1._line_key`; `magnitude="gaussian"` option added; `run_line.py` defaults to arm C + gaussian; three-line screen built (`code/vcc6/screen.py`).
- **Observation (baseline screen):** G is too large on every line, with calibration slopes of 0.67 (RPE1), 0.18 (H1) and 0.22 (K562). G loses to the slab on NMAE and MSE in effect space, but wins on ranking (AUROC). End to end on RPE1 it wins by +0.32.
- **Decision:** option (a): keep the rule, do not go back to the slab, fit the scale per line from source lines (arm E, by the knocked-down gene's expression). **Arm M** (isotonic magnitude curve) runs as its own arm against C+G. If both win, stacking comes later.
- **Observation:** the research brief (shared by Jose). Final deadline **Nov 5, 23:59 UTC**, lines D/E/F, CLI ≥ 0.3.0 needed, last entry counts, scores hidden until Nov 19. Rules risk in our stack: the generator's DE-significance shaping, the use of NTC cells, and pooling published deltas.
- **Decision:** S1 (CLI and final data) yes; **S2 (rules question) and S3 (incumbent submission) skipped**.
- **Action (S1):** CLI upgraded to 0.3.0; `controls-final` in `data_final/extracted/` (D/E/F, 18,400 controls each, one list of 300 targets).
- **Observation:** final-target coverage is 299 of 300 (ANG uncovered). The covering sources are Marson CD4, Replogle K562 (267), Orion pseudobulk (160 / 152), KOLF (50), H1 (28) and Tahoe. **RPE1, HepG2 and Jurkat cover none.**
- **Observation:** our RPE1 eval's generic baseline is raw MSE **0.664** (scaled 0), and its replicate is −0.024 (scaled 1). C+G's raw 0.867 would scale to −0.30, which clamps to 0. To score above 0, raw MSE must fall below 0.664.
- **Decision (MSE):** go after MSE. Predict the target line's shared response as well as the knockdown-specific part, fix the scale (Step 1), and predict on the scorer's pooled profile. Over 200 teams beat the baseline in public.
- **Decision:** track **R² and β** on every regression: `beta` (factor to scale predictions by), `r2_fit` and `r2_raw`. Report them in-sample (the data a model or curve was fitted on, i.e. its training set) and out-of-sample (the held-out line, i.e. validation), to catch overfitting.
- **Decision:** rebuild the dev screen around the final targets. K562 is scored on its 267 final targets, H1 on its 28, and Orion pseudobulk HCT116 / HEK293T as effect-space lines. This replaces the RPE1/H1/K562-Eval-B screen after Step 1.
- **Observation:** end-to-end C+G baseline, 3 salts on the RPE1 300: **0.1216 / 0.1234 / 0.1220** (range 0.0018).
- **Observation:** Orion per-cell data covers only the 300 validation A–C targets (plus 20,000 controls). The other about 16,000 knockdowns are per-knockdown sums only, with no DE labels. For the final targets that leaves about 118 per line in cells and 160 / 152 in pseudobulk; the two sets are disjoint.
- **Decision (09:35):** re-stream Orion from Hugging Face and keep per-cell data for the final targets as well; same QC as other lines; then rebuild the screen.
- **Action (09:50):** `scripts/orion_stream.py` written. It reproduces the old A–C file exactly on batch 1: same 270 cells, labels, counts and gene panel (symbol mapping, 18,106 measured genes). A first launch nearly exhausted vcc2027's memory (whole parquet batches in memory) and was stopped. The reader now streams in record batches. Full streams (A–C ∪ final, about 471 targets, 20,000 hashed controls) started at 09:53.
- **Observation (Step 1, old three-line screen, 09:49):**
  - **E (scale by the knocked-down gene's expression):** the fitted curve is nearly flat (0.12–0.24, higher only in the top expression bin), so it is effectively a global scale fitted on source lines.
    - NMAE improves on all 3 lines (RPE1 −0.031, H1 −0.462, K562 −0.085). MSE improves on H1 (−0.028) and K562 (−0.014), RPE1 within noise. Sign unchanged.
    - Calls fall to 0 / 0 / 6, because they follow the magnitude.
    - Out-of-sample β: RPE1 5.36, H1 0.79, K562 1.09. Right for H1 and K562, five times too small for RPE1.
  - **M (isotonic curve):** NMAE and MSE improve, but sign falls on all 3 lines (−0.04 to −0.12): the curve is about −0.05 at zero and pushes small predictions negative. **M fails the keep rule.**
  - **R²:** r²_fit ≤ 0.02 in and out of sample, for every variant. The gene-level fit is small, but there is no overfitting gap.
- **Proposal (pending Jose):**
  - E's magnitude with G's call set (the generator takes the two separately).
  - Run arm L (per-line amplitude) despite E's win, because RPE1's β of 5.36 shows a line-level difference.
  - Re-decide on the rebuilt screen.
- **Observation (10:06):** the first Orion HCT116 rebuild failed at the table build (433 knockdowns with no cells). Cause: `write_replogle_qc.gene_id_lookup` reads the raw K562 and H1 files, which were archived off the VM on 10-19; it skipped them silently and left every target ID blank.
- **Action:** the lookup now falls back to the QC'd K562 and H1 files (same gene table) and refuses to run if empty. It reproduces v1's 300 of 300 target IDs and resolves all 471 new targets. HCT116 was re-run from `write_replogle_qc`; HEK293T ran with the fix.
- **Observation:** not every source is cell-based. Orion pseudobulk ×2 (per-knockdown sums) and Marson CD4 ×3 (published replicate tables) have no cells and no DE labels. Tahoe's labels are flagged `significance_unknown`.
- **Decision (Jose):** after the rebuild, delete the old Orion single-cell version. The new tables keep the same names in `lib/effects_v5` so the model reads only them.
- **Decisions (Jose, after the rebuild plan):**
  - End-to-end evals: **HCT116-final and K562-final** (two lines, two labs; HEK293T only if the two disagree).
  - **E's magnitude with G's call set:** yes.
  - **5a ablation** (C+G without the five non-cell tables): yes, and it runs **before** Step 1's arms, so E, M and L are fitted on the sources that will be kept.
  - **Keep rule on the final screen:** K562, HCT116 and HEK293T count; H1 is reported only.
  - **The reference model is C+G × one β** fitted on held-out source lines. Any curve, per-line amplitude or shared-response term is a second parameter, kept only if it beats that one number on the rebuilt screen.
- **Action and observation (10:52): Orion rebuild complete.**
  - Streamed: HCT116 2.08M cells read, HEK293T 4.53M, keeping the A–C ∪ final targets and 20,000 hashed controls.
  - Cell QC: HCT116 kept 97,249 of 97,864 cells. The knockdown QC is the existing table QC.
  - Single-cell tables: `orion_hct116` **434** knockdowns (was 274) and `orion_hek293t` **418** (was 266), both `all_pass`.
  - Pseudobulk complements: **16,058** and **15,520**. Single-cell plus pseudobulk equals the QC-kept totals exactly (16,492 and 15,938).
  - Final targets held per cell: **278 (HCT116) and 271 (HEK293T)**; 0 left in pseudobulk.
  - The pseudobulk builder's own check on the 469 knockdowns held both ways: pseudobulk vs per-cell LFC Pearson median **0.990** (p10 0.975), median |Δ| 0.016.
  - Old v1 Orion files deleted (Jose's request). The new cell files are backed up as `gs://vcc2026/orion/*_cells_targets_v2.h5ad`; the 10-08 raw A–C files remain in the bucket.
- **Observation (11:00):** the screen's final set had HCT116 and HEK293T pointed at the pseudobulk tables, which hold no final targets after the rebuild. It was fixed (single-cell tables, 278 / 271 targets) and the run restarted. The final lists are K562 267, H1 28, HCT116 278, HEK293T 271.
- **Action:** final-target end-to-end evals started with the frozen harness (`code/vcc6/final_evals.py`): HCT116final (278 knockdowns, 45,452 real cells, 9,940 + 9,940 controls), with K562final queued behind it.
- **Observation (11:12): final screen, C+G baseline.**
  - The final targets are weak knockdowns in the lines that cover them: median true DE genes are 3 (K562), 1 (HCT116), 1 (HEK293T), against 272 in H1 (28 targets).
  - C+G overshoots: β = 0.071 / 0.055 / 0.028, with H1 at 0.127. r = 0.025 / 0.014 / 0.013 (H1 0.099).
  - NMAE is 1.12 to 1.41 on every line (worse than zero), and C+G calls 110–279 genes where 1–3 are true.
  - Scorer-space pb_nmse is 1.35 / 1.39 / 1.86 (H1 1.60), against 0.99 / 0.99 / 0.99 (H1 0.92) for the mean response.
  - AUROC 0.73–0.78 (H1 0.61).
- **Observation (5a, without the five non-cell tables):** NMAE and MSE improve on all 4 lines (P = 1), but sign on true DE worsens on all 4 (−0.006 to −0.069), and so do `sign_z3` and AUROC. **Not kept by the rule.** Reading: the NMAE and MSE gain is magnitude (fewer sources), which a fitted β absorbs; the sign and AUROC loss is direction the non-cell tables carry. Proposal: keep them (pending Jose).
- **Implication:** a single β fitted on K562 and Orion will be about 0.03–0.07, close to predicting zero, while H1 (the only Arc line, nearest to D/E/F) has strong effects for the same targets. The scale for D/E/F is the key open question (arm L).
- **Check (11:58):** the final-screen C+G baseline used the rebuilt single-cell Orion tables: `orion_hct116` 434 and `orion_hek293t` 418 knockdowns, with md5 identical to vcc2026's v2 build (synced at 10:55). Non-cell tables are kept (the 5a proposal).
- **Action (12:00): single β and arms E and M on the final screen** (`step_scale.py --set final`).
  - **β (the reference):** C+G × β, where β = Σxy / Σx² over the pooled held-out source-line pairs (no intercept, the target's background excluded).
  - **E and M:** compared against β and against C+G.
  - **Calls:** every variant keeps C+G's call set (Jose: E's magnitude with G's calls), so sign and AUROC are unchanged for β and E.
  - **Scorer-space pb_nmse:** added per variant (`screen.extra_metrics(lfc_all=…)`).
  - **Output:** one json per line in `results/stage1/step1_scale/final/`.
  - **Memory:** one process holds 14.8 GB after loading the tables, so the four lines run in sequence, not in parallel.
  - **Smoke test** (H1, 3 calibration knockdowns per line): β = 0.27; H1 NMAE 1.41 → 1.03, pb_nmse 1.60 → 1.01 (mean response 0.92).
- **Observation (12:30): salts vs knockdown resampling** (RPE1dev C+G, three salts, per-knockdown `results.csv`).
  - Per-knockdown values are identical across salts for 8 of 9 metrics: per-knockdown r = 1.000 for pds, MSE, distance, nmae, Jaccard and yield.
  - Only `direction_reach` moves: salt range 0.006 raw, per-knockdown r = 0.991. That metric accounts for the whole 0.0018 scaled spread (0.1216 / 0.1234 / 0.1220).
  - The bootstrap SD over knockdowns is 3–100× the salt range: pds 0.014, yield 0.012, reach 0.019, nmae 0.005, MSE 0.0008.
  - Implication: one salt plus a paired bootstrap over knockdowns measures the noise that matters. Three salts only measure generator noise, which is negligible. Pending Jose.
- **Observation (12:30): an H1 end-to-end eval already exists.** Harness Eval A (`lib/harness/evalA`, built 10-25): H1 with the whole line held out, 292 knockdowns, 400 cells each, 18,400 controls. That is the D/E/F format exactly, and it includes all 28 H1 final targets. HCT116final and K562final are both weak-response lines (median 1–3 true DE genes); H1 is Arc's own data and strong (272). Proposal: Eval A is the primary end-to-end eval and HCT116final the weak-regime check; drop the queued K562final. Pending Jose.
- **Decisions (Jose, 12:40):**
  - **End-to-end pair:** H1 (harness Eval A, strong response) and HCT116final (weak response), always reported side by side.
  - **No salt:** end-to-end runs use salt "" only. Uncertainty comes from a bootstrap over knockdowns.
  - **Step 1** goes ahead on the rebuilt final screen: single β, then E and M on top of C+G.
- **Action (12:45): `code/vcc6/boot_score.py`.** Resamples cell-eval2's per-knockdown `results.csv` (2,000 resamples, seed 2026) and rebuilds the aggregate as cell-eval2 does (means; the MSE ratio of sums). It then scales with cell-eval2's own `score_metrics` against the eval's bundle, with the anchors fixed. Runs on the same eval are compared paired.
  - The full-set rebuild reproduces the harness's `summary.json` exactly (checked on every run).
  - On RPE1 C+G: overall 0.1216, 95% CI [0.088, 0.157]. Paired salt s1 − none = +0.0018 [+0.0005, +0.0031]. So the generator's draw is a small but consistent shift, about a tenth of the noise from which knockdowns are scored.
- **Action (12:45): C+G end-to-end baseline queued** (`results/stage1/2026-10-27_rpe1/logs/e2e_final/run_CG.sh`), behind the HCT116final build:
  - H1 on Eval A (`results/target_lists/evalA_h1_292.txt`, the H1 tables held out);
  - then HCT116 on HCT116final;
  - salt "" only; outputs in `results/stage1/e2e_final/CG/`.
- **Observation (12:22): K562 final screen, β/E/M.** β = 0.195 (calibration on the 7 non-K562 lines).
  - **β vs C+G:** NMAE 1.121 → 1.004 (−0.117, CI [−0.126, −0.105]); MSE −0.016; better in every tercile. Sign and AUROC are unchanged (same calls). Scorer-space pb_nmse 1.346 → 1.024, but the mean response is 0.988.
  - **E = β:** the E curve is flat (0.18–0.25), and its difference from β is 0.000 on every metric.
  - **M vs β:** NMAE +0.028 worse, MSE equal, sign −0.149 (CI [−0.207, −0.085]). M fails.
  - **Out of sample:** K562's own slope is still 0.36 after β, so K562 wants about a third of the β fitted on the other lines.
- **Observation (12:38): HCT116 final screen, β/E/M.** β = 0.225 (calibration on the 7 non-HCT116 lines). Same pattern as K562.
  - **β vs C+G:** NMAE 1.151 → 1.006 (−0.145, CI [−0.154, −0.136]); MSE 0.0895 → 0.0737, against 0.0735 for predicting zero; better in every tercile.
  - **Scorer space:** pb_nmse 1.392 → 1.037, against 0.989 for the mean response.
  - **E = β:** the E curve is flat, 0.22–0.32.
  - **M vs β:** NMAE +0.029, sign −0.098 (CI [−0.157, −0.026]). M fails.
  - **Out of sample:** the slope after β is 0.25, so HCT116 is still overshot about 4×.
- **Observation (12:54): HEK293T final screen, β/E/M.** β = 0.223. Same pattern again.
  - **β vs C+G:** NMAE 1.338 → 1.017 (−0.321); MSE 0.0714 → 0.0414, against 0.0402 for predicting zero; better in every tercile.
  - **Scorer space:** pb_nmse 1.857 → 1.109, against 0.986 for the mean response.
  - **E = β** (flat curve, 0.21–0.30).
  - **M vs β:** NMAE +0.045, sign −0.126. M fails.
  - **Out of sample:** the slope after β is 0.12, so HEK293T is still overshot about 8×.
- **Status:** the HCT116final build is still in `prep-real-bundle` (1 h 19 min, at the bundle's DE passes). It is slower than estimated.
- **Observation (13:02): H1 final screen, β/E/M** (28 targets, reported only). β = 0.205 (calibration on the 7 non-H1 lines).
  - **β vs C+G:** NMAE 1.406 → 1.018 (−0.389); MSE 0.0417 → 0.0248, against 0.0249 for predicting zero.
  - **Scorer space:** pb_nmse 1.598 → 1.001, against 0.922 for the mean response. In the strong tercile it is 0.967.
  - **E ≈ β:** NMAE −0.001, CI [−0.004, +0.006].
  - **M vs β:** NMAE +0.132, MSE +0.003. M fails.
  - **Out of sample:** the slope after β is 0.62, so H1 is overshot about 1.6×. Even the strong line wants less than β.
- **Step 1 result on the final screen (all four lines):**
  - **β is stable:** 0.195 / 0.225 / 0.223 / 0.205, whichever line is held out.
  - **β beats C+G** on NMAE and MSE on all 4 lines and in every tercile. Sign, z3 and AUROC are unchanged by construction: same calls, positive scale.
  - **E is identical to β** (its curve is flat everywhere). **M worsens sign** on all lines. **Neither E nor M is kept.**
  - **C+G × β stays the reference.** By the letter of the keep rule, β cannot "improve sign", since a positive scale cannot change it. It is the reference model by Jose's decision, not a kept change.
- **Observation: every target line wants a smaller scale than β.** The optimal scale of C+G on the target is K562 0.071, HCT116 0.055, HEK293T 0.028, H1 0.127, against β ≈ 0.2. β is fitted on 200 random knockdowns per calibration line, but the final targets are weaker and less predictable (r = 0.01–0.10) than random knockdowns.
  - Proposal: β_final, one number fitted on the held-out source lines restricted to the final-target knockdowns. It uses source data only, no D/E/F truth.
- **Observation: no variant beats the mean response in scorer space** (pb_nmse 1.00–1.11 against 0.92–0.99). The scorer's baseline holds the line's shared response, which our per-knockdown deltas lack. This is the case for the shared + specific step of the MSE block.
- **Observation (13:30): R² on the final screen.** Figures in `results/stage1/step1_scale/final/figures/` (`plot_step_scale.py`).
  - **Gene-level r², out of sample** (pooled over the final targets; the same for C+G, β and E): K562 0.0006, HCT116 0.0002, HEK293T 0.0002, H1 0.0098. M: 0.0009 / 0.0002 / 0.0003 / 0.012.
  - **In sample** (held-out source lines): 0.0005–0.0007 for every variant. There is no gap between in and out of sample, so no overfitting; but there is almost no fit to overfit.
  - **Noise ceiling** (1 − mean se² / var(true lfc), same genes; `r2_ceiling.json`): K562 0.136, HCT116 0.123, HEK293T 0.143, H1 0.606. We reach 0.1–1.6% of the r² the truth's noise allows. H1 does best (1.6%).
  - **R² with the scale counted** (1 − MSE / var): C+G is −0.11 / −0.06 / −0.20 / −0.45; C+G × β is −0.005 / −0.005 / −0.011 / +0.006. Shrinking moves us onto "predict nothing" and no further.
- **Decision (Jose, 14:00):** set aside the H1 end-to-end run and finish HCT116final as soon as possible. Add the "Current state" section at the top of this report. Start the shrinkage and source-dependence study (3 folds).
- **Action (14:05):** the C+G end-to-end queue was stopped (its H1 step never ran).
- **Observation (14:10): study Step 0 inventory, first findings.**
  - There is no control-cell similarity measure in the codebase.
  - The tables carry no dataset or study ID beyond their names. By name prefix, only feng (21 tables), jiang (30) and tahoe (50) have 3 or more lines; everything else falls into "other" under the study's rule.
  - Marson (3 tables) has no cells: no median cell count and no per-cell control profile.
- **Decisions (Jose, 14:15), shrinkage study:**
  - **Source unit:** the 117 tables, excluded by cell line. Several drugs on one target within a Tahoe table are first combined into one row with `combine_within`.
  - **Knockdowns:** the final-target lists everywhere, for the held-out line F and for the training targets U.
  - **Harness metrics:** only on existing evals, HCT116final and H1 Eval A. The K562 fold has none, and HEK293T is effect space only.
  - **Similarity:** cosine of each table's log1p mean control profile (`log1p_bulk_50k_full`) over the 2,000 genes most variable across lines, among genes every line measures.
  - **Marson:** its log cell count is imputed with the median across all sources, so it gets a t̂ and stays in Stage A. This is flagged in the study report.
- **Action (14:26): shrinkage study build started** (`code/vcc6/shrink_study.py build`). It does one load of the 117 tables and writes:
  - the inventory, the SE distribution and the similarity matrix;
  - per-target arrays for K562, HCT116, HEK293T and H1;
  - 12 baseline (C+G) runs, one per (target, removed lines) pair the folds need.
  - `analyze` (`shrink_analyze.py`) runs Steps 0.4–5 from that cache. The study report is `results/stage1/2026-10-27_shrinkage_study/report.md`.
- **Decisions (Jose, 14:35), shrinkage study, before analyze:**
  - **ρ_w / ρ_b and Stage A's R_k:** the dataset is the study (table-name prefix), whatever its number of lines. Different studies are always between-dataset. The "≥3 lines, else other" pooling applies only to the 3b fixed effects.
  - **Missing t̂:** a table whose 3b dataset has no rows in the fold's fit gets t̂ at the reference level (similarity and cell count only). It is flagged and its share of (k, g) reported; the run does not stop.
  - **Tahoe:** report what it contains (done in the study report).
  - **Harness approved:** baseline and Stage C end to end on HCT116final and H1 Eval A after analyze, each with the study's fold pool. That is 4 runs.
- **Observation (14:40): Tahoe in the study.**
  - Tahoe holds single-target drugs at their top dose against DMSO, mapped to each drug's annotated target gene: 141 drugs on 94 targets per table, with 2–6 drugs on 29 of those targets.
  - Only one Tahoe target is on any final list: MAPK7 (HCT116 and HEK293T; 48 of 50 tables hold it).
  - So Tahoe cannot enter 3b (≥10 shared knockdowns), and its t̂ will be at the reference level.
- **Observation (14:31):** study build at the per-target arrays. Similarity uses 5,961 genes that every table measures, top 2,000 by variance; off-diagonal cosine ranges 0.47–1.0 (median 0.91). The K562 array is 267 knockdowns × 18,533 genes with 1,485 source rows. Marson's imputed log cell count is 5.99 (≈400 cells).
- **Observation (14:45): Tahoe and the study's similarity** (Jose's question; `code/vcc6/similarity_check.py`).
  - 50 of the 117 tables are Tahoe.
  - Recomputing the 2,000-gene HVG set without Tahoe: both sets are drawn from the same 5,961 genes every table measures, and they overlap in 1,744 genes (87%).
  - Each fold target's similarity ranking of its 116 sources is almost unchanged: Spearman 0.997 (K562), 0.998 (HCT116), 0.997 (HEK293T), 0.994 (H1); 0.994–0.999 on non-Tahoe sources only.
  - Tahoe does not drive the similarity.
- **Observation (14:48): HCT116final end to end, C+G (screen pool: HCT116 held out; salt ""). Scaled overall 0.190.**
  - **pds:** 0.655 (raw 0.628).
  - **direction yield:** 0.915 (raw 0.540).
  - **direction reach:** 0.741 (raw 0.667).
  - **lfc_nmae:** 0.401 (raw 0.747).
  - **MSE:** 0, clamped (raw 2.90, far worse than baseline).
  - **sig Jaccard:** −1.57 (raw 0.015). C+G requests a median of 132 calls, where the truth has about 1.
  - Reading: on this weak line the direction metrics carry the score. MSE (oversized effects) and Jaccard (far too many calls) are the two losses, which is exactly what β and a smaller call set address. Bootstrap CI pending (`boot_score.py`).
- **Observation (14:55): HCT116 C+G end to end**, bootstrap over knockdowns: 0.190, 95% CI [0.144, 0.235].
- **Observation (14:51): shrinkage study, analyze done** (full write-up in `results/stage1/2026-10-27_shrinkage_study/report.md`).
  - **Headline:** on the final targets, no predictor beats predicting zero on LFC MSE by more than 0.7%. Stage C is lowest on every line, but only matches zero: −0.08% K562, −0.08% HCT116, −0.10% HEK293T, −0.72% H1.
  - **The gain is shrinkage:** B_k ≈ 0.3 × c ≈ 0.14–0.17.
  - **Stage A's similarity-modelled weights are worse than the baseline's** on all 4 lines.
  - **The sources are nearly independent:** ρ_b ≈ 0.001–0.006, ρ_w ≈ 0.04; Marson's 3 conditions are the only correlated group (0.30).
  - **3b is unreliable:** 11–22 rows from about 8 tables, the similarity coefficient flips sign across folds, and the cell-count effect is carried by Marson (imputed cells, t floored because its SE² exceeds its error).
  - **The calibration curve is linear;** slopes are 0.02–0.11.
  - **The final config's c (0.256) is not comparable to the fold c (0.14–0.17):** in Step 5, HCT116 and HEK293T are sources for each other.
- **Observation (14:52): study end-to-end inputs.**
  - The study's baseline (HCT116 and HEK293T both held out) requests a median of 464 calls on HCT116, against 132 for the screen's C+G (HCT116 held out only); H1 requests 488.
  - Stage C's median |lfc| is 0.003, against 0.04–0.06 for the baseline.
  - Both keep the baseline's calls, so the generator will shape 464–488 DE genes per knockdown at near-zero means.
- **Observation (15:00): truth on each line's final list** (the tables' `de_sig` uses the harness's rule: cell-eval2 vcc2026 Wilcoxon, BH p_adj < 0.05, own gene excluded, 400 cells at 20k depth; RMS over the tested genes).

  | Line | True DE per knockdown, mean / median (p90) | RMS true LFC per knockdown, mean / median |
  |---|---|---|
  | K562 (267) | 13.1 / 3 (36) | 0.299 / 0.239 |
  | HCT116 (278) | 10.7 / 1 (13) | 0.436 / 0.270 |
  | HEK293T (271) | 25.7 / 1 (16) | 0.284 / 0.200 |
  | H1 (28) | 849 / 272.5 (2,554) | 0.175 / 0.151 |

  - The harness's own calls on the real eval cells agree: HCT116final 10.6 / 1; Eval A on the 28 final targets 851 / 278, and on all 292 knockdowns 1,048 / 218.
  - **Eval A does not use HCT116final's list.** It is the whole H1 line (292 knockdowns, 28 of them D/E/F targets); the two lists share 27 knockdowns.
- **Decisions (Jose, 15:00):**
  - **Deciding metric:** the end-to-end harness's scaled overall decides which stage is recommended; LFC MSE is reported only.
  - **Study metrics:** add the screen's effect-space metrics per stage (NMAE, r, AUROC, sign, slope through origin, calls vs true calls by tercile), using the screen's code.
  - **B collapse check:** report B_k's distribution and the share below 0.1.
  - **Report:** add a section comparing effect-space and harness rankings.
  - **Another harness eval:** to be built after the truth statistics above have been reported.
  - **Magnitude sweep (start now):** C+G × m, m ∈ {1, 0.5, 0.25, 0.14, 0.07, 0.03}, end to end on HCT116final and H1 Eval A, screen setup.
- **Action (15:00): magnitude sweep started** (`code/vcc6/magnitude_sweep.py`; `results/stage1/2026-10-27_rpe1/logs/sweep_*.sh`).
  - **Estimate:** 12 runs (11 new; HCT116 m = 1 is the 14:48 run), about 4.5–5 h wall, ending around 20:00 UTC.
  - **Method:** calls are recomputed on m × lfc by the model's own rule. At m = 1 the recomputation reproduces the model's calls exactly. HCT116 mean calls per knockdown are 196 / 15 / 1.3 / 0.2 / 0 / 0 for m = 1 / 0.5 / 0.25 / 0.14 / 0.07 / 0.03.
- **Action (15:02):** study analyze re-run with the screen metrics and the B collapse check added.
- **Decision (Jose, 15:05): stop everything and archive the study (move, never delete).**
    - `code/vcc6/shrink_study.py` (12 KB) and `shrink_analyze.py` (40 KB).
  - **Left in place:** the frozen generator (sha256 422f8a3d…, both VMs), `scripts/harness/`, the 117 effect tables, `stage1.py` (C+G).
  - **Also left in place (not on the list):**
    - `study_e2e.py` and `similarity_check.py` (both import the archived `shrink_study.py`, so they no longer run);
    - `magnitude_sweep.py` and `results/stage1/magnitude_sweep/` (HCT116 inputs only, nothing generated);
    - `step_scale.py` and `plot_step_scale.py`;
    - the step-scale logs.
  - **Kept as measurements of the data, not study conclusions** (in this log):
    - true DE per knockdown, HCT116 median 1 and H1 median 272;
    - HCT116 C+G end to end, 0.190 [0.144, 0.235];
    - per-source r with the target, 0.001–0.07.
- **Decision (Jose, 2026-10-28 ~05:25 UTC):** **C+G is the D/E/F incumbent.** On Eval A the Model 2 port (v3's pooling on effects_v5, C+G's calls) scored 0.184 and C+G 0.219; paired difference +0.034 [+0.020, +0.048] for C+G. Run C+G on D/E/F and **upload it** (explicit OK). Write-up: `results/stage1/2026-10-27_incumbent/report.md`; run: `results/stage1/2026-10-28_DEF/`.
- **Submitted (2026-10-28 10:07 UTC): the D/E/F final entry `JFCghVNyHwfCOWhnHG6M` (LHY), C+G.**
  - `vcc prep` exit 0: 360,000 cells × 18,533 genes, D / E / F 120,000 each, targets verified.
- **Scorer-space study (2026-10-28), consolidated** in `results/stage1/2026-10-28_scorer_space/report.md`.
  - The metric is reproduced exactly. The noise credit is not the gap; the plug-in distance is. A noise-free proxy matches end-to-end MSE to 0.0006.
  - **Recommendation pending Jose:** C+G LFC × 0.5 (m* = 0.5 on every fold). Eval A MSE 0 → +0.12, PDS −0.03; expected overall about +0.015 on Eval A.

---

# 3. 2026-10-27 — Shrinkage and source-dependence study, 3 folds (archived)
*Originally `results/archive/2026-10-27/2026-10-27_shrinkage_study/report.md` (title: "Shrinkage and source-dependence study (3 folds) — report"). Its tables, logs and outputs stay in `results/archive/2026-10-27/2026-10-27_shrinkage_study/`.*

Status (2026-10-27 14:25 UTC): Step 0 inventory done; build running (`results/archive/2026-10-27/code/vcc6/shrink_study.py build`), then
`analyze`. Results sections below are filled in as the steps finish.

Code: `results/archive/2026-10-27/code/vcc6/shrink_study.py` (build: one table load → cache) and `results/archive/2026-10-27/code/vcc6/shrink_analyze.py` (Steps 0.4–5).
Outputs: this folder: `results/archive/2026-10-27/2026-10-27_shrinkage_study/tables/*.parquet`, `results/archive/2026-10-27/2026-10-27_shrinkage_study/figures/*.png`, `results/archive/2026-10-27/2026-10-27_shrinkage_study/results.json`, `results/archive/2026-10-27/2026-10-27_shrinkage_study/config_final.json`, `results/archive/2026-10-27/2026-10-27_shrinkage_study/stageC_<line>.npz`;
`results/archive/2026-10-27/2026-10-27_shrinkage_study/cache/` (per-target arrays and baseline predictions; copied to `gs://vcc2026/shrinkage_study/`).

### Decisions taken before the build (Jose, 2026-10-27 14:15)

| Question | Decision |
|---|---|
| What is a source *s* | The 117 tables. Leakage exclusion is by cell line (obs.background, with the two agreed aliases). Several Tahoe drugs on one target in one table are combined into one row first (`stage1.combine_within`, multiplier 1). |
| Knockdowns | The final-target lists everywhere: K562 267, HCT116 278, HEK293T 271, H1 28. For H1 also Eval A's 292 knockdowns, used only for its harness run. |
| Harness metrics | Existing evals only: HCT116final (HCT116/293T fold) and H1 Eval A (H1 fold). The K562 fold has no harness metrics; HEK293T is effect space only. |
| Similarity | The codebase has none, so the spec's default with one change of profile: cosine between each table's log1p mean control profile (`log1p_bulk_50k_full`, the log of mean CPM at 50k depth; mean-of-log per cell does not exist for every table), over the 2,000 genes most variable across the 117 tables among genes every table measures. |
| Marson (no cells) | Its log median cell count is **imputed with the median over all sources**, so it gets a t̂ and stays in Stage A. Flagged here and in the config. |

Choices fixed by the codebase rather than the spec (stated so they can be challenged):
- A source value counts where its LFC and SE are finite and |LFC| ≤ 6. That is the epsilon-floor artifact rule every stage-1 pool already applies.
- The target's (k, g) set is its `de_tested` genes (5-CPM gate) with a finite truth, own gene excluded: the screen's scoring set. Truth is the target's own effect table.
- "Current precision-weighted transfer model" = stage-1 arm C with the Gaussian shrink (C+G), no β.
- Harness runs (baseline and Stage C) keep the baseline's DE call set. The spec does not say how Stage C calls; magnitude arms have kept the baseline calls since 10-27.
- Similarity of a source to a target = similarity between the source table and the target's own table (its controls only).
- Dataset ID: the tables carry none beyond their names, so the name's study prefix is used.

### Changes before analyze (Jose, 2026-10-27 14:35)

| Item | Change |
|---|---|
| Dataset for ρ_w / ρ_b (and R_k in Stage A) | The **study** (table-name prefix), whatever its number of lines. Pairs from different studies are always between-dataset. The "≥ 3 lines, else other" pooling is kept **only** for the 3b fixed effects. |
| Missing t̂ | A source table in S_k whose 3b dataset has no rows in that fold's fit gets t̂ from the 3b model **at the reference level** (similarity and cell count only). It is flagged, and the share of predicted (k, g) whose weights involve such tables is reported. The run is not stopped. |
| Tahoe | Reported below. |
| Harness | **Approved:** baseline and Stage C end to end on HCT116final and H1 Eval A after analyze. That is 4 runs, each with the study's fold-specific source pool (HCT116: HCT116 and HEK293T both held out; H1: H1 held out). |

**Study mapping (table-name prefix → tables; cell lines; 3b label):**

| Study | Tables | Lines | 3b label |
|---|---|---|---|
| tahoe | 50 | 50 | tahoe |
| jiang | 30 | 6 | jiang |
| feng | 21 | 20 | feng |
| adamson | 3 | 1 (K562) | other |
| marson | 3 | 1 (CD4 T) | other |
| orion | 4 | 2 (HCT116, HEK293T) | other |
| replogle | 2 | 2 (K562, RPE1) | other |
| nadig | 2 | 2 (HepG2, Jurkat) | other |
| arc | 1 | 1 (H1) | other |
| kolf | 1 | 1 (KOLF2.1J) | other |

**Tahoe:**
- **Perturbation type:** small-molecule drugs (Tahoe-100M), each at its highest dose, against a plate-matched DMSO vehicle control. Each drug is single-target and is mapped to its annotated target gene; that gene is the "knockdown" ID. The tables carry the note "drug effects are not on-target verified".
- **Size:** 141 drug conditions on 94 target genes per table (two small tables have 18 and 16 drugs).
- **Several drugs per target:** 29 targets have 2–6 drugs. In this study each table's drugs on one target are combined into one row (`combine_within`).
- **Depth:** about 1.7k UMI per cell, 10× below the scorer, so DE labels are flagged `significance_unknown`.
- **Overlap with the final lists:** the 94 Tahoe targets overlap the final lists in **one** gene, MAPK7, held by 48 of the 50 tables. That gives K562 0, HCT116 1, HEK293T 1, H1 0.
  - So no Tahoe table can meet Step 3's ≥ 10 shared knockdowns with any training target. Tahoe gets no 3b rows; its t̂ will be at the reference level.
  - Tahoe enters only MAPK7's prediction in HCT116 and HEK293T.
- **Tahoe's share of the tables:** 50 of the 117 tables (43%) are Tahoe, so the tables are dominated by Tahoe in number but not in coverage of the final targets.
- **Does Tahoe drive the similarity?** No (`code/vcc6/similarity_check.py`; `results/archive/2026-10-27/2026-10-27_shrinkage_study/tables/similarity_tahoe_check.{parquet,json}`). The HVG set was recomputed over the 67 non-Tahoe tables only:
  - **Genes every table measures:** 5,961 in both cases. The non-Tahoe tables set that limit.
  - **The two 2,000-gene sets overlap in 1,744 genes (87%).**
  - **Spearman correlation of each fold target's similarity to its 116 source tables, all-tables genes vs non-Tahoe genes:**

    | Target | All sources | Non-Tahoe sources only |
    |---|---|---|
    | K562 | 0.997 | 0.994 |
    | HCT116 | 0.998 | 0.999 |
    | HEK293T | 0.997 | 0.999 |
    | H1 | 0.994 | 0.998 |

  So the ranking of sources by similarity is essentially unchanged; Tahoe tables shift about 13% of the HVGs but not the ordering.

### Step 0 — inventory

#### 0.1 Where things are
| Item | Location |
|---|---|
| Transfer model | `code/vcc6/stage1.py` (`predict_line`, arm C, magnitude `gaussian`) |
| Control-cell similarity | **None in the codebase** → the spec's default (see decisions) |
| Effect tables | `lib/effects_v5/` and `lib/effects_v5_extra/`: 117 tables (`stage1.discover_tables`), each with LFC and SE layers, a `_baseline.csv` (control profile, `de_tested`) and a control-quantile sidecar |
| Fold definitions | `code/vcc6/screen.py` `FINAL` (line → truth table) and `results/target_lists/final_<line>.txt` (D/E/F targets ∩ table) |
| Eval harness | `scripts/harness/` (frozen), driven by `code/vcc6/run_line.py score`; evals `lib/harness/evalA` (H1) and `results/stage1/final_evals/harness/lib/evalHCT116final` |

#### How the baseline sets its weights (C+G)
1. **Sources:** every table whose cell line is not the target's.
2. **One vote per cell line:** a line's tables (and a Tahoe table's several drugs on one target) are combined by a DerSimonian–Laird random-effects mean, with weights 1 / (SE² × m).
   - m = 1 for single-cell with non-targeting controls; 4 for rest-of-line reference or CRISPR knockout; 16 for Marson replicate tables, Tahoe drugs and Orion pseudobulk.
   - Non-single-cell SE is floored at 0.25.
   - A line's variance is floored at its best single row's, so repeat screens add no precision.
3. **Pool across lines:** the votes are averaged per gene with weights 1/SE² (DL random effects). That gives mu, var_fe = 1/Σw, and τ² (gene-specific DL, floored at a global τ² fitted on knockdowns seen in at least 2 lines).
4. **Gaussian shrink:** lfc = mu × v0 / (v0 + var_fe + τ²). v0 is the mean excess variance of effects with |z| > 2 (about 0.5–0.76).
   - The global τ² has sat at its 1e-4 floor on every screen.
   - So the shrink factor is about 0.98 for a single-cell vote with SE 0.1, and about 0.33 only for the ×16 sources.
5. **Weights are never correlation-aware:** two tables of one line count once (step 2), but two lines from one study count as independent. Nothing uses control similarity.

### Summary

**On the final targets, none of the four predictors (baseline, A, B, C) beats predicting zero on LFC MSE by more than 0.7%.**
- Stage C has the lowest MSE on every held-out line, but it reaches "predict nothing" and no further: −0.08% (K562), −0.08% (HCT116), −0.10% (HEK293T), −0.72% (H1).
- The baseline is 6–49% worse than zero. Almost all of that is removed by **shrinkage** (Stage B ≈ ×0.3, then c ≈ ×0.14–0.17), not by better weights.
- **Stage A's correlation-aware, similarity-modelled weights make MSE slightly worse** than the baseline's on all four lines.
- **The reason is that the sources barely agree with each other or with the target on these knockdowns:** source–source error correlation net of the target is ρ_b ≈ 0.001–0.006 and ρ_w ≈ 0.04. Per-source r with the target is 0.001–0.07.
- **The transfer-variance model (3b) is fitted on 11–22 source–target pairs drawn from about 8 tables per fold, and its similarity coefficient changes sign across folds.** The final config's 3b is not reliable for unseen lines.
- **Harness metrics** (baseline vs C on HCT116final and H1 Eval A) are pending.

### Direct answers

| Question | Answer |
|---|---|
| **Is the pooled slope well below 1?** | Yes, far below. y on p through the origin: K562 0.053, HCT116 0.034, HEK293T 0.018, H1 0.113. Huber is lower still (0.036 / 0.022 / 0.016 / 0.063). OLS intercepts are −0.022 / −0.034 / −0.023 / −0.004. |
| **Does the slope rise with source count?** | In K562 yes, monotonically: 0.017 (1 source) → 0.029 → 0.030 → 0.066 (5–9) → 0.073 (10–29, 13 knockdowns). In H1 broadly (0.05 → 0.09 → 0.03 → 0.15). In HCT116 and HEK293T no clear trend (0.01–0.05). Their 10+ buckets come from only 14 knockdowns and are ≈ 0. |
| **Linear or S-shaped?** | **Linear on all four lines.** The through-origin slopes below and above the 90th percentile of \|p\| have overlapping knockdown-bootstrap CIs: K562 0.049 vs 0.055, HCT116 0.038 vs 0.031, HEK293T 0.025 vs 0.014, H1 0.109 vs 0.115. |
| **How dependent are the sources?** | By 3c, barely: ρ_w = 0.050 / 0.036 / 0.040 and ρ_b = 0.006 / 0.001 / 0.004 (folds K562 / HCT116-293T / H1), so n_eff ≈ S (S = 30 → n_eff 25–29). Marson's three conditions (ρ ≈ 0.30) are the only correlated group. The 3d plateau implies ρ = 0.65 / 0.87 / 0.35 / 0.23 (K562 / HCT116 / HEK293T / H1), but on only 13 / 14 / 14 / 1 knockdowns, and it measures something else (see "Undermines the design"). |
| **Does transfer variance track similarity, cell count or dataset?** | Not reliably. 3b's R² is 0.72 / 0.79 / 0.99, but on 19 / 11 / 22 rows from about 8 tables, with 6–9 rows floored. The similarity coefficient is −36.7 / −34.3 / +13.1 (sign flip). Cell count (−4.4 / −4.0 / −3.7) and feng (−4.7 / −2.6 / −4.3) are stable in sign, but the cell-count effect is carried by Marson: imputed ~400 cells, and t floored because its SE² exceeds its whole error. Raw t against cell count shows no trend apart from Marson. |
| **How variable is β_k, and does agreement predict it?** | Very variable around a small mean. β_k p5–p95: K562 −0.04…0.23 (median 0.035), HCT116 −0.03…0.21 (0.022), HEK293T −0.03…0.12 (0.016), H1 0.006…0.40 (0.047). Agreement predicts it weakly: Spearman β_k vs τ̂²_k = 0.27 / 0.28 / 0.25 / 0.04, and B_k vs β_k = 0.31 / 0.24 / 0.21 / 0.03. β_k falls slightly with strength (−0.19 / −0.16 / −0.13 / −0.06). |

### Step 0 — results

#### 0.4 Coverage (sources = every table except the fold's lines)

| Fold line | Knockdowns listed | With truth and ≥ 1 source | (k, g) with 0 sources | 1 | 2 | 3–4 | 5–9 | 10–29 | 30+ |
|---|---|---|---|---|---|---|---|---|---|
| K562 | 267 | 267 | 24,659 | 26,070 | 70,605 | 319,415 | 1,514,071 | 95,778 | 0 |
| HCT116 | 278 | 276 | 55,482 | 91,790 | 70,056 | 1,906,437 | 618,480 | 116,286 | 10,161 |
| HEK293T | 271 | 269 | 62,107 | 85,565 | 80,897 | 1,992,348 | 645,448 | 121,527 | 10,733 |
| H1 | 28 | 28 | 1,100 | 2,970 | 5,824 | 51,059 | 231,006 | 9,573 | 0 |

Most (k, g) have 3–9 sources: Marson ×3 plus K562 / H1 / Orion / KOLF / Feng genome-wide. RPE1, HepG2 and Jurkat hold none of the final targets, and Tahoe holds one (MAPK7).

#### 0.5 SE distribution (per study, medians of per-table quantiles; `results/archive/2026-10-27/2026-10-27_shrinkage_study/tables/se_distribution.parquet` per table)
No table has a zero SE. NaN SE where the LFC is present is common in some tables (Orion, Feng, Tahoe, KOLF): genes with too few cells for an SE.

| Study | Tables | Entries with LFC | SE NaN | SE p05 / p50 / p95 |
|---|---|---|---|---|
| replogle | 2 | 87.3M | 21.6k | 0.083 / 0.235 / 0.555 |
| nadig | 2 | 32.9M | 23.7k | 0.103 / 0.287 / 0.677 |
| arc | 1 | 5.4M | 985k | 0.037 / 0.118 / 1.160 |
| adamson | 3 | 1.6M | 566k | 0.051 / 0.181 / 1.461 |
| orion | 4 | 587.2M | 146.4M | 0.068 / 0.206 / 1.448 |
| marson | 3 | 569.4M | 0 | 0.058 / 0.368 / 1.262 |
| kolf | 1 | 29.2M | 8.0M | 0.139 / 0.437 / 1.448 |
| feng | 21 | 152.2M | 41.4M | 0.106 / 0.340 / 1.512 |
| jiang | 30 | 23.7M | 5.1M | 0.091 / 0.312 / 1.540 |
| tahoe | 50 | 123.9M | 35.2M | 0.127 / 0.406 / 1.470 |

### Step 1 — baseline (pooled LFC MSE over the target's pairs)

| Line | MSE | strong | mid | weak | Predict zero |
|---|---|---|---|---|---|
| K562 | 0.1599 | 0.1489 | 0.1184 | 0.2124 | 0.1439 |
| HCT116 | 0.3986 | 0.3604 | 0.5078 | 0.3272 | 0.3742 |
| HEK293T | 0.1804 | 0.2120 | 0.1748 | 0.1548 | 0.1549 |
| H1 | 0.0514 | 0.0587 | 0.0452 | 0.0505 | 0.0345 |

Terciles are by strength_k, the RMS of the baseline p on the target, fixed per target. Harness metrics are below.

### Step 2 — calibration diagnostics (figures `step2a_hexbin_*`, `step2b_binned_*`, `step2e_4B_*`)

| Line | OLS₀ | Huber₀ | OLS slope / intercept | Slope \|p\| ≤ p90 / > p90 | Shape | Top-10% knockdowns' share of Σp² |
|---|---|---|---|---|---|---|
| K562 | 0.053 | 0.036 | 0.050 / −0.022 | 0.049 / 0.055 | linear | 26% |
| HCT116 | 0.034 | 0.022 | 0.028 / −0.034 | 0.038 / 0.031 | linear | 19% |
| HEK293T | 0.018 | 0.016 | 0.015 / −0.023 | 0.025 / 0.014 | linear | 20% |
| H1 | 0.113 | 0.063 | 0.114 / −0.004 | 0.109 / 0.115 | linear | 18% |

2c (slope by |S_kg| bucket) is in `results/archive/2026-10-27/2026-10-27_shrinkage_study/tables/step2c_buckets_*.parquet` and in the direct answers above. The negative intercepts mean the targets' mean LFC over these genes is below zero; the baseline has no such offset.

### Step 3 — source error statistics

**3a.** Few sources qualify: at least 20,000 shared pairs and at least 10 shared knockdowns *with a training target's final list*. That is 19 / 11 / 22 source–target rows per fold, from Marson ×3, Arc H1, KOLF2.1J, Feng genome-wide ×2, Orion ×2 and Replogle K562.
- **Scale:** per-source through-origin slopes of y on x are tiny, 0.003–0.06, and r is 0.001–0.07. No source transfers with a usable scale.
- **Marson has t < 0 in every row:** mean se² (0.07–0.23) exceeds V. Its replicate-table SEs overstate its error. All Marson rows are floored, which accounts for 6 / 6 / 9 of the floored rows.
- Tables: `results/archive/2026-10-27/2026-10-27_shrinkage_study/tables/step3a_sources_<fold>.parquet`.

**3b** (`results/archive/2026-10-27/2026-10-27_shrinkage_study/figures/step3b_<fold>.png`):

| Fold | n | R² | intercept | similarity | log cells | feng | Residual mean by target |
|---|---|---|---|---|---|---|---|
| K562 | 19 | 0.72 | 52.0 (19.0) | −36.7 (17.4) | −4.39 (0.83) | −4.65 (1.47) | H1 −0.30, HCT116 +0.15, HEK293T +0.07 |
| HCT116/293T | 11 | 0.79 | 46.5 (21.0) | −34.3 (21.4) | −4.01 (0.97) | −2.62 (2.38) | H1 −0.21, K562 +0.12 |
| H1 | 22 | 0.99 | 1.7 (4.1) | +13.1 (3.5) | −3.72 (0.19) | −4.25 (0.26) | ≈ 0 |
| Final (Step 5) | 33 | 0.68 | 28.6 | −13.0 | −4.10 | −4.85 | — |

The fit sees only "other" and "feng" rows. Jiang and Tahoe get their t̂ at the reference level (Jose 14:35; flagged).

**3c** (`results/archive/2026-10-27/2026-10-27_shrinkage_study/figures/step3c_<fold>.png`, `results/archive/2026-10-27/2026-10-27_shrinkage_study/tables/step3c_pairs_<fold>.parquet`):

| Fold | Source pairs | ρ_w | ρ_b | Within study | n_eff at S = 2 / 5 / 10 / 30 |
|---|---|---|---|---|---|
| K562 | 473 | 0.050 | 0.006 | feng 0.043, jiang 0.027, marson 0.304, orion 0.009 | 1.99 / 4.88 / 9.47 / 25.4 |
| HCT116/293T | 225 | 0.036 | 0.001 | feng 0.024, jiang 0.017, marson 0.304 | 2.00 / 4.99 / 9.94 / 29.4 |
| H1 | 771 | 0.040 | 0.004 | feng 0.036, jiang 0.022, marson 0.296, orion 0.013 | 1.99 / 4.92 / 9.65 / 26.9 |

No ρ was undefined or outside [−1, 1].

**3d** (`results/archive/2026-10-27/2026-10-27_shrinkage_study/figures/step3d_averaging.png`), equal-weight averaging over knockdowns with at least 10 sources:

| Line | Knockdowns | MSE at S = 1 → 10 | A (plateau) | B | Implied ρ = A / (A + B) |
|---|---|---|---|---|---|
| K562 | 13 | 0.277 → 0.189 | 0.181 | 0.096 | 0.65 |
| HCT116 | 14 | 0.625 → 0.548 | 0.555 | 0.085 | 0.87 |
| HEK293T | 14 | 0.269 → 0.105 | 0.097 | 0.180 | 0.35 |
| H1 | 1 | 0.231 → 0.062 | 0.054 | 0.184 | 0.23 |

### Step 4 — revised predictor (pooled LFC MSE on F; tercile columns strong / mid / weak)

| Line | Baseline | A | B | C | Zero | Slope y on p: base / A / B / C |
|---|---|---|---|---|---|---|
| K562 | 0.1599 (0.149 / 0.118 / 0.212) | 0.1601 | 0.1452 | **0.1438** (0.121 / 0.105 / 0.205) | 0.1439 | 0.05 / 0.09 / 0.25 / 1.76 |
| HCT116 | 0.3986 (0.360 / 0.508 / 0.327) | 0.4050 | 0.3786 | **0.3738** (0.323 / 0.484 / 0.314) | 0.3742 | 0.03 / 0.10 / 0.22 / 1.59 |
| HEK293T | 0.1804 (0.212 / 0.175 / 0.155) | 0.1923 | 0.1617 | **0.1548** (0.173 / 0.151 / 0.141) | 0.1549 | 0.02 / 0.06 / 0.13 / 0.97 |
| H1 | 0.0514 (0.059 / 0.045 / 0.051) | 0.0628 | 0.0357 | **0.0343** (0.037 / 0.029 / 0.037) | 0.0345 | 0.11 / 0.08 / 0.32 / 1.82 |
| Mean across folds* | 0.1669 | 0.1738 | 0.1503 | 0.1475 | 0.1477 | |

\*HCT116/293T counts once, as the mean of its two lines.

Per fold:

| Fold | c | c per training target | c per tercile (training targets) | B_k p10–p90 on F | Weights clipped | t̂ at reference level |
|---|---|---|---|---|---|---|
| K562 | 0.141 | HCT116 0.175, HEK293T 0.108, H1 0.279 | 0.10–0.33 | 0.15–0.45 | 1 / 1,485 | 2.6% of F's (k, g) (15 Jiang tables) |
| HCT116/293T | 0.136 | K562 0.131, H1 0.165 | 0.05–0.44 | 0.23–0.57 | 3 / 1,318 and 3 / 1,294 | 3.3% / 3.4% (66 tables, Jiang + Tahoe) |
| H1 | 0.174 | K562 0.244, HCT116 0.214, HEK293T 0.129 | 0.11–0.29 | 0.16–0.39 | 16 / 2,565 | 3.5% (6 tables) |

Single-source knockdowns on F: K562 2, HCT116 9, HEK293T 8, H1 0. B_k vs the oracle β_k: Spearman 0.31 / 0.24 / 0.21 / 0.03 (`results/archive/2026-10-27/2026-10-27_shrinkage_study/figures/step2e_4B_*`). After C, the binned conditional mean (`results/archive/2026-10-27/2026-10-27_shrinkage_study/figures/step4_binned_C_*`) is flat around zero with a small positive slope.

### Harness metrics (baseline and Stage C)

(pending: 4 runs)

### Stability across folds

| Parameter | K562 | HCT116/293T | H1 | Final (Step 5) | Reliable for unseen lines? |
|---|---|---|---|---|---|
| c | 0.141 | 0.136 | 0.174 | **0.256** | Stable across folds. **The final c is not comparable:** in Step 5 each target excludes only its own line, so HCT116 is predicted with HEK293T (same lab) as a source and vice versa. That raises the agreement, and with it c. An unseen line has no sibling, so the fold value (~0.14–0.17) is the relevant one. |
| ρ_w / ρ_b | 0.050 / 0.006 | 0.036 / 0.001 | 0.040 / 0.004 | 0.042 / 0.002 | Stable, and ≈ independence. |
| 3b similarity | −36.7 | −34.3 | +13.1 | −13.0 | **No:** the sign flips, SEs are large, and n is small. |
| 3b log cells | −4.39 | −4.01 | −3.72 | −4.10 | Stable, but carried by the imputed and floored Marson rows. |
| 3b feng | −4.65 | −2.62 | −4.25 | −4.85 | Two tables (Feng genome-wide, one pooled line). |
| Decile medians of B | 0.39–0.52 | 0.37–0.55 | 0.26–0.47 | 0.26–0.39 | Moderately stable. |

### Recommendation for the 3 unseen lines

(pending harness; effect space alone says the following)
- **In effect space, C is the best of the four on every fold**, but only by matching predict-zero.
- **Its gain over the baseline is entirely shrinkage:** B_k (≈ 0.3) × c (≈ 0.15).
- **Stage A's similarity-modelled weights should not be used as built.** They are fitted from an unstable 3b and are worse than the baseline's weights on all four lines.

### What I would have done differently / what undermines the design
1. **"Final lists everywhere" starves Step 3.** Only about 8 tables share at least 10 final-target knockdowns with a training target. Feng's 19 targeted lines, Jiang, Tahoe, Nadig and RPE1 never enter 3a/3b, so 3b is a 4-parameter fit on 11–22 rows. Fitting Step 3 on the training targets' *full* tables (K562 9,245 knockdowns, H1 300, Orion 434/418) would make 3b estimable. The metrics would stay on the final lists.
2. **3b's cell-count effect is a Marson artefact.** Marson has no cells (its count is imputed at ~400), and its reported SEs exceed its total error, so its t is negative and floored. The cell-count coefficient mostly separates Marson from the rest.
3. **The 3c formula is not an error correlation.** Its numerator, mean(x_s x_s') − ½(mean(y x_s) + mean(y x_s')), and V_s = mean(x_s²) − mean(y x_s) both omit mean(y²) and the doubled cross term. When y is barely predictable (r ≈ 0.01), ρ reduces to the correlation between the two sources' x, i.e. source–source agreement.
   - The 3d plateau, by contrast, includes everything common to all sources: the target's own measurement noise and its idiosyncratic response.
   - So the two disagree (ρ_b 0.004 against an implied 0.23–0.87) because they measure different things, not because of a bug. With so little shared signal, the plateau is mostly the truth's own variance.
4. **One fold is a single study:** HCT116/293T is Orion only, and as a training target it is one lab's assay.
5. **The 3d estimates are thin:** few knockdowns reach 10 sources (13 / 14 / 14 / 1).
6. **Scale differences are not the problem; agreement is.** Per-source slopes are all small and similar (0.003–0.06), so per-source scale correction (3a) would not help.

---

### Task (as given, 2026-10-27 14:00)

## Task: Shrinkage and source-dependence study for the zero-shot transfer model (3-fold)

### Goal
Measure how well our current precision-weighted transfer model is calibrated, how dependent the source cell lines are, and whether three changes reduce LFC MSE:
(A) correlation-aware, similarity-modelled source weights
(B) per-knockdown shrinkage from source agreement
(C) one global calibration constant
Evaluate with 3 held-out folds: K562, HCT116/293T, H1. Then refit on all three for the unseen lines.

### Folds and leakage rules
- Folds are K562, HCT116/293T and H1, using the fold definitions already in the codebase. If a fold contains more than one line, every line in it is treated as the held-out line together.
- In fold F, the held-out line's perturbation data is used nowhere: not as a target, and not as a source in any parameter fit. Only its control cells are used (for similarity).
- In fold F, parameters are fit using the other two fold lines as training targets. When a training target U is the target, the source pool is all lines except U and F.
- When F itself is predicted, the source pool is all lines except F.
- RPE1 is an ordinary source line throughout.

### Ground rules
- Build exactly what is specified. Do not add "use X only if it scores better" logic.
- If something cannot be built as specified on the data, stop that step and report why. Do not improvise a substitute.
- Do not modify the cell generator or the eval harness. Both are frozen.
- Write outputs to a new folder in the project on the VM: report.md, figures/, and tables as parquet. Large intermediates go to the GCS bucket.

### Notation
- T = target line, k = knockdown, g = downstream gene, s = source line.
- x_skg = source LFC, se_skg = its SE. Both come from the effect tables, which hold LFC and SE as separate arrays.
- y_kg = target true LFC.
- S_kg = the sources that measure (k, g). S_k = the sources that measure knockdown k.
- Coverage is asymmetric across the 117 lines, so S_kg varies.
- p_kg = baseline prediction from the current precision-weighted transfer model.
- strength_k = RMS over genes of the baseline p_kg on the target. Define strong/mid/weak terciles once per target, from its baseline, and keep them fixed for every metric on that target.
- SE floor: wherever se² is used, floor it at the 5th percentile of that source's se².

### Step 0: Inventory
1. Locate the current transfer model and how it sets weights; any control-cell similarity measure; the effect tables; the fold definitions; the eval harness. Summarise the baseline weighting in the report.
2. Similarity: use the existing control-cell similarity measure. If the codebase has none, use cosine similarity between mean log-normalised control profiles over the top 2,000 HVGs shared across lines, and say so.
3. Dataset/study ID for each of the 117 lines: take it from table metadata and report the mapping. Datasets with fewer than 3 lines are grouped as "other" wherever dataset is used below.
4. Coverage for each fold line as target (sources = all lines except that fold): number of knockdowns that have truth and at least 1 source; number of (k, g) pairs in each |S_kg| bucket: {1, 2, 3–4, 5–9, 10–29, 30+}.
5. Per-source SE distribution, plus counts of zero or NaN SE.

### Step 1: Baseline
For each fold F, run the current model zero-shot on F. Report per fold: LFC MSE overall and by strength tercile; eval harness metrics.

### Step 2: Calibration diagnostics (no model changes; run on each fold line's baseline)
2a. Pooled slope of y on p: through the origin with OLS and with Huber (HuberRegressor, fit_intercept=False); also OLS with an intercept, reporting the intercept. Figure: hexbin log-density of y vs p, with the y=x line and the fitted line.
2b. Binned conditional mean: 50 quantile bins of p; mean y (95% CI) against mean p per bin. Separate through-origin slopes for |p| below and above the 90th percentile of |p|. State whether the curve is linear or S-shaped.
2c. Through-origin slope within each |S_kg| bucket, with counts.
2d. Leverage: the share of Σp² from the top 10% of knockdowns by strength.
2e. Per-knockdown slope β_k (OLS through origin over genes): distribution; scatter against strength_k and against agreement τ̂²_k (Step 4B); Spearman correlations.

### Step 3: Source error statistics
For every (fold F, training target U): 3 folds × 2 training targets = 6 runs; source pool = all lines except U and F. Only source–target pairs sharing ≥ 20,000 (k, g) pairs and ≥ 10 knockdowns. All means over shared (k, g).
3a. Per source: shared pairs; through-origin slope of y on x_s (scale diagnostic only); Pearson r; mean se²; V_s = mean(x_s²) − mean(y·x_s); t_s = V_s − mean(se_s²), floored at 1% of the median positive t_s (report how many floored).
3b. Per fold, OLS over source–target pairs pooled over the two training targets: log t_s,U ~ similarity(s, U) + log(median cells per perturbation in s) + dataset fixed effects; no target fixed effect. Report coefficients and R²; observed vs fitted coloured by target; residual mean by target; raw t_s against median cell count. This model supplies t̂_s,T from control-cell similarity.
3c. Error correlation net of the target, for source pairs sharing ≥ 20,000 pairs with each other and with U: ρ_ss' = [mean(x_s·x_s') − ½(mean(y·x_s) + mean(y·x_s'))] / sqrt(V_s·V_s'). Per fold, Fisher-z mean pooled over the two training targets, within datasets (ρ_w) and between (ρ_b). Figure: ρ against source–source control similarity. n_eff = S / (1 + (S−1)ρ_b) for S = 2, 5, 10, 30.
3d. Empirical averaging curve (descriptive; once per fold line as target, sources = all lines except that fold): knockdowns with ≥ 10 sources; for S ∈ {1, 2, 3, 5, 8, 10}, 20 random S-subsets, equal average, MSE against y; fit MSE = A + B/S; report plateau A and implied ρ; compare with 3c.

### Step 4: Revised predictor (per fold; fit on the two training targets, apply to F with sources = all lines except F; cumulative stages, metrics after each)
Stage A (weights), per knockdown over S_k: d_s = median over genes of floored se²_skg + t̂_s,T (3b); R_k with ρ_w for same-dataset pairs, ρ_b otherwise; Σ_k = D^½ R_k D^½; w ∝ Σ_k⁻¹·1; clip negatives to 0 and normalise (report fraction clipped); renormalise over present sources per gene (no re-solve); p^A_kg = Σ_s w_s x_skg.
Stage B (per-knockdown shrinkage from agreement; sources only, computed on the target predicted): τ̂²_k = pair-count-weighted mean over source pairs of mean_g(x_s·x_s'); σ̂²_k = mean over sources of mean_g(x_s²) − τ̂²_k; ρ̄_k = mean off-diagonal of R_k; n_eff,k = |S_k| / (1 + (|S_k|−1)ρ̄_k); B_k = clip(τ̂²_k / (τ̂²_k + σ̂²_k / n_eff,k), 0, 1); single-source knockdowns get the median B_k of their p^A strength decile, from the fold's two training targets; p^B = B_k · p^A. On F: Spearman and scatter of B_k against the oracle β_k (2e).
Stage C (global calibration): c = OLS through-origin slope of y on p^B pooled over the two training targets (each built with A–B, sources excluding it and F); p^C = c · p^B on F. Diagnostics only: c per training target, c within each strength tercile. Repeat 2b's figure for p^C on F.
Metrics on F for baseline, A, B, C: LFC MSE overall and by tercile; pooled slope of y on the prediction; eval harness metrics for the baseline and Stage C only.

### Step 5: Final fit for the unseen lines
Refit 3b, ρ_w, ρ_b, the decile medians and c with all three fold lines as training targets (each target: sources = all lines except that target). Save as a single frozen config. No predictions for the unseen lines in this task.

### Step 6: Report
Summary paragraph leading with the most important finding; direct answers (slope well below 1? rises with source count? linear or S-shaped? source dependence ρ_w, ρ_b, n_eff, plateau? transfer variance vs similarity, cell count, dataset? β_k variability and whether agreement predicts it); stage-wise metric tables per fold and the mean across folds; stability across folds (c, ρ_w, ρ_b, 3b coefficients), flagging parameters too unstable for the final config; recommendation for the 3 unseen lines; what I would have done differently and anything in the data that undermines the design.

---

# 4. 2026-10-27 to 10-28 — D/E/F incumbent: Model 2 (effects_v5, no floor) vs C+G on Eval A
*Originally `results/stage1/2026-10-27_incumbent/report.md` (title: "D/E/F incumbent: Model 2 (effects_v5, no floor) vs C+G on Eval A — report"). Its tables, logs and outputs stay in `results/stage1/2026-10-27_incumbent/`.*

Task (Jose, 2026-10-27 15:10): port the leaderboard Model 2 recipe to effects_v5 without the detection floor, score
it against C+G end to end on Eval A, recommend one, then (after Jose's choice) prepare the D/E/F file. Nothing uploaded.

### Step 1 — the Model 2 recipe (2026-10-27 15:20)

#### Two corrections to the task's premises
1. **No entry scored 0.145.** +0.145 is the leaderboard *gain* from v2 to v3 (−0.054 → +0.091).
   - The best entry is **v3: +0.091**, rank 749, submitted 2026-10-13. It is Model 2 with top-2,000 calls and a detection floor, from `vcc_2026_3`.
   - The `vcc_2026_2` recipe itself (v2, gated Model 2) scored −0.054.
2. **The leaderboard Model 2 never read effects_v4.**
   - effects_v4 was built on 10-18, after every Model 2 entry (v2 10-1x, v3 10-13, v3a 10-14, v4 10-15).
   - v3 read the panel source tables `lib/sources/<line>/effects_panel.h5ad` (from `scripts/build_panel_sources.py`: 575 targets × 18,533 genes, layers lfc / se).
   - **Those panel tables no longer exist:** archived to the bucket on 10-19, and deleted from the bucket on 10-21 (user-approved cleanup; `ARCHIVE_MANIFEST.csv` rows now point to deleted objects).
   - The effects_v4 tables do survive in `gs://vcc2026/vm/lib/effects_v4/`, but no submission used them.

#### The recipe (v3, entry rank 749)
- **Code:** notebook `vcc_2026_3.ipynb` at commit `c208334` ("the vcc_2026_3 that ranked 749"). Extracted to `results/stage1/2026-10-27_incumbent/recipe/`: the notebook, `vcc2_meta.py` (Model 2), `vcc3_generate.py` (call rule, floor, generator) and `vcc2_common.py`, exactly as written by that commit. The current `code/vcc_2026_3c/` copies were edited later.
- **Selected arm:** confirmed from `gs://vcc2026/lib/vcc_2026_3/val_final/selected_generator.json`: `model2_top2000_floor`, array `re`. Settings: `call=topn`, `top_n=2000`, `detect_floor=true`, `z_req=5.0`, `alpha=1.0`, `max_lfc=8`, `min_effect=0`, `z=0`, `n_cells=400`, `composition=projection`, `kd_clip=[0.01, 0.9]`, `default_kd=0.3`, `seed=2026`. H1 validation scaled overall was 0.1177.
- **Sources:** one per cell line, 8 lines: HCT116, HEK293T, K562_gwps, CD4_Stim48hr (the only T-cell vote), H1, RPE1, HepG2, Jurkat. H1 was held out for validation; for the leaderboard all 8 voted.
- **Pooling (Model 2, `vcc2_meta.py`):** random-effects meta-analysis.
  - Per gene, DerSimonian–Laird τ² is pooled over every target seen in at least 2 backgrounds.
  - Prediction: re = Σ w·lfc / Σ w, with w = 1 / (SE² + τ²) and re_se = sqrt(1/Σw).
  - **No shrinkage:** the average is the prediction. The knocked-down gene's own column never enters.
- **Ranking and top N (`call_top_n`):** per knockdown, keep the 2,000 genes with the largest |re| / re_se among candidates.
  - **Candidates:** genes above 5 CPM in the line's own controls, finite and nonzero, and never the knocked-down gene.
  - **All other genes are set to 0.**
- **Finishing (`calibrate`):**
  1. × α (1.0).
  2. **Detection floor:** each called gene becomes sign(Δ) × max(|Δ|, floor_g). floor_g is the smallest |log2 change| a 400-cell group would show against the line's controls at z = 5, from the controls' per-cell CPM mean and variance (bisection; unreachable → 8).
  3. Clip to ±8, zero outside the 5-CPM gate, zero the own gene. The knocked-down gene is applied separately by measured guide efficiency (`kd_clip`, `default_kd`).
- **Leaderboard cells:** `vcc3_generate.py submit`, the vcc_2026_2a generator (composition projection, metadata cell counts and library sizes, state-dependent overdispersed induction). **This is not the frozen `code/cellgen/cell_generator.py`,** which was written and frozen later (10-26). It then went through `vcc prep` and `vcc submit`.

#### Leaderboard entries (overall / pds / nmae / fidelity / reach / Jaccard; MSE always 0)

| Entry | Recipe | Overall |
|---|---|---|
| v1 | — | −0.270 |
| v2 | `vcc_2026_2`, gated Model 2 (\|Δ\| > 3 SE and > 0.1, α 2) | −0.054 |
| **v3** | **`vcc_2026_3`, Model 2 top-2,000 + floor (z 5, α 1)** | **+0.091** / 0.416 / 0.048 / −0.088 / 0.192 / −0.024 |
| v3a | Model 2 top-1,000 + floor | +0.087 |
| v4 | `vcc_2026_3c`, every gene at 0.75 × Δ, no calling | +0.070 |
| v5 | transformer, no PIE, every gene + floor | +0.064 |
| v6 | transformer with PIE | +0.057 |
| v7 | transformer Δ̄ + 2(Δ − Δ̄) + floor | +0.078 |

So the entry to port is **v3, with N = 2,000.**

#### What this means for Steps 2–3 (needs Jose)
- **Step 2:** the v3 → v5 table mapping is per line and exists for every line:

  | v3 source | v5 table |
  |---|---|
  | HCT116 | `orion_hct116` + `orion_hct116_pseudobulk` |
  | HEK293T | `orion_hek293t` + `orion_hek293t_pseudobulk` |
  | K562_gwps | `replogle_k562_gw` |
  | CD4_Stim48hr | `marson_cd4_stim48hr` |
  | H1 | `arc_h1` |
  | RPE1 | `replogle_rpe1` |
  | HepG2 | `nadig_hepg2` |
  | Jurkat | `nadig_jurkat` |

  For Orion, the single-cell and pseudobulk tables split one line's knockdowns disjointly, so both are needed to cover what the one panel table covered.
- **Step 3's "v4 version (as submitted)" cannot be built as specified.** The tables v3 submitted with are deleted, and effects_v4 was never submitted.
- **v3 set every non-called gene to 0.** The task asks for "predicted LFC for all genes, plus the call set". Keeping non-called genes' LFC would be a third change.

#### Decisions after Step 1 (Jose, 15:30)
- Port v3, but do not hard-code N and do not set a floor.
- No v4 comparison: the effects_v4 tables will not be recovered for it.
- **Revised task (supersedes the top-N version):** the two arms differ only in pooling.
  - Model 2's pooling on effects_v5.
  - **C+G's call rule:** BH on |LFC| / the target's 400-cell measurement SE.
  - Non-called genes handled as C+G handles them.
  - Same generation path.
  - Step 3 compares the two models before any cells are made. If they are nearly identical (r > 0.95, similar call counts), Step 4 is cancelled.
  - If the paired CI includes 0, recommend C+G.
- The first port (top-N, zero non-called) was stopped mid table load at 15:34; nothing from it is used.

#### Step 1, addendum: how C+G hands non-called genes to the generator
**Kept at their predicted LFC.**
- `stage1.prepare_generator_inputs` (called by `gen_cells` for every knockdown) zeroes only genes **no source measured** (n_sources = 0); those are never called.
- Every other gene keeps the shrunk posterior mean, called or not. `de_sig` only marks which genes the generator shapes to test significant.
- `p_sign` only orders a saved `reach_order` that the frozen generator never reads.
- So generating C+G through `run_line --lfc-npz` from its own arrays is the same as its model path. Both arms use `--lfc-npz`.

### Step 2 — the port (`code/vcc6/model2_port.py`)

| v3 source (`lib/sources/<line>`) | v5 table(s) |
|---|---|
| HCT116 | `orion_hct116` + `orion_hct116_pseudobulk` (disjoint knockdowns of one line) |
| HEK293T | `orion_hek293t` + `orion_hek293t_pseudobulk` |
| K562_gwps | `replogle_k562_gw` |
| CD4_Stim48hr | `marson_cd4_stim48hr` |
| H1 | `arc_h1` (held out on Eval A) |
| RPE1 | `replogle_rpe1` |
| HepG2 | `nadig_hepg2` |
| Jurkat | `nadig_jurkat` |

Every line has a v5 table.

**Kept from the recipe, verbatim:** `tau2_dl` and `combine` from `vcc2_meta.py` at c208334, with:
- one vote per line;
- DL τ² pooled over the panel's targets seen in ≥ 2 lines (own column out, floored at 0, median fallback);
- w = 1/(SE² + τ²);
- **no shrinkage.**

**Changed:**
- **Calls:** `stage1.predictive_calls` with C+G's measurement SE for the same controls (`se_measurement` from C+G's own model run), over the controls' 5-CPM gate, own gene excluded.
- **Non-called genes:** keep their Model 2 LFC. Unmeasured genes (no line) are 0.

**One analogue:** v3 pooled τ² over its 575-target panel (its challenge + validation targets). The port pools over the targets it predicts: D/E/F 300 ∪ Eval A 292 = 564.

**Where this differs from how the leaderboard cells were made:**
- Cells come from the frozen generator (`code/cellgen`, 10-26), not `vcc3_generate.py submit`. That one differed in composition projection, metadata library sizes, the knockdown applied by measured guide efficiency, and the detection floor.
- Calls come from C+G's BH rule, not the top 2,000.

#### Step 2 — run (15:35–15:46, exit 0 for both jobs)
- **C+G model:** as `run_line predict` builds it, on Eval A's context controls. Only `arc_h1` dropped (the only H1 table). SE curve r = 0.923; prior τ² = 1e-4 (the floor), v0 = 0.541. Median 486 calls.
- **Model 2 port:** 7 lines vote, H1 held out. Panel rows used per table:

  | Table | Rows |
  |---|---|
  | `orion_hct116` + pseudobulk | 287 + 226 |
  | `orion_hek293t` + pseudobulk | 278 + 222 |
  | `replogle_k562_gw` | 501 |
  | `marson_cd4_stim48hr` | 492 |
  | `replogle_rpe1` | 84 |
  | `nadig_hepg2` | 64 |
  | `nadig_jurkat` | 86 |

  τ² was estimated for 17,712 genes (4,605 at 0); the median is 0.0030, 30× C+G's floored 1e-4.

### Step 3 — pre-generation comparison (Eval A, 292 knockdowns, H1 held out; `results/stage1/2026-10-27_incumbent/evalA/step3_compare.json`)

| | Model 2 port | C+G |
|---|---|---|
| Calls per knockdown, median / mean (truth: 218 / 1,048) | 249.5 / 388.6 | 486.5 / 551.4 |
| Median \|LFC\| of called genes | 0.459 | 0.315 |
| Median \|LFC\|, all measured genes | 0.075 | 0.062 |

- **Per-knockdown Pearson r between the two models' LFCs:** median 0.285 (IQR 0.210–0.364) over genes any source measured; 0.353 (0.262–0.488) over measured genes in the 5-CPM gate.
- **Call-set Jaccard:** median 0.085 (IQR 0.049–0.177).
- **Not nearly identical** (r ≪ 0.95, call counts differ about 2×), so Step 4 proceeds.
- The port calls about half as many genes at larger |LFC|. The two pools differ in sources: 7 lines vs 116 tables, including Feng's 20 iPSC lines, KOLF, Jiang, Adamson, three Marson conditions and Tahoe. They also differ in τ² (0.0030 vs 1e-4) and in shrinkage (none vs Gaussian).

### Step 4 — end to end on Eval A

#### Step 4 — results (Eval A, 292 knockdowns, H1 held out, salt "", frozen generator → frozen harness → `boot_score.py`)

**Overall (95% CI, knockdown bootstrap):** Model 2 port **0.1842 [0.1703, 0.1982]**; C+G **0.2186 [0.2035, 0.2329]**.

| Component | Model 2 port, scaled (raw) | C+G, scaled (raw) | Paired diff M2 − C+G [95% CI] | P(M2 better) |
|---|---|---|---|---|
| pds_cosine | 0.663 (0.805) | 0.723 (0.834) | −0.061 [−0.113, −0.003] | 0.019 |
| lfc_nmae | −0.037 (0.971) | 0.040 (0.927) | −0.077 [−0.099, −0.057] | 0.000 |
| direction fidelity | 0.357 (0.360) | 0.380 (0.378) | −0.023 [−0.047, −0.000] | 0.025 |
| direction reach | 0.047 (0.164) | 0.070 (0.184) | −0.023 [−0.048, +0.003] | 0.036 |
| sig Jaccard | 0.075 (0.043) | 0.097 (0.052) | −0.022 [−0.046, −0.004] | 0.005 |
| MSE (capped norm) | 0 (2.129) | 0 (1.496) | 0 (both clamped) | — |
| **Overall** | **0.184** | **0.219** | **−0.034 [−0.048, −0.020]** | **0.000** |

**Calls per knockdown, median / mean** (reference, the harness on the real cells: 218 / 1,048):

| | Requested | Realised by the generator | Called by the harness on the predicted cells |
|---|---|---|---|
| Model 2 port | 249.5 / 388.6 | 271 / 469.4 | 275.5 / 476.1 |
| C+G | 486.5 / 551.4 | 513.5 / 581.4 | 520 / 585.8 |

**MSE:** the bundle's baseline raw MSE is 0.967 and its replicate raw MSE is 0.0055.
- Model 2 port: raw 2.129.
- C+G: raw 1.496.
- Both are worse than the baseline, so both clamp to 0.
- cell-eval2's own `score_metrics` still returns 0 with `clamp_low=-inf`, because its MSE policy floors there. The linear unclamped (u − b)/(r − b) is **−1.21** (Model 2 port) and **−0.55** (C+G).

**Reading:**
- C+G is ahead on every component. Every paired CI excludes 0 except direction reach (−0.023 [−0.048, +0.003], P(M2 better) 0.036).
- The port's call count is closer to the reference median (250 vs C+G's 487, truth 218). It still loses Jaccard, because its calls are less often the right genes.
- Its larger effects (median called |LFC| 0.46 vs 0.32) cost lfc_nmae and MSE.

### Step 5 — recommendation
**C+G for D/E/F.** Its paired advantage over the Model 2 port on Eval A is +0.034, 95% CI [+0.020, +0.048], which excludes 0. It is ahead on every component.
- It also has end-to-end results on RPE1 (+0.122) and HCT116final (0.190 [0.144, 0.235]).
- **Caveats:**
  - H1 has misranked candidates before (v4 and v5 on the leaderboard). Both of those were call-everything recipes against a top-N one. Here the two arms share the call rule and generator and differ only in pooling, the kind of difference H1 ranked correctly in v2 → v3.
  - Model 2's leaderboard evidence (v3, +0.091) was earned with a different generator, top-2,000 calls and a detection floor. None of those are in either arm here.

**What I would have done differently:**
1. **Checked the recipe's history before porting.** The task's premises (an entry at 0.145; effects_v4 as the input) took a step to correct. The top-N port was started under the old task and stopped mid-load.
3. **Separated pooling from source set.** The port pools 7 lines; C+G pools 116 tables (Feng, KOLF, Jiang, Adamson, three Marson conditions, Tahoe). Their predictions correlate only 0.29, so this comparison cannot say whether C+G wins on its random-effects shrinkage or on its extra sources. A C+G arm restricted to Model 2's 7 lines would separate the two. It would only matter if the margin needed explaining, not for the choice.

**Stopped here, waiting for Jose's choice before Step 6** (the D/E/F file).

---

# 5. 2026-10-28 — D/E/F final file, C+G (submitted 10:07 UTC)
*Originally `results/stage1/2026-10-28_DEF/report.md` (title: "D/E/F final file — C+G (2026-10-28)"). Its tables, logs and outputs stay in `results/stage1/2026-10-28_DEF/`.*

**Decision (Jose, 2026-10-28 ~05:25 UTC):** C+G is the D/E/F model (Eval A: +0.034 [+0.020, +0.048] over the Model 2 port; `results/stage1/2026-10-27_incumbent/report.md`). Run it on D/E/F **and upload it.**

**Model:** `code/vcc6/stage1.py`, arm C with the Gaussian shrink, through `run_line.py predict`.
- Each context's own controls (`data_final/extracted/context_<X>.h5ad`, 18,400 each).
- The 300 final targets (`results/target_lists/final_DEF_300.txt`, from `pert_counts.csv`).
- All 117 tables vote: D/E/F match no source line.
- Salt "", 400 cells per knockdown, the frozen generator (sha256 422f8a3d…).

**File:** `code/vcc6/def_file.py` strips the controls, sets `context`, and concatenates to 360,000 cells × 18,533 genes. It then runs `vcc prep` with the final `gene_names.csv` and `pert_counts.csv` (contexts D, E, F), checking the exit code. After that comes `vcc submit`.

### Log
- 05:27: context D started. Controls: 18,400, median panel library 20,236, 9,633 gated genes; 117 tables loading.
- **06:25: restarted with model and generation in separate processes** (`results/stage1/2026-10-28_DEF/logs/run_def2.sh`), as the Eval A runs did.

### Submitted (2026-10-28 10:07 UTC)
- **Merge:** 2 min, on disk.
- **`vcc prep`:** 09:47–10:03, exit 0. 360,000 cells × 18,533 genes (float32); contexts D / E / F = 120,000 each; targets verified against the official list; counts preserved. Output `results/stage1/2026-10-28_DEF/final/prediction_DEF.vcc`, 3.19 GB.
- **`vcc submit -m LHY`:** exit 0, 3.0 GiB uploaded in 4 min 06 s. **Entry `JFCghVNyHwfCOWhnHG6M`**, final phase ("counts toward final scoring"). Server scoring job `vcc-metrics-prod-jfcghvnyhwfcowhnhg6m`. Scores are hidden until the winners are announced.
- **Model:** C+G (stage1 arm C + Gaussian shrink), all 117 tables, each context's own controls, salt "", the frozen generator. Median calls requested / realised: D 242.5 / 242.5, E 278.5 / 279, F 245 / 246.5.

---

# 6. 2026-10-28 to 10-30 — Scorer-space study (pseudobulk MSE, proxies, configs 1 and 1′, knockdown strength)
*Originally `results/stage1/2026-10-28_scorer_space/report.md` (title: "Scorer-space MSE study (pseudobulk only) — report"). Its tables, logs and outputs stay in `results/stage1/2026-10-28_scorer_space/`.*

## Part 1 — Steps 0–8 (task: Jose, 2026-10-28 ~06:05 UTC)

### Consolidated findings (Steps 0–8, 2026-10-28; detail in the step sections below)

**Recommendation (awaiting Jose's approval): arm (a), m = 0.5.** C+G's predicted LFC × 0.5 on every gene, with C+G's own call set (call rule untouched); κ off, w_arc = 1; the Flex weighting w is not fitted (blocked).
- **Fit:** m* = 0.5 on every left-out fold, so the spread max/min = 1.0 and m qualifies.
- **Why:** it is the only arm that takes Eval A's MSE above zero. Noise-free raw 1.496 → 0.853 against a 0.967 baseline, scaled 0 → +0.118. The cost is PDS 0.723 → 0.694.
- **Expected effect** from these two components alone:

  | Line | Overall change |
  |---|---|
  | Eval A | **+0.015** |
  | HCT116 | +0.014 |
  | HEK293T | −0.015 |
  | K562 | −0.014 |

  The four DE components need the end-to-end run. **The gain is about the size of the paired CI, so it may not pass the adoption rule.** On the weak lines the PDS loss (−0.09 to −0.25) is the main risk.

**Findings, in order of what they mean for the score:**
1. **The metric is fully understood and reproduced.** cell_eval2 0.18.0's MSE, rebuilt term by term and run through `mse_study.score_generated`, matches the harness to 9e-16 (HCT116 C+G, 2.896) and exactly (Eval A C+G, 1.496). PDS matches too.
   - The trace is right, with three refinements: per-knockdown averaging over \|G_p\|; the per-row credit cap; and a total budget that scales every credit by B/claim when it binds.
2. **The noise credit is not the problem.** The budget never binds for C+G. Restoring all lost credit would lower raw MSE by only 0.024 (Eval A) and 0.036 (HCT116). **The gap is in the plug-in distance.**
3. **A noise-free pseudobulk proxy predicts end-to-end MSE to 0.0006** (H1: 1.4956 vs 1.4962). Arms can be screened in seconds, without generating cells.
4. **Where C+G's error sits:**
   - **Non-called genes carry 76–89% of it.**
   - **Its largest LFCs sit on barely expressed genes:** K562's top 5% have a median target control CPM of 0.
   - **On Eval A, 6% of realised DE calls (median 29 per knockdown) were never requested.** 77% of those are in the top three expression deciles: kept small shifts on highly expressed genes turn significant.
5. **No single scaling removes more than 2% of the error** (R0² ≤ 0.02 in LFC, J²-weighted LFC and Δ space).
   - The calibration curve is a straight line with a small constant slope (Δ space: H1 ~0.14, others 0.02–0.07), not an S-curve.
   - The slope roughly doubles to quadruples where 4+ sources agree.
   - Beating the baseline on H1 needs R0² ≥ 0.047, which Δ-space arms never reach. Arm (a) gets there through the called genes and the shrink, not the non-called Δ.
6. **Scorer-space Δ pooling (arm c, m* 0.1) is the strongest arm overall but not on Eval A:**

   | Line | Unclamped scaled MSE (arm c, m* 0.1) |
   |---|---|
   | K562 | +0.18 |
   | HCT116 | +0.93 |
   | HEK293T | −0.43 |
   | H1 | −0.04 |

   It also keeps PDS higher, but leaves Eval A's MSE clamped at 0.
7. **The shared response is not the bottleneck:** arm (d)'s ceiling, the true mean Δ, gains ≤ 0.03. **The renormalisation offset** is under 1% of the error; κ is noise.
8. **Arc H1 upweighting hurts** where H1 has a say (Arc covers 28 of 300 D/E/F targets): ×20 turns +0.55…+2.2 into −0.3…−2.6 on the knockdowns it holds. H1's large effects transfer as magnitude error.
9. **Data limits:**
   - **Marson has no full profiles,** so 24–26% of the HCT116+HEK293T fold's tested (k, g) have no Δ source.
   - **Chemistry is recorded nowhere** (Jiang's note "Parse split-pool" excepted), so Step 7's family weighting is blocked pending a mapping from Jose.
   - **Orion's per-cell tables carry a stale "A–C only" note**; they actually cover 278 and 271 of the D/E/F targets.
   - **K562 and HEK293T references were built for this study** (harness cells + cell-eval2 MSE/PDS-only anchors). Their anchors are cell-eval2's own code; only the metric list was narrowed.
10. **Holding out both Orion lines together** (the study's fold) makes C+G much worse on them (HCT116 raw 7.0 vs 2.90 with HEK293T voting). For both lines, predicting the shared response alone beats C+G.

**What I would have done differently:**
- Used the noise-free proxy from the start.
- Checked metadata (chemistry, Marson profiles) before designing steps that needed it.
- Kept in mind that all arm numbers are noise-free with Ĉ_pred = 1.13 × C0, and the budget binds for arm (a) at m ≤ 0.5, so only the end-to-end run settles the recommendation.

**Next, after approval:** arm (a) m = 0.5 against the incumbent, end to end on Eval A and HCT116final, paired bootstrap, all six components. One job per VM.

---

Task: Jose, 2026-10-28 ~06:05 UTC. Rules: CPU and pseudobulk only through Step 5; no cell generation before Jose approves Step 6.

### Step 0 — the MSE path in the installed cell_eval2 0.18.0 (code reading, 06:05–06:15)
Paths are relative to `~/venv/lib/python3.12/site-packages/cell_eval2/`.

| Context item | Verdict | Where |
|---|---|---|
| Pseudobulk = group **sum** of counts, scaled to 50,000, then log1p | **Confirmed.** `bulk_lognorm_means`: `log1p(TS × P / Σ_g P)` on the per-group count sums (`_grouped_sums`); TS = `bulk_target_sum: 50000.0`. | `prep.py:329`, `prep.py:353–361`; `configs/vcc2026.yaml:60` |
| Numerator = ‖pred − real‖² with own target dropped, − C_real, − credit | **Confirmed, with two refinements.** (1) Per row it is divided by \|G_p\| (G − 1 when the target resolves to a measured gene), so each row is a per-gene mean. (2) The own-target coordinate is dropped from the plug-in distance only; both corrections stay whole. | `metrics/delta.py:959` (docstring), `:276` `_drop_on_target`, `:877–879` (`unbiased_sq_dist` → `_drop_on_target` → `_per_gene`) |
| C_real = jackknife of real cells | **Confirmed.** Delete-1 jackknife in bulk_lognorm space: `C_p = (n−1)/n Σ_g Σ_i (v_ig − v̄_g)²` with `v_ig = log1p(TS (P_pg − y_ig)/r_i)`. | `moments.py:160` `jackknife_correction` |
| credit = min(Ĉ_pred, C_real) per knockdown | **Confirmed**, with `PRED_TRACE_CAP_K = 1.0` (#247): `pred_term = min(pred_tn, k × real_tn)`. | `metrics/delta.py:43`, `:699` |
| Total credit capped at B, the across-knockdown centred spread of predicted means | **Confirmed, and refined.** B = Σ_g Σ_{p∈R_g} w_p (μ_pred[p,g] − wmean_{R_g})², with w_p = 1/\|G_p\| and each row's own target gene left out (#348). If B < Σ_p w_p × capped claim_p, every row's credit is scaled by **r = B / claim** (proportional allocation), not capped one by one. No tolerance multiplier. The budget is formed from the whole predicted panel. | `metrics/delta.py:157` `_across_pert_budget`; `:204` `_whole_panel_budget`; `:858–861` (`ratio = budget / claim; if ratio < 1: pred_term *= ratio`) |
| Denominator = ‖real − real_ctrl‖² − C_real − C_ctrl | **Confirmed** (`expr_distance_unbiased`, same own-gene drop and per-gene averaging). | `metrics/delta.py:1122` `distance_unbiased` |
| Ratio of sums; scaled vs baseline/replicate; clamped [0, 1] | **Confirmed.** `expr_mse_unbiased_capped_norm` is `agg="ratio_of_sums"` of `expr_mse_unbiased_capped` over `expr_distance_unbiased`. Scoring: `direction="lower", anchor=0.0, penalty="boxcox", clamp_low=0.0, clamp_high=1.0`; profiles include `vcc2026`. | `catalog.py:345–352` (components `:299`) |
| The vcc2026-test-1 config uses this metric | **Confirmed locally, not server-side.** `configs/vcc2026.yaml:10` sets `metrics: vcc2026`, and `expr_mse_unbiased_capped_norm` is in that profile; our harness runs report it. The final phase's `manifest.json` names the panel `vcc2026-test-1`, but the server's scoring config is not shipped. The vcc 0.3.0 CLI only reads `panel_id` from the server (`vcc/season.py:37,69`). | as cited |

**Correction to an earlier note of mine (10-28 05:20):** passing `clamp_low=-inf` to `score_metrics` still returned 0 for MSE. The metric's own Scoring has a Box-Cox penalty with `clamp_low=0.0`, so "unclamped" values in this study are computed directly as (u − b)/(r − b).

- **What we use:** the 10-21 MSE study, `analysis/mse_metric/mse_study.py` (with `scripts/v5/pb.py`). It **does implement both bounds** through cell_eval2's own internals:
  - per row, `min(C_pred, C_real)` (#247);
  - the #348 budget via `_across_pert_budget` / `_row_weights`, with `ratio = min(1, budget / claim)` (`mse_study.py:141–153`).

#### How C+G passes non-called genes to the generator
**Kept at their predicted LFC.** `stage1.prepare_generator_inputs` zeroes only genes no source measured (n_sources = 0); those are never called. Every other gene carries the shrunk posterior mean; `de_sig` only marks which genes get shaped to test significant.

#### How the generator turns LFCs into counts (`code/cellgen/cell_generator.py`, frozen, sha256 422f8a3d…)
- **Means:** mean CPM per gene = (control_mean_cpm + 1e-9) × 2^lfc − 1e-9 (`target_means`, line 168). This is the v5 LFC inverted, so the generated cells' mean CPM equals the requested one up to rounding.
- **Counts:** counts = round(cpm × library_size / 1e6), with largest-remainder rounding, so each gene's total count is exact.
  - `run_line` passes `library_size` = the controls' median panel library (`run_line.py:148`), e.g. 12,133 for HCT116final and 20,236 for context D.
- **Per-cell totals are not renormalized.** Each cell is built at `library_size` × (Σ requested CPM / 1e6). The generator's own note: "the scorer renormalizes each cell, so the realized means are the requested ones times 1e6 / sum(means)", and the realised composition is recorded in `uns['cellgen']['composition']`.
  - So when a prediction's means don't sum to 1e6, every gene's realised level shifts by the same factor after the scorer's per-cell (and per-group, CP50K) normalisation.
  - That shared shift is the "renormalization offset" Step 4 measures.

#### Reproduction (06:44, `code/vcc6/scorer_space.py check`; `results/stage1/2026-10-28_scorer_space/tables/step0_reproduce_HCT116.json`)
`scorer_space.Ref.parts` rebuilds the metric term by term from cell_eval2's own primitives. The real side comes from the eval's harness cache (`real_cache/*.moments.npz`). The predicted side comes from the cells (`prep.pseudobulk_bulk_lognorm_with_moments` + `moments.correction_for`). The budget is `delta._across_pert_budget` with weights 1/\|G_p\|.

| HCT116 C+G (14:48 run, 278 knockdowns) | Value |
|---|---|
| Rebuilt raw MSE | **2.8961771873171123** |
| Harness raw MSE | 2.896177187317113 |
| \|difference\| | **8.9e-16** (reproduces to 1e-6: yes) |
| Budget B vs claim Σ w·min(Ĉ_pred, C_real) | 0.550 vs 0.426, so r = 1 (does not bind) |
| Knockdowns whose credit the per-row cap limits (Ĉ_pred > C_real) | 16 of 278 |
| Bundle baseline / replicate raw MSE | 1.852 / 0.844 |
| Scaled MSE, unclamped (u − b)/(r − b) | −1.036 |

**Step 0 is complete.**

### Decision before Steps 1–5 (Jose, ~06:50)
**K562 and HEK293T get pseudobulk real sides built from their QC'd cells**, since only H1 (Eval A) and HCT116final have harness bundles.
- **Cells:** the harness's own build step (`final_evals.py build K562final|HEK293Tfinal --cells-only`): seeded 400 cells per knockdown, 20k-depth thinning, context/scoring control halves, on each line's final list (K562 267, HEK293T 271). No DE bundle.
- **Real side:** pseudobulk and jackknife in scorer space, computed with cell-eval2's own functions.
- **Anchors:** cell-eval2's own baseline and replicate-anchor code if it can be restricted to the MSE metric (no DE). Otherwise my own computation, labelled as such.
- **Scheduling:** runs after the D/E/F upload (the K562 cell file is 6.6 GB).

### Amended task (Jose, 2026-10-28 ~08:25): Steps 1–7
Every MSE number goes through `analysis/mse_metric/mse_study.py:score_generated`, with `code/vcc6/scorer_study.py:EvalPB` as its eval object (pb.py's attribute set on any harness eval; anchors and PDS via cell-eval2's own `_replicate_entries` / `score_one` / `discrimination_score`).

### Step 1 — reproduction through mse_study.py (`results/stage1/2026-10-28_scorer_space/tables/step1_reproduce.json`)

| | HCT116 C+G | Eval A C+G |
|---|---|---|
| raw MSE, rebuilt vs harness | 2.8961771873171123 vs 2.896177187317113 (diff 8.9e-16) | 1.4961872607668685 vs 1.4961872607668685 (diff 0) |
| PDS raw, rebuilt vs harness | 0.62784718 vs 0.62784718 (diff 1.1e-16) | 0.83353340 (diff 0) |
| knockdowns with Ĉ_pred > C_real (#247 cap) | 16 / 278 | **249 / 292** |
| budget B / claim (#348) | 1.29 (does not bind) | 3.91 (does not bind) |
| baseline / replicate raw MSE | 1.852 / 0.844 | 0.967 / 0.0055 |
| scaled MSE (clamped / unclamped) | 0 / −1.036 | 0 / −0.550 |

On HCT116 the credit actually granted, Σ min(Ĉ_pred, C_real)/|G|, is 6.3× the summed denominator: HCT116's real signal is tiny next to its sampling noise. On Eval A it is 0.27×.

### Step 2 — Δ per source (`scorer_study.py step2`; `results/stage1/2026-10-28_scorer_space/delta/<table>.npz`, `results/stage1/2026-10-28_scorer_space/tables/step2_*`)
Δ_src,k,g = `log1p_bulk_50k_full`[knockdown row] − `log1p_bulk_50k_full`[that source's control, `<source>_baseline.csv`].
- Written as separate files, one per table, for every knockdown any fold or D/E/F needs. The existing layers are untouched.
- Every table's gene order matches the panel.

**1. Which tables have a Δ.** 110 of 117 tables have both a finite knockdown full profile and a control profile.
- **No Δ: Marson CD4 ×3** (rest, stim 8 h, stim 48 h). They have a control profile, but the knockdown full profile is empty (replicate-level tables, no cells). They hold 492–502 of the needed knockdowns, including 278 of the 300 D/E/F targets.
- **Not needed:** `adamson_k562_10x001`, `adamson_k562_10x005`, `feng_targeted_fiaj_3` and `feng_targeted_pipw_5` hold none of the needed knockdowns.
- **Tahoe (50) and Jiang (30) do have Δ.** Their tables carry full profiles.

**2. (k, g) with no Δ source once the fold's lines are removed:**

| Fold | Line | Knockdowns | Tested (k, g) without Δ | Knockdowns with no Δ at all |
|---|---|---|---|---|
| K562 | K562 | 267 | 60,510 / 2,050,827 (3.0%) | 5 |
| HCT116+HEK293T | HCT116 | 278 | **685,872 / 2,868,960 (23.9%)** | 20 |
| HCT116+HEK293T | HEK293T | 271 | **774,830 / 2,998,886 (25.8%)** | 20 |
| H1 | H1 (Eval A) | 292 | 46,311 / 3,144,840 (1.5%) | 2 |

"Tested" means the target table's 5-CPM gate. Over all panel genes it's 3.4% / 40% / 40% / 2.6%. The HCT116+HEK293T hole is Marson: with both Orion lines out, Marson was the main remaining source for the final targets, and it has no Δ.

**3. Rest-of-line reference:** `feng_gw_fitness` and `feng_gw_nonfitness` (`qc_replogle_umi.REFERENCE = {'FengGW_fitness': 'rest', 'FengGW_nonfitness': 'rest'}`, `scripts/qc_replogle_umi.py:112`; the tables' obs `reference = rest_of_line`). No other table.

**4. Orion coverage of the 300 D/E/F targets:** `orion_hct116` **278**, `orion_hek293t` **271** (the per-cell tables since the 10-27 rebuild); the two pseudobulk complements 0.
- ⚠️ Both per-cell tables still carry the note "Orion: the A–C targets subset only" in `uns`. It is stale since the rebuild; the content covers the final targets.

**5. Panel genes measured, and chemistry:**
- **Genes measured** (median per study): Replogle 0.43 (K562 0.41, RPE1 0.45) and Nadig 0.47. Those releases ship only their expressed genes, about 8k. Everything else is 0.92–0.98: Arc 0.98, Orion 0.98, Feng 0.98, KOLF 0.98, Tahoe 0.98, Jiang 0.92, Marson 0.96, Adamson 0.96.
- **Chemistry is not recorded in the tables' metadata,** except Jiang (`depth_note`: "Parse split-pool, shallow (~3.9k UMI per cell)"). Tahoe's note gives depth only (~1.7k UMI per cell). Step 7 needs a chemistry label; I'll check the cell files' metadata there and stop that step if it isn't recorded, rather than guess.

**6. Sanity, H1 SIN3B → SERPINB9:**

| Quantity | Value |
|---|---|
| lfc | 1.468; × ln 2 = **1.018** |
| control mean CPM | 73.7 |
| full profile: knockdown / control | 2.424 / 1.582 |
| **Δ** | **0.841** |

The expectation was Δ ≈ 0.87. That's log((1 + v·2^lfc)/(1 + v)) with v = 73.7 × 0.05 = 3.68. The small gap is the full profile's own depth normalisation. So Δ is about 0.83× lfc·ln 2 here: log1p compresses moderately expressed genes.

### Step 7 — blocked for the chemistry families (08:42)
**No chemistry is recorded in the metadata:** not in the 117 effect tables' `uns`/`obs`, and not in any `lib/qc_cells/*_qc.h5ad` (`uns` searched recursively for chem/tech/platform/assay/kit/protocol/flex/library_type).
- The only hint is Jiang's depth note, "Parse split-pool".
- Per the rules this part stops rather than guess. It needs a source → chemistry mapping from Jose (e.g. from the publications).
- The Arc-first sub-arm (Arc H1's weight × w_arc) needs no chemistry and will run once Step 5 has the arm (c) m*.

### Progress (08:42)
Chain:
1. Step 4 (noise credit).
2. K562/HEK293T real caches + MSE/PDS-only anchors (HCT116final's exact cell-eval2 commands, metric list narrowed).
3. Per-fold C+G + Δ pool (one table load).
4. Steps 3, 3b, 5, 6.

The K562final cells are built (267 knockdowns, 48,020 cells, 18,400 + 18,400 controls); HEK293Tfinal is building.

### Step 4 — noise-credit accounting (`scorer_study.py step4`; `results/stage1/2026-10-28_scorer_space/tables/step4.json`, `results/stage1/2026-10-28_scorer_space/tables/step4_per_knockdown.parquet`)
- **Cells:** the already generated C+G cells (Eval A, 10-28; HCT116, 10-27 14:25).
- **C0:** cell_eval2's delete-1 jackknife on 400 control cells, 50 draws from each eval's scoring controls.
- **Units:** "credit" and "Q" are in score units, Σ_p w_p × term with w_p = 1/|G_p|, as the scorer forms them. The budget numbers come from cell_eval2's `_across_pert_budget` / `_row_weights`.

| | Eval A (H1, 292) | HCT116 (278) |
|---|---|---|
| C0 (mean ± sd over 50 draws) | 14.12 ± 0.28 | 24.85 ± 1.25 |
| Ĉ_pred median / C_real median | 15.71 / 14.29 | 28.56 / 68.64 |
| share of knockdowns with Ĉ_pred > C_real (#247 cap applies) | **85.3%** | 5.8% |
| Ĉ_pred / C0: median (p10–p90) | **1.113** (1.100–1.135) | **1.150** (1.125–1.170) |
| claim before cap → Q after cap | 0.248 → 0.228 | 0.428 → 0.426 |
| budget B | 0.893 (Q < B: **does not bind**) | 0.550 (**does not bind**) |
| credit lost to the cap / to the budget | 0.020 / 0 | 0.0025 / 0 |
| Σ denominator (for scale) | 0.834 | 0.067 |
| raw MSE actual | 1.4962 | 2.8962 |
| raw with the budget loss restored (Q − min(B, Q) granted) | 1.4962 (diff 0) | 2.8962 (diff 4e-16) |
| raw with the cap loss also restored (full Ĉ_pred credited) | 1.4723 (**−0.024**) | 2.8598 (**−0.036**) |

**Reading:**
- The budget never binds, so restoring the lost credit changes nothing.
- Even removing the per-knockdown cap would lower raw MSE by only 0.024 (Eval A) and 0.036 (HCT116). Beating the baseline on Eval A needs raw < 0.967.
- **The MSE gap is in the plug-in distance, not in the noise credit.**
- On HCT116 the real side's own correction (C_real median 68.6, from its smaller real groups) is 2.8× C0. The corrections total 6× the denominator there, which is why the HCT116 scale is so fragile.
- Step 5 uses Ĉ_pred = C0 × **1.130**, the median Ĉ_pred/C0 pooled over both lines' 570 knockdowns (`scorer_arms.ratio_from_step4`).

### K562 and HEK293T reference sides (08:36–08:56)
- **Cells:** `final_evals.py build K562final|HEK293Tfinal --cells-only`, the harness's own build step. Seeded 400 cells per knockdown, 20k-depth thinning (factor 1.0 for both), context and scoring control halves.
  - K562final: 267 knockdowns, 48,020 real cells, 18,400 + 18,400 controls.
  - HEK293Tfinal: 271 knockdowns.
- **Real side and anchors:** cell-eval2 itself (`results/stage1/2026-10-28_scorer_space/logs/scorer_refs.sh`): HCT116final's exact `prep-cache` / `baseline` / `prep-real-bundle` commands (vcc2026 preset, pdex, counts, 5 anchor splits). The only change is `--set metrics=[expr_mse_unbiased_capped, expr_distance_unbiased, expr_mse_unbiased_capped_norm, pds_cosine]`, so no DE is computed (cache 48 s instead of hours). Every exit code was 0.

| Fold line | Baseline raw MSE | Replicate raw MSE | Baseline PDS | Replicate PDS | Source |
|---|---|---|---|---|---|
| K562 | 1.095 | −0.076 | 0.515 | 0.869 | built today (MSE/PDS only) |
| HCT116 | 1.852 | 0.844 | — | — | HCT116final bundle |
| HEK293T | 1.141 | −0.059 | 0.501 | 0.768 | built today (MSE/PDS only) |
| H1 (Eval A) | 0.967 | 0.0055 | — | — | Eval A bundle |

### Per-fold models (`scorer_study.py folds`; `results/stage1/2026-10-28_scorer_space/cache/fold_<line>.npz`)
- **One table load for all folds.** Per fold, C+G is rebuilt exactly as `run_line`/`predict_line` builds it: arm C votes, SE curve, prior, Gaussian shrink, on each eval's context controls, with the fold's lines' tables dropped. The Δ pool uses the same prior, SE curve and votes, with each source's LFC replaced by its Δ and duplicates collapsed as `load_source` does.
- **Identity check:** the H1 fold's C+G equals the incumbent Eval A C+G model (`results/stage1/2026-10-27_incumbent/evalA/CG/model`) exactly: max |ΔLFC| 0, calls identical, measurement SE identical.

| Fold | Tables voting / with Δ | Prior τ² / v0 | Line: median calls, share of (k, g) with a Δ source |
|---|---|---|---|
| K562 | 108 / 103 | 1e-4 / 0.783 | K562: 113, 0.79 |
| HCT116+HEK293T | 113 / 106 | 1e-4 / 0.499 | HCT116: 461, 0.52; HEK293T: 592, 0.52 |
| H1 | 116 / 109 | 1e-4 / 0.541 | H1: 486, 0.81 |

### Step 3 — where the MSE error sits (`scorer_arms.py step3`, `step3b`; `results/stage1/2026-10-28_scorer_space/tables/step3.json`, `results/stage1/2026-10-28_scorer_space/tables/step3b.json`)
**Setup.**
- C+G per fold (fold pool), as a **noise-free pseudobulk**: the generator's mean inversion (control + 1e-9)·2^lfc − 1e-9, then the scorer's log1p(5e4·μ/Σμ). Ĉ_pred = C0 × 1.1295 (Step 4).
- Scored by `mse_study.score_generated`.
- **Proxy check:** H1 noise-free C+G raw **1.4956** against the generated cells' **1.4962**. The pseudobulk proxy tracks end-to-end MSE to 0.0006.

**Raw MSE (unclamped scaled) and anchors:**

| Fold line | C+G | (i) predict control | (ii) shared response s̄ alone | Baseline / replicate raw | PDS scaled, C+G |
|---|---|---|---|---|---|
| K562 | **1.057 (+0.033)** | 1.216 (−0.103) | 1.132 (−0.032) | 1.095 / −0.076 | 0.511 |
| HCT116 | 7.005 (−5.11) | 1.810 (+0.041) | **1.341 (+0.507)** | 1.852 / 0.844 | 0.528 |
| HEK293T | 5.320 (−3.48) | 1.653 (−0.427) | 1.136 (+0.004) | 1.141 / −0.059 | 0.149 |
| H1 (Eval A) | 1.496 (−0.550) | 1.014 (−0.049) | 1.005 (−0.040) | 0.967 / 0.0055 | 0.724 |

- Clamped scaled MSE is max(0, unclamped): only K562 C+G (0.033) and HCT116's control/s̄ arms are above 0.
- s̄ = the mean over the target's knockdowns of the pooled Δ (Step 5d).
- In the HCT116+HEK293T fold both Orion lines are held out, so C+G there is much worse than the screen-setup HCT116 run (2.90, where HEK293T voted). Predicting nothing beats it.

**Split of C+G's plug-in error ‖P − R‖² (own gene out):**

| Fold line | Called genes | Non-called genes | Ungated genes | Strong / mid / weak tercile | (C_real + credit) / plug-in |
|---|---|---|---|---|---|
| K562 | 11% | **89%** | 0% | 39 / 31 / 30% | 0.92 |
| HCT116 | 18% | **82%** | 5% | 39 / 35 / 26% | 0.78 |
| HEK293T | 24% | **76%** | 4% | 41 / 34 / 25% | 0.61 |
| H1 | 23% | **78%** | 5% | 48 / 30 / 23% | 0.28 |

By target control-CPM decile (gated genes, low → high), the error is spread fairly evenly. It peaks in the middle deciles for HCT116/HEK293T/H1 and in the low deciles for K562 (`results/stage1/2026-10-28_scorer_space/tables/step3.json`).

**Top 5% of (k, g) by |C+G LFC| (step3b):** median control CPM in the target vs in the sources that hold k:

| | K562 | HCT116 | HEK293T | H1 |
|---|---|---|---|---|
| Target | **0.0** (only 7% in K562's gate) | 15.6 | 20.7 | 4.7 |
| Sources | 0.56 | 18.1 | 18.2 | 5.9 |

C+G's largest effects sit on barely expressed genes. In K562 they are genes K562 doesn't express, so they carry almost no weight in scorer space.

**Realised DE calls not in the requested set** (harness DE on the generated cells; p_adj < 0.05, own gene out):

| | Realised calls | Not requested | Share | Per knockdown, median / mean | Where (target control-CPM decile) |
|---|---|---|---|---|---|
| HCT116 C+G | 53,219 | 147 | 0.3% | 0 / 0.5 | low deciles (0–2: 122) |
| Eval A C+G | 171,056 | 10,650 | **6.2%** | 29 / 37 | **top deciles: 7–9 = 8,236 (77%)** |

On H1 the non-called genes keep their small predicted shifts, and at high expression 400 cells make those shifts significant against the real controls. That's the generator's "not DE" shaping being judged against a different control set.

### Step 5 — magnitude, four arms (`scorer_arms.py step5`; `results/stage1/2026-10-28_scorer_space/tables/step5_grid.parquet`, `results/stage1/2026-10-28_scorer_space/tables/step5.json`)
**Setup.**
- Noise-free pseudobulks scored by `mse_study.score_generated`, with Ĉ_pred = C0 × 1.1295.
- Grid m ∈ {0, 0.1, 0.2, 0.3, 0.5, 0.75, 1}; arm (d) also a ∈ {0, 0.5, 1, 1.5}.
- **The arms:**
  - (a) m × LFC on all genes.
  - (b) m × LFC on non-called genes only.
  - (c) non-called genes in Δ space: ctrl_profile + m·Δ_pool; called genes at C+G's LFC; genes with no Δ source at control.
  - (d) like (c), with Δ = a·s̄ + m·(Δ_pool − s̄).

**Unclamped scaled MSE over the grid** (PDS scaled in brackets at m = 0 / 0.5 / 1):

| Arm | Line | m = 0 | 0.1 | 0.2 | 0.3 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|---|---|---|
| a | K562 | −0.103 | −0.068 | −0.035 | −0.002 | **+0.058** | **+0.127** | +0.033 |
| a | HCT116 | +0.041 | +0.105 | +0.166 | +0.224 | **+0.332** | −0.041 | −5.11 |
| a | HEK293T | −0.427 | −0.405 | −0.384 | −0.365 | −0.330 | −1.27 | −3.48 |
| a | H1 | −0.049 | −0.015 | +0.019 | +0.052 | **+0.118** | −0.141 | −0.550 |
| c | K562 | +0.168 | +0.184 | +0.194 | **+0.196** | +0.179 | +0.119 | −1.14 |
| c | HCT116 | +0.903 | +0.933 | +0.960 | **+0.983** | +0.218 | −2.14 | −5.45 |
| c | HEK293T | −0.406 | −0.425 | −0.501 | −0.634 | −1.07 | −1.92 | −3.12 |
| c | H1 | −0.053 | −0.040 | −0.045 | −0.069 | −0.168 | −0.391 | −0.720 |

- **PDS:**
  - Arm (a): H1 0.013 / 0.694 / 0.723; K562 −0.04 / 0.40 / 0.51; HCT116 −0.12 / 0.28 / 0.53.
  - Arm (c) keeps PDS high: H1 0.72 / 0.78 / 0.76; HCT116 0.74 / 0.92 / 0.95; K562 0.52 / 0.64 / 0.63.
- **Arm (b)** (m on non-called only, called kept) is worse than (a) everywhere except K562 at m = 0.75 (+0.154). Full table in `step5_grid.parquet`.
- **The budget B binds** for arm (a) at every m ≤ 0.5 on all four lines, and for (b)/(c) on 2 lines at m ≤ 0.3. Shrinking all effects shrinks the across-knockdown spread that funds the noise credit. It never binds at m = 1.

**m* fitted on the other two folds (pooled numerator) and applied to the third:**

| Arm | Fold K562: m* → K562 | Fold HCT116+HEK293T: m* → HCT116 / HEK293T | Fold H1: m* → H1 (Eval A) |
|---|---|---|---|
| a | 0.5 → +0.058 (raw 1.027), PDS 0.400 | 0.5 → +0.332 / −0.330, PDS 0.282 / 0.061 | 0.5 → **+0.118 (raw 0.853)**, PDS 0.694 |
| b | 0.1 → +0.057 | 0.2 → +0.398 / −0.915 | 0.1 → −0.051 |
| c | 0.1 → **+0.184** (raw 0.879), PDS 0.576 | 0.1 → **+0.933** / −0.425, PDS 0.791 / 0.416 | 0.1 → −0.040 (raw 1.005), PDS 0.761 |
| d | (a 0.5, m 0.1) → +0.172 | (0.5, 0.1) → +0.910 / −0.434 | (0, 0.1) → −0.042 |
| d, ceiling (true mean Δ in place of s̄) | (0.5, 0.1) → +0.203 | (1.0, 0.1) → +0.916 / −0.435 | (0.5, 0.1) → −0.017 |

- **m* is perfectly stable across folds within an arm:** (a) 0.5, 0.5, 0.5; (c) 0.1, 0.1, 0.1. So max/min = 1 for both.
- **Only arm (a) at m* = 0.5 takes Eval A above zero:** raw 0.853 against a 0.967 baseline. Its cost is PDS, 0.694 against C+G's 0.723.
- **The ceiling barely helps arm (d):** knowing the true mean shift instead of s̄ gains ≤ 0.03. So the shared response is not where the Δ-space error is.

**Calibration** (gated non-called genes for Δ; true LFC from the line's table for LFC space):

| Line | β LFC unweighted (R0²) | β LFC J²-weighted (R0²) | β Δ space (R0²) | R0² needed to beat the baseline: 1 − baseline_raw / control_raw |
|---|---|---|---|---|
| K562 | 0.053 (0.000) | 0.056 (0.001) | 0.046 (0.001) | 0.099 |
| HCT116 | 0.033 (0.000) | 0.031 (0.000) | 0.022 (0.000) | −0.023 (control already beats it) |
| HEK293T | 0.018 (0.000) | 0.021 (0.000) | 0.016 (0.000) | 0.31 |
| H1 | 0.154 (0.014) | 0.156 (0.020) | 0.143 (0.013) | 0.047 |

- **No single scaling can remove more than 2% of the error on any line** (R0² ≤ 0.02). Arms (c)/(d) in Δ space never reach the R0² needed: H1 0.013 vs 0.047, K562 0.001 vs 0.099.
- So arm (a)'s gain on H1 does not come from the non-called genes' Δ. It comes from shrinking the called genes' large LFCs, and from the budget, which binds at m = 0.5.
- **Binned calibration in Δ space** (deciles of |Δ_pool|): **a straight line with a small, constant slope.** It's not flat near 0 and steeper in the tails. Per-decile β: H1 0.14–0.19, K562 0.04–0.10, HCT116 0.02–0.14, HEK293T 0.01–0.06. The true Δ also carries a **constant negative offset** in every bin (about −0.002 to −0.003 for K562/HCT116/HEK293T, about −0.001 for H1): a renormalisation shift (Step 6).
- **β in Δ space by strength tercile** (strong / mid / weak): K562 0.027 / 0.056 / 0.091; HCT116 0.021 / 0.024 / 0.021; HEK293T 0.019 / 0.015 / 0.010; H1 0.177 / 0.097 / 0.133.
- **β in Δ space by number of sources** (1 / 2–3 / 4+): K562 0.023 / 0.042 / **0.125**; HCT116 0.013 / 0.036 / 0.081; HEK293T 0.011 / 0.030 / 0.005; H1 0.063 / 0.055 / **0.232**. Agreement among 4+ sources roughly doubles to quadruples the usable slope.

### Step 6 — renormalisation offset (`scorer_arms.py step6`; `results/stage1/2026-10-28_scorer_space/tables/step6.json`)
**Setup.**
- **Arm:** the Step 5 arm with the best pooled unclamped MSE over the held-out lines, **arm (c) at m* = 0.1** (mean +0.163, against +0.045 for arm (a)).
- **Offset:** per knockdown, the median scorer-space delta over gated non-called genes, P − control bulk for the prediction and R − real control for the truth.
- **Share:** n_non-called × (gap of the medians)² over the plug-in total.
- **κ:** one constant per knockdown, iterated through the renormalisation (5 passes, residual |median| ≤ 0.0022).

| Line | Median delta, predicted (p10 / p50 / p90) | True (p10 / p50 / p90) | Offset share of the plug-in error | MSE unclamped: before → with κ | PDS scaled: before → with κ |
|---|---|---|---|---|---|
| K562 | −0.009 / −0.003 / +0.004 | −0.009 / −0.001 / +0.006 | 0.73% | +0.184 → +0.165 | 0.576 → 0.604 |
| HCT116 | +0.001 / +0.003 / +0.005 | −0.010 / −0.002 / +0.006 | 0.71% | +0.933 → +0.929 | 0.791 → 0.768 |
| HEK293T | +0.001 / +0.003 / +0.004 | −0.011 / −0.003 / +0.006 | 0.95% | −0.425 → −0.402 | 0.416 → 0.408 |
| H1 | −0.002 / +0.001 / +0.004 | −0.005 / −0.000 / +0.006 | 0.53% | −0.040 → −0.036 | 0.761 → 0.763 |

**The renormalisation offset is under 1% of the plug-in error on every line**, and κ moves unclamped MSE by at most 0.02, in both directions. It is not a lever. The constant negative offset in the binned true Δ (Step 5) is real but tiny next to the per-gene error.

### Step 7 — Arc-first sub-arm (`scorer_study.py arcw`, `scorer_arms.py step7arc`; `results/stage1/2026-10-28_scorer_space/tables/step7_arc.*`)
**Arc H1 covers 28 of the 300 D/E/F targets.**
- **Setup:** arm (c) at its m* (0.1). Arc H1's pooling weight × w_arc in both pools (C+G's LFC and the Δ pool): the H1 background's vote variance divided by w_arc inside the same random-effects pool, with w (Flex) = 1. Folds: K562 and HCT116+HEK293T, where H1 is a source.
- **Columns:** unclamped scaled MSE / PDS scaled / J²-weighted r in LFC space / J²-weighted r in Δ space.

| Line, scope | w_arc 1 | 2 | 5 | 20 |
|---|---|---|---|---|
| K562, all 267 | +0.184 / 0.576 / 0.034 / 0.045 | +0.188 / 0.581 | +0.190 / 0.582 | +0.191 / 0.575 / 0.034 / 0.045 |
| K562, the 27 Arc holds | **+0.554** / 0.655 / 0.075 / 0.100 | +0.425 | +0.163 | **−0.319** / 0.638 / 0.064 / 0.099 |
| HCT116, all 278 | +0.933 / 0.791 / 0.022 / 0.033 | +0.940 | +0.946 | +0.948 / 0.764 |
| HCT116, the 27 Arc holds | **+2.213** / 1.32\* / 0.066 / 0.092 | +1.712 | +0.861 | **−0.361** / 1.05\* / 0.056 / 0.096 |
| HEK293T, all 271 | −0.425 / 0.416 / 0.019 / 0.029 | −0.479 | −0.568 | −0.692 / 0.410 |
| HEK293T, the 25 Arc holds | **+0.544** / 0.528 / 0.062 / 0.065 | −0.093 | −1.133 | **−2.592** / 0.464 / 0.051 / 0.073 |

\* PDS scaled on a subset can exceed 1, because the anchors are whole-panel values.

**Reading:**
- Upweighting Arc H1 makes MSE worse exactly where H1 has a say (the 25–27 knockdowns it holds), and leaves J²-weighted r flat or lower.
- Over all knockdowns the effect is under ±0.02, except HEK293T, which gets worse.
- H1's large, consistent effects transfer as magnitude error, not signal. **w_arc stays at 1.**

**The chemistry-family part of Step 7 remains blocked:** no chemistry metadata exists (see above).

### Step 8 — recommendation (stop for Jose's approval)
**Config: arm (a), m = 0.5. κ off, w (Flex) = 1 (blocked, not fitted), w_arc = 1.**
- **What it is:** C+G's predicted LFC × 0.5 on every gene, called and non-called, handed to the frozen generator.
- **Calls unchanged:** C+G's own call set; the call rule is untouched.
- **Fitted on training folds only:** m* = 0.5 when fitted on any two folds and applied to the third (K562, HCT116+HEK293T, H1). **Spread max/min = 1.0 ≤ 2, so m qualifies.**

**Why this arm and not (c):**
- Only arm (a) takes **Eval A's MSE above zero**: noise-free raw 1.496 → **0.853**, below the 0.967 baseline, so scaled +0.118.
- Arm (c) (m* 0.1) is better on K562 and HCT116 and keeps PDS higher, but leaves Eval A at −0.040, clamped to 0. It fails the adoption rule's first condition.
- (b) and (d) are dominated, (d)'s ceiling shows the shared response isn't the bottleneck, κ is noise, and w_arc > 1 hurts.

**Expected change against C+G, noise-free scorer-space estimates.** Only MSE and PDS can be estimated without cells; the four DE components need the end-to-end run.

| Fold line | MSE scaled (clamped): C+G → config | PDS scaled: C+G → config | Δ overall from these two (÷6) |
|---|---|---|---|
| H1 / Eval A | 0 → **+0.118** (raw 1.496 → 0.853) | 0.723 → 0.694 (−0.029) | **+0.015** |
| HCT116 (fold pool) | 0 → +0.332 | 0.528 → 0.282 (−0.246) | +0.014 |
| HEK293T (fold pool) | 0 → 0 (unclamped −3.48 → −0.33) | 0.149 → 0.061 (−0.088) | −0.015 |
| K562 | +0.033 → +0.058 | 0.511 → 0.400 (−0.111) | −0.014 |

**Against the adoption rule** (Eval A: MSE scaled > 0, every other component > 0, overall above C+G with paired CI lower bound > −0.01; HCT116 overall within its CI):
- **MSE > 0 is expected** (+0.118).
- **PDS stays positive** (0.694).
- **The other four** (nmae, fidelity, reach, Jaccard) are unknown until cells exist. They are positive for C+G today, and the call set is unchanged.
- **The overall gain expected from MSE + PDS is about +0.015**, on the order of the paired CI width. **It may not pass.**
- **The PDS cost is the main risk:** shrinking every effect pulls every knockdown's profile toward control and costs discrimination. It's −0.03 on H1, but −0.11 to −0.25 on the weak lines.

**After approval:** run this config against the incumbent end to end on Eval A and HCT116final (generated through `run_line --lfc-npz` with C+G's calls and 0.5 × its LFC), with the paired bootstrap and all six components. That's one job per VM, on whichever VM is idle.

**What I would have done differently:**
1. **Used the noise-free pseudobulk proxy from the start.** It reproduces end-to-end MSE to 0.0006 (H1) and costs seconds per arm instead of hours. Most of this week's end-to-end runs could have been screened this way first.
2. **Checked the metadata first.** Chemistry isn't recorded anywhere, so Step 7's main part was blocked. Marson's missing full profiles opened a 24–26% Δ hole in the HCT116+HEK293T fold.
4. **One caveat on the evidence:** every number here is noise-free with Ĉ_pred = 1.13 × C0. The budget B binds for arm (a) at m ≤ 0.5, so the realised credit in generated cells could differ. Only the end-to-end run settles it.

**Stopped here, waiting for Jose's approval.**

---

## Part 2 — Steps 9–15 (task: Jose, 2026-10-28 ~17:30 UTC)
**Goal:** Eval A MSE scaled > 0, every other component > 0, overall above C+G (0.219). HCT116 and HEK293T are guardrails only.

**Infrastructure:**
- **K562 and HEK293T needed DE references for the four DE proxies.** Their Steps 3–8 references were MSE/PDS-only.
- They are being rebuilt with the harness's own `finish` (real cache, ceiling, generic baseline, 5 anchor splits) on the existing cells: `final_evals.py build <X>final --finish-only`.
- The MSE/PDS-only cache, baseline and bundle moved to `eval<X>final/msepds_2026-10-28/`. Nothing was deleted.
- I stopped the first attempt (17:41–17:49) after 8 min, so no table load ran alongside it. It restarts after the Step 10 fold builds (`results/stage1/2026-10-28_scorer_space/logs/step10_folds.sh`).

### Step 9 — proxies for all six components (`scorer_proxy.py step9`; `results/stage1/2026-10-28_scorer_space/tables/step9*.{json,csv,parquet}`)
**Proxy:**
- **Means:** the requested LFC → the generator's own mean inversion, μ = (m0 + 1e-9)·2^lfc − 1e-9. m0 is the run's context-control CPM.
- **MSE and PDS:** P = log1p(5e4·μ/Σμ), as in Steps 3–5, with Ĉ_pred = C0 × 1.1295.
- **DE components:** a synthetic pred DE table on the real DE table's (target, feature) rows, passed through cell-eval2's own `prep_de_side → rank_de_side → assemble_prepared_de → run.dispatch_de_metrics` (the path `cell-eval2 run` takes), against the eval's cached real DE.
  - **LFC:** pdex's log2((mean_pred + 1e-9)/(mean_ctrl + 1e-9)). mean_pred is μ renormalised to 1e6, as the scorer's per-cell normalisation does. mean_ctrl is the real scoring controls' arithmetic mean CPM (`control_source: real`).
  - **Significance:** p = 1e-3/(1 + 1e4·s) for requested calls and 0.5 + 0.5/(1 + 1e4·s) otherwise, with s = |P − real control bulk|. So the requested set is exactly the significant set, ranked by effect size in scorer space.
- **Scaling:** cell-eval2's `score_one` on each bundle's replicate entries. Overall = mean of the six.

**Plumbing checks (all pass):**
- The same aggregation and scaling, applied to each run's own per-knockdown `results.csv`, reproduce its `summary.json` (every component and overall within 7e-16).
- cell-eval2's DE dispatch, re-run on each run's own `de_pred.parquet`, reproduces `results.csv` per knockdown (max |diff| 1.2e-15).

**Proxy vs actual (scaled; raw in brackets):**

| Component | Eval A C+G: proxy / actual | Eval A Model 2 port | HCT116 C+G (screen setup) |
|---|---|---|---|
| MSE | 0 / 0 (1.4956 / 1.4962) | 0 / 0 (2.128 / 2.129) | 0 / 0 (2.987 / 2.896) |
| PDS | 0.7235 / 0.7234 | 0.6627 / 0.6629 | 0.6542 / 0.6549 |
| NMAE | 0.0402 / 0.0404 | −0.0369 / −0.0368 | 0.4004 / 0.4013 |
| Direction fidelity | 0.3711 / 0.3802 | 0.3429 / 0.3571 | 0.9155 / 0.9146 |
| Direction reach | 0.0739 / 0.0701 | 0.0485 / 0.0468 | 0.7443 / 0.7414 |
| Jaccard | 0.0953 / 0.0975 | 0.0689 / 0.0754 | −1.5714 / −1.5705 |
| **Overall** | **0.2173 / 0.2186 (−0.0013)** | **0.1810 / 0.1842 (−0.0032)** | **0.1905 / 0.1903 (+0.0002)** |

- **The expected overall tracks the actual one to within 0.003 on all three runs, and it ranks C+G above the Model 2 port by the same margin** (0.036 vs 0.034 actual).
- The weakest proxy is direction fidelity, 0.009–0.014 low on Eval A. The realised calls add true-direction calls the proxy doesn't request (below).
- On HCT116 the raw MSE proxy is 0.09 high (2.987 vs 2.896). It stays clamped at 0 either way: the denominator is tiny.

**Realised minus requested calls** (realised = harness `de_pred`, p_adj < 0.05, own gene out). By the target's context-control CPM decile (gated genes, 0 = lowest; −1 = ungated):

| Decile | Eval A C+G: realised − requested per knockdown (not requested / requested not realised, total) | Eval A M2 | HCT116 C+G |
|---|---|---|---|
| ungated | +0.2 (63 / 0) | +0.8 (223 / 0) | +0.1 (19 / 0) |
| 0 | +0.6 (202 / 13) | +5.1 (1,596 / 103) | +0.1 (33 / 6) |
| 1 | +1.4 (402 / 0) | +5.0 | +0.2 |
| 2 | +0.6 | +2.0 | −0.1 |
| 3 | −0.1 | +2.2 | −0.2 |
| 4 | +0.5 | +4.6 | −0.2 |
| 5 | +2.7 (804 / 22) | +7.0 | −0.2 |
| 6 | +1.8 | +8.4 | −0.4 |
| 7 | +6.1 (1,801 / 32) | +12.0 | −0.6 |
| 8 | +8.4 (2,505 / 39) | +15.0 | −0.8 |
| 9 | **+12.2 (3,930 / 363)** | **+25.4 (8,233 / 827)** | −2.9 (5 / 807) |
| **Per knockdown, median** | **+28** (requested 486.5 → realised 520) | **+20** (249.5 → 275.5) | **−3** (132.5 → 125.5) |

- **On H1 the generator over-calls**, mostly in the top expression deciles, where small kept shifts turn significant against the real controls (as Step 3b found).
- **On HCT116 it under-calls slightly**, again in the top decile: about 3 requested calls per knockdown aren't realised.
- **Requested-but-not-realised is rare on H1** (median 0 per knockdown).
- Every arm from here is ranked by this expected overall.

### Step 10 — source filters × magnitude arms (`scorer_part2.py step10`, `step10fit`; `results/stage1/2026-10-28_scorer_space/tables/step10*.{json,csv,parquet}`)

**Carried: F0 (all tables), arm (a) at m = 0.5.** On Eval A it is the only qualifying family: MSE scaled 0.118 and every component > 0. Expected overall is **0.2283, against 0.2173 for C+G (+0.011)**. The full fit matches the provisional carry.

**Eval A (H1 held out), scaled:**

| Filter | Arm | Overall | MSE | PDS | NMAE | Fidelity | Reach | Jaccard | Qualifies |
|---|---|---|---|---|---|---|---|---|---|
| F0 | C+G | 0.217 | 0 | 0.723 | 0.040 | 0.371 | 0.074 | 0.095 | no (MSE 0) |
| F0 | **a 0.5** | **0.228** | **0.118** | 0.694 | 0.020 | 0.371 | 0.071 | 0.095 | **yes** |
| F0 | c 0.1 | 0.202 | 0 | 0.761 | −0.051 | 0.371 | 0.038 | 0.095 | no |
| F0 | e (LOFO: m_c 1, m 0.2) | 0.210 | 0 | 0.772 | −0.038 | 0.371 | 0.060 | 0.095 | no |
| F1 | C+G | 0.196 | 0 | 0.664 | 0.022 | 0.355 | 0.041 | 0.096 | no |
| F1 | a 0.5 | 0.209 | 0.107 | 0.648 | 0.012 | 0.355 | 0.035 | 0.096 | yes |
| F1 | c 0.1 | 0.190 | 0 | 0.720 | −0.057 | 0.355 | 0.026 | 0.096 | no |
| F1 | e | 0.195 | 0 | 0.726 | −0.045 | 0.355 | 0.042 | 0.096 | no |
| F2 | C+G | 0.165 | 0 | 0.568 | 0.001 | 0.301 | 0.025 | 0.091 | no |
| F2 | a 0.5 | 0.179 | 0.100 | 0.554 | 0.007 | 0.301 | 0.020 | 0.091 | yes |
| F2 | c 0.1 | 0.163 | 0 | 0.639 | −0.066 | 0.301 | 0.015 | 0.091 | no |
| F2 | e | 0.168 | 0 | 0.645 | −0.055 | 0.301 | 0.027 | 0.091 | no |

**Every line, overall (scaled), held out:**

| Filter / arm | K562 | HCT116 | HEK293T | H1 (Eval A) |
|---|---|---|---|---|
| F0 C+G | 0.146 | 0.124 | 0.106 | 0.217 |
| F0 a 0.5 | 0.115 | 0.117 | 0.072 | **0.228** |
| F0 c 0.1 | 0.158 | 0.319 | 0.139 | 0.202 |
| F0 e (LOFO) | 0.165 | 0.312 | 0.161 | 0.210 |
| F1 C+G / a / c / e | 0.168 / 0.114 / 0.158 / 0.167 | 0.102 / 0.101 / 0.320 / 0.315 | 0.088 / 0.070 / 0.139 / 0.164 | 0.196 / 0.209 / 0.190 / 0.195 |
| F2 C+G / a / c / e | 0.085 / 0.066 / 0.096 / 0.095 | 0.136 / 0.138 / 0.332 / 0.300 | 0.103 / 0.090 / 0.172 / 0.173 | 0.165 / 0.179 / 0.163 / 0.168 |

The full six-component rows for every line are in `results/stage1/2026-10-28_scorer_space/tables/step10_heldout.csv`.

**What the table says:**
- **Arm (a) wins only on H1.** It scales every Δ by 0.5, which buys MSE (0 → 0.118 on H1) but costs PDS on every line: H1 −0.03, K562 −0.11, HCT116 −0.25, HEK293T −0.09. On H1 the MSE gain outweighs that loss; on the other three lines it doesn't.
- **Arms (c) and (e) are the reverse.** They are much better on the guardrail lines (HCT116 0.32 vs 0.12), because PDS rises and HCT116's MSE clears. On H1 they lose NMAE (−0.04 to −0.05): the source Δ on non-called genes is wrong-sized often enough to hurt the LFC error. MSE stays clamped at 0.
- **Arm (e) never qualifies on H1 within the specified grid.** NMAE is negative in all 36 grid cells. LOFO puts non-called m at 0.2, the top of the {0, 0.1, 0.2} grid, on every fold and filter, and m_c at 0.75 (K562 and HCT116+HEK293T folds) or 1 (H1 fold). NMAE improves as m rises (F0, m_c 0.5: −0.041, −0.027, −0.014), so the grid edge matters. See the side check below.
- **Fidelity, reach and Jaccard barely move between arms.** All arms within a filter use the same call set (C+G's), so these three components change only through filter-dependent calls. Step 12 is where they are tested.
- **J²-weighted r in Δ space:**
  - H1: 0.156 (C+G), 0.155 (a, a pure rescale), 0.128 (c), 0.136 (e).
  - The other lines are 0.02–0.04 for every arm, and c/e go negative on HEK293T (−0.02).
  - The extra source Δ that c/e add on non-called genes lowers r on H1, which is why NMAE falls.
- **Filters hurt, everywhere it matters.** F1 costs Eval A 0.02 and F2 costs 0.05, under every arm. The cleaner sources don't make up for the coverage lost:

| Filter | Tables | H1: share of (k, g) with no source | H1 knockdowns with no source | D/E/F targets with no source |
|---|---|---|---|---|
| F0 | 117 | 0.009 | 1 | 1 (ANG) |
| F1 | 27 | 0.227 | 21 | 2 (ANG, RAB5A) |
| F2 | 27 (rows: efficiency ≤ log2 0.2, known) | 0.354 | 64 | 90 |

  Under F2, K562 loses 61% of (k, g) pairs (159 knockdowns with no source), and C+G's median K562 call count drops to 0.
- C+G's own calls change with the filter, because its posterior pools the kept sources. H1's median requested calls per knockdown are 486.5 under F0, 589 under F1 and 503.5 under F2.

**Side check: arm (e) past the grid edge** (`scorer_part2.py step10ext`, `step10extfit`; `results/stage1/2026-10-28_scorer_space/tables/step10ext*.{json,parquet}`). Non-called m ∈ {0.3, 0.5, 0.75, 1}, F0, all m_c. This is outside the task's grid, so it does not enter the carry rule.

H1, expected overall (NMAE in brackets):

| m_c \ m | 0.2 | 0.3 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| 0.25 | 0.223 | 0.227 | 0.234 (+0.004) | 0.215 (MSE 0) | 0.208 |
| 0.5 | 0.231 (−0.014) | 0.238 (−0.001) | **0.239 (+0.020, MSE 0.095, all > 0)** | 0.222 (MSE 0) | 0.215 |
| 0.75 | 0.234 | 0.234 | 0.223 (MSE 0) | 0.224 | 0.217 |
| 1 | 0.210 | 0.215 | 0.219 | 0.221 | 0.216 |

- **In sample, H1 peaks at m_c 0.5 and m 0.5 (0.239).** That beats arm (a)'s 0.228, and every component is positive.
- **Held out, the peak is not reached.** The K562 and HCT116+HEK293T folds prefer m_c 0.75 (or 1) and m 0.3. LOFO on the extended grid gives H1 m_c 1, m 0.3: 0.215, with MSE 0 and NMAE −0.025, so it doesn't qualify.
- **Arm (a)'s m = 0.5 isn't held out either.** It was carried from Part 1 and was not refit here.
- **Step 11 needs no further extension of (e).** Its arms (f–i) fix the called genes at the carried 0.5 × LFC and fit only the non-called side, which is the m_c 0.5 row above.

### Step 11 — adding signal on the carried config (`scorer_s11.py`, `scorer_step11.py`; `results/stage1/2026-10-28_scorer_space/tables/step11*.{json,parquet}`)
Run on the provisional carry that Step 10's fit then confirmed.

The arms keep the called genes at the carried 0.5 × LFC and change only the non-called side, with parameters fit by LOFO.

**Carried: arm (i), confidence-weighted source Δ.** Non-called Δ_pred = 0.5 × c × Δ_pool, with c = z_a²/(z_a² + 1) and c₁ = 1 for single-source genes.
- **Eval A: 0.2436**, against 0.2283 for the Step 10 carry and 0.2173 for C+G.
- Components: MSE 0.113, PDS 0.808, NMAE **0.0010**, fidelity 0.371, reach 0.073, Jaccard 0.095. r_J2_delta is 0.156.

**Held out, overall (scaled):**

| Arm (H1-fold parameters) | K562 | HCT116 | HEK293T | H1 (Eval A) | Eval A MSE / NMAE | Qualifies |
|---|---|---|---|---|---|---|
| Step 10 carry (a 0.5) | 0.115 | 0.117 | 0.072 | 0.228 | 0.118 / 0.020 | yes |
| f: strength terciles (0.2, 0.2, 0.2) | 0.153 | 0.275 | 0.142 | 0.231 | 0.084 / −0.014 | no |
| g: source PCs (K 50, m 0.3) | 0.101 | 0.235 | 0.122 | 0.208 | 0.075 / −0.025 | no |
| h: partners (n 0 → no partners, m 0.2) | 0.153 | 0.275 | 0.142 | 0.231 | 0.084 / −0.014 | no |
| **i: confidence (m 0.5, κ 1, c₁ 1)** | **0.156** | **0.292** | **0.149** | **0.244** | **0.113 / 0.001** | **yes** |
| C+G (Step 10, for reference) | 0.146 | 0.124 | 0.106 | 0.217 | 0 / 0.040 | no |

- **Arm (i) also beats C+G and the Step 10 carry held out on the other three lines.** It is the first config in this study that does.
- **Fragility, which Step 14 has to weigh:**
  - **NMAE clears zero by only 0.001.**
  - **All three parameters sit on a grid edge, and LOFO picks the same values on every fold:**
    - m = 0.5 is the top of {0.1, 0.2, 0.3, 0.5}. The Step 10 side check shows m 0.75 drives MSE to 0 at m_c 0.5, so going higher is unlikely to help.
    - κ = 1 is the bottom of {1, 4, 16}. No shrinkage (κ → 0) is arm (e) at (0.5, 0.5): 0.239. So the agreement shrinkage is worth +0.005, and its optimum lies between κ 0 and 4.
    - c₁ = 1 is the top of {0, 0.5, 1}, meaning single-source genes are not down-weighted.
  - **Fold spreads:** m is 0.5 on all folds and κ is 1 on all folds. c₁ is 1, 0.5 and 1 (max/min 2, the limit).

H1 grid for arm (i), overall (NMAE):

| κ, c₁ \\ m | 0.1 | 0.2 | 0.3 | 0.5 |
|---|---|---|---|---|
| 1, 1 | 0.217 | 0.228 | 0.234 (−0.014) | **0.244 (+0.001)** |
| 1, 0 | 0.216 | 0.226 | 0.232 | 0.241 (−0.001) |
| 4, 1 | 0.215 | 0.221 | 0.228 | 0.237 (−0.014) |
| 16, 1 | 0.209 | 0.216 | 0.219 | 0.225 (−0.027) |

#### S1 — strength
Strength = J²-weighted ||Δ_k||² over gated genes, own target out.

Spearman between lines on shared knockdowns:

| Pair | Shared | ρ |
|---|---|---|
| HCT116 ~ HEK293T | 260 | **0.83** |
| K562 ~ HCT116 | 249 | 0.11 |
| K562 ~ HEK293T | 242 | 0.05 |
| K562 ~ H1 | 27 | 0.34 |
| HCT116 ~ H1 | 27 | 0.11 |
| HEK293T ~ H1 | 25 | −0.04 |

- **Only the two Orion screens agree on strength.**
- **Pooled predicted vs true strength:** K562 0.02, HCT116 0.07, HEK293T 0.03, H1 0.22.
- **Arm (f) collapses to a uniform m.** LOFO picks 0.2 for every tercile on every fold, so it is arm (e) at m_c 0.5, m 0.2. Strength prediction adds nothing.

#### S2 — PCs of the source lines' Δ (fold-excluded)
Energy share of the line's true Δ captured by K PCs:

| Line | Source PCs, K 5 / 10 / 20 / 50 | Own top K (ceiling), K 50 | Random (floor), K 50 |
|---|---|---|---|
| K562 | 0.016 / 0.027 / 0.036 / 0.047 | 0.43 | 0.006 |
| HCT116 | 0.035 / 0.050 / 0.066 / 0.080 | 0.44 | 0.005 |
| HEK293T | 0.065 / 0.091 / 0.115 / 0.135 | 0.52 | 0.005 |
| H1 | 0.043 / 0.087 / 0.125 / 0.213 | 0.75 | 0.005 |

- **The shared subspace is real, but it holds only a fifth of H1's energy at K 50.** That is 40× random but well under H1's own ceiling.
- **Arm (g):** LOFO picks K 50 (the top of the grid) on every fold, with m 0.5 / 0.2 / 0.3. Projecting first lowers H1 to 0.208 and drops reach to 0.011. Projection throws away the per-gene signal that reach and NMAE need.

#### S3 — complex partners (sources only)
- **D/E/F targets:** 0/300 have a partner at cosine > 0.5. 267 of them are in K562 GW, where the median best cosine is 0.115.
- **Arm (h):** LOFO picks n = 0 on every fold. Every partner setting lowers H1 (n 3: 0.201–0.208; n 10: 0.194–0.202, against 0.231 at n 0).
- **The 30 strongest H1 knockdowns** are almost all chromatin and transcription machinery, plus H1's own pluripotency factors:

| Group | Knockdowns in H1's top 30 |
|---|---|
| m6A methyltransferase complex | METTL3, METTL14 (each the other's K562 partner) |
| BAF/PBAF | SMARCA4, SMARCE1 |
| TFIID | TAF11, TAF13 |
| EMSY complex | SIN3B, PHF12 |
| bHLH | MAX, USF2 |
| Pluripotency TFs | PRDM14, SOX2, SALL4 |
| Others | MED19, KDM1A, EP300, DNMT1, DOT1L, INO80B, SUPT4H1, RNF40, EIF3H |
| Sterol synthesis (my reading; no HGNC group) | FDPS, NSDHL |

  HGNC groups shared within the top 30 are listed in `results/stage1/2026-10-28_scorer_space/tables/step11.json` → `h1_top30_shared_groups`. Only METTL3/METTL14 are both in the top 30 and each other's partners.

#### S4 — per-(k, g) confidence
β and R0² are the slope and uncentred R² of true on pooled Δ, in Δ space.

β by decile, low → high (Spearman of β over the ten deciles in brackets):

| Line | Agreement z (iii) | C+G z = \|LFC\|/SE | Source tier 1 / 2–3 / 4+ |
|---|---|---|---|
| K562 | 0.016 → 0.051, peak 0.063 at d7 (0.89) | 0.023 → 0.057 (0.90) | 0.023 / 0.042 / 0.125 |
| HCT116 | 0.017 → 0.052, peak 0.058 at d8 (0.96) | 0.013 → 0.024, flat from d3 (0.50) | 0.013 / 0.036 / 0.081 |
| HEK293T | 0.010 → 0.038 (0.85) | 0.000 → 0.017, flat (0.60) | 0.011 / 0.030 / **0.005** |
| H1 | 0.051 → 0.272, strictly rising (1.00) | 0.043 → 0.169 (1.00) | 0.063 / 0.055 / 0.232 |

- **Agreement z is the only measure whose β rises on every line.** It rises 3–5× from bottom to top decile. It is strictly monotone only on H1; elsewhere it has one or two dips.
- C+G's z is flat on the two Orion lines. Source tier breaks on HEK293T.
- R0² stays small even in the top agreement decile (H1 0.050, the others ≤ 0.003).
- Median sign agreement among sources is 0.67 (K562, H1) and 1.0 (HCT116, HEK293T).
- This is why arm (i) works. Agreement z ranks which non-called source Δ to trust. The fit prefers light shrinkage (κ 1) at a larger m (0.5) over heavy shrinkage.

### Step 12 — call rule on the carried magnitudes (`scorer_step12.py step12`, `step12fit`; `results/stage1/2026-10-28_scorer_space/tables/step12*.{json,csv,parquet}`)
Every variant uses arm (i)'s magnitudes, so MSE, PDS and NMAE are the same in each row. K1 at m_call 1 reproduces K0's call set exactly on all four lines (0 mismatches).

**Eval A call precision per decile, first.** These are the realised C+G calls from Step 9: 171,056 calls, 34,016 true, precision 0.20.

| Decile | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|---|
| Precision | **0.26** | **0.21** | **0.20** | **0.18** | 0.16 | 0.15 | 0.15 | 0.15 | 0.18 | **0.26** |

The median is 0.182, so K3 keeps deciles 0–3 and 9.

**Carried: K0 (C+G's calls).** On Eval A it stays at 0.2436, unchanged from Step 11. No held-out variant beats it.

**The fold design can't fit a call rule.** The truth call set per knockdown has a median of 218 on H1 but 3 on K562 and 1 on HCT116 and HEK293T. So when H1 is held out, both training folds reward calling almost nothing:
- **K1:** LOFO picks m_call 0.15 for the H1 fold. That gives 0 calls and H1 0.154.
- **K4:** LOFO picks N_med 50, giving H1 0.202.
- **On the other lines** the same truth makes K2_4 (≥ 4 same-sign sources, almost no calls) the best variant: K562 0.204–0.220, HCT116 0.439, HEK293T 0.220. Those gains come from Jaccard moving from about −1.2 to about +0.08 by calling nothing. That is the near-noise denominator the task warns about.

**Held out, overall (variant at the H1 fold's parameters):**

| Variant | K562 | HCT116 | HEK293T | H1 (Eval A) | H1 median calls (truth 218) | H1 fidelity / reach / Jaccard |
|---|---|---|---|---|---|---|
| **K0** | 0.156 | 0.292 | 0.149 | **0.244** | 486.5 | 0.371 / 0.073 / 0.095 |
| K1 (m_call 0.15 on H1 and K562; 1 on HCT116 and HEK293T) | 0.183 | 0.292 | 0.149 | 0.154 | 0 | −0.100 / 0.081 / 0.022 |
| K2, ≥ 2 sources | 0.197 | 0.455 | 0.220 | 0.154 | 0 | −0.102 / 0.081 / 0.022 |
| K2, ≥ 4 sources | 0.204 | 0.439 | 0.220 | 0.154 | 0 | −0.103 / 0.081 / 0.022 |
| K3 (deciles 0–3, 9) | 0.183 | 0.290 | 0.148 | 0.154 | 0 | −0.100 / 0.081 / 0.022 |
| K4 (N_med 800 on K562 and HCT116+HEK293T; 50 on H1) | 0.144 | 0.302 | 0.154 | 0.202 | 50 | 0.191 / 0.078 / 0.019 |

**H1 in sample (not held out, not eligible for the carry):**

| Variant | Median calls | Overall | Fidelity | Reach | Jaccard |
|---|---|---|---|---|---|
| K0 = K1 at m_call 1 | 486.5 | 0.244 | 0.371 | 0.073 | 0.095 |
| K2 ≥ 2 / ≥ 4 at m_call 1 | 380.5 / 39.5 | 0.238 / 0.191 | 0.347 / 0.077 | 0.077 / 0.086 | 0.085 / 0.060 |
| K3 at m_call 1 | 198.5 | 0.223 | 0.290 | 0.064 | 0.062 |
| K1 at m_call 0.5 / 0.25 / 0.15 | 16 / 1 / 0 | 0.188 / 0.158 / 0.154 | 0.110 / −0.067 / −0.100 | 0.081 | 0.015 / 0.009 / 0.022 |
| K4, N_med 100 / 200 / 400 / **800** | 100 / 200 / 400 / 800 | 0.214 / 0.228 / 0.242 / **0.260** | 0.252 / 0.314 / 0.379 / **0.449** | 0.077 / 0.075 / 0.071 / 0.075 | 0.034 / 0.055 / 0.082 / **0.114** |

- **On H1, more calls score better, all the way to the top of the grid.** K4 at 800 calls per knockdown (about 3.7× truth) gives 0.260, against 0.244 for K0. Fidelity and Jaccard both keep rising.
- **This matches Step 9:** the generator's +28 leaked calls per knockdown raised C+G's actual Eval A score above the proxy (fidelity +0.009, Jaccard +0.002, overall +0.001).
- **Every restriction lowers H1:** source support (K2), dropping the low-precision deciles (K3), and scaling the test (K1).
- **K4's N_med = 800 is a grid edge, and an H1-only fit.** No other line can confirm it, so it does not enter the carry rule.

### Step 13 — generator check (`scorer_step12.py step13`; `results/stage1/2026-10-28_scorer_space/tables/step13.json`) — measurement only, run ahead of Steps 10–12 because it needs no carried config
**Cells:** the existing Eval A C+G (10-28) and HCT116 C+G (screen setup) cells, against each eval's real cells. Realised calls are taken from the harness `de_pred` (p_adj < 0.05, own gene out). Truth is the eval's real DE table at p_adj < 0.05.

**1. Per-cell library size (5 / 25 / 50 / 75 / 95%):**

| | Generated knockdown cells | Real knockdown cells | Real scoring controls |
|---|---|---|---|
| Eval A | 19,639 / 19,871 / **20,026** / 20,185 / 20,473 | 10,655 / 15,576 / **19,662** / 24,698 / 34,332 | 11,015 / 15,951 / 19,989 / 25,050 / 34,595 |
| HCT116 | 11,699 / 11,945 / **12,125** / 12,322 / 12,668 | 5,668 / 8,892 / **12,188** / 16,878 / 27,162 | 5,733 / 8,861 / 12,145 / 16,732 / 26,782 |

- **The medians match** (generated/real knockdown 1.019 on Eval A, 0.995 on HCT116).
- **The spread does not.** The generator builds every knockdown cell at the controls' median library × Σ requested CPM / 1e6, so its depth CV is about 1%, against about 35–45% in real cells.

**2. Per-gene zero fraction, generated vs real cells of the same knockdown, by the target's control-CPM decile:**

| Decile | Eval A: mean gen − real (share \|diff\| > 0.05) | HCT116 |
|---|---|---|
| 0 (lowest) | +0.010 (7%) | +0.016 (12%) |
| 2 | +0.011 (18%) | +0.028 (29%) |
| 4 | +0.006 (19%) | +0.034 (38%) |
| 6 | −0.004 (13%) | +0.030 (37%) |
| 8 | −0.005 (4%) | +0.014 (23%) |
| 9 (highest) | −0.002 (1%) | +0.013 (9%) |

On HCT116 generated cells have slightly more zeros in the middle deciles. On Eval A the zero fractions are close everywhere, and closest in the top deciles.

**3. Realised calls that were not requested ("leaked"):**

| | Eval A C+G | HCT116 C+G |
|---|---|---|
| Leaked calls (per knockdown, median) | **10,650 (28)** | 147 (0) |
| … of which also false against the truth | 8,928 (84%) | 147 (100%) |
| … in control-CPM deciles 7–9 | **8,236 (77%)** | 5 |
| … in deciles 0–2 | 781 | 122 |
| \|requested scorer-space Δ\| of leaked calls, 10 / 50 / 90% | 0.082 / 0.128 / 0.171 | 0.116 / 0.153 / 0.181 |
| … with \|zero-fraction gen − real\| > 0.05 | **1,160 (11%)** | 78 (53%) |
| Requested calls not realised | 596 (median 0) | 1,532 (median 3) |
| All realised calls false against the truth | 137,040 of 171,056 | 52,355 of 53,219 |

**Reading:**
- **On Eval A the leaked calls are not a zero-fraction (dropout) effect.** 89% of them sit on genes whose zero fractions match the real cells to within 0.05, and 77% are in the top three expression deciles, where zeros barely occur.
  - They carry real requested shifts (median |Δ| 0.13 in scorer space) that the call rule left uncalled.
  - The generator's near-constant depth removes the between-cell depth variance real cells have. At 400 cells per knockdown, the Mann–Whitney test against the real controls then has more power than it would on real data, so those kept shifts turn significant.
- **On HCT116 the few leaked calls (147) are low-expression genes.** Half of them have a zero-fraction mismatch, where the generator over-produces zeros in the low and middle deciles.

**How many could a depth-matched generator plausibly remove?**
- **Upper bound on Eval A:** the 8,236 leaked calls in deciles 7–9 (77%, median about 22 per knockdown), plus at most the 1,160 with zero mismatches. That's about 9,400 of 10,650, of which about 7,900 are also false against the truth.
- **HCT116:** at most the 78 with zero mismatches.
- This is a plausibility bound, not a measured effect. Real-like depth spread would also cut power on the requested calls, so some requested calls (596 not realised today on Eval A) could drop out too. Only a generator run would measure the net. **The generator is not changed here.**
- **For scale:** the false realised calls on Eval A (137,040) are dominated by the requested set itself (about 128,000 requested calls that the truth does not hold), not by leakage.

### Step 14 — recommendation (stop for Jose's approval)

**Two configs, by held-out expected Eval A overall.** Both use all source tables (F0) and C+G's calls (K0):

| | Config 1: F0 / arm (i) / K0 | Config 2: F0 / arm (a) 0.5 / K0 | Incumbent C+G |
|---|---|---|---|
| Called genes | 0.5 × C+G LFC | 0.5 × C+G LFC | C+G LFC |
| Non-called genes | 0.5 × c × Δ_pool, c = z_a²/(z_a² + 1); single-source c₁ = 1 | 0.5 × C+G LFC | C+G LFC |
| **Expected Eval A overall** | **0.2436 (+0.026)** | **0.2283 (+0.011)** | 0.2173 |

Six proxy components (scaled) on every line, held out:

| Line | Config | Overall | MSE | PDS | NMAE | Fidelity | Reach | Jaccard |
|---|---|---|---|---|---|---|---|---|
| **H1 (Eval A)** | 1 | **0.244** | 0.113 | 0.808 | **0.001** | 0.371 | 0.073 | 0.095 |
| | 2 | 0.228 | 0.118 | 0.694 | 0.020 | 0.371 | 0.071 | 0.095 |
| | C+G | 0.217 | 0 | 0.723 | 0.040 | 0.371 | 0.074 | 0.095 |
| K562 | 1 | 0.156 | 0.140 | 0.622 | 0.101 | 0.843 | 0.346 | −1.117 |
| | 2 | 0.115 | 0.058 | 0.400 | 0.118 | 0.843 | 0.386 | −1.117 |
| | C+G | 0.146 | 0.033 | 0.511 | 0.221 | 0.843 | 0.383 | −1.117 |
| HCT116 (guardrail) | 1 | 0.292 | 0.848 | 0.780 | 0.207 | 0.871 | 0.645 | −1.598 |
| | 2 | 0.117 | 0.332 | 0.282 | 0.216 | 0.871 | 0.598 | −1.598 |
| | C+G | 0.124 | 0 | 0.528 | 0.343 | 0.871 | 0.599 | −1.598 |
| HEK293T (guardrail) | 1 | 0.149 | 0.123 | 0.365 | 0.151 | 0.925 | 0.585 | −1.252 |
| | 2 | 0.072 | 0 | 0.061 | 0.159 | 0.925 | 0.537 | −1.252 |
| | C+G | 0.106 | 0 | 0.149 | 0.225 | 0.925 | 0.591 | −1.252 |

The negative Jaccard on the three non-H1 lines comes from their near-empty truth (median 1–3 calls per knockdown against C+G's 113–592). It is the same for every config, because the calls are C+G's.

**Parameter spreads across folds:**

| Parameter | Folds (K562 / HCT116+HEK293T / H1) | Max/min | Qualifies (m-type ≤ 2) | Note |
|---|---|---|---|---|
| Config 1, called scale | 0.5 fixed | — | — | Carried from Part 1, not refit |
| Config 1, m | 0.5 / 0.5 / 0.5 | 1 | yes | Top of the grid |
| Config 1, κ | 1 / 1 / 1 | 1 | yes | Bottom of the grid; κ → 0 is worse (arm (e) 0.239) |
| Config 1, c₁ | 1 / 0.5 / 1 | 2 | yes, at the limit | |
| Config 2, m | 0.5 fixed | — | — | Part 1's choice, not refit in Part 2 |

**Why these two:**
- Config 1 is the best held-out config that qualifies, and it is above C+G on all four lines.
- Config 2 is the next best qualifying config. It works by a different mechanism (a pure rescale), and its NMAE margin is 0.020 against config 1's 0.001.
- **The main risk to config 1 under Step 15's rule is NMAE.** Step 9's proxy matched actual NMAE to within 0.0002 on Eval A (0.0402 vs 0.0404; −0.0369 vs −0.0368), so +0.001 should survive. But it is thin, and the bootstrap may straddle zero.

**Considered and not recommended:** both of these are H1-fitted, in-sample only.
- **Arm (e) at m_c 0.5, m 0.5:** 0.239. LOFO doesn't reach it, and config 1 beats it held out.
- **Config 1 with K4 at 800 calls per knockdown:** 0.260. N = 800 sits at the grid edge, and no other line can confirm it (Step 12). Scoring it in Step 15 would grade it on the same H1 data it was picked on. If Jose wants the extra calls tested, the honest version is a two-fold split of H1's knockdowns to fit N, run on vcc2027 with the proxy only (about 10 minutes).

**Does Step 13 support a generator change? No.**
- **The leaked calls help Eval A on net.** Step 9's actual minus proxy is +0.0013 overall for C+G (fidelity +0.009, Jaccard +0.002, reach −0.004) and +0.003 for Model 2.
- **More calls help too.** Step 12 shows H1's fidelity and Jaccard rising with call count all the way to 800 per knockdown.
- **So matching depth would likely cost score.** A depth-matched generator would remove up to about 9,400 leaked calls, mostly false. On the evidence here that would lower fidelity and Jaccard rather than raise them.
- **HCT116's leakage is negligible** (147 calls).
- The net is unmeasured without a generator run, and the generator stays frozen.

**What I would have done differently:**
1. **Wider grids from the start.** Most fitted parameters landed on a grid edge:
   - arm (e)'s m;
   - all three of arm (i)'s parameters;
   - arm (g)'s K;
   - K4's N.

   A proxy row costs seconds, so wider grids would have been cheap. The Step 10 side check fixed this for (e) only.
2. **Check each fold's call truth before Step 12.** Only H1 has a usable truth call set (median 218, against 1–3 elsewhere). Fold-level LOFO was bound to pick "call almost nothing" for H1. A knockdown-level split within H1 would have given a fair fit of the call rule.
3. **Refit arm (a)'s m inside the same LOFO** as (e) and (i), rather than carrying 0.5 from Part 1. Config 2's 0.228 is not held out in the same sense as config 1's.

**Step 15 plan, after approval:**
- Config 1 and config 2, in parallel.
- Each runs on Eval A and HCT116final (screen setup, HEK293T voting) against the incumbent C+G scores.
- Report: paired bootstrap of all six components, scaled and raw; realised vs requested calls; the adoption rule.
- Stop and report before building any D/E/F file.

## Step 16 — wider magnitude grid and H1 knockdown-split call budgets (task: Jose, 2026-10-29)
**Order of work.** Jose's note put Step 16 after Step 15, but Step 15 had not been approved or started. Asked at 07:30 UTC, Jose chose Step 16 first.

Noise-free proxy only; no cells were generated.

### A — magnitude (`step16 --line`, `step16fit`; `results/stage1/2026-10-28_scorer_space/tables/step16{.json,_heldout.csv,_rows_*.parquet,_curve_*.json,_zbeta_*.json}`)
- **Arms:**
  - **(i):** config 1's form. Called genes s_c × LFC; non-called genes m × c × Δ_pool, with c = z_a²/(z_a² + κ), c₁ for a single source, and 0 with no source.
  - **(j):** like (i), but the called genes are also multiplied by c (c₁ for one source or none).
  - **(k):** like (i), but c is a monotone curve of z_a, fixed at 0 at z = 0. For each held-out fold, the curve is fitted on the two training folds:
    - each line's per-decile β is normalised to its top decile;
    - folds are weighted equally, and the fit is isotonic in the decile's median z;
    - the curve is 1 past the top decile.

    Every line is scored with each fold's curve, so a fold's fit never sees its own β.
  - **(kj):** the curve applied on both sides, as in (j).
- **Grid:** m {0.5, 0.6, 0.75}, κ {0.25, 0.5, 1, 2}, s_c {0.35, 0.5, 0.65, 0.8}, c₁ {0.5, 1}. That is 336 rows per line.
- **Check:** config 1 (i: m 0.5, κ 1, s_c 0.5, c₁ 1) reproduces its Step 11 row on all four lines (|diff| 0).

**Held out (LOFO, three folds), scaled:**

| Arm | H1-fold parameters | K562 | HCT116 | HEK293T | **H1 (Eval A)** | H1: MSE / PDS / NMAE / fidelity / reach / Jaccard | Qualifies |
|---|---|---|---|---|---|---|---|
| i | m 0.6, κ 0.25, s_c 0.8, c₁ 0.5 | 0.170 | 0.320 | 0.166 | 0.224 | **0** / 0.790 / 0.012 / 0.371 / 0.077 / 0.095 | no (MSE 0) |
| j | m 0.75, κ 0.25, s_c 0.8, c₁ 0.5 | 0.167 | 0.307 | 0.155 | 0.229 | 0.009 / 0.789 / 0.030 / 0.371 / 0.082 / 0.095 | yes |
| k | m 0.6, s_c 0.8, c₁ 0.5 | 0.167 | 0.320 | 0.166 | 0.224 | **0** / 0.791 / 0.014 / 0.371 / 0.076 / 0.095 | no (MSE 0) |
| **kj** | m 0.6, s_c 0.8, c₁ 1 | 0.170 | 0.301 | 0.150 | **0.237** | 0.051 / 0.798 / 0.021 / 0.371 / 0.082 / 0.095 | **yes** |
| *Config 1 (Step 11, LOFO)* | *m 0.5, κ 1, s_c 0.5, c₁ 1* | *0.156* | *0.292* | *0.149* | ***0.244*** | *0.113 / 0.808 / 0.001 / 0.371 / 0.073 / 0.095* | *yes* |

- **By the Step 10 rule within Step 16, kj is carried (0.2365).** But **config 1 still has the higher held-out Eval A score (0.2436)**, and it is a point inside the new (i) grid. Widening the grid made the held-out fit worse on H1.
- **Why the fit can't transfer:** the training folds (K562, HCT116, HEK293T) keep MSE > 0 up to s_c 0.8 and score best there. On H1, s_c 0.8 cuts MSE to 0.07 or less, and to 0 at the LOFO picks. The H1 fold therefore inherits a called-gene scale that H1 can't take.
- **That's why LOFO passes over config 1.** On the H1 fold's training lines config 1 averages 0.191, against 0.217 for the LOFO pick.

**Fold spreads** (K562 / HCT116+HEK293T / H1 fold); **edge** = the value lies on a grid boundary:

| Arm | m | κ | s_c | c₁ |
|---|---|---|---|---|
| i | 0.75 **edge** / 0.6 / 0.6 (1.25) | 2 **edge** / 1 / 0.25 **edge** (**8**, fails ≤ 2) | 0.65 / 0.65 / 0.8 **edge** (1.23) | 0.5 **edge** ×3 (1) |
| j | 0.5 **edge** / 0.75 **edge** / 0.75 **edge** (1.5) | 0.25 **edge** / 1 / 0.25 **edge** (**4**, fails) | 0.65 / 0.8 **edge** / 0.8 **edge** (1.23) | 1 **edge** / 0.5 **edge** / 0.5 **edge** (2) |
| k | 0.5 **edge** / 0.5 **edge** / 0.6 (1.2) | — | 0.65 / 0.65 / 0.8 **edge** (1.23) | 1 **edge** / 0.5 **edge** / 0.5 **edge** (2) |
| kj | 0.6 / 0.5 **edge** / 0.6 (1.2) | — | 0.65 / 0.8 **edge** / 0.8 **edge** (1.23) | 1 **edge** / 0.5 **edge** / 1 **edge** (2) |

- κ is unstable across folds (8× and 4×).
- Every arm puts s_c at 0.8, the top of the grid, on at least one fold.
- c₁ is binary in this grid, so it is always on an edge.

**H1 in sample (not eligible):**
- **Best per arm:**
  - i: 0.250 (m 0.75, κ 1, s_c 0.5, c₁ 0.5);
  - j: 0.253;
  - k: 0.249;
  - kj: 0.255.
- **The (i) surface around config 1 is flat:** for κ 1, s_c 0.35–0.65 and every m, overall ranges 0.235–0.250. s_c 0.8 is where it falls (0.226–0.233).
- **Config 1 sits on an NMAE edge.** At config 1 NMAE is +0.001, and each of κ 2, c₁ 0.5 or s_c 0.35 makes it ≤ 0 (−0.006, −0.0001, −0.007). Raising m to 0.6 gives +0.008 (overall 0.247, MSE 0.127), and m 0.75 gives +0.019. That safety margin, though, is an H1-only pick.

**The learned curve (k)** is flatter than config 1's c at κ 1. At z ≈ 0.5 it gives 0.37–0.48, against 0.2. At z ≈ 1 it gives 0.59–0.65, against 0.5. It reaches 1 at z ≈ 2.5–4.4. It behaves like κ ≈ 0.25–0.5, and on its own it doesn't beat the κ form held out.

### B — calls (`step16calls`; `results/stage1/2026-10-28_scorer_space/tables/step16_calls{,_config1,_config1_ext}.json`, `…_rows*.parquet`)
- **Halves.** H1's 292 knockdowns are split into two halves of 146 (seed 2026). Knockdowns are paired by rank of true DE count, and one of each pair goes to each half. Median truth calls per knockdown are 218.5 and 218.
- **Fitting.** Each variant's budget is fitted on one half and scored on the other, then the halves are swapped. Halves use the full panel's scaling, as `boot_score` does, so compare variants within a half rather than across halves.
- **Magnitudes.** The magnitudes stay fixed; only the call set changes. I ran B on both:
  - the Step 16 carry (kj), as the task specifies;
  - config 1, because it still has the best held-out magnitude (A above).
- **Variants.**
  - K4: per knockdown, the top N tested genes by C+G z.
  - K5: the top N by agreement z among genes with ≥ 2 sources, then the rest by C+G z.
- **Extension.** K4 and K5 at 1,200 won on both halves, the edge of {400, 600, 800, 1200}. So I added a side check at {1,600, 2,400, 3,200, 4,800} on config 1. H1 has a median of 10,747 tested genes per knockdown.

**Config 1 magnitudes; H1 held-out halves (half 0 / half 1), scaled:**

| Variant | Median calls | Overall | Fidelity | Reach | Jaccard | HCT116 overall | HEK293T overall |
|---|---|---|---|---|---|---|---|
| **K0** | 495.5 / 447 | **0.264 / 0.224** | 0.371 / 0.371 | 0.102 / 0.044 | 0.111 / 0.079 | 0.300 | 0.154 |
| K4 400 | 400 | 0.259 / 0.222 | 0.382 / 0.371 | 0.100 / 0.044 | 0.075 / 0.073 | 0.299 | 0.142 |
| K4 600 | 600 | 0.269 / 0.232 | 0.420 / 0.410 | 0.104 / 0.044 | 0.094 / 0.090 | 0.294 | 0.141 |
| K4 800 | 800 | 0.276 / 0.239 | 0.448 / 0.438 | 0.103 / 0.044 | 0.108 / 0.103 | 0.292 | 0.141 |
| **K4 1200** | 1,200 | **0.286 / 0.250** | 0.489 / 0.480 | 0.100 / 0.048 | 0.128 / 0.122 | 0.291 | 0.140 |
| K5 400 | 400 | 0.249 / 0.209 | 0.370 / 0.360 | 0.106 / 0.030 | 0.020 / 0.017 | 0.299 | 0.139 |
| K5 800 | 800 | 0.265 / 0.226 | 0.440 / 0.430 | 0.105 / 0.036 | 0.050 / 0.046 | 0.293 | 0.138 |
| K5 1200 | 1,200 | 0.277 / 0.238 | 0.485 / 0.473 | 0.107 / 0.041 | 0.072 / 0.068 | 0.291 | 0.138 |
| *K4 1600 (ext)* | 1,600 | *0.294 / 0.258* | 0.520 / 0.510 | 0.104 / 0.051 | 0.144 / 0.138 | 0.291 | 0.138 |
| *K4 2400 (ext)* | 2,400 | *0.306 / 0.269* | 0.560 / 0.550 | 0.113 / 0.055 | 0.165 / 0.160 | 0.288 | 0.138 |
| *K4 3200 (ext)* | 3,200 | *0.311 / 0.274* | 0.580 / 0.567 | 0.110 / 0.058 | 0.178 / 0.174 | 0.289 | 0.138 |
| *K4 4800 (ext)* | 4,800 | *0.315 / 0.278* | 0.593 / 0.579 | 0.107 / 0.054 | 0.194 / 0.190 | 0.290 | 0.139 |

**Held-out fit:**
- **On config 1**, K4 picks 1,200 on both halves (spread 1, but on the grid edge). Its held-out mean is **0.268, against 0.244 for K0 (+0.024)**. K5 also picks 1,200, with a held-out mean of 0.257. Choosing across all variants picks K4 1,200 on both halves.
- **On the kj carry**, the result is the same: K4 1,200 on both halves, 0.260 against 0.236 for K0 (+0.024). K5 gives 0.250.

**Reading:**
- **Ranking by C+G z beats agreement z at every budget.** The gap is all Jaccard: K5's genes overlap the truth set far less. Fidelity and reach are about equal.
- **On H1, overall rises monotonically with the budget, all the way to 4,800 calls per knockdown** (45% of tested genes; held-out mean 0.297). The gain per step shrinks: about +0.004 from 3,200 to 4,800. Fidelity and Jaccard keep rising.
  - The scorer rewards calling far more genes than the truth holds: about 22× truth at 4,800.
- **The guardrails don't collapse.** HCT116 loses at most 0.012 and HEK293T at most 0.016 at any budget, because their Jaccard is already at about −1.3 to −1.6 from near-empty truth.
  - If D/E/F's truth is sparse, mass calling costs little. If it is rich like H1's, mass calling gains. The incentive is one-sided.
- **Config 1's NMAE is not robust to which knockdowns are scored:** +0.008 on half 0, −0.006 on half 1. Its +0.001 panel mean will very likely straddle zero in Step 15's paired bootstrap, which would fail the "every component > 0" rule.

### Step 16 recommendation (stop for Jose's approval before any cells are generated)
**Magnitude: keep config 1** (F0, arm i: called genes 0.5 × LFC, non-called 0.5 × z²/(z² + 1) × Δ_pool, single-source c₁ = 1). Held-out Eval A is 0.2436. Nothing in Step 16 beats it held out; the formal Step 16 carry, kj, scores 0.2365.

| Parameter | Folds (K562 / HCT116+HEK293T / H1) | Max/min | Grid position (Step 16 grid) |
|---|---|---|---|
| s_c (called scale) | 0.5 fixed | — | Interior; H1 is flat over 0.35–0.65 |
| m | 0.5 / 0.5 / 0.5 | 1 | Bottom edge of {0.5, 0.6, 0.75}; Step 11's grid edge was 0.5 from above |
| κ | 1 / 1 / 1 | 1 | Interior of {0.25, 0.5, 1, 2} |
| c₁ | 1 / 0.5 / 1 | 2 | Binary grid, always an edge |

- **Risk: NMAE (+0.001; the halves give +0.008 and −0.006).**
- **The obvious fix is m 0.6:** NMAE +0.008, overall 0.247 and MSE 0.127 on H1, with both halves positive to be checked. But it is an H1-only pick. Held out, the folds choose otherwise, because they pull s_c up rather than m.

**Call variant: K4 (top N by C+G z) at N = 1,200 per knockdown.**
- It is the held-out choice on both halves (spread 1) and adds about +0.024 on H1.
- It costs the guardrails 0.009 (HCT116) and 0.013 (HEK293T).
- **N is on the edge of the specified grid**, and the extension shows more calls keep scoring higher.

**Two things for Jose to decide before cells:**
1. **How far to push N.** The proxy rewards calling up to at least 45% of tested genes. Two reasons I would not go past 1,200 without a cell check:
   - **Generator realism:** the generator reaches calls by reshaping cell spread at a fixed mean. Step 9 validated the proxy at about 500 calls per knockdown, not 1,200–4,800, and many of the added genes carry near-zero requested shifts.
   - **Rules:** the 10-27 brief rates DE-significance shaping as our highest rules risk. Requesting thousands of significant genes per knockdown that the model has no evidence for leans hard on exactly that.
2. **Whether to accept the NMAE risk on config 1, or use m 0.6.** The alternative is to keep config 1 and let Step 15's bootstrap decide.

**Proposed Step 15, revised, still needing approval:** generate cells on Eval A and HCT116final for:
- config 1 + K0;
- config 1 + K4 1,200.

**What I would have done differently:**
- **Put s_c in the LOFO grid from Step 10 on.** The fold disagreement on the called scale is what blocks transfer, and it was visible only once s_c moved.
- **Run budgets out to "all tested genes" at once.** The response is monotone, so a short grid only produced an edge.
- **Check the NMAE margin under knockdown resampling before recommending config 1 in Step 14.** A half split takes seconds and would have shown +0.001 to be fragile.
- **Unit-test the fit summary before launching.** The `spreads()` bug cost one rerun of the fit (2 minutes; no rows lost).

## Step 15 — proxy checks, knockdown strength, end to end (task: Jose, 2026-10-30)
Outputs (tables, logs) in `results/stage1/2026-10-29_step15/`. Code: `code/vcc6/step15.py` (new module), run from `vcc_2026_8.ipynb`'s Step 15 section. The frozen generator, the harness and config 1's frozen files are imported and not changed.

**Rules:**
- 15a and 15b: noise-free proxy only, no cells.
- 15c: cells, only after Jose's go.
- No D/E/F file, no upload.

**Definitions used:**
- **Config 1′** (F0): called genes 0.5 × C+G LFC; non-called 0.6 × c × Δ_pool, c = z_a²/(z_a² + 1); single-source c₁ = 1; no source → control. That is config 1 with m 0.5 → 0.6.
- **K0:** C+G's calls. **K4 1,200:** per knockdown, the top 1,200 tested genes by C+G z. The magnitude split (LFC side vs Δ side) stays K0's; only the call set changes, as in Step 16 B.
- **H1 halves:** `scorer_step16.halves` (seed 2026, pairs by rank of true DE count). Halves are scored with the full panel's scaling, as `boot_score` does.

### 15b.1 — direction reach in cell-eval2 0.18.0
Paths relative to `~/venv/lib/python3.12/site-packages/cell_eval2/`.

**The scored member** is `de_wilcoxon_direction_reach_raw` = `de_direction_reach(universe="adjudicated", corrected=False)` (`catalog.py:649–651`), with the vcc2026 DE settings `p_adj_threshold 0.05`, `min_abs_log2fc 0.0`, `nan_lfc_policy mask` (`configs/vcc2026.yaml:43–48`).

**Definition**, per knockdown t (`metrics/direction.py:862`, body `:1014–1056`):

> reach(t) = k\*(t) / N_conf(t)

- **N_conf(t)** (`_reference_stats`, `:344–350`): the number of **reference-significant genes**, real `p_adj < 0.05`, with the target's own gene removed (`_exclude_on_target`). It counts significance only, not signs. N_conf = 0 → NaN (`:1052–1053`).
- **Ranking pool** (`_purity_curve`, `:470–471`): only reference-significant genes (`universe == "adjudicated"` → `filter(p_adj_real < alpha)`), own gene out, inner-joined with the prediction's rows.
- **Order** (`:508–511`): by the **prediction's** significance first (`p_adj_pred < 0.05` descending), then pred `p_adj`, pred `p_value`, pred `|log2FC|` descending, feature name.
- **Depth k** (`:513–520`): the running count of **adjudicable** pairs, i.e. real log2FC non-null, non-NaN, non-zero (`in_denom`, `:186`).
- **Match** (`:195–199`): adjudicable, the prediction has a direction (pred log2FC non-zero), and **sign(pred log2FC) = sign(real log2FC)**.
- **Purity(k)** = matches / k (`:521–523`). **k\*** is the **deepest** k with purity ≥ **P0 = 0.9** (`REACH_PURITY_FLOOR`, `:29`; `p0 = REACH_PURITY_FLOOR`, `:1043`; `_k_star`, `:553–562`). If no prefix reaches 0.9, k\* = 0.

**What reach uses:**

| Quantity | Used? |
|---|---|
| Genes | **All reference-DE genes** of the knockdown (real p_adj < 0.05, own gene out). Not the predicted call set: the prediction's significance only **orders** them. |
| Signs | Yes: the prediction's sign against the reference's sign, gene by gene. |
| Magnitudes | Only as a tie-break inside the ranking (pred \|log2FC\| after pred p_adj and p_value). The reference's magnitude is never used. |
| Thresholds | Reference p_adj < 0.05 (pool and N_conf); pred p_adj < 0.05 (first sort key); purity floor 0.9. No LFC floor. |

**What follows for the strength question:**
- Reach rewards putting the right sign on as many of the knockdown's **true** DE genes as possible, ahead of any wrong sign. Calling more genes helps only if the extra genes are true DE genes with the right sign, ranked before the first run of errors that drops purity below 0.9.
- A knockdown with many true DE genes has a large N_conf. Its reach rises only if the prediction calls and correctly signs a large share of them.

#### Addendum to 15b.1 — lfc_nmae, sig Jaccard, and p-value order in reach
**lfc_nmae** = `de_lfc_nmae` (`metrics/de.py:587`), registered with no override (`catalog.py:507–508`), so `min_gate_size = 10`. Per knockdown, nmae = mean |lfc_pred − lfc_real| / mean |lfc_real| over the gate (`:765–767`, `:785`).
- **Genes (the gate):** real `p_adj < 0.05` with a finite real log2FC, own gene removed (`:717–725`). The prediction's significance is **not used**: a gate gene missing from the prediction, or with a non-finite predicted LFC, counts as LFC 0 (`:737–740`, `:759–764`).
- **Knockdowns:** a knockdown is omitted if its gate is empty (`:773–774`), has **fewer than 10 genes** (`:779–781`), or has a zero or non-finite denominator (`:782–784`). All three are decided on the real side alone, so the same knockdowns are omitted for every submission. From cell-eval2's warnings in `results/stage1/2026-10-29_step15/logs/step15a.log`, NMAE is averaged over:
  - Eval A (H1): 216 of 292 knockdowns (6 empty gates, 70 with < 10 genes);
  - K562: about 78 of 267 (70 empty, 119 small);
  - HCT116: about 42 of 278 (92 empty, 144 small);
  - HEK293T: about 36 of 271.

  So the guardrail NMAEs rest on few knockdowns. The per-line counts are read from the warnings' order in the log, so treat them as approximate.
- **So:** call variants (K0 vs K4) cannot move NMAE. Only the magnitudes do, as the 15a table shows.

**sig Jaccard** = `de_sig_jaccard` (`metrics/de.py:176`; `catalog.py:669–670`, mean over knockdowns, no fill).
- Sets: R = real `p_adj < 0.05`, P = pred `p_adj < 0.05`, unique (target, feature), own gene removed from **both** sides (`:234–241`).
- **Both empty:** union 0 → **1.0** (`:250–256`). It is not masked; every knockdown gets a value.
- **One empty:** tp = 0 and union > 0 → **0.0**.
- That explains the guardrail lines' Jaccard (−1.1 to −1.6 scaled): their truth is mostly empty (median 1–3 genes). A knockdown with an empty truth scores 1 only if nothing is called for it, and every config here calls 113+ genes per knockdown.

**Does reach depend on the predicted p-values within the called set? Yes, through their order, not their values.**
- The ranking pool is every reference-DE gene, called or not. The sort is the called flag (`p_adj_pred < 0.05`), then pred `p_adj`, pred `p_value`, pred \|log2FC\| (`direction.py:506–511`). So within the called set, genes are ordered by their predicted p-values.
- The values matter only through that order and the 0.05 cut. Purity is evaluated at every prefix, and k\* is the **deepest** prefix with purity ≥ 0.9 (`:553`, `:562`). So reordering the called genes can move k\*. Only the purity at exactly k = \|called ∩ reference-DE\| is order-free (it equals precision there).
- **Uncalled reference-DE genes count too.** They rank after the called ones by their own predicted p-values, and extend k\* if purity stays ≥ 0.9. So reach is not limited to the call set.
- **In the proxy**, called genes get p = 1e-3/(1 + 1e4·s) and the rest 0.5 + 0.5/(1 + 1e4·s), with s = |P − control|. Within each group, the order is by scorer-space effect size. On real generated cells the order is the Wilcoxon p-values', which the generator shapes.

### 15a — proxy check of the planned configs (03:02–03:10 UTC; `results/stage1/2026-10-29_step15/tables/step15a.csv`, `results/stage1/2026-10-29_step15/logs/step15a.log`)
**Check:** config 1 + K0 and config 1 + K4 1,200 reproduce Step 16 B's rows on both H1 halves (|diff| 0).

**H1 (Eval A), scaled.** Halves use the full panel's scaling.

| Config | Calls | Part | Overall | MSE | PDS | **NMAE** | Fidelity | Reach | Jaccard | Calls (median) |
|---|---|---|---|---|---|---|---|---|---|---|
| **1′** | **K0** | full | **0.2470** | 0.127 | 0.804 | **+0.0081** | 0.371 | 0.077 | 0.095 | 486.5 |
| | | half 0 | 0.2679 | 0.145 | 0.857 | **+0.0147** | 0.371 | 0.109 | 0.111 | 495.5 |
| | | half 1 | 0.2261 | 0.109 | 0.751 | **+0.0014** | 0.371 | 0.046 | 0.079 | 447 |
| **1′** | **K4 1,200** | full | **0.2711** | 0.127 | 0.804 | **+0.0081** | 0.484 | 0.078 | 0.125 | 1,200 |
| | | half 0 | 0.2901 | 0.145 | 0.857 | **+0.0147** | 0.489 | 0.107 | 0.128 | 1,200 |
| | | half 1 | 0.2519 | 0.109 | 0.751 | **+0.0014** | 0.479 | 0.049 | 0.122 | 1,200 |
| 1 (ref.) | K0 | full / h0 / h1 | 0.2436 / 0.2636 / 0.2235 | 0.113 | 0.808 | +0.0010 / +0.0076 / **−0.0057** | 0.371 | 0.073 | 0.095 | 486.5 |
| 1 (ref.) | K4 1,200 | full / h0 / h1 | 0.2676 / 0.2857 / 0.2495 | 0.113 | 0.808 | same as K0 | 0.485 | 0.074 | 0.125 | 1,200 |
| C+G | K0 | full / h0 / h1 | 0.2173 / 0.2300 / 0.2046 | 0 | 0.724 | +0.040 / +0.052 / +0.029 | 0.371 | 0.074 | 0.095 | 486.5 |

**Guardrail lines, full panel** (fold pools, held out as before; note HCT116 here has both Orion lines held out, not the screen setup):

| Line | 1′ + K0 | 1′ + K4 1,200 | 1 + K0 | 1 + K4 1,200 | C+G |
|---|---|---|---|---|---|
| K562 | 0.158 | 0.140 | 0.156 | 0.139 | 0.146 |
| HCT116 | 0.311 | 0.301 | 0.300 | 0.291 | 0.124 |
| HEK293T | 0.154 | 0.141 | 0.154 | 0.140 | 0.106 |

Full six-component rows for every line are in `results/stage1/2026-10-29_step15/tables/step15a.csv`. On the guardrail lines every component except Jaccard is positive for config 1′.

**Reading:**
- **Config 1′ fixes config 1's NMAE sign on both halves**: +0.0147 and +0.0014, against +0.0076 and −0.0057. Half 1's margin is still thin, so a paired bootstrap on cells could straddle zero.
- **Config 1′ is at least as good as config 1 everywhere**: H1 +0.003 overall (MSE +0.014, PDS −0.004), HCT116 +0.011, K562 +0.002, HEK293T equal.
- **K4 1,200 adds +0.024 on H1** (fidelity +0.113, Jaccard +0.030). It costs the guardrails 0.009–0.018, all through Jaccard and a little fidelity. NMAE, MSE and PDS are unchanged, since they don't read the calls.
- **Every Eval A component is above zero for both planned configs, on the full panel and on each half.**

### 15b.2 — oracle ceiling (03:10–03:21; `results/stage1/2026-10-29_step15/tables/grid_oracle_p*.parquet`, `results/stage1/2026-10-29_step15/tables/step15b.json`)
**Setup** (Eval A, H1 halves; fit on one half, score on the other, both ways; base = config 1′ + K4 1,200):
- **True strength** = true DE count per knockdown (real p_adj < 0.05, own gene out). Terciles are taken over the whole panel; the halves are balanced by this count.
- **O1:** called and non-called magnitudes × s_k, with s_k per tercile from {0, 0.25, 0.5, 1, 1.5}³ (125 cells), on K4 1,200's calls.
- **O2:** N_k = clip(λ × true DE count, 0, 4,000), λ ∈ {0.5, 1, 2}, top N_k tested genes by C+G z, on config 1′ magnitudes.
- **O1 + O2:** the joint grid (375 cells).

**Held-out mean over the two directions (scaled):**

| Arm | Picked (both halves) | Overall | MSE | PDS | NMAE | Fidelity | Reach | Jaccard | Δ overall vs base |
|---|---|---|---|---|---|---|---|---|---|
| Base: 1′ + K4 1,200 | — | 0.2710 | 0.126 | 0.804 | +0.0081 | 0.484 | 0.078 | 0.125 | — |
| O1 (magnitude) | s = 0.5 / 0.5 / **1.5** (weak / mid / strong) | 0.2811 | 0.193 | 0.819 | **−0.0002** | 0.486 | 0.064 | 0.125 | **+0.010** |
| O2 (calls) | λ = **2** (median 436 calls) | 0.3151 | 0.126 | 0.804 | +0.0081 | 0.599 | 0.081 | 0.271 | **+0.044** |
| O1 + O2 | s = 0.5 / 0.5 / 1.5, λ = 2 | **0.3262** | 0.193 | 0.819 | −0.0002 | 0.600 | 0.074 | 0.271 | **+0.055** |

O2 per λ (the score on each half; not a fit):

| λ | Median calls | Overall, half 0 / half 1 | Fidelity | Jaccard |
|---|---|---|---|---|
| 0.5 | 109 | 0.263 / 0.225 | 0.24 | 0.20 |
| 1 | 218 | 0.329 / 0.290 | 0.59 | 0.25 |
| 2 | 436 | 0.337 / 0.293 | 0.60 | 0.27 |
| (K4 1,200, flat) | 1,200 | 0.290 / 0.252 | 0.48 | 0.12–0.13 |

**Reading:**
- **A perfect strength predictor is worth about +0.044 through the calls and +0.010 through the magnitudes; +0.055 together.**
- **The call gain comes from allocation, not volume.** O2 at λ = 2 asks for a median of 436 calls, a third of K4's 1,200, and still adds Jaccard +0.146 and fidelity +0.115. Even λ = 1 (218 calls, exactly the truth count) beats K4 1,200 by +0.04. Weak knockdowns stop drawing 1,200 false calls each; strong ones get up to 4,000.
  - *Correction (2026-10-30, Step 17):* the **median** is a third of K4's, but the **mean** is not: O2 λ 2 makes a mean of 1,390 calls per knockdown (λ 1: 956). The conclusion holds, though. Against flat K4 at the same mean volume, the oracle still gains +0.038 to +0.044 (section 9).
- **λ = 2 is on the grid edge**, but λ 1 → 2 is worth only +0.004–0.008, so the curve is flattening.
- **The magnitude gain is not usable as is.** O1 wants the strong tercile at 1.5× (grid edge) and the weak and mid at 0.5×. That raises MSE (+0.067) and PDS (+0.015), but **pushes NMAE to −0.0002** (−0.0054 on one direction) and costs reach −0.015. Under the adoption rule ("every component > 0"), O1 fails even with perfect strength.

### 15b.3 — strength predictors (03:02–03:04; `results/stage1/2026-10-29_step15/tables/step15b_predictors{.json,_H1.csv,_K562.csv}`)
Spearman ρ with the true DE count. Sources and the target's own controls only; the line's own source tables are excluded (H1: 111 tables; K562: 103, so the jiang K562 tables drop out).

| Predictor | H1 full | H1 half 0 | H1 half 1 | **H1 held-out (mean of halves)** | K562 (sources without K562) |
|---|---|---|---|---|---|
| (a) pooled predicted strength, J²-weighted ‖Δ_k‖² (Step 11 S1) | 0.224 | 0.293 | 0.153 | 0.223 | 0.054 |
| **(b) mean source DE count (n_real)** | **0.642** | 0.719 | 0.563 | **0.641** | **0.408** |
| (c) share of genes with agreement z > panel p90 (H1 2.89, K562 4.31) | −0.293 | −0.315 | −0.265 | −0.290 | −0.061 |
| (d) target CPM in the line's controls | 0.075 | −0.031 | 0.179 | 0.074 | 0.103 |
| (e) HGNC group: chromatin remodelling / TF / Mediator / TFIID / m6A | 0.247 | 0.301 | 0.193 | 0.247 | 0.117 |
| (f) ridge of (a)–(e), fit on the other half | 0.634 | 0.709 | 0.565 | 0.637 | 0.388 (fit on all H1) |

**Definitions used:**
- **(b):** per table (the line's own excluded; tables with `de_labels_known` False or no finite `n_real` skipped), the mean `n_real` over its rows for the knockdown; then the mean over the tables that hold it. 276 of 292 H1 knockdowns have a value; the rest get the median in 15b.4.
- **(c):** among each knockdown's tested genes with an agreement z (≥ 2 source lines), the share above the panel-wide 90th percentile.
- **(e):** HGNC has no single "transcription factor" group, so the rule is a keyword match on group names (`step15.HGNC_RULES`):
  - chromatin remodelling: BAF, PBAF, GBAF, INO80, NuRD, SRCAP, B-WICH, chromatin accessibility complex (8 groups);
  - TF: groups named "… transcription factor …" (not "General transcription factor"), homeoboxes, bHLH, bZIP / basic leucine zipper, forkhead, nuclear receptor subfamilies, SRY-box (69 groups); C2H2 zinc fingers are **not** included;
  - Mediator complex; TFIID (GTF IID subunits, TBP, TAF1 family); m6A methyltransferase complex.
  - H1 hits: 9 chromatin, 16 TF, 9 Mediator, 2 TFIID, 2 m6A. Single-group ρ on H1: 0.07–0.14.
- **(f):** ridge (α = 1) on standardised log1p(a), log1p(b), c, log1p(d), e, target log1p(true count). β is stable across halves: (b) 1.60 / 1.32, (c) −0.34 / −0.32, (e) 0.41 / 0.25, (a) and (d) near 0.

**Reading:**
- **(b) is the only strong predictor:** held-out ρ 0.64 on H1 and 0.41 on K562. A knockdown that moves many genes in other lines tends to move many in H1. The ridge adds nothing to it.
- **(a), the Δ-energy predictor Step 11 used, is weak** (0.22). That is why Step 11's strength arm (f) collapsed to a uniform m.
- **(c) is negative:** knockdowns where many genes agree strongly across sources are, if anything, weaker in H1. My guess is that strong cross-source agreement marks a shared, generic response (stress, proliferation) rather than a knockdown-specific one. I haven't tested that.
- **Best held-out predictor: (b).** It goes to 15b.4.

#### Addendum O3 — the model's own expected DE count (`results/stage1/2026-10-29_step15/tables/grid_o3_fit*_p*.parquet`, `step15b_o3_counts.json`)
- **Budget:** N_k = round(Σ_g P̂(DE | C+G z)) over the knockdown's tested genes, ranked by C+G z.
- **The map:** P̂ is fitted on the fit half's tested (k, g): 200 quantile bins of z, isotonic in the bin medians. It is applied to the other half, on config 1′ magnitudes.
- **O1 + O3** also uses the oracle's true-strength terciles.

| Arm (held-out mean) | Median calls | Overall | MSE | PDS | NMAE | Fidelity | Reach | Jaccard | Δ vs base |
|---|---|---|---|---|---|---|---|---|---|
| Base: 1′ + K4 1,200 | 1,200 | 0.2710 | 0.126 | 0.804 | +0.0081 | 0.484 | 0.078 | 0.125 | — |
| **O2** (true count, λ 2) | 436 | 0.3151 | 0.126 | 0.804 | +0.0081 | 0.599 | 0.081 | 0.271 | **+0.044** |
| **O3** (expected count) | 1,050 | 0.2681 | 0.126 | 0.804 | +0.0081 | 0.471 | 0.079 | 0.120 | **−0.003** |
| O1 + O3 (true terciles) | 1,050 | 0.2786 | 0.193 | 0.819 | −0.0002 | 0.473 | 0.068 | 0.120 | +0.008 |

- **O3 does not help:** its N_k barely tracks the truth. Held-out Spearman of N_k with the true count is 0.15 and 0.29, against 0.64 for predictor (b). C+G's z is inflated almost evenly across knockdowns, so the expected count comes out near 1,050 for nearly every knockdown, close to a flat budget, and slightly worse than K4's 1,200.
- So "how many genes the model is confident about" is not a usable strength signal. The source DE count is.

### 15b.4 — predictor (b) in place of true strength (`results/stage1/2026-10-29_step15/tables/grid_pred_b_fit*_p*.parquet`, `grid_o3p_b_*`, `grid_o3s_b_*`)
**Setup:**
- **Terciles:** by (b)'s rank, which needs no truth.
- **Predicted DE count for O2:** (b)'s panel percentile mapped onto the **fit half's** true-count quantiles. So N_k is fitted on one half and scored on the other.
- **O3 with (b)**, run two ways:
  - **O1 + O3**, with (b)'s terciles;
  - **O3s**, with the z → P(DE) map fitted separately within each (b) tercile (alone, and with O1).

**Held-out mean (scaled):**

| Arm | Picked | Median calls | Overall | MSE | PDS | NMAE | Fidelity | Reach | Jaccard | Δ vs base | Share of the oracle's gain |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Base: 1′ + K4 1,200 | — | 1,200 | 0.2710 | 0.126 | 0.804 | +0.0081 | 0.484 | 0.078 | 0.125 | — | — |
| O1 (b) | s = 0.5 / **1.5** / 1 (weak / mid / strong), both halves | 1,200 | 0.2747 | 0.140 | 0.817 | +0.0052 | 0.485 | 0.076 | 0.125 | **+0.004** | 37% of O1's +0.010 |
| O2 (b) | λ 2, both halves | 437 | 0.2659 | 0.126 | 0.804 | +0.0081 | 0.417 | 0.087 | 0.153 | **−0.005** | −12% of O2's +0.044 |
| O1 + O2 (b) | as above | 437 | 0.2698 | 0.140 | 0.817 | +0.0052 | 0.417 | 0.087 | 0.153 | −0.001 | −2% of +0.055 |
| O1 + O3 (b terciles) | s 0.5 / 1.5 / 1 | 1,050 | 0.2720 | 0.140 | 0.817 | +0.0052 | 0.472 | 0.078 | 0.120 | +0.001 | 2% of +0.055 |
| **O3s** (map per (b) tercile) | — | 980 / 1,199 | 0.2743 | 0.126 | 0.804 | +0.0081 | 0.484 | 0.082 | **0.142** | **+0.003** | 8% of O2's +0.044 |
| **O1 + O3s** | s 0.5/1.5/1 · 0.5/0.5/1.5 | 980 / 1,199 | **0.2757** | 0.134 | 0.811 | +0.0074 | 0.485 | 0.076 | 0.142 | **+0.005** | 9% of +0.055 |

(The shares for O3 rows are against O2 and O1 + O2. `step15b.json`'s automatic share for the O3s rows divides by the pooled O3 "oracle", whose own gain is −0.003, so that ratio is meaningless.)

**Side check, outside the task's grid (`lamext`; `results/stage1/2026-10-29_step15/tables/side_lambda.parquet`): O2 at λ 3, 4, 6.** λ = 2 sat on the grid edge, and with predicted counts it also cuts the volume to a third of K4's.

| λ | Median calls | Oracle: half 0 / half 1 | (b), fitted on the other half: half 0 / half 1 | (b) Δ vs base (mean) |
|---|---|---|---|---|
| 2 (grid) | 437 | 0.337 / 0.293 | 0.289 / 0.243 | −0.005 |
| 3 | 655 | 0.336 / 0.293 | 0.297 / 0.252 | +0.003 |
| 4 | 873 | 0.333 / 0.291 | 0.302 / 0.256 | +0.008 |
| 6 | 1,310 | 0.331 / 0.290 | **0.307 / 0.262** | **+0.013** |
| base (K4 1,200, flat) | 1,200 | 0.290 / 0.252 | 0.290 / 0.252 | — |

- **The oracle peaks at λ 2–3 and then falls slowly.** With perfect strength, about 2–3× the truth count is right.
- **With (b), more volume keeps helping up to λ 6**, which is still an edge. At K4's volume, (b)'s allocation beats flat K4 by **+0.010 and +0.017 on the two halves**. That is about 30% of the oracle's call gain. It comes from Jaccard (+0.05) and fidelity (+0.02).
- **Caveat for D/E/F:** this needs a DE-count scale for the target line. Here it comes from H1's own true-count quantiles (fitted on the other half). D/E/F have no truth, so the shape would have to be borrowed from H1 or from the sources. That is an assumption the H1 proxy can't test.
- **Correction (2026-10-30, Step 17, section 9):** "at K4's volume" above means the same **median** calls. An allocation rule concentrates calls on strong knockdowns, so its **mean** (total) volume is much higher: (b) at λ 6 makes a mean of ~1,910 calls per knockdown against flat K4's 1,200. Against flat K4 at the same mean volume, (b)'s gain is **−0.002 / +0.003**, about zero. The +0.013 was mostly volume. The oracle's call gain survives mean matching (+0.04 at λ 1–2).

**Reading of 15b.4:**
- **A source-only strength predictor recovers little within the task's grid:** at most +0.005 (O1 + O3s), with every component still positive.
- **Its magnitude half (O1) recovers 37% of a small oracle gain** (+0.004 of +0.010) and keeps NMAE positive (+0.005). The fit wants the *mid* tercile at 1.5×, not the strong one, which suggests (b)'s terciles don't line up with true strength at the top.
- **Its call half helps only at K4's volume or above (side check, +0.013).** At the oracle's preferred volume it loses, because a wrong budget costs more than a flat one.
  - *Correction (2026-10-30, Step 17):* matched on mean volume rather than median, this +0.013 is about zero (section 9). The "+0.013" in the stop summary below has the same caveat.

### Stop — summary and recommendation for 15c (2026-10-30 04:15 UTC)

**15a:**
- Config 1′ scores **0.2470 (K0)** and **0.2711 (K4 1,200)** on Eval A, against C+G's 0.2173.
- Every component is positive on the full panel and on both halves.
- **NMAE is now +0.0147 / +0.0014 on the halves** (config 1: +0.0076 / −0.0057), positive on both but still thin on half 1.
- The guardrail lines are at or above config 1. K4 costs them 0.009–0.018.

**15b.1:** reach = k\*/N_conf over **all reference-DE genes**, ordered by the prediction's significance. It reads signs only; magnitude is a tie-break. Purity floor 0.9, deepest qualifying prefix.
- The order of the predicted p-values within the called set matters; their values don't, beyond the 0.05 cut.
- NMAE ignores the calls. It drops knockdowns with < 10 real-DE genes (Eval A keeps 216 of 292).
- Jaccard scores ∅/∅ = 1 and one-side-empty = 0, with no mask.

**15b.2 (oracle ceiling over 1′ + K4 1,200):**
- Calls by true count (O2) **+0.044**.
- Magnitudes by true tercile (O1) **+0.010**, but NMAE turns negative.
- Both together **+0.055**.

**15b.3:** only **(b), the mean source DE count**, predicts strength: ρ 0.64 held out on H1 and 0.41 on K562. Ridge adds nothing; Step 11's Δ-energy measure is 0.22.

**15b.4:**
- Within the grid, (b) recovers at most **+0.005** (O1 + O3s), about 9% of the oracle's gain.
- **Off-grid, (b)-allocated calls at K4's volume give +0.013**, still rising at the λ edge.
- O3, the model's own expected count, doesn't track strength (ρ 0.15–0.29) and gives −0.003.

**Recommendation: keep the default pair for 15c.**
- **config 1′ + K0 and config 1′ + K4 1,200.**
- **Why not swap in a strength variant:**
  - Their in-grid proxy gains (+0.003 to +0.005) are below the paired-bootstrap resolution we've seen on Eval A (CI half-widths around 0.01–0.015).
  - The larger +0.013 comes from an off-grid λ at an edge, with a D/E/F count scale that can't be checked.
- **What the pair settles on cells:**
  - (1) Does config 1′'s thin NMAE (+0.0014 on half 1) survive generation?
  - (2) Does the proxy's +0.024 for 1,200 flat calls survive the generator's spread-shaping (Step 9 validated the proxy only at ~500 calls)?
- If K4 1,200 holds up on cells, (b)-allocated budgets at the same volume are the natural next proxy step. They would add calls only where the sources say the knockdown is strong.
- **Not recommended:** O1-style magnitude scaling by strength. Even the oracle version turns NMAE negative.

**What I would have done differently:**
1. **Put λ past 2 in the O2 grid from the start.** λ = 2 was an edge for the oracle and for (b), and the useful result for (b) is at λ 6. A proxy row costs seconds; I repeated Step 16's lesson here.
2. **Matched call volume when comparing allocation rules.** O2 with predicted counts at λ 2 mixed "where the calls go" with "how many calls", and lost on volume. The side check had to separate the two afterwards.
3. **Checked a predictor's own count against the truth before running its grid.** O3's N_k has ρ 0.15–0.29 with the truth, which predicts its result without 250 proxy rows.

**Waiting for Jose's go before 15c.**

---

# 7. 2026-10-30 — Illustrations: profiles, changes and DE calls (H1, HCT116, context D)
Jose asked for these before 15c. Code: `code/vcc6/illustrations.py`. Figures and `stats.json` are in `results/stage1/2026-10-30_illustrations/`.

**Data:** the scorer's own real side, i.e. the harness real-cache pseudobulk, log1p of CP50K.
- H1 = Eval A: 292 knockdowns, 400 cells each.
- HCT116 = HCT116final: 278 knockdowns, 146–320 cells each.
- Context D: its 18,400 control cells, pseudobulked the same way. D has no knockdown truth.
- Truth = each eval's real DE table, p_adj < 0.05, own gene out.
- Calls = the incumbent C+G end-to-end runs as the scorer saw them (harness `de_pred`, p_adj < 0.05).

## Figure 1 — control vs knockdown profiles
![log1p profiles](../2026-10-30_illustrations/fig1_log1p_profiles.png)

- **Grey:** the control profile, one dot per gene, with genes ranked by it.
- **Orange:** every knockdown's profile on the same genes (a random quarter of the 5.4 M dots, to keep the cloud readable).
- **Reading:**
  - Knockdown profiles hug the control curve. On the scale of the whole expression range, knockdowns move genes very little; the effects are fine structure on top of a shared profile.
  - HCT116's cloud is visibly wider than H1's. That is mostly sampling noise: its knockdowns have 146–320 cells against H1's 400. It is not larger effects (figure 2).
  - D's control curve has the same shape. More genes sit at zero in D (2,975) than in the H1 controls (1,422) or HCT116's (2,088).

## Figure 2 — how far, and which way, knockdowns move genes
![delta](../2026-10-30_illustrations/fig2_delta.png)

**Unit: Δ = knockdown log1p − control log1p**, in the scorer's space, not log fold change.
- Δ is what the MSE metric scores, and it shows the size of a change where it matters.
- LFC blows up on barely expressed genes: a ratio of two tiny CPMs. That would put the most extreme points on the genes that carry the least weight (Step 3b found C+G's largest LFCs on genes with control CPM ≈ 0).
- Gene order is as in figure 1; both Δ panels share one y-scale.
- Right-hand panels: true DE genes per knockdown, sorted, up above and down below, on a symlog axis.

**Reading:**
- **H1 has rich truth:** 306,169 true DE (knockdown, gene) pairs, 46% up and 54% down. The median knockdown has 218 true DE genes (mean 1,049, so a few strong ones dominate), and only 6 of 292 have none.
- **H1's effects are small:** true DE genes move by a median |Δ| of 0.14, against 0.014 for the rest. DE colour starts near gene rank 7,700 because the real DE test only covers genes expressed in the controls.
- **HCT116's truth is almost empty:** 2,943 true DE pairs in total, a median of 1 per knockdown, and 92 of 278 knockdowns with none. Its grey cloud is mostly noise (median |Δ| 0.024 for non-DE genes, larger than H1's).
- The curved grey band down to about −0.4 between ranks 5,000 and 9,000 is low-expression genes that read zero in knockdowns with few cells (Δ = −control). That is noise too.
- **Both lines split roughly evenly between up and down.** Direction has to be predicted per gene; there is no "mostly down" shortcut.

## Figure 3 — DE calls against the truth (incumbent C+G)
![DE calls](../2026-10-30_illustrations/fig3_de_calls.png)

**Layout:**
- Per knockdown, sorted by true DE count (strongest left).
- Above 0, what was called: correct (true DE, same sign), wrong direction (true DE, opposite sign), false call (not true DE).
- Below 0, missed (true DE, not called).
- Right: totals over all knockdowns.

| | Called | Correct | Wrong direction | False call | Missed | Precision | Recall | Sign right on hits | Median per knockdown: called / true DE / correct / false / missed |
|---|---|---|---|---|---|---|---|---|---|
| **H1 (Eval A)** | 171,056 | 26,260 | 7,756 | 137,040 | 272,153 | 0.20 | 0.11 | 77% | 520 / 218 / 16.5 / 384.5 / 188 |
| **HCT116** | 53,219 | 788 | 76 | 52,355 | 2,079 | 0.02 | 0.29 | 91% | 126 / 1 / 0 / 124 / 1 |

**Reading:**
- **The call count doesn't follow the truth.** On H1, the strongest knockdowns (left) have thousands of true DE genes, and we call a few hundred of them, so most are missed. The weakest (right) have almost none, and still get several hundred calls each, all false. That is the pattern Step 15's oracle O2 exploits (+0.044): move calls from weak knockdowns to strong ones.
- **Where we do hit, the sign is mostly right** (77% on H1, 91% on HCT116). The wrong-direction slice is small next to the false calls and the misses.
- **On HCT116 nearly every call is false** (98%), because the truth is nearly empty. This is the near-noise-truth situation behind the guardrail lines' negative Jaccard. Calling less would only help there.

## Figures 4–5 — change against control expression (Jose's follow-up, 2026-10-30)
**What they are:** one file per line and unit, from the same harness real side as figures 1–2.
- **x:** the gene's control log1p (pseudobulk CP50K, the scorer's space), as a value rather than a rank.
- **y:** every knockdown × gene pair's Δ (figure 4) or log2 fold change (figure 5).
- **The LFC:** from the same pseudobulk, log2((kd + ε) / (ctrl + ε)) with ε = 0.01 CP50K, so it covers all 18,533 genes. The scorer's Wilcoxon table has LFCs only for the genes it tests (10,778 on H1, 10,320 on HCT116).
- **Colours:** true DE up (red) and down (blue) as in figure 2; tested but not DE (grey, 15% shown); not in the scorer's DE table (light grey, 15% shown).
- **Code:** `code/vcc6/illustrations.py --only fig4 fig5`. The numbers below are in `stats.json` under `fig4_fig5`.

![H1, delta against control expression](../2026-10-30_illustrations/fig4_H1_delta_vs_control.png)
![HCT116, delta against control expression](../2026-10-30_illustrations/fig4_HCT116_delta_vs_control.png)
![H1, log2 fold change against control expression](../2026-10-30_illustrations/fig5_H1_lfc_vs_control.png)
![HCT116, log2 fold change against control expression](../2026-10-30_illustrations/fig5_HCT116_lfc_vs_control.png)

**Reading:**
- **The two units weight expression levels in opposite ways.**
  - **Δ:** on H1 the 99th percentile of |Δ| is 0.33 among the top tenth of tested genes by control expression, and 0.11 in the bottom tenth.
  - **LFC:** the order flips: 0.51 at the top, 1.46 at the bottom.
  - log1p compresses low counts, so a doubling of a barely-expressed gene is a small Δ but a full unit of LFC.
- **LFC fans out at low expression.** For genes outside the scorer's DE table the 99th percentile of |LFC| is 1.5 (H1) and 2.8 (HCT116), almost all of it noise. That is why figure 2 used Δ.
- **The scorer's DE gate is visible as a vertical edge** at x ≈ 0.22, which is 5 CPM. Left of it everything is light grey: those genes are never tested for DE, so only MSE and PDS see them.
- **H1's true DE pairs sit off the zero line at every expression level.** Their median |Δ| is 0.14 and median |LFC| 0.29, against 0.016 and 0.088 over all pairs. The real changes are small, and present across the whole range.
- **HCT116 has almost no true DE** (2,943 pairs). Its grey cloud is mostly sampling noise, which is wider than H1's because HCT116's knockdowns have fewer cells.
- **The diagonal lower edge** in each panel is a knockdown with zero counts for a gene: there Δ = −control, and LFC = log2(ε / (control + ε)).

## Figure 6 — the same knockdown in two lines (Jose's follow-up, 2026-10-30)
**What it is:** 10 of the 27 knockdowns in both Eval A (H1) and HCT116final, chosen to span H1's true DE count and to include the two with the most HCT116 DE genes (SHPRH, ACLY).
- **Axes:** x = Δ in H1 and y = Δ in HCT116, for all 18,533 genes. Each line's Δ is knockdown log1p minus its own controls' log1p, in pseudobulk CP50K (the scorer's space, as MSE and PDS see it).
- **Colours:** true DE in both lines (orange), in H1 only (blue), in HCT116 only (aqua), in neither (grey).
- **The knocked-down gene:** a star, pinned to the panel edge with its values printed when off-scale. It is left out of the statistics, as the scorer leaves it out.
- **Files:** a grid `fig6_cross_line_delta_grid.png` and one file per knockdown `fig6_<target>_H1_vs_HCT116_delta.png`; numbers in `stats.json` → `fig6`.

![Same knockdown, H1 vs HCT116 delta](../2026-10-30_illustrations/fig6_cross_line_delta_grid.png)

| Knockdown | true DE H1 / HCT116 / both | Pearson r | cosine | slope (HCT116 on H1) | ‖Δ_HCT116‖ / ‖Δ_H1‖ | target Δ H1 / HCT116 |
|---|---|---|---|---|---|---|
| STUB1 | 2 / 2 / 0 | 0.01 | 0.01 | 0.04 | 3.17 | −1.84 / −1.18 |
| TARBP2 | 26 / 1 / 0 | 0.01 | 0.01 | 0.01 | 1.42 | −1.19 / −0.61 |
| SHPRH | 37 / 175 / 3 | 0.08 | 0.09 | 0.20 | 2.38 | −1.20 / −0.72 |
| PLAGL2 | 100 / 1 / 0 | 0.01 | 0.01 | 0.02 | 1.57 | −1.42 / −0.51 |
| SLIRP | 251 / 11 / 3 | 0.03 | 0.03 | 0.07 | 2.52 | −1.72 / −1.31 |
| STAT6 | 323 / 13 / 5 | 0.05 | 0.05 | 0.06 | 1.12 | −0.22 / −0.55 |
| MED13 | 808 / 1 / 0 | −0.01 | −0.01 | −0.01 | 1.12 | −1.58 / −0.63 |
| ACLY | 1,247 / 290 / 41 | 0.12 | 0.12 | 0.17 | 1.43 | −2.39 / −1.02 |
| BRPF1 | 2,418 / 20 / 6 | 0.06 | 0.06 | 0.06 | 0.95 | −0.79 / −0.39 |
| SIN3B | 5,324 / 26 / 22 | 0.21 | 0.21 | 0.20 | 0.93 | −1.51 / −0.90 |

**Reading:**
- **The same knockdown barely transfers between these two lines in Δ.** The cross-line correlation over all genes is 0.01–0.21; the best case is SIN3B, the strongest H1 knockdown.
  - That doesn't mean the biology differs completely. Most genes' Δ is sampling noise, and the noise is independent between lines, which caps the correlation.
  - Where both lines call a gene DE (orange), the sign agrees: 3/3 for SHPRH and SLIRP, 5/5 for STAT6, 34/41 for ACLY, 20/22 for SIN3B, and 4/6 for BRPF1.
- **Magnitudes don't transfer either.** The slope of HCT116 on H1 is 0.0–0.2, and the clouds are near-perpendicular: H1's DE genes (blue) spread along x at y ≈ 0, HCT116's (aqua) along y at x ≈ 0. A model copying H1's Δ onto HCT116 would put the wrong genes in the wrong places for MSE. For PDS, which ranks knockdowns by cosine, cosines of 0.01–0.2 between lines show how little of one line's direction survives into the other.
- **HCT116's cloud is wider** (norm ratio up to 3.2× on weak knockdowns). That's noise from its smaller cells-per-knockdown, not a bigger effect, as figure 2 showed.
- **The knocked-down gene itself is the one shared signal:** it's down in both lines for all 10, usually by more in H1 (e.g. ACLY −2.39 vs −1.02).
- **The horizontal streak at y = 0** is genes with no HCT116 counts in either the knockdown or the controls, so their Δ is exactly 0 there.
- **What this says about the approach:** our Stage-1 model averages Δ across source lines. These panels show why the shrunk, agreement-weighted pool (config 1′) helps MSE, while copying any single line's magnitudes would not.

---

# 8. 2026-10-30 — Step 15c: config 1′ on cells (Eval A, HCT116final)
*Task: Step 15c, with the per-tercile split from the addendum. Go from Jose, 2026-10-30 ~06:05 UTC. Code: `code/vcc6/step15c.py` and `step15c.sh`. Outputs: `results/stage1/2026-10-29_step15/c/<config>/<line>/`, comparison tables in `c/tables/` (`step15c_scores.csv`, `step15c_paired.csv`, `step15c.json`).*

**Checks:**
- **Model:** on Eval A, the model code with m 0.5 reproduces v7's freeze-check input exactly (max |ΔLFC| 0). Calls, LFC and SE equal the scorer-space fold cache exactly.
- **Harness:** every run's rebuilt aggregate reproduces its `summary.json` (|diff| < 1e-14), and cell-eval2's DE dispatch re-run on its own `de_pred` reproduces `results.csv`.
- **Generator:** it reports a median of 1,201 realised calls for K4's 1,200 requested.

## 1. Six components, scaled and raw, with paired bootstrap CIs
2,000 knockdown resamples (`boot_score`), with each eval's own bundle scaling.

**Scores (scaled, with 95% CI; raw in brackets):**

| Line | Run | **Overall** | MSE | PDS | NMAE | Fidelity | Reach | Jaccard |
|---|---|---|---|---|---|---|---|---|
| Eval A | C+G (incumbent) | **0.2186** [0.204, 0.233] | 0 (1.496) | 0.723 (0.834) | 0.040 (0.927) | 0.380 (0.378) | 0.070 (0.184) | 0.097 (0.052) |
| Eval A | 1′ + K0 | **0.2426** [0.227, 0.258] | 0.127 (0.845) | 0.795 (0.867) | +0.0077 (0.946) | 0.372 (0.371) | 0.060 (0.175) | 0.094 (0.050) |
| Eval A | 1′ + K4 1,200 | **0.2638** [0.250, 0.278] | 0.127 (0.845) | 0.795 (0.867) | +0.0078 (0.946) | 0.483 (0.457) | 0.048 (0.165) | 0.123 (0.062) |
| HCT116final | C+G (incumbent) | **0.1903** [0.144, 0.235] | 0 (2.896) | 0.655 (0.628) | 0.401 (0.747) | 0.915 (0.540) | 0.741 (0.668) | −1.571 (0.0145) |
| HCT116final | 1′ + K0 | **0.3617** [0.303, 0.403] | 1.0 (0.825) | 0.860 (0.662) | 0.224 (0.859) | 0.918 (0.542) | 0.738 (0.666) | −1.570 (0.0146) |
| HCT116final | 1′ + K4 1,200 | **0.3476** [0.291, 0.386] | 1.0 (0.825) | 0.860 (0.662) | 0.224 (0.859) | 0.899 (0.531) | 0.731 (0.662) | −1.628 (0.0038) |

**Paired differences** (scaled; 95% CI):

| Line | Comparison | **Overall** | MSE | PDS | NMAE | Fidelity | Reach | Jaccard |
|---|---|---|---|---|---|---|---|---|
| Eval A | K0 − C+G | **+0.024 [+0.016, +0.033]** | +0.127 [+0.104, +0.150] | +0.072 [+0.038, +0.105] | −0.033 [−0.045, −0.020] | −0.009 [−0.012, −0.005] | −0.010 [−0.031, +0.013] | −0.003 [−0.005, −0.002] |
| Eval A | K4 − C+G | **+0.045 [+0.035, +0.055]** | +0.127 [+0.104, +0.150] | +0.072 [+0.038, +0.105] | −0.033 [−0.045, −0.020] | +0.102 [+0.080, +0.126] | −0.022 [−0.047, +0.005] | +0.025 [+0.001, +0.044] |
| Eval A | K4 − K0 | **+0.021 [+0.015, +0.028]** | 0.000 | 0.000 | +0.000 | +0.111 [+0.088, +0.134] | −0.012 [−0.028, +0.002] | +0.029 [+0.004, +0.048] |
| HCT116final | K0 − C+G | **+0.171 [+0.127, +0.206]** | +1.000 [+0.853, +1.000] | +0.206 [+0.031, +0.380] | −0.178 [−0.212, −0.144] | +0.003 [+0.001, +0.006] | −0.003 [−0.091, +0.089] | +0.001 [−0.001, +0.002] |
| HCT116final | K4 − C+G | **+0.157 [+0.114, +0.191]** | +1.000 [+0.853, +1.000] | +0.206 [+0.031, +0.380] | −0.178 [−0.212, −0.144] | −0.016 [−0.028, −0.002] | −0.010 [−0.096, +0.073] | −0.058 [−0.071, −0.046] |
| HCT116final | K4 − K0 | **−0.014 [−0.022, −0.007]** | 0.000 | 0.000 | 0.000 | −0.019 [−0.031, −0.006] | −0.007 [−0.050, +0.031] | −0.058 [−0.072, −0.047] |

**Reading:**
- **Config 1′'s magnitudes:** +0.024 on Eval A and +0.171 on HCT116final over C+G, entirely through MSE and PDS.
  - **NMAE drops against C+G** (−0.033 on Eval A, −0.178 on HCT116) but stays positive: +0.0077 on Eval A, with a CI of [−0.004, +0.019] that crosses 0.
  - **HCT116's MSE clamps at 1.0** (raw 0.825 against C+G's 2.896).
- **K4 1,200 adds +0.021 on Eval A** (fidelity +0.111, Jaccard +0.029, reach −0.012, ns) but **costs −0.014 on HCT116final**, almost all of it Jaccard (−0.058). On a line with a median of 1 true DE gene per knockdown, 1,200 calls are almost all false.
- **K4 changes only the calls:** MSE, PDS and NMAE are identical between K0 and K4 to 1e-3, as the proxy said.

## 1b. By true-strength tercile (Eval A; terciles by true DE count; paired CIs within each tercile, full-panel scaling)
Terciles: weak 0–26 true DE genes (median 5), mid 26–862 (median 215), strong 870–7,343 (median 2,490); 97 / 97 / 98 knockdowns.

| Tercile | Comparison | **Overall** | MSE | PDS | NMAE | Fidelity | Reach | Jaccard |
|---|---|---|---|---|---|---|---|---|
| weak | K0 − C+G | +0.193 [+0.070, +0.216] | +0.971 [+0.254, +1.0] | +0.196 [+0.125, +0.273] | −0.016 [−0.058, +0.028] | +0.002 [−0.004, +0.008] | +0.005 [−0.068, +0.073] | +0.001 [+0.000, +0.002] |
| weak | K4 − C+G | +0.177 [+0.052, +0.210] | +0.929 [+0.213, +1.0] | +0.196 | −0.016 | +0.016 [−0.007, +0.042] | −0.021 [−0.104, +0.055] | −0.043 [−0.106, −0.004] |
| weak | K4 − K0 | **−0.016 [−0.026, +0.000]** | −0.042 [−0.046, 0.000] | 0 | 0 | +0.014 [−0.008, +0.039] | −0.026 [−0.068, +0.014] | −0.044 [−0.107, −0.005] |
| mid | K0 − C+G | +0.011 [−0.004, +0.032] | +0.067 [0.000, +0.182] | +0.040 [−0.009, +0.087] | −0.007 [−0.028, +0.014] | −0.013 [−0.019, −0.007] | −0.019 [−0.033, −0.007] | −0.001 [−0.003, +0.001] |
| mid | K4 − C+G | +0.026 [+0.007, +0.050] | +0.066 | +0.040 | −0.007 | +0.080 [+0.040, +0.119] | −0.019 [−0.033, −0.007] | −0.003 [−0.011, +0.005] |
| mid | K4 − K0 | **+0.015 [+0.007, +0.023]** | −0.001 | 0 | 0 | +0.092 [+0.051, +0.133] | −0.001 [−0.011, +0.011] | −0.002 [−0.011, +0.006] |
| strong | K0 − C+G | **+0.000 [−0.009, +0.008]** | +0.121 [+0.099, +0.136] | −0.020 [−0.064, +0.021] | **−0.061 [−0.075, −0.049]** | −0.015 [−0.018, −0.012] | −0.015 [−0.022, −0.008] | −0.010 [−0.013, −0.008] |
| strong | K4 − C+G | +0.058 [+0.045, +0.070] | +0.122 | −0.020 | −0.061 | +0.210 [+0.170, +0.250] | −0.026 [−0.036, −0.016] | +0.121 [+0.098, +0.145] |
| strong | K4 − K0 | **+0.057 [+0.048, +0.067]** | +0.001 | 0 | 0 | +0.224 [+0.185, +0.264] | −0.011 [−0.019, −0.004] | +0.131 [+0.108, +0.156] |

Overall per tercile (C+G / K0 / K4): weak 0.226 / 0.419 / 0.403; mid 0.214 / 0.225 / 0.240; strong 0.220 / 0.220 / 0.278.

- **Config 1′'s gain over C+G is a weak-knockdown gain.**
  - **Weak tercile, +0.19:** MSE goes from 0 to ~1.0 there, because shrinking towards the control profile is right when little changes.
  - **Strong tercile, +0.000:** MSE +0.12 is cancelled by **NMAE −0.061**. Halving the called LFCs hurts where the true LFCs are large.
  - This is where Step 15's O1 oracle wanted a larger strong-tercile magnitude. A strength-aware magnitude remains the open lever, but (b) predicts strength too poorly to drive it (section 9).
- **K4's gain is a strong-knockdown gain:** +0.057 strong, +0.015 mid, −0.016 weak. That is the same split Step 17's proxy showed, and it is why it loses on HCT116, where almost every knockdown is weak.

## 2. Requested vs realised calls (realised = harness `de_pred` p_adj < 0.05 on the generated cells, own gene out)

| Line | Config | requested (median / total) | realised (median / total) | requested, not realised | leaked (realised, not requested) |
|---|---|---|---|---|---|
| Eval A | K0 | 486.5 / 161,002 | 482 / 159,127 | 2,547 (1.6%) | 672 |
| Eval A | K4 1,200 | 1,200 / 349,200 | 1,177 / 341,274 | 8,583 (2.5%) | 657 |
| HCT116final | K0 | 132.5 / 54,604 | 124 / 53,158 | 1,606 (2.9%) | 160 |
| HCT116final | K4 1,200 | 1,200 / 333,600 | 1,148 / 318,763 | 15,044 (4.5%) | 207 |

**By control-CPM decile of the target line's gated genes** (`<config>/<line>/calls_by_decile.csv`). Eval A, requested-not-realised / leaked, from decile 0 (lowest) to 9:
- **K0:** 13/470, 7/101, 38/22, 5/6, 22/1, 72/0, 242/2, 369/2, 591/0, 1,188/15. Ungated genes: 53 leaked.
- **K4:** 68/483, 82/84, 139/10, 269/3, 485/1, 857/0, 916/0, 1,189/1, 1,780/0, 2,798/11. Ungated genes: 64 leaked.

- **Failures to realise grow with expression.** Decile 9 accounts for 33–47% of them, and up to 5.6% of decile 9's calls on HCT116 K4.
- **Leaks sit almost all in the lowest decile and the ungated genes**, the same pattern as Step 9.
- **K4 loses 2.5–4.5% of its requested calls,** against 1.6–2.9% for K0.

## 3. The calls only K4 makes (K4 \ K0)

| Line | pairs (median per knockdown) | median \|requested Δ\| | median \|real Δ\| | true-DE share | gen / real cell SD, K4 run: median [IQR] | same genes in the K0 run (uncalled) |
|---|---|---|---|---|---|---|
| Eval A | 199,191 (713) | **0.025** | 0.030 | 13% | 1.03 **[0.54, 1.93]** | 1.00 [0.96, 1.04] |
| HCT116final | 279,996 (1,068) | **0.033** | 0.049 | 0.17% | 0.92 **[0.62, 1.24]** | 1.02 [0.97, 1.07] |

Reference rows on Eval A: C+G's own calls have a median |requested Δ| of 0.092, a true-DE share of 20%, and an SD ratio of 1.00 [0.88, 1.26] in the K4 run or 1.01 [0.83, 1.58] in the K0 run. Uncalled gated genes have 0.009, 9%, and 1.00 [0.97, 1.04].

- **SD:** cell-level standard deviation of log1p(CP10k) within each knockdown, generated against real, for up to 60 sampled genes per knockdown.
- **The K4-only calls request almost no shift** (|Δ| 0.025, a quarter of C+G's calls), yet the generator must make them significant. It does that by reshaping the cell-level spread at a fixed mean.
- **The cost is realism.**
  - Half of those genes end up more than 1.9× or less than 0.54× the real cell-level SD on Eval A.
  - Left uncalled (the K0 run), the same genes sit at 0.96–1.04, as realistic as any uncalled gene.
  - This is the DE-significance shaping that the rules brief rated the highest risk, made measurable.
- **On HCT116 only 0.17% of K4's extra calls are true DE.**

## 4. Proxy minus actual

| Line | Config | 15a proxy (fold cache) | proxy on the generator's input | actual (cells) | gap (input − actual) | largest component gap |
|---|---|---|---|---|---|---|
| Eval A | K0 | 0.2470 | 0.2456 | 0.2426 | +0.0030 | reach +0.018 |
| Eval A | K4 1,200 | 0.2711 | 0.2697 | 0.2638 | +0.0058 | reach +0.031 |
| HCT116final | K0 | — | 0.3679 | 0.3617 | +0.0062 | reach +0.042 |
| HCT116final | K4 1,200 | — | 0.3555 | 0.3476 | +0.0079 | reach +0.048 |

(HCT116final's 15a fold dropped HEK293T as well, so it has no comparable 15a row.)

- **The proxy tracks the cells to within 0.003–0.008, and all of the gap is reach.** Every other component is within 0.003.
- **The K4 − K0 difference:** proxy +0.0241, cells +0.0213.
- **The 0.0014 between 15a and the generator's input is a plumbing difference.**
  - The cause: genes with no source (0.9% of gated pairs, in every knockdown) reach the generator as "no change from the control's mean per-cell CPM" (`prepare_generator_inputs` zeroes their LFC). Config 1′ means them to be the control's pooled profile.
  - Keeping config 1′'s value lifts the proxy from 0.2456 to 0.2470 (PDS 0.795 → 0.804).
  - **That is a free +0.0014,** available in the model's `input.npz` (set `n_sources` ≥ 1 for those genes) without touching frozen code. v7's config 1 has the same issue; the freeze check compared unmasked means.
  - Not applied: it would change v7's frozen config.

## 5. Adoption rule (Eval A: MSE > 0, every other component > 0, overall above C+G with paired CI lower bound > −0.01; HCT116final: overall within its CI)

| Config | MSE > 0 | others > 0 | overall > C+G | CI lo > −0.01 | HCT116 within CI | **Result** |
|---|---|---|---|---|---|---|
| 1′ + K0 | 0.127 ✓ | ✓ (min NMAE +0.0077) | +0.024 ✓ | +0.016 ✓ | +0.171 [+0.127, +0.206] ✓ | **PASS** |
| 1′ + K4 1,200 | 0.127 ✓ | ✓ (min NMAE +0.0078) | +0.045 ✓ | +0.035 ✓ | +0.157 [+0.114, +0.191] ✓ | **PASS** |

- **How I applied the HCT116 test:** the paired CI against C+G must contain 0 or lie above it. Both lie well above.
- **Both pass on the point estimate for NMAE**, but NMAE's own CI crosses 0 ([−0.004, +0.019]).

## Recommendation (stopped before any D/E/F file)
**Both configs pass. The choice between them turns on what D/E/F look like, and on rules risk.**
- **Config 1′ + K4 1,200:**
  - Best on Eval A, +0.045 over C+G.
  - Worse than K0 on HCT116 (−0.014), and the costs land on weak knockdowns.
  - Its 713 extra calls per knockdown carry almost no requested shift, and their cell-level SD is off by up to 2× (IQR 0.54–1.93). That is the most visible form of DE-significance shaping.
- **Config 1′ + K0:**
  - +0.024 on Eval A, +0.171 on HCT116; better than K4 on the weak-responding line.
  - No extra shaping beyond C+G's own calls.
- **Expected value on D/E/F:** if they respond like H1, K4 gains about +0.02; if like HCT116 or K562, it costs about −0.014. Step 17 found that D/E/F's targets have source DE counts like the guardrail lines' (median 4 vs H1's 22), which leans weak, though it measures the targets, not the lines.
- **My recommendation: config 1′ + K0** for the D/E/F file, with the free +0.0014 no-source fix if you allow a change to the input step. Choose K4 1,200 only if you accept the shaping risk for an Eval A-type upside.
- **Your decision either way.**

## What I would have done differently
3. **Computed the proxy on the exact generator input from the start.** The masking of no-source genes was hidden by the freeze check, which compared means without that step. It's a small but free 0.0014.
4. **Split by strength in Step 16 already.** The tercile split shows config 1′'s gain is all weak-knockdown MSE, and K4's is all strong-knockdown calls. That would have pointed at strength-aware magnitudes earlier.

---

# 9. 2026-10-30 — Step 17: call allocation by predictor (b), volume-matched (proxy)
*Task: Jose, 2026-10-30, queued behind 15c. Noise-free proxy only, no cells. Code: `code/vcc6/step17.py`, `code/vcc6/step17.sh`. Tables and logs: `results/stage1/2026-10-30_step17/`.*

- **Checks:** the flat K4 1,200 rows reproduce 15a exactly (H1 halves 0.2901 / 0.2519; K562 0.140, HCT116 0.301, HEK293T 0.141).

**Setup:**
- **Magnitudes:** config 1′ for every row. Calls are the top N_k tested genes by C+G z.
- **(b)-allocation:** N_k = clip(λ × q(pct_k), 0, 4,000). pct_k is (b)'s percentile in the line's panel, and q is the true DE-count quantile function of the fit H1 half (item 1) or of all of H1 (item 2). Guardrail lines use H1's q.
- **Every row has a flat K4 partner** at the same **median** calls on the part it is scored on, as the task asks.
- **Halves:** Step 16's (seed 2026). MSE, PDS and NMAE don't read the calls, so within a half they are the same for every row: half 0 MSE 0.145, PDS 0.857, NMAE +0.0147; half 1 MSE 0.109, PDS 0.751, NMAE +0.0014.

## The volume caveat (read first)
An allocation rule moves calls from weak knockdowns to strong ones. At the same **median** it therefore makes many more calls in total. (b) at λ 4 has a median of 873 calls per knockdown but a **mean of 1,727**, against 873 for flat K4. The proxy rewards volume, so a median-matched comparison mixes allocation with volume.

So I added a side check, `meanmatch` → `tables/meanmatch.csv`: flat K4 at each row's **mean** calls. On H1 these are exact rows; on the guardrails they are interpolated linearly between flat rows. Both comparisons are reported below.

## 1. (b)-allocation, fitted on one H1 half, scored on the other

| λ | fit → scored | median / mean calls | (b) overall | flat, median-matched | Δ | flat, mean-matched | **Δ** | Δ fidelity / reach / Jaccard (median-matched) |
|---|---|---|---|---|---|---|---|---|
| 4 | 0 → 1 | 874 / 1,667 | 0.2561 | 0.2437 | +0.012 | 0.2615 | **−0.005** | +0.022 / +0.008 / +0.045 |
| 6 | 0 → 1 | 1,311 / 1,866 | 0.2619 | 0.2544 | +0.008 | 0.2642 | **−0.002** | +0.005 / +0.007 / +0.033 |
| 8 | 0 → 1 | 1,748 / 1,991 | 0.2652 | 0.2625 | +0.003 | 0.2661 | **−0.001** | −0.010 / +0.004 / +0.022 |
| 10 | 0 → 1 | 2,185 / 2,091 | 0.2674 | 0.2687 | −0.001 | 0.2674 | **0.000** | −0.022 / +0.001 / +0.013 |
| 14 | 0 → 1 | 3,059 / 2,239 | 0.2705 | 0.2759 | −0.005 | 0.2691 | **+0.001** | −0.031 / −0.001 / −0.001 |
| 4 | 1 → 0 | 872 / 1,787 | 0.3016 | 0.2823 | +0.019 | 0.3017 | **−0.000** | +0.033 / +0.014 / +0.069 |
| 6 | 1 → 0 | 1,308 / 1,951 | 0.3071 | 0.2925 | +0.015 | 0.3039 | **+0.003** | +0.019 / +0.014 / +0.054 |
| 8 | 1 → 0 | 1,744 / 2,064 | 0.3105 | 0.3011 | +0.009 | 0.3056 | **+0.005** | +0.007 / +0.008 / +0.042 |
| 10 | 1 → 0 | 2,180 / 2,157 | 0.3133 | 0.3076 | +0.006 | 0.3074 | **+0.006** | −0.000 / +0.002 / +0.032 |
| 14 | 1 → 0 | 3,052 / 2,303 | 0.3147 | 0.3143 | +0.000 | 0.3093 | **+0.005** | −0.019 / +0.004 / +0.017 |

Reference: flat K4 1,200 scores 0.2901 / 0.2519 (half 0 / half 1).

**Fitting λ on the fit half:**
- **By overall, the obvious rule:** it picks λ 14 on both halves, the top of the grid (5/5), because the overall rises with volume for every rule. Held out: 0.2705 / 0.3147. That is −0.005 / +0.000 against median-matched flat.
- **By gain over median-matched flat (the allocation question):** it picks λ 4 on both halves, the bottom of the grid (1/5). Held out: +0.012 / +0.019 against median-matched flat, but **−0.005 / −0.000 against mean-matched flat**.

**Guardrails** (fit half 0; fit half 1 agrees within 0.001). Δ is overall minus flat, at matched median / at matched mean (interpolated):

| λ | median / mean calls | K562 | Δ | HCT116 | Δ | HEK293T | Δ |
|---|---|---|---|---|---|---|---|
| 4 | 874 / 1,729 | 0.1404 | −0.001 / +0.003 | 0.3022 | −0.002 / +0.002 | 0.1415 | −0.001 / +0.002 |
| 6 | 1,311 / 1,912 | 0.1398 | −0.000 / +0.002 | 0.3034 | +0.003 / +0.003 | 0.1442 | +0.004 / +0.005 |
| 8 | 1,748 / 2,030 | 0.1398 | +0.002 / +0.002 | 0.3043 | +0.004 / +0.004 | 0.1421 | +0.003 / +0.003 |
| 14 | 3,059 / 2,274 | 0.1387 | +0.001 / +0.002 | 0.3023 | +0.002 / +0.002 | 0.1406 | +0.001 / +0.002 |

Flat K4 1,200 on the guardrails: 0.1404 / 0.3014 / 0.1410.
- **The guardrails barely care.** Every row is within ±0.004 of flat K4 1,200, at up to about twice the calls. Their true counts have a median of 1–3, so nearly every added call is false whichever rule places it, and Jaccard is already deeply negative.

## 1b. Oracle: N_k = clip(λ × true count, 0, 4,000)

| λ | median / mean calls | half 0 / half 1 | Δ vs flat, median-matched (h0 / h1) | Δ vs flat, mean-matched (h0 / h1) |
|---|---|---|---|---|
| 1 | 218 / 956 | 0.3294 / 0.2899 | **+0.078 / +0.077** | **+0.044 / +0.044** |
| 2 | 437 / 1,390 | 0.3368 / 0.2934 | +0.072 / +0.067 | +0.042 / +0.038 |
| 3 | 655 / 1,599 | 0.3357 / 0.2930 | +0.061 / +0.057 | +0.037 / +0.033 |
| 4 | 873 / 1,731 | 0.3334 / 0.2905 | +0.051 / +0.047 | +0.033 / +0.028 |
| 6 | 1,310 / 1,914 | 0.3310 / 0.2897 | +0.038 / +0.035 | +0.028 / +0.025 |
| 8 | 1,746 / 2,032 | 0.3303 / 0.2880 | +0.029 / +0.026 | +0.025 / +0.021 |
| 10 | 2,183 / 2,131 | 0.3284 / 0.2879 | +0.021 / +0.019 | +0.021 / +0.020 |
| 14 | 3,056 / 2,279 | 0.3285 / 0.2869 | +0.014 / +0.011 | +0.020 / +0.017 |

**Share of the oracle's gain that (b) recovers:**

| (b) λ | median-matched, oracle at the same λ (scored h1 / h0) | median-matched, oracle's best λ (λ 1: +0.077 / +0.078) | mean-matched, oracle at the same λ |
|---|---|---|---|
| 4 | 27% / 38% | 16% / 25% | −19% / −1% |
| 6 | 21% / 38% | 10% / 19% | −9% / +11% |
| 8 | 11% / 32% | 3% / 12% | −4% / +19% |
| 10 | −7% / 27% | −2% / 7% | 0% / +28% |
| 14 | −49% / 3% | −7% / 1% | +8% / +27% |

- **The oracle's allocation is real.** It survives mean matching (+0.02 to +0.044), and its best λ is 1, the bottom of the grid: exactly the true count.
- **(b)'s allocation mostly isn't.** Mean-matched, it recovers −19% to +28% of the oracle's gain, centred near zero. The best case (λ 10–14 on half 0) is about +0.006.

## 1c. Predictor (b) checks
**Spearman ρ with the true DE count** (each line's own tables excluded):
- **H1:** full 0.64 (n = 276 with a value); half 0 0.72; half 1 0.56.
- **K562:** 0.41 (262 with a value).
- **HCT116:** 0.50 (277).
- **HEK293T:** 0.41 (270).

**Cross-fit calibration on H1** (`tables/calibration_H1.csv`). Predicted counts come from the other half's quantiles. Deciles are of the predicted count; the values are medians, averaged over the two directions:

| decile | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| predicted | 2 | 4 | 8 | 43 | 147 | 327 | 829 | 1,681 | 2,596 | 4,674 |
| true | 6 | 8 | 25 | 147 | 411 | 708 | 404 | 1,747 | 1,866 | 1,969 |

- **The order is roughly right, but the shape is too steep.**
- **Deciles 1–6: too low,** by about 2–3.5× (decile 4: 43 vs 147).
- **Deciles 8–10: too high,** by up to 2.4× (4,674 vs 1,969). The truth flattens near 1,900 in the top three deciles.
- That is what a rank-to-quantile map does at ρ 0.64: it spreads predictions over the full true range, while the truth within a predicted decile regresses towards the middle.
- Decile 7 dips (404). The interquartile ranges within deciles are wide (for example 20–2,555 in decile 6).

## 2. Volume-free variant (q = H1's full true-count distribution, as D/E/F would borrow)

| λ | half | median / mean calls | overall | flat, median-matched | Δ | flat, mean-matched | Δ |
|---|---|---|---|---|---|---|---|
| 4 | 0 / 1 | 872 / 1,788, 1,666 | 0.3015 / 0.2560 | 0.2823 / 0.2436 | +0.019 / +0.012 | 0.3017 / 0.2615 | **−0.000 / −0.005** |
| 6 | 0 / 1 | 1,308 / 1,954, 1,864 | 0.3069 / 0.2617 | 0.2925 / 0.2543 | +0.014 / +0.007 | 0.3040 / 0.2642 | **+0.003 / −0.002** |
| 8 | 0 / 1 | 1,744 / 2,065, 1,988 | 0.3105 / 0.2651 | 0.3011 / 0.2624 | +0.009 / +0.003 | 0.3063 / 0.2660 | **+0.004 / −0.001** |

- **Within 0.001 of item 1** in every cell. The two halves' count distributions are balanced, so borrowing the whole of H1 changes nothing on H1.
- **The quantile shape is in-sample here** (H1's own counts), so this is not a held-out estimate of the borrowing step.
- **Guardrails** are the same as item 1 to 0.001: K562 0.1404 / 0.1399 / 0.1393, HCT116 0.3024 / 0.3034 / 0.3044, HEK293T 0.1414 / 0.1443 / 0.1422 at λ 4 / 6 / 8.

## 3. Per tercile of (b): (b)-allocation against flat K4, held out
Each knockdown is scored from the direction that holds it out. The flat partner is median-matched. CIs come from 1,000 bootstrap resamples within the tercile, with full-panel scaling. MSE, PDS and NMAE are identical between the arms.

**Gain-fitted λ 4 (the "best variant" by allocation gain; flat N 874 / 872):**

| (b) tercile (true DE median) | Jaccard Δ | fidelity Δ | reach Δ | overall Δ |
|---|---|---|---|---|
| weak (8) | +0.008 [−0.008, +0.027] | **−0.155 [−0.216, −0.099]** | +0.012 [+0.002, +0.023] | **−0.022 [−0.035, −0.011]** |
| mid (404) | +0.008 [−0.008, +0.026] | +0.001 [−0.033, +0.034] | +0.005 [−0.002, +0.013] | +0.002 [−0.006, +0.011] |
| strong (1,353) | **+0.153 [+0.111, +0.196]** | **+0.235 [+0.186, +0.288]** | +0.014 [+0.004, +0.027] | **+0.067 [+0.052, +0.082]** |

**Overall-fitted λ 14 (flat N 3,059 / 3,052):**

| (b) tercile | Jaccard Δ | fidelity Δ | reach Δ | overall Δ |
|---|---|---|---|---|
| weak | +0.009 [−0.005, +0.021] | −0.082 [−0.123, −0.043] | +0.001 [−0.014, +0.016] | −0.012 [−0.021, −0.004] |
| mid | −0.008 [−0.023, +0.004] | −0.020 [−0.037, −0.004] | +0.003 [−0.001, +0.009] | −0.004 [−0.009, +0.001] |
| strong | +0.024 [+0.015, +0.033] | +0.027 [+0.014, +0.041] | +0.001 [−0.004, +0.006] | +0.009 [+0.005, +0.012] |

- **What (b) does:**
  - It buys Jaccard and fidelity on the strong tercile, because those knockdowns get up to 4,000 calls instead of 873.
  - It pays in fidelity on the weak tercile. There it cuts calls to a few dozen, and fidelity, a yield over the true-DE genes, collapses.
  - The strong gain is mostly volume: that tercile gets about 3–4× the flat partner's calls, which is why the effect vanishes when the total is matched.

## 4. (b) for the 300 D/E/F targets (all 117 sources, none held out; `tables/b_DEF.csv`)
- **Coverage:** 298 of 300 targets have a value. The **2 without a source DE count are ANG and RAB5A**; ANG is the one target with no source at all.
- **Tables per target:** median 3, 90th percentile 5, maximum 50.

| | n | q10 | q25 | median | q75 | q90 | max | share 0 |
|---|---|---|---|---|---|---|---|---|
| D/E/F (all sources) | 298 | 0 | 1.0 | **4.0** | 23.0 | 105 | 1,777 | 11% |
| H1 (H1 tables out) | 276 | 0 | 1.0 | **21.9** | 113 | 324 | 3,109 | 17% |
| K562 / HCT116 / HEK293T (own out) | 262–277 | 0 | 0.5–1.0 | 3.4–4.5 | 19–31 | 97–123 | 1,238–1,813 | — |

- **The D/E/F targets look like the guardrail lines' targets, not H1's.** (b) is about 5× lower at the median and 3× lower at q90. Its values are source DE counts, which says something about the targets; it doesn't say how strongly D/E/F themselves respond.
- **The percentile map discards that level anyway.** D/E/F budgets would get H1's shape whatever their (b) values. On H1-like lines the shape is right; on guardrail-like lines (true median 1–3) it gives ~1,300 median / ~1,900 mean calls, which the proxy treats as neutral (section 1).

## Recommendation for one confirming cell run
**Flat K4, no allocation rule.**
- **Why not (b):** matched on total call volume, (b)-allocation is worth −0.005 to +0.006 per half on H1 (mean ≈ 0), and ±0.005 on the guardrails. Its apparent gain in the median-matched comparison is extra calls.
- **The parameter is N.**
  - N = 1,200 is Step 16's fit: chosen on both halves (spread 0), at the top of {400, 600, 800, 1,200} (grid position 4/4).
  - Step 16's side check shows the proxy still rising at 4,800.
  - 15c's run already confirms config 1′ + K4 1,200 on cells: Eval A 0.2638 against C+G's 0.2186, paired +0.045 [+0.035, +0.055]; full write-up in section 8.
- **So no new cell run is needed for allocation.**
  - If you want one more cell run, the question worth spending it on is volume: flat K4 at N ≈ 1,900, the total that (b) at λ 6 would have used.
  - On the proxy that adds +0.014 / +0.012 over 1,200 (mean-matched flat rows at 1,951 / 1,866: 0.3039 / 0.2642 vs 0.2901 / 0.2519). Guardrails are flat.
  - Going past 1,200 is your call: generator realism (see 15c item 3) and the DE-shaping rules risk both grow with N.
- **If an allocation rule is wanted anyway:**
  - The least-bad choice is the volume-free (b) map at **λ 6**: grid position 2/5. The allocation-gain fit picked λ 4 (1/5, an edge) on both halves, and the overall fit picked λ 14 (5/5, an edge).
  - Its held-out gain is +0.003 / −0.002 at matched mean volume, and +0.014 / +0.007 at matched median.
  - Expect nothing beyond what flat K4 at the same total gives.

## What I would have done differently
1. **Matched total (mean) call volume from the start, in Step 15's side check as well.** Median matching was the stated spec, but any rule that concentrates calls moves the mean, and the proxy rewards volume. I only caught it here, and it reverses Step 15's "+0.013 at K4's volume" (corrected in place in section 6).
2. **Stated the λ-fitting criterion in advance.** "Best overall on the fit half" just picks the largest volume (λ 14). The allocation question needs "best gain over volume-matched flat". I added the second fit after seeing the first.
3. **Put λ < 1 in the oracle grid.** The oracle's best λ (1) sits at the bottom edge.

---

# 10. 2026-10-30 — Step 18: D/E/F file for config 1′ + K0 (not uploaded)
*Task: Jose, 2026-10-30, after 15c (Jose's choice of config). Code: `code/vcc6/step18.py`, `step18_gen.sh`, `step18_final.sh`, `step18_prep.sh`; the model is `step15c.model`, refactored to take the lines and output folder, with 15c's behaviour unchanged. Outputs: `results/stage1/2026-10-30_step18/`.*

**The file:** `results/stage1/2026-10-30_step18/final/prediction_DEF.vcc`, 3,179,929,600 bytes (3.18 GB; md5 `ec8ac09a8c94545094c54790b0e9b6a0`). The merged h5ad it was built from is `final/prediction_DEF.h5ad` (15.56 GB). **Not uploaded.**

**Model:**
- Config 1′ + K0 with all 117 source tables, nothing held out. Called genes (C+G's calls) get 0.5 × C+G LFC. Non-called genes get 0.6 × c × Δ_pool, with c = z_a² / (z_a² + 1), single-source c₁ = 1, and no source → control.
- Each context D, E and F uses its own 18,400 controls, the 300 final targets and salt "". The generator's input step is unchanged, as in 15c.
- **Checks against the uploaded C+G file's inputs** (`results/stage1/2026-10-28_DEF/<x>/input.npz`), on all three contexts:
  - Calls and source counts are bit-identical.
  - C+G's LFC differs by 0.
  - Called LFC = 0.5 × C+G's exactly (max |diff| 0).

## Before generation (per context, from the saved arrays as the generator sees them)

| Context | Calls / knockdown, median (mean) | Uploaded C+G file | median \|LFC\| of called genes: 1′ + K0 / C+G | median \|LFC\|, uncalled: 1′ + K0 / C+G |
|---|---|---|---|---|
| D | 242.5 (352.9) | 242.5 (352.9) | 0.158 / 0.315 | 0.119 / 0.065 |
| E | 278.5 (380.8) | 278.5 (380.8) | 0.151 / 0.303 | 0.062 / 0.065 |
| F | 245.0 (355.9) | 245.0 (355.9) | 0.147 / 0.294 | 0.070 / 0.065 |

**Why D's uncalled shifts are larger:**
- Config 1′ targets the controls' pooled profile, the scorer's mean, rather than the mean per-cell CPM the generator starts from.
- D's controls vary more in library size, so that baseline offset is larger: median |log2(pooled / mean CPM)| on gated genes is 0.060, against 0.021 for E and 0.025 for F.
- Without the offset, the uncalled shifts are the same in all three contexts (0.012–0.013).

## After generation

| Context | Cells | Knockdowns × cells | Genes | Controls | Calls / knockdown, requested: median (mean) | Realised (generator's own test): median (mean) | recall / precision (median) |
|---|---|---|---|---|---|---|---|
| D | 120,000 | 300 × 400 | 18,533 | 0 | 242.5 (352.9) | 244.0 (354.6) | 1.000 / 0.995 |
| E | 120,000 | 300 × 400 | 18,533 | 0 | 278.5 (380.8) | 279.0 (383.4) | 1.000 / 0.995 |
| F | 120,000 | 300 × 400 | 18,533 | 0 | 245.0 (355.9) | 245.5 (358.6) | 1.000 / 0.995 |

- **Merge checks** (`tables/final.json`): shape 360,000 × 18,533, genes in the final `gene_names.csv` order, unique obs names, each context's 300 targets equal the final list, 400 cells per knockdown, no controls.

**`vcc prep`** (`logs/vcc_prep.log`; exit 0):
```
✓ prep complete
  cells:         360000
  genes:         18533  (encoding: float32)
  perturbations: column 'target_gene' → 'target_gene'
  contexts:      column 'context' → D=120000, E=120000, F=120000
  targets:       verified against the official list
  normalization: counts-preserved
```
- **Its warning:** 1,942,114,131 stored entries (5,395 per cell), which may need about 34 GiB of RAM.

## What I would have done differently

## Integrity check and upload (2026-10-30)
**Integrity check** (`code/vcc6/check18.py` → `results/stage1/2026-10-30_step18/tables/check.json`):
- **Full scans:** every file was read in full: the three generator outputs, the merged h5ad, the .vcc's own payload (decompressed to a temp copy, then deleted) and the uploaded C+G file. All counts are whole numbers; none is negative or missing.
- **Blocks:** 900 knockdown blocks of exactly 400 contiguous cells, and no controls.
- **Identity:** per-knockdown count sums are identical (max |diff| 0) from generator output to merged file, and from merged file to .vcc payload. The payload has the same block order and nnz (1,942,114,131).
- **Fidelity:** generated vs requested per knockdown, r ≥ 0.9999996 for the profiles and median 0.9996 for the changes (minimum 0.973, RAB5A, which requests almost no change). A one-gene shift drops r to about 0.005.
- **Against the uploaded C+G file:** the same layout, with changes 0.38–0.43× as large and a median correlation of 0.67–0.80, as config 1′ intends.
- **Recall 0:** only the one knockdown per context with zero calls requested.

**Upload (Jose's explicit OK, 2026-10-30):**
- `vcc submit results/stage1/2026-10-30_step18/final/prediction_DEF.vcc -m LHY`, 14:02:40–14:06:52 UTC. It uploaded 3.0 GiB in 4 min 01 s and exited 0.
- **Entry `ZGdlukwIpzkCZuC4M1Ix`**, FINAL mode (counts toward final scoring); scoring job `vcc-metrics-prod-zgdlukwipzkczuc4m1ix`. Log: `results/stage1/2026-10-30_step18/logs/vcc_submit.log`.
- As the last entry, it replaces C+G's `JFCghVNyHwfCOWhnHG6M` (10-28) for final scoring. Scores are hidden until the winners are announced.

**Stuck entry, cancelled and resubmitted (2026-10-30):**
- **Stuck:** entry `ZGdlukwIpzkCZuC4M1Ix` stayed at `launching` for 46 min (`past_window` true, `holds_slot` true, no error), and no confirmation email arrived.
- **Cancelled** on Jose's call: `vcc cancel --yes`, 14:49:45. The CLI reported "stopped the queued scoring job before it started running", so the job had been queued and never ran. The entry is now `failed` ("Cancelled by submitter"), which doesn't count against the daily limit.
- **Resubmitted** the same file (md5 `ec8ac09a…`, unchanged): 14:50:13–14:54:37, 3.0 GiB in 4 min 12 s, exit 0. **New entry `RTXvB7VtfATH5NQLobww`**, FINAL mode; scoring job `vcc-metrics-prod-rtxvb7vtfath5nqlobww`.
- Logs: `results/stage1/2026-10-30_step18/logs/vcc_cancel.log`, `vcc_submit_attempt2.log`, `vcc_status_watch*.log`.
- **15:20 UTC:** `RTXvB7VtfATH5NQLobww` is also still at `launching` with `past_window` true and `holds_slot` true, and no error. It is the same server-side stall as the first entry (job queued, never started), not the file: the 10-28 C+G file of the same size scored. Jose chose to leave it in place and stop polling. A support message to Arc, with both entry and job IDs, has been drafted for him to send.

---

# 11. 2026-11-01 — Step 19: studies before the next D/E/F file
*Task: Jose, 2026-11-01, with two supplements the same day (Part E replaced, adding the program regression E-reg; source-support gate rows added to Part D). No cells were generated, no new cells were scored, no D/E/F file was built and nothing was uploaded. Code: `code/vcc6/step19.py` (all parts and their unit tests), `step19.sh` (the chain), `step19_wait.sh` (waiter: chain log, tmux liveness, kernel OOM), `step19_queue2.sh` (supplement chain), `step19_fetch.sh` (copy from vcc2027). Tables, logs and the runtime estimate: `results/stage1/2026-11-01_step19/{tables,logs}`. Frozen code (generator, harness, `run_line.py`, `boot_score.py`, config 1′'s files, v7, Step 18's outputs) was imported only.*

## Part 0 — checks before any new result (`tables/part0.json`, `logs/part0.log`): all pass

**1. Combining the 15c runs.** The combined panel takes each knockdown's DE components (NMAE, fidelity, reach, Jaccard) from `results.csv` of the run that supplies it. It recomputes MSE and PDS on the combined pseudobulk through `mse_study.score_generated` (each knockdown's Ĉ_pred from its own run; the #348 budget formed on the combined panel) and `EvalPB.pds_rows`, and scales through `boot_score`'s path.

| Line | Endpoint | Rebuilt overall | 15c `summary.json` | Max \|diff\|, six components + overall | Per-knockdown max \|diff\| (MSE num / den / PDS) |
|---|---|---|---|---|---|
| Eval A | K0 everywhere | 0.242586 | 0.242586 | 5.2e-16 | 1e-16 / 1e-16 / 1e-16 |
| Eval A | K4 1,200 everywhere | 0.263839 | 0.263839 | 3.8e-16 | 1e-16 / 1e-16 / 1e-16 |
| HCT116final | K0 everywhere | 0.361650 | 0.361650 | 4.4e-16 | 1e-16 / 1e-16 / 1e-16 |
| HCT116final | K4 1,200 everywhere | 0.347623 | 0.347623 | 3.3e-16 | 1e-16 / 1e-16 / 1e-16 |

- **The two runs' predicted pseudobulk means are bit-identical** (max |P_K0 − P_K4| = 0 on both evals); only Ĉ_pred and the DE tables differ.

**2. Proxy.** Config 1′ + K0 on the generator's input reproduces 15c §4 exactly (0.245618; |diff| 0). On the fold cache it reproduces 15a: full 0.247034, halves 0.267867 / 0.226062 (|diff| ≤ 6e-17).

**Bootstrap scaling.** `boot_score.Run.scaled` costs 0.34 s a call, which would have made Part A's CIs take ~4.5 h. `step19.Fast` vectorises the same aggregate and maps it through the bundle's own `score_one` entries. Unit test: it matches `Run.scaled` to 4.4e-16 on 15 random resamples and a 20-knockdown slice, and reproduces `boot_score`'s resample sequence (seed 2026).

## Part A — hybrid calls from the existing 15c cells (`tables/partA_rows.csv`, `partA_paired.csv`, `partA.json`)
**Rule:** knockdown k takes K4 1,200's outputs if b_k ≥ T, else K0's. b is predictor (b) with the line's own tables out (H1: `step15b_predictors_H1.csv`, 276 of 292 with a value; HCT116: Step 17's `b_HCT116.csv`, 277 of 278). Knockdowns without a value take K0. "Calls" are the requested calls per knockdown (own gene out), from each supplying run.

**Eval A, full panel (scaled; MSE, PDS and NMAE move by ≤ 0.0002 across T, so only the call components are shown):**

| Rule | Share on K4 | Calls median / mean | **Overall** | Fidelity | Reach | Jaccard | Paired vs K0 (2,000 resamples) |
|---|---|---|---|---|---|---|---|
| K4 everywhere | 1.00 | 1,200 / 1,196 | **0.2638** | 0.483 | 0.048 | 0.123 | +0.0213 [+0.0145, +0.0276] |
| T 1 | 0.75 | 1,200 / 1,007 | 0.2623 | 0.460 | 0.051 | 0.133 | +0.0197 [+0.0144, +0.0249] |
| T 2 | 0.67 | 1,200 / 956 | 0.2621 | 0.460 | 0.050 | 0.134 | +0.0195 [+0.0147, +0.0246] |
| T 5 | 0.60 | 1,200 / 906 | 0.2627 | 0.458 | 0.055 | 0.134 | +0.0201 [+0.0155, +0.0251] |
| T 10 | 0.55 | 1,200 / 869 | 0.2610 | 0.450 | 0.055 | 0.132 | +0.0184 [+0.0141, +0.0231] |
| T 20 | 0.48 | 1,200 / 819 | 0.2583 | 0.437 | 0.056 | 0.127 | +0.0157 [+0.0116, +0.0200] |
| T 35 | 0.41 | 839 / 761 | 0.2553 | 0.427 | 0.057 | 0.119 | +0.0128 [+0.0087, +0.0169] |
| T 50 | 0.37 | 710 / 732 | 0.2537 | 0.420 | 0.057 | 0.116 | +0.0111 [+0.0073, +0.0150] |
| T 100 | 0.25 | 580 / 657 | 0.2493 | 0.402 | 0.057 | 0.108 | +0.0067 [+0.0036, +0.0100] |
| T 200 | 0.17 | 523 / 598 | 0.2451 | 0.384 | 0.058 | 0.099 | +0.0025 [+0.0005, +0.0046] |
| T 500 | 0.07 | 491 / 551 | 0.2426 | 0.373 | 0.060 | 0.093 | +0.0000 [−0.0012, +0.0011] |
| K0 everywhere | 0 | 487 / 551 | **0.2426** | 0.372 | 0.060 | 0.094 | — |

**HCT116final (full panel):**

| Rule | Share on K4 | Calls median / mean | **Overall** | Fidelity | Reach | Jaccard | Paired vs K0 |
|---|---|---|---|---|---|---|---|
| K4 everywhere | 1.00 | 1,200 / 1,200 | 0.3476 | 0.899 | 0.731 | −1.628 | −0.0140 [−0.0220, −0.0067] |
| T 1 | 0.76 | 1,200 / 968 | 0.3476 | 0.897 | 0.728 | −1.624 | −0.0141 [−0.0219, −0.0071] |
| T 2 | 0.68 | 1,200 / 882 | 0.3491 | 0.899 | 0.734 | −1.622 | −0.0125 [−0.0198, −0.0059] |
| T 5 | 0.49 | 778 / 684 | 0.3546 | 0.903 | 0.753 | −1.613 | −0.0070 [−0.0124, −0.0017] |
| T 10 | 0.38 | 271 / 577 | 0.3561 | 0.907 | 0.751 | −1.605 | −0.0056 [−0.0109, −0.0006] |
| T 20 | 0.29 | 204 / 478 | 0.3571 | 0.910 | 0.748 | −1.600 | −0.0046 [−0.0095, +0.0001] |
| T 35 | 0.23 | 173 / 424 | 0.3571 | 0.911 | 0.742 | −1.594 | −0.0046 [−0.0094, −0.0004] |
| T 50 | 0.18 | 160 / 372 | 0.3583 | 0.914 | 0.743 | −1.591 | −0.0034 [−0.0080, +0.0006] |
| T 100 | 0.12 | 146 / 311 | 0.3602 | 0.916 | 0.741 | −1.580 | −0.0014 [−0.0050, +0.0018] |
| T 200 | 0.06 | 138 / 255 | 0.3622 | 0.918 | 0.746 | −1.575 | +0.0006 [−0.0007, +0.0024] |
| T 500 | 0.03 | 135 / 224 | 0.3618 | 0.916 | 0.743 | −1.572 | +0.0001 [−0.0004, +0.0007] |
| K0 everywhere | 0 | 133 / 196 | 0.3617 | 0.918 | 0.738 | −1.570 | — |

Eval A halves (overall, half 0 / half 1; full-panel scaling): K4 everywhere 0.2812 / 0.2463; T 1 0.2828 / 0.2416; **T 2 0.2838** / 0.2403; T 5 0.2835 / 0.2418; T 10 0.2822 / 0.2397; T 20 0.2792 / 0.2373; T 50 0.2759 / 0.2314; T 100 0.2724 / 0.2261; T 200 0.2653 / 0.2247; K0 0.2628 / 0.2222.

**Fit (pick T on one half by overall; ties to the higher T; score the other half):**

| Fit half | T picked | Grid position (12 rules, K4-everywhere … K0-everywhere) | Held-out overall | K0 on that half | **Held-out gain** | Share on K4 (held out) | HCT116final at that T (paired vs K0) |
|---|---|---|---|---|---|---|---|
| 0 | **T 2** | 3/12 (interior) | 0.2403 (half 1) | 0.2222 | **+0.0181** (fidelity +0.082, Jaccard +0.040, reach −0.013) | 0.68 | 0.3491 (−0.0125 [−0.0198, −0.0059]) |
| 1 | **K4 everywhere** | 1/12 (edge) | 0.2812 (half 0) | 0.2628 | **+0.0184** (fidelity +0.116, Jaccard +0.016, reach −0.021) | 1.00 | 0.3476 (−0.0140 [−0.0220, −0.0067]) |

**At the picked rules, Eval A by true-strength tercile** (15c §1b terciles; paired CIs within the tercile, full-panel scaling):

| Tercile (true DE median) | T 2: share on K4 | T 2: overall vs K0 | K4 everywhere: overall vs K0 |
|---|---|---|---|
| weak (5) | 0.31 | −0.0100 [−0.0135, −0.0008] | −0.0162 [−0.0262, +0.0004] |
| mid (215) | 0.75 | +0.0081 [+0.0016, +0.0148] | +0.0148 [+0.0068, +0.0232] |
| strong (2,490) | 0.96 | +0.0533 [+0.0439, +0.0631] | +0.0574 [+0.0481, +0.0670] |

**Share of targets with b ≥ T:**

| T | 1 | 2 | 5 | 10 | 20 | 35 | 50 | 100 | 200 | 500 |
|---|---|---|---|---|---|---|---|---|---|---|
| D/E/F's 300 (all 117 sources; 298 with b) | 0.78 | 0.66 | 0.46 | 0.34 | 0.27 | 0.20 | 0.17 | 0.10 | 0.06 | 0.02 |
| H1 panel (292; 276 with b) | 0.75 | 0.67 | 0.60 | 0.55 | 0.48 | 0.41 | 0.37 | 0.25 | 0.17 | 0.07 |
| HCT116 panel (278; 277 with b) | 0.76 | 0.68 | 0.49 | 0.38 | 0.29 | 0.23 | 0.18 | 0.12 | 0.06 | 0.03 |

**Criterion** (Eval A held-out gain over K0 positive on both halves; HCT116final overall within K0's CI [0.303, 0.402]): **met** in both directions (+0.018 / +0.018; HCT116final 0.349 / 0.348).

**Reading:**
- **The criterion is met, but the hybrid adds nothing over plain K4 on Eval A.** One direction picks K4 everywhere, and the other picks T 2, which puts 67% of knockdowns on K4; full-panel overall rises almost monotonically as T falls.
- **On HCT116 the hybrid only stops costing at T ≥ 100,** where Eval A's gain is down to +0.007. At the picked rules HCT116 still pays −0.013 to −0.014 (paired CI below 0). The criterion passes only because K0's own CI on HCT116 is 0.10 wide.
- **(b) does sort knockdowns by strength:** at T 2 the weak tercile has 31% on K4 and the strong 96%. But dropping K4 on the weak tercile only halves its loss (−0.010 vs −0.016), and the Eval A total falls with the call volume, as Step 17 found for allocation.
- **D/E/F's (b) distribution is HCT116's, not H1's** (b ≥ 5: 46% / 49% / 60%). Part B shows that the guardrails' sparse truth is partly a cell-count effect D/E/F won't have.

## Part B — power and depth (`tables/partB.json`, `partB_kd_<eval>.csv`, `controls_<set>.npz`)
**1. Knockdowns:** real cells per knockdown (the eval's `real.h5ad`), true DE count (real p_adj < 0.05, own gene out, the eval's real DE table) and (b) (own line out; K562 from `step15b_predictors_K562.csv`).

| Eval | Knockdowns (with b) | Cells per knockdown: median [q10, q90] (min–max) | True DE median | Spearman true ~ cells | **Partial, given log1p(b)** | Spearman true ~ b | Spearman cells ~ b |
|---|---|---|---|---|---|---|---|
| Eval A (H1) | 292 (276) | **400** [352, 400] (138–400) | 218 | −0.35 | **−0.27** | 0.64 | −0.23 |
| HCT116final | 278 (277) | **144** [59, 310] (30–400) | 1 | +0.26 | **+0.23** | 0.50 | +0.12 |
| K562final | 267 (262) | **176** [80, 280] (35–400) | 3 | +0.12 | **+0.14** | 0.41 | −0.03 |
| HEK293Tfinal | 271 (270) | **188** [66, 389] (30–400) | 1 | +0.31 | **+0.30** | 0.41 | +0.10 |

Median true DE count by cell-count tercile (low / mid / high) within each (b) tercile:

| Eval | (b) weak | (b) mid | (b) strong |
|---|---|---|---|
| Eval A | 13.5 / 8 / 7 | 984 / 381 / 306 | 2,210 / 1,653 / 1,098 |
| HCT116final | 0 / 0 / 1 | 1 / 1 / 4 | 3 / 4 / 8 |
| K562final | 0 / 2 / 1 | 2 / 3 / 6 | 11 / 25 / 6 |
| HEK293Tfinal | 0 / 1 / 0 | 0 / 2 / 2 | 1 / 3 / 6 |

- **Eval A is capped at 400 cells, and strong knockdowns lose cells** (partial ρ −0.27): depletion, not power.
- **On the three guardrail lines, more cells means more true DE genes at the same (b)** (partial +0.14 to +0.30; within the strong-(b) tercile, the median doubles from the lowest to the highest cell tercile on HCT116 and HEK293T). Their sparse truth is partly a power limit, at a median of 144–188 cells.

**2. Control cells** (all 18,400 or all available; UMI = per-cell total over the 18,533 panel genes):

| Set | Cells | UMI q5 / q25 / **q50** / q75 / q95 | UMI CV | Genes > 5 CPM | **Baseline offset** |
|---|---|---|---|---|---|
| D | 18,400 | 6,436 / 13,960 / **20,236** / 27,387 / 39,868 | 0.47 | 9,633 | **0.060** |
| E | 18,400 | 7,282 / 13,869 / **19,681** / 26,887 / 38,715 | 0.46 | 9,743 | **0.021** |
| F | 18,400 | 7,029 / 13,647 / **19,845** / 27,047 / 38,970 | 0.46 | 9,537 | **0.025** |
| Eval A context | 18,400 | 10,794 / 15,966 / **20,001** / 25,098 / 34,602 | 0.35 | 10,778 | 0.010 |
| Eval A scoring | 18,400 | 11,015 / 15,951 / **19,989** / 25,050 / 34,595 | 0.35 | 10,778 | 0.010 |
| HCT116final context / scoring | 9,940 / 9,940 | 5,669 / 8,797 / **12,133** / 16,909 / 26,846 (context) | 0.51 | 10,326 / 10,320 | 0.021 / 0.021 |
| K562final context / scoring | 18,400 / 18,400 | 3,174 / 5,804 / **7,915** / 10,055 / 14,158 (context) | 0.42 | 7,681 / 7,681 | 0.015 / 0.015 |
| HEK293Tfinal context / scoring | 9,942 / 9,943 | 6,868 / 11,021 / **15,171** / 20,880 / 35,592 (context) | 0.56 | 11,066 / 11,070 | 0.038 / 0.040 |

**Check:** D 0.0602, E 0.0213, F 0.0249 reproduce Step 18's 0.060 / 0.021 / 0.025 (pass).

**3. Final-phase truth, as found.** `manifest.json`: `cells_per_pert` 400, `n_constructs` 300; per context D/E/F: `n_perturbations` 300, `control_cells` 18,400, `n_ntc_ids` 46, `ground_truth_cells` 138,400. `pert_counts.csv` has one column (`target_gene`, 300 rows) and no counts. 138,400 = 18,400 + 300 × 400, so if `ground_truth_cells` includes the controls (my reading), the truth has **400 cells per knockdown**, like Eval A's cap.

**Where D, E and F sit:**
- **Cells per knockdown:** at Eval A's level (400), against 144–188 on the guardrails.
- **Depth:** median UMI 19.7–20.2k, the same as Eval A's 20.0k. That is 1.3–2.5× the guardrails' 7.9–15.2k.
- **Spread:** library-size CV 0.46–0.47 is wider than Eval A's 0.35, like HCT116 / K562 / HEK293T.
- **Offset:** E and F (0.021–0.025) are at HCT116's level. D (0.060) is the largest of any set, 6× Eval A's.
- **So by power, D/E/F's truth should be denser than the guardrails' at the same knockdown strength.** Their sparse truth is not a prior for D/E/F on depth grounds. Target strength ((b), Part A) still looks like the guardrails'.

## Part C — fingerprint of D/E/F (`tables/partC_markers.csv`, `partC_spearman.csv`, `partC_marker_similarity.csv`, `partC.json`)
**Controls on disk:** A, B and C (`data/context_{A,B,C}.h5ad`, 18,400 each), D/E/F, and the Arc H1 controls (Eval A context, plus `arc_h1_baseline.csv`). XIST and RPS4Y1 are not on the 18,533-gene panel.

**Mean CPM (detection fraction) of the deciding markers** (full table of all 52 markers × all profiles in `partC_markers.csv`):

| Marker | D | E | F | Arc H1 | A | B | C |
|---|---|---|---|---|---|---|---|
| POU5F1 | 4.8 (0.09) | 0.3 | 1.8 | **3,087 (1.00)** | 0.1 | 2.4 | 6.0 |
| NANOG / L1TD1 / DPPA4 | 0.4 / 0 / 0 | 0 | 0 | 168 / 1,514 / 649 | 0 | 0 | 0 |
| SPI1 (PU.1) | **886 (0.995)** | 0 | 0 | 0.7 | 0 | 0.2 | 0.1 |
| PTPRC (CD45) | **27 (0.38)** | 0 | 0 | 0 | 21 | 0 | 0 |
| CD3E / CD3D / LCK / TAL1 | 0 / 2 / 0 / 0 | 0 | 0 | 0 / 0 / 73 / 0 | **760 / 939 / 431 / 149** | 0 | 0 |
| GATA1 / HBG2 | 1 / 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| VIM | 3,133 | **31,826** | **27,498** | 351 | 1,612 | 6,052 | 2 |
| COL1A1 / COL1A2 | 1 / 0 | 317 / 405 | **1,973 / 687** | 196 / 132 | 4 / 0 | 358 / 0 | 10 / 0 |
| THY1 | 0 | 0.4 | **542 (0.98)** | 477 | 78 | 2 | 0 |
| NES / PAX6 | 0 / 0 | **51 (0.52) / 9.5 (0.16)** | 4 / 0 | 503 / 0 | 1 / 11 | 39 / 34 | 1 / 0 |
| KRT8 / KRT18 / KRT19 | 13 / 14 / 0 | 21 / 5 / 0 | 67 / **257** / 19 | 533 / 749 / 3 | 4 / 0 / 0 | **2,506 / 835 / 140** | 1,427 / 667 / 4 |
| EPCAM / CDH1 | 0 / 0 | 2 / 0 | 1 / 0 | 402 / 380 | 0 | 2 / 0 | 310 / 127 |
| DDX3Y / KDM5D / UTY | **157 / 72 / 72** | 0 | 0 | 122 / 39 / 15 | 12 / 51 / 33 | 0 | 95 / 40 / 40 |
| CDKN2A | 0 | 0 | **941 (0.998)** | 0.4 | 0 | 2,480 | 54 |
| TP53 / MDM2 / CDKN1A | 372 / 253 / 260 | 85 / 180 / 122 | 519 / 138 / 311 | 461 / 184 / 6 | 22 / 140 / 2 | 418 / 130 / 140 | 238 / 181 / 32 |
| MKI67 / TOP2A | 217 / 203 | 248 / 272 | 522 / 475 | 175 / 331 | 408 / 607 | 362 / 398 | 274 / 277 |

**Working identities (my reading):**
- **D: a myeloid (monocyte/macrophage-like) leukaemia line, male, CDKN2A silent.**
  - Evidence: SPI1 at 886 CPM in 99.5% of cells and PTPRC+.
  - It is not erythroid like K562 (no GATA1 or HBG) and not T-lineage (no CD3).
  - Y genes are expressed, TP53 is expressed, and no epithelial or pluripotency genes are present.
- **E: a mesenchymal line, Y-negative, CDKN2A silent, TP53 low.**
  - Evidence: VIM 31.8k CPM and collagens.
  - Nestin is detected in 52% of cells and PAX6 in 16%, so it could be glial or neural-crest-like.
  - No keratins.
- **F: a fibroblast-like mesenchymal line with simple-epithelial keratins, Y-negative, CDKN2A high** (941 CPM in 99.8% of cells, which usually means the RB pathway is off).
  - Evidence: VIM 27.5k CPM, COL1A1 2.0k, THY1+, KRT18+.
- **None is pluripotent and none resembles H1:** POU5F1 ≤ 5 CPM against H1's 3,087, and NANOG/L1TD1/DPPA4 are absent.
  - A looks like Jurkat (CD3/LCK/TAL1, male); C is epithelial and male.

**Correlations.** Whole-transcriptome Spearman of log1p mean CPM over each context's gated genes ranks every Arc context first (D: A 0.61 … E 0.54; E: F 0.81, B 0.76; F: E 0.81, B 0.80), then Arc H1 (0.57–0.61), then every non-Arc source below 0.47. **This mostly measures platform (Arc Flex against the rest),** as earlier work found.

**Marker-only similarity** (Pearson of log1p mean CPM over the 45 panel markers; my addition, to read biology rather than platform):
- **F:** Tahoe HS-578T 0.92, RPMI-7951 0.92, B 0.84, HOP62 0.80, LOX-IMVI 0.79, SW 1088 0.79.
- **E:** Tahoe RPMI-7951 0.84, HS-578T 0.81, MIA PaCa-2 0.81, HOP62 0.80, hTERT-HPNE 0.80.
- **D:** best 0.72 (Tahoe A549 / A-427 / H4). No source is myeloid; the highest SPI1 in any source baseline is K562's 49 CPM.

**Biologically closest source lines:**
- **E and F:** Tahoe's mesenchymal lines (HS-578T, RPMI-7951, SW 1088, HOP62). These enter the model only through drug rows, not knockdowns.
- **D:** nothing close. By lineage, the haematopoietic suspension lines K562, Jurkat and CD4 T. CD4 T (Marson) has no control CPM in its baselines.
- **Among the genetic screens,** no line is mesenchymal or myeloid.

No model change was made.

## Part D — called-gene scale by expression tier (proxy)
**Base:** the uploaded file's config on Eval A: config 1′ + K0 on the generator's input (15c `config1p_K0/EvalA/input.npz`, no-source genes zeroed), 0.2456. The guardrails use the fold pools (fold caches, held out as before), config 1′ + K0 as in 15a: K562 0.1582, HCT116 0.3108, HEK293T 0.1536. Every line's base row reproduces its reference to ≤ 2e-16.

### D1 — where the uploaded config's error sits (Eval A; `tables/partD1.json`)
Deciles are of H1's context-control mean CPM over its 10,778 gated genes (decile −1 = ungated). "Num share" is each category's mean over knockdowns of Σ_cat |pred − real| / Σ_gate |real|; these sum to the raw NMAE (0.9461, 216 knockdowns, which reproduces the proxy exactly). "Den share" is the category's share of Σ |real|. A ratio num/den near 1 means the prediction there is about zero.

| Category | NMAE num share | Den share | num / den | Gate (k, g) pairs |
|---|---|---|---|---|
| called | 0.0995 | 0.1230 | **0.81** | 32,380 |
| uncalled | **0.8466** | 0.8770 | 0.97 | 273,474 |
| pred \|LFC\| < 0.1 | 0.7406 | 0.7566 | 0.98 | 250,538 |
| 0.1–0.25 | 0.1447 | 0.1666 | 0.87 | 45,253 |
| 0.25–0.5 | 0.0448 | 0.0576 | 0.78 | 8,837 |
| 0.5–1 | 0.0132 | 0.0178 | 0.74 | 1,126 |
| > 1 | 0.0028 | 0.0014 | **2.0** | 100 |
| decile 0 / 4 / 9, called | 0.0005 / 0.0058 / 0.0385 | 0.0007 / 0.0075 / 0.0451 | 0.71 / 0.77 / 0.85 | 114 / 1,412 / 14,560 |
| decile 0 / 4 / 9, uncalled | 0.1119 / 0.0788 / 0.0940 | 0.1157 / 0.0817 / 0.0967 | 0.97 / 0.96 / 0.97 | 18,059 / 24,595 / 44,328 |

Sign accuracy on true DE genes, called / uncalled, by decile 0 → 9:
- **Called:** 0.83, 0.88, 0.82, 0.80, 0.79, 0.79, 0.78, 0.77, 0.78, 0.75.
- **Uncalled:** 0.60–0.62 throughout.

MSE plug-in error (scorer space, own gene out, raw-MSE units; total 1.390):
- **Called genes:** 15.4% of it, rising with expression from 0.07% (decile 0) to 4.4% (decile 9).
- **Uncalled gated genes:** 75.7%, spread evenly (5.8–8.9% per decile).
- **Ungated genes:** 8.9%.

**Reading:**
- **The called genes are the part of the prediction that works.** They carry 12% of NMAE's |real| mass, at error ratio 0.81 and sign accuracy 0.75–0.88; their share of MSE grows with expression.
- **85% of NMAE's numerator sits on uncalled genes** whose prediction is essentially zero. The called-gene scale can only move the remaining ~10% (and ~15% of the MSE plug-in).
- **The few called genes with predicted |LFC| > 1 overshoot**, at twice the real magnitude, but they are 100 pairs.

### D2 — tier grid (`tables/partD_rows_<line>.parquet`, `partD.json`)
**Grid:**
- **Tiers:** terciles of each line's own context-control mean CPM among its gated (> 5 CPM) genes; H1: 3,592 / 3,593 / 3,593 genes, and no called gene is ungated.
- **Scales:** called gated genes get s_tier × C+G LFC, with each s ∈ {0.35, 0.5, 0.65, 0.8, 1.0}.
- **Cap:** after scaling, |LFC| is capped at {none, 2, 1}.
- **Unchanged:** uncalled genes stay config 1′'s values, and calls stay K0.

That is 375 rows × (H1 full / half 0 / half 1, K562, HCT116, HEK293T). It ran as four single-threaded workers, 15:44–17:35 (H1 was the slowest line, at 19–21 s a row).

**Fit (pick on one half among rows with all six components > 0 on that half; score the other):**

| Fit half | Eligible rows | Pick (s_low / s_mid / s_high, cap) | Grid positions | Held-out overall vs base | Gain: overall / MSE / PDS / NMAE / fidelity / reach / Jaccard | Guardrails at the pick (K562 / HCT116 / HEK293T, vs base) |
|---|---|---|---|---|---|---|
| 0 | 300 | 0.65 / 0.65 / 0.5, cap 2 | 3/5, 3/5, 2/5 (interior); cap 2/3 (interior) | 0.2263 vs 0.2247 (half 1) | **+0.0017** / +0.007 / +0.002 / +0.001 / 0.000 / −0.000 / 0 | +0.0035 / +0.0049 / +0.0036 |
| 1 | 81 | 0.8 / 0.65 / 0.5, cap 1 | 4/5, 3/5, 2/5 (interior); cap 3/3 (**edge**) | 0.2680 vs 0.2664 (half 0) | **+0.0016** / +0.011 / −0.003 / +0.002 / 0 / −0.001 / 0 | +0.0040 / +0.0047 / +0.0034 |

Only 81 rows qualify on half 1, because there NMAE is ≤ 0 for most rows that raise s_high.

**In-sample H1 surface (full panel; base 0.2456):**
- **Best row per cap:** none and 2 tie at (0.65, 0.65, 0.5) with 0.2474; cap 1 gives (0.8, 0.65, 0.5) with 0.2473. **A cap at 2 changes nothing,** and a cap at 1 changes the fourth decimal.
- **Marginal mean overall over all rows with that value (no cap):**
  - s_low 0.35 → 1.0: 0.2346, 0.2349, **0.2350**, 0.2348, 0.2341 (flat);
  - s_mid: 0.2370, **0.2373**, 0.2367, 0.2343, 0.2283;
  - s_high: 0.2416, **0.2430**, 0.2393, 0.2304, 0.2193.
- **Marginal mean NMAE:** s_high 0.5 is the only value with clearly positive NMAE (+0.0072; s_high 1.0 −0.0149).
- **So the high-expression tier wants exactly config 1′'s 0.5,** and the low and mid tiers want slightly more (0.65), worth ~+0.002.

**Criterion** (held-out gain ≥ +0.005 on both halves, NMAE > 0 on both, no guardrail more than 0.005 worse): **not met.** The gains are +0.0017 / +0.0016; NMAE stays positive and the guardrails gain +0.003 to +0.005.

Reach is unchanged by the tier rows (Δ ≤ 0.001), and the proxy overstates reach by 0.02–0.05 compared with cells (15c §4).

### Addendum — source-support gate rows (`tables/partDgate_rows.csv`, `partDgate.json`; run separately, 15:54–16:04 UTC, alongside the D grid)
**Support value** for a called gene in knockdown k: the median, over the source backgrounds that hold k (the line's fold lines excluded), of the background's control mean CPM (`<source>_baseline.csv`, its tables averaged first).
- **"Holds" means a genetic row for k.** Tahoe's drug rows do carry a `target_gene` (the drug's annotated target, e.g. AZD-8055 → MTOR), but they are not knockdowns of it; I left them out, as the Part E supplement also does.
- **Marson's baselines carry no control CPM,** so CD4 T never contributes a value.
- A called gene with no value (no holding background measures it) counts as below every threshold: 0 on H1, 51 on K562, 2,156 on HCT116, 2,384 on HEK293T.

**1. How much each threshold touches** (share of called (k, g) pairs; share of the called genes' total |LFC| at the base):

| Line | Called pairs | < 5 CPM: pairs / \|LFC\| | < 20 CPM: pairs / \|LFC\| |
|---|---|---|---|
| H1 | 161,002 | 5.5% / 8.5% | 18.1% / 25.1% |
| K562 | 67,699 | 11.3% / 18.5% | 22.4% / 34.5% |
| HCT116 | 154,463 | 3.2% / 4.4% | 11.3% / 15.4% |
| HEK293T | 186,890 | 3.5% / 4.9% | 12.2% / 16.9% |

The touched genes do carry larger called LFCs than average, which is the premise.

**2. Where they sit in D1 (Eval A):**

| Set | NMAE num share | Den share | num / den | Sign accuracy on true DE | Share in H1 deciles 0–2 / 3–6 / 7–9 |
|---|---|---|---|---|---|
| all called | 0.0995 | 0.1230 | 0.81 | 0.769 | 3% / 28% / 70% |
| support < 5 | 0.0125 | 0.0145 | 0.86 | 0.755 | 16% / 41% / 43% |
| support < 20 | 0.0279 | 0.0331 | 0.84 | 0.731 | 11% / 49% / 40% |

- **The low-support genes are only slightly worse than other called genes:** error ratio 0.84–0.86 against 0.81, sign accuracy 0.73–0.76 against 0.77.
- **H1 expresses them well:** 40–43% sit in its top three deciles. They are genes the sources barely express but the target line does; the LFC is not inflated noise on a gene H1 lacks.

**3. Scores** (scaled; Δ against the base on the same line and part):

| Line / part | Base | G1 (< 5, × 0.5) | G2 (< 5, × 0) | G3 (< 20, × 0.5) | G4 (< 20, × 0) |
|---|---|---|---|---|---|
| H1 full | 0.2456 | −0.0002 | −0.0011 | −0.0012 | −0.0057 |
| H1 half 0 | 0.2664 | −0.0009 | −0.0020 | −0.0015 | −0.0072 |
| H1 half 1 | 0.2247 | +0.0006 | −0.0001 | −0.0009 | −0.0041 |
| H1 NMAE, half 0 / half 1 | 0.0138 / 0.0015 | 0.0121 / 0.0013 | 0.0094 / **−0.0013** | 0.0108 / **−0.0008** | 0.0046 / **−0.0069** |
| K562 | 0.1582 | +0.0031 | −0.0003 | +0.0008 | −0.0082 |
| HCT116 | 0.3108 | −0.0030 | −0.0057 | −0.0057 | −0.0184 |
| HEK293T | 0.1536 | +0.0002 | +0.0056 | −0.0020 | −0.0066 |

On H1 every gate row lowers MSE (−0.005 to −0.022) and NMAE and raises reach a little (+0.003 to +0.005), and the proxy overstates reach by 0.02–0.05 compared with cells (15c §4). Full six-component rows are in the CSV.

**Criterion** (gain ≥ +0.005 on both halves, NMAE > 0 on both, no guardrail more than 0.005 worse): **no row meets it.** G1 is the mildest: −0.0009 / +0.0006, guardrails −0.003 to +0.003. G2–G4 turn half 1's NMAE negative, and G2–G4 cost HCT116 0.006–0.018.

**Reading:** the hypothesis that large called LFCs on source-unexpressed genes are inflated is not borne out on H1. Halving or zeroing them removes real signal (MSE falls), so I would not use any gate row.

## Part E — programs (PCA on source knockdowns), the ceiling, then the regression (`tables/partE.json`, `partE_ceiling.csv`, `partE_<basis>_{basis.npz,rows.parquet}`, `partE_DEF_target_scores.csv`; `logs/partE.log`)
**Which text ran.** Jose's supplement replaced Part E before Part E had started. Nothing was built under the earlier text, so nothing was kept or rebuilt: every basis and ceiling below is from the supplement's spec.

**Build.** Rows are every perturbation row of every source table except Arc H1, with Tahoe's drug rows and Marson's lfc/se included. Entries are LFC × z²/(z² + 1), 0 at ≤ 5 CPM in the table's own controls and where non-finite. Each table is centred per gene over its measured entries. Rows are weighted 1/√n_b per background, with all Tahoe tables as one background. Tables are loaded one at a time into a preallocated float32 matrix (E1: 97,763 × 18,533, 7.2 GB; process peak 15.5 GB). Randomized SVD, 100 components (20 oversamples, 5 power iterations, seed 2026). E1−K562, E1−Orion and E5 drop whole backgrounds, so each is exactly E1 built without those tables (same per-table centring, same weights for the others); they were cut from E1's matrix.

| Basis | Rows | Tables | Largest backgrounds (rows) | Cum. variance K 5 / 10 / 20 / 50 / 100 | K at 25% |
|---|---|---|---|---|---|
| **E1** | 97,763 | 116 | CD4 T 32,040; HCT116 16,492; HEK293T 15,938; K562 9,567; Tahoe 6,808; Jurkat 2,054; HipSci 2,047; RPE1 1,968 | 0.055 / 0.078 / 0.112 / 0.157 / **0.196** | **not reached by 100** (K 100 used) |
| E1−K562 (fold basis) | 88,196 | 107 | as E1 without K562 | 0.056 / 0.079 / 0.113 / 0.158 / 0.197 | not reached (K 100) |
| E1−Orion (fold basis) | 65,333 | 112 | as E1 without HCT116, HEK293T | 0.057 / 0.079 / 0.114 / 0.159 / 0.199 | not reached (K 100) |
| E2 (n_real ≥ 50) | 13,418 | 107 | Tahoe 6,588; K562 1,521; RPE1 1,383; KOLF 1,222; Jurkat 993; HepG2 837 | 0.105 / 0.146 / 0.193 / 0.279 / 0.366 | **39** |
| E3 (20 CPM gate) | 97,763 | 116 | as E1 | 0.079 / 0.110 / 0.148 / 0.207 / 0.250 | 100 (just) |
| E4 (Δ space) | 65,723 | 113 | as E1 without Marson (no full profiles) | 0.074 / 0.106 / 0.148 / 0.211 / 0.253 | **96** |
| E5 (E1 without Tahoe) | 90,955 | 66 | as E1 without Tahoe | 0.056 / 0.078 / 0.113 / 0.155 / 0.196 | not reached (K 100) |

Per-background row counts for every basis are in `partE.json`. Marson (CD4 T) contributes no rows to E2: it has no `n_real`, and neither do the Orion pseudobulk tables.

**Tahoe's share of E1's top-10 components' row-score energy:** PC1 0.04, PC2 0.005, PC3 0.02, **PC4 0.61**, PC5 0.001, **PC6 0.48**, PC7 0.14, PC8 0.07, PC9 0.06, PC10 0.01.

### Ceiling on H1 (Eval A, 292 knockdowns; captured energy = 1 − Σ‖x − x̂‖² / Σ‖x‖² over the group)
(a) LFC, all gated (tested) genes; (b) LFC, the knockdown's true-DE genes with (a)'s coefficients; (c) Δ space (LFC fit converted gene by gene with H1's control CP50K; E4 fits Δ directly).

| Basis, K | (a) all / weak / mid / **strong** | (b) all / weak / mid / **strong** | (c) all / weak / mid / **strong** |
|---|---|---|---|
| E1, 5 | 0.019 / 0.002 / 0.007 / **0.024** | 0.030 / 0.002 / 0.012 / **0.032** | 0.022 / 0.000 / 0.009 / **0.027** |
| E1, 10 | 0.039 / 0.005 / 0.015 / **0.049** | 0.061 / 0.004 / 0.024 / **0.064** | 0.054 / 0.004 / 0.023 / **0.065** |
| E1, 20 | 0.057 / 0.007 / 0.023 / **0.072** | 0.092 / 0.006 / 0.039 / **0.097** | 0.076 / 0.005 / 0.032 / **0.091** |
| E1, 50 | 0.170 / 0.018 / 0.088 / 0.212 | 0.274 / 0.017 / 0.194 / 0.283 | 0.221 / 0.018 / 0.129 / 0.259 |
| E1, 100 (K "25%") | 0.197 / 0.026 / 0.106 / 0.243 | 0.313 / 0.029 / 0.229 / 0.323 | 0.252 / 0.021 / 0.149 / 0.295 |
| E5, 20 / 100 | 0.058 / 0.194 | 0.092 / 0.308 | 0.079 / 0.248 |
| E2, 20 / 39 | 0.117 / 0.147 (strong 0.149 / 0.185) | 0.195 / 0.242 (strong 0.204 / 0.252) | 0.157 / 0.194 (strong 0.188 / 0.231) |
| E3, 20 / 100 | 0.055 / 0.152 | 0.093 / 0.243 | 0.090 / 0.240 |
| E4 (Δ), 20 / 96 | — | 0.177 / 0.453 | 0.116 / 0.319 (strong 0.137 / 0.369) |
| Random orthonormal, 20 / 50 (mean of 5) | 0.002 / 0.005 | 0.003 / 0.007 | 0.001 / 0.001 |
| H1's own top-K PCs (in sample), 5 / 10 / 20 / 50 | 0.30 / 0.42 / 0.53 / 0.69 (strong 0.38 / 0.53 / 0.68 / 0.86) | 0.46 / 0.62 / 0.76 / 0.90 | 0.34 / 0.47 / 0.59 / 0.74 |
| S2 (cited; source PCs of Δ, fold-excluded), 5 / 10 / 20 / 50 | — | — | 0.043 / 0.087 / 0.125 / 0.213 (all) |
| **Config 1′'s actual Eval A prediction** (15c generator input) | **−0.109** / −0.355 / −0.342 / −0.019 | **0.042** / 0.116 / −0.089 / **0.050** | **−0.040** / −0.610 / −0.290 / **0.066** |

**Ceiling criterion** (E1, strong tercile, K ≤ 20: clearly above config 1′ on (b) and (c), and well above S2's 4–21%): **not met.**
- **(b):** 0.032 / 0.064 / 0.097 at K 5 / 10 / 20, against config 1′'s 0.050. Above it only from K 10, and about 2× at K 20.
- **(c):** 0.027 / 0.065 / 0.091, against 0.066. Level at K 10, +0.025 at K 20.
- **Against S2:** at the bottom of S2's range, not well above it. On all knockdowns at K 20, (c) is 0.076 against S2's 0.125.

**Reading:**
- **The source programs are real but thin.** They are 30–40× a random basis, yet E1 holds under 10% of H1's strong-knockdown energy at K 20 and ~25–30% even at K 100, against 68% (K 20) for H1's own PCs.
- **Restricting to strong source rows (E2) doubles the ceiling at the same K** (strong (b) 0.20 at K 20). That is the best of the LFC bases.
- **Δ space (E4) captures more of the DE genes' energy** (0.18 at K 20, 0.45 at K 96).
- **Tahoe adds almost nothing to the ceiling** (E1 vs E5 within 0.005).
- **Config 1′'s own prediction is worse than zero on all gated genes** in LFC (−0.11) and Δ (−0.04). Only on the true-DE genes of strong knockdowns is it positive (0.05 / 0.07), which matches D1: the uncalled background, not the calls, costs it.

### Programs (E1's top 10; labels are **my reading**)

| PC | Var. share | + side (top genes) | − side (top genes) | Largest row scores | Label (my reading) |
|---|---|---|---|---|---|
| 1 | 1.7% | CDKN1A, TP53INP1, DDIT3, BBC3, ZMAT3, MDM2, RRM2B, SESN2, PHLDA3, FDXR | HIST1H1B, CLSPN, RRM2, CCNA2, MKI67, MCM10, KIF11, CDC20, TOP2A, AURKB | + HepG2 RPL17/RPL36/RPL23A/RPS19, PSMB5 (RPE1, HepG2) | **p53 activation / arrest vs proliferation (G2/M, S)** |
| 2 | 1.2% | ZNF524, STN1, TSEN2, SERPINE1, PSMB9, NFKB2, ICAM3 | MT1H, CAMK4, THBS4, CASP2, NNAT, PNMT | Feng fiaj_3 (4 rows), pipw_5 (13), tolg_4, oikd_2 | small-table artefact (no coherent program) |
| 3 | 1.0% | MT1E, CISH, DDIT3, FOXD1, SLIT3 | ADGRB2, TSEN2, STARD10, SETBP1, COL5A1 | Feng fiaj_3 / pipw_5; Tahoe AURKB inhibitors | small-table artefact, partly AURKB |
| 4 | 0.9% | IFI6, ISG15, POLD4, MT1E, ENO2 | TTK, PRPF3, FOSL1, CCNL1, CLK1, HK2, MSMO1 | − Tahoe AURKB (8 lines); + Feng, TGFBR1/2 (Jiang) | **AURKB-inhibitor response (Tahoe)** vs interferon-like |
| 5 | 0.8% | EFNB1, ZNF44, DDIT3, FOXD1, SOX21, IRX1 | ST8SIA2, LOXL1, GAS7, NMI, MECOM, NPY | Feng fiaj_3 / pipw_5 / oikd_2 / tolg_4 | small-table artefact |
| 6 | 0.6% | KLHL24, CREBRF, HBP1, CBLB, PDCD4, RBL2, EZH1, PDK1, APAF1 | NOP56, SRM, HSPE1, NOP16, NME1, MRTO4, EBNA1BP2, ATP5MC1, UBA52, COX6B1 | − Tahoe AURKB; Feng pipw_5 | **quiescence vs ribosome biogenesis / mitochondrial translation (MYC-like growth)** |
| 7 | 0.4% | MT1H, MT1G, MT2A, TRIM22, IFITM1, IFITM3, TXNIP | DDIT3, PPP1R15A, TRIB3, AREG, HSPH1, BRIX1, LSG1 | + Jiang IFNB USP18 (A549, BXPC3); − HepG2 EIF2S1, HSPA9, EIF3B | **metallothionein / interferon vs integrated stress response (ATF4)** |
| 8 | 0.4% | MT1H, HIST1H2BG, POLE2, OPTN | PKIB, LEFTY2, CTHRC1, PREX1 | Feng pipw_5; HepG2 EIF2S1/EIF2S3 | mixed; small-table dominated |
| 9 | 0.4% | ID2, NHEJ1, TRIM22, NRG1 | MSI1, FOXO4, GLI3, DGAT2 | Feng pipw_5; HepG2 splicing (DDX41, SF3B3, CACTIN) | mixed; small-table dominated |
| 10 | 0.4% | EPHB1, HDAC4, MET, FOSL1 | PLAUR, CIC, DYRK1B, HEG1 | Feng pipw_5 / fiaj_3; HepG2 EIF2S1/S3 | mixed; small-table dominated |

- **Equal weight per background lets tiny backgrounds steer the components.** Feng's fiaj_3 (4 rows) and pipw_5 (13 rows) each carry the same total weight as K562's 9,567 rows. Their rows have weight 1/2 and 1/3.6 against K562's 1/98, so PCs 2, 3, 5 and 8–10 are largely those few rows. Only PC1, PC6 and PC7 (and the Tahoe AURKB axis in PC4) read as biology.

**The 300 D/E/F targets on E1's PCs 1–5** (mean over the genetic tables that hold each target, per table first; `partE_DEF_target_scores.csv`). 299 have a score. The 90th percentiles of |score| over all genetic rows are 1.66 / 0.41 / 0.41 / 1.50 / 0.33. **6 targets exceed one of them:**

| Target | PC1 | PC4 | Other | (b) | HGNC group |
|---|---|---|---|---|---|
| DPY30 | −2.73 | −3.78 | PC2 −0.52, PC5 −1.32 | 1,191 | WRAD complex; DPY30 domain containing |
| PAPOLA | +1.72 | +0.51 | — | 127 | — |
| LTO1 | +0.79 | **+1.58** | — | 457 | — |
| PHF23 | −0.48 | −1.76 | PC3 +0.47 | 17.7 | PHD finger proteins |
| GPATCH1 | +0.33 | −0.96 | PC5 +0.44 | 25.5 | G-patch domain containing; spliceosomal Bact / C complex |
| LYSMD4 | −0.33 | −0.45 | PC5 −0.48 | 0.33 | LysM domain containing |

Averaging over tables pulls targets towards 0, so few clear the row-level p90.

**HGNC keyword counts among the 300.** The rules are exact HGNC group names (`step19.HGNC_E`):
- histones = "H1/H2A/H2B/H3/H4 histones";
- replication = "MCM family", "minichromosome maintenance 2-7 complex", "GINS complex", "DNA polymerases";
- mitosis = "Constitutive centromere associated network", "MIS12 kinetochore complex", "NDC80 kinetochore complex", "Anaphase promoting complex", "Cyclins";
- TP53 pathway = groups whose name contains p53/TP53 ("p53 transcription factor family", "TP53 regulating kinase family"). HGNC has no "p53 pathway" group.

Results: **histones 0; replication 1 (REV1, DNA polymerases); mitosis 1 (CCNE2, cyclins); TP53 0.** The D/E/F panel is not built around the proliferation / p53 machinery that E1's one clean program (PC1) describes.

### E-reg — program regression, scored on the proxy (`tables/partEreg.json`, `partEreg_rows.csv`, `partEreg_programs.csv`, `partEreg_DEF_Pprog.npz`; 18:02–18:19 UTC)
**Setup:**
- **Bases:** E1 for H1; E1−K562 for K562; E1−Orion for HCT116 and HEK293T.
- **ŝ_k:** per knockdown, each source background that holds k casts one vote: its tables' rows averaged, then projected (unweighted, centred, confidence-weighted, as in the basis). Votes are then averaged over backgrounds, with the fold's own lines excluded and ŝ = 0 where no background holds k. 291 of 292 H1 knockdowns have votes (K562 267/267, HCT116 276/278, HEK293T 269/271).
- **"Holds" means a genetic row.** The supplement says Tahoe's drug rows have no target gene, but on disk they carry the drug's annotated target. I left them out of the votes, as intended, which also keeps the standing rule that Tahoe drugs are never averaged.
- **Truth s_k:** H1's true LFC (gated genes) projected on the same K vectors by least squares.
- **β_j:** one slope per program through the origin. For H1, fitted on one half and applied to the other, both ways. For the guardrails, fitted on all 292 H1 knockdowns under the fold's basis, with the fold's lines removed from H1's votes too.
- **P_prog** = Σ_j β_j ŝ_j v_j on the line's gated genes.
- **K grid:** {5, 10, 20, 100}. K at 25% is not reached within 100 components on any of the three bases, so 100 stands in for it.

**The prediction itself, before any arm.** Per program, held out:

| K 5, program | β (fit h0 / h1) | r(ŝ, s) held out | R² held out |
|---|---|---|---|
| 1 (p53 vs proliferation) | **−0.32 / −0.40** | −0.18 / −0.05 | 0.10 / 0.04 |
| 2 | 0.32 / 0.70 | 0.24 / 0.14 | 0.03 / −0.01 |
| 3 | 0.91 / 0.48 | 0.14 / 0.31 | 0.01 / 0.07 |
| 4 (AURKB axis) | 1.39 / 1.34 | 0.34 / 0.49 | 0.12 / 0.23 |
| 5 | 0.37 / −0.14 | −0.05 / 0.15 | −0.04 / −0.01 |

- **At K 100:** β > 0 for 98% / 93% of programs, median held-out r 0.23 / 0.25, and R² > 0 for 66 / 85 of the 100.
- **PC1's β is negative in both directions:** H1 moves along the p53/proliferation axis opposite to the sources' average.

**Held-out captured energy of P_prog alone** (cross-fitted: each half from the other half's β), against the ceiling (E1, perfect coefficients) and config 1′:

| | (a) all / weak / mid / strong | (b) all / **strong** | (c) all / weak / **strong** |
|---|---|---|---|
| P_prog, K 5 | 0.002 / −0.002 / −0.002 / 0.003 | 0.003 / **0.003** | 0.003 / −0.006 / **0.005** |
| P_prog, K 10 | 0.003 / −0.007 / −0.003 / 0.006 | 0.006 / **0.006** | 0.004 / −0.018 / **0.008** |
| P_prog, K 20 | 0.004 / −0.012 / −0.005 / 0.009 | 0.010 / **0.010** | 0.006 / −0.031 / **0.013** |
| P_prog, K 100 | 0.017 / −0.040 / −0.027 / 0.035 | 0.043 / **0.045** | 0.022 / −0.088 / **0.048** |
| Ceiling E1, K 20 / 100 | 0.057 / 0.197 | 0.092 / 0.323 (strong 0.097 / 0.323) | 0.076 / 0.252 (strong 0.091 / 0.295) |
| Config 1′ | −0.109 | 0.042 / 0.050 | −0.040 / −0.610 / 0.066 |

- **The regression recovers about a tenth of the ceiling at K ≤ 20 and ~14% at K 100.**
- **On weak knockdowns it adds error,** because any non-zero prediction there costs energy.

**Proxy arms** (calls K0; P_prog enters as an LFC on m0, as config 1′'s called values do):
- **R1:** uncalled gated genes = m × P_prog.
- **R2:** config 1′ + γ × P_prog on all gated genes.

The table gives overall on the fit half / the held-out half, using β fitted on the fit half (base: half 0 0.2664, half 1 0.2247).

| Arm, scale | K 5 | K 10 | K 20 | K 100 |
|---|---|---|---|---|
| R1, 0.25 – 1.0 (fit h0 → h1) | 0.188–0.189 / 0.166–0.167 | 0.189–0.192 / 0.166–0.167 | 0.192–0.197 / 0.169–0.172 | 0.202–0.219 / 0.177–0.190 |
| R2, 0.25 (fit h0 → h1; fit h1 → h0) | 0.2667 / 0.2251; 0.2251 / 0.2666 | 0.2659 / 0.2238; 0.2248 / 0.2656 | 0.2667 / 0.2270; 0.2251 / 0.2667 | 0.2716 / 0.2303; 0.2291 / 0.2696 |
| R2, 0.5 | 0.2657 / 0.2254; 0.2249 / 0.2656 | 0.2632 / 0.2245; 0.2244 / 0.2651 | 0.2652 / 0.2261; 0.2261 / 0.2669 | **0.2723 / 0.2340**; 0.2329 / 0.2719 |
| R2, 0.75 | 0.2651 / 0.2257; 0.2253 / 0.2650 | 0.2608 / 0.2230; 0.2244 / 0.2637 | 0.2629 / 0.2245; 0.2257 / 0.2656 | 0.2702 / 0.2304; **0.2360 / 0.2715** |
| R2, 1.0 | 0.2653 / 0.2256; 0.2245 / 0.2651 | 0.2586 / 0.2217; 0.2238 / 0.2625 | 0.2623 / 0.2233; 0.2271 / 0.2646 | 0.2661 / 0.2260; 0.2339 / 0.2697 |

**R1 collapses** (~0.18 on H1, 0.04–0.06 on every guardrail): replacing config 1′'s uncalled values removes the pooled-profile baseline they carry.

**Fit and guardrails:**

| Fit half | Pick | Grid positions | Held-out overall vs config 1′ | Gain: overall / MSE / PDS / NMAE / fidelity / reach / Jaccard | K562 / HCT116 / HEK293T vs base |
|---|---|---|---|---|---|
| 0 | R2, γ 0.5, K 100 | arm 2/2 (edge), γ 2/4 (interior), K 4/4 (**edge**) | 0.2340 vs 0.2247 (half 1) | **+0.0094** / +0.036 / −0.021 / +0.020 / +0.000 / +0.022 / 0 | +0.0035 / **−0.0117** / +0.0030 |
| 1 | R2, γ 0.75, K 100 | arm 2/2 (edge), γ 3/4 (interior), K 4/4 (**edge**) | 0.2715 vs 0.2664 (half 0) | **+0.0051** / +0.048 / −0.042 / +0.020 / +0.000 / +0.004 / 0 | +0.0056 / **−0.0174** / +0.0029 |

- Every held-out component is above 0.
- The proxy overstates reach by 0.02–0.05 compared with cells (15c §4). About 40% of fit 0's gain is reach (+0.022 of the components' +0.056 sum), so on cells it may shrink to about the MSE + NMAE − PDS part (~+0.006).
- HCT116's loss is all PDS (0.844 → 0.765 / 0.738).

**Criterion** (≥ +0.005 held out on both halves, every component > 0, no guardrail more than 0.005 worse): **not met.** The two held-out gains pass (+0.0094 / +0.0051, the second only just), but HCT116 is 0.012–0.017 worse.

**D/E/F.** P_prog for the 300 targets per context (all sources including Arc H1 as one more background, E1 basis, β fitted on all of H1, each context's own > 5 CPM gate; K 5 / 10 / 20 / 100) is saved in `tables/partEreg_DEF_Pprog.npz`. 299 targets have votes, with a median of 4 backgrounds. Nothing was built from it.

## Part F — recommendation for the next D/E/F file
**One spec, the same for D, E and F.** Nothing in Step 19's results separates the three contexts on calls or magnitudes:
- The targets, and so (b), are shared.
- None of the contexts has a close source line (Part C).
- D's larger baseline offset (0.060) is already handled by config 1′'s pooled-profile target.

| Item | D | E | F | Evidence |
|---|---|---|---|---|
| **Call rule** | **Hybrid, T 5:** K4 1,200 where b ≥ 5, else K0 (46% of the 300 targets on K4; ANG and RAB5A, with no b, on K0) | same | same | Part A (cells) |
| **Called-gene magnitudes** | **config 1′** (0.5 × C+G LFC); no tier row, no gate row | same | same | Part D, gate rows (proxy) |
| **Uncalled genes** | **config 1′** (no E-reg arm) | same | same | E-reg (proxy) |
| **No-source fix** | **include:** `n_sources` ≥ 1 for gated genes with no source, so they keep config 1′'s pooled-profile value | same | same | 15c §4; here +0.0013 (Eval A), +0.0005 (HCT116final), proxy |

**Why T 5 and not the fit's pick.**
- **The held-out fit picked T 2 (3/12, interior) and K4 everywhere (1/12, edge).** Both pass Part A's criterion. On Eval A, T 5 scores the same as those picks within noise:
  - full +0.0201 [+0.0155, +0.0251] against T 2 +0.0195 and K4 +0.0213;
  - halves 0.2835 / 0.2418 against T 2 0.2838 / 0.2403.
- **On HCT116final it costs half as much:** −0.0070 [−0.0124, −0.0017], against −0.0125 at T 2 and −0.0140 at K4.
- **This choice uses the guardrail,** so it is a judgment, not a held-out fit.
- **Part B argues D/E/F's truth will be denser than the guardrails' at the same strength** (400 cells per knockdown and Eval A's depth, against 144–188 cells on the guardrails). That favours more K4. Their targets' (b), which decides who gets K4, looks like HCT116's (b ≥ 5: 46% against 49%).
- **The rules risk from 15c §3 still applies** to the targets that go to K4: ~700 extra calls per knockdown with almost no requested shift, and generated cell SD off by up to 2×. **Whether to accept it is your call.** Without it, the spec is the uploaded file plus the no-source fix.

**Expected change against the uploaded file** (config 1′ + K0, md5 `ec8ac09a…`; Eval A 0.2426, HCT116final 0.3617 on cells):

| Component of the change | Eval A | HCT116final | Source |
|---|---|---|---|
| Hybrid T 5 calls | **+0.0201** [+0.0155, +0.0251] | **−0.0070** [−0.0124, −0.0017] | **cells** (Part A, paired, 2,000 resamples) |
| No-source fix | +0.0013 | +0.0005 | proxy (generator input) |
| **Expected total** | **≈ +0.021** (→ ≈ 0.264) | **≈ −0.007** (→ ≈ 0.355) | cells + proxy |

**Do the parts combine?** On the proxy (`tables/partF_rows.csv`), as one row next to each part alone (H1 full / half 0 / half 1; fold pools full):

| Row | H1 full | H1 halves | K562 | HCT116 (fold) | HEK293T |
|---|---|---|---|---|---|
| base (uploaded config) | 0.2456 | 0.2664 / 0.2247 | 0.1582 | 0.3108 | 0.1536 |
| no-source fix alone | 0.2469 (+0.0013) | 0.2678 / 0.2259 | 0.1582 (0) | 0.3108 (0) | 0.1536 (0) |
| hybrid T 5 alone | 0.2669 (+0.0213) | 0.2878 / 0.2458 | 0.1471 (−0.011) | 0.3017 (−0.009) | 0.1530 (−0.001) |
| **hybrid T 5 + fix (the spec)** | **0.2682 (+0.0226)** | 0.2892 / 0.2471 | 0.1471 (−0.011) | 0.3017 (−0.009) | 0.1530 (−0.001) |
| hybrid T 2 + fix (the fit's pick, for comparison) | 0.2684 (+0.0228) | 0.2898 / 0.2468 | 0.1437 (−0.015) | 0.3001 (−0.011) | 0.1492 (−0.004) |

- **The two parts add exactly** (+0.0213 + 0.0013 = +0.0226): the fix moves MSE/PDS, the calls move fidelity/Jaccard.
- **The fold pools can't show the fix,** because their inputs already carry config 1′'s intended value. Its HCT116final number comes from that line's generator input (0.3679 → 0.3684).
- **The proxy's hybrid gain on H1 (+0.0213) is close to the cells' (+0.0201).**
- **Its guardrail costs are larger than HCT116final's on cells.** The fold pool holds out both Orion lines.

**Cell confirmation run:** not needed for the call rule.
- **The hybrid is a cell result.** The generator works knockdown by knockdown, so a hybrid's cells are the two 15c runs' cells. Part 0 reproduced both endpoints to < 1e-15, and Part A recomputed the panel-wide MSE budget.
- **The no-source fix changes ~0.9% of gated inputs** and, per your note, needs no test.
- **The K4 targets may generate a little slower** (1,200 shaped genes); unmeasured.

**Does Part E justify a regression?**
- **The ceiling criterion is not met:** E1 holds < 10% of strong-knockdown energy at K ≤ 20, at the bottom of S2's range.
- **E-reg misses its criterion only on HCT116** (PDS −0.08 / −0.11). Its Eval A gain (+0.009 / +0.005 held out) needs K 100, a grid edge, and ~40% of it is reach, which the proxy overstates.
- **So: not for the next file.** Pursuing it before 5 Nov would take:
  - a proxy round to explain HCT116's PDS loss and to try E2 (whose ceiling is 2× E1's), ~1–1.5 h;
  - a D/E/F build, ~2 h.

  That is ~7 h of compute plus review, feasible by 3–4 Nov but with little margin for a second upload. P_prog for D/E/F is saved if you want it.

## What I would have done differently
2. **Raised the background weighting before running Part E.** Equal weight per background gives a 4-row Feng table the same total weight as K562's 9,567 rows, and those few rows steer 6 of E1's top 10 components. A minimum background size (or capping weight at, say, 1/√50) is the obvious variant; I reported the effect rather than changing the spec.
3. **Checked the Tahoe `target_gene` premise against the data before coding E-reg and the gate rows.** The supplement says drug rows have no target gene; on disk they carry the drug's annotated target. I excluded them as intended, but should have flagged it to you before running.

5. **Avoided a pandas method name as a column** (`assign`), which cost Part A one crash.
6. **Compared Part A's thresholds against K4 at matched mean calls too.** As in Step 17, Eval A's overall falls with the call volume as T rises, so most of the T sweep is volume, not allocation. The criterion compares against K0, which hides that.

## Follow-up (Jose, 2026-11-02): MSE by true-strength tercile, and the programs in the mid/strong residual
*Code: `step19.py mseterc`, `mseterc2`. Tables: `results/stage1/2026-11-01_step19/tables/mse_terciles.json`, `mse_terciles_capture.csv`, `mse_terciles_corrections.csv`; logs `mseterc*.log`. Eval A only; terciles as 15c §1b (rank thirds of the true DE count: weak median 5, mid 215, strong 2,490). Each tercile's MSE is the cell-eval2 ratio of sums over its knockdowns, scaled with the **full panel's** anchors (as boot_score does), so these are not per-tercile leaderboard numbers. No cells generated.*

**1. MSE by tercile (scaled; unclamped and raw in brackets).** The uploaded file's full-panel raw MSE reproduces its `summary.json` (|diff| 1e-16).

| Prediction | All | Weak | Mid | Strong |
|---|---|---|---|---|
| Bundle baseline (anchor) | 0 (0.967) | 0 (−1.58; 2.481) | 0 (−0.04; 1.002) | 0.026 (0.942) |
| C+G (cells) | 0 (−0.55; 1.496) | 0 (−20.4; 20.60) | 0 (−2.66; 3.52) | 0 (−0.04; 1.004) |
| **Config 1′ + K0, uploaded (cells)** | **0.127** (0.845) | **0.971** (0.034) | **0.067** (0.902) | **0.121** (0.850) |
| Config 1′ + K0 (proxy) | 0.127 | 1.0 (+1.02) | 0.070 | 0.120 |
| E-reg R2 γ 0.5 K 100, cross-fitted (proxy) | 0.165 | 1.0 (+1.37) | 0.111 | 0.154 |
| E-reg R2 γ 0.75 K 100, cross-fitted (proxy) | 0.163 | 1.0 (+1.13) | 0.058 | 0.161 |

**Where the panel's MSE comes from:**

| Share of… | Weak | Mid | Strong |
|---|---|---|---|
| the denominator Σ den (real effect energy) | **1.2%** | 9.8% | **88.9%** |
| the uploaded config's numerator | 0.0% | 10.5% | 89.5% |

- **Your recollection holds:** the weak tercile is at 0.971 scaled on cells, essentially solved.
- **But it carries 1% of the denominator, so the panel MSE is the strong tercile's MSE.** Mid adds a tenth; any MSE gain has to come from strong knockdowns.
- **E-reg's MSE gain is mostly mid/strong** (strong 0.121 → 0.154–0.161), but it costs PDS (Part E).

**2. What the uploaded config's residual (P − R, scorer space, own gene out) is made of:**

| Tercile | Called genes | Uncalled gated | Ungated | Tercile-mean residual (a shared shift) | Panel-mean residual |
|---|---|---|---|---|---|
| weak | 23% | 68% | 10% | 2.7% | 9.0% |
| mid | 19% | 71% | 11% | 4.5% | 5.6% |
| strong | 13% | **78%** | 8% | **7.7%** | 1.6% |

The strong tercile's error is on the genes we don't call, and part of it is a shift shared by strong knockdowns that the panel-wide mean doesn't contain (7.7% against 1.6%).

**3. How much of that residual the source programs hold** (least-squares projection on the top-K basis vectors with **perfect** per-knockdown coefficients, i.e. a ceiling; the MSE it would leave in brackets, scaled, cells' Ĉ_pred, budget recomputed):

| Basis | K | Captured: weak / mid / **strong** | MSE all / mid / **strong** if removed |
|---|---|---|---|
| E4 (Δ space) | 5 | 1.5% / 2.9% / **4.6%** | 0.218 / 0.265 / **0.187** |
| E4 | 20 | 2.4% / 5.8% / **10.5%** | 0.286 / 0.346 / **0.253** |
| E4 | 96 | 5.3% / 15.4% / **29.8%** | 0.510 / 0.615 / **0.469** |
| E1 (LFC) | 20 / 100 | 1.9% / 3.9% / 7.0%; 4.1% / 11.5% / 22.8% | 0.245 / 0.291 / 0.214; 0.427 / 0.506 / 0.391 |
| Random | 50 | 0.3% each | 0.133 / 0.084 / 0.125 |

- **No single source program carries much.** The ten E4 programs with the most mid+strong residual hold 0.7–1.4% each.
  - The largest is E4's p53-vs-proliferation axis (+ CDKN1A, MDM2, ZMAT3, PHLDA3; − MKI67, TOP2A, CENPF, CCNB1; top rows are RPE1 ribosomal knockdowns). Its mean coefficient is negative in 62% of mid/strong knockdowns: truth is more p53-active / less proliferative than predicted.
  - Most of the rest are again driven by the small Feng tables (Part E's weighting caveat). None of their coefficients track (b) (|Spearman| ≤ 0.23).

**4. H1's own residual programs** (SVD of the mid+strong residual on one H1 half, applied to the other half's mid+strong knockdowns):

| Fit half | Held-out residual captured, K 1 / 2 / 5 / 10 / 20 (oracle coefficients) |
|---|---|
| 0 | 6.7% / 11.4% / **21.4%** / 26.7% / 30.6% |
| 1 | 7.4% / 13.1% / **22.4%** / 29.1% / 33.6% |

- **The residual is far more concentrated than the source programs suggest:** 5 of H1's own directions hold ~22% of held-out strong/mid residual, against 4.6% for E4's top 5.
- **The leading ones have a consistent sign:** PC1's held-out coefficient has the same sign in 77% / 90% of knockdowns.
- **They are only partly in the sources:** 34–57% of each PC's norm² lies in E4's 96-dimensional span (a random 96-dimensional subspace would hold 0.5%), and 20–48% in E1's. Max |cos| with any single E4 vector is 0.18–0.35.

The PCs, as **my reading** of their top genes (they repeat across the two halves):

| Program | Up in truth relative to prediction | Down in truth relative to prediction | Label (my reading) |
|---|---|---|---|
| PC1 (both halves) | LEFTY1, LEFTY2, NODAL / CER1, H1F0, CDH2, FST, SPP1, FZD2, CHD7 | DNMT3B, DPPA2, MT1X, GAL, SCGB3A2, SLC7A5 | **exit from pluripotency / Nodal–Lefty differentiation response** |
| PC2 (half 0) | p53 targets (CDKN1A, BBC3, PHLDA3), FN1, VIM, CAV1 | RRM2, PCLAF, DUT, CDCA7, ODC1, NOP16 | **p53 / arrest vs S-phase proliferation** |
| PC4 (half 0); PC2 (half 1, − side) | CHAC1, PCK2, STC2, DDIT4, ASS1, SLC7A11, PSAT1, SESN2, ATF5, INHBE, VEGFA | — | **integrated stress / amino-acid response (ATF4)** |
| PC3 / PC4 (half 1) | ACSS2, MVD, HMGCS1, ACAT2, PCSK9, FADS2, INSIG1, ELOVL6, ERG28 | — | **sterol synthesis (SREBP)** |

The direction ("up in truth") follows from the sign of the held-out mean coefficient on P − R. For PCs 3–4 the sign varies by half, so read those as an axis.

**5. Do corrections learned on one half lower the other half's MSE?** Uploaded cells, held-out half, scaled (fit half 0 → half 1 / fit half 1 → half 0):

| Correction (subtracted from the prediction) | All | Mid | Strong |
|---|---|---|---|
| none (uploaded) | 0.108 / 0.144 | 0.000 / 0.168 | 0.115 / 0.127 |
| fit half's tercile-mean residual, assigned by **true** tercile | 0.207 / 0.229 | 0.119 / 0.330 | 0.196 / 0.193 |
| the same shift, assigned by **(b)** tercile (no truth on the target) | **0.167 / 0.198** | 0.000 / 0.194 | **0.178 / 0.184** |
| top-1 own residual PC, fit-half mean coefficient by (b) tercile | 0.134 / 0.190 | 0.000 / 0.186 | 0.143 / 0.174 |
| top-K own PCs with oracle coefficients, K 5 / 20 (ceiling) | 0.420 / 0.459; 0.540 / 0.599 | 0.337 / 0.541; 0.462 / 0.731 | 0.412 / 0.428; 0.532 / 0.565 |

- **A single shared shift for strong and mid knockdowns, learned on half the panel and assigned by (b), raises held-out MSE by +0.055–0.06 scaled.** That is ≈ +0.01 overall, as much as E-reg's best arm, from one vector per tercile.
- **Placed by the true tercile it is worth ~+0.09.** The difference is (b)'s misassignment, which on half 1 also zeroes the mid tercile.
- **Per-knockdown coefficients on five such programs would be worth ~+0.3 (a ceiling).** That is where most of the reducible mid/strong MSE is.

**What this means for D/E/F:**
- **The biggest and most consistent direction is H1-specific.** Strong H1 knockdowns push cells out of pluripotency (LEFTY/NODAL/H1F0 up, DNMT3B/DPPA2 down), and D, E and F are not pluripotent (Part C). That shift cannot be transferred as is.
- **The others are generic:** p53/arrest vs proliferation, the ATF4 stress response and sterol synthesis. Each exists in the source lines (E4/E1 span 34–57%), so a "shared strong-knockdown response" built from the sources could carry them.
- **The natural next proxy test (not run):**
  - estimate each source line's own mean response of its strong knockdowns (by its own true DE count) minus its weak ones, in scorer/Δ space;
  - pool across lines, one vote per background;
  - add γ × that vector to knockdowns by (b) tercile;
  - score on H1 halves + guardrails, against this table's +0.055 that H1's own shift gets.

  If the source-built shift recovers a good part of it on H1 without hurting HCT116's PDS, it is a cheap, transferable MSE lever.
- **The cost to watch is PDS.** A shift common to many knockdowns makes them more alike, which is what cost E-reg on HCT116.

---

# 12. 2026-11-02 — Step 20: lab effects, the strong-knockdown shift and program regression v2
*Task: Jose, 2026-11-02.*
- **Not done:** no cells were generated, no new cells were scored, no D/E/F file was built and nothing was uploaded.
- **Code:** `code/vcc6/step20.py` (chain scripts `step20*.sh` alongside).
- **Outputs:** `results/stage1/2026-11-02_step20/{tables,logs}`.
- **Frozen inputs** were imported or read only: the generator, the harness, `run_line.py`, `boot_score.py`, config 1′'s files, v7 and the Step 18 / 19 outputs.

**Purpose:** four choices for the next D/E/F file:
1. the call rule;
2. a shift for strong knockdowns;
3. Arc H1's weight for the targets it holds;
4. a program-regression term on the uncalled genes.

**The hypothesis under test:** lab protocol shapes knockdown responses more than biology does. If so, D/E/F's truth, which comes from Arc's protocol, should resemble H1's.

**Unit tests** (`step20.py test`, `testF`; synthetic slices, run before each chain):
- the (b) cut points and terciles;
- the Δ-space addition and its round trip;
- pair cosines (masked, and with both targets' genes out);
- duplicate collapsing, pair classes and line / lab maps;
- the table-bootstrap regression (it recovers known coefficients);
- weighted votes and agreement z;
- the strong-minus-weak vector, S with exclusions, the projection;
- `pick_fit`;
- the weighting fix;
- weighted ridge, which matches duplicated rows, with CV and constant columns;
- the design matrix, per-row least squares, background scores, the residual correlation, energy sums, the artefact counts, HGNC mapping and gene features.

**Summary:**
- **Lab against biology (A):** in 513 source pairs, sharing a lab is the only sizeable term (+0.10 in pair cosine on strong knockdowns), but its CI touches 0 [−0.012, +0.130].
  - Sharing a line or a lineage adds ≤ +0.02.
  - **H1 agrees with Weissman's lines as well as they agree with each other.** It agrees no better with the one other-lab iPSC screen that shares many strong knockdowns with it.
  - The data can't test Arc's own protocol, because Arc has one line.
- **Same-lab upweighting (B):** misses its criterion.
  - A same-lab line is worth +0.06 to +0.08 as an ordinary voter.
  - Upweighting it lowers overall and PDS on both Orion analogues.
  - The K562 analogue can't test it: the upweighted lines hold none of K562's targets.
  - **Arc H1's weight stays at 1.**
- **Strong-knockdown shift (C, D):** misses its criterion (+0.0034 / +0.0060 held out).
  - What gain there is comes from pushing weak knockdowns away from S (PDS), not from shifting strong ones.
  - The follow-up's H1-own shift gains MSE but loses more PDS on the proxy (−0.002 / −0.006).
  - S_H1 is 71% the pluripotency program. Without it, it does nothing on the guardrails.
- **Program regression v2 (F):**
  - **The weighting fix works.** E4w's strong-tercile ceiling at K 20 is 0.231, against E1's 0.091, and the small-table artefacts are gone.
  - **Ridge on separate source lines doubles Step 19's prediction.** DepMap adds a little, gene features alone beat their null, and the co-expression prior and the gene embeddings carry almost nothing.
  - **On the proxy, R3** (the term only where b ≥ 5) **with E4w, M3, DepMap, K 100 and γ 1 gains +0.025 / +0.023 held out, and +0.007 to +0.011 on every guardrail.**
  - **F5's criterion is still not met.** The fit rule picked the pooled M3p in one direction, and it fails held out.
  - ARCHS4 was not used (see the addendum).
- **Recommendation (G):**
  - **Spec:** hybrid T 5 + config 1′ + no-source fix (Step 19's spec), unchanged. Expected Eval A ≈ +0.021, HCT116final ≈ −0.007.
  - **The F5 program term** (+0.021 more on the proxy) **is the one candidate worth a cell check** before the deadline. That's your call.

## Part 0 — checks (`tables/part0.json`, `logs/part0.log`)
Both pass.

| Check | Got | Step 19 / follow-up | max \|diff\| |
|---|---|---|---|
| 1. Config 1′ + K0 + no-source fix, H1 generator input: full / half 0 / half 1 | 0.2469 / 0.2678 / 0.2259 | 0.2469 / 0.2678 / 0.2259 (Part F `nosource_fix` row) | 5.6e-17 |
| 1. Same row on the fold pools: K562 / HCT116 / HEK293T | 0.1582 / 0.3108 / 0.1536 | same | ≤ 8.3e-17 |
| 2. H1-own tercile shift on the uploaded cells, by (b): held-out MSE, fit half 0 / 1 | 0.108 → 0.167 / 0.144 → 0.198 | 0.108 → 0.167 / 0.144 → 0.198 | < 1e-9 |

## Part A — lab against biology in the source tables (`tables/partA.json`, `partA_pairs.csv`, `partA_labs.csv`)

### A1. Lab map (67 genetic tables; Tahoe out)
| Lab | Tables | Lines | Sure? |
|---|---|---|---|
| Weissman | 7: Replogle K562 GW, Replogle RPE1, Adamson K562 ×3, Nadig HepG2, Nadig Jurkat | K562, RPE1, HepG2, Jurkat | sure |
| Xaira (X-Atlas/Orion) | 4: Orion HCT116, HEK293T and their pseudobulk complements | HCT116, HEK293T | fairly sure (Orion = X-Atlas/Orion) |
| Arc | 1: Arc H1 | H1 | sure |
| Feng | 21: Feng genome-wide (fitness, non-fitness; HipSci donors pooled) and 19 targeted tables (HipSci donor clones) | 10 donors + HipSci pooled | **unsure of the lab's name**: our metadata records the study only. The tables share one lab and protocol |
| KOLF2.1J study | 1: KOLF2.1J CRISPRi (`KOLF_Strong_Perturbations`) | kolf | **unsure**: the producing lab isn't recorded in our download notes, so it is treated as its own lab |
| Jiang | 30: 6 lines × 5 cytokine conditions, CRISPRko, Parse split-pool | A549, BXPC3, HAP1, HT29, K562, MCF7 | **unsure** it is the Satija lab |
| Marson | 3: CD4 T rest / stim 8 h / stim 48 h (replicate tables) | CD4 T | sure |

**Line and lineage rules** (stated before the run):
- Feng clones of one donor count as one line, so eipl_1 and eipl_3 are the same line.
- KOLF2.1J (a gene-corrected subclone of HPSI0114i-kolf_2) and Feng's kolf_2 and kolf_3 are all line "kolf".
- Lineages:
  - **pluripotent:** H1, HipSci and KOLF;
  - **blood:** K562, Jurkat, CD4 T;
  - **colorectal:** HCT116, HT29;
  - **one lineage each:** HEK293T, HepG2, RPE1, A549, BXPC3 and MCF7. HAP1 is counted as its own lineage, not as blood.

### A2–A3. Pairs and measures
- **Scope:** 696 pairs of genetic tables share at least 20 knockdowns.
- **Marson's 183 pairs drop out of every cosine.** Its baseline has no control CPM (replicate tables, no cells), so the "> 5 CPM in both controls" gate can't be formed. That leaves 513 pairs with a cosine. I left them out rather than invent a gate.
- **Measures:** each pair's cosine is the mean over its shared knockdowns of the per-knockdown cosine.
- **"DE ≥ 20":** knockdowns whose source DE count (`n_real`) is at least 20 in both tables.

| Class | Pairs (with a cosine) | Knockdowns (DE ≥ 20) | Conf.-weighted LFC cosine, all knockdowns | **Conf.-weighted, DE ≥ 20** | Δ cosine, DE ≥ 20 | Spearman of `n_real` |
|---|---|---|---|---|---|---|
| same lab, same line | 11 (8) | 33,300 (198) | 0.061 [0.055, 0.070] | **0.243** [0.233, 0.253] | 0.325 | 0.66 |
| same lab, same lineage, different line | 163 (163) | 47,912 (4,842) | 0.061 [0.040, 0.071] | **0.240** [0.228, 0.252] | 0.319 | 0.59 |
| same lab, different line | 84 (84) | 26,912 (4,964) | 0.033 [0.026, 0.057] | **0.117** [0.045, 0.172] | 0.188 | 0.42 |
| different lab, same line | 7 (7) | 410 (36) | 0.029 [0.027, 0.038] | **0.094** (3 pairs) | 0.130 | 0.11 |
| different lab, same lineage | 48 (24) | 27,185 (241) | 0.038 [0.031, 0.047] | **0.101** [0.097, 0.127] | 0.143 | 0.02 |
| different lab, different lineage | 383 (227) | 148,680 (3,849) | 0.033 [0.018, 0.052] | **0.134** [0.108, 0.143] | 0.168 | 0.20 |

Medians [IQR] over pairs. All-knockdown cosines are small, because most shared knockdowns are weak in at least one table.

**The classes are dominated by single labs, so they must be read lab by lab:**
- **Feng:** 161 of the 163 same-lab, same-lineage pairs are Feng–Feng (0.240). Feng tables agree as well across donors as within one donor (same donor 0.235, 7 pairs).
- **Jiang:** 75 of the 84 same-lab, different-line pairs are Jiang–Jiang (0.090; shallow CRISPRko at ~3.9k UMI).
- **Weissman, 10 pairs, 740–1,041 strong shared knockdowns each:**
  - most cross-line pairs: 0.17–0.21 (K562–Jurkat 0.201, K562–RPE1 0.171, HepG2–Jurkat 0.171, Jurkat–RPE1 0.182);
  - the two exceptions: HepG2–RPE1 0.338, and Adamson–Replogle K562 (same line) 0.367.

**Key contrasts:**
- **Same lab / different line (0.117, 70 pairs) against different lab / same line (0.094, 3 pairs).** The second class is tiny:
  - KOLF2.1J against Feng's kolf_2 and kolf_3: 0.094 / 0.108;
  - Replogle K562 against Jiang K562: only one condition has strong shared knockdowns, at 0.039.
- **H1 against other labs' pluripotent tables (0.205, 2 pairs) against H1 and non-pluripotent tables (0.169, 6 pairs).**
  - The H1 pairs that carry weight are the ones with many strong shared knockdowns:
    - KOLF2.1J 0.170 (77 strong knockdowns);
    - Replogle K562 0.195 (91);
    - Replogle RPE1 0.184 (73);
    - Nadig Jurkat 0.189 (60);
    - Nadig HepG2 0.155 (35).
  - The pluripotent edge comes from Feng GW fitness: 0.239, on only 6 strong knockdowns.
  - **H1 agrees with Weissman's lines as well as Weissman's lines agree with each other, and no better with the one other-lab iPSC screen that shares many strong knockdowns with it.**

**Regression** of pair cosine (confidence-weighted, DE ≥ 20; 305 pairs, 55 tables, Marson out) on same lab, same line, same lineage and log median cells of both tables. 95% CIs are from a bootstrap over tables (1,000 draws):

| Term | β | 95% CI |
|---|---|---|
| same lab | **+0.096** | [−0.012, +0.130] |
| same line | +0.010 | [−0.030, +0.100] |
| same lineage | +0.019 | [−0.015, +0.119] |
| log cells, smaller table | −0.063 | [−0.103, −0.002] |
| log cells, larger table | −0.001 | [−0.072, +0.035] |

- **Same lab is the only sizeable term,** at about +0.10, but its CI touches 0.
- **Line and lineage add about +0.01–0.02,** with wide CIs.
- **The cells terms have the "wrong" sign.** The best-powered tables (Jiang and Adamson at 400 cells, shallow) agree least. Cells per knockdown is a poor proxy for power across protocols that differ ~5× in depth.
- On all knockdowns (513 pairs), every term is within ±0.02 and every CI includes 0.

**Reading (mine):**
- In the source data, sharing a lab and protocol moves agreement more than sharing a line or a lineage does, but by an amount the bootstrap can't separate from 0. Lineage on its own does nothing across labs: different lab / same lineage 0.101, against 0.134 for different lineage.
- **For D/E/F:** the evidence that Arc's protocol, more than biology, sets the response is weak in direction and unmeasurable for Arc itself, because Arc has only one line in the sources. What the data do say is that H1 transfers to other lines (and from them) about as well as Weissman's lines transfer among themselves.

## Part B — same-lab upweighting on two analogues (`tables/partB.json`, `partB_rows.csv`, `partB_pools.json`, `partB_pool_<line>_w<w>.npz`)
**Pools** (one table load, 117 tables, 6 min; then 0.09–0.11 s per knockdown for both pools):
- **HCT116:** HCT116's tables held out, so HEK293T votes (115 tables). HEK293T (Orion single-cell + pseudobulk) holds 260 of HCT116final's 278 targets.
- **HEK293T:** HEK293T's tables held out (115 tables). HCT116 holds 260 of its 271 targets.
- **K562:** every K562 table held out (108 tables), which is exactly the K562 fold pool.
  - **None of K562final's 267 targets is held by RPE1, HepG2 or Jurkat,** the three upweighted Weissman lines. So **w has no effect on K562**: every row is identical to the fold pool.
  - The "held" split for K562 is therefore empty, and its row in `partB_rows.csv` (−1.54) is meaningless.
  - (The HepG2 background also includes Tahoe's HepG2/C3A drug rows, which stage1 aliases into it; it holds none of the K562 targets either.)

**Weighting:** as Step 7's w_arc. The upweighted backgrounds' vote variance is divided by w, in both C+G's LFC pool and the Δ pool. Config 1′'s agreement z is recomputed from the new pool's unweighted Δ votes.

**Check:**
- The K562 pool at w 1 reproduces `fold_K562.npz` and `s11_K562.npz` exactly (LFC, calls, Δ pool, n_bg and z all identical).
- Its proxy score reproduces 15a (0.158234, |diff| 2.8e-17).

**Results** (proxy, config 1′ + K0; the fold pools carry the no-source fix already):

| Target, w | 1 | 2 | 5 | 20 | Fold pool (both Orion lines out), for reference |
|---|---|---|---|---|---|
| **HCT116** overall | **0.3684** | 0.3533 | 0.3419 | 0.3257 | 0.3108 |
| PDS | 0.864 | 0.831 | 0.791 | 0.731 | 0.844 |
| MSE / NMAE / fidelity / reach / Jaccard at w 1 | 1.000 / 0.223 / 0.916 / 0.780 / −1.571 | | | | 0.876 / 0.221 / 0.871 / 0.651 / −1.598 |
| knockdowns HEK293T holds (260) | 0.3687 | 0.3526 | 0.3405 | 0.3155 | |
| knockdowns it doesn't hold (18) | 0.3623 | 0.3623 | 0.3041 | 0.1957 | |
| **HEK293T** overall | **0.2335** | **0.2359** | 0.2085 | 0.1890 | 0.1536 |
| PDS | 0.625 | 0.607 | 0.581 | 0.565 | 0.424 |
| knockdowns HCT116 holds (260) | 0.2377 | 0.2381 | 0.2061 | 0.1906 | |
| knockdowns it doesn't hold (11) | −0.878 | −0.852 | −0.815 | −0.815 | |
| **K562** overall | 0.1582 | 0.1582 | 0.1582 | 0.1582 | 0.1582 |

- Subsets are scaled with the full panel's anchors, so the small "doesn't hold" subsets are noisy.
- Their MSE still moves with w, because cell-eval2's MSE budget is formed on the whole panel.

**Criterion** (overall rises with w, without PDS falling, on all three targets): **not met.**
- K562 is flat (untestable).
- HCT116 falls from w 1.
- HEK293T rises by +0.002 from w 1 to w 2, then falls.
- PDS falls with every step of w on both Orion lines.

**What the analogue says for Arc H1 and D/E/F:**
- **A same-lab source is very valuable as an ordinary voter.** Adding HEK293T to HCT116's pool, or HCT116 to HEK293T's, gains +0.058 / +0.080 over the pools without the Orion lines, and HCT116's MSE reaches the ceiling (1.000).
- **Upweighting it beyond its pooled precision hurts,** and mostly through PDS: it makes knockdowns look more like the source line's and less like each other.
- This matches Step 7's w_arc result.
- **So Arc H1's weight should stay at 1 for the 28 D/E/F targets it holds.**

**Arc H1's coverage of D/E/F:** it holds **28 of the 300 targets**. Over the targets × each context's gated genes, it has a finite LFC in **9.2%** of the pairs in all three contexts (D 266,084 / (300 × 9,633); E 9.2%; F 9.2%).

## Part C — a source-derived shift for strong knockdowns (`tables/partC.json`, `partC_rows.csv`, `partC_S.npz`, `partC_lines.npz`)
**Building S:**
- **Lines:** 33 genetic lines with full profiles, keyed by stage1 background, one vote each.
  - Left out: Marson CD4 T (no full profiles) and the Orion pseudobulk rows (no `n_real`).
  - **Feng's iPSC lines are 20 of the 33 votes** (19 donor clones plus HipSci pooled). Jiang has 5 more lines of its own (its K562 tables join the K562 line), Weissman 4 (K562, RPE1, HepG2, Jurkat), Orion 2, KOLF2.1J 1 and Arc H1 1.
  - Per line, knockdowns run from 4 (Feng fiaj_3) and 13 (pipw_5) to 9,294 (K562). So the smallest lines' "terciles" hold one to four knockdowns.
- **Each line's vector:** within a line, a knockdown's Δ is averaged over the line's tables first, using each table's own > 5 CPM gate and with the knockdown's own gene left out. The vector is the mean over the strong tercile minus the mean over the weak tercile of the line's own `n_real`.
- **S:** the mean of the line vectors over the lines that gate the gene.
- **Fold exclusion:** each scored line uses S without its own lines (H1: without h1; K562: without k562; HCT116 / HEK293T: without both Orion lines). D/E/F use all lines.

**Applying it:**
- Δ_pred + a_t × S on every gene where S is defined, except each knockdown's own gene. It is added to config 1′'s means in log1p CP50K and converted back the way config 1′'s Δ values are.
- **(b) cut points** (absolute):
  - half 0: 3.0 / 74.2;
  - half 1: 4.0 / 58.3;
  - all of H1: 3.0 / 67.7 (used for the full-panel surface and the guardrail grid).
- **H1 (b) terciles** (weak / mid / strong) under all-H1 cuts: 103 / 97 / 92. A knockdown with no b goes to weak.

**Fit** (on one H1 half: the highest overall among rows with all six components > 0; 38 and 32 of the 40 rows are eligible):

| Fit half → held out | Pick (a_weak, a_mid, a_strong) | Grid positions | Held-out overall vs base | Gain: overall / MSE / **PDS** / NMAE / fidelity / reach / Jaccard | K562 / HCT116 / HEK293T vs base (with the fit half's cuts) |
|---|---|---|---|---|---|
| 0 → 1 | **−0.25, 0, 0.5** | 1/2 (**edge**), 1/4 (**edge**), 3/5 (interior) | 0.2293 vs 0.2259 | **+0.0034** / +0.0054 / **+0.0078** / +0.0062 / 0 / +0.0012 / 0 | −0.0012 / +0.0069 / +0.0067 |
| 1 → 0 | **−0.25, 0, 0.25** | 1/2 (**edge**), 1/4 (**edge**), 2/5 (interior) | 0.2738 vs 0.2678 | **+0.0060** / +0.0029 / **+0.0304** / +0.0042 / 0 / −0.0017 / 0 | −0.0000 / +0.0001 / +0.0078 |

- Every held-out component is above 0, which is the criterion's test. Fidelity and reach move by at most −0.002.
- **PDS rises in both directions** (and on the guardrails): +0.002 to +0.050.
- Reach moves by ±0.002, so the proxy's reach bias (0.02–0.05 above cells, 15c §4) plays almost no part here.

**Criterion** (≥ +0.005 held out on both halves, every component > 0, no guardrail more than 0.005 worse): **not met.** The gain on the fit-0 direction is +0.0034; everything else passes.

**H1 full-panel surface** (all-H1 cuts; base 0.2469):
- **Best row:** (−0.25, 0, 0.25) at 0.2518 (+0.0049), which is also the best with all six components > 0.
- **Where the gain comes from:** almost all of it is a_weak.
  - (−0.25, 0, 0) alone gives 0.2509 (+0.0040), with PDS 0.803 → 0.829.
  - a_strong adds +0.001 at 0.25 and nothing at 0.5.
  - **A shift in the direction S points, on mid or strong knockdowns, does not help.**
- **Marginal means of overall:**
  - a_weak −0.25 / 0: 0.2379 / 0.2339;
  - a_mid 0 / 0.25 / 0.5 / 1: 0.2425 / 0.2399 / 0.2359 / 0.2253;
  - a_strong 0 / 0.25 / 0.5 / 1 / 1.5: 0.2430 / 0.2430 / 0.2415 / 0.2326 / 0.2194.
- **Large a_mid and a_strong raise MSE but cost more PDS:**
  - mean MSE falls 0.118 → 0.003 over a_strong 0 → 1.5;
  - mean PDS falls 0.784 → 0.750.
- **The guardrail grid** (all-H1 cuts) is the same shape. Best rows:
  - K562 (0, 0, 0.25) +0.0006;
  - HCT116 (0, 0, 0.5) +0.0097;
  - HEK293T (−0.25, 0, 0.5) +0.0069.
  - Every guardrail loses with a_mid > 0.

**Reference: the follow-up's H1-own shift, on the proxy** (the fit half's true-tercile mean residual, mid and strong, assigned by (b) tercile; the same halves):
- fit 0 → half 1: 0.2238 vs 0.2259 (**−0.0021**); MSE 0.109 → 0.166, but **PDS 0.750 → 0.563**.
- fit 1 → half 0: 0.2614 vs 0.2678 (**−0.0064**); MSE 0.144 → 0.194, but **PDS 0.857 → 0.667**.
- Its MSE gains (as on cells: 0.108 → 0.167, 0.144 → 0.198) are real, but **a shared shift on many knockdowns costs PDS more than it gains MSE.** This is the risk the section 11 follow-up named.

**How S relates to H1:**
- **Cosine with H1's strong-tercile mean residual** (truth − prediction), over H1's 10,746 gated genes where S is defined: **0.30**, against both the proxy base and the uploaded cells.
- **S's extremes are lineage-specific genes,** each gated in only one or two lines, which therefore carry the whole mean there:

| S (all lines, as for D/E/F) | Top 30 genes |
|---|---|
| **Up** | IGFBP1, SAA4, ANGPTL8, TM4SF5, SPINK1, PRAP1, EGR3, GZMA, VTN, ANKRD65, ITIH3, HCST, SLA, F10, TTR, HEPN1, TRAC, IL1B, ALB, CBLN4, KISS1, NR0B2, SCN3A, ALAS2, CLSTN2, UGT2B10, RRAD, HBA2, CCL20, GDNF |
| **Down** | TMEM151A, IGLL1, FABP1, CTSG, CAPS, CXCL9, HAMP, PGC, UPK3A, MTTP, GIMAP8, OVCH1, TMEM255A, GBP5, CHRD, OAS2, CDH19, PTX3, CPB2, RBP5, CCN3, TMEM253, FOXC2, TRPA1, C20orf204, GBP4, TSHR, PIGR, TNFRSF4, SLC5A11 |

  - Liver (HepG2): ALB, TTR, HAMP, FABP1, MTTP.
  - T cell (Jurkat): GZMA, TRAC, HCST.
  - Erythroid (K562): HBA2, ALAS2.
  - Interferon (Jiang).
  - **The rule "mean over lines that measure the gene" gives a gene one line's value wherever only that line gates it.**
  - The H1-excluded S has the same top genes.
  - **I would have restricted S to genes gated in at least half the lines.** I report it as specified.

## Part D — H1-derived shift with pluripotency removed (`tables/partD.json`, `partD_rows.csv`, `partD_S_H1_clean.npz`)
**This part can't be validated on an Arc line other than H1: none has truth.**

**Vectors:**
- **S_H1:** H1's strong-tercile minus weak-tercile mean residual, truth − config 1′ (proxy base, no-source fix), over all 292 knockdowns.
- **The pluripotency PC1:** the follow-up's residual programs reproduce exactly (each half's top five PCs, top 15 genes on both sides). I project out one vector: PC1 of the uploaded cells' mid + strong residual over the whole panel, whose |cos| with the two halves' PC1s is 0.93 and 0.87. Its genes:
  - +: H2AFY2, H1F0, CDH2, PHC2, LEFTY2, CER1;
  - −: GAL, MT1X, DNMT3B, DPPA2, COL1A1, SLC7A5.
- **PC1 holds 71% of S_H1's energy.**

**Energy shares** ((v·u)² / ‖v‖² for the follow-up's three generic programs: half 0 PC2, half 0 PC4, half 1 PC4):

| Vector | p53 / arrest vs S-phase | amino-acid stress | sterol synthesis | pluripotency PC1 |
|---|---|---|---|---|
| S_H1 | 0.319 | 0.001 | 0.112 | **0.713** |
| **S_H1,clean** | **0.493** | 0.000 | **0.281** | 0 |
| S, all lines | 0.014 | 0.000 | 0.011 | 0.009 |
| S without H1 | 0.010 | 0.000 | 0.009 | 0.004 |

**Cosines** over H1's gated genes:
- S_H1,clean with S: 0.31 (all lines), 0.27 (without H1).
- S_H1 with S: 0.40 / 0.31.
- **The three named programs hold < 3% of S's energy,** because S's energy sits on the lineage genes above.

**S_H1,clean, top 30:**
- **Up:** SERPINB9, LEFTY1, PLAAT3, H1F0, LEFTY2, NODAL, TMEM64, MMP24, C9orf135, PHLDA3, CDKN1A, SPP1, PLAAT5, GNA14, BBC3, PMAIP1, KLKB1, GRPR, RASGRF2, FN1, CHAC1, TAC1, PLK2, SESN1, KDR, AASS, HERC5, DHRS2, ZMAT3, PCDH18.
- **Down:** DGKK, DUSP6, NEFM, SFRP1, ID3, PCLAF, NPR1, AMOTL1, OPN3, SEMA5B, IDH2, PIR, DUT, RRM2, GNG4, SLC38A5, IGDCC3, KCNQ2, IVNS1ABP, HIST1H1A, CENPV, BRIP1, TMSB15A, LRP3, DERA, GINS2, LMO3, BTBD2, BTBD17, ARMCX2.
- **My reading:** p53 targets (CDKN1A, BBC3, PMAIP1, PHLDA3, ZMAT3, SESN1) up; S-phase (RRM2, PCLAF, DUT, GINS2, HIST1H1A) down.
- **Some lineage signal survives the projection:** LEFTY1/2, NODAL, DUSP6, SFRP1, ID3. One PC does not hold all of H1's differentiation response.

**S_H1,clean on the guardrails** (by (b) tercile, at Part C's picks, with each pick's fit-half cuts; vs base):

| Pick (from fit half) | K562 | HCT116 | HEK293T |
|---|---|---|---|
| (−0.25, 0, 0.5) (half 0) | +0.0005 (PDS +0.001) | **−0.0038** (PDS −0.017) | +0.0004 |
| (−0.25, 0, 0.25) (half 1) | +0.0017 (PDS −0.002) | −0.0008 (PDS −0.007) | +0.0001 |

On lines that aren't pluripotent, H1's strong-knockdown residual (pluripotency removed) neither helps nor clearly hurts. HCT116's PDS falls a little.

## Part F — program regression v2

`step20.py` was re-synced (md5-checked) twice before the parts that use it started.

### F1. Bases (`tables/partF1.json`, `partF1_<basis>_<fold>_basis.npz`, `_rows.parquet`)
**Weighting fix:** row weight 1/√max(n_b, 50), with all Tahoe tables as one background. A background with n_b ≥ 50 rows then carries total squared weight 1, and a smaller one carries n_b/50.

| Basis | Rows (tables) | Backgrounds < 50 rows: count, weight share | Cum. variance K 5 / 10 / 20 / 50 / 100 | K at 25% |
|---|---|---|---|---|
| **E2w** (E2's rows: n_real ≥ 50; conf.-weighted LFC) | 13,418 (107) | 18 (A549, HAP1, HCT116, HEK293T, HipSci pooled, 13 Feng donor clones), **0.417** | 0.089 / 0.129 / 0.172 / 0.237 / 0.305 | **59** (E2: 39) |
| **E4w** (E4's rows: all, full-profile tables; Δ) | 65,723 (113) | 2 (Feng fiaj_3, pipw_5), 0.011 | 0.072 / 0.096 / 0.126 / 0.168 / 0.209 | **not reached** (E4: 96) → 100 |
| **E24w** (n_real ≥ 50 rows; Δ entries) | 13,418 (107) | same 18 as E2w, 0.417 | 0.129 / 0.189 / 0.252 / 0.337 / 0.405 | **20** |

- **Fold bases** (rows / K at 25%):
  - **E2w:** −K562 11,897 / 56; −Orion 13,373 / 58.
  - **E4w:** −K562 56,156 / none; −Orion 33,293 / none.
  - **E24w:** −K562 11,897 / 19; −Orion 13,373 / 20.
- **E2's and E24w's rows are the same 13,418:** every row with n_real ≥ 50 comes from a table with full profiles.
- **Rows per background:**
  - **E2w / E24w:** Tahoe 6,588, K562 1,521, RPE1 1,383, KOLF2.1J 1,222, Jurkat 993, HepG2 837; the other 24 backgrounds hold ≤ 104 each.
  - **E4w:** HCT116 16,492, HEK293T 15,938, K562 9,567, Tahoe 6,808, then 1,600–2,050 for Jurkat, HipSci pooled, RPE1, HepG2 and KOLF2.1J.
- **The weight is spread evenly:** each background with ≥ 50 rows holds 1/20.5 = 0.049 of E2w's weight and 1/31 = 0.032 of E4w's.

**Did the weighting fix remove Step 19's small-table artefacts?** Per component: the share of row-score energy from backgrounds with < 50 rows, and how many of the 12 largest row scores they hold:

| Component | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| Step 19 E1 | 0.00/0 | **0.98/8** | **0.96/9** | 0.05/3 | **0.98/10** | 0.08/7 | 0.09/8 | **0.75/11** | **0.77/10** | **0.92/12** |
| Step 19 E4 | 0.00/0 | 0.09/10 | 0.04/3 | **0.91/8** | 0.34/11 | **0.71/10** | 0.11/7 | **0.77/12** | 0.23/10 | 0.30/8 |
| **E4w** | 0.00/0 | 0.00/1 | 0.01/0 | 0.00/0 | 0.01/1 | 0.00/0 | 0.00/0 | 0.01/0 | 0.00/0 | 0.01/0 |
| Step 19 E2 | 0.83/12 | 0.68/12 | 0.46/12 | 0.53/12 | 0.69/12 | 0.69/12 | 0.55/12 | 0.61/9 | 0.58/10 | 0.37/12 |
| E2w | 0.06/1 | 0.51/10 | 0.58/9 | 0.65/6 | 0.33/9 | 0.42/9 | 0.42/4 | 0.36/9 | 0.31/9 | 0.30/7 |
| E24w | 0.12/6 | 0.50/11 | 0.45/8 | 0.54/9 | 0.43/9 | 0.39/9 | 0.42/5 | 0.24/7 | 0.21/9 | 0.38/10 |

- **In E4w the artefacts are gone.** Each small background now carries weight in proportion to its rows, and only two of E4w's backgrounds are small.
- **In E2w and E24w they are reduced but not gone,** and they change character. No single 4-row table steers a component any more. Instead, the n_real ≥ 50 filter leaves 13 Feng donor clones (and HipSci pooled) with 7–70 rows each, and together they hold 42% of the weight. Components 2–6 are then Feng iPSC knockdowns: POU5F1, the PAF1 complex (PAF1, CTR9, RTF1), CNOT1, SRRT, SUPT5H.

**Top components, my reading (marked as such):**

*E4w (Δ):*

| PC | Var. | Top genes + / − | Largest rows | My reading |
|---|---|---|---|---|
| 1 | 0.027 | CDKN1A MDM2 ZMAT3 DDIT3 PHLDA3 RRM2B GADD45A / MKI67 TOP2A CENPF RRM2 HMGB2 CCNB1 | RPE1 ribosome-biogenesis knockdowns (NOL11, RPL14, WDR12) / Jiang receptor knockdowns | p53 arrest vs proliferation |
| 2 | 0.015 | NUP58 USP3 HBP1 ZNF292 TTK KDM6A / IFI6 ISG15 MIF MT2A HSPE1 NOP56 ATP5MC1 | Tahoe AURKB inhibitors / Tahoe RPS6KB1 inhibitors | long-transcript / nuclear vs translation-mitochondria (drug-driven) |
| 3 | 0.013 | BRIX1 NOP58 ODC1 SNRPD1 NME1 / SYNE2 UTRN DOCK1 RNF213 VPS13C | Tahoe AURKB / Feng pipw_5 SNRPN, oikd_2 UFD1 | MYC-like biogenesis vs long genes |
| 4 | 0.009 | NOP56 SLC7A5 HSPH1 ODC1 / MT1G MT1H TXNIP MT2A DDIT4 TOP2A | Feng CNOT1 / Jiang USP18, Feng tolg_4 | growth vs metallothionein / TXNIP stress |
| 5 | 0.007 | DDIT3 TXNIP TRIB3 IFRD1 GDF15 / HMGCS1 LDLR SCD IDI1 DHCR24 LEFTY1/2 | HepG2 EIF2S1 / EIF2S3 / EIF3 / Jiang USP18 | ISR (ATF4) vs sterol synthesis |
| 6 | 0.006 | MT2A CENPF PRPF3 INTS6 / GALK1 TP53I3 ACADVL FAH SAA4 | HepG2 splicing knockdowns / HepG2 RPL/RPS | splicing vs ribosome (HepG2) |
| 7 | 0.005 | POLR2A HSPA1A/B HMGCS1 IDI1 / MT1G MT1H IFIT2/3 DDX58 | Feng PAF1/CTR9 / Jiang USP18 | heat shock + sterol vs interferon |
| 8 | 0.005 | WARS PARP14 IFIT3 XAF1 IFIH1 / ZMAT3 CDKN1A MDM2 FDXR | Jiang IFNB USP18 / Jiang STAT1 / IFNGR1 | interferon vs p53 |
| 9 | 0.005 | SAA4 FOSL1 GDF15 CDKN1A TP53I3 / DDIT4 SLC7A11 PSAT1 PHGDH ATF4 PCK2 | HepG2 RPL / Jiang IFNGR1 | p53 vs amino-acid stress (ATF4) |
| 10 | 0.004 | LEFTY1/2 GDF15 TAGLN VIM SERPINE1 / STAT1 IFIT2/3 STAT2 | Jiang BXPC3 IFN pathway | mesenchymal / TGF-β vs interferon |

*E2w (LFC):*
- PC1 (0.025) is again p53 vs proliferation.
- PCs 2–4 (0.024 / 0.016 / 0.015) are pluripotency and the exit from it: POU5F1, GDF3, IDO1 vs LEFTY1/2, CER1, PKDCC. They are driven by Feng POU5F1 and PAF1-complex knockdowns.
- PC5 is interferon vs the integrated stress response (Jiang USP18 against Feng CNOT1).
- PC7 is interferon.
- PC8 is sterol synthesis (PCSK9, MVD, INSIG1, MSMO1).
- PC10 is sterol vs ISR.

*E24w (Δ, n_real ≥ 50):*
- PC1 (0.035) is p53 vs proliferation.
- PCs 2–6 are pluripotency, the exit from it, and Feng-specific programs (PAF1, POU5F1, CNOT1).
- PCs 7, 8 and 10 are interferon.

### F2. Ceiling on H1 (`tables/partF2_ceiling.csv`)
Perfect per-knockdown coefficients on the top K components, over each knockdown's gated genes. Captured energy is 1 − error/energy, in all knockdowns and by true-strength tercile.

**Scale:**
- The LFC bases report (a) LFC on gated genes, (b) LFC on true-DE genes, and (c) Δ, where the LFC fit is converted gene by gene.
- **The Δ bases fit Δ directly and report (c) only.**

| Basis, K | (a) LFC all / strong | (b) LFC all / **strong** | (c) Δ all / mid / **strong** |
|---|---|---|---|
| E2w, 5 / 10 / 20 / 50 / **59** (25%) | 0.054 / 0.070 · 0.084 / 0.108 · 0.132 / 0.168 · 0.179 / 0.224 · 0.185 / 0.231 | 0.092 / **0.098** · 0.138 / **0.147** · 0.218 / **0.228** · 0.290 / **0.301** · 0.299 / **0.310** | 0.072 / 0.025 / **0.088** · 0.109 / 0.041 / **0.132** · 0.174 / 0.084 / **0.207** · 0.236 / 0.128 / **0.278** · 0.243 / 0.134 / **0.286** |
| E4w, 5 / 10 / 20 / 50 / **100** (cap) | — | — | 0.061 / 0.026 / **0.073** · 0.104 / 0.054 / **0.123** · 0.198 / 0.119 / **0.231** · 0.278 / 0.172 / **0.323** · 0.322 / 0.204 / **0.373** |
| E24w, 5 / 10 / **20** (25%) / 50 | — | — | 0.075 / 0.029 / **0.091** · 0.116 / 0.047 / **0.140** · 0.184 / 0.092 / **0.219** · 0.253 / 0.145 / **0.296** |
| Step 19 E1, 5 / 10 / 20 / 50 / 100 | 0.019 / 0.024 · 0.039 / 0.049 · 0.057 / 0.072 · 0.170 / 0.212 · 0.197 / 0.243 | 0.030 / 0.032 · 0.061 / 0.064 · 0.092 / 0.097 · 0.274 / 0.283 · 0.313 / 0.323 | 0.022 / 0.009 / 0.027 · 0.054 / 0.023 / 0.065 · 0.076 / 0.032 / 0.091 · 0.221 / 0.129 / 0.259 · 0.252 / 0.149 / 0.295 |
| Step 19 E2, 5 / 10 / 20 / 39 (25%) | 0.050 / 0.065 · 0.082 / 0.106 · 0.117 / 0.149 · 0.147 / 0.185 | 0.087 / 0.092 · 0.137 / 0.145 · 0.195 / 0.204 · 0.242 / 0.252 | 0.064 / 0.023 / 0.078 · 0.109 / 0.041 / 0.133 · 0.157 / 0.075 / 0.188 · 0.194 / 0.097 / 0.231 |
| Step 19 E4, 5 / 10 / 20 / 50 / 96 (25%) | — | — | 0.051 / 0.021 / 0.061 · 0.074 / 0.037 / 0.087 · 0.116 / 0.063 / 0.137 · 0.254 / 0.155 / 0.295 · 0.319 / 0.202 / 0.369 |
| config 1′ (the base), for scale | −0.109 / −0.019 | 0.042 / 0.050 | −0.040 / −0.290 / 0.066 |

**What the ceiling changes:**
- **At K 20, the weighting fix more than doubles the strong-tercile ceiling** compared with Step 19's E1: E2w (b) 0.228 against 0.097, and E4w (c) 0.231 against 0.091. E4w also beats E4 (0.137).
- **This clears Step 19's ceiling criterion,** which E1 missed: E1 had less than 10% of strong-knockdown energy at K ≤ 20, while E2w, E4w and E24w are all above 20% at K 20.
- **It is still a ceiling with perfect coefficients.** F4 measures how much of it a prediction from the sources recovers.
- **Weak knockdowns stay at ≤ 4%** in every basis.

### F3. Gene features (`tables/partF3*.json`, `partF3_depmap.npz`, `partF3_string.npz`)
| Set | Release, URL | File, size, md5 | Built |
|---|---|---|---|
| **DepMap** | **DepMap 24Q4 Public** (figshare article 27993248, 2024-12-10), `https://ndownloader.figshare.com/files/51064667` | `CRISPRGeneEffect.csv`, 428,678,699 B, md5 `6edf7ade09b9b34199210b559d4745d3` (matches figshare's) | 17,915 genes × 1,178 lines; 32 PCs + mean gene effect |
| **STRING** | **v12.0**, `https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz` and `.../protein.info.v12.0/9606.protein.info.v12.0.txt.gz` | 83,164,437 B, md5 `5b7c91d02926e723ad30e81ab9d009a4`; 1,970,090 B, md5 `1ca9efcefcc68b5cc3b037a335cb870b` (STRING publishes no md5) | 16,040 genes, 234,368 edges with score ≥ 700; 32 eigenvectors + log(1 + degree) |
| **ARCHS4** | `human_gene_v2.latest.h5` (last modified 2026-07-07) | **62.3 GB** | **not downloaded** (62.3 GB). Nothing was downloaded, and M3c / M4c are not run. **Abandoned by Jose's decision (see the addendum at the end of this section)** |

**Two deviations, flagged:**
1. **DepMap isn't the latest release.** The portal's download API, which serves 26Q1, answers with a bot-verification page. On figshare, DepMap's newest release carrying `CRISPRGeneEffect.csv` is 24Q4. The 2026 figshare items are "Chronos parameters (Public 26Q1)": the raw Chronos `gene_effect.csv`, not the release file. The spec allows the portal or figshare, so I took the named file from figshare: two releases older, the same pipeline.
2. **STRING's spectral embedding is on the largest connected component** (15,726 of 16,040 genes; 125 components). On a disconnected graph, D^(−1/2) A D^(−1/2) has eigenvalue 1 once per component, so the whole graph's top 32 eigenvectors would just be indicators of small components. The other 314 genes keep log(1 + degree) and get the missing indicator. The adjacency is binary at score ≥ 700, one edge per gene pair (maximum score). The 32 eigenvalues run from 1.000 to 0.891.

**Mapping:** both sets are mapped to approved HGNC symbols. DepMap is mapped by Entrez ID, then by symbol; STRING by preferred name. Previous and alias symbols are used when they are unique. Targets are mapped the same way (H1F0 → H1-0, for example).

| | DepMap | STRING |
|---|---|---|
| Coverage: H1 / D/E/F / K562 / HCT116 / HEK293T targets | 289 of 292 / 300 of 300 / 267 / 278 / 271 (all) | 279 / 258 / 232 / 236 / 235 |
| Variance in the 32 components | 0.894 (**PC1 alone 0.856:** the common-essentiality axis) | eigenvalues 1.000 → 0.891 |
| Spearman of the extra feature with H1's true DE count / with (b) | mean gene effect **−0.47** / −0.46 (essential genes give stronger knockdowns) | log degree +0.28 / +0.29 |
| Within-complex cosine (10 HGNC complexes) vs random pairs | random 0.22. Exosome 0.97, NADH core 0.88, L / S ribosome 0.75 / 0.67, Integrator 0.70, CCT 0.63, Mediator 0.46, proteasome 0.36, COP9 0.25, **BAF 0.00** | random 0.35. Every complex 0.94–1.00 |
| Overlap: cosine of the two sets' gene-to-gene similarities (3,000 of 14,943 shared genes) | 0.33 (Pearson 0.22) | |

- **Both sets pass the complex check.**
- **STRING passes it by construction:** complexes are dense cliques in STRING.
- **DepMap's cosines are dominated by the essentiality axis.** That inflates random pairs (0.22), and the non-essential BAF complex scores 0.
- **The two sets overlap only weakly.**

### F3d. Do the embeddings carry response information? (no fitting; `tables/partF3d.json`)
**Method:**
- Pairs are every two knockdowns on a truth line whose targets both have features.
- Statistic: the Spearman between the pair's embedding cosine (32 components) and the cosine of the two knockdowns' true Δ (the line's gated genes, both targets out).
- **Null:** the gene-to-embedding mapping shuffled across the line's knockdowns, 100 times.
- **Baseline:** the same Spearman with the cosine of the two targets' source scores (E4w, 100 components; backgrounds averaged; the fold's own lines out), which is what M1 / M2 use.

| Line | Source scores (baseline): all / mid+strong | DepMap: all (null pct) / mid+strong (pct) | STRING: all (pct) / mid+strong (pct) |
|---|---|---|---|
| H1 | 0.042 / **0.083** | **0.025** (100th) / 0.028 (100th) | **0.024** (95th) / **0.065** (100th) |
| K562final | 0.008 / 0.008 | 0.014 (92nd) / 0.001 (60th) | 0.024 (97th) / 0.031 (96th) |
| HCT116final | 0.010 / −0.001 | 0.006 (83rd) / 0.004 (73rd) | 0.007 (81st) / 0.001 (61st) |
| HEK293Tfinal | 0.008 / −0.001 | −0.002 (44th) / 0.001 (58th) | 0.012 (82nd) / 0.009 (84th) |

- **On H1 both feature sets carry response information beyond chance,** but very little. The Spearman is 0.02–0.07 over ~18–40k pairs, about half of what the source scores carry.
- **On the guardrails** (sparser truth, fewer cells) neither feature set nor the source scores rise clearly above noise.

### F4. Coefficient prediction (`tables/partF4.json`, `partF4_rows.csv`, `partF4_preds.pkl`)
**Setup:**
- **Truth:** H1's program scores (LFC for LFC bases, Δ for Δ bases; least squares over each knockdown's gated genes).
- **Source scores ŝ:** every genetic row of every table (Tahoe drugs and Arc H1 out), unweighted and centred per table over all its rows. Each source line's tables are averaged first, then the result is projected. The fold's own lines are out.
- **R² convention:** held-out R² per program is 1 − SSE / Σy², relative to predicting 0 (Step 19's convention). An intercept alone therefore scores above 0, which is why the M4 null sits at ~0.02–0.07.
- **Ridge:** standardised features, intercept, and α by 5-fold CV within the training set.
- **M2 / M3 features** span 30–33 source backgrounds (33 in LFC space, 32 in Δ; minus the fold's own lines), each with its score and a missing indicator.
- **Pooled training (M2p, M3p):** each line's knockdowns are weighted by their signal fraction max(0, D − C)/D, rescaled so the lines contribute equally (mean weight 1, so α means the same with and without weights).

**Signal weights** (quantiles 10 / 50 / 90%; knockdowns above 0.1):

| Line | 10% / 50% / 90% | Knockdowns above 0.1 |
|---|---|---|
| H1 | 0.19 / 0.52 / 0.89 | 288 of 292 |
| K562final | 0.02 / 0.09 / 0.27 | 126 of 267 |
| HCT116final | 0.00 / 0.05 / 0.20 | 88 of 278 |
| HEK293Tfinal | 0.02 / 0.11 / 0.37 | 143 of 271 |

**α's grid position:** **76% of program fits choose α = 1,000, the top edge of the grid**, and the rest choose 100. No fit chooses 0.1–10.
- M2: 55% at the edge.
- M4: 73–76%.
- Pooled M3p: 98%.

The data want more shrinkage than the grid allows, so these fits are regularisation-bound. **A grid reaching 10⁴–10⁵ is the first thing I'd add.**

**Held-out (cross-fitted H1 halves) median R² over programs, and the captured energy of P_prog alone in the strong tercile (c, Δ scale):**

| Basis | K | M1 (Step 19) | M2 | M3 DepMap / STRING | M4 DepMap / STRING | **M4 null** (20 perms) | M2p | M3p DepMap / STRING |
|---|---|---|---|---|---|---|---|---|
| E2w | 10 | 0.053 · 0.024 | 0.196 · 0.041 | 0.241 · 0.042 / 0.207 · 0.041 | 0.108 · 0.019 / 0.107 · 0.018 | 0.047 · 0.013 / 0.043 · 0.012 | 0.184 · 0.044 | 0.242 · 0.046 / 0.234 · 0.048 |
| E2w | 20 | 0.054 · 0.038 | 0.126 · 0.058 | 0.210 · 0.060 / 0.163 · 0.057 | 0.092 · 0.028 / 0.054 · 0.021 | 0.029 · 0.015 / 0.023 · 0.015 | 0.137 · 0.059 | 0.233 · 0.065 / 0.189 · 0.065 |
| E2w | 59 | 0.055 · 0.050 | 0.117 · 0.061 | 0.139 · 0.074 / 0.122 · 0.063 | 0.074 · 0.037 / 0.048 · 0.028 | 0.023 · 0.021 | 0.112 · 0.038 | 0.127 · 0.037 / 0.127 · −0.057 |
| E4w | 10 | 0.053 · 0.015 | 0.167 · 0.037 | 0.152 · 0.035 / 0.182 · 0.035 | 0.110 · 0.020 / 0.090 · 0.018 | 0.048 · 0.014 / 0.041 · 0.014 | 0.249 · 0.040 | 0.259 · 0.042 / 0.256 · 0.042 |
| E4w | 20 | 0.070 · 0.035 | 0.174 · 0.070 | 0.205 · 0.072 / 0.201 · 0.069 | 0.119 · 0.035 / 0.070 · 0.030 | 0.043 · 0.022 / 0.036 · 0.021 | 0.219 · 0.069 | 0.261 · **0.077** / 0.242 · 0.073 |
| E4w | 100 | 0.074 · 0.058 | 0.124 · 0.063 | 0.140 · **0.080** / 0.127 · 0.061 | 0.063 · 0.051 / 0.043 · 0.039 | 0.018 · 0.031 | 0.126 · **−0.284** | 0.169 · −0.104 / 0.148 · −0.501 |
| E24w | 10 | 0.042 · 0.016 | 0.195 · 0.039 | 0.220 · 0.042 / 0.187 · 0.041 | 0.098 · 0.019 / 0.074 · 0.019 | 0.045 · 0.013 / 0.041 · 0.013 | 0.248 · 0.044 | **0.300** · 0.048 / 0.279 · 0.049 |
| E24w | 20 | 0.070 · 0.031 | 0.184 · 0.063 | 0.215 · 0.065 / 0.177 · 0.065 | 0.088 · 0.029 / 0.056 · 0.023 | 0.032 · 0.016 / 0.030 · 0.016 | 0.220 · 0.065 | 0.259 · 0.073 / 0.248 · 0.072 |
| Step 19 E2 | 20 | 0.059 · 0.030 | 0.141 · 0.049 | 0.190 · 0.054 / 0.136 · 0.053 | 0.077 · 0.024 / 0.056 · 0.020 | 0.028 · 0.015 | 0.084 · 0.046 | 0.133 · 0.054 / 0.114 · 0.055 |
| Step 19 E4 | 20 | 0.058 · 0.017 | 0.153 · **−0.146** | 0.196 · −0.150 / 0.182 · −0.147 | 0.135 · 0.023 / 0.099 · 0.021 | 0.034 · 0.015 | 0.210 · −0.197 | 0.246 · −0.568 / 0.247 · −0.405 |

Each cell is median R² · c_strong.

**Reference, the same strong-tercile measure (c):**

| | K 20 | K 100 |
|---|---|---|
| Step 19's E-reg (E1, M1 with one averaged vote) | 0.013 | 0.048 |
| Config 1′ itself | 0.066 (overall) | |
| Perfect-coefficient ceilings: E2w / E4w / E24w | 0.207 / 0.231 / 0.219 | E4w 0.373 |

**What each ingredient adds** (E4w and E24w, K 10–20):
1. **Ridge on separate source lines (M2) against one averaged vote (M1):** median R² 0.05–0.07 → 0.17–0.20; strong-tercile capture ×2 (0.035 → 0.070 at E4w K 20). **The biggest single step.**
2. **Gene features on top (M3):**
   - R² +0.01–0.04, capture +0.000–0.003.
   - DepMap slightly ahead of STRING at K 20; STRING ahead at K 10 on E4w.
   - The inner-CV choice picks DepMap in most directions.
3. **Gene features alone (M4) against the permutation null:**
   - R² 0.06–0.12 against 0.03–0.05;
   - capture 0.018–0.035 against 0.013–0.022.
   - **Gene function alone predicts some of the response, beyond chance, at about half of what the source lines do.**
4. **Pooled truth lines (M2p, M3p):**
   - At K ≤ 20: R² +0.03–0.08 over H1-only training, capture +0.000–0.008.
   - **At high K pooling breaks** (E4w K 100: −0.28 to −0.50). The guardrails' labels on low-variance components are noise.
5. **Share of the ceiling recovered** (strong, c):
   - M1: 0.11–0.18;
   - M2 / M3: 0.27–0.31;
   - M3p: 0.30–0.35 at K ≤ 20.

   **A third of the perfect-coefficient ceiling is the most any method gets.**
6. **Step 19's E4 with source scores fails:**
   - The median R² is positive, but a few weak knockdowns get huge predictions (weak-tercile c −178).
   - Their source scores are extreme on E4's small-table artefact components (PCs 4, 6, 8), and ridge extrapolates from them.
   - M1 and M4 are unaffected.
   - The weighting fix removes the problem in E4w.

### F5. Proxy arms (`tables/partF5.json`, `partF5_rows.csv`, `partF5_DEF_Pprog_*.npz`)
**Setup:**
- **Calls:** K0. **Base:** config 1′ + K0 + the no-source fix.
- **Grid:** bases E2w / E4w / E24w; methods M1, M2, M3, M2p, M3p; K ∈ {10, 20, K at 25%}; arms R2 (all knockdowns) and R3 (b ≥ 5 only); γ ∈ {0.25, 0.5, 0.75, 1}.
- **How P_prog is added:** in Δ space for E4w / E24w, converted like config 1′'s Δ values. For E2w it is added as an LFC.
- **Rows:** H1 has 320 rows per fit direction, not 360, because E24w's K at 25% is 20, the same as K 20. K562 has 360 (its fold basis reaches 25% at K 19); HCT116 and HEK293T have 320.
- **Which M3 / M3p variant:** the inner-CV choice on the training set picked DepMap in 24 of 32 H1 slots. STRING won the other 8, all in the fit-1 direction at K 10–20. On the guardrails DepMap won 48 of 50 (STRING: E4w K 10 for HCT116 and HEK293T).
- **Run:** two 4-thread workers took 43 min (H1 3.4 s a row, guardrails 2.2–2.9 s).

**Fit** (on one H1 half: the highest overall among rows with all six components > 0 there; 320 of 320 rows eligible in each direction):

| Fit half → held out | Pick | Grid positions: basis / method / K / arm / γ | Held-out overall vs base | Gain: overall / MSE / **PDS** / NMAE / fidelity / reach / Jaccard | Reach share of component gain |
|---|---|---|---|---|---|
| 0 → 1 | **E4w, M3 (DepMap), K 100, R3, γ 1** | 2/3 · 3/5 · **3/3 edge** · **2/2 edge** · **4/4 edge** | 0.2510 vs 0.2259 | **+0.0251** / +0.071 / **−0.061** / +0.086 / +0.001 / +0.054 / 0 | 0.36 |
| 1 → 0 | **E4w, M3p (DepMap), K 100, R3, γ 1** | 2/3 · **5/5 edge** · **3/3 edge** · **2/2 edge** · **4/4 edge** | 0.2482 vs 0.2678 | **−0.0197** / **−0.144 (MSE scaled 0)** / **−0.051** / +0.024 / +0.001 / +0.053 / 0 | — |

**Guardrails at the picks** (the fold's basis, K 100):

| Pick | K562 | HCT116 | HEK293T |
|---|---|---|---|
| Fit 0 (M3) | +0.0068 (PDS −0.036) | **+0.0107** (PDS **+0.014**) | +0.0108 (PDS +0.030) |
| Fit 1 (M3p) | +0.0050 (PDS −0.052) | +0.0050 (PDS **−0.037**) | +0.0128 (PDS +0.023) |

**Criterion** (≥ +0.005 held out on both halves, every component > 0, no guardrail more than 0.005 worse): **not met.**
- The fit-1 direction picks the pooled M3p at K 100, which is in-sample on its own half. It fails held out, with MSE collapsing to 0 scaled, as F4's cross-fitted M3p at K 100 predicted (strong-tercile c −0.10).
- The guardrails pass in both directions.

**Why the procedure picked a failing row:** the fit-half rows use the model trained on that half, so they are in-sample. Ridge with ~70–100 features (and pooled lines) overfits its own half far more than Step 19's one-slope M1 did. The spec's fit rule therefore favours the most flexible method. **Cross-fitted fit-half rows** (inner folds within the fit half) **would have been the right selection rule.** I report the procedure as specified.

**The same configuration for every method** (E4w, K 100, R3, γ 1; not a selection, all reported):

| Method | Held-out gain, fit 0 → 1 / fit 1 → 0 | PDS change on H1 | K562 / HCT116 / HEK293T |
|---|---|---|---|
| M1 (Step 19's slope on the averaged vote) | +0.011 / +0.018 | −0.026 / −0.004 | +0.006 / −0.003 / +0.006 |
| M2 (ridge on source lines) | +0.020 / +0.008 | −0.060 / −0.075 | +0.001 / −0.001 / +0.014 |
| **M3 (ridge + DepMap)** | **+0.025 / +0.023** | −0.061 / −0.062 | **+0.007 / +0.011 / +0.011** |
| M2p (pooled) | +0.023 / **−0.031** | −0.051 / −0.060 | −0.001 / −0.004 / +0.008 |
| M3p (pooled + DepMap) | +0.028 / **−0.020** | −0.052 / −0.051 | +0.005 / +0.005 / +0.013 |

**What each ingredient adds:**
- **The arm matters most.**
  - R3 (only knockdowns with b ≥ 5) beats R2 (all knockdowns) in almost every cell.
  - At E4w K 100 γ 0.5, R2 loses held out with every method except M1 and pooled M3p in one direction (M2 −0.002 / −0.023), and costs 0.15–0.25 in PDS.
  - **Adding programs to weak knockdowns makes them look alike.**
- **The new basis and the R3 arm alone, with Step 19's M1,** give +0.011 / +0.018. That's above Step 19's E-reg pick (+0.009 / +0.005).
- **Ridge on separate source lines plus DepMap (M3)** adds about +0.01 on top. It is the only method positive in both directions and on every guardrail at this configuration.
- **Pooling truth lines (M2p, M3p) helps at K ≤ 20** (E4w K 20 R3 γ 1: M3p +0.017 / +0.016) and breaks at K 100.

**Best held-out row per method** (any configuration; selection on the held-out half, so descriptive only):

| Method | Fit 0 → 1 | Fit 1 → 0 |
|---|---|---|
| M1 | +0.012 | +0.018 |
| M2 | +0.021 | +0.014 |
| M3 | +0.025 | +0.023 |
| M2p | +0.023 | +0.014 (E4w K 20) |
| M3p | +0.028 | +0.017 (E24w K 20) |

All are at R3; most are at E4w K 100.

**For comparison, Step 19's E-reg pick on the same halves:** R2 γ 0.5 / 0.75, E1 K 100: +0.009 / +0.005, PDS −0.02 / −0.04, and HCT116 −0.012 / −0.017.

**F4's held-out captured energy for every M3 / M4 variant and the null** is in F4's table above:
- DepMap and STRING are within 0.01 of each other everywhere;
- both M4s are 1.5–3× their permutation null;
- ARCHS4 was not built.

**Proxy caveat:** the proxy overstates reach by 0.02–0.05 compared with cells (15c §4). In the fit-0 pick, reach is +0.054 of the components' +0.150 sum (36%). If cells showed no reach gain at all, the held-out gain would shrink from +0.025 to about +0.016.

**D/E/F:** P_prog for the 300 targets in each context, for both picks. All sources are used, with Arc H1 as one more background. The fit uses all of H1 (the pooled pick: H1 + the three guardrail lines, with Arc H1's column available to the guardrail rows). The scores are masked to each context's > 5 CPM gate.
- Files: `tables/partF5_DEF_Pprog_E4w_M3DepMap_K100.npz` and `..._M3pDepMap_K100.npz`.
- 298 of 300 targets have scores, with a median of 3 backgrounds.
- **Nothing was built from them.**
- (For M3, an H1-only fit gives the Arc H1 column no weight, because H1's own rows can't use it.)

### F6. Co-expression prior from each line's own controls (`tables/partF6.json`, `partF6_rows.csv`, `partF6_C_{D,E,F}.npz`, `partF6_prior_<line>.npz`)
**Build:**
- Metacells: k = min(600, cells / 30). D, E, F, H1 and K562 have 18,400 cells, so 600 metacells. HCT116 and HEK293T have ~9,940 cells, so 331.
- Gated genes: D / E / F 9,633 / 9,743 / 9,537; H1 10,778.
- **Sanity check on D / E / F:** mitochondrial genes correlate 0.55 / 0.43 / 0.17 with each other, against 0.01 for random gene pairs. So C picks up modules.

**F6a (no fitting):** the per-knockdown Pearson of true Δ with −C_k. Null: C's row for a random gated gene, 100 draws per knockdown.

| Line | Knockdowns with a prior | Median r, all (null; percentile) | Weak / mid / **strong** | True DE genes in the target's top 100 co-expressed (null) |
|---|---|---|---|---|
| H1 | 289 / 292 | 0.003 (0.000; 78th) | −0.003 / −0.003 / **0.017** (98th) | 0.024 (0.013) |
| K562final | 229 / 267 | −0.006 (0.000; 3rd) | 0.004 / −0.008 / −0.009 | 0.011 (0.015) |
| HCT116final | 270 / 278 | −0.012 (0.002; 6th) | −0.023 / −0.011 / −0.012 | 0.010 (0.011) |
| HEK293Tfinal | 261 / 271 | −0.000 (0.004; 34th) | 0.010 / −0.017 / 0.008 | 0.011 (0.014) |

- **A knockdown's response is unrelated to its target's co-expression partners in the same line's controls.**
- The only signal above the null is H1's strong tercile: r 0.017, with DE genes enriched about 2× in the top 100.
- On the three guardrails the correlation is at or below the null.

**F6b (proxy):**
- λ_t, fitted per (b) tercile through the origin on the fit half, is tiny and changes sign between halves:
  - fit half 0: −0.004 / +0.013 / +0.003;
  - fit half 1: −0.001 / −0.011 / +0.023;
  - all of H1: −0.003 / +0.004 / +0.009.
- **Both directions pick γ 0.5** (2/4, interior; the grid includes γ 0 = base):

| Fit half → held out | Held-out overall vs base | Gain: overall / MSE / **PDS** / NMAE / fidelity / reach / Jaccard | K562 / HCT116 / HEK293T vs base |
|---|---|---|---|
| 0 → 1 | 0.2263 vs 0.2259 | **+0.0004** / −0.0001 / **−0.0022** / +0.0026 / 0 / +0.0020 / 0 | +0.0002 / −0.0043 / −0.0006 |
| 1 → 0 | 0.2671 vs 0.2678 | **−0.0007** / −0.0003 / **−0.0021** / +0.0011 / 0 / −0.0030 / 0 | +0.0002 / −0.0043 / −0.0006 |

**Criterion** (≥ +0.005 held out on both halves, every component > 0, no guardrail more than 0.005 worse): **not met.** The gains are about 0.

**D/E/F:** C (float16) and −C_k for the 300 targets are saved in `partF6_C_{D,E,F}.npz`. 258–300 of the targets have a prior, depending on the context. Nothing was built from them.
## Part G — recommendation (`tables/partG_rows.csv`, `partG_specs.json`)
**One spec, the same for D, E and F.**
- Step 20 found nothing that separates the contexts.
- The parts that depend on context (each context's gate, its own controls) enter through the existing pipeline.

| Item | D | E | F | Evidence |
|---|---|---|---|---|
| **Calls** | **hybrid T 5** (K4 1,200 where b ≥ 5, else K0) | same | same | Step 19 Part A (cells). Parts A and B give no reason to move to K0 or K4 everywhere |
| **Shift** | **none** | none | none | Part C misses its criterion (+0.003 / +0.006). Part D can't be validated and does nothing on the guardrails |
| **Arc H1 weight** (its 28 targets) | **1** | 1 | 1 | Part B: upweighting a same-lab source lowers overall and PDS on both Orion analogues |
| **Program regression** | **none in the file as built; candidate for a cell check:** E4w, M3 (DepMap), K 100, R3, γ 1 | same | same | F5 misses its criterion (fit-1 direction). The fit-0 pick is positive held out and on every guardrail, and this configuration is positive in both directions for M3 |
| **Own-control co-expression prior** | **none** | none | none | F6: no relation to responses; F6b gain ≈ 0 |
| **No-source fix** | **included** | included | included | 15c §4; Step 19 Part F |

**Why calls stay at hybrid T 5:**
- Part A doesn't show D/E/F's truth will be more H1-like than any other line's: H1 transfers to Weissman's lines about as well as they transfer among themselves.
- Part B shows a same-lab source is worth a lot as an ordinary voter, which says nothing about call density.
- The case for denser calls on D/E/F is still Step 19 Part B's: Arc's depth (400 cells, ~20k UMI) is the same as H1's. That case was already behind T 5.
- **K4 everywhere** gains +0.0012 more on Eval A cells but costs HCT116final twice as much (Step 19 Part A).
- **The K4 rules risk** (15c §3) **remains your call.**

**Combination rows on the proxy** (H1 full / half 0 / half 1; fold pools full):

| Row | H1 full | H1 halves | K562 | HCT116 (fold) | HEK293T | H1 PDS |
|---|---|---|---|---|---|---|
| live entry (config 1′ + K0, no fix; Step 19 Part F) | 0.2456 | 0.2664 / 0.2247 | 0.1582 | 0.3108 | 0.1536 | 0.795 |
| no-source fix alone | 0.2469 | 0.2678 / 0.2259 | 0.1582 | 0.3108 | 0.1536 | 0.803 |
| hybrid T 5 alone (Step 19) | 0.2669 | 0.2878 / 0.2458 | 0.1471 | 0.3017 | 0.1530 | |
| **recommended spec: hybrid T 5 + fix** | **0.2682** | 0.2892 / 0.2471 | 0.1471 | 0.3017 | 0.1530 | 0.803 |
| fix + Part C shift (−0.25, 0, 0.25) | 0.2518 | 0.2737 / 0.2297 | 0.1578 | 0.3127 | 0.1597 | 0.824 |
| fix + F5 term (E4w M3 DepMap K 100 R3 γ 1) | 0.2709 | 0.2907 / 0.2510 | 0.1650 | 0.3215 | 0.1644 | 0.742 |
| spec + shift | 0.2728 | 0.2948 / 0.2507 | 0.1458 | 0.3042 | 0.1591 | 0.824 |
| **spec + F5 term** | **0.2893** | 0.3095 / 0.2690 | **0.1517** | **0.3147** | **0.1651** | 0.742 |
| spec + shift + F5 term | 0.2912 | 0.3127 / 0.2695 | 0.1505 | 0.3143 | 0.1695 | 0.764 |

- **On H1 the F5 term is cross-fitted:** each half's knockdowns get the model trained on the other half. Its configuration was chosen on half 0, so half 1 (0.2510 / 0.2690) is the clean held-out number.
- **The parts nearly add:**
  - spec +0.0213 and F5 term +0.0240 from the fix-only base give +0.0424 together, against a sum of +0.0453;
  - the shift adds +0.002–0.005 on top.
- **What the F5 term does on top of the spec:**
  - it brings the guardrail losses from the hybrid calls (K562 −0.011, HCT116 −0.009 against the fold pools) back to K562 −0.007 and HCT116 +0.004, and HEK293T to +0.011;
  - **it moves H1's MSE 0.127 → 0.196 and NMAE 0.008 → 0.092, and costs PDS 0.803 → 0.742.**

**Expected change against the live entry** (Eval A 0.2426, HCT116final 0.3617 on cells):

| Change | Eval A | HCT116final | Source |
|---|---|---|---|
| hybrid T 5 calls | +0.0201 [+0.0155, +0.0251] | −0.0070 [−0.0124, −0.0017] | **cells** (Step 19 Part A) |
| no-source fix | +0.0013 | +0.0005 | proxy (generator input) |
| **recommended spec, total** | **≈ +0.021** (→ ≈ 0.264) | **≈ −0.007** (→ ≈ 0.355) | cells + proxy |
| *if the F5 term is added after a cell check* | *+0.021 more on the proxy (H1 0.2682 → 0.2893); ~+0.013 if its reach gain doesn't hold on cells* | *proxy only on the HCT116 fold pool: +0.013 over the spec; HCT116final's generator input not scored* | proxy |

**Is a cell check needed?**
- **For the recommended spec: no.** The calls are a cell result, and the fix is small and was tested on the proxy (Step 19).
  - **Build:** ≈ 2 h (Step 18's route).
- **For the F5 term: yes, before it goes in a file.**
  - It changes the means of the knockdowns with b ≥ 5 on every gated gene.
  - It may change which uncalled genes turn significant, which the proxy tracks least well.
  - A third of its proxy gain is reach.
  - **Then:** a D/E/F build (~2 h), from `partF5_DEF_Pprog_E4w_M3DepMap_K100.npz`.
  - **Timing:** both would fit before the 5 Nov deadline (today is 2 Nov), with room for one upload of each file. **Whether to spend that time is your call.**

## What I would have done differently
2. **Checked Part B's K562 analogue before running it.** RPE1, HepG2 and Jurkat are essential-gene screens and hold none of K562final's 267 targets. A K562 analogue needs targets those screens hold.
3. **Raised F5's selection rule before the run.** The fit-half rows are in-sample for each model. That was harmless for Step 19's one-slope M1, but it makes the rule pick the most flexible method (pooled ridge at K 100), and F4's cross-fitted table had already shown that row failing. Selecting on cross-fitted inner folds within the fit half would avoid it.
4. **Extended the ridge α grid.** 76% of the program fits chose α = 1,000, the grid's top edge.
5. **Flagged two things in Part C's S before running it:**
   - "One vote per line" gives Feng's iPSC lines 20 of 33 votes.
   - "Mean over the lines that measure the gene" lets a single line define lineage-specific genes, which then dominate S's extremes (ALB, TTR, GZMA, HBA2).
   - A minimum line count per gene (half the lines, say) and a vote per lab or lineage would have made S the generic response it was meant to be.
6. **Sized Part A's memory before writing it.** Float32 for 90k rows × 4 layers was never going to fit, and float16 was the plan from attempt 2 (16 s lost).
8. **Estimated from per-item timings, not whole-part analogies.** Part B ran 3× faster than estimated (0.1 s per knockdown) and F1 / F6 ran 3× faster, while F4 ran slower than my upper bound.

## Addendum (Jose, 2026-11-02): decisions after the write-up
**ARCHS4 abandoned (Jose's decision).** The 62.3 GB ARCHS4 file was not downloaded, and the ARCHS4 feature set is dropped, so M3c and M4c will not be run. The grounds:
- **The three sets are close kin.** ARCHS4 is another gene-level embedding of the knocked-down gene, built the same way as DepMap's and STRING's.
- **The two we built carried little.**
  - DepMap carried weak information:
    - knockdown strength: Spearman −0.47 with H1's true DE count, and −0.25 after allowing for (b);
    - response similarity on H1: Spearman 0.025, about half what the source scores carry;
    - about +0.01 of the best F5 term (M3 against M2 at E4w, K 100, R3, γ 1).
  - **STRING added nothing useful** beyond the complex structure it holds by construction.
- **The expected gain from a third set of this kind is very small,** and the time is better spent before the 5 Nov deadline.

**DepMap is accepted as run** (24Q4 release, via figshare). Jose is satisfied with it, so it is not re-fetched from the 26Q1 portal.

**Part C rerun with the fixed shift vector (Jose, 2026-11-02; `step20.py partC2`, `tables/partC2.json`, `partC2_rows.csv`, `partC2_S.npz`, `logs/partC2.log`, 09:10–09:31 UTC).**

*Why S was flawed:*
- **Single-line genes dominated it.**
  - "The mean over the lines that measure the gene" let a gene measured in one line take that line's value alone.
  - 68% of S's energy sat on 1,531 genes defined by only 1–2 lines (liver genes from HepG2, T-cell genes from Jurkat, globins from K562).
  - Part C added S on every gene where it was defined, and 82% of its energy was on genes H1 barely expresses.
- **Most of the votes were noise.**
  - One vote per line gave Feng's iPSC lines 20 of the 33 votes.
  - Their median `n_real` is 0–2, so their "strong minus weak" compares almost-null knockdowns.

*The fix* (rules set from the source tables alone, before any scoring):
1. **Keep a line only if** it has ≥ 30 knockdowns per `n_real` tercile and its strong tercile's median `n_real` is ≥ 30.
   - The filter is on the strong tercile, not the line's overall median, because K562's overall median is 1 while its strong third has 41.
   - **13 lines pass:** A549, BXPC3, HAP1, HT29, MCF7 (Jiang), K562, RPE1, HepG2, Jurkat (Weissman), KOLF2.1J, Feng jejf_3 and oikd_5, and H1.
   - **Dropped:** the other 18 Feng lines, both Orion lines (strong-tercile median 7–8) and the tiny tables.
   - H1 is left out when scoring H1, and the fold exclusion now changes nothing for HCT116 / HEK293T.
   - Jiang now holds 5 of the 12 votes that build H1's S.
2. **Define S only on genes that at least half of the voting lines measure:** 10,552–10,829 genes.
3. **Apply S only on the scored line's own gated genes.**

The grid, the (b) cut points and the fit rule are unchanged from Part C.

*What S looks like now.* Its extremes are a generic program:
- **up:** p53 / stress (CDKN1A, DDIT3, GDF15, ZMAT3, MDM2, TP53INP1, PHLDA3, GADD45A, BBC3, TP53I3, SQSTM1, TXNIP);
- **down:** proliferation (MKI67, TOP2A, CENPF, RRM2, PCLAF, MCM3/4/5/7, GINS2, CCNA2, CCNB1, CDC20, histones).

This is close to E4w's PC1. The cosine with H1's strong-tercile residual (H1's gated genes) is 0.29, against 0.30 before.

| | Part C (as specified) | **Part C, fixed** |
|---|---|---|
| Pick, fit half 0 → 1 | (−0.25, 0, 0.5) | (−0.25, 0, 0.5); grid 1/2 **edge**, 1/4 **edge**, 3/5 |
| held-out gain: overall / MSE / PDS / NMAE / reach | +0.0034 / +0.005 / +0.008 / +0.006 / +0.001 | **+0.0040** / +0.014 / −0.001 / +0.008 / +0.003 |
| guardrails K562 / HCT116 / HEK293T | −0.001 / +0.007 / +0.007 | +0.002 / +0.009 / +0.007 |
| Pick, fit half 1 → 0 | (−0.25, 0, 0.25) | (−0.25, 0, 0.25); grid 1/2 **edge**, 1/4 **edge**, 2/5 |
| held-out gain: overall / MSE / PDS / NMAE / reach | +0.0060 / +0.003 / +0.030 / +0.004 / −0.002 | **+0.0088** / +0.009 / +0.033 / +0.007 / +0.004 |
| guardrails | −0.000 / +0.000 / +0.008 | −0.000 / +0.002 / +0.010 |
| H1 full-panel best (base 0.2469) | (−0.25, 0, 0.25) 0.2518 | (−0.25, 0, 0.25) **0.2538** |
| a_weak −0.25 alone | 0.2509 (+0.0040; PDS +0.026) | 0.2528 (**+0.0059**; PDS +0.029) |
| a_strong alone: 0.25 / 0.5 / 1.0 | +0.0006 / −0.0008 / −0.0067 | +0.0006 / −0.0005 / −0.0068 |

**Criterion** (≥ +0.005 held out on both halves, every component > 0, no guardrail more than 0.005 worse): **still not met.**
- The gain in the fit-0 direction is +0.0040.
- Every held-out component is above 0, and no guardrail is worse than −0.0002.

**Reading (mine):**
- **The fix cleaned S and lifted every number a little,** but the conclusion stands: **shifting strong or mid knockdowns along S does not help.** a_strong ≥ 0.5 and every a_mid > 0 lower overall, because the PDS cost outweighs the MSE gain.
- **All the gain comes from a_weak = −0.25:**
  - it moves weak knockdowns' predictions against the generic p53-up / proliferation-down program;
  - +0.006 on H1, positive in both held-out directions, and +0.002 to +0.010 on the guardrails.
  - My reading is that config 1′'s weak-knockdown predictions inherit some of the pooled stress response from strong source knockdowns, which the truth for weak knockdowns doesn't show. Removing it makes weak knockdowns more distinct (PDS).
- **−0.25 is the edge of the grid.**
- **The next test, if wanted, is a weak-only arm:** a_weak ∈ {−0.25, −0.5, −0.75, −1}, mid and strong at 0, fitted the same way. It would take about 10 min on the proxy, and would need a cell check before going into a file.
- **Recommendation for the next file is unchanged:** no shift.

**Decision (Jose, 2026-11-02): take the fixed Part C shift's gain.** Jose is inclined to include it in the final spec despite the missed criterion. The weak-only arm (a_weak beyond −0.25) is **not** to be run.

**Submission checklist on the live entry's cells (Jose's eight scenarios; `code/vcc6/checklist20.py`, `results/stage1/2026-11-02_step20/checklist/checklist.json`).**
- **What was checked:** config 1′ + K0's generator outputs for D / E / F (Step 18, 300 × 400 cells each, with each context's 18,400 real controls), and the same config's Eval A cells (15c) against Eval A's real cells.
- **Panel note:** the panel has no RPL/RPS genes, so the ribosomal share is 0 on both sides.

| # | Scenario | Result | Verdict |
|---|---|---|---|
| 1 | Rounding to integers | Gated genes (> 5 CPM) whose knockdown mean is exactly 0: median 0 per knockdown, max 1–2, in 6–11 knockdowns per context. Eval A real cells: max 1, 2 knockdowns | pass |
| 2 | Unmeasured genes | Gated (knockdown, gene) pairs with no source: 0.8–1.0% of pairs. Their predicted mean vs control: median log2 −0.02 (D) / −0.09 (E) / −0.07 (F), 90th percentile of \|log2\| 0.13–0.19 (measured pairs: median \|log2\| 0.07–0.09). 1–4 gated genes have no source for any target | pass (E sits ~6% low) |
| 3 | Total-count distortion | Top-50 control genes' share of counts: 0.156 vs 0.155 (D), 0.149 vs 0.149 (E), 0.121 vs 0.123 (F), 0.0809 vs 0.0810 (Eval A; real knockdown cells 0.0806). MT share within 1% of controls everywhere. Per-knockdown top-50 share, 5–95%: ±3% | pass |
| 4 | Cells too uniform → too many calls | Eval A significant genes per knockdown (cell-eval's Wilcoxon): **ours median 483 / 90th percentile 1,082; real 219 / 3,335**. Spearman across knockdowns 0.09 | **partial fail:** total volume is near the real mean, but it is spread evenly. Weak knockdowns are over-called about 2×, strong ones under-called about 3×. This is the call allocation problem (Steps 15b, 17), not cell spread |
| 5 | Cells too noisy | Within-knockdown Fano factor at matched mean, ours vs reference (Eval A real knockdown cells, or the context's controls): low-expression genes (0.1–1 counts per cell) **0–26% too dispersed** (Eval A +6% / 0%; D/E/F +7–26%); high-expression genes (3–30 counts) **28–54% under-dispersed** (Eval A 1.49 vs 3.16 at 10–30 counts per cell; D 3.49 vs 6.72) | **mixed:** the low-gene excess fits Step 4's finding that 85% of Eval A knockdowns hit the noise-credit cap (Ĉ_pred > C_real), worth at most raw −0.024 MSE |
| 6 | Cells reused | 0 duplicate cells within or across knockdowns, 0 identical to a control | pass |
| 7 | Depth too low | Median total counts: 20,259 / 19,560 / 19,705 (D/E/F) vs controls 20,236 / 19,681 / 19,845; Eval A 19,991 vs real 19,662 | pass |
| 8 | Depth inconsistent | **The opposite problem: depth is far too uniform.** CV of totals in knockdown cells 0.032 / 0.041 / 0.035 (D/E/F) and 0.014 (Eval A), against 0.47 / 0.46 / 0.46 in the controls and **0.36 in Eval A's real knockdown cells** | **flag:** a frozen-generator property (the "exact pseudobulk totals" path). It narrows high-expression genes' spread (item 5) and makes the Wilcoxon test sharper per call. It has been in every scored cell run, so all our Eval A numbers include it. Changing it would need a generator change and a full cell re-validation |

# 13. 2026-11-02 — Step 21: final proxy round and the cell run of row 4
*Task: Jose, 2026-11-02 (consolidated prompt; run 2026-11-03, 07:15–09:41 UTC).*
- **Not done:** no D/E/F file, no upload, nothing from Part H.
- **Code:** `code/vcc6/step21.py` (chain scripts `step21*.sh` alongside).
- **Outputs:** `results/stage1/2026-11-02_step21/{tables,logs,cells}` (`matrix/` holds C7's small scratch files).
- **Frozen inputs** were imported or read only: the generator, the harness, `run_line.py`, `boot_score.py`, config 1′'s files, v7 and the Step 18–20 outputs.
- **Kept outputs:** nothing from Step 21 had run before this prompt, so nothing was kept or skipped. The prompt was revised mid-run (C3's coupled row, C7d, the C8 additions, C11). Every part was built to the revised version, because none of them had started by then.
- **The MSE simulator task** was not run (cancelled).

**Row 4** = spec (hybrid T 5 calls + config 1′ + no-source fix) + the fixed Part C shift + the F5 term.
- **H1:** F5 term cross-fitted (each half from the model trained on the other half); the shift uses the other half's (b) cut points.
- **Fold pools:** F5 from all of H1 under the fold's basis; shift cut points all-H1 (3.0 / 67.7).
- **HCT116final, screen setup** (`HCT116s` in the tables): 15c's generator input (HCT116's own tables out, HEK293T votes). F5 with sources excluding HCT116's tables, the E4w−Orion basis and β from all of H1 (`tables/f5_HCT116s.npz`); shift cut points all-H1.
  - **Check:** the same code with both Orion lines out rebuilds Step 20's HCT116 fold prediction exactly (0 after casting to float32, as Step 20 stored it).
  - 277 of 278 targets have a source background, and 260 have HEK293T.

**Summary:**
- **Track 1, cells:** row 4 scores **0.2895 [0.2695, 0.3097] on Eval A**. That is **+0.0469 [+0.032, +0.061] over the live entry** and +0.0268 [+0.013, +0.040] over the hybrid T 5 cells (Step 19).
  - **HCT116final 0.3621:** +0.0005 [−0.022, +0.032] against the live entry, and +0.0075 [−0.013, +0.038] against hybrid T 5.
  - **Proxy minus cells** is +0.003 overall on both lines; reach accounts for +0.032 / +0.022 of it.
- **Track 2, proxy:** **no new part met its criterion** (B, C3, C4, C6, C7c both arms, C7d-b, C9c), so **Arm N = none**.
  - C4 came closest: it gained on both H1 halves but cost every guardrail 0.006–0.029.
  - The diagnostics (C2, C5, C7a, C7d-a, C8, C9b, C10, C11) found no usable new signal.
  - Neither C8 nor C9b should replace (b) in the switch.
- **Recommendation:** File 1 = **row 4** (cell evidence above); File 2 = spec + fixed shift (fallback). Both await your go (Part H).

**Deviations, flagged:**
3. **Three of my own bugs stopped a part and were fixed, then resumed:**
   - f5screen's first check compared float64 against Step 20's float32-stored prediction (diff 5.9e-8). The check is now an exact float32 comparison, which passes.
   - C7c's summary failed when a pick was the base row (γ = 0, outside the grid). The rows were kept; picks of the base row are now labelled "base row (row 4), outside the grid".
   - C11 failed on `HipSci_pooled` (pooled donors, not a Cellosaurus line), then on RPE1's ID (`hTERT-RPE1`). I also fixed the TP53 parser: Cellosaurus' text format holds variants in `CC   Sequence variation:` lines.

## 13.1 Track 1: row 4 on cells

**Input checks** (`tables/t1inputs.json`): both pass.

| Check | Eval A | HCT116final (screen) |
|---|---|---|
| calls = hybrid T 5 sets (15c K4 input where b ≥ 5, K0 input elsewhere), bit for bit | yes (176 of 292 knockdowns on K4) | yes (135 of 278) |
| F5 term and shift removed: lfc = 15c K0 input bit for bit; mu = spec mu exactly | yes / yes | yes / yes |
| proxy of that input = spec proxy | 0.268192 = 0.268192 (diff 0) | 0.360018 = 0.360018 (diff 0) |
| proxy of the row-4 input vs row 4 from mu | 0.292798 vs 0.292798 (4.9e-11) | 0.365145 vs 0.365145 (4.1e-12) |
| entries changed / n_sources fixed (gated no-source genes) | 2,951,504 / 28,364 | 2,468,374 / 5,923 |

**Six components** (cells; 95% CIs from 2,000 knockdown resamples; the rebuild reproduces the harness summary to 9e-16):

| Eval A | MSE | PDS | NMAE | fidelity | reach | Jaccard | overall |
|---|---|---|---|---|---|---|---|
| scaled | 0.187 [0.124, 0.240] | 0.766 [0.708, 0.820] | 0.094 [0.072, 0.116] | 0.469 [0.446, 0.492] | 0.079 [0.036, 0.126] | 0.142 [0.119, 0.169] | **0.2895** [0.2695, 0.3097] |
| raw | 0.787 [0.736, 0.848] | 0.854 [0.826, 0.879] | 0.897 [0.884, 0.910] | 0.447 [0.429, 0.465] | 0.191 [0.155, 0.232] | 0.070 [0.060, 0.081] | |
| **HCT116final** scaled | 1.000 [0.952, 1.000] | 0.903 [0.697, 1.117] | 0.250 [0.212, 0.289] | 0.903 [0.889, 0.917] | 0.729 [0.618, 0.838] | −1.612 [−1.622, −1.601] | **0.3621** [0.3174, 0.4024] |
| raw | 0.730 [0.538, 0.892] | 0.669 [0.635, 0.704] | 0.843 [0.818, 0.866] | 0.534 [0.526, 0.541] | 0.661 [0.606, 0.716] | 0.007 [0.005, 0.009] | |

**Paired differences** (same resamples):

| | MSE | PDS | NMAE | fidelity | reach | Jaccard | **overall** |
|---|---|---|---|---|---|---|---|
| Eval A vs **live entry** (15c config 1′ + K0) | +0.061 [+0.009, +0.099] | −0.029 [−0.068, +0.010] | +0.086 [+0.070, +0.103] | +0.097 [+0.076, +0.119] | +0.019 [−0.004, +0.041] | +0.048 [+0.035, +0.061] | **+0.0469 [+0.032, +0.061]**, P(better) 1.00 |
| Eval A vs **hybrid T 5** (Step 19 Part A) | +0.061 [+0.009, +0.099] | −0.029 [−0.068, +0.010] | +0.086 [+0.070, +0.103] | +0.011 [+0.008, +0.014] | +0.024 [+0.002, +0.047] | +0.008 [+0.004, +0.012] | **+0.0268 [+0.013, +0.040]** |
| HCT116final vs live entry | 0.000 [0.000, +0.117] | +0.042 [−0.069, +0.160] | +0.026 [+0.017, +0.037] | −0.015 [−0.024, −0.004] | −0.009 [−0.080, +0.064] | −0.042 [−0.054, −0.031] | **+0.0005 [−0.022, +0.032]**, P 0.61 |
| HCT116final vs hybrid T 5 | 0.000 [0.000, +0.119] | +0.042 [−0.069, +0.160] | +0.026 [+0.017, +0.036] | +0.000 [−0.002, +0.002] | −0.024 [−0.093, +0.043] | +0.001 [0.000, +0.002] | **+0.0075 [−0.013, +0.038]** |

- **MSE, PDS and NMAE are the same against both references.** The K0 and K4 15c runs have the same pseudobulk means (Step 19 Part 0), so the hybrid panel's mean-based components equal K0's.
- **Against hybrid T 5,** the row-4 gain is therefore the F5 term + shift + no-source fix on the means: MSE +0.061 and NMAE +0.086, at a PDS cost of −0.029 (CI includes 0). There is a little more on the DE components.
- **H1 halves on cells:** 0.3105 / 0.2684, each +0.048 / +0.046 over the live entry. Half 1 is the clean held-out half for the F5 configuration (chosen on half 0 in Step 20).

**Proxy minus cells** (row 4; the proxy rows are Part A's):

| | MSE | PDS | NMAE | fidelity | reach | Jaccard | overall |
|---|---|---|---|---|---|---|---|
| Eval A | −0.003 | −0.000 | −0.000 | −0.004 | **+0.032** | −0.005 | **+0.003** |
| HCT116final | 0.000 | −0.001 | −0.001 | −0.002 | **+0.022** | −0.000 | **+0.003** |

The proxy's reach bias (+0.02–0.05, 15c §4) holds here too. Everything else agrees to ±0.005.

**By true-strength tercile, Eval A** (15c §1b terciles; paired CIs within the tercile, full-panel scaling):

| Tercile (n) | overall | MSE | PDS | vs live entry | vs hybrid T 5 |
|---|---|---|---|---|---|
| weak (97) | 0.419 | 0.980 | 0.541 | −0.001 [−0.076, +0.034] | +0.006 [−0.069, +0.040] |
| mid (97) | 0.249 | 0.059 | 0.771 | +0.024 [+0.008, +0.044] | +0.017 [+0.003, +0.034] |
| **strong (98)** | 0.307 | 0.190 | 0.985 | **+0.087 [+0.067, +0.108]** | **+0.034 [+0.017, +0.050]** |

**Requested against realised calls** (`tables/t1_calls_per_kd_*.csv`, `t1_calls_by_decile_*.csv`; realised = harness `de_pred` p_adj < 0.05, own gene out):

| | requested median / mean | realised median / mean | requested-not-realised / realised-not-requested (total) | Spearman with true DE count: realised / requested |
|---|---|---|---|---|
| Eval A | 1,200 / 906 | 1,180 / 913 | 4,843 / 7,070 | **0.62** / 0.56 (live entry's cells: 0.09) |
| HCT116final | 778 / 684 | 738 / 659 | 7,350 / 296 | 0.40 / 0.41 |

The hybrid calls fix the allocation problem the checklist flagged (item 4): realised call counts now follow knockdown strength.

## 13.2 Track 2: the proxy round
*Proxy note, for every result below: the proxy overstates reach by 0.02–0.05 compared with cells (15c §4; +0.032 / +0.022 for row 4 in 13.1). Each pick's reach share of the component gain is given.*

### Part 0 — checks (`tables/part0.json`): all pass
| Check | Got | Step 20 | max \|diff\| |
|---|---|---|---|
| spec, H1 full / half 0 / half 1 | 0.2682 / 0.2892 / 0.2471 | Part G, same | 2.8e-17 |
| spec + F5 term, H1 | 0.2893 / 0.3095 / 0.2690 | same | 5.6e-17 |
| spec and spec + F5 on K562 / HCT116 / HEK293T | 0.1471, 0.1517 / 0.3017, 0.3147 / 0.1530, 0.1651 | same | ≤ 5.6e-17 |
| fixed shift (−0.25, 0, 0.25), fix-only base, all-H1 cuts | 0.25378 (0.2768 / 0.2306) | C2: 0.25378 | 1.2e-12 |

### Part A — the four combinations (`tables/partA_rows.csv`, `partA.json`)
| Line | Row | overall | MSE | **PDS** | NMAE | fidelity | reach | Jaccard | calls / kd median / mean |
|---|---|---|---|---|---|---|---|---|---|
| H1 full (half 0 / half 1) | spec | 0.2682 (0.2892 / 0.2471) | 0.127 | **0.803** | 0.008 | 0.459 | 0.075 | 0.137 | 1,200 / 906 |
| | spec + shift | 0.2748 (0.2975 / 0.2519) | 0.135 | **0.825** | 0.014 | 0.459 | 0.079 | 0.137 | same |
| | spec + F5 | 0.2893 (0.3095 / 0.2690) | 0.196 | **0.742** | 0.092 | 0.465 | 0.105 | 0.137 | same |
| | **row 4** | **0.2928 (0.3145 / 0.2709)** | 0.185 | **0.766** | 0.094 | 0.465 | 0.110 | 0.137 | same |
| K562 | spec / +shift / +F5 / **row 4** | 0.1471 / 0.1477 / 0.1517 / **0.1532** | 0.138 / 0.138 / 0.172 / 0.179 | 0.627 / 0.629 / 0.591 / 0.597 | | | | | 575 / 670 |
| HCT116 (fold) | | 0.3017 / 0.3061 / 0.3147 / **0.3153** | 0.876 / 0.901 / 0.921 / 0.963 | 0.844 / 0.869 / 0.858 / 0.837 | | | | | 1,200 / 853 |
| HEK293T | | 0.1530 / 0.1608 / 0.1651 / **0.1695** | 0.132 / 0.147 / 0.131 / 0.154 | 0.424 / 0.445 / 0.454 / 0.448 | | | | | 1,200 / 916 |
| HCT116final, screen | | 0.3600 / **0.3517** / 0.3795 / **0.3651** | 1.000 throughout | 0.864 / 0.864 / 0.932 / 0.901 | 0.223 / 0.227 / 0.245 / 0.249 | | 0.784 / **0.730** / 0.810 / **0.751** | | 778 / 684 |

- **On H1 the parts add:** the shift +0.0066 and the F5 term +0.0211 over the spec give +0.0246 together. The shift restores 0.024 of the F5 term's PDS cost.
- **On every fold pool, row 4 is the best of the four rows.**
- **HCT116final in the screen setup is the exception: the shift costs 0.008 there** (spec 0.3600 → 0.3517; spec + F5 0.3795 → row 4 0.3651), all of it through reach (−0.054 / −0.059).
  - On cells, row 4 still lands at the live entry's level (13.1), so the loss doesn't show beyond noise there.
  - **Its reach is the component the proxy overstates.**

### Part B — shared-component removal (`tables/partB.json`, `partB_rows.csv`): criterion not met
- **Picks:** η 0 (row 4 itself) in both directions; grid position 1/5, edge. Held-out gain 0.
- **H1 full panel:** η 0 / 0.25 / 0.5 / 0.75 / 1 → 0.2928 / 0.2908 / 0.2862 / 0.2820 / 0.2772.
  - PDS 0.766 → 0.742 and MSE 0.185 → 0.151.
  - **Removing the mean predicted Δ lowers PDS as well as MSE.** The shared component in row 4's predictions is not what holds PDS down.

### Part C1 — C+G's pooled quantities (`tables/partC1.json`, `C1_<context>.npz`, `C1_uncalled_{H1,HCT116s}.npz`; one table load, 07:40–07:59)
- **Saved arrays:** mu, var_fe, τ² and τ²_untruncated were not saved anywhere, so the model was rerun once, as `predict_line` builds C+G, for the 8 contexts. The voting sets are the 5 distinct exclusions (H1; K562; both Orion lines; HCT116 only; none).
- **Checks, every context: shrunk LFC, calls, n_sources and SE_meas reproduce the saved arrays bit for bit** (fold caches, 15c HCT116final input, Step 18 D/E/F inputs).
- **For C3's coupled row,** the Δ pool and agreement z were rebuilt for H1 and HCT116s the way 15c builds them; `config_lfc` with K0's calls reproduces both 15c inputs bit for bit.

**τ²_g** (C+G's gene-specific DL τ², floored at the global τ²):
- **The global τ² is 1e-4 in every voting set, and 77–94% of estimable (≥ 2 sources) gated pairs sit on that floor.** Quantiles of estimable pairs on H1: 50th 1e-4, 90th 0.0022, 99th 0.23 (D/E/F: 90th 0.028, 99th 0.43).
- **Fallback** (< 2 sources → the line's median, which is 1e-4 everywhere): 4.1% of gated, own-out pairs on H1; 4.9% K562; 4.8% HCT116s; 3.8–5.1% D/E/F; **28% / 30% on the HCT116 / HEK293T fold pools** (both Orion lines out).

### Part C2 — uncertainty-aware calls, diagnostic (`tables/partC2.json`)
**The z_u denominator is var_fe.** On H1, its share of SE_meas² + var_fe + τ² is 0.89–0.94 in every control-CPM decile; SE_meas² is 0.02–0.08 and τ² 0.01–0.05. So z_u ≈ |mu| / √var_fe, the pooled fixed-effect z.

| Line | Rule | calls / kd median / mean (10th–90th) | Spearman with true DE count / with (b) | precision / recall / sign accuracy |
|---|---|---|---|---|
| H1 | requested K0 | 487 / 551 (53–1,098) | **0.09** / 0.16 | 0.20 / 0.11 / 0.77 |
| | hybrid T 5 (row 4) | 1,200 / 906 (168–1,200) | **0.56** / 0.70 | 0.21 / 0.18 / 0.76 |
| | z_u BH 0.01 / 0.05 / 0.2 | 0 / 3.4, 6.7, 15.0 (0–10, 16, 31) | 0.51, 0.56, 0.59 / 0.70–0.78 | 0.74–0.64 / 0.002–0.009 / 0.96–0.93 |
| K562 | requested K0 / hybrid | 113 / 254 · 575 / 670 | 0.11 / −0.05 · 0.32 / 0.64 | 0.014 / 0.28 / 0.80 · 0.008 / 0.41 / 0.79 |
| | z_u BH 0.05 | 0 / 1.6 | 0.31 / 0.53 | 0.14 / 0.016 / 0.98 |
| HCT116 (fold) | K0 / hybrid / z_u 0.05 | 461 / 556 · 1,200 / 853 · 0 / 4.6 | 0.11 · 0.27 · 0.42 | precision 0.006 · 0.004 · 0.06 |
| HEK293T | same | 592 / 690 · 1,200 / 916 · 0 / 5.3 | 0.10 · 0.22 · 0.35 | 0.006 · 0.006 · 0.045 |
| HCT116final, screen | same | 133 / 196 · 778 / 684 · 0 / 1.6 | −0.00 · 0.41 · 0.43 | 0.016 · 0.007 · 0.16 |
| D / E / F | requested K0 | 243 / 353 · 279 / 381 · 245 / 356 | — / −0.07 · −0.11 · −0.11 | — |
| | z_u BH 0.05 | 0 / 1.6–1.9 | — / 0.49–0.51 | — |

- **Requested K0's call count is unrelated to knockdown strength** (Spearman 0.09 on H1, −0.00 on HCT116final), as on Eval A's realised cells (0.09). The hybrid rule lifts that to 0.56 / 0.41.
- **BH on z_u calls almost nothing** (median 0 per knockdown everywhere): its few calls are precise and rank knockdowns well, but recall is < 1%.
- By true tercile on H1 (mean calls): K0 524 / 545 / 586; hybrid 639 / 923 / 1,153; z_u 0.2 1.5 / 8.6 / 34.5.

### Part C3 — calls only (`tables/partC3_rows.csv`, `partC3_matched.csv`, `partC3.json`): criterion not met
Magnitudes are fixed at row 4; only the call set changes. The extra coupled row lets magnitudes follow the call rule.

| Rule (H1 full) | overall | fidelity | reach | Jaccard | calls / kd mean | Spearman calls ~ true DE |
|---|---|---|---|---|---|---|
| **R-ref** (hybrid T 5) | **0.2928** | 0.465 | 0.110 | 0.137 | 906 | 0.56 |
| R-u 0.01 / 0.05 / 0.2 | 0.191 / 0.193 / 0.195 | **−0.088 / −0.083 / −0.075** | 0.161 | 0.03–0.04 | 3 / 7 / 15 | 0.50 / 0.56 / 0.59 |
| R-uK4 | 0.2920 | 0.465 | 0.129 | 0.113 | 906 | 0.56 |
| R-uu 0.01 / 0.05 / 0.2 | 0.257 / 0.257 / 0.258 | 0.24 | 0.135 | 0.12 | 723 | 0.65 |
| R-low 50 / 100 / 200 / 300 | 0.2867 / 0.2913 / 0.2938 / **0.2943** | 0.435–0.484 | 0.112–0.116 | 0.125–0.126 | 743–841 | 0.64 |
| R-u 0.05, **coupled** | 0.2067 | −0.083 | 0.191 | 0.032 | 7 | — (MSE 0.242 vs 0.185: dropping calls moves the magnitudes to config 1′'s uncalled values) |

- **Picks:** R-low 300 in both directions (12/12, edge). Held-out +0.0024 / +0.0007 (reach share 0.03 / 0.45).
- **Guardrails:** K562 −0.002, HCT116 +0.002, **HEK293T −0.012**, HCT116final −0.003.
- **Matched mean call volume** (rules above R-ref on H1 full: R-low 200 and 300). R-ref with the b ≥ 5 budget N set to match:
  - H1: 0.2896 (N 1,028) against R-low 200's 0.2938, and 0.2909 (N 1,094) against R-low 300's 0.2943. **R-low keeps +0.003–0.004 on H1 at matched volume.**
  - Guardrails at matched volume: R-ref is ahead on K562 (+0.003) and HEK293T (+0.011), and about level on HCT116 and HCT116final.
- **On the guardrails, R-u's tiny call sets score far higher on the proxy** (HEK293T 0.29–0.31, HCT116final 0.53–0.55, against 0.17 / 0.37). Their sparse truth makes "call almost nothing" look good. I read this as a volume artefact of the proxy, not a signal: on H1 the same rules lose 0.10.

### Part C4 — uncertainty-aware magnitudes (`tables/partC4_rows.csv`, `partC4.json`): criterion not met (guardrails)
"Called" here means config 1′'s magnitude split (K0's calls). Those genes get s × C+G LFC × z_u² / (z_u² + κ); everything else is as row 4.

| H1 full | MSE | PDS | NMAE | overall | weak / mid / strong |
|---|---|---|---|---|---|
| row 4 | 0.185 | 0.766 | 0.094 | 0.2928 | 0.418 / 0.254 / 0.314 |
| κ 1: s 0.5 / 0.75 / 1 | 0.221 / 0.162 / 0.070 | 0.762 / 0.775 / 0.780 | 0.092 / 0.092 / 0.082 | 0.2997 / 0.2904 / 0.2738 | |
| κ 4: s 0.5 / 0.75 / 1 | **0.245** / 0.227 / 0.194 | 0.749 / 0.766 / 0.776 | 0.085 / 0.091 / 0.091 | **0.3032 / 0.3027** / 0.2981 | κ 4 s 0.5: 0.422 / **0.301** / 0.313 |

- **Picks:** (κ 4, s 0.75) from half 0 (κ 2/2 edge, s interior) and (κ 4, s 0.5) from half 1 (both edges).
- **Held out:** +0.0069 / +0.0125, almost all MSE (+0.046 / +0.054); reach share 0.01 / 0.36.
- **Guardrails: −0.006 to −0.008 at the first pick, and −0.016 to −0.029 at the second.** The criterion fails on the guardrails.
- **Where the gain sits:** shrinking called genes with low z_u helps H1's mid tercile (MSE 0.053 → 0.333) and does nothing for the strong tercile.
- **For scale:** Step 16's arms j / kj (0.229 / 0.237 against 0.244) also lost on H1. C4 wins on H1 but loses every other line.

### Part C5 — up against down, and sign calibration (diagnostic, H1, row 4; `tables/partC5.json`)
1. **Plug-in error split.** Strong tercile: error and real-change energy are spread across the CPM bands in proportion (20–100 CPM holds 47% of each; > 100 CPM 26% / 28%; ungated 9% / 7%).
   - The excess is in **wrong-sign cells:** true down / predicted up (20–100 CPM) holds 10.7% of the error but 4.9% of the energy; true up / predicted down 7.2% vs 3.5%.
   - **Mid tercile:** wrong-sign cells together hold about 37% of the error against about 16% of the energy.
2. **Sign accuracy on true DE genes:**
   - called predicted-up **0.72**, called predicted-down **0.81**; uncalled 0.63 / 0.73.
   - By z_u decile it rises steadily from 0.63 to 0.83 (top decile, z_u > 1.55). By (b) tercile 0.53 / 0.67 / 0.71.
3. **Calibration slope** (Δ space, through the origin, called genes), against the 2p − 1 implied by sign accuracy:

| | 5–20 CPM | 20–100 CPM | > 100 CPM |
|---|---|---|---|
| predicted up | 0.69 (implied 0.54) | 0.48 (0.50) | **0.40** (0.40) |
| predicted down | 0.38 (0.61) | 0.39 (0.57) | 0.55 (0.63) |

- The slopes are about what sign accuracy alone allows (a shrink of 0.4–0.7).
- **Predicted-up genes at > 100 CPM are both the least sign-reliable and the most over-predicted.** This motivated C6's s_up.

### Part C6 — asymmetric called-gene scale (`tables/partC6_rows.csv`, `partC6.json`): criterion not met
- **Picks:** (s_up 0.5, s_down 0.35, c_floor 20) in both directions, a corner of the grid (all three edges). Held-out +0.0002 / +0.0007.
- **Guardrails:** K562 −0.011, HCT116 −0.015, HEK293T −0.003, HCT116final −0.012.
- **Marginal means, H1 full** (overall / MSE / PDS):
  - s_up 0.5 → 1.0: 0.2915 → 0.2703; MSE 0.181 → 0.082; PDS flat.
  - s_down 0.35 / 0.5 / 0.65: 0.2843 / 0.2838 / 0.2796; MSE 0.163 → 0.111; PDS 0.751 → 0.773.
  - **Larger scales cost MSE, smaller ones cost PDS.** No asymmetry beats config 1′'s 0.5 on both.
- **Compared with C4:** C4's held-out H1 gain was 10× larger, but both lose the guardrails.

### Part C7 — co-response modules (`tables/partC7.json`, `partC7a_per_kd.csv`, `partC7b_rows.csv`, `partC7c.json`, `c7build.json`)
**Matrix:** E4w exactly as Step 20 builds it: 65,723 rows; the rows match Step 20's and project onto its basis with the saved scores (relative diff 2e-6).
- **c_X is computed per target** with the scored line's tables dropped (full / −K562 / −Orion / −HCT116), rows targeting X dropped, only tables that gate X, row weight w².
- **Rows kept per target:** median 65,715 (full); 56,153 (−K562); 33,290 (−Orion); 49,227 (−HCT116). 10th percentile 18,000–51,000.
- **Low rank against full** (median Pearson over the line's gate): K 20 0.85–0.87; K 100 0.92–0.93.

**C7a** (no fitting; per-knockdown Pearson of true Δ with −c_X):

| Line | median r (null; percentile) | strong tercile | true DE in top 100 co-movers (null) |
|---|---|---|---|
| H1 | −0.004 (−0.000; 12th) | −0.003 (38th) | 0.034 (0.013) |
| K562 | −0.005 (0.000; 8th) | −0.009 (2nd) | 0.026 (0.018) |
| HCT116 | −0.008 (0.001; 5th) | −0.013 (6th) | 0.019 (0.014) |
| HEK293T | +0.007 (0.002; 84th) | +0.010 (79th) | 0.039 (0.015) |
| HCT116final, screen | −0.010 (−0.000; 3rd) | −0.016 (4th) | 0.019 (0.014) |

- **For scale:** F6 (own-control co-expression) H1 median r 0.003. F3d's source-score baseline on H1 is a pair Spearman of 0.042 (mid + strong 0.083).
- **Response co-variation with X does not predict X's knockdown.** The correlation is at or below the null on 4 of 5 lines. True DE genes are enriched 1.4–2.6× among X's co-movers, but without the predicted sign.

**C7b** (Step 20 F4's setup; held-out median R² · strong-tercile captured energy). M2, M3 DepMap and M4 DepMap reproduce Step 20's F4 exactly.

| K | M2 | M3 (DepMap) | **M3e** | M4 (DepMap) | **M4e** | M4e null (20 shuffles) |
|---|---|---|---|---|---|---|
| 20 | 0.174 · 0.070 | 0.205 · 0.072 | 0.178 · 0.066 | 0.119 · 0.035 | 0.047 · 0.022 | 0.046 · 0.022 (95th 0.051 · 0.023) |
| 100 | 0.124 · 0.063 | 0.140 · 0.080 | 0.130 · 0.030 | 0.063 · 0.051 | 0.046 · 0.037 | 0.014 · 0.030 (95th 0.027 · 0.033) |

The target's own loadings carry little response information. At K 20, M4e equals its null. At K 100 it beats the null (R² 0.046 against a 95th percentile of 0.027), but stays below DepMap's M4 (0.063 · 0.051). Added to the source-line ridge, M3e is below M3 at both K (strong-tercile capture at K 100: 0.030 against 0.080).

**C7c** (proxy arms on row 4): criterion not met for either.
- **Co-movement arm:** picks γ 0 (row 4) and γ 0.5 (1/3, edge); held-out ±0.0000. λ_t per (b) tercile in `partC7c_lambda.json`.
- **F5 with M3e in place of M3:** held-out −0.026 / −0.001. Guardrails −0.004 to −0.010. H1 full 0.2789 (MSE 0.145, PDS 0.754, NMAE 0.071; row 4: 0.185 / 0.766 / 0.094). By tercile 0.419 / 0.251 / 0.293 vs row 4's 0.418 / 0.254 / 0.314.

### Part C7d — X's column as a dose-response regression (`tables/partC7d.json`, `partC7d_info_<line>.csv`, `partC7db_rows.csv`)
- **Rows:** C7's exclusions plus |z_X| ≥ 3 in the row's own table. Median kept per target: 1,021 (H1), 688 (K562), 541–565 (Orion lines). Slopes are per source line (≥ 3 rows), pooled one vote per line.
- **No prediction (< 20 rows):** 0% of H1's targets, 0.4% K562, 0.7% HCT116 / HEK293T. For b < 5: 0–1.4%. D/E/F 0.3% (0.6% for b < 5).
- **Reliability ratio:** median 0.92; none clipped. SE of Δ_X ≈ |Δ_X| / |z_X|, since the LFC → Δ factor cancels.
- **Premise:** true own-gene change vs d_X (80% knockdown):
  - H1 median −1.51 (IQR −1.87 to −1.14) vs −1.13: ratio 1.37, so H1's knockdowns are stronger than 80%;
  - HCT116 −0.56 vs −0.80: ratio 0.83.
- **Extrapolation:** median |d_X| ÷ median |Δ_X| in the kept rows = 3.9 (H1), 2.3–2.6 (guardrails).
- **C7d-a** (median per-knockdown Pearson; capped = uncapped, because capping rescales each knockdown's prediction by one positive factor and leaves r unchanged):

| Line | all | weak / mid / strong | b < 5 | C7a on the same knockdowns (null) |
|---|---|---|---|---|
| H1 | −0.005 | −0.003 / −0.006 / −0.010 | −0.003 | −0.004 (−0.000) |
| K562 | −0.004 | −0.005 / −0.003 / −0.002 | −0.006 | −0.004 (0.001) |
| HCT116 | +0.001 | −0.002 / +0.006 / −0.004 | +0.003 | −0.007 (−0.000) |
| HEK293T | −0.000 | +0.003 / +0.003 / −0.003 | +0.003 | +0.008 (0.002) |

- **C7d-b** (row 4 + γ Δ̂ on uncalled gated genes; every knockdown with a prediction): criterion not met.
  - **Picks: γ 0 (row 4) in both directions.**
  - Even γ 0.25 drives MSE to 0 on every line (H1 0.2928 → 0.1949; b < 5 0.361 → 0.179; b ≥ 5 0.292 → 0.208), and NMAE turns negative. The term's magnitude is far too large, which is the extrapolation above.
- **Caveat:** these rows show genes that move *with* X, not *because of* X. Rows where X fell because its upstream regulator was knocked down also move that regulator's other targets.

### Part C8 — strength from the target line's expressed downstream genes, plus b_own and DepMap (measure only; `tables/partC8.json`, `C8_counts.csv`)
Spearman with the true DE count (partial given (b); partial given b_expr5 for b_own and DepMap):

| Predictor | H1 | K562 | HCT116 (fold) | HEK293T | HCT116final screen |
|---|---|---|---|---|---|
| **(b)** | **0.642** | **0.408** | **0.497** | **0.409** | **0.497** |
| (b) recomputed, z > 3 (b_z3) | 0.621 (+0.22) | 0.324 (+0.00) | 0.465 (+0.05) | 0.318 (−0.11) | 0.432 (+0.06) |
| b_expr5 | 0.619 (+0.23) | 0.315 (+0.01) | 0.461 (+0.05) | 0.316 (−0.11) | 0.428 (+0.07) |
| b_expr100 | 0.613 (+0.22) | 0.330 (+0.06) | 0.453 (+0.07) | 0.309 (−0.07) | 0.428 (+0.08) |
| b_own (target's control CPM) | 0.075 (+0.02; +0.03) | 0.103 (+0.06; +0.05) | 0.128 (+0.04; +0.05) | 0.064 (+0.02; −0.02) | 0.128 (+0.04; +0.05) |
| DepMap mean gene effect | −0.470 (−0.25; −0.22) | −0.181 (−0.16; −0.17) | −0.195 (−0.14; −0.11) | −0.115 (−0.04; −0.04) | −0.195 (−0.14; −0.14) |

- **The expression filter adds nothing.** b_expr5 ≈ b_z3 everywhere, and both are below (b) itself.
- On H1, a z > 3 count adds +0.22 beyond (b); on the guardrails it adds −0.11 to +0.07.
- **b_own** is nearly uninformative. **DepMap** (essential genes knock down harder) adds −0.14 to −0.25 beyond (b) on three of four lines.
- **Medians, H1 vs D/E/F:** (b) 21.9 vs 4.0; b_expr5 62.6 vs 45–46; b_z3 73.6 vs 55; b_own 109 vs 56–64 CPM. Within-(b)-tercile Spearmans are in the json.
- **Recommendation: none of these should replace (b) in the hybrid switch.** (b) has the highest Spearman on every line. DepMap's partial is the only consistent extra signal, but (b) ≥ 5 is a switch, not a regression; adding DepMap to it would need its own cell test.

### Part C9 — knockdown clusters (`tables/partC9.json`, `partC9_clusters.pkl`, `partC9c_rows.csv`, `partC9c.json`)
- **Embeddings:** E24w (strong rows, n_real ≥ 50, K 20) and E4w (all rows, K 100), from Step 20's saved row scores, per target averaged per background then over backgrounds.
- **Bases:** fold bases for the fold pools; the full basis with HCT116's rows dropped for HCT116final screen.
- **Clusters:** cosine k-means (unit vectors, 10 starts, seed 2026).

**C9a, sanity** (share of member pairs in one cluster, full variant; random pairs 0.037 / 0.014 / 0.030 / 0.010):

| Complex (members, E24w / E4w) | E24w k 50 / 150 | E4w k 50 / 150 |
|---|---|---|
| ribosome L (46 / 54) | 0.62 / 0.37 | 0.63 / 0.40 |
| ribosome S (31 / 35) | 0.82 / 0.51 | 0.78 / 0.78 |
| proteasome (29 / 45) | 0.52 / 0.37 | 0.47 / 0.44 |
| Mediator (28 / 33) | 0.28 / 0.22 | 0.47 / 0.36 |
| mitochondrial ribosome (76 / 82) | 0.14 / 0.10 | 0.30 / 0.28 |
| CCT (7 / 9) | 1.00 / 1.00 | 0.58 / 0.58 |
| exosome (11 / 12) | 0.51 / 0.65 | 0.24 / 0.55 |
| Integrator (16 / 16) | 0.12 / 0.05 | 0.09 / 0.08 |
| BAF (9 / 25) | 0.06 / 0.06 | 0.02 / 0.01 |
| TFIID (13 / 14) | 0.18 / 0.32 | 0.41 / 0.73 |
| spliceosome (77 / 97; HGNC Spliceosomal A–P + Sm) | 0.17 / 0.08 | 0.37 / 0.15 |
| NADH dehydrogenase (19 / 37; core + supernumerary) | 0.15 / 0.10 | 0.06 / 0.03 |

- Clusters are 5–100× random for the essential machines (ribosome, CCT, proteasome, exosome, TFIID).
- They are near random for BAF, Integrator and complex I.

**C9b, strength prior:**
- **Routes:**
  - E24w: own embedding for 135 of 292 H1 targets, STRING for 125, none for 32; on the guardrails only 51–65 own, 142–153 STRING, 62–74 none.
  - E4w: 250–290 own, the rest STRING.
- **b_clust's partial Spearman given (b) is −0.15 to +0.08** on every line, embedding and k.
- Raw Spearmans: E4w k 50 0.46 (H1) and 0.13–0.23 (guardrails), all explained by (b). E24w −0.23 to +0.04.
- **For b < 5:** −0.20 to +0.15.
- **D/E/F:** b_clust median 1 (E4w) / 90–98 (E24w).
- **No knockdown-cluster strength prior beyond (b).**

**C9c, cluster-mean response** (E24w, k 150, declared before the run; the mates' mean pooled raw source Δ, the scored line's tables and the target out; uncalled gated genes; b ≥ 5): criterion not met.
- **Picks: γ 0 (row 4) in both directions.**
- **Reference:** Step 11's partner arm (h) lowered H1 (0.201–0.208 against 0.231). S3: 0 of 300 D/E/F targets had a K562 partner above cosine 0.5.

### Part C10 — cascades (measure only; `tables/partC10.json`, `partC10_H1.csv`)
- **X candidates:** genes with an own source knockdown at (b) ≥ 50: 2,361, of which 2,254 have pooled data. Pooled z from C1 (mu / √var_fe); pooled Δ = the per-background mean of raw source Δ, averaged over backgrounds.
- **H1:** 158 of 176 knockdowns with b ≥ 5 have at least one X (median 8, 90th percentile 64, max 419).
  - **Pearson of T_A with the residual** (truth − row 4): median −0.033 against a random-X null of −0.038 (81st percentile).
  - **With the truth:** 0.046 against 0.037 (99th percentile).
- **D/E/F:** 271 of 300 targets have at least one X (134 with b ≥ 5; context D's pooled z).
- **Reading:** a cascade term from downstream knockdowns tracks the truth a little better than random X sets do, but carries nothing row 4 is missing.

### Part C11 — the p53 program (measure only; `tables/partC11.json`)
- **Gene set:** Fischer 2017 (Oncogene 36:3943, doi:10.1038/onc.2016.502), Table 1, "Top 116 genes identified as activated p53 targets in at least 6 out of 16 genome-wide data sets". Parsed from Europe PMC's full text (PMC5511239) into `lib/downloads/step21_2026-11-03/fischer2017_core116.json` (116 genes, all on the panel).
- **TP53 status** (Cellosaurus API; raw entries cached in `lib/downloads/step21_2026-11-03/cellosaurus/`):

| Source line | Cellosaurus | TP53 | Mean Δ of the set: RPL/RPS knockdowns (n) / strong n_real tercile |
|---|---|---|---|
| RPE1 | hTERT-RPE1, CVCL_4388 | no variant listed (hTERT-immortalised) | **0.452** (63) / 0.286 |
| HepG2 | Hep-G2, CVCL_0027 | no variant listed | **0.374** (71) / 0.220 |
| HCT116 | HCT 116, CVCL_0291 | "None_reported" (wild type) | 0.133 (62) / 0.003 |
| HEK293T | HEK293T, CVCL_0063 | no variant; **transformants adenovirus 5 (E1B) and SV40 large T** | 0.082 (71) / −0.000 |
| K562 | K-562, CVCL_0004 | **p.Gln136fs*13** (null) | 0.051 (75) / 0.016 |
| Jurkat | Jurkat, CVCL_0065 | **p.Arg196Ter** | 0.050 (74) / 0.034 |
| KOLF2.1J | KOLF2.1J, CVCL_B5P3 | no variant listed | 0.098 (4) / 0.038 |
| Feng HipSci clones (19) | HPSI…-<clone> entries | no variant listed | 0.04–0.13 (4–8 each) / −0.02–0.03 |
| HipSci pooled | (not a line) | not looked up | 0.030 (5) / −0.011 |
| A549 / MCF7 | CVCL_0023 / CVCL_0031 | "None_reported" | no RPL/RPS knockdowns / 0.006 / −0.005 |
| BXPC3 / HT29 | CVCL_0186 / CVCL_0320 | **p.Tyr220Cys / p.Arg273His** | none / −0.005 / −0.002 |
| HAP1 | HAP-1, CVCL_5G07 | no TP53 variant listed in the entry | none / −0.006 |

- **The set rises 0.37–0.45 on ribosome knockdowns in RPE1 and HepG2** (p53 intact, no viral oncoproteins).
- It rises 0.05 in K562 and Jurkat (TP53-mutant). HCT116 (wild type) is in between at 0.13, and HEK293T (p53 bound by viral proteins) at 0.08.
- **Arc H1 and HCT116final** (row 4):
  - **neither panel holds an RPL/RPS knockdown** (0 of 292 / 278), so the RPL/RPS comparison can't be made;
  - the set holds 1.5% of H1's strong-tercile error and 1.7% of its real-change energy (0.6% of genes). Strong-tercile mean Δ of the set: predicted 0.036 vs true 0.053 on H1; 0.017 vs 0.002 on HCT116final.
  - **The p53 program is a small part of the error on the scored lines.**
- **Basal CPM** (MDM2 / CDKN1A / TP53 / CDKN2A):
  - D 253 / 259 / 372 / 0; E 180 / 122 / 85 / 0.004; F 138 / 311 / 519 / 941;
  - H1 184 / 6 / 461 / 0.4; HCT116final 127 / 205 / 47 / 207; K562 57 / 0 / 0 / 0; HEK293T 129 / 29 / 62 / 540;
  - per source line in the json.
  - **This is weak evidence of p53 status.** D/E/F all express TP53 mRNA and CDKN1A well above H1's level. E's TP53 is lowest; F looks fibroblast-like (CDKN2A high).

### Part D — recommendation (`tables/partD.json`)
- **Parts that met their criterion: none** (B, C3, C4, C6, C7c co-movement, C7c F5-M3e, C7d-b, C9c). No combination row to score; **Arm N = none.**
- **C8 or C9b in place of (b):** no (C8, C9b above).

**Candidates against the live entry** (live Eval A 0.2426, HCT116final 0.3617 on cells):

| Candidate | Eval A | HCT116final | Evidence |
|---|---|---|---|
| **Row 4** (spec + fixed shift + F5 term + no-source fix) | **+0.0469 [+0.032, +0.061]** (→ 0.2895) | **+0.0005 [−0.022, +0.032]** (→ 0.3621) | **cells (Track 1)**; proxy agrees to +0.003 |
| Spec + fixed shift (File 2, fallback) | ≈ +0.028: hybrid T 5 +0.0201 (cells) + fix +0.0013 + shift +0.0066 (proxy) | ≈ −0.015: hybrid T 5 −0.0070 (cells) + fix +0.0005 + shift −0.0083 (proxy, all reach) | cells + proxy |
| Spec alone (Step 19) | ≈ +0.021 | ≈ −0.007 | cells + proxy |

**Recommendation:**
- **File 1 = row 4.** It is the only candidate with a held-out cell gain, and its HCT116final cost against the live entry is nil.
- File 2 as the prompt specifies. On the proxy its shift hurts HCT116final (through reach, the proxy's weakest component).

## 13.3 What I would have done differently
2. **Built C3's guardrail reading in from the start.** On the guardrails' sparse truth, rules that call almost nothing score far higher on the proxy (R-u: HCT116final 0.55). Any call rule with a very different volume needs the matched-volume comparison on every line, not only when it beats R-ref on H1.
3. **Sized C7d's dose before running the arm.** The 3.9× extrapolation and H1's own-change ratio (1.37) were both visible before C7d-b's proxy rows. A dose fitted to the source rows' own range (the capped version) was the sensible arm to score; I scored the uncapped one as specified.
4. **Written C11's Cellosaurus parser against a cached entry first.** Two restarts (pooled donors, RPE1's ID, the `CC Sequence variation` format) cost about 4 min.
6. **Labelled base-row picks in `fit_report` from the start.** Half of the arms picked γ 0, and C7c's summary crashed on it.

# 14. 2026-11-03 — Step 22: the line's own programs, line-vote sign agreement, two quick checks, and reach by ordering
*Task: Jose, 2026-11-03, with three follow-ups during the run: the gene-set matching report, P1 item 1 plus P2b, and the pre-approved E1. Run 2026-11-03, 13:01–15:34 UTC.*
- **Not done:** no D/E/F file and no upload. E1 did not run (R3 failed check 3), so the only cells generated were R3's 20-knockdown slice.
- **Code:**
  - `code/vcc6/step22.py`;
  - `step22.sh` (the chain), `step22_e1.sh` (the E1 / P2b queue), `step22_stage.sh` / `step22_sync.sh` (copies);
  - `gen_v8/cell_generator_v8.py` and `gen_v8/gen_cells_v8.py`, which are copies; the frozen files were not edited. `code/cellgen/cell_generator.py` sha256 is still 422f8a3d…, read-only.
- **Outputs:** `results/stage1/2026-11-03_step22/{tables,logs}`.
- **Frozen inputs** were imported or read only: the generator, the harness, `run_line.py`, `boot_score.py`, config 1′'s files, v7 and the Step 18–21 outputs.

**Summary:**
- **No part met its criterion.** That covers P2, P2b, S3, Q2 (including M1b), and the report-only criteria of P1 and Part R. **Arm N2 = none; row 4 stays File 1.**
  - **Expected Eval A:** unchanged, row 4's cells at 0.2895.
- **P1 (control-derived directions):**
  - H1's own control set adds +4.4 points of mid + strong residual energy beyond E4w K 100, 20× its random match. HEK293T adds +4.2.
  - K562 (+2.9) and HCT116 (+2.4) fall just short of 3 points, so the report-only criterion is not met.
  - The directions are line-specific: other lines' controls capture 5–7% of H1's residual, against 14% for H1's own.
- **P2 / P2b (activating the pathway programs from the sources):**
  - Every P2 row lost on every line, and the base was picked.
  - P2b, with one pathway (cholesterol) and joint ridge slopes, moved nothing: held-out −0.0002 / 0.0000.
  - The sources do not predict the pathway coefficients: P1's Pearson is |r| < 0.2 except cholesterol on H1, which is 0.39.
- **S (line votes):**
  - The lines agree *more* closely than their SEs imply (median Cochran's Q/df 0.16–0.48).
  - p̃ predicts sign about as well as z_u (AUC 0.56–0.65), with little added information.
  - S3's calibrated-weight arm gained only +0.001 on both held-out halves.
- **R (reach by ordering):**
  - **The ceiling is large:** truth-sign-first ordering would add about +0.08 scaled reach per H1 half on the proxy.
  - **No key gets near it:** p_cal captures +1% / −22% of that on H1.
  - **On cells, the order is the whole gap:** re-ordering the proxy by the cells' realised p-values reproduces Track 1's cell reach (0.078 against 0.079). So the +0.032 proxy-minus-cells reach gap is all ordering, and none of it sign.
  - **Generator v8 (rank-graded spread):** means bit for bit, MSE/PDS/NMAE unchanged. But it moved 1.9% of the called set (Jaccard 0.981 < 0.99), so **R3 failed check 3. Part R stopped, and E1 was not run.**

**Deviations, flagged:**
2. **The Tirosh S/G2M lists:** my first URL (a scanpy_usage mirror) returned 404. I used Seurat's own `data/cc.genes.rda` (md5 e84661f9…) read with Rscript: 43 S + 54 G2M genes.
4. **Gene-set matching (Jose's follow-up):** each set uses direct panel matches, or HGNC previous/alias mapping only where that matches more genes.
   - The rule was in place before P0 used the sets. Table below.
5. **HCT116final's E4w basis:** no E4w basis with only HCT116's rows dropped exists. P1 used E4w−Orion (HCT116 and HEK293T both out), the basis row 4's F5 term uses for HCT116final.
6. **S3's normalisation:**
   - The applied weight m^ρ, not m itself, is normalised to mean 1 within each class. That keeps the "class average unchanged" intent at ρ 0.5.
   - It is normalised separately in each H1 half (each scored panel). Pairs with no vote keep weight 1.
8. **Q1's per-tercile changes** on the fold pools needed the two C4 picks rescored with tercile parts. These were the same rows, not new arms.
9. **P2b choices (mine, logged in `partP2b.json`):**
   - P2b's pathway selection used P1's coefficient on v (basis H, as P1 specifies). That was fixed before reading the numbers; v⊥ would select the same pathway.
   - The ridge is fitted on all 10 retained pathway scores together; the term uses the selected pathway's joint slope only.
   - Metacell weights are scaled to mean 1, so λ = 1 sits on the usual 200-row scale.

## 14.1 Part 0 — checks (`tables/part0.json`, `row4_rows.csv`): pass

**Row 4 reproduces Step 21 Part A exactly** (max |diff| ≤ 8e-17):

| line | H1 (halves) | K562 | HCT116 | HEK293T | HCT116final |
|---|---|---|---|---|---|
| overall | 0.2928 (0.3145 / 0.2709) | 0.1532 | 0.3153 | 0.1695 | 0.3651 |

**Single-cell controls:** every line has raw integer counts on disk. The gate recomputed from them (> 5 CPM) equals the line's gate on every line, so no line is left out of P1 or P2. HCT116final uses HCT116's controls, the same file (identical gate).

| line | file | control cells | median UMI per cell | gated genes |
|---|---|---|---|---|
| H1 | `lib/harness/evalA/context.h5ad` | 18,400 | 20,001 | 10,778 |
| K562 (Replogle NT) | `final_evals/harness/lib/evalK562final/context.h5ad` | 18,400 | 7,915 | 7,681 |
| HCT116 (Orion NT) | `…/evalHCT116final/context.h5ad` | 9,940 | 12,133 | 10,326 |
| HEK293T (Orion NT) | `…/evalHEK293Tfinal/context.h5ad` | 9,942 | 15,171 | 11,066 |
| D / E / F | `data_final/extracted/context_{D,E,F}.h5ad` | 18,400 each | 20,236 / 19,681 / 19,845 | 9,633 / 9,743 / 9,537 |

## 14.2 Part P — the line's own programs

### P0 — control bases (`tables/P0_<line>.{npz,json}`, `partP0.json`, `gene_sets.json`)
**Gene-set matching:** "direct" = on the panel as is; "mapped" adds HGNC previous/alias symbols with a unique target.

| set | size | direct | mapped | used |
|---|---|---|---|---|
| S (Tirosh 2016) | 43 | 42 | 43 | mapped (MLF1IP → CENPU) |
| G2M (Tirosh 2016) | 54 | 52 | 54 | mapped (FAM64A → PIMREG, HN1 → JPT1) |
| E2F_TARGETS / G2M_CHECKPOINT / P53_PATHWAY | 200 / 200 / 200 | 193 / 192 / 191 | same | direct |
| UNFOLDED_PROTEIN_RESPONSE / CHOLESTEROL_HOMEOSTASIS | 113 / 74 | 107 / 74 | same | direct |
| MTORC1_SIGNALING / MYC_TARGETS_V1 / INTERFERON_ALPHA_RESPONSE | 200 / 200 / 97 | 189 / 169 / 95 | same | direct |
| TNFA_SIGNALING_VIA_NFKB / HYPOXIA | 200 / 200 | 198 / 195 | same | direct |

The Hallmark genes that don't match are simply not among the 18,533 panel genes (for example GAPDH, LDHA and 31 RPL/RPS genes in MYC_TARGETS_V1). Mapping adds none.

**Metacells, stability and phases.** Every line has M = 200 metacells. The pathway programs v_{c,P} are very stable across random halves of the cells: min r_P is 0.90 on H1 and 0.62–0.97 on the other lines. **All 10 pathways are kept on every line** (rule r_P ≥ 0.5).

| line | cells / metacell (median) | B overlap K 5 / K 20 | phases G1 / S / G2M | r_A (S−G1, G2M−G1) |
|---|---|---|---|---|
| H1 | 92 | 0.92 / 0.52 | **0.03%** / 48% / 52% | — (A empty) |
| K562 | 94 | 0.90 / 0.43 | 10.8% / 48% / 41% | 0.96 / 0.97 |
| HCT116 | 45 | 0.86 / 0.39 | 10.3% / 32% / 58% | 0.95 / 0.96 |
| HEK293T | 52 | 0.84 / 0.33 | **0.02%** / 47% / 53% | — (A empty) |
| D | 92 | 0.97 / 0.54 | 30% / 48% / 22% | 1.00 / 1.00 |
| E | 93 | 0.95 / 0.62 | 35% / 34% / 30% | 0.99 / 0.99 |
| F | 92 | 0.96 / 0.65 | 1.0% / 51% / 48% | 0.93 / 0.94 |

- **Basis A is empty on H1 and HEK293T.** scanpy's `score_genes_cell_cycle` calls a cell G1 only when both its S and G2M scores are negative, and nearly every control cell in those lines scores above zero on one of them. A G1 pseudobulk needs at least 10 cells. This is what the specified method gives; I did not change it.
- **The pathway scores vary much more in D and E's controls than in H1's.** var(s_P) for E2F is 0.117 / 0.103 against H1's 0.001, and G2M is 0.074 / 0.080 against 0.003. In F the pathways are flatter (0.014–0.026).
- **The programs differ by line.** Cosines of v_{c,P} over shared gated genes:
  - D/E/F against H1: −0.19 to 0.59, median about 0.3.
  - D/E/F against K562: higher, 0.35–0.84 for most pathways (a few negative: cholesterol −0.23 on D; UPR −0.13 and mTORC1 −0.12 on F).
  - H1 against K562: 0.27 median.

### P1 — do control-derived directions hold the missing response? (oracle ceiling; `tables/partP1.json`)
Captured energy (least squares, own gene out), mid + strong knockdowns. "C" = A ∪ the 10 pathway programs ∪ B K 20 (30–32 dimensions).

| line | target | B 5 / 20 / 50 | A | H (10) | C | E4w K 20 / 100 | E4w 100 + C | C's increment (random match) | CV own PCs K 5 / 20 |
|---|---|---|---|---|---|---|---|---|---|
| H1 | truth | 0.048 / 0.168 / 0.202 | 0 | 0.115 | 0.189 | 0.213 / 0.345 | 0.390 | +0.045 (+0.002) | 0.272 / 0.407 |
| H1 | residual | 0.038 / 0.125 / 0.150 | 0 | 0.085 | 0.140 | 0.159 / 0.278 | 0.322 | **+0.044** (+0.002) | 0.161 / 0.270 |
| K562 | residual | 0.036 / 0.060 / 0.069 | 0.012 | 0.048 | 0.065 | 0.050 / 0.077 | 0.106 | +0.029 (+0.004) | 0.013 / 0.034 |
| HCT116 | residual | 0.073 / 0.091 / 0.098 | 0.023 | 0.080 | 0.095 | 0.081 / 0.111 | 0.135 | +0.024 (+0.003) | 0.057 / 0.067 |
| HEK293T | residual | 0.148 / 0.159 / 0.165 | 0 | 0.151 | 0.162 | 0.130 / 0.165 | 0.207 | **+0.042** (+0.002) | 0.133 / 0.145 |
| HCT116final | residual | 0.074 / 0.092 / 0.099 | 0.023 | 0.080 | 0.096 | 0.077 / 0.106 | 0.132 | +0.026 (+0.003) | 0.057 / 0.067 |

- Random orthonormal bases of the same dimensions capture 0.000–0.006. Every control basis sits 20–50× above random.
- **Whose controls (H1's residual, mid + strong):** H1's own C 0.140 against D 0.059, E 0.061, F 0.065, K562 0.052, HCT116 0.065, HEK293T 0.052. The directions are line-specific.
- **Reference (section 11):** E4 (Δ) 0.177 / 0.453 at K 20 / 96 on (b); H1's own top PCs in sample 0.30 / 0.53 at K 5 / 20 (strong 0.38 / 0.68). The cross-validated own PCs here (0.16 / 0.27 on the residual) are the out-of-sample version.
- **Criterion (report only), not met.** C adds ≥ 3 points and ≥ 5× random on H1 (+4.4 points, 22×) and on HEK293T (+4.2). K562 (+2.9) and HCT116 (+2.4) fall short, so only one fold pool passes, against the two required.

**Per pathway (item 1, Jose's follow-up):** for mid + strong knockdowns, each knockdown's joint least-squares coefficient of its residual on the 10 v's.
- **r(a):** Pearson of the coefficient with the source activity a_{k,P}.
- **sign:** sign consistency, the share of knockdowns with the majority sign.
- **ρ(b):** Spearman of the coefficient with (b).

| pathway | H1 r(a) / sign / ρ(b) | K562 | HCT116 | HEK293T |
|---|---|---|---|---|
| E2F_TARGETS | +0.02 / 0.67 / −0.06 | +0.10 / 0.55 / +0.06 | −0.06 / 0.51 / −0.24 | −0.12 / 0.57 / +0.32 |
| G2M_CHECKPOINT | −0.06 / 0.66 / +0.07 | −0.09 / 0.54 / −0.08 | +0.01 / 0.60 / +0.16 | −0.08 / 0.61 / −0.41 |
| P53_PATHWAY | −0.16 / 0.57 / −0.25 | +0.07 / 0.57 / +0.07 | −0.05 / 0.56 / −0.26 | −0.17 / 0.52 / −0.52 |
| UNFOLDED_PROTEIN_RESPONSE | +0.01 / 0.67 / +0.09 | −0.06 / 0.52 / +0.02 | −0.09 / 0.61 / −0.21 | −0.04 / 0.65 / +0.36 |
| CHOLESTEROL_HOMEOSTASIS | **+0.39** / 0.51 / +0.02 | −0.05 / 0.54 / +0.09 | +0.04 / 0.51 / +0.05 | +0.02 / 0.61 / −0.03 |
| MTORC1_SIGNALING | −0.02 / 0.57 / −0.14 | +0.10 / 0.51 / −0.06 | −0.10 / 0.63 / +0.18 | −0.03 / 0.55 / −0.13 |
| MYC_TARGETS_V1 | −0.02 / 0.61 / +0.14 | −0.28 / 0.55 / +0.13 | +0.03 / 0.56 / +0.12 | +0.02 / 0.55 / −0.17 |
| INTERFERON_ALPHA_RESPONSE | +0.10 / 0.56 / +0.18 | +0.03 / 0.57 / −0.16 | −0.16 / 0.72 / −0.26 | −0.06 / 0.66 / −0.13 |
| TNFA_SIGNALING_VIA_NFKB | −0.16 / 0.54 / −0.14 | −0.04 / 0.54 / −0.06 | −0.01 / 0.52 / −0.24 | −0.08 / 0.51 / −0.04 |
| HYPOXIA | −0.01 / 0.56 / +0.21 | +0.05 / 0.63 / +0.06 | −0.01 / 0.55 / +0.21 | +0.09 / 0.56 / +0.36 |

- The v⊥ version (set genes zeroed) gives the same picture (H1 cholesterol +0.31); it is in `partP1.json`.
- **Collinearity of H1's 10 v's:** signed median pairwise r −0.08, but median |r| 0.48, and 47% of pairs have |r| > 0.5.
  - Block 1 (P53, cholesterol, mTORC1, IFN-α, TNF-α, hypoxia) correlates +0.66 to +0.94 internally.
  - E2F, UPR and MYC form a second block, anti-correlated with the first.
- **Sterol targets (lookup):** HGNC has no "sterol biosynthesis" gene group.
  - H1 has 7 targets in Hallmark cholesterol (ACAT2 b 2.5, CTNNB1 16.0, FDPS 13.5, **HMGCR 1,333**, **HMGCS1 150**, NSDHL 28.0, TNFRSF12A 0.0), plus DHCR24 and LBR from the nearest HGNC groups (dehydrocholesterol reductases; b 0.0 each).
  - D/E/F has 4 (ATXN2 10.1, **LDLR 15.3**, PDK3 1.3, S100A11 5.7), and none in the HGNC groups.
  - H1's cholesterol r of 0.39 plausibly rests on its few core sterol knockdowns; D/E/F has almost none.

### P2 — pathway-anchored programs, activated by the sources (`tables/partP2_rows.csv`, `partP2.json`): criterion not met
Term = γ Σ_P a_{k,P} v⊥_{c,P} on b ≥ 5 knockdowns' uncalled gated genes. The source activity a uses sources with the scored line's tables out (`srcdelta`: 59 tables, 2,454 target × background pairs, 560 of 564 needed targets).

| line | row 4 | raw γ 0.25 / 0.5 / 1 | centred γ 0.25 / 0.5 / 1 |
|---|---|---|---|
| H1 full | 0.2928 | 0.2807 / 0.2525 / 0.2324 | 0.2849 / 0.2591 / 0.2401 |
| K562 | 0.1532 | 0.1486 / 0.1422 / 0.1302 | 0.1488 / 0.1448 / 0.1392 |
| HCT116 | 0.3153 | 0.2658 / 0.2304 / 0.0479 | 0.2719 / 0.2333 / 0.0462 |
| HEK293T | 0.1695 | 0.1059 / 0.0925 / 0.0637 | 0.0927 / 0.0706 / 0.0349 |
| HCT116final | 0.3651 | 0.3367 / 0.3161 / 0.1786 | 0.3366 / 0.3119 / 0.2460 |

- Both halves pick the base row ("base"); the gain is 0.
- The term is large (RMS 0.03–0.10 log1p over b ≥ 5 knockdowns), and every component falls with γ. H1's MSE is clamped to 0 from γ 0.5.
- **D/E/F activity, median [q10, q90]:** the spread is about half H1's for E2F, G2M, MYC and mTORC1 (for example E2F: H1 −0.008 [−0.045, +0.007] against D/E/F −0.004 [−0.017, +0.010]), and similar elsewhere. 298 of 300 targets have a source.

### P2b — joint ridge slopes, the pathways with P1 r > 0.2 (Jose's follow-up; `tables/partP2b_rows.csv`, `partP2b.json`): criterion not met
- **Setup:**
  - Selected pathway: CHOLESTEROL_HOMEOSTASIS only.
  - The slopes are joint (ridge, λ = 1, all 10 standardised scores). The term uses the centred activity, γ ∈ {0.05, 0.1, 0.15}, and is otherwise as P2.
- **Check:** the rebuilt metacells reproduce P0's separate slopes (max |diff| ≤ 4.7e-7).
- **Joint against separate cholesterol slope:** correlation 0.40 on H1, 0.74 on HCT116, 0.55 on HEK293T and 0.09 on K562.
- **Term RMS:** 0.007–0.015.

| | fit half 0 → half 1 | fit half 1 → half 0 |
|---|---|---|
| pick (grid) | γ 0.15 (edge) | γ 0.05 (edge) |
| held-out gain | −0.0002 | −0.0000 |
| fold pools K562 / HCT116 / HEK293T | +0.0007 / −0.0047 / +0.0001 | +0.0001 / −0.0005 / −0.0002 |
| HCT116final (reported) | −0.0002 | −0.0000 |

## 14.3 Part S — line-vote sign agreement

### S1 — the votes (`tables/partS1.json`, `S1_<context>.npz`)
The C+G rerun reproduces Step 21 C1's mu, var_fe, τ², shrunk LFC and n **bit for bit on all 8 contexts**.

| context | votes 1 / 2 / 3 / ≥4 (gated, own out) | Cochran's Q/df, n ≥ 3: median (q25–q75, q90) | Q/df > 2 | τ² at the floor (n ≥ 2) |
|---|---|---|---|---|
| H1 | 3% / 8% / 18% / 70% | 0.16 (0.04–0.54, 1.08) | 2.6% | 89% |
| K562 | 4% / 11% / 52% / 32% | 0.35 (0.11–0.84, 1.54) | 5.8% | 81% |
| HCT116 | 26% / 45% / 16% / 11% | 0.26 (0.09–0.69, 1.29) | 4.1% | 94% |
| HEK293T | 27% / 44% / 16% / 11% | 0.26 (0.09–0.68, 1.29) | 4.1% | 94% |
| HCT116final | 5% / 26% / 43% / 26% | 0.33 (0.11–0.82, 1.52) | 5.7% | 85% |
| D / E / F | 3–4% / 7–8% / 23–25% / 63–64% | 0.48 (0.19–0.97, 1.64) | 6.3% | 77% |

**Reading:** if the lines agreed within noise, Q/df would centre near 1. It sits at 0.16–0.48, so the lines agree more closely than their SEs imply. The τ² floor is not hiding disagreement; the per-vote SEs are, if anything, inflated.

### S2 — does agreement predict sign accuracy? (`tables/partS2.json`)
Sign accuracy of row 4 on true DE genes (real p_adj < 0.05, gated, own out), by p̃ (votes ≥ 2):

| line | 0.5–0.6 | 0.6–0.7 | 0.7–0.8 | 0.8–0.9 | 0.9–1.0 | 1 vote |
|---|---|---|---|---|---|---|
| H1 | 0.663 (84,835) | 0.689 (83,120) | 0.741 (47,551) | 0.768 (30,382) | 0.934 (2,993) | 0.599 (4,804) |
| K562 | 0.638 (520) | 0.641 (984) | 0.693 (251) | 0.809 (854) | 1.00 (23) | 0.507 (469) |
| HCT116 | 0.589 (538) | 0.722 (425) | 0.812 (579) | 0.867 (315) | 1.00 (11) | 0.593 (920) |
| HEK293T | 0.504 (1,575) | 0.599 (1,025) | 0.648 (1,283) | 0.706 (657) | 1.00 (13) | 0.560 (1,182) |
| HCT116final | 0.663 (692) | 0.743 (576) | 0.790 (699) | 0.902 (624) | 0.882 (17) | 0.563 (87) |

For comparison, Step 21 C5's z_u deciles on H1 run 0.632 → 0.825.

AUC for "sign correct" ("both" = logistic on z_u and p̃; cross-fitted on the H1 halves, all-H1 fit for the others):

| line | all: z_u / p̃ / both | called | uncalled |
|---|---|---|---|
| H1 | 0.572 / 0.557 / 0.575 | 0.604 / 0.619 / 0.632 | 0.551 / 0.536 / 0.552 |
| K562 | 0.635 / 0.598 / 0.635 | 0.672 / 0.675 / 0.696 | 0.528 / 0.505 / 0.514 |
| HCT116 | 0.645 / 0.617 / 0.660 | 0.692 / 0.672 / 0.717 | 0.545 / 0.560 / 0.561 |
| HEK293T | 0.544 / 0.566 / 0.557 | 0.599 / 0.632 / 0.633 | 0.511 / 0.535 / 0.520 |
| HCT116final | 0.699 / 0.647 / 0.709 | 0.734 / 0.681 / 0.754 | 0.573 / 0.566 / 0.580 |

- **By predicted direction:** sign accuracy is higher for genes predicted down than up on every line. Examples: H1 called 0.81 against 0.72, and HEK293T called 0.78 against **0.49**. Called genes predicted up on HEK293T are a coin flip.

### S3 — magnitude × calibrated (2p − 1) (`tables/partS3_rows.csv`, `partS3.json`): criterion not met
Isotonic p_cal per vote class (1, 2, ≥ 3), cross-fitted within the fit half (2 inner folds paired by true-DE rank).

| fit half | pick (grid) | held-out (row 4) | gain | MSE | PDS | NMAE | reach | reach share | held-out terciles weak / mid / strong |
|---|---|---|---|---|---|---|---|---|---|
| 0 | ρ 1, called (edge / edge) | 0.2720 (0.2709) | **+0.0010** | −0.0017 | +0.0096 | −0.0076 | +0.0059 | 97% | +0.049 / 0 / +0.002 |
| 1 | ρ 0.5, both (edge / edge) | 0.3160 (0.3145) | **+0.0014** | −0.0050 | +0.0074 | +0.0057 | +0.0007 | 8% | +0.002 / −0.005 / +0.003 |

- **Fold pools at the picks:** K562 +0.0165 / +0.0094, HCT116 +0.0081 / +0.0115, HEK293T −0.0012 / −0.0042. HCT116final (reported) +0.0100 / +0.0076.
- **Against the criterion:** the gain is below +0.005 on both halves. All components are > 0 and the fold pools are within 0.005.
- **Weights:** at most 0.02% of pairs get m = 0 (calibrated sign no better than a coin flip).

### S4 — reach order (proxy, measure only; `tables/partS4_rows.csv`, `partS4.json`)
Called genes keep row 4's calls and magnitudes; only their proxy p-order changes. Scaled reach (raw):

| order | H1 half 0 | H1 half 1 | K562 | HCT116 | HEK293T |
|---|---|---|---|---|---|
| (i) p_cal | 0.1367 (0.241) | 0.0687 (0.183) | 0.382 | 0.622 | 0.555 |
| (ii) z_u | 0.1297 (0.235) | 0.0670 (0.181) | 0.372 | 0.624 | 0.552 |
| (iii) effect size (proxy) | 0.1359 (0.240) | 0.0848 (0.196) | 0.357 | 0.610 | 0.549 |
| (iv) oracle | **0.2160** (0.308) | **0.1585** (0.259) | 0.456 | 0.660 | 0.584 |
| (i) / (ii) share of the oracle's gain | +1% / −8% | −22% / −24% | 26% / 16% | 25% / 28% | 16% / 9% |

- The oracle's gain on H1 sits in the mid and strong terciles. For example, half 1 mid: −0.011 under (iii), against +0.095 under (iv).
- **Key for R3: p_cal.** It gives the higher held-out reach than z_u on both halves.
- **Note:** on cells, the order comes from the generator's realised Wilcoxon p-values. The frozen generator does not read a requested order: `stage1.prepare_generator_inputs` builds `reach_order` and nothing passes it on. So these are ceilings.

## 14.4 Part Q — two quick checks

### Q1 — why C4 loses the guardrails (`tables/partQ1.json`)
Change against row 4, full panel (scaled):

| line | pick | overall | MSE | PDS | NMAE | fidelity | reach | Jaccard | terciles weak / mid / strong |
|---|---|---|---|---|---|---|---|---|---|
| K562 | κ 4, s 0.75 | −0.0075 | −0.0067 | +0.0049 | −0.0151 | −0.0117 | −0.0165 | 0 | −0.011 / −0.005 / −0.011 |
| K562 | κ 4, s 0.5 | −0.0209 | −0.0290 | −0.0204 | −0.0405 | −0.0113 | −0.0244 | 0 | −0.024 / −0.023 / −0.022 |
| HCT116 | κ 4, s 0.75 | −0.0058 | −0.0124 | +0.0051 | −0.0224 | −0.0046 | −0.0003 | 0 | −0.001 / +0.001 / −0.002 |
| HCT116 | κ 4, s 0.5 | −0.0293 | −0.0578 | −0.0400 | −0.0624 | −0.0018 | −0.0138 | 0 | −0.002 / −0.007 / −0.027 |
| HEK293T | κ 4, s 0.75 | −0.0063 | −0.0039 | −0.0002 | −0.0229 | −0.0013 | −0.0097 | 0 | +0.004 / −0.024 / +0.002 |
| HEK293T | κ 4, s 0.5 | −0.0162 | −0.0162 | −0.0193 | −0.0427 | +0.0003 | −0.0192 | 0 | +0.000 / −0.043 / −0.003 |
| HCT116final | κ 4, s 0.75 | −0.0071 | 0 (clamped) | −0.0127 | −0.0106 | +0.0019 | −0.0214 | 0 | +0.005 / −0.007 / −0.007 |
| HCT116final | κ 4, s 0.5 | −0.0199 | 0 (clamped) | −0.0476 | −0.0504 | +0.0026 | −0.0240 | 0 | −0.021 / −0.009 / −0.027 |

**Anchors** (raw → scaled; baseline anchor, replicate anchor; clamped?):

| line | MSE | PDS | NMAE | fidelity | reach | Jaccard |
|---|---|---|---|---|---|---|
| K562 proxy, row 4 | 0.885 → 0.179 (1.095, −0.076) | 0.727 → 0.597 (0.515, 0.869) | 0.935 → 0.125 (0.994, 0.522) | 0.528 → 0.831 (0.032, 0.628) | 0.487 → 0.357 (0.280, 0.860) | 0.010 → −1.170 (0.243, 0.443) |
| HCT116 proxy, row 4 | 0.881 → 0.963 (1.852, 0.844) | 0.658 → 0.837 (0.520, 0.685) | 0.847 → 0.243 (1.001, 0.369) | 0.513 → 0.866 (0.024, 0.588) | 0.601 → 0.610 (0.296, 0.797) | 0.004 → −1.627 (0.306, 0.491) |
| HEK293T proxy, row 4 | 0.956 → 0.154 (1.141, −0.059) | 0.621 → 0.448 (0.501, 0.768) | 0.878 → 0.197 (0.995, 0.399) | 0.508 → 0.928 (0.035, 0.545) | 0.582 → 0.549 (0.321, 0.796) | 0.026 → −1.259 (0.299, 0.515) |
| HCT116final proxy, row 4 | 0.729 → **1.000 clamped** (1.852, 0.844) | 0.668 → 0.901 | 0.843 → 0.249 | 0.533 → 0.902 | 0.672 → 0.751 | 0.007 → −1.612 |
| Eval A cells, row 4 (Track 1) | 0.787 → 0.187 (0.967, 0.006) | 0.854 → 0.766 (0.494, 0.963) | 0.897 → 0.094 (0.950, 0.382) | 0.447 → 0.469 (0.082, 0.860) | 0.191 → 0.079 (0.124, 0.976) | 0.070 → 0.142 (0.012, 0.420) |
| Eval A cells, live entry | 0.845 → 0.127 | 0.867 → 0.795 | 0.946 → 0.008 | 0.371 → 0.372 | 0.175 → 0.060 | 0.050 → 0.094 |
| HCT116final cells, row 4 | 0.730 → **1.000 clamped** | 0.669 → 0.903 | 0.843 → 0.250 | 0.534 → 0.903 | 0.661 → 0.729 | 0.007 → −1.612 |
| HCT116final cells, live entry | 0.825 → **1.000 clamped** | 0.662 → 0.860 | 0.859 → 0.224 | 0.542 → 0.918 | 0.666 → 0.738 | 0.015 → −1.570 |

**Reading:**
- **What C4's loss consists of.** Shrinking called genes by z_u loses NMAE (−0.011 to −0.062) and reach on every guardrail. It also loses MSE wherever MSE isn't clamped.
- **Why the scaled losses are large.** The fold pools' NMAE and reach anchors are narrow, so small raw losses scale up (NMAE spans 0.47–0.63 raw; the HCT116 reach span is 0.50).
- **Where nothing can move.** HCT116final's MSE is clamped at 1 for every model. Its Jaccard sits far below its baseline anchor (−1.6 scaled) for row 4, the live entry and C4 alike.

### Q2 (M1) — DepMap in the shift (`tables/partQ2_rows.csv`, `partQ2.json`): criterion not met
- **Coverage:** z_k was found for 289 of H1's 292 targets (H1 mean effect −0.234, SD 0.365).

| δ | H1 full (half 0 / half 1) | K562 | HCT116 | HEK293T | HCT116final |
|---|---|---|---|---|---|
| 0 (row 4) | 0.2928 (0.3145 / 0.2709) | 0.1532 | 0.3153 | 0.1695 | 0.3651 |
| 0.1 | 0.2920 (0.3136 / 0.2702) | 0.1521 | 0.3139 | 0.1735 | 0.3645 |
| 0.2 | 0.2903 (0.3117 / 0.2687) | 0.1489 | 0.3077 | 0.1774 | 0.3599 |
| 0.3 | 0.2862 (0.3080 / 0.2642) | 0.1468 | 0.3033 | 0.1776 | 0.3600 |
| M1b (rank-average terciles) | 0.2931 (0.3146 / 0.2714) | 0.1567 | 0.3053 | 0.1682 | 0.3571 |

- Both halves pick δ = 0 ("base").
- **M1b:**
  - Held-out +0.0005 / +0.0001; 89 H1 knockdowns change group.
  - Fold pools: K562 +0.0035, HCT116 **−0.0101**, HEK293T −0.0014; HCT116final −0.0080. It fails on HCT116.

## 14.5 Part R — reach, by ordering the generator's p-values

### R1 — ceiling: S4 above

### R2 — the generator, and the Track 1 cells (`tables/partR2.json`)
**1. Code reading** (`code/cellgen/cell_generator.py`):
- **The spread parameter.** Called genes are shaped by `realize_gene` (143–165). Each starts from a control-shaped base draw (`base_draw`, 96–102), scaled to the target mean, then pulled toward that mean: `_spread(b, mu, c)` = mu + c (y − mu) (105–108). The spread parameter is **c**.
- **How c is chosen.** For a called gene the candidates are c ∈ `SPREADS` = (1, 0.8, 0.6, 0.45, 0.3, 0.2, 0.12, 0.07, 0.03, 0) (line 49, used at 156), then fewer and fewer carrier cells (`_concentrate`). **The first candidate whose test rejects is kept** (161–162).
  - The grid is **global**. The chosen c is **per gene**, found by search.
  - The test threshold p_sig = 0.5 × α × K / M (187) is **per knockdown** (K = its number of calls).
- **Consequence:** a called gene keeps c = 1 (the plain scaled control shape) whenever that already rejects, and is tightened only when it would not. So its realised p is set by its effect size and expression, not by our confidence in it (item 3).
- **Uncalled genes** take the plain scaled control shape (c = 1). If that still rejects, the change goes on the top cells, or is clipped off the top (158).
- **Untested genes** (≤ 5 CPM) carry only the mean (190–193).
- **The final re-check** re-tests every gene on each cell's own library size and rebuilds wrong-side genes once with tighter margins (199–205).

**2. Achievable p by expression** (Eval A cells, row 4's called genes; −log10 p from the harness's `de_pred`):

| control CPM | called pairs | −log10 p: q10 / median / q90 / max | p_adj < 0.05 | control zero fraction | smallest p possible, 400 cells vs 18,400 controls (up / down) |
|---|---|---|---|---|---|
| 5–20 | 10,907 | 2.8 / 4.6 / 9.0 / 72 | 99.0% | 0.80 | < 1e-300 / 1e-22 |
| 20–100 | 97,811 | 2.8 / 4.8 / 27.5 / 62 | 98.6% | 0.41 | 5e-275 / 4e-98 |
| > 100 | 155,685 | 2.8 / 5.0 / 62.8 / 254 | 97.9% | 0.06 | 2e-257 / 1e-229 |

Every band has orders of magnitude of headroom, so even low-count genes can be ordered. The median called gene sits at p ≈ 1e-5.

**3. The realised order:**
- **Median per-knockdown Spearman of realised p rank** (more significant first) within called genes, 276 knockdowns:
  - with |requested LFC| −0.09;
  - with control CPM +0.11;
  - with z_u −0.10;
  - with p_cal +0.06.
- **The order carries no confidence information.**

**Split of the reach gap** (row 4; scaled reach):

| part | proxy (effect-size order) | proxy with the cells' realised order | cells | gap: total = order + rest |
|---|---|---|---|---|
| full | 0.1104 | 0.0781 | 0.0788 | +0.0315 = **+0.0323** + (−0.0008) |
| half 0 | 0.1359 | 0.1117 | 0.0985 | +0.0374 = +0.0242 + 0.0132 |
| half 1 | 0.0848 | 0.0444 | 0.0592 | +0.0256 = +0.0404 + (−0.0148) |

**On the full panel, Track 1's proxy-minus-cells reach gap is entirely the order.** The proxy's signs reproduce the cells' reach once the cells' order is used.

### R3 — generator v8 (`gen_v8/`; `tables/R3_ranks.json`, `partR3_gen.json`, `partR3.json`): **check 3 failed; Part R stopped**
- **Mapping:**
  - Each called gene's spread search starts at c_top = 1 − r and continues down the frozen grid. Here r ∈ [0, 1] is the within-knockdown rank by p_cal, with ties broken by |requested LFC|.
  - The least confident called gene in each knockdown (r = 0) keeps the frozen search.
  - Uncalled genes keep c = 1.
  - **Tightest value: c = 0, the frozen grid's own last value.** At c = 0 every cell sits at the target mean before integerisation, and every c in [0, 1] keeps the mean exactly.
- **Rank arrays saved:** Eval A (p_cal cross-fitted from the other half), HCT116final (all-H1 p_cal), D/E/F (row 4's hybrid T 5 calls from C1's arrays, all-H1 p_cal; 139 of 300 on K4). These are `R3/rank_*.npz`.
- **Slice:** 20 Eval A knockdowns, 10 strong and 10 mid by true tercile (seed 2026), salt "".
  - The 20 were generated three ways (19 min): frozen, v8 without rank, v8 with rank.
  - v8 with rank went into Track 1's file in place of those knockdowns, and that file was scored by the frozen harness (20 min).

| check | result | pass rule | pass |
|---|---|---|---|
| (1) v8 without rank = frozen counts | bit for bit, 20 / 20 (frozen rerun = Track 1 cells, bit for bit) | exact | **yes** |
| (2) pseudobulk log1p CP50K | max \|diff\| 0, p99 0 | — | — |
| (2) MSE / PDS / NMAE, slice-scaled | +0.00002 / 0.0000 / +0.00007 | within 0.002 | **yes** |
| (3) realised called set (p_adj < 0.05, own out), pooled Jaccard | **0.981**; +261 gained, −154 lost | ≥ 0.99 | **no** |
| (3) fidelity / Jaccard on the slice | 0.4237 → 0.4225 / 0.2142 → 0.2144 | — | — |
| (4) Spearman of realised p with requested rank, 5–20 / 20–100 / > 100 CPM | 0.25 / 0.33 / 0.15 (frozen 0.06 / 0.02 / −0.01) | > 0.5 above 20 CPM | no |
| (5) reach on the slice | −0.0898 → −0.0871 (raw 0.0476 → 0.0499) | — | — |

**What moved under check (3):**
- Most knockdowns keep a called-set Jaccard of 0.98–0.995.
- The exceptions are CAST (+84 / −8; 0.921), RRM2 (+17 / −37; 0.957), VCL (0.965) and INSIG1 (0.972).
- **Cause:** tighter spreads change the cells' library sizes. The frozen re-check (199–205) then flips borderline genes on both sides of the BH cut-off, including uncalled ones.
- Tightening does order the p-values (Spearman 0.15–0.33 against about 0), but less than 0.5.
- As specified, I stopped Part R here and did not tune it.

## 14.6 Part D — recommendation (`tables/partD.json`)
- **Candidates:** P2, P2b, S3, Q2 and Q2-M1b all fail their criterion. There is nothing to combine.
- **Arm N2: none.** File 1 stays **row 4**.
  - **Expected Eval A** = row 4's cells, **0.2895** (cell evidence, Step 21 Track 1).
  - The largest proxy gain seen was S3's +0.001 per H1 half. That is below the 0.02–0.05 proxy reach bias and below the criterion.
- **Part R:**
  - **Ceilings:** the oracle adds +0.080 / +0.074 scaled reach on the H1 halves. The chosen key (p_cal) captures +1% / −22% of it, and z_u −8% / −24%.
  - **Criterion (report only): not met.** The oracle passes +0.03; the key does not capture half.
  - **R2:** Track 1's proxy-minus-cells reach gap is all ordering (+0.032 of +0.032 on the full panel).
  - **R3 failed check 3.** **I do not recommend cell-testing v8 now.** The proxy shows no ordering gain from any key we have, and v8 moves the called sets past tolerance.
- **P1, my reading:**
  - The controls do hold missing response, about 4 points of H1's mid + strong residual beyond E4w K 100, 20× random.
  - That response is line-specific (other lines' controls capture less than half as much).
  - The pathway coefficients that would place it are not predictable from the sources (|r| < 0.2 except H1 cholesterol, where few targets drive it), so it is a ceiling we can't yet use.
- **S2, my reading:**
  - Line agreement does predict sign: accuracy climbs from about 0.6 to 0.8–0.9 across the p̃ bins.
  - But it is nearly redundant with z_u (AUC 0.56–0.65 each; combined at most +0.015, about +0.03 on called genes).
  - The lines agree better than their SEs imply, so disagreement isn't the missing signal.

## 14.7 Part E — cell runs
- **E1** (row 4 on v8, pre-approved on condition that R3 passed every check) **did not run.** The queue (`step22_e1.sh`) read `partR3.json`, found check 3 failed, and skipped it (`logs/chain_e1.log` 15:30:43).
- **E2** (Arm N2 on the frozen generator) has nothing to run, because Arm N2 = none.

## 14.8 What I would have done differently
1. **Phase assignment.** With scanpy's sign rule, H1 and HEK293T have no G1 cells, so basis A is empty there. I would have checked the phase proportions on a slice first and proposed a within-line rule (for example the lowest-score tercile as G1) before the run.
2. **P2's term scale.** The ten separate slopes stack collinear programs (median |r| 0.48), and the term's RMS was 0.03–0.10 log1p. I did not compare that with the residual's size before running the grid; a one-minute check would have shown whether γ ≥ 0.25 was in range. Jose's P2b fixed the collinearity, but with one pathway left there was little to activate.
3. **R3's design.** Grading the spread changes each cell's library size, and the frozen re-check then moves borderline calls. A v8 that also re-balanced libraries, or kept the re-check on the frozen library sizes, would have been needed for check 3. Also, S4 had already shown that no key we have orders well on H1. I would have said so before R3 rather than after.
4. **Memory.** R3gen's hybrid file was built in memory from the full Track 1 prediction (29.8 GB) without a logged estimate. It should have been streamed by knockdown, with an estimate logged first.
5. **The write-up.** I started it only after the last result. Drafting it while S and Q ran would have saved about 40 minutes.

# 15. 2026-11-03 — Row 4 for D/E/F: build, checks and submission
*Task: Jose, 2026-11-03: "build out row 4 now; we will submit it tonight". Once the file is ready, run the final checks, submit and document it here. This is Step 21's Part H, File 1.*
- **Code:**
  - `code/vcc6/row4_def.py` (`inputs` / `final` / `check` / `check_vcc`);
- **Outputs:** `results/stage1/2026-11-03_row4_DEF/{D,E,F,final,tables,logs}`.
- **Frozen code, used unchanged:** generator, `run_line.py`, `vcc prep`.

**Submitted: entry `3P7TGy6GVvDqx9slS8GX`.**
- **Upload:** model name LHY, `is_final` true. 3.0 GiB uploaded in 4 min 10 s; exit 0 at 18:12:17 UTC.
- **Server scoring job:** `vcc-metrics-prod-3p7tgy6gvvdqx9sls8gx`.
- **Status:** `launching` until 18:15, `scoring` from 18:16, **`published` at 18:25:57 UTC**, with no error (`logs/vcc_status.log`, `logs/vcc_status_final.json`). Partition test, panel `vcc2026-test-1`, anchors testD/E/F-r1. Scores stay hidden until the winners are announced.
- **The file:** `final/prediction_DEF.vcc`, md5 `27b8e247793080fca78c2d1764283398`, 3,168,481,280 B.
- **It replaces** the 10-30 entry `RTXvB7VtfATH5NQLobww` (config 1′ + K0), which was `published` at the time of this submission. The last entry is the one that counts.

## 15.1 The generator input (`tables/inputs.json`)
**Row 4** = hybrid T 5 calls (K4 1,200 where b ≥ 5, else K0) + config 1′ magnitudes + the no-source fix + the fixed shift + the F5 term.
- **Fixed shift:** (−0.25, 0, +0.25) × C2's `S_DEF` on each context's gate, by (b) tercile, with the all-H1 cuts 3.0 / 67.7.
- **F5 term:** `partF5_DEF_Pprog_E4w_M3DepMap_K100.npz`, md5 `3165cf2d…` (Step 20 F5's D/E/F prediction). It applies on b ≥ 5 knockdowns.
- **Base arrays:** Step 18's config 1′ + K0 inputs (the 10-30 entry's, all 117 tables, each context's own 18,400 controls). The construction is exactly `step21.input_arrays`.

| Check | Result |
|---|---|
| (1) the same builder on HCT116final (screen) reproduces Step 21 Track 1's generator input | lfc, calls, n_sources and K4 assignment identical, bit for bit |
| (2) terms removed: lfc = Step 18's input; calls = K4 where b ≥ 5, K0 elsewhere | yes, on D, E and F |
| (3) the F5 file's P_x = ŝ × E4w-full basis on the gate | max \|diff\| 3e-8 (float32) |

| context | on K4 | (b) terciles weak / mid / strong | requested calls, median (mean) | live entry | entries changed | n_sources fixed |
|---|---|---|---|---|---|---|
| D | 139 / 300 | 129 / 128 / 43 | 808 (751) | 243 (353) | 2.44 M | 29,382 |
| E | 139 | 129 / 128 / 43 | 870 (771) | 279 (381) | 2.49 M | 26,655 |
| F | 139 | 129 / 128 / 43 | 906 (759) | 245 (356) | 2.46 M | 23,525 |

The added term's RMS on gated genes is 0.011–0.012 log1p, and the F5 term's on b ≥ 5 knockdowns 0.010–0.011.

## 15.2 Generation and merge
- **Generation**: `run_line.py predict --lfc-npz`, salt "", 400 cells, 8 workers × 1 thread.
  - D 16:23–16:43, E 16:43–17:02, F 17:02–17:22.
  - Realised against requested calls, median: D 811 / 808, E 875 / 870, F 913 / 906. Median recall 1.0, precision 0.998.
  - ANG (the one uncovered target) has 0 calls requested and realised.
- **Fetch** (28 files, md5 both sides) at 17:22–17:24.
- **Merge** (`row4_def.py final`, 17:24–17:38): straight from the generator outputs, controls stripped, `context` set, `concat_on_disk`. The per-context submission files were temporary and deleted after the merge.
  - Result: **360,000 × 18,533**; D / E / F 120,000 cells each = 300 knockdowns × 400; target lists match; 0 controls; obs names unique; genes = the final `gene_names.csv`.

## 15.3 Final checks
**Before `vcc prep`** (`tables/check_gen.json`, 17:38–17:44): PASS.
- **Generator reports:**
  - Each has 300 rows in list order, with requested calls and n_unmeasured equal to the input. Seeds are unique.
  - The largest mean-count error is 0.00125. Realised minus requested is −39 to +69 calls, median +1.
- **Every file** (three generator outputs, merged, the live entry's three outputs): counts integer, non-negative and finite.
- **Identity:** per-knockdown count sums are exactly equal between each generator output and the merged file (max |diff| 0), with no controls in the merged file.
- **Fidelity to the requested means:**
  - Pseudobulk level r ≥ 0.999999. Δ r median 0.9997 (minimum 0.994, on ANG).
  - The gene-shift check: r at offset 0 is 0.9997; at ±1 it is ≤ 0.01.
  - The controls are identical to the live entry's.
- **Against the live entry's cells (Step 18):**
  - Median per-knockdown Δ correlation is 0.95 (D / E / F 0.950 / 0.955 / 0.952). Median Δ norm ratio, row 4 over live, is 1.03 / 1.05 / 1.05.
  - The minimum (−0.50) is **ANG**. The live entry predicts no change there (no source, Δ = 0, so its cells are noise). Row 4 gives it the weak-tercile shift via the no-source fix: n_sources 1 on its gated genes, about 15.5k genes changed, no calls.
  - That is row 4 as defined, and the construction Track 1 reproduced. Over the inputs, no other knockdown correlates below 0.59 with the live entry (lowest: ATXN2, ZNF526, RBM27).

**`vcc prep`** (17:44–17:59): exit 0.
- Targets verified against the official list; contexts D / E / F at 120,000 each; normalization counts-preserved.
- The .vcc meta matches the live entry's .vcc (CLI 0.3.0, panel `vcc2026-test-1`, the same perts fingerprint) except nnz: 1,942,625,295 against 1,942,114,131.

**After `vcc prep`** (`tables/check_vcc.json`, 18:02–18:06): PASS.
- The decompressed .vcc payload has the same 900 groups in the same order as the merged file, exactly equal per-knockdown count sums, and the same nnz.
- All counts integer, non-negative and finite.

**`vcc whoami`:** phase final, contexts D / E / F, panel `vcc2026-test-1`, deadline Nov 05 23:59 UTC, "ready to submit".

## 15.4 Notes
- **m0 convention:** the builder uses each context's controls-CSV `control_mean_cpm` as m0, as Step 21's row 4 does. It differs from the input file's m0 by CSV round-off only (at most 5e-14 relative); check (1) still reproduces Track 1 bit for bit.
- **`checklist20.py`** is hard-wired to Step 18's files. Its role here was covered by `check` and `check_vcc` (Step 18's `check18.py` logic, applied to this build).

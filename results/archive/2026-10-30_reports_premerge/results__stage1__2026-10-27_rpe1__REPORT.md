# Stage 1 — report (2026-10-27; started as the RPE1 dev run)

## Current state (updated 2026-10-27 14:10 UTC)

The RPE1 dev run (Steps 0–3, below) is finished. Validation has since moved to the D/E/F final target list. The dated log at the end holds every observation and decision.

### Validation setup (since 2026-10-27 11:00)
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

### Baselines on the new setup

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

### Running now
- **HCT116final end to end:** built (14:28); C+G scored 0.190 [0.144, 0.235].
- **Nothing is running** (15:05, Jose: all study, harness and sweep jobs stopped; study outputs archived in `results/archive/2026-10-27/`).

---

# RPE1 dev run (Steps 0–3)

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

## Step 0 — one copy of the model, a fixed list

**(a) The move.** The notebook's model cells went into `code/vcc6/stage1.py` verbatim, and the notebook imports it.
- **Changes:** `knockdown_seed(name, salt="")` (the default reproduces the old seeds), and `calibration_stats` written to `summary["fit"]` for the shrunk and the unshrunk predictions.
- **`calibration_stats`:** this function did not exist anywhere on vcc2026: not in the deleted `stage1_model.py`, the notebook or the bytecode. It is new. It reports Pearson r, slope of true on predicted (calibration) and slope of predicted on true (transmission), pooled over each knockdown's tested genes.
- **Self-test:** passes, including new checks for the salt, `calibration_stats` and the Step 2 arms.
- **Identity:** `predictions.npz` on 20 RPE1 knockdowns is **byte-identical** before and after the move (`step0/identity.json`; numpy writes fixed zip timestamps). The "after" run used the later, refactored `stage1.py` (arms, votes, `load_line`), so arm A is unchanged by everything since.

**(b) Archive.** `stage1_model.py` had already been deleted from `~/vcc`.
- My snapshot of it, and its stale `.pyc`, are in `archive/2026-10-26/` with a README.
- `stage1_rpe1.ipynb` did not exist on vcc2026.
- The only file that mentioned either was `vcc_2026_6.ipynb`, in markdown, and it has since been rewritten.

**(c) The list.** All 1,968 knockdowns in the `replogle_rpe1` table have kept cells, so the 300 were drawn from all 1,968. Per knockdown (`step1/knockdowns_300.csv`):
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

**The eval** (`harness/lib/evalRPE1dev`) was built with the harness's own `build_cells` and `finish`. `code/vcc6/rpe1_eval.py` only points `common.LIB` and `score.RES` at this folder.

---

## Step 1 — oracle diagnosis

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

**Per knockdown, arm C** (`step1/armC_per_knockdown.csv`), in scorer space over the genes RPE1 measures, own gene excluded:
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

## Step 2 — how sources vote (effect space, the 300)

**Structure** (`step2/backgrounds.csv`, `tables.csv`, `structure.json`):
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

## Step 3 — end to end (`run_line.py`: controls only → model → frozen generator → frozen harness)

**Controls-only profile against the sidecar** (`step3/profile_check.json`):
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
- **Run time:** about 25 minutes per 300-knockdown arm per salt to generate, and about 11 minutes to score, on vcc2026 alone.

**`vcc prep` correction.** The first Step 3 runs called `vcc prep` on the harness file and it failed every time. I had not checked its exit code.
- **The format needs two things:** no control cells in the file, and a `context` column covering A, B and C in one file.
- **Fix:** `run_line.py` now writes `prediction.h5ad` (scored by the harness, controls kept) and, with `--prep` only, `submission.h5ad` (no controls; `--context A|B|C`, or an A/B/C split of a dev line's knockdowns as a format-only check).
- **Result:** G's file passes (120,000 cells, 18,533 genes, counts preserved).
- **Opt-in:** `vcc prep` does not affect the scores and is only for files to upload, as the plan says.

---

## Magnitude test (no refit; arrays from the Step 3 arm B model run)

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

## Recommendation

**Adopt arm C with the Gaussian shrink (decided).** It is the first stage-1 model above the generic baseline, and the gain is far outside the salt noise.
- **Weakest metrics:** lfc NMAE (−0.58) and direction reach (−0.24), together with the weak tercile's MSE.
- **What would move them:** a calibrated magnitude (G is still about 1.5× too large, calibration slope 0.67), and a gate on knockdowns that do nothing in the target line.

## What I would have done differently

- **Tested the premise first.** A ten-minute oracle-only run before building Step 1 would have shown the oracle scores 2.2, and redirected Step 1 before it was built.
- **Acted on the first sign of the shrinkage problem.** Step 0's 20-knockdown run already showed π at its floor and NMAE 0.9997. The magnitude test belonged straight after Step 0, not after Step 3.
- **Checked every exit code.** I only checked `vcc prep`'s after Step 3. Its failures sat in 12 run records unnoticed.
- **Changed one thing at a time between Steps 1 and 3.** Steps 1 and 3 differ in two ways at once: library size (20k against the line's own about 8.5k) and the controls used (the sidecar's 11,309 against the context half's 5,654). So their oracle difference (2.24 against 2.10, Jaccard 0.54 against 0.34) can't be attributed to either.
- **Harmonised background names before Step 2,** since one-vote-per-background and same-line exclusion both depend on them.
- **Run the two lanes on two machines.** Running Step 2 alongside Step 1 on 4 CPUs slowed generation by about 2×. vcc2027 would have halved the wall time.

---

## Study plan (agreed 2026-10-27)

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

## Study plan v2 (2026-10-27, after the external research brief)

The brief changes the frame more than the method. **The final deadline is Nov 5, 23:59 UTC (9 days).** The final test lines are D, E and F, released Oct 21–22. Two submissions are allowed per UTC day, only the last counts, and scores stay hidden until Nov 19.

Two gaps on this machine:
- **The CLI is too old.** The `vcc` CLI is **0.2.2**; the final phase needs **0.3.0 or newer**.
- **The final data isn't here.** `controls-final` has not been downloaded, so `data/` holds only the validation A/B/C contexts and `pert_counts.csv`.

Without these, nothing we build can be submitted. The plan therefore gains a submission track ahead of the model steps.

### What the brief confirms in our setup
- **Scorer:** cell-eval2 vcc2026, with the six metrics as implemented, scaled from baseline 0 to replicate 1. MSE is a ratio of sums with the noise credit capped, and pds is cosine with all 300 targets removed. Our harness runs that package (0.18.0) unchanged.
- **The MSE baseline sits at raw 0.986–0.992, and strong entries reach raw 0.70–0.85.** C+G is at raw 0.867 on RPE1, though RPE1 isn't D/E/F.
- **Tahoe is low value** (CC0, drug data; useful mostly for basal states). That matches Step 2.
- **PIE and text embeddings add little.** That matches the decision to drop PIE.
- **Shrinkage and amplitude can't be tuned on D/E/F;** they must come from cross-validation or a rule learned on source lines. That is exactly Step 1, arms E and M.

### Rules risk in our current stack (needs Jose's decision)
1. **The frozen generator shapes significance.** For every gene in the predicted `de_sig`, it pulls the control-shaped cells toward the mean until the Wilcoxon test rejects. The brief lists "DE-significance shaping: narrow or reshape cell spread at a fixed mean so chosen genes test significant" as the **highest rules risk** of all, unanswered by staff. It is also the source of our Fidelity and Jaccard gains.
2. **Building cells from the target line's control distribution** (resampled or quantile-based NTC cells): grey, no staff ruling. Peers say "almost everybody" does it.
3. **Precision-weighted pooling of published deltas, plus shrinkage.** The 9/17 rule allows training or fine-tuning on published results, but not pasting them in or revising predictions with them. A learned model (fitted scale and weights, cross-validated) is on the safer side. A plain weighted average of public deltas is closer to "pasting". Arms E and M, and any learned weighting, move us toward "learned".
- **Suggested action:** one written question to help@virtualcellchallenge.org covering 1–3. Keep a fallback generator mode that doesn't shape significance ready. That has to be a new module with its own oracle test, because `cell_generator.py` is frozen.

### Ordered plan

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

### Update 2026-10-27 08:10 UTC

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

## Log of observations and decisions (UTC)

Kept in time order. "Decision" means Jose decided it; "Observation" is a measured result; "Action" is something done.

### 2026-10-26 — source tables
- **Decision:** Jiang genes with two symbols on one Ensembl ID keep the panel's symbol and drop the alias; never sum them.
- **Decision:** Jiang's blank Ensembl IDs are filled from the panel by symbol (617 genes per file).
- **Action:** `write_replogle_qc` now maps `''` to `unresolved:<symbol>` for new contexts. QC had pooled 2–4 blank-ID targets into one knockdown in all 30 Jiang sources. QC and the 30 tables were redone, all passing.
- **Observation:** Tahoe 50/50 and Jiang 30/30 tables built and passing. `qc_effects_extra.py`: 80 tables, 0 failing.
- **Decision:** a source never votes on a target of its own cell line (`obs.background`), whatever its stem.
- **Decision:** Tahoe drugs that share a target are separate observations (×16), never averaged.
- **Decision:** small Tahoe tables are kept.
- **Action:** `stage1_model.py` moved into `vcc_2026_6.ipynb`. The source list was checked: 117 of 117 tables. Orion pseudobulk falls to the ×8 default.

### 2026-10-27 — stage-1 RPE1 dev run (Steps 0–3)
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
- **Observation:** Step 1's first launch was OOM-killed at 08:21 while loading tables next to a cell-eval2 scoring. It was re-queued behind the scoring chain and started at 08:49. Lesson: never overlap a full table load with scoring on vcc2026.
- **Decision (09:15):** vcc2027 may be used when it suits the job.
- **Action (09:24–09:40):** effect tables, code and lists synced to vcc2027 through `gs://vcc2026/sync_stage1/`. The final-target screen was queued there, then cancelled on Jose's request (no runs until the order is agreed). vcc2027 was stopped, then restarted for the Orion stream.
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
  - vcc2026 and vcc2027 both discover 117 tables, with identical Orion checksums.
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
- **Check (11:58):** the final-screen C+G baseline (vcc2027, 11:00–11:11) used the rebuilt single-cell Orion tables: `orion_hct116` 434 and `orion_hek293t` 418 knockdowns, with md5 identical to vcc2026's v2 build (synced at 10:55). Non-cell tables are kept (the 5a proposal).
- **Action (12:00): single β and arms E and M on the final screen** (`step_scale.py --set final`, vcc2027).
  - **β (the reference):** C+G × β, where β = Σxy / Σx² over the pooled held-out source-line pairs (no intercept, the target's background excluded).
  - **E and M:** compared against β and against C+G.
  - **Calls:** every variant keeps C+G's call set (Jose: E's magnitude with G's calls), so sign and AUROC are unchanged for β and E.
  - **Scorer-space pb_nmse:** added per variant (`screen.extra_metrics(lfc_all=…)`).
  - **Output:** one json per line in `results/stage1/step1_scale/final/`.
  - **Memory:** one process holds 14.8 GB after loading the tables, so the four lines run in sequence, not in parallel.
  - **Smoke test** (H1, 3 calibration knockdowns per line): β = 0.27; H1 NMAE 1.41 → 1.03, pb_nmse 1.60 → 1.01 (mean response 0.92).
- **Status (12:00):** vcc2026 is in the HCT116final build's `prep-real-bundle` step (since 11:05); K562final is queued behind it.
- **Observation (12:30): salts vs knockdown resampling** (RPE1dev C+G, three salts, per-knockdown `results.csv`).
  - Per-knockdown values are identical across salts for 8 of 9 metrics: per-knockdown r = 1.000 for pds, MSE, distance, nmae, Jaccard and yield.
  - Only `direction_reach` moves: salt range 0.006 raw, per-knockdown r = 0.991. That metric accounts for the whole 0.0018 scaled spread (0.1216 / 0.1234 / 0.1220).
  - The bootstrap SD over knockdowns is 3–100× the salt range: pds 0.014, yield 0.012, reach 0.019, nmae 0.005, MSE 0.0008.
  - Implication: one salt plus a paired bootstrap over knockdowns measures the noise that matters. Three salts only measure generator noise, which is negligible. Pending Jose.
- **Observation (12:30): an H1 end-to-end eval already exists.** Harness Eval A (`lib/harness/evalA`, built 10-25): H1 with the whole line held out, 292 knockdowns, 400 cells each, 18,400 controls. That is the D/E/F format exactly, and it includes all 28 H1 final targets. HCT116final and K562final are both weak-response lines (median 1–3 true DE genes); H1 is Arc's own data and strong (272). Proposal: Eval A is the primary end-to-end eval and HCT116final the weak-regime check; drop the queued K562final. Pending Jose.
- **Decisions (Jose, 12:40):**
  - **K562final cancelled** (its queued waiter was killed; the build never started).
  - **End-to-end pair:** H1 (harness Eval A, strong response) and HCT116final (weak response), always reported side by side.
  - **No salt:** end-to-end runs use salt "" only. Uncertainty comes from a bootstrap over knockdowns.
  - **Step 1** goes ahead on the rebuilt final screen: single β, then E and M on top of C+G.
- **Action (12:45): `code/vcc6/boot_score.py`.** Resamples cell-eval2's per-knockdown `results.csv` (2,000 resamples, seed 2026) and rebuilds the aggregate as cell-eval2 does (means; the MSE ratio of sums). It then scales with cell-eval2's own `score_metrics` against the eval's bundle, with the anchors fixed. Runs on the same eval are compared paired.
  - The full-set rebuild reproduces the harness's `summary.json` exactly (checked on every run).
  - On RPE1 C+G: overall 0.1216, 95% CI [0.088, 0.157]. Paired salt s1 − none = +0.0018 [+0.0005, +0.0031]. So the generator's draw is a small but consistent shift, about a tenth of the noise from which knockdowns are scored.
- **Action (12:45): C+G end-to-end baseline queued on vcc2026** (`logs/e2e_final/run_CG.sh`), behind the HCT116final build:
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
- **Action (13:03):** results copied to vcc2026 (`results/stage1/step1_scale/final/`); vcc2027 stopped (idle).
- **Observation (13:30): R² on the final screen.** Figures in `results/stage1/step1_scale/final/figures/` (`plot_step_scale.py`).
  - **Gene-level r², out of sample** (pooled over the final targets; the same for C+G, β and E): K562 0.0006, HCT116 0.0002, HEK293T 0.0002, H1 0.0098. M: 0.0009 / 0.0002 / 0.0003 / 0.012.
  - **In sample** (held-out source lines): 0.0005–0.0007 for every variant. There is no gap between in and out of sample, so no overfitting; but there is almost no fit to overfit.
  - **Noise ceiling** (1 − mean se² / var(true lfc), same genes; `r2_ceiling.json`): K562 0.136, HCT116 0.123, HEK293T 0.143, H1 0.606. We reach 0.1–1.6% of the r² the truth's noise allows. H1 does best (1.6%).
  - **R² with the scale counted** (1 − MSE / var): C+G is −0.11 / −0.06 / −0.20 / −0.45; C+G × β is −0.005 / −0.005 / −0.011 / +0.006. Shrinking moves us onto "predict nothing" and no further.
- **Decision (Jose, 14:00):** set aside the H1 end-to-end run and finish HCT116final as soon as possible. Add the "Current state" section at the top of this report. Start the shrinkage and source-dependence study (3 folds); vcc2027 may be used.
- **Action (14:05):** the C+G end-to-end queue was stopped (its H1 step never ran). The HCT116 C+G prediction runs on vcc2027 (8 CPUs) while the bundle finishes on vcc2026; scoring follows on vcc2026.
- **Observation (14:03): vcc2027 held a stale, pre-freeze `cell_generator.py`** (sha256 0ba855a1…, 10-26 10:26), not the frozen 422f8a3d…. It also had an older `build_effects_v5.py`.
  - Neither was used by today's vcc2027 runs: the screens and Step 1 are effect space only, with no generator.
  - Both were replaced from vcc2026; the stale copies are kept in `archive/2026-10-27_vcc2027_stale/`. The generator is now read-only and hash-identical on both VMs.
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
- **Action (14:25):** HCT116 C+G prediction done on vcc2027 (278 knockdowns, median 129 calls realised of 132.5 requested). It was copied to vcc2026, and scoring is queued behind the HCT116final bundle (`logs/e2e_final/score_HCT116_CG.sh`).
- **Action (14:26): shrinkage study build started on vcc2027** (`code/vcc6/shrink_study.py build`). It does one load of the 117 tables and writes:
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
- **Check (14:40):** vcc2027 already had a copy of Eval A from 10-25; copying vcc2026's context over it was harmless, since both builds have identical split files and build records.
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
- **Action (15:00): magnitude sweep started** (`code/vcc6/magnitude_sweep.py`; `logs/sweep_*.sh`).
  - **Estimate:** 12 runs (11 new; HCT116 m = 1 is the 14:48 run), about 4.5–5 h wall, ending around 20:00 UTC.
  - **Method:** calls are recomputed on m × lfc by the model's own rule. At m = 1 the recomputation reproduces the model's calls exactly. HCT116 mean calls per knockdown are 196 / 15 / 1.3 / 0.2 / 0 / 0 for m = 1 / 0.5 / 0.25 / 0.14 / 0.07 / 0.03.
  - **Where it runs:** H1 runs are scored on vcc2027, whose Eval A is byte-identical to vcc2026's (all 18 files' md5); HCT116 runs are scored on vcc2026 after the study's 4 scorings.
- **Action (15:02):** study analyze re-run on vcc2027 with the screen metrics and the B collapse check added.
- **Decision (Jose, 15:05): stop everything and archive the study (move, never delete).**
  - **Stopped**, both VMs confirmed idle:
    - vcc2026: the study's copy-and-score queue (`study_score.sh`; nothing had been scored yet) and the sweep's HCT116 score queue (`sweep_score_hct116.sh`; waiting).
    - vcc2027: the study chain (`study_chain.sh`), mid-generation of HCT116 baseline cells (about 210/278 knockdowns, never scored); the analyze re-run with the added screen metrics (3 min in, not finished); the sweep chain (`sweep_chain.sh`) and H1 score queue (`sweep_score_h1.sh`), both still waiting, nothing generated.
  - **Moved to `results/archive/2026-10-27/`, on both VMs:**
    - `2026-10-27_shrinkage_study/`: 2.4 MB on vcc2026 (report, tables, figures, results.json, config); 2.4 GB on vcc2027 (plus `cache/` and the partial `e2e/`).
    - `step1_scale/` (arms β / E / M, v1 and final): 1.3 GB on vcc2026, 1.6 GB on vcc2027.
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
  - Details and the memory workarounds (model and generation in separate processes; an on-disk merge; temporary swap) are in `results/stage1/2026-10-28_DEF/report.md`. Final scores are hidden until the winners are announced.
- **Scorer-space study (2026-10-28), consolidated** in `results/stage1/2026-10-28_scorer_space/report.md`.
  - The metric is reproduced exactly. The noise credit is not the gap; the plug-in distance is. A noise-free proxy matches end-to-end MSE to 0.0006.
  - **Recommendation pending Jose:** C+G LFC × 0.5 (m* = 0.5 on every fold). Eval A MSE 0 → +0.12, PDS −0.03; expected overall about +0.015 on Eval A.

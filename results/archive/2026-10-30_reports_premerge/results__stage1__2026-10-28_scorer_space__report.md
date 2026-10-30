# Scorer-space MSE study (pseudobulk only) — report

## Consolidated findings (Steps 0–8, 2026-10-28; detail in the step sections below)

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
- Built the D/E/F file with an on-disk merge.
- Kept in mind that all arm numbers are noise-free with Ĉ_pred = 1.13 × C0, and the budget binds for arm (a) at m ≤ 0.5, so only the end-to-end run settles the recommendation.

**Next, after approval:** arm (a) m = 0.5 against the incumbent, end to end on Eval A and HCT116final, paired bootstrap, all six components. One job per VM.

---

Task: Jose, 2026-10-28 ~06:05 UTC. Rules: CPU and pseudobulk only through Step 5; no cell generation before Jose approves Step 6; never alongside scoring on one VM.

**Scheduling note:** the D/E/F final file (C+G, Jose's choice, to be uploaded) was generating on vcc2026 when this task arrived. It held about 26 GB of 31 GB, leaving 4 GB free. Step 0's code reading was done at once. Everything that loads cells or pseudobulks waits until the D/E/F file is prepped and uploaded, so the submission is never put at risk.

## Step 0 — the MSE path in the installed cell_eval2 0.18.0 (code reading, 06:05–06:15)
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

### `~/h1_decomp`
- **Not on vcc2026.** Memory places it on vcc2027 (`~/h1_decomp.py`, run 10-20), which is off. Its memory note says it scored "noise-free pseudobulks with the scorer's own formulas", so it would not have needed the credit or the budget.
- **What we use instead:** vcc2026 has the 10-21 MSE study, `analysis/mse_metric/mse_study.py` (with `scripts/v5/pb.py`). It **does implement both bounds** through cell_eval2's own internals:
  - per row, `min(C_pred, C_real)` (#247);
  - the #348 budget via `_across_pert_budget` / `_row_weights`, with `ratio = min(1, budget / claim)` (`mse_study.py:141–153`).
- This study builds on that code rather than starting vcc2027 to read `h1_decomp.py`.
- **Reproducing the HCT116 C+G raw MSE (2.896177):** pending. It needs the prediction (1.1 GB) and the reference loaded, which waits for the D/E/F job.

### How C+G passes non-called genes to the generator
**Kept at their predicted LFC.** `stage1.prepare_generator_inputs` zeroes only genes no source measured (n_sources = 0); those are never called. Every other gene carries the shrunk posterior mean; `de_sig` only marks which genes get shaped to test significant.

### How the generator turns LFCs into counts (`code/cellgen/cell_generator.py`, frozen, sha256 422f8a3d…)
- **Means:** mean CPM per gene = (control_mean_cpm + 1e-9) × 2^lfc − 1e-9 (`target_means`, line 168). This is the v5 LFC inverted, so the generated cells' mean CPM equals the requested one up to rounding.
- **Counts:** counts = round(cpm × library_size / 1e6), with largest-remainder rounding, so each gene's total count is exact.
  - `run_line` passes `library_size` = the controls' median panel library (`run_line.py:148`), e.g. 12,133 for HCT116final and 20,236 for context D.
- **Per-cell totals are not renormalized.** Each cell is built at `library_size` × (Σ requested CPM / 1e6). The generator's own note: "the scorer renormalizes each cell, so the realized means are the requested ones times 1e6 / sum(means)", and the realised composition is recorded in `uns['cellgen']['composition']`.
  - So when a prediction's means don't sum to 1e6, every gene's realised level shifts by the same factor after the scorer's per-cell (and per-group, CP50K) normalisation.
  - That shared shift is the "renormalization offset" Step 4 measures.

### Reproduction (06:44, `code/vcc6/scorer_space.py check`; `tables/step0_reproduce_HCT116.json`)
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

## Decision before Steps 1–5 (Jose, ~06:50)
**K562 and HEK293T get pseudobulk real sides built from their QC'd cells**, since only H1 (Eval A) and HCT116final have harness bundles.
- **Cells:** the harness's own build step (`final_evals.py build K562final|HEK293Tfinal --cells-only`): seeded 400 cells per knockdown, 20k-depth thinning, context/scoring control halves, on each line's final list (K562 267, HEK293T 271). No DE bundle.
- **Real side:** pseudobulk and jackknife in scorer space, computed with cell-eval2's own functions.
- **Anchors:** cell-eval2's own baseline and replicate-anchor code if it can be restricted to the MSE metric (no DE). Otherwise my own computation, labelled as such.
- **Scheduling:** runs after the D/E/F upload (the K562 cell file is 6.6 GB).

## Amended task (Jose, 2026-10-28 ~08:25): Steps 1–7 on vcc2027, vcc2026 stays with the D/E/F upload
vcc2027 was started for this; its stale Orion QC'd cells (10-18, A–C only) were moved to `archive/2026-10-28_vcc2027_stale_qc/` and replaced by vcc2026's 10-27 v2 files (md5-checked). Every MSE number goes through `analysis/mse_metric/mse_study.py:score_generated`, with `code/vcc6/scorer_study.py:EvalPB` as its eval object (pb.py's attribute set on any harness eval; anchors and PDS via cell-eval2's own `_replicate_entries` / `score_one` / `discrimination_score`).

## Step 1 — reproduction through mse_study.py (`tables/step1_reproduce.json`)

| | HCT116 C+G | Eval A C+G |
|---|---|---|
| raw MSE, rebuilt vs harness | 2.8961771873171123 vs 2.896177187317113 (diff 8.9e-16) | 1.4961872607668685 vs 1.4961872607668685 (diff 0) |
| PDS raw, rebuilt vs harness | 0.62784718 vs 0.62784718 (diff 1.1e-16) | 0.83353340 (diff 0) |
| knockdowns with Ĉ_pred > C_real (#247 cap) | 16 / 278 | **249 / 292** |
| budget B / claim (#348) | 1.29 (does not bind) | 3.91 (does not bind) |
| baseline / replicate raw MSE | 1.852 / 0.844 | 0.967 / 0.0055 |
| scaled MSE (clamped / unclamped) | 0 / −1.036 | 0 / −0.550 |

On HCT116 the credit actually granted, Σ min(Ĉ_pred, C_real)/|G|, is 6.3× the summed denominator: HCT116's real signal is tiny next to its sampling noise. On Eval A it is 0.27×.

## Step 2 — Δ per source (`scorer_study.py step2`, vcc2027; `delta/<table>.npz`, `tables/step2_*`)
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

## Step 7 — blocked for the chemistry families (08:42)
**No chemistry is recorded in the metadata:** not in the 117 effect tables' `uns`/`obs`, and not in any `lib/qc_cells/*_qc.h5ad` (`uns` searched recursively for chem/tech/platform/assay/kit/protocol/flex/library_type).
- The only hint is Jiang's depth note, "Parse split-pool".
- Per the rules this part stops rather than guess. It needs a source → chemistry mapping from Jose (e.g. from the publications).
- The Arc-first sub-arm (Arc H1's weight × w_arc) needs no chemistry and will run once Step 5 has the arm (c) m*.

## Progress (08:42)
vcc2027 chain:
1. Step 4 (noise credit).
2. K562/HEK293T real caches + MSE/PDS-only anchors (HCT116final's exact cell-eval2 commands, metric list narrowed).
3. Per-fold C+G + Δ pool (one table load).
4. Steps 3, 3b, 5, 6.

The K562final cells are built (267 knockdowns, 48,020 cells, 18,400 + 18,400 controls); HEK293Tfinal is building.

## Step 4 — noise-credit accounting (`scorer_study.py step4`; `tables/step4.json`, `tables/step4_per_knockdown.parquet`)
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

## K562 and HEK293T reference sides (vcc2027, 08:36–08:56)
- **Cells:** `final_evals.py build K562final|HEK293Tfinal --cells-only`, the harness's own build step. Seeded 400 cells per knockdown, 20k-depth thinning (factor 1.0 for both), context and scoring control halves.
  - K562final: 267 knockdowns, 48,020 real cells, 18,400 + 18,400 controls.
  - HEK293Tfinal: 271 knockdowns.
- **Real side and anchors:** cell-eval2 itself (`logs/scorer_refs.sh`): HCT116final's exact `prep-cache` / `baseline` / `prep-real-bundle` commands (vcc2026 preset, pdex, counts, 5 anchor splits). The only change is `--set metrics=[expr_mse_unbiased_capped, expr_distance_unbiased, expr_mse_unbiased_capped_norm, pds_cosine]`, so no DE is computed (cache 48 s instead of hours). Every exit code was 0.

| Fold line | Baseline raw MSE | Replicate raw MSE | Baseline PDS | Replicate PDS | Source |
|---|---|---|---|---|---|
| K562 | 1.095 | −0.076 | 0.515 | 0.869 | built today (MSE/PDS only) |
| HCT116 | 1.852 | 0.844 | — | — | HCT116final bundle |
| HEK293T | 1.141 | −0.059 | 0.501 | 0.768 | built today (MSE/PDS only) |
| H1 (Eval A) | 0.967 | 0.0055 | — | — | Eval A bundle |

## Per-fold models (`scorer_study.py folds`, vcc2027, 08:56–09:05; `cache/fold_<line>.npz`)
- **One table load for all folds.** Per fold, C+G is rebuilt exactly as `run_line`/`predict_line` builds it: arm C votes, SE curve, prior, Gaussian shrink, on each eval's context controls, with the fold's lines' tables dropped. The Δ pool uses the same prior, SE curve and votes, with each source's LFC replaced by its Δ and duplicates collapsed as `load_source` does.
- **Identity check:** the H1 fold's C+G equals the incumbent Eval A C+G model (`2026-10-27_incumbent/evalA/CG/model`) exactly: max |ΔLFC| 0, calls identical, measurement SE identical.

| Fold | Tables voting / with Δ | Prior τ² / v0 | Line: median calls, share of (k, g) with a Δ source |
|---|---|---|---|
| K562 | 108 / 103 | 1e-4 / 0.783 | K562: 113, 0.79 |
| HCT116+HEK293T | 113 / 106 | 1e-4 / 0.499 | HCT116: 461, 0.52; HEK293T: 592, 0.52 |
| H1 | 116 / 109 | 1e-4 / 0.541 | H1: 486, 0.81 |

## Step 3 — where the MSE error sits (`scorer_arms.py step3`, `step3b`; `tables/step3.json`, `tables/step3b.json`)
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

By target control-CPM decile (gated genes, low → high), the error is spread fairly evenly. It peaks in the middle deciles for HCT116/HEK293T/H1 and in the low deciles for K562 (`tables/step3.json`).

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

## Step 5 — magnitude, four arms (`scorer_arms.py step5`; `tables/step5_grid.parquet`, `tables/step5.json`)
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

## Step 6 — renormalisation offset (`scorer_arms.py step6`; `tables/step6.json`)
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

## Step 7 — Arc-first sub-arm (`scorer_study.py arcw`, `scorer_arms.py step7arc`; `tables/step7_arc.*`)
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

## Step 8 — recommendation (stop for Jose's approval)
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
3. **Built the D/E/F file with an on-disk merge.** The in-memory combine hit 30 GB plus 15 GB swap and had to be stopped.
4. **One caveat on the evidence:** every number here is noise-free with Ĉ_pred = 1.13 × C0. The budget B binds for arm (a) at m ≤ 0.5, so the realised credit in generated cells could differ. Only the end-to-end run settles it.

**Stopped here, waiting for Jose's approval.**

---

# Part 2 — Steps 9–15 (task: Jose, 2026-10-28 ~17:30 UTC)
**Rules:** vcc2027 for Steps 9–14 (noise-free proxy, existing cells only); vcc2026 free for Step 15 after approval; never alongside scoring on one VM; report after each step.
**Goal:** Eval A MSE scaled > 0, every other component > 0, overall above C+G (0.219). HCT116 and HEK293T are guardrails only.

**Infrastructure:**
- vcc2027 started 17:38. The first start attempt hit a zone STOCKOUT; the second worked.
- **K562 and HEK293T needed DE references for the four DE proxies.** Their Steps 3–8 references were MSE/PDS-only.
- They are being rebuilt with the harness's own `finish` (real cache, ceiling, generic baseline, 5 anchor splits) on the existing cells: `final_evals.py build <X>final --finish-only`.
- The MSE/PDS-only cache, baseline and bundle moved to `eval<X>final/msepds_2026-10-28/`. Nothing was deleted.
- I stopped the first attempt (17:41–17:49) after 8 min, so no table load ran alongside it. It restarts after the Step 10 fold builds (`logs/step10_folds.sh`).

## Step 9 — proxies for all six components (`scorer_proxy.py step9`; `tables/step9*.{json,csv,parquet}`)
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

## Step 10 — source filters × magnitude arms (`scorer_part2.py step10`, `step10fit`; `tables/step10*.{json,csv,parquet}`)
Run on vcc2027: H1 and HCT116 at 17:50–18:05, then K562 and HEK293T at 21:13–21:16, after their DE references were built (K562final 18:00–19:40, HEK293Tfinal 19:40–21:13, `logs/de_refs.log`). All exit codes were 0.

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

The full six-component rows for every line are in `tables/step10_heldout.csv`.

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

**Side check: arm (e) past the grid edge** (`scorer_part2.py step10ext`, `step10extfit`; `tables/step10ext*.{json,parquet}`, vcc2027 21:20–21:24). Non-called m ∈ {0.3, 0.5, 0.75, 1}, F0, all m_c. This is outside the task's grid, so it does not enter the carry rule.

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


## Step 11 — adding signal on the carried config (`scorer_s11.py`, `scorer_step11.py`; `tables/step11*.{json,parquet}`)
Run on vcc2027: H1 and HCT116 at 18:42–18:51, on the provisional carry that Step 10's fit then confirmed. K562 and HEK293T ran 21:16–21:23. The first `step11fit` failed (exit 1) because the HGNC table was only on vcc2026. I copied it over (md5 checked) and reran the fit at 21:24 (`logs/step11_12_resume.sh`); it exited 0.

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

### S1 — strength
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

### S2 — PCs of the source lines' Δ (fold-excluded)
Energy share of the line's true Δ captured by K PCs:

| Line | Source PCs, K 5 / 10 / 20 / 50 | Own top K (ceiling), K 50 | Random (floor), K 50 |
|---|---|---|---|
| K562 | 0.016 / 0.027 / 0.036 / 0.047 | 0.43 | 0.006 |
| HCT116 | 0.035 / 0.050 / 0.066 / 0.080 | 0.44 | 0.005 |
| HEK293T | 0.065 / 0.091 / 0.115 / 0.135 | 0.52 | 0.005 |
| H1 | 0.043 / 0.087 / 0.125 / 0.213 | 0.75 | 0.005 |

- **The shared subspace is real, but it holds only a fifth of H1's energy at K 50.** That is 40× random but well under H1's own ceiling.
- **Arm (g):** LOFO picks K 50 (the top of the grid) on every fold, with m 0.5 / 0.2 / 0.3. Projecting first lowers H1 to 0.208 and drops reach to 0.011. Projection throws away the per-gene signal that reach and NMAE need.

### S3 — complex partners (sources only)
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

  HGNC groups shared within the top 30 are listed in `tables/step11.json` → `h1_top30_shared_groups`. Only METTL3/METTL14 are both in the top 30 and each other's partners.

### S4 — per-(k, g) confidence
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

## Step 12 — call rule on the carried magnitudes (`scorer_step12.py step12`, `step12fit`; `tables/step12*.{json,csv,parquet}`)
Run on vcc2027 21:24–21:27, on the four lines in parallel. Every variant uses arm (i)'s magnitudes, so MSE, PDS and NMAE are the same in each row. K1 at m_call 1 reproduces K0's call set exactly on all four lines (0 mismatches).

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

## Step 13 — generator check (`scorer_step12.py step13`; `tables/step13.json`) — measurement only, run ahead of Steps 10–12 because it needs no carried config
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

## Step 14 — recommendation (stop for Jose's approval)
vcc2027 was stopped at 21:31 UTC once every table and log had been copied to vcc2026 (status TERMINATED). The code is identical on both VMs (md5 checked).

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
4. **Stage every input on the run VM before chaining.** The HGNC table was missing on vcc2027, which cost one restart of Step 11's fit (about 1 minute; no results lost).

**Step 15 plan, after approval:**
- Config 1 on vcc2026 and config 2 on vcc2027, one config per VM, in parallel.
- Each runs on Eval A and HCT116final (screen setup, HEK293T voting) against the incumbent C+G scores.
- Report: paired bootstrap of all six components, scaled and raw; realised vs requested calls; the adoption rule.
- Stop and report before building any D/E/F file.

# Step 16 — wider magnitude grid and H1 knockdown-split call budgets (task: Jose, 2026-10-29)
**Order of work.** Jose's note put Step 16 on vcc2027 after Step 15, but Step 15 had not been approved or started. Asked at 07:30 UTC, Jose chose Step 16 first.

vcc2027 runs (`scorer_step16.py`, `logs/step16_*.{sh,log}`):
- **Stage A:** 07:32–08:09.
- **Fit:** the first run failed on my bug in `spreads()`, which was passed the whole fits dict. I fixed it and reran at 08:09; no rows were lost.
- **Stage B:** 08:10 onwards.

Noise-free proxy only; no cells were generated.

## A — magnitude (`step16 --line`, `step16fit`; `tables/step16{.json,_heldout.csv,_rows_*.parquet,_curve_*.json,_zbeta_*.json}`)
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

## B — calls (`step16calls`; `tables/step16_calls{,_config1,_config1_ext}.json`, `…_rows*.parquet`)
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

## Step 16 recommendation (stop for Jose's approval before any cells are generated)
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

Config 2 is dropped as superseded. That pair measures the realised call gain on real cells directly, as one config per VM.

**What I would have done differently:**
- **Put s_c in the LOFO grid from Step 10 on.** The fold disagreement on the called scale is what blocks transfer, and it was visible only once s_c moved.
- **Run budgets out to "all tested genes" at once.** The response is monotone, so a short grid only produced an edge.
- **Check the NMAE margin under knockdown resampling before recommending config 1 in Step 14.** A half split takes seconds and would have shown +0.001 to be fragile.
- **Unit-test the fit summary before launching.** The `spreads()` bug cost one rerun of the fit (2 minutes; no rows lost).

# Step 15 — proxy checks, knockdown strength, end to end (task: Jose, 2026-10-30)
Outputs (tables, logs) in `results/stage1/2026-10-29_step15/`. Code: `code/vcc6/step15.py` (new module), run from `vcc_2026_8.ipynb`'s Step 15 section. The frozen generator, the harness and config 1's frozen files are imported and not changed.

**Rules:**
- 15a and 15b: noise-free proxy only, no cells, on vcc2027.
- 15c: cells, one config per VM, never next to a table load or a scoring on the same VM; only after Jose's go.
- No D/E/F file, no upload.

**Definitions used:**
- **Config 1′** (F0): called genes 0.5 × C+G LFC; non-called 0.6 × c × Δ_pool, c = z_a²/(z_a² + 1); single-source c₁ = 1; no source → control. That is config 1 with m 0.5 → 0.6.
- **K0:** C+G's calls. **K4 1,200:** per knockdown, the top 1,200 tested genes by C+G z. The magnitude split (LFC side vs Δ side) stays K0's; only the call set changes, as in Step 16 B.
- **H1 halves:** `scorer_step16.halves` (seed 2026, pairs by rank of true DE count). Halves are scored with the full panel's scaling, as `boot_score` does.

## 15b.1 — direction reach in cell-eval2 0.18.0
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

### Addendum to 15b.1 — lfc_nmae, sig Jaccard, and p-value order in reach
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

## 15a — proxy check of the planned configs (vcc2027, 03:02–03:10 UTC; `results/stage1/2026-10-29_step15/tables/step15a.csv`, `results/stage1/2026-10-29_step15/logs/step15a.log`)
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

## 15b.2 — oracle ceiling (vcc2027, 03:10–03:21; `results/stage1/2026-10-29_step15/tables/grid_oracle_p*.parquet`, `results/stage1/2026-10-29_step15/tables/step15b.json`)
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
- **λ = 2 is on the grid edge**, but λ 1 → 2 is worth only +0.004–0.008, so the curve is flattening.
- **The magnitude gain is not usable as is.** O1 wants the strong tercile at 1.5× (grid edge) and the weak and mid at 0.5×. That raises MSE (+0.067) and PDS (+0.015), but **pushes NMAE to −0.0002** (−0.0054 on one direction) and costs reach −0.015. Under the adoption rule ("every component > 0"), O1 fails even with perfect strength.

## 15b.3 — strength predictors (vcc2027, 03:02–03:04; `results/stage1/2026-10-29_step15/tables/step15b_predictors{.json,_H1.csv,_K562.csv}`)
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

### Addendum O3 — the model's own expected DE count (`results/stage1/2026-10-29_step15/tables/grid_o3_fit*_p*.parquet`, `step15b_o3_counts.json`)
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

## 15b.4 — predictor (b) in place of true strength (`results/stage1/2026-10-29_step15/tables/grid_pred_b_fit*_p*.parquet`, `grid_o3p_b_*`, `grid_o3s_b_*`)
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

**Reading of 15b.4:**
- **A source-only strength predictor recovers little within the task's grid:** at most +0.005 (O1 + O3s), with every component still positive.
- **Its magnitude half (O1) recovers 37% of a small oracle gain** (+0.004 of +0.010) and keeps NMAE positive (+0.005). The fit wants the *mid* tercile at 1.5×, not the strong one, which suggests (b)'s terciles don't line up with true strength at the top.
- **Its call half helps only at K4's volume or above (side check, +0.013).** At the oracle's preferred volume it loses, because a wrong budget costs more than a flat one.

## Stop — summary and recommendation for 15c (2026-10-30 04:15 UTC)
vcc2027 ran 03:01–04:09 UTC and is **stopped** (TERMINATED). Every table and log is copied to vcc2026 (md5 checked). No cells were generated, no D/E/F file was built, nothing was uploaded.

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
- **vcc2026: config 1′ + K0. vcc2027: config 1′ + K4 1,200.**
- **Why not swap in a strength variant:**
  - Their in-grid proxy gains (+0.003 to +0.005) are below the paired-bootstrap resolution we've seen on Eval A (CI half-widths around 0.01–0.015).
  - The larger +0.013 comes from an off-grid λ at an edge, with a D/E/F count scale that can't be checked.
- **What the pair settles on cells:**
  - (1) Does config 1′'s thin NMAE (+0.0014 on half 1) survive generation?
  - (2) Does the proxy's +0.024 for 1,200 flat calls survive the generator's spread-shaping (Step 9 validated the proxy only at ~500 calls)?
- If K4 1,200 holds up on cells, (b)-allocated budgets at the same volume are the natural next proxy step. They would add calls only where the sources say the knockdown is strong.
- **Not recommended:** O1-style magnitude scaling by strength. Even the oracle version turns NMAE negative.

**15c runtime estimate (after your go), per VM, the two VMs in parallel:**

| Stage | vcc2026 (4 CPUs) | vcc2027 (8 CPUs, 8 workers) |
|---|---|---|
| Table load + model, both lines (v7 cells) | ~13 min (freeze check: 9 min load + 1.5 min per line) | ~13 min |
| Generate Eval A (292 × 400 cells) | ~56 min (incumbent run) | ~30–35 min |
| Generate HCT116final (278) | ~45 min | ~25–30 min |
| Score Eval A (harness) | ~50 min | ~40–50 min |
| Score HCT116final | ~20 min | ~20 min |
| boot_score ×2 + call / SD diagnostics (items 2–4, terciles) | ~15 min | ~15 min |
| **Total** | **≈ 3.3 h** | **≈ 2.5–3 h** |

- K4's 1,200 shaped genes per knockdown may slow generation; I haven't measured that, so add up to 30 min on vcc2027.
- Steps run strictly in sequence on each VM: the model's table load finishes and its process exits before generation, and no scoring overlaps a table load.
- vcc2027 can hit a STOCKOUT again on start (this morning's start took 6 tries, about 6 minutes).
- vcc2027 needs the incumbent C+G harness runs and the Eval A / HCT116final references. I'll check those before launching and copy any that are missing (md5 checked).

**What I would have done differently:**
1. **Put λ past 2 in the O2 grid from the start.** λ = 2 was an edge for the oracle and for (b), and the useful result for (b) is at λ 6. A proxy row costs seconds; I repeated Step 16's lesson here.
2. **Matched call volume when comparing allocation rules.** O2 with predicted counts at λ 2 mixed "where the calls go" with "how many calls", and lost on volume. The side check had to separate the two afterwards.
3. **Checked a predictor's own count against the truth before running its grid.** O3's N_k has ρ 0.15–0.29 with the truth, which predicts its result without 250 proxy rows.
4. **Set thread limits from the start.** The first oracle batch ran 6 multithreaded processes at load ~38 on 8 vCPUs. The single-threaded chains afterwards were clearly faster per row.

**Waiting for Jose's go before 15c.**

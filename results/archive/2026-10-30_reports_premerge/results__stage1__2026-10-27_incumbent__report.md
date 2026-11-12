# D/E/F incumbent: Model 2 (effects_v5, no floor) vs C+G on Eval A — report

Task (Jose, 2026-10-27 15:10): port the leaderboard Model 2 recipe to effects_v5 without the detection floor, score
it against C+G end to end on Eval A, recommend one, then (after Jose's choice) prepare the D/E/F file. vcc2026 only;
nothing uploaded.

## Step 1 — the Model 2 recipe (2026-10-27 15:20)

### Two corrections to the task's premises
1. **No entry scored 0.145.** +0.145 is the leaderboard *gain* from v2 to v3 (−0.054 → +0.091).
   - The best entry is **v3: +0.091**, rank 749, submitted 2026-10-13. It is Model 2 with top-2,000 calls and a detection floor, from `vcc_2026_3`.
   - The `vcc_2026_2` recipe itself (v2, gated Model 2) scored −0.054.
2. **The leaderboard Model 2 never read effects_v4.**
   - effects_v4 was built on 10-18, after every Model 2 entry (v2 10-1x, v3 10-13, v3a 10-14, v4 10-15).
   - v3 read the panel source tables `lib/sources/<line>/effects_panel.h5ad` (from `scripts/build_panel_sources.py`: 575 targets × 18,533 genes, layers lfc / se).
   - **Those panel tables no longer exist:** archived to the bucket on 10-19, and deleted from the bucket on 10-21 (user-approved cleanup; `ARCHIVE_MANIFEST.csv` rows now point to deleted objects).
   - The effects_v4 tables do survive in `gs://vcc2026/vm/lib/effects_v4/`, but no submission used them.

### The recipe (v3, entry rank 749)
- **Code:** notebook `vcc_2026_3.ipynb` at commit `c208334` ("the vcc_2026_3 that ranked 749"). Extracted to `recipe/`: the notebook, `vcc2_meta.py` (Model 2), `vcc3_generate.py` (call rule, floor, generator) and `vcc2_common.py`, exactly as written by that commit. The current `code/vcc_2026_3c/` copies were edited later.
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

### Leaderboard entries (overall / pds / nmae / fidelity / reach / Jaccard; MSE always 0)

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

### What this means for Steps 2–3 (needs Jose)
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

### Decisions after Step 1 (Jose, 15:30)
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

### Step 1, addendum: how C+G hands non-called genes to the generator
**Kept at their predicted LFC.**
- `stage1.prepare_generator_inputs` (called by `gen_cells` for every knockdown) zeroes only genes **no source measured** (n_sources = 0); those are never called.
- Every other gene keeps the shrunk posterior mean, called or not. `de_sig` only marks which genes the generator shapes to test significant.
- `p_sign` only orders a saved `reach_order` that the frozen generator never reads.
- So generating C+G through `run_line --lfc-npz` from its own arrays is the same as its model path. Both arms use `--lfc-npz`.

## Step 2 — the port (`code/vcc6/model2_port.py`)

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

### Step 2 — run (15:35–15:46, exit 0 for both jobs)
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

## Step 3 — pre-generation comparison (Eval A, 292 knockdowns, H1 held out; `evalA/step3_compare.json`)

| | Model 2 port | C+G |
|---|---|---|
| Calls per knockdown, median / mean (truth: 218 / 1,048) | 249.5 / 388.6 | 486.5 / 551.4 |
| Median \|LFC\| of called genes | 0.459 | 0.315 |
| Median \|LFC\|, all measured genes | 0.075 | 0.062 |

- **Per-knockdown Pearson r between the two models' LFCs:** median 0.285 (IQR 0.210–0.364) over genes any source measured; 0.353 (0.262–0.488) over measured genes in the 5-CPM gate.
- **Call-set Jaccard:** median 0.085 (IQR 0.049–0.177).
- **Not nearly identical** (r ≪ 0.95, call counts differ about 2×), so Step 4 proceeds.
- The port calls about half as many genes at larger |LFC|. The two pools differ in sources: 7 lines vs 116 tables, including Feng's 20 iPSC lines, KOLF, Jiang, Adamson, three Marson conditions and Tahoe. They also differ in τ² (0.0030 vs 1e-4) and in shrinkage (none vs Gaussian).

## Step 4 — end to end on Eval A
**Runtime estimate (stated 15:52): about 2.5 h on vcc2026, one job at a time.** Generate 2 × ~25 min, then score 2 × ~45 min, then a ~10 min paired bootstrap. Expected done around 18:30 UTC.
- **Revised estimate (16:03):** generation runs at about 10 s per knockdown on vcc2026's 4 CPUs (50 of 292 in 590 s), so about 55 min per arm. That makes Step 4 about 3.5 h in total, done around 19:20 UTC.
- **Step 4 stopped (17:55): scoring the Model 2 port on Eval A was killed out of memory.** cell-eval2 reached 27 GB RSS on vcc2026 (31 GB, no swap; VS Code holds about 4 GB). Both arms' cells are generated and kept (exit 0: Model 2 58 min, median 249.5 requested / 271 realised calls; C+G 59 min). Nothing was scored.
- **Action (2026-10-28 03:35 UTC):** 16 GB swap added on vcc2026 (`/swapfile_eval`, not persistent; to be removed after scoring), Jose's choice. The killed run's partial output was moved to `harness_oom_partial/`, and scoring restarted (`logs/step4_score.sh`).

### Step 4 — results (Eval A, 292 knockdowns, H1 held out, salt "", frozen generator → frozen harness → `boot_score.py`)
Scoring ran with 16 GB swap: peak about 28 GB RSS plus 5.5 GB in swap. Model 2 port 03:35–04:27, C+G 04:27–05:18, bootstrap 05:18–05:21; every exit code 0.

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

## Step 5 — recommendation
**C+G for D/E/F.** Its paired advantage over the Model 2 port on Eval A is +0.034, 95% CI [+0.020, +0.048], which excludes 0. It is ahead on every component.
- It also has end-to-end results on RPE1 (+0.122) and HCT116final (0.190 [0.144, 0.235]).
- **Caveats:**
  - H1 has misranked candidates before (v4 and v5 on the leaderboard). Both of those were call-everything recipes against a top-N one. Here the two arms share the call rule and generator and differ only in pooling, the kind of difference H1 ranked correctly in v2 → v3.
  - Model 2's leaderboard evidence (v3, +0.091) was earned with a different generator, top-2,000 calls and a detection floor. None of those are in either arm here.

**What I would have done differently:**
1. **Checked the recipe's history before porting.** The task's premises (an entry at 0.145; effects_v4 as the input) took a step to correct. The top-N port was started under the old task and stopped mid-load.
2. **Checked memory before scoring Eval A on vcc2026.** HCT116final (55k cells) scored there, but Eval A (135k cells) needs about 33 GB. The first scoring was killed after 58 min of generation, and the run then waited overnight for a decision.
3. **Separated pooling from source set.** The port pools 7 lines; C+G pools 116 tables (Feng, KOLF, Jiang, Adamson, three Marson conditions, Tahoe). Their predictions correlate only 0.29, so this comparison cannot say whether C+G wins on its random-effects shrinkage or on its extra sources. A C+G arm restricted to Model 2's 7 lines would separate the two. It would only matter if the margin needed explaining, not for the choice.

**Stopped here, waiting for Jose's choice before Step 6** (the D/E/F file).
- **Swap removed (2026-10-28 ~05:25 UTC):** `swapoff` and `rm /swapfile_eval`, exit 0; vcc2026 is back to no swap.

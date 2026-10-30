# Shrinkage and source-dependence study (3 folds) — report

Status (2026-10-27 14:25 UTC): Step 0 inventory done; build running on vcc2027 (`code/vcc6/shrink_study.py build`), then
`analyze`. Results sections below are filled in as the steps finish.

Code: `code/vcc6/shrink_study.py` (build: one table load → cache) and `code/vcc6/shrink_analyze.py` (Steps 0.4–5).
Outputs: this folder: `tables/*.parquet`, `figures/*.png`, `results.json`, `config_final.json`, `stageC_<line>.npz`;
`cache/` (per-target arrays and baseline predictions; copied to `gs://vcc2026/shrinkage_study/`).

## Decisions taken before the build (Jose, 2026-10-27 14:15)

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

## Changes before analyze (Jose, 2026-10-27 14:35)

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
- **Does Tahoe drive the similarity?** No (`code/vcc6/similarity_check.py`; `tables/similarity_tahoe_check.{parquet,json}`). The HVG set was recomputed over the 67 non-Tahoe tables only:
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

## Step 0 — inventory

### 0.1 Where things are
| Item | Location |
|---|---|
| Transfer model | `code/vcc6/stage1.py` (`predict_line`, arm C, magnitude `gaussian`) |
| Control-cell similarity | **None in the codebase** → the spec's default (see decisions) |
| Effect tables | `lib/effects_v5/` and `lib/effects_v5_extra/`: 117 tables (`stage1.discover_tables`), each with LFC and SE layers, a `_baseline.csv` (control profile, `de_tested`) and a control-quantile sidecar |
| Fold definitions | `code/vcc6/screen.py` `FINAL` (line → truth table) and `results/target_lists/final_<line>.txt` (D/E/F targets ∩ table) |
| Eval harness | `scripts/harness/` (frozen), driven by `code/vcc6/run_line.py score`; evals `lib/harness/evalA` (H1) and `results/stage1/final_evals/harness/lib/evalHCT116final` |

### How the baseline sets its weights (C+G)
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

## Summary

**On the final targets, none of the four predictors (baseline, A, B, C) beats predicting zero on LFC MSE by more than 0.7%.**
- Stage C has the lowest MSE on every held-out line, but it reaches "predict nothing" and no further: −0.08% (K562), −0.08% (HCT116), −0.10% (HEK293T), −0.72% (H1).
- The baseline is 6–49% worse than zero. Almost all of that is removed by **shrinkage** (Stage B ≈ ×0.3, then c ≈ ×0.14–0.17), not by better weights.
- **Stage A's correlation-aware, similarity-modelled weights make MSE slightly worse** than the baseline's on all four lines.
- **The reason is that the sources barely agree with each other or with the target on these knockdowns:** source–source error correlation net of the target is ρ_b ≈ 0.001–0.006 and ρ_w ≈ 0.04. Per-source r with the target is 0.001–0.07.
- **The transfer-variance model (3b) is fitted on 11–22 source–target pairs drawn from about 8 tables per fold, and its similarity coefficient changes sign across folds.** The final config's 3b is not reliable for unseen lines.
- **Harness metrics** (baseline vs C on HCT116final and H1 Eval A) are pending.

## Direct answers

| Question | Answer |
|---|---|
| **Is the pooled slope well below 1?** | Yes, far below. y on p through the origin: K562 0.053, HCT116 0.034, HEK293T 0.018, H1 0.113. Huber is lower still (0.036 / 0.022 / 0.016 / 0.063). OLS intercepts are −0.022 / −0.034 / −0.023 / −0.004. |
| **Does the slope rise with source count?** | In K562 yes, monotonically: 0.017 (1 source) → 0.029 → 0.030 → 0.066 (5–9) → 0.073 (10–29, 13 knockdowns). In H1 broadly (0.05 → 0.09 → 0.03 → 0.15). In HCT116 and HEK293T no clear trend (0.01–0.05). Their 10+ buckets come from only 14 knockdowns and are ≈ 0. |
| **Linear or S-shaped?** | **Linear on all four lines.** The through-origin slopes below and above the 90th percentile of \|p\| have overlapping knockdown-bootstrap CIs: K562 0.049 vs 0.055, HCT116 0.038 vs 0.031, HEK293T 0.025 vs 0.014, H1 0.109 vs 0.115. |
| **How dependent are the sources?** | By 3c, barely: ρ_w = 0.050 / 0.036 / 0.040 and ρ_b = 0.006 / 0.001 / 0.004 (folds K562 / HCT116-293T / H1), so n_eff ≈ S (S = 30 → n_eff 25–29). Marson's three conditions (ρ ≈ 0.30) are the only correlated group. The 3d plateau implies ρ = 0.65 / 0.87 / 0.35 / 0.23 (K562 / HCT116 / HEK293T / H1), but on only 13 / 14 / 14 / 1 knockdowns, and it measures something else (see "Undermines the design"). |
| **Does transfer variance track similarity, cell count or dataset?** | Not reliably. 3b's R² is 0.72 / 0.79 / 0.99, but on 19 / 11 / 22 rows from about 8 tables, with 6–9 rows floored. The similarity coefficient is −36.7 / −34.3 / +13.1 (sign flip). Cell count (−4.4 / −4.0 / −3.7) and feng (−4.7 / −2.6 / −4.3) are stable in sign, but the cell-count effect is carried by Marson: imputed ~400 cells, and t floored because its SE² exceeds its whole error. Raw t against cell count shows no trend apart from Marson. |
| **How variable is β_k, and does agreement predict it?** | Very variable around a small mean. β_k p5–p95: K562 −0.04…0.23 (median 0.035), HCT116 −0.03…0.21 (0.022), HEK293T −0.03…0.12 (0.016), H1 0.006…0.40 (0.047). Agreement predicts it weakly: Spearman β_k vs τ̂²_k = 0.27 / 0.28 / 0.25 / 0.04, and B_k vs β_k = 0.31 / 0.24 / 0.21 / 0.03. β_k falls slightly with strength (−0.19 / −0.16 / −0.13 / −0.06). |

## Step 0 — results

### 0.4 Coverage (sources = every table except the fold's lines)

| Fold line | Knockdowns listed | With truth and ≥ 1 source | (k, g) with 0 sources | 1 | 2 | 3–4 | 5–9 | 10–29 | 30+ |
|---|---|---|---|---|---|---|---|---|---|
| K562 | 267 | 267 | 24,659 | 26,070 | 70,605 | 319,415 | 1,514,071 | 95,778 | 0 |
| HCT116 | 278 | 276 | 55,482 | 91,790 | 70,056 | 1,906,437 | 618,480 | 116,286 | 10,161 |
| HEK293T | 271 | 269 | 62,107 | 85,565 | 80,897 | 1,992,348 | 645,448 | 121,527 | 10,733 |
| H1 | 28 | 28 | 1,100 | 2,970 | 5,824 | 51,059 | 231,006 | 9,573 | 0 |

Most (k, g) have 3–9 sources: Marson ×3 plus K562 / H1 / Orion / KOLF / Feng genome-wide. RPE1, HepG2 and Jurkat hold none of the final targets, and Tahoe holds one (MAPK7).

### 0.5 SE distribution (per study, medians of per-table quantiles; `tables/se_distribution.parquet` per table)
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

## Step 1 — baseline (pooled LFC MSE over the target's pairs)

| Line | MSE | strong | mid | weak | Predict zero |
|---|---|---|---|---|---|
| K562 | 0.1599 | 0.1489 | 0.1184 | 0.2124 | 0.1439 |
| HCT116 | 0.3986 | 0.3604 | 0.5078 | 0.3272 | 0.3742 |
| HEK293T | 0.1804 | 0.2120 | 0.1748 | 0.1548 | 0.1549 |
| H1 | 0.0514 | 0.0587 | 0.0452 | 0.0505 | 0.0345 |

Terciles are by strength_k, the RMS of the baseline p on the target, fixed per target. Harness metrics are below.

## Step 2 — calibration diagnostics (figures `step2a_hexbin_*`, `step2b_binned_*`, `step2e_4B_*`)

| Line | OLS₀ | Huber₀ | OLS slope / intercept | Slope \|p\| ≤ p90 / > p90 | Shape | Top-10% knockdowns' share of Σp² |
|---|---|---|---|---|---|---|
| K562 | 0.053 | 0.036 | 0.050 / −0.022 | 0.049 / 0.055 | linear | 26% |
| HCT116 | 0.034 | 0.022 | 0.028 / −0.034 | 0.038 / 0.031 | linear | 19% |
| HEK293T | 0.018 | 0.016 | 0.015 / −0.023 | 0.025 / 0.014 | linear | 20% |
| H1 | 0.113 | 0.063 | 0.114 / −0.004 | 0.109 / 0.115 | linear | 18% |

2c (slope by |S_kg| bucket) is in `tables/step2c_buckets_*.parquet` and in the direct answers above. The negative intercepts mean the targets' mean LFC over these genes is below zero; the baseline has no such offset.

## Step 3 — source error statistics

**3a.** Few sources qualify: at least 20,000 shared pairs and at least 10 shared knockdowns *with a training target's final list*. That is 19 / 11 / 22 source–target rows per fold, from Marson ×3, Arc H1, KOLF2.1J, Feng genome-wide ×2, Orion ×2 and Replogle K562.
- **Scale:** per-source through-origin slopes of y on x are tiny, 0.003–0.06, and r is 0.001–0.07. No source transfers with a usable scale.
- **Marson has t < 0 in every row:** mean se² (0.07–0.23) exceeds V. Its replicate-table SEs overstate its error. All Marson rows are floored, which accounts for 6 / 6 / 9 of the floored rows.
- Tables: `tables/step3a_sources_<fold>.parquet`.

**3b** (`figures/step3b_<fold>.png`):

| Fold | n | R² | intercept | similarity | log cells | feng | Residual mean by target |
|---|---|---|---|---|---|---|---|
| K562 | 19 | 0.72 | 52.0 (19.0) | −36.7 (17.4) | −4.39 (0.83) | −4.65 (1.47) | H1 −0.30, HCT116 +0.15, HEK293T +0.07 |
| HCT116/293T | 11 | 0.79 | 46.5 (21.0) | −34.3 (21.4) | −4.01 (0.97) | −2.62 (2.38) | H1 −0.21, K562 +0.12 |
| H1 | 22 | 0.99 | 1.7 (4.1) | +13.1 (3.5) | −3.72 (0.19) | −4.25 (0.26) | ≈ 0 |
| Final (Step 5) | 33 | 0.68 | 28.6 | −13.0 | −4.10 | −4.85 | — |

The fit sees only "other" and "feng" rows. Jiang and Tahoe get their t̂ at the reference level (Jose 14:35; flagged).

**3c** (`figures/step3c_<fold>.png`, `tables/step3c_pairs_<fold>.parquet`):

| Fold | Source pairs | ρ_w | ρ_b | Within study | n_eff at S = 2 / 5 / 10 / 30 |
|---|---|---|---|---|---|
| K562 | 473 | 0.050 | 0.006 | feng 0.043, jiang 0.027, marson 0.304, orion 0.009 | 1.99 / 4.88 / 9.47 / 25.4 |
| HCT116/293T | 225 | 0.036 | 0.001 | feng 0.024, jiang 0.017, marson 0.304 | 2.00 / 4.99 / 9.94 / 29.4 |
| H1 | 771 | 0.040 | 0.004 | feng 0.036, jiang 0.022, marson 0.296, orion 0.013 | 1.99 / 4.92 / 9.65 / 26.9 |

No ρ was undefined or outside [−1, 1].

**3d** (`figures/step3d_averaging.png`), equal-weight averaging over knockdowns with at least 10 sources:

| Line | Knockdowns | MSE at S = 1 → 10 | A (plateau) | B | Implied ρ = A / (A + B) |
|---|---|---|---|---|---|
| K562 | 13 | 0.277 → 0.189 | 0.181 | 0.096 | 0.65 |
| HCT116 | 14 | 0.625 → 0.548 | 0.555 | 0.085 | 0.87 |
| HEK293T | 14 | 0.269 → 0.105 | 0.097 | 0.180 | 0.35 |
| H1 | 1 | 0.231 → 0.062 | 0.054 | 0.184 | 0.23 |

## Step 4 — revised predictor (pooled LFC MSE on F; tercile columns strong / mid / weak)

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

Single-source knockdowns on F: K562 2, HCT116 9, HEK293T 8, H1 0. B_k vs the oracle β_k: Spearman 0.31 / 0.24 / 0.21 / 0.03 (`figures/step2e_4B_*`). After C, the binned conditional mean (`figures/step4_binned_C_*`) is flat around zero with a small positive slope.

## Harness metrics (baseline and Stage C)

(pending: 4 runs, generated on vcc2027 and scored on vcc2026)

## Stability across folds

| Parameter | K562 | HCT116/293T | H1 | Final (Step 5) | Reliable for unseen lines? |
|---|---|---|---|---|---|
| c | 0.141 | 0.136 | 0.174 | **0.256** | Stable across folds. **The final c is not comparable:** in Step 5 each target excludes only its own line, so HCT116 is predicted with HEK293T (same lab) as a source and vice versa. That raises the agreement, and with it c. An unseen line has no sibling, so the fold value (~0.14–0.17) is the relevant one. |
| ρ_w / ρ_b | 0.050 / 0.006 | 0.036 / 0.001 | 0.040 / 0.004 | 0.042 / 0.002 | Stable, and ≈ independence. |
| 3b similarity | −36.7 | −34.3 | +13.1 | −13.0 | **No:** the sign flips, SEs are large, and n is small. |
| 3b log cells | −4.39 | −4.01 | −3.72 | −4.10 | Stable, but carried by the imputed and floored Marson rows. |
| 3b feng | −4.65 | −2.62 | −4.25 | −4.85 | Two tables (Feng genome-wide, one pooled line). |
| Decile medians of B | 0.39–0.52 | 0.37–0.55 | 0.26–0.47 | 0.26–0.39 | Moderately stable. |

## Recommendation for the 3 unseen lines

(pending harness; effect space alone says the following)
- **In effect space, C is the best of the four on every fold**, but only by matching predict-zero.
- **Its gain over the baseline is entirely shrinkage:** B_k (≈ 0.3) × c (≈ 0.15).
- **Stage A's similarity-modelled weights should not be used as built.** They are fitted from an unstable 3b and are worse than the baseline's weights on all four lines.

## What I would have done differently / what undermines the design
1. **"Final lists everywhere" starves Step 3.** Only about 8 tables share at least 10 final-target knockdowns with a training target. Feng's 19 targeted lines, Jiang, Tahoe, Nadig and RPE1 never enter 3a/3b, so 3b is a 4-parameter fit on 11–22 rows. Fitting Step 3 on the training targets' *full* tables (K562 9,245 knockdowns, H1 300, Orion 434/418) would make 3b estimable. The metrics would stay on the final lists.
2. **3b's cell-count effect is a Marson artefact.** Marson has no cells (its count is imputed at ~400), and its reported SEs exceed its total error, so its t is negative and floored. The cell-count coefficient mostly separates Marson from the rest.
3. **The 3c formula is not an error correlation.** Its numerator, mean(x_s x_s') − ½(mean(y x_s) + mean(y x_s')), and V_s = mean(x_s²) − mean(y x_s) both omit mean(y²) and the doubled cross term. When y is barely predictable (r ≈ 0.01), ρ reduces to the correlation between the two sources' x, i.e. source–source agreement.
   - The 3d plateau, by contrast, includes everything common to all sources: the target's own measurement noise and its idiosyncratic response.
   - So the two disagree (ρ_b 0.004 against an implied 0.23–0.87) because they measure different things, not because of a bug. With so little shared signal, the plateau is mostly the truth's own variance.
4. **One fold is a single study:** HCT116/293T is Orion only, and as a training target it is one lab's assay.
5. **The 3d estimates are thin:** few knockdowns reach 10 sources (13 / 14 / 14 / 1).
6. **Scale differences are not the problem; agreement is.** Per-source slopes are all small and similar (0.003–0.06), so per-source scale correction (3a) would not help.

---

## Task (as given, 2026-10-27 14:00)

# Task: Shrinkage and source-dependence study for the zero-shot transfer model (3-fold)

## Goal
Measure how well our current precision-weighted transfer model is calibrated, how dependent the source cell lines are, and whether three changes reduce LFC MSE:
(A) correlation-aware, similarity-modelled source weights
(B) per-knockdown shrinkage from source agreement
(C) one global calibration constant
Evaluate with 3 held-out folds: K562, HCT116/293T, H1. Then refit on all three for the unseen lines.

## Folds and leakage rules
- Folds are K562, HCT116/293T and H1, using the fold definitions already in the codebase. If a fold contains more than one line, every line in it is treated as the held-out line together.
- In fold F, the held-out line's perturbation data is used nowhere: not as a target, and not as a source in any parameter fit. Only its control cells are used (for similarity).
- In fold F, parameters are fit using the other two fold lines as training targets. When a training target U is the target, the source pool is all lines except U and F.
- When F itself is predicted, the source pool is all lines except F.
- RPE1 is an ordinary source line throughout.

## Ground rules
- Build exactly what is specified. Do not add "use X only if it scores better" logic.
- If something cannot be built as specified on the data, stop that step and report why. Do not improvise a substitute.
- Do not modify the cell generator or the eval harness. Both are frozen.
- Write outputs to a new folder in the project on the VM: report.md, figures/, and tables as parquet. Large intermediates go to the GCS bucket.

## Notation
- T = target line, k = knockdown, g = downstream gene, s = source line.
- x_skg = source LFC, se_skg = its SE. Both come from the effect tables, which hold LFC and SE as separate arrays.
- y_kg = target true LFC.
- S_kg = the sources that measure (k, g). S_k = the sources that measure knockdown k.
- Coverage is asymmetric across the 117 lines, so S_kg varies.
- p_kg = baseline prediction from the current precision-weighted transfer model.
- strength_k = RMS over genes of the baseline p_kg on the target. Define strong/mid/weak terciles once per target, from its baseline, and keep them fixed for every metric on that target.
- SE floor: wherever se² is used, floor it at the 5th percentile of that source's se².

## Step 0: Inventory
1. Locate the current transfer model and how it sets weights; any control-cell similarity measure; the effect tables; the fold definitions; the eval harness. Summarise the baseline weighting in the report.
2. Similarity: use the existing control-cell similarity measure. If the codebase has none, use cosine similarity between mean log-normalised control profiles over the top 2,000 HVGs shared across lines, and say so.
3. Dataset/study ID for each of the 117 lines: take it from table metadata and report the mapping. Datasets with fewer than 3 lines are grouped as "other" wherever dataset is used below.
4. Coverage for each fold line as target (sources = all lines except that fold): number of knockdowns that have truth and at least 1 source; number of (k, g) pairs in each |S_kg| bucket: {1, 2, 3–4, 5–9, 10–29, 30+}.
5. Per-source SE distribution, plus counts of zero or NaN SE.

## Step 1: Baseline
For each fold F, run the current model zero-shot on F. Report per fold: LFC MSE overall and by strength tercile; eval harness metrics.

## Step 2: Calibration diagnostics (no model changes; run on each fold line's baseline)
2a. Pooled slope of y on p: through the origin with OLS and with Huber (HuberRegressor, fit_intercept=False); also OLS with an intercept, reporting the intercept. Figure: hexbin log-density of y vs p, with the y=x line and the fitted line.
2b. Binned conditional mean: 50 quantile bins of p; mean y (95% CI) against mean p per bin. Separate through-origin slopes for |p| below and above the 90th percentile of |p|. State whether the curve is linear or S-shaped.
2c. Through-origin slope within each |S_kg| bucket, with counts.
2d. Leverage: the share of Σp² from the top 10% of knockdowns by strength.
2e. Per-knockdown slope β_k (OLS through origin over genes): distribution; scatter against strength_k and against agreement τ̂²_k (Step 4B); Spearman correlations.

## Step 3: Source error statistics
For every (fold F, training target U): 3 folds × 2 training targets = 6 runs; source pool = all lines except U and F. Only source–target pairs sharing ≥ 20,000 (k, g) pairs and ≥ 10 knockdowns. All means over shared (k, g).
3a. Per source: shared pairs; through-origin slope of y on x_s (scale diagnostic only); Pearson r; mean se²; V_s = mean(x_s²) − mean(y·x_s); t_s = V_s − mean(se_s²), floored at 1% of the median positive t_s (report how many floored).
3b. Per fold, OLS over source–target pairs pooled over the two training targets: log t_s,U ~ similarity(s, U) + log(median cells per perturbation in s) + dataset fixed effects; no target fixed effect. Report coefficients and R²; observed vs fitted coloured by target; residual mean by target; raw t_s against median cell count. This model supplies t̂_s,T from control-cell similarity.
3c. Error correlation net of the target, for source pairs sharing ≥ 20,000 pairs with each other and with U: ρ_ss' = [mean(x_s·x_s') − ½(mean(y·x_s) + mean(y·x_s'))] / sqrt(V_s·V_s'). Per fold, Fisher-z mean pooled over the two training targets, within datasets (ρ_w) and between (ρ_b). Figure: ρ against source–source control similarity. n_eff = S / (1 + (S−1)ρ_b) for S = 2, 5, 10, 30.
3d. Empirical averaging curve (descriptive; once per fold line as target, sources = all lines except that fold): knockdowns with ≥ 10 sources; for S ∈ {1, 2, 3, 5, 8, 10}, 20 random S-subsets, equal average, MSE against y; fit MSE = A + B/S; report plateau A and implied ρ; compare with 3c.

## Step 4: Revised predictor (per fold; fit on the two training targets, apply to F with sources = all lines except F; cumulative stages, metrics after each)
Stage A (weights), per knockdown over S_k: d_s = median over genes of floored se²_skg + t̂_s,T (3b); R_k with ρ_w for same-dataset pairs, ρ_b otherwise; Σ_k = D^½ R_k D^½; w ∝ Σ_k⁻¹·1; clip negatives to 0 and normalise (report fraction clipped); renormalise over present sources per gene (no re-solve); p^A_kg = Σ_s w_s x_skg.
Stage B (per-knockdown shrinkage from agreement; sources only, computed on the target predicted): τ̂²_k = pair-count-weighted mean over source pairs of mean_g(x_s·x_s'); σ̂²_k = mean over sources of mean_g(x_s²) − τ̂²_k; ρ̄_k = mean off-diagonal of R_k; n_eff,k = |S_k| / (1 + (|S_k|−1)ρ̄_k); B_k = clip(τ̂²_k / (τ̂²_k + σ̂²_k / n_eff,k), 0, 1); single-source knockdowns get the median B_k of their p^A strength decile, from the fold's two training targets; p^B = B_k · p^A. On F: Spearman and scatter of B_k against the oracle β_k (2e).
Stage C (global calibration): c = OLS through-origin slope of y on p^B pooled over the two training targets (each built with A–B, sources excluding it and F); p^C = c · p^B on F. Diagnostics only: c per training target, c within each strength tercile. Repeat 2b's figure for p^C on F.
Metrics on F for baseline, A, B, C: LFC MSE overall and by tercile; pooled slope of y on the prediction; eval harness metrics for the baseline and Stage C only.

## Step 5: Final fit for the unseen lines
Refit 3b, ρ_w, ρ_b, the decile medians and c with all three fold lines as training targets (each target: sources = all lines except that target). Save as a single frozen config. No predictions for the unseen lines in this task.

## Step 6: Report
Summary paragraph leading with the most important finding; direct answers (slope well below 1? rises with source count? linear or S-shaped? source dependence ρ_w, ρ_b, n_eff, plateau? transfer variance vs similarity, cell count, dataset? β_k variability and whether agreement predicts it); stage-wise metric tables per fold and the mean across folds; stability across folds (c, ρ_w, ρ_b, 3b coefficients), flagging parameters too unstable for the final config; recommendation for the 3 unseen lines; what I would have done differently and anything in the data that undermines the design.

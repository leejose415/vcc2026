# Evaluation: the scorer replica, the proxy and the pass criterion

## The official scorer, rebuilt (`analysis/mse_metric/mse_study.py`, `code/vcc6/scorer_study.py`)
- **Six metrics:** MSE (noise-corrected, capped, ratio of sums), PDS (cosine discrimination), LFC NMAE, direction fidelity, direction reach and DE Jaccard. Each is scaled between the eval's baseline (0) and replicate (1) anchors.
- **Reproduced against the harness:**
  - MSE rebuilt term by term matches to **9e-16** (HCT116) and exactly (Eval A).
  - PDS matches too.
  - The trace found three details beyond the public description: per-knockdown averaging over |G_p|, a per-row credit cap, and a total budget that scales every credit.

## The noise-free proxy (`code/vcc6/scorer_proxy.py`)
**What it computes:** all six components straight from a requested (LFC, call set), without generating cells.
- **Means:** the generator's own mean inversion.
- **MSE and PDS:** from the scorer's pseudobulk.
- **DE metrics:** from a synthetic DE table pushed through `cell-eval2`'s own DE path.

**Accuracy:** it predicts end-to-end MSE to 0.0006. Against real cells it overstates reach by 0.02–0.05. Step 22 traced that gap to the cells' gene ordering, not to signs.

**Speed:** 2–7 s per arm, against about 3 h to generate and score cells. Almost every idea was screened this way first.

## Evaluation lines
| Line | Role |
|---|---|
| **Eval A (H1, 292 knockdowns, 400 cells each)** | Main held-out line: the H1 tables are excluded from the sources. Scored on real cells by the official harness |
| HCT116final (278 knockdowns) | Second end-to-end line (screen setup: HCT116's own tables out) |
| K562, HCT116, HEK293T fold pools | Proxy guardrails: each line's tables held out |

Bundles are built by `scripts/harness/` from the real cells (`real_cache`, baseline and replicate anchors).

## How a change is accepted (written before each run)
- **Halves:** H1's knockdowns are split into two halves balanced by true DE count (seed 2026).
  - Every grid is **picked on one half and scored on the other**, in both directions.
  - Any model fitted inside a configuration is cross-fitted within the fit half (2 inner folds).
- **Criterion:**
  - a held-out gain of at least **+0.005 on both halves**;
  - every component above 0;
  - **no guardrail line more than 0.005 worse**.
- **Grid reporting:** each pick's position in its grid is reported (edge or interior), and a pick of the baseline row is labelled as such.
- **Before anything ships:** a cell run on Eval A and HCT116final, compared with the current entry by a **paired bootstrap** over knockdowns (2,000 resamples, `code/vcc6/boot_score.py`).
- **Volume:** when a rule moves calls between knockdowns, compare at matched *mean* call volume. The proxy rewards volume, and a median-matched comparison once produced a false gain.

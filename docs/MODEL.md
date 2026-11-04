# The model: from 117 screens to one prediction per knockdown
A knockdown's response in an unseen cell line is predicted by transferring what other lines showed for the same knockdown. The pieces were added one at a time, each with held-out evidence (see RESULTS.md).

## 1. Stage 1: precision-weighted transfer (`code/vcc6/stage1.py`, arm C + Gaussian shrink = "C+G")
1. **Sources:** all 117 tables, minus every table from the target's own cell line (matched on `obs.background`).
2. **One vote per cell line.**
   - A line's rows for the knockdown (several screens, conditions or drugs) are combined by a random-effects average with weights 1 / (SE² × multiplier).
   - Multipliers: single-cell ×1; knockout or rest-of-line reference ×4; replicate tables (Marson), drugs and Orion pseudobulk ×16.
   - A line's variance is never below its best single row's, so repeat screens of one line add no precision.
3. **Pool across lines:** a precision-weighted average of the line votes (DerSimonian–Laird random effects) gives μ, its sampling variance var_fe, and the between-line variance τ² (floored at a global 10⁻⁴).
4. **Shrink:** LFC = μ × v0 / (v0 + var_fe + τ²), a Gaussian prior.
5. **Calls (K0):** BH on |LFC| / the target line's 400-cell measurement SE (from the line's own controls).

**Bug found and fixed on 27 Oct:** an earlier slab-and-spike prior multiplied effects by about 0.003. Removing it took the RPE1 hold-out from −0.20 to +0.12.

## 2. Config 1′: magnitudes in the scorer's space (`code/vcc6/step15c.py`, frozen in `CONFIG1_FROZEN.md` as config 1, then 1′)
- **Called genes:** 0.5 × C+G's LFC.
- **Uncalled genes:** the line's control profile + 0.6 × c × Δ_pool.
  - Δ_pool is the sources' knockdown-minus-control pseudobulk (log1p CP50K), pooled with the same stage-1 weights.
  - c = z²/(z² + 1), where z is the agreement z across source lines; c = 1 with a single source and 0 with none.
- **Why:** MSE and PDS read the pseudobulk, not the LFC. About 80% of C+G's error sat on uncalled genes, and its largest LFCs sat on barely expressed genes.

## 3. Calls: hybrid T5 (`code/vcc6/step19.py`)
- **K4:** the top 1,200 tested genes per knockdown by C+G's z (|LFC| / SE). The scorer's DE metrics reward calling more genes, up to a point.
- **(b), predicted knockdown strength:** the mean `n_real` of the target across source tables (the target line's own tables excluded). It predicts the true DE count (held-out ρ 0.64 on H1).
- **Hybrid T5:** K4 where b ≥ 5, else K0. Weak knockdowns keep a short call list, strong ones get 1,200.

## 4. The no-source fix
A gated gene that no source measured is passed to the generator as measured with no change (n_sources ≥ 1). It then lands on the pooled control profile the scorer uses, not on the mean per-cell CPM (+0.0014).

## 5. The fixed strength shift (`code/vcc6/step20.py`, part C2)
- **S:** a shared "strong-minus-weak" response vector, from 13 source lines with enough strong and weak knockdowns.
  - It uses genes measured by at least half of those lines.
  - It mostly carries p53 up and proliferation down.
- **How it's applied:** on the target line's gated genes, as (−0.25, 0, +0.25) × S by (b) tercile (cuts 3.0 / 67.7 from H1).

## 6. The F5 program term (`code/vcc6/step20.py`, part F; `step21.py`)
- **Basis E4w:** an SVD of all genetic source rows in Δ space.
  - Each table's > 5 CPM gate; per-table centring.
  - Row weight 1/√max(n_b, 50), so tiny backgrounds can't steer components.
  - Tahoe counted as one background; Arc H1 out.
  - K = 100 components.
- **Regression M3:** for each knockdown, predict the target line's component scores from each source line's scores (ridge, one feature block per source line) plus the target gene's DepMap features. β is fitted on H1.
- **Term:** ŝ × basis on the gated genes, added to Δ for knockdowns with b ≥ 5 only (γ = 1).

## 7. Row 4: the final entry
Row 4 = sections 1 to 6 together:
- hybrid T5 calls + config 1′ magnitudes;
- the no-source fix, the fixed shift and the F5 term;
- the frozen generator.

| | Eval A (held-out H1 cells) |
|---|---|
| Row 4 | 0.2895 |
| Live entry (config 1′ + K0) | +0.047 [+0.032, +0.061] |

The D/E/F inputs are built by `code/vcc6/row4_def.py inputs`, which also checks that it reproduces the evaluated Eval A input bit for bit.

## Earlier models (validation phase, see RESULTS.md)
- **Model 2** (`code/vcc_2026_3c/vcc2_meta.py`): random-effects meta-analysis on 8 cell lines with top-N calls and a detection floor. Best leaderboard entry +0.091.
- **Transformer** (`code/transformer/`): gene-gene attention over pretrained gene embeddings (ESM2, STRING, DepMap, text). Better on held-out genes, worse on the leaderboard (PDS collapse).

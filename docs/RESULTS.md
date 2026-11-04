# Results
Full detail, with every table, is in `results/stage1/2026-10-28_scorer_space/report.md` (15 dated studies). This page is the summary.

## Validation phase: public leaderboard (lines A–C; MSE was 0 for every entry)

| Entry | Model | Overall | PDS | NMAE | Fidelity | Reach | Jaccard |
|---|---|---|---|---|---|---|---|
| v1 | first pipeline | −0.270 | 0.088 | 0.023 | −1.688 | 0.029 | −0.074 |
| v2 | Model 2 (random-effects meta-analysis), gated | −0.054 | 0.350 | 0.024 | −0.805 | 0.083 | 0.026 |
| **v3** | **Model 2 + top-2,000 calls + detection floor** | **+0.091** | 0.416 | 0.048 | −0.088 | 0.192 | −0.024 |
| v3a | same, top-1,000 | +0.087 | 0.420 | 0.076 | −0.121 | 0.180 | −0.032 |
| v4 | every gene at 0.75 × Δ, no calls | +0.070 | 0.403 | 0.119 | −0.241 | 0.149 | −0.011 |
| v5 | transformer, every gene + floor | +0.064 | 0.214 | 0.058 | +0.012 | 0.110 | −0.010 |
| v6 | transformer + PIE evidence | +0.057 | 0.197 | 0.046 | +0.009 | 0.098 | −0.009 |
| v7 | transformer, Δ̄ + 2(Δ − Δ̄) + floor | +0.078 | 0.292 | 0.062 | +0.001 | 0.123 | −0.011 |

## Final phase: held-out H1 cells (Eval A), official harness, paired bootstrap

| Model | Eval A | Change | HCT116final |
|---|---|---|---|
| Model 2 port | 0.184 | — | — |
| C+G | 0.219 | +0.034 [+0.020, +0.048] | 0.190 |
| Config 1′ + K0 | 0.243 | +0.024 [+0.016, +0.033] | 0.362 |
| Config 1′ + K4 1,200 | 0.264 | +0.045 vs C+G | 0.348 |
| Hybrid T5 | 0.263 | +0.020 vs K0 | 0.355 |
| **Row 4 (final entry)** | **0.2895** | **+0.047 [+0.032, +0.061] vs config 1′ + K0** | 0.362 |

**Row 4's Eval A components:** MSE 0.187, PDS 0.766, NMAE 0.094, fidelity 0.469, reach 0.079, Jaccard 0.142.

**Final-phase entries:**
- 28 Oct: C+G.
- 30 Oct: config 1′ + K0.
- **3 Nov: row 4 (the counting entry).**

Test scores are hidden until the challenge announces winners.

## Tested and rejected
Each idea was held to the criterion in EVALUATION.md.
- **Transformer:** leaderboard 0.064–0.078 against v3's 0.091 (PDS collapse).
- **Allocating calls by predicted strength:** about 0 at matched mean volume (an artefact we caught).
- **Expression tiers, support gates, same-lab upweighting:** no gain.
- **Uncertainty-aware magnitudes:** H1 +0.007 / +0.013, but guardrails −0.006 to −0.029.
- **Cell-line programs from control cells** (pathway-anchored, joint ridge): every row lost, or ~0.
- **Line-agreement weights:** +0.001.
- **DepMap essentiality in the shift:** no gain.
- **Generator v8 (ordering p-values for reach):** failed its call-set check.

## What we learned about the biology
- **Hidden test lines from control cells alone:** D myeloid, E mesenchymal, F fibroblast-like.
- **H1's knockdown response is low-rank but line-specific:** its own top-20 components hold 68% of strong-knockdown energy; components shared across sources hold 7–14%. The leading H1 component is exit from pluripotency.
- **Ribosomal-protein knockdowns trigger p53** in TP53-wild-type lines (0.37–0.45), but not in TP53-mutant K562 or Jurkat (0.05).
- **Lines agree more closely than their SEs imply** (Cochran's Q/df 0.16–0.48). Lab protocol matters little (+0.10).
- **The proxy-minus-cells reach gap (0.032) is entirely gene ordering.** A truth-sign-first ordering would add about +0.08 reach: a measured target for future work.

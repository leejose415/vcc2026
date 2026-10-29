# Config 1 — frozen 2026-10-29

`CONFIG s_c 0.5, m 0.5, kappa 1.0, c1 1.0`

- v7 data loader + model cells: `sha256 cells 4cff99dde2f37c974a8ef72d38756db5651b500e28113a54da61cfed5cb29c38`
- `code/vcc6/stage1.py`: `sha256 stage1 1021148db09fd48137072d373c9a4593a472e6f11e12e0409f6bb64041ad9af2`
- `code/vcc6/run_line.py`: `sha256 run_line 39d8143a1b010fccd59ae95d8366718390af41c24c62577591fbc0ba2bc78de8`
- the 117 effect tables in `lib/effects_v5` and `lib/effects_v5_extra`, sorted stems joined by newlines:
  `sha256 tables 6c3fcff7e1e93457ecbb03849091f67d537cd0e5e84bf04f33f53baba20affba`

These are checked by `tests/test_config1.py`. The cells are hashed by source, so running the notebook or toggling the
run flags in its setup cell (`RUN_EVAL`, `RUN_UPLOAD`, `SUBMIT`, `MAKE_SWAP`) does not break the freeze. Changing
`CONFIG`, the data-loader or model cell, `stage1.py`, `run_line.py` or the table set does. A new model goes in a new
module or notebook (`vcc_2026_8.ipynb`), never in place.

## What config 1 is
C+G (stage-1 arm C + Gaussian shrink, all 117 tables, the line's own tables excluded), with new magnitudes handed to
the frozen generator:
- **Called genes** (C+G's calls, K0): 0.5 × C+G's LFC.
- **Non-called genes:** control profile + 0.5 × c × Δ_pool. Δ_pool is the sources' knockdown-minus-control
  pseudobulk, pooled by the stage-1 machinery.
- **c:** z²/(z² + 1) for genes with ≥ 2 source lines, where z = |mean| / (SD / √n). c = 1 for one source line, and 0
  for none.

## Evidence (Eval A, H1 held out, scaled)
Config 1 has not been run end to end. Its score is the noise-free proxy, which matched C+G's end-to-end overall to
0.0013.

| Metric | C+G, end-to-end run | C+G, proxy | Config 1, proxy (held out) |
|---|---|---|---|
| MSE | 0 (raw 1.496) | 0 | 0.113 |
| PDS | 0.723 | 0.724 | 0.808 |
| NMAE | 0.040 | 0.040 | 0.001 |
| Direction fidelity | 0.380 | 0.371 | 0.371 |
| Direction reach | 0.070 | 0.074 | 0.073 |
| Jaccard | 0.097 | 0.095 | 0.095 |
| **Overall** | **0.219** | 0.217 | **0.244** |

Sources: `results/stage1/2026-10-28_scorer_space/tables/step11_heldout.csv` (arm i_lofo, H1) and `step9.json`
(EvalA_CG); `results/v7/config1_record.json`.

## Freeze check (2026-10-29, `results/v7/freeze_check/`)
The proxy's 0.244 came from the scorer-space study's code and caches (`cache/fold_H1.npz`, `cache/s11_H1.npz`, copied
from vcc2027, md5 checked). The v7 notebook rebuilds config 1 with its own code. To confirm the two agree, v7's
unchanged setup, data-loader and model cells were run for Eval A (`run_v7_model.py`) and compared
(`compare.py` → `compare.json`):

- Calls: identical (0 of 292 × 18,533 differ).
- Generator means on gated genes: max relative difference 1.1e-7.
- Called-gene LFC: exactly 0.5 × C+G's.
- Proxy overall from v7's input: 0.24361700, against 0.24361700 from the caches. All six components agree to 1e-10.

## Open risk
NMAE +0.001 is fragile: +0.008 and −0.006 on the two H1 knockdown halves (scorer-space Step 16). The real end-to-end
Eval A score comes from v7's section 4.

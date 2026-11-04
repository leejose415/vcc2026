# Virtual Cell Challenge 2026: knockdown-response prediction

Our entry to Arc Institute's [Virtual Cell Challenge 2026](https://virtualcellchallenge.org).
- **The task:** predict, for a cell line never seen in training, how the expression of 18,533 genes responds to CRISPRi knockdown of each of 300 genes. Then deliver 360,000 synthetic single cells (3 lines × 300 knockdowns × 400 cells) whose statistics carry that prediction.
- **Scoring:** six metrics computed by `cell-eval2` (MSE, PDS, LFC NMAE, direction fidelity, direction reach, DE Jaccard), each scaled between a baseline and a replicate.

## What is in this repository

| Part | Where | Docs |
|---|---|---|
| Data sources and downloads | `scripts/download_*`, `scripts/build_raw_cells.py`, `scripts/convert_new_contexts.py`, `scripts/orion_stream.py` | [docs/DATA_SOURCES.md](docs/DATA_SOURCES.md) |
| Cell QC and perturbation QC | `scripts/qc_replogle_umi.py` → `qc_replogle_apply.py` → `write_replogle_qc.py` → `qc_perturbations.py` (`qc_table_perturbations.py` for table-only sources) | [docs/CELL_QC.md](docs/CELL_QC.md) |
| Effect tables (117 sources, scorer regime) | `scripts/build_effects_v5*.py`, `scripts/build_effects_extra.py` | [docs/EFFECT_TABLES.md](docs/EFFECT_TABLES.md) |
| Transfer model (stage 1, config 1′, the final "row 4") | `code/vcc6/stage1.py`, `run_line.py`, `step15c.py`, `step19.py`–`step21.py`, `row4_def.py` | [docs/MODEL.md](docs/MODEL.md) |
| Cell generator (frozen) | `code/cellgen/cell_generator.py`, `code/vcc6/gen_cells.py`; experimental v8 in `code/vcc6/gen_v8/` | [docs/GENERATOR.md](docs/GENERATOR.md) |
| Scorer replica, noise-free proxy, evaluation | `code/vcc6/scorer_*.py`, `boot_score.py`, `analysis/mse_metric/`, `scripts/harness/` | [docs/EVALUATION.md](docs/EVALUATION.md) |
| Every study, with results | `code/vcc6/step15*.py`–`step22.py`, `results/stage1/*/tables/` | [docs/RESULTS.md](docs/RESULTS.md); full log in [results/stage1/2026-10-28_scorer_space/report.md](results/stage1/2026-10-28_scorer_space/report.md) |
| Earlier models (validation phase) | `code/vcc_2026_3c/` (meta-analysis + top-N calls), `code/transformer/` (attention model), `analysis/v6lin/` | [docs/RESULTS.md](docs/RESULTS.md) |
| Tests | `tests/` (frozen generator hash, config-1 freeze, builders, transformer) | |

## Results in one table
Held-out H1 cells (Eval A, 292 knockdowns, the official harness, paired bootstrap over knockdowns). Test-line scores are hidden until the challenge announces winners.

| Model | Eval A overall |
|---|---|
| Model 2 port (random-effects meta-analysis) | 0.184 |
| C+G (stage-1 transfer + Gaussian shrink) | 0.219 |
| Config 1′ + K0 | 0.243 |
| Config 1′ + K4 (1,200 calls) | 0.264 |
| **Row 4: final entry** (hybrid calls + strength shift + DepMap program term) | **0.290** (+0.047 [+0.032, +0.061] over config 1′) |

## Reproducing
- **Environment:** Python 3.12; `pip install -r requirements.lock` (cell-eval2 0.18.0, the `vcc` CLI 0.3.0, anndata, scanpy, polars, torch CPU).
- **Paths:** the code expects the repository at `~/vcc`, so clone or symlink it there. Data goes under `~/vcc/lib/` and `~/vcc/data*/`, none of which is in git.
  - Three files are frozen and pinned by hash tests: `code/cellgen/cell_generator.py`, `code/vcc6/stage1.py` and `code/vcc6/run_line.py`. That is why the root path is not configurable.
- **Order:**
  1. Download the sources (DATA_SOURCES).
  2. Cell QC and perturbation QC (CELL_QC).
  3. Effect tables (EFFECT_TABLES).
  4. Evaluation bundles (`scripts/harness/`).
  5. Model and generator (`run_line.py predict`).
  6. Score (`run_line.py score`, `boot_score.py`).
- **The final entry** was built by `code/vcc6/row4_def.py` (inputs, merge, checks) with `row4_gen.sh` / `row4_final.sh`. Those two shell scripts assume our two-VM setup and are kept as a record.
- **Memory:** scoring a 135k-cell eval with `cell-eval2` needs about 33 GB of RAM, and `vcc prep` on a 360k-cell file about 34 GB plus about 19 GB of temporary disk. Use swap on a 32 GB machine.

## Licence
Code: MIT (see `LICENSE`). Data belong to their original authors (see DATA_SOURCES); none is redistributed here.

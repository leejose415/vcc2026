# cell_generator.py — frozen 2026-10-26

`sha256 422f8a3d8a04f2f07c8039212d79155f1bb114e217d1b59357d3aad42333a196`

This hash is checked by `tests/test_cellgen.py::test_generator_is_frozen`. The file is read-only. A new or changed
generator goes in a new module with its own oracle test. It never replaces this one in place.

## What the generator does

It turns one knockdown's predicted v5 LFC and significance call into integer-count cells. Each gene's mean CPM is
exact. Its rank test against the v5 control quantiles lands on the requested side at BH FDR 0.05. The cell shape is
otherwise free, because no metric reads it. The module docstring has the details.

## Oracle evidence

The test is `scripts/cellgen_oracle_test.py`. Results are in `results/cellgen_tests/`, on vcc2026 and in full on
vcc2027.

**Setup:**
- 10 H1 knockdowns from Eval A, seed 2026: ANKZF1, C1QBP, CXCL12, INSIG1, IRF3, JAZF1, MED13, PIAS1, PMS1, TAF11.
- Each one's cells are split 50/50.
- The truth comes from half 1, using the v5 definitions against the v5 H1 controls.
- Every arm is scored against half 2 plus Eval A's 18,400 scoring controls.
- Scoring is official cell-eval2 (vcc2026, pdex), with a bundle built on that reference.

| arm | scaled overall | MSE | sig. Jaccard |
|---|---|---|---|
| replicate (real half 1) | 1.0971 | 1.000 | 1.114 |
| **oracle_nhalf_pb** (n = half 1's cells; mean target = half 1's pseudobulk) | **1.0957** | 1.000 | 1.107 |
| oracle_nhalf (n = half 1's cells; mean CPM target) | 1.0734 | 0.885 | 1.108 |
| oracle, 400 cells, mean CPM target | 0.9642 | 0.401 | 0.935 |
| null (lfc = 0, no calls) | −0.1389 | 0 | −0.029 |
| control-paste | −0.1272 | 0 | −0.029 |

## Decision

**The confirming oracle is oracle_nhalf_pb.** Its gap to the replicate is 0.0014, inside the 0.05 criterion. The
null arm sits at the control-paste floor.

**The 400-cell / mean-CPM arm is not the success criterion** when the call set and the mean estimator come from a
different cell count than the generated cells. Here both came from 200-cell halves, which causes two artefacts:
- **Sig. Jaccard.** A 200-cell call set realized in 400 cells adds calls. These are small decreases, mostly in
  low-expression genes, that become significant at twice the cells whatever the cell shape.
- **MSE.** The mean per-cell CPM of 200 cells is a noisier estimate than their pseudobulk, which is what the MSE
  metric aggregates. The generator reproduced its input exactly, so that noise is the whole excess.

In production the v5 calls are made at 400 cells, the same as the generated cells.

**Attribution.** Any later score gap on Eval A or B belongs to the Stage-1 model, not to this generator.

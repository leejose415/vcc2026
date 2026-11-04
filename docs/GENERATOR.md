# The cell generator (frozen)
`code/cellgen/cell_generator.py`, frozen on 2026-10-26. Its sha256 (422f8a3d…) is checked by `tests/test_cellgen.py`.

**What it does:** turns one knockdown's predicted LFC vector and call set into 400 integer-count cells.
- **Means are exact:** each gene's mean CPM is (control mean CPM + ε) · 2^LFC − ε. Each gene's total count is preserved exactly by largest-remainder rounding.
- **Calls are realised:**
  - Every called gene's rank test against the stored control quantiles rejects at BH 5%.
  - Every other tested gene does not reject.
- **How called genes are shaped:** the control shape is pulled toward the mean, with spread c from (1, 0.8, … 0). The first c that rejects is kept.
- **Uncalled genes** keep the control shape. If it still rejects, the change is carried by the top cells, or clipped off the top.
- **Deterministic:** gene g draws from `default_rng([seed, g])`, and the seed comes from the knockdown name.

## Evidence that it loses nothing on its own
- **Oracle:** feed it the true effects of 10 H1 knockdowns, scored against a held-out half of the real cells (`scripts/cellgen_oracle_test.py`). It reaches **1.0957 against the real replicate's 1.0971 (gap 0.0014)**.
- **Exactness:** generated pseudobulk against requested has slope 0.998. In the final build, r ≥ 0.999999 on every knockdown.

So any score gap belongs to the model, not the generator.

## Generator v8 (experiment, not adopted)
- **Code:** `code/vcc6/gen_v8/`, a copy of the frozen file plus an optional rank input.
- **Idea:** start each called gene's spread search at c = 1 − rank, so the most confident genes get the smallest p-values. Direction reach ranks genes by p.
- **Result:** means stayed bit-identical, but 1.9% of the realised call sets moved (Jaccard 0.981 < 0.99). It failed its own check, so it was stopped.

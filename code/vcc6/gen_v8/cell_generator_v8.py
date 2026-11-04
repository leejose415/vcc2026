"""Cell generator v8 (Step 22 Part R3, 2026-11-03): a COPY of the frozen code/cellgen/cell_generator.py (not edited) with
one optional input, rank (per gene, in [0, 1] for called genes, 1 = most confident; NaN / None elsewhere). A called
gene's spread search starts at c_top = 1 - rank * (1 - C_MIN) instead of 1 and continues down the frozen SPREADS
below c_top, then the frozen concentrate fallback. rank 0 (each knockdown's least confident called gene) and
rank=None reproduce the frozen generator bit for bit (same candidate list, same draws). C_MIN = 0.0 is the frozen
search's own tightest spread: every cell at the target mean before integerisation, so the mean is unchanged.
Everything below is the frozen text.

Cell generator: one knockdown's predicted LFC and significance call -> the cells the scorer grades.

The vcc2026 metrics read two things from a knockdown's cells: each gene's mean CPM (pds, mse, lfc_nmae and the sign
the direction metrics use) and which genes the Wilcoxon test calls (sig_jaccard, direction fidelity / reach). The
generator realizes exactly those two:
  mean   mean CPM of the generated cells = (control_mean_cpm + 1e-9) * 2**lfc - 1e-9, the v5 LFC inverted (so a gene
         absent from the controls keeps its level). Each gene's total count is exact, so the mean is within 1 / n_cells
         of a count.
  ranks  every gene in de_sig: the scorer's test (two-sided Mann-Whitney U on log1p(CPM), normal approximation, tie +
         continuity corrected) against the stored control quantiles rejects at BH FDR 0.05 over the gated genes; every
         other gene: it does not reject.
Shape is otherwise free (no metric reads it), so each gene's cells start from a control-shaped base:
  base     n_cells stratified draws from the gene's control distribution (the quantile sidecar), at library_size;
  not DE   the base scaled to the target mean (the control shape). If that still rejects, the mean change is carried by
           the top fraction w of the cells only (largest w that does not reject), which moves the rank sum least;
  DE       the base scaled to the target mean, then pulled toward it, x = mean + c * (x - mean): the largest c in
           (0, 1] that rejects (c = 1 first, the plain scaled control shape). This is the per-gene search the test needs.
  Genes outside the sidecar (control mean <= 5 CPM: the scorer never tests them) carry the mean only: their total
  count is spread over the cells at random (multinomial).
Counts: counts = round(cpm * library_size / 1e6), library_size = 20,000 UMI (the scorer regime). The rounding keeps each
  gene's total exact (largest remainder: every cell gets the floor or the ceiling); plain rounding would erase every
  gene below 25 CPM, half a count at 20k.
Check: the finished cells are re-tested as the scorer computes CPM (each cell's own row sum) against the quantiles,
  with BH over the gated genes; genes on the wrong side are rebuilt once with tighter margins. The report is in
  uns['cellgen'] (requested vs realized calls, the composition sum).
Composition: the scorer renormalizes each cell, so the realized means are the requested ones times 1e6 / sum(means).
  An exact v5 LFC vector sums to 1e6 by construction; a model's may not (uns['cellgen']['composition']).
Deterministic given the seed: gene g draws from default_rng([seed, g]).

FROZEN 2026-10-26 (code/cellgen/FROZEN.md; the file hash is checked by tests/test_cellgen.py). Confirming oracle:
oracle_nhalf_pb (scripts/cellgen_oracle_test.py --extra) scored 1.0957 against the real replicate's 1.0971 on 10 H1
knockdowns (gap 0.0014). The 400-cell / mean-CPM arm (0.9642) is not the success criterion: its call set and mean
estimator came from 200-cell halves, a different cell count from the generated cells. Any later score gap on Eval A/B
belongs to the Stage-1 model, not to this generator. Do not change this file; a new generator is a new module.
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import false_discovery_control
from scipy.special import ndtr

EPS, FDR = 1e-9, 0.05
SIG_MARGIN, NONSIG_P = 0.5, 0.1           # DE genes aim at p <= 0.5 x the BH cut-off; other genes at p > 0.1
SPREADS = (1.0, 0.8, 0.6, 0.45, 0.3, 0.2, 0.12, 0.07, 0.03, 0.0)
TOPS = (1.0, 0.5, 0.25, 0.12, 0.06, 0.03, 0.015, 0.0)        # 0.0 -> a single cell
C_MIN = 0.0                               # v8: the tightest spread (SPREADS' last value)


class ControlQuantiles:
    """The v5 sidecar (<source>_control_quantiles.npz): per gated panel gene, the sorted unique nonzero log1p(CPM) of
    the sampled controls, their multiplicities and the zeros. Loaded eagerly; per-gene pieces cached for the test."""
    def __init__(self, path):
        q = dict(np.load(path))
        self.path = str(path); self.col = q['panel_col'].astype(np.int64); self.ptr = q['ptr'].astype(np.int64)
        self.values = q['values']; self.counts = q['counts'].astype(np.int64); self.zeros = q['zeros'].astype(np.int64)
        self.n0 = int(q['n_controls'])
        self.index = {int(c): i for i, c in enumerate(self.col)}

    def gene(self, i):
        a, b = self.ptr[i], self.ptr[i + 1]
        v, c = self.values[a:b], self.counts[a:b]
        cf = c.astype(float)
        return dict(v=v, c=c, below=np.concatenate([[0], np.cumsum(c)]), tie_c=float((cf ** 3 - cf).sum()),
                    nz0=int(self.zeros[i]), cum=(self.zeros[i] + np.cumsum(c)) / self.n0)


def rank_p(x, g, n0):
    """Two-sided MWU p of the values x (log1p CPM) against one gene's stored control quantiles: the scorer's test
    (normal approximation, tie + continuity corrected; build_effects_v5.quantile_p with the control sums cached)."""
    n1 = len(x); xs, p = np.unique(x[x > 0], return_counts=True); z1 = n1 - int(p.sum())
    v, cnts, below, nz0 = g['v'], g['c'], g['below'], g['nz0']
    i = np.searchsorted(v, xs); j = np.minimum(i, len(v) - 1)
    ce = np.where((i < len(v)) & (v[j] == xs), cnts[j], 0).astype(float)
    U1 = 0.5 * nz0 * z1 + float((p * (nz0 + below[i] + 0.5 * ce)).sum())
    t = ce + p
    tie = g['tie_c'] + float(((t ** 3 - t) - (ce ** 3 - ce)).sum()) + (nz0 + z1) ** 3 - (nz0 + z1)
    n = n0 + n1; s_sq = n1 * n0 / 12.0 * ((n + 1.0) - tie / (n * (n - 1.0)))
    if s_sq <= 0: return 1.0
    z = (max(U1, n1 * n0 - U1) - n1 * n0 / 2.0 - 0.5) / np.sqrt(s_sq)
    return float(min(max(2.0 * ndtr(-z), 0.0), 1.0))


def integerize(y, rng):
    """Floor or ceiling per cell, the gene's total = round(sum y) (largest remainders; ties broken at random)."""
    f = np.floor(y); k = int(round(float(y.sum()) - float(f.sum())))
    if k > 0:
        order = np.lexsort((rng.random(len(y)), -(y - f)))
        f[order[:k]] += 1
    return f.astype(np.int64)


def base_draw(g, n, scale, rng):
    """n stratified draws from the control distribution (zeros, then the sorted values), in count space, shuffled."""
    u = (np.arange(n) + rng.random(n)) / n
    zero_share = g['nz0'] / (g['nz0'] + g['below'][-1])
    k = np.minimum(np.searchsorted(g['cum'], u, side='right'), len(g['v']) - 1)    # first value whose cumulative share > u
    cpm = np.where(u < zero_share, 0.0, np.expm1(g['v'][k]))
    return rng.permutation(cpm * scale / 1e6)


def _spread(b, mu, c):
    mb = b.mean()
    y = b * (mu / mb) if mb > 0 else np.full_like(b, mu)
    return mu + c * (y - mu)


def _tops(b, mu, w):
    """The mean change carried by the top ceil(w n) cells only (None when they cannot absorb a decrease)."""
    n = len(b); m = max(1, int(np.ceil(w * n))); top = np.argsort(-b, kind='stable')[:m]
    y = b.astype(float).copy(); d = n * (mu - b.mean())
    if d >= 0:
        y[top] += d / m
    else:
        s = y[top].sum()
        if s < -d: return None
        y[top] *= 1 + d / s
    return y


def _concentrate(n, mu, w, rng):
    """The gene's whole mass on ceil(w n) random cells, equal shares; every other cell at zero."""
    m = max(1, int(np.ceil(w * n))); y = np.zeros(n)
    y[rng.permutation(n)[:m]] = n * mu / m
    return y


def _clip_top(b, mu):
    """A decrease taken off the highest values only: every cell above a level T is cut to T (water-filling), which
    moves the rank sum least, the control distribution being flat in its upper tail. None for an increase."""
    d = len(b) * (b.mean() - mu)
    if d <= 0: return None
    lo, hi = 0.0, float(b.max())
    for _ in range(60):
        t = (lo + hi) / 2
        lo, hi = (t, hi) if np.maximum(b - t, 0).sum() > d else (lo, t)
    return np.minimum(b, hi)


def spreads_from(c_top):
    """v8: the spread candidates starting at c_top (the frozen SPREADS when c_top >= 1)."""
    return SPREADS if c_top >= 1.0 else (float(c_top),) + tuple(c for c in SPREADS if c < c_top)


def realize_gene(g, n0, mu, sig, p_sig, p_nonsig, scale, rng, n, c_top=1.0):
    """One tested gene: integer counts whose mean is mu (counts) and whose test is on the requested side, judged with
    every cell at library `scale` (the finished cells are re-tested on their own row sums by check()). Candidate shapes
    in order of preference (the first one on the requested side is kept):
      DE      the control shape pulled toward the mean (c = 1 ... 0), then the mass carried by fewer and fewer cells
              (the rest at zero: the way a low-expressed gene's decrease shows in ranks);
      not DE  the control shape; then, for an increase, the extra mass on the top cells only, for a decrease, the
              highest values clipped (water-filling)."""
    b = base_draw(g, n, scale, rng)
    to_x = lambda k: np.log1p(k * 1e6 / scale)
    ok = (lambda p: p <= p_sig) if sig else (lambda p: p > p_nonsig)
    best, best_key = None, None
    if sig:
        tried = [_spread(b, mu, c) for c in spreads_from(c_top)] + [_concentrate(n, mu, w, rng) for w in TOPS[1:]]
    else:
        tried = [_spread(b, mu, 1.0)] + ([_tops(b, mu, w) for w in TOPS] if mu >= b.mean() else [_clip_top(b, mu)])
    for y in tried:
        if y is None: continue
        k = integerize(np.maximum(y, 0), rng); p = rank_p(to_x(k), g, n0)
        if ok(p): return k, p, True
        key = -p if sig else p
        if best_key is None or key > best_key: best, best_key = (k, p), key
    return best[0], best[1], False


def target_means(lfc, control_mean_cpm):
    """The v5 LFC inverted: mean CPM = (m0 + eps) 2^lfc - eps; NaN lfc -> the control mean; NaN m0 -> 0."""
    m0 = np.nan_to_num(np.asarray(control_mean_cpm, float), nan=0.0)
    l = np.nan_to_num(np.asarray(lfc, float), nan=0.0, posinf=60.0, neginf=-60.0)
    return np.maximum((m0 + EPS) * np.exp2(np.clip(l, -60, 60)) - EPS, 0.0)


def generate_cells(lfc, de_sig, control_quantiles, control_mean_cpm, n_cells=400, library_size=20000, seed=0,
                   target='perturbed', genes=None, alpha=FDR, rank=None):
    """lfc: (n_genes,) log2 fold change vs control mean CPM (NaN = unmeasured: control-like); de_sig: (n_genes,) bool;
    control_quantiles: ControlQuantiles or the sidecar path; control_mean_cpm: (n_genes,) (the baseline CSV's
    control_mean_cpm). Returns an AnnData of n_cells x n_genes integer counts, obs['target'] = obs['target_gene'] =
    target, uns['cellgen'] = the realization report. Genes are panel columns in order (genes = their names)."""
    cq = control_quantiles if isinstance(control_quantiles, ControlQuantiles) else ControlQuantiles(control_quantiles)
    mu_cpm = target_means(lfc, control_mean_cpm); G = len(mu_cpm); n = int(n_cells); L = float(library_size)
    sig = np.zeros(G, bool) if de_sig is None else np.asarray(de_sig, bool).copy()
    if len(sig) != G: raise ValueError('lfc, de_sig and control_mean_cpm must have one value per gene')
    mu = mu_cpm * L / 1e6                                                  # target mean count per cell
    rk = None if rank is None else np.asarray(rank, float)
    ctop = lambda col: 1.0 if (rk is None or not sig[col] or not np.isfinite(rk[col])) else 1.0 - float(np.clip(rk[col], 0, 1)) * (1.0 - C_MIN)
    tested = np.array([c for c in cq.col if c < G]); K = int(sig[tested].sum()); M = len(tested)
    p_sig = SIG_MARGIN * alpha * max(K, 1) / max(M, 1)
    X = np.zeros((n, G), np.int64)
    gated = np.zeros(G, bool); gated[tested] = True
    for col in np.flatnonzero(~gated):                                     # untested: the mean only
        if mu[col] <= 0: continue
        rng = np.random.default_rng([seed, int(col)])
        X[:, col] = rng.multinomial(int(round(mu[col] * n)), np.full(n, 1.0 / n))
    unmet = 0
    for col in tested:
        rng = np.random.default_rng([seed, int(col)])
        k, _, met = realize_gene(cq.gene(cq.index[int(col)]), cq.n0, mu[col], bool(sig[col]), p_sig, NONSIG_P, L, rng, n, ctop(int(col)))
        X[:, col] = k; unmet += not met
    report = check(X, tested, sig, cq, alpha)
    redo = [int(c) for c in tested[report.pop('_wrong')]]
    for col in redo:                                                       # one rebuild with tighter margins
        rng = np.random.default_rng([seed, int(col), 1])
        X[:, col] = realize_gene(cq.gene(cq.index[col]), cq.n0, mu[col], bool(sig[col]), 0.1 * p_sig, 0.5, L, rng, n, ctop(col))[0]
    if redo:
        report = check(X, tested, sig, cq, alpha); report.pop('_wrong')
    lib = X.sum(1)
    report.update(first_pass_unmet=int(unmet), rebuilt_genes=len(redo), composition=float(mu_cpm.sum() / 1e6), library_median=float(np.median(lib)),
                  library_size=L, n_cells=n, seed=int(seed), quantiles=cq.path,
                  mean_count_max_abs_error=float(np.abs(X.mean(0) - mu).max()))
    obs = pd.DataFrame(dict(target=target, target_gene=target), index=[f'{target}_{i:04d}' for i in range(n)])
    var = pd.DataFrame(index=pd.Index(genes if genes is not None else [f'g{j}' for j in range(G)], name='gene'))
    a = ad.AnnData(X=sp.csr_matrix(X.astype(np.int32)), obs=obs, var=var)
    a.uns['cellgen'] = report
    return a


def check(X, tested, sig, cq, alpha=FDR):
    """The finished cells' calls, as the scorer computes them (CPM over each cell's own row sum) against the stored
    controls, BH over the gated genes: requested vs realized."""
    lib = X.sum(1).astype(float); f = 1e6 / np.where(lib > 0, lib, 1)
    p = np.array([rank_p(np.log1p(X[:, c] * f), cq.gene(cq.index[int(c)]), cq.n0) for c in tested])
    called = false_discovery_control(p, method='bh') < alpha if len(p) else np.zeros(0, bool)
    want = sig[tested]
    tp = int((called & want).sum())
    return dict(gated_genes=int(len(tested)), sig_requested=int(want.sum()), sig_realized=int(called.sum()),
                sig_true_positive=tp, sig_recall=tp / max(int(want.sum()), 1), sig_precision=tp / max(int(called.sum()), 1),
                _wrong=np.flatnonzero(called != want))

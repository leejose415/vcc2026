"""Pin each gene's total count over one knockdown's generated cells to the generator's expected total.

The scorer's mse and pds use log1p CPM of each knockdown's summed counts (pseudobulk), so random sampling noise in those
totals is pure error; its Wilcoxon DE metrics compare our cells with the real controls cell by cell, so per-cell
variation must stay control-like. pin_totals fixes the totals and keeps the per-cell shape: surplus genes are thinned
binomially, shortfalls are added as a multinomial over cells weighted towards cells that already express the gene
(w_c ~ X[c, g] + alpha * s_c, so zeros can still rise), and a final hypergeometric pass removes any exact remainder.
"""
import numpy as np


def pin_totals(X, expected_totals, rng, alpha=0.05, max_iter=4):
    """X: one knockdown's counts (cells x genes, integer); expected_totals: float per gene; rng: numpy Generator.
    Returns int64 counts whose column sums equal max(0, round(expected_totals)) exactly."""
    X = np.array(X, dtype=np.int64, copy=True)
    T = np.maximum(0, np.rint(np.asarray(expected_totals, float))).astype(np.int64)
    if T.shape != (X.shape[1],): raise ValueError(f'expected_totals has shape {T.shape}, X has {X.shape[1]} genes')
    lib = X.sum(1).astype(float); mean = lib.mean()
    s = lib / mean if mean > 0 else np.ones(X.shape[0])
    for _ in range(max_iter):
        S = X.sum(0); D = T - S
        if not D.any(): break
        over = np.flatnonzero(D < 0)
        if len(over):                                             # binomial thinning never goes negative
            X[:, over] = rng.binomial(X[:, over], (T[over] / S[over])[None, :])
        S = X.sum(0); D = T - S
        under = np.flatnonzero(D > 0)
        if len(under):
            W = X[:, under].astype(float) + alpha * s[:, None]
            W /= W.sum(0, keepdims=True)
            X[:, under] += rng.multinomial(D[under], W.T).T
    S = X.sum(0)
    for g in np.flatnonzero(S > T):                               # exact: remove the remaining surplus without replacement
        X[:, g] -= rng.multivariate_hypergeometric(X[:, g], int(S[g] - T[g]))
    if not np.array_equal(X.sum(0), T) or X.min(initial=0) < 0:
        raise AssertionError('pin_totals did not reach the targets exactly')
    return X

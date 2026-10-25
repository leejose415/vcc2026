"""Unit tests for the v5 builder's scorer-matching pieces: the MWU kernel against scipy (tie-heavy counts, all-zero and
constant genes) and the CPM against scanpy's normalize_total, bit for bit.
Run from ~/vcc: python -m pytest tests/test_effects_v5.py -q"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from build_effects_v5 import control_index, cpm_csr, mwu_vs_control, quantile_p           # noqa: E402


def counts(rng, n, G):
    lib = rng.choice([800, 1000, 1000, 1200], n)                   # repeated libraries -> nonzero ties within groups
    mu = np.exp(rng.uniform(np.log(1e-3), np.log(5), G)) / 1000
    X = rng.poisson(lib[:, None] * mu[None, :])
    X[:, 0] = 0                                                     # all-zero gene
    return X


def test_mwu_matches_scipy():
    rng = np.random.default_rng(1)
    G = 300; ctrl = counts(rng, 2000, G); pert = counts(rng, 150, G)
    pert[:, 1:20] = rng.poisson(3, (150, 19))                       # shifted genes
    full = np.vstack([ctrl, pert])
    csr = sp.csr_matrix(full)
    C, _ = cpm_csr(csr.data, csr.indices, csr.indptr, G, np.ones(G, bool))
    C = C.tocsr(); Cc, Pc = C[:2000], C[2000:].toarray()
    cols = np.arange(G)
    ptr, vals, cnz, tie = control_index(Cc, cols)
    p = mwu_vs_control(np.ascontiguousarray(np.log1p(Pc).T), ptr, vals, cnz, tie, 2000)
    # On log1p(CPM), the values pdex ranks. (On raw CPM the libraries here, 800 / 1000 / 1200, make coincidental equal
    # ratios such as 4/800 = 5/1000 that log1p merges or splits at the last bit, so the tie terms differ slightly; real
    # libraries rarely do this -- the build's own scipy check runs on CPM.)
    L, Lc = np.log1p(Pc), np.log1p(Cc.toarray())
    ref = mannwhitneyu(L, Lc, use_continuity=True, alternative='two-sided', method='asymptotic', axis=0).pvalue
    ref = np.where(np.isnan(ref), 1.0, ref)                          # scipy: NaN for a constant gene; pdex / gpudge: 1
    assert np.allclose(p, ref, rtol=1e-10, atol=1e-300)
    from numba_mwu import mannwhitneyu_sparse                       # pdex's own kernel
    assert np.array_equal(p, np.asarray(mannwhitneyu_sparse(sp.csr_matrix(L), sp.csr_matrix(Lc)).pvalue))
    assert p[0] == 1.0 and np.median(p[1:20]) < 1e-6


def test_quantiles_reproduce_kernel():
    """The sidecar form (unique sorted nonzero control values + multiplicities + zeros) gives the kernel's p exactly."""
    rng = np.random.default_rng(3)
    G = 200; full = np.vstack([counts(rng, 1500, G), counts(rng, 120, G)]); csr = sp.csr_matrix(full)
    C, _ = cpm_csr(csr.data, csr.indices, csr.indptr, G, np.ones(G, bool)); C = C.tocsr()
    Cc, L = C[:1500], np.log1p(C[1500:].toarray())
    cols = np.flatnonzero(np.asarray(Cc.mean(0)).ravel() > 0)
    ptr, vals, cnz, tie = control_index(Cc, cols)
    p = mwu_vs_control(np.ascontiguousarray(L[:, cols].T), ptr, vals, cnz, tie, 1500)
    # equal to rounding: the kernel's 2 * 0.5 * erfc and quantile_p's erfc differ only in the subnormal range (< 1e-308)
    q = np.array([quantile_p(L[:, cols[g]], *np.unique(vals[ptr[g]:ptr[g + 1]], return_counts=True), int(cnz[g]), 1500)
                  for g in range(len(cols))])
    normal = p > 1e-300
    assert np.array_equal(q[normal], p[normal]) and (q[~normal] < 1e-300).all()


def test_cpm_matches_scanpy():
    import anndata as ad
    import scanpy as sc
    rng = np.random.default_rng(2)
    X = sp.csr_matrix(counts(rng, 500, 400).astype(np.float64))
    mine, _ = cpm_csr(X.data, X.indices, X.indptr, X.shape[1], np.ones(X.shape[1], bool))
    a = ad.AnnData(X=X.copy()); sc.pp.normalize_total(a, target_sum=1e6)
    assert np.array_equal(mine.toarray(), a.X.toarray())            # identical floats, so identical ties

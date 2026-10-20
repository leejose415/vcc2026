"""Unit test for pin_totals on synthetic NB counts (400 cells x 18,500 genes).
Run from ~/vcc: python -m pytest tests/test_exact_totals.py -q -s"""
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code/vcc_2026_3c'))
from vcc_exact_totals import pin_totals                                  # noqa: E402


def synthetic(rng, n_cells=400, n_genes=18500, dispersion=2.0):
    means = np.exp(rng.uniform(np.log(1e-3), np.log(50), n_genes))         # counts per cell, ~1e-3 to ~50
    lib = rng.lognormal(0, 0.3, n_cells); lib /= lib.mean()
    mu = lib[:, None] * means[None, :]
    lam = rng.gamma(1 / dispersion, mu * dispersion)                       # gamma-Poisson = NB, dispersion ~2
    return rng.poisson(lam).astype(np.int64), mu.sum(0)


def log_cpm(totals):
    return np.log1p(totals / max(totals.sum(), 1) * 5e4)


def test_pin_totals():
    rng = np.random.default_rng(2026)
    X, expected = synthetic(rng)
    T = np.maximum(0, np.rint(expected)).astype(np.int64)
    t0 = time.time(); Y = pin_totals(X, expected, np.random.default_rng(7)); sec = time.time() - t0
    assert Y.dtype == np.int64 and Y.min() >= 0
    assert np.array_equal(Y.sum(0), T), 'column sums must equal round(expected) exactly'
    mid = (expected / X.shape[0] > 0.5) & (expected / X.shape[0] < 20)    # mid-expression genes
    vr = np.median(Y[:, mid].var(0) / np.maximum(X[:, mid].var(0), 1e-12))
    det_old, det_new = (X[:, mid] > 0).mean(0), (Y[:, mid] > 0).mean(0)
    mse_old = np.mean((log_cpm(X.sum(0)) - log_cpm(expected)) ** 2)
    mse_new = np.mean((log_cpm(Y.sum(0)) - log_cpm(expected)) ** 2)
    print(f'\n  runtime {sec:.2f} s | variance ratio (mid genes, median) {vr:.3f} | detection {det_old.mean():.4f} -> '
          f'{det_new.mean():.4f} | pseudobulk log1p-CPM MSE {mse_old:.3e} -> {mse_new:.3e} ({mse_old / max(mse_new, 1e-300):.0f}x)')
    assert 0.9 < vr < 1.1
    assert abs(det_new.mean() - det_old.mean()) < 0.01
    assert mse_new < mse_old / 100
    assert sec < 2.0

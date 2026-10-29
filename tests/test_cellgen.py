"""Cell generator (code/cellgen/cell_generator.py) on a synthetic control sidecar: its rank test is the v5 one, means
are exact, calls are realized, output is deterministic given the seed."""
import sys
from pathlib import Path

import numpy as np
import pytest

V = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(V / 'code/cellgen')); sys.path.insert(0, str(V / 'scripts'))
import cell_generator as C                                                                  # noqa: E402


@pytest.fixture(scope='module')
def sidecar(tmp_path_factory):
    """Controls: 3,000 cells x 200 genes (NB counts) plus an untested remainder that brings the library to ~20k."""
    rng = np.random.default_rng(7); n0, G = 3000, 200
    mean = np.exp(rng.uniform(np.log(0.3), np.log(40), G)); r = rng.uniform(1, 8, G)
    counts = rng.negative_binomial(r, r / (r + mean), size=(n0, G))
    lib = counts.sum(1) + rng.integers(15000, 25000, n0)
    x = np.log1p(counts * 1e6 / lib[:, None])
    ptr, vals, cnts, zeros = [0], [], [], []
    for g in range(G):
        v, c = np.unique(x[x[:, g] > 0, g], return_counts=True)
        vals.append(v); cnts.append(c.astype(np.int32)); zeros.append(int((x[:, g] == 0).sum())); ptr.append(ptr[-1] + len(v))
    path = tmp_path_factory.mktemp('cq') / 'syn_control_quantiles.npz'
    np.savez_compressed(path, panel_col=np.arange(G, dtype=np.int32), ptr=np.array(ptr, np.int64), values=np.concatenate(vals),
                        counts=np.concatenate(cnts), zeros=np.array(zeros, np.int32), n_controls=np.int64(n0))
    m0 = np.expm1(x).mean(0)
    m0 = np.append(m0, 1e6 - m0.sum())                  # the untested remainder: the means sum to 1e6, as on the panel
    return path, m0, x


def test_rank_p_is_the_v5_test(sidecar):
    import build_effects_v5 as B
    path, m0, _ = sidecar; cq = C.ControlQuantiles(path); rng = np.random.default_rng(1)
    q = dict(np.load(path))
    for i in rng.choice(len(cq.col), 25, replace=False):
        x = np.log1p(rng.poisson(rng.uniform(0.1, 5), 400) * 50.0)
        a, b = q['ptr'][i], q['ptr'][i + 1]
        ref = B.quantile_p(x, q['values'][a:b], q['counts'][a:b], int(q['zeros'][i]), int(q['n_controls']))
        assert C.rank_p(x, cq.gene(i), cq.n0) == pytest.approx(ref, rel=1e-12, abs=1e-300)


def test_means_calls_and_determinism(sidecar):
    path, m0, _ = sidecar; rng = np.random.default_rng(2); G = len(m0)
    sig = np.zeros(G, bool); sig[rng.choice(G - 1, 40, replace=False)] = True     # DE genes move, the others barely
    lfc = np.where(sig, rng.choice([-1, 1], G) * rng.normal(1.0, 0.3, G), rng.normal(0, 0.08, G))
    mu = C.target_means(lfc, m0); lfc[-1] = np.log2((1e6 - mu[:-1].sum()) / m0[-1])   # a real LFC vector sums to 1e6
    a = C.generate_cells(lfc, sig, path, m0, n_cells=400, library_size=20000, seed=11)
    b = C.generate_cells(lfc, sig, path, m0, n_cells=400, library_size=20000, seed=11)
    assert (a.X != b.X).nnz == 0                                                  # deterministic given the seed
    X = a.X.toarray(); mu = C.target_means(lfc, m0) * 20000 / 1e6
    assert np.abs(X.mean(0) - mu).max() <= 1 / 400 + 1e-9                          # each gene's total is exact
    assert X.dtype.kind in 'iu' and X.min() >= 0
    r = a.uns['cellgen']
    assert r['sig_recall'] >= 0.95 and r['sig_precision'] >= 0.95, r


def test_null_input_calls_nothing(sidecar):
    path, m0, _ = sidecar; G = len(m0)
    a = C.generate_cells(np.zeros(G), np.zeros(G, bool), path, m0, seed=3)
    assert a.uns['cellgen']['sig_realized'] == 0


def test_generator_is_frozen():
    """cell_generator.py is frozen (code/cellgen/FROZEN.md): its bytes must match the recorded hash."""
    import hashlib, re
    record = (V / 'code/cellgen/FROZEN.md').read_text()
    want = re.search(r'sha256 ([0-9a-f]{64})', record).group(1)
    assert hashlib.sha256((V / 'code/cellgen/cell_generator.py').read_bytes()).hexdigest() == want

"""Shared helpers for the VCC 2026 meta-analysis pipeline.

Imported by:
  vcc2_meta.py      : read_names, save_json
  vcc2_generate.py  : seed_for, save_json, expression_gate, int64_indptr, validate_counts, alignment, moments
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

CONTROL = "non-targeting"


def read_names(path, col):
    """One column of a CSV as a numpy string array (first column if `col` is absent)."""
    df = pd.read_csv(path)
    return np.asarray(df[col] if col in df.columns else df.iloc[:, 0], dtype=str)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(f"not JSON serialisable: {type(o)}")


def save_json(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=1, default=_json_default))


def seed_for(seed, *keys):
    """Deterministic seed from (seed, keys...): stable across runs and machines (unlike Python's hash())."""
    h = hashlib.sha256(repr((int(seed),) + tuple(str(k) for k in keys)).encode()).digest()
    return int.from_bytes(h[:8], "little")


def expression_gate(mean_cpm):
    """Strict scorer boundary, using arithmetic mean per-cell CPM in real controls."""
    mean = np.asarray(mean_cpm, dtype=float)
    return np.isfinite(mean) & (mean > 5.0)


def int64_indptr(matrix):
    """CSR with an int64 row pointer.

    An h5ad X that is written once and then grown with append() keeps the dtype of its
    FIRST block's indptr. scipy uses int32 for small blocks, which caps the finished file at
    2**31 stored values (merged H1 has ~3.5e9; a 360k-cell submission is ~2.2e9).
    """
    m = matrix.tocsr()
    m.indptr = m.indptr.astype(np.int64)
    return m


def validate_counts(matrix, label="counts"):
    values = matrix.data if sp.issparse(matrix) else np.asarray(matrix)
    if not (np.isfinite(values).all() and (values >= 0).all() and (values == np.floor(values)).all()):
        raise ValueError(f"{label}: expected finite, nonnegative integer raw counts")


def alignment(source_names, genes):
    """For each of OUR genes, its column in the source (-1 = the source did not measure it)."""
    pos = {g: i for i, g in enumerate(np.asarray(source_names, dtype=str))}
    idx = np.array([pos.get(g, -1) for g in np.asarray(genes, dtype=str)], dtype=np.int64)
    if (idx >= 0).sum() == 0:
        raise ValueError("no genes in common with the provided gene list")
    return idx


def _matrix(a, layer):
    if layer is None:
        return a.X
    if layer in a.layers:
        return a.layers[layer]
    raise KeyError(f"layer {layer!r} not found")


def moments(a, rows, idx, layer=None, block=512):
    """Per-gene statistics over `rows` (cells), streamed `block` cells at a time. Library size = sum over OUR genes.
    Returns arrays aligned to our gene list (NaN where the source lacks the gene)."""
    X = _matrix(a, layer); G = len(idx); ok = idx >= 0; cols = idx[ok]
    s_cnt = np.zeros(ok.sum()); s_cpm = np.zeros(ok.sum()); s_cpm2 = np.zeros(ok.sum()); s_det = np.zeros(ok.sum())
    libs = []
    for b in range(0, len(rows), block):
        r = np.asarray(rows[b:b + block])
        M = X[r]
        M = (M.tocsr() if sp.issparse(M) else sp.csr_matrix(np.asarray(M)))[:, cols].astype(np.float64)
        validate_counts(M)
        lib = np.asarray(M.sum(1)).ravel(); libs.append(lib)
        C = sp.diags(1e6 / np.maximum(lib, 1)) @ M
        s_cnt += np.asarray(M.sum(0)).ravel(); s_cpm += np.asarray(C.sum(0)).ravel()
        s_cpm2 += np.asarray(C.multiply(C).sum(0)).ravel(); s_det += np.asarray((M > 0).sum(0)).ravel()
    n = len(rows); lib = np.concatenate(libs) if libs else np.array([0.0])

    def full(v):
        out = np.full(G, np.nan, dtype=np.float32); out[ok] = v; return out
    if n == 0:
        raise ValueError("no cells to profile")
    mean = s_cpm / n; var = np.maximum(s_cpm2 / n - mean ** 2, 0)
    pooled = s_cnt / max(s_cnt.sum(), 1) * 1e6
    stats = dict(mean_cpm=np.where(ok, 0.0, np.nan), pooled_cpm=full(pooled), bulk_cpm=full(pooled), log_mean=full(np.log1p(pooled)),
                 detection=full(s_det / n), var_cpm=full(var), se_mean=full(np.sqrt(var / n)),
                 n_cells=np.int64(n), median_library=np.float64(np.median(lib)), measured=ok.copy(),
                 legacy_baseline=np.bool_(False))
    stats["mean_cpm"][ok] = mean  # float64 preserves the strict 5-CPM boundary
    return stats
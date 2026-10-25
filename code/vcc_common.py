"""Shared helpers for the VCC 2026 v3 pipeline.

Imported by:
  vcc_context_profile_v3.py : read_names, alignment, select_controls, moments, save_profile
  vcc_transfer_model_v3.py  : Effects, PROGRAMS, load_profile, read_names, save_json, seed_for, shrink
  vcc_build_effects_v3.py   : everything below (writes the formats the other two read)

On-disk formats (one folder per public context, written by vcc_build_effects_v3.py)
  <dir>/lfc.npy              knockdowns x genes, float16 (memory-mapped), log2 fold change vs that context's controls
  <dir>/se.npy               same shape, standard error of the LFC (NaN = unknown)
  <dir>/targets.npy          knockdown names (row order)
  <dir>/genes.npy            gene names (column order)
  <dir>/meta.json            how the numbers were made (se_kind, lfc definition, legacy flag, ...)
  <dir>/control_profile.npz  baseline profile of the context's controls (same format as a profile .npz)
  <dir>/context.json         suggested name / background / platform for contexts.json
  A context is "legacy" when its SE is missing (no honest uncertainty); loading it then needs allow_legacy=True.

Profile .npz (written by save_profile): genes, meta (json string), and per-gene arrays aligned to `genes`:
  mean_cpm   mean of per-cell CPM (the scorer's > 5 CPM gene gate uses this)
  pooled_cpm ratio-of-sums CPM (sum of counts / total * 1e6), comparable across single-cell AND pseudobulk contexts
  log_mean   log1p(pooled_cpm)  <- the baseline the transfer model compares across contexts
  detection, var_cpm, se_mean (NaN when not computable, e.g. pseudobulk-only contexts)
  plus program_names / program_levels (mean log1p pooled CPM over each program's genes), n_cells, legacy_baseline.
  Genes a context did not measure are NaN.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

CONTROL = "non-targeting"

# Program gene sets. 'prefix:' entries are expanded against the gene list (ribosomal proteins, mitochondrial genes).
PROGRAMS = {
    "cc_S": "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 HELLS RFC2 RPA2 NASP RAD51AP1 GMNN "
            "WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 "
            "CHAF1B BRIP1 E2F8".split(),
    "cc_G2M": "HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 SMC4 CCNB2 CKAP2L "
              "CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 "
              "CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA".split(),
    "stress": "HSPA1A HSPA1B HSPA8 HSPH1 DNAJB1 HSP90AA1 FOS JUN EGR1 ATF3 IER2 DUSP1 DDIT3 ATF4".split(),
    "interferon": "ISG15 IFI6 IFIT1 IFIT2 IFIT3 MX1 OAS1 IFI44L STAT1 IRF7 IFITM1 BST2".split(),
    "hypoxia_glyc": "VEGFA CA9 BNIP3 PGK1 LDHA SLC2A1 ENO1 ALDOA NDRG1 P4HA1 PFKFB3 HK2".split(),
    "biogenesis": "NPM1 NCL NOP56 NOP58 FBL RRP9 PA2G4 CAD ODC1 DKC1 NOP16 POLR1B".split(),
    "ribosome": ["prefix:RPL", "prefix:RPS"],
    "mito": ["prefix:MT-"],
}
_PREFIX_EXCLUDE = ("RPL3L", "RPS4Y")


# ---------------------------------------------------------------- small utilities
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


SE_EPS = 1e-6


def valid_se(se):
    se = np.asarray(se, dtype=float)
    return np.isfinite(se) & (se > SE_EPS)


def estimate_prior(lfc, se, floor=0.01):
    """Normal-prior variance: observed variance minus mean measurement variance."""
    lfc, se = np.asarray(lfc, float), np.asarray(se, float)
    good = np.isfinite(lfc) & valid_se(se)
    if good.sum() < 2:
        return floor
    return float(max(np.var(lfc[good], ddof=1)-np.mean(se[good]**2), floor))


def shrink(lfc, se, cap=8.0, prior_variance=None, missing_se_factor=0.5):
    """Normal-prior posterior mean; zero/tiny/nonfinite SE has missing-SE fallback.

    Pass a context's prior_variance for consistent shrinkage across gene chunks.
    The fallback row estimate is intended only for standalone diagnostics.
    """
    lfc = np.asarray(lfc, dtype=float)
    se = np.asarray(se, dtype=float) if se is not None else np.full_like(lfc, np.nan)
    prior = estimate_prior(lfc, se) if prior_variance is None else float(prior_variance)
    if not np.isfinite(prior) or prior <= 0:
        raise ValueError("prior variance must be finite and positive")
    factor = np.full_like(lfc, missing_se_factor)
    good = valid_se(se)
    factor[good] = prior/(prior+se[good]**2)
    return np.where(np.isfinite(lfc), np.clip(lfc*factor, -cap, cap), np.nan).astype(np.float32)


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


# ---------------------------------------------------------------- gene alignment and programs
def alignment(source_names, genes):
    """For each of OUR genes, its column in the source (-1 = the source did not measure it)."""
    pos = {g: i for i, g in enumerate(np.asarray(source_names, dtype=str))}
    idx = np.array([pos.get(g, -1) for g in np.asarray(genes, dtype=str)], dtype=np.int64)
    if (idx >= 0).sum() == 0:
        raise ValueError("no genes in common with the provided gene list")
    return idx


def program_gene_index(genes):
    """{program: indices into genes}."""
    genes = np.asarray(genes, dtype=str); pos = {g: i for i, g in enumerate(genes)}; out = {}
    for name, members in PROGRAMS.items():
        idx = []
        for m in members:
            if m.startswith("prefix:"):
                p = m.split(":", 1)[1]
                idx += [i for i, g in enumerate(genes) if g.startswith(p) and not g.startswith(_PREFIX_EXCLUDE)]
            elif m in pos:
                idx.append(pos[m])
        out[name] = np.array(sorted(set(idx)), dtype=np.int64)
    return out


def program_levels(pooled_cpm, genes):
    """Mean log1p(pooled CPM) over each program's measured genes (NaN if < 3 measured)."""
    lv = []
    for name, idx in program_gene_index(genes).items():
        v = np.log1p(pooled_cpm[idx]) if len(idx) else np.array([])
        v = v[np.isfinite(v)]
        lv.append(float(v.mean()) if len(v) >= 3 else np.nan)
    return np.array(list(PROGRAMS)), np.array(lv, dtype=np.float32)


# ---------------------------------------------------------------- profiles from cells
def select_controls(a, obs_col, control=CONTROL, all_controls=False, max_cells=20000, seed=0):
    """Row indices of control cells. If `obs_col` is absent every row is treated as a control (e.g. context_A.h5ad).
    Without all_controls, at most `max_cells` controls (random, seeded) are used."""
    if obs_col in a.obs.columns:
        rows = np.flatnonzero(a.obs[obs_col].astype(str).values == control)
    else:
        rows = np.arange(a.n_obs)
    if len(rows) == 0:
        raise ValueError(f"no rows with {obs_col} == {control!r}")
    if not all_controls and len(rows) > max_cells:
        rows = np.sort(np.random.default_rng(seed).choice(rows, max_cells, replace=False))
    return np.sort(rows)


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
    stats = dict(mean_cpm=full(mean), pooled_cpm=full(pooled), bulk_cpm=full(pooled), log_mean=full(np.log1p(pooled)),
                 detection=full(s_det / n), var_cpm=full(var), se_mean=full(np.sqrt(var / n)),
                 n_cells=np.int64(n), median_library=np.float64(np.median(lib)), measured=ok.copy(),
                 legacy_baseline=np.bool_(False))
    return stats


def finalise_profile(stats, genes):
    """Add program levels and the scorer gate flag."""
    names, lv = program_levels(stats["pooled_cpm"], genes)
    stats["program_names"] = names; stats["program_levels"] = lv
    stats["scored"] = np.nan_to_num(stats["mean_cpm"], nan=0) > 5
    return stats


def save_profile(output, genes, stats, meta=None):
    """Write a profile .npz (program levels added if missing)."""
    stats = dict(stats)
    if "pooled_cpm" in stats:
        stats["bulk_cpm"] = np.asarray(stats["pooled_cpm"]).copy()
        stats["baseline_convention"] = np.array("ratio_of_sums_cpm")
    if "program_levels" not in stats:
        stats = finalise_profile(stats, genes)
    out = Path(output); out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, genes=np.asarray(genes, dtype=str), meta=json.dumps(meta or {}, default=_json_default), **stats)
    return out


def load_profile(path, genes):
    """Profile re-aligned to `genes`. Per-gene arrays: NaN (or False) where missing. Returns a plain dict."""
    z = np.load(path, allow_pickle=False)
    pg = z["genes"].astype(str); src = alignment(pg, genes); ok = src >= 0
    out = {}
    for k in z.files:
        if k in ("genes", "meta"):
            continue
        v = z[k]
        if v.ndim == 1 and len(v) == len(pg) and k not in ("program_names", "program_levels", "baseline_program_levels", "baseline_program_coverage", "control_obs_names"):
            if v.dtype == bool:
                w = np.zeros(len(genes), bool); w[ok] = v[src[ok]]
            else:
                w = np.full(len(genes), np.nan, dtype=np.float32); w[ok] = v[src[ok]]
            out[k] = w
        else:
            out[k] = v.item() if v.ndim == 0 else v
    out["genes"] = np.asarray(genes, dtype=str)
    out["meta"] = json.loads(str(z["meta"])) if "meta" in z.files else {}
    out["legacy_baseline"] = bool(out.get("legacy_baseline", False))
    if "program_levels" not in out:
        out["program_names"], out["program_levels"] = program_levels(np.expm1(out["log_mean"]), genes)
    out["program_levels"] = np.nan_to_num(np.asarray(out["program_levels"], dtype=np.float32), nan=0.0)
    for k in ("log_mean", "mean_cpm"):
        if k not in out:
            raise KeyError(f"profile {path} lacks {k}")
    return out


# ---------------------------------------------------------------- effect tables
class Effects:
    """One public context's effect table: LFC and SE per (knockdown, gene), memory-mapped from disk.
    get(kd, gene_idx) -> (lfc, se) float32 arrays for OUR gene indices (NaN = not measured)."""

    def __init__(self, path, genes, allow_legacy=False):
        p = Path(path); self.path = p
        self.meta = json.loads((p / "meta.json").read_text()) if (p / "meta.json").exists() else {}
        self.legacy = bool(self.meta.get("legacy", False))
        if self.legacy and not allow_legacy:
            raise ValueError(f"{p} is a legacy effect table (no usable SE); pass --allow-legacy-effects to use it")
        self.lfc = np.load(p / "lfc.npy", mmap_mode="r")
        self.se = np.load(p / "se.npy", mmap_mode="r") if (p / "se.npy").exists() else None
        self.targets = [str(t) for t in np.load(p / "targets.npy")]
        self.pos = {t: i for i, t in enumerate(self.targets)}
        self.col = alignment(np.load(p / "genes.npy").astype(str), genes)    # our gene -> stored column
        self.n_genes = len(genes)

    @property
    def prior_variance(self):
        if not hasattr(self, "_prior_variance"):
            n = total = square = noise = 0.0
            if self.se is not None:
                for start in range(0, len(self.targets), 64):
                    l = np.asarray(self.lfc[start:start+64], float)
                    se = np.asarray(self.se[start:start+64], float)
                    good = np.isfinite(l) & valid_se(se)
                    v, e = l[good], se[good]
                    n += len(v); total += v.sum(); square += np.dot(v, v); noise += np.dot(e, e)
            self._prior_variance = max((square-total*total/n)/max(n-1, 1)-noise/n, 0.01) if n > 1 else 0.01
        return self._prior_variance

    def prior_for(self, kd, floor=0.01):
        """Normal-prior variance for ONE knockdown: var(LFC) - mean(SE^2) over that row's genes.

        A table-wide prior is dominated by the many weak knockdowns, so it shrank the strong
        ones' clearly real effects by ~30-45%. Per row, a weak knockdown still sits at the floor
        (identical to the table prior, which is itself at the floor for Orion/T cells), while a
        strong one keeps its signal. Estimated once from the FULL row, so every gene chunk of a
        knockdown is shrunk consistently. Tables without SE fall back to the table prior.
        """
        cache = self.__dict__.setdefault("_row_priors", {})
        if kd not in cache:
            r = self.pos.get(kd)
            if r is None or self.se is None:
                cache[kd] = self.prior_variance
            else:
                cache[kd] = estimate_prior(np.asarray(self.lfc[r], float), np.asarray(self.se[r], float), floor=floor)
        return cache[kd]

    def get(self, kd, gi):
        gi = np.asarray(gi, dtype=np.int64).ravel()
        lfc = np.full(len(gi), np.nan, np.float32); se = np.full(len(gi), np.nan, np.float32)
        r = self.pos.get(kd)
        if r is None:
            return lfc, se
        c = self.col[gi]; ok = c >= 0
        if ok.any():
            row = np.asarray(self.lfc[r], dtype=np.float32); lfc[ok] = row[c[ok]]
            if self.se is not None:
                srow = np.asarray(self.se[r], dtype=np.float32); se[ok] = srow[c[ok]]
        return lfc, se

    def __repr__(self):
        return f"Effects({self.path.name}: {len(self.targets):,} knockdowns, legacy={self.legacy})"


def write_effects(out, genes, targets, rows_iter, meta, n_rows):
    """Stream (k, lfc_row, se_row) into an effect-table folder. Rows arrive in `targets` order."""
    from numpy.lib.format import open_memmap
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    L = open_memmap(out / "lfc.npy", mode="w+", dtype=np.float16, shape=(n_rows, len(genes)))
    S = open_memmap(out / "se.npy", mode="w+", dtype=np.float16, shape=(n_rows, len(genes)))
    done = 0
    for k, lfc, se in rows_iter:
        L[k] = np.clip(lfc, -60000, 60000); S[k] = np.clip(se, 0, 60000); done += 1
    L.flush(); S.flush(); del L, S
    np.save(out / "targets.npy", np.asarray(targets, dtype=str)); np.save(out / "genes.npy", np.asarray(genes, dtype=str))
    save_json(out / "meta.json", dict(meta, n_targets=int(n_rows), n_genes=int(len(genes))))
    assert done == n_rows, (done, n_rows)
    return out

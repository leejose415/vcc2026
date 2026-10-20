"""Generate counts with projected composition, metadata cell counts/library sizes,
state-dependent overdispersed induction and reproducible common control sampling.

Finishing: a DE call, then one global alpha, an optional detection floor, the > 5 CPM gate in the line's own controls,
and the knocked-down gene applied separately. Two calling rules:
  z     keep |delta| > z* x SE_pred and |delta| > a minimum size (the vcc_2026_2 submission). With gate_from (e.g.
        're'), the keep-mask comes from that array and its SE (|re| > z* x re_se) and is applied to the selected
        array's values, so a ridge output without its own SE can borrow Model 2's significance call
  topn  per knockdown keep the N genes with the largest |delta| / SE_pred among genes > 5 CPM (optionally only genes
        whose backgrounds agree in sign, or measured in >= min_bg backgrounds); every other gene is 0. With gate_from, that array's |delta| / SE ranks the
        genes and the selected array supplies the values (e.g. Model 2 picks, the corrected ridge output is used)
Detection floor: each called gene is raised to the smallest log2 change a 400-cell test against the line's controls
would detect at z_req, so that the genes we call actually show up as DE in the generated cells.

Pseudo-model `nochange`: zero effects AND no knockdown (the target gene is left at its control level), the true
no-change arm; the `zero` array instead keeps the measured knockdown.

`benchmark` scores arms on the exam (H1-held-out validation) with official cell-eval2 (vcc2026), streaming every
file one target at a time. The real side, the ceiling and the competition bundle (generic baseline + replicate anchor)
are built once, and can be reused by a later stage (--reference). Every arm is scored with the leaderboard formula
(scaled, clamped, averaged: scaled_scores.csv) and selection uses that overall; raw metrics and ranks are kept.
`oracle` is a separate generator ceiling check and never selects submission knobs. With --variants it scores several
arms on the same real cells (method `oracle` or `nochange`, any alpha and cell rule); with --scale it also builds the
competition bundle on those cells, so the ceiling is reported on the leaderboard scale.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import shlex
import subprocess
import shutil
import copy
from vcc2_common import int64_indptr, seed_for, validate_counts, save_json, expression_gate, moments, alignment

try:
    from anndata.io import read_elem, write_elem, sparse_dataset
except ImportError:  # older anndata
    from anndata.experimental import read_elem, write_elem, sparse_dataset

VCC2026_METRICS = {
    'pds_cosine': 'max', 'expr_mse_unbiased_capped_norm': 'min',
    'de_wilcoxon_lfc_nmae': 'min', 'de_wilcoxon_direction_fidelity_yield_raw': 'max',
    'de_wilcoxon_direction_reach_raw': 'max', 'de_wilcoxon_sig_jaccard': 'max',
}
CONTROL = "non-targeting"
SCHEMA = "vcc_2026_3_ridge_v1"
CC_GENES = ("MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 HELLS RFC2 RPA2 NASP RAD51AP1 GMNN "
            "WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 "
            "CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 SMC4 "
            "CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 CDC20 TTK CDC25C KIF2C RANGAP1 "
            "NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA").split()


class Cfg:
    def __init__(self, a):
        self.alpha = getattr(a, "alpha", 1.0); self.max_lfc = 8.0; self.min_lfc = 0.0
        self.kd_clip = (0.01, 0.9); self.default_kd = 0.30; self.n_cells = getattr(a, "n_cells", 400)
        self.composition = "projection"; self.comp_thr = getattr(a, "comp_thr", 0.02)
        self.pi = getattr(a, "pi", 0.0); self.up_w = getattr(a, "up_w", 0.75)
        self.tilt = getattr(a, "tilt", False); self.match_pool = getattr(a, "match_pool", False) or self.tilt
        self.seed = getattr(a, "seed", 2026)
        self.z = getattr(a, "z", 0.0); self.min_effect = getattr(a, "min_effect", 0.0)   # z call; 0 and 0 = off
        self.call = getattr(a, "call", "z"); self.top_n = int(getattr(a, "top_n", 1000))
        self.sign_rule = bool(getattr(a, "sign_rule", False)); self.detect_floor = bool(getattr(a, "detect_floor", False))
        self.z_req = float(getattr(a, "z_req", 4.0))
        self.min_bg = int(getattr(a, "min_bg", 1))                # top-N: only genes measured in >= min_bg backgrounds
        self.gate_from = getattr(a, "gate_from", "") or ""
        self.metadata = getattr(a, "metadata", True)
        self.library_size = None
        self.exact_totals = bool(getattr(a, "exact_totals", False))   # pin pseudobulk totals (vcc_exact_totals.py)
        if (not 0 <= self.pi < 1 or not 0 <= self.up_w <= 1 or self.n_cells < 1 or self.alpha < 0
                or self.z < 0 or self.min_effect < 0 or self.call not in ("z", "topn") or self.top_n < 1 or self.z_req <= 0 or self.min_bg < 1
                or (self.gate_from and self.call == "z" and self.z <= 0)):
            raise ValueError("invalid generator settings")

    def as_dict(self):
        # min_bg is left out at its default, so settings written before it existed still describe the same arm.
        return {k: v for k, v in vars(self).items() if not (k == "min_bg" and v == 1) and not (k == "exact_totals" and not v)}


# ---------------------------------------------------------------- inputs
def read_list(path, col):
    return pd.read_csv(path)[col].astype(str).tolist()


def aligned_X(a, genes):
    """cells x our genes, CSR float32 (missing genes = 0, extra genes dropped)."""
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    names = a.var_names.astype(str).tolist()
    if names == list(genes):
        return X.astype(np.float32)
    pos = {g: i for i, g in enumerate(genes)}
    src = np.array([i for i, g in enumerate(names) if g in pos]); dst = np.array([pos[names[i]] for i in src])
    X = X[:, src].tocoo()
    print(f"  re-indexed genes: {len(src):,}/{len(genes):,} of ours present")
    return sp.csr_matrix((X.data.astype(np.float32), (X.row, dst[X.col])), shape=(X.shape[0], len(genes)))


class Line:
    """Control cells of one context plus the statistics the generator needs."""
    def __init__(self, X, guide, genes):
        validate_counts(X, "control cells")
        self.X = X.tocsr().astype(np.float32); self.genes = genes; self.n = X.shape[0]
        if not self.n or self.X.sum() <= 0: raise ValueError("controls must contain cells and positive counts")
        self.lib = np.asarray(self.X.sum(1)).ravel()
        self.cpm_cells = sp.diags(1e6 / np.maximum(self.lib.astype(np.float64), 1)) @ self.X
        self.cpm = np.asarray(self.cpm_cells.mean(0)).ravel()           # mean per-cell CPM (scorer gate uses > 5)
        self.cpm_var = np.maximum(np.asarray(self.cpm_cells.multiply(self.cpm_cells).mean(0)).ravel() - self.cpm ** 2, 0)
        self.n_reference = self.n                                       # control cells the scorer compares against
        self.share = np.asarray(self.X.sum(0)).ravel() / self.X.sum()   # pooled share of molecules
        self.gate = expression_gate(self.cpm)
        count_mean = np.asarray(self.X.mean(0)).ravel()
        count_var = np.maximum(np.asarray(self.X.multiply(self.X).mean(0)).ravel()-count_mean**2, 0)
        self.overdispersion = np.where(count_mean > 0, np.clip((count_var-count_mean)/np.maximum(count_mean**2, 1e-12), 0, 10), 1.0)
        g = np.asarray(guide).astype(str) if guide is not None else np.zeros(self.n, str)
        self.groups = [np.flatnonzero(g == u) for u in np.unique(g)]
        gp = {x: i for i, x in enumerate(genes)}
        self.cc_idx = np.array([gp[x] for x in CC_GENES if x in gp and self.cpm[gp[x]] > 1], int)
        if len(self.cc_idx) >= 5:
            L = self.cpm_cells[:, self.cc_idx].toarray()
            self.cyc = np.log2(L + 1).mean(1)                            # per-cell cycling level (log2 units)
        else:
            self.cyc = None

    @staticmethod
    def from_h5ad(path, genes, control_col="target_gene"):
        a = ad.read_h5ad(path)
        validate_counts(a.X, str(path))
        if control_col and control_col in a.obs:
            a = a[a.obs[control_col].astype(str).values == CONTROL].copy()
        guide = next((a.obs[c].values for c in ("ntc_id", "guide_id", "sgID_AB", "guide") if c in a.obs), None)
        return Line(aligned_X(a, genes), guide, genes)


def load_deltas(path, targets, genes, array="re", fill=0.0):
    """meta prediction dir (genes.npy, targets.npy, <array>.npy) or .npz (delta [+targets, genes]) -> T x G.

    Missing targets/genes and non-finite values become `fill`: 0 (no change) for deltas, inf (unknown) for SEs.
    """
    p = Path(path)
    if p.is_dir():
        D = np.load(p / f"{array}.npy", mmap_mode="r"); t = np.load(p / "targets.npy").astype(str)
        g = np.load(p / "genes.npy").astype(str)
    else:
        z = np.load(p, allow_pickle=False); D = z["delta"]
        t = z["targets"].astype(str) if "targets" in z else np.array(targets)
        g = z["genes"].astype(str) if "genes" in z else np.array(genes)
    out = np.full((len(targets), len(genes)), fill, np.float32)
    tp = {x: i for i, x in enumerate(t)}; gp = {x: i for i, x in enumerate(g)}
    gsrc = np.array([gp.get(x, -1) for x in genes]); ok = gsrc >= 0
    missing = 0
    for i, x in enumerate(targets):
        if x in tp:
            out[i, ok] = np.asarray(D[tp[x]])[gsrc[ok]]
        else:
            missing += 1
    if missing:
        print(f"  WARNING {path}: {missing} targets missing -> no change predicted for them")
    out[~np.isfinite(out)] = fill
    return out


def load_info(path, targets, genes, model, gate_from=""):
    """The calling inputs that exist next to a model's table: se, z (|delta| / SE), n_bg and sign_agree, plus the
    gating array and its SE when another array decides which genes are called (gate_from)."""
    p = Path(path); info = {}
    for key, name, fill in [("se", model + "_se", np.inf), ("z", model + "_z", np.nan), ("n_bg", "n_bg", np.nan),
                            ("sign_agree", "sign_agree", np.nan)]:
        if p.is_dir() and (p / f"{name}.npy").exists(): info[key] = load_deltas(path, targets, genes, name, fill=fill)
    if gate_from:
        for key, name, fill in [("gate_delta", gate_from, 0.0), ("gate_se", gate_from + "_se", np.inf)]:
            if not (p / f"{name}.npy").exists(): raise ValueError(f"gate_from={gate_from!r} needs {name}.npy")
            info[key] = load_deltas(path, targets, genes, name, fill=fill)
    return info


def is_nochange(model):
    return model == "nochange"


def load_effects(path, targets, genes, model, cfg):
    """Delta table of one model, plus the calling inputs (SE, z, n_bg, sign_agree, gate) that exist for it.
    `nochange` is all zeros (and the caller skips the knockdown)."""
    if is_nochange(model): return np.zeros((len(targets), len(genes)), np.float32), {}
    D = load_deltas(path, targets, genes, model)
    info = load_info(path, targets, genes, model, cfg.gate_from)
    if cfg.z > 0 and not cfg.gate_from and "se" not in info: raise ValueError(f"{model}: the z call needs {model}_se.npy")
    return D, info


def kd_from_sources(targets, genes, orion_paths=(), tcell_path=None, tcell_condition="Stim48hr", pseudo=0.5):
    """Knocked-down gene's own log2 drop per target, averaged across sources (M1 rule). NaN if no source has it."""
    gp = {g: i for i, g in enumerate(genes)}; vals = [[] for _ in targets]
    for path in orion_paths:
        with h5py.File(path, "r") as f:
            obs = read_elem(f["obs"]); names = list(obs.index.astype(str)); n = obs["n_cells"].values.astype(float)
            var = read_elem(f["var"]); vg = list(var.index.astype(str))
            vpos = {g: i for i, g in enumerate(vg)}
            c = names.index(CONTROL)
            want = [(i, t) for i, t in enumerate(targets) if t in names and t in vpos]
            if not want:
                continue
            rows = sorted({c} | {names.index(t) for _, t in want})
            S = f["layers/sum_cpm"][rows, :]; rpos = {r: k for k, r in enumerate(rows)}
            ctrl = S[rpos[c]] / n[c]
            for i, t in want:
                r = names.index(t); col = vpos[t]
                vals[i].append(np.log2((S[rpos[r], col] / n[r] + pseudo) / (ctrl[col] + pseudo)))
    if tcell_path:
        t = ad.read_h5ad(tcell_path)
        if "condition" in t.obs:
            t = t[t.obs["condition"].astype(str) == tcell_condition]
        tg = t.obs["target_gene"].astype(str).values; tv = {g: i for i, g in enumerate(t.var_names.astype(str))}
        X = np.asarray(t.X)
        for i, x in enumerate(targets):
            r = np.flatnonzero(tg == x)
            if len(r) and x in tv and np.isfinite(X[r[0], tv[x]]):
                vals[i].append(float(X[r[0], tv[x]]))
    return np.array([np.mean(v) if v else np.nan for v in vals])


def kd_residual(log2_drop, cfg):
    r = np.where(np.isfinite(log2_drop), 2.0 ** np.nan_to_num(log2_drop), cfg.default_kd)
    return np.clip(r, *cfg.kd_clip)


# ---------------------------------------------------------------- core generator
def call_top_n(D, line, own, cfg, info):
    """Per knockdown, keep the top_n genes by |delta| / SE among genes > 5 CPM (never the knocked-down gene).
    With sign_rule, a gene measured in >= 2 backgrounds is a candidate only if they all agree in sign.
    With gate_from, the gating array's |delta| / SE ranks the genes; the values stay those of D."""
    if cfg.gate_from:
        with np.errstate(divide="ignore", invalid="ignore"): z = np.abs(info["gate_delta"]) / info["gate_se"]
    else:
        z = info.get("z")
        if z is None:
            if "se" not in info: raise ValueError("the top-N call needs re_z or the prediction SE")
            with np.errstate(divide="ignore", invalid="ignore"): z = np.abs(D) / info["se"]
    cand = line.gate[None, :] & np.isfinite(D) & (D != 0) & np.isfinite(z)
    rows = np.flatnonzero(own >= 0); cand[rows, own[rows]] = False
    if cfg.sign_rule:
        if "n_bg" not in info or "sign_agree" not in info: raise ValueError("the sign rule needs n_bg and sign_agree")
        with np.errstate(invalid="ignore"):
            cand &= ~((info["n_bg"] >= 2) & (info["sign_agree"] < 1))
    if cfg.min_bg > 1:
        if "n_bg" not in info: raise ValueError("min_bg needs n_bg")
        with np.errstate(invalid="ignore"):
            cand &= info["n_bg"] >= cfg.min_bg
    score = np.where(cand, z, -np.inf)
    keep = np.zeros_like(cand)
    n = min(cfg.top_n, D.shape[1])
    top = np.argpartition(-score, n - 1, axis=1)[:, :n]
    np.put_along_axis(keep, top, True, axis=1)
    return np.where(keep & cand, D, 0).astype(np.float32)


def detection_floor(line, cfg, iters=40):
    """Per gene, the smallest |log2 change| (up, down) at which a cfg.n_cells-cell group would differ from the line's
    controls at z_req: m |f - 1| / sqrt(v f^2 / n_cells + v / n_reference) >= z_req, with f = 2^lfc and m, v the
    per-cell CPM mean and variance of the controls (the knockdown keeps the controls' CV). Unreachable -> max_lfc."""
    key = (cfg.n_cells, line.n_reference, cfg.z_req, cfg.max_lfc)
    cache = getattr(line, "_floor", {})
    if key in cache: return cache[key]
    m, v = line.cpm.astype(float), line.cpm_var.astype(float)
    def z_of(l):
        f = np.exp2(l)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(v > 0, m * np.abs(f - 1) / np.sqrt(v * f ** 2 / cfg.n_cells + v / line.n_reference), np.inf)
    floors = []
    for end in (cfg.max_lfc, -cfg.max_lfc):                      # up, then down
        lo, hi = np.zeros_like(m), np.full_like(m, end)
        reachable = z_of(hi) >= cfg.z_req
        for _ in range(iters):                                    # bisection on the log2 scale between 0 and the cap
            mid = (lo + hi) / 2; ok = z_of(mid) >= cfg.z_req
            hi = np.where(ok, mid, hi); lo = np.where(ok, lo, mid)
        floors.append(np.where(reachable, np.abs(hi), cfg.max_lfc).astype(np.float32))
    cache[key] = tuple(floors); line._floor = cache
    return cache[key]


def calibrate(D, line, targets, cfg, se=None, info=None):
    """DE call, then the single global alpha, the optional detection floor and a binary responder-gene gate.

    z call: an effect stays at full size if |delta| > min_effect and, when z* > 0, |delta| > z* x SE_pred; else 0.
    topn call: see call_top_n. Unknown SE (inf/NaN) is never convincing. No shrinkage of retained effects.
    Detection floor: each called gene becomes sign(delta) x max(alpha |delta|, floor_g).
    The +/-8 generator limit is a numerical safety cap, not an estimated prior.
    On-target suppression is handled separately by measured guide efficiency.
    """
    info = dict(info or {})
    if se is not None: info.setdefault("se", se)
    D = np.asarray(D, np.float32)
    gp = {g: i for i, g in enumerate(line.genes)}
    own = np.array([gp.get(t, -1) for t in targets])
    if cfg.call == "topn":
        D = call_top_n(D, line, own, cfg, info)
    elif cfg.gate_from:
        # Model 2-gated: genes Model 2 calls (|re| > z* x re_se), with the selected array's values.
        with np.errstate(invalid="ignore"):
            keep = (np.abs(info["gate_delta"]) > cfg.z * np.asarray(info["gate_se"], np.float32)) & (np.abs(D) > cfg.min_effect)
        D = np.where(keep, D, 0)
    elif cfg.z > 0 or cfg.min_effect > 0:
        size = np.abs(D); keep = size > cfg.min_effect
        if cfg.z > 0:
            if "se" not in info: raise ValueError("the SE gate (z* > 0) needs the prediction SE")
            with np.errstate(invalid="ignore"):
                keep &= size > cfg.z * np.asarray(info["se"], np.float32)
        D = np.where(keep, D, 0)
    D = (cfg.alpha * D).astype(np.float32)
    if cfg.detect_floor:
        up, down = detection_floor(line, cfg)
        floor = np.where(D > 0, up[None, :], down[None, :])
        D = np.where(D != 0, np.sign(D) * np.maximum(np.abs(D), floor), 0)
    D = np.clip(D, -cfg.max_lfc, cfg.max_lfc).astype(np.float32)
    D[:, ~line.gate] = 0.0
    rows = np.flatnonzero(own >= 0); D[rows, own[rows]] = 0.0
    return D


def pick_cells(line, d, cfg, rng):
    """Step 3 (+ optional tilt): which control cells become this knockdown's 400 cells."""
    n = cfg.n_cells
    if cfg.tilt and line.cyc is not None:
        dcc = float(np.mean(d[line.cc_idx]))                       # predicted cycling shift (log2)
        v = line.cyc - line.cyc.mean(); s = max(v.std(), 1e-6)
        if abs(dcc) > 0.02:
            lo, hi = -30.0, 30.0
            for _ in range(50):                                     # bisection: weighted mean shift == dcc
                b = (lo + hi) / 2; w = np.exp(np.clip(b * v / s, -50, 50)); m = (w * v).sum() / w.sum()
                lo, hi = (b, hi) if m < dcc else (lo, b)
            w = np.exp(np.clip(b * v / s, -50, 50)); w /= w.sum()
            k = min(n, int((w > 0).sum()))
            return rng.choice(line.n, k, replace=False, p=w) if k == n else rng.choice(line.n, n, replace=True, p=w)
    if n > line.n:
        order = np.concatenate([rng.permutation(g) for g in line.groups])
        return np.concatenate([order, rng.choice(line.n, n-line.n, replace=True)])
    per = int(np.ceil(n / len(line.groups)))
    pick = np.concatenate([rng.choice(v, size=min(per, len(v)), replace=False) for v in line.groups])
    if len(pick) >= n:
        return rng.choice(pick, n, replace=False)
    rest = np.setdiff1d(np.arange(line.n), pick)
    return np.concatenate([pick, rng.choice(rest, n - len(pick), replace=len(rest) < n - len(pick))])


def project_factors(share, delta, j=None, kd=None):
    """Project all-gene factors onto the composition simplex without an LFC floor.

    The realized expectation is log2(projected_factor), not the unconstrained
    delta. Zero-baseline induction uses the same 0.5-CPM pseudocount as labels.
    """
    fac = np.exp2(np.asarray(delta, float))
    if j is not None and kd is not None: fac[j] = kd
    proposed = np.asarray(share, float)*fac
    zero_induced = (share == 0) & (fac > 1)
    proposed[zero_induced] = 0.5e-6*(fac[zero_induced]-1)
    mass = float(proposed.sum())
    if not np.isfinite(mass) or mass <= 0: raise ValueError("invalid projected composition")
    projected = proposed/mass
    factors = np.divide(projected, share, out=np.ones_like(projected), where=share > 0)
    return factors, projected, mass


def generate_one(line, d, j, kd, cfg, rng=None, pick=None, line_key="line", target_key="target"):
    """Integer counts, with stable per-gene streams shared between perturbations."""
    if rng is None: rng = np.random.default_rng(seed_for(cfg.seed, line_key, target_key))
    if pick is None:
        # Same base pool per line; target order does not advance a line-level stream.
        pick = pick_cells(line, d, cfg, np.random.default_rng(seed_for(cfg.seed, line_key, "base_controls")))
    B = line.X[pick].toarray().astype(np.int64); d = np.asarray(d, float).copy()
    if cfg.match_pool:
        sampled = np.asarray(line.cpm_cells[pick].mean(0)).ravel()
        deviation = np.log2((sampled+0.5)/(line.cpm+0.5))
        mask = line.gate.copy()
        if j is not None: mask[j] = False
        d[mask] = np.clip(d[mask]-deviation[mask], -cfg.max_lfc, cfg.max_lfc)
    factors, projected, mass = project_factors(line.share, d, j, kd)
    # Library size is an independent whole-cell scale after composition projection.
    lib_scale = 1.0
    if cfg.library_size is not None:
        lib_scale = float(cfg.library_size)/max(float(np.median(line.lib[pick])), 1)
    positive = line.share > 0
    effective_pi = min(cfg.pi, float(factors[positive].min())) if positive.any() else 0.0
    esc_rng = np.random.default_rng(seed_for(cfg.seed, line_key, "escapers"))
    escapers = esc_rng.random(len(pick)) < effective_pi
    cols = np.flatnonzero((np.abs(factors*lib_scale-1) > 1e-12) | ((line.share == 0) & (projected > 0)))
    E = B.astype(np.float64) if cfg.exact_totals else None    # expected counts; untouched genes keep their base counts
    for g in cols:
        # Fresh per-gene seed gives common random draws without target-order dependence.
        grng = np.random.default_rng(seed_for(cfg.seed, line_key, "gene", line.genes[g]))
        factor = factors[g]
        inside = max((factor-effective_pi)/(1-effective_pi), 0)
        F = np.where(escapers, 1.0, inside)*lib_scale
        original = B[:, g].copy()
        values = grng.binomial(original, np.minimum(F, 1.0))
        up = np.maximum(F-1, 0)
        base = line.lib[pick]*line.share[g]
        if line.share[g] == 0 and projected[g] > 0:
            # Induced genes absent from controls must be able to gain counts.
            intensity = line.lib[pick]*projected[g]*lib_scale
        else:
            expected = cfg.up_w*original+(1-cfg.up_w)*base
            intensity = up*expected
        if E is not None:                                         # binomial mean + gamma-Poisson mean
            E[:, g] = original*np.minimum(F, 1.0) + np.maximum(intensity, 0)
        if np.any(intensity > 0):
            dispersion = float(line.overdispersion[g])
            if dispersion > 1e-8:
                intensity = grng.gamma(1/dispersion, intensity*dispersion)
            values += grng.poisson(intensity)
        B[:, g] = values
    pin = {}
    if E is not None:
        from vcc_exact_totals import pin_totals
        expected_totals = E.sum(0); realized = B.sum(0); big = np.rint(expected_totals) > 50
        pin = dict(bias_ratio_median=float(np.median(realized[big]/expected_totals[big])) if big.any() else float('nan'),
                   genes_T_gt_50=int(big.sum()), pre_pin_abs_err=float(np.abs(realized-np.rint(expected_totals)).sum()))
        B = pin_totals(B, expected_totals, rng); pin['_expected_totals'] = expected_totals   # float64: float32 rounds x.5 totals 1 off
        pin['post_pin_max_abs_err'] = float(np.abs(B.sum(0)-np.maximum(0, np.rint(expected_totals))).max())
    if B.max(initial=0) > np.iinfo(np.int32).max: raise ValueError("count exceeds int32 range")
    out = sp.csr_matrix(B.astype(np.int32)); out.eliminate_zeros(); validate_counts(out)
    planned = np.log2((projected*1e6+0.5)/(line.share*1e6+0.5))
    return out, dict(net_share_change=mass-1, projection_normalizer=mass, n_changed=int(len(cols)),
                     effective_pi=effective_pi, library_scale=lib_scale, **pin), planned


def target_metadata(path, line, fallback_n, enabled=True):
    """Optional CSV metadata; supports shared or per-context rows/column suffixes."""
    df = pd.read_csv(path)
    if 'target_gene' not in df: raise ValueError("pert_counts.csv needs target_gene")
    if 'context' in df: df = df[df.context.astype(str) == str(line)]
    if df.target_gene.duplicated().any(): raise ValueError("duplicate target rows; supply context labels")
    result = {}
    for row in df.to_dict('records'):
        n, umi = fallback_n, None
        if enabled:
            for field, default in [('n_cells', fallback_n), ('median_umi_per_cell', None)]:
                value = row.get(field+'_'+str(line), row.get(field, default))
                if value is None or pd.isna(value): value = default
                if value is not None and (not np.isfinite(value) or value <= 0): raise ValueError(f"invalid {field} for {row['target_gene']}")
                if field == 'n_cells':
                    if float(value) != int(value): raise ValueError("n_cells must be integer")
                    n = int(value)
                else: umi = None if value is None else float(value)
        result[str(row['target_gene'])] = dict(n_cells=n, median_umi_per_cell=umi)
    if not result: raise ValueError(f"no requested targets for context {line}")
    return result


def check_output(path, expected_shape, block=512):
    a = ad.read_h5ad(path, backed='r')
    try:
        if a.shape != expected_shape: raise ValueError(f"output shape {a.shape} != {expected_shape}")
        stored = 0
        for start in range(0, a.n_obs, block):
            matrix = a.X[start:start+block]; validate_counts(matrix, f"output rows {start}:{start+block}")
            if np.any(np.asarray(matrix.sum(1)).ravel() > 1e6): raise ValueError('cell exceeds 1,000,000 count cap')
            stored += matrix.nnz if sp.issparse(matrix) else matrix.size
        if a.n_obs > 400000 or stored > 4750000000: raise ValueError('submission exceeds cell or stored-entry cap')
    finally: a.file.close()
    return True


def answer_free_check(cells, line, planned, j):
    pooled = np.asarray(cells.sum(0)).ravel(); pooled = pooled/max(pooled.sum(), 1)*1e6
    realized = np.log2((pooled+0.5)/(line.share*1e6+0.5))
    mask = np.isfinite(realized) & (line.share > 0)
    if j is not None: mask[j] = False
    return dict(projected_effect_mae=float(np.abs(realized[mask]-planned[mask]).mean()),
                n_cells=cells.shape[0], median_umi_per_cell=float(np.median(np.asarray(cells.sum(1)).ravel())))


def prediction_ontarget(path, targets):
    p = Path(path)
    if not p.is_dir() or not (p/'ontarget_lfc.npy').exists(): return np.full(len(targets), np.nan)
    names = np.load(p/'targets.npy').astype(str); values = np.load(p/'ontarget_lfc.npy'); pos = dict(zip(names, values))
    return np.asarray([pos.get(t, np.nan) for t in targets])


def cmd_submit(a):
    cfg = Cfg(a)
    if a.settings:
        settings = json.loads(Path(a.settings).read_text())
        if settings.get('schema') != SCHEMA: raise ValueError('use generator settings from this notebook')
        chosen = settings['settings']
        for key in ['alpha', 'z', 'min_effect', 'call', 'top_n', 'sign_rule', 'detect_floor', 'z_req', 'gate_from', 'min_bg',
                    'pi', 'up_w', 'n_cells', 'metadata', 'match_pool', 'tilt']:
            if key in chosen: setattr(cfg, key, chosen[key])
        if a.delta_array == 'selected': a.delta_array = settings['delta_array']
    if a.delta_array == 'selected': raise ValueError('--delta-array selected needs --settings')
    data = Path(a.data); genes = read_list(data/'gene_names.csv', 'gene_name')
    gp = {g: i for i, g in enumerate(genes)}; deltas = dict(x.split('=', 1) for x in a.delta)
    specifications = {line: target_metadata(data/'pert_counts.csv', line, cfg.n_cells, cfg.metadata) for line in a.lines}
    sets = [set(s) for s in specifications.values()]
    if any(s != sets[0] for s in sets): raise ValueError('all submitted contexts must share the exact perturbation panel')
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True); path = out/'prediction.h5ad'
    if path.exists(): raise FileExistsError('choose a new output directory to preserve the previous submission')
    for line_specs in specifications.values():
        if any(s['n_cells'] != a.required_n_cells for s in line_specs.values()):
            raise ValueError(f'2026 submissions require exactly {a.required_n_cells} cells per target; optional CSV metadata must respect that contract')
    frames = [pd.DataFrame(dict(target_gene=t, context=line), index=[f'{line}_{t}_{i:05d}' for i in range(spec['n_cells'])])
              for line in a.lines for t, spec in specifications[line].items()]
    obs = pd.concat(frames); obs['target_gene'] = obs.target_gene.astype('category'); obs['context'] = obs.context.astype('category')
    report = dict(settings=cfg.as_dict(), delta_array=a.delta_array, lines={})
    with h5py.File(path, 'w') as f:
        write_elem(f, 'obs', obs); write_elem(f, 'var', pd.DataFrame(index=pd.Index(genes, name='gene_name')))
        f.attrs['encoding-type'] = 'anndata'; f.attrs['encoding-version'] = '0.1.0'; ds = None
        for name in a.lines:
            if name not in deltas: raise ValueError(f'no predictions for {name}')
            line = Line.from_h5ad(data/f'context_{name}.h5ad', genes)
            targets = list(specifications[name]); D, info = load_effects(deltas[name], targets, genes, a.delta_array, cfg)
            D = calibrate(D, line, targets, cfg, info=info); del info
            own = prediction_ontarget(deltas[name], targets)
            if a.kd_orion or a.kd_tcell:
                legacy_own = kd_from_sources(targets, genes, a.kd_orion or (), a.kd_tcell, a.kd_tcell_condition)
                own = np.where(np.isfinite(own), own, legacy_own)
            kd = [None] * len(targets) if is_nochange(a.delta_array) else kd_residual(own, cfg); checks = []; diags = []
            # Pick max count once; different n_cells use prefixes of this same pool.
            pool_cfg = copy.copy(cfg); pool_cfg.n_cells = max(s['n_cells'] for s in specifications[name].values())
            common = pick_cells(line, np.zeros(len(genes)), pool_cfg, np.random.default_rng(seed_for(cfg.seed, name, 'base_controls')))
            global _SUBMIT_CTX
            _SUBMIT_CTX = (line, D, gp, kd, cfg, specifications[name], name, common, targets)
            if getattr(a, 'workers', 1) > 1:
                import multiprocessing as mp
                pool = mp.get_context('fork').Pool(a.workers); results = pool.imap(_generate_task, range(len(targets)), chunksize=2)
            else:
                pool = None; results = map(_generate_task, range(len(targets)))
            expected = []
            for i, (t, (cells, diag, planned)) in enumerate(zip(targets, results)):
                if '_expected_totals' in diag: expected.append(diag.pop('_expected_totals'))
                diags.append(dict(target=t, **diag))
                if i < a.checks: checks.append(dict(target=t, **answer_free_check(cells, line, planned, gp.get(t))))
                # float32 holds these integer counts exactly and matches vcc prep's output encoding,
                # so prep packages without a second full-matrix copy (int32 needs ~37 GB at 360k cells).
                if cells.data.max(initial=0) >= 2**24: raise ValueError('count too large for exact float32')
                cells = cells.astype(np.float32)
                if ds is None: write_elem(f, 'X', int64_indptr(cells)); ds = sparse_dataset(f['X'])
                else: ds.append(cells)
            if pool is not None: pool.close(); pool.join()
            if expected: np.savez(out/f'expected_totals_{name}.npz', targets=np.asarray(targets, dtype=str), expected=np.stack(expected))
            pd.DataFrame(checks).to_csv(out/f'checks_{name}.csv', index=False)
            pd.DataFrame(diags).to_csv(out/f'projection_{name}.csv', index=False)
            report['lines'][name] = dict(n_targets=len(targets), n_cells=sum(s['n_cells'] for s in specifications[name].values()),
                metadata=specifications[name], measured_ontarget=int(np.isfinite(own).sum()))
            print(f"generated context {name}", flush=True)
    report['integer_counts'] = check_output(path, (len(obs), len(genes)))
    save_json(out/'generator_report.json', report); print(f'WROTE {path}; checked every output row')


_SUBMIT_CTX = None


def _generate_task(i):
    """One knockdown of cmd_submit's current line (module level so a forked worker pool can run it)."""
    line, D, gp, kd, cfg, spec, name, common, targets = _SUBMIT_CTX
    t = targets[i]
    local = copy.copy(cfg); local.n_cells = spec[t]['n_cells']; local.library_size = spec[t]['median_umi_per_cell']
    pick = None if local.tilt else common[:local.n_cells]
    return generate_one(line, D[i], gp.get(t), kd[i], local, pick=pick, line_key=name, target_key=t)


class StreamingH5AD:
    """An .h5ad whose obs is known up front and whose X arrives one block at a time.

    Nothing but the current block is held in memory. The row pointer is int64 from the first
    block, so the finished file may exceed 2**31 stored values. Every block is cast to the first
    block's dtype, and the row total is checked against obs on close.
    """
    def __init__(self, path, obs, genes):
        self.path = Path(path); self.path.parent.mkdir(parents=True, exist_ok=True)
        self.n = len(obs); self.rows = 0; self.ds = None; self.dtype = None
        self.f = h5py.File(self.path, 'w')
        self.f.attrs['encoding-type'] = 'anndata'; self.f.attrs['encoding-version'] = '0.1.0'
        write_elem(self.f, 'obs', obs); write_elem(self.f, 'var', pd.DataFrame(index=pd.Index(genes, name='gene_name')))

    def append(self, block):
        block = block.tocsr() if sp.issparse(block) else sp.csr_matrix(np.asarray(block))
        if self.ds is None:
            self.dtype = block.dtype
            write_elem(self.f, 'X', int64_indptr(block)); self.ds = sparse_dataset(self.f['X'])
        else:
            self.ds.append(block.astype(self.dtype, copy=False))
        self.rows += block.shape[0]

    def close(self):
        self.f.close()
        if self.rows != self.n: raise ValueError(f'{self.path}: wrote {self.rows} rows for {self.n} obs entries')


def labelled_obs(labels, prefix):
    labels = np.asarray(labels, dtype=str)
    return pd.DataFrame({'target_gene': pd.Categorical(labels)}, index=[f'{prefix}{i:07d}' for i in range(len(labels))])


def h1_controls(h, labels, genes, seed):
    """The benchmark's seeded control split: the generator half (a Line) and the scorer's half (rows). The gate and
    the detection floor's reference size come from the scorer's half, which is what the evaluator compares against."""
    ci = np.flatnonzero(labels == CONTROL)
    if len(ci) < 2: raise ValueError('H1 needs controls')
    ci = np.random.default_rng(seed_for(seed, 'H1', 'controls')).permutation(ci)
    train_rows, score_rows = ci[:len(ci) // 2], ci[len(ci) // 2:]
    small = h[np.sort(train_rows)].to_memory(); validate_counts(small.X, 'H1 raw counts')
    line = Line(aligned_X(small, genes), None, genes); del small
    # The scorer uses its REAL control half. Match that exact gate, even when
    # the generator half happens to put a near-threshold gene on the other side.
    line.gate = expression_gate(moments(h, np.sort(score_rows), alignment(h.var_names.astype(str), genes))['mean_cpm'])
    line.n_reference = len(score_rows)
    return line, train_rows, score_rows


def reference_targets(real_path):
    """Targets and their real cell counts in an existing scoring reference (a previous stage's real.h5ad)."""
    a = ad.read_h5ad(real_path, backed='r')
    try: counts = a.obs['target_gene'].astype(str).value_counts()
    finally: a.file.close()
    counts = counts.drop(CONTROL)
    return sorted(counts.index), {t: int(n) for t, n in counts.items()}


def stream_h1(a, genes, real_path, line_pseudo=0.5, block=4096):
    """Write the scoring reference straight from the backed H1 file, one target at a time.

    Controls are split into two disjoint, seeded halves: the first half becomes the generator's
    input (held in memory, it is sampled for every target); the second half is appended to the
    reference after the perturbed cells. In oracle mode each target's cells are split the same
    way as before: the estimation half only produces that target's effect vector, the other half
    is appended to the reference. Peak memory is the control half plus one target.

    Returns (Line, targets, {target: reference cell count}, oracle effects or None, oracle own LFCs).
    """
    h = ad.read_h5ad(a.cells, backed='r')
    try:
        labels = h.obs[a.obs_col].astype(str).values
        requested = sorted(set(read_list(a.targets, 'target_gene'))) if a.targets else sorted(set(labels) - {CONTROL})
        counts = pd.Series(labels).value_counts()
        targets = [t for t in requested if counts.get(t, 0) >= a.min_cells and t != CONTROL]
        if a.n_targets and len(targets) > a.n_targets:
            targets = sorted(np.random.default_rng(seed_for(a.seed, 'H1', 'panel')).choice(targets, a.n_targets, replace=False).tolist())
        if not targets: raise ValueError('no eligible H1 targets')

        def subset(rows):
            rows = np.sort(np.asarray(rows))
            small = h[rows].to_memory(); validate_counts(small.X, 'H1 raw counts')
            return aligned_X(small, genes)

        line, train_rows, score_rows = h1_controls(h, labels, genes, a.seed)
        save_json(Path(real_path).parent / 'control_split.json', dict(seed=a.seed,
                  generator_controls=len(train_rows), scoring_controls=len(score_rows),
                  gate_source='scoring real controls', genes_above_5_cpm=int(line.gate.sum())))
        plan = {}
        for t in targets:
            rows = np.flatnonzero(labels == t)
            if a.mode == 'oracle':
                rows = np.random.default_rng(seed_for(a.seed, 'H1', t, 'split')).permutation(rows)
                plan[t] = (rows[:len(rows) // 2], rows[len(rows) // 2:])
            else:
                plan[t] = (None, rows)
        real_counts = {t: len(ref) for t, (_, ref) in plan.items()}
        obs_labels = [t for t in targets for _ in range(real_counts[t])] + [CONTROL] * len(score_rows)
        writer = StreamingH5AD(real_path, labelled_obs(obs_labels, 'real_'), genes)
        oracle, own = ([], []) if a.mode == 'oracle' else (None, None)
        gp = {g: i for i, g in enumerate(genes)}
        try:
            for i, t in enumerate(targets):
                estimation, reference = plan[t]
                if estimation is not None:
                    x = subset(estimation); lib = np.asarray(x.sum(1)).ravel()
                    mean = np.asarray((sp.diags(1e6 / np.maximum(lib, 1)) @ x).mean(0)).ravel()
                    d = np.log2((mean + line_pseudo) / (line.cpm + line_pseudo)); oracle.append(d)
                    own.append(d[gp[t]] if t in gp else np.nan); del x
                writer.append(subset(reference))
                print(f'  reference: {i + 1}/{len(targets)} targets', end='\r', flush=True)
            sorted_score = np.sort(score_rows)
            for start in range(0, len(sorted_score), block):
                writer.append(subset(sorted_score[start:start + block]))
        finally:
            writer.close()
        print(f'\n  wrote {real_path}: {writer.rows:,} cells ({len(score_rows):,} scoring controls)', flush=True)
        return (line, targets, real_counts, None if oracle is None else np.asarray(oracle),
                None if own is None else np.asarray(own))
    finally:
        h.file.close()


def write_arm(path, line, effects, kd, targets, n_cells, common, cfg, line_key='H1', block=4096):
    """Generate each target's cells and append them immediately, then the generator's control cells.

    cell-eval2 requires the predicted file to carry the control label even though vcc2026 scores
    every perturbation against the REAL controls (`control_source: real`), so the generator's own
    input controls are appended as before; they are not used as the scoring reference.
    """
    gp = {g: i for i, g in enumerate(line.genes)}
    labels = [t for t in targets for _ in range(n_cells[t])] + [CONTROL] * line.n
    writer = StreamingH5AD(path, labelled_obs(labels, 'pred_'), line.genes)
    try:
        for i, t in enumerate(targets):
            local = copy.copy(cfg); local.n_cells = n_cells[t]
            cells = generate_one(line, effects[i], gp.get(t), kd[i], local, pick=common[:local.n_cells],
                                 line_key=line_key, target_key=t)[0]
            writer.append(cells)
        for start in range(0, line.n, block):
            writer.append(line.X[start:start + block].astype(np.int32))
    finally:
        writer.close()


class OfficialEvaluator:
    """Runs official cell-eval2 against one reference file, and scores arms with the leaderboard formula.

    The reference never changes between arms, so its expensive side (real DE tests, pseudobulks, noise corrections)
    is cached once with `prep-cache --side real` and reused by every arm through `--cache-real`. The ceiling and the
    competition bundle (cell-eval2's generic-response baseline as the 0 end, a real-vs-real replicate anchor as the
    1 end) depend on the real data only, so they are built once too, and a later stage can point at them
    (`reference=` a previous benchmark folder) instead of rebuilding them.
    """
    def __init__(self, a, reference, out, reuse=None):
        self.binary = shutil.which(a.cell_eval)
        if self.binary is None: raise RuntimeError('install official cell-eval2 (pip install -r requirements-eval.txt) before running the H1 benchmark')
        help_text = subprocess.run([self.binary, 'run', '--help'], check=True, text=True, capture_output=True).stdout
        if '--ceiling' not in help_text or '--cache-real' not in help_text or '--write-degenes' not in help_text:
            raise RuntimeError('installed cell-eval2 lacks --ceiling/--cache-real/--write-degenes; install requirements-eval.txt')
        self.a, self.reference, self.out = a, Path(reference), Path(out)
        home = Path(reuse) if reuse else self.out
        self.cache, self.bundle = home / 'real_cache', home / 'bundle'
        self.common = ['--pert-col', 'target_gene', '--control', CONTROL, '--profile', a.eval_profile,
                       '--input-type', 'counts', '--preset', 'vcc2026', '--set', 'de.backend=pdex',
                       '--set', f'num_threads={a.num_threads}']
        self.invocations = []

    def _run(self, command):
        print('  $ ' + ' '.join(shlex.quote(str(x)) for x in command), flush=True)
        subprocess.run([str(x) for x in command], check=True); self.invocations.append([str(x) for x in command])

    def prepare(self):
        # One-sided cache warm: the expression comparator must be named. vcc2026 on counts resolves
        # to bulk_lognorm (cell_eval2.competition.competition_comparator()).
        self._run([self.binary, 'prep-cache', '--side', 'real', '--adata', self.reference,
                   '--cache-real', self.cache, '--comparator', 'bulk_lognorm', *self.common])
        self._run([self.binary, 'run', '-ar', self.reference, '-o', self.out / 'ceiling', '--ceiling',
                   '--cache-real', self.cache, *self.common])

    def build_bundle(self):
        """The leaderboard's two ends for this reference: generic baseline (0) and replicate anchor (1)."""
        base = self.out / 'baseline'
        self._run([self.binary, 'baseline', '-ar', self.reference, '-o', base, '--save-pred', base / 'baseline_pred.h5ad',
                   '--cache-real', self.cache, *self.common])
        self._run([self.binary, 'prep-real-bundle', '--real', self.reference, '--baseline', base / 'baseline_pred.h5ad',
                   '-o', self.bundle, '--anchor-splits', str(self.a.anchor_splits), '--cache-real', self.cache, *self.common])
        (base / 'baseline_pred.h5ad').unlink()                  # large; the bundle keeps only aggregates

    def score(self, prediction, output):
        self._run([self.binary, 'run', '-ap', prediction, '-ar', self.reference, '-o', output,
                   '--cache-real', self.cache, '--write-degenes', *self.common])
        self.record(Path(output) / 'invocation.json')

    def scale(self, output):
        """The leaderboard number for one scored arm: cell-eval2 score against the bundle (from_replicate)."""
        output = Path(output)
        self._run([self.binary, 'score', '--user-agg', output / 'agg_results.csv', '--real-bundle', self.bundle,
                   '--user-meta', output / 'run_meta.json', '-o', output / 'scaled.csv'])
        return read_scaled(output / 'scaled.csv')

    def record(self, path):
        import importlib.metadata
        versions = {}
        for package in ['cell-eval2', 'pdex', 'numpy', 'scipy', 'anndata']:
            try: versions[package] = importlib.metadata.version(package)
            except importlib.metadata.PackageNotFoundError: versions[package] = 'not installed in generator environment'
        save_json(path, dict(commands=self.invocations, versions=versions))


def read_scaled(path):
    """A `cell-eval2 score` CSV -> {metric: leaderboard value}, with 'overall' = avg_score (from_replicate)."""
    frame = pd.read_csv(path).set_index('metric')
    values = frame['from_replicate'].astype(float).to_dict()
    values['overall'] = values.pop('avg_score')
    return values


def per_target_table(arm_dir, alpha=0.05):
    """Per knockdown: our called and the real DE genes (scorer's p_adj < alpha, own gene excluded), their overlap,
    precision against chance (random genes hit real DE genes at n_real / scored genes), direction accuracy, recall."""
    pred = pd.read_parquet(Path(arm_dir) / 'de_pred.parquet'); real = pd.read_parquet(Path(arm_dir) / 'de_real.parquet')
    rows = []
    for t, r in real.groupby('target'):
        p = pred[pred.target == t]
        P = set(p.feature[(p.p_adj < alpha) & (p.feature != t)]); R = set(r.feature[(r.p_adj < alpha) & (r.feature != t)])
        both = P & R; same = 0
        if both:
            lp = p.set_index('feature').log2_fold_change; lr = r.set_index('feature').log2_fold_change
            same = int(sum(np.sign(lp[x]) == np.sign(lr[x]) for x in both))
        scored = r.feature.nunique(); precision = len(both) / len(P) if P else np.nan
        rows.append(dict(target=t, n_called=len(P), n_real=len(R), overlap=len(both), same_direction=same,
                         precision=precision, chance_precision=len(R) / scored if scored else np.nan,
                         enrichment=precision * scored / len(R) if R else np.nan,
                         direction_accuracy=same / len(both) if both else np.nan, recall=len(both) / len(R) if R else np.nan))
    return pd.DataFrame(rows).set_index('target')


def per_target_summary(table):
    return dict(median_n_called=float(table.n_called.median()), median_n_real=float(table.n_real.median()),
                median_precision=float(table.precision.median()), median_chance_precision=float(table.chance_precision.median()),
                median_enrichment=float(table.enrichment.median()), direction_accuracy=float(table.same_direction.sum() / max(table.overlap.sum(), 1)),
                median_recall=float(table.recall.median()))


def metric_rank_table(outputs, directions_path):
    """Rank only explicitly requested metric families, with explicit directions.

    The evaluator emits diagnostics as well as scoring metrics. Never average
    their ranks blindly. Supply metric,direction ('min'/'max') from current rules.
    """
    directions = pd.read_csv(directions_path)
    if not {'metric', 'direction'} <= set(directions) or directions.metric.duplicated().any():
        raise ValueError('directions CSV needs unique metric,direction columns')
    values = []
    for label, path in outputs.items():
        df = pd.read_csv(Path(path)/'agg_results.csv')
        if not {'metric', 'value'} <= set(df):
            if {'metric', 'score'} <= set(df): df = df.rename(columns={'score': 'value'})
            elif set(directions.metric) <= set(df.columns):
                # cell-eval2 writes a wide describe-style frame; use only its mean row.
                mean_row = np.zeros(len(df), bool)
                for col in df.columns:
                    if col not in set(directions.metric): mean_row |= df[col].astype(str).eq('mean').to_numpy()
                if not mean_row.any() and len(df) == 1: mean_row[:] = True
                if mean_row.sum() != 1: raise ValueError('official aggregate CSV needs exactly one mean row')
                df = pd.DataFrame({'metric': directions.metric, 'value': [df.loc[mean_row, m].iloc[0] for m in directions.metric]})
            else: raise ValueError(f'unrecognized official aggregate schema: {df.columns.tolist()}')
        for metric, direction in directions[['metric', 'direction']].itertuples(index=False, name=None):
            v = df.loc[df.metric == metric, 'value']
            if len(v) != 1 or not np.isfinite(v.iloc[0]): raise ValueError(f'{label}: metric {metric} missing/nonfinite')
            if direction not in ['min', 'max']: raise ValueError('metric direction must be min/max')
            values.append(dict(variant=label, metric=metric, value=float(v.iloc[0]), direction=direction))
    frame = pd.DataFrame(values); frame['rank'] = np.nan
    for metric in directions.metric:
        use = frame.metric == metric
        frame.loc[use, 'rank'] = frame.loc[use, 'value'].rank(ascending=frame.loc[use, 'direction'].iloc[0] == 'min', method='average')
    return frame, frame.groupby('variant')['rank'].mean().sort_values(kind='stable')


ARM_KEYS = ['method', 'call', 'top_n', 'sign_rule', 'min_bg', 'detect_floor', 'z_req', 'gate_from', 'alpha', 'z', 'min_effect', 'pi', 'up_w', 'cell_rule']


def benchmark_variants(a):
    """The arms to score: an explicit JSON list (--variants) or the product of the grid flags."""
    if a.variants:
        variants = json.loads(Path(a.variants).read_text())
        for v in variants:
            unknown = set(v) - set(ARM_KEYS) - {'label'}
            if unknown or 'method' not in v: raise ValueError(f'bad variant {v}: needs method, unknown keys {sorted(unknown)}')
            v.setdefault('cell_rule', 'fixed')
        return variants
    return [dict(method=m, call=a.call, top_n=a.top_n, sign_rule=a.sign_rule, detect_floor=a.detect_floor, z_req=a.z_req,
                 alpha=alpha, z=z, min_effect=me, pi=pi, up_w=w, cell_rule=rule)
            for m in a.delta_array for z in a.z_grid for me in a.min_effect_grid for alpha in a.alpha_grid
            for pi in a.pi_grid for w in a.up_w_grid for rule in a.cell_rules]


def cmd_benchmark(a):
    """Official H1 benchmark, streamed: reference and each arm go to disk one target at a time.

    Every arm is scored raw (agg_results.csv), with the leaderboard formula (scaled.csv) and per knockdown
    (per_target.csv). Selection uses the scaled overall. Oracle mode never writes selected generator settings.
    """
    genes = read_list(a.genes, 'gene_name')
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    if a.mode == 'transfer':
        provenance = json.loads((Path(a.predictions)/'prediction_report.json').read_text())
        if provenance.get('heldout_background') != 'H1': raise ValueError('generator tuning requires the H1 validation: predictions with H1 held out')
        if provenance.get('schema') != SCHEMA: raise ValueError('rebuild the validation predictions with vcc3_ridge.py')
        if 'H1' in provenance.get('source_backgrounds', []): raise ValueError('H1 leaked into the validation sources')
    if a.mode == 'transfer':
        variants = benchmark_variants(a)
    elif a.variants:
        # Oracle arms on the same real cells: the oracle effects and/or true no change, at any alpha and cell rule.
        variants = benchmark_variants(a)
        if {v['method'] for v in variants} - {'oracle', 'nochange'}: raise ValueError('oracle arms use method oracle or nochange')
    else:
        variants = [dict(method='oracle', alpha=a.alpha, pi=a.pi, up_w=a.up_w, cell_rule='matched')]
    if a.mode == 'oracle' and (a.reference or a.resume):
        raise ValueError('the oracle estimates effects from its own split of the real cells: no --reference or --resume')
    if a.mode == 'transfer':
        unknown = {v['method'] for v in variants} - set(provenance['models']) - {'nochange'}
        if unknown: raise ValueError(f'unknown models {sorted(unknown)}; the exam has {provenance["models"]}')
        # Fail before the (hours-long) reference build: an SE gate needs an SE table, its own or the gating array's.
        for v in variants:
            se_from = v.get('gate_from') or v['method']
            if float(v.get('z', 0) or 0) > 0 and v.get('call', 'z') == 'z' and not (Path(a.predictions) / f'{se_from}_se.npy').exists():
                raise ValueError(f"{v.get('label', v['method'])}: z* > 0 needs {se_from}_se.npy; ridge outputs have none, so use z = 0 "
                                 "or gate_from='re'")
    if a.reference:
        # A previous stage's reference: same real cells, cache and bundle; rebuild only the generator's controls.
        ref = Path(a.reference); real_path = ref / 'real.h5ad'
        if not (ref / 'bundle').exists(): raise ValueError(f'{ref} has no competition bundle')
        targets, real_counts = reference_targets(real_path)
        h = ad.read_h5ad(a.cells, backed='r')
        try: line = h1_controls(h, h.obs[a.obs_col].astype(str).values, genes, a.seed)[0]
        finally: h.file.close()
        oracle_effects = oracle_own = None
        evaluator = OfficialEvaluator(a, real_path, out, reuse=ref)
        save_json(out / 'reference.json', dict(reused=str(ref.resolve()), targets=len(targets)))
    elif a.resume and all((out / x).exists() for x in ('real.h5ad', 'real_cache/manifest.json', 'ceiling/ceiling_agg.csv', 'control_split.json')):
        # An interrupted run of this folder: keep its real cells, cache and ceiling; build only what is missing.
        real_path = out / 'real.h5ad'
        if json.loads((out / 'control_split.json').read_text()).get('seed') != a.seed: raise ValueError('--resume with a different --seed')
        targets, real_counts = reference_targets(real_path)
        h = ad.read_h5ad(a.cells, backed='r')
        try: line = h1_controls(h, h.obs[a.obs_col].astype(str).values, genes, a.seed)[0]
        finally: h.file.close()
        oracle_effects = oracle_own = None
        evaluator = OfficialEvaluator(a, real_path, out)
        if a.mode == 'transfer' and not (out / 'bundle/manifest.json').exists():
            for partial in (out / 'baseline', out / 'bundle'):           # left by the interrupted run
                if partial.exists(): shutil.rmtree(partial)
            evaluator.build_bundle()
        print(f'resumed {out}: {len(targets)} targets, reference reused', flush=True)
    else:
        real_path = out/'real.h5ad'
        line, targets, real_counts, oracle_effects, oracle_own = stream_h1(a, genes, real_path)
        evaluator = OfficialEvaluator(a, real_path, out); evaluator.prepare()
        if a.mode == 'transfer' or a.scale: evaluator.build_bundle()
    print(f'{len(variants)} arms to score on {len(targets)} targets', flush=True)
    if a.mode == 'transfer':
        kd = kd_residual(prediction_ontarget(a.predictions, targets), Cfg(a))
        effects = {}
        for m, gate in sorted({(v['method'], v.get('gate_from', '') or '') for v in variants}):
            effects[m, gate] = ((np.zeros((len(targets), len(genes)), np.float32), {}) if is_nochange(m) else
                                (load_deltas(a.predictions, targets, genes, m), load_info(a.predictions, targets, genes, m, gate)))
    else:
        # Ceiling-only oracle inputs: per-cell arithmetic CPM from each target's estimation half.
        effects = {('oracle', ''): (oracle_effects, {}), ('nochange', ''): (np.zeros_like(oracle_effects), {})}
        kd = kd_residual(oracle_own, Cfg(a))
    outputs, settings, scaled, summaries = {}, {}, {}, {}
    maximum = max([a.n_cells] + list(real_counts.values()))
    pool_cfg = Cfg(a); pool_cfg.n_cells = maximum
    common = pick_cells(line, np.zeros(len(genes)), pool_cfg, np.random.default_rng(seed_for(a.seed, 'H1', 'base_controls')))
    for index, variant in enumerate(variants):
        method = variant['method']; effects_array, info = effects[method, variant.get('gate_from', '') or '']
        label = variant.get('label') or f'{method}_{index:03d}'
        cfg = Cfg(argparse.Namespace(**dict(vars(a), **{k: v for k, v in variant.items() if k not in ('method', 'cell_rule', 'label')})))
        calibrated = calibrate(effects_array, line, targets, cfg, info=info)
        n_cells = {t: (real_counts[t] if variant['cell_rule'] == 'matched' else a.n_cells) for t in targets}
        pred_path = out/f'{label}.h5ad'
        write_arm(pred_path, line, calibrated, [None] * len(targets) if is_nochange(method) else kd, targets, n_cells, common, cfg)
        del calibrated
        evaluator.score(pred_path, out/label)
        if not a.keep_arm_files: pred_path.unlink()   # ~2 GB each; the scores are in out/label
        outputs[label] = str(out/label); settings[label] = dict(variant)
        save_json(out/f'{label}.json', settings[label])
        table = per_target_table(out/label); table.to_csv(out/label/'per_target.csv')
        summaries[label] = per_target_summary(table)
        if a.mode == 'transfer' or a.scale: scaled[label] = evaluator.scale(out/label)
        print(f'official evaluation complete: {label}' + (f" (scaled overall {scaled[label]['overall']:.4f})" if label in scaled else ''), flush=True)
    pd.DataFrame(summaries).T.to_csv(out/'per_target_summary.csv')
    save_json(out/'benchmark.json', dict(mode=a.mode, outputs=outputs, settings=settings, n_targets=len(targets),
              ceiling=str(evaluator.out / 'ceiling') if not a.reference else str(Path(a.reference) / 'ceiling'),
              real_cache=str(evaluator.cache), bundle=str(evaluator.bundle),
              note='H1 is generator calibration; final evidence requires unseen lines. Oracle is a generator-only ceiling check.'))
    if a.metric_directions or a.eval_profile == 'vcc2026':
        directions = a.metric_directions or str(out/'metric_directions.csv')
        if a.metric_directions is None:
            pd.DataFrame({'metric': list(VCC2026_METRICS), 'direction': list(VCC2026_METRICS.values())}).to_csv(directions, index=False)
        frame, ranks = metric_rank_table(outputs, directions)
        frame.to_csv(out/'metric_ranks.csv', index=False); ranks.to_csv(out/'aggregate_rank.csv')
    if a.mode == 'oracle' and scaled:
        table = pd.DataFrame(scaled).T.sort_values('overall', ascending=False)
        table.index.name = 'arm'; table.to_csv(out/'scaled_scores.csv')
    if a.mode == 'transfer':
        table = pd.DataFrame(scaled).T.sort_values('overall', ascending=False)
        table.index.name = 'arm'; table.to_csv(out/'scaled_scores.csv')
        eligible = [label for label in table.index if settings[label]['cell_rule'] == 'fixed']
        if not eligible: raise ValueError('include fixed cell rule to select submission-compatible generator settings')
        best = settings[eligible[0]]
        chosen = {k: v for k, v in best.items() if k not in ('method', 'cell_rule', 'label')}
        margin = {}
        if a.margin_over:
            if a.margin_over not in table.index: raise ValueError(f'--margin-over {a.margin_over!r} is not an arm here')
            margin = dict(margin_over=a.margin_over, margin_over_overall=float(table.loc[a.margin_over, 'overall']),
                          margin=float(table.loc[eligible[0], 'overall'] - table.loc[a.margin_over, 'overall']))
        save_json(out/'selected_generator.json', dict(settings=dict(Cfg(argparse.Namespace(**dict(vars(a), **chosen))).as_dict() | chosen, n_cells=a.n_cells),
                   cell_rule='fixed', n_cells=a.n_cells, delta_array=best['method'], arm=eligible[0],
                   scaled_overall=float(table.loc[eligible[0], 'overall']), **margin,
                   selection='scaled overall (leaderboard formula)', schema=SCHEMA,
                   note='H1 is calibration, not an untouched test'))


def add_switches(p):
    p.add_argument('--alpha', type=float, default=1.0); p.add_argument('--n-cells', type=int, default=400)
    p.add_argument('--pi', type=float, default=0); p.add_argument('--up-w', type=float, default=0.75)
    p.add_argument('--metadata', action=argparse.BooleanOptionalAction, default=True)
    p.add_argument('--z', type=float, default=0); p.add_argument('--min-effect', type=float, default=0)
    p.add_argument('--call', choices=['z', 'topn'], default='z'); p.add_argument('--top-n', type=int, default=1000)
    p.add_argument('--sign-rule', action='store_true'); p.add_argument('--detect-floor', action='store_true')
    p.add_argument('--z-req', type=float, default=4.0); p.add_argument('--gate-from', default='')
    p.add_argument('--min-bg', type=int, default=1)
    p.add_argument('--match-pool', action='store_true'); p.add_argument('--tilt', action='store_true')
    p.add_argument('--seed', type=int, default=2026)
    p.add_argument('--exact-totals', action='store_true', help='pin each gene\'s pseudobulk total per knockdown to its expected value')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('submit'); add_switches(s)
    for key in ['data', 'out']: s.add_argument('--'+key, required=True)
    s.add_argument('--lines', nargs='+', default=['A', 'B', 'C']); s.add_argument('--delta', action='append', required=True)
    s.add_argument('--required-n-cells', type=int, default=400)
    s.add_argument('--delta-array', default='selected'); s.add_argument('--settings')
    s.add_argument('--kd-orion', nargs='*'); s.add_argument('--kd-tcell'); s.add_argument('--kd-tcell-condition', default='Stim48hr')
    s.add_argument('--checks', type=int, default=20)
    s.add_argument('--workers', type=int, default=1, help='processes generating knockdowns in parallel (output identical)')
    for command, mode in [('benchmark', 'transfer'), ('oracle', 'oracle')]:
        b = sub.add_parser(command); add_switches(b); b.set_defaults(mode=mode)
        for key in ['cells', 'genes', 'out']: b.add_argument('--'+key, required=True)
        b.add_argument('--predictions', required=mode == 'transfer'); b.add_argument('--targets')
        b.add_argument('--obs-col', default='target_gene'); b.add_argument('--delta-array', nargs='+', default=['re'])
        b.add_argument('--min-cells', type=int, default=100); b.add_argument('--n-targets', type=int, default=0)
        b.add_argument('--alpha-grid', nargs='+', type=float, default=[0.75, 1.0, 1.5, 2.0])
        b.add_argument('--z-grid', nargs='+', type=float, default=[2.0, 3.0, 4.0])
        b.add_argument('--min-effect-grid', nargs='+', type=float, default=[0.1, 0.2])
        b.add_argument('--pi-grid', nargs='+', type=float, default=[0.0])
        b.add_argument('--up-w-grid', nargs='+', type=float, default=[0.75])
        b.add_argument('--cell-rules', nargs='+', choices=['matched', 'fixed'], default=['fixed'])
        b.add_argument('--cell-eval', default='cell-eval2'); b.add_argument('--eval-profile', default='vcc2026')
        b.add_argument('--num-threads', type=int, default=4)
        b.add_argument('--keep-arm-files', action='store_true', help='keep each scored arm .h5ad (deleted by default)')
        b.add_argument('--variants', help='JSON list of arms (method, call, top_n, sign_rule, detect_floor, z_req, alpha, z, min_effect, ...)')
        b.add_argument('--reference', help='a previous benchmark folder: reuse its real.h5ad, real_cache, ceiling and bundle')
        b.add_argument('--anchor-splits', type=int, default=5, help='replicate-anchor splits for the competition bundle')
        b.add_argument('--margin-over', help="an arm label; selected_generator.json records the winner's scaled-overall margin over it")
        b.add_argument('--scale', action='store_true', help='oracle: also build the competition bundle on its real cells and scale every arm')
        b.add_argument('--resume', action='store_true', help="reuse this folder's real cells, cache, ceiling (and bundle, if complete) from an interrupted run")
        b.add_argument('--metric-directions', help='CSV metric,direction for current official scoring metrics only')
    a = p.parse_args()
    if a.cmd == 'submit': cmd_submit(a)
    else: cmd_benchmark(a)
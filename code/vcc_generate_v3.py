"""Generate counts with projected composition, metadata cell counts/library sizes,
state-dependent overdispersed induction and reproducible common control sampling.

`benchmark` tunes on LOBO-H1 transfer predictions with official cell-eval2 (vcc2026), streaming every
file one target at a time; the real side is cached once and the ceiling computed once.
`oracle` is a separate generator ceiling check and never selects submission knobs.
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
from vcc_common import int64_indptr, seed_for, validate_counts, shrink, save_json

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
CC_GENES = ("MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 HELLS RFC2 RPA2 NASP RAD51AP1 GMNN "
            "WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 "
            "CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 SMC4 "
            "CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 CDC20 TTK CDC25C KIF2C RANGAP1 "
            "NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA").split()


class Cfg:
    def __init__(self, a):
        self.alpha = getattr(a, "alpha", 1.0); self.max_lfc = 4.0; self.min_lfc = 0.0
        self.kd_clip = (0.01, 0.9); self.default_kd = 0.30; self.n_cells = getattr(a, "n_cells", 400)
        self.composition = "projection"; self.comp_thr = getattr(a, "comp_thr", 0.02)
        self.pi = getattr(a, "pi", 0.0); self.up_w = getattr(a, "up_w", 0.75)
        self.tilt = getattr(a, "tilt", False); self.match_pool = getattr(a, "match_pool", False) or self.tilt
        self.seed = getattr(a, "seed", 2026)
        self.metadata = getattr(a, "metadata", True)
        self.library_size = None
        if not 0 <= self.pi < 1 or not 0 <= self.up_w <= 1 or self.n_cells < 1 or self.alpha < 0:
            raise ValueError("invalid generator settings")

    def as_dict(self):
        return {k: v for k, v in vars(self).items()}


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
        self.cpm_cells = sp.diags(1e6 / np.maximum(self.lib, 1)) @ self.X
        self.cpm = np.asarray(self.cpm_cells.mean(0)).ravel()           # mean per-cell CPM (scorer gate uses > 5)
        self.share = np.asarray(self.X.sum(0)).ravel() / self.X.sum()   # pooled share of molecules
        self.gate = self.cpm > 5
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


def load_deltas(path, targets, genes, array="model"):
    """v3 prediction dir (genes.npy, targets.npy, <array>.npy) or .npz (delta [+targets, genes]) -> T x G, NaN -> 0."""
    p = Path(path)
    if p.is_dir():
        D = np.load(p / f"{array}.npy", mmap_mode="r"); t = np.load(p / "targets.npy").astype(str)
        g = np.load(p / "genes.npy").astype(str)
    else:
        z = np.load(p, allow_pickle=False); D = z["delta"]
        t = z["targets"].astype(str) if "targets" in z else np.array(targets)
        g = z["genes"].astype(str) if "genes" in z else np.array(genes)
    out = np.zeros((len(targets), len(genes)), np.float32)
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
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


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
def calibrate(D, line, targets, cfg):
    """Steps 0-1: alpha, clip, floor, 5-CPM mask; zero the knocked-down gene (its factor is handled separately)."""
    D = np.clip(cfg.alpha * D, -cfg.max_lfc, cfg.max_lfc).astype(np.float32)
    D[(D < 0) & ~line.gate[None, :]] = 0.0
    gp = {g: i for i, g in enumerate(line.genes)}
    for i, t in enumerate(targets):
        if t in gp:
            D[i, gp[t]] = 0.0
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
        if np.any(intensity > 0):
            dispersion = float(line.overdispersion[g])
            if dispersion > 1e-8:
                intensity = grng.gamma(1/dispersion, intensity*dispersion)
            values += grng.poisson(intensity)
        B[:, g] = values
    if B.max(initial=0) > np.iinfo(np.int32).max: raise ValueError("count exceeds int32 range")
    out = sp.csr_matrix(B.astype(np.int32)); out.eliminate_zeros(); validate_counts(out)
    planned = np.log2((projected*1e6+0.5)/(line.share*1e6+0.5))
    return out, dict(net_share_change=mass-1, projection_normalizer=mass, n_changed=int(len(cols)),
                     effective_pi=effective_pi, library_scale=lib_scale), planned


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
        chosen = json.loads(Path(a.settings).read_text())['settings']
        for key in ['alpha', 'pi', 'up_w', 'n_cells', 'metadata', 'match_pool', 'tilt']:
            if key in chosen: setattr(cfg, key, chosen[key])
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
            targets = list(specifications[name]); D = calibrate(load_deltas(deltas[name], targets, genes, a.delta_array), line, targets, cfg)
            own = prediction_ontarget(deltas[name], targets)
            if a.kd_orion or a.kd_tcell:
                legacy_own = kd_from_sources(targets, genes, a.kd_orion or (), a.kd_tcell, a.kd_tcell_condition)
                own = np.where(np.isfinite(own), own, legacy_own)
            kd = kd_residual(own, cfg); checks = []; diags = []
            # Pick max count once; different n_cells use prefixes of this same pool.
            pool_cfg = copy.copy(cfg); pool_cfg.n_cells = max(s['n_cells'] for s in specifications[name].values())
            common = pick_cells(line, np.zeros(len(genes)), pool_cfg, np.random.default_rng(seed_for(cfg.seed, name, 'base_controls')))
            for i, t in enumerate(targets):
                local = copy.copy(cfg); local.n_cells = specifications[name][t]['n_cells']; local.library_size = specifications[name][t]['median_umi_per_cell']
                pick = None if local.tilt else common[:local.n_cells]
                cells, diag, planned = generate_one(line, D[i], gp.get(t), kd[i], local, pick=pick, line_key=name, target_key=t)
                diags.append(dict(target=t, **diag))
                if i < a.checks: checks.append(dict(target=t, **answer_free_check(cells, line, planned, gp.get(t))))
                # float32 holds these integer counts exactly and matches vcc prep's output encoding,
                # so prep packages without a second full-matrix copy (int32 needs ~37 GB at 360k cells).
                if cells.data.max(initial=0) >= 2**24: raise ValueError('count too large for exact float32')
                cells = cells.astype(np.float32)
                if ds is None: write_elem(f, 'X', int64_indptr(cells)); ds = sparse_dataset(f['X'])
                else: ds.append(cells)
            pd.DataFrame(checks).to_csv(out/f'checks_{name}.csv', index=False)
            pd.DataFrame(diags).to_csv(out/f'projection_{name}.csv', index=False)
            report['lines'][name] = dict(n_targets=len(targets), n_cells=sum(s['n_cells'] for s in specifications[name].values()),
                metadata=specifications[name], measured_ontarget=int(np.isfinite(own).sum()))
            print(f"generated context {name}", flush=True)
    report['integer_counts'] = check_output(path, (len(obs), len(genes)))
    save_json(out/'generator_report.json', report); print(f'WROTE {path}; checked every output row')


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
        if a.n_targets: targets = targets[:a.n_targets]
        if not targets: raise ValueError('no eligible H1 targets')

        def subset(rows):
            rows = np.sort(np.asarray(rows))
            small = h[rows].to_memory(); validate_counts(small.X, 'H1 raw counts')
            return aligned_X(small, genes)

        ci = np.flatnonzero(labels == CONTROL)
        if len(ci) < 2: raise ValueError('H1 needs controls')
        ci = np.random.default_rng(seed_for(a.seed, 'H1', 'controls')).permutation(ci)
        train_rows, score_rows = ci[:len(ci) // 2], ci[len(ci) // 2:]
        line = Line(subset(train_rows), None, genes)
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
    """Runs official cell-eval2 against one reference file.

    The reference never changes between arms, so its expensive side (real DE tests,
    pseudobulks, noise corrections) is cached once with `prep-cache --side real` and reused by
    every arm through `--cache-real`. The ceiling depends on the real data only, so it is computed
    once in real-only mode instead of inside every arm.
    """
    def __init__(self, a, reference, out):
        self.binary = shutil.which(a.cell_eval)
        if self.binary is None: raise RuntimeError('install official cell-eval2 (pip install -r requirements-eval.txt) before running the H1 benchmark')
        help_text = subprocess.run([self.binary, 'run', '--help'], check=True, text=True, capture_output=True).stdout
        if '--ceiling' not in help_text or '--cache-real' not in help_text:
            raise RuntimeError('installed cell-eval2 lacks --ceiling/--cache-real; install requirements-eval.txt')
        self.a, self.reference, self.out = a, Path(reference), Path(out)
        self.cache = self.out / 'real_cache'
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

    def score(self, prediction, output):
        self._run([self.binary, 'run', '-ap', prediction, '-ar', self.reference, '-o', output,
                   '--cache-real', self.cache, *self.common])
        self.record(Path(output) / 'invocation.json')

    def record(self, path):
        import importlib.metadata
        versions = {}
        for package in ['cell-eval2', 'pdex', 'numpy', 'scipy', 'anndata']:
            try: versions[package] = importlib.metadata.version(package)
            except importlib.metadata.PackageNotFoundError: versions[package] = 'not installed in generator environment'
        save_json(path, dict(commands=self.invocations, versions=versions))


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


def cmd_benchmark(a):
    """Official H1 benchmark, streamed: reference and each arm go to disk one target at a time.

    Oracle mode never writes selected generator settings.
    """
    genes = read_list(a.genes, 'gene_name')
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    if a.mode == 'transfer':
        provenance = json.loads((Path(a.predictions)/'prediction_report.json').read_text())
        if provenance.get('heldout_background') != 'H1': raise ValueError('generator tuning requires LOBO-H1 predictions (H1 excluded)')
        if a.delta_array == 'selected': a.delta_array = provenance['selected_method']
    real_path = out/'real.h5ad'
    line, targets, real_counts, oracle_effects, oracle_own = stream_h1(a, genes, real_path)
    evaluator = OfficialEvaluator(a, real_path, out); evaluator.prepare()
    if a.mode == 'transfer':
        methods = list(dict.fromkeys([a.delta_array, 'model', 'source_mean']))
        kd = kd_residual(prediction_ontarget(a.predictions, targets), Cfg(a))
        effects = {method: load_deltas(a.predictions, targets, genes, method) for method in methods}
        variants = [dict(alpha=alpha, pi=pi, up_w=w, cell_rule=rule)
                    for alpha in a.alpha_grid for pi in a.pi_grid for w in a.up_w_grid for rule in a.cell_rules]
    else:
        # Ceiling-only oracle inputs: per-cell arithmetic CPM from each target's estimation half.
        effects = {'oracle': oracle_effects}; kd = kd_residual(oracle_own, Cfg(a))
        variants = [dict(alpha=a.alpha, pi=a.pi, up_w=a.up_w, cell_rule='matched')]
    outputs, settings = {}, {}
    maximum = max([a.n_cells] + list(real_counts.values()))
    pool_cfg = Cfg(a); pool_cfg.n_cells = maximum
    common = pick_cells(line, np.zeros(len(genes)), pool_cfg, np.random.default_rng(seed_for(a.seed, 'H1', 'base_controls')))
    for method, effects_array in effects.items():
        for index, variant in enumerate(variants):
            label = f'{method}_{index:03d}'
            cfg = Cfg(argparse.Namespace(**dict(vars(a), **{k: v for k, v in variant.items() if k != 'cell_rule'})))
            calibrated = calibrate(effects_array, line, targets, cfg)
            n_cells = {t: (real_counts[t] if variant['cell_rule'] == 'matched' else a.n_cells) for t in targets}
            pred_path = out/f'{label}.h5ad'
            write_arm(pred_path, line, calibrated, kd, targets, n_cells, common, cfg)
            del calibrated
            evaluator.score(pred_path, out/label)
            outputs[label] = str(out/label); settings[label] = dict(method=method, **variant)
            save_json(out/f'{label}.json', settings[label])
            print(f'official evaluation complete: {label}', flush=True)
    save_json(out/'benchmark.json', dict(mode=a.mode, outputs=outputs, settings=settings, n_targets=len(targets),
              ceiling=str(out/'ceiling'), real_cache=str(evaluator.cache),
              note='H1 is generator calibration; final evidence requires unseen lines. Oracle is a generator-only ceiling check.'))
    if a.metric_directions or a.eval_profile == 'vcc2026':
        directions = a.metric_directions or str(out/'metric_directions.csv')
        if a.metric_directions is None:
            pd.DataFrame({'metric': list(VCC2026_METRICS), 'direction': list(VCC2026_METRICS.values())}).to_csv(directions, index=False)
        frame, ranks = metric_rank_table(outputs, directions)
        frame.to_csv(out/'metric_ranks.csv', index=False); ranks.to_csv(out/'aggregate_rank.csv')
        if a.mode == 'transfer':
            # Compare model/source_mean officially, but tune knobs on the transfer
            # method selected by multi-background LOBO, not an H1-only method winner.
            eligible = [label for label in ranks.index if settings[str(label)]['method'] == a.delta_array and settings[str(label)]['cell_rule'] == 'fixed']
            if not eligible: raise ValueError('include fixed cell rule to select submission-compatible generator settings')
            best = settings[str(eligible[0])]
            save_json(out/'selected_generator.json', dict(settings=dict({k:v for k,v in best.items() if k not in ['method', 'cell_rule']}, n_cells=a.n_cells),
                       cell_rule='fixed', n_cells=a.n_cells, diagnostic_best_method=best['method'],
                       delta_array='selected', selection='official metric average rank on LOBO-H1 transfer predictions',
                       note='method remains controlled by multi-background LOBO; no per-line settings'))


def add_switches(p):
    p.add_argument('--alpha', type=float, default=1.0); p.add_argument('--n-cells', type=int, default=400)
    p.add_argument('--pi', type=float, default=0); p.add_argument('--up-w', type=float, default=0.75)
    p.add_argument('--metadata', action=argparse.BooleanOptionalAction, default=True)
    p.add_argument('--match-pool', action='store_true'); p.add_argument('--tilt', action='store_true')
    p.add_argument('--seed', type=int, default=2026)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('submit'); add_switches(s)
    for key in ['data', 'out']: s.add_argument('--'+key, required=True)
    s.add_argument('--lines', nargs='+', default=['A', 'B', 'C']); s.add_argument('--delta', action='append', required=True)
    s.add_argument('--required-n-cells', type=int, default=400)
    s.add_argument('--delta-array', default='selected'); s.add_argument('--settings')
    s.add_argument('--kd-orion', nargs='*'); s.add_argument('--kd-tcell'); s.add_argument('--kd-tcell-condition', default='Stim48hr')
    s.add_argument('--checks', type=int, default=20)
    for command, mode in [('benchmark', 'transfer'), ('oracle', 'oracle')]:
        b = sub.add_parser(command); add_switches(b); b.set_defaults(mode=mode)
        for key in ['cells', 'genes', 'out']: b.add_argument('--'+key, required=True)
        b.add_argument('--predictions', required=mode == 'transfer'); b.add_argument('--targets')
        b.add_argument('--obs-col', default='target_gene'); b.add_argument('--delta-array', default='selected')
        b.add_argument('--min-cells', type=int, default=100); b.add_argument('--n-targets', type=int, default=0)
        b.add_argument('--alpha-grid', nargs='+', type=float, default=[0.5, 0.75, 1.0])
        b.add_argument('--pi-grid', nargs='+', type=float, default=[0.0, 0.15])
        b.add_argument('--up-w-grid', nargs='+', type=float, default=[0.5, 0.75, 1.0])
        b.add_argument('--cell-rules', nargs='+', choices=['matched', 'fixed'], default=['fixed'])
        b.add_argument('--cell-eval', default='cell-eval2'); b.add_argument('--eval-profile', default='vcc2026')
        b.add_argument('--num-threads', type=int, default=4)
        b.add_argument('--metric-directions', help='CSV metric,direction for current official scoring metrics only')
    a = p.parse_args()
    if a.cmd == 'submit': cmd_submit(a)
    else: cmd_benchmark(a)

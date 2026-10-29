"""Stage-1 transfer model for VCC 2026 effect tables (effects_v5).

Predicts, for a target cell line, each (knockdown, gene):

  lfc_unshrunk     precision-weighted mean across other lines (weights 1/SE^2)
  v_model          uncertainty of that transferred effect (NOT measurement noise)
  se_measurement   noise of a 400-cell experiment in the *target* line
  lfc              shrunk magnitude used for NMAE / MSE / the generator
  de_sig           genes the reference is predicted to call (BH on
                   |shrunk LFC| / measurement SE; count varies by knockdown)
  p_sign           P(sign is right), for reach ranking

Control quantiles are not predicted. They are a property of the target line's
own controls (the v5 sidecar) and stay an input to the frozen cell generator.
This model only predicts the first moment and which genes are called.

Two uncertainties, never collapsed into one number
-----------------------------------------------
Measurement SE  — how noisy the gene is in a 400-cell target-line experiment.
  Depends on the target line's control CPM and the cell count. The reference
  calls gene g when |true LFC| / SE_measurement clears BH. Fitted as
  log(SE) = a + b log(control_mean_cpm + 1) + c log(n_cells), then evaluated
  at n_cells = 400 on the target line's baseline.

Model uncertainty — how far a new line's true effect can sit from the pooled
  estimate. Fixed-effect variance of the pool (1/SE^2 weights) plus a
  between-line tau^2. Tau^2 is gene-specific when lines disagree, and is never
  smaller than the typical between-line scatter, so a copied effect is not
  treated as exact.

Most genes do not move. A single Gaussian prior fitted on the non-null genes
is too wide to shrink those nulls, and NMAE then loses to predicting zero.
The magnitude handed to the generator is therefore the Gaussian shrinkage the
metrics ask for, multiplied by the posterior probability the gene is one of
the genes that actually moves:

  P(slab | data) * prediction * v0 / (v0 + v_model)

v0 is the variance of effects that already look non-null. The slab weight is
the empirical-Bayes probability under a point-null versus N(0, v0), with
v_model as the likelihood variance. Training is the precision-weighted
random-effects likelihood (each source weighted by 1/SE^2, times a modality
multiplier). No fixed top-N. The call count falls out of BH, so it varies by
knockdown.

Calls use |shrunk LFC| / measurement SE (what the reference divides by).
Model uncertainty enters the call through the shrinkage, not by inflating the
denominator a second time. Reach is ranked by P(sign correct) under the same
posterior: null genes stay near 1/2, confident slabs approach 1.

Modality multipliers (applied to source variance, so noisy tables vote less):
  single_cell, non-targeting reference     x1
  single_cell, rest_of_line (Feng GW)      x4
  CRISPRko                                 x4
  replicate_table (Marson) and drug        x16
Non-single-cell SEs are also floored at 0.25 because they are not 400-cell
Wilcoxon SEs. |LFC| > 6 is dropped (epsilon-floor artifacts).

The target line is never a source for its own knockdowns, and neither is any
other table of the same cell line (obs.background): holding out
replogle_k562_gw also drops jiang_*_k562 and adamson_k562_*. Priors (tau^2,
the shrinkage variance, the SE curve) are fit on the remaining lines only.

Duplicate target symbols in one table are precision-pooled, except in drug
tables: each Tahoe drug on a shared target is its own x16 observation, so
drugs that disagree raise tau^2 instead of averaging to a falsely precise
value. n_sources counts observations, so a Tahoe line can add more than one.

Locked run choices (do not silently change)
------------------------------------------
  sources     all of effects_v5 and effects_v5_extra; noisy tables are
              down-weighted by the multipliers above, not dropped
  knockdowns  scored lists only. Default, no flag needed:
              arc_h1            -> results/target_lists/h1_targets_300.txt
              replogle_k562_gw  -> results/target_lists/k562_heldout_evalB_100.txt
              looked up from the working directory, then from ~/vcc
  generator   the shrunk posterior mean, the `lfc` array. `lfc_unshrunk`
              is saved for comparison and must not be what the generator sees

Module form of vcc_2026_6.ipynb's model cells (moved verbatim 2026-10-27; the notebook imports
it). Self-test: ~/venv/bin/python code/vcc6/stage1.py
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

EPS = 1e-9
LFC_CAP = 6.0
SE_MIN = 1e-3
SE_FLOOR_NON_SC = 0.25
FDR = 0.05
N_CELLS_REF = 400
V0_FLOOR = 0.05
V0_CAP = 4.0
TAU2_FLOOR = 1e-4

# Default scored lists. Relative to the working directory, then to ~/vcc.
DEFAULT_LISTS = {
    "arc_h1": "results/target_lists/h1_targets_300.txt",
    "replogle_k562_gw": "results/target_lists/k562_heldout_evalB_100.txt",
}


def resolve_knockdowns_file(target, given):
    """Scored-list path for a known target. Explicit --knockdowns-file wins."""
    if given:
        path = Path(given).expanduser()
        if not path.exists():
            raise FileNotFoundError(path)
        return path
    rel = DEFAULT_LISTS.get(target)
    if rel is None:
        return None
    candidates = [Path(rel), Path.home() / "vcc" / rel]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(
        f"scored list for {target} not found. Looked for {rel} in the working "
        f"directory and in ~/vcc. Pass --knockdowns-file if it lives elsewhere."
    )
MULTIPLIER = {
    ("single_cell", "non-targeting"): 1.0,
    ("single_cell", "rest_of_line"): 4.0,
    ("CRISPRko", "non-targeting"): 4.0,
    ("CRISPRko", "rest_of_line"): 4.0,
    ("replicate_table", "non-targeting"): 16.0,
    ("replicate_table", "rest_of_line"): 16.0,
    ("drug", "non-targeting"): 16.0,
    ("drug", "rest_of_line"): 16.0,
}


def multiplier_for(modality, reference):
    mod = modality or "single_cell"
    ref = reference or "non-targeting"
    if (mod, ref) in MULTIPLIER:
        return MULTIPLIER[(mod, ref)]
    if mod in ("replicate_table", "drug"):
        return 16.0
    if mod == "CRISPRko" or ref == "rest_of_line":
        return 4.0
    return 8.0


def norm_sf(z):
    """P(Z > z) for standard normal Z. z is any shape.

    math.erfc, not np.erfc: this numpy build does not expose erfc.
    """
    z = np.asarray(z, dtype=np.float64)
    scaled = z / np.sqrt(2.0)
    flat = np.fromiter((math.erfc(float(v)) for v in scaled.ravel()), dtype=np.float64, count=z.size)
    return (0.5 * flat).reshape(z.shape)


def norm_cdf(z):
    return 1.0 - norm_sf(z)


def bh_reject(p, fdr=FDR):
    """Benjamini-Hochberg. p may contain NaN (not tested). Returns a bool mask."""
    p = np.asarray(p, dtype=np.float64)
    out = np.zeros(p.shape, dtype=bool)
    tested = np.isfinite(p)
    pv = p[tested]
    m = int(pv.size)
    if m == 0:
        return out
    order = np.argsort(pv, kind="mergesort")
    thresh = fdr * (np.arange(1, m + 1, dtype=np.float64) / m)
    below = np.flatnonzero(pv[order] <= thresh)
    if below.size == 0:
        return out
    keep = np.zeros(m, dtype=bool)
    keep[order[: int(below[-1]) + 1]] = True
    out[np.flatnonzero(tested)[keep]] = True
    return out


def auroc(scores, labels):
    """AUROC via the rank-sum formula. NaN if either class is empty."""
    scores = np.asarray(scores, dtype=np.float64)
    labels = np.asarray(labels, dtype=bool)
    ok = np.isfinite(scores)
    scores, labels = scores[ok], labels[ok]
    n_pos = int(labels.sum())
    n_neg = int((~labels).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(scores.size, dtype=np.float64)
    # average ranks for ties
    sorted_s = scores[order]
    i = 0
    while i < sorted_s.size:
        j = i + 1
        while j < sorted_s.size and sorted_s[j] == sorted_s[i]:
            j += 1
        ranks[order[i:j]] = 0.5 * (i + j - 1) + 1.0  # average of 1-based ranks
        i = j
    sum_pos = float(ranks[labels].sum())
    return (sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


class Source:
    """One v5 effect table, aligned later onto the prediction panel."""

    def __init__(self, name, background, modality, reference, genes, knockdowns,
                 lfc, se, n_cells, mult):
        self.name = str(name)
        self.background = str(background)
        self.modality = str(modality)
        self.reference = str(reference)
        self.genes = np.asarray(genes, dtype=object)
        self.knockdowns = np.asarray(knockdowns, dtype=object)
        self.lfc = np.asarray(lfc, dtype=np.float32)
        self.se = np.asarray(se, dtype=np.float32)
        self.n_cells = np.asarray(n_cells, dtype=np.float32)
        self.mult = float(mult)
        if self.lfc.shape != self.se.shape:
            raise ValueError(f"{name}: lfc/se shape mismatch")
        if self.lfc.shape[0] != self.knockdowns.shape[0]:
            raise ValueError(f"{name}: knockdown axis mismatch")
        self.rows_of = {}
        for i, k in enumerate(self.knockdowns):
            # one row per symbol after load_source, except drug tables: every drug on a shared
            # target stays its own observation of that gene
            self.rows_of.setdefault(str(k), []).append(i)
        self.col_of = {str(g): i for i, g in enumerate(self.genes)}


class SECalibration:
    """log(SE) = a + b log(control_cpm + 1) + c log(n_cells)."""

    def __init__(self, a, b, c, n, r):
        self.a = float(a)
        self.b = float(b)
        self.c = float(c)
        self.n = int(n)
        self.r = float(r)

    def predict(self, control_cpm, n_cells=N_CELLS_REF):
        m = np.asarray(control_cpm, dtype=np.float64)
        log_se = self.a + self.b * np.log(np.maximum(m, 0.0) + 1.0) + self.c * np.log(float(n_cells))
        return np.clip(np.exp(log_se), SE_MIN, 5.0)

    def as_dict(self):
        return dict(a=self.a, b=self.b, c=self.c, n=self.n, r=self.r)


class Prior:
    def __init__(self, tau2, v0, pi, n_tau, n_v0):
        self.tau2 = float(tau2)
        self.v0 = float(v0)
        self.pi = float(pi)
        self.n_tau = int(n_tau)
        self.n_v0 = int(n_v0)

    def as_dict(self):
        return dict(tau2=self.tau2, v0=self.v0, pi=self.pi, n_tau=self.n_tau, n_v0=self.n_v0)


def _collapse_duplicate_knockdowns(symbols, lfc, se, n_cells):
    """Precision-pool rows that share a target symbol. Returns uniques in first-seen order.
    Not applied to drug tables (see load_source)."""
    symbols = np.asarray(symbols, dtype=object)
    if len(set(map(str, symbols))) == len(symbols):
        return symbols, lfc, se, n_cells
    order = {}
    groups = []
    for i, s in enumerate(symbols):
        key = str(s)
        if key not in order:
            order[key] = len(groups)
            groups.append([i])
        else:
            groups[order[key]].append(i)
    out_s, out_y, out_e, out_n = [], [], [], []
    for idxs in groups:
        if len(idxs) == 1:
            i = idxs[0]
            out_s.append(symbols[i])
            out_y.append(lfc[i])
            out_e.append(se[i])
            out_n.append(n_cells[i])
            continue
        Y = np.stack([lfc[i] for i in idxs]).astype(np.float64)
        E = np.stack([se[i] for i in idxs]).astype(np.float64)
        ok = np.isfinite(Y) & np.isfinite(E) & (E > 0)
        w = np.where(ok, 1.0 / np.maximum(E, SE_MIN) ** 2, 0.0)
        sw = w.sum(0)
        mu = np.divide((w * np.where(ok, Y, 0.0)).sum(0), sw, out=np.full(Y.shape[1], np.nan), where=sw > 0)
        var = np.divide(1.0, sw, out=np.full(Y.shape[1], np.nan), where=sw > 0)
        out_s.append(symbols[idxs[0]])
        out_y.append(mu.astype(np.float32))
        out_e.append(np.sqrt(var).astype(np.float32))
        out_n.append(float(np.nanmedian([n_cells[i] for i in idxs])))
    return (np.asarray(out_s, dtype=object), np.stack(out_y), np.stack(out_e),
            np.asarray(out_n, dtype=np.float32))


def pool_effects(Y, SE, mults, tau2_global):
    """Precision-weighted pool across sources. Weights are 1 / (SE^2 * multiplier).

    Returns mu, var_fe, tau2, tau_untrunc, n.
    var_fe is the sampling variance of mu as an estimate of the shared effect.
    tau2 is the between-line variance used for a new line: gene-specific
    DerSimonian-Laird, floored at tau2_global so a tight pair of lines is not
    treated as more transferable than the typical gene. tau_untrunc is the
    untruncated DL moment (can be negative), for estimating tau2_global.
    Genes with no finite source have var_fe = inf and n = 0.
    """
    Y = np.asarray(Y, dtype=np.float64)
    SE = np.asarray(SE, dtype=np.float64)
    S, G = Y.shape
    finite = np.isfinite(Y) & np.isfinite(SE) & (SE > 0) & (np.abs(Y) <= LFC_CAP)
    var = (np.maximum(SE, SE_MIN) ** 2) * np.asarray(mults, dtype=np.float64)[:, None]
    w = np.where(finite, 1.0 / var, 0.0)
    sw = w.sum(0)
    n = finite.sum(0).astype(np.int32)
    mu = np.divide((w * np.where(finite, Y, 0.0)).sum(0), sw, out=np.zeros(G), where=sw > 0)
    var_fe = np.divide(1.0, sw, out=np.full(G, np.inf), where=sw > 0)
    resid = np.where(finite, Y - mu, 0.0)
    Q = (w * resid ** 2).sum(0)
    C = sw - np.divide((w ** 2).sum(0), sw, out=np.zeros(G), where=sw > 0)
    tau_untrunc = np.full(G, np.nan)
    ok = (n >= 2) & (C > 1e-12)
    tau_untrunc[ok] = (Q[ok] - (n[ok] - 1.0)) / C[ok]
    tau_dl = np.where(ok, np.maximum(0.0, tau_untrunc), 0.0)
    tau2 = np.maximum(tau_dl, float(tau2_global))
    tau2 = np.where(n >= 1, tau2, np.nan)
    return mu, var_fe, tau2, tau_untrunc, n


def combine_within(Y, SE, mults):
    """One background's observations (rows: tables, conditions, drugs) -> one vote per gene.

    DerSimonian-Laird random effects over the rows (weights 1 / (SE^2 * multiplier), pool_effects with no
    tau^2 floor): variance = var_fe + tau^2_DL, so rows that disagree (an agonist and an antagonist) inflate it.
    Floored per gene at the smallest single row's variance, so repeat screens of one line add no precision.
    Returns (mu, se, n_rows); se is NaN where no row is finite."""
    Y = np.asarray(Y, dtype=np.float64)
    SE = np.asarray(SE, dtype=np.float64)
    m = np.asarray(mults, dtype=np.float64)[:, None]
    mu, var_fe, tau2, _, n = pool_effects(Y, SE, m[:, 0], tau2_global=0.0)
    finite = np.isfinite(Y) & np.isfinite(SE) & (SE > 0) & (np.abs(Y) <= LFC_CAP)
    var_row = np.where(finite, np.maximum(SE, SE_MIN) ** 2 * m, np.inf)
    v = np.maximum(var_fe + np.where(np.isfinite(tau2), tau2, 0.0), var_row.min(0))
    se = np.where(n > 0, np.sqrt(v), np.nan)
    return np.where(n > 0, mu, np.nan), se, n


# Step 2 arms B / C: Orion's pseudobulk tables are not in MULTIPLIER (x8 by default in arm A).
MULTIPLIER_VOTES = {**MULTIPLIER, ("pseudobulk_table", "non-targeting"): 16.0,
                    ("pseudobulk_table", "rest_of_line"): 16.0}


class BackgroundVotes:
    """One vote per background (obs.background, case-insensitive) for arms B and C.

    gather(k, genes) combines every row of every table of a background that holds knockdown k with
    combine_within, and returns one (Y, SE) row per background, multiplier 1 (each row's own multiplier is
    already inside its background's variance). pool_effects across these rows is the existing random-effects
    pool, so n_sources counts backgrounds."""

    def __init__(self, sources):
        self.groups = {}
        for s in sources:
            self.groups.setdefault(_line_key(s.background), []).append(s)
        self.owners = {}
        for bg, ss in self.groups.items():
            for s in ss:
                for k in s.rows_of:
                    lst = self.owners.setdefault(k, [])
                    if bg not in lst:
                        lst.append(bg)

    @staticmethod
    def mult(src):
        return MULTIPLIER_VOTES.get((src.modality, src.reference), src.mult)

    def gather(self, knockdown, genes):
        Y, SE, mults = [], [], []
        for bg in self.owners.get(str(knockdown), ()):
            ys, es, ms = [], [], []
            for src in self.groups[bg]:
                for row in src.rows_of.get(str(knockdown), ()):
                    y, e = _align_row(src, row, genes)
                    ys.append(y)
                    es.append(e)
                    ms.append(self.mult(src))
            mu, se, _ = combine_within(np.stack(ys), np.stack(es), np.array(ms))
            Y.append(mu)
            SE.append(se)
            mults.append(1.0)
        return Y, SE, mults


def mixture_posterior(mu, var_fe, tau2, v0, pi):
    """Shrink a pooled effect toward 0 with a null-vs-slab prior.

    v_model = var_fe + tau2 is the uncertainty of the transferred effect.
    Slab N(0, v0), point null at 0, prior weight `pi` on the slab.
    Returned LFC is P(slab|data) * mu * v0 / (v0 + v_model).
    """
    mu = np.asarray(mu, dtype=np.float64)
    var_fe = np.asarray(var_fe, dtype=np.float64)
    tau2 = np.asarray(tau2, dtype=np.float64)
    known = np.isfinite(var_fe) & np.isfinite(mu) & (var_fe >= 0)
    tau = np.where(np.isfinite(tau2) & (tau2 > 0), tau2, 0.0)
    v_model = np.where(known, var_fe + tau, np.inf)
    v = np.where(known, np.maximum(v_model, 1e-8), 1.0)
    mu_slab = mu * (v0 / (v0 + v))
    v_slab = (v0 * v) / (v0 + v)
    # log likelihood ratio of slab vs point null, clipped so tiny SEs cannot overflow
    log_lr = (math.log(pi) - math.log(1.0 - pi)
              - 0.5 * np.log1p(v0 / v)
              + 0.5 * (mu * mu) * (v0 / (v * (v0 + v))))
    log_lr = np.clip(log_lr, -40.0, 40.0)
    p_slab = 1.0 / (1.0 + np.exp(-log_lr))
    p_slab = np.where(known, p_slab, 0.0)
    mu_post = np.where(known, p_slab * mu_slab, 0.0)
    second = p_slab * (v_slab + mu_slab ** 2)
    v_post = np.where(known, np.maximum(second - mu_post ** 2, 0.0), v0)
    sd_slab = np.sqrt(np.maximum(v_slab, TAU2_FLOOR))
    p_pos = p_slab * norm_cdf(mu_slab / sd_slab) + (1.0 - p_slab) * 0.5
    p_sign = np.where(mu_post >= 0.0, p_pos, 1.0 - p_pos)
    p_sign = np.where(known, np.clip(p_sign, 0.5, 1.0), 0.5)
    return mu_post, v_post, v_model, p_sign, p_slab


def predictive_calls(mu_post, se_meas, tested, fdr=FDR):
    """BH on z = shrunk LFC / measurement SE.

    The reference standardizes by measurement SE, not by our transfer uncertainty.
    Transfer uncertainty is already in the shrunk numerator: null genes are near
    0 and drop out of the call set, so the count varies by knockdown.
    """
    se = np.maximum(np.asarray(se_meas, dtype=np.float64), SE_MIN)
    z = np.divide(np.asarray(mu_post, dtype=np.float64), se, out=np.zeros_like(se), where=se > 0)
    p = np.full(z.shape, np.nan)
    use = np.asarray(tested, dtype=bool) & np.isfinite(z)
    p[use] = np.minimum(1.0, 2.0 * norm_sf(np.abs(z[use])))
    called = bh_reject(p, fdr=fdr)
    return called, p


def _align_row(source, row, genes):
    """One knockdown row reindexed onto `genes` (NaN where the source lacks the gene)."""
    idx = np.array([source.col_of.get(str(g), -1) for g in genes], dtype=np.int32)
    y = np.full(len(genes), np.nan, dtype=np.float64)
    e = np.full(len(genes), np.nan, dtype=np.float64)
    ok = idx >= 0
    y[ok] = source.lfc[row, idx[ok]]
    e[ok] = source.se[row, idx[ok]]
    if source.modality != "single_cell":
        e = np.where(np.isfinite(e), np.maximum(e, SE_FLOOR_NON_SC), e)
    return y, e


def _prepare_se(se, modality):
    out = np.array(se, dtype=np.float64, copy=True)
    if modality != "single_cell":
        out = np.where(np.isfinite(out), np.maximum(out, SE_FLOOR_NON_SC), out)
    return out


def fit_se_calibration(sources, basals, max_rows=25, max_genes=2500, seed=2026):
    """Fit the measurement-SE curve on source lines. `basals` maps source name -> (genes, cpm)."""
    rng = np.random.default_rng(seed)
    rows_x, rows_y = [], []
    for src in sources:
        if src.modality != "single_cell":
            continue  # only real 400-cell-style SEs identify the measurement curve
        basal = basals.get(src.name)
        if basal is None:
            continue
        b_genes, cpm = basal
        col = {str(g): i for i, g in enumerate(b_genes)}
        take_rows = np.arange(src.lfc.shape[0])
        if take_rows.size > max_rows:
            take_rows = rng.choice(take_rows, max_rows, replace=False)
        for r in take_rows:
            se = src.se[r].astype(np.float64)
            n = float(src.n_cells[r]) if np.isfinite(src.n_cells[r]) and src.n_cells[r] >= 20 else np.nan
            if not np.isfinite(n):
                continue
            m0 = np.array([cpm[col[str(g)]] if str(g) in col else np.nan for g in src.genes], dtype=np.float64)
            ok = np.isfinite(se) & (se > SE_MIN) & (se < 5) & np.isfinite(m0) & (m0 > 0)
            ix = np.flatnonzero(ok)
            if ix.size > max_genes:
                ix = rng.choice(ix, max_genes, replace=False)
            if ix.size == 0:
                continue
            x = np.column_stack([
                np.ones(ix.size),
                np.log(m0[ix] + 1.0),
                np.full(ix.size, np.log(n)),
            ])
            rows_x.append(x)
            rows_y.append(np.log(se[ix]))
    if not rows_x:
        # depth-only fallback: SE proportional to 1/sqrt(n), mild expression slope
        return SECalibration(a=np.log(0.8), b=-0.25, c=-0.5, n=0, r=np.nan)
    X = np.vstack(rows_x)
    y = np.concatenate(rows_y)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coef
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return SECalibration(a=coef[0], b=coef[1], c=coef[2], n=int(y.size), r=r)


def _fit_pi(mu, v_model, v0):
    """One-dimensional MLE of the slab weight. Subsampled; grid on (0, 0.6]."""
    mu = np.asarray(mu, dtype=np.float64)
    v = np.asarray(v_model, dtype=np.float64)
    ok = np.isfinite(mu) & np.isfinite(v) & (v > 0)
    mu, v = mu[ok], v[ok]
    if mu.size == 0:
        return 0.05
    if mu.size > 20000:
        rng = np.random.default_rng(0)
        take = rng.choice(mu.size, 20000, replace=False)
        mu, v = mu[take], v[take]
    v = np.maximum(v, 1e-8)
    # log LR without the pi term: slab density minus null density
    log_lr0 = -0.5 * np.log1p(v0 / v) + 0.5 * (mu * mu) * (v0 / (v * (v0 + v)))
    best_pi, best_ll = 0.05, -np.inf
    for pi in np.linspace(0.005, 0.60, 40):
        log_lr = np.clip(math.log(pi) - math.log(1.0 - pi) + log_lr0, -40.0, 40.0)
        # log(pi m_s + (1-pi) m_n) = log m_n + softplus(log LR)
        # constant log m_n doesn't change the argmax, but softplus(log LR) does
        ll = float((math.log(1.0 - pi) + np.logaddexp(0.0, log_lr)).sum())
        if ll > best_ll:
            best_ll, best_pi = ll, float(pi)
    return best_pi


def fit_prior(sources, max_knockdowns=150, seed=2026, votes=None):
    """tau^2, slab variance v0, and slab weight pi, from knockdowns seen in >= 2 lines.

    votes: a BackgroundVotes (arms B / C): knockdowns seen in >= 2 backgrounds, one vote per background,
    so the moments come from cross-background pairs only. None: every table row (arm A).

    v0 is the mean of (mu^2 - var_fe) over genes with |mu| / sqrt(var_fe) > 2.
    Fitting it on every gene would set v0 near 0 and over-shrink the real effects.
    tau^2 is a trimmed mean of the untruncated DL moment. The median is ~0, because
    measurement noise hides between-line scatter for the typical null gene, and a
    zero tau^2 would treat a copied effect as exact.
    pi is the marginal MLE of how many genes are in the slab.
    """
    rng = np.random.default_rng(seed)
    owners = {}
    if votes is not None:
        owners = votes.owners
    else:
        for src in sources:
            for k in src.rows_of:
                owners.setdefault(k, []).append(src)
    shared = [k for k, ss in owners.items() if len(ss) >= 2]
    if not shared:
        return Prior(tau2=0.05, v0=0.25, pi=0.05, n_tau=0, n_v0=0)
    if len(shared) > max_knockdowns:
        shared = list(rng.choice(np.array(shared, dtype=object), max_knockdowns, replace=False))
    tau_chunks = []
    signal = []
    mu_keep = []
    v_keep = []
    panel = sources[0].genes
    for k in shared:
        ss = owners[str(k)]
        Y, SE, mults = [], [], []
        if votes is not None:
            Y, SE, mults = votes.gather(k, panel)
        else:
            for src in ss:
                for row in src.rows_of[str(k)]:
                    y, e = _align_row(src, row, panel)
                    Y.append(y)
                    SE.append(e)
                    mults.append(src.mult)
        mu, var_fe, tau2, tau_untrunc, n = pool_effects(
            np.stack(Y), np.stack(SE), np.array(mults), tau2_global=0.0)
        use_tau = (n >= 2) & np.isfinite(tau_untrunc)
        if use_tau.any():
            tau_chunks.append(tau_untrunc[use_tau])
        z = np.divide(np.abs(mu), np.sqrt(np.maximum(var_fe, TAU2_FLOOR)),
                      out=np.zeros_like(mu), where=np.isfinite(var_fe))
        sig = (n >= 2) & (z > 2.0) & np.isfinite(mu) & np.isfinite(var_fe)
        if sig.any():
            signal.append(np.clip(mu[sig] ** 2 - var_fe[sig], 0.0, 25.0))
        use = (n >= 2) & np.isfinite(mu) & np.isfinite(var_fe)
        if use.any():
            mu_keep.append(mu[use])
            # v_model used at prediction time floors tau at the global estimate.
            # pi is fit once that floor is known; stash var_fe and add it below.
            v_keep.append(var_fe[use] + np.where(np.isfinite(tau2[use]), tau2[use], 0.0))
    if tau_chunks:
        raw = np.concatenate(tau_chunks)
        lo, hi = np.quantile(raw, [0.10, 0.90])
        mid = raw[(raw >= lo) & (raw <= hi)]
        tau2 = float(mid.mean()) if mid.size else float(np.median(raw))
        n_tau = int(raw.size)
    else:
        tau2, n_tau = 0.05, 0
    tau2 = float(np.clip(tau2, TAU2_FLOOR, 4.0))
    if signal:
        v0 = float(np.mean(np.concatenate(signal)))
        n_v0 = int(sum(len(s) for s in signal))
    else:
        v0, n_v0 = 0.25, 0
    v0 = float(np.clip(v0, V0_FLOOR, V0_CAP))
    if mu_keep:
        # raise the likelihood variance to the global floor used in prediction
        mus = np.concatenate(mu_keep)
        vs = np.maximum(np.concatenate(v_keep), tau2)
        pi = _fit_pi(mus, vs, v0)
    else:
        pi = 0.05
    return Prior(tau2=tau2, v0=v0, pi=pi, n_tau=n_tau, n_v0=n_v0)


def gaussian_posterior(mu, var_fe, tau2, v0):
    """The Gaussian shrink without the slab weight (magnitude "gaussian", G): mu x v0 / (v0 + v_model), v_model =
    var_fe + tau2 as in mixture_posterior. Returns (lfc, v_post, v_model, p_sign); p_sign from N(lfc, v_post)."""
    mu = np.asarray(mu, dtype=np.float64)
    var_fe = np.asarray(var_fe, dtype=np.float64)
    tau2 = np.asarray(tau2, dtype=np.float64)
    known = np.isfinite(var_fe) & np.isfinite(mu) & (var_fe >= 0)
    tau = np.where(np.isfinite(tau2) & (tau2 > 0), tau2, 0.0)
    v_model = np.where(known, var_fe + tau, np.inf)
    v = np.where(known, np.maximum(v_model, 1e-8), 1.0)
    lfc = np.where(known, mu * (v0 / (v0 + v)), 0.0)
    v_post = np.where(known, (v0 * v) / (v0 + v), v0)
    p_sign = np.where(known, norm_cdf(np.abs(lfc) / np.sqrt(np.maximum(v_post, TAU2_FLOOR))), 0.5)
    return lfc, v_post, v_model, p_sign


def predict_knockdown(sources, knockdown, genes, basal_cpm, tested, calibration, prior,
                      drop_gene=None, votes=None, magnitude="slab"):
    """Predict one knockdown onto `genes`. `drop_gene` (the knockdown symbol) is not called.
    votes: a BackgroundVotes (arms B / C) in place of one row per table row (arm A).
    magnitude: "slab" (P(slab) x Gaussian shrink, vcc_2026_6) or "gaussian" (the Gaussian shrink only, G: the
    stage-1 baseline with arm C since 2026-10-27). The calls follow the lfc either way."""
    Y, SE, mults = [], [], []
    if votes is not None:
        Y, SE, mults = votes.gather(knockdown, genes)
    else:
        for src in sources:
            for row in src.rows_of.get(str(knockdown), ()):
                y, e = _align_row(src, row, genes)
                Y.append(y)
                SE.append(e)
                mults.append(src.mult)
    G = len(genes)
    if not Y:
        mu = np.zeros(G)
        var_fe = np.full(G, np.inf)
        tau2 = np.full(G, np.nan)
        n_src = np.zeros(G, dtype=np.int32)
    else:
        mu, var_fe, tau2, _, n_src = pool_effects(
            np.stack(Y), np.stack(SE), np.asarray(mults), prior.tau2)
    tau2_used = np.where(np.isfinite(tau2), tau2, prior.tau2)
    if magnitude == "gaussian":
        mu_post, v_post, v_model, p_sign = gaussian_posterior(mu, var_fe, tau2_used, prior.v0)
    elif magnitude == "slab":
        mu_post, v_post, v_model, p_sign, _p_slab = mixture_posterior(mu, var_fe, tau2_used, prior.v0, prior.pi)
    else:
        raise ValueError(f"magnitude must be slab or gaussian, got {magnitude!r}")
    se_meas = calibration.predict(basal_cpm, N_CELLS_REF)
    tested = np.asarray(tested, dtype=bool).copy()
    if drop_gene is not None:
        for i, g in enumerate(genes):
            if str(g) == str(drop_gene):
                tested[i] = False
                break
    called, p = predictive_calls(mu_post, se_meas, tested)
    return dict(lfc=mu_post.astype(np.float32), lfc_unshrunk=mu.astype(np.float32),
                v_model=v_model.astype(np.float32), v_post=v_post.astype(np.float32),
                se_measurement=se_meas.astype(np.float32), de_sig=called, p_value=p.astype(np.float32),
                p_sign=p_sign.astype(np.float32), n_sources=n_src.astype(np.int16))


def score_against_truth(pred, truth_lfc, truth_se, truth_sig, tested):
    """Effect-space scores for one knockdown. truth_* aligned to the panel."""
    y = np.asarray(truth_lfc, dtype=np.float64)
    se = np.asarray(truth_se, dtype=np.float64)
    sig = np.asarray(truth_sig, dtype=bool)
    ok = np.isfinite(y) & np.asarray(tested, dtype=bool)
    out = dict(n=int(ok.sum()), n_true=int((sig & ok).sum()), n_called=int(pred["de_sig"][ok].sum()) if ok.any() else 0)
    if ok.sum() == 0:
        return out
    err = pred["lfc"].astype(np.float64) - y
    raw = pred["lfc_unshrunk"].astype(np.float64) - y
    denom = np.abs(y[ok]).sum() + EPS
    out.update(
        mae=float(np.abs(err[ok]).mean()),
        nmae=float(np.abs(err[ok]).sum() / denom),
        mse=float((err[ok] ** 2).mean()),
        mae_unshrunk=float(np.abs(raw[ok]).mean()),
        nmae_unshrunk=float(np.abs(raw[ok]).sum() / denom),
        mse_unshrunk=float((raw[ok] ** 2).mean()),
        nmae_zero=float(np.abs(y[ok]).sum() / denom),
        mse_zero=float((y[ok] ** 2).mean()),
    )
    if (sig & ok).any():
        s = sig & ok
        out["sign_acc_true_de"] = float((np.sign(pred["lfc"][s]) == np.sign(y[s])).mean())
        out["sign_acc_unshrunk"] = float((np.sign(pred["lfc_unshrunk"][s]) == np.sign(y[s])).mean())
    # calibration of the predictive distribution against the observed LFC
    sd = np.sqrt(np.maximum(pred["v_post"].astype(np.float64), 0.0) + np.maximum(se, 0.0) ** 2)
    both = ok & np.isfinite(se) & (sd > 0)
    if both.any():
        out["coverage_95"] = float((np.abs(pred["lfc"][both] - y[both]) <= 1.96 * sd[both]).mean())
    # discrimination: |posterior| / measurement SE vs the reference call
    score = np.abs(pred["lfc"].astype(np.float64)) / np.maximum(pred["se_measurement"].astype(np.float64), SE_MIN)
    if ok.any() and (sig & ok).any() and (~sig & ok).any():
        out["auroc"] = float(auroc(score[ok], sig[ok]))
    return out


def calibration_stats(pred, truth, tested):
    """Fit of predicted to true LFC, pooled over knockdowns and over each one's tested genes with a
    finite truth. pred: (K, G); truth: {row: true lfc (G,)}. pearson_r; slope_true_on_pred
    (calibration: 1 = right scale, < 1 = overconfident); slope_pred_on_true (transmission:
    1 = full size, < 1 = attenuated). Pass the shrunk lfc or lfc_unshrunk."""
    xs, ys = [], []
    gate = np.asarray(tested, dtype=bool)
    for i, y in truth.items():
        y = np.asarray(y, dtype=np.float64)
        x = np.asarray(pred[i], dtype=np.float64)
        ok = gate & np.isfinite(y) & np.isfinite(x)
        xs.append(x[ok])
        ys.append(y[ok])
    x = np.concatenate(xs) if xs else np.zeros(0)
    y = np.concatenate(ys) if ys else np.zeros(0)
    out = dict(n=int(x.size), pearson_r=None, slope_true_on_pred=None, slope_pred_on_true=None,
               beta=None, r2_fit=None, r2_raw=None)
    if x.size < 3:
        return out
    vx, vy = float(np.var(x)), float(np.var(y))
    cov = float(np.mean((x - x.mean()) * (y - y.mean())))
    if vx > 0 and vy > 0:
        r = cov / math.sqrt(vx * vy)
        out.update(pearson_r=r, slope_true_on_pred=cov / vx, slope_pred_on_true=cov / vy,
                   # beta: the factor to scale the predictions by (true regressed on predicted); r2_fit: the share of
                   # the true lfc's variance explained once rescaled by beta; r2_raw: the share explained by the
                   # predictions as they are (1 - SS_res / SS_tot; negative when worse than predicting the mean)
                   beta=cov / vx, r2_fit=r * r, r2_raw=1.0 - float(np.mean((y - x) ** 2)) / vy)
    return out


def _nanmean(xs):
    v = np.array([x for x in xs if x is not None and np.isfinite(x)], dtype=np.float64)
    return float(v.mean()) if v.size else None


def _nanmedian(xs):
    v = np.array([x for x in xs if x is not None and np.isfinite(x)], dtype=np.float64)
    return float(np.median(v)) if v.size else None


def summarize_scores(rows):
    keys = ["mae", "nmae", "mse", "nmae_unshrunk", "mse_unshrunk", "nmae_zero", "mse_zero",
            "sign_acc_true_de", "coverage_95", "auroc", "n_called", "n_true"]
    out = {k: _nanmedian([r.get(k) for r in rows]) for k in keys}
    # call-count tracking uses the per-knockdown vectors, not the median
    called = np.array([r.get("n_called", np.nan) for r in rows], dtype=np.float64)
    truth = np.array([r.get("n_true", np.nan) for r in rows], dtype=np.float64)
    m = np.isfinite(called) & np.isfinite(truth)
    if m.sum() >= 3 and truth[m].std() > 0 and called[m].std() > 0:
        out["call_count_pearson"] = float(np.corrcoef(called[m], truth[m])[0, 1])
    else:
        out["call_count_pearson"] = None
    out["knockdowns"] = len(rows)
    return out


# ---------------------------------------------------------------- I/O

def _as_str_array(values):
    return np.array(["" if v is None else str(v) for v in values], dtype=object)


def table_meta(obs, uns, stem):
    """(modality, reference, background) as the loader reads them from a table's obs (and uns)."""
    uns = dict(uns.get("effects_v5", {}))
    if "modality" in obs.columns:
        mods = set(map(str, obs["modality"].unique()))
        modality = mods.pop() if len(mods) == 1 else "single_cell"
    else:
        modality = "single_cell"
    if "reference" in obs.columns:
        refs = set(map(str, obs["reference"].unique()))
        reference = refs.pop() if len(refs) == 1 else "non-targeting"
    else:
        reference = str(uns.get("reference", "non-targeting"))
    if "de_labels_known" in obs.columns and not bool(np.all(obs["de_labels_known"].to_numpy())):
        if modality == "single_cell":
            modality = "replicate_table"
    background = str(obs["background"].iloc[0]) if "background" in obs.columns else stem
    return modality, reference, background


def load_source(path):
    """Read one effects_v5 h5ad. Requires anndata (present on the analysis VM)."""
    import anndata as ad
    path = Path(path)
    a = ad.read_h5ad(path)
    try:
        if "lfc" not in a.layers or "se" not in a.layers:
            raise ValueError(f"{path.name}: missing lfc/se layers")
        genes = _as_str_array(a.var_names)
        if "target_gene" not in a.obs.columns:
            raise ValueError(f"{path.name}: obs has no target_gene")
        knockdowns = _as_str_array(a.obs["target_gene"].tolist())
        modality, reference, background = table_meta(a.obs, a.uns, path.stem)
        n_cells = (a.obs["n_cells_used"].to_numpy(dtype=np.float64)
                   if "n_cells_used" in a.obs.columns else np.full(a.n_obs, np.nan))
        src = Source(
            name=path.stem, background=background, modality=modality, reference=reference,
            genes=genes, knockdowns=knockdowns,
            lfc=np.asarray(a.layers["lfc"], dtype=np.float32),
            se=np.asarray(a.layers["se"], dtype=np.float32),
            n_cells=n_cells.astype(np.float32),
            mult=multiplier_for(modality, reference),
        )
        if modality == "drug":
            # Tahoe: drugs sharing a target are different drugs, not replicates (an agonist and an
            # antagonist would cancel). Each stays its own x16 observation; their disagreement
            # goes into the between-observation tau^2 instead of a falsely small pooled SE.
            return src
        symbols, lfc, se, n_cells = _collapse_duplicate_knockdowns(
            src.knockdowns, src.lfc, src.se, src.n_cells)
        if len(symbols) != src.lfc.shape[0]:
            src = Source(src.name, src.background, src.modality, src.reference, src.genes,
                         symbols, lfc, se, n_cells, src.mult)
        return src
    finally:
        if getattr(a, "file", None) is not None:
            a.file.close()


def load_baseline(path, genes):
    """control_mean_cpm and de_tested aligned to `genes`. Missing file -> None."""
    path = Path(path)
    if not path.exists():
        return None
    import pandas as pd
    df = pd.read_csv(path, index_col=0)
    cpm_col = "control_mean_cpm" if "control_mean_cpm" in df.columns else df.columns[0]
    cpm = df[cpm_col]
    tested = df["de_tested"] if "de_tested" in df.columns else pd.Series(True, index=df.index)
    cpm_map = {str(i): float(v) for i, v in cpm.items() if np.isfinite(v)}
    tested_map = {str(i): bool(v) for i, v in tested.items()}
    out_cpm = np.array([cpm_map.get(str(g), np.nan) for g in genes], dtype=np.float64)
    out_tested = np.array([tested_map.get(str(g), False) for g in genes], dtype=bool)
    # a gene with no de_tested column is testable when it is expressed
    if "de_tested" not in df.columns:
        out_tested = np.isfinite(out_cpm) & (out_cpm > 5.0)
    else:
        # still require the 5-CPM gate even if the column says tested
        out_tested = out_tested & np.isfinite(out_cpm) & (out_cpm > 5.0)
    return out_cpm, out_tested


def discover_tables(dirs):
    found = []
    for d in dirs:
        root = Path(d)
        if not root.exists():
            raise FileNotFoundError(root)
        files = sorted(root.glob("*.h5ad")) + sorted(root.glob("*/*.h5ad"))
        for p in files:
            if p.name.endswith("_qc.h5ad"):
                continue
            found.append(p)
    # unique, stable
    seen, out = set(), []
    for p in found:
        key = str(p.resolve())
        if key not in seen:
            seen.add(key)
            out.append(p)
    return out


def truth_signature(path, knockdowns, genes):
    """Target-line LFC, SE and de_sig for scoring, aligned to genes. None if absent."""
    import anndata as ad
    a = ad.read_h5ad(path)
    try:
        if "lfc" not in a.layers:
            return None
        kds = _as_str_array(a.obs["target_gene"].tolist())
        row = {str(k): i for i, k in enumerate(kds)}
        col_index = {str(g): i for i, g in enumerate(_as_str_array(a.var_names))}
        cols = np.array([col_index.get(str(g), -1) for g in genes], dtype=np.int32)
        sig = a.layers["de_sig"].tocsr() if hasattr(a.layers["de_sig"], "tocsr") else None
        out = {}
        for k in knockdowns:
            i = row.get(str(k))
            if i is None:
                continue
            y = np.full(len(genes), np.nan)
            e = np.full(len(genes), np.nan)
            ok = cols >= 0
            y[ok] = np.asarray(a.layers["lfc"][i], dtype=np.float64)[cols[ok]]
            e[ok] = np.asarray(a.layers["se"][i], dtype=np.float64)[cols[ok]]
            called = np.zeros(len(genes), dtype=bool)
            if sig is not None:
                src_cols = sig[i].indices
                keep = cols >= 0
                pos = {int(c): j for j, c in enumerate(cols) if keep[j]}
                called[[pos[int(c)] for c in src_cols if int(c) in pos]] = True
            elif "de_sig" in a.layers:
                dense = np.asarray(a.layers["de_sig"][i]).ravel().astype(bool)
                called[ok] = dense[cols[ok]]
            out[str(k)] = (y, e, called)
        return out
    finally:
        if getattr(a, "file", None) is not None:
            a.file.close()


# One cell line under two names (Jose, 2026-10-27): Tahoe's HT-29 is Jiang's HT29; HepG2/C3A is a HepG2 subclone,
# so it never votes on a HepG2 hold-out. The only merges; kolf_2 / kolf_3 / KOLF2.1J and the Feng clone pairs stay apart.
BACKGROUND_ALIASES = {"ht-29": "ht29", "hepg2/c3a": "hepg2"}


def _line_key(background):
    key = str(background).strip().lower()
    return BACKGROUND_ALIASES.get(key, key)


def split_sources(sources, target):
    """(target source, the sources that may vote, names dropped as the target's cell line).

    A source is dropped when its background is the target's cell line, not only when its stem
    is the target: while replogle_k562_gw is held out, jiang_*_k562 and adamson_k562_* are
    K562 knockdowns of the held-out line and do not vote. They still vote for other targets.
    """
    target_src = next(s for s in sources if s.name == target)
    line = _line_key(target_src.background)
    others, same_line = [], []
    for s in sources:
        if s.name == target:
            continue
        (same_line if _line_key(s.background) == line else others).append(s)
    return target_src, others, [s.name for s in same_line]


def load_line(effects_dirs, target, crispri_only=False):
    """Load every table once, hold out the target line (and its cell line), read the baselines. The result
    feeds predict_line, so several arms can run on one load."""
    tables = discover_tables(effects_dirs)
    if not tables:
        raise FileNotFoundError(f"no h5ad tables under {effects_dirs}")
    by_stem = {p.stem: p for p in tables}
    if target not in by_stem:
        raise FileNotFoundError(f"target {target} not among: {sorted(by_stem)[:12]} ...")
    return line_for_target(load_tables(tables, crispri_only), by_stem, target)


def load_tables(tables, crispri_only=False):
    """Every table as a Source (the expensive part: about 15 GB for all 117)."""
    print(f"loading {len(tables)} tables", flush=True)
    sources = []
    for p in tables:
        src = load_source(p)
        if crispri_only and not (src.modality == "single_cell" and src.reference == "non-targeting"):
            continue
        sources.append(src)
        print(f"  {src.name}: {src.lfc.shape[0]} knockdowns, modality={src.modality}, "
              f"reference={src.reference}, variance_x{src.mult:g}", flush=True)
    return sources


def line_for_target(sources, by_stem, target):
    """load_line's hold-out on already-loaded sources: several targets share one load."""
    target_src, others, same_line = split_sources(sources, target)
    if same_line:
        print(f"dropped {len(same_line)} sources of the target's cell line ({target_src.background}): "
              f"{', '.join(same_line)}", flush=True)
    if not others:
        raise RuntimeError("no source lines left after excluding the target")
    genes = target_src.genes
    basal_path = by_stem[target].with_name(by_stem[target].stem + "_baseline.csv")
    loaded = load_baseline(basal_path, genes)
    if loaded is None:
        raise FileNotFoundError(f"target baseline not found: {basal_path}")
    basal_cpm, tested = loaded
    # calibration uses each source's own baseline, so the curve sees that line's expression
    basals = {}
    for src in others:
        bp = by_stem[src.name].with_name(src.name + "_baseline.csv") if src.name in by_stem else None
        if bp is not None and bp.exists():
            got = load_baseline(bp, src.genes)
            if got is not None:
                basals[src.name] = (src.genes, got[0])
    return dict(by_stem=by_stem, sources=sources, target=target, target_src=target_src, others=others,
                same_line=same_line, genes=genes, basal_cpm=basal_cpm, tested=tested, basals=basals)


def load_line_for_controls(effects_dirs, background, genes, basal_cpm, tested, crispri_only=False):
    """load_line for a line known only from its control cells (no effects table of its own): basal_cpm and the
    tested gate come from those controls. Every table whose background is `background` is dropped, as
    split_sources drops the target's cell line."""
    tables = discover_tables(effects_dirs)
    by_stem = {p.stem: p for p in tables}
    print(f"loading {len(tables)} tables", flush=True)
    sources = []
    for p in tables:
        src = load_source(p)
        if crispri_only and not (src.modality == "single_cell" and src.reference == "non-targeting"):
            continue
        sources.append(src)
    line = _line_key(background)
    others = [s for s in sources if _line_key(s.background) != line]
    same_line = [s.name for s in sources if _line_key(s.background) == line]
    if same_line:
        print(f"dropped {len(same_line)} sources of the target's cell line ({background}): {', '.join(same_line)}",
              flush=True)
    basals = {}
    for src in others:
        bp = by_stem[src.name].with_name(src.name + "_baseline.csv")
        if bp.exists():
            got = load_baseline(bp, src.genes)
            if got is not None:
                basals[src.name] = (src.genes, got[0])
    return dict(by_stem=by_stem, sources=sources, target=None, target_src=None, others=others, same_line=same_line,
                genes=np.asarray(genes, dtype=object), basal_cpm=np.asarray(basal_cpm, dtype=np.float64),
                tested=np.asarray(tested, dtype=bool), basals=basals)


def arm_votes(others, arm):
    """(votes for prediction, votes for fit_prior). A: (None, None), every table row. B: the genetic tables
    (no drug table), one vote per background, for both. C: every table for prediction, genetic only for the prior."""
    if arm == "A":
        return None, None
    genetic = [s for s in others if s.modality != "drug"]
    prior_votes = BackgroundVotes(genetic)
    return (prior_votes if arm == "B" else BackgroundVotes(others)), prior_votes


def run_predict(effects_dirs, target, out_dir, crispri_only=False, knockdowns_file=None, arm="A", magnitude="slab"):
    return predict_line(load_line(effects_dirs, target, crispri_only), out_dir, knockdowns_file, arm, magnitude)


def predict_line(line, out_dir, knockdowns_file=None, arm="A", magnitude="slab"):
    """Fit and predict on a loaded line. arm A: every table row votes (vcc_2026_6). B: one vote per background,
    drug tables left out of pooling and of fit_prior. C: B plus the drug tables, one vote per background
    (each Tahoe line's drugs combined by combine_within); fit_prior on the genetic tables only, as in B.
    magnitude: "slab" (vcc_2026_6) or "gaussian" (no slab weight); the baseline since 2026-10-27 is arm C +
    gaussian."""
    by_stem, target, target_src, others = line["by_stem"], line["target"], line["target_src"], line["others"]
    same_line, genes, basal_cpm, tested, basals = (line["same_line"], line["genes"], line["basal_cpm"],
                                                   line["tested"], line["basals"])
    if arm not in ("A", "B", "C"):
        raise ValueError(f"arm must be A, B or C, got {arm!r}")
    votes, prior_votes = arm_votes(others, arm)
    print("fitting measurement-SE curve and transfer prior on the other lines", flush=True)
    calibration = fit_se_calibration(others, basals)
    prior = fit_prior(others) if arm == "A" else fit_prior(others, votes=prior_votes)
    print(f"  SE curve: log(SE) = {calibration.a:.3f} + {calibration.b:.3f}*log(cpm+1) "
          f"+ {calibration.c:.3f}*log(n)   n={calibration.n} r={calibration.r:.3f}", flush=True)
    print(f"  prior: tau2={prior.tau2:.4f} (from {prior.n_tau} genes)  "
          f"v0={prior.v0:.3f} (from {prior.n_v0} non-null genes)  pi={prior.pi:.3f}", flush=True)
    if knockdowns_file:
        want = [ln.strip() for ln in Path(knockdowns_file).read_text().splitlines() if ln.strip()]
    else:
        want = [str(k) for k in target_src.knockdowns]
    truth = truth_signature(by_stem[target], want, genes) if target in by_stem else None
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    K, G = len(want), len(genes)
    lfc = np.zeros((K, G), np.float32)
    lfc_raw = np.zeros((K, G), np.float32)
    se_m = np.zeros((K, G), np.float32)
    se_model = np.zeros((K, G), np.float32)
    p_sign = np.full((K, G), 0.5, np.float32)
    called = np.zeros((K, G), bool)
    n_src = np.zeros((K, G), np.int16)
    scores = []
    kd_scores = []
    kd_rows = []
    for i, k in enumerate(want):
        pred = predict_knockdown(others, k, genes, basal_cpm, tested, calibration, prior, drop_gene=k, votes=votes,
                                 magnitude=magnitude)
        lfc[i] = pred["lfc"]
        lfc_raw[i] = pred["lfc_unshrunk"]
        se_m[i] = pred["se_measurement"]
        # store sd, not variance, so the file is readable next to measurement SE
        vm = pred["v_model"].astype(np.float64)
        se_model[i] = np.where(np.isfinite(vm), np.sqrt(np.maximum(vm, 0.0)), np.nan).astype(np.float32)
        p_sign[i] = pred["p_sign"]
        called[i] = pred["de_sig"]
        n_src[i] = pred["n_sources"]
        rec = dict(knockdown=k, n_sources_median=float(np.median(pred["n_sources"])),
                   n_called=int(pred["de_sig"].sum()))
        if truth is not None and k in truth:
            y, e, sig = truth[k]
            sc = score_against_truth(pred, y, e, sig, tested)
            scores.append(sc)
            kd_scores.append(dict(knockdown=k, **sc))
            rec["n_true"] = sc.get("n_true")
            rec["nmae"] = sc.get("nmae")
        kd_rows.append(rec)
        if (i + 1) % 50 == 0 or i + 1 == K:
            print(f"  {i + 1}/{K} {k}: sources(median) {rec['n_sources_median']:.0f}, "
                  f"called {rec['n_called']}", flush=True)
    summary = dict(target=target, arm=arm, magnitude=magnitude, knockdowns=K, genes=G,
                   n_sources_counts="backgrounds" if votes is not None else "table rows",
                   backgrounds=None if votes is None else {bg: [s.name for s in ss] for bg, ss in votes.groups.items()},
                   sources=[dict(name=s.name, modality=s.modality, reference=s.reference, mult=s.mult,
                                 knockdowns=int(s.lfc.shape[0])) for s in others],
                   dropped_same_line=same_line,
                   calibration=calibration.as_dict(), prior=prior.as_dict(),
                   effect_space=summarize_scores(scores) if scores else None,
                   fit=_fit_stats(want, truth, lfc, lfc_raw, tested),
                   generator_lfc="shrunk posterior mean (lfc layer). lfc_unshrunk is the unshrunk pool.",
                   note=("Control quantiles are not in this file. Pass the target line's "
                         "v5 *_control_quantiles.npz plus these lfc and de_sig arrays to the frozen generator."))
    np.savez_compressed(
        out_dir / "predictions.npz",
        knockdowns=np.array(want, dtype=object), genes=np.array(list(map(str, genes)), dtype=object),
        lfc=lfc, lfc_unshrunk=lfc_raw, se_measurement=se_m, se_model=se_model,
        p_sign=p_sign, de_sig=called, n_sources=n_src,
    )
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")
    _write_kd_csv(out_dir / "per_knockdown.csv", kd_rows)
    (out_dir / "per_knockdown_scores.json").write_text(json.dumps(kd_scores, default=str) + "\n")
    print(json.dumps({k: summary[k] for k in ("target", "knockdowns", "calibration", "prior", "effect_space")},
                     indent=1, default=str), flush=True)
    return summary


def _fit_stats(want, truth, lfc, lfc_raw, tested):
    """summary["fit"]: calibration_stats of the shrunk and the unshrunk prediction (None without truth)."""
    if not truth:
        return None
    rows = {i: truth[str(k)][0] for i, k in enumerate(want) if str(k) in truth}
    return dict(shrunk=calibration_stats(lfc, rows, tested), unshrunk=calibration_stats(lfc_raw, rows, tested))


def _write_kd_csv(path, rows):
    keys = ["knockdown", "n_sources_median", "n_called", "n_true", "nmae"]
    lines = [",".join(keys)]
    for r in rows:
        lines.append(",".join("" if r.get(k) is None else str(r.get(k)) for k in keys))
    Path(path).write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------- generator handoff

def prepare_generator_inputs(lfc, de_sig, n_sources, p_sign):
    """What the frozen generator sees for one knockdown (steps 5 and 6).

    lfc is the shrunk posterior mean, de_sig the call set. A gene no source measured
    (n_sources == 0) is no change: LFC 0 and not called, so the generator draws it from the
    target line's own control quantiles. reach_order (gene indices by P(sign correct), most
    confident first) is saved only: the generator has no ranking input.
    """
    unmeasured = np.asarray(n_sources) <= 0
    out_lfc = np.where(unmeasured, 0.0, np.asarray(lfc, dtype=np.float64)).astype(np.float32)
    out_sig = np.asarray(de_sig, dtype=bool) & ~unmeasured
    reach = np.argsort(-np.asarray(p_sign, dtype=np.float64), kind="stable").astype(np.int32)
    return dict(lfc=out_lfc, de_sig=out_sig, reach_order=reach,
                n_unmeasured=int(unmeasured.sum()), n_called=int(out_sig.sum()))


def knockdown_seed(name, salt=""):
    """Generator seed from the knockdown name (not its row), so a rerun draws the same cells
    and two knockdowns never share one control draw. `salt` gives an independent draw
    (sha256(salt + name)); the default "" reproduces the unsalted seeds."""
    import hashlib
    return int.from_bytes(hashlib.sha256((str(salt) + str(name)).encode()).digest()[:4], "little")


# ---------------------------------------------------------------- synthetic self-test

def _synthetic(seed=7):
    rng = np.random.default_rng(seed)
    G, K = 350, 30
    names = ["L1", "L2", "L3", "L4"]
    genes = np.array([f"g{i}" for i in range(G)], dtype=object)
    kds = np.array([f"k{i}" for i in range(K)], dtype=object)
    m0 = np.exp(rng.normal(2.2, 1.1, size=(len(names), G)))
    shared = rng.normal(0, 0.7, size=(K, G)) * (rng.random((K, G)) < 0.10)
    sources, truth = [], {}
    for li, name in enumerate(names):
        shift = rng.normal(0, 0.15, size=(K, G))
        true = shared + shift
        n_cells = rng.integers(60, 380, size=K).astype(np.float32)
        se = (0.55 / np.sqrt(m0[li] + 1.0)) * np.sqrt(400.0 / n_cells)[:, None]
        se = se * rng.uniform(0.85, 1.15, size=(K, G))
        y = true + rng.normal(0, 1, size=(K, G)) * se
        if name == "L4":
            y[:, 280:] = np.nan
            se[:, 280:] = np.nan
        truth[name] = true
        sources.append(Source(name, name, "single_cell", "non-targeting", genes, kds,
                              y.astype(np.float32), se.astype(np.float32), n_cells, mult=1.0))
    return sources, truth, m0, genes, kds


def self_test():
    sources, truth, m0, genes, kds = _synthetic()
    rng_t = np.random.default_rng(11)
    target, others = sources[-1], sources[:-1]
    basals = {s.name: (genes, m0[i]) for i, s in enumerate(sources)}
    cal = fit_se_calibration(others, basals, max_rows=30, max_genes=350, seed=1)
    prior = fit_prior(others, max_knockdowns=30, seed=1)
    # expression slope should be negative, depth slope negative
    assert cal.b < 0, cal.as_dict()
    assert cal.c < 0, cal.as_dict()
    tested = m0[-1] > 5.0
    rows = []
    copy_abs, ours_abs = [], []
    for i, k in enumerate(kds):
        pred = predict_knockdown(others, k, genes, m0[-1], tested, cal, prior, drop_gene=None)
        # reference calls from the true effect at the target measurement SE
        se_ref = cal.predict(m0[-1], N_CELLS_REF)
        z = np.abs(truth["L4"][i]) / np.maximum(se_ref, SE_MIN)
        p = np.full(z.shape, np.nan)
        p[tested] = np.minimum(1.0, 2.0 * norm_sf(z[tested]))
        sig = bh_reject(p)
        sc = score_against_truth(pred, truth["L4"][i], se_ref, sig, tested)
        rows.append(sc)
        ok = np.isfinite(truth["L4"][i]) & tested
        copy_abs.append(np.abs(others[0].lfc[i].astype(np.float64) - truth["L4"][i])[ok].mean())
        ours_abs.append(sc["mae"])
        # no measured source gene that is all-NaN should be called off the panel expression gate
        assert not pred["de_sig"][~tested].any()
    summary = summarize_scores(rows)
    mae_copy = float(np.mean(copy_abs))
    mae_ours = float(np.mean(ours_abs))
    print(json.dumps(dict(calibration=cal.as_dict(), prior=prior.as_dict(),
                          effect_space=summary, mae_single_source_copy=mae_copy, mae_model=mae_ours),
                     indent=1), flush=True)
    assert mae_ours < mae_copy, (mae_ours, mae_copy)
    assert summary["nmae"] < summary["nmae_zero"]
    assert summary["coverage_95"] is not None and 0.85 <= summary["coverage_95"] <= 0.99
    assert summary["auroc"] is not None and summary["auroc"] > 0.75
    # disagreement shrinks toward 0 and does not call; agreement does
    g = np.array(["x"], dtype=object)
    hi = Source("A", "A", "single_cell", "non-targeting", g, np.array(["k"]),
                np.array([[3.0]], np.float32), np.array([[0.05]], np.float32),
                np.array([400], np.float32), 1.0)
    lo = Source("B", "B", "single_cell", "non-targeting", g, np.array(["k"]),
                np.array([[-3.0]], np.float32), np.array([[0.05]], np.float32),
                np.array([400], np.float32), 1.0)
    agree_a = Source("A", "A", "single_cell", "non-targeting", g, np.array(["k"]),
                     np.array([[1.70]], np.float32), np.array([[0.05]], np.float32),
                     np.array([400], np.float32), 1.0)
    agree_b = Source("B", "B", "single_cell", "non-targeting", g, np.array(["k"]),
                     np.array([[1.55]], np.float32), np.array([[0.05]], np.float32),
                     np.array([400], np.float32), 1.0)
    basal = np.array([50.0])
    gate = np.array([True])
    tiny_prior = Prior(tau2=prior.tau2, v0=prior.v0, pi=prior.pi, n_tau=1, n_v0=1)
    clash = predict_knockdown([hi, lo], "k", g, basal, gate, cal, tiny_prior)
    same = predict_knockdown([agree_a, agree_b], "k", g, basal, gate, cal, tiny_prior)
    assert abs(float(clash["lfc"][0])) < 0.5, clash["lfc"]
    assert not bool(clash["de_sig"][0])
    assert bool(same["de_sig"][0]) and float(same["lfc"][0]) > 0.5, same["lfc"]
    # empty sources -> no call, zero effect
    empty = predict_knockdown([], "missing", g, basal, gate, cal, tiny_prior)
    assert float(empty["lfc"][0]) == 0.0 and not bool(empty["de_sig"][0])
    # drug table: two drugs on one target stay two observations; an agonist and an antagonist
    # do not average to a precise value but shrink to no call
    drugs = Source("T", "T", "drug", "non-targeting", g, np.array(["k", "k"]),
                   np.array([[3.0], [-3.0]], np.float32), np.array([[0.05], [0.05]], np.float32),
                   np.array([400, 400], np.float32), multiplier_for("drug", "non-targeting"))
    assert drugs.rows_of["k"] == [0, 1]
    two = predict_knockdown([drugs], "k", g, basal, gate, cal, tiny_prior)
    assert int(two["n_sources"][0]) == 2, two["n_sources"]
    assert abs(float(two["lfc"][0])) < 0.5 and not bool(two["de_sig"][0]), two["lfc"]
    # same cell line as the target does not vote; other lines do
    mk = lambda name, bg: Source(name, bg, "single_cell", "non-targeting", g, np.array(["k"]),
                                 np.array([[1.0]], np.float32), np.array([[0.1]], np.float32),
                                 np.array([400], np.float32), 1.0)
    tgt, voters, dropped = split_sources(
        [mk("replogle_k562_gw", "K562"), mk("jiang_ifnb_k562", "K562"), mk("replogle_rpe1", "RPE1")],
        "replogle_k562_gw")
    assert tgt.name == "replogle_k562_gw" and [s.name for s in voters] == ["replogle_rpe1"], voters
    assert dropped == ["jiang_ifnb_k562"], dropped
    # Step 2, within-background rule: a second identical table of the same background adds no precision
    one_row = Source("K1", "K562", "single_cell", "non-targeting", g, np.array(["k"]),
                     np.array([[1.2]], np.float32), np.array([[0.2]], np.float32), np.array([400], np.float32), 1.0)
    twin = Source("K2", "K562", "single_cell", "non-targeting", g, np.array(["k"]),
                  np.array([[1.2]], np.float32), np.array([[0.2]], np.float32), np.array([400], np.float32), 1.0)
    Y1, S1, _ = BackgroundVotes([one_row]).gather("k", g)
    Y2, S2, _ = BackgroundVotes([one_row, twin]).gather("k", g)
    assert len(Y2) == 1 and abs(S2[0][0] - S1[0][0]) < 1e-12 and abs(Y2[0][0] - 1.2) < 1e-6, (S1, S2)
    # an agonist and an antagonist on one target in one Tahoe line: large variance, no call
    vt = BackgroundVotes([drugs])
    Yt, St, _ = vt.gather("k", g)
    assert St[0][0] > 2.0, St
    vote = predict_knockdown([drugs], "k", g, basal, gate, cal, tiny_prior, votes=vt)
    assert not bool(vote["de_sig"][0]) and abs(float(vote["lfc"][0])) < 0.5 and int(vote["n_sources"][0]) == 1, vote
    # arm B: a drug table changes neither the predictions nor the prior
    drug_tab = Source("tahoe_x", "T", "drug", "non-targeting", genes, kds,
                      rng_t.normal(0, 2, size=(len(kds), len(genes))).astype(np.float32),
                      np.full((len(kds), len(genes)), 0.1, np.float32), np.full(len(kds), 400, np.float32),
                      multiplier_for("drug", "non-targeting"))
    outs = []
    for pool in (others, others + [drug_tab]):
        v_pred, v_prior = arm_votes(pool, "B")
        pr = fit_prior(pool, max_knockdowns=30, seed=1, votes=v_prior)
        p0 = predict_knockdown(pool, kds[0], genes, m0[-1], tested, cal, pr, votes=v_pred)
        outs.append((pr.as_dict(), p0["lfc"], p0["de_sig"]))
    assert outs[0][0] == outs[1][0], (outs[0][0], outs[1][0])
    assert np.array_equal(outs[0][1], outs[1][1]) and np.array_equal(outs[0][2], outs[1][2])
    # magnitude "gaussian": mu x v0 / (v0 + v_model), the slab weight not applied; the agreeing pair keeps its size
    gm = predict_knockdown([agree_a, agree_b], "k", g, basal, gate, cal, tiny_prior, magnitude="gaussian")
    sl = predict_knockdown([agree_a, agree_b], "k", g, basal, gate, cal, tiny_prior)
    mu_ab, var_ab, tau_ab, _, _ = pool_effects(np.array([[1.70], [1.55]]), np.array([[0.05], [0.05]]), np.ones(2), tiny_prior.tau2)
    v_ab = var_ab + np.maximum(tau_ab, 0)
    assert abs(float(gm["lfc"][0]) - float(mu_ab[0] * tiny_prior.v0 / (tiny_prior.v0 + v_ab[0]))) < 1e-6, gm["lfc"]
    assert float(gm["lfc"][0]) >= float(sl["lfc"][0]) - 1e-9 and bool(gm["de_sig"][0])
    # background aliases: HT-29 is HT29, HepG2/C3A is HepG2; KOLF2.1J and kolf_2 stay apart
    assert _line_key("HT-29") == _line_key("ht29") and _line_key("HepG2/C3A") == _line_key("HepG2")
    assert _line_key("KOLF2.1J") != _line_key("kolf_2")
    # handoff: an unmeasured gene goes to the generator as no change, never called
    h = prepare_generator_inputs(np.array([1.5, 0.7, -2.0]), np.array([True, True, False]),
                                 np.array([3, 0, 2]), np.array([0.9, 0.99, 0.6]))
    assert h["lfc"].tolist() == [1.5, 0.0, -2.0] and h["de_sig"].tolist() == [True, False, False], h
    assert h["n_unmeasured"] == 1 and h["n_called"] == 1 and h["reach_order"].tolist() == [1, 0, 2], h
    # seed: fixed per name, different across names
    assert knockdown_seed("TP53") == knockdown_seed("TP53") != knockdown_seed("MYC")
    # salt: "" is the unsalted seed; another salt is another draw
    import hashlib
    assert knockdown_seed("TP53") == int.from_bytes(hashlib.sha256(b"TP53").digest()[:4], "little")
    assert knockdown_seed("TP53", "") == knockdown_seed("TP53") != knockdown_seed("TP53", "s1")
    # calibration_stats: a perfect prediction has r = 1 and both slopes 1; half-size has transmission 0.5
    y = np.array([[1.0, -2.0, 0.5, np.nan]])
    cs = calibration_stats(y, {0: y[0]}, np.array([True, True, True, True]))
    assert cs["n"] == 3 and abs(cs["pearson_r"] - 1) < 1e-12 and abs(cs["slope_true_on_pred"] - 1) < 1e-12, cs
    half = calibration_stats(0.5 * y, {0: y[0]}, np.ones(4, bool))
    assert abs(half["slope_pred_on_true"] - 0.5) < 1e-12 and abs(half["slope_true_on_pred"] - 2) < 1e-12, half
    assert abs(half["beta"] - 2) < 1e-12 and abs(half["r2_fit"] - 1) < 1e-12 and half["r2_raw"] < 1, half
    assert abs(cs["r2_raw"] - 1) < 1e-12 and abs(cs["beta"] - 1) < 1e-12, cs
    # BH: only the first of three p-values rejects at 0.05
    mask = bh_reject(np.array([0.001, 0.04, 0.5]))
    assert mask.tolist() == [True, False, False]
    print("self-test passed", flush=True)


if __name__ == "__main__":
    self_test()

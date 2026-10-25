"""Per-source effect tables v5 (lib/effects_v5/): the quantities cell-eval2 scores (vcc2026 profile, rule_version 3),
in the regime the scorer sees. Never touches lib/effects_v4/ or the QC outputs.

Scorer definitions (cell_eval2 0.18 de_compute.compute_de + configs/vcc2026.yaml; checked against pdex 0.3 / numba_mwu
and gpudge 0.9.1, which the competition path uses -- all three run the same test):
  CPM      per cell, counts / (library / 1e6), library = the cell's sum over the genes in the matrix (scanpy
           normalize_total, target_sum 1e6). Here: the source's genes that map to the 18,533-gene panel.
  LFC      log2((mean CPM of the knockdown's cells + 1e-9) / (mean CPM of the control cells + 1e-9)), arithmetic means.
  test     two-sided Mann-Whitney U, normal approximation, tie-corrected, continuity-corrected (scipy asymptotic
           defaults), on log1p(CPM) (rank-identical to CPM; it is pdex's input, so the tie structure is the scorer's).
           A gene whose combined sample is constant gets p = 1.
  gate     a gene is tested only when the CONTROL mean CPM > 5 (strict; rule_version 3: control alone decides).
  FDR      Benjamini-Hochberg per knockdown over the gated genes; significant = p_adj < 0.05. The target's own gene is
           in the BH universe (as in the scorer) and is then removed from the significant set (as the direction metrics do).
  profile  log1p(5e4 * summed counts / summed total) per group (bulk_target_sum 5e4, natural log).

Scorer regime (the leaderboard's held-out contexts A-C: 18,400 control cells each, median 20k UMI, q10-q90 ~9k-35k):
  1. depth: every cell of the source is binomially thinned by ONE factor f = min(1, 20,000 / median control n_umi),
     drawn per cell from default_rng([2026, row in the QC'd file]) so any cell can be re-drawn independently. A single
     factor moves the source's median to 20k and keeps its depth spread, as the contexts have one; sources already
     below 20k (all but H1) are left as they are (no upsampling: their DE power stays below the scorer's).
  2. cells: per kept knockdown, a random 400 of its kept cells (all when fewer); controls: a random 18,400 (all when
     fewer). Seeded (2026).
  LFC, SE, the gate, the test and BH are all computed on these sampled, thinned cells.
  SE       delta method on the sampled cells: sqrt(var1 / n1 / (m1 + eps)^2 + var0 / n0 / (m0 + eps)^2) / ln 2, sample
           variances of per-cell CPM. +inf where either mean is 0 (the LFC sits on the epsilon floor, about +-30, and
           has no finite SE). NaN where the source did not measure the gene.
  profiles two per knockdown (and for the controls), log1p(5e4 * summed counts / summed total) over the panel genes:
           log1p_bulk_50k_matched (PRIMARY): the same sampled, thinned cells as LFC and the test, so it carries the
             sampling variance of the single ~400-cell / ~20k-UMI draw the scorer sees;
           log1p_bulk_50k_full: all kept cells (all controls), raw counts, not thinned (the best estimate of the
             expected profile).
  library  CPM, the gate and both profiles use ONLY the source's genes that map onto the 18,533-gene panel; genes off
           the panel never enter a library or an output layer (the build reports how far the panel-only library is
           from the whole-matrix one). The one exception is obs efficiency, the target's own LFC, which is read in the
           source's genes (with the panel library) when the target is off the panel.

Control quantiles (sidecar <source>_control_quantiles.npz, cell sources): per gated panel gene, the sorted unique
nonzero log1p(CPM) values of the sampled control cells (float64, bit-identical to what the scorer ranks), their
multiplicities and the zero count. With these, a predicted knockdown (cell values, or a first moment turned into cell
values) can be rank-tested against the controls exactly as the scorer does, without re-reading cells. de_sig / n_real
stay the ready-made labels; the SE layer is a modelling aid, not the significance source.

Cell sources (K562, RPE1, HepG2, Jurkat, H1, and the Orion HCT116 / HEK293T subsets we hold: the 300 A-C targets plus
controls): lib/qc_cells/<source>_qc.h5ad (already cell-QC'd), with the existing perturbation QC: guide_qc.csv (kept
constructs) and knockdown_qc.csv (kept knockdowns); Orion has no constructs, so its knockdown QC is the table QC
(results/qc/orion_<line>_table/knockdown_qc.csv), restricted to the targets in the subset.
Marson CD4 T (Rest, Stim8hr, Stim48hr): replicate-level tables only (lib/tcell/). No single cells, so NO Wilcoxon
labels: de_sig is empty, n_real is NaN, obs de_labels_known = False, modality = 'replicate_table',
significance_unknown = True (and the same in uns), so training can treat them like the Tahoe modality flag. Their LFC
stays the table's replicate-paired log2 ratio with pseudocount 0.5: the per-knockdown means needed for epsilon = 1e-9
are not in the tables we hold (documented gap; lfc_type says so). No per-cell ranks are invented. Cell sources carry
modality = 'single_cell', significance_unknown = False.

Output lib/effects_v5/<source>.h5ad: rows = kept knockdowns, id '<source>|<background>|<condition>|<target_gene_id>';
obs source, background, condition, target_gene, target_gene_id, n_cells_kept, n_cells_used, efficiency (official LFC
of the target, NaN unless measurable), efficiency_known, n_real, de_labels_known, modality, significance_unknown,
lfc_type, se_type; var = the panel; layers 'lfc', 'se', 'log1p_bulk_50k_matched', 'log1p_bulk_50k_full' (float32, NaN
unmeasured), 'de_sig' (sparse bool). Plus <source>_baseline.csv (per panel gene: control mean CPM and detection rate
in the sampled regime, the full-depth control mean CPM, both control profiles, de_tested = the gate),
<source>_control_quantiles.npz (cell sources) and <source>_build.json (summary, checks, notes).

Validation (cell sources). The cells are re-read from the file and re-thinned independently of the build's batches,
then (a) cell_eval2's own compute_de (pdex backend, vcc2026 parameters) and (b) scipy.stats.mannwhitneyu +
false_discovery_control are run on them, and (c) a plain numpy rank sum against the stored control quantiles (read back
from the sidecar) re-derives de_sig.
  fast (every source): 5 random knockdowns x 20 random measured genes.
  --full-validate N: N random knockdowns (all when fewer) x every measured gene; reports the max |LFC| difference, the
    knockdowns with a nonzero significant-set symmetric difference and the total genes in it, for (a), (b) and (c).
Checks: |LFC - stored| < 1e-5, symmetric differences 0, epsilon 1e-9, n_cells_used <= 400, the panel-library check,
and the count check (per-cell integer counts, one row per cell: Orion's files are cells, not sums).

Usage: ~/venv/bin/python scripts/build_effects_v5.py --source K562|RPE1|HepG2|Jurkat|H1|HCT116|HEK293T|CD4_Rest|CD4_Stim8hr|CD4_Stim48hr
       [--full-validate 30]
"""
import argparse
import json
import math
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numba as nb
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem
from scipy.stats import false_discovery_control

sys.path.insert(0, str(Path(__file__).parent))
from build_effects_v4 import panel                                                        # noqa: E402
from qc_replogle_umi import SOURCES                                                        # noqa: E402

V = Path.home() / 'vcc'
OUT = V / 'lib/effects_v5'
EPS, MIN_CPM, FDR, BULK_TS = 1e-9, 5.0, 0.05, 5e4
DEPTH, N_CELLS, N_CTRL, SEED = 20000, 400, 18400, 2026
CHUNK = 50000
CELL_SOURCES = {'K562': 'replogle_k562_gw', 'RPE1': 'replogle_rpe1', 'HepG2': 'nadig_hepg2', 'Jurkat': 'nadig_jurkat',
                'H1': 'arc_h1', 'HCT116': 'orion_hct116', 'HEK293T': 'orion_hek293t'}
ORION = ('HCT116', 'HEK293T')
TABLES = {'CD4_Rest': ('marson_cd4_rest', 'CD4_T', 'Rest'), 'CD4_Stim8hr': ('marson_cd4_stim8hr', 'CD4_T', 'Stim8hr'),
          'CD4_Stim48hr': ('marson_cd4_stim48hr', 'CD4_T', 'Stim48hr')}
LFC_TYPE, SE_TYPE = 'official_cpm_1e-9', 'downsampled_400_20k_delta'
LIB_TOL = 0.01                     # panel-only vs whole-matrix library: documented when the medians differ by more


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


# ---------------------------------------------------------------- sampling regime
def thin_row(counts, row, f):
    """The one depth draw for a cell: binomial thinning by f, seeded by the cell's row in the QC'd file."""
    if f >= 1.0:
        return counts
    return np.random.default_rng([SEED, int(row)]).binomial(counts, f).astype(counts.dtype)


def read_rows(X, rows, f):
    """Rows (sorted) of an h5 CSR group, thinned: (data, indices, indptr). Reads only the slices it needs."""
    ip = X['indptr']; data, idx, ptr = [], [], [0]
    for r in rows:
        a, b = int(ip[r]), int(ip[r + 1])
        d = thin_row(X['data'][a:b], r, f); i = X['indices'][a:b]; k = d > 0
        data.append(d[k]); idx.append(i[k]); ptr.append(ptr[-1] + int(k.sum()))
    return np.concatenate(data), np.concatenate(idx), np.array(ptr)


# ---------------------------------------------------------------- the scorer's test, exactly
@nb.njit(cache=True)
def _ndtr(x):
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


@nb.njit(parallel=True, cache=True)
def mwu_vs_control(PT, cptr, cvals, cnz, ctie, n0):
    """Two-sided MWU p per gene: knockdown cells (rows of PT: genes x n1, log1p CPM) vs the controls (per gene: sorted
    nonzero values cvals[cptr[g]:cptr[g+1]], zero count cnz[g], tie term sum(t^3 - t) over the nonzero control values
    ctie[g]). U1 = sum over knockdown cells of (#controls below + 0.5 #controls equal): the cross-pair count, identical to
    R1 - n1(n1+1)/2 with mid-ranks. Variance, tie and continuity terms as numba_mwu._sparse_mwu_column (pdex)."""
    G, n1 = PT.shape
    out = np.ones(G)
    for g in nb.prange(G):
        v = np.sort(PT[g])
        c = cvals[cptr[g]:cptr[g + 1]]
        nz0 = cnz[g]
        z1 = 0
        while z1 < n1 and v[z1] == 0.0:
            z1 += 1
        U1 = 0.5 * nz0 * z1
        tie = ctie[g]
        i = z1
        while i < n1:
            j = i
            while j + 1 < n1 and v[j + 1] == v[i]:
                j += 1
            p = j - i + 1
            lo = np.searchsorted(c, v[i], 'left'); hi = np.searchsorted(c, v[i], 'right'); ce = hi - lo
            U1 += p * (nz0 + lo + 0.5 * ce)
            t = ce + p
            tie += (t * t * t - t) - (ce * ce * ce - ce)
            i = j + 1
        zt = float(nz0 + z1)
        tie += zt * zt * zt - zt
        n = n0 + n1
        U2 = n1 * n0 - U1
        U = max(U1, U2)
        s_sq = n1 * n0 / 12.0 * ((n + 1.0) - tie / (float(n) * (n - 1.0)))
        if s_sq <= 0.0:
            out[g] = 1.0
            continue
        pv = 2.0 * _ndtr(-(U - n1 * n0 / 2.0 - 0.5) / math.sqrt(s_sq))
        out[g] = min(max(pv, 0.0), 1.0)
    return out


def control_index(C, cols):
    """Per tested gene (columns `cols` of the control CPM CSR): sorted nonzero log1p values, zero count, tie term."""
    Cc = C[:, cols].tocsc(); Cc.sort_indices()
    vals = np.log1p(Cc.data.astype(np.float64)); ptr = Cc.indptr.astype(np.int64)
    for g in range(len(cols)):
        vals[ptr[g]:ptr[g + 1]].sort()
    nnz = np.diff(ptr); cnz = (C.shape[0] - nnz).astype(np.int64)
    tie = np.zeros(len(cols))
    for g in range(len(cols)):
        s = vals[ptr[g]:ptr[g + 1]]
        if len(s):
            _, t = np.unique(s, return_counts=True); t = t.astype(float); tie[g] = (t ** 3 - t).sum()
    return ptr, vals, cnz, tie


def save_quantiles(source, tested, to_panel, ptr, cvals, cnz, n0):
    """Sidecar: per gated gene, the sorted unique nonzero log1p(CPM) control values, their multiplicities, the zeros."""
    uptr, uv, uc = [0], [], []
    for g in range(len(tested)):
        v, c = np.unique(cvals[ptr[g]:ptr[g + 1]], return_counts=True)
        uv.append(v); uc.append(c.astype(np.int32)); uptr.append(uptr[-1] + len(v))
    path = OUT / f'{source}_control_quantiles.npz'
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, panel_col=to_panel[tested].astype(np.int32), ptr=np.array(uptr, np.int64),
                        values=np.concatenate(uv) if uv else np.zeros(0), counts=np.concatenate(uc) if uc else np.zeros(0, np.int32),
                        zeros=cnz.astype(np.int32), n_controls=np.int64(n0),
                        description=np.array('per gated panel gene g (panel_col[g]): values[ptr[g]:ptr[g+1]] = sorted unique '
                                             'nonzero log1p(CPM) of the sampled controls (float64), counts = their '
                                             'multiplicities, zeros[g] = controls at 0; n_controls in all'))
    return path, int(uptr[-1])


def quantile_p(x, vals, cnts, nz0, n0):
    """Rank-sum p for one gene from the stored control quantiles alone (plain numpy; independent of the numba kernel):
    two-sided MWU, normal approximation, tie + continuity corrected, as the scorer. A gated gene always has nonzero
    controls (mean CPM > 5), so vals is never empty."""
    n1 = len(x); xs, p = np.unique(x[x > 0], return_counts=True); z1 = n1 - int(p.sum())
    below = np.concatenate([[0], np.cumsum(cnts)])
    i = np.searchsorted(vals, xs); j = np.minimum(i, len(vals) - 1)
    ce = np.where((i < len(vals)) & (vals[j] == xs), cnts[j], 0).astype(float)
    U1 = 0.5 * nz0 * z1 + float((p * (nz0 + below[i] + 0.5 * ce)).sum())
    c = cnts.astype(float); t = ce + p
    tie = float((c ** 3 - c).sum() + ((t ** 3 - t) - (ce ** 3 - ce)).sum()) + (nz0 + z1) ** 3 - (nz0 + z1)
    n = n0 + n1; s_sq = n1 * n0 / 12.0 * ((n + 1.0) - tie / (n * (n - 1.0)))
    if s_sq <= 0: return 1.0
    z = (max(U1, n1 * n0 - U1) - n1 * n0 / 2.0 - 0.5) / math.sqrt(s_sq)
    return float(min(max(math.erfc(z / math.sqrt(2.0)), 0.0), 1.0))       # 2 * ndtr(-z) = erfc(z / sqrt 2)


def bh_sig(p):
    """Indices with BH p_adj < FDR (scipy's false_discovery_control, as the scorer)."""
    return np.flatnonzero(false_discovery_control(p, method='bh') < FDR) if len(p) else np.zeros(0, int)


def cpm_csr(data, indices, indptr, n_vars, panel_mask):
    """CSR of per-cell CPM, library = the row's sum over the panel-mapped genes. The same float operation as scanpy's
    normalize_total (x / (lib / 1e6)), so values, and so ties, are bit-identical to what the scorer ranks."""
    data = data.astype(np.float64); n = len(indptr) - 1
    lib = np.bincount(np.repeat(np.arange(n), np.diff(indptr)), weights=data * panel_mask[indices], minlength=n)
    factor = lib / 1e6
    out = data / np.repeat(np.where(factor > 0, factor, 1.0), np.diff(indptr))
    return sp.csr_matrix((out, indices, indptr), shape=(n, n_vars)), lib


# ---------------------------------------------------------------- cell sources
def source_genes(dataset, var_src, var):
    """Map the source's genes to panel columns (-1 = off the panel or not measured)."""
    if dataset in ORION:                                   # Orion subset: var = the panel symbols, 'measured' mask
        col = {s: i for i, s in enumerate(var.index)}
        to_panel = np.array([col.get(s, -1) for s in var_src.index.astype(str)])
        to_panel[~var_src['measured'].to_numpy(bool)] = -1
        ids = var.gene_id.to_numpy()[np.maximum(to_panel, 0)]
        ids = np.where(to_panel >= 0, ids, var_src.index.astype(str))
        return to_panel, np.array([i if i else s for i, s in zip(ids, var_src.index.astype(str))])
    ids = var_src.index.astype(str).to_numpy() if 'gene_name' in var_src else var_src['gene_id'].astype(str).to_numpy()
    pid = {g: i for i, g in enumerate(var.gene_id) if g}
    return np.array([pid.get(g, -1) for g in ids]), ids


def kept_cells(dataset, obs):
    """Kept knockdowns (DataFrame) and, per QC'd cell, its kept-knockdown row (-1 = not kept), from the existing QC."""
    _, qc_dir, _, _ = SOURCES[dataset]
    tid = obs['target_gene_id'].astype(str).to_numpy(); ctrl = obs['is_control'].to_numpy(bool)
    if dataset in ORION:
        kq = V / f'results/qc/orion_{dataset.lower()}_table/knockdown_qc.csv'
        qc = pd.read_csv(kq); qc['target_gene'] = qc['target_gene'].astype(str)
        present = set(obs.loc[~ctrl, 'target_gene'].astype(str))
        kept = qc[qc.knockdown_keep.astype(bool) & qc.target_gene.isin(present)].reset_index(drop=True)
        sym2id = dict(zip(obs['target_gene'].astype(str), tid))
        kept['target_gene_id'] = [sym2id[t] for t in kept.target_gene]
        row_of = {t: i for i, t in enumerate(kept.target_gene_id)}
        kd = np.array([-1 if c else row_of.get(t, -1) for c, t in zip(ctrl, tid)])
        return kept, kd, dict(knockdown_qc=str(kq), guide_qc='none (Orion: no constructs; table knockdown QC)')
    gq, kq = qc_dir / 'guide_qc.csv', qc_dir / 'knockdown_qc.csv'
    guides = pd.read_csv(gq); kds = pd.read_csv(kq)
    kds['target_gene_id'] = kds['target_gene_id'].astype(str); kds['target_gene'] = kds['target_gene'].astype(str)
    kept = kds[kds.knockdown_keep.astype(bool)].reset_index(drop=True)
    # cell -> kept construct -> kept knockdown, through guide_qc (as v4 pooled them; read-through targets such as
    # NEDD8-MDP1 have no Ensembl ID, so the obs target_gene_id cannot be the key)
    row_of = {t: i for i, t in enumerate(kept.target_gene_id)}
    g = guides[guides.keep.astype(bool)]
    guide_kd = {gid: row_of.get(t, -1) for gid, t in zip(g.guide_id.astype(str), g.target_gene_id.astype(str))}
    guide = obs['guide_id'].astype(str).to_numpy()
    kd = np.array([-1 if c else guide_kd.get(gg, -1) for c, gg in zip(ctrl, guide)])
    return kept, kd, dict(knockdown_qc=str(kq), guide_qc=str(gq))


def sample(kd, ctrl, n_kd, rng):
    """Sampled rows: per knockdown up to N_CELLS of its kept cells, and up to N_CTRL controls (sorted row lists)."""
    order = np.argsort(kd, kind='stable'); ks = kd[order]
    starts = np.searchsorted(ks, np.arange(n_kd)); ends = np.searchsorted(ks, np.arange(n_kd), 'right')
    picked = []
    for k in range(n_kd):
        r = order[starts[k]:ends[k]]
        picked.append(np.sort(rng.choice(r, N_CELLS, replace=False)) if len(r) > N_CELLS else np.sort(r))
    c = np.flatnonzero(ctrl)
    c = np.sort(rng.choice(c, N_CTRL, replace=False)) if len(c) > N_CTRL else c
    return picked, c, (ends - starts)


def stream(path, f, want, accumulate=None):
    """One pass over the file: thinned CSR (uint32 counts, uint16 indices) of the rows with want[row] True, in row
    order. Buffers are preallocated at the unthinned size and trimmed, so memory never holds two copies."""
    with h5py.File(path, 'r') as fh:
        X = fh['X']; ip = X['indptr'][:]; n_obs = len(ip) - 1; n_vars = int(X.attrs['shape'][1])
        rows_out = np.flatnonzero(want); ub = int(np.diff(ip)[rows_out].sum())
        data = np.empty(ub, np.uint32); indices = np.empty(ub, np.uint16); indptr = np.zeros(len(rows_out) + 1, np.int64)
        t0 = time.time(); w = 0; ri = 0
        for s in range(0, n_obs, CHUNK):
            e = min(s + CHUNK, n_obs); a, b = int(ip[s]), int(ip[e])
            D = X['data'][a:b]; I = X['indices'][a:b]; P = ip[s:e + 1] - a
            if accumulate is not None:
                accumulate(s, e, D, I, P, n_vars)
            for r in np.flatnonzero(want[s:e]) + s:
                lo, hi = int(ip[r] - a), int(ip[r + 1] - a)
                d = thin_row(D[lo:hi], r, f); k = d > 0; m = int(k.sum())
                data[w:w + m] = d[k]; indices[w:w + m] = I[lo:hi][k]; w += m
                ri += 1; indptr[ri] = w
            if (s // CHUNK) % 10 == 0 or e == n_obs:
                log(f'  pass: {e:,} / {n_obs:,} cells ({time.time() - t0:.0f}s)')
    return rows_out, data[:w], indices[:w], indptr, n_vars


def build_cells(dataset, batch_nnz):
    source = CELL_SOURCES[dataset]; _, qc_dir, path, title = SOURCES[dataset]
    background = dataset
    with h5py.File(path, 'r') as fh:
        obs = read_elem(fh['obs']); var_src = read_elem(fh['var']); row_nnz = np.diff(fh['X']['indptr'][:])
    var, _ = panel()
    to_panel, src_ids = source_genes(dataset, var_src, var); on = to_panel >= 0
    assert len(np.unique(to_panel[on])) == int(on.sum()), "duplicate panel mapping"
    assert len(var_src) < 65536, f"gene index overflow: {len(var_src)}"      # stream() holds gene indices as uint16
    kept, kd, qc_files = kept_cells(dataset, obs)
    ctrl = obs['is_control'].to_numpy(bool); n_umi = obs['n_umi'].to_numpy(float)
    f = min(1.0, DEPTH / float(np.median(n_umi[ctrl])))
    rng = np.random.default_rng(SEED)
    picked, ctrl_rows, n_kept = sample(kd, ctrl, len(kept), rng)
    if dataset not in ORION and not (n_kept == kept.n_cells_kept.to_numpy()).all():
        sys.exit(f'error: kept cells per knockdown differ from knockdown_qc n_cells_kept for '
                 f'{int((n_kept != kept.n_cells_kept.to_numpy()).sum())} knockdowns')
    if (n_kept == 0).any():
        sys.exit(f'error: {int((n_kept == 0).sum())} kept knockdowns have no cells in {path}')
    log(f'{source}: {len(kept):,} kept knockdowns, {int(n_kept.sum()):,} kept cells, {ctrl.sum():,} controls '
        f'({len(ctrl_rows):,} sampled); median control n_umi {np.median(n_umi[ctrl]):.0f} -> thinning factor {f:.4f}; '
        f'{int(on.sum()):,} of {len(on):,} source genes on the panel')
    # batches of knockdowns by sampled nnz
    kd_nnz = np.array([row_nnz[r].sum() * min(1.0, f * 1.05) for r in picked])
    batches, cur, acc = [], [], 0.0
    for k in range(len(kept)):
        if cur and acc + kd_nnz[k] > batch_nnz:
            batches.append(cur); cur, acc = [], 0.0
        cur.append(k); acc += kd_nnz[k]
    batches.append(cur)
    log(f'  {len(batches)} pass(es) over the file (batch budget {batch_nnz:.2e} nnz)')

    G = len(var); n_src = len(src_ids); nk = len(kept)
    lfc = np.full((nk, G), np.nan, np.float32); se = np.full_like(lfc, np.nan)
    bulk_full = np.full_like(lfc, np.nan); bulk_matched = np.full_like(lfc, np.nan)
    sig_rows, sig_cols = [], []
    eff = np.full(nk, np.nan); n_real = np.zeros(nk, int); n_used = np.array([len(r) for r in picked])
    # pass-1 accumulators: raw-count sums per kept knockdown (all kept cells) and controls; full-depth control CPM;
    # per-cell raw row sums and the count dtype checks (the count check)
    raw_sum = np.zeros((nk + 1, n_src)); cpm0_full = np.zeros(n_src)
    rowsum_raw = np.zeros(len(obs)); count_flags = dict(dtype=None, non_integer=0, negative=0)
    grp_all = np.where(ctrl, nk, kd)
    panel_cols = np.flatnonzero(on)

    def accumulate(s, e, D, I, P, n_vars):
        count_flags['dtype'] = str(D.dtype)
        if D.dtype.kind == 'f': count_flags['non_integer'] += int((D != np.floor(D)).sum())
        count_flags['negative'] += int((D < 0).sum())
        rowsum_raw[s:e] = np.bincount(np.repeat(np.arange(e - s), np.diff(P)), weights=D, minlength=e - s)
        g = grp_all[s:e]; r = np.flatnonzero(g >= 0)
        if not len(r): return
        M = sp.csr_matrix((D.astype(np.float64), I, P), shape=(e - s, n_vars))[r]
        Gm = sp.csr_matrix((np.ones(len(r)), (g[r], np.arange(len(r)))), shape=(nk + 1, len(r)))
        raw_sum[:] += (Gm @ M).toarray()
        c = np.flatnonzero(g[r] == nk)
        if len(c):
            Mc = M[c]; lib = np.asarray(Mc[:, panel_cols].sum(1)).ravel()
            cpm0_full[:] += np.asarray((sp.diags(1e6 / lib) @ Mc).sum(0)).ravel()

    target_src = np.array([{g: i for i, g in enumerate(src_ids)}.get(t, -1) for t in kept.target_gene_id]) \
        if dataset not in ORION else np.array([{s: i for i, s in enumerate(var_src.index.astype(str))}.get(t, -1)
                                               for t in kept.target_gene])
    target_col = np.where(target_src >= 0, to_panel[np.maximum(target_src, 0)], -1)
    measurable = kept.measurable.astype(bool).to_numpy()
    ctrl_stats = None
    for bi, batch in enumerate(batches):
        want = np.zeros(len(obs), bool)
        for k in batch: want[picked[k]] = True
        if bi == 0: want[ctrl_rows] = True
        rows, data, indices, indptr, n_vars = stream(path, f, want, accumulate if bi == 0 else None)
        assert n_vars < 65536, f"gene index overflow: {n_vars}"
        pos = {r: i for i, r in enumerate(rows)}
        if bi == 0:
            # count check, right after the first pass (every row's raw total is known by now), so a file that is not
            # per-cell integer counts -- e.g. pre-aggregated pseudobulk sums -- fails in minutes, not after all batches
            count_check = check_counts(path, obs, rowsum_raw, n_umi, count_flags)
            ci = np.array([pos[r] for r in ctrl_rows])
            d0, i0, p0 = sub_csr(data, indices, indptr, ci)
            C, lib0 = cpm_csr(d0, i0, p0, n_vars, on)
            lib0_matrix = np.bincount(np.repeat(np.arange(len(ci)), np.diff(p0)), weights=d0, minlength=len(ci))
            raw0 = np.bincount(i0, weights=d0, minlength=n_vars)               # thinned counts, sampled controls
            n0 = C.shape[0]
            s1 = np.asarray(C.sum(0)).ravel(); s2 = np.asarray(C.multiply(C).sum(0)).ravel()
            m0 = s1 / n0; v0 = np.maximum(s2 - s1 ** 2 / n0, 0) / (n0 - 1)
            det0 = np.diff(C.tocsc().indptr) / n0
            tested = np.flatnonzero(on & (m0 > MIN_CPM))                   # source columns the scorer would test
            ptr, cvals, cnz, ctie = control_index(C, tested)
            q_path, q_unique = save_quantiles(source, tested, to_panel, ptr, cvals, cnz, n0)
            ctrl_stats = dict(n0=n0, m0=m0, v0=v0, det0=det0, tested=tested, lib_median=float(np.median(lib0)),
                              lib_matrix_median=float(np.median(lib0_matrix)), raw0=raw0)
            log(f'  controls: {n0:,} cells, median panel library {np.median(lib0):.0f} (whole matrix '
                f'{np.median(lib0_matrix):.0f}), {len(tested):,} genes pass the {MIN_CPM:g}-CPM gate; '
                f'{q_unique:,} unique control values -> {q_path.name}')
            del C, d0, i0
        t0 = time.time()
        for j, k in enumerate(batch):
            ri = np.array([pos[r] for r in picked[k]])
            d1, i1, p1 = sub_csr(data, indices, indptr, ri)
            raw1 = np.bincount(i1, weights=d1, minlength=n_vars)               # thinned counts, sampled cells
            bulk_matched[k, to_panel[on]] = np.log1p(BULK_TS * raw1[on] / max(raw1[on].sum(), 1))
            M, _ = cpm_csr(d1, i1, p1, n_vars, on)
            n1 = M.shape[0]; Md = M.toarray()
            m1 = Md.mean(0); v1 = Md.var(0, ddof=1)
            l = np.log2((m1 + EPS) / (m0 + EPS))
            e_ = np.sqrt(v1 / n1 / (m1 + EPS) ** 2 + v0 / n0 / (m0 + EPS) ** 2) / np.log(2)
            e_[(m1 == 0) | (m0 == 0)] = np.inf
            lfc[k, to_panel[on]] = l[on]; se[k, to_panel[on]] = e_[on]
            if target_src[k] >= 0 and measurable[k]: eff[k] = l[target_src[k]]
            PT = np.ascontiguousarray(np.log1p(Md[:, tested]).T)
            p = mwu_vs_control(PT, ptr, cvals, cnz, ctie, n0)
            s = tested[bh_sig(p)]; s = s[s != target_src[k]]
            n_real[k] = len(s); sig_rows.append(np.full(len(s), k)); sig_cols.append(to_panel[s])
            if (j + 1) % 500 == 0 or j + 1 == len(batch):
                log(f'  batch {bi + 1}/{len(batches)}: {j + 1:,} / {len(batch):,} knockdowns ({time.time() - t0:.0f}s)')
        del data, indices
    # full-depth profiles: log1p(5e4 * sum / total) over the panel-mapped genes, all kept cells, raw counts
    tot = raw_sum[:, on].sum(1, keepdims=True)
    prof = np.log1p(BULK_TS * raw_sum[:, on] / np.where(tot > 0, tot, 1))
    bulk_full[:, to_panel[on]] = prof[:-1]
    raw0 = ctrl_stats['raw0']; prof0_matched = np.log1p(BULK_TS * raw0[on] / max(raw0[on].sum(), 1))
    sig = sp.csr_matrix((np.ones(sum(len(r) for r in sig_rows), bool),
                         (np.concatenate(sig_rows), np.concatenate(sig_cols))), shape=(nk, G))
    base = pd.DataFrame(dict(control_mean_cpm=np.nan, detection_rate=np.nan, control_mean_cpm_fulldepth=np.nan,
                             log1p_bulk_50k_matched=np.nan, log1p_bulk_50k_full=np.nan, de_tested=False), index=var.index)
    pc = to_panel[on]
    base.iloc[pc, 0] = ctrl_stats['m0'][on]; base.iloc[pc, 1] = ctrl_stats['det0'][on]
    base.iloc[pc, 2] = (cpm0_full / ctrl.sum())[on]; base.iloc[pc, 3] = prof0_matched; base.iloc[pc, 4] = prof[-1]
    base.iloc[to_panel[ctrl_stats['tested']], 5] = True
    obs_out = pd.DataFrame(dict(target_gene=kept.target_gene.astype(str), target_gene_id=kept.target_gene_id.astype(str),
                                n_cells_kept=n_kept.astype(int), n_cells_used=n_used, efficiency=eff,
                                efficiency_known=measurable, n_real=n_real.astype(float), de_labels_known=True,
                                modality='single_cell', significance_unknown=False, lfc_type=LFC_TYPE, se_type=SE_TYPE))
    # library check: panel-only vs whole-matrix library of the sampled, thinned controls
    lp, lm = ctrl_stats['lib_median'], ctrl_stats['lib_matrix_median']; rel = lm / lp - 1 if lp else float('nan')
    library_check = dict(median_panel_library_after_thinning=lp, median_matrix_library_after_thinning=lm,
                         median_n_umi_times_f=float(np.median(n_umi[ctrl_rows]) * f), relative_difference=rel,
                         source_genes=int(n_src), source_genes_on_panel=int(on.sum()),
                         within_1pct=bool(abs(rel) <= LIB_TOL),
                         statement=(f'panel-only library matches the whole matrix within {LIB_TOL:.0%} ({rel:+.2%}): it is the '
                                    'matrix the scorer would see' if abs(rel) <= LIB_TOL else
                                    f'panel-only library is {rel:+.2%} from the whole QC\'d matrix ({int(n_src - on.sum()):,} '
                                    'source genes are off the panel or unmeasured); kept as the intended definition (the scorer '
                                    'normalizes over the panel), documented here'))
    meta = dict(cell_file=str(path), **qc_files, thinning_factor=f, median_control_n_umi=float(np.median(n_umi[ctrl])),
                median_control_panel_library_after_thinning=lp, control_cells=int(ctrl.sum()),
                control_cells_used=int(ctrl_stats['n0']), source_genes=int(n_src), source_genes_on_panel=int(on.sum()),
                de_tested_genes=int(len(ctrl_stats['tested'])), de_labels='wilcoxon_bh', modality='single_cell',
                significance_unknown=False, control_quantiles=q_path.name,
                wilcoxon_reconstruction=('possible: rank a predicted knockdown (cell values, or a predicted first moment turned '
                                         'into cell values) against ' + q_path.name + ' (sorted unique nonzero log1p CPM of '
                                         'the sampled controls, multiplicities, zeros, per gated gene) with the scorer\'s '
                                         'two-sided tie + continuity corrected MWU, then BH over the gated genes'),
                note=('H1: training table only; the scoring harness keeps the unfiltered real H1 cells' if dataset == 'H1'
                      else 'Orion: the A-C targets subset only' if dataset in ORION else ''))
    a = write(source, background, 'none', var, obs_out, lfc, se, bulk_matched, bulk_full, sig, base, meta)
    return a, dict(picked=picked, ctrl_rows=ctrl_rows, f=f, to_panel=to_panel, src_ids=src_ids, target_src=target_src,
                   kept=kept, path=path, quantiles=q_path, library_check=library_check, count_check=count_check)


MAX_ROW_TOTAL, ROW_MEDIAN_RANGE = 1e6, (500, 100_000)      # a single cell; pseudobulk sums run far above both


def check_counts(path, obs, rowsum_raw, n_umi, count_flags):
    """Count check: per-cell integer counts, one row per cell (not pre-aggregated sums). Exits on failure."""
    ratio = rowsum_raw / np.maximum(n_umi, 1)
    med, mx = float(np.median(rowsum_raw)), float(rowsum_raw.max())
    c = dict(dtype=count_flags['dtype'], integer_dtype=bool(np.dtype(count_flags['dtype']).kind in 'iu'),
             non_integer_values=count_flags['non_integer'], negative_values=count_flags['negative'],
             rows=int(len(obs)), unique_cell_ids=int(obs.index.nunique()),
             median_row_total=med, max_row_total=mx,
             max_row_total_below_1e6=bool(mx < MAX_ROW_TOTAL),
             median_row_total_single_cell_range=bool(ROW_MEDIAN_RANGE[0] <= med <= ROW_MEDIAN_RANGE[1]),
             single_cell_bounds=f'max row total < {MAX_ROW_TOTAL:.0e}, median row total in {ROW_MEDIAN_RANGE}',
             share_rows_total_equals_n_umi=float(np.mean(rowsum_raw == n_umi)),
             median_row_total_over_n_umi=float(np.median(ratio)))
    c['per_cell_integer_counts'] = bool(c['integer_dtype'] and count_flags['non_integer'] == 0 and count_flags['negative'] == 0
                                        and c['unique_cell_ids'] == len(obs) and c['max_row_total_below_1e6']
                                        and c['median_row_total_single_cell_range'])
    if not c['per_cell_integer_counts']:
        sys.exit(f'error: {path} does not hold per-cell integer counts: {c}')
    log(f'  count check passed: {c["dtype"]}, median row total {med:.0f}, max {mx:.0f}')
    return c


def sub_csr(data, indices, indptr, rows):
    """Rows of a (data, indices, indptr) CSR as a new (data, indices, indptr)."""
    starts, ends = indptr[rows], indptr[rows + 1]
    sel = np.concatenate([np.arange(a, b) for a, b in zip(starts, ends)]) if len(rows) else np.zeros(0, int)
    return data[sel], indices[sel].astype(np.int32), np.concatenate([[0], np.cumsum(ends - starts)])


# ---------------------------------------------------------------- writing
def write(source, background, condition, var, obs, lfc, se, bulk_matched, bulk_full, sig, base, meta):
    obs = obs.copy(); obs.insert(0, 'condition', condition); obs.insert(0, 'background', background); obs.insert(0, 'source', source)
    obs.index = pd.Index([f'{source}|{background}|{condition}|{t}' for t in obs.target_gene_id], name='row_id')
    for c in ('source', 'background', 'condition', 'modality', 'lfc_type', 'se_type'): obs[c] = obs[c].astype('category')
    layers = dict(lfc=lfc.astype(np.float32, copy=False), se=se.astype(np.float32, copy=False),
                  log1p_bulk_50k_matched=bulk_matched.astype(np.float32, copy=False),
                  log1p_bulk_50k_full=bulk_full.astype(np.float32, copy=False), de_sig=sig)
    a = ad.AnnData(obs=obs, var=var.copy(), layers=layers)
    a.uns['effects_v5'] = dict(meta, built=time.strftime('%Y-%m-%d %H:%M'), epsilon=EPS, fdr=FDR, p_adj_rule='p_adj < 0.05',
                               filter_gene_min_cpm_cell=MIN_CPM, gate='control mean CPM > 5 (rule_version 3)',
                               bulk_target_sum=BULK_TS, target_depth=DEPTH, max_cells=N_CELLS, max_control_cells=N_CTRL,
                               seed=SEED, scorer='cell_eval2 0.18 vcc2026 preset (pdex / gpudge MWU)',
                               test='two-sided Mann-Whitney U, normal approx, tie + continuity corrected, on log1p(CPM)',
                               target_excluded_from_de_sig=True, library='sum over the source genes that map onto the panel',
                               primary_profile='log1p_bulk_50k_matched',
                               profiles=('log1p_bulk_50k_matched: the sampled, thinned cells used for lfc and the test, so it '
                                         'carries the sampling variance of the single ~400-cell / ~20k-UMI draw the scorer '
                                         'sees (PRIMARY); log1p_bulk_50k_full: all kept cells, raw counts (the expected '
                                         'profile)'),
                               se_role='modelling aid only (delta method); significance comes from de_sig / the control quantiles')
    OUT.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(OUT / f'{source}.h5ad', compression='gzip')
    base.to_csv(OUT / f'{source}_baseline.csv')
    log(f'wrote {OUT / f"{source}.h5ad"}: {a.n_obs:,} knockdowns x {a.n_vars:,} panel genes, '
        f'{int(np.isfinite(lfc).any(0).sum()):,} measured')
    return a


# ---------------------------------------------------------------- validation (cell sources)
def validate(dataset, a, ctx, n_fast=5, n_genes=20, n_full=0):
    """Independent recomputation: re-read and re-thin the sampled cells of random knockdowns and the sampled controls
    straight from the file (one read for both modes), then (a) cell_eval2's compute_de (pdex), (b) scipy MWU + BH and
    (c) a numpy rank sum against the stored control quantiles. fast: the first n_fast knockdowns x n_genes random
    measured genes; full (n_full > 0): n_full knockdowns (all when fewer) x every measured gene."""
    from cell_eval2.de_compute import compute_de
    from scipy.stats import mannwhitneyu
    rng = np.random.default_rng(SEED + 1)
    kept, picked, f = ctx['kept'], ctx['picked'], ctx['f']
    to_panel, src_ids, target_src = ctx['to_panel'], ctx['src_ids'], ctx['target_src']
    on = np.flatnonzero(to_panel >= 0)
    draw = rng.permutation(len(kept))[:min(max(n_fast, n_full), len(kept))]
    fast = set(draw[:n_fast].tolist()); ks = np.sort(draw)
    q = dict(np.load(ctx['quantiles']))                                # eager: every array read once, not per access
    pos_on = {c: j for j, c in enumerate(to_panel[on])}                  # panel col -> position in `on`
    q_on = np.array([pos_on[c] for c in q['panel_col']])
    rows = np.concatenate([picked[k] for k in ks] + [ctx['ctrl_rows']])
    labels = np.concatenate([[kept.target_gene_id[k]] * len(picked[k]) for k in ks] + [['non-targeting'] * len(ctx['ctrl_rows'])])
    order = np.argsort(rows); rows, labels = rows[order], labels[order]
    t0 = time.time()
    with h5py.File(ctx['path'], 'r', rdcc_nbytes=1 << 30) as fh:
        X = fh['X']; n_vars = int(X.attrs['shape'][1])
        d, i, p = read_rows(X, rows, f)
    M = sp.csr_matrix((d.astype(np.float64), i, p), shape=(len(rows), n_vars))[:, on]
    M.sort_indices()
    log(f'  validation: re-read {len(rows):,} cells ({time.time() - t0:.0f}s)')
    cells = ad.AnnData(X=M, obs=pd.DataFrame(dict(target=labels), index=[str(r) for r in rows]),
                       var=pd.DataFrame(index=[f'g{j}' for j in range(len(on))]))
    df = compute_de(cells, backend='pdex', groupby='target', reference='non-targeting', mean_calc='arithmetic',
                    epsilon=EPS, input_type='counts', target_sum=1e6, clip_value=None, filter_gene_min_cpm_cell=MIN_CPM,
                    fdr_scope='per_pert', threads=-1).to_pandas()
    col_of = {s: j for j, s in enumerate(cells.var_names)}
    # scipy on the same cells
    lib = np.asarray(M.sum(1)).ravel()
    cpm = sp.csr_matrix((M.data / np.repeat(lib / 1e6, np.diff(M.indptr)), M.indices, M.indptr), shape=M.shape)
    ctrl = labels == 'non-targeting'; C = cpm[np.flatnonzero(ctrl)]
    m0 = np.asarray(C.mean(0)).ravel(); gated = np.flatnonzero(m0 > MIN_CPM)
    Cd = C[:, gated].toarray()
    res = []
    lfc_tab = a.layers['lfc']; sig_tab = a.layers['de_sig'].tocsr()
    for k in ks:
        t = kept.target_gene_id[k]; rk = np.flatnonzero(labels == t); Pk = cpm[rk]
        m1 = np.asarray(Pk.mean(0)).ravel()
        direct = np.log2((m1 + EPS) / (m0 + EPS))
        stored = lfc_tab[k, to_panel[on]].astype(float)                      # every measured gene, in `on` order
        sub = df[df.target == t]
        ce_lfc = sub.set_index('feature').log2_fold_change.reindex(cells.var_names).to_numpy(float)
        genes = rng.choice(len(on), n_genes, replace=False)                 # the fast check's 20 genes
        d_all, d_ce = np.abs(direct - stored), np.abs(ce_lfc - stored)
        tgt = target_src[k]; drop = {to_panel[tgt]} if tgt >= 0 else set()
        ce_sig = {to_panel[on[col_of[g]]] for g in sub.loc[sub.p_adj < FDR, 'feature']} - drop
        Pd = Pk[:, gated].toarray()
        pv = np.ones(len(gated))
        for s in range(0, len(gated), 2000):
            # on log1p(CPM), what the build and pdex rank (raw CPM is not rank-identical in float64)
            r_ = mannwhitneyu(np.log1p(Pd[:, s:s + 2000]), np.log1p(Cd[:, s:s + 2000]),
                              use_continuity=True, alternative='two-sided', method='asymptotic', axis=0).pvalue
            pv[s:s + 2000] = np.where(np.isnan(r_), 1.0, r_)
        sp_sig = {to_panel[on[gated[j]]] for j in np.flatnonzero(false_discovery_control(pv, method='bh') < FDR)} - drop
        # (c) the stored control quantiles alone, against this knockdown's re-read cells
        Lk = Pk[:, q_on].tocsc(); qp = np.ones(len(q_on))
        for g in range(len(q_on)):
            x = np.zeros(len(rk)); a_, b_ = Lk.indptr[g], Lk.indptr[g + 1]; x[Lk.indices[a_:b_]] = np.log1p(Lk.data[a_:b_])
            v, c = q['values'][q['ptr'][g]:q['ptr'][g + 1]], q['counts'][q['ptr'][g]:q['ptr'][g + 1]]
            qp[g] = quantile_p(x, v, c, int(q['zeros'][g]), int(q['n_controls']))
        q_sig = {int(q['panel_col'][j]) for j in np.flatnonzero(false_discovery_control(qp, method='bh') < FDR)} - drop
        stored_sig = set(sig_tab[k].indices.tolist())
        # cell_eval2's table holds only the gated genes, so its LFC is compared on those (direct: every measured gene)
        ce20 = d_ce[genes][np.isfinite(d_ce[genes])]
        res.append(dict(knockdown=str(kept.target_gene[k]), fast=k in fast, n_cells=int(len(picked[k])), n_stored=len(stored_sig),
                        lfc_max_abs_diff_20genes=float(max(d_all[genes].max(), ce20.max() if len(ce20) else 0.0)),
                        lfc_max_abs_diff_direct=float(d_all.max()), lfc_max_abs_diff_cell_eval2=float(np.nanmax(d_ce)),
                        symdiff_cell_eval2=len(stored_sig ^ ce_sig), symdiff_scipy=len(stored_sig ^ sp_sig),
                        symdiff_quantiles=len(stored_sig ^ q_sig)))
    fr = [r for r in res if r['fast']]
    out = dict(fast=dict(knockdowns=fr, lfc_max_abs_diff=max(r['lfc_max_abs_diff_20genes'] for r in fr),
                         lfc_pass=bool(max(r['lfc_max_abs_diff_20genes'] for r in fr) < 1e-5),
                         sig_pass=bool(all(r['symdiff_cell_eval2'] == 0 and r['symdiff_scipy'] == 0 for r in fr)),
                         quantiles_pass=bool(all(r['symdiff_quantiles'] == 0 for r in fr))),
               gated_genes_validation=int(len(gated)), gated_genes_build=int(a.uns['effects_v5']['de_tested_genes']),
               gated_genes_quantile_file=int(len(q_on)), gated_sets_identical=bool(set(on[gated]) == set(on[q_on])))
    if n_full:
        out['full'] = dict(
            knockdowns_checked=len(res),
            genes_per_knockdown='every measured gene (direct recomputation); the gated genes (cell_eval2, scipy, quantiles)',
            lfc_max_abs_diff_direct=max(r['lfc_max_abs_diff_direct'] for r in res),
            lfc_max_abs_diff_cell_eval2=max(r['lfc_max_abs_diff_cell_eval2'] for r in res),
            **{f'{w}_knockdowns_with_symdiff': sum(r[f'symdiff_{w}'] > 0 for r in res) for w in ('cell_eval2', 'scipy', 'quantiles')},
            **{f'{w}_total_symdiff_genes': sum(r[f'symdiff_{w}'] for r in res) for w in ('cell_eval2', 'scipy', 'quantiles')},
            per_knockdown=res)
        out['full']['pass'] = bool(max(out['full']['lfc_max_abs_diff_direct'], out['full']['lfc_max_abs_diff_cell_eval2']) < 1e-5
                                   and all(out['full'][f'{w}_total_symdiff_genes'] == 0 for w in ('cell_eval2', 'scipy', 'quantiles')))
    return out


# ---------------------------------------------------------------- Marson CD4 T tables
def build_table(name):
    source, background, condition = TABLES[name]
    var, _ = panel(); sym_col = {s: i for i, s in enumerate(var.index)}
    kq = V / 'results/qc' / source / 'knockdown_qc.csv'
    qc = pd.read_csv(kq); qc['target_gene'] = qc['target_gene'].astype(str)
    kept = qc[qc.knockdown_keep.astype(bool)].reset_index(drop=True)
    path = V / f'lib/tcell/TCELL_{condition}_effects.h5ad'
    G = len(var); lfc = np.full((len(kept), G), np.nan, np.float32); se = np.full_like(lfc, np.nan)
    with h5py.File(path, 'r') as fh:
        obs = read_elem(fh['obs']); vv = read_elem(fh['var']); pos = {t: i for i, t in enumerate(obs.index.astype(str))}
        cols = np.array([sym_col.get(s, -1) for s in vv.index.astype(str)]); meas = vv['measured'].to_numpy(bool) & (cols >= 0)
        rows = np.array([pos[t] for t in kept.target_gene]); o = np.argsort(rows); back = np.argsort(o)
        Xl = fh['X'][rows[o]][back]; Xs = fh['layers/se'][rows[o]][back]
    lfc[:, cols[meas]] = Xl[:, meas]; se[:, cols[meas]] = Xs[:, meas]
    target_col = np.array([sym_col.get(t, -1) for t in kept.target_gene])
    measurable = kept.measurable.astype(bool).to_numpy()
    eff = np.array([lfc[i, j] if (j >= 0 and measurable[i]) else np.nan for i, j in enumerate(target_col)])
    tid = [var.gene_id.iloc[sym_col[t]] if t in sym_col and var.gene_id.iloc[sym_col[t]] else f'symbol:{t}' for t in kept.target_gene]
    obs_out = pd.DataFrame(dict(target_gene=kept.target_gene, target_gene_id=tid, n_cells_kept=kept.n_cells.astype(float),
                                n_cells_used=np.nan, efficiency=eff, efficiency_known=measurable, n_real=np.nan,
                                de_labels_known=False, modality='replicate_table', significance_unknown=True,
                                lfc_type='marson_replicate_paired_pc0.5_NOT_OFFICIAL',
                                se_type='marson_donor_moderated_table_NOT_DOWNSAMPLED'))
    base = pd.DataFrame(dict(control_mean_cpm=np.nan, detection_rate=np.nan, control_mean_cpm_fulldepth=np.nan,
                             log1p_bulk_50k_matched=np.nan, log1p_bulk_50k_full=np.nan, de_tested=False), index=var.index)
    ctrl_cpm = vv['ctrl_cpm'].to_numpy(float)
    base.iloc[cols[meas], 2] = ctrl_cpm[meas]
    base.iloc[cols[meas], 4] = np.log1p(BULK_TS * ctrl_cpm[meas] / ctrl_cpm[meas].sum())
    meta = dict(table=str(path), knockdown_qc=str(kq), de_labels='absent', modality='replicate_table',
                significance_unknown=True, control_quantiles='none (no single cells)',
                wilcoxon_reconstruction='not possible: no single cells, no per-cell ranks (none are invented)',
                method='replicate-level table only: no single cells, so no Wilcoxon / BH labels (de_sig empty, n_real NaN, '
                       'de_labels_known False, modality replicate_table, significance_unknown True; treat significance '
                       'as unknown, like the Tahoe modality flag). lfc = the '
                       "table's replicate-paired log2 ratio of pseudobulk CPM with pseudocount 0.5 (NOT the official "
                       'epsilon 1e-9 LFC: the per-knockdown means are not in the tables we hold); se = the table SE '
                       '(between-replicate, limma-style moderated), not floored, not in the downsampled regime. '
                       'Profiles: per knockdown absent (NaN, both layers); the baseline carries a full-depth control '
                       "profile from the table's pooled control CPM (approximate); no matched profile.")
    sig = sp.csr_matrix((len(kept), G), dtype=bool)
    nan = np.full_like(lfc, np.nan)
    a = write(source, background, condition, var, obs_out, lfc, se, nan, nan, sig, base, meta)
    return a


# ---------------------------------------------------------------- main
if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--source', required=True, choices=list(CELL_SOURCES) + list(TABLES))
    p.add_argument('--batch-nnz', type=float, default=1.2e9, help='sampled nonzeros held per pass (4-6 bytes each)')
    p.add_argument('--full-validate', type=int, default=0, metavar='N',
                   help='also validate N random knockdowns (all when fewer) on every measured gene')
    args = p.parse_args(); t0 = time.time()
    if args.source in TABLES:
        a = build_table(args.source); name = TABLES[args.source][0]
        checks = dict(epsilon=None, de_labels='absent', de_labels_known_flag_set=bool((~a.obs.de_labels_known).all()),
                      modality_flag_set=bool((a.obs.modality == 'replicate_table').all() and a.obs.significance_unknown.all()),
                      statement='Marson CD4 T: no Wilcoxon labels (replicate-level tables only); lfc is the pseudocount-0.5 '
                                'replicate-paired ratio, flagged NOT_OFFICIAL')
        notes = dict(library_check='not applicable (no cells; the table CPMs are the source\'s own)',
                     control_quantiles='not stored (no single cells)',
                     modality="obs and uns: modality = 'replicate_table', significance_unknown = True, de_labels_known = False",
                     profiles='both layers NaN per knockdown (no per-knockdown means); baseline log1p_bulk_50k_full from the '
                              'pooled control CPM (approximate)',
                     validation='not applicable (no cells)', orion_path='not applicable')
    else:
        a, ctx = build_cells(args.source, args.batch_nnz); name = CELL_SOURCES[args.source]
        log('validating' + (f' (fast + full on {args.full_validate} knockdowns)' if args.full_validate else ' (fast)'))
        checks = validate(args.source, a, ctx, n_full=args.full_validate)
        checks.update(epsilon=EPS, epsilon_is_1e9=a.uns['effects_v5']['epsilon'] == 1e-9,
                      max_n_cells_used=int(a.obs.n_cells_used.max()), n_cells_used_le_400=bool(a.obs.n_cells_used.max() <= N_CELLS),
                      target_depth=DEPTH, thinning_factor=ctx['f'], library_check=ctx['library_check'],
                      count_check=ctx['count_check'])
        fv = checks.get('full')
        notes = dict(
            library_check=ctx['library_check']['statement'],
            control_quantiles=f"{ctx['quantiles'].name}: {checks['gated_genes_quantile_file']:,} gated genes; reproduces the "
                              f"stored de_sig on the validated knockdowns: {checks['fast']['quantiles_pass'] and (fv is None or fv['quantiles_total_symdiff_genes'] == 0)}",
            modality="obs and uns: modality = 'single_cell', significance_unknown = False",
            profiles='log1p_bulk_50k_matched (primary: the sampled, thinned cells) and log1p_bulk_50k_full (all kept cells)',
            validation=(f"fast 5 x 20: lfc {checks['fast']['lfc_pass']}, sig {checks['fast']['sig_pass']}; " +
                        (f"full {fv['knockdowns_checked']} knockdowns: max |dLFC| {max(fv['lfc_max_abs_diff_direct'], fv['lfc_max_abs_diff_cell_eval2']):.2e}, "
                         f"knockdowns with symdiff cell_eval2 / scipy / quantiles {fv['cell_eval2_knockdowns_with_symdiff']} / "
                         f"{fv['scipy_knockdowns_with_symdiff']} / {fv['quantiles_knockdowns_with_symdiff']}, genes "
                         f"{fv['cell_eval2_total_symdiff_genes']} / {fv['scipy_total_symdiff_genes']} / {fv['quantiles_total_symdiff_genes']}"
                         if fv else 'full: not run for this source')),
            orion_path=(f"QC'd file holds per-cell integer counts ({ctx['count_check']['dtype']}, {ctx['count_check']['rows']:,} rows "
                        f"= cells, median row total {ctx['count_check']['median_row_total']:.0f}); the same thinning + sampling + "
                        'Wilcoxon path as every cell source' if args.source in ORION else 'not an Orion source'))
    if args.source in TABLES:
        required = dict(de_labels_known_flag_set=checks['de_labels_known_flag_set'], modality_flag_set=checks['modality_flag_set'])
    else:
        required = dict(fast_lfc_pass=checks['fast']['lfc_pass'], fast_sig_pass=checks['fast']['sig_pass'],
                        fast_quantiles_pass=checks['fast']['quantiles_pass'], gated_sets_identical=checks['gated_sets_identical'],
                        epsilon_is_1e9=checks['epsilon_is_1e9'], n_cells_used_le_400=checks['n_cells_used_le_400'])
        if args.full_validate: required['full_pass'] = checks['full']['pass']
    failing = [k for k, v in required.items() if not v]
    all_pass = not failing
    lfc = a.layers['lfc']
    summary = dict(source=name, rows=int(a.n_obs), measured_genes=int(np.isfinite(lfc).any(0).sum()),
                   de_tested_genes=a.uns['effects_v5'].get('de_tested_genes'),
                   median_n_real=float(a.obs.n_real.median()) if a.obs.n_real.notna().any() else None,
                   p90_n_real=float(a.obs.n_real.quantile(0.9)) if a.obs.n_real.notna().any() else None,
                   median_efficiency=float(np.nanmedian(a.obs.efficiency)) if a.obs.efficiency.notna().any() else None,
                   lfc_on_epsilon_floor_share=float(np.mean(np.abs(lfc[np.isfinite(lfc)]) > 25)) if np.isfinite(lfc).any() else None,
                   all_pass=all_pass, required_checks=required, failing_checks=failing,
                   notes=notes, checks=checks, minutes=round((time.time() - t0) / 60, 1))
    (OUT / f'{name}_build.json').write_text(json.dumps(summary, indent=1, default=str) + '\n')
    print(json.dumps(summary, indent=1, default=str))
    if not all_pass:
        print(f'VALIDATION FAILED for {name}: {", ".join(failing)}', file=sys.stderr, flush=True)
        sys.exit(1)

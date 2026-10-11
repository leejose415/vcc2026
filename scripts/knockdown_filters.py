"""Per-knockdown, per-dataset QC filters: enough cells, target expressed in controls, knockdown worked.

For every knockdown in every training dataset:

  1. cells          n_cells >= 100
  2. expressed      the knocked-down gene is visible in that dataset's own control cells:
                      umi: >= 0.1 UMI per control cell    (sum of target counts / control cells)
                      cpm: >= 5 CPM                       (mean per-cell CPM, the evaluator's convention)
                    --expr-rule picks which one decides keep/drop (umi / cpm / both / either).
  3. knockdown      on-target log2 fold change <= -1 (at least a 50% drop) AND adjusted p < 0.05 within the
                    dataset (Benjamini-Hochberg over all of its knockdowns with a measured target).

On-target LFC and test, by how the dataset is stored:
  sums  (Orion HCT116/HEK293T, H1)  log2((mean CPM + 0.5) / (control mean CPM + 0.5)), the pipeline's label
                                    convention; Welch t-test on per-cell CPM from the stored moments.
  cells (Replogle K562 GWPS)        same LFC and Welch test, streamed from raw counts; CPM uses obs UMI_count
                                    as the library (the matrix holds only 8,248 of the measured genes).
  tcell (Marson CD4+ T cells)       the table's own replicate-paired LFC / SE; t-test with (n_reps - 1) + d0
                                    degrees of freedom (limma-style moderated SE). No UMI counts were
                                    released, so the umi rule falls back to cpm here.

A target the dataset did not measure fails rules 2 and 3: its knockdown cannot be verified.
Replogle knockdowns are grouped by target gene, pooling the gene's guides/transcripts.

Usage:  ~/venv/bin/python scripts/knockdown_filters.py [--datasets HCT116 K562_gwps ...] [--expr-rule umi]
Writes (lib/filters/):
  <dataset>.h5ad    kept knockdowns (rows) x genes (columns), X = log2 fold change vs controls;
                    obs = per-knockdown QC numbers, var = control expression of each gene
  _cache/           raw per-knockdown numbers; changing a threshold re-uses them (seconds).
                             The Replogle pass streams the 61 GiB matrix once (~7 min); --recompute redoes it.
"""
import argparse
import ast
import json
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import stats

try:
    from anndata.io import read_elem, sparse_dataset
except ImportError:
    from anndata.experimental import read_elem, sparse_dataset

V = Path.home() / 'vcc'
LIB = V / 'lib'
CONTROL = 'non-targeting'
PC = 0.5                   # CPM pseudocount for LFCs, same as the pipeline's labels
RULES = ('umi', 'cpm', 'both', 'either')
GENES = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy()   # challenge gene list

# name -> (storage kind, file, extra)
DATASETS = {
    'HCT116':       ('sums',  LIB / 'orion/HCT116_pseudobulk.h5ad', {}),
    'HEK293T':      ('sums',  LIB / 'orion/HEK293T_pseudobulk.h5ad', {}),
    'H1':           ('sums',  LIB / 'h1/h1_pseudobulk.h5ad', {}),
    'CD4_Stim48hr': ('tcell', LIB / 'tcell/TCELL_Stim48hr_effects.h5ad', {}),
    'CD4_Stim8hr':  ('tcell', LIB / 'tcell/TCELL_Stim8hr_effects.h5ad', {}),
    'CD4_Rest':     ('tcell', LIB / 'tcell/TCELL_Rest_effects.h5ad', {}),
    'K562_gwps':    ('cells', LIB / 'replogle/K562_gwps_raw_singlecell_01.h5ad',
                     dict(target_col='gene', target_id_col='gene_id', library_col='UMI_count')),
}
TESTS = dict(sums='Welch t on per-cell CPM', cells='Welch t on per-cell CPM',
             tcell='moderated t on replicate SE')


# ---------------------------------------------------------------- tests
def welch(m1, v1, n1, m0, v0, n0):
    """Two-sided Welch t-test from group means and sample variances (vectorised). No spread at all -> NaN."""
    se2 = v1 / n1 + v0 / n0
    with np.errstate(divide='ignore', invalid='ignore'):
        t = (m1 - m0) / np.sqrt(se2)
        df = se2 ** 2 / ((v1 / n1) ** 2 / np.maximum(n1 - 1, 1) + (v0 / n0) ** 2 / np.maximum(n0 - 1, 1))
        p = 2 * stats.t.sf(np.abs(t), df)
    return np.where(se2 > 0, p, np.nan)


def bh(p):
    """Benjamini-Hochberg adjusted p-values; NaN stays NaN and is not counted as a test."""
    p = np.asarray(p, float); out = np.full_like(p, np.nan); ok = np.isfinite(p)
    if ok.any():
        q = p[ok]; order = np.argsort(q); m = len(q)
        adj = np.minimum.accumulate((q[order] * m / np.arange(1, m + 1))[::-1])[::-1]
        res = np.empty(m); res[order] = np.minimum(adj, 1); out[ok] = res
    return out


def moments_lfc_p(n, s_cpm, s_cpm2, n0, s0_cpm, s0_cpm2):
    """On-target LFC and Welch p from per-group sums of per-cell CPM and CPM^2."""
    m1, m0 = s_cpm / n, s0_cpm / n0
    v1 = np.maximum(s_cpm2 - n * m1 ** 2, 0) / np.maximum(n - 1, 1)
    v0 = np.maximum(s0_cpm2 - n0 * m0 ** 2, 0) / max(n0 - 1, 1)
    return np.log2((m1 + PC) / (m0 + PC)), welch(m1, v1, n, m0, v0, n0)


def lfc(mean_cpm, ctrl_cpm):
    return np.log2((mean_cpm + PC) / (ctrl_cpm + PC))


def relvar(mean, var_of_mean, unit):
    """Relative variance of (mean + PC), with the variance of the mean floored at Poisson noise plus one count
    (`unit` = CPM one count adds to the mean), as in vcc_build_effects_v3.mean_variance: an all-zero group
    must not look precise."""
    return np.maximum(var_of_mean, (np.maximum(mean, 0) + unit) * unit) / (mean + PC) ** 2


def se_lfc(kd_relvar, ctrl_relvar):
    """Delta-method SE of the log2 fold change."""
    return np.sqrt(kd_relvar + ctrl_relvar) / np.log(2)


# ---------------------------------------------------------------- reading
def layer_reader(f, key):
    """(start, stop) or sorted row array -> dense block, for dense or sparse h5ad matrices."""
    node = f[key]
    if isinstance(node, h5py.Group):
        ds = sparse_dataset(node)
        return lambda rows: ds[rows].toarray()
    return lambda rows: node[rows]


def target_columns(names, var):
    """Each row's own-target column in var (-1 = not measured)."""
    measured = var['measured'].to_numpy(bool) if 'measured' in var else np.ones(len(var), bool)
    pos = {g: i for i, g in enumerate(var.index.astype(str))}
    col = np.array([pos.get(t, -1) for t in names])
    col[(col >= 0) & ~measured[np.maximum(col, 0)]] = -1
    return col, measured


def gene_columns(source_genes, measured):
    """Challenge genes this dataset measured, and their columns in the source."""
    pos = {}
    for i, g in enumerate(np.asarray(source_genes, str)):
        if measured[i]: pos.setdefault(g, i)
    keep = np.array([g in pos for g in GENES])
    return GENES[keep], np.array([pos[g] for g in GENES[keep]], dtype=np.int64)


def chunks(rows, size):
    for a in range(0, len(rows), size): yield rows[a:a + size]


# ---------------------------------------------------------------- per-knockdown numbers (cached)
def stats_sums(name, path, chunk=512, **_):
    """Orion / H1: one row per knockdown holding sums over its cells, plus one aggregated control row."""
    with h5py.File(path, 'r') as f:
        obs, var = read_elem(f['obs']), read_elem(f['var'])
        names = obs.index.astype(str).to_numpy(); n = obs['n_cells'].to_numpy(float)
        col, _ = target_columns(names, var)
        ci = int(np.flatnonzero(names == CONTROL)[0])
        L = {k: layer_reader(f, f'layers/{k}') for k in ('sum_counts', 'sum_cpm', 'sum_cpm2')}
        ctrl = {k: np.asarray(r(slice(ci, ci + 1)), float)[0] for k, r in L.items()}
        own = {k: np.full(len(names), np.nan) for k in L}           # each row's value at its own target gene
        for a in range(0, len(names), chunk):
            b = min(a + chunk, len(names)); c = col[a:b]; ok = c >= 0
            for k, r in L.items():
                block = np.asarray(r(slice(a, b)), float)
                own[k][a:b][ok] = block[np.flatnonzero(ok), c[ok]]
            print(f'  {name}: {b:,}/{len(names):,}', end='\r', flush=True)
    print()
    n0 = n[ci]; tc = np.maximum(col, 0)
    lfc_, p = moments_lfc_p(n, own['sum_cpm'], own['sum_cpm2'], n0, ctrl['sum_cpm'][tc], ctrl['sum_cpm2'][tc])
    out = pd.DataFrame(dict(
        target=names, n_cells=n, n_control_cells=n0, target_measured=col >= 0,
        ctrl_umi_per_cell=ctrl['sum_counts'][tc] / n0, ctrl_cpm=ctrl['sum_cpm'][tc] / n0,
        ontarget_lfc=lfc_, pvalue=p))
    out.loc[~out.target_measured, ['ctrl_umi_per_cell', 'ctrl_cpm', 'ontarget_lfc', 'pvalue']] = np.nan
    return out[names != CONTROL].reset_index(drop=True)


def stats_tcell(name, path, **_):
    """Marson T cells: an effect table (no cells). Uses its own on-target LFC/SE and control CPM."""
    with h5py.File(path, 'r') as f:
        obs, var = read_elem(f['obs']), read_elem(f['var'])
        settings = ast.literal_eval(str(read_elem(f['uns/settings'])))
    d0 = float(settings.get('d0', 0))
    names = obs.index.astype(str).to_numpy()
    col, _ = target_columns(names, var)
    lfc_ = obs['target_lfc'].to_numpy(float); se = obs['target_se'].to_numpy(float)
    has = (col >= 0) & np.isfinite(lfc_) & np.isfinite(se) & (se > 0)
    df = obs['n_reps'].to_numpy(float) - 1 + d0
    p = np.where(has, 2 * stats.t.sf(np.abs(lfc_ / np.where(has, se, 1)), np.maximum(df, 1)), np.nan)
    ctrl_cpm = np.where(col >= 0, var['ctrl_cpm'].to_numpy(float)[np.maximum(col, 0)], np.nan)
    out = pd.DataFrame(dict(
        target=names, n_cells=obs['n_cells'].to_numpy(float), n_control_cells=np.nan, target_measured=has,
        ctrl_umi_per_cell=np.nan, ctrl_cpm=np.where(has, ctrl_cpm, np.nan), ontarget_lfc=np.where(has, lfc_, np.nan),
        pvalue=p))
    return out[names != CONTROL].reset_index(drop=True)


def stats_cells(name, path, target_col, target_id_col, library_col, block=2048, **_):
    """Raw cells (Replogle): stream the matrix once. Keeps control moments, each cell's own-target count, and
    every group's per-cell CPM mean and variance over all genes (cached for the effects matrix and its SE)."""
    with h5py.File(path, 'r') as f:
        obs, var = read_elem(f['obs']), read_elem(f['var'])
        X = f['X']
        if isinstance(X, h5py.Group): raise ValueError('expected a dense X; add a sparse path for this dataset')
        labels = obs[target_col].astype(str).to_numpy()
        lib = obs[library_col].to_numpy(float)
        # match the target by Ensembl id (var index), falling back to the symbol
        by_id = {g: i for i, g in enumerate(var.index.astype(str))}
        by_name = {}
        for i, g in enumerate(var['gene_name'].astype(str)): by_name.setdefault(g, i)
        targets, code = np.unique(labels, return_inverse=True)
        first = pd.Series(np.arange(len(labels))).groupby(code).first().to_numpy()
        ids = obs[target_id_col].astype(str).to_numpy()[first]
        tcol = np.array([by_id.get(i, by_name.get(t, -1)) for t, i in zip(targets, ids)])
        cell_col = tcol[code]
        is_ctrl = labels == CONTROL
        G = X.shape[1]; c_cnt = np.zeros(G)
        g_cpm = np.zeros((len(targets), G)); g_cpm2 = np.zeros((len(targets), G))   # per group: sum CPM, CPM^2
        t0 = time.time()
        for a in range(0, X.shape[0], block):
            b = min(a + block, X.shape[0]); x = X[a:b].astype(np.float64)
            cc = x * (1e6 / np.maximum(lib[a:b], 1))[:, None]
            g = code[a:b]; order = np.argsort(g, kind='stable'); gs = g[order]
            starts = np.r_[0, np.flatnonzero(gs[1:] != gs[:-1]) + 1]
            co = cc[order]
            g_cpm[gs[starts]] += np.add.reduceat(co, starts, axis=0)
            g_cpm2[gs[starts]] += np.add.reduceat(co * co, starts, axis=0)
            ctl = is_ctrl[a:b]
            if ctl.any(): c_cnt += x[ctl].sum(0)
            if (a // block) % 50 == 0:
                rate = b / max(time.time() - t0, 1e-9)
                print(f'  {name}: {b:,}/{X.shape[0]:,} cells, ~{(X.shape[0] - b) / rate / 60:.0f} min left',
                      end='\r', flush=True)
        source_genes = var['gene_name'].astype(str).to_numpy()
    print()
    n = np.bincount(code, minlength=len(targets)).astype(float)
    ci = int(np.flatnonzero(targets == CONTROL)[0]); n0 = n[ci]; tc = np.maximum(tcol, 0)
    kd = np.arange(len(targets)) != ci
    own = np.arange(len(targets))
    lfc_, p = moments_lfc_p(n, g_cpm[own, tc], g_cpm2[own, tc], n0, g_cpm[ci][tc], g_cpm2[ci][tc])
    out = pd.DataFrame(dict(
        target=targets, n_cells=n, n_control_cells=n0, target_measured=tcol >= 0,
        ctrl_umi_per_cell=c_cnt[tc] / n0, ctrl_cpm=g_cpm[ci][tc] / n0,
        ontarget_lfc=lfc_, pvalue=p))
    out.loc[~out.target_measured, ['ctrl_umi_per_cell', 'ctrl_cpm', 'ontarget_lfc', 'pvalue']] = np.nan
    # per group and gene: mean and (population) variance of per-cell CPM, for the effects matrix and its SE
    mean = g_cpm / n[:, None]; g_cpm2 /= n[:, None]; g_cpm2 -= mean ** 2; var_ = np.maximum(g_cpm2, 0)
    unit = np.bincount(code, weights=1e6 / np.maximum(lib, 1), minlength=len(targets)) / n / n
    np.savez(matrix_cache(name), targets=targets[kd].astype(str), genes=source_genes.astype(str), n_cells=n[kd],
             mean_cpm=mean[kd].astype(np.float32), var_cpm=var_[kd].astype(np.float32), unit=unit[kd],
             ctrl_cpm=mean[ci], ctrl_var=var_[ci], ctrl_unit=unit[ci], ctrl_n=n0, ctrl_umi=c_cnt / n0)
    return out[kd].reset_index(drop=True)


BUILDERS = dict(sums=stats_sums, tcell=stats_tcell, cells=stats_cells)


def matrix_cache(name):
    return CACHE / f'{name}_mean_cpm.npz'


# ---------------------------------------------------------------- kept knockdowns x genes
def effects_sums(name, path, targets, chunk=256, **_):
    with h5py.File(path, 'r') as f:
        obs, var = read_elem(f['obs']), read_elem(f['var'])
        names = obs.index.astype(str).to_numpy(); n = obs['n_cells'].to_numpy(float)
        _, measured = target_columns(names, var)
        genes, src = gene_columns(var.index, measured)
        pos = {t: i for i, t in enumerate(names)}; ci = pos[CONTROL]
        cpm, cpm2, cnt = (layer_reader(f, f'layers/{k}') for k in ('sum_cpm', 'sum_cpm2', 'sum_counts'))

        def moments(r, nn):
            """rows -> mean CPM, variance of that mean, CPM of one count (from the group's total counts)"""
            counts = np.asarray(cnt(r), float)[:, src]
            m = np.asarray(cpm(r), float)[:, src] / nn
            v = np.maximum(np.asarray(cpm2(r), float)[:, src] / nn - m ** 2, 0) / nn
            return m, v, 1e6 / np.maximum(counts.sum(1, keepdims=True), 1), counts

        c_cpm, c_v, c_u, c_cnt = (x[0] for x in moments(slice(ci, ci + 1), n[ci]))
        c_rv = relvar(c_cpm, c_v, c_u)
        rows = np.array(sorted(pos[t] for t in targets), dtype=np.int64)
        X = np.empty((len(rows), len(genes)), np.float32); S = np.empty_like(X)
        for k, r in enumerate(chunks(rows, chunk)):
            m, v, u, _ = moments(r, n[r, None]); out = slice(k * chunk, k * chunk + len(r))
            X[out] = lfc(m, c_cpm); S[out] = se_lfc(relvar(m, v, u), c_rv)
    return names[rows], genes, X, S, dict(ctrl_cpm=c_cpm, ctrl_umi_per_cell=c_cnt / n[ci])


def effects_tcell(name, path, targets, chunk=512, **_):
    with h5py.File(path, 'r') as f:
        obs, var = read_elem(f['obs']), read_elem(f['var'])
        names = obs.index.astype(str).to_numpy()
        _, measured = target_columns(names, var)
        genes, src = gene_columns(var.index, measured)
        pos = {t: i for i, t in enumerate(names)}
        rows = np.array(sorted(pos[t] for t in targets), dtype=np.int64)
        rx, rs = layer_reader(f, 'X'), layer_reader(f, 'layers/se')
        X = np.empty((len(rows), len(genes)), np.float32); S = np.empty_like(X)
        for k, r in enumerate(chunks(rows, chunk)):
            out = slice(k * chunk, k * chunk + len(r))
            X[out] = np.asarray(rx(r), np.float32)[:, src]; S[out] = np.asarray(rs(r), np.float32)[:, src]
        c_cpm = var['ctrl_cpm'].to_numpy(float)[src]
    return names[rows], genes, X, S, dict(ctrl_cpm=c_cpm)


def effects_cells(name, path, targets, **_):
    z = np.load(matrix_cache(name), allow_pickle=False)
    genes, src = gene_columns(z['genes'], np.ones(len(z['genes']), bool))
    pos = {t: i for i, t in enumerate(z['targets'])}
    rows = np.array(sorted(pos[t] for t in targets), dtype=np.int64)
    c_cpm = z['ctrl_cpm'][src]
    c_rv = relvar(c_cpm, z['ctrl_var'][src] / z['ctrl_n'], z['ctrl_unit'])
    m = z['mean_cpm'][rows][:, src].astype(float)
    v = z['var_cpm'][rows][:, src].astype(float) / z['n_cells'][rows, None]
    X = lfc(m, c_cpm).astype(np.float32)
    S = se_lfc(relvar(m, v, z['unit'][rows, None]), c_rv).astype(np.float32)
    return z['targets'][rows], genes, X, S, dict(ctrl_cpm=c_cpm, ctrl_umi_per_cell=z['ctrl_umi'][src])


EFFECTS = dict(sums=effects_sums, tcell=effects_tcell, cells=effects_cells)


# ---------------------------------------------------------------- the three filters
def apply_filters(s, a):
    """Adds padj, status (keep/drop) and failed (which filters, comma-separated). Returns the table and
    keep masks under every expression rule (for the summary)."""
    s = s.copy()
    s['padj'] = bh(s['pvalue'])
    measured = s['target_measured'].to_numpy(bool)
    cells = (s['n_cells'] >= a.min_cells).to_numpy()
    cpm = (s['ctrl_cpm'] >= a.min_cpm).to_numpy()
    has_umi = s['ctrl_umi_per_cell'].notna().any()
    umi = (s['ctrl_umi_per_cell'] >= a.min_umi).to_numpy() if has_umi else cpm   # no UMI counts: cpm decides
    expr = dict(umi=umi, cpm=cpm, both=umi & cpm, either=umi | cpm)
    strong = (s['ontarget_lfc'] <= a.max_lfc).to_numpy()
    signif = (s['padj'] < a.alpha).to_numpy()
    keep = {r: cells & expr[r] & strong & signif for r in RULES}
    checks = [('cells', ~cells), ('not measured', ~measured), ('not expressed', measured & ~expr[a.expr_rule]),
              ('weak knockdown', measured & ~strong), ('not significant', measured & ~signif)]
    s['status'] = np.where(keep[a.expr_rule], 'keep', 'drop')
    s['failed'] = [', '.join(label for label, bad in checks if bad[i]) for i in range(len(s))]
    return s, keep, dict(checks), has_umi


COLUMNS = ['target', 'challenge_target', 'status', 'failed', 'n_cells', 'ctrl_umi_per_cell', 'ctrl_cpm',
           'ontarget_lfc', 'padj']


def write_dataset(name, kind, path, extra, s, a):
    table = s[COLUMNS].dropna(axis=1, how='all')                   # e.g. T cells have no UMI column
    kept = s.loc[s.status == 'keep', 'target'].astype(str).to_numpy()
    if not len(kept): return
    rows, genes, X, S, gene_info = EFFECTS[kind](name, path, kept, **extra)
    obs = table.set_index('target').loc[rows, ['challenge_target', 'n_cells', 'ontarget_lfc', 'padj']]
    var = pd.DataFrame(gene_info, index=pd.Index(genes, name='gene'))
    adata = ad.AnnData(X=X, layers=dict(se=S), obs=obs, var=var, uns=dict(
        dataset=name, source_file=str(path),
        X=f'log2((knockdown mean per-cell CPM + {PC}) / (control mean per-cell CPM + {PC}))'
          if kind != 'tcell' else 'replicate-paired log2 pseudobulk CPM ratio from the source table (pseudocount 0.5)',
        se='standard error of X: ' + ('between-donor, limma-style moderated (source table)' if kind == 'tcell' else
            'delta method from the cell-to-cell spread of per-cell CPM, floored at one Poisson count'),
        filters=json.dumps(settings(a))))
    adata.write_h5ad(a.out / f'{name}.h5ad', compression='gzip')
    print(f'{name}.h5ad: {X.shape[0]:,} knockdowns x {X.shape[1]:,} genes')


def settings(a):
    return dict(min_cells=a.min_cells, min_umi_per_control_cell=a.min_umi, min_ctrl_cpm=a.min_cpm,
                expr_rule=a.expr_rule, max_ontarget_lfc=a.max_lfc, max_padj=a.alpha)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--datasets', nargs='+', default=list(DATASETS), choices=list(DATASETS))
    ap.add_argument('--out', type=Path, default=LIB / 'filters')
    ap.add_argument('--min-cells', type=int, default=100)
    ap.add_argument('--min-umi', type=float, default=0.1, help='UMI per control cell')
    ap.add_argument('--min-cpm', type=float, default=5.0, help='mean per-cell CPM in controls')
    ap.add_argument('--expr-rule', choices=RULES, default='umi', help='which expression rule decides keep/drop')
    ap.add_argument('--max-lfc', type=float, default=-1.0, help='on-target log2FC must be <= this')
    ap.add_argument('--alpha', type=float, default=0.05)
    ap.add_argument('--recompute', action='store_true', help='ignore cached numbers')
    a = ap.parse_args()

    global CACHE
    CACHE = a.out / '_cache'; CACHE.mkdir(parents=True, exist_ok=True)
    panel = set(pd.read_csv(V / 'data/pert_counts.csv')['target_gene'].astype(str))
    rows = {}
    for name in a.datasets:
        kind, path, extra = DATASETS[name]
        cache = CACHE / f'{name}.csv'
        fresh = cache.exists() and not a.recompute and (
            kind != 'cells' or (matrix_cache(name).exists() and 'var_cpm' in np.load(matrix_cache(name)).files))
        if fresh:
            s = pd.read_csv(cache, keep_default_na=False, na_values=[''])
        elif not path.exists():
            print(f'{name}: missing {path}, skipped'); continue
        else:
            print(f'{name}: computing from {path}')
            t0 = time.time(); s = BUILDERS[kind](name, path, **extra); s.to_csv(cache, index=False)
            print(f'{name}: {len(s):,} knockdowns in {time.time() - t0:.0f} s')
        s['challenge_target'] = s['target'].astype(str).isin(panel)
        s, keep, failed, has_umi = apply_filters(s, a)
        write_dataset(name, kind, path, extra, s, a)
        kept = s.status == 'keep'
        rows[name] = dict(
            knockdowns=len(s), kept=int(kept.sum()), kept_challenge_targets=int((kept & s.challenge_target).sum()),
            **{f'failed_{k.replace(" ", "_")}': int(v.sum()) for k, v in failed.items()},
            **{f'kept_if_{r}': int(m.sum()) for r, m in keep.items()},
            control_cells=s['n_control_cells'].iloc[0] if s['n_control_cells'].notna().any() else None,
            umi_rule=('umi counts' if has_umi else 'no UMI counts: cpm used'), test=TESTS[kind])
    if not rows: raise SystemExit('no datasets available')
    summary = pd.DataFrame(rows).T.rename_axis('dataset')
    pd.set_option('display.width', 240); pd.set_option('display.max_columns', 30)
    print('\n' + summary.drop(columns=['umi_rule', 'test']).to_string())
    print(f'\nwrote {a.out}/<dataset>.h5ad')


if __name__ == '__main__':
    main()

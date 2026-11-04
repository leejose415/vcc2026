"""Per-source effect tables v4 (lib/effects_v4/), on the challenge's 18,533-gene panel. Never touches lib/effects/.

Decisions: on-target threshold -0.5 and MIN_CELLS = 30 (applied upstream: guide_qc.csv / knockdown_qc.csv);
LFC = log2((knockdown mean per-cell CPM + 0.5) / (control mean per-cell CPM + 0.5)), CPM = counts / n_umi * 1e6, against the
source's own controls (all lanes); SE = delta method from the per-cell spread, Poisson-floored, floored at 0.05; no SE
scaling. Genes a source did not measure are NaN (never 0). The target's own column is kept.

Cell sources (K562, RPE1, HepG2, Jurkat, H1; --cells): from lib/qc_cells/<source>_qc.h5ad plus the QC tables (required).
  1. per construct and gene: n cells, sum CPM, sum CPM^2, total depth (and the same for controls); saved to
     results/qc/<source>/construct_totals.npz and reused when present
  2. kept constructs (guide_qc keep) pooled into kept knockdowns (knockdown_qc knockdown_keep); lfc, se from the pooled totals
Table sources (--table): X-Atlas/Orion HCT116 / HEK293T from per-knockdown sums of per-cell CPM, CPM^2 and counts
(lib/orion/<line>_pseudobulk.h5ad; per-cell CPM over the full library, same formulas as cells) and the Marson CD4 T-cell
tables (lib/tcell/TCELL_<condition>_effects.h5ad: replicate-paired log2 ratio of pseudobulk CPM, pseudocount 0.5; SE from
between-replicate variance moderated toward a Poisson term, limma-style, d0 = 4; not DESeq2, not lfcShrink-ed), with the
knockdown QC of scripts/qc_table_perturbations.py.
n_de_bh: per knockdown, genes significant at BH FDR 5% (two-sided normal p from z = lfc / se), measured genes, target excluded.

Output lib/effects_v4/<source>.h5ad: rows = kept knockdowns, id '<source>|<background>|<condition>|<target_gene_id>';
obs source, background, condition, target_gene, target_gene_id, n_cells, efficiency (target's own lfc in the source's
genes, also when the target is off the panel; NaN unless the target is measurable), efficiency_known, n_de_bh, lfc_type, se_type; var = the panel (index symbol, gene_id = Ensembl ID or
'' when no source maps it); layers 'lfc', 'se' (float32, NaN where unmeasured); no X.
Plus lib/effects_v4/<source>_baseline.csv: per panel gene, control mean CPM and detection rate (NaN if unmeasured).

Usage: ~/venv/bin/python scripts/build_effects_v4.py --cells HepG2 | --table HCT116|HEK293T|CD4_Rest|CD4_Stim8hr|CD4_Stim48hr
"""
import argparse
import json
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).parent))
from qc_replogle_umi import SOURCES                                                        # noqa: E402

V = Path.home() / 'vcc'
OUT = V / 'lib/effects_v4'
PSEUDO, SE_FLOOR, FDR, CHUNK = 0.5, 0.05, 0.05, 50000
CELL_SOURCES = {'K562': 'replogle_k562_gw', 'RPE1': 'replogle_rpe1', 'HepG2': 'nadig_hepg2', 'Jurkat': 'nadig_jurkat', 'H1': 'arc_h1'}
TABLES = {'HCT116': ('orion_hct116', 'HCT116', 'none'), 'HEK293T': ('orion_hek293t', 'HEK293T', 'none'),
          'CD4_Rest': ('marson_cd4_rest', 'CD4_T', 'Rest'), 'CD4_Stim8hr': ('marson_cd4_stim8hr', 'CD4_T', 'Stim8hr'),
          'CD4_Stim48hr': ('marson_cd4_stim48hr', 'CD4_T', 'Stim48hr')}


# ---------------------------------------------------------------- panel
def panel():
    """The 18,533-gene panel: symbol -> Ensembl ID from the cell sources' gene tables (H1 first). '' when none maps it."""
    symbols = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy()
    sym2id, conflicts = {}, []
    for ds in ('H1', 'K562', 'RPE1', 'HepG2', 'Jurkat'):
        p = SOURCES[ds][2]
        if not p.exists(): continue
        with h5py.File(p, 'r') as f: v = read_elem(f['var'])
        pairs = zip(v['gene_name'].astype(str), v.index.astype(str)) if 'gene_name' in v else zip(v.index.astype(str), v['gene_id'].astype(str))
        for s, g in pairs:
            if s in sym2id and sym2id[s] != g: conflicts.append((s, sym2id[s], g, ds))
            sym2id.setdefault(s, g)
    var = pd.DataFrame(dict(gene_id=[sym2id.get(s, '') for s in symbols]), index=pd.Index(symbols, name='gene'))
    return var, conflicts


# ---------------------------------------------------------------- shared maths
def stats(s1, s2, n, depth):
    """Mean per-cell CPM and the relative variance of (mean + 0.5) with the Poisson floor, per row x gene."""
    nn = np.maximum(n, 1)[:, None]; mean = s1 / nn
    var_mean = np.maximum(s2 / nn - mean ** 2, 0) * nn / np.maximum(nn - 1, 1) / nn
    u = 1e6 / np.maximum(depth, 1)[:, None]
    return mean, np.maximum(var_mean, (np.maximum(mean, 0) + u) * u) / (mean + PSEUDO) ** 2


def bh_count(lfc, se, exclude):
    """Per row: measured genes (target excluded) significant at BH FDR 5%, two-sided normal p from lfc / se."""
    out = np.zeros(len(lfc), int)
    for i in range(len(lfc)):
        ok = np.isfinite(lfc[i]) & np.isfinite(se[i])
        if exclude[i] >= 0: ok[exclude[i]] = False
        p = 2 * norm.sf(np.abs(lfc[i, ok] / se[i, ok]))
        if not len(p): continue
        q = np.sort(p) * len(p) / np.arange(1, len(p) + 1)
        passed = np.flatnonzero(q <= FDR)
        out[i] = passed[-1] + 1 if len(passed) else 0
    return out


def write(source, background, condition, var, kd, lfc, se, baseline, meta):
    """kd: DataFrame with target_gene, target_gene_id, n_cells, efficiency, efficiency_known, n_de_bh, lfc_type, se_type."""
    obs = kd.copy(); obs.insert(0, 'condition', condition); obs.insert(0, 'background', background); obs.insert(0, 'source', source)
    obs.index = pd.Index([f'{source}|{background}|{condition}|{t}' for t in obs.target_gene_id], name='row_id')
    for c in ('source', 'background', 'condition', 'lfc_type', 'se_type'): obs[c] = obs[c].astype('category')
    a = ad.AnnData(obs=obs, var=var.copy(), layers=dict(lfc=lfc.astype(np.float32), se=se.astype(np.float32)))
    a.uns['effects_v4'] = dict(meta, built=time.strftime('%Y-%m-%d %H:%M'), se_floor=SE_FLOOR, pseudocount=PSEUDO, fdr=FDR)
    OUT.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(OUT / f'{source}.h5ad', compression='gzip')
    baseline.to_csv(OUT / f'{source}_baseline.csv')
    print(f'wrote {OUT / f"{source}.h5ad"}: {a.n_obs:,} knockdowns x {a.n_vars:,} panel genes, '
          f'{int(np.isfinite(lfc).any(0).sum()):,} measured', flush=True)
    return a


# ---------------------------------------------------------------- cell sources
def construct_totals(dataset, guides):
    """Per construct (every one in guide_qc) and controls: n, sum CPM, sum CPM^2, depth, control detection; cached."""
    _, qc_dir, path, _ = SOURCES[dataset]
    cache = qc_dir / 'construct_totals.npz'
    if cache.exists():
        z = np.load(cache, allow_pickle=False)
        if int(z['mtime_ns']) == path.stat().st_mtime_ns and list(z['guides'].astype(str)) == list(guides):
            print(f'  reusing {cache}', flush=True); return {k: z[k] for k in z.files}
    gidx = {g: i for i, g in enumerate(guides)}; t0 = time.time()
    with h5py.File(path, 'r') as f:
        obs = read_elem(f['obs']); X = f['X']; n_obs = len(obs); n_vars = int(X.attrs['shape'][1])
        ctrl = obs['is_control'].to_numpy(bool); lib = obs['n_umi'].to_numpy(float)
        grp = np.where(ctrl, len(guides), [gidx.get(g, -1) for g in obs['guide_id'].astype(str)])
        k = len(guides) + 1
        s1 = np.zeros((k, n_vars)); s2 = np.zeros((k, n_vars)); n = np.zeros(k); depth = np.zeros(k); detect = np.zeros(n_vars)
        ip = X['indptr']
        for s in range(0, n_obs, CHUNK):
            e = min(s + CHUNK, n_obs); g = grp[s:e]; r = np.flatnonzero(g >= 0)
            if len(r):
                p = ip[s:e + 1]; a, b = int(p[0]), int(p[-1])
                M = sp.csr_matrix((X['data'][a:b], X['indices'][a:b], p - a), shape=(e - s, n_vars))[r]
                c = g[r] == len(guides)
                if c.any(): detect += np.asarray((M[np.flatnonzero(c)] > 0).sum(0)).ravel()
                M = sp.diags(1e6 / lib[s:e][r]) @ M.astype(np.float64)
                G = sp.csr_matrix((np.ones(len(r)), (g[r], np.arange(len(r)))), shape=(k, len(r)))
                s1 += (G @ M).toarray(); s2 += (G @ M.multiply(M)).toarray()
                n += np.bincount(g[r], minlength=k); depth += np.bincount(g[r], weights=lib[s:e][r], minlength=k)
            print(f'  totals: {e:,} / {n_obs:,} cells ({time.time() - t0:.0f}s)', flush=True)
    tot = dict(guides=np.array(guides, dtype=str), s1=s1, s2=s2, n=n, depth=depth, detect=detect, mtime_ns=np.int64(path.stat().st_mtime_ns))
    np.savez(cache, **tot)
    return tot


def build_cells(dataset):
    source = CELL_SOURCES[dataset]; _, qc_dir, path, title = SOURCES[dataset]
    gq, kq = qc_dir / 'guide_qc.csv', qc_dir / 'knockdown_qc.csv'
    missing = [str(p) for p in (gq, kq, path) if not p.exists()]
    if missing: sys.exit(f'error: missing input(s): {missing}')
    guides = pd.read_csv(gq); kds = pd.read_csv(kq)
    guides['guide_id'] = guides['guide_id'].astype(str); guides['target_gene_id'] = guides['target_gene_id'].astype(str)
    kds['target_gene_id'] = kds['target_gene_id'].astype(str)
    tot = construct_totals(dataset, guides.guide_id.tolist())
    with h5py.File(path, 'r') as f: var_src = read_elem(f['var'])
    src_ids = var_src.index.astype(str).to_numpy() if 'gene_name' in var_src else var_src['gene_id'].astype(str).to_numpy()
    var, _ = panel(); pid = {g: i for i, g in enumerate(var.gene_id) if g}
    to_panel = np.array([pid.get(g, -1) for g in src_ids]); ok = to_panel >= 0
    kept = kds[kds.knockdown_keep.astype(bool)].reset_index(drop=True)
    kg = guides[guides.keep.astype(bool) & guides.target_gene_id.isin(set(kept.target_gene_id))]
    row_of = {t: i for i, t in enumerate(kept.target_gene_id)}
    P = sp.csr_matrix((np.ones(len(kg)), ([row_of[t] for t in kg.target_gene_id], kg.index.to_numpy())), shape=(len(kept), len(guides)))
    s1, s2, n, depth = P @ tot['s1'][:-1], P @ tot['s2'][:-1], P @ tot['n'][:-1], P @ tot['depth'][:-1]
    m1, rv1 = stats(s1, s2, n, depth)
    m0, rv0 = stats(tot['s1'][-1:], tot['s2'][-1:], tot['n'][-1:], tot['depth'][-1:])
    lfc_src = np.log2((m1 + PSEUDO) / (m0 + PSEUDO)); se_src = np.maximum(np.sqrt(rv1 + rv0) / np.log(2), SE_FLOOR)
    G = len(var); lfc = np.full((len(kept), G), np.nan); se = np.full_like(lfc, np.nan)
    lfc[:, to_panel[ok]] = lfc_src[:, ok]; se[:, to_panel[ok]] = se_src[:, ok]
    target_col = np.array([pid.get(t, -1) for t in kept.target_gene_id])
    src_col = {g: i for i, g in enumerate(src_ids)}
    target_src = np.array([src_col.get(t, -1) for t in kept.target_gene_id])           # the target in the source's own genes
    measurable = kept.measurable.astype(bool).to_numpy()                               # (it may be off the panel)
    eff = np.array([lfc_src[i, j] if (j >= 0 and measurable[i]) else np.nan for i, j in enumerate(target_src)])
    kd = pd.DataFrame(dict(target_gene=kept.target_gene.astype(str), target_gene_id=kept.target_gene_id, n_cells=n.astype(int),
                           efficiency=eff, efficiency_known=measurable, n_de_bh=bh_count(lfc, se, target_col),
                           lfc_type='mean_cpm_ratio', se_type='delta_cells'))
    base = pd.DataFrame(dict(control_mean_cpm=np.nan, detection_rate=np.nan), index=var.index)
    base.iloc[to_panel[ok], 0] = m0[0, ok]; base.iloc[to_panel[ok], 1] = (tot['detect'] / tot['n'][-1])[ok]
    meta = dict(cell_file=str(path), guide_qc=str(gq), knockdown_qc=str(kq), control_cells=int(tot['n'][-1]),
                source_genes=int(len(src_ids)), source_genes_on_panel=int(ok.sum()),
                note='H1: training table only; it is excluded from training in validation runs and the scoring harness keeps the unfiltered real H1 cells' if dataset == 'H1' else '')
    a = write(source, dataset, 'none', var, kd, lfc, se, base, meta)
    checks = cell_checks(dataset, a, kept, kg, kds, src_ids, pid)
    return a, checks


def cell_checks(dataset, a, kept, kg, kds, src_ids, pid, n_kd=5, n_genes=20):
    """3: 5 random knockdowns x 20 random measured genes, lfc recomputed straight from the filtered cells (1e-6)."""
    _, qc_dir, path, _ = SOURCES[dataset]
    rng = np.random.default_rng(2026)
    pick = rng.choice(len(kept), min(n_kd, len(kept)), replace=False)
    measured_src = np.flatnonzero(np.array([pid.get(g, -1) for g in src_ids]) >= 0)
    genes = rng.choice(measured_src, n_genes, replace=False)
    want = {kept.target_gene_id[i]: set(kg.loc[kg.target_gene_id == kept.target_gene_id[i], 'guide_id']) for i in pick}
    sums = {t: [0.0 * genes, 0] for t in want}; csum = np.zeros(n_genes); cn = 0
    with h5py.File(path, 'r') as f:
        obs = read_elem(f['obs']); X = f['X']; ip = X['indptr']; n_vars = int(X.attrs['shape'][1])
        ctrl = obs['is_control'].to_numpy(bool); lib = obs['n_umi'].to_numpy(float)
        guide = obs['guide_id'].astype(str).to_numpy(); tid = obs['target_gene_id'].astype(str).to_numpy()
        for s in range(0, len(obs), CHUNK):
            e = min(s + CHUNK, len(obs)); p = ip[s:e + 1]; a0, b0 = int(p[0]), int(p[-1])
            M = sp.csr_matrix((X['data'][a0:b0], X['indices'][a0:b0], p - a0), shape=(e - s, n_vars))[:, genes].toarray() * (1e6 / lib[s:e])[:, None]
            cm = ctrl[s:e]; csum += M[cm].sum(0); cn += int(cm.sum())
            for t, gs in want.items():
                m = (tid[s:e] == t) & np.isin(guide[s:e], list(gs))
                if m.any(): sums[t][0] = sums[t][0] + M[m].sum(0); sums[t][1] += int(m.sum())
    worst = 0.0
    for t in want:
        direct = np.log2((sums[t][0] / sums[t][1] + PSEUDO) / (csum / cn + PSEUDO))
        cols = [pid[src_ids[g]] for g in genes]
        table = a.layers['lfc'][a.obs.target_gene_id.astype(str).to_numpy() == t][0, cols].astype(float)
        worst = max(worst, float(np.nanmax(np.abs(direct - table))))
    eff_qc = kds.set_index('target_gene_id').loc[a.obs.target_gene_id.astype(str), 'efficiency'].to_numpy(float)
    eff_tab = a.obs['efficiency'].to_numpy(float)
    both = np.isfinite(eff_qc) | np.isfinite(eff_tab)
    return dict(check3_max_abs_diff=worst, check3_pass=worst < 1e-6,
                check4_max_abs_diff=float(np.nanmax(np.abs(eff_qc[both] - eff_tab[both]))) if both.any() else 0.0,
                check4_same_nan=bool((np.isfinite(eff_qc) == np.isfinite(eff_tab)).all()))


# ---------------------------------------------------------------- table sources
def build_table(name):
    source, background, condition = TABLES[name]
    var, _ = panel(); sym_col = {s: i for i, s in enumerate(var.index)}
    qc_name = {'HCT116': 'orion_hct116_table', 'HEK293T': 'orion_hek293t_table'}.get(name, source)
    kq = V / 'results/qc' / qc_name / 'knockdown_qc.csv'
    if not kq.exists(): sys.exit(f'error: missing {kq}; run scripts/qc_table_perturbations.py first')
    qc = pd.read_csv(kq); qc['target_gene'] = qc['target_gene'].astype(str)
    kept = qc[qc.knockdown_keep.astype(bool)].reset_index(drop=True)
    G = len(var); lfc = np.full((len(kept), G), np.nan); se = np.full_like(lfc, np.nan)
    if name in ('HCT116', 'HEK293T'):
        path = V / f'lib/orion/{name}_pseudobulk.h5ad'
        with h5py.File(path, 'r') as f:
            obs = read_elem(f['obs']); vv = read_elem(f['var']); names = obs.index.astype(str); pos = {t: i for i, t in enumerate(names)}
            ncell = obs['n_cells'].to_numpy(float); cols = np.array([sym_col.get(s, -1) for s in vv.index.astype(str)])
            meas = vv['measured'].to_numpy(bool) & (cols >= 0)
            L = {k: f['layers'][k] for k in ('sum_cpm', 'sum_cpm2', 'sum_counts', 'n_expr')}
            def moments(rows):
                rows = np.asarray(rows); order = np.argsort(rows); r = rows[order]
                vals = {k: L[k][r].astype(float) for k in L}; back = np.argsort(order)
                vals = {k: v[back] for k, v in vals.items()}
                depth = vals['sum_counts'].sum(1)
                return stats(vals['sum_cpm'], vals['sum_cpm2'], ncell[rows], depth) + (vals,)
            m0, rv0, c0 = moments([pos['non-targeting']])
            for s in range(0, len(kept), 256):
                rows = [pos[t] for t in kept.target_gene[s:s + 256]]
                m1, rv1, _ = moments(rows)
                l = np.log2((m1 + PSEUDO) / (m0 + PSEUDO)); e = np.maximum(np.sqrt(rv1 + rv0) / np.log(2), SE_FLOOR)
                lfc[s:s + len(rows)][:, cols[meas]] = l[:, meas]; se[s:s + len(rows)][:, cols[meas]] = e[:, meas]
        base = pd.DataFrame(dict(control_mean_cpm=np.nan, detection_rate=np.nan), index=var.index)
        nc = ncell[pos['non-targeting']]
        base.iloc[cols[meas], 0] = (c0['sum_cpm'][0] / nc)[meas]; base.iloc[cols[meas], 1] = (c0['n_expr'][0] / nc)[meas]
        lfc_type, se_type = 'mean_cpm_ratio', 'orion_cell_moments_delta'
        meta = dict(table=str(path), method='per-knockdown sums of per-cell CPM (full-library depth), CPM^2, counts, cells expressing; '
                    'lfc and se with the same formulas as the cell sources', control_cells=int(nc))
    else:
        cond = TABLES[name][2]; path = V / f'lib/tcell/TCELL_{cond}_effects.h5ad'
        with h5py.File(path, 'r') as f:
            obs = read_elem(f['obs']); vv = read_elem(f['var']); names = obs.index.astype(str); pos = {t: i for i, t in enumerate(names)}
            cols = np.array([sym_col.get(s, -1) for s in vv.index.astype(str)]); meas = vv['measured'].to_numpy(bool) & (cols >= 0)
            rows = np.array([pos[t] for t in kept.target_gene]); order = np.argsort(rows)
            Xl = f['X'][np.sort(rows)][np.argsort(order)]; Xs = f['layers/se'][np.sort(rows)][np.argsort(order)]
        lfc[:, cols[meas]] = Xl[:, meas]; se[:, cols[meas]] = np.maximum(Xs[:, meas], SE_FLOOR)
        base = pd.DataFrame(dict(control_mean_cpm=np.nan, detection_rate=np.nan), index=var.index)
        base.iloc[cols[meas], 0] = vv['ctrl_cpm'].to_numpy(float)[meas]
        lfc_type, se_type = 'sum_ratio', 'donor_moderated'
        meta = dict(table=str(path), method='replicate-paired log2 ratio of pseudobulk CPM (pseudocount 0.5); SE from between-replicate '
                    'variance moderated toward a Poisson term (limma-style, d0 = 4, source floor 0.02, here floored at 0.05); not DESeq2, '
                    'not shrunk', detection_rate='not available (no per-cell data)')
    tid = [var.gene_id.iloc[sym_col[t]] if t in sym_col and var.gene_id.iloc[sym_col[t]] else f'symbol:{t}' for t in kept.target_gene]
    target_col = np.array([sym_col.get(t, -1) for t in kept.target_gene])
    measurable = kept.measurable.astype(bool).to_numpy()
    eff = np.array([lfc[i, j] if (j >= 0 and measurable[i]) else np.nan for i, j in enumerate(target_col)])
    kd = pd.DataFrame(dict(target_gene=kept.target_gene, target_gene_id=tid, n_cells=kept.n_cells.astype(float), efficiency=eff,
                           efficiency_known=measurable, n_de_bh=bh_count(lfc, se, target_col), lfc_type=lfc_type, se_type=se_type))
    meta.update(knockdown_qc=str(kq))
    a = write(source, background, condition, var, kd, lfc, se, base, meta)
    eff_qc = kept.ontarget_lfc.to_numpy(float, copy=True); eff_qc[~measurable] = np.nan
    both = np.isfinite(eff_qc) | np.isfinite(eff)
    return a, dict(check4_max_abs_diff=float(np.nanmax(np.abs(eff_qc[both] - eff[both]))) if both.any() else 0.0,
                   check4_same_nan=bool((np.isfinite(eff_qc) == np.isfinite(eff)).all()), knockdown_qc_kept=int(len(kept)))


def common_checks(a, expected_rows):
    lfc, se = a.layers['lfc'], a.layers['se']
    return dict(check1_rows=int(a.n_obs), check1_expected=int(expected_rows), check1_pass=int(a.n_obs) == int(expected_rows),
                check2_all_nan_rows=int((~np.isfinite(lfc).any(1)).sum()), check2_nan_mismatch=int((np.isfinite(lfc) != np.isfinite(se)).sum()))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--cells', choices=list(CELL_SOURCES)); g.add_argument('--table', choices=list(TABLES))
    args = p.parse_args(); t0 = time.time()
    if args.cells:
        a, checks = build_cells(args.cells)
        kq = pd.read_csv(SOURCES[args.cells][1] / 'knockdown_qc.csv'); expected = int(kq.knockdown_keep.sum()); name = CELL_SOURCES[args.cells]
    else:
        a, checks = build_table(args.table); expected = checks['knockdown_qc_kept']; name = TABLES[args.table][0]
    checks.update(common_checks(a, expected))
    summary = dict(source=name, rows=int(a.n_obs), measured_genes=int(np.isfinite(a.layers['lfc']).any(0).sum()),
                   median_efficiency=float(np.nanmedian(a.obs.efficiency)) if a.obs.efficiency.notna().any() else None,
                   n_de_bh_median=float(a.obs.n_de_bh.median()), n_de_bh_p90=float(a.obs.n_de_bh.quantile(0.9)),
                   checks=checks, minutes=round((time.time() - t0) / 60, 1))
    (OUT / f'{name}_build.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))

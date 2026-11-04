"""v5 effect tables for the rest of X-Atlas/Orion: every knockdown the table QC keeps (results/qc/orion_<line>_table/
knockdown_qc.csv) that the single-cell v5 table (lib/effects_v5/orion_<line>.h5ad, the A-C targets subset) does not hold,
from the per-knockdown pseudobulk sums (lib/orion/<LINE>_pseudobulk.h5ad: n_cells, sum_counts, sum_cpm, sum_cpm2, n_expr on
the panel genes). Output lib/effects_v5/orion_<line>_pseudobulk.h5ad (+ _baseline.csv, _build.json), v5 schema, written
by build_effects_v5.write (imported, not modified).

  LFC      log2((m1 + 1e-9) / (m0 + 1e-9)), mean CPM vs non-targeting, v5's formula. The source's CPM is over each cell's
           whole library (~29% of counts are off the panel); v5's is over the panel genes, so each group's mean is divided by
           its mean panel share s = sum_g(sum_cpm) / (n * 1e6) -- exact when every cell of the group has the same share.
  SE       delta method from sum_cpm2 (sample variance of per-cell CPM, same rescaling), n capped at the scorer's 400
           knockdown / 18,400 control cells. Full-depth variance: sums cannot be thinned to 20k UMI.
  de_sig   none: no per-cell values, no Wilcoxon (significance_unknown = True, de_labels_known = False, n_real NaN).
  profiles log1p_bulk_50k_full from sum_counts over the measured panel genes; log1p_bulk_50k_matched NaN (no sampling).
  The pseudobulk carries no cell QC: all cells of each knockdown and of non-targeting are in the sums.
Check (always run): the knockdowns held as single cells (lib/qc_cells/orion_<line>_qc.h5ad, the A-C subset; QC'd cells,
full depth, panel library -- the exact v5 quantities) vs this approximation on the same knockdowns, with and without the
panel rescaling: control mean CPM and LFC on the gated genes.

Usage: ~/venv/bin/python scripts/build_effects_v5_orion.py --lines HCT116 HEK293T
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

sys.path.insert(0, str(Path(__file__).parent))
import build_effects_v5 as B                                                             # noqa: E402

V = Path.home() / 'vcc'
LFC_TYPE, SE_TYPE = 'official_cpm_1e-9_pseudobulk_panel_rescaled', 'pseudobulk_fulldepth_delta_n400cap'
CHUNK = 2000


def group_stats(f, rows, n):
    """Per row of the pseudobulk: panel-rescaled mean and variance of per-cell CPM, raw count sums, detection counts,
    and the mean panel share s; plus the unrescaled (whole-library) mean, for the check."""
    order = np.argsort(rows); r = np.asarray(rows)[order]
    L = {k: f['layers'][k][r.tolist(), :].astype(np.float64) for k in ('sum_cpm', 'sum_cpm2', 'sum_counts', 'n_expr')}
    inv = np.empty_like(order); inv[order] = np.arange(len(order))
    L = {k: v[inv] for k, v in L.items()}
    n = np.asarray(n, float)[:, None]
    s = L['sum_cpm'].sum(1, keepdims=True) / (n * 1e6)
    m_full = L['sum_cpm'] / n
    var_full = np.maximum(L['sum_cpm2'] - L['sum_cpm'] ** 2 / n, 0) / np.maximum(n - 1, 1)
    return dict(m=m_full / s, v=var_full / s ** 2, m_full=m_full, counts=L['sum_counts'], n_expr=L['n_expr'], s=s.ravel())


def check(line, f, idx, meas, var, gated_cols):
    """Exact v5 quantities from the A-C single cells vs the pseudobulk approximation, same knockdowns."""
    p = V / f'lib/qc_cells/orion_{line.lower()}_qc.h5ad'
    with h5py.File(p, 'r') as h:
        o = read_elem(h['obs']); X = sp.csr_matrix((h['X/data'][:], h['X/indices'][:], h['X/indptr'][:]),
                                                    shape=(len(o), h['X'].attrs['shape'][1]))
    lib = np.asarray(X.sum(1)).ravel()
    cpm = sp.diags(1e6 / np.where(lib > 0, lib, 1)) @ X
    ctrl = o.is_control.to_numpy(bool); tg = o.target_gene.astype(str).to_numpy()
    targets = sorted(t for t in set(tg[~ctrl]) if t in idx)
    m0_sc = np.asarray(cpm[ctrl].mean(0)).ravel()
    pb = group_stats(f, [idx['non-targeting']] + [idx[t] for t in targets],
                     f['obs']['n_cells'][:][[idx['non-targeting']] + [idx[t] for t in targets]])
    g = gated_cols
    out = dict(knockdowns=len(targets), single_cell_controls=int(ctrl.sum()), gated_genes=int(len(g)))
    for tag, M in (('rescaled', pb['m']), ('whole_library', pb['m_full'])):
        r0 = np.log2((M[0, g] + B.EPS) / (m0_sc[g] + B.EPS))
        d, rr = [], []
        for j, t in enumerate(targets):
            m1_sc = np.asarray(cpm[tg == t].mean(0)).ravel()
            l_sc = np.log2((m1_sc[g] + B.EPS) / (m0_sc[g] + B.EPS)); l_pb = np.log2((M[j + 1, g] + B.EPS) / (M[0, g] + B.EPS))
            ok = (m1_sc[g] > 0) & (M[j + 1, g] > 0)
            d.append(np.median(np.abs(l_pb[ok] - l_sc[ok]))); rr.append(np.corrcoef(l_pb[ok], l_sc[ok])[0, 1])
        out[tag] = dict(control_mean_cpm_median_abs_log2_ratio=float(np.median(np.abs(r0))),
                        control_mean_cpm_median_log2_ratio=float(np.median(r0)),
                        lfc_median_abs_diff=float(np.median(d)), lfc_p90_abs_diff=float(np.quantile(d, 0.9)),
                        lfc_pearson_median=float(np.median(rr)), lfc_pearson_p10=float(np.quantile(rr, 0.1)))
    out['note'] = ('single cells: the QC\'d A-C subset (cell QC applied; its own control subset); pseudobulk: all cells. '
                   'Differences mix the panel rescaling with cell QC and control sampling.')
    return out


def build(line):
    t0 = time.time(); l = line.lower(); name = f'orion_{l}_pseudobulk'
    if (B.OUT / f'{name}.h5ad').exists(): sys.exit(f'{B.OUT / name}.h5ad exists; refusing to overwrite')
    sc = ad.read_h5ad(B.OUT / f'orion_{l}.h5ad', backed='r'); var = sc.var.to_memory() if hasattr(sc.var, 'to_memory') else sc.var.copy()
    held = set(sc.obs.target_gene.astype(str)); sc.file.close()
    kq_path = V / f'results/qc/orion_{l}_table/knockdown_qc.csv'
    kq = pd.read_csv(kq_path); kq['target_gene'] = kq.target_gene.astype(str)
    keep = kq.knockdown_keep.astype(bool)
    kept = kq[keep & ~kq.target_gene.isin(held) & (kq.target_gene != 'non-targeting')].reset_index(drop=True)
    path = V / f'lib/orion/{line}_pseudobulk.h5ad'
    with h5py.File(path, 'r') as f:
        o = read_elem(f['obs']); pv = read_elem(f['var'])
        assert np.array_equal(pv.index.astype(str).to_numpy(), var.index.astype(str).to_numpy()), 'pseudobulk genes != panel'
        meas = pv['measured'].to_numpy(bool); idx = {k: i for i, k in enumerate(o.index.astype(str))}
        missing = [t for t in kept.target_gene if t not in idx]
        if missing: sys.exit(f'{line}: {len(missing)} kept knockdowns missing from the pseudobulk: {missing[:5]}')
        n_all = o.n_cells.to_numpy(float); ci = idx['non-targeting']; n0 = n_all[ci]
        c = group_stats(f, [ci], [n0]); m0, v0 = c['m'][0], c['v'][0]
        gated = np.flatnonzero(meas & (m0 > B.MIN_CPM))
        B.log(f'{name}: {len(kept):,} knockdowns (table-QC kept {int(keep.sum()):,}, {len(held)} held as single cells), '
              f'{int(n0):,} non-targeting cells (panel share {c["s"][0]:.3f}), {int(meas.sum()):,} measured genes, '
              f'{len(gated):,} pass the {B.MIN_CPM:g}-CPM gate')
        G = len(var); nk = len(kept)
        lfc = np.full((nk, G), np.nan, np.float32); se = np.full_like(lfc, np.nan)
        bulk_full = np.full_like(lfc, np.nan); bulk_matched = np.full_like(lfc, np.nan)
        n1 = n_all[[idx[t] for t in kept.target_gene]]; share = np.zeros(nk)
        n0c = min(n0, B.N_CTRL)
        for s0 in range(0, nk, CHUNK):
            rows = [idx[t] for t in kept.target_gene[s0:s0 + CHUNK]]; st = group_stats(f, rows, n1[s0:s0 + CHUNK])
            m1, v1 = st['m'][:, meas], st['v'][:, meas]; n1c = np.minimum(n1[s0:s0 + CHUNK], B.N_CELLS)[:, None]
            l_ = np.log2((m1 + B.EPS) / (m0[meas] + B.EPS))
            e_ = np.sqrt(v1 / n1c / (m1 + B.EPS) ** 2 + v0[meas] / n0c / (m0[meas] + B.EPS) ** 2) / np.log(2)
            e_[(m1 == 0) | (m0[meas] == 0)[None, :]] = np.inf
            e_[n1[s0:s0 + CHUNK] < 2] = np.nan
            lfc[s0:s0 + CHUNK][:, meas] = l_; se[s0:s0 + CHUNK][:, meas] = e_
            cnt = st['counts'][:, meas]; tot = cnt.sum(1, keepdims=True)
            bulk_full[s0:s0 + CHUNK][:, meas] = np.log1p(B.BULK_TS * cnt / np.where(tot > 0, tot, 1))
            share[s0:s0 + CHUNK] = st['s']
            B.log(f'  {min(s0 + CHUNK, nk):,} / {nk:,} knockdowns')
        cnt0 = c['counts'][0, meas]
        chk = check(line, f, idx, meas, var, gated)
    col = {s: j for j, s in enumerate(var.index.astype(str))}; gid = var.gene_id.astype(str).to_numpy()
    tcol = np.array([col.get(t, -1) for t in kept.target_gene])
    tid = np.array([gid[j] if j >= 0 and gid[j] else f'unresolved:{t}' for t, j in zip(kept.target_gene, tcol)])
    if len(set(tid)) < len(tid):                                  # two symbols on one ID: keep both rows apart
        dup = pd.Series(tid).duplicated(keep=False).to_numpy()
        tid = np.where(dup, [f'{i}|{t}' for i, t in zip(tid, kept.target_gene)], tid)
    measurable = kept.measurable.astype(bool).to_numpy()
    eff = np.array([float(lfc[k, j]) if j >= 0 and measurable[k] and meas[j] else np.nan for k, j in enumerate(tcol)])
    obs = pd.DataFrame(dict(target_gene=kept.target_gene.to_numpy(), target_gene_id=tid, n_cells_kept=n1.astype(int),
                            n_cells_used=n1.astype(int), efficiency=eff, efficiency_known=measurable, n_real=np.nan,
                            de_labels_known=False, modality='pseudobulk_table', significance_unknown=True,
                            lfc_type=LFC_TYPE, se_type=SE_TYPE, reference='non-targeting'))
    base = pd.DataFrame(dict(control_mean_cpm=np.nan, detection_rate=np.nan, control_mean_cpm_fulldepth=np.nan,
                             log1p_bulk_50k_matched=np.nan, log1p_bulk_50k_full=np.nan, de_tested=False), index=var.index)
    mc = np.flatnonzero(meas)
    base.iloc[mc, 0] = m0[meas]; base.iloc[mc, 1] = c['n_expr'][0, meas] / n0; base.iloc[mc, 2] = m0[meas]
    base.iloc[mc, 4] = np.log1p(B.BULK_TS * cnt0 / max(cnt0.sum(), 1)); base.iloc[gated, 5] = True
    sig = sp.csr_matrix((nk, G), dtype=bool)
    meta = dict(cell_file=str(path), knockdown_qc=str(kq_path), guide_qc='none (Orion: pseudobulk table, table knockdown QC)',
                method='pseudobulk sums -> v5 LFC formula; CPM rescaled from whole-library to panel library by each group\'s '
                       'mean panel share; SE by the delta method from sum_cpm2, n capped at 400 / 18,400; no Wilcoxon',
                rows_excluded_held_as_single_cells=len(held & set(kq.target_gene[keep])), control_cells=int(n0),
                control_cells_used=int(n0), control_panel_share=float(c['s'][0]),
                knockdown_panel_share_median=float(np.median(share)), source_genes_on_panel=int(meas.sum()),
                de_tested_genes=int(len(gated)), de_labels='none', modality='pseudobulk_table', significance_unknown=True,
                reference='non-targeting', control_quantiles='not stored (no single cells)', thinning_factor='none (sums cannot be thinned)',
                note='Orion: the table-QC kept knockdowns not held as single cells; companion of orion_<line>.h5ad')
    a = B.write(name, line, 'none', var, obs, lfc, se, bulk_matched, bulk_full, sig, base, meta)
    summary = dict(source=name, rows=int(a.n_obs), measured_genes=int(meas.sum()), de_tested_genes=int(len(gated)),
                   median_efficiency=float(np.nanmedian(eff)), median_n_cells=float(np.median(n1)),
                   lfc_on_epsilon_floor_share=float(np.mean(np.abs(lfc[np.isfinite(lfc)]) > 25)),
                   control_panel_share=float(c['s'][0]), knockdown_panel_share_median=float(np.median(share)),
                   check=chk, minutes=round((time.time() - t0) / 60, 1))
    (B.OUT / f'{name}_build.json').write_text(json.dumps(summary, indent=1, default=str) + '\n')
    print(json.dumps(summary, indent=1, default=str), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--lines', nargs='+', choices=['HCT116', 'HEK293T'], default=['HCT116', 'HEK293T'])
    for line in p.parse_args().lines: build(line)

"""Baselines v2: control features on one shared library definition (lib/effects_v4/*_baseline_v2.csv). Never overwrites
older baselines and never touches lfc / se.

S = the panel genes measured by every effects_v4 source (6,013). Per control cell, library_S = sum of counts over S
and CPM_S = counts / library_S * 1e6 (for every measured panel gene). Per gene:
  mean_cpm_S, log1p_cpm_S = log1p(mean CPM_S), detection_rate (share of control cells with count > 0),
  pct_rank_S = percentile rank of mean CPM_S among the S genes of that line (0-1; S genes rank among themselves, other
  genes are placed on the S scale), in_S, measured. NaN where the line did not measure the gene.
Lines with cells: K562, RPE1, HepG2, Jurkat, H1 (lib/qc_cells), HCT116 and HEK293T (their per-cell controls in
lib/qc_cells), challenge lines A-C (data/context_<X>.h5ad). CD4 T (Rest, Stim8hr, Stim48hr) have no cells: mean CPM_S is
the pseudobulk control CPM rescaled so the S genes sum to 1e6 (a ratio of sums); detection_rate NaN, detection_missing 1.
Also writes lib/effects_v4/baseline_v2.json: per line cells, median library_S, the median v2 / panel-only CPM ratio.
Usage: python scripts/build_baseline_v2.py [--lines K562 RPE1 ...]
"""
import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

sys.path.insert(0, str(Path(__file__).parent))
from build_effects_v4 import CELL_SOURCES, OUT, TABLES, V, panel             # noqa: E402
from qc_replogle_umi import SOURCES                                          # noqa: E402

CHUNK = 20_000


def shared_genes(genes):
    """S: measured (finite control CPM in the current baseline) by every effects_v4 source."""
    ok = np.ones(len(genes), bool)
    for f in sorted(OUT.glob('*_baseline.csv')):
        ok &= pd.read_csv(f, index_col=0)['control_mean_cpm'].reindex(genes).notna().to_numpy()
    return ok


def cell_baseline(path, to_panel, ctrl, in_S):
    """Mean per-cell CPM_S and detection over control cells; to_panel maps source columns to panel columns (-1 = off)."""
    G = len(in_S); ok = to_panel >= 0; lib_cols = ok & in_S[np.where(ok, to_panel, 0)]
    keep = sp.diags(ok.astype(float)).tocsr(); libsel = sp.diags(lib_cols.astype(float)).tocsr()
    s1 = np.zeros(G); det = np.zeros(G); n = 0; libs = []
    with h5py.File(path, 'r') as f:
        X = f['X']; ip = X['indptr']; n_obs, n_vars = (int(x) for x in X.attrs['shape'])
        for s in range(0, n_obs, CHUNK):
            e = min(s + CHUNK, n_obs); r = np.flatnonzero(ctrl[s:e])
            if not len(r): continue
            p = ip[s:e + 1]; a, b = int(p[0]), int(p[-1])
            M = sp.csr_matrix((X['data'][a:b], X['indices'][a:b], p - a), shape=(e - s, n_vars))[r].astype(np.float64)
            lib = np.asarray((M @ libsel).sum(1)).ravel(); libs.append(lib)
            C = sp.diags(1e6 / np.maximum(lib, 1)) @ (M @ keep)
            np.add.at(s1, to_panel[ok], np.asarray(C.sum(0)).ravel()[ok])
            np.add.at(det, to_panel[ok], np.asarray((M > 0).sum(0)).ravel()[ok]); n += len(r)
    mean = np.full(G, np.nan); rate = np.full(G, np.nan); cols = np.unique(to_panel[ok])
    mean[cols] = s1[cols] / n; rate[cols] = det[cols] / n
    libs = np.concatenate(libs)
    return mean, rate, dict(control_cells=n, library_S_genes=int(lib_cols.sum()), median_library_S=float(np.median(libs)),
                            zero_library_cells=int((libs == 0).sum()))


def features(genes, mean, rate, in_S, missing_detection=False):
    measured = np.isfinite(mean)
    ref = np.sort(mean[in_S & measured])
    pct = np.full(len(genes), np.nan)
    pct[measured] = np.searchsorted(ref, mean[measured], side='right') / len(ref)
    return pd.DataFrame(dict(mean_cpm_S=mean, log1p_cpm_S=np.log1p(mean), detection_rate=rate,
                             detection_missing=int(missing_detection), pct_rank_S=pct, in_S=in_S.astype(int),
                             measured=measured.astype(int)), index=pd.Index(genes, name='gene'))


def save(name, df, info, summary):
    path = OUT / f'{name}_baseline_v2.csv'; df.to_csv(path)
    old = OUT / f'{name}_baseline_panel.csv'
    if old.exists():
        o = pd.read_csv(old, index_col=0)['control_mean_cpm'].reindex(df.index).to_numpy(float)
        m = df['mean_cpm_S'].to_numpy(float); both = np.isfinite(o) & np.isfinite(m) & (o > 0)
        info['median_ratio_v2_over_panel'] = float(np.median(m[both] / o[both]))
    info['detection_rate_present'] = bool(df['detection_rate'].notna().any())
    info['measured_genes'] = int(df['measured'].sum())
    summary[name] = info
    print(f'{name:22s} {json.dumps(info)}', flush=True)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--lines', nargs='*'); a = ap.parse_args()
    var, _ = panel(); genes = var.index.to_numpy(); G = len(genes)
    pid = {g: i for i, g in enumerate(var.gene_id) if g}; sym = {s: i for i, s in enumerate(genes)}
    in_S = shared_genes(genes); print(f'S = {int(in_S.sum())} genes', flush=True)
    want = set(a.lines) if a.lines else None
    summ_path = OUT / 'baseline_v2.json'
    summary = json.loads(summ_path.read_text()) if summ_path.exists() and want else {}
    summary['_S_genes'] = int(in_S.sum())
    for ds, source in CELL_SOURCES.items():
        if want and ds not in want: continue
        path = SOURCES[ds][2]
        with h5py.File(path, 'r') as f:
            v = read_elem(f['var']); ctrl = read_elem(f['obs'])['is_control'].to_numpy(bool)
        src_ids = v.index.astype(str).to_numpy() if 'gene_name' in v else v['gene_id'].astype(str).to_numpy()
        mean, rate, info = cell_baseline(path, np.array([pid.get(g, -1) for g in src_ids]), ctrl, in_S)
        save(source, features(genes, mean, rate, in_S), dict(info, method=f'{path.name}: control cells'), summary)
    for line in ('HCT116', 'HEK293T'):
        if want and line not in want: continue
        source = TABLES[line][0]; path = V / f'lib/qc_cells/{source}_qc.h5ad'
        with h5py.File(path, 'r') as f:
            v = read_elem(f['var']); ctrl = read_elem(f['obs'])['is_control'].to_numpy(bool)
        to_panel = np.array([sym.get(s, -1) for s in v.index.astype(str)]); to_panel[~v['measured'].to_numpy(bool)] = -1
        mean, rate, info = cell_baseline(path, to_panel, ctrl, in_S)
        save(source, features(genes, mean, rate, in_S), dict(info, method=f'{path.name}: subsampled control cells'), summary)
    for name, (source, _, cond) in TABLES.items():
        if not name.startswith('CD4') or (want and 'CD4' not in want): continue
        with h5py.File(V / f'lib/tcell/TCELL_{cond}_effects.h5ad', 'r') as f:
            v = read_elem(f['var'])
        cols = np.array([sym.get(s, -1) for s in v.index.astype(str)]); meas = v['measured'].to_numpy(bool) & (cols >= 0)
        cpm = np.full(G, np.nan); cpm[cols[meas]] = v['ctrl_cpm'].to_numpy(float)[meas]
        lib = np.nansum(cpm[in_S]); mean = cpm / lib * 1e6
        save(source, features(genes, mean, np.full(G, np.nan), in_S, missing_detection=True),
             dict(control_cells=None, library_S_genes=int((in_S & np.isfinite(cpm)).sum()),
                  method='pseudobulk control CPM rescaled to sum 1e6 over S (ratio of sums); no per-cell data'), summary)
    for c in 'ABC':
        if want and c not in want: continue
        path = V / f'data/context_{c}.h5ad'
        with h5py.File(path, 'r') as f:
            v = read_elem(f['var']); ctrl = read_elem(f['obs'])['target_gene'].astype(str).to_numpy() == 'non-targeting'
        assert list(v.index.astype(str)) == list(genes)
        mean, rate, info = cell_baseline(path, np.arange(G), ctrl, in_S)
        save(f'challenge_{c}', features(genes, mean, rate, in_S), dict(info, method=f'{path.name}: control cells'), summary)
    summ_path.write_text(json.dumps(summary, indent=1) + '\n')
    print(f'wrote {OUT}/*_baseline_v2.csv and {summ_path.name}')


if __name__ == '__main__':
    main()

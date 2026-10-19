"""Control baselines with a panel-only library size, for every line (lib/effects_v4/*_baseline_panel.csv). Never overwrites
the current *_baseline.csv tables.

Per cell, CPM = count / (sum of the cell's counts over the source genes that map onto the 18,533-gene panel) * 1e6; the
baseline is the mean of that over control cells, plus the share of controls expressing the gene. This is the convention
the challenge lines A-C force (their files hold panel genes only), applied to every line. A source that releases fewer
genes (Replogle K562: 8,248) still has a smaller panel library than A-C; json records the genes in each library.
- K562, RPE1, HepG2, Jurkat, H1: control cells of lib/qc_cells/<source>_qc.h5ad, genes mapped by Ensembl ID as in
  build_effects_v4.py.
- HCT116, HEK293T: control cells of lib/qc_cells/orion_<line>_qc.h5ad (the 20,000 subsampled controls; the current
  baseline used the pseudobulk table's full-library sums instead), panel symbols.
- CD4 T (Rest, Stim8hr, Stim48hr): no per-cell data; the pseudobulk control CPM over measured panel genes is rescaled to
  sum to 1e6 (a ratio of sums, not a mean of per-cell ratios); detection rate NaN.
- Challenge lines A-C: data/context_<X>.h5ad, written as challenge_<X>_baseline_panel.csv.
Also writes lib/effects_v4/baseline_panel.json (cells, library genes, agreement with the current baseline).
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from build_effects_v4 import CELL_SOURCES, OUT, TABLES, V, panel             # noqa: E402
from qc_replogle_umi import SOURCES                                          # noqa: E402

CHUNK = 20_000


def control_baseline(path, to_panel, ctrl, n_panel):
    """Mean per-cell panel CPM and detection over the control rows of a CSR h5ad; to_panel: source column -> panel column."""
    ok = to_panel >= 0; s1 = np.zeros(n_panel); det = np.zeros(n_panel); n = 0; libs = []
    with h5py.File(path, 'r') as f:
        X = f['X']; ip = X['indptr']; n_obs, n_vars = (int(x) for x in X.attrs['shape'])
        keep = sp.diags(ok.astype(float)).tocsr()
        for s in range(0, n_obs, CHUNK):
            e = min(s + CHUNK, n_obs); r = np.flatnonzero(ctrl[s:e])
            if not len(r): continue
            p = ip[s:e + 1]; a, b = int(p[0]), int(p[-1])
            M = sp.csr_matrix((X['data'][a:b], X['indices'][a:b], p - a), shape=(e - s, n_vars))[r].astype(np.float64) @ keep
            lib = np.asarray(M.sum(1)).ravel(); libs.append(lib)
            C = sp.diags(1e6 / np.maximum(lib, 1)) @ M
            np.add.at(s1, to_panel[ok], np.asarray(C.sum(0)).ravel()[ok])
            np.add.at(det, to_panel[ok], np.asarray((M > 0).sum(0)).ravel()[ok]); n += len(r)
    base = np.full(n_panel, np.nan); rate = np.full(n_panel, np.nan)
    cols = np.unique(to_panel[ok]); base[cols] = s1[cols] / n; rate[cols] = det[cols] / n
    libs = np.concatenate(libs)
    return base, rate, dict(control_cells=n, library_genes=int(ok.sum()), median_panel_library=float(np.median(libs)),
                            zero_library_cells=int((libs == 0).sum()))


def write(name, genes, base, rate, info, summary, current=None):
    df = pd.DataFrame(dict(control_mean_cpm=base, detection_rate=rate), index=pd.Index(genes, name='gene'))
    path = OUT / f'{name}_baseline_panel.csv'; assert not path.name.endswith('_baseline.csv'); df.to_csv(path)
    if current is not None:
        old = pd.read_csv(OUT / f'{current}_baseline.csv', index_col=0)['control_mean_cpm'].reindex(genes).to_numpy(float)
        both = np.isfinite(old) & np.isfinite(base) & ((old > 0) | (base > 0))
        info['vs_current'] = dict(genes=int(both.sum()), spearman=float(spearmanr(old[both], base[both])[0]),
                                  median_ratio_new_over_old=float(np.median(base[both & (old > 0)] / old[both & (old > 0)])),
                                  measured_same=bool((np.isfinite(old) == np.isfinite(base)).all()))
    info['measured_panel_genes'] = int(np.isfinite(base).sum())
    summary[name] = info
    print(f"{name:22s} {json.dumps({k: v for k, v in info.items() if k != 'method'})}", flush=True)


def main():
    var, _ = panel(); genes = var.index.to_numpy(); G = len(genes)
    pid = {g: i for i, g in enumerate(var.gene_id) if g}; sym = {s: i for i, s in enumerate(genes)}
    summary = {}
    for ds, source in CELL_SOURCES.items():
        path = SOURCES[ds][2]
        with h5py.File(path, 'r') as f:
            v = read_elem(f['var']); ctrl = read_elem(f['obs'])['is_control'].to_numpy(bool)
        src_ids = v.index.astype(str).to_numpy() if 'gene_name' in v else v['gene_id'].astype(str).to_numpy()
        to_panel = np.array([pid.get(g, -1) for g in src_ids])
        base, rate, info = control_baseline(path, to_panel, ctrl, G)
        info['method'] = f'{path.name}: control cells, panel genes by Ensembl ID'
        write(source, genes, base, rate, info, summary, current=source)
    for line in ('HCT116', 'HEK293T'):
        source = TABLES[line][0]; path = V / f'lib/qc_cells/{source}_qc.h5ad'
        with h5py.File(path, 'r') as f:
            v = read_elem(f['var']); ctrl = read_elem(f['obs'])['is_control'].to_numpy(bool)
        to_panel = np.array([sym.get(s, -1) for s in v.index.astype(str)]); to_panel[~v['measured'].to_numpy(bool)] = -1
        base, rate, info = control_baseline(path, to_panel, ctrl, G)
        info['method'] = f'{path.name}: subsampled control cells, measured panel genes'
        write(source, genes, base, rate, info, summary, current=source)
    for name, (source, _, cond) in TABLES.items():
        if not name.startswith('CD4'): continue
        with h5py.File(V / f'lib/tcell/TCELL_{cond}_effects.h5ad', 'r') as f:
            v = read_elem(f['var'])
        cols = np.array([sym.get(s, -1) for s in v.index.astype(str)]); meas = v['measured'].to_numpy(bool) & (cols >= 0)
        cpm = v['ctrl_cpm'].to_numpy(float); base = np.full(G, np.nan)
        base[cols[meas]] = cpm[meas] / cpm[meas].sum() * 1e6
        info = dict(control_cells=None, library_genes=int(meas.sum()), method='pseudobulk control CPM rescaled to the '
                    'measured panel genes (ratio of sums); no per-cell data')
        write(source, genes, base, np.full(G, np.nan), info, summary, current=source)
    for c in 'ABC':
        path = V / f'data/context_{c}.h5ad'
        with h5py.File(path, 'r') as f:
            v = read_elem(f['var']); ctrl = read_elem(f['obs'])['target_gene'].astype(str).to_numpy() == 'non-targeting'
        assert list(v.index.astype(str)) == list(genes)
        base, rate, info = control_baseline(path, np.arange(G), ctrl, G)
        info['method'] = f'{path.name}: control cells, all panel genes'
        write(f'challenge_{c}', genes, base, rate, info, summary)
    (OUT / 'baseline_panel.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(f'wrote {OUT}/*_baseline_panel.csv and baseline_panel.json')


if __name__ == '__main__':
    main()

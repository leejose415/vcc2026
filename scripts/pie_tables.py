"""PIE predictions -> tables in the effects_v4 layout: lib/pie_tables/pie_<line>.h5ad (layers lfc_pred, p_de; panel
columns; NaN where PIE gave no prediction), plus lib/pie_tables/pie_<line>.npz (targets, genes, lfc, p_de) for
scripts/k562_benchmark.py --pie. Rows: the line's effects_v4 knockdowns (K562, RPE1, HepG2, Jurkat) or our 300 targets
(A-C), in that order; a knockdown PIE did not return stays NaN. Never made for HCT116, HEK293T or H1 (PIE saw them).
Usage: python scripts/pie_tables.py --lines K562 A ...
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from anndata.io import read_elem

V = Path.home() / 'vcc'
OUT = V / 'lib/pie_tables'
SOURCE = {'K562': 'replogle_k562_gw', 'RPE1': 'replogle_rpe1', 'HepG2': 'nadig_hepg2', 'Jurkat': 'nadig_jurkat'}
LEAKED = {'HCT116', 'HEK293T', 'H1'}


def build(line):
    if line in LEAKED: raise SystemExit(f'{line}: PIE is leaked; no PIE table')
    with h5py.File(V / 'lib/effects_v4/replogle_k562_gw.h5ad', 'r') as f:
        pvar = read_elem(f['var'])
    genes = pvar.index.astype(str).to_numpy(); gpos = {g: i for i, g in enumerate(genes)}
    t = pq.read_table(V / f'lib/pie/predictions/{line}.parquet')
    meta = json.loads(t.schema.metadata[b'pie']); pie_genes = next(iter(meta['genes'].values()))
    cols = np.array([gpos.get(g, -1) for g in pie_genes]); ok = cols >= 0
    df = t.to_pandas()
    if line in SOURCE:
        with h5py.File(V / f'lib/effects_v4/{SOURCE[line]}.h5ad', 'r') as f:
            kobs = read_elem(f['obs'])
        obs = pd.DataFrame(dict(target_gene=kobs['target_gene'].astype(str).to_numpy(),
                                target_gene_id=kobs['target_gene_id'].astype(str).to_numpy()), index=kobs.index.astype(str))
    else:
        tg = sorted(pd.read_csv(V / 'data/pert_counts.csv')['target_gene'].astype(str))
        obs = pd.DataFrame(dict(target_gene=tg, target_gene_id=[pvar['gene_id'].iloc[gpos[g]] if g in gpos else '' for g in tg]),
                           index=[f'pie_{line}|{g}' for g in tg])
    row_of = {p: i for i, p in enumerate(df['perturbation'].astype(str))}
    L = np.full((len(obs), len(genes)), np.nan, np.float32); P = np.full_like(L, np.nan)
    hit = 0
    for i, g in enumerate(obs['target_gene']):
        j = row_of.get(g)
        if j is None: continue
        L[i, cols[ok]] = np.asarray(df['lfc_pred'].iloc[j], np.float32)[ok]
        P[i, cols[ok]] = np.asarray(df['p_de'].iloc[j], np.float32)[ok]; hit += 1
    obs['source'] = 'pie'; obs['background'] = line; obs['predicted'] = np.isfinite(L).any(1)
    a = ad.AnnData(obs=obs, var=pvar.copy(), layers=dict(lfc_pred=L, p_de=P))
    a.uns['pie'] = dict(model='arcinstitute/PIE_replogle_xdataset (best_auprc)', parquet=f'lib/pie/predictions/{line}.parquet',
                        genes_predicted=int(ok.sum()), rows_predicted=int(hit),
                        note='control-only context; log1p CP10k; lfc_pred is PIE log2 fold change; p_de its DE probability')
    OUT.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(OUT / f'pie_{line}.h5ad', compression='gzip')
    np.savez(OUT / f'pie_{line}.npz', targets=obs['target_gene'].to_numpy().astype(str), genes=genes, lfc=L, p_de=P)
    info = dict(line=line, rows=len(obs), rows_predicted=hit, genes_predicted=int(ok.sum()), pie_rows=len(df))
    print(json.dumps(info), flush=True)
    return info


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--lines', nargs='+', required=True); a = ap.parse_args()
    for line in a.lines: build(line)


if __name__ == '__main__':
    main()

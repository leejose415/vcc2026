"""Data inventory for vcc_2026_4, built from the files: every line in lib/effects_v4/ plus the challenge lines A-C.

Per effects_v4 table: source, background, condition, role in the transformer (target: trained and validated on;
validation: H1, a target only in final mode; evidence-only: CD4 T), knockdowns (rows), panel genes measured, how many
of our 300 targets it has as knockdowns, and PIE status (clean: not in PIE's training/validation; leaked: seen by PIE;
n/a: no PIE table is made). Challenge lines: their three closest lines by basal_sim (results/analysis/basal_sim.csv).
Writes results/vcc_2026_4/data_inventory.csv and challenge_neighbours.csv; prints both.
"""
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

V = Path.home() / 'vcc'
E = V / 'lib/effects_v4'
OUT = V / 'results/vcc_2026_4'
ROLE = {'K562': 'target', 'RPE1': 'target', 'HepG2': 'target', 'Jurkat': 'target', 'HCT116': 'target',
        'HEK293T': 'target', 'H1': 'validation (target in final mode)', 'CD4_T': 'evidence-only'}
PIE = {'K562': 'clean', 'RPE1': 'clean', 'HepG2': 'clean', 'Jurkat': 'clean', 'HCT116': 'leaked', 'HEK293T': 'leaked',
       'H1': 'leaked', 'CD4_T': 'n/a (no PIE table)'}
BASAL_NAME = {'replogle_k562_gw': 'K562', 'replogle_rpe1': 'RPE1', 'nadig_hepg2': 'HepG2', 'nadig_jurkat': 'Jurkat',
              'orion_hct116': 'HCT116', 'orion_hek293t': 'HEK293T', 'marson_cd4_stim48hr': 'CD4T_Stim48hr', 'arc_h1': 'H1'}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    targets = set(pd.read_csv(V / 'data/pert_counts.csv')['target_gene'].astype(str))
    rows = []
    for p in sorted(E.glob('*.h5ad')):
        with h5py.File(p, 'r') as f:
            obs = read_elem(f['obs']); lfc = f['layers/lfc']
            measured = np.zeros(lfc.shape[1], bool)
            for s in range(0, lfc.shape[0], 2000):
                measured |= np.isfinite(lfc[s:s + 2000]).any(0)
        bg = str(obs['background'].iloc[0])
        rows.append(dict(source=p.stem, background=bg, condition=str(obs['condition'].iloc[0]), role=ROLE.get(bg, '?'),
                         knockdowns=len(obs), panel_genes_measured=int(measured.sum()),
                         our_300_targets=int(obs['target_gene'].astype(str).isin(targets).sum()),
                         median_n_de_bh=float(obs['n_de_bh'].median()), lfc_type=str(obs['lfc_type'].iloc[0]),
                         se_type=str(obs['se_type'].iloc[0]), pie=PIE.get(bg, '?')))
    inv = pd.DataFrame(rows)
    inv.to_csv(OUT / 'data_inventory.csv', index=False)
    basal = pd.read_csv(V / 'results/analysis/basal_sim.csv', index_col=0)
    train = [c for c in basal.columns if not c.startswith('line_')]
    near = []
    for c in 'ABC':
        s = basal.loc[f'line_{c}', train].sort_values(ascending=False)
        near.append(dict(line=c, **{f'closest_{i + 1}': f'{s.index[i]} ({s.iloc[i]:.2f})' for i in range(3)}))
    near = pd.DataFrame(near); near.to_csv(OUT / 'challenge_neighbours.csv', index=False)
    pd.set_option('display.width', 250)
    print(inv.to_string(index=False)); print(); print(near.to_string(index=False))


if __name__ == '__main__':
    main()

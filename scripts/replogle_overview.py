"""Quick look at the Replogle K562 genome-wide Perturb-seq file without loading the 61 GiB matrix.

Opens the .h5ad in backed mode (metadata in RAM, expression matrix stays on disk), prints the
structure, cells per knockdown gene, guide and QC summaries, and reads a small random sample of
rows to show what the expression values look like.

Usage:  ~/venv/bin/python scripts/replogle_overview.py [path.h5ad]
Writes: lib/replogle/cells_per_gene.csv
"""
import resource
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

PATH = Path(sys.argv[1] if len(sys.argv) > 1 else
            Path.home() / 'vcc/lib/replogle/K562_gwps_raw_singlecell_01.h5ad')
CONTROL = 'non-targeting'
SAMPLE_ROWS = 2000  # ~66 MB dense; enough to judge value range and sparsity

pd.set_option('display.width', 160)
pd.set_option('display.max_columns', 20)


def section(title):
    print(f'\n{"=" * 8} {title} {"=" * (70 - len(title))}')


a = ad.read_h5ad(PATH, backed='r')
obs, var = a.obs, a.var

section('object')
print(a)

section('obs (one row per cell), first 5')
print(obs.head())

section('var (one row per measured gene), first 5')
print(var.head())

section('knockdowns')
counts = obs['gene'].value_counts()
n_ctrl = int(counts.get(CONTROL, 0))
per_gene = counts.drop(CONTROL, errors='ignore')
per_gene = per_gene[per_gene > 0]
print(f'cells total            : {a.n_obs:,}')
print(f'control cells ({CONTROL}): {n_ctrl:,} ({n_ctrl / a.n_obs:.1%})')
print(f'knockdown genes        : {len(per_gene):,}')
print(f'measured genes (var)   : {a.n_vars:,}')

section('cells per knockdown gene')
print(per_gene.describe(percentiles=[.05, .25, .5, .75, .95]).round(1).to_string())
bins = [0, 10, 25, 50, 100, 200, 500, np.inf]
labels = ['<10', '10-24', '25-49', '50-99', '100-199', '200-499', '500+']
print('\nhistogram (genes per bin):')
print(pd.cut(per_gene, bins=bins, labels=labels, right=False).value_counts().sort_index().to_string())
print('\nmost cells:')
print(per_gene.head(10).to_string())
print('\nfewest cells:')
print(per_gene.tail(10).to_string())

section('guides')
guides_per_gene = obs.groupby('gene', observed=True)['sgID_AB'].nunique()
print('distinct sgRNA pairs per gene:')
print(guides_per_gene.drop(CONTROL, errors='ignore').value_counts().sort_index().to_string())
print(f'control guides: {guides_per_gene.get(CONTROL, 0)}')
tx = obs.groupby('gene', observed=True)['transcript'].nunique().drop(CONTROL, errors='ignore')
print(f'genes targeted at >1 transcript (TSS): {(tx > 1).sum()}')

section('knockdown targets measured in var?')
measured = set(var['gene_name'].astype(str))
in_var = per_gene.index.astype(str).isin(measured)
print(f'{in_var.sum():,} of {len(per_gene):,} knockdown genes are also measured, '
      'so their own knockdown efficiency can be checked directly')

section('batches (gem_group)')
gg = obs['gem_group'].value_counts()
print(f'{len(gg)} GEM groups, cells per group: min {gg.min():,}, median {int(gg.median()):,}, max {gg.max():,}')

section('per-cell QC')
print(obs[['UMI_count', 'mitopercent', 'z_gemgroup_UMI', 'core_scale_factor']]
      .describe(percentiles=[.01, .05, .5, .95, .99]).round(3).to_string())

section(f'expression values ({SAMPLE_ROWS:,} random cells)')
rng = np.random.default_rng(0)
rows = np.sort(rng.choice(a.n_obs, SAMPLE_ROWS, replace=False))
X = np.asarray(a.X[rows])
print(f'dtype {X.dtype}, min {X.min():g}, max {X.max():g}, mean {X.mean():.3f}')
print(f'fraction zeros: {(X == 0).mean():.1%}  (dense on disk; '
      f'as sparse this file would be roughly {(X != 0).mean() * 61.1 * 2:.0f} GiB)')
print(f'all integers (raw counts?): {bool(np.all(X == np.round(X)))}')
print(f'per-cell sum in sample: median {np.median(X.sum(1)):,.0f} '
      f'vs obs UMI_count median {obs["UMI_count"].iloc[rows].median():,.0f}')
print('\nfirst 5 sampled cells x first 8 genes:')
print(pd.DataFrame(X[:5, :8], index=obs.index[rows[:5]], columns=var['gene_name'].iloc[:8]))

out = PATH.parent / 'cells_per_gene.csv'
per_gene.rename('n_cells').rename_axis('gene').to_csv(out)
print(f'\nwrote {out}')
print(f'peak RAM: {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.1f} GB')

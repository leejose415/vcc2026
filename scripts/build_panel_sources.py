"""Panel tables: every source on one shape, ready for modelling.

For each cleaned source (lib/filters/<name>.h5ad, built by knockdown_filters.py):

  rows     the panel targets in a fixed order: our 300 (data/pert_counts.csv order), then H1's other targets
  columns  the 18,533 genes of data/gene_names.csv, in that order
  layers   lfc (log2 fold change vs that source's controls) and se (its standard error); X is empty

  NaN vs 0:
    whole row NaN    the knockdown failed filters 1-3 there, or the source never screened it
    column NaN       gene not measurable in that source: not measured, or < 1 CPM in its controls (filter 4)
    0                a real "no change"; values are copied, never filled

  obs   panel (ours / H1 / both), status (kept / which filters failed / not screened), n_cells, ontarget_lfc, padj
  var   ctrl_cpm (and ctrl_umi_per_cell where the source has UMI counts); NaN = not measured

Coverage table: one row per target, each source's status, n_sources, n_backgrounds (the three T-cell tables are
one background), and decision = "no source: needs a decision" where nothing survived. Those targets have no
measured effect anywhere; the only thing we know is that the target itself should drop.

After writing, every table is re-read and checked against the rules above.

Usage:  ~/venv/bin/python scripts/build_panel_sources.py [--min-cpm 1]
Writes: lib/sources/<NAME>/effects_panel.h5ad, lib/sources/coverage.csv
"""
import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import knockdown_filters as kf                                       # noqa: E402  (filters + cached numbers)

V = Path.home() / 'vcc'
LIB = V / 'lib'
GENES = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy()
BACKGROUND = dict(HCT116='HCT116', HEK293T='HEK293T', K562_gwps='K562', H1='H1',
                  CD4_Stim48hr='CD4_T', CD4_Stim8hr='CD4_T', CD4_Rest='CD4_T')


def panel_targets():
    ours = pd.read_csv(V / 'data/pert_counts.csv')['target_gene'].astype(str).tolist()
    h1 = ad.read_h5ad(LIB / 'h1/h1_pseudobulk.h5ad', backed='r')
    h1_targets = [t for t in h1.obs_names.astype(str) if t != kf.CONTROL]; h1.file.close()
    extra = sorted(set(h1_targets) - set(ours))
    targets = ours + extra
    label = ['both' if t in set(h1_targets) else 'ours' for t in ours] + ['H1'] * len(extra)
    return np.array(targets), np.array(label)


def filter_status(name, settings):
    """Every screened knockdown's QC numbers and keep/drop reason, from the filter step's cache."""
    s = pd.read_csv(LIB / 'filters/_cache' / f'{name}.csv', keep_default_na=False, na_values=[''])
    s['challenge_target'] = False
    a = argparse.Namespace(min_cells=settings['min_cells'], min_umi=settings['min_umi_per_control_cell'],
                           min_cpm=settings['min_ctrl_cpm'], expr_rule=settings['expr_rule'],
                           max_lfc=settings['max_ontarget_lfc'], alpha=settings['max_padj'])
    s, *_ = kf.apply_filters(s, a)
    s['status'] = np.where(s['status'] == 'keep', 'kept', 'failed: ' + s['failed'])
    return s.drop_duplicates('target').set_index('target')


def build(name, targets, label, a):
    clean = ad.read_h5ad(LIB / 'filters' / f'{name}.h5ad')
    settings = json.loads(clean.uns['filters'])
    qc = filter_status(name, settings)

    gpos = {g: i for i, g in enumerate(clean.var_names.astype(str))}
    cols = np.array([gpos.get(g, -1) for g in GENES])               # -1: not measured
    ctrl = np.full(len(GENES), np.nan); ok = cols >= 0
    ctrl[ok] = clean.var['ctrl_cpm'].to_numpy(float)[cols[ok]]
    measurable = ok & (ctrl >= a.min_cpm)                             # filter 4

    rpos = {t: i for i, t in enumerate(clean.obs_names.astype(str))}
    rows = np.array([rpos.get(t, -1) for t in targets])               # -1: failed or never screened
    lfc = np.full((len(targets), len(GENES)), np.nan, np.float32); se = lfc.copy()
    kept = rows >= 0
    src_rows, dst_cols = rows[kept], np.flatnonzero(measurable)
    lfc[np.ix_(kept, dst_cols)] = np.asarray(clean.X)[np.ix_(src_rows, cols[measurable])]
    se[np.ix_(kept, dst_cols)] = np.asarray(clean.layers['se'])[np.ix_(src_rows, cols[measurable])]

    status = np.array([qc.loc[t, 'status'] if t in qc.index else 'not screened' for t in targets], dtype=object)
    status[kept] = 'kept'                                             # the cleaned file is the authority on keep
    obs = pd.DataFrame(dict(panel=label, status=status), index=pd.Index(targets, name='target'))
    for c in ('n_cells', 'ontarget_lfc', 'padj'):
        obs[c] = [qc.loc[t, c] if t in qc.index else np.nan for t in targets]
    var = pd.DataFrame(dict(ctrl_cpm=ctrl), index=pd.Index(GENES, name='gene'))
    if 'ctrl_umi_per_cell' in clean.var:
        umi = np.full(len(GENES), np.nan); umi[ok] = clean.var['ctrl_umi_per_cell'].to_numpy(float)[cols[ok]]
        var['ctrl_umi_per_cell'] = umi

    out = ad.AnnData(X=None, obs=obs, var=var, layers=dict(lfc=lfc, se=se), uns=dict(
        source=name, background=BACKGROUND[name], built_from=str(LIB / 'filters' / f'{name}.h5ad'),
        lfc=clean.uns['X'], se=clean.uns['se'], knockdown_filters=clean.uns['filters'],
        gene_filter=f'NaN where the gene is not measured or ctrl_cpm < {a.min_cpm:g}',
        row_filter='NaN rows: knockdown failed filters 1-3 (obs.status says which) or was not screened'))
    folder = a.out / name; folder.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(folder / 'effects_panel.h5ad', compression='gzip')
    return obs['status'], dict(measured=int(ok.sum()), measurable=int(measurable.sum()),
                               dropped_below_cpm=int((ok & ~measurable).sum()))


def check(name, targets, a):
    """Re-read the written table and test every rule. Returns one row of results."""
    t = ad.read_h5ad(a.out / name / 'effects_panel.h5ad')
    lfc, se = t.layers['lfc'], t.layers['se']
    kept = (t.obs['status'] == 'kept').to_numpy()
    measurable = (t.var['ctrl_cpm'] >= a.min_cpm).to_numpy()
    finite = np.isfinite(lfc)
    r = dict(
        shape_ok=t.shape == (len(targets), len(GENES)),
        gene_order_ok=bool((t.var_names.to_numpy() == GENES).all()),
        target_order_ok=bool((t.obs_names.to_numpy() == targets).all()),
        layers_ok=sorted(t.layers.keys()) == ['lfc', 'se'],
        failed_rows_all_nan=bool(~finite[~kept].any()),
        unmeasurable_cols_all_nan=bool(~finite[:, ~measurable].any()),
        kept_x_measurable_all_finite=bool(finite[np.ix_(kept, measurable)].all()),
        se_matches_lfc=bool((np.isfinite(se) == finite).all()),
        se_positive=bool((se[finite] > 0).all()),
        rows_kept=int(kept.sum()), genes_measurable=int(measurable.sum()),
        exact_zero_lfc=int((lfc == 0).sum()),
    )
    clean = ad.read_h5ad(LIB / 'filters' / f'{name}.h5ad', backed='r')  # values copied unchanged?
    pick = [x for x in t.obs_names[kept][:20]]
    gi = np.flatnonzero(measurable)[:500]
    cpos = {g: i for i, g in enumerate(clean.var_names)}
    src = np.asarray(clean[pick].X)[:, [cpos[g] for g in GENES[gi]]]
    r['values_match_clean'] = bool(np.array_equal(src, t[pick].layers['lfc'][:, gi]))
    clean.file.close()
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sources', nargs='+', default=list(BACKGROUND), choices=list(BACKGROUND))
    ap.add_argument('--min-cpm', type=float, default=1.0, help='filter 4: genes below this control CPM -> NaN')
    ap.add_argument('--out', type=Path, default=LIB / 'sources')
    a = ap.parse_args()

    targets, label = panel_targets()
    print(f'panel: {len(targets)} targets ({(label != "H1").sum()} ours, {(label == "H1").sum()} H1-only) '
          f'x {len(GENES):,} genes')
    cov = pd.DataFrame(dict(panel=label), index=pd.Index(targets, name='target'))
    genes, checks = {}, {}
    for name in a.sources:
        if not (LIB / 'filters' / f'{name}.h5ad').exists():
            print(f'{name}: no cleaned file, skipped'); continue
        cov[name], genes[name] = build(name, targets, label, a)
        checks[name] = check(name, targets, a)
        print(f'{name}: {int((cov[name] == "kept").sum())} targets kept, '
              f'{genes[name]["measurable"]:,} measurable genes')

    srcs = list(genes)
    kept = cov[srcs] == 'kept'
    cov['n_sources'] = kept.sum(axis=1)
    cov['n_backgrounds'] = kept.T.groupby([BACKGROUND[s] for s in srcs]).any().sum(axis=0)
    cov['decision'] = np.where(cov['n_sources'] == 0, 'no source: needs a decision', '')
    cov.to_csv(a.out / 'coverage.csv')

    pd.set_option('display.width', 240); pd.set_option('display.max_columns', 30)
    print('\ngenes per source (of 18,533):')
    print(pd.DataFrame(genes).T.to_string())
    print('\nchecks:')
    print(pd.DataFrame(checks).T.to_string())
    print('\ntargets by number of backgrounds with a surviving knockdown:')
    print(pd.crosstab(cov['n_backgrounds'], cov['panel']).to_string())
    none = cov.index[cov['n_sources'] == 0]
    print(f'\n{len(none)} targets with no source (need a decision), '
          f'{int((cov.loc[none, "panel"] != "H1").sum())} of them ours: {", ".join(none[cov.loc[none, "panel"] != "H1"])}')
    print(f'\nwrote {a.out}/<NAME>/effects_panel.h5ad and {a.out / "coverage.csv"}')


if __name__ == '__main__':
    main()

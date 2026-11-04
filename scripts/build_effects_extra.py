"""Effect tables for the PIE-relevant extra sources, with the v5 builder's own math and checks, into lib/effects_v5_extra/
(never lib/effects_v5/: these are not scored-task tables). build_effects_v5.py is imported, not modified: its output
folder and, for Tahoe, its kept-cell grouping are redirected for the duration of a build.

  jiang   Jiang 2025 pathway Perturb-seq, one source per cell line x pathway (Jiang_<PATHWAY>_<LINE>): build_cells on the
          QC'd file with the standard knockdown QC (guide_qc keep -> knockdown_qc keep: on-target LFC <= -0.5 and p < 0.05,
          unmeasurable targets kept unverified), unchanged; rows tagged modality='CRISPRko', background = the line,
          condition = the pathway. Parse split-pool depth (~3.9k UMI): SEs and Wilcoxon calls are noisier (build JSON).
  tahoe   Tahoe-100M, one source per cell line (Tahoe_<line>): rows = single-target drugs at their highest dose (>= 30
          QC'd cells), each a perturbation of its annotated target gene, against the line's DMSO_TF (the controls).
          LFC epsilon 1e-9, panel-only library, <= 400 cells / <= 18,400 controls, Wilcoxon + BH as v5. Tagged
          modality='drug', significance_unknown=True (~1.7k UMI per cell, 10x below the scorer; drug effects are not
          on-target-verified). <source>_baseline.csv = the line's DMSO_TF basal profile (context background).
Output lib/effects_v5_extra/<kind>/<source>.h5ad (+ _baseline.csv, _control_quantiles.npz, _build.json), schema as v5.

Usage: ~/venv/bin/python scripts/build_effects_extra.py --kind jiang|tahoe [--sources KEY ...]
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
import qc_replogle_umi as Q                                                              # noqa: E402

V = Path.home() / 'vcc'
OUT = V / 'lib/effects_v5_extra'
MIN_CELLS = 30


def tahoe_kept(dataset, obs):
    """Tahoe grouping for build_cells: one row per drug (guide_id = drug_<dose>uM) with >= 30 QC'd cells; target_gene_id
    is a unique 'ENSG...|drug' label (so validation keeps drugs apart); the target is resolved back afterwards."""
    ctrl = obs['is_control'].to_numpy(bool); guide = obs['guide_id'].astype(str).to_numpy()
    tid = obs['target_gene_id'].astype(str).to_numpy(); tg = obs['target_gene'].astype(str).to_numpy()
    d = pd.DataFrame(dict(guide=guide, tid=tid, tg=tg))[~ctrl]
    n = d.groupby('guide').size(); first = d.groupby('guide')[['tid', 'tg']].first()
    rows = n[n >= MIN_CELLS].index
    kept = pd.DataFrame(dict(target_gene=first.loc[rows, 'tg'].to_numpy(), target_gene_id=[f'{first.loc[g, "tid"]}|{g}' for g in rows],
                             drug_dose=list(rows), gene_id=first.loc[rows, 'tid'].to_numpy(), n_cells_kept=n[rows].to_numpy().astype(int),
                             measurable=True))
    row_of = {g: i for i, g in enumerate(rows)}
    kd = np.array([-1 if c else row_of.get(g, -1) for c, g in zip(ctrl, guide)])
    skipped = n[n < MIN_CELLS]
    tahoe_kept.skipped = dict(zip(skipped.index, skipped.to_numpy().astype(int).tolist()))
    return kept.reset_index(drop=True), kd, dict(knockdown_qc='none (drugs: grouping by drug at its highest dose, >= 30 QC\'d cells)',
                                                 guide_qc='none (drugs have no on-target rule)')


def retag(path, kind, key, kept):
    """Rewrite obs (modality, background, condition, row ids) and, for drugs, resolve the target: efficiency = its LFC,
    removed from de_sig / n_real (build_cells could not, the row label being 'ENSG...|drug')."""
    a = ad.read_h5ad(path)
    _, _, bg = Q.NEW_CONTEXTS[key]; src = a.obs['source'].astype(str).iloc[0]
    o = a.obs.copy()
    if kind == 'jiang':
        cond = key.split('_')[1]
        o['modality'] = 'CRISPRko'; o['significance_unknown'] = False
        o['lfc_type'] = B.LFC_TYPE; o['se_type'] = B.SE_TYPE
    else:
        cond_ = kept['drug_dose'].to_numpy(); gid = kept['gene_id'].to_numpy()
        col = {g: j for j, g in enumerate(a.var['gene_id'].astype(str)) if g}
        tcol = np.array([col.get(g, -1) for g in gid])
        sig = a.layers['de_sig'].tocsr().tolil()
        for i, j in enumerate(tcol):
            if j >= 0 and sig[i, j]: sig[i, j] = False
        sig = sig.tocsr(); sig.eliminate_zeros(); a.layers['de_sig'] = sig
        o['n_real'] = np.asarray(sig.sum(1)).ravel().astype(float)
        lfc = a.layers['lfc']
        o['efficiency'] = [float(lfc[i, j]) if j >= 0 else np.nan for i, j in enumerate(tcol)]
        o['efficiency_known'] = tcol >= 0
        o['target_gene_id'] = gid; o['drug_dose'] = cond_
        o['modality'] = 'drug'; o['significance_unknown'] = True
        cond = None
    o['background'] = bg
    o['condition'] = cond if cond is not None else o['drug_dose']
    o.index = pd.Index([f'{src}|{bg}|{c}|{t}' for c, t in zip(o['condition'].astype(str), o['target_gene_id'].astype(str))], name='row_id')
    if not o.index.is_unique: raise ValueError(f'{key}: duplicate row ids after retagging')
    for c in ('source', 'background', 'condition', 'modality', 'lfc_type', 'se_type'): o[c] = o[c].astype(str).astype('category')
    a.obs = o
    meta = dict(a.uns['effects_v5'])
    meta.update(modality='CRISPRko' if kind == 'jiang' else 'drug', significance_unknown=kind == 'tahoe',
                table_family='lib/effects_v5_extra: PIE-relevant extra source, NOT a scored-task table',
                depth_note=('Parse split-pool, shallow (~3.9k UMI per cell): SEs and Wilcoxon calls noisier than the CRISPRi sources'
                            if kind == 'jiang' else '~1.7k UMI per cell (10x below the scorer): labels computed as v5 but flagged '
                            'significance_unknown; drug effects are not on-target verified'))
    if kind == 'tahoe':
        meta.update(perturbation='single-target drug at its highest dose, as a perturbation of its annotated target gene',
                    controls='DMSO_TF (plate-matched vehicle)', basal_profile=f'{src}_baseline.csv (the DMSO_TF profile of {bg})',
                    conditions_skipped_lt_30_cells=json.dumps(tahoe_kept.skipped))
    a.uns['effects_v5'] = meta
    a.write_h5ad(path, compression='gzip')
    return a


def build(kind, key, batch_nnz):
    name = Q.NEW_CONTEXTS[key][0]
    B.OUT = OUT / kind; B.CELL_SOURCES[key] = name
    orig = B.kept_cells
    if kind == 'tahoe': B.kept_cells = tahoe_kept
    try:
        t0 = time.time()
        a, ctx = B.build_cells(key, batch_nnz)
        kept = ctx['kept']
        checks = B.validate(key, a, ctx)
    finally:
        B.kept_cells = orig
    path = B.OUT / f'{name}.h5ad'
    a = retag(path, kind, key, kept)
    _, _, path_cells, _ = Q.SOURCES[key]
    with h5py.File(path_cells, 'r') as f: var_src = read_elem(f['var'])
    to_panel, _ = B.source_genes(key, var_src, B.panel()[0])
    lfc = a.layers['lfc']
    required = dict(fast_lfc_pass=checks['fast']['lfc_pass'], fast_sig_pass=checks['fast']['sig_pass'],
                    fast_quantiles_pass=checks['fast']['quantiles_pass'], gated_sets_identical=checks['gated_sets_identical'],
                    epsilon_is_1e9=a.uns['effects_v5']['epsilon'] == 1e-9, n_cells_used_le_400=bool(a.obs.n_cells_used.max() <= B.N_CELLS))
    failing = [k for k, v in required.items() if not v]
    summary = dict(source=name, kind=kind, modality=a.uns['effects_v5']['modality'], rows=int(a.n_obs),
                   measured_genes=int(np.isfinite(lfc).any(0).sum()), source_genes=int(len(var_src)),
                   source_genes_without_panel_mapping=int((to_panel < 0).sum()),
                   source_genes_without_ensembl=int((var_src['gene_id'].astype(str) == '').sum()),
                   targets_off_panel=int((~a.obs['efficiency_known'].astype(bool)).sum()) if kind == 'tahoe' else None,
                   de_tested_genes=a.uns['effects_v5'].get('de_tested_genes'), median_n_real=float(a.obs.n_real.median()),
                   median_cells_used=float(a.obs.n_cells_used.median()), thinning_factor=ctx['f'],
                   median_control_panel_library=a.uns['effects_v5']['median_control_panel_library_after_thinning'],
                   conditions_skipped_lt_30_cells=(tahoe_kept.skipped if kind == 'tahoe' else 'see results/qc/<source>/knockdown_qc.csv'),
                   depth_note=a.uns['effects_v5']['depth_note'], all_pass=not failing, failing_checks=failing,
                   library_check=ctx['library_check'], count_check=ctx['count_check'], checks=checks,
                   minutes=round((time.time() - t0) / 60, 1))
    (B.OUT / f'{name}_build.json').write_text(json.dumps(summary, indent=1, default=str) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k not in ('checks', 'library_check', 'count_check', 'conditions_skipped_lt_30_cells')},
                     indent=1, default=str), flush=True)
    return not failing


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--kind', required=True, choices=['jiang', 'tahoe']); p.add_argument('--sources', nargs='*')
    p.add_argument('--batch-nnz', type=float, default=1.2e9)
    a = p.parse_args()
    keys = a.sources or [k for k in Q.NEW_CONTEXTS if k.startswith('Jiang_' if a.kind == 'jiang' else 'Tahoe_')]
    bad = []
    for k in keys:
        if not Q.SOURCES[k][2].exists(): print(f'{k}: no QC\'d file, skipped', flush=True); continue
        if not build(a.kind, k, a.batch_nnz): bad.append(k)
    if bad: print('VALIDATION FAILED:', bad, file=sys.stderr); sys.exit(1)

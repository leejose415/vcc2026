"""Checks on a leaderboard file (generator prediction.h5ad + its .vcc).

  python filecheck.py pbulk <prediction.h5ad> <out.npz>          per-knockdown summed counts (line, target order of obs)
  python filecheck.py check <new run dir> <old pbulk.npz> <old .vcc> <pred dir with A/ B/ C/> <array>
"""
import json
import sys
import tarfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

V = Path.home() / 'vcc'


def blocks(path):
    with h5py.File(path, 'r') as f:
        obs = read_elem(f['obs'])
    key = obs['context'].astype(str) + '|' + obs['target_gene'].astype(str)
    change = np.flatnonzero(np.r_[True, key.values[1:] != key.values[:-1], True])
    return obs, [(key.values[a], a, b) for a, b in zip(change[:-1], change[1:])]


def pbulk(path, out):
    obs, bl = blocks(path)
    with h5py.File(path, 'r') as f:
        X = f['X']; G = X.attrs['shape'][1] if 'shape' in X.attrs else read_elem(f['var']).shape[0]
        indptr = X['indptr'][:]; S = np.zeros((len(bl), G)); lib_max = 0.0; mn = np.inf; nonint = 0
        for i, (k, a, b) in enumerate(bl):
            s, e = indptr[a], indptr[b]; d = X['data'][s:e].astype(np.float64); j = X['indices'][s:e]
            S[i] = np.bincount(j, weights=d, minlength=G)
            rows = np.repeat(np.arange(b - a), np.diff(indptr[a:b + 1]))
            lib_max = max(lib_max, np.bincount(rows, weights=d, minlength=b - a).max()); mn = min(mn, d.min(initial=0))
            nonint += int((d != np.round(d)).sum())
    np.savez(out, keys=np.array([k for k, _, _ in bl]), sums=S, max_cell_total=lib_max, min_value=mn, nonint=nonint)
    print(f'{path}: {len(bl)} blocks, max counts per cell {lib_max:.0f}, min {mn}, non-integer {nonint}')


def logcpm(S):
    return np.log1p(S / S.sum(1, keepdims=True) * 1e6)


def vcc_meta(path):
    with tarfile.open(path) as t:
        return json.loads(t.extractfile('meta.json').read()), {m.name: m.size for m in t.getmembers()}


def check(run, old_npz, old_vcc, pred, array):
    sys.path.insert(0, str(V / 'code/vcc_2026_3c'))
    import argparse
    import anndata as ad
    from scipy.stats import mannwhitneyu
    from vcc3_generate import Cfg, Line, calibrate, generate_one, kd_residual, load_deltas, pick_cells, prediction_ontarget, read_list
    from vcc2_common import seed_for
    run = Path(run); res = {}
    genes = read_list(V / 'data/gene_names.csv', 'gene_name'); G = len(genes); gp = {g: i for i, g in enumerate(genes)}
    targets = read_list(V / 'data/pert_counts.csv', 'target_gene')
    # 1. structure vs the 0.0638 file (generator h5ad obs/var; .vcc container)
    with h5py.File(run / 'prediction.h5ad', 'r') as f:
        obs = read_elem(f['obs']); var = read_elem(f['var']); shape = (len(obs), len(var))
    old = np.load(old_npz)
    new_npz = run / 'pbulk.npz'
    if not new_npz.exists(): pbulk(run / 'prediction.h5ad', new_npz)
    new = np.load(new_npz)
    ref_obs_keys = old['keys']
    om, osz = vcc_meta(old_vcc); nm, nsz = vcc_meta(run / 'prediction.vcc')
    res['1_structure'] = dict(n_obs=shape[0], n_vars=shape[1], gene_order_matches_gene_names=bool(list(var.index.astype(str)) == genes),
                              obs_columns=list(obs.columns), target_col_present='target_gene' in obs,
                              block_keys_identical_to_old=bool(np.array_equal(new['keys'], ref_obs_keys)),
                              targets_per_line_match_pert_counts=all(sorted(obs.loc[obs.context == l, 'target_gene'].astype(str).unique()) == sorted(targets) for l in 'ABC'),
                              cells_per_target=sorted(obs.groupby(['context', 'target_gene'], observed=True).size().unique().tolist()),
                              non_targeting_label_in_h5ad=bool((obs.target_gene.astype(str) == 'non-targeting').any()),
                              vcc_members=list(nsz), old_vcc_members=list(osz), vcc_meta=nm, old_vcc_meta=om,
                              integer_counts=int(new['nonint']) == 0, min_value=float(new['min_value']),
                              max_counts_per_cell=float(new['max_cell_total']))
    res['1_structure']['PASS'] = bool(res['1_structure']['gene_order_matches_gene_names'] and res['1_structure']['block_keys_identical_to_old']
                                      and res['1_structure']['targets_per_line_match_pert_counts'] and res['1_structure']['cells_per_target'] == [400]
                                      and res['1_structure']['integer_counts'] and res['1_structure']['min_value'] >= 0
                                      and res['1_structure']['max_counts_per_cell'] <= 1e6 and list(nsz) == list(osz)
                                      and {k: v for k, v in nm.items() if k != 'nnz'} == {k: v for k, v in om.items() if k != 'nnz'})
    # 2 + 3. realised vs intended pseudobulk (30 knockdowns per line, same seeds) and the null test, per line
    settings = json.loads((run / 'generator_report.json').read_text())['settings']
    cfg = Cfg(argparse.Namespace(**settings))
    res['2_realised_vs_intended'] = {}; res['3_null'] = {}
    keys = list(new['keys']); kpos = {k: i for i, k in enumerate(keys)}
    for name in 'ABC':
        line = Line.from_h5ad(V / f'data/context_{name}.h5ad', genes)
        D = calibrate(load_deltas(Path(pred) / name, targets, genes, array), line, targets, cfg)
        kd = kd_residual(prediction_ontarget(Path(pred) / name, targets), cfg)
        common = pick_cells(line, np.zeros(G), cfg, np.random.default_rng(seed_for(cfg.seed, name, 'base_controls')))
        ratios = []; same = []
        for i in np.random.default_rng(0).choice(len(targets), 30, replace=False):
            t = targets[i]; c = Cfg(argparse.Namespace(**dict(settings, exact_totals=True)))
            X, diag, _ = generate_one(line, D[i], gp.get(t), kd[i], c, pick=common[:400], line_key=name, target_key=t)
            ratios.append(diag['bias_ratio_median'])
            # the exact-totals call draws the same counts before pinning: rebuild the unpinned sum to tie it to the file
            Xu, _, _ = generate_one(line, D[i], gp.get(t), kd[i], cfg, pick=common[:400], line_key=name, target_key=t)
            same.append(bool(np.array_equal(np.asarray(Xu.sum(0)).ravel(), new['sums'][kpos[f'{name}|{t}']])))
        res['2_realised_vs_intended'][name] = dict(median_of_per_kd_median_ratio=float(np.median(ratios)), min=float(np.min(ratios)),
                                                   max=float(np.max(ratios)), file_matches_regeneration=f'{sum(same)}/{len(same)}',
                                                   PASS=bool(0.97 <= np.min(ratios) and np.max(ratios) <= 1.03 and all(same)))
        X0, _, _ = generate_one(line, np.zeros(G), None, None, cfg, pick=common[:400], line_key=name, target_key='__null__')
        ctrl = line.cpm_cells; keep = line.cpm > 5
        lib0 = np.asarray(X0.sum(1)).ravel(); P0 = (X0.multiply(1e6 / lib0[:, None])).tocsc()[:, keep].toarray()
        C0 = ctrl.tocsc()[:, keep].toarray()
        p = mannwhitneyu(P0, C0, axis=0).pvalue; p = np.where(np.isfinite(p), p, 1)
        from scipy.stats import false_discovery_control
        q = false_discovery_control(p, method='bh')
        res['3_null'][name] = dict(genes_tested=int(keep.sum()), significant_bh05=int((q < 0.05).sum()), PASS=bool((q < 0.05).sum() <= 5))
        print(name, res['2_realised_vs_intended'][name], res['3_null'][name], flush=True)
    # 4. diff vs the previous file (per-knockdown log1p CPM pseudobulks)
    A, B = logcpm(old['sums']), logcpm(new['sums']); diff = np.abs(B - A)
    lines = np.array([k.split('|')[0] for k in keys])
    res['4_diff_vs_previous'] = {l: dict(mean_abs_dlog1p_cpm=float(diff[lines == l].mean()),
                                         genes_moved_gt_0_1_median_per_kd=float(np.median((diff[lines == l] > 0.1).sum(1))),
                                         genes_moved_gt_0_1_iqr=[float(x) for x in np.percentile((diff[lines == l] > 0.1).sum(1), [25, 75])])
                                 for l in 'ABC'}
    (run / 'checks_v5.json').write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))


if __name__ == '__main__':
    if sys.argv[1] == 'pbulk': pbulk(sys.argv[2], sys.argv[3])
    else: check(*sys.argv[2:7])

"""Build one eval of the local harness: split definition, context + reference h5ads, real cache, ceiling and bundle.

  ~/venv/bin/python scripts/harness/build.py --eval A          # H1, raw unfiltered file, all knockdowns
  ~/venv/bin/python scripts/harness/build.py --eval B          # K562, 100 held-out kept knockdowns
  ~/venv/bin/python scripts/harness/build.py --eval A_shared   # after A (and B's QC tables): A restricted to K562-kept

Writes lib/harness/eval<name>/: split.json (knockdowns, cells per knockdown, seeds, thinning factor, control split),
split_rows.npz (source rows of every reference / context cell), context.h5ad (context controls), real.h5ad (held-out
knockdowns' real cells, then the scoring controls), real_cache/, ceiling/, bundle/ (cell-eval2 vcc2026, pdex).
Cells: per knockdown a seeded 400 of its real cells (all when fewer); controls: a seeded 2 x 18,400 (all when fewer),
split into two seeded halves (context / scoring); each cell binomially thinned by f = min(1, 20,000 / median control
panel library), seeded per source row; genes reindexed to the 18,533-gene panel. Eval A takes H1 from the RAW file
(lib/h1/H1_full_2025.h5ad: every cell of the 300 knockdowns, no cell or knockdown QC; the QC'd H1 file is training data
only). Eval B takes K562 from lib/qc_cells/replogle_k562_gw_qc.h5ad, kept cells only (guide_qc keep -> knockdown_qc keep).
"""
import argparse
import hashlib
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

sys.path.insert(0, str(Path(__file__).parent))
from common import (V, LIB, CONTROL, SEED, DEPTH, N_CELLS, N_CTRL, G, eval_dir, panel_genes, evaluator, thin,  # noqa: E402
                    read_rows, write_h5ad, log, save_json, seed_for, read_json)

H1_RAW = V / 'lib/h1/H1_full_2025.h5ad'
K562 = V / 'lib/qc_cells/replogle_k562_gw_qc.h5ad'
K562_QC = V / 'results/qc/replogle_k562_gw'
MIN_CELLS = 100                     # vcc3_generate benchmark's default eligibility for a held-out knockdown
N_B = 100


def panel_map(var_src, kind):
    """Source gene -> panel column. H1: by symbol (its index). K562: Ensembl ID -> panel symbol through the v4/v5
    panel table (lib/harness/panel.csv, symbol + gene_id), the same mapping the v5 effect tables use."""
    genes = panel_genes(); col = {g: i for i, g in enumerate(genes)}
    if kind == 'symbol':
        return np.array([col.get(s, -1) for s in var_src.index.astype(str)]), genes
    pt = pd.read_csv(LIB / 'panel.csv', keep_default_na=False)
    id2col = {g: col[s] for s, g in zip(pt.gene, pt.gene_id) if g}
    return np.array([id2col.get(g, -1) for g in var_src.index.astype(str)]), genes


def k562_kept(obs):
    """K562 kept knockdowns and, per QC'd cell, its kept-knockdown symbol ('' = not kept): cell -> kept guide ->
    kept knockdown, as build_effects_v5 does (read-through targets have no Ensembl ID, so guides are the key)."""
    guides = pd.read_csv(K562_QC / 'guide_qc.csv'); kds = pd.read_csv(K562_QC / 'knockdown_qc.csv')
    kept = kds[kds.knockdown_keep.astype(bool)]
    sym_of = dict(zip(kept.target_gene_id.astype(str), kept.target_gene.astype(str)))
    g = guides[guides.keep.astype(bool)]
    guide_kd = {gid: sym_of.get(t, '') for gid, t in zip(g.guide_id.astype(str), g.target_gene_id.astype(str))}
    ctrl = obs['is_control'].to_numpy(bool)
    label = np.array(['' if c else guide_kd.get(x, '') for c, x in zip(ctrl, obs['guide_id'].astype(str))], dtype=object)
    return set(kept.target_gene.astype(str)), label, ctrl


def h1_raw_targets():
    with h5py.File(H1_RAW, 'r') as fh:
        lab = read_elem(fh['obs'])['target_gene'].astype(str)
    vc = lab.value_counts()
    return sorted(t for t, n in vc.items() if t != CONTROL and n >= MIN_CELLS)


def sample_rows(groups, label, rng_key):
    """Per knockdown a seeded N_CELLS of its rows (all when fewer); controls a seeded 2 x N_CTRL, split in two halves."""
    picked = {}
    for t in groups:
        r = np.flatnonzero(label == t)
        rng = np.random.default_rng(seed_for(SEED, 'harness', rng_key, t))
        picked[t] = np.sort(rng.choice(r, N_CELLS, replace=False)) if len(r) > N_CELLS else r
    c = np.flatnonzero(label == CONTROL)
    rng = np.random.default_rng(seed_for(SEED, 'harness', rng_key, 'controls'))
    c = rng.permutation(c)[:2 * N_CTRL]
    half = len(c) // 2
    return picked, np.sort(c[:half]), np.sort(c[half:])              # context half, scoring half


def build_cells(name, path, label, targets, to_panel, genes, extra_obs, rng_key):
    """Read, thin, write context.h5ad and real.h5ad; return the split record."""
    d = eval_dir(name); d.mkdir(parents=True, exist_ok=True)
    picked, ctx_rows, score_rows = sample_rows(targets, label, rng_key)
    pert_rows = np.concatenate([picked[t] for t in targets])
    all_rows = np.unique(np.concatenate([pert_rows, ctx_rows, score_rows]))
    log(f'eval{name}: {len(targets)} knockdowns, {len(pert_rows):,} real cells, {len(ctx_rows):,} context + '
        f'{len(score_rows):,} scoring controls; reading {len(all_rows):,} rows of {path}')
    M = read_rows(path, all_rows, to_panel, len(genes))
    pos = {r: i for i, r in enumerate(all_rows)}
    lib = np.asarray(M.sum(1)).ravel()
    ctrl_lib = lib[[pos[r] for r in np.concatenate([ctx_rows, score_rows])]]
    f = min(1.0, DEPTH / float(np.median(ctrl_lib)))
    log(f'  median control panel library {np.median(ctrl_lib):.0f} -> thinning factor {f:.4f}')
    if f < 1.0:
        for i, r in enumerate(all_rows):
            a, b = M.indptr[i], M.indptr[i + 1]; M.data[a:b] = thin(M.data[a:b], r, f)
        M.eliminate_zeros()
    lib_after = np.asarray(M.sum(1)).ravel()

    def obs_for(rows, labels):
        o = pd.DataFrame(dict(target=labels, target_gene=labels, source_row=rows), index=[f'{name}_{r}' for r in rows])
        for k, v in extra_obs.items(): o[k] = np.asarray(v)[rows]
        return o
    ri = np.array([pos[r] for r in ctx_rows])
    write_h5ad(d / 'context.h5ad', M[ri], obs_for(ctx_rows, [CONTROL] * len(ctx_rows)), genes)
    ref_rows = np.concatenate([pert_rows, score_rows])
    ref_labels = [t for t in targets for _ in picked[t]] + [CONTROL] * len(score_rows)
    write_h5ad(d / 'real.h5ad', M[np.array([pos[r] for r in ref_rows])], obs_for(ref_rows, ref_labels), genes)
    np.savez_compressed(d / 'split_rows.npz', context=ctx_rows, scoring_controls=score_rows,
                        **{f'kd:{t}': picked[t] for t in targets})
    med = lambda rows: float(np.median(lib_after[[pos[r] for r in rows]]))
    return dict(knockdowns=list(targets), n_knockdowns=len(targets), cells_per_knockdown={t: int(len(picked[t])) for t in targets},
                context_controls=int(len(ctx_rows)), scoring_controls=int(len(score_rows)), thinning_factor=f,
                median_panel_library_before=float(np.median(ctrl_lib)), median_panel_library_after_controls=med(np.concatenate([ctx_rows, score_rows])),
                median_panel_library_after_perturbed=med(pert_rows), seed=SEED, seed_keys=f"seed_for({SEED}, 'harness', '{rng_key}', <knockdown | 'controls'>)",
                thinning_seed=f'default_rng([{SEED}, source row])', depth_target=DEPTH, max_cells=N_CELLS, max_controls_per_half=N_CTRL,
                genes_on_panel=int((to_panel >= 0).sum()), source=str(path))


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''): h.update(b)
    return h.hexdigest()


def finish(name, split, threads):
    d = eval_dir(name)
    split['generator'] = {p.name: sha(p) for p in sorted((V / 'code/vcc_2026_3c').glob('*.py'))}
    save_json(d / 'split.json', split)
    for sub in ('real_cache', 'ceiling', 'baseline', 'bundle'):
        if (d / sub).exists(): shutil.rmtree(d / sub)
    ev = evaluator(name, threads=threads)
    log('  cell-eval2: real cache + ceiling'); ev.prepare()
    log('  cell-eval2: competition bundle (generic baseline + replicate anchor)'); ev.build_bundle()
    save_json(d / 'build_invocations.json', dict(commands=ev.invocations))
    log(f'eval{name} ready: {d}')


def build_A(threads):
    with h5py.File(H1_RAW, 'r') as fh:
        obs = read_elem(fh['obs']); var = read_elem(fh['var'])
    to_panel, genes = panel_map(var, 'symbol')
    label = obs['target_gene'].astype(str).to_numpy()
    targets = h1_raw_targets()
    split = build_cells('A', H1_RAW, label, targets, to_panel, genes,
                        dict(guide_id=obs['guide_id'].astype(str).to_numpy(), batch=obs['batch'].astype(str).to_numpy()), 'A')
    split.update(line='H1', rule='raw unfiltered file (no cell or knockdown QC); the QC\'d arc_h1_qc.h5ad is training data only',
                 eligibility=f'every H1 knockdown with >= {MIN_CELLS} raw cells')
    finish('A', split, threads)


def build_B(threads):
    with h5py.File(K562, 'r') as fh:
        obs = read_elem(fh['obs']); var = read_elem(fh['var'])
    to_panel, genes = panel_map(var, 'ensembl')
    kept, label, ctrl = k562_kept(obs)
    label = np.where(ctrl, CONTROL, label)
    shared = sorted(kept & set(h1_raw_targets()))
    rng = np.random.default_rng(seed_for(SEED, 'harness', 'B', 'targets'))
    if len(shared) >= N_B:
        targets = sorted(rng.choice(shared, N_B, replace=False).tolist()); fill = []
    else:
        fill = sorted(rng.choice(sorted(kept - set(shared)), N_B - len(shared), replace=False).tolist()); targets = sorted(shared + fill)
    split = build_cells('B', K562, label, targets, to_panel, genes,
                        dict(guide_id=obs['guide_id'].astype(str).to_numpy(), batch=obs['batch'].astype(str).to_numpy()), 'B')
    split.update(line='K562', rule='QC\'d cells of kept constructs of kept knockdowns (guide_qc.csv -> knockdown_qc.csv)',
                 shared_with_h1_available=len(shared), filled_at_random=fill,
                 target_seed=f"seed_for({SEED}, 'harness', 'B', 'targets')")
    finish('B', split, threads)


def build_A_shared(threads):
    """Eval A's cells for the knockdowns K562 also keeps: same reference rows, same context, its own cache and bundle."""
    import anndata as ad
    src = eval_dir('A'); d = eval_dir('A_shared'); d.mkdir(parents=True, exist_ok=True)
    with h5py.File(K562, 'r') as fh: obs = read_elem(fh['obs'])
    kept, _, _ = k562_kept(obs)
    a = ad.read_h5ad(src / 'real.h5ad')
    lab = a.obs['target_gene'].astype(str).values
    shared = sorted(set(lab) - {CONTROL} & kept)
    keep = np.isin(lab, shared + [CONTROL])
    write_h5ad(d / 'real.h5ad', a.X[np.flatnonzero(keep)], a.obs.iloc[np.flatnonzero(keep)][['target', 'target_gene', 'source_row', 'guide_id', 'batch']],
               a.var_names.astype(str).tolist())
    del a
    shutil.copy(src / 'context.h5ad', d / 'context.h5ad')
    split = read_json(src / 'split.json')
    split.update(knockdowns=shared, n_knockdowns=len(shared), cells_per_knockdown={t: split['cells_per_knockdown'][t] for t in shared},
                 subset_of='evalA', subset_rule='H1 knockdowns that K562 also keeps (knockdown_qc.csv keep)')
    finish('A_shared', split, threads)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--eval', required=True, choices=['A', 'B', 'A_shared'])
    p.add_argument('--threads', type=int, default=8)
    a = p.parse_args()
    dict(A=build_A, B=build_B, A_shared=build_A_shared)[a.eval](a.threads)

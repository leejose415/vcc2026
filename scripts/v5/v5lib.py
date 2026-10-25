"""vcc_2026_4 v5 queue (21 Oct): output arms on dev_nopie's zero-shot H1 predictions and on final_nopie's A-C predictions.

Notation: delta_k = predicted log2FC of knockdown k; dbar = the line's mean delta over its (300) knockdowns, per gene;
specific part = delta_k - dbar. Arms: delta'_k = beta * dbar + s * (delta_k - dbar) on the run's 5,000 genes, optionally
plus an extension to the other panel genes (beta x the mean over source lines of each line's mean knockdown lfc,
genes measured in >= 3 source lines, else 0). Everything downstream (gate, own gene, composition) is the generator's.
"""
import json
import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

V = Path.home() / 'vcc'
RUN = V / 'lib/vcc_2026_4/runs/dev_nopie'
FINAL = V / 'lib/vcc_2026_4/runs/final_nopie'
H1_SRC_DIR = RUN / 'h1_pred_g3'                         # the scored H1 predictions (tf_top5000 = every run gene)
H1_DIR = RUN / 'h1_pred_v5'
FINAL_SRC_DIR = FINAL / 'pred_top5000'
EFFECTS = V / 'lib/effects_v4'
# effects_v4 file per source line; H1 (arc_h1) is a source only for the challenge lines A-C
SOURCE_FILES = dict(K562='replogle_k562_gw', RPE1='replogle_rpe1', HepG2='nadig_hepg2', Jurkat='nadig_jurkat',
                    HCT116='orion_hct116', HEK293T='orion_hek293t', CD4_Rest='marson_cd4_rest',
                    CD4_Stim8hr='marson_cd4_stim8hr', CD4_Stim48hr='marson_cd4_stim48hr', H1='arc_h1')


def panel_genes():
    return pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy()


def arm(delta, s=1.0, beta=1.0, ext=None, per_gene=None):
    """delta: (K, n_run) run-gene predictions of one line. Returns run-gene values and the extension (or None).
    per_gene: optional (n_run,) multiplier on the specific part (Part D); s may be a (K,) vector (Part C)."""
    dbar = delta.mean(0, keepdims=True)
    spec = delta - dbar
    s = np.asarray(s, np.float64).reshape(-1, 1) if np.ndim(s) else s
    if per_gene is not None: spec = spec * per_gene[None, :]
    return (beta * dbar + s * spec).astype(np.float32)


def full_array(run_vals, gcols, G, ext_vals=None):
    """(K, G) panel array: run genes from run_vals, other genes from ext_vals (a (G,) vector, same for every knockdown)
    or 0."""
    D = np.zeros((run_vals.shape[0], G), np.float32)
    if ext_vals is not None:
        out = np.ones(G, bool); out[gcols] = False
        D[:, out] = ext_vals[out][None, :]
    D[:, gcols] = run_vals
    return D


def source_line_means(sources, genes, min_lines=3):
    """Per panel gene: the mean over source lines of that line's mean knockdown lfc (NaN-aware), kept where >= min_lines
    lines measure the gene, else 0. Also returns the per-line means and the count."""
    G = len(genes); gp = {g: i for i, g in enumerate(genes)}; per = []
    for src in sources:
        with h5py.File(EFFECTS / f'{SOURCE_FILES[src]}.h5ad', 'r') as f:
            var = read_elem(f['var']); vg = var.index.astype(str)
            n_rows, Gs = f['layers/lfc'].shape; s1 = np.zeros(Gs); n = np.zeros(Gs)
            for a in range(0, n_rows, 2000):
                L = f['layers/lfc'][a:a + 2000].astype(np.float64); ok = np.isfinite(L)
                s1 += np.where(ok, L, 0).sum(0); n += ok.sum(0)
        m = np.full(G, np.nan)
        with np.errstate(invalid='ignore', divide='ignore'): mm = s1 / n
        idx = np.array([gp.get(g, -1) for g in vg]); ok = (idx >= 0) & (n >= 30)
        m[idx[ok]] = mm[ok]; per.append(m); print(f'  line mean: {src} ({int(ok.sum())} genes)', flush=True)
    P = np.stack(per); cnt = np.isfinite(P).sum(0)
    with np.errstate(invalid='ignore'): mean = np.nanmean(P, 0)
    return np.where(cnt >= min_lines, np.nan_to_num(mean), 0.0), P, cnt


def write_dir(dst, src, arrays, heldout, sources):
    """A prediction folder for vcc3_generate.py: genes / targets / ontarget copied from src, the given arrays, report."""
    dst = Path(dst); dst.mkdir(parents=True, exist_ok=True)
    for f in ('genes.npy', 'targets.npy', 'ontarget_lfc.npy'): shutil.copy(Path(src) / f, dst / f)
    for name, D in arrays.items(): np.save(dst / f'{name}.npy', D.astype(np.float32))
    rep = dst / 'prediction_report.json'
    old = json.loads(rep.read_text())['models'] if rep.exists() else []
    models = old + [m for m in arrays if m not in old]
    rep.write_text(json.dumps(dict(schema='vcc_2026_3_ridge_v1', heldout_background=heldout, source_backgrounds=sources,
                                   models=models, note='v5 output arms: beta * dbar + s * (delta - dbar) (+ extension)'), indent=1))
    return models

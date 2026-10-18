"""H1 SE-calibration anomaly: which genes pass |z| > 3 most often across the n = 300 control-only draws? Read only.

Replays scripts/se_calibration.py exactly (same seed, same draw order: 200 draws at n = 30, 200 at n = 100, then the 200
at n = 300), and for the n = 300 draws counts per gene how often |z| > 3 (and how often it also has |lfc| > 0.5, i.e.
counts towards the noise-aware total). For the most frequent genes: the max per-cell count among H1 controls vs its 99.9th
percentile, the share of controls expressing it, its control mean CPM, and a gene-type flag (mitochondrial, ribosomal,
hemoglobin, stress / immediate-early). Changes nothing; writes results/qc/arc_h1/se_anomaly_genes.csv.
"""
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

sys.path.insert(0, str(Path(__file__).parent))
from qc_replogle_umi import SOURCES                                                        # noqa: E402
from se_calibration import CHUNK, LFC_MIN, PSEUDO, REPS, SE_FLOOR, SEED, SIZES, Z_MIN, stats   # noqa: E402

TOP = 25
FLAGS = [('mitochondrial', r'^MT-'), ('ribosomal', r'^RP[LS]\d|^RPLP|^RPSA'), ('hemoglobin', r'^HB[ABDEGMQZ]\d?$'),
         ('stress / immediate-early', r'^HSP[ABDEH]|^HSP90|^DNAJ[AB]|^FOS|^JUN|^EGR[1-3]$|^ATF3$|^IER[235]|^DDIT3$|^MT1|^MT2A$|^HSPA')]


def flag(symbol):
    hits = [name for name, rx in FLAGS if re.search(rx, symbol)]
    return ', '.join(hits) if hits else ''


def main():
    _, out, path, _ = SOURCES['H1']
    parts, libs, rows_ctrl = [], [], []
    with h5py.File(path, 'r') as f:
        obs = read_elem(f['obs']); var = read_elem(f['var']); X = f['X']; n_obs = len(obs); n_vars = int(X.attrs['shape'][1])
        ctrl = obs['is_control'].to_numpy(bool); lib = obs['n_umi'].to_numpy(float); ip = X['indptr']
        for s in range(0, n_obs, CHUNK):
            e = min(s + CHUNK, n_obs); cm = np.flatnonzero(ctrl[s:e])
            if len(cm):
                p = ip[s:e + 1]; a, b = int(p[0]), int(p[-1])
                M = sp.csr_matrix((X['data'][a:b], X['indices'][a:b], p - a), shape=(e - s, n_vars))[cm]
                parts.append((sp.diags(1e6 / lib[s:e][cm]) @ M).astype(np.float32)); libs.append(lib[s:e][cm]); rows_ctrl.append(s + cm)
    C = sp.vstack(parts).tocsr(); del parts; depth = np.concatenate(libs); n_ctrl = C.shape[0]   # float32 CPM: the memory limit
    total1 = np.asarray(C.sum(0, dtype=np.float64)).ravel(); total2 = np.asarray(C.multiply(C).sum(0, dtype=np.float64)).ravel(); total_depth = depth.sum()
    measured = np.asarray((C > 0).sum(0)).ravel() > 0
    rng = np.random.default_rng(SEED)
    hit_z = np.zeros(n_vars, int); hit_both = np.zeros(n_vars, int)
    for n in SIZES:
        for _ in range(REPS):
            pick = rng.choice(n_ctrl, n, replace=False)
            if n != 300: continue
            G = C[pick].astype(np.float64); g1 = np.asarray(G.sum(0)).ravel(); g2 = np.asarray(G.multiply(G).sum(0)).ravel(); gd = depth[pick].sum()
            m1, rv1 = stats(g1, g2, n, gd); m0, rv0 = stats(total1 - g1, total2 - g2, n_ctrl - n, total_depth - gd)
            lfc = np.log2((m1 + PSEUDO) / (m0 + PSEUDO)); se = np.maximum(np.sqrt(rv1 + rv0) / np.log(2), SE_FLOOR)
            big = (np.abs(lfc / se) > Z_MIN) & measured
            hit_z += big; hit_both += big & (np.abs(lfc) > LFC_MIN)
    symbols = var.index.astype(str).to_numpy()
    order = np.argsort(-hit_both * 1000 - hit_z)[:TOP]
    del C
    keep_rows = np.concatenate(rows_ctrl); K = np.zeros((n_ctrl, len(order)), np.float32); k = 0
    with h5py.File(path, 'r') as f:                                   # second pass: raw counts of the top genes only
        X = f['X']; ip = X['indptr']
        for s in range(0, n_obs, CHUNK):
            e = min(s + CHUNK, n_obs); cm = keep_rows[(keep_rows >= s) & (keep_rows < e)] - s
            if len(cm):
                p = ip[s:e + 1]; a, b = int(p[0]), int(p[-1])
                M = sp.csr_matrix((X['data'][a:b], X['indices'][a:b], p - a), shape=(e - s, n_vars))[cm][:, order]
                K[k:k + len(cm)] = M.toarray(); k += len(cm)
    rows = []
    for jj, j in enumerate(order):
        col = K[:, jj]
        rows.append(dict(gene=symbols[j], draws_abs_z_gt3=int(hit_z[j]), draws_noise_aware=int(hit_both[j]), of_draws=REPS,
                         max_count=float(col.max()), p999_count=float(np.quantile(col, 0.999)),
                         max_over_p999=float(col.max() / max(np.quantile(col, 0.999), 1)), cells_expressing=float((col > 0).mean()),
                         control_mean_cpm=float(total1[j] / n_ctrl), flag=flag(symbols[j])))
    t = pd.DataFrame(rows)
    t.to_csv(out / 'se_anomaly_genes.csv', index=False)
    expected = 0.0027 * REPS
    print(f'H1 controls: {n_ctrl:,}; n = 300 draws: {REPS}; a gene with a correct SE passes |z| > 3 in ~{expected:.1f} draws')
    print(f'genes passing |z| > 3 in >= 10 of {REPS} draws: {int((hit_z >= 10).sum())}; in >= 50: {int((hit_z >= 50).sum())}')
    pd.set_option('display.width', 220); print(t.round(3).to_string(index=False))
    print(f'wrote {out / "se_anomaly_genes.csv"}')


if __name__ == '__main__':
    main()

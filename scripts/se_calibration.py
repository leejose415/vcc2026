"""Control-only noise and SE calibration: is the effect tables' SE right? Read only; writes a JSON and a CSV.

For a cell-QC'd source (lib/qc_cells/<source>_qc.h5ad), only its control cells ('non-targeting'):
for each group size n in {30, 100, 300}, 200 times: draw n random controls as a pseudo-knockdown; the remaining
controls are its comparison group. lfc and se exactly as scripts/build_effect_tables.py:
  per-cell CPM = counts / n_umi * 1e6;  lfc = log2((mean CPM + 0.5) / (comparison mean CPM + 0.5))
  se = sqrt(relvar_group + relvar_comparison) / ln 2, relvar = max(var / n, (mean + u) * u) / (mean + 0.5)^2,
       u = 1e6 / (the group's total depth), var = sample variance of the per-cell CPM; floored at 0.05
over the measured genes (>= 1 count in this source's controls). Recorded: z = lfc / se for every measured gene, and the
noise-aware downstream count (|lfc| > 0.5 and |z| > 3).

Report per n: SD of z (about 1 if the SE is right), share of genes with |z| > 3 (about 0.27% if Gaussian and right),
median and mean noise-aware count (the noise floor).
Writes results/qc/<source>/se_calibration.json and se_calibration_reps.csv.
Usage: ~/venv/bin/python scripts/se_calibration.py --dataset K562|RPE1|HepG2|Jurkat|H1
"""
import argparse
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

sys.path.insert(0, str(Path(__file__).parent))
from qc_replogle_umi import SOURCES                                                        # noqa: E402

SIZES, REPS, SEED, PSEUDO, SE_FLOOR, LFC_MIN, Z_MIN, CHUNK = (30, 100, 300), 200, 2026, 0.5, 0.05, 0.5, 3.0, 50000


def stats(s1, s2, n, depth):
    mean = s1 / n; var = np.maximum(s2 / n - mean ** 2, 0) * n / max(n - 1, 1)
    u = 1e6 / max(depth, 1)
    return mean, np.maximum(var / n, (np.maximum(mean, 0) + u) * u) / (mean + PSEUDO) ** 2


def main(dataset):
    _, out, path, title = SOURCES[dataset]
    t0 = time.time(); parts = []; libs = []
    with h5py.File(path, 'r') as f:
        obs = read_elem(f['obs']); X = f['X']; n_obs = len(obs); n_vars = int(X.attrs['shape'][1])
        ctrl = obs['is_control'].to_numpy(bool); lib = obs['n_umi'].to_numpy(float); ip = X['indptr']
        for s in range(0, n_obs, CHUNK):
            e = min(s + CHUNK, n_obs); cm = np.flatnonzero(ctrl[s:e])
            if len(cm):
                p = ip[s:e + 1]; a, b = int(p[0]), int(p[-1])
                M = sp.csr_matrix((X['data'][a:b], X['indices'][a:b], p - a), shape=(e - s, n_vars))[cm]
                parts.append(sp.diags(1e6 / lib[s:e][cm]) @ M.astype(np.float64)); libs.append(lib[s:e][cm])
    C = sp.vstack(parts).tocsr(); depth = np.concatenate(libs); n_ctrl = C.shape[0]
    total1 = np.asarray(C.sum(0)).ravel(); total2 = np.asarray(C.multiply(C).sum(0)).ravel(); total_depth = depth.sum()
    measured = np.asarray((C > 0).sum(0)).ravel() > 0
    print(f'{title}: {n_ctrl:,} controls, {int(measured.sum()):,} of {n_vars:,} genes measured ({time.time() - t0:.0f}s to load)', flush=True)
    rng = np.random.default_rng(SEED); rows, summary = [], {}
    for n in SIZES:
        sd_acc = []; big = 0; total = 0
        for rep in range(REPS):
            pick = rng.choice(n_ctrl, n, replace=False); G = C[pick]
            g1 = np.asarray(G.sum(0)).ravel(); g2 = np.asarray(G.multiply(G).sum(0)).ravel(); gd = depth[pick].sum()
            m1, rv1 = stats(g1, g2, n, gd)
            m0, rv0 = stats(total1 - g1, total2 - g2, n_ctrl - n, total_depth - gd)
            lfc = np.log2((m1 + PSEUDO) / (m0 + PSEUDO))[measured]
            se = np.maximum(np.sqrt(rv1 + rv0) / np.log(2), SE_FLOOR)[measured]
            z = lfc / se
            noise_aware = int(((np.abs(lfc) > LFC_MIN) & (np.abs(z) > Z_MIN)).sum())
            sd_acc.append(z); big += int((np.abs(z) > Z_MIN).sum()); total += z.size
            rows.append(dict(n=n, rep=rep, sd_z=float(z.std()), share_abs_z_gt3=float((np.abs(z) > Z_MIN).mean()), noise_aware=noise_aware,
                             se_at_floor=float((se == SE_FLOOR).mean())))
        zs = np.concatenate(sd_acc); r = pd.DataFrame([x for x in rows if x['n'] == n])
        summary[n] = dict(sd_z=float(zs.std()), share_abs_z_gt3=float(big / total), median_noise_aware=float(r.noise_aware.median()),
                          mean_noise_aware=float(r.noise_aware.mean()), p95_noise_aware=float(r.noise_aware.quantile(0.95)),
                          share_se_at_floor=float(r.se_at_floor.mean()), sd_z_reps_median=float(r.sd_z.median()))
        print(f'  n={n}: SD(z) {summary[n]["sd_z"]:.3f}; |z| > 3: {summary[n]["share_abs_z_gt3"]:.3%}; noise-aware count median '
              f'{summary[n]["median_noise_aware"]:.0f}, mean {summary[n]["mean_noise_aware"]:.1f} ({time.time() - t0:.0f}s)', flush=True)
    pd.DataFrame(rows).to_csv(out / 'se_calibration_reps.csv', index=False)
    (out / 'se_calibration.json').write_text(json.dumps(dict(source=out.name, controls=n_ctrl, measured_genes=int(measured.sum()),
                                                             reps=REPS, seed=SEED, sizes=list(SIZES), results={str(k): v for k, v in summary.items()}), indent=1) + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--dataset', required=True, choices=list(SOURCES))
    main(p.parse_args().dataset)

"""F1 A-C prediction folders: Model 2 v4 (H1 as a source) -> scorer space on each line's control share -> the full-H1
calibration (sign logistic + k_q) -> log2 factors; plus the called-gene masks (top 2,000 by re_z among the line's
genes > 5 CPM, own gene out). Arrays f1_cal, f1_konly, called_top2000 in lib/vcc_2026_5/final/f1/<line>/.

  cd ~/vcc && PYTHONPATH=scripts/v5:analysis/model7:analysis/final ~/venv/bin/python analysis/final/f1_abc.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from calib import features, groups, to_scorer
from m7_score import to_log2_factor
from v5lib import write_dir

V = Path.home() / 'vcc'
F = V / 'analysis/final'
sys.path.insert(0, str(V / 'code/vcc_2026_3c'))
from vcc3_generate import Line, read_list  # noqa: E402

ABCSRC = V / 'lib/vcc_2026_5/final_nopie_pred_top5000'
genes = read_list(V / 'data/gene_names.csv', 'gene_name'); gp = {g: i for i, g in enumerate(genes)}
z = np.load(F / 'm2v4_abc.npz', allow_pickle=False); tg = z['targets'].astype(str)
cf = np.load(F / 'calib_full.npz')
out = V / 'lib/vcc_2026_5/final/f1'; rep = {}
for c in 'ABC':
    assert list(np.load(ABCSRC / c / 'targets.npy').astype(str)) == list(tg)
    line = Line.from_h5ad(V / f'data/context_{c}.h5ad', genes)
    share = line.share; q = groups(np.log1p(5e4 * share), share > 0)
    pred = to_scorer(z['re'], share)
    p = np.full(pred.shape, 0.5, np.float32)
    for i in range(len(tg)):
        X = features(z['sign_agree'][i], z['n_src'][i], z['re_z'][i], q)
        p[i] = np.where(pred[i] != 0, 1 / (1 + np.exp(-(X @ cf['coef'] + cf['intercept'][0]))), 0.5)
    cal = cf['k_cal'][q][None] * pred * (2 * p - 1); konly = cf['k_only'][q][None] * pred
    zr = np.where(line.gate[None], z['re_z'], -1)
    for i, t in enumerate(tg):
        if t in gp: zr[i, gp[t]] = -1
    called = np.zeros(pred.shape, bool); top = np.argpartition(-zr, 1999, axis=1)[:, :2000]
    called[np.arange(len(tg))[:, None], top] = True; called &= zr > 0
    write_dir(out / c, ABCSRC / c, {'f1_cal': to_log2_factor(cal, share), 'f1_konly': to_log2_factor(konly, share), 'called_top2000': called},
              'none (final predictions)', ['K562', 'RPE1', 'HepG2', 'Jurkat', 'HCT116', 'HEK293T', 'CD4_Stim48hr', 'H1'])
    rep[c] = dict(called_median=float(np.median(called.sum(1))), mean_abs_cal_expressed=float(np.abs(cal[:, line.gate]).mean()),
                  mean_abs_konly_expressed=float(np.abs(konly[:, line.gate]).mean()))
base = json.loads((V / 'results/vcc_2026_4_final_nopie_top5000/generator_report.json').read_text())['settings']
s = dict(base, call='z', z=0.0, min_effect=0.0, alpha=1.0, detect_floor=False, match_pool=True)
for arr in ('f1_cal', 'f1_konly'):
    (out / f'settings_{arr}.json').write_text(json.dumps(dict(schema='vcc_2026_3_ridge_v1', settings=s, delta_array=arr), indent=1))
print(json.dumps(rep, indent=1))

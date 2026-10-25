"""F1 emission check on 10 H1 knockdowns (f1_konly array, match_pool + exact_totals, narrowing m in {1, 0.5, 0.25}):
achieved vs intended (a) pooled profile sum X_g / sum L = control pool share x 2^D and (b) per-cell normalised mean
mean_i(X_ig / L_i) = control per-cell mean x 2^D, on the genes > 5 CPM in the controls (own gene out).

  cd ~/vcc/code/vcc_2026_3c && PYTHONPATH=~/vcc/scripts/v5:~/vcc/analysis/final ~/venv/bin/python ~/vcc/analysis/final/f1_check.py
"""
import argparse
import copy
import json
from pathlib import Path

import numpy as np

import gen_f1
from pb import H1

V = Path.home() / 'vcc'; F = V / 'analysis/final'
G = gen_f1.G
h = H1(); gen_f1.PRED_DIRS['*'] = str(F / 'h1_pred')
D = np.load(F / 'h1_pred/f1_konly.npy'); tg = list(np.load(F / 'h1_pred/targets.npy').astype(str)); pos = {t: i for i, t in enumerate(tg)}
cfg0 = G.Cfg(argparse.Namespace(**dict(vars(h.cfg), match_pool=True, exact_totals=True)))
common = G.pick_cells(h.line, np.zeros(h.G), h.cfg, np.random.default_rng(G.seed_for(2026, 'H1', 'base_controls')))
pick10 = np.random.default_rng(1).choice(h.T, 10, replace=False)
Dc = G.calibrate(D[[pos[t] for t in h.targets]], h.line, h.targets, cfg0)
share = h.line.share; cell_mean = h.line.cpm / 1e6
out = []
for m in (1.0, 0.5, 0.25):
    cfg = copy.copy(cfg0); cfg.narrow_m = m; cfg.called_from = 'called_top2000' if m < 1 else ''
    for i in pick10:
        t = h.targets[i]; j = h.own[i] if h.own[i] >= 0 else None
        X, diag, planned = gen_f1.generate_one(h.line, Dc[i], j, h.kd[i], cfg, pick=common[:400], line_key='H1', target_key=t)
        X = X.toarray().astype(np.float64); L = X.sum(1)
        keep = h.line.gate.copy()
        if j is not None: keep[j] = False
        f = np.exp2(Dc[i].astype(np.float64))
        pooled = X.sum(0) / L.sum(); intended_p = share * f; intended_p /= intended_p[share > 0].sum()   # composition-projected
        per_cell = (X / L[:, None]).mean(0); intended_c = cell_mean * f
        e_p = np.abs(np.log2((pooled[keep] + 1e-12) / (intended_p[keep] + 1e-12)))
        e_c = np.abs(np.log2((per_cell[keep] + 1e-12) / (intended_c[keep] + 1e-12)))
        sd_ratio = np.nanmedian(((X / L[:, None])[:, keep].std(0) + 1e-12) / (h.line.cpm_cells[common[:400]][:, keep].toarray().std(0) / 1e6 + 1e-12))
        out.append(dict(m=m, target=t, pooled_median_abs_log2_err=float(np.median(e_p)), pooled_max=float(e_p.max()),
                        percell_median_abs_log2_err=float(np.median(e_c)), percell_max=float(e_c.max()),
                        narrowed=int(diag.get('narrowed_genes', 0)), median_sd_ratio_vs_controls=float(sd_ratio)))
import pandas as pd
df = pd.DataFrame(out); df.to_csv(F / 'f1_emission_check.csv', index=False)
print(df.groupby('m')[['pooled_median_abs_log2_err', 'pooled_max', 'percell_median_abs_log2_err', 'percell_max', 'narrowed', 'median_sd_ratio_vs_controls']].median().round(4).to_string())

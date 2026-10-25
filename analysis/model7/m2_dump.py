"""Model 2 on effects_v4 for H1 (reference arm): random-effects mean over the seven Model 2 lines (weight 1 / (se^2 +
tau2_g), tau2 DerSimonian-Laird, own gene excluded), log2 mean-CPM-ratio space, as analysis/v6lin/v6lin.py builds it.

  cd ~/vcc && PYTHONPATH=analysis/v6lin:scripts/v5 ~/venv/bin/python analysis/model7/m2_dump.py
"""
import numpy as np
import pandas as pd
from v6lin import LINES, Line, V, gather, tau2

genes = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy()
lines = {k: Line(k, s, genes) for k, s in LINES.items()}
with np.errstate(all='ignore'):
    h1 = Line('H1', 'arc_h1', genes)
panel = np.ones(len(genes), bool)
for l in lines.values(): panel &= l.measured
t2 = np.nan_to_num(tau2(lines, genes), nan=0.0)
g = gather(list(h1.tg), list(lines), lines, {s: 1.0 for s in lines}, panel, np.zeros((panel.sum(), 1), np.float32), t2, genes)
np.savez(V / 'analysis/model7/h1_model2_effects_v4.npz', tg=h1.tg, re=np.nan_to_num(g['re']).astype(np.float32))
print('saved', g['re'].shape, int(g['has'].sum()))

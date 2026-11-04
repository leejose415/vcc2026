"""Prediction folders + generator settings for a v5 leaderboard file from final_nopie's A-C predictions.

  python build_file.py <out dir> --s S [--beta B] [--ext] [--floor] [--topn N]
Writes <out>/{A,B,C}/ (array 'v5', plus genes / targets / ontarget_lfc from final_nopie/pred_top5000) and
<out>/settings.json: the 0.0638 file's generator settings with only detect_floor changed (off unless --floor).
--topn: keep per knockdown only the N run genes with the highest p_de (after the arm), as the H1 top-N arms.
--ext: genes outside the 5,000 get beta x the per-gene mean over source lines (H1 included) of each line's mean
knockdown lfc, >= 3 lines (cache from the H1 VM: line_means_abcsources.npz).
"""
import argparse
import json
from pathlib import Path

import numpy as np

from a1b_arrays import topn
from v5lib import arm, full_array, panel_genes, write_dir

V = Path.home() / 'vcc'
SRC = V / 'lib/vcc_2026_5/final_nopie_pred_top5000'      # copy of lib/vcc_2026_4/runs/final_nopie/pred_top5000
OLD = V / 'results/vcc_2026_4_final_nopie_top5000/generator_report.json'

ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--s', type=float, required=True)
ap.add_argument('--beta', type=float, default=1.0); ap.add_argument('--ext', action='store_true'); ap.add_argument('--floor', action='store_true')
ap.add_argument('--topn', type=int)
ap.add_argument('--line-means', default=str(V / 'lib/vcc_2026_5/line_means_abcsources.npz'))
a = ap.parse_args(); out = Path(a.out)
genes = panel_genes(); G = len(genes)
ext = a.beta * np.load(a.line_means)['mean'] if a.ext else None
summary = {}
for line in 'ABC':
    z = np.load(SRC / line / 'raw_predictions.npz'); delta, gcols = z['delta'].astype(np.float64), z['gene_cols']
    assert np.array_equal(z['targets'].astype(str), np.load(SRC / line / 'targets.npy').astype(str))
    old = np.load(SRC / line / 'tf_top5000.npy')
    assert np.abs(old[:, gcols] - delta).max() < 1e-6                     # the 0.0638 file's array = delta on run genes
    vals = arm(delta, s=a.s, beta=a.beta)
    if a.topn: vals = topn(vals, z['p'], a.topn)
    D = full_array(vals, gcols, G, ext)
    rep = json.loads((SRC / line / 'prediction_report.json').read_text())
    write_dir(out / line, SRC / line, {'v5': D}, rep['heldout_background'], rep['source_backgrounds'])
    summary[line] = dict(mean_abs_change_vs_0638_array=float(np.abs(D - old).mean()), nonzero_genes=int((D != 0).any(0).sum()))
settings = json.loads(OLD.read_text())['settings']; settings['detect_floor'] = bool(a.floor)
(out / 'settings.json').write_text(json.dumps(dict(schema='vcc_2026_3_ridge_v1', settings=settings, delta_array='v5',
                                                   arm=dict(s=a.s, beta=a.beta, ext=a.ext, floor=a.floor, topn=a.topn)), indent=1))
print(json.dumps(dict(arm=vars(a), lines=summary), indent=1))

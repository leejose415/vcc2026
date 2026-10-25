"""F0: Model 2 v4 + the v3 (+0.091 leaderboard) generator recipe + totals matching, floor on.

  cd ~/vcc && PYTHONPATH=scripts/v5 ~/venv/bin/python analysis/final/f0.py
Recipe of the v3 file (vcc_2026_3, commit c208334, generator_report in gs://vcc2026/vm/results/vcc_2026_3_ridge):
call topn, top_n 2000 by re_z among genes > 5 CPM, detection floor with z_req 5.0, alpha 1, up_w 0.75, pi 0.
Added here: match_pool (each knockdown's picked control cells corrected to the full control pool) and exact_totals
(pseudobulk totals pinned to the expected values); both are CLI/settings switches of vcc3_generate.py.
Writes the H1 prediction folder (analysis/final/h1_pred, arrays m2v4 / m2v4_z / m2v4_se), its variants file, and the
A-C prediction folders + settings (lib/vcc_2026_5/final/f0/).
"""
import json
from pathlib import Path

import numpy as np

from v5lib import write_dir

V = Path.home() / 'vcc'
F = V / 'analysis/final'
H1SRC = V / 'lib/vcc_2026_4/runs/dev_nopie/h1_pred_g3'
ABCSRC = V / 'lib/vcc_2026_5/final_nopie_pred_top5000'
SOURCES = ['K562', 'RPE1', 'HepG2', 'Jurkat', 'HCT116', 'HEK293T', 'CD4_Stim48hr']

z = np.load(F / 'm2v4_h1.npz', allow_pickle=False)
assert list(z['targets'].astype(str)) == list(np.load(H1SRC / 'targets.npy').astype(str))
arr = {'m2v4': z['re'], 'm2v4_z': z['re_z'], 'm2v4_se': np.where(np.isfinite(z['re_se']), z['re_se'], 1e6).astype(np.float32)}
write_dir(F / 'h1_pred', H1SRC, arr, 'H1', SOURCES)
v3 = dict(call='topn', top_n=2000, detect_floor=True, z_req=5.0, alpha=1.0, z=0.0, min_effect=0.0)
(F / 'variants_f0.json').write_text(json.dumps([dict(v3, method='m2v4', label='f0_m2v4_top2000_floor_z5')], indent=1))
za = np.load(F / 'm2v4_abc.npz', allow_pickle=False); tg = za['targets'].astype(str)
out = V / 'lib/vcc_2026_5/final/f0'
for c in 'ABC':
    assert list(np.load(ABCSRC / c / 'targets.npy').astype(str)) == list(tg)
    write_dir(out / c, ABCSRC / c, {'m2v4': za['re'], 'm2v4_z': za['re_z'],
                                    'm2v4_se': np.where(np.isfinite(za['re_se']), za['re_se'], 1e6).astype(np.float32)},
              'none (final predictions)', SOURCES + ['H1'])
rep = json.loads((V / 'results/vcc_2026_4_final_nopie_top5000/generator_report.json').read_text())['settings']
s = dict(rep, call='topn', top_n=2000, detect_floor=True, z_req=5.0, alpha=1.0, z=0.0, min_effect=0.0, match_pool=True)
(out / 'settings.json').write_text(json.dumps(dict(schema='vcc_2026_3_ridge_v1', settings=s, delta_array='m2v4'), indent=1))
print('written', F / 'h1_pred', out)

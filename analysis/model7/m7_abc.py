"""Model 7 A-C prediction folders + generator settings (floor on = FILE 1 settings, floor off).

  cd ~/vcc && PYTHONPATH=scripts/v5 ~/venv/bin/python analysis/model7/m7_abc.py [variant]   (default model7)
Scorer-space predictions (m7_fit.py final) -> log2 factors on each line's control share (all its control cells, the
submission generator's pool); targets and ontarget_lfc from final_nopie's folders, as FILE 1.
"""
import json
import sys
from pathlib import Path

import numpy as np

from m7_score import to_log2_factor
from v5lib import write_dir

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc_2026_3c'))
from vcc3_generate import Line, read_list  # noqa: E402

variant = sys.argv[1] if len(sys.argv) > 1 else 'model7'
genes = read_list(V / 'data/gene_names.csv', 'gene_name')
z = np.load(V / 'analysis/model7/abc_predictions_scorer_space.npz', allow_pickle=True)
tg = z['tg'].astype(str)
SRC = V / 'lib/vcc_2026_4/runs/final_nopie/pred_top5000'
out = V / f'lib/vcc_2026_5/m7_pred/{variant}'
for c in 'ABC':
    assert list(np.load(SRC / c / 'targets.npy').astype(str)) == list(tg)
    line = Line.from_h5ad(V / f'data/context_{c}.h5ad', genes)
    D = to_log2_factor(z[f'{c}_{variant}'], line.share)
    write_dir(out / c, SRC / c, {'m7': D}, 'none (final predictions)', ['K562', 'RPE1', 'HepG2', 'Jurkat', 'HCT116', 'HEK293T', 'CD4_Rest', 'CD4_Stim8hr', 'CD4_Stim48hr', 'H1'])
    print(c, 'mean |log2 factor| on expressed genes', float(np.abs(D[:, line.gate]).mean()))
rep = json.loads((V / 'results/vcc_2026_4_final_nopie_top5000/generator_report.json').read_text())['settings']
for floor in (True, False):
    s = dict(rep, detect_floor=floor)
    (out / f'settings_{"floor" if floor else "nofloor"}.json').write_text(json.dumps(dict(schema='vcc_2026_3_ridge_v1', settings=s, delta_array='m7'), indent=1))
print('wrote', out)

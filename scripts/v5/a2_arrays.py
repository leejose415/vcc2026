"""A2 cell-based check: s = 2, beta = 1.5 (no extension), detection floor; against arm 7 (s = 2, beta = 1, floor)."""
import json
import numpy as np
from v5lib import H1_SRC_DIR, H1_DIR, RUN, arm, full_array, panel_genes, write_dir

genes = panel_genes(); z = np.load(H1_SRC_DIR / 'raw_predictions.npz')
sources = json.loads((H1_SRC_DIR / 'prediction_report.json').read_text())['source_backgrounds']
arrays = {'v5_s2_b1.5': full_array(arm(z['delta'].astype(np.float64), s=2, beta=1.5), z['gene_cols'], len(genes))}
print(write_dir(H1_DIR, H1_SRC_DIR, arrays, 'H1', sources))
v = [dict(method='v5_s2_b1.5', call='z', z=0.0, min_effect=0.0, alpha=1.0, detect_floor=True, label='a2_s2_b1.5_floor')]
(RUN / 'h1_variants_v5a2.json').write_text(json.dumps(v, indent=1))

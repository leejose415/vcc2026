"""A1: H1 arrays s = 1, 1.5, 2 (no shared change) from dev_nopie's saved zero-shot H1 predictions, plus variants files."""
import json
import numpy as np
from v5lib import H1_SRC_DIR, H1_DIR, RUN, arm, full_array, panel_genes, write_dir

genes = panel_genes(); G = len(genes)
z = np.load(H1_SRC_DIR / 'raw_predictions.npz'); delta, gcols = z['delta'].astype(np.float64), z['gene_cols']
old = np.load(H1_SRC_DIR / 'tf_top5000.npy')
assert np.abs(old[:, gcols] - delta).max() < 1e-6 and not old[:, np.setdiff1d(np.arange(G), gcols)].any()
sources = json.loads((H1_SRC_DIR / 'prediction_report.json').read_text())['source_backgrounds']
assert 'H1' not in sources
arrays = {f'v5_s{s:g}': full_array(arm(delta, s=s), gcols, G) for s in (1, 1.5, 2)}
assert np.abs(arrays['v5_s1'] - old).max() < 1e-5
print(write_dir(H1_DIR, H1_SRC_DIR, arrays, 'H1', sources))
base = dict(call='z', z=0.0, min_effect=0.0, alpha=1.0)
p1 = [dict(base, method='v5_s1', detect_floor=True, label='a1_arm1_s1_floor'),
      dict(base, method='v5_s1', detect_floor=False, label='a1_arm2_s1')]
p2 = [dict(base, method='v5_s1.5', detect_floor=False, label='a1_arm3_s1.5'),
      dict(base, method='v5_s2', detect_floor=False, label='a1_arm4_s2')]
for i, v in enumerate((p1, p2), 1): (RUN / f'h1_variants_v5a1_p{i}.json').write_text(json.dumps(v, indent=1))
print('ok', {k: float(np.abs(v).mean()) for k, v in arrays.items()})

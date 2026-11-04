"""Identity check for a build_effects_v5.py change: rebuild one existing source into a scratch folder with the current
builder (through build_effects_v5_new.py's runner) and compare every output with the table already in lib/effects_v5/:
obs, var, every layer (NaN-aware, exact), uns (except the build time), the baseline CSV and the control-quantiles npz.
Usage: ~/venv/bin/python scripts/check_v5_identity.py <dataset key> <scratch dir>     (exit 1 on any difference)
"""
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import build_effects_v5 as B                                                             # noqa: E402
import build_effects_v5_new as N                                                         # noqa: E402
import qc_replogle_umi as Q                                                              # noqa: E402

key, tmp = sys.argv[1], Path(sys.argv[2]); name = Q.NEW_CONTEXTS[key][0]; ref_dir = B.OUT
tmp.mkdir(parents=True, exist_ok=True); B.OUT = tmp
assert N.build(key, 1.2e9)
a, b = ad.read_h5ad(ref_dir / f'{name}.h5ad'), ad.read_h5ad(tmp / f'{name}.h5ad')
diff = []
if not a.obs.equals(b.obs): diff.append('obs')
if not a.var.equals(b.var): diff.append('var')
for l in a.layers:
    x, y = a.layers[l], b.layers[l]
    if hasattr(x, 'toarray'): x, y = x.toarray(), y.toarray()
    if not np.array_equal(x, y, equal_nan=x.dtype.kind == 'f'): diff.append(f'layer {l}')
ua, ub = (json.dumps({k: v for k, v in u['effects_v5'].items() if k != 'built'}, sort_keys=True, default=str) for u in (a.uns, b.uns))
if ua != ub: diff.append('uns')
if not pd.read_csv(ref_dir / f'{name}_baseline.csv').equals(pd.read_csv(tmp / f'{name}_baseline.csv')): diff.append('baseline csv')
qa, qb = np.load(ref_dir / f'{name}_control_quantiles.npz'), np.load(tmp / f'{name}_control_quantiles.npz')
if set(qa) != set(qb) or any(not np.array_equal(qa[k], qb[k]) for k in qa): diff.append('control quantiles')
print(f'IDENTITY {"PASS" if not diff else "FAIL"}: {name} rebuilt with the current builder vs lib/effects_v5/' +
      ('' if not diff else f'; differs in {diff}'), flush=True)
sys.exit(1 if diff else 0)

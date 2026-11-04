"""A1 extra arms (addendum): top-N run genes per knockdown by the DE-probability head p_de (the transformer's |z|
analogue; it has no SE), delta values, everything else 0; scored with the detection floor. Plus s = 2 + floor."""
import json
import numpy as np
from v5lib import H1_SRC_DIR, H1_DIR, RUN, full_array, panel_genes, write_dir


def topn(delta, p, n):
    out = np.zeros_like(delta); idx = np.argpartition(-p, n - 1, axis=1)[:, :n]; r = np.arange(len(delta))[:, None]
    out[r, idx] = delta[r, idx]; return out


if __name__ == '__main__':
    genes = panel_genes(); G = len(genes)
    z = np.load(H1_SRC_DIR / 'raw_predictions.npz'); delta, p, gcols = z['delta'], z['p'], z['gene_cols']
    sources = json.loads((H1_SRC_DIR / 'prediction_report.json').read_text())['source_backgrounds']
    arrays = {f'v5_top{n}': full_array(topn(delta, p, n), gcols, G) for n in (2000, 1000)}
    assert np.abs(arrays['v5_top2000'] - np.load(H1_SRC_DIR / 'tf_top2000.npy')).max() < 1e-6   # same rule as the g3 arm
    print(write_dir(H1_DIR, H1_SRC_DIR, arrays, 'H1', sources))
    base = dict(call='z', z=0.0, min_effect=0.0, alpha=1.0, detect_floor=True)
    p3 = [dict(base, method='v5_top2000', label='a1_arm5_top2000_floor'), dict(base, method='v5_s2', label='a1_arm7_s2_floor')]
    p4 = [dict(base, method='v5_top1000', label='a1_arm6_top1000_floor')]
    for i, v in ((3, p3), (4, p4)): (RUN / f'h1_variants_v5a1_p{i}.json').write_text(json.dumps(v, indent=1))
    print('ok')

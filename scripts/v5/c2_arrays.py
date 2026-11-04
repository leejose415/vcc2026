"""C2 cell-based check: per-knockdown s by within-line thirds of predicted strength (sum of p_de over run genes, no H1
data): s_lo, s_mid, s_hi = 3, 2, 3 (the pseudobulk best by mse_scaled + pds_scaled), beta 1, detection floor."""
import json
import numpy as np
from v5lib import H1_SRC_DIR, H1_DIR, RUN, arm, full_array, panel_genes, write_dir


def thirds(p):
    x = p.sum(1); return np.argsort(np.argsort(x)) * 3 // len(x)


if __name__ == '__main__':
    genes = panel_genes(); z = np.load(H1_SRC_DIR / 'raw_predictions.npz')
    sources = json.loads((H1_SRC_DIR / 'prediction_report.json').read_text())['source_backgrounds']
    s = np.array([3.0, 2.0, 3.0])[thirds(z['p'])]
    arrays = {'v5_c2_323': full_array(arm(z['delta'].astype(np.float64), s=s), z['gene_cols'], len(genes))}
    print(write_dir(H1_DIR, H1_SRC_DIR, arrays, 'H1', sources))
    v = [dict(method='v5_c2_323', call='z', z=0.0, min_effect=0.0, alpha=1.0, detect_floor=True, label='c2_s323_floor')]
    (RUN / 'h1_variants_v5c2.json').write_text(json.dumps(v, indent=1))

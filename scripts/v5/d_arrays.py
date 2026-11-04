"""D cell-based check: per-gene amplification of the specific part, s0 (1 + alpha r_g), s0 = 2, alpha = 2 (the
pseudobulk best by mse_scaled + pds_scaled); r_g = rank of training-line signal variance (lfc variance minus mean se^2,
six training lines, no H1) over the run genes, scaled to [0, 1]; beta 1, detection floor."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from v5lib import H1_SRC_DIR, H1_DIR, RUN, arm, full_array, panel_genes, write_dir


def gene_rank(gcols):
    v = np.load(Path.home() / 'h1_decomp/train_signal_var.npy')[gcols]
    rk = pd.Series(v).rank(method='average').to_numpy().copy(); rk[~np.isfinite(v)] = 1
    return (rk - rk.min()) / (rk.max() - rk.min())


if __name__ == '__main__':
    genes = panel_genes(); z = np.load(H1_SRC_DIR / 'raw_predictions.npz')
    sources = json.loads((H1_SRC_DIR / 'prediction_report.json').read_text())['source_backgrounds']
    r = gene_rank(z['gene_cols'])
    arrays = {'v5_d_s2_a2': full_array(arm(z['delta'].astype(np.float64), s=2, per_gene=1 + 2 * r), z['gene_cols'], len(genes))}
    print(write_dir(H1_DIR, H1_SRC_DIR, arrays, 'H1', sources))
    v = [dict(method='v5_d_s2_a2', call='z', z=0.0, min_effect=0.0, alpha=1.0, detect_floor=True, label='d_s2_a2_floor')]
    (RUN / 'h1_variants_v5d.json').write_text(json.dumps(v, indent=1))

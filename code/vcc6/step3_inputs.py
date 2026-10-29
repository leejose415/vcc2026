"""Step 3 inputs and checks.

  inputs   results/stage1/2026-10-27_rpe1/step3/inputs/{zero,oracle}.npz for run_line.py --lfc-npz (the 300 list; oracle =
           the replogle_rpe1 table's lfc and de_sig, n_sources = the genes it measured; zero = lfc 0, nothing called)
  profile  the controls-only profile (run_line.control_profile) against the replogle_rpe1 sidecar files: rebuilt from ALL
           RPE1 controls (the table's own control set: should match) and from the eval's context half (what run_line
           is given), max absolute differences of the quantiles and the baseline

  ~/venv/bin/python code/vcc6/step3_inputs.py inputs|profile
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
import rpe1_eval as E                                                                      # noqa: E402
import step1_oracle as O                                                                   # noqa: E402
import run_line as R                                                                       # noqa: E402

OUT = E.OUT / 'step3'


def inputs():
    d = OUT / 'inputs'; d.mkdir(parents=True, exist_ok=True)
    tg = O.targets(); genes = E.H.panel_genes()
    lfc, sig, _, _ = O.table_rows(tg, genes)
    meas = np.isfinite(lfc)
    np.savez_compressed(d / 'oracle.npz', knockdowns=np.array(tg, dtype=object), genes=np.array(genes, dtype=object),
                        lfc=lfc.astype(np.float32), de_sig=sig, n_sources=meas.astype(np.int16))
    np.savez_compressed(d / 'zero.npz', knockdowns=np.array(tg, dtype=object), genes=np.array(genes, dtype=object),
                        lfc=np.zeros_like(lfc, np.float32), de_sig=np.zeros_like(sig), n_sources=meas.astype(np.int16))
    print(f'wrote {d}/oracle.npz and zero.npz ({len(tg)} knockdowns)')


def compare(prof_dir, tag):
    import cell_generator as C
    side = C.ControlQuantiles(O.QUANT); mine = C.ControlQuantiles(prof_dir / 'controls_control_quantiles.npz')
    b0 = pd.read_csv(O.BASE, index_col=0); b1 = pd.read_csv(prof_dir / 'controls_baseline.csv', index_col=0)
    m_a, m_b = b0.control_mean_cpm.to_numpy(float), b1.control_mean_cpm.to_numpy(float)
    ok = np.isfinite(m_a)
    out = dict(source=tag, n_controls_sidecar=side.n0, n_controls_rebuilt=mine.n0,
               baseline_control_mean_cpm_max_abs_diff=float(np.nanmax(np.abs(m_a[ok] - m_b[ok]))),
               baseline_control_mean_cpm_max_rel_diff=float(np.nanmax(np.abs(m_a[ok] - m_b[ok]) / np.maximum(m_a[ok], 1e-9))),
               baseline_unmeasured_genes_rebuilt_nonzero=int((m_b[~ok] != 0).sum()),
               gated_sidecar=int(len(side.col)), gated_rebuilt=int(len(mine.col)),
               gated_symmetric_difference=int(len(set(side.col) ^ set(mine.col))))
    shared = sorted(set(side.col) & set(mine.col))
    if side.n0 == mine.n0:
        diffs = [0.0]; exact = True
        for c in shared:
            a, b = side.gene(side.index[c]), mine.gene(mine.index[c])
            same = len(a['v']) == len(b['v']) and a['nz0'] == b['nz0'] and np.array_equal(a['c'], b['c'])
            exact &= bool(same)
            if len(a['v']) == len(b['v']): diffs.append(float(np.abs(a['v'] - b['v']).max()) if len(a['v']) else 0.0)
        out.update(quantiles_identical_structure=exact, quantile_values_max_abs_diff=max(diffs))
    else:
        probs = np.linspace(0.01, 0.99, 99); diffs = []
        for c in shared:
            qa = [_q(side.gene(side.index[c]), p) for p in probs]; qb = [_q(mine.gene(mine.index[c]), p) for p in probs]
            diffs.append(float(np.max(np.abs(np.array(qa) - np.array(qb)))))
        out.update(quantile_percentile_max_abs_diff_log1p_cpm=float(np.max(diffs)),
                   quantile_percentile_median_of_gene_max=float(np.median(diffs)))
    return out


def _q(g, p):
    """The p-quantile of log1p(CPM) from one gene's stored distribution (zeros first)."""
    n = g['nz0'] + g['below'][-1]; k = p * n
    if k < g['nz0']: return 0.0
    i = int(np.searchsorted(g['below'][1:], k - g['nz0'], side='right'))
    return float(g['v'][min(i, len(g['v']) - 1)])


def profile():
    import h5py
    from anndata.io import read_elem
    import build as HB
    genes = E.H.panel_genes()
    obs, var = E.cells_obs()
    to_panel, _ = HB.panel_map(var, 'ensembl')
    ctrl_rows = np.flatnonzero(obs['is_control'].to_numpy(bool))
    X = E.H.read_rows(E.CELLS, ctrl_rows, to_panel, len(genes))
    res = []
    d_all = OUT / 'profile_all_controls'
    R.control_profile(X.tocsr(), d_all)
    res.append(compare(d_all, f'all {len(ctrl_rows)} RPE1 controls'))
    Xc = R.read_controls(E.H.eval_dir(E.NAME) / 'context.h5ad', genes)
    d_ctx = OUT / 'profile_context_half'
    R.control_profile(Xc, d_ctx)
    res.append(compare(d_ctx, f'eval context half ({Xc.shape[0]} controls)'))
    (OUT / 'profile_check.json').write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    {'inputs': inputs, 'profile': profile}[sys.argv[1]]()

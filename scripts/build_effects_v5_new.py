"""v5 effect tables for the new CRISPRi sources (KOLF2.1J, Adamson, Feng targeted, Feng genome-wide) into lib/effects_v5/,
beside the first v5 build, whose files this never rewrites. build_effects_v5.py is imported, not modified: the source is
registered for the duration of the build and its rows get the background from qc_replogle_umi.NEW_CONTEXTS (Adamson =
K562, Feng = the donor line) and a reference flag: obs/uns reference = 'non-targeting' (the marked controls) or
'rest_of_line' (Feng GW: the QC's rest-of-line reference; needs build_cells support, refused until then).
Same checks, build JSON and exit rule as build_effects_v5.py's main.

Usage: ~/venv/bin/python scripts/build_effects_v5_new.py [--sources KEY ...] [--batch-nnz 1.2e9]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import build_effects_v5 as B                                                             # noqa: E402
import qc_replogle_umi as Q                                                              # noqa: E402

KEYS = [k for k in Q.NEW_CONTEXTS if not k.startswith(('Tahoe_', 'Jiang_'))]


def build(key, batch_nnz):
    name = Q.NEW_CONTEXTS[key][0]; ref = 'rest_of_line' if Q.REFERENCE.get(key) == 'rest' else 'non-targeting'
    if ref == 'rest_of_line' and not getattr(B, 'REST_REFERENCE', False):
        sys.exit(f'{key}: rest-of-line reference not supported by build_effects_v5.build_cells yet')
    B.CELL_SOURCES[key] = name
    if (B.OUT / f'{name}.h5ad').exists(): sys.exit(f'{key}: {B.OUT / name}.h5ad exists; refusing to overwrite')
    orig_write = B.write

    def write(source, background, condition, var, obs, *rest, meta):
        obs = obs.copy(); obs['reference'] = ref
        return orig_write(source, Q.BACKGROUND[key], condition, var, obs, *rest,
                          dict(meta, reference=ref, background_from='qc_replogle_umi.NEW_CONTEXTS'))
    B.write = lambda *a: write(*a[:-1], meta=a[-1])
    try:
        t0 = time.time()
        a, ctx = B.build_cells(key, batch_nnz)
        B.log('validating (fast)')
        checks = B.validate(key, a, ctx)
    finally:
        B.write = orig_write
    checks.update(epsilon=B.EPS, epsilon_is_1e9=a.uns['effects_v5']['epsilon'] == 1e-9,
                  max_n_cells_used=int(a.obs.n_cells_used.max()), n_cells_used_le_400=bool(a.obs.n_cells_used.max() <= B.N_CELLS),
                  target_depth=B.DEPTH, thinning_factor=ctx['f'], library_check=ctx['library_check'], count_check=ctx['count_check'])
    required = dict(fast_lfc_pass=checks['fast']['lfc_pass'], fast_sig_pass=checks['fast']['sig_pass'],
                    fast_quantiles_pass=checks['fast']['quantiles_pass'], gated_sets_identical=checks['gated_sets_identical'],
                    epsilon_is_1e9=checks['epsilon_is_1e9'], n_cells_used_le_400=checks['n_cells_used_le_400'])
    failing = [k for k, v in required.items() if not v]
    lfc = a.layers['lfc']
    summary = dict(source=name, key=key, background=Q.BACKGROUND[key], reference=ref, rows=int(a.n_obs),
                   measured_genes=int(np.isfinite(lfc).any(0).sum()), de_tested_genes=a.uns['effects_v5'].get('de_tested_genes'),
                   median_n_real=float(a.obs.n_real.median()), p90_n_real=float(a.obs.n_real.quantile(0.9)),
                   median_efficiency=float(np.nanmedian(a.obs.efficiency)) if a.obs.efficiency.notna().any() else None,
                   lfc_on_epsilon_floor_share=float(np.mean(np.abs(lfc[np.isfinite(lfc)]) > 25)),
                   all_pass=not failing, required_checks=required, failing_checks=failing,
                   notes=dict(library_check=ctx['library_check']['statement'],
                              modality="obs and uns: modality = 'single_cell', significance_unknown = False",
                              runner='scripts/build_effects_v5_new.py (build_effects_v5.py imported unchanged)'),
                   checks=checks, minutes=round((time.time() - t0) / 60, 1))
    (B.OUT / f'{name}_build.json').write_text(json.dumps(summary, indent=1, default=str) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k not in ('checks', 'notes')}, indent=1, default=str), flush=True)
    return not failing


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sources', nargs='*', choices=KEYS); p.add_argument('--batch-nnz', type=float, default=1.2e9)
    args = p.parse_args()
    bad = [k for k in (args.sources or KEYS) if not build(k, args.batch_nnz)]
    if bad: print('VALIDATION FAILED:', bad, file=sys.stderr); sys.exit(1)

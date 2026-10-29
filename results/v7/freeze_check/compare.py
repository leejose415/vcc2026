"""Freeze check, part 2: v7's EvalA input (freeze_check/EvalA/input.npz) against the proxy's config 1 on H1
(cache/fold_H1.npz + cache/s11_H1.npz via scorer_step16.Mag), then v7's input scored by the same proxy.
Pass: same calls, same means to float32 precision, proxy overall within 0.001 of 0.24362."""
import json
import sys
from pathlib import Path

import numpy as np

V = Path.home() / 'vcc'; sys.path.insert(0, str(V / 'code/vcc6'))
import scorer_arms as SA                                                                   # noqa: E402
import scorer_proxy as SP                                                                  # noqa: E402
import scorer_step16 as S16                                                                # noqa: E402

OUT = V / 'results/v7/freeze_check'
G = S16.Mag('H1'); fo = G.fo
mu_px = G.build('i', S16.CONFIG1['params']['H1'])
z = np.load(OUT / 'EvalA/input.npz', allow_pickle=True)
assert list(z['knockdowns']) == list(fo.h.targets) and list(z['genes']) == list(fo.h.genes)
sig_v7 = z['de_sig'] & (z['n_sources'] > 0)
mu_v7 = SA.mu_from_lfc(z['lfc'].astype(np.float64), fo.m0)
ok = fo.ok
rel = np.abs(mu_v7 - mu_px) / np.maximum(np.abs(mu_px), 1e-3)
px = G.ref.proxy(mu=mu_px, sig=fo.sig); v7 = G.ref.proxy(mu=mu_v7, sig=sig_v7)
res = dict(
    calls_identical=bool((sig_v7 == fo.sig).all()), calls_differ=int((sig_v7 != fo.sig).sum()),
    mu_rel_diff_max_gated=float(rel[ok].max()), mu_rel_diff_p999_gated=float(np.quantile(rel[ok], 0.999)),
    called_lfc_max_abs_diff=float(np.abs(z['lfc'][fo.sig] - 0.5 * fo.lfc[fo.sig]).max()),
    proxy_overall_from_caches=px['overall'], proxy_overall_from_v7=v7['overall'],
    scaled_from_v7={SP.SHORT[m]: v7['scaled'][m] for m in SP.SIX},
    scaled_from_caches={SP.SHORT[m]: px['scaled'][m] for m in SP.SIX})
res['pass'] = res['calls_identical'] and abs(res['proxy_overall_from_v7'] - 0.24361699565446907) < 1e-3
(OUT / 'compare.json').write_text(json.dumps(res, indent=1) + '\n')
print(json.dumps(res, indent=1))

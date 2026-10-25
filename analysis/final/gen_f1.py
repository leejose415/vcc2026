"""vcc3_generate.py with the F1 'decoupled emission': unchanged generator plus, for each knockdown's CALLED genes, the
cell-to-cell spread of the normalised expression narrowed around its pooled value by a factor m (pooled value fixed):
  r_i = y_i / l_i,  r_bar = sum y / sum l,  y_i' = l_i (r_bar + m (r_i - r_bar)),  stochastic rounding, then the
  exact-totals pin restores every gene's pooled total.
m <= 1 only narrows: the spread is never wider than what the generator emits from the control cells.
Switches (benchmark variants or submit CLI): narrow_m (float, 1 = off) and called_from (name of a boolean
<name>.npy array next to the predictions, targets x genes: the genes to narrow). Everything else is vcc3_generate.py.

  python analysis/final/gen_f1.py benchmark|submit <vcc3_generate.py arguments> [--narrow-m M --called-from NAME]
"""
import sys
from pathlib import Path

import numpy as np

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc_2026_3c'))
import vcc3_generate as G  # noqa: E402

CALLED = {}            # name -> (targets -> row, bool array)
_orig_init = G.Cfg.__init__


def _init(self, a):
    _orig_init(self, a)
    self.narrow_m = float(getattr(a, 'narrow_m', 1.0) or 1.0)
    self.called_from = getattr(a, 'called_from', '') or ''


G.Cfg.__init__ = _init
G.ARM_KEYS.extend(['narrow_m', 'called_from'])


def load_called(pred_dir, name):
    key = (str(pred_dir), name)
    if key not in CALLED:
        p = Path(pred_dir); t = np.load(p / 'targets.npy').astype(str)
        CALLED[key] = ({x: i for i, x in enumerate(t)}, np.load(p / f'{name}.npy', mmap_mode='r'))
    return CALLED[key]


PRED_DIRS = {}         # line_key -> prediction folder (set by the entry point)


def narrow(B, cols, m, rng):
    if not len(cols) or m >= 1: return B
    l = B.sum(1).astype(np.float64); l = np.maximum(l, 1)
    Y = B[:, cols].astype(np.float64); rbar = Y.sum(0) / l.sum()
    target = l[:, None] * (rbar[None, :] + m * (Y / l[:, None] - rbar[None, :]))
    target = np.maximum(target, 0)
    low = np.floor(target); new = low + (rng.random(target.shape) < (target - low))
    B = B.copy(); B[:, cols] = new.astype(B.dtype)
    return B


_orig_generate_one = G.generate_one


def generate_one(line, d, j, kd, cfg, rng=None, pick=None, line_key='line', target_key='target'):
    m = getattr(cfg, 'narrow_m', 1.0); name = getattr(cfg, 'called_from', '')
    if m >= 1 or not name:
        return _orig_generate_one(line, d, j, kd, cfg, rng=rng, pick=pick, line_key=line_key, target_key=target_key)
    # run the original with exact totals off to get the cells, then narrow and pin to the original's expected totals
    import copy
    c2 = copy.copy(cfg); c2.exact_totals = True
    out, diag, planned = _orig_generate_one(line, d, j, kd, c2, rng=rng, pick=pick, line_key=line_key, target_key=target_key)
    expected = diag.pop('_expected_totals')
    pos, arr = load_called(PRED_DIRS.get(line_key, PRED_DIRS.get('*')), name)
    row = pos.get(target_key, -1)
    if row < 0: return out, diag, planned
    cols = np.flatnonzero(np.asarray(arr[row]).astype(bool) & line.gate)
    if j is not None: cols = cols[cols != j]
    B = out.toarray().astype(np.int64)
    nrng = np.random.default_rng(G.seed_for(cfg.seed, line_key, target_key, 'narrow'))
    B = narrow(B, cols, m, nrng)
    from vcc_exact_totals import pin_totals
    B = pin_totals(B, expected, nrng)
    import scipy.sparse as sp
    res = sp.csr_matrix(B.astype(np.int32)); res.eliminate_zeros(); G.validate_counts(res)
    diag['narrowed_genes'] = int(len(cols)); diag['narrow_m'] = m
    if getattr(cfg, 'exact_totals', False): diag['_expected_totals'] = expected
    return res, diag, planned


G.generate_one = generate_one


def main():
    argv = sys.argv[1:]
    # extract our switches for the submit path (benchmark takes them per variant)
    m, name = 1.0, ''
    if '--narrow-m' in argv:
        i = argv.index('--narrow-m'); m = float(argv[i + 1]); del argv[i:i + 2]
    if '--called-from' in argv:
        i = argv.index('--called-from'); name = argv[i + 1]; del argv[i:i + 2]
    sys.argv = [sys.argv[0]] + argv
    import argparse
    p = argparse.ArgumentParser(); sub = p.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('submit'); G.add_switches(s)
    for key in ['data', 'out']: s.add_argument('--' + key, required=True)
    s.add_argument('--lines', nargs='+', default=['A', 'B', 'C']); s.add_argument('--delta', action='append', required=True)
    s.add_argument('--required-n-cells', type=int, default=400)
    s.add_argument('--delta-array', default='selected'); s.add_argument('--settings')
    s.add_argument('--kd-orion', nargs='*'); s.add_argument('--kd-tcell'); s.add_argument('--kd-tcell-condition', default='Stim48hr')
    s.add_argument('--checks', type=int, default=20); s.add_argument('--workers', type=int, default=1)
    b = sub.add_parser('benchmark'); G.add_switches(b); b.set_defaults(mode='transfer')
    for key in ['cells', 'genes', 'out']: b.add_argument('--' + key, required=True)
    b.add_argument('--predictions', required=True); b.add_argument('--targets')
    b.add_argument('--obs-col', default='target_gene'); b.add_argument('--guide-col', default='')
    b.add_argument('--delta-array', nargs='+', default=['re'])
    b.add_argument('--min-cells', type=int, default=100); b.add_argument('--n-targets', type=int, default=0)
    for k, d in (('--alpha-grid', [1.0]), ('--z-grid', [0.0]), ('--min-effect-grid', [0.0]), ('--pi-grid', [0.0]), ('--up-w-grid', [0.75])):
        b.add_argument(k, nargs='+', type=float, default=d)
    b.add_argument('--cell-rules', nargs='+', default=['fixed']); b.add_argument('--cell-eval', default='cell-eval2')
    b.add_argument('--eval-profile', default='vcc2026'); b.add_argument('--num-threads', type=int, default=4)
    b.add_argument('--keep-arm-files', action='store_true'); b.add_argument('--variants'); b.add_argument('--reference')
    b.add_argument('--anchor-splits', type=int, default=5); b.add_argument('--margin-over'); b.add_argument('--scale', action='store_true')
    b.add_argument('--resume', action='store_true'); b.add_argument('--metric-directions')
    a = p.parse_args()
    a.narrow_m = m; a.called_from = name
    if a.cmd == 'benchmark':
        PRED_DIRS['*'] = a.predictions; G.cmd_benchmark(a)
    else:
        for x in a.delta:
            k, v = x.split('=', 1); PRED_DIRS[k] = v
        G.cmd_submit(a)


if __name__ == '__main__':
    main()

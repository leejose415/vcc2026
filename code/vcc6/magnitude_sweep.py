"""Magnitude sweep (Jose 2026-10-27; separate from the shrinkage study): the screen-setup C+G predictions with every
predicted LFC multiplied by m in {1, 0.5, 0.25, 0.14, 0.07, 0.03}, end to end on HCT116final and H1 Eval A (salt "").

The C+G model run is run_line's own (controls-only path; the line's tables held out):
  HCT116  results/stage1/e2e_final/CG/HCT116/model  (scored 14:48, = m 1)
  H1      results/stage1/magnitude_sweep/H1/m1/model (run_line predict on Eval A's context controls)
For m < 1, the calls are recomputed on m x lfc by the model's own rule (stage1.predictive_calls: BH on lfc /
measurement SE over the controls' 5-CPM gate, own gene excluded), so the call count moves with m as it would if the
model itself predicted smaller effects. The m = 1 recomputation must reproduce the model's calls exactly (checked).

  inputs  ~/venv/bin/python code/vcc6/magnitude_sweep.py inputs HCT116|H1   -> <sweep>/<line>/inputs/m<m>.npz
  report  ~/venv/bin/python code/vcc6/magnitude_sweep.py report             -> <sweep>/sweep.csv, figures/sweep.png
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
import stage1 as S                                                                         # noqa: E402

OUT = V / 'results/stage1/magnitude_sweep'
MS = [1, 0.5, 0.25, 0.14, 0.07, 0.03]
MODEL = {'HCT116': V / 'results/stage1/e2e_final/CG/HCT116', 'H1': OUT / 'H1/m1'}


def mtag(m):
    return f'm{m:g}'


def inputs(line):
    d = MODEL[line]
    z = np.load(d / 'model/predictions.npz', allow_pickle=True)
    base = pd.read_csv(d / 'controls/controls_baseline.csv', index_col=0)
    genes = [str(g) for g in z['genes']]
    gate = base['de_tested'].reindex(genes).fillna(False).to_numpy(bool)
    col = {g: j for j, g in enumerate(genes)}
    lfc = z['lfc'].astype(np.float64); se = z['se_measurement'].astype(np.float64)
    out = OUT / line / 'inputs'; out.mkdir(parents=True, exist_ok=True)
    for m in MS:
        sig = np.zeros(lfc.shape, bool)
        for i, k in enumerate(z['knockdowns']):
            t = gate.copy()
            if str(k) in col: t[col[str(k)]] = False
            sig[i], _ = S.predictive_calls(m * lfc[i], se[i], t)
        if m == 1:
            same = (sig == z['de_sig']).all()
            print(f'{line}: m 1 recomputed calls reproduce the model\'s: {same}')
            assert same, 'call recomputation does not reproduce the model'
        np.savez_compressed(out / f'{mtag(m)}.npz', knockdowns=z['knockdowns'], genes=z['genes'],
                            lfc=(m * lfc).astype(np.float32), de_sig=sig, n_sources=z['n_sources'])
        print(f'{line} m {m:g}: calls per knockdown mean {sig.sum(1).mean():.1f}, median {np.median(sig.sum(1)):.0f}')


def report():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = []
    for line in MODEL:
        for m in MS:
            d = (MODEL[line] if (line == 'HCT116' and m == 1) else OUT / line / mtag(m)) / 'salt_none'
            ev = 'evalHCT116final' if line == 'HCT116' else 'evalA'
            sp = d / 'harness' / ev / 'summary.json'
            if not sp.exists():
                continue
            s = json.loads(sp.read_text())
            npz = (OUT / line / 'inputs' / f'{mtag(m)}.npz')
            calls = float(np.load(npz, allow_pickle=True)['de_sig'].sum(1).mean()) if npz.exists() else np.nan
            rows.append(dict(line=line, m=m, calls_mean=calls, overall=s['scaled']['overall'],
                             **{f'scaled:{k}': v for k, v in s['scaled'].items() if k != 'overall'},
                             **{f'raw:{k}': v for k, v in s['raw'].items()}))
    d = pd.DataFrame(rows); d.to_csv(OUT / 'sweep.csv', index=False)
    comps = ['overall'] + sorted(c for c in d.columns if c.startswith('scaled:')) + sorted(c for c in d.columns if c.startswith('raw:')) + ['calls_mean']
    n = len(comps); cols = 4; rws = int(np.ceil(n / cols))
    fig, axs = plt.subplots(rws, cols, figsize=(4.2 * cols, 3.2 * rws))
    col = {'HCT116': '#eb6834', 'H1': '#2a78d6'}
    for ax, c in zip(axs.ravel(), comps):
        for line, g in d.groupby('line'):
            g = g.sort_values('m')
            ax.plot(g.m, g[c], 'o-', color=col[line], lw=2, label=line)
            if c == 'overall':
                best = g.loc[g.overall.idxmax()]
                ax.scatter([best.m], [best.overall], s=160, facecolor='none', edgecolor=col[line], lw=2)
        ax.set_xscale('log'); ax.set_title(c, loc='left', fontsize=9); ax.set_xlabel('m')
        ax.spines[['top', 'right']].set_visible(False); ax.grid(color='#e6e5e0', lw=0.7)
    for ax in axs.ravel()[n:]:
        ax.axis('off')
    axs.ravel()[0].legend(frameon=False, fontsize=8)
    fig.suptitle('Magnitude sweep: C+G LFC × m, end to end (circle: best scaled overall per line)', x=0.01, ha='left')
    fig.tight_layout(); (OUT / 'figures').mkdir(exist_ok=True); fig.savefig(OUT / 'figures/sweep.png', dpi=140)
    print(d[['line', 'm', 'calls_mean', 'overall']].to_string())


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['inputs', 'report'])
    p.add_argument('line', nargs='?')
    a = p.parse_args()
    inputs(a.line) if a.cmd == 'inputs' else report()

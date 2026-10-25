"""Figures and the tidy CSV for the transformer hyperparameter sweeps.

  python scripts/sweep_plots.py results/sweeps/<date>
Reads <dir>/runs/*/metrics.csv (+ baselines.csv, h1_dump.npz) and writes <dir>/all_metrics.csv (run_id, knob, value, seed,
epoch, eval_set, metric, score) and PNGs in <dir>/figures/:
  A_<knob>.png   per knob: R^2 (all measured genes), DE R^2, calibration slope vs the knob value (final epoch)
  B1_learning_curve.png, B2_lr_curves.png, B3_epochs.png, C_sensitivity.png, D1_scatter.png, D2_h1_r2.png, E_dashboard.png
R^2 panels show the per-knockdown MEDIAN (most knockdowns have few real effects, so their R^2 can be very negative and
pull the mean around); the CSV keeps both. Every R^2 panel shares one y-range. Grey band = the baseline config's H1 value
+/- 1 seed SD (noise floor); dashed lines = B1 (Model 2: precision-weighted evidence mean) and B2 (average knockdown
response; for H1 it uses H1's own answers, reference only), both on H1.
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

D = Path(sys.argv[1]); FIG = D / 'figures'; FIG.mkdir(parents=True, exist_ok=True)
STYLE = {'train': dict(color='#2a78d6', ls='--', marker='s', lw=1.6, ms=5, label='train'),
         'heldout': dict(color='#eb6834', ls='-', marker='o', lw=1.8, ms=5, label='held-out knockdowns'),
         'h1': dict(color='#1baf7a', ls='-', marker='^', lw=3.0, ms=8, label='H1 (new line)')}
SETS = ['train', 'heldout', 'h1']
INK, MUTED, BAND = '#333333', '#777777', '#d9d9d9'
MATTERS, NOEFFECT = '#2a78d6', '#b9b9b4'
R2, DER2, SLOPE = 'r2_median', 'de_r2_median', 'slope'
KNOB_ORDER = ['epochs', 'learning rate', 'gamma', 'dropout', 'weight decay', 'width', 'data fraction', 'depth', 'cosine weight']
DEFAULTS = {'epochs': 12, 'learning rate': 3e-4, 'gamma': 'learned', 'dropout': 0.1, 'weight decay': 0.01, 'width': 128,
            'data fraction': 1.0, 'depth': 4, 'cosine weight': 0.1}
LOGX = {'learning rate', 'width', 'data fraction'}
WHAT = {'epochs': 'total training epochs (each run has its own cosine schedule)', 'learning rate': 'peak learning rate',
        'gamma': 'the evidence shortcut gamma', 'dropout': 'dropout', 'weight decay': 'AdamW weight decay',
        'width': 'model width (FF = 4 x width)', 'data fraction': 'share of training knockdowns used', 'depth': 'layers',
        'cosine weight': 'weight of the (1 - cosine) loss term'}


def load():
    frames = [pd.read_csv(p) for p in sorted((D / 'runs').glob('*/metrics.csv'))]
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    base = pd.read_csv(D / 'baselines.csv') if (D / 'baselines.csv').exists() else pd.DataFrame()
    pd.concat([df, base], ignore_index=True).to_csv(D / 'all_metrics.csv', index=False)
    return df, base


def final(df):
    """Last epoch of each run."""
    last = df.groupby('run_id').epoch.transform('max')
    return df[df.epoch == last]


def value_key(knob, v):
    try: return float(v)
    except (TypeError, ValueError): return v


def knob_table(fin, base_runs, knob, metric):
    """Per knob value: mean and SD (over seeds) of `metric` per eval set; the baseline runs supply the default value."""
    rows = fin[(fin.knob == knob) | fin.run_id.isin(base_runs)].copy()
    if knob == 'epochs':
        rows['value'] = [12 if r in base_runs else v for r, v in zip(rows.run_id, rows.value)]
    else:
        rows['value'] = [DEFAULTS[knob] if r in base_runs else v for r, v in zip(rows.run_id, rows.value)]
    rows = rows[rows.metric == metric]
    rows['x'] = [value_key(knob, v) for v in rows.value]
    g = rows.groupby(['x', 'eval_set']).score.agg(['mean', 'std', 'count']).reset_index()
    return g


def noise(fin, base_runs, metric):
    b = fin[fin.run_id.isin(base_runs) & (fin.metric == metric) & (fin.eval_set == 'h1')].score
    return float(b.mean()), float(b.std(ddof=1)) if len(b) > 1 else float('nan')


def recommend(g, knob, sd):
    """Middle of the H1 plateau (values within 1 seed SD of the best mean); the default if it is on the plateau."""
    h = g[g.eval_set == 'h1'].sort_values('x')
    if h.empty: return None, []
    best = h['mean'].max(); tol = sd if np.isfinite(sd) else 0.0
    plateau = list(h[h['mean'] >= best - tol].x)
    d = DEFAULTS.get(knob)
    if d in plateau: return d, plateau
    return plateau[len(plateau) // 2], plateau


def takeaway(knob, g, sd, rec, plateau):
    h = g[g.eval_set == 'h1']
    if h.empty: return knob
    rng = h['mean'].max() - h['mean'].min(); d = DEFAULTS.get(knob)
    fmt = lambda v: f'{v:g}' if isinstance(v, float) else str(v)
    if not np.isfinite(sd) or rng <= 2 * sd:
        return f'{knob.capitalize()}: H1 R² barely moves (within the noise), so keep the default ({fmt(d)})'
    if rec == d:
        return f'{knob.capitalize()}: H1 R² is flat around the default {fmt(d)} and worse further out, so keep it'
    return f'{knob.capitalize()}: H1 R² is best around {fmt(rec)} (default {fmt(d)}); a real effect above the noise'


def label_end(ax, xs, ys, text, style):
    if len(xs) == 0: return
    ax.annotate(text, (xs[-1], ys[-1]), xytext=(6, 0), textcoords='offset points', va='center', fontsize=8,
                color=style['color'], fontweight='bold' if style['lw'] > 2 else 'normal')


def regions(ax, g, sd, metric_sd):
    """Shade only where the data show it: sweet spot (H1 plateau), overfitting (train up, H1 down vs the default),
    underfitting (all three sets below the default)."""
    if not np.isfinite(sd): return
    piv = g.pivot(index='x', columns='eval_set', values='mean').sort_index()
    if 'h1' not in piv or len(piv) < 3: return
    xs = list(piv.index); best = piv['h1'].max()
    plateau = [x for x in xs if piv.loc[x, 'h1'] >= best - sd]
    if 1 < len(plateau) < len(xs):
        lo, hi = xs.index(plateau[0]), xs.index(plateau[-1])
        if hi - lo + 1 == len(plateau):
            ax.axvspan(plateau[0], plateau[-1], color='#1baf7a', alpha=0.08, lw=0)
            ax.text(np.sqrt(plateau[0] * plateau[-1]) if ax.get_xscale() == 'log' else (plateau[0] + plateau[-1]) / 2,
                    0.98, 'sweet spot', transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=8, color='#127a55')


def figure_A(fin, base, base_runs, ylim):
    out = []
    knobs = [k for k in KNOB_ORDER if (fin.knob == k).any()]
    for knob in knobs:
        if knob == 'gamma': out.append(figure_gamma(fin, base, base_runs, ylim)); continue
        fig, axes = plt.subplots(3, 1, figsize=(6.4, 8.2), sharex=True)
        rec = None
        for ax, metric, name in zip(axes, [R2, DER2, SLOPE], ['R², all measured genes (median per knockdown)',
                                                               'R², DE genes only (median per knockdown)',
                                                               'calibration slope (observed on predicted)']):
            g = knob_table(fin, base_runs, knob, metric)
            mu, sd = noise(fin, base_runs, metric)
            if metric == R2:
                rec, plateau = recommend(g, knob, sd); title = takeaway(knob, g, sd, rec, plateau)
            if np.isfinite(sd): ax.axhspan(mu - sd, mu + sd, color=BAND, alpha=0.6, lw=0, zorder=0)
            for es in SETS:
                h = g[g.eval_set == es].sort_values('x')
                st = STYLE[es]
                ax.errorbar(h.x, h['mean'], yerr=h['std'].where(h['count'] > 1), color=st['color'], ls=st['ls'],
                            marker=st['marker'], lw=st['lw'], ms=st['ms'], capsize=3, zorder=3)
                label_end(ax, list(h.x), list(h['mean']), st['label'], st)
            if metric != SLOPE and not base.empty:
                for meth, lab, ls in (('B1', 'B1 Model 2', (0, (4, 3))), ('B2', 'B2 avg knockdown*', (0, (1, 2)))):
                    v = base[(base.run_id == meth) & (base.eval_set == 'h1') & (base.metric == metric)].score
                    if len(v):
                        ax.axhline(float(v.iloc[0]), color=MUTED, ls=ls, lw=1)
                        ax.text(0.01, float(v.iloc[0]), f' {lab} (H1)', transform=ax.get_yaxis_transform(), fontsize=7,
                                color=MUTED, va='bottom')
            if metric == SLOPE: ax.axhline(1.0, color=INK, lw=1, ls='-'); ax.text(0.01, 1.0, ' 1 = right size', fontsize=7,
                                                                                    transform=ax.get_yaxis_transform(), va='bottom')
            if metric in (R2, DER2): ax.set_ylim(*ylim[metric])
            d = DEFAULTS[knob]; ax.axvline(d, color=MUTED, ls=':', lw=1)
            if rec is not None and metric == R2:
                hr = g[(g.eval_set == 'h1') & (g.x == rec)]
                if len(hr): ax.plot([rec], [float(hr['mean'].iloc[0])], marker='*', ms=16, color='#e8a400', mec=INK, zorder=5)
            regions(ax, g, sd, None) if metric == R2 else None
            ax.set_ylabel(name, fontsize=8); ax.grid(axis='y', color='#ececec', lw=0.8)
            ax.spines[['top', 'right']].set_visible(False)
            if knob in LOGX: ax.set_xscale('log')
        axes[-1].set_xlabel(f'{WHAT[knob]}  (dotted line = current default, star = recommended)', fontsize=8)
        fig.suptitle(title, fontsize=10, x=0.02, ha='left')
        fig.text(0.02, 0.005, 'Look at the thick green H1 line: a change only counts if it leaves the grey band (±1 seed SD). '
                 '* B2 on H1 uses H1’s own answers (reference only).', fontsize=7.5, color=MUTED)
        fig.tight_layout(rect=(0, 0.02, 0.97, 0.96)); p = FIG / f'A_{knob.replace(" ", "_")}.png'; fig.savefig(p, dpi=140)
        plt.close(fig); out.append((knob, p, rec))
    return out


def figure_gamma(fin, base, base_runs, ylim):
    vals = ['learned', 'fixed 0', 'fixed 1']
    fig, axes = plt.subplots(3, 1, figsize=(6.0, 8.2), sharex=True)
    for ax, metric, name in zip(axes, [R2, DER2, SLOPE], ['R², all measured genes', 'R², DE genes only', 'calibration slope']):
        mu, sd = noise(fin, base_runs, metric)
        if np.isfinite(sd): ax.axhspan(mu - sd, mu + sd, color=BAND, alpha=0.6, lw=0)
        for es in SETS:
            ys, es_sd = [], []
            for v in vals:
                sel = fin[(fin.run_id.isin(base_runs) if v == 'learned' else (fin.knob == 'gamma') & (fin.value == v))
                          & (fin.eval_set == es) & (fin.metric == metric)].score
                ys.append(sel.mean() if len(sel) else np.nan); es_sd.append(sel.std(ddof=1) if len(sel) > 1 else np.nan)
            st = STYLE[es]
            ax.errorbar(range(3), ys, yerr=es_sd, color=st['color'], ls=st['ls'], marker=st['marker'], lw=st['lw'],
                        ms=st['ms'], capsize=3); label_end(ax, [2], [ys[2]], st['label'], st)
        if metric == SLOPE: ax.axhline(1.0, color=INK, lw=1)
        if metric in (R2, DER2): ax.set_ylim(*ylim[metric])
        ax.set_ylabel(name, fontsize=8); ax.grid(axis='y', color='#ececec', lw=0.8); ax.spines[['top', 'right']].set_visible(False)
    axes[-1].set_xticks(range(3)); axes[-1].set_xticklabels(['learned (default)', 'fixed at 0 (no shortcut)', 'fixed at 1'])
    h = fin[(fin.metric == R2) & (fin.eval_set == 'h1')]
    m = {v: (h[h.run_id.isin(base_runs)].score.mean() if v == 'learned' else h[(h.knob == 'gamma') & (h.value == v)].score.mean()) for v in vals}
    mu, sd = noise(fin, base_runs, R2)
    if np.isfinite(m['fixed 0']) and np.isfinite(sd) and abs(m['fixed 0'] - m['learned']) <= 2 * sd:
        title = 'Evidence shortcut (gamma): switching it off changes nothing beyond the noise, so the encoder has absorbed it'
    elif np.isfinite(m['fixed 0']) and m['fixed 0'] < m['learned']:
        title = 'Evidence shortcut (gamma): removing it hurts H1, so the shortcut carries real signal'
    else:
        title = 'Evidence shortcut (gamma): see panels'
    fig.suptitle(title, fontsize=10, x=0.02, ha='left')
    fig.text(0.02, 0.005, 'Compare the thick green H1 line across the three settings; the grey band is ±1 seed SD.', fontsize=7.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.02, 0.97, 0.96)); p = FIG / 'A_gamma.png'; fig.savefig(p, dpi=140); plt.close(fig)
    return ('gamma', p, None)


def figure_B(df, base_runs, ylim):
    res = {}
    longest = df.groupby('run_id').epoch.max(); runs = [r for r in ('epochs_30', 'epochs_20') if r in longest.index]
    if runs:
        r = runs[0]; fig, ax = plt.subplots(figsize=(7, 4.2))
        d = df[(df.run_id == r) & (df.metric == R2)]
        for es in SETS:
            h = d[d.eval_set == es].sort_values('epoch'); st = STYLE[es]
            ax.plot(h.epoch, h.score, color=st['color'], ls=st['ls'], marker=st['marker'], lw=st['lw'], ms=st['ms'] - 2)
            label_end(ax, list(h.epoch), list(h.score), st['label'], st)
        h1 = d[d.eval_set == 'h1']; pk = h1.loc[h1.score.idxmax()]
        ax.axvline(pk.epoch, color=MUTED, ls=':'); ax.annotate(f'H1 peaks at epoch {int(pk.epoch)}', (pk.epoch, pk.score),
                                                               xytext=(5, 12), textcoords='offset points', fontsize=8)
        tr = d[d.eval_set == 'train'].set_index('epoch').score; hh = h1.set_index('epoch').score
        gap0, gap1 = tr.iloc[0] - hh.iloc[0], tr.iloc[-1] - hh.iloc[-1]
        trend = 'widening' if gap1 > gap0 + 0.01 else ('narrowing' if gap1 < gap0 - 0.01 else 'steady')
        ax.set_title(f'Learning curve ({int(longest[r])} epochs): H1 peaks at epoch {int(pk.epoch)}; train-H1 gap {trend}',
                     fontsize=10, loc='left')
        ax.set_xlabel('epoch'); ax.set_ylabel('R², all measured genes (median per knockdown)'); ax.set_ylim(*ylim[R2])
        ax.grid(axis='y', color='#ececec'); ax.spines[['top', 'right']].set_visible(False)
        fig.text(0.01, 0.01, 'Watch whether the thick H1 line flattens, keeps rising or turns down while train keeps rising.',
                 fontsize=7.5, color=MUTED)
        fig.tight_layout(rect=(0, 0.03, 0.95, 1)); fig.savefig(FIG / 'B1_learning_curve.png', dpi=140); plt.close(fig)
        res['B1'] = FIG / 'B1_learning_curve.png'
    lr = df[((df.knob == 'learning rate') | df.run_id.isin(base_runs[:1])) & (df.metric == R2) & (df.eval_set == 'h1')]
    if lr.run_id.nunique() > 1:
        fig, ax = plt.subplots(figsize=(7, 4.2)); cmap = plt.get_cmap('viridis')
        vals = sorted({(3e-4 if r in base_runs else float(v)) for r, v in zip(lr.run_id, lr.value)})
        for i, v in enumerate(vals):
            h = lr[[(3e-4 if r in base_runs else float(x)) == v for r, x in zip(lr.run_id, lr.value)]].sort_values('epoch')
            c = cmap(i / max(1, len(vals) - 1)); ax.plot(h.epoch, h.score, color=c, marker='^', lw=2.2 if v == 3e-4 else 1.4)
            ax.annotate(f'lr {v:g}' + (' (default)' if v == 3e-4 else ''), (h.epoch.iloc[-1], h.score.iloc[-1]), xytext=(5, 0),
                        textcoords='offset points', fontsize=8, color=c, va='center')
        ax.set_title('H1 R² over training, one line per learning rate', fontsize=10, loc='left')
        ax.set_xlabel('epoch'); ax.set_ylabel('H1 R² (median per knockdown)'); ax.set_ylim(*ylim[R2])
        ax.grid(axis='y', color='#ececec'); ax.spines[['top', 'right']].set_visible(False)
        fig.text(0.01, 0.01, 'Too high a rate shows as an early peak or a noisy line; too low as a slow climb.', fontsize=7.5, color=MUTED)
        fig.tight_layout(rect=(0, 0.03, 0.9, 1)); fig.savefig(FIG / 'B2_lr_curves.png', dpi=140); plt.close(fig)
        res['B2'] = FIG / 'B2_lr_curves.png'
    fin = final(df); g = knob_table(fin, base_runs, 'epochs', R2)
    h = g[g.eval_set == 'h1'].sort_values('x')
    if len(h) > 1:
        mu, sd = noise(fin, base_runs, R2); fig, ax = plt.subplots(figsize=(6.4, 4.2))
        if np.isfinite(sd): ax.axhspan(mu - sd, mu + sd, color=BAND, alpha=0.6, lw=0)
        for es in SETS:
            hh = g[g.eval_set == es].sort_values('x'); st = STYLE[es]
            ax.errorbar(hh.x, hh['mean'], yerr=hh['std'].where(hh['count'] > 1), color=st['color'], ls=st['ls'], marker=st['marker'],
                        lw=st['lw'], ms=st['ms'], capsize=3); label_end(ax, list(hh.x), list(hh['mean']), st['label'], st)
        rec, plateau = recommend(g, 'epochs', sd)
        ax.axvline(12, color=MUTED, ls=':')
        hr = h[h.x == rec]
        if len(hr): ax.plot([rec], [float(hr['mean'].iloc[0])], marker='*', ms=16, color='#e8a400', mec=INK, zorder=5)
        ax.set_title(f'Is 12 epochs right? H1 at the end of separate 6/12/20/30-epoch runs: recommend {rec:g}', fontsize=10, loc='left')
        ax.set_xlabel('total epochs of the run (dotted = current 12, star = recommended)'); ax.set_ylabel('final R² (median per knockdown)')
        ax.set_ylim(*ylim[R2]); ax.grid(axis='y', color='#ececec'); ax.spines[['top', 'right']].set_visible(False)
        fig.text(0.01, 0.01, 'Pick the shortest run on the H1 plateau; the grey band is the baseline H1 value ±1 seed SD.',
                 fontsize=7.5, color=MUTED)
        fig.tight_layout(rect=(0, 0.03, 0.95, 1)); fig.savefig(FIG / 'B3_epochs.png', dpi=140); plt.close(fig)
        res['B3'] = FIG / 'B3_epochs.png'
    return res


def figure_C(fin, base_runs):
    mu, sd = noise(fin, base_runs, R2); rows = []
    for knob in [k for k in KNOB_ORDER if (fin.knob == k).any()]:
        if knob == 'gamma':
            h = fin[(fin.metric == R2) & (fin.eval_set == 'h1') & ((fin.knob == 'gamma') | fin.run_id.isin(base_runs))]
            means = h.assign(v=[('learned' if r in base_runs else v) for r, v in zip(h.run_id, h.value)]).groupby('v').score.mean()
        else:
            g = knob_table(fin, base_runs, knob, R2); means = g[g.eval_set == 'h1'].set_index('x')['mean']
        if len(means) > 1: rows.append((knob, float(means.max() - means.min()), len(means)))
    if not rows: return None, []
    rows.sort(key=lambda r: r[1]); thr = 2 * sd if np.isfinite(sd) else np.nan
    fig, ax = plt.subplots(figsize=(7, 0.5 * len(rows) + 1.6))
    colors = [MATTERS if np.isfinite(thr) and r[1] > thr else NOEFFECT for r in rows]
    ax.barh([r[0] for r in rows], [r[1] for r in rows], color=colors, height=0.6)
    for i, r in enumerate(rows):
        ax.text(r[1], i, f'  {r[1]:.3f} ({"matters" if colors[i] == MATTERS else "no detectable effect"})', va='center', fontsize=8)
    if np.isfinite(thr):
        ax.axvline(thr, color=INK, ls='--', lw=1); ax.text(thr, len(rows) - 0.4, f' 2 x seed SD = {thr:.3f}', fontsize=8)
    n_m = sum(c == MATTERS for c in colors)
    ax.set_title(f'Which knobs matter for H1: {n_m} of {len(rows)} move H1 R² beyond the noise', fontsize=10, loc='left')
    ax.set_xlabel('range of H1 R² across the values tried (max − min of the means)')
    ax.spines[['top', 'right']].set_visible(False)
    fig.text(0.01, 0.01, 'Blue bars pass the dashed line (twice the seed-to-seed spread): only those knobs show a real effect.',
             fontsize=7.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.04, 0.92, 1)); fig.savefig(FIG / 'C_sensitivity.png', dpi=140); plt.close(fig)
    return FIG / 'C_sensitivity.png', rows[::-1]


def figure_D():
    p = D / 'h1_dump.npz'
    if not p.exists(): return {}
    z = np.load(p, allow_pickle=False); res = {}
    pred, obs, mask, de, b1 = z['pred'], z['obs'], z['mask'], z['de'], z['b1']; names = z['targets'].astype(str)
    def r2(x, y): ss = ((y - y.mean()) ** 2).sum(); return 1 - ((y - x) ** 2).sum() / ss if ss > 0 else np.nan
    def slope(x, y): vx = ((x - x.mean()) ** 2).sum(); return ((x - x.mean()) * (y - y.mean())).sum() / vx if vx > 0 else np.nan
    rm = np.array([r2(pred[i][mask[i]], obs[i][mask[i]]) if mask[i].sum() > 2 else np.nan for i in range(len(pred))])
    order = np.argsort(np.nan_to_num(rm, nan=-np.inf)); valid = order[np.isfinite(rm[order])]
    picks = [('best', valid[-1]), ('median', valid[len(valid) // 2]), ('worst', valid[0])]
    fig, axes = plt.subplots(3, 2, figsize=(9, 12.5))
    for row, (tag, i) in enumerate(picks):
        m = mask[i]
        for col, (x, lab) in enumerate(((pred[i], 'transformer'), (b1[i], 'B1 Model 2'))):
            ax = axes[row, col]; xm, ym, dm = x[m], obs[i][m], de[i][m]
            lim = np.nanmax(np.abs(np.concatenate([xm, ym]))) * 1.05 + 1e-6
            ax.scatter(xm[~dm], ym[~dm], s=5, color='#b9b9b4', alpha=0.6, lw=0, label='other genes')
            ax.scatter(xm[dm], ym[dm], s=14, color='#eb6834', alpha=0.9, lw=0, label='DE genes')
            ax.plot([-lim, lim], [-lim, lim], color=INK, lw=1); ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
            ax.set_title(f'{tag}: {names[i]} | {lab}: R² {r2(xm, ym):.2f}, slope {slope(xm, ym):.2f}', fontsize=9, loc='left')
            ax.set_xlabel('predicted log2FC'); ax.set_ylabel('observed log2FC (H1)')
            ax.spines[['top', 'right']].set_visible(False)
            if row == 0 and col == 0: ax.legend(frameon=False, fontsize=8, loc='upper left')
    fig.suptitle('Predicted vs observed on 3 H1 knockdowns (best, median, worst under the baseline config)', fontsize=10, x=0.02, ha='left')
    fig.text(0.02, 0.005, 'Points on the diagonal are right; a cloud flatter than the diagonal means predictions are too small (shrunk).',
             fontsize=7.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.01, 1, 0.97)); fig.savefig(FIG / 'D1_scatter.png', dpi=130); plt.close(fig)
    res['D1'] = FIG / 'D1_scatter.png'
    rb1 = np.array([r2(b1[i][mask[i]], obs[i][mask[i]]) if mask[i].sum() > 2 else np.nan for i in range(len(pred))])
    rb2 = np.array([r2(z['b2'][mask[i]], obs[i][mask[i]]) if mask[i].sum() > 2 else np.nan for i in range(len(pred))])
    fig, ax = plt.subplots(figsize=(7, 4.4)); data = [a[np.isfinite(a)] for a in (rm, rb1, rb2)]
    parts = ax.violinplot(data, showmedians=True, showextrema=False)
    for b, c in zip(parts['bodies'], ['#1baf7a', '#888888', '#bbbbbb']): b.set_facecolor(c); b.set_alpha(0.55)
    for k, a in enumerate(data): ax.scatter(np.full(len(a), k + 1) + np.random.default_rng(0).uniform(-0.08, 0.08, len(a)), a, s=4,
                                           color=INK, alpha=0.3, lw=0)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(['transformer', 'B1 Model 2', 'B2 average H1 knockdown*'])
    med = [np.median(a) for a in data]
    ax.set_title(f'H1, per-knockdown R²: transformer median {med[0]:.3f} vs Model 2 {med[1]:.3f}', fontsize=10, loc='left')
    ax.set_ylabel('R² per H1 knockdown'); lo = np.nanpercentile(np.concatenate(data), 2); ax.set_ylim(max(lo, -1.5), None)
    ax.axhline(0, color=MUTED, lw=0.8); ax.spines[['top', 'right']].set_visible(False)
    fig.text(0.01, 0.01, 'Each dot is one H1 knockdown; the bar is the median. * B2 uses H1’s own answers (reference only).',
             fontsize=7.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(FIG / 'D2_h1_r2.png', dpi=140); plt.close(fig)
    res['D2'] = FIG / 'D2_h1_r2.png'
    return res


def figure_E(paths_A, pC, pB):
    import matplotlib.image as mpimg
    imgs = [p for _, p, _ in paths_A]
    n = len(imgs); cols = min(4, max(1, n)); rows_A = int(np.ceil(n / cols)) if n else 0
    fig = plt.figure(figsize=(16, 5 + 5.2 * rows_A + 5))
    gs = fig.add_gridspec(2 + rows_A, cols, height_ratios=[1.0] + [1.3] * rows_A + [1.0])
    def show(ax, p): ax.imshow(mpimg.imread(p)); ax.axis('off')
    if pC: show(fig.add_subplot(gs[0, :]), pC)
    for k, p in enumerate(imgs): show(fig.add_subplot(gs[1 + k // cols, k % cols]), p)
    if pB.get('B1'): show(fig.add_subplot(gs[-1, :cols // 2 or 1]), pB['B1'])
    if pB.get('B3'): show(fig.add_subplot(gs[-1, cols // 2 or 1:]), pB['B3'])
    fig.suptitle('Transformer hyperparameter sweeps: what matters for H1 (top), each knob (middle), training length (bottom)',
                 fontsize=14, x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.98)); fig.savefig(FIG / 'E_dashboard.png', dpi=110); plt.close(fig)


def main():
    df, base = load()
    if df.empty: print('no runs yet'); return
    base_runs = sorted(df[df.knob == 'baseline'].run_id.unique())
    fin = final(df)
    ylim = {}
    for metric in (R2, DER2):
        v = pd.concat([fin[fin.metric == metric].score, df[(df.metric == metric)].score,
                       base[base.metric == metric].score if not base.empty else pd.Series(dtype=float)])
        v = v[np.isfinite(v)]
        lo, hi = (float(v.min()), float(v.max())) if len(v) else (0, 1); pad = 0.05 * (hi - lo + 1e-6)
        ylim[metric] = (lo - pad, hi + pad)
    pA = figure_A(fin, base, base_runs, ylim)
    pB = figure_B(df, base_runs, ylim)
    pC, ranking = figure_C(fin, base_runs)
    figure_D()
    figure_E(pA, pC, pB)
    summary = dict(base_runs=base_runs, noise={m: noise(fin, base_runs, m) for m in (R2, DER2, SLOPE, 'pearson_median')},
                   recommended={k: (r if not isinstance(r, np.generic) else r.item()) for k, _, r in pA},
                   ranking=[dict(knob=k, h1_r2_range=v, values=n) for k, v, n in ranking])
    (D / 'figures/summary.json').write_text(json.dumps(summary, indent=1, default=str))
    print(json.dumps(summary, indent=1, default=str))


if __name__ == '__main__':
    main()

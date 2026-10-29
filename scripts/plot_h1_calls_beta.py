"""H1 (Eval A, line held out) figures for the incumbent C+G: which genes are called, and true vs predicted Delta.

Input: the npz from scripts/h1_calls_beta_prep.py <out_dir> (truth / realised calls from the C+G end-to-end harness run, p_adj < 0.05, own gene out;
predicted Delta = the C+G mean-inversion proxy in log1p CP50K, true Delta = held-out H1 pseudobulk minus control).

  python scripts/plot_h1_calls_beta.py <h1_arrays.npz> <out_dir>
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                                           # noqa: E402
import numpy as np                                                                       # noqa: E402

HIT, TRUE, CALL, NONE = '#2a78d6', '#eb6834', '#1baf7a', '#c9c8c2'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e6e5e0'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'figure.facecolor': '#fcfcfb',
                     'axes.facecolor': '#fcfcfb'})

z = np.load(sys.argv[1], allow_pickle=True); out = Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
T, R, ok, tgt = z['truth'], z['real'], z['ok'], list(z['targets'])
m0, gate, lt = z['m0'], z['gate'], z['lfc_true']
pdel, tdel = z['pred_delta'].astype(np.float64), z['true_delta'].astype(np.float64)
T, R = T & ok, R & ok
EX = 'SIN3B'

# ---------- figure 1: calls ----------
fig = plt.figure(figsize=(15, 10.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.38, wspace=0.22)
ax = fig.add_subplot(gs[0, :])
nt, nr, nh = T.sum(1), R.sum(1), (T & R).sum(1)
o = np.argsort(nt, kind='stable'); x = np.arange(len(o))
for v, c, lab in ((nt, TRUE, 'true DE genes (real H1 cells)'), (nr, CALL, "genes the model's cells call DE"),
                  (nh, HIT, 'called and truly DE (hits)')):
    ax.plot(x, v[o], 'o', ms=3.2, color=c, label=lab, mec='white', mew=0.3)
ax.set_yscale('symlog', linthresh=1); ax.set_ylim(-0.3, 1.2e4)
ax.set_xlim(-3, len(o) + 2); ax.set_xlabel('H1 knockdowns, sorted by number of true DE genes')
ax.set_ylabel('genes per knockdown')
ax.set_title(f'A. Per knockdown: truth vs calls. Median true DE {np.median(nt):.0f}, median model calls {np.median(nr):.0f}; '
             f'pooled recall {nh.sum() / nt.sum():.0%}, precision {nh.sum() / nr.sum():.0%}', loc='left', color=INK)
ax.legend(loc='lower right', frameon=False)
i = tgt.index(EX); j = int(np.flatnonzero(o == i)[0])
ax.annotate(EX, (j, nt[i]), (j - 40, nt[i] * 2.2), color=INK2, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.6))

ax = fig.add_subplot(gs[1, 0])
lm = np.log10(m0 + 1e-9); qs = np.quantile(lm[gate], np.linspace(0, 1, 11))
dec = np.clip(np.searchsorted(qs, lm, side='right') - 1, 0, 9)
eb = [0, 0.25, 0.5, 1, 2, np.inf]; elab = ['<0.25', '0.25–0.5', '0.5–1', '1–2', '≥2']
alt = np.abs(np.nan_to_num(lt))
rec = np.full((5, 10), np.nan); cnt = np.zeros((5, 10), int)
for a in range(5):
    sa = T & (alt >= eb[a]) & (alt < eb[a + 1])
    for d in range(10):
        s = sa & (dec == d)[None, :]; cnt[a, d] = s.sum()
        if cnt[a, d] >= 30: rec[a, d] = (s & R).sum() / cnt[a, d]
im = ax.imshow(rec, cmap='Blues', vmin=0, vmax=0.5, aspect='auto', origin='lower')
for a in range(5):
    for d in range(10):
        if np.isfinite(rec[a, d]):
            ax.text(d, a, f'{rec[a, d]:.0%}', ha='center', va='center', fontsize=8.5,
                    color='white' if rec[a, d] > 0.3 else INK)
ax.set_xticks(range(10), [f'{k + 1}' for k in range(10)]); ax.set_yticks(range(5), elab); ax.grid(False)
ax.set_xlabel('control expression decile (1 = lowest, gated genes)'); ax.set_ylabel('|true log2 fold change|')
ax.set_title('B. Share of true DE genes the model calls (recall)', loc='left', color=INK)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.outline.set_visible(False)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))

ax = fig.add_subplot(gs[1, 1])
g = gate & np.isfinite(lt[i]) & ok[i]; y = lt[i]
groups = [(g & ~T[i] & ~R[i], NONE, 'neither', 2, 0.35), (g & T[i] & ~R[i], TRUE, 'true DE, missed', 4, 0.6),
          (g & ~T[i] & R[i], CALL, 'called, not true DE', 6, 0.85), (g & T[i] & R[i], HIT, 'hit', 6, 0.85)]
for s, c, lab, ms, al in groups:
    ax.scatter(lm[s], y[s], s=ms, color=c, alpha=al, lw=0, label=f'{lab} ({s.sum():,})')
ax.axhline(0, color=INK2, lw=0.6)
ax.set_xlabel('control mean CPM (log10)'); ax.set_ylabel('true log2 fold change (real H1 cells)')
ax.set_title(f'C. One knockdown ({EX}): which genes are called', loc='left', color=INK)
ax.legend(loc='upper right', frameon=False, markerscale=2.5, fontsize=9)
fig.suptitle('H1 held out (Eval A), incumbent C+G end-to-end run: which genes get called DE (p_adj < 0.05, own gene excluded)',
             x=0.01, ha='left', fontsize=13, color=INK)
fig.savefig(out / 'h1_calls.png', dpi=150, bbox_inches='tight'); fig.savefig(out / 'h1_calls.pdf', bbox_inches='tight')

# ---------- figure 2: true vs predicted Delta ----------
def beta(x, y):
    return float((x * y).sum() / (x * x).sum())

Q = z['req'] & ok
xa, ya = pdel[ok], tdel[ok]
fits = {'all': beta(xa, ya), 'called': beta(pdel[Q], tdel[Q]), 'not_called': beta(pdel[ok & ~Q], tdel[ok & ~Q])}
r_all = float(np.corrcoef(xa, ya)[0, 1])
bk = np.array([beta(pdel[k][ok[k]], tdel[k][ok[k]]) if (pdel[k][ok[k]] ** 2).sum() > 0 else np.nan for k in range(len(tgt))])

fig, (a1, a2) = plt.subplots(1, 2, figsize=(15, 6.6), gridspec_kw=dict(width_ratios=[1.1, 1], wspace=0.22))
hb = a1.hexbin(xa, ya, gridsize=140, bins='log', cmap='Blues', mincnt=1, extent=(-0.5, 0.55, -1.0, 1.05), lw=0)
cb = fig.colorbar(hb, ax=a1, fraction=0.04, pad=0.02); cb.set_label('(knockdown, gene) pairs'); cb.outline.set_visible(False)
xx = np.array([-0.5, 0.55])
a1.plot(xx, xx, '--', color=INK2, lw=1.2, label='perfect: true = predicted')
a1.plot(xx, fits['all'] * xx, color=INK, lw=2, label=f"fit, all genes: β = {fits['all']:.3f}")
a1.plot(xx, fits['called'] * xx, color=TRUE, lw=2, label=f"fit, called genes: β = {fits['called']:.3f}")
a1.plot(xx, fits['not_called'] * xx, color=CALL, lw=2, label=f"fit, not-called genes: β = {fits['not_called']:.3f}")
a1.set_xlim(-0.5, 0.55); a1.set_ylim(-1.0, 1.05); a1.grid(False)
a1.set_xlabel('predicted Δ (C+G, log1p CP50K)'); a1.set_ylabel('true Δ (held-out H1 cells, log1p CP50K)')
a1.set_title(f'A. 292 knockdowns × gated genes ({ok.sum() / 1e6:.2f}M pairs), r = {r_all:.3f}\n'
             'β = no-intercept OLS slope', loc='left', color=INK)
a1.legend(loc='upper left', frameon=False, fontsize=9)

good = np.isfinite(bk)
a2.scatter(np.maximum(nt[good], 0.8), bk[good], s=16, color=HIT, alpha=0.75, lw=0)
a2.axhline(fits['all'], color=INK, lw=1.2, label=f"pooled β = {fits['all']:.3f}")
a2.axhline(np.nanmedian(bk), color=INK2, lw=1.2, ls=':', label=f'median per-knockdown β = {np.nanmedian(bk):.3f}')
a2.axhline(1, color=INK2, lw=1, ls='--', label='β = 1 (magnitude right)')
a2.set_xscale('log'); a2.set_xlabel('true DE genes in the knockdown (0 drawn at 0.8)')
a2.set_ylabel('β, true Δ on predicted Δ (one knockdown)')
for name in (EX, 'DNMT1'):
    k = tgt.index(name); a2.annotate(name, (max(nt[k], 0.8), bk[k]), (6, 4), textcoords='offset points', color=INK2, fontsize=9)
a2.set_title('B. β per knockdown\nthe signal sits in the strong knockdowns', loc='left', color=INK)
a2.legend(loc='upper left', frameon=False, fontsize=9)
fig.suptitle('H1 held out (Eval A): true Δ regressed on the incumbent C+G predicted Δ', x=0.01, ha='left', fontsize=13, color=INK)
fig.savefig(out / 'h1_delta_beta.png', dpi=150, bbox_inches='tight'); fig.savefig(out / 'h1_delta_beta.pdf', bbox_inches='tight')

summ = dict(model='incumbent C+G (stage-1 arm C + Gaussian shrink), Eval A, H1 line held out',
            calls='realised calls from the end-to-end harness run (de_pred p_adj < 0.05), own gene out, gated genes',
            delta='predicted = mean-inversion proxy (bulk of mu_from_lfc), true = H1 pseudobulk - control, log1p CP50K',
            median_true_de=float(np.median(nt)), median_calls=float(np.median(nr)), pooled_recall=float(nh.sum() / nt.sum()),
            pooled_precision=float(nh.sum() / nr.sum()), beta=fits, r_all=r_all,
            beta_per_knockdown_quantiles={q: float(np.nanpercentile(bk, q)) for q in (5, 25, 50, 75, 95)},
            beta_called_note='called = requested de_sig (the K0 call set config 1 also uses)')
(out / 'summary.json').write_text(json.dumps(summ, indent=1) + '\n')
print(json.dumps(summ, indent=1))

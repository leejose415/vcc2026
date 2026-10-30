"""Illustrations for Jose (2026-10-30), static PNGs in results/stage1/2026-10-30_illustrations/.

  fig1  log1p profiles: per line, genes ranked by the control profile; grey = control, orange = every knockdown's
        profile. H1 (Eval A) and HCT116final from the harness's own real side (pseudobulk log1p of CP50K, the scorer's
        space); context D from its 18,400 control cells the same way (no knockdowns).
  fig2  change: Delta = knockdown log1p - control log1p (the scorer's space), same gene order; true DE up (red) / down
        (blue) / not DE (grey) from the eval's real DE table (p_adj < 0.05, own gene out).
  fig3  DE calls of the incumbent C+G end-to-end runs (realised: harness de_pred p_adj < 0.05, own gene out) against the
        truth, per knockdown sorted by true DE count: correct (true DE, same sign), wrong direction, false call (not
        true DE), missed (true DE, not called); plus totals.

  ~/venv/bin/python code/vcc6/illustrations.py
"""
import glob
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                                            # noqa: E402
import numpy as np                                                                         # noqa: E402
import pandas as pd                                                                        # noqa: E402

V = Path.home() / 'vcc'
OUT = V / 'results/stage1/2026-10-30_illustrations'
GENES = np.array(pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str))
CTRL, BTS = 'non-targeting', 5e4
LINES = {
    'H1 (Eval A)': dict(eval=V / 'lib/harness/evalA',
                        run=V / 'results/stage1/2026-10-27_incumbent/evalA/CG/e2e/salt_none/harness/evalA'),
    'HCT116 (HCT116final)': dict(eval=V / 'results/stage1/final_evals/harness/lib/evalHCT116final',
                                 run=V / 'results/stage1/e2e_final/CG/HCT116/salt_none/harness/evalHCT116final'),
}
# reference palette (dataviz skill, references/palette.md, light mode)
GREY, ORANGE, BLUE, AQUA, RED = '#8f8d87', '#eb6834', '#2a78d6', '#1baf7a', '#e34948'
MISSED = '#c3c2b7'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'axes.labelcolor': INK2, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.axisbelow': True,
                     'figure.facecolor': '#fcfcfb', 'axes.facecolor': '#fcfcfb', 'savefig.facecolor': '#fcfcfb'})


def real_side(d):
    z = np.load(glob.glob(str(d / 'real_cache' / '*.moments.npz'))[0], allow_pickle=True)
    perts = z['perts'].astype(str); means = z['means'].astype(np.float64)
    i = list(perts).index(CTRL); kd = [j for j in range(len(perts)) if j != i]
    return perts[kd], means[i], means[kd], int(z['counts'][i])


def truth(d, targets):
    t = pd.read_parquet(glob.glob(str(d / 'real_cache' / 'de_wilcoxon_table-*.parquet'))[0],
                        columns=['target', 'feature', 'log2_fold_change', 'p_adj'])
    t = t[(t.p_adj < 0.05) & (t.feature != t.target)]
    tp = {k: i for i, k in enumerate(targets)}; gp = {g: j for j, g in enumerate(GENES)}
    S = np.zeros((len(targets), len(GENES)), np.int8)                          # +1 up, -1 down, 0 not DE
    S[t.target.map(tp).to_numpy(), t.feature.map(gp).to_numpy()] = np.sign(t.log2_fold_change.to_numpy()).astype(np.int8)
    return S


def calls(run, targets):
    t = pd.read_parquet(run / 'de_pred.parquet', columns=['target', 'feature', 'log2_fold_change', 'p_adj'])
    t = t[(t.p_adj < 0.05) & (t.feature != t.target)]
    tp = {k: i for i, k in enumerate(targets)}; gp = {g: j for j, g in enumerate(GENES)}
    C = np.zeros((len(targets), len(GENES)), np.int8)
    C[t.target.map(tp).to_numpy(), t.feature.map(gp).to_numpy()] = np.where(t.log2_fold_change.to_numpy() >= 0, 1, -1)
    return C


def context_d():
    import anndata as ad
    a = ad.read_h5ad(V / 'data_final/extracted/context_D.h5ad')
    assert list(a.var_names.astype(str)) == list(GENES)
    s = np.asarray(a.X.sum(0)).ravel().astype(np.float64)
    return np.log1p(BTS * s / s.sum()), a.n_obs


def fig1(data, dprof, dn):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), sharey=True)
    rng = np.random.default_rng(2026)
    for ax, (name, d) in zip(axes, data.items()):
        o = np.argsort(d['ctrl'], kind='stable'); x = np.arange(len(o))
        K = d['kd'][:, o]
        xs = np.broadcast_to(x, K.shape).ravel(); ys = K.ravel()
        keep = rng.random(len(xs)) < 0.25                                     # a quarter of the (knockdown, gene) dots
        ax.scatter(xs[keep], ys[keep], s=0.15, c=ORANGE, alpha=0.06, linewidths=0, rasterized=True)
        ax.scatter(x, d['ctrl'][o], s=0.6, c=GREY, linewidths=0, rasterized=True)
        ax.set_title(f'{name}\n{K.shape[0]} knockdowns × 18,533 genes; controls: {d["n_ctrl"]:,} cells', color=INK, fontsize=10.5, loc='left')
        ax.set_xlabel('genes, ranked by control log1p')
    o = np.argsort(dprof, kind='stable')
    axes[2].scatter(np.arange(len(o)), dprof[o], s=0.6, c=GREY, linewidths=0, rasterized=True)
    axes[2].set_title(f'Context D (final test)\ncontrols only: {dn:,} cells; no knockdown truth', color=INK, fontsize=10.5, loc='left')
    axes[2].set_xlabel('genes, ranked by control log1p')
    axes[0].set_ylabel('log1p(pseudobulk CP50K), the scorer\'s space')
    h = [plt.Line2D([], [], ls='', marker='o', ms=6, mfc=GREY, mec='none', label='control profile (one dot per gene)'),
         plt.Line2D([], [], ls='', marker='o', ms=6, mfc=ORANGE, mec='none', label='knockdown profiles (one dot per knockdown × gene; 25% shown)')]
    fig.legend(handles=h, loc='upper left', frameon=False, ncol=2, bbox_to_anchor=(0.005, 0.955))
    fig.suptitle('Expression profiles: control vs knockdown', x=0.01, ha='left', color=INK, fontsize=13, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.89)); fig.savefig(OUT / 'fig1_log1p_profiles.png', dpi=170); plt.close(fig)


def fig2(data):
    fig, axes = plt.subplots(2, 2, figsize=(15, 8.6), gridspec_kw=dict(width_ratios=[3.2, 1]))
    stats, lims = {}, []
    for r, (name, d) in enumerate(data.items()):
        ax, axh = axes[r]
        o = np.argsort(d['ctrl'], kind='stable'); x = np.arange(len(o))
        D = (d['kd'] - d['ctrl'][None, :])[:, o]; S = d['truth'][:, o]
        xs = np.broadcast_to(x, D.shape)
        n = S == 0
        rng = np.random.default_rng(2026); keep = n & (rng.random(D.shape) < 0.15)
        ax.scatter(xs[keep], D[keep], s=0.15, c=GREY, alpha=0.08, linewidths=0, rasterized=True)
        for sgn, col in ((1, RED), (-1, BLUE)):
            m = S == sgn
            ax.scatter(xs[m], D[m], s=0.5, c=col, alpha=0.35, linewidths=0, rasterized=True)
        ax.axhline(0, color=INK2, lw=0.8)
        lims.append(np.quantile(np.abs(D), 0.9999))
        nde = (S != 0).sum(1)
        stats[name] = dict(true_de_total=int((S != 0).sum()), up=int((S == 1).sum()), down=int((S == -1).sum()),
                           true_de_per_kd_median=float(np.median(nde)), true_de_per_kd_mean=float(nde.mean()),
                           kd_with_0_true_de=int((nde == 0).sum()),
                           abs_delta_de_median=float(np.median(np.abs(D[S != 0]))) if (S != 0).any() else None,
                           abs_delta_notde_median=float(np.median(np.abs(D[S == 0]))))
        st = stats[name]
        ax.set_title(f'{name}: Δ = knockdown − control, every knockdown × gene\n'
                     f'true DE: {st["true_de_total"]:,} ({st["up"]:,} up, {st["down"]:,} down); median {st["true_de_per_kd_median"]:.0f} per knockdown, '
                     f'{st["kd_with_0_true_de"]} of {len(nde)} knockdowns with none', color=INK, fontsize=10.5, loc='left')
        ax.set_ylabel('Δ log1p (scorer\'s space)'); ax.set_xlabel('genes, ranked by control log1p (as in figure 1)')
        # right: per-knockdown true DE count, sorted, split up / down
        up = (S == 1).sum(1); dn = (S == -1).sum(1); oo = np.argsort(up + dn)[::-1]; xk = np.arange(len(oo))
        axh.bar(xk, up[oo], width=1.0, color=RED, linewidth=0, label='up')
        axh.bar(xk, -dn[oo], width=1.0, color=BLUE, linewidth=0, label='down')
        axh.axhline(0, color=INK2, lw=0.8)
        axh.set_title('true DE genes per knockdown\n(sorted; up above, down below)', color=INK, fontsize=10.5, loc='left')
        axh.set_xlabel('knockdowns'); axh.set_yscale('symlog', linthresh=10)
    for ax in axes[:, 0]:
        ax.set_ylim(-max(lims) * 1.05, max(lims) * 1.05)                      # one scale for both lines
    h = [plt.Line2D([], [], ls='', marker='o', ms=6, mfc=RED, mec='none', label='true DE, up'),
         plt.Line2D([], [], ls='', marker='o', ms=6, mfc=BLUE, mec='none', label='true DE, down'),
         plt.Line2D([], [], ls='', marker='o', ms=6, mfc=GREY, mec='none', label='not DE (15% shown)')]
    fig.legend(handles=h, loc='upper left', frameon=False, ncol=3, bbox_to_anchor=(0.005, 0.965))
    fig.suptitle('How far, and which way, knockdowns move genes', x=0.01, ha='left', color=INK, fontsize=13, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(OUT / 'fig2_delta.png', dpi=170); plt.close(fig)
    return stats


def fig3(data):
    fig, axes = plt.subplots(2, 2, figsize=(15, 8.6), gridspec_kw=dict(width_ratios=[3.2, 1]))
    stats = {}
    for r, (name, d) in enumerate(data.items()):
        ax, axt = axes[r]; S, C = d['truth'], d['calls']
        cor = ((C != 0) & (S != 0) & (C == S)).sum(1); wrong = ((C != 0) & (S != 0) & (C != S)).sum(1)
        false = ((C != 0) & (S == 0)).sum(1); miss = ((C == 0) & (S != 0)).sum(1)
        o = np.argsort((S != 0).sum(1), kind='stable')[::-1]; x = np.arange(len(o))
        ax.bar(x, cor[o], width=1.0, color=BLUE, linewidth=0)
        ax.bar(x, wrong[o], bottom=cor[o], width=1.0, color=AQUA, linewidth=0)
        ax.bar(x, false[o], bottom=(cor + wrong)[o], width=1.0, color=ORANGE, linewidth=0)
        ax.bar(x, -miss[o], width=1.0, color=MISSED, linewidth=0)
        ax.axhline(0, color=INK2, lw=0.8)
        ax.set_ylabel('genes per knockdown\n(called above 0, missed below)')
        ax.set_xlabel('knockdowns, sorted by true DE count (strongest left)')
        tot = dict(correct=int(cor.sum()), wrong_direction=int(wrong.sum()), false_call=int(false.sum()), missed=int(miss.sum()))
        tot['called'] = tot['correct'] + tot['wrong_direction'] + tot['false_call']; tot['true_de'] = tot['correct'] + tot['wrong_direction'] + tot['missed']
        tot['precision'] = (tot['correct'] + tot['wrong_direction']) / max(tot['called'], 1)
        tot['recall'] = (tot['correct'] + tot['wrong_direction']) / max(tot['true_de'], 1)
        tot['sign_accuracy_on_hits'] = tot['correct'] / max(tot['correct'] + tot['wrong_direction'], 1)
        tot['per_kd_median'] = dict(called=float(np.median(cor + wrong + false)), true_de=float(np.median((S != 0).sum(1))),
                                    correct=float(np.median(cor)), false_call=float(np.median(false)), missed=float(np.median(miss)))
        stats[name] = tot
        ax.set_title(f'{name}: C+G end-to-end calls (realised, p_adj < 0.05), per knockdown\n'
                     f'precision {tot["precision"]:.2f}, recall {tot["recall"]:.2f}, sign right on {tot["sign_accuracy_on_hits"]:.0%} of hits; '
                     f'median per knockdown: {tot["per_kd_median"]["called"]:.0f} called, {tot["per_kd_median"]["true_de"]:.0f} true DE',
                     color=INK, fontsize=10.5, loc='left')
        # right: totals
        labs = ['correct', 'wrong direction', 'false call', 'missed']; vals = [tot['correct'], tot['wrong_direction'], tot['false_call'], tot['missed']]
        cols = [BLUE, AQUA, ORANGE, MISSED]
        yy = np.arange(4)[::-1]
        axt.barh(yy, vals, color=cols, height=0.62)
        for y_, v in zip(yy, vals):
            axt.text(v, y_, f'  {v:,}', va='center', ha='left', color=INK, fontsize=9.5)
        axt.set_yticks(yy); axt.set_yticklabels(labs); axt.set_xlim(0, max(vals) * 1.45)
        axt.grid(axis='y', visible=False); axt.set_xlabel('genes, all knockdowns')
        axt.set_title('totals', color=INK, fontsize=10.5, loc='left')
    h = [plt.Rectangle((0, 0), 1, 1, color=c, label=l) for c, l in
         ((BLUE, 'correct: true DE, same direction'), (AQUA, 'wrong direction: true DE, opposite sign'),
          (ORANGE, 'false call: called, not true DE'), (MISSED, 'missed: true DE, not called'))]
    fig.legend(handles=h, loc='upper left', frameon=False, ncol=4, bbox_to_anchor=(0.005, 0.965))
    fig.suptitle('DE calls against the truth (incumbent C+G)', x=0.01, ha='left', color=INK, fontsize=13, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(OUT / 'fig3_de_calls.png', dpi=170); plt.close(fig)
    return stats


EPS = 0.01                                     # CP50K pseudo-count for fig5's LFC
LIGHT = '#d9d8d2'


def tested_genes(d):
    """Genes in the scorer's real DE table (its own expression gate)."""
    t = pd.read_parquet(glob.glob(str(d / 'real_cache' / 'de_wilcoxon_table-*.parquet'))[0], columns=['feature'])
    return np.isin(GENES, t.feature.unique())


def short(name):
    return name.split(' ')[0]


def fig45(data, evals):
    stats = {}
    xmax = max(d['ctrl'].max() for d in data.values()) * 1.02
    lim_d = max(np.quantile(np.abs(d['kd'] - d['ctrl'][None, :]), 0.9999) for d in data.values()) * 1.05
    lfcs = {n: np.log2((np.expm1(d['kd']) + EPS) / (np.expm1(d['ctrl'])[None, :] + EPS)) for n, d in data.items()}
    lim_l = max(np.quantile(np.abs(v), 0.9999) for v in lfcs.values()) * 1.05
    for name, d in data.items():
        tg = tested_genes(evals[name]); S = d['truth']; x = np.broadcast_to(d['ctrl'], S.shape)
        rng = np.random.default_rng(2026); sub = rng.random(S.shape) < 0.15
        for tag, Y, lim, ylab, ttl in (
                ('fig4', d['kd'] - d['ctrl'][None, :], lim_d, "Δ log1p (knockdown − control, scorer's space)",
                 'Δ against control expression'),
                ('fig5', lfcs[name], lim_l, f'log2 fold change (pseudobulk, ε = {EPS} CP50K)',
                 'log2 fold change against control expression')):
            fig, ax = plt.subplots(figsize=(11, 6.2))
            untested = ~tg[None, :] & (S == 0)
            ax.scatter(x[untested & sub], Y[untested & sub], s=0.15, c=LIGHT, alpha=0.25, linewidths=0, rasterized=True)
            nd = tg[None, :] & (S == 0) & sub
            ax.scatter(x[nd], Y[nd], s=0.15, c=GREY, alpha=0.10, linewidths=0, rasterized=True)
            for sgn, col in ((1, RED), (-1, BLUE)):
                m = S == sgn
                ax.scatter(x[m], Y[m], s=0.5, c=col, alpha=0.35, linewidths=0, rasterized=True)
            ax.axhline(0, color=INK2, lw=0.8)
            ax.set_xlim(-0.1, xmax); ax.set_ylim(-lim, lim)
            ax.set_xlabel("control expression: log1p(pseudobulk CP50K), the scorer's space")
            ax.set_ylabel(ylab)
            nz = int((d['ctrl'] == 0).sum())
            ax.set_title(f'{name}: {S.shape[0]} knockdowns × 18,533 genes; {int(tg.sum()):,} genes in the scorer\'s DE table, '
                         f'{nz:,} with zero control expression (x = 0)', color=INK, fontsize=10, loc='left')
            h = [plt.Line2D([], [], ls='', marker='o', ms=6, mfc=RED, mec='none', label='true DE, up'),
                 plt.Line2D([], [], ls='', marker='o', ms=6, mfc=BLUE, mec='none', label='true DE, down'),
                 plt.Line2D([], [], ls='', marker='o', ms=6, mfc=GREY, mec='none', label='not DE, tested (15% shown)'),
                 plt.Line2D([], [], ls='', marker='o', ms=6, mfc=LIGHT, mec='none', label='not in the DE table (15% shown)')]
            fig.legend(handles=h, loc='upper left', frameon=False, ncol=4, bbox_to_anchor=(0.005, 0.935))
            fig.suptitle(f'{short(name)} — {ttl}', x=0.01, ha='left', color=INK, fontsize=13, y=0.995)
            fig.tight_layout(rect=(0, 0, 1, 0.88))
            fn = OUT / f'{tag}_{short(name)}_{"delta" if tag == "fig4" else "lfc"}_vs_control.png'
            fig.savefig(fn, dpi=170); plt.close(fig)
            lo = d['ctrl'] < np.quantile(d['ctrl'][tg], 0.1)                # the lowest-expressed tenth of tested genes
            stats[f'{tag}_{short(name)}'] = dict(
                file=fn.name, abs_median_all=float(np.median(np.abs(Y))),
                abs_median_true_de=float(np.median(np.abs(Y[S != 0]))) if (S != 0).any() else None,
                abs_q99_low_expr=float(np.quantile(np.abs(Y[:, lo]), 0.99)),
                abs_q99_high_expr=float(np.quantile(np.abs(Y[:, d['ctrl'] > np.quantile(d['ctrl'][tg], 0.9)]), 0.99)),
                abs_q99_untested=float(np.quantile(np.abs(Y[:, ~tg]), 0.99)))
    return stats


FIG6_KD = ['STUB1', 'TARBP2', 'SHPRH', 'PLAGL2', 'SLIRP', 'STAT6', 'MED13', 'ACLY', 'BRPF1', 'SIN3B']
PURPLE = '#8e5bd0'


def fig6(data):
    (a, da), (b, db) = list(data.items())
    stats = {}
    panels = []
    for t in FIG6_KD:
        i, j = list(da['targets']).index(t), list(db['targets']).index(t)
        x = da['kd'][i] - da['ctrl']; y = db['kd'][j] - db['ctrl']
        sa, sb = da['truth'][i] != 0, db['truth'][j] != 0
        own = GENES == t; keep = ~own
        r = float(np.corrcoef(x[keep], y[keep])[0, 1])
        cos = float(x[keep] @ y[keep] / (np.linalg.norm(x[keep]) * np.linalg.norm(y[keep])))
        slope = float(x[keep] @ y[keep] / (x[keep] @ x[keep]))
        both = sa & sb
        st = dict(true_de_H1=int(sa.sum()), true_de_HCT116=int(sb.sum()), true_de_both=int(both.sum()),
                  pearson_r=r, cosine=cos, slope_through_origin=slope,
                  norm_ratio_HCT116_over_H1=float(np.linalg.norm(y[keep]) / np.linalg.norm(x[keep])),
                  sign_agreement_on_both_DE=float((np.sign(x[both]) == np.sign(y[both])).mean()) if both.any() else None,
                  target_delta_H1=float(x[own][0]) if own.any() else None, target_delta_HCT116=float(y[own][0]) if own.any() else None)
        stats[t] = st; panels.append((t, x, y, sa, sb, own, st))

    def draw(ax, t, x, y, sa, sb, own, st, small=False):
        nn = ~sa & ~sb & ~own
        ax.scatter(x[nn], y[nn], s=1.2 if small else 2, c=GREY, alpha=0.25, linewidths=0, rasterized=True)
        for m, col in ((sa & ~sb, BLUE), (sb & ~sa, AQUA), (sa & sb, ORANGE)):
            ax.scatter(x[m & ~own], y[m & ~own], s=3 if small else 5, c=col, alpha=0.7, linewidths=0, rasterized=True)
        lim = max(np.quantile(np.abs(np.r_[x[~own], y[~own]]), 0.9999), 0.05) * 1.1
        if own.any():                                   # the target: pinned to the edge when off-scale, values printed
            tx, ty = float(x[own][0]), float(y[own][0]); cx, cy = np.clip([tx, ty], -0.97 * lim, 0.97 * lim)
            off = (cx, cy) != (tx, ty)
            ax.scatter([cx], [cy], s=40 if small else 70, marker='*', c=INK, linewidths=0, zorder=5, clip_on=False)
            ax.text(0.98, 0.03, f'★ target Δ: H1 {tx:+.2f}, HCT116 {ty:+.2f}' + (' (off-scale, pinned to edge)' if off else ''),
                    transform=ax.transAxes, ha='right', va='bottom', fontsize=7.5 if small else 9, color=INK2)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
        ax.plot([-lim, lim], [-lim, lim], color=INK2, lw=0.7, ls='--', zorder=1)
        ax.axhline(0, color=GRID, lw=0.8, zorder=0); ax.axvline(0, color=GRID, lw=0.8, zorder=0)
        ax.set_title(f'{t}  (true DE: H1 {st["true_de_H1"]:,}, HCT116 {st["true_de_HCT116"]:,})\n'
                     f'r {st["pearson_r"]:.2f} · cosine {st["cosine"]:.2f} · slope {st["slope_through_origin"]:.2f}',
                     color=INK, fontsize=9 if small else 10.5, loc='left')

    handles = [plt.Line2D([], [], ls='', marker='o', ms=6, mfc=ORANGE, mec='none', label='true DE in both'),
               plt.Line2D([], [], ls='', marker='o', ms=6, mfc=BLUE, mec='none', label='true DE in H1 only'),
               plt.Line2D([], [], ls='', marker='o', ms=6, mfc=AQUA, mec='none', label='true DE in HCT116 only'),
               plt.Line2D([], [], ls='', marker='o', ms=6, mfc=GREY, mec='none', label='not DE in either'),
               plt.Line2D([], [], ls='', marker='*', ms=9, mfc=INK, mec='none', label='the knocked-down gene'),
               plt.Line2D([], [], ls='--', color=INK2, label='y = x')]
    fig, axes = plt.subplots(2, 5, figsize=(22, 10))
    for ax, p_ in zip(axes.ravel(), panels):
        draw(ax, *p_, small=True)
    for ax in axes[1]:
        ax.set_xlabel('Δ in H1 (Eval A)')
    for ax in axes[:, 0]:
        ax.set_ylabel('Δ in HCT116 (HCT116final)')
    fig.legend(handles=handles, loc='upper left', frameon=False, ncol=6, bbox_to_anchor=(0.005, 0.955))
    fig.suptitle('Same knockdown, two lines: Δ log1p (knockdown − own controls, scorer\'s space), 18,533 genes per panel; '
                 'ordered by H1 true DE count', x=0.01, ha='left', color=INK, fontsize=13, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.92)); fig.savefig(OUT / 'fig6_cross_line_delta_grid.png', dpi=150); plt.close(fig)
    for p_ in panels:
        fig, ax = plt.subplots(figsize=(7.2, 7.6))
        draw(ax, *p_)
        ax.set_xlabel('Δ in H1 (Eval A)'); ax.set_ylabel('Δ in HCT116 (HCT116final)')
        fig.legend(handles=handles, loc='upper left', frameon=False, ncol=3, bbox_to_anchor=(0.01, 0.995), fontsize=8.5)
        fig.tight_layout(rect=(0, 0, 1, 0.91)); fig.savefig(OUT / f'fig6_{p_[0]}_H1_vs_HCT116_delta.png', dpi=150); plt.close(fig)
    return stats


def main(only=None):
    OUT.mkdir(parents=True, exist_ok=True)
    data = {}
    for name, s in LINES.items():
        targets, ctrl, kd, n_ctrl = real_side(s['eval'])
        data[name] = dict(targets=targets, ctrl=ctrl, kd=kd, n_ctrl=n_ctrl, truth=truth(s['eval'], targets),
                          calls=calls(s['run'], targets))
    if only:
        st = json.loads((OUT / 'stats.json').read_text()) if (OUT / 'stats.json').exists() else {}
        if 'fig4' in only or 'fig5' in only:
            st['fig4_fig5'] = fig45(data, {n: s['eval'] for n, s in LINES.items()}); print(json.dumps(st['fig4_fig5'], indent=1))
        if 'fig6' in only:
            st['fig6'] = fig6(data); print(json.dumps(st['fig6'], indent=1))
        (OUT / 'stats.json').write_text(json.dumps(st, indent=1) + '\n')
        return
    dprof, dn = context_d()
    fig1(data, dprof, dn)
    st = dict(fig2=fig2(data), fig3=fig3(data), fig4_fig5=fig45(data, {n: s['eval'] for n, s in LINES.items()}), fig6=fig6(data),
              fig1={**{k: dict(knockdowns=int(v['kd'].shape[0]), control_cells=v['n_ctrl'],
                               genes_zero_in_control=int((v['ctrl'] == 0).sum())) for k, v in data.items()},
                    'Context D': dict(control_cells=dn, genes_zero_in_control=int((dprof == 0).sum()))})
    (OUT / 'stats.json').write_text(json.dumps(st, indent=1) + '\n')
    print(json.dumps(st, indent=1))


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument('--only', nargs='*')
    main(ap.parse_args().only)

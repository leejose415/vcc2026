"""Figures for Step 1 on the final screen (results/stage1/step1_scale/final/<line>.json, from step_scale.py --set final).

  fig 1  step1_final_metrics.png  four metrics x four lines; C+G, C+G x beta, E, M as dots (median over knockdowns with
         its bootstrap 95% CI; pb_nmse is a line-level ratio, no CI); the gray bar is the reference each metric is
         judged against (predict nothing / coin flip / the scorer's mean-response baseline).
  fig 3  step1_final_r2.png      left: out-of-sample r^2 (gene-level, pooled over the final targets; the same for C+G,
         beta and E, a positive scale cannot change it) against the noise ceiling, the share of the truth's variance
         that is not measurement noise (1 - mean(se^2) / var(lfc), same genes; cached in r2_ceiling.json); right:
         R^2 = 1 - MSE / var(true lfc), in and out of sample (0 = as good as predicting the mean, about predicting 0).
  fig 2  step1_final_scale.png   left: beta fitted on the held-out source lines vs the scale the target line itself
         wants (C+G's out-of-sample slope); right: the E curves (scale by the knocked-down gene's expression) with beta.

  ~/venv/bin/python code/vcc6/plot_step_scale.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                                           # noqa: E402
import numpy as np                                                                        # noqa: E402

V = Path.home() / 'vcc'
IN = V / 'results/stage1/step1_scale/final'
OUT = IN / 'figures'
LINES = ['K562', 'HCT116', 'HEK293T', 'H1']
VARIANTS = [('C+G', 'C+G', '#2a78d6', 'o'), ('beta', 'C+G × β', '#eb6834', 's'),
            ('E', 'E', '#1baf7a', '^'), ('M', 'M', '#eda100', 'D')]
INK, INK2, MUTED, GRID, SURFACE, REF = '#0b0b0b', '#52514e', '#8a8984', '#e6e5e0', '#fcfcfb', '#c9c8c2'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK2, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.titlecolor': INK, 'figure.facecolor': SURFACE,
                     'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE})


def load():
    return {n: json.loads((IN / f'{n}.json').read_text()) for n in LINES}


def style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='y', color=GRID, lw=0.8); ax.set_axisbelow(True)
    ax.tick_params(length=0)


def metrics_figure(R):
    panels = [('nmae', 'NMAE (lower is better)', 'predict nothing = 1'),
              ('mse', 'MSE, effect space (lower is better)', 'predict nothing'),
              ('sign_acc_true_de', 'Sign accuracy on true DE genes (higher is better)', 'coin flip = 0.5'),
              ('pb_nmse', 'pb_nmse, scorer space (lower is better)', "scorer's mean-response baseline")]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.2))
    x = np.arange(len(LINES)); off = np.linspace(-0.27, 0.27, len(VARIANTS))
    for ax, (m, title, ref_label) in zip(axes.ravel(), panels):
        for j, n in enumerate(LINES):
            cg = R[n]['summary']['C+G']
            ref = {'nmae': 1.0, 'mse': cg['all']['mse_zero'], 'sign_acc_true_de': 0.5,
                   'pb_nmse': cg['scorer_space']['all']['pb_nmse_meanresp']}[m]
            ax.hlines(ref, j - 0.42, j + 0.42, color=REF, lw=5, zorder=1,
                      label=ref_label if j == 0 else None)
            for k, (v, lab, col, mk) in enumerate(VARIANTS):
                s = R[n]['summary'][v]
                if m == 'pb_nmse':
                    y, ci = s['scorer_space']['all']['pb_nmse'], None
                else:
                    y, ci = s['all'][m], s['all'][f'{m}_ci']
                if ci is not None:
                    ax.vlines(x[j] + off[k], ci[0], ci[1], color=col, lw=2, zorder=2)
                ax.scatter(x[j] + off[k], y, s=64, marker=mk, color=col, edgecolor=SURFACE, linewidth=2, zorder=3,
                           label=lab if j == 0 else None)
        ax.set_xticks(x, [f'{n}\n(reported only)' if n == 'H1' else n for n in LINES])
        ax.set_title(title, loc='left', fontsize=11)
        ax.set_xlim(-0.6, len(LINES) - 0.4)
        style(ax)
        ax.legend(frameon=False, fontsize=8.5, loc='best', handletextpad=0.3, labelcolor=INK2)
    betas = ', '.join(f'{n} {R[n]["curves"]["beta"]:.3f}' for n in LINES)
    fig.suptitle('Step 1 on the final screen: C+G vs one β, E and M', x=0.01, ha='left', fontsize=13, color=INK)
    fig.text(0.01, 0.935, f'Median over the final-target knockdowns, 95% bootstrap CI; pb_nmse is a ratio over the '
             f'line. β fitted on the held-out source lines: {betas}.', fontsize=9, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(OUT / 'step1_final_metrics.png', dpi=150)


def scale_figure(R):
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw=dict(width_ratios=[1, 1.25]))
    x = np.arange(len(LINES))
    fitted = [R[n]['curves']['beta'] for n in LINES]
    wanted = [R[n]['fit_out_of_sample_target']['C+G']['beta'] for n in LINES]
    a.bar(x - 0.19, fitted, 0.36, color='#2a78d6', label='β fitted on held-out source lines')
    a.bar(x + 0.19, wanted, 0.36, color='#eb6834', label="scale the target line wants (C+G's own slope)")
    for i in range(len(LINES)):
        a.text(x[i] - 0.19, fitted[i] + 0.004, f'{fitted[i]:.3f}', ha='center', va='bottom', fontsize=8.5, color=INK2)
        a.text(x[i] + 0.19, wanted[i] + 0.004, f'{wanted[i]:.3f}', ha='center', va='bottom', fontsize=8.5, color=INK2)
    a.set_xticks(x, LINES); a.set_ylabel('scale applied to C+G')
    a.set_title('Every line wants less than the fitted β', loc='left', fontsize=11)
    a.legend(frameon=False, fontsize=8.5, loc='upper left', labelcolor=INK2); a.set_ylim(0, max(fitted) * 1.35)
    style(a)
    cols = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']
    for n, col in zip(LINES, cols):
        e = R[n]['curves']['E']
        b.plot(e['expr_log1p_cpm'], e['scale'], color=col, lw=2, label=f'{n}: E curve')
        b.hlines(R[n]['curves']['beta'], 0, 8, color=col, lw=1, linestyles=(0, (4, 3)))
    b.set_xlabel('log1p CPM of the knocked-down gene in the control cells'); b.set_ylabel('scale')
    b.set_title('E stays within about 0.03 of β (dashed) for most knockdowns', loc='left', fontsize=11)
    b.set_ylim(0, None); b.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc='upper left')
    style(b)
    fig.tight_layout()
    fig.savefig(OUT / 'step1_final_scale.png', dpi=150)


def r2_ceiling():
    """1 - mean(se^2) / var(true lfc) over the tested genes of each line's final targets (own gene excluded)."""
    path = IN / 'r2_ceiling.json'
    if path.exists():
        return json.loads(path.read_text())
    import sys
    import anndata as ad
    sys.path.insert(0, str(V / 'code/vcc6'))
    import screen as SC
    import stage1 as S
    out = {}
    for n in LINES:
        t = V / f'lib/effects_v5/{SC.FINAL[n]}.h5ad'
        a = ad.read_h5ad(t, backed='r')
        genes = np.array(a.var_names.astype(str), dtype=object)
        _, gate = S.load_baseline(t.with_name(t.stem + '_baseline.csv'), genes)
        row = {k: i for i, k in enumerate(a.obs['target_gene'].astype(str))}; col = {g: j for j, g in enumerate(genes)}
        ys, ses = [], []
        for k in (V / f'results/target_lists/final_{n}.txt').read_text().split():
            y = np.asarray(a.layers['lfc'][row[k]], float); se = np.asarray(a.layers['se'][row[k]], float)
            ok = gate & np.isfinite(y) & np.isfinite(se)
            if k in col: ok[col[k]] = False
            ys.append(y[ok]); ses.append(se[ok])
        y, se = np.concatenate(ys), np.concatenate(ses)
        out[n] = dict(pairs=int(y.size), var_lfc=float(y.var()), mean_se2=float(np.mean(se ** 2)),
                      ceiling=float(max(0.0, 1 - np.mean(se ** 2) / y.var())))
    path.write_text(json.dumps(out, indent=1) + '\n')
    return out


def r2_figure(R):
    ceil = r2_ceiling()
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 5.8))
    x = np.arange(len(LINES))
    c = [ceil[n]['ceiling'] for n in LINES]
    g = [R[n]['summary']['C+G']['fit']['shrunk']['r2_fit'] for n in LINES]
    m = [R[n]['summary']['M']['fit']['shrunk']['r2_fit'] for n in LINES]
    a.bar(x - 0.27, c, 0.26, color=REF, label='noise ceiling (signal share of the truth)')
    a.bar(x, g, 0.26, color='#2a78d6', label='C+G = C+G × β = E')
    a.bar(x + 0.27, m, 0.26, color='#eda100', label='M')
    for i in range(len(LINES)):
        for dx, v in ((-0.27, c[i]), (0, g[i]), (0.27, m[i])):
            a.text(x[i] + dx, v * 1.15, f'{v:.4f}' if v < 0.01 else f'{v:.3f}', ha='center', va='bottom', fontsize=7.5,
                   color=INK2, rotation=90)
    a.set_yscale('log'); a.set_ylim(1e-4, 3)
    a.set_xticks(x, [f'{n}\n(reported only)' if n == 'H1' else n for n in LINES]); a.set_ylabel('r² (log scale)')
    a.set_title('Out of sample r² vs what the truth\'s noise allows', loc='left', fontsize=11)
    a.legend(frameon=False, fontsize=8.5, loc='upper left', labelcolor=INK2, ncol=1)
    style(a)
    off = np.linspace(-0.33, 0.27, len(VARIANTS))
    b.axhline(0, color=REF, lw=5, zorder=1, label='0 = as good as predicting the mean')
    for k, (v, lab, col, mk) in enumerate(VARIANTS):
        key = 'G' if v == 'C+G' else v
        ins = [R[n]['fit_in_sample_source_lines'][key]['r2_raw'] for n in LINES]
        oos = [R[n]['summary'][v]['fit']['shrunk']['r2_raw'] for n in LINES]
        b.scatter(x + off[k], oos, s=64, marker=mk, color=col, edgecolor=SURFACE, linewidth=2, zorder=3,
                  label=f'{lab}: out of sample (target line)')
        b.scatter(x + off[k] + 0.07, ins, s=40, marker=mk, facecolor='none', edgecolor=col, linewidth=1.5, zorder=3,
                  label=f'{lab}: in sample (held-out source lines)')
    b.set_xticks(x, [f'{n}\n(reported only)' if n == 'H1' else n for n in LINES])
    b.set_ylabel('R² = 1 − MSE / var(true lfc)')
    b.set_title('R² with the scale counted: shrinking only reaches about 0', loc='left', fontsize=11)
    b.legend(frameon=False, fontsize=7.5, loc='upper center', bbox_to_anchor=(0.5, -0.16), labelcolor=INK2, ncol=2)
    b.set_ylim(-0.5, 0.05)
    style(b)
    fig.tight_layout()
    fig.savefig(OUT / 'step1_final_r2.png', dpi=150)


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    R = load()
    metrics_figure(R)
    scale_figure(R)
    r2_figure(R)
    print(OUT)

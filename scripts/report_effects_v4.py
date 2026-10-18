"""Report on the v4 effect tables (lib/effects_v4/): per-source summary, scale plot, per-gene coverage. Read only.

- Per source: rows, measured panel genes, median efficiency, median and 90th-percentile n_de_bh, checks (<source>_build.json).
- results/qc/effects_v4_scale.png: |lfc| and se per source within 4 control-expression bins (control mean CPM < 1, 1-10,
  10-100, >= 100, from each source's baseline table), to spot sources that are off in scale.
- results/qc/effects_v4_coverage.csv: per panel gene, the number of sources (files) and backgrounds that measured it;
  genes measured by 0 or 1 source are listed in the report.
Writes results/qc/effects_v4_report.json.
"""
import json
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
DIR = V / 'lib/effects_v4'
OUT = V / 'results/qc'
BINS = [(0, 1, '< 1 CPM'), (1, 10, '1–10'), (10, 100, '10–100'), (100, np.inf, '≥ 100')]
SAMPLE = 400_000


def main():
    files = sorted(p for p in DIR.glob('*.h5ad'))
    rng = np.random.default_rng(2026); summary, coverage, scale = [], {}, []
    background = {}
    for p in files:
        name = p.stem; a = ad.read_h5ad(p)
        build = json.loads((DIR / f'{name}_build.json').read_text()) if (DIR / f'{name}_build.json').exists() else {}
        base = pd.read_csv(DIR / f'{name}_baseline.csv', index_col=0)
        lfc, se = a.layers['lfc'], a.layers['se']
        measured = np.isfinite(lfc).any(0); coverage[name] = measured; background[name] = str(a.obs['background'].iloc[0])
        summary.append(dict(source=name, background=background[name], rows=a.n_obs, measured_genes=int(measured.sum()),
                            median_efficiency=float(np.nanmedian(a.obs.efficiency)) if a.obs.efficiency.notna().any() else None,
                            n_de_bh_median=float(a.obs.n_de_bh.median()), n_de_bh_p90=float(a.obs.n_de_bh.quantile(0.9)),
                            lfc_type=str(a.obs.lfc_type.iloc[0]), se_type=str(a.obs.se_type.iloc[0]), checks=build.get('checks', {})))
        cpm = base['control_mean_cpm'].reindex(a.var_names).to_numpy(float)
        for lo, hi, label in BINS:
            cols = np.flatnonzero(measured & (cpm >= lo) & (cpm < hi))
            if not len(cols): continue
            vals_l = np.abs(lfc[:, cols]).ravel(); vals_s = se[:, cols].ravel(); ok = np.isfinite(vals_l)
            idx = np.flatnonzero(ok); idx = rng.choice(idx, min(SAMPLE, len(idx)), replace=False)
            scale.append(dict(source=name, bin=label, lfc=vals_l[idx], se=vals_s[idx], genes=len(cols)))
        del a
    S = pd.DataFrame(summary)
    # coverage
    panel = pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).to_numpy()
    cov = pd.DataFrame(coverage, index=panel)
    by_bg = pd.DataFrame({bg: cov[[n for n in cov.columns if background[n] == bg]].any(axis=1) for bg in sorted(set(background.values()))})
    cov_out = pd.DataFrame(dict(n_sources=cov.sum(1), n_backgrounds=by_bg.sum(1)))
    cov_out.to_csv(OUT / 'effects_v4_coverage.csv')
    # plot
    names = [n for n in S.source]
    fig, axes = plt.subplots(2, 1, figsize=(max(10, 1.1 * len(names) * len(BINS) / 2), 8), sharex=True)
    colors = ['#c7dcef', '#7fb2dd', '#3a7bbf', '#1b4b80']
    for ax, key, ylabel in ((axes[0], 'lfc', '|lfc|'), (axes[1], 'se', 'se')):
        pos, data, cols = [], [], []
        for i, n in enumerate(names):
            for k, (_, _, label) in enumerate(BINS):
                d = [x for x in scale if x['source'] == n and x['bin'] == label]
                if not d: continue
                pos.append(i * (len(BINS) + 1) + k); data.append(d[0][key]); cols.append(colors[k])
        bp = ax.boxplot(data, positions=pos, widths=0.8, showfliers=False, patch_artist=True)
        for patch, c in zip(bp['boxes'], cols): patch.set_facecolor(c)
        ax.set_yscale('log'); ax.set_ylabel(ylabel + ' (log scale)')
    axes[1].set_xticks([i * (len(BINS) + 1) + 1.5 for i in range(len(names))]); axes[1].set_xticklabels(names, rotation=30, ha='right')
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colors]
    axes[0].legend(handles, [b[2] for b in BINS], title='control mean CPM', ncol=4, frameon=False, loc='upper left')
    axes[0].set_title('effects v4: |lfc| and se per source within control-expression bins (boxes: quartiles; whiskers 1.5 IQR)')
    fig.tight_layout(); fig.savefig(OUT / 'effects_v4_scale.png', dpi=150); plt.close(fig)
    scale_table = pd.DataFrame([dict(source=x['source'], bin=x['bin'], genes=x['genes'], median_abs_lfc=float(np.median(x['lfc'])),
                                     median_se=float(np.median(x['se']))) for x in scale])
    report = dict(sources=S.drop(columns='checks').to_dict('records'), checks={r['source']: r['checks'] for r in summary},
                  scale=scale_table.to_dict('records'),
                  coverage=dict(panel_genes=len(panel), by_sources=cov_out.n_sources.value_counts().sort_index().to_dict(),
                                by_backgrounds=cov_out.n_backgrounds.value_counts().sort_index().to_dict(),
                                genes_0_sources=cov_out.index[cov_out.n_sources == 0].tolist(),
                                genes_1_source=cov_out.index[cov_out.n_sources == 1].tolist()))
    (OUT / 'effects_v4_report.json').write_text(json.dumps(report, indent=1, default=str) + '\n')
    pd.set_option('display.width', 220)
    print(S.drop(columns='checks').to_string(index=False))
    print('\nscale (medians):'); print(scale_table.pivot(index='source', columns='bin', values='median_abs_lfc').round(3).to_string())
    print(scale_table.pivot(index='source', columns='bin', values='median_se').round(3).to_string())
    print('\ncoverage, panel genes by number of sources measuring them:', report['coverage']['by_sources'])
    print('by number of backgrounds:', report['coverage']['by_backgrounds'])
    print(f"measured by 0 sources: {len(report['coverage']['genes_0_sources'])}; by 1 source: {len(report['coverage']['genes_1_source'])}")


if __name__ == '__main__':
    main()

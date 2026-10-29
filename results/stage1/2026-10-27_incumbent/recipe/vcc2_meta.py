"""Random-effects meta-analysis of the public panel sources: precision-weighted average (Model 2).

Inputs: lib/sources/<NAME>/effects_panel.h5ad from scripts/build_panel_sources.py (575 targets x 18,533 genes,
layers lfc / se, NaN = not available). One source per background: the T cells enter as CD4_Stim48hr only.

  fit      per-gene tau2: DerSimonian-Laird, pooled over every target with >= 2 backgrounds -> model dir
  predict  one delta table per line, in the folder format vcc2_generate.py reads:
             re.npy, re_se.npy         weight = 1/(SE^2 + tau2); prediction = sum(weight x lfc) / sum(weight);
                                       SE = sqrt(1 / sum(weight))
             n_backgrounds.npy         backgrounds behind each value; 0 = no information (delta stored as NaN)
             ontarget_lfc.npy          mean of the sources' ontarget_lfc: the knocked-down gene, applied separately
  exam     the held-out background's kept targets that at least one other background also kept
  explain  one target x gene: each background's lfc, SE, weight and share, then the prediction

The knocked-down gene's own column never enters tau2 or the averages. No shrinkage: the average is the prediction.

Example: python vcc2_meta.py fit --sources-dir ~/vcc/lib/sources --genes gene_names.csv --holdout H1 --out models/exam
         python vcc2_meta.py predict --model models/exam --targets exam_targets.csv --lines H1 --out exam
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from vcc2_common import read_names, save_json

SCHEMA = 'vcc_2026_2_meta_v1'
SOURCES = ['HCT116', 'HEK293T', 'K562_gwps', 'CD4_Stim48hr', 'H1', 'RPE1', 'HepG2', 'Jurkat']   # step 1: Stim48hr is the one T-cell vote
MODELS = ('re',)


def background_of(root, name):
    a = ad.read_h5ad(Path(root) / name / 'effects_panel.h5ad', backed='r')
    try: return str(a.uns['background']), a.obs.copy()
    finally: a.file.close()


def included(root, names, holdout=None):
    """Sources that vote, in order; refuses two sources from one background."""
    chosen, seen = [], {}
    for name in names:
        background, _ = background_of(root, name)
        if background == holdout: continue
        if background in seen: raise ValueError(f'{name} and {seen[background]} are both {background}: one vote per background')
        seen[background] = name; chosen.append((name, background))
    if not chosen: raise ValueError('no source backgrounds remain')
    return chosen


class Panel:
    """Included sources stacked as sources x targets x genes (lfc, se), plus each source's on-target LFC."""
    def __init__(self, root, sources, genes):
        self.names = [n for n, _ in sources]; self.backgrounds = [b for _, b in sources]
        lfc, se, own = [], [], []
        for name in self.names:
            a = ad.read_h5ad(Path(root) / name / 'effects_panel.h5ad')
            if not np.array_equal(a.var_names.astype(str), genes): raise ValueError(f'{name}: gene order differs from gene_names.csv')
            targets = a.obs_names.astype(str).to_numpy()
            if lfc and not np.array_equal(targets, self.targets): raise ValueError(f'{name}: target order differs between panels')
            self.targets = targets
            kept = (a.obs['status'] == 'kept').to_numpy()
            lfc.append(np.asarray(a.layers['lfc'], np.float32)); se.append(np.asarray(a.layers['se'], np.float32))
            own.append(np.where(kept, a.obs['ontarget_lfc'].to_numpy(float), np.nan))
        self.lfc, self.se, self.ontarget = np.stack(lfc), np.stack(se), np.stack(own)
        gp = {g: i for i, g in enumerate(genes)}
        self.own_col = np.array([gp.get(t, -1) for t in self.targets])
        self.row = {t: i for i, t in enumerate(self.targets)}


def tau2_dl(lfc, se, own_col):
    """Per-gene between-background variance, DerSimonian-Laird pooled over targets.

    For each target with k >= 2 backgrounds: Q = sum w (y - mu_FE)^2 with w = 1/SE^2 (noise alone predicts k - 1)
    and C = sum w - sum w^2 / sum w.  tau2 = max(0, (sum Q - sum (k - 1)) / sum C), floored at 0.
    Returns tau2 (NaN where no target ever had two backgrounds) and the pooled degrees of freedom.
    """
    _, T, G = lfc.shape
    Q, df, C = np.zeros(G), np.zeros(G), np.zeros(G)
    for t in range(T):
        y, v = lfc[:, t].astype(np.float64), se[:, t].astype(np.float64) ** 2
        ok = np.isfinite(y) & np.isfinite(v) & (v > 0)
        if own_col[t] >= 0: ok[:, own_col[t]] = False
        k = ok.sum(0); use = k >= 2
        if not use.any(): continue
        w = np.where(ok, 1 / np.where(ok, v, 1), 0); y = np.where(ok, y, 0)
        sw = w.sum(0); mu = np.divide((w * y).sum(0), sw, out=np.zeros(G), where=sw > 0)
        Q[use] += (w * (y - mu) ** 2).sum(0)[use]; df[use] += k[use] - 1
        C[use] += (sw - np.divide((w ** 2).sum(0), sw, out=np.zeros(G), where=sw > 0))[use]
    tau2 = np.full(G, np.nan); est = C > 0
    tau2[est] = np.maximum(0, (Q[est] - df[est]) / C[est])
    return tau2, df


def combine(lfc, se, tau2):
    """Sources x rows x genes -> (re, re_se, k) per row x gene. k == 0 gives NaN (no information)."""
    y = lfc.astype(np.float64); total = se.astype(np.float64) ** 2 + tau2   # SE^2 + tau2
    ok = np.isfinite(y) & np.isfinite(total)
    k = ok.sum(0); y = np.where(ok, y, 0)
    w = np.where(ok, 1 / np.where(ok, total, 1), 0); sw = w.sum(0)
    with np.errstate(divide='ignore', invalid='ignore'):
        out = ((w * y).sum(0) / sw, np.sqrt(1 / sw))
    return tuple(np.where(k > 0, x, np.nan).astype(np.float32) for x in out) + (k.astype(np.uint8),)


def fit(sources_dir, genes_csv, output, holdout=None, names=SOURCES):
    genes = read_names(genes_csv, 'gene_name'); sources = included(sources_dir, names, holdout)
    panel = Panel(sources_dir, sources, genes)
    tau2, df = tau2_dl(panel.lfc, panel.se, panel.own_col)
    est = np.isfinite(tau2); fallback = float(np.median(tau2[est]))
    measured_in = np.isfinite(panel.lfc).any(1).sum(0)                  # backgrounds measuring each gene at all
    out = Path(output); out.mkdir(parents=True, exist_ok=True)
    np.save(out / 'genes.npy', genes); np.save(out / 'tau2.npy', tau2.astype(np.float32)); np.save(out / 'tau2_df.npy', df)
    q = np.quantile(tau2[est], [0.1, 0.25, 0.5, 0.75, 0.9])
    model = dict(schema=SCHEMA, sources_dir=str(Path(sources_dir).resolve()), genes_csv=str(Path(genes_csv).resolve()),
                 sources=[n for n, _ in sources], backgrounds=[b for _, b in sources], heldout_background=holdout,
                 tau2_estimator='DerSimonian-Laird, pooled over targets with >= 2 backgrounds, own column excluded, floored at 0',
                 tau2_fallback=fallback, tau2_fallback_rule='median tau2 for genes no target measured in 2 backgrounds',
                 genes_estimated=int(est.sum()), genes_fallback=int((~est & (measured_in > 0)).sum()),
                 genes_no_data=int((measured_in == 0).sum()), genes_tau2_zero=int((tau2[est] == 0).sum()),
                 tau2_quantiles=dict(zip(['10%', '25%', '50%', '75%', '90%'], q.round(5).tolist())))
    save_json(out / 'model.json', model)
    print(f'fit: {len(sources)} backgrounds {model["backgrounds"]} (held out: {holdout})')
    print(f'  tau2 estimated for {model["genes_estimated"]:,} genes ({model["genes_tau2_zero"]:,} at 0), '
          f'median {fallback:.4f}; {model["genes_fallback"]:,} genes measured in one background use the median; '
          f'{model["genes_no_data"]:,} genes have no data anywhere')
    return model


def load_model(model_dir):
    m = Path(model_dir); model = json.loads((m / 'model.json').read_text())
    if model.get('schema') != SCHEMA: raise ValueError('incompatible model; fit with this notebook')
    genes = np.load(m / 'genes.npy').astype(str); tau2 = np.load(m / 'tau2.npy').astype(np.float64)
    used = np.where(np.isfinite(tau2), tau2, model['tau2_fallback'])
    panel = Panel(model['sources_dir'], list(zip(model['sources'], model['backgrounds'])), genes)
    return model, genes, used, panel


def predict(model_dir, targets_csv, lines, output, chunk=64):
    model, genes, tau2, panel = load_model(model_dir)
    targets = read_names(targets_csv, 'target_gene')
    if not len(targets): raise ValueError('empty perturbation panel')
    out = Path(output)
    for line in lines:   # an earlier meta prediction is replaced (it is deterministic); anything else is kept
        d = out / line; earlier = d / 'prediction_report.json'
        if d.exists() and any(d.iterdir()) and not (earlier.exists() and json.loads(earlier.read_text()).get('schema') == SCHEMA):
            raise ValueError(f'choose an empty prediction directory: {d}')
        if d.exists():
            for old in [*d.glob('*.npy'), d / 'prediction_report.json']: old.unlink(missing_ok=True)
    T, G = len(targets), len(genes)
    arrays = {name: np.full((T, G), np.nan, np.float32) for name in ['re', 're_se']}
    count = np.zeros((T, G), np.uint8); ontarget = np.full(T, np.nan, np.float32)
    rows = np.array([panel.row.get(t, -1) for t in targets]); missing = targets[rows < 0].tolist()
    gp = {g: i for i, g in enumerate(genes)}
    for start in range(0, T, chunk):
        idx = np.arange(start, min(start + chunk, T)); idx = idx[rows[idx] >= 0]
        if not len(idx): continue
        lfc, se = panel.lfc[:, rows[idx]].copy(), panel.se[:, rows[idx]]
        for j, i in enumerate(idx):
            if targets[i] in gp: lfc[:, j, gp[targets[i]]] = np.nan   # step 5: own gene is not averaged
        arrays['re'][idx], arrays['re_se'][idx], k = combine(lfc, se, tau2)
        count[idx] = k
        own = panel.ontarget[:, rows[idx]]; n = np.isfinite(own).sum(0)
        ontarget[idx] = np.divide(np.nan_to_num(own).sum(0), n, out=np.full(len(idx), np.nan), where=n > 0)
    voters = [[b for b, kept in zip(panel.backgrounds, np.isfinite(panel.ontarget[:, r])) if kept] if r >= 0 else [] for r in rows]
    info = (count > 0).sum(1)
    report = dict(schema=SCHEMA, models=list(MODELS), heldout_background=model['heldout_background'],
                  source_backgrounds=model['backgrounds'], source_names=model['sources'], model_dir=str(Path(model_dir).resolve()),
                  targets=str(Path(targets_csv).resolve()), weighting='1 / (SE^2 + tau2)', shrinkage='none', ontarget='mean of the backgrounds\' ontarget_lfc, kept rows only',
                  no_information='NaN in re (n_backgrounds 0); the generator predicts no change there',
                  not_in_panel=missing, no_source_targets=[t for t, v in zip(targets.tolist(), voters) if not v],
                  coverage=[dict(target=t, backgrounds=v, genes_with_information=int(n)) for t, v, n in zip(targets.tolist(), voters, info)])
    for line in lines:
        d = out / line; d.mkdir(parents=True, exist_ok=True)
        np.save(d / 'genes.npy', genes); np.save(d / 'targets.npy', targets)
        for name, value in arrays.items(): np.save(d / f'{name}.npy', value)
        np.save(d / 'n_backgrounds.npy', count); np.save(d / 'ontarget_lfc.npy', ontarget)
        save_json(d / 'prediction_report.json', dict(report, line=line))
    by_k = pd.Series([len(v) for v in voters]).value_counts().sort_index()
    print(f'predict: {T} targets x {G:,} genes -> {", ".join(str(out / l) for l in lines)}')
    print('  targets by number of backgrounds: ' + ', '.join(f'{k}: {n}' for k, n in by_k.items()))
    print(f'  no source at all (all genes "no information"): {report["no_source_targets"] or "none"}')
    if missing: print(f'  not in the panel: {missing}')
    return report


def exam(sources_dir, holdout, output, names=SOURCES):
    status = {}
    for name in names:
        background, obs = background_of(sources_dir, name)
        status.setdefault(background, []).append(obs['status'].eq('kept'))
    if holdout not in status: raise ValueError(f'no source for {holdout}')
    held = pd.concat(status.pop(holdout), axis=1).any(axis=1)
    others = pd.concat([pd.concat(v, axis=1).any(axis=1) for v in status.values()], axis=1).sum(axis=1)
    targets = held.index[held & (others >= 1)].tolist()
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({'target_gene': targets}).to_csv(output, index=False)
    print(f'exam: {len(targets)} {holdout} targets also kept in >= 1 other background -> {output}')
    return targets


def explain(model_dir, target, gene):
    """The worked example for one target x gene, from the real panel."""
    if target == gene: raise ValueError('the knocked-down gene is not averaged; see ontarget_lfc')
    model, genes, tau2, panel = load_model(model_dir)
    r, g = panel.row[target], int(np.flatnonzero(genes == gene)[0])
    t2 = tau2[g]; y = panel.lfc[:, r, g].astype(float); se = panel.se[:, r, g].astype(float)
    frame = pd.DataFrame(dict(background=panel.backgrounds, lfc=y, se=se)).set_index('background')
    frame['se2_plus_tau2'] = se ** 2 + t2; frame['weight'] = 1 / frame.se2_plus_tau2
    frame['share'] = frame.weight / frame.weight.sum()
    re, re_se, k = (x[0, 0] for x in combine(frame.lfc.to_numpy()[:, None, None],
                                                                 frame.se.to_numpy()[:, None, None], np.array([t2])))
    summary = dict(target=target, gene=gene, tau2=float(t2), tau2_estimated=bool(np.isfinite(np.load(Path(model_dir) / 'tau2.npy')[g])),
                   n_backgrounds=int(k), re=float(re), re_se=float(re_se))
    return frame, summary


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    f = sub.add_parser('fit')
    for k in ['sources-dir', 'genes', 'out']: f.add_argument('--' + k, required=True)
    f.add_argument('--holdout'); f.add_argument('--sources', nargs='+', default=SOURCES)
    pr = sub.add_parser('predict')
    for k in ['model', 'targets', 'out']: pr.add_argument('--' + k, required=True)
    pr.add_argument('--lines', nargs='+', required=True); pr.add_argument('--chunk', type=int, default=64)
    e = sub.add_parser('exam')
    for k in ['sources-dir', 'holdout', 'out']: e.add_argument('--' + k, required=True)
    e.add_argument('--sources', nargs='+', default=SOURCES)
    x = sub.add_parser('explain')
    for k in ['model', 'target', 'gene']: x.add_argument('--' + k, required=True)
    a = p.parse_args()
    if a.cmd == 'fit': fit(a.sources_dir, a.genes, a.out, a.holdout, a.sources)
    elif a.cmd == 'predict': predict(a.model, a.targets, a.lines, a.out, a.chunk)
    elif a.cmd == 'exam': exam(a.sources_dir, a.holdout, a.out, a.sources)
    else:
        frame, summary = explain(a.model, a.target, a.gene)
        print(frame.to_string(float_format='%.4f')); print(json.dumps(summary, indent=2))
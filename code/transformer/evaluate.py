"""Effect-space evaluation and predictions for a trained vcc_2026_4 run.

  PYTHONPATH=code python -m transformer.evaluate effect --run lib/vcc_2026_4/runs/dev
  PYTHONPATH=code python -m transformer.evaluate h1-pipeline --run ... --taus 0.3 0.5 0.7 0.9
  PYTHONPATH=code python -m transformer.evaluate predict --run ... --lines A B C --out <dir>
Effect space, per knockdown then the median: Pearson of lfc over measured genes, Pearson over truly DE genes (>= 3),
sign accuracy on truly DE genes, AUPRC of p_de for truly DE genes; model next to the baseline (precision-weighted mean of
the evidence lfc; score |z|). Truly DE = BH 5% in the target line over its measured genes (own gene excluded).
Held-out genes: the run's validation rows (every target line). H1: all 300 H1 knockdowns, zero-shot (validation mode: H1
is neither evidence nor target).
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
from anndata.io import read_elem
from sklearn.metrics import average_precision_score

from . import config as C
from .data import LINE, LINES, Bank, de_full, kd_keys
from .model import VCCTransformer

V = Path.home() / 'vcc'


def load_run(run, device='cuda', ckpt='best.pt'):
    run = Path(run); cfg = json.loads((run / 'config.json').read_text())
    bank = Bank(cfg, device=device)
    assert np.array_equal(np.load(run / 'gene_cols.npy'), bank.gcols), 'gene set differs from the run'
    model = VCCTransformer(cfg, len(bank.genes)).to(bank.device)
    state = torch.load(run / ckpt, map_location=bank.device); model.load_state_dict(state['model']); model.eval()
    return cfg, bank, model


@torch.no_grad()
def predict(model, bank, kd, line, trow, bs=32, bf16=True):
    """delta_hat, p_de, baseline pw_mean, baseline z: (N, G) float32 numpy."""
    if bank.device.type == 'cpu': bs = min(bs, 8)
    out = {k: [] for k in ('delta', 'p', 'pw_mean', 'pw_z')}
    for s in range(0, len(kd), bs):
        b = bank.batch(kd=kd[s:s + bs], line=line[s:s + bs], trow=trow[s:s + bs])
        with torch.autocast('cuda', dtype=torch.bfloat16, enabled=bf16 and bank.device.type == 'cuda'):
            o = model(b)
        out['delta'].append(o['delta'].float().cpu()); out['p'].append(torch.sigmoid(o['p_logit'].float()).cpu())
        out['pw_mean'].append(b['pw_mean'].cpu()); out['pw_z'].append(b['pw_z'].cpu())
    return {k: torch.cat(v).numpy() for k, v in out.items()}


def kd_metrics(true, de, measured, pred, score):
    m = measured; t = true[m]; d = de[m]; p = pred[m]; s = score[m]
    r = dict(n_genes=int(m.sum()), n_de=int(d.sum()))
    r['r_all'] = float(np.corrcoef(t, p)[0, 1]) if m.sum() > 2 and np.std(p) > 0 and np.std(t) > 0 else np.nan
    r['r_de'] = float(np.corrcoef(t[d], p[d])[0, 1]) if d.sum() >= 3 and np.std(p[d]) > 0 else np.nan
    r['sign_de'] = float((np.sign(p[d]) == np.sign(t[d])).mean()) if d.any() else np.nan
    r['auprc'] = float(average_precision_score(d, s)) if d.any() and not d.all() else np.nan
    r['auprc_random'] = float(d.mean()) if m.any() else np.nan
    return r


def table(rows_true, preds, info):
    """Per-knockdown metrics for model and baseline; rows_true: list of (lfc, de, measured) arrays."""
    recs = []
    for i, (t, d, m) in enumerate(rows_true):
        base = dict(info[i])
        recs.append(dict(base, method='model', **kd_metrics(t, d, m, preds['delta'][i], preds['p'][i])))
        recs.append(dict(base, method='baseline', **kd_metrics(t, d, m, preds['pw_mean'][i], np.abs(preds['pw_z'][i]))))
    return pd.DataFrame(recs)


def summary(df, strong=50):
    cols = ['r_all', 'r_de', 'sign_de', 'auprc', 'auprc_random']
    out = {}
    for meth, g in df.groupby('method'):
        out[(meth, 'all')] = g[cols].median()
        s = g[g.n_de_bh >= strong]
        if len(s): out[(meth, f'n_de_bh>={strong}')] = s[cols].median()
    t = pd.DataFrame(out).T
    t['knockdowns'] = [len(df[df.method == m]) if k == 'all' else int((df[df.method == m].n_de_bh >= strong).sum())
                       for m, k in t.index]
    return t


def heldout_effect(model, bank, bs=32):
    rows = bank.val_rows
    kd, line, trow = (x[torch.as_tensor(rows, device=bank.device)] for x in bank._rows_t)
    preds = predict(model, bank, kd, line, trow, bs)
    rows_true, info = [], []
    for s in range(0, len(rows), bs):
        b = bank.batch(rows=rows[s:s + bs])
        lt = b['lfc_t'].cpu().numpy(); de = b['de'].cpu().numpy(); lm = b['loss_mask'].cpu().numpy()
        for j in range(len(lt)):
            r = rows[s + j]; ln = LINES[bank.row_line[r]]
            rows_true.append((lt[j], de[j], lm[j]))
            info.append(dict(line=ln, target=bank.kd_vocab[bank.row_kd[r]], n_de_bh=int(de[j].sum())))
    return table(rows_true, preds, info)


def h1_truth(bank):
    """H1's own table on the run's genes (read from disk; never part of a validation-mode bank)."""
    E = Path(bank.cfg['paths']['effects'])
    with h5py.File(E / 'arc_h1.h5ad', 'r') as f:
        obs = read_elem(f['obs']); L = f['layers/lfc'][:][:, bank.gcols]; S = f['layers/se'][:][:, bank.gcols]
    keys = kd_keys(obs); col = {g: i for i, g in enumerate(bank.genes)}
    kcols = np.array([col.get(s, -1) for s in obs['target_gene'].astype(str)])
    de = de_full(E, 'arc_h1', bank.cfg['paths']['cache'], kcols)[:, bank.gcols]
    own = bank.gcols[None, :] == kcols[:, None]
    measured = np.isfinite(L) & np.isfinite(S) & ~own
    return obs, keys, np.nan_to_num(L), de & measured, measured


def h1_effect(model, bank, bs=32):
    obs, keys, L, de, measured = h1_truth(bank)
    kd, line, trow = bank.rows_for('H1', keys)
    preds = predict(model, bank, kd, line, trow, bs)
    info = [dict(line='H1', target=k, target_gene=g, n_de_bh=int(n)) for k, g, n in
            zip(keys, obs['target_gene'].astype(str), de.sum(1))]
    return table(list(zip(L, de, measured)), preds, info), preds, obs


SCHEMA = 'vcc_2026_3_ridge_v1'          # the prediction-folder schema vcc3_generate.py checks


def write_pred_dir(out, targets, genes, gcols, delta, p, ontarget, taus, sources, heldout, extra=None, topn=None):
    """A prediction folder for vcc3_generate.py: one full-panel array per calling threshold (p_de > tau keeps delta_hat,
    every other gene and every gene outside the run's set is 0 = no change), plus ontarget_lfc and the report."""
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    G = len(genes); models = []
    np.save(out / 'genes.npy', np.asarray(genes).astype(str)); np.save(out / 'targets.npy', np.asarray(targets).astype(str))
    np.save(out / 'ontarget_lfc.npy', np.asarray(ontarget, np.float32))
    for tau in taus:
        D = np.zeros((len(targets), G), np.float32)
        D[:, gcols] = np.where(p > tau, delta, 0.0)
        name = f'tf_tau{tau:g}'; np.save(out / f'{name}.npy', D); models.append(name)
    for n in topn or []:                                   # per knockdown, the n genes with the highest p_de
        D = np.zeros((len(targets), G), np.float32); k = min(n, p.shape[1])
        idx = np.argpartition(-p, k - 1, axis=1)[:, :k]; rows = np.arange(len(targets))[:, None]
        D[rows, gcols[idx]] = delta[rows, idx]
        name = f'tf_top{n}'; np.save(out / f'{name}.npy', D); models.append(name)
    for name, arr in (extra or {}).items():
        D = np.zeros((len(targets), G), np.float32); D[:, gcols] = arr; np.save(out / f'{name}.npy', D); models.append(name)
    report = dict(schema=SCHEMA, heldout_background=heldout, source_backgrounds=sources, models=models,
                  note='vcc_2026_4 transformer; calling p_de > tau, uncalled genes and genes outside the run set = no change')
    (out / 'prediction_report.json').write_text(json.dumps(report, indent=1))
    return models


def h1_pipeline(cfg, bank, model, run, taus, n_targets=75, floor=False, topn=None, tag='', baseline=True, alphas=(1.0,)):
    """Zero-shot H1 predictions -> existing generator -> official scorer on the unfiltered H1 cells, scored on the same
    75-target reference (lib/vcc_2026_3/val_final) as the earlier H1 results."""
    import subprocess
    if bank.mode != 'validation': raise ValueError('the H1 pipeline needs a validation-mode run (H1 unseen)')
    obs, keys, L, de, measured = h1_truth(bank)
    kd, line, trow = bank.rows_for('H1', keys)
    preds = predict(model, bank, kd, line, trow)
    eff = torch.cat([bank.batch(kd=kd[s:s + 64], line=line[s:s + 64], trow=trow[s:s + 64])['kd_feats'][:, -2:].cpu()
                     for s in range(0, len(kd), 64)]).numpy()
    ontarget = np.where(eff[:, 1] > 0, eff[:, 0], np.nan)
    pdir = Path(run) / f'h1_pred{tag}'; pdir.mkdir(parents=True, exist_ok=True)
    np.savez(pdir / 'raw_predictions.npz', delta=preds['delta'], p=preds['p'], pw_mean=preds['pw_mean'], gene_cols=bank.gcols)
    sources = [l for l in bank.loaded if l != 'H1']
    models = write_pred_dir(pdir, obs['target_gene'].astype(str).to_numpy(), bank.genes, bank.gcols, preds['delta'],
                            preds['p'], ontarget, taus, sources, 'H1',
                            extra=dict(baseline_pw=preds['pw_mean']) if baseline else None, topn=topn)
    variants = [dict(method=m, call='z', z=0.0, min_effect=0.0, alpha=al, detect_floor=floor,
                     label=m + (f'_a{al:g}' if al != 1 else '') + ('_floor' if floor else '')) for m in models for al in alphas]
    vpath = Path(run) / f'h1_variants{tag}.json'; vpath.write_text(json.dumps(variants, indent=1))
    cmd = ['python', 'vcc3_generate.py', 'benchmark', '--cells', str(V / 'lib/h1/H1_full_2025.h5ad'),
           '--genes', str(V / 'data/gene_names.csv'), '--targets', str(V / 'lib/vcc_2026_3/validation/targets.csv'),
           '--predictions', str(pdir), '--eval-profile', 'vcc2026', '--n-targets', str(n_targets),
           '--reference', str(V / 'lib/vcc_2026_3/val_final'), '--variants', str(vpath), '--out', str(Path(run) / f'h1_bench{tag}')]
    import gc, sys
    cmd[0] = sys.executable
    # The generator + scorer run for about an hour on CPU: give the GPU back first (other runs share it).
    for t in list(bank.__dict__):
        if isinstance(getattr(bank, t), (torch.Tensor, list, dict)): setattr(bank, t, None)
    model.cpu(); gc.collect(); torch.cuda.empty_cache()
    print(' '.join(cmd), flush=True)
    import os
    env = dict(os.environ, PATH=str(Path(sys.executable).parent) + os.pathsep + os.environ.get('PATH', ''))   # cell-eval2
    subprocess.run(cmd, cwd=V / 'code/vcc_2026_3c', check=True, env=env)
    scaled = pd.read_csv(Path(run) / f'h1_bench{tag}/scaled_scores.csv', index_col=0)
    print(scaled.round(4).to_string())
    return scaled


def challenge_rows(bank, line):
    """Rows for our 300 targets (data/pert_counts.csv order) in challenge line A, B or C."""
    tg = pd.read_csv(V / 'data/pert_counts.csv')['target_gene'].astype(str).tolist()
    col = {g: i for i, g in enumerate(bank.genes)}
    keys = [bank.panel_ids[col[g]] if bank.panel_ids[col[g]].startswith('ENSG') else f'symbol:{g}' for g in tg]
    return tg, bank.rows_for(line, keys)


def efficiency(bank, kd, line, trow):
    eff = torch.cat([bank.batch(kd=kd[s:s + 64], line=line[s:s + 64], trow=trow[s:s + 64])['kd_feats'][:, -2:].cpu()
                     for s in range(0, len(kd), 64)]).numpy()
    return np.where(eff[:, 1] > 0, eff[:, 0], np.nan)


def predict_lines(bank, model, out, taus, lines=('A', 'B', 'C'), pie_only=False, topn=None):
    """Prediction folders for the challenge lines: out/<line>/ (tf_tau<t>.npy or pie_tau<t>.npy, ontarget_lfc.npy).
    pie_only: PIE lfc_pred / p_de from lib/pie_tables instead of the transformer (the fallback submission)."""
    import h5py
    out = Path(out); summary = {}
    for line in lines:
        targets, (kd, ln, tr) = challenge_rows(bank, line)
        ontarget = efficiency(bank, kd, ln, tr)
        if pie_only:
            with h5py.File(Path(bank.cfg['paths']['pie_tables']) / f'pie_{line}.h5ad', 'r') as f:
                obs = read_elem(f['obs']); L = f['layers/lfc_pred'][:]; P = f['layers/p_de'][:]
            pos = {g: i for i, g in enumerate(obs['target_gene'].astype(str))}
            idx = [pos[t] for t in targets]
            delta = np.nan_to_num(L[idx][:, bank.gcols]); p = np.nan_to_num(P[idx][:, bank.gcols])
            prefix = 'pie'
        else:
            preds = predict(model, bank, kd, ln, tr); delta, p = preds['delta'], preds['p']; prefix = 'tf'
        models = write_pred_dir(out / line, targets, bank.genes, bank.gcols, delta, p, ontarget, taus,
                                [l for l in bank.loaded], 'none (final predictions)', topn=topn)
        np.savez(out / line / 'raw_predictions.npz', delta=delta, p=p, gene_cols=bank.gcols, targets=np.asarray(targets, dtype=str))
        if prefix == 'pie':
            for m in models:                                   # name the arrays after their source
                (out / line / f'{m}.npy').rename(out / line / f'{m.replace("tf_", "pie_")}.npy')
            models = [m.replace('tf_', 'pie_') for m in models]
            rep = json.loads((out / line / 'prediction_report.json').read_text()); rep['models'] = models
            (out / line / 'prediction_report.json').write_text(json.dumps(rep, indent=1))
        calls = {f'tau{t:g}': float((p > t).sum(1).mean()) for t in taus}
        summary[line] = dict(targets=len(targets), ontarget_known=int(np.isfinite(ontarget).sum()),
                             mean_calls_per_knockdown=calls, models=models)
        print(line, json.dumps(summary[line]), flush=True)
    (out / 'predict_summary.json').write_text(json.dumps(summary, indent=1))
    return summary


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('cmd', choices=['effect', 'h1-pipeline', 'predict-lines'])
    ap.add_argument('--run'); ap.add_argument('--ckpt', default='best.pt'); ap.add_argument('--out')
    ap.add_argument('--pie-only', action='store_true', help='predict-lines: PIE tables instead of a run (fallback)')
    ap.add_argument('--config', help='predict-lines --pie-only: bank config (final mode)')
    ap.add_argument('--taus', nargs='*', type=float, default=[0.3, 0.5, 0.7, 0.9])
    ap.add_argument('--topn', nargs='*', type=int, help='h1-pipeline: also top-N genes by p_de per knockdown')
    ap.add_argument('--floor', action='store_true', help='h1-pipeline: detection floor on called genes')
    ap.add_argument('--tag', default='', help='h1-pipeline: suffix for h1_pred / h1_bench folders')
    ap.add_argument('--no-baseline', action='store_true')
    ap.add_argument('--alphas', nargs='+', type=float, default=[1.0], help='h1-pipeline: generator alpha per arm')
    ap.add_argument('--device', default='cuda'); a = ap.parse_args()
    if a.cmd == 'predict-lines' and a.pie_only:
        cfg = C.load(a.config, mode='final', pie='on', gene_set='full'); bank = Bank(cfg, device='cpu')   # every panel gene PIE predicts
        predict_lines(bank, None, a.out, a.taus, pie_only=True, topn=a.topn); return
    cfg, bank, model = load_run(a.run, device=a.device, ckpt=a.ckpt); run = Path(a.run).resolve()   # absolute: the generator runs from its own folder
    if a.cmd == 'h1-pipeline':
        h1_pipeline(cfg, bank, model, run, a.taus, floor=a.floor, topn=a.topn, tag=a.tag, baseline=not a.no_baseline,
                    alphas=a.alphas); return
    if a.cmd == 'predict-lines':
        predict_lines(bank, model, a.out or run / 'pred_lines', a.taus, topn=a.topn); return
    if a.cmd == 'effect':
        held = heldout_effect(model, bank); held.to_csv(run / 'effect_heldout.csv', index=False)
        h1, _, _ = h1_effect(model, bank); h1.to_csv(run / 'effect_h1.csv', index=False)
        pd.set_option('display.width', 200)
        print('held-out genes (all target lines):'); print(summary(held).round(3).to_string())
        for ln, g in held[held.method == 'model'].groupby('line'):
            print(f'  {ln}: model r_all {g.r_all.median():.3f}, baseline {held[(held.method == "baseline") & (held.line == ln)].r_all.median():.3f} ({len(g)} knockdowns)')
        print('\nH1, all 300 knockdowns (zero-shot):'); print(summary(h1).round(3).to_string())
        res = dict(heldout={f'{m}|{k}': v.to_dict() for (m, k), v in summary(held).iterrows()},
                   h1={f'{m}|{k}': v.to_dict() for (m, k), v in summary(h1).iterrows()})
        (run / 'effect_summary.json').write_text(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()

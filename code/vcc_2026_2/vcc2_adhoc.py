"""Ad hoc diagnostics after the vcc_2026_2 submission (step 0 of the review). Nothing here changes the submission.

  arms  regenerate four arms on the exam's 25 H1 targets, reusing the exam's real file and real-side cache, and score
        each with cell-eval2, keeping both sides' DE tables (--write-degenes):
          null_shared       delta = 0 and no knockdown; every target is the one shared 400-cell pool (pure pool artefact)
          null_independent  delta = 0 and no knockdown; a fresh 400-cell draw per target
          submitted         the submitted settings (exam_b's selection)
          z2                z* = 2, min 0.1, alpha 1 (stage A's best arm on fid and Jaccard)
        The shared-pool arms rebuild the benchmark's cells exactly (same control split, pool and seeds), so 'submitted'
        must score as its exam_b arm did.

Functions used by the notebook cells:
  de_table(arm)        per-target predicted vs real DE genes (scorer's p_adj < 0.05, own gene excluded)
  generic_share(...)   how much of H1's real DE lands on genes DE in > 30% of knockdowns
  generic_transfer()   correlation of each source's average response with H1's
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import stats

import vcc2_generate as g
from vcc2_common import alignment, expression_gate, moments, seed_for, save_json

try:
    from anndata.io import read_elem
except ImportError:  # older anndata
    from anndata.experimental import read_elem

CONTROL = 'non-targeting'
SEED = 2026
ARMS = ['null_shared', 'null_independent', 'submitted', 'z2']
SOURCES = ['HCT116', 'HEK293T', 'K562_gwps', 'CD4_Stim48hr']


# ---------------------------------------------------------------- arms
def exam_line(cells, genes, seed=SEED):
    """The benchmark's generator controls and its scoring-half gate, rebuilt as stream_h1 does."""
    h = ad.read_h5ad(cells, backed='r')
    try:
        labels = h.obs['target_gene'].astype(str).values
        ci = np.random.default_rng(seed_for(seed, 'H1', 'controls')).permutation(np.flatnonzero(labels == CONTROL))
        train, score = ci[:len(ci) // 2], ci[len(ci) // 2:]
        small = h[np.sort(train)].to_memory()
        line = g.Line(g.aligned_X(small, genes), None, genes); del small
        line.gate = expression_gate(moments(h, np.sort(score), alignment(h.var_names.astype(str), genes))['mean_cpm'])
    finally:
        h.file.close()
    return line


def exam_targets(real):
    a = ad.read_h5ad(real, backed='r')
    try: counts = a.obs['target_gene'].astype(str).value_counts()
    finally: a.file.close()
    counts = counts.drop(CONTROL)
    return sorted(counts.index), counts.to_dict()


def write_arm(path, line, targets, effects, kd, picks, cfg, own):
    labels = [t for t in targets for _ in range(cfg.n_cells)] + [CONTROL] * line.n
    writer = g.StreamingH5AD(path, g.labelled_obs(labels, 'pred_'), line.genes)
    try:
        for i, t in enumerate(targets):
            writer.append(g.generate_one(line, effects[i], own[i], None if kd is None else kd[i], cfg,
                                         pick=picks[i], line_key='H1', target_key=t)[0])
        for start in range(0, line.n, 4096):
            writer.append(line.X[start:start + 4096].astype(np.int32))
    finally:
        writer.close()


def score(binary, pred, real, cache, out, threads=4):
    command = [binary, 'run', '-ap', pred, '-ar', real, '-o', out, '--cache-real', cache, '--pert-col', 'target_gene',
               '--control', CONTROL, '--profile', 'vcc2026', '--input-type', 'counts', '--preset', 'vcc2026',
               '--set', 'de.backend=pdex', '--set', f'num_threads={threads}', '--write-degenes']
    print('  $ ' + ' '.join(map(str, command)), flush=True)
    subprocess.run([str(x) for x in command], check=True)


def cmd_arms(a):
    genes = g.read_list(a.genes, 'gene_name'); exam = Path(a.exam); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    targets, real_counts = exam_targets(exam / 'real.h5ad')
    print(f'{len(targets)} exam targets; rebuilding the generator controls from {a.cells}', flush=True)
    line = exam_line(a.cells, genes)
    selected = json.loads((exam / 'selected_generator.json').read_text())['settings']
    settings = dict(submitted=dict(alpha=selected['alpha'], z=selected['z'], min_effect=selected['min_effect']),
                    z2=dict(alpha=1.0, z=2.0, min_effect=0.1))
    base = dict(pi=0.0, up_w=0.75, n_cells=a.n_cells, seed=SEED)
    pool_cfg = g.Cfg(argparse.Namespace(**base)); pool_cfg.n_cells = max([a.n_cells] + list(real_counts.values()))
    common = g.pick_cells(line, np.zeros(len(genes)), pool_cfg, np.random.default_rng(seed_for(SEED, 'H1', 'base_controls')))
    shared = [common[:a.n_cells]] * len(targets)
    cfg0 = g.Cfg(argparse.Namespace(**base))
    independent = [g.pick_cells(line, np.zeros(len(genes)), cfg0, np.random.default_rng(seed_for(SEED, 'H1', 'independent', t)))
                   for t in targets]
    gp = {x: i for i, x in enumerate(genes)}; own = [gp.get(t) for t in targets]
    zeros = np.zeros((len(targets), len(genes)), np.float32)
    D = g.load_deltas(a.predictions, targets, genes, 're'); se = g.load_deltas(a.predictions, targets, genes, 're_se', fill=np.inf)
    plan = dict(null_shared=(zeros, None, shared, [None] * len(targets), cfg0),
                null_independent=(zeros, None, independent, [None] * len(targets), cfg0))
    for name, s in settings.items():
        cfg = g.Cfg(argparse.Namespace(**base, **s))
        plan[name] = (g.calibrate(D, line, targets, cfg, se), g.kd_residual(g.prediction_ontarget(a.predictions, targets), cfg),
                      shared, own, cfg)
    for name in a.arms:
        effects, kd, picks, owns, cfg = plan[name]
        path = out / f'{name}.h5ad'
        write_arm(path, line, targets, effects, kd, picks, cfg, owns)
        score(a.cell_eval, path, exam / 'real.h5ad', exam / 'real_cache', out / name, a.num_threads)
        path.unlink()                                             # ~2 GB; the scores and DE tables are kept
        save_json(out / name / 'arm.json', dict(arm=name, settings={k: v for k, v in vars(cfg).items() if k in ['alpha', 'z', 'min_effect']},
                                                pool='independent' if name == 'null_independent' else 'shared', targets=targets))
        print(f'arm complete: {name}', flush=True)


# ---------------------------------------------------------------- analysis
def de_table(arm_dir, alpha=0.05):
    """Per target: predicted and real DE genes (p_adj < alpha, own gene excluded), overlap, precision, recall, and the
    chance levels: random genes would hit real DE genes at n_real / scored genes, and get the direction right half the time."""
    pred = pd.read_parquet(Path(arm_dir) / 'de_pred.parquet'); real = pd.read_parquet(Path(arm_dir) / 'de_real.parquet')
    rows = []
    for t, r in real.groupby('target'):
        p = pred[pred.target == t]
        P = set(p.feature[(p.p_adj < alpha) & (p.feature != t)]); R = set(r.feature[(r.p_adj < alpha) & (r.feature != t)])
        both = P & R
        same_sign = 0
        if both:
            lp = p.set_index('feature').log2_fold_change; lr = r.set_index('feature').log2_fold_change
            same_sign = int(sum(np.sign(lp[x]) == np.sign(lr[x]) for x in both))
        scored = r.feature.nunique()
        precision = len(both) / len(P) if P else np.nan
        rows.append(dict(target=t, n_pred=len(P), n_real=len(R), overlap=len(both), same_direction=same_sign,
                         precision=precision, chance_precision=len(R) / scored, enrichment=precision * scored / len(R) if R else np.nan,
                         direction_agreement=same_sign / len(both) if both else np.nan, recall=len(both) / len(R) if R else np.nan))
    return pd.DataFrame(rows).set_index('target')


def generic_from_de(de, alpha=0.05, share=0.3):
    """From a scorer DE table: genes DE in > `share` of targets, and the fraction of all DE calls on them."""
    calls = de[(de.p_adj < alpha) & (de.feature != de.target)]
    n_targets = de.target.nunique()
    freq = calls.groupby('feature').target.nunique() / n_targets
    generic = set(freq.index[freq > share])
    return dict(targets=n_targets, de_calls=len(calls), generic_genes=len(generic),
                calls_on_generic=int(calls.feature.isin(generic).sum()),
                share_of_calls_on_generic=float(calls.feature.isin(generic).mean()) if len(calls) else np.nan), generic, freq


def effect_table(path, min_cpm):
    with h5py.File(path, 'r') as f:
        obs, var = read_elem(f['obs']), read_elem(f['var'])
        keep = np.flatnonzero(var['ctrl_cpm'].to_numpy(float) > min_cpm)
        lfc = f['X'][:][:, keep]; se = f['layers/se'][:][:, keep]
    genes = var.index.astype(str).to_numpy()[keep]; targets = obs.index.astype(str).to_numpy()
    gp = {x: i for i, x in enumerate(genes)}
    for i, t in enumerate(targets):                                   # own gene is not a downstream response
        if t in gp: lfc[i, gp[t]] = np.nan
    return targets, genes, lfc, se


def generic_from_effects(path, alpha=0.05, share=0.3, min_cpm=5.0):
    """Same question on a pseudobulk effect table: z = lfc / SE, BH within each knockdown (like the scorer's per-pert FDR)."""
    targets, genes, lfc, se = effect_table(path, min_cpm)
    with np.errstate(invalid='ignore', divide='ignore'):
        p = 2 * stats.norm.sf(np.abs(lfc / se))
    de = np.zeros_like(lfc, bool)
    for i in range(len(targets)):
        ok = np.isfinite(p[i])
        if ok.any(): de[i, np.flatnonzero(ok)] = stats.false_discovery_control(p[i, ok]) < alpha
    freq = de.mean(0); generic = freq > share
    calls = int(de.sum())
    return dict(targets=len(targets), genes=len(genes), de_calls=calls, generic_genes=int(generic.sum()),
                calls_on_generic=int(de[:, generic].sum()),
                share_of_calls_on_generic=float(de[:, generic].sum() / calls) if calls else np.nan), pd.Series(freq, index=genes)


def average_response(path, min_cpm=1.0):
    targets, genes, lfc, _ = effect_table(path, min_cpm)
    return pd.Series(np.nanmean(lfc, 0), index=genes), len(targets)


def generic_transfer(filters_dir, h1_generic_freq, real_de=None, top=200):
    """Each source's average knockdown response vs H1's: Pearson r, and overlap of the strongest genes."""
    h1, n_h1 = average_response(Path(filters_dir) / 'H1.h5ad', min_cpm=5.0)
    references = {'H1 effect table (all kept knockdowns)': h1}
    if real_de is not None:
        references['H1 scorer DE (25 exam targets)'] = real_de[real_de.feature != real_de.target].groupby('feature').log2_fold_change.mean()
    sources = {name: average_response(Path(filters_dir) / f'{name}.h5ad')[0] for name in SOURCES}
    sources['mean of the 4 sources'] = pd.concat(sources.values(), axis=1).mean(axis=1)
    generic = set(h1_generic_freq.index[h1_generic_freq > 0.3])
    rows = []
    for ref_name, ref in references.items():
        for name, s in sources.items():
            common = ref.index.intersection(s.index); x, y = s[common], ref[common]
            ok = np.isfinite(x) & np.isfinite(y); x, y = x[ok], y[ok]
            top_ref = set(y.abs().nlargest(top).index); top_src = set(x.abs().nlargest(top).index)
            src_on_generic = set(x.abs().nlargest(len(generic)).index) & generic if generic else set()
            rows.append(dict(reference=ref_name, source=name, genes=len(x), pearson=float(np.corrcoef(x, y)[0, 1]),
                             spearman=float(stats.spearmanr(x, y).statistic), **{f'top{top}_overlap': len(top_ref & top_src) / top},
                             same_sign_on_top=float(np.mean(np.sign(x[list(top_ref & top_src)]) == np.sign(y[list(top_ref & top_src)])))
                             if top_ref & top_src else np.nan,
                             h1_generic_recovered=len(src_on_generic) / len(generic) if generic else np.nan,
                             source_mean_abs=float(x.abs().mean()), h1_mean_abs=float(y.abs().mean())))
    return pd.DataFrame(rows).set_index(['reference', 'source'])


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('arms')
    for key in ['cells', 'genes', 'exam', 'predictions', 'out']: s.add_argument('--' + key, required=True)
    s.add_argument('--arms', nargs='+', default=ARMS, choices=ARMS); s.add_argument('--n-cells', type=int, default=400)
    s.add_argument('--cell-eval', default='cell-eval2'); s.add_argument('--num-threads', type=int, default=4)
    a = p.parse_args()
    cmd_arms(a)

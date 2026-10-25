"""PCA programs + ridge transfer, on top of the Model 2 baseline (vcc_2026_3).

Inputs
  lib/filters/<NAME>.h5ad  every QC-passed knockdown of each source (X = lfc, layers/se, var ctrl_cpm): training rows
  lib/sources panels       only to refit each training destination's Model 2 tau2 without that destination
  a Model 2 prediction dir (vcc2_meta.py predict) for the line(s) being predicted, built with the same holdout

fit
  genes     one fixed gene set: measured (>= 1 CPM) in HCT116, HEK293T and K562, in gene_names.csv order
  programs  uncentred PCA of every knockdown of the program sources (Orion, K562, RPE1, HepG2, Jurkat) on that set.
            Each row's own knocked-down gene is set to 0 (its drop is handled separately by ontarget_lfc), so every
            row uses the same columns. The newer lines miss 6-10% of the set; those entries are filled iteratively
            from the low-rank fit (EM-PCA) instead of shrinking the gene set to what every line measured
  scores    every source row -> k numbers by least squares on the genes that source measured
  pairs     the destination d rotates over every included line; inputs are the other lines' scores for the same
            knockdown (0 + missing flag where a line did not keep it; d's own slot is always missing):
              (a) programs   target = d's program scores
              (b) corrected  target = scores of d's residual  R = measured d - Model 2 built without d
                             (tau2 refit without d); Model 2's own scores are an extra input
  ridge     one small ridge per program and variant; penalty from a fixed grid by 5-fold CV grouped by knockdown
            (all pairs of a knockdown share a fold): this measures transfer to unseen knockdowns in known lines,
            not to unseen lines. The H1 run in the notebook is validation, not an untouched exam.

predict
  (a) programs_k<k>  = loadings x predicted scores on the program genes; NaN elsewhere (no prediction)
  (b) corrected_k<k> = Model 2 + loadings x predicted residual scores on the program genes; Model 2 elsewhere
  Rows with all slots present (never seen in training, where the destination's slot is always missing) are
  predicted with each slot masked in turn, then averaged. For corrected, each mask also rebuilds Model 2 without
  that slot (tau2 refit too), exactly as training built its baseline, and the output averages
  (Model 2 without the slot + correction). Other rows are unchanged.
  <variant>_k<k>_spread.npy stores the program-space SD across masks (NaN on unaffected rows, not an SE).

Choice of k: besides the fixed grid, k80 = the smallest k whose programs hold 80% of the variance ABOVE THE NOISE.
  The noise spectrum comes from a table of pure measurement noise (each entry drawn with its own SE), so the
  comparison respects that some genes are far noisier than others; components above it carry signal.
  The knocked-down gene's own column is NaN (ontarget_lfc from Model 2 applies), and targets that no line kept stay
  NaN ("no information"). No SE is written for the ridge outputs: the generator's SE gate is off for them (z* = 0).
  Nothing about the destination enters, so the Delta tables for A, B and C are identical.

Example: python vcc3_ridge.py fit --filters-dir ~/vcc/lib/filters --sources-dir ~/vcc/lib/sources \
             --genes gene_names.csv --holdout H1 --out models/validation
         python vcc3_ridge.py predict --model models/validation --baseline validation/model2/H1 --lines H1 --out validation/pred
"""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

try:
    from anndata.io import read_elem
except ImportError:  # older anndata
    from anndata.experimental import read_elem

import vcc2_meta as meta
from vcc2_common import read_names, save_json, seed_for

SCHEMA = 'vcc_2026_3_ridge_v1'
BACKGROUND = dict(HCT116='HCT116', HEK293T='HEK293T', K562_gwps='K562', CD4_Stim48hr='CD4_T', H1='H1',
                  RPE1='RPE1', HepG2='HepG2', Jurkat='Jurkat')      # Replogle RPE1; Nadig HepG2 / Jurkat
LINES = list(BACKGROUND)                                   # one source per background: T cells = Stim48hr only
GENE_SET_SOURCES = ['HCT116', 'HEK293T', 'K562_gwps']      # the program genes: measured in all three
PROGRAM_SOURCES = ['HCT116', 'HEK293T', 'K562_gwps', 'RPE1', 'HepG2', 'Jurkat']   # every non-H1 cell line learns programs
K_GRID = (5, 10, 20)                                       # 20 is the ceiling: 21 components clear the noise edge
PENALTIES = (1e-4, 1e-3, 1e-2, 1e-1, 1.0)                  # fixed grid, chosen per program by grouped CV
SIGNAL_SHARE = 0.8                                         # k80: share of the above-noise variance the programs keep
FOLDS = 5
MIN_CPM = 1.0                                              # same gene filter as the panel tables
SEED = 2026


# ---------------------------------------------------------------- data
def check_backgrounds(sources_dir):
    for name in LINES:
        b, _ = meta.background_of(sources_dir, name)
        if b != BACKGROUND[name]:
            raise ValueError(f'{name}: panel says {b!r}, vcc3 says {BACKGROUND[name]!r}')


def var_of(filters_dir, name):
    with h5py.File(Path(filters_dir) / f'{name}.h5ad', 'r') as f:
        return read_elem(f['var']), read_elem(f['obs'])


def program_genes(filters_dir, genes, names=GENE_SET_SOURCES):
    """Challenge genes measured (>= MIN_CPM in controls) in every program source, in gene list order."""
    keep = np.ones(len(genes), bool)
    for name in names:
        var, _ = var_of(filters_dir, name)
        cpm = dict(zip(var.index.astype(str), var['ctrl_cpm'].to_numpy(float)))
        keep &= np.array([cpm.get(g, np.nan) >= MIN_CPM for g in genes])
    return np.asarray(genes)[keep]


class Line:
    """One source's QC-passed knockdowns on the program genes: lfc and se (NaN = not measured there)."""
    def __init__(self, filters_dir, name, genes, block=512, with_se=True):
        self.name, self.background = name, BACKGROUND[name]
        var, obs = var_of(filters_dir, name)
        self.targets = obs.index.astype(str).to_numpy()
        pos = {g: i for i, g in enumerate(var.index.astype(str))}
        cols = np.array([pos.get(g, -1) for g in genes]); cpm = np.full(len(genes), np.nan)
        ok = cols >= 0; cpm[ok] = var['ctrl_cpm'].to_numpy(float)[cols[ok]]
        self.measured = ok & (cpm >= MIN_CPM)
        src, dst = cols[self.measured], np.flatnonzero(self.measured)
        shape = (len(self.targets), len(genes))
        self.lfc = np.full(shape, np.nan, np.float32); self.se = np.full(shape, np.nan, np.float32) if with_se else None
        with h5py.File(Path(filters_dir) / f'{name}.h5ad', 'r') as f:   # row blocks: never the whole matrix at once
            for start in range(0, shape[0], block):
                stop = min(start + block, shape[0])
                self.lfc[start:stop, dst] = f['X'][start:stop][:, src]
                if with_se: self.se[start:stop, dst] = f['layers/se'][start:stop][:, src]
        gp = {g: i for i, g in enumerate(genes)}
        self.own = np.array([gp.get(t, -1) for t in self.targets])
        self.row = {t: i for i, t in enumerate(self.targets)}


def own_set(values, own, fill):
    """Copy with each row's own knocked-down gene set to `fill` (0 for programs/scores, NaN for averaging)."""
    out = np.array(values, np.float32, copy=True); r = np.flatnonzero(own >= 0)
    out[r, own[r]] = fill
    return out


# ---------------------------------------------------------------- programs and scores
def fit_programs(lines, k, iters=25, tol=1e-4):
    """Uncentred PCA: the top-k right singular vectors of the stacked knockdown x gene table (own gene = 0).
    Genes a line did not measure are missing, not zero: they are filled from the current rank-k fit and the SVD is
    redone (EM-PCA) until the fill stops changing. With nothing missing this is one plain SVD."""
    from sklearn.utils.extmath import randomized_svd
    Y = np.vstack([own_set(l.lfc, l.own, 0.0) for l in lines])
    observed = np.isfinite(Y); filled = np.where(observed, Y, 0).astype(np.float32)
    for _ in range(iters if not observed.all() else 1):
        u, s, vt = randomized_svd(filled, n_components=k, n_iter=7, random_state=SEED)
        if observed.all(): break
        fit = ((u * s) @ vt).astype(np.float32)
        change = float(np.sqrt(np.mean((fit[~observed] - filled[~observed]) ** 2)))
        filled[~observed] = fit[~observed]
        if change < tol: break
    total = float(np.square(filled, dtype=np.float64).sum())
    return vt.T.astype(np.float32), s, total, filled


def scores(values, W, own=None):
    """Least-squares program scores of each row on the genes it has (own gene counted as 0). NaN rows -> NaN."""
    Y = own_set(values, own, 0.0) if own is not None else np.asarray(values, np.float32)
    obs = np.isfinite(Y); out = np.full((len(Y), W.shape[1]), np.nan, np.float32)
    if not len(Y): return out
    _, group = np.unique(np.packbits(obs, axis=1), axis=0, return_inverse=True)
    for g in np.unique(group):
        r = np.flatnonzero(group.ravel() == g); m = obs[r[0]]
        if m.sum() < W.shape[1]: continue                   # too few genes to place the row
        Wm = W[m].astype(np.float64)
        out[r] = np.linalg.solve(Wm.T @ Wm, Wm.T @ Y[np.ix_(r, m)].T.astype(np.float64)).T
    return out


# ---------------------------------------------------------------- Model 2 baseline without the destination
def model2_tau2(sources_dir, all_genes, excluded):
    """Model 2's per-gene tau2, refit with the excluded backgrounds removed (same rule as vcc2_meta.fit)."""
    sources = [(n, b) for n, b in meta.included(sources_dir, LINES) if b not in excluded]
    panel = meta.Panel(sources_dir, sources, all_genes)
    tau2, _ = meta.tau2_dl(panel.lfc, panel.se, panel.own_col)
    est = np.isfinite(tau2)
    return np.where(est, tau2, np.median(tau2[est]))


def model2_rows(targets, own, sources, tau2, chunk=1024):
    """Model 2 prediction for `targets` (own gene column per target, -1 = none) from `sources` only."""
    T, G = len(targets), len(tau2)
    out = np.full((T, G), np.nan, np.float32)
    rows = [np.array([o.row.get(t, -1) for t in targets]) for o in sources]
    for start in range(0, T, chunk):
        idx = np.arange(start, min(start + chunk, T))
        lfc = np.full((len(sources), len(idx), G), np.nan, np.float32); se = lfc.copy()
        for s, (o, r) in enumerate(zip(sources, rows)):
            have = r[idx] >= 0
            lfc[s, have] = o.lfc[r[idx][have]]; se[s, have] = o.se[r[idx][have]]
        o_ = own[idx]; j = np.flatnonzero(o_ >= 0)
        lfc[:, j, o_[j]] = np.nan
        out[idx] = meta.combine(lfc, se, tau2)[0]
    return out


def model2_baseline(dest, others, tau2, chunk=1024):
    """Model 2 prediction of every knockdown of `dest` from the other lines only (own gene NaN)."""
    return model2_rows(dest.targets, dest.own, others, tau2, chunk)


def signal_k(lines, share=SIGNAL_SHARE, seed=SEED, filled=None):
    """k80: the smallest k whose components hold `share` of the variance above the noise (parallel analysis).

    The noise spectrum comes from a table of pure measurement noise of the same shape, each entry drawn N(0, SE^2)
    with its own SE (own gene 0 in both). A component carries signal when its eigenvalue exceeds the LARGEST noise
    eigenvalue (the noise edge; comparing rank by rank is biased, since signal components push the data's noise bulk
    down the ranks). k80 accumulates the signal components' excess over the noise edge until `share` is reached."""
    Y = np.vstack([own_set(l.lfc, l.own, 0.0) for l in lines]).astype(np.float64)
    if filled is not None: Y = filled.astype(np.float64)                      # EM-filled entries carry no noise
    sd = np.vstack([own_set(l.se, l.own, 0.0) for l in lines]).astype(np.float64)
    sd[~np.isfinite(sd)] = 0
    data = np.clip(np.linalg.eigvalsh(Y.T @ Y)[::-1], 0, None)
    noise_table = np.random.default_rng(seed_for(seed, 'signal_k')).standard_normal(Y.shape) * sd; del Y, sd
    noise = np.clip(np.linalg.eigvalsh(noise_table.T @ noise_table)[::-1], 0, None); del noise_table
    edge = float(noise[0])
    n_signal = int((data > edge).sum())
    excess = data[:n_signal] - edge
    cumulative = np.cumsum(excess) / excess.sum() if n_signal else np.array([1.0])
    k80 = int(np.searchsorted(cumulative, share) + 1) if n_signal else 1
    return dict(k80=k80, share=share, n_signal_components=n_signal, noise_edge=edge, noise_share_of_total=float(noise.sum() / data.sum()),
                signal_share_kept={int(k): float(cumulative[min(k, len(cumulative)) - 1]) for k in (*K_GRID, k80)})


# ---------------------------------------------------------------- ridge
def ridge_fit(X, y, penalty):
    """Ridge on standardised inputs; the intercept is not penalised. Loss: mean squared error + penalty x |b|^2."""
    mu = X.mean(0); sd = X.std(0); sd[sd == 0] = 1
    Z = (X - mu) / sd; n = len(y)
    coef = np.linalg.solve(Z.T @ Z / n + penalty * np.eye(X.shape[1]), Z.T @ (y - y.mean()) / n)
    return dict(mu=mu, sd=sd, coef=coef, intercept=float(y.mean()), penalty=float(penalty))


def ridge_predict(fit, X):
    return ((X - fit['mu']) / fit['sd']) @ fit['coef'] + fit['intercept']


def folds_for(targets, folds=FOLDS):
    """Fold of each pair by its knockdown: every pair of one knockdown lands in the same fold."""
    return np.array([seed_for(SEED, 'fold', t) % folds for t in targets])


def fit_program(X, y, fold, penalties=PENALTIES):
    """Pick the penalty by grouped CV, refit on all pairs. Returns (fit, cv table, out-of-fold predictions)."""
    cv, oof = {}, {}
    for lam in penalties:
        pred = np.empty(len(y))
        for f in np.unique(fold):
            test = fold == f
            pred[test] = ridge_predict(ridge_fit(X[~test], y[~test], lam), X[test])
        cv[lam] = float(np.mean((pred - y) ** 2)); oof[lam] = pred
    best = min(cv, key=cv.get)
    return ridge_fit(X, y, best), cv, oof[best]


def features(slot_scores, slot_missing, j, extra=None):
    """Inputs of program j: each slot's score (0 when missing) and missing flag, plus an optional extra column."""
    cols = [np.where(slot_missing, 0.0, np.nan_to_num(slot_scores[:, :, j])), slot_missing.astype(float)]
    if extra is not None: cols.append(extra[:, [j]])
    return np.hstack(cols)


# ---------------------------------------------------------------- fit
def slot_inputs(targets, lines, W, cache, exclude=None):
    """Scores of every slot line for `targets`: (n x slots x k) and missing flags (n x slots)."""
    n, k = len(targets), W.shape[1]
    S = np.zeros((n, len(lines), k), np.float32); missing = np.ones((n, len(lines)), bool)
    for s, line in enumerate(lines):
        if line.name == exclude: continue
        r = np.array([line.row.get(t, -1) for t in targets]); have = r >= 0
        sc = cache[line.name][r[have]]
        ok = np.isfinite(sc).all(1); idx = np.flatnonzero(have)[ok]
        S[idx, s] = sc[ok]; missing[idx, s] = False
    return S, missing


def fit(filters_dir, sources_dir, genes_csv, output, holdout=None, k_grid=K_GRID, signal_share=SIGNAL_SHARE):
    check_backgrounds(sources_dir)
    all_genes = read_names(genes_csv, 'gene_name')
    names = [n for n in LINES if BACKGROUND[n] != holdout]
    pgenes = program_genes(filters_dir, all_genes)
    print(f'fit: lines {names} (held out: {holdout}); {len(pgenes):,} program genes', flush=True)
    lines = [Line(filters_dir, n, pgenes) for n in names]
    program_lines = [l for l in lines if l.name in PROGRAM_SOURCES]
    print(f"  program sources: {[l.name for l in program_lines]}", flush=True)
    W_all, sing, total, filled = fit_programs(program_lines, max(k_grid))
    signal = None
    if signal_share:
        signal = signal_k(program_lines, signal_share, filled=filled)
        print(f"  k80: {signal['k80']} programs hold {signal_share:.0%} of the above-noise variance "
              f"({signal['n_signal_components']} components above noise; noise = {signal['noise_share_of_total']:.0%} of the total)", flush=True)
        if signal['k80'] < max(k_grid) and signal['k80'] not in k_grid:   # never above the grid's ceiling
            k_grid = tuple(sorted(set(k_grid) | {signal['k80']}))
    del filled
    explained = (sing ** 2 / total).tolist()
    print('  variance explained by the top programs: ' + ', '.join(f'k={k}: {sum(explained[:k]):.1%}' for k in k_grid), flush=True)

    # Destination rows and the Model 2 baseline built without each destination.
    gpos = {g: i for i, g in enumerate(all_genes)}; pidx = np.array([gpos[g] for g in pgenes])
    pairs = []
    for d in lines:
        others = [o for o in lines if o is not d]
        kept_elsewhere = np.array([any(t in o.row for o in others) for t in d.targets])
        if not kept_elsewhere.any(): continue
        tau2 = model2_tau2(sources_dir, all_genes, {holdout, d.background} - {None})[pidx]
        base = model2_baseline(d, others, tau2)
        rows = np.flatnonzero(kept_elsewhere)
        pairs.append(dict(dest=d, rows=rows, base=base[rows]))
        print(f'  destination {d.name}: {len(rows):,} knockdowns kept in >= 1 other line', flush=True)
    for l in lines: l.se = None                                  # only the baselines needed it

    out = Path(output); out.mkdir(parents=True, exist_ok=True)
    np.save(out / 'program_genes.npy', pgenes); np.save(out / 'loadings.npy', W_all)
    report, params = [], {}
    for k in k_grid:
        W = W_all[:, :k]
        cache = {l.name: scores(l.lfc, W, l.own) for l in lines}
        X_rows = {'programs': [], 'corrected': []}; Y_rows = {'programs': [], 'corrected': []}; meta_rows = []
        for p in pairs:
            d, rows = p['dest'], p['rows']; targets = d.targets[rows]
            S, missing = slot_inputs(targets, lines, W, cache, exclude=d.name)
            target_scores = cache[d.name][rows]
            resid = own_set(d.lfc[rows] - p['base'], d.own[rows], 0.0)
            base_scores = scores(p['base'], W, d.own[rows])
            resid_scores = scores(resid, W)
            ok = np.isfinite(target_scores).all(1) & np.isfinite(resid_scores).all(1) & np.isfinite(base_scores).all(1)
            meta_rows.append(dict(dest=d, rows=rows[ok], targets=targets[ok], base=p['base'][ok]))
            X_rows['programs'].append((S[ok], missing[ok], None)); Y_rows['programs'].append(target_scores[ok])
            X_rows['corrected'].append((S[ok], missing[ok], base_scores[ok])); Y_rows['corrected'].append(resid_scores[ok])
        targets_all = np.concatenate([m['targets'] for m in meta_rows]); fold = folds_for(targets_all)
        params[k] = {}
        for variant in ['programs', 'corrected']:
            S = np.concatenate([x[0] for x in X_rows[variant]]); M = np.concatenate([x[1] for x in X_rows[variant]])
            extra = None if variant == 'programs' else np.concatenate([x[2] for x in X_rows[variant]])
            Y = np.concatenate(Y_rows[variant]); oof = np.empty_like(Y); fits = []
            for j in range(k):
                f, cv, oof[:, j] = fit_program(features(S, M, j, extra), Y[:, j], fold)
                fits.append(f)
            params[k][variant] = fits
            report += cv_report(k, variant, oof, Y, W, meta_rows)
        print(f'  k={k}: {len(targets_all):,} training pairs, {len(np.unique(targets_all)):,} knockdowns', flush=True)
    save_params(out / 'ridge.npz', params)
    frame = pd.DataFrame(report); frame.to_csv(out / 'cv_report.csv', index=False)
    model = dict(schema=SCHEMA, filters_dir=str(Path(filters_dir).resolve()), sources_dir=str(Path(sources_dir).resolve()),
                 genes_csv=str(Path(genes_csv).resolve()), heldout_background=holdout, lines=names,
                 backgrounds=[BACKGROUND[n] for n in names], program_sources=[l.name for l in program_lines],
                 gene_set_sources=GENE_SET_SOURCES, k_grid=list(k_grid),
                 penalties=list(PENALTIES), folds=FOLDS, n_program_genes=len(pgenes),
                 variance_explained={k: float(sum(explained[:k])) for k in k_grid},
                 own_gene='set to 0 in programs and scores; handled by ontarget_lfc',
                 cv='grouped by knockdown: a selection estimate for unseen knockdowns in known lines, not a transfer estimate',
                 signal_k=signal)
    save_json(out / 'model.json', model)
    summary = frame[frame.destination == 'all'].pivot(index='k', columns='variant', values='gene_mse')
    print('  CV gene-space MSE on held-out knockdowns (all destinations):'); print(summary.round(5).to_string())
    return model


def cv_report(k, variant, oof, Y, W, meta_rows):
    """Gene-space error of out-of-fold predictions vs Model 2 and vs predicting no change."""
    out, start, totals = [], 0, dict(n=0, model2=0.0, zero=0.0, pred=0.0)
    for m in meta_rows:
        n = len(m['rows']); d = m['dest']; s = oof[start:start + n]; start += n
        measured = own_set(d.lfc[m['rows']], d.own[m['rows']], np.nan)
        rebuilt = s @ W.T
        pred = rebuilt if variant == 'programs' else m['base'] + rebuilt
        ok = np.isfinite(measured) & np.isfinite(pred) & np.isfinite(m['base'])
        e = dict(n=int(ok.sum()), model2=float(((measured - m['base'])[ok] ** 2).sum()),
                 zero=float((measured[ok] ** 2).sum()), pred=float(((measured - pred)[ok] ** 2).sum()))
        for key in totals: totals[key] += e[key]
        out.append(dict(k=k, variant=variant, destination=d.name, pairs=n, gene_mse=e['pred'] / e['n'],
                        model2_mse=e['model2'] / e['n'], no_change_mse=e['zero'] / e['n']))
    score_r2 = 1 - ((oof - Y) ** 2).sum(0) / np.maximum(((Y - Y.mean(0)) ** 2).sum(0), 1e-12)
    out.append(dict(k=k, variant=variant, destination='all', pairs=int(len(Y)), gene_mse=totals['pred'] / totals['n'],
                    model2_mse=totals['model2'] / totals['n'], no_change_mse=totals['zero'] / totals['n'],
                    score_r2_mean=float(score_r2.mean())))
    return out


def save_params(path, params):
    flat = {}
    for k, variants in params.items():
        for variant, fits in variants.items():
            for j, f in enumerate(fits):
                for key, value in f.items(): flat[f'k{k}/{variant}/{j}/{key}'] = np.asarray(value)
    np.savez(path, **flat)


def load_params(path, k, variant):
    z = np.load(path); fits = []
    for j in range(k):
        fits.append({key: z[f'k{k}/{variant}/{j}/{key}'] for key in ['mu', 'sd', 'coef', 'intercept', 'penalty']})
    return fits


# ---------------------------------------------------------------- predict
def predict_scores(fits, S, missing, extra):
    return np.stack([ridge_predict(f, features(S, missing, j, extra)) for j, f in enumerate(fits)], 1)


def leave_one_slot_out(fits, S, missing, extra, extra_by_slot=None):
    """Average masked predictions only for full rows; spread is in program space.

    Spread has one row per target, with NaN for rows not subject to averaging. For corrected, `extra_by_slot[s]`
    holds the full rows' Model 2 scores rebuilt WITHOUT slot s (as in training); without it, `extra` is reused.
    """
    before = predict_scores(fits, S, missing, extra)
    after = before.copy()
    full = ~missing.any(1)
    spread = np.full(before.shape, np.nan, dtype=before.dtype)
    if full.any():
        loo = []
        for slot in range(S.shape[1]):
            m = missing[full].copy(); m[:, slot] = True
            X = S[full].copy(); X[:, slot] = 0
            e = None if extra is None else (extra_by_slot[slot] if extra_by_slot is not None else extra[full])
            loo.append(predict_scores(fits, X, m, e))
        loo = np.stack(loo)
        after[full] = loo.mean(0)
        spread[full] = loo.std(0)
    return before, after, spread, full


def mean_abs_delta(delta):
    finite = np.isfinite(delta)
    return float(np.abs(delta[finite]).mean(dtype=np.float64)) if finite.any() else None


def predict(model_dir, baseline_dir, lines_out, output):
    m = Path(model_dir); model = json.loads((m / 'model.json').read_text())
    check_backgrounds(model['sources_dir'])
    if model.get('schema') != SCHEMA: raise ValueError('incompatible model; fit with this notebook')
    b = Path(baseline_dir); report2 = json.loads((b / 'prediction_report.json').read_text())
    if report2.get('heldout_background') != model['heldout_background']:
        raise ValueError(f'baseline held out {report2.get("heldout_background")}, model held out {model["heldout_background"]}')
    pgenes = np.load(m / 'program_genes.npy').astype(str); W_all = np.load(m / 'loadings.npy')
    genes = np.load(b / 'genes.npy').astype(str); targets = np.load(b / 'targets.npy').astype(str)
    re = np.load(b / 're.npy'); gpos = {g: i for i, g in enumerate(genes)}; pidx = np.array([gpos[g] for g in pgenes])
    lines = [Line(model['filters_dir'], n, pgenes) for n in model['lines']]
    own_p = np.array([{g: i for i, g in enumerate(pgenes)}.get(t, -1) for t in targets])
    own_all = np.array([gpos.get(t, -1) for t in targets])
    out = Path(output)
    for line in lines_out:
        d = out / line; earlier = d / 'prediction_report.json'
        if d.exists() and any(d.iterdir()) and not (earlier.exists() and json.loads(earlier.read_text()).get('schema') == SCHEMA):
            raise ValueError(f'choose an empty prediction directory: {d}')
    arrays, coverage, spreads, averaging = {}, {}, {}, {}
    # Rows with every slot present: Model 2 rebuilt without each slot in turn (tau2 refit), for the corrected masks.
    present = np.array([[t in l.row for l in lines] for t in targets])
    full_rows = np.flatnonzero(present.all(1))
    all_genes = read_names(model['genes_csv'], 'gene_name')
    without_slot = []
    if len(full_rows):
        for s_, line in enumerate(lines):
            tau2 = model2_tau2(model['sources_dir'], all_genes, {model['heldout_background'], line.background} - {None})
            tau2 = tau2[np.array([{g: i for i, g in enumerate(all_genes)}[g] for g in pgenes])]
            without_slot.append(model2_rows(targets[full_rows], own_p[full_rows], [l for l in lines if l is not line], tau2))
    for k in model['k_grid']:
        W = W_all[:, :k]
        cache = {l.name: scores(l.lfc, W, l.own) for l in lines}
        S, missing = slot_inputs(targets, lines, W, cache)
        none = missing.all(1)                                    # no line kept this knockdown: no information
        base_scores = scores(re[:, pidx], W, own_p)
        for variant in ['programs', 'corrected']:
            fits = load_params(m / 'ridge.npz', k, variant)
            extra = None if variant == 'programs' else np.nan_to_num(base_scores)
            extra_by_slot = None
            full = ~missing.any(1)
            sel = np.searchsorted(full_rows, np.flatnonzero(full))      # full rows are a subset of the rows with every line
            if variant == 'corrected' and full.any():
                assert np.isin(np.flatnonzero(full), full_rows).all(), 'every full row needs its per-slot baselines'
                extra_by_slot = [np.nan_to_num(scores(b[sel], W, own_p[full])) for b in without_slot]
            before, s, spread, full = leave_one_slot_out(fits, S, missing, extra, extra_by_slot)
            rebuilt = (s @ W.T).astype(np.float32)
            if variant == 'programs':
                delta = np.full(re.shape, np.nan, np.float32); delta[:, pidx] = rebuilt
            else:
                delta = re.copy(); delta[:, pidx] = re[:, pidx] + rebuilt   # NaN baseline stays NaN
                if full.any():
                    # each mask predicts the residual of Model 2 WITHOUT that slot: average (that baseline + correction)
                    with warnings.catch_warnings():              # genes no remaining line measured stay NaN
                        warnings.simplefilter('ignore', RuntimeWarning)
                        base_avg = np.nanmean(np.stack([b[sel] for b in without_slot]), 0)
                    delta[np.ix_(np.flatnonzero(full), pidx)] = base_avg + rebuilt[full]
            delta[none] = np.nan
            r = np.flatnonzero(own_all >= 0); delta[r, own_all[r]] = np.nan
            name = f'{variant}_k{k}'
            arrays[name] = delta
            spreads[f'{name}_spread'] = spread
            prior = delta[full].copy()
            raw_rebuilt = (before[full] @ W.T).astype(np.float32)
            prior[:, pidx] = raw_rebuilt if variant == 'programs' else re[full][:, pidx] + raw_rebuilt
            prior = own_set(prior, own_all[full], np.nan)
            before_abs, after_abs = mean_abs_delta(prior), mean_abs_delta(delta[full])
            # The correction alone (programs x predicted scores), own gene excluded: for corrected, |delta| is mostly Model 2.
            corr_before = mean_abs_delta(own_set(raw_rebuilt, own_p[full], np.nan)) if full.any() else None
            corr_after = mean_abs_delta(own_set(rebuilt[full], own_p[full], np.nan)) if full.any() else None
            averaging[name] = dict(affected_rows=int(full.sum()), affected_targets=targets[full].tolist(),
                                   mean_abs_delta_before=before_abs, mean_abs_delta_after=after_abs,
                                   mean_abs_correction_before=corr_before, mean_abs_correction_after=corr_after,
                                   after_smaller=None if corr_before is None else corr_after < corr_before,
                                   corrected_baseline='Model 2 rebuilt without each masked slot' if variant == 'corrected' else None)
            print(f'  {name}: leave-one-slot-out on {full.sum()} rows; mean |correction| {corr_before} -> {corr_after}', flush=True)
        coverage[k] = int((~none).sum())
    models = ['re'] + list(arrays) + ['zero']
    report = dict(schema=SCHEMA, models=models, heldout_background=model['heldout_background'],
                  source_backgrounds=model['backgrounds'], model_dir=str(m.resolve()), baseline_dir=str(b.resolve()),
                  re='Model 2 (vcc2_meta.py), copied unchanged; re_se is its SE',
                  zero='all-zero effects with the same shape as re; the generator still applies the knockdown',
                  programs='loadings x ridge-predicted program scores on the program genes; NaN elsewhere',
                  corrected='Model 2 + loadings x ridge-predicted residual scores on the program genes; Model 2 elsewhere',
                  se='none for the ridge outputs: use the generator with --z-grid 0',
                  leave_one_slot_out=averaging,
                  averaging_metric='mean absolute effect over finite output genes of affected rows; own gene excluded',
                  spread='population SD across masked predictions in program space; targets x k, NaN on unaffected rows; not an SE',
                  n_program_genes=len(pgenes), targets_with_any_line=coverage[model['k_grid'][0]],
                  no_information=[t for t, n in zip(targets.tolist(), missing.all(1)) if n])
    for line in lines_out:
        d = out / line; d.mkdir(parents=True, exist_ok=True)
        for old in [*d.glob('*.npy'), d / 'prediction_report.json']: old.unlink(missing_ok=True)
        for name in ['genes', 'targets', 're', 're_se', 'n_backgrounds', 'ontarget_lfc']:
            np.save(d / f'{name}.npy', np.load(b / f'{name}.npy'))
        for name, value in arrays.items(): np.save(d / f'{name}.npy', value)
        for name, value in spreads.items(): np.save(d / f'{name}.npy', value)
        np.save(d / 'zero.npy', np.zeros_like(re))
        save_json(d / 'prediction_report.json', dict(report, line=line))
    print(f'predict: {len(targets)} targets -> {", ".join(str(out / l) for l in lines_out)}; models {models}')
    print(f'  {report["targets_with_any_line"]} targets kept in >= 1 line; no information: {report["no_information"] or "none"}')
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    f = sub.add_parser('fit')
    for key in ['filters-dir', 'sources-dir', 'genes', 'out']: f.add_argument('--' + key, required=True)
    f.add_argument('--holdout'); f.add_argument('--k-grid', nargs='+', type=int, default=list(K_GRID))
    f.add_argument('--signal-share', type=float, default=SIGNAL_SHARE, help='k80: share of the above-noise variance; 0 = grid only')
    pr = sub.add_parser('predict')
    for key in ['model', 'baseline', 'out']: pr.add_argument('--' + key, required=True)
    pr.add_argument('--lines', nargs='+', required=True)
    a = p.parse_args()
    if a.cmd == 'fit': fit(a.filters_dir, a.sources_dir, a.genes, a.out, a.holdout, tuple(a.k_grid), a.signal_share)
    else: predict(a.model, a.baseline, a.lines, a.out)

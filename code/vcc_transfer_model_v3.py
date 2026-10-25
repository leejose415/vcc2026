"""Context transfer with normal-prior shrinkage and nested biological LOBO selection.

Tune ridge penalty, expression gates and similarity temperature ONLY on held-out
backgrounds. Final refit uses all eligible targets. source_mean is the default
unless ridge wins >=3 of the four original background folds. Extra Replogle
backgrounds expand coverage and can supply additional LOBO diagnostics.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from scipy.special import expit
from vcc_common import (Effects, PROGRAMS, estimate_prior, load_profile, program_levels, read_names,
                        save_json, seed_for, shrink, valid_se)

SRC_NAMES = ['source_lfc', 'source_se', 'se_missing', 'source_shrunk_lfc',
             'source_baseline', 'baseline_difference', 'source_knockdown_gene_baseline',
             'knockdown_gene_baseline_difference', 'source_ontarget_lfc', 'ontarget_missing',
             'legacy_effect', 'context_pearson', 'agreement', 'legacy_source_baseline',
             'gated_source_effect', 'output_expression_gate', 'knockdown_expression_gate',
             'knockdown_detection_gate', 'source_scale'] + ['program_difference_'+p for p in PROGRAMS]
TGT_NAMES = ['target_baseline', 'target_knockdown_gene_baseline', 'knockdown_gene_missing',
             'legacy_target_baseline'] + ['target_program_'+p for p in PROGRAMS]
DEFAULT_KNOBS = dict(temperature=0.5, gate_center=1.8, gate_scale=1.0,
                     detection_center=0.05, detection_scale=0.05)


def transfer_profile(profile, mode='bulk', declared_convention=None, label='profile'):
    p = dict(profile)
    convention = declared_convention or str(p.get('baseline_convention', ''))
    if mode == 'bulk':
        if 'pooled_cpm' in p: base = p['pooled_cpm']
        elif 'bulk_cpm' in p: base = p['bulk_cpm']
        elif convention == 'ratio_of_sums_cpm': base = p['mean_cpm']
        else: raise ValueError(f'{label}: no verified pooled CPM baseline')
    elif mode == 'per_cell':
        if bool(p.get('legacy_baseline', False)) and convention != 'arithmetic_mean_per_cell_cpm':
            raise ValueError(f'{label}: no verified per-cell baseline')
        base = p['mean_cpm']
    else: raise ValueError('baseline mode must be bulk or per_cell')
    base = np.asarray(base, float)
    if not np.isfinite(base).any() or np.any(base[np.isfinite(base)] < 0) or np.isinf(base).any():
        raise ValueError(f'{label}: invalid baseline')
    p['baseline_cpm'] = base; p['baseline_log'] = np.log1p(base)
    _, levels = program_levels(base, p['genes'])
    p['baseline_program_levels'] = np.nan_to_num(levels, nan=0)
    p['baseline_mode'] = mode
    return p


def centred_similarity(a, b):
    valid = np.isfinite(a) & np.isfinite(b)
    if valid.sum() < 2: return 0.0
    x, y = np.asarray(a[valid], float).copy(), np.asarray(b[valid], float).copy()
    x -= x.mean(); y -= y.mean(); den = np.linalg.norm(x)*np.linalg.norm(y)
    return float(np.clip(np.dot(x, y)/den, -1, 1)) if den > 1e-12 else 0.0


class Context:
    def __init__(self, spec, genes, allow_legacy=False, baseline_mode='bulk'):
        self.name = spec['name']; self.background = spec['background']; self.platform = spec.get('platform', 'unknown')
        self.priority = spec.get('source_priority', 0); self.spec = spec
        self.effects = Effects(spec['effects'], genes, allow_legacy)
        p = load_profile(spec.get('profile', str(Path(spec['effects'])/'control_profile.npz')), genes)
        self.profile = transfer_profile(p, baseline_mode, spec.get('baseline_convention'), self.name)


def load_contexts(config, genes, excluded=(), allow_legacy=False, baseline_mode='bulk'):
    # Exclude before opening effects/profiles; all states of a held-out line stay out.
    return [Context(c, genes, allow_legacy, baseline_mode) for c in config['contexts'] if c['background'] not in excluded]


def choose_sources(contexts, target_background=None):
    groups = {}
    for c in contexts:
        if c.background == target_background: continue
        if c.background not in groups or (c.priority, c.name) > (groups[c.background].priority, groups[c.background].name):
            groups[c.background] = c
    return sorted(groups.values(), key=lambda c: c.name)


def shrunk(c, kd, gi, cap):
    l, se = c.effects.get(kd, gi)
    return shrink(l, se, cap=cap, prior_variance=c.effects.prior_for(kd))


def generic_signatures(contexts, targets, genes, cap):
    out = {}; gi = np.arange(len(genes))
    for c in contexts:
        total, count = np.zeros(len(genes)), np.zeros(len(genes))
        for kd in targets:
            if kd not in c.effects.pos: continue
            sh = shrunk(c, kd, gi, cap); ok = np.isfinite(sh)
            total[ok] += sh[ok]; count[ok] += 1
        out[c.name] = np.divide(total, count, out=np.zeros_like(total), where=count > 0).astype(np.float32)
    return out


def features(target_profile, target_platform, sources, kd, gi, genes, generics, cap=8.0,
             line_features=True, knobs=None, scales=None):
    """Only source effects + untreated target context; no target perturbation labels."""
    gi = np.asarray(gi, int); gp = {g: i for i, g in enumerate(genes)}; ki = gp.get(kd)
    knobs = dict(DEFAULT_KNOBS, **(knobs or {})); scales = scales or {}
    if knobs['temperature'] <= 0 or knobs['gate_scale'] <= 0 or knobs['detection_scale'] <= 0:
        raise ValueError('temperatures and gate widths must be positive')
    if 'baseline_mode' not in target_profile: target_profile = transfer_profile(target_profile)
    tb = target_profile['baseline_log']; missing = ki is None or not np.isfinite(tb[ki])
    tk = 0.0 if missing else tb[ki]
    # Unknown target expression/detection passes through; unknown is not measured off.
    kgate = 1.0 if missing else (0.0 if target_profile['baseline_cpm'][ki] == 0 else float(expit((tk-knobs['gate_center'])/knobs['gate_scale'])))
    detection = target_profile.get('detection', np.full(len(genes), np.nan))
    det = np.nan if ki is None else detection[ki]
    dgate = 1.0 if not np.isfinite(det) else (0.0 if det == 0 else float(expit((det-knobs['detection_center'])/knobs['detection_scale'])))
    ogate = np.where(target_profile['baseline_cpm'][gi] == 0, 0.0, expit((tb[gi]-knobs['gate_center'])/knobs['gate_scale']))
    tp = target_profile['baseline_program_levels'] if line_features else np.zeros(len(PROGRAMS))
    tf = np.column_stack([tb[gi], np.full(len(gi), tk), np.full(len(gi), missing),
                          np.full(len(gi), bool(target_profile.get('legacy_baseline', False))),
                          np.broadcast_to(tp, (len(gi), len(tp)))])
    sv, sf, masks, gv, correlations = [], [], [], [], []
    for c in sources:
        l, se = c.effects.get(kd, gi); sb = c.profile['baseline_log']
        sk = sb[ki] if ki is not None and np.isfinite(sb[ki]) else 0.0
        own, _ = c.effects.get(kd, [ki]) if ki is not None and kd in c.effects.pos else (np.array([np.nan]), None)
        valid = np.isfinite(l) & np.isfinite(sb[gi]) & np.isfinite(tb[gi])
        sh = shrunk(c, kd, gi, cap); source_scale = scales.get(c.name, 1.0)
        # Down effects vanish as the responder turns off; source-supported induction remains available.
        gated = sh*source_scale*np.where(sh < 0, ogate, 1.0)*kgate*dgate
        corr = centred_similarity(tb, sb); correlations.append(corr)
        pdiff = tp-c.profile['baseline_program_levels'] if line_features else np.zeros(len(PROGRAMS))
        cols = [np.clip(l, -cap, cap), np.where(valid_se(se), np.clip(se, 0, cap), 0), ~valid_se(se), sh,
                sb[gi], tb[gi]-sb[gi], np.full(len(gi), sk), np.full(len(gi), tk-sk),
                np.full(len(gi), np.clip(own[0], -cap, cap)), np.full(len(gi), not np.isfinite(own[0])),
                np.full(len(gi), c.effects.legacy), np.full(len(gi), corr), np.zeros(len(gi)),
                np.full(len(gi), bool(c.profile.get('legacy_baseline', False))), gated, ogate,
                np.full(len(gi), kgate), np.full(len(gi), dgate), np.full(len(gi), source_scale),
                np.broadcast_to(pdiff, (len(gi), len(tp)))]
        sf.append(np.nan_to_num(np.column_stack(cols), nan=0, posinf=cap, neginf=-cap))
        sv.append(np.nan_to_num(gated, nan=0)); masks.append(valid); gv.append(generics[c.name][gi])
    if not sources: raise ValueError('no permitted sources')
    sf, sv, mask = np.stack(sf, 1), np.stack(sv, 1), np.stack(masks, 1)
    for j in range(len(sources)):
        other = mask.copy(); other[:, j] = False
        votes = other & (np.abs(sv) > 0.05); same = np.sign(sv) == np.sign(sv[:, [j]])
        sf[:, j, SRC_NAMES.index('agreement')] = np.divide((votes & same).sum(1), votes.sum(1),
                                                    out=np.full(len(gi), 0.5), where=votes.sum(1) > 0)
    logw = np.asarray(correlations)/knobs['temperature']; logw -= logw.max()
    weight = mask*np.exp(logw)[None, :]
    weight /= np.maximum(weight.sum(1, keepdims=True), 1e-30)
    generic = (np.stack(gv, 1)*weight).sum(1)
    return dict(sv=sv.astype(np.float32), sf=sf.astype(np.float32), mask=mask,
                weights=weight.astype(np.float32), tf=np.nan_to_num(tf).astype(np.float32), generic=generic.astype(np.float32))


def pad(f, n):
    d = n-f['sv'].shape[1]
    if d < 0: raise ValueError('invalid source padding')
    f = dict(f)
    for k in ['sv', 'mask', 'weights']: f[k] = np.pad(f[k], ((0, 0), (0, d)))
    f['sf'] = np.pad(f['sf'], ((0, 0), (0, d), (0, 0)))
    return f


def pooled(f):
    w = f['weights'][:, :, None]; mean = (f['sf']*w).sum(1)
    sd = np.sqrt(((f['sf']-mean[:, None, :])**2*w).sum(1)); base = (f['sv']*f['weights']).sum(1)
    ii = [SRC_NAMES.index('baseline_difference'), SRC_NAMES.index('knockdown_gene_baseline_difference')]
    inter = (f['sv'][:, :, None]*f['sf'][:, :, ii]*w).sum(1)
    return np.column_stack([mean, sd, f['tf'], inter, f['generic'], f['mask'].sum(1)]).astype(np.float32), base


def balanced_label_weights(label_se, task, contexts, noise_floor=0.25, precision_cap=10.0, enabled=True):
    if noise_floor <= 0 or precision_cap <= 0: raise ValueError('invalid label weighting')
    w = np.zeros(len(task)); active = np.unique(task)
    states = {c.background: sum(contexts[i].background == c.background for i in active) for c in contexts}
    for ti in active:
        use = task == ti; se = np.asarray(label_se[use], float); good = valid_se(se)
        precision = np.ones(len(se))
        if enabled and good.any():
            precision[good] = np.minimum(precision_cap, 1/(se[good]**2+noise_floor**2))
            precision[~good] = np.median(precision[good])
        w[use] = precision/precision.sum()/states[contexts[ti].background]
    return (w/w.mean()).astype(np.float32)


def training_targets(contexts, seed, max_kd):
    overlap = set()
    for c in contexts:
        others = set().union(*(set(s.effects.targets) for s in contexts if s.background != c.background))
        overlap |= set(c.effects.targets) & others
    kd = np.array(sorted(overlap)); np.random.default_rng(seed).shuffle(kd)
    if max_kd: kd = kd[:max_kd]
    if len(kd) < 5: raise ValueError('need >=5 overlapping successful targets in training backgrounds')
    return kd


def build_rows(contexts, targets, genes, generics, seed=0, per_kd=128, dropout=0.25, cap=8,
               line_features=True, knobs=None, scales=None, label_contexts=None,
               label_noise_floor=0.25, label_weight_cap=10.0, label_weighting=True):
    labels_contexts = contexts if label_contexts is None else label_contexts
    blocks, labels, errors, tasks = [], [], [], []
    smax = len({c.background for c in contexts})
    for ti, target in enumerate(labels_contexts):
        # Include low-baseline genes; the prediction gate handles down/up asymmetrically.
        eligible = np.flatnonzero(np.isfinite(target.profile['baseline_cpm']))
        for kd in targets:
            if kd not in target.effects.pos: continue
            src = choose_sources([c for c in contexts if kd in c.effects.pos], target.background)
            if not src: continue
            rng = np.random.default_rng(seed_for(seed, target.name, kd))
            keep = rng.random(len(src)) >= dropout
            if not keep.any(): keep[rng.integers(len(src))] = True
            src = [c for c, yes in zip(src, keep) if yes]
            gi = rng.choice(eligible[genes[eligible] != kd], min(per_kd, int((genes[eligible] != kd).sum())), replace=False)
            y, se = target.effects.get(kd, gi)
            f = features(target.profile, target.platform, src, kd, gi, genes, generics, cap, line_features, knobs, scales)
            valid = np.isfinite(y) & f['mask'].any(1)
            if valid.any():
                blocks.append(pad({k: v[valid] for k, v in f.items()}, smax)); labels.append(y[valid])
                errors.append(se[valid]); tasks.append(np.full(valid.sum(), ti))
    if not blocks: raise ValueError('no aligned training/validation examples')
    f = {k: np.concatenate([b[k] for b in blocks]) for k in blocks[0]}
    y, task = np.concatenate(labels), np.concatenate(tasks)
    w = balanced_label_weights(np.concatenate(errors), task, labels_contexts, label_noise_floor, label_weight_cap, label_weighting)
    return f, y, w


def fit_source_scales(contexts, targets, genes, cap, seed=0, budget=128):
    """One regularized global scale per source table, learned on training backgrounds only.

    Average moments across target backgrounds before fitting so cell/row count
    cannot make one background dominate. Shrink each scale toward 1 and bound it.
    """
    scales = {}; gi_all = np.arange(len(genes))
    for source in contexts:
        moments_by_bg = {}
        for target in contexts:
            if target.background == source.background: continue
            xx = xy = n = 0.0
            for kd in targets:
                if kd not in source.effects.pos or kd not in target.effects.pos: continue
                rng = np.random.default_rng(seed_for(seed, 'scale', target.name, kd))
                gi = rng.choice(gi_all[genes != kd], min(budget, int((genes != kd).sum())), replace=False)
                x = shrunk(source, kd, gi, cap); y, _ = target.effects.get(kd, gi)
                good = np.isfinite(x) & np.isfinite(y)
                xx += np.dot(x[good], x[good]); xy += np.dot(x[good], y[good]); n += good.sum()
            if n: moments_by_bg.setdefault(target.background, []).append((xx/n, xy/n))
        pairs = [np.mean(v, axis=0) for v in moments_by_bg.values()]
        if pairs:
            xx, xy = np.mean(pairs, axis=0)
            scales[source.name] = float(np.clip((xy+0.25)/(xx+0.25), 0.25, 2.0))
        else: scales[source.name] = 1.0
    return scales


def fit_ridge(f, y, w, alpha):
    x, base = pooled(f); scaler = StandardScaler().fit(x, sample_weight=w); active = scaler.var_ > 1e-8
    z = np.clip(scaler.transform(x)*active, -5, 5)
    # Gate the whole residual path (including its intercept), so an off target
    # cannot acquire effects through ungated context covariates.
    gate = f['sf'][:, 0, SRC_NAMES.index('knockdown_expression_gate')]*f['sf'][:, 0, SRC_NAMES.index('knockdown_detection_gate')]
    z = np.column_stack([z*gate[:, None], gate])
    model = Ridge(alpha=float(alpha), fit_intercept=False).fit(z, y-base, sample_weight=w)
    return dict(kind='ridge', estimator=model, scaler=scaler, active=active, alpha=float(alpha))


def predict(model, f):
    x, base = pooled(f)
    if model['kind'] == 'source_mean': return base
    z = np.clip(model['scaler'].transform(x)*model['active'], -5, 5)
    gate = f['sf'][:, 0, SRC_NAMES.index('knockdown_expression_gate')]*f['sf'][:, 0, SRC_NAMES.index('knockdown_detection_gate')]
    return base+model['estimator'].predict(np.column_stack([z*gate[:, None], gate]))


def lobo_select(contexts, cfg, genes, seed, max_kd, per_kd, cap, row_options):
    grid = cfg.get('transfer_grid', [DEFAULT_KNOBS]); alphas = cfg.get('ridge_alphas', [10, 100, 1000, 10000])
    backgrounds = sorted({c.background for c in contexts})
    if len(backgrounds) < 3: raise ValueError('LOBO tuning needs >=3 available backgrounds after exclusions')
    reports = []
    for heldout in backgrounds:
        training = [c for c in contexts if c.background != heldout]
        validation = choose_sources([c for c in contexts if c.background == heldout])
        kd = training_targets(training, seed, max_kd)
        generics = generic_signatures(training, kd, genes, cap)
        scales = fit_source_scales(training, kd, genes, cap, seed, per_kd)
        for index, knobs in enumerate(grid):
            f, y, w = build_rows(training, kd, genes, generics, seed, per_kd, 0.25, cap, knobs=knobs, scales=scales, **row_options)
            fv, yv, wv = build_rows(training, kd, genes, generics, seed, per_kd, 0, cap,
                                     knobs=knobs, scales=scales, label_contexts=validation, **row_options)
            base = pooled(fv)[1]
            baseline_mae = float(np.average(np.abs(base-yv), weights=wv))
            for alpha in alphas:
                model = fit_ridge(f, y, w, alpha)
                # All learned quantities come from training; validation supplies labels only here.
                ridge_mae = float(np.average(np.abs(predict(model, fv)-yv), weights=wv))
                reports.append(dict(background=heldout, knobs_index=index, ridge_alpha=float(alpha),
                                    source_mean_mae=baseline_mae, ridge_mae=ridge_mae, n_rows=len(yv)))
        print(f'LOBO tuning: {heldout} complete', flush=True)
    frame = pd.DataFrame(reports)
    # Knobs are selected for the robust baseline first; ridge then competes against it.
    chosen_knobs = int(frame.groupby('knobs_index').source_mean_mae.mean().idxmin())
    selected = frame[frame.knobs_index == chosen_knobs]
    alpha = float(selected.groupby('ridge_alpha').ridge_mae.mean().idxmin())
    folds = selected[selected.ridge_alpha == alpha].copy()
    anchors = cfg.get('selection_backgrounds', ['CD4_T', 'HCT116', 'HEK293T', 'H1'])
    original = folds[folds.background.isin(anchors)]
    required = max(3, int(np.ceil(0.75*len(original))))
    wins = int((original.ridge_mae < original.source_mean_mae).sum())
    approved = len(original) >= 3 and wins >= required
    return dict(knobs=grid[chosen_knobs], ridge_alpha=alpha, selected_method='model' if approved else 'source_mean',
                original_background_wins=wins, required_wins=required, original_folds=len(original),
                folds=folds.to_dict('records'), history=reports)


def train(config_path, output, holdout=None, kind='ridge', allow_legacy=False, seed=0,
          max_kd=3000, per_kd=128, **overrides):
    cfg = json.loads(Path(config_path).read_text()); genes = read_names(cfg['genes'], 'gene_name')
    if holdout and holdout not in {c['background'] for c in cfg['contexts']}: raise ValueError('unknown holdout background')
    if kind not in ['ridge', 'source_mean']: raise ValueError('reviewed workflow supports ridge/source_mean with LOBO selection')
    if overrides.get('signal_fraction') not in [None, 0]: raise ValueError('use uniform sampling (signal_fraction=0) to avoid unweighted enrichment')
    mode = overrides.get('baseline_mode') or cfg.get('baseline_mode', 'bulk')
    contexts = load_contexts(cfg, genes, [holdout] if holdout else [], allow_legacy, mode)
    cap = float(cfg.get('effect_cap', 8)); row_options = {}
    for key, default in [('line_features', True), ('label_noise_floor', 0.25), ('label_weight_cap', 10), ('label_weighting', True)]:
        row_options[key] = cfg.get(key, default) if overrides.get(key) is None else overrides[key]
    selection = lobo_select(contexts, cfg, genes, seed, max_kd, per_kd, cap, row_options)
    kd = training_targets(contexts, seed, max_kd)
    generics = generic_signatures(contexts, kd, genes, cap); scales = fit_source_scales(contexts, kd, genes, cap, seed, per_kd)
    f, y, w = build_rows(contexts, kd, genes, generics, seed, per_kd, 0.25, cap,
                         knobs=selection['knobs'], scales=scales, **row_options)
    model = fit_ridge(f, y, w, selection['ridge_alpha']) if kind == 'ridge' else dict(kind='source_mean')
    if kind == 'source_mean': selection['selected_method'] = 'source_mean'
    artifact = dict(review_revision=3, model=model, genes=genes, contexts=[c.spec for c in contexts],
                    excluded_background=holdout, generics=generics, cap=cap, seed=seed, allow_legacy=allow_legacy,
                    max_trained_sources=int(f['mask'].sum(1).max()), line_features=row_options['line_features'],
                    baseline_mode=mode, scales=scales, knobs=selection['knobs'], train_targets=kd.tolist(),
                    selected_method=selection['selected_method'], neighbours=cfg.get('neighbours'),
                    neighbour_k=int(cfg.get('neighbour_k', 5)), feature_names=dict(source=SRC_NAMES, target=TGT_NAMES))
    out = Path(output); out.mkdir(parents=True, exist_ok=True)
    joblib.dump(artifact, out/'model.joblib', compress=3)
    save_json(out/'selection.json', {k: v for k, v in selection.items() if k != 'history'})
    save_json(out/'training_report.json', dict(selection=selection, rows=len(y), training_targets=len(kd),
               training_backgrounds=sorted({c.background for c in contexts}), heldout_background=holdout,
               source_scales=scales, prior_variances={c.name: c.effects.prior_variance for c in contexts},
               shrinkage='normal prior estimated per knockdown row (floor 0.01); table prior reported for reference',
               sampling='uniform; all finite baseline genes including low-expression responders',
               model_selection='nested LOBO only; all targets refit after tuning',
               labels={c.name: c.effects.meta for c in contexts}))
    print(f"LOBO selection: {selection['selected_method']}; ridge wins {selection['original_background_wins']}/{selection['original_folds']} original folds")
    return artifact


class FunctionalNeighbours:
    def __init__(self, path, k=5):
        self.k = k; self.adj = {}
        if k < 1: raise ValueError('neighbour_k must be positive')
        if path:
            df = pd.read_csv(path, sep=None, engine='python')
            if not {'gene_a', 'gene_b', 'score'} <= set(df): raise ValueError('neighbours need gene_a,gene_b,score columns')
            for a, b, score in df[['gene_a', 'gene_b', 'score']].itertuples(index=False, name=None):
                if not np.isfinite(score) or score <= 0 or a == b: continue
                for x, y in [(str(a), str(b)), (str(b), str(a))]:
                    self.adj.setdefault(x, {})[y] = max(float(score), self.adj.get(x, {}).get(y, 0))

    def available(self, target, covered):
        pairs = [(g, score) for g, score in self.adj.get(target, {}).items() if g in covered]
        return sorted(pairs, key=lambda x: (-x[1], x[0]))[:self.k]


class BorrowedEffects:
    """Effect table that fills uncovered targets from functional neighbours.

    Each neighbour's row still contains that neighbour's OWN on-target knockdown. That
    column is not a downstream effect of the borrowed target, and it is scored for it,
    so it is masked before averaging; the remaining neighbours fill that gene in.
    """
    def __init__(self, original, graph, genes):
        self.original = original; self.graph = graph
        self.gene_pos = {str(g): i for i, g in enumerate(genes)}
        self.pos = original.pos; self.targets = original.targets; self.meta = original.meta; self.legacy = original.legacy
        self.prior_variance = original.prior_variance

    def neighbours(self, kd): return self.graph.available(kd, self.original.pos)

    def prior_for(self, kd, floor=0.01):
        """Measured knockdowns keep their own row prior; a borrowed profile gets one from its averaged row."""
        if kd in self.pos: return self.original.prior_for(kd, floor)
        cache = self.__dict__.setdefault('_row_priors', {})
        if kd not in cache:
            l, se = self.get(kd, np.arange(self.original.n_genes))
            cache[kd] = estimate_prior(l, se, floor=floor) if np.isfinite(se).any() else self.prior_variance
        return cache[kd]

    def get(self, kd, gi):
        if kd in self.pos: return self.original.get(kd, gi)
        neighbours = self.neighbours(kd)
        if not neighbours: return self.original.get(kd, gi)
        pairs = [self.original.get(g, gi) for g, _ in neighbours]
        l, se = np.stack([x[0] for x in pairs]), np.stack([x[1] for x in pairs])
        columns = np.asarray(gi, dtype=np.int64).ravel()
        for row, (g, _) in enumerate(neighbours):
            own = self.gene_pos.get(g)
            if own is not None: l[row, columns == own] = np.nan
        w = np.asarray([score for _, score in neighbours])[:, None]*np.isfinite(l)
        den = w.sum(0); w /= np.maximum(den, 1e-30)
        avg = (np.nan_to_num(l)*w).sum(0)
        # Conservative mixture uncertainty includes disagreement; unknown donor SE stays missing.
        variance = (w*(np.nan_to_num(se)**2+(np.nan_to_num(l)-avg)**2)).sum(0)
        valid = den > 0; avg[~valid] = np.nan
        result_se = np.sqrt(variance); result_se[~valid | ((w > 0) & ~valid_se(se)).any(0)] = np.nan
        return avg.astype(np.float32), result_se.astype(np.float32)


def holdout_panels(config_path, output, max_targets=1000, seed=0):
    """Per-context diagnostic panels for the LOBO holdout folds.

    A fold's model never sees the held-out background, so prediction raises on any target
    no OTHER background measured. Keep only those directly covered targets, and subsample
    genome-wide panels (seeded per context) so each fold stays a few hundred MB.
    """
    cfg = json.loads(Path(config_path).read_text())
    measured = {c['name']: set(np.load(Path(c['effects'])/'targets.npy').astype(str)) for c in cfg['contexts']}
    out = Path(output); out.mkdir(parents=True, exist_ok=True); summary = {}
    for c in cfg['contexts']:
        others = set().union(*(measured[o['name']] for o in cfg['contexts'] if o['background'] != c['background']))
        covered = sorted(measured[c['name']] & others); n_covered = len(covered)
        if max_targets and n_covered > max_targets:
            rng = np.random.default_rng(seed_for(seed, 'panel', c['name']))
            covered = sorted(rng.choice(covered, max_targets, replace=False).tolist())
        if not covered: print(f"WARNING {c['name']}: no targets covered by other backgrounds; its fold cannot be evaluated")
        pd.DataFrame({'target_gene': covered}).to_csv(out/f"{c['name']}.csv", index=False)
        summary[c['name']] = dict(background=c['background'], retained_targets=len(measured[c['name']]),
                                  covered_by_other_backgrounds=n_covered, panel=len(covered))
    save_json(out/'panels_summary.json', dict(max_targets=max_targets, seed=seed, contexts=summary))
    print(json.dumps(summary, indent=1))
    return summary


def coverage(config_path, targets_path, output, profiles=None):
    cfg = json.loads(Path(config_path).read_text()); requested = read_names(targets_path, 'target_gene')
    graph = FunctionalNeighbours(cfg.get('neighbours'), int(cfg.get('neighbour_k', 5)))
    sources = {s['name']: set(np.load(Path(s['effects'])/'targets.npy').astype(str)) for s in cfg['contexts']}
    union = set().union(*sources.values()); records = []
    for kd in requested:
        direct = [name for name, covered in sources.items() if kd in covered]
        neighbours = graph.available(kd, union)
        records.append(dict(target_gene=kd, direct=bool(direct), sources=';'.join(direct),
                            borrowed=bool(not direct and neighbours), neighbours=';'.join(g for g, _ in neighbours)))
    frame = pd.DataFrame(records); out = Path(output); out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out/'coverage.csv', index=False)
    summary = dict(n_targets=len(frame), direct_fraction=float(frame.direct.mean()),
                   resolved_fraction=float((frame.direct | frame.borrowed).mean()),
                   unresolved=frame.loc[~frame.direct & ~frame.borrowed, 'target_gene'].tolist(),
                   by_source={name: int(sum(kd in s for kd in requested)) for name, s in sources.items()},
                   note='direct coverage after guide QC; optional profiles add per-line measured target expression')
    if profiles:
        genes = read_names(cfg['genes'], 'gene_name'); gp = {g: i for i, g in enumerate(genes)}
        summary['lines'] = {}
        for spec in profiles:
            line, path = spec.split('=', 1); p = load_profile(path, genes)
            summary['lines'][line] = dict(n_targets=len(requested), direct_fraction=summary['direct_fraction'],
                 measured_target_fraction=float(np.mean([k in gp and np.isfinite(p['mean_cpm'][gp[k]]) for k in requested])))
    save_json(out/'coverage_summary.json', summary); print(json.dumps(summary, indent=1))


def infer(model_dir, profile_path, targets_csv, output, platform='unknown', chunk=2048, allow_unsupported=False):
    a = joblib.load(Path(model_dir)/'model.joblib')
    if a.get('review_revision') != 3 or a.get('feature_names') != dict(source=SRC_NAMES, target=TGT_NAMES): raise ValueError('retrain: prior checkpoint schema is incompatible')
    genes = a['genes']; targets = read_names(targets_csv, 'target_gene'); profile = transfer_profile(load_profile(profile_path, genes), a['baseline_mode'])
    contexts = [Context(s, genes, a['allow_legacy'], a['baseline_mode']) for s in a['contexts']]
    graph = FunctionalNeighbours(a.get('neighbours'), a.get('neighbour_k', 5))
    for c in contexts: c.effects = BorrowedEffects(c.effects, graph, genes)
    usable_for = {}
    for kd in targets:
        direct = [c for c in contexts if kd in c.effects.pos]
        candidates = direct or [c for c in contexts if c.effects.neighbours(kd)]
        usable_for[kd] = choose_sources(candidates)
    unresolved = [kd for kd in targets if not usable_for[kd]]
    if unresolved and not allow_unsupported:
        raise ValueError(f'{len(unresolved)} unsupported targets: {unresolved[:20]}. Load Replogle or configure functional neighbours; explicit --allow-unsupported permits documented zeros.')
    out = Path(output)
    if out.exists() and any(out.iterdir()): raise ValueError('choose an empty prediction directory')
    out.mkdir(parents=True, exist_ok=True); np.save(out/'genes.npy', genes); np.save(out/'targets.npy', targets)
    from numpy.lib.format import open_memmap
    arrays = {name: open_memmap(out/f'{name}.npy', mode='w+', dtype='float32', shape=(len(targets), len(genes)))
              for name in ['model', 'zero', 'source_mean', 'generic', 'selected']}
    counts = open_memmap(out/'source_count.npy', mode='w+', dtype='uint16', shape=(len(targets), len(genes)))
    diagnostics = []; ontarget = np.full(len(targets), np.nan, np.float32); gp = {g: i for i, g in enumerate(genes)}
    for i, kd in enumerate(targets):
        usable = usable_for[kd]; ncovered = 0
        for start in range(0, len(genes), chunk):
            gi = np.arange(start, min(start+chunk, len(genes)))
            if not usable:
                for arr in arrays.values(): arr[i, gi] = 0
                counts[i, gi] = 0; continue
            f = features(profile, platform, usable, kd, gi, genes, a['generics'], a['cap'], a['line_features'], a['knobs'], a['scales'])
            number = f['mask'].sum(1); base = pooled(f)[1]
            model = predict(a['model'], f)
            values = dict(zero=np.zeros(len(gi)), source_mean=base, generic=f['generic'], model=model)
            values['selected'] = values[a['selected_method']]
            for name, v in values.items(): arrays[name][i, gi] = np.where(number > 0, np.clip(v, -a['cap'], a['cap']), 0)
            counts[i, gi] = number; ncovered += int((number > 0).sum())
        # Direct successful source on-target effects, similarity-weighted; borrowers do not invent guide efficiency.
        ki = gp.get(kd); own_sources = [c for c in usable if kd in c.effects.pos]
        if ki is not None and own_sources:
            values = np.array([c.effects.get(kd, [ki])[0][0] for c in own_sources])
            corr = np.array([centred_similarity(profile['baseline_log'], c.profile['baseline_log']) for c in own_sources])
            good = np.isfinite(values) & (values < -0.1)
            if good.any():
                weights = np.exp((corr[good]-corr[good].max())/a['knobs']['temperature'])
                ontarget[i] = np.average(values[good], weights=weights)
        diagnostics.append(dict(target=kd, covered_genes=ncovered, source_names=[c.name for c in usable],
                    borrowed={c.name: c.effects.neighbours(kd) for c in usable if kd not in c.effects.pos}))
    for arr in arrays.values(): arr.flush()
    counts.flush(); np.save(out/'ontarget_lfc.npy', ontarget)
    save_json(out/'prediction_report.json', dict(selected_method=a['selected_method'], heldout_background=a['excluded_background'],
             model_dir=str(Path(model_dir).resolve()), target_profile=str(Path(profile_path).resolve()),
             source_names=[c.name for c in contexts], coverage=diagnostics, unresolved=unresolved,
             knobs=a['knobs'], source_scales=a['scales'], allow_unsupported=allow_unsupported))


def evaluate_tables(pred_dir, truth_dir, profile_path, output, allow_legacy=False):
    """LOBO effect MAE diagnostics only. Cell generation is scored with official cell-eval."""
    p = Path(pred_dir); genes = np.load(p/'genes.npy'); targets = np.load(p/'targets.npy')
    truth = Effects(truth_dir, genes, allow_legacy); profile = load_profile(profile_path, genes)
    arrays = {name: np.load(p/f'{name}.npy', mmap_mode='r') for name in ['model', 'source_mean', 'selected', 'zero']}
    records = []
    for i, kd in enumerate(targets):
        y, _ = truth.get(kd, np.arange(len(genes)))
        # Include induced genes as well as baseline-expressed genes.
        valid = np.isfinite(y) & np.isfinite(profile['mean_cpm']) & (genes != kd)
        for name, arr in arrays.items():
            if valid.any(): records.append(dict(target=kd, method=name, mae=float(np.abs(arr[i, valid]-y[valid]).mean()),
                                               mse=float(((arr[i, valid]-y[valid])**2).mean())))
    frame = pd.DataFrame(records)
    if frame.empty: raise ValueError('no truth-aligned targets')
    out = Path(output); out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out/'per_knockdown.csv', index=False); frame.groupby('method')[['mae', 'mse']].mean().to_csv(out/'summary.csv')
    save_json(out/'limitations.json', dict(note='Effect diagnostics only. Official generator scoring uses cell-eval --ceiling.'))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='cmd', required=True)
    t = sub.add_parser('train')
    for key in ['config', 'out']: t.add_argument('--'+key, required=True)
    t.add_argument('--holdout'); t.add_argument('--kind', choices=['ridge', 'source_mean'], default='ridge')
    t.add_argument('--allow-legacy-effects', action='store_true'); t.add_argument('--seed', type=int, default=0)
    t.add_argument('--max-knockdowns', type=int, default=3000); t.add_argument('--genes-per-knockdown', type=int, default=128)
    t.add_argument('--line-features', action=argparse.BooleanOptionalAction, default=None)
    t.add_argument('--label-weighting', action=argparse.BooleanOptionalAction, default=None)
    t.add_argument('--label-noise-floor', type=float); t.add_argument('--label-weight-cap', type=float)
    t.add_argument('--signal-fraction', type=float, default=0.0)
    t.add_argument('--baseline-mode', choices=['bulk', 'per_cell'])
    i = sub.add_parser('predict')
    for key in ['model', 'profile', 'targets', 'out']: i.add_argument('--'+key, required=True)
    i.add_argument('--platform', default='unknown'); i.add_argument('--allow-unsupported', action='store_true')
    pn = sub.add_parser('panels')
    for key in ['config', 'out']: pn.add_argument('--'+key, required=True)
    pn.add_argument('--max-targets', type=int, default=1000); pn.add_argument('--seed', type=int, default=0)
    c = sub.add_parser('coverage')
    for key in ['config', 'targets', 'out']: c.add_argument('--'+key, required=True)
    c.add_argument('--profile', action='append')
    e = sub.add_parser('evaluate-tables')
    for key in ['predictions', 'truth', 'profile', 'out']: e.add_argument('--'+key, required=True)
    e.add_argument('--allow-legacy-effects', action='store_true')
    a = p.parse_args()
    if a.cmd == 'train':
        train(a.config, a.out, a.holdout, a.kind, a.allow_legacy_effects, a.seed, a.max_knockdowns, a.genes_per_knockdown,
              line_features=a.line_features, label_weighting=a.label_weighting, label_noise_floor=a.label_noise_floor,
              label_weight_cap=a.label_weight_cap, baseline_mode=a.baseline_mode, signal_fraction=a.signal_fraction)
    elif a.cmd == 'predict': infer(a.model, a.profile, a.targets, a.out, a.platform, allow_unsupported=a.allow_unsupported)
    elif a.cmd == 'coverage': coverage(a.config, a.targets, a.out, a.profile)
    elif a.cmd == 'panels': holdout_panels(a.config, a.out, a.max_targets, a.seed)
    else: evaluate_tables(a.predictions, a.truth, a.profile, a.out, a.allow_legacy_effects)

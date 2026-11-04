"""Scorer-space MSE study (Jose 2026-10-28; results/stage1/2026-10-28_scorer_space/), on vcc2027.

Every MSE number goes through analysis/mse_metric/mse_study.py's `score_generated(h, P, jk_pred)` (cell_eval2's own
_across_pert_budget / _row_weights, the #247 cap and the #348 budget), with `h` an eval in scorer space:

  EvalPB   one eval's real side from its harness cache (real_cache/*.moments.npz: per-group bulk_lognorm means and
           delete-1 jackknife, built by cell-eval2), its anchors from the bundle (cell_eval2.score._replicate_entries,
           scored with cell_eval2.scoring.score_one, as scripts/v5/pb.py does), PDS from
           cell_eval2.metrics.discrimination.discrimination_score (pb.py's pds_rows, verbatim).

  step1    ~/venv/bin/python code/vcc6/scorer_study.py step1      reproduce the harness raw MSE of HCT116 C+G and Eval A C+G
  step2    ~/venv/bin/python code/vcc6/scorer_study.py step2      Delta per source (no model change): full knockdown
           profile - that source's control profile (log1p_bulk_50k_full), written to OUT/delta/<table>.npz for every
           knockdown any fold or D/E/F needs; the existing layers are untouched. Writes tables/step2_*.
  step4    ~/venv/bin/python code/vcc6/scorer_study.py step4      noise-credit accounting on the generated C+G cells (Eval A,
           HCT116): C_pred, C_real, C0 (control jackknife at n = 400, 50 bootstrap draws from the eval's scoring controls),
           the #247 cap, Q vs B (#348), credit lost, and the raw MSE with the lost credit restored.
  folds    ~/venv/bin/python code/vcc6/scorer_study.py folds      one table load; per fold and line, C+G exactly as
           run_line builds it (stage1 arm C + gaussian on the eval's context controls, the fold's lines' tables
           dropped) and Delta pooled with C+G's weights (same prior, SE curve, votes, shrink; each source's LFC
           replaced by its Delta). Writes cache/fold_<line>.npz.
  arcw     ~/venv/bin/python code/vcc6/scorer_study.py arcw       Step 7's Arc-first sub-arm: on the K562 and HCT116+HEK293T
           folds (H1 is a source there), Arc H1's pooling weight x w_arc in {2, 5, 20} for both pools (C+G's LFC and the
           Delta pool): the H1 background vote's variance / w_arc inside the same random-effects pool. Writes
           cache/fold_<line>_arcw<w>.npz (w = 1 is cache/fold_<line>.npz).
"""
import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
for p in (V / 'analysis/mse_metric', V / 'scripts/v5', V / 'code/vcc6', V / 'code/vcc_2026_3c'):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
OUT = V / 'results/stage1/2026-10-28_scorer_space'
CTRL, BTS = 'non-targeting', 50000.0
EVALS = {'H1': V / 'lib/harness/evalA', 'HCT116': V / 'results/stage1/final_evals/harness/lib/evalHCT116final'}
RUNS = {'H1': V / 'results/stage1/2026-10-27_incumbent/evalA/CG/e2e/salt_none',
        'HCT116': V / 'results/stage1/e2e_final/CG/HCT116/salt_none'}


class EvalPB:
    """The attribute set mse_study.score_generated reads (pb.py's H1 harness, generalised to any eval)."""
    MSE, PDS = 'expr_mse_unbiased_capped_norm', 'pds_cosine'

    def __init__(self, eval_dir, bundle=True):
        from cell_eval2.scoring import score_one
        self._score_one = score_one
        d = Path(eval_dir); self.dir = d
        z = np.load(glob.glob(str(d / 'real_cache' / '*.moments.npz'))[0], allow_pickle=True)
        self.perts = z['perts'].astype(str); self.means = z['means'].astype(np.float64)
        jk = z['jk'].astype(np.float64); idx = {p: i for i, p in enumerate(self.perts)}
        self.genes = np.array(pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str))
        G = self.G = len(self.genes); gp = {g: j for j, g in enumerate(self.genes)}
        self.targets = sorted(p for p in self.perts if p != CTRL); T = self.T = len(self.targets); self.rows = np.arange(T)
        r = [idx[t] for t in self.targets]
        self.R = self.means[r]; self.Creal = jk[r]; self.n_real = z['counts'][r]
        self.ctrl = self.means[idx[CTRL]]; self.Cctrl = float(jk[idx[CTRL]])
        self.own = np.array([gp.get(t, -1) for t in self.targets]); self.nrow = np.where(self.own >= 0, G - 1, G).astype(float)
        dr = self.R - self.ctrl[None]
        self.den = (np.einsum('ij,ij->i', dr, dr) - self.Creal - self.Cctrl
                    - self.own_sq(self.R, np.broadcast_to(self.ctrl, self.R.shape))) / self.nrow
        if bundle:
            from cell_eval2.score import _replicate_entries
            anchor = pd.read_parquet(d / 'bundle/anchor_agg.parquet')
            base = pd.read_csv(d / 'bundle/baseline_agg.csv').set_index('statistic').loc['mean'].to_dict()
            self.entries = _replicate_entries(base, anchor); rep = anchor.set_index('metric')['replicate'].to_dict()
            self.bm, self.rm = float(base[self.MSE]), float(rep[self.MSE])
            self.bp, self.rp = float(base[self.PDS]), float(rep[self.PDS])

    def own_sq(self, A, B):
        return np.array([(A[i, self.own[i]] - B[i, self.own[i]]) ** 2 if self.own[i] >= 0 else 0.0 for i in self.rows])

    def scale(self, metric, u):
        e = self.entries[metric]; return float(self._score_one(float(u), e.base, e.scoring))

    def pds_rows(self, P):
        from cell_eval2.metrics.discrimination import discrimination_score
        out = discrimination_score(pred_bulk=(np.array(self.targets), P), real_bulk=(self.perts, self.means), distance='cosine',
                                   rank_denominator='n-1', tie_policy='midrank', exclude_target_gene=True,
                                   exclusion_scope='panel', control_source='real', genes=list(self.genes), control=CTRL)
        return np.array([out[t] for t in self.targets])


def pred_cells(pred_h5ad, targets):
    """Predicted bulk_lognorm means and jackknife from generated cells (cell_eval2's prep + moments)."""
    import anndata as ad
    from cell_eval2.moments import correction_for
    from cell_eval2.prep import pseudobulk_bulk_lognorm_with_moments
    a = ad.read_h5ad(pred_h5ad)
    col = 'target_gene' if 'target_gene' in a.obs else 'target'
    perts, means, mom = pseudobulk_bulk_lognorm_with_moments(a, col, bulk_target_sum=BTS)
    C = np.asarray(correction_for(mom, perts, means, comparator='bulk_lognorm'), np.float64)
    idx = {str(p): i for i, p in enumerate(perts)}
    rows = [idx[t] for t in targets]
    return np.asarray(means, np.float64)[rows], C[rows]


def step1(a):
    import mse_study as M
    res = {}
    for line in ('HCT116', 'H1'):
        h = EvalPB(EVALS[line])
        P, C = pred_cells(RUNS[line] / 'prediction.h5ad', h.targets)
        np.savez_compressed(OUT / f'tables/pred_cells_CG_{line}.npz', targets=np.array(h.targets), P=P, C=C)
        s, _ = M.score_generated(h, P, C)
        ev = 'evalHCT116final' if line == 'HCT116' else 'evalA'
        harness = json.loads((RUNS[line] / 'harness' / ev / 'summary.json').read_text())
        res[line] = dict(knockdowns=h.T, mse_study=s, harness_raw=harness['raw'][h.MSE],
                         harness_pds_raw=harness['raw'][h.PDS], abs_diff=abs(s['mse_raw'] - harness['raw'][h.MSE]),
                         reproduces_1e6=bool(abs(s['mse_raw'] - harness['raw'][h.MSE]) < 1e-6),
                         pds_abs_diff=abs(s['pds_raw'] - harness['raw'][h.PDS]), baseline_raw=h.bm, replicate_raw=h.rm)
        print(line, json.dumps(res[line], default=float), flush=True)
    (OUT / 'tables/step1_reproduce.json').write_text(json.dumps(res, indent=1, default=float) + '\n')


def control_C0(eval_dir, n=400, draws=50, seed=2026):
    """Delete-1 jackknife (cell_eval2.moments.jackknife_correction) of n control cells drawn from the eval's scoring
    controls (real.h5ad), `draws` times."""
    import anndata as ad
    from cell_eval2.moments import jackknife_correction
    a = ad.read_h5ad(Path(eval_dir) / 'real.h5ad', backed='r')
    col = 'target_gene' if 'target_gene' in a.obs else 'target'
    ctrl = np.flatnonzero(a.obs[col].astype(str).to_numpy() == CTRL)
    X = a[ctrl].to_memory().X
    X = X.tocsr() if hasattr(X, 'tocsr') else X
    rng = np.random.default_rng(seed); out = []
    for _ in range(draws):
        pick = np.sort(rng.choice(X.shape[0], n, replace=False))
        out.append(float(np.asarray(jackknife_correction(X[pick].astype(np.float64), np.zeros(n, np.int64), 1, BTS)).ravel()[0]))
    return np.array(out), int(len(ctrl))


def step4(a):
    import mse_study as M
    from cell_eval2.metrics.delta import _across_pert_budget, _row_weights
    res, per = {}, []
    for line in ('H1', 'HCT116'):
        h = EvalPB(EVALS[line])
        z = np.load(OUT / f'tables/pred_cells_CG_{line}.npz', allow_pickle=True)
        assert list(z['targets']) == h.targets
        P, C = z['P'], z['C']
        C0s, n_ctrl = control_C0(EVALS[line])
        C0 = float(C0s.mean())
        excl = {i: int(h.own[i]) for i in h.rows if h.own[i] >= 0}
        w = _row_weights(h.T, h.G, excl)
        capped = np.minimum(C, h.Creal)
        Q = float((capped * w).sum()); Qfull = float((C * w).sum())
        B = _across_pert_budget(P, excl, w)
        actual, _ = M.score_generated(h, P, C)
        # all lost credit restored: the #348 part (Q - min(B, Q)) and, separately, the #247 cap too
        d = P - h.R; plug = np.einsum('ij,ij->i', d, d) - h.own_sq(P, h.R)
        raw_budget_restored = float(((plug - capped - h.Creal) / h.nrow).sum() / h.den.sum())
        raw_all_restored = float(((plug - C - h.Creal) / h.nrow).sum() / h.den.sum())
        res[line] = dict(knockdowns=h.T, C0_mean=C0, C0_sd=float(C0s.std()), C0_draws=len(C0s), scoring_controls=n_ctrl,
                         C_pred_median=float(np.median(C)), C_real_median=float(np.median(h.Creal)),
                         share_Cpred_gt_Creal=float(np.mean(C > h.Creal)),
                         ratio_Cpred_over_C0_median=float(np.median(C / C0)),
                         ratio_Cpred_over_C0_q=np.quantile(C / C0, [0.1, 0.25, 0.75, 0.9]).tolist(),
                         Q_claim_after_cap=Q, claim_before_cap=Qfull, B=B, budget_binds=bool(B < Q),
                         credit_lost_to_cap=Qfull - Q, credit_lost_to_budget=max(0.0, Q - B),
                         credit_granted=min(B, Q), den_sum=float(h.den.sum()),
                         raw_actual=actual['mse_raw'], raw_budget_loss_restored=raw_budget_restored,
                         raw_cap_and_budget_restored=raw_all_restored,
                         diff_budget_restored=raw_budget_restored - actual['mse_raw'],
                         diff_all_restored=raw_all_restored - actual['mse_raw'])
        per.append(pd.DataFrame(dict(line=line, target=h.targets, C_pred=C, C_real=h.Creal, C0=C0, ratio=C / C0,
                                     n_real=h.n_real)))
        print(line, json.dumps(res[line], default=float), flush=True)
    pd.concat(per).to_parquet(OUT / 'tables/step4_per_knockdown.parquet')
    (OUT / 'tables/step4.json').write_text(json.dumps(res, indent=1, default=float) + '\n')


EVAL_DIRS = dict(EVALS, K562=V / 'results/stage1/final_evals/harness/lib/evalK562final',
                 HEK293T=V / 'results/stage1/final_evals/harness/lib/evalHEK293Tfinal')


F1_DROP = ('tahoe', 'jiang', 'marson', 'adamson')
EFF_MAX = float(np.log2(0.2))                                   # F2: on-target knockdown >= 80%


def table_kept(name, modality, reference, filt):
    """Step 10's source filters. F0 every table. F1 single-cell genetic CRISPRi only: Tahoe, Jiang, Marson, Adamson and
    the rest-of-line (Feng GW) tables out; also anything not single-cell (the Orion pseudobulk complements). F2 = F1."""
    if filt == 'F0':
        return True
    return (name.split('_')[0] not in F1_DROP and reference != 'rest_of_line' and modality == 'single_cell')


def row_mask(path, filt):
    """Table rows a filter keeps (None = all). F2: efficiency known and <= log2(0.2); unknown efficiency dropped."""
    if filt != 'F2':
        return None
    import h5py
    from anndata.io import read_elem
    with h5py.File(path, 'r') as fh:
        o = read_elem(fh['obs'])
    return (o['efficiency_known'].astype(bool) & (o['efficiency'].astype(float) <= EFF_MAX)).to_numpy()


def load_source_rows(path, mask):
    """stage1.load_source on a subset of the table's rows (F2), duplicates collapsed exactly as it collapses them."""
    import anndata as ad
    import stage1 as S
    a = ad.read_h5ad(path)
    try:
        sub = np.flatnonzero(mask)
        modality, reference, background = S.table_meta(a.obs, a.uns, path.stem)
        n = (a.obs['n_cells_used'].to_numpy(np.float64)[sub] if 'n_cells_used' in a.obs.columns else np.full(len(sub), np.nan))
        src = S.Source(name=path.stem, background=background, modality=modality, reference=reference,
                       genes=S._as_str_array(a.var_names), knockdowns=S._as_str_array(a.obs['target_gene'].iloc[sub].tolist()),
                       lfc=np.asarray(a.layers['lfc'][sub], np.float32), se=np.asarray(a.layers['se'][sub], np.float32),
                       n_cells=n.astype(np.float32), mult=S.multiplier_for(modality, reference))
        if modality == 'drug':
            return src
        kd, lfc, se, n = S._collapse_duplicate_knockdowns(src.knockdowns, src.lfc, src.se, src.n_cells)
        return S.Source(src.name, src.background, src.modality, src.reference, src.genes, kd, lfc, se, n, src.mult)
    finally:
        if getattr(a, 'file', None) is not None:
            a.file.close()


def filter_plan(tables, filt):
    """Which tables (and, for F2, which rows) a filter keeps, from the obs alone; and the knockdowns they hold."""
    import h5py
    import stage1 as S
    from anndata.io import read_elem
    keep, masks, held = [], {}, set()
    for p in tables:
        with h5py.File(p, 'r') as fh:
            obs = read_elem(fh['obs'])
        mod, ref, _ = S.table_meta(obs, {}, p.stem)
        if table_kept(p.stem, mod, ref, filt):
            m = row_mask(p, filt); keep.append(p); masks[p.stem] = m
            kd = obs['target_gene'].astype(str).to_numpy()
            held |= set(kd if m is None else kd[m])
    return keep, masks, held


def load_filtered(keep, masks, filt):
    """The sources a filter keeps: whole tables (F0/F1) or the kept rows (F2)."""
    import stage1 as S
    if filt != 'F2':
        return S.load_tables(keep)
    return [load_source_rows(p, masks[p.stem]) for p in keep if masks[p.stem].any()]


def delta_sources(others, genes, masks=None):
    """Delta copies of the voting sources that have a Delta file: lfc := Delta rows, se := the raw table's se rows,
    duplicates collapsed exactly as stage1.load_source collapses them (drug tables keep every row). masks: per table,
    the rows a source filter keeps (Step 10's F2)."""
    import h5py
    import stage1 as S
    from anndata.io import read_elem
    out = []
    for src in others:
        f_ = OUT / 'delta' / f'{src.name}.npz'
        if not f_.exists():
            continue
        z = np.load(f_, allow_pickle=True); kd = z['knockdowns'].astype(str); D = z['delta']
        path = next(p for p in S.discover_tables([str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')]) if p.stem == src.name)
        with h5py.File(path, 'r') as fh:
            allk = read_elem(fh['obs'])['target_gene'].astype(str).to_numpy()
            need = set(kd); rows = sorted(i for i, k in enumerate(allk) if k in need)
            se = fh['layers/se'][rows]
        assert list(allk[rows]) == list(kd), src.name
        m = (masks or {}).get(src.name)
        if m is not None:
            j = np.array([r for r, i in enumerate(rows) if m[i]], dtype=int)
            if len(j) == 0:
                continue
            kd, D, se = kd[j], D[j], se[j]
        n = np.full(len(kd), np.nan, np.float32)
        if src.modality != 'drug':
            kd, D, se, n = S._collapse_duplicate_knockdowns(kd, D, se, n)
        out.append(S.Source(src.name, src.background, src.modality, src.reference, src.genes, kd,
                            np.asarray(D, np.float32), np.asarray(se, np.float32), n, src.mult))
    return out


def folds(a):
    import stage1 as S
    import screen as SC
    import run_line as RL
    filt = getattr(a, 'filter', 'F0'); suffix = '' if filt == 'F0' else f'_{filt}'
    genes = RL.panel_genes()
    tables = S.discover_tables(SC.DIRS); by_stem = {p.stem: p for p in tables}
    print(f'loading {len(tables)} tables, filter {filt}', flush=True)
    keep, masks, held = filter_plan(tables, filt)
    (OUT / 'cache').mkdir(exist_ok=True)
    defs = list(pd.read_csv(V / 'data_final/extracted/pert_counts.csv')['target_gene'].astype(str))
    cov = dict(filter=filt, tables=[p.stem for p in keep], n_tables=len(keep),
               rows_kept={k: (int(m.sum()), int(m.size)) for k, m in masks.items() if m is not None},
               def_targets_no_source=sorted(set(defs) - held), def_targets_with_source=len(set(defs) & held))
    (OUT / 'tables' / f'step10_coverage_{filt}.json').write_text(json.dumps(cov, indent=1) + '\n')
    print(f'{filt}: {len(keep)} tables; D/E/F targets with a source {cov["def_targets_with_source"]} / {len(defs)}', flush=True)
    if getattr(a, 'coverage_only', False):
        return
    sources = load_filtered(keep, masks, filt)
    for f_, fl in FOLD_LINES.items():
        drop = {LINE_KEY[x] for x in fl}
        others = [s for s in sources if S._line_key(s.background) not in drop]
        basals = {}
        for src in others:
            got = S.load_baseline(by_stem[src.name].with_name(src.name + '_baseline.csv'), src.genes)
            if got is not None: basals[src.name] = (src.genes, got[0])
        votes, prior_votes = S.arm_votes(others, 'C')
        calibration = S.fit_se_calibration(others, basals)
        prior = S.fit_prior(others, votes=prior_votes)
        dsrc = delta_sources(others, genes, masks); dvotes = S.BackgroundVotes(dsrc)
        print(f'fold {f_}: {len(others)} tables vote, {len(dsrc)} with Delta; prior {prior.as_dict()}', flush=True)
        for line in fl:
            e = EvalPB(EVAL_DIRS[line], bundle=False)
            out = OUT / 'cache' / f'controls_{line}'; out.mkdir(parents=True, exist_ok=True)
            X = RL.read_controls(EVAL_DIRS[line] / 'context.h5ad', genes)
            m0, gate, _, _ = RL.control_profile(X, out)
            raw = np.asarray(X.sum(0)).ravel().astype(np.float64)
            ctrl_profile = np.log1p(BTS * raw / raw.sum())
            K, G = len(e.targets), len(genes)
            A = {k: np.zeros((K, G), np.float32) for k in ('lfc', 'se', 'dpool')}
            sig = np.zeros((K, G), bool); ns = np.zeros((K, G), np.int16); dns = np.zeros((K, G), np.int16)
            for i, k in enumerate(e.targets):
                pr = S.predict_knockdown(others, k, np.asarray(genes, dtype=object), m0, gate, calibration, prior,
                                         drop_gene=k, votes=votes, magnitude='gaussian')
                A['lfc'][i] = pr['lfc']; A['se'][i] = pr['se_measurement']; sig[i] = pr['de_sig']; ns[i] = pr['n_sources']
                pd_ = S.predict_knockdown(None, k, np.asarray(genes, dtype=object), m0, gate, calibration, prior,
                                          drop_gene=k, votes=dvotes, magnitude='gaussian')
                A['dpool'][i] = pd_['lfc']; dns[i] = pd_['n_sources']
            np.savez_compressed(OUT / 'cache' / f'fold_{line}{suffix}.npz', targets=np.array(e.targets), genes=np.array(genes),
                                lfc=A['lfc'], de_sig=sig, se_measurement=A['se'], n_sources=ns, dpool=A['dpool'],
                                dpool_n=dns, m0=m0, gate=gate, ctrl_profile=ctrl_profile,
                                prior=json.dumps(prior.as_dict()), tables=json.dumps([s.name for s in others]))
            print(f'  {line}: {K} knockdowns, calls median {np.median(sig.sum(1)):.0f}; Delta measured share '
                  f'{(dns > 0).mean():.3f}', flush=True)


class WeightedVotes:
    """stage1.BackgroundVotes with one background's vote variance divided by w (its pooling weight x w)."""
    def __init__(self, votes, background, w):
        self.v, self.bg, self.w = votes, background, float(w)
        self.groups, self.owners = votes.groups, votes.owners

    def gather(self, knockdown, genes):
        Y, SE, mults = self.v.gather(knockdown, genes)
        bgs = list(self.v.owners.get(str(knockdown), ()))
        SE = [se / np.sqrt(self.w) if bg == self.bg else se for bg, se in zip(bgs, SE)]
        return Y, SE, mults


def arcw(a):
    import stage1 as S
    import screen as SC
    import run_line as RL
    genes = RL.panel_genes(); G_ = np.asarray(genes, dtype=object)
    tables = S.discover_tables(SC.DIRS); by_stem = {p.stem: p for p in tables}
    sources = S.load_tables(tables)
    for f_ in ('K562', 'HCT116+HEK293T'):
        fl = FOLD_LINES[f_]; drop = {LINE_KEY[x] for x in fl}
        others = [s for s in sources if S._line_key(s.background) not in drop]
        basals = {}
        for src in others:
            got = S.load_baseline(by_stem[src.name].with_name(src.name + '_baseline.csv'), src.genes)
            if got is not None: basals[src.name] = (src.genes, got[0])
        votes, prior_votes = S.arm_votes(others, 'C')
        calibration = S.fit_se_calibration(others, basals); prior = S.fit_prior(others, votes=prior_votes)
        dvotes = S.BackgroundVotes(delta_sources(others, genes))
        for line in fl:
            base = np.load(OUT / 'cache' / f'fold_{line}.npz', allow_pickle=True)
            m0, gate = base['m0'], base['gate']; targets = list(base['targets'])
            for w in (2, 5, 20):
                V1, V2 = WeightedVotes(votes, 'h1', w), WeightedVotes(dvotes, 'h1', w)
                K, G = len(targets), len(genes)
                lfc = np.zeros((K, G), np.float32); dp = np.zeros((K, G), np.float32); sig = np.zeros((K, G), bool)
                ns = np.zeros((K, G), np.int16); dn = np.zeros((K, G), np.int16)
                for i, k in enumerate(targets):
                    pr = S.predict_knockdown(None, k, G_, m0, gate, calibration, prior, drop_gene=k, votes=V1, magnitude='gaussian')
                    lfc[i], sig[i], ns[i] = pr['lfc'], pr['de_sig'], pr['n_sources']
                    pd_ = S.predict_knockdown(None, k, G_, m0, gate, calibration, prior, drop_gene=k, votes=V2, magnitude='gaussian')
                    dp[i], dn[i] = pd_['lfc'], pd_['n_sources']
                out = dict(base); out.update(lfc=lfc, de_sig=sig, n_sources=ns, dpool=dp, dpool_n=dn, arc_w=w)
                np.savez_compressed(OUT / 'cache' / f'fold_{line}_arcw{w}.npz', **out)
                print(f'{f_} {line} w_arc {w}: calls median {np.median(sig.sum(1)):.0f}', flush=True)


FOLD_LINES = {'K562': ['K562'], 'HCT116+HEK293T': ['HCT116', 'HEK293T'], 'H1': ['H1']}
LINE_TABLE = {'K562': 'replogle_k562_gw', 'HCT116': 'orion_hct116', 'HEK293T': 'orion_hek293t', 'H1': 'arc_h1'}
LINE_KEY = {'K562': 'k562', 'HCT116': 'hct116', 'HEK293T': 'hek293t', 'H1': 'h1'}


def knockdown_lists():
    L = {l: [k for k in (V / f'results/target_lists/final_{l}.txt').read_text().split() if k] for l in LINE_TABLE}
    L['H1_evalA'] = [k for k in (V / 'results/target_lists/evalA_h1_292.txt').read_text().split() if k]
    L['DEF'] = list(pd.read_csv(V / 'data_final/extracted/pert_counts.csv')['target_gene'].astype(str))
    return L


def step2(a):
    import h5py
    import stage1 as S
    import screen as SC
    from anndata.io import read_elem
    lists = knockdown_lists(); need = set().union(*map(set, lists.values()))
    genes = list(pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str))
    (OUT / 'delta').mkdir(parents=True, exist_ok=True)
    rows, cov = [], {}
    for p in S.discover_tables(SC.DIRS):
        name = p.stem
        base = pd.read_csv(p.with_name(name + '_baseline.csv'), index_col=0)
        ctrl = base['log1p_bulk_50k_full'].reindex(genes).to_numpy(np.float64)
        with h5py.File(p, 'r') as f:
            obs = read_elem(f['obs']); var = np.array(read_elem(f['var']).index.astype(str))
            assert list(var) == genes, f'{name}: gene order'
            kd = obs['target_gene'].astype(str).to_numpy()
            pick = [i for i, k in enumerate(kd) if k in need]
            full = f['layers/log1p_bulk_50k_full'][sorted(pick)] if pick else np.zeros((0, len(genes)))
            lfc_any = np.zeros(len(genes), bool)
            for s0 in range(0, len(kd), 2000):
                lfc_any |= np.isfinite(f['layers/lfc'][s0:s0 + 2000]).any(0)
            mod, ref, bg = S.table_meta(obs, {}, name)
        pick = sorted(pick)
        full = np.asarray(full, np.float64)
        ok_ctrl = bool(np.isfinite(ctrl).any()); ok_full = bool(np.isfinite(full).any()) if full.size else None
        usable = ok_ctrl and bool(ok_full)
        if usable:
            D = (full - ctrl[None, :]).astype(np.float32)
            np.savez_compressed(OUT / 'delta' / f'{name}.npz', knockdowns=kd[pick], genes=np.array(genes), delta=D)
            cov[name] = set(kd[pick])
        u = dict(read_elem(h5py.File(p, 'r')['uns']['effects_v5']))
        rows.append(dict(table=name, background=bg, line_key=S._line_key(bg), modality=mod, reference=ref,
                         rest_of_line=bool(ref == 'rest_of_line'), knockdowns_needed=len(pick),
                         control_full_profile=ok_ctrl, knockdown_full_profile=ok_full, delta_usable=usable,
                         panel_genes_measured_share=float(lfc_any.mean()),
                         chemistry_metadata=str(u.get('depth_note', '') or u.get('method', '') or 'not recorded')[:120],
                         def_targets_held=len(set(kd) & set(lists['DEF']))))
        print(name, usable, len(pick), flush=True)
    t = pd.DataFrame(rows); t.to_parquet(OUT / 'tables/step2_sources.parquet')
    # per fold: (k, g) of the fold lines' knockdowns with no Delta source once the fold's lines are removed
    fold = []
    tables = {r['table']: r for r in rows}
    for f_, fl in FOLD_LINES.items():
        drop = {LINE_KEY[x] for x in fl}
        for line in fl:
            ks = lists['H1_evalA'] if line == 'H1' else lists[line]
            tb = LINE_TABLE[line]
            base = pd.read_csv(next(p for p in S.discover_tables(SC.DIRS) if p.stem == tb).with_name(tb + '_baseline.csv'), index_col=0)
            gate = base['de_tested'].reindex(genes).fillna(False).to_numpy(bool)
            srcs = [n for n in cov if tables[n]['line_key'] not in drop]
            n_no_all = n_no_gate = 0; kd_no = 0
            have = np.zeros((len(ks), len(genes)), bool)
            for n in srcs:
                z = np.load(OUT / 'delta' / f'{n}.npz', allow_pickle=True); pos = {k: i for i, k in enumerate(z['knockdowns'])}
                for i, k in enumerate(ks):
                    if k in pos: have[i] |= np.isfinite(z['delta'][pos[k]])
            fold.append(dict(fold=f_, line=line, knockdowns=len(ks), kg_all=int(have.size), kg_no_delta_all=int((~have).sum()),
                             kg_gated=int(gate.sum()) * len(ks), kg_no_delta_gated=int((~have[:, gate]).sum()),
                             knockdowns_no_delta_at_all=int((~have.any(1)).sum())))
    pd.DataFrame(fold).to_parquet(OUT / 'tables/step2_fold_coverage.parquet')
    # sanity: H1 SIN3B -> SERPINB9
    p = next(p for p in S.discover_tables(SC.DIRS) if p.stem == 'arc_h1')
    with h5py.File(p, 'r') as f:
        obs = read_elem(f['obs']); var = list(read_elem(f['var']).index.astype(str))
        i = list(obs['target_gene'].astype(str)).index('SIN3B'); j = var.index('SERPINB9')
        lfc = float(f['layers/lfc'][i, j]); fk = float(f['layers/log1p_bulk_50k_full'][i, j])
    b = pd.read_csv(p.with_name('arc_h1_baseline.csv'), index_col=0).loc['SERPINB9']
    san = dict(lfc=lfc, lfc_ln2=lfc * np.log(2), control_mean_cpm=float(b['control_mean_cpm']), full_knockdown=fk,
               full_control=float(b['log1p_bulk_50k_full']), delta=fk - float(b['log1p_bulk_50k_full']))
    (OUT / 'tables/step2_sanity.json').write_text(json.dumps(san, indent=1) + '\n')
    print(t.drop(columns=['chemistry_metadata']).to_string()); print(pd.DataFrame(fold).to_string()); print(san)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    for c in ('step1', 'step2', 'step4', 'folds', 'arcw'): sub.add_parser(c)
    sub.choices['folds'].add_argument('--filter', default='F0', choices=['F0', 'F1', 'F2'])
    sub.choices['folds'].add_argument('--coverage-only', action='store_true')
    a = p.parse_args()
    {'step1': step1, 'step2': step2, 'step4': step4, 'folds': folds, 'arcw': arcw}[a.cmd](a)

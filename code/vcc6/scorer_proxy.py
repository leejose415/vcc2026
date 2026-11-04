"""Scorer-space study part 2 (Jose 2026-10-28, results/stage1/2026-10-28_scorer_space/): noise-free proxies for all six
vcc2026 components, so arms can be ranked by expected overall without generating cells.

Per eval (`Ref`), from a requested (LFC, call set) and the target's context-control CPM m0:
  mu       the generator's own mean inversion, (m0 + 1e-9) 2^lfc - 1e-9 (cell_generator.target_means)
  MSE/PDS  P = log1p(5e4 mu / sum mu), C_pred = C0 x 1.1295 (Step 4), through mse_study.score_generated (Steps 3-5)
  DE       a synthetic pred DE table on the real DE table's (target, feature) rows, scored by cell-eval2's own
           prep_de_side -> rank_de_side -> assemble_prepared_de -> run.dispatch_de_metrics (the path cell-eval2 run
           takes), against the eval's cached real DE (real_cache/de_wilcoxon_table-*.parquet):
             log2_fold_change  pdex's log2((mean_pred + 1e-9) / (mean_ctrl + 1e-9)), mean_pred = mu renormalised to
                               1e6 (the scorer's per-cell normalisation), mean_ctrl = the real scoring controls'
                               arithmetic mean CPM (control_source: real)
             p_adj = p_value   requested call -> 1e-3 / (1 + 1e4 s), else 0.5 + 0.5 / (1 + 1e4 s), s = |P - real
                               control bulk| (so the requested set is exactly the significant set, ranked by size)
Scaled values: cell_eval2.scoring.score_one on the bundle's replicate entries (EvalPB.scale); overall = their mean.

  ~/venv/bin/python code/vcc6/scorer_proxy.py step9
"""
import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
import scorer_study as SS                                                                  # noqa: E402
import scorer_arms as SA                                                                   # noqa: E402

OUT = SS.OUT
MSE, PDS = SS.EvalPB.MSE, SS.EvalPB.PDS
DE = ['de_wilcoxon_lfc_nmae', 'de_wilcoxon_direction_fidelity_yield_raw', 'de_wilcoxon_direction_reach_raw',
      'de_wilcoxon_sig_jaccard']
SIX = [MSE, PDS] + DE
SHORT = {MSE: 'mse', PDS: 'pds', DE[0]: 'nmae', DE[1]: 'fidelity', DE[2]: 'reach', DE[3]: 'jaccard'}


def cfg_vcc2026():
    from cell_eval2.config import EvalConfig
    return EvalConfig.from_preset('vcc2026')


class Ref:
    """One eval's real side for all six components."""

    def __init__(self, eval_dir, ratio=None, name=None):
        import polars as pl
        from cell_eval2.de import prep_de_side, rank_de_side, resolve_target_genes
        self.dir = Path(eval_dir); self.name = name or self.dir.name
        self.h = h = SS.EvalPB(eval_dir)
        self.cfg = c = cfg_vcc2026()
        tab = pl.read_parquet(glob.glob(str(self.dir / 'real_cache' / 'de_wilcoxon_table-*.parquet'))[0])
        self.real_df, self.real_perts = prep_de_side(tab, name='real', sort_by=c.de.sort_by,
                                                     nan_lfc_policy=c.de.nan_lfc_policy, min_abs_log2fc=c.de.min_abs_log2fc)
        assert self.real_perts == h.targets, f'{self.name}: real DE targets differ from the moments targets'
        self.real_rank = rank_de_side(self.real_df, sort_by=c.de.sort_by, p_adj_threshold=c.de.p_adj_threshold)
        self.res = resolve_target_genes(self.real_df, self.real_perts, target_gene_map=c.target_gene_map)
        gp = {g: j for j, g in enumerate(h.genes)}; tp = {t: i for i, t in enumerate(h.targets)}
        self.rows_t = self.real_df['target'].to_numpy(); self.rows_f = self.real_df['feature'].to_numpy()
        self.ri = np.array([tp[t] for t in self.rows_t]); self.rg = np.array([gp[f] for f in self.rows_f])
        self.ctrl_cpm = self._ctrl_cpm()
        self.Cpred = None
        if ratio is not None:
            C0s, _ = SS.control_C0(eval_dir); self.C0 = float(C0s.mean()); self.Cpred = np.full(h.T, self.C0 * ratio)

    def _ctrl_cpm(self):
        """Arithmetic mean over the real scoring controls of per-cell CPM (pdex's reference mean), cached."""
        f = OUT / 'cache' / f'realctrl_cpm_{self.name}.npy'
        if f.exists():
            return np.load(f)
        import anndata as ad
        a = ad.read_h5ad(self.dir / 'real.h5ad', backed='r')
        col = 'target_gene' if 'target_gene' in a.obs else 'target'
        ctrl = np.flatnonzero(a.obs[col].astype(str).to_numpy() == SS.CTRL)
        X = a[ctrl].to_memory().X; X = X.tocsr() if hasattr(X, 'tocsr') else X
        tot = np.asarray(X.sum(1)).ravel()
        v = np.asarray(X.multiply(1e6 / tot[:, None]).mean(0)).ravel() if hasattr(X, 'multiply') else (X * (1e6 / tot[:, None])).mean(0)
        assert list(a.var_names.astype(str)) == list(self.h.genes)
        f.parent.mkdir(exist_ok=True); np.save(f, v)
        return v

    # ------------------------------------------------------------------ DE
    def de_metrics(self, pred_df):
        """cell-eval2's four DE components per knockdown (h.targets order) for a pred DE table."""
        from cell_eval2.de import assemble_prepared_de, prep_de_side, rank_de_side
        from cell_eval2.run import dispatch_de_metrics
        c = self.cfg
        pdf, pperts = prep_de_side(pred_df, name='pred', sort_by=c.de.sort_by, nan_lfc_policy=c.de.nan_lfc_policy,
                                   min_abs_log2fc=c.de.min_abs_log2fc)
        prank = rank_de_side(pdf, sort_by=c.de.sort_by, p_adj_threshold=c.de.p_adj_threshold)
        prep = assemble_prepared_de(self.real_rank, self.real_perts, prank, pperts, control=c.control, sort_by=c.de.sort_by,
                                    p_adj_threshold=c.de.p_adj_threshold, real_df=self.real_df, pred_df=pdf,
                                    target_resolution=self.res)
        rows = pd.DataFrame(dispatch_de_metrics(DE, prep, c))
        w = rows.pivot(index='perturbation', columns='metric', values='value').reindex(self.h.targets)
        return {m: w[m].to_numpy(np.float64) for m in DE}

    def pred_table(self, mu, sig, P):
        """The synthetic pred DE table for noise-free means mu (T x G, CPM scale), call set sig and scorer bulk P."""
        import polars as pl
        cpm = mu * (1e6 / mu.sum(1, keepdims=True))
        lfc = np.log2((cpm[self.ri, self.rg] + 1e-9) / (self.ctrl_cpm[self.rg] + 1e-9))
        s = np.abs(P[self.ri, self.rg] - self.h.ctrl[self.rg])
        called = sig[self.ri, self.rg]
        p = np.where(called, 1e-3 / (1 + 1e4 * s), 0.5 + 0.5 / (1 + 1e4 * s))
        return pl.DataFrame(dict(target=self.rows_t, feature=self.rows_f, log2_fold_change=lfc, p_value=p, p_adj=p))

    # ------------------------------------------------------------------ all six
    def proxy(self, lfc=None, sig=None, m0=None, mu=None, per_kd=False):
        """Six raw + scaled components and the overall for a requested prediction (lfc, sig on h.targets x panel,
        m0 = the target's context-control CPM), or directly from noise-free means mu."""
        import mse_study as M
        h = self.h
        if mu is None:
            mu = SA.mu_from_lfc(lfc, m0)
        P = SA.bulk(mu)
        s, num = M.score_generated(h, P, self.Cpred)
        pds = h.pds_rows(P)
        de = self.de_metrics(self.pred_table(mu, sig, P))
        kd = {MSE + ':num': num, MSE + ':den': h.den, PDS: pds, **de}
        return self.summarise(kd, per_kd)

    def summarise(self, kd, per_kd=False, rows=None):
        """Aggregate per-knockdown values as cell-eval2 does (MSE ratio of sums, others means), then scale."""
        rows = np.arange(self.h.T) if rows is None else rows
        raw = {MSE: float(kd[MSE + ':num'][rows].sum() / kd[MSE + ':den'][rows].sum())}
        for m in [PDS] + DE:
            raw[m] = float(np.nanmean(kd[m][rows]))
        sc = {m: self.h.scale(m, raw[m]) for m in SIX}
        out = dict(raw=raw, scaled=sc, overall=float(np.mean([sc[m] for m in SIX])),
                   mse_unclamped=(raw[MSE] - self.h.bm) / (self.h.rm - self.h.bm))
        if per_kd:
            out['per_kd'] = kd
        return out


def requested(npz):
    """(knockdowns, lfc, sig) as the generator sees them (stage1.prepare_generator_inputs: unmeasured genes 0, uncalled)."""
    z = np.load(npz, allow_pickle=True)
    un = z['n_sources'] <= 0
    return list(map(str, z['knockdowns'])), np.where(un, 0.0, z['lfc'].astype(np.float64)), z['de_sig'] & ~un


INC = V / 'results/stage1/2026-10-27_incumbent/evalA'
RUNS9 = {
    'EvalA_CG': dict(line='H1', eval=SS.EVALS['H1'], harness=INC / 'CG/e2e/salt_none/harness/evalA',
                     inp=INC / 'CG/input.npz', controls=INC / 'CG/e2e/controls/controls_baseline.csv'),
    'EvalA_M2': dict(line='H1', eval=SS.EVALS['H1'], harness=INC / 'M2/e2e/salt_none/harness/evalA',
                     inp=INC / 'M2/input.npz', controls=INC / 'M2/e2e/controls/controls_baseline.csv'),
    'HCT116_CG': dict(line='HCT116', eval=SS.EVALS['HCT116'],
                      harness=V / 'results/stage1/e2e_final/CG/HCT116/salt_none/harness/evalHCT116final',
                      inp=V / 'results/stage1/e2e_final/CG/HCT116/model/predictions.npz',
                      controls=V / 'results/stage1/e2e_final/CG/HCT116/controls/controls_baseline.csv'),
}


def actual(ref, hdir):
    """The run's own per-knockdown values (results.csv) through the same aggregation + scaling, checked against
    summary.json; and cell-eval2's DE dispatch re-run on the run's own de_pred.parquet, checked against results.csv."""
    import polars as pl
    summ = json.loads((hdir / 'summary.json').read_text())
    w = pd.read_csv(hdir / 'results.csv').pivot(index='perturbation', columns='metric', values='value').reindex(ref.h.targets)
    kd = {MSE + ':num': w['expr_mse_unbiased_capped'].to_numpy(), MSE + ':den': w['expr_distance_unbiased'].to_numpy(),
          PDS: w[PDS].to_numpy(), **{m: w[m].to_numpy() for m in DE}}
    out = ref.summarise(kd, per_kd=True)
    chk = {m: abs(out['scaled'][m] - summ['scaled'][m]) for m in SIX}
    chk['overall'] = abs(out['overall'] - summ['scaled']['overall'])
    redo = ref.de_metrics(pl.read_parquet(hdir / 'de_pred.parquet'))
    plumb = {m: float(np.nanmax(np.abs(redo[m] - kd[m]))) for m in DE}
    return out, dict(scaled_abs_diff_vs_summary=chk, de_dispatch_max_abs_diff_vs_results_csv=plumb)


def call_diff(ref, hdir, sig, controls):
    """Realised (harness de_pred, p_adj < 0.05, own gene out) minus requested calls, per knockdown and by the target's
    context-control CPM decile (gated genes; decile -1 = ungated)."""
    de = pd.read_parquet(hdir / 'de_pred.parquet', columns=['target', 'feature', 'p_adj'])
    de = de[(de.p_adj < 0.05) & (de.feature != de.target)]
    h = ref.h; gp = {g: j for j, g in enumerate(h.genes)}; tp = {t: i for i, t in enumerate(h.targets)}
    real = np.zeros((h.T, h.G), bool); real[de.target.map(tp).to_numpy(), de.feature.map(gp).to_numpy()] = True
    req = sig.copy()
    for i in h.rows:
        if h.own[i] >= 0: req[i, h.own[i]] = False
    cb = pd.read_csv(controls, index_col=0)
    m0 = cb['control_mean_cpm'].reindex(h.genes).to_numpy(float); gate = cb['de_tested'].reindex(h.genes).fillna(False).to_numpy(bool)
    gi = np.flatnonzero(gate); edges = np.quantile(m0[gi], np.linspace(0.1, 0.9, 9))
    dec = np.where(gate, np.searchsorted(edges, m0, side='right'), -1)
    per = pd.DataFrame(dict(target=h.targets, requested=req.sum(1), realised=real.sum(1), both=(req & real).sum(1),
                            realised_not_requested=(real & ~req).sum(1), requested_not_realised=(req & ~real).sum(1)))
    per['realised_minus_requested'] = per.realised - per.requested
    by = []
    for d in range(-1, 10):
        g = dec == d
        by.append(dict(decile=d, genes=int(g.sum()), requested=int(req[:, g].sum()), realised=int(real[:, g].sum()),
                       realised_not_requested=int((real & ~req)[:, g].sum()), requested_not_realised=int((req & ~real)[:, g].sum()),
                       per_kd_realised_minus_requested_mean=float((real[:, g].sum(1) - req[:, g].sum(1)).mean())))
    return per, pd.DataFrame(by)


def step9(a):
    ratio = SA.ratio_from_step4()
    refs, res, rows, dec = {}, {}, [], []
    for run, r in RUNS9.items():
        key = str(r['eval'])
        if key not in refs:
            refs[key] = Ref(r['eval'], ratio, name=r['line'])
        ref = refs[key]
        kds, lfc, sig = requested(r['inp'])
        assert kds == ref.h.targets, run
        m0 = pd.read_csv(r['controls'], index_col=0)['control_mean_cpm'].reindex(ref.h.genes).to_numpy(float)
        prox = ref.proxy(lfc, sig, m0)
        act, checks = actual(ref, r['harness'])
        per, by = call_diff(ref, r['harness'], sig, r['controls'])
        per.insert(0, 'run', run); by.insert(0, 'run', run); dec.append(by)
        per.to_parquet(OUT / f'tables/step9_calls_{run}.parquet')
        res[run] = dict(proxy={k: prox[k] for k in ('raw', 'scaled', 'overall', 'mse_unclamped')},
                        actual={k: act[k] for k in ('raw', 'scaled', 'overall', 'mse_unclamped')}, checks=checks,
                        calls=dict(requested_median=float(per.requested.median()), realised_median=float(per.realised.median()),
                                   realised_minus_requested_median=float(per.realised_minus_requested.median()),
                                   realised_minus_requested_mean=float(per.realised_minus_requested.mean()),
                                   realised_not_requested_median=float(per.realised_not_requested.median()),
                                   requested_not_realised_median=float(per.requested_not_realised.median())))
        for m in SIX + ['overall']:
            pv = prox['overall'] if m == 'overall' else prox['scaled'][m]; av = act['overall'] if m == 'overall' else act['scaled'][m]
            rows.append(dict(run=run, component=SHORT.get(m, m), proxy_raw=prox['raw'].get(m), actual_raw=act['raw'].get(m),
                             proxy_scaled=pv, actual_scaled=av, scaled_gap=pv - av))
        print(run, json.dumps(res[run], default=float), flush=True)
    pd.DataFrame(rows).to_csv(OUT / 'tables/step9_proxy_vs_actual.csv', index=False)
    pd.concat(dec).to_csv(OUT / 'tables/step9_calls_by_decile.csv', index=False)
    (OUT / 'tables/step9.json').write_text(json.dumps(dict(ratio=ratio, runs=res), indent=1, default=float) + '\n')
    print(pd.DataFrame(rows).to_string())


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['step9'])
    a = p.parse_args()
    {'step9': step9}[a.cmd](a)

"""Step 1 analysis: per-knockdown decomposition of the scorer's MSE for arms A-C, arm D (no cells), the gene-group
and shared-error tables, the effect-size terciles, the Step 0c knockdown table, and the alignment / units / seed checks.

Scorer space = cell-eval2's own pseudobulk (streaming_bulk.inmem_pseudobulk, bulk_lognorm: log1p(5e4 x group sum /
group total)) with its delete-1 jackknife moments. The MSE pieces come from cell-eval2's metrics.delta helpers in the
order _numerator applies them (the value is checked against delta.mse_unbiased_capped and the harness's results.csv).

  ~/venv/bin/python code/vcc6/step1_analyze.py
"""
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
import rpe1_eval as E                                                                      # noqa: E402
import step1_oracle as O                                                                   # noqa: E402
sys.path.insert(0, str(V / 'code/cellgen'))
import cell_generator as C                                                                 # noqa: E402
from cell_eval2.metrics import delta as D                                                  # noqa: E402
from cell_eval2.moments import unbiased_sq_dist                                            # noqa: E402
from cell_eval2.streaming_bulk import inmem_pseudobulk                                     # noqa: E402

OUT = O.OUT
CTRL = E.H.CONTROL
TS = 5e4
NORM = 'bulk_lognorm'


def bulk_of(path):
    a = ad.read_h5ad(path)
    a.X = a.X.tocsr()
    b, m = inmem_pseudobulk(a, pert_col='target_gene', norms=[NORM], target_sum=None, with_moments=True,
                            per_gene_jk=True)
    perts, means = b[NORM]
    return (np.asarray(perts).astype(str), np.asarray(means, np.float64)), m[NORM], a.var_names.astype(str).tolist()


def mse_parts(pb, pm, rb, rm, genes):
    """_numerator(cap=True), step by step, with cell-eval2's helpers. Per knockdown, in the metric's units (/|G_p|)."""
    got = D._aligned_pair(pb, rb, pm, rm, None, CTRL, 'expr_mse_unbiased_capped', comparator=NORM)
    keep, rows, pred_means, real_means, pred_tn, real_tn, real_index, pred_perts = got
    G = pred_means.shape[1]
    labels = [str(pred_perts[i]) for i in keep]
    excl = D._exclusion_cols(labels, genes, target_gene_map=None,
                             gate_labels=[p for p in real_index if p != CTRL], who='expr_mse_unbiased_capped', n_genes=G)
    pred_sel = pred_means[keep]
    c_raw = pred_tn[keep]
    c_cap = np.minimum(c_raw, D.PRED_TRACE_CAP_K * real_tn[rows])
    row_w = D._row_weights(keep.size, G, excl)
    claim = float((c_cap * row_w).sum())
    budget = float(D._across_pert_budget(pred_sel, excl, row_w))
    rho = min(1.0, budget / claim) if claim > 0 else 1.0
    c_used = c_cap * rho

    def per_gene(pt, rt):
        s = unbiased_sq_dist(pred_sel, real_means[rows], pt, rt)
        s, n_row = D._drop_on_target(s, pred_means, keep, real_means, rows, G, excl)
        return D._per_gene(s, n_row), n_row
    value, n_row = per_gene(c_used, real_tn[rows])
    raw, _ = per_gene(np.zeros_like(c_used), np.zeros_like(c_used))
    t = pd.DataFrame(dict(knockdown=labels, raw_sq=raw, c_pred_jackknife=c_raw / n_row, c_pred_capped=c_cap / n_row,
                          c_pred_used=c_used / n_row, c_real=real_tn[rows] / n_row, value=value, genes_scored=n_row))
    return t.set_index('knockdown'), dict(rho=rho, budget=budget, claim=claim, cap_k=float(D.PRED_TRACE_CAP_K),
                                          cap_fired=int((c_cap < c_raw).sum()))


def denominators(rb, rm, genes):
    d = D.distance_unbiased(comparator=NORM, real_bulk=rb, real_moments=rm, control=CTRL, genes=genes)
    return pd.Series(d, name='denominator')


def main():
    tg = O.targets(); genes = E.H.panel_genes(); G = len(genes)
    d = E.H.eval_dir(E.NAME)
    rb, rm, rgenes = bulk_of(d / 'real.h5ad')
    assert rgenes == genes
    rperts, rmeans = rb; ci = int(np.flatnonzero(rperts == CTRL)[0]); rrow = {p: i for i, p in enumerate(rperts)}
    ctrl_pb = rmeans[ci]
    den = denominators(rb, rm, genes)
    lfc, sig, se, tobs = O.table_rows(tg, genes)
    measured = np.isfinite(lfc)
    on_panel = measured.any(0)                                  # genes the RPE1 table measures
    m0 = O.control_mean(genes)
    gene_col = {g: j for j, g in enumerate(genes)}
    own = np.array([gene_col.get(t, -1) for t in tg])

    # ---- Step 0c table: cells, de_sig count, tercile
    obs, _ = E.cells_obs()
    _, label, _ = E.rpe1_kept(obs)
    counts = pd.Series(label[label != '']).value_counts()
    split = json.loads((d / 'split.json').read_text())
    kd = pd.DataFrame(dict(knockdown=tg, real_cells=[int(counts[t]) for t in tg],
                           reference_cells=[int(split['cells_per_knockdown'][t]) for t in tg],
                           table_cells_used=tobs['n_cells_used'].to_numpy(int), table_de_sig=sig.sum(1)))
    kd['denominator'] = den.reindex(tg).to_numpy()
    q = kd.denominator.rank(method='first', ascending=False)
    kd['tercile'] = np.where(q <= 100, 'strong', np.where(q <= 200, 'mid', 'weak'))
    kd = kd.set_index('knockdown')
    kd.to_csv(OUT / 'knockdowns_300.csv')
    terc = kd['tercile']

    # ---- arms A-C: harness values + decomposition
    res, parts_all = {}, {}
    for arm in ('A', 'B', 'C'):
        s = json.loads((E.H.RES / f'step1_arm{arm}' / f'eval{E.NAME}' / 'summary.json').read_text())
        rc = pd.read_csv(E.H.RES / f'step1_arm{arm}' / f'eval{E.NAME}' / 'results.csv')
        rc = rc.pivot_table(index='perturbation', columns='metric', values='value')
        pb, pm, pgenes = bulk_of(OUT / f'arm{arm}_pred.h5ad')
        assert pgenes == genes
        parts, glob = mse_parts(pb, pm, rb, rm, genes)
        check = D.mse_unbiased_capped(comparator=NORM, pred_bulk=pb, real_bulk=rb, pred_moments=pm, real_moments=rm,
                                      control=CTRL, genes=genes)
        parts['denominator'] = den.reindex(parts.index)
        parts['harness_capped'] = rc['expr_mse_unbiased_capped'].reindex(parts.index)
        parts['harness_denominator'] = rc['expr_distance_unbiased'].reindex(parts.index)
        parts['ratio'] = parts.value / parts.denominator
        parts['tercile'] = terc.reindex(parts.index)
        # scorer-space deltas
        pperts, pmeans = pb; prow = {p: i for i, p in enumerate(pperts)}
        rr, slope_t, slope_d = [], [], []
        for t in parts.index:
            g = on_panel.copy()
            if gene_col.get(t, -1) >= 0: g[gene_col[t]] = False
            gen = pmeans[prow[t]] - ctrl_pb; real = rmeans[rrow[t]] - ctrl_pb
            rr.append(float(np.corrcoef(gen[g], real[g])[0, 1]) if gen[g].std() > 0 and real[g].std() > 0 else np.nan)
        parts['pearson_gen_real'] = rr
        parts.to_csv(OUT / f'arm{arm}_per_knockdown.csv')
        parts_all[arm] = (parts, pb)
        res[arm] = dict(scaled=s['scaled'], raw=s['raw'], numerator_global=glob,
                        recomputed_vs_harness=dict(
                            max_abs_diff_capped=float(np.nanmax(np.abs(parts.value - parts.harness_capped))),
                            max_abs_diff_delta_fn=float(max(abs(check[k] - parts.value[k]) for k in parts.index)),
                            max_abs_diff_denominator=float(np.nanmax(np.abs(parts.denominator - parts.harness_denominator)))),
                        norm_ratio_of_sums=float(parts.value.sum() / parts.denominator.sum()),
                        sums=dict(raw_sq=float(parts.raw_sq.sum()), c_pred_used=float(parts.c_pred_used.sum()),
                                  c_real=float(parts.c_real.sum()), value=float(parts.value.sum()),
                                  denominator=float(parts.denominator.sum())),
                        median_pearson_gen_real=float(np.nanmedian(parts.pearson_gen_real)))

    # ---- input delta (the oracle's intended scorer-space change) and arm D
    def ideal_pb(lfc_row):
        m1 = C.target_means(lfc_row, m0)                                    # the generator's own inversion
        m1 = np.where(np.isfinite(m0), m1, 0.0)
        return np.log1p(TS * m1 / m1.sum())
    ideal_ctrl = ideal_pb(np.zeros(G))
    ideal = np.vstack([ideal_pb(np.where(measured[i], lfc[i], np.nan)) for i in range(len(tg))])
    table_pb = None
    t_all = ad.read_h5ad(V / 'lib/effects_v5/replogle_rpe1.h5ad', backed='r')
    trow = {s: i for i, s in enumerate(t_all.obs['target_gene'].astype(str))}
    table_pb = np.vstack([np.nan_to_num(np.asarray(t_all.layers['log1p_bulk_50k_matched'][trow[t]], float), nan=0.0) for t in tg])
    base = pd.read_csv(O.BASE, index_col=0)
    table_ctrl_pb = np.nan_to_num(base['log1p_bulk_50k_matched'].to_numpy(float), nan=0.0)

    def no_cell_arm(P, tag):
        """A prediction with no cells: pseudobulk rows P (300 x G), no sampling correction (C_pred = 0)."""
        perts = np.array(list(tg) + [CTRL]); means = np.vstack([P, ideal_ctrl])
        order = np.argsort(perts); pbx = (perts[order], means[order])
        zero = type(rm)(perts=pbx[0], counts=np.full(len(perts), 2.0), sumsq=np.zeros(len(perts)),
                        jk=np.zeros(len(perts)), jk_gene=np.zeros((len(perts), G)))
        parts, glob = mse_parts(pbx, zero, rb, rm, genes)
        parts['denominator'] = den.reindex(parts.index); parts['tercile'] = terc.reindex(parts.index)
        parts.to_csv(OUT / f'arm{tag}_per_knockdown.csv')
        return dict(raw_sq_sum=float(parts.raw_sq.sum()), value_sum=float(parts.value.sum()),
                    norm=float(parts.value.sum() / parts.denominator.sum()),
                    norm_clamped=float(min(1.0, max(0.0, parts.value.sum() / parts.denominator.sum()))),
                    by_tercile={k: dict(norm=float(g.value.sum() / g.denominator.sum()), raw_sq=float(g.raw_sq.sum()))
                                for k, g in parts.groupby('tercile')}), parts
    res['D'], partsD = no_cell_arm(ideal, 'D')
    res['D2_table_bulk_matched'], partsD2 = no_cell_arm(table_pb, 'D2')

    # ---- arm C: transmission and definition slopes, gene-group splits, shared error
    pc, pbc = parts_all['C']
    pperts, pmeans = pbc; prow = {p: i for i, p in enumerate(pperts)}
    rows = []
    for i, t in enumerate(tg):
        g = on_panel.copy()
        if own[i] >= 0: g[own[i]] = False
        inp = ideal[i] - ideal_ctrl; gen = pmeans[prow[t]] - ctrl_pb; real = rmeans[rrow[t]] - ctrl_pb
        v = np.var(inp[g])
        rows.append(dict(knockdown=t, slope_gen_on_input=float(np.cov(gen[g], inp[g])[0, 1] / v) if v > 0 else np.nan,
                         slope_real_on_input=float(np.cov(real[g], inp[g])[0, 1] / v) if v > 0 else np.nan,
                         pearson_input_real=float(np.corrcoef(inp[g], real[g])[0, 1]) if v > 0 else np.nan))
    sl = pd.DataFrame(rows).set_index('knockdown')
    pc = pc.join(sl)
    pc.to_csv(OUT / 'armC_per_knockdown.csv')

    m0dec = pd.Series(np.nan, index=range(G))
    meas_cols = np.flatnonzero(on_panel & np.isfinite(m0) & (m0 > 0))
    m0dec[meas_cols] = pd.qcut(m0[meas_cols], 10, labels=False) + 1
    def splits(arm_pb, arm_parts, tag):
        pp, pmn = arm_pb; prw = {p: i for i, p in enumerate(pp)}
        R = np.vstack([pmn[prw[t]] - rmeans[rrow[t]] for t in tg])          # residual pred - real, scorer space
        for i in range(len(tg)):
            if own[i] >= 0: R[i, own[i]] = np.nan
        E2 = R ** 2; tot = np.nansum(E2)
        groups = {'inside RPE1 table panel': on_panel, 'outside RPE1 table panel': ~on_panel}
        out = dict(total_raw_sq=float(tot), panel={k: float(np.nansum(E2[:, m]) / tot) for k, m in groups.items()})
        out['expression_decile'] = {int(k): float(np.nansum(E2[:, (m0dec == k).to_numpy()]) / tot) for k in range(1, 11)}
        called = sig & on_panel[None, :]
        out['called_vs_uncalled'] = dict(called=float(np.nansum(E2[called]) / tot), uncalled=float(np.nansum(E2[~called]) / tot))
        rbar = np.nanmean(R, axis=0)
        shared = float(np.nansum(np.where(np.isfinite(R), rbar[None, :] ** 2, 0.0)) / tot)
        out['shared_error_share'] = shared
        top = np.argsort(-np.nan_to_num(rbar ** 2))[:15]
        out['shared_top_genes'] = [dict(gene=genes[j], mean_residual=float(rbar[j]), control_mean_cpm=float(m0[j]) if np.isfinite(m0[j]) else None,
                                        in_table=bool(on_panel[j])) for j in top]
        out['by_tercile'] = {}
        for k in ('strong', 'mid', 'weak'):
            idx = np.array([terc[t] == k for t in tg])
            Ek = E2[idx]; tk = np.nansum(Ek); Rk = R[idx]; rb_k = np.nanmean(Rk, 0)
            out['by_tercile'][k] = dict(
                total_raw_sq=float(tk),
                panel={n: float(np.nansum(Ek[:, m]) / tk) for n, m in groups.items()},
                expression_decile={int(j): float(np.nansum(Ek[:, (m0dec == j).to_numpy()]) / tk) for j in range(1, 11)},
                called=float(np.nansum(Ek[called[idx]]) / tk),
                shared_error_share=float(np.nansum(np.where(np.isfinite(Rk), rb_k[None, :] ** 2, 0.0)) / tk))
        return out
    res['C_splits'] = splits(pbc, pc, 'C')
    res['B_splits'] = splits(parts_all['B'][1], parts_all['B'][0], 'B')
    def tercile_table(p):
        return {k: dict(norm=float(g.value.sum() / g.denominator.sum()), raw_sq=float(g.raw_sq.sum()),
                        c_pred_used=float(g.c_pred_used.sum()), c_real=float(g.c_real.sum()),
                        denominator=float(g.denominator.sum()), median_pearson=float(np.nanmedian(g.pearson_gen_real))
                        if 'pearson_gen_real' in g else None, knockdowns=len(g))
                for k, g in p.groupby('tercile')}
    res['by_tercile'] = {a: tercile_table(parts_all[a][0]) for a in ('A', 'B', 'C')}
    res['C_slopes'] = dict(median_slope_gen_on_input=float(np.nanmedian(pc.slope_gen_on_input)),
                           median_slope_real_on_input=float(np.nanmedian(pc.slope_real_on_input)),
                           median_pearson_input_real=float(np.nanmedian(pc.pearson_input_real)),
                           median_pearson_gen_real=float(np.nanmedian(pc.pearson_gen_real)),
                           by_tercile={k: dict(slope_gen_on_input=float(np.nanmedian(g.slope_gen_on_input)),
                                               slope_real_on_input=float(np.nanmedian(g.slope_real_on_input)),
                                               pearson_gen_real=float(np.nanmedian(g.pearson_gen_real)))
                                       for k, g in pc.groupby('tercile')})

    # ---- checks
    qz = np.load(O.QUANT)
    base = pd.read_csv(O.BASE, index_col=0)
    real_var = ad.read_h5ad(d / 'real.h5ad', backed='r').var_names.astype(str).tolist()
    predC_var = ad.read_h5ad(OUT / 'armC_pred.h5ad', backed='r').var_names.astype(str).tolist()
    t_var = t_all.var_names.astype(str).tolist()
    qcols = qz['panel_col'].astype(int)
    gen_rep = pd.read_csv(OUT / 'armC_pred_generator.csv')
    res['checks'] = dict(
        gene_order_mismatches=dict(table_lfc_vs_panel=int(sum(a != b for a, b in zip(t_var, genes))),
                                   baseline_csv_vs_panel=int(sum(a != b for a, b in zip(base.index.astype(str), genes))),
                                   scorer_reference_vs_panel=int(sum(a != b for a, b in zip(real_var, genes))),
                                   prediction_vs_panel=int(sum(a != b for a, b in zip(predC_var, genes))),
                                   lengths=dict(table=len(t_var), baseline=len(base), reference=len(real_var), panel=G),
                                   quantile_cols_outside_table_measured=int((~on_panel[qcols]).sum()),
                                   quantile_cols_not_gated_in_baseline=int((~base['de_tested'].to_numpy(bool)[qcols]).sum()),
                                   baseline_gated_not_in_quantiles=int(base['de_tested'].sum() - len(qcols))),
        seeds=dict(distinct=int(gen_rep.seed.nunique()), knockdowns=len(gen_rep)),
        generator=dict(median_sig_requested=float(gen_rep.sig_requested.median()), median_sig_realized=float(gen_rep.sig_realized.median()),
                       median_library=float(gen_rep.library_median.median()), median_composition=float(gen_rep.composition.median())),
        depth=dict(reference_median_panel_library=split.get('median_panel_library_after_perturbed'),
                   reference_controls_median_panel_library=split.get('median_panel_library_after_controls'),
                   generator_library_size=20000))
    (OUT / 'step1.json').write_text(json.dumps(res, indent=1, default=str) + '\n')
    print(json.dumps({k: res[k] for k in res if k not in ('C_splits', 'B_splits')}, indent=1, default=str))


if __name__ == '__main__':
    main()

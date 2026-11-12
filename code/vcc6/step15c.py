"""Step 15c (Jose 2026-10-30): config 1' on cells, one config per VM, through the frozen generator and harness.
Outputs in results/stage1/2026-10-29_step15/c/<config>/<line>/; write-up in the scorer-space report.md, Step 15.

  config1p_K0       (vcc2026)  config 1' magnitudes, C+G's calls
  config1p_K4_1200  (vcc2027)  config 1' magnitudes, per knockdown the top 1,200 tested genes by C+G z = |LFC| / SE
Config 1' = vcc_2026_7's config 1 with the non-called scale m 0.6 (CONFIG1P). The magnitude split is always C+G's calls (K0), as
in the proxy (step15.py): K4 changes only the call set handed to the generator.

  model   --config C              vcc_2026_7's data-loader + model cells (code copied, vcc_2026_7 untouched) with CONFIG1P, for EvalA
                                  and HCT116final (own line's tables out; HEK293T votes for HCT116, the screen setup)
                                  -> <line>/input.npz (lfc, de_sig = this config's calls, n_sources; extras: de_sig_K0,
                                  de_sig_K4, cg_lfc, se). Checks on EvalA: m 0.5 reproduces vcc_2026_7's freeze-check input; calls,
                                  LFC and SE equal the scorer-space fold cache (diag then checks the proxy on input.npz
                                  against 15a: 0.2470 K0, 0.2711 K4 1,200).
  (generate / score: run_line.py predict --lfc-npz / run_line.py score, unchanged)
  boot    --config C              boot_score.py of this config against the incumbent C+G run, per line
  diag    --config C              items 2-4 for this config: calls requested vs realised by control-CPM decile
                                  (scorer_proxy.call_diff), requested-not-realised, leaked; item 3 (median |requested
                                  Delta| and generated / real cell-level SD ratio for the calls only in K4, with K0 calls
                                  and uncalled genes as reference rows); proxy (scorer_proxy.Ref.proxy on input.npz) minus
                                  actual (the harness run)
  compare                         on vcc2026 with both configs' harness folders: paired CIs (vs C+G, K4 vs K0), all six
                                  components scaled and raw, by true-strength tercile on Eval A, the adoption rule
"""
import argparse
import gc
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
for p in (V / 'code/vcc6', V / 'scripts', V / 'scripts/harness', V / 'code/cellgen'):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

OUT = V / 'results/stage1/2026-10-29_step15/c'
CONFIG1P = dict(s_c=0.5, m=0.6, kappa=1.0, c1=1.0)
CONFIG1 = dict(s_c=0.5, m=0.5, kappa=1.0, c1=1.0)          # frozen config 1, EvalA check only
N_K4 = 1200
CONFIGS = ('config1p_K0', 'config1p_K4_1200')
FE = V / 'results/stage1/final_evals/harness/lib'
LINES = {
    'EvalA': dict(controls=V / 'lib/harness/evalA/context.h5ad', knockdowns=V / 'results/target_lists/evalA_h1_292.txt',
                  background='H1', reference=V / 'lib/harness/evalA/real.h5ad', proxy_line='H1',
                  incumbent=V / 'results/stage1/2026-10-27_incumbent/evalA/CG/e2e/salt_none/harness/evalA',
                  incumbent_input=V / 'results/stage1/2026-10-27_incumbent/evalA/CG/input.npz'),
    'HCT116final': dict(controls=FE / 'evalHCT116final/context.h5ad', knockdowns=V / 'results/target_lists/final_HCT116.txt',
                        background='HCT116', reference=FE / 'evalHCT116final/real.h5ad', proxy_line='HCT116',
                        incumbent=V / 'results/stage1/e2e_final/CG/HCT116/salt_none/harness/evalHCT116final',
                        incumbent_input=V / 'results/stage1/e2e_final/CG/HCT116/model/predictions.npz'),
}
COMPS = ['mse', 'pds', 'nmae', 'fidelity', 'reach', 'jaccard']
TERC = ('weak', 'mid', 'strong')


def log(msg):
    print(f'[{time.strftime("%H:%M:%S", time.gmtime())}] {msg}', flush=True)


def kds(line):
    return [k for k in Path(LINES[line]['knockdowns']).read_text().split() if k]


def harness_dir(cfg, line):
    return OUT / cfg / line / 'e2e/salt_none/harness' / LINES[line]['reference'].parent.name


def config_lfc(cfg, cg_lfc, cg_sig, dpool, n_bg, z, m0, cprof):
    """vcc_2026_7's config1() with the configuration passed in."""
    z2 = np.nan_to_num(z.astype(float)) ** 2
    c = np.where((n_bg >= 2) & np.isfinite(z), z2 / (z2 + cfg['kappa']), np.where(n_bg == 1, cfg['c1'], 0.0))
    mu = np.maximum(20.0 * np.expm1(cprof + cfg['m'] * c * dpool), 0.0)
    return np.where(cg_sig, cfg['s_c'] * cg_lfc, np.log2((mu + 1e-9) / (m0 + 1e-9))).astype(np.float32)


def k4_calls(cg_lfc, se, n_src, gate, own, n=N_K4):
    """Top n tested genes by C+G z (scorer_step16.call_sets / step15.cgz_order + top_n)."""
    import stage1 as S
    lfc = np.where(n_src <= 0, 0.0, cg_lfc.astype(np.float64))
    ok = gate & (n_src > 0) & (lfc != 0)
    if own >= 0:
        ok[own] = False
    zc = np.where(ok, np.abs(lfc) / np.maximum(se, S.SE_MIN), -1.0)
    order = np.argsort(-zc, kind='stable')
    sig = np.zeros(len(lfc), bool); sig[order[:int(min(n, ok.sum()))]] = True
    return sig


# ---------------------------------------------------------------- model (vcc_2026_7 cells 1-2, config passed in)
def model(a, lines=None, out=None, extra_checks=None):
    """lines: {name: dict(controls, knockdowns, background)}; background None = no table held out (Step 18's D/E/F).
    extra_checks(line, targets, B, A, n_src) -> dict, asserted like the EvalA checks."""
    import h5py
    from anndata.io import read_elem
    import stage1 as S
    import run_line as RL
    cfg_name = a.config; assert cfg_name in CONFIGS
    lines = LINES if lines is None else lines; out = OUT if out is None else out
    kds_ = lambda l: [k for k in Path(lines[l]['knockdowns']).read_text().split() if k]
    GENES = RL.panel_genes(); G_ = np.asarray(GENES, dtype=object); gpos = {g: j for j, g in enumerate(GENES)}
    todo = [l for l in lines if not (out / cfg_name / l / 'input.npz').exists()]
    log(f'{cfg_name}: model for {todo}; config {CONFIG1P}')
    if not todo:
        return
    TABLES = S.discover_tables([str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')])
    BY_STEM = {p.stem: p for p in TABLES}
    SOURCES = S.load_tables(TABLES)
    need = set().union(*(kds_(l) for l in todo))
    DSRC = []
    for src in SOURCES:
        p = BY_STEM[src.name]
        ctrl = pd.read_csv(p.with_name(src.name + '_baseline.csv'), index_col=0)['log1p_bulk_50k_full'].reindex(GENES).to_numpy(float)
        with h5py.File(p, 'r') as f:
            kd = read_elem(f['obs'])['target_gene'].astype(str).to_numpy()
            rows = [i for i, k in enumerate(kd) if k in need]
            if not rows or not np.isfinite(ctrl).any():
                continue
            full, se = np.asarray(f['layers/log1p_bulk_50k_full'][rows], float), np.asarray(f['layers/se'][rows], np.float32)
        if not np.isfinite(full).any():
            continue
        D, k_, n = (full - ctrl).astype(np.float32), kd[rows], np.full(len(rows), np.nan, np.float32)
        if src.modality != 'drug':
            k_, D, se, n = S._collapse_duplicate_knockdowns(k_, D, se, n)
        DSRC.append(S.Source(src.name, src.background, src.modality, src.reference, src.genes, k_, D, se, n, src.mult))
    log(f'{len(SOURCES)} tables; {len(DSRC)} with Delta rows')
    for l in todo:
        d = out / cfg_name / l; d.mkdir(parents=True, exist_ok=True)
        X = RL.read_controls(lines[l]['controls'], GENES)
        m0, gate, _, _ = RL.control_profile(X, d / 'controls')
        raw = np.asarray(X.sum(0)).ravel().astype(float); cprof = np.log1p(5e4 * raw / raw.sum()); del X
        targets = kds_(l)
        bg = lines[l]['background']
        others = list(SOURCES) if bg is None else [s for s in SOURCES if S._line_key(s.background) != S._line_key(bg)]
        basals = {}
        for s in others:
            bp = BY_STEM[s.name].with_name(s.name + '_baseline.csv')
            got = S.load_baseline(bp, s.genes) if bp.exists() else None
            if got is not None: basals[s.name] = (s.genes, got[0])
        votes, prior_votes = S.arm_votes(others, 'C')
        calib, prior = S.fit_se_calibration(others, basals), S.fit_prior(others, votes=prior_votes)
        dvotes = S.BackgroundVotes([s for s in DSRC if s.name in {o.name for o in others}])
        K, G = len(targets), len(GENES)
        A = {k: np.zeros((K, G), np.float32) for k in ('lfc', 'lfc_m05', 'cg_lfc', 'se')}
        B = {k: np.zeros((K, G), bool) for k in ('K0', 'K4')}
        n_src = np.zeros((K, G), np.int16)
        for i, k in enumerate(targets):
            cg = S.predict_knockdown(others, k, G_, m0, gate, calib, prior, drop_gene=k, votes=votes, magnitude='gaussian')
            dp = S.predict_knockdown(None, k, G_, m0, gate, calib, prior, drop_gene=k, votes=dvotes, magnitude='gaussian')
            dpool = np.where(dp['n_sources'] > 0, dp['lfc'], 0.0)
            Y, _, _ = dvotes.gather(k, G_)
            n_bg, z = np.zeros(G, int), np.full(G, np.nan)
            if Y:
                Y = np.stack(Y); fin = np.isfinite(Y); n_bg = fin.sum(0)
                mean = np.where(fin, Y, 0.0).sum(0) / np.maximum(n_bg, 1)
                sd = np.sqrt(np.where(fin, (Y - mean) ** 2, 0.0).sum(0) / np.maximum(n_bg - 1, 1))
                z = np.where(n_bg >= 2, np.abs(mean) / np.maximum(sd / np.sqrt(np.maximum(n_bg, 1)), 1e-12), np.nan)
            A['lfc'][i] = config_lfc(CONFIG1P, cg['lfc'], cg['de_sig'], dpool, n_bg, z, m0, cprof)
            A['lfc_m05'][i] = config_lfc(CONFIG1, cg['lfc'], cg['de_sig'], dpool, n_bg, z, m0, cprof)
            A['cg_lfc'][i], A['se'][i] = cg['lfc'], cg['se_measurement']
            B['K0'][i], n_src[i] = cg['de_sig'], cg['n_sources']
            B['K4'][i] = k4_calls(cg['lfc'], cg['se_measurement'], cg['n_sources'], gate, gpos.get(k, -1))
        sig = B['K0'] if cfg_name.endswith('K0') else B['K4']
        checks = {}
        if l == 'EvalA':
            fz = V / 'results/7/freeze_check/EvalA/input.npz'
            if fz.exists():
                z7 = np.load(fz, allow_pickle=True); assert list(map(str, z7['knockdowns'])) == targets
                checks['v7_m05_lfc_max_abs_diff'] = float(np.abs(z7['lfc'] - A['lfc_m05']).max())
                checks['v7_calls_identical'] = bool((z7['de_sig'] == B['K0']).all())
            fc = V / 'results/stage1/2026-10-28_scorer_space/cache/fold_H1.npz'
            if fc.exists():
                f = np.load(fc, allow_pickle=True); r = [list(map(str, f['targets'])).index(t) for t in targets]
                checks['fold_cache_calls_identical'] = bool((f['de_sig'][r] == B['K0']).all())
                checks['fold_cache_lfc_max_abs_diff'] = float(np.abs(f['lfc'][r] - A['cg_lfc']).max())
                checks['fold_cache_se_max_abs_diff'] = float(np.abs(f['se_measurement'][r] - A['se']).max())
        if extra_checks is not None:
            checks.update(extra_checks(l, targets, B, A, n_src))
        if checks:
            log(f'{l} checks {checks}')
            assert all(v is True or (not isinstance(v, bool) and v < 1e-5) for v in checks.values()), checks
        np.savez_compressed(d / 'input.npz', knockdowns=np.array(targets, dtype=object), genes=G_, lfc=A['lfc'],
                            de_sig=sig, n_sources=n_src, de_sig_K0=B['K0'], de_sig_K4=B['K4'], cg_lfc=A['cg_lfc'],
                            se=A['se'], m0=m0, gate=gate, cprof=cprof)
        (d / 'model.json').write_text(json.dumps(dict(config=cfg_name, params=CONFIG1P, tables_voting=len(others),
                                                      calls_K0_median=float(np.median(B['K0'].sum(1))),
                                                      calls_K4_median=float(np.median(B['K4'].sum(1))),
                                                      calls_this_median=float(np.median(sig.sum(1))), checks=checks), indent=1) + '\n')
        log(f'{l}: {K} knockdowns, {len(others)} tables vote; calls K0 median {np.median(B["K0"].sum(1)):.0f}, '
            f'K4 {np.median(B["K4"].sum(1)):.0f}')
    del SOURCES, DSRC; gc.collect()


# ---------------------------------------------------------------- boot vs C+G (per VM)
def boot(a):
    import boot_score as BS
    for l, s in LINES.items():
        out = OUT / a.config / l / 'boot_vs_CG.json'
        if out.exists():
            continue
        res = BS.boot([BS.Run(s['incumbent']), BS.Run(harness_dir(a.config, l))])
        out.write_text(json.dumps(res, indent=1) + '\n')
        c = next(iter(res['compare'].values()))['overall']
        log(f'{a.config} {l}: overall vs C+G {c["diff"]:+.4f} [{c["ci"][0]:+.4f}, {c["ci"][1]:+.4f}]')


# ---------------------------------------------------------------- items 2-4 (per VM)
def cell_sd(path, targets, pairs, label_col='target_gene'):
    """Cell-level SD of log1p(CP10k) for (knockdown, gene) pairs: pairs = {knockdown: gene indices}."""
    import anndata as ad
    a = ad.read_h5ad(path, backed='r')
    lab = a.obs[label_col].astype(str).to_numpy()
    out = {}
    for k in targets:
        g = pairs.get(k)
        if g is None or not len(g):
            continue
        rows = np.flatnonzero(lab == k)
        X = a.X[rows[0]:rows[-1] + 1] if (np.diff(rows) == 1).all() else a[rows].to_memory().X
        X = X.tocsr() if hasattr(X, 'tocsr') else X
        tot = np.asarray(X.sum(1)).ravel().astype(float)
        Y = np.asarray(X[:, g].todense() if hasattr(X, 'todense') else X[:, g], float)
        Y = np.log1p(1e4 * Y / np.maximum(tot, 1)[:, None])
        out[k] = Y.std(0, ddof=1)
    a.file.close()
    return out


def diag(a):
    import scorer_proxy as SP
    import scorer_arms as SA
    cfg = a.config
    ratio = SA.ratio_from_step4()
    res = {}
    for l, s in LINES.items():
        d = OUT / cfg / l; hd = harness_dir(cfg, l)
        ref = SP.Ref(s['reference'].parent, ratio, name=s['proxy_line'])
        z = np.load(d / 'input.npz', allow_pickle=True)
        names = list(map(str, z['knockdowns'])); r = [names.index(t) for t in ref.h.targets]
        un = z['n_sources'][r] <= 0
        lfc = np.where(un, 0.0, z['lfc'][r].astype(np.float64)); sig = z['de_sig'][r] & ~un
        controls = d / 'e2e/controls/controls_baseline.csv'
        m0 = pd.read_csv(controls, index_col=0)['control_mean_cpm'].reindex(ref.h.genes).to_numpy(float)
        prox = ref.proxy(lfc, sig, m0)
        act, checks = SP.actual(ref, hd)
        per, by = SP.call_diff(ref, hd, sig, controls)
        per.to_parquet(d / 'calls_per_kd.parquet'); by.to_csv(d / 'calls_by_decile.csv', index=False)
        gap = {m: dict(proxy=prox['scaled'][m], actual=act['scaled'][m], gap=prox['scaled'][m] - act['scaled'][m],
                       proxy_raw=prox['raw'][m], actual_raw=act['raw'][m]) for m in SP.SIX}
        gap['overall'] = dict(proxy=prox['overall'], actual=act['overall'], gap=prox['overall'] - act['overall'])
        # item 3: requested Delta (scorer bulk space) and the cell-level SD ratio, by gene set
        k0 = z['de_sig_K0'][r] & ~un; k4 = z['de_sig_K4'][r] & ~un
        own = np.zeros_like(k0)
        for i in ref.h.rows:
            if ref.h.own[i] >= 0: own[i, ref.h.own[i]] = True
        gate = pd.read_csv(controls, index_col=0)['de_tested'].reindex(ref.h.genes).fillna(False).to_numpy(bool)
        sets = {'K4_only': k4 & ~k0 & ~own, 'K0_called': k0 & ~own, 'uncalled_gated': gate[None, :] & ~k0 & ~k4 & ~own & ~un}
        rng = np.random.default_rng(2026)
        P = SA.bulk(SA.mu_from_lfc(lfc, m0)); dreq = P - ref.h.ctrl[None, :]
        T = SP_truth(ref)
        item3 = {}
        for nm, M in sets.items():
            pairs = {}
            for i, t in enumerate(ref.h.targets):
                g = np.flatnonzero(M[i])
                if len(g) > 60:
                    g = np.sort(rng.choice(g, 60, replace=False))       # 60 genes per knockdown for the SD ratio
                pairs[t] = g
            gen = cell_sd(d / 'e2e/salt_none/prediction.h5ad', ref.h.targets, pairs)
            real = cell_sd(s['reference'], ref.h.targets, pairs)
            ratio_sd, dq = [], []
            for i, t in enumerate(ref.h.targets):
                if t in gen and t in real:
                    ok = real[t] > 0
                    ratio_sd.append(gen[t][ok] / real[t][ok])
            ratio_sd = np.concatenate(ratio_sd) if ratio_sd else np.array([np.nan])
            item3[nm] = dict(pairs=int(M.sum()), per_kd_median=float(np.median(M.sum(1))),
                             abs_requested_delta_median=float(np.median(np.abs(dreq[M]))) if M.any() else None,
                             abs_real_delta_median=float(np.median(np.abs((ref.h.R - ref.h.ctrl[None, :])[M]))) if M.any() else None,
                             true_de_share=float(T[M].mean()) if M.any() else None,
                             sd_ratio_pairs=int(np.isfinite(ratio_sd).sum()), sd_ratio_median=float(np.nanmedian(ratio_sd)),
                             sd_ratio_q25=float(np.nanquantile(ratio_sd, 0.25)), sd_ratio_q75=float(np.nanquantile(ratio_sd, 0.75)))
            log(f'{cfg} {l} item3 {nm}: {item3[nm]}')
        res[l] = dict(proxy_minus_actual=gap, checks=checks, item3=item3,
                      calls=dict(requested_median=float(per.requested.median()), realised_median=float(per.realised.median()),
                                 requested_not_realised_total=int(per.requested_not_realised.sum()),
                                 realised_not_requested_total=int(per.realised_not_requested.sum()),
                                 requested_total=int(per.requested.sum()), realised_total=int(per.realised.sum()),
                                 requested_not_realised_median=float(per.requested_not_realised.median()),
                                 realised_not_requested_median=float(per.realised_not_requested.median())))
        log(f'{cfg} {l}: proxy {prox["overall"]:.4f} actual {act["overall"]:.4f}; calls {res[l]["calls"]}')
        del ref; gc.collect()
    (OUT / cfg / 'diag.json').write_text(json.dumps(res, indent=1, default=float) + '\n')


def SP_truth(ref):
    import scorer_step12 as S12
    return S12.truth_calls(ref)


# ---------------------------------------------------------------- compare (vcc2026, both configs)
def compare(a):
    import boot_score as BS
    import scorer_proxy as SP
    import scorer_arms as SA
    tab = OUT / 'tables'; tab.mkdir(parents=True, exist_ok=True)
    rows, terc_rows, res = [], [], {}
    for l, s in LINES.items():
        runs = {'CG': BS.Run(s['incumbent']), **{c: BS.Run(harness_dir(c, l)) for c in CONFIGS}}
        common = runs['CG'].per_kd.index
        for r in runs.values():
            common = common.intersection(r.per_kd.index)
        parts = {'all': common}
        if l == 'EvalA':
            ref = SP.Ref(s['reference'].parent, SA.ratio_from_step4(), name='H1')
            n_true = pd.Series(SP_truth(ref).sum(1), index=ref.h.targets).reindex(common)
            rk = n_true.rank(method='first').to_numpy(); n = len(rk)
            tt = np.where(rk <= n / 3, 0, np.where(rk <= 2 * n / 3, 1, 2))
            for j, t in enumerate(TERC):
                parts[t] = common[tt == j]
            res['EvalA_terciles'] = {t: dict(knockdowns=int((tt == j).sum()), true_de_min=int(n_true[tt == j].min()),
                                             true_de_median=float(n_true[tt == j].median()), true_de_max=int(n_true[tt == j].max()))
                                     for j, t in enumerate(TERC)}
        pairs = [('config1p_K0', 'CG'), ('config1p_K4_1200', 'CG'), ('config1p_K4_1200', 'config1p_K0')]
        for part, idx in parts.items():
            rng = np.random.default_rng(BS.SEED)
            draws = {k: [] for k in runs}
            for _ in range(BS.N_BOOT):
                ii = idx[rng.integers(0, len(idx), len(idx))]
                for k, r in runs.items():
                    draws[k].append(r.scaled(ii))
            F = {k: pd.DataFrame(v) for k, v in draws.items()}
            sh = lambda d_: {SP.SHORT.get(m, m): v for m, v in d_.items()}
            pt = {k: sh(r.scaled(idx)) for k, r in runs.items()}
            raw = {k: sh(r.aggregate(idx).to_dicts()[0]) for k, r in runs.items()}
            F = {k: f.rename(columns=SP.SHORT) for k, f in F.items()}
            for k in runs:
                for m in F[k].columns:
                    rows.append(dict(line=l, part=part, knockdowns=len(idx), run=k, metric=m, scaled=pt[k][m],
                                     ci_lo=F[k][m].quantile(0.025), ci_hi=F[k][m].quantile(0.975), raw=raw[k].get(m)))
            for c, b in pairs:
                dd = F[c] - F[b]
                for m in dd.columns:
                    terc_rows.append(dict(line=l, part=part, knockdowns=len(idx), cand=c, ref=b, metric=m,
                                          diff=pt[c][m] - pt[b][m], ci_lo=dd[m].quantile(0.025), ci_hi=dd[m].quantile(0.975),
                                          p_better=(dd[m] > 0).mean(),
                                          raw_diff=(raw[c][m] - raw[b][m]) if m in raw[c] else None))
            log(f'{l} {part}: done ({len(idx)} knockdowns)')
        for k, r in runs.items():
            r.check()
    P = pd.DataFrame(rows); D = pd.DataFrame(terc_rows)
    P.to_csv(tab / 'step15c_scores.csv', index=False); D.to_csv(tab / 'step15c_paired.csv', index=False)
    # adoption rule, per candidate
    adopt = {}
    for c in CONFIGS:
        e = P[(P.line == 'EvalA') & (P.part == 'all') & (P.run == c)].set_index('metric').scaled
        d = D[(D.line == 'EvalA') & (D.part == 'all') & (D.cand == c) & (D.ref == 'CG')].set_index('metric')
        h = D[(D.line == 'HCT116final') & (D.part == 'all') & (D.cand == c) & (D.ref == 'CG')].set_index('metric')
        rule = dict(mse_pos=bool(e['mse'] > 0), others_pos=bool((e[[m for m in COMPS if m != 'mse']] > 0).all()),
                    overall_above_CG=bool(d.loc['overall', 'diff'] > 0), overall_ci_lo_gt_m001=bool(d.loc['overall', 'ci_lo'] > -0.01),
                    hct_within_ci=bool(h.loc['overall', 'ci_lo'] <= 0 <= h.loc['overall', 'ci_hi'] or h.loc['overall', 'ci_lo'] > 0))
        rule['pass'] = all(rule.values())
        adopt[c] = rule
    res['adoption'] = adopt
    (tab / 'step15c.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(json.dumps(res, indent=1, default=float))
    with pd.option_context('display.width', 250, 'display.max_rows', 500):
        print(D[D.metric.isin(COMPS + ['overall'])].round(4).to_string())



if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['model', 'boot', 'diag', 'compare'])
    p.add_argument('--config', choices=CONFIGS)
    a = p.parse_args()
    {'model': model, 'boot': boot, 'diag': diag, 'compare': compare}[a.cmd](a)

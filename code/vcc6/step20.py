"""Step 20 (Jose 2026-11-02): lab effects, the strong-knockdown shift and program regression v2, before the next D/E/F
file. CPU only; no cells generated, no harness scoring of new cells, no D/E/F file, no upload. Write-up: section 12 of
results/stage1/2026-10-28_scorer_space/report.md. Outputs: results/stage1/2026-11-02_step20/{tables,logs}. Frozen code
(generator, harness, run_line, boot_score, config 1' files, v7) and the Step 18 / 19 outputs are imported or read only.

vcc2026 (step20.sh):
  part0   checks: config 1' + K0 + no-source fix on the generator input reproduces Step 19 Part F's row; the section 11
          follow-up's H1-own tercile shift (cells) reproduces its held-out MSE rows
  partA   lab against biology in the source tables: lab map, pair classes, cosines (LFC raw / confidence-weighted /
          Delta), n_real Spearman, cells; per-class summaries and a regression with table-bootstrap CIs
  partB   same-lab upweighting on HCT116 / HEK293T / K562 pools (one table load), w in {1, 2, 5, 20}, proxy scored
  partC   source-derived strong-minus-weak shift S by (b) tercile; 40-row grid; fit on one H1 half, score the other
  partD   S_H1 and S_H1,clean (pluripotency PC1 projected out): cosines, program shares, guardrails at Part C's picks
  partG   the recommended spec's combination row on the proxy
vcc2027 (step20_prog.sh): partF1 .. partF6 (program regression v2; see the partF* docstrings)
  test    unit tests on small slices (no table load, no proxy)
"""
import argparse
import gc
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
for p in (V / 'code/vcc6', V / 'scripts', V / 'scripts/harness', V / 'code/cellgen', V / 'analysis/mse_metric'):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import step19 as S19                                                                       # noqa: E402
from step19 import (log, rss_gb, avail_gb, mem_before, mem_after, save_json, grid_pos, Fast, COMPS, TERC,  # noqa: E402
                    tercile_of, ref_for, b_for, DLine, C15, D_LINES, OWN_BG)

OUT = V / 'results/stage1/2026-11-02_step20'
TAB, LOG = OUT / 'tables', OUT / 'logs'
S19T = S19.TAB
SS_OUT = V / 'results/stage1/2026-10-28_scorer_space'
SEED = 2026
B_LINE = {'H1': 'EvalA', 'K562': 'K562final', 'HCT116': 'HCT116final', 'HEK293T': 'HEK293Tfinal'}
NOSRC = {'nosource_fix': True}
DEF_TARGETS = lambda: pd.read_csv(V / 'data_final/extracted/pert_counts.csv')['target_gene'].astype(str).tolist()


# ================================================================ shared
def base_mu(D):
    """Config 1' + K0 with the no-source fix (Step 19 Part F's combo_mu); the fold pools already carry it."""
    return S19.combo_mu(D, NOSRC)


def prof(mu):
    """Generator means (CPM scale) -> log1p CP50K profile (the inverse of scorer_arms.mu_from_profile)."""
    return np.log1p(np.maximum(mu, 0.0) / 20.0)


def mu_add_delta(mu, add, own):
    """config 1''s means with a Delta-space term added (log1p CP50K), converted back as config 1''s Delta-space values
    are (mu = 20 (exp(profile) - 1), floored at 0). NaN entries of add and each knockdown's own gene are left alone."""
    import scorer_arms as SA
    A = np.where(np.isfinite(add), add, 0.0)
    for i, g in enumerate(own):
        if g >= 0:
            A[i, g] = 0.0
    return SA.mu_from_profile(prof(mu) + A)


def score_rows(D, mu, sig=None, **tag):
    """Proxy all six components for one prediction, summarised on each of D.parts."""
    out = D.ref.proxy(mu=mu, sig=D.sig if sig is None else sig, per_kd=True)
    rows = []
    for part, rr in D.parts.items():
        sm = D.ref.summarise(out['per_kd'], rows=rr)
        r = dict(line=D.line, part=part, overall=sm['overall'], **tag)
        for m, k in Fast.KEYS.items():
            r[m] = sm['scaled'][k]
        rows.append(r)
    return rows, out


def b_cuts(b):
    """(b) tercile cut points from a set of knockdowns' finite b: the 1/3 and 2/3 quantiles."""
    b = np.asarray(b, float); b = b[np.isfinite(b)]
    return np.quantile(b, [1 / 3, 2 / 3])


def b_terc(b, cuts):
    """Absolute (b) tercile: 0 below the first cut, 2 at or above the second; no b -> weak (0)."""
    b = np.asarray(b, float)
    t = np.where(b >= cuts[1], 2, np.where(b >= cuts[0], 1, 0))
    return np.where(np.isfinite(b), t, 0)


def own_index(h):
    return np.asarray(h.own)


def pick_fit(R, fit_part, held_part, key, base_sel):
    """Highest overall on fit_part among rows with all six components > 0 there; that row on held_part; the base."""
    fit = R[R.part == fit_part]
    elig = fit[(fit[COMPS] > 0).all(1)]
    if not len(elig):
        return None
    pk = elig.sort_values('overall', ascending=False).iloc[0]
    sel = np.ones(len(R), bool)
    for k in key:
        sel &= (R[k] == pk[k]).to_numpy()
    ho = R[sel & (R.part == held_part).to_numpy()].iloc[0]
    b0 = R[base_sel(R) & (R.part == held_part)].iloc[0]
    return pk, ho, b0, len(elig)


def comp_dict(r):
    return {c: float(r[c]) for c in COMPS + ['overall']}


# ================================================================ Part 0: checks
def part0(a):
    import mse_study as M
    import scorer_step12 as S12
    import scorer_step16 as S16
    res = {}
    # (1) config 1' + K0 + no-source fix on the generator's input = Step 19 Part F's 'nosource_fix' row (H1 and folds)
    F = pd.read_csv(S19T / 'partF_rows.csv')
    res['check1'] = {}
    for line in D_LINES:
        D = DLine(line)
        t0 = time.time(); rows, _ = score_rows(D, base_mu(D)); dt = time.time() - t0
        want = F[(F.line == line) & (F.row == 'nosource_fix')].set_index('part')
        diffs = {r['part']: abs(r['overall'] - want.loc[r['part'], 'overall']) for r in rows}
        res['check1'][line] = dict(got={r['part']: r['overall'] for r in rows}, want=want['overall'].to_dict(), max_diff=max(diffs.values()),
                                   seconds_per_row=dt, threads=os.environ.get('OMP_NUM_THREADS'))
        log(f'part0 check1 {line}: {res["check1"][line]["got"]} vs Step 19 {res["check1"][line]["want"]}; max |diff| {max(diffs.values()):.1e}; '
            f'one proxy row {dt:.1f} s at OMP_NUM_THREADS={os.environ.get("OMP_NUM_THREADS")}')
        del D; gc.collect()
    ok1 = all(v['max_diff'] < 1e-9 for v in res['check1'].values())
    # (2) the follow-up's H1-own tercile-mean shift, learned on one half (true terciles), assigned by (b) tercile
    ref = ref_for('EvalA'); h = ref.h
    T = S12.truth_calls(ref); terc = tercile_of(T.sum(1)); halves = S16.halves(T)
    rows = own_shift_cells(h, terc, halves)
    C = pd.read_csv(S19T / 'mse_terciles_corrections.csv')
    res['check2'] = []
    for r in rows:
        w = C[(C.arm == r['arm']) & (C.fit_half == r['fit_half']) & (C.part == r['part'])].iloc[0]
        res['check2'].append(dict(**r, want=float(w.scaled), diff=abs(r['scaled'] - float(w.scaled))))
    for r in res['check2']:
        if r['part'] == 'all':
            log(f'part0 check2 fit half {r["fit_half"]} {r["arm"]}: held-out MSE {r["scaled"]:.4f} vs follow-up {r["want"]:.4f}')
    ok2 = all(r['diff'] < 1e-9 for r in res['check2'])
    res['pass'] = dict(check1=bool(ok1), check2=bool(ok2))
    save_json(TAB / 'part0.json', res)
    log(f'part0: check1 {"PASS" if ok1 else "FAIL"}, check2 {"PASS" if ok2 else "FAIL"}')
    if not (ok1 and ok2):
        sys.exit(3)


def own_shift_cells(h, terc, halves):
    """Section 11 follow-up (mseterc2), arms 'uploaded' and the (b)-assigned tercile-mean shift, on the uploaded cells."""
    import mse_study as M
    b = b_for('EvalA', h.targets); bt = np.full(h.T, 0)
    ok = np.isfinite(b); bt[ok] = tercile_of(b[ok])
    z = np.load(S19T / 'predside_K0_EvalA.npz'); P, C = z['P'], z['C']
    Rr = P - h.R
    for i in h.rows:
        if h.own[i] >= 0: Rr[i, h.own[i]] = 0.0
    rows = []
    for f in (0, 1):
        fit, ho = halves[f], halves[1 - f]
        for name in ('uploaded (no correction)', 'tercile-mean shift (mid, strong), assigned by (b) tercile'):
            Rh = np.zeros_like(Rr)
            if name != 'uploaded (no correction)':
                for j in (1, 2):
                    Rh[np.intersect1d(ho, np.flatnonzero(bt == j))] = Rr[np.intersect1d(fit, np.flatnonzero(terc == j))].mean(0)
            _, num = M.score_generated(h, P - Rh, C)
            for t in ('all',) + TERC:
                rr = ho if t == 'all' else np.intersect1d(ho, np.flatnonzero(terc == TERC.index(t)))
                u = num[rr].sum() / h.den[rr].sum()
                rows.append(dict(arm=name, fit_half=f, part=t, scaled=h.scale(Fast.KEYS['mse'], u)))
    return rows


# ================================================================ Part A: lab against biology
# Lab of each genetic table (table-name prefix), from the publication / provenance; flagged where unsure.
LAB_RULES = [('replogle_', 'Weissman'), ('adamson_', 'Weissman'), ('nadig_', 'Weissman'), ('orion_', 'Xaira (X-Atlas/Orion)'),
             ('arc_', 'Arc'), ('feng_', 'Feng'), ('kolf_', 'KOLF2.1J study'), ('jiang_', 'Jiang (Satija)'), ('marson_', 'Marson')]
LAB_FLAGS = {
    'Weissman': 'Replogle 2022 (K562 GW, RPE1) and Adamson 2016 are Weissman-lab; Nadig 2025 (HepG2, Jurkat) is Weissman-lab with '
                'the Replogle/Nadig authors: sure.',
    'Xaira (X-Atlas/Orion)': 'Orion = X-Atlas/Orion (Xaira Therapeutics, 2025), both lines one platform; the pseudobulk tables are '
                             'the same screens\' complement rows: fairly sure.',
    'Arc': 'Arc Institute (the VCC training screen): sure.',
    'Feng': 'Feng et al. CRISPRi Perturb-seq across HipSci iPSC donors (targeted and genome-wide); our metadata records the study, '
            'not the lab: UNSURE of the lab name, sure the tables share one lab and protocol.',
    'KOLF2.1J study': 'KOLF_Strong_Perturbations.h5ad (KRAB-ZIM3-dCas9 KOLF2.1J CRISPRi); the producing lab is not recorded in '
                      'our download notes: UNSURE (treated as its own lab).',
    'Jiang (Satija)': 'Jiang et al. (Parse split-pool CRISPRko, 6 lines x 5 cytokines); I believe Satija lab (NYGC): UNSURE.',
    'Marson': 'Marson lab CD4 T (replicate tables, no cells): sure.',
}
FENG_DONOR = lambda bg: bg.split('_')[0] if '_' in bg and bg != 'HipSci_pooled' else bg
LINEAGE = {'H1': 'pluripotent', 'HipSci_pooled': 'pluripotent', 'KOLF2.1J': 'pluripotent',
           'K562': 'blood', 'Jurkat': 'blood', 'CD4_T': 'blood', 'HAP1': 'HAP1 (CML-derived, adherent)',
           'HCT116': 'colorectal', 'HT29': 'colorectal', 'HEK293T': 'kidney', 'HepG2': 'liver', 'RPE1': 'retina',
           'A549': 'lung', 'BXPC3': 'pancreas', 'MCF7': 'breast'}
LINE_FLAGS = {'kolf': 'KOLF2.1J is a gene-corrected subclone of HPSI0114i-kolf_2; Feng kolf_2 / kolf_3 are clones of the same '
                      'donor. All three are line "kolf" here.',
              'feng clones': 'Feng clones of one donor (eipl_1 / eipl_3, ...) are the same line (one donor genome).',
              'HipSci_pooled': 'Feng genome-wide tables pool many HipSci donors: own line, lineage pluripotent.',
              'HAP1': 'HAP1 (near-haploid, from KBM-7 CML) is its own lineage, not counted as blood.'}
PAIR_CLASSES = ['same lab, same line', 'same lab, same lineage, different line', 'same lab, different line',
                'different lab, same line', 'different lab, same lineage', 'different lab, different lineage']


def lab_of(stem):
    return next(l for p, l in LAB_RULES if stem.startswith(p))


def line_of(bg):
    """Line identity for 'same line': Feng clones -> donor; KOLF2.1J -> kolf."""
    if bg == 'KOLF2.1J':
        return 'kolf'
    if bg in LINEAGE:
        return bg
    return FENG_DONOR(bg)


def lineage_of(bg):
    return LINEAGE.get(bg, 'pluripotent')                       # Feng donor clones are HipSci iPSC lines


def pair_class(lab1, line1, lin1, lab2, line2, lin2):
    if lab1 == lab2:
        if line1 == line2:
            return PAIR_CLASSES[0]
        return PAIR_CLASSES[1] if lin1 == lin2 else PAIR_CLASSES[2]
    if line1 == line2:
        return PAIR_CLASSES[3]
    return PAIR_CLASSES[4] if lin1 == lin2 else PAIR_CLASSES[5]


def genetic_tables():
    """Every genetic (non-drug) source table: (path, stem, obs)."""
    import h5py
    import stage1 as S
    from anndata.io import read_elem
    out = []
    for p in S.discover_tables([str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')]):
        if p.stem.startswith('tahoe_'):
            continue
        with h5py.File(p, 'r') as f:
            o = read_elem(f['obs'])
        if (o['modality'].astype(str) == 'drug').all():
            continue
        out.append((p, p.stem, o))
    return out


def collapse(kd, X):
    """Average duplicate knockdown rows (NaN-aware); returns (unique kd, rows)."""
    u, inv = np.unique(kd, return_inverse=True)
    if len(u) == len(kd):
        o = np.argsort(kd); return kd[o], X[o]
    S_ = np.zeros((len(u), X.shape[1]), np.float64); C_ = np.zeros_like(S_)
    np.add.at(S_, inv, np.nan_to_num(X)); np.add.at(C_, inv, np.isfinite(X))
    return u, np.where(C_ > 0, S_ / np.maximum(C_, 1), np.nan).astype(np.float32)


def read_rows(p, o, want, layers, dtype=np.float32):
    """Rows of a table for the wanted knockdowns (duplicates averaged), per layer, stored as dtype, chunked h5 reads.
    layer 'delta' = log1p_bulk_50k_full - the table's control profile; 'wlfc' = LFC z^2 / (z^2 + 1)."""
    import h5py
    import run_line as RL
    genes = RL.panel_genes()
    kd_all = o['target_gene'].astype(str).to_numpy()
    idx = np.flatnonzero(np.isin(kd_all, list(want)))
    b = pd.read_csv(p.with_name(p.stem + '_baseline.csv'), index_col=0).reindex(genes)
    out = {}
    with h5py.File(p, 'r') as f:
        for lay in layers:
            src = {'lfc': 'lfc', 'se': 'se', 'delta': 'log1p_bulk_50k_full'}[lay]
            if lay == 'delta' and not np.isfinite(f['layers/log1p_bulk_50k_full'][:min(3, len(kd_all))]).any():
                continue
            blocks = []
            for c0 in range(0, len(idx), 2000):
                kk = idx[c0:c0 + 2000]
                if not len(kk):
                    continue
                x = f[f'layers/{src}'][kk[0]:kk[-1] + 1][kk - kk[0]].astype(np.float32)
                if lay == 'delta':
                    x = x - b['log1p_bulk_50k_full'].to_numpy(np.float32)[None, :]
                blocks.append(x)
            out[lay] = np.concatenate(blocks) if blocks else np.zeros((0, len(genes)), np.float32)
    if 'lfc' in out and 'se' in out:
        z2 = (out['lfc'] / out['se']) ** 2
        out['wlfc'] = (out['lfc'] * z2 / (z2 + 1.0)).astype(np.float32)
    kd = kd_all[idx]
    res = {}
    for lay, X in out.items():
        if lay == 'se':
            continue
        u, Y = collapse(kd, X); res[lay] = Y.astype(dtype)
    nr = o['n_real'].astype(float).to_numpy()[idx]; u, nrc = collapse(kd, nr[:, None])
    res['kd'] = u; res['n_real'] = nrc[:, 0]
    res['gate'] = (b['control_mean_cpm'].to_numpy(float) > 5.0)
    return res


def pair_cos(A, B, M):
    """Per-row cosine of A and B over mask M (finite in both and gated in both); NaN where a row has no energy."""
    A = np.where(M, A, 0.0); B = np.where(M, B, 0.0)
    num = (A * B).sum(1, dtype=np.float64)
    den = np.sqrt((A * A).sum(1, dtype=np.float64) * (B * B).sum(1, dtype=np.float64))
    return np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)


def partA(a):
    from scipy.stats import spearmanr
    T = genetic_tables()
    meta = []
    for p, stem, o in T:
        bg = str(o['background'].iloc[0])
        nc = o['n_cells_used'].astype(float).to_numpy() if 'n_cells_used' in o else np.full(len(o), np.nan)
        meta.append(dict(table=stem, lab=lab_of(stem), background=bg, line=line_of(bg), lineage=lineage_of(bg),
                         modality=str(o['modality'].iloc[0]), rows=len(o), knockdowns=int(o['target_gene'].astype(str).nunique()),
                         median_cells=float(np.nanmedian(nc)) if np.isfinite(nc).any() else np.nan,
                         median_n_real=float(np.nanmedian(o['n_real'].astype(float))) if np.isfinite(o['n_real'].astype(float)).any() else np.nan))
    meta = pd.DataFrame(meta); meta.to_csv(TAB / 'partA_labs.csv', index=False)
    kds = {stem: set(o['target_gene'].astype(str)) for _, stem, o in T}
    stems = [s for _, s, _ in T]
    pairs = [(i, j) for i in range(len(T)) for j in range(i + 1, len(T)) if len(kds[stems[i]] & kds[stems[j]]) >= 20]
    want = {s: set() for s in stems}
    for i, j in pairs:
        sh = kds[stems[i]] & kds[stems[j]]; want[stems[i]] |= sh; want[stems[j]] |= sh
    nrow = sum(len(v) for v in want.values())
    log(f'A: {len(T)} genetic tables, {len(pairs)} pairs sharing >= 20 knockdowns, {nrow} table rows to read')
    mem_before(f'A rows ({nrow} x 18080 x 3 layers float16, plus the largest table\'s float64 read buffers)', nrow * 18080 * 3 * 2 / 1e9 + 4)
    data = {}
    for p, stem, o in T:
        if want[stem]:
            data[stem] = read_rows(p, o, want[stem], ('lfc', 'se', 'delta'), dtype=np.float16)
    mem_after('A rows')
    rows = []
    for i, j in pairs:
        s1, s2 = stems[i], stems[j]; d1, d2 = data[s1], data[s2]
        sh = np.array(sorted(kds[s1] & kds[s2]))
        i1 = np.searchsorted(d1['kd'], sh); i2 = np.searchsorted(d2['kd'], sh)
        g = d1['gate'] & d2['gate']
        m1, m2 = meta.set_index('table').loc[s1], meta.set_index('table').loc[s2]
        r = dict(table1=s1, table2=s2, lab1=m1.lab, lab2=m2.lab, line1=m1.line, line2=m2.line, lineage1=m1.lineage, lineage2=m2.lineage,
                 cls=pair_class(m1.lab, m1.line, m1.lineage, m2.lab, m2.line, m2.lineage), shared=len(sh),
                 cells1=m1.median_cells, cells2=m2.median_cells, genes=int(g.sum()))
        nr1, nr2 = d1['n_real'][i1], d2['n_real'][i2]
        strong = (nr1 >= 20) & (nr2 >= 20)
        r['shared_de20'] = int(strong.sum())
        ok = np.isfinite(nr1) & np.isfinite(nr2)
        r['spearman_n_real'] = float(spearmanr(nr1[ok], nr2[ok]).statistic) if ok.sum() >= 5 and np.ptp(nr1[ok]) > 0 and np.ptp(nr2[ok]) > 0 else np.nan
        for lay, nm in (('lfc', 'raw'), ('wlfc', 'cw'), ('delta', 'delta')):
            if lay not in d1 or lay not in d2:
                r[f'cos_{nm}_all'] = r[f'cos_{nm}_de20'] = np.nan; continue
            A, B = d1[lay][i1].astype(np.float32), d2[lay][i2].astype(np.float32)
            M = np.isfinite(A) & np.isfinite(B) & g[None, :]
            c = pair_cos(A, B, M)
            r[f'cos_{nm}_all'] = float(np.nanmean(c)) if np.isfinite(c).any() else np.nan
            r[f'cos_{nm}_de20'] = float(np.nanmean(c[strong])) if np.isfinite(c[strong]).sum() >= 3 else np.nan
        rows.append(r)
    P = pd.DataFrame(rows); P.to_csv(TAB / 'partA_pairs.csv', index=False)
    del data; gc.collect()
    res = partA_summary(P, meta)
    save_json(TAB / 'partA.json', res)
    print(json.dumps({k: v for k, v in res.items() if k != 'labs'}, indent=1, default=float))


A_MEASURES = ['cos_raw_all', 'cos_cw_all', 'cos_raw_de20', 'cos_cw_de20', 'cos_delta_all', 'cos_delta_de20', 'spearman_n_real']


def iqr_summary(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if not len(x):
        return dict(n=0)
    q = np.quantile(x, [0.25, 0.5, 0.75])
    return dict(n=int(len(x)), median=float(q[1]), q25=float(q[0]), q75=float(q[2]))


def partA_summary(P, meta, n_boot=1000):
    res = dict(labs=meta.to_dict('records'), lab_flags=LAB_FLAGS, line_flags=LINE_FLAGS, classes={}, contrasts={}, regression={})
    for c in PAIR_CLASSES:
        d = P[P.cls == c]
        res['classes'][c] = dict(pairs=int(len(d)), knockdowns=int(d.shared.sum()), knockdowns_de20=int(d.shared_de20.sum()),
                                 **{m: iqr_summary(d[m]) for m in A_MEASURES})
    sl_dl = P[P.cls == 'same lab, different line']; dl_sl = P[P.cls == 'different lab, same line']
    res['contrasts']['same lab / different line vs different lab / same line'] = {
        m: dict(same_lab_diff_line=iqr_summary(sl_dl[m]), diff_lab_same_line=iqr_summary(dl_sl[m])) for m in A_MEASURES}
    h1 = P[(P.table1 == 'arc_h1') | (P.table2 == 'arc_h1')]
    other_lin = np.where(h1.table1 == 'arc_h1', h1.lineage2, h1.lineage1)
    res['contrasts']['H1 vs pluripotent tables of other labs, vs H1 vs non-pluripotent'] = {
        m: dict(pluripotent=iqr_summary(h1[other_lin == 'pluripotent'][m]), non_pluripotent=iqr_summary(h1[other_lin != 'pluripotent'][m]))
        for m in A_MEASURES}
    res['h1_pairs'] = h1[['table1', 'table2', 'cls', 'shared', 'shared_de20'] + A_MEASURES].to_dict('records')
    for y in ('cos_cw_all', 'cos_cw_de20', 'cos_delta_all'):
        res['regression'][y] = pair_regression(P, y, n_boot)
    return res


def reg_design(P):
    c1, c2 = np.log(P.cells1.to_numpy(float)), np.log(P.cells2.to_numpy(float))
    X = np.column_stack([np.ones(len(P)), (P.lab1 == P.lab2).to_numpy(float), (P.line1 == P.line2).to_numpy(float),
                         (P.lineage1 == P.lineage2).to_numpy(float), np.minimum(c1, c2), np.maximum(c1, c2)])
    return X, ['intercept', 'same_lab', 'same_line', 'same_lineage', 'log_cells_smaller', 'log_cells_larger']


def pair_regression(P, y, n_boot=1000):
    """OLS of pair cosine on {same lab, same line, same lineage, log cells of both tables (smaller, larger)}; 95% CIs from a
    bootstrap over tables (draw tables with replacement; a pair enters once per drawn copy of each of its two tables).
    Pairs without a cell count (Marson) or without y are left out."""
    d = P[np.isfinite(P[y]) & np.isfinite(P.cells1) & np.isfinite(P.cells2)].reset_index(drop=True)
    X, names = reg_design(d); Y = d[y].to_numpy(float)
    beta = np.linalg.lstsq(X, Y, rcond=None)[0]
    tabs = np.array(sorted(set(d.table1) | set(d.table2))); ti = {t: k for k, t in enumerate(tabs)}
    a1 = d.table1.map(ti).to_numpy(); a2 = d.table2.map(ti).to_numpy()
    rng = np.random.default_rng(SEED); B = []
    for _ in range(n_boot):
        cnt = np.bincount(rng.integers(0, len(tabs), len(tabs)), minlength=len(tabs))
        w = cnt[a1] * cnt[a2]
        if (w > 0).sum() < X.shape[1] + 2:
            continue
        sw = np.sqrt(w)
        try:
            B.append(np.linalg.lstsq(X * sw[:, None], Y * sw, rcond=None)[0])
        except np.linalg.LinAlgError:
            continue
    B = np.array(B)
    return dict(n_pairs=int(len(d)), n_tables=int(len(tabs)), n_boot=int(len(B)),
                coef={n: dict(beta=float(beta[k]), ci95=[float(np.nanquantile(B[:, k], 0.025)), float(np.nanquantile(B[:, k], 0.975))])
                      for k, n in enumerate(names)})


# ================================================================ Part B: same-lab upweighting (one table load)
B_POOLS = {'K562': dict(drop={'k562'}, up={'rpe1', 'hepg2', 'jurkat'}),
           'HCT116': dict(drop={'hct116'}, up={'hek293t'}),
           'HEK293T': dict(drop={'hek293t'}, up={'hct116'})}
B_W = [1, 2, 5, 20]


class MultiWeightedVotes:
    """scorer_study.WeightedVotes for a set of backgrounds: each one's vote variance / w (SE / sqrt w)."""

    def __init__(self, votes, backgrounds, w):
        self.v, self.bgs, self.w = votes, set(backgrounds), float(w)
        self.groups, self.owners = votes.groups, votes.owners

    def gather(self, knockdown, genes):
        Y, SE, mults = self.v.gather(knockdown, genes)
        bgs = list(self.v.owners.get(str(knockdown), ()))
        SE = [se / np.sqrt(self.w) if bg in self.bgs else se for bg, se in zip(bgs, SE)]
        return Y, SE, mults


def agreement(dvotes, k, genes):
    """Step 11's per-(k, g) agreement over the unweighted Delta background votes: n_bg and z = |mean| / (SD / sqrt n)."""
    Y, _, _ = dvotes.gather(k, genes)
    G = len(genes)
    if not Y:
        return np.zeros(G, np.int16), np.full(G, np.nan, np.float32)
    Y = np.stack(Y); fin = np.isfinite(Y); nn = fin.sum(0)
    Yz = np.where(fin, Y, 0.0); mean = Yz.sum(0) / np.maximum(nn, 1)
    sd = np.sqrt(np.where(fin, (Y - mean) ** 2, 0.0).sum(0) / np.maximum(nn - 1, 1))
    z = np.where(nn >= 2, np.abs(mean) / np.maximum(sd / np.sqrt(np.maximum(nn, 1)), 1e-12), np.nan)
    return nn.astype(np.int16), z.astype(np.float32)


def partB(a):
    """One table load; per pool: votes, calibration, prior, Delta votes; per w the C+G LFC pool and the Delta pool with
    the upweighted backgrounds' vote variance / w; the unweighted agreement z. Writes tables/partB_pool_<line>_w<w>.npz."""
    import stage1 as S
    import screen as SC
    import run_line as RL
    import scorer_study as SS
    genes = RL.panel_genes(); G_ = np.asarray(genes, dtype=object)
    tables = S.discover_tables(SC.DIRS); by_stem = {p.stem: p for p in tables}
    mem_before(f'B table load ({len(tables)} tables)', 20.0)
    t0 = time.time(); sources = S.load_tables(tables)
    mem_after(f'B table load ({time.time() - t0:.0f} s)')
    info = {}
    for line, pdef in B_POOLS.items():
        others = [s for s in sources if S._line_key(s.background) not in pdef['drop']]
        basals = {}
        for src in others:
            got = S.load_baseline(by_stem[src.name].with_name(src.name + '_baseline.csv'), src.genes)
            if got is not None: basals[src.name] = (src.genes, got[0])
        votes, prior_votes = S.arm_votes(others, 'C')
        calibration = S.fit_se_calibration(others, basals); prior = S.fit_prior(others, votes=prior_votes)
        dvotes = S.BackgroundVotes(SS.delta_sources(others, genes))
        fc = np.load(SS_OUT / 'cache' / f'fold_{line}.npz', allow_pickle=True)          # targets, m0, gate, ctrl profile
        m0, gate, targets = fc['m0'], fc['gate'], list(fc['targets'])
        held = np.array([any(bg in pdef['up'] for bg in prior_votes.owners.get(k, ())) for k in targets])     # genetic rows
        held_any = np.array([any(bg in pdef['up'] for bg in votes.owners.get(k, ())) for k in targets])     # incl. drug rows
        dheld = np.array([any(bg in pdef['up'] for bg in dvotes.owners.get(k, ())) for k in targets])
        info[line] = dict(tables_voting=len(others), prior=prior.as_dict(), upweighted=sorted(pdef['up']), dropped=sorted(pdef['drop']),
                          upweighted_present=sorted(set(pdef['up']) & set(votes.groups)), held_by_upweighted=int(held.sum()),
                          held_by_upweighted_incl_drug_rows=int(held_any.sum()),
                          held_by_upweighted_delta=int(dheld.sum()), knockdowns=len(targets),
                          upweighted_tables={bg: [s.name for s in votes.groups.get(bg, [])] for bg in pdef['up']})
        log(f'B pool {line}: {len(others)} tables vote; upweighted {info[line]["upweighted_tables"]}; '
            f'{held.sum()} / {len(targets)} knockdowns held by them; prior {prior.as_dict()}')
        K, G = len(targets), len(genes)
        nbg = np.zeros((K, G), np.int16); za = np.full((K, G), np.nan, np.float32)
        for i, k in enumerate(targets):
            nbg[i], za[i] = agreement(dvotes, k, G_)
        for w in B_W:
            V1, V2 = MultiWeightedVotes(votes, pdef['up'], w), MultiWeightedVotes(dvotes, pdef['up'], w)
            lfc = np.zeros((K, G), np.float32); dp = np.zeros((K, G), np.float32); sig = np.zeros((K, G), bool)
            se = np.zeros((K, G), np.float32); ns = np.zeros((K, G), np.int16); dn = np.zeros((K, G), np.int16)
            t1 = time.time()
            for i, k in enumerate(targets):
                pr = S.predict_knockdown(None, k, G_, m0, gate, calibration, prior, drop_gene=k, votes=V1, magnitude='gaussian')
                lfc[i], sig[i], ns[i], se[i] = pr['lfc'], pr['de_sig'], pr['n_sources'], pr['se_measurement']
                pd_ = S.predict_knockdown(None, k, G_, m0, gate, calibration, prior, drop_gene=k, votes=V2, magnitude='gaussian')
                dp[i], dn[i] = pd_['lfc'], pd_['n_sources']
                if i == 9:
                    log(f'B {line} w {w}: {(time.time() - t1) / 10:.2f} s per knockdown (both pools)')
            np.savez_compressed(TAB / f'partB_pool_{line}_w{w}.npz', targets=np.array(targets), genes=np.array(genes), lfc=lfc, de_sig=sig,
                                se_measurement=se, n_sources=ns, dpool=dp, dpool_n=dn, m0=m0, gate=gate, ctrl_profile=fc['ctrl_profile'],
                                prior=json.dumps(prior.as_dict()), n_bg=nbg, z_agree=za, held=held, w=w)
            log(f'B {line} w {w}: done in {time.time() - t1:.0f} s; calls median {np.median(sig.sum(1)):.0f}')
            if line == 'K562' and w == 1:                     # the K562 pool at w 1 is the K562 fold pool exactly
                s11 = np.load(SS_OUT / 'cache' / 's11_K562.npz')
                chk = dict(lfc=float(np.abs(lfc - fc['lfc']).max()), de_sig=int((sig != fc['de_sig']).sum()),
                           dpool=float(np.abs(np.where(dn > 0, dp, 0) - np.where(fc['dpool_n'] > 0, fc['dpool'], 0)).max()),
                           n_bg=int((nbg != s11['n_bg']).sum()),
                           z_agree_rel=float(np.nanmax(np.abs(np.nan_to_num(s11['z_agree'].astype(np.float64)) - np.nan_to_num(za.astype(np.float64)))
                                                      / np.maximum(1.0, np.abs(np.nan_to_num(s11['z_agree'].astype(np.float64)))))),
                           z_agree_nan_mismatch=int((np.isfinite(s11['z_agree']) != np.isfinite(za)).sum()))
                info['K562_w1_vs_fold_cache'] = chk
                log(f'B check: K562 w 1 vs fold_K562 / s11_K562: {chk}')
                if chk['lfc'] > 1e-5 or chk['de_sig'] or chk['dpool'] > 1e-5 or chk['n_bg'] or chk['z_agree_rel'] > 1e-4 or chk['z_agree_nan_mismatch']:
                    save_json(TAB / 'partB_pools.json', info); log('B CHECK FAILED: K562 w 1 does not reproduce the fold pool'); sys.exit(4)
        save_json(TAB / 'partB_pools.json', info)
        del votes, prior_votes, dvotes, calibration; gc.collect()
    # Arc H1's coverage of D/E/F: targets held, and the share of (target, gated gene) pairs where Arc H1 has a finite LFC
    h1 = next(s for s in sources if s.name == 'arc_h1')
    tg = DEF_TARGETS(); held = [t for t in tg if t in h1.rows_of]
    gi = {g: j for j, g in enumerate(h1.genes)}; jj = np.array([gi.get(g, -1) for g in genes])
    cov = {}
    for x in ('D', 'E', 'F'):
        gate = S19.control_stats(x)['gate']; n = 0
        for t in held:
            v = np.full(len(genes), np.nan); v[jj >= 0] = h1.lfc[h1.rows_of[t][0], jj[jj >= 0]]
            n += int((np.isfinite(v) & gate).sum())
        cov[x] = dict(gated_genes=int(gate.sum()), arc_pairs=n, share_of_pairs=float(n / (len(tg) * gate.sum())))
    info['DEF_arc_h1'] = dict(targets=len(tg), held=len(held), held_targets=held, **cov)
    save_json(TAB / 'partB_pools.json', info)
    log(f'B: Arc H1 holds {len(held)} of {len(tg)} D/E/F targets; pair shares {[cov[x]["share_of_pairs"] for x in "DEF"]}')


class BLine:
    """A Part B pool scored as a fold pool: P2.Fold + S16.Mag rebuilt from tables/partB_pool_<line>_w<w>.npz, config 1'."""

    def __init__(self, line, w):
        import scorer_part2 as P2
        import scorer_arms as SA
        import step15
        import scorer_study as SS
        self.line = line; self.ref = P2.refs_for(line); h = self.ref.h
        z = np.load(TAB / f'partB_pool_{line}_w{w}.npz', allow_pickle=True)
        assert list(z['targets']) == h.targets
        un = z['n_sources'] <= 0
        lfc = np.where(un, 0.0, z['lfc'].astype(np.float64)); self.sig = z['de_sig'] & ~un
        self.m0 = z['m0'].astype(np.float64); cprof = z['ctrl_profile'].astype(np.float64)
        dn = z['dpool_n']; dpool = np.where(dn > 0, z['dpool'].astype(np.float64), 0.0)
        n = z['n_bg']; za = z['z_agree'].astype(np.float64); has2 = (n >= 2) & np.isfinite(za)
        p = step15.CONFIGS['config1p']
        z2 = np.nan_to_num(za) ** 2
        c_non = np.where(has2, z2 / (z2 + p['kappa']), np.where(n == 1, p['c1'], 0.0))
        self.mu = np.where(self.sig, SA.mu_from_lfc(p['s_c'] * lfc, self.m0),
                           SA.mu_from_profile(cprof[None, :] + np.where(dn > 0, p['m'] * c_non * dpool, 0.0)))
        self.held = z['held']
        self.parts = {'full': np.arange(h.T), 'held_by_upweighted': np.flatnonzero(self.held), 'not_held': np.flatnonzero(~self.held)}


def partBscore(a):
    rows = []
    for line in B_POOLS:
        for w in B_W:
            D = BLine(line, w)
            rr, _ = score_rows(D, D.mu, w=w)
            rows += rr
            if line == 'K562' and w == 1:
                want = S19.d_base_reference('K562'); d = abs(rr[0]['overall'] - want)
                log(f'B score check: K562 w 1 {rr[0]["overall"]:.6f} vs the fold pool (15a) {want:.6f} |diff| {d:.1e}')
                if d > 1e-6:
                    log('B SCORE CHECK FAILED'); sys.exit(4)
            log(f'B score {line} w {w}: ' + ', '.join(f'{r["part"]} {r["overall"]:.4f} (PDS {r["pds"]:.3f})' for r in rr))
            del D; gc.collect()
        pd.DataFrame(rows).to_csv(TAB / 'partB_rows.csv', index=False)
    R = pd.DataFrame(rows)
    # the fold pools at w 1 for reference (K562's equals the K562 pool at w 1)
    ref = {l: S19.d_base_reference(l) for l in B_POOLS}
    res = dict(rows=R.to_dict('records'), fold_pool_base=ref, criterion={})
    ok = {}
    for line in B_POOLS:
        f = R[(R.line == line) & (R.part == 'full')].sort_values('w')
        ov, pds = f.overall.to_numpy(), f.pds.to_numpy()
        ok[line] = dict(overall_rises_with_w=bool(np.all(np.diff(ov) > 0)), pds_never_falls=bool(np.all(np.diff(pds) >= 0)),
                        overall_by_w=dict(zip(f.w.astype(int), ov)), pds_by_w=dict(zip(f.w.astype(int), pds)))
    res['criterion'] = dict(per_target=ok, met=bool(all(v['overall_rises_with_w'] and v['pds_never_falls'] for v in ok.values())))
    save_json(TAB / 'partB.json', res)
    with pd.option_context('display.width', 250, 'display.max_rows', 200):
        print(R.round(4).to_string())
    print(json.dumps(res['criterion'], indent=1, default=float))


# ================================================================ Part C: source-derived strong-minus-weak shift
C_AW, C_AM, C_AS = [-0.25, 0.0], [0.0, 0.25, 0.5, 1.0], [0.0, 0.25, 0.5, 1.0, 1.5]
C_GRID = [(w, m, s) for w in C_AW for m in C_AM for s in C_AS]
C_EXCL = {'H1': ('h1',), 'K562': ('k562',), 'HCT116': ('hct116', 'hek293t'), 'HEK293T': ('hct116', 'hek293t'), 'DEF': ()}


def strong_minus_weak(X, gated, n_real):
    """One line's vector: mean Delta over its strong knockdowns (n_real tercile 3) minus the mean over its weak ones
    (tercile 1), per gene over the knockdowns that measure it (gated, finite). NaN where either side has none."""
    t = tercile_of(n_real)
    out = []
    for j in (2, 0):
        Y = X[t == j]; M = np.isfinite(Y) & gated[t == j]
        out.append(np.where(M.sum(0) > 0, np.where(M, Y, 0.0).sum(0) / np.maximum(M.sum(0), 1), np.nan))
    return out[0] - out[1]


def line_vectors():
    """Per genetic full-profile source line (stage1 background key): tables averaged per knockdown first (Delta and
    n_real, each table's own > 5 CPM gate, own gene out), then the line's strong-minus-weak vector by its own n_real
    terciles. Rows without n_real (the Orion pseudobulk complements) are left out. Cached in tables/partC_lines.npz."""
    import stage1 as S
    import run_line as RL
    f = TAB / 'partC_lines.npz'
    if f.exists():
        z = np.load(f, allow_pickle=True); return z['lines'].tolist(), z['vec'], json.loads(str(z['info']))
    genes = RL.panel_genes(); gp = {g: j for j, g in enumerate(genes)}
    by_line = {}
    for p, stem, o in genetic_tables():
        by_line.setdefault(S._line_key(str(o['background'].iloc[0])), []).append((p, stem, o))
    lines, vecs, info = [], [], {}
    for ln, tabs in sorted(by_line.items()):
        acc = {}
        used = []
        for p, stem, o in tabs:
            nr = o['n_real'].astype(float).to_numpy()
            want = set(o['target_gene'].astype(str).to_numpy()[np.isfinite(nr)])
            if not want:
                continue
            d = read_rows(p, o, want, ('delta',))
            if 'delta' not in d:
                continue
            X = d['delta'].copy(); X[:, ~d['gate']] = np.nan
            for i, k in enumerate(d['kd']):
                if k in gp: X[i, gp[k]] = np.nan
                acc.setdefault(k, []).append((X[i], d['n_real'][i]))
            used.append(stem)
        if not acc:
            info[ln] = dict(tables=[s for _, s, _ in tabs], used=[], note='no full profile or no n_real: left out'); continue
        kd = sorted(acc)
        Xl = np.stack([np.nanmean(np.stack([r for r, _ in acc[k]]), 0) if len(acc[k]) > 1 else acc[k][0][0] for k in kd])
        nr = np.array([np.nanmean([n for _, n in acc[k]]) for k in kd])
        v = strong_minus_weak(Xl, np.isfinite(Xl), nr)
        lines.append(ln); vecs.append(v.astype(np.float32))
        info[ln] = dict(tables=[s for _, s, _ in tabs], used=used, knockdowns=len(kd), per_tercile=int(len(kd) // 3),
                        median_n_real=float(np.median(nr)), genes_defined=int(np.isfinite(v).sum()))
        log(f'C line {ln}: {len(kd)} knockdowns from {used}; vector defined on {np.isfinite(v).sum()} genes')
        del Xl, acc; gc.collect()
    vec = np.stack(vecs)
    np.savez_compressed(f, lines=np.array(lines), vec=vec, info=json.dumps(info), genes=np.array(genes))
    return lines, vec, info


def shift_S(lines, vec, exclude):
    """S: the mean over lines (excluded lines out) of their vectors, per gene over the lines that define it."""
    keep = np.array([l not in exclude for l in lines])
    V_ = vec[keep]; M = np.isfinite(V_)
    with np.errstate(invalid='ignore'):
        return np.where(M.sum(0) > 0, np.where(M, V_, 0.0).sum(0) / np.maximum(M.sum(0), 1), np.nan)


def cosine(x, y, m=None):
    ok = np.isfinite(x) & np.isfinite(y) & (True if m is None else m)
    return float((x[ok] * y[ok]).sum() / np.sqrt((x[ok] ** 2).sum() * (y[ok] ** 2).sum()))


def mean_residual(D, mu, rows):
    """Mean over rows of true minus predicted scorer-space Delta (R - bulk(mu)), own gene 0."""
    import scorer_arms as SA
    h = D.ref.h; Rr = h.R - SA.bulk(mu)
    for i in h.rows:
        if h.own[i] >= 0: Rr[i, h.own[i]] = 0.0
    return Rr[rows].mean(0), Rr


def c_rows(D, mu0, S, b, cuts, grid, tag):
    t = b_terc(b, cuts); own = own_index(D.ref.h); rows = []
    for aw, am, as_ in grid:
        av = np.array([aw, am, as_])[t]
        mu = mu0 if (aw, am, as_) == (0.0, 0.0, 0.0) else mu_add_delta(mu0, av[:, None] * S[None, :], own)
        rr, _ = score_rows(D, mu, a_weak=aw, a_mid=am, a_strong=as_, cuts=tag)
        rows += rr
    return rows


C2_MIN = 30


def line_strength():
    """Per source line (stage1 background key; genetic, n_real known): knockdowns, the size of its n_real terciles, and
    the median n_real of its strong and weak terciles (n_real averaged over the line's tables first, as line_vectors)."""
    import stage1 as S
    by = {}
    for p, stem, o in genetic_tables():
        nr = o['n_real'].astype(float).to_numpy(); kd = o['target_gene'].astype(str).to_numpy(); ok = np.isfinite(nr)
        for k, n in zip(kd[ok], nr[ok]):
            by.setdefault(S._line_key(str(o['background'].iloc[0])), {}).setdefault(k, []).append(n)
    rows = []
    for l, d in sorted(by.items()):
        nr = np.array([np.mean(v) for v in d.values()]); t = tercile_of(nr)
        rows.append(dict(line=l, knockdowns=len(nr), per_tercile=int((t == 2).sum()), strong_median_n_real=float(np.median(nr[t == 2])),
                         weak_median_n_real=float(np.median(nr[t == 0]))))
    return pd.DataFrame(rows).set_index('line')


def shift_S2(lines, vec, exclude, keep_lines):
    """S over the kept lines (excluded lines out), only on genes that at least half of those lines define."""
    keep = np.array([l in keep_lines and l not in exclude for l in lines])
    V_ = vec[keep]; M = np.isfinite(V_); need = int(np.ceil(keep.sum() / 2))
    with np.errstate(invalid='ignore'):
        return np.where(M.sum(0) >= need, np.where(M, V_, 0.0).sum(0) / np.maximum(M.sum(0), 1), np.nan)


def partC2(a):
    partC(a, 'C2')


def partC(a, tag='C'):
    """tag C: S as specified. tag C2 (Jose 2026-11-02, the fix): only lines with >= 30 knockdowns per n_real tercile and
    a strong-tercile median n_real >= 30; S only on genes measured by at least half the voting lines; applied only on
    the scored line's gated genes."""
    import scorer_step12 as S12
    import scorer_step16 as S16
    import run_line as RL
    genes = np.array(RL.panel_genes())
    lines, vec, info = line_vectors()
    res = dict(lines=info, lines_used=lines, excluded={k: list(v) for k, v in C_EXCL.items()}, fits={}, guardrails={}, surface={})
    if tag == 'C':
        Ss = {k: shift_S(lines, vec, ex) for k, ex in C_EXCL.items()}
    else:
        st = line_strength()
        keep = [l for l in lines if l in st.index and st.loc[l, 'per_tercile'] >= C2_MIN and st.loc[l, 'strong_median_n_real'] >= C2_MIN]
        res['line_filter'] = dict(rule=f'>= {C2_MIN} knockdowns per tercile and strong-tercile median n_real >= {C2_MIN}', kept=keep,
                                  dropped=[l for l in lines if l not in keep], strength=st.reset_index().to_dict('records'))
        Ss = {k: shift_S2(lines, vec, ex, keep) for k, ex in C_EXCL.items()}
        log(f'C2: {len(keep)} lines kept {keep}; S defined on ' + ', '.join(f'{k} {np.isfinite(v).sum()}' for k, v in Ss.items()) + ' genes')
    use = (lambda S_, D_: S_) if tag == 'C' else (lambda S_, D_: np.where(D_.gate, S_, np.nan))
    np.savez_compressed(TAB / f'part{tag}_S.npz', genes=genes, **{f'S_{k}': v.astype(np.float32) for k, v in Ss.items()})
    f_rows = TAB / f'part{tag}_rows.csv'
    rows = pd.read_csv(f_rows).to_dict('records') if f_rows.exists() else []
    done = {(r['line'], r['cuts'], r['a_weak'], r['a_mid'], r['a_strong']) for r in rows}
    D = DLine('H1'); h = D.ref.h; mu0 = base_mu(D)
    T = S12.truth_calls(D.ref); terc = tercile_of(T.sum(1)); halves = S16.halves(T)
    bH = b_for('EvalA', h.targets)
    cuts = {'half0': b_cuts(bH[halves[0]]), 'half1': b_cuts(bH[halves[1]]), 'full': b_cuts(bH)}
    res['cuts'] = {k: v.tolist() for k, v in cuts.items()}
    res['b_tercile_counts_H1'] = {k: np.bincount(b_terc(bH, v), minlength=3).tolist() for k, v in cuts.items()}
    t0 = time.time()
    for ck in ('half0', 'half1', 'full'):
        todo = [g for g in C_GRID if ('H1', ck) + g not in done]
        for g in todo:
            rows += c_rows(D, mu0, use(Ss['H1'], D), bH, cuts[ck], [g], ck)
            if len(rows) % 30 == 0:
                pd.DataFrame(rows).to_csv(f_rows, index=False)
        pd.DataFrame(rows).to_csv(f_rows, index=False)
        log(f'C H1 cuts {ck}: done ({(time.time() - t0) / max(len(rows) // 3, 1):.1f} s per row so far)')
    # cosines and top genes
    gate = D.gate
    m_strong, Rr = mean_residual(D, mu0, np.flatnonzero(terc == 2))
    z = np.load(S19T / 'predside_K0_EvalA.npz'); Rc = h.R - z['P']
    for i in h.rows:
        if h.own[i] >= 0: Rc[i, h.own[i]] = 0.0
    res['cos_S_H1excl_vs_H1_strong_mean_residual'] = dict(proxy_base=cosine(Ss['H1'], m_strong, gate),
                                                          uploaded_cells=cosine(Ss['H1'], Rc[terc == 2].mean(0), gate),
                                                          genes=int((np.isfinite(Ss['H1']) & gate).sum()))
    E_ = np.nan_to_num(Ss['H1']) ** 2
    res['S_H1_energy_share_off_H1_gate'] = float(E_[~gate].sum() / E_.sum())      # C2 applies S on the gate only
    for k in ('H1', 'DEF'):
        v = np.where(np.isfinite(Ss[k]), Ss[k], 0.0); o = np.argsort(v)
        res[f'top30_S_{k}'] = dict(up=[(genes[i], float(v[i])) for i in o[::-1][:30]], down=[(genes[i], float(v[i])) for i in o[:30]])
    # the follow-up's H1-own shift on the proxy base (fit-half true-tercile mean residual, mid and strong, by (b) tercile)
    bt = np.full(h.T, 0); ok = np.isfinite(bH); bt[ok] = tercile_of(bH[ok])
    own_rows = []
    for f in (0, 1):
        fit, ho = halves[f], halves[1 - f]; add = np.zeros_like(Rr)
        for j in (1, 2):
            add[np.intersect1d(ho, np.flatnonzero(bt == j))] = Rr[np.intersect1d(fit, np.flatnonzero(terc == j))].mean(0)
        rr, _ = score_rows(D, mu_add_delta(mu0, add, own_index(h)), arm='H1-own shift (follow-up), proxy', fit_half=f)
        own_rows += rr
    res['H1_own_shift_proxy'] = own_rows
    del D; gc.collect()
    R = pd.DataFrame(rows)
    H = R[R.line == 'H1']
    isbase = lambda d: (d.a_weak == 0) & (d.a_mid == 0) & (d.a_strong == 0)
    key = ['a_weak', 'a_mid', 'a_strong']
    picks = {}
    for f in (0, 1):
        got = pick_fit(H[H.cuts == f'half{f}'], f'half{f}', f'half{1 - f}', key, isbase)
        if got is None:
            res['fits'][f'fit_half{f}'] = dict(picked=None); continue
        pk, ho, b0, ne = got
        picks[f] = tuple(float(pk[k]) for k in key)
        res['fits'][f'fit_half{f}'] = dict(picked=dict(zip(key, picks[f])), eligible=ne, fit_overall=float(pk.overall),
                                          grid_position={k: grid_pos(g, pk[k]) for k, g in zip(key, (C_AW, C_AM, C_AS))},
                                          held_out=comp_dict(ho), base_held_out=comp_dict(b0),
                                          gain={c: float(ho[c] - b0[c]) for c in COMPS + ['overall']})
        log(f'C fit half{f}: pick {picks[f]} held-out {ho.overall:.4f} vs base {b0.overall:.4f}')
    F_ = H[(H.cuts == 'full') & (H.part == 'full')]
    res['surface'] = dict(base=comp_dict(F_[isbase(F_)].iloc[0]), best=F_.sort_values('overall').iloc[-1][key + COMPS + ['overall']].to_dict(),
                          best_all_positive=(lambda d: d.sort_values('overall').iloc[-1][key + COMPS + ['overall']].to_dict() if len(d) else None)(F_[(F_[COMPS] > 0).all(1)]),
                          marginal_mean_overall={k: F_.groupby(k).overall.mean().to_dict() for k in key},
                          marginal_mean_pds={k: F_.groupby(k).pds.mean().to_dict() for k in key},
                          marginal_mean_mse={k: F_.groupby(k).mse.mean().to_dict() for k in key})
    # guardrails: the full grid with all-H1 cut points, and the picks with their fit half's cut points
    for line in D_LINES[1:]:
        D = DLine(line); mu0 = base_mu(D); b = b_for(B_LINE[line], D.ref.h.targets)
        todo = [g for g in C_GRID if (line, 'full') + g not in done]
        for g in todo:
            rows += c_rows(D, mu0, use(Ss[line], D), b, cuts['full'], [g], 'full')
        for f, g in picks.items():
            if (line, f'half{f}') + g not in done:
                rows += c_rows(D, mu0, use(Ss[line], D), b, cuts[f'half{f}'], [(0.0, 0.0, 0.0), g], f'half{f}')
        pd.DataFrame(rows).to_csv(f_rows, index=False)
        res['b_tercile_counts'] = res.get('b_tercile_counts', {}) | {line: np.bincount(b_terc(b, cuts['full']), minlength=3).tolist()}
        log(f'C {line}: rows done'); del D; gc.collect()
    R = pd.DataFrame(rows)
    for f, g in picks.items():
        gr = {}
        for line in D_LINES[1:]:
            L = R[(R.line == line) & (R.cuts == f'half{f}')]
            x = L[(L.a_weak == g[0]) & (L.a_mid == g[1]) & (L.a_strong == g[2])].iloc[0]; y = L[isbase(L)].iloc[0]
            gr[line] = dict(pick=comp_dict(x), base=comp_dict(y), diff=float(x.overall - y.overall))
        res['fits'][f'fit_half{f}']['guardrails'] = gr
    for line in D_LINES[1:]:
        L = R[(R.line == line) & (R.cuts == 'full')]
        res['guardrails'][line] = dict(base=comp_dict(L[isbase(L)].iloc[0]), best=L.sort_values('overall').iloc[-1][key + ['overall', 'pds']].to_dict(),
                                       marginal_mean_overall={k: L.groupby(k).overall.mean().to_dict() for k in key})
    fv = [v for v in res['fits'].values() if v.get('picked')]
    res['criterion'] = dict(both_picked=len(fv) == 2,
                            heldout_gain_ge_0005_both=bool(len(fv) == 2 and all(v['gain']['overall'] >= 0.005 for v in fv)),
                            all_components_pos_both=bool(len(fv) == 2 and all(all(v['held_out'][c] > 0 for c in COMPS) for v in fv)),
                            guardrails_within_0005=bool(len(fv) == 2 and all(g['diff'] >= -0.005 for v in fv for g in v['guardrails'].values())))
    res['criterion']['met'] = all(res['criterion'].values())
    res['picks'] = {f: list(g) for f, g in picks.items()}
    save_json(TAB / f'part{tag}.json', res)
    print(json.dumps({k: v for k, v in res.items() if k not in ('lines',)}, indent=1, default=float))


# ================================================================ Part D: H1-derived shift, pluripotency removed
D_PROGRAMS = {   # the follow-up's residual PCs (fit half, PC) that carry the three named generic programs
    'p53 / arrest vs S-phase': (0, 2), 'amino-acid stress': (0, 4), 'sterol synthesis': (1, 4)}


def residual_pcs(Rr, rows, n=5):
    _, s, W = np.linalg.svd(Rr[rows], full_matrices=False)
    return W[:n], s


def unit(v):
    v = np.where(np.isfinite(v), v, 0.0); return v / np.linalg.norm(v)


def partD(a):
    import scorer_step12 as S12
    import scorer_step16 as S16
    import run_line as RL
    genes = np.array(RL.panel_genes())
    D = DLine('H1'); h = D.ref.h; mu0 = base_mu(D)
    T = S12.truth_calls(D.ref); terc = tercile_of(T.sum(1)); halves = S16.halves(T)
    _, Rr = mean_residual(D, mu0, np.arange(h.T))
    S_H1 = Rr[terc == 2].mean(0) - Rr[terc == 0].mean(0)
    # the follow-up's residual programs, from the uploaded cells (P - R), mid + strong, as mseterc computed them
    z = np.load(S19T / 'predside_K0_EvalA.npz'); Rc = z['P'] - h.R
    for i in h.rows:
        if h.own[i] >= 0: Rc[i, h.own[i]] = 0.0
    ms = terc >= 1
    J = json.loads((S19T / 'mse_terciles.json').read_text())['own_pcs']
    W_half, chk = {}, {}
    for f in (0, 1):
        W, _ = residual_pcs(Rc, np.intersect1d(halves[f], np.flatnonzero(ms)))
        W_half[f] = W
        for c in range(5):
            o = np.argsort(W[c]); want = J[f'fit_half{f}']['pcs'][c]
            chk[f'half{f}_pc{c + 1}'] = bool(genes[o[::-1][:15]].tolist() == want['pos_genes'] and genes[o[:15]].tolist() == want['neg_genes'])
    W_all, _ = residual_pcs(Rc, np.flatnonzero(ms))
    pc1 = unit(W_all[0])
    res = dict(reproduces_followup_pcs=chk, pc1_full_vs_half_pcs1={f'half{f}': float(abs(pc1 @ W_half[f][0])) for f in (0, 1)})
    log(f'D: follow-up PCs reproduced {chk}; full-panel PC1 |cos| with the halves\' PC1 {res["pc1_full_vs_half_pcs1"]}')
    if not all(chk.values()):
        save_json(TAB / 'partD.json', res); log('D CHECK FAILED: the follow-up residual PCs do not reproduce'); sys.exit(5)
    o = np.argsort(pc1); res['pc1_genes'] = dict(pos=genes[o[::-1][:15]].tolist(), neg=genes[o[:15]].tolist())
    S_clean = S_H1 - (S_H1 @ pc1) * pc1
    np.savez_compressed(TAB / 'partD_S_H1_clean.npz', genes=genes, S_H1=S_H1.astype(np.float32), S_H1_clean=S_clean.astype(np.float32),
                        pc1=pc1.astype(np.float32))
    Sz = np.load(TAB / 'partC_S.npz')
    lines, vec, _ = line_vectors()
    S_all = shift_S(lines, vec, ()); S_noH1 = Sz['S_H1']
    gate = D.gate
    res['pc1_share_of_S_H1_energy'] = float((S_H1 @ pc1) ** 2 / (S_H1 @ S_H1))
    res['cos'] = dict(S_H1clean_vs_S_all_lines=cosine(S_clean, S_all, gate), S_H1clean_vs_S_without_H1=cosine(S_clean, S_noH1, gate),
                      S_H1_vs_S_all_lines=cosine(S_H1, S_all, gate), S_H1_vs_S_without_H1=cosine(S_H1, S_noH1, gate))
    progs = {nm: unit(W_half[f][c - 1]) for nm, (f, c) in D_PROGRAMS.items()}
    res['program_vectors'] = {nm: dict(fit_half=f, pc=c, pos=genes[np.argsort(W_half[f][c - 1])[::-1][:10]].tolist(),
                                       neg=genes[np.argsort(W_half[f][c - 1])[:10]].tolist()) for nm, (f, c) in D_PROGRAMS.items()}
    res['energy_share'] = {}
    for vn, v in (('S_H1', S_H1), ('S_H1_clean', S_clean), ('S_all_lines', S_all), ('S_without_H1', S_noH1)):
        v0 = np.where(np.isfinite(v), v, 0.0)
        res['energy_share'][vn] = {nm: float((v0 @ u) ** 2 / (v0 @ v0)) for nm, u in progs.items()}
        res['energy_share'][vn]['pluripotency PC1'] = float((v0 @ pc1) ** 2 / (v0 @ v0))
    o = np.argsort(S_clean)
    res['top30_S_H1_clean'] = dict(up=[(genes[i], float(S_clean[i])) for i in o[::-1][:30]], down=[(genes[i], float(S_clean[i])) for i in o[:30]])
    del D; gc.collect()
    # guardrails: S_H1,clean by (b) tercile at Part C's picks (each with its fit half's cut points)
    C = json.loads((TAB / 'partC.json').read_text())
    bH = b_for('EvalA', h.targets); res['guardrails'] = {}
    rows = []
    for f, g in C['picks'].items():
        cuts = b_cuts(bH[halves[int(f)]])
        for line in D_LINES[1:]:
            Dl = DLine(line); m0 = base_mu(Dl); b = b_for(B_LINE[line], Dl.ref.h.targets)
            rr = c_rows(Dl, m0, S_clean, b, cuts, [(0.0, 0.0, 0.0), tuple(g)], f'half{f}')
            for r in rr: r['pick_from_fit_half'] = int(f)
            rows += rr
            x, y = rr[1], rr[0]
            res['guardrails'].setdefault(f'pick_fit_half{f}', {})[line] = dict(a=g, pick={c: x[c] for c in COMPS + ['overall']},
                                                                              base={c: y[c] for c in COMPS + ['overall']}, diff=x['overall'] - y['overall'])
            log(f'D guardrail {line} pick {g}: {x["overall"]:.4f} vs base {y["overall"]:.4f} (PDS {x["pds"]:.3f} vs {y["pds"]:.3f})')
            del Dl; gc.collect()
    pd.DataFrame(rows).to_csv(TAB / 'partD_rows.csv', index=False)
    save_json(TAB / 'partD.json', res)
    print(json.dumps(res, indent=1, default=float))



# ================================================================ Part F (vcc2027): program regression v2
F_BASES = {'E2w': ('lfc', True), 'E4w': ('delta', False), 'E24w': ('delta', True)}       # (space, n_real >= 50 rows)
F_FOLDS = {'full': (), 'K562': ('K562',), 'Orion': ('HCT116', 'HEK293T')}
FOLD_OF = {'H1': 'full', 'K562': 'K562', 'HCT116': 'Orion', 'HEK293T': 'Orion'}
OWN_BG_F = {'H1': ('H1',), 'K562': ('K562',), 'HCT116': ('HCT116', 'HEK293T'), 'HEK293T': ('HCT116', 'HEK293T')}
MIN_BG = 50
DL = V / 'lib/downloads/step20_2026-11-02'
COMPLEXES = ['Proteasome', 'Mediator complex', 'BAF complex subunits', 'L ribosomal proteins', 'S ribosomal proteins',
             'Integrator complex', 'Exosome complex', 'COP9 signalosome', 'Chaperonin containing TCP1',
             'NADH:ubiquinone oxidoreductase core subunits']


def f_weight(n_b):
    """Step 20's weighting fix: row weight 1 / sqrt(max(n_b, 50))."""
    return 1.0 / np.sqrt(np.maximum(n_b, MIN_BG))


def f_build(variant):
    """Basis matrix (weighted rows) and row metadata for E2w / E4w / E24w: Step 19's e_build rows and entries (5 CPM gate
    in each table's own controls, per-table centring over measured entries), Tahoe as one background, Arc H1 out, with
    row weight 1 / sqrt(max(n_b, 50))."""
    import h5py
    import run_line as RL
    from anndata.io import read_elem
    space, filt = F_BASES[variant]
    genes = RL.panel_genes(); G = len(genes)
    plan = []
    for p, stem, bg, tahoe in S19.e_tables():
        with h5py.File(p, 'r') as f:
            o = read_elem(f['obs'])
            if space == 'delta' and not np.isfinite(f['layers/log1p_bulk_50k_full'][:min(5, len(o))]).any():
                continue
        nr = o['n_real'].astype(float).to_numpy()
        keep = (np.isfinite(nr) & (nr >= 50)) if filt else np.ones(len(o), bool)
        if keep.any():
            plan.append((p, stem, bg, tahoe, o, np.flatnonzero(keep)))
    n_bg = pd.Series({}, dtype=float)
    for _, _, bg, _, _, k in plan:
        n_bg[bg] = n_bg.get(bg, 0) + len(k)
    n_rows = int(n_bg.sum())
    mem_before(f'{variant} basis {n_rows} x {G}', n_rows * G * 4 / 1e9 + 3)
    M = np.zeros((n_rows, G), np.float32); meta = []; r0 = 0
    for p, stem, bg, tahoe, o, k in plan:
        X, Mm = table_entries(p, k, space)
        mean = np.where(Mm, X, 0.0).sum(0) / np.maximum(Mm.sum(0), 1)
        w = f_weight(n_bg[bg])
        M[r0:r0 + len(k)] = np.where(Mm, X - mean[None, :], 0.0) * w; r0 += len(k)
        meta.append(pd.DataFrame(dict(table=stem, background=bg, tahoe=tahoe, modality=o['modality'].astype(str).to_numpy()[k],
                                      target=o['target_gene'].astype(str).to_numpy()[k], n_real=o['n_real'].astype(float).to_numpy()[k], weight=w)))
        del X, Mm
    mem_after(f'{variant} basis')
    return M, pd.concat(meta, ignore_index=True), n_bg


def table_entries(p, k, space, gate_cpm=5.0):
    """Raw entries (before centring) for rows k of a table, and the measured mask (Step 19's e_entries)."""
    import h5py
    import run_line as RL
    genes = RL.panel_genes()
    b = pd.read_csv(p.with_name(p.stem + '_baseline.csv'), index_col=0).reindex(genes)
    gated = b['control_mean_cpm'].to_numpy(float) > gate_cpm
    ctrl = b['log1p_bulk_50k_full'].to_numpy(float) if space == 'delta' else None
    blocks, meass = [], []
    with h5py.File(p, 'r') as f:
        for c0 in range(0, len(k), 2000):
            kk = k[c0:c0 + 2000]; sl = slice(kk[0], kk[-1] + 1)
            if space == 'lfc':
                x, m = S19.e_entries(f['layers/lfc'][sl].astype(np.float64)[kk - kk[0]], f['layers/se'][sl].astype(np.float64)[kk - kk[0]], gated, 'lfc')
            else:
                x, m = S19.e_entries(None, None, gated, 'delta', f['layers/log1p_bulk_50k_full'][sl].astype(np.float64)[kk - kk[0]], ctrl)
            blocks.append(x); meass.append(m)
    return np.concatenate(blocks), np.concatenate(meass)


def f_source_rows(space):
    """Every genetic row (Tahoe drugs and Arc H1 out) in a basis space, unweighted, centred per table over all its rows:
    the rows projected for the source scores s_hat."""
    import h5py
    import run_line as RL
    from anndata.io import read_elem
    plan = []
    for p, stem, bg, tahoe in S19.e_tables():
        if tahoe:
            continue
        with h5py.File(p, 'r') as f:
            o = read_elem(f['obs'])
            if space == 'delta' and not np.isfinite(f['layers/log1p_bulk_50k_full'][:min(5, len(o))]).any():
                continue
        g = np.flatnonzero(o['modality'].astype(str).to_numpy() != 'drug')
        if len(g):
            plan.append((p, stem, bg, o, g))
    n = sum(len(g) for *_, g in plan); G = len(RL.panel_genes())
    mem_before(f'source rows {space} {n} x {G}', n * G * 4 / 1e9 + 3)
    M = np.zeros((n, G), np.float32); meta = []; r0 = 0
    for p, stem, bg, o, g in plan:
        X, Mm = table_entries(p, g, space)
        mean = np.where(Mm, X, 0.0).sum(0) / np.maximum(Mm.sum(0), 1)
        M[r0:r0 + len(g)] = np.where(Mm, X - mean[None, :], 0.0); r0 += len(g)
        meta.append(pd.DataFrame(dict(table=stem, background=bg, target=o['target_gene'].astype(str).to_numpy()[g])))
    mem_after(f'source rows {space}')
    return M, pd.concat(meta, ignore_index=True)


def arc_rows(space):
    """Arc H1's rows in a basis space (own 5 CPM gate, centred over its rows), averaged per target."""
    import h5py
    from anndata.io import read_elem
    p = V / 'lib/effects_v5/arc_h1.h5ad'
    with h5py.File(p, 'r') as f:
        o = read_elem(f['obs'])
    X, Mm = table_entries(p, np.arange(len(o)), space)
    X = np.where(Mm, X - (np.where(Mm, X, 0.0).sum(0) / np.maximum(Mm.sum(0), 1))[None, :], 0.0)
    u, Y = collapse(o['target_gene'].astype(str).to_numpy(), X)
    return u, Y


def basis_file(b, fold='full'):
    if b in ('E2', 'E4'):
        return S19T / f'partE_{b}_basis.npz'
    return TAB / f'partF1_{b}_{fold}_basis.npz'


def k25_of(b, fold='full'):
    z = np.load(basis_file(b, fold)); cum = z['cum']
    return int(np.searchsorted(cum, 0.25) + 1) if cum[-1] >= 0.25 else None


def small_bg_artefacts(US, meta, n_bg, top=10):
    """Per top component: share of its row-score energy from backgrounds with < 50 rows, and how many of its 12 largest
    |row scores| come from them."""
    small = meta.background.map(lambda b: n_bg.get(b, 0) < MIN_BG).to_numpy()
    out = []
    for c in range(top):
        e = US[:, c].astype(np.float64) ** 2
        o = np.argsort(-np.abs(US[:, c]))[:12]
        out.append(dict(comp=c + 1, small_bg_energy_share=float(e[small].sum() / e.sum()), small_bg_in_top12=int(small[o].sum())))
    return out


def partF1(a):
    import run_line as RL
    genes = np.array(RL.panel_genes())
    res = dict(bases={}, step19_artefacts={})
    for variant, (space, filt) in F_BASES.items():
        M, meta, n_bg = f_build(variant)
        for fold, drop in F_FOLDS.items():
            keep = ~meta.background.isin(drop).to_numpy()
            Mf = M[keep] if drop else M
            mf = meta[keep].reset_index(drop=True)
            US, s, Vt, cum, k25 = S19.e_svd(Mf, f'{variant}-{fold}')
            np.savez_compressed(basis_file(variant, fold), Vt=Vt, s=s, cum=cum, US=US)
            mf.to_parquet(TAB / f'partF1_{variant}_{fold}_rows.parquet')
            if fold != 'full':
                res['bases'][variant][f'fold_{fold}'] = dict(rows=int(len(mf)), K_at_25pct=k25, cum_share={K: float(cum[K - 1]) for K in (5, 10, 20, 50, 100)})
                del Mf; gc.collect(); continue
            ew = (meta.weight.to_numpy() ** 2); tot = ew.sum()
            nb = meta.groupby('background').size()
            small = [b for b, n in nb.items() if n < MIN_BG]
            comps = []
            for c in range(10):
                o = np.argsort(Vt[c]); sc = US[:, c]
                comps.append(dict(comp=c + 1, var_share=float(cum[c] - (cum[c - 1] if c else 0.0)),
                                  pos_genes=genes[o[::-1][:20]].tolist(), neg_genes=genes[o[:20]].tolist(),
                                  top_pos_rows=[f'{mf.table[i]}:{mf.target[i]}' for i in np.argsort(-sc)[:10]],
                                  top_neg_rows=[f'{mf.table[i]}:{mf.target[i]}' for i in np.argsort(sc)[:10]]))
            res['bases'][variant] = dict(space=space, n_real_ge_50=filt, rows=int(len(meta)), tables=int(meta.table.nunique()),
                                         rows_per_background={k: int(v) for k, v in nb.items()},
                                         backgrounds_under_50=small, weight_share_under_50=float(ew[meta.background.isin(small)].sum() / tot),
                                         weight_share_by_background={k: float(ew[(meta.background == k).to_numpy()].sum() / tot) for k in nb.index},
                                         cum_share={K: float(cum[K - 1]) for K in (5, 10, 20, 50, 100)}, K_at_25pct=k25,
                                         components=comps, artefacts=small_bg_artefacts(US, meta, nb))
            log(f'F1 {variant}: {len(meta)} rows; under-50 backgrounds {small} hold {res["bases"][variant]["weight_share_under_50"]:.4f} of the weight; '
                f'K at 25% {k25}')
            save_json(TAB / 'partF1.json', res)
        del M; gc.collect()
    for b in ('E1', 'E2', 'E4'):
        z = np.load(S19T / f'partE_{b}_basis.npz'); m = pd.read_parquet(S19T / f'partE_{b}_rows.parquet')
        res['step19_artefacts'][b] = small_bg_artefacts(z['US'], m, m.groupby('background').size())
    save_json(TAB / 'partF1.json', res)
    # s_hat inputs: every genetic row projected on each basis (and fold basis), Step 19's E2 / E4 for comparison
    for space in ('lfc', 'delta'):
        Ms, ms = f_source_rows(space)
        ms.to_parquet(TAB / f'partF1_srcrows_{space}.parquet')
        targets = {}
        for b, (sp, _) in list(F_BASES.items()) + [('E2', ('lfc', 0)), ('E4', ('delta', 0))]:
            if sp != space:
                continue
            for fold in (F_FOLDS if b in F_BASES else ['full']):
                Vt = np.load(basis_file(b, fold))['Vt'].astype(np.float32)
                RS = np.empty((len(ms), Vt.shape[0]), np.float32)
                for r0 in range(0, len(ms), 10000):
                    RS[r0:r0 + 10000] = Ms[r0:r0 + 10000] @ Vt.T
                np.save(TAB / f'partF1_RS_{b}_{fold}.npy', RS); del RS
        au, AX = arc_rows(space)
        for b, (sp, _) in list(F_BASES.items()) + [('E2', ('lfc', 0)), ('E4', ('delta', 0))]:
            if sp == space:
                Vt = np.load(basis_file(b, 'full'))['Vt'].astype(np.float32)
                np.savez_compressed(TAB / f'partF1_arc_{b}.npz', targets=au, scores=(AX @ Vt.T).astype(np.float32))
        del Ms; gc.collect(); log(f'F1 source rows {space}: projected')
    print(json.dumps(res, indent=1, default=float)[:20000])


def partF2(a):
    tr = S19.e_truth(); rows = []
    for b, (space, _) in F_BASES.items():
        Vt = np.load(basis_file(b))['Vt']; k25 = k25_of(b)
        for K in sorted(set([5, 10, 20, 50, min(k25 or 100, 100)])):
            sums = S19.e_ceiling(Vt, K, tr, space)
            r = dict(basis=b, K=K, K_is_25pct=K == min(k25 or 100, 100), scale='LFC (a, b), Delta (c)' if space == 'lfc' else 'Delta (c) only')
            for m in ('abc' if space == 'lfc' else 'c'):
                r.update({f'{m}_{g}': v for g, v in S19.captured(sums[m], tr['terc']).items()})
            rows.append(r); log(f'F2 {b} K {K}: {r}')
    C = pd.read_csv(S19T / 'partE_ceiling.csv')
    old = C[C.basis.isin(['E1', 'E2', 'E4', 'config1p_prediction'])].copy()
    old['scale'] = np.where(old.basis == 'E4', 'Delta (c); Step 19 (b) on Delta not shown', 'LFC (a, b), Delta (c)')
    old.loc[old.basis == 'E4', [c for c in old.columns if c.startswith(('a_', 'b_'))]] = np.nan
    R = pd.concat([pd.DataFrame(rows), old.assign(step19=True)], ignore_index=True)
    R.to_csv(TAB / 'partF2_ceiling.csv', index=False)
    with pd.option_context('display.width', 250, 'display.max_rows', 200):
        print(R.round(4).to_string())


# ---------------------------------------------------------------- F3: gene features
def hgnc_maps():
    t = pd.read_csv(V / 'lib/downloads/hgnc_2026-10-28/hgnc_complete_set.txt', sep='\t',
                    usecols=['symbol', 'prev_symbol', 'alias_symbol', 'entrez_id', 'gene_group'], dtype=str)
    appr = set(t.symbol)
    prev, alias = {}, {}
    for r in t.itertuples():
        for col, d in (('prev_symbol', prev), ('alias_symbol', alias)):
            for x in str(getattr(r, col)).split('|'):
                if x and x != 'nan':
                    d.setdefault(x, set()).add(r.symbol)
    ent = {str(e).split('.')[0]: s for e, s in zip(t.entrez_id, t.symbol) if isinstance(e, str)}
    groups = {r.symbol: [g for g in str(r.gene_group).split('|') if g and g != 'nan'] for r in t.itertuples()}
    return dict(appr=appr, prev={k: next(iter(v)) for k, v in prev.items() if len(v) == 1},
                alias={k: next(iter(v)) for k, v in alias.items() if len(v) == 1}, entrez=ent, groups=groups)


def to_hgnc(sym, H):
    """A symbol in its approved HGNC form: itself if approved, else a unique previous symbol, else a unique alias."""
    if sym in H['appr']:
        return sym
    return H['prev'].get(sym) or H['alias'].get(sym)


def download(url, dest, md5=None):
    """curl to dest (skipped if present with the right md5); logs URL, size and md5."""
    import hashlib
    import subprocess
    dest = Path(dest); dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        log(f'download {url} -> {dest}')
        r = subprocess.run(['curl', '-sSfL', '--retry', '3', '-o', str(dest) + '.part', url])
        if r.returncode != 0:
            raise RuntimeError(f'download failed ({r.returncode}): {url}')
        os.replace(str(dest) + '.part', dest)
    h = hashlib.md5()
    with open(dest, 'rb') as f:
        for ch in iter(lambda: f.read(1 << 24), b''):
            h.update(ch)
    got = h.hexdigest()
    info = dict(url=url, file=str(dest), bytes=dest.stat().st_size, md5=got, md5_published=md5, md5_ok=(md5 is None or md5 == got))
    log(f'download ok: {info}')
    if md5 and md5 != got:
        raise RuntimeError(f'md5 mismatch for {dest}: {got} vs published {md5}')
    return info


def figshare_file(article, name):
    import urllib.request
    a = json.loads(urllib.request.urlopen(f'https://api.figshare.com/v2/articles/{article}', timeout=60).read())
    f = next(x for x in a['files'] if x['name'] == name)
    return a['title'], a['published_date'], f['download_url'], f.get('computed_md5') or f.get('supplied_md5'), f['size']


def partF3a(a):
    """DepMap CRISPRGeneEffect.csv: the portal's download API answers with a bot-verification page (checked 2026-11-02),
    so the file comes from DepMap's figshare release that carries it: DepMap 24Q4 Public (article 27993248)."""
    from sklearn.decomposition import PCA
    title, date, url, md5, size = figshare_file(27993248, 'CRISPRGeneEffect.csv')
    info = download(url, DL / 'depmap_24Q4/CRISPRGeneEffect.csv', md5); info.update(release=title, published=date)
    H = hgnc_maps()
    E = pd.read_csv(DL / 'depmap_24Q4/CRISPRGeneEffect.csv', index_col=0)
    cols = []
    for c in E.columns:
        sym, ent = c.split(' (')[0], c.split('(')[-1].rstrip(')')
        cols.append(H['entrez'].get(ent) or to_hgnc(sym, H))
    E.columns = cols; E = E.loc[:, [c is not None for c in E.columns]]
    G = E.T.groupby(level=0).mean()                                         # genes x cell lines, duplicates averaged
    mean_effect = G.mean(1)
    X = G.T.fillna(G.mean(1)).T.to_numpy(float)                              # missing: the gene's mean
    X = (X - X.mean(0)) / X.std(0)                                           # standardise each cell line
    pca = PCA(32, random_state=SEED).fit(X)
    emb = pca.transform(X)
    np.savez_compressed(TAB / 'partF3_depmap.npz', genes=np.array(G.index), emb=emb.astype(np.float32), extra=mean_effect.to_numpy(np.float32),
                        var_share=pca.explained_variance_ratio_)
    out = dict(download=info, genes=int(len(G)), cell_lines=int(G.shape[1]), unmapped_columns=int(sum(c is None for c in cols)),
               var_share_32=float(pca.explained_variance_ratio_.sum()), var_share_each=pca.explained_variance_ratio_.tolist())
    save_json(TAB / 'partF3a.json', out); log(f'F3a: {out["genes"]} genes x {out["cell_lines"]} lines; 32 PCs hold {out["var_share_32"]:.3f}')


def partF3b(a):
    import gzip
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse.linalg import eigsh
    base = 'https://stringdb-downloads.org/download'
    i1 = download(f'{base}/protein.links.v12.0/9606.protein.links.v12.0.txt.gz', DL / 'string_v12/9606.protein.links.v12.0.txt.gz')
    i2 = download(f'{base}/protein.info.v12.0/9606.protein.info.v12.0.txt.gz', DL / 'string_v12/9606.protein.info.v12.0.txt.gz')
    H = hgnc_maps()
    info = pd.read_csv(DL / 'string_v12/9606.protein.info.v12.0.txt.gz', sep='\t', usecols=[0, 1])
    pmap = {p: to_hgnc(n, H) for p, n in zip(info.iloc[:, 0], info.iloc[:, 1])}
    L = pd.read_csv(DL / 'string_v12/9606.protein.links.v12.0.txt.gz', sep=' ')
    L = L[L.combined_score >= 700]
    L['g1'] = L.protein1.map(pmap); L['g2'] = L.protein2.map(pmap)
    L = L.dropna(subset=['g1', 'g2']); L = L[L.g1 != L.g2]
    a_, b_ = np.minimum(L.g1, L.g2), np.maximum(L.g1, L.g2)
    E = pd.DataFrame(dict(a=a_, b=b_, s=L.combined_score.to_numpy())).groupby(['a', 'b']).s.max().reset_index()  # one edge per gene pair, max score
    genes = np.array(sorted(set(E.a) | set(E.b))); gi = {g: i for i, g in enumerate(genes)}
    r, c = E.a.map(gi).to_numpy(), E.b.map(gi).to_numpy()
    A = sp.coo_matrix((np.ones(len(r)), (r, c)), shape=(len(genes),) * 2).tocsr(); A = ((A + A.T) > 0).astype(np.float64)
    deg = np.asarray(A.sum(1)).ravel()
    ncomp, lab = connected_components(A, directed=False)
    big = lab == np.bincount(lab).argmax()
    Ab = A[big][:, big]; d = np.asarray(Ab.sum(1)).ravel(); Dm = sp.diags(1 / np.sqrt(d))
    N = Dm @ Ab @ Dm
    vals, vecs = eigsh(N, k=32, which='LA', tol=1e-8, maxiter=20000)
    o = np.argsort(-vals); vals, vecs = vals[o], vecs[:, o]
    emb = np.zeros((len(genes), 32), np.float32); emb[big] = vecs
    np.savez_compressed(TAB / 'partF3_string.npz', genes=genes, emb=emb, extra=np.log1p(deg).astype(np.float32), in_lcc=big, eigenvalues=vals)
    out = dict(download=[i1, i2], release='STRING v12.0 (9606)', edges_ge_700=int(len(E)), genes=int(len(genes)), components=int(ncomp),
               largest_component=int(big.sum()), eigenvalues=vals.tolist(),
               note='Spectral embedding on the largest connected component: the normalised adjacency of a graph with many '
                    'components has eigenvalue 1 once per component, so the top 32 eigenvectors of the whole graph are indicators '
                    'of small components. Genes outside it keep log(1 + degree) and get zeros plus the missing indicator.')
    save_json(TAB / 'partF3b.json', out); log(f'F3b: {len(genes)} genes, {len(E)} edges >= 700, LCC {big.sum()}; eigenvalues {vals[:5].round(3)} .. {vals[-1]:.3f}')


def partF3c(a):
    """ARCHS4: stop unless free disk >= file size + 20 GB (checked before any download)."""
    import shutil
    import subprocess
    url = 'https://s3.dev.maayanlab.cloud/archs4/files/human_gene_v2.latest.h5'
    hd = subprocess.run(['curl', '-sIL', url], capture_output=True, text=True).stdout
    size = [int(l.split(':')[1]) for l in hd.lower().splitlines() if l.startswith('content-length')][-1]
    free = shutil.disk_usage(V).free
    out = dict(url=url, file_bytes=size, free_bytes=free, needed_bytes=size + 20e9, room=bool(free >= size + 20e9),
               last_modified=[l.split(':', 1)[1].strip() for l in hd.splitlines() if l.lower().startswith('last-modified')][-1:])
    log(f'F3c disk check: file {size / 1e9:.1f} GB, free {free / 1e9:.1f} GB, need {(size + 20e9) / 1e9:.1f} GB -> {"room" if out["room"] else "NO ROOM: F3c stops"}')
    if not out['room']:
        out['status'] = 'stopped at the disk check: no download; M3c and M4c are not run'
    save_json(TAB / 'partF3c.json', out)
    if out['room']:
        raise NotImplementedError('F3c processing is not reached on vcc2027 (42 GB free for a 62 GB file)')


def feature_sets():
    """Available gene-feature sets: name -> (approved-symbol index, 32 components, extra feature)."""
    out = {}
    for nm, f in (('DepMap', 'partF3_depmap.npz'), ('STRING', 'partF3_string.npz'), ('ARCHS4', 'partF3_archs4.npz')):
        if (TAB / f).exists():
            z = np.load(TAB / f, allow_pickle=True)
            emb = z['emb'].astype(np.float64)
            ok = np.ones(len(z['genes']), bool) if nm != 'STRING' else z['in_lcc']
            out[nm] = dict(idx={g: i for i, g in enumerate(z['genes'].astype(str))}, emb=emb, extra=z['extra'].astype(np.float64), ok=ok)
    return out


def gene_features(targets, fs, H):
    """(T x 34) features for knockdown targets: 32 components, the extra feature, and a missing indicator."""
    X = np.zeros((len(targets), 34)); cos = np.zeros((len(targets), 32)); has = np.zeros(len(targets), bool)
    for i, t in enumerate(targets):
        g = to_hgnc(t, H); j = fs['idx'].get(g) if g else None
        if j is None:
            X[i, 33] = 1.0; continue
        X[i, 32] = fs['extra'][j]
        if fs['ok'][j]:
            X[i, :32] = fs['emb'][j]; cos[i] = fs['emb'][j]; has[i] = True
        else:
            X[i, 33] = 1.0
    return X, cos, has


def line_targets():
    out = {}
    for line, ev in B_LINE.items():
        out[line] = ref_for(ev).h.targets if line == 'H1' else list(np.load(SS_OUT / 'cache' / f'fold_{line}.npz', allow_pickle=True)['targets'])
    out['DEF'] = DEF_TARGETS()
    return out


def partF3sum(a):
    """Coverage, spectra, Spearman with H1's true DE count and (b), complex sanity check, cross-set similarity overlap."""
    import scorer_step12 as S12
    from scipy.stats import spearmanr
    H = hgnc_maps(); fs = feature_sets(); LT = line_targets()
    ref = ref_for('EvalA'); nt = S12.truth_calls(ref).sum(1); bH = b_for('EvalA', ref.h.targets)
    res = dict(sets=list(fs), coverage={}, spearman={}, complexes={}, overlap={})
    for nm, f in fs.items():
        res['coverage'][nm] = {l: dict(n=len(t), covered=int(gene_features(t, f, H)[2].sum())) for l, t in LT.items()}
        X, _, has = gene_features(ref.h.targets, f, H)
        x = X[:, 32]; ok = X[:, 33] == 0
        okb = ok & np.isfinite(bH)
        res['spearman'][nm] = dict(feature='mean gene effect' if nm == 'DepMap' else ('log(1 + degree)' if nm == 'STRING' else 'mean log CPM'),
                                   vs_true_DE_count=float(spearmanr(x[ok], nt[ok]).statistic), n=int(ok.sum()),
                                   vs_b=float(spearmanr(x[okb], bH[okb]).statistic), n_b=int(okb.sum()))
        # complexes: mean within-complex cosine of the 32 components against random pairs
        U = f['emb'] / np.maximum(np.linalg.norm(f['emb'], axis=1, keepdims=True), 1e-12)
        rng = np.random.default_rng(SEED); allg = np.flatnonzero(f['ok'])
        i1, i2 = rng.choice(allg, 20000), rng.choice(allg, 20000); keep = i1 != i2
        rnd = float(np.mean((U[i1[keep]] * U[i2[keep]]).sum(1)))
        cx = {}
        for c in COMPLEXES:
            mem = [f['idx'][g] for g, gr in H['groups'].items() if c in gr and g in f['idx'] and f['ok'][f['idx'][g]]]
            if len(mem) >= 3:
                S_ = U[mem] @ U[mem].T; n = len(mem)
                cx[c] = dict(members=n, mean_cos=float((S_.sum() - n) / (n * (n - 1))))
        res['complexes'][nm] = dict(random_pairs_mean_cos=rnd, complexes=cx)
        log(f'F3 {nm}: coverage {res["coverage"][nm]}; Spearman {res["spearman"][nm]}; random cos {rnd:.3f}, complexes ' +
            ', '.join(f'{k} {v["mean_cos"]:.2f}' for k, v in cx.items()))
    names = list(fs); rng = np.random.default_rng(SEED)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            A_, B_ = fs[names[i]], fs[names[j]]
            com = sorted(g for g in A_['idx'] if g in B_['idx'] and A_['ok'][A_['idx'][g]] and B_['ok'][B_['idx'][g]])
            sub = rng.choice(com, min(3000, len(com)), replace=False)
            ua = A_['emb'][[A_['idx'][g] for g in sub]]; ub = B_['emb'][[B_['idx'][g] for g in sub]]
            ua /= np.linalg.norm(ua, axis=1, keepdims=True); ub /= np.linalg.norm(ub, axis=1, keepdims=True)
            Sa, Sb = ua @ ua.T, ub @ ub.T; iu = np.triu_indices(len(sub), 1)
            res['overlap'][f'{names[i]} vs {names[j]}'] = dict(common_genes=len(com), sampled=len(sub), cosine_of_similarities=cosine(Sa[iu], Sb[iu]),
                                                               pearson_of_similarities=float(np.corrcoef(Sa[iu], Sb[iu])[0, 1]))
    save_json(TAB / 'partF3.json', res)
    print(json.dumps(res, indent=1, default=float))



# ---------------------------------------------------------------- truth lines, source scores (shared by F3d, F4, F5, F6)
def truth_of(line):
    """A truth line from its eval reference: true LFC on the real DE table's tested genes (own gene out), true Delta,
    true-DE mask, H1-style terciles of the true DE count, signal fraction max(0, D - C) / D (D = |R - control|^2,
    C = the real side's jackknife correction). The proxy Ref is dropped afterwards to save memory."""
    import scorer_step12 as S12
    ev = B_LINE[line]; ref = ref_for(ev); h = ref.h
    X = np.full((h.T, h.G), np.nan)
    rl = ref.real_df['log2_fold_change'].to_numpy().astype(float); rp = ref.real_df['p_adj'].to_numpy().astype(float)
    X[ref.ri, ref.rg] = rl
    DE = np.zeros((h.T, h.G), bool); DE[ref.ri, ref.rg] = (rp < 0.05) & np.isfinite(rl)
    for i in h.rows:
        if h.own[i] >= 0: X[i, h.own[i]] = np.nan; DE[i, h.own[i]] = False
    Dt = h.R - h.ctrl[None, :]
    Dk = (Dt ** 2).sum(1)
    nt = S12.truth_calls(ref).sum(1)
    out = dict(line=line, targets=list(h.targets), X=X, gated=np.isfinite(X), DE=DE, Delta=Dt, v=np.expm1(h.ctrl), terc=tercile_of(nt), n_true=nt,
               signal=np.maximum(0.0, Dk - h.Creal) / Dk, own=np.asarray(h.own))
    S19._REFS.pop(ev, None); gc.collect()
    return out


def bg_scores(RS, meta, targets, exclude, K, backgrounds):
    """Per knockdown x background program scores (T x B x K, NaN = background does not hold it): each background's
    tables averaged first; genetic rows only; excluded backgrounds left NaN."""
    ok = (~meta.background.isin(exclude) & meta.target.isin(targets)).to_numpy()
    d = pd.DataFrame(RS[ok, :K].astype(np.float64)); d['table'] = meta.table.to_numpy()[ok]
    d['background'] = meta.background.to_numpy()[ok]; d['target'] = meta.target.to_numpy()[ok]
    per = d.groupby(['target', 'background', 'table']).mean(numeric_only=True).groupby(level=['target', 'background']).mean()
    A = np.full((len(targets), len(backgrounds), K), np.nan)
    ti = {t: i for i, t in enumerate(targets)}; bi = {b: j for j, b in enumerate(backgrounds)}
    for (t, b), row in per.iterrows():
        if b in bi:
            A[ti[t], bi[b]] = row.to_numpy()
    return A


def avg_scores(A):
    """M1's input: the mean over the backgrounds that hold the knockdown (0 where none)."""
    n = np.isfinite(A).any(2).sum(1)
    return np.where(n[:, None] > 0, np.nansum(A, 1) / np.maximum(n, 1)[:, None], 0.0), n


# ---------------------------------------------------------------- F6: co-expression prior from the target line's own controls
F6_CTRL = {'H1': S19.EVAL['EvalA'] / 'context.h5ad', 'HCT116': S19.EVAL['HCT116final'] / 'context.h5ad',
           'K562': S19.EVAL['K562final'] / 'context.h5ad', 'HEK293T': S19.EVAL['HEK293Tfinal'] / 'context.h5ad',
           'D': V / 'data_final/extracted/context_D.h5ad', 'E': V / 'data_final/extracted/context_E.h5ad', 'F': V / 'data_final/extracted/context_F.h5ad'}
F6_GAMMA = [0.5, 1.0, 1.5]


def line_gate(line):
    """The line's gated genes (> 5 CPM), as everywhere else: H1 the 15c controls baseline, guardrails the fold cache,
    D/E/F Step 19's control stats."""
    import run_line as RL
    if line == 'H1':
        cb = pd.read_csv(C15 / 'config1p_K0/EvalA/e2e/controls/controls_baseline.csv', index_col=0)
        return cb['de_tested'].reindex(RL.panel_genes()).fillna(False).to_numpy(bool)
    if line in ('D', 'E', 'F'):
        return np.load(S19T / f'controls_{line}.npz')['gate'].astype(bool)
    return np.load(SS_OUT / 'cache' / f'fold_{line}.npz')['gate'].astype(bool)


def coexp(line):
    """C over the line's gated genes from its control cells: log1p CP10k -> top 50 PCs -> k-means (k = min(600,
    cells // 30), seed 2026) metacells; per metacell the mean over its cells of log1p CPM (gated genes) and of log library
    size; each gene regressed on the log library size; Pearson correlation of the residuals."""
    import scipy.sparse as sp
    import run_line as RL
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    genes = RL.panel_genes(); gate = line_gate(line); gi = np.flatnonzero(gate)
    path = F6_CTRL[line]
    import h5py
    with h5py.File(path, 'r') as f:
        nnz = f['X/data'].shape[0]
    mem_before(f'F6 controls {line}', nnz * 12 / 1e9 + 3)
    X = RL.read_controls(path, genes)
    n = X.shape[0]; lib = np.asarray(X.sum(1)).ravel().astype(np.float64)
    Xg = X[:, gi].astype(np.float64)
    L10 = sp.diags(1e4 / lib) @ Xg; L10.data = np.log1p(L10.data)
    mem_before(f'F6 dense {line} {n} x {len(gi)}', n * len(gi) * 8 * 2 / 1e9 + 1)
    D10 = L10.toarray()
    pcs = PCA(50, svd_solver='randomized', random_state=SEED).fit_transform(D10); del D10
    k = int(min(600, n // 30))
    lab = KMeans(n_clusters=k, n_init=1, random_state=SEED).fit_predict(pcs)
    Lc = sp.diags(1e6 / lib) @ Xg; Lc.data = np.log1p(Lc.data)
    Z = sp.csr_matrix((np.ones(n), (lab, np.arange(n))), shape=(k, n)); cnt = np.asarray(Z.sum(1)).ravel()
    Mc = np.asarray((Z @ Lc).todense()) / cnt[:, None]
    ll = (Z @ np.log(lib)) / cnt
    C, sd = resid_corr(Mc, ll)
    info = dict(line=line, cells=int(n), metacells=k, cells_per_metacell_median=float(np.median(cnt)), gated_genes=int(len(gi)),
                zero_variance_genes=int((~np.isfinite(sd)).sum()))
    log(f'F6 C {line}: {info}')
    del X, Xg, L10, Lc, Mc; gc.collect()
    return C, gi, info


def resid_corr(Mc, ll):
    """Each gene (column of Mc, metacells x genes) regressed on the metacells' mean log library size (with intercept);
    Pearson correlation between the genes' residuals (0 for a constant gene)."""
    k = len(Mc); Xd = np.column_stack([np.ones(k), ll])
    Rz = Mc - Xd @ np.linalg.lstsq(Xd, Mc, rcond=None)[0]
    Rz -= Rz.mean(0); sd = Rz.std(0); sd = np.where(sd > 1e-12, sd, np.inf)
    Zs = Rz / sd
    return (Zs.T @ Zs) / k, sd


def prior_rows(C, gi, targets, G):
    """-C_k on the panel (T x G, 0 off the gate and at the target), and whether the knockdown has a prior."""
    import run_line as RL
    gp = {g: j for j, g in enumerate(RL.panel_genes())}; pos = {j: i for i, j in enumerate(gi)}
    out = np.zeros((len(targets), G)); has = np.zeros(len(targets), bool)
    for i, t in enumerate(targets):
        j = gp.get(t)
        if j is None or j not in pos:
            continue
        row = -C[pos[j]].astype(np.float64); row[pos[j]] = 0.0
        out[i, gi] = row; has[i] = True
    return out, has


def f6_diag(tr, C, gi, n_null=100):
    """Per knockdown Pearson of its true Delta (gated, own out) with -C_k; nulls from random gated genes' rows; share of
    its true DE genes among the target's top-100 co-expressed genes (largest C)."""
    rng = np.random.default_rng(SEED); G = tr['Delta'].shape[1]
    Pk, has = prior_rows(C, gi, tr['targets'], G)
    pos_of = {j: i for i, j in enumerate(gi)}
    r_obs = np.full(len(has), np.nan); r_null = np.full((len(has), n_null), np.nan)
    de_obs = np.full(len(has), np.nan); de_null = np.full((len(has), n_null), np.nan)
    for i in np.flatnonzero(has):
        own = tr['own'][i]; m = np.zeros(G, bool); m[gi] = True
        if own >= 0: m[own] = False
        y = tr['Delta'][i, m]; r_obs[i] = np.corrcoef(y, Pk[i, m])[0, 1]
        de = tr['DE'][i, m]; mg = np.flatnonzero(m)
        top = mg[np.argsort(Pk[i, m])[:100]]                                   # largest C = most negative -C
        if de.sum():
            de_obs[i] = float(np.isin(np.flatnonzero(tr['DE'][i] & m), top).sum() / de.sum())
        draws = rng.choice([j for j in gi if j != own], n_null, replace=False)
        for d, j in enumerate(draws):
            row = -C[pos_of[j]].astype(np.float64); full = np.zeros(G); full[gi] = row
            if own >= 0: full[own] = 0.0
            r_null[i, d] = np.corrcoef(y, full[m])[0, 1]
            if de.sum():
                topn = mg[np.argsort(full[m])[:100]]
                de_null[i, d] = float(np.isin(np.flatnonzero(tr['DE'][i] & m), topn).sum() / de.sum())
    out = dict(knockdowns_with_prior=int(has.sum()), knockdowns=int(len(has)))
    for nm, rr in (('all', has),) + tuple((t, has & (tr['terc'] == j)) for j, t in enumerate(TERC)):
        null_meds = np.nanmedian(r_null[rr], 0)
        out[nm] = dict(n=int(rr.sum()), median_r=float(np.nanmedian(r_obs[rr])), null_median_r=float(np.nanmedian(null_meds)),
                       null_median_r_q95=float(np.nanquantile(null_meds, 0.95)),
                       percentile_vs_null=float((null_meds < np.nanmedian(r_obs[rr])).mean() * 100),
                       de_in_top100=float(np.nanmean(de_obs[rr])), de_in_top100_null=float(np.nanmean(de_null[rr])))
    return out, Pk, has


def partF6(a):
    """F6a diagnostics on the four truth lines; F6b proxy arm (lambda per (b) tercile on the fit half; gamma grid);
    C for D, E and F saved (float16) with -C_k for the 300 targets."""
    import scorer_step12 as S12
    import scorer_step16 as S16
    res = dict(info={}, diag={}, fits={}, guardrails={}, lambda_={})
    tgt = DEF_TARGETS()
    for x in ('D', 'E', 'F'):
        C, gi, info = coexp(x); res['info'][x] = info
        Pk, has = prior_rows(C, gi, tgt, len(line_gate(x)))
        np.savez_compressed(TAB / f'partF6_C_{x}.npz', C=C.astype(np.float16), gated_index=gi, targets=np.array(tgt),
                            minus_C_k=Pk.astype(np.float32), has_prior=has)
        log(f'F6 {x}: {has.sum()} of {len(tgt)} targets have a prior'); del C; gc.collect()
    save_json(TAB / 'partF6.json', res)
    Pks = {}
    for line in D_LINES:
        C, gi, info = coexp(line); res['info'][line] = info
        tr = truth_of(line)
        d, Pk, has = f6_diag(tr, C, gi); res['diag'][line] = d
        Pks[line] = (Pk, has, tr)
        np.savez_compressed(TAB / f'partF6_prior_{line}.npz', minus_C_k=Pk.astype(np.float32), has_prior=has, targets=np.array(tr['targets']))
        log(f'F6a {line}: {d["all"]}')
        del C; gc.collect(); save_json(TAB / 'partF6.json', res)
    # F6b: lambda_t through the origin on the fit half's knockdowns of (b) tercile t (gated, own out)
    Pk, has, tr = Pks['H1']
    D = DLine('H1'); h = D.ref.h; mu0 = base_mu(D)
    T_ = S12.truth_calls(D.ref); halves = S16.halves(T_); bH = b_for('EvalA', h.targets)
    ok = np.broadcast_to(line_gate('H1')[None, :], Pk.shape).copy()
    for i in h.rows:
        if h.own[i] >= 0: ok[i, h.own[i]] = False

    def lam(rows, cuts):
        t = b_terc(bH, cuts); out = np.zeros(3)
        for j in range(3):
            rr = np.intersect1d(rows, np.flatnonzero((t == j) & has))
            x = np.where(ok[rr], Pk[rr], 0.0); y = np.where(ok[rr], tr['Delta'][rr], 0.0)
            out[j] = (x * y).sum() / max((x * x).sum(), 1e-30)
        return out
    rows = []
    for f in (0, 1):
        cuts = b_cuts(bH[halves[f]]); L = lam(halves[f], cuts); res['lambda_'][f'fit_half{f}'] = dict(cuts=cuts.tolist(), lam=L.tolist())
        t = b_terc(bH, cuts)
        rr, _ = score_rows(D, mu0, gamma=0.0, fit_half=f); rows += rr
        for g in F6_GAMMA:
            add = g * L[t][:, None] * np.where(ok & has[:, None], Pk, 0.0)
            rr, _ = score_rows(D, mu_add_delta(mu0, add, own_index(h)), gamma=g, fit_half=f); rows += rr
        log(f'F6b fit half {f}: lambda {L.round(4)}')
    cuts_all = b_cuts(bH); L_all = lam(np.arange(h.T), cuts_all); res['lambda_']['all_H1'] = dict(cuts=cuts_all.tolist(), lam=L_all.tolist())
    del D; gc.collect()
    for line in D_LINES[1:]:
        Pk, has, tr = Pks[line]
        D = DLine(line); mu0 = base_mu(D); b = b_for(B_LINE[line], D.ref.h.targets); t = b_terc(b, cuts_all)
        okl = np.broadcast_to(line_gate(line)[None, :], Pk.shape).copy()
        for i in D.ref.h.rows:
            if D.ref.h.own[i] >= 0: okl[i, D.ref.h.own[i]] = False
        rr, _ = score_rows(D, mu0, gamma=0.0, fit_half=-1); rows += rr
        for g in F6_GAMMA:
            add = g * L_all[t][:, None] * np.where(okl & has[:, None], Pk, 0.0)
            rr, _ = score_rows(D, mu_add_delta(mu0, add, own_index(D.ref.h)), gamma=g, fit_half=-1); rows += rr
        del D; gc.collect()
    R = pd.DataFrame(rows); R.to_csv(TAB / 'partF6_rows.csv', index=False)
    H = R[R.line == 'H1']
    for f in (0, 1):
        got = pick_fit(H[H.fit_half == f], f'half{f}', f'half{1 - f}', ['gamma'], lambda d: d.gamma == 0)
        if got is None:
            res['fits'][f'fit_half{f}'] = dict(picked=None); continue
        pk, ho, b0, ne = got
        g = {}
        for line in D_LINES[1:]:
            Lr = R[(R.line == line)]
            x, y = Lr[Lr.gamma == pk.gamma].iloc[0], Lr[Lr.gamma == 0].iloc[0]
            g[line] = dict(pick=comp_dict(x), base=comp_dict(y), diff=float(x.overall - y.overall))
        res['fits'][f'fit_half{f}'] = dict(picked=float(pk.gamma), eligible=ne, grid_position=grid_pos([0.0] + F6_GAMMA, pk.gamma),
                                          held_out=comp_dict(ho), base_held_out=comp_dict(b0),
                                          gain={c: float(ho[c] - b0[c]) for c in COMPS + ['overall']}, guardrails=g)
    fv = [v for v in res['fits'].values() if v.get('picked') is not None]
    res['criterion'] = dict(both_picked=len(fv) == 2, heldout_gain_ge_0005_both=bool(len(fv) == 2 and all(v['gain']['overall'] >= 0.005 for v in fv)),
                            all_components_pos_both=bool(len(fv) == 2 and all(all(v['held_out'][c] > 0 for c in COMPS) for v in fv)),
                            guardrails_within_0005=bool(len(fv) == 2 and all(x['diff'] >= -0.005 for v in fv for x in v['guardrails'].values())))
    res['criterion']['met'] = all(res['criterion'].values())
    save_json(TAB / 'partF6.json', res)
    print(json.dumps({k: v for k, v in res.items()}, indent=1, default=float))


# ---------------------------------------------------------------- F3d: do the embeddings carry response information?
def pair_cos_both_out(Dg, own):
    """Cosine between every two rows of Dg (gated genes already masked), each pair with both rows' own genes out."""
    Dg = Dg.copy(); T = len(Dg)
    for i in range(T):
        if own[i] >= 0: Dg[i, own[i]] = 0.0
    P = Dg @ Dg.T; nrm = (Dg ** 2).sum(1)
    Q = np.where(np.asarray(own)[None, :] >= 0, Dg[:, np.maximum(own, 0)], 0.0)             # Q[i, j] = Dg_i at own_j
    ni = nrm[:, None] - Q ** 2; nj = nrm[None, :] - (Q ** 2).T
    return P / np.sqrt(np.maximum(ni, 1e-30) * np.maximum(nj, 1e-30))


def partF3d(a, n_perm=100):
    from scipy.stats import spearmanr
    H = hgnc_maps(); fs = feature_sets(); res = {}
    meta = pd.read_parquet(TAB / 'partF1_srcrows_delta.parquet'); bgs = sorted(meta.background.unique())
    rng = np.random.default_rng(SEED)
    for line in D_LINES:
        tr = truth_of(line); T = len(tr['targets']); gate = line_gate(line)
        Dg = np.where(gate[None, :], tr['Delta'], 0.0)
        own = tr['own']
        # pair cosine of true Delta over gated genes with both targets out
        N = pair_cos_both_out(Dg, own)
        iu = np.triu_indices(T, 1); dcos = N[iu]
        ms = tr['terc'] >= 1; msp = ms[iu[0]] & ms[iu[1]]
        out = {}
        RS = np.load(TAB / f'partF1_RS_E4w_{FOLD_OF[line]}.npy')
        A = bg_scores(RS, meta, tr['targets'], OWN_BG_F[line], 100, bgs); sh, nb = avg_scores(A)
        okb = nb > 0; U = sh / np.maximum(np.linalg.norm(sh, axis=1, keepdims=True), 1e-12); scos = (U @ U.T)[iu]
        pb = okb[iu[0]] & okb[iu[1]]
        out['source_scores_E4w_K100'] = dict(pairs=int(pb.sum()), spearman=float(spearmanr(scos[pb], dcos[pb]).statistic),
                                             pairs_mid_strong=int((pb & msp).sum()), spearman_mid_strong=float(spearmanr(scos[pb & msp], dcos[pb & msp]).statistic))
        for nm, f in fs.items():
            _, E, has = gene_features(tr['targets'], f, H)
            U = E / np.maximum(np.linalg.norm(E, axis=1, keepdims=True), 1e-12)
            pe = has[iu[0]] & has[iu[1]]
            ecos = (U @ U.T)[iu]
            obs = spearmanr(ecos[pe], dcos[pe]).statistic; obs_ms = spearmanr(ecos[pe & msp], dcos[pe & msp]).statistic
            hi = np.flatnonzero(has); nul, nul_ms = [], []
            for _ in range(n_perm):
                Up = U.copy(); Up[hi] = U[rng.permutation(hi)]
                pc = (Up @ Up.T)[iu]
                nul.append(spearmanr(pc[pe], dcos[pe]).statistic); nul_ms.append(spearmanr(pc[pe & msp], dcos[pe & msp]).statistic)
            nul, nul_ms = np.array(nul), np.array(nul_ms)
            out[nm] = dict(knockdowns_with_features=int(has.sum()), pairs=int(pe.sum()), spearman=float(obs), null_mean=float(nul.mean()), null_sd=float(nul.std()),
                           percentile=float((nul < obs).mean() * 100), pairs_mid_strong=int((pe & msp).sum()), spearman_mid_strong=float(obs_ms),
                           null_mean_mid_strong=float(nul_ms.mean()), percentile_mid_strong=float((nul_ms < obs_ms).mean() * 100))
        res[line] = out; log(f'F3d {line}: {out}')
        save_json(TAB / 'partF3d.json', res)
    print(json.dumps(res, indent=1, default=float))



# ---------------------------------------------------------------- F4: coefficient prediction
ALPHAS = [0.1, 1.0, 10.0, 100.0, 1000.0]
F4_BASES = ['E2w', 'E4w', 'E24w', 'E2', 'E4']
SPACE = {'E2w': 'lfc', 'E4w': 'delta', 'E24w': 'delta', 'E2': 'lfc', 'E4': 'delta'}
POOL_TRAIN = {'H1': ['HCT116', 'K562', 'HEK293T'], 'K562': ['H1', 'HCT116', 'HEK293T'], 'HCT116': ['H1', 'K562'], 'HEK293T': ['H1', 'K562']}
N_NULL = 20


def program_truth(tr, Vt, Ks, space):
    """Least-squares program scores of a truth line's true values (LFC or Delta) over its gated genes, for nested K
    (normal equations on the top max(K) rows, sub-blocks per K)."""
    X = tr['X'] if space == 'lfc' else tr['Delta']
    Km = max(Ks); B = Vt[:Km].astype(np.float64); out = {K: np.zeros((len(X), K)) for K in Ks}
    for i in range(len(X)):
        g = tr['gated'][i]; Bg = B[:, g]; Gm = Bg @ Bg.T; r = Bg @ X[i, g]
        for K in Ks:
            out[K][i] = np.linalg.solve(Gm[:K, :K], r[:K])
    return out


def std_design(Xtr, w):
    """Weighted standardisation per program and column: (mean, sd) with sd = inf for constant columns."""
    ws = w / w.sum()
    mu = np.einsum('n,jnp->jp', ws, Xtr)
    sd = np.sqrt(np.einsum('n,jnp->jp', ws, (Xtr - mu[:, None, :]) ** 2))
    return mu, np.where(sd > 1e-12, sd, np.inf)


def ridge_solve(Xtr, Ytr, w, alpha):
    """Weighted ridge with intercept on standardised columns, per program. Xtr (J, n, p), Ytr (n, J), alpha (J,)."""
    mu, sd = std_design(Xtr, w); Z = (Xtr - mu[:, None, :]) / sd[:, None, :]
    ws = w / w.sum(); ym = ws @ Ytr; Yc = Ytr - ym[None, :]
    A = np.einsum('jnp,n,jnq->jpq', Z, w, Z) + np.asarray(alpha)[:, None, None] * np.eye(Z.shape[2])[None]
    b = np.einsum('jnp,n,nj->jp', Z, w, Yc)
    beta = np.linalg.solve(A, b[..., None])[..., 0]
    return dict(mu=mu, sd=sd, beta=beta, ym=ym)


def ridge_predict(m, X):
    Z = (X - m['mu'][:, None, :]) / m['sd'][:, None, :]
    return (np.einsum('jnp,jp->nj', Z, m['beta']) + m['ym'][None, :])


def ridge_cv(Xtr, Ytr, w, alphas=ALPHAS, n_fold=5, seed=SEED):
    """alpha per program by 5-fold CV (weighted SSE) within the training set; refit on all of it. cv_r2 = 1 - SSE / sum w y^2
    at the chosen alpha (relative to predicting 0, the convention of every R^2 here)."""
    n, J = Ytr.shape
    fold = np.random.default_rng(seed).permutation(n) % n_fold
    sse = np.zeros((len(alphas), J)); tot = np.zeros(J)
    for f in range(n_fold):
        tr, va = fold != f, fold == f
        for ai, al in enumerate(alphas):
            m = ridge_solve(Xtr[:, tr], Ytr[tr], w[tr], np.full(J, al))
            e = Ytr[va] - ridge_predict(m, Xtr[:, va])
            sse[ai] += (w[va][:, None] * e ** 2).sum(0)
        tot += (w[va][:, None] * Ytr[va] ** 2).sum(0)
    best = sse.argmin(0); alpha = np.array(alphas)[best]
    m = ridge_solve(Xtr, Ytr, w, alpha); m.update(alpha=alpha, cv_r2=1 - sse[best, np.arange(J)] / np.maximum(tot, 1e-30))
    return m


def build_X(method, A, GF):
    """Per-program design (J, n, p). A: source scores (n, B, K) NaN = missing; GF: (n, 34) gene features or None."""
    n, B, K = A.shape
    parts = []
    if method != 'M4':
        S_ = np.nan_to_num(np.transpose(A, (2, 0, 1)))                              # (K, n, B)
        miss = np.broadcast_to(np.isnan(A[:, :, 0]).astype(float)[None], (K, n, B))
        parts += [S_, miss]
    if GF is not None:
        parts.append(np.broadcast_to(GF[None], (K, n, GF.shape[1])))
    return np.concatenate(parts, axis=2)


class FCtx:
    """Shared Part F state: truth lines, gene features, source rows meta, bases."""

    def __init__(self, lines=D_LINES):
        import scorer_step16 as S16
        import scorer_step12 as S12
        self.H = hgnc_maps(); self.fs = {k: v for k, v in feature_sets().items() if k in ('DepMap', 'STRING', 'ARCHS4')}
        self.tr = {l: truth_of(l) for l in lines}
        ref = ref_for('EvalA'); self.halves = S16.halves(S12.truth_calls(ref)); S19._REFS.pop('EvalA', None); gc.collect()
        self.meta = {sp: pd.read_parquet(TAB / f'partF1_srcrows_{sp}.parquet') for sp in ('lfc', 'delta')}
        self.bgs = {sp: sorted(m.background.unique()) for sp, m in self.meta.items()}
        self.GF = {(l, nm): gene_features(self.tr[l]['targets'], f, self.H)[0] for l in lines for nm, f in self.fs.items()}
        self._rs, self._vt, self._lab, self._sc = {}, {}, {}, {}

    def vt(self, b, fold):
        if (b, fold) not in self._vt:
            self._vt[(b, fold)] = np.load(basis_file(b, fold))['Vt'].astype(np.float64)
        return self._vt[(b, fold)]

    def rs(self, b, fold):
        if (b, fold) not in self._rs:
            self._rs[(b, fold)] = np.load(TAB / f'partF1_RS_{b}_{fold}.npy')
        return self._rs[(b, fold)]

    def labels(self, line, b, fold, Ks):
        key = (line, b, fold, tuple(Ks))
        if key not in self._lab:
            self._lab[key] = program_truth(self.tr[line], self.vt(b, fold), Ks, SPACE[b])
        return self._lab[key]

    def scores(self, line, b, fold, K, exclude):
        key = (line, b, fold, tuple(sorted(exclude)))
        if key not in self._sc:
            sp = SPACE[b]; RS = self.rs(b, fold)
            self._sc[key] = bg_scores(RS, self.meta[sp], self.tr[line]['targets'], exclude, RS.shape[1], self.bgs[sp])
        return self._sc[key][:, :, :K]

    def weights(self, line):
        s = self.tr[line]['signal']; return s / s.sum()


def fit_predict(ctx, method, b, fold, K, Ks, train, pred_A, pred_GF, fset=None, gf_override=None):
    """Fit one method on training blocks and predict program scores for pred rows.
    train: list of (line, rows, exclude, weighted); gf_override: line -> gene features to use in place of the line's own
    (the M4 permutation null). Returns (pred (m, K), model info)."""
    As, Ys, Gs, Ws = [], [], [], []
    for line, rows, excl, weighted in train:
        As.append(ctx.scores(line, b, fold, K, excl)[rows]); Ys.append(ctx.labels(line, b, fold, Ks)[K][rows])
        Gs.append(((gf_override or {}).get(line, ctx.GF[(line, fset)] if fset else None))[rows] if fset else None)
        w = ctx.weights(line)[rows] if weighted else np.ones(len(rows))
        Ws.append(w / w.sum() if weighted else w)
    A = np.concatenate(As); Y = np.concatenate(Ys); W = np.concatenate(Ws)
    W = W * len(W) / W.sum()                       # mean weight 1, so alpha means the same with or without weights
    G = np.concatenate(Gs) if fset else None
    if method == 'M1':
        sh, _ = avg_scores(A); beta = S19.fit_beta(Y, sh, np.arange(len(Y)))
        ph, _ = avg_scores(pred_A)
        return ph * beta[None, :], dict(beta=beta)
    X = build_X(method, A, G); m = ridge_cv(X, Y, W)
    return ridge_predict(m, build_X(method, pred_A, pred_GF)), dict(alpha=m['alpha'], cv_r2=m['cv_r2'])


def energy_delta(tr, P):
    T = len(P); out = np.zeros((T, 2))
    for i in range(T):
        g = tr['gated'][i]; d = tr['Delta'][i, g]; out[i] = [((d - P[i, g]) ** 2).sum(), (d ** 2).sum()]
    return out


def partF4(a):
    """Held-out (H1 halves, cross-fitted) R^2 and captured energy of every method, and the program-score predictions F5
    needs (H1 per fit direction; guardrails from all of H1 or the pooled lines, under the fold's basis)."""
    import pickle
    ctx = FCtx(); trH = ctx.tr['H1']; hv = ctx.halves; T = len(trH['targets'])
    C2 = pd.read_csv(TAB / 'partF2_ceiling.csv')
    fsets = list(ctx.fs)
    rows, preds, choice = [], {}, {}
    for b in F4_BASES:
        k25 = k25_of(b); Ks = sorted(set([5, 10, 20, min(k25 or 100, 100)]))
        Vt = ctx.vt(b, 'full'); lab = ctx.labels('H1', b, 'full', Ks)
        for K in Ks:
            AH = ctx.scores('H1', b, 'full', K, ('H1',))
            meth = {'M1': (None, False), 'M2': (None, False)}
            meth.update({f'M3{nm}': (nm, False) for nm in fsets}); meth.update({f'M4{nm}': (nm, False) for nm in fsets})
            meth.update({'M2p': (None, True)}); meth.update({f'M3p{nm}': (nm, True) for nm in fsets})
            cf = {}
            for mname, (fset, pooled) in meth.items():
                base_m = mname[:2] if not mname.startswith(('M2p', 'M3p')) else mname[:2]
                cf[mname] = np.zeros((T, K)); info = {}
                for f in (0, 1):
                    train = [('H1', hv[f], ('H1',), pooled)]
                    if pooled:
                        train += [(l, np.arange(len(ctx.tr[l]['targets'])), OWN_BG_F[l], True) for l in POOL_TRAIN['H1']]
                    pred, inf = fit_predict(ctx, base_m, b, 'full', K, Ks, train, AH, ctx.GF[('H1', fset)] if fset else None, fset)
                    cf[mname][hv[1 - f]] = pred[hv[1 - f]]
                    info[f] = inf
                    if b in F_BASES and K in (10, 20, min(k25 or 100, 100)):
                        preds[(b, K, mname, 'H1', f)] = pred.astype(np.float32)
                y = lab[K]
                r2 = 1 - ((y - cf[mname]) ** 2).sum(0) / (y ** 2).sum(0)
                r2_half = {f: 1 - ((y[hv[1 - f]] - cf[mname][hv[1 - f]]) ** 2).sum(0) / (y[hv[1 - f]] ** 2).sum(0) for f in (0, 1)}
                P = cf[mname] @ Vt[:K]
                if SPACE[b] == 'lfc':
                    es = S19.energy_sums(trH, P); cap = {f'{m}_{g}': v for m in 'abc' for g, v in S19.captured(es[m], trH['terc']).items()}
                else:
                    cap = {f'c_{g}': v for g, v in S19.captured(energy_delta(trH, P), trH['terc']).items()}
                cl = C2[(C2.basis == b) & (C2.K == K)] if b in F_BASES else C2[(C2.basis == b) & (C2.K == K) & (C2.step19 == True)]
                ceil = {c: float(cl[c].iloc[0]) for c in cap if len(cl) and c in cl and np.isfinite(cl[c].iloc[0])}
                row = dict(basis=b, K=K, K_is_25pct=K == min(k25 or 100, 100), method=mname, median_r2_heldout=float(np.median(r2)),
                           median_r2_heldout_half0=float(np.median(r2_half[1])), median_r2_heldout_half1=float(np.median(r2_half[0])),
                           r2_heldout_first5=r2[:5].round(3).tolist(), **cap,
                           **{f'share_of_ceiling_{c}': (cap[c] / ceil[c] if ceil.get(c) else np.nan) for c in cap if c.endswith(('all', 'strong'))})
                for f in (0, 1):
                    if 'cv_r2' in info[f]:
                        row[f'cv_r2_median_fit{f}'] = float(np.median(info[f]['cv_r2']))
                        row[f'alpha_counts_fit{f}'] = {str(al): int((info[f]['alpha'] == al).sum()) for al in ALPHAS}
                rows.append(row)
            # the M4 permutation null (gene-to-embedding mapping shuffled; 20 permutations per set)
            rng = np.random.default_rng(SEED)
            for nm in fsets:
                caps = []
                for _ in range(N_NULL):
                    perm = rng.permutation(T); cfn = np.zeros((T, K))
                    for f in (0, 1):
                        Gp = ctx.GF[('H1', nm)][perm]
                        pred, _ = fit_predict(ctx, 'M4', b, 'full', K, Ks, [('H1', hv[f], ('H1',), False)], AH, Gp, nm, gf_override={'H1': Gp})
                        cfn[hv[1 - f]] = pred[hv[1 - f]]
                    P = cfn @ Vt[:K]
                    caps.append(S19.captured(energy_delta(trH, P) if SPACE[b] == 'delta' else S19.energy_sums(trH, P)['c'], trH['terc']))
                    caps[-1]['r2'] = float(np.median(1 - ((lab[K] - cfn) ** 2).sum(0) / (lab[K] ** 2).sum(0)))
                rows.append(dict(basis=b, K=K, method=f'M4{nm}_null', null_draws=N_NULL, median_r2_heldout=float(np.mean([c['r2'] for c in caps])),
                                 median_r2_heldout_null_q95=float(np.quantile([c['r2'] for c in caps], 0.95)),
                                 **{f'c_{g}': float(np.mean([c[g] for c in caps])) for g in ('all',) + TERC},
                                 c_strong_null_q95=float(np.quantile([c['strong'] for c in caps], 0.95))))
            log(f'F4 {b} K {K}: ' + ', '.join(f'{r["method"]} R2 {r["median_r2_heldout"]:.3f} c_strong {r.get("c_strong", np.nan):.3f}' for r in rows if r['basis'] == b and r['K'] == K))
            pd.DataFrame(rows).to_csv(TAB / 'partF4_rows.csv', index=False)
        # guardrail predictions (F_BASES only): fit on all of H1 (or the pooled lines) under the fold's basis
        if b not in F_BASES:
            continue
        for line in D_LINES[1:]:
            fold = FOLD_OF[line]; k25f = k25_of(b, fold); Kf = sorted(set([10, 20, min(k25f or 100, 100)]))
            Ksf = sorted(set(Kf + [5]))
            for K in Kf:
                Ag = ctx.scores(line, b, fold, K, OWN_BG_F[line])
                excl = OWN_BG_F[line]
                for mname, fset, pooled in [('M1', None, False), ('M2', None, False)] + [(f'M3{nm}', nm, False) for nm in fsets] + \
                                            [('M2p', None, True)] + [(f'M3p{nm}', nm, True) for nm in fsets]:
                    train = [('H1', np.arange(T), tuple(set(('H1',) + excl)), pooled)]
                    if pooled:
                        train += [(l, np.arange(len(ctx.tr[l]['targets'])), tuple(set(OWN_BG_F[l] + excl)), True) for l in POOL_TRAIN[line]]
                    pred, inf = fit_predict(ctx, mname[:2], b, fold, K, Ksf, train, Ag, ctx.GF[(line, fset)] if fset else None, fset)
                    preds[(b, K, mname, line, -1)] = pred.astype(np.float32)
                    if 'cv_r2' in inf:
                        choice[(b, K, mname, line, -1)] = float(np.median(inf['cv_r2']))
            log(f'F4 {b}: guardrail predictions for {line} done')
        # M3 / M3p variant choice per fit direction: highest median inner-CV R^2 on the training set
    R = pd.DataFrame(rows); R.to_csv(TAB / 'partF4_rows.csv', index=False)
    sel = {}
    for b in F_BASES:
        for K in sorted({k[1] for k in preds if k[0] == b}):
            for f in (0, 1):
                for fam in ('M3', 'M3p'):
                    c = {nm: R[(R.basis == b) & (R.K == K) & (R.method == f'{fam}{nm}')][f'cv_r2_median_fit{f}'] for nm in fsets}
                    c = {nm: float(v.iloc[0]) for nm, v in c.items() if len(v)}
                    if c: sel[(b, K, fam, 'H1', f)] = max(c, key=c.get)
            for line in D_LINES[1:]:
                for fam in ('M3', 'M3p'):
                    c = {nm: choice.get((b, K, f'{fam}{nm}', line, -1)) for nm in fsets}
                    c = {k: v for k, v in c.items() if v is not None}
                    if c: sel[(b, K, fam, line, -1)] = max(c, key=c.get)
    with open(TAB / 'partF4_preds.pkl', 'wb') as fh:
        pickle.dump(dict(preds=preds, selection=sel, fsets=fsets), fh)
    save_json(TAB / 'partF4.json', dict(selection={str(k): v for k, v in sel.items()}, fsets=fsets,
                                        weights={l: dict(q=np.quantile(ctx.tr[l]['signal'], [0, .1, .25, .5, .75, .9, 1]).tolist(),
                                                         above_0p1=int((ctx.tr[l]['signal'] > 0.1).sum()), n=len(ctx.tr[l]['signal'])) for l in D_LINES}))
    with pd.option_context('display.width', 300, 'display.max_rows', 500, 'display.max_columns', 40):
        print(R.drop(columns=[c for c in R.columns if c.startswith('alpha_counts') or c == 'r2_heldout_first5']).round(3).to_string())


# ---------------------------------------------------------------- F5: proxy arms
F5_METHODS = ['M1', 'M2', 'M3', 'M2p', 'M3p']
F5_ARMS = [(arm, g) for arm in ('R2', 'R3') for g in (0.25, 0.5, 0.75, 1.0)]


def f5_P(ctx_vt, b, fold, K, shat, gate):
    return (shat.astype(np.float64) @ ctx_vt[:K]) * gate[None, :]


def f5_mu(D, mu0, lfc0, b, P, arm, g, bvals):
    import scorer_arms as SA
    on = np.ones(len(P), bool) if arm == 'R2' else (np.nan_to_num(bvals, nan=-np.inf) >= 5)
    add = np.where(on[:, None], g * P, 0.0)
    if SPACE[b] == 'lfc':
        return SA.mu_from_lfc(lfc0 + add, D.m0)
    return mu_add_delta(mu0, add, own_index(D.ref.h))


def partF5(a):
    """Worker 0: H1 (both fit directions); worker 1: K562, HCT116, HEK293T. Resumable."""
    import pickle
    z = pickle.load(open(TAB / 'partF4_preds.pkl', 'rb')); preds, sel = z['preds'], z['selection']
    lines = ['H1'] if a.worker == 0 else D_LINES[1:]
    f_rows = TAB / f'partF5_rows_w{a.worker}.csv'
    rows = pd.read_csv(f_rows).to_dict('records') if f_rows.exists() else []
    done = {(r['line'], r['fit_half'], r['basis'], r['method'], r['K'], r['arm'], r['gamma']) for r in rows}
    for line in lines:
        D = DLine(line); mu0 = base_mu(D); lfc0 = np.log2((mu0 + 1e-9) / (D.m0[None, :] + 1e-9))
        import scorer_arms as SA
        chk = float(np.abs(SA.mu_from_lfc(lfc0, D.m0) - mu0).max() / max(mu0.max(), 1))
        log(f'F5 {line}: LFC round trip of the base max rel diff {chk:.1e}')
        bvals = b_for(B_LINE[line], D.ref.h.targets)
        fold = FOLD_OF[line]
        if not any(r['line'] == line and r['method'] == 'base' for r in rows):
            rr, _ = score_rows(D, mu0, fit_half=-1, basis='base', method='base', K=0, arm='base', gamma=0.0, variant=''); rows += rr
        dirs = (0, 1) if line == 'H1' else (-1,)
        t0 = time.time(); n = 0
        for f in dirs:
            for b in F_BASES:
                Vt = np.load(basis_file(b, fold))['Vt'].astype(np.float64)
                for K in sorted({k[1] for k in preds if k[0] == b and k[3] == line}):
                    for m in F5_METHODS:
                        key_m = m
                        if m in ('M3', 'M3p'):
                            v = sel.get((b, K, m, line, f)); key_m = f'{m}{v}'
                        P = f5_P(Vt, b, fold, K, preds[(b, K, key_m, line, f)], D.gate)
                        for arm, g in F5_ARMS:
                            if (line, f, b, m, K, arm, g) in done:
                                continue
                            rr, _ = score_rows(D, f5_mu(D, mu0, lfc0, b, P, arm, g, bvals), fit_half=f, basis=b, method=m, K=K, arm=arm, gamma=g, variant=key_m)
                            rows += rr; n += 1
                            if n % 10 == 0:
                                pd.DataFrame(rows).to_csv(f_rows, index=False)
                                log(f'F5 {line}: {n} rows, {(time.time() - t0) / n:.1f} s per row')
        pd.DataFrame(rows).to_csv(f_rows, index=False)
        log(f'F5 {line}: done'); del D; gc.collect()


def partF5time(a):
    """One proxy row at the thread setting the workers use, before splitting."""
    D = DLine('H1'); mu0 = base_mu(D)
    for _ in range(2):
        t0 = time.time(); score_rows(D, mu0); log(f'F5 timing: one H1 proxy row {time.time() - t0:.1f} s at OMP_NUM_THREADS={os.environ.get("OMP_NUM_THREADS")}')



def scores_targets(ctx, targets, b, exclude, with_arc):
    """Source scores (n x B' x 100) for arbitrary targets; with_arc adds Arc H1 as one more background column 'H1'."""
    sp = SPACE[b]; RS = ctx.rs(b, 'full')
    A = bg_scores(RS, ctx.meta[sp], targets, exclude, RS.shape[1], ctx.bgs[sp])
    if not with_arc:
        return A
    z = np.load(TAB / f'partF1_arc_{b}.npz', allow_pickle=True); ai = {t: i for i, t in enumerate(z['targets'].astype(str))}
    col = np.full((len(targets), 1, RS.shape[1]), np.nan)
    if 'H1' not in exclude:
        for i, t in enumerate(targets):
            if t in ai: col[i, 0] = z['scores'][ai[t]]
    return np.concatenate([A, col], axis=1)


def def_pprog(ctx, b, m, v, K):
    """D/E/F P_prog for one picked (basis, method, K): all sources with Arc H1 as one more background; fit on all of H1
    (pooled methods: all of H1 plus the three guardrail lines, Arc H1's column available to the guardrail rows)."""
    tg = DEF_TARGETS(); fset = v[len(m):] if m in ('M3', 'M3p') else None
    pooled = m in ('M2p', 'M3p'); base_m = m[:2]
    train = [('H1', np.arange(len(ctx.tr['H1']['targets'])), ('H1',))] + ([(l, np.arange(len(ctx.tr[l]['targets'])), OWN_BG_F[l]) for l in D_LINES[1:]] if pooled else [])
    Ks = sorted(set([5, 10, 20, min(k25_of(b) or 100, 100)]))
    As, Ys, Gs, Ws = [], [], [], []
    for line, rows, excl in train:
        As.append(scores_targets(ctx, ctx.tr[line]['targets'], b, excl, True)[rows, :, :K]); Ys.append(ctx.labels(line, b, 'full', Ks)[K][rows])
        Gs.append(ctx.GF[(line, fset)][rows] if fset else None)
        w = ctx.weights(line)[rows] if pooled else np.ones(len(rows)); Ws.append(w / w.sum() if pooled else w)
    A = np.concatenate(As); Y = np.concatenate(Ys); W = np.concatenate(Ws); W = W * len(W) / W.sum(); G = np.concatenate(Gs) if fset else None
    Ad = scores_targets(ctx, tg, b, (), True)[:, :, :K]
    Gd = gene_features(tg, ctx.fs[fset], ctx.H)[0] if fset else None
    if base_m == 'M1':
        sh, _ = avg_scores(A); beta = S19.fit_beta(Y, sh, np.arange(len(Y))); pred = avg_scores(Ad)[0] * beta[None, :]
    else:
        mod = ridge_cv(build_X(base_m, A, G), Y, W); pred = ridge_predict(mod, build_X(base_m, Ad, Gd))
    Vt = ctx.vt(b, 'full')
    out = dict(targets=np.array(tg), shat=pred.astype(np.float32), n_backgrounds=np.isfinite(Ad[:, :, 0]).sum(1))
    for x in ('D', 'E', 'F'):
        out[f'P_{x}'] = ((pred @ Vt[:K]) * line_gate(x)[None, :]).astype(np.float32)
    return out


def partF5sum(a):
    import pickle
    R = pd.concat([pd.read_csv(TAB / f'partF5_rows_w{w}.csv') for w in (0, 1)], ignore_index=True)
    R.to_csv(TAB / 'partF5_rows.csv', index=False)
    key = ['basis', 'method', 'K', 'arm', 'gamma']
    n_arms = R[(R.method != 'base') & (R.part == 'full')].groupby(['line', 'fit_half']).size().to_dict()
    res = dict(rows_per_line={str(k): int(v) for k, v in n_arms.items()}, fits={}, best_per_method={}, guardrails_base={})
    H = R[R.line == 'H1']; base = H[H.method == 'base'].set_index('part')
    for l in D_LINES[1:]:
        res['guardrails_base'][l] = comp_dict(R[(R.line == l) & (R.method == 'base')].iloc[0])
    picks = {}
    for f in (0, 1):
        Hf = H[(H.fit_half == f)]
        fit = Hf[Hf.part == f'half{f}']; elig = fit[(fit[COMPS] > 0).all(1)]
        if not len(elig):
            res['fits'][f'fit_half{f}'] = dict(picked=None); continue
        pk = elig.sort_values('overall', ascending=False).iloc[0]
        sel = np.ones(len(Hf), bool)
        for k in key: sel &= (Hf[k] == pk[k]).to_numpy()
        ho = Hf[sel & (Hf.part == f'half{1 - f}').to_numpy()].iloc[0]; b0 = base.loc[f'half{1 - f}']
        gain = {c: float(ho[c] - b0[c]) for c in COMPS + ['overall']}
        Ks_b = sorted(H[(H.basis == pk.basis)].K.unique())
        gr = {}
        for l in D_LINES[1:]:
            L = R[(R.line == l) & (R.basis == pk.basis) & (R.method == pk.method) & (R.arm == pk.arm) & (R.gamma == pk.gamma)]
            Kl = sorted(L.K.unique()); Kg = Kl[Ks_b.index(pk.K)] if pk.K in Ks_b and len(Kl) == len(Ks_b) else pk.K
            x = L[L.K == Kg].iloc[0]; y = R[(R.line == l) & (R.method == 'base')].iloc[0]
            gr[l] = dict(K=int(Kg), pick=comp_dict(x), base=comp_dict(y), diff=float(x.overall - y.overall), pds_diff=float(x.pds - y.pds))
        sum_c = sum(gain[c] for c in COMPS)
        res['fits'][f'fit_half{f}'] = dict(
            picked={k: (pk[k] if isinstance(pk[k], str) else float(pk[k])) for k in key} | dict(variant=pk.variant), eligible=int(len(elig)),
            grid_position=dict(basis=grid_pos(list(F_BASES), pk.basis), method=grid_pos(F5_METHODS, pk.method), K=grid_pos(Ks_b, pk.K),
                               arm=grid_pos(['R2', 'R3'], pk.arm), gamma=grid_pos([0.25, 0.5, 0.75, 1.0], pk.gamma)),
            held_out=comp_dict(ho), base_held_out=comp_dict(b0), gain=gain, reach_share_of_component_gain=float(gain['reach'] / sum_c) if sum_c else np.nan,
            guardrails=gr)
        picks[f] = (pk.basis, pk.method, pk.variant, int(pk.K))
        log(f'F5 fit half{f}: pick {res["fits"][f"fit_half{f}"]["picked"]} held-out {ho.overall:.4f} vs base {b0.overall:.4f}; guardrails ' +
            ', '.join(f'{l} {v["diff"]:+.4f}' for l, v in gr.items()))
        for m in F5_METHODS:
            ho_m = Hf[(Hf.method == m) & (Hf.part == f'half{1 - f}')]
            if len(ho_m):
                bm = ho_m.sort_values('overall').iloc[-1]
                res['best_per_method'].setdefault(m, {})[f'heldout_half{1 - f}'] = dict(config={k: (bm[k] if isinstance(bm[k], str) else float(bm[k])) for k in key} | dict(variant=bm.variant),
                                                                                     overall=float(bm.overall), gain=float(bm.overall - b0.overall), pds=float(bm.pds))
    fv = [v for v in res['fits'].values() if v.get('picked')]
    res['criterion'] = dict(both_picked=len(fv) == 2, heldout_gain_ge_0005_both=bool(len(fv) == 2 and all(v['gain']['overall'] >= 0.005 for v in fv)),
                            all_components_pos_both=bool(len(fv) == 2 and all(all(v['held_out'][c] > 0 for c in COMPS) for v in fv)),
                            guardrails_within_0005=bool(len(fv) == 2 and all(g['diff'] >= -0.005 for v in fv for g in v['guardrails'].values())))
    res['criterion']['met'] = all(res['criterion'].values())
    E = json.loads((S19T / 'partEreg.json').read_text())
    res['step19_ereg_fits'] = {k: dict(picked=v.get('picked'), held_out=v.get('held_out'), gain=v.get('gain')) for k, v in E['fits'].items()}
    F4 = pd.read_csv(TAB / 'partF4_rows.csv')
    res['F4_variants_at_picks'] = {}
    for f, (b, m, v, K) in picks.items():
        d = F4[(F4.basis == b) & (F4.K == K)]
        res['F4_variants_at_picks'][f'fit_half{f}'] = d[['method', 'median_r2_heldout'] + [c for c in d.columns if c.startswith('c_')]].to_dict('records')
    # D/E/F P_prog for each distinct pick
    ctx = FCtx()
    for f, pk in picks.items():
        if f == 1 and picks.get(0) == pk:
            continue
        b, m, v, K = pk; out = def_pprog(ctx, b, m, v, K)
        np.savez_compressed(TAB / f'partF5_DEF_Pprog_{b}_{v}_K{K}.npz', **out)
        res.setdefault('DEF', {})[f'{b}_{v}_K{K}'] = dict(file=f'tables/partF5_DEF_Pprog_{b}_{v}_K{K}.npz', from_fit_half=f,
                                                       targets_with_scores=int((out['n_backgrounds'] > 0).sum()), median_backgrounds=float(np.median(out['n_backgrounds'])))
    save_json(TAB / 'partF5.json', res)
    print(json.dumps(res, indent=1, default=float))



# ================================================================ Part G: the recommended spec's combination row (proxy)
def g_parts(D, spec, ctx):
    """(mu, sig) for a spec on one line. Keys: calls ('K0' | 'hybrid' with T | 'K4'); shift {'vec': 'C' | 'D', 'a': [aw, am, as],
    'cuts': 'full' | 'half0' | 'half1'}; prog {'basis', 'method', 'variant', 'K', 'arm', 'gamma'} (H1 cross-fitted: each
    half from the other half's fit; guardrails: the fold's fit at the matching K); coexp {'gamma'} (lambda per (b) tercile
    from all of H1). The no-source fix is always in (base_mu)."""
    import scorer_step16 as S16
    import scorer_step12 as S12
    import scorer_arms as SA
    h = D.ref.h; line = D.line; mu0 = base_mu(D); sig = D.sig
    bvals = b_for(B_LINE[line], h.targets)
    if spec.get('calls') in ('K4', 'hybrid'):
        if line == 'H1':
            z = np.load(C15 / 'config1p_K0/EvalA/input.npz', allow_pickle=True); k4 = z['de_sig_K4'] & ~(z['n_sources'] <= 0)
        else:
            k4 = S16.call_sets(D.G, [1200])['K4_1200']
        on = np.ones(h.T, bool) if spec['calls'] == 'K4' else (np.nan_to_num(bvals, nan=-np.inf) >= spec['T'])
        sig = np.where(on[:, None], k4, D.sig)
    H1ref = ref_for('EvalA') if line == 'H1' else None
    hv = S16.halves(S12.truth_calls(D.ref)) if line == 'H1' else None
    bH = b_for('EvalA', (H1ref or ref_for('EvalA')).h.targets)
    halvesH = hv if hv is not None else S16.halves(S12.truth_calls(ref_for('EvalA')))
    add = np.zeros((h.T, h.G)); mu = mu0
    pr = spec.get('prog')
    if pr:
        b, K, v = pr['basis'], int(pr['K']), pr['variant']
        if line == 'H1':
            sh = np.zeros((h.T, K))
            for f in (0, 1):
                sh[hv[1 - f]] = ctx['preds'][(b, K, v, 'H1', f)][hv[1 - f]]
            Vt = np.load(basis_file(b, 'full'))['Vt'].astype(np.float64)
        else:
            fold = FOLD_OF[line]; Kl = sorted({k[1] for k in ctx['preds'] if k[0] == b and k[3] == line})
            KH = sorted({k[1] for k in ctx['preds'] if k[0] == b and k[3] == 'H1'})
            Kf = Kl[KH.index(K)] if K in KH and len(Kl) == len(KH) else K
            sh = ctx['preds'][(b, Kf, v, line, -1)]; Vt = np.load(basis_file(b, fold))['Vt'].astype(np.float64); K = Kf
        P = f5_P(Vt, b, None, K, sh, D.gate)
        on = np.ones(h.T, bool) if pr['arm'] == 'R2' else (np.nan_to_num(bvals, nan=-np.inf) >= 5)
        P = np.where(on[:, None], pr['gamma'] * P, 0.0)
        if SPACE[b] == 'lfc':
            mu = SA.mu_from_lfc(np.log2((mu0 + 1e-9) / (D.m0[None, :] + 1e-9)) + P, D.m0)
        else:
            add += P
    sh_ = spec.get('shift')
    if sh_:
        cuts = b_cuts(bH if sh_['cuts'] == 'full' else bH[halvesH[int(sh_['cuts'][-1])]])
        t = b_terc(bvals, cuts)
        if sh_['vec'] == 'C':
            S = np.load(TAB / 'partC_S.npz')[f'S_{line}']
        else:
            S = np.load(TAB / 'partD_S_H1_clean.npz')['S_H1_clean']
        add += np.array(sh_['a'])[t][:, None] * np.where(np.isfinite(S), S, 0.0)[None, :]
    co = spec.get('coexp')
    if co:
        z = np.load(TAB / f'partF6_prior_{line}.npz'); F6 = json.loads((TAB / 'partF6.json').read_text())
        lam = np.array(F6['lambda_']['all_H1']['lam']); t = b_terc(bvals, np.array(F6['lambda_']['all_H1']['cuts']))
        ok = np.broadcast_to(D.gate[None, :], (h.T, h.G)).copy()
        for i in h.rows:
            if h.own[i] >= 0: ok[i, h.own[i]] = False
        add += co['gamma'] * lam[t][:, None] * np.where(ok & z['has_prior'][:, None], z['minus_C_k'], 0.0)
    if np.any(add != 0):
        mu = mu_add_delta(mu, add, own_index(h))
    return mu, sig


def partG(a):
    """Score a list of specs (JSON list in --spec) on H1 (full, halves) and the three fold pools: base, each spec."""
    import pickle
    specs = json.loads(a.spec) if a.spec not in ('{}', '') else json.loads((TAB / 'partG_specs.json').read_text())
    z = pickle.load(open(TAB / 'partF4_preds.pkl', 'rb')) if (TAB / 'partF4_preds.pkl').exists() else dict(preds={})
    ctx = dict(preds=z['preds'])
    f = TAB / 'partG_rows.csv'; rows = pd.read_csv(f).to_dict('records') if f.exists() else []
    for line in D_LINES:
        D = DLine(line)
        for sp in [dict(name='base')] + specs:
            if any(r['line'] == line and r['row'] == sp['name'] for r in rows):
                continue
            mu, sig = g_parts(D, sp, ctx)
            rr, _ = score_rows(D, mu, sig, row=sp['name'], spec=json.dumps(sp)); rows += rr
            log(f'G {line} {sp["name"]}: ' + ', '.join(f'{r["part"]} {r["overall"]:.4f}' for r in rr))
            ref19 = {'base': 'nosource_fix', 'hybrid T5 + fix (spec)': 'hybrid_T5+nosource_fix'}.get(sp['name'])
            if ref19:                                     # must reproduce Step 19 Part F's row
                F19 = pd.read_csv(S19T / 'partF_rows.csv'); w = F19[(F19.line == line) & (F19.row == ref19)].set_index('part')['overall']
                d = max(abs(r['overall'] - w[r['part']]) for r in rr)
                log(f'G check {line} {sp["name"]} vs Step 19 {ref19}: max |diff| {d:.1e}')
                if d > 1e-9:
                    log('G CHECK FAILED'); sys.exit(6)
            pd.DataFrame(rows).to_csv(f, index=False)
        del D; gc.collect()


# ================================================================ unit tests (small slices; no table load, no proxy)
def test(a):
    rng = np.random.default_rng(0)
    # (b) cut points and absolute terciles
    b = np.array([np.nan, 0, 1, 2, 3, 4, 5, 6, 7, 8.0]); c = b_cuts(b)
    assert np.allclose(c, np.quantile(np.arange(9.0), [1 / 3, 2 / 3])), c
    assert b_terc(b, c).tolist() == [0, 0, 0, 0, 1, 1, 1, 2, 2, 2], b_terc(b, c)
    assert b_terc(np.array([100.0, -1, np.nan]), c).tolist() == [2, 0, 0]
    # Delta-space addition: zero term is the identity; a term round-trips through the profile; own gene untouched
    import scorer_arms as SA
    mu = rng.gamma(2, 50, (3, 6)); add = np.zeros((3, 6))
    assert np.allclose(mu_add_delta(mu, add, np.array([-1, 2, 0])), mu, rtol=1e-12)
    add = rng.normal(0, 0.2, (3, 6)); m2 = mu_add_delta(mu, add, np.array([-1, 2, 0]))
    d = prof(m2) - prof(mu); A = add.copy(); A[1, 2] = 0; A[2, 0] = 0
    assert np.allclose(np.where(prof(mu) + A > 0, d, A), A, atol=1e-12) and np.isclose(m2[1, 2], mu[1, 2], rtol=1e-12) and np.isclose(m2[2, 0], mu[2, 0], rtol=1e-12)
    assert np.allclose(SA.mu_from_profile(prof(mu)), mu)
    # pair cosine over a mask, duplicates collapse, pair classes
    A_ = np.array([[1.0, 2, 3], [1, 0, 0]]); B_ = np.array([[2.0, 4, 100], [0, 1, 0]]); M = np.array([[1, 1, 0], [1, 1, 1]], bool)
    cc = pair_cos(A_, B_, M); assert np.isclose(cc[0], 1) and np.isclose(cc[1], 0), cc
    u, X = collapse(np.array(['b', 'a', 'b']), np.array([[1.0, np.nan], [5, 5], [3, 2]], np.float32))
    assert u.tolist() == ['a', 'b'] and np.allclose(X, [[5, 5], [2, 2]]), (u, X)
    assert pair_class('W', 'K562', 'blood', 'W', 'RPE1', 'retina') == 'same lab, different line'
    assert pair_class('W', 'K562', 'blood', 'W', 'Jurkat', 'blood') == 'same lab, same lineage, different line'
    assert pair_class('W', 'K562', 'blood', 'J', 'K562', 'blood') == 'different lab, same line'
    assert pair_class('A', 'H1', 'pluripotent', 'F', 'eipl', 'pluripotent') == 'different lab, same lineage'
    assert pair_class('A', 'H1', 'pluripotent', 'W', 'K562', 'blood') == 'different lab, different lineage'
    assert pair_class('F', 'eipl', 'pluripotent', 'F', 'eipl', 'pluripotent') == 'same lab, same line'
    assert line_of('eipl_3') == 'eipl' and line_of('KOLF2.1J') == 'kolf' and line_of('kolf_2') == 'kolf' and line_of('HipSci_pooled') == 'HipSci_pooled'
    assert lineage_of('eipl_3') == 'pluripotent' and lineage_of('K562') == 'blood' and lab_of('nadig_jurkat') == 'Weissman'
    # regression recovers known coefficients; the table bootstrap returns finite CIs
    n = 400; tabs = [f't{i}' for i in range(30)]
    P = pd.DataFrame(dict(table1=rng.choice(tabs, n), table2=rng.choice(tabs, n), lab1=rng.choice(['x', 'y'], n), lab2=rng.choice(['x', 'y'], n),
                          line1=rng.choice(['p', 'q', 'r'], n), line2=rng.choice(['p', 'q', 'r'], n), cells1=rng.uniform(20, 400, n), cells2=rng.uniform(20, 400, n)))
    P['lineage1'] = P.line1.map({'p': 'L', 'q': 'L', 'r': 'M'}); P['lineage2'] = P.line2.map({'p': 'L', 'q': 'L', 'r': 'M'})
    X, names = reg_design(P); truth = np.array([0.1, 0.2, 0.05, 0.03, 0.02, 0.01])
    P['y'] = X @ truth + rng.normal(0, 0.01, n)
    r = pair_regression(P, 'y', n_boot=200)
    assert all(abs(r['coef'][nm]['beta'] - t) < 0.01 for nm, t in zip(names, truth)), r
    assert all(r['coef'][nm]['ci95'][0] < r['coef'][nm]['beta'] < r['coef'][nm]['ci95'][1] for nm in names), r
    # weighted votes: only the named backgrounds' SE change, by 1/sqrt(w)
    class _V:
        groups = {'a': [], 'b': []}; owners = {'k': ['a', 'b']}
        def gather(self, k, g): return [np.ones(2), np.ones(2)], [np.full(2, 2.0), np.full(2, 2.0)], [1.0, 1.0]
    Y, SE, _ = MultiWeightedVotes(_V(), {'b'}, 4).gather('k', None)
    assert np.allclose(SE[0], 2) and np.allclose(SE[1], 1)
    class _D:
        def gather(self, k, g): return ([np.array([1.0, 1, np.nan]), np.array([3.0, -1, 2])], None, None) if k == 'k' else ([], None, None)
    nn, zz = agreement(_D(), 'k', ['g1', 'g2', 'g3'])
    assert nn.tolist() == [2, 2, 1] and np.isclose(zz[0], 2 / (np.sqrt(2) / np.sqrt(2))) and np.isclose(zz[1], 0) and np.isnan(zz[2]), (nn, zz)
    assert agreement(_D(), 'z', ['g1'])[0].tolist() == [0]
    # strong-minus-weak, S with exclusions, projection
    X = np.array([[1.0, 0], [2, 0], [3, 0], [10, 1], [11, np.nan], [12, 3]]); nr = np.array([0, 1, 2, 50, 60, 70.0])
    t = tercile_of(nr); assert t.tolist() == [0, 0, 1, 1, 2, 2]
    v = strong_minus_weak(X, np.isfinite(X), nr); assert np.allclose(v, [11.5 - 1.5, 3 - 0]), v
    vec = np.array([[1.0, np.nan], [3, 4], [5, 6]]); assert np.allclose(shift_S(['a', 'b', 'c'], vec, ('c',)), [2, 4])
    assert np.allclose(shift_S(['a', 'b', 'c'], vec, ()), [3, 5])
    vec4 = np.array([[1.0, np.nan, 5], [3, np.nan, np.nan], [5, 2, np.nan], [7, 4, 9]])
    s2 = shift_S2(['a', 'b', 'c', 'd'], vec4, ('d',), ['a', 'b', 'c', 'd'])        # 3 voting lines, need 2
    assert np.allclose(s2[0], 3) and np.isnan(s2[1]) and np.isnan(s2[2]), s2
    s3 = shift_S2(['a', 'b', 'c', 'd'], vec4, (), ['a', 'c', 'd'])                 # 3 kept lines, need 2
    assert np.allclose(s3, [13 / 3, 3, 7]), s3
    w = unit(rng.normal(size=8)); s_ = rng.normal(size=8); sc = s_ - (s_ @ w) * w; assert abs(sc @ w) < 1e-12
    # pick_fit: eligibility (all six > 0 on the fit part), held-out row, base row
    R = pd.DataFrame([dict(part=p, x=x, overall=o, **{c: (v if c == 'pds' else 0.1) for c in COMPS})
                      for p, x, o, v in [('h0', 0, .2, .1), ('h0', 1, .3, -.1), ('h0', 2, .25, .1), ('h1', 0, .1, .1), ('h1', 1, .5, .1), ('h1', 2, .15, .1)]])
    pk, ho, b0, ne = pick_fit(R, 'h0', 'h1', ['x'], lambda d: d.x == 0)
    assert pk.x == 2 and ho.overall == .15 and b0.overall == .1 and ne == 2
    log('test: all unit tests pass')



def testF(a):
    rng = np.random.default_rng(1)
    assert np.allclose(f_weight(np.array([10, 50, 100.0])), [1 / np.sqrt(50), 1 / np.sqrt(50), 0.1])
    # ridge: recovers coefficients at low noise, picks a small alpha; weights equal duplicated rows
    J, n, p = 3, 300, 4
    X = rng.normal(size=(J, n, p)) * np.array([1, 2, 0.5, 3])[None, None, :] + 1.0
    B = rng.normal(size=(J, p)); Y = np.einsum('jnp,jp->nj', X, B) + 0.5 + rng.normal(0, 0.01, (n, J))
    m = ridge_cv(X, Y, np.ones(n))
    pr = ridge_predict(m, X); assert np.abs(pr - Y).max() < 0.1 and (m['alpha'] <= 1).all() and (m['cv_r2'] > 0.99).all(), (m['alpha'], m['cv_r2'])
    w = rng.integers(1, 3, n).astype(float); idx = np.repeat(np.arange(n), w.astype(int))
    m1 = ridge_solve(X, Y, w, np.full(J, 10.0)); m2 = ridge_solve(X[:, idx], Y[idx], np.ones(len(idx)), np.full(J, 10.0 * len(idx) / w.sum() * w.sum() / len(idx)))
    assert np.allclose(ridge_predict(m1, X), ridge_predict(m2, X), atol=1e-8), 'weighted ridge differs from duplicated rows'
    Xc = X.copy(); Xc[:, :, 1] = 7.0; mc = ridge_solve(Xc, Y, np.ones(n), np.ones(J)); assert np.isfinite(ridge_predict(mc, Xc)).all()
    # build_X: source scores, missing indicators, gene features
    A = np.array([[[1.0, 2], [np.nan, np.nan]], [[3, 4], [5, 6]]])        # (n=2, B=2, K=2)
    Xm = build_X('M2', A, None); assert Xm.shape == (2, 2, 4) and np.allclose(Xm[0], [[1, 0, 0, 1], [3, 5, 0, 0]]) and np.allclose(Xm[1], [[2, 0, 0, 1], [4, 6, 0, 0]])
    GF = np.array([[9.0, 1], [8, 0]]); assert build_X('M3', A, GF).shape == (2, 2, 6) and np.allclose(build_X('M4', A, GF)[1], GF)
    sh, nb = avg_scores(A); assert np.allclose(sh, [[1, 2], [4, 5]]) and nb.tolist() == [1, 2]
    # program_truth equals per-row least squares
    Vt = np.linalg.qr(rng.normal(size=(30, 6)))[0].T; Xt = rng.normal(size=(5, 30)); gated = rng.random((5, 30)) > 0.3
    tr = dict(X=np.where(gated, Xt, np.nan), gated=gated, Delta=Xt)
    pt = program_truth(tr, Vt, [2, 6], 'lfc')
    for K in (2, 6):
        assert np.allclose(pt[K], S19.e_coefs(Vt, K, tr['X'], gated)), K
    assert np.allclose(program_truth(tr, Vt, [3], 'delta')[3], S19.e_coefs(Vt, 3, np.where(gated, Xt, np.nan), gated))
    # bg_scores: tables averaged within a background first; exclusions; missing NaN
    mt = pd.DataFrame(dict(background=['P', 'P', 'Q', 'H1'], table=['p1', 'p2', 'q1', 'h'], target=['a', 'a', 'a', 'a']))
    RSt = np.array([[1.0, 0], [3, 0], [10, 2], [1000, 0]], np.float32)
    Ab = bg_scores(RSt, mt, ['a', 'z'], ('H1',), 2, ['P', 'Q', 'H1'])
    assert np.allclose(Ab[0, 0], [2, 0]) and np.allclose(Ab[0, 1], [10, 2]) and np.isnan(Ab[0, 2]).all() and np.isnan(Ab[1]).all()
    # pair cosine with both own genes out equals the brute force
    Dg = rng.normal(size=(4, 7)); own = np.array([0, 3, -1, 3])
    N = pair_cos_both_out(Dg, own)
    for i in range(4):
        for j in range(i + 1, 4):
            m_ = np.ones(7, bool); m_[[g for g in (own[i], own[j]) if g >= 0]] = False
            x, y = Dg[i, m_], Dg[j, m_]; assert np.isclose(N[i, j], x @ y / np.sqrt((x @ x) * (y @ y))), (i, j)
    # residual correlation: library size regressed out, unit diagonal, constant gene 0
    ll = rng.normal(size=50); g0 = rng.normal(size=50); Mc = np.column_stack([g0 + 3 * ll, g0 + 0.1 * rng.normal(size=50) - 2 * ll, np.ones(50)])
    Cc, sd_ = resid_corr(Mc, ll); assert np.isclose(Cc[0, 0], 1) and Cc[0, 1] > 0.95 and Cc[2, 2] == 0 and not np.isfinite(sd_[2])
    # energy_delta, small-background artefacts, HGNC mapping
    trd = dict(gated=np.ones((2, 3), bool), Delta=np.array([[1.0, 0, 0], [0, 2, 0]]))
    assert np.allclose(energy_delta(trd, np.array([[1.0, 0, 0], [0, 1, 0]])), [[0, 1], [1, 4]])
    US = np.array([[10.0], [0.1], [0.1]]); mm = pd.DataFrame(dict(background=['s', 'b', 'b']))
    sa = small_bg_artefacts(np.repeat(US, 10, 1), mm, {'s': 4, 'b': 2000}, top=1); assert sa[0]['small_bg_in_top12'] == 1 and sa[0]['small_bg_energy_share'] > 0.99
    Hm = dict(appr={'H1-0', 'TP53'}, prev={'H1F0': 'H1-0'}, alias={'P53': 'TP53'})
    assert to_hgnc('H1F0', Hm) == 'H1-0' and to_hgnc('P53', Hm) == 'TP53' and to_hgnc('TP53', Hm) == 'TP53' and to_hgnc('XX', Hm) is None
    fsd = dict(idx={'TP53': 0, 'H1-0': 1}, emb=np.array([[1.0] * 32, [2.0] * 32]), extra=np.array([5.0, 6]), ok=np.array([True, False]))
    Gx, Ec, hs = gene_features(['P53', 'H1F0', 'ZZ'], fsd, Hm)
    assert hs.tolist() == [True, False, False] and Gx[0, 32] == 5 and Gx[1, 32] == 6 and Gx[1, 33] == 1 and Gx[2, 33] == 1 and (Gx[1, :32] == 0).all()
    log('testF: all Part F unit tests pass')


CMDS = {'test': test, 'testF': testF, 'part0': part0, 'partA': partA, 'partB': partB, 'partBscore': partBscore, 'partC': partC, 'partC2': partC2, 'partD': partD,
        'partF1': partF1, 'partF2': partF2, 'partF3a': partF3a, 'partF3b': partF3b, 'partF3c': partF3c, 'partF3sum': partF3sum,
        'partF3d': partF3d, 'partF4': partF4, 'partF5time': partF5time, 'partF5': partF5, 'partF5sum': partF5sum, 'partF6': partF6, 'partG': partG}

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd')
    p.add_argument('--line'); p.add_argument('--worker', type=int, default=-1); p.add_argument('--spec', default='{}')
    a = p.parse_args()
    log(f'step20 {a.cmd} start; OMP_NUM_THREADS={os.environ.get("OMP_NUM_THREADS")} host {os.uname().nodename.split(".")[0]}')
    CMDS[a.cmd](a)
    log(f'step20 {a.cmd} done; process peak {rss_gb():.1f} GB')

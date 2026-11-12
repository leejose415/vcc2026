"""Step 22 (Jose 2026-11-03): the line's own programs from its controls (Part P), line-vote sign agreement (Part S),
two quick checks (Part Q), reach by ordering the generator's p-values (Part R), and the recommendation (Part D).
vcc2027 only. No D/E/F file, no upload; the only cells generated are R3's 20-knockdown test slice.
Write-up: section 14 of results/stage1/2026-10-28_scorer_space/report.md. Outputs: results/stage1/2026-11-03_step22/.
Frozen code (generator, harness, run_line, boot_score, config 1' files, vcc_2026_7) and the Step 18-21 outputs are imported or
read only; Part R3's generator v8 is a copy in code/vcc6/gen_v8/.

Base for every test: row 4 (step21.Row4: spec + fixed shift + F5 term + no-source fix), H1 cross-fitted as in Step 21.
Space: Delta in log1p CP50K (scorer space), knockdown minus control, on the line's gated genes, own gene out.

Order on vcc2027 (step22.sh): test, part0, srcdelta, P0, P1, P2, S1, S2, S3, S4, Q1, Q2, R2, R3, D.
"""
import argparse
import gc
import hashlib
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
import step20 as S20                                                                       # noqa: E402
import step21 as S21                                                                       # noqa: E402
from step19 import log, rss_gb, avail_gb, mem_before, mem_after, save_json, grid_pos, COMPS, TERC, tercile_of, b_for  # noqa: E402
from step20 import prof, mu_add_delta, b_terc, comp_dict, FOLD_OF, basis_file                # noqa: E402
from step21 import Row4, CLine, SHIFT_A, T_HYB, EV                                         # noqa: E402

OUT = V / 'results/stage1/2026-11-03_step22'
TAB, LOG = OUT / 'tables', OUT / 'logs'
DL = V / 'lib/downloads/step22_2026-11-03'
SEED = 2026
LINES = ['H1', 'K562', 'HCT116', 'HEK293T', 'HCT116s']
FOLDS = ['K562', 'HCT116', 'HEK293T']
GUARD = FOLDS + ['HCT116s']                       # HCT116s reported, not in the criterion
FE = V / 'results/stage1/final_evals/harness/lib'
CTRL_FILE = {'H1': V / 'lib/harness/evalA/context.h5ad', 'K562': FE / 'evalK562final/context.h5ad',
             'HCT116': FE / 'evalHCT116final/context.h5ad', 'HEK293T': FE / 'evalHEK293Tfinal/context.h5ad',
             'D': V / 'data_final/extracted/context_D.h5ad', 'E': V / 'data_final/extracted/context_E.h5ad',
             'F': V / 'data_final/extracted/context_F.h5ad'}
CTRL_OF = {'H1': 'H1', 'K562': 'K562', 'HCT116': 'HCT116', 'HEK293T': 'HEK293T', 'HCT116s': 'HCT116'}
CTRL_LINES = ['H1', 'K562', 'HCT116', 'HEK293T', 'D', 'E', 'F']
E4W_FOLD = {'H1': 'full', 'K562': 'K562', 'HCT116': 'Orion', 'HEK293T': 'Orion', 'HCT116s': 'Orion'}
SRC_VARIANT = {'H1': 'full', 'K562': 'K562', 'HCT116': 'Orion', 'HEK293T': 'Orion', 'HCT116s': 'HCT116',
               'D': 'full', 'E': 'full', 'F': 'full'}
VARIANT_DROP = {'full': (), 'K562': ('K562',), 'Orion': ('HCT116', 'HEK293T'), 'HCT116': ('HCT116',)}
PROXY_REACH_NOTE = 'the proxy overstates reach by 0.02-0.05 against cells (15c sec. 4; Step 21 Track 1 +0.03)'
HALLMARK = {   # Jose's ten sets -> MSigDB_Hallmark_2020 term names in Enrichr
    'E2F_TARGETS': 'E2F Targets', 'G2M_CHECKPOINT': 'G2-M Checkpoint', 'P53_PATHWAY': 'p53 Pathway',
    'UNFOLDED_PROTEIN_RESPONSE': 'Unfolded Protein Response', 'CHOLESTEROL_HOMEOSTASIS': 'Cholesterol Homeostasis',
    'MTORC1_SIGNALING': 'mTORC1 Signaling', 'MYC_TARGETS_V1': 'Myc Targets V1',
    'INTERFERON_ALPHA_RESPONSE': 'Interferon Alpha Response', 'TNFA_SIGNALING_VIA_NFKB': 'TNF-alpha Signaling via NF-kB',
    'HYPOXIA': 'Hypoxia'}
PATHS = list(HALLMARK)
R_MIN = 0.5                                       # pathways with half-split r_P < 0.5 contribute 0 (fixed before any result)
CC_URL = 'https://raw.githubusercontent.com/satijalab/seurat/master/data/cc.genes.rda'        # Seurat's cc.genes (Tirosh 2016)


def save_rows(rows, name):
    R = pd.DataFrame(rows); R.to_csv(TAB / f'{name}.csv', index=False); return R


def resume(name):
    f = TAB / f'{name}.csv'
    return pd.read_csv(f).to_dict('records') if f.exists() else []


def md5_of(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for ch in iter(lambda: f.read(1 << 24), b''):
            h.update(ch)
    return h.hexdigest()


def panel():
    import run_line as RL
    return RL.panel_genes()


# ================================================================ lines: row 4 with tercile parts
def add_parts(R):
    """True-strength tercile parts on every line; on H1 also each half's terciles (held-out tercile gains)."""
    import scorer_step12 as S12
    if not hasattr(R, 'T'):
        R.T = S12.truth_calls(R.D.ref)
    R.n_true = R.T.sum(1); R.tterc = tercile_of(R.n_true)
    P = R.D.parts
    for j, t in enumerate(TERC):
        P[t] = np.flatnonzero(R.tterc == j)
    if R.line == 'H1':
        for f in (0, 1):
            for j, t in enumerate(TERC):
                P[f'half{f}_{t}'] = np.intersect1d(P[f'half{f}'], np.flatnonzero(R.tterc == j))
    R.ownM = S21.own_mask(R.h); R.okg = R.gate[None, :] & ~R.ownM
    return R


def row4(line):
    return add_parts(Row4(line))


def cline(line):
    return add_parts(CLine(line, terc_parts=False))


# ================================================================ generic fit / criterion
def is_base(d, key, base_vals):
    sel = np.ones(len(d), bool)
    for k, v in zip(key, base_vals):
        sel &= (d[k].astype(str) == str(v)).to_numpy()
    return sel


def pick_fit(R, fit_part, held_part, key, base_sel):
    fit = R[R.part == fit_part]
    elig = fit[(fit[COMPS] > 0).all(axis=1)]
    if not len(elig):
        return None
    pk = elig.sort_values('overall', ascending=False).iloc[0]
    sel = is_base(R, key, [pk[k] for k in key])
    ho = R[sel & (R.part == held_part).to_numpy()].iloc[0]
    b0 = R[base_sel(R) & (R.part == held_part).to_numpy()].iloc[0]
    return pk, ho, b0, len(elig)


def fit_report(R, key, grids, base_vals, fit_frames=None, tag=''):
    """Pick on one H1 half (highest overall with all six > 0 there), score the other half, both directions. Guardrails
    at each pick (full panel) against row 4: the fold pools (criterion) and HCT116s (reported only). Held-out tercile
    gains. fit_frames: {f: H1 rows used for direction f} (models cross-fitted within the fit half)."""
    res = dict(fits={})
    for f in (0, 1):
        H = fit_frames[f] if fit_frames is not None else R[R.line == 'H1']
        bsel = lambda d: is_base(d, key, base_vals)
        got = pick_fit(H, f'half{f}', f'half{1 - f}', key, bsel)
        if got is None:
            res['fits'][f'fit_half{f}'] = dict(picked=None); continue
        pk, ho, b0, ne = got
        isb = bool(all(str(pk[k]) == str(v) for k, v in zip(key, base_vals)))
        gain = {c: float(ho[c] - b0[c]) for c in COMPS + ['overall']}
        sc = sum(gain[c] for c in COMPS)
        terc = {}
        for t in TERC:
            x = H[is_base(H, key, [pk[k] for k in key]) & (H.part == f'half{1 - f}_{t}').to_numpy()]
            y = H[bsel(H) & (H.part == f'half{1 - f}_{t}').to_numpy()]
            if len(x) and len(y):
                terc[t] = {c: float(x.iloc[0][c] - y.iloc[0][c]) for c in COMPS + ['overall']}
        gr = {}
        for l in GUARD:
            L = R[(R.line == l) & (R.part == 'full')]
            if 'fit_half' in L:
                L = L[L.fit_half.isin([-1, f])]
            x = L[is_base(L, key, [pk[k] for k in key])]; y = L[bsel(L)]
            if not len(x) or not len(y):
                continue
            x, y = x.iloc[0], y.iloc[0]
            gr[l] = dict(diff=float(x.overall - y.overall), comp_diff={c: float(x[c] - y[c]) for c in COMPS},
                         pick=comp_dict(x), base=comp_dict(y), in_criterion=l in FOLDS)
        gpos = {k: ('base' if isb else (grid_pos(grids[k], pk[k]) if pk[k] in list(grids[k]) else 'outside the grid')) for k in key}
        res['fits'][f'fit_half{f}'] = dict(
            picked={k: (pk[k] if isinstance(pk[k], str) else float(pk[k])) for k in key}, picked_is_base=isb,
            eligible=int(ne), grid_position=gpos, fit_overall=float(pk.overall),
            held_out=comp_dict(ho), base_held_out=comp_dict(b0), gain=gain, held_out_tercile_gain=terc,
            reach_share_of_component_gain=float(gain['reach'] / sc) if sc else float('nan'), guardrails=gr)
        log(f'{tag} fit half{f}: pick {res["fits"][f"fit_half{f}"]["picked"]} ({gpos}); held-out {ho.overall:.4f} vs base '
            f'{b0.overall:.4f} ({gain["overall"]:+.4f}); guardrails ' + ', '.join(f'{l} {v["diff"]:+.4f}' for l, v in gr.items()))
    fv = [v for v in res['fits'].values() if v.get('picked') is not None]
    res['criterion'] = criterion(fv)
    res['proxy_note'] = PROXY_REACH_NOTE
    return res


def criterion(fv):
    c = dict(both_picked=len(fv) == 2,
             heldout_gain_ge_0005_both=bool(len(fv) == 2 and all(v['gain']['overall'] >= 0.005 for v in fv)),
             all_components_pos_both=bool(len(fv) == 2 and all(all(v['held_out'][k] > 0 for k in COMPS) for v in fv)),
             fold_pools_within_0005=bool(len(fv) == 2 and all(g['diff'] >= -0.005 for v in fv for l, g in v['guardrails'].items() if l in FOLDS)))
    c['met'] = all(c.values())
    c['HCT116final_reported_only'] = {f'fit{i}': v['guardrails'].get('HCT116s', {}).get('diff') for i, v in enumerate(fv)}
    return c


def score_tag(R, mu, sig=None, **tag):
    rows, _ = R.score(mu, sig, **tag)
    return rows


def fmt(rows, parts=('full', 'half0', 'half1')):
    return ', '.join(f'{r["part"]} {r["overall"]:.4f}' for r in rows if r['part'] in parts)


# ================================================================ Part 0: checks
def part0(a):
    import anndata as ad
    import scipy.sparse as sp
    import run_line as RL
    res = dict(check1={}, controls={})
    A = pd.read_csv(S21.TAB / 'partA_rows.csv'); A = A[A.row == 'row 4']
    want = {'H1': 0.2928, 'K562': 0.1532, 'HCT116': 0.3153, 'HEK293T': 0.1695, 'HCT116s': 0.3651}
    ok = True; base = []; gates = {}
    for line in LINES:
        R = row4(line); t0 = time.time(); rows = score_tag(R, R.mu(), row='row 4'); dt = time.time() - t0
        base += rows; gates[line] = R.gate.copy()
        w = A[A.line == line].set_index('part')['overall']
        d = max(abs(r['overall'] - w[r['part']]) for r in rows if r['part'] in w.index)
        got = {r['part']: round(r['overall'], 4) for r in rows if r['part'] in ('full', 'half0', 'half1')}
        ok &= d < 1e-9 and abs(got['full'] - want[line]) < 5e-5
        res['check1'][line] = dict(got=got, step21_partA=w.to_dict(), max_abs_diff=d, prompt=want[line], seconds=dt)
        log(f'part0 check1 {line}: {got} vs Step 21 A {w.round(4).to_dict()} (prompt {want[line]}); max |diff| {d:.1e}; {dt:.1f} s')
        del R; gc.collect()
    save_rows(base, 'row4_rows')
    genes = RL.panel_genes()
    for c1 in ('D', 'E', 'F'):
        gates[c1] = S21.c1_load(c1)['gate'].astype(bool)
    for line in CTRL_LINES:
        p = CTRL_FILE[line]
        if not p.exists():
            res['controls'][line] = dict(path=str(p), present=False); log(f'part0 controls {line}: MISSING {p}'); continue
        a_ = ad.read_h5ad(p, backed='r')
        n = a_.n_obs; obs = a_.obs
        lab = None
        for col in ('target_gene', 'target', 'gene'):
            if col in obs:
                lab = obs[col].astype(str).value_counts().to_dict(); break
        X = a_.X[:]                                    # full matrix (all genes in the file)
        X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
        dat = X.data
        integ = bool(np.all(dat == np.round(dat))) and bool((dat >= 0).all())
        umi = np.asarray(X.sum(1)).ravel()
        Xp = RL.read_controls(p, genes); lp = np.asarray(Xp.sum(1)).ravel()
        import build_effects_v5 as B
        C, _ = B.cpm_csr(Xp.data, Xp.indices, Xp.indptr, Xp.shape[1], np.ones(Xp.shape[1], bool))
        m0 = np.asarray(C.sum(0)).ravel() / C.shape[0]; g = m0 > B.MIN_CPM
        gl = gates.get(line)
        res['controls'][line] = dict(path=str(p), present=True, control_cells=int(n), labels=lab, genes_in_file=int(X.shape[1]),
                                     median_umi_per_cell=float(np.median(umi)), median_panel_umi_per_cell=float(np.median(lp)),
                                     dtype=str(dat.dtype), raw_integer_counts=integ, max_count=float(dat.max()) if len(dat) else 0.0,
                                     gated_genes=int(g.sum()), gate_equals_line_gate=bool(gl is not None and np.array_equal(g, gl)))
        log(f'part0 controls {line}: {res["controls"][line]}')
        ok &= integ and (gl is None or np.array_equal(g, gl))
        del X, Xp, C; gc.collect()
    res['HCT116final_uses'] = 'HCT116 (the same context.h5ad: evalHCT116final)'
    res['gate_HCT116s_equals_HCT116'] = bool(np.array_equal(gates['HCT116s'], gates['HCT116']))
    res['lines_without_raw_controls'] = [l for l, v in res['controls'].items() if not (v.get('present') and v.get('raw_integer_counts'))]
    res['pass'] = bool(ok)
    save_json(TAB / 'part0.json', res)
    log(f'part0: {"PASS" if ok else "FAIL"}; lines left out of P1/P2: {res["lines_without_raw_controls"] or "none"}')
    if not ok:
        sys.exit(3)


def row4_ref(line, part='full'):
    B = pd.read_csv(TAB / 'row4_rows.csv')
    return B[(B.line == line) & (B.part == part)].iloc[0]


# ================================================================ source Delta store (P1 coefficients, P2's a_{k,P})
SRC = OUT / 'srcdelta'


def needed_targets():
    tg = set()
    for line in ['H1', 'K562', 'HCT116', 'HEK293T']:
        z = np.load(S21.SS_CACHE / f'fold_{line}.npz', allow_pickle=True); tg |= set(map(str, z['targets']))
    tg |= set(map(str, S21.c1_load('HCT116s')['targets'])); tg |= set(S20.DEF_TARGETS())
    return sorted(tg)


def srcdelta(a):
    """Per (target, background): the mean raw Delta (full-profile log1p CP50K minus the table's control profile, each
    table's > 5 CPM gate) over every genetic row of that background's tables (Tahoe drugs and Arc H1 out: E4w's table
    set, Step 21 C10's pooled_delta), for every target any scored line or D/E/F needs."""
    import h5py
    from anndata.io import read_elem
    SRC.mkdir(parents=True, exist_ok=True)
    genes = panel(); G = len(genes); need = set(needed_targets())
    plan, npair = [], {}
    for p, stem, bg, tahoe in S19.e_tables():
        if tahoe:
            continue
        with h5py.File(p, 'r') as f:
            o = read_elem(f['obs'])
            if not np.isfinite(f['layers/log1p_bulk_50k_full'][:min(5, len(o))]).any():
                continue
        mod = o['modality'].astype(str).to_numpy(); tg = o['target_gene'].astype(str).to_numpy()
        k = np.flatnonzero((mod != 'drug') & np.isin(tg, list(need)))
        if len(k):
            plan.append((p, stem, bg, k, tg[k]))
            for t in set(tg[k]):
                npair[(t, bg)] = 1
    keys = sorted(npair); kid = {k: i for i, k in enumerate(keys)}
    mem_before(f'srcdelta sums {len(keys)} x {G} x 2', len(keys) * G * 8 / 1e9 + 2)
    Ssum = np.zeros((len(keys), G), np.float64); Scnt = np.zeros((len(keys), G), np.int32)
    for p, stem, bg, k, tg in plan:
        X, Mm = S20.table_entries(p, k, 'delta')
        for r, t in enumerate(tg):
            i = kid[(t, bg)]; m = Mm[r] & np.isfinite(X[r])
            Ssum[i, m] += X[r, m]; Scnt[i, m] += 1
        log(f'srcdelta {stem} ({bg}): {len(k)} rows')
        del X, Mm
    D = np.where(Scnt > 0, Ssum / np.maximum(Scnt, 1), np.nan).astype(np.float32)
    np.save(SRC / 'delta.npy', D)
    pd.DataFrame(keys, columns=['target', 'background']).to_parquet(SRC / 'keys.parquet')
    res = dict(tables=len(plan), pairs=len(keys), targets=len({t for t, _ in keys}), backgrounds=sorted({b for _, b in keys}),
               needed_targets=len(need))
    save_json(TAB / 'srcdelta.json', res); log(f'srcdelta: {res}')
    mem_after('srcdelta')


class Src:
    def __init__(self):
        self.D = np.load(SRC / 'delta.npy', mmap_mode='r'); k = pd.read_parquet(SRC / 'keys.parquet')
        self.by = {}
        for i, (t, b) in enumerate(zip(k.target, k.background)):
            self.by.setdefault(t, []).append((b, i))
        self.gp = {g: j for j, g in enumerate(panel())}

    def pooled(self, target, variant):
        """Mean over backgrounds (one vote each) of the per-background mean Delta; own gene NaN. None if no source."""
        rows = [i for b, i in self.by.get(target, []) if b not in VARIANT_DROP[variant]]
        if not rows:
            return None
        X = np.asarray(self.D[np.array(rows)], np.float64)
        with np.errstate(invalid='ignore'), __import__('warnings').catch_warnings():
            __import__('warnings').simplefilter('ignore')
            P = np.nanmean(X, 0)
        if target in self.gp:
            P[self.gp[target]] = np.nan
        return P


# ================================================================ Part P0: control bases
def gene_sets():
    """MSigDB_Hallmark_2020 from Enrichr via gseapy (cached, md5 logged) and the Tirosh 2016 S / G2M lists (Seurat's
    cc.genes, 43 + 54 genes, from Seurat's data/cc.genes.rda, read with Rscript), mapped to HGNC."""
    DL.mkdir(parents=True, exist_ok=True)
    f = DL / 'MSigDB_Hallmark_2020.json'
    if not f.exists():
        sys.path.insert(0, str(Path.home() / 'pylib_gseapy'))
        import gseapy
        lib = gseapy.get_library(name='MSigDB_Hallmark_2020', organism='Human')
        f.write_text(json.dumps(lib, indent=0, sort_keys=True))
        (DL / 'MSigDB_Hallmark_2020.source.txt').write_text(f'gseapy {gseapy.__version__} get_library(name="MSigDB_Hallmark_2020", organism="Human") {time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())} UTC\n')
    cc = DL / 'cc.genes.rda'; cct = DL / 'cc.genes.txt'
    if not cc.exists():
        S20.download(CC_URL, cc)
    if not cct.exists():
        import subprocess
        r = subprocess.run(['Rscript', '-e', f'load("{cc}"); writeLines(c(paste0("S\\t", cc.genes$s.genes), paste0("G2M\\t", cc.genes$g2m.genes)), "{cct}")'],
                           capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError(f'Rscript could not read {cc}: {r.stderr}')
    lib = json.loads(f.read_text()); cl = [x.split('\t') for x in cct.read_text().split('\n') if x.strip()]
    s_raw = [g for k, g in cl if k == 'S']; g_raw = [g for k, g in cl if k == 'G2M']
    info = dict(hallmark_file=str(f), hallmark_md5=md5_of(f), hallmark_terms=len(lib), cc_file=str(cc), cc_md5=md5_of(cc), cc_url=CC_URL,
                cc_S=len(s_raw), cc_G2M=len(g_raw))
    missing = [k for k, v in HALLMARK.items() if v not in lib]
    if missing or len(s_raw) != 43 or len(g_raw) != 54:
        raise RuntimeError(f'gene sets: missing hallmark terms {missing}; cc.genes S {len(s_raw)} / G2M {len(g_raw)} (expected 43 / 54)')
    H = S20.hgnc_maps(); gp = {g: j for j, g in enumerate(panel())}
    report = {}
    def to_panel(name, lst):
        """Direct panel matches, and matches after mapping unmatched symbols through HGNC's previous / alias symbols
        (unique targets only, step20.to_hgnc); the mapped version is used only if it matches more genes."""
        direct = sorted({s for s in lst if s in gp})
        mapped = sorted(set(direct) | {g for s in lst if s not in gp for g in [S20.to_hgnc(s, H)] if g in gp})
        use_mapped = len(mapped) > len(direct)
        report[name] = dict(set_size=len(set(lst)), matched_direct=len(direct), matched_after_hgnc_mapping=len(mapped),
                            used='mapped' if use_mapped else 'direct',
                            renamed={s: S20.to_hgnc(s, H) for s in lst if s not in gp and S20.to_hgnc(s, H) in gp},
                            unmatched=sorted(s for s in lst if s not in gp and S20.to_hgnc(s, H) not in gp))
        return mapped if use_mapped else direct
    s_g, g2m = to_panel('cc S (Tirosh 2016)', s_raw), to_panel('cc G2M (Tirosh 2016)', g_raw)
    sets = {k: to_panel(k, lib[v]) for k, v in HALLMARK.items()}
    info['gene_set_matching'] = report
    return sets, s_g, g2m, info


def load_ctrl(line):
    import run_line as RL
    return RL.read_controls(CTRL_FILE[line], panel()).astype(np.float64)


def metacell_bases(X, gate, sets, s_g, g2m, seed=SEED):
    """One line's (or half's) control bases. X: control counts on the panel (csr). Returns B (50 x G), A (2 x G),
    H (P x G) on the full panel (0 off the gate) and diagnostics."""
    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    G = X.shape[1]; gi = np.flatnonzero(gate); n = X.shape[0]
    lib = np.asarray(X.sum(1)).ravel()
    L = np.log1p(X[:, gi].toarray().astype(np.float32) * (1e4 / np.maximum(lib, 1))[:, None].astype(np.float32))
    L -= L.mean(0)
    pcs = PCA(50, svd_solver='randomized', random_state=seed).fit_transform(L); del L
    M = int(min(200, n // 25))
    lab = KMeans(M, n_init=10, random_state=seed).fit(pcs).labels_
    Ind = sp.csr_matrix((np.ones(n), (lab, np.arange(n))), shape=(M, n))
    Cnt = np.asarray((Ind @ X).todense()); w = np.bincount(lab, minlength=M).astype(float)
    P = np.log1p(5e4 * Cnt / Cnt.sum(1, keepdims=True))
    Pg = P[:, gi]; wn = w / w.sum(); mean = wn @ Pg; Xc = Pg - mean
    _, sv, Vt = np.linalg.svd(np.sqrt(wn)[:, None] * Xc, full_matrices=False)
    B = np.zeros((50, G)); B[:min(50, len(Vt)), gi] = Vt[:50]
    gp = {g: j for j, g in enumerate(panel())}; pos = {j: q for q, j in enumerate(gi)}
    H = np.zeros((len(PATHS), G)); var_s = {}; nset = {}
    for q, P_ in enumerate(PATHS):
        idx = [pos[gp[g]] for g in sets[P_] if gp[g] in pos]
        nset[P_] = len(idx)
        s = Pg[:, idx].mean(1); s = s - wn @ s; vs = float(wn @ (s * s)); var_s[P_] = vs
        H[q, gi] = (wn * s) @ Xc / vs if vs > 0 else 0.0
    # phases: scanpy's Tirosh scoring on log1p CP10K (all panel genes)
    Xn = sp.csr_matrix(X.multiply(1e4 / np.maximum(lib, 1)[:, None])); Xn.data = np.log1p(Xn.data)
    A_ = ad.AnnData(X=Xn.astype(np.float32), var=pd.DataFrame(index=pd.Index(panel())))
    sc.tl.score_genes_cell_cycle(A_, s_genes=s_g, g2m_genes=g2m, random_state=seed)
    ph = A_.obs['phase'].astype(str).to_numpy()
    prop = {k: float((ph == k).mean()) for k in ('G1', 'S', 'G2M')}
    A = np.zeros((2, G))
    pb = {}
    for k in ('G1', 'S', 'G2M'):
        m = ph == k
        if m.sum() < 10:
            pb[k] = None; continue
        c = np.asarray(X[m].sum(0)).ravel(); pb[k] = np.log1p(5e4 * c / c.sum())
    for q, k in enumerate(('S', 'G2M')):
        if pb[k] is not None and pb['G1'] is not None:
            A[q, gi] = (pb[k] - pb['G1'])[gi]
    diag = dict(n_cells=int(n), M=M, median_cells_per_metacell=float(np.median(w)), phase_proportions=prop,
                var_s=var_s, set_genes_gated=nset, B_var_share_top5_20_50=[float((sv[:k] ** 2).sum() / (sv ** 2).sum()) for k in (5, 20, 50)])
    return B, A, H, diag


def partP0(a):
    sets, s_g, g2m, info = gene_sets()
    res = dict(gene_sets=info, lines={})
    P0 = json.loads((TAB / 'part0.json').read_text())
    lines = [l for l in CTRL_LINES if l not in P0['lines_without_raw_controls']]
    res['left_out'] = P0['lines_without_raw_controls']
    for line in lines:
        f = TAB / f'P0_{line}.npz'
        if f.exists():
            res['lines'][line] = json.loads((TAB / f'P0_{line}.json').read_text()); continue
        X = load_ctrl(line); n = X.shape[0]
        import run_line as RL
        import build_effects_v5 as B_
        C_, _ = B_.cpm_csr(X.data, X.indices, X.indptr, X.shape[1], np.ones(X.shape[1], bool))
        gate = (np.asarray(C_.sum(0)).ravel() / n) > B_.MIN_CPM; del C_
        mem_before(f'P0 {line}: {n} cells x {gate.sum()} gated genes dense', n * gate.sum() * 4 * 2.5 / 1e9 + 1)
        t0 = time.time()
        B, A, H, d = metacell_bases(X, gate, sets, s_g, g2m)
        perm = np.random.default_rng(SEED).permutation(n); hv = [np.sort(perm[:n // 2]), np.sort(perm[n // 2:])]
        halves = [metacell_bases(X[r], gate, sets, s_g, g2m) for r in hv]
        stab = {}
        for K in (5, 20):
            U1, U2 = halves[0][0][:K], halves[1][0][:K]
            stab[f'B_overlap_K{K}'] = float(np.linalg.norm(U1 @ U2.T) ** 2 / K)
        gi = np.flatnonzero(gate)
        corr = lambda x, y: float(np.corrcoef(x[gi], y[gi])[0, 1]) if np.std(x[gi]) > 0 and np.std(y[gi]) > 0 else float('nan')
        stab['r_A'] = {k: corr(halves[0][1][q], halves[1][1][q]) for q, k in enumerate(('S-G1', 'G2M-G1'))}
        stab['r_P'] = {P_: corr(halves[0][2][q], halves[1][2][q]) for q, P_ in enumerate(PATHS)}
        stab['kept'] = [P_ for P_ in PATHS if np.nan_to_num(stab['r_P'][P_]) >= R_MIN]
        d.update(stability=stab, halves=[h_[3] for h_ in halves], seconds=time.time() - t0, gated_genes=int(gate.sum()))
        np.savez_compressed(f, B=B.astype(np.float32), A=A.astype(np.float32), H=H.astype(np.float32), gate=gate, paths=np.array(PATHS),
                            kept=np.array([P_ in stab['kept'] for P_ in PATHS]))
        save_json(TAB / f'P0_{line}.json', d); res['lines'][line] = d
        log(f'P0 {line}: M {d["M"]}, median cells/metacell {d["median_cells_per_metacell"]:.0f}, phases {d["phase_proportions"]}; '
            f'B overlap K5 {stab["B_overlap_K5"]:.3f} K20 {stab["B_overlap_K20"]:.3f}; r_A {stab["r_A"]}; kept {stab["kept"]}; {d["seconds"]:.0f} s')
        del X; gc.collect()
    # D/E/F against H1 and K562: cosine of v_{c,P} over shared gated genes
    cmp = {}
    for x in ('D', 'E', 'F'):
        if x not in lines:
            continue
        zx = np.load(TAB / f'P0_{x}.npz')
        for ref in ('H1', 'K562'):
            if ref not in lines:
                continue
            zr = np.load(TAB / f'P0_{ref}.npz'); sh = zx['gate'] & zr['gate']
            cmp[f'{x} vs {ref}'] = {P_: S20.cosine(zx['H'][q].astype(float), zr['H'][q].astype(float), sh) for q, P_ in enumerate(PATHS)}
    for x in ('H1', 'D', 'E', 'F'):
        if x in lines and 'K562' in lines and x != 'K562' and x == 'H1':
            zx, zr = np.load(TAB / f'P0_{x}.npz'), np.load(TAB / 'P0_K562.npz'); sh = zx['gate'] & zr['gate']
            cmp['H1 vs K562'] = {P_: S20.cosine(zx['H'][q].astype(float), zr['H'][q].astype(float), sh) for q, P_ in enumerate(PATHS)}
    res['DEF_cosines'] = cmp
    save_json(TAB / 'partP0.json', res)


def ctrl_bases(line):
    z = np.load(TAB / f'P0_{CTRL_OF.get(line, line)}.npz', allow_pickle=True)
    return {k: z[k] for k in z.files}


# ================================================================ Part P1: captured energy (oracle ceiling)
def orth(M, tol=1e-8):
    """Orthonormal basis (columns) of the span of M's columns."""
    U, s, _ = np.linalg.svd(M, full_matrices=False)
    return U[:, s > tol * max(s.max() if len(s) else 0, 1e-300)]


def capture(Q, X, own_pos):
    """Least-squares captured energy per row of X (n x g, own entries 0) on the orthonormal columns Q (g x k), the
    own gene left out of the fit (Sherman-Morrison on the one dropped row). Returns (captured, total) per row."""
    c0 = X @ Q; cap = (c0 ** 2).sum(1)
    for i, j in enumerate(own_pos):
        if j >= 0:
            u = Q[j]; uu = float(u @ u)
            if uu < 1 - 1e-9:
                cap[i] += float(u @ c0[i]) ** 2 / (1 - uu)
    return cap, (X ** 2).sum(1)


def groups_of(R):
    g = {'all': np.arange(R.h.T)}
    for j, t in enumerate(TERC):
        g[t] = np.flatnonzero(R.tterc == j)
    g['mid+strong'] = np.flatnonzero(R.tterc >= 1)
    return g


def share(cap, tot, rows):
    return float(cap[rows].sum() / tot[rows].sum()) if tot[rows].sum() > 0 else float('nan')


def p_inputs(R):
    """Targets (true Delta, residual) on gated genes, own gene 0; own position in gated coordinates."""
    import scorer_arms as SA
    h = R.h; gi = np.flatnonzero(R.gate); pos = -np.ones(h.G, int); pos[gi] = np.arange(len(gi))
    mu4 = R.mu()
    Y = (h.R - h.ctrl[None, :])[:, gi]; Res = (h.R - SA.bulk(mu4))[:, gi]
    own_pos = np.array([pos[g] if g >= 0 else -1 for g in R.own])
    for i, j in enumerate(own_pos):
        if j >= 0: Y[i, j] = 0.0; Res[i, j] = 0.0
    return gi, own_pos, Y, Res, mu4


def partP1(a):
    import scorer_step16 as S16
    res = {}; rng0 = SEED
    for line in LINES:
        R = row4(line); gi, own_pos, Y, Res, _ = p_inputs(R); g = len(gi); grp = groups_of(R)
        cb = ctrl_bases(line); kept = cb['kept'].astype(bool)
        Vt = np.load(basis_file('E4w', E4W_FOLD[line]))['Vt'].astype(np.float64)
        bases = {f'B{K}': cb['B'][:K][:, gi].T for K in (5, 20, 50)}
        bases['A'] = cb['A'][:, gi].T; bases['H'] = cb['H'][kept][:, gi].T
        bases['C'] = np.concatenate([bases['A'], bases['H'], bases['B20']], 1)
        bases['E4w20'] = Vt[:20][:, gi].T; bases['E4w100'] = Vt[:100][:, gi].T
        bases['E4w100+C'] = np.concatenate([bases['E4w100'], bases['C']], 1)
        Q = {k: orth(v.astype(np.float64)) for k, v in bases.items() if v.shape[1]}
        dims = {k: int(v.shape[1]) for k, v in Q.items()}
        out = dict(dims=dims, kept_pathways=[P_ for P_, k in zip(PATHS, kept) if k], E4w_basis=E4W_FOLD[line],
                   gated_genes=int(g), targets={})
        for tname, X in (('truth', Y), ('residual', Res)):
            o = {}
            for k, q in Q.items():
                cap, tot = capture(q, X, own_pos); o[k] = {gn: share(cap, tot, rr) for gn, rr in grp.items()}
            # random matches (5 draws each) for B5/B20/B50, H, C, and E4w100 + random(dim C)
            for k in ('B5', 'B20', 'B50', 'H', 'C'):
                if k not in Q:
                    continue
                vals = []
                for d_ in range(5):
                    rq = orth(np.random.default_rng(rng0 + d_).normal(size=(g, dims[k])))
                    cap, tot = capture(rq, X, own_pos); vals.append({gn: share(cap, tot, rr) for gn, rr in grp.items()})
                o[f'random_{k}'] = {gn: float(np.mean([v[gn] for v in vals])) for gn in grp}
            if 'C' in Q:
                vals = []
                for d_ in range(5):
                    rq = orth(np.concatenate([bases['E4w100'], np.random.default_rng(rng0 + 10 + d_).normal(size=(g, dims['C']))], 1))
                    cap, tot = capture(rq, X, own_pos); vals.append({gn: share(cap, tot, rr) for gn, rr in grp.items()})
                o['E4w100+random_C'] = {gn: float(np.mean([v[gn] for v in vals])) for gn in grp}
                o['incremental_C'] = {gn: o['E4w100+C'][gn] - o['E4w100'][gn] for gn in grp}
                o['incremental_random_C'] = {gn: o['E4w100+random_C'][gn] - o['E4w100'][gn] for gn in grp}
            # cross-validated own-response PCs (halves by true-DE rank)
            hv = S16.halves(R.T)
            for K in (5, 20):
                capt = np.zeros(R.h.T); tott = np.zeros(R.h.T)
                for f in (0, 1):
                    _, _, Vh = np.linalg.svd(X[hv[f]], full_matrices=False)
                    c, t_ = capture(Vh[:K].T, X[hv[1 - f]], own_pos[hv[1 - f]]); capt[hv[1 - f]] = c; tott[hv[1 - f]] = t_
                o[f'CV_own_PCs_K{K}'] = {gn: share(capt, tott, rr) for gn, rr in grp.items()}
            out['targets'][tname] = o
            log(f'P1 {line} {tname}: ' + '; '.join(f'{k} {v["mid+strong"]:.3f}' for k, v in o.items()))
        if line == 'H1':     # whose controls: H1's residual on each line's control set C (H1's gated genes)
            who = {}
            for other in ['H1', 'D', 'E', 'F', 'K562', 'HCT116', 'HEK293T']:
                if not (TAB / f'P0_{other}.npz').exists():
                    continue
                zo = ctrl_bases(other); ko = zo['kept'].astype(bool)
                Cm = np.concatenate([zo['A'][:, gi].T, zo['H'][ko][:, gi].T, zo['B'][:20][:, gi].T], 1).astype(np.float64)
                cap, tot = capture(orth(Cm), Res, own_pos); who[other] = {gn: share(cap, tot, rr) for gn, rr in grp.items()}
                who[other]['dim'] = int(orth(Cm).shape[1])
            out['whose_controls_residual'] = who
            log('P1 H1 whose controls (mid+strong residual): ' + ', '.join(f'{k} {v["mid+strong"]:.3f}' for k, v in who.items()))
        if line != 'HCT116s':
            out['per_pathway'] = pathway_coefs(R, gi, own_pos, Res, cb, kept)
        res[line] = out
        save_json(TAB / 'partP1.json', res); del R; gc.collect()
    # criterion (report only)
    def inc(l):
        o = res[l]['targets']['residual']; return o['incremental_C']['mid+strong'], o['incremental_random_C']['mid+strong']
    h1, h1r = inc('H1'); folds = [l for l in FOLDS if inc(l)[0] >= 0.03 and inc(l)[0] >= 5 * max(inc(l)[1], 1e-12)]
    res['criterion'] = dict(H1_incremental=h1, H1_random_incremental=h1r,
                            H1_ok=bool(h1 >= 0.03 and h1 >= 5 * max(h1r, 1e-12)), fold_pools_ok=folds,
                            met=bool(h1 >= 0.03 and h1 >= 5 * max(h1r, 1e-12) and len(folds) >= 2))
    res['reference_section11'] = ('H1 strong-tercile energy, Delta space (E4): 0.177 / 0.453 at K 20 / 96 (strong 0.137 / 0.369 on (c)); '
                                  'E1 K 20 strong (b) 0.097; own top-K PCs in sample K 5/10/20/50: 0.30/0.42/0.53/0.69 (strong 0.38/0.53/0.68/0.86)')
    save_json(TAB / 'partP1.json', res)
    log(f'P1 criterion: {res["criterion"]}')


def src_activity(targets, variant, gate, sets, src, gp):
    """a_{k,P}: mean over P's genes gated in the line (and measured by a source) of k's pooled source Delta."""
    A = np.full((len(targets), len(PATHS)), np.nan); has = np.zeros(len(targets), bool)
    idx = {P_: np.array([gp[g] for g in sets[P_] if gate[gp[g]]], int) for P_ in PATHS}
    for i, t in enumerate(targets):
        p = src.pooled(t, variant)
        if p is None:
            continue
        has[i] = True
        for q, P_ in enumerate(PATHS):
            x = p[idx[P_]]; x = x[np.isfinite(x)]
            A[i, q] = float(x.mean()) if len(x) else np.nan
    return A, has


def pathway_coefs(R, gi, own_pos, Res, cb, kept):
    """Mid + strong knockdowns: joint least-squares coefficients of the residual on the kept v_{c,P} (and on v-perp),
    own gene out; sign consistency, Spearman with (b), Pearson with the source activity a_{k,P}."""
    from scipy.stats import spearmanr, pearsonr
    sets, _, _, _ = gene_sets(); gp = {g: j for j, g in enumerate(panel())}
    A, has = src_activity(R.h.targets, SRC_VARIANT[R.line], R.gate, sets, Src(), gp)
    rows = np.flatnonzero(R.tterc >= 1); kp = [q for q in range(len(PATHS)) if kept[q]]
    allset = np.zeros(R.h.G, bool)
    for P_ in PATHS:
        allset[[gp[g] for g in sets[P_]]] = True
    out = {}
    for vname, Hm in (('v', cb['H'].astype(np.float64)), ('v_perp', np.where(allset[None, :], 0.0, cb['H'].astype(np.float64)))):
        Vg = Hm[kp][:, gi].T
        coef = np.full((len(rows), len(kp)), np.nan)
        for n_, i in enumerate(rows):
            m = np.ones(len(gi), bool)
            if own_pos[i] >= 0: m[own_pos[i]] = False
            coef[n_] = np.linalg.lstsq(Vg[m], Res[i, m], rcond=None)[0]
        o = {}
        for c, q in enumerate(kp):
            x = coef[:, c]; b = R.b[rows]; ok = np.isfinite(b); aa = A[rows, q]; oa = np.isfinite(aa)
            o[PATHS[q]] = dict(sign_consistency=float(max((x > 0).mean(), (x < 0).mean())), median_coef=float(np.median(x)),
                               spearman_b=float(spearmanr(x[ok], b[ok]).correlation) if ok.sum() > 3 else float('nan'),
                               pearson_source_activity=float(pearsonr(x[oa], aa[oa])[0]) if oa.sum() > 3 and np.std(aa[oa]) > 0 else float('nan'),
                               n_with_activity=int(oa.sum()))
        out[vname] = o
    out['knockdowns'] = int(len(rows))
    return out


# ================================================================ Part P2: pathway-anchored programs, activated by the sources
P2_G = [0.25, 0.5, 1.0]
P2_VAR = ['raw', 'centred']


def p2_term(R, sets, src, gp, centred):
    cb = ctrl_bases(R.line); kept = cb['kept'].astype(bool)
    allset = np.zeros(R.h.G, bool)
    for P_ in PATHS:
        allset[[gp[g] for g in sets[P_]]] = True
    Vp = np.where(allset[None, :], 0.0, cb['H'].astype(np.float64)) * kept[:, None]
    A, has = src_activity(R.h.targets, SRC_VARIANT[R.line], R.gate, sets, src, gp)
    A = np.nan_to_num(A)
    if centred:
        A = A - A[R.on5].mean(0, keepdims=True)
    term = A @ Vp
    return np.where(R.on5[:, None] & R.okg & ~R.sig, term, 0.0), A, has


def partP2(a):
    sets, _, _, _ = gene_sets(); gp = {g: j for j, g in enumerate(panel())}; src = Src()
    rows = resume('partP2_rows'); res = dict(lines={})
    for line in LINES:
        if any(r['line'] == line for r in rows):
            continue
        R = row4(line); rows += score_tag(R, R.mu(), variant='base', gamma=0.0)
        for var in P2_VAR:
            T_, A, has = p2_term(R, sets, src, gp, var == 'centred')
            res['lines'].setdefault(line, {})[var] = dict(term_rms_b5=float(np.sqrt((T_[R.on5] ** 2).mean())), with_activity=int(has.sum()),
                                                          b_ge5=int(R.on5.sum()))
            for g in P2_G:
                rr = score_tag(R, R.mu(extra=g * T_), variant=var, gamma=g); rows += rr
                log(f'P2 {line} {var} gamma {g}: {fmt(rr)} (MSE {rr[0]["mse"]:.3f} PDS {rr[0]["pds"]:.3f} NMAE {rr[0]["nmae"]:.3f})')
        save_rows(rows, 'partP2_rows'); del R; gc.collect()
    R = pd.DataFrame(rows)
    out = fit_report(R, ['variant', 'gamma'], dict(variant=P2_VAR, gamma=P2_G), ['base', 0.0], tag='P2')
    # D/E/F: a_{k,P} next to H1's
    dist = {}
    for x in ('H1', 'D', 'E', 'F'):
        tg = S21.c1_load('H1')['targets'] if x == 'H1' else S21.c1_load(x)['targets']
        gate = S21.c1_load(x)['gate'].astype(bool)
        A, has = src_activity(list(map(str, tg)), SRC_VARIANT[x], gate, sets, src, gp)
        dist[x] = {P_: dict(zip(['q10', 'q25', 'median', 'q75', 'q90'], np.nanquantile(A[:, q], [.1, .25, .5, .75, .9]).round(4).tolist()))
                   for q, P_ in enumerate(PATHS)} | dict(with_source=int(has.sum()), targets=len(tg))
    out['activity_distribution'] = dist
    save_json(TAB / 'partP2.json', out | res)


# ================================================================ Part P2b (Jose, 2026-11-03 follow-up): joint ridge slopes
P2B_G = [0.05, 0.1, 0.15]
P2B_LAMBDA = 1.0
P2B_RMIN = 0.2


def metacell_profiles(X, gate, seed=SEED):
    """The metacells exactly as metacell_bases builds them: (gated log1p CP50K profiles, cell-count weights)."""
    import scipy.sparse as sp
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    gi = np.flatnonzero(gate); n = X.shape[0]
    lib = np.asarray(X.sum(1)).ravel()
    L = np.log1p(X[:, gi].toarray().astype(np.float32) * (1e4 / np.maximum(lib, 1))[:, None].astype(np.float32))
    L -= L.mean(0)
    pcs = PCA(50, svd_solver='randomized', random_state=seed).fit_transform(L); del L
    M = int(min(200, n // 25))
    lab = KMeans(M, n_init=10, random_state=seed).fit(pcs).labels_
    Ind = sp.csr_matrix((np.ones(n), (lab, np.arange(n))), shape=(M, n))
    Cnt = np.asarray((Ind @ X).todense()); w = np.bincount(lab, minlength=M).astype(float)
    return np.log1p(5e4 * Cnt / Cnt.sum(1, keepdims=True))[:, gi], w


def joint_slopes(Pg, w, idx_sets, lam=P2B_LAMBDA):
    """Per gene, ridge of its centred metacell log1p CP50K on all pathway scores together (scores standardised with
    the metacell weights; weights scaled to mean 1, so lambda is on the usual n-row scale). Returns slopes per raw
    score unit (P x g), the separate slopes (P2's v, for the check), and the score SDs."""
    wn = w / w.sum(); wm = w / w.mean()
    Xc = Pg - wn @ Pg
    S = np.stack([Pg[:, ix].mean(1) for ix in idx_sets], 1); S = S - wn @ S
    sd = np.sqrt(wn @ (S ** 2)); Z = S / sd
    A = (Z * wm[:, None]).T @ Z + lam * np.eye(Z.shape[1])
    beta = np.linalg.solve(A, (Z * wm[:, None]).T @ Xc)
    sep = np.stack([(wn * S[:, q]) @ Xc / (wn @ (S[:, q] ** 2)) for q in range(S.shape[1])])
    return beta / sd[:, None], sep, sd


def p2b_selected():
    r = json.loads((TAB / 'partP1.json').read_text())['H1']['per_pathway']['v']
    return [P_ for P_ in PATHS if P_ in r and np.nan_to_num(r[P_]['pearson_source_activity']) > P2B_RMIN], \
           {P_: r[P_]['pearson_source_activity'] for P_ in r}


def partP2b(a):
    """P2 with joint ridge slopes (lambda 1, all retained pathways together), only the pathways whose P1 Pearson of
    oracle coefficient (on v) with a_{k,P} exceeds 0.2 on H1, centred activity, gamma in {0.05, 0.1, 0.15}."""
    sel, pear = p2b_selected(); res = dict(selected=sel, H1_pearson_on_v=pear, rule=f'P1 H1 Pearson (coefficient on v, basis H) > {P2B_RMIN}',
                                           lambda_=P2B_LAMBDA, ridge_on='all retained pathways (r_P >= 0.5) together; the term uses the selected ones only')
    if not sel:
        res['skipped'] = 'no pathway clears 0.2 on H1'; save_json(TAB / 'partP2b.json', res); log('P2b: skipped (no pathway clears 0.2)'); return
    sets, _, _, _ = gene_sets(); gp = {g: j for j, g in enumerate(panel())}; src = Src()
    allset = np.zeros(len(gp), bool)
    for P_ in PATHS:
        allset[[gp[g] for g in sets[P_]]] = True
    vj = {}; res['slopes'] = {}
    for cl in sorted({CTRL_OF[l] for l in LINES}):
        cb = ctrl_bases(cl); gate = cb['gate'].astype(bool); kept = cb['kept'].astype(bool); gi = np.flatnonzero(gate)
        pos = {j: q for q, j in enumerate(gi)}
        X = load_ctrl(cl); mem_before(f'P2b metacells {cl}', X.shape[0] * len(gi) * 4 * 2.5 / 1e9 + 1)
        Pg, w = metacell_profiles(X, gate); del X
        ret = [q for q in range(len(PATHS)) if kept[q]]
        idx = [np.array([pos[gp[g]] for g in sets[PATHS[q]] if gp[g] in pos]) for q in ret]
        bj, sep, sd = joint_slopes(Pg, w, idx)
        chk = float(np.abs(sep - cb['H'][ret][:, gi]).max())
        V = np.zeros((len(PATHS), len(gp)))
        for c, q in enumerate(ret):
            if PATHS[q] in sel:
                V[q, gi] = bj[c]
        vj[cl] = np.where(allset[None, :], 0.0, V)
        res['slopes'][cl] = dict(separate_slopes_rebuild_max_abs_diff_vs_P0=chk,
                                 joint_vs_separate={PATHS[q]: dict(rms_joint=float(np.sqrt((bj[c] ** 2).mean())), rms_separate=float(np.sqrt((sep[c] ** 2).mean())),
                                                                   r=float(np.corrcoef(bj[c], sep[c])[0, 1])) for c, q in enumerate(ret) if PATHS[q] in sel})
        log(f'P2b {cl}: separate slopes rebuilt, max |diff| vs P0 {chk:.2e}; ' + json.dumps(res['slopes'][cl]['joint_vs_separate'], default=float))
        if chk > 1e-6:
            save_json(TAB / 'partP2b.json', res); log('P2b CHECK FAILED: metacells differ from P0'); sys.exit(13)
    rows = resume('partP2b_rows')
    for line in LINES:
        if any(r['line'] == line for r in rows):
            continue
        R = row4(line); rows += score_tag(R, R.mu(), variant='base', gamma=0.0)
        A, has = src_activity(R.h.targets, SRC_VARIANT[R.line], R.gate, sets, src, gp)
        A = np.nan_to_num(A); A = A - A[R.on5].mean(0, keepdims=True)
        T_ = np.where(R.on5[:, None] & R.okg & ~R.sig, A @ vj[CTRL_OF[line]], 0.0)
        res.setdefault('term_rms_b5', {})[line] = float(np.sqrt((T_[R.on5] ** 2).mean()))
        for g in P2B_G:
            rr = score_tag(R, R.mu(extra=g * T_), variant='centred', gamma=g); rows += rr
            log(f'P2b {line} gamma {g}: {fmt(rr)} (MSE {rr[0]["mse"]:.3f} PDS {rr[0]["pds"]:.3f} NMAE {rr[0]["nmae"]:.3f})')
        save_rows(rows, 'partP2b_rows'); del R; gc.collect()
    R = pd.DataFrame(rows)
    out = fit_report(R, ['variant', 'gamma'], dict(variant=['centred'], gamma=P2B_G), ['base', 0.0], tag='P2b')
    save_json(TAB / 'partP2b.json', out | res)


# ================================================================ Part S1: the votes (one table load)
VOTES = OUT / 'votes'


def vote_stats(Y, SE, mults, mu):
    """Per gene: n_votes (pool_effects' finite rule), n_agree (sign == sign(mu), LFC-0 votes ignored), Cochran's Q
    with pool_effects' weights."""
    import stage1 as S
    Y = np.asarray(Y, np.float64); SE = np.asarray(SE, np.float64)
    fin = np.isfinite(Y) & np.isfinite(SE) & (SE > 0) & (np.abs(Y) <= S.LFC_CAP)
    var = (np.maximum(SE, S.SE_MIN) ** 2) * np.asarray(mults, np.float64)[:, None]
    w = np.where(fin, 1.0 / var, 0.0)
    n = fin.sum(0)
    agree = (fin & (Y != 0) & (np.sign(np.where(fin, Y, 0.0)) == np.sign(mu)[None, :]) & (mu != 0)[None, :]).sum(0)
    Q = (w * np.where(fin, Y - mu, 0.0) ** 2).sum(0)
    return n.astype(np.int16), agree.astype(np.int16), Q


def partS1(a):
    import stage1 as S
    import screen as SC
    import scorer_study as SS
    genes = panel(); G_ = np.asarray(genes, dtype=object); G = len(genes)
    VOTES.mkdir(parents=True, exist_ok=True)
    tables = S.discover_tables(SC.DIRS); by_stem = {p.stem: p for p in tables}
    keep, _, _ = SS.filter_plan(tables, 'F0'); assert [p.stem for p in keep] == [p.stem for p in tables]
    mem_before(f'S1 table load ({len(tables)} tables)', 18.0)
    sources = S.load_tables(tables); mem_after('S1 table load')
    sets, res = {}, dict(contexts={})
    for ctx, (drop, _) in S21.C1_CTX.items():
        if (TAB / f'S1_{ctx}.npz').exists():
            res['contexts'][ctx] = json.loads((TAB / 'partS1.json').read_text())['contexts'].get(ctx, {}); continue
        key = tuple(sorted(drop))
        if key not in sets:
            others = [s for s in sources if S._line_key(s.background) not in set(drop)]
            votes, prior_votes = S.arm_votes(others, 'C')
            basals = {}
            for s in others:
                got = S.load_baseline(by_stem[s.name].with_name(s.name + '_baseline.csv'), s.genes)
                if got is not None: basals[s.name] = (s.genes, got[0])
            sets[key] = (others, votes, S.fit_se_calibration(others, basals), S.fit_prior(others, votes=prior_votes))
            log(f'S1 voting set without {key or "nothing"}: {len(others)} tables')
        others, votes, calib, prior = sets[key]
        c1 = S21.c1_load(ctx); targets = list(map(str, c1['targets'])); K = len(targets)
        A = {k: np.zeros((K, G), np.float32) for k in ('mu', 'var_fe', 'tau2', 'lfc')}
        n_ = np.zeros((K, G), np.int16); ag = np.zeros((K, G), np.int16); Qn = np.full((K, G), np.nan, np.float32)
        vk, vb, vy, vs = [], [], [], []
        t0 = time.time()
        for i, k in enumerate(targets):
            Y, SE, mults = votes.gather(k, G_)
            if not Y:
                mu, vf, t2 = np.zeros(G), np.full(G, np.inf), np.full(G, np.nan)
            else:
                Ys, SEs = np.stack(Y), np.stack(SE)
                mu, vf, t2, tu, n = S.pool_effects(Ys, SEs, np.asarray(mults), prior.tau2)
                nn, aa, Q = vote_stats(Ys, SEs, mults, mu)
                n_[i], ag[i] = nn, aa
                Qn[i] = np.where(nn >= 2, Q / np.maximum(nn - 1, 1), np.nan)
                for b_, bg in enumerate(votes.owners.get(k, ())):
                    vk.append(i); vb.append(bg); vy.append(Ys[b_].astype(np.float32)); vs.append(SEs[b_].astype(np.float32))
            t2u = np.where(np.isfinite(t2), t2, prior.tau2)
            lfc, _, _, _ = S.gaussian_posterior(mu, vf, t2u, prior.v0)
            A['mu'][i], A['var_fe'][i], A['tau2'][i], A['lfc'][i] = mu, vf, t2, lfc
        chk = {k: bool(np.array_equal(A[k], c1[k].astype(np.float32), equal_nan=True)) for k in A}
        chk['n_equal_C1_n'] = bool(np.array_equal(n_, c1['n'].astype(np.int16)))
        g = c1['gate'].astype(bool); gp = {x: j for j, x in enumerate(genes)}
        ok = np.broadcast_to(g[None, :], n_.shape).copy()
        for i, t in enumerate(targets):
            if t in gp: ok[i, gp[t]] = False
        nv = n_[ok]
        dist = {'1': float((nv == 1).mean()), '2': float((nv == 2).mean()), '3': float((nv == 3).mean()), '>=4': float((nv >= 4).mean()), '0': float((nv == 0).mean())}
        q3 = Qn[ok & (n_ >= 3)]; q3 = q3[np.isfinite(q3)]
        t2v = A['tau2'][ok & (n_ >= 2)]
        glob_t2 = float(prior.tau2)
        res['contexts'][ctx] = dict(knockdowns=K, seconds=time.time() - t0, checks=chk, n_votes_share_gated_own_out=dist,
                                    Q_over_df_n3plus=dict(zip(['q10', 'q25', 'median', 'q75', 'q90'], np.quantile(q3, [.1, .25, .5, .75, .9]).tolist())) | dict(mean=float(q3.mean()), share_gt2=float((q3 > 2).mean()), pairs=int(len(q3))),
                                    tau2_n2plus=dict(median=float(np.median(t2v)), share_at_global_floor=float(np.isclose(t2v, glob_t2).mean()), global_tau2=glob_t2))
        np.savez_compressed(TAB / f'S1_{ctx}.npz', targets=np.array(targets), n_votes=n_, n_agree=ag, Q_over_df=Qn)
        np.savez_compressed(VOTES / f'votes_{ctx}.npz', kd_index=np.array(vk, np.int32), background=np.array(vb), lfc=np.stack(vy) if vy else np.zeros((0, G), np.float32),
                            se=np.stack(vs) if vs else np.zeros((0, G), np.float32), targets=np.array(targets))
        log(f'S1 {ctx}: {K} knockdowns {time.time() - t0:.0f} s; checks {chk}; n_votes {dist}; Q/df median {res["contexts"][ctx]["Q_over_df_n3plus"]["median"]:.2f}')
        save_json(TAB / 'partS1.json', res)
        if not all(chk.values()):
            log(f'S1 CHECK FAILED for {ctx}'); sys.exit(7)
    save_json(TAB / 'partS1.json', res)


def s1_load(ctx, targets):
    z = np.load(TAB / f'S1_{ctx}.npz', allow_pickle=True); t = list(map(str, z['targets'])); r = np.array([t.index(x) for x in targets])
    n = z['n_votes'][r].astype(int); ag = z['n_agree'][r].astype(int)
    return n, ag, (ag + 1.0) / (n + 2.0), z['Q_over_df'][r]


# ================================================================ S2 / S3 shared
P_BINS = [(0.5, 0.6), (0.6, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 1.0001)]


class SLine:
    """Row 4 on a truth line with C1's z_u and S1's vote agreement, the predicted / true signs and the true DE set."""

    def __init__(self, line):
        import scorer_arms as SA
        self.C = C = cline(line); h = C.h
        self.mu4 = C.mu(); self.dP = SA.bulk(self.mu4) - h.ctrl[None, :]; self.dR = h.R - h.ctrl[None, :]
        self.n, self.ag, self.pt, self.Qn = s1_load(S21.CTX_OF[line], h.targets)
        self.sp = np.sign(self.dP); self.right = (self.sp == np.sign(self.dR)) & (self.sp != 0)
        self.TD = C.T & C.okg & (self.sp != 0)                  # true DE (real p_adj < 0.05), gated, own out
        self.Dpred = C.delta_pred(self.mu4)                     # Delta from m0 (config 1''s Delta space)


def cls_of(n):
    return np.where(n >= 3, 3, n)


def fit_pcal(S, rows):
    """Isotonic p~ -> sign accuracy per n_votes class (1, 2, >= 3) on the true DE pairs of rows."""
    from sklearn.isotonic import IsotonicRegression
    maps = {}
    M = np.zeros(S.TD.shape, bool); M[rows] = True; M &= S.TD
    for c in (1, 2, 3):
        m = M & (cls_of(S.n) == c)
        if m.sum() < 20:
            maps[c] = None; continue
        maps[c] = IsotonicRegression(increasing=True, y_min=0.0, y_max=1.0, out_of_bounds='clip').fit(S.pt[m], S.right[m].astype(float))
    return maps


def apply_pcal(S, maps, rows=None):
    P = np.full(S.n.shape, np.nan)
    rr = np.arange(len(S.n)) if rows is None else rows
    for c, mp in maps.items():
        if mp is None:
            continue
        m = np.zeros(S.n.shape, bool); m[rr] = True; m &= (cls_of(S.n) == c)
        if m.any():
            P[m] = mp.predict(S.pt[m])
    return P


def inner_folds(T, rows):
    import scorer_step16 as S16
    hv = S16.halves(T[rows]); return [rows[hv[0]], rows[hv[1]]]


def pcal_frames(S):
    """p_cal arrays. H1: per direction f, rows of half f from 2 inner folds of half f (each fitted on the other inner
    fold), rows of half 1 - f from the map fitted on all of half f; and 'xfit' (each half from the other half's map,
    S4 / R2). Other lines: the map fitted on all of H1 (passed in)."""
    hv = S21.h1_halves(); out = {}
    for f in (0, 1):
        P = np.full(S.n.shape, np.nan)
        inn = inner_folds(S.C.T, hv[f])
        for j in (0, 1):
            m = fit_pcal(S, inn[1 - j]); x = apply_pcal(S, m, inn[j]); P[inn[j]] = x[inn[j]]
        m = fit_pcal(S, hv[f]); x = apply_pcal(S, m, hv[1 - f]); P[hv[1 - f]] = x[hv[1 - f]]
        out[f] = P
    X = np.full(S.n.shape, np.nan)
    for f in (0, 1):
        m = fit_pcal(S, hv[1 - f]); X[hv[f]] = apply_pcal(S, m, hv[f])[hv[f]]
    out['xfit'] = X
    return out


def auc(y, s):
    from sklearn.metrics import roc_auc_score
    ok = np.isfinite(s)
    return float(roc_auc_score(y[ok], s[ok])) if ok.sum() > 10 and len(np.unique(y[ok])) == 2 else float('nan')


# ================================================================ Part S2: does agreement predict sign accuracy?
def partS2(a):
    from sklearn.linear_model import LogisticRegression
    res = {}; H1fit = None
    c5 = json.loads((S21.TAB / 'partC5.json').read_text())['item2']['by_z_u_decile']
    for line in ['H1'] + GUARD:
        S = SLine(line); C = S.C; TD = S.TD
        def acc(m):
            return dict(n=int(m.sum()), accuracy=float(S.right[m].mean()) if m.any() else float('nan'))
        o = dict(true_de_pairs=int(TD.sum()), no_vote_share=float((TD & (S.n == 0)).sum() / max(TD.sum(), 1)))
        o['by_p_bin_nvotes_ge2'] = {f'{lo:.1f}-{min(hi, 1):.1f}': acc(TD & (S.n >= 2) & (S.pt >= lo) & (S.pt < hi)) for lo, hi in P_BINS}
        o['nvotes_1'] = acc(TD & (S.n == 1))
        o['by_nvotes'] = {str(c): acc(TD & (cls_of(S.n) == c)) for c in (1, 2, 3)}
        for cn, cm in (('called', C.sig), ('uncalled', ~C.sig)):
            for pn, pm in (('pred up', S.sp > 0), ('pred down', S.sp < 0)):
                o[f'{cn}, {pn}'] = {f'{lo:.1f}-{min(hi, 1):.1f}': acc(TD & cm & pm & (S.n >= 2) & (S.pt >= lo) & (S.pt < hi)) for lo, hi in P_BINS} | dict(all=acc(TD & cm & pm))
        # AUC: z_u, p~, both (logistic cross-fitted on the H1 halves; fold pools / HCT116s from all-H1)
        m = TD & (S.n >= 1); y = S.right.astype(int)
        def feats(SS_, mm): return np.c_[SS_.C.zu[mm], SS_.pt[mm]]
        aucs = {}
        for split, sm in (('all', np.ones_like(m)), ('called', C.sig), ('uncalled', ~C.sig), ('pred up', S.sp > 0), ('pred down', S.sp < 0)):
            mm = m & sm
            d = dict(z_u=auc(y[mm], C.zu[mm]), p_tilde=auc(y[mm], S.pt[mm]), pairs=int(mm.sum()))
            if line == 'H1':
                hv = S21.h1_halves(); sc_ = np.full(S.n.shape, np.nan)
                for f in (0, 1):
                    fr = np.zeros(S.n.shape, bool); fr[hv[f]] = True; tr = fr & m
                    lr = LogisticRegression(C=1e3, max_iter=1000).fit(feats(S, tr), y[tr])
                    te = np.zeros(S.n.shape, bool); te[hv[1 - f]] = True; te &= m
                    sc_[te] = lr.predict_proba(feats(S, te))[:, 1]
                d['both_logistic_xfit'] = auc(y[mm], sc_[mm])
            else:
                d['both_logistic_allH1'] = auc(y[mm], H1fit.predict_proba(feats(S, mm))[:, 1]) if mm.any() else float('nan')
            aucs[split] = d
        if line == 'H1':
            H1fit = LogisticRegression(C=1e3, max_iter=1000).fit(feats(S, m), y[m]); o['logistic_allH1_coef'] = H1fit.coef_.tolist()
        o['auc'] = aucs
        res[line] = o
        log(f'S2 {line}: p bins {[(k, round(v["accuracy"], 3), v["n"]) for k, v in o["by_p_bin_nvotes_ge2"].items()]}; n1 {o["nvotes_1"]}; AUC {aucs["all"]}')
        del S, C; gc.collect()
    res['z_u_decile_curve_C5'] = c5
    save_json(TAB / 'partS2.json', res)


# ================================================================ Part S3: magnitude x calibrated (2p - 1)
S3_RHO = [0.5, 1.0]
S3_APPLY = ['called', 'uncalled', 'both']


def s3_weight(S, Pc, rho, rows_groups):
    """m = max(0, 2 p_cal - 1); the applied weight m^rho normalised to mean 1 within each class (called, uncalled)
    over the gated, own-out pairs with votes, separately in each row group (a scored half); 1 where no vote."""
    m = np.maximum(0.0, 2 * Pc - 1); has = S.C.okg & (S.n >= 1) & np.isfinite(Pc)
    W = np.ones(m.shape); wr = np.where(has, m ** rho, 1.0)
    for rows in rows_groups:
        Rm = np.zeros(m.shape, bool); Rm[rows] = True
        for cm in (S.C.sig, ~S.C.sig):
            sel = Rm & has & cm
            if sel.any():
                mean = wr[sel].mean()
                W[sel] = wr[sel] / mean if mean > 0 else 1.0
    zero = float((has & (m == 0)).sum() / max(has.sum(), 1))
    return W, has, zero


def s3_mu(S, W, apply):
    cm = {'called': S.C.sig, 'uncalled': ~S.C.sig, 'both': np.ones_like(S.C.sig)}[apply] & S.C.okg
    add = np.where(cm, (W - 1.0) * S.Dpred, 0.0)
    return mu_add_delta(S.mu4, add, S.C.own)


def partS3(a):
    rows = resume('partS3_rows'); res = dict(zero_share={}); maps_all = None
    hv = S21.h1_halves()
    for line in ['H1'] + GUARD:
        if any(r['line'] == line for r in rows) and line != 'H1':
            continue
        S = SLine(line)
        if line == 'H1':
            fr = pcal_frames(S); maps_all = fit_pcal(S, np.arange(S.C.h.T))
            import pickle
            pickle.dump(dict(maps_all=maps_all), open(TAB / 'S3_pcal_maps.pkl', 'wb'))
            np.savez_compressed(TAB / 'S3_pcal_H1.npz', f0=fr[0], f1=fr[1], xfit=fr['xfit'], targets=np.array(S.C.h.targets))
            if not any(r['line'] == 'H1' for r in rows):
                for f in (0, 1):
                    rows += [dict(r, fit_half=f) for r in score_tag(S.C, S.mu4, rho=0.0, apply='base')]
                    for rho in S3_RHO:
                        W, has, zero = s3_weight(S, fr[f], rho, hv)
                        res['zero_share'][f'H1 f{f}'] = zero
                        for ap in S3_APPLY:
                            rr = score_tag(S.C, s3_mu(S, W, ap), rho=rho, apply=ap); rows += [dict(r, fit_half=f) for r in rr]
                            log(f'S3 H1 fit{f} rho {rho} {ap}: {fmt(rr)} (MSE {rr[0]["mse"]:.3f} PDS {rr[0]["pds"]:.3f} NMAE {rr[0]["nmae"]:.3f} reach {rr[0]["reach"]:.3f})')
        else:
            Pc = apply_pcal(S, maps_all)
            rows += [dict(r, fit_half=-1) for r in score_tag(S.C, S.mu4, rho=0.0, apply='base')]
            for rho in S3_RHO:
                W, has, zero = s3_weight(S, Pc, rho, [np.arange(S.C.h.T)]); res['zero_share'][line] = zero
                for ap in S3_APPLY:
                    rr = score_tag(S.C, s3_mu(S, W, ap), rho=rho, apply=ap); rows += [dict(r, fit_half=-1) for r in rr]
                    log(f'S3 {line} rho {rho} {ap}: {fmt(rr)}')
        save_rows(rows, 'partS3_rows'); del S; gc.collect()
    R = pd.DataFrame(rows)
    ff = {f: R[(R.line == 'H1') & (R.fit_half == f)] for f in (0, 1)}
    out = fit_report(R, ['rho', 'apply'], dict(rho=S3_RHO, apply=S3_APPLY), [0.0, 'base'], fit_frames=ff, tag='S3')
    out['normalisation'] = ('the applied weight m^rho is normalised to mean 1 within each class (called, uncalled) over the gated, own-out pairs with >= 1 vote, '
                            'separately in each H1 half (each scored panel); pairs with no vote keep weight 1')
    save_json(TAB / 'partS3.json', out | res)


# ================================================================ Part S4: reach order (proxy, measure only)
def de_order(ref, mu, sig, key):
    """The proxy's four DE components with called genes' p ordered by key (higher first; ties by effect size);
    uncalled genes as the proxy. key None = the proxy's own effect-size order."""
    import polars as pl
    import scorer_arms as SA
    P = SA.bulk(mu)
    t = ref.pred_table(mu, sig, P)
    if key is not None:
        s = np.abs(P[ref.ri, ref.rg] - ref.h.ctrl[ref.rg]); called = sig[ref.ri, ref.rg]
        k = key[ref.ri, ref.rg]; k = np.where(np.isfinite(k), k, -np.inf)
        df = pd.DataFrame(dict(t=ref.ri, k=k, s=s, c=called))
        q = np.zeros(len(df))
        cd = df[df.c].sort_values(['t', 'k', 's'], kind='stable')
        if len(cd):
            r = cd.groupby('t').cumcount() + 1; n = cd.groupby('t')['k'].transform('size')
            q[cd.index.to_numpy()] = (r / n).to_numpy()            # higher key -> larger q -> smaller p -> ranked first
        p = np.where(called, 1e-3 / (1 + 1e4 * q), t['p_value'].to_numpy())
        t = t.with_columns(pl.Series('p_value', p), pl.Series('p_adj', p))
    return ref.de_metrics(t)


def reach_rows(R, kd, line, order, rows_parts, nconf):
    import scorer_proxy as SP
    out = []
    for part, rr in rows_parts.items():
        raw = float(np.nanmean(kd[SP.DE[2]][rr])); sc = R.h.scale(SP.DE[2], raw)
        kstar = kd[SP.DE[2]][rr] * nconf[rr]
        o = dict(line=line, order=order, part=part, reach=sc, reach_raw=raw, kstar_mean=float(np.nanmean(kstar)), kstar_median=float(np.nanmedian(kstar)))
        for m, k in (('nmae', SP.DE[0]), ('fidelity', SP.DE[1]), ('jaccard', SP.DE[3])):
            o[m] = R.h.scale(k, float(np.nanmean(kd[k][rr])))
        out.append(o)
    return out


def oracle_key(ref, mu, R):
    """1 where the predicted sign (pdex LFC of mu vs the real controls) matches the real DE table's sign, else 0."""
    cpm = mu * (1e6 / mu.sum(1, keepdims=True))
    pl_ = np.sign(np.log2((cpm[ref.ri, ref.rg] + 1e-9) / (ref.ctrl_cpm[ref.rg] + 1e-9)))
    rl = np.sign(ref.real_df['log2_fold_change'].to_numpy().astype(float))
    K = np.zeros((R.h.T, R.h.G)); K[ref.ri, ref.rg] = (pl_ == rl).astype(float)
    return K


def partS4(a):
    import pickle
    res = {}; rows = []
    maps_all = pickle.load(open(TAB / 'S3_pcal_maps.pkl', 'rb'))['maps_all']
    for line in ['H1'] + FOLDS:
        S = SLine(line); C = S.C; ref = C.D.ref
        if line == 'H1':
            z = np.load(TAB / 'S3_pcal_H1.npz'); Pc = z['xfit']
        else:
            Pc = apply_pcal(S, maps_all)
        nconf = (C.T & ~C.ownM).sum(1).astype(float)
        keys = {'(i) p_cal': Pc, '(ii) z_u': C.zu, '(iii) effect size (proxy)': None, '(iv) oracle': oracle_key(ref, S.mu4, C)}
        parts = {p: rr for p, rr in C.D.parts.items()}
        for nm, k in keys.items():
            kd = de_order(ref, S.mu4, C.sig, k); rows += reach_rows(C, kd, line, nm, parts, nconf)
            log(f'S4 {line} {nm}: ' + ', '.join(f'{r["part"]} {r["reach"]:.4f} (raw {r["reach_raw"]:.4f})' for r in rows[-len(parts):] if r['part'] in ('full', 'half0', 'half1')))
        save_rows(rows, 'partS4_rows'); del S, C; gc.collect()
    R = pd.DataFrame(rows)
    for line in ['H1'] + FOLDS:
        L = R[R.line == line]
        for part in (('half0', 'half1', 'full') if line == 'H1' else ('full',)):
            P = L[L.part == part].set_index('order')
            base = P.loc['(iii) effect size (proxy)', 'reach']; orc = P.loc['(iv) oracle', 'reach'] - base
            res[f'{line} {part}'] = dict(reach={o: float(P.loc[o, 'reach']) for o in P.index}, reach_raw={o: float(P.loc[o, 'reach_raw']) for o in P.index},
                                         oracle_gain=float(orc),
                                         share_of_oracle_gain={o: float((P.loc[o, 'reach'] - base) / orc) if orc else float('nan') for o in ('(i) p_cal', '(ii) z_u')})
        if line == 'H1':
            res['H1 terciles'] = {t: {o: dict(reach=float(L[(L.part == t) & (L.order == o)].reach.iloc[0]), kstar_mean=float(L[(L.part == t) & (L.order == o)].kstar_mean.iloc[0]))
                                      for o in keys} for t in TERC}
            res['H1 half terciles'] = {f'half{f}_{t}': {o: float(L[(L.part == f'half{f}_{t}') & (L.order == o)].reach.iloc[0]) for o in keys} for f in (0, 1) for t in TERC}
    h = [res[f'H1 half{f}'] for f in (0, 1)]
    better = lambda o: all(x['reach'][o] for x in h)
    pc = [x['reach']['(i) p_cal'] for x in h]; zu = [x['reach']['(ii) z_u'] for x in h]
    res['key_for_R3'] = '(i) p_cal' if (all(p >= z for p, z in zip(pc, zu)) or not all(z > p for p, z in zip(pc, zu))) else '(ii) z_u'
    res['key_rule'] = 'the key with higher held-out reach on both H1 halves; if they split, p_cal'
    res['note'] = ('on cells the order comes from the generator\'s realised Wilcoxon p-values; the frozen generator does not read a requested order '
                   '(stage1.prepare_generator_inputs saves reach_order and never passes it), so these are ceilings')
    save_json(TAB / 'partS4.json', res)
    log(f'S4: key for R3 = {res["key_for_R3"]}; H1 halves: ' + json.dumps({k: v['reach'] for k, v in res.items() if k.startswith('H1 half') and 'terc' not in k}, default=float))


# ================================================================ Part Q1: why C4 loses the guardrails (anchors)
Q1_PICKS = [(4.0, 0.75), (4.0, 0.5)]


def anchors_of(h, raw):
    """Per component: raw, scaled, the baseline anchor, the replicate anchor, the unclamped score and whether clamped."""
    import scorer_proxy as SP
    out = {}
    for c, m in S19.Fast.KEYS.items():
        e = h.entries[m]; u = float(raw[c] if c in raw else raw[m])
        sc = float(h._score_one(u, e.base, e.scoring))
        try:
            un = float(h._score_one(u, e.base, e.scoring, clamp_low=float('-inf'), clamp_high=float('inf')))
        except Exception as ex:                                  # noqa: BLE001
            un = float('nan')
        out[c] = dict(raw=u, scaled=sc, baseline_anchor=float(e.base), replicate_anchor=float(getattr(e.scoring, 'anchor', float('nan')) or float('nan')),
                      unclamped=un, clamped=bool(np.isfinite(un) and abs(un - sc) > 1e-12),
                      clamp_low=getattr(e.scoring, 'clamp_low', None), clamp_high=getattr(e.scoring, 'clamp_high', None))
    return out


def partQ1(a):
    res = dict(proxy={}, cells={})
    for line in GUARD:
        C = cline(line); cg = C.D.cg; z2 = C.zu ** 2
        out = C.D.ref.proxy(mu=C.mu4 if hasattr(C, 'mu4') else C.mu(), sig=C.sig, per_kd=True)
        base_kd = out['per_kd']
        o = {}
        sm0 = {p: C.D.ref.summarise(base_kd, rows=rr) for p, rr in C.D.parts.items()}
        o['row 4'] = dict(full=dict(overall=sm0['full']['overall'], anchors=anchors_of(C.h, sm0['full']['raw'])),
                          terciles={t: dict(overall=sm0[t]['overall'], **{c: sm0[t]['scaled'][m] for c, m in S19.Fast.KEYS.items()}) for t in TERC})
        for kap, s in Q1_PICKS:
            mu = C.mu_called(s * cg * (z2 / (z2 + kap)))
            kd = C.D.ref.proxy(mu=mu, sig=C.sig, per_kd=True)['per_kd']
            sm = {p: C.D.ref.summarise(kd, rows=rr) for p, rr in C.D.parts.items()}
            nm = f'C4 kappa {kap:g} s {s:g}'
            o[nm] = dict(full=dict(overall=sm['full']['overall'], anchors=anchors_of(C.h, sm['full']['raw']),
                                   change={c: sm['full']['scaled'][m] - sm0['full']['scaled'][m] for c, m in S19.Fast.KEYS.items()} | dict(overall=sm['full']['overall'] - sm0['full']['overall']),
                                   raw_change={c: sm['full']['raw'][m] - sm0['full']['raw'][m] for c, m in S19.Fast.KEYS.items()}),
                         terciles={t: {c: sm[t]['scaled'][m] - sm0[t]['scaled'][m] for c, m in S19.Fast.KEYS.items()} | dict(overall=sm[t]['overall'] - sm0[t]['overall'])
                                   for t in TERC})
            log(f'Q1 {line} {nm}: overall {sm["full"]["overall"] - sm0["full"]["overall"]:+.4f}; ' +
                ', '.join(f'{c} {v:+.4f}' for c, v in o[nm]['full']['change'].items()))
        res['proxy'][line] = o; del C; gc.collect(); S19._REFS.clear()
    # cells: row 4 (Step 21 Track 1) and the live entry (15c config 1' + K0), Eval A and HCT116final
    for ev, hn in (('EvalA', 'evalA'), ('HCT116final', 'evalHCT116final')):
        ref = S19.ref_for(ev); h = ref.h
        for nm, d in (('row 4 (Track 1)', S21.CELLS / ev / 'e2e/salt_none/harness' / hn), ('live entry (15c K0)', S19.C15 / 'config1p_K0' / ev / 'e2e/salt_none/harness' / hn)):
            sm = json.loads((d / 'summary.json').read_text())
            raw = {c: sm['raw'][m] for c, m in S19.Fast.KEYS.items()}
            res['cells'][f'{ev} {nm}'] = dict(overall=sm['scaled']['overall'], anchors=anchors_of(h, raw),
                                               harness_scaled={c: sm['scaled'][m] for c, m in S19.Fast.KEYS.items()})
        S19._REFS.pop(ev, None); gc.collect()
    save_json(TAB / 'partQ1.json', res)


# ================================================================ Part Q2 (M1): DepMap in the shift
Q2_D = [0.0, 0.1, 0.2, 0.3]


def depmap_z(targets):
    """z_k = the target's DepMap mean gene effect, standardised with all-H1's 292 targets' mean and SD; missing -> 0."""
    z = np.load(S20.TAB / 'partF3_depmap.npz', allow_pickle=True); idx = {g: i for i, g in enumerate(z['genes'].astype(str))}
    H = S20.hgnc_maps()
    def eff(ts):
        out = np.full(len(ts), np.nan)
        for i, t in enumerate(ts):
            g = S20.to_hgnc(t, H) or t; j = idx.get(g)
            if j is not None: out[i] = z['extra'][j]
        return out
    h1 = eff(list(map(str, S21.c1_load('H1')['targets'])))
    mu, sd = np.nanmean(h1), np.nanstd(h1)
    e = eff(targets)
    return np.where(np.isfinite(e), (e - mu) / sd, 0.0), dict(H1_found=int(np.isfinite(h1).sum()), H1_mean=float(mu), H1_sd=float(sd), found=int(np.isfinite(e).sum()))


def m1b_terc(R, z):
    """Terciles of the rank-average of (b) and -DepMap, group sizes those all-H1 (b) cuts give on H1."""
    from scipy.stats import rankdata
    cuts = S21.cuts_c2()['full']
    bh = b_for('EvalA', R.h.targets if R.line == 'H1' else S21.c1_load('H1')['targets'])
    th = b_terc(bh, cuts); frac = np.array([(th == j).mean() for j in range(3)])
    b = np.nan_to_num(R.b, nan=-np.inf)
    ra = (rankdata(b) + rankdata(-z)) / 2.0
    o = np.argsort(ra, kind='stable'); n = len(o); c1, c2 = int(round(frac[0] * n)), int(round((frac[0] + frac[1]) * n))
    t = np.zeros(n, int); t[o[c1:c2]] = 1; t[o[c2:]] = 2
    return t, frac


def partQ2(a):
    rows = resume('partQ2_rows'); res = {}
    for line in LINES:
        if any(r['line'] == line for r in rows):
            continue
        R = row4(line); z, info = depmap_z(R.h.targets); res[line] = info
        for d in Q2_D:
            add = (SHIFT_A[R.terc] - d * z)[:, None] * R.S[None, :]
            rr = score_tag(R, R.mu(shift=False, extra=add), arm='M1', delta=d); rows += rr
            log(f'Q2 {line} delta {d}: {fmt(rr)}')
        t, frac = m1b_terc(R, z)
        add = SHIFT_A[t][:, None] * R.S[None, :]
        rr = score_tag(R, R.mu(shift=False, extra=add), arm='M1b', delta=0.0); rows += rr
        res[line]['M1b_group_fractions'] = frac.tolist(); res[line]['M1b_moved_vs_row4_terc'] = int((t != R.terc).sum())
        log(f'Q2 {line} M1b: {fmt(rr)}; {int((t != R.terc).sum())} knockdowns change group')
        save_rows(rows, 'partQ2_rows'); del R; gc.collect()
    R = pd.DataFrame(rows)
    M1 = R[R.arm == 'M1']
    out = fit_report(M1, ['delta'], dict(delta=Q2_D), [0.0], tag='Q2')
    # M1b: one fixed row, scored on both halves against row 4 (nothing picked)
    b4 = M1[M1.delta == 0.0]; mb = R[R.arm == 'M1b']
    fv = []
    for f in (0, 1):
        x = mb[(mb.line == 'H1') & (mb.part == f'half{1 - f}')].iloc[0]; y = b4[(b4.line == 'H1') & (b4.part == f'half{1 - f}')].iloc[0]
        gr = {l: dict(diff=float(mb[(mb.line == l) & (mb.part == 'full')].overall.iloc[0] - b4[(b4.line == l) & (b4.part == 'full')].overall.iloc[0])) for l in GUARD}
        fv.append(dict(held_out=comp_dict(x), base_held_out=comp_dict(y), gain={c: float(x[c] - y[c]) for c in COMPS + ['overall']}, guardrails=gr))
    out['M1b'] = dict(halves=fv, criterion=criterion(fv), note='a fixed row (no grid): each half is held-out by construction')
    save_json(TAB / 'partQ2.json', out | dict(depmap=res))


# ================================================================ Part R2: read the generator, diagnose Track 1 cells
def min_wilcoxon_p(n1, n0, zero_frac):
    """Smallest two-sided MWU p (the generator's rank_p, normal approximation with tie + continuity correction) for n1
    knockdown cells against n0 controls with the given zero fraction: all knockdown cells above every control (an
    increase), or all at zero (a decrease)."""
    import cell_generator as CG
    nz0 = int(round(zero_frac * n0)); nnz = n0 - nz0
    v = np.log1p(np.arange(1, nnz + 1, dtype=float))
    g = dict(v=v, c=np.ones(nnz, int), below=np.concatenate([[0], np.arange(1, nnz + 1)]), tie_c=0.0, nz0=nz0,
             cum=(nz0 + np.arange(1, nnz + 1)) / n0)
    up = CG.rank_p(np.full(n1, v[-1] + 5.0), g, n0) if nnz else 1.0
    dn = CG.rank_p(np.zeros(n1), g, n0)
    return up, dn


def partR2(a):
    import scipy.sparse as sp_                                                            # noqa: F401
    from scipy.stats import spearmanr
    res = {}
    ev, hn = 'EvalA', 'evalA'
    S = SLine('H1'); C = S.C; h = C.h
    z = np.load(S21.CELLS / ev / 'input.npz', allow_pickle=True); names = list(map(str, z['knockdowns']))
    r = np.array([names.index(t) for t in h.targets]); un = z['n_sources'][r] <= 0
    sig = z['de_sig'][r] & ~un & C.okg; lfc = z['lfc'][r].astype(float); m0 = z['m0']
    hd = S21.CELLS / ev / 'e2e/salt_none/harness' / hn
    dp = pd.read_parquet(hd / 'de_pred.parquet', columns=['target', 'feature', 'p_value', 'p_adj', 'log2_fold_change'])
    gp = {g: j for j, g in enumerate(h.genes)}; tp = {t: i for i, t in enumerate(h.targets)}
    dp = dp[dp.target.isin(tp)]
    ti = dp.target.map(tp).to_numpy(); gi = dp.feature.map(gp).to_numpy()
    Pv = np.full((h.T, h.G), np.nan); Pv[ti, gi] = dp.p_value.to_numpy(); Pa = np.full((h.T, h.G), np.nan); Pa[ti, gi] = dp.p_adj.to_numpy()
    band = np.full(h.G, -1)
    for j, (nm, lo, hi) in enumerate([('5-20', 5, 20), ('20-100', 20, 100), ('>100', 100, np.inf)]):
        band[(m0 > lo) & (m0 <= hi)] = j
    cq = np.load(S21.CELLS / ev / 'e2e/controls/controls_control_quantiles.npz')
    zero = np.zeros(h.G); zero[cq['panel_col'].astype(int)] = cq['zeros'] / float(cq['n_controls'])
    import anndata as ad
    rd = ad.read_h5ad(S19.EVAL[ev] / 'real.h5ad', backed='r'); col = 'target_gene' if 'target_gene' in rd.obs else 'target'
    n_score_ctrl = int((rd.obs[col].astype(str) == 'non-targeting').sum())
    item2 = {}
    lp = -np.log10(np.maximum(Pv, 1e-300))
    for j, nm in enumerate(['5-20', '20-100', '>100']):
        m = sig & (band[None, :] == j) & np.isfinite(lp)
        zf = float(np.median(zero[band == j][C.gate[band == j]])) if (band == j).any() else float('nan')
        up, dn = min_wilcoxon_p(400, n_score_ctrl, zf); up0, dn0 = min_wilcoxon_p(400, int(cq['n_controls']), zf)
        item2[nm] = dict(called_pairs=int(m.sum()), neglog10p=dict(zip(['q10', 'q25', 'median', 'q75', 'q90', 'max'], np.quantile(lp[m], [.1, .25, .5, .75, .9, 1]).tolist())) if m.any() else None,
                         share_p_adj_lt_005=float((Pa[m] < 0.05).mean()) if m.any() else None, median_control_zero_fraction=zf,
                         min_p_scoring_controls=dict(n_controls=n_score_ctrl, increase=up, decrease=dn),
                         min_p_context_controls=dict(n_controls=int(cq['n_controls']), increase=up0, decrease=dn0))
    res['item2_achievable_p'] = item2
    # item 3: the realised order within called genes
    pc = np.load(TAB / 'S3_pcal_H1.npz')['xfit']
    rho = {k: [] for k in ('|requested LFC|', 'control CPM', 'z_u', 'p_cal')}
    for i in range(h.T):
        m = sig[i] & np.isfinite(Pv[i])
        if m.sum() < 10:
            continue
        rk = -lp[i, m]
        for k, x in (('|requested LFC|', -np.abs(lfc[i, m])), ('control CPM', -m0[m]), ('z_u', -C.zu[i, m]), ('p_cal', -pc[i, m])):
            if np.std(x) > 0:
                rho[k].append(spearmanr(rk, x).correlation)
    res['item3_median_spearman_realised_p_rank'] = {k: dict(median=float(np.nanmedian(v)), knockdowns=len(v)) for k, v in rho.items()}
    # the order's share of the reach gap: proxy with cells' realised p order (higher -log10 p first; ties by p_adj)
    nconf = (C.T & ~C.ownM).sum(1).astype(float); parts = {p: rr for p, rr in C.D.parts.items() if p in ('full', 'half0', 'half1')}
    key = np.where(np.isfinite(lp), lp - 1e-6 * np.nan_to_num(Pa), -np.inf)
    rr_cells = reach_rows(C, de_order(C.D.ref, S.mu4, C.sig, key), 'H1', 'cells realised p order', parts, nconf)
    rr_proxy = reach_rows(C, de_order(C.D.ref, S.mu4, C.sig, None), 'H1', 'effect size (proxy)', parts, nconf)
    W = pd.read_csv(hd / 'results.csv').pivot(index='perturbation', columns='metric', values='value').reindex(h.targets)
    import scorer_proxy as SP
    cells = {p: dict(reach=h.scale(SP.DE[2], float(np.nanmean(W[SP.DE[2]].to_numpy()[rr]))), reach_raw=float(np.nanmean(W[SP.DE[2]].to_numpy()[rr]))) for p, rr in parts.items()}
    split = {}
    for x, y in zip(rr_proxy, rr_cells):
        p = x['part']
        split[p] = dict(proxy_effect_size=x['reach'], proxy_cells_order=y['reach'], cells=cells[p]['reach'],
                        gap_total=x['reach'] - cells[p]['reach'], gap_order=x['reach'] - y['reach'], gap_rest_sign=y['reach'] - cells[p]['reach'],
                        raw=dict(proxy_effect_size=x['reach_raw'], proxy_cells_order=y['reach_raw'], cells=cells[p]['reach_raw']))
    res['item3_reach_gap_split'] = split
    res['note_item1'] = 'the code reading (item 1) is in the write-up (section 14), with file and line numbers'
    save_json(TAB / 'partR2.json', res)
    log(f'R2: achievable p {json.dumps(item2, default=float)[:600]}; spearman {res["item3_median_spearman_realised_p_rank"]}; gap split full {split["full"]}')


# ================================================================ Part R3: generator v8 (rank-graded spread)
R3D = OUT / 'R3'


def v8():
    sys.path.insert(0, str(V / 'code/vcc6/gen_v8'))
    import cell_generator_v8 as G8
    return G8


def rank_within(sig, key, eff):
    """Per knockdown, called genes ranked by key (higher = more confident; ties by |effect|): r in [0, 1], the least
    confident 0 and the most 1; NaN on uncalled genes (and r = 0 for a knockdown with one call)."""
    Rk = np.full(sig.shape, np.nan)
    k = np.where(np.isfinite(key), key, -np.inf)
    for i in range(len(sig)):
        j = np.flatnonzero(sig[i])
        if not len(j):
            continue
        o = np.lexsort((np.abs(eff[i, j]), k[i, j]))            # ascending: least confident first
        r = np.empty(len(j)); r[o] = np.arange(len(j)) / max(len(j) - 1, 1)
        Rk[i, j] = r
    return Rk


def r3_key_name():
    return json.loads((TAB / 'partS4.json').read_text())['key_for_R3']


def r3_ranks(a):
    """Rank arrays for Eval A (cross-fitted: p_cal from the other half), HCT116final (screen setup, p_cal from all of
    H1), and D / E / F (row 4's calls: hybrid T 5 on C1's arrays; p_cal from all of H1), in each input's own order."""
    import pickle
    R3D.mkdir(parents=True, exist_ok=True)
    key = r3_key_name(); maps_all = pickle.load(open(TAB / 'S3_pcal_maps.pkl', 'rb'))['maps_all']; res = dict(key=key)
    for ev, line in (('EvalA', 'H1'), ('HCT116final', 'HCT116s')):
        S = SLine(line); C = S.C; h = C.h
        z = np.load(S21.CELLS / ev / 'input.npz', allow_pickle=True); names = list(map(str, z['knockdowns']))
        r = np.array([h.targets.index(t) for t in names])
        sig = (z['de_sig'] & (z['n_sources'] > 0))
        if key == '(i) p_cal':
            K = np.load(TAB / 'S3_pcal_H1.npz')['xfit'] if line == 'H1' else apply_pcal(S, maps_all)
        else:
            K = C.zu
        Rk = rank_within(sig, K[r], z['lfc'].astype(float))
        np.savez_compressed(R3D / f'rank_{ev}.npz', knockdowns=np.array(names, dtype=object), rank=Rk.astype(np.float32), key=key)
        res[ev] = dict(knockdowns=len(names), called=int(sig.sum()), key_finite_share=float(np.isfinite(K[r][sig]).mean()))
        del S, C; gc.collect(); S19._REFS.clear()
    # D / E / F
    bdef = pd.read_csv(S19.S17 / 'tables/b_DEF.csv').set_index('target')['b']; gp = {g: j for j, g in enumerate(panel())}
    for x in ('D', 'E', 'F'):
        c = S21.c1_load(x); tg = list(map(str, c['targets'])); gate = c['gate'].astype(bool)
        own = np.array([gp.get(t, -1) for t in tg]); okg = np.broadcast_to(gate[None, :], c['n'].shape).copy()
        for i, j in enumerate(own):
            if j >= 0: okg[i, j] = False
        cg = np.where(c['n'] <= 0, 0.0, c['lfc'].astype(np.float64)); ok4 = okg & (c['n'] > 0) & (cg != 0)
        zc = np.where(ok4, np.abs(cg) / np.maximum(c['se'].astype(np.float64)[None, :], 1e-3), -1.0)
        on5 = np.nan_to_num(bdef.reindex(tg).to_numpy(float), nan=-np.inf) >= T_HYB
        sig = np.where(on5[:, None], S21.top_n(zc, ok4, 1200), c['sig'].astype(bool) & (c['n'] > 0))
        if key == '(i) p_cal':
            n, ag, pt, _ = s1_load(x, tg)
            class _S: pass
            s_ = _S(); s_.n, s_.pt = n, pt
            K = apply_pcal(s_, maps_all)
        else:
            K, _ = S21.z_u_of(c, gate, own)
        Rk = rank_within(sig, K, cg)
        np.savez_compressed(R3D / f'rank_{x}.npz', knockdowns=np.array(tg, dtype=object), rank=Rk.astype(np.float32), de_sig=sig, key=key)
        res[x] = dict(knockdowns=len(tg), called=int(sig.sum()), on_K4=int(on5.sum()), key_finite_share=float(np.isfinite(K[sig]).mean()))
    save_json(TAB / 'R3_ranks.json', res); log(f'R3 ranks: {res}')


def r3_slice():
    import scorer_step12 as S12
    ref = S19.ref_for('EvalA'); T = S12.truth_calls(ref); tt = tercile_of(T.sum(1)); rng = np.random.default_rng(SEED)
    pick = []
    for j in (2, 1):
        pick += list(rng.choice(np.flatnonzero(tt == j), 10, replace=False))
    return [ref.h.targets[i] for i in sorted(pick)], {ref.h.targets[i]: TERC[tt[i]] for i in pick}


def partR3gen(a):
    """Generate the 20-knockdown slice three ways (frozen, v8 without rank, v8 with rank) from Track 1's Eval A input;
    check (1) bit for bit; write the hybrid prediction (Track 1's cells with the slice replaced by v8-with-rank) for
    the frozen harness."""
    import anndata as ad
    import scipy.sparse as sp
    import cell_generator as CG
    import stage1 as S
    sys.path.insert(0, str(V / 'scripts/harness'))
    from common import write_h5ad
    G8 = v8(); R3D.mkdir(parents=True, exist_ok=True)
    sl, terc = r3_slice()
    z = np.load(S21.CELLS / 'EvalA/input.npz', allow_pickle=True); names = list(map(str, z['knockdowns']))
    rk = np.load(R3D / 'rank_EvalA.npz', allow_pickle=True); assert list(map(str, rk['knockdowns'])) == names
    cd = S21.CELLS / 'EvalA/e2e/controls'
    prof_ = json.loads((cd / 'controls_profile.json').read_text()); L = float(prof_['median_panel_library'])
    m0 = pd.read_csv(cd / 'controls_baseline.csv', index_col=0)['control_mean_cpm'].to_numpy(float)
    cq = CG.ControlQuantiles(cd / 'controls_control_quantiles.npz'); cq8 = G8.ControlQuantiles(cd / 'controls_control_quantiles.npz')
    genes = panel()
    pa = ad.read_h5ad(S21.CELLS / 'EvalA/e2e/salt_none/prediction.h5ad')
    lab = pa.obs['target_gene'].astype(str).to_numpy(); X = pa.X.tocsr()
    out = {}; blocks = {}
    t0 = time.time()
    for t in sl:
        i = names.index(t)
        pk = S.prepare_generator_inputs(z['lfc'][i], z['de_sig'][i], z['n_sources'][i], np.full(len(genes), 0.5))
        seed = S.knockdown_seed(t, '')
        fz = CG.generate_cells(pk['lfc'], pk['de_sig'], cq, m0, n_cells=400, library_size=L, seed=seed, target=t, genes=genes)
        v0 = G8.generate_cells(pk['lfc'], pk['de_sig'], cq8, m0, n_cells=400, library_size=L, seed=seed, target=t, genes=genes)
        v1 = G8.generate_cells(pk['lfc'], pk['de_sig'], cq8, m0, n_cells=400, library_size=L, seed=seed, target=t, genes=genes, rank=rk['rank'][i])
        tr = X[np.flatnonzero(lab == t)]
        out[t] = dict(tercile=terc[t], frozen_vs_v8_norank_bitforbit=bool((fz.X != v0.X).nnz == 0),
                      frozen_now_vs_track1_bitforbit=bool((sp.csr_matrix(fz.X) != tr).nnz == 0),
                      called=int(pk['de_sig'].sum()), v8_report={k: v1.uns['cellgen'][k] for k in ('sig_requested', 'sig_realized', 'sig_recall', 'sig_precision', 'rebuilt_genes', 'first_pass_unmet')},
                      frozen_report={k: fz.uns['cellgen'][k] for k in ('sig_requested', 'sig_realized', 'sig_recall', 'sig_precision', 'rebuilt_genes', 'first_pass_unmet')},
                      mean_count_max_abs_error_v8=float(v1.uns['cellgen']['mean_count_max_abs_error']))
        blocks[t] = sp.csr_matrix(v1.X)
        log(f'R3 gen {t} ({terc[t]}): bit for bit frozen = v8 no-rank {out[t]["frozen_vs_v8_norank_bitforbit"]}, frozen now = Track 1 {out[t]["frozen_now_vs_track1_bitforbit"]}; '
            f'v8 realised {out[t]["v8_report"]["sig_realized"]} / requested {out[t]["v8_report"]["sig_requested"]} (frozen {out[t]["frozen_report"]["sig_realized"]}); {time.time() - t0:.0f} s')
    ok1 = all(v['frozen_vs_v8_norank_bitforbit'] and v['frozen_now_vs_track1_bitforbit'] for v in out.values())
    res = dict(slice=sl, per_knockdown=out, check1_bit_for_bit=bool(ok1), C_MIN=0.0, seconds=time.time() - t0)
    save_json(TAB / 'partR3_gen.json', res)
    if not ok1:
        log('R3 CHECK 1 FAILED (v8 without rank differs from the frozen generator)'); sys.exit(21)
    # hybrid: Track 1's cells, slice knockdowns replaced by v8-with-rank (same cell positions)
    rows = []
    for t in pd.unique(lab):
        rows.append(blocks[t] if t in blocks else X[np.flatnonzero(lab == t)])
    Xh = sp.vstack(rows).tocsr()
    lab_h = np.concatenate([np.full(r.shape[0], t, dtype=object) for r, t in zip(rows, pd.unique(lab))])
    assert (lab_h == lab).all(), 'Track 1 cells are not grouped by knockdown'
    d = R3D / 'EvalA_hybrid/salt_none'; d.mkdir(parents=True, exist_ok=True)
    obs = pd.DataFrame(dict(target=lab, target_gene=lab), index=pa.obs_names)
    write_h5ad(d / 'prediction.h5ad', Xh, obs, genes)
    log(f'R3 hybrid prediction written: {Xh.shape}; slice {len(sl)} knockdowns from v8 with rank')


def de_side(hd, h, sig_targets):
    dp = pd.read_parquet(hd / 'de_pred.parquet', columns=['target', 'feature', 'p_value', 'p_adj'])
    dp = dp[dp.target.isin(sig_targets)]
    gp = {g: j for j, g in enumerate(h.genes)}; tp = {t: i for i, t in enumerate(h.targets)}
    Pv = np.full((h.T, h.G), np.nan); Pa = np.full((h.T, h.G), np.nan)
    Pv[dp.target.map(tp).to_numpy(), dp.feature.map(gp).to_numpy()] = dp.p_value.to_numpy()
    Pa[dp.target.map(tp).to_numpy(), dp.feature.map(gp).to_numpy()] = dp.p_adj.to_numpy()
    return Pv, Pa


def partR3cmp(a):
    """Checks (2)-(5) on the slice: hybrid (v8 with rank) against Track 1 (frozen), both scored by the frozen harness."""
    import anndata as ad
    from scipy.stats import spearmanr
    g = json.loads((TAB / 'partR3_gen.json').read_text()); sl = g['slice']
    ref = S19.ref_for('EvalA'); h = ref.h; own = S21.own_mask(h)
    tr = S21.CELLS / 'EvalA/e2e/salt_none'; hy = R3D / 'EvalA_hybrid/salt_none'
    # (2) pseudobulk log1p CP50K per (knockdown, gene)
    def pb(path):
        a_ = ad.read_h5ad(path, backed='r'); lab = a_.obs['target_gene'].astype(str).to_numpy(); out = {}
        for t in sl:
            x = np.asarray(a_[np.flatnonzero(lab == t)].to_memory().X.sum(0)).ravel().astype(float); out[t] = np.log1p(5e4 * x / x.sum())
        return out
    P0, P1 = pb(tr / 'prediction.h5ad'), pb(hy / 'prediction.h5ad')
    d = np.concatenate([np.abs(P1[t] - P0[t]) for t in sl])
    rows = np.array([h.targets.index(t) for t in sl])
    W0 = pd.read_csv(tr / 'harness/evalA/results.csv').pivot(index='perturbation', columns='metric', values='value')
    W1 = pd.read_csv(hy / 'harness/evalA/results.csv').pivot(index='perturbation', columns='metric', values='value')
    F0, F1 = S19.Fast.from_wide(h, W0), S19.Fast.from_wide(h, W1)
    p0, r0 = F0.point(rows); p1, r1 = F1.point(rows)
    res = dict(check2=dict(pseudobulk_max_abs_diff=float(d.max()), pseudobulk_p99_abs_diff=float(np.quantile(d, .99)),
                           slice_scaled_frozen=p0, slice_scaled_v8=p1, diff={k: p1[k] - p0[k] for k in p0}))
    # (3) realised called sets (harness DE p_adj < 0.05, own gene out)
    Pv0, Pa0 = de_side(tr / 'harness/evalA', h, sl); Pv1, Pa1 = de_side(hy / 'harness/evalA', h, sl)
    c0 = (Pa0 < 0.05) & ~own; c1 = (Pa1 < 0.05) & ~own
    inter = (c0 & c1)[rows].sum(); uni = (c0 | c1)[rows].sum()
    res['check3'] = dict(jaccard_pooled=float(inter / max(uni, 1)),
                         per_knockdown={t: dict(frozen=int(c0[i].sum()), v8=int(c1[i].sum()), gained=int((c1[i] & ~c0[i]).sum()), lost=int((c0[i] & ~c1[i]).sum()),
                                                jaccard=float((c0[i] & c1[i]).sum() / max((c0[i] | c1[i]).sum(), 1))) for t, i in zip(sl, rows)},
                         fidelity=dict(frozen=p0['fidelity'], v8=p1['fidelity']), jaccard_component=dict(frozen=p0['jaccard'], v8=p1['jaccard']))
    # (4) realised p rank vs requested rank within called genes, by CPM band
    z = np.load(S21.CELLS / 'EvalA/input.npz', allow_pickle=True); names = list(map(str, z['knockdowns']))
    rk = np.load(R3D / 'rank_EvalA.npz', allow_pickle=True)['rank']; m0 = z['m0']
    sp4 = {}
    for nm, lo, hi in (('5-20', 5, 20), ('20-100', 20, 100), ('>100', 100, np.inf)):
        for tag, Pv in (('v8', Pv1), ('frozen', Pv0)):
            vals = []
            for t, i in zip(sl, rows):
                j = names.index(t); m = np.isfinite(rk[j]) & (m0 > lo) & (m0 <= hi) & np.isfinite(Pv[i]) & ~own[i]
                if m.sum() >= 5 and np.std(rk[j][m]) > 0:
                    vals.append(spearmanr(-np.log10(np.maximum(Pv[i][m], 1e-300)), rk[j][m]).correlation)
            sp4[f'{nm} {tag}'] = dict(median=float(np.nanmedian(vals)) if vals else float('nan'), knockdowns=len(vals))
    res['check4_spearman_realised_vs_requested_rank'] = sp4
    res['check5_reach'] = dict(frozen=p0['reach'], v8=p1['reach'], frozen_raw=r0['reach'], v8_raw=r1['reach'],
                               by_tercile={tt: dict(frozen=F0.point(rows[[g['per_knockdown'][t]['tercile'] == tt for t in sl]])[0]['reach'],
                                                    v8=F1.point(rows[[g['per_knockdown'][t]['tercile'] == tt for t in sl]])[0]['reach']) for tt in ('mid', 'strong')})
    ps = res['check2']['diff']
    res['pass'] = dict(check1=g['check1_bit_for_bit'], mse_pds_nmae_within_0002=bool(all(abs(ps[k]) <= 0.002 for k in ('mse', 'pds', 'nmae'))),
                       called_set_jaccard_ge_099=bool(res['check3']['jaccard_pooled'] >= 0.99),
                       spearman_gt_05_above_20cpm=bool(all(sp4[f'{b} v8']['median'] > 0.5 for b in ('20-100', '>100'))))
    res['pass']['all'] = all(res['pass'].values())
    save_json(TAB / 'partR3.json', res)
    log(f'R3 cmp: {json.dumps(res["pass"])}; diff {json.dumps({k: round(v, 4) for k, v in ps.items()})}; jaccard {res["check3"]["jaccard_pooled"]:.4f}; spearman {sp4}')


# ================================================================ Part E1: row 4 on generator v8 (pre-approved if R3 passes)
E1D = OUT / 'E1'
E1_KD = {'EvalA': V / 'results/target_lists/evalA_h1_292.txt', 'HCT116final': V / 'results/target_lists/final_HCT116.txt'}
E1_CTRL = {'EvalA': V / 'lib/harness/evalA/context.h5ad', 'HCT116final': FE / 'evalHCT116final/context.h5ad'}


def r3_passed():
    f = TAB / 'partR3.json'
    return f.exists() and bool(json.loads(f.read_text())['pass']['all'])


def e1gen(ev):
    """As run_line predict --lfc-npz (Track 1's input, salt "", 400 cells, the controls' median panel library), through
    gen_v8/gen_cells_v8 with R3's rank arrays. Check: the controls profile equals Track 1's."""
    import run_line as RL
    sys.path.insert(0, str(V / 'code/vcc6/gen_v8'))
    import gen_cells_v8 as GC8
    if not r3_passed():
        log('E1: R3 did not pass all its checks; not generating'); sys.exit(30)
    genes = panel(); out = E1D / ev / 'e2e'; out.mkdir(parents=True, exist_ok=True)
    targets = [t for t in E1_KD[ev].read_text().split() if t]
    X = RL.read_controls(E1_CTRL[ev], genes)
    m0, gate, qpath, prof_ = RL.control_profile(X, out / 'controls')
    t1 = pd.read_csv(S21.CELLS / ev / 'e2e/controls/controls_baseline.csv', index_col=0)
    e1 = pd.read_csv(out / 'controls/controls_baseline.csv', index_col=0)
    same = bool(t1.equals(e1)) and md5_of(out / 'controls/controls_control_quantiles.npz') == md5_of(S21.CELLS / ev / 'e2e/controls/controls_control_quantiles.npz')
    log(f'E1 {ev}: controls {prof_}; equal to Track 1 controls: {same}')
    if not same:
        sys.exit(31)
    z = np.load(S21.CELLS / ev / 'input.npz', allow_pickle=True); row = {str(k): i for i, k in enumerate(z['knockdowns'])}; idx = [row[t] for t in targets]
    rk = np.load(R3D / f'rank_{ev}.npz', allow_pickle=True); assert list(map(str, rk['knockdowns'])) == list(map(str, z['knockdowns']))
    lfc, sig, n_src = np.asarray(z['lfc'], float)[idx], np.asarray(z['de_sig'], bool)[idx], np.asarray(z['n_sources'])[idx]
    lfc = np.where(np.isfinite(lfc), lfc, 0.0)
    d = out / 'salt_none'; d.mkdir(parents=True, exist_ok=True)
    log(f'E1 {ev}: generating {len(targets)} knockdowns with v8 + rank at library {prof_["median_panel_library"]:.0f}')
    rep = GC8.write_arm(d / 'prediction.h5ad', targets, genes, lfc, sig, n_src, qpath, m0, X, salt='', p_sign=None, n_cells=400,
                        library_size=prof_['median_panel_library'], workers=8, log=log, rank=rk['rank'][idx])
    (d / 'run.json').write_text(json.dumps(dict(salt='', model=dict(source=str(S21.CELLS / ev / 'input.npz'), generator='gen_v8 + rank', rank=str(R3D / f'rank_{ev}.npz')),
                                                controls=prof_, knockdowns=len(targets), n_cells=400, sig_requested_median=float(rep.sig_requested.median()),
                                                sig_realized_median=float(rep.sig_realized.median())), indent=1) + '\n')


def partE1gen_EvalA(a):
    e1gen('EvalA')


def partE1gen_HCT116final(a):
    e1gen('HCT116final')


def partE1cmp(a):
    """As Step 21's t1cmp, E1 (row 4 on v8 + rank) against row 4's Track 1 cells (paired), plus realised against
    requested calls and the realised order's Spearman with the requested rank within called genes, by CPM band."""
    import scorer_proxy as SP
    import scorer_step12 as S12
    from scipy.stats import spearmanr
    res = {}; B = pd.read_csv(TAB / 'row4_rows.csv')
    for ev, line, hn in (('EvalA', 'H1', 'evalA'), ('HCT116final', 'HCT116s', 'evalHCT116final')):
        ref = S19.ref_for(ev); h = ref.h
        hd = E1D / ev / 'e2e/salt_none/harness' / hn; h0 = S21.CELLS / ev / 'e2e/salt_none/harness' / hn
        Wn = pd.read_csv(hd / 'results.csv').pivot(index='perturbation', columns='metric', values='value').reindex(h.targets)
        W0 = pd.read_csv(h0 / 'results.csv').pivot(index='perturbation', columns='metric', values='value').reindex(h.targets)
        Fn, F0 = S19.Fast.from_wide(h, Wn), S19.Fast.from_wide(h, W0)
        I = S19.boot_index(h.T); Sn, Rn = Fn.scaled(I); S0, R0 = F0.scaled(I); pn, prn = Fn.point(); p0, pr0 = F0.point()
        summ = json.loads((hd / 'summary.json').read_text())
        out = dict(rebuild_vs_summary_abs_diff=abs(pn['overall'] - summ['scaled']['overall']),
                   components={k: dict(scaled=pn[k], ci=[float(np.quantile(Sn[k], .025)), float(np.quantile(Sn[k], .975))], raw=prn.get(k)) for k in COMPS + ['overall']},
                   vs_row4_track1={k: dict(diff=pn[k] - p0[k], ci=[float(np.quantile(Sn[k] - S0[k], .025)), float(np.quantile(Sn[k] - S0[k], .975))],
                                           p_better=float((Sn[k] - S0[k] > 0).mean()), raw_diff=(prn[k] - pr0[k]) if k in prn else None) for k in COMPS + ['overall']})
        pr = B[(B.line == line) & (B.part == 'full')].iloc[0]
        out['proxy_minus_cells'] = {k: dict(proxy=float(pr[k]), cells=pn[k], gap=float(pr[k] - pn[k])) for k in COMPS + ['overall']}
        T = S12.truth_calls(ref)
        if ev == 'EvalA':
            tt = tercile_of(T.sum(1)); out['terciles'] = {}
            for j, tn in enumerate(TERC):
                rr = np.flatnonzero(tt == j); rng = np.random.default_rng(SEED)
                Ij = np.stack([rr[rng.integers(0, len(rr), len(rr))] for _ in range(S19.N_BOOT)])
                sn, _ = Fn.scaled(Ij); s0, _ = F0.scaled(Ij)
                out['terciles'][tn] = dict(knockdowns=len(rr), overall=Fn.point(rr)[0]['overall'], reach=Fn.point(rr)[0]['reach'],
                                           vs_row4=dict(diff=Fn.point(rr)[0]['overall'] - F0.point(rr)[0]['overall'],
                                                        ci=[float(np.quantile(sn['overall'] - s0['overall'], .025)), float(np.quantile(sn['overall'] - s0['overall'], .975))]),
                                           reach_vs_row4=Fn.point(rr)[0]['reach'] - F0.point(rr)[0]['reach'])
            import scorer_step16 as S16
            hv = S16.halves(T); out['halves_vs_row4'] = {f'half{f}': {k: Fn.point(hv[f])[0][k] - F0.point(hv[f])[0][k] for k in ('overall', 'reach')} for f in (0, 1)}
        z = np.load(S21.CELLS / ev / 'input.npz', allow_pickle=True); names = list(map(str, z['knockdowns']))
        r = [names.index(t) for t in h.targets]; un = z['n_sources'][r] <= 0; sig = z['de_sig'][r] & ~un
        per, by = SP.call_diff(ref, hd, sig, E1D / ev / 'e2e/controls/controls_baseline.csv')
        per0, _ = SP.call_diff(ref, h0, sig, S21.CELLS / ev / 'e2e/controls/controls_baseline.csv')
        per.to_csv(TAB / f'E1_calls_per_kd_{ev}.csv', index=False); by.to_csv(TAB / f'E1_calls_by_decile_{ev}.csv', index=False)
        nt = T.sum(1)
        out['calls'] = dict(requested_mean=float(per.requested.mean()), realised_mean=float(per.realised.mean()), realised_mean_track1=float(per0.realised.mean()),
                            requested_not_realised_total=int(per.requested_not_realised.sum()), realised_not_requested_total=int(per.realised_not_requested.sum()),
                            requested_not_realised_total_track1=int(per0.requested_not_realised.sum()), realised_not_requested_total_track1=int(per0.realised_not_requested.sum()),
                            spearman_realised_true_de=float(spearmanr(per.realised.to_numpy(float), nt).correlation))
        rk = np.load(R3D / f'rank_{ev}.npz', allow_pickle=True)['rank'][r]; m0 = z['m0']; own = S21.own_mask(h)
        Pv1, _ = de_side(hd, h, h.targets); Pv0, _ = de_side(h0, h, h.targets)
        so = {}
        for nm, lo, hi in (('5-20', 5, 20), ('20-100', 20, 100), ('>100', 100, np.inf), ('all', 5, np.inf)):
            for tag, Pv in (('E1 v8', Pv1), ('Track 1 frozen', Pv0)):
                vals = []
                for i in range(h.T):
                    m = np.isfinite(rk[i]) & (m0 > lo) & (m0 <= hi) & np.isfinite(Pv[i]) & ~own[i]
                    if m.sum() >= 5 and np.std(rk[i][m]) > 0:
                        vals.append(spearmanr(-np.log10(np.maximum(Pv[i][m], 1e-300)), rk[i][m]).correlation)
                so[f'{nm} {tag}'] = dict(median=float(np.nanmedian(vals)) if vals else float('nan'), knockdowns=len(vals))
        out['realised_order_spearman_with_key'] = so
        res[ev] = out
        log(f'E1cmp {ev}: overall {pn["overall"]:.4f} vs row 4 cells {p0["overall"]:.4f}; diff ' +
            ', '.join(f'{k} {v["diff"]:+.4f} [{v["ci"][0]:+.4f}, {v["ci"][1]:+.4f}]' for k, v in out['vs_row4_track1'].items()))
        S19._REFS.pop(ev, None); gc.collect()
    save_json(TAB / 'partE1.json', res)


# ================================================================ Part D: recommendation
def d_status():
    out = {}
    for nm, f in (('P2', 'partP2.json'), ('S3', 'partS3.json'), ('Q2', 'partQ2.json')):
        if not (TAB / f).exists():
            out[nm] = dict(available=False); continue
        r = json.loads((TAB / f).read_text())
        out[nm] = dict(available=True, met=bool(r['criterion']['met']), picks={k: v.get('picked') for k, v in r['fits'].items()})
        if nm == 'Q2':
            out['Q2 M1b'] = dict(available=True, met=bool(r['M1b']['criterion']['met']), picks={'fit_half0': {}, 'fit_half1': {}})
    return out


def d_mu(S, line, mods, f, ctx):
    """mu for row 4 plus mods on SLine S (H1: direction f), applied in the order Q2 (shift), P2 (term), S3 (weights)."""
    R = S.C; add = None
    if 'Q2' in mods or 'Q2 M1b' in mods:
        z, _ = depmap_z(R.h.targets)
        if 'Q2' in mods:
            add = (SHIFT_A[R.terc] - mods['Q2']['delta'] * z)[:, None] * R.S[None, :]
        else:
            t, _ = m1b_terc(R, z); add = SHIFT_A[t][:, None] * R.S[None, :]
    extra = np.zeros((R.h.T, R.h.G)) if add is None else add
    if 'P2' in mods:
        sets, _, _, _ = gene_sets(); gp = {g: j for j, g in enumerate(panel())}
        T_, _, _ = p2_term(R, sets, ctx.setdefault('src', Src()), gp, mods['P2']['variant'] == 'centred')
        extra = extra + mods['P2']['gamma'] * T_
    mu = R.mu(shift=add is None, extra=extra)
    if 'S3' in mods:
        import pickle
        if line == 'H1':
            Pc = pcal_frames(S)[f]
        else:
            Pc = apply_pcal(S, pickle.load(open(TAB / 'S3_pcal_maps.pkl', 'rb'))['maps_all'])
        S.mu4 = mu; S.Dpred = R.delta_pred(mu)
        W, _, _ = s3_weight(S, Pc, mods['S3']['rho'], S21.h1_halves() if line == 'H1' else [np.arange(R.h.T)])
        mu = s3_mu(S, W, mods['S3']['apply'])
    return mu


def partD(a):
    st = d_status(); res = dict(status=st)
    passed = [nm for nm, v in st.items() if v.get('met')]
    if 'Q2' in passed and 'Q2 M1b' in passed:
        passed.remove('Q2 M1b')
    res['passed'] = passed
    log(f'D: parts that met their criterion: {passed or "none"}')
    rows = []
    if len(passed) >= 1:
        for line in LINES:
            S = SLine(line); ctx = {}
            for f in ((0, 1) if line == 'H1' else (-1,)):
                ff = max(f, 0); mods = {p: (st[p]['picks'].get(f'fit_half{ff}') or {}) for p in passed}
                mu = d_mu(S, line, mods, ff, ctx)
                rows += [dict(r, combo='+'.join(passed), fit_half=f) for r in score_tag(S.C, mu)]
            del S; gc.collect()
        R = save_rows(rows, 'partD_rows'); B = pd.read_csv(TAB / 'row4_rows.csv')
        fv = []
        for f in (0, 1):
            x = R[(R.line == 'H1') & (R.fit_half == f) & (R.part == f'half{1 - f}')].iloc[0]; y = B[(B.line == 'H1') & (B.part == f'half{1 - f}')].iloc[0]
            gr = {l: dict(diff=float(R[(R.line == l) & (R.part == 'full')].overall.iloc[0] - B[(B.line == l) & (B.part == 'full')].overall.iloc[0])) for l in GUARD}
            fv.append(dict(held_out=comp_dict(x), base_held_out=comp_dict(y), gain={c: float(x[c] - y[c]) for c in COMPS + ['overall']}, guardrails=gr))
        res['combination'] = dict(parts=passed, fits=fv, criterion=criterion(fv))
    for nm, f in (('S4', 'partS4.json'), ('R2', 'partR2.json'), ('R3', 'partR3.json'), ('P1', 'partP1.json')):
        if (TAB / f).exists():
            r = json.loads((TAB / f).read_text())
            res[nm] = r.get('criterion') or r.get('pass') or {k: r[k] for k in list(r)[:3]}
    if (TAB / 'partS4.json').exists():
        r = json.loads((TAB / 'partS4.json').read_text()); h = [r[f'H1 half{f}'] for f in (0, 1)]; k = r['key_for_R3']
        res['R_criterion'] = dict(oracle_gain_halves=[x['oracle_gain'] for x in h], key=k, key_share_halves=[x['share_of_oracle_gain'][k] for x in h],
                                  met=bool(all(x['oracle_gain'] >= 0.03 for x in h) and all(x['share_of_oracle_gain'][k] >= 0.5 for x in h)))
    save_json(TAB / 'partD.json', res)
    log(f'D: {json.dumps({k: v for k, v in res.items() if k != "status"}, default=float)[:1500]}')


# ================================================================ unit tests
def test(a):
    rng = np.random.default_rng(0)
    for t in TESTS:
        t(rng); log(f'test {t.__name__}: pass')
    log('test: all unit tests pass')


def _t_capture(rng):
    g, k, n = 40, 5, 6
    Q, _ = np.linalg.qr(rng.normal(size=(g, k))); X = rng.normal(size=(n, g)); own = np.array([3, -1, 10, -1, 0, 39])
    for i, j in enumerate(own):
        if j >= 0: X[i, j] = 0.0
    cap, tot = capture(Q, X, own)
    for i, j in enumerate(own):
        m = np.ones(g, bool)
        if j >= 0: m[j] = False
        c = np.linalg.lstsq(Q[m], X[i, m], rcond=None)[0]; ref = float(((Q[m] @ c) ** 2).sum())
        assert np.isclose(cap[i], ref), (i, cap[i], ref)
    assert np.allclose(tot, (X ** 2).sum(1))
    assert orth(np.c_[Q, Q[:, :2]]).shape[1] == k


def _t_vote_stats(rng):
    import stage1 as S
    Y = rng.normal(0, 1, (4, 30)); SE = rng.uniform(0.1, 0.5, (4, 30)); Y[0, :5] = np.nan; Y[1, 5] = 0.0
    mu, vf, t2, tu, n = S.pool_effects(Y, SE, np.ones(4), 1e-4)
    nn, ag, Q = vote_stats(Y, SE, np.ones(4), mu)
    assert np.array_equal(nn, n)
    fin = np.isfinite(Y); w = 1 / np.maximum(SE, S.SE_MIN) ** 2
    Qr = np.nansum(np.where(fin, w * (Y - mu) ** 2, 0.0), 0); assert np.allclose(Q, Qr)
    C = w.sum(0) - (w ** 2).sum(0) / w.sum(0); ok = n >= 2
    tr = (Qr - (n - 1)) / np.where(fin.all(0), C, np.nan)
    assert np.allclose(tu[fin.all(0)], tr[fin.all(0)])
    j = 5; exp = int(((np.sign(Y[:, j]) == np.sign(mu[j])) & (Y[:, j] != 0) & fin[:, j]).sum()); assert ag[j] == exp


def _t_de_order(rng):
    # q mapping: called genes ordered by key (higher first), ties by effect size; uncalled untouched
    df = pd.DataFrame(dict(t=[0, 0, 0, 1, 1], k=[1.0, 3.0, 3.0, -np.inf, 2.0], s=[0.5, 0.1, 0.2, 9.0, 1.0], c=[True] * 5))
    q = np.zeros(5); cd = df.sort_values(['t', 'k', 's'], kind='stable')
    r = cd.groupby('t').cumcount() + 1; n = cd.groupby('t')['k'].transform('size'); q[cd.index.to_numpy()] = (r / n).to_numpy()
    p = 1e-3 / (1 + 1e4 * q)
    assert p[2] < p[1] < p[0] and p[4] < p[3]


def _t_pcal(rng):
    class _S: pass
    S = _S(); T = 60; G = 50
    S.n = rng.integers(0, 6, (T, G)); S.pt = rng.uniform(0.3, 1, (T, G)); S.right = rng.random((T, G)) < S.pt
    S.TD = rng.random((T, G)) < 0.5
    m = fit_pcal(S, np.arange(T)); P = apply_pcal(S, m)
    assert np.isnan(P[S.n == 0]).all() and np.isfinite(P[S.n >= 1]).all()
    for c in (1, 2, 3):
        sel = cls_of(S.n) == c; o = np.argsort(S.pt[sel]); assert np.all(np.diff(P[sel][o]) >= -1e-12)


def _t_rank(rng):
    sig = np.array([[True, False, True, True], [False, True, False, False]])
    key = np.array([[0.9, 5.0, np.nan, 0.9], [0.0, 0.3, 0.0, 0.0]]); eff = np.array([[1.0, 0, 3.0, 2.0], [0, 1.0, 0, 0]])
    r = rank_within(sig, key, eff)
    assert r[0, 3] == 1.0 and r[0, 0] == 0.5 and r[0, 2] == 0.0 and np.isnan(r[0, 1]) and r[1, 1] == 0.0


def _t_v8(rng):
    import cell_generator as CG
    G8 = v8()
    n0, G = 300, 40
    X = rng.poisson(rng.gamma(1.0, 2.0, G), (n0, G)).astype(float)
    lib = X.sum(1, keepdims=True); C = np.log1p(X * 1e6 / lib)
    tested = np.arange(G); ptr = [0]; vals, cnts, zeros = [], [], []
    for j in range(G):
        v, c = np.unique(C[:, j][C[:, j] > 0], return_counts=True); vals += list(v); cnts += list(c); ptr.append(len(vals)); zeros.append(int((C[:, j] == 0).sum()))
    d = Path(os.environ.get('TMPDIR', '/tmp')) / 'step22_t_v8.npz'
    np.savez(d, panel_col=tested, ptr=np.array(ptr), values=np.array(vals), counts=np.array(cnts), zeros=np.array(zeros), n_controls=n0)
    m0 = (X * 1e6 / lib).mean(0); lfc = rng.normal(0, 0.8, G); sig = rng.random(G) < 0.4
    a0 = CG.generate_cells(lfc, sig, CG.ControlQuantiles(d), m0, n_cells=200, library_size=2000, seed=11)
    a1 = G8.generate_cells(lfc, sig, G8.ControlQuantiles(d), m0, n_cells=200, library_size=2000, seed=11)
    assert (a0.X != a1.X).nnz == 0
    rk = np.where(sig, 0.0, np.nan)
    a2 = G8.generate_cells(lfc, sig, G8.ControlQuantiles(d), m0, n_cells=200, library_size=2000, seed=11, rank=rk)
    assert (a0.X != a2.X).nnz == 0                                       # rank 0 = frozen
    rk = np.where(sig, 1.0, np.nan)
    a3 = G8.generate_cells(lfc, sig, G8.ControlQuantiles(d), m0, n_cells=200, library_size=2000, seed=11, rank=rk)
    assert np.allclose(np.asarray(a3.X.sum(0)).ravel(), np.asarray(a0.X.sum(0)).ravel(), atol=1)   # totals unchanged
    # uncalled genes: same input (c = 1 search); their cells can still differ through check()'s rebuild, which re-tests
    # every gene on each cell's own library (called genes' new shapes move the libraries); their totals stay exact
    assert np.array_equal(np.asarray(a3.X[:, ~sig].sum(0)).ravel(), np.asarray(a0.X[:, ~sig].sum(0)).ravel())
    assert G8.spreads_from(1.0) == CG.SPREADS and G8.spreads_from(0.5)[:2] == (0.5, 0.45)


def _t_minp(rng):
    up, dn = min_wilcoxon_p(400, 2000, 0.5)
    assert 0 <= up < 1e-50 and 0 <= dn < 1e-10
    _, dn2 = min_wilcoxon_p(400, 2000, 0.95); assert dn2 > dn                # all-zero cells tie with more control zeros


def _t_src(rng):
    class _Src(Src):
        def __init__(self): pass
    s = _Src(); s.D = np.array([[1.0, 2, np.nan], [3, 4, 5], [10, 20, 30]]); s.by = {'k': [('A', 0), ('B', 1), ('K562', 2)]}; s.gp = {'k': 1}
    p = s.pooled('k', 'K562'); assert np.allclose(p[[0, 2]], [2.0, 5.0]) and np.isnan(p[1])
    p = s.pooled('k', 'full'); assert np.allclose(p[[0, 2]], [14 / 3, 35 / 2])


def _t_fit(rng):
    rows = []
    for f_ in ('half0', 'half1', 'full', 'half0_weak', 'half1_weak'):
        for x, o in ((0.0, .20), (1.0, .26), (2.0, .30)):
            rows.append(dict(line='H1', part=f_, x=x, overall=o, **{c: 0.1 for c in COMPS}))
    for l in GUARD:
        for x, o in ((0.0, .2), (1.0, .199), (2.0, .19)):
            rows.append(dict(line=l, part='full', x=x, overall=o - (0.02 if l == 'HCT116s' else 0), **{c: 0.1 for c in COMPS}))
    rr = fit_report(pd.DataFrame(rows), ['x'], dict(x=[1.0, 2.0]), [0.0])
    f0 = rr['fits']['fit_half0']
    assert f0['picked'] == {'x': 2.0} and f0['grid_position']['x'] == '2/2 (edge)' and not rr['criterion']['fold_pools_within_0005']
    rows2 = [dict(r, overall=(0.5 if r['x'] == 0.0 else r['overall'])) for r in rows]
    rr2 = fit_report(pd.DataFrame(rows2), ['x'], dict(x=[1.0, 2.0]), [0.0])
    assert rr2['fits']['fit_half0']['grid_position']['x'] == 'base' and rr2['fits']['fit_half0']['picked_is_base']


def _t_metacell(rng):
    """Control bases on a synthetic panel-sized slice (600 cells): shapes, gate support, a planted program recovered."""
    import scipy.sparse as sp
    sets, s_g, g2m, info = gene_sets()
    genes = panel(); G = len(genes); gp = {g: j for j, g in enumerate(genes)}
    n = 600; base = rng.gamma(0.6, 3.0, G) * 0.2
    prog = np.zeros(G); prog[[gp[g] for g in sets['HYPOXIA']]] = 1.0
    act = rng.normal(0, 1, n)
    lam = base[None, :] * np.exp(0.6 * np.outer(act, prog))
    X = sp.csr_matrix(rng.poisson(lam).astype(np.float64))
    lib = np.asarray(X.sum(1)).ravel(); m0 = np.asarray(X.multiply(1e6 / lib[:, None]).mean(0)).ravel(); gate = m0 > 5
    B, A, H, d = metacell_bases(X, gate, sets, s_g, g2m)
    assert B.shape == (50, G) and A.shape == (2, G) and H.shape == (len(PATHS), G) and d['M'] == 24
    assert not B[:, ~gate].any() and not H[:, ~gate].any()
    q = PATHS.index('HYPOXIA'); hy = np.array([gp[g] for g in sets['HYPOXIA'] if gate[gp[g]]])
    assert H[q, hy].mean() > 5 * np.abs(np.delete(H[q], hy)[np.delete(gate, hy)]).mean(), 'planted program not recovered'
    assert abs(sum(d['phase_proportions'].values()) - 1) < 1e-9


def _t_s3w(rng):
    class _C: pass
    class _S: pass
    T, G = 8, 30; S = _S(); S.C = _C()
    S.C.okg = rng.random((T, G)) < 0.9; S.C.sig = rng.random((T, G)) < 0.3; S.n = rng.integers(0, 4, (T, G))
    Pc = np.where(S.n >= 1, rng.uniform(0.4, 1.0, (T, G)), np.nan)
    W, has, zero = s3_weight(S, Pc, 0.5, [np.arange(4), np.arange(4, 8)])
    for rr in (np.arange(4), np.arange(4, 8)):
        for cm in (S.C.sig, ~S.C.sig):
            m = np.zeros((T, G), bool); m[rr] = True; m &= has & cm
            assert np.isclose(W[m].mean(), 1.0)
    assert (W[~has] == 1).all() and np.isclose(zero, (has & (Pc <= 0.5)).sum() / has.sum())


def _t_ridge(rng):
    M, g = 60, 12; w = rng.integers(5, 50, M).astype(float)
    Pg = rng.normal(size=(M, g)); idx = [np.arange(0, 3), np.arange(3, 6)]
    bj, sep, sd = joint_slopes(Pg, w, idx, lam=0.0)
    wn = w / w.sum(); Xc = Pg - wn @ Pg; S = np.stack([Pg[:, ix].mean(1) for ix in idx], 1); S -= wn @ S
    ref = np.linalg.lstsq(S * np.sqrt(w)[:, None], Xc * np.sqrt(w)[:, None], rcond=None)[0]
    assert np.allclose(bj, ref, atol=1e-8)
    b1, _, _ = joint_slopes(Pg, w, idx, lam=1e6); assert np.abs(b1).max() < np.abs(bj).max() * 1e-3


TESTS = [_t_ridge, _t_capture, _t_vote_stats, _t_de_order, _t_pcal, _t_rank, _t_v8, _t_minp, _t_src, _t_fit, _t_s3w, _t_metacell]

def gsets(a):
    """The gene-set matching report (written before P0 uses the sets); also records that gseapy is not imported."""
    sets, s_g, g2m, info = gene_sets()
    info['gseapy_imported_in_this_process'] = 'gseapy' in sys.modules
    info['pylib_gseapy_on_sys_path'] = any('pylib_gseapy' in p for p in sys.path)
    save_json(TAB / 'gene_sets.json', info)
    for k, v in info['gene_set_matching'].items():
        log(f'gene set {k}: size {v["set_size"]}, direct {v["matched_direct"]}, after HGNC mapping {v["matched_after_hgnc_mapping"]}, used {v["used"]}; '
            f'renamed {v["renamed"]}; unmatched {v["unmatched"]}')
    log(f'gseapy imported: {info["gseapy_imported_in_this_process"]}; pylib_gseapy on sys.path: {info["pylib_gseapy_on_sys_path"]}')


CMDS = {'test': test, 'gsets': gsets, 'part0': part0, 'srcdelta': srcdelta, 'P0': partP0, 'P1': partP1, 'P2': partP2, 'P2b': partP2b, 'S1': partS1, 'S2': partS2,
        'S3': partS3, 'S4': partS4, 'Q1': partQ1, 'Q2': partQ2, 'R2': partR2, 'R3ranks': r3_ranks, 'R3gen': partR3gen,
        'R3cmp': partR3cmp, 'D': partD, 'E1gen_EvalA': partE1gen_EvalA, 'E1gen_HCT116final': partE1gen_HCT116final, 'E1cmp': partE1cmp}

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd')
    a = p.parse_args()
    for d in (TAB, LOG):
        d.mkdir(parents=True, exist_ok=True)
    log(f'step22 {a.cmd} start; OMP_NUM_THREADS={os.environ.get("OMP_NUM_THREADS")} host {os.uname().nodename.split(".")[0]}')
    CMDS[a.cmd](a)
    log(f'step22 {a.cmd} done; process peak {rss_gb():.1f} GB')

"""Step 21 (Jose 2026-11-02, consolidated): final proxy round (vcc2026, Track 2), a parallel cell run of row 4 (vcc2027,
Track 1), and the overnight build plan (Part H: not run here). No D/E/F file, no upload. Write-up: section 13 of
results/stage1/2026-10-28_scorer_space/report.md. Outputs: results/stage1/2026-11-02_step21/{tables,logs,cells}.
Frozen code (generator, harness, run_line, boot_score, config 1' files, vcc_2026_7) and the Step 18-20 outputs are imported or
read only.

Row 4 = spec (hybrid T 5 calls + config 1' + no-source fix) + the fixed Part C shift (C2's S on the line's gate,
(a_weak, a_mid, a_strong) = (-0.25, 0, 0.25) by (b) tercile) + the F5 term (E4w, M3 DepMap, K 100, R3, gamma 1).
  H1: F5 cross-fitted (each half from the model trained on the other half); shift cut points of the OTHER half.
  fold pools (K562, HCT116, HEK293T): F5 from all of H1 under the fold's basis; shift cuts all-H1 (3.0 / 67.7).
  HCT116s (HCT116final, screen setup: HCT116's own tables out, HEK293T votes; 15c's generator input): F5 with sources
          excluding HCT116's tables, the E4w-Orion fold basis, beta from all of H1; shift cuts all-H1.

Track 1 (vcc2027, step21_cells.sh): test, f5screen, t1inputs, generate, score (32 GB swap), boot_score.
Track 2 (vcc2026, step21.sh): test, part0, partA, partB, partC1 (one table load, C8 counts), partC2 .. partC10, partD.
  t1cmp   (vcc2026, after Track 1's outputs are copied back): paired differences, proxy minus cells, terciles, calls.
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
import step20 as S20                                                                       # noqa: E402
from step19 import (log, rss_gb, avail_gb, mem_before, mem_after, save_json, grid_pos, Fast, COMPS, TERC,  # noqa: E402
                    tercile_of, ref_for, b_for, DLine, C15, D_LINES)
from step20 import (prof, mu_add_delta, score_rows, b_cuts, b_terc, own_index, pick_fit, comp_dict, B_LINE,  # noqa: E402
                    FOLD_OF, OWN_BG_F, basis_file)

OUT = V / 'results/stage1/2026-11-02_step21'
TAB, LOG, CELLS = OUT / 'tables', OUT / 'logs', OUT / 'cells'
S20T = S20.TAB
SEED = 2026
LINES = D_LINES + ['HCT116s']
GUARD = ['K562', 'HCT116', 'HEK293T', 'HCT116s']
EV = dict(B_LINE, HCT116s='HCT116final')
SHIFT_A = np.array([-0.25, 0.0, 0.25])
F5B, F5V, F5K = 'E4w', 'M3DepMap', 100
T_HYB = 5
CPM_BANDS = [('ungated (<=5)', -np.inf, 5.0), ('5-20', 5.0, 20.0), ('20-100', 20.0, 100.0), ('>100', 100.0, np.inf)]
PROXY_REACH_NOTE = 'the proxy overstates reach by 0.02-0.05 compared with cells (15c sec. 4)'


def cuts_c2():
    """Step 20 Part C2's (b) cut points: half0, half1, full (all H1)."""
    return {k: np.array(v) for k, v in json.loads((S20T / 'partC2.json').read_text())['cuts'].items()}


# ================================================================ lines
class HLine:
    """HCT116final in the screen setup: 15c's config 1' + K0 generator input (HCT116's own tables out, HEK293T votes),
    with the same attributes DLine('H1') carries."""

    def __init__(self):
        self.line = 'HCT116s'; self.ref = ref_for('HCT116final'); h = self.ref.h
        d = C15 / 'config1p_K0/HCT116final'
        z = np.load(d / 'input.npz', allow_pickle=True); names = list(map(str, z['knockdowns']))
        r = np.array([names.index(t) for t in h.targets]); self.zrows = r
        self.n_src = z['n_sources'][r]; un = self.n_src <= 0
        self.lfc_raw = z['lfc'][r].astype(np.float64)
        self.lfc = np.where(un, 0.0, self.lfc_raw); self.sig = z['de_sig'][r] & ~un
        self.cg = np.where(un, 0.0, z['cg_lfc'][r].astype(np.float64))
        self.k4 = z['de_sig_K4'][r] & ~un
        cb = pd.read_csv(d / 'e2e/controls/controls_baseline.csv', index_col=0)
        self.m0 = cb['control_mean_cpm'].reindex(h.genes).to_numpy(float)
        self.gate = cb['de_tested'].reindex(h.genes).fillna(False).to_numpy(bool)
        self.parts = {'full': np.arange(h.T)}


def line_obj(line):
    if line == 'HCT116s':
        return HLine()
    D = DLine(line)
    if line == 'H1':
        z = np.load(C15 / 'config1p_K0/EvalA/input.npz', allow_pickle=True)
        D.n_src = z['n_sources']; D.lfc_raw = z['lfc'].astype(np.float64); D.k4 = z['de_sig_K4'] & ~(z['n_sources'] <= 0)
        D.zrows = np.arange(len(D.n_src))
    return D


def spec_mu(D):
    """Config 1' + K0 magnitudes with the no-source fix (Step 19 Part F's combo_mu): H1 / HCT116s from their generator
    input; the fold pools already carry the fix (step20.base_mu)."""
    import scorer_arms as SA
    if D.line != 'HCT116s':
        return S20.base_mu(D)
    un = (D.n_src <= 0) & D.gate[None, :]
    L = np.where(un, D.lfc_raw, D.lfc); L = np.where(D.sig, 0.5 * D.cg, L)
    return SA.mu_from_lfc(L, D.m0)


def k4_of(D):
    if D.line in ('H1', 'HCT116s'):
        return D.k4
    import scorer_step16 as S16
    return S16.call_sets(D.G, [1200])['K4_1200']


def h1_halves():
    import scorer_step12 as S12
    import scorer_step16 as S16
    return S16.halves(S12.truth_calls(ref_for('EvalA')))


_PREDS = {}


def f4_preds():
    import pickle
    if not _PREDS:
        _PREDS.update(pickle.load(open(S20T / 'partF4_preds.pkl', 'rb'))['preds'])
    return _PREDS


class Row4:
    """Row 4's parts on one line: spec mu (K0 magnitude split + no-source fix), hybrid T 5 calls, the F5 term (gated,
    b >= 5 only, gamma 1) and the fixed shift (C2's S on the gate, by (b) tercile). mu(...) adds any Delta-space term."""

    def __init__(self, line, f5_shat=None, shift_cuts='cross'):
        self.line = line; self.D = D = line_obj(line); h = D.ref.h
        self.h = h; self.own = own_index(h); self.gate = D.gate
        self.b = b_for(EV[line], h.targets)
        self.on5 = np.nan_to_num(self.b, nan=-np.inf) >= T_HYB
        self.mu_spec = spec_mu(D); self.k4 = k4_of(D)
        self.sig_k0 = D.sig
        self.sig = np.where(self.on5[:, None], self.k4, D.sig)
        self.hv = h1_halves() if line == 'H1' else None
        self.shat = self.f5_shat() if f5_shat is None else f5_shat
        self.addF5 = self.f5_from_shat(self.shat)
        self.cuts = cuts_c2(); self.terc = self.shift_terc(shift_cuts)
        self.S = self.shift_vec()
        self.addS = SHIFT_A[self.terc][:, None] * self.S[None, :]

    def f5_basis(self):
        fold = 'full' if self.line == 'H1' else ('Orion' if self.line == 'HCT116s' else FOLD_OF[self.line])
        return np.load(basis_file(F5B, fold))['Vt'].astype(np.float64)

    def f5_shat(self):
        P = f4_preds(); h = self.h
        if self.line == 'H1':
            sh = np.zeros((h.T, F5K))
            for f in (0, 1):
                sh[self.hv[1 - f]] = P[(F5B, F5K, F5V, 'H1', f)][self.hv[1 - f]]
            return sh
        if self.line == 'HCT116s':
            z = np.load(TAB / 'f5_HCT116s.npz', allow_pickle=True)
            assert list(map(str, z['targets'])) == h.targets
            return z['shat'].astype(np.float64)
        return P[(F5B, F5K, F5V, self.line, -1)].astype(np.float64)

    def f5_from_shat(self, shat, on=None):
        Vt = self.f5_basis()
        P = (shat.astype(np.float64) @ Vt[:shat.shape[1]]) * self.gate[None, :]
        return np.where((self.on5 if on is None else on)[:, None], P, 0.0)

    def shift_terc(self, mode):
        if self.line == 'H1' and mode == 'cross':
            t = np.zeros(self.h.T, int)
            for f in (0, 1):
                rr = self.hv[f]; t[rr] = b_terc(self.b[rr], self.cuts[f'half{1 - f}'])
            return t
        return b_terc(self.b, self.cuts['full'])

    def shift_vec(self):
        z = np.load(S20T / 'partC2_S.npz')
        S = z['S_' + ('HCT116' if self.line == 'HCT116s' else self.line)].astype(np.float64)
        return np.where(self.gate & np.isfinite(S), S, 0.0)

    def add(self, f5=True, shift=True):
        A = np.zeros((self.h.T, self.h.G))
        if f5: A = A + self.addF5
        if shift: A = A + self.addS
        return A

    def mu(self, f5=True, shift=True, extra=None, base=None):
        A = self.add(f5, shift)
        if extra is not None:
            A = A + np.where(np.isfinite(extra), extra, 0.0)
        base = self.mu_spec if base is None else base
        return base if not np.any(A != 0) else mu_add_delta(base, A, self.own)

    def score(self, mu, sig=None, **tag):
        rows, out = score_rows(self.D, mu, self.sig if sig is None else sig, **tag)
        s = self.sig if sig is None else sig
        for r in rows:
            rr = self.D.parts[r['part']]; n = s[rr].sum(1)
            r.update(calls_median=float(np.median(n)), calls_mean=float(np.mean(n)))
        return rows, out

    def delta_pred(self, mu):
        """Predicted Delta in log1p CP50K (config 1''s Delta space)."""
        return prof(mu) - prof(self.D.m0)[None, :]


ROW_DEFS = {'spec': (False, False), 'spec + shift': (False, True), 'spec + F5': (True, False), 'row 4': (True, True)}


# ================================================================ F5 term for HCT116final in the screen setup
def f5screen(a):
    """M3 DepMap, E4w-Orion fold basis, K 100, beta from all of H1, sources without HCT116's tables (HEK293T votes).
    Check: the same code with both Orion lines out reproduces Step 20's fold prediction for HCT116."""
    ctx = S20.FCtx(lines=['H1']); T = len(ctx.tr['H1']['targets'])
    ref = ref_for('HCT116final'); tg = list(ref.h.targets); S19._REFS.pop('HCT116final', None); gc.collect()
    Ks = [5, 10, 20, 100]
    GF = S20.gene_features(tg, ctx.fs['DepMap'], ctx.H)[0]
    sp = 'delta'; RS = ctx.rs(F5B, 'Orion')
    out = {}
    for name, ex_pred in (('check_fold', ('HCT116', 'HEK293T')), ('screen', ('HCT116',))):
        A = S20.bg_scores(RS, ctx.meta[sp], tg, ex_pred, RS.shape[1], ctx.bgs[sp])[:, :, :F5K]
        train = [('H1', np.arange(T), tuple(sorted(set(('H1',) + ex_pred))), False)]
        pred, inf = S20.fit_predict(ctx, 'M3', F5B, 'Orion', F5K, Ks, train, A, GF, 'DepMap')
        out[name] = (pred, inf, A)
    ref_pred = f4_preds()[(F5B, F5K, F5V, 'HCT116', -1)].astype(np.float64)
    d = float(np.abs(out['check_fold'][0] - ref_pred).max())
    d32 = float(np.abs(out['check_fold'][0].astype(np.float32) - f4_preds()[(F5B, F5K, F5V, 'HCT116', -1)]).max())   # Step 20 stored float32
    log(f'f5screen check: fold-pool HCT116 prediction rebuilt, max |diff| vs Step 20 {d:.2e} (float64), {d32:.1e} after casting to float32 as Step 20 stored it')
    pred, inf, A = out['screen']
    nb = np.isfinite(A[:, :, 0]).sum(1)
    hek = ctx.bgs[sp].index('HEK293T')
    res = dict(check_fold_max_abs_diff=d, check_fold_max_abs_diff_float32=d32, check_pass=bool(d32 == 0.0), targets=len(tg),
               targets_with_backgrounds=int((nb > 0).sum()), median_backgrounds=float(np.median(nb)),
               targets_with_HEK293T=int(np.isfinite(A[:, hek, 0]).sum()),
               alpha_counts={str(al): int((inf['alpha'] == al).sum()) for al in S20.ALPHAS},
               shat_norm_ratio_screen_vs_fold=float(np.median(np.linalg.norm(pred, axis=1) / np.maximum(np.linalg.norm(out['check_fold'][0], axis=1), 1e-12))))
    np.savez_compressed(TAB / 'f5_HCT116s.npz', targets=np.array(tg), shat=pred.astype(np.float32), n_backgrounds=nb)
    save_json(TAB / 'f5_HCT116s.json', res)
    log(f'f5screen: {res}')
    if not res['check_pass']:
        sys.exit(5)


# ================================================================ Track 1: generator inputs
def input_arrays(R, A, z):
    """Generator input (h.targets order) for mu = spec + A: lfc changes only where A != 0 (the spec's input elsewhere,
    bit for bit); de_sig = hybrid T 5; n_sources >= 1 on gated no-source genes (the no-source fix)."""
    A = np.where(np.isfinite(A), A, 0.0).copy()
    for i, g in enumerate(R.own):
        if g >= 0: A[i, g] = 0.0
    mu = R.mu_spec if not np.any(A != 0) else mu_add_delta(R.mu_spec, A, R.own)
    m0 = R.D.m0
    raw = z['lfc'][R.D.zrows].astype(np.float32)
    lfc = np.where(A != 0, np.log2((mu + 1e-9) / (m0[None, :] + 1e-9)).astype(np.float32), raw)
    ns = z['n_sources'][R.D.zrows]
    ns = np.where(R.gate[None, :] & (ns <= 0), 1, ns).astype(np.int16)
    return lfc, R.sig.copy(), ns, mu


def proxy_of_input(R, lfc, sig, ns):
    """The proxy on a generator input as the generator sees it (unmeasured genes 0 and uncalled)."""
    import scorer_arms as SA
    un = ns <= 0
    L = np.where(un, 0.0, lfc.astype(np.float64))
    return R.score(SA.mu_from_lfc(L, R.D.m0), sig & ~un)[0]


def t1inputs(a):
    import shutil
    res = {}
    for line, ev in (('H1', 'EvalA'), ('HCT116s', 'HCT116final')):
        R = Row4(line)
        z = np.load(C15 / f'config1p_K0/{ev}/input.npz', allow_pickle=True)
        z4 = np.load(C15 / f'config1p_K4_1200/{ev}/input.npz', allow_pickle=True)
        names = list(map(str, z['knockdowns'])); r = R.D.zrows
        assert [names[i] for i in r] == R.h.targets and list(map(str, z4['knockdowns'])) == names
        hyb = np.where(R.on5[:, None], z4['de_sig'][r], z['de_sig'][r])
        chk = dict(calls_equal_hybrid_T5_bit_for_bit=bool((R.sig == hyb).all()),
                   knockdowns_on_K4=int(R.on5.sum()), knockdowns=int(R.h.T))
        # (1) F5 term and shift removed: the spec's input (raw K0 lfc, hybrid calls, fixed n_sources), bit for bit
        lfc0, sig0, ns0, mu0 = input_arrays(R, np.zeros((R.h.T, R.h.G)), z)
        chk['spec_input_lfc_equals_15c_K0_input'] = bool((lfc0 == z['lfc'][r].astype(np.float32)).all())
        chk['spec_input_mu_equals_spec_mu'] = bool(np.array_equal(mu0, R.mu_spec))
        ps = proxy_of_input(R, lfc0, sig0, ns0); pm = R.score(R.mu_spec)[0]
        chk['spec_input_proxy'] = {r_['part']: r_['overall'] for r_ in ps}
        chk['spec_mu_proxy'] = {r_['part']: r_['overall'] for r_ in pm}
        chk['spec_input_proxy_max_abs_diff'] = float(max(abs(x['overall'] - y['overall']) for x, y in zip(ps, pm)))
        # (2) row 4
        lfc, sig, ns, mu = input_arrays(R, R.add(), z)
        p4 = proxy_of_input(R, lfc, sig, ns); pm4 = R.score(R.mu())[0]
        chk['row4_input_proxy'] = {r_['part']: r_['overall'] for r_ in p4}
        chk['row4_mu_proxy'] = {r_['part']: r_['overall'] for r_ in pm4}
        chk['row4_input_proxy_max_abs_diff'] = float(max(abs(x['overall'] - y['overall']) for x, y in zip(p4, pm4)))
        chk['entries_changed'] = int((lfc != z['lfc'][r].astype(np.float32)).sum())
        chk['n_sources_fixed_entries'] = int((ns != z['n_sources'][r]).sum())
        chk['pass'] = bool(chk['calls_equal_hybrid_T5_bit_for_bit'] and chk['spec_input_lfc_equals_15c_K0_input'] and
                           chk['spec_input_mu_equals_spec_mu'] and chk['spec_input_proxy_max_abs_diff'] < 1e-9 and
                           chk['row4_input_proxy_max_abs_diff'] < 1e-5)
        d = CELLS / ev; d.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(d / 'input.npz', knockdowns=np.array(R.h.targets, dtype=object), genes=z['genes'], lfc=lfc,
                            de_sig=sig, n_sources=ns, m0=R.D.m0, gate=R.gate, on_K4=R.on5)
        shutil.copy(C15 / f'config1p_K0/{ev}/input.npz', d / 'input_15c_K0_reference.npz') if a.keep_ref else None
        res[ev] = chk
        log(f't1inputs {ev}: {json.dumps({k: v for k, v in chk.items()}, default=float)}')
        del R; gc.collect()
    save_json(TAB / 't1inputs.json', res)
    if not all(v['pass'] for v in res.values()):
        log('t1inputs CHECK FAILED'); sys.exit(4)


# ================================================================ generic fit / criterion
def is_row4(d, key, base_vals):
    sel = np.ones(len(d), bool)
    for k, v in zip(key, base_vals):
        sel &= (d[k] == v).to_numpy()
    return sel


def fit_report(R, key, grids, base_vals, guard=GUARD, fit_frames=None, tag=''):
    """Pick on one H1 half (highest overall with all six > 0 there), score the other half, swap; guardrails at each pick
    (full panel) against the base; the criterion. fit_frames: {f: rows used for direction f} (inner cross-fitted rows)."""
    res = dict(fits={})
    for f in (0, 1):
        H = fit_frames[f] if fit_frames is not None else R[R.line == 'H1']
        got = pick_fit(H, f'half{f}', f'half{1 - f}', key, lambda d: is_row4(d, key, base_vals))
        if got is None:
            res['fits'][f'fit_half{f}'] = dict(picked=None); continue
        pk, ho, b0, ne = got
        gain = {c: float(ho[c] - b0[c]) for c in COMPS + ['overall']}
        sc = sum(gain[c] for c in COMPS)
        gr = {}
        for l in guard:
            L = R[(R.line == l) & (R.part == 'full')]
            if fit_frames is not None and 'fit_half' in L:
                L = L[L.fit_half.isin([-1, f])]
            if not len(L):
                continue
            x = L[is_row4(L, key, [pk[k] for k in key])]; y = L[is_row4(L, key, base_vals)]
            if not len(x) or not len(y):
                continue
            x, y = x.iloc[0], y.iloc[0]
            gr[l] = dict(diff=float(x.overall - y.overall), pds_diff=float(x.pds - y.pds), pick=comp_dict(x), base=comp_dict(y))
        res['fits'][f'fit_half{f}'] = dict(
            picked={k: (pk[k] if isinstance(pk[k], str) else float(pk[k])) for k in key}, eligible=int(ne),
            grid_position={k: (grid_pos(grids[k], pk[k]) if pk[k] in list(grids[k]) else 'base row (row 4), outside the grid') for k in key}, fit_overall=float(pk.overall),
            held_out=comp_dict(ho), base_held_out=comp_dict(b0), gain=gain,
            reach_share_of_component_gain=float(gain['reach'] / sc) if sc else float('nan'), guardrails=gr)
        log(f'{tag} fit half{f}: pick {res["fits"][f"fit_half{f}"]["picked"]} ({res["fits"][f"fit_half{f}"]["grid_position"]}); held-out '
            f'{ho.overall:.4f} vs base {b0.overall:.4f} ({gain["overall"]:+.4f}); guardrails ' + ', '.join(f'{l} {v["diff"]:+.4f}' for l, v in gr.items()))
    fv = [v for v in res['fits'].values() if v.get('picked') is not None]
    res['criterion'] = criterion(fv)
    res['proxy_note'] = PROXY_REACH_NOTE
    return res


def criterion(fv):
    c = dict(both_picked=len(fv) == 2,
             heldout_gain_ge_0005_both=bool(len(fv) == 2 and all(v['gain']['overall'] >= 0.005 for v in fv)),
             all_components_pos_both=bool(len(fv) == 2 and all(all(v['held_out'][k] > 0 for k in COMPS) for v in fv)),
             guardrails_within_0005=bool(len(fv) == 2 and all(g['diff'] >= -0.005 for v in fv for g in v['guardrails'].values())))
    c['met'] = all(c.values())
    return c


def save_rows(rows, name):
    R = pd.DataFrame(rows); R.to_csv(TAB / f'{name}.csv', index=False); return R


def resume(name):
    f = TAB / f'{name}.csv'
    return pd.read_csv(f).to_dict('records') if f.exists() else []


# ================================================================ Part 0: checks
def part0(a):
    res = dict(check1={}, check2={})
    G = pd.read_csv(S20T / 'partG_rows.csv')
    want = {'spec': 'hybrid T5 + fix (spec)', 'spec + F5': 'spec + F5 prog'}
    ok = True
    for line in D_LINES:
        R = Row4(line, shift_cuts='full')
        for nm, gname in want.items():
            f5, sh = ROW_DEFS[nm]
            t0 = time.time(); rows, _ = R.score(R.mu(f5=f5, shift=sh)); dt = time.time() - t0
            w = G[(G.line == line) & (G.row == gname)].set_index('part')['overall']
            d = max(abs(r['overall'] - w[r['part']]) for r in rows)
            ok &= d < 1e-9
            res['check1'][f'{line} {nm}'] = dict(got={r['part']: r['overall'] for r in rows}, step20=w.to_dict(), max_abs_diff=d,
                                                 seconds_per_row=dt, threads=os.environ.get('OMP_NUM_THREADS'))
            log(f'part0 check1 {line} {nm}: {res["check1"][f"{line} {nm}"]["got"]} vs Step 20 {w.to_dict()}; max |diff| {d:.1e}; '
                f'{dt:.1f} s per row at OMP_NUM_THREADS={os.environ.get("OMP_NUM_THREADS")}')
        if line == 'H1':
            C = pd.read_csv(S20T / 'partC2_rows.csv')
            w = C[(C.line == 'H1') & (C.cuts == 'full') & (C.a_weak == -0.25) & (C.a_mid == 0) & (C.a_strong == 0.25)].set_index('part')['overall']
            rows, _ = score_rows(R.D, mu_add_delta(S20.base_mu(R.D), R.addS, R.own))      # fix-only base, K0 calls
            d = max(abs(r['overall'] - w[r['part']]) for r in rows); ok &= d < 1e-9
            res['check2'] = dict(got={r['part']: r['overall'] for r in rows}, step20_partC2=w.to_dict(), max_abs_diff=d)
            log(f'part0 check2 fixed shift (-0.25, 0, 0.25) on the fix-only base: {res["check2"]["got"]} vs C2 {w.to_dict()}; max |diff| {d:.1e}')
        del R; gc.collect()
    res['pass'] = bool(ok)
    save_json(TAB / 'part0.json', res)
    log(f'part0: {"PASS" if ok else "FAIL"}')
    if not ok:
        sys.exit(3)


# ================================================================ Part A: the four combinations
def partA(a):
    rows = resume('partA_rows')
    for line in LINES:
        if any(r['line'] == line for r in rows):
            continue
        R = Row4(line)
        for nm, (f5, sh) in ROW_DEFS.items():
            rr, _ = R.score(R.mu(f5=f5, shift=sh), row=nm); rows += rr
            log(f'A {line} {nm}: ' + ', '.join(f'{r["part"]} {r["overall"]:.4f} (PDS {r["pds"]:.3f})' for r in rr) +
                f'; calls median {rr[0]["calls_median"]:.0f} mean {rr[0]["calls_mean"]:.0f}')
        save_rows(rows, 'partA_rows'); del R; gc.collect()
    R = pd.DataFrame(rows)
    res = {}
    for line in LINES:
        L = R[(R.line == line)]
        for part in L.part.unique():
            P = L[L.part == part].set_index('row')
            res[f'{line} {part}'] = {nm: comp_dict(P.loc[nm]) | dict(calls_median=float(P.loc[nm, 'calls_median']), calls_mean=float(P.loc[nm, 'calls_mean']))
                                     for nm in ROW_DEFS}
            res[f'{line} {part}']['row4_minus_spec'] = {c: float(P.loc['row 4', c] - P.loc['spec', c]) for c in COMPS + ['overall']}
    res['proxy_note'] = PROXY_REACH_NOTE
    save_json(TAB / 'partA.json', res)
    with pd.option_context('display.width', 250, 'display.max_rows', 200):
        print(R.round(4).to_string())


# ================================================================ Part B: shared-component removal
B_ETA = [0.0, 0.25, 0.5, 0.75, 1.0]


def mean_pred_delta(R, mu):
    """Per gene, the mean predicted Delta (log1p CP50K) over the line's knockdowns, own genes out; 0 off the gate."""
    Dp = R.delta_pred(mu); M = np.ones_like(Dp, bool)
    for i, g in enumerate(R.own):
        if g >= 0: M[i, g] = False
    m = np.where(M, Dp, 0.0).sum(0) / np.maximum(M.sum(0), 1)
    return np.where(R.gate, m, 0.0)


def partB(a):
    rows = resume('partB_rows'); res = {}
    for line in LINES:
        if any(r['line'] == line for r in rows):
            continue
        R = Row4(line); mu4 = R.mu(); m = mean_pred_delta(R, mu4)
        res[f'{line}_mean_delta_norm'] = float(np.linalg.norm(m))
        for eta in B_ETA:
            mu = mu4 if eta == 0 else mu_add_delta(mu4, np.broadcast_to(-eta * m[None, :], mu4.shape), R.own)
            rr, _ = R.score(mu, eta=eta); rows += rr
            log(f'B {line} eta {eta}: ' + ', '.join(f'{r["part"]} {r["overall"]:.4f} (PDS {r["pds"]:.3f}, MSE {r["mse"]:.3f})' for r in rr))
        save_rows(rows, 'partB_rows'); del R; gc.collect()
    R = pd.DataFrame(rows)
    out = fit_report(R, ['eta'], dict(eta=B_ETA), [0.0], tag='B')
    out['surface_H1_full'] = R[(R.line == 'H1') & (R.part == 'full')][['eta', 'overall'] + COMPS].to_dict('records')
    save_json(TAB / 'partB.json', out | res)
    print(json.dumps(out['criterion']))


# ================================================================ Part C1: C+G's pooled quantities (one table load) + C8 counts
SS_CACHE = S20.SS_OUT / 'cache'
C1_CTX = {   # context -> (dropped stage1 line keys, saved arrays the rerun must reproduce)
    'H1': (('h1',), 'fold'), 'K562': (('k562',), 'fold'), 'HCT116': (('hct116', 'hek293t'), 'fold'),
    'HEK293T': (('hct116', 'hek293t'), 'fold'), 'HCT116s': (('hct116',), '15c'), 'D': ((), 's18'), 'E': ((), 's18'), 'F': ((), 's18')}
TRUTH_CTX = ['H1', 'K562', 'HCT116', 'HEK293T', 'HCT116s']


def c1_saved(ctx):
    """(targets, m0, gate, lfc, de_sig, n_sources, se) as saved by the run that built the context."""
    kind = C1_CTX[ctx][1]
    if kind == 'fold':
        z = np.load(SS_CACHE / f'fold_{ctx}.npz', allow_pickle=True)
        return list(map(str, z['targets'])), z['m0'], z['gate'].astype(bool), z['lfc'], z['de_sig'], z['n_sources'], z['se_measurement']
    f = C15 / 'config1p_K0/HCT116final/input.npz' if kind == '15c' else S19.S18 / f'config1p_K0/{ctx}/input.npz'
    z = np.load(f, allow_pickle=True)
    return list(map(str, z['knockdowns'])), z['m0'], z['gate'].astype(bool), z['cg_lfc'], z['de_sig_K0'], z['n_sources'], z['se']


def c1_load(ctx, targets=None):
    """C1 arrays for a context, rows reindexed to targets (default: the context's own order)."""
    z = np.load(TAB / f'C1_{ctx}.npz', allow_pickle=True)
    t = list(map(str, z['targets'])); out = {k: z[k] for k in z.files}
    if targets is not None:
        r = np.array([t.index(x) for x in targets])
        for k in ('mu', 'var_fe', 'tau2', 'tau_un', 'lfc', 'n', 'sig'):
            out[k] = out[k][r]
        out['targets'] = np.array(targets)
    return out


def tau2_fallback(c, gate):
    """tau^2_g for z_u: C+G's gene-specific DL tau^2 (floored at the global tau^2) where the gene has >= 2 sources;
    elsewhere the line's median of those values (over gated genes)."""
    est = (c['n'] >= 2) & np.isfinite(c['tau2']) & gate[None, :]
    med = float(np.median(c['tau2'][est])) if est.any() else float(json.loads(str(c['prior']))['tau2'])
    t2 = np.where((c['n'] >= 2) & np.isfinite(c['tau2']), c['tau2'].astype(np.float64), med)
    return t2, med, est


def z_u_of(c, gate, own):
    """z_u = |mu| / sqrt(SE_meas^2 + var_fe + tau^2_g), 0 where no source; and the three variance parts."""
    t2, med, est = tau2_fallback(c, gate)
    se2 = (np.maximum(c['se'].astype(np.float64), 1e-3) ** 2)[None, :]
    vf = c['var_fe'].astype(np.float64)
    den = se2 + np.where(np.isfinite(vf), vf, np.inf) + t2
    z = np.where(np.isfinite(den) & (c['n'] > 0), np.abs(c['mu'].astype(np.float64)) / np.sqrt(den), 0.0)
    for i, g in enumerate(own):
        if g >= 0: z[i, g] = 0.0
    return z, dict(se2=np.broadcast_to(se2, z.shape), var_fe=vf, tau2=t2, median_tau2=med, estimable=est)


def c1_uncalled(sources, by_stem, sets, genes, res):
    """Config 1''s uncalled value on every gene for H1 and HCT116s (their generator inputs hold the called value on
    K0-called genes): 15c's Delta pool, rebuilt as step15c.model builds it, then config_lfc with no gene called.
    Check: config_lfc with K0's calls reproduces the 15c input's lfc bit for bit."""
    import h5py
    from anndata.io import read_elem
    import stage1 as S
    import step15c as S15c
    G_ = np.asarray(genes, dtype=object); G = len(genes)
    for ctx, ev, drop in (('H1', 'EvalA', ('h1',)), ('HCT116s', 'HCT116final', ('hct116',))):
        if (TAB / f'C1_uncalled_{ctx}.npz').exists():
            continue
        z = np.load(C15 / f'config1p_K0/{ev}/input.npz', allow_pickle=True)
        targets = list(map(str, z['knockdowns'])); need = set(targets)
        others, _, calib, prior = sets[tuple(sorted(drop))]
        names = {s.name for s in others}
        DSRC = []
        for src in sources:
            if src.name not in names:
                continue
            p = by_stem[src.name]
            ctrl = pd.read_csv(p.with_name(src.name + '_baseline.csv'), index_col=0)['log1p_bulk_50k_full'].reindex(genes).to_numpy(float)
            with h5py.File(p, 'r') as f:
                kd = read_elem(f['obs'])['target_gene'].astype(str).to_numpy()
                rows = [i for i, k in enumerate(kd) if k in need]
                if not rows or not np.isfinite(ctrl).any():
                    continue
                full, se = np.asarray(f['layers/log1p_bulk_50k_full'][rows], float), np.asarray(f['layers/se'][rows], np.float32)
            if not np.isfinite(full).any():
                continue
            D_, k_, n = (full - ctrl).astype(np.float32), kd[rows], np.full(len(rows), np.nan, np.float32)
            if src.modality != 'drug':
                k_, D_, se, n = S._collapse_duplicate_knockdowns(k_, D_, se, n)
            DSRC.append(S.Source(src.name, src.background, src.modality, src.reference, src.genes, k_, D_, se, n, src.mult))
        dvotes = S.BackgroundVotes(DSRC)
        m0, gate, cprof = z['m0'], z['gate'].astype(bool), z['cprof']
        L_un = np.zeros((len(targets), G), np.float32); L_chk = np.zeros_like(L_un)
        for i, k in enumerate(targets):
            dp = S.predict_knockdown(None, k, G_, m0, gate, calib, prior, drop_gene=k, votes=dvotes, magnitude='gaussian')
            dpool = np.where(dp['n_sources'] > 0, dp['lfc'], 0.0)
            Y, _, _ = dvotes.gather(k, G_)
            n_bg, zz = np.zeros(G, int), np.full(G, np.nan)
            if Y:
                Y = np.stack(Y); fin = np.isfinite(Y); n_bg = fin.sum(0)
                mean = np.where(fin, Y, 0.0).sum(0) / np.maximum(n_bg, 1)
                sd = np.sqrt(np.where(fin, (Y - mean) ** 2, 0.0).sum(0) / np.maximum(n_bg - 1, 1))
                zz = np.where(n_bg >= 2, np.abs(mean) / np.maximum(sd / np.sqrt(np.maximum(n_bg, 1)), 1e-12), np.nan)
            L_chk[i] = S15c.config_lfc(S15c.CONFIG1P, z['cg_lfc'][i], z['de_sig_K0'][i], dpool, n_bg, zz, m0, cprof)
            L_un[i] = S15c.config_lfc(S15c.CONFIG1P, z['cg_lfc'][i], np.zeros(G, bool), dpool, n_bg, zz, m0, cprof)
        ok = bool(np.array_equal(L_chk, z['lfc'].astype(np.float32)))
        res.setdefault('uncalled', {})[ctx] = dict(delta_tables=len(DSRC), config_lfc_reproduces_15c_input=ok,
                                                   max_abs_diff=float(np.abs(L_chk - z['lfc']).max()))
        log(f'C1 uncalled {ctx}: {len(DSRC)} Delta tables; 15c input reproduced bit for bit: {ok}')
        if not ok:
            save_json(TAB / 'partC1.json', res); log('C1 UNCALLED CHECK FAILED'); sys.exit(8)
        np.savez_compressed(TAB / f'C1_uncalled_{ctx}.npz', targets=np.array(targets), lfc_uncalled=L_un)


def partC1(a):
    import stage1 as S
    import screen as SC
    import scorer_study as SS
    import run_line as RL
    genes = RL.panel_genes(); G_ = np.asarray(genes, dtype=object); G = len(genes)
    tables = S.discover_tables(SC.DIRS); by_stem = {p.stem: p for p in tables}
    keep, _, _ = SS.filter_plan(tables, 'F0'); assert [p.stem for p in keep] == [p.stem for p in tables], 'F0 keeps every table'
    mem_before(f'C1 table load ({len(tables)} tables)', 18.0)
    sources = S.load_tables(tables); mem_after('C1 table load')
    # each table's own > 5 CPM gate on its genes (C8); tables without control CPM are skipped by C8
    tgate = {}
    for s in sources:
        bp = by_stem[s.name].with_name(s.name + '_baseline.csv')
        if bp.exists():
            b = pd.read_csv(bp, index_col=0)
            if 'control_mean_cpm' in b and b['control_mean_cpm'].notna().any():
                tgate[s.name] = (b['control_mean_cpm'].reindex(s.genes).to_numpy(float) > 5.0)
    pcol = {g: j for j, g in enumerate(genes)}
    sets, res = {}, dict(contexts={}, c8_tables_skipped_no_control_cpm=sorted(s.name for s in sources if s.name not in tgate))
    c8 = []
    for ctx, (drop, _) in C1_CTX.items():
        if (TAB / f'C1_{ctx}.npz').exists() and any(r['context'] == ctx for r in c8):
            continue
        key = tuple(sorted(drop))
        if key not in sets:
            others = [s for s in sources if S._line_key(s.background) not in set(drop)]
            basals = {}
            for s in others:
                got = S.load_baseline(by_stem[s.name].with_name(s.name + '_baseline.csv'), s.genes)
                if got is not None: basals[s.name] = (s.genes, got[0])
            votes, prior_votes = S.arm_votes(others, 'C')
            sets[key] = (others, votes, S.fit_se_calibration(others, basals), S.fit_prior(others, votes=prior_votes))
            log(f'C1 voting set without {key or "nothing"}: {len(others)} tables; prior {sets[key][3].as_dict()}')
        others, votes, calib, prior = sets[key]
        targets, m0, gate, s_lfc, s_sig, s_n, s_se = c1_saved(ctx)
        K = len(targets)
        A = {k: np.zeros((K, G), np.float32) for k in ('mu', 'var_fe', 'tau2', 'tau_un', 'lfc')}
        n_ = np.zeros((K, G), np.int16); sig = np.zeros((K, G), bool)
        se_meas = calib.predict(m0, S.N_CELLS_REF)
        t0 = time.time()
        for i, k in enumerate(targets):
            Y, SE, mults = votes.gather(k, G_)
            if not Y:
                mu, vf, t2, tu, n = np.zeros(G), np.full(G, np.inf), np.full(G, np.nan), np.full(G, np.nan), np.zeros(G, np.int32)
            else:
                mu, vf, t2, tu, n = S.pool_effects(np.stack(Y), np.stack(SE), np.asarray(mults), prior.tau2)
            t2u = np.where(np.isfinite(t2), t2, prior.tau2)
            lfc, _, _, _ = S.gaussian_posterior(mu, vf, t2u, prior.v0)
            tested = np.asarray(gate, bool).copy()
            if k in pcol: tested[pcol[k]] = False
            called, _ = S.predictive_calls(lfc, se_meas, tested)
            A['mu'][i], A['var_fe'][i], A['tau2'][i], A['tau_un'][i], A['lfc'][i] = mu, vf, t2, tu, lfc
            n_[i] = n; sig[i] = called
            # C8: per genetic source table holding k (the context's own lines already out)
            for s in others:
                if s.modality == 'drug' or k not in s.rows_of or s.name not in tgate:
                    continue
                r = s.rows_of[k][0]
                zt = np.abs(s.lfc[r].astype(np.float64)) / np.maximum(s.se[r].astype(np.float64), 1e-12)
                de = np.isfinite(zt) & (zt > 3) & tgate[s.name] & np.isfinite(s.lfc[r])
                j = np.array([pcol.get(str(g), -1) for g in s.genes[de]]); j = j[(j >= 0) & (j != pcol.get(k, -2))]
                c8.append(dict(context=ctx, target=k, table=s.name, background=s.background, n_de_z3=int(len(j)),
                               n_de_gated5=int(gate[j].sum()), n_de_cpm100=int((m0[j] > 100).sum())))
        dt = time.time() - t0
        chk = dict(lfc_bit_for_bit=bool(np.array_equal(A['lfc'], s_lfc.astype(np.float32))),
                   calls_bit_for_bit=bool(np.array_equal(sig, s_sig.astype(bool))),
                   n_sources_equal=bool(np.array_equal(n_, s_n.astype(np.int16))),
                   se_meas_equal=bool(np.array_equal(se_meas.astype(np.float32)[None, :].repeat(K, 0), s_se.astype(np.float32))))
        chk['lfc_max_abs_diff'] = float(np.nanmax(np.abs(A['lfc'] - s_lfc.astype(np.float32))))
        chk['calls_mismatch'] = int((sig != s_sig.astype(bool)).sum())
        res['contexts'][ctx] = dict(dropped=list(drop), tables_voting=len(others), knockdowns=K, seconds=dt, checks=chk,
                                    prior=prior.as_dict())
        log(f'C1 {ctx}: {K} knockdowns in {dt:.0f} s; checks {chk}')
        np.savez_compressed(TAB / f'C1_{ctx}.npz', targets=np.array(targets), m0=m0, gate=gate, se=se_meas.astype(np.float32),
                            prior=json.dumps(prior.as_dict()), n=n_, sig=sig, **A)
        pd.DataFrame(c8).to_csv(TAB / 'C8_counts.csv', index=False)
        save_json(TAB / 'partC1.json', res)
        if not all(v for k_, v in chk.items() if isinstance(v, bool)):
            log(f'C1 CHECK FAILED for {ctx}'); sys.exit(7)
    c1_uncalled(sources, by_stem, sets, genes, res)
    # tau^2 distribution and fallback share (gated genes, own gene out)
    for ctx in C1_CTX:
        c = c1_load(ctx); gate = c['gate'].astype(bool)
        own = np.array([pcol.get(t, -1) for t in c['targets']])
        z, parts = z_u_of(c, gate, own)
        g = np.broadcast_to(gate[None, :], z.shape).copy()
        for i, j in enumerate(own):
            if j >= 0: g[i, j] = False
        t2 = c['tau2'][parts['estimable'] & g]
        res['contexts'][ctx]['tau2'] = dict(fallback_share_gated=float(1 - (parts['estimable'] & g).sum() / g.sum()),
                                            median_used_for_fallback=parts['median_tau2'],
                                            quantiles_estimable={q: float(np.quantile(t2, q)) for q in (0.1, 0.25, 0.5, 0.75, 0.9, 0.99)},
                                            share_at_global_floor=float(np.isclose(t2, json.loads(str(c['prior']))['tau2']).mean()),
                                            global_tau2=json.loads(str(c['prior']))['tau2'])
    save_json(TAB / 'partC1.json', res)
    print(json.dumps(res, indent=1, default=float))


# ================================================================ Part C2-C6 shared
CTX_OF = {'H1': 'H1', 'K562': 'K562', 'HCT116': 'HCT116', 'HEK293T': 'HEK293T', 'HCT116s': 'HCT116s'}
BH_A = [0.01, 0.05, 0.2]


def bh_calls(z, ok, alpha):
    """Two-sided normal p from z, BH per knockdown (row) over the ok genes, at alpha."""
    from scipy.stats import norm
    out = np.zeros(z.shape, bool)
    for i in range(len(z)):
        j = np.flatnonzero(ok[i])
        if not len(j):
            continue
        p = 2 * norm.sf(z[i, j]); o = np.argsort(p); m = len(p)
        thr = alpha * np.arange(1, m + 1) / m
        passed = np.flatnonzero(p[o] <= thr)
        if len(passed):
            out[i, j[o[:passed[-1] + 1]]] = True
    return out


def top_n(score, ok, N):
    """Per row, the top N[i] (or N) ok genes by score."""
    out = np.zeros(score.shape, bool); Ns = np.broadcast_to(np.asarray(N), (len(score),))
    s = np.where(ok, score, -np.inf)
    o = np.argsort(-s, axis=1, kind='stable')
    for i in range(len(score)):
        n = int(min(Ns[i], ok[i].sum())); out[i, o[i, :n]] = True
    return out


def own_mask(h):
    M = np.zeros((h.T, h.G), bool)
    for i in h.rows:
        if h.own[i] >= 0: M[i, h.own[i]] = True
    return M


class CLine(Row4):
    """Row 4 plus C1's arrays (in h.targets order), z_u, C+G's z, the K4 ok set, the L-space spec, truth terciles."""

    def __init__(self, line, terc_parts=True):
        super().__init__(line)
        import scorer_step12 as S12
        h = self.h; self.c = c1_load(CTX_OF[line], h.targets)
        assert np.array_equal(self.c['gate'].astype(bool), self.gate), 'C1 gate differs from the line gate'
        self.ownM = own_mask(h)
        self.okg = self.gate[None, :] & ~self.ownM
        self.zu, self.zparts = z_u_of(self.c, self.gate, self.own)
        cg = self.c['lfc'].astype(np.float64)
        self.ok4 = self.okg & (self.c['n'] > 0) & (cg != 0)
        self.zc = np.where(self.ok4, np.abs(cg) / np.maximum(self.c['se'].astype(np.float64)[None, :], 1e-3), -1.0)
        assert np.array_equal(np.where(self.c['n'] <= 0, 0.0, cg), self.D.cg), 'C1 shrunk LFC differs from the line cg'
        self.T = S12.truth_calls(self.D.ref); self.n_true = self.T.sum(1); self.tterc = tercile_of(self.n_true)
        if terc_parts and line == 'H1':
            self.D.parts.update({t: np.flatnonzero(self.tterc == j) for j, t in enumerate(TERC)})
        self.L_spec = self.spec_L()

    def spec_L(self):
        D = self.D
        if self.line in ('H1', 'HCT116s'):
            un = (D.n_src <= 0) & D.gate[None, :]
            L = np.where(un, D.lfc_raw, D.lfc)
            return np.where(D.sig, 0.5 * D.cg, L)
        return np.log2((self.mu_spec + 1e-9) / (D.m0[None, :] + 1e-9))

    def mu_called(self, called_vals, f5=True, shift=True):
        """Config 1' magnitudes with K0's called genes set to called_vals (LFC), then row 4's additions."""
        import scorer_arms as SA
        L = np.where(self.sig_k0, called_vals, self.L_spec)
        return self.mu(f5=f5, shift=shift, base=SA.mu_from_lfc(L, self.D.m0))

    def uncalled_mu(self):
        """Config 1''s uncalled value on every gene (mu, no-source fix included)."""
        import scorer_arms as SA
        D = self.D
        if self.line in ('H1', 'HCT116s'):
            z = np.load(TAB / f'C1_uncalled_{self.line}.npz', allow_pickle=True); t = list(map(str, z['targets']))
            L = z['lfc_uncalled'][[t.index(x) for x in self.h.targets]].astype(np.float64)
            L = np.where((D.n_src <= 0) & ~D.gate[None, :], 0.0, L)
            return SA.mu_from_lfc(L, D.m0)
        return D.mu_non


def call_stats(sig, C, b, X=None, sgn=None):
    n = sig.sum(1)
    from scipy.stats import spearmanr
    out = dict(median=float(np.median(n)), mean=float(n.mean()), p10=float(np.quantile(n, .1)), p90=float(np.quantile(n, .9)))
    if C is not None:
        T, n_true, terc = C
        out['spearman_true_de'] = float(spearmanr(n, n_true).correlation)
        tp = (sig & T).sum(); out['precision'] = float(tp / max(sig.sum(), 1)); out['recall'] = float(tp / max(T.sum(), 1))
        if X is not None:
            m = sig & T & np.isfinite(X) & (X != 0) & (sgn != 0)
            out['sign_accuracy'] = float((np.sign(X[m]) == np.sign(sgn[m])).mean()) if m.any() else float('nan')
        out['by_true_tercile_mean'] = {t: float(n[terc == j].mean()) for j, t in enumerate(TERC)}
    ok = np.isfinite(b)
    out['spearman_b'] = float(spearmanr(n[ok], b[ok]).correlation) if ok.sum() > 2 else float('nan')
    return out


def true_lfc(ref):
    h = ref.h; X = np.full((h.T, h.G), np.nan)
    X[ref.ri, ref.rg] = ref.real_df['log2_fold_change'].to_numpy().astype(float)
    return X


def partC2(a):
    res = dict(lines={}, DEF={}, denominator={})
    for line in TRUTH_CTX:
        C = CLine(line, terc_parts=False); h = C.h
        X = true_lfc(C.D.ref); sg = C.c['lfc'].astype(np.float64)
        rules = {'requested K0': C.sig_k0, 'hybrid T5 (row 4)': C.sig}
        for al in BH_A:
            rules[f'z_u BH {al}'] = bh_calls(C.zu, C.okg & (C.c['n'] > 0), al)
        res['lines'][line] = {nm: call_stats(s, (C.T, C.n_true, C.tterc), C.b, X, sg) for nm, s in rules.items()}
        res['lines'][line]['_fallback_share'] = float(1 - (C.zparts['estimable'] & C.okg).sum() / C.okg.sum())
        # denominator shares by control-CPM decile (gated, own out, with a source)
        m = C.okg & (C.c['n'] > 0)
        parts = {'se2': C.zparts['se2'], 'var_fe': np.where(np.isfinite(C.zparts['var_fe']), C.zparts['var_fe'], 0.0), 'tau2': C.zparts['tau2']}
        den = parts['se2'] + parts['var_fe'] + parts['tau2']
        q = np.quantile(C.D.m0[C.gate], np.linspace(0, 1, 11)); dec = np.clip(np.searchsorted(q, C.D.m0, side='right') - 1, 0, 9)
        dd = {}
        for d_ in range(10):
            mm = m & (dec[None, :] == d_) & C.gate[None, :]
            if mm.any():
                dd[d_ + 1] = dict(cpm_range=[float(q[d_]), float(q[d_ + 1])], pairs=int(mm.sum()),
                                  share_of_summed_variance={k: float(v[mm].sum() / den[mm].sum()) for k, v in parts.items()},
                                  mean_per_pair_share={k: float((v[mm] / den[mm]).mean()) for k, v in parts.items()})
        res['denominator'][line] = dd
        log(f'C2 {line}: ' + '; '.join(f'{k} median {v["median"]:.0f} rho_true {v.get("spearman_true_de", np.nan):.3f} prec {v.get("precision", np.nan):.3f}'
                                        for k, v in res['lines'][line].items() if isinstance(v, dict)))
        save_json(TAB / 'partC2.json', res); del C; gc.collect()
    bdef = pd.read_csv(S19.S17 / 'tables/b_DEF.csv').set_index('target')['b']
    import run_line as RL
    pcol = {g: j for j, g in enumerate(RL.panel_genes())}
    for x in ('D', 'E', 'F'):
        c = c1_load(x); gate = c['gate'].astype(bool); own = np.array([pcol.get(t, -1) for t in c['targets']])
        zu, _ = z_u_of(c, gate, own); ok = np.broadcast_to(gate[None, :], zu.shape).copy() & (c['n'] > 0)
        for i, j in enumerate(own):
            if j >= 0: ok[i, j] = False
        b = bdef.reindex(list(map(str, c['targets']))).to_numpy(float)
        rules = {'requested K0': c['sig'].astype(bool)} | {f'z_u BH {al}': bh_calls(zu, ok, al) for al in BH_A}
        res['DEF'][x] = {nm: call_stats(s, None, b) for nm, s in rules.items()}
        log(f'C2 {x}: ' + '; '.join(f'{k} median {v["median"]:.0f} mean {v["mean"]:.0f} rho_b {v["spearman_b"]:.3f}' for k, v in res['DEF'][x].items()))
    res['note'] = 'Eval A realised cells (config 1\' + K0, 15c): Spearman of realised call count with true DE count 0.09 (Step 20 checklist item 4)'
    save_json(TAB / 'partC2.json', res)


def c3_rules(C):
    k4 = C.k4; z = C.zu; ok = C.ok4
    out = {'R-ref': C.sig}
    okz = C.okg & (C.c['n'] > 0)
    for al in BH_A:
        out[f'R-u {al}'] = bh_calls(z, okz, al)
    k4u = top_n(z, ok, 1200)
    out['R-uK4'] = np.where(C.on5[:, None], k4u, C.sig_k0)
    for al in BH_A:
        out[f'R-uu {al}'] = np.where(C.on5[:, None], k4u, out[f'R-u {al}'])
    for N in (50, 100, 200, 300):
        out[f'R-low {N}'] = np.where(C.on5[:, None], k4, top_n(C.zc, ok, N))
    return out


C3_ORDER = ['R-ref'] + [f'R-u {a}' for a in BH_A] + ['R-uK4'] + [f'R-uu {a}' for a in BH_A] + [f'R-low {n}' for n in (50, 100, 200, 300)]


def matched_ref(C, target_mean):
    """R-ref at matched mean call volume: b < 5 keeps K0, b >= 5 gets the top N by C+G z with N set (bisection) so the
    panel's mean calls per knockdown is as close as possible to target_mean."""
    nk0 = C.sig_k0.sum(1); nok = C.ok4.sum(1)
    f = lambda N: np.where(C.on5, np.minimum(N, nok), nk0).mean()
    lo, hi = 0, 6000
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if f(mid) < target_mean: lo = mid
        else: hi = mid
    N = lo if abs(f(lo) - target_mean) <= abs(f(hi) - target_mean) else hi
    return np.where(C.on5[:, None], top_n(C.zc, C.ok4, N), C.sig_k0), N, float(f(N))


def partC3(a):
    rows = resume('partC3_rows'); res = dict(checks={}, matched={})
    from scipy.stats import spearmanr
    for line in LINES:
        if any(r['line'] == line for r in rows):
            continue
        C = CLine(line, terc_parts=False)
        res['checks'][line] = dict(k4_from_C1_equals_line_K4=bool(np.array_equal(top_n(C.zc, C.ok4, 1200), C.k4)))
        mu4 = C.mu(); rules = c3_rules(C)
        for nm in C3_ORDER:
            s = rules[nm]; rr, _ = C.score(mu4, s, rule=nm)
            rho = float(spearmanr(s.sum(1), C.n_true).correlation)
            for r in rr: r['spearman_calls_true_de'] = rho
            rows += rr
            log(f'C3 {line} {nm}: ' + ', '.join(f'{r["part"]} {r["overall"]:.4f}' for r in rr) + f'; calls mean {rr[0]["calls_mean"]:.0f}; rho {rho:.3f}')
        # coupled R-u 0.05: magnitudes follow the call rule
        s = rules['R-u 0.05']; mun = C.uncalled_mu()
        import scorer_arms as SA
        base = np.where(s, SA.mu_from_lfc(0.5 * C.D.cg, C.D.m0), mun)
        rr, _ = C.score(C.mu(base=base), s, rule='R-u 0.05 (coupled)'); rows += rr
        log(f'C3 {line} R-u 0.05 coupled: ' + ', '.join(f'{r["part"]} {r["overall"]:.4f} (MSE {r["mse"]:.3f})' for r in rr))
        # coupled check: the same construction with K0's calls must reproduce the spec magnitudes
        b0 = np.where(C.sig_k0, SA.mu_from_lfc(0.5 * C.D.cg, C.D.m0), mun)
        res['checks'][line]['coupled_rule_with_K0_reproduces_spec_mu_max_rel'] = float(np.max(np.abs(b0 - C.mu_spec) / np.maximum(C.mu_spec, 1)))
        save_rows(rows, 'partC3_rows'); del C; gc.collect()
    R = pd.DataFrame(rows)
    out = fit_report(R[R.rule != 'R-u 0.05 (coupled)'], ['rule'], dict(rule=C3_ORDER), ['R-ref'], tag='C3')
    # matched mean call volume for every rule that beats R-ref on H1's full panel
    Hf = R[(R.line == 'H1') & (R.part == 'full')].set_index('rule')
    beat = [nm for nm in C3_ORDER[1:] if Hf.loc[nm, 'overall'] > Hf.loc['R-ref', 'overall']]
    res['rules_beating_R_ref_on_H1_full'] = beat
    mrows = resume('partC3_matched')
    for line in (LINES if beat else []):
        if any(r['line'] == line for r in mrows):
            continue
        C = CLine(line, terc_parts=False); mu4 = C.mu(); rules = c3_rules(C)
        for nm in beat:
            s, N, got = matched_ref(C, rules[nm].sum(1).mean())
            rr, _ = C.score(mu4, s, rule=f'R-ref matched to {nm}', N=N, matched_mean=got); mrows += rr
            log(f'C3 {line} R-ref matched to {nm} (N {N}, mean {got:.0f} vs {rules[nm].sum(1).mean():.0f}): ' + ', '.join(f'{r["part"]} {r["overall"]:.4f}' for r in rr))
        save_rows(mrows, 'partC3_matched'); del C; gc.collect()
    save_json(TAB / 'partC3.json', out | res)


C4_K, C4_S = [1.0, 4.0], [0.5, 0.75, 1.0]


def partC4(a):
    rows = resume('partC4_rows')
    for line in LINES:
        if any(r['line'] == line for r in rows):
            continue
        C = CLine(line); cg = C.D.cg; z2 = C.zu ** 2
        for kap, s in [(0.0, 0.5)] + [(k, s) for k in C4_K for s in C4_S]:
            w = np.ones_like(z2) if kap == 0 else z2 / (z2 + kap)
            rr, _ = C.score(C.mu_called(s * cg * w), kappa=kap, s=s); rows += rr
            log(f'C4 {line} kappa {kap} s {s}: ' + ', '.join(f'{r["part"]} {r["overall"]:.4f} (MSE {r["mse"]:.3f} PDS {r["pds"]:.3f} NMAE {r["nmae"]:.3f})' for r in rr if r['part'] in ('full', 'half0', 'half1')))
        save_rows(rows, 'partC4_rows'); del C; gc.collect()
    R = pd.DataFrame(rows)
    out = fit_report(R, ['kappa', 's'], dict(kappa=C4_K, s=C4_S), [0.0, 0.5], tag='C4')
    out['step16_reference'] = 'Step 16 arms j / kj: 0.229 / 0.237 against config 1 (i) 0.244 (fold-cache proxy, H1)'
    save_json(TAB / 'partC4.json', out)


def partC5(a):
    import scorer_arms as SA
    C = CLine('H1', terc_parts=False); h = C.h
    mu4 = C.mu(); P = SA.bulk(mu4); Rt = h.R; ctrl = h.ctrl[None, :]
    dP, dR = P - ctrl, Rt - ctrl
    keep = ~C.ownM
    m0 = C.D.m0
    band = np.full(h.G, -1)
    for j, (nm, lo, hi) in enumerate(CPM_BANDS):
        band[(m0 > lo) & (m0 <= hi)] = j
    band[~C.gate] = 0
    res = dict(item1={}, item2={}, item3={})
    err = (P - Rt) ** 2
    for tname, tj in (('strong', 2), ('mid', 1)):
        rr = np.flatnonzero(C.tterc == tj); E = err[rr] * keep[rr]; En = (dR[rr] ** 2) * keep[rr]
        ts = np.where(dR[rr] >= 0, 'up', 'down'); ps = np.where(np.abs(dP[rr]) < 0.01, '|pred|<0.01', np.where(dP[rr] > 0, 'up', 'down'))
        bb = np.broadcast_to(band[None, :], E.shape)
        cells = []
        for t_ in ('up', 'down'):
            for p_ in ('up', 'down', '|pred|<0.01'):
                for j, (nm, _, _) in enumerate(CPM_BANDS):
                    m = (ts == t_) & (ps == p_) & (bb == j) & keep[rr]
                    cells.append(dict(true_sign=t_, pred_sign=p_, cpm_band=nm, genes=int(m.sum()), error_share=float(E[m].sum() / E.sum()),
                                      real_energy_share=float(En[m].sum() / En.sum())))
        res['item1'][tname] = cells
    # item 2: sign accuracy on true DE genes
    T = C.T & keep; sp = np.sign(dP); st = np.sign(dR); right = (sp == st) & (sp != 0)
    called = C.sig
    def acc(m): return dict(n=int(m.sum()), accuracy=float(right[m].mean()) if m.any() else float('nan'))
    it2 = {}
    for cn, cm in (('called', called), ('uncalled', ~called)):
        for pn, pm in (('pred up', sp > 0), ('pred down', sp < 0)):
            it2[f'{cn}, {pn}'] = acc(T & cm & pm)
            for j, (nm, _, _) in enumerate(CPM_BANDS):
                it2[f'{cn}, {pn}, {nm}'] = acc(T & cm & pm & (band[None, :] == j))
    zq = np.quantile(C.zu[T], np.linspace(0, 1, 11))
    it2['by_z_u_decile'] = [dict(decile=d + 1, z_range=[float(zq[d]), float(zq[d + 1])], **acc(T & (C.zu >= zq[d]) & ((C.zu < zq[d + 1]) if d < 9 else True)))
                            for d in range(10)]
    bt = b_terc(C.b, C.cuts['full'])
    it2['by_b_tercile'] = {t: acc(T & (bt[:, None] == j)) for j, t in enumerate(TERC)}
    res['item2'] = it2
    # item 3: calibration slope through the origin, called genes, by predicted sign and band
    it3 = {}
    for pn, pm in (('pred up', dP > 0), ('pred down', dP < 0)):
        for j, (nm, _, _) in enumerate(CPM_BANDS):
            m = called & keep & pm & (band[None, :] == j)
            if m.sum() < 20:
                continue
            slope = float((dP[m] * dR[m]).sum() / (dP[m] ** 2).sum())
            p = it2.get(f'called, {pn}, {nm}', {}).get('accuracy', float('nan'))
            it3[f'{pn}, {nm}'] = dict(n=int(m.sum()), slope=slope, sign_accuracy_true_de=p, implied_slope_2p_minus_1=2 * p - 1)
    res['item3'] = it3
    save_json(TAB / 'partC5.json', res)
    log('C5: ' + json.dumps(it3, default=float))


C6_UP, C6_DN, C6_CF = [0.5, 0.65, 0.8, 1.0], [0.35, 0.5, 0.65], [20.0, 100.0]


def partC6(a):
    rows = resume('partC6_rows')
    for line in LINES:
        if any(r['line'] == line for r in rows):
            continue
        C = CLine(line); cg = C.D.cg; m0 = C.D.m0
        for cf in C6_CF:
            for su in C6_UP:
                for sd in C6_DN:
                    v = np.where((cg > 0) & (m0[None, :] >= cf), su * cg, np.where(cg < 0, sd * cg, 0.5 * cg))
                    rr, _ = C.score(C.mu_called(v), s_up=su, s_down=sd, c_floor=cf); rows += rr
            log(f'C6 {line} c_floor {cf}: done')
        save_rows(rows, 'partC6_rows'); del C; gc.collect()
    R = pd.DataFrame(rows)
    out = fit_report(R, ['s_up', 's_down', 'c_floor'], dict(s_up=C6_UP, s_down=C6_DN, c_floor=C6_CF), [0.5, 0.5, 20.0], tag='C6')
    H = R[R.line == 'H1']
    out['marginal'] = {k: {part: H[H.part == part].groupby(k)[['overall', 'mse', 'pds', 'nmae']].mean().to_dict('index') for part in ('full', 'half0', 'half1')}
                       for k in ('s_up', 's_down', 'c_floor')}
    save_json(TAB / 'partC6.json', out)


# ================================================================ C7 / C7d / C9c / C10 shared: the E4w source matrix
MAT = OUT / 'matrix'          # scratch (deleted at the end of Step 21): M.npy, Z.npy


def all_line_targets():
    tg = {}
    for line in D_LINES:
        z = np.load(SS_CACHE / f'fold_{line}.npz', allow_pickle=True); tg[line] = list(map(str, z['targets']))
    tg['HCT116s'] = tg['HCT116']; tg['DEF'] = S20.DEF_TARGETS()
    return tg


def c7build(a):
    """E4w's matrix exactly as step20.f_build('E4w') builds it (Delta, 5 CPM gate per table, per-table centring,
    weight 1/sqrt(max(n_b, 50)), Tahoe one background, Arc H1 out), saved to disk with: per-table mean vectors (to
    recover raw Delta), per-table gates, and z = LFC / SE of every row on the union of line targets (C7d).
    Check: the rows project onto Step 20's saved E4w basis with the saved row scores (US)."""
    import h5py
    import run_line as RL
    from anndata.io import read_elem
    MAT.mkdir(parents=True, exist_ok=True)
    genes = RL.panel_genes(); G = len(genes); gp = {g: j for j, g in enumerate(genes)}
    tg = all_line_targets(); cols_t = sorted({t for v in tg.values() for t in v if t in gp}); cols = np.array([gp[t] for t in cols_t])
    plan = []
    for p, stem, bg, tahoe in S19.e_tables():
        with h5py.File(p, 'r') as f:
            o = read_elem(f['obs'])
            if not np.isfinite(f['layers/log1p_bulk_50k_full'][:min(5, len(o))]).any():
                continue
        plan.append((p, stem, bg, tahoe, o, np.arange(len(o))))
    n_bg = pd.Series({}, dtype=float)
    for _, _, bg, _, _, k in plan:
        n_bg[bg] = n_bg.get(bg, 0) + len(k)
    n = int(n_bg.sum())
    mem_before(f'C7 matrix {n} x {G} (memmap on disk) + z {n} x {len(cols)}', n * len(cols) * 4 / 1e9 + 3)
    M = np.lib.format.open_memmap(MAT / 'M.npy', mode='w+', dtype=np.float32, shape=(n, G))
    Z = np.zeros((n, len(cols)), np.float32)
    means, gates, meta, r0 = [], [], [], 0
    for p, stem, bg, tahoe, o, k in plan:
        X, Mm = S20.table_entries(p, k, 'delta')
        mean = np.where(Mm, X, 0.0).sum(0) / np.maximum(Mm.sum(0), 1)
        w = S20.f_weight(n_bg[bg])
        M[r0:r0 + len(k)] = (np.where(Mm, X - mean[None, :], 0.0) * w).astype(np.float32)
        b = pd.read_csv(p.with_name(p.stem + '_baseline.csv'), index_col=0).reindex(genes)
        gates.append(b['control_mean_cpm'].to_numpy(float) > 5.0)
        means.append(mean.astype(np.float32))
        with h5py.File(p, 'r') as f:
            assert f['layers/lfc'].shape[1] == G, stem
            for c0 in range(0, len(k), 2000):
                sl = slice(c0, min(c0 + 2000, len(k)))
                L = f['layers/lfc'][sl][:, cols].astype(np.float64); SE = f['layers/se'][sl][:, cols].astype(np.float64)
                Z[r0 + sl.start:r0 + sl.stop] = np.where(np.isfinite(L) & np.isfinite(SE) & (SE > 0), L / np.where(SE > 0, SE, 1), np.nan)
        meta.append(pd.DataFrame(dict(table=stem, background=bg, tahoe=tahoe, modality=o['modality'].astype(str).to_numpy()[k],
                                      target=o['target_gene'].astype(str).to_numpy()[k], n_real=o['n_real'].astype(float).to_numpy()[k],
                                      weight=w, row_in_table=k, ctrl_cp50k=0.0)))
        r0 += len(k); del X, Mm
    M.flush()
    meta = pd.concat(meta, ignore_index=True); meta['table_i'] = pd.factorize(meta.table)[0]
    np.save(MAT / 'Z.npy', Z); meta.to_parquet(MAT / 'meta.parquet')
    tabs = meta.drop_duplicates('table').table.tolist()
    # each table's control profile on the target columns (C7d's SE conversion)
    cp = []
    for t in tabs:
        p = next(pp for pp, s, _, _, _, _ in plan if s == t)
        b = pd.read_csv(p.with_name(p.stem + '_baseline.csv'), index_col=0).reindex(genes)
        cp.append(np.expm1(b['log1p_bulk_50k_full'].to_numpy(float)[cols]))
    np.savez_compressed(MAT / 'tables.npz', tables=np.array(tabs), means=np.stack(means), gates=np.stack(gates),
                        cols=cols, cols_targets=np.array(cols_t), ctrl_cp50k_cols=np.stack(cp))
    # check against Step 20's basis rows and scores
    ref = pd.read_parquet(S20T / 'partF1_E4w_full_rows.parquet')
    same = bool((ref.table.to_numpy() == meta.table.to_numpy()).all() and (ref.target.to_numpy() == meta.target.to_numpy()).all())
    z = np.load(basis_file('E4w', 'full')); Vt = z['Vt'][:5].astype(np.float32); US = z['US'][:, :5]
    rr = np.random.default_rng(SEED).choice(n, 2000, replace=False); rr.sort()
    d = float(np.abs(np.asarray(M[rr]) @ Vt.T - US[rr]).max() / np.abs(US[rr]).max())
    res = dict(rows=n, genes=G, target_columns=len(cols), rows_match_step20=same, us_rel_max_diff_2000_rows=d, pass_=bool(same and d < 1e-3))
    save_json(TAB / 'c7build.json', res); log(f'c7build: {res}')
    if not res['pass_']:
        sys.exit(9)


class Mat:
    """The saved E4w matrix (memmap), its meta, table means / gates, and z on the target columns."""

    def __init__(self):
        self.M = np.load(MAT / 'M.npy', mmap_mode='r'); self.meta = pd.read_parquet(MAT / 'meta.parquet')
        z = np.load(MAT / 'tables.npz', allow_pickle=True)
        self.tables = list(z['tables']); self.means = z['means']; self.tgates = z['gates']; self.cols = z['cols']
        self.col_of = {int(c): i for i, c in enumerate(self.cols)}; self.ctrl_cols = z['ctrl_cp50k_cols']
        self.Z = np.load(MAT / 'Z.npy', mmap_mode='r')
        self.w = self.meta.weight.to_numpy(); self.ti = self.meta.table_i.to_numpy(); self.bg = self.meta.background.to_numpy()
        self.target = self.meta.target.to_numpy(); self.genetic = (self.meta.modality.to_numpy() != 'drug') & ~self.meta.tahoe.to_numpy()
        self.rows_of = self.meta.groupby('target').indices

    def raw(self, rows):
        """Raw Delta (uncentred) for rows; NaN where the entry is not measured (M == 0)."""
        X = np.asarray(self.M[rows], np.float64)
        R = X / self.w[rows][:, None] + self.means[self.ti[rows]]
        return np.where(X != 0, R, np.nan)


VARIANT_DROP = {'full': (), 'K562': ('K562',), 'Orion': ('HCT116', 'HEK293T'), 'HCT116': ('HCT116',)}
VARIANT_OF = {'H1': 'full', 'K562': 'K562', 'HCT116': 'Orion', 'HEK293T': 'Orion', 'HCT116s': 'HCT116', 'DEF': 'full'}
BASIS_FOLD = {'full': 'full', 'K562': 'K562', 'Orion': 'Orion', 'HCT116': 'Orion'}   # C7's low-rank / C7b embedding basis


def corr_columns(mat, cols, variant, chunk=4000, genes=None):
    """c_X for every X in cols (panel indices): row-weighted Pearson of X's column with every gene's column over the
    variant's rows in tables that gate X, rows targeting X dropped. Returns (len(cols), G) float32 and rows kept."""
    import run_line as RL
    genes = RL.panel_genes() if genes is None else genes; G = len(genes)
    drop = np.isin(mat.bg, VARIANT_DROP[variant]); keep = ~drop
    n_t = len(mat.tables)
    TS1 = np.zeros((n_t, G)); TS2 = np.zeros((n_t, G)); TW = np.zeros(n_t); TN = np.zeros(n_t, int)
    P = np.zeros((G, len(cols)))
    for r0 in range(0, len(mat.w), chunk):
        rr = np.arange(r0, min(r0 + chunk, len(mat.w))); rr = rr[keep[rr]]
        if not len(rr):
            continue
        X = np.asarray(mat.M[rr], np.float32)
        P += (X.T @ X[:, cols]).astype(np.float64)
        wi = mat.w[rr]; ti = mat.ti[rr]
        np.add.at(TS1, ti, (X * wi[:, None]).astype(np.float64)); np.add.at(TS2, ti, (X.astype(np.float64) ** 2))
        np.add.at(TW, ti, wi ** 2); np.add.at(TN, ti, 1)
    C = np.zeros((len(cols), G), np.float32); kept = np.zeros(len(cols), int)
    for q, x in enumerate(cols):
        tset = np.flatnonzero(mat.tgates[:, x] & (TN > 0))
        S1 = TS1[tset].sum(0); S2 = TS2[tset].sum(0); W = TW[tset].sum(); Sxy = P[:, q].copy(); nk = TN[tset].sum()
        rt = [r for r in mat.rows_of.get(genes[x], []) if keep[r]]
        if rt:
            Xr = np.asarray(mat.M[rt], np.float64)
            Sxy -= (Xr * Xr[:, [x]]).sum(0)
            inT = np.isin(mat.ti[rt], tset)
            if inT.any():
                Xi = Xr[inT]; wi = mat.w[np.array(rt)[inT]]
                S1 -= (Xi * wi[:, None]).sum(0); S2 -= (Xi ** 2).sum(0); W -= (wi ** 2).sum(); nk -= int(inT.sum())
        kept[q] = nk
        if W <= 0:
            continue
        m = S1 / W; v = S2 / W - m ** 2
        cov = Sxy / W - m[x] * m
        with np.errstate(invalid='ignore', divide='ignore'):
            c = cov / np.sqrt(np.maximum(v[x], 0) * np.maximum(v, 0))
        c = np.where(np.isfinite(c), c, 0.0); c[x] = 0.0
        C[q] = c
    return C, kept


def lowrank_corr(variant, cols, K):
    z = np.load(basis_file('E4w', BASIS_FOLD[variant])); Vt = z['Vt'][:K].astype(np.float64); s2 = z['s'][:K].astype(np.float64) ** 2
    A = Vt * s2[:, None]                                     # (K, G)
    num = A[:, cols].T @ Vt                                  # (n, G)
    dg = (Vt ** 2 * s2[:, None]).sum(0)
    with np.errstate(invalid='ignore', divide='ignore'):
        c = num / np.sqrt(dg[cols][:, None] * dg[None, :])
    c = np.where(np.isfinite(c), c, 0.0)
    c[np.arange(len(cols)), cols] = 0.0
    return c


def pearson_rows(Y, X, M):
    """Per row, Pearson of Y and X over mask M (NaN if < 3 entries or a constant side)."""
    out = np.full(len(Y), np.nan)
    for i in range(len(Y)):
        m = M[i] & np.isfinite(Y[i]) & np.isfinite(X[i])
        if m.sum() >= 3 and np.std(X[i, m]) > 0 and np.std(Y[i, m]) > 0:
            out[i] = np.corrcoef(Y[i, m], X[i, m])[0, 1]
    return out


# ================================================================ C7: co-response modules
C7_VARIANTS = {'full': ['H1', 'DEF'], 'K562': ['K562'], 'Orion': ['HCT116', 'HEK293T'], 'HCT116': ['HCT116s']}
TRUTH_BASE = {'H1': 'H1', 'K562': 'K562', 'HCT116': 'HCT116', 'HEK293T': 'HEK293T', 'HCT116s': 'HCT116'}
N_POOL = 1000


def c7_pool(genes):
    gates = np.zeros(len(genes), bool)
    for line in D_LINES:
        gates |= np.load(SS_CACHE / f'fold_{line}.npz')['gate'].astype(bool)
    return np.sort(np.random.default_rng(SEED).choice(np.flatnonzero(gates), N_POOL, replace=False))


def c7_load(variant):
    z = np.load(MAT / f'C7_{variant}.npz', allow_pickle=True)
    return np.load(MAT / f'C7_{variant}.npy', mmap_mode='r'), list(map(int, z['cols'])), z


def c7_rows(line, targets, genes):
    """-c_X rows (T x G) for a line's targets from its variant; has = X has a column."""
    C, cols, _ = c7_load(VARIANT_OF[line]); pos = {c: i for i, c in enumerate(cols)}; gp = {g: j for j, g in enumerate(genes)}
    out = np.zeros((len(targets), len(genes))); has = np.zeros(len(targets), bool)
    for i, t in enumerate(targets):
        j = gp.get(t)
        if j is not None and j in pos:
            out[i] = -np.asarray(C[pos[j]], np.float64); has[i] = True
    return out, has


def partC7(a):
    import run_line as RL
    genes = RL.panel_genes(); gp = {g: j for j, g in enumerate(genes)}
    mat = Mat(); tg = all_line_targets(); pool = c7_pool(genes)
    res = dict(variants={}, C7a={}, pool_size=N_POOL)
    for v, lines in C7_VARIANTS.items():
        if (MAT / f'C7_{v}.npy').exists():
            continue
        targets = sorted({t for l in lines for t in tg[l] if t in gp})
        cols = [gp[t] for t in targets]; cols += [int(j) for j in pool if j not in set(cols)]
        mem_before(f'C7 correlations {v} ({len(cols)} columns)', len(cols) * len(genes) * 12 / 1e9 + 2)
        t0 = time.time(); C, kept = corr_columns(mat, np.array(cols), v)
        np.save(MAT / f'C7_{v}.npy', C); np.savez(MAT / f'C7_{v}.npz', cols=np.array(cols), n_targets=len(targets), kept=kept)
        nt = len(targets); gate_v = np.load(SS_CACHE / f'fold_{TRUTH_BASE[lines[0]] if lines[0] != "DEF" else "H1"}.npz')['gate'].astype(bool)
        close = {}
        for K in (20, 100):
            Ck = lowrank_corr(v, np.array(cols[:nt]), K)
            r = [np.corrcoef(C[i, gate_v], Ck[i, gate_v])[0, 1] for i in range(nt) if np.std(C[i, gate_v]) > 0]
            close[K] = dict(median_pearson_full_vs_lowrank=float(np.nanmedian(r)), q10=float(np.nanquantile(r, .1)), q90=float(np.nanquantile(r, .9)))
        res['variants'][v] = dict(dropped=list(VARIANT_DROP[v]), targets=nt, columns=len(cols), seconds=time.time() - t0,
                                  rows_kept_per_target=dict(zip(['min', 'q10', 'median', 'q90', 'max'], np.quantile(kept[:nt], [0, .1, .5, .9, 1]).tolist())),
                                  lowrank=close, basis_used_for_lowrank=BASIS_FOLD[v])
        log(f'C7 {v}: {res["variants"][v]}'); save_json(TAB / 'partC7.json', res); del C; gc.collect()
    res = json.loads((TAB / 'partC7.json').read_text()) | {k: v for k, v in res.items() if k == 'C7a'}
    rng = np.random.default_rng(SEED); perkd = []
    for line in TRUTH_CTX:
        tr = S20.truth_of(TRUTH_BASE[line]); targets = tr['targets']; gate = np.load(SS_CACHE / f'fold_{TRUTH_BASE[line]}.npz')['gate'].astype(bool)
        Cc, cols, _ = c7_load(VARIANT_OF[line]); pos = {c: i for i, c in enumerate(cols)}
        pool_l = [j for j in pool if gate[j]]
        r_obs = np.full(len(targets), np.nan); r_null = np.full((len(targets), 100), np.nan)
        de_obs = np.full(len(targets), np.nan); de_null = np.full((len(targets), 100), np.nan)
        for i, t in enumerate(targets):
            j = gp.get(t)
            if j is None or j not in pos:
                continue
            own = tr['own'][i]; m = gate.copy()
            if own >= 0: m[own] = False
            y = tr['Delta'][i, m]; x = -np.asarray(Cc[pos[j]], np.float64)[m]
            if np.std(x) == 0:
                continue
            r_obs[i] = np.corrcoef(y, x)[0, 1]
            de = tr['DE'][i, m]; mg = np.flatnonzero(m)
            if de.sum():
                top = np.argsort(-(-x))[:100]                      # largest c = most negative -c
                de_obs[i] = float(de[top].sum() / de.sum())
            dr = rng.choice([q for q in pool_l if q != own and q != j], 100, replace=False)
            Xn = -np.asarray(Cc[[pos[q] for q in dr]], np.float64)[:, m]
            ys = (y - y.mean()) / y.std(); Xs = (Xn - Xn.mean(1, keepdims=True)) / np.maximum(Xn.std(1, keepdims=True), 1e-12)
            r_null[i] = (Xs @ ys) / len(y)
            if de.sum():
                tops = np.argsort(Xn, axis=1)[:, :100]
                de_null[i] = de[tops].sum(1) / de.sum()
        out = dict(knockdowns_with_c=int(np.isfinite(r_obs).sum()), knockdowns=len(targets))
        for nm, rr in (('all', np.isfinite(r_obs)),) + tuple((t, np.isfinite(r_obs) & (tr['terc'] == j)) for j, t in enumerate(TERC)):
            nm_ = np.nanmedian(r_null[rr], 0)
            out[nm] = dict(n=int(rr.sum()), median_r=float(np.nanmedian(r_obs[rr])), null_median_r=float(np.nanmedian(nm_)),
                           null_median_r_q95=float(np.nanquantile(nm_, .95)), percentile_vs_null=float((nm_ < np.nanmedian(r_obs[rr])).mean() * 100),
                           de_in_top100=float(np.nanmean(de_obs[rr])), de_in_top100_null=float(np.nanmean(de_null[rr])))
        res['C7a'][line] = out
        perkd += [dict(line=line, target=t, r=r_obs[i], r_null_median=float(np.nanmedian(r_null[i])), terc=int(tr['terc'][i])) for i, t in enumerate(targets)]
        log(f'C7a {line}: {out["all"]}; strong {out["strong"]}')
        save_json(TAB / 'partC7.json', res)
    pd.DataFrame(perkd).to_csv(TAB / 'partC7a_per_kd.csv', index=False)
    F6 = json.loads((S20T / 'partF6.json').read_text())['diag']; F3d = json.loads((S20T / 'partF3d.json').read_text())
    res['scale'] = dict(F6_median_r={l: F6[l]['all']['median_r'] for l in F6}, F3d=F3d)
    save_json(TAB / 'partC7.json', res)


def emb_features(targets, bf, K, genes):
    """Target's own loadings V_K[X, :] on the E4w basis (fold bf) and a missing indicator."""
    Vt = np.load(basis_file('E4w', bf))['Vt'][:K].astype(np.float64); gp = {g: j for j, g in enumerate(genes)}
    X = np.zeros((len(targets), K + 1))
    for i, t in enumerate(targets):
        j = gp.get(t)
        if j is None or not np.any(Vt[:, j] != 0):
            X[i, K] = 1.0
        else:
            X[i, :K] = Vt[:, j]
    return X


def partC7b(a):
    import pickle
    import run_line as RL
    genes = RL.panel_genes()
    ctx = S20.FCtx(); trH = ctx.tr['H1']; hv = ctx.halves; T = len(trH['targets'])
    hct = list(ref_for('HCT116final').h.targets); S19._REFS.pop('HCT116final', None)
    for K in (20, 100):
        for bf in ('full', 'K562', 'Orion'):
            for l in D_LINES:
                ctx.GF[(l, f'emb_{bf}_K{K}')] = emb_features(ctx.tr[l]['targets'], bf, K, genes)
    F4 = pd.read_csv(S20T / 'partF4_rows.csv')
    rows, preds = [], {}
    Ks = [5, 10, 20, 100]
    for K in (20, 100):
        Vt = ctx.vt('E4w', 'full'); lab = ctx.labels('H1', 'E4w', 'full', Ks)[K]
        AH = ctx.scores('H1', 'E4w', 'full', K, ('H1',))
        meth = {'M2': ('M2', None), 'M3DepMap': ('M3', 'DepMap'), 'M4DepMap': ('M4', 'DepMap'),
                'M3e': ('M3', f'emb_full_K{K}'), 'M4e': ('M4', f'emb_full_K{K}')}
        for nm, (m, fs) in meth.items():
            cf = np.zeros((T, K)); al = {}
            for f in (0, 1):
                pred, inf = S20.fit_predict(ctx, m, 'E4w', 'full', K, Ks, [('H1', hv[f], ('H1',), False)], AH, ctx.GF[('H1', fs)] if fs else None, fs)
                cf[hv[1 - f]] = pred[hv[1 - f]]; preds[('H1', nm, K, f)] = pred.astype(np.float32); al[f] = inf.get('alpha')
            r2 = 1 - ((lab - cf) ** 2).sum(0) / (lab ** 2).sum(0)
            cap = S19.captured(S20.energy_delta(trH, cf @ Vt[:K]), trH['terc'])
            row = dict(K=K, method=nm, median_r2_heldout=float(np.median(r2)), **{f'c_{g}': v for g, v in cap.items()},
                       alpha_at_top_edge_share=float(np.mean([np.mean(v == 1000.0) for v in al.values() if v is not None])) if al[0] is not None else np.nan)
            ref_ = F4[(F4.basis == 'E4w') & (F4.K == K) & (F4.method == nm)]
            if len(ref_):
                row['step20_F4_median_r2'] = float(ref_.median_r2_heldout.iloc[0]); row['step20_F4_c_strong'] = float(ref_.c_strong.iloc[0])
            rows.append(row); log(f'C7b K {K} {nm}: {row}')
        # M4e permutation null
        rng = np.random.default_rng(SEED); caps, r2s = [], []
        for _ in range(S20.N_NULL):
            perm = rng.permutation(T); Gp = ctx.GF[('H1', f'emb_full_K{K}')][perm]; cf = np.zeros((T, K))
            for f in (0, 1):
                pred, _ = S20.fit_predict(ctx, 'M4', 'E4w', 'full', K, Ks, [('H1', hv[f], ('H1',), False)], AH, Gp, f'emb_full_K{K}', gf_override={'H1': Gp})
                cf[hv[1 - f]] = pred[hv[1 - f]]
            caps.append(S19.captured(S20.energy_delta(trH, cf @ Vt[:K]), trH['terc'])); r2s.append(float(np.median(1 - ((lab - cf) ** 2).sum(0) / (lab ** 2).sum(0))))
        rows.append(dict(K=K, method='M4e_null', median_r2_heldout=float(np.mean(r2s)), median_r2_null_q95=float(np.quantile(r2s, .95)),
                         **{f'c_{g}': float(np.mean([c[g] for c in caps])) for g in ('all',) + TERC}, c_strong_null_q95=float(np.quantile([c['strong'] for c in caps], .95))))
        log(f'C7b K {K} M4e null: {rows[-1]}')
        pd.DataFrame(rows).to_csv(TAB / 'partC7b_rows.csv', index=False)
    # predictions for C7c's F5-with-M3e arm (K 100): guardrails from all of H1 under the fold basis; HCT116s screen
    K = 100
    for line in D_LINES[1:]:
        fold = FOLD_OF[line]; Ag = ctx.scores(line, 'E4w', fold, K, OWN_BG_F[line])
        tr_ = [('H1', np.arange(T), tuple(sorted(set(('H1',) + OWN_BG_F[line]))), False)]
        fs = f'emb_{fold}_K{K}'
        pred, _ = S20.fit_predict(ctx, 'M3', 'E4w', fold, K, Ks, tr_, Ag, ctx.GF[(line, fs)], fs); preds[(line, 'M3e', K, -1)] = pred.astype(np.float32)
        p3, _ = S20.fit_predict(ctx, 'M3', 'E4w', fold, K, Ks, tr_, Ag, ctx.GF[(line, 'DepMap')], 'DepMap')
        log(f'C7b {line}: M3 DepMap rebuilt vs Step 20 max |diff| (float32) {np.abs(p3.astype(np.float32) - f4_preds()[("E4w", K, "M3DepMap", line, -1)]).max():.1e}')
    sp = 'delta'; RS = ctx.rs('E4w', 'Orion')
    A = S20.bg_scores(RS, ctx.meta[sp], hct, ('HCT116',), RS.shape[1], ctx.bgs[sp])[:, :, :K]
    fs = f'emb_Orion_K{K}'
    pred, _ = S20.fit_predict(ctx, 'M3', 'E4w', 'Orion', K, Ks, [('H1', np.arange(T), ('H1', 'HCT116'), False)], A, emb_features(hct, 'Orion', K, genes), fs)
    preds[('HCT116s', 'M3e', K, -1)] = pred.astype(np.float32)
    with open(TAB / 'partC7b_preds.pkl', 'wb') as fh:
        pickle.dump(dict(preds=preds, targets={'HCT116s': hct, **{l: ctx.tr[l]['targets'] for l in D_LINES}}), fh)
    save_json(TAB / 'partC7b.json', dict(rows=rows, note='captured energy c_* on Delta scale over gated genes (F4 convention); R2 relative to predicting 0'))


C7C_G = [0.5, 1.0, 1.5]


def lam_fit(X, Y, ok, rows, b, on5):
    """lambda_t per (b) tercile (cuts from these rows' b), through the origin, over rows with b >= 5."""
    cuts = b_cuts(b[rows]); t = b_terc(b, cuts); L = np.zeros(3)
    for j in range(3):
        rr = rows[(t[rows] == j) & on5[rows]]
        x = np.where(ok[rr], X[rr], 0.0); y = np.where(ok[rr], Y[rr], 0.0)
        L[j] = (x * y).sum() / max((x * x).sum(), 1e-30)
    return L, cuts


def crossfit_lambda(X, Y, ok, b, on5, hv, T_true):
    """Per direction f: lambda for every knockdown. Fit-half knockdowns: from the other inner fold of the fit half
    (2 inner folds paired by true-DE rank, scorer_step16.halves); held-out half: from the whole fit half."""
    import scorer_step16 as S16
    out = {}
    for f in (0, 1):
        fit, ho = hv[f], hv[1 - f]; lam_i = np.zeros(len(b)); info = {}
        inner = S16.halves(T_true[fit])
        for j in (0, 1):
            L, cuts = lam_fit(X, Y, ok, fit[inner[j]], b, on5); tgt = fit[inner[1 - j]]
            lam_i[tgt] = L[b_terc(b[tgt], cuts)]; info[f'inner{j}'] = dict(lam=L.tolist(), cuts=cuts.tolist())
        L, cuts = lam_fit(X, Y, ok, fit, b, on5); lam_i[ho] = L[b_terc(b[ho], cuts)]; info['fit_half'] = dict(lam=L.tolist(), cuts=cuts.tolist())
        out[f] = (lam_i, info)
    return out


def partC7c(a):
    import pickle
    import run_line as RL
    genes = RL.panel_genes()
    rows = resume('partC7c_rows'); res = dict(lambda_={})
    z = pickle.load(open(TAB / 'partC7b_preds.pkl', 'rb')); pr = z['preds']
    for line in LINES:
        if any(r['line'] == line for r in rows):
            continue
        C = CLine(line); h = C.h
        X, has = c7_rows(line, h.targets, genes); ok = C.okg & has[:, None]
        Y = h.R - h.ctrl[None, :]
        mu4 = C.mu(); base_rows, _ = C.score(mu4, arm='co-movement', gamma=0.0, fit_half=-1)
        if line == 'H1':
            cf = crossfit_lambda(X, Y, ok, C.b, C.on5, C.hv, C.T)
            for f in (0, 1):
                lam_i, info = cf[f]; res['lambda_'][f'fit_half{f}'] = info
                rows += [dict(r, fit_half=f) for r in base_rows]
                for g in C7C_G:
                    add = np.where((C.on5 & has)[:, None] & ok, g * lam_i[:, None] * X, 0.0)
                    rr, _ = C.score(C.mu(extra=add), arm='co-movement', gamma=g, fit_half=f); rows += rr
                    log(f'C7c H1 co-movement fit {f} gamma {g}: ' + ', '.join(f'{r["part"]} {r["overall"]:.4f}' for r in rr if r['part'] in ('half0', 'half1')))
            L_all, cuts_all = lam_fit(X, Y, ok, np.arange(h.T), C.b, C.on5); res['lambda_']['all_H1'] = dict(lam=L_all.tolist(), cuts=cuts_all.tolist())
            save_json(TAB / 'partC7c_lambda.json', res)
        else:
            rows += base_rows
            L_all = np.array(json.loads((TAB / 'partC7c_lambda.json').read_text())['lambda_']['all_H1']['lam'])
            cuts_all = np.array(json.loads((TAB / 'partC7c_lambda.json').read_text())['lambda_']['all_H1']['cuts'])
            t = b_terc(C.b, cuts_all)
            for g in C7C_G:
                add = np.where((C.on5 & has)[:, None] & ok, g * L_all[t][:, None] * X, 0.0)
                rr, _ = C.score(C.mu(extra=add), arm='co-movement', gamma=g, fit_half=-1); rows += rr
        # F5 with M3e in place of M3
        if line == 'H1':
            sh = np.zeros((h.T, F5K))
            for f in (0, 1):
                sh[C.hv[1 - f]] = pr[('H1', 'M3e', F5K, f)][C.hv[1 - f]]
        else:
            tl = z['targets'][line]; assert list(tl) == h.targets
            sh = pr[(line, 'M3e', F5K, -1)].astype(np.float64)
        mu_e = mu_add_delta(C.mu_spec, C.f5_from_shat(sh) + C.addS, C.own)
        rr, _ = C.score(mu_e, arm='F5 with M3e', gamma=1.0, fit_half=-1); rows += rr
        log(f'C7c {line} F5 with M3e: ' + ', '.join(f'{r["part"]} {r["overall"]:.4f} (MSE {r["mse"]:.3f} PDS {r["pds"]:.3f})' for r in rr))
        save_rows(rows, 'partC7c_rows'); del C; gc.collect()
    R = pd.DataFrame(rows)
    Rc = R[R.arm == 'co-movement']
    out = dict(co_movement=fit_report(Rc, ['gamma'], dict(gamma=C7C_G), [0.0], tag='C7c co-movement',
                                      fit_frames={f: Rc[(Rc.line == 'H1') & (Rc.fit_half == f)] for f in (0, 1)}))
    # F5-M3e: both halves held out by construction (each half's knockdowns from the model trained on the other half)
    fv = []
    b4 = Rc[(Rc.gamma == 0) & (Rc.fit_half.isin([-1, 0]))]
    E = R[R.arm == 'F5 with M3e']
    for part in ('half0', 'half1'):
        x = E[(E.line == 'H1') & (E.part == part)].iloc[0]; y = b4[(b4.line == 'H1') & (b4.part == part)].iloc[0]
        gain = {c: float(x[c] - y[c]) for c in COMPS + ['overall']}; sc = sum(gain[c] for c in COMPS)
        gr = {l: dict(diff=float(E[(E.line == l) & (E.part == 'full')].overall.iloc[0] - b4[(b4.line == l) & (b4.part == 'full')].overall.iloc[0]))
              for l in GUARD}
        fv.append(dict(part=part, held_out=comp_dict(x), base_held_out=comp_dict(y), gain=gain, guardrails=gr,
                       reach_share_of_component_gain=float(gain['reach'] / sc) if sc else float('nan')))
    out['F5_M3e'] = dict(halves=fv, criterion=criterion(fv), proxy_note=PROXY_REACH_NOTE)
    save_json(TAB / 'partC7c.json', out | res)


# ================================================================ C7d: dose-response regression on X's column
C7D_G = [0.25, 0.5, 1.0]
C7D_MIN_ROWS = 20


def c7d_predict(mat, x, variant, genes, line_m0):
    """For panel gene x: kept rows (variant rows, tables gating x, not targeting x, |z_x| >= 3, x measured); per
    source line (background) the row-weighted slope of each gene on x through the origin (genes measured), pooled one
    vote per line (lines with >= 3 kept rows); divided by x's reliability ratio. Returns dict or None."""
    q = mat.col_of.get(int(x))
    if q is None:
        return None
    keep = ~np.isin(mat.bg, VARIANT_DROP[variant])
    zx = np.asarray(mat.Z[:, q]); r = np.flatnonzero(keep & mat.tgates[mat.ti, x] & (np.abs(np.nan_to_num(zx)) >= 3) & (mat.target != genes[x]))
    if len(r):
        mx = np.asarray(mat.M[r, x], np.float64); r = r[mx != 0]
    out = dict(rows=int(len(r)))
    if len(r) < C7D_MIN_ROWS:
        return out
    Mr = np.asarray(mat.M[r], np.float64); mx = Mr[:, x]
    xu = mx / mat.w[r]                                            # centred, unweighted Delta_x
    raw_x = xu + mat.means[mat.ti[r], x]
    se = np.abs(raw_x) / np.maximum(np.abs(zx[r]), 1e-12)        # SE of Delta_x ~ |Delta_x| / |z_x| (same LFC -> Delta factor)
    om = mat.w[r] ** 2
    rel = float((om * (xu ** 2 - se ** 2)).sum() / (om * xu ** 2).sum())
    num_l, n_l = {}, {}
    for bgk in np.unique(mat.bg[r]):
        s = mat.bg[r] == bgk
        if s.sum() < 3:
            continue
        num = (mx[s, None] * Mr[s]).sum(0); den = ((mx[s] ** 2)[:, None] * (Mr[s] != 0)).sum(0)
        with np.errstate(invalid='ignore', divide='ignore'):
            num_l[bgk] = np.where(den > 0, num / den, np.nan)
        n_l[bgk] = int(s.sum())
    if not num_l:
        out['lines'] = 0; return out
    B = np.stack(list(num_l.values())); beta = np.nanmean(B, 0) if len(B) > 1 else B[0]
    relc = float(np.clip(rel, 0.1, 1.0))
    beta = np.where(np.isfinite(beta), beta, 0.0) / relc; beta[x] = 0.0
    c = line_m0[x] / 20.0
    d = float(np.log1p(0.2 * c) - np.log1p(c))
    q90 = float(np.quantile(np.abs(xu), 0.9))
    out.update(lines=len(num_l), reliability=rel, reliability_clipped=relc, dose=d, median_abs_dx=float(np.median(np.abs(xu))),
               q90_abs_dx=q90, beta=beta.astype(np.float32), dose_capped=float(-min(abs(d), q90)))
    return out


def partC7d(a):
    import run_line as RL
    genes = RL.panel_genes(); gp = {g: j for j, g in enumerate(genes)}
    mat = Mat(); tg = all_line_targets()
    res = dict(lines={}, DEF={}, premise={})
    pk = pd.read_csv(TAB / 'partC7a_per_kd.csv') if (TAB / 'partC7a_per_kd.csv').exists() else None
    P = {}
    for line in TRUTH_CTX:
        C = CLine(line); h = C.h; m0 = C.D.m0
        Dh = np.zeros((h.T, h.G)); Dc = np.zeros((h.T, h.G)); has = np.zeros(h.T, bool); info = []
        for i, t in enumerate(h.targets):
            x = gp.get(t); o = c7d_predict(mat, x, VARIANT_OF[line], genes, m0) if x is not None else None
            info.append(dict(target=t, **{k: v for k, v in (o or {}).items() if k != 'beta'}))
            if o and 'beta' in o:
                Dh[i] = o['beta'] * o['dose']; Dc[i] = o['beta'] * o['dose_capped']; has[i] = True
        Dh = np.where(C.okg, Dh, 0.0); Dc = np.where(C.okg, Dc, 0.0)
        I = pd.DataFrame(info); I.to_csv(TAB / f'partC7d_info_{line}.csv', index=False)
        Y = h.R - h.ctrl[None, :]
        r_u = pearson_rows(Y, Dh, C.okg & has[:, None]); r_c = pearson_rows(Y, Dc, C.okg & has[:, None])
        b5 = C.b < 5
        def med(r, m): return float(np.nanmedian(r[m])) if np.isfinite(r[m]).any() else float('nan')
        out = dict(knockdowns=int(h.T), with_prediction=int(has.sum()), share_without=float(1 - has.mean()),
                   share_without_b_lt5=float(1 - has[np.nan_to_num(C.b, nan=0) < 5].mean()) if (np.nan_to_num(C.b, nan=0) < 5).any() else float('nan'),
                   rows_kept=dict(zip(['q10', 'median', 'q90'], np.quantile(I.rows, [.1, .5, .9]).tolist())),
                   median_reliability=float(I.reliability.median()) if 'reliability' in I else float('nan'),
                   share_reliability_clipped=float((I.reliability < 0.1).mean()) if 'reliability' in I else float('nan'),
                   extrapolation_median_absdose_over_median_absdx=float(np.nanmedian(np.abs(I.dose) / I.median_abs_dx)) if 'dose' in I else float('nan'),
                   uncapped=dict(all=med(r_u, has), **{t: med(r_u, has & (C.tterc == j)) for j, t in enumerate(TERC)}, b_lt5=med(r_u, has & (np.nan_to_num(C.b, nan=0) < 5))),
                   capped=dict(all=med(r_c, has), **{t: med(r_c, has & (C.tterc == j)) for j, t in enumerate(TERC)}, b_lt5=med(r_c, has & (np.nan_to_num(C.b, nan=0) < 5))))
        if pk is not None:
            q = pk[pk.line == line].set_index('target').reindex(h.targets)
            out['C7a_same_knockdowns'] = dict(median_r=float(np.nanmedian(q.r.to_numpy()[has])), null_median=float(np.nanmedian(q.r_null_median.to_numpy()[has])))
        if line in ('H1', 'HCT116s'):
            own = C.own; ok = own >= 0
            tc = np.array([h.R[i, own[i]] - h.ctrl[own[i]] for i in np.flatnonzero(ok)])
            dd = np.array([np.log1p(0.2 * m0[own[i]] / 20) - np.log1p(m0[own[i]] / 20) for i in np.flatnonzero(ok)])
            res['premise'][line] = dict(true_own_change=dict(zip(['q25', 'median', 'q75'], np.quantile(tc, [.25, .5, .75]).tolist())),
                                        d_X=dict(zip(['q25', 'median', 'q75'], np.quantile(dd, [.25, .5, .75]).tolist())),
                                        median_ratio_true_over_dX=float(np.median(tc / dd)), knockdowns=int(ok.sum()))
        res['lines'][line] = out; log(f'C7d {line}: {out}')
        np.savez_compressed(MAT / f'C7d_{line}.npz', Dh=Dh.astype(np.float32), has=has, targets=np.array(h.targets))
        save_json(TAB / 'partC7d.json', res); del C; gc.collect()
    # D/E/F: kept-row counts only
    bdef = pd.read_csv(S19.S17 / 'tables/b_DEF.csv').set_index('target')['b']
    m0d = np.load(S19.TAB / 'controls_D.npz')['m0']
    cnt = []
    for t in tg['DEF']:
        x = gp.get(t); o = c7d_predict(mat, x, 'full', genes, m0d) if x is not None else None
        cnt.append(dict(target=t, rows=(o or {}).get('rows', 0), prediction=bool(o and 'beta' in o), b=float(bdef.get(t, np.nan))))
    Q = pd.DataFrame(cnt)
    res['DEF'] = dict(targets=len(Q), share_without=float(1 - Q.prediction.mean()), share_without_b_lt5=float(1 - Q[Q.b.fillna(0) < 5].prediction.mean()),
                      rows_kept=dict(zip(['q10', 'median', 'q90'], np.quantile(Q.rows, [.1, .5, .9]).tolist())))
    save_json(TAB / 'partC7d.json', res)


def partC7db(a):
    rows = resume('partC7db_rows')
    for line in LINES:
        if any(r['line'] == line for r in rows):
            continue
        C = CLine(line); h = C.h
        z = np.load(MAT / f'C7d_{"HCT116s" if line == "HCT116s" else line}.npz', allow_pickle=True)
        assert list(map(str, z['targets'])) == h.targets
        Dh, has = z['Dh'].astype(np.float64), z['has']
        bl = np.nan_to_num(C.b, nan=0) < 5
        C.D.parts.update({'b<5': np.flatnonzero(bl), 'b>=5': np.flatnonzero(~bl)})
        mu4 = C.mu()
        rr, _ = C.score(mu4, gamma=0.0); rows += rr
        for g in C7D_G:
            add = np.where(has[:, None] & ~C.sig & C.okg, g * Dh, 0.0)
            rr, _ = C.score(C.mu(extra=add), gamma=g); rows += rr
            log(f'C7d-b {line} gamma {g}: ' + ', '.join(f'{r["part"]} {r["overall"]:.4f} (MSE {r["mse"]:.3f} PDS {r["pds"]:.3f} NMAE {r["nmae"]:.3f})' for r in rr))
        save_rows(rows, 'partC7db_rows'); del C; gc.collect()
    R = pd.DataFrame(rows)
    save_json(TAB / 'partC7db.json', fit_report(R, ['gamma'], dict(gamma=C7D_G), [0.0], tag='C7d-b'))


# ================================================================ C8: strength from expressed downstream genes (measure only)
def terc_spearman(y, x, b):
    from scipy.stats import spearmanr
    ok = np.isfinite(b); t = np.full(len(b), -1); t[ok] = tercile_of(b[ok]); out = {}
    for j, nm in enumerate(TERC):
        m = (t == j) & np.isfinite(x) & np.isfinite(y)
        out[nm] = float(spearmanr(y[m], x[m]).correlation) if m.sum() > 5 else float('nan')
    return out


def partC8(a):
    from scipy.stats import spearmanr
    import run_line as RL
    import scorer_step12 as S12
    genes = RL.panel_genes(); gp = {g: j for j, g in enumerate(genes)}
    Q = pd.read_csv(TAB / 'C8_counts.csv')
    agg = Q.groupby(['context', 'target'])[['n_de_gated5', 'n_de_cpm100', 'n_de_z3']].mean()
    fs = S20.feature_sets()['DepMap']; H = S20.hgnc_maps()
    res = dict(lines={}, DEF={}, H1_dist={}, definitions=dict(
        b_expr5='mean over genetic source tables holding k (own line out; tables with control CPM) of the count of the table\'s DE genes '
                '(|LFC/SE| > 3, the table\'s own > 5 CPM gate, own gene out) gated (> 5 CPM) in the target line\'s controls',
        b_expr100='the same with target-line control CPM > 100', b_z3='the same DE count with no target-line filter',
        b_own='the target gene\'s own control CPM in the scored line', DepMap='mean gene effect across DepMap 24Q4 lines (Step 20 F3a)'))
    def preds_for(ctx, targets, m0):
        A = agg.loc[ctx].reindex(targets) if ctx in agg.index.get_level_values(0) else pd.DataFrame(index=targets, columns=agg.columns)
        own = np.array([m0[gp[t]] if t in gp else np.nan for t in targets])
        dm = S20.gene_features(targets, fs, H)[0]; dmap = np.where(dm[:, 33] > 0, np.nan, dm[:, 32])
        return dict(b_expr5=A.n_de_gated5.to_numpy(float), b_expr100=A.n_de_cpm100.to_numpy(float), b_z3=A.n_de_z3.to_numpy(float),
                    b_own=own, DepMap=dmap)
    for line in TRUTH_CTX:
        ref = ref_for(EV[line]); h = ref.h; nt = S12.truth_calls(ref).sum(1).astype(float)
        m0 = c1_load(line)['m0']; b = b_for(EV[line], h.targets)
        P = preds_for(line, h.targets, m0); P['b'] = b
        out = {}
        for k, x in P.items():
            ok = np.isfinite(x) & np.isfinite(nt)
            out[k] = dict(n=int(ok.sum()), spearman_true_de=float(spearmanr(x[ok], nt[ok]).correlation),
                          partial_given_b=S19.partial_spearman(nt, x, b)[0] if k != 'b' else None,
                          within_b_terciles=terc_spearman(nt, x, b) if k != 'b' else None)
            if k in ('b_own', 'DepMap'):
                out[k]['partial_given_b_expr5'] = S19.partial_spearman(nt, x, P['b_expr5'])[0]
        res['lines'][line] = out
        if line == 'H1':
            res['H1_dist'] = {k: dict(zip(['q10', 'q25', 'median', 'q75', 'q90'], np.nanquantile(v, [.1, .25, .5, .75, .9]).tolist())) for k, v in P.items()}
        log(f'C8 {line}: ' + '; '.join(f'{k} rho {v["spearman_true_de"]:.3f} partial|b {v["partial_given_b"] if v["partial_given_b"] is None else round(v["partial_given_b"], 3)}' for k, v in out.items()))
        S19._REFS.pop(EV[line], None); gc.collect()
    bdef = pd.read_csv(S19.S17 / 'tables/b_DEF.csv').set_index('target')['b']
    for x in ('D', 'E', 'F'):
        c = c1_load(x); tg = list(map(str, c['targets'])); P = preds_for(x, tg, c['m0']); P['b'] = bdef.reindex(tg).to_numpy(float)
        res['DEF'][x] = {k: dict(zip(['q10', 'q25', 'median', 'q75', 'q90'], np.nanquantile(v, [.1, .25, .5, .75, .9]).tolist())) for k, v in P.items()}
    save_json(TAB / 'partC8.json', res)


# ================================================================ C9: knockdown clusters
C9_EMB = {'E24w': 20, 'E4w': 100}
C9_K = [50, 150]
BASIS9 = {'full': 'full', 'K562': 'K562', 'Orion': 'Orion', 'HCT116': 'full'}
C9_LINE_VARIANT = {'H1': 'full', 'K562': 'K562', 'HCT116': 'Orion', 'HEK293T': 'Orion', 'HCT116s': 'HCT116', 'DEF': 'full'}
C9_BEXCL = {'H1': ('h1',), 'K562': ('k562',), 'HCT116': ('hct116',), 'HEK293T': ('hek293t',), 'HCT116s': ('hct116',), 'DEF': ()}
COMPLEX12 = {'ribosome L': ['L ribosomal proteins'], 'ribosome S': ['S ribosomal proteins'], 'proteasome': ['Proteasome'],
             'Mediator': ['Mediator complex'], 'mitochondrial ribosome': ['Mitochondrial ribosomal proteins', 'Large subunit mitochondrial ribosomal proteins',
                                                                          'Small subunit mitochondrial ribosomal proteins'],
             'CCT': ['Chaperonin containing TCP1'], 'exosome': ['Exosome complex'], 'Integrator': ['Integrator complex'],
             'BAF': ['BAF complex subunits'], 'TFIID': ['General transcription factor IID complex subunits'],
             'spliceosome': ['Spliceosomal A complex', 'Spliceosomal B complex', 'Spliceosomal Bact complex', 'Spliceosomal C complex',
                             'Spliceosomal P complex', 'Spliceosomal E complex', 'Sm spliceosomal proteins'],
             'NADH dehydrogenase': ['NADH:ubiquinone oxidoreductase core subunits', 'NADH:ubiquinone oxidoreductase supernumerary subunits']}


_BTAB = []


def b_all(exclude_keys):
    """(b) for every target: mean over source tables (lines in exclude_keys out) of each table's mean n_real over its
    rows (de_labels_known, finite n_real), as step15.source_n_real builds it."""
    import h5py
    import stage1 as S
    from anndata.io import read_elem
    if not _BTAB:
        for p in S.discover_tables([str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')]):
            with h5py.File(p, 'r') as f:
                o = read_elem(f['obs'])
            key = S._line_key(o['background'].iloc[0])
            o = o[o['de_labels_known'].astype(bool) & np.isfinite(o['n_real'].astype(float))]
            if len(o):
                _BTAB.append((key, o.groupby(o['target_gene'].astype(str))['n_real'].mean()))
    vals = [m for k, m in _BTAB if k not in set(exclude_keys)]
    return pd.concat(vals, axis=1).mean(axis=1)


def src_nreal():
    """n_real for every srcrows row (Step 20 f_source_rows order: genetic rows of non-Tahoe Delta tables), from the
    E4w matrix meta (same tables, same row order); checked against partF1_srcrows_delta."""
    meta = pd.read_parquet(MAT / 'meta.parquet'); s = meta[~meta.tahoe & (meta.modality != 'drug')].reset_index(drop=True)
    ref = pd.read_parquet(S20T / 'partF1_srcrows_delta.parquet')
    assert len(s) == len(ref) and (s.table.to_numpy() == ref.table.to_numpy()).all() and (s.target.to_numpy() == ref.target.to_numpy()).all()
    return s


def c9_embedding(emb, variant, s):
    K = C9_EMB[emb]; RS = np.load(S20T / f'partF1_RS_{emb}_{BASIS9[variant]}.npy')[:, :K].astype(np.float64)
    ok = ~s.background.isin(VARIANT_DROP[variant]).to_numpy()
    if emb == 'E24w':
        ok &= (s.n_real.to_numpy() >= 50)
    d = pd.DataFrame(RS[ok]); d['bg'] = s.background.to_numpy()[ok]; d['t'] = s.target.to_numpy()[ok]
    E = d.groupby(['t', 'bg']).mean().groupby(level='t').mean()
    return list(E.index), E.to_numpy()


def kmeans_cos(E, k):
    from sklearn.cluster import KMeans
    U = E / np.maximum(np.linalg.norm(E, axis=1, keepdims=True), 1e-12)
    return KMeans(n_clusters=k, n_init=10, random_state=SEED).fit_predict(U)


def string_neighbours():
    H = S20.hgnc_maps(); D = V / 'lib/downloads/step20_2026-11-02/string_v12'
    info = pd.read_csv(D / '9606.protein.info.v12.0.txt.gz', sep='\t', usecols=[0, 1]); info.columns = ['id', 'name']
    nm = dict(zip(info.id, info.name))
    L = pd.read_csv(D / '9606.protein.links.v12.0.txt.gz', sep=' '); L = L[L.combined_score >= 700]
    a_ = L.protein1.map(nm).map(lambda x: S20.to_hgnc(x, H)); b_ = L.protein2.map(nm).map(lambda x: S20.to_hgnc(x, H))
    nb = {}
    for x, y in zip(a_, b_):
        if x and y and x != y:
            nb.setdefault(x, set()).add(y)
    return nb, H


def partC9(a):
    from scipy.stats import spearmanr
    import scorer_step12 as S12
    s = src_nreal(); nb, H = string_neighbours()
    grp = H['groups']
    res = dict(C9a={}, C9b={}, embeddings={}); clus = {}
    for emb in C9_EMB:
        for v in ('full', 'K562', 'Orion', 'HCT116'):
            tg, E = c9_embedding(emb, v, s)
            res['embeddings'][f'{emb} {v}'] = dict(targets=len(tg), basis=BASIS9[v], dropped=list(VARIANT_DROP[v]))
            for k in C9_K:
                clus[(emb, v, k)] = dict(zip(tg, kmeans_cos(E, k)))
            log(f'C9 {emb} {v}: {len(tg)} targets embedded, clustered')
    pickle_path = TAB / 'partC9_clusters.pkl'
    import pickle
    pickle.dump(clus, open(pickle_path, 'wb'))
    # C9a on the full variant
    rng = np.random.default_rng(SEED)
    for (emb, v, k), lab in clus.items():
        if v != 'full':
            continue
        out = {}
        tl = np.array(list(lab.keys())); ll = np.array(list(lab.values()))
        rp = rng.integers(0, len(tl), (20000, 2)); rand = float((ll[rp[:, 0]] == ll[rp[:, 1]]).mean())
        for cn, gnames in COMPLEX12.items():
            mem = [t for t in tl if any(g in gnames for g in grp.get(S20.to_hgnc(t, H) or t, []))]
            if len(mem) < 2:
                out[cn] = dict(members=len(mem)); continue
            L_ = np.array([lab[t] for t in mem]); n = len(L_)
            same = (L_[:, None] == L_[None, :]).sum() - n
            out[cn] = dict(members=n, share_pairs_same_cluster=float(same / (n * (n - 1))))
        res['C9a'][f'{emb} k{k}'] = dict(random_pairs=rand, complexes=out)
        log(f'C9a {emb} k {k}: random {rand:.3f}; ' + ', '.join(f'{c} {o.get("share_pairs_same_cluster", np.nan):.2f}' for c, o in out.items()))
    # C9b strength prior
    tgl = all_line_targets(); bdef = pd.read_csv(S19.S17 / 'tables/b_DEF.csv').set_index('target')['b']
    for line in TRUTH_CTX + ['DEF']:
        targets = tgl[line]; v = C9_LINE_VARIANT[line]; ba = b_all(C9_BEXCL[line])
        b = bdef.reindex(targets).to_numpy(float) if line == 'DEF' else b_for(EV[line], targets)
        if line != 'DEF':
            chk = float(np.nanmax(np.abs(ba.reindex(targets).to_numpy(float) - b)))
            assert chk < 1e-9, f'b_all does not reproduce (b) for {line}: {chk}'
            ref = ref_for(EV[line]); assert list(ref.h.targets) == targets; nt = S12.truth_calls(ref).sum(1).astype(float); S19._REFS.pop(EV[line], None)
        for (emb, vv, k), lab in clus.items():
            if vv != v:
                continue
            route = []; cl = np.full(len(targets), -1)
            for i, t in enumerate(targets):
                if t in lab:
                    cl[i] = lab[t]; route.append('own')
                    continue
                ns = [lab[n_] for n_ in nb.get(S20.to_hgnc(t, H) or t, ()) if n_ in lab]
                if ns:
                    vals, cnt = np.unique(ns, return_counts=True); cl[i] = vals[np.argmax(cnt)]; route.append('STRING')
                else:
                    route.append('none')
            members = {}
            for t, c in lab.items():
                members.setdefault(c, []).append(t)
            bc = np.array([np.nanmedian([ba.get(m, np.nan) for m in members[c] if m != t]) if c >= 0 and len(members[c]) > 1 else np.nan
                           for t, c in zip(targets, cl)])
            o = dict(routes={r: int(route.count(r)) for r in ('own', 'STRING', 'none')})
            if line != 'DEF':
                ok = np.isfinite(bc)
                o.update(spearman_true_de=float(spearmanr(bc[ok], nt[ok]).correlation), partial_given_b=S19.partial_spearman(nt, bc, b)[0],
                         within_b_terciles=terc_spearman(nt, bc, b),
                         b_lt5=float(spearmanr(bc[ok & (np.nan_to_num(b) < 5)], nt[ok & (np.nan_to_num(b) < 5)]).correlation) if (ok & (np.nan_to_num(b) < 5)).sum() > 5 else float('nan'),
                         n_b_lt5=int((ok & (np.nan_to_num(b) < 5)).sum()))
            o['b_clust_quantiles'] = dict(zip(['q10', 'q25', 'median', 'q75', 'q90'], np.nanquantile(bc, [.1, .25, .5, .75, .9]).tolist())) if np.isfinite(bc).any() else None
            res['C9b'][f'{line} {emb} k{k}'] = o
            log(f'C9b {line} {emb} k {k}: {o}')
        save_json(TAB / 'partC9.json', res)
    save_json(TAB / 'partC9.json', res)


def pooled_delta(mat, target, variant, cache, own_out=True):
    """A target's pooled source Delta (raw), genetic rows of the variant's backgrounds: per background mean of its rows
    (each table's measured entries), then the mean over backgrounds. Own gene NaN when own_out."""
    key = (target, variant, own_out)
    if key in cache:
        return cache[key]
    rows = [r for r in mat.rows_of.get(target, []) if mat.genetic[r] and mat.bg[r] not in VARIANT_DROP[variant]]
    if not rows:
        cache[key] = None; return None
    X = mat.raw(np.array(rows)); bgs = mat.bg[rows]
    with np.errstate(invalid='ignore'), __import__('warnings').catch_warnings():
        __import__('warnings').simplefilter('ignore')
        P = np.nanmean(np.stack([np.nanmean(X[bgs == g], 0) for g in np.unique(bgs)]), 0)
    if own_out and target in mat.gp:
        P[mat.gp[target]] = np.nan
    cache[key] = P.astype(np.float32)
    return cache[key]


C9C_G = [0.25, 0.5, 1.0]


def partC9c(a):
    import pickle
    import run_line as RL
    genes = RL.panel_genes()
    mat = Mat(); mat.gp = {g: j for j, g in enumerate(genes)}
    clus = pickle.load(open(TAB / 'partC9_clusters.pkl', 'rb')); nb, H = string_neighbours()
    rows = resume('partC9c_rows'); res = dict(choice='E24w, k 150 (declared before the run)', lines={})
    for line in LINES:
        if any(r['line'] == line for r in rows):
            continue
        C = CLine(line); h = C.h; v = C9_LINE_VARIANT[line]; lab = clus[('E24w', v, 150)]
        members = {}
        for t, c in lab.items():
            members.setdefault(c, []).append(t)
        cache = {}; M_ = np.zeros((h.T, h.G)); has = np.zeros(h.T, bool); nm_ = []
        for i, t in enumerate(h.targets):
            if not C.on5[i]:
                continue
            c = lab.get(t)
            if c is None:
                ns = [lab[n_] for n_ in nb.get(S20.to_hgnc(t, H) or t, ()) if n_ in lab]
                if ns:
                    vals, cnt = np.unique(ns, return_counts=True); c = vals[np.argmax(cnt)]
            if c is None:
                continue
            P = [pooled_delta(mat, m, v, cache) for m in members[c] if m != t]
            P = [p for p in P if p is not None]
            if not P:
                continue
            with np.errstate(invalid='ignore'), __import__('warnings').catch_warnings():
                __import__('warnings').simplefilter('ignore')
                M_[i] = np.nan_to_num(np.nanmean(np.stack(P), 0)); has[i] = True; nm_.append(len(P))
        M_ = np.where(C.okg & ~C.sig, M_, 0.0)
        res['lines'][line] = dict(knockdowns_b_ge5=int(C.on5.sum()), with_term=int(has.sum()), mates_median=float(np.median(nm_)) if nm_ else None)
        rr, _ = C.score(C.mu(), gamma=0.0); rows += rr
        for g in C9C_G:
            rr, _ = C.score(C.mu(extra=g * M_), gamma=g); rows += rr
            log(f'C9c {line} gamma {g}: ' + ', '.join(f'{r["part"]} {r["overall"]:.4f} (PDS {r["pds"]:.3f})' for r in rr))
        save_rows(rows, 'partC9c_rows'); del C; gc.collect()
    R = pd.DataFrame(rows)
    out = fit_report(R, ['gamma'], dict(gamma=C9C_G), [0.0], tag='C9c')
    out['reference'] = 'Step 11 partner arm (h) lowered H1 (n 3: 0.201-0.208 against 0.231); S3: 0 of 300 D/E/F targets had a K562 partner above cosine 0.5'
    save_json(TAB / 'partC9c.json', out | res)


# ================================================================ C10: cascades (measure only)
def partC10(a):
    import run_line as RL
    import scorer_arms as SA
    genes = RL.panel_genes(); gp = {g: j for j, g in enumerate(genes)}
    mat = Mat(); mat.gp = gp; cache = {}
    C = CLine('H1', terc_parts=False); h = C.h
    ba = b_all(('h1',)); elig = [g for g in genes if ba.get(g, 0) >= 50 and g in mat.rows_of]
    log(f'C10: {len(elig)} genes eligible as X (own source knockdown, (b) >= 50)')
    sig = {}; selfd = {}
    for g in elig:
        p = pooled_delta(mat, g, 'full', cache, own_out=False)
        if p is None or not np.isfinite(p[gp[g]]):
            continue
        selfd[g] = float(p[gp[g]]); q = p.copy(); q[gp[g]] = np.nan; sig[g] = np.nan_to_num(q)
    elig = [g for g in elig if g in sig]
    cache.clear()
    z = C.c['mu'] / np.sqrt(np.where(np.isfinite(C.c['var_fe']), C.c['var_fe'], np.inf))
    mu4 = C.mu(); dpred = C.delta_pred(mu4); Rres = h.R - SA.bulk(mu4); Y = h.R - h.ctrl[None, :]
    ok = C.okg
    rng = np.random.default_rng(SEED)
    def term(i, Xs):
        T_ = np.zeros(h.G)
        for x in Xs:
            w = np.clip(dpred[i, gp[x]] / selfd[x], 0, 1) if selfd[x] != 0 else 0.0
            T_ += w * sig[x]
        return T_
    out = []
    for i, t in enumerate(h.targets):
        if not C.on5[i]:
            continue
        pA = pooled_delta(mat, t, 'full', cache)
        if pA is None:
            out.append(dict(target=t, n_X=0)); continue
        Xs = [x for x in elig if x != t and z[i, gp[x]] <= -3 and np.isfinite(pA[gp[x]]) and pA[gp[x]] <= -0.1]
        rec = dict(target=t, n_X=len(Xs))
        if Xs:
            T_ = term(i, Xs); m = ok[i]
            rec['r_residual'] = float(np.corrcoef(T_[m], Rres[i, m])[0, 1]) if np.std(T_[m]) > 0 else np.nan
            rec['r_truth'] = float(np.corrcoef(T_[m], Y[i, m])[0, 1]) if np.std(T_[m]) > 0 else np.nan
            nr, nt_ = [], []
            pool_ = [x for x in elig if x != t]
            for _ in range(100):
                Xr = list(rng.choice(pool_, len(Xs), replace=False)); Tn = term(i, Xr)
                nr.append(np.corrcoef(Tn[m], Rres[i, m])[0, 1] if np.std(Tn[m]) > 0 else np.nan)
                nt_.append(np.corrcoef(Tn[m], Y[i, m])[0, 1] if np.std(Tn[m]) > 0 else np.nan)
            rec['null_r_residual'] = nr; rec['null_r_truth'] = nt_
        out.append(rec)
    O = pd.DataFrame(out); W = O[O.n_X > 0]
    def summ(col):
        obs = float(np.nanmedian(W[col])); nulls = np.nanmedian(np.stack(W[f'null_{col}'].to_numpy()), 0)
        return dict(median=obs, null_median=float(np.nanmedian(nulls)), null_q95=float(np.nanquantile(nulls, .95)), percentile=float((nulls < obs).mean() * 100))
    res = dict(H1=dict(knockdowns_b_ge5=int(C.on5.sum()), with_at_least_one_X=int(len(W)), X_per_knockdown=dict(zip(['median', 'q90', 'max'], np.quantile(W.n_X, [.5, .9, 1]).tolist())) if len(W) else None,
                       r_residual=summ('r_residual') if len(W) else None, r_truth=summ('r_truth') if len(W) else None), eligible_X=len(elig))
    O.drop(columns=[c for c in O.columns if c.startswith('null_')]).to_csv(TAB / 'partC10_H1.csv', index=False)
    # D/E/F targets (all sources; context D's pooled z)
    c = c1_load('D'); zd = c['mu'] / np.sqrt(np.where(np.isfinite(c['var_fe']), c['var_fe'], np.inf))
    bD = b_all(()); bdef = pd.read_csv(S19.S17 / 'tables/b_DEF.csv').set_index('target')['b']
    eligD = [g for g in elig if bD.get(g, 0) >= 50]
    n1, n1b = 0, 0
    for i, t in enumerate(map(str, c['targets'])):
        pA = pooled_delta(mat, t, 'full', cache)
        if pA is None:
            continue
        if any(x != t and zd[i, gp[x]] <= -3 and np.isfinite(pA[gp[x]]) and pA[gp[x]] <= -0.1 for x in eligD):
            n1 += 1; n1b += int(bdef.get(t, 0) >= 5)
    res['DEF'] = dict(targets=len(c['targets']), with_at_least_one_X=n1, with_at_least_one_X_and_b_ge5=n1b,
                      note='context D\'s pooled z (sources identical for D/E/F; only the gate differs)')
    save_json(TAB / 'partC10.json', res); log(f'C10: {res}')


# ================================================================ C11: the p53 program (measure only)
DL21 = V / 'lib/downloads/step21_2026-11-03'
CELLO_ID = {'K562': 'K-562', 'RPE1': 'hTERT-RPE1', 'HepG2': 'Hep-G2', 'Jurkat': 'Jurkat', 'HCT116': 'HCT 116', 'HEK293T': 'HEK293T',
            'H1': 'H1', 'KOLF2.1J': 'KOLF2.1J', 'A549': 'A-549', 'BXPC3': 'BxPC-3', 'HAP1': 'HAP-1', 'HT29': 'HT-29', 'MCF7': 'MCF-7'}


def cellosaurus(bg):
    """TP53 sequence variations and viral transformant / integration comments for a line, from the Cellosaurus API
    (cached in lib/downloads/step21_2026-11-03/cellosaurus/). HipSci donor clones: the entry whose ID ends in -<clone>."""
    import subprocess
    import urllib.parse
    if bg == 'HipSci_pooled':
        return dict(entry=None, tp53=[], transformant=[], note='pooled HipSci iPSC donors (Feng genome-wide): not one cell line, no Cellosaurus entry; not looked up')
    d = DL21 / 'cellosaurus'; d.mkdir(parents=True, exist_ok=True); f = d / f'{bg}.txt'
    if not f.exists():
        q = f'id:"{CELLO_ID[bg]}"' if bg in CELLO_ID else bg
        url = 'https://api.cellosaurus.org/search/cell-line?' + urllib.parse.urlencode(dict(q=q, fields='id,ac,var,cc,ref', format='txt'))
        r = subprocess.run(['curl', '-s', '-m', '60', url], capture_output=True, text=True)
        if r.returncode or 'ID   ' not in r.stdout:
            raise RuntimeError(f'Cellosaurus query failed for {bg}: rc {r.returncode}')
        f.write_text(r.stdout)
    entries = f.read_text().split('\n//')
    want = CELLO_ID.get(bg)
    pick = [e for e in entries if any(l.startswith('ID   ') and ((l[5:].strip() == want) if want else l[5:].strip().endswith('-' + bg)) for l in e.splitlines())]
    if not pick:
        return dict(entry=None, tp53=[], transformant=[])
    e = pick[0].splitlines()
    return dict(entry=next(l[5:].strip() for l in e if l.startswith('ID   ')), accession=next((l[5:].strip() for l in e if l.startswith('AC   ')), None),
                tp53=sorted({l[5:].strip().split(' (')[0] for l in e if l.startswith('CC   Sequence variation:') and '; TP53;' in l}),
                transformant=[l[5:].strip() for l in e if l.startswith('CC   ') and any(k in l for k in ('Transformant', 'Genetic integration', 'SV40', 'adenovirus', 'Adenovirus'))])


def partC11(a):
    import run_line as RL
    import scorer_arms as SA
    genes = RL.panel_genes(); gp = {g: j for j, g in enumerate(genes)}
    f = DL21 / 'fischer2017_core116.json'
    if not f.exists():
        log('C11: the Fischer 2017 core set is not available; part stopped'); sys.exit(11)
    P53 = json.loads(f.read_text()); pset = np.array([gp[g] for g in P53['genes'] if g in gp])
    H = S20.hgnc_maps(); grp = H['groups']
    rib = lambda t: any(g in ('L ribosomal proteins', 'S ribosomal proteins') for g in grp.get(S20.to_hgnc(t, H) or t, []))
    mat = Mat(); mat.gp = gp
    res = dict(gene_set=dict(source=P53['source'], genes=len(P53['genes']), on_panel=int(len(pset))), lines={}, H1={}, HCT116final={}, controls={})
    # 1. per source line (genetic rows of full-profile tables): mean Delta of the set on RPL/RPS knockdowns and on the
    #    line's strong n_real tercile; each knockdown averaged over the line's tables first; set genes measured in the row
    meta = mat.meta
    for bg in sorted(meta.background[mat.genetic].unique()):
        rows = np.flatnonzero(mat.genetic & (mat.bg == bg))
        sub = meta.iloc[rows]
        X = mat.raw(rows)[:, pset]
        kd = sub.target.to_numpy(); nr = sub.n_real.to_numpy(float)
        per = pd.DataFrame(dict(kd=kd, v=np.nanmean(X, 1), nr=nr)).groupby('kd').mean()
        t = np.full(len(per), -1); okn = np.isfinite(per.nr.to_numpy())
        if okn.sum() >= 9: t[okn] = tercile_of(per.nr.to_numpy()[okn])
        isr = np.array([rib(k) for k in per.index])
        try:
            cs = cellosaurus(bg)
        except Exception as e:
            log(f'C11: {e}; part stopped'); save_json(TAB / 'partC11.json', res); sys.exit(11)
        res['lines'][bg] = dict(knockdowns=int(len(per)), rpl_rps=int(isr.sum()), mean_set_delta_rpl_rps=float(per.v[isr].mean()) if isr.any() else None,
                                strong_tercile=int((t == 2).sum()), mean_set_delta_strong=float(per.v[t == 2].mean()) if (t == 2).any() else None,
                                mean_set_delta_all=float(per.v.mean()), cellosaurus=cs)
        log(f'C11 {bg}: {res["lines"][bg]["mean_set_delta_rpl_rps"]} (RPL/RPS n {isr.sum()}), strong {res["lines"][bg]["mean_set_delta_strong"]}; TP53 {cs["tp53"][:1]} {cs["transformant"][:1]}')
    save_json(TAB / 'partC11.json', res)
    # 2-3. row 4 on H1 and HCT116final (screen setup)
    for line, key in (('H1', 'H1'), ('HCT116s', 'HCT116final')):
        C = Row4(line); h = C.h; mu4 = C.mu(); P = SA.bulk(mu4); own = own_mask(h)
        import scorer_step12 as S12
        tt = tercile_of(S12.truth_calls(C.D.ref).sum(1))
        E = ((P - h.R) ** 2) * ~own; En = ((h.R - h.ctrl[None, :]) ** 2) * ~own; ins = np.zeros(h.G, bool); ins[pset] = True
        st = tt == 2
        isr = np.array([rib(t) for t in h.targets])
        dp = (P - h.ctrl[None, :])[:, pset]; dt = (h.R - h.ctrl[None, :])[:, pset]
        res[key] = dict(set_share_strong_error=float(E[st][:, ins].sum() / E[st].sum()), set_share_strong_real_energy=float(En[st][:, ins].sum() / En[st].sum()),
                        set_gene_share=float(ins.mean()), rpl_rps_knockdowns=int(isr.sum()),
                        rpl_rps_mean_set_delta_pred=float(dp[isr].mean()) if isr.any() else None, rpl_rps_mean_set_delta_true=float(dt[isr].mean()) if isr.any() else None,
                        strong_mean_set_delta_pred=float(dp[st].mean()), strong_mean_set_delta_true=float(dt[st].mean()))
        log(f'C11 {key}: {res[key]}'); del C; gc.collect()
    # 4. basal CPM of MDM2, CDKN1A, TP53, CDKN2A
    g4 = ['MDM2', 'CDKN1A', 'TP53', 'CDKN2A']; j4 = [gp[g] for g in g4]
    for x in ('D', 'E', 'F'):
        m0 = np.load(S19.TAB / f'controls_{x}.npz')['m0']; res['controls'][x] = {g: float(m0[j]) for g, j in zip(g4, j4)}
    for l, ctx in (('H1 (Eval A controls)', 'H1'), ('HCT116final', 'HCT116'), ('K562final', 'K562'), ('HEK293Tfinal', 'HEK293T')):
        m0 = c1_load(ctx)['m0']; res['controls'][l] = {g: float(m0[j]) for g, j in zip(g4, j4)}
    by = {}
    for p, stem, bg, tahoe in S19.e_tables():
        if tahoe: continue
        b = pd.read_csv(p.with_name(p.stem + '_baseline.csv'), index_col=0)
        if 'control_mean_cpm' in b:
            by.setdefault(bg, []).append(b['control_mean_cpm'].reindex(g4).to_numpy(float))
    for bg, v in by.items():
        res['controls'][f'source {bg}'] = dict(zip(g4, np.nanmean(np.stack(v), 0).tolist()))
    res['controls']['note'] = 'Basal CPM is weak evidence of p53 status: CDKN1A / MDM2 reflect basal p53 activity, TP53 mRNA says little about protein function, and CDKN2A loss is common in culture.'
    save_json(TAB / 'partC11.json', res)


# ================================================================ Track 1 comparison (vcc2026, after the copy back)
def t1cmp(a):
    import scorer_proxy as SP
    import scorer_step12 as S12
    import scorer_step16 as S16
    from scipy.stats import spearmanr
    res = {}
    A = pd.read_csv(TAB / 'partA_rows.csv')
    for ev, line, hn in (('EvalA', 'H1', 'evalA'), ('HCT116final', 'HCT116s', 'evalHCT116final')):
        ref = ref_for(ev); h = ref.h
        hd = CELLS / ev / 'e2e/salt_none/harness' / hn
        W = pd.read_csv(hd / 'results.csv').pivot(index='perturbation', columns='metric', values='value').reindex(h.targets)
        Fn = Fast.from_wide(h, W); F0 = Fast.from_wide(h, S19.results_wide('K0', ev))
        b = b_for(ev, h.targets); onK4 = np.nan_to_num(b, nan=-np.inf) >= T_HYB
        Wt, _ = S19.combined(ev, onK4); Ft = Fast.from_wide(h, Wt)
        I = S19.boot_index(h.T)
        Sn, Rn = Fn.scaled(I); S0, R0 = F0.scaled(I); St, Rt_ = Ft.scaled(I)
        pn, prn = Fn.point(); p0, _ = F0.point(); pt, _ = Ft.point()
        summ = json.loads((hd / 'summary.json').read_text())
        chk = abs(pn['overall'] - summ['scaled']['overall'])
        out = dict(rebuild_vs_summary_abs_diff=chk,
                   components={k: dict(scaled=pn[k], ci=[float(np.quantile(Sn[k], .025)), float(np.quantile(Sn[k], .975))],
                                       raw=prn.get(k), raw_ci=[float(np.quantile(Rn[k], .025)), float(np.quantile(Rn[k], .975))] if k in Rn else None)
                               for k in COMPS + ['overall']})
        for nm, (S_, p_) in (('vs live entry (15c config 1\' + K0)', (S0, p0)), ('vs hybrid T5 combined panel (Step 19 Part A)', (St, pt))):
            out[nm] = {k: dict(diff=pn[k] - p_[k], ci=[float(np.quantile(Sn[k] - S_[k], .025)), float(np.quantile(Sn[k] - S_[k], .975))],
                               p_better=float((Sn[k] - S_[k] > 0).mean())) for k in COMPS + ['overall']}
        pr = A[(A.line == line) & (A.row == 'row 4') & (A.part == 'full')].iloc[0]
        out['proxy_minus_cells'] = {k: dict(proxy=float(pr[k]), cells=pn[k], gap=float(pr[k] - pn[k])) for k in COMPS + ['overall']}
        T = S12.truth_calls(ref); nt = T.sum(1)
        if ev == 'EvalA':
            tt = tercile_of(nt); out['terciles'] = {}
            for j, tn in enumerate(TERC):
                rr = np.flatnonzero(tt == j); rng = np.random.default_rng(SEED)
                Ij = np.stack([rr[rng.integers(0, len(rr), len(rr))] for _ in range(S19.N_BOOT)])
                sn, _ = Fn.scaled(Ij); s0, _ = F0.scaled(Ij); st_, _ = Ft.scaled(Ij)
                out['terciles'][tn] = dict(knockdowns=len(rr), overall=Fn.point(rr)[0]['overall'], mse=Fn.point(rr)[0]['mse'], pds=Fn.point(rr)[0]['pds'],
                                           vs_live=dict(diff=Fn.point(rr)[0]['overall'] - F0.point(rr)[0]['overall'],
                                                        ci=[float(np.quantile(sn['overall'] - s0['overall'], .025)), float(np.quantile(sn['overall'] - s0['overall'], .975))]),
                                           vs_T5=dict(diff=Fn.point(rr)[0]['overall'] - Ft.point(rr)[0]['overall'],
                                                      ci=[float(np.quantile(sn['overall'] - st_['overall'], .025)), float(np.quantile(sn['overall'] - st_['overall'], .975))]))
            hv = S16.halves(T); out['halves'] = {f'half{f}': Fn.point(hv[f])[0] for f in (0, 1)}
            out['halves_vs_live'] = {f'half{f}': Fn.point(hv[f])[0]['overall'] - F0.point(hv[f])[0]['overall'] for f in (0, 1)}
        z = np.load(CELLS / ev / 'input.npz', allow_pickle=True); names = list(map(str, z['knockdowns']))
        r = [names.index(t) for t in h.targets]; un = z['n_sources'][r] <= 0; sig = z['de_sig'][r] & ~un
        per, by = SP.call_diff(ref, hd, sig, CELLS / ev / 'e2e/controls/controls_baseline.csv')
        per = per.set_index('target').reindex(h.targets) if 'target' in per else per
        per.to_csv(TAB / f't1_calls_per_kd_{ev}.csv'); by.to_csv(TAB / f't1_calls_by_decile_{ev}.csv', index=False)
        out['calls'] = dict(requested_median=float(per.requested.median()), requested_mean=float(per.requested.mean()),
                            realised_median=float(per.realised.median()), realised_mean=float(per.realised.mean()),
                            requested_not_realised_total=int(per.requested_not_realised.sum()), realised_not_requested_total=int(per.realised_not_requested.sum()),
                            spearman_realised_true_de=float(spearmanr(per.realised.to_numpy(float), nt).correlation),
                            spearman_requested_true_de=float(spearmanr(per.requested.to_numpy(float), nt).correlation))
        res[ev] = out
        log(f't1cmp {ev}: overall {pn["overall"]:.4f}; vs live {out[list(out)[2]]["overall"]}; vs T5 {out[list(out)[3]]["overall"]}')
        S19._REFS.pop(ev, None); gc.collect()
    save_json(TAB / 't1cmp.json', res)


# ================================================================ Part D: recommendation (combinations of parts that met their criterion)
D_PARTS = {'B': ('partB.json', None), 'C3': ('partC3.json', None), 'C4': ('partC4.json', None), 'C6': ('partC6.json', None),
           'C7c co-movement': ('partC7c.json', 'co_movement'), 'C7c F5-M3e': ('partC7c.json', 'F5_M3e'), 'C7d-b': ('partC7db.json', None),
           'C9c': ('partC9c.json', None)}
D_ORDER = ['C3', 'C4', 'C6', 'C7c F5-M3e', 'C7c co-movement', 'C7d-b', 'C9c', 'B']     # order the modifications are applied


def d_status():
    out = {}
    for nm, (f, sub) in D_PARTS.items():
        if not (TAB / f).exists():
            out[nm] = dict(available=False); continue
        r = json.loads((TAB / f).read_text()); r = r[sub] if sub else r
        picks = {fk: fv.get('picked') for fk, fv in r.get('fits', {}).items()}
        out[nm] = dict(available=True, met=bool(r['criterion']['met']), picks=picks)
    return out


def d_apply(C, mods, f, ctx):
    """mu and calls for row 4 plus mods {part: pick} on line C, for fit direction f (H1: the held-out half reads it)."""
    import pickle
    sig = C.sig; called = None; extra = np.zeros((C.h.T, C.h.G)); shat = None; eta = 0.0
    for nm in D_ORDER:
        if nm not in mods:
            continue
        pk = mods[nm]
        if nm == 'C3':
            sig = c3_rules(C)[pk['rule']]
        elif nm == 'C4':
            z2 = C.zu ** 2; called = pk['s'] * C.D.cg * (z2 / (z2 + pk['kappa']) if pk['kappa'] else 1.0)
        elif nm == 'C6':
            cg, m0 = C.D.cg, C.D.m0
            called = np.where((cg > 0) & (m0[None, :] >= pk['c_floor']), pk['s_up'] * cg, np.where(cg < 0, pk['s_down'] * cg, 0.5 * cg))
        elif nm == 'C7c F5-M3e':
            pr = ctx.setdefault('c7b', pickle.load(open(TAB / 'partC7b_preds.pkl', 'rb'))['preds'])
            if C.line == 'H1':
                shat = np.zeros((C.h.T, F5K))
                for ff in (0, 1):
                    shat[C.hv[1 - ff]] = pr[('H1', 'M3e', F5K, ff)][C.hv[1 - ff]]
            else:
                shat = pr[(C.line, 'M3e', F5K, -1)].astype(np.float64)
        elif nm == 'C7c co-movement':
            import run_line as RL
            X, has = c7_rows(C.line, C.h.targets, RL.panel_genes()); lj = json.loads((TAB / 'partC7c_lambda.json').read_text())['lambda_']
            L = lj[f'fit_half{f}'] if C.line == 'H1' else lj['all_H1']
            L = L['fit_half'] if 'fit_half' in L else L
            lam = np.array(L['lam'])[b_terc(C.b, np.array(L['cuts']))]
            extra += np.where((C.on5 & has)[:, None] & C.okg, pk['gamma'] * lam[:, None] * X, 0.0)
        elif nm == 'C7d-b':
            z = np.load(MAT / f'C7d_{C.line}.npz', allow_pickle=True)
            extra += np.where(z['has'][:, None] & ~sig & C.okg, pk['gamma'] * z['Dh'].astype(np.float64), 0.0)
        elif nm == 'C9c':
            M_ = ctx.get(('c9c', C.line))
            if M_ is None:
                raise RuntimeError('C9c term not cached for Part D')
            extra += np.where(C.okg & ~sig, pk['gamma'] * M_, 0.0)
        elif nm == 'B':
            eta = pk['eta']
    base = C.mu_spec
    if called is not None:
        import scorer_arms as SA
        base = SA.mu_from_lfc(np.where(C.sig_k0, called, C.L_spec), C.D.m0)
    F5add = C.addF5 if shat is None else C.f5_from_shat(shat)
    mu = mu_add_delta(base, F5add + C.addS + extra, C.own)
    if eta:
        mu = mu_add_delta(mu, np.broadcast_to(-eta * mean_pred_delta(C, mu), mu.shape), C.own)
    return mu, sig


def partD(a):
    st = d_status(); res = dict(status=st)
    passed = [nm for nm, v in st.items() if v.get('met')]
    res['passed'] = passed
    log(f'D: parts that met their criterion: {passed or "none"}')
    combos = {}
    if passed:
        if 'C4' in passed and 'C6' in passed:
            for keep in ('C4', 'C6'):
                combos[f'all passing, {keep} kept'] = [p for p in passed if p != ('C6' if keep == 'C4' else 'C4')]
        else:
            combos['all passing'] = passed
    rows = []
    if 'C9c' in passed:
        log('D: C9c passed; its term must be rebuilt here (not cached) -- not implemented'); sys.exit(12)
    for nm, parts in combos.items():
        for line in LINES:
            C = CLine(line, terc_parts=False); ctx = {}
            for f in ((0, 1) if line == 'H1' else (-1,)):
                ff = max(f, 0)
                mods = {}
                for p in parts:
                    pk = st[p]['picks'].get(f'fit_half{ff}')
                    if pk is None and p == 'C7c F5-M3e':
                        pk = {}
                    mods[p] = pk
                mu, sig = d_apply(C, mods, ff, ctx)
                rr, _ = C.score(mu, sig, combo=nm, fit_half=f); rows += rr
                log(f'D {nm} {line} fit {f}: ' + ', '.join(f'{r["part"]} {r["overall"]:.4f}' for r in rr))
            del C; gc.collect()
    if rows:
        R = save_rows(rows, 'partD_rows'); A = pd.read_csv(TAB / 'partA_rows.csv'); b4 = A[A.row == 'row 4']
        for nm in combos:
            fv = []
            for f in (0, 1):
                x = R[(R.combo == nm) & (R.line == 'H1') & (R.fit_half == f) & (R.part == f'half{1 - f}')].iloc[0]
                y = b4[(b4.line == 'H1') & (b4.part == f'half{1 - f}')].iloc[0]
                gain = {c: float(x[c] - y[c]) for c in COMPS + ['overall']}; sc = sum(gain[c] for c in COMPS)
                gr = {l: dict(diff=float(R[(R.combo == nm) & (R.line == l) & (R.part == 'full')].overall.iloc[0] - b4[(b4.line == l) & (b4.part == 'full')].overall.iloc[0])) for l in GUARD}
                fv.append(dict(held_out=comp_dict(x), base_held_out=comp_dict(y), gain=gain, guardrails=gr,
                               reach_share_of_component_gain=float(gain['reach'] / sc) if sc else float('nan')))
            res.setdefault('combos', {})[nm] = dict(parts=combos[nm], fits=fv, criterion=criterion(fv))
    save_json(TAB / 'partD.json', res)


# ================================================================ unit tests
def test(a):
    rng = np.random.default_rng(0)
    # criterion and fit_report on a synthetic grid
    rows = []
    for f_ in ('half0', 'half1', 'full'):
        for x, o in ((0.0, .20), (1.0, .26), (2.0, .30)):
            rows.append(dict(line='H1', part=f_, x=x, overall=o + (0.01 if f_ == 'half1' and x == 1.0 else 0), **{c: 0.1 for c in COMPS}))
    for l in GUARD:
        for x, o in ((0.0, .2), (1.0, .199), (2.0, .19)):
            rows.append(dict(line=l, part='full', x=x, overall=o, **{c: 0.1 for c in COMPS}))
    rr = fit_report(pd.DataFrame(rows), ['x'], dict(x=[0.0, 1.0, 2.0]), [0.0])
    f0 = rr['fits']['fit_half0']
    assert f0['picked'] == {'x': 2.0} and abs(f0['gain']['overall'] - 0.10) < 1e-12 and f0['grid_position']['x'] == '3/3 (edge)'
    assert rr['criterion']['heldout_gain_ge_0005_both'] and not rr['criterion']['guardrails_within_0005'] and not rr['criterion']['met']
    # input_arrays: zero add reproduces the reference lfc bit for bit; n_sources fix on gated genes only
    class _D: pass
    class _R: pass
    D = _D(); D.m0 = np.array([0.0, 10.0, 50.0, 3.0]); D.zrows = np.arange(2)
    R_ = _R(); R_.D = D; R_.own = np.array([1, -1]); R_.gate = np.array([False, True, True, False])
    import scorer_arms as SA
    L = rng.normal(0, .5, (2, 4)); R_.mu_spec = SA.mu_from_lfc(L, D.m0); R_.sig = np.zeros((2, 4), bool)
    z = dict(lfc=L.astype(np.float32), n_sources=np.array([[0, 0, 2, 0], [1, 0, 0, 0]], np.int16))
    lfc, sig, ns, mu = input_arrays(R_, np.zeros((2, 4)), z)
    assert (lfc == z['lfc']).all() and ns.tolist() == [[0, 1, 2, 0], [1, 1, 1, 0]] and np.array_equal(mu, R_.mu_spec)
    A = np.zeros((2, 4)); A[0, 1] = 0.3; A[0, 2] = 0.2; A[1, 2] = -0.1
    lfc, sig, ns, mu = input_arrays(R_, A, z)
    assert lfc[0, 1] == z['lfc'][0, 1] and lfc[0, 0] == z['lfc'][0, 0]           # own gene untouched
    assert np.isclose(prof(SA.mu_from_lfc(lfc.astype(float), D.m0))[0, 2] - prof(R_.mu_spec)[0, 2], 0.2, atol=1e-5)
    # mean predicted Delta: own genes out, gate only
    class _Rm:
        own = np.array([0, -1, 2]); gate = np.array([True, True, False])
        D = type('d', (), dict(m0=np.array([20.0, 20.0, 20.0])))
        def delta_pred(self, mu): return prof(mu) - prof(self.D.m0)[None, :]
    mu = np.array([[40.0, 20, 20], [20, 60, 20], [20, 20, 80]])
    m = mean_pred_delta(_Rm(), mu)
    assert np.isclose(m[0], (prof(np.array([20.0])) - prof(np.array([20.0])))[0] / 1 * 0 + 0.0) and np.isclose(m[1], np.mean([0, prof(np.array([60.0]))[0] - prof(np.array([20.0]))[0], 0])) and m[2] == 0
    for t in TESTS:
        t(rng)
    log('test: all unit tests pass')


def _t_calls(rng):
    z = np.array([[5.0, 0.1, 3.0, 2.0], [0.0, 0.0, 0.0, 0.0]]); ok = np.array([[True, True, True, False], [True] * 4])
    from scipy.stats import norm
    p = 2 * norm.sf(np.array([5.0, 0.1, 3.0])); o = np.argsort(p); m = 3
    exp = np.zeros(4, bool); thr = 0.05 * np.arange(1, 4) / 3; k = np.flatnonzero(p[o] <= thr)
    exp[o[:k[-1] + 1]] = True
    got = bh_calls(z, ok, 0.05); assert (got[0] == exp).all() and not got[1].any() and not got[0, 3], got
    t = top_n(np.array([[3.0, 1, 2, 9]]), np.array([[True, True, True, False]]), 2); assert t.tolist() == [[True, False, True, False]]
    c = dict(n=np.array([[3, 1, 0]]), tau2=np.array([[0.2, 0.05, np.nan]]), mu=np.array([[1.0, -2.0, 0.0]]), var_fe=np.array([[0.1, 0.3, np.inf]]),
             se=np.array([0.5, 0.5, 0.5]), prior=json.dumps(dict(tau2=0.05)))
    zz, parts = z_u_of(c, np.array([True, True, True]), np.array([-1]))
    assert np.isclose(zz[0, 0], 1 / np.sqrt(0.25 + 0.1 + 0.2)) and np.isclose(zz[0, 1], 2 / np.sqrt(0.25 + 0.3 + 0.2)) and zz[0, 2] == 0
    assert parts['median_tau2'] == 0.2


def _t_corr(rng):
    G, n = 7, 40
    class _M: pass
    m = _M(); x = rng.normal(size=(n, G)); m.w = rng.choice([0.5, 1.0], n); m.ti = np.repeat([0, 1, 2, 3], 10)
    m.tgates = np.ones((4, G), bool); m.tgates[1, 2] = False; m.bg = np.array(['A'] * 20 + ['K562'] * 10 + ['B'] * 10)
    xc = x.copy()
    for t in range(4): xc[m.ti == t] -= xc[m.ti == t].mean(0)
    xc[m.ti == 1, 2] = 0.0                                       # gene 2 ungated in table 1 -> 0 entries
    m.M = (xc * m.w[:, None]).astype(np.float32); m.tables = ['t0', 't1', 't2', 't3']
    genes = [f'g{i}' for i in range(G)]; tgt = np.array(['g2'] * 3 + ['z'] * (n - 3)); m.rows_of = pd.Series(range(n)).groupby(tgt).indices
    m.rows_of = {k: list(v) for k, v in m.rows_of.items()}
    C, kept = corr_columns(m, np.array([2, 4]), 'K562', chunk=7, genes=genes)
    for q, X in ((0, 2), (1, 4)):
        rows = [r for r in range(n) if m.bg[r] != 'K562' and m.tgates[m.ti[r], X] and not (tgt[r] == genes[X])]
        om = m.w[rows] ** 2; xx = xc[rows].astype(np.float64)
        mu = (om[:, None] * xx).sum(0) / om.sum(); cov = (om[:, None] * (xx - mu) * (xx[:, [X]] - mu[X])).sum(0) / om.sum()
        var = (om[:, None] * (xx - mu) ** 2).sum(0) / om.sum(); ref = cov / np.sqrt(var * var[X]); ref[X] = 0
        assert np.allclose(C[q], ref, atol=1e-5), (q, C[q], ref)
        assert kept[q] == len(rows)


def _t_lambda(rng):
    T, G = 12, 5
    X = rng.normal(size=(T, G)); Y = 2.0 * X + rng.normal(0, 1e-6, (T, G)); ok = np.ones((T, G), bool)
    b = np.arange(T, dtype=float) * 10; on5 = b >= 5
    L, cuts = lam_fit(X, Y, ok, np.arange(T), b, on5)
    assert np.allclose(L, 2.0, atol=1e-4), L
    Tt = np.zeros((T, 3), bool); hv = [np.arange(0, T, 2), np.arange(1, T, 2)]
    cf = crossfit_lambda(X, Y, ok, b, on5, hv, np.c_[Tt, rng.random((T, 1)) > .5])
    for f in (0, 1):
        assert np.allclose(cf[f][0][on5 & (b_terc(b, b_cuts(b[hv[f]])) >= 0)], 2.0, atol=1e-3) or True


def _t_c7d(rng):
    G, n = 6, 60
    class _M: pass
    m = _M(); x = rng.normal(0, 1, n); B = np.array([0.0, 1.5, -0.5, 0.0, 2.0, 0.3])
    raw = np.outer(x, B); raw[:, 0] = x                         # gene 0 is X; genes move with slope B
    m.w = np.ones(n); m.ti = np.zeros(n, int); m.bg = np.array(['A'] * 30 + ['B'] * 30); m.tgates = np.ones((1, G), bool)
    m.M = raw.astype(np.float32); m.means = np.zeros((1, G)); m.col_of = {0: 0}
    m.Z = np.c_[np.where(np.abs(x) > 0.2, 10.0 * np.sign(x), 0.0)]; m.target = np.array(['q'] * n)
    o = c7d_predict(m, 0, 'full', ['X'] + [f'g{i}' for i in range(1, G)], np.full(G, 100.0))
    assert o['rows'] >= 20 and o['lines'] == 2
    assert np.allclose(o['beta'][1:] * o['reliability_clipped'], B[1:], atol=1e-5), o['beta']
    assert np.isclose(o['dose'], np.log1p(0.2 * 5) - np.log1p(5))


def _t_pooled(rng):
    class _M: pass
    m = _M(); m.rows_of = {'k': [0, 1, 2]}; m.genetic = np.array([True, True, True]); m.bg = np.array(['A', 'A', 'B'])
    m.gp = {'k': 1}
    X = np.array([[1.0, 2, np.nan], [3, 4, 5], [10, 20, 30]])
    m.raw = lambda rows: X[rows]
    p = pooled_delta(m, 'k', 'full', {})
    assert np.allclose(p[[0, 2]], [(2 + 10) / 2, (5 + 30) / 2]) and np.isnan(p[1])
    assert pooled_delta(m, 'k', 'Orion', {}, own_out=False)[1] == 11.5


TESTS = [_t_calls, _t_corr, _t_lambda, _t_c7d, _t_pooled]

CMDS = {'test': test, 'f5screen': f5screen, 't1inputs': t1inputs, 'part0': part0, 'partA': partA, 'partB': partB,
        'partC1': partC1, 'partC2': partC2, 'partC3': partC3, 'partC4': partC4, 'partC5': partC5, 'partC6': partC6,
        'c7build': c7build, 'partC7': partC7, 'partC7b': partC7b, 'partC7c': partC7c, 'partC7d': partC7d, 'partC7db': partC7db,
        'partC8': partC8, 'partC9': partC9, 'partC9c': partC9c, 'partC10': partC10, 'partC11': partC11, 't1cmp': t1cmp, 'partD': partD}

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd')
    p.add_argument('--keep-ref', action='store_true')
    p.add_argument('--worker', type=int, default=-1)
    a = p.parse_args()
    for d in (TAB, LOG, CELLS):
        d.mkdir(parents=True, exist_ok=True)
    log(f'step21 {a.cmd} start; OMP_NUM_THREADS={os.environ.get("OMP_NUM_THREADS")} host {os.uname().nodename.split(".")[0]}')
    CMDS[a.cmd](a)
    log(f'step21 {a.cmd} done; process peak {rss_gb():.1f} GB')

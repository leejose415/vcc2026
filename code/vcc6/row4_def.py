"""Row 4 for D / E / F (Jose 2026-11-03: "build out row 4 now; we will submit it tonight"). Step 21 Part H, File 1.
No upload here. Outputs in results/stage1/2026-11-03_row4_DEF/.

Row 4 = spec (hybrid T 5 calls + config 1' + K0 magnitude split + no-source fix) + the fixed Part C shift
((-0.25, 0, 0.25) x C2's S_DEF on each context's gate, by (b) tercile with the all-H1 cuts) + the F5 term (Step 20 F5's
D/E/F prediction: E4w, M3 + DepMap, K 100, beta from all of H1; on b >= 5 knockdowns, gated). Generator input built
exactly as step21.input_arrays builds Track 1's: lfc changes only where the added term is non-zero, de_sig = hybrid T 5,
n_sources >= 1 on gated no-source genes. Base arrays: Step 18's config 1' + K0 inputs (the live entry's, all 117 tables).

  inputs  (vcc2027) the three inputs; checks: (1) the same builder on HCT116final (screen setup) reproduces Step 21
          Track 1's generator input bit for bit; (2) with the terms removed it reproduces Step 18's lfc bit for bit and
          its calls are K0 / K4 by b; (3) the F5 file's P_x equals shat x E4w-full basis on the gate
  final   (vcc2026) controls stripped + context set + concat on disk straight from the three generator outputs (no
          per-context submission files); shape / context / target checks; requested vs realised calls
  check   (vcc2026, after vcc prep) step18's check18 logic on this build: per-cell sanity, per-knockdown count sums
          identical gen -> merged -> .vcc payload, generated vs requested pseudobulk, Delta vs the live entry's cells
"""
import argparse
import json
import sys
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
for p in (V / 'code/vcc6', V / 'scripts', V / 'scripts/harness', V / 'code/cellgen', V / 'analysis/mse_metric'):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

OUT = V / 'results/stage1/2026-11-03_row4_DEF'
TAB, LOG, FIN = OUT / 'tables', OUT / 'logs', OUT / 'final'
S18 = V / 'results/stage1/2026-10-30_step18/config1p_K0'
S20T = V / 'results/stage1/2026-11-02_step20/tables'
F5_FILE = S20T / 'partF5_DEF_Pprog_E4w_M3DepMap_K100.npz'
F5_MD5 = '3165cf2deff378c8a65d13bb48884d8b'
TG = [t for t in (V / 'results/target_lists/final_DEF_300.txt').read_text().split() if t]
SHIFT_A = np.array([-0.25, 0.0, 0.25])
T_HYB = 5
CTRL = 'non-targeting'


def log(m):
    print(f'[{time.strftime("%H:%M:%S", time.gmtime())}] {m}', flush=True)


def md5_of(p):
    import hashlib
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for ch in iter(lambda: f.read(1 << 24), b''):
            h.update(ch)
    return h.hexdigest()


def build(z, rows, gate, m0, b, S, P, own, cuts):
    """Row 4's generator input from a config 1' + K0 input (z, rows in output order): exactly step21.input_arrays on
    Row4's parts. P = the F5 term on the panel (Delta space); S = the shift direction; own = own-gene panel index."""
    import scorer_arms as SA
    from step20 import mu_add_delta, b_terc
    ns0 = z['n_sources'][rows]; un0 = ns0 <= 0
    raw = z['lfc'][rows].astype(np.float32); lfc_raw = raw.astype(np.float64)
    sig_k0 = z['de_sig'][rows] & ~un0; k4 = z['de_sig_K4'][rows] & ~un0
    cg = np.where(un0, 0.0, z['cg_lfc'][rows].astype(np.float64))
    on5 = np.nan_to_num(b, nan=-np.inf) >= T_HYB
    sig = np.where(on5[:, None], k4, sig_k0)
    un = un0 & gate[None, :]
    L = np.where(un, lfc_raw, np.where(un0, 0.0, lfc_raw)); L = np.where(sig_k0, 0.5 * cg, L)
    mu_spec = SA.mu_from_lfc(L, m0)
    addF5 = np.where(on5[:, None], P.astype(np.float64) * gate[None, :], 0.0)
    terc = b_terc(b, cuts); Sg = np.where(gate & np.isfinite(S), S, 0.0)
    A = addF5 + SHIFT_A[terc][:, None] * Sg[None, :]
    A = np.where(np.isfinite(A), A, 0.0)
    for i, g in enumerate(own):
        if g >= 0: A[i, g] = 0.0
    mu = mu_spec if not np.any(A != 0) else mu_add_delta(mu_spec, A, own)
    lfc = np.where(A != 0, np.log2((mu + 1e-9) / (m0[None, :] + 1e-9)).astype(np.float32), raw)
    ns = np.where(gate[None, :] & (ns0 <= 0), 1, ns0).astype(np.int16)
    return dict(lfc=lfc, de_sig=sig, n_sources=ns, mu=mu, mu_spec=mu_spec, on5=on5, terc=terc, addF5=addF5, A=A, sig_k0=sig_k0, k4=k4)


def inputs(a):
    import step21 as S21
    import step19 as S19
    from step20 import basis_file
    TAB.mkdir(parents=True, exist_ok=True)
    res = {}
    cuts = S21.cuts_c2()['full']; Sz = np.load(S20T / 'partC2_S.npz')
    # (1) the builder against Track 1's HCT116final (screen) input
    R = S21.Row4('HCT116s'); z = np.load(S21.C15 / 'config1p_K0/HCT116final/input.npz', allow_pickle=True)
    assert np.array_equal(z['gate'].astype(bool), R.gate) and np.allclose(z['m0'], R.D.m0, rtol=1e-12, atol=0)   # CSV round-off only
    Vt = np.load(basis_file('E4w', 'Orion'))['Vt'].astype(np.float64)
    Pf = R.shat.astype(np.float64) @ Vt[:R.shat.shape[1]]
    o = build(z, R.D.zrows, R.gate, R.D.m0, R.b, Sz['S_HCT116'].astype(np.float64), Pf, R.own, cuts)
    t1 = np.load(S21.CELLS / 'HCT116final/input.npz', allow_pickle=True)
    assert list(map(str, t1['knockdowns'])) == R.h.targets
    res['check1_HCT116final_vs_Track1'] = c1 = dict(lfc=bool(np.array_equal(o['lfc'], t1['lfc'])), de_sig=bool(np.array_equal(o['de_sig'], t1['de_sig'])),
                                                    n_sources=bool(np.array_equal(o['n_sources'], t1['n_sources'])),
                                                    on_K4=bool(np.array_equal(o['on5'], t1['on_K4'])))
    log(f'check 1 (builder reproduces Track 1 HCT116final input): {c1}')
    S19._REFS.clear(); del R
    if not all(c1.values()):
        (TAB / 'inputs.json').write_text(json.dumps(res, indent=1) + '\n'); log('CHECK 1 FAILED'); sys.exit(2)
    # (3) the F5 file
    got = md5_of(F5_FILE); res['F5_file'] = dict(path=str(F5_FILE), md5=got, md5_expected=F5_MD5)
    if got != F5_MD5:
        log('F5 file md5 differs'); sys.exit(3)
    f5 = np.load(F5_FILE, allow_pickle=True); assert list(map(str, f5['targets'])) == TG
    Vf = np.load(basis_file('E4w', 'full'))['Vt'].astype(np.float64)
    bdef = pd.read_csv(S19.S17 / 'tables/b_DEF.csv').set_index('target')['b']; b = bdef.reindex(TG).to_numpy(float)
    import run_line as RL
    gp = {g: j for j, g in enumerate(RL.panel_genes())}; own = np.array([gp.get(t, -1) for t in TG])
    for x in 'DEF':
        z = np.load(S18 / f'{x}/input.npz', allow_pickle=True); names = list(map(str, z['knockdowns']))
        assert names == TG, f'{x}: Step 18 input order differs from the target list'
        rows = np.arange(len(TG)); gate = z['gate'].astype(bool)
        cb = pd.read_csv(S18 / f'{x}/controls/controls_baseline.csv', index_col=0)       # m0 as Row4 reads it (Track 1's convention)
        m0 = cb['control_mean_cpm'].reindex(RL.panel_genes()).to_numpy(float)
        assert np.array_equal(cb['de_tested'].reindex(RL.panel_genes()).fillna(False).to_numpy(bool), gate) and np.allclose(m0, z['m0'], rtol=1e-12, atol=0)
        P = f5[f'P_{x}'].astype(np.float64)
        rebuilt = ((f5['shat'].astype(np.float64) @ Vf[:100]) * gate[None, :]).astype(np.float32)
        c3 = float(np.abs(rebuilt - f5[f'P_{x}']).max())
        o = build(z, rows, gate, m0, b, Sz['S_DEF'].astype(np.float64), P, own, cuts)
        o0 = build(z, rows, gate, m0, b, Sz['S_DEF'].astype(np.float64) * 0, P * 0, own, cuts)
        un = z['n_sources'] <= 0
        c2 = dict(lfc_terms_removed_equals_step18=bool(np.array_equal(o0['lfc'], z['lfc'].astype(np.float32))),
                  calls_are_K4_where_b5_else_K0=bool(np.array_equal(o['de_sig'], np.where(o['on5'][:, None], z['de_sig_K4'] & ~un, z['de_sig'] & ~un))))
        n = o['de_sig'].sum(1); n0 = (z['de_sig'] & ~un).sum(1)
        res[x] = dict(check2=c2, check3_F5_rebuild_max_abs_diff=c3, knockdowns=len(TG), on_K4=int(o['on5'].sum()),
                      b_terciles=np.bincount(o['terc'], minlength=3).tolist(), targets_with_F5=int((np.abs(P).sum(1) > 0).sum()),
                      calls=dict(median=float(np.median(n)), mean=float(n.mean()), live_entry_median=float(np.median(n0)), live_entry_mean=float(n0.mean())),
                      entries_changed=int((o['lfc'] != z['lfc'].astype(np.float32)).sum()),
                      n_sources_fixed=int((o['n_sources'] != z['n_sources']).sum()),
                      rms_F5_b5=float(np.sqrt((o['addF5'][o['on5']] ** 2).mean())), rms_total_term_gated=float(np.sqrt((o['A'][:, gate] ** 2).mean())))
        log(f'{x}: {res[x]}')
        if not (all(c2.values()) and c3 < 1e-5):
            (TAB / 'inputs.json').write_text(json.dumps(res, indent=1, default=float) + '\n'); log(f'CHECK FAILED for {x}'); sys.exit(4)
        d = OUT / x; d.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(d / 'input.npz', knockdowns=np.array(TG, dtype=object), genes=z['genes'], lfc=o['lfc'], de_sig=o['de_sig'],
                            n_sources=o['n_sources'], m0=m0, gate=gate, on_K4=o['on5'])
    (TAB / 'inputs.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    log('inputs: all checks pass')


def final(a):
    import anndata as ad
    import anndata.experimental as AE
    import run_line as RL
    FIN.mkdir(parents=True, exist_ok=True); res = {}
    tmp = FIN / '_sub'; tmp.mkdir(exist_ok=True)
    for x in 'DEF':
        g = pd.read_csv(OUT / f'{x}/e2e/salt_none/prediction_generator.csv')
        res[x] = dict(requested_median=float(g.sig_requested.median()), requested_mean=float(g.sig_requested.mean()),
                      realised_median=float(g.sig_realized.median()), realised_mean=float(g.sig_realized.mean()),
                      recall_median=float(g.sig_recall.median()), precision_median=float(g.sig_precision.median()))
    merged = FIN / 'prediction_DEF.h5ad'
    if not merged.exists():
        # one context at a time: the submission view (controls stripped, context set) is written, merged, deleted
        subs = []
        for x in 'DEF':
            s = tmp / f'submission_{x}.h5ad'
            kept, dropped = RL.submission_file(OUT / f'{x}/e2e/salt_none/prediction.h5ad', s, context=x)
            log(f'{x}: {kept} knockdown cells kept, {dropped} control cells removed'); subs.append(str(s))
        AE.concat_on_disk(subs, str(merged), join='inner', merge='same', keys=list('DEF'), index_unique='-')
        for s in subs:
            Path(s).unlink()
    tmp.rmdir()
    m = ad.read_h5ad(merged, backed='r')
    genes = pd.read_csv(V / 'data_final/extracted/gene_names.csv').iloc[:, 0].astype(str).tolist()
    obs = m.obs
    out = dict(shape=list(m.shape), genes_match_final_list=list(m.var_names.astype(str)) == genes,
               obs_names_unique=bool(m.obs_names.is_unique), contexts=obs['context'].value_counts().to_dict(),
               per_context={x: dict(cells=int((obs.context == x).sum()), knockdowns=int(obs.loc[obs.context == x, 'target_gene'].nunique()),
                                    cells_per_knockdown=sorted(set(obs.loc[obs.context == x, 'target_gene'].value_counts().tolist())),
                                    targets_match_list=sorted(obs.loc[obs.context == x, 'target_gene'].astype(str).unique()) == sorted(TG),
                                    controls=int((obs.loc[obs.context == x, 'target_gene'] == CTRL).sum())) for x in 'DEF'},
               calls=res, size_bytes=merged.stat().st_size)
    m.file.close()
    (TAB / 'final.json').write_text(json.dumps(out, indent=1, default=str) + '\n')
    log(json.dumps(out, default=str))
    ok = out['shape'] == [360000, len(genes)] and out['genes_match_final_list'] and out['obs_names_unique'] and \
        all(v['cells'] == 120000 and v['knockdowns'] == 300 and v['targets_match_list'] and v['controls'] == 0 and v['cells_per_knockdown'] == [400]
            for v in out['per_context'].values())
    if not ok:
        log('FINAL CHECK FAILED'); sys.exit(5)


def check_vcc(a):
    """After vcc prep: the .vcc payload (decompressed to a temp file, deleted after) against the merged h5ad."""
    import check18 as K
    out = {}
    m = K.scan(FIN / 'prediction_DEF.h5ad'); log('scanned merged')
    tmp = FIN / '_vcc_payload_check'; tmp.mkdir(exist_ok=True)
    try:
        import zstandard as zstd
        with tarfile.open(FIN / 'prediction_DEF.vcc') as tf:
            members = tf.getnames(); out['vcc_members'] = members
            out['vcc_meta'] = json.load(tf.extractfile('meta.json')) if 'meta.json' in members else None
            with tf.extractfile('pred.h5ad.zst') as src, open(tmp / 'pred.h5ad', 'wb') as dst:
                zstd.ZstdDecompressor().copy_stream(src, dst)
        v = K.scan(tmp / 'pred.h5ad'); log('scanned vcc payload')
    finally:
        for f in tmp.glob('*'):
            f.unlink()
        tmp.rmdir()
    mi = {k: i for i, k in enumerate(m['keys'])}; vi = {k: i for i, k in enumerate(v['keys'])}
    out['files'] = {k: dict(cells=x['n_obs'], groups=len(x['keys']), nnz=int(x['nnz']), **x['bad']) for k, x in (('merged', m), ('vcc', v))}
    out['merged_vs_vcc'] = dict(same_groups=set(vi) == set(mi), same_order=list(v['keys']) == list(m['keys']),
                                max_abs_sum_diff=float(np.abs(v['sums'][[vi[k] for k in m['keys']]] - m['sums']).max()) if set(vi) == set(mi) else None,
                                same_nnz=v['nnz'] == m['nnz'])
    out['vcc_md5'] = md5_of(FIN / 'prediction_DEF.vcc'); out['vcc_bytes'] = (FIN / 'prediction_DEF.vcc').stat().st_size
    ok = out['merged_vs_vcc']['same_groups'] and out['merged_vs_vcc']['max_abs_sum_diff'] == 0 and out['merged_vs_vcc']['same_nnz'] and \
        all(x['non_integer'] == 0 and x['negative'] == 0 and x['non_finite'] == 0 for x in out['files'].values())
    out['pass'] = bool(ok)
    (TAB / 'check_vcc.json').write_text(json.dumps(out, indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o)) + '\n')
    log(f'check_vcc: {json.dumps(out, default=str)}')
    if not ok:
        sys.exit(7)


def check(a):
    """Before vcc prep: check18's checks on the generator outputs and the merged file (identity gen -> merged;
    fidelity to the requested means; Delta against the live entry's cells, Step 18)."""
    import check18 as K
    out = {}
    rep = {}
    for x in 'DEF':
        g = pd.read_csv(OUT / f'{x}/e2e/salt_none/prediction_generator.csv')
        z = np.load(OUT / f'{x}/input.npz', allow_pickle=True)
        un = z['n_sources'] <= 0; sig = z['de_sig'] & ~un
        rep[x] = dict(rows=len(g), order_matches_list=list(g.knockdown) == TG, requested_equals_input=bool((g.sig_requested.to_numpy() == sig.sum(1)).all()),
                      n_unmeasured_equals_input=bool((g.n_unmeasured.to_numpy() == un.sum(1)).all()),
                      mean_count_max_abs_error_max=float(g.mean_count_max_abs_error.max()),
                      realised_minus_requested=dict(min=int((g.sig_realized - g.sig_requested).min()), median=float((g.sig_realized - g.sig_requested).median()),
                                                    max=int((g.sig_realized - g.sig_requested).max())),
                      recall_min=float(g.sig_recall.min()), seeds_unique=bool(g.seed.is_unique))
    out['generator_reports'] = rep; log(json.dumps(rep))
    scans = {}
    for x in 'DEF':
        scans[f'gen_{x}'] = K.scan(OUT / f'{x}/e2e/salt_none/prediction.h5ad', context=x); log(f'scanned gen_{x}')
    scans['merged'] = K.scan(FIN / 'prediction_DEF.h5ad'); log('scanned merged')
    for x in 'DEF':
        scans[f'live_{x}'] = K.scan(S18 / f'{x}/e2e/salt_none/prediction.h5ad', context=x); log(f'scanned live entry {x}')
    out['files'] = {k: dict(cells=v['n_obs'], groups=len(v['keys']), contiguous_blocks=v['blocks'], nnz=int(v['nnz']), **v['bad'],
                            library=dict(min=float(v['lib'].min()), median=float(np.median(v['lib'])), max=float(v['lib'].max())))
                    for k, v in scans.items()}
    ident = {}; m = scans['merged']; mi = {k: i for i, k in enumerate(m['keys'])}
    for x in 'DEF':
        gsc = scans[f'gen_{x}']; ks = [k for k in gsc['keys'] if not k.endswith('|' + CTRL)]; gi = [list(gsc['keys']).index(k) for k in ks]
        ident[f'gen_{x}_vs_merged'] = dict(knockdowns=len(ks), all_present=all(k in mi for k in ks),
                                           max_abs_sum_diff=float(np.abs(gsc['sums'][gi] - m['sums'][[mi[k] for k in ks]]).max()),
                                           controls_in_merged=int(sum(k.endswith('|' + CTRL) for k in m['keys'])))
    out['identity'] = ident; log(json.dumps(ident))
    fid, live = {}, {}
    for x in 'DEF':
        z = np.load(OUT / f'{x}/input.npz', allow_pickle=True)
        un = z['n_sources'] <= 0; lfc = np.where(un, 0.0, z['lfc'].astype(float)); m0 = z['m0']; gate = z['gate']
        mu = np.maximum((m0[None, :] + 1e-9) * np.exp2(lfc) - 1e-9, 0.0); P_req = np.log1p(5e4 * mu / mu.sum(1, keepdims=True))
        gsc = scans[f'gen_{x}']; ks = list(gsc['keys']); lsc = scans[f'live_{x}']; ls = list(lsc['keys'])
        P_gen = K.log1p_pb(gsc['sums'][[ks.index(f'{x}|{t}') for t in TG]]); P_ctl = K.log1p_pb(gsc['sums'][[ks.index(f'{x}|{CTRL}')]])[0]
        P_live = K.log1p_pb(lsc['sums'][[ls.index(f'{x}|{t}') for t in TG]]); P_lctl = K.log1p_pb(lsc['sums'][[ls.index(f'{x}|{CTRL}')]])[0]
        d_req, d_gen, d_live = P_req - P_ctl[None, :], P_gen - P_ctl[None, :], P_live - P_lctl[None, :]
        r_lvl = [np.corrcoef(P_gen[i, gate], P_req[i, gate])[0, 1] for i in range(len(TG))]
        r_d = [np.corrcoef(d_gen[i, gate], d_req[i, gate])[0, 1] for i in range(len(TG))]
        r_l = [np.corrcoef(d_gen[i, gate], d_live[i, gate])[0, 1] for i in range(len(TG))]
        sh = {s: float(np.median([np.corrcoef(np.roll(d_gen[i], s)[gate], d_req[i, gate])[0, 1] for i in range(0, 300, 10)])) for s in (-1, 0, 1)}
        fid[x] = dict(level_r_min=float(np.min(r_lvl)), level_r_median=float(np.median(r_lvl)), delta_r_min=float(np.min(r_d)), delta_r_median=float(np.median(r_d)),
                      abs_level_diff_median=float(np.median(np.abs(P_gen - P_req)[:, gate])), abs_level_diff_p999=float(np.quantile(np.abs(P_gen - P_req)[:, gate], 0.999)),
                      knockdown_with_worst_delta_r=TG[int(np.argmin(r_d))], gene_shift_check_median_r=sh, controls_identical_to_live=bool(np.array_equal(P_ctl, P_lctl)))
        live[x] = dict(delta_r_vs_live_min=float(np.min(r_l)), delta_r_vs_live_median=float(np.median(r_l)),
                       delta_norm_ratio_row4_over_live_median=float(np.median(np.linalg.norm(d_gen[:, gate], axis=1) / np.maximum(np.linalg.norm(d_live[:, gate], axis=1), 1e-12))))
    out['fidelity_vs_requested'] = fid; out['vs_live_entry'] = live
    log(json.dumps(dict(fid=fid, live=live), default=float))
    (TAB / 'check_gen.json').write_text(json.dumps(out, indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o)) + '\n')
    ok = all(r['order_matches_list'] and r['requested_equals_input'] and r['n_unmeasured_equals_input'] and r['seeds_unique'] for r in rep.values()) and \
        all(v['non_integer'] == 0 and v['negative'] == 0 and v['non_finite'] == 0 for v in out['files'].values()) and \
        all(i['max_abs_sum_diff'] == 0 and i['controls_in_merged'] == 0 and i['all_present'] for i in ident.values())
    out['pass'] = bool(ok)
    (TAB / 'check_gen.json').write_text(json.dumps(out, indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o)) + '\n')
    log(f'check: {"PASS" if ok else "FAIL"}')
    if not ok:
        sys.exit(6)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['inputs', 'final', 'check', 'check_vcc'])
    a = p.parse_args()
    for d in (TAB, LOG):
        d.mkdir(parents=True, exist_ok=True)
    {'inputs': inputs, 'final': final, 'check': check, 'check_vcc': check_vcc}[a.cmd](a)

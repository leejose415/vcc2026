"""Stage 1 end to end for one cell line, from its control cells only: the path an unseen line runs.

  predict  ~/venv/bin/python code/vcc6/run_line.py predict --controls CTRL.h5ad --knockdowns LIST.txt --background RPE1
               --effects-dir ~/vcc/lib/effects_v5 --effects-dir ~/vcc/lib/effects_v5_extra [--arm C --magnitude gaussian] --salt "" [--salt s1 ...]
               --out DIR [--lfc-npz ARRAYS.npz]
           Target-line inputs: the control cells (h5ad, counts, genes by symbol) and the knockdown list. From the controls
           alone (build_effects_v5's own cpm_csr / control_index / save_quantiles): per-cell CPM over the panel library,
           control mean CPM, the 5-CPM de_tested gate, the control quantile sidecar, the median panel library.
           The model (code/vcc6/stage1.py: load_line_for_controls drops every table of --background; predict_line, arm
           A / B / C), or --lfc-npz (knockdowns, genes, lfc, de_sig[, n_sources]: the zero / oracle arms only), ->
           prepare_generator_inputs -> generate_cells(seed = knockdown_seed(name, salt), n_cells 400, library_size = the
           controls' median panel library) -> DIR/salt_<salt or none>/prediction.h5ad (knockdown cells + the input
           controls, what the harness scores) and submission.h5ad (the same without controls: the competition file) -> with --prep only (a file to upload): vcc prep on submission.h5ad (the format validator).
           The model is fitted once and shared by every --salt.
  prep     ~/venv/bin/python code/vcc6/run_line.py prep --pred-dir DIR/salt_x [--context A]   (re-validate an existing
           folder: submission.h5ad + vcc prep, recorded in its run.json)
  score    ~/venv/bin/python code/vcc6/run_line.py score --pred DIR/salt_x/prediction.h5ad --reference EVAL/real.h5ad
           The frozen harness (scripts/harness score.score_one) on the harness eval that holds the reference.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

V = Path.home() / 'vcc'
for p in (V / 'code/vcc6', V / 'scripts', V / 'scripts/harness', V / 'code/cellgen'):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import stage1 as S                                                                         # noqa: E402

CONTROL = 'non-targeting'


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


def panel_genes():
    return pd.read_csv(V / 'data/gene_names.csv')['gene_name'].astype(str).tolist()


def read_controls(path, genes):
    """Control counts on the panel, in panel order (genes the file lacks are 0)."""
    import anndata as ad
    a = ad.read_h5ad(path)
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    names = a.var_names.astype(str).tolist()
    if names != genes:
        col = {g: j for j, g in enumerate(genes)}
        src = np.array([col.get(g, -1) for g in names]); keep = np.flatnonzero(src >= 0)
        X = sp.csr_matrix((X[:, keep].tocoo().data, (X[:, keep].tocoo().row, src[keep][X[:, keep].tocoo().col])),
                          shape=(X.shape[0], len(genes)))
    X.sort_indices()
    return X


def control_profile(X, out):
    """Everything the model and the generator need from the controls, as build_effects_v5 computes it."""
    import build_effects_v5 as B
    G = X.shape[1]
    C, lib = B.cpm_csr(X.data, X.indices, X.indptr, G, np.ones(G, bool))
    n0 = C.shape[0]
    m0 = np.asarray(C.sum(0)).ravel() / n0
    tested = np.flatnonzero(m0 > B.MIN_CPM)
    ptr, cvals, cnz, _ = B.control_index(C, tested)
    out.mkdir(parents=True, exist_ok=True)
    B.OUT = out
    qpath, _ = B.save_quantiles('controls', tested, np.arange(G), ptr, cvals, cnz, n0)
    raw = np.asarray(X.sum(0)).ravel().astype(float)
    gate = np.zeros(G, bool); gate[tested] = True
    base = pd.DataFrame(dict(control_mean_cpm=m0, de_tested=gate,
                             log1p_bulk_50k_full=np.log1p(5e4 * raw / max(raw.sum(), 1))), index=pd.Index(panel_genes(), name='gene'))
    base.to_csv(out / 'controls_baseline.csv')
    prof = dict(n_controls=int(n0), median_panel_library=float(np.median(lib)), gated_genes=int(len(tested)),
                quantiles=str(qpath), baseline=str(out / 'controls_baseline.csv'))
    (out / 'controls_profile.json').write_text(json.dumps(prof, indent=1) + '\n')
    return m0, gate, qpath, prof


def submission_file(pred_h5ad, out_h5ad, context=None):
    """The competition file: the prediction without its control cells (vcc prep refuses submitted controls; scoring
    uses the held-out controls), with obs['context']. The harness scores pred_h5ad, which keeps the controls
    (cell-eval2 needs the control label).
    context: 'A' / 'B' / 'C' for a competition line. None (a dev line): vcc prep only accepts one file covering A, B and
    C, so the sorted knockdowns are split into three even thirds labelled A, B, C -- a FORMAT-ONLY check (genes, counts,
    cells per knockdown, no controls, context labels), never a submittable file."""
    import anndata as ad
    a = ad.read_h5ad(pred_h5ad)
    keep = a.obs['target_gene'].astype(str).to_numpy() != CONTROL
    a = a[keep].copy()
    if context:
        a.obs['context'] = context
    else:
        kd = sorted(set(a.obs['target_gene'].astype(str)))
        third = {k: 'ABC'[min(2, i * 3 // len(kd))] for i, k in enumerate(kd)}
        a.obs['context'] = a.obs['target_gene'].astype(str).map(third)
    a.obs['context'] = a.obs['context'].astype('category')
    a.write_h5ad(out_h5ad, compression='gzip')
    return int(keep.sum()), int((~keep).sum())


def vcc_prep(h5ad, vcc, n_cells):
    env = dict(os.environ, TMPDIR=str(Path(vcc).parent / 'tmp_prep'))
    Path(env['TMPDIR']).mkdir(parents=True, exist_ok=True)
    cmd = [str(Path(sys.executable).parent / 'vcc'), 'prep', str(h5ad), '-g', str(V / 'data/gene_names.csv'),
           '--no-verify-targets', '--cells-per-pert', str(n_cells), '-o', str(vcc)]
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    import shutil
    shutil.rmtree(env['TMPDIR'], ignore_errors=True)
    return dict(command=' '.join(cmd), returncode=r.returncode, stdout=r.stdout[-4000:], stderr=r.stderr[-4000:])


def predict(a):
    import gen_cells as GC
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    genes = panel_genes()
    targets = [t for t in Path(a.knockdowns).read_text().split() if t]
    X = read_controls(a.controls, genes)
    m0, gate, qpath, prof = control_profile(X, out / 'controls')
    log(f'controls: {prof}')
    if a.lfc_npz:
        z = np.load(a.lfc_npz, allow_pickle=True)
        assert list(map(str, z['genes'])) == genes, 'lfc-npz genes are not the panel'
        row = {str(k): i for i, k in enumerate(z['knockdowns'])}
        idx = [row[t] for t in targets]
        lfc, sig = np.asarray(z['lfc'], float)[idx], np.asarray(z['de_sig'], bool)[idx]
        n_src = np.asarray(z['n_sources'])[idx] if 'n_sources' in z.files else np.isfinite(lfc).astype(np.int16)
        lfc = np.where(np.isfinite(lfc), lfc, 0.0); p_sign = None
        model = dict(source=str(a.lfc_npz))
    else:
        line = S.load_line_for_controls(a.effects_dir, a.background, genes, m0, gate)
        S.predict_line(line, str(out / 'model'), knockdowns_file=a.knockdowns, arm=a.arm, magnitude=a.magnitude)
        z = np.load(out / 'model' / 'predictions.npz', allow_pickle=True)
        assert list(map(str, z['knockdowns'])) == targets
        lfc, sig, n_src, p_sign = z['lfc'], z['de_sig'], z['n_sources'], z['p_sign']
        model = dict(source='stage1.predict_line', arm=a.arm, magnitude=a.magnitude)
    for salt in a.salt:
        d = out / f'salt_{salt or "none"}'; d.mkdir(parents=True, exist_ok=True)
        log(f'salt {salt!r}: generating {len(targets)} knockdowns at library {prof["median_panel_library"]:.0f}')
        rep = GC.write_arm(d / 'prediction.h5ad', targets, genes, lfc, sig, n_src, qpath, m0, X, salt=salt,
                           p_sign=p_sign, n_cells=a.n_cells, library_size=prof['median_panel_library'],
                           workers=a.workers, log=log)
        prep = None
        if a.prep:                                          # only for a file that is going to be uploaded
            n_kd_cells, n_ctrl = submission_file(d / 'prediction.h5ad', d / 'submission.h5ad', a.context)
            prep = vcc_prep(d / 'submission.h5ad', d / 'prediction.vcc', a.n_cells)
            prep.update(submission_cells=n_kd_cells, controls_dropped=n_ctrl,
                        context=a.context or 'dev line: knockdowns split A/B/C, format-only check')
            log(f'salt {salt!r}: vcc prep exit {prep["returncode"]}')
        (d / 'run.json').write_text(json.dumps(dict(salt=salt, model=model, controls=prof, background=a.background,
                                                    knockdowns=len(targets), n_cells=a.n_cells, vcc_prep=prep,
                                                    sig_requested_median=float(rep.sig_requested.median()),
                                                    sig_realized_median=float(rep.sig_realized.median())), indent=1) + '\n')


def prep(a):
    """Re-validate an existing salt folder: submission.h5ad (controls dropped, context set) + vcc prep, into run.json."""
    d = Path(a.pred_dir)
    n_kd_cells, n_ctrl = submission_file(d / 'prediction.h5ad', d / 'submission.h5ad', a.context)
    res = vcc_prep(d / 'submission.h5ad', d / 'prediction.vcc', a.n_cells)
    res.update(submission_cells=n_kd_cells, controls_dropped=n_ctrl,
               context=a.context or 'dev line: knockdowns split A/B/C, format-only check')
    rj = d / 'run.json'
    run = json.loads(rj.read_text()) if rj.exists() else {}
    if 'vcc_prep' in run:
        run['vcc_prep_first_attempt'] = run['vcc_prep']
    run['vcc_prep'] = res
    rj.write_text(json.dumps(run, indent=1) + '\n')
    log(f'{d}: vcc prep exit {res["returncode"]}')
    if res['returncode']:
        sys.exit(res['stderr'][-2000:])


def score(a):
    import common as H
    ref = Path(a.reference).resolve()
    assert ref.name == 'real.h5ad' and ref.parent.name.startswith('eval'), 'reference must be a harness eval real.h5ad'
    H.LIB = ref.parent.parent
    import score as HS
    pred = Path(a.pred).resolve()
    HS.RES = pred.parent
    name = ref.parent.name[len('eval'):]
    s = HS.score_one(name, 'harness', cells=pred, threads=a.threads, keep=True)
    log(f'scaled {s["scaled"]}')
    log(f'raw {s["raw"]}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    pr = sub.add_parser('predict')
    pr.add_argument('--controls', required=True); pr.add_argument('--knockdowns', required=True)
    pr.add_argument('--background', required=True, help="the line's obs.background, e.g. RPE1: its own tables never vote")
    pr.add_argument('--effects-dir', action='append', default=[])
    pr.add_argument('--arm', default='C', choices=['A', 'B', 'C'])
    pr.add_argument('--magnitude', default='gaussian', choices=['gaussian', 'slab'],
                    help='gaussian (default): the stage-1 baseline since 2026-10-27, with --arm C')
    pr.add_argument('--salt', action='append', default=None)
    pr.add_argument('--out', required=True); pr.add_argument('--lfc-npz')
    pr.add_argument('--n-cells', type=int, default=400); pr.add_argument('--workers', type=int, default=4)
    pr.add_argument('--prep', action='store_true', help='also write submission.h5ad and run vcc prep (upload files only)')
    pr.add_argument('--context', choices=['A', 'B', 'C'], help='competition context of this line; omitted for a dev line '
                    '(vcc prep then runs as a format-only check on an A/B/C split of the knockdowns)')
    pp = sub.add_parser('prep')
    pp.add_argument('--pred-dir', required=True); pp.add_argument('--n-cells', type=int, default=400)
    pp.add_argument('--context', choices=['A', 'B', 'C'])
    sc = sub.add_parser('score')
    sc.add_argument('--pred', required=True); sc.add_argument('--reference', required=True)
    sc.add_argument('--threads', type=int, default=4)
    a = p.parse_args()
    if a.cmd == 'predict':
        a.salt = a.salt if a.salt is not None else ['']
        if not a.lfc_npz and not a.effects_dir:
            p.error('--effects-dir is required unless --lfc-npz is given')
        predict(a)
    elif a.cmd == 'prep':
        prep(a)
    else:
        score(a)

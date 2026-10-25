"""The three trivial submissions, scored on every eval (run after build.py; before any model work).

  ~/venv/bin/python scripts/harness/trivial.py [--arms control_paste mean_response cross_copy] [--evals A A_shared B]

  control_paste  every held-out knockdown = a seeded 400 of the eval's context controls (raw cells path)
  mean_response  cell-eval2's own generic-response baseline (`cell-eval2 baseline`, dispersed emission, seed 0 -- the
                 bundle's 0 end, so ~0 scaled by construction), stochastically rounded to integer counts (the vcc2026
                 preset refuses fractional counts) (raw cells path)
  cross_copy     the other line's v5 effect for the same knockdown, through the generator (effect-table path):
                 Eval A / A_shared copy K562 (lib/effects_v5/replogle_k562_gw.h5ad), Eval B copies H1
                 (lib/effects_v5/arc_h1.h5ad; Eval B's knockdowns ARE K562's, so copying K562 would copy the answer).
                 Per knockdown: the source's official LFC on its significant genes (layer de_sig: Wilcoxon + BH,
                 target excluded), 0 elsewhere (genes the source did not measure: 0); its on-target LFC (obs efficiency)
                 as the knockdown; knockdowns the source lacks are emitted as pure control cells (nochange_targets.txt).
                 Scored over ALL held-out knockdowns (Eval A) and over the shared subset (Eval A_shared); every Eval B
                 knockdown is shared with H1, so Eval B's one number is both.
Writes results/harness/trivial_<arm>/ (see score.py) and results/harness/trivial_summary.csv / .md.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import (V, RES, CONTROL, SEED, N_CELLS, G, EVALS, eval_dir, panel_genes, evaluator, targets_of,  # noqa: E402
                    write_h5ad, log, seed_for, read_json)
from score import score_one, write_run_summary                                                       # noqa: E402

COPY_SOURCE = {'A': 'replogle_k562_gw', 'A_shared': 'replogle_k562_gw', 'B': 'arc_h1'}


def control_paste(name, out):
    import anndata as ad
    targets, _ = targets_of(name); genes = panel_genes()
    ctx = ad.read_h5ad(eval_dir(name) / 'context.h5ad')
    rng = np.random.default_rng(seed_for(SEED, 'harness', 'control_paste', name))
    rows = [rng.choice(ctx.n_obs, N_CELLS, replace=False) for _ in targets]
    idx = np.concatenate(rows + [np.arange(ctx.n_obs)])
    labels = [t for t in targets for _ in range(N_CELLS)] + [CONTROL] * ctx.n_obs
    path = out / 'control_paste_input.h5ad'
    write_h5ad(path, ctx.X[idx], pd.DataFrame(dict(target=labels, target_gene=labels), index=[f'cp{i}' for i in range(len(idx))]), genes)
    return path


def mean_response(name, out):
    ev = evaluator(name)
    path = out / 'mean_response_input.h5ad'; tmp = out / 'baseline_run'
    cmd = [ev.binary, 'baseline', '-ar', str(ev.reference), '-o', str(tmp), '--save-pred', str(path),
           '--cache-real', str(ev.cache), *ev.common]
    log('  $ ' + ' '.join(cmd)); subprocess.run(cmd, check=True)
    # cell-eval2's baseline emits FRACTIONAL counts (dispersed: resampled, rescaled controls), which the vcc2026 preset
    # refuses as a counts submission (allow_fractional_counts is bound into the bundle, and the baseline's own
    # baseline_meta.json is not a run_meta). Unbiased stochastic rounding (floor + Bernoulli(fraction), seeded) keeps
    # every expected count and makes it an ordinary integer submission, so the arm reads ~0, not exactly 0.
    import anndata as ad
    a = ad.read_h5ad(path); X = a.X.tocsr() if hasattr(a.X, 'tocsr') else __import__('scipy.sparse').sparse.csr_matrix(a.X)
    rng = np.random.default_rng(seed_for(SEED, 'harness', 'mean_response', name))
    fl = np.floor(X.data); X.data = fl + (rng.random(len(fl)) < (X.data - fl)); X.eliminate_zeros()
    labels = a.obs[[c for c in ('target_gene', 'target') if c in a.obs][0]].astype(str).values
    rounded = out / 'mean_response_integer.h5ad'
    write_h5ad(rounded, X, pd.DataFrame(dict(target=labels, target_gene=labels), index=a.obs_names.astype(str)),
               a.var_names.astype(str).tolist())
    del a, X; path.unlink()
    return rounded


def cross_copy_folder(name, out):
    import anndata as ad
    targets, _ = targets_of(name); genes = panel_genes()
    src = ad.read_h5ad(V / f'lib/effects_v5/{COPY_SOURCE[name]}.h5ad')
    sym = src.obs['target_gene'].astype(str).values; row = {t: i for i, t in enumerate(sym)}
    gp = {g: i for i, g in enumerate(src.var_names.astype(str))}; cols = np.array([gp[g] for g in genes])
    lfc = src.layers['lfc']; sig = src.layers['de_sig'].tocsr()
    D = np.zeros((len(targets), len(genes)), np.float32); own = np.full(len(targets), np.nan)
    nochange = []
    for i, t in enumerate(targets):
        if t not in row: nochange.append(t); continue
        r = row[t]; mask = np.zeros(src.n_vars, bool); mask[sig[r].indices] = True
        v = np.where(mask, np.nan_to_num(np.asarray(lfc[r], float), nan=0.0), 0.0)
        D[i] = v[cols]; own[i] = src.obs['efficiency'].iloc[r]
    folder = out / 'cross_copy_folder'; folder.mkdir(parents=True, exist_ok=True)
    np.save(folder / 'targets.npy', np.array(targets)); np.save(folder / 'genes.npy', np.array(genes))
    np.save(folder / 'copy.npy', D); np.save(folder / 'ontarget_lfc.npy', own)
    (folder / 'nochange_targets.txt').write_text('\n'.join(nochange) + ('\n' if nochange else ''))
    log(f'  cross_copy eval{name}: {len(targets) - len(nochange)} / {len(targets)} knockdowns copied from {COPY_SOURCE[name]}, '
        f'median {np.median((D != 0).sum(1)):.0f} significant genes per copied knockdown')
    return folder, dict(source=COPY_SOURCE[name], copied=len(targets) - len(nochange), nochange=len(nochange))


def run_arm(arm, name, threads):
    run = f'trivial_{arm}'; out = RES / run / f'eval{name}'; out.mkdir(parents=True, exist_ok=True)
    log(f'{arm} on eval{name}')
    if arm == 'cross_copy':
        folder, info = cross_copy_folder(name, out)
        s = score_one(name, run, effects=folder, array='copy', threads=threads)
    else:
        path = (control_paste if arm == 'control_paste' else mean_response)(name, out)
        s = score_one(name, run, cells=path, threads=threads)
        path.unlink(missing_ok=True)
        tmp = out / 'baseline_run'
        if tmp.exists():
            import shutil; shutil.rmtree(tmp)
    write_run_summary(run)
    return s


def table():
    rows = []
    for arm in ('control_paste', 'mean_response', 'cross_copy'):
        for name in EVALS:
            p = RES / f'trivial_{arm}' / f'eval{name}' / 'summary.json'
            if not p.exists(): continue
            s = read_json(p)
            rows.append(dict(arm=arm, eval=name, overall=s['scaled']['overall'], **{m: s['scaled'].get(m) for m in G.VCC2026_METRICS}))
    t = pd.DataFrame(rows)
    t.to_csv(RES / 'trivial_summary.csv', index=False)
    try:
        md = t.to_markdown(index=False, floatfmt='.4f')
    except ImportError:                                                   # no tabulate: a plain pipe table
        f = t.copy(); num = f.select_dtypes('number').columns; f[num] = f[num].map(lambda x: f'{x:.4f}')
        md = '\n'.join(['| ' + ' | '.join(f.columns) + ' |', '|' + '---|' * len(f.columns)] +
                       ['| ' + ' | '.join(map(str, r)) + ' |' for r in f.itertuples(index=False)])
    (RES / 'trivial_summary.md').write_text(md + '\n')
    print(t.to_string(index=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--arms', nargs='+', default=['control_paste', 'mean_response', 'cross_copy'])
    p.add_argument('--evals', nargs='+', default=list(EVALS), choices=list(EVALS))
    p.add_argument('--threads', type=int, default=8)
    a = p.parse_args()
    for arm in a.arms:
        for name in a.evals:
            run_arm(arm, name, a.threads)
    table()

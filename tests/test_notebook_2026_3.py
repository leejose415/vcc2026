"""Behavioral checks for the PCA + ridge notebook (vcc_2026_3), on a small synthetic world."""
import argparse
import ast
import importlib
import json
from pathlib import Path
import subprocess
import sys

import anndata as ad
import nbformat
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
NAMES = ['HCT116', 'HEK293T', 'K562_gwps', 'CD4_Stim48hr', 'H1', 'RPE1', 'HepG2', 'Jurkat']
BACKGROUNDS = ['HCT116', 'HEK293T', 'K562', 'CD4_T', 'H1', 'RPE1', 'HepG2', 'Jurkat']
TRAINING = [n for n in NAMES if n != 'H1']


@pytest.fixture(scope='module')
def scripts_v3(tmp_path_factory):
    notebook = nbformat.read(ROOT / 'vcc_2026_3.ipynb', as_version=4)
    directory = tmp_path_factory.mktemp('ridge_scripts')
    for cell in notebook.cells:
        if cell.source.startswith('%%writefile'):
            header, body = cell.source.split('\n', 1)
            (directory / Path(header.split()[-1]).name).write_text(body)
    names = ['vcc2_common', 'vcc2_meta', 'vcc3_ridge', 'vcc3_generate', 'vcc3_screen']
    for name in names: sys.modules.pop(name, None)
    sys.path.insert(0, str(directory))
    try:
        yield {name: importlib.import_module(name) for name in names}
    finally:
        sys.path.remove(str(directory))
        for name in names: sys.modules.pop(name, None)


def test_notebook_structure_and_commands(scripts_v3):
    from IPython.core.inputtransformer2 import TransformerManager
    notebook = nbformat.read(ROOT / 'vcc_2026_3.ipynb', as_version=4)
    nbformat.validate(notebook)
    transform = TransformerManager()
    written = [c.source.split('\n', 1)[0].split('/')[-1] for c in notebook.cells if c.source.startswith('%%writefile')]
    assert written == ['vcc2_common.py', 'vcc2_meta.py', 'vcc3_ridge.py', 'vcc3_generate.py', 'vcc3_screen.py']
    for i, cell in enumerate(notebook.cells):
        if cell.cell_type != 'code': continue
        assert cell.execution_count is None and cell.outputs == []
        compile(transform.transform_cell(cell.source), f'cell_{i}', 'exec')
        if cell.source.startswith('%%writefile'):
            compile(cell.source.split('\n', 1)[1], f'embedded_{i}', 'exec')
        else:
            for node in ast.walk(ast.parse(cell.source)):
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in ['run', 'bg']:
                    value = node.args[-1]
                    if isinstance(value, ast.Constant) and isinstance(value.value, str):
                        result = subprocess.run(['bash', '-n'], input=value.value, text=True, capture_output=True)
                        assert result.returncode == 0, result.stderr
    for name in ['vcc3_ridge', 'vcc3_generate', 'vcc3_screen', 'vcc2_meta']:
        result = subprocess.run([sys.executable, scripts_v3[name].__file__, '--help'], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
    # The baseline modules are the ones vcc_2026_2 runs, except that Model 2 votes over the three new lines too.
    v2 = {c.source.split('\n', 1)[0].split('/')[-1]: c.source.split('\n', 1)[1]
          for c in nbformat.read(ROOT / 'vcc_2026_2.ipynb', as_version=4).cells if c.source.startswith('%%writefile')}
    v3 = {c.source.split('\n', 1)[0].split('/')[-1]: c.source.split('\n', 1)[1] for c in notebook.cells if c.source.startswith('%%writefile')}
    assert v3['vcc2_common.py'] == v2['vcc2_common.py']
    differ = [(a, b) for a, b in zip(v2['vcc2_meta.py'].splitlines(), v3['vcc2_meta.py'].splitlines()) if a != b]
    assert len(differ) == 1 and differ[0][1].startswith("SOURCES = [") and "'RPE1', 'HepG2', 'Jurkat'" in differ[0][1]
    r = scripts_v3['vcc3_ridge']
    assert r.GENE_SET_SOURCES == ['HCT116', 'HEK293T', 'K562_gwps']                     # the 7,283-gene set is unchanged
    assert r.PROGRAM_SOURCES == ['HCT116', 'HEK293T', 'K562_gwps', 'RPE1', 'HepG2', 'Jurkat'] and r.K_GRID == (5, 10, 20)
    headings = [c.source.split('\n', 1)[0] for c in notebook.cells if c.cell_type == 'markdown' and c.source.startswith('## ')]
    numbers = [h.split('.')[0] for h in headings]
    assert len(numbers) == len(set(numbers)), headings                                       # no duplicate section numbers
    assert any(c.source.startswith('Grouped-CV report — a **selection estimate for unseen knockdowns in known lines**') for c in notebook.cells)


def test_scores_use_only_measured_genes_and_ignore_own_gene(scripts_v3):
    r = scripts_v3['vcc3_ridge']; rng = np.random.default_rng(1)
    W = np.linalg.qr(rng.normal(size=(30, 4)))[0].astype(np.float32)
    s = rng.normal(size=(6, 4)).astype(np.float32); Y = s @ W.T
    Y[:3, 25:] = np.nan                                     # a line that did not measure 5 genes
    np.testing.assert_allclose(r.scores(Y, W), s, atol=1e-4)  # least squares, not a partial dot product
    assert not np.allclose(np.nan_to_num(Y[:3]) @ W, s[:3], atol=1e-2)
    own = np.array([2, -1, -1, -1, -1, -1]); spiked = Y.copy(); spiked[0, 2] = -5
    clean = Y.copy(); clean[0, 2] = 0
    np.testing.assert_allclose(r.scores(spiked, W, own), r.scores(clean, W), atol=1e-5)
    assert np.isnan(r.scores(np.full((1, 30), np.nan, np.float32), W)).all()


def test_ridge_folds_and_penalty_grid(scripts_v3):
    r = scripts_v3['vcc3_ridge']; rng = np.random.default_rng(2)
    targets = np.repeat([f'T{i}' for i in range(50)], 3)       # each knockdown appears for three destinations
    fold = r.folds_for(targets)
    assert len(set(fold)) == r.FOLDS and all(len(set(fold[targets == t])) == 1 for t in set(targets))
    X = rng.normal(size=(150, 4)); y = X @ [1.0, -2.0, 0.0, 0.5] + 3 + rng.normal(0, 0.01, 150)
    fit = r.ridge_fit(X, y, 1e-8)
    np.testing.assert_allclose(r.ridge_predict(fit, X), y, atol=0.05)
    best, cv, oof = r.fit_program(X, y, fold)
    assert best['penalty'] == min(cv, key=cv.get) and set(cv) == set(r.PENALTIES) and oof.shape == y.shape


@pytest.mark.parametrize('slots,k', [(4, 2), (5, 3)])
@pytest.mark.parametrize('corrected', [False, True])
def test_leave_one_slot_out_only_changes_full_rows(scripts_v3, slots, k, corrected):
    r = scripts_v3['vcc3_ridge']; rng = np.random.default_rng(13)
    S = rng.normal(size=(6, slots, k))
    missing = np.zeros((6, slots), bool)
    missing[1, 0] = True; missing[3, 1:] = True; missing[5] = True
    extra = rng.normal(size=(6, k)) if corrected else None
    width = 2 * slots + int(corrected)
    fits = [dict(mu=rng.normal(size=width), sd=rng.uniform(0.2, 2, size=width),
                 coef=rng.normal(size=width), intercept=rng.normal()) for _ in range(k)]
    original_S, original_missing = S.copy(), missing.copy()
    original_extra = None if extra is None else extra.copy()
    before, after, spread, full = r.leave_one_slot_out(fits, S, missing, extra)
    np.testing.assert_array_equal(full, [True, False, True, False, True, False])
    np.testing.assert_array_equal(after[~full], before[~full])
    assert np.isnan(spread[~full]).all()
    for row in np.flatnonzero(full):
        expected = []
        for slot in range(slots):
            x = S[row].copy(); x[slot] = 0
            flags = np.zeros(slots); flags[slot] = 1
            pred = []
            for j, fit in enumerate(fits):
                features = np.concatenate([x[:, j], flags, [] if extra is None else extra[row, j:j+1]])
                pred.append(((features - fit['mu']) / fit['sd']) @ fit['coef'] + fit['intercept'])
            expected.append(pred)
        np.testing.assert_allclose(after[row], np.mean(expected, axis=0))
        np.testing.assert_allclose(spread[row], np.std(expected, axis=0))
    np.testing.assert_array_equal(S, original_S)
    np.testing.assert_array_equal(missing, original_missing)
    if corrected: np.testing.assert_array_equal(extra, original_extra)
    if corrected:
        # With per-slot baselines, mask s uses the full rows' Model 2 scores rebuilt WITHOUT slot s.
        by_slot = [rng.normal(size=(int(full.sum()), k)) for _ in range(slots)]
        _, after_slot, _, _ = r.leave_one_slot_out(fits, S, missing, extra, by_slot)
        for n, row in enumerate(np.flatnonzero(full)):
            expected = []
            for slot in range(slots):
                x = S[row].copy(); x[slot] = 0; flags = np.zeros(slots); flags[slot] = 1
                expected.append([((np.concatenate([x[:, j], flags, by_slot[slot][n, j:j + 1]]) - f['mu']) / f['sd']) @ f['coef'] + f['intercept']
                                 for j, f in enumerate(fits)])
            np.testing.assert_allclose(after_slot[row], np.mean(expected, axis=0))
    missing[:, 0] = True
    before, after, spread, full = r.leave_one_slot_out(fits, S, missing, extra)
    assert not full.any() and np.isnan(spread).all()
    np.testing.assert_array_equal(after, before)


def make_world(root, rng, n_genes=40, n_kd=260, k_true=3, h1_shift=0.0):
    """Two programs every line shares, T cells seeing them at 40% strength, and a third program whose size is
    different in every line (it does not transfer); K562 misses 5 genes and T cells 3 inside the program set."""
    genes = [f'G{i}' for i in range(n_genes)]
    pd.DataFrame({'gene_name': genes}).to_csv(root / 'genes.csv', index=False)
    targets = [f'G{i}' for i in range(8)] + [f'T{i}' for i in range(n_kd - 8)]
    W = np.linalg.qr(rng.normal(size=(n_genes, k_true)))[0] * 3
    u = rng.normal(size=(n_kd, k_true))
    keep_rate = dict(HCT116=0.8, HEK293T=0.7, K562_gwps=0.8, CD4_Stim48hr=0.8, H1=0.3, RPE1=0.4, HepG2=0.4, Jurkat=0.5)
    low = dict(K562_gwps=range(35, 40), CD4_Stim48hr=range(10, 13))
    for name, background in zip(NAMES, BACKGROUNDS):
        mix = np.array([1.0, 1.0, rng.normal(0, 1.5)]) * (0.4 if name == 'CD4_Stim48hr' else 1.0)
        lfc = (u * mix) @ W.T + rng.normal(0, 0.05, (n_kd, n_genes)) + (h1_shift if name == 'H1' else 0)
        for i in range(8): lfc[i, i] = -2.5                  # knocked-down genes drop hard
        se = np.full_like(lfc, 0.05); cpm = np.full(n_genes, 50.0); cpm[list(low.get(name, []))] = 0.5
        kept = rng.random(n_kd) < keep_rate[name]; kept[:8] = True
        obs = pd.DataFrame(dict(challenge_target=False, n_cells=200, ontarget_lfc=-2.5, padj=0.01), index=np.array(targets)[kept])
        ad.AnnData(lfc[kept].astype(np.float32), obs=obs, var=pd.DataFrame(dict(ctrl_cpm=cpm), index=genes),
                   layers=dict(se=se[kept].astype(np.float32))).write_h5ad(root / 'filters' / f'{name}.h5ad')
        panel_lfc = np.where(kept[:, None] & (cpm >= 1)[None], lfc, np.nan).astype(np.float32)
        panel_se = np.where(np.isfinite(panel_lfc), se, np.nan).astype(np.float32)
        pobs = pd.DataFrame(dict(panel='ours', status=np.where(kept, 'kept', 'not screened'), n_cells=200.0,
                                 ontarget_lfc=np.where(kept, -2.5, np.nan), padj=0.01), index=pd.Index(targets, name='target'))
        (root / 'sources' / name).mkdir(parents=True, exist_ok=True)
        ad.AnnData(X=None, obs=pobs, var=pd.DataFrame(dict(ctrl_cpm=cpm), index=pd.Index(genes, name='gene')),
                   layers=dict(lfc=panel_lfc, se=panel_se), uns=dict(background=background)).write_h5ad(root / 'sources' / name / 'effects_panel.h5ad')
    return genes, targets


@pytest.fixture(scope='module')
def world(scripts_v3, tmp_path_factory):
    root = tmp_path_factory.mktemp('world'); (root / 'filters').mkdir()
    genes, targets = make_world(root, np.random.default_rng(0))
    meta, r = scripts_v3['vcc2_meta'], scripts_v3['vcc3_ridge']
    pd.DataFrame({'target_gene': ['G1', 'T5', 'T6', 'NOWHERE']}).to_csv(root / 'targets.csv', index=False)
    meta.fit(root / 'sources', root / 'genes.csv', root / 'm2', 'H1', meta.SOURCES)
    meta.predict(root / 'm2', root / 'targets.csv', ['H1'], root / 'm2pred')
    r.fit(root / 'filters', root / 'sources', root / 'genes.csv', root / 'ridge', 'H1', (2, 3), 0)   # grid only here; k80 has its own test
    report = r.predict(root / 'ridge', root / 'm2pred/H1', ['H1'], root / 'pred')
    return dict(root=root, genes=genes, targets=targets, report=report)


def test_fit_learns_transfer_and_corrects_model2(world):
    root = world['root']; model = json.loads((root / 'ridge/model.json').read_text())
    assert model['lines'] == TRAINING and model['heldout_background'] == 'H1'
    assert model['n_program_genes'] == 35                     # K562 misses G35-G39; T cells' gaps stay inside the set
    cv = pd.read_csv(root / 'ridge/cv_report.csv'); total = cv[(cv.destination == 'all') & (cv.k == 3)].set_index('variant')
    # Averaging keeps the weak T-cell copy and the non-transferring program; the rotating ridge learns to rescale
    # and drop them. It has no destination features, so it cannot learn any one line's own mixing.
    assert total.loc['programs', 'gene_mse'] < 0.5 * total.loc['programs', 'no_change_mse']
    assert total.loc['corrected', 'gene_mse'] < 0.8 * total.loc['corrected', 'model2_mse']
    assert set(cv.destination) == set(TRAINING) | {'all'}


def test_predict_outputs_and_gene_rules(scripts_v3, world):
    root, genes = world['root'], world['genes']; d = root / 'pred/H1'
    report = json.loads((d / 'prediction_report.json').read_text())
    assert report['schema'] == 'vcc_2026_3_ridge_v1' and report['heldout_background'] == 'H1' and 'H1' not in report['source_backgrounds']
    assert report['models'] == ['re', 'programs_k2', 'corrected_k2', 'programs_k3', 'corrected_k3', 'zero']
    np.testing.assert_array_equal(np.load(d / 'zero.npy'), np.zeros_like(np.load(d / 're.npy')))
    assert not (d / 'programs_k3_se.npy').exists() and not (d / 'corrected_k3_se.npy').exists()
    for name in report['models'][1:-1]:
        spread = np.load(d / f'{name}_spread.npy')
        assert spread.shape == (4, int(name.split('_k')[1]))
        affected = np.isfinite(spread).all(axis=1)
        check = report['leave_one_slot_out'][name]
        assert check['affected_rows'] == int(affected.sum())
        assert check['affected_targets'] == np.load(d / 'targets.npy')[affected].tolist()
        assert check['mean_abs_delta_after'] == scripts_v3['vcc3_ridge'].mean_abs_delta(np.load(d / f'{name}.npy')[affected])
        assert check['after_smaller'] == (check['mean_abs_correction_after'] < check['mean_abs_correction_before'])
        assert (check['corrected_baseline'] is not None) == name.startswith('corrected')
    for name in ['re', 're_se', 'n_backgrounds', 'ontarget_lfc', 'targets', 'genes']:
        np.testing.assert_array_equal(np.load(d / f'{name}.npy'), np.load(root / 'm2pred/H1' / f'{name}.npy'))
    re, prog, corr = (np.load(d / f'{x}.npy') for x in ['re', 'programs_k3', 'corrected_k3'])
    outside = np.arange(35, 40)
    assert np.isnan(prog[:, outside]).all()                   # programs only: no prediction outside the set
    np.testing.assert_array_equal(corr[1:3, outside], re[1:3, outside])   # corrected: Model 2 kept there
    assert np.isnan(prog[0, 1]) and np.isnan(corr[0, 1])      # G1's own column: ontarget_lfc handles it
    assert np.isfinite(prog[1:3, :35]).all() and not np.allclose(corr[1:3, :35], re[1:3, :35])
    assert np.isnan(prog[3]).all() and np.isnan(corr[3]).all() and report['no_information'] == ['NOWHERE']


def test_holdout_and_destination_never_leak(scripts_v3, world, tmp_path):
    r = scripts_v3['vcc3_ridge']; root = world['root']
    # Shifting every H1 effect changes nothing about the H1-held-out fit.
    (tmp_path / 'filters').mkdir()
    make_world(tmp_path, np.random.default_rng(0), h1_shift=4.0)
    r.fit(tmp_path / 'filters', tmp_path / 'sources', tmp_path / 'genes.csv', tmp_path / 'ridge', 'H1', (2, 3), 0)
    a, b = np.load(root / 'ridge/ridge.npz'), np.load(tmp_path / 'ridge/ridge.npz')
    assert set(a.files) == set(b.files) and all(np.array_equal(a[f], b[f]) for f in a.files)
    # A destination's Model 2 baseline and its input slots never see its own values.
    pgenes = np.load(root / 'ridge/program_genes.npy').astype(str)
    lines = [r.Line(root / 'filters', n, pgenes) for n in NAMES[:4]]
    dest, others = lines[0], lines[1:]
    tau2 = np.full(len(pgenes), 0.01)
    before = r.model2_baseline(dest, others, tau2)
    dest.lfc += 9.0
    np.testing.assert_array_equal(before, r.model2_baseline(dest, others, tau2))
    W = np.load(root / 'ridge/loadings.npy')[:, :3]
    cache = {l.name: r.scores(l.lfc, W, l.own) for l in lines}
    S, missing = r.slot_inputs(dest.targets, lines, W, cache, exclude=dest.name)
    assert missing[:, 0].all() and (S[:, 0] == 0).all()


def test_background_mismatch_rejected_before_fit_and_predict(scripts_v3, world, monkeypatch, tmp_path):
    r = scripts_v3['vcc3_ridge']; root = world['root']
    r.check_backgrounds(root / 'sources')
    monkeypatch.setitem(r.BACKGROUND, 'K562_gwps', 'incorrect label')
    message = "K562_gwps: panel says 'K562', vcc3 says 'incorrect label'"
    with pytest.raises(ValueError, match=message):
        r.fit(tmp_path / 'missing_filters', root / 'sources', root / 'genes.csv', tmp_path / 'fit')
    with pytest.raises(ValueError, match=message):
        r.predict(root / 'ridge', root / 'm2pred/H1', ['H1'], tmp_path / 'pred')
    assert not (tmp_path / 'fit').exists() and not (tmp_path / 'pred').exists()


def test_generator_minimum_size_without_se_gate(scripts_v3):
    c, g = scripts_v3['vcc2_common'], scripts_v3['vcc3_generate']
    line = g.Line(sp.csr_matrix([[100, 100, 100, 100]] * 2), None, ['KD', 'G1', 'G2', 'LOW'])
    line.gate = c.expression_gate([10, 10, 10, 1])
    D = np.array([[-2.0, 0.6, 0.15, 0.6]])
    np.testing.assert_allclose(g.calibrate(D, line, ['KD'], g.Cfg(argparse.Namespace(alpha=2.0, min_effect=0.2))), [[0, 1.2, 0, 0]])
    with pytest.raises(ValueError, match='needs the prediction SE'):
        g.calibrate(D, line, ['KD'], g.Cfg(argparse.Namespace(z=2.0)))
    se = np.array([[0.1, 0.4, 0.01, 0.01]])
    np.testing.assert_allclose(g.calibrate(D, line, ['KD'], g.Cfg(argparse.Namespace(z=2.0, min_effect=0.1)), se), [[0, 0, 0.15, 0]])


def test_screen_metrics_and_paired_bootstrap(scripts_v3):
    s = scripts_v3['vcc3_screen']
    truth = np.tile(np.arange(1, 81, dtype=float), (3, 1))
    pred = truth[:, ::-1].copy(); pred[1] = 0; pred[2] = truth[2]
    eligible = np.ones(truth.shape, bool); sig = eligible.copy(); sig[2] = False
    metrics = s.target_metrics(pred, truth, eligible, sig)
    np.testing.assert_allclose(metrics[:, 1], [-1, 0, 1])
    np.testing.assert_allclose(metrics[:, 2], [0.4, 0, 1])               # an all-zero prediction earns no overlap
    assert metrics[0, 3] == 1 and np.isnan(metrics[1, 3]) and np.isnan(metrics[2, 3])   # accuracy only where predicted
    np.testing.assert_allclose(metrics[:2, 4], [1, 0]); assert np.isnan(metrics[2, 4])  # coverage kept apart
    sparse = np.zeros_like(truth); sparse[0, -10:] = truth[0, -10:]           # 10 nonzero genes, all in the true top 50
    one = s.target_metrics(sparse[:1], truth[:1], eligible[:1], sig[:1])
    assert one[0, 2] == pytest.approx(10 / 50) and one[0, 3] == 1 and one[0, 4] == pytest.approx(10 / 80)
    np.testing.assert_allclose(metrics[:, 0], [np.mean((np.arange(1, 81) * 2 - 81) ** 2), np.mean(truth[1] ** 2), 0])
    reference = metrics.copy(); reference[:, 0] += [1, 2, 3]
    draws = np.array([[0, 0, 0], [0, 1, 2], [2, 2, 2], [1, 1, 1]])
    result = s.paired_summary(metrics, reference, draws)
    assert result['mse_diff'] == -2 and result['direction_accuracy_n_targets'] == 1 and result['direction_coverage_n_targets'] == 2
    np.testing.assert_allclose([result['mse_diff_ci_low'], result['mse_diff_ci_high']], np.percentile([-1, -2, -3, -2], [2.5, 97.5]))
    assert result['direction_coverage_diff_ci_low'] == result['direction_coverage_diff_ci_high'] == 0


def test_screen_alignment_grid_and_stop_check(scripts_v3, tmp_path):
    s = scripts_v3['vcc3_screen']
    genes = np.array(['KD0', 'KD1', 'KD2', 'LOW', 'FIVE', 'NAN'] + [f'G{i}' for i in range(58)])
    targets = np.array(['KD0', 'KD1', 'KD2', 'DROP'])
    truth = np.tile(np.linspace(-1, 1, len(genes)), (len(targets), 1)).astype(np.float32)
    for i in range(3): truth[i, i] = -99
    truth[:, 3:5] = 99; truth[:, 5] = np.nan
    se = np.full_like(truth, 0.05)
    order_t = np.array([2, 3, 0, 1]); order_g = np.arange(len(genes))[::-1]
    obs = pd.DataFrame({'status': ['kept', 'kept', 'kept', 'failed: cells']}, index=targets)
    ad.AnnData(X=None, obs=obs.iloc[order_t], var=pd.DataFrame(index=genes[order_g]),
               layers=dict(lfc=truth[np.ix_(order_t, order_g)], se=se[np.ix_(order_t, order_g)])).write_h5ad(tmp_path / 'truth.h5ad')
    cpm = np.full(len(genes), 10.0); cpm[3:5] = [4.9, 5.0]
    ad.AnnData(np.zeros((1, len(genes)-1)), var=pd.DataFrame({'ctrl_cpm': cpm[:-1][::-1]}, index=genes[:-1][::-1])).write_h5ad(tmp_path / 'filters.h5ad')
    pd.DataFrame({'target_gene': ['KD1', 'DROP', 'KD0', 'KD2']}).to_csv(tmp_path / 'targets.csv', index=False)
    p = tmp_path / 'pred'; p.mkdir()
    np.save(p / 'genes.npy', genes); np.save(p / 'targets.npy', targets)
    np.save(p / 're.npy', truth); np.save(p / 'programs_k2.npy', truth * 0.8)
    np.save(p / 'corrected_k2.npy', truth * 1.1); np.save(p / 'zero.npy', np.zeros_like(truth))
    np.save(p / 'corrected_k2_spread.npy', np.ones((4, 2)))  # diagnostic, not a candidate
    (p / 'prediction_report.json').write_text(json.dumps(dict(heldout_background='H1', source_backgrounds=['K562'])))
    args = (p, tmp_path / 'truth.h5ad', tmp_path / 'filters.h5ad', tmp_path / 'targets.csv')
    g, t, measured, _, eligible, _, gate, info = s.load_inputs(*args)
    np.testing.assert_array_equal(g, genes); np.testing.assert_array_equal(t, ['KD1', 'KD0', 'KD2'])
    np.testing.assert_array_equal(measured, truth[[1, 0, 2]])
    assert not gate[3:5].any() and not gate[-1]
    assert not eligible[:, 5].any() and not eligible[[0, 1, 2], [1, 0, 2]].any()
    assert info['excluded_not_kept'] == ['DROP'] and info['genes_absent_from_filter'] == 1
    report = s.screen(*args, tmp_path / 'out', resamples=1000)
    frame = pd.read_csv(tmp_path / 'out/screen.csv')
    assert len(frame) == 4 * 5 * 3 and set(frame.delta_array) == {'re', 'programs_k2', 'corrected_k2', 'zero'}
    assert set(frame.alpha) == set(s.ALPHAS) and set(frame.min_effect) == set(s.MINIMUMS)
    baseline = frame.set_index('arm').loc['re_a1_m0']
    assert baseline.mse == 0 and baseline.mse_diff == baseline.mse_diff_ci_low == baseline.mse_diff_ci_high == 0
    assert report['n_targets'] == 3 and report['zero_check_passed']      # every shortlisted arm beats zero on Pearson and top-50
    selected = json.loads((tmp_path / 'out/shortlist.json').read_text())
    assert len(selected['arms']) in [3, 4] and selected['ready_for_benchmark']
    assert all(a['alpha'] == 1.0 for a in selected['arms'])                 # the screen chooses array and minimum size, not alpha
    zero = frame.set_index('arm').loc['zero_a1_m0']
    assert zero.pearson == 0 and zero.top50_overlap == 0 and zero.direction_coverage == 0
    np.testing.assert_array_equal(np.load(tmp_path / 'out/bootstrap_indices.npy'), np.random.default_rng(2026).integers(0, 3, size=(1000, 3)))
    # A missing truth name is a hard alignment error, not a silent zero or a dropped target.
    pd.DataFrame({'target_gene': ['ABSENT']}).to_csv(tmp_path / 'targets.csv', index=False)
    with pytest.raises(ValueError, match='missing'):
        s.load_inputs(*args)


def test_signal_k_counts_programs_above_noise(scripts_v3):
    """k80 sees through heteroscedastic noise: 3 true programs among 200 genes, noise SD 0.1-1 per gene."""
    from types import SimpleNamespace
    r = scripts_v3['vcc3_ridge']; rng = np.random.default_rng(5)
    N, G = 900, 200
    W = np.linalg.qr(rng.normal(size=(G, 3)))[0]
    signal = (rng.normal(size=(N, 3)) * [6.0, 4.0, 3.0]) @ W.T
    sd = np.tile(rng.uniform(0.1, 1.0, G), (N, 1))
    lines = [SimpleNamespace(lfc=(signal + rng.normal(size=(N, G)) * sd).astype(np.float32), se=sd.astype(np.float32), own=np.full(N, -1))]
    result = r.signal_k(lines, 0.8)
    assert 3 <= result['n_signal_components'] <= 12 and result['k80'] in (2, 3)
    assert 0.2 < result['noise_share_of_total'] < 0.9


def test_em_pca_recovers_programs_with_missing_genes(scripts_v3):
    """Lines that miss some genes still teach the programs: EM-PCA beats treating the gaps as zero effects."""
    from types import SimpleNamespace
    r = scripts_v3['vcc3_ridge']; rng = np.random.default_rng(9)
    G, k = 60, 3
    W = np.linalg.qr(rng.normal(size=(G, k)))[0]
    def line(n, missing):
        Y = (rng.normal(size=(n, k)) * [5.0, 3.0, 2.0]) @ W.T + rng.normal(0, 0.05, (n, G))
        Y[:, missing] = np.nan
        return SimpleNamespace(lfc=Y.astype(np.float32), own=np.full(n, -1))
    lines = [line(80, []), line(300, rng.choice(G, 12, replace=False))]           # the bigger line misses 20% of genes
    V, _, _, filled = r.fit_programs(lines, k)
    angle = lambda A: np.linalg.norm(W - A @ (A.T @ W))                          # distance from the true program space
    zero = [SimpleNamespace(lfc=np.nan_to_num(l.lfc), own=l.own) for l in lines]
    V0 = r.fit_programs(zero, k)[0]
    assert angle(V) < 0.05 and angle(V) < 0.5 * angle(V0)
    assert np.isfinite(filled).all()


def test_gate_from_and_nochange(scripts_v3, tmp_path):
    c, g = scripts_v3['vcc2_common'], scripts_v3['vcc3_generate']
    line = g.Line(sp.csr_matrix([[100, 100, 100, 100]] * 2), None, ['KD', 'G1', 'G2', 'G3'])
    line.gate = c.expression_gate([10, 10, 10, 10])
    D = np.array([[-2.0, 0.6, 0.6, 0.05]])                                # e.g. a ridge array without its own SE
    info = dict(gate_delta=np.array([[-2.0, 0.9, 0.1, 0.9]]), gate_se=np.array([[0.1, 0.2, 0.2, 0.2]]))
    cfg = g.Cfg(argparse.Namespace(gate_from='re', z=3.0, min_effect=0.1, alpha=2.0))
    # G1: Model 2 calls it -> ridge value x alpha; G2: Model 2 does not; G3: called but below the minimum size.
    np.testing.assert_allclose(g.calibrate(D, line, ['KD'], cfg, info=info), [[0, 1.2, 0, 0]])
    with pytest.raises(ValueError, match='invalid generator settings'):
        g.Cfg(argparse.Namespace(gate_from='re', z=0.0))
    genes, targets = ['KD', 'G1'], ['KD']
    D0, info0 = g.load_effects(tmp_path, targets, genes, 'nochange', g.Cfg(argparse.Namespace()))
    assert not D0.any() and info0 == {} and g.is_nochange('nochange') and not g.is_nochange('zero')


def test_official_finalists_benchmark(scripts_v3, tmp_path):
    """Real CPU evaluator on a small synthetic H1: leaderboard-scored finalists, Model 2 gating, zero and nochange."""
    binary = Path(sys.executable).parent / 'cell-eval2'
    if not binary.exists(): pytest.skip('official cell-eval2 is not installed')
    g = scripts_v3['vcc3_generate']; genes = np.array(['K0', 'K1', 'K2'] + [f'G{i}' for i in range(37)])
    rng = np.random.default_rng(2); baseline = rng.integers(30, 100, len(genes))
    controls = rng.poisson(baseline, (160, len(genes))); groups = {}
    for i, kd in enumerate(genes[:3]):
        factors = np.ones(len(genes)); factors[3 + i::3] = 3; factors[i] = 0.2
        groups[kd] = rng.poisson(baseline * factors, (100, len(genes)))
    labels = ['non-targeting'] * len(controls) + [t for t, x in groups.items() for _ in range(len(x))]
    ad.AnnData(sp.csr_matrix(np.vstack([controls] + list(groups.values()))),
               obs=pd.DataFrame({'target_gene': labels}, index=[str(i) for i in range(len(labels))]),
               var=pd.DataFrame(index=genes)).write_h5ad(tmp_path / 'h1.h5ad')
    pd.DataFrame({'gene_name': genes}).to_csv(tmp_path / 'genes.csv', index=False)
    pd.DataFrame({'target_gene': list(groups)}).to_csv(tmp_path / 'targets.csv', index=False)
    effects = np.zeros((3, len(genes)))
    for i in range(3): effects[i, 3 + i::3] = 0.75; effects[i, i] = -1.5
    pred = tmp_path / 'pred'; pred.mkdir()
    for name, value in [('genes', genes), ('targets', list(groups)), ('re', effects), ('re_se', np.full(effects.shape, 0.1)),
                        ('programs_k3', effects * 0.8), ('zero', np.zeros_like(effects)), ('ontarget_lfc', [-1.5] * 3)]:
        np.save(pred / f'{name}.npy', np.asarray(value))
    (pred / 'prediction_report.json').write_text(json.dumps(dict(schema='vcc_2026_3_ridge_v1', models=['re', 'programs_k3', 'zero'],
        heldout_background='H1', source_backgrounds=['HCT116', 'HEK293T'])))
    arms = [dict(method='programs_k3', call='z', z=0.0, min_effect=0.1, alpha=2.0, label='programs_alpha2'),
            dict(method='re', call='z', z=3.0, min_effect=0.1, alpha=2.0, label='gated_model2'),
            dict(method='programs_k3', call='z', gate_from='re', z=3.0, min_effect=0.1, alpha=2.0, label='programs_gated_by_model2'),
            dict(method='zero', call='z', alpha=1.0, label='ontarget_only'),
            dict(method='nochange', call='z', alpha=1.0, label='true_nochange')]
    (tmp_path / 'arms.json').write_text(json.dumps(arms))
    common = ['benchmark', '--cells', str(tmp_path / 'h1.h5ad'), '--genes', str(tmp_path / 'genes.csv'), '--targets', str(tmp_path / 'targets.csv'),
              '--predictions', str(pred), '--cell-eval', str(binary), '--num-threads', '1', '--n-cells', '100', '--min-cells', '50',
              '--anchor-splits', '2']
    (tmp_path / 'bad.json').write_text(json.dumps([dict(method='programs_k3', call='z', z=2.0, alpha=1.0)]))
    refused = subprocess.run([sys.executable, g.__file__, *common, '--variants', str(tmp_path / 'bad.json'), '--out', str(tmp_path / 'bad')],
                             capture_output=True, text=True, timeout=120)
    assert refused.returncode != 0 and 'needs programs_k3_se.npy' in refused.stderr and not (tmp_path / 'bad/real.h5ad').exists()
    out = tmp_path / 'final'
    result = subprocess.run([sys.executable, g.__file__, *common, '--variants', str(tmp_path / 'arms.json'), '--margin-over', 'gated_model2',
                             '--keep-arm-files', '--out', str(out)], capture_output=True, text=True, timeout=900)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-3000:]
    scaled = pd.read_csv(out / 'scaled_scores.csv', index_col='arm')
    assert set(scaled.index) == {a['label'] for a in arms} and scaled['overall'].notna().all() and scaled['overall'].is_monotonic_decreasing
    for arm in scaled.index: assert (out / arm / 'per_target.csv').exists() and (out / arm / 'scaled.csv').exists()
    selected = json.loads((out / 'selected_generator.json').read_text())
    assert selected['arm'] == scaled.index[0] and selected['selection'] == 'scaled overall (leaderboard formula)'
    assert selected['margin_over'] == 'gated_model2'
    assert selected['margin'] == pytest.approx(scaled['overall'].iloc[0] - scaled.loc['gated_model2', 'overall'])
    # zero keeps the knockdown, nochange does not: the target gene of the nochange arm stays at its control level.
    def target_mean(arm, target):
        a = ad.read_h5ad(out / f'{arm}.h5ad'); x = a[a.obs.target_gene == target].X
        return float(np.asarray(x[:, list(genes).index(target)].mean()))
    control = float(controls[:, 0].mean())
    assert target_mean('ontarget_only', 'K0') < 0.5 * control and target_mean('true_nochange', 'K0') > 0.7 * control
    # --resume reuses this folder's real cells, cache, ceiling and bundle: nothing is rebuilt.
    stamps = {x: (out / x).stat().st_mtime for x in ('real.h5ad', 'bundle/manifest.json', 'ceiling/ceiling_agg.csv')}
    (tmp_path / 'more.json').write_text(json.dumps([dict(method='re', call='topn', top_n=6, detect_floor=True, z_req=5.0, alpha=1.0, label='top6')]))
    result = subprocess.run([sys.executable, g.__file__, *common, '--variants', str(tmp_path / 'more.json'), '--resume', '--out', str(out)],
                            capture_output=True, text=True, timeout=900)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-3000:]
    assert 'reference reused' in result.stdout and {x: (out / x).stat().st_mtime for x in stamps} == stamps
    assert pd.read_csv(out / 'scaled_scores.csv', index_col='arm').index.tolist() == ['top6']

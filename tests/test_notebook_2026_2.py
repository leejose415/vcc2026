"""Behavioral checks for scripts embedded in the random-effects meta-analysis notebook."""
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


@pytest.fixture(scope='session')
def scripts_v2(tmp_path_factory):
    notebook = nbformat.read(ROOT / 'vcc_2026_2.ipynb', as_version=4)
    directory = tmp_path_factory.mktemp('meta_scripts')
    for cell in notebook.cells:
        if cell.source.startswith('%%writefile'):
            header, body = cell.source.split('\n', 1)
            (directory / Path(header.split()[-1]).name).write_text(body)
    sys.path.insert(0, str(directory))
    return {name: importlib.import_module('vcc2_' + name) for name in ['common', 'meta', 'generate']}


def test_notebook_structure_and_commands(scripts_v2):
    from IPython.core.inputtransformer2 import TransformerManager
    notebook = nbformat.read(ROOT / 'vcc_2026_2.ipynb', as_version=4)
    nbformat.validate(notebook)
    transform = TransformerManager()
    written = [c.source.split('\n', 1)[0].split('/')[-1] for c in notebook.cells if c.source.startswith('%%writefile')]
    assert written == ['vcc2_common.py', 'vcc2_meta.py', 'vcc2_generate.py']
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
    for name, module in scripts_v2.items():
        if name == 'common': continue
        result = subprocess.run([sys.executable, module.__file__, '--help'], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr


def test_worked_example(scripts_v2):
    """Step 4's example: K562 sits out, tau2 = 0.04, weights 20 / 18.4 / 9.8 (the spec rounds SE^2 first)."""
    m = scripts_v2['meta']
    lfc = np.array([-0.80, -0.50, np.nan, -1.10])[:, None, None]
    se = np.array([0.10, 0.12, np.nan, 0.25])[:, None, None]
    re, re_se, k = (x[0, 0] for x in m.combine(lfc, se, np.array([0.04])))
    w = 1 / (np.array([0.10, 0.12, 0.25]) ** 2 + 0.04)
    np.testing.assert_allclose(w, [20.0, 18.5, 9.7], atol=0.15)
    assert k == 3
    assert re == pytest.approx(-0.75, abs=0.005) and re_se == pytest.approx(0.14, abs=0.005)
    assert re == pytest.approx((w * [-0.8, -0.5, -1.1]).sum() / w.sum()) and re_se == pytest.approx(np.sqrt(1 / w.sum()))
    # Large tau2 flattens the weights toward the plain average; tau2 = 0 is inverse-variance weighting.
    assert m.combine(lfc, se, np.array([100.0]))[0][0, 0] == pytest.approx(-0.8, abs=1e-3)
    inverse = 1 / np.array([0.10, 0.12, 0.25]) ** 2
    assert m.combine(lfc, se, np.array([0.0]))[0][0, 0] == pytest.approx((inverse * [-0.8, -0.5, -1.1]).sum() / inverse.sum())
    # A gene no background measured is "no information": NaN, count 0.
    nothing = m.combine(np.full((2, 1, 1), np.nan), np.full((2, 1, 1), np.nan), np.array([0.01]))
    assert np.isnan(nothing[0]).all() and np.isnan(nothing[1]).all() and nothing[2][0, 0] == 0


def test_tau2_recovers_heterogeneity_and_floors_at_zero(scripts_v2):
    m = scripts_v2['meta']; rng = np.random.default_rng(0)
    S, T = 4, 3000; se = np.full((S, T, 2), 0.2)
    truth = rng.normal(0, 0.5, (1, T, 2))
    lfc = truth + rng.normal(0, np.sqrt(0.09), (S, T, 2)) * [1, 0] + rng.normal(0, 0.2, (S, T, 2))
    lfc[1:, ::2, :] = np.nan                                  # targets with one background add nothing
    tau2, df = m.tau2_dl(lfc, se, np.full(T, -1))
    assert tau2[0] == pytest.approx(0.09, rel=0.15) and tau2[1] == pytest.approx(0, abs=0.004)
    assert df[0] == (S - 1) * (T // 2)
    # The knocked-down gene's own column never enters tau2.
    wild = lfc.copy(); wild[:, :, 1] = rng.normal(0, 5, (S, T)); own = np.full(T, 1)
    assert m.tau2_dl(wild, se, own)[1][1] == 0 and np.isnan(m.tau2_dl(wild, se, own)[0][1])


def write_panel(root, name, background, targets, genes, lfc, se, status, ontarget):
    obs = pd.DataFrame(dict(panel='ours', status=status, n_cells=200.0, ontarget_lfc=ontarget, padj=0.01),
                       index=pd.Index(targets, name='target'))
    a = ad.AnnData(X=None, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name='gene')),
                   layers=dict(lfc=np.asarray(lfc, np.float32), se=np.asarray(se, np.float32)), uns=dict(background=background))
    (root / name).mkdir(parents=True, exist_ok=True); a.write_h5ad(root / name / 'effects_panel.h5ad')


def make_panels(tmp_path, h1_value=7.0):
    genes = ['KD', 'G1', 'G2', 'G3']; targets = ['KD', 'T2']; root = tmp_path / 'sources'; tmp_path.mkdir(exist_ok=True)
    pd.DataFrame({'gene_name': genes}).to_csv(tmp_path / 'genes.csv', index=False)
    nan = np.nan
    write_panel(root, 'A', 'BG_A', targets, genes, [[-2, 1.0, 0.2, nan], [nan] * 4], [[0.1, 0.1, 0.2, nan], [nan] * 4],
                ['kept', 'failed: weak knockdown'], [-2.0, -0.2])
    write_panel(root, 'B', 'BG_B', targets, genes, [[-3, 0.0, nan, nan], [0, 0.5, 0.5, nan]], [[0.1, 0.2, nan, nan], [0.1] * 3 + [nan]],
                ['kept', 'kept'], [-1.0, -1.5])
    write_panel(root, 'H1', 'H1', targets, genes, [[h1_value] * 4] * 2, [[0.01] * 4] * 2, ['kept', 'kept'], [-9.0, -9.0])
    write_panel(root, 'B_rest', 'BG_B', targets, genes, [[5] * 4] * 2, [[0.1] * 4] * 2, ['kept', 'kept'], [-1.0, -1.0])
    pd.DataFrame({'target_gene': ['KD', 'MISSING']}).to_csv(tmp_path / 'targets.csv', index=False)
    return root, genes


def test_fit_predict_one_vote_per_background_and_holdout(scripts_v2, tmp_path):
    m = scripts_v2['meta']; root, genes = make_panels(tmp_path)
    with pytest.raises(ValueError, match='one vote per background'):
        m.fit(root, tmp_path / 'genes.csv', tmp_path / 'bad', 'H1', ['A', 'B', 'B_rest', 'H1'])
    model = m.fit(root, tmp_path / 'genes.csv', tmp_path / 'model', 'H1', ['A', 'B', 'H1'])
    assert model['backgrounds'] == ['BG_A', 'BG_B'] and model['heldout_background'] == 'H1'
    report = m.predict(tmp_path / 'model', tmp_path / 'targets.csv', ['A', 'B'], tmp_path / 'pred')
    out = tmp_path / 'pred/A'
    re, count = (np.load(out / f'{x}.npy') for x in ['re', 'n_backgrounds'])
    assert sorted(p.name for p in out.iterdir()) == ['genes.npy', 'n_backgrounds.npy', 'ontarget_lfc.npy', 'prediction_report.json',
                                                     're.npy', 're_se.npy', 'targets.npy']
    tau2 = np.load(tmp_path / 'model/tau2.npy')
    assert np.isnan(tau2[0]) and np.isnan(tau2[3])             # own column excluded / never measured twice
    np.testing.assert_array_equal(count, [[0, 2, 1, 0], [0, 0, 0, 0]])
    assert np.isnan(re[0, 0])                                  # own gene: not averaged
    t2 = tau2[1]; w = 1 / (np.array([0.1, 0.2]) ** 2 + t2)
    assert re[0, 1] == pytest.approx((w * [1.0, 0.0]).sum() / w.sum())   # weighted, not (1 + 0) / 2
    assert re[0, 2] == pytest.approx(0.2)                      # B has no value: sits out
    assert np.isnan(re[0, 3]) and np.isnan(re[1]).all()        # no information anywhere; MISSING not in the panel
    np.testing.assert_allclose(np.load(out / 'ontarget_lfc.npy'), [-1.5, np.nan])
    assert report['not_in_panel'] == ['MISSING'] and report['no_source_targets'] == ['MISSING']
    assert report['heldout_background'] == 'H1' and 'H1' not in report['source_backgrounds']
    for name in ['re', 're_se', 'n_backgrounds', 'ontarget_lfc']:
        np.testing.assert_array_equal(np.load(out / f'{name}.npy'), np.load(tmp_path / 'pred/B' / f'{name}.npy'))
    m.predict(tmp_path / 'model', tmp_path / 'targets.csv', ['A'], tmp_path / 'pred')   # re-running replaces it
    np.testing.assert_array_equal(np.load(out / 're.npy'), re)
    (tmp_path / 'pred/other').mkdir(); (tmp_path / 'pred/other/keep.txt').write_text('not ours')
    with pytest.raises(ValueError, match='empty prediction directory'):
        m.predict(tmp_path / 'model', tmp_path / 'targets.csv', ['other'], tmp_path / 'pred')
    # The held-out background's values never reach tau2 or the averages.
    root2, _ = make_panels(tmp_path / 'other', h1_value=-50.0)
    m.fit(root2, tmp_path / 'genes.csv', tmp_path / 'model2', 'H1', ['A', 'B', 'H1'])
    np.testing.assert_array_equal(np.load(tmp_path / 'model2/tau2.npy'), tau2)
    assert m.exam(root, 'H1', tmp_path / 'exam.csv', ['A', 'B', 'H1']) == ['KD', 'T2']


def test_de_call_alpha_gate_and_own_column(scripts_v2):
    c, g = scripts_v2['common'], scripts_v2['generate']
    line = g.Line(sp.csr_matrix([[100, 100, 100, 100, 100]] * 2), None, ['KD', 'G1', 'G2', 'G3', 'LOW'])
    line.gate = c.expression_gate([10, 10, 10, 10, 1])
    D = np.array([[-2.0, 0.6, 0.15, 0.6, 0.6]]); se = np.array([[0.1, 0.2, 0.01, np.inf, 0.01]])
    cfg = g.Cfg(argparse.Namespace(alpha=2.0, z=2.0, min_effect=0.2))
    # G1 passes both and is doubled at full size; G2 too small; G3 unknown SE; LOW gated; KD handled separately.
    np.testing.assert_allclose(g.calibrate(D, line, ['KD'], cfg, se), [[0, 1.2, 0, 0, 0]])
    cfg.z = 4.0
    np.testing.assert_allclose(g.calibrate(D, line, ['KD'], cfg, se), [[0, 0, 0, 0, 0]])
    off = g.Cfg(argparse.Namespace(alpha=1.0))
    np.testing.assert_allclose(g.calibrate(D, line, ['KD'], off), [[0, 0.6, 0.15, 0.6, 0]])
    with pytest.raises(ValueError, match='needs the prediction SE'):
        g.calibrate(D, line, ['KD'], g.Cfg(argparse.Namespace(z=2.0)))


def test_binary_gate_uses_per_cell_cpm(scripts_v2):
    c, g = scripts_v2['common'], scripts_v2['generate']
    np.testing.assert_array_equal(c.expression_gate([0, 4.9, 5, np.nextafter(5., 6.), 6, np.nan]),
                                  [False, False, False, True, True, False])
    line = g.Line(sp.csr_matrix([[1, 999], [100, 999900]]), None, ['G', 'KD'])
    assert line.cpm[0] == pytest.approx(550)  # mean per cell, not pooled ~101


def write_prediction(pred, genes, targets, effects, se, ontarget, heldout='H1'):
    pred.mkdir()
    for name, value in [('genes', genes), ('targets', targets), ('re', effects), ('re_se', se),
                        ('ontarget_lfc', ontarget)]:
        np.save(pred / f'{name}.npy', np.asarray(value))
    (pred / 'prediction_report.json').write_text(json.dumps(dict(schema='vcc_2026_2_meta_v1', models=['re'],
        heldout_background=heldout, source_backgrounds=['BG_A', 'BG_B'])))


def test_count_generation_and_submission(scripts_v2, tmp_path):
    g, c = scripts_v2['generate'], scripts_v2['common']
    genes = ['KD', 'G']; pd.DataFrame({'gene_name': genes}).to_csv(tmp_path / 'gene_names.csv', index=False)
    pd.DataFrame({'target_gene': ['KD'], 'n_cells': [400]}).to_csv(tmp_path / 'pert_counts.csv', index=False)
    control = ad.AnnData(sp.csr_matrix([[30, 20]]*30), obs=pd.DataFrame({'target_gene': ['non-targeting']*30},
                         index=[str(i) for i in range(30)]), var=pd.DataFrame(index=genes))
    control.write_h5ad(tmp_path / 'context_A.h5ad')
    write_prediction(tmp_path / 'pred', genes, ['KD'], [[np.nan, 0.5]], [[np.nan, 0.1]], [-1], heldout=None)
    (tmp_path / 'settings.json').write_text(json.dumps(dict(schema='vcc_2026_2_meta_v1', delta_array='re',
        settings=dict(alpha=1.5, z=2.0, min_effect=0.1, pi=0, up_w=0.75, n_cells=400))))
    command = [sys.executable, g.__file__, 'submit', '--data', str(tmp_path), '--lines', 'A', '--delta', f'A={tmp_path / "pred"}',
               '--delta-array', 'selected', '--settings', str(tmp_path / 'settings.json'), '--out', str(tmp_path / 'output')]
    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    output = ad.read_h5ad(tmp_path / 'output/prediction.h5ad')
    assert output.shape == (400, 2)
    c.validate_counts(output.X)
    assert output.obs.context.unique().tolist() == ['A']
    report = json.loads((tmp_path / 'output/generator_report.json').read_text())
    assert report['delta_array'] == 're' and report['settings']['z'] == 2.0 and report['settings']['alpha'] == 1.5


def test_official_exam_benchmark(scripts_v2, tmp_path):
    """Run the real CPU evaluator on a small synthetic reference, never real H1."""
    binary = Path(sys.executable).parent / 'cell-eval2'
    if not binary.exists(): pytest.skip('official cell-eval2 is not installed')
    g = scripts_v2['generate']; genes = np.array(['K0', 'K1', 'K2'] + [f'G{i}' for i in range(37)])
    rng = np.random.default_rng(2); baseline = rng.integers(30, 100, len(genes))
    controls = rng.poisson(baseline, (120, len(genes))); groups = {}
    for i, kd in enumerate(genes[:3]):
        factors = np.ones(len(genes)); factors[3+i::3] = 3; factors[i] = 0.2
        groups[kd] = rng.poisson(baseline*factors, (100, len(genes)))
    labels = ['non-targeting']*len(controls) + [t for t, x in groups.items() for _ in range(len(x))]
    ad.AnnData(sp.csr_matrix(np.vstack([controls] + list(groups.values()))),
        obs=pd.DataFrame({'target_gene': labels}, index=[str(i) for i in range(len(labels))]),
        var=pd.DataFrame(index=genes)).write_h5ad(tmp_path / 'h1.h5ad')
    pd.DataFrame({'gene_name': genes}).to_csv(tmp_path / 'genes.csv', index=False)
    pd.DataFrame({'target_gene': list(groups)}).to_csv(tmp_path / 'targets.csv', index=False)
    # Independent synthetic source predictions: no perturbation outcomes used here.
    effects = np.zeros((3, len(genes)))
    for i in range(3): effects[i, 3+i::3] = 0.75; effects[i, i] = -1.5
    write_prediction(tmp_path / 'pred', genes, list(groups), effects, np.full(effects.shape, 0.1), [-1.5]*3)
    out = tmp_path / 'benchmark'
    command = [sys.executable, g.__file__, 'benchmark', '--cells', str(tmp_path / 'h1.h5ad'),
        '--genes', str(tmp_path / 'genes.csv'), '--targets', str(tmp_path / 'targets.csv'), '--predictions', str(tmp_path / 'pred'),
        '--out', str(out), '--cell-eval', str(binary), '--num-threads', '1', '--n-cells', '100', '--min-cells', '50',
        '--delta-array', 're', '--z-grid', '2', '--min-effect-grid', '0.1', '--alpha-grid', '1',
        '--pi-grid', '0', '--up-w-grid', '0.75', '--cell-rules', 'fixed']
    result = subprocess.run(command, capture_output=True, text=True, timeout=300)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (out / 're_000/agg_results.csv').exists() and not (out / 're_000.h5ad').exists()
    assert (out / 'ceiling/ceiling_agg.csv').exists()
    ranks = pd.read_csv(out / 'metric_ranks.csv')
    assert len(ranks) == 6 and ranks.value.notna().all()
    assert json.loads((out / 're_000.json').read_text()) == dict(method='re', alpha=1, z=2, min_effect=0.1, pi=0, up_w=0.75, cell_rule='fixed')
    selected = json.loads((out / 'selected_generator.json').read_text())
    assert selected['schema'] == 'vcc_2026_2_meta_v1' and selected['delta_array'] == 're'
    assert selected['settings']['alpha'] == 1 and selected['settings']['z'] == 2
    split = json.loads((out / 'control_split.json').read_text())
    assert split['gate_source'] == 'scoring real controls'
    assert split['generator_controls'] == split['scoring_controls'] == 60
    # H1 among the sources is a leak, not an exam.
    report = json.loads((tmp_path / 'pred/prediction_report.json').read_text()); report['source_backgrounds'].append('H1')
    (tmp_path / 'pred/prediction_report.json').write_text(json.dumps(report))
    result = subprocess.run(command[:-1] + ['fixed', '--out', str(tmp_path / 'leak')], capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and 'H1 leaked' in result.stderr


def test_h1_expression_gate_uses_scoring_half(scripts_v2, tmp_path):
    c, g = scripts_v2['common'], scripts_v2['generate']
    seed = 2026
    order = np.random.default_rng(c.seed_for(seed, 'H1', 'controls')).permutation(6)
    values = np.tile([10, 999990], (8, 1))
    values[order[:3]] = [0, 1000000]  # generation controls below 5 CPM
    values[order[3:]] = [6, 999994]   # scoring controls above 5 CPM
    labels = ['non-targeting']*6 + ['KD']*2
    ad.AnnData(sp.csr_matrix(values), obs=pd.DataFrame({'target_gene': labels}, index=[str(i) for i in range(8)]),
               var=pd.DataFrame(index=['G', 'KD'])).write_h5ad(tmp_path / 'h1.h5ad')
    opts = argparse.Namespace(cells=str(tmp_path / 'h1.h5ad'), obs_col='target_gene', targets=None,
                              min_cells=2, n_targets=0, mode='transfer', seed=seed)
    line, *_ = g.stream_h1(opts, ['G', 'KD'], tmp_path / 'real.h5ad')
    assert line.cpm[0] == 0
    assert line.gate[0]
    np.testing.assert_allclose(g.calibrate(np.array([[1., 0]]), line, ['KD'], g.Cfg(argparse.Namespace())), [[1, 0]])

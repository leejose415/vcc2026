"""Behavioral checks for scripts embedded in the new averaging notebook."""
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
    directory = tmp_path_factory.mktemp('averaging_scripts')
    for cell in notebook.cells:
        if cell.source.startswith('%%writefile'):
            header, body = cell.source.split('\n', 1)
            (directory / Path(header.split()[-1]).name).write_text(body)
    sys.path.insert(0, str(directory))
    return {name: importlib.import_module('vcc2_' + name)
            for name in ['common', 'build_effects', 'context_profile', 'transfer', 'generate']}


def test_notebook_structure_and_commands(scripts_v2):
    from IPython.core.inputtransformer2 import TransformerManager
    notebook = nbformat.read(ROOT / 'vcc_2026_2.ipynb', as_version=4)
    nbformat.validate(notebook)
    transform = TransformerManager()
    deferred = next(i for i, c in enumerate(notebook.cells) if c.source.startswith('### Deferred:'))
    assert notebook.cells[deferred + 1].cell_type == 'code'
    assert notebook.cells[deferred + 1].source == ''
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
    common = scripts_v2['common']; transfer = scripts_v2['transfer']
    assert not hasattr(common, 'shrink') and not hasattr(common.Effects, 'prior_for')
    for name in ['lobo_select', 'fit_source_scales', 'fit_ridge']:
        assert not hasattr(transfer, name)


def test_filter_retains_hits_at_full_size_and_excludes_unknown(scripts_v2):
    c = scripts_v2['common']
    # Strict >2 boundary, both signs, measured zeros, and unreliable SEs.
    result = c.significance_filter([1, -1, 0.5, -0.5, 0.49, 0, 1, 1, 1, 1, np.nan],
                                  [0.3, 0.3, 0.25, 0.25, 0.25, 0.25, 0, 1e-9, np.nan, -1, 0.2])
    np.testing.assert_allclose(result[:6], [1, -1, 0, 0, 0, 0])
    assert np.isnan(result[6:]).all()
    # Adding 18,532 null genes cannot change the amplitude of the one real hit.
    lfc = np.zeros(18533); lfc[-1] = 1
    filtered = c.significance_filter(lfc, np.full(18533, 0.3))
    assert filtered[-1] == 1 and np.count_nonzero(filtered) == 1


def make_contexts(scripts, tmp_path):
    c, t = scripts['common'], scripts['transfer']
    genes = np.array(['KD', 'G1', 'G2', 'G3'])
    pd.DataFrame({'gene_name': genes}).to_csv(tmp_path / 'genes.csv', index=False)
    specs = []
    for name, bg, values, errors, priority in [
        ('A', 'BG_A', [-1, 1, 1, 0.25], [0.25]*4, 0),
        ('A_extra', 'BG_A', [-1, 6, 6, 6], [0.25]*4, -1),
        ('B', 'BG_B', [-1, 3, 0, 1], [0.25, 0.25, np.nan, 0.25], 0),
        ('H1', 'H1', [-1, 7, 7, 7], [0.25]*4, 0),
    ]:
        d = tmp_path / name
        c.write_effects(d, genes, ['KD'], iter([(0, np.array(values), np.array(errors))]), {'legacy': False}, 1)
        c.save_profile(d / 'control_profile.npz', genes, dict(mean_cpm=np.array([0, 5, 10, 20.]),
            pooled_cpm=np.array([0, 8, 8, 19.]), log_mean=np.log1p([0, 8, 8, 19.]),
            detection=np.array([0, 0.1, 0.5, 0.8])))
        specs.append(dict(name=name, background=bg, source_priority=priority, platform='test',
                          effects=str(d), profile=str(d / 'control_profile.npz')))
    config = dict(schema=t.SCHEMA, genes=str(tmp_path / 'genes.csv'), contexts=specs)
    (tmp_path / 'config.json').write_text(json.dumps(config))
    pd.DataFrame({'target_gene': ['KD']}).to_csv(tmp_path / 'targets.csv', index=False)
    return genes, config


def test_average_is_equal_per_background_and_chunk_independent(scripts_v2, tmp_path):
    t = scripts_v2['transfer']; _, config = make_contexts(scripts_v2, tmp_path)
    # The excluded background is genuinely unreadable; preparing must not open it.
    config['contexts'][-1]['effects'] = '/nonexistent/heldout_effects'
    config['contexts'][-1]['profile'] = '/nonexistent/heldout_profile'
    (tmp_path / 'config.json').write_text(json.dumps(config))
    artifact = t.prepare(tmp_path / 'config.json', tmp_path / 'model', holdout='H1')
    assert all(c['background'] != 'H1' for c in artifact['contexts'])
    for chunk in [1, 4]:
        out = tmp_path / f'pred_{chunk}'
        t.infer(tmp_path / 'model', tmp_path / 'A/control_profile.npz', tmp_path / 'targets.csv', out, chunk=chunk)
        # G1: two independent backgrounds, not three tables. G2: B unknown, not zero.
        # G3: A measured non-hit is zero and B hit is one, hence 0.5 (no second filter).
        np.testing.assert_allclose(np.load(out / 'source_mean.npy'), [[-1, 2, 1, 0.5]])
        np.testing.assert_array_equal(np.load(out / 'source_count.npy'), [[2, 2, 1, 2]])
        np.testing.assert_array_equal(np.load(out / 'selected.npy'), np.load(out / 'source_mean.npy'))
        report = json.loads((out / 'prediction_report.json').read_text())
        assert report['heldout_background'] == 'H1'
        assert report['coverage'][0]['source_names'] == ['A', 'B']
        assert report['coverage'][0]['background_weights'] == {'BG_A': 0.5, 'BG_B': 0.5}
        assert report['weighting'] == 'equal_background'


def test_neighbour_donors_filtered_once_and_own_column_masked(scripts_v2, tmp_path):
    c, t = scripts_v2['common'], scripts_v2['transfer']
    genes = np.array(['N1', 'N2', 'G'])
    c.write_effects(tmp_path / 'effects', genes, ['N1', 'N2'],
        iter([(0, np.array([-4., 1, 1]), np.full(3, 0.25)),
              (1, np.array([2., -4, 0.25]), np.full(3, 0.25))]), {'legacy': False}, 2)
    pd.DataFrame({'gene_a': ['NEW', 'NEW'], 'gene_b': ['N1', 'N2'], 'score': [1, 1]}).to_csv(tmp_path / 'edges.csv', index=False)
    effects = t.FilteredEffects(c.Effects(tmp_path / 'effects', genes), t.FunctionalNeighbours(tmp_path / 'edges.csv'), genes)
    # At G the one retained effect averages with a measured zero to 0.5; it must
    # not be filtered a second time using uncertainty from the borrowed mixture.
    np.testing.assert_allclose(effects.get('NEW', np.arange(3)), [2, 1, 0.5])


def test_binary_gate_uses_per_cell_cpm_and_does_not_gate_by_knockdown(scripts_v2):
    c, g = scripts_v2['common'], scripts_v2['generate']
    np.testing.assert_array_equal(c.expression_gate([0, 4.9, 5, np.nextafter(5., 6.), 6, np.nan]),
                                  [False, False, False, True, True, False])
    line = g.Line(sp.csr_matrix([[1, 999], [100, 999900]]), None, ['G', 'KD'])
    assert line.cpm[0] == pytest.approx(550)  # mean per cell, not pooled ~101
    line.genes = ['KD', 'G1', 'G2', 'G3']
    line.gate = c.expression_gate([0, 5, 5.001, 100])
    cfg = g.Cfg(argparse.Namespace(alpha=1))
    delta = np.array([[-2, 3, -3, 2], [-1, -3, 3, -2.]])
    np.testing.assert_allclose(g.calibrate(delta, line, ['KD', 'KD'], cfg), [[0, 0, -3, 2], [0, 0, 3, -2]])


def test_unresolved_targets_require_explicit_fallback(scripts_v2, tmp_path):
    t = scripts_v2['transfer']; make_contexts(scripts_v2, tmp_path)
    t.prepare(tmp_path / 'config.json', tmp_path / 'model', holdout='H1')
    pd.DataFrame({'target_gene': ['MISSING']}).to_csv(tmp_path / 'missing.csv', index=False)
    with pytest.raises(ValueError, match='unsupported targets'):
        t.infer(tmp_path / 'model', tmp_path / 'A/control_profile.npz', tmp_path / 'missing.csv', tmp_path / 'strict')
    t.infer(tmp_path / 'model', tmp_path / 'A/control_profile.npz', tmp_path / 'missing.csv', tmp_path / 'fallback', allow_unsupported=True)
    np.testing.assert_array_equal(np.load(tmp_path / 'fallback/selected.npy'), np.zeros((1, 4)))
    report = json.loads((tmp_path / 'fallback/prediction_report.json').read_text())
    assert report['unresolved'] == ['MISSING']


def test_real_control_profile_preserves_baselines(scripts_v2, tmp_path):
    c, p = scripts_v2['common'], scripts_v2['context_profile']
    genes = ['G', 'KD']; counts = np.array([[1, 999], [100, 999900]])
    pd.DataFrame({'gene_name': genes}).to_csv(tmp_path / 'genes.csv', index=False)
    ad.AnnData(sp.csr_matrix(counts), obs=pd.DataFrame({'target_gene': ['non-targeting']*2}, index=['c1', 'c2']),
               var=pd.DataFrame(index=genes)).write_h5ad(tmp_path / 'controls.h5ad')
    p.profile(tmp_path / 'controls.h5ad', tmp_path / 'genes.csv', tmp_path / 'profile.npz', all_controls=True)
    result = c.load_profile(tmp_path / 'profile.npz', genes)
    np.testing.assert_allclose(result['mean_cpm'], [550, 999450])
    np.testing.assert_allclose(result['pooled_cpm'], counts.sum(0)/counts.sum()*1e6)


def test_count_generation_and_submission(scripts_v2, tmp_path):
    g, c = scripts_v2['generate'], scripts_v2['common']
    genes = ['KD', 'G']; pd.DataFrame({'gene_name': genes}).to_csv(tmp_path / 'gene_names.csv', index=False)
    pd.DataFrame({'target_gene': ['KD'], 'n_cells': [400]}).to_csv(tmp_path / 'pert_counts.csv', index=False)
    control = ad.AnnData(sp.csr_matrix([[30, 20]]*30), obs=pd.DataFrame({'target_gene': ['non-targeting']*30},
                         index=[str(i) for i in range(30)]), var=pd.DataFrame(index=genes))
    control.write_h5ad(tmp_path / 'context_A.h5ad')
    pred = tmp_path / 'pred'; pred.mkdir()
    for name, value in [('genes', genes), ('targets', ['KD']), ('source_mean', [[-1, 0.5]]), ('ontarget_lfc', [-1])]:
        np.save(pred / f'{name}.npy', value)
    command = [sys.executable, g.__file__, 'submit', '--data', str(tmp_path), '--lines', 'A',
               '--delta', f'A={pred}', '--delta-array', 'source_mean', '--out', str(tmp_path / 'output')]
    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    output = ad.read_h5ad(tmp_path / 'output/prediction.h5ad')
    assert output.shape == (400, 2)
    c.validate_counts(output.X)
    assert output.obs.context.unique().tolist() == ['A']


def test_official_h1_benchmark(scripts_v2, tmp_path):
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
    pred = tmp_path / 'pred'; pred.mkdir()
    # Independent synthetic source predictions: no perturbation outcomes used here.
    effects = np.zeros((3, len(genes)))
    for i in range(3): effects[i, 3+i::3] = 0.75; effects[i, i] = -1.5
    for name, value in [('genes', genes), ('targets', list(groups)), ('source_mean', effects), ('ontarget_lfc', [-1.5]*3)]:
        np.save(pred / f'{name}.npy', value)
    (pred / 'prediction_report.json').write_text(json.dumps(dict(schema='vcc_2026_2_average_v1',
        heldout_background='H1', selected_method='source_mean', source_backgrounds=['BG_A', 'BG_B'])))
    out = tmp_path / 'benchmark'
    command = [sys.executable, g.__file__, 'benchmark', '--cells', str(tmp_path / 'h1.h5ad'),
        '--genes', str(tmp_path / 'genes.csv'), '--targets', str(tmp_path / 'targets.csv'), '--predictions', str(pred),
        '--out', str(out), '--cell-eval', str(binary), '--num-threads', '1', '--n-cells', '100', '--min-cells', '50',
        '--alpha-grid', '1', '--pi-grid', '0', '--up-w-grid', '0.75', '--cell-rules', 'fixed']
    result = subprocess.run(command, capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (out / 'source_mean_000/agg_results.csv').exists()
    assert (out / 'ceiling/ceiling_agg.csv').exists()
    ranks = pd.read_csv(out / 'metric_ranks.csv')
    assert len(ranks) == 6 and ranks.value.notna().all()
    assert not (out / 'model_000.h5ad').exists()
    selected = json.loads((out / 'selected_generator.json').read_text())
    assert selected['schema'] == 'vcc_2026_2_average_v1'
    assert selected['settings']['alpha'] == 1
    split = json.loads((out / 'control_split.json').read_text())
    assert split['gate_source'] == 'scoring real controls'
    assert split['generator_controls'] == split['scoring_controls'] == 60


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

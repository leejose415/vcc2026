"""Regression tests execute scripts extracted from the notebook, which is the source of truth."""
import argparse
import importlib
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import pytest


@pytest.fixture(scope='session')
def scripts(tmp_path_factory):
    notebook = json.loads((Path(__file__).resolve().parents[1]/'vcc_2026.ipynb').read_text())
    directory = tmp_path_factory.mktemp('notebook_scripts')
    for cell in notebook['cells']:
        source = ''.join(cell['source'])
        if source.startswith('%%writefile'):
            header, body = source.split('\n', 1)
            (directory/header.rsplit('/', 1)[-1]).write_text(body)
    sys.path.insert(0, str(directory))
    return {name: importlib.import_module(name) for name in ['vcc_common', 'vcc_build_effects_v3',
            'vcc_context_profile_v3', 'vcc_transfer_model_v3', 'vcc_generate_v3', 'vcc_merge_h1']}


def test_profile_program_interface(scripts, tmp_path):
    c = scripts['vcc_common']; t = scripts['vcc_transfer_model_v3']; p = scripts['vcc_context_profile_v3']
    genes = np.array(['MCM5', 'PCNA', 'TYMS', 'RPL1'])
    pd.DataFrame({'gene_name': genes}).to_csv(tmp_path/'genes.csv', index=False)
    counts = np.array([[1, 2, 3, 4], [4, 4, 4, 8], [8, 8, 8, 16]])
    h = ad.AnnData(sp.csr_matrix(counts), obs=pd.DataFrame({'target_gene': ['non-targeting']*3}),
                   var=pd.DataFrame(index=genes))
    h.write_h5ad(tmp_path/'controls.h5ad')
    p.profile(tmp_path/'controls.h5ad', tmp_path/'genes.csv', tmp_path/'profile.npz', all_controls=True)
    profile = c.load_profile(tmp_path/'profile.npz', genes); selected = t.transfer_profile(profile)
    expected = counts.sum(0)/counts.sum()*1e6
    np.testing.assert_allclose(selected['baseline_cpm'], expected, rtol=1e-6)
    np.testing.assert_allclose(profile['bulk_cpm'], profile['pooled_cpm'])
    assert selected['baseline_program_levels'].dtype.kind == 'f'
    assert selected['baseline_program_levels'][0] > 0
    aligned = c.load_profile(tmp_path/'profile.npz', genes[::-1])
    np.testing.assert_allclose(aligned['pooled_cpm'], expected[::-1], rtol=1e-6)


def test_prior_and_zero_se(scripts):
    c = scripts['vcc_common']; t = scripts['vcc_transfer_model_v3']
    lfc = np.array([1., 1., 1., 1.]); se = np.array([1., 0., 1e-9, np.nan])
    np.testing.assert_allclose(c.shrink(lfc, se, prior_variance=1), [0.5]*4)
    contexts = [argparse.Namespace(background='A')]
    weights = t.balanced_label_weights(np.array([1., 0., np.nan, 1e-9]), np.zeros(4, int), contexts)
    np.testing.assert_allclose(weights, np.ones(4))
    assert c.estimate_prior(np.array([-2., 2.]), np.array([1., 1.])) == pytest.approx(7.)


def test_projection_and_low_baseline_induction(scripts):
    g = scripts['vcc_generate_v3']
    share = np.array([0.2, 0.8, 0.0]); delta = np.array([-1., 0., 4.])
    factors, composition, normalizer = g.project_factors(share, delta)
    assert composition.sum() == pytest.approx(1)
    np.testing.assert_allclose(factors[:2], np.exp2(delta[:2])/normalizer)
    assert composition[2] > 0
    line = g.Line(sp.csr_matrix(np.tile([4000, 16000, 0], (100, 1))), None, ['A', 'B', 'ATF3'])
    cfg = g.Cfg(argparse.Namespace(n_cells=100, up_w=0.75))
    calibrated = g.calibrate(np.array([[-1., 0, 4.]]), line, ['other'], cfg)
    assert calibrated[0, 2] == 4
    counts, _, _ = g.generate_one(line, calibrated[0], None, 0.3, cfg)
    assert counts[:, 2].sum() > 0
    scripts['vcc_common'].validate_counts(counts)


def test_low_expression_down_gate(scripts):
    g = scripts['vcc_generate_v3']
    line = g.Line(sp.csr_matrix([[100, 0]]*20), None, ['G', 'ATF3'])
    cfg = g.Cfg(argparse.Namespace())
    result = g.calibrate(np.array([[0., -3.], [0., 3.]]), line, ['X', 'Y'], cfg)
    assert result[0, 1] == 0 and result[1, 1] == 3


def test_generator_order_and_common_draws(scripts):
    g = scripts['vcc_generate_v3']; rng = np.random.default_rng(3)
    line = g.Line(sp.csr_matrix(rng.poisson([30, 40, 20], (40, 3))), None, ['A', 'B', 'C'])
    cfg = g.Cfg(argparse.Namespace(n_cells=20, pi=0.15))
    effects = {'first': np.array([0.4, -0.3, 0.]), 'second': np.array([-0.2, 0.2, 0.1])}
    first = {t: g.generate_one(line, d, None, 0.3, cfg, line_key='A', target_key=t)[0] for t, d in effects.items()}
    second = {t: g.generate_one(line, effects[t], None, 0.3, cfg, line_key='A', target_key=t)[0] for t in reversed(effects)}
    for target in first: assert (first[target] != second[target]).nnz == 0
    x = g.generate_one(line, effects['first'], None, 0.3, cfg, line_key='A', target_key='different')[0]
    assert (first['first'] != x).nnz == 0


def test_escaper_floor_and_metadata(scripts, tmp_path):
    g = scripts['vcc_generate_v3']
    line = g.Line(sp.csr_matrix([[100, 100]]*50), None, ['A', 'B'])
    cfg = g.Cfg(argparse.Namespace(n_cells=40, pi=0.5))
    _, diagnostics, _ = g.generate_one(line, np.zeros(2), 0, 0.02, cfg)
    assert diagnostics['effective_pi'] < 0.05
    pd.DataFrame({'target_gene': ['A', 'B'], 'n_cells': [400, 400], 'median_umi_per_cell': [1000, 2000]}).to_csv(tmp_path/'metadata.csv', index=False)
    result = g.target_metadata(tmp_path/'metadata.csv', 'A', 400)
    assert result['B'] == {'n_cells': 400, 'median_umi_per_cell': 2000}
    assert g.target_metadata(tmp_path/'metadata.csv', 'A', 30, False)['B']['n_cells'] == 30


def test_every_output_row_validated(scripts, tmp_path):
    g = scripts['vcc_generate_v3']
    values = np.ones((2200, 2)); values[2199, 1] = 1.25
    h = ad.AnnData(sp.csr_matrix(values)); h.write_h5ad(tmp_path/'bad.h5ad')
    with pytest.raises(ValueError, match='integer raw counts'):
        g.check_output(tmp_path/'bad.h5ad', (2200, 2))


def test_failed_guides_filtered_before_selection(scripts, tmp_path):
    builder = scripts['vcc_build_effects_v3']; c = scripts['vcc_common']
    genes = np.array(['MCM5', 'PCNA', 'TYMS', 'RPL1'])
    pd.DataFrame({'gene_name': genes}).to_csv(tmp_path/'genes.csv', index=False)
    blocks = [np.tile([10, 10, 10, 10], (20, 1)), np.tile([10, 10, 10, 10], (60, 1)), np.tile([1, 10, 10, 10], (30, 1))]
    labels = ['non-targeting']*20+['MCM5']*90
    guides = ['ntc']*20+['failed']*60+['good']*30
    h = ad.AnnData(sp.csr_matrix(np.vstack(blocks)), obs=pd.DataFrame({'target_gene': labels, 'guide': guides}), var=pd.DataFrame(index=genes))
    h.write_h5ad(tmp_path/'raw.h5ad')
    opts = argparse.Namespace(genes=tmp_path/'genes.csv', cells=tmp_path/'raw.h5ad', out=tmp_path/'effects', gene_col=None,
             obs_col='target_gene', control='non-targeting', layer=None, target_split=None, guide_col='guide',
             name='replogle_test', background='K562', platform='3prime', kind='replogle', priority=0,
             min_cells=20, min_ontarget_lfc=-0.1, block=16, pseudocount=0.5)
    builder.cmd_cells(opts)
    qc = pd.read_csv(tmp_path/'effects/guide_qc.csv')
    assert not qc.loc[qc.guide == 'failed', 'passed'].iloc[0]
    assert qc.loc[qc.guide == 'good', 'passed'].iloc[0]
    table = c.Effects(tmp_path/'effects', genes)
    observed = table.get('MCM5', np.arange(4))[0]
    expected = np.log2((np.array([1, 10, 10, 10])/31*1e6+0.5)/(250000+0.5))
    np.testing.assert_allclose(observed, expected, atol=0.002)
    assert table.meta['label_estimator'] == 'per_cell_cpm' and table.meta['label_pseudocount'] == 0.5


def make_contexts(scripts, tmp_path, backgrounds=('CD4_T', 'HCT116', 'HEK293T', 'H1')):
    c = scripts['vcc_common']; rng = np.random.default_rng(14)
    genes = np.array(['K0', 'K1', 'K2', 'K3', 'K4', 'MCM5', 'PCNA', 'TYMS', 'ATF3', 'RPL1'])
    pd.DataFrame({'gene_name': genes}).to_csv(tmp_path/'genes.csv', index=False)
    contexts = []
    for index, bg in enumerate(backgrounds):
        directory = tmp_path/bg
        effect = rng.normal(0, 0.5, (5, len(genes)))
        for i in range(5): effect[i, i] = -1.5
        error = np.full_like(effect, 0.25)
        c.write_effects(directory, genes, genes[:5], ((i, effect[i], error[i]) for i in range(5)), {'legacy': False}, 5)
        pooled = rng.uniform(50, 100, len(genes)); pooled = pooled/pooled.sum()*1e6
        c.save_profile(directory/'control_profile.npz', genes, dict(mean_cpm=pooled, pooled_cpm=pooled,
                       log_mean=np.log1p(pooled), detection=np.full(len(genes), 0.8), legacy_baseline=np.bool_(False)))
        contexts.append(dict(name=bg, background=bg, platform='flex', effects=str(directory), profile=str(directory/'control_profile.npz')))
    cfg = dict(genes=str(tmp_path/'genes.csv'), contexts=contexts, ridge_alphas=[10, 100], transfer_grid=[dict(temperature=0.5, gate_center=1.8)],
               selection_backgrounds=list(backgrounds), effect_cap=8)
    (tmp_path/'config.json').write_text(json.dumps(cfg))
    return genes, cfg


def test_nested_lobo_and_refit_exclude_background(scripts, tmp_path):
    t = scripts['vcc_transfer_model_v3']; genes, cfg = make_contexts(scripts, tmp_path)
    first = t.train(tmp_path/'config.json', tmp_path/'model', holdout='H1', max_kd=0, per_kd=8)
    assert 'H1' not in [c['background'] for c in first['contexts']]
    assert set(first['train_targets']) == set(genes[:5])
    report = json.loads((tmp_path/'model/training_report.json').read_text())
    assert {fold['background'] for fold in report['selection']['history']} == {'CD4_T', 'HCT116', 'HEK293T'}
    assert report['selection']['original_folds'] == 3
    # Make the held-out files unreadable/unavailable: training must still be identical.
    for spec in cfg['contexts']:
        if spec['background'] == 'H1':
            spec['effects'] = '/nonexistent/heldout_effects'
            spec['profile'] = '/nonexistent/heldout_profile'
    (tmp_path/'config2.json').write_text(json.dumps(cfg))
    second = t.train(tmp_path/'config2.json', tmp_path/'model2', holdout='H1', max_kd=0, per_kd=8)
    np.testing.assert_allclose(first['model']['estimator'].coef_, second['model']['estimator'].coef_)
    assert first['scales'] == second['scales']


def test_neighbours_and_unresolved_targets(scripts, tmp_path):
    t = scripts['vcc_transfer_model_v3']; c = scripts['vcc_common']; genes, cfg = make_contexts(scripts, tmp_path)
    path = tmp_path/'neighbours.csv'
    pd.DataFrame({'gene_a': ['NEW', 'NEW'], 'gene_b': ['K0', 'K1'], 'score': [3, 1]}).to_csv(path, index=False)
    graph = t.FunctionalNeighbours(path, 2)
    original = c.Effects(tmp_path/'CD4_T', genes); borrowed = t.BorrowedEffects(original, graph, genes)
    result, errors = borrowed.get('NEW', np.arange(len(genes)))
    k0, k1 = (original.get(g, np.arange(len(genes)))[0] for g in ['K0', 'K1'])
    expected = k0*0.75+k1*0.25
    # Each neighbour's own knockdown is masked; the other neighbour fills that gene in.
    expected[0], expected[1] = k1[0], k0[1]
    np.testing.assert_allclose(result, expected, rtol=1e-6)
    assert result[0] > -1 and result[1] > -1      # the -1.5 on-target drops did not leak
    assert np.isfinite(errors).all()
    # Single-neighbour borrowing leaves that neighbour's own column unknown, not -1.5.
    single = t.BorrowedEffects(original, t.FunctionalNeighbours(path, 1), genes).get('NEW', np.arange(len(genes)))[0]
    assert np.isnan(single[0]) and np.isfinite(single[1:]).all()
    cfg['neighbours'] = str(path); (tmp_path/'config.json').write_text(json.dumps(cfg))
    t.train(tmp_path/'config.json', tmp_path/'model', max_kd=0, per_kd=8)
    pd.DataFrame({'target_gene': ['NEW']}).to_csv(tmp_path/'targets.csv', index=False)
    t.infer(tmp_path/'model', tmp_path/'CD4_T/control_profile.npz', tmp_path/'targets.csv', tmp_path/'prediction')
    report = json.loads((tmp_path/'prediction/prediction_report.json').read_text())
    assert report['coverage'][0]['borrowed']
    assert np.load(tmp_path/'prediction/source_count.npy').min() > 0
    pd.DataFrame({'target_gene': ['UNSUPPORTED']}).to_csv(tmp_path/'targets2.csv', index=False)
    with pytest.raises(ValueError, match='unsupported targets'):
        t.infer(tmp_path/'model', tmp_path/'CD4_T/control_profile.npz', tmp_path/'targets2.csv', tmp_path/'bad_prediction')


def test_gate_off_and_similarity(scripts, tmp_path):
    t = scripts['vcc_transfer_model_v3']; genes, cfg = make_contexts(scripts, tmp_path)
    contexts = t.load_contexts(cfg, genes); profile = dict(contexts[0].profile)
    profile['baseline_cpm'] = profile['baseline_cpm'].copy(); profile['baseline_log'] = profile['baseline_log'].copy()
    profile['baseline_cpm'][0] = profile['baseline_log'][0] = 0
    generics = t.generic_signatures(contexts[1:], genes[:5], genes, 8)
    f = t.features(profile, 'flex', contexts[1:], 'K0', np.arange(len(genes)), genes, generics)
    np.testing.assert_allclose(f['sv'], 0)
    np.testing.assert_allclose(f['weights'].sum(1), 1)
    assert t.centred_similarity(np.arange(5.), np.arange(5.)) == pytest.approx(1)


def test_official_aggregate_ranks(scripts, tmp_path):
    g = scripts['vcc_generate_v3']; metrics = g.VCC2026_METRICS
    paths = {}
    for name, benefit in [('good', 1), ('bad', 0)]:
        directory = tmp_path/name; directory.mkdir(); paths[name] = str(directory)
        row = {'statistic': 'mean'}
        row.update({m: (benefit if d == 'max' else 1-benefit) for m, d in metrics.items()})
        pd.DataFrame([row]).to_csv(directory/'agg_results.csv', index=False)
    directions = tmp_path/'directions.csv'
    pd.DataFrame({'metric': list(metrics), 'direction': list(metrics.values())}).to_csv(directions, index=False)
    frame, ranks = g.metric_rank_table(paths, directions)
    assert len(frame) == 12 and ranks.index[0] == 'good'
    assert ranks['good'] == 1 and ranks['bad'] == 2


def test_submit_metadata_counts_and_validation(scripts, tmp_path):
    g = scripts['vcc_generate_v3']
    genes = ['G', 'ATF3']; pd.DataFrame({'gene_name': genes}).to_csv(tmp_path/'gene_names.csv', index=False)
    pd.DataFrame({'target_gene': ['G'], 'n_cells': [8], 'median_umi_per_cell': [100]}).to_csv(tmp_path/'pert_counts.csv', index=False)
    controls = ad.AnnData(sp.csr_matrix([[30, 20]]*20), obs=pd.DataFrame({'target_gene': ['non-targeting']*20}), var=pd.DataFrame(index=genes))
    controls.write_h5ad(tmp_path/'context_A.h5ad')
    pred = tmp_path/'pred'; pred.mkdir(); np.save(pred/'genes.npy', genes); np.save(pred/'targets.npy', ['G'])
    np.save(pred/'selected.npy', [[-1., 0.2]]); np.save(pred/'ontarget_lfc.npy', [-1.])
    opts = argparse.Namespace(data=str(tmp_path), lines=['A'], delta=[f'A={pred}'], delta_array='selected', out=str(tmp_path/'result'),
               checks=1, settings=None, kd_orion=None, kd_tcell=None, kd_tcell_condition='Stim48hr', required_n_cells=400)
    with pytest.raises(ValueError, match='exactly 400'): g.cmd_submit(opts)
    opts.required_n_cells = 8
    g.cmd_submit(opts)
    h = ad.read_h5ad(tmp_path/'result/prediction.h5ad')
    assert h.shape == (8, 2)
    assert set(h.obs.context.astype(str)) == {'A'}
    assert json.loads((tmp_path/'result/generator_report.json').read_text())['integer_counts']


def test_notebook_and_cli_syntax(scripts):
    import ast
    import nbformat
    from IPython.core.inputtransformer2 import TransformerManager
    path = Path(__file__).resolve().parents[1]/'vcc_2026.ipynb'
    notebook = nbformat.read(path, as_version=4); nbformat.validate(notebook)
    transform = TransformerManager()
    for cell in notebook.cells:
        if cell.cell_type == 'code': ast.parse(transform.transform_cell(cell.source))


def test_sums_guide_qc_and_baseline_compatibility(scripts, tmp_path):
    builder = scripts['vcc_build_effects_v3']; c = scripts['vcc_common']; t = scripts['vcc_transfer_model_v3']
    genes = np.array(['MCM5', 'PCNA', 'TYMS', 'RPL1']); sizes = np.array([20, 60, 30])
    rows = np.array([[10., 10, 10, 10], [10., 10, 10, 10], [1., 10, 10, 10]])
    cpm = rows/rows.sum(1, keepdims=True)*1e6
    obs = pd.DataFrame({'target_gene': ['non-targeting', 'MCM5', 'MCM5'], 'n_cells': sizes}, index=['ctrl', 'bad', 'good'])
    h = ad.AnnData(rows, obs=obs, var=pd.DataFrame(index=genes),
        layers={'sum_counts': rows*sizes[:, None], 'sum_cpm': cpm*sizes[:, None],
                'sum_cpm2': cpm**2*sizes[:, None], 'n_expr': np.broadcast_to(sizes[:, None], rows.shape).copy()})
    h.write_h5ad(tmp_path/'sums.h5ad'); pd.DataFrame({'gene_name': genes}).to_csv(tmp_path/'genes.csv', index=False)
    opts = argparse.Namespace(genes=tmp_path/'genes.csv', pseudobulk=tmp_path/'sums.h5ad', out=tmp_path/'effects', control=None,
            min_cells=20, chunk=2, lfc_estimator='per_cell', min_ontarget_lfc=-0.1,
            name='Orion', background='HCT116', platform='3prime', priority=0)
    builder.cmd_sums(opts, 'orion')
    table = c.Effects(tmp_path/'effects', genes)
    expected = np.log2((cpm[2]+0.5)/(cpm[0]+0.5))
    np.testing.assert_allclose(table.get('MCM5', np.arange(4))[0], expected, atol=0.002)
    profile = t.transfer_profile(c.load_profile(tmp_path/'effects/control_profile.npz', genes))
    np.testing.assert_allclose(profile['baseline_cpm'], cpm[0])
    qc = pd.read_csv(tmp_path/'effects/guide_qc.csv')
    assert qc.passed.tolist() == [False, True]
    assert table.meta['estimator_matches_official'] is False


def test_shell_cells_and_pipeline_commands_parse():
    import ast
    import subprocess
    nb = json.loads((Path(__file__).resolve().parents[1]/'vcc_2026.ipynb').read_text())
    for cell in nb['cells']:
        source = ''.join(cell['source'])
        if cell['cell_type'] != 'code': continue
        if source.startswith('%%bash'):
            command = source.split('\n', 1)[1]
            assert '%%bash' not in command
            subprocess.run(['bash', '-n'], input=command, text=True, check=True)
        elif not source.startswith('%%'):
            for node in ast.walk(ast.parse(source)):
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in ['run', 'bg']:
                    arg = node.args[-1]
                    if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                        subprocess.run(['bash', '-n'], input=arg.value, text=True, check=True)


def test_official_cli_invocation(scripts, tmp_path, monkeypatch):
    g = scripts['vcc_generate_v3']; calls = []
    monkeypatch.setattr(g.shutil, 'which', lambda binary: '/tmp/cell-eval2')
    def run(command, **kwargs):
        calls.append(command)
        return argparse.Namespace(stdout='--ceiling --preset --profile --set --control --input-type')
    monkeypatch.setattr(g.subprocess, 'run', run)
    opts = argparse.Namespace(cell_eval='cell-eval2', eval_profile='vcc2026', num_threads=4)
    g.run_official(tmp_path/'pred.h5ad', tmp_path/'real.h5ad', tmp_path/'metrics', opts)
    command = calls[-1]
    assert command[command.index('--preset')+1] == 'vcc2026'
    assert command[command.index('--profile')+1] == 'vcc2026'
    assert 'de.backend=pdex' in command and 'num_threads=4' in command
    assert '--ceiling' in command and '--control-pert' not in command
    assert (tmp_path/'metrics/invocation.json').exists()


def test_actual_official_h1_benchmark(scripts, tmp_path):
    """Optional CPU integration test: actual official metrics, not mocked proxies."""
    import subprocess
    binary = Path(sys.executable).parent/'cell-eval2'
    if not binary.exists(): pytest.skip('install requirements-eval.txt for the official evaluator integration test')
    genes = np.array(['K0', 'K1', 'K2']+[f'G{i}' for i in range(37)])
    rng = np.random.default_rng(2); baseline = rng.integers(30, 100, len(genes))
    controls = rng.poisson(baseline, (120, len(genes))); groups = {}
    for i, kd in enumerate(genes[:3]):
        factors = np.ones(len(genes)); factors[3+i::3] = 3; factors[i] = 0.2
        groups[kd] = rng.poisson(baseline*factors, (100, len(genes)))
    labels = ['non-targeting']*len(controls)+[t for t, x in groups.items() for _ in range(len(x))]
    h = ad.AnnData(sp.csr_matrix(np.vstack([controls]+list(groups.values()))),
        obs=pd.DataFrame({'target_gene': labels}, index=[str(i) for i in range(len(labels))]), var=pd.DataFrame(index=genes))
    h.write_h5ad(tmp_path/'h1.h5ad'); pd.DataFrame({'gene_name': genes}).to_csv(tmp_path/'genes.csv', index=False)
    pd.DataFrame({'target_gene': list(groups)}).to_csv(tmp_path/'targets.csv', index=False)
    def mean_cpm(x): return (x/np.maximum(x.sum(1, keepdims=True), 1)*1e6).mean(0)
    effects = np.stack([np.log2((mean_cpm(x)+1e-9)/(mean_cpm(controls)+1e-9)) for x in groups.values()])
    pred = tmp_path/'pred'; pred.mkdir(); np.save(pred/'genes.npy', genes); np.save(pred/'targets.npy', list(groups))
    np.save(pred/'model.npy', effects*0.8); np.save(pred/'source_mean.npy', effects*0.6)
    np.save(pred/'ontarget_lfc.npy', [effects[i, i] for i in range(3)])
    (pred/'prediction_report.json').write_text(json.dumps({'heldout_background': 'H1', 'selected_method': 'model'}))
    output = tmp_path/'benchmark'
    command = [sys.executable, scripts['vcc_generate_v3'].__file__, 'benchmark',
        '--cells', str(tmp_path/'h1.h5ad'), '--genes', str(tmp_path/'genes.csv'), '--targets', str(tmp_path/'targets.csv'),
        '--predictions', str(pred), '--delta-array', 'selected', '--out', str(output), '--cell-eval', str(binary),
        '--num-threads', '1', '--n-cells', '100', '--min-cells', '50', '--alpha-grid', '1', '--pi-grid', '0',
        '--up-w-grid', '0.75', '--cell-rules', 'fixed']
    result = subprocess.run(command, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout+result.stderr
    for method in ['model', 'source_mean']:
        assert (output/f'{method}_000/agg_results.csv').exists()
        assert (output/f'{method}_000/ceiling_agg.csv').exists()
    ranks = pd.read_csv(output/'metric_ranks.csv')
    assert len(ranks) == 12 and ranks.value.notna().all()
    selected = json.loads((output/'selected_generator.json').read_text())
    assert selected['diagnostic_best_method'] == 'model'
    assert selected['settings']['n_cells'] == 100


def test_low_count_labels_stay_bounded(scripts, tmp_path):
    """Zero-count groups on low-count genes used to give LFC near -30 with SE ~0.5."""
    builder = scripts['vcc_build_effects_v3']; c = scripts['vcc_common']; rng = np.random.default_rng(0)
    genes = np.array(['KD']+[f'G{i}' for i in range(1, 80)])
    pd.DataFrame({'gene_name': genes}).to_csv(tmp_path/'genes.csv', index=False)
    # Shallow libraries (~1.2k counts) on purpose: one count is ~800 CPM, so 20 genes at
    # ~2.5 CPM are mostly all-zero in 100 cells. Only the SE floor keeps those labels honest.
    mu = np.r_[np.full(60, 20.0), np.full(20, 0.003)]
    controls = rng.poisson(mu, (3000, len(genes))); perturbed = rng.poisson(mu, (100, len(genes)))
    perturbed[:, 0] = rng.poisson(4, 100); perturbed[:, 1] = rng.poisson(60, 100)    # knockdown; real induction
    h = ad.AnnData(sp.csr_matrix(np.vstack([controls, perturbed])), var=pd.DataFrame(index=genes),
                   obs=pd.DataFrame({'target_gene': ['non-targeting']*3000+['KD']*100}, index=[str(i) for i in range(3100)]))
    h.write_h5ad(tmp_path/'raw.h5ad')
    opts = argparse.Namespace(genes=tmp_path/'genes.csv', cells=tmp_path/'raw.h5ad', out=tmp_path/'effects', gene_col=None,
             obs_col='target_gene', control='non-targeting', layer=None, target_split=None, guide_col=None, name='low',
             background='K562', platform='3prime', kind='replogle', priority=0, min_cells=20, min_ontarget_lfc=-0.1,
             block=512, pseudocount=0.5)
    builder.cmd_cells(opts)
    table = c.Effects(tmp_path/'effects', genes); lfc, se = table.get('KD', np.arange(len(genes)))
    assert np.abs(lfc[60:]).max() < 5                        # bounded by the pseudocount (was ~-30)
    assert se[60:].min() > 1                                 # all-zero groups are uncertain, not precise
    shrunk = c.shrink(lfc, se, prior_variance=table.prior_variance)
    assert np.abs(shrunk[60:]).max() < 0.5                   # ... so shrinkage pulls them toward 0
    assert shrunk[0] < -1 and shrunk[1] > 1                  # the knockdown and genuine induction survive
    assert table.prior_variance < 2                          # shrinkage prior is no longer inflated (was ~130)


def test_holdout_panels_only_covered_targets(scripts, tmp_path):
    t = scripts['vcc_transfer_model_v3']; c = scripts['vcc_common']; genes, cfg = make_contexts(scripts, tmp_path)
    rng = np.random.default_rng(5); effect = rng.normal(0, 0.5, (6, len(genes)))
    c.write_effects(tmp_path/'H1', genes, list(genes[:5])+['ATF3'], ((i, effect[i], np.full(len(genes), 0.25)) for i in range(6)),
                    {'legacy': False}, 6)                                        # ATF3 only measured in H1
    summary = t.holdout_panels(tmp_path/'config.json', tmp_path/'panels', max_targets=0)
    assert summary['H1']['retained_targets'] == 6 and summary['H1']['panel'] == 5
    assert 'ATF3' not in pd.read_csv(tmp_path/'panels/H1.csv').target_gene.tolist()
    t.train(tmp_path/'config.json', tmp_path/'model', holdout='H1', max_kd=0, per_kd=8)
    t.infer(tmp_path/'model', tmp_path/'H1/control_profile.npz', tmp_path/'panels/H1.csv', tmp_path/'pred')   # no unsupported-target error
    capped = t.holdout_panels(tmp_path/'config.json', tmp_path/'capped', max_targets=3)
    again = t.holdout_panels(tmp_path/'config.json', tmp_path/'capped2', max_targets=3)
    assert capped['H1']['panel'] == 3
    assert pd.read_csv(tmp_path/'capped/H1.csv').equals(pd.read_csv(tmp_path/'capped2/H1.csv'))


def test_merge_h1_deduplicates_controls(scripts, tmp_path):
    m = scripts['vcc_merge_h1']; c = scripts['vcc_common']; rng = np.random.default_rng(1)
    genes = [f'G{i}' for i in range(6)]; controls = rng.poisson(5, (7, 6)); paths, expected = [], [controls]
    for k, target in enumerate(['A', 'B', 'C']):
        perturbed = rng.poisson(5, (3+k, 6)); expected.append(perturbed)
        obs = pd.DataFrame({'target_gene': ['non-targeting']*7+[target]*len(perturbed), 'batch': f'b{k}'},
                           index=[f'ctrl{i}' for i in range(7)]+[f'{target}{i}' for i in range(len(perturbed))])
        path = tmp_path/f'split{k}.h5ad'; paths.append(str(path))
        ad.AnnData(sp.csr_matrix(np.vstack([controls, perturbed]).astype(np.float32)), obs=obs, var=pd.DataFrame(index=genes)).write_h5ad(path)
    report = m.merge(paths, tmp_path/'H1_full.h5ad', block=2)                 # small blocks exercise append
    merged = ad.read_h5ad(tmp_path/'H1_full.h5ad')
    assert report['n_controls'] == 7 and merged.n_obs == 7+3+4+5 and report['n_targets'] == 3
    assert (merged.obs.target_gene == 'non-targeting').sum() == 7 and merged.obs_names.is_unique
    np.testing.assert_array_equal(merged.X.toarray(), np.vstack(expected))
    assert list(merged.var_names) == genes and all(x['controls_match_first'] for x in report['inputs'])
    c.validate_counts(merged.X)
    with pytest.raises(FileExistsError): m.merge(paths, tmp_path/'H1_full.h5ad')

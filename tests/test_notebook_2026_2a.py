"""Behavioral checks for vcc_2026_2a: leaderboard-scored tuning, top-N calling, detection floor, diagnostics."""
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
MODULES = ['vcc2_common', 'vcc2_meta', 'vcc2_generate', 'vcc2_adhoc']


@pytest.fixture(scope='module')
def v2a(tmp_path_factory):
    notebook = nbformat.read(ROOT / 'vcc_2026_2a.ipynb', as_version=4)
    directory = tmp_path_factory.mktemp('scripts_2a')
    for cell in notebook.cells:
        if cell.source.startswith('%%writefile'):
            header, body = cell.source.split('\n', 1)
            (directory / Path(header.split()[-1]).name).write_text(body)
    for name in MODULES: sys.modules.pop(name, None)
    sys.path.insert(0, str(directory))
    try:
        yield {name.split('_')[1]: importlib.import_module(name) for name in MODULES}
    finally:
        sys.path.remove(str(directory))
        for name in MODULES: sys.modules.pop(name, None)


def test_notebook_structure_and_switches(v2a):
    from IPython.core.inputtransformer2 import TransformerManager
    notebook = nbformat.read(ROOT / 'vcc_2026_2a.ipynb', as_version=4)
    nbformat.validate(notebook)
    transform = TransformerManager()
    written = [c.source.split('\n', 1)[0] for c in notebook.cells if c.source.startswith('%%writefile')]
    assert written == [f'%%writefile ~/vcc/code/vcc_2026_2a/{n}.py' for n in MODULES]
    adhoc = next(i for i, c in enumerate(notebook.cells) if c.source.startswith('## Ad hoc tests'))
    for i, cell in enumerate(notebook.cells):
        if cell.cell_type != 'code': continue
        if i < adhoc: assert cell.execution_count is None and cell.outputs == []   # only the ad hoc checks carry outputs
        assert "'vcc_2026_2'" not in cell.source and 'code/vcc_2026_2/' not in cell.source   # never writes into vcc_2026_2
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
    switches = [l for c in notebook.cells for l in c.source.splitlines() if l.startswith(('SUBMIT_TO_LEADERBOARD =', 'BACK_UP_TO_BUCKET =', 'DOWNLOAD_INPUTS ='))]
    assert sorted(switches) == ['BACK_UP_TO_BUCKET = False', 'DOWNLOAD_INPUTS = False', 'SUBMIT_TO_LEADERBOARD = False']
    assert any("RUN_NAME = 'vcc_2026_2a_topn'" in c.source for c in notebook.cells)
    for name in ['meta', 'generate']:
        result = subprocess.run([sys.executable, v2a[name].__file__, '--help'], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
    assert v2a['meta'].SCHEMA == v2a['generate'].SCHEMA == 'vcc_2026_2_meta_v2'


# ---------------------------------------------------------------- edit 3
def write_panel(root, name, background, targets, genes, lfc, se):
    obs = pd.DataFrame(dict(panel='ours', status='kept', n_cells=200.0, ontarget_lfc=-2.0, padj=0.01), index=pd.Index(targets, name='target'))
    (root / name).mkdir(parents=True, exist_ok=True)
    ad.AnnData(X=None, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name='gene')),
               layers=dict(lfc=np.asarray(lfc, np.float32), se=np.asarray(se, np.float32)), uns=dict(background=background)).write_h5ad(root / name / 'effects_panel.h5ad')


def test_meta_writes_calling_inputs(v2a, tmp_path):
    m = v2a['meta']; nan = np.nan
    genes = ['KD', 'G1', 'G2', 'G3', 'G4']; targets = ['KD']
    pd.DataFrame({'gene_name': genes}).to_csv(tmp_path / 'genes.csv', index=False)
    pd.DataFrame({'target_gene': targets}).to_csv(tmp_path / 'targets.csv', index=False)
    # G1: three backgrounds agree; G2: one disagrees; G3: one background; G4: none.
    write_panel(tmp_path / 'src', 'HCT116', 'HCT116', targets, genes, [[-2, 1.0, 1.0, 0.5, nan]], [[0.1, 0.1, 0.1, 0.2, nan]])
    write_panel(tmp_path / 'src', 'HEK293T', 'HEK293T', targets, genes, [[-2, 0.8, 0.9, nan, nan]], [[0.1, 0.1, 0.1, nan, nan]])
    write_panel(tmp_path / 'src', 'K562_gwps', 'K562', targets, genes, [[-2, 1.2, -0.3, nan, nan]], [[0.1, 0.1, 0.1, nan, nan]])
    m.fit(tmp_path / 'src', tmp_path / 'genes.csv', tmp_path / 'model', None, ['HCT116', 'HEK293T', 'K562_gwps'])
    m.predict(tmp_path / 'model', tmp_path / 'targets.csv', ['A'], tmp_path / 'pred')
    d = tmp_path / 'pred/A'; re, se, z, n_bg, agree = (np.load(d / f'{x}.npy') for x in ['re', 're_se', 're_z', 'n_bg', 'sign_agree'])
    for x in [z, n_bg, agree]: assert x.shape == re.shape
    np.testing.assert_allclose(z[0, 1:4], np.abs(re[0, 1:4]) / se[0, 1:4], rtol=1e-6)
    np.testing.assert_array_equal(n_bg[0, 1:4], [3, 3, 1])
    np.testing.assert_allclose(agree[0, 1:4], [1, 2 / 3, 1])
    assert np.isnan(z[0, [0, 4]]).all() and np.isnan(n_bg[0, [0, 4]]).all() and np.isnan(agree[0, [0, 4]]).all()   # own gene, G4


# ---------------------------------------------------------------- edits 4 and 5
def line_with(g, c, genes, cpm_cells):
    """A Line built from integer counts proportional to the given per-cell values."""
    X = sp.csr_matrix(np.rint(np.asarray(cpm_cells, float)))
    return g.Line(X, None, genes)


def test_top_n_calling(v2a):
    g, c = v2a['generate'], v2a['common']
    genes = ['KD'] + [f'G{i}' for i in range(9)]
    line = line_with(g, c, genes, np.full((4, 10), 1e5))
    line.gate = c.expression_gate([100] * 8 + [4, 100])             # G7 sits at 4 CPM: never called
    D = np.array([[-3.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 5.0, 0.2],
                  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1]])
    info = dict(z=np.abs(D) * 10, n_bg=np.full(D.shape, 3.0), sign_agree=np.ones(D.shape))
    info['sign_agree'][0, 1] = 2 / 3                                    # G0: backgrounds disagree in sign
    cfg = g.Cfg(argparse.Namespace(call='topn', top_n=3, sign_rule=True))
    out = g.calibrate(D, line, ['KD', 'G0'], cfg, info=info)
    assert (out[0] != 0).sum() == 3 and (out[1] != 0).sum() == 1      # exactly N, or fewer when fewer candidates
    np.testing.assert_allclose(out[0], [0, 0, 0.8, 0.7, 0.6, 0, 0, 0, 0, 0])   # KD own, G0 sign, G7 gate
    assert out[1, 1] == 0                                                # target G0's own column
    loose = g.calibrate(D, line, ['KD', 'G0'], g.Cfg(argparse.Namespace(call='topn', top_n=3, sign_rule=False)), info=info)
    np.testing.assert_allclose(loose[0], [0, 0.9, 0.8, 0.7, 0, 0, 0, 0, 0, 0])
    # A gene measured in one background is never dropped by the sign rule.
    info['n_bg'][0, 1] = 1
    np.testing.assert_allclose(g.calibrate(D, line, ['KD', 'G0'], cfg, info=info)[0], [0, 0.9, 0.8, 0.7, 0, 0, 0, 0, 0, 0])
    with pytest.raises(ValueError, match='sign rule needs'):
        g.calibrate(D, line, ['KD', 'G0'], cfg, info=dict(z=info['z']))
    info['n_bg'][0, 1] = 3; info['sign_agree'][0, 1] = np.nan
    assert g.calibrate(D, line, ['KD', 'G0'], cfg, info=info)[0, 1] == 0  # missing agreement is not unanimity
    line.gate = c.expression_gate([100] * 8 + [5, 100])
    assert g.calibrate(D, line, ['KD', 'G0'], cfg, info=info)[0, 8] == 0  # strict boundary


def test_detection_floor(v2a):
    g, c = v2a['generate'], v2a['common']
    rng = np.random.default_rng(0)
    genes = ['KD', 'HIGH', 'LOW', 'BIG', 'OFF']
    cpm = np.column_stack([np.full(200, 50.0), rng.normal(400, 40, 200), rng.normal(20, 15, 200).clip(0), rng.normal(400, 40, 200),
                           rng.normal(400, 40, 200)])
    line = line_with(g, c, genes, cpm); line.gate = np.ones(5, bool)
    cfg = g.Cfg(argparse.Namespace(call='z', min_effect=0.01, detect_floor=True, z_req=4.0, alpha=1.0))
    up, down = g.detection_floor(line, cfg)
    m, v = line.cpm, line.cpm_var
    z_at = lambda l: m * np.abs(np.exp2(l) - 1) / np.sqrt(v * np.exp2(l) ** 2 / cfg.n_cells + v / line.n_reference)
    np.testing.assert_allclose(z_at(up)[1:], 4.0, rtol=1e-3)          # the floor is exactly where the test reaches z_req
    assert up[2] > up[1]                                                # a noisier, lower gene needs a bigger change
    D = np.array([[0.0, 0.02, -0.02, 3.0, 0.0]])
    out = g.calibrate(D, line, ['KD'], cfg)
    assert out[0, 1] == pytest.approx(up[1]) and out[0, 2] == pytest.approx(-down[2])   # raised to the floor, sign kept
    assert out[0, 3] == pytest.approx(3.0) and out[0, 4] == 0                         # never shrunk; uncalled stays 0
    no_floor = g.calibrate(D, line, ['KD'], g.Cfg(argparse.Namespace(call='z', min_effect=0.01, alpha=1.0)))
    np.testing.assert_allclose(no_floor, [[0, 0.02, -0.02, 3.0, 0]])
    cfg.alpha = 0
    zero_alpha = g.calibrate(D, line, ['KD'], cfg)
    np.testing.assert_allclose(zero_alpha[0, 1:4], [up[1], -down[2], up[3]])
    assert zero_alpha[0, 4] == 0


# ---------------------------------------------------------------- edit 6
def test_per_target_table(v2a, tmp_path):
    g = v2a['generate']; features = [f'G{i}' for i in range(10)]
    real = pd.DataFrame(dict(target='T', feature=features, log2_fold_change=[1.0] * 5 + [0.0] * 5, p_value=0.0,
                             p_adj=[0.01] * 5 + [0.9] * 5, abs_log2_fold_change=0.0))
    pred = real.copy(); pred['p_adj'] = [0.01, 0.01, 0.9, 0.9, 0.9, 0.01, 0.9, 0.9, 0.9, 0.9]
    pred['log2_fold_change'] = [1.0, -1.0] + [0.0] * 8
    real.to_parquet(tmp_path / 'de_real.parquet'); pred.to_parquet(tmp_path / 'de_pred.parquet')
    row = g.per_target_table(tmp_path).loc['T']
    assert (row.n_called, row.n_real, row.overlap, row.same_direction) == (3, 5, 2, 1)
    assert row.precision == pytest.approx(2 / 3) and row.chance_precision == pytest.approx(0.5)
    assert row.enrichment == pytest.approx(4 / 3) and row.direction_accuracy == pytest.approx(0.5) and row.recall == pytest.approx(0.4)
    summary = g.per_target_summary(g.per_target_table(tmp_path))
    assert summary['median_direction_accuracy'] == 0.5


def test_scaled_score_requires_all_six_metrics_and_their_mean(v2a, tmp_path):
    g = v2a['generate']; metrics = list(g.VCC2026_METRICS)
    values = [0.5, 0, -0.2, -0.3, -0.1, 0.1]
    frame = pd.DataFrame({'metric': metrics + ['avg_score'], 'from_replicate': values + [np.mean(values)]})
    path = tmp_path / 'scaled.csv'; frame.to_csv(path, index=False)
    score = g.read_scaled(path)
    assert score['overall'] == pytest.approx(np.mean(values)) and score[metrics[3]] == -0.3
    frame.iloc[:-2].to_csv(path, index=False)
    with pytest.raises(ValueError, match='six competition metrics'): g.read_scaled(path)
    frame.loc[6, 'from_replicate'] = 0.9; frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match='mean of the six'): g.read_scaled(path)
    frame.loc[0, 'from_replicate'] = np.nan; frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match='non-finite'): g.read_scaled(path)


def test_gene_split_cross_validation(v2a):
    d = v2a['adhoc']; rng = np.random.default_rng(7)
    x = rng.normal(size=(100, 5)); y = x @ [1, -2, 3, 0, 1] + 0.5
    np.testing.assert_allclose(d.r2_with_cv(y, x), [1, 1])
    assert np.isnan(d.r2_with_cv(np.ones(100), x)).all()


def test_wilcoxon_blocks_preserve_full_gene_tests(v2a):
    from scipy import stats
    rng = np.random.default_rng(3)
    pred, real = rng.poisson(2, (20, 13)), rng.poisson(3, (30, 13))
    p = v2a['adhoc'].wilcoxon_blocks(sp.csc_matrix(pred), sp.csc_matrix(real), block=4)
    np.testing.assert_array_equal(p, stats.mannwhitneyu(pred, real, axis=0).pvalue)


def test_tune_chains_stages_on_the_same_reference(v2a, tmp_path, monkeypatch):
    g = v2a['generate']; submitted = tmp_path / 'submitted.json'
    submitted.write_text(json.dumps(dict(z=3, min_effect=0.1, alpha=2, pi=0, up_w=0.75)))
    stages = []
    def benchmark(a):
        arms = json.loads(Path(a.variants).read_text()); stages.append((a, arms))
        if len(stages) == 1:
            out = Path(a.out); out.mkdir()
            by_name = {x['label']: x for x in arms}
            (out / 'benchmark.json').write_text(json.dumps(dict(settings=by_name)))
            # Reference wins overall, but stage B refines the best top-N candidate.
            pd.DataFrame({'arm': ['reference_v2_submitted', 'topn500_floor1'], 'overall': [0.9, 0.7]}).to_csv(out / 'scaled_scores.csv', index=False)
    monkeypatch.setattr(g, 'cmd_benchmark', benchmark)
    g.cmd_tune(argparse.Namespace(out=str(tmp_path / 'run'), submitted_settings=str(submitted), z_req=5, seed=2026, n_targets=50))
    first, second = stages
    assert len(first[1]) == 9 and len(second[1]) == 7
    assert second[0].reference == first[0].out and first[0].seed == second[0].seed == 2026
    assert first[0].n_targets == second[0].n_targets == 50
    assert first[1][-1] == second[1][-1]
    assert {(a['alpha'], a['sign_rule']) for a in second[1][:-1]} == {(alpha, sign) for alpha in [1, 2, 3] for sign in [True, False]}
    assert all(a['top_n'] == 500 and a['detect_floor'] and a['z_req'] == 5 for a in second[1][:-1])


def test_submission_with_top_n_settings(v2a, tmp_path):
    g, c = v2a['generate'], v2a['common']
    genes = ['KD', 'G', 'H']; pd.DataFrame({'gene_name': genes}).to_csv(tmp_path / 'gene_names.csv', index=False)
    pd.DataFrame({'target_gene': ['KD'], 'n_cells': [400]}).to_csv(tmp_path / 'pert_counts.csv', index=False)
    ad.AnnData(sp.csr_matrix([[30, 20, 25]] * 30), obs=pd.DataFrame({'target_gene': ['non-targeting'] * 30}, index=[str(i) for i in range(30)]),
               var=pd.DataFrame(index=genes)).write_h5ad(tmp_path / 'context_A.h5ad')
    pred = tmp_path / 'pred'; pred.mkdir()
    for name, value in [('genes', genes), ('targets', ['KD']), ('re', [[np.nan, 0.5, 0.1]]), ('re_se', [[np.nan, 0.1, 0.1]]),
                        ('re_z', [[np.nan, 5, 1]]), ('n_bg', [[np.nan, 2, 2]]), ('sign_agree', [[np.nan, 1, 1]]), ('ontarget_lfc', [-1])]:
        np.save(pred / f'{name}.npy', np.asarray(value))
    (tmp_path / 'settings.json').write_text(json.dumps(dict(schema='vcc_2026_2_meta_v2', delta_array='re',
        settings=dict(call='topn', top_n=1, sign_rule=True, detect_floor=True, z_req=4.0, alpha=1.0, z=0, min_effect=0, pi=0, up_w=0.75, n_cells=400))))
    result = subprocess.run([sys.executable, g.__file__, 'submit', '--data', str(tmp_path), '--lines', 'A', '--delta', f'A={pred}',
                             '--delta-array', 'selected', '--settings', str(tmp_path / 'settings.json'), '--out', str(tmp_path / 'out')],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads((tmp_path / 'out/generator_report.json').read_text())
    assert report['settings']['call'] == 'topn' and report['settings']['top_n'] == 1 and report['settings']['detect_floor'] is True
    c.validate_counts(ad.read_h5ad(tmp_path / 'out/prediction.h5ad').X)


# ---------------------------------------------------------------- edits 1, 2, 7: the benchmark on the real evaluator
def test_benchmark_scaled_selection_and_reference_reuse(v2a, tmp_path):
    """Real CPU evaluator on a small synthetic H1: bundle, scaled scores, per-target files, selection, stage reuse."""
    binary = Path(sys.executable).parent / 'cell-eval2'
    if not binary.exists(): pytest.skip('official cell-eval2 is not installed')
    g = v2a['generate']; genes = np.array(['K0', 'K1', 'K2'] + [f'G{i}' for i in range(37)])
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
                        ('re_z', np.abs(effects) / 0.1), ('n_bg', np.full(effects.shape, 2.0)), ('sign_agree', np.ones(effects.shape)),
                        ('ontarget_lfc', [-1.5] * 3)]:
        np.save(pred / f'{name}.npy', np.asarray(value))
    (pred / 'prediction_report.json').write_text(json.dumps(dict(schema='vcc_2026_2_meta_v2', models=['re'], heldout_background='H1',
                                                                 source_backgrounds=['HCT116', 'HEK293T'])))
    arms = [dict(method='re', call='topn', top_n=12, sign_rule=True, detect_floor=True, alpha=1.0, label='topn12_floor1'),
            dict(method='re', call='z', z=3.0, min_effect=0.1, alpha=2.0, label='reference')]
    (tmp_path / 'a.json').write_text(json.dumps(arms))
    common = ['--cells', str(tmp_path / 'h1.h5ad'), '--genes', str(tmp_path / 'genes.csv'), '--targets', str(tmp_path / 'targets.csv'),
              '--predictions', str(pred), '--cell-eval', str(binary), '--num-threads', '1', '--n-cells', '100', '--min-cells', '50',
              '--anchor-splits', '5']
    run = lambda *extra: subprocess.run([sys.executable, g.__file__, 'benchmark', *common, *extra], capture_output=True, text=True, timeout=900)
    result = run('--variants', str(tmp_path / 'a.json'), '--out', str(tmp_path / 'stage_a'))
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-3000:]
    a = tmp_path / 'stage_a'
    assert (a / 'bundle/manifest.json').exists() and not (a / 'baseline/baseline_pred.h5ad').exists()
    assert (a / 'anchor/anchor_agg.parquet').exists()
    assert json.loads((a / 'bundle/manifest.json').read_text())['rule_digest']
    scaled = pd.read_csv(a / 'scaled_scores.csv', index_col='arm')
    assert set(scaled.index) == {'topn12_floor1', 'reference'} and 'overall' in scaled and scaled['overall'].notna().all()
    assert scaled['overall'].is_monotonic_decreasing
    for arm in scaled.index:
        assert (a / arm / 'scaled.csv').exists() and (a / arm / 'per_target.csv').exists() and not (a / f'{arm}.h5ad').exists()
    assert len(pd.read_csv(a / 'per_target_summary.csv')) == 2
    selected = json.loads((a / 'selected_generator.json').read_text())
    assert selected['arm'] == scaled.index[0] and selected['selection'] == 'scaled overall (leaderboard formula)'
    assert selected['schema'] == 'vcc_2026_2_meta_v2' and {'call', 'top_n', 'detect_floor', 'sign_rule', 'z_req', 'alpha'} <= set(selected['settings'])
    # Stage B reuses stage A's reference and bundle: no new real file, cache or bundle.
    (tmp_path / 'b.json').write_text(json.dumps([dict(arms[0], alpha=2.0, label='topn12_floor1_alpha2'), arms[1]]))
    result = run('--variants', str(tmp_path / 'b.json'), '--reference', str(a), '--out', str(tmp_path / 'stage_b'))
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-3000:]
    b = tmp_path / 'stage_b'
    assert not (b / 'real.h5ad').exists() and not (b / 'bundle').exists() and not (b / 'real_cache').exists()
    assert json.loads((b / 'reference_info.json').read_text())['targets'] == 3
    assert set(pd.read_csv(b / 'scaled_scores.csv', index_col='arm').index) == {'topn12_floor1_alpha2', 'reference'}
    assert (b / 'reference/reused_from.json').exists()
    pd.testing.assert_frame_equal(pd.read_csv(a / 'reference/per_target.csv'), pd.read_csv(b / 'reference/per_target.csv'))

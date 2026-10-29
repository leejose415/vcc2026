"""Freeze check, part 1: run vcc_2026_v7.ipynb's own setup, data-loader and model cells (unchanged source) for EvalA only,
writing results/v7/freeze_check/EvalA/input.npz instead of results/v7/EvalA/."""
import json
from pathlib import Path

V = Path.home() / 'vcc'
cells = [''.join(c['source']) for c in json.load(open(V / 'vcc_2026_v7.ipynb'))['cells'] if c['cell_type'] == 'code']
setup, data, model = cells[0], cells[1], cells[2]
assert 'CONFIG = dict(' in setup and 'S.load_tables' in data and 'def config1' in model
g = {'__name__': '__main__'}
exec(compile(setup, 'v7_setup', 'exec'), g)
g['OUT'] = V / 'results/v7/freeze_check'; g['RUN'] = ['EvalA']; g['todo'] = ['EvalA']
(g['OUT'] / 'EvalA').mkdir(parents=True, exist_ok=True)
exec(compile(data, 'v7_data', 'exec'), g)
exec(compile(model, 'v7_model', 'exec'), g)

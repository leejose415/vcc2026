"""Run a queue of transformer sweep runs, one after another on the GPU, from a plan file that may change while it runs.

  python scripts/sweep_queue.py results/sweeps/<date>/plan.json
The plan is re-read before every run, so tiers can be cut, reordered or extended without stopping the queue. Each entry:
  {"run_id": "lr_1e-4", "knob": "learning rate", "value": 1e-4, "seed": 2026, "set": {"train.lr": 1e-4}, "epochs": 12}
A run is skipped when <out>/runs/<run_id>/done.json exists; a failed run is logged and the queue moves on.
Validation mode only (the base config is code/transformer/configs/dev.yaml with PIE evidence off, no early stopping).
Logs: <out>/logs/<run_id>.log; per-epoch metrics: <out>/runs/<run_id>/metrics.csv.
"""
import json
import subprocess
import sys
import time
from pathlib import Path

V = Path.home() / 'vcc'


def main():
    plan_path = Path(sys.argv[1]).resolve(); out = plan_path.parent
    (out / 'logs').mkdir(parents=True, exist_ok=True)
    tried = set()
    while True:
        plan = json.loads(plan_path.read_text())
        todo = [r for r in plan['runs'] if r['run_id'] not in tried and not (out / 'runs' / r['run_id'] / 'done.json').exists()]
        if not todo:
            print(f'{time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())} UTC queue empty', flush=True); break
        r = todo[0]; tried.add(r['run_id'])
        sets = dict(plan.get('base_set', {}), **r.get('set', {}))
        cmd = [sys.executable, '-m', 'transformer.train', 'train', '--config', str(V / plan.get('base_config', 'code/transformer/configs/dev.yaml')),
               '--name', r['run_id'], '--runs-dir', str(out / 'runs'), '--epochs', str(r.get('epochs', plan.get('epochs', 12))),
               '--set', f'seed={r.get("seed", 2026)}', '--set', 'sweep.enabled=true', '--set', 'train.early_stop=false',
               '--set', f'sweep.knob={r["knob"]}', '--set', f'sweep.value={r["value"]}']
        for k, v in sets.items():
            cmd += ['--set', f'{k}={v}']
        t0 = time.time()
        print(f'{time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())} UTC start {r["run_id"]}', flush=True)
        with open(out / 'logs' / f'{r["run_id"]}.log', 'w') as log:
            rc = subprocess.run(cmd, cwd=V, env=dict(__import__('os').environ, PYTHONPATH=str(V / 'code')), stdout=log,
                                stderr=subprocess.STDOUT).returncode
        print(f'{time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())} UTC end {r["run_id"]} rc={rc} '
              f'({(time.time() - t0) / 60:.1f} min)', flush=True)


if __name__ == '__main__':
    main()

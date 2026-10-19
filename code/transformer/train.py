"""Train / smoke-test / benchmark the vcc_2026_4 transformer.

  PYTHONPATH=code python -m transformer.train train --config code/transformer/configs/dev.yaml [--name N] [--pie on|off|auto]
  PYTHONPATH=code python -m transformer.train smoke --config code/transformer/configs/tiny.yaml
  PYTHONPATH=code python -m transformer.train bench --config code/transformer/configs/dev.yaml
Run from ~/vcc. Training: AdamW, linear warm-up then cosine, grad clip, bf16 autocast, per-epoch CSV log, best checkpoint
by held-out-gene loss with early stopping. Runs go to lib/vcc_2026_4/runs/<name>/.
"""
import argparse
import csv
import json
import math
import time
from pathlib import Path

import numpy as np
import torch

from . import config as C
from .data import Bank
from .losses import losses
from .model import VCCTransformer


def seed_all(seed):
    torch.manual_seed(seed); np.random.seed(seed); torch.cuda.manual_seed_all(seed)


def lr_at(step, total, warmup, base):
    if step < warmup:
        return base * (step + 1) / warmup
    p = min(1.0, (step - warmup) / max(1, total - warmup))
    return base * 0.5 * (1 + math.cos(math.pi * p))


def step_batch(model, bank, rows, cfg, opt=None, train=False, gen=None):
    b = bank.batch(rows=rows, train=train, generator=gen)
    with torch.autocast('cuda', dtype=torch.bfloat16, enabled=cfg['train']['bf16']):
        out = model(b)
    loss, parts = losses(out, b, cfg)
    if train:
        opt.zero_grad(set_to_none=True); loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), cfg['train']['grad_clip']); opt.step()
    return parts


@torch.no_grad()
def evaluate(model, bank, rows, cfg, max_batches=None):
    model.eval(); bs = cfg['train']['batch_size']; acc = {}; n = 0
    for i, s in enumerate(range(0, len(rows), bs)):
        if max_batches and i >= max_batches: break
        p = step_batch(model, bank, rows[s:s + bs], cfg)
        for k, v in p.items(): acc[k] = acc.get(k, 0.0) + v * p['measured']
        n += p['measured']
    model.train()
    return {k: v / max(n, 1) for k, v in acc.items() if k != 'measured'}


def make_model(cfg, bank):
    return VCCTransformer(cfg, len(bank.genes)).to(bank.device)


def train(cfg, bank=None, log=print):
    seed_all(cfg['seed'])
    bank = bank or Bank(cfg)
    model = make_model(cfg, bank); tc = cfg['train']
    opt = torch.optim.AdamW(model.parameters(), lr=tc['lr'], weight_decay=tc['weight_decay'])
    run = Path(cfg['paths']['runs']) / cfg['name']; run.mkdir(parents=True, exist_ok=True)
    (run / 'config.json').write_text(json.dumps(cfg, indent=1))
    np.save(run / 'gene_cols.npy', bank.gcols)
    (run / 'heldout_genes.json').write_text(json.dumps(sorted(bank.heldout_keys)))
    bs = tc['batch_size']; spe = tc['steps_per_epoch'] or math.ceil(len(bank.train_rows) / bs)
    total = spe * tc['epochs']; gen = torch.Generator(device=bank.device); gen.manual_seed(cfg['seed'])
    rng = np.random.default_rng(cfg['seed'])
    n_params = sum(p.numel() for p in model.parameters())
    log(f'run {run} | params {n_params / 1e6:.2f} M | {spe} steps/epoch x {tc["epochs"]} epochs | batch {bs}')
    best, bad, step = float('inf'), 0, 0
    with open(run / 'log.csv', 'w', newline='') as fh:
        w = None
        for epoch in range(1, tc['epochs'] + 1):
            t0 = time.time(); order = rng.permutation(bank.train_rows); acc = {}; nb = 0
            for i in range(spe):
                rows = order[(i * bs) % len(order):(i * bs) % len(order) + bs]
                if len(rows) < 2: continue
                for g in opt.param_groups: g['lr'] = lr_at(step, total, tc['warmup'], tc['lr'])
                p = step_batch(model, bank, rows, cfg, opt, train=True, gen=gen); step += 1; nb += 1
                for k, v in p.items(): acc[k] = acc.get(k, 0.0) + v
                if not math.isfinite(p['loss']): raise FloatingPointError(f'non-finite loss at step {step}: {p}')
            tr = {f'train_{k}': v / max(nb, 1) for k, v in acc.items() if k not in ('measured', 'pos_weight')}
            va = {f'val_{k}': v for k, v in evaluate(model, bank, bank.val_rows, cfg, tc['eval_batches']).items()
                  if k != 'pos_weight'}
            rec = dict(epoch=epoch, step=step, lr=opt.param_groups[0]['lr'], gamma=float(model.gamma), seconds=time.time() - t0,
                       **tr, **va)
            if w is None:
                w = csv.DictWriter(fh, fieldnames=list(rec)); w.writeheader()
            w.writerow(rec); fh.flush()
            improved = va['val_loss'] < best
            if improved:
                best, bad = va['val_loss'], 0
                torch.save(dict(model=model.state_dict(), epoch=epoch, val_loss=best), run / 'best.pt')
            else:
                bad += 1
            log(f"epoch {epoch:2d} | train {tr['train_loss']:.4f} (huber {tr['train_huber']:.4f} cos {tr['train_cosine']:.4f} "
                f"bce {tr['train_bce']:.4f}) | held-out {va['val_loss']:.4f} | gamma {float(model.gamma):.3f} | "
                f"{rec['seconds']:.0f}s{' *' if improved else ''}")
            if bad >= tc['patience']:
                log(f'early stop: no held-out improvement for {bad} epochs'); break
    torch.save(dict(model=model.state_dict(), epoch=epoch), run / 'last.pt')
    (run / 'done.json').write_text(json.dumps(dict(best_val_loss=best, epochs_run=epoch, steps=step)))
    return model, bank, run


def smoke(cfg):
    """Overfit 64 examples (tiny model, no dropout / noise), then 200 real training steps."""
    seed_all(cfg['seed']); bank = Bank(cfg); res = {}
    sc = dict(dict(overfit_steps=600, overfit_lr=1e-3, overfit_rows=64), **cfg.get('smoke', {}))
    model = make_model(cfg, bank); opt = torch.optim.AdamW(model.parameters(), lr=sc['overfit_lr'], weight_decay=0)
    rows = bank.train_rows[np.random.default_rng(0).choice(len(bank.train_rows), sc['overfit_rows'], replace=False)]
    first = None; n = sc['overfit_steps']
    for i in range(n):
        p = step_batch(model, bank, rows, cfg, opt, train=True)
        first = first or p
        if i % max(1, n // 6) == 0 or i == n - 1: print(f'  overfit step {i:3d}: loss {p["loss"]:.4f} (huber {p["huber"]:.4f} cos {p["cosine"]:.4f} bce {p["bce"]:.4f})', flush=True)
    res['overfit'] = dict(first=first, last=p, ratio=p['loss'] / first['loss'])
    seed_all(cfg['seed']); model = make_model(cfg, bank)
    opt = torch.optim.AdamW(model.parameters(), lr=cfg['train']['lr'], weight_decay=cfg['train']['weight_decay'])
    rng = np.random.default_rng(1); bs = cfg['train']['batch_size']; hist = []; t0 = time.time()
    for i in range(200):
        for g in opt.param_groups: g['lr'] = lr_at(i, 200, 20, cfg['train']['lr'])
        p = step_batch(model, bank, rng.choice(bank.train_rows, bs, replace=False), cfg, opt, train=True); hist.append(p['loss'])
        if i % 50 == 0 or i == 199: print(f'  real step {i:3d}: loss {p["loss"]:.4f}', flush=True)
    res['real'] = dict(first10=float(np.mean(hist[:10])), last10=float(np.mean(hist[-10:])), seconds=time.time() - t0,
                       finite=bool(np.all(np.isfinite(hist))))
    res['heldout'] = evaluate(model, bank, bank.val_rows, cfg, max_batches=20)
    print(json.dumps(res, indent=1, default=float))
    return res


def bench(cfg, sizes=(500, 5000, 18533), max_bs=256):
    """Peak GPU memory and time per training step at several gene counts and batch sizes; max batch that fits."""
    out = []
    for G in sizes:
        c = C._merge(cfg, dict(gene_set='full' if G >= 18533 else 'dev', n_genes=G))
        bank = Bank(c, verbose=False); torch.cuda.synchronize(); base = torch.cuda.memory_allocated()
        model = make_model(c, bank); opt = torch.optim.AdamW(model.parameters(), lr=1e-4)
        bs, best = 1, None
        while bs <= max_bs:
            try:
                rows = bank.train_rows[:bs]
                step_batch(model, bank, rows, c, opt, train=True)                       # warm-up
                torch.cuda.synchronize(); torch.cuda.reset_peak_memory_stats(); t = time.time()
                for _ in range(3): step_batch(model, bank, rows, c, opt, train=True)
                torch.cuda.synchronize(); sec = (time.time() - t) / 3
                peak = (torch.cuda.max_memory_allocated() - base) / 2**30
                best = dict(genes=len(bank.gcols), batch=bs, sec_per_step=sec, peak_gib=peak,
                            data_on_gpu_gib=base / 2**30,
                            epoch_minutes=sec * math.ceil(len(bank.train_rows) / bs) / 60)
                print(f"  genes {len(bank.gcols):5d} batch {bs:3d}: {sec * 1000:7.1f} ms/step, peak {peak:5.2f} GiB "
                      f"(+ data {base / 2**30:.2f} GiB), epoch {best['epoch_minutes']:.1f} min", flush=True)
                bs *= 2
            except torch.OutOfMemoryError:
                print(f'  genes {len(bank.gcols):5d} batch {bs:3d}: out of memory', flush=True); break
            finally:
                opt.zero_grad(set_to_none=True); torch.cuda.empty_cache()
        out.append(best); del model, opt, bank; torch.cuda.empty_cache()
    print(json.dumps(out, indent=1))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('cmd', choices=['train', 'smoke', 'bench'])
    ap.add_argument('--config'); ap.add_argument('--name'); ap.add_argument('--pie'); ap.add_argument('--mode')
    ap.add_argument('--gene-features'); ap.add_argument('--epochs', type=int); ap.add_argument('--n-genes', type=int)
    ap.add_argument('--gene-set'); ap.add_argument('--out')
    a = ap.parse_args()
    over = {k: v for k, v in dict(name=a.name, pie=a.pie, mode=a.mode, gene_features=a.gene_features, n_genes=a.n_genes,
                                  gene_set=a.gene_set).items() if v is not None}
    if a.epochs: over['train'] = dict(epochs=a.epochs)
    cfg = C.load(a.config, **over)
    print(f"MODE: {cfg['mode']} | name {cfg['name']} | genes {cfg['gene_set']} ({cfg['n_genes']}) | PIE {cfg['pie']} | "
          f"gene features {cfg['gene_features']}", flush=True)
    res = dict(train=lambda: str(train(cfg)[2]), smoke=lambda: smoke(cfg), bench=lambda: bench(cfg))[a.cmd]()
    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True); Path(a.out).write_text(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()

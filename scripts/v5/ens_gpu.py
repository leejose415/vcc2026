"""Item 3 predictions (GPU, no training): dev_nopie zero-shot H1 predictions with no source masked and with each loaded
evidence source masked in turn (TTA), and from last.pt (checkpoint averaging: only best.pt and last.pt were saved).

  cd ~/vcc && PYTHONPATH=code ~/venv/bin/python ~/v5/ens_gpu.py
"""
from pathlib import Path

import numpy as np
import torch

from transformer.data import SLOT
from transformer.evaluate import h1_truth, load_run, predict

RUN = Path.home() / 'vcc/lib/vcc_2026_4/runs/dev_nopie'
out = {}
for ckpt in ('best.pt', 'last.pt'):
    cfg, bank, model = load_run(RUN, device='cuda', ckpt=ckpt)
    obs, keys, *_ = h1_truth(bank)
    kd, line, trow = bank.rows_for('H1', keys)
    tag = ckpt.split('.')[0]
    out[f'{tag}_none'] = predict(model, bank, kd, line, trow)['delta']
    if ckpt == 'best.pt':
        full = bank.loaded_slots.clone()
        for src in bank.loaded:
            bank.loaded_slots = full[full != SLOT[src]]
            out[f'best_mask_{src}'] = predict(model, bank, kd, line, trow)['delta']
            print('masked', src, float(np.abs(out[f'best_mask_{src}'] - out['best_none']).mean()), flush=True)
        bank.loaded_slots = full
    del model, bank; torch.cuda.empty_cache()
saved = np.load(RUN / 'h1_pred_g3/raw_predictions.npz')
print('best_none vs saved delta max |diff|:', float(np.abs(out['best_none'] - saved['delta']).max()))
np.savez(Path.home() / 'v5/ens_preds.npz', gene_cols=saved['gene_cols'], **out)
print('wrote', sorted(out))

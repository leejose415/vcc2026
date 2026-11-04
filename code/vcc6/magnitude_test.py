"""Magnitude test on the arrays already saved (no refit): the Step 3 arm B model run on the RPE1 dev 300.

  S  the shrunk posterior mean (lfc), as run          U  lfc_unshrunk, no shrinkage
  G  Gaussian shrink only: lfc_unshrunk x v0 / (v0 + v_model), v_model = se_model^2, the slab weight not applied
  H  0.5 x lfc_unshrunk (the factor the unshrunk calibration slope names)

Calls for every variant: stage1.predictive_calls (BH on |lfc| / se_measurement over the run's tested gate), the model's
own rule. Effect space against the replogle_rpe1 table on its tested genes (as Step 2): r and both slopes
(calibration_stats), NMAE / sign accuracy / AUROC / call count (score_against_truth, summarize_scores). The best of G
and H, chosen before the end-to-end run by: higher pooled r, then calibration slope closer to 1, is written as an
--lfc-npz for run_line.py.

  ~/venv/bin/python code/vcc6/magnitude_test.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
import stage1 as S                                                                         # noqa: E402

R = V / 'results/stage1/2026-10-27_rpe1'
RUN = R / 'step3/armB'
OUT = R / 'magnitude'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    z = np.load(RUN / 'model/predictions.npz', allow_pickle=True)
    summ = json.loads((RUN / 'model/summary.json').read_text())
    v0 = float(summ['prior']['v0'])
    kds = [str(k) for k in z['knockdowns']]; genes = [str(g) for g in z['genes']]
    gate_run = pd.read_csv(RUN / 'controls/controls_baseline.csv', index_col=0)['de_tested'].to_numpy(bool)
    tab = V / 'lib/effects_v5/replogle_rpe1.h5ad'
    _, gate_tab = S.load_baseline(V / 'lib/effects_v5/replogle_rpe1_baseline.csv', np.array(genes, dtype=object))
    truth = S.truth_signature(tab, kds, np.array(genes, dtype=object))
    mu = z['lfc_unshrunk'].astype(np.float64)
    v_model = np.where(np.isfinite(z['se_model']), z['se_model'].astype(np.float64) ** 2, np.inf)
    variants = dict(S=z['lfc'].astype(np.float64), U=mu,
                    G=np.where(np.isfinite(v_model), mu * v0 / (v0 + v_model), 0.0), H=0.5 * mu)
    se_meas = z['se_measurement']
    res, calls = {}, {}
    for name, L in variants.items():
        sig = np.zeros(L.shape, bool) if name != 'S' else z['de_sig']
        if name != 'S':
            for i in range(len(kds)):
                sig[i], _ = S.predictive_calls(L[i], se_meas[i], gate_run & (z['n_sources'][i] > 0))
        calls[name] = sig
        rows = []
        for i, k in enumerate(kds):
            y, e, tsig = truth[k]
            pred = dict(lfc=L[i].astype(np.float32), lfc_unshrunk=z['lfc_unshrunk'][i], de_sig=sig[i],
                        v_post=np.zeros(len(genes), np.float32), se_measurement=se_meas[i])
            rows.append(S.score_against_truth(pred, y, e, tsig, gate_tab))
        es = S.summarize_scores(rows)
        fit = S.calibration_stats(L, {i: truth[k][0] for i, k in enumerate(kds)}, gate_tab)
        res[name] = dict(pearson_r=fit['pearson_r'], calibration_slope=fit['slope_true_on_pred'],
                         transmission_slope=fit['slope_pred_on_true'], nmae=es['nmae'], nmae_zero=es['nmae_zero'],
                         mse=es['mse'], mse_zero=es['mse_zero'], sign_acc_true_de=es['sign_acc_true_de'], auroc=es['auroc'],
                         median_n_called=es['n_called'], median_n_true=es['n_true'], call_count_pearson=es['call_count_pearson'],
                         median_abs_lfc_called=float(np.median(np.abs(L[sig]))) if sig.any() else None)
    t = pd.DataFrame(res).T
    t.to_csv(OUT / 'effect_space.csv')
    print(t.to_string())
    gh = t.loc[['G', 'H']].copy()
    gh['slope_gap'] = (gh.calibration_slope - 1).abs()
    best = gh.sort_values(['pearson_r', 'slope_gap'], ascending=[False, True]).index[0]
    np.savez_compressed(OUT / f'{best}.npz', knockdowns=np.array(kds, dtype=object), genes=np.array(genes, dtype=object),
                        lfc=variants[best].astype(np.float32), de_sig=calls[best], n_sources=z['n_sources'])
    (OUT / 'choice.json').write_text(json.dumps(dict(best=best, rule='higher pooled r, then calibration slope closer to 1',
                                                     v0=v0, source=str(RUN / 'model/predictions.npz')), indent=1) + '\n')
    print('best of G and H:', best)


if __name__ == '__main__':
    main()

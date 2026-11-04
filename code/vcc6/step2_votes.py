"""Step 2: how sources vote. Effect space only, RPE1 held out, the dev 300 list. One load of the 117 tables.

  structure  per background: tables, rows, modalities, multipliers, votes for RPE1 / for H1; per knockdown: Tahoe
             rows and Tahoe's share of the precision weight (arm A weights 1 / (SE^2 x multiplier), median over genes
             and over the bottom expression decile of RPE1's controls); how many of fit_prior's 150 sampled knockdowns
             have a Tahoe row; the Orion tables' backgrounds and their overlap with the 300.
  arms       stage1.predict_line with arm A (vcc_2026_6), B (one vote per background, no drug tables), C (B + Tahoe,
             one vote per background); the metrics on the 300 and on the knockdowns any Tahoe table holds.

  ~/venv/bin/python code/vcc6/step2_votes.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
import stage1 as S                                                                           # noqa: E402

OUT = V / 'results/stage1/2026-10-27_rpe1/step2'
DIRS = [str(V / 'lib/effects_v5'), str(V / 'lib/effects_v5_extra')]
LIST = V / 'results/target_lists/rpe1_dev_300.txt'


def structure(line, want):
    others, genes = line['others'], line['genes']
    rows = []
    for s in line['sources']:
        rows.append(dict(background=S._line_key(s.background), table=s.name, modality=s.modality, reference=s.reference,
                         mult_A=s.mult, mult_BC=S.BackgroundVotes.mult(s), rows=int(s.lfc.shape[0]),
                         rows_in_300=int(sum(len(s.rows_of.get(k, ())) for k in want))))
    t = pd.DataFrame(rows)
    by = t.groupby('background').agg(tables=('table', 'count'), observations=('rows', 'sum'),
                                     observations_in_300=('rows_in_300', 'sum'),
                                     modality=('modality', lambda x: '+'.join(sorted(set(x)))),
                                     multiplier_A=('mult_A', lambda x: '/'.join(f'{v:g}' for v in sorted(set(x)))),
                                     multiplier_BC=('mult_BC', lambda x: '/'.join(f'{v:g}' for v in sorted(set(x)))),
                                     tables_list=('table', lambda x: ' '.join(sorted(x)))).reset_index()
    by['votes_for_rpe1'] = by.background != S._line_key(line['target_src'].background)
    by['votes_for_h1'] = by.background != 'h1'
    t.to_csv(OUT / 'tables.csv', index=False)
    by.to_csv(OUT / 'backgrounds.csv', index=False)
    # Tahoe exposure under arm A's pooling
    tahoe = [s for s in others if s.modality == 'drug']
    m0 = line['basal_cpm']
    expr = np.isfinite(m0) & (m0 > 0)
    bottom = expr & (m0 <= np.nanquantile(m0[expr], 0.10))
    per = []
    for k in want:
        w_all = np.zeros(len(genes)); w_t = np.zeros(len(genes)); n_t = 0
        for s in others:
            for r in s.rows_of.get(k, ()):
                y, e = S._align_row(s, r, genes)
                ok = np.isfinite(y) & np.isfinite(e) & (e > 0) & (np.abs(y) <= S.LFC_CAP)
                w = np.where(ok, 1.0 / (np.maximum(e, S.SE_MIN) ** 2 * s.mult), 0.0)
                w_all += w
                if s.modality == 'drug':
                    w_t += w; n_t += 1
        has = w_all > 0
        share = np.divide(w_t, w_all, out=np.full(len(genes), np.nan), where=has)
        per.append(dict(knockdown=k, tahoe_rows=n_t, tahoe_tables=sum(1 for s in tahoe if k in s.rows_of),
                        tahoe_share_median=float(np.nanmedian(share)) if has.any() else None,
                        tahoe_share_bottom_decile=float(np.nanmedian(share[bottom & has])) if (bottom & has).any() else None,
                        genes_with_weight=int(has.sum())))
    per = pd.DataFrame(per)
    per.to_csv(OUT / 'tahoe_per_knockdown.csv', index=False)
    # fit_prior's 150 under arm A (same owners, rng and draw as stage1.fit_prior)
    owners = {}
    for s in others:
        for k in s.rows_of:
            owners.setdefault(k, []).append(s)
    shared = [k for k, ss in owners.items() if len(ss) >= 2]
    rng = np.random.default_rng(2026)
    drawn = list(rng.choice(np.array(shared, dtype=object), 150, replace=False)) if len(shared) > 150 else shared
    prior_tahoe = [k for k in drawn if any(s.modality == 'drug' for s in owners[k])]
    # Orion
    orion = t[t.table.str.startswith('orion')][['table', 'background', 'modality', 'mult_A', 'mult_BC', 'rows', 'rows_in_300']]
    rep = dict(n_backgrounds=int(len(by)), fit_prior_A=dict(shared_knockdowns=len(shared), drawn=len(drawn),
                                                            drawn_with_tahoe=len(prior_tahoe), examples=prior_tahoe[:10]),
               tahoe_knockdowns_in_300=int((per.tahoe_rows > 0).sum()),
               tahoe_share_median_over_touched=float(per.loc[per.tahoe_rows > 0, 'tahoe_share_median'].median())
               if (per.tahoe_rows > 0).any() else None,
               tahoe_share_bottom_decile_over_touched=float(per.loc[per.tahoe_rows > 0, 'tahoe_share_bottom_decile'].median())
               if (per.tahoe_rows > 0).any() else None,
               orion=orion.to_dict('records'))
    (OUT / 'structure.json').write_text(json.dumps(rep, indent=1, default=str) + '\n')
    return per, rep


def arm_metrics(arm_dir, subset, line, want):
    summ = json.loads((arm_dir / 'summary.json').read_text())
    kd = json.loads((arm_dir / 'per_knockdown_scores.json').read_text())
    pred = np.load(arm_dir / 'predictions.npz', allow_pickle=True)
    truth = S.truth_signature(line['by_stem'][line['target']], want, line['genes'])
    out = {}
    for name, keep in (('all_300', set(want)), ('tahoe_touched', set(subset))):
        rows = [r for r in kd if r['knockdown'] in keep]
        es = S.summarize_scores(rows) if rows else None
        idx = [i for i, k in enumerate(want) if k in keep and k in truth]
        tr = {i: truth[want[i]][0] for i in idx}
        fit = dict(shrunk=S.calibration_stats(pred['lfc'], tr, line['tested']),
                   unshrunk=S.calibration_stats(pred['lfc_unshrunk'], tr, line['tested']))
        n_called = [r['n_called'] for r in rows]; n_true = [r['n_true'] for r in rows]
        out[name] = dict(knockdowns=len(rows), effect_space=es, fit=fit,
                         median_n_called=float(np.median(n_called)) if rows else None,
                         median_n_true=float(np.median(n_true)) if rows else None)
    out['prior'] = summ['prior']
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    want = [k for k in LIST.read_text().split() if k]
    line = S.load_line(DIRS, 'replogle_rpe1')
    per, rep = structure(line, want)
    print(json.dumps(rep, indent=1, default=str)[:3000], flush=True)
    subset = per.loc[per.tahoe_rows > 0, 'knockdown'].tolist()
    res = {}
    for arm in ('A', 'B', 'C'):
        d = OUT / f'arm{arm}'
        S.predict_line(line, str(d), knockdowns_file=str(LIST), arm=arm)
        res[arm] = arm_metrics(d, subset, line, want)
        (OUT / 'arms.json').write_text(json.dumps(res, indent=1, default=str) + '\n')
    print(json.dumps(res, indent=1, default=str), flush=True)


if __name__ == '__main__':
    main()

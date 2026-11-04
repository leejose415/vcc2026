"""How much Tahoe drives the study's similarity (Jose 2026-10-27). The similarity of shrink_study.build, computed twice:
  all     the 2,000 genes most variable across all 117 tables (among genes every table measures)
  noTahoe the same over the 67 non-Tahoe tables only (genes every non-Tahoe table measures)
Profiles: each table's log1p_bulk_50k_full. Reports the overlap of the two 2,000-gene sets and, per fold target, the
Spearman correlation of its similarity to every source table under the two gene sets (all sources scored by both; also
non-Tahoe sources only). Writes <study>/tables/similarity_tahoe_check.parquet and prints a summary.

  ~/venv/bin/python code/vcc6/similarity_check.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
import stage1 as S                                                                         # noqa: E402
from shrink_study import DIRS, LINES, OUT                                                  # noqa: E402


def hv_genes(P, cols):
    common = P.index[P[cols].notna().all(axis=1)]
    return common, P.loc[common, cols].var(axis=1).sort_values(ascending=False).index[:2000]


def cosine_to(P, genes, target):
    M = P.loc[genes]
    ok = M.notna().all(axis=0)
    Mn = M.loc[:, ok] / np.linalg.norm(M.loc[:, ok].to_numpy(float), axis=0)
    return pd.Series(Mn.T.to_numpy(float) @ Mn[target].to_numpy(float), index=Mn.columns)


def main():
    by = {p.stem: p for p in S.discover_tables(DIRS)}
    P = pd.DataFrame({n: pd.read_csv(p.with_name(n + '_baseline.csv'), index_col=0)['log1p_bulk_50k_full'] for n, p in by.items()})
    tahoe = [n for n in P.columns if n.startswith('tahoe')]
    rest = [n for n in P.columns if n not in tahoe]
    c_all, g_all = hv_genes(P, list(P.columns))
    c_rest, g_rest = hv_genes(P, rest)
    overlap = len(set(g_all) & set(g_rest))
    rows, summ = [], dict(tables=len(P.columns), tahoe_tables=len(tahoe), common_genes_all=len(c_all),
                          common_genes_noTahoe=len(c_rest), hvg_overlap=overlap, hvg_overlap_share=overlap / 2000, targets={})
    for line, t in LINES.items():
        a = cosine_to(P, g_all, t).drop(t); b = cosine_to(P, g_rest, t).drop(t)
        both = a.index.intersection(b.index)
        nt = [s for s in both if not s.startswith('tahoe')]
        rho_all = spearmanr(a[both], b[both]).statistic; rho_nt = spearmanr(a[nt], b[nt]).statistic
        summ['targets'][line] = dict(sources_scored_both=len(both), spearman_all_sources=float(rho_all),
                                     spearman_non_tahoe_sources=float(rho_nt),
                                     sources_dropped_noTahoe_genes=sorted(set(a.index) - set(b.index)))
        rows += [dict(target=line, source=s, sim_all=float(a[s]), sim_noTahoe=float(b.get(s, np.nan))) for s in a.index]
    pd.DataFrame(rows).to_parquet(OUT / 'tables/similarity_tahoe_check.parquet')
    (OUT / 'tables/similarity_tahoe_check.json').write_text(json.dumps(summ, indent=1) + '\n')
    print(json.dumps(summ, indent=1))


if __name__ == '__main__':
    main()

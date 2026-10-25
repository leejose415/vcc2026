"""Reconcile completed QC outputs and write one inventory row per cell line."""
import json,re
from pathlib import Path
import numpy as np,pandas as pd
from qc_new_context import ROOT,OUT,BAD
records=[];dataset_rows=[];details=[]
registry=pd.read_csv(ROOT/'lib/ipsc/feng_lfc/00_raw_data_id_mapping.csv');full={}
for text in registry.lines_used.dropna():
 for name in str(text).split(';'):
  if '-' in name:full[name.split('-',1)[1]]=name
for p in sorted(OUT.glob('*/summary.json')):
 s=json.loads(p.read_text());source=p.parent.name
 if source.startswith('Feng'):study='Feng'
 elif source.startswith('Jiang'):study='Jiang'
 elif source.startswith('KOLF'):study='Doctor/KOLF'
 else:study=re.match(r'^([A-Za-z]+\d{4})',source).group(1)
 for r in s['rows']:
  records.append(dict(**r,study=study,full_cell_line=full.get(r['cell_line'],r['cell_line']),qc_status='complete'))
 d=pd.read_csv(p.parent/'guide_qc.csv')
 if len(d):
  for (line,status),g in d.groupby(['cell_line','efficacy_status']):details.append(dict(source=source,cell_line=line,efficacy_status=status,construct_contexts=len(g),cells=int(g.n_cells.sum())))
if records:
 per=pd.DataFrame(records);per.to_csv(OUT/'inventory_by_dataset_and_cell_line.csv',index=False)
 rows=[]
 for line,g in per.groupby('cell_line',sort=True):
  rows.append(dict(cell_line=line,full_cell_line=g.full_cell_line.iloc[0],studies=g.study.nunique(),expression_matrices=g.source.nunique(),cell_records=int(g.cells.sum()),cells_kept=int(g.cells_kept.sum()),controls_kept=int(g.controls_kept.sum()),target_condition_groups=int(g.targets_total.sum()),target_condition_groups_kept=int(g.targets_kept.sum()),study_names='; '.join(sorted(g.study.unique())),sources='; '.join(g.source)))
 table=pd.DataFrame(rows);table.to_csv(OUT/'inventory_by_cell_line.csv',index=False)
 pd.DataFrame(details).to_csv(OUT/'efficacy_status_by_cell_line.csv',index=False)
 text='# Newly downloaded datasets: inventory and QC\n\n'
 text+=f'Completed expression matrices: **{per.source.nunique()} / 22**. Studies represented: **{per.study.nunique()} / 10**. Known cell-line identifiers represented: **{sum(table.cell_line != 'unassigned_cell_line')} / 54**, plus an explicitly labeled unassigned row where needed. Known matrix × cell-line records: **{sum(per.cell_line != 'unassigned_cell_line')} / 121**.\n\n'
 text+='Cell counts are cell records across source files, not a claim of deduplicated cells across studies. Target-condition groups are counted separately by source and condition; KO and CRISPRa results are not knockdown-efficacy claims.\n\n'
 text+='| Cell line | New expression matrices | Cell records | Cells passing QC | Studies |\n|---|---:|---:|---:|---|\n'
 for r in rows:text+=f"| {r['full_cell_line']} | {r['expression_matrices']} | {r['cell_records']:,} | {r['cells_kept']:,} | {r['study_names']} |\n"
 text+='\nSee `METHOD.md` for thresholds, applicability, and limitations. Per-source folders contain cell masks, guide QC and knockdown QC.\n'
 (OUT/'REPORT.md').write_text(text)
 print('COMPLETED',per.source.nunique(),'matrices;',len(table),'cell lines;',len(per),'matrix-cell-line rows;',int(per.cells.sum()),'cell records')
t=pd.read_parquet(ROOT/'lib/tahoe/metadata/cell_line_metadata.parquet');t=t[['cell_name','Cell_ID_DepMap','Cell_ID_Cellosaur','Organ']].drop_duplicates();t['expression_downloaded']=False;t['cell_qc_status']='not_applicable_metadata_only';t.to_csv(OUT/'tahoe_metadata_only_cell_lines.csv',index=False)
print('Tahoe metadata-only cell lines',t.cell_name.nunique())

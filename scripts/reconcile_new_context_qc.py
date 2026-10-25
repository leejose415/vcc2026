"""Reconcile QC masks with every original construct, including groups with zero retained cells."""
import json
from pathlib import Path
import numpy as np,pandas as pd
from qc_new_context import ROOT,OUT,BAD,dump
invalid=BAD|{'*'}
for p in sorted(OUT.glob('*/summary.json')):
 dest=p.parent;source=dest.name;summary=json.loads(p.read_text())
 if summary.get('reconciled'):continue
 o=pd.read_csv(dest/'cell_qc.csv.gz',keep_default_na=False,low_memory=False)
 if 'source_row' not in o:o.insert(0,'source_row',np.arange(len(o)))
 keys=['cell_line','condition','guide','target'];eligible=(~o.is_control)&(~o.target.isin(invalid))
 original=o.loc[eligible].groupby(keys,sort=False,dropna=False).agg(n_cells_original=('cell_id','size'),n_cells_after_cell_qc=('cell_keep','sum')).reset_index()
 try:g=pd.read_csv(dest/'guide_qc.csv',keep_default_na=False)
 except pd.errors.EmptyDataError:g=pd.DataFrame(columns=['cell_line','condition','guide_id','target_gene'])
 g=g.rename(columns={'guide_id':'guide','target_gene':'target'});g=g[~g.target.isin(invalid)]
 g=original.merge(g,on=keys,how='left',validate='one_to_one');missing=g.n_cells.isna() if 'n_cells' in g else pd.Series(True,index=g.index)
 if missing.any():
  assert (g.loc[missing,'n_cells_after_cell_qc']==0).all(),source
  g.loc[missing,'efficacy_status']='no_cells_after_cell_qc';g.loc[missing,'guide_keep']=False;g.loc[missing,'transcript_effect_pass']=False
  g.loc[missing,'modality']=summary['rows'][0]['modality'];g.loc[missing,'source']=source
 g['n_cells']=g.n_cells_after_cell_qc.astype(int);g['guide_keep']=g.guide_keep.fillna(False).astype(bool)
 ctrl=o[o.cell_keep&o.is_control].groupby(['cell_line','condition']).size()
 g['n_controls']=[int(ctrl.get((l,c),0)) for l,c in zip(g.cell_line,g.condition)]
 rows=[]
 for (line,condition,target),a in g.groupby(['cell_line','condition','target'],sort=False):
  kept=a[a.guide_keep];n=int(kept.n_cells.sum());lfc=pd.to_numeric(kept.ontarget_lfc,errors='coerce');cm=pd.to_numeric(kept.ctrl_cpm,errors='coerce');valid=np.isfinite(lfc)&np.isfinite(cm)&(kept.n_cells>0)
  eff=np.nan
  if valid.any():
   mean=np.average(np.exp2(lfc[valid])*(cm[valid]+.5)-.5,weights=kept.loc[valid,'n_cells']);ref=np.average(cm[valid],weights=kept.loc[valid,'n_cells']);eff=float(np.log2((mean+.5)/(ref+.5)))
  rows.append(dict(source=source,cell_line=line,condition=condition,target_gene=target,n_cells_original=int(a.n_cells_original.sum()),n_cells_after_cell_qc=int(a.n_cells.sum()),n_constructs=len(a),n_constructs_kept=int(a.guide_keep.sum()),n_cells_kept=n,measurable=bool(a.efficacy_status.eq('measurable').any()),efficiency=eff,knockdown_keep=n>=30,modality=summary['rows'][0]['modality'],unverifiable_constructs=int(a.efficacy_status.isin(['no_matched_controls','target_not_measured_or_non_single_gene']).sum())))
 kd=pd.DataFrame(rows)
 guidekeys=pd.MultiIndex.from_frame(g.loc[g.guide_keep,keys]);okguide=pd.MultiIndex.from_frame(o[keys]).isin(guidekeys)
 kepttargets=pd.MultiIndex.from_frame(kd.loc[kd.knockdown_keep,['cell_line','condition','target_gene']].rename(columns={'target_gene':'target'})) if len(kd) else pd.MultiIndex.from_tuples([],names=['cell_line','condition','target'])
 oktarget=pd.MultiIndex.from_frame(o[['cell_line','condition','target']]).isin(kepttargets)
 o['analysis_keep']=o.cell_keep & (o.is_control | (okguide & oktarget & eligible))
 o['perturbation_usable']=o.analysis_keep & ~o.is_control
 g.rename(columns={'guide':'guide_id','target':'target_gene'}).to_csv(dest/'guide_qc.csv',index=False);kd.to_csv(dest/'knockdown_qc.csv',index=False);o.to_csv(dest/'cell_qc.csv.gz',index=False)
 for row in summary['rows']:
  a=o[o.cell_line==row['cell_line']];k=kd[kd.cell_line==row['cell_line']] if len(kd) else kd
  row.update(targets_total=len(k),targets_kept=int(k.knockdown_keep.sum()) if len(k) else 0,analysis_cells=int(a.analysis_keep.sum()),analysis_perturbed_cells=int(a.perturbation_usable.sum()),unassigned_cells=int(a.target.isin(invalid).sum()))
 summary['reconciled']=True;summary['mask_note']='cell_keep is RNA-depth QC; analysis_keep additionally requires explicit control or usable assigned perturbation. source_row records original source row/column order.';dump(p,summary)
 print('RECONCILED',source,len(o),'cells',len(g),'construct-context rows',flush=True)

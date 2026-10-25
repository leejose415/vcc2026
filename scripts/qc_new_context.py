"""Read-only QC for newly downloaded context data; outputs masks/tables, never count copies.
Uses existing VCC UMI/MAD and on-target thresholds; stratifies biological contexts.
"""
import argparse,json,re,time,hashlib
from pathlib import Path
import numpy as np,pandas as pd,h5py,scipy.sparse as sp
from scipy.stats import mannwhitneyu, norm
from anndata.io import read_elem
from qc_replogle_umi import cutoffs
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results/qc/new_context'
BAD={'','nan','None','NA','unassigned','unknown'}
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
def mw_reference(v):
 values,counts=np.unique(v,return_counts=True);counts=counts.astype(float)
 return np.sort(v),values,counts,float(np.sum(counts**3-counts))
def mw_cached(x,reference,greater=False):
 sorted_y,unique_y,y_counts,tie_base=reference;n1=len(x);n2=len(sorted_y)
 left=np.searchsorted(sorted_y,x,side='left');right=np.searchsorted(sorted_y,x,side='right');u=float(np.sum((left+right)*.5))
 ux,nx=np.unique(x,return_counts=True);ix=np.searchsorted(unique_y,ux);ny=np.zeros(len(ux));ok=ix<len(unique_y);ok[ok]=unique_y[ix[ok]]==ux[ok];ny[ok]=y_counts[ix[ok]]
 ties=tie_base+np.sum((ny+nx)**3-(ny+nx)-(ny**3-ny));n=n1+n2
 variance=n1*n2/12*((n+1)-ties/(n*(n-1)))
 if variance<=0:return 1.0
 z=(u-n1*n2/2+(-.5 if greater else .5))/np.sqrt(variance)
 return float(norm.sf(z) if greater else norm.cdf(z))

def cell_qc(o):
 o=o.copy();o['cell_keep']=True;o['cell_reason']='';o['reference']='';cut=[]
 for context,idx in o.groupby(['cell_line','condition'],sort=False,observed=True).indices.items():
  g=o.iloc[idx];u=g.n_umi.to_numpy(float);x=np.log10(np.maximum(u,0)+1);ctrl=g.is_control.to_numpy(bool)
  valid=np.isfinite(x);ref=ctrl & valid;fallback=not ref.any()
  if fallback:ref=valid
  if not ref.any():o.loc[o.index[idx],'cell_keep']=False;o.loc[o.index[idx],'cell_reason']='invalid_counts';continue
  whole=cutoffs(x[ref]);bins=np.linspace(x[valid].min(),x[valid].max()+1e-9,81);hist,_=np.histogram(x[valid],bins=bins);share=hist[hist>0][0]/max(hist.max(),1);sharp=bool(share>=.01)
  low=np.zeros(len(g));high=np.zeros(len(g))
  for batch,ii in g.groupby('batch',sort=False,observed=True).indices.items():
   rr=ii[ref[ii]];c=cutoffs(x[rr]) if len(rr)>=50 else whole
   low[ii]=c['low_umi'];high[ii]=c['high_umi'];cut.append(dict(cell_line=context[0],condition=context[1],batch=batch,n_controls=int(ctrl[ii].sum()),fallback=len(rr)<50 or fallback,sharp_edge=sharp,**c))
  initial=(u<low)|(u>high)|~valid;rate=initial[ref].mean();flagged=set()
  for pert,ii in g.groupby('perturbation',sort=False,observed=True).indices.items():
   if not ctrl[ii].any() and len(ii)>=30 and initial[ii].mean()>max(2*rate,rate+.05):
    flagged.add(pert)
    if not sharp:low[ii]=cutoffs(x[ii][valid[ii]])['low_umi'] if valid[ii].any() else 200
  reason=np.select([~valid,u<200,u>high,(u<low)&(not sharp)],['invalid_counts','floor','high','low'],'')
  o.loc[o.index[idx],'cell_keep']=reason=='';o.loc[o.index[idx],'cell_reason']=reason
  o.loc[o.index[idx],'reference']='all_cells_no_controls' if fallback else 'context_controls'
 pd.DataFrame(cut).to_csv(OUT/o.source.iloc[0]/'cell_cutoffs.csv',index=False)
 return o

def run_qc(source,o,counts,genes,modality,notes=None):
 """counts: sparse cells x targeted genes only; total RNA UMIs supplied in o."""
 dest=OUT/source;dest.mkdir(parents=True,exist_ok=True);o=o.reset_index(drop=True);o['source']=source
 o=cell_qc(o);gene_col={str(g):i for i,g in enumerate(genes)};counts=counts.tocsc();rows=[];pooled=[]
 for context,idx in o.groupby(['cell_line','condition'],sort=False,observed=True).indices.items():
  selected=idx[o.cell_keep.to_numpy()[idx]];ctrl=selected[o.is_control.to_numpy()[selected]];cache={}
  for (guide,target),ix in o.iloc[selected].groupby(['guide','target'],sort=False,observed=True).groups.items():
   ix=np.asarray(list(ix),int)
   if o.loc[ix,'is_control'].all():continue
   if target in BAD:continue
   j=gene_col.get(str(target));n=len(ix);status='';lfc=p=cm=np.nan
   if not len(ctrl):status='no_matched_controls'
   elif j is None:status='target_not_measured_or_non_single_gene'
   else:
    if j not in cache:
     column=counts[:,j];control_values=column[ctrl].toarray().ravel()*1e6/o.n_umi.to_numpy()[ctrl];cache[j]=(column,control_values,mw_reference(control_values))
    column,control_values,reference=cache[j];values=column[ix].toarray().ravel()*1e6/o.n_umi.to_numpy()[ix];cm=float(control_values.mean());lfc=float(np.log2((values.mean()+.5)/(cm+.5)));p=mw_cached(values,reference,greater=modality=='CRISPRa')
    status='measurable' if cm>=5 else 'low_control_expression'
   passes=bool(status=='measurable' and ((lfc>=.5) if modality=='CRISPRa' else (lfc<=-.5)) and p<.05)
   # Transcript reduction is not a reliable functional KO criterion. Report it, don't discard KO cells on that basis.
   applicable=modality in ['CRISPRi','CRISPRa'];keep=bool(status not in ['no_matched_controls','target_not_measured_or_non_single_gene'] and (not applicable or status=='low_control_expression' or passes))
   rows.append(dict(source=source,cell_line=context[0],condition=context[1],guide_id=guide,target_gene=target,n_cells=n,n_controls=len(ctrl),ctrl_cpm=cm,ontarget_lfc=lfc,p_one_sided=p,transcript_effect_pass=passes,efficacy_status=status,efficacy_filter_applied=applicable,guide_keep=keep,modality=modality))
  # don't retain whole-context dense vectors across contexts
  del cache
 guides=pd.DataFrame(rows)
 if len(guides):
  for (line,condition,target),g in guides.groupby(['cell_line','condition','target_gene'],sort=False):
   kept=g[g.guide_keep];n=int(kept.n_cells.sum());meas=g.efficacy_status.eq('measurable').any()
   pooled.append(dict(source=source,cell_line=line,condition=condition,target_gene=target,n_constructs=len(g),n_constructs_kept=len(kept),n_cells=int(g.n_cells.sum()),n_cells_kept=n,measurable=bool(meas),knockdown_keep=bool(n>=30),modality=modality,interpretation='transcript_effect_diagnostic_only' if modality=='CRISPRko' else 'on_target_filter',unverifiable_constructs=int(g.efficacy_status.isin(['no_matched_controls','target_not_measured_or_non_single_gene']).sum())))
 guides.to_csv(dest/'guide_qc.csv',index=False);kd=pd.DataFrame(pooled);kd.to_csv(dest/'knockdown_qc.csv',index=False)
 o.to_csv(dest/'cell_qc.csv.gz',index=False)
 summary=[]
 for line,g in o.groupby('cell_line',sort=False):
  k=kd[kd.cell_line==line] if len(kd) else kd
  summary.append(dict(source=source,cell_line=line,cells=len(g),cells_kept=int(g.cell_keep.sum()),controls=int(g.is_control.sum()),controls_kept=int((g.is_control & g.cell_keep).sum()),conditions=g.condition.nunique(),targets_total=len(k),targets_kept=int(k.knockdown_keep.sum()) if len(k) else 0,modality=modality))
 dump(dest/'summary.json',dict(rows=summary,notes=notes or [],parameters=dict(umi_floor=200,low_mad=3,high_mad=4,mad_scale=1.4826,mad_floor=.05,min_batch_controls=50,min_cells=30,min_ctrl_cpm=5,ontarget_lfc=.5 if modality=='CRISPRa' else -.5,p_max=.05),outputs='cell masks and guide/target tables; original counts unchanged'))
 print(source,summary,flush=True);return summary

def sc_adapter(path):
 source=path.stem;dest=OUT/source;dest.mkdir(parents=True,exist_ok=True);notes=[]
 with h5py.File(path,'r') as f:
  obs=read_elem(f['obs']);var=read_elem(f['var']);symbols=var.index.astype(str).to_numpy();o=pd.DataFrame({'cell_id':obs.index.astype(str)},index=obs.index)
  kolf=source.startswith('KOLF');modality='CRISPRi' if kolf or source.startswith(('Adamson','Gasperini')) else ('CRISPRa' if source.startswith('Norman') else 'CRISPRko')
  o['cell_line']='KOLF2.1J' if kolf else ('melanoma_2686' if source.startswith('Frangieh') else obs.cell_line.astype(str).str.replace(' cells','',regex=False))
  o['condition']=obs['perturbation_2'].astype(str) if 'perturbation_2' in obs else 'none'
  batch=next((c for c in (['channel','batch'] if kolf else ['gemgroup','sample','replicate']) if c in obs),None);o['batch']=obs[batch].astype(str) if batch else 'not_recorded'
  pert=obs.perturbation.astype(str);o['perturbation']=pert;o['is_control']=pert.str.lower().eq('control')
  o['guide']=obs['guide_id'].astype(str) if 'guide_id' in obs else pert
  if kolf:o['target']=obs.gene_target.astype(str);o['is_control']=o.target.eq('NTC');o['guide']=pert
  elif source.startswith('Adamson'):
   control=pert.str.contains(r'neg_ctrl|Gal4|^6[23]\(mod\)',regex=True,case=False);o['is_control']=control
   o['target']=pert.str.split('_').str[0].replace({'PERK':'EIF2AK3','IRE1':'ERN1'})
   multi=pert.str.contains('ATF6_PERK|ATF6_IRE1|PERK_IRE1');o.loc[multi,'target']='combination:'+pert[multi]
   notes.append('Controls: published negative constructs 62/63(mod), Gal4, and neg_ctrl; PERK=EIF2AK3, IRE1=ERN1. Combination constructs excluded from single-gene efficiency.')
  elif source.startswith('Dixit') or source.startswith('DatlingerBock2017'):o['target']=obs.target.astype(str)
  elif source.startswith('DatlingerBock2021'):o['target']=pert.str.replace(r'_\d+$','',regex=True)
  else:o['target']=pert
  if source.startswith('Gasperini'):
   notes.append('Enhancer and multi-target constructs cannot be assessed as single-gene knockdowns. Unassigned cells are not controls.')
  o.loc[o.is_control,'target']='control'
  invalid=o.target.isin(BAD)|o.perturbation.isin(BAD);o.loc[invalid,'target']='unassigned'
  if source.startswith('Dixit'):
   # high-MOI target metadata can be scalar even when guide assignment is multi; avoid conflating constructs.
   mult=obs.get('nperts',pd.Series(1,index=obs.index)).fillna(0).to_numpy()>1;o.loc[mult & ~o.is_control,'target']='combination:'+o.loc[mult & ~o.is_control,'target']
  # Ensembl-indexed at-scale data need symbol -> ID lookup from local scPerturb var tables.
  if source.endswith('atscale'):
   lookup={}
   for other in sorted((ROOT/'lib/scperturb').glob('*.h5ad')):
    if other==path:continue
    with h5py.File(other,'r') as h:
     v=read_elem(h['var']);c=next((c for c in ['ensembl_id','gene_id','ensemble_id'] if c in v),None)
     if c:lookup.update(dict(zip(v[c].astype(str),v.index.astype(str))))
   symbols=np.array([lookup.get(s,s) for s in symbols])
  wanted=set(o.loc[~o.is_control,'target'])-BAD;cols=np.flatnonzero(np.isin(symbols,list(wanted)));X=f['layers/counts'] if kolf else f['X'];sums=np.zeros(len(o),float);ngenes=np.zeros(len(o),np.int32);mito=np.zeros(len(o),float);selected=[];selected_ids=[];bad=0
  enc=X.attrs.get('encoding-type','')
  if enc=='csc_matrix':
   ip=X['indptr'][:]
   for a in range(0,len(symbols),128):
    b=min(a+128,len(symbols));lo,hi=int(ip[a]),int(ip[b]);data=X['data'][lo:hi];ii=X['indices'][lo:hi];bad+=int(np.count_nonzero((data<0)|~np.isfinite(data)|(data!=np.rint(data))))
    sums+=np.bincount(ii,weights=data,minlength=len(o));ngenes+=np.bincount(ii[data>0],minlength=len(o)).astype(np.int32)
    block=sp.csc_matrix((data,ii,ip[a:b+1]-lo),shape=(len(o),b-a));mi=np.flatnonzero(np.char.startswith(symbols[a:b].astype(str),'MT-'))
    if len(mi):mito+=np.asarray(block[:,mi].sum(1)).ravel()
    cc=cols[(cols>=a)&(cols<b)]
    if len(cc):selected.append(block[:,cc-a]);selected_ids.extend(cc)
  elif enc=='csr_matrix':
   ip=X['indptr'][:]
   for a in range(0,len(o),5000):
    b=min(a+5000,len(o));lo,hi=int(ip[a]),int(ip[b]);data=X['data'][lo:hi];ii=X['indices'][lo:hi];bad+=int(np.count_nonzero((data<0)|~np.isfinite(data)|(data!=np.rint(data))))
    block=sp.csr_matrix((data,ii,ip[a:b+1]-lo),shape=(b-a,len(symbols)));sums[a:b]=np.asarray(block.sum(1)).ravel();ngenes[a:b]=np.asarray((block>0).sum(1)).ravel();mi=np.flatnonzero(np.char.startswith(symbols.astype(str),'MT-'));mito[a:b]=np.asarray(block[:,mi].sum(1)).ravel();selected.append(block[:,cols])
   selected_ids=cols.tolist()
  else:raise ValueError('Unsupported raw-count encoding '+str(enc))
  if bad:raise ValueError(f'{source}: {bad} non-integer/nonfinite/negative counts')
  C=(sp.vstack(selected).tocsc() if enc=='csr_matrix' else sp.hstack(selected).tocsc()) if selected else sp.csc_matrix((len(o),0))
  o['n_umi']=sums;o['n_genes']=ngenes;o['percent_mito']=np.divide(mito*100,sums,out=np.zeros(len(o)),where=sums>0)
  notes.append('Raw counts validated nonnegative/integer; RNA depth recomputed from raw matrix. Mitochondrial fraction reported, no new cutoff imposed.')
 return run_qc(source,o,C,symbols[selected_ids],modality,notes)

def main():
 p=argparse.ArgumentParser();p.add_argument('--source');args=p.parse_args();paths=sorted((ROOT/'lib/scperturb').glob('*.h5ad'))+list((ROOT/'lib/ipsc/kolf').glob('*.h5ad'))
 for path in paths:
  if args.source and args.source not in path.stem:continue
  if (OUT/path.stem/'summary.json').exists():continue
  print('START',path,flush=True);sc_adapter(path)
if __name__=='__main__':main()

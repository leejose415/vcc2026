import gzip,subprocess,re,json,sys
from pathlib import Path
import numpy as np,pandas as pd,scipy.sparse as sp
from qc_new_context import ROOT,OUT,run_qc,BAD
for screen in (sys.argv[1:] or ['GenomeWideScreen_FitnessGenes','GenomeWideScreen_NonFitnessGenes','TargetedScreen']):
 source='Feng_'+screen;tmp=OUT/(source+'_extract');tmp.mkdir(parents=True,exist_ok=True)
 if (OUT/source/'summary.json').exists():continue
 m=pd.read_csv(ROOT/'lib/ipsc/feng_counts'/f'{screen}_Cell-Metadata.tsv.gz',sep='\t',dtype=str)
 # Published RNA columns omit MP- or PC-/day prefixes, while metadata retains them.
 norm=m.Cell_ID.str.replace(r'^MP-','',regex=True).str.replace(r'^PC-','',regex=True)
 if screen=='TargetedScreen':norm=norm.str.replace(r'-D3_','_',regex=True)
 assert not norm.duplicated().any();m.index=norm
 rna=ROOT/'lib/ipsc/feng_counts'/f'{screen}_RNA-UMI-Counts.csv.gz'
 with gzip.open(rna,'rt') as f:columns=f.readline().strip().split(',')[1:]
 assert len(columns)==len(set(columns));missing=set(columns)-set(m.index)
 if len(missing)>1:raise ValueError(f'{len(missing)} RNA columns missing metadata; examples {list(missing)[:5]}')
 mapping_audit={'rna_columns':len(columns),'metadata_rows':len(m),'rna_columns_missing_metadata':sorted(missing),'metadata_rows_absent_from_rna':len(set(m.index)-set(columns))}
 (tmp/'barcode_mapping_audit.json').write_text(json.dumps(mapping_audit,indent=2))
 m=m.reindex(columns)
 for cell in missing:m.loc[cell,['Cell_ID','Batch','Guide_Call','Cell_Line']]=[cell,'metadata_not_available','unassigned','unassigned_cell_line']
 print(screen,len(m),'aligned RNA cells; missing metadata',len(missing),flush=True)
 guide=m.Guide_Call;target=guide.str.rsplit('_',n=1).str[0];control=target.eq('NonTarget');target=target.mask(control,'control');target=target.mask(guide.eq('unassigned'),'unassigned')
 o=pd.DataFrame({'cell_id':m.Cell_ID.values,'cell_line':m.Cell_Line.values,'condition':m.Batch.str.extract(r'(Day\d+)',expand=False).fillna('unspecified').values,'batch':m.Batch.values,'guide':guide.values,'target':target.values,'perturbation':target.values,'is_control':control.values})
 if not (tmp/'complete.txt').exists():
  with open(tmp/'setup.tsv','w') as f:
   for t,c in zip(target,control):f.write(str(t)+'\t'+str(int(c))+'\n')
  subprocess.run(['/tmp/vcc-stream-feng-qc',str(rna),str(tmp/'setup.tsv'),str(tmp)],check=True)
 genes=(tmp/'genes.txt').read_text().splitlines();measured=np.fromfile(tmp/'measured.bin',dtype=np.uint32).astype(bool);trip=np.fromfile(tmp/'evidence.bin',dtype=[('row','<u4'),('col','<u4'),('value','<u4')]);C=sp.coo_matrix((trip['value'],(trip['row'],trip['col'])),shape=(len(o),len(genes))).tocsc();del trip
 o['n_umi']=np.fromfile(tmp/'totals.bin',dtype=np.float64);o['n_genes']=np.fromfile(tmp/'ngenes.bin',dtype=np.uint32);mt=np.fromfile(tmp/'mito.bin',dtype=np.float64);o['percent_mito']=np.divide(mt*100,o.n_umi.to_numpy(),out=np.zeros(len(o)),where=o.n_umi.to_numpy()>0)
 # Unmeasured target columns must remain distinguishable from genuinely measured all-zero genes.
 run_qc(source,o,C[:,measured],np.asarray(genes)[measured],'CRISPRi',['Raw gene-by-cell counts streamed without full matrix conversion. RNA/metadata barcode mapping validated. Explicit NonTarget controls; unassigned cells are not controls. Source LFC and guide-count files are supporting measurements, not additional independent cell datasets.'])

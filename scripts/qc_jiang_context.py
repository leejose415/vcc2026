import subprocess,sys,json,shutil
from pathlib import Path
import pandas as pd,numpy as np,scipy.sparse as sp
from qc_new_context import ROOT,OUT,run_qc
for p in sorted((ROOT/'lib/jiang').glob('Seurat_object_*Perturb_seq.rds'),key=lambda p:p.stat().st_size):
 source='Jiang_'+p.name.split('_')[2]
 if (OUT/source/'summary.json').exists():continue
 tmp=OUT/(source+'_extract');tmp.mkdir(exist_ok=True)
 if not (tmp/'shape.txt').exists():subprocess.run(['Rscript',str(ROOT/'scripts/extract_jiang_qc.R'),str(p),str(tmp)],check=True)
 m=pd.read_csv(tmp/'metadata.csv.gz',index_col=0,low_memory=False);genes=(tmp/'genes.txt').read_text().splitlines();shape=tuple(map(int,(tmp/'shape.txt').read_text().split()))
 C=sp.csc_matrix((np.fromfile(tmp/'data.bin',dtype='<f8'),np.fromfile(tmp/'indices.bin',dtype='<i4'),np.fromfile(tmp/'indptr.bin',dtype='<i4')),shape=shape)
 o=pd.DataFrame({'cell_id':m.index.astype(str),'cell_line':m.cell_type.astype(str).values,'condition':m.pathway.astype(str).values,'batch':m.Batch_info.astype(str).values,'guide':m.guide.astype(str).values,'target':m.gene.astype(str).values,'perturbation':m.gene.astype(str).values,'n_umi':m.qc_n_umi.values,'n_genes':m.qc_n_genes.values,'percent_mito':m['percent.mito'].values})
 o['is_control']=o.target.eq('NT');o.loc[o.is_control,'target']='control'
 run_qc(source,o,C,genes,'CRISPRi',['Source Seurat RNA counts; non-targeting control=NT; separate cell lines and stimulated pathways; original counts unchanged.'])
 # Keep metadata for reproducibility, remove only temporary extracted numeric buffers.
 for f in ['data.bin','indices.bin','indptr.bin']:(tmp/f).unlink()

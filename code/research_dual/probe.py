"""Ten-knockdown G2/G3 feasibility check using only frozen predictions and controls."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='2'
import sys,pickle,json,argparse,time
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'code/vcc_2026_3c'))
import vcc3_generate as gen
from generator import DualGenerator

def inputs():
    line,*_=pickle.load(open(ROOT/'lib/research_mse/h1_generator_controls.pkl','rb'))
    pdir=ROOT/'lib/vcc_2026_4/runs/dev_nopie/h1_pred_g3'
    z=np.load(pdir/'raw_predictions.npz')
    ts=np.load(pdir/'targets.npy').astype(str)
    genes=pd.read_csv(ROOT/'data/gene_names.csv').gene_name.tolist();gp={g:i for i,g in enumerate(genes)}
    settings=json.load(open(pdir.parent/'h1_bench_g3/selected_generator.json'))['settings']
    settings['detect_floor']=False;c=gen.Cfg(argparse.Namespace(**settings))
    D=np.zeros((len(ts),len(genes)));D[:,z['gene_cols']]=z['delta']
    D=gen.calibrate(D,line,ts,c)
    kd=gen.kd_residual(gen.prediction_ontarget(pdir,ts),c)
    pick=gen.pick_cells(line,np.zeros(len(genes)),c,np.random.default_rng(1729))[:400]
    # Randomly assign templates to depth blocks; the moments, not sampled depths, determine the covariance.
    template=line.cpm_cells[pick].toarray()/1e6
    engine=DualGenerator(line.cpm/1e6,line.share,line.cpm_var/1e12,template)
    return line,ts,genes,gp,D,kd,engine

if __name__=='__main__':
    line,ts,genes,gp,D,kd,engine=inputs();records=[]
    out=ROOT/'results/research_dual';out.mkdir(exist_ok=True,parents=True)
    n=int(sys.argv[1]) if len(sys.argv)>1 else 10
    for i in range(n):
        f,q,_=gen.project_factors(line.share,D[i],gp[ts[i]],kd[i])
        for m in [1.,.5,.25]:
            t=time.time();X,stats=engine.generate(f,D[i]!=0,m,pooled=q)
            stats.update(target=ts[i],seconds=time.time()-t)
            records.append(stats);print(stats,flush=True)
            if i==0:np.savez_compressed(out/f'probe_first_m{m}.npz',X=X,f=f,q=q,called=D[i]!=0)
            pd.DataFrame(records).to_csv(out/'g2_probe.csv',index=False)

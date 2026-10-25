"""Isolated G6/G7 fits. Reads Model 7 tables; writes only results/research_dual/holdouts."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='2'
import json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path.home()/'vcc'
sys.path.insert(0,str(ROOT/'analysis/model7'))
import m7_fit as original
OUT=ROOT/'results/research_dual/holdouts'
OUT.mkdir(parents=True,exist_ok=True)
CHALLENGE=['K562','HCT116','HEK293T','CD4_Rest','CD4_Stim8hr','CD4_Stim48hr']

def eligible(source,destination):
    return source!=destination and not (source.startswith('CD4_') and destination.startswith('CD4_'))

def main():
    genes=original.genes_list();gp={g:i for i,g in enumerate(genes)}
    tabs={n:original.Tab(n) for n in original.TRAIN+['H1']}
    cache={}
    def panel_stats(source,destination,matched):
        key=(source,destination,matched)
        if key in cache:return cache[key]
        s=tabs[source]
        if not matched:return s.R,s.varR
        wanted=set(tabs[destination].row)
        rows=np.array([i for t,i in s.row.items() if t in wanted])
        if not len(rows):
            ans=(np.full(len(genes),np.nan),np.full(len(genes),np.nan))
        else:
            D=s.D[rows].astype(float);v=s.vD[rows].astype(float)
            valid=np.isfinite(D)&np.isfinite(v)&s.measured[None]
            K=valid.sum(0)
            R=np.divide(np.where(valid,D,0).sum(0),K,out=np.full(len(genes),np.nan),where=K>0)
            independent=np.maximum(v-s.var0[None],0)
            vr=np.divide(np.where(valid,independent,0).sum(0),K*K,
                         out=np.full(len(genes),np.nan),where=K>0)+s.var0
            ans=R,vr
        cache[key]=ans;return ans
    def shared(dest,sources,matched):
        num=np.zeros(len(genes));den=num.copy()
        for source in sources:
            R,v=panel_stats(source,dest,matched)
            ok=np.isfinite(R)&np.isfinite(v)&(v>0)
            w=np.divide(1,v,out=np.zeros_like(v),where=ok)
            num+=w*np.nan_to_num(R);den+=w
        return np.divide(num,den,out=np.zeros_like(num),where=den>0)
    def features(dest,rows,sources,matched):
        d=tabs[dest];allowed=[s for s in sources if eligible(s,dest)]
        F=np.zeros((1+len(sources),len(rows),len(genes)),np.float32)
        F[0]=shared(dest,allowed,matched)[None]
        for si,s in enumerate(sources):
            if s not in allowed:continue
            positions=np.array([tabs[s].row.get(t,-1) for t in d.tg[rows]])
            ok=positions>=0
            if ok.any():
                R,_=panel_stats(s,dest,matched)
                F[1+si,ok]=np.nan_to_num(tabs[s].D[positions[ok]].astype(np.float32)-R.astype(np.float32))
        return F
    records=[]
    for tag,sources,matched in [('all_family',original.TRAIN,False),('challenge_family',CHALLENGE,False),
                               ('challenge_matched_se',CHALLENGE,True)]:
        done=OUT/f'{tag}.npz'
        if done.exists():continue
        nf=1+len(sources);acc=dict(G=np.zeros((7,nf,nf)),Xy=np.zeros((7,nf)),yy=np.zeros(7),n=np.zeros(7))
        # Same destination set in all comparisons; only the available sources change.
        for dest in original.TRAIN:
            started=time.time();d=tabs[dest];q=original.groups_of(d.S0,d.measured)
            rows=np.flatnonzero(d.w>0)
            for start in range(0,len(rows),96):
                rr=rows[start:start+96];y=d.D[rr].astype(np.float32)
                mask=np.isfinite(y)&d.measured[None]
                for j,r in enumerate(rr):
                    if d.tg[r] in gp:mask[j,gp[d.tg[r]]]=False
                sw=np.sqrt(d.w[rr]).astype(np.float32)[:,None]
                F=features(dest,rr,sources,matched)*sw[None]*mask[None]
                y=np.where(mask,y,0)*sw
                for k in range(7):
                    keep=q==k;A=F[:,:,keep].reshape(nf,-1).astype(float);b=y[:,keep].ravel().astype(float)
                    acc['G'][k]+=A@A.T;acc['Xy'][k]+=A@b;acc['yy'][k]+=b@b;acc['n'][k]+=mask[:,keep].sum()
            print(tag,dest,'seconds',round(time.time()-started,1),flush=True)
        coef,info=original.solve(acc)
        d=tabs['H1'];q=original.groups_of(d.S0,d.measured)
        P=np.empty(d.D.shape,np.float32)
        for start in range(0,len(d.tg),96):
            rr=np.arange(start,min(start+96,len(d.tg)))
            F=features('H1',rr,sources,matched)
            P[rr]=np.einsum('frg,gf->rg',F,coef[q])
        np.savez_compressed(done,tg=d.tg.astype(str),delta=P,coef=coef,sources=np.array(sources),
                            S0=d.S0,real_delta=d.D,var_delta=d.vD,control_var=d.var0,w=d.w,measured=d.measured)
        (OUT/f'{tag}.json').write_text(json.dumps(dict(tag=tag,sources=sources,matched=matched,fit=info,
                a_by_group=coef[:,0].tolist(),cd4_family_heldout=True),indent=2))
        # Shared-response diagnostic: fit only one oracle multiplier for diagnosis, never a prediction arm.
        rh=shared('H1',sources,matched);DD=d.D.astype(float);mask=np.isfinite(DD)&d.measured[None]&(d.w[:,None]>0)
        for i,t in enumerate(d.tg):
            if t in gp:mask[i,gp[t]]=False
        vv=d.vD.astype(float);den=np.where(mask,DD**2-vv,0).sum()
        a=np.where(mask,DD*rh,0).sum()/np.where(mask,rh[None]**2,0).sum()
        def raw(x):return float(np.where(mask,(x-DD)**2-(vv-d.var0),0).sum()/den)
        records.append(dict(tag=tag,shared_a_oracle=float(a),shared_raw=raw(rh),shared_raw_oracle_a=raw(a*rh),
                            fit_a_by_group=coef[:,0].tolist()))
        (OUT/'shared_diagnostics.json').write_text(json.dumps(records,indent=2))
        print('FINISHED',tag,records[-1],flush=True)

if __name__=='__main__':main()

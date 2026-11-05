from pathlib import Path
import os,time,json
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control
B=Path(__file__).resolve().parents[1];P=B;O=B;R=B/'results'
m=pd.read_csv(P/'data/baseline_master.csv');d=m[(m.Months_from_Baseline==0)&m.ALPS_Mean.notna()&m.WholeBrain_SUVR_AAL3.notna()].copy();roi=[c for c in d if c.startswith('AAL3_') and d[c].std()>1e-10];assert len(roi)==134
extra=['Education_years','log1p_WMH','SBP_preindex','Hypertension_documented','Diabetes_documented','Dyslipidemia_documented','History_source']
d=d.dropna(subset=extra).copy()
ii,jj=np.triu_indices(len(roi),1);cut=d.loc[d.Group.isin(['MCI','AD']),'ALPS_Mean'].median();rng=np.random.default_rng(20261020)
def corr(x):
 mask=np.isfinite(x).astype(float);a=np.nan_to_num(x);n=mask.T@mask;s=a.T@mask;s2=(a*a).T@mask;c=a.T@a-s*s.T/n;v=s2-s*s/n
 return np.clip((c/np.sqrt(v*v.T))[ii,jj],-.999999999,.999999999)
def run(name,x,hi,groups=None):
 start=time.time();obs=np.arctanh(corr(x[hi]))-np.arctanh(corr(x[~hi]));counts=np.zeros(len(ii),dtype=int);Bperm=9999
 for b in range(Bperm):
  h=hi.copy()
  if groups is None:h=rng.permutation(h)
  else:
   for g in np.unique(groups):idx=np.flatnonzero(groups==g);h[idx]=rng.permutation(h[idx])
  z=np.arctanh(corr(x[h]))-np.arctanh(corr(x[~h]));counts+=np.abs(z)>=np.abs(obs)
  if (b+1)%2000==0:print(name,'permutation',b+1,'seconds',round(time.time()-start),flush=True)
 p=(1+counts)/(Bperm+1);q=false_discovery_control(p)
 a=np.flatnonzero(hi);bb=np.flatnonzero(~hi);boots=[]
 for b in range(500):boots.append(np.arctanh(corr(x[rng.choice(a,len(a),replace=True)]))-np.arctanh(corr(x[rng.choice(bb,len(bb),replace=True)])))
 low,high=np.quantile(np.asarray(boots),[.025,.975],axis=0)
 result=pd.DataFrame({'ROI1':np.array(roi)[ii],'ROI2':np.array(roi)[jj],'Delta_Fisher_z':obs,'Bootstrap_CI_low':low,'Bootstrap_CI_high':high,'Permutation_P':p,'BH_q_all_8911':q,'Permutations':Bperm,'Bootstrap_draws':500})
 result.sort_values('Permutation_P').to_csv(R/f'{name}_8911_participant_permutation.csv',index=False)
 return {'Analysis':name,'N':len(x),'N_high':len(a),'N_low':len(bb),'Edges':len(q),'BH_significant':int((q<.05).sum()),'Min_P':float(p.min()),'Min_q':float(q.min()),'seconds':time.time()-start}
z=d[d.Group.isin(['MCI','AD'])];res=[run('Vascular_unadjusted',z[roi].to_numpy(),(z.ALPS_Mean>=cut).to_numpy())]
cc=d.dropna(subset=roi+['Age','Sex','FA_min','Group']).copy();sel=pd.read_csv(P/'results/Vascular_S12_equalised_subject_selection.csv');assert list(sel.Subject)==list(cc.Subject)
design=pd.DataFrame({'Intercept':1,'Age':cc.Age,'Sex_F':(cc.Sex=='F').astype(int),'FA_min':cc.FA_min,'MCI':(cc.Group=='MCI').astype(int),'AD':(cc.Group=='AD').astype(int)})
design=pd.concat([design.reset_index(drop=True),cc[extra[:-1]].reset_index(drop=True),pd.get_dummies(cc.History_source,drop_first=True,dtype=float).reset_index(drop=True)],axis=1)
xx=cc[roi].to_numpy();des=design.to_numpy(dtype=float);resid=xx-des@np.linalg.lstsq(des,xx,rcond=None)[0];use=(sel.In_High_matched|sel.In_Low_matched).to_numpy()
res.append(run('Vascular_diagnosis_matched_residualised',resid[use],sel.loc[use,'In_High_matched'].to_numpy(),cc.loc[use,'Group'].to_numpy()))
pd.DataFrame(res).to_csv(R/'Vascular_Edge_permutation_summary.csv',index=False);print(res,flush=True)

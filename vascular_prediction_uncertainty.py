from pathlib import Path
import pandas as pd,numpy as np,json
O=Path(__file__).resolve().parents[1];P=O;R=O/'results'
x=pd.read_csv(R/'Vascular_Nested_OOF_predictions_PRIVATE.csv');m=pd.read_csv(P/'data/baseline_master.csv').set_index('Subject_ID');ids=sorted(x.Subject_ID.unique());idx={v:i for i,v in enumerate(ids)};rng=np.random.default_rng(20261020)
groups=[np.flatnonzero(m.loc[ids,'Group'].to_numpy()==g) for g in ['CN','MCI','AD']];draws=np.zeros((2000,len(ids)))
for b in range(2000):
 ix=np.concatenate([rng.choice(z,len(z),replace=True) for z in groups]);draws[b]=np.bincount(ix,minlength=len(ids))
res=[];drawsave=[];pooled=[]
for alg in ['RandomForest','XGBoost']:
 vals={}
 for fs in ['Clinical','Plus_ALPS']:
  z=x[(x.Algorithm==alg)&(x.Features==fs)];R2=[];BS=[];rm=[]
  for _,a in z.groupby(['Repeat','Fold']):
   ix=np.array([idx[v] for v in a.Subject_ID]);y=a.Observed.to_numpy();p=a.Predicted.to_numpy();w=draws[:,ix];s=w.sum(1);mean=(w@y)/s;ss=(w@(y*y))-s*mean*mean;rss=w@((y-p)**2)
   BS.append(1-rss/ss);R2.append(1-np.sum((y-p)**2)/np.sum((y-y.mean())**2));rm.append(np.sqrt(rss/s))
  bv=np.mean(BS,axis=0);vals[fs]=bv
  for metric,point,boot in [('Mean_test_fold_R2',np.mean(R2),bv),('Mean_test_fold_RMSE',z.groupby(['Repeat','Fold']).apply(lambda a:np.sqrt(np.mean((a.Observed-a.Predicted)**2)),include_groups=False).mean(),np.mean(rm,axis=0))]:
   lo,hi=np.quantile(boot,[.025,.975]);res.append({'Algorithm':alg,'Features':fs,'Metric':metric,'N':len(ids),'Estimate':point,'CI_low':lo,'CI_high':hi,'Bootstrap':2000,'Uncertainty':'Participant resampling conditional on saved outer-fold predictions; tuning/fitting not rerun'})
  a=z.groupby('Subject_ID').agg(Observed=('Observed','first'),Predicted=('Predicted','mean'));y=a.Observed.to_numpy();p=a.Predicted.to_numpy();coef=np.linalg.lstsq(np.column_stack([np.ones(len(p)),p]),y,rcond=None)[0];pooled.append({'Algorithm':alg,'Features':fs,'N':len(ids),'Pooled_R2_averaged_predictions':1-np.sum((y-p)**2)/np.sum((y-y.mean())**2),'RMSE_averaged_predictions':np.sqrt(np.mean((y-p)**2)),'Calibration_intercept':coef[0],'Calibration_slope':coef[1]})
 delta=vals['Plus_ALPS']-vals['Clinical'];point=np.mean(x[(x.Algorithm==alg)&(x.Features=='Plus_ALPS')].groupby(['Repeat','Fold']).apply(lambda a:1-np.sum((a.Observed-a.Predicted)**2)/np.sum((a.Observed-a.Observed.mean())**2),include_groups=False))-np.mean(x[(x.Algorithm==alg)&(x.Features=='Clinical')].groupby(['Repeat','Fold']).apply(lambda a:1-np.sum((a.Observed-a.Predicted)**2)/np.sum((a.Observed-a.Observed.mean())**2),include_groups=False))
 lo,hi=np.quantile(delta,[.025,.975]);res.append({'Algorithm':alg,'Features':'Plus_ALPS_minus_Clinical','Metric':'Delta_R2','N':len(ids),'Estimate':point,'CI_low':lo,'CI_high':hi,'Bootstrap':2000,'Uncertainty':'Participant resampling conditional on saved outer-fold predictions; tuning/fitting not rerun'})
pd.DataFrame(res).to_csv(R/'Vascular_Nested_performance_conditional_intervals.csv',index=False);pd.DataFrame(pooled).to_csv(R/'Vascular_Nested_pooled_metrics_and_calibration.csv',index=False)
print(pd.DataFrame(res)[['Algorithm','Features','Metric','Estimate','CI_low','CI_high']].to_string(index=False))

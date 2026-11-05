from pathlib import Path
import pandas as pd, numpy as np, json
B=Path(__file__).resolve().parents[1]; O=B
P=B
m=pd.read_csv(P/'data/baseline_master.csv'); dx=pd.read_csv(B/'data/DXSUM_26May2026.csv',low_memory=False)
dx['date']=pd.to_datetime(dx.EXAMDATE,errors='coerce');dx['code']=pd.to_numeric(dx.DIAGNOSIS,errors='coerce')
dx=dx[dx.date.notna()&dx.code.isin([1,2,3])].copy()
rr=[];summary=[]
for g in ['CN','MCI']:
 base=m[m.Group.eq(g)]; endpoints=[2,3] if g=='CN' else [3]
 subset=dx[dx.RID.isin(base.RID)]
 assert subset.groupby(['RID','date']).code.nunique().max()==1
 subset=subset.drop_duplicates(['RID','date']).sort_values('date')
 for _,r in base.iterrows():
  z=subset[subset.RID.eq(r.RID)];date=pd.Timestamp(r.DTI_DATE);pre=z[z.date.le(date)];post=z[z.date.gt(date)]
  last=pre.iloc[-1] if len(pre) else None
  excluded=last is not None and last.code in endpoints
  decision='latest endpoint' if excluded else 'no follow-up' if not len(post) else 'included'
  ev=post[post.code.isin(endpoints)];end=ev.date.min() if len(ev) else post.date.max()
  rr.append(dict(RID=r.RID,Group=g,Decision=decision,Any_prior_endpoint=pre.code.isin(endpoints).any(),Latest_prior_endpoint=excluded,Event=int(len(ev)>0) if decision=='included' else np.nan,Time=(end-date).days/365.25 if decision=='included' else np.nan))
 a=pd.DataFrame(rr);a=a[a.Group.eq(g)];inc=a[a.Decision.eq('included')]
 old=pd.read_csv(P/'data'/f'survival_{g}.csv').set_index('RID');new=inc.set_index('RID')
 assert set(new.index)==set(old.index),(g,'risk membership mismatch')
 assert np.allclose(new.loc[old.index,'Time'],old.time_years)
 assert np.array_equal(new.loc[old.index,'Event'],old.conversion)
 summary.append(dict(Group=g,Source_N=len(base),Excluded_latest_preindex_endpoint=int(a.Latest_prior_endpoint.sum()),Excluded_no_followup=int(a.Decision.eq('no follow-up').sum()),Retained_N=len(inc),Events=int(inc.Event.sum()),Any_prior_endpoint_N=int(a.Any_prior_endpoint.sum()),Any_vs_latest_discordant=int(a.Any_prior_endpoint.ne(a.Latest_prior_endpoint).sum()),Membership_events_times_match=True))
pd.DataFrame(rr).to_csv(O/'data/private/Survival_rule_verification_PRIVATE.csv',index=False)
pd.DataFrame(summary).to_csv(O/'results/Survival_eligibility_rule_verification.csv',index=False)
paths=pd.read_csv(O/'results/Complete_path_coefficients_and_intervals.csv')
paths['Reporting_role']=np.where(paths.Path.eq('ab'),'Primary exploratory indirect inference: pointwise bootstrap CI','Auxiliary OLS model summary; c/c_prime primary uncertainty is the bootstrap CI in S3')
paths['Interval_method']=paths.Interval_method.str.replace('P not inferred','P not reported; inference uses the interval',regex=False)
paths.to_csv(O/'results/Complete_path_coefficients_and_intervals.csv',index=False)
print(pd.DataFrame(summary).to_string(index=False))

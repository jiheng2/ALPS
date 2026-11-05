from pathlib import Path
import pandas as pd,numpy as np
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
B=Path(__file__).resolve().parents[1];O=B
d=pd.read_csv(B/'data/structural_PET216.csv');c=pd.read_csv(B/'data/regional_CSF_long.csv');cw=c.pivot(index='FDG_VISIT',columns='ROI',values='roi_CSF_mean').add_prefix('CSF__');d=d.merge(cw,left_on='PET_visit_id',right_index=True,how='left')
rows=[]
for roi in [c for c in d if c.startswith('AAL3_')]:
 for fa in ['FA_mean','FA_min']:
  csf='CSF__'+roi;z=d.dropna(subset=[roi,csf,fa,'ALPS_Mean','Age','Sex','Group','logLV','logICV','DTI_PET_gap_days']).copy()
  if len(z)<30 or z[roi].std()<1e-10:continue
  z['outcome']=z[roi];z['CSF']=z[csf]
  for model,add in [('Base',''),('LV_ICV','+logLV+logICV'),('LV_ICV_CSF','+logLV+logICV+CSF')]:
   model_input=smf.ols('outcome~ALPS_Mean+Age+C(Sex)+C(Group)+DTI_PET_gap_days+'+fa+add,z)
   if np.linalg.matrix_rank(model_input.exog)<model_input.exog.shape[1]:continue
   f=model_input.fit(cov_type='HC3');ci=f.conf_int().loc['ALPS_Mean'];sd=z.ALPS_Mean.std()/z.outcome.std();rows.append(dict(ROI=roi,FA=fa,Model=model,N=len(z),Beta=f.params.ALPS_Mean,CI_low=ci.iloc[0],CI_high=ci.iloc[1],Standardized_beta=f.params.ALPS_Mean*sd,Standardized_CI_low=ci.iloc[0]*sd,Standardized_CI_high=ci.iloc[1]*sd,P=f.pvalues.ALPS_Mean))
z=pd.DataFrame(rows);z['q_BH140']=z.groupby(['FA','Model']).P.transform(lambda p:multipletests(np.r_[p,np.ones(140-len(p))],method='fdr_bh')[1][:len(p)]);z.to_csv(O/'results/FA_swap_with_regional_CSF.csv',index=False);print(z.groupby(['FA','Model']).q_BH140.agg(['size',lambda p:(p<.05).sum()]).to_string())

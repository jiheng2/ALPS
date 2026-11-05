from pathlib import Path
import json, itertools, re
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.outliers_influence import variance_inflation_factor

B=Path(__file__).resolve().parents[1];O=B
P=B; R=O/'results'; A=O/'logs'
m=pd.read_csv(P/'data/baseline_master.csv'); an=pd.read_csv(O/'data/private/anatomy_linked_PRIVATE.csv')
def save(d,n): d.to_csv(R/n,index=False,encoding='utf-8-sig')
rows=[]
for h in ['Mean','Left','Right']:
    for g,k in itertools.combinations(['CN','MCI','AD'],2):
        x=m.loc[m.Group==g,'ALPS_'+h]; y=m.loc[m.Group==k,'ALPS_'+h]
        pooled=np.sqrt(((len(x)-1)*x.var(ddof=1)+(len(y)-1)*y.var(ddof=1))/(len(x)+len(y)-2))
        t=stats.ttest_ind(x,y,equal_var=False)
        rows.append(dict(Hemisphere=h,Group1=g,Group2=k,N1=len(x),N2=len(y),Mean1=x.mean(),Mean2=y.mean(),Cohens_d=(x.mean()-y.mean())/pooled,P=t.pvalue))
d=pd.DataFrame(rows);d['q_three_pairs']=d.groupby('Hemisphere').P.transform(lambda z:multipletests(z,method='fdr_bh')[1]); save(d,'Baseline_pairwise_Welch_Cohens_d.csv')
# Descriptive whole-brain FDG contrasts: restore original B with exact, named tests.
pet=m[(m.Months_from_Baseline==0)&m.WholeBrain_SUVR_AAL3.notna()].copy();assert len(pet)==pd.read_csv(P/"data/structural_PET216.csv").shape[0]
rows=[]
for g,k in itertools.combinations(['CN','MCI','AD'],2):
    x=pet.loc[pet.Group==g,'WholeBrain_SUVR_AAL3']; y=pet.loc[pet.Group==k,'WholeBrain_SUVR_AAL3']; t=stats.ttest_ind(x,y,equal_var=False)
    rows.append(dict(Group1=g,Group2=k,N1=len(x),N2=len(y),Mean1=x.mean(),Mean2=y.mean(),P=t.pvalue))
d=pd.DataFrame(rows);d['q']=multipletests(d.P,method='fdr_bh')[1]; save(d,'Whole_brain_FDG_group_contrasts.csv')
# Explicit target membership and available-ROI weights. This audit does not imply prospective target selection.
roi=[c for c in m if c.startswith('AAL3_')];target=[c for c in roi if re.search('Caudate|Putamen|Pallidum|Thal|N_Acc',c)]
assert len(target)==38
save(pd.DataFrame({'ROI':target,'Candidate_weight':1/38,'Rule':'Unweighted mean of available finite ROI SUVR values; renormalized per participant'}),'Prediction_target_all38_ROIs.csv')
nvalid=pet[target].notna().sum(axis=1)
save(pd.DataFrame({'Valid_ROIs':nvalid.value_counts().index,'Participants':nvalid.value_counts().values}),'Prediction_target_availability.csv')
# Index scans versus first PASS date; fixed Group originates from the selected-index master record.
long=pd.read_csv(P/'data/longitudinal_with_site_protocol_proxy.csv');long['date']=pd.to_datetime(long.Acq_Date,errors='coerce')
lp=long[long.QC_Status.eq('PASS')&long.date.notna()]; earliest=lp.groupby('Subject').date.min()
idx=pd.to_datetime(m.DTI_DATE); earliest=m.Subject.map(earliest)
save(pd.DataFrame({'Subject':m.Subject,'Index_date':idx,'First_PASS_DTI_date':earliest,'Difference_days':(idx-earliest).dt.days}),'Index_first_available_date_PRIVATE.csv')
# Summarize only date links actually retained in prepared datasets.
tim=[]
def ts(label,v,rule):
    v=pd.Series(v).dropna();tim.append(dict(Modality=label,N=len(v),Signed_median=v.median(),Signed_Q1=v.quantile(.25),Signed_Q3=v.quantile(.75),Signed_min=v.min(),Signed_max=v.max(),Absolute_median=v.abs().median(),Absolute_Q1=v.abs().quantile(.25),Absolute_Q3=v.abs().quantile(.75),Rule=rule))
ts('FDG-PET',(pd.to_datetime(pet.FDG_Date)-pd.to_datetime(pet.DTI_DATE)).dt.days,'PET minus DTI, signed; eligible absolute interval <=183 d')
ts('Structural MRI',(pd.to_datetime(an.FS_DATE)-pd.to_datetime(an.DTI_DATE)).dt.days,'Exact examination-date match; valid LV/ICV N=459')
ex=pd.read_csv(P/'data/expression_metadata_pre_RIN.csv');ex['DTI_minus_blood_years']=pd.to_datetime(ex.DTI_DATE).dt.year-ex.YearOfCollection
ts('Blood expression calendar-year difference',ex.DTI_minus_blood_years,'DTI calendar year minus blood collection year; not exact days')
save(pd.DataFrame(tim),'Modality_timing_available.csv')
save(pd.DataFrame({'Maximum_absolute_year_difference':[1,2], 'Participants':[sum(ex.DTI_minus_blood_years.abs()<=k) for k in [1,2]],'Total':[len(ex)]*2}),'Expression_timing_subsets.csv')
# Same-case FA specification check and VIF after adding ventricular anatomy.
ap=an[an.Subject_ID.isin(pet.Subject_ID)].copy();ap=ap.merge(pet[['Subject_ID']+roi+['DTI_FDG_gap_days']],on='Subject_ID',suffixes=('_drop',''))
controls='Age+C(Sex)+C(Group)+DTI_FDG_gap_days'
rows=[]; vifs=[]
for fa in ['FA_mean','FA_min']:
    for c in roi:
        z=ap.dropna(subset=[c,'ALPS_Mean','Age','Sex','Group','DTI_FDG_gap_days',fa,'logLV','logICV']).copy()
        if len(z)<30 or z[c].std()<1e-10:continue
        z['outcome']=z[c]
        for model,add in [('Base',''),('LV_ICV','+logLV+logICV')]:
            f=smf.ols('outcome~ALPS_Mean+'+controls+'+'+fa+add,z).fit(cov_type='HC3');ci=f.conf_int().loc['ALPS_Mean'];sdscale=z.ALPS_Mean.std()/z.outcome.std()
            rows.append(dict(ROI=c,FA=fa,Model=model,N=len(z),Beta=f.params['ALPS_Mean'],CI_low=ci.iloc[0],CI_high=ci.iloc[1],Standardized_beta=f.params['ALPS_Mean']*sdscale,Standardized_CI_low=ci.iloc[0]*sdscale,Standardized_CI_high=ci.iloc[1]*sdscale,P=f.pvalues['ALPS_Mean'],Max_Cooks_distance=f.get_influence().cooks_distance[0].max()))
            if c==roi[0]:
                for i,term in enumerate(f.model.exog_names):
                    if term!='Intercept':vifs.append(dict(FA=fa,Model=model,Term=term,N=len(z),VIF=variance_inflation_factor(f.model.exog,i)))
fa_result=pd.DataFrame(rows);d=fa_result;d['q_BH140']=d.groupby(['FA','Model']).P.transform(lambda p:multipletests(np.r_[p,np.ones(140-len(p))],method='fdr_bh')[1][:len(p)])
save(d,'FA_swap_structural_FDG_HC3.csv');save(pd.DataFrame(vifs),'Structural_predictor_VIF.csv')
# Shared tensor components: numerator and denominator are not independent biomarkers.
an['ALPS_numerator']=an[['M_L_Dxx_Assoc','M_L_Dxx_Proj','M_R_Dxx_Assoc','M_R_Dxx_Proj']].mean(axis=1,skipna=False)
an['ALPS_denominator']=an[['M_L_Dyy_Proj','M_L_Dzz_Assoc','M_R_Dyy_Proj','M_R_Dzz_Assoc']].mean(axis=1,skipna=False)
rows=[]
for component in ['ALPS_numerator','ALPS_denominator']:
    z=an.dropna(subset=[component,'logLV','logICV','FA_mean','Age','Sex','Group']).copy()
    X=pd.get_dummies(z[['Age','Sex','Group','FA_mean','logICV']],drop_first=True,dtype=float)
    for k in ['Age','FA_mean','logICV']:X[k]=stats.rankdata(X[k])
    X=sm.add_constant(X).to_numpy(float)
    def residual(x):return x-X@np.linalg.lstsq(X,x,rcond=None)[0]
    a=residual(stats.rankdata(z[component]));b=residual(stats.rankdata(z.logLV));r=np.corrcoef(a,b)[0,1];df=len(z)-np.linalg.matrix_rank(X)-1;p=2*stats.t.sf(abs(r*np.sqrt(df/(1-r*r))),df)
    rows.append(dict(Component=component,N=len(z),Target='logLV',rho=r,P=p))
d=pd.DataFrame(rows);d['q']=multipletests(d.P,method='fdr_bh')[1];save(d,'Tensor_components_ventricular_relationship.csv')
# Follow-up coverage selection, not evidence of known dropout causes.
vs=pd.read_csv(P/'results/S9_subject_visit_followup.csv');z=m.merge(vs[['Subject','N_visits']],on='Subject');rows=[]
for variable in ['Age','ALPS_Mean','MMSE_Score','CDR_SumBox']:
    a=z.loc[z.N_visits==1,variable].dropna();b=z.loc[z.N_visits>=2,variable].dropna();t=stats.ttest_ind(a,b,equal_var=False)
    rows.append(dict(Variable=variable,Single_N=len(a),Single_mean=a.mean(),Single_SD=a.std(),Repeated_N=len(b),Repeated_mean=b.mean(),Repeated_SD=b.std(),Welch_P=t.pvalue))
save(pd.DataFrame(rows),'Followup_coverage_baseline_comparison.csv')
summary=dict(Index_equal_first_PASS=int((idx==earliest).sum()),Index_different_first_PASS=int((idx!=earliest).sum()),Expression_year_gap=ex.DTI_minus_blood_years.describe().to_dict(),Expression_abs_gt1=int((ex.DTI_minus_blood_years.abs()>1).sum()),Expression_abs_gt2=int((ex.DTI_minus_blood_years.abs()>2).sum()),Protocol_switch_proxy_subjects=int((vs.Protocol_families>1).sum()),FA_swap_significant=fa_result.groupby(['FA','Model']).q_BH140.apply(lambda z:int((z<.05).sum())).rename_axis(['FA','Model']).reset_index(name='N_significant').to_dict('records'),Target_thalamic_candidates=sum('Thal' in c for c in target),Target_valid_ROIs=nvalid.value_counts().to_dict())
(A/'additional_checks_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False));print(pd.DataFrame(vifs).query('Term == "ALPS_Mean"').to_string(index=False))


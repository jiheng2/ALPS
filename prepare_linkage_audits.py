"""Read-only source linkage audit. No historical output or participant table is overwritten."""
from pathlib import Path
import json
import numpy as np,pandas as pd
q=Path(__file__).resolve().parents[1];out=q/'logs';private=q/'results/private';private.mkdir(exist_ok=True)
fs=pd.read_csv(q/'data/UCSFFSX7_20Oct2026.csv',usecols=['RID','EXAMDATE','RUNDATE','STATUS','VENTQC','ST37SV','ST96SV','ST10CV'],low_memory=False)
fs['FS_DATE']=pd.to_datetime(fs.EXAMDATE,errors='coerce');fs['Run_date']=pd.to_datetime(fs.RUNDATE,errors='coerce')
fs['LV_Total_mm3']=pd.to_numeric(fs.ST37SV,errors='coerce')+pd.to_numeric(fs.ST96SV,errors='coerce');fs['ICV_mm3']=pd.to_numeric(fs.ST10CV,errors='coerce');fs['Completed']=fs.STATUS.astype(str).str.lower().eq('complete');fs['Source_row']=fs.index
valid=fs[(fs.LV_Total_mm3>0)&(fs.ICV_mm3>0)]
def choose(na):return valid.sort_values(['RID','FS_DATE','Completed','Run_date'],kind='stable',na_position=na).drop_duplicates(['RID','FS_DATE'],keep='last').set_index(['RID','FS_DATE'])
old=choose('last');fixed=choose('first');changed=old.Source_row.ne(fixed.Source_row);fixed.reset_index().to_csv(private/'FreeSurfer_selected_source_records_PRIVATE.csv',index=False)
m=pd.read_csv(q/'data/baseline_master.csv');src=pd.read_csv(q/'data/FDG_AAL3_SUVR_v8.csv');src['RID']=pd.to_numeric(src.Subject_ID.str.extract(r'_S_(\d+)')[0]);src['PET_DATE']=pd.to_datetime(src.Subject_ID.str.extract(r'(\d{4}-\d{2}-\d{2})$')[0]);pet=m[m.WholeBrain_SUVR_AAL3.notna()&m.Months_from_Baseline.eq(0)]
roi=[x for x in m if x.startswith('AAL3_')];rows=[]
for _,r in pet.iterrows():
 z=src[src.RID.eq(r.RID)];hits=[]
 for k,s in z.iterrows():
  if np.allclose(r[roi].to_numpy(float),s[roi].to_numpy(float),equal_nan=True,atol=1e-10,rtol=1e-10) and np.isclose(r.WholeBrain_SUVR_AAL3,s.WholeBrain_SUVR_AAL3,atol=1e-10,rtol=1e-10):hits.append(k)
 if len(hits)!=1:raise ValueError('PET ROI vector has no unique dated source record')
 s=src.loc[hits[0]];gap=(s.PET_DATE-pd.Timestamp(r.DTI_DATE)).days
 if abs(gap)>183 or s.PET_DATE!=pd.Timestamp(r.FDG_Date) or abs(gap)!=r.DTI_FDG_gap_days:raise ValueError('Date/window/value provenance mismatch')
 rows.append(dict(RID=r.RID,Subject_ID=r.Subject_ID,PET_source_record_ID=s.Subject_ID,PET_source_date=s.PET_DATE,PET_minus_DTI_days=gap,PET_absolute_gap_days=abs(gap)))
d=pd.DataFrame(rows);d.to_csv(private/'PET_source_record_and_intervals_PRIVATE.csv',index=False)
report=dict(FreeSurfer_records=len(fs),Selections_changed_by_NaT_sort=int(changed.sum()),PET_records=len(d),All_unique_ROI_vectors_and_dates=True,Within_183_days=int(d.PET_absolute_gap_days.le(183).sum()),PET_before_DTI=int(d.PET_minus_DTI_days.lt(0).sum()),Original_model_interval='Absolute days',Image_UID_available=False,Dynamic_frames_independently_verified=False)
(out/'Linkage_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(report)

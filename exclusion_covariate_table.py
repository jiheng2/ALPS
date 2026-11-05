"""Summarize covariate availability in the filtered baseline cohort."""
from pathlib import Path
import pandas as pd

root = Path(__file__).resolve().parents[1]
baseline = pd.read_csv(root/'data/baseline_master.csv', low_memory=False)
columns = [
    'WMH_mL', 'WMH_mL_365', 'ADNI1_adjusted_WMH_separate', 'UCD_parenchyma_GM_WM_WMH_mL',
    'UCD_cerebrum_brain_fraction', 'Education_years', 'SBP_preindex', 'Cardiovascular_history_explicit',
    'Smoking_history_explicit', 'Preindex_history_record_available', 'Hypertension_documented',
    'Diabetes_documented', 'Dyslipidemia_documented', 'Stroke_TIA_documented', 'Coronary_disease_documented',
]
rows = []
for group, cohort in [('All', baseline), *baseline.groupby('Group')]:
    for name in columns:
        values = pd.to_numeric(cohort[name], errors='coerce').dropna()
        rows.append(dict(Group=group, Measurement=name, Population_N=len(cohort),
                         Available_N=len(values), Missing_N=len(cohort)-len(values),
                         Median=values.median(), Min=values.min(), Max=values.max()))
pd.DataFrame(rows).to_csv(root/'results/New_covariate_availability.csv', index=False, encoding='utf-8-sig')

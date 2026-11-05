import os
"""Check cohort exclusion, fresh figure provenance, and absence of old cached inputs."""
import csv
import hashlib
import json
from pathlib import Path
import pandas as pd

root = Path(__file__).resolve().parents[1]
def table(name):
    return pd.read_csv(root / name, encoding='utf-8-sig', low_memory=False)

m = table('data/baseline_master.csv')
l = table('data/longitudinal_with_site_protocol_proxy.csv')
pet = table('data/structural_PET216.csv')
flags = table('logs/Exact_exclusion_visits_PRIVATE.csv')
index = table('logs/Excluded_baseline_index_PRIVATE.csv')
source = Path(os.environ.get('ALPS_RAW_DATA',str(root/'Numeric_Source')))
base = pd.read_csv(source/'baseline_master.csv',encoding='utf-8-sig')
assert len(m) == 467 and m.PTID.is_unique
assert len(pet) == 210 and (m.WholeBrain_SUVR_AAL3.notna() & m.Months_from_Baseline.eq(0)).sum() == 210
assert m[['ALPS_Left','ALPS_Right','ALPS_Mean']].max().max() <= 2.5
assert len(set(base.PTID) - set(m.PTID)) == 3
assert not set(l.Subject) & (set(base.PTID)-set(m.PTID))
assert flags.groupby('Modality').size().to_dict() == {'DTI': 7, 'PET': 6}
assert len(index) == 3
assert set(index.Subject) == set(base.PTID) - set(m.PTID)
assert not set(flags.loc[flags.Modality.eq('DTI'),'Exact_visit']) & set(l.Subject_ID)
assert not set(flags.loc[flags.Modality.eq('PET'),'Exact_visit']) & set(pet.PET_visit_id)
assert set(l.Subject).issubset(pd.read_csv(source/'longitudinal_with_site_protocol_proxy.csv',encoding='utf-8-sig').Subject)
figures = sorted((root/'figures').glob('Figure*.pdf'))
expected = {f'Figure{i}' for i in range(1,6)} | {f'FigureS{i}' for i in range(1,20)}
assert {p.stem for p in figures} == expected
assert all(p.read_bytes()[:5] == b'%PDF-' and p.stat().st_size > 1000 for p in figures)
provenance = table('logs/plot_table_provenance.csv')
for r in provenance.itertuples():
    path = root/'results'/r.Table
    assert path.is_file() and hashlib.md5(path.read_bytes()).hexdigest() == r.MD5, r.Table
summary = dict(Baseline=len(m),PET=len(pet),Longitudinal_input=len(l),
               Figures=len(figures),Figure_table_sources=len(provenance),
               Generated_result_csvs=len(list((root/'results').rglob('*.csv'))))
(root/'logs/Exclusion_validation.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(summary)




"""Select analysis cohorts from finalized post-correction numerical source tables."""
import os
from pathlib import Path
import json
import argparse
import runpy
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = Path(os.environ.get('ALPS_RAW_DATA',str(ROOT/'Numeric_Source')))
OUT = ROOT / 'data'


def read(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False, encoding='utf-8-sig')


def write(table, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(dest, index=False, encoding='utf-8-sig')


def regroup_survival(table):
    alps = pd.to_numeric(table['ALPS_Mean'], errors='raise')
    table['ALPS_Group'] = alps.ge(alps.median()).map({True: 'High ALPS', False: 'Low ALPS'})


def apply_exclusions(policy='primary'):
    runpy.run_path(str(ROOT/'python/verify_numeric_manifest.py'),run_name='__main__')
    OUT.mkdir(parents=True, exist_ok=True)
    flagged = read(SRC / 'QC_11_Cases_PRIVATE.csv')
    baseline = read(SRC / 'baseline_master.csv')
    dti_flags = set(flagged.loc[flagged.Modality == 'DTI', 'Visit_key'])
    pet_flags = set(flagged.loc[flagged.Modality == 'PET', 'Visit_key'])
    extreme = baseline[['ALPS_Left','ALPS_Right','ALPS_Mean']].apply(pd.to_numeric).gt(2.5).any(axis=1)
    assert len(dti_flags) == 6 and len(pet_flags) == 5 and extreme.sum() == 1
    extreme_keys = set(baseline.loc[extreme, 'Subject_ID']) if policy=='primary' else set()
    exact_dti = dti_flags | extreme_keys if policy!='none' else set()
    lost = baseline[baseline.Subject_ID.isin(exact_dti)]
    expected={'primary':(467,210,3),'technical':(468,211,2),'none':(470,216,0)}[policy]
    assert len(lost) == expected[2]
    lost_ptid, lost_rid = set(lost.PTID), set(lost.RID)
    lost_pairs = lost.loc[lost.FDG_Date.ne(''), 'PTID'] + '_' + lost.loc[lost.FDG_Date.ne(''), 'FDG_Date']
    pet_remove = pet_flags | set(lost_pairs) if policy!='none' else set()
    pet_ptid = {v.rsplit('_',1)[0] for v in pet_remove}

    contract=json.loads((ROOT/'fixed_resources/numeric_input_contract.json').read_text(encoding='utf8'))
    audit = []
    for relative_name in contract['model_inputs']:
        relative = Path(relative_name)
        source = SRC / relative
        name = source.name
        dest = OUT / relative
        dest.parent.mkdir(parents=True,exist_ok=True)
        if name in {'DXSUM_26May2026.csv','UCSFFSX7_20Oct2026.csv','AAL3_voxel_counts.csv','expression_probe_info.csv','previously_displayed_GO_terms.csv'}:
            dest.write_bytes(source.read_bytes())
            continue
        table = read(source)
        initial = len(table)
        if name == 'expression_matched_matrix.csv':
            metadata = read(SRC/'expression_metadata_pre_RIN.csv')
            drop = metadata.loc[metadata.RID.isin(lost_rid),'SampleID'].tolist()
            table = table.drop(columns=drop)
        else:
            keep = pd.Series(True, index=table.index)
            for col in ('Subject','PTID'):
                if col in table:
                    keep &= ~table[col].isin(lost_ptid)
            if 'RID' in table:
                keep &= ~table.RID.isin(lost_rid)
            if 'Subject_ID' in table:
                sid = table.Subject_ID
                keep &= ~sid.str.rsplit('_',n=1).str[0].isin(lost_ptid) & ~sid.isin(lost_ptid | exact_dti)
            if 'Subject' in table and 'Acq_Date' in table:
                keep &= ~(table.Subject + '_' + table.Acq_Date).isin(exact_dti)
            for col in ('FDG_VISIT','PET_visit_id'):
                if col in table:
                    keep &= ~table[col].isin(pet_remove)
            if name == 'FDG_AAL3_SUVR_v8.csv':
                keep &= ~table.Subject_ID.isin(pet_remove)
            if name == 'structural_PET216.csv':
                keep &= ~table.Subject.isin(pet_ptid)
            table = table.loc[keep].copy()
            if name in {'survival_CN.csv', 'survival_MCI.csv'}:
                regroup_survival(table)
            if name in {'baseline_master.csv','baseline_with_dated_amyloid.csv','BBB_analysis.csv','PRS_analysis.csv','anatomy_linked_PRIVATE.csv'} and 'FDG_Date' in table:
                flagged_pet = (table.Subject + '_' + table.FDG_Date).isin(pet_remove)
                fields = [c for c in table if c.startswith('AAL3_') or c.startswith('WholeBrain_SUVR') or c in
                          {'FDG_Date','DTI_FDG_gap_days','FDG_DATE','PET_visit_id','DTI_PET_gap_days'}]
                table.loc[flagged_pet,fields] = ''
            if name in {'ALPS_all_visits_corrected.csv','longitudinal_with_site_protocol_proxy.csv'} and 'N_Scans_per_Subject' in table:
                table['N_Scans_per_Subject'] = table.groupby('Subject').Subject.transform('size')
        write(table, dest)
        audit.append(dict(File=str(relative),Action='cohort policy applied',Before=initial,After=len(table)))
    write(pd.DataFrame(audit), ROOT/'logs'/'Exclusion_input_audit.csv')
    write(lost[['Subject','RID','Subject_ID','Group','FDG_Date']],ROOT/'logs'/'Excluded_baseline_index_PRIVATE.csv')
    write(pd.DataFrame([
      *[{'Modality':'DTI','Exact_visit':k,'Reason':'Figure 1 >2.5 baseline' if k in extreme_keys else 'QC physical-value anomaly'} for k in sorted(exact_dti)],
      *[{'Modality':'PET','Exact_visit':k,'Reason':'QC reference coverage' if k in pet_flags else 'DTI baseline pair removed'} for k in sorted(pet_remove)]
    ],columns=['Modality','Exact_visit','Reason']),ROOT/'logs'/'Exact_exclusion_visits_PRIVATE.csv')
    assert len(read(OUT/'baseline_master.csv')) == expected[0]
    assert len(read(OUT/'structural_PET216.csv')) == expected[1]
    for filename, key in [('baseline_master.csv','Subject_ID'),('longitudinal_with_site_protocol_proxy.csv','Subject_ID')]:
        r = read(OUT/filename)
        assert not r[key].isin(exact_dti).any() and not r.Subject.isin(lost_ptid).any()
    result = {'policy':policy,'baseline':expected[0],'PET':expected[1],'removed_index_participants':len(lost),
              'longitudinal_source_rows':len(read(OUT/'longitudinal_with_site_protocol_proxy.csv')),
              'source_DTI_technical_flags':len(dti_flags),'source_PET_technical_flags':len(pet_flags),
              'applied_DTI_scan_keys':len(exact_dti),'applied_PET_scan_keys':len(pet_remove)}
    (ROOT/'logs'/'Exclusion_inputs_summary.json').write_text(json.dumps(result, indent=2),encoding='utf-8')
    print(result)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--policy',choices=['primary','technical','none'],default='primary')
    apply_exclusions(parser.parse_args().policy)


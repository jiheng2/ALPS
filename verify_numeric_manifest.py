"""Require the authoritative frozen corrected numerical source, without image-correction programs."""
from pathlib import Path
import os,csv,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
source=Path(os.environ.get('ALPS_RAW_DATA',str(ROOT/'Numeric_Source')))
contract=json.loads((ROOT/'fixed_resources/numeric_input_contract.json').read_text(encoding='utf8'))
expected=set(contract['model_inputs']+contract['exclusion_rules']+contract['flow_only_inputs'])
with (ROOT/'fixed_resources/EXPECTED_NUMERIC_SOURCE_SHA256.csv').open(encoding='utf-8-sig',newline='')as f:rows=list(csv.DictReader(f))
if {r['File']for r in rows}!=expected:raise ValueError('Input contract and locked source manifest differ')
checked=[]
for r in rows:
 p=source/r['File']
 if not p.is_file():raise FileNotFoundError('Required corrected numerical source is missing: '+str(p))
 h=hashlib.sha256()
 with p.open('rb')as f:
  for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
 if p.stat().st_size!=int(r['Bytes'])or h.hexdigest()!=r['SHA256']:raise ValueError('Authoritative corrected numerical source mismatch: '+r['File'])
 checked.append({'File':r['File'],'SHA256_matches':True,'Role':r['Role']})
(ROOT/'logs').mkdir(exist_ok=True)
with (ROOT/'logs/Numeric_source_verification.csv').open('w',encoding='utf8',newline='')as f:
 w=csv.DictWriter(f,fieldnames=checked[0]);w.writeheader();w.writerows(checked)
print('CORRECTED_NUMERIC_SOURCE_VERIFIED',len(checked),'files; integrity check only, not upstream image-correction validation')

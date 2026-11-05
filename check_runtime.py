import importlib,json,sys
from pathlib import Path
parent=Path(__file__).resolve().parents[1]
root=parent if(parent/'project_paths.json').exists()else parent/'Local_Reproduction'
names=['numpy','pandas','scipy','statsmodels','sklearn','networkx','nibabel']
rows=[]
for name in names:
    try:
        module=importlib.import_module(name)
        rows.append({'Module':name,'Installed':True,'Version':str(getattr(module,'__version__','unknown'))})
    except Exception as exc:
        rows.append({'Module':name,'Installed':False,'Error':str(exc)})
(root/'logs/Python_runtime_check.json').write_text(json.dumps({'Python':sys.version,'Modules':rows},indent=2),encoding='utf-8')
for row in rows:print(row)
if not all(row['Installed']for row in rows):raise SystemExit(1)
print('PYTHON_RUNTIME_OK')

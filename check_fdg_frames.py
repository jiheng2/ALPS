"""Select a verified complete uptake window; frame indices are zero-based.
TSV: frame_index, start_seconds_post_injection, duration_seconds, registered.
Output is start-index and frame-count for fslroi; unequal-duration frames are rejected.
"""
import csv,sys,math
def select(path,start=1800.,end=3600.):
 with open(path,encoding='utf-8-sig',newline='')as f:rows=list(csv.DictReader(f,delimiter='\t'))
 if not rows:raise ValueError('Empty timing table')
 data=[]
 for r in rows:
  i=int(r['frame_index']);t=float(r['start_seconds_post_injection']);d=float(r['duration_seconds'])
  if not all(math.isfinite(x)for x in [t,d])or d<=0:raise ValueError('Invalid frame timing')
  if r['registered'].strip().lower()not in ['true','1','yes']:raise ValueError('Frame registration is not documented')
  data.append((i,t,d))
 data.sort();ids=[r[0]for r in data]
 if ids!=list(range(len(data))):raise ValueError('Frame indices must be unique and contiguous from zero')
 chosen=[r for r in data if r[1]>=start and r[1]+r[2]<=end]
 if not chosen or abs(chosen[0][1]-start)>1e-6 or abs(chosen[-1][1]+chosen[-1][2]-end)>1e-6:raise ValueError('Frames do not exactly cover the requested window')
 if any(abs(a[1]+a[2]-b[1])>1e-6 for a,b in zip(chosen,chosen[1:])):raise ValueError('Frame overlap or gap')
 if any(abs(r[2]-chosen[0][2])>1e-6 for r in chosen):raise ValueError('Unequal durations require explicit duration-weighted averaging')
 return chosen[0][0],len(chosen),len(data)
if __name__=='__main__':
 a=select(sys.argv[1],float(sys.argv[2]),float(sys.argv[3]));print(*a)

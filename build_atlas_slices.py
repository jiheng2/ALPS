"""Build standard atlas/background slices from separately licensed NIfTI resources."""
from pathlib import Path
import argparse,json,numpy as np,nibabel as nib
from nibabel.processing import resample_from_to
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--atlas',required=True);ap.add_argument('--background',required=True);ap.add_argument('--out',default=str(Path(__file__).resolve().parents[1]/'assets'));a=ap.parse_args()
 out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
 atlas=nib.as_closest_canonical(nib.load(a.atlas));bg=resample_from_to(nib.load(a.background),atlas,order=1)
 x=atlas.get_fdata();y=bg.get_fdata();planes=[0,-18,10]
 for axis,mm in enumerate(planes):
  point=np.zeros(3);point[axis]=mm;k=int(round(nib.affines.apply_affine(np.linalg.inv(atlas.affine),point)[axis]))
  if not 0<=k<x.shape[axis]:raise ValueError('Requested plane is outside the atlas')
  np.savetxt(out/f'atlas_slice{axis+1}.csv',np.take(x,k,axis=axis),delimiter=',',fmt='%.0f')
  np.savetxt(out/f'background_slice{axis+1}.csv',np.take(y,k,axis=axis),delimiter=',',fmt='%.6f')
 nib.save(atlas,out/'AAL3v1.nii.gz')
 (out/'slice_provenance.json').write_text(json.dumps({'Planes_mm':planes,'Orientation':'Closest canonical RAS; sagittal, coronal, axial','Atlas_affine':atlas.affine.tolist(),'Background':'Provided standard template linearly resampled to atlas'},indent=2))
if __name__=='__main__':main()

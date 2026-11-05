#!/usr/bin/env bash
# Prospective preparation only; does not validate historical ADNI images.
set -euo pipefail
input=${1:?input NIfTI required}
output_dir=${2:?new empty output directory required}
timing=${3:-}
: "${ALPS_PYTHON:=python3}"
: "${WINDOW_START_SECONDS:=1800}"
: "${WINDOW_END_SECONDS:=3600}"
if [[ -e "$output_dir" ]]; then echo 'Refusing to reuse existing outputs' >&2; exit 1; fi
[[ -f "$input" ]] || exit 1
for cmd in fslval fslroi fslmaths; do command -v "$cmd" >/dev/null; done
n=$(fslval "$input" dim4 | tr -d '[:space:]')
[[ "$n" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid frame count' >&2; exit 1; }
mkdir -p "$output_dir"
trap 'echo "FDG static preparation failed" >&2' ERR
if [[ "$n" == 1 ]]; then
  : "${STATIC_IMAGE_PROVENANCE:?Document the ADNI 3D image product and uptake window}"
  cp "$input" "$output_dir/fdg_static.nii.gz"
  printf '3D source: %s\n' "$STATIC_IMAGE_PROVENANCE" > "$output_dir/frame_selection.log"
else
  [[ -f "$timing" ]] || { echo '4D input requires frame timing and registration documentation' >&2; exit 1; }
  script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
  selection=$("$ALPS_PYTHON" "$script_dir/check_fdg_frames.py" "$timing" "$WINDOW_START_SECONDS" "$WINDOW_END_SECONDS")
  read -r first count documented_n <<< "$selection"
  [[ "$documented_n" == "$n" ]] || { echo 'Timing rows do not match NIfTI frames' >&2; exit 1; }
  fslroi "$input" "$output_dir/selected_frames.nii.gz" "$first" "$count"
  fslmaths "$output_dir/selected_frames.nii.gz" -Tmean "$output_dir/fdg_static.nii.gz"
  cp "$timing" "$output_dir/frame_timing.tsv"
  printf 'Window: %s–%s seconds; first frame %s; count %s\n' "$WINDOW_START_SECONDS" "$WINDOW_END_SECONDS" "$first" "$count" > "$output_dir/frame_selection.log"
fi
[[ -s "$output_dir/fdg_static.nii.gz" ]] || exit 1
printf 'Static preparation complete\n' >> "$output_dir/frame_selection.log"

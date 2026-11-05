#!/usr/bin/env bash
# Prospective PET processing, requiring a separately checked static image and frame_selection.log.
# Not run against this study's raw dynamic images; raw-image review is deferred.
set -euo pipefail

FDG_DIR="${1:?用法: bash preprocess_fdg_v5.sh <FDG_subject_dir>}"
T1_BASE_DIR="${T1_BASE_DIR:?Set the T1 input directory}"
REF_MASK="${REF_MASK:?Set the validated pons-vermis mask}"

MNI_T1_2mm="$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz"
MNI_T1_2mm_brain="$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz"
BRAIN_MASK="$FSLDIR/data/standard/MNI152_T1_2mm_brain_mask.nii.gz"

[ ! -f "$REF_MASK" ] && { echo "❌ 缺 $REF_MASK"; exit 1; }

SUBJ_DATE=$(basename "$FDG_DIR")
SUBJ_ID=$(echo "$SUBJ_DATE" | sed 's/_[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}$//')
FDG_DATE=$(echo "$SUBJ_DATE" | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}$')

LOG="$FDG_DIR/preprocess.log"
echo "==========================================" > "$LOG"
echo "FDG v5: $FDG_DIR" >> "$LOG"
echo "Subject: $SUBJ_ID, FDG date: $FDG_DATE" >> "$LOG"
echo "Started: $(date)" >> "$LOG"

# ── 1. 配对 T1 ──
BEST_T1=""
BEST_DIFF=99999
for T1_CAND in "$T1_BASE_DIR/${SUBJ_ID}"_*; do
    [ -d "$T1_CAND" ] || continue
    T1_DATE=$(basename "$T1_CAND" | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}$')
    [ -z "$T1_DATE" ] && continue
    [ ! -f "$T1_CAND/T1.nii.gz" ] && continue
    DIFF=$(( ($(date -d "$FDG_DATE" +%s) - $(date -d "$T1_DATE" +%s)) / 86400 ))
    DIFF=${DIFF#-}
    [ "$DIFF" -lt "$BEST_DIFF" ] && { BEST_DIFF=$DIFF; BEST_T1="$T1_CAND/T1.nii.gz"; }
done

[ -z "$BEST_T1" ] && { echo "❌ 找不到 T1" | tee -a "$LOG"; exit 1; }
echo "✓ T1: $BEST_T1 (差 ${BEST_DIFF}天)" | tee -a "$LOG"

cd "$FDG_DIR"

# Static FDG must be prepared separately in a new output directory with timing provenance.
[[ -s fdg_static.nii.gz && -s frame_selection.log ]] || { echo 'Run prepare_fdg_static_checked.sh first' >&2; exit 1; }
# ── 3. BET ──
if [ ! -f T1_brain.nii.gz ]; then
    echo "[2/7] T1 BET..." | tee -a "$LOG"
    bet "$BEST_T1" T1_brain.nii.gz -f 0.35 -B
fi

# ── 4. FLIRT T1→MNI ──
if [ ! -f T1_2MNI_linear.mat ]; then
    echo "[3/7] T1→MNI FLIRT (12 DOF)..." | tee -a "$LOG"
    flirt -in T1_brain.nii.gz -ref "$MNI_T1_2mm_brain" \
          -out T1_2MNI_linear.nii.gz -omat T1_2MNI_linear.mat \
          -dof 12 -cost corratio
fi

# ── 5. FNIRT T1→MNI ──
if [ ! -f T1_2MNI_warp.nii.gz ]; then
    echo "[4/7] T1→MNI FNIRT (~5 min)..." | tee -a "$LOG"
    fnirt --in="$BEST_T1" --ref="$MNI_T1_2mm" \
          --aff=T1_2MNI_linear.mat \
          --cout=T1_2MNI_warp.nii.gz \
          --iout=T1_2MNI_nonlin.nii.gz \
          --config=T1_2_MNI152_2mm
fi

# ── 6. FDG→T1 ──
if [ ! -f fdg_2T1.nii.gz ]; then
    echo "[5/7] FDG→T1 FLIRT (6 DOF)..." | tee -a "$LOG"
    flirt -in fdg_static.nii.gz -ref T1_brain.nii.gz \
          -out fdg_2T1.nii.gz -omat fdg_2T1.mat \
          -dof 6 -cost mutualinfo
fi

# ── 7. FDG→MNI ──
if [ ! -f fdg_MNI.nii.gz ]; then
    echo "[6/7] FDG→MNI applywarp..." | tee -a "$LOG"
    applywarp --in=fdg_static.nii.gz --ref="$MNI_T1_2mm" \
              --warp=T1_2MNI_warp.nii.gz \
              --premat=fdg_2T1.mat \
              --out=fdg_MNI.nii.gz \
              --interp=trilinear
fi

# ── 8. SUVR ──
if [ ! -f fdg_SUVR.nii.gz ]; then
    echo "[7/7] SUVR (÷ Pons+Vermis)..." | tee -a "$LOG"
    REF=$(fslstats fdg_MNI.nii.gz -k "$REF_MASK" -M | tr -d '[:space:]')
    echo "  Ref 均值 = $REF" | tee -a "$LOG"
    
    if awk -v r="$REF" 'BEGIN { exit !(r ~ /^[+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ && r+0>0) }'; then
        fslmaths fdg_MNI.nii.gz -div "$REF" fdg_SUVR.nii.gz
        SUVR=$(fslstats fdg_SUVR.nii.gz -k "$BRAIN_MASK" -M | awk '{printf "%.3f", $1}')
        echo "  ✓ 全脑 SUVR = $SUVR" | tee -a "$LOG"
        
        bad=$(awk "BEGIN {print ($SUVR<0.7||$SUVR>1.4)?1:0}")
        [ "$bad" -eq 1 ] && echo "  ⚠ SUVR 异常 ($SUVR)" | tee -a "$LOG"
    else
        echo "Invalid reference mean" | tee -a "$LOG"; exit 1
    fi
fi

[[ -s fdg_SUVR.nii.gz ]] || { echo "Missing SUVR output" >&2; exit 1; }
echo "✓ Done: $(date)" | tee -a "$LOG"

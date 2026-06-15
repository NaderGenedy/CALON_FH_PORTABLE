#!/bin/bash
# ==============================================================================
# CALON-2 PIPELINE: STEP 07c — EXTRACT HIGH-QUALITY CARDIAC MRI IDPs
# ==============================================================================
# PURPOSE: Extract Category 157 (Bai et al. CNN pipeline) cardiac MRI
#          imaging-derived phenotypes — superior quality to Category 133
#          (scanner-automated) already extracted in p22420-p22426.
#
#          Also extracts imaging visit date (p53_i2) to verify that MRI
#          measurements occurred BEFORE the index ASCVD event.
#
# KEY FIELDS:
#   p24105  LV myocardial mass (strong ASCVD predictor — LVH)
#   p24103  LV ejection fraction (systolic function)
#   p24100  LV end-diastolic volume (cardiac remodeling)
#   p24120  Ascending aorta distensibility (arterial stiffness)
#   p24110  LA maximum volume (atrial remodeling)
#   p24181  Global longitudinal strain (subclinical dysfunction)
#   p53_i2  Imaging visit date (for temporal filtering)
#
# NOTE: hsCRP is ALREADY in base dataset as 'crp' (field 30710,
#       1,619/1,623 non-missing). No new extraction needed.
#
# NOTE: Category 133 (p22420-p22426) has QUALITY ISSUES — Siemens InlineVF
#       scanner-automated, no expert QC. Category 157 uses validated CNN
#       pipeline with ~81,000 participants.
#
# RUN ON: UKB-RAP Platform (JupyterLab Bash terminal)
# AUTHOR: Dr Nader Genedy
# DATE:   February 2026
# ==============================================================================

RECORD="project-J6K175jJZ01XppV5477pkYvJ:record-J6K32f8JgZ4JX4gYF39zjBQz"

set -e

extract_with_retry() {
  local output="$1"
  local fields="$2"
  local max_retries=3

  for attempt in $(seq 1 $max_retries); do
    echo "  Attempt ${attempt}/${max_retries}..."
    rm -f "$output"
    if dx extract_dataset "$RECORD" \
         --fields "$fields" \
         --output "$output" \
         --delimiter "," 2>&1; then
      local lines=$(wc -l < "$output")
      if [ "$lines" -gt 1 ]; then
        echo "  ✓ OK: $output ($lines lines)"
        return 0
      else
        echo "  ✗ Empty result"
      fi
    else
      echo "  ✗ Failed attempt $attempt"
    fi
    sleep 3
  done
  echo "  ✗✗ FAILED after $max_retries attempts: $output"
  return 1
}

echo ""
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║  CALON-2: EXTRACT HIGH-QUALITY CARDIAC MRI IDPs                    ║"
echo "║  Category 157 (Bai et al. CNN Pipeline, ~81K participants)         ║"
echo "║  + Imaging visit date for temporal filtering                       ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

# ==============================================================================
# BATCH 1: IMAGING VISIT DATE (for temporal filtering)
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 1: Imaging Visit Date (p53_i2 = first imaging visit)"
echo "  This is essential: MRI must precede ASCVD event!"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

# Try both with and without instance suffix
IMG_DATE_SUCCESS=0
for fmt in "participant.p53_i2" "participant.p53_i2_a0"; do
  if extract_with_retry calon_extra_imaging_date.csv "participant.eid,$fmt"; then
    non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_imaging_date.csv)
    echo "  Non-empty imaging dates: $non_empty"
    if [ "$non_empty" -gt 0 ]; then
      IMG_DATE_SUCCESS=1
      break
    fi
  fi
done

# Also get repeat imaging date (instance 3) if available
echo ""
echo "  Also trying repeat imaging date (p53_i3)..."
for fmt in "participant.p53_i3" "participant.p53_i3_a0"; do
  if extract_with_retry calon_extra_imaging_date_repeat.csv "participant.eid,$fmt" 2>/dev/null; then
    non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_imaging_date_repeat.csv)
    if [ "$non_empty" -gt 0 ]; then
      echo "  ✓ Repeat imaging dates: $non_empty non-empty"
      break
    fi
  fi
done
echo ""

# ==============================================================================
# BATCH 2: CORE CARDIAC MRI IDPs — Category 157 (Bai et al.)
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 2: Category 157 — LV Function & Mass (most important)"
echo "  p24100  LV end-diastolic volume"
echo "  p24101  LV end-systolic volume"
echo "  p24102  LV stroke volume"
echo "  p24103  LV ejection fraction"
echo "  p24104  LV cardiac output"
echo "  p24105  LV myocardial mass ★ (LVH = strong ASCVD predictor)"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

# Try without instance first (Category 157 may use different instancing)
LV_SUCCESS=0
for suffix in "" "_i2" "_i0"; do
  LV_FIELDS="participant.eid"
  LV_FIELDS="$LV_FIELDS,participant.p24100${suffix}"
  LV_FIELDS="$LV_FIELDS,participant.p24101${suffix}"
  LV_FIELDS="$LV_FIELDS,participant.p24102${suffix}"
  LV_FIELDS="$LV_FIELDS,participant.p24103${suffix}"
  LV_FIELDS="$LV_FIELDS,participant.p24104${suffix}"
  LV_FIELDS="$LV_FIELDS,participant.p24105${suffix}"

  echo "  Trying suffix '$suffix'..."
  if extract_with_retry calon_extra_mri_lv.csv "$LV_FIELDS"; then
    non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_mri_lv.csv)
    echo "  Non-empty LV values: $non_empty"
    if [ "$non_empty" -gt 0 ]; then
      LV_SUCCESS=1
      break
    fi
  fi
done

if [ $LV_SUCCESS -eq 0 ]; then
  echo "  ✗ Category 157 LV fields not dispensed."
  echo "  Will use Category 133 (p22420-p22426) as fallback."
fi
echo ""

# ==============================================================================
# BATCH 3: AORTIC DISTENSIBILITY — Category 157
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 3: Category 157 — Aortic Distensibility (arterial stiffness)"
echo "  p24118  Ascending aorta max area"
echo "  p24119  Ascending aorta min area"
echo "  p24120  Ascending aorta distensibility ★"
echo "  p24121  Descending aorta max area"
echo "  p24122  Descending aorta min area"
echo "  p24123  Descending aorta distensibility"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

AORTA_SUCCESS=0
for suffix in "" "_i2" "_i0"; do
  AO_FIELDS="participant.eid"
  AO_FIELDS="$AO_FIELDS,participant.p24118${suffix}"
  AO_FIELDS="$AO_FIELDS,participant.p24119${suffix}"
  AO_FIELDS="$AO_FIELDS,participant.p24120${suffix}"
  AO_FIELDS="$AO_FIELDS,participant.p24121${suffix}"
  AO_FIELDS="$AO_FIELDS,participant.p24122${suffix}"
  AO_FIELDS="$AO_FIELDS,participant.p24123${suffix}"

  echo "  Trying suffix '$suffix'..."
  if extract_with_retry calon_extra_mri_aorta.csv "$AO_FIELDS"; then
    non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_mri_aorta.csv)
    echo "  Non-empty aortic values: $non_empty"
    if [ "$non_empty" -gt 0 ]; then
      AORTA_SUCCESS=1
      break
    fi
  fi
done
echo ""

# ==============================================================================
# BATCH 4: LEFT ATRIAL VOLUMES + RV — Category 157
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 4: Category 157 — LA Volumes + RV Function"
echo "  p24110  LA maximum volume ★ (stroke/AF risk)"
echo "  p24111  LA minimum volume"
echo "  p24112  LA stroke volume"
echo "  p24113  LA ejection fraction"
echo "  p24106  RV end-diastolic volume"
echo "  p24109  RV ejection fraction"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

LA_SUCCESS=0
for suffix in "" "_i2" "_i0"; do
  LA_FIELDS="participant.eid"
  LA_FIELDS="$LA_FIELDS,participant.p24110${suffix}"
  LA_FIELDS="$LA_FIELDS,participant.p24111${suffix}"
  LA_FIELDS="$LA_FIELDS,participant.p24112${suffix}"
  LA_FIELDS="$LA_FIELDS,participant.p24113${suffix}"
  LA_FIELDS="$LA_FIELDS,participant.p24106${suffix}"
  LA_FIELDS="$LA_FIELDS,participant.p24109${suffix}"

  echo "  Trying suffix '$suffix'..."
  if extract_with_retry calon_extra_mri_la_rv.csv "$LA_FIELDS"; then
    non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_mri_la_rv.csv)
    echo "  Non-empty LA/RV values: $non_empty"
    if [ "$non_empty" -gt 0 ]; then
      LA_SUCCESS=1
      break
    fi
  fi
done
echo ""

# ==============================================================================
# BATCH 5: MYOCARDIAL STRAIN — Category 157
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 5: Category 157 — Myocardial Strain (subclinical dysfunction)"
echo "  p24140  LV mean wall thickness"
echo "  p24157  Circumferential strain (global)"
echo "  p24174  Radial strain (global)"
echo "  p24181  Longitudinal strain (global) ★"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

STRAIN_SUCCESS=0
for suffix in "" "_i2" "_i0"; do
  ST_FIELDS="participant.eid"
  ST_FIELDS="$ST_FIELDS,participant.p24140${suffix}"
  ST_FIELDS="$ST_FIELDS,participant.p24157${suffix}"
  ST_FIELDS="$ST_FIELDS,participant.p24174${suffix}"
  ST_FIELDS="$ST_FIELDS,participant.p24181${suffix}"

  echo "  Trying suffix '$suffix'..."
  if extract_with_retry calon_extra_mri_strain.csv "$ST_FIELDS"; then
    non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_mri_strain.csv)
    echo "  Non-empty strain values: $non_empty"
    if [ "$non_empty" -gt 0 ]; then
      STRAIN_SUCCESS=1
      break
    fi
  fi
done
echo ""

# ==============================================================================
# BATCH 6: CATEGORY 162 (ALTERNATIVE) — Petersen et al.
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 6: Category 162 (Petersen et al.) — Alternative if 157 fails"
echo "  p31063  LV mass (alternative)"
echo "  p31060  LV ejection fraction (alternative)"
echo "  p31075  LA biplanar max volume"
echo "  p31085  Pericardial fat area"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

CAT162_SUCCESS=0
for suffix in "" "_i2" "_i0"; do
  C162_FIELDS="participant.eid"
  C162_FIELDS="$C162_FIELDS,participant.p31063${suffix}"
  C162_FIELDS="$C162_FIELDS,participant.p31060${suffix}"
  C162_FIELDS="$C162_FIELDS,participant.p31075${suffix}"
  C162_FIELDS="$C162_FIELDS,participant.p31085${suffix}"

  echo "  Trying suffix '$suffix'..."
  if extract_with_retry calon_extra_mri_cat162.csv "$C162_FIELDS"; then
    non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_mri_cat162.csv)
    echo "  Non-empty Cat 162 values: $non_empty"
    if [ "$non_empty" -gt 0 ]; then
      CAT162_SUCCESS=1
      break
    fi
  fi
done
echo ""

# ==============================================================================
# UPLOAD ALL FILES
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "UPLOADING ALL FILES"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

for f in calon_extra_imaging_date.csv \
         calon_extra_imaging_date_repeat.csv \
         calon_extra_mri_lv.csv \
         calon_extra_mri_aorta.csv \
         calon_extra_mri_la_rv.csv \
         calon_extra_mri_strain.csv \
         calon_extra_mri_cat162.csv; do
  if [ -f "$f" ] && [ -s "$f" ]; then
    lines=$(wc -l < "$f")
    if [ "$lines" -gt 1 ]; then
      dx upload "$f" --destination "/" 2>/dev/null || true
      echo "  ✓ Uploaded: $f ($lines lines)"
    fi
  else
    echo "  ✗ Skipped: $f (missing or empty)"
  fi
done

echo ""
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║  MRI IDP EXTRACTION COMPLETE                                       ║"
echo "╠══════════════════════════════════════════════════════════════════════╣"

if [ $LV_SUCCESS -eq 1 ]; then
  echo "║  ✓ LV Function & Mass (Cat 157): SUCCESS                         ║"
else
  echo "║  ✗ LV Function (Cat 157): Failed — will use Cat 133 fallback     ║"
fi

if [ $AORTA_SUCCESS -eq 1 ]; then
  echo "║  ✓ Aortic Distensibility (Cat 157): SUCCESS                      ║"
else
  echo "║  ✗ Aortic Distensibility: Not available                          ║"
fi

if [ $LA_SUCCESS -eq 1 ]; then
  echo "║  ✓ LA Volumes + RV (Cat 157): SUCCESS                            ║"
else
  echo "║  ✗ LA Volumes: Not available                                     ║"
fi

if [ $STRAIN_SUCCESS -eq 1 ]; then
  echo "║  ✓ Myocardial Strain (Cat 157): SUCCESS                          ║"
else
  echo "║  ✗ Myocardial Strain: Not available                              ║"
fi

if [ $CAT162_SUCCESS -eq 1 ]; then
  echo "║  ✓ Category 162 (Petersen): SUCCESS                              ║"
fi

if [ $IMG_DATE_SUCCESS -eq 1 ]; then
  echo "║  ✓ Imaging Visit Date: Extracted for temporal filtering          ║"
else
  echo "║  ✗ Imaging Visit Date: Not extracted (will use assessment date)  ║"
fi

echo "║                                                                    ║"
echo "║  REMINDER: hsCRP already in base dataset as 'crp' field           ║"
echo "║            (p30710, 1619/1623 non-missing)                         ║"
echo "║                                                                    ║"
echo "╠══════════════════════════════════════════════════════════════════════╣"
echo "║  NEXT STEPS:                                                       ║"
echo "║  1. Download all calon_extra_mri_*.csv + imaging_date.csv          ║"
echo "║  2. Place in calon_ukb_pipeline/                                   ║"
echo "║  3. Run 05_CALON2_develop.R (updated for MRI + hsCRP + temporal)   ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"

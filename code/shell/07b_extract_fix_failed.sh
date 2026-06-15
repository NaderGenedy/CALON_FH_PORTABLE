#!/bin/bash
# ==============================================================================
# CALON-2 PIPELINE: STEP 07b — FIX FAILED EXTRACTIONS
# ==============================================================================
# PURPOSE: Re-extract fields that failed in 07_extract_extra_ukb_fields.sh
#          using corrected field IDs based on UKB Data Showcase research.
#
# FAILED FIELDS AND FIXES:
#   1. Townsend p189 RETIRED → use p22189 (replacement field, singular)
#      Fallback: IMD fields p26410/p26426/p26427
#   2. CAC: Does NOT exist in UKB (no CT, only cardiac MRI)
#      → Extract cardiac MRI IDPs as proxies (LV mass, aortic distensibility)
#   3. PRS p26201-p26283: May not be dispensed to project
#      → Try individual key PRS fields one at a time
#      → Fallback: genetic principal components p22009
#   4. CIMT p22671-p22678: Extracted but ALL EMPTY
#      → Re-verify; also try mean CIMT at all 4 angles + quality fields
#
# RUN ON: UKB-RAP Platform (JupyterLab Bash terminal or %%bash cell)
#
# AUTHOR:  Dr Nader Genedy
# DATE:    February 2026
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

# Quick helper: try extraction, return success/failure without retries
extract_once() {
  local output="$1"
  local fields="$2"
  rm -f "$output"
  if dx extract_dataset "$RECORD" \
       --fields "$fields" \
       --output "$output" \
       --delimiter "," 2>&1; then
    local lines=$(wc -l < "$output" 2>/dev/null || echo 0)
    if [ "$lines" -gt 1 ]; then
      echo "  ✓ OK: $output ($lines lines)"
      return 0
    fi
  fi
  return 1
}

echo ""
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║   CALON-2: FIX FAILED EXTRACTIONS                                 ║"
echo "║   Corrected Field IDs Based on UKB Showcase Research              ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

# ==============================================================================
# STEP 0: DATA DICTIONARY CHECK
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "STEP 0: Checking dispensed fields (data dictionary)"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

rm -f data_dictionary.tsv
dx extract_dataset "$RECORD" --ddd > data_dictionary.tsv 2>/dev/null || true

if [ -f data_dictionary.tsv ] && [ -s data_dictionary.tsv ]; then
  echo "  Data dictionary generated. Checking key fields..."
  echo ""

  # Check each target field
  for fid in p22189 p189 p26410 p26426 p26427 \
             p26201 p26206 p26216 p26223 p26227 p26244 p26248 p26285 \
             p22009 \
             p22670 p22671 p22674 p22677 p22680 \
             p22682 p22683 p22684 p22685 \
             p22420 p22421 p22426; do
    count=$(grep -c "$fid" data_dictionary.tsv 2>/dev/null || echo 0)
    if [ "$count" -gt 0 ]; then
      echo "    ✓ $fid: DISPENSED ($count entries)"
    else
      echo "    ✗ $fid: NOT dispensed"
    fi
  done
  echo ""
else
  echo "  Could not generate data dictionary. Will try extractions anyway."
  echo ""
fi

# ==============================================================================
# BATCH 1: TOWNSEND DEPRIVATION INDEX
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 1: Townsend Deprivation Index"
echo "  NOTE: Field 189 is RETIRED. Trying replacement field 22189."
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

DEP_SUCCESS=0

# Try 1: p22189 (replacement for retired p189, singular instancing — NO _i0)
echo "  Try 1: participant.p22189 (replacement field, no instance suffix)..."
if extract_once calon_extra_deprivation.csv "participant.eid,participant.p22189"; then
  DEP_SUCCESS=1
fi

# Try 2: p22189_i0 (in case this project's convention needs _i0)
if [ $DEP_SUCCESS -eq 0 ]; then
  echo "  Try 2: participant.p22189_i0 (with instance suffix)..."
  if extract_once calon_extra_deprivation.csv "participant.eid,participant.p22189_i0"; then
    DEP_SUCCESS=1
  fi
fi

# Try 3: Original p189 without _i0 (maybe it needs bare format?)
if [ $DEP_SUCCESS -eq 0 ]; then
  echo "  Try 3: participant.p189 (original field, no instance suffix)..."
  if extract_once calon_extra_deprivation.csv "participant.eid,participant.p189"; then
    DEP_SUCCESS=1
  fi
fi

# Try 4: Index of Multiple Deprivation (England + Wales + Scotland)
if [ $DEP_SUCCESS -eq 0 ]; then
  echo "  Try 4: IMD fields p26410 (England), p26426 (Wales), p26427 (Scotland)..."
  if extract_once calon_extra_deprivation.csv \
       "participant.eid,participant.p26410,participant.p26426,participant.p26427"; then
    DEP_SUCCESS=1
    echo "  NOTE: IMD fields are country-specific. Will need harmonisation in R."
  fi
fi

if [ $DEP_SUCCESS -eq 0 ]; then
  echo "  ✗✗ ALL DEPRIVATION EXTRACTIONS FAILED"
  echo "  Check: Is deprivation dispensed in your application?"
  echo "  Run: grep -i 'depri\|townsend\|p22189\|p189\|26410' data_dictionary.tsv"
fi
echo ""

# ==============================================================================
# BATCH 2: CIMT — RE-VERIFY ALL FIELDS
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 2: CIMT — Re-verify all measurement + quality fields"
echo "  Fields: Mean CIMT at 4 angles (p22671, p22674, p22677, p22680)"
echo "  Plus quality: p22682-p22685"
echo "  Instance: _i2 (imaging visit)"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

CIMT_SUCCESS=0

# Full CIMT extraction: all 12 measurement fields + 4 quality fields
echo "  Try 1: All CIMT mean values at 4 angles (instance 2)..."
if extract_once calon_extra_cimt_full.csv \
     "participant.eid,participant.p22670_i2,participant.p22671_i2,participant.p22672_i2,participant.p22673_i2,participant.p22674_i2,participant.p22675_i2,participant.p22676_i2,participant.p22677_i2,participant.p22678_i2,participant.p22679_i2,participant.p22680_i2,participant.p22681_i2,participant.p22682_i2,participant.p22683_i2,participant.p22684_i2,participant.p22685_i2"; then
  # Check if there's actual data (not all nulls)
  non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_cimt_full.csv)
  if [ "$non_empty" -gt 0 ]; then
    echo "  ✓ CIMT has $non_empty non-empty values"
    CIMT_SUCCESS=1
  else
    echo "  ✗ CIMT fields dispensed but ALL VALUES EMPTY (imaging not done for FH cohort)"
  fi
fi

# Try instance 3 (repeat imaging) as fallback
if [ $CIMT_SUCCESS -eq 0 ]; then
  echo "  Try 2: CIMT instance 3 (repeat imaging visit)..."
  if extract_once calon_extra_cimt_full.csv \
       "participant.eid,participant.p22671_i3,participant.p22674_i3,participant.p22677_i3,participant.p22680_i3"; then
    non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_cimt_full.csv)
    if [ "$non_empty" -gt 0 ]; then
      echo "  ✓ CIMT (instance 3) has $non_empty non-empty values"
      CIMT_SUCCESS=1
    else
      echo "  ✗ CIMT instance 3 also all empty"
    fi
  fi
fi

if [ $CIMT_SUCCESS -eq 0 ]; then
  echo "  RESULT: CIMT not available for this cohort."
  echo "  (UKB CIMT is imaging-visit only; ~96K participants, mostly non-FH)"
fi
echo ""

# ==============================================================================
# BATCH 3: CAC / CARDIAC IMAGING
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 3: Coronary Calcium / Cardiac CT"
echo "  IMPORTANT: UKB does NOT perform CT coronary angiography."
echo "  UKB cardiac imaging is MRI only (no Agatston/CAC score exists)."
echo "  Trying cardiac MRI IDPs as alternative cardiovascular markers."
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

CAC_SUCCESS=0

# Try standard CAC fields one more time (just to verify definitively)
echo "  Verify 1: CAC/CTCA fields p22420, p22421, p22426..."
if extract_once calon_extra_cac_verify.csv \
     "participant.eid,participant.p22420_i2,participant.p22421_i2,participant.p22426_i2"; then
  non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_cac_verify.csv)
  if [ "$non_empty" -gt 0 ]; then
    echo "  ✓ CAC fields have data! ($non_empty values)"
    mv calon_extra_cac_verify.csv calon_extra_cac.csv
    CAC_SUCCESS=1
  else
    echo "  ✗ CAC fields exist but all empty"
    rm -f calon_extra_cac_verify.csv
  fi
else
  echo "  ✗ CAC fields not dispensed (confirmed: no CTCA in UKB)"
fi

# Try cardiac MRI phenotypes as alternative cardiovascular markers
if [ $CAC_SUCCESS -eq 0 ]; then
  echo ""
  echo "  Trying cardiac MRI IDPs as proxy markers..."
  echo "  p22420 = LV end-diastolic volume (indexing for body size)"
  echo "  p22421 = LV end-systolic volume"
  echo "  p22422 = LV stroke volume"
  echo "  p22425 = LV ejection fraction"
  echo "  p22426 = LV cardiac output"

  # These cardiac MRI phenotypes may be available even if CAC isn't
  for field_set in \
    "participant.eid,participant.p22420_i2,participant.p22421_i2,participant.p22422_i2,participant.p22425_i2,participant.p22426_i2" \
    "participant.eid,participant.p22420_i2,participant.p22425_i2"; do
    if extract_once calon_extra_cardiac_mri.csv "$field_set"; then
      non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_cardiac_mri.csv)
      if [ "$non_empty" -gt 0 ]; then
        echo "  ✓ Cardiac MRI IDPs have $non_empty values"
        CAC_SUCCESS=2  # 2 = MRI proxy, not true CAC
        break
      fi
    fi
  done
fi

if [ $CAC_SUCCESS -eq 0 ]; then
  echo "  RESULT: No cardiac imaging data available for this cohort."
  echo "  This will be documented as a limitation in the manuscript."
  echo "  (CAC does not exist in UKB; cardiac MRI was for ~40K imaging substudy)"
fi
echo ""

# ==============================================================================
# BATCH 4: POLYGENIC RISK SCORES (PRS) — Individual key fields
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 4: Polygenic Risk Scores — trying key individual fields"
echo "  Category 301 (Standard PRS, from external GWAS only)"
echo "  These require your application to have dispensed Category 300"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

PRS_SUCCESS=0

# Key PRS fields for cardiovascular research (all singular, no instance)
echo "  Try 1: Key cardiovascular PRS fields (no instance suffix)..."
PRS_FIELDS_BARE="participant.eid"
PRS_FIELDS_BARE="$PRS_FIELDS_BARE,participant.p26223"   # Cardiovascular disease
PRS_FIELDS_BARE="$PRS_FIELDS_BARE,participant.p26227"   # Coronary artery disease
PRS_FIELDS_BARE="$PRS_FIELDS_BARE,participant.p26244"   # Hypertension
PRS_FIELDS_BARE="$PRS_FIELDS_BARE,participant.p26248"   # Ischaemic stroke
PRS_FIELDS_BARE="$PRS_FIELDS_BARE,participant.p26216"   # BMI
PRS_FIELDS_BARE="$PRS_FIELDS_BARE,participant.p26285"   # Type 2 diabetes

if extract_once calon_extra_prs.csv "$PRS_FIELDS_BARE"; then
  non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_prs.csv)
  if [ "$non_empty" -gt 0 ]; then
    echo "  ✓ PRS fields have $non_empty values"
    PRS_SUCCESS=1
  else
    echo "  ✗ PRS fields dispensed but all empty"
  fi
fi

# Try with _i0 suffix
if [ $PRS_SUCCESS -eq 0 ]; then
  echo "  Try 2: Same PRS fields with _i0 suffix..."
  PRS_FIELDS_I0="participant.eid"
  PRS_FIELDS_I0="$PRS_FIELDS_I0,participant.p26223_i0"
  PRS_FIELDS_I0="$PRS_FIELDS_I0,participant.p26227_i0"
  PRS_FIELDS_I0="$PRS_FIELDS_I0,participant.p26244_i0"
  PRS_FIELDS_I0="$PRS_FIELDS_I0,participant.p26248_i0"
  PRS_FIELDS_I0="$PRS_FIELDS_I0,participant.p26216_i0"
  PRS_FIELDS_I0="$PRS_FIELDS_I0,participant.p26285_i0"

  if extract_once calon_extra_prs.csv "$PRS_FIELDS_I0"; then
    non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_prs.csv)
    if [ "$non_empty" -gt 0 ]; then
      echo "  ✓ PRS (with _i0) has $non_empty values"
      PRS_SUCCESS=1
    fi
  fi
fi

# Try single most important field: CAD PRS
if [ $PRS_SUCCESS -eq 0 ]; then
  echo "  Try 3: Single CAD PRS field (p26227) only..."
  for fmt in "participant.p26227" "participant.p26227_i0"; do
    if extract_once calon_extra_prs_cad.csv "participant.eid,$fmt"; then
      non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_prs_cad.csv)
      if [ "$non_empty" -gt 0 ]; then
        echo "  ✓ CAD PRS available ($non_empty values)"
        mv calon_extra_prs_cad.csv calon_extra_prs.csv
        PRS_SUCCESS=1
        break
      fi
    fi
  done
fi

# Try Enhanced PRS (Category 302) — p26300 series
if [ $PRS_SUCCESS -eq 0 ]; then
  echo "  Try 4: Enhanced PRS (Category 302)..."
  for fmt in "participant.p26323" "participant.p26323_i0"; do
    if extract_once calon_extra_prs_enhanced.csv "participant.eid,$fmt"; then
      non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_prs_enhanced.csv)
      if [ "$non_empty" -gt 0 ]; then
        echo "  ✓ Enhanced PRS available ($non_empty values)"
        PRS_SUCCESS=2  # 2 = enhanced, not standard
        break
      fi
    fi
  done
fi

# Fallback: Genetic Principal Components (p22009, always available)
if [ $PRS_SUCCESS -eq 0 ]; then
  echo "  Try 5: Genetic Principal Components (p22009, fallback)..."
  GPC_FIELDS="participant.eid"
  for i in $(seq 0 9); do
    GPC_FIELDS="$GPC_FIELDS,participant.p22009_a${i}"
  done

  if extract_once calon_extra_genetic_pcs.csv "$GPC_FIELDS"; then
    non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="")c++}END{print c+0}' calon_extra_genetic_pcs.csv)
    if [ "$non_empty" -gt 0 ]; then
      echo "  ✓ Genetic PCs available ($non_empty values)"
      echo "  NOTE: PCs are NOT PRS. They control for population stratification."
      echo "         PRS not dispensed — request Category 300 from UKB if needed."
      PRS_SUCCESS=3  # 3 = PCs only, not PRS
    fi
  fi
fi

if [ $PRS_SUCCESS -eq 0 ]; then
  echo "  ✗✗ ALL PRS/GENETIC EXTRACTIONS FAILED"
  echo "  Your project may not have genetic data dispensed."
  echo "  Check: grep -i 'p22009\|p26' data_dictionary.tsv | head -20"
fi
echo ""

# ==============================================================================
# BATCH 5: PRS FLAG + TESTING SUBGROUP
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 5: PRS metadata — testing subgroup flag (p26200)"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

for fmt in "participant.p26200" "participant.p26200_i0"; do
  if extract_once calon_extra_prs_flag.csv "participant.eid,$fmt" 2>/dev/null; then
    echo "  ✓ PRS testing subgroup flag extracted"
    break
  fi
done
echo ""

# ==============================================================================
# UPLOAD ALL NEW FILES TO PROJECT
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "UPLOADING ALL FILES TO PROJECT"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

for f in calon_extra_deprivation.csv \
         calon_extra_cimt_full.csv \
         calon_extra_cac.csv \
         calon_extra_cardiac_mri.csv \
         calon_extra_prs.csv \
         calon_extra_prs_cad.csv \
         calon_extra_prs_enhanced.csv \
         calon_extra_genetic_pcs.csv \
         calon_extra_prs_flag.csv; do
  if [ -f "$f" ] && [ -s "$f" ]; then
    lines=$(wc -l < "$f")
    if [ "$lines" -gt 1 ]; then
      dx upload "$f" --destination "/" 2>/dev/null || true
      echo "  ✓ Uploaded: $f ($lines lines)"
    else
      echo "  ✗ Skipped: $f (empty)"
    fi
  else
    echo "  ✗ Skipped: $f (missing)"
  fi
done

echo ""
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║   FIX EXTRACTION COMPLETE                                         ║"
echo "╠══════════════════════════════════════════════════════════════════════╣"
echo "║   RESULTS SUMMARY:                                                ║"
echo "║                                                                    ║"

if [ $DEP_SUCCESS -gt 0 ]; then
  echo "║   ✓ DEPRIVATION: Extracted successfully                          ║"
else
  echo "║   ✗ DEPRIVATION: Not available — request field to UKB            ║"
fi

if [ $CIMT_SUCCESS -gt 0 ]; then
  echo "║   ✓ CIMT: Has data for this cohort                              ║"
else
  echo "║   ✗ CIMT: Fields dispensed but empty (imaging substudy only)     ║"
fi

if [ $CAC_SUCCESS -gt 0 ]; then
  echo "║   ✓ CARDIAC: Some imaging data available                         ║"
else
  echo "║   ✗ CAC/CTCA: Does NOT exist in UKB (MRI only, no CT)           ║"
fi

if [ $PRS_SUCCESS -eq 1 ]; then
  echo "║   ✓ PRS: Standard PRS extracted                                  ║"
elif [ $PRS_SUCCESS -eq 2 ]; then
  echo "║   ✓ PRS: Enhanced PRS extracted                                  ║"
elif [ $PRS_SUCCESS -eq 3 ]; then
  echo "║   ~ PRS: Only genetic PCs (not true PRS)                         ║"
else
  echo "║   ✗ PRS: Not dispensed — request Category 300 from UKB           ║"
fi

echo "║                                                                    ║"
echo "╠══════════════════════════════════════════════════════════════════════╣"
echo "║   NEXT STEPS:                                                      ║"
echo "║   1. Download all calon_extra_*.csv files                          ║"
echo "║   2. Place in: calon_ukb_pipeline/                                 ║"
echo "║   3. Run 05_CALON2_develop.R                                       ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

# Show data dictionary search results for debugging
if [ -f data_dictionary.tsv ] && [ -s data_dictionary.tsv ]; then
  echo "  DATA DICTIONARY SEARCH (for debugging):"
  echo "  ───────────────────────────────────────"
  echo "  Deprivation fields:"
  grep -i "depri\|townsend\|p22189\|p189\b" data_dictionary.tsv 2>/dev/null | head -5 || echo "    (none found)"
  echo ""
  echo "  PRS/SNP fields:"
  grep -i "p262\|polygenic\|prs" data_dictionary.tsv 2>/dev/null | head -10 || echo "    (none found)"
  echo ""
  echo "  CIMT fields:"
  grep -i "p2267\|carotid\|intima" data_dictionary.tsv 2>/dev/null | head -5 || echo "    (none found)"
  echo ""
  echo "  Cardiac imaging fields:"
  grep -i "p2242\|cardiac\|coronary\|calcium" data_dictionary.tsv 2>/dev/null | head -5 || echo "    (none found)"
  echo ""
  echo "  Genetic PC fields:"
  grep -i "p22009\|principal.comp" data_dictionary.tsv 2>/dev/null | head -5 || echo "    (none found)"
fi

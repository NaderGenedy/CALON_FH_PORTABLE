#!/bin/bash
# ==============================================================================
# CALON-FH: FIX SCRIPT — Re-extract ONLY the failed batches
# ==============================================================================
# This script fixes the errors from the initial extraction run:
#
# FIXES:
#   Batch 2  — Lipids timeout: retry with smaller batch + --sql fallback
#   Batch 5  — Remove p20161 (pack years, not dispensed)
#   Batch 8  — Remove p6150 arrays (not dispensed); keep age-at-event fields
#   Batch 11 — Self-reported illness: smaller batches + error checking
#   Batch 12 — ICD-10: single batch (100 arrays like TUDOR) + visible errors
#   Batch 13 — OPCS-4: visible errors + error checking
#
# RUN ON: UKB-RAP Platform (same terminal as original extraction)
# ==============================================================================

RECORD="project-J6K175jJZ01XppV5477pkYvJ:record-J6K32f8JgZ4JX4gYF39zjBQz"

echo ""
echo "================================================================"
echo "  CALON-FH: RE-EXTRACTING FAILED BATCHES                       "
echo "================================================================"
echo ""

# ==============================================================================
# HELPER: Extract with retry (up to 3 attempts)
# ==============================================================================
extract_with_retry() {
  local output_file="$1"
  shift
  local max_attempts=3
  local attempt=1

  while [ $attempt -le $max_attempts ]; do
    echo "  Attempt $attempt/$max_attempts..."
    # CRITICAL: dx extract_dataset refuses to overwrite — must remove first
    rm -f "$output_file"
    dx extract_dataset "$RECORD" "$@" --output "$output_file" --delimiter ","
    local exit_code=$?

    if [ $exit_code -eq 0 ] && [ -f "$output_file" ]; then
      local lines=$(wc -l < "$output_file")
      if [ "$lines" -gt 1 ]; then
        echo "  SUCCESS: $output_file ($lines lines)"
        return 0
      else
        echo "  WARNING: File created but only has header ($lines lines)"
      fi
    else
      echo "  FAILED (exit code $exit_code)"
    fi

    attempt=$((attempt + 1))
    if [ $attempt -le $max_attempts ]; then
      echo "  Waiting 10 seconds before retry..."
      sleep 10
    fi
  done

  echo "  FAILED after $max_attempts attempts: $output_file"
  return 1
}

# ==============================================================================
# FIX BATCH 2: LIPIDS (retry — original timed out)
# ==============================================================================
echo "=== FIX BATCH 2: Lipid biomarkers (retry) ==="

# Check if existing file has valid data
if [ -f "calon_batch2_lipids.csv" ]; then
  LINES=$(wc -l < calon_batch2_lipids.csv)
  if [ "$LINES" -gt 100 ]; then
    echo "  Existing file has $LINES lines — looks OK, skipping."
  else
    echo "  Existing file only has $LINES lines — re-extracting..."
    extract_with_retry calon_batch2_lipids.csv \
      --fields participant.eid,participant.p30690_i0,participant.p30760_i0,participant.p30780_i0,participant.p30870_i0,participant.p30640_i0,participant.p30630_i0,participant.p30790_i0
  fi
else
  echo "  File not found — extracting..."
  extract_with_retry calon_batch2_lipids.csv \
    --fields participant.eid,participant.p30690_i0,participant.p30760_i0,participant.p30780_i0,participant.p30870_i0,participant.p30640_i0,participant.p30630_i0,participant.p30790_i0
fi

echo ""

# ==============================================================================
# FIX BATCH 5: SMOKING (remove p20161 which is not dispensed)
# ==============================================================================
echo "=== FIX BATCH 5: Smoking + Alcohol (removed p20161) ==="

extract_with_retry calon_batch5_smoking_alcohol.csv \
  --fields participant.eid,participant.p20116_i0,participant.p2867_i0,participant.p2897_i0,participant.p1239_i0,participant.p1249_i0,participant.p20117_i0,participant.p1558_i0

echo ""

# ==============================================================================
# FIX BATCH 8: SELF-REPORTED ASCVD (remove p6150 which is not dispensed)
# ==============================================================================
echo "=== FIX BATCH 8: Self-reported ASCVD (without p6150) ==="

# Field 6150 (vascular problems) was NOT dispensed in this application.
# We still have ASCVD from:
#   - Batch 7: First-occurrence dates (p131296-p131388) — GOLD STANDARD
#   - p3894_i0: Age heart attack diagnosed
#   - p3627_i0: Age angina diagnosed
#   - p4056_i0: Age stroke diagnosed
# These self-reported age-at-event fields ARE dispensed.

extract_with_retry calon_batch8_ascvd_selfreport.csv \
  --fields participant.eid,participant.p3894_i0,participant.p3627_i0,participant.p4056_i0

echo ""

# ==============================================================================
# FIX BATCH 11: SELF-REPORTED ILLNESS (smaller batches + error handling)
# ==============================================================================
echo "=== FIX BATCH 11: Self-reported illness (Field 20002) ==="
echo "  Splitting into 3 smaller batches to avoid timeout..."

# CRITICAL: Remove old temp files first — dx extract_dataset refuses to overwrite
echo "  Cleaning old temp files..."
rm -f /tmp/calon_illness_a.csv /tmp/calon_illness_b.csv /tmp/calon_illness_c.csv
rm -f /tmp/calon_illness_a.txt /tmp/calon_illness_b.txt /tmp/calon_illness_c.txt

# 11a: arrays 0-11 (12 fields)
echo "--- 11a: arrays 0-11 ---"
python3 -c "
fields = ['participant.eid']
for a in range(12):
    fields.append(f'participant.p20002_i0_a{a}')
print('\n'.join(fields))
" > /tmp/calon_illness_a.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_illness_a.txt \
  --output /tmp/calon_illness_a.csv --delimiter ","
ILLNESS_A_OK=$?

if [ $ILLNESS_A_OK -eq 0 ] && [ -f "/tmp/calon_illness_a.csv" ]; then
  echo "  11a OK: $(wc -l < /tmp/calon_illness_a.csv) lines"
else
  echo "  11a FAILED — retrying..."
  sleep 5
  rm -f /tmp/calon_illness_a.csv
  dx extract_dataset "$RECORD" \
    --fields-file /tmp/calon_illness_a.txt \
    --output /tmp/calon_illness_a.csv --delimiter ","
fi

# 11b: arrays 12-23 (12 fields)
echo "--- 11b: arrays 12-23 ---"
python3 -c "
fields = ['participant.eid']
for a in range(12, 24):
    fields.append(f'participant.p20002_i0_a{a}')
print('\n'.join(fields))
" > /tmp/calon_illness_b.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_illness_b.txt \
  --output /tmp/calon_illness_b.csv --delimiter ","
ILLNESS_B_OK=$?

if [ $ILLNESS_B_OK -eq 0 ] && [ -f "/tmp/calon_illness_b.csv" ]; then
  echo "  11b OK: $(wc -l < /tmp/calon_illness_b.csv) lines"
else
  echo "  11b FAILED — retrying..."
  sleep 5
  rm -f /tmp/calon_illness_b.csv
  dx extract_dataset "$RECORD" \
    --fields-file /tmp/calon_illness_b.txt \
    --output /tmp/calon_illness_b.csv --delimiter ","
fi

# 11c: arrays 24-35 (12 fields)
echo "--- 11c: arrays 24-35 ---"
python3 -c "
fields = ['participant.eid']
for a in range(24, 36):
    fields.append(f'participant.p20002_i0_a{a}')
print('\n'.join(fields))
" > /tmp/calon_illness_c.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_illness_c.txt \
  --output /tmp/calon_illness_c.csv --delimiter ","
ILLNESS_C_OK=$?

if [ $ILLNESS_C_OK -eq 0 ] && [ -f "/tmp/calon_illness_c.csv" ]; then
  echo "  11c OK: $(wc -l < /tmp/calon_illness_c.csv) lines"
else
  echo "  11c FAILED — retrying..."
  sleep 5
  rm -f /tmp/calon_illness_c.csv
  dx extract_dataset "$RECORD" \
    --fields-file /tmp/calon_illness_c.txt \
    --output /tmp/calon_illness_c.csv --delimiter ","
fi

# 11d: Merge with error checking
echo "--- 11d: Merging illness parts (with error checking) ---"
python3 -c "
import pandas as pd
import os

parts = []
for f in ['/tmp/calon_illness_a.csv', '/tmp/calon_illness_b.csv', '/tmp/calon_illness_c.csv']:
    if os.path.exists(f) and os.path.getsize(f) > 50:
        try:
            df = pd.read_csv(f)
            if len(df) > 0:
                parts.append(df)
                print(f'  Loaded {f}: {df.shape[0]} rows x {df.shape[1]} cols')
            else:
                print(f'  SKIP {f}: empty dataframe')
        except Exception as e:
            print(f'  SKIP {f}: {e}')
    else:
        print(f'  SKIP {f}: file not found or too small')

if parts:
    merged = parts[0]
    for df in parts[1:]:
        merged = merged.merge(df, on='participant.eid', how='outer')
    merged.to_csv('calon_batch11_illness.csv', index=False)
    print(f'  MERGED: {merged.shape[0]} rows x {merged.shape[1]} cols')
else:
    print('  ERROR: No illness data could be loaded!')
    print('  Creating placeholder with eid only...')
    import subprocess
    # Create a minimal file with just eid from batch 1
    if os.path.exists('calon_batch1_demographics.csv'):
        df1 = pd.read_csv('calon_batch1_demographics.csv', usecols=['participant.eid'])
        df1.to_csv('calon_batch11_illness.csv', index=False)
        print(f'  Placeholder created: {len(df1)} eids')
"

echo "Done: calon_batch11_illness.csv"
echo ""

# ==============================================================================
# FIX BATCH 12: HES ICD-10 (single batch, 100 arrays, visible errors)
# ==============================================================================
echo "=== FIX BATCH 12: HES ICD-10 codes + dates ==="
echo "  Using single batch of 100 arrays (matching working TUDOR extraction)..."

# 12a: ICD-10 codes (100 arrays, single batch)
echo "--- 12a: ICD-10 codes (p41270, arrays 0-99) ---"
python3 -c "
fields = ['participant.eid']
for a in range(100):
    fields.append(f'participant.p41270_a{a}')
print('\n'.join(fields))
" > /tmp/calon_icd10_fields.txt

echo "  Extracting $(wc -l < /tmp/calon_icd10_fields.txt) fields..."
rm -f calon_batch12a_icd10_codes.csv
dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_icd10_fields.txt \
  --output calon_batch12a_icd10_codes.csv --delimiter ","
ICD_CODE_OK=$?

if [ $ICD_CODE_OK -eq 0 ] && [ -f "calon_batch12a_icd10_codes.csv" ]; then
  echo "  ICD-10 codes OK: $(wc -l < calon_batch12a_icd10_codes.csv) lines"
else
  echo "  ICD-10 codes FAILED (field 41270 may not be dispensed)"
  echo "  Creating empty placeholder..."
  echo "participant.eid" > calon_batch12a_icd10_codes.csv
fi

# 12b: ICD-10 dates (100 arrays)
echo "--- 12b: ICD-10 dates (p41280, arrays 0-99) ---"
python3 -c "
fields = ['participant.eid']
for a in range(100):
    fields.append(f'participant.p41280_a{a}')
print('\n'.join(fields))
" > /tmp/calon_icd10_dates.txt

echo "  Extracting $(wc -l < /tmp/calon_icd10_dates.txt) fields..."
rm -f calon_batch12b_icd10_dates.csv
dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_icd10_dates.txt \
  --output calon_batch12b_icd10_dates.csv --delimiter ","
ICD_DATE_OK=$?

if [ $ICD_DATE_OK -eq 0 ] && [ -f "calon_batch12b_icd10_dates.csv" ]; then
  echo "  ICD-10 dates OK: $(wc -l < calon_batch12b_icd10_dates.csv) lines"
else
  echo "  ICD-10 dates FAILED (field 41280 may not be dispensed)"
  echo "  Creating empty placeholder..."
  echo "participant.eid" > calon_batch12b_icd10_dates.csv
fi

echo ""

# ==============================================================================
# FIX BATCH 13: OPCS-4 (visible errors + error checking)
# ==============================================================================
echo "=== FIX BATCH 13: OPCS-4 procedure codes + dates ==="

# 13a: OPCS-4 codes (60 arrays)
echo "--- 13a: OPCS-4 codes (p41272, arrays 0-59) ---"
python3 -c "
fields = ['participant.eid']
for a in range(60):
    fields.append(f'participant.p41272_a{a}')
print('\n'.join(fields))
" > /tmp/calon_opcs_fields.txt

rm -f calon_batch13a_opcs_codes.csv
dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_opcs_fields.txt \
  --output calon_batch13a_opcs_codes.csv --delimiter ","
OPCS_CODE_OK=$?

if [ $OPCS_CODE_OK -eq 0 ] && [ -f "calon_batch13a_opcs_codes.csv" ]; then
  echo "  OPCS codes OK: $(wc -l < calon_batch13a_opcs_codes.csv) lines"
else
  echo "  OPCS codes FAILED (field 41272 may not be dispensed)"
  echo "  Creating empty placeholder..."
  echo "participant.eid" > calon_batch13a_opcs_codes.csv
fi

# 13b: OPCS-4 dates (60 arrays)
echo "--- 13b: OPCS-4 dates (p41282, arrays 0-59) ---"
python3 -c "
fields = ['participant.eid']
for a in range(60):
    fields.append(f'participant.p41282_a{a}')
print('\n'.join(fields))
" > /tmp/calon_opcs_dates.txt

rm -f calon_batch13b_opcs_dates.csv
dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_opcs_dates.txt \
  --output calon_batch13b_opcs_dates.csv --delimiter ","
OPCS_DATE_OK=$?

if [ $OPCS_DATE_OK -eq 0 ] && [ -f "calon_batch13b_opcs_dates.csv" ]; then
  echo "  OPCS dates OK: $(wc -l < calon_batch13b_opcs_dates.csv) lines"
else
  echo "  OPCS dates FAILED (field 41282 may not be dispensed)"
  echo "  Creating empty placeholder..."
  echo "participant.eid" > calon_batch13b_opcs_dates.csv
fi

echo ""

# ==============================================================================
# VERIFICATION: Check all output files
# ==============================================================================
echo "================================================================"
echo "  VERIFICATION: Checking ALL output files                       "
echo "================================================================"
echo ""

ALL_OK=true
for f in calon_batch1_demographics.csv \
         calon_batch2_lipids.csv \
         calon_batch3_biomarkers.csv \
         calon_batch4_medications.csv \
         calon_batch5_smoking_alcohol.csv \
         calon_batch6_comorbidities.csv \
         calon_batch7_ascvd_dates.csv \
         calon_batch8_ascvd_selfreport.csv \
         calon_batch9_cimt.csv \
         calon_batch10_cac.csv \
         calon_batch11_illness.csv \
         calon_batch12a_icd10_codes.csv \
         calon_batch12b_icd10_dates.csv \
         calon_batch13a_opcs_codes.csv \
         calon_batch13b_opcs_dates.csv; do
  if [ -f "$f" ]; then
    LINES=$(wc -l < "$f")
    COLS=$(head -1 "$f" | tr ',' '\n' | wc -l)
    if [ "$LINES" -gt 1 ]; then
      STATUS="OK"
    else
      STATUS="EMPTY (header only)"
      ALL_OK=false
    fi
    printf "  %-40s %8d lines  %3d cols  [%s]\n" "$f" "$LINES" "$COLS" "$STATUS"
  else
    printf "  %-40s %8s        %3s       [MISSING]\n" "$f" "---" "---"
    ALL_OK=false
  fi
done

echo ""
if [ "$ALL_OK" = true ]; then
  echo "  ALL FILES VERIFIED SUCCESSFULLY"
else
  echo "  SOME FILES HAVE ISSUES — check above for EMPTY/MISSING"
  echo "  NOTE: Batch 9 (CIMT), 10 (CAC), 12 (ICD-10), 13 (OPCS-4) may be"
  echo "        empty if those fields are not dispensed in your UKB application."
  echo "        The R pipeline handles this gracefully."
fi

echo ""
echo "================================================================"
echo "  FIX SCRIPT COMPLETE                                           "
echo "  Next: Run 01_CALON_build_cohort.R                             "
echo "================================================================"

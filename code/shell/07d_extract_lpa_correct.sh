#!/bin/bash
# ==============================================================================
# CALON-2 PIPELINE: STEP 07d — RE-EXTRACT Lp(a) WITH CORRECT FIELD MAPPING
# ==============================================================================
# PURPOSE: The initial extraction of p30900 (Lipoprotein(a)) in
#          ukb_reviewer_apob_lpa.csv returned corrupted/mis-encoded values
#          (range 75-2922, ceiling effect at 2922, 100% > 125 nmol/L).
#          This is clearly artefactual — UKB p30900 should have median ~20
#          nmol/L in the general population.
#
#          This script re-extracts Lp(a) using correct field references,
#          including both instances (baseline i0, repeat i1) and the
#          immunoturbidimetric assay results.
#
# UKB FIELD REFERENCE:
#   p30900  Lipoprotein(a) — Immunoturbidimetric, nmol/L
#           Randox Biosciences Ltd, County Antrim, UK
#           Category 17518 (Blood biochemistry)
#           ~390,000 participants with baseline data
#
#   Typical distribution: median ~20, IQR 7-75, range 3.8-189+
#   Clinical thresholds: >50 nmol/L elevated, >125 nmol/L high risk (ESC 2019)
#
# ALSO EXTRACTS:
#   p30790  Lipoprotein(a) in mg/dL (alternative unit, some UKB versions)
#   p30890  Apolipoprotein B (g/L) — for ApoB/Lp(a) co-analysis
#
# RUN ON: UKB-RAP Platform (JupyterLab Bash terminal)
# AUTHOR: Dr Nader Genedy
# DATE:   April 2026
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
        echo "  OK: $output ($lines lines)"
        return 0
      else
        echo "  Empty result"
      fi
    else
      echo "  Failed attempt $attempt"
    fi
    sleep 3
  done
  echo "  FAILED after $max_retries attempts: $output"
  return 1
}

echo ""
echo "================================================================="
echo "  CALON-2: RE-EXTRACT Lp(a) — CORRECT FIELD MAPPING"
echo "  UKB Field p30900 (Immunoturbidimetric, nmol/L)"
echo "================================================================="
echo ""

# ==============================================================================
# EXTRACTION 1: Lp(a) nmol/L — both instances
# ==============================================================================
echo "--- Extraction 1: Lp(a) in nmol/L (p30900) ---"
echo "  Expected: ~390,000 non-null at baseline"
echo ""

LPA_SUCCESS=0
for suffix in "_i0" "" "_i0_a0"; do
  LPA_FIELDS="participant.eid"
  LPA_FIELDS="$LPA_FIELDS,participant.p30900${suffix}"

  # Also try repeat visit
  LPA_FIELDS_FULL="$LPA_FIELDS"
  repeat_suffix="${suffix/i0/i1}"
  if [ "$suffix" != "$repeat_suffix" ]; then
    LPA_FIELDS_FULL="$LPA_FIELDS_FULL,participant.p30900${repeat_suffix}"
  fi

  echo "  Trying suffix '$suffix'..."
  if extract_with_retry calon_extra_lpa_nmol.csv "$LPA_FIELDS_FULL"; then
    # Validate: check value range
    non_empty=$(awk -F',' 'NR>1{for(i=2;i<=NF;i++)if($i!="" && $i+0>0)c++}END{print c+0}' calon_extra_lpa_nmol.csv)
    echo "  Non-empty Lp(a) values: $non_empty"

    if [ "$non_empty" -gt 1000 ]; then
      # Check median is in expected range (should be ~20, NOT ~2900)
      median=$(awk -F',' 'NR>1 && $2!="" && $2+0>0{a[NR]=$2+0; n++} END{
        asort(a); print a[int(n/2)]}' calon_extra_lpa_nmol.csv)
      echo "  Approximate median: $median"

      if [ "$(echo "$median < 500" | bc -l)" -eq 1 ]; then
        echo "  VALIDATION PASSED: median < 500 nmol/L (expected ~20)"
        LPA_SUCCESS=1
        break
      else
        echo "  WARNING: median $median > 500 — likely still mis-encoded"
        echo "  Trying next suffix..."
      fi
    fi
  fi
done

if [ $LPA_SUCCESS -eq 0 ]; then
  echo ""
  echo "  All suffixes failed for p30900. Trying alternative field p30790..."
  # p30790 is Lp(a) in mg/dL in some UKB versions
  for suffix in "_i0" "" "_i0_a0"; do
    ALT_FIELDS="participant.eid,participant.p30790${suffix}"
    echo "  Trying p30790 with suffix '$suffix'..."
    if extract_with_retry calon_extra_lpa_mgdl.csv "$ALT_FIELDS"; then
      non_empty=$(awk -F',' 'NR>1 && $2!="" && $2+0>0{c++}END{print c+0}' calon_extra_lpa_mgdl.csv)
      echo "  Non-empty Lp(a) mg/dL values: $non_empty"
      if [ "$non_empty" -gt 1000 ]; then
        echo "  SUCCESS with p30790 (mg/dL). Convert: nmol/L = mg/dL * 2.15"
        LPA_SUCCESS=2
        break
      fi
    fi
  done
fi

echo ""

# ==============================================================================
# EXTRACTION 2: ApoB (p30890) — fresh extraction for consistency
# ==============================================================================
echo "--- Extraction 2: ApoB (p30890, g/L) ---"
echo ""

APOB_SUCCESS=0
for suffix in "_i0" "" "_i0_a0"; do
  APOB_FIELDS="participant.eid,participant.p30890${suffix}"
  echo "  Trying suffix '$suffix'..."
  if extract_with_retry calon_extra_apob.csv "$APOB_FIELDS"; then
    non_empty=$(awk -F',' 'NR>1 && $2!="" && $2+0>0{c++}END{print c+0}' calon_extra_apob.csv)
    echo "  Non-empty ApoB values: $non_empty"
    if [ "$non_empty" -gt 1000 ]; then
      APOB_SUCCESS=1
      break
    fi
  fi
done

echo ""

# ==============================================================================
# UPLOAD
# ==============================================================================
echo "--- Uploading files ---"
for f in calon_extra_lpa_nmol.csv calon_extra_lpa_mgdl.csv calon_extra_apob.csv; do
  if [ -f "$f" ] && [ -s "$f" ]; then
    lines=$(wc -l < "$f")
    if [ "$lines" -gt 1 ]; then
      dx upload "$f" --destination "/" 2>/dev/null || true
      echo "  Uploaded: $f ($lines lines)"
    fi
  fi
done

echo ""
echo "================================================================="
echo "  EXTRACTION COMPLETE"
echo "================================================================="
echo ""

if [ $LPA_SUCCESS -eq 1 ]; then
  echo "  Lp(a) (nmol/L, p30900): SUCCESS"
elif [ $LPA_SUCCESS -eq 2 ]; then
  echo "  Lp(a) (mg/dL, p30790): SUCCESS — needs conversion (* 2.15)"
else
  echo "  Lp(a): FAILED — check field availability on RAP"
  echo "  DEBUG: Try manual extraction on RAP:"
  echo "    dx extract_dataset $RECORD \\"
  echo "      --fields 'participant.eid,participant.p30900_i0' \\"
  echo "      --output test_lpa.csv --delimiter ','"
  echo "    head -5 test_lpa.csv"
  echo "    awk -F',' 'NR>1 && \$2!=\"\"{c++}END{print c}' test_lpa.csv"
fi

if [ $APOB_SUCCESS -eq 1 ]; then
  echo "  ApoB (g/L, p30890): SUCCESS"
else
  echo "  ApoB: FAILED"
fi

echo ""
echo "  NEXT: Download files, place in calon_ukb_pipeline/"
echo "        Update 53_aortic_phenotype_panel.py to use new Lp(a) data"
echo "================================================================="

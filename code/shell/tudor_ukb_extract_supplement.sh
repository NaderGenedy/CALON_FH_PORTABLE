#!/usr/bin/env bash
##############################################################################
# TUDOR — Supplementary UKB Field Extraction
# Flaws: 5 (ethnicity), 6 (T2DM/BMI/alcohol), 7 (ASCVD), 10 (longit. LDL)
#
# ENVIRONMENT: UKB Research Analysis Platform (RAP) — JupyterLab terminal
#              OR Swiss Army Knife app
#
# HOW TO RUN:
#   1.  Open RAP JupyterLab → File → New → Terminal
#   2.  Upload this file: dx upload tudor_ukb_extract_supplement.sh
#   3.  bash tudor_ukb_extract_supplement.sh
#   4.  dx download tudor_supplement_fields.csv
#   5.  On local machine: Rscript tudor_supplement_process.R
#
# OUTPUT:  tudor_supplement_fields.csv   (one row per participant)
##############################################################################

set -euo pipefail

# ── CONFIGURE THESE TWO LINES ─────────────────────────────────────────────────
# Find your dataset record ID in RAP: Data → Datasets → click your dispensation
# It looks like: record-XXXXXXXXXXXXXXXXXXXXXXXX
DATASET_ID="record-XXXXXXXXXXXXXXXXXXXXXXXX"   # ← REPLACE WITH YOUR DATASET ID

# Path to your restricted lipid-clinic EID list on RAP
# (upload with: dx upload ukb_lipid_clinic_eids.txt)
EID_FILE="ukb_lipid_clinic_eids.txt"   # one EID per line, no header

OUT="tudor_supplement_fields.csv"
# ─────────────────────────────────────────────────────────────────────────────

echo "=================================================="
echo " TUDOR UKB Supplementary Extraction"
echo " Dataset : $DATASET_ID"
echo " EID file: $EID_FILE"
echo " Output  : $OUT"
echo "=================================================="

##############################################################################
# FIELD LIST
# ──────────────────────────────────────────────────────────────────────────
#  F5  — ETHNICITY
#    participant.p21000_i0   Ethnic background (instance 0)
#
#  F6  — TRG SHIELD CONFOUNDERS
#    participant.p2443_i0    Diabetes diagnosed by doctor (1=Yes, 0=No)
#    participant.p21001_i0   BMI kg/m² (baseline)
#    participant.p21001_i1   BMI kg/m² (repeat 1)
#    participant.p21001_i2   BMI kg/m² (repeat 2)
#    participant.p21001_i3   BMI kg/m² (repeat 3)
#    participant.p1558_i0    Alcohol intake frequency
#    participant.p20116_i0   Smoking status
#
#  F7  — ASCVD EVENTS (hospital admissions)
#    participant.p41270       ICD-10 diagnoses  (array, all admissions)
#    participant.p41280       Date of diagnoses (array, paired with p41270)
#    participant.p53_i0       Date of baseline assessment
#    participant.p34_i0       Year of birth
#    participant.p52_i0       Month of birth
#    participant.p31_i0       Sex (1=Male, 0=Female)
#
#  F10 — LONGITUDINAL LDL (multiple assessment instances)
#    participant.p30780_i0   LDL direct mmol/L  (baseline ~2006-10)
#    participant.p30780_i1   LDL direct mmol/L  (repeat 1 ~2012-13)
#    participant.p30780_i2   LDL direct mmol/L  (imaging ~2014+)
#    participant.p30780_i3   LDL direct mmol/L  (repeat imaging)
#    participant.p30690_i0   Total cholesterol  (baseline)
#    participant.p30690_i1   Total cholesterol  (repeat 1)
#    participant.p30690_i2   Total cholesterol  (imaging)
#    participant.p30760_i0   HDL cholesterol    (baseline)
#    participant.p30760_i1   HDL cholesterol    (repeat 1)
#    participant.p30760_i2   HDL cholesterol    (imaging)
#    participant.p30870_i0   Triglycerides      (baseline)
#    participant.p30870_i1   Triglycerides      (repeat 1)
#    participant.p30870_i2   Triglycerides      (imaging)
#    participant.p30890_i0   ApoB               (baseline)
#    participant.p30900_i0   Lipoprotein(a)     (baseline)
##############################################################################

FIELDS="participant.eid,\
participant.p21000_i0,\
participant.p2443_i0,\
participant.p21001_i0,participant.p21001_i1,participant.p21001_i2,participant.p21001_i3,\
participant.p1558_i0,\
participant.p20116_i0,\
participant.p41270,\
participant.p41280,\
participant.p53_i0,\
participant.p34_i0,\
participant.p52_i0,\
participant.p31_i0,\
participant.p30780_i0,participant.p30780_i1,participant.p30780_i2,participant.p30780_i3,\
participant.p30690_i0,participant.p30690_i1,participant.p30690_i2,\
participant.p30760_i0,participant.p30760_i1,participant.p30760_i2,\
participant.p30870_i0,participant.p30870_i1,participant.p30870_i2,\
participant.p30890_i0,\
participant.p30900_i0"

# Remove whitespace from field string
FIELDS=$(echo "$FIELDS" | tr -d ' \t\n')

echo ""
echo "Extracting ${#FIELDS} characters of field definitions..."
echo ""

##############################################################################
# STEP 1 — Extract via dx extract_dataset
##############################################################################

# Option A: extract all participants in dataset, then filter by EID in R
# (use this if --entity-filter is not available in your RAP version)

dx extract_dataset "$DATASET_ID" \
  --fields "$FIELDS" \
  --delimiter "," \
  --output "$OUT"

echo "✓ Raw extraction complete → $OUT"
echo "  Rows: $(wc -l < "$OUT") (including header)"

##############################################################################
# STEP 2 — Filter to your restricted cohort EIDs (if EID file provided)
##############################################################################

if [[ -f "$EID_FILE" ]]; then
  echo ""
  echo "Filtering to restricted lipid-clinic cohort ($EID_FILE)..."

  # Use awk to filter: keep header + rows where column 1 (eid) is in EID_FILE
  awk -F',' 'NR==FNR{eids[$1]=1; next}
             FNR==1{print; next}
             ($1 in eids){print}' \
    "$EID_FILE" "$OUT" > "${OUT%.csv}_filtered.csv"

  N_FILTERED=$(wc -l < "${OUT%.csv}_filtered.csv")
  echo "✓ Filtered file: ${OUT%.csv}_filtered.csv"
  echo "  Rows after filter: $((N_FILTERED - 1)) participants"

  # Use filtered file as main output
  mv "${OUT%.csv}_filtered.csv" "$OUT"
else
  echo "No EID filter file found — keeping all participants."
  echo "(Upload ukb_lipid_clinic_eids.txt to filter to your cohort)"
fi

##############################################################################
# STEP 3 — Upload output back to RAP project storage
##############################################################################

dx upload "$OUT" --destination "/"
echo "✓ Uploaded to RAP: /$OUT"

echo ""
echo "=================================================="
echo " NEXT STEPS"
echo "=================================================="
echo "  On your local machine:"
echo "    dx download $OUT"
echo "    Rscript tudor_supplement_process.R"
echo ""
echo "  Field → variable mapping:"
echo "    p21000_i0  → ethnicity_group (5 groups)"
echo "    p2443_i0   → t2dm (0/1)"
echo "    p21001_*   → bmi (most recent non-missing)"
echo "    p1558_i0   → alcohol_cat (Heavy/Moderate/Light)"
echo "    p41270/80  → premature_ascvd_55, age_first_ascvd"
echo "    p30780_*   → ldl_baseline, ldl_mean, chol_years_ukb"
echo "    p30890_i0  → apob_measured (g/L)"
echo "    p30900_i0  → lpa_nmol_L"
echo "=================================================="

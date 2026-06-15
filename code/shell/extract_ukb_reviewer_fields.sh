#!/bin/bash
# =============================================================================
# CALON-FH: UK Biobank RAP Field Extraction for Reviewer Response
# =============================================================================
# Purpose: Extract missing fields from UKB RAP to address reviewer flaws
# Author:  Dr Nader Genedy / Claude
# Date:    2026-03-29
#
# INSTRUCTIONS:
# 1. Log into UKB RAP: https://ukbiobank.dnanexus.com
# 2. Open a JupyterLab or Cloud Workstation terminal
# 3. Upload this script: dx upload extract_ukb_reviewer_fields.sh
# 4. Run: bash extract_ukb_reviewer_fields.sh
# 5. Download outputs: dx download ukb_reviewer_*.csv
#
# IMPORTANT: Run each extraction separately if memory is limited.
# Each block can be run independently.
# =============================================================================

set -euo pipefail

# --- Configuration ---
# Update this to your UKB dataset ID (find via: dx find data --name "*.dataset" --project)
DATASET_ID="app*.dataset"  # Replace with actual dataset ID
PROJECT=$(dx env --bash | grep DX_PROJECT_CONTEXT_ID | cut -d= -f2 | tr -d '"')
OUTPUT_DIR="./ukb_reviewer_extracts"
mkdir -p "$OUTPUT_DIR"

echo "=================================================="
echo "CALON-FH: UKB RAP Field Extraction"
echo "Project: $PROJECT"
echo "Output:  $OUTPUT_DIR"
echo "=================================================="

# Find the dataset
echo "[1/6] Finding dataset..."
DATASET=$(dx find data --name "*.dataset" --project "$PROJECT" --brief | head -1)
if [ -z "$DATASET" ]; then
    echo "ERROR: No dataset found. Check your project ID."
    echo "Try: dx find data --name '*.dataset'"
    exit 1
fi
echo "Dataset: $DATASET"

# =============================================================================
# EXTRACTION 1: Demographics & Covariates (Flaw 5 + 6)
# Fields: ethnicity, alcohol, birth year/month, baseline date, BMI longitudinal
# =============================================================================
echo ""
echo "[2/6] Extracting demographics & covariates..."

dx extract_dataset "$DATASET" \
    --fields \
        "participant.eid" \
        "participant.p21000" \
        "participant.p21000_i1" \
        "participant.p21000_i2" \
        "participant.p1558_i0" \
        "participant.p1558_i1" \
        "participant.p1558_i2" \
        "participant.p34" \
        "participant.p52" \
        "participant.p53_i0" \
        "participant.p53_i1" \
        "participant.p53_i2" \
        "participant.p53_i3" \
        "participant.p21001_i0" \
        "participant.p21001_i1" \
        "participant.p21001_i2" \
        "participant.p21001_i3" \
    --delimiter "," \
    --output "$OUTPUT_DIR/ukb_reviewer_demographics.csv"

echo "  -> Saved: ukb_reviewer_demographics.csv"

# =============================================================================
# EXTRACTION 2: ICD-10 Diagnoses + Dates (Flaw 7)
# Fields: p41270 (ICD-10 array), p41280 (date array)
# NOTE: These are array fields — each participant can have multiple diagnoses
# =============================================================================
echo ""
echo "[3/6] Extracting ICD-10 diagnoses and dates..."

# ICD-10 codes (array field — up to ~200 entries per participant)
# Extract first 50 array indices (covers vast majority of participants)
ICD_FIELDS="participant.eid"
for i in $(seq 0 49); do
    ICD_FIELDS="$ICD_FIELDS,participant.p41270_a${i}"
done

dx extract_dataset "$DATASET" \
    --fields "$ICD_FIELDS" \
    --delimiter "," \
    --output "$OUTPUT_DIR/ukb_reviewer_icd10_codes.csv"

echo "  -> Saved: ukb_reviewer_icd10_codes.csv"

# ICD-10 dates (matching array)
DATE_FIELDS="participant.eid"
for i in $(seq 0 49); do
    DATE_FIELDS="$DATE_FIELDS,participant.p41280_a${i}"
done

dx extract_dataset "$DATASET" \
    --fields "$DATE_FIELDS" \
    --delimiter "," \
    --output "$OUTPUT_DIR/ukb_reviewer_icd10_dates.csv"

echo "  -> Saved: ukb_reviewer_icd10_dates.csv"

# =============================================================================
# EXTRACTION 3: Longitudinal Lipids (Flaw 10 — Cholesterol-Years)
# Fields: LDL, TC, HDL, TG across all instances (i0-i3)
# =============================================================================
echo ""
echo "[4/6] Extracting longitudinal lipids..."

dx extract_dataset "$DATASET" \
    --fields \
        "participant.eid" \
        "participant.p30780_i0" \
        "participant.p30780_i1" \
        "participant.p30780_i2" \
        "participant.p30780_i3" \
        "participant.p30690_i0" \
        "participant.p30690_i1" \
        "participant.p30690_i2" \
        "participant.p30690_i3" \
        "participant.p30760_i0" \
        "participant.p30760_i1" \
        "participant.p30760_i2" \
        "participant.p30760_i3" \
        "participant.p30870_i0" \
        "participant.p30870_i1" \
        "participant.p30870_i2" \
        "participant.p30870_i3" \
    --delimiter "," \
    --output "$OUTPUT_DIR/ukb_reviewer_longitudinal_lipids.csv"

echo "  -> Saved: ukb_reviewer_longitudinal_lipids.csv"

# =============================================================================
# EXTRACTION 4: ApoB and Lp(a) (Flaw 12)
# p30890 = ApoB (g/L), p30900 = Lp(a) immunoassay (nmol/L)
# =============================================================================
echo ""
echo "[5/6] Extracting ApoB and Lp(a)..."

dx extract_dataset "$DATASET" \
    --fields \
        "participant.eid" \
        "participant.p30890_i0" \
        "participant.p30890_i1" \
        "participant.p30890_i2" \
        "participant.p30890_i3" \
        "participant.p30900_i0" \
        "participant.p30900_i1" \
        "participant.p30900_i2" \
        "participant.p30900_i3" \
    --delimiter "," \
    --output "$OUTPUT_DIR/ukb_reviewer_apob_lpa.csv"

echo "  -> Saved: ukb_reviewer_apob_lpa.csv"

# =============================================================================
# EXTRACTION 5: Smoking status (all instances, for completeness)
# =============================================================================
echo ""
echo "[6/6] Extracting smoking (all instances)..."

dx extract_dataset "$DATASET" \
    --fields \
        "participant.eid" \
        "participant.p20116_i0" \
        "participant.p20116_i1" \
        "participant.p20116_i2" \
        "participant.p20116_i3" \
    --delimiter "," \
    --output "$OUTPUT_DIR/ukb_reviewer_smoking.csv"

echo "  -> Saved: ukb_reviewer_smoking.csv"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=================================================="
echo "EXTRACTION COMPLETE"
echo "=================================================="
echo ""
echo "Files generated in $OUTPUT_DIR/:"
ls -lh "$OUTPUT_DIR/"
echo ""
echo "To download to your local machine:"
echo "  dx download -r $OUTPUT_DIR/"
echo ""
echo "Field mapping for reviewer response:"
echo "  Flaw 5  (Ethnicity)         -> ukb_reviewer_demographics.csv (p21000)"
echo "  Flaw 6  (Covariates)        -> ukb_reviewer_demographics.csv (p1558 alcohol, p21001 BMI)"
echo "  Flaw 7  (Age/Events)        -> ukb_reviewer_demographics.csv (p34, p52, p53)"
echo "                              -> ukb_reviewer_icd10_codes.csv (p41270)"
echo "                              -> ukb_reviewer_icd10_dates.csv (p41280)"
echo "  Flaw 10 (Cholesterol-years) -> ukb_reviewer_longitudinal_lipids.csv"
echo "  Flaw 12 (ApoB/Lp(a))       -> ukb_reviewer_apob_lpa.csv (p30890, p30900)"
echo ""
echo "NEXT STEPS:"
echo "  1. Download CSVs to D:/calon_ukb_pipeline/"
echo "  2. Run the merging script: python 42_merge_reviewer_fields.py"
echo "=================================================="

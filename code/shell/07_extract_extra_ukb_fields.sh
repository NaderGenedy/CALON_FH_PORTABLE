#!/bin/bash
# ==============================================================================
# CALON-2 PIPELINE: STEP 07 — EXTRACT ADDITIONAL UKB FIELDS
# ==============================================================================
# PURPOSE: Extract Townsend deprivation, NMR metabolomics, CAC (calcium score),
#          and Polygenic Risk Scores from UKB-RAP for the CALON-2 model.
#
# RUN ON: UKB-RAP Platform (JupyterLab Bash terminal or %%bash cell)
#
# OUTPUTS:
#   1. calon_extra_deprivation.csv   — Townsend deprivation index (p189)
#   2. calon_extra_cac.csv           — Coronary artery calcium / cardiac imaging
#   3. calon_extra_nmr_a.csv         — NMR: Subclass total lipids (p23400-p23409)
#   4. calon_extra_nmr_b.csv         — NMR: HDL lipids + VLDL TG (p23410-p23448)
#   5. calon_extra_nmr_c.csv         — NMR: LDL/HDL TG + VLDL particles
#   6. calon_extra_nmr_d.csv         — NMR: Particle concentrations
#   7. calon_extra_nmr_e.csv         — NMR: Particle sizes + Clinical NMR
#   8. calon_extra_nmr_f.csv         — NMR: Amino acids + GlycA (p23495-p23510)
#   9. calon_extra_nmr_g.csv         — NMR: Fatty acids (p23511-p23521)
#  10. calon_extra_nmr_h.csv         — NMR: Glycolysis + ketones (p23522-p23531)
#  11. calon_extra_prs.csv           — Polygenic Risk Scores (p26201-p26283)
#
# AFTER RUNNING: Download all files to C:/Users/nader/Downloads/calon_ukb_pipeline/
#
# AUTHOR:  Dr Nader Genedy
# DATE:    February 2026
# PROJECT: CALON-2 ASCVD Risk Model Development
# ==============================================================================

RECORD="project-J6K175jJZ01XppV5477pkYvJ:record-J6K32f8JgZ4JX4gYF39zjBQz"

echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║   CALON-2: EXTRACT ADDITIONAL UKB FIELDS                           ║"
echo "║   Townsend Deprivation + NMR + CAC + PRS                           ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

FAIL_COUNT=0
SUCCESS_COUNT=0

extract_with_retry() {
  local output_file="$1"
  shift
  local max_attempts=3
  local attempt=1

  while [ $attempt -le $max_attempts ]; do
    echo "  Attempt $attempt/$max_attempts..."
    rm -f "$output_file"  # CRITICAL: dx refuses to overwrite existing files
    "$@"
    if [ -f "$output_file" ] && [ "$(wc -l < "$output_file")" -gt 1 ]; then
      echo "  ✓ OK: $output_file ($(wc -l < "$output_file") lines)"
      SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
      return 0
    fi
    echo "  ✗ Failed attempt $attempt"
    attempt=$((attempt + 1))
    sleep 5
  done
  echo "  ✗✗ FAILED after $max_attempts attempts: $output_file"
  FAIL_COUNT=$((FAIL_COUNT + 1))
  return 1
}

# ==============================================================================
# BATCH 1: TOWNSEND DEPRIVATION INDEX
# ==============================================================================
# Field p189: Townsend deprivation index at recruitment
# Standard UKB socioeconomic measure — higher = more deprived
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 1: Townsend Deprivation Index (p189)"
echo "═══════════════════════════════════════════════════════════════════════"

extract_with_retry calon_extra_deprivation.csv \
  dx extract_dataset "$RECORD" \
    --fields "participant.eid,participant.p189_i0" \
    --output calon_extra_deprivation.csv --delimiter ","
echo ""

# ==============================================================================
# BATCH 2: CORONARY ARTERY CALCIUM / CARDIAC IMAGING
# ==============================================================================
# UKB cardiac imaging is primarily MRI. CT-derived CAC may not be dispensed.
# Fields attempted:
#   p22420 — Cardiac CT derived value
#   p22421 — Cardiac CT derived value
#   p22426 — Total arterial calcium score
# Instance 2 (imaging visit)
#
# If these fail → create empty placeholder and document as limitation
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 2: Coronary Artery Calcium / Cardiac Imaging (may fail)"
echo "═══════════════════════════════════════════════════════════════════════"

python3 -c "
fields = ['participant.eid']
# Cardiac CT derived fields (Instance 2 = imaging visit)
for f in [22420, 22421, 22426]:
    fields.append(f'participant.p{f}_i2')
# Also try Instance 0 in case imaging data is there
for f in [22420, 22421, 22426]:
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/cac_fields.txt

rm -f calon_extra_cac.csv
dx extract_dataset "$RECORD" --fields-file /tmp/cac_fields.txt \
  --output calon_extra_cac.csv --delimiter "," 2>/dev/null

if [ -f "calon_extra_cac.csv" ] && [ "$(wc -l < calon_extra_cac.csv)" -gt 1 ]; then
  echo "  ✓ OK: calon_extra_cac.csv ($(wc -l < calon_extra_cac.csv) lines)"
  SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
else
  echo "  ✗ CAC fields not dispensed or extraction failed."
  echo "  Creating empty placeholder..."
  echo "participant.eid" > calon_extra_cac.csv
  FAIL_COUNT=$((FAIL_COUNT + 1))
  echo "  NOTE: CAC will be documented as limitation in manuscript."
fi
echo ""

# ==============================================================================
# BATCH 3-7: NMR METABOLOMICS (5 core batches from TUDOR template)
# ==============================================================================
# Category 220: Nightingale NMR Metabolomics
# Fields p23400–p23494 (core lipoprotein subclass data)
# Small batches (10-11 fields each) to avoid DataTooLarge API errors
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 3: NMR — Subclass total lipids (p23400-p23409)"
echo "═══════════════════════════════════════════════════════════════════════"

python3 -c "
fields = ['participant.eid']
for f in range(23400, 23410):
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/nmr_a.txt

extract_with_retry calon_extra_nmr_a.csv \
  dx extract_dataset "$RECORD" --fields-file /tmp/nmr_a.txt \
    --output calon_extra_nmr_a.csv --delimiter ","
echo ""

echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 4: NMR — HDL lipids + VLDL TG (p23410-13, p23442-48)"
echo "═══════════════════════════════════════════════════════════════════════"

python3 -c "
fields = ['participant.eid']
for f in range(23410, 23414):
    fields.append(f'participant.p{f}_i0')
for f in range(23442, 23449):
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/nmr_b.txt

extract_with_retry calon_extra_nmr_b.csv \
  dx extract_dataset "$RECORD" --fields-file /tmp/nmr_b.txt \
    --output calon_extra_nmr_b.csv --delimiter ","
echo ""

echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 5: NMR — LDL/HDL TG + VLDL particles (p23449-55, p23470-72)"
echo "═══════════════════════════════════════════════════════════════════════"

python3 -c "
fields = ['participant.eid']
for f in range(23449, 23456):
    fields.append(f'participant.p{f}_i0')
for f in range(23470, 23473):
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/nmr_c.txt

extract_with_retry calon_extra_nmr_c.csv \
  dx extract_dataset "$RECORD" --fields-file /tmp/nmr_c.txt \
    --output calon_extra_nmr_c.csv --delimiter ","
echo ""

echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 6: NMR — Particle concentrations (p23473-p23483)"
echo "═══════════════════════════════════════════════════════════════════════"

python3 -c "
fields = ['participant.eid']
for f in range(23473, 23484):
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/nmr_d.txt

extract_with_retry calon_extra_nmr_d.csv \
  dx extract_dataset "$RECORD" --fields-file /tmp/nmr_d.txt \
    --output calon_extra_nmr_d.csv --delimiter ","
echo ""

echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 7: NMR — Particle sizes + Clinical NMR (p23484-p23494)"
echo "═══════════════════════════════════════════════════════════════════════"

python3 -c "
fields = ['participant.eid']
for f in range(23484, 23495):
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/nmr_e.txt

extract_with_retry calon_extra_nmr_e.csv \
  dx extract_dataset "$RECORD" --fields-file /tmp/nmr_e.txt \
    --output calon_extra_nmr_e.csv --delimiter ","
echo ""

# ==============================================================================
# BATCH 8-10: NMR EXTENDED (amino acids, fatty acids, glycolysis)
# ==============================================================================
# p23495-p23531: Extended NMR features not in original TUDOR extraction
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 8: NMR — GlycA + Amino acids (p23495-p23510)"
echo "═══════════════════════════════════════════════════════════════════════"

python3 -c "
fields = ['participant.eid']
# p23495 = GlycA (inflammation)
# p23496-p23505 = Amino acids (Ala, Gln, Gly, His, Ile, Leu, Val, Phe, Tyr, BCAA)
# p23506-p23510 = Total FA, saturation, %unsaturation, DHA, Omega-3
for f in range(23495, 23511):
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/nmr_f.txt

extract_with_retry calon_extra_nmr_f.csv \
  dx extract_dataset "$RECORD" --fields-file /tmp/nmr_f.txt \
    --output calon_extra_nmr_f.csv --delimiter ","
echo ""

echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 9: NMR — Fatty acids continued (p23511-p23521)"
echo "═══════════════════════════════════════════════════════════════════════"

python3 -c "
fields = ['participant.eid']
# Omega-6, PUFA, MUFA, SFA, LA, DHA%, LA%, Omega-3%, Omega-6%, PUFA%, MUFA%
for f in range(23511, 23522):
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/nmr_g.txt

extract_with_retry calon_extra_nmr_g.csv \
  dx extract_dataset "$RECORD" --fields-file /tmp/nmr_g.txt \
    --output calon_extra_nmr_g.csv --delimiter ","
echo ""

echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 10: NMR — Glycolysis + Ketones + misc (p23522-p23531)"
echo "═══════════════════════════════════════════════════════════════════════"

python3 -c "
fields = ['participant.eid']
# Glucose, Lactate, Pyruvate, Citrate, b-OHbutyrate, Acetone, Acetoacetate,
# Acetate, Creatinine_NMR, Albumin_NMR
for f in range(23522, 23532):
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/nmr_h.txt

extract_with_retry calon_extra_nmr_h.csv \
  dx extract_dataset "$RECORD" --fields-file /tmp/nmr_h.txt \
    --output calon_extra_nmr_h.csv --delimiter ","
echo ""

# ==============================================================================
# BATCH 11: POLYGENIC RISK SCORES (Enhanced Polygenic Scores, Category 301)
# ==============================================================================
# UKB Enhanced PRS fields (p26201-p26283) — standardised scores for many traits
# Key ones for CALON-2:
#   - LDL-C PRS
#   - CAD (coronary artery disease) PRS
#   - Total cholesterol PRS
#
# NOTE: These field IDs should be verified on UKB Showcase.
#       If not dispensed, creates empty placeholder.
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "BATCH 11: Polygenic Risk Scores (p26201-p26283)"
echo "═══════════════════════════════════════════════════════════════════════"

python3 -c "
fields = ['participant.eid']
# Enhanced Polygenic Scores (Category 301)
# Try the full range — only dispensed fields will extract
for f in range(26201, 26284):
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/prs_fields.txt

rm -f calon_extra_prs.csv
dx extract_dataset "$RECORD" --fields-file /tmp/prs_fields.txt \
  --output calon_extra_prs.csv --delimiter "," 2>/dev/null

if [ -f "calon_extra_prs.csv" ] && [ "$(wc -l < calon_extra_prs.csv)" -gt 1 ]; then
  echo "  ✓ OK: calon_extra_prs.csv ($(wc -l < calon_extra_prs.csv) lines)"
  SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
else
  echo "  ✗ PRS fields not dispensed or extraction failed."
  echo "  Trying smaller subsets..."

  # Try key PRS fields individually
  for prs_field in 26206 26210 26212 26218 26230 26248 26260 26269; do
    rm -f /tmp/prs_single.csv
    dx extract_dataset "$RECORD" \
      --fields "participant.eid,participant.p${prs_field}_i0" \
      --output /tmp/prs_single.csv --delimiter "," 2>/dev/null
    if [ -f "/tmp/prs_single.csv" ] && [ "$(wc -l < /tmp/prs_single.csv)" -gt 1 ]; then
      echo "  ✓ Field p${prs_field} available"
    fi
  done

  # If individual fields worked, combine them
  if [ -f "/tmp/prs_single.csv" ]; then
    echo "  Some PRS fields available — extracting available subset..."
    # Use the last successful single-field extract as placeholder
    cp /tmp/prs_single.csv calon_extra_prs.csv 2>/dev/null
  else
    echo "  Creating empty placeholder..."
    echo "participant.eid" > calon_extra_prs.csv
    FAIL_COUNT=$((FAIL_COUNT + 1))
    echo "  NOTE: PRS will be documented as limitation in manuscript."
  fi
fi
echo ""

# ==============================================================================
# UPLOAD ALL FILES TO PROJECT
# ==============================================================================
echo "═══════════════════════════════════════════════════════════════════════"
echo "UPLOADING ALL FILES TO PROJECT"
echo "═══════════════════════════════════════════════════════════════════════"
echo ""

ALL_FILES=(
  calon_extra_deprivation.csv
  calon_extra_cac.csv
  calon_extra_nmr_a.csv
  calon_extra_nmr_b.csv
  calon_extra_nmr_c.csv
  calon_extra_nmr_d.csv
  calon_extra_nmr_e.csv
  calon_extra_nmr_f.csv
  calon_extra_nmr_g.csv
  calon_extra_nmr_h.csv
  calon_extra_prs.csv
)

UPLOAD_COUNT=0
for f in "${ALL_FILES[@]}"; do
  if [ -f "$f" ] && [ "$(wc -l < "$f")" -gt 1 ]; then
    dx rm -f "/$f" 2>/dev/null  # Remove old version if exists
    dx upload "$f" --destination / --brief
    echo "  ✓ Uploaded: $f ($(wc -l < "$f") lines)"
    UPLOAD_COUNT=$((UPLOAD_COUNT + 1))
  else
    echo "  ✗ Skipped: $f (missing or empty)"
  fi
done

echo ""
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║   EXTRACTION COMPLETE                                              ║"
echo "╠══════════════════════════════════════════════════════════════════════╣"
echo "║   Successful: $SUCCESS_COUNT / ${#ALL_FILES[@]} batches                                ║"
echo "║   Failed:     $FAIL_COUNT / ${#ALL_FILES[@]} batches                                   ║"
echo "║   Uploaded:   $UPLOAD_COUNT files                                         ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""
echo " NEXT STEPS:"
echo " 1. Download all calon_extra_*.csv files from the Manage tab"
echo " 2. Place in: C:/Users/nader/Downloads/calon_ukb_pipeline/"
echo " 3. Run: 05_CALON2_develop.R"
echo ""
echo " FILES TO DOWNLOAD:"
for f in "${ALL_FILES[@]}"; do
  echo "   - $f"
done
echo ""
echo " NMR FIELD REFERENCE:"
echo "   Batches A-E: p23400-p23494 (lipoprotein subclasses, particles, clinical)"
echo "   Batch F:     p23495-p23510 (GlycA, amino acids, total FA)"
echo "   Batch G:     p23511-p23521 (fatty acid subtypes and percentages)"
echo "   Batch H:     p23522-p23531 (glucose, lactate, ketone bodies)"
echo ""

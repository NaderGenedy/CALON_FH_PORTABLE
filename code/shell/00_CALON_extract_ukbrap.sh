#!/bin/bash
# ==============================================================================
# CALON-FH PIPELINE: STEP 00 — UKB-RAP DATA EXTRACTION
# ==============================================================================
# PURPOSE: Extract ALL fields needed for CALON external validation from UKB-RAP.
#          CALON predicts ASCVD in genetically confirmed FH patients.
#          This is SEPARATE from TUDOR (which diagnoses FH).
#
# COMPARATOR: SAFEHEART-RE / SAFEHEART-UK calculator
#
# RUN ON: UKB-RAP Platform (JupyterLab Bash terminal or %%bash cell)
#
# OUTPUTS: 12 CSV files uploaded to project root for download.
#
# AFTER RUNNING: Download all files to C:/Users/nader/Downloads/calon_ukb_pipeline/data/
#
# AUTHOR: Dr Nader Genedy
# DATE:   February 2026
# PROJECT: CALON-FH External Validation (Nature-calibre manuscript)
# ==============================================================================

RECORD="project-J6K175jJZ01XppV5477pkYvJ:record-J6K32f8JgZ4JX4gYF39zjBQz"

echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║   CALON-FH: UK BIOBANK EXTRACTION FOR EXTERNAL VALIDATION          ║"
echo "║   Predicting ASCVD in Genetically Confirmed FH                     ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

# ==============================================================================
# BATCH 1: CORE DEMOGRAPHICS + ANTHROPOMETRICS + BLOOD PRESSURE
# ==============================================================================
# Fields:
#   p31        — Sex (0=Female, 1=Male)
#   p21022     — Age at recruitment (years)
#   p34        — Year of birth
#   p52        — Month of birth
#   p53_i0     — Date of attending assessment centre (Visit 0)
#   p53_i1     — Date of attending assessment centre (Visit 1, if repeat)
#   p21001_i0  — BMI (kg/m²) Visit 0
#   p21001_i1  — BMI (kg/m²) Visit 1
#   p48_i0     — Waist circumference (cm) Visit 0
#   p48_i1     — Waist circumference (cm) Visit 1
#   p49_i0     — Hip circumference (cm) Visit 0
#   p4080_i0_a0 — Systolic BP automated reading 1
#   p4080_i0_a1 — Systolic BP automated reading 2
#   p4079_i0_a0 — Diastolic BP automated reading 1
#   p4079_i0_a1 — Diastolic BP automated reading 2
#   p22009_a1-a10 — Genetic principal components (population structure)
#   p22006     — Genetic ethnic grouping (White British confirmed)
#   p22001     — Genetic sex
# ==============================================================================
echo "=== BATCH 1: Core demographics + anthropometrics + BP ==="
dx extract_dataset "$RECORD" \
  --fields participant.eid,\
participant.p31,\
participant.p21022,\
participant.p34,\
participant.p52,\
participant.p53_i0,\
participant.p53_i1,\
participant.p21001_i0,\
participant.p21001_i1,\
participant.p48_i0,\
participant.p48_i1,\
participant.p49_i0,\
participant.p4080_i0_a0,\
participant.p4080_i0_a1,\
participant.p4079_i0_a0,\
participant.p4079_i0_a1,\
participant.p22006,\
participant.p22001 \
  --output calon_batch1_demographics.csv --delimiter ","
echo "Done: calon_batch1_demographics.csv"

# ==============================================================================
# BATCH 2: LIPID BIOMARKERS + ApoB + Lp(a) + ApoA1
# ==============================================================================
# Fields:
#   p30690_i0  — Total cholesterol (mmol/L)
#   p30760_i0  — HDL cholesterol (mmol/L)
#   p30780_i0  — LDL direct (mmol/L)
#   p30870_i0  — Triglycerides (mmol/L)
#   p30640_i0  — Apolipoprotein B (g/L) ← CRITICAL for CALON ApoB/LDL ratio
#   p30630_i0  — Apolipoprotein A (g/L) ← For inverse ApoA1
#   p30790_i0  — Lipoprotein(a) (nmol/L) ← CRITICAL for CALON Lp(a) binary
# ==============================================================================
echo "=== BATCH 2: Lipid biomarkers (TC, HDL, LDL, TG, ApoB, ApoA, Lp(a)) ==="
dx extract_dataset "$RECORD" \
  --fields participant.eid,\
participant.p30690_i0,\
participant.p30760_i0,\
participant.p30780_i0,\
participant.p30870_i0,\
participant.p30640_i0,\
participant.p30630_i0,\
participant.p30790_i0 \
  --output calon_batch2_lipids.csv --delimiter ","
echo "Done: calon_batch2_lipids.csv"

# ==============================================================================
# BATCH 3: OTHER BIOMARKERS (Glucose, HbA1c, CRP, Liver, Renal)
# ==============================================================================
# Fields:
#   p30750_i0  — HbA1c (mmol/mol) ← Diabetes assessment
#   p30740_i0  — Glucose (mmol/L)
#   p30710_i0  — C-reactive protein (mg/L) ← Inflammation
#   p30620_i0  — ALT (U/L) ← Liver function
#   p30650_i0  — AST (U/L)
#   p30680_i0  — Calcium (mmol/L)
#   p30700_i0  — Creatinine (umol/L) ← Renal function / eGFR
#   p30720_i0  — Cystatin C (mg/L) ← Better eGFR
# ==============================================================================
echo "=== BATCH 3: Other biomarkers (HbA1c, Glucose, CRP, Liver, Renal) ==="
dx extract_dataset "$RECORD" \
  --fields participant.eid,\
participant.p30750_i0,\
participant.p30740_i0,\
participant.p30710_i0,\
participant.p30620_i0,\
participant.p30650_i0,\
participant.p30680_i0,\
participant.p30700_i0,\
participant.p30720_i0 \
  --output calon_batch3_biomarkers.csv --delimiter ","
echo "Done: calon_batch3_biomarkers.csv"

# ==============================================================================
# BATCH 4: MEDICATIONS — Visit 0 (arrays 0-47)
# ==============================================================================
# Field 20003: Treatment/medication code (touchscreen)
# Split into 3 parts to avoid DataTooLarge
#
# CRITICAL STATIN CODES for CALON reverse-engineering:
#   1140861958 = Simvastatin
#   1140888648 = Atorvastatin
#   1140910632 = Rosuvastatin
#   1141146234 = Pravastatin
#   1141192414 = Fluvastatin
#   1140861922 = Lovastatin
#   1141146138 = Ezetimibe
#   1140888594 = Fenofibrate
#   1140861924 = Bezafibrate
#
# NOTE: UKB does NOT record statin DOSE — must assume moderate intensity
#       unless evidence of high-intensity statin name
# ==============================================================================

# 4a. Medications arrays 0-15
echo "=== BATCH 4a: Medications V0 arrays 0-15 ==="
python3 -c "
fields = ['participant.eid']
for a in range(16):
    fields.append(f'participant.p20003_i0_a{a}')
print('\n'.join(fields))
" > /tmp/calon_meds_a.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_meds_a.txt \
  --output /tmp/calon_meds_v0_a.csv --delimiter ","
echo "Done: /tmp/calon_meds_v0_a.csv"

# 4b. Medications arrays 16-31
echo "=== BATCH 4b: Medications V0 arrays 16-31 ==="
python3 -c "
fields = ['participant.eid']
for a in range(16, 32):
    fields.append(f'participant.p20003_i0_a{a}')
print('\n'.join(fields))
" > /tmp/calon_meds_b.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_meds_b.txt \
  --output /tmp/calon_meds_v0_b.csv --delimiter ","
echo "Done: /tmp/calon_meds_v0_b.csv"

# 4c. Medications arrays 32-47
echo "=== BATCH 4c: Medications V0 arrays 32-47 ==="
python3 -c "
fields = ['participant.eid']
for a in range(32, 48):
    fields.append(f'participant.p20003_i0_a{a}')
print('\n'.join(fields))
" > /tmp/calon_meds_c.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_meds_c.txt \
  --output /tmp/calon_meds_v0_c.csv --delimiter ","
echo "Done: /tmp/calon_meds_v0_c.csv"

# 4d. Merge medication parts
echo "=== BATCH 4d: Merging medication files ==="
python3 -c "
import pandas as pd
a = pd.read_csv('/tmp/calon_meds_v0_a.csv')
b = pd.read_csv('/tmp/calon_meds_v0_b.csv')
c = pd.read_csv('/tmp/calon_meds_v0_c.csv')
merged = a.merge(b, on='participant.eid', how='outer').merge(c, on='participant.eid', how='outer')
merged.to_csv('calon_batch4_medications.csv', index=False)
print(f'Merged: {merged.shape[0]} rows x {merged.shape[1]} cols')
"
echo "Done: calon_batch4_medications.csv"

# ==============================================================================
# BATCH 5: SMOKING + ALCOHOL
# ==============================================================================
# Fields:
#   p20116_i0  — Smoking status (0=Never, 1=Previous, 2=Current)
#   p20161     — Pack years of smoking
#   p2867_i0   — Age started smoking in current smokers
#   p2897_i0   — Age stopped smoking
#   p1239_i0   — Current tobacco smoking (touchscreen)
#   p1249_i0   — Past tobacco smoking
#   p20117_i0  — Alcohol drinker status
#   p1558_i0   — Alcohol intake frequency
# ==============================================================================
echo "=== BATCH 5: Smoking + Alcohol ==="
dx extract_dataset "$RECORD" \
  --fields participant.eid,\
participant.p20116_i0,\
participant.p20161,\
participant.p2867_i0,\
participant.p2897_i0,\
participant.p1239_i0,\
participant.p1249_i0,\
participant.p20117_i0,\
participant.p1558_i0 \
  --output calon_batch5_smoking_alcohol.csv --delimiter ","
echo "Done: calon_batch5_smoking_alcohol.csv"

# ==============================================================================
# BATCH 6: DIABETES + HYPERTENSION + COMORBIDITIES
# ==============================================================================
# Fields:
#   p2443_i0   — Diabetes diagnosed by doctor (0=No, 1=Yes, -1=DK, -3=Prefer not)
#   p6177_i0   — Medication for BP/cholesterol/diabetes (Male)
#   p6153_i0   — Medication for BP/cholesterol/diabetes (Female)
#   p2966_i0   — Age high blood pressure diagnosed
#   p4080_i0_a0 — Systolic BP (already in Batch 1 but needed for HTN flag)
# ==============================================================================
echo "=== BATCH 6: Diabetes + Hypertension ==="
dx extract_dataset "$RECORD" \
  --fields participant.eid,\
participant.p2443_i0,\
participant.p6177_i0,\
participant.p6153_i0,\
participant.p2966_i0 \
  --output calon_batch6_comorbidities.csv --delimiter ","
echo "Done: calon_batch6_comorbidities.csv"

# ==============================================================================
# BATCH 7: ASCVD OUTCOMES — FIRST OCCURRENCE DATES (ICD-10)
# ==============================================================================
# These are the GOLD STANDARD outcome fields for time-to-event analysis.
# Each field = date when ICD-10 code FIRST appeared in ANY data source
# (HES, death registry, primary care, self-report).
#
# ASCVD COMPOSITE:
#   p131296 — Date I20 first reported (Angina pectoris)
#   p131298 — Date I21 first reported (Acute MI)
#   p131300 — Date I22 first reported (Subsequent MI)
#   p131306 — Date I25 first reported (Chronic IHD)
#   p131364 — Date I63 first reported (Cerebral infarction / ischaemic stroke)
#   p131366 — Date I64 first reported (Stroke NOS)
#   p131380 — Date I70 first reported (Atherosclerosis)
#   p131382 — Date I71 first reported (Aortic aneurysm and dissection)
#   p131386 — Date I73 first reported (Other PVD)
#   p131388 — Date I74 first reported (Arterial embolism & thrombosis)
#
# DEATH:
#   p40000_i0 — Date of death (instance 0)
#   p40001_i0 — Primary cause of death ICD-10
#   p40007_i0 — Age at death
#
# SELF-REPORTED ASCVD:
#   p6150_i0   — Vascular/heart problems diagnosed (arrays)
#   p3894_i0   — Age heart attack diagnosed
#   p3627_i0   — Age angina diagnosed
#   p4056_i0   — Age stroke diagnosed
# ==============================================================================
echo "=== BATCH 7: ASCVD first-occurrence dates ==="
dx extract_dataset "$RECORD" \
  --fields participant.eid,\
participant.p131296,\
participant.p131298,\
participant.p131300,\
participant.p131306,\
participant.p131364,\
participant.p131366,\
participant.p131380,\
participant.p131382,\
participant.p131386,\
participant.p131388,\
participant.p40000_i0,\
participant.p40001_i0,\
participant.p40007_i0 \
  --output calon_batch7_ascvd_dates.csv --delimiter ","
echo "Done: calon_batch7_ascvd_dates.csv"

# ==============================================================================
# BATCH 8: SELF-REPORTED ASCVD HISTORY + AGE AT EVENT
# ==============================================================================
echo "=== BATCH 8: Self-reported ASCVD + age at event ==="
python3 -c "
fields = ['participant.eid']
# p6150 vascular problems (4 arrays)
for a in range(4):
    fields.append(f'participant.p6150_i0_a{a}')
# Age at MI, angina, stroke
fields += ['participant.p3894_i0', 'participant.p3627_i0', 'participant.p4056_i0']
print('\n'.join(fields))
" > /tmp/calon_ascvd_sr.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_ascvd_sr.txt \
  --output calon_batch8_ascvd_selfreport.csv --delimiter ","
echo "Done: calon_batch8_ascvd_selfreport.csv"

# ==============================================================================
# BATCH 9: CAROTID IMT (Imaging Subset ~100k)
# ==============================================================================
# Carotid intima-media thickness — available for imaging participants only.
# Mean IMT across multiple angles for left and right carotid arteries.
#
# Fields 22671-22680: Mean carotid IMT at various angles
# Field 22681-22694: Max carotid IMT at various angles
#
# NOTE: These fields may not be dispensed. Extraction may fail.
#       If so, CIMT analysis will be excluded (acknowledged limitation).
# ==============================================================================
echo "=== BATCH 9: Carotid IMT (imaging subset, may fail if not dispensed) ==="
dx extract_dataset "$RECORD" \
  --fields participant.eid,\
participant.p22671_i2,\
participant.p22672_i2,\
participant.p22673_i2,\
participant.p22674_i2,\
participant.p22675_i2,\
participant.p22676_i2,\
participant.p22677_i2,\
participant.p22678_i2 \
  --output calon_batch9_cimt.csv --delimiter "," 2>/dev/null

if [ -f "calon_batch9_cimt.csv" ] && [ "$(wc -l < calon_batch9_cimt.csv)" -gt 1 ]; then
  echo "Done: calon_batch9_cimt.csv"
else
  echo "NOTE: Carotid IMT fields not dispensed. Creating empty placeholder."
  echo "participant.eid" > calon_batch9_cimt.csv
fi

# ==============================================================================
# BATCH 10: CORONARY ARTERY CALCIUM (CAC) SCORE
# ==============================================================================
# UKB Cardiac Imaging includes cardiac MRI. CAC scoring from CT is NOT part
# of the standard UKB imaging protocol. However, derived imaging variables
# may include coronary calcium from dedicated CT sub-studies.
#
# Check Category 157 (Imaging) and derived fields.
# Field 22420-22421: may contain cardiac CT derived values
# Field 22426: Total arterial calcium score (if available)
#
# NOTE: CAC is LIMITED in UKB. This extraction may fail.
# ==============================================================================
echo "=== BATCH 10: Coronary Calcium Score (limited availability) ==="
# Attempt extraction — may fail if fields not dispensed
dx extract_dataset "$RECORD" \
  --fields participant.eid,\
participant.p22420_i2,\
participant.p22421_i2,\
participant.p22426_i2 \
  --output calon_batch10_cac.csv --delimiter "," 2>/dev/null

if [ -f "calon_batch10_cac.csv" ] && [ "$(wc -l < calon_batch10_cac.csv)" -gt 1 ]; then
  echo "Done: calon_batch10_cac.csv"
else
  echo "NOTE: CAC fields not dispensed or unavailable. Creating empty placeholder."
  echo "NOTE: UKB imaging is primarily cardiac MRI, NOT CT for calcium scoring."
  echo "NOTE: CAC analysis will be noted as limitation in manuscript."
  echo "participant.eid" > calon_batch10_cac.csv
fi

# ==============================================================================
# BATCH 11: SELF-REPORTED ILLNESSES (for comorbidity phenotyping)
# ==============================================================================
# Field 20002: Non-cancer illness codes (arrays 0-35)
# Critical codes for CALON:
#   1065 = Hypertension
#   1075 = Heart attack / MI
#   1066 = Heart failure
#   1067 = Angina
#   1081 = Stroke
#   1082 = TIA
#   1083 = Peripheral vascular disease
#   1223 = Type 2 diabetes
#   1220 = Diabetes (unspecified)
#   1226 = Hypothyroidism
#   1261 = DVT
#   1473 = High cholesterol
#   1074 = AF / atrial fibrillation
# ==============================================================================
echo "=== BATCH 11a: Self-reported illness arrays 0-17 ==="
python3 -c "
fields = ['participant.eid']
for a in range(18):
    fields.append(f'participant.p20002_i0_a{a}')
print('\n'.join(fields))
" > /tmp/calon_illness_a.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_illness_a.txt \
  --output /tmp/calon_illness_a.csv --delimiter ","

echo "=== BATCH 11b: Self-reported illness arrays 18-35 ==="
python3 -c "
fields = ['participant.eid']
for a in range(18, 36):
    fields.append(f'participant.p20002_i0_a{a}')
print('\n'.join(fields))
" > /tmp/calon_illness_b.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_illness_b.txt \
  --output /tmp/calon_illness_b.csv --delimiter ","

# Merge
echo "=== BATCH 11c: Merging illness parts ==="
python3 -c "
import pandas as pd
a = pd.read_csv('/tmp/calon_illness_a.csv')
b = pd.read_csv('/tmp/calon_illness_b.csv')
merged = a.merge(b, on='participant.eid', how='outer')
merged.to_csv('calon_batch11_illness.csv', index=False)
print(f'Merged: {merged.shape[0]} rows x {merged.shape[1]} cols')
"
echo "Done: calon_batch11_illness.csv"

# ==============================================================================
# BATCH 12: HES ICD-10 DIAGNOSES (Hospital Episode Statistics)
# ==============================================================================
# Field 41270: Diagnoses - ICD10 (up to 250+ arrays)
# This provides HOSPITAL-CODED ASCVD events, the gold standard for outcomes.
# Extract first 150 arrays (covers >99% of patients).
#
# CRITICAL ICD-10 codes for ASCVD composite:
#   I20   — Angina pectoris
#   I21   — Acute myocardial infarction
#   I22   — Subsequent MI
#   I23   — Complications of acute MI
#   I24   — Other acute IHD
#   I25   — Chronic IHD
#   I63   — Cerebral infarction
#   I64   — Stroke, not specified
#   I65   — Occlusion/stenosis precerebral arteries
#   I66   — Occlusion/stenosis cerebral arteries
#   I70   — Atherosclerosis
#   I71   — Aortic aneurysm
#   I73   — Other PVD
#   I74   — Arterial embolism
#
# Also need corresponding dates:
# Field 41280: Dates of diagnoses (same array structure as 41270)
# ==============================================================================

# 12a. ICD-10 codes (arrays 0-74)
echo "=== BATCH 12a: HES ICD-10 codes arrays 0-74 ==="
python3 -c "
fields = ['participant.eid']
for a in range(75):
    fields.append(f'participant.p41270_a{a}')
print('\n'.join(fields))
" > /tmp/calon_icd10_a.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_icd10_a.txt \
  --output /tmp/calon_icd10_a.csv --delimiter "," 2>/dev/null

# 12b. ICD-10 codes (arrays 75-149)
echo "=== BATCH 12b: HES ICD-10 codes arrays 75-149 ==="
python3 -c "
fields = ['participant.eid']
for a in range(75, 150):
    fields.append(f'participant.p41270_a{a}')
print('\n'.join(fields))
" > /tmp/calon_icd10_b.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_icd10_b.txt \
  --output /tmp/calon_icd10_b.csv --delimiter "," 2>/dev/null

# 12c. ICD-10 dates (arrays 0-74)
echo "=== BATCH 12c: HES diagnosis dates arrays 0-74 ==="
python3 -c "
fields = ['participant.eid']
for a in range(75):
    fields.append(f'participant.p41280_a{a}')
print('\n'.join(fields))
" > /tmp/calon_dates_a.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_dates_a.txt \
  --output /tmp/calon_dates_a.csv --delimiter "," 2>/dev/null

# 12d. ICD-10 dates (arrays 75-149)
echo "=== BATCH 12d: HES diagnosis dates arrays 75-149 ==="
python3 -c "
fields = ['participant.eid']
for a in range(75, 150):
    fields.append(f'participant.p41280_a{a}')
print('\n'.join(fields))
" > /tmp/calon_dates_b.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_dates_b.txt \
  --output /tmp/calon_dates_b.csv --delimiter "," 2>/dev/null

# 12e. Merge ICD-10 and dates
echo "=== BATCH 12e: Merging HES ICD-10 and dates ==="
python3 -c "
import pandas as pd
import os

# ICD-10 codes
icd_files = ['/tmp/calon_icd10_a.csv', '/tmp/calon_icd10_b.csv']
icd_dfs = []
for f in icd_files:
    if os.path.exists(f):
        icd_dfs.append(pd.read_csv(f))
if icd_dfs:
    icd_merged = icd_dfs[0]
    for df in icd_dfs[1:]:
        icd_merged = icd_merged.merge(df, on='participant.eid', how='outer')
    icd_merged.to_csv('calon_batch12a_icd10_codes.csv', index=False)
    print(f'ICD-10 codes: {icd_merged.shape[0]} rows x {icd_merged.shape[1]} cols')
else:
    print('WARNING: No ICD-10 code files found')

# Dates
date_files = ['/tmp/calon_dates_a.csv', '/tmp/calon_dates_b.csv']
date_dfs = []
for f in date_files:
    if os.path.exists(f):
        date_dfs.append(pd.read_csv(f))
if date_dfs:
    dates_merged = date_dfs[0]
    for df in date_dfs[1:]:
        dates_merged = dates_merged.merge(df, on='participant.eid', how='outer')
    dates_merged.to_csv('calon_batch12b_icd10_dates.csv', index=False)
    print(f'ICD-10 dates: {dates_merged.shape[0]} rows x {dates_merged.shape[1]} cols')
else:
    print('WARNING: No ICD-10 date files found')
"
echo "Done: calon_batch12a_icd10_codes.csv, calon_batch12b_icd10_dates.csv"

# ==============================================================================
# BATCH 13: OPCS-4 PROCEDURES (Revascularisation)
# ==============================================================================
# Field 41272: Operative procedures - OPCS4 (up to 120+ arrays)
# Field 41282: Dates of operative procedures
#
# CRITICAL OPCS-4 codes for ASCVD:
#   K40-K46 — CABG
#   K49, K50, K75 — PCI / coronary angioplasty
#   L29-L35 — Carotid endarterectomy
#   L51-L63 — Peripheral arterial procedures
# ==============================================================================

echo "=== BATCH 13a: OPCS-4 procedure codes ==="
python3 -c "
fields = ['participant.eid']
for a in range(60):
    fields.append(f'participant.p41272_a{a}')
print('\n'.join(fields))
" > /tmp/calon_opcs_a.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_opcs_a.txt \
  --output /tmp/calon_opcs_a.csv --delimiter "," 2>/dev/null

echo "=== BATCH 13b: OPCS-4 procedure dates ==="
python3 -c "
fields = ['participant.eid']
for a in range(60):
    fields.append(f'participant.p41282_a{a}')
print('\n'.join(fields))
" > /tmp/calon_opcs_dates.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/calon_opcs_dates.txt \
  --output /tmp/calon_opcs_dates.csv --delimiter "," 2>/dev/null

# Merge OPCS
echo "=== BATCH 13c: Merging OPCS-4 ==="
python3 -c "
import pandas as pd, os
for pair in [('calon_batch13a_opcs_codes.csv', '/tmp/calon_opcs_a.csv'),
             ('calon_batch13b_opcs_dates.csv', '/tmp/calon_opcs_dates.csv')]:
    out, src = pair
    if os.path.exists(src):
        df = pd.read_csv(src)
        df.to_csv(out, index=False)
        print(f'{out}: {df.shape[0]} rows x {df.shape[1]} cols')
    else:
        print(f'WARNING: {src} not found')
"
echo "Done: calon_batch13a_opcs_codes.csv, calon_batch13b_opcs_dates.csv"

# ==============================================================================
# UPLOAD ALL OUTPUTS
# ==============================================================================
echo ""
echo "=== Uploading all CALON extraction files ==="
UPLOAD_FILES=(
  calon_batch1_demographics.csv
  calon_batch2_lipids.csv
  calon_batch3_biomarkers.csv
  calon_batch4_medications.csv
  calon_batch5_smoking_alcohol.csv
  calon_batch6_comorbidities.csv
  calon_batch7_ascvd_dates.csv
  calon_batch8_ascvd_selfreport.csv
  calon_batch9_cimt.csv
  calon_batch10_cac.csv
  calon_batch11_illness.csv
  calon_batch12a_icd10_codes.csv
  calon_batch12b_icd10_dates.csv
  calon_batch13a_opcs_codes.csv
  calon_batch13b_opcs_dates.csv
)

for f in "${UPLOAD_FILES[@]}"; do
  if [ -f "$f" ] && [ "$(wc -l < "$f")" -gt 1 ]; then
    dx upload "$f" --destination /calon_data/ --brief
    echo "  Uploaded: $f"
  else
    echo "  Skipped: $f (empty or not found)"
  fi
done

echo ""
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║                    ALL EXTRACTIONS COMPLETE                        ║"
echo "║                                                                    ║"
echo "║  Download these files from the Manage tab → /calon_data/           ║"
echo "║  Place in: C:/Users/nader/Downloads/calon_ukb_pipeline/data/       ║"
echo "║                                                                    ║"
echo "║  Then run: 01_CALON_build_cohort.R                                 ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"

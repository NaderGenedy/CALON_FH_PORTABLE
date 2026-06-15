"""
CALON Clinical Risk Model - Raw Data Inventory
================================================
Locates and verifies every raw UKB file needed for the clinical CALON-FH
risk model (lipid profile, ratios, prior ASCVD, FH, sex, smoking, T2DM,
HTN, imaging -> incident ASCVD).

Critically verifies the ASCVD outcome source: ukb_dates_mace.csv is
MISNAMED (contains I10-I15 hypertension first-occurrence, not ASCVD).
Correct ASCVD must come from HES ICD-10 array.
"""
import os
import pandas as pd

# Candidate raw files per feature domain
CANDIDATES = {
    'lipids_assay':   [r'D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_longitudinal_lipids.csv',
                       r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'],
    'apob_lpa':       [r'D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_apob_lpa.csv'],
    'icd10_hes':      [r'D:/Projects/CALON_AlphaFold_Rebuild/ukb_icd10_full.csv',
                       r'D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_icd10_codes.csv'],
    'fh_carriers':    [r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv'],
    'sex':            [r'D:/Projects/CALON_AlphaFold_Rebuild/Paper3_ASCVD_Prediction/calon_sex.csv'],
    'smoking_bp':     [r'D:/Projects/CALON_AlphaFold_Rebuild/Paper3_ASCVD_Prediction/paper3_smoking_bp.csv',
                       r'D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_smoking.csv'],
    'hba1c_t2dm':     [r'D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_hba1c.csv'],
    'medications':    [r'D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_medications.csv'],
    'carotid_imt':    [r'D:/Projects/CALON_AlphaFold_Rebuild/ukb_carotid_imt.csv',
                       r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carotid_imt.csv'],
    'cardiac_mri':    [r'D:/Projects/CALON_AlphaFold_Rebuild/ukb_cardiac_mri.csv',
                       r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_cardiac_mri.csv'],
    'recruitment':    [r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_recruitment_dates.csv',
                       r'D:/Projects/CALON_AlphaFold_Rebuild/ukb_recruitment_dates.csv'],
    'mace_MISNAMED':  [r'D:/Projects/CALON_AlphaFold_Rebuild/ukb_dates_mace.csv'],
}

print('='*72)
print('CALON CLINICAL RISK MODEL - RAW DATA INVENTORY')
print('='*72)

found = {}
for domain, paths in CANDIDATES.items():
    print(f'\n[{domain}]')
    hit = None
    for p in paths:
        if os.path.exists(p):
            sz = os.path.getsize(p) / 1e6
            try:
                head = pd.read_csv(p, nrows=2)
                ncol = len(head.columns)
                cols = [c.replace('participant.', '') for c in head.columns]
                print(f'  FOUND  {p}')
                print(f'         {sz:.1f} MB, {ncol} columns')
                print(f'         cols: {cols[:18]}{"..." if ncol>18 else ""}')
                if hit is None:
                    hit = p
            except Exception as e:
                print(f'  ERROR  {p}: {e}')
        else:
            print(f'  absent {p}')
    found[domain] = hit

print('\n' + '='*72)
print('RESOLVED SOURCES')
print('='*72)
for d, p in found.items():
    status = 'OK' if p else 'MISSING'
    print(f'  {d:18s} {status:8s} {p or ""}')

# Verify the MACE misnaming explicitly
print('\n' + '='*72)
print('ASCVD OUTCOME SOURCE VERIFICATION')
print('='*72)
mace = found.get('mace_MISNAMED')
if mace:
    df = pd.read_csv(mace, nrows=5000)
    foc_cols = [c for c in df.columns if 'p131' in c]
    print(f'  ukb_dates_mace.csv first-occurrence fields: {[c.replace("participant.","") for c in foc_cols]}')
    print(f'  -> These are p131286-p131296 = ICD-10 I10-I15 HYPERTENSION family.')
    print(f'  -> DO NOT use as ASCVD outcome. Verdict: MISNAMED FILE.')
    for c in foc_cols:
        nn = df[c].notna().sum()
        print(f'     {c.replace("participant.",""):14s} populated in {nn}/5000 sampled rows ({nn/50:.1f}%)')
icd = found.get('icd10_hes')
if icd:
    print(f'\n  CORRECT ASCVD source: {icd}')
    print(f'  -> Define ASCVD from HES ICD-10 array: I20-I25, I63, I64, I70, I73, I74')
    print(f'  -> EXCLUDE I10-I15 (hypertension), I35 (aortic stenosis)')

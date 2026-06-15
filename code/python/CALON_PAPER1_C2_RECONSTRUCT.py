"""
CALON Paper 1 — Reconstruct the expanded 775 statin cohort (claim C2)
======================================================================
The carrier CSV's on_statin = primary-field 374. Manuscript expands to 775
using self-reported cholesterol medication across p6153/p6177 instances 0-2.

UKB coding:
  p6153 (female) / p6177 (male) — medication for cholesterol, BP, diabetes
  Value 1 = "Cholesterol lowering medication"
  Arrays may hold value 1 in any element.

Method:
  1. Load 3,544 carriers
  2. Load medication arrays for those eids
  3. on_statin_expanded = TRUE if value 1 appears in ANY p6153/p6177 cell
  4. Re-derive SSS<->residual LDL on the expanded cohort
"""
import pandas as pd
import numpy as np
from scipy import stats

CARRIERS = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_full_prs_sss.csv'
MEDS = r'D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_medications.csv'

print('CALON Paper 1 C2 — reconstruct expanded statin cohort')
print('='*64)

car = pd.read_csv(CARRIERS)
car = car.dropna(subset=['eid'])
car['eid'] = car['eid'].astype(int)
print(f'Carriers loaded: {len(car)}')

meds = pd.read_csv(MEDS)
meds = meds.rename(columns={'participant.eid':'eid'})
med_cols = [c for c in meds.columns if 'p6153' in c or 'p6177' in c]
meds = meds.dropna(subset=['eid'])
meds['eid'] = meds['eid'].astype(int)

# Cholesterol-lowering medication = value 1 in any p6153/p6177 cell.
# UKB arrays are pipe- or comma-delimited strings; value 1 = cholesterol med.
def has_chol_med(row):
    # UKB arrays formatted as '[1,2]' '[1]' '[-7]'. Cholesterol med = value 1.
    for c in med_cols:
        v = row[c]
        if pd.isna(v):
            continue
        tokens = str(v).strip('[]').replace('|', ',').replace(';', ',').split(',')
        if '1' in [t.strip() for t in tokens]:
            return True
    return False

meds_carrier = meds[meds['eid'].isin(set(car['eid']))].copy()
print(f'Medication rows matched to carriers: {len(meds_carrier)}')
meds_carrier['on_statin_expanded'] = meds_carrier.apply(has_chol_med, axis=1).astype(int)
n_expanded = meds_carrier['on_statin_expanded'].sum()
print(f'Expanded statin-treated carriers: {n_expanded}')

# Merge back
car = car.merge(meds_carrier[['eid','on_statin_expanded']], on='eid', how='left')
car['on_statin_expanded'] = car['on_statin_expanded'].fillna(0).astype(int)

# Re-derive C2 on expanded cohort
treated = car[car['on_statin_expanded']==1].dropna(subset=['sss_v3','ldl'])
print()
print(f'Expanded cohort complete-case (sss_v3 & ldl): n={len(treated)}')
if len(treated) > 50:
    rho, p = stats.spearmanr(treated['sss_v3'], treated['ldl'])
    print(f'  SSS v3 <-> residual LDL:  rho={rho:.4f}  p={p:.4f}')
    print(f'  Manuscript claim:         rho=0.083  p=0.020  n=775')
    # PRS-adjusted
    adj = treated.dropna(subset=['ldl_prs'])
    if len(adj) > 50:
        import statsmodels.api as sm
        X = sm.add_constant(adj['ldl_prs'])
        resid_ldl = sm.OLS(adj['ldl'], X).fit().resid
        rho_adj, p_adj = stats.spearmanr(adj['sss_v3'], resid_ldl)
        print(f'  PRS-adjusted (residualised LDL): rho={rho_adj:.4f}  p={p_adj:.4f}  n={len(adj)}')

# Also check: position-level SSS may be the f_* layers; how many have ANY sss
treated_any = car[car['on_statin_expanded']==1]
n_sss = treated_any['sss_v3'].notna().sum()
n_ldl = treated_any['ldl'].notna().sum()
print()
print(f'Expanded cohort total (incl. missing SSS): n={len(treated_any)}')
print(f'  with sss_v3 non-null:  {n_sss}')
print(f'  with ldl non-null:     {n_ldl}')

#!/usr/bin/env python3
"""
72_detailed_thresholds_correlations.py
========================================
1. Lp(a) threshold for EACH vascular outcome separately
2. Full ApoB / LDL / Lp(a) correlation matrix
3. ApoB/LDL ratio vs Lp(a) detailed analysis
4. Outcome-specific Lp(a) quintile gradients

Author: Dr Nader Genedy
Date:   April 2026
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings('ignore')
import os

ANALYSIS = r"C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis"
DATA = r"D:/CALON_AF3_PROJECT/data"

print("=" * 70)
print("DETAILED THRESHOLDS + FULL CORRELATION MATRIX")
print("=" * 70)

# Load
lpa = pd.read_csv(f"{DATA}/ukb_wide/calon_lpa_CORRECT.csv")
lpa.columns = [c.replace('participant.', '') for c in lpa.columns]
lpa['eid'] = lpa['eid'].astype(str)
lpa['lpa'] = pd.to_numeric(lpa['p30790_i0'], errors='coerce')

apob = pd.read_csv(f"{DATA}/ukb_wide/calon_apob_CORRECT.csv")
apob.columns = [c.replace('participant.', '') for c in apob.columns]
apob['eid'] = apob['eid'].astype(str)
apob['apob'] = pd.to_numeric(apob['p30640_i0'], errors='coerce')

lipids = pd.read_csv(f"{DATA}/ukb_wide/ukb_reviewer_longitudinal_lipids.csv")
lipids.columns = [c.replace('participant.', '') for c in lipids.columns]
lipids['eid'] = lipids['eid'].astype(str)
lipids['ldl_direct'] = pd.to_numeric(lipids['p30780_i0'], errors='coerce')
for col in lipids.columns:
    if '30690' in col and '_i0' in col:
        lipids['tc'] = pd.to_numeric(lipids[col], errors='coerce')
    if '30760' in col and '_i0' in col:
        lipids['hdl'] = pd.to_numeric(lipids[col], errors='coerce')
    if '30870' in col and '_i0' in col:
        lipids['trig'] = pd.to_numeric(lipids[col], errors='coerce')
lipids['ldl_calc'] = lipids['tc'] - lipids['hdl'] - (lipids['trig'] / 2.2)

icd = pd.read_csv(f"{DATA}/ukb_wide/calon_batch_icd10.csv",
                   usecols=['participant.eid', 'participant.p41270'], dtype=str)
icd.columns = ['eid', 'icd10_raw']
icd['ascvd'] = icd['icd10_raw'].str.contains('"I2[0-5]', na=False, regex=True).astype(int)
icd['mi'] = icd['icd10_raw'].str.contains('"I21|"I22', na=False, regex=True).astype(int)
icd['stroke'] = icd['icd10_raw'].str.contains('"I6[3-4]', na=False, regex=True).astype(int)
icd['tia'] = icd['icd10_raw'].str.contains('"G45', na=False, regex=True).astype(int)
icd['pvd'] = icd['icd10_raw'].str.contains('"I7[0-4]', na=False, regex=True).astype(int)
icd['ao_stenosis'] = icd['icd10_raw'].str.contains('"I35', na=False, regex=True).astype(int)
icd['carotid'] = icd['icd10_raw'].str.contains('"I65', na=False, regex=True).astype(int)

carriers = pd.read_csv(os.path.join(ANALYSIS, "ukb_carriers_annotated.csv"), dtype={'eid': str})
patients = carriers.sort_values('sss', ascending=False).drop_duplicates('eid')
fh_eids = set(patients['eid'])

# Merge
ukb = lpa[['eid', 'lpa']].merge(apob[['eid', 'apob']], on='eid', how='left')
ukb = ukb.merge(lipids[['eid', 'ldl_direct', 'ldl_calc', 'tc', 'hdl', 'trig']], on='eid', how='left')
ukb = ukb.merge(icd[['eid', 'ascvd', 'mi', 'stroke', 'tia', 'pvd', 'ao_stenosis', 'carotid']], on='eid', how='left')
ukb['is_fh'] = ukb['eid'].isin(fh_eids).astype(int)
ukb = ukb[ukb['lpa'].notna()].copy()

ukb['apob_ldl_direct'] = np.where(ukb['ldl_direct'] > 0, ukb['apob'] / ukb['ldl_direct'], np.nan)
ukb['apob_ldl_calc'] = np.where(ukb['ldl_calc'] > 0, ukb['apob'] / ukb['ldl_calc'], np.nan)

fh = ukb[ukb['is_fh'] == 1].copy()
nonfh = ukb[ukb['is_fh'] == 0].copy()

# Wales
ar = pd.read_csv(f"{DATA}/fh_cohort/calon_ukb_analysis_ready.csv")
ar['lpa_val'] = pd.to_numeric(ar['lpa'], errors='coerce')
ar['ascvd'] = pd.to_numeric(ar['ascvd_combined'], errors='coerce')

# ============================================================================
# PART 1: Lp(a) THRESHOLD FOR EACH OUTCOME — NON-FH
# ============================================================================
print(f"\n{'='*70}")
print("PART 1: Lp(a) THRESHOLD PER OUTCOME (Non-FH, n={})".format(len(nonfh)))
print("=" * 70)

outcomes = [('mi', 'MI'), ('ao_stenosis', 'Aortic Stenosis'), ('stroke', 'Stroke'),
            ('pvd', 'PVD'), ('carotid', 'Carotid'), ('tia', 'TIA'), ('ascvd', 'ASCVD')]

for col, name in outcomes:
    events = nonfh[col].sum()
    if events < 100:
        continue

    print(f"\n  --- {name} (n events={events:.0f}) ---")
    print(f"  {'Threshold':>12s}  {'Rate above':>10s}  {'Rate below':>10s}  {'RR':>6s}  {'P':>12s}")

    best_rr = 0
    best_thresh = 0
    for t in [25, 30, 40, 50, 75, 100, 125, 150]:
        hi = nonfh[nonfh['lpa'] > t]
        lo = nonfh[nonfh['lpa'] <= t]
        if len(hi) >= 100 and len(lo) >= 100:
            r_hi = 100 * hi[col].mean()
            r_lo = 100 * lo[col].mean()
            rr = r_hi / r_lo if r_lo > 0 else 0
            table = np.array([[int(hi[col].sum()), len(hi) - int(hi[col].sum())],
                              [int(lo[col].sum()), len(lo) - int(lo[col].sum())]])
            chi2, p, _, _ = stats.chi2_contingency(table)
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            print(f"  >{t:>3d} nmol/L  {r_hi:>9.2f}%  {r_lo:>9.2f}%  {rr:>6.2f}  {p:>12.2e} {sig}")

            if rr > best_rr:
                best_rr = rr
                best_thresh = t

    print(f"  Best threshold: >{best_thresh} nmol/L (RR={best_rr:.2f})")

# ============================================================================
# PART 1B: SAME FOR FH (Wales — better powered)
# ============================================================================
print(f"\n{'='*70}")
print("PART 1B: Lp(a) THRESHOLD PER OUTCOME (Wales FH)")
print("=" * 70)

ar_lpa = ar[ar['lpa_val'].notna()].copy()
print(f"  n={len(ar_lpa)}, ASCVD events={ar_lpa['ascvd'].sum():.0f}")

for t in [30, 50, 75, 100, 125]:
    hi = ar_lpa[ar_lpa['lpa_val'] > t]
    lo = ar_lpa[ar_lpa['lpa_val'] <= t]
    if len(hi) >= 10 and len(lo) >= 10:
        r_hi = 100 * hi['ascvd'].mean()
        r_lo = 100 * lo['ascvd'].mean()
        rr = r_hi / r_lo if r_lo > 0 else 0
        table = np.array([[int(hi['ascvd'].sum()), len(hi) - int(hi['ascvd'].sum())],
                          [int(lo['ascvd'].sum()), len(lo) - int(lo['ascvd'].sum())]])
        _, p, _, _ = stats.chi2_contingency(table)
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        print(f"  >{t:>3d}: ASCVD {r_hi:.1f}% vs {r_lo:.1f}%, RR={rr:.2f}, P={p:.4f} {sig}")

# ============================================================================
# PART 2: FULL CORRELATION MATRIX — ApoB, LDL, Lp(a), ratio
# ============================================================================
print(f"\n{'='*70}")
print("PART 2: FULL CORRELATION MATRIX")
print("=" * 70)

corr_vars = [
    ('lpa', 'Lp(a)'),
    ('apob', 'ApoB'),
    ('ldl_direct', 'LDL direct'),
    ('ldl_calc', 'LDL Friedewald'),
    ('tc', 'Total cholesterol'),
    ('hdl', 'HDL-C'),
    ('trig', 'Triglycerides'),
    ('apob_ldl_direct', 'ApoB/LDL ratio (direct)'),
    ('apob_ldl_calc', 'ApoB/LDL ratio (calc)'),
]

for group_label, df_g in [("FH (UKB)", fh), ("Non-FH (UKB)", nonfh)]:
    print(f"\n  --- {group_label} ---")
    print(f"  {'':>25s}", end='')
    for _, name in corr_vars[:6]:
        print(f"  {name[:8]:>8s}", end='')
    print()

    for col1, name1 in corr_vars:
        print(f"  {name1:>25s}", end='')
        for col2, name2 in corr_vars[:6]:
            v = df_g[[col1, col2]].dropna()
            if len(v) >= 50 and col1 != col2:
                rho, p = stats.spearmanr(v[col1], v[col2])
                sig = '*' if p < 0.05 else ''
                print(f"  {rho:>7.3f}{sig}", end='')
            elif col1 == col2:
                print(f"  {'1.000':>8s}", end='')
            else:
                print(f"  {'N/A':>8s}", end='')
        print()

# ============================================================================
# PART 3: DETAILED ApoB/LDL vs Lp(a)
# ============================================================================
print(f"\n{'='*70}")
print("PART 3: ApoB/LDL DISCORDANCE vs Lp(a) — DETAILED")
print("=" * 70)

for group_label, df_g in [("FH", fh), ("Non-FH", nonfh)]:
    sub = df_g[['lpa', 'apob', 'ldl_direct', 'ldl_calc', 'apob_ldl_direct', 'apob_ldl_calc']].dropna(subset=['lpa', 'apob'])
    if len(sub) < 50:
        continue

    print(f"\n  --- {group_label} (n={len(sub)}) ---")

    # Lp(a) quintiles vs each lipid measure
    sub['lpa_q'] = pd.qcut(sub['lpa'], 5, labels=['Q1', 'Q2', 'Q3', 'Q4', 'Q5'])

    for metric, name in [('apob', 'ApoB'), ('ldl_direct', 'LDL direct'),
                          ('ldl_calc', 'LDL calc'), ('apob_ldl_direct', 'ApoB/LDL ratio')]:
        vals = []
        for q in ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']:
            v = sub[sub['lpa_q'] == q][metric].dropna()
            vals.append(v.median())

        # Trend
        rho, p = stats.spearmanr(range(5), vals)
        gradient = f"{vals[0]:.3f} -> {vals[-1]:.3f}"
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        print(f"    {name:>15s} by Lp(a) Q: {gradient}, trend rho={rho:.3f}, P={p:.4f} {sig}")

    # Does Lp(a) correction change LDL?
    sub_lpa = sub[sub['lpa'].notna()].copy()
    sub_lpa['lpa_chol'] = (sub_lpa['lpa'] / 2.5) * 0.30 / 38.7
    sub_lpa['ldl_corr_direct'] = sub_lpa['ldl_direct'] - sub_lpa['lpa_chol']
    sub_lpa['ldl_corr_calc'] = sub_lpa['ldl_calc'] - sub_lpa['lpa_chol']

    print(f"\n    Lp(a) cholesterol correction ({group_label}):")
    print(f"      Mean Lp(a)-C: {sub_lpa['lpa_chol'].mean():.3f} mmol/L")
    print(f"      Direct LDL: {sub_lpa['ldl_direct'].median():.2f} -> corrected: {sub_lpa['ldl_corr_direct'].median():.2f}")
    print(f"      Calc LDL:   {sub_lpa['ldl_calc'].median():.2f} -> corrected: {sub_lpa['ldl_corr_calc'].median():.2f}")

    # Does correction change ApoB/LDL ratio?
    sub_lpa['ratio_corr'] = np.where(sub_lpa['ldl_corr_direct'] > 0,
                                      sub_lpa['apob'] / sub_lpa['ldl_corr_direct'], np.nan)
    valid_ratio = sub_lpa[sub_lpa['ratio_corr'].notna() & sub_lpa['apob_ldl_direct'].notna()]
    if len(valid_ratio) >= 50:
        print(f"      ApoB/LDL ratio uncorrected: {valid_ratio['apob_ldl_direct'].median():.3f}")
        print(f"      ApoB/LDL ratio Lp(a)-corrected: {valid_ratio['ratio_corr'].median():.3f}")
        diff = valid_ratio['ratio_corr'].median() - valid_ratio['apob_ldl_direct'].median()
        print(f"      Change: {diff:+.3f} ({'increases' if diff > 0 else 'decreases'} discordance)")

# ============================================================================
# PART 4: Lp(a) QUINTILE x EACH OUTCOME — WITH ABSOLUTE NUMBERS
# ============================================================================
print(f"\n{'='*70}")
print("PART 4: Lp(a) QUINTILE GRADIENTS — DETAILED (Non-FH)")
print("=" * 70)

nonfh_q = nonfh.copy()
nonfh_q['lpa_q'] = pd.qcut(nonfh_q['lpa'], 5, labels=['Q1', 'Q2', 'Q3', 'Q4', 'Q5'])

for col, name in outcomes:
    events = nonfh[col].sum()
    if events < 100:
        continue

    print(f"\n  --- {name} ---")
    print(f"  {'Quintile':>10s}  {'Lp(a) range':>15s}  {'N':>8s}  {'Events':>7s}  {'Rate%':>7s}  {'RR vs Q1':>9s}")

    q1_rate = None
    for q in ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']:
        sub = nonfh_q[nonfh_q['lpa_q'] == q]
        n_ev = sub[col].sum()
        rate = 100 * sub[col].mean()
        lo = sub['lpa'].min()
        hi = sub['lpa'].max()
        if q1_rate is None:
            q1_rate = rate
        rr = rate / q1_rate if q1_rate > 0 else 0
        print(f"  {q:>10s}  {lo:>6.0f}-{hi:>6.0f}  {len(sub):>8d}  {n_ev:>7.0f}  {rate:>7.2f}  {rr:>9.2f}")

    # Q5 vs Q1 formal test
    q1 = nonfh_q[nonfh_q['lpa_q'] == 'Q1']
    q5 = nonfh_q[nonfh_q['lpa_q'] == 'Q5']
    table = np.array([[int(q5[col].sum()), len(q5) - int(q5[col].sum())],
                      [int(q1[col].sum()), len(q1) - int(q1[col].sum())]])
    chi2, p, _, _ = stats.chi2_contingency(table)
    rr_q5q1 = (q5[col].mean()) / (q1[col].mean()) if q1[col].mean() > 0 else 0
    print(f"  Q5 vs Q1: RR={rr_q5q1:.2f}, P={p:.2e}")

# ============================================================================
# PART 5: WALES FH — FULL LIPID CORRELATIONS WITH Lp(a)
# ============================================================================
print(f"\n{'='*70}")
print("PART 5: WALES FH — FULL CORRELATIONS WITH Lp(a)")
print("=" * 70)

wales_vars = ['ldl', 'tc', 'hdl', 'trig', 'apob', 'apob_ldl_ratio', 'bmi', 'crp']
for col in wales_vars:
    if col in ar.columns:
        v = ar[['lpa_val', col]].dropna()
        v[col] = pd.to_numeric(v[col], errors='coerce')
        v = v.dropna()
        if len(v) >= 30:
            rho, p = stats.spearmanr(v['lpa_val'], v[col])
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            print(f"  Lp(a) vs {col:>15s}: rho={rho:.3f}, P={p:.4f} {sig} (n={len(v)})")

print("\n" + "=" * 70)
print("COMPLETE")
print("=" * 70)

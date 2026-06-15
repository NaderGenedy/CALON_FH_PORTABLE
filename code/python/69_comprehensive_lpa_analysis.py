#!/usr/bin/env python3
"""
69_comprehensive_lpa_analysis.py
==================================
Comprehensive Lp(a) analysis covering ALL remaining questions:

PART 1: Lp(a) threshold for ASCVD risk — FH vs matched non-FH
PART 2: Risk factor correlations (smoking, T2DM, BMI, TRG, HTN) × Lp(a)
PART 3: ApoB/LDL discordance vs Lp(a)
PART 4: Vascular outcomes (stroke, TIA, PVD, CAC) by Lp(a)
PART 5: LPA gene mutations in UKB

Author: Dr Nader Genedy
Date:   April 2026
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
import warnings
warnings.filterwarnings('ignore')
import os

ANALYSIS = r"C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis"
DATA = r"D:/CALON_AF3_PROJECT/data"

print("=" * 70)
print("COMPREHENSIVE Lp(a) ANALYSIS")
print("=" * 70)

# ============================================================================
# LOAD ALL DATA
# ============================================================================
print("\n[LOADING DATA]...")

# Lp(a)
lpa = pd.read_csv(f"{DATA}/ukb_wide/calon_lpa_CORRECT.csv")
lpa.columns = [c.replace('participant.', '') for c in lpa.columns]
lpa['eid'] = lpa['eid'].astype(str)
lpa['lpa'] = pd.to_numeric(lpa['p30790_i0'], errors='coerce')

# ApoB
apob = pd.read_csv(f"{DATA}/ukb_wide/calon_apob_CORRECT.csv")
apob.columns = [c.replace('participant.', '') for c in apob.columns]
apob['eid'] = apob['eid'].astype(str)
apob['apob'] = pd.to_numeric(apob['p30640_i0'], errors='coerce')

# Sex
sex = pd.read_csv(f"{DATA}/ukb_wide/calon_sex.csv")
sex.columns = [c.replace('participant.', '') for c in sex.columns]
sex.rename(columns={'p31': 'sex'}, inplace=True)
sex['eid'] = sex['eid'].astype(str)

# Demographics
demo = pd.read_csv(f"{DATA}/ukb_wide/ukb_reviewer_demographics.csv")
demo.columns = [c.replace('participant.', '') for c in demo.columns]
demo['eid'] = demo['eid'].astype(str)
demo['age'] = 2010 - pd.to_numeric(demo['p34'], errors='coerce')
demo['bmi'] = pd.to_numeric(demo['p21001_i0'], errors='coerce')

# Lipids
lipids = pd.read_csv(f"{DATA}/ukb_wide/ukb_reviewer_longitudinal_lipids.csv")
lipids.columns = [c.replace('participant.', '') for c in lipids.columns]
lipids['eid'] = lipids['eid'].astype(str)
lipids['ldl_direct'] = pd.to_numeric(lipids['p30780_i0'], errors='coerce')
# Calculated LDL columns
for col in lipids.columns:
    if '30690' in col:  # Total cholesterol
        lipids['tc'] = pd.to_numeric(lipids[col], errors='coerce')
    if '30760' in col:  # HDL
        lipids['hdl'] = pd.to_numeric(lipids[col], errors='coerce')
    if '30870' in col:  # Triglycerides
        lipids['trig'] = pd.to_numeric(lipids[col], errors='coerce')

# Friedewald calculated LDL
lipids['ldl_calc'] = lipids['tc'] - lipids['hdl'] - (lipids['trig'] / 2.2)

# Smoking
smoke = pd.read_csv(f"{DATA}/ukb_wide/ukb_reviewer_smoking.csv")
smoke.columns = [c.replace('participant.', '') for c in smoke.columns]
smoke['eid'] = smoke['eid'].astype(str)
smoke['smoking_status'] = pd.to_numeric(smoke['p20116_i0'], errors='coerce')
smoke['current_smoker'] = (smoke['smoking_status'] == 2).astype(int)
smoke['ever_smoker'] = (smoke['smoking_status'].isin([1, 2])).astype(int)

# ICD-10
icd = pd.read_csv(f"{DATA}/ukb_wide/calon_batch_icd10.csv",
                   usecols=['participant.eid', 'participant.p41270'], dtype=str)
icd.columns = ['eid', 'icd10_raw']

# Parse multiple outcomes
icd['ascvd'] = icd['icd10_raw'].str.contains('"I2[0-5]', na=False, regex=True).astype(int)
icd['mi'] = icd['icd10_raw'].str.contains('"I21|"I22', na=False, regex=True).astype(int)
icd['stroke'] = icd['icd10_raw'].str.contains('"I6[3-4]', na=False, regex=True).astype(int)
icd['tia'] = icd['icd10_raw'].str.contains('"G45', na=False, regex=True).astype(int)
icd['pvd'] = icd['icd10_raw'].str.contains('"I7[0-4]', na=False, regex=True).astype(int)
icd['ao_stenosis'] = icd['icd10_raw'].str.contains('"I35', na=False, regex=True).astype(int)
icd['htn'] = icd['icd10_raw'].str.contains('"I10|"I11|"I12|"I13|"I15', na=False, regex=True).astype(int)
icd['t2dm'] = icd['icd10_raw'].str.contains('"E11', na=False, regex=True).astype(int)
icd['carotid'] = icd['icd10_raw'].str.contains('"I65', na=False, regex=True).astype(int)

# FH carriers
carriers = pd.read_csv(os.path.join(ANALYSIS, "ukb_carriers_annotated.csv"), dtype={'eid': str})
patients = carriers.sort_values('sss', ascending=False).drop_duplicates('eid')
fh_eids = set(patients['eid'])

# Merge everything for ALL UKB
ukb = lpa[['eid', 'lpa']].copy()
for d, cols in [(apob, ['eid', 'apob']), (sex, ['eid', 'sex']),
                (demo, ['eid', 'age', 'bmi']),
                (lipids, ['eid', 'ldl_direct', 'ldl_calc', 'tc', 'hdl', 'trig']),
                (smoke, ['eid', 'current_smoker', 'ever_smoker']),
                (icd, ['eid', 'ascvd', 'mi', 'stroke', 'tia', 'pvd', 'ao_stenosis',
                       'htn', 't2dm', 'carotid'])]:
    ukb = ukb.merge(d[cols], on='eid', how='left')

ukb['is_fh'] = ukb['eid'].isin(fh_eids).astype(int)
ukb = ukb[ukb['lpa'].notna()].copy()

# ApoB/LDL ratios
ukb['apob_ldl_direct'] = np.where(ukb['ldl_direct'] > 0, ukb['apob'] / ukb['ldl_direct'], np.nan)
ukb['apob_ldl_calc'] = np.where(ukb['ldl_calc'] > 0, ukb['apob'] / ukb['ldl_calc'], np.nan)
ukb['discordant'] = (ukb['apob_ldl_direct'] >= 0.31).astype(int)

fh = ukb[ukb['is_fh'] == 1].copy()
nonfh = ukb[ukb['is_fh'] == 0].copy()

print(f"  FH with Lp(a): {len(fh)}")
print(f"  Non-FH with Lp(a): {len(nonfh)}")

# Also load Wales FH
ar = pd.read_csv(f"{DATA}/fh_cohort/calon_ukb_analysis_ready.csv")
ar['lpa_val'] = pd.to_numeric(ar['lpa'], errors='coerce')
ar['ascvd'] = pd.to_numeric(ar['ascvd_combined'], errors='coerce')
ar_lpa = ar[ar['lpa_val'].notna()].copy()
print(f"  Wales FH with Lp(a): {len(ar_lpa)}")

# ============================================================================
# PART 1: Lp(a) THRESHOLD FOR ASCVD RISK
# ============================================================================
print(f"\n{'='*70}")
print("PART 1: AT WHAT Lp(a) DOES ASCVD RISK BECOME ELEVATED?")
print("=" * 70)

for label, df_a in [("FH (UKB)", fh), ("Non-FH (UKB)", nonfh), ("FH (Wales)", ar_lpa)]:
    lpa_col = 'lpa' if label != "FH (Wales)" else 'lpa_val'
    ascvd_col = 'ascvd'

    sub = df_a[[lpa_col, ascvd_col]].dropna()
    if len(sub) < 50 or sub[ascvd_col].sum() < 10:
        continue

    print(f"\n  --- {label} (n={len(sub)}, events={sub[ascvd_col].sum():.0f}) ---")

    # Decile analysis
    sub['decile'] = pd.qcut(sub[lpa_col], 10, labels=False, duplicates='drop')
    print(f"  {'Decile':>8s}  {'Lp(a) range':>15s}  {'N':>6s}  {'Events':>7s}  {'Rate%':>7s}")

    baseline_rate = None
    threshold_found = False
    for d in sorted(sub['decile'].unique()):
        grp = sub[sub['decile'] == d]
        rate = 100 * grp[ascvd_col].mean()
        lo = grp[lpa_col].min()
        hi = grp[lpa_col].max()
        if baseline_rate is None:
            baseline_rate = rate

        elevated = ""
        if rate > baseline_rate * 1.5 and not threshold_found:
            elevated = " <-- 1.5x baseline"
            threshold_found = True

        print(f"  {d:>8d}  {lo:>6.0f}-{hi:>6.0f}  {len(grp):>6d}  {grp[ascvd_col].sum():>7.0f}  {rate:>7.1f}{elevated}")

    # Continuous threshold analysis
    print(f"\n  Threshold analysis:")
    for t in [25, 30, 50, 75, 100, 125, 150]:
        hi = sub[sub[lpa_col] > t]
        lo = sub[sub[lpa_col] <= t]
        if len(hi) >= 20 and len(lo) >= 20:
            r_hi = 100 * hi[ascvd_col].mean()
            r_lo = 100 * lo[ascvd_col].mean()
            rr = r_hi / r_lo if r_lo > 0 else 0
            table = np.array([[int(hi[ascvd_col].sum()), len(hi) - int(hi[ascvd_col].sum())],
                              [int(lo[ascvd_col].sum()), len(lo) - int(lo[ascvd_col].sum())]])
            chi2, p, _, _ = stats.chi2_contingency(table)
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            print(f"    >{t:>3d}: {r_hi:.1f}% vs {r_lo:.1f}%, RR={rr:.2f}, P={p:.4f} {sig}")

# ============================================================================
# PART 2: RISK FACTOR CORRELATIONS WITH Lp(a)
# ============================================================================
print(f"\n{'='*70}")
print("PART 2: RISK FACTORS × Lp(a)")
print("=" * 70)

risk_factors = [
    ('current_smoker', 'Current smoker', 'binary'),
    ('ever_smoker', 'Ever smoker', 'binary'),
    ('t2dm', 'T2DM (ICD E11)', 'binary'),
    ('htn', 'Hypertension (ICD I10-I15)', 'binary'),
    ('bmi', 'BMI', 'continuous'),
    ('trig', 'Triglycerides', 'continuous'),
    ('hdl', 'HDL-C', 'continuous'),
]

for group_label, df_g in [("FH", fh), ("Non-FH", nonfh)]:
    print(f"\n  --- {group_label} ---")
    for col, name, rtype in risk_factors:
        if col not in df_g.columns:
            continue
        valid = df_g[['lpa', col]].dropna()
        if len(valid) < 50:
            continue

        if rtype == 'binary':
            yes = valid[valid[col] == 1]['lpa']
            no = valid[valid[col] == 0]['lpa']
            if len(yes) >= 10 and len(no) >= 10:
                u, p = stats.mannwhitneyu(yes, no)
                d_cohen = (yes.mean() - no.mean()) / np.sqrt((yes.std()**2 + no.std()**2) / 2)
                print(f"    {name:>25s}: {col}=1 Lp(a)={yes.median():.1f} vs {col}=0 Lp(a)={no.median():.1f}, "
                      f"P={p:.4f}, d={d_cohen:.3f}")
        else:
            rho, p = stats.spearmanr(valid['lpa'], valid[col])
            print(f"    {name:>25s}: rho={rho:.3f}, P={p:.4f}")

# ============================================================================
# PART 3: ApoB/LDL DISCORDANCE vs Lp(a)
# ============================================================================
print(f"\n{'='*70}")
print("PART 3: ApoB/LDL DISCORDANCE × Lp(a)")
print("=" * 70)

for group_label, df_g in [("FH", fh), ("Non-FH", nonfh)]:
    sub = df_g[['lpa', 'apob', 'ldl_direct', 'ldl_calc', 'apob_ldl_direct', 'apob_ldl_calc', 'discordant']].dropna(subset=['lpa', 'apob_ldl_direct'])
    if len(sub) < 50:
        continue

    print(f"\n  --- {group_label} (n={len(sub)}) ---")

    # Correlations
    for x, y, name in [('lpa', 'apob', 'Lp(a) vs ApoB'),
                        ('lpa', 'ldl_direct', 'Lp(a) vs LDL direct'),
                        ('lpa', 'ldl_calc', 'Lp(a) vs LDL calculated'),
                        ('lpa', 'apob_ldl_direct', 'Lp(a) vs ApoB/LDL ratio (direct)'),
                        ('lpa', 'apob_ldl_calc', 'Lp(a) vs ApoB/LDL ratio (calc)')]:
        v = sub[[x, y]].dropna()
        if len(v) >= 50:
            rho, p = stats.spearmanr(v[x], v[y])
            print(f"    {name:>35s}: rho={rho:.3f}, P={p:.2e}")

    # Discordant vs concordant Lp(a)
    disc = sub[sub['discordant'] == 1]['lpa']
    conc = sub[sub['discordant'] == 0]['lpa']
    if len(disc) >= 20 and len(conc) >= 20:
        u, p = stats.mannwhitneyu(disc, conc)
        print(f"    {'Discordant vs concordant Lp(a)':>35s}: {disc.median():.1f} vs {conc.median():.1f}, P={p:.4f}")

# ============================================================================
# PART 4: VASCULAR OUTCOMES BY Lp(a)
# ============================================================================
print(f"\n{'='*70}")
print("PART 4: VASCULAR OUTCOMES BY Lp(a)")
print("=" * 70)

outcomes = [
    ('ascvd', 'ASCVD (I20-I25)'),
    ('mi', 'MI (I21-I22)'),
    ('stroke', 'Stroke (I63-I64)'),
    ('tia', 'TIA (G45)'),
    ('pvd', 'PVD (I70-I74)'),
    ('ao_stenosis', 'Aortic stenosis (I35)'),
    ('carotid', 'Carotid disease (I65)'),
]

for group_label, df_g in [("FH", fh), ("Non-FH", nonfh)]:
    print(f"\n  --- {group_label} ---")
    print(f"  {'Outcome':>25s}  {'N events':>9s}  {'Rate%':>6s}  {'Lp(a) yes':>10s}  {'Lp(a) no':>9s}  {'P':>10s}  {'AUC':>5s}")
    print("  " + "-" * 85)

    for col, name in outcomes:
        if col not in df_g.columns:
            continue
        valid = df_g[['lpa', col]].dropna()
        n_events = valid[col].sum()
        rate = 100 * valid[col].mean()

        if n_events >= 10:
            yes_lpa = valid[valid[col] == 1]['lpa'].median()
            no_lpa = valid[valid[col] == 0]['lpa'].median()
            u, p = stats.mannwhitneyu(valid[valid[col]==1]['lpa'], valid[valid[col]==0]['lpa'])
            auc = roc_auc_score(valid[col], valid['lpa'])
            print(f"  {name:>25s}  {n_events:>9.0f}  {rate:>6.1f}  {yes_lpa:>10.1f}  {no_lpa:>9.1f}  {p:>10.2e}  {auc:>5.3f}")
        else:
            print(f"  {name:>25s}  {n_events:>9.0f}  {rate:>6.1f}  insufficient events")

# ============================================================================
# PART 5: SPECIFIC OUTCOME ANALYSIS — Lp(a) QUINTILES
# ============================================================================
print(f"\n{'='*70}")
print("PART 5: Lp(a) QUINTILE × EACH OUTCOME (Non-FH)")
print("=" * 70)

nonfh_q = nonfh.copy()
nonfh_q['lpa_q'] = pd.qcut(nonfh_q['lpa'], 5, labels=['Q1', 'Q2', 'Q3', 'Q4', 'Q5'])

print(f"\n  {'Outcome':>20s}", end='')
for q in ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']:
    print(f"  {q:>7s}", end='')
print(f"  {'Trend P':>10s}")
print("  " + "-" * 70)

for col, name in outcomes:
    if col not in nonfh_q.columns:
        continue
    rates = []
    for q in ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']:
        sub = nonfh_q[nonfh_q['lpa_q'] == q]
        r = 100 * sub[col].mean()
        rates.append(r)

    # Trend test
    rho, p = stats.spearmanr(range(5), rates)
    print(f"  {name:>20s}", end='')
    for r in rates:
        print(f"  {r:>6.2f}%", end='')
    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
    print(f"  {p:>10.4f} {sig}")

# ============================================================================
# WALES FH: DETAILED RISK FACTOR ANALYSIS
# ============================================================================
print(f"\n{'='*70}")
print("WALES FH: RISK FACTORS × Lp(a)")
print("=" * 70)

wales_rf = [
    ('smoking_binary', 'Smoking'),
    ('diabetes', 'Diabetes'),
    ('hypertension', 'Hypertension'),
]

for col, name in wales_rf:
    if col in ar_lpa.columns:
        val = pd.to_numeric(ar_lpa[col], errors='coerce')
        valid = ar_lpa[val.notna()].copy()
        valid['rf'] = val[valid.index]
        yes = valid[valid['rf'] == 1]['lpa_val']
        no = valid[valid['rf'] == 0]['lpa_val']
        if len(yes) >= 5 and len(no) >= 10:
            u, p = stats.mannwhitneyu(yes, no)
            print(f"  {name}: {col}=1 Lp(a)={yes.median():.1f} (n={len(yes)}) vs "
                  f"{col}=0 Lp(a)={no.median():.1f} (n={len(no)}), P={p:.4f}")

# BMI, TRG correlations
for col, name in [('bmi', 'BMI'), ('trig', 'Triglycerides')]:
    if col in ar_lpa.columns:
        val = pd.to_numeric(ar_lpa[col], errors='coerce')
        valid = ar_lpa[val.notna() & ar_lpa['lpa_val'].notna()]
        if len(valid) >= 50:
            rho, p = stats.spearmanr(valid['lpa_val'], val[valid.index])
            print(f"  Lp(a) vs {name}: rho={rho:.3f}, P={p:.4f}")

# ApoB/LDL discordance in Wales
if 'apob_ldl_ratio' in ar_lpa.columns and 'lpa_val' in ar_lpa.columns:
    valid = ar_lpa[ar_lpa['apob_ldl_ratio'].notna() & ar_lpa['lpa_val'].notna()]
    rho, p = stats.spearmanr(valid['lpa_val'], valid['apob_ldl_ratio'])
    print(f"  Lp(a) vs ApoB/LDL ratio (Wales): rho={rho:.3f}, P={p:.4f}")

print("\n" + "=" * 70)
print("COMPLETE")
print("=" * 70)

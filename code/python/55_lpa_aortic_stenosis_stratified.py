#!/usr/bin/env python3
"""
55_lpa_aortic_stenosis_stratified.py
=====================================
Lp(a) vs aortic stenosis in non-FH UKB participants.
Age-stratified, sex-stratified, and interaction analysis.
"""

import pandas as pd
import numpy as np
from scipy import stats
import json
import warnings
warnings.filterwarnings('ignore')

# Load data
print("Loading data...")
lpa = pd.read_csv('D:/calon_lpa_CORRECT.csv')
lpa.columns = [c.replace('participant.', '') for c in lpa.columns]
lpa['eid'] = lpa['eid'].astype(str)
lpa['lpa'] = pd.to_numeric(lpa['p30790_i0'], errors='coerce')

icd = pd.read_csv('D:/calon_backup_data/calon_batch_icd10.csv',
                   usecols=['participant.eid', 'participant.p41270'], dtype=str)
icd.columns = ['eid', 'icd10_raw']

def parse_icd(raw):
    if pd.isna(raw) or raw == '':
        return []
    try:
        codes = json.loads(raw)
        if isinstance(codes, list):
            return [str(c).strip().strip('"') for c in codes]
    except:
        pass
    return [c.strip().strip('"') for c in str(raw).strip('[]"').split(',') if c.strip()]

# Direct string search — avoids memory-heavy JSON parsing
icd['ao_stenosis'] = icd['icd10_raw'].str.contains('"I35', na=False).astype(int)

# Demographics
demo = pd.read_csv('D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_demographics.csv',
                    usecols=['participant.eid', 'participant.p34'],
                    dtype={'participant.eid': str})
demo.columns = [c.replace('participant.', '') for c in demo.columns]
demo['birth_year'] = pd.to_numeric(demo['p34'], errors='coerce')
demo['age_approx'] = 2010 - demo['birth_year']

# Sex (UKB field p31: 0=Female, 1=Male)
sex = pd.read_csv('D:/calon_sex.csv')
sex.columns = [c.replace('participant.', '') for c in sex.columns]
sex.rename(columns={'p31': 'sex'}, inplace=True)
sex['eid'] = sex['eid'].astype(str)
print(f"  Sex: {sex['sex'].value_counts().to_dict()} (0=F, 1=M)")

# FH carriers to exclude
carriers = pd.read_csv(
    'C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis/ukb_carriers_annotated.csv',
    usecols=['eid'], dtype={'eid': str})
fh_eids = set(carriers['eid'].unique())

# Merge
df = lpa[['eid', 'lpa']].merge(icd[['eid', 'ao_stenosis']], on='eid', how='inner')
df = df.merge(demo[['eid', 'age_approx']], on='eid', how='left')
df = df.merge(sex, on='eid', how='left')
df = df[~df['eid'].isin(fh_eids)]
has_lpa = df[df['lpa'].notna()].copy()

print("=" * 70)
print("Lp(a) AND AORTIC STENOSIS - AGE & SEX STRATIFIED")
print("Non-FH UKB participants with Lp(a) data")
print("=" * 70)
print(f"Total: {len(has_lpa):,}")
print(f"Sex: {has_lpa['sex'].value_counts().to_dict()}")
print(f"Age: {has_lpa['age_approx'].min():.0f}-{has_lpa['age_approx'].max():.0f}")
print()

# =====================================================================
# SEX-STRATIFIED
# =====================================================================
print("=" * 70)
print("SEX-STRATIFIED ANALYSIS")
print("=" * 70)

for sex_val, sex_label in [(1, 'Male'), (0, 'Female')]:
    sub = has_lpa[has_lpa['sex'] == sex_val]
    if len(sub) < 100:
        continue

    yes = sub[sub['ao_stenosis'] == 1]['lpa']
    no = sub[sub['ao_stenosis'] == 0]['lpa']

    if len(yes) < 5:
        continue
    u, p = stats.mannwhitneyu(yes, no, alternative='two-sided')
    d = (yes.mean() - no.mean()) / np.sqrt((yes.std()**2 + no.std()**2) / 2)

    print(f"\n{sex_label} (n={len(sub):,}, stenosis={len(yes):,} [{100*len(yes)/len(sub):.2f}%]):")
    print(f"  Stenosis:    median Lp(a) = {yes.median():.1f}, mean = {yes.mean():.1f} nmol/L")
    print(f"  No stenosis: median Lp(a) = {no.median():.1f}, mean = {no.mean():.1f} nmol/L")
    print(f"  P = {p:.2e}, Cohen d = {d:.3f}")

    # Quintile gradient
    sub = sub.copy()
    sub['lpa_q'] = pd.qcut(sub['lpa'], 5, labels=['Q1', 'Q2', 'Q3', 'Q4', 'Q5'])
    print(f"  Quintile gradient:")
    for q in ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']:
        qsub = sub[sub['lpa_q'] == q]
        rate = 100 * qsub['ao_stenosis'].mean()
        print(f"    {q} ({qsub['lpa'].min():.1f}-{qsub['lpa'].max():.1f}): {rate:.2f}%")

    q_rates = sub.groupby('lpa_q')['ao_stenosis'].mean().values
    rho, p_t = stats.spearmanr(range(5), q_rates)
    print(f"  Trend: rho={rho:.3f}, P={p_t:.4f}")

    for thresh in [50, 75, 125]:
        hi_t = sub[sub['lpa'] > thresh]
        lo_t = sub[sub['lpa'] <= thresh]
        r_hi = 100 * hi_t['ao_stenosis'].mean()
        r_lo = 100 * lo_t['ao_stenosis'].mean()
        rr = r_hi / r_lo if r_lo > 0 else float('inf')
        table = np.array([[hi_t['ao_stenosis'].sum(), len(hi_t) - hi_t['ao_stenosis'].sum()],
                          [lo_t['ao_stenosis'].sum(), len(lo_t) - lo_t['ao_stenosis'].sum()]])
        chi2, p_c, _, _ = stats.chi2_contingency(table)
        print(f"  Lp(a)>{thresh}: {r_hi:.2f}% vs {r_lo:.2f}%, RR={rr:.2f}, P={p_c:.2e}")

# =====================================================================
# AGE-STRATIFIED
# =====================================================================
print()
print("=" * 70)
print("AGE-STRATIFIED ANALYSIS")
print("=" * 70)

age_bins = [(40, 50, '40-49'), (50, 55, '50-54'), (55, 60, '55-59'),
            (60, 65, '60-64'), (65, 70, '65-69'), (70, 80, '70+')]

header = f"{'Age':>8s}  {'N':>8s}  {'Stenosis':>8s}  {'Rate%':>6s}  {'Lpa(st)':>8s}  {'Lpa(no)':>8s}  {'P':>10s}  {'d':>6s}"
print(f"\n{header}")
print("-" * 75)

for lo, hi, label in age_bins:
    sub = has_lpa[(has_lpa['age_approx'] >= lo) & (has_lpa['age_approx'] < hi)]
    if len(sub) < 100:
        continue

    yes = sub[sub['ao_stenosis'] == 1]['lpa']
    no = sub[sub['ao_stenosis'] == 0]['lpa']
    n_st = len(yes)
    rate = 100 * n_st / len(sub)

    if len(yes) >= 10 and len(no) >= 10:
        u, p = stats.mannwhitneyu(yes, no, alternative='two-sided')
        d = (yes.mean() - no.mean()) / np.sqrt((yes.std()**2 + no.std()**2) / 2)
        p_str = f"{p:.2e}"
        d_str = f"{d:.3f}"
    else:
        p_str = "N/A"
        d_str = "N/A"

    y_med = f"{yes.median():.1f}" if len(yes) > 0 else "N/A"
    n_med = f"{no.median():.1f}" if len(no) > 0 else "N/A"

    print(f"{label:>8s}  {len(sub):>8,}  {n_st:>8,}  {rate:>6.2f}  {y_med:>8s}  {n_med:>8s}  {p_str:>10s}  {d_str:>6s}")

# =====================================================================
# AGE x Lp(a) INTERACTION
# =====================================================================
print()
print("=" * 70)
print("AGE x Lp(a) INTERACTION - STENOSIS RATES")
print("=" * 70)

age_med = has_lpa['age_approx'].median()
lpa_med = has_lpa['lpa'].median()
print(f"Age median: {age_med:.0f}, Lp(a) median: {lpa_med:.1f} nmol/L")
print()

has_lpa['age_group'] = np.where(has_lpa['age_approx'] >= age_med, 'Older', 'Younger')
has_lpa['lpa_group'] = np.where(has_lpa['lpa'] > lpa_med, 'High Lp(a)', 'Low Lp(a)')

print(f"{'Group':>25s}  {'N':>8s}  {'Stenosis':>8s}  {'Rate%':>7s}")
print("-" * 55)
for age_g in ['Younger', 'Older']:
    for lpa_g in ['Low Lp(a)', 'High Lp(a)']:
        sub = has_lpa[(has_lpa['age_group'] == age_g) & (has_lpa['lpa_group'] == lpa_g)]
        n_st = sub['ao_stenosis'].sum()
        rate = 100 * n_st / len(sub)
        print(f"{age_g + ' + ' + lpa_g:>25s}  {len(sub):>8,}  {n_st:>8,}  {rate:>7.2f}")

# Fold gradient
young_lo = has_lpa[(has_lpa['age_group'] == 'Younger') & (has_lpa['lpa_group'] == 'Low Lp(a)')]['ao_stenosis'].mean()
old_hi = has_lpa[(has_lpa['age_group'] == 'Older') & (has_lpa['lpa_group'] == 'High Lp(a)')]['ao_stenosis'].mean()
if young_lo > 0:
    print(f"\nOlder+High Lp(a) vs Younger+Low Lp(a): {old_hi/young_lo:.1f}-fold gradient")

# =====================================================================
# Lp(a) QUINTILE x AGE TERTILE HEATMAP
# =====================================================================
print()
print("=" * 70)
print("Lp(a) QUINTILES x AGE TERTILES - STENOSIS RATE HEATMAP")
print("=" * 70)

has_lpa['lpa_q'] = pd.qcut(has_lpa['lpa'], 5, labels=['Q1', 'Q2', 'Q3', 'Q4', 'Q5'])
has_lpa['age_t'] = pd.qcut(has_lpa['age_approx'], 3, labels=['Young', 'Middle', 'Old'])

print(f"\n{'':>8s}", end='')
for q in ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']:
    print(f"  {q:>7s}", end='')
print()
print("-" * 50)

for at in ['Young', 'Middle', 'Old']:
    print(f"{at:>8s}", end='')
    for q in ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']:
        sub = has_lpa[(has_lpa['age_t'] == at) & (has_lpa['lpa_q'] == q)]
        rate = 100 * sub['ao_stenosis'].mean() if len(sub) > 0 else 0
        print(f"  {rate:>6.2f}%", end='')
    print()

# Max gradient
young_q1 = has_lpa[(has_lpa['age_t'] == 'Young') & (has_lpa['lpa_q'] == 'Q1')]
old_q5 = has_lpa[(has_lpa['age_t'] == 'Old') & (has_lpa['lpa_q'] == 'Q5')]
r_min = young_q1['ao_stenosis'].mean() * 100
r_max = old_q5['ao_stenosis'].mean() * 100
if r_min > 0:
    print(f"\nMax gradient: Old Q5 ({r_max:.2f}%) vs Young Q1 ({r_min:.2f}%) = {r_max/r_min:.1f}-fold")

# =====================================================================
# FH vs NON-FH COMPARISON
# =====================================================================
print()
print("=" * 70)
print("COMPARISON: FH CARRIERS vs NON-FH POPULATION")
print("=" * 70)

fh_lpa = lpa[lpa['eid'].isin(fh_eids) & lpa['lpa'].notna()][['eid', 'lpa']]
fh_icd = icd[icd['eid'].isin(fh_eids)][['eid', 'ao_stenosis']]
fh_m = fh_lpa.merge(fh_icd, on='eid', how='inner')

for label, dataset in [('Non-FH (n=372,830)', has_lpa), ('FH carriers', fh_m)]:
    yes_v = dataset[dataset['ao_stenosis'] == 1]['lpa']
    no_v = dataset[dataset['ao_stenosis'] == 0]['lpa']
    n_st = len(yes_v)
    rate = 100 * n_st / len(dataset) if len(dataset) > 0 else 0

    if len(yes_v) >= 5 and len(no_v) >= 10:
        u, p = stats.mannwhitneyu(yes_v, no_v, alternative='two-sided')
        p_str = f"{p:.2e}"
    else:
        p_str = "N/A (n too small)"

    y_med = f"{yes_v.median():.1f}" if len(yes_v) > 0 else "N/A"
    n_med = f"{no_v.median():.1f}"

    print(f"\n{label}:")
    print(f"  N={len(dataset):,}, stenosis={n_st} ({rate:.2f}%)")
    print(f"  Lp(a) with stenosis: median {y_med} nmol/L")
    print(f"  Lp(a) without:       median {n_med} nmol/L")
    print(f"  P = {p_str}")

# =====================================================================
# LOGISTIC REGRESSION — age + Lp(a) + interaction
# =====================================================================
print()
print("=" * 70)
print("LOGISTIC REGRESSION (age-adjusted)")
print("=" * 70)

from sklearn.linear_model import LogisticRegression

reg_df = has_lpa[['lpa', 'age_approx', 'sex', 'ao_stenosis']].dropna().copy()
print(f"N for regression: {len(reg_df):,}")

reg_df['lpa_std'] = (reg_df['lpa'] - reg_df['lpa'].mean()) / reg_df['lpa'].std()
reg_df['age_std'] = (reg_df['age_approx'] - reg_df['age_approx'].mean()) / reg_df['age_approx'].std()
reg_df['age_x_lpa'] = reg_df['age_std'] * reg_df['lpa_std']
reg_df['sex_x_lpa'] = reg_df['sex'] * reg_df['lpa_std']

for model_name, features in [
    ('Model 1: Lp(a) only', ['lpa_std']),
    ('Model 2: Age + sex + Lp(a)', ['age_std', 'sex', 'lpa_std']),
    ('Model 3: + age x Lp(a)', ['age_std', 'sex', 'lpa_std', 'age_x_lpa']),
    ('Model 4: + sex x Lp(a)', ['age_std', 'sex', 'lpa_std', 'age_x_lpa', 'sex_x_lpa']),
]:
    X = reg_df[features].values
    y = reg_df['ao_stenosis'].values
    lr = LogisticRegression(max_iter=1000, solver='lbfgs')
    lr.fit(X, y)

    print(f"\n{model_name}:")
    for feat, coef in zip(features, lr.coef_[0]):
        odds = np.exp(coef)
        print(f"  {feat:>12s}: OR = {odds:.3f}  (beta={coef:.4f})")

print()
print("=" * 70)
print("COMPLETE")
print("=" * 70)

#!/usr/bin/env python3
"""
73_lpa_special_populations.py
==============================
1. LPA genotype proxy analysis using Lp(a) levels (90% heritable)
2. Lp(a) in pregnancy-related conditions (UKB has pregnancy history)
3. Lp(a) by menopausal status / HRT use
4. Medication effects on Lp(a) (statins, HRT, oral contraceptives)
5. Lp(a) by age decade (proxy for menopause effect)

Author: Dr Nader Genedy
Date:   April 2026
"""

import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')
import os

DATA = r"D:/CALON_AF3_PROJECT/data"
ANALYSIS = r"C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis"

print("=" * 70)
print("Lp(a) IN SPECIAL POPULATIONS + DRUG EFFECTS")
print("=" * 70)

# Load
lpa = pd.read_csv(f"{DATA}/ukb_wide/calon_lpa_CORRECT.csv")
lpa.columns = [c.replace('participant.', '') for c in lpa.columns]
lpa['eid'] = lpa['eid'].astype(str)
lpa['lpa'] = pd.to_numeric(lpa['p30790_i0'], errors='coerce')

sex = pd.read_csv(f"{DATA}/ukb_wide/calon_sex.csv")
sex.columns = [c.replace('participant.', '') for c in sex.columns]
sex.rename(columns={'p31': 'sex'}, inplace=True)
sex['eid'] = sex['eid'].astype(str)

demo = pd.read_csv(f"{DATA}/ukb_wide/ukb_reviewer_demographics.csv")
demo.columns = [c.replace('participant.', '') for c in demo.columns]
demo['eid'] = demo['eid'].astype(str)
demo['age'] = 2010 - pd.to_numeric(demo['p34'], errors='coerce')

icd = pd.read_csv(f"{DATA}/ukb_wide/calon_batch_icd10.csv",
                   usecols=['participant.eid', 'participant.p41270'], dtype=str)
icd.columns = ['eid', 'icd10_raw']
# Pregnancy complications
icd['preeclampsia'] = icd['icd10_raw'].str.contains('"O1[0-6]', na=False, regex=True).astype(int)
icd['gest_diabetes'] = icd['icd10_raw'].str.contains('"O24', na=False, regex=True).astype(int)
icd['pregnancy_any'] = icd['icd10_raw'].str.contains('"O', na=False, regex=True).astype(int)
icd['menopause'] = icd['icd10_raw'].str.contains('"N95|"E28', na=False, regex=True).astype(int)
icd['ascvd'] = icd['icd10_raw'].str.contains('"I2[0-5]', na=False, regex=True).astype(int)

carriers = pd.read_csv(os.path.join(ANALYSIS, "ukb_carriers_annotated.csv"), dtype={'eid': str})
fh_eids = set(carriers.sort_values('sss', ascending=False).drop_duplicates('eid')['eid'])

# Merge
df = lpa[['eid', 'lpa']].merge(sex[['eid', 'sex']], on='eid', how='left')
df = df.merge(demo[['eid', 'age']], on='eid', how='left')
df = df.merge(icd[['eid', 'preeclampsia', 'gest_diabetes', 'pregnancy_any',
                    'menopause', 'ascvd']], on='eid', how='left')
df['is_fh'] = df['eid'].isin(fh_eids).astype(int)
df = df[df['lpa'].notna()].copy()

women = df[df['sex'] == 0].copy()
men = df[df['sex'] == 1].copy()

print(f"Total with Lp(a): {len(df)}")
print(f"Women: {len(women)}, Men: {len(men)}")

# ============================================================================
# PART 1: LPA GENOTYPE PROXY — Lp(a) LEVEL CATEGORIES
# ============================================================================
print(f"\n{'='*70}")
print("PART 1: LPA GENOTYPE PROXY (Lp(a) is 90% heritable)")
print("=" * 70)

# Classify by Lp(a) level as proxy for LPA genotype
df['lpa_geno'] = 'Normal (<30)'
df.loc[df['lpa'] >= 30, 'lpa_geno'] = 'Elevated (30-75)'
df.loc[df['lpa'] >= 75, 'lpa_geno'] = 'High (75-125)'
df.loc[df['lpa'] >= 125, 'lpa_geno'] = 'Very High (>125)'

print(f"\n  LPA genotype proxy distribution:")
for g in ['Normal (<30)', 'Elevated (30-75)', 'High (75-125)', 'Very High (>125)']:
    sub = df[df['lpa_geno'] == g]
    ascvd_r = 100 * sub['ascvd'].mean()
    print(f"    {g:>20s}: n={len(sub):>7,} ({100*len(sub)/len(df):.1f}%), "
          f"median Lp(a)={sub['lpa'].median():.0f}, ASCVD={ascvd_r:.1f}%")

# In FH vs non-FH
print(f"\n  By FH status:")
for fh_val, lab in [(1, 'FH'), (0, 'Non-FH')]:
    sub = df[df['is_fh'] == fh_val]
    for g in ['Normal (<30)', 'Elevated (30-75)', 'High (75-125)', 'Very High (>125)']:
        n = (sub['lpa_geno'] == g).sum()
        pct = 100 * n / len(sub) if len(sub) > 0 else 0
        print(f"    {lab} {g}: {n} ({pct:.1f}%)")

# ============================================================================
# PART 2: Lp(a) BY SEX AND AGE DECADE
# ============================================================================
print(f"\n{'='*70}")
print("PART 2: Lp(a) BY SEX AND AGE DECADE")
print("(Proxy for menopausal effect)")
print("=" * 70)

print(f"\n  {'Age':>6s}  {'Women n':>8s}  {'Women Lpa':>10s}  {'Men n':>8s}  {'Men Lpa':>10s}  {'W-M diff':>9s}  {'P':>10s}")
print("  " + "-" * 70)

for lo, hi, lab in [(40, 45, '40-44'), (45, 50, '45-49'), (50, 55, '50-54'),
                     (55, 60, '55-59'), (60, 65, '60-64'), (65, 70, '65-69')]:
    w = women[(women['age'] >= lo) & (women['age'] < hi)]['lpa']
    m = men[(men['age'] >= lo) & (men['age'] < hi)]['lpa']
    if len(w) >= 50 and len(m) >= 50:
        u, p = stats.mannwhitneyu(w, m)
        diff = w.median() - m.median()
        print(f"  {lab:>6s}  {len(w):>8,}  {w.median():>10.1f}  {len(m):>8,}  {m.median():>10.1f}  {diff:>+9.1f}  {p:>10.4f}")

# The menopause effect: Lp(a) typically rises 20-50% after menopause
# Compare pre vs post menopausal age groups in women
pre = women[(women['age'] >= 40) & (women['age'] < 50)]['lpa']
post = women[(women['age'] >= 55) & (women['age'] < 65)]['lpa']
if len(pre) >= 100 and len(post) >= 100:
    u, p = stats.mannwhitneyu(pre, post)
    pct_change = 100 * (post.median() - pre.median()) / pre.median()
    print(f"\n  Pre-menopausal (40-49) vs Post-menopausal (55-64):")
    print(f"    Pre:  n={len(pre):,}, median={pre.median():.1f}")
    print(f"    Post: n={len(post):,}, median={post.median():.1f}")
    print(f"    Change: {pct_change:+.1f}%, P={p:.4f}")

# ============================================================================
# PART 3: PREGNANCY COMPLICATIONS AND Lp(a)
# ============================================================================
print(f"\n{'='*70}")
print("PART 3: Lp(a) AND PREGNANCY COMPLICATIONS")
print("=" * 70)

for cond, name in [('preeclampsia', 'Pre-eclampsia (O10-O16)'),
                    ('gest_diabetes', 'Gestational diabetes (O24)'),
                    ('pregnancy_any', 'Any pregnancy code (O)')]:
    yes = women[women[cond] == 1]['lpa']
    no = women[women[cond] == 0]['lpa']
    if len(yes) >= 20:
        u, p = stats.mannwhitneyu(yes, no)
        d = (yes.mean() - no.mean()) / np.sqrt((yes.std()**2 + no.std()**2) / 2)
        print(f"  {name}:")
        print(f"    With: n={len(yes):,}, median Lp(a)={yes.median():.1f}")
        print(f"    Without: n={len(no):,}, median Lp(a)={no.median():.1f}")
        print(f"    P={p:.4f}, d={d:.3f}")

# Pre-eclampsia + high Lp(a) -> ASCVD risk?
women_pe = women[women['preeclampsia'].notna()].copy()
women_pe['high_lpa'] = (women_pe['lpa'] > 50).astype(int)
if women_pe['preeclampsia'].sum() >= 20:
    print(f"\n  Pre-eclampsia x Lp(a) -> ASCVD:")
    for pe in [0, 1]:
        for lpa_hi in [0, 1]:
            sub = women_pe[(women_pe['preeclampsia'] == pe) & (women_pe['high_lpa'] == lpa_hi)]
            if len(sub) >= 20:
                r = 100 * sub['ascvd'].mean()
                pe_lab = 'Pre-eclampsia' if pe else 'No pre-eclampsia'
                lpa_lab = 'High Lp(a)' if lpa_hi else 'Low Lp(a)'
                print(f"    {pe_lab} + {lpa_lab}: ASCVD={r:.1f}% (n={len(sub):,})")

# ============================================================================
# PART 4: MEDICATION EFFECTS ON Lp(a)
# ============================================================================
print(f"\n{'='*70}")
print("PART 4: KNOWN DRUG EFFECTS ON Lp(a) (Literature)")
print("=" * 70)

print("""
  DRUGS THAT RAISE Lp(a):
    - Niacin: paradoxically raises Lp(a) in some patients
    - Growth hormone: increases Lp(a) 50-100%
    - Isotretinoin (Accutane): reported to increase Lp(a)
    - Androgens/testosterone: can increase Lp(a)

  DRUGS THAT LOWER Lp(a):
    - PCSK9 inhibitors (evolocumab, alirocumab): -20-30%
    - Pelacarsen (ASO targeting LPA mRNA): -35-80%
    - Olpasiran (siRNA targeting LPA mRNA): -70-101%
    - Lepodisiran (siRNA): -55-97%
    - Muvalaplin (oral small molecule, KIV-7/8 target): -45-65%
    - Lipoprotein apheresis: -60-75% (acute, rebounds)
    - Estrogen/HRT: -20-40% (clinically significant)
    - Tamoxifen: -20%

  DRUGS WITH NO CONSISTENT EFFECT ON Lp(a):
    - Statins: variable, no net effect (some increase, some decrease)
    - Ezetimibe: no significant effect
    - Fibrates: minimal effect
    - SGLT2 inhibitors: no effect
    - GLP-1 agonists: no consistent effect
    - Metformin: no effect

  CRITICAL FOR FH PATIENTS:
    - HRT (estrogen) lowers Lp(a) 20-40% — this is clinically relevant
      for post-menopausal FH women with high Lp(a)
    - PCSK9i lowers both LDL and Lp(a) — dual benefit
    - Statins do NOT lower Lp(a) — our Saturation Threshold explains why
""")

# Can we see the HRT effect in UKB data?
# HRT use might be in medication data
# For now, use age as proxy: post-menopausal women may have higher Lp(a)

# ============================================================================
# PART 5: Lp(a) IN FH WOMEN — MENOPAUSAL EFFECT
# ============================================================================
print(f"{'='*70}")
print("PART 5: Lp(a) IN FH WOMEN BY AGE")
print("=" * 70)

fh_women = women[women['is_fh'] == 1].copy()
nonfh_women = women[women['is_fh'] == 0].copy()

print(f"  FH women with Lp(a): {len(fh_women)}")
print(f"  Non-FH women with Lp(a): {len(nonfh_women)}")

for lab, w_df in [("FH Women", fh_women), ("Non-FH Women", nonfh_women)]:
    print(f"\n  --- {lab} ---")
    for lo, hi, age_lab in [(40, 50, '40-49 (pre-meno)'), (50, 55, '50-54 (peri)'),
                             (55, 65, '55-64 (post)'), (65, 75, '65-74 (late post)')]:
        sub = w_df[(w_df['age'] >= lo) & (w_df['age'] < hi)]
        if len(sub) >= 10:
            print(f"    {age_lab:>25s}: n={len(sub):>5}, median Lp(a)={sub['lpa'].median():.1f}, "
                  f"ASCVD={100*sub['ascvd'].mean():.1f}%")

# Menopause effect in FH women
fh_pre = fh_women[(fh_women['age'] >= 40) & (fh_women['age'] < 50)]['lpa']
fh_post = fh_women[(fh_women['age'] >= 55) & (fh_women['age'] < 65)]['lpa']
if len(fh_pre) >= 10 and len(fh_post) >= 10:
    u, p = stats.mannwhitneyu(fh_pre, fh_post)
    pct = 100 * (fh_post.median() - fh_pre.median()) / fh_pre.median()
    print(f"\n  FH Women menopause effect:")
    print(f"    Pre (40-49): median={fh_pre.median():.1f} (n={len(fh_pre)})")
    print(f"    Post (55-64): median={fh_post.median():.1f} (n={len(fh_post)})")
    print(f"    Change: {pct:+.1f}%, P={p:.4f}")

print("\n" + "=" * 70)
print("COMPLETE")
print("=" * 70)

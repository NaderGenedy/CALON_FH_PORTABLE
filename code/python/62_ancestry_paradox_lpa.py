#!/usr/bin/env python3
"""
62_ancestry_paradox_lpa.py
===========================
The Ancestry-Lp(a)-ASCVD Paradox in FH:
  Black carriers: highest Lp(a) (65 nmol/L) but LOWEST ASCVD (6.6%)
  White carriers: lowest Lp(a) (21 nmol/L) but HIGHEST ASCVD (13.8%)

Is this real biology or confounding? Test with age, sex, treatment adjustment.

Author: Dr Nader Genedy
Date:   April 2026
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.linear_model import LogisticRegression
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

ANALYSIS = r"C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis"
FIGURES = os.path.join(ANALYSIS, "figures")

print("=" * 70)
print("ANCESTRY-Lp(a)-ASCVD PARADOX IN FH CARRIERS")
print("=" * 70)

# Load and merge
carriers = pd.read_csv(os.path.join(ANALYSIS, "ukb_carriers_annotated.csv"), dtype={'eid': str})
patients = carriers.sort_values('sss', ascending=False).drop_duplicates('eid')

lpa = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/calon_lpa_CORRECT.csv")
lpa.columns = [c.replace('participant.', '') for c in lpa.columns]
lpa['eid'] = lpa['eid'].astype(str)
lpa['lpa'] = pd.to_numeric(lpa['p30790_i0'], errors='coerce')

sex = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/calon_sex.csv")
sex.columns = [c.replace('participant.', '') for c in sex.columns]
sex.rename(columns={'p31': 'sex'}, inplace=True)
sex['eid'] = sex['eid'].astype(str)

demo = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/ukb_reviewer_demographics.csv")
demo.columns = [c.replace('participant.', '') for c in demo.columns]
demo['eid'] = demo['eid'].astype(str)
demo['age'] = 2010 - pd.to_numeric(demo['p34'], errors='coerce')
demo['ethnicity_raw'] = pd.to_numeric(demo['p21000_i0'], errors='coerce')
demo['bmi'] = pd.to_numeric(demo['p21001_i0'], errors='coerce')
demo['ethnicity'] = 'Other'
demo.loc[demo['ethnicity_raw'].isin([1001, 1002, 1003]), 'ethnicity'] = 'White'
demo.loc[demo['ethnicity_raw'].isin([3001, 3002, 3003, 3004]), 'ethnicity'] = 'Asian'
demo.loc[demo['ethnicity_raw'].isin([4001, 4002, 4003]), 'ethnicity'] = 'Black'

icd = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/calon_batch_icd10.csv",
                   usecols=['participant.eid', 'participant.p41270'], dtype=str)
icd.columns = ['eid', 'icd10_raw']
icd['ascvd'] = icd['icd10_raw'].str.contains('"I2[0-5]', na=False, regex=True).astype(int)

smoke = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/ukb_reviewer_smoking.csv")
smoke.columns = [c.replace('participant.', '') for c in smoke.columns]
smoke['eid'] = smoke['eid'].astype(str)
smoke['smoking'] = (pd.to_numeric(smoke['p20116_i0'], errors='coerce') == 2).astype(int)

ldl_df = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/ukb_reviewer_longitudinal_lipids.csv",
                      usecols=['participant.eid', 'participant.p30780_i0'])
ldl_df.columns = [c.replace('participant.', '') for c in ldl_df.columns]
ldl_df['eid'] = ldl_df['eid'].astype(str)
ldl_df['ldl'] = pd.to_numeric(ldl_df['p30780_i0'], errors='coerce')

# Merge all
df = patients[['eid', 'ldlr_domain', 'gene', 'consequence']].copy()
for d, cols in [(lpa, ['eid', 'lpa']), (sex, ['eid', 'sex']),
                (demo, ['eid', 'age', 'ethnicity', 'bmi']),
                (icd, ['eid', 'ascvd']), (smoke, ['eid', 'smoking']),
                (ldl_df, ['eid', 'ldl'])]:
    df = df.merge(d[cols], on='eid', how='left')

fh = df[df['lpa'].notna()].copy()
print(f"FH carriers with Lp(a): {len(fh)}")

# =====================================================================
# 1. THE RAW PARADOX
# =====================================================================
print("\n" + "=" * 70)
print("1. THE RAW PARADOX")
print("=" * 70)

print(f"\n  {'Ethnicity':>10s}  {'N':>6s}  {'Med Lpa':>8s}  {'Mean Lpa':>9s}  {'ASCVD%':>7s}  {'Med Age':>8s}  {'%Male':>6s}  {'Med LDL':>8s}")
print("  " + "-" * 75)
for eth in ['White', 'Black', 'Asian', 'Other']:
    sub = fh[fh['ethnicity'] == eth]
    if len(sub) >= 10:
        print(f"  {eth:>10s}  {len(sub):>6d}  {sub['lpa'].median():>8.1f}  {sub['lpa'].mean():>9.1f}  "
              f"{100*sub['ascvd'].mean():>7.1f}  {sub['age'].median():>8.0f}  "
              f"{100*sub['sex'].mean():>5.0f}%  {sub['ldl'].median():>8.2f}")

# =====================================================================
# 2. AGE EXPLAINS IT?
# =====================================================================
print("\n" + "=" * 70)
print("2. IS AGE THE EXPLANATION?")
print("=" * 70)

white = fh[fh['ethnicity'] == 'White']
black = fh[fh['ethnicity'] == 'Black']

print(f"\n  White age: median={white['age'].median():.0f}, mean={white['age'].mean():.1f}")
print(f"  Black age: median={black['age'].median():.0f}, mean={black['age'].mean():.1f}")
u, p = stats.mannwhitneyu(white['age'].dropna(), black['age'].dropna())
print(f"  Age difference P={p:.4f}")

# Age-stratified ASCVD by ethnicity
print("\n  Age-stratified ASCVD rates:")
for lo, hi, lab in [(40, 55, '<55'), (55, 65, '55-64'), (65, 80, '65+')]:
    for eth in ['White', 'Black']:
        sub = fh[(fh['ethnicity'] == eth) & (fh['age'] >= lo) & (fh['age'] < hi)]
        if len(sub) >= 10:
            r = 100 * sub['ascvd'].mean()
            print(f"    {eth:>6s} age {lab}: n={len(sub):>4d}, ASCVD={r:.1f}%, Lp(a)={sub['lpa'].median():.1f}")

# =====================================================================
# 3. MUTATION SPECTRUM DIFFERS?
# =====================================================================
print("\n" + "=" * 70)
print("3. DO BLACK AND WHITE CARRIERS HAVE DIFFERENT MUTATIONS?")
print("=" * 70)

for eth in ['White', 'Black']:
    sub = fh[fh['ethnicity'] == eth]
    cons = sub['consequence'].dropna()
    if len(cons) >= 10:
        miss_pct = 100 * cons.str.contains('missense', na=False).mean()
        lof_pct = 100 * cons.str.contains('frameshift|stop_gained', na=False).mean()
        splice_pct = 100 * cons.str.contains('splice', na=False).mean()
        print(f"  {eth}: missense={miss_pct:.1f}%, LoF={lof_pct:.1f}%, splice={splice_pct:.1f}%")

# Domain distribution
print("\n  Domain distribution:")
for eth in ['White', 'Black']:
    sub = fh[(fh['ethnicity'] == eth) & fh['ldlr_domain'].notna()]
    if len(sub) >= 10:
        top = sub['ldlr_domain'].value_counts().head(5)
        print(f"  {eth} (n={len(sub)}):")
        for dom, n in top.items():
            short = dom.replace('Ligand-binding ', 'LB-')
            print(f"    {short}: {n} ({100*n/len(sub):.1f}%)")

# =====================================================================
# 4. ADJUSTED LOGISTIC REGRESSION
# =====================================================================
print("\n" + "=" * 70)
print("4. ADJUSTED ASCVD ~ ETHNICITY + CONFOUNDERS")
print("=" * 70)

reg = fh[['ascvd', 'ethnicity', 'age', 'sex', 'lpa', 'ldl', 'bmi', 'smoking']].dropna().copy()
reg['is_black'] = (reg['ethnicity'] == 'Black').astype(int)
reg['is_asian'] = (reg['ethnicity'] == 'Asian').astype(int)
reg['log_lpa'] = np.log1p(reg['lpa'])

models = [
    ('Ethnicity only', ['is_black', 'is_asian']),
    ('+ Age + Sex', ['is_black', 'is_asian', 'age', 'sex']),
    ('+ Age + Sex + LDL', ['is_black', 'is_asian', 'age', 'sex', 'ldl']),
    ('+ Age + Sex + LDL + Lp(a)', ['is_black', 'is_asian', 'age', 'sex', 'ldl', 'log_lpa']),
    ('Full (+ BMI + smoking)', ['is_black', 'is_asian', 'age', 'sex', 'ldl', 'log_lpa', 'bmi', 'smoking']),
]

print(f"\n  n={len(reg)}, ASCVD events={reg['ascvd'].sum()}")
print(f"\n  {'Model':>30s}  {'Black OR':>9s}  {'Asian OR':>9s}  {'Lp(a) OR':>9s}")
print("  " + "-" * 65)

for name, feats in models:
    X = reg[feats].values
    y = reg['ascvd'].values
    lr = LogisticRegression(max_iter=1000, solver='lbfgs')
    lr.fit(X, y)

    black_or = np.exp(lr.coef_[0][0])
    asian_or = np.exp(lr.coef_[0][1])
    lpa_or = np.exp(lr.coef_[0][feats.index('log_lpa')]) if 'log_lpa' in feats else None
    lpa_str = f"{lpa_or:.3f}" if lpa_or else "N/A"
    print(f"  {name:>30s}  {black_or:>9.3f}  {asian_or:>9.3f}  {lpa_str:>9s}")

# =====================================================================
# 5. Lp(a) PREDICTIVE VALUE BY ETHNICITY
# =====================================================================
print("\n" + "=" * 70)
print("5. Lp(a) PREDICTS ASCVD DIFFERENTLY BY ETHNICITY?")
print("=" * 70)

for eth in ['White', 'Black', 'Asian']:
    sub = fh[(fh['ethnicity'] == eth) & fh['ascvd'].notna()]
    if len(sub) >= 30 and sub['ascvd'].sum() >= 5:
        rho, p = stats.spearmanr(sub['lpa'], sub['ascvd'])
        try:
            from sklearn.metrics import roc_auc_score
            auc = roc_auc_score(sub['ascvd'], sub['lpa'])
            print(f"  {eth:>6s}: n={len(sub)}, ASCVD={sub['ascvd'].sum()}, "
                  f"rho={rho:.3f}, P={p:.4f}, AUC={auc:.3f}")
        except:
            print(f"  {eth:>6s}: n={len(sub)}, rho={rho:.3f}, P={p:.4f}")

# =====================================================================
# 6. NON-FH COMPARISON
# =====================================================================
print("\n" + "=" * 70)
print("6. NON-FH POPULATION: SAME PARADOX?")
print("=" * 70)

# Check if non-FH also shows paradox
nonfh = icd[~icd['eid'].isin(set(patients['eid']))].copy()
nonfh = nonfh.merge(lpa[['eid', 'lpa']], on='eid', how='inner')
nonfh = nonfh.merge(demo[['eid', 'ethnicity', 'age']], on='eid', how='left')
nonfh = nonfh[nonfh['lpa'].notna()]

for eth in ['White', 'Black', 'Asian']:
    sub = nonfh[nonfh['ethnicity'] == eth]
    if len(sub) >= 100:
        print(f"  {eth:>6s}: n={len(sub):>6,}, Lp(a)={sub['lpa'].median():.1f}, "
              f"ASCVD={100*sub['ascvd'].mean():.1f}%, age={sub['age'].median():.0f}")

# =====================================================================
# FIGURE
# =====================================================================
print(f"\n{'='*70}")
print("Generating figure...")

fig, axes = plt.subplots(2, 2, figsize=(14, 11))
RED = '#b2182b'
BLUE = '#2166ac'
GREEN = '#1b7837'
ORANGE = '#ef8a62'

# Panel A: Lp(a) by ethnicity (boxplot)
ax = axes[0, 0]
eth_data = []
eth_labels = []
eth_cols = []
for eth, col in [('White', BLUE), ('Black', RED), ('Asian', ORANGE)]:
    vals = fh[fh['ethnicity'] == eth]['lpa'].dropna()
    if len(vals) >= 10:
        eth_data.append(vals.values)
        eth_labels.append(f"{eth}\n(n={len(vals)})")
        eth_cols.append(col)

if eth_data:
    bp = ax.boxplot(eth_data, labels=eth_labels, patch_artist=True,
                    medianprops=dict(color='black', linewidth=2), showfliers=False)
    for patch, c in zip(bp['boxes'], eth_cols):
        patch.set_facecolor(c)
        patch.set_alpha(0.7)
    ax.set_ylabel('Lp(a) (nmol/L)', fontsize=10)
    ax.set_title('A. Lp(a) by Ethnicity (FH Carriers)', fontsize=12, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: ASCVD rate by ethnicity
ax = axes[0, 1]
eth_rates = []
for eth in ['White', 'Black', 'Asian']:
    sub = fh[fh['ethnicity'] == eth]
    if len(sub) >= 10:
        eth_rates.append((eth, 100 * sub['ascvd'].mean(), len(sub)))

if eth_rates:
    x_pos = range(len(eth_rates))
    bars = ax.bar(x_pos, [r[1] for r in eth_rates],
                  color=[BLUE, RED, ORANGE][:len(eth_rates)],
                  edgecolor='black', linewidth=0.5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"{r[0]}\n(n={r[2]})" for r in eth_rates])
    ax.set_ylabel('ASCVD Rate (%)', fontsize=10)
    ax.set_title('B. ASCVD by Ethnicity (FH Carriers)', fontsize=12, fontweight='bold')
    for bar, r in zip(bars, eth_rates):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                f'{r[1]:.1f}%', ha='center', fontsize=10, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: THE PARADOX — Lp(a) vs ASCVD scatter by ethnicity
ax = axes[1, 0]
for eth, col, marker in [('White', BLUE, 'o'), ('Black', RED, 's'), ('Asian', ORANGE, '^')]:
    sub = fh[fh['ethnicity'] == eth]
    if len(sub) >= 10:
        lpa_med = sub['lpa'].median()
        ascvd_r = 100 * sub['ascvd'].mean()
        ax.scatter(lpa_med, ascvd_r, c=col, s=len(sub)/2, marker=marker,
                   edgecolors='black', linewidth=1, label=f'{eth} (n={len(sub)})', zorder=5)

ax.set_xlabel('Median Lp(a) (nmol/L)', fontsize=10)
ax.set_ylabel('ASCVD Rate (%)', fontsize=10)
ax.set_title('C. THE PARADOX: Highest Lp(a) = Lowest ASCVD', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)
# Add arrow showing paradox
ax.annotate('', xy=(65, 7), xytext=(22, 14),
            arrowprops=dict(arrowstyle='->', color='red', lw=2))
ax.text(40, 11, 'Expected\ndirection', fontsize=8, color='red', ha='center')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Age distribution by ethnicity
ax = axes[1, 1]
for eth, col in [('White', BLUE), ('Black', RED), ('Asian', ORANGE)]:
    vals = fh[fh['ethnicity'] == eth]['age'].dropna()
    if len(vals) >= 10:
        ax.hist(vals, bins=20, alpha=0.5, color=col, label=f'{eth} (med={vals.median():.0f})',
                density=True, edgecolor='white', linewidth=0.3)

ax.set_xlabel('Age at Recruitment', fontsize=10)
ax.set_ylabel('Density', fontsize=10)
ax.set_title('D. Age Distribution by Ethnicity', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout(pad=2.0)
fig_path = os.path.join(FIGURES, "Figure_Ancestry_Paradox.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {fig_path}")

print("\n" + "=" * 70)
print("COMPLETE")
print("=" * 70)

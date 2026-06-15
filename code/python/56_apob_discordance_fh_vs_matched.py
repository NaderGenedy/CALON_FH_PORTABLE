#!/usr/bin/env python3
"""
56_apob_discordance_fh_vs_matched.py
======================================
ApoB/LDL-C discordance and ASCVD risk: FH carriers vs propensity-matched
non-FH controls from UK Biobank.

Discordance = ApoB/LDL ratio >= 0.31 (small, dense, particle-rich LDL).

Propensity matching: 1:5 nearest-neighbour on age, sex, BMI, smoking.
NOT matching on LDL/ApoB (exposure variables).

Author: Dr Nader Genedy
Date:   April 2026
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")
FIGURES = os.path.join(ANALYSIS, "figures")
os.makedirs(FIGURES, exist_ok=True)

DISC_THRESH = 0.31  # ApoB/LDL threshold for discordance
ASCVD_CODES = ['I20', 'I21', 'I22', 'I23', 'I24', 'I25',
               'I63', 'I64', 'I65', 'I66', 'I70', 'I71', 'I73', 'I74']

print("=" * 70)
print("ApoB/LDL DISCORDANCE x ASCVD: FH vs PROPENSITY-MATCHED NON-FH")
print("=" * 70)

# ============================================================================
# STEP 1: BUILD NON-FH UKB COHORT
# ============================================================================
print("\n[1/6] Building non-FH UKB cohort...")

# ApoB
apob = pd.read_csv("D:/calon_apob_CORRECT.csv")
apob.columns = [c.replace('participant.', '') for c in apob.columns]
apob['eid'] = apob['eid'].astype(str)
apob['apob'] = pd.to_numeric(apob['p30640_i0'], errors='coerce')
print(f"  ApoB: {apob['apob'].notna().sum():,}")

# LDL
ldl = pd.read_csv("D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_longitudinal_lipids.csv")
ldl.columns = [c.replace('participant.', '') for c in ldl.columns]
ldl['eid'] = ldl['eid'].astype(str)
ldl['ldl'] = pd.to_numeric(ldl['p30780_i0'], errors='coerce')
print(f"  LDL: {ldl['ldl'].notna().sum():,}")

# Sex
sex = pd.read_csv("D:/calon_sex.csv")
sex.columns = [c.replace('participant.', '') for c in sex.columns]
sex.rename(columns={'p31': 'sex'}, inplace=True)
sex['eid'] = sex['eid'].astype(str)

# Age + BMI
demo = pd.read_csv("D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_demographics.csv")
demo.columns = [c.replace('participant.', '') for c in demo.columns]
demo['eid'] = demo['eid'].astype(str)
demo['birth_year'] = pd.to_numeric(demo['p34'], errors='coerce')
demo['age'] = 2010 - demo['birth_year']
demo['bmi'] = pd.to_numeric(demo['p21001_i0'], errors='coerce')

# Smoking
smoke = pd.read_csv("D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_smoking.csv")
smoke.columns = [c.replace('participant.', '') for c in smoke.columns]
smoke['eid'] = smoke['eid'].astype(str)
smoke['smoking'] = pd.to_numeric(smoke['p20116_i0'], errors='coerce')
smoke['smoking_binary'] = (smoke['smoking'] == 2).astype(int)  # current smoker

# ASCVD from ICD-10
print("  Parsing ICD-10 for ASCVD...")
icd = pd.read_csv("D:/calon_backup_data/calon_batch_icd10.csv",
                   usecols=['participant.eid', 'participant.p41270'], dtype=str)
icd.columns = ['eid', 'icd10_raw']
# Fast string search for ASCVD codes
icd['ascvd'] = 0
for code in ASCVD_CODES:
    icd['ascvd'] = icd['ascvd'] | icd['icd10_raw'].str.contains(f'"{code}', na=False).astype(int)
icd['ascvd'] = icd['ascvd'].clip(0, 1)
print(f"  ASCVD events in UKB: {icd['ascvd'].sum():,} ({100*icd['ascvd'].mean():.1f}%)")

# FH carrier eids to exclude
carriers = pd.read_csv(os.path.join(ANALYSIS, "ukb_carriers_annotated.csv"),
                        usecols=['eid'], dtype={'eid': str})
fh_eids = set(carriers['eid'].unique())

# Merge all
ukb = apob[['eid', 'apob']].copy()
for df_merge, cols in [(ldl, ['eid', 'ldl']),
                        (sex, ['eid', 'sex']),
                        (demo, ['eid', 'age', 'bmi']),
                        (smoke, ['eid', 'smoking_binary']),
                        (icd, ['eid', 'ascvd'])]:
    ukb = ukb.merge(df_merge[cols], on='eid', how='left')

# Compute ApoB/LDL ratio
ukb['apob_ldl_ratio'] = np.where(
    ukb['ldl'].notna() & (ukb['ldl'] > 0) & ukb['apob'].notna(),
    ukb['apob'] / ukb['ldl'], np.nan)

# Exclude FH carriers
ukb_nonfh = ukb[~ukb['eid'].isin(fh_eids)].copy()

# Complete cases for matching
match_cols = ['apob_ldl_ratio', 'age', 'sex', 'bmi', 'smoking_binary', 'ascvd']
ukb_complete = ukb_nonfh.dropna(subset=match_cols).copy()
ukb_complete['discordant'] = (ukb_complete['apob_ldl_ratio'] >= DISC_THRESH).astype(int)

print(f"\n  Non-FH complete cases: {len(ukb_complete):,}")
print(f"  ApoB/LDL ratio: median={ukb_complete['apob_ldl_ratio'].median():.3f}")
print(f"  Discordant (>={DISC_THRESH}): {ukb_complete['discordant'].sum():,} "
      f"({100*ukb_complete['discordant'].mean():.1f}%)")
print(f"  ASCVD: {ukb_complete['ascvd'].sum():,} ({100*ukb_complete['ascvd'].mean():.1f}%)")

# ============================================================================
# STEP 2: BUILD FH COHORT
# ============================================================================
print("\n[2/6] Building FH cohort...")

fh = pd.read_csv("D:/CALON_FH_BACKUP_FULL/data/calon_ukb_analysis_ready.csv",
                  dtype={'eid': str})
# Standardise column names
fh_cols = {
    'apob_ldl_ratio': 'apob_ldl_ratio',
    'age': 'age',
    'sex': 'sex',
    'bmi': 'bmi',
    'smoking_binary': 'smoking_binary',
    'ascvd_combined': 'ascvd',
    'apob': 'apob',
    'ldl': 'ldl',
    'hypertension': 'hypertension',
    'diabetes': 'diabetes',
    'on_statin': 'on_statin',
}
fh_df = fh.rename(columns=fh_cols)[list(fh_cols.values()) + ['eid']].copy()
fh_df = fh_df.dropna(subset=['apob_ldl_ratio', 'age', 'sex', 'bmi', 'smoking_binary', 'ascvd'])
fh_df['discordant'] = (fh_df['apob_ldl_ratio'] >= DISC_THRESH).astype(int)

print(f"  FH complete cases: {len(fh_df):,}")
print(f"  ApoB/LDL ratio: median={fh_df['apob_ldl_ratio'].median():.3f}")
print(f"  Discordant (>={DISC_THRESH}): {fh_df['discordant'].sum():,} "
      f"({100*fh_df['discordant'].mean():.1f}%)")
print(f"  ASCVD: {fh_df['ascvd'].sum():,} ({100*fh_df['ascvd'].mean():.1f}%)")

# ============================================================================
# STEP 3: PROPENSITY SCORE MATCHING
# ============================================================================
print("\n[3/6] Propensity score matching (1:5)...")

# Combine for PS model
fh_df['is_fh'] = 1
ukb_complete['is_fh'] = 0

# PS model: P(FH) ~ age + sex + BMI + smoking
ps_features = ['age', 'sex', 'bmi', 'smoking_binary']

# Fit PS model
combined = pd.concat([fh_df[ps_features + ['is_fh', 'eid']],
                       ukb_complete[ps_features + ['is_fh', 'eid']]], ignore_index=True)
combined = combined.dropna(subset=ps_features)

X_ps = combined[ps_features].values
y_ps = combined['is_fh'].values

lr_ps = LogisticRegression(max_iter=1000, solver='lbfgs')
lr_ps.fit(X_ps, y_ps)
combined['ps'] = lr_ps.predict_proba(X_ps)[:, 1]

fh_ps = combined[combined['is_fh'] == 1].copy()
nonfh_ps = combined[combined['is_fh'] == 0].copy()

print(f"  FH propensity scores: median={fh_ps['ps'].median():.4f}, "
      f"range={fh_ps['ps'].min():.4f}-{fh_ps['ps'].max():.4f}")
print(f"  Non-FH PS: median={nonfh_ps['ps'].median():.4f}")

# 1:5 nearest neighbour matching with caliper
caliper = 0.2 * combined['ps'].std()
print(f"  Caliper: {caliper:.4f}")

nn = NearestNeighbors(n_neighbors=5, metric='euclidean')
nn.fit(nonfh_ps[['ps']].values)
distances, indices = nn.kneighbors(fh_ps[['ps']].values)

matched_nonfh_idx = set()
matched_pairs = []

for i, (dists, idxs) in enumerate(zip(distances, indices)):
    for d, j in zip(dists, idxs):
        if d <= caliper and j not in matched_nonfh_idx:
            matched_nonfh_idx.add(j)
            matched_pairs.append((fh_ps.index[i], nonfh_ps.index[int(j)]))

matched_nonfh_eids = set(nonfh_ps.iloc[list(matched_nonfh_idx)]['eid'])
matched_nonfh = ukb_complete[ukb_complete['eid'].isin(matched_nonfh_eids)].copy()

print(f"  Matched non-FH: {len(matched_nonfh):,} (ratio {len(matched_nonfh)/len(fh_df):.1f}:1)")

# Covariate balance
print(f"\n  Covariate balance (SMD):")
print(f"  {'Variable':>15s}  {'FH':>8s}  {'Matched':>8s}  {'Unmatched':>10s}  {'SMD_m':>6s}  {'SMD_u':>6s}")
print("  " + "-" * 65)
for var in ps_features:
    fh_mean = fh_df[var].mean()
    m_mean = matched_nonfh[var].mean()
    u_mean = ukb_complete[var].mean()
    fh_std = fh_df[var].std()
    pooled_m = np.sqrt((fh_df[var].std()**2 + matched_nonfh[var].std()**2) / 2)
    pooled_u = np.sqrt((fh_df[var].std()**2 + ukb_complete[var].std()**2) / 2)
    smd_m = abs(fh_mean - m_mean) / pooled_m if pooled_m > 0 else 0
    smd_u = abs(fh_mean - u_mean) / pooled_u if pooled_u > 0 else 0
    bal = "OK" if smd_m < 0.1 else "CHECK"
    print(f"  {var:>15s}  {fh_mean:>8.2f}  {m_mean:>8.2f}  {u_mean:>10.2f}  {smd_m:>6.3f}  {smd_u:>6.3f}  {bal}")

# ============================================================================
# STEP 4: DISCORDANCE x ASCVD ANALYSIS
# ============================================================================
print("\n[4/6] Discordance x ASCVD analysis...")

for label, df_anal in [("FH carriers", fh_df), ("Matched non-FH", matched_nonfh),
                        ("Unmatched non-FH", ukb_complete)]:
    disc = df_anal[df_anal['discordant'] == 1]
    conc = df_anal[df_anal['discordant'] == 0]

    r_disc = 100 * disc['ascvd'].mean() if len(disc) > 0 else 0
    r_conc = 100 * conc['ascvd'].mean() if len(conc) > 0 else 0

    # Odds ratio
    a = disc['ascvd'].sum()
    b = len(disc) - a
    c = conc['ascvd'].sum()
    d_val = len(conc) - c

    if min(a, b, c, d_val) > 0:
        odds, p_fish = stats.fisher_exact([[a, b], [c, d_val]])
        # 95% CI for OR
        log_or = np.log(odds)
        se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d_val)
        ci_lo = np.exp(log_or - 1.96 * se_log_or)
        ci_hi = np.exp(log_or + 1.96 * se_log_or)
        or_str = f"OR={odds:.2f} (95% CI {ci_lo:.2f}-{ci_hi:.2f}), P={p_fish:.4f}"
    else:
        or_str = "insufficient events"

    print(f"\n  {label} (n={len(df_anal):,}):")
    print(f"    Discordant (ratio>={DISC_THRESH}): n={len(disc):,} ({100*len(disc)/len(df_anal):.1f}%)")
    print(f"    ASCVD discordant:  {r_disc:.1f}%")
    print(f"    ASCVD concordant:  {r_conc:.1f}%")
    print(f"    {or_str}")

# Interaction test: discordance x FH status
print("\n  --- Interaction: Discordance x FH Status ---")
pooled = pd.concat([fh_df.assign(is_fh=1), matched_nonfh.assign(is_fh=0)], ignore_index=True)
pooled['disc_x_fh'] = pooled['discordant'] * pooled['is_fh']

X_int = pooled[['discordant', 'is_fh', 'disc_x_fh', 'age', 'sex']].values
y_int = pooled['ascvd'].values
lr_int = LogisticRegression(max_iter=1000, solver='lbfgs')
lr_int.fit(X_int, y_int)

feat_names = ['discordant', 'is_fh', 'disc_x_fh', 'age', 'sex']
print(f"  Logistic regression (FH + matched non-FH pooled):")
for feat, coef in zip(feat_names, lr_int.coef_[0]):
    print(f"    {feat:>15s}: OR={np.exp(coef):.3f} (beta={coef:.4f})")

# Continuous ApoB/LDL ratio as predictor
print("\n  --- Continuous ApoB/LDL ratio vs ASCVD ---")
for label, df_anal in [("FH", fh_df), ("Matched non-FH", matched_nonfh)]:
    sub = df_anal[['apob_ldl_ratio', 'ascvd']].dropna()
    if len(sub) < 50:
        continue
    X_c = sub[['apob_ldl_ratio']].values
    y_c = sub['ascvd'].values
    lr_c = LogisticRegression(max_iter=1000, solver='lbfgs')
    lr_c.fit(X_c, y_c)
    or_per_01 = np.exp(lr_c.coef_[0][0] * 0.1)  # OR per 0.1 unit increase
    # Point-biserial correlation
    rho, p = stats.pointbiserialr(y_c, sub['apob_ldl_ratio'].values)
    print(f"  {label}: OR per 0.1 ratio increase = {or_per_01:.3f}, "
          f"correlation rho={rho:.3f}, P={p:.4f}")

# ============================================================================
# STEP 5: DOMAIN-SPECIFIC DISCORDANCE (FH ONLY)
# ============================================================================
print("\n[5/6] Domain-specific discordance (FH carriers)...")

# Merge FH analysis-ready with carrier annotations for domain
fh_full = fh.copy()
fh_full['eid'] = fh_full['eid'].astype(str)
carriers_full = pd.read_csv(os.path.join(ANALYSIS, "ukb_carriers_annotated.csv"),
                             dtype={'eid': str})
carriers_dom = carriers_full.sort_values('sss', ascending=False).drop_duplicates('eid')
fh_dom = fh_full.merge(carriers_dom[['eid', 'ldlr_domain', 'mech_simple', 'sss', 'foldx_ddg']],
                        on='eid', how='inner')

if 'apob_ldl_ratio' in fh_dom.columns and 'ldlr_domain' in fh_dom.columns:
    dom_disc = fh_dom.groupby('ldlr_domain').agg(
        n=('eid', 'count'),
        mean_ratio=('apob_ldl_ratio', 'mean'),
        median_ratio=('apob_ldl_ratio', 'median'),
        pct_discordant=('apob_ldl_ratio', lambda x: 100 * (x >= DISC_THRESH).mean()),
        ascvd_rate=('ascvd_combined', lambda x: 100 * x.mean()),
    ).reset_index()
    dom_disc = dom_disc[dom_disc['n'] >= 10].sort_values('median_ratio', ascending=False)

    print(f"\n  {'Domain':>25s}  {'N':>5s}  {'Med Ratio':>10s}  {'Disc%':>6s}  {'ASCVD%':>7s}")
    print("  " + "-" * 60)
    for _, r in dom_disc.iterrows():
        print(f"  {r['ldlr_domain']:>25s}  {int(r['n']):>5d}  {r['median_ratio']:>10.3f}  "
              f"{r['pct_discordant']:>6.1f}  {r['ascvd_rate']:>7.1f}")

    dom_disc.to_csv(os.path.join(ANALYSIS, "apob_discordance_by_domain.csv"), index=False)

    # ASCVD by discordance within domains
    print(f"\n  ASCVD rate by discordance status within top domains:")
    for dom in dom_disc['ldlr_domain'].head(8):
        sub = fh_dom[fh_dom['ldlr_domain'] == dom]
        disc_sub = sub[sub['apob_ldl_ratio'] >= DISC_THRESH]
        conc_sub = sub[sub['apob_ldl_ratio'] < DISC_THRESH]
        if len(disc_sub) >= 5 and len(conc_sub) >= 5:
            r_d = 100 * disc_sub['ascvd_combined'].mean()
            r_c = 100 * conc_sub['ascvd_combined'].mean()
            short = dom.replace('Ligand-binding ', 'LB-').replace('EGF-like ', 'EGF-')
            print(f"    {short:>20s}: disc={r_d:.1f}% (n={len(disc_sub)}) "
                  f"vs conc={r_c:.1f}% (n={len(conc_sub)})")

# ============================================================================
# STEP 6: PUBLICATION FIGURE
# ============================================================================
print("\n[6/6] Generating publication figure...")

RED = '#b2182b'
BLUE = '#2166ac'
ORANGE = '#ef8a62'
LIGHT_BLUE = '#d1e5f0'

fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# Panel A: ApoB vs LDL scatter
ax = axes[0, 0]
# Sample non-FH for visibility
sample_nonfh = matched_nonfh.sample(min(5000, len(matched_nonfh)), random_state=42)
ax.scatter(sample_nonfh['ldl'], sample_nonfh['apob'], c=BLUE, alpha=0.15, s=5,
           label=f'Matched non-FH (n={len(matched_nonfh):,})', rasterized=True)
if 'ldl' in fh_df.columns:
    ax.scatter(fh_df['ldl'], fh_df['apob'], c=RED, alpha=0.5, s=15,
               label=f'FH carriers (n={len(fh_df):,})')
# Discordance line
x_line = np.linspace(0.5, 10, 100)
ax.plot(x_line, DISC_THRESH * x_line, '--', color='black', linewidth=1,
        label=f'Discordance line (ratio={DISC_THRESH})')
ax.set_xlabel('LDL-C (mmol/L)', fontsize=11, fontfamily='Arial')
ax.set_ylabel('ApoB (g/L)', fontsize=11, fontfamily='Arial')
ax.set_title('A. ApoB vs LDL-C', fontsize=12, fontweight='bold', fontfamily='Arial')
ax.legend(fontsize=7, loc='upper left', framealpha=0.8)
ax.set_xlim(0, 10)
ax.set_ylim(0, 2.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: ASCVD by discordance — grouped bar
ax = axes[0, 1]
groups = []
for label, df_b in [("FH", fh_df), ("Matched\nnon-FH", matched_nonfh)]:
    disc_r = 100 * df_b[df_b['discordant'] == 1]['ascvd'].mean()
    conc_r = 100 * df_b[df_b['discordant'] == 0]['ascvd'].mean()
    groups.append((label, conc_r, disc_r))

x_pos = np.arange(len(groups))
width = 0.35
conc_vals = [g[1] for g in groups]
disc_vals = [g[2] for g in groups]
bars1 = ax.bar(x_pos - width/2, conc_vals, width, label='Concordant', color=BLUE, edgecolor='black', linewidth=0.5)
bars2 = ax.bar(x_pos + width/2, disc_vals, width, label='Discordant', color=RED, edgecolor='black', linewidth=0.5)

ax.set_ylabel('ASCVD Rate (%)', fontsize=11, fontfamily='Arial')
ax.set_title('B. ASCVD by Discordance Status', fontsize=12, fontweight='bold', fontfamily='Arial')
ax.set_xticks(x_pos)
ax.set_xticklabels([g[0] for g in groups], fontsize=10)
ax.legend(fontsize=8)

# Add value labels
for bar_set in [bars1, bars2]:
    for bar in bar_set:
        h = bar.get_height()
        ax.annotate(f'{h:.1f}%', xy=(bar.get_x() + bar.get_width()/2, h),
                    xytext=(0, 3), textcoords='offset points', ha='center', fontsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: ApoB/LDL ratio distribution — FH vs non-FH
ax = axes[1, 0]
bins = np.linspace(0, 0.8, 40)
ax.hist(matched_nonfh['apob_ldl_ratio'].dropna(), bins=bins, density=True, alpha=0.6,
        color=BLUE, label=f'Matched non-FH', edgecolor='white', linewidth=0.3)
ax.hist(fh_df['apob_ldl_ratio'].dropna(), bins=bins, density=True, alpha=0.6,
        color=RED, label=f'FH carriers', edgecolor='white', linewidth=0.3)
ax.axvline(x=DISC_THRESH, color='black', linestyle='--', linewidth=1.5,
           label=f'Discordance threshold ({DISC_THRESH})')
ax.set_xlabel('ApoB/LDL-C Ratio', fontsize=11, fontfamily='Arial')
ax.set_ylabel('Density', fontsize=11, fontfamily='Arial')
ax.set_title('C. ApoB/LDL Ratio Distribution', fontsize=12, fontweight='bold', fontfamily='Arial')
ax.legend(fontsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Domain-specific discordance (FH)
ax = axes[1, 1]
if len(dom_disc) > 0:
    dd = dom_disc.sort_values('median_ratio')
    colours_d = [RED if r['pct_discordant'] > 50 else ORANGE if r['pct_discordant'] > 30 else BLUE
                 for _, r in dd.iterrows()]
    bars_d = ax.barh(range(len(dd)), dd['median_ratio'], color=colours_d,
                     edgecolor='black', linewidth=0.5)
    ax.axvline(x=DISC_THRESH, color='black', linestyle='--', linewidth=1)
    ax.set_yticks(range(len(dd)))
    ax.set_yticklabels([f"{r['ldlr_domain'].replace('Ligand-binding ','LB-')}\n"
                        f"(n={int(r['n'])}, ASCVD {r['ascvd_rate']:.0f}%)"
                        for _, r in dd.iterrows()], fontsize=7)
    ax.set_xlabel('Median ApoB/LDL-C Ratio', fontsize=11, fontfamily='Arial')
    ax.set_title('D. ApoB/LDL by LDLR Domain', fontsize=12, fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout(pad=2.0)
fig_path = os.path.join(FIGURES, "Figure_ApoB_Discordance_FH_vs_Matched.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {fig_path}")

# Save summary
summary = pd.DataFrame([
    {'group': 'FH', 'n': len(fh_df), 'disc_pct': 100*fh_df['discordant'].mean(),
     'ascvd_disc': 100*fh_df[fh_df['discordant']==1]['ascvd'].mean(),
     'ascvd_conc': 100*fh_df[fh_df['discordant']==0]['ascvd'].mean(),
     'median_ratio': fh_df['apob_ldl_ratio'].median()},
    {'group': 'Matched_nonFH', 'n': len(matched_nonfh), 'disc_pct': 100*matched_nonfh['discordant'].mean(),
     'ascvd_disc': 100*matched_nonfh[matched_nonfh['discordant']==1]['ascvd'].mean(),
     'ascvd_conc': 100*matched_nonfh[matched_nonfh['discordant']==0]['ascvd'].mean(),
     'median_ratio': matched_nonfh['apob_ldl_ratio'].median()},
])
summary.to_csv(os.path.join(ANALYSIS, "apob_discordance_fh_vs_matched.csv"), index=False)

print("\n" + "=" * 70)
print("COMPLETE")
print("=" * 70)

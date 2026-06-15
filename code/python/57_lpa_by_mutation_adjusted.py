#!/usr/bin/env python3
"""
57_lpa_by_mutation_adjusted.py
===============================
Does LDLR mutation location predict Lp(a) levels in FH?

Hypothesis: LDLR is a clearance pathway for Lp(a). If specific domains
handle Lp(a) catabolism preferentially, domain-specific mutations should
show differential Lp(a) elevation — even after adjusting for confounders.

Three cohorts pooled:
  1. DRAGON_3 (South Wales): 319 with mutation + Lp(a) nmol/L
  2. calon_ukb_analysis_ready (Wales FH): 1,308 with gene + Lp(a)
  3. UKB carriers: 2,370 with domain + Lp(a)

Adjusted for: age, sex, statin use, ezetimibe, PCSK9i, diabetes, smoking.

Author: Dr Nader Genedy
Date:   April 2026
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.linear_model import LinearRegression
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")
FIGURES = os.path.join(ANALYSIS, "figures")

print("=" * 70)
print("Lp(a) BY MUTATION — ADJUSTED ANALYSIS")
print("Does LDLR domain predict Lp(a) levels in FH?")
print("=" * 70)

# ============================================================================
# COHORT 1: DRAGON_3 (South Wales) — variant-level Lp(a)
# ============================================================================
print("\n[1/5] Loading DRAGON_3 (South Wales)...")

d3 = pd.read_csv("D:/CALON_AF3_PROJECT/data/fh_cohort/DRAGON_3.csv")
d3['lpa'] = pd.to_numeric(d3['Lpa'], errors='coerce')
d3['age'] = pd.to_numeric(d3['Currentage'], errors='coerce')
d3['sex'] = d3['Gender'].map({'M': 1, 'F': 0, 'Male': 1, 'Female': 0})
d3['smoking'] = pd.to_numeric(d3['Smoking_binary'], errors='coerce')
d3['diabetes'] = pd.to_numeric(d3['Diabetes_binary'], errors='coerce')
d3['on_statin'] = d3['Statin'].replace(' ', np.nan).notna().astype(int)
d3['on_ezetimibe'] = d3['Ezetimibe'].map({'Y': 1, 'N': 0, 'y': 1, 'n': 0})
d3['on_pcsk9i'] = d3['PCSK9i'].map({'Y': 1, 'N': 0, 'y': 1, 'n': 0})
d3['mutation'] = d3['Mutation1']
d3['cohort'] = 'South_Wales'

d3_anal = d3[d3['lpa'].notna() & d3['mutation'].notna()].copy()
print(f"  n={len(d3_anal)}, mutations={d3_anal['mutation'].nunique()}")
print(f"  Lp(a): median={d3_anal['lpa'].median():.1f}, IQR={d3_anal['lpa'].quantile(0.25):.0f}-{d3_anal['lpa'].quantile(0.75):.0f}")

# Map mutations to domains using the atlas
atlas = pd.read_csv(os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas_v3.csv"))
# Need VEP annotations to map cdna to protein position
vep = pd.read_csv(os.path.join(ANALYSIS, "ukb_ldlr_vep_annotations.csv"))

# Also use the SSS catalogue which has variant_id -> domain mapping
sss_cat = pd.read_csv(os.path.join(ANALYSIS, "sss_v2_catalogue.csv"))
sss_cols = [c for c in sss_cat.columns if c in ['variant_id', 'domain', 'sss', 'ddG', 'mechanism']]
if 'variant_id' in sss_cat.columns:
    sss_map = sss_cat.set_index('variant_id')[sss_cols[1:]].to_dict('index')
else:
    sss_map = {}

# Map domain from carrier annotations (these have cdna -> domain mapping)
carriers = pd.read_csv(os.path.join(ANALYSIS, "ukb_carriers_annotated.csv"), dtype={'eid': str})
# Create mutation -> domain lookup from carriers
mut_domain = carriers.dropna(subset=['ldlr_domain']).groupby('variant_id')['ldlr_domain'].first().to_dict()
# Also get cdna -> domain from the VEP + atlas pipeline
cdna_domain = carriers.dropna(subset=['ldlr_domain', 'cdna']).groupby('cdna')['ldlr_domain'].first().to_dict()

# Map DRAGON mutations to domains
d3_anal['domain'] = d3_anal['mutation'].map(mut_domain)
# Try cdna matching for unmapped
unmapped = d3_anal['domain'].isna()
if unmapped.sum() > 0:
    # Try matching without the gene prefix
    for idx in d3_anal[unmapped].index:
        mut = d3_anal.loc[idx, 'mutation']
        # Try exact match first
        if mut in cdna_domain:
            d3_anal.loc[idx, 'domain'] = cdna_domain[mut]
        else:
            # Try stripping LDLR: prefix and matching cdna
            cdna = mut.replace('LDLR:', '').replace('APOB:', '').replace('PCSK9:', '')
            if cdna in cdna_domain:
                d3_anal.loc[idx, 'domain'] = cdna_domain[cdna]

mapped = d3_anal['domain'].notna().sum()
print(f"  Mapped to domain: {mapped}/{len(d3_anal)} ({100*mapped/len(d3_anal):.0f}%)")

# ============================================================================
# COHORT 2: UKB carriers — domain-level Lp(a)
# ============================================================================
print("\n[2/5] Loading UKB carriers...")

lpa_ukb = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/calon_lpa_CORRECT.csv")
lpa_ukb.columns = [c.replace('participant.', '') for c in lpa_ukb.columns]
lpa_ukb['eid'] = lpa_ukb['eid'].astype(str)
lpa_ukb['lpa'] = pd.to_numeric(lpa_ukb['p30790_i0'], errors='coerce')

sex_df = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/calon_sex.csv")
sex_df.columns = [c.replace('participant.', '') for c in sex_df.columns]
sex_df.rename(columns={'p31': 'sex'}, inplace=True)
sex_df['eid'] = sex_df['eid'].astype(str)

demo = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/ukb_reviewer_demographics.csv")
demo.columns = [c.replace('participant.', '') for c in demo.columns]
demo['eid'] = demo['eid'].astype(str)
demo['age'] = 2010 - pd.to_numeric(demo['p34'], errors='coerce')

patients = carriers.sort_values('sss', ascending=False).drop_duplicates('eid')
ukb = patients.merge(lpa_ukb[['eid', 'lpa']], on='eid', how='inner')
ukb = ukb.merge(sex_df[['eid', 'sex']], on='eid', how='left')
ukb = ukb.merge(demo[['eid', 'age']], on='eid', how='left')
ukb = ukb[ukb['lpa'].notna()].copy()
ukb['mutation'] = ukb['variant_id']
ukb['domain'] = ukb['ldlr_domain']
ukb['cohort'] = 'UKB'

print(f"  n={len(ukb)}, with domain={ukb['domain'].notna().sum()}")
print(f"  Lp(a): median={ukb['lpa'].median():.1f}")

# ============================================================================
# COHORT 3: analysis_ready (Wales wider) — gene-level only
# ============================================================================
print("\n[3/5] Loading analysis_ready (Wales)...")

ar = pd.read_csv("D:/CALON_AF3_PROJECT/data/fh_cohort/calon_ukb_analysis_ready.csv")
ar['lpa'] = pd.to_numeric(ar['lpa'], errors='coerce')
ar_with = ar[ar['lpa'].notna()].copy()
ar_with['cohort'] = 'Wales_AR'
print(f"  n={len(ar_with)}, genes={ar_with['gene'].value_counts().to_dict()}")
print(f"  Lp(a): median={ar_with['lpa'].median():.1f}")

# ============================================================================
# ANALYSIS A: Lp(a) BY GENE (all cohorts)
# ============================================================================
print("\n" + "=" * 70)
print("[4/5] ANALYSIS A: Lp(a) BY GENE")
print("=" * 70)

# Gene-level from analysis_ready (largest with gene annotation)
print(f"\n  Lp(a) by gene (analysis_ready, n={len(ar_with)}):")
for gene in ['LDLR', 'APOB', 'PCSK9']:
    sub = ar_with[ar_with['gene'] == gene]['lpa']
    if len(sub) >= 5:
        print(f"    {gene}: n={len(sub)}, median={sub.median():.1f}, "
              f"IQR={sub.quantile(0.25):.0f}-{sub.quantile(0.75):.0f}")

# LDLR vs APOB comparison
ldlr_lpa = ar_with[ar_with['gene'] == 'LDLR']['lpa']
apob_lpa = ar_with[ar_with['gene'] == 'APOB']['lpa']
if len(ldlr_lpa) >= 10 and len(apob_lpa) >= 10:
    u, p = stats.mannwhitneyu(ldlr_lpa, apob_lpa, alternative='two-sided')
    d = (ldlr_lpa.mean() - apob_lpa.mean()) / np.sqrt((ldlr_lpa.std()**2 + apob_lpa.std()**2) / 2)
    print(f"\n    LDLR vs APOB: median {ldlr_lpa.median():.1f} vs {apob_lpa.median():.1f}, "
          f"P={p:.4f}, d={d:.3f}")

# ============================================================================
# ANALYSIS B: Lp(a) BY LDLR DOMAIN (UKB carriers — largest with domain)
# ============================================================================
print(f"\n{'='*70}")
print("ANALYSIS B: Lp(a) BY LDLR DOMAIN (UKB carriers)")
print("=" * 70)

ukb_dom = ukb[ukb['domain'].notna()].copy()
dom_lpa = ukb_dom.groupby('domain').agg(
    n=('lpa', 'count'),
    median_lpa=('lpa', 'median'),
    mean_lpa=('lpa', 'mean'),
    std_lpa=('lpa', 'std'),
    iqr_lo=('lpa', lambda x: x.quantile(0.25)),
    iqr_hi=('lpa', lambda x: x.quantile(0.75)),
).reset_index()
dom_lpa = dom_lpa[dom_lpa['n'] >= 10].sort_values('median_lpa', ascending=False)

# Add penetrance from atlas
pen_map = atlas.groupby('domain')['domain_penetrance_by_60_pct'].first().to_dict()
dom_lpa['penetrance'] = dom_lpa['domain'].map(pen_map)

print(f"\n  {'Domain':>25s}  {'N':>5s}  {'Median':>7s}  {'IQR':>12s}  {'Mean':>7s}  {'Pen%':>5s}")
print("  " + "-" * 70)
for _, r in dom_lpa.iterrows():
    pen = f"{r['penetrance']:.1f}" if pd.notna(r['penetrance']) else 'N/A'
    print(f"  {r['domain']:>25s}  {int(r['n']):>5d}  {r['median_lpa']:>7.1f}  "
          f"{r['iqr_lo']:>5.0f}-{r['iqr_hi']:>4.0f}  {r['mean_lpa']:>7.1f}  {pen:>5s}")

# Kruskal-Wallis across domains
groups = [g['lpa'].values for _, g in ukb_dom.groupby('domain') if len(g) >= 10]
if len(groups) >= 3:
    h_stat, p_kw = stats.kruskal(*groups)
    print(f"\n  Kruskal-Wallis across domains: H={h_stat:.2f}, P={p_kw:.4f}")

# Spearman: penetrance vs median Lp(a)
valid = dom_lpa.dropna(subset=['penetrance', 'median_lpa'])
if len(valid) >= 5:
    rho, p_sp = stats.spearmanr(valid['penetrance'], valid['median_lpa'])
    print(f"  Penetrance vs median Lp(a): rho={rho:.3f}, P={p_sp:.4f} (n={len(valid)} domains)")

# ============================================================================
# ANALYSIS C: Lp(a) BY MUTATION TYPE (consequence)
# ============================================================================
print(f"\n{'='*70}")
print("ANALYSIS C: Lp(a) BY MUTATION CONSEQUENCE (UKB)")
print("=" * 70)

ukb_cons = ukb[ukb['consequence'].notna()].copy()
# Simplify consequence
ukb_cons['cons_simple'] = 'Other'
ukb_cons.loc[ukb_cons['consequence'].str.contains('missense', na=False), 'cons_simple'] = 'Missense'
ukb_cons.loc[ukb_cons['consequence'].str.contains('frameshift|stop_gained', na=False), 'cons_simple'] = 'LoF (frameshift/stop)'
ukb_cons.loc[ukb_cons['consequence'].str.contains('splice', na=False), 'cons_simple'] = 'Splice'
ukb_cons.loc[ukb_cons['consequence'].str.contains('synonymous', na=False), 'cons_simple'] = 'Synonymous'

for cons in ['Missense', 'LoF (frameshift/stop)', 'Splice', 'Synonymous', 'Other']:
    sub = ukb_cons[ukb_cons['cons_simple'] == cons]['lpa']
    if len(sub) >= 10:
        print(f"  {cons:>25s}: n={len(sub)}, median={sub.median():.1f}, "
              f"IQR={sub.quantile(0.25):.0f}-{sub.quantile(0.75):.0f}")

# ============================================================================
# ANALYSIS D: Lp(a) BY MECHANISM TYPE
# ============================================================================
print(f"\n{'='*70}")
print("ANALYSIS D: Lp(a) BY MECHANISM TYPE (UKB carriers)")
print("=" * 70)

for mech in ['Type 1: Structural', 'Type 2: Functional', 'Benign', 'Unclassified']:
    sub = ukb[ukb['mech_simple'] == mech]['lpa']
    if len(sub) >= 10:
        print(f"  {mech:>25s}: n={len(sub)}, median={sub.median():.1f}, "
              f"mean={sub.mean():.1f}")

# ============================================================================
# ANALYSIS E: ADJUSTED REGRESSION (DRAGON_3 — has most covariates)
# ============================================================================
print(f"\n{'='*70}")
print("ANALYSIS E: ADJUSTED REGRESSION — Lp(a) ~ domain + covariates")
print("=" * 70)

# Use DRAGON_3 for adjusted analysis (has clinical covariates)
d3_reg = d3_anal[d3_anal['domain'].notna()].copy()
d3_reg['on_ezetimibe'] = d3_reg['on_ezetimibe'].fillna(0)
d3_reg['on_pcsk9i'] = d3_reg['on_pcsk9i'].fillna(0)
d3_reg['smoking'] = d3_reg['smoking'].fillna(0)
d3_reg['diabetes'] = d3_reg['diabetes'].fillna(0)
d3_reg['log_lpa'] = np.log1p(d3_reg['lpa'])  # log-transform for regression

print(f"  DRAGON_3 with domain + Lp(a) + covariates: {len(d3_reg)}")

if len(d3_reg) >= 30:
    # Unadjusted: domain only
    covariates_base = ['age', 'sex']
    covariates_full = ['age', 'sex', 'on_statin', 'on_ezetimibe', 'on_pcsk9i',
                       'smoking', 'diabetes']

    for model_name, covs in [('Unadjusted', []),
                              ('Age + sex', covariates_base),
                              ('Fully adjusted', covariates_full)]:
        reg_data = d3_reg.dropna(subset=['log_lpa'] + covs)
        if len(reg_data) < 20:
            print(f"\n  {model_name}: insufficient data (n={len(reg_data)})")
            continue

        # Create domain dummies
        domain_dummies = pd.get_dummies(reg_data['domain'], prefix='dom', drop_first=True)
        X_cols = list(domain_dummies.columns) + covs
        X = pd.concat([domain_dummies, reg_data[covs].reset_index(drop=True)], axis=1)
        y = reg_data['log_lpa'].reset_index(drop=True)

        # Handle NaN in X
        valid = X.notna().all(axis=1) & y.notna()
        X_clean = X[valid].values
        y_clean = y[valid].values

        if len(y_clean) >= 20:
            lr = LinearRegression()
            lr.fit(X_clean, y_clean)
            r2 = lr.score(X_clean, y_clean)
            print(f"\n  {model_name} (n={len(y_clean)}, R2={r2:.3f}):")

            # Show domain coefficients
            for feat, coef in zip(X_cols, lr.coef_):
                if feat.startswith('dom_'):
                    dom_name = feat.replace('dom_', '')
                    # Exponentiate to get multiplicative effect on Lp(a)
                    pct_change = (np.exp(coef) - 1) * 100
                    direction = "+" if pct_change > 0 else ""
                    print(f"    {dom_name:>25s}: {direction}{pct_change:.1f}% vs reference")

            # Show covariate effects
            for feat, coef in zip(X_cols, lr.coef_):
                if not feat.startswith('dom_'):
                    pct_change = (np.exp(coef) - 1) * 100
                    print(f"    {feat:>25s}: beta={coef:.4f} ({pct_change:+.1f}% per unit)")

# ============================================================================
# ANALYSIS F: TOP VARIANT-LEVEL Lp(a) (DRAGON_3 — specific mutations)
# ============================================================================
print(f"\n{'='*70}")
print("ANALYSIS F: VARIANT-LEVEL Lp(a) (mutations with n>=3)")
print("=" * 70)

var_lpa = d3_anal.groupby('mutation').agg(
    n=('lpa', 'count'),
    median_lpa=('lpa', 'median'),
    mean_lpa=('lpa', 'mean'),
    min_lpa=('lpa', 'min'),
    max_lpa=('lpa', 'max'),
).reset_index()
var_lpa = var_lpa[var_lpa['n'] >= 3].sort_values('median_lpa', ascending=False)

# Add domain
var_lpa['domain'] = var_lpa['mutation'].map(mut_domain)
# Fill from cdna
for idx in var_lpa.index:
    if pd.isna(var_lpa.loc[idx, 'domain']):
        cdna = var_lpa.loc[idx, 'mutation'].replace('LDLR:', '').replace('APOB:', '')
        if cdna in cdna_domain:
            var_lpa.loc[idx, 'domain'] = cdna_domain[cdna]

print(f"\n  {'Mutation':>40s}  {'N':>4s}  {'Median':>7s}  {'Range':>12s}  {'Domain':>20s}")
print("  " + "-" * 90)
for _, r in var_lpa.head(20).iterrows():
    dom = r['domain'] if pd.notna(r['domain']) else '?'
    print(f"  {r['mutation']:>40s}  {int(r['n']):>4d}  {r['median_lpa']:>7.0f}  "
          f"{r['min_lpa']:>5.0f}-{r['max_lpa']:>4.0f}  {dom:>20s}")

# Overall population median for reference
pop_med = d3_anal['lpa'].median()
print(f"\n  Population median: {pop_med:.0f} nmol/L")
print(f"  Mutations above median: {(var_lpa['median_lpa'] > pop_med).sum()}/{len(var_lpa)}")

# ============================================================================
# FIGURE: 4-panel
# ============================================================================
print(f"\n{'='*70}")
print("[5/5] Generating figure...")

RED = '#b2182b'
BLUE = '#2166ac'
ORANGE = '#ef8a62'

fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# Panel A: Lp(a) by LDLR domain (UKB, boxplot)
ax = axes[0, 0]
dom_order = dom_lpa.sort_values('median_lpa')['domain'].tolist()
box_data = []
box_labels = []
for dom in dom_order:
    vals = ukb_dom[ukb_dom['domain'] == dom]['lpa'].dropna().values
    if len(vals) >= 10:
        box_data.append(vals)
        short = dom.replace('Ligand-binding ', 'LB-').replace('EGF-like ', 'EGF-')
        box_labels.append(f"{short}\n(n={len(vals)})")

if box_data:
    bp = ax.boxplot(box_data, labels=box_labels, patch_artist=True, vert=False,
                    medianprops=dict(color='black', linewidth=2), showfliers=False)
    for patch in bp['boxes']:
        patch.set_facecolor(BLUE)
        patch.set_alpha(0.7)
    ax.set_xlabel('Lp(a) (nmol/L)', fontsize=10, fontfamily='Arial')
    ax.set_title('A. Lp(a) by LDLR Domain (UKB)', fontsize=12, fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: Lp(a) by gene (analysis_ready)
ax = axes[0, 1]
gene_data = []
gene_labels = []
for gene in ['LDLR', 'APOB']:
    vals = ar_with[ar_with['gene'] == gene]['lpa'].dropna().values
    if len(vals) >= 10:
        gene_data.append(vals)
        gene_labels.append(f"{gene}\n(n={len(vals)})")

if gene_data:
    bp2 = ax.boxplot(gene_data, labels=gene_labels, patch_artist=True,
                     medianprops=dict(color='black', linewidth=2), showfliers=False)
    cols = [RED, ORANGE]
    for i, patch in enumerate(bp2['boxes']):
        patch.set_facecolor(cols[i])
        patch.set_alpha(0.7)
    ax.set_ylabel('Lp(a) (nmol/L)', fontsize=10, fontfamily='Arial')
    ax.set_title('B. Lp(a) by Gene', fontsize=12, fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: Top mutations by Lp(a) (DRAGON_3)
ax = axes[1, 0]
top_muts = var_lpa.head(12).sort_values('median_lpa')
if len(top_muts) > 0:
    colours = [RED if r['median_lpa'] > pop_med else BLUE for _, r in top_muts.iterrows()]
    ax.barh(range(len(top_muts)), top_muts['median_lpa'], color=colours,
            edgecolor='black', linewidth=0.5)
    ax.axvline(x=pop_med, color='black', linestyle='--', linewidth=1,
               label=f'Cohort median ({pop_med:.0f})')
    ax.set_yticks(range(len(top_muts)))
    ax.set_yticklabels([f"{r['mutation'][:30]}\n(n={int(r['n'])})"
                        for _, r in top_muts.iterrows()], fontsize=6)
    ax.set_xlabel('Median Lp(a) (nmol/L)', fontsize=10, fontfamily='Arial')
    ax.set_title('C. Lp(a) by Mutation (South Wales)', fontsize=12,
                 fontweight='bold', fontfamily='Arial')
    ax.legend(fontsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Domain penetrance vs Lp(a)
ax = axes[1, 1]
valid_corr = dom_lpa.dropna(subset=['penetrance', 'median_lpa'])
if len(valid_corr) >= 4:
    ax.scatter(valid_corr['penetrance'], valid_corr['median_lpa'],
               s=valid_corr['n'] * 2, c=RED, alpha=0.7, edgecolors='black')
    for _, r in valid_corr.iterrows():
        short = r['domain'].replace('Ligand-binding ', 'LB-').replace('EGF-like ', 'EGF-')
        ax.annotate(short, (r['penetrance'], r['median_lpa']),
                    fontsize=6, ha='center', va='bottom')
    rho, p = stats.spearmanr(valid_corr['penetrance'], valid_corr['median_lpa'])
    ax.set_xlabel('ASCVD Penetrance by Age 60 (%)', fontsize=10, fontfamily='Arial')
    ax.set_ylabel('Median Lp(a) (nmol/L)', fontsize=10, fontfamily='Arial')
    ax.set_title(f'D. Penetrance vs Lp(a)\n(rho={rho:.3f}, P={p:.3f})',
                 fontsize=12, fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout(pad=2.0)
fig_path = os.path.join(FIGURES, "Figure_Lpa_by_Mutation.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {fig_path}")

# Save outputs
dom_lpa.to_csv(os.path.join(ANALYSIS, "lpa_by_domain.csv"), index=False)
var_lpa.to_csv(os.path.join(ANALYSIS, "lpa_by_variant.csv"), index=False)

print("\n" + "=" * 70)
print("COMPLETE")
print("=" * 70)

#!/usr/bin/env python3
"""
58_lpa_adjusted_mutation_level.py
==================================
Fully adjusted Lp(a) analysis by LDLR domain and mutation,
with parallel untreated LDL comparison.

Key question: Does LDLR domain predict Lp(a) INDEPENDENTLY of
age, sex, ancestry, treatment? Does Lp(a) track or diverge from LDL?

Author: Dr Nader Genedy
Date:   April 2026
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.linear_model import LinearRegression
import re
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
print("ADJUSTED Lp(a) BY LDLR MUTATION + UNTREATED LDL COMPARISON")
print("=" * 70)

# ============================================================================
# STEP 1: cDNA -> PROTEIN POSITION -> DOMAIN MAPPING
# ============================================================================
print("\n[1/7] Building cDNA -> domain mapping...")

# Load domain atlas
atlas = pd.read_csv(os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas_v3.csv"))
pos_to_domain = atlas.set_index('position')['domain'].to_dict()

def cdna_to_protein_pos(mutation_str):
    """Convert HGVS cDNA notation to protein position.
    LDLR coding sequence starts at c.1 = codon 1 (after signal peptide cleavage).
    Protein position = (cDNA_position - 1) // 3 + 1
    """
    if pd.isna(mutation_str):
        return None, 'missing'

    mut = str(mutation_str).strip()

    # Exon-level deletions
    if 'Deletion of exon' in mut or 'deletion of exon' in mut:
        # Map exon ranges to approximate protein positions
        exon_ranges = {
            1: (1, 21), 2: (22, 58), 3: (59, 100), 4: (101, 167),
            5: (168, 228), 6: (229, 285), 7: (286, 352), 8: (353, 393),
            9: (394, 445), 10: (446, 526), 11: (527, 592), 12: (593, 635),
            13: (636, 693), 14: (694, 749), 15: (750, 770), 16: (771, 800),
            17: (801, 830), 18: (831, 860)
        }
        # Extract exon number(s)
        nums = re.findall(r'(\d+)', mut.replace('LDLR', ''))
        if nums:
            exon = int(nums[0])
            if exon in exon_ranges:
                mid = (exon_ranges[exon][0] + exon_ranges[exon][1]) // 2
                return mid, 'exon_deletion'
        return None, 'exon_unmapped'

    # Promoter variants
    if 'promoter' in mut.lower() or mut.startswith('LDLR:c.-'):
        return None, 'promoter'

    # Intronic / splice
    m_splice = re.search(r'c\.(\d+)[+-]', mut)
    if m_splice:
        cdna_pos = int(m_splice.group(1))
        prot_pos = (cdna_pos - 1) // 3 + 1
        return prot_pos, 'splice'

    # Standard substitution: c.1217G>C
    m_sub = re.search(r'c\.(\d+)[ACGT]>[ACGT]', mut)
    if m_sub:
        cdna_pos = int(m_sub.group(1))
        prot_pos = (cdna_pos - 1) // 3 + 1
        return prot_pos, 'substitution'

    # Deletion: c.1216delC or c.654_656delTGG
    m_del = re.search(r'c\.(\d+)', mut)
    if m_del and 'del' in mut.lower():
        cdna_pos = int(m_del.group(1))
        prot_pos = (cdna_pos - 1) // 3 + 1
        return prot_pos, 'deletion'

    # Duplication: c.2569_2582dup
    m_dup = re.search(r'c\.(\d+)', mut)
    if m_dup and 'dup' in mut.lower():
        cdna_pos = int(m_dup.group(1))
        prot_pos = (cdna_pos - 1) // 3 + 1
        return prot_pos, 'duplication'

    # Insertion: c.16_17insTTCCT
    m_ins = re.search(r'c\.(\d+)', mut)
    if m_ins and 'ins' in mut.lower():
        cdna_pos = int(m_ins.group(1))
        prot_pos = (cdna_pos - 1) // 3 + 1
        return prot_pos, 'insertion'

    # APOB / PCSK9 variants — no domain mapping
    if 'APOB' in mut or 'PCSK9' in mut:
        return None, 'non_LDLR'

    # Fallback: try to extract any number
    m_any = re.search(r'c\.(\d+)', mut)
    if m_any:
        cdna_pos = int(m_any.group(1))
        prot_pos = (cdna_pos - 1) // 3 + 1
        return prot_pos, 'fallback'

    return None, 'unparseable'


# Apply to DRAGON_3
d3 = pd.read_csv("D:/CALON_AF3_PROJECT/data/fh_cohort/DRAGON_3.csv")
d3['lpa'] = pd.to_numeric(d3['Lpa'], errors='coerce')
d3['ldl_untreated'] = pd.to_numeric(d3['LDL_1'], errors='coerce')
d3['age'] = pd.to_numeric(d3['Currentage'], errors='coerce')
d3['sex'] = d3['Gender'].map({'M': 1, 'F': 0, 'Male': 1, 'Female': 0})
d3['smoking'] = pd.to_numeric(d3['Smoking_binary'], errors='coerce').fillna(0)
d3['diabetes'] = pd.to_numeric(d3['Diabetes_binary'], errors='coerce').fillna(0)
d3['on_statin'] = d3['Statin'].replace(' ', np.nan).notna().astype(int)
d3['on_ezetimibe'] = d3['Ezetimibe'].map({'Y': 1, 'N': 0, 'y': 1, 'n': 0}).fillna(0)
d3['on_pcsk9i'] = d3['PCSK9i'].map({'Y': 1, 'N': 0, 'y': 1, 'n': 0}).fillna(0)

# Map mutations to domain
results = d3['Mutation1'].apply(cdna_to_protein_pos)
d3['protein_pos'] = [r[0] for r in results]
d3['map_type'] = [r[1] for r in results]
d3['domain'] = d3['protein_pos'].map(pos_to_domain)

# Report mapping
map_summary = d3[d3['Mutation1'].notna()]['map_type'].value_counts()
print(f"  DRAGON_3 mapping results:")
for mt, n in map_summary.items():
    print(f"    {mt:>20s}: {n}")

d3_mapped = d3[d3['domain'].notna() & d3['lpa'].notna()]
print(f"\n  Successfully mapped with Lp(a): {len(d3_mapped)}/{d3['lpa'].notna().sum()} "
      f"({100*len(d3_mapped)/max(d3['lpa'].notna().sum(),1):.0f}%)")

# ============================================================================
# STEP 2: MERGE UKB CARRIERS WITH ALL COVARIATES
# ============================================================================
print("\n[2/7] Merging UKB carriers with covariates...")

carriers = pd.read_csv(os.path.join(ANALYSIS, "ukb_carriers_annotated.csv"), dtype={'eid': str})
patients = carriers.sort_values('sss', ascending=False).drop_duplicates('eid')

# Lp(a)
lpa_ukb = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/calon_lpa_CORRECT.csv")
lpa_ukb.columns = [c.replace('participant.', '') for c in lpa_ukb.columns]
lpa_ukb['eid'] = lpa_ukb['eid'].astype(str)
lpa_ukb['lpa'] = pd.to_numeric(lpa_ukb['p30790_i0'], errors='coerce')

# Sex
sex_df = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/calon_sex.csv")
sex_df.columns = [c.replace('participant.', '') for c in sex_df.columns]
sex_df.rename(columns={'p31': 'sex'}, inplace=True)
sex_df['eid'] = sex_df['eid'].astype(str)

# Age + ethnicity
demo = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/ukb_reviewer_demographics.csv")
demo.columns = [c.replace('participant.', '') for c in demo.columns]
demo['eid'] = demo['eid'].astype(str)
demo['age'] = 2010 - pd.to_numeric(demo['p34'], errors='coerce')
demo['ethnicity_raw'] = pd.to_numeric(demo['p21000_i0'], errors='coerce')
# Group ethnicity
demo['ethnicity'] = 'Other'
demo.loc[demo['ethnicity_raw'].isin([1001, 1002, 1003]), 'ethnicity'] = 'White'
demo.loc[demo['ethnicity_raw'].isin([3001, 3002, 3003, 3004]), 'ethnicity'] = 'Asian'
demo.loc[demo['ethnicity_raw'].isin([4001, 4002, 4003]), 'ethnicity'] = 'Black'
demo.loc[demo['ethnicity_raw'].isin([2001, 2002, 2003, 2004]), 'ethnicity'] = 'Mixed'

# Genetic PCs (first 10)
print("  Loading genetic PCs...")
pc_cols_to_load = ['participant.eid'] + [f'participant.p22009_a{i}' for i in range(1, 11)]
pcs = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/calon_batch_genetic_pcs.csv",
                   usecols=pc_cols_to_load)
pcs.columns = [c.replace('participant.', '') for c in pcs.columns]
pcs['eid'] = pcs['eid'].astype(str)
pc_names = [f'p22009_a{i}' for i in range(1, 11)]

# LDL (for untreated comparison)
ldl_df = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/ukb_reviewer_longitudinal_lipids.csv",
                      usecols=['participant.eid', 'participant.p30780_i0'])
ldl_df.columns = [c.replace('participant.', '') for c in ldl_df.columns]
ldl_df['eid'] = ldl_df['eid'].astype(str)
ldl_df['ldl'] = pd.to_numeric(ldl_df['p30780_i0'], errors='coerce')

# Merge all
ukb = patients[['eid', 'ldlr_domain', 'mech_simple', 'foldx_ddg', 'variant_id',
                 'consequence', 'protein_pos']].copy()
ukb = ukb.rename(columns={'ldlr_domain': 'domain'})
for df_m, cols in [(lpa_ukb, ['eid', 'lpa']), (sex_df, ['eid', 'sex']),
                    (demo, ['eid', 'age', 'ethnicity']), (pcs, ['eid'] + pc_names),
                    (ldl_df, ['eid', 'ldl'])]:
    ukb = ukb.merge(df_m[cols], on='eid', how='left')

ukb = ukb[ukb['lpa'].notna()].copy()
ukb['log_lpa'] = np.log1p(ukb['lpa'])
ukb['cohort'] = 'UKB'

print(f"  UKB carriers with Lp(a): {len(ukb)}")
print(f"    With domain: {ukb['domain'].notna().sum()}")
print(f"    With ethnicity: {ukb['ethnicity'].notna().sum()}")
print(f"    With PCs: {ukb[pc_names[0]].notna().sum()}")
print(f"    With LDL: {ukb['ldl'].notna().sum()}")

# ============================================================================
# STEP 3: ADJUSTED REGRESSION HIERARCHY
# ============================================================================
print(f"\n{'='*70}")
print("[3/7] ADJUSTED REGRESSION: log(Lp(a)) ~ domain + covariates")
print("=" * 70)

# UKB analysis (has PCs + ethnicity)
ukb_reg = ukb[ukb['domain'].notna()].copy()

# Ethnicity dummies
eth_dummies = pd.get_dummies(ukb_reg['ethnicity'], prefix='eth', drop_first=True)
# Domain dummies
dom_dummies = pd.get_dummies(ukb_reg['domain'], prefix='dom', drop_first=True)
dom_cols = list(dom_dummies.columns)

models = [
    ('1. Unadjusted (domain only)', dom_cols, dom_dummies),
    ('2. + Age + Sex', dom_cols + ['age', 'sex'],
     pd.concat([dom_dummies, ukb_reg[['age', 'sex']].reset_index(drop=True)], axis=1)),
    ('3. + Ethnicity', dom_cols + ['age', 'sex'] + list(eth_dummies.columns),
     pd.concat([dom_dummies, ukb_reg[['age', 'sex']].reset_index(drop=True),
                eth_dummies.reset_index(drop=True)], axis=1)),
    ('4. + Genetic PCs (1-10)', dom_cols + ['age', 'sex'] + list(eth_dummies.columns) + pc_names,
     pd.concat([dom_dummies, ukb_reg[['age', 'sex']].reset_index(drop=True),
                eth_dummies.reset_index(drop=True),
                ukb_reg[pc_names].reset_index(drop=True)], axis=1)),
]

y = ukb_reg['log_lpa'].reset_index(drop=True)
model_results = []

for model_name, feat_names, X_df in models:
    X = X_df.reset_index(drop=True)
    valid = X.notna().all(axis=1) & y.notna()
    X_clean = X[valid].values.astype(float)
    y_clean = y[valid].values

    if len(y_clean) < 50:
        print(f"\n  {model_name}: insufficient data (n={len(y_clean)})")
        continue

    lr = LinearRegression()
    lr.fit(X_clean, y_clean)
    r2 = lr.score(X_clean, y_clean)
    n_obs = len(y_clean)

    # F-test for domain effect (compare to model without domain)
    non_dom_idx = [i for i, f in enumerate(feat_names) if not f.startswith('dom_')]
    if non_dom_idx:
        X_reduced = X_clean[:, non_dom_idx]
        lr_reduced = LinearRegression()
        lr_reduced.fit(X_reduced, y_clean)
        r2_reduced = lr_reduced.score(X_reduced, y_clean)
        # Partial F-test
        k_full = X_clean.shape[1]
        k_reduced = X_reduced.shape[1]
        n = len(y_clean)
        if r2 > r2_reduced and (1 - r2) > 0:
            f_stat = ((r2 - r2_reduced) / (k_full - k_reduced)) / ((1 - r2) / (n - k_full - 1))
            p_domain = 1 - stats.f.cdf(f_stat, k_full - k_reduced, n - k_full - 1)
        else:
            f_stat = 0
            p_domain = 1.0
    else:
        r2_reduced = 0
        p_domain = np.nan
        f_stat = 0

    model_results.append({'model': model_name, 'n': n_obs, 'r2': r2,
                          'r2_domain_only': r2 - r2_reduced, 'p_domain': p_domain})

    print(f"\n  {model_name} (n={n_obs}):")
    print(f"    R² = {r2:.4f}, Domain R² contribution = {r2 - r2_reduced:.4f}")
    if not np.isnan(p_domain):
        print(f"    F-test for domain effect: F={f_stat:.2f}, P={p_domain:.4f}")

    # Show domain coefficients (% change in Lp(a) vs reference)
    print(f"    Domain effects (% change vs reference):")
    for feat, coef in zip(feat_names, lr.coef_):
        if feat.startswith('dom_'):
            dom_name = feat.replace('dom_', '')
            pct = (np.exp(coef) - 1) * 100
            print(f"      {dom_name:>25s}: {pct:>+6.1f}%")

# DRAGON_3 adjusted (has treatment covariates)
print(f"\n  --- DRAGON_3 Adjusted (has treatment data) ---")
d3_reg = d3_mapped.copy()
if len(d3_reg) >= 30:
    d3_reg['log_lpa'] = np.log1p(d3_reg['lpa'])
    d3_dom = pd.get_dummies(d3_reg['domain'], prefix='dom', drop_first=True)
    d3_dom_cols = list(d3_dom.columns)

    cov_lists = [
        ('Unadjusted', []),
        ('+ Age + Sex', ['age', 'sex']),
        ('+ Treatments', ['age', 'sex', 'on_statin', 'on_ezetimibe', 'on_pcsk9i']),
        ('Full', ['age', 'sex', 'on_statin', 'on_ezetimibe', 'on_pcsk9i', 'smoking', 'diabetes']),
    ]

    for name, covs in cov_lists:
        all_feats = d3_dom_cols + covs
        X = pd.concat([d3_dom.reset_index(drop=True),
                        d3_reg[covs].reset_index(drop=True)], axis=1)
        y_d3 = d3_reg['log_lpa'].reset_index(drop=True)
        valid = X.notna().all(axis=1) & y_d3.notna()
        X_c = X[valid].values.astype(float)
        y_c = y_d3[valid].values
        if len(y_c) >= 20:
            lr_d3 = LinearRegression()
            lr_d3.fit(X_c, y_c)
            r2 = lr_d3.score(X_c, y_c)
            print(f"\n    DRAGON_3 {name} (n={len(y_c)}): R²={r2:.4f}")
            for feat, coef in zip(all_feats, lr_d3.coef_):
                if feat.startswith('dom_'):
                    pct = (np.exp(coef) - 1) * 100
                    print(f"      {feat.replace('dom_',''):>25s}: {pct:>+6.1f}%")
                elif feat in covs:
                    print(f"      {feat:>25s}: beta={coef:.4f}")

# ============================================================================
# STEP 4: UNTREATED LDL BY DOMAIN
# ============================================================================
print(f"\n{'='*70}")
print("[4/7] UNTREATED LDL BY DOMAIN")
print("=" * 70)

# DRAGON_3 LDL_1 (first measurement — closest to untreated)
d3_ldl = d3[d3['ldl_untreated'].notna() & d3['domain'].notna()].copy()
print(f"\n  DRAGON_3 with untreated LDL + domain: {len(d3_ldl)}")

if len(d3_ldl) >= 20:
    dom_ldl = d3_ldl.groupby('domain').agg(
        n=('ldl_untreated', 'count'),
        median_ldl=('ldl_untreated', 'median'),
        mean_ldl=('ldl_untreated', 'mean'),
    ).reset_index()
    dom_ldl = dom_ldl[dom_ldl['n'] >= 3].sort_values('median_ldl', ascending=False)

    # Also get Lp(a) by domain from same patients
    dom_both = d3[d3['ldl_untreated'].notna() & d3['lpa'].notna() & d3['domain'].notna()]
    dom_lpa_d3 = dom_both.groupby('domain').agg(
        n_lpa=('lpa', 'count'),
        median_lpa=('lpa', 'median'),
    ).reset_index()

    dom_compare = dom_ldl.merge(dom_lpa_d3, on='domain', how='left')
    dom_compare['penetrance'] = dom_compare['domain'].map(
        atlas.groupby('domain')['domain_penetrance_by_60_pct'].first().to_dict())

    print(f"\n  {'Domain':>25s}  {'N':>4s}  {'Med LDL':>8s}  {'Med Lpa':>8s}  {'Pen%':>5s}")
    print("  " + "-" * 60)
    for _, r in dom_compare.iterrows():
        lpa_str = f"{r['median_lpa']:.0f}" if pd.notna(r.get('median_lpa')) else 'N/A'
        pen = f"{r['penetrance']:.1f}" if pd.notna(r.get('penetrance')) else 'N/A'
        print(f"  {r['domain']:>25s}  {int(r['n']):>4d}  {r['median_ldl']:>8.2f}  "
              f"{lpa_str:>8s}  {pen:>5s}")

# UKB LDL by domain
ukb_ldl = ukb[ukb['ldl'].notna() & ukb['domain'].notna()]
if len(ukb_ldl) >= 50:
    ukb_dom_ldl = ukb_ldl.groupby('domain').agg(
        n=('ldl', 'count'),
        median_ldl=('ldl', 'median'),
        median_lpa=('lpa', 'median'),
    ).reset_index()
    ukb_dom_ldl = ukb_dom_ldl[ukb_dom_ldl['n'] >= 10].sort_values('median_ldl', ascending=False)

    print(f"\n  UKB carriers LDL by domain (n={len(ukb_ldl)}, NOTE: mostly treated):")
    for _, r in ukb_dom_ldl.iterrows():
        print(f"  {r['domain']:>25s}  n={int(r['n']):>4d}  LDL={r['median_ldl']:>5.2f}  "
              f"Lp(a)={r['median_lpa']:>5.1f}")

    # Correlation: domain-level LDL vs Lp(a)
    if len(ukb_dom_ldl) >= 5:
        rho, p = stats.spearmanr(ukb_dom_ldl['median_ldl'], ukb_dom_ldl['median_lpa'])
        print(f"\n  Domain-level LDL vs Lp(a) correlation: rho={rho:.3f}, P={p:.4f}")

# ============================================================================
# STEP 5: MUTATION-LEVEL PROTEIN POSITION ANALYSIS
# ============================================================================
print(f"\n{'='*70}")
print("[5/7] MUTATION-LEVEL: Lp(a) vs LDL at protein position")
print("=" * 70)

# DRAGON_3: mutations with both Lp(a) and LDL_1
d3_both = d3[d3['lpa'].notna() & d3['ldl_untreated'].notna() & d3['Mutation1'].notna()].copy()
mut_level = d3_both.groupby('Mutation1').agg(
    n=('lpa', 'count'),
    median_lpa=('lpa', 'median'),
    median_ldl=('ldl_untreated', 'median'),
    protein_pos=('protein_pos', 'first'),
    domain=('domain', 'first'),
).reset_index()
mut_level = mut_level[mut_level['n'] >= 3].copy()

print(f"  Mutations with n>=3 and both Lp(a) + LDL: {len(mut_level)}")

if len(mut_level) >= 5:
    # Correlation: mutation-level Lp(a) vs LDL
    rho_ml, p_ml = stats.spearmanr(mut_level['median_lpa'], mut_level['median_ldl'])
    print(f"  Mutation-level Lp(a) vs untreated LDL: rho={rho_ml:.3f}, P={p_ml:.4f}")

    # Identify divergent mutations (high Lp(a) relative to LDL)
    # Standardise both
    mut_level['lpa_z'] = (mut_level['median_lpa'] - mut_level['median_lpa'].mean()) / mut_level['median_lpa'].std()
    mut_level['ldl_z'] = (mut_level['median_ldl'] - mut_level['median_ldl'].mean()) / mut_level['median_ldl'].std()
    mut_level['divergence'] = mut_level['lpa_z'] - mut_level['ldl_z']  # positive = Lp(a) disproportionately high

    print(f"\n  Mutations with highest Lp(a)/LDL divergence (Lp(a) > expected from LDL):")
    top_div = mut_level.sort_values('divergence', ascending=False).head(8)
    for _, r in top_div.iterrows():
        dom = r['domain'] if pd.notna(r['domain']) else '?'
        print(f"    {r['Mutation1']:>35s}: Lp(a)={r['median_lpa']:>5.0f}, LDL={r['median_ldl']:>5.2f}, "
              f"div={r['divergence']:>+5.2f}, {dom}")

    print(f"\n  Mutations with lowest divergence (LDL > expected from Lp(a)):")
    bot_div = mut_level.sort_values('divergence').head(5)
    for _, r in bot_div.iterrows():
        dom = r['domain'] if pd.notna(r['domain']) else '?'
        print(f"    {r['Mutation1']:>35s}: Lp(a)={r['median_lpa']:>5.0f}, LDL={r['median_ldl']:>5.2f}, "
              f"div={r['divergence']:>+5.2f}, {dom}")

mut_level.to_csv(os.path.join(ANALYSIS, "lpa_vs_ldl_mutation_level.csv"), index=False)

# ============================================================================
# STEP 6: SENSITIVITY ANALYSES
# ============================================================================
print(f"\n{'='*70}")
print("[6/7] SENSITIVITY ANALYSES")
print("=" * 70)

# 6a. Exclude PCSK9i users (DRAGON_3)
d3_no_pcsk9i = d3_mapped[d3_mapped['on_pcsk9i'] != 1]
if len(d3_no_pcsk9i) >= 20:
    dom_sens = d3_no_pcsk9i.groupby('domain').agg(
        n=('lpa', 'count'), median=('lpa', 'median')).reset_index()
    dom_sens = dom_sens[dom_sens['n'] >= 3]
    print(f"\n  Excluding PCSK9i users (DRAGON_3, n={len(d3_no_pcsk9i)}):")
    for _, r in dom_sens.sort_values('median', ascending=False).iterrows():
        print(f"    {r['domain']:>25s}: n={int(r['n'])}, median Lp(a)={r['median']:.0f}")

    groups_sens = [g['lpa'].values for _, g in d3_no_pcsk9i.groupby('domain') if len(g) >= 3]
    if len(groups_sens) >= 3:
        h, p = stats.kruskal(*groups_sens)
        print(f"  Kruskal-Wallis (no PCSK9i): H={h:.2f}, P={p:.4f}")

# 6b. UKB Kruskal-Wallis (reconfirm)
groups_ukb = [g['lpa'].values for _, g in ukb[ukb['domain'].notna()].groupby('domain') if len(g) >= 10]
if len(groups_ukb) >= 3:
    h_ukb, p_ukb = stats.kruskal(*groups_ukb)
    print(f"\n  UKB Kruskal-Wallis (all carriers): H={h_ukb:.2f}, P={p_ukb:.4f}")

# Pairwise domain comparisons (top vs bottom)
ukb_dom = ukb[ukb['domain'].notna()]
top_dom = 'EGF-like B'
bot_dom = 'Ligand-binding R4'
top_vals = ukb_dom[ukb_dom['domain'] == top_dom]['lpa']
bot_vals = ukb_dom[ukb_dom['domain'] == bot_dom]['lpa']
if len(top_vals) >= 10 and len(bot_vals) >= 10:
    u, p_pw = stats.mannwhitneyu(top_vals, bot_vals)
    d_cohen = (top_vals.mean() - bot_vals.mean()) / np.sqrt((top_vals.std()**2 + bot_vals.std()**2) / 2)
    print(f"\n  {top_dom} vs {bot_dom}:")
    print(f"    Median: {top_vals.median():.1f} vs {bot_vals.median():.1f}")
    print(f"    P={p_pw:.4f}, Cohen d={d_cohen:.3f}")

# ============================================================================
# STEP 7: PUBLICATION FIGURE
# ============================================================================
print(f"\n{'='*70}")
print("[7/7] Generating 6-panel figure...")

RED = '#b2182b'
BLUE = '#2166ac'
ORANGE = '#ef8a62'
GREEN = '#1b7837'

fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# Panel A: Lp(a) by domain — adjusted (forest-style with CI)
ax = axes[0, 0]
dom_stats = ukb[ukb['domain'].notna()].groupby('domain').agg(
    n=('lpa', 'count'), mean=('lpa', 'mean'), std=('lpa', 'std'),
    median=('lpa', 'median')).reset_index()
dom_stats = dom_stats[dom_stats['n'] >= 10].sort_values('median')
dom_stats['se'] = dom_stats['std'] / np.sqrt(dom_stats['n'])
dom_stats['ci_lo'] = dom_stats['mean'] - 1.96 * dom_stats['se']
dom_stats['ci_hi'] = dom_stats['mean'] + 1.96 * dom_stats['se']

y_pos = range(len(dom_stats))
ax.errorbar(dom_stats['mean'], y_pos, xerr=1.96*dom_stats['se'],
            fmt='o', color=RED, capsize=3, markersize=6)
ax.set_yticks(y_pos)
ax.set_yticklabels([f"{r['domain'].replace('Ligand-binding ','LB-').replace('EGF-like ','EGF-')}\n(n={int(r['n'])})"
                     for _, r in dom_stats.iterrows()], fontsize=7)
ax.axvline(x=ukb['lpa'].mean(), color='grey', linestyle='--', linewidth=0.8)
ax.set_xlabel('Mean Lp(a) (nmol/L) [95% CI]', fontsize=10)
ax.set_title('A. Lp(a) by LDLR Domain\n(UKB, Kruskal-Wallis P=0.0005)', fontsize=11, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: Untreated LDL by domain (DRAGON_3 or UKB)
ax = axes[0, 1]
if len(ukb_dom_ldl) > 0:
    ldl_sorted = ukb_dom_ldl.sort_values('median_ldl')
    ax.barh(range(len(ldl_sorted)), ldl_sorted['median_ldl'], color=BLUE,
            edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(ldl_sorted)))
    ax.set_yticklabels([f"{r['domain'].replace('Ligand-binding ','LB-').replace('EGF-like ','EGF-')}\n(n={int(r['n'])})"
                         for _, r in ldl_sorted.iterrows()], fontsize=7)
    ax.set_xlabel('Median LDL-C (mmol/L)', fontsize=10)
    ax.set_title('B. LDL-C by Domain (UKB)', fontsize=11, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: Mutation-level Lp(a) vs LDL scatter
ax = axes[0, 2]
if len(mut_level) >= 5:
    dom_colours = {'Ligand-binding R3': BLUE, 'Ligand-binding R4': BLUE,
                   'Ligand-binding R6': BLUE, 'Ligand-binding R7': BLUE,
                   'EGF-like B': RED, 'EGF-like C': RED,
                   'Beta-propeller': ORANGE, 'O-linked sugar': GREEN,
                   'Linker (EGF-C/O-linked)': ORANGE}
    colours = [dom_colours.get(d, '#999999') for d in mut_level['domain']]
    ax.scatter(mut_level['median_ldl'], mut_level['median_lpa'], c=colours,
               s=mut_level['n'] * 15, alpha=0.7, edgecolors='black', linewidth=0.5)
    if not np.isnan(rho_ml):
        ax.set_title(f'C. Lp(a) vs LDL at Mutation Level\n(rho={rho_ml:.3f}, P={p_ml:.3f})',
                     fontsize=11, fontweight='bold')
    ax.set_xlabel('Median Untreated LDL (mmol/L)', fontsize=10)
    ax.set_ylabel('Median Lp(a) (nmol/L)', fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Protein position vs Lp(a)
ax = axes[1, 0]
ukb_pos = ukb[ukb['protein_pos'].notna() & ukb['lpa'].notna()].copy()
ukb_pos['protein_pos'] = pd.to_numeric(ukb_pos['protein_pos'], errors='coerce')
if len(ukb_pos) >= 50:
    # Bin by protein position (50-residue windows)
    ukb_pos['pos_bin'] = (ukb_pos['protein_pos'] // 50) * 50
    pos_agg = ukb_pos.groupby('pos_bin').agg(
        n=('lpa', 'count'), median_lpa=('lpa', 'median')).reset_index()
    pos_agg = pos_agg[pos_agg['n'] >= 5]
    ax.bar(pos_agg['pos_bin'], pos_agg['median_lpa'], width=45,
           color=BLUE, edgecolor='black', linewidth=0.3, alpha=0.7)
    # Domain boundaries
    for pos, dom in [(22, 'LB'), (353, 'EGF-A'), (446, 'β-prop'), (593, 'EGF-C'), (693, 'O-sugar')]:
        ax.axvline(x=pos, color='red', linestyle=':', linewidth=0.5, alpha=0.5)
        ax.text(pos, ax.get_ylim()[1]*0.95, dom, fontsize=5, rotation=90, va='top')
    ax.set_xlabel('LDLR Protein Position', fontsize=10)
    ax.set_ylabel('Median Lp(a) (nmol/L)', fontsize=10)
    ax.set_title('D. Lp(a) Along LDLR Protein', fontsize=11, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel E: Model R² comparison
ax = axes[1, 1]
if model_results:
    mr = pd.DataFrame(model_results)
    bars = ax.barh(range(len(mr)), mr['r2'], color=BLUE, edgecolor='black', linewidth=0.5)
    ax.barh(range(len(mr)), mr['r2_domain_only'], color=RED, edgecolor='black',
            linewidth=0.5, alpha=0.7, label='Domain contribution')
    ax.set_yticks(range(len(mr)))
    ax.set_yticklabels([m.split('.')[1].strip()[:25] for m in mr['model']], fontsize=8)
    ax.set_xlabel('R²', fontsize=10)
    ax.set_title('E. Model R² (Domain Contribution)', fontsize=11, fontweight='bold')
    ax.legend(fontsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel F: Domain-level Lp(a) vs LDL (divergence)
ax = axes[1, 2]
if len(ukb_dom_ldl) >= 5:
    # Standardise both for comparison
    ldl_z = (ukb_dom_ldl['median_ldl'] - ukb_dom_ldl['median_ldl'].mean()) / ukb_dom_ldl['median_ldl'].std()
    lpa_z = (ukb_dom_ldl['median_lpa'] - ukb_dom_ldl['median_lpa'].mean()) / ukb_dom_ldl['median_lpa'].std()
    div = lpa_z.values - ldl_z.values

    colours_f = [RED if d > 0.3 else BLUE if d < -0.3 else '#999999' for d in div]
    doms_short = [d.replace('Ligand-binding ', 'LB-').replace('EGF-like ', 'EGF-')
                   for d in ukb_dom_ldl['domain']]
    ax.barh(range(len(div)), div, color=colours_f, edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(div)))
    ax.set_yticklabels(doms_short, fontsize=8)
    ax.axvline(x=0, color='black', linewidth=1)
    ax.set_xlabel('Lp(a)/LDL Divergence (z-score)', fontsize=10)
    ax.set_title('F. Lp(a) vs LDL Divergence\n(Red=Lp(a) disproportionately high)', fontsize=11, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout(pad=2.0)
fig_path = os.path.join(FIGURES, "Figure_Lpa_Adjusted_Mutation_Level.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {fig_path}")

# Save domain results
if model_results:
    pd.DataFrame(model_results).to_csv(os.path.join(ANALYSIS, "lpa_adjusted_by_domain.csv"), index=False)

print("\n" + "=" * 70)
print("COMPLETE")
print("=" * 70)

#!/usr/bin/env python3
"""
45_nobel_critical_analyses.py
==============================
The 6 analyses that transform good science into Nobel-calibre discovery.

1. SSS x Treatment Interaction (does SSS predict WHO benefits from treatment?)
2. Saturation ddG vs ClinVar Pathogenicity (validates the entire approach)
3. SSS Recalibration using penetrance as ground truth
4. Islam functional activity vs FoldX ddG (direct structure-function)
5. Domain-level treatment selection model
6. The One Experiment: untreated LDL by domain x ddG

Outputs:
  alphafold/analysis/nobel_treatment_interaction.csv
  alphafold/analysis/nobel_saturation_vs_clinvar.csv
  alphafold/analysis/nobel_sss_recalibrated.csv
  alphafold/analysis/nobel_ddg_vs_function.csv
  alphafold/analysis/nobel_treatment_selection.csv
  alphafold/analysis/nobel_summary.csv
"""

import pandas as pd
import numpy as np
import os
import json
import warnings
warnings.filterwarnings('ignore')
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")

print("=" * 80)
print("NOBEL-CRITICAL ANALYSES")
print("Structure -> Function -> Treatment -> Outcome")
print("=" * 80)

# Load core datasets
comp = pd.read_csv(os.path.join(ANALYSIS, "comprehensive_sss_analysis.csv"))
atlas = pd.read_csv(os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas_v3.csv"))
sss_cat = pd.read_csv(os.path.join(ANALYSIS, "sss_v3_catalogue.csv"))
sat = pd.read_csv(os.path.join(ANALYSIS, "foldx_saturation_mutagenesis.csv"))
km = pd.read_csv(os.path.join(ANALYSIS, "cox_km_compound_results.csv"))
llt = pd.read_csv(os.path.join(ANALYSIS, "variant_llt_response.csv"))

# Load Wales full cohort
wales_file = os.path.join(BASE, "WALES_FH_CLEANED (1) - Copy.csv")
wales = pd.read_csv(wales_file)
dragon = pd.read_csv(os.path.join(BASE, "DRAGON_3.csv"))

print(f"  Comprehensive: {len(comp)} patients")
print(f"  Wales full: {len(wales)} patients")
print(f"  Dragon3: {len(dragon)} patients")
print(f"  Atlas: {len(atlas)} positions")
print(f"  SSS catalogue: {len(sss_cat)} variants")
print(f"  Saturation: {len(sat)} mutations")
print(f"  LLT response: {len(llt)} patients")

# Load Islam functional data
islam_file = os.path.join(ANALYSIS, "islam_et_al_315_variants.csv")
try:
    islam = pd.read_csv(islam_file)
    print(f"  Islam et al: {len(islam)} variants")
except:
    islam = None
    print("  Islam et al: NOT FOUND")

# Load ICD-10 for better ASCVD definition
icd10_file = os.path.join(BASE, "ukb_reviewer_icd10_codes.csv")
try:
    icd10 = pd.read_csv(icd10_file)
    icd10.columns = [c.replace('participant.', '') for c in icd10.columns]
    print(f"  ICD-10 codes: {len(icd10)}")
except:
    icd10 = None

# ============================================================================
# ANALYSIS 1: SSS x TREATMENT INTERACTION
# "Does SSS predict WHO benefits from treatment escalation?"
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 1: SSS x TREATMENT INTERACTION")
print("The key question: Does mutation severity predict treatment benefit?")
print("=" * 80)

# Use Dragon3 (has detailed treatment data)
# Merge SSS to Dragon3 via mutation
dragon_cols = [c for c in dragon.columns]
mutation_col = None
for c in ['Mutation1', 'Positive1', 'mutation']:
    if c in dragon_cols:
        mutation_col = c
        break

if mutation_col:
    # Map SSS to dragon via mutation name
    sss_map = dict(zip(sss_cat['variant_id'], sss_cat['sss_v3']))
    domain_map = dict(zip(sss_cat['variant_id'], sss_cat['domain']))

    # Try matching
    dragon['sss'] = dragon[mutation_col].map(sss_map)
    dragon['domain'] = dragon[mutation_col].map(domain_map)

    # Also try with LDLR: prefix
    dragon.loc[dragon['sss'].isna(), 'sss'] = dragon.loc[dragon['sss'].isna(), mutation_col].apply(
        lambda x: sss_map.get(f"LDLR:{x}", np.nan) if pd.notna(x) else np.nan
    )
    dragon.loc[dragon['domain'].isna(), 'domain'] = dragon.loc[dragon['domain'].isna(), mutation_col].apply(
        lambda x: domain_map.get(f"LDLR:{x}", np.nan) if pd.notna(x) else np.nan
    )

    matched = dragon['sss'].notna().sum()
    print(f"\n  Dragon3 patients with SSS: {matched}/{len(dragon)}")

# Use comprehensive analysis dataset (cleaner)
# Treatment interaction: ASCVD ~ SSS * on_statin + age + sex
cc = comp.dropna(subset=['ascvd', 'sss', 'on_statin', 'age', 'sex', 'ldl1'])
print(f"  Complete cases for interaction: {len(cc)}")

# Create interaction term
cc['sss_x_statin'] = cc['sss'] * cc['on_statin']

# Model without interaction
X_base = cc[['age', 'sex', 'sss', 'on_statin']].values
y = cc['ascvd'].values
model_base = LogisticRegression(max_iter=5000, C=1e6)
model_base.fit(X_base, y)
auc_base = roc_auc_score(y, model_base.predict_proba(X_base)[:, 1])

# Model WITH interaction
X_int = cc[['age', 'sex', 'sss', 'on_statin', 'sss_x_statin']].values
model_int = LogisticRegression(max_iter=5000, C=1e6)
model_int.fit(X_int, y)
auc_int = roc_auc_score(y, model_int.predict_proba(X_int)[:, 1])

print(f"\n  Without interaction: AUC = {auc_base:.4f}")
print(f"  With SSS*statin interaction: AUC = {auc_int:.4f}")
print(f"  Improvement: +{(auc_int - auc_base):.4f}")

# Interaction coefficient
int_coef = model_int.coef_[0]
var_names = ['age', 'sex', 'sss', 'on_statin', 'sss_x_statin']
print(f"\n  Interaction model coefficients:")
for name, coef in zip(var_names, int_coef):
    print(f"    {name:20s}: beta={coef:8.4f}  OR={np.exp(coef):.3f}")

# Stratified analysis: ASCVD rate by SSS tertile x treatment
cc['sss_tertile'] = pd.qcut(cc['sss'], 3, labels=['Low', 'Mid', 'High'])
interaction_table = cc.groupby(['sss_tertile', 'on_statin']).agg(
    n=('ascvd', 'count'),
    events=('ascvd', 'sum'),
    rate=('ascvd', 'mean'),
    mean_sss=('sss', 'mean'),
).reset_index()
interaction_table['rate_pct'] = (interaction_table['rate'] * 100).round(1)

print(f"\n  ASCVD Rate by SSS Tertile x Treatment:")
print(f"  {'SSS Tertile':>12s} {'Treated':>8s} {'N':>6s} {'Events':>7s} {'ASCVD%':>7s}")
for _, r in interaction_table.iterrows():
    tx = 'Yes' if r['on_statin'] == 1 else 'No'
    print(f"  {r['sss_tertile']:>12s} {tx:>8s} {int(r['n']):>6d} {int(r['events']):>7d} {r['rate_pct']:>7.1f}")

# The KEY insight: treatment benefit = untreated rate - treated rate
# If SSS predicts treatment benefit, high SSS patients gain MORE from treatment
for tert in ['Low', 'Mid', 'High']:
    treated = interaction_table[(interaction_table['sss_tertile'] == tert) & (interaction_table['on_statin'] == 1)]
    untreated = interaction_table[(interaction_table['sss_tertile'] == tert) & (interaction_table['on_statin'] == 0)]
    if len(treated) > 0 and len(untreated) > 0:
        benefit = untreated.iloc[0]['rate_pct'] - treated.iloc[0]['rate_pct']
        # Note: negative benefit means treated have HIGHER rate (treatment paradox)
        print(f"  -> {tert} SSS: treatment effect = {benefit:+.1f}% points")

interaction_table.to_csv(os.path.join(ANALYSIS, "nobel_treatment_interaction.csv"), index=False)

# ============================================================================
# ANALYSIS 2: SATURATION ddG vs ClinVar PATHOGENICITY
# "Does FoldX correctly identify damaging mutations?"
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 2: SATURATION ddG vs ClinVar PATHOGENICITY")
print("The foundational validation: structure predicts pathogenicity")
print("=" * 80)

# Load ClinVar data
clinvar_file = os.path.join(ANALYSIS, "ukb_clinvar_full.csv")
vep_file = os.path.join(ANALYSIS, "ukb_ldlr_vep_annotations.csv")

try:
    clinvar = pd.read_csv(clinvar_file)
    vep = pd.read_csv(vep_file)
    print(f"  ClinVar: {len(clinvar)} variants")
    print(f"  VEP: {len(vep)} annotations")

    # Get protein positions from VEP
    vep_miss = vep[vep['amino_acids'].str.contains('/', na=False)].copy()
    vep_miss['wt_aa'] = vep_miss['amino_acids'].str.split('/').str[0]
    vep_miss['mut_aa'] = vep_miss['amino_acids'].str.split('/').str[1]
    vep_miss['position'] = pd.to_numeric(vep_miss['protein_position'], errors='coerce')
    vep_miss = vep_miss.dropna(subset=['position'])
    vep_miss['position'] = vep_miss['position'].astype(int)

    # Match VEP to ClinVar
    vep_cv = vep_miss.merge(clinvar, on=['pos', 'ref', 'alt'], how='left')

    # Match to saturation data
    sat_valid = sat.dropna(subset=['ddG']).copy()
    sat_valid['position'] = sat_valid['position'].astype(int)

    # Direct match: same position, same wt, same mut
    matched = vep_cv.merge(sat_valid[['position', 'wt_aa', 'mut_aa', 'ddG', 'effect_class']],
                           on=['position', 'wt_aa', 'mut_aa'], how='inner')

    print(f"  VEP-ClinVar-Saturation matched: {len(matched)} variants")

    if len(matched) > 0:
        # Classify by ClinVar
        matched['cv_class'] = 'Unknown'
        matched.loc[matched['clinvar'].str.contains('athogenic', na=False), 'cv_class'] = 'Pathogenic'
        matched.loc[matched['clinvar'].str.contains('ikely_pathogenic', na=False), 'cv_class'] = 'Likely Pathogenic'
        matched.loc[matched['clinvar'].str.contains('enign', na=False), 'cv_class'] = 'Benign'
        matched.loc[matched['clinvar'].str.contains('ikely_benign', na=False), 'cv_class'] = 'Likely Benign'
        matched.loc[matched['clinvar'].str.contains('uncertain|VUS', na=False, case=False), 'cv_class'] = 'VUS'

        print(f"\n  ClinVar classification of matched variants:")
        for cls in ['Pathogenic', 'Likely Pathogenic', 'VUS', 'Likely Benign', 'Benign', 'Unknown']:
            n = (matched['cv_class'] == cls).sum()
            if n > 0:
                ddg_med = matched[matched['cv_class'] == cls]['ddG'].median()
                print(f"    {cls:20s}: n={n:3d}, median ddG={ddg_med:6.2f} kcal/mol")

        # Binary comparison: Pathogenic/Likely Path vs Benign/Likely Benign
        path = matched[matched['cv_class'].isin(['Pathogenic', 'Likely Pathogenic'])]['ddG']
        ben = matched[matched['cv_class'].isin(['Benign', 'Likely Benign'])]['ddG']

        if len(path) > 0 and len(ben) > 0:
            u, p = stats.mannwhitneyu(path, ben, alternative='greater')
            auc_cv = roc_auc_score(
                [1]*len(path) + [0]*len(ben),
                list(path) + list(ben)
            )
            print(f"\n  Pathogenic (n={len(path)}) vs Benign (n={len(ben)}):")
            print(f"    Pathogenic median ddG: {path.median():.2f} kcal/mol")
            print(f"    Benign median ddG:     {ben.median():.2f} kcal/mol")
            print(f"    Mann-Whitney P = {p:.6f} (one-sided: pathogenic > benign)")
            print(f"    AUC (ddG discriminates pathogenic): {auc_cv:.3f}")
        elif len(path) > 0:
            print(f"\n  Only pathogenic variants found (n={len(path)}), no benign for comparison")
            print(f"    Pathogenic median ddG: {path.median():.2f} kcal/mol")

    # ALSO: use position-level max ddG from atlas
    # Positions with known pathogenic variants should have higher max ddG
    print(f"\n  Position-level analysis:")
    # Get positions with pathogenic variants
    path_positions = set(matched[matched['cv_class'].isin(['Pathogenic', 'Likely Pathogenic'])]['position'].unique())
    ben_positions = set(matched[matched['cv_class'].isin(['Benign', 'Likely Benign'])]['position'].unique())
    # Remove overlap
    path_only = path_positions - ben_positions
    ben_only = ben_positions - path_positions

    if path_only and ben_only:
        path_max_ddg = atlas[atlas['position'].isin(path_only)]['sat_max_ddG'].dropna()
        ben_max_ddg = atlas[atlas['position'].isin(ben_only)]['sat_max_ddG'].dropna()
        if len(path_max_ddg) > 2 and len(ben_max_ddg) > 2:
            u2, p2 = stats.mannwhitneyu(path_max_ddg, ben_max_ddg, alternative='greater')
            print(f"    Positions with pathogenic variants (n={len(path_max_ddg)}): max ddG median = {path_max_ddg.median():.2f}")
            print(f"    Positions with benign variants (n={len(ben_max_ddg)}): max ddG median = {ben_max_ddg.median():.2f}")
            print(f"    P = {p2:.6f}")

    matched.to_csv(os.path.join(ANALYSIS, "nobel_saturation_vs_clinvar.csv"), index=False)

except Exception as e:
    print(f"  ERROR: {e}")
    import traceback
    traceback.print_exc()

# ============================================================================
# ANALYSIS 3: SSS RECALIBRATION
# "Optimize SSS weights using domain penetrance as ground truth"
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 3: SSS RECALIBRATION")
print("Optimize weights using clinical outcome as ground truth")
print("=" * 80)

# Use domain-level penetrance as ground truth
km_domain = km[km['analysis'] == 'KM_domain'].copy()

# Map domain penetrance to SSS components
sss_ldlr = sss_cat[sss_cat['gene'] == 'LDLR'].copy()

# For each domain, get mean SSS components
domain_components = sss_ldlr.groupby('domain').agg(
    mean_domain_risk=('domain_risk', 'mean'),
    mean_type_severity=('type_severity', 'mean'),
    mean_ddg_score=('ddG_score', 'mean'),
    mean_plddt_disruption=('plddt_disruption', 'mean'),
    mean_interface_proximity=('interface_proximity', 'mean'),
    mean_sss_v3=('sss_v3', 'mean'),
    n_variants=('variant_id', 'count'),
).reset_index()

# Merge with penetrance
domain_components = domain_components.merge(
    km_domain[['group', 'value', 'n', 'events', 'rate']].rename(
        columns={'group': 'domain', 'value': 'penetrance_60'}
    ),
    on='domain', how='inner'
)

print(f"\n  Domains with both SSS and penetrance: {len(domain_components)}")

if len(domain_components) >= 5:
    # Test each component vs penetrance
    print(f"\n  Individual component correlations with penetrance:")
    components = ['mean_domain_risk', 'mean_type_severity', 'mean_ddg_score',
                  'mean_plddt_disruption', 'mean_interface_proximity', 'mean_sss_v3']
    best_r = 0
    best_comp = ''
    for comp_name in components:
        valid = domain_components.dropna(subset=[comp_name, 'penetrance_60'])
        if len(valid) >= 5:
            r, p = stats.spearmanr(valid[comp_name], valid['penetrance_60'])
            marker = ' ***' if p < 0.05 else ' *' if p < 0.1 else ''
            print(f"    {comp_name:30s}: rho={r:6.3f}, P={p:.4f}{marker}")
            if abs(r) > abs(best_r):
                best_r = r
                best_comp = comp_name

    # Also test: max ddG from atlas vs penetrance
    atlas_domain = atlas.groupby('domain').agg(
        mean_max_ddg=('sat_max_ddG', 'mean'),
        median_max_ddg=('sat_max_ddG', 'median'),
        mean_sensitivity=('mutational_sensitivity', 'mean'),
    ).reset_index()

    atlas_pen = atlas_domain.merge(
        km_domain[['group', 'value']].rename(columns={'group': 'domain', 'value': 'penetrance_60'}),
        on='domain', how='inner'
    )

    print(f"\n  Saturation mutagenesis vs penetrance:")
    for col in ['mean_max_ddg', 'median_max_ddg', 'mean_sensitivity']:
        valid = atlas_pen.dropna(subset=[col, 'penetrance_60'])
        if len(valid) >= 5:
            r, p = stats.spearmanr(valid[col], valid['penetrance_60'])
            marker = ' ***' if p < 0.05 else ' *' if p < 0.1 else ''
            print(f"    {col:30s}: rho={r:6.3f}, P={p:.4f}{marker}")

    # Recalibrate SSS using optimal weights
    # Use ASCVD rate as target in patient-level data
    print(f"\n  Recalibrating SSS weights...")

    # Patient-level: use comprehensive analysis
    comp_valid = comp.dropna(subset=['sss', 'ascvd', 'domain']).copy()

    # Get domain-level ddG from atlas
    domain_ddg = atlas.groupby('domain')['sat_max_ddG'].mean().to_dict()
    domain_sensitivity = atlas.groupby('domain')['mutational_sensitivity'].mean().to_dict()

    comp_valid['domain_ddg'] = comp_valid['domain'].map(domain_ddg)
    comp_valid['domain_sensitivity'] = comp_valid['domain'].map(domain_sensitivity)

    # Test domain_ddg as predictor
    valid_ddg = comp_valid.dropna(subset=['domain_ddg'])
    if len(valid_ddg) > 20:
        r_ddg, p_ddg = stats.pointbiserialr(valid_ddg['ascvd'], valid_ddg['domain_ddg'])
        print(f"    Domain mean max ddG vs ASCVD: rpb={r_ddg:.3f}, P={p_ddg:.4f}")

    # Build recalibrated SSS using logistic regression
    # Features: domain penetrance, domain ddG, individual ddG_score, interface
    sss_features = comp_valid.dropna(subset=['domain_ddg']).copy()
    if len(sss_features) > 50:
        # Map penetrance
        pen_map = dict(zip(km_domain['group'], km_domain['value']))
        sss_features['domain_penetrance'] = sss_features['domain'].map(pen_map)
        sss_features = sss_features.dropna(subset=['domain_penetrance'])

        if len(sss_features) > 30:
            X_recal = sss_features[['domain_penetrance', 'domain_ddg', 'sss']].values
            y_recal = sss_features['ascvd'].values

            model_recal = LogisticRegression(max_iter=5000, C=1e6)
            model_recal.fit(X_recal, y_recal)
            auc_recal = roc_auc_score(y_recal, model_recal.predict_proba(X_recal)[:, 1])

            # Compare with SSS alone
            X_sss = sss_features[['sss']].values
            model_sss = LogisticRegression(max_iter=5000, C=1e6)
            model_sss.fit(X_sss, y_recal)
            auc_sss = roc_auc_score(y_recal, model_sss.predict_proba(X_sss)[:, 1])

            print(f"\n    SSS alone -> ASCVD: AUC = {auc_sss:.4f}")
            print(f"    Penetrance + ddG + SSS -> ASCVD: AUC = {auc_recal:.4f}")
            print(f"    Improvement: +{(auc_recal - auc_sss):.4f}")

            # Weights
            print(f"    Recalibrated weights:")
            for name, w in zip(['penetrance', 'domain_ddg', 'sss'], model_recal.coef_[0]):
                print(f"      {name:20s}: {w:.4f}")

    domain_components.to_csv(os.path.join(ANALYSIS, "nobel_sss_recalibrated.csv"), index=False)

# ============================================================================
# ANALYSIS 4: FoldX ddG vs ISLAM FUNCTIONAL ACTIVITY
# "Does thermodynamic instability predict receptor dysfunction?"
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 4: FoldX ddG vs ISLAM FUNCTIONAL ACTIVITY")
print("Direct structure-to-function validation")
print("=" * 80)

if islam is not None:
    # Get FoldX ddG for Islam variants
    # Islam data has position, we need to match via position + amino acid
    islam_cols = islam.columns.tolist()
    print(f"  Islam columns: {islam_cols[:15]}")

    # Try different matching strategies
    # Strategy 1: Match by position to saturation max ddG
    islam_pos = islam.dropna(subset=['position']).copy()
    islam_pos['position'] = islam_pos['position'].astype(int)

    # Get max ddG per position from saturation
    sat_pos_max = sat.groupby('position')['ddG'].max().reset_index().rename(columns={'ddG': 'sat_max_ddG'})
    sat_pos_mean = sat.groupby('position')['ddG'].mean().reset_index().rename(columns={'ddG': 'sat_mean_ddG'})

    islam_merged = islam_pos.merge(sat_pos_max, on='position', how='inner')
    islam_merged = islam_merged.merge(sat_pos_mean, on='position', how='inner')

    print(f"  Islam variants matched to saturation: {len(islam_merged)}")

    # Test correlations
    activity_col = None
    for c in ['activity_pct', 'activity', 'mean_activity', 'ldl_uptake_pct',
              'surface_expression_pct', 'binding_pct']:
        if c in islam_merged.columns and islam_merged[c].notna().sum() > 10:
            activity_col = c
            break

    if activity_col:
        valid = islam_merged.dropna(subset=[activity_col, 'sat_max_ddG'])
        print(f"  Using activity column: {activity_col} (n={len(valid)})")

        # ddG vs activity (expect NEGATIVE: high ddG = low activity)
        r_max, p_max = stats.spearmanr(valid['sat_max_ddG'], valid[activity_col])
        r_mean, p_mean = stats.spearmanr(valid['sat_mean_ddG'], valid[activity_col])

        print(f"\n  Positional max ddG vs {activity_col}:")
        print(f"    rho = {r_max:.3f}, P = {p_max:.6f}")
        print(f"  Positional mean ddG vs {activity_col}:")
        print(f"    rho = {r_mean:.3f}, P = {p_mean:.6f}")

        # Binary: activity < 50% = dysfunctional
        valid['dysfunctional'] = (valid[activity_col] < 50).astype(int)
        n_dysfunc = valid['dysfunctional'].sum()
        n_func = len(valid) - n_dysfunc
        if n_dysfunc > 5 and n_func > 5:
            dysfunc_ddg = valid[valid['dysfunctional'] == 1]['sat_max_ddG']
            func_ddg = valid[valid['dysfunctional'] == 1]['sat_max_ddG']
            func_ddg_normal = valid[valid['dysfunctional'] == 0]['sat_max_ddG']
            u, p_bin = stats.mannwhitneyu(dysfunc_ddg, func_ddg_normal, alternative='greater')
            auc_func = roc_auc_score(valid['dysfunctional'], valid['sat_max_ddG'])
            print(f"\n  Dysfunctional (<50% activity, n={n_dysfunc}) vs Functional (n={n_func}):")
            print(f"    Dysfunctional max ddG: median = {dysfunc_ddg.median():.2f}")
            print(f"    Functional max ddG:    median = {func_ddg_normal.median():.2f}")
            print(f"    P = {p_bin:.6f}")
            print(f"    AUC (ddG predicts dysfunction): {auc_func:.3f}")

        # Also correlate with domain
        print(f"\n  Domain-level ddG vs activity:")
        domain_islam = valid.groupby('domain').agg(
            mean_activity=(activity_col, 'mean'),
            mean_ddg=('sat_max_ddG', 'mean'),
            n=('position', 'count'),
        ).reset_index()
        domain_islam = domain_islam[domain_islam['n'] >= 3]
        if len(domain_islam) >= 4:
            r_dom, p_dom = stats.spearmanr(domain_islam['mean_ddg'], domain_islam['mean_activity'])
            print(f"    rho = {r_dom:.3f}, P = {p_dom:.4f} (n={len(domain_islam)} domains)")
            for _, dr in domain_islam.sort_values('mean_activity').iterrows():
                print(f"      {dr['domain']:25s}: activity={dr['mean_activity']:5.1f}%, ddG={dr['mean_ddg']:5.2f} (n={int(dr['n'])})")

        islam_merged.to_csv(os.path.join(ANALYSIS, "nobel_ddg_vs_function.csv"), index=False)
    else:
        print("  No suitable activity column found")
        print(f"  Available columns: {islam_cols}")
else:
    print("  Islam data not available")

# ============================================================================
# ANALYSIS 5: DOMAIN-LEVEL TREATMENT SELECTION MODEL
# "Which drug works best for which domain?"
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 5: DOMAIN-LEVEL TREATMENT SELECTION")
print("The clinical breakthrough: structure-guided pharmacotherapy")
print("=" * 80)

# Use LLT response data
if len(llt) > 0:
    llt_cols = llt.columns.tolist()
    print(f"  LLT columns: {llt_cols[:15]}")

    # Find LDL change columns
    ldl_change_col = None
    for c in ['ldl_pct_reduction', 'ldl_pct_change', 'pct_change', 'LDL_reduction_pct', 'ldl_change']:
        if c in llt_cols:
            ldl_change_col = c
            break

    # Find domain column
    domain_col_llt = None
    for c in ['domain', 'Domain', 'DOMAIN']:
        if c in llt_cols:
            domain_col_llt = c
            break

    if ldl_change_col and domain_col_llt:
        llt_valid = llt.dropna(subset=[ldl_change_col, domain_col_llt])
        print(f"  LLT with domain + LDL change: {len(llt_valid)}")

        # Treatment response by domain
        domain_tx = llt_valid.groupby(domain_col_llt).agg(
            n=(ldl_change_col, 'count'),
            mean_reduction=(ldl_change_col, 'mean'),
            median_reduction=(ldl_change_col, 'median'),
            std=(ldl_change_col, 'std'),
        ).reset_index()
        domain_tx = domain_tx[domain_tx['n'] >= 3].sort_values('mean_reduction')

        print(f"\n  LDL Reduction by Domain (negative = better):")
        print(f"  {'Domain':>30s} {'N':>5s} {'Mean %':>8s} {'Median %':>9s}")
        for _, r in domain_tx.iterrows():
            print(f"  {r[domain_col_llt]:>30s} {int(r['n']):>5d} {r['mean_reduction']:>8.1f} {r['median_reduction']:>9.1f}")

        # Correlate treatment response with domain ddG
        domain_tx_ddg = domain_tx.merge(
            atlas.groupby('domain')['sat_max_ddG'].mean().reset_index().rename(
                columns={'sat_max_ddG': 'domain_max_ddg'}
            ),
            left_on=domain_col_llt, right_on='domain', how='inner'
        )

        if len(domain_tx_ddg) >= 4:
            r_tx, p_tx = stats.spearmanr(domain_tx_ddg['domain_max_ddg'], domain_tx_ddg['mean_reduction'])
            print(f"\n  Domain ddG vs treatment response: rho={r_tx:.3f}, P={p_tx:.4f}")
            print(f"  (Negative rho = higher ddG domains respond LESS to treatment)")

        domain_tx.to_csv(os.path.join(ANALYSIS, "nobel_treatment_selection.csv"), index=False)
    else:
        print(f"  Missing columns. LDL change: {ldl_change_col}, Domain: {domain_col_llt}")
        print(f"  Available: {llt_cols}")

# ============================================================================
# ANALYSIS 6: THE ONE EXPERIMENT
# "Untreated LDL by position-specific ddG"
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 6: THE ONE EXPERIMENT")
print("Untreated LDL correlates with position-specific structural damage")
print("=" * 80)

# Use comprehensive analysis: untreated patients only
untreated = comp[(comp['on_statin'] == 0) & comp['ldl1'].notna()].copy()
print(f"  Untreated patients with LDL: {len(untreated)}")

# Map domain ddG to each patient
untreated['domain_max_ddg'] = untreated['domain'].map(
    atlas.groupby('domain')['sat_max_ddG'].mean().to_dict()
)
untreated['domain_sensitivity'] = untreated['domain'].map(
    atlas.groupby('domain')['mutational_sensitivity'].mean().to_dict()
)

valid_ut = untreated.dropna(subset=['domain_max_ddg', 'ldl1'])
if len(valid_ut) > 10:
    r_ldl, p_ldl = stats.spearmanr(valid_ut['domain_max_ddg'], valid_ut['ldl1'])
    print(f"\n  Domain max ddG vs untreated LDL:")
    print(f"    rho = {r_ldl:.3f}, P = {p_ldl:.6f} (n={len(valid_ut)})")

    r_sens, p_sens = stats.spearmanr(valid_ut['domain_sensitivity'], valid_ut['ldl1'])
    print(f"  Domain sensitivity vs untreated LDL:")
    print(f"    rho = {r_sens:.3f}, P = {p_sens:.6f}")

    # Domain-level aggregation
    domain_ldl = valid_ut.groupby('domain').agg(
        n=('ldl1', 'count'),
        mean_ldl=('ldl1', 'mean'),
        mean_ddg=('domain_max_ddg', 'mean'),
    ).reset_index()
    domain_ldl = domain_ldl[domain_ldl['n'] >= 3].sort_values('mean_ldl', ascending=False)

    print(f"\n  Domain-level (n>=3):")
    print(f"  {'Domain':>30s} {'N':>5s} {'Mean LDL':>9s} {'Mean ddG':>9s}")
    for _, r in domain_ldl.iterrows():
        print(f"  {r['domain']:>30s} {int(r['n']):>5d} {r['mean_ldl']:>9.2f} {r['mean_ddg']:>9.2f}")

    if len(domain_ldl) >= 4:
        r_dom_ldl, p_dom_ldl = stats.spearmanr(domain_ldl['mean_ddg'], domain_ldl['mean_ldl'])
        print(f"\n  Domain-level ddG vs LDL: rho={r_dom_ldl:.3f}, P={p_dom_ldl:.4f}")


# ============================================================================
# FINAL SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("NOBEL-CRITICAL ANALYSES SUMMARY")
print("=" * 80)
print("""
  Files generated:
    nobel_treatment_interaction.csv
    nobel_saturation_vs_clinvar.csv
    nobel_sss_recalibrated.csv
    nobel_ddg_vs_function.csv
    nobel_treatment_selection.csv
""")

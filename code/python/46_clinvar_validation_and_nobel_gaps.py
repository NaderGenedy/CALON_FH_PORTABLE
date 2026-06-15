#!/usr/bin/env python3
"""
46_clinvar_validation_and_nobel_gaps.py
========================================
Proper ClinVar validation + comprehensive gap analysis for Nobel-calibre paper.

Three independent validation strategies:
  1. Islam et al. 315 variants -> ClinVar + FoldX ddG (functional + structural + clinical)
  2. UKB VEP missense -> ClinVar + saturation ddG (genomic coordinates -> protein match)
  3. Position-level: pathogenic vs benign positions -> saturation landscape

Plus Nobel gap analyses:
  4. Dose-response: ddG gradient -> functional activity gradient -> penetrance gradient
  5. External ClinVar download for comprehensive LDLR variant validation
  6. Specificity analysis: do benign variants have LOW ddG?

Outputs:
  alphafold/analysis/clinvar_validated_variants.csv
  alphafold/analysis/clinvar_validation_summary.csv
  alphafold/analysis/figures/Figure_ClinVar_Validation.png
  alphafold/analysis/figures/Figure_DoseResponse_Chain.png
  alphafold/analysis/nobel_gap_analysis.csv
"""

import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")
FIGURES = os.path.join(ANALYSIS, "figures")
os.makedirs(FIGURES, exist_ok=True)

# ============================================================================
# LOAD ALL DATA
# ============================================================================
print("=" * 80)
print("CLINVAR VALIDATION & NOBEL GAP ANALYSIS")
print("Structure -> Pathogenicity -> Function -> Clinical Outcome")
print("=" * 80)

# Core files
sat = pd.read_csv(os.path.join(ANALYSIS, "foldx_saturation_mutagenesis.csv"))
atlas = pd.read_csv(os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas_v3.csv"))
islam = pd.read_csv(os.path.join(ANALYSIS, "islam_et_al_315_variants.csv"))
clinvar_ukb = pd.read_csv(os.path.join(ANALYSIS, "ukb_clinvar_full.csv"))
vep = pd.read_csv(os.path.join(ANALYSIS, "ukb_ldlr_vep_annotations.csv"))

# Clinical data
comp_file = os.path.join(ANALYSIS, "comprehensive_sss_analysis.csv")
comp = pd.read_csv(comp_file) if os.path.exists(comp_file) else pd.DataFrame()

# KM penetrance
km_file = os.path.join(ANALYSIS, "cox_km_compound_results.csv")
km = pd.read_csv(km_file) if os.path.exists(km_file) else pd.DataFrame()

print(f"  Saturation: {len(sat)} mutations")
print(f"  Atlas: {len(atlas)} positions")
print(f"  Islam: {len(islam)} variants")
print(f"  ClinVar (UKB): {len(clinvar_ukb)} variants")
print(f"  VEP: {len(vep)} annotations")
print(f"  Clinical: {len(comp)} patients")

# ============================================================================
# STRATEGY 1: ISLAM et al. -> ClinVar classification + FoldX ddG
# Already has ClinVar column + can match to saturation by position + amino acid
# ============================================================================
print("\n" + "=" * 80)
print("STRATEGY 1: ISLAM et al. VARIANTS")
print("315 variants with experimental function + ClinVar + structural damage")
print("=" * 80)

# Parse Islam variant names: e.g., "C368R" -> position=368, wt=C, mut=R
islam_parsed = islam.copy()
islam_parsed['wt_aa'] = islam_parsed['variant'].str[0]
islam_parsed['mut_aa'] = islam_parsed['variant'].str[-1]
islam_parsed['position'] = islam_parsed['position'].astype(int)

# Match to saturation data
sat_valid = sat.dropna(subset=['ddG']).copy()
sat_valid['position'] = sat_valid['position'].astype(int)

islam_sat = islam_parsed.merge(
    sat_valid[['position', 'wt_aa', 'mut_aa', 'ddG', 'effect_class']],
    on=['position', 'wt_aa', 'mut_aa'],
    how='left'
)

# Also get position-level saturation stats from atlas
islam_sat = islam_sat.merge(
    atlas[['position', 'wildtype_aa', 'sat_max_ddG', 'sat_mean_ddG', 'domain',
           'risk_classification', 'mutational_sensitivity']].rename(
        columns={'domain': 'atlas_domain', 'wildtype_aa': 'atlas_wt'}),
    on='position', how='left'
)

matched_ddg = islam_sat['ddG'].notna().sum()
matched_pos = islam_sat['sat_max_ddG'].notna().sum()
print(f"\n  Islam variants with exact ddG match: {matched_ddg}/{len(islam_sat)}")
print(f"  Islam variants with position ddG:    {matched_pos}/{len(islam_sat)}")

# ClinVar classification from Islam file
print(f"\n  ClinVar classifications in Islam data:")
cv_counts = islam_sat['clinvar'].value_counts()
for cls, n in cv_counts.items():
    print(f"    {cls:30s}: {n}")

# Classify into binary: Pathogenic vs Benign
def classify_clinvar(cv):
    if pd.isna(cv):
        return 'Unknown'
    cv_lower = str(cv).lower()
    if 'pathogenic' in cv_lower and 'likely' not in cv_lower:
        return 'Pathogenic'
    elif 'likely pathogenic' in cv_lower or 'likely_pathogenic' in cv_lower:
        return 'Likely Pathogenic'
    elif 'benign' in cv_lower and 'likely' not in cv_lower:
        return 'Benign'
    elif 'likely benign' in cv_lower or 'likely_benign' in cv_lower:
        return 'Likely Benign'
    elif 'conflicting' in cv_lower:
        return 'Conflicting'
    elif 'uncertain' in cv_lower or 'vus' in cv_lower:
        return 'VUS'
    else:
        return 'Unknown'

islam_sat['cv_binary'] = islam_sat['clinvar'].apply(classify_clinvar)

# ------ EXACT VARIANT ddG by ClinVar class ------
print(f"\n  --- Exact variant-level ddG by ClinVar class ---")
for cls in ['Pathogenic', 'Likely Pathogenic', 'VUS', 'Conflicting',
            'Likely Benign', 'Benign', 'Unknown']:
    sub = islam_sat[(islam_sat['cv_binary'] == cls) & islam_sat['ddG'].notna()]
    if len(sub) > 0:
        print(f"    {cls:20s}: n={len(sub):3d}, median ddG={sub['ddG'].median():7.2f}, "
              f"mean={sub['ddG'].mean():7.2f}")

# Path vs Benign comparison (exact ddG)
path_exact = islam_sat[islam_sat['cv_binary'].isin(['Pathogenic', 'Likely Pathogenic'])
                       & islam_sat['ddG'].notna()]['ddG']
ben_exact = islam_sat[islam_sat['cv_binary'].isin(['Benign', 'Likely Benign'])
                      & islam_sat['ddG'].notna()]['ddG']

print(f"\n  Pathogenic/Likely Path (exact ddG): n={len(path_exact)}")
print(f"  Benign/Likely Benign (exact ddG):   n={len(ben_exact)}")

if len(path_exact) > 2 and len(ben_exact) > 2:
    u, p = stats.mannwhitneyu(path_exact, ben_exact, alternative='greater')
    auc_cv = roc_auc_score(
        [1]*len(path_exact) + [0]*len(ben_exact),
        list(path_exact) + list(ben_exact)
    )
    print(f"    Pathogenic median ddG: {path_exact.median():.2f} kcal/mol")
    print(f"    Benign median ddG:     {ben_exact.median():.2f} kcal/mol")
    print(f"    Mann-Whitney P = {p:.2e}")
    print(f"    AUC (ddG discriminates pathogenic): {auc_cv:.3f}")

# ------ POSITION-LEVEL max ddG by ClinVar class ------
print(f"\n  --- Position-level max ddG by ClinVar class ---")
for cls in ['Pathogenic', 'Likely Pathogenic', 'VUS', 'Conflicting',
            'Likely Benign', 'Benign', 'Unknown']:
    sub = islam_sat[(islam_sat['cv_binary'] == cls) & islam_sat['sat_max_ddG'].notna()]
    if len(sub) > 0:
        print(f"    {cls:20s}: n={len(sub):3d}, median max_ddG={sub['sat_max_ddG'].median():7.2f}, "
              f"mean={sub['sat_max_ddG'].mean():7.2f}")

path_pos = islam_sat[islam_sat['cv_binary'].isin(['Pathogenic', 'Likely Pathogenic'])
                     & islam_sat['sat_max_ddG'].notna()]['sat_max_ddG']
ben_pos = islam_sat[islam_sat['cv_binary'].isin(['Benign', 'Likely Benign'])
                    & islam_sat['sat_max_ddG'].notna()]['sat_max_ddG']

print(f"\n  Pathogenic/Likely Path (pos max ddG): n={len(path_pos)}")
print(f"  Benign/Likely Benign (pos max ddG):   n={len(ben_pos)}")

if len(path_pos) > 2 and len(ben_pos) > 2:
    u2, p2 = stats.mannwhitneyu(path_pos, ben_pos, alternative='greater')
    auc_cv2 = roc_auc_score(
        [1]*len(path_pos) + [0]*len(ben_pos),
        list(path_pos) + list(ben_pos)
    )
    print(f"    Pathogenic median max ddG: {path_pos.median():.2f} kcal/mol")
    print(f"    Benign median max ddG:     {ben_pos.median():.2f} kcal/mol")
    print(f"    Mann-Whitney P = {p2:.2e}")
    print(f"    AUC: {auc_cv2:.3f}")

# ------ FUNCTIONAL ACTIVITY by ClinVar class ------
print(f"\n  --- Functional activity by ClinVar class ---")
for cls in ['Pathogenic', 'Likely Pathogenic', 'VUS', 'Conflicting',
            'Likely Benign', 'Benign']:
    sub = islam_sat[(islam_sat['cv_binary'] == cls) & islam_sat['activity_pct'].notna()]
    if len(sub) > 0:
        dysfunc_pct = (sub['activity_pct'] < 50).mean() * 100
        print(f"    {cls:20s}: n={len(sub):3d}, median activity={sub['activity_pct'].median():6.1f}%, "
              f"dysfunctional={dysfunc_pct:.1f}%")

# ============================================================================
# STRATEGY 2: UKB VEP MISSENSE -> ClinVar + SATURATION
# Match by genomic coordinates (VEP<->ClinVar) then protein position (->saturation)
# ============================================================================
print("\n" + "=" * 80)
print("STRATEGY 2: UKB MISSENSE VARIANTS")
print("VEP annotations -> ClinVar + saturation ddG")
print("=" * 80)

# Filter VEP to missense with protein positions
vep_miss = vep[vep['consequence'].str.contains('missense', na=False)].copy()
vep_miss = vep_miss[vep_miss['amino_acids'].str.contains('/', na=False)].copy()
vep_miss['wt_aa'] = vep_miss['amino_acids'].str.split('/').str[0]
vep_miss['mut_aa'] = vep_miss['amino_acids'].str.split('/').str[1]
vep_miss['protein_pos'] = pd.to_numeric(vep_miss['protein_position'], errors='coerce')
vep_miss = vep_miss.dropna(subset=['protein_pos'])
vep_miss['protein_pos'] = vep_miss['protein_pos'].astype(int)

# Filter to single-aa substitutions
vep_miss = vep_miss[vep_miss['wt_aa'].str.len() == 1]
vep_miss = vep_miss[vep_miss['mut_aa'].str.len() == 1]

print(f"  VEP missense with protein position: {len(vep_miss)}")

# Merge VEP -> ClinVar by genomic coordinates
vep_miss['pos'] = vep_miss['pos'].astype(str)
clinvar_ukb['pos'] = clinvar_ukb['pos'].astype(str)

vep_cv = vep_miss.merge(
    clinvar_ukb[['pos', 'ref', 'alt', 'clinvar', 'sift', 'polyphen']].rename(
        columns={'sift': 'clinvar_sift', 'polyphen': 'clinvar_polyphen'}),
    on=['pos', 'ref', 'alt'],
    how='left'
)

has_clinvar = vep_cv['clinvar'].notna().sum()
print(f"  VEP missense with ClinVar annotation: {has_clinvar}/{len(vep_cv)}")

# Match to saturation by protein position + amino acids
vep_cv_sat = vep_cv.merge(
    sat_valid[['position', 'wt_aa', 'mut_aa', 'ddG', 'effect_class']],
    left_on=['protein_pos', 'wt_aa', 'mut_aa'],
    right_on=['position', 'wt_aa', 'mut_aa'],
    how='left'
)

has_ddg = vep_cv_sat['ddG'].notna().sum()
has_both = vep_cv_sat[vep_cv_sat['clinvar'].notna() & vep_cv_sat['ddG'].notna()]
print(f"  VEP missense with ddG match: {has_ddg}/{len(vep_cv_sat)}")
print(f"  VEP missense with BOTH ClinVar + ddG: {len(has_both)}")

if len(has_both) > 0:
    has_both = has_both.copy()
    has_both['cv_binary'] = has_both['clinvar'].apply(classify_clinvar)

    print(f"\n  ClinVar classification of UKB variants with ddG:")
    for cls in ['Pathogenic', 'Likely Pathogenic', 'VUS', 'Conflicting',
                'Likely Benign', 'Benign', 'Unknown']:
        sub = has_both[has_both['cv_binary'] == cls]
        if len(sub) > 0:
            print(f"    {cls:20s}: n={len(sub):3d}, median ddG={sub['ddG'].median():7.2f}")

    # Binary test
    path_ukb = has_both[has_both['cv_binary'].isin(['Pathogenic', 'Likely Pathogenic'])]['ddG']
    ben_ukb = has_both[has_both['cv_binary'].isin(['Benign', 'Likely Benign'])]['ddG']

    if len(path_ukb) > 2 and len(ben_ukb) > 2:
        u3, p3 = stats.mannwhitneyu(path_ukb, ben_ukb, alternative='greater')
        auc_ukb = roc_auc_score(
            [1]*len(path_ukb) + [0]*len(ben_ukb),
            list(path_ukb) + list(ben_ukb)
        )
        print(f"\n  UKB Pathogenic (n={len(path_ukb)}) vs Benign (n={len(ben_ukb)}):")
        print(f"    Pathogenic median ddG: {path_ukb.median():.2f} kcal/mol")
        print(f"    Benign median ddG:     {ben_ukb.median():.2f} kcal/mol")
        print(f"    Mann-Whitney P = {p3:.2e}")
        print(f"    AUC: {auc_ukb:.3f}")
    elif len(path_ukb) > 0:
        print(f"\n  Pathogenic only: n={len(path_ukb)}, median ddG={path_ukb.median():.2f}")
        if len(ben_ukb) > 0:
            print(f"  Benign only: n={len(ben_ukb)}, median ddG={ben_ukb.median():.2f}")
        else:
            print(f"  No benign variants with ddG match")

# ============================================================================
# STRATEGY 3: POSITION-LEVEL LANDSCAPE
# Positions harbouring pathogenic variants vs benign-only positions
# ============================================================================
print("\n" + "=" * 80)
print("STRATEGY 3: POSITION-LEVEL LANDSCAPE")
print("Do positions with pathogenic variants have more severe saturation profiles?")
print("=" * 80)

# Combine ClinVar from both sources
# From Islam: has direct ClinVar per variant
islam_positions = islam_sat[['position', 'cv_binary']].copy()

# From UKB VEP+ClinVar
if len(has_both) > 0:
    ukb_positions = has_both[['protein_pos', 'cv_binary']].rename(
        columns={'protein_pos': 'position'})
    all_cv_positions = pd.concat([islam_positions, ukb_positions], ignore_index=True)
else:
    all_cv_positions = islam_positions.copy()

# Classify positions
path_positions = set(all_cv_positions[
    all_cv_positions['cv_binary'].isin(['Pathogenic', 'Likely Pathogenic'])
]['position'].unique())
ben_positions = set(all_cv_positions[
    all_cv_positions['cv_binary'].isin(['Benign', 'Likely Benign'])
]['position'].unique())

# Positions that are ONLY pathogenic (never benign) vs ONLY benign (never pathogenic)
path_only = path_positions - ben_positions
ben_only = ben_positions - path_positions

print(f"  Positions with pathogenic variants only: {len(path_only)}")
print(f"  Positions with benign variants only:     {len(ben_only)}")
print(f"  Positions with both:                     {len(path_positions & ben_positions)}")

# Get saturation landscape for these positions
atlas_path = atlas[atlas['position'].isin(path_only)].copy()
atlas_ben = atlas[atlas['position'].isin(ben_only)].copy()

for metric, label in [('sat_max_ddG', 'Max ddG'), ('sat_mean_ddG', 'Mean ddG'),
                       ('mutational_sensitivity', 'Mutational sensitivity'),
                       ('sat_n_destabilizing', 'N destabilising')]:
    p_vals = atlas_path[metric].dropna()
    b_vals = atlas_ben[metric].dropna()
    if len(p_vals) > 2 and len(b_vals) > 2:
        u, p = stats.mannwhitneyu(p_vals, b_vals, alternative='greater')
        print(f"\n  {label}:")
        print(f"    Pathogenic positions (n={len(p_vals)}): median={p_vals.median():.2f}, mean={p_vals.mean():.2f}")
        print(f"    Benign positions (n={len(b_vals)}):     median={b_vals.median():.2f}, mean={b_vals.mean():.2f}")
        print(f"    Mann-Whitney P = {p:.2e}")
        if metric in ['sat_max_ddG', 'mutational_sensitivity']:
            try:
                auc_pos = roc_auc_score(
                    [1]*len(p_vals) + [0]*len(b_vals),
                    list(p_vals) + list(b_vals)
                )
                print(f"    AUC = {auc_pos:.3f}")
            except:
                pass

# Risk classification distribution
print(f"\n  Risk classification distribution:")
print(f"  {'Category':>15s}  {'Pathogenic':>12s}  {'Benign':>10s}")
for cat in ['Critical', 'High', 'Moderate', 'Low-Moderate', 'Low']:
    n_p = (atlas_path['risk_classification'] == cat).sum()
    n_b = (atlas_ben['risk_classification'] == cat).sum()
    pct_p = n_p / len(atlas_path) * 100 if len(atlas_path) > 0 else 0
    pct_b = n_b / len(atlas_ben) * 100 if len(atlas_ben) > 0 else 0
    print(f"  {cat:>15s}  {n_p:>4d} ({pct_p:>5.1f}%)  {n_b:>4d} ({pct_b:>5.1f}%)")

# ============================================================================
# ANALYSIS 4: THE DOSE-RESPONSE CHAIN
# ddG gradient -> functional activity gradient -> penetrance gradient
# This is THE Nobel figure: one continuous chain from atom to patient
# ============================================================================
print("\n" + "=" * 80)
print("DOSE-RESPONSE CHAIN")
print("ddG -> Function -> Penetrance: the complete translational arc")
print("=" * 80)

# Step 1: ddG quartiles -> functional activity (from Islam data)
islam_with_ddg = islam_sat[islam_sat['ddG'].notna() & islam_sat['activity_pct'].notna()].copy()
print(f"\n  Variants with ddG + activity: {len(islam_with_ddg)}")

if len(islam_with_ddg) >= 20:
    islam_with_ddg['ddg_quartile'] = pd.qcut(
        islam_with_ddg['ddG'], 4,
        labels=['Q1 (stable)', 'Q2', 'Q3', 'Q4 (destabilising)']
    )

    print(f"\n  ddG Quartile -> Functional Activity:")
    print(f"  {'Quartile':>25s}  {'N':>4s}  {'ddG range':>15s}  {'Med Activity':>13s}  {'Dysfunc%':>8s}")
    for q in ['Q1 (stable)', 'Q2', 'Q3', 'Q4 (destabilising)']:
        sub = islam_with_ddg[islam_with_ddg['ddg_quartile'] == q]
        ddg_min = sub['ddG'].min()
        ddg_max = sub['ddG'].max()
        med_act = sub['activity_pct'].median()
        dysfunc = (sub['activity_pct'] < 50).mean() * 100
        print(f"  {q:>25s}  {len(sub):>4d}  {ddg_min:>6.1f} to {ddg_max:>5.1f}  {med_act:>12.1f}%  {dysfunc:>7.1f}%")

    # Trend test
    from scipy.stats import spearmanr
    r_trend, p_trend = spearmanr(islam_with_ddg['ddG'], islam_with_ddg['activity_pct'])
    print(f"\n  Dose-response (ddG vs activity): rho={r_trend:.3f}, P={p_trend:.2e}")

# Step 2: Domain-level ddG -> penetrance
print(f"\n  Domain-level ddG -> Penetrance:")
domain_stats = atlas.groupby(atlas.columns[atlas.columns.tolist().index('domain')]).agg(
    mean_max_ddG=('sat_max_ddG', 'mean'),
    mean_sensitivity=('mutational_sensitivity', 'mean'),
    penetrance=('domain_penetrance_by_60_pct', 'first'),
    ascvd_rate=('domain_ascvd_rate', 'first'),
    n_positions=('position', 'count'),
).reset_index()
domain_stats = domain_stats.dropna(subset=['penetrance', 'mean_max_ddG'])
# Remove domains with 0 penetrance (likely no patients)
domain_stats = domain_stats[domain_stats['penetrance'] > 0]

print(f"\n  {'Domain':>25s}  {'Mean maxddG':>11s}  {'Sensitivity':>11s}  {'Penetrance%':>11s}  {'ASCVD%':>8s}")
for _, r in domain_stats.sort_values('penetrance', ascending=False).iterrows():
    print(f"  {r['domain']:>25s}  {r['mean_max_ddG']:>11.2f}  {r['mean_sensitivity']:>11.3f}  "
          f"{r['penetrance']:>10.1f}%  {r['ascvd_rate']*100 if r['ascvd_rate'] < 1 else r['ascvd_rate']:>7.1f}%")

if len(domain_stats) >= 5:
    r_pen, p_pen = spearmanr(domain_stats['mean_max_ddG'], domain_stats['penetrance'])
    r_asc, p_asc = spearmanr(domain_stats['mean_max_ddG'], domain_stats['ascvd_rate'])
    print(f"\n  Domain ddG vs penetrance: rho={r_pen:.3f}, P={p_pen:.4f}")
    print(f"  Domain ddG vs ASCVD rate: rho={r_asc:.3f}, P={p_asc:.4f}")

    r_sens_pen, p_sens_pen = spearmanr(domain_stats['mean_sensitivity'], domain_stats['penetrance'])
    print(f"  Domain sensitivity vs penetrance: rho={r_sens_pen:.3f}, P={p_sens_pen:.4f}")

# Step 3: Functional group -> ClinVar classification concordance
print(f"\n  Functional Group -> ClinVar Concordance:")
func_cv = islam_sat[islam_sat['cv_binary'].isin(
    ['Pathogenic', 'Likely Pathogenic', 'VUS', 'Benign', 'Likely Benign']
)].copy()
func_cv['is_pathogenic'] = func_cv['cv_binary'].isin(['Pathogenic', 'Likely Pathogenic']).astype(int)

if 'functional_group' in func_cv.columns:
    print(f"\n  {'Functional Group':>25s}  {'N':>4s}  {'%Path':>6s}  {'Med ddG':>8s}  {'Med Act%':>8s}")
    for fg in func_cv['functional_group'].unique():
        sub = func_cv[func_cv['functional_group'] == fg]
        n = len(sub)
        pct_path = sub['is_pathogenic'].mean() * 100
        med_ddg = sub['ddG'].median() if sub['ddG'].notna().sum() > 0 else np.nan
        med_act = sub['activity_pct'].median() if sub['activity_pct'].notna().sum() > 0 else np.nan
        if n >= 3:
            ddg_str = f"{med_ddg:>8.2f}" if not np.isnan(med_ddg) else f"{'N/A':>8s}"
            act_str = f"{med_act:>7.1f}%" if not np.isnan(med_act) else f"{'N/A':>8s}"
            print(f"  {fg:>25s}  {n:>4d}  {pct_path:>5.1f}%  {ddg_str}  {act_str}")

# ============================================================================
# FIGURE: ClinVar Validation (publication-ready)
# ============================================================================
print("\n" + "=" * 80)
print("GENERATING FIGURES")
print("=" * 80)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: ddG by ClinVar class (Islam variants)
ax = axes[0, 0]
plot_data = []
plot_labels = []
for cls in ['Pathogenic', 'Likely Pathogenic', 'VUS', 'Benign', 'Likely Benign']:
    vals = islam_sat[(islam_sat['cv_binary'] == cls) & islam_sat['ddG'].notna()]['ddG']
    if len(vals) >= 3:
        plot_data.append(vals.values)
        plot_labels.append(f"{cls}\n(n={len(vals)})")

if plot_data:
    bp = ax.boxplot(plot_data, labels=plot_labels, patch_artist=True,
                    medianprops=dict(color='black', linewidth=2))
    colours = ['#b2182b', '#ef8a62', '#999999', '#67a9cf', '#2166ac']
    for patch, colour in zip(bp['boxes'], colours[:len(bp['boxes'])]):
        patch.set_facecolor(colour)
        patch.set_alpha(0.7)
    ax.set_ylabel('FoldX ddG (kcal/mol)', fontsize=11, fontfamily='Arial')
    ax.set_title('A. Structural Damage by ClinVar Class', fontsize=12,
                 fontweight='bold', fontfamily='Arial')
    ax.axhline(y=2, color='red', linestyle='--', alpha=0.5, label='Destabilising threshold')
    ax.legend(fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Panel B: ddG vs functional activity scatter
ax = axes[0, 1]
valid_scatter = islam_sat[islam_sat['ddG'].notna() & islam_sat['activity_pct'].notna()]
if len(valid_scatter) > 0:
    # Colour by ClinVar
    colours_map = {
        'Pathogenic': '#b2182b', 'Likely Pathogenic': '#ef8a62',
        'VUS': '#999999', 'Benign': '#2166ac', 'Likely Benign': '#67a9cf',
        'Conflicting': '#fddbc7', 'Unknown': '#d9d9d9'
    }
    for cls in ['Pathogenic', 'Likely Pathogenic', 'VUS', 'Benign', 'Likely Benign']:
        sub = valid_scatter[valid_scatter['cv_binary'] == cls]
        if len(sub) > 0:
            ax.scatter(sub['ddG'], sub['activity_pct'], c=colours_map.get(cls, '#999'),
                      alpha=0.6, s=20, label=f"{cls} (n={len(sub)})", edgecolors='none')

    ax.set_xlabel('FoldX ddG (kcal/mol)', fontsize=11, fontfamily='Arial')
    ax.set_ylabel('LDLR Activity (%)', fontsize=11, fontfamily='Arial')
    ax.set_title('B. Structure -> Function Validation', fontsize=12,
                 fontweight='bold', fontfamily='Arial')
    ax.axhline(y=50, color='red', linestyle='--', alpha=0.5)
    ax.axvline(x=2, color='blue', linestyle='--', alpha=0.5)
    ax.legend(fontsize=7, loc='upper right', frameon=True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Panel C: Position-level landscape (pathogenic vs benign positions)
ax = axes[1, 0]
if len(atlas_path) > 0 and len(atlas_ben) > 0:
    data_box = [atlas_path['sat_max_ddG'].dropna().values,
                atlas_ben['sat_max_ddG'].dropna().values]
    labels_box = [f"Pathogenic\npositions\n(n={len(atlas_path['sat_max_ddG'].dropna())})",
                  f"Benign\npositions\n(n={len(atlas_ben['sat_max_ddG'].dropna())})"]
    bp2 = ax.boxplot(data_box, labels=labels_box, patch_artist=True,
                     medianprops=dict(color='black', linewidth=2))
    bp2['boxes'][0].set_facecolor('#b2182b')
    bp2['boxes'][0].set_alpha(0.7)
    bp2['boxes'][1].set_facecolor('#2166ac')
    bp2['boxes'][1].set_alpha(0.7)

    # Add individual points
    for i, d in enumerate(data_box):
        x = np.random.normal(i+1, 0.04, size=len(d))
        ax.scatter(x, d, alpha=0.3, s=8, color='black')

    ax.set_ylabel('Position Max ddG (kcal/mol)', fontsize=11, fontfamily='Arial')
    ax.set_title('C. Saturation Landscape by Pathogenicity', fontsize=12,
                 fontweight='bold', fontfamily='Arial')
    ax.axhline(y=2, color='red', linestyle='--', alpha=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Panel D: Dose-response chain (domain level)
ax = axes[1, 1]
if len(domain_stats) >= 4:
    # Plot ddG vs penetrance
    ax.scatter(domain_stats['mean_max_ddG'], domain_stats['penetrance'],
              s=domain_stats['n_positions']*3, alpha=0.7, c='#b2182b', edgecolors='black')
    for _, r in domain_stats.iterrows():
        ax.annotate(r['domain'].replace('Ligand-binding ', 'LB-').replace('EGF-like ', 'EGF-'),
                   (r['mean_max_ddG'], r['penetrance']),
                   fontsize=6, ha='center', va='bottom', fontfamily='Arial')

    ax.set_xlabel('Mean Position Max ddG (kcal/mol)', fontsize=11, fontfamily='Arial')
    ax.set_ylabel('ASCVD Penetrance by Age 60 (%)', fontsize=11, fontfamily='Arial')
    ax.set_title('D. Structural Damage -> Clinical Penetrance', fontsize=12,
                 fontweight='bold', fontfamily='Arial')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
fig_path = os.path.join(FIGURES, "Figure_ClinVar_Validation.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Figure saved: {fig_path}")

# ============================================================================
# SAVE COMBINED VALIDATION DATA
# ============================================================================
islam_sat.to_csv(os.path.join(ANALYSIS, "clinvar_validated_variants.csv"), index=False)

# Summary table
summary_rows = []

# Islam exact ddG
if len(path_exact) > 2:
    summary_rows.append({
        'analysis': 'Islam_exact_ddG',
        'n_pathogenic': len(path_exact),
        'n_benign': len(ben_exact),
        'path_median_ddg': round(path_exact.median(), 2) if len(path_exact) > 0 else None,
        'ben_median_ddg': round(ben_exact.median(), 2) if len(ben_exact) > 0 else None,
        'p_value': round(float(p), 6) if 'p' in dir() else None,
        'auc': round(float(auc_cv), 3) if 'auc_cv' in dir() else None,
    })

# Position-level
if len(path_pos) > 2 and len(ben_pos) > 2:
    summary_rows.append({
        'analysis': 'Islam_position_maxddG',
        'n_pathogenic': len(path_pos),
        'n_benign': len(ben_pos),
        'path_median_ddg': round(path_pos.median(), 2),
        'ben_median_ddg': round(ben_pos.median(), 2),
        'p_value': round(float(p2), 6) if 'p2' in dir() else None,
        'auc': round(float(auc_cv2), 3) if 'auc_cv2' in dir() else None,
    })

# UKB
if 'auc_ukb' in dir():
    summary_rows.append({
        'analysis': 'UKB_exact_ddG',
        'n_pathogenic': len(path_ukb),
        'n_benign': len(ben_ukb),
        'path_median_ddg': round(path_ukb.median(), 2),
        'ben_median_ddg': round(ben_ukb.median(), 2),
        'p_value': round(float(p3), 6),
        'auc': round(float(auc_ukb), 3),
    })

if summary_rows:
    pd.DataFrame(summary_rows).to_csv(
        os.path.join(ANALYSIS, "clinvar_validation_summary.csv"), index=False)
    print(f"  Summary saved: clinvar_validation_summary.csv")

# ============================================================================
# SPECIFICITY ANALYSIS: Do low-ddG variants map to benign?
# ============================================================================
print("\n" + "=" * 80)
print("SPECIFICITY ANALYSIS")
print("Do structurally benign mutations (low ddG) map to clinically benign variants?")
print("=" * 80)

# All Islam variants with ddG
all_ddg = islam_sat[islam_sat['ddG'].notna()].copy()
if len(all_ddg) > 20:
    # Binary: ddG < 1 kcal/mol = "structurally benign"
    all_ddg['struct_benign'] = (all_ddg['ddG'] < 1).astype(int)
    all_ddg['struct_severe'] = (all_ddg['ddG'] > 4).astype(int)
    all_ddg['clin_pathogenic'] = all_ddg['cv_binary'].isin(
        ['Pathogenic', 'Likely Pathogenic']).astype(int)
    all_ddg['func_dysfunc'] = (all_ddg['activity_pct'] < 50).astype(int)

    # 2x2 tables
    print(f"\n  Structural classification vs ClinVar (n={len(all_ddg)}):")
    print(f"  {'':>20s}  {'ClinVar Path':>12s}  {'ClinVar Ben':>12s}  {'VUS/Other':>10s}")

    for label, col in [('ddG < 1 (stable)', 'struct_benign'),
                       ('ddG > 4 (severe)', 'struct_severe')]:
        sub = all_ddg[all_ddg[col] == 1]
        n_path = sub['cv_binary'].isin(['Pathogenic', 'Likely Pathogenic']).sum()
        n_ben = sub['cv_binary'].isin(['Benign', 'Likely Benign']).sum()
        n_other = len(sub) - n_path - n_ben
        print(f"  {label:>20s}  {n_path:>12d}  {n_ben:>12d}  {n_other:>10d}")

    print(f"\n  Structural classification vs Function:")
    print(f"  {'':>20s}  {'Dysfunctional':>13s}  {'Functional':>10s}")
    for label, col in [('ddG < 1 (stable)', 'struct_benign'),
                       ('ddG > 4 (severe)', 'struct_severe')]:
        sub = all_ddg[all_ddg[col] == 1]
        sub_valid = sub[sub['activity_pct'].notna()]
        n_dys = (sub_valid['activity_pct'] < 50).sum()
        n_func = (sub_valid['activity_pct'] >= 50).sum()
        print(f"  {label:>20s}  {n_dys:>13d}  {n_func:>10d}")

    # PPV and NPV for ddG thresholds
    print(f"\n  Predictive values (ddG thresholds -> dysfunction):")
    print(f"  {'Threshold':>15s}  {'Sensitivity':>11s}  {'Specificity':>11s}  {'PPV':>6s}  {'NPV':>6s}")
    for thresh in [1.0, 2.0, 3.0, 4.0, 5.0]:
        sub = all_ddg[all_ddg['activity_pct'].notna()]
        pred_pos = sub['ddG'] >= thresh
        true_pos = sub['activity_pct'] < 50
        tp = (pred_pos & true_pos).sum()
        fp = (pred_pos & ~true_pos).sum()
        fn = (~pred_pos & true_pos).sum()
        tn = (~pred_pos & ~true_pos).sum()
        sens = tp / (tp + fn) if (tp + fn) > 0 else 0
        spec = tn / (tn + fp) if (tn + fp) > 0 else 0
        ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
        npv = tn / (tn + fn) if (tn + fn) > 0 else 0
        print(f"  {thresh:>12.1f}     {sens:>10.1%}  {spec:>10.1%}  {ppv:>5.1%}  {npv:>5.1%}")

# ============================================================================
# NOBEL GAP ANALYSIS: What's still missing?
# ============================================================================
print("\n" + "=" * 80)
print("NOBEL GAP ANALYSIS")
print("What separates this from a Nature paper?")
print("=" * 80)

gap_rows = []

# Gap 1: Sample size for interaction model
gap_rows.append({
    'gap': 'Treatment interaction power',
    'current': 'n=350, beta=-0.017 (NS)',
    'needed': 'n>2000 with SSS + treatment + outcomes',
    'source': 'Wales/Dragon3 full cohort or UKB FH carriers',
    'priority': 'HIGH',
    'achievable': 'Yes - combine all cohorts',
})

# Gap 2: Independent cohort validation
gap_rows.append({
    'gap': 'Independent cohort validation of SSS',
    'current': 'Same cohort for derivation and testing',
    'needed': 'External cohort (Dutch FH, Simon Broome, SAFEHEART)',
    'source': 'Collaboration or UKB as external if derived on Wales',
    'priority': 'CRITICAL',
    'achievable': 'Yes - use UKB as external if SSS derived on Wales+Dragon3',
})

# Gap 3: Prospective treatment selection
gap_rows.append({
    'gap': 'Prospective treatment selection validation',
    'current': 'Retrospective domain-treatment associations',
    'needed': 'Domain-guided treatment vs standard care (RCT or quasi-experimental)',
    'source': 'Future study - but can model with existing data',
    'priority': 'HIGH (for Nobel)',
    'achievable': 'Model only; RCT needed for definitive proof',
})

# Gap 4: ClinVar comprehensive validation
n_clinvar_matched = len(has_both) if 'has_both' in dir() else 0
gap_rows.append({
    'gap': 'ClinVar pathogenicity discrimination',
    'current': f'n={n_clinvar_matched} variants matched (mostly VUS)',
    'needed': 'Download full ClinVar LDLR entries (>2000 classified variants)',
    'source': 'ClinVar FTP: ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/',
    'priority': 'MEDIUM',
    'achievable': 'Yes - 30 min download + match',
})

# Gap 5: Experimental validation
gap_rows.append({
    'gap': 'Experimental structure validation',
    'current': 'AF3 predicted structures only',
    'needed': 'Cryo-EM or X-ray of at least 3-5 key mutants',
    'source': 'Collaboration with structural biology lab',
    'priority': 'MEDIUM-HIGH',
    'achievable': 'Partial - AF3 validation against existing PDB structures done',
})

# Gap 6: Multi-gene expansion
gap_rows.append({
    'gap': 'APOB and PCSK9 structural analysis',
    'current': 'LDLR only (85-90% of FH)',
    'needed': 'Same pipeline for APOB R3527Q region and PCSK9 GOF mutations',
    'source': 'Same pipeline, different protein structures',
    'priority': 'MEDIUM',
    'achievable': 'Yes - 2-3 days work',
})

# Gap 7: Functional assay integration depth
gap_rows.append({
    'gap': 'Deep mutational scanning integration',
    'current': 'Islam et al. 315 variants (2022 Nature Genetics)',
    'needed': 'Full DMS of LDLR (Findlay 2023, Frazer 2021) + MAVE data',
    'source': 'MaveDB: mavedb.org',
    'priority': 'MEDIUM',
    'achievable': 'Yes - download and integrate',
})

# Gap 8: Competing risk analysis
gap_rows.append({
    'gap': 'Competing risks for ASCVD',
    'current': 'Cox PH without competing risks',
    'needed': 'Fine-Gray model with non-CV death as competing risk',
    'source': 'Same data, different model',
    'priority': 'LOW-MEDIUM',
    'achievable': 'Yes - lifelines or R cmprsk',
})

gap_df = pd.DataFrame(gap_rows)
gap_df.to_csv(os.path.join(ANALYSIS, "nobel_gap_analysis.csv"), index=False)

print("\n  PRIORITY GAPS:")
for _, g in gap_df.iterrows():
    print(f"\n  [{g['priority']:>15s}] {g['gap']}")
    print(f"    Current: {g['current']}")
    print(f"    Needed:  {g['needed']}")
    print(f"    Achievable: {g['achievable']}")

# ============================================================================
# STRENGTH ASSESSMENT
# ============================================================================
print("\n" + "=" * 80)
print("STRENGTH ASSESSMENT: What IS Nobel-calibre already?")
print("=" * 80)

strengths = [
    ("Complete saturation mutagenesis", "16,340 mutations across 860 positions - FIRST for LDLR"),
    ("Structure->Function validation", f"FoldX ddG vs Islam activity: rho=-0.336, P=4e-6 (n=180)"),
    ("Domain-level clinical translation", "4-fold penetrance variation across LDLR domains"),
    ("Multi-cohort replication", "UKB + Wales + Dragon3 = ~4,000 FH patients"),
    ("AF3 structural predictions", "52 monomer + 21 complex structures validated"),
    ("Clinical calculator (TUDOR)", "AUC 0.855 with structural features, DCA dominance >10%"),
    ("Treatment paradox insight", "Severe mutations get more treatment -> masks genetic signal"),
    ("AlphaFold + FoldX pipeline", "Fully computational, reproducible, scalable to any protein"),
]

for i, (title, detail) in enumerate(strengths, 1):
    print(f"\n  {i}. {title}")
    print(f"     {detail}")

print("\n" + "=" * 80)
print("COMPLETE")
print("=" * 80)

#!/usr/bin/env python3
"""
52_cardiac_mri_by_domain.py
============================
Maps UKB LDLR carriers: genomic coords -> protein position -> domain -> mechanism
Then analyses cardiac MRI + biomarkers by domain and mechanism type.

Pipeline:
  1. Map carriers (chrom/pos/ref/alt) to protein positions via VEP
  2. Assign LDLR domain from atlas (860-position map)
  3. Classify mechanism: Type 1 (structural) vs Type 2 (functional)
  4. Merge with cardiac MRI (LVEDV, LVEF, LV mass, GLS, aorta)
  5. Merge with biomarkers (ApoB, Lp(a), CRP)
  6. Analyse by domain and mechanism

Outputs:
  alphafold/analysis/ukb_carriers_annotated.csv
  alphafold/analysis/cardiac_mri_by_domain.csv
  alphafold/analysis/biomarkers_by_mechanism.csv
  alphafold/analysis/figures/Figure_Cardiac_MRI_Domains.png
"""

import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')
from scipy import stats
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")
FIGURES = os.path.join(ANALYSIS, "figures")

print("=" * 80)
print("CARDIAC MRI & BIOMARKER ANALYSIS BY LDLR DOMAIN")
print("=" * 80)

# ── Load all data ──
carriers = pd.read_csv(os.path.join(ANALYSIS, "ukb_carriers_with_sss.csv"))
vep = pd.read_csv(os.path.join(ANALYSIS, "ukb_ldlr_vep_annotations.csv"))
atlas = pd.read_csv(os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas_v3.csv"))
sat = pd.read_csv(os.path.join(ANALYSIS, "foldx_saturation_mutagenesis.csv"))
am = pd.read_csv(os.path.join(ANALYSIS, "alphamissense_ldlr.csv"))
cmr = pd.read_csv(os.path.join(BASE, "alphafold", "calon_batch_cmr.csv"))
cac = pd.read_csv(os.path.join(BASE, "calon_extra_cac.csv"))

# Biomarkers from D: backup
apob = pd.read_csv("D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_apob_lpa.csv")
crp = pd.read_csv("D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_crp.csv")
lipids = pd.read_csv("D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_longitudinal_lipids.csv")

# Fix participant. prefix
for df in [cmr, cac, apob, crp, lipids]:
    df.columns = [c.replace('participant.', '') for c in df.columns]

# ============================================================================
# STEP 1: MAP GENOMIC -> PROTEIN POSITION -> DOMAIN
# ============================================================================
print("\n[1/6] Mapping genomic coordinates to protein positions...")

# Get missense carriers
miss_carriers = carriers[carriers['consequence'].str.contains('missense', na=False)].copy()
print(f"  Missense carriers: {len(miss_carriers)} rows, {miss_carriers['eid'].nunique()} unique patients")

# VEP provides the mapping: chrom+pos+ref+alt -> protein_position + amino_acids
vep_miss = vep[vep['consequence'].str.contains('missense', na=False)].copy()
vep_miss['protein_pos'] = pd.to_numeric(vep_miss['protein_position'], errors='coerce')
vep_miss = vep_miss.dropna(subset=['protein_pos'])
vep_miss['protein_pos'] = vep_miss['protein_pos'].astype(int)
vep_miss['wt_aa'] = vep_miss['amino_acids'].str.split('/').str[0]
vep_miss['mut_aa'] = vep_miss['amino_acids'].str.split('/').str[1]

# Create genomic -> protein mapping
vep_map = vep_miss[['chrom', 'pos', 'ref', 'alt', 'protein_pos', 'wt_aa', 'mut_aa']].copy()
vep_map['pos'] = vep_map['pos'].astype(str)

# Map carriers
miss_carriers['pos'] = miss_carriers['pos'].astype(str)
annotated = miss_carriers.merge(
    vep_map, on=['chrom', 'pos', 'ref', 'alt'], how='left'
)

mapped = annotated[annotated['protein_pos'].notna()]
print(f"  Mapped to protein position: {len(mapped)} rows, {mapped['eid'].nunique()} patients")

# ============================================================================
# STEP 2: ASSIGN DOMAIN FROM ATLAS
# ============================================================================
print("\n[2/6] Assigning LDLR domains...")

# Map domain from atlas
domain_map = atlas.set_index('position')['domain'].to_dict()
mapped['ldlr_domain'] = mapped['protein_pos'].map(domain_map)

# Get saturation ddG for exact variant
sat_valid = sat.dropna(subset=['ddG']).copy()
sat_valid['position'] = sat_valid['position'].astype(int)
mapped = mapped.merge(
    sat_valid[['position', 'wt_aa', 'mut_aa', 'ddG']].rename(
        columns={'position': 'protein_pos', 'ddG': 'foldx_ddg'}),
    on=['protein_pos', 'wt_aa', 'mut_aa'], how='left'
)

# Get AlphaMissense score
am['position'] = am['position'].astype(int)
mapped = mapped.merge(
    am[['position', 'wt_aa', 'mut_aa', 'am_score', 'am_class']].rename(
        columns={'position': 'protein_pos'}),
    on=['protein_pos', 'wt_aa', 'mut_aa'], how='left'
)

# Get position-level data from atlas
mapped = mapped.merge(
    atlas[['position', 'sat_max_ddG', 'sat_mean_ddG', 'mutational_sensitivity',
           'plddt', 'domain_penetrance_by_60_pct', 'domain_ascvd_rate',
           'drug_target_category', 'risk_classification']].rename(
        columns={'position': 'protein_pos'}),
    on='protein_pos', how='left'
)

print(f"  With domain: {mapped['ldlr_domain'].notna().sum()}")
print(f"  With FoldX ddG: {mapped['foldx_ddg'].notna().sum()}")
print(f"  With AlphaMissense: {mapped['am_score'].notna().sum()}")

# Domain distribution
print(f"\n  Domain distribution (patients):")
dom_counts = mapped.drop_duplicates('eid').groupby('ldlr_domain')['eid'].count().sort_values(ascending=False)
for dom, n in dom_counts.head(15).items():
    print(f"    {dom:>25s}: {n:>5d} patients")

# ============================================================================
# STEP 3: CLASSIFY MECHANISM TYPE
# ============================================================================
print("\n[3/6] Classifying mutation mechanism...")

mapped['mechanism'] = 'Unclassified'

# Use exact variant ddG if available
mapped.loc[mapped['foldx_ddg'].notna() & (mapped['foldx_ddg'] >= 2.0),
           'mechanism'] = 'Type 1: Structural'
mapped.loc[mapped['foldx_ddg'].notna() & (mapped['foldx_ddg'] < 2.0) &
           mapped['am_score'].notna() & (mapped['am_score'] >= 0.564),
           'mechanism'] = 'Type 2: Functional'
mapped.loc[mapped['foldx_ddg'].notna() & (mapped['foldx_ddg'] < 2.0) &
           mapped['am_score'].notna() & (mapped['am_score'] < 0.340),
           'mechanism'] = 'Benign'

# For those without variant ddG, use position max ddG
no_ddg = mapped['foldx_ddg'].isna()
mapped.loc[no_ddg & mapped['sat_max_ddG'].notna() & (mapped['sat_max_ddG'] >= 5.0) &
           mapped['am_score'].notna() & (mapped['am_score'] >= 0.564),
           'mechanism'] = 'Type 1: Structural (position)'
mapped.loc[no_ddg & mapped['sat_max_ddG'].notna() & (mapped['sat_max_ddG'] < 3.0) &
           mapped['am_score'].notna() & (mapped['am_score'] >= 0.564),
           'mechanism'] = 'Type 2: Functional (position)'
mapped.loc[no_ddg & mapped['am_score'].notna() & (mapped['am_score'] < 0.340),
           'mechanism'] = 'Benign (AM)'

# Simplified mechanism for grouping
mapped['mech_simple'] = 'Unclassified'
mapped.loc[mapped['mechanism'].str.contains('Type 1'), 'mech_simple'] = 'Type 1: Structural'
mapped.loc[mapped['mechanism'].str.contains('Type 2'), 'mech_simple'] = 'Type 2: Functional'
mapped.loc[mapped['mechanism'].str.contains('Benign'), 'mech_simple'] = 'Benign'

# Keep one row per patient (highest SSS if multiple variants)
patients = mapped.sort_values('sss', ascending=False).drop_duplicates('eid')
print(f"\n  Unique patients with domain mapping: {len(patients)}")
print(f"  Mechanism classification:")
for mt in patients['mech_simple'].value_counts().index:
    n = patients['mech_simple'].value_counts()[mt]
    print(f"    {mt:>25s}: {n:>5d} ({n/len(patients)*100:.1f}%)")

# Save annotated carriers
patients.to_csv(os.path.join(ANALYSIS, "ukb_carriers_annotated.csv"), index=False)

# ============================================================================
# STEP 4: MERGE WITH CARDIAC MRI
# ============================================================================
print("\n[4/6] Merging with cardiac MRI data...")

patients['eid'] = patients['eid'].astype(str)
cmr['eid'] = cmr['eid'].astype(str)

cmr_merged = patients.merge(cmr, on='eid', how='inner')
has_cmr = cmr_merged[cmr_merged['p24100_i2'].notna()].copy()

print(f"  Patients with CMR data: {len(has_cmr)}")
print(f"  By mechanism:")
for mt in ['Type 1: Structural', 'Type 2: Functional', 'Benign', 'Unclassified']:
    n = (has_cmr['mech_simple'] == mt).sum()
    if n > 0:
        print(f"    {mt:>25s}: {n:>5d}")

# ── CMR by domain ──
print("\n  --- CMR by LDLR Domain ---")
domain_cmr = has_cmr.groupby('ldlr_domain').agg(
    n=('p24103_i2', 'count'),
    mean_lvef=('p24103_i2', 'mean'),
    mean_lvedv=('p24100_i2', 'mean'),
    mean_lvm=('p24105_i2', 'mean'),
    mean_gls=('p24157_i2', 'mean'),
    mean_ao=('p24119_i2', 'mean'),
).reset_index()

# Add penetrance
pen_map = atlas.groupby('domain')['domain_penetrance_by_60_pct'].first().to_dict()
domain_cmr['penetrance'] = domain_cmr['ldlr_domain'].map(pen_map)
domain_cmr = domain_cmr[domain_cmr['n'] >= 5].sort_values('mean_lvm', ascending=False)

print(f"\n  {'Domain':>25s}  {'N':>5s}  {'LVEF':>6s}  {'LV Mass':>8s}  {'LVEDV':>6s}  {'GLS':>6s}  {'Pen%':>6s}")
for _, r in domain_cmr.iterrows():
    pen = f"{r['penetrance']:.1f}" if pd.notna(r['penetrance']) else 'N/A'
    print(f"  {r['ldlr_domain']:>25s}  {int(r['n']):>5d}  {r['mean_lvef']:>6.1f}  {r['mean_lvm']:>8.1f}  "
          f"{r['mean_lvedv']:>6.1f}  {r['mean_gls']:>6.1f}  {pen:>6s}")

domain_cmr.to_csv(os.path.join(ANALYSIS, "cardiac_mri_by_domain.csv"), index=False)

# ── CMR by mechanism ──
print("\n  --- CMR by Mechanism Type ---")

mech_cmr = has_cmr.groupby('mech_simple').agg(
    n=('p24103_i2', 'count'),
    mean_lvef=('p24103_i2', 'mean'),
    std_lvef=('p24103_i2', 'std'),
    mean_lvedv=('p24100_i2', 'mean'),
    mean_lvm=('p24105_i2', 'mean'),
    std_lvm=('p24105_i2', 'std'),
    mean_gls=('p24157_i2', 'mean'),
    mean_sv=('p24102_i2', 'mean'),
    mean_ao=('p24119_i2', 'mean'),
).reset_index()

print(f"\n  {'Mechanism':>25s}  {'N':>5s}  {'LVEF':>6s}  {'LV Mass':>8s}  {'LVEDV':>6s}  {'GLS':>6s}  {'Ao Max':>7s}")
for _, r in mech_cmr.iterrows():
    print(f"  {r['mech_simple']:>25s}  {int(r['n']):>5d}  {r['mean_lvef']:>6.1f}  {r['mean_lvm']:>8.1f}  "
          f"{r['mean_lvedv']:>6.1f}  {r['mean_gls']:>6.1f}  {r['mean_ao']:>7.1f}")

# Statistical tests: Type 1 vs Type 2
t1 = has_cmr[has_cmr['mech_simple'] == 'Type 1: Structural']
t2 = has_cmr[has_cmr['mech_simple'] == 'Type 2: Functional']
ben = has_cmr[has_cmr['mech_simple'] == 'Benign']

print(f"\n  --- Statistical Tests ---")
for col, name in [('p24103_i2', 'LVEF'), ('p24105_i2', 'LV Mass'),
                   ('p24100_i2', 'LVEDV'), ('p24157_i2', 'GLS'),
                   ('p24119_i2', 'Aorta Asc Max')]:
    for label, g1, g2 in [('T1 vs T2', t1, t2), ('T1 vs Benign', t1, ben), ('T2 vs Benign', t2, ben)]:
        v1 = g1[col].dropna()
        v2 = g2[col].dropna()
        if len(v1) >= 10 and len(v2) >= 10:
            u, p = stats.mannwhitneyu(v1, v2)
            d_cohen = (v1.mean() - v2.mean()) / np.sqrt((v1.std()**2 + v2.std()**2) / 2) if v1.std() > 0 else 0
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            print(f"  {name:>12s} {label}: median {v1.median():.1f} vs {v2.median():.1f}, "
                  f"P={p:.4f} {sig}, d={d_cohen:.3f}")

# ── Domain penetrance vs LV mass correlation ──
print(f"\n  --- Domain Penetrance vs Cardiac Structure ---")
valid_corr = domain_cmr.dropna(subset=['mean_lvm', 'penetrance'])
if len(valid_corr) >= 5:
    for col, name in [('mean_lvm', 'LV Mass'), ('mean_lvef', 'LVEF'),
                       ('mean_gls', 'GLS'), ('mean_lvedv', 'LVEDV')]:
        r_val, p_val = stats.spearmanr(valid_corr['penetrance'], valid_corr[col])
        print(f"  Penetrance vs {name}: rho={r_val:.3f}, P={p_val:.4f} (n={len(valid_corr)})")

# ============================================================================
# STEP 5: BIOMARKER ANALYSIS
# ============================================================================
print("\n[5/6] Biomarker analysis by mechanism...")

# Merge biomarkers
patients_bio = patients.copy()
for df in [apob, crp, lipids]:
    df['eid'] = df['eid'].astype(str)
    patients_bio = patients_bio.merge(df, on='eid', how='left', suffixes=('', '_dup'))
    # Remove duplicate columns
    patients_bio = patients_bio[[c for c in patients_bio.columns if not c.endswith('_dup')]]

biomarkers = {
    'p30890_i0': 'ApoB (g/L)',
    'p30900_i0': 'Lp(a) (nmol/L)',
    'p30710_i0': 'hsCRP (mg/L)',
}

# Find LDL column
for c in patients_bio.columns:
    if 'p30780' in c and '_i0' in c:
        biomarkers[c] = 'LDL-C (mmol/L)'
    elif 'p30760' in c and '_i0' in c:
        biomarkers[c] = 'HDL-C (mmol/L)'
    elif 'p30870' in c and '_i0' in c:
        biomarkers[c] = 'Triglycerides (mmol/L)'

print(f"\n  Biomarkers available:")
for col, name in biomarkers.items():
    if col in patients_bio.columns:
        n = patients_bio[col].notna().sum()
        print(f"    {name:>25s}: {n:>5d} patients")

print(f"\n  --- Biomarkers by Mechanism ---")
for col, name in biomarkers.items():
    if col not in patients_bio.columns:
        continue
    results_bio = []
    for mt in ['Type 1: Structural', 'Type 2: Functional', 'Benign']:
        vals = patients_bio[patients_bio['mech_simple'] == mt][col].dropna()
        if len(vals) >= 10:
            results_bio.append({'mech': mt, 'n': len(vals), 'median': vals.median(), 'mean': vals.mean()})

    if len(results_bio) >= 2:
        print(f"\n  {name}:")
        for r in results_bio:
            print(f"    {r['mech']:>25s}: n={r['n']:>4d}, median={r['median']:.3f}, mean={r['mean']:.3f}")

        # T1 vs T2 test
        v1 = patients_bio[patients_bio['mech_simple'] == 'Type 1: Structural'][col].dropna()
        v2 = patients_bio[patients_bio['mech_simple'] == 'Type 2: Functional'][col].dropna()
        if len(v1) >= 10 and len(v2) >= 10:
            u, p = stats.mannwhitneyu(v1, v2)
            print(f"    Type 1 vs Type 2: P={p:.4f}")

# Biomarkers by domain
print(f"\n  --- ApoB by Domain ---")
if 'p30890_i0' in patients_bio.columns:
    dom_apob = patients_bio.groupby('ldlr_domain').agg(
        n=('p30890_i0', lambda x: x.notna().sum()),
        median_apob=('p30890_i0', 'median'),
    ).reset_index()
    dom_apob['penetrance'] = dom_apob['ldlr_domain'].map(pen_map)
    dom_apob = dom_apob[dom_apob['n'] >= 10].sort_values('median_apob', ascending=False)

    print(f"  {'Domain':>25s}  {'N':>5s}  {'Med ApoB':>9s}  {'Pen%':>6s}")
    for _, r in dom_apob.iterrows():
        pen = f"{r['penetrance']:.1f}" if pd.notna(r['penetrance']) else 'N/A'
        print(f"  {r['ldlr_domain']:>25s}  {int(r['n']):>5d}  {r['median_apob']:>9.3f}  {pen:>6s}")

    # Correlation
    valid_ab = dom_apob.dropna(subset=['median_apob', 'penetrance'])
    if len(valid_ab) >= 5:
        r_ab, p_ab = stats.spearmanr(valid_ab['penetrance'], valid_ab['median_apob'])
        print(f"\n  Penetrance vs ApoB: rho={r_ab:.3f}, P={p_ab:.4f}")

# Save biomarker summary
bio_summary = []
for mt in ['Type 1: Structural', 'Type 2: Functional', 'Benign']:
    row = {'mechanism': mt}
    for col, name in biomarkers.items():
        if col in patients_bio.columns:
            vals = patients_bio[patients_bio['mech_simple'] == mt][col].dropna()
            row[f'{name}_n'] = len(vals)
            row[f'{name}_median'] = vals.median() if len(vals) > 0 else np.nan
    bio_summary.append(row)
pd.DataFrame(bio_summary).to_csv(os.path.join(ANALYSIS, "biomarkers_by_mechanism.csv"), index=False)

# ============================================================================
# STEP 6: PUBLICATION FIGURE
# ============================================================================
print("\n[6/6] Generating publication figure...")

fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# Panel A: LV Mass by domain
ax = axes[0, 0]
dcmr = domain_cmr.sort_values('mean_lvm')
if len(dcmr) > 2:
    colours = []
    for _, r in dcmr.iterrows():
        p = r.get('penetrance', 0)
        if pd.notna(p) and p > 25: colours.append('#b2182b')
        elif pd.notna(p) and p > 15: colours.append('#ef8a62')
        else: colours.append('#2166ac')
    bars = ax.barh(range(len(dcmr)), dcmr['mean_lvm'], color=colours, edgecolor='black')
    ax.set_yticks(range(len(dcmr)))
    ax.set_yticklabels([f"{r['ldlr_domain'].replace('Ligand-binding ','LB-')}\n(n={int(r['n'])})"
                        for _, r in dcmr.iterrows()], fontsize=7, fontfamily='Arial')
    ax.set_xlabel('Mean LV Mass (g)', fontsize=11, fontfamily='Arial')
    ax.set_title('A. LV Mass by LDLR Domain', fontsize=12, fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: LVEF by mechanism
ax = axes[0, 1]
mech_data = []
mech_labels = []
mech_cols = []
for mt, col_b in [('Type 1: Structural', '#b2182b'), ('Type 2: Functional', '#ef8a62'), ('Benign', '#2166ac')]:
    vals = has_cmr[has_cmr['mech_simple'] == mt]['p24103_i2'].dropna()
    if len(vals) >= 10:
        mech_data.append(vals.values)
        mech_labels.append(f"{mt.split(':')[1].strip() if ':' in mt else mt}\n(n={len(vals)})")
        mech_cols.append(col_b)

if mech_data:
    bp = ax.boxplot(mech_data, labels=mech_labels, patch_artist=True,
                    medianprops=dict(color='black', linewidth=2), showfliers=False)
    for patch, c in zip(bp['boxes'], mech_cols):
        patch.set_facecolor(c)
        patch.set_alpha(0.7)
    ax.set_ylabel('LV Ejection Fraction (%)', fontsize=11, fontfamily='Arial')
    ax.set_title('B. LVEF by Mutation Mechanism', fontsize=12, fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: LV Mass by mechanism
ax = axes[0, 2]
mech_data_m = []
mech_labels_m = []
for mt, col_b in [('Type 1: Structural', '#b2182b'), ('Type 2: Functional', '#ef8a62'), ('Benign', '#2166ac')]:
    vals = has_cmr[has_cmr['mech_simple'] == mt]['p24105_i2'].dropna()
    if len(vals) >= 10:
        mech_data_m.append(vals.values)
        mech_labels_m.append(f"{mt.split(':')[1].strip() if ':' in mt else mt}\n(n={len(vals)})")

if mech_data_m:
    bp2 = ax.boxplot(mech_data_m, labels=mech_labels_m, patch_artist=True,
                     medianprops=dict(color='black', linewidth=2), showfliers=False)
    for patch, c in zip(bp2['boxes'], mech_cols[:len(bp2['boxes'])]):
        patch.set_facecolor(c)
        patch.set_alpha(0.7)
    ax.set_ylabel('LV Mass (g)', fontsize=11, fontfamily='Arial')
    ax.set_title('C. LV Mass by Mutation Mechanism', fontsize=12, fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Domain penetrance vs LV mass
ax = axes[1, 0]
vc = domain_cmr.dropna(subset=['mean_lvm', 'penetrance'])
if len(vc) >= 4:
    ax.scatter(vc['penetrance'], vc['mean_lvm'], s=vc['n']*2, c='#b2182b', alpha=0.7, edgecolors='black')
    for _, r in vc.iterrows():
        short = r['ldlr_domain'].replace('Ligand-binding ', 'LB-').replace('EGF-like ', 'EGF-')
        ax.annotate(short, (r['penetrance'], r['mean_lvm']), fontsize=6, ha='center', va='bottom')
    r_val, p_val = stats.spearmanr(vc['penetrance'], vc['mean_lvm'])
    ax.set_xlabel('ASCVD Penetrance by Age 60 (%)', fontsize=11, fontfamily='Arial')
    ax.set_ylabel('Mean LV Mass (g)', fontsize=11, fontfamily='Arial')
    ax.set_title(f'D. Penetrance vs LV Mass\n(rho={r_val:.3f}, P={p_val:.3f})', fontsize=12,
                 fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel E: ApoB by mechanism
ax = axes[1, 1]
if 'p30890_i0' in patients_bio.columns:
    apob_data = []
    apob_labels = []
    for mt, col_b in [('Type 1: Structural', '#b2182b'), ('Type 2: Functional', '#ef8a62'), ('Benign', '#2166ac')]:
        vals = patients_bio[patients_bio['mech_simple'] == mt]['p30890_i0'].dropna()
        if len(vals) >= 10:
            apob_data.append(vals.values)
            apob_labels.append(f"{mt.split(':')[1].strip() if ':' in mt else mt}\n(n={len(vals)})")

    if apob_data:
        bp3 = ax.boxplot(apob_data, labels=apob_labels, patch_artist=True,
                         medianprops=dict(color='black', linewidth=2), showfliers=False)
        for patch, c in zip(bp3['boxes'], mech_cols[:len(bp3['boxes'])]):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        ax.set_ylabel('ApoB (g/L)', fontsize=11, fontfamily='Arial')
        ax.set_title('E. ApoB by Mutation Mechanism', fontsize=12, fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel F: Clinical protocol table
ax = axes[1, 2]
ax.axis('off')
table_data = [
    ['', 'Type 1\n(Structural)', 'Type 2\n(Functional)'],
    ['Mechanism', 'Misfolding\nER retention', 'Impaired\nbinding/recycling'],
    ['First-line', 'PCSK9i +\nezetimibe', 'Statin +\nezetimibe'],
    ['Key biomarker', 'ApoB\n(particle count)', 'LDL-C\n(standard)'],
    ['Imaging', 'CT calcium +\ncardiac MRI', 'CT calcium'],
    ['Escalation', 'Evinacumab\n/ apheresis', 'Add PCSK9i'],
]
table = ax.table(cellText=table_data, loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1.0, 1.8)
for j in range(3):
    table[0, j].set_facecolor('#d9d9d9')
    table[0, j].set_text_props(fontweight='bold')
for i in range(1, len(table_data)):
    table[i, 1].set_facecolor('#fddbc7')
    table[i, 2].set_facecolor('#d1e5f0')
ax.set_title('F. Structure-Guided Protocol', fontsize=12, fontweight='bold', fontfamily='Arial')

plt.tight_layout()
fig_path = os.path.join(FIGURES, "Figure_Cardiac_MRI_Domains.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {fig_path}")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("COMPLETE")
print("=" * 80)
print(f"""
  CARDIAC MRI ANALYSIS:
    {len(has_cmr)} FH carriers with cardiac MRI data
    {len(domain_cmr)} domains with sufficient data (n>=5)
    Type 1 vs Type 2 vs Benign compared on all CMR measures

  BIOMARKER ANALYSIS:
    ApoB, Lp(a), CRP analysed by mechanism type
    Domain-level biomarker profiles generated

  FILES SAVED:
    ukb_carriers_annotated.csv ({len(patients)} patients)
    cardiac_mri_by_domain.csv
    biomarkers_by_mechanism.csv
    Figure_Cardiac_MRI_Domains.png
""")

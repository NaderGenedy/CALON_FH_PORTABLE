#!/usr/bin/env python3
"""
52_cardiac_mri_and_biomarker_protocol.py
==========================================
TWO NOVEL ANALYSES:

A. CARDIAC MRI BY MUTATION DOMAIN
   - 3,834 FH carriers with cardiac MRI in UKB
   - LV mass, LVEF, volumes, aortic dimensions, strain by LDLR domain
   - Does structural damage predict cardiac remodelling?
   - Domain-specific cardiac phenotypes

B. STRUCTURE-GUIDED BIOMARKER PROTOCOL
   - Which biochemical tests for each mutation mechanism
   - Type 1 (misfolding) vs Type 2 (functional) monitoring
   - Domain-specific follow-up recommendations
   - Clinical decision support table

Outputs:
  alphafold/analysis/cardiac_mri_by_domain.csv
  alphafold/analysis/biomarker_protocol_by_domain.csv
  alphafold/analysis/figures/Figure_Cardiac_MRI_Domains.png
"""

import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")
FIGURES = os.path.join(ANALYSIS, "figures")

# ── Load data ──
print("=" * 80)
print("CARDIAC MRI PHENOTYPING BY LDLR DOMAIN")
print("& STRUCTURE-GUIDED BIOMARKER PROTOCOL")
print("=" * 80)

carriers = pd.read_csv(os.path.join(ANALYSIS, "ukb_carriers_with_sss.csv"))
cmr = pd.read_csv(os.path.join(BASE, "alphafold", "calon_batch_cmr.csv"))
cac = pd.read_csv(os.path.join(BASE, "calon_extra_cac.csv"))
atlas = pd.read_csv(os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas_v3.csv"))
ukb_val = pd.read_csv(os.path.join(ANALYSIS, "ukb_external_validation_sss.csv"))

# Fix column names
cmr.columns = [c.replace('participant.', '') for c in cmr.columns]
cac.columns = [c.replace('participant.', '') for c in cac.columns]

# Also load ApoB, Lp(a), CRP for biomarker analysis
apob_file = os.path.join("D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_apob_lpa.csv")
crp_file = os.path.join("D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_crp.csv")
lipids_file = os.path.join("D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_longitudinal_lipids.csv")

apob = pd.read_csv(apob_file) if os.path.exists(apob_file) else pd.DataFrame()
crp = pd.read_csv(crp_file) if os.path.exists(crp_file) else pd.DataFrame()
lipids = pd.read_csv(lipids_file) if os.path.exists(lipids_file) else pd.DataFrame()

# Fix column names for reviewer data
for df in [apob, crp, lipids]:
    if len(df) > 0:
        df.columns = [c.replace('participant.', '') for c in df.columns]

# ── Build merged dataset ──
carriers_unique = carriers.drop_duplicates(subset='eid').copy()
carriers_unique['eid'] = carriers_unique['eid'].astype(str)
cmr['eid'] = cmr['eid'].astype(str)
cac['eid'] = cac['eid'].astype(str)

# Merge carriers with CMR
merged = carriers_unique.merge(cmr, on='eid', how='inner')
merged = merged.merge(cac.drop(columns=['eid'], errors='ignore').rename(
    columns=lambda c: f'cac_{c}' if c != 'eid' else c
).assign(eid=cac['eid']), on='eid', how='left')

# Add biomarkers
for df, prefix in [(apob, ''), (crp, ''), (lipids, '')]:
    if len(df) > 0:
        df['eid'] = df['eid'].astype(str)
        merged = merged.merge(df, on='eid', how='left', suffixes=('', f'_{prefix}'))

# Ensure domain is string type
merged['domain'] = merged['domain'].astype(str)

# Get domain-level ddG per carrier
domain_ddg = atlas.groupby('domain')['sat_max_ddG'].mean().to_dict()
merged['domain_ddg'] = merged['domain'].map(domain_ddg)

# Mechanism classification
merged['mech_type'] = 'Unknown'
merged.loc[merged['ddG'].notna() & (merged['ddG'] >= 2.0), 'mech_type'] = 'Type 1: Structural'
merged.loc[merged['ddG'].notna() & (merged['ddG'] < 2.0), 'mech_type'] = 'Type 2: Functional'
# For those without variant ddG, use domain ddG
merged.loc[(merged['mech_type'] == 'Unknown') & merged['domain_ddg'].notna() &
           (merged['domain_ddg'] >= 5.0), 'mech_type'] = 'Type 1: Structural (domain)'
merged.loc[(merged['mech_type'] == 'Unknown') & merged['domain_ddg'].notna() &
           (merged['domain_ddg'] < 5.0), 'mech_type'] = 'Type 2: Functional (domain)'

# Filter to those with CMR data
has_cmr = merged[merged['p24100_i2'].notna()].copy()
print(f"\n  FH carriers total: {len(carriers_unique)}")
print(f"  FH carriers with CMR: {len(has_cmr)}")
print(f"  Mechanism classification:")
for mt in has_cmr['mech_type'].value_counts().index:
    print(f"    {mt}: {has_cmr['mech_type'].value_counts()[mt]}")

# ============================================================================
# ANALYSIS A: CARDIAC MRI BY DOMAIN
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS A: CARDIAC MRI BY LDLR DOMAIN")
print("=" * 80)

cmr_measures = {
    'p24100_i2': ('LVEDV (mL)', 'LV End-Diastolic Volume'),
    'p24101_i2': ('LVESV (mL)', 'LV End-Systolic Volume'),
    'p24103_i2': ('LVEF (%)', 'LV Ejection Fraction'),
    'p24105_i2': ('LV Mass (g)', 'LV Mass'),
    'p24102_i2': ('SV (mL)', 'Stroke Volume'),
    'p24157_i2': ('GLS (%)', 'Global Longitudinal Strain'),
    'p24118_i2': ('Ao Asc Min', 'Ascending Aorta Min Area'),
    'p24119_i2': ('Ao Asc Max', 'Ascending Aorta Max Area'),
}

# ── A1: Domain-level CMR ──
print("\n  --- A1: CMR by Domain ---")

domain_cmr = has_cmr.groupby('domain').agg(
    n=('p24103_i2', 'count'),
    mean_lvef=('p24103_i2', 'mean'),
    mean_lvedv=('p24100_i2', 'mean'),
    mean_lvm=('p24105_i2', 'mean'),
    mean_gls=('p24157_i2', 'mean'),
    mean_ao=('p24119_i2', 'mean'),
    ascvd_rate=('sss', lambda x: np.nan),  # placeholder
).reset_index()

# Add domain penetrance
domain_pen = atlas.groupby('domain')['domain_penetrance_by_60_pct'].first().to_dict()
domain_cmr['penetrance'] = domain_cmr['domain'].map(domain_pen)
domain_cmr = domain_cmr[domain_cmr['n'] >= 20].sort_values('mean_lvm', ascending=False)

print(f"\n  {'Domain':>25s}  {'N':>5s}  {'LVEF%':>6s}  {'LV Mass':>8s}  {'LVEDV':>6s}  {'GLS%':>6s}  {'Pen%':>5s}")
for _, r in domain_cmr.iterrows():
    pen_str = f"{r['penetrance']:.1f}" if not np.isnan(r.get('penetrance', np.nan)) else 'N/A'
    print(f"  {r['domain']:>25s}  {int(r['n']):>5d}  {r['mean_lvef']:>6.1f}  {r['mean_lvm']:>8.1f}  "
          f"{r['mean_lvedv']:>6.1f}  {r['mean_gls']:>6.1f}  {pen_str:>5s}")

# ── A2: CMR by mechanism type ──
print("\n  --- A2: CMR by Mechanism Type ---")

mech_cmr = has_cmr.groupby('mech_type').agg(
    n=('p24103_i2', 'count'),
    mean_lvef=('p24103_i2', 'mean'),
    std_lvef=('p24103_i2', 'std'),
    mean_lvedv=('p24100_i2', 'mean'),
    mean_lvesv=('p24101_i2', 'mean'),
    mean_lvm=('p24105_i2', 'mean'),
    mean_gls=('p24157_i2', 'mean'),
    mean_ao_max=('p24119_i2', 'mean'),
).reset_index()

print(f"\n  {'Mechanism':>30s}  {'N':>5s}  {'LVEF':>6s}  {'LVEDV':>6s}  {'LV Mass':>8s}  {'GLS':>6s}")
for _, r in mech_cmr.iterrows():
    print(f"  {r['mech_type']:>30s}  {int(r['n']):>5d}  {r['mean_lvef']:>6.1f}  {r['mean_lvedv']:>6.1f}  "
          f"{r['mean_lvm']:>8.1f}  {r['mean_gls']:>6.1f}")

# Statistical tests
for measure, desc in [('p24103_i2', 'LVEF'), ('p24105_i2', 'LV Mass'),
                       ('p24100_i2', 'LVEDV'), ('p24157_i2', 'GLS')]:
    groups = []
    labels = []
    for mt in ['Type 1: Structural', 'Type 2: Functional',
               'Type 1: Structural (domain)', 'Type 2: Functional (domain)']:
        vals = has_cmr[has_cmr['mech_type'] == mt][measure].dropna()
        if len(vals) >= 10:
            groups.append(vals)
            labels.append(mt)

    if len(groups) >= 2:
        # Combine Type 1s and Type 2s
        t1_all = pd.concat([has_cmr[has_cmr['mech_type'].str.contains('Type 1')][measure].dropna()])
        t2_all = pd.concat([has_cmr[has_cmr['mech_type'].str.contains('Type 2')][measure].dropna()])
        if len(t1_all) >= 10 and len(t2_all) >= 10:
            u, p = stats.mannwhitneyu(t1_all, t2_all)
            print(f"\n  {desc}: Type 1 (n={len(t1_all)}, median={t1_all.median():.1f}) vs "
                  f"Type 2 (n={len(t2_all)}, median={t2_all.median():.1f}), P={p:.4f}")

# ── A3: Correlation: domain penetrance vs CMR measures ──
print("\n  --- A3: Domain Penetrance vs Cardiac Structure ---")

if len(domain_cmr) >= 5:
    for measure, desc in [('mean_lvm', 'LV Mass'), ('mean_lvef', 'LVEF'),
                           ('mean_gls', 'GLS'), ('mean_lvedv', 'LVEDV')]:
        valid = domain_cmr.dropna(subset=[measure, 'penetrance'])
        if len(valid) >= 5:
            r, p = stats.spearmanr(valid['penetrance'], valid[measure])
            print(f"  Penetrance vs {desc}: rho={r:.3f}, P={p:.4f} (n={len(valid)} domains)")

# ── A4: SSS correlation with CMR ──
print("\n  --- A4: SSS vs Cardiac Structure (individual level) ---")

for measure, desc in [('p24103_i2', 'LVEF'), ('p24105_i2', 'LV Mass'),
                       ('p24100_i2', 'LVEDV'), ('p24157_i2', 'GLS')]:
    valid = has_cmr.dropna(subset=['sss', measure])
    if len(valid) >= 50:
        r, p = stats.spearmanr(valid['sss'], valid[measure])
        print(f"  SSS vs {desc}: rho={r:.3f}, P={p:.4f} (n={len(valid)})")

# Save domain CMR
domain_cmr.to_csv(os.path.join(ANALYSIS, "cardiac_mri_by_domain.csv"), index=False)

# ============================================================================
# ANALYSIS B: BIOMARKER PROTOCOL BY DOMAIN
# ============================================================================
print("\n\n" + "=" * 80)
print("ANALYSIS B: BIOMARKER LEVELS BY DOMAIN")
print("Which biochemical tests differ by mutation mechanism?")
print("=" * 80)

# Merge carriers with biomarker data
bio = carriers_unique.copy()
bio['eid'] = bio['eid'].astype(str)

# ApoB and Lp(a)
if len(apob) > 0:
    apob['eid'] = apob['eid'].astype(str)
    bio = bio.merge(apob, on='eid', how='left')

# CRP
if len(crp) > 0:
    crp['eid'] = crp['eid'].astype(str)
    bio = bio.merge(crp, on='eid', how='left')

# Lipids
if len(lipids) > 0:
    lipids['eid'] = lipids['eid'].astype(str)
    bio = bio.merge(lipids, on='eid', how='left')

bio['domain_ddg'] = bio['domain'].map(domain_ddg)
bio['mech_type'] = 'Unknown'
bio.loc[bio['ddG'].notna() & (bio['ddG'] >= 2.0), 'mech_type'] = 'Type 1'
bio.loc[bio['ddG'].notna() & (bio['ddG'] < 2.0), 'mech_type'] = 'Type 2'
bio.loc[(bio['mech_type'] == 'Unknown') & bio['domain_ddg'].notna() &
        (bio['domain_ddg'] >= 5.0), 'mech_type'] = 'Type 1'
bio.loc[(bio['mech_type'] == 'Unknown') & bio['domain_ddg'].notna() &
        (bio['domain_ddg'] < 5.0), 'mech_type'] = 'Type 2'

# Biomarker columns to check
bio_cols = {}
for col in bio.columns:
    col_lower = col.lower()
    if 'p30890' in col: bio_cols[col] = 'ApoB (g/L)'
    elif 'p30900' in col and '_i0' in col: bio_cols[col] = 'Lp(a) (nmol/L)'
    elif 'p30710' in col and '_i0' in col: bio_cols[col] = 'CRP (mg/L)'
    elif 'p30780' in col and '_i0' in col: bio_cols[col] = 'LDL-C (mmol/L)'
    elif 'p30760' in col and '_i0' in col: bio_cols[col] = 'HDL-C (mmol/L)'
    elif 'p30870' in col and '_i0' in col: bio_cols[col] = 'Triglycerides (mmol/L)'

print(f"\n  Biomarker data available: {len(bio_cols)} measures")
for col, desc in bio_cols.items():
    n = bio[col].notna().sum()
    print(f"    {desc:>25s} ({col}): {n:>6d} carriers with data")

# ── B1: Biomarkers by mechanism type ──
print("\n  --- B1: Biomarkers by Mechanism ---")

for col, desc in bio_cols.items():
    t1 = bio[(bio['mech_type'] == 'Type 1') & bio[col].notna()][col]
    t2 = bio[(bio['mech_type'] == 'Type 2') & bio[col].notna()][col]

    if len(t1) >= 20 and len(t2) >= 20:
        u, p = stats.mannwhitneyu(t1, t2)
        print(f"\n  {desc}:")
        print(f"    Type 1 (n={len(t1)}): median={t1.median():.2f}, mean={t1.mean():.2f}")
        print(f"    Type 2 (n={len(t2)}): median={t2.median():.2f}, mean={t2.mean():.2f}")
        print(f"    P={p:.4f}")

# ── B2: Biomarkers by domain ──
print("\n  --- B2: Key Biomarkers by Domain ---")

for col, desc in list(bio_cols.items())[:3]:  # Top 3 biomarkers
    domain_bio = bio.groupby('domain').agg(
        n=(col, lambda x: x.notna().sum()),
        median=(col, 'median'),
        mean=(col, 'mean'),
    ).reset_index()
    domain_bio = domain_bio[domain_bio['n'] >= 20].sort_values('median', ascending=False)

    if len(domain_bio) >= 3:
        print(f"\n  {desc} by domain:")
        print(f"  {'Domain':>25s}  {'N':>5s}  {'Median':>8s}  {'Mean':>8s}")
        for _, r in domain_bio.iterrows():
            print(f"  {r['domain']:>25s}  {int(r['n']):>5d}  {r['median']:>8.2f}  {r['mean']:>8.2f}")

# ============================================================================
# CLINICAL PROTOCOL TABLE
# ============================================================================
print("\n\n" + "=" * 80)
print("STRUCTURE-GUIDED CLINICAL MONITORING PROTOCOL")
print("=" * 80)

protocol = []

# Type 1: Structural destabilisation
protocol.append({
    'mechanism': 'Type 1: Structural Destabilisation',
    'brown_goldstein': 'Class 2 (transport-defective)',
    'molecular': 'Protein misfolding -> ER retention -> no surface receptor',
    'baseline_labs': 'LDL-C, ApoB, Lp(a), hsCRP, LFTs, CK',
    'monitoring_labs': 'LDL-C + ApoB q3mo (first year), then q6mo',
    'additional_tests': 'ApoB (better than LDL for particle count), Lp(a) (independent risk), PCSK9 free/total',
    'imaging': 'CT calcium score at diagnosis, repeat q5yr. Cardiac MRI if LV mass elevated.',
    'cardiac_mri_rationale': 'LV mass may be elevated from chronic pressure overload + atherosclerosis',
    'first_line': 'PCSK9i + ezetimibe (bypass LDLR pathway)',
    'escalation': 'Evinacumab (ANGPTL3i) if LDL >3.5 despite max therapy',
    'avoid': 'Statin monotherapy trial (futile if no functional receptor)',
    'follow_up_interval': 'q3mo until target, then q6mo',
    'target': 'LDL-C <1.4 mmol/L (ESC very high risk) or >50% reduction',
})

protocol.append({
    'mechanism': 'Type 2: Functional Disruption',
    'brown_goldstein': 'Class 3/4/5 (binding/internalisation/recycling)',
    'molecular': 'Protein folds and reaches surface but function impaired',
    'baseline_labs': 'LDL-C, ApoB, Lp(a), hsCRP, LFTs, CK, HbA1c',
    'monitoring_labs': 'LDL-C q3mo, ApoB q6mo',
    'additional_tests': 'sdLDL or LDL particle number (NMR), oxLDL if available',
    'imaging': 'CT calcium score at diagnosis. Carotid IMT optional.',
    'cardiac_mri_rationale': 'Monitor for subclinical cardiomyopathy if prolonged LDL exposure',
    'first_line': 'High-intensity statin + ezetimibe',
    'escalation': 'Add PCSK9i if LDL >1.8 after 3 months',
    'avoid': 'Delay in treatment escalation',
    'follow_up_interval': 'q3mo first year, then q6-12mo if at target',
    'target': 'LDL-C <1.4 mmol/L or >50% reduction',
})

# Domain-specific additions
domain_protocols = {
    'EGF-like A': {
        'special': 'PCSK9 binding site mutation - MONITOR PCSK9i RESPONSE CLOSELY',
        'additional': 'Measure free PCSK9 levels pre/post PCSK9i. Consider switching PCSK9i class if poor response.',
        'rationale': 'EGF-A is where PCSK9 binds LDLR. Mutations here may alter drug-target interaction.',
    },
    'Ligand-binding R4': {
        'special': 'HIGHEST RISK DOMAIN (39.6% penetrance)',
        'additional': 'Urgent CV risk assessment. Family cascade screening. CT calcium immediately.',
        'rationale': 'R4-R5 are the essential LDL binding repeats. Near-null phenotype expected.',
    },
    'Ligand-binding R5': {
        'special': 'PRIMARY LDL BINDING SITE',
        'additional': 'ApoB particle count essential. Consider LDL apheresis if refractory.',
        'rationale': 'Critical binding repeat. Complete loss of LDL binding likely.',
    },
    'Beta-propeller': {
        'special': 'RECYCLING DOMAIN - pH-dependent LDL release',
        'additional': 'Monitor LDL response to statin carefully - recycling defect may retain partial function.',
        'rationale': 'Receptor cannot release LDL at endosomal pH. PCSK9i blocks degradation of impaired receptor.',
    },
    'Ligand-binding R2': {
        'special': 'LOWEST RISK DOMAIN (8.7% penetrance)',
        'additional': 'Standard follow-up may suffice. Genetic counselling reassurance.',
        'rationale': 'R2 is auxiliary binding repeat. Partial function often retained.',
    },
}

print("\n  GENERAL PROTOCOL BY MECHANISM:")
for p in protocol:
    print(f"\n  {p['mechanism']}:")
    print(f"    B&G Class:        {p['brown_goldstein']}")
    print(f"    Molecular:        {p['molecular']}")
    print(f"    Baseline labs:    {p['baseline_labs']}")
    print(f"    Monitoring:       {p['monitoring_labs']}")
    print(f"    Additional:       {p['additional_tests']}")
    print(f"    Imaging:          {p['imaging']}")
    print(f"    First-line Rx:    {p['first_line']}")
    print(f"    Escalation:       {p['escalation']}")
    print(f"    Avoid:            {p['avoid']}")
    print(f"    Target:           {p['target']}")

print("\n\n  DOMAIN-SPECIFIC ADDITIONS:")
for dom, info in domain_protocols.items():
    print(f"\n  {dom}: {info['special']}")
    print(f"    Additional:  {info['additional']}")
    print(f"    Rationale:   {info['rationale']}")

# Save protocol
pd.DataFrame(protocol).to_csv(os.path.join(ANALYSIS, "biomarker_protocol_by_domain.csv"), index=False)

# ============================================================================
# FIGURE
# ============================================================================
print("\n" + "=" * 80)
print("GENERATING FIGURES")
print("=" * 80)

fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# Panel A: LV Mass by domain
ax = axes[0, 0]
if len(domain_cmr) > 3:
    dcmr = domain_cmr.sort_values('mean_lvm')
    colours_dom = []
    for _, r in dcmr.iterrows():
        pen = r.get('penetrance', 0)
        if pen and pen > 25: colours_dom.append('#b2182b')
        elif pen and pen > 15: colours_dom.append('#ef8a62')
        else: colours_dom.append('#2166ac')
    bars = ax.barh(range(len(dcmr)), dcmr['mean_lvm'], color=colours_dom, edgecolor='black')
    ax.set_yticks(range(len(dcmr)))
    ax.set_yticklabels([f"{r['domain']}\n(n={int(r['n'])})" for _, r in dcmr.iterrows()],
                       fontsize=8, fontfamily='Arial')
    ax.set_xlabel('Mean LV Mass (g)', fontsize=11, fontfamily='Arial')
    ax.set_title('A. LV Mass by LDLR Domain', fontsize=12, fontweight='bold', fontfamily='Arial')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Panel B: LVEF by mechanism
ax = axes[0, 1]
mech_groups = []
mech_labels_p = []
mech_colours = []
for mt, col in [('Type 1: Structural', '#b2182b'), ('Type 1: Structural (domain)', '#ef8a62'),
                ('Type 2: Functional', '#2166ac'), ('Type 2: Functional (domain)', '#67a9cf')]:
    vals = has_cmr[has_cmr['mech_type'] == mt]['p24103_i2'].dropna()
    if len(vals) >= 20:
        mech_groups.append(vals.values)
        mech_labels_p.append(f"{mt.split(':')[1].strip()}\n(n={len(vals)})")
        mech_colours.append(col)

if mech_groups:
    bp = ax.boxplot(mech_groups, labels=mech_labels_p, patch_artist=True,
                    medianprops=dict(color='black', linewidth=2), showfliers=False)
    for patch, col in zip(bp['boxes'], mech_colours):
        patch.set_facecolor(col)
        patch.set_alpha(0.7)
    ax.set_ylabel('LV Ejection Fraction (%)', fontsize=11, fontfamily='Arial')
    ax.set_title('B. LVEF by Mutation Mechanism', fontsize=12, fontweight='bold', fontfamily='Arial')
    ax.tick_params(axis='x', labelsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: Domain penetrance vs LV mass
ax = axes[1, 0]
valid_corr = domain_cmr.dropna(subset=['mean_lvm', 'penetrance'])
if len(valid_corr) >= 4:
    ax.scatter(valid_corr['penetrance'], valid_corr['mean_lvm'], s=valid_corr['n']/2,
              c='#b2182b', alpha=0.7, edgecolors='black')
    for _, r in valid_corr.iterrows():
        ax.annotate(r['domain'].replace('Ligand-binding ', 'LB-').replace('EGF-like ', 'EGF-'),
                   (r['penetrance'], r['mean_lvm']), fontsize=6, ha='center', va='bottom')
    r_val, p_val = stats.spearmanr(valid_corr['penetrance'], valid_corr['mean_lvm'])
    ax.set_xlabel('ASCVD Penetrance by Age 60 (%)', fontsize=11, fontfamily='Arial')
    ax.set_ylabel('Mean LV Mass (g)', fontsize=11, fontfamily='Arial')
    ax.set_title(f'C. Penetrance vs LV Mass\n(rho={r_val:.3f}, P={p_val:.3f})',
                 fontsize=12, fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Clinical protocol summary
ax = axes[1, 1]
ax.axis('off')
table_data = [
    ['Measure', 'Type 1\n(Structural)', 'Type 2\n(Functional)'],
    ['Mechanism', 'Misfolding\nER retention', 'Impaired binding\nor recycling'],
    ['First-line Rx', 'PCSK9i +\nezetimibe', 'Statin +\nezetimibe'],
    ['Key biomarker', 'ApoB\n(particle count)', 'LDL-C\n(standard)'],
    ['Imaging', 'CT calcium +\ncardiac MRI', 'CT calcium'],
    ['Monitoring', 'q3 months', 'q3-6 months'],
    ['Escalation', 'Evinacumab\n/ apheresis', 'Add PCSK9i\nif not at target'],
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
ax.set_title('D. Structure-Guided Clinical Protocol', fontsize=12,
             fontweight='bold', fontfamily='Arial')

plt.tight_layout()
fig_path = os.path.join(FIGURES, "Figure_Cardiac_MRI_Domains.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {fig_path}")

print("\n" + "=" * 80)
print("COMPLETE")
print("=" * 80)

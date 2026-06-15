#!/usr/bin/env python3
"""
53_aortic_phenotype_panel.py
=============================
Full aortic phenotype panel by LDLR domain and mutation mechanism.

Integrates:
  1. Aortic distensibility (CMR p24118-p24123)
  2. Aortic valve sclerosis (ICD-10 I35.x from HES)
  3. Coronary artery calcium score (CT Category 162)
  4. Lipoprotein(a) — independent driver of aortic valve calcification

Stratified by LDLR domain and Type 1/Type 2/Benign mechanism.
Lp(a) interaction analysis: does Lp(a) amplify aortic disease risk
differentially by domain or mechanism?

Outputs:
  alphafold/analysis/aortic_phenotype_by_domain.csv
  alphafold/analysis/aortic_phenotype_by_mechanism.csv
  alphafold/analysis/aortic_sclerosis_prevalence.csv
  alphafold/analysis/figures/Figure_Aortic_Phenotype_Panel.png

Author: Dr Nader Genedy
Date:   April 2026
"""

import pandas as pd
import numpy as np
import os
import json
import warnings
warnings.filterwarnings('ignore')
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")
FIGURES = os.path.join(ANALYSIS, "figures")
BACKUP = r"D:/calon_backup_data"
os.makedirs(FIGURES, exist_ok=True)

print("=" * 80)
print("AORTIC PHENOTYPE PANEL BY LDLR DOMAIN")
print("=" * 80)

# ============================================================================
# STEP 1: LOAD ALL DATA
# ============================================================================
print("\n[1/6] Loading data...")

# 1a. Annotated carriers (from script 52)
carriers = pd.read_csv(os.path.join(ANALYSIS, "ukb_carriers_annotated.csv"))
carriers['eid'] = carriers['eid'].astype(str)
# Keep one row per patient (highest SSS)
patients = carriers.sort_values('sss', ascending=False).drop_duplicates('eid')
print(f"  Annotated carriers: {len(patients)} unique patients")

# 1b. Domain atlas for penetrance
atlas = pd.read_csv(os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas_v3.csv"))
pen_map = atlas.groupby('domain')['domain_penetrance_by_60_pct'].first().to_dict()

# 1c. LV function data
lv = pd.read_csv(os.path.join(BACKUP, "calon_extra_mri_lv (1).csv"))
lv.columns = [c.replace('participant.', '') for c in lv.columns]
lv['eid'] = lv['eid'].astype(str)
print(f"  LV data: {lv['p24100_i2'].notna().sum()} with LVEDV")

# 1d. Aortic CMR data (all 6 fields)
aorta = pd.read_csv(os.path.join(BACKUP, "calon_extra_mri_aorta (1).csv"))
aorta.columns = [c.replace('participant.', '') for c in aorta.columns]
aorta['eid'] = aorta['eid'].astype(str)
print(f"  Aortic CMR: {aorta['p24120_i2'].notna().sum()} with ascending distensibility")

# 1e. Strain data
strain = pd.read_csv(os.path.join(BACKUP, "calon_extra_mri_strain.csv"))
strain.columns = [c.replace('participant.', '') for c in strain.columns]
strain['eid'] = strain['eid'].astype(str)
print(f"  Strain: {strain['p24157_i2'].notna().sum()} with circumferential strain")

# 1f. Category 162 (CACS)
cat162 = pd.read_csv(os.path.join(BACKUP, "calon_extra_mri_cat162.csv"))
cat162.columns = [c.replace('participant.', '') for c in cat162.columns]
cat162['eid'] = cat162['eid'].astype(str)
print(f"  Cat 162 (CACS): {cat162['p31063_i2'].notna().sum()} with LV mass (alt)")

# 1g. Lp(a) and ApoB biomarkers — CORRECT UKB FIELDS
# FIELD MAPPING (verified 2026-04-05 on UKB-RAP):
#   p30790 = Lipoprotein A (Lp(a), nmol/L) — immunoturbidimetric
#   p30640 = Apolipoprotein B (g/L)
#   p30900 = Number of proteins measured (proteomics) — NOT Lp(a)!
#   p30890 = NOT listed — was wrong field in previous extraction
lpa_raw = pd.read_csv("D:/calon_lpa_CORRECT.csv")
lpa_raw.columns = [c.replace('participant.', '') for c in lpa_raw.columns]
lpa_raw['eid'] = lpa_raw['eid'].astype(str)
lpa_raw['lpa'] = pd.to_numeric(lpa_raw['p30790_i0'], errors='coerce')
print(f"  Lp(a) (UKB p30790): {lpa_raw['lpa'].notna().sum()} with data "
      f"(median={lpa_raw['lpa'].median():.1f} nmol/L)")

apob_raw = pd.read_csv("D:/calon_apob_CORRECT.csv")
apob_raw.columns = [c.replace('participant.', '') for c in apob_raw.columns]
apob_raw['eid'] = apob_raw['eid'].astype(str)
apob_raw['apob'] = pd.to_numeric(apob_raw['p30640_i0'], errors='coerce')
print(f"  ApoB (UKB p30640): {apob_raw['apob'].notna().sum()} with data "
      f"(median={apob_raw['apob'].median():.2f} g/L)")

# Combine
lpa_apob = lpa_raw[['eid', 'lpa']].merge(apob_raw[['eid', 'apob']], on='eid', how='outer')

# 1h. ICD-10 codes — parse JSON array in p41270
print("  Loading ICD-10 codes (this may take a moment)...")
icd = pd.read_csv(os.path.join(BACKUP, "calon_batch_icd10.csv"),
                   usecols=['participant.eid', 'participant.p41270'],
                   dtype=str)
icd.columns = ['eid', 'icd10_raw']


def parse_icd10(raw):
    """Parse the JSON-like ICD-10 array from UKB p41270."""
    if pd.isna(raw) or raw == '':
        return []
    try:
        codes = json.loads(raw)
        if isinstance(codes, list):
            return [str(c).strip().strip('"') for c in codes]
    except (json.JSONDecodeError, TypeError):
        pass
    # Fallback: try splitting on common delimiters
    raw = str(raw).strip('[]"')
    return [c.strip().strip('"') for c in raw.split(',') if c.strip()]


icd['icd10_codes'] = icd['icd10_raw'].apply(parse_icd10)

# Create binary aortic valve phenotypes
icd['ao_sclerosis'] = icd['icd10_codes'].apply(
    lambda codes: int(any(c.startswith('I35') for c in codes)))
icd['ao_stenosis'] = icd['icd10_codes'].apply(
    lambda codes: int(any(c in ('I350', 'I352') for c in codes)))
icd['ao_regurg'] = icd['icd10_codes'].apply(
    lambda codes: int(any(c == 'I351' for c in codes)))

n_scl = icd['ao_sclerosis'].sum()
n_sten = icd['ao_stenosis'].sum()
n_reg = icd['ao_regurg'].sum()
print(f"  ICD-10 aortic phenotypes in full UKB:")
print(f"    Any I35 (sclerosis/stenosis): {n_scl} ({100*n_scl/len(icd):.2f}%)")
print(f"    I350/I352 (stenosis):         {n_sten}")
print(f"    I351 (regurgitation):          {n_reg}")

icd_pheno = icd[['eid', 'ao_sclerosis', 'ao_stenosis', 'ao_regurg']].copy()

# ── Merge everything onto carriers ──
print("\n  Merging all data onto carriers...")
df = patients.copy()
for data, name in [(lv, 'LV'), (aorta, 'Aorta'), (strain, 'Strain'),
                    (cat162, 'Cat162'), (lpa_apob[['eid', 'lpa', 'apob']], 'Lpa'),
                    (icd_pheno, 'ICD')]:
    df = df.merge(data, on='eid', how='left', suffixes=('', f'_dup_{name}'))
    # Remove duplicate columns
    df = df[[c for c in df.columns if '_dup_' not in c]]

print(f"  Merged dataset: {len(df)} patients")
print(f"    With aortic distensibility: {df['p24120_i2'].notna().sum()}")
print(f"    With aortic sclerosis code: {df['ao_sclerosis'].sum()}")
print(f"    With Lp(a):                 {df['lpa'].notna().sum()}")
print(f"    With ApoB:                  {df['apob'].notna().sum()}")
print(f"    With LV mass:               {df['p24105_i2'].notna().sum()}")

# ============================================================================
# STEP 2: CREATE DERIVED AORTIC PHENOTYPES
# ============================================================================
print("\n[2/6] Creating aortic phenotype variables...")

# Rename for clarity
df['ao_asc_min'] = df['p24118_i2']       # Ascending aorta min area (mm2)
df['ao_asc_max'] = df['p24119_i2']       # Ascending aorta max area (mm2)
df['ao_asc_dist'] = df['p24120_i2']      # Ascending aorta distensibility
df['ao_desc_min'] = df['p24121_i2']      # Descending aorta min area
df['ao_desc_max'] = df['p24122_i2']      # Descending aorta max area
df['ao_desc_dist'] = df['p24123_i2']     # Descending aorta distensibility
df['lvef'] = df['p24103_i2']
df['lvm'] = df['p24105_i2']
df['lvedv'] = df['p24100_i2']
df['gls'] = df['p24157_i2']

# Derived: aortic pulsatility index (fractional area change)
# = (max_area - min_area) / max_area
df['ao_asc_pulsatility'] = np.where(
    df['ao_asc_max'].notna() & (df['ao_asc_max'] > 0),
    (df['ao_asc_max'] - df['ao_asc_min']) / df['ao_asc_max'],
    np.nan
)
df['ao_desc_pulsatility'] = np.where(
    df['ao_desc_max'].notna() & (df['ao_desc_max'] > 0),
    (df['ao_desc_max'] - df['ao_desc_min']) / df['ao_desc_max'],
    np.nan
)

# CACS from Cat 162
df['cacs'] = pd.to_numeric(df.get('p31063_i2'), errors='coerce')

# Lp(a) thresholds (nmol/L): >125 = high risk (ESC 2019), >75 = elevated
df['lpa_high'] = (df['lpa'] > 125).astype(int).where(df['lpa'].notna())
df['lpa_elevated'] = (df['lpa'] > 75).astype(int).where(df['lpa'].notna())

for var, name in [('ao_asc_dist', 'Asc distensibility'),
                  ('ao_desc_dist', 'Desc distensibility'),
                  ('ao_asc_pulsatility', 'Asc pulsatility'),
                  ('ao_sclerosis', 'Aortic sclerosis (ICD)'),
                  ('cacs', 'CACS')]:
    if var in df.columns:
        n = df[var].notna().sum() if var != 'ao_sclerosis' else df[var].sum()
        print(f"  {name:>30s}: {n}")

# ============================================================================
# STEP 3: ANALYSIS BY DOMAIN
# ============================================================================
print("\n[3/6] Aortic phenotypes by LDLR domain...")

has_ao = df[df['ao_asc_dist'].notna()].copy()
print(f"  Carriers with aortic distensibility: {len(has_ao)}")

# Domain-level aggregation
domain_ao = has_ao.groupby('ldlr_domain').agg(
    n=('ao_asc_dist', 'count'),
    mean_asc_dist=('ao_asc_dist', 'mean'),
    std_asc_dist=('ao_asc_dist', 'std'),
    median_asc_dist=('ao_asc_dist', 'median'),
    mean_desc_dist=('ao_desc_dist', 'mean'),
    mean_asc_max=('ao_asc_max', 'mean'),
    mean_asc_puls=('ao_asc_pulsatility', 'mean'),
    mean_lvef=('lvef', 'mean'),
    mean_lvm=('lvm', 'mean'),
    mean_gls=('gls', 'mean'),
).reset_index()

# Add sclerosis prevalence and Lp(a) by domain (use ALL carriers)
scl_by_dom = df.groupby('ldlr_domain').agg(
    n_total=('eid', 'count'),
    n_sclerosis=('ao_sclerosis', 'sum'),
    n_lpa=('lpa', lambda x: x.notna().sum()),
    median_lpa=('lpa', 'median'),
    mean_lpa=('lpa', 'mean'),
    n_lpa_high=('lpa_high', 'sum'),
    median_apob=('apob', 'median'),
    mean_apob=('apob', 'mean'),
).reset_index()
scl_by_dom['sclerosis_pct'] = 100 * scl_by_dom['n_sclerosis'] / scl_by_dom['n_total']
scl_by_dom['lpa_high_pct'] = np.where(
    scl_by_dom['n_lpa'] > 0,
    100 * scl_by_dom['n_lpa_high'] / scl_by_dom['n_lpa'], np.nan)

domain_ao = domain_ao.merge(
    scl_by_dom[['ldlr_domain', 'n_total', 'n_sclerosis', 'sclerosis_pct',
                'n_lpa', 'median_lpa', 'mean_lpa', 'lpa_high_pct', 'median_apob']],
    on='ldlr_domain', how='left')

# Add CACS mean by domain
cacs_dom = df[df['cacs'].notna()].groupby('ldlr_domain')['cacs'].agg(
    ['count', 'mean', 'median']).reset_index()
cacs_dom.columns = ['ldlr_domain', 'n_cacs', 'mean_cacs', 'median_cacs']
domain_ao = domain_ao.merge(cacs_dom, on='ldlr_domain', how='left')

# Add penetrance
domain_ao['penetrance'] = domain_ao['ldlr_domain'].map(pen_map)

# Filter to domains with enough data
domain_ao = domain_ao[domain_ao['n'] >= 5].sort_values('mean_asc_dist')

print(f"\n  {'Domain':>25s}  {'N':>4s}  {'Asc Dist':>9s}  {'Scl%':>5s}  {'Lp(a)':>7s}  {'Lpa>125%':>8s}  {'ApoB':>6s}  {'Pen%':>5s}")
print("  " + "-" * 85)
for _, r in domain_ao.iterrows():
    pen = f"{r['penetrance']:.1f}" if pd.notna(r['penetrance']) else 'N/A'
    scl = f"{r['sclerosis_pct']:.1f}" if pd.notna(r['sclerosis_pct']) else 'N/A'
    lpa_m = f"{r['median_lpa']:.0f}" if pd.notna(r['median_lpa']) else 'N/A'
    lpa_h = f"{r['lpa_high_pct']:.1f}" if pd.notna(r['lpa_high_pct']) else 'N/A'
    apob_m = f"{r['median_apob']:.2f}" if pd.notna(r['median_apob']) else 'N/A'
    print(f"  {r['ldlr_domain']:>25s}  {int(r['n']):>4d}  {r['mean_asc_dist']:>9.2f}  "
          f"{scl:>5s}  {lpa_m:>7s}  {lpa_h:>8s}  {apob_m:>6s}  {pen:>5s}")

domain_ao.to_csv(os.path.join(ANALYSIS, "aortic_phenotype_by_domain.csv"), index=False)

# ============================================================================
# STEP 4: ANALYSIS BY MECHANISM TYPE
# ============================================================================
print("\n[4/6] Aortic phenotypes by mechanism type...")

mech_groups = ['Type 1: Structural', 'Type 2: Functional', 'Benign']
mech_results = []

for mt in mech_groups:
    sub = df[df['mech_simple'] == mt]
    row = {
        'mechanism': mt,
        'n_total': len(sub),
        'n_ao_dist': sub['ao_asc_dist'].notna().sum(),
        'mean_asc_dist': sub['ao_asc_dist'].mean(),
        'std_asc_dist': sub['ao_asc_dist'].std(),
        'median_asc_dist': sub['ao_asc_dist'].median(),
        'mean_desc_dist': sub['ao_desc_dist'].mean(),
        'mean_asc_puls': sub['ao_asc_pulsatility'].mean(),
        'n_sclerosis': sub['ao_sclerosis'].sum(),
        'sclerosis_pct': 100 * sub['ao_sclerosis'].sum() / len(sub) if len(sub) > 0 else 0,
        'n_lpa': sub['lpa'].notna().sum(),
        'median_lpa': sub['lpa'].median(),
        'mean_lpa': sub['lpa'].mean(),
        'lpa_high_pct': 100 * sub['lpa_high'].sum() / sub['lpa'].notna().sum() if sub['lpa'].notna().sum() > 0 else 0,
        'median_apob': sub['apob'].median(),
        'mean_lvm': sub['lvm'].mean(),
        'mean_lvef': sub['lvef'].mean(),
        'mean_cacs': sub['cacs'].mean(),
    }
    mech_results.append(row)

mech_df = pd.DataFrame(mech_results)
print(f"\n  {'Mechanism':>25s}  {'N':>5s}  {'Asc Dist':>9s}  {'Scl%':>5s}  {'Lp(a)':>7s}  {'Lpa>125%':>8s}  {'ApoB':>6s}  {'LVM':>6s}")
print("  " + "-" * 85)
for _, r in mech_df.iterrows():
    lpa_m = f"{r['median_lpa']:.0f}" if pd.notna(r['median_lpa']) else 'N/A'
    apob_m = f"{r['median_apob']:.2f}" if pd.notna(r['median_apob']) else 'N/A'
    print(f"  {r['mechanism']:>25s}  {int(r['n_total']):>5d}  {r['mean_asc_dist']:>9.2f}  "
          f"{r['sclerosis_pct']:>5.1f}  {lpa_m:>7s}  {r['lpa_high_pct']:>8.1f}  "
          f"{apob_m:>6s}  {r['mean_lvm']:>6.1f}")

# Statistical tests
print(f"\n  --- Statistical Tests ---")
t1 = df[df['mech_simple'] == 'Type 1: Structural']
t2 = df[df['mech_simple'] == 'Type 2: Functional']
ben = df[df['mech_simple'] == 'Benign']

continuous_vars = [
    ('ao_asc_dist', 'Asc Distensibility'),
    ('ao_desc_dist', 'Desc Distensibility'),
    ('ao_asc_pulsatility', 'Asc Pulsatility'),
    ('lpa', 'Lp(a)'),
    ('apob', 'ApoB'),
    ('lvm', 'LV Mass'),
    ('ao_asc_max', 'Asc Aorta Max Area'),
]

for col, name in continuous_vars:
    for label, g1, g2 in [('T1 vs T2', t1, t2), ('T1 vs Benign', t1, ben), ('T2 vs Benign', t2, ben)]:
        v1 = g1[col].dropna()
        v2 = g2[col].dropna()
        if len(v1) >= 10 and len(v2) >= 10:
            u, p = stats.mannwhitneyu(v1, v2, alternative='two-sided')
            pooled_std = np.sqrt((v1.std()**2 + v2.std()**2) / 2)
            d_cohen = (v1.mean() - v2.mean()) / pooled_std if pooled_std > 0 else 0
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            print(f"  {name:>20s} {label}: {v1.median():.2f} vs {v2.median():.2f}, "
                  f"P={p:.4f} {sig}, d={d_cohen:.3f}")

# Aortic sclerosis: Fisher's exact test
print(f"\n  Aortic sclerosis (Fisher's exact):")
for label, g1, g2 in [('T1 vs T2', t1, t2), ('T1 vs Benign', t1, ben), ('T2 vs Benign', t2, ben)]:
    a = g1['ao_sclerosis'].sum()
    b = len(g1) - a
    c = g2['ao_sclerosis'].sum()
    d_val = len(g2) - c
    if a + c > 0:
        table = np.array([[a, b], [c, d_val]])
        odds, p_fish = stats.fisher_exact(table)
        sig = '***' if p_fish < 0.001 else '**' if p_fish < 0.01 else '*' if p_fish < 0.05 else ''
        r1 = 100 * a / (a + b) if (a + b) > 0 else 0
        r2 = 100 * c / (c + d_val) if (c + d_val) > 0 else 0
        print(f"    {label}: {r1:.1f}% vs {r2:.1f}%, OR={odds:.2f}, P={p_fish:.4f} {sig}")

mech_df.to_csv(os.path.join(ANALYSIS, "aortic_phenotype_by_mechanism.csv"), index=False)

# ============================================================================
# STEP 5: CORRELATION ANALYSIS
# ============================================================================
print("\n[5/6] Correlation analyses...")

# 5a. Domain-level: penetrance vs aortic distensibility
valid = domain_ao.dropna(subset=['mean_asc_dist', 'penetrance'])
if len(valid) >= 5:
    for col, name in [('mean_asc_dist', 'Asc Distensibility'),
                       ('mean_desc_dist', 'Desc Distensibility'),
                       ('sclerosis_pct', 'Sclerosis %'),
                       ('mean_lvm', 'LV Mass')]:
        v = valid.dropna(subset=[col])
        if len(v) >= 5:
            rho, p = stats.spearmanr(v['penetrance'], v[col])
            print(f"  Domain penetrance vs {name}: rho={rho:.3f}, P={p:.4f} (n={len(v)} domains)")

# 5b. Individual-level correlations
print(f"\n  Individual-level correlations:")
has_both = df[df['ao_asc_dist'].notna() & df['sss'].notna()]
if len(has_both) >= 20:
    rho, p = stats.spearmanr(has_both['sss'], has_both['ao_asc_dist'])
    print(f"  SSS vs Asc Distensibility: rho={rho:.3f}, P={p:.4f} (n={len(has_both)})")

# FoldX ddG vs distensibility
has_ddg = df[df['ao_asc_dist'].notna() & df['foldx_ddg'].notna()]
if len(has_ddg) >= 20:
    rho, p = stats.spearmanr(has_ddg['foldx_ddg'], has_ddg['ao_asc_dist'])
    print(f"  FoldX ddG vs Asc Distensibility: rho={rho:.3f}, P={p:.4f} (n={len(has_ddg)})")

# Distensibility vs LV mass (structure-function)
has_lvm = df[df['ao_asc_dist'].notna() & df['lvm'].notna()]
if len(has_lvm) >= 20:
    rho, p = stats.spearmanr(has_lvm['ao_asc_dist'], has_lvm['lvm'])
    print(f"  Asc Distensibility vs LV Mass: rho={rho:.3f}, P={p:.4f} (n={len(has_lvm)})")

# 5c. Lp(a)-specific correlations
print(f"\n  --- Lp(a) Correlations ---")
has_lpa = df[df['lpa'].notna()]
print(f"  Carriers with Lp(a): {len(has_lpa)}")

# Lp(a) vs distensibility
lpa_ao = df[df['lpa'].notna() & df['ao_asc_dist'].notna()]
if len(lpa_ao) >= 20:
    rho, p = stats.spearmanr(lpa_ao['lpa'], lpa_ao['ao_asc_dist'])
    print(f"  Lp(a) vs Asc Distensibility: rho={rho:.3f}, P={p:.4f} (n={len(lpa_ao)})")

# Lp(a) vs aortic sclerosis (point-biserial)
if has_lpa['ao_sclerosis'].sum() >= 5:
    lpa_scl = has_lpa[has_lpa['ao_sclerosis'] == 1]['lpa']
    lpa_noscl = has_lpa[has_lpa['ao_sclerosis'] == 0]['lpa']
    if len(lpa_scl) >= 5 and len(lpa_noscl) >= 10:
        u, p = stats.mannwhitneyu(lpa_scl, lpa_noscl, alternative='two-sided')
        print(f"  Lp(a) in sclerosis vs no sclerosis: {lpa_scl.median():.0f} vs {lpa_noscl.median():.0f} nmol/L, P={p:.4f}")

# Lp(a) vs SSS
lpa_sss = df[df['lpa'].notna() & df['sss'].notna()]
if len(lpa_sss) >= 20:
    rho, p = stats.spearmanr(lpa_sss['sss'], lpa_sss['lpa'])
    print(f"  SSS vs Lp(a): rho={rho:.3f}, P={p:.4f} (n={len(lpa_sss)})")

# Lp(a) vs LV mass
lpa_lvm = df[df['lpa'].notna() & df['lvm'].notna()]
if len(lpa_lvm) >= 20:
    rho, p = stats.spearmanr(lpa_lvm['lpa'], lpa_lvm['lvm'])
    print(f"  Lp(a) vs LV Mass: rho={rho:.3f}, P={p:.4f} (n={len(lpa_lvm)})")

# 5d. Lp(a) interaction: high vs low Lp(a) x mechanism
print(f"\n  --- Lp(a) x Mechanism Interaction ---")
for mt in mech_groups:
    sub = df[(df['mech_simple'] == mt) & df['lpa'].notna()]
    if len(sub) >= 10:
        hi = sub[sub['lpa'] > 125]
        lo = sub[sub['lpa'] <= 125]
        scl_hi = hi['ao_sclerosis'].mean() * 100 if len(hi) > 0 else 0
        scl_lo = lo['ao_sclerosis'].mean() * 100 if len(lo) > 0 else 0
        short = mt.split(':')[1].strip() if ':' in mt else mt
        print(f"  {short}: Lp(a)>125 sclerosis={scl_hi:.1f}% (n={len(hi)}), "
              f"Lp(a)<=125 sclerosis={scl_lo:.1f}% (n={len(lo)})")

# 5e. Domain-level: Lp(a) vs penetrance
print(f"\n  --- Domain-level Lp(a) ---")
dom_lpa = scl_by_dom[scl_by_dom['n_lpa'] >= 5].dropna(subset=['median_lpa'])
dom_lpa['penetrance'] = dom_lpa['ldlr_domain'].map(pen_map)
valid_lpa = dom_lpa.dropna(subset=['penetrance'])
if len(valid_lpa) >= 5:
    rho, p = stats.spearmanr(valid_lpa['penetrance'], valid_lpa['median_lpa'])
    print(f"  Domain penetrance vs median Lp(a): rho={rho:.3f}, P={p:.4f} (n={len(valid_lpa)} domains)")

print(f"\n  {'Domain':>25s}  {'N Lpa':>6s}  {'Med Lpa':>8s}  {'Lpa>125%':>8s}  {'Scl%':>5s}  {'Pen%':>5s}")
for _, r in dom_lpa.sort_values('median_lpa', ascending=False).iterrows():
    pen = f"{r['penetrance']:.1f}" if pd.notna(r['penetrance']) else 'N/A'
    print(f"  {r['ldlr_domain']:>25s}  {int(r['n_lpa']):>6d}  {r['median_lpa']:>8.0f}  "
          f"{r['lpa_high_pct']:>8.1f}  {r['sclerosis_pct']:>5.1f}  {pen:>5s}")

# Sclerosis prevalence table
print(f"\n  --- Aortic Sclerosis Prevalence by Domain ---")
scl_out = scl_by_dom[scl_by_dom['n_total'] >= 20].sort_values('sclerosis_pct', ascending=False)
scl_out['penetrance'] = scl_out['ldlr_domain'].map(pen_map)
print(f"  {'Domain':>25s}  {'N':>5s}  {'Scl':>4s}  {'%':>5s}  {'Pen%':>5s}")
for _, r in scl_out.iterrows():
    pen = f"{r['penetrance']:.1f}" if pd.notna(r['penetrance']) else 'N/A'
    print(f"  {r['ldlr_domain']:>25s}  {int(r['n_total']):>5d}  {int(r['n_sclerosis']):>4d}  "
          f"{r['sclerosis_pct']:>5.1f}  {pen:>5s}")

scl_out.to_csv(os.path.join(ANALYSIS, "aortic_sclerosis_prevalence.csv"), index=False)

# ============================================================================
# STEP 6: PUBLICATION FIGURE (6-panel)
# ============================================================================
print("\n[6/6] Generating publication figure...")

fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# Colour scheme
RED = '#b2182b'
ORANGE = '#ef8a62'
BLUE = '#2166ac'
LIGHT_BLUE = '#d1e5f0'

# ── Panel A: Aortic distensibility by domain ──
ax = axes[0, 0]
dcmr = domain_ao.sort_values('mean_asc_dist')
if len(dcmr) > 2:
    colours = []
    for _, r in dcmr.iterrows():
        p = r.get('penetrance', 0)
        if pd.notna(p) and p > 25: colours.append(RED)
        elif pd.notna(p) and p > 15: colours.append(ORANGE)
        else: colours.append(BLUE)
    bars = ax.barh(range(len(dcmr)), dcmr['mean_asc_dist'], color=colours,
                   edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(dcmr)))
    ax.set_yticklabels([f"{r['ldlr_domain'].replace('Ligand-binding ','LB-')}\n(n={int(r['n'])})"
                        for _, r in dcmr.iterrows()], fontsize=7, fontfamily='Arial')
    # Error bars (SEM)
    sems = dcmr['std_asc_dist'] / np.sqrt(dcmr['n'])
    ax.errorbar(dcmr['mean_asc_dist'].values, range(len(dcmr)), xerr=sems.values,
                fmt='none', ecolor='black', capsize=2, linewidth=0.8)
    ax.set_xlabel('Mean Asc. Aortic Distensibility\n(10$^{-3}$ mmHg$^{-1}$)', fontsize=10, fontfamily='Arial')
    ax.set_title('A. Aortic Distensibility by Domain', fontsize=12, fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ── Panel B: Aortic sclerosis prevalence by domain ──
ax = axes[0, 1]
scl_plot = scl_by_dom[scl_by_dom['n_total'] >= 20].sort_values('sclerosis_pct', ascending=False).head(12)
if len(scl_plot) > 0:
    scl_plot['penetrance'] = scl_plot['ldlr_domain'].map(pen_map)
    colours_b = []
    for _, r in scl_plot.iterrows():
        p = r.get('penetrance', 0)
        if pd.notna(p) and p > 25: colours_b.append(RED)
        elif pd.notna(p) and p > 15: colours_b.append(ORANGE)
        else: colours_b.append(BLUE)
    bars_b = ax.bar(range(len(scl_plot)), scl_plot['sclerosis_pct'], color=colours_b,
                    edgecolor='black', linewidth=0.5)
    ax.set_xticks(range(len(scl_plot)))
    ax.set_xticklabels([d.replace('Ligand-binding ', 'LB-').replace('EGF-like ', 'EGF-')
                         for d in scl_plot['ldlr_domain']], rotation=45, ha='right', fontsize=7)
    ax.set_ylabel('Aortic Sclerosis (%)', fontsize=10, fontfamily='Arial')
    ax.set_title('B. Aortic Sclerosis by Domain', fontsize=12, fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ── Panel C: Aortic distensibility by mechanism (boxplot) ──
ax = axes[0, 2]
mech_data = []
mech_labels = []
mech_cols = []
for mt, col_b in [('Type 1: Structural', RED), ('Type 2: Functional', ORANGE), ('Benign', BLUE)]:
    vals = df[df['mech_simple'] == mt]['ao_asc_dist'].dropna()
    if len(vals) >= 10:
        mech_data.append(vals.values)
        short = mt.split(':')[1].strip() if ':' in mt else mt
        mech_labels.append(f"{short}\n(n={len(vals)})")
        mech_cols.append(col_b)

if mech_data:
    bp = ax.boxplot(mech_data, labels=mech_labels, patch_artist=True,
                    medianprops=dict(color='black', linewidth=2), showfliers=False,
                    whiskerprops=dict(linewidth=1), capprops=dict(linewidth=1))
    for patch, c in zip(bp['boxes'], mech_cols):
        patch.set_facecolor(c)
        patch.set_alpha(0.7)
    ax.set_ylabel('Asc. Aortic Distensibility\n(10$^{-3}$ mmHg$^{-1}$)', fontsize=10, fontfamily='Arial')
    ax.set_title('C. Distensibility by Mechanism', fontsize=12, fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ── Panel D: Penetrance vs distensibility (scatter) ──
ax = axes[1, 0]
vc = domain_ao.dropna(subset=['mean_asc_dist', 'penetrance'])
if len(vc) >= 4:
    ax.scatter(vc['penetrance'], vc['mean_asc_dist'], s=vc['n'] * 3,
               c=RED, alpha=0.7, edgecolors='black', linewidth=0.5)
    for _, r in vc.iterrows():
        short = r['ldlr_domain'].replace('Ligand-binding ', 'LB-').replace('EGF-like ', 'EGF-')
        ax.annotate(short, (r['penetrance'], r['mean_asc_dist']),
                    fontsize=6, ha='center', va='bottom', fontfamily='Arial')
    rho_d, p_d = stats.spearmanr(vc['penetrance'], vc['mean_asc_dist'])
    # Trend line
    z = np.polyfit(vc['penetrance'], vc['mean_asc_dist'], 1)
    x_line = np.linspace(vc['penetrance'].min(), vc['penetrance'].max(), 50)
    ax.plot(x_line, np.polyval(z, x_line), '--', color='grey', linewidth=1, alpha=0.7)
    ax.set_xlabel('ASCVD Penetrance by Age 60 (%)', fontsize=10, fontfamily='Arial')
    ax.set_ylabel('Mean Asc. Distensibility', fontsize=10, fontfamily='Arial')
    ax.set_title(f'D. Penetrance vs Distensibility\n(rho={rho_d:.3f}, P={p_d:.3f})',
                 fontsize=12, fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ── Panel E: Lp(a) by mechanism (boxplot) ──
ax = axes[1, 1]
lpa_data = []
lpa_labels = []
lpa_cols_e = []
for mt, col_b in [('Type 1: Structural', RED), ('Type 2: Functional', ORANGE), ('Benign', BLUE)]:
    vals = df[df['mech_simple'] == mt]['lpa'].dropna()
    if len(vals) >= 5:
        lpa_data.append(vals.values)
        short = mt.split(':')[1].strip() if ':' in mt else mt
        lpa_labels.append(f"{short}\n(n={len(vals)})")
        lpa_cols_e.append(col_b)

if lpa_data:
    bp2 = ax.boxplot(lpa_data, labels=lpa_labels, patch_artist=True,
                     medianprops=dict(color='black', linewidth=2), showfliers=False,
                     whiskerprops=dict(linewidth=1), capprops=dict(linewidth=1))
    for patch, c in zip(bp2['boxes'], lpa_cols_e):
        patch.set_facecolor(c)
        patch.set_alpha(0.7)
    ax.axhline(y=125, color='red', linestyle='--', linewidth=0.8, alpha=0.6, label='High risk (>125)')
    ax.axhline(y=75, color='orange', linestyle='--', linewidth=0.8, alpha=0.6, label='Elevated (>75)')
    ax.set_ylabel('Lp(a) (nmol/L)', fontsize=10, fontfamily='Arial')
    ax.set_title('E. Lp(a) by Mutation Mechanism', fontsize=12, fontweight='bold', fontfamily='Arial')
    ax.legend(fontsize=7, loc='upper right', framealpha=0.8)
else:
    # Fallback: Lp(a) vs aortic sclerosis across all carriers
    lpa_scl_v = df[df['lpa'].notna()].copy()
    if len(lpa_scl_v) >= 20:
        scl_yes = lpa_scl_v[lpa_scl_v['ao_sclerosis'] == 1]['lpa'].dropna()
        scl_no = lpa_scl_v[lpa_scl_v['ao_sclerosis'] == 0]['lpa'].dropna()
        bp_data = [scl_no.values, scl_yes.values] if len(scl_yes) >= 3 else [scl_no.values]
        bp_labs = [f'No Sclerosis\n(n={len(scl_no)})', f'Sclerosis\n(n={len(scl_yes)})'] if len(scl_yes) >= 3 else [f'No Sclerosis\n(n={len(scl_no)})']
        bp2 = ax.boxplot(bp_data, labels=bp_labs, patch_artist=True,
                         medianprops=dict(color='black', linewidth=2), showfliers=False)
        cols_bp = [BLUE, RED]
        for i, patch in enumerate(bp2['boxes']):
            patch.set_facecolor(cols_bp[i])
            patch.set_alpha(0.7)
        ax.axhline(y=125, color='red', linestyle='--', linewidth=0.8, alpha=0.6)
        ax.set_ylabel('Lp(a) (nmol/L)', fontsize=10, fontfamily='Arial')
        ax.set_title('E. Lp(a) by Aortic Sclerosis', fontsize=12, fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ── Panel F: Summary heatmap ──
ax = axes[1, 2]
heatmap_data = domain_ao[domain_ao['n'] >= 5].copy()
if len(heatmap_data) >= 3:
    # Normalise each column to 0-1 for heatmap
    heat_cols = ['mean_asc_dist', 'sclerosis_pct', 'median_lpa', 'lpa_high_pct', 'mean_lvm']
    heat_labels = ['Asc Dist', 'Sclerosis %', 'Lp(a)', 'Lp(a)>125 %', 'LV Mass']

    heat_matrix = []
    for col in heat_cols:
        vals = heatmap_data[col].values.astype(float)
        if np.nanmax(vals) != np.nanmin(vals):
            normed = (vals - np.nanmin(vals)) / (np.nanmax(vals) - np.nanmin(vals))
        else:
            normed = np.zeros_like(vals)
        heat_matrix.append(normed)

    heat_matrix = np.array(heat_matrix)
    domain_labels = [d.replace('Ligand-binding ', 'LB-').replace('EGF-like ', 'EGF-')
                     for d in heatmap_data['ldlr_domain']]

    im = ax.imshow(heat_matrix, cmap='RdYlBu_r', aspect='auto')
    ax.set_xticks(range(len(domain_labels)))
    ax.set_xticklabels(domain_labels, rotation=45, ha='right', fontsize=7)
    ax.set_yticks(range(len(heat_labels)))
    ax.set_yticklabels(heat_labels, fontsize=9, fontfamily='Arial')
    ax.set_title('F. Multi-Modal Aortic Summary', fontsize=12, fontweight='bold', fontfamily='Arial')

    # Annotate cells with actual values
    for i, col in enumerate(heat_cols):
        for j, val in enumerate(heatmap_data[col].values):
            if pd.notna(val):
                fmt = f"{val:.1f}" if col != 'sclerosis_pct' else f"{val:.1f}%"
                ax.text(j, i, fmt, ha='center', va='center', fontsize=6,
                        color='white' if heat_matrix[i, j] > 0.6 else 'black')

    plt.colorbar(im, ax=ax, shrink=0.6, label='Normalised')

plt.tight_layout(pad=2.0)
fig_path = os.path.join(FIGURES, "Figure_Aortic_Phenotype_Panel.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {fig_path}")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("COMPLETE")
print("=" * 80)
n_with_ao = df['ao_asc_dist'].notna().sum()
n_with_scl = df['ao_sclerosis'].sum()
n_domains = len(domain_ao)
print(f"""
  AORTIC PHENOTYPE PANEL:
    {len(df)} FH carriers total
    {n_with_ao} with aortic distensibility data
    {n_with_scl} with aortic sclerosis (ICD-10 I35)
    {n_domains} domains with sufficient data (n>=5)

  ANALYSES:
    1. Aortic distensibility by domain (ascending + descending)
    2. Aortic sclerosis prevalence by domain
    3. Mechanism comparison (T1 vs T2 vs Benign)
    4. Correlations: penetrance, SSS, ddG vs distensibility
    5. Structure-function: distensibility vs LV mass

  FILES SAVED:
    aortic_phenotype_by_domain.csv
    aortic_phenotype_by_mechanism.csv
    aortic_sclerosis_prevalence.csv
    Figure_Aortic_Phenotype_Panel.png
""")

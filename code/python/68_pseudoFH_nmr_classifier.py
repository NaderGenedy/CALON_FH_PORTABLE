#!/usr/bin/env python3
"""
68_pseudoFH_nmr_classifier.py
===============================
Nature Gap 4: Pseudo-FH Classifier

Can NMR metabolomics distinguish "true monogenic FH" from
"Lp(a)-driven phenocopy" (pseudo-FH)?

Approach:
  1. Define pseudo-FH: high Lp(a) + LDL inflated by Lp(a)-cholesterol
  2. Compare NMR metabolomic profiles: true FH vs pseudo-FH vs non-FH
  3. Identify discriminating metabolites
  4. Build classifier and test AUC

Key metabolites from Nightingale NMR platform:
  - Lipoprotein subclasses (VLDL, IDL, LDL, HDL particle sizes)
  - ApoB, ApoA1
  - Fatty acids (SFA, MUFA, PUFA, omega-3, omega-6)
  - GlycA (inflammation)
  - Amino acids
  - Glycolysis intermediates

Author: Dr Nader Genedy
Date:   April 2026
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

ANALYSIS = r"C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis"
FIGURES = os.path.join(ANALYSIS, "figures")
DATA = r"D:/CALON_AF3_PROJECT/data"

print("=" * 70)
print("GAP 4: PSEUDO-FH NMR CLASSIFIER")
print("Can metabolomics distinguish true FH from Lp(a)-phenocopy?")
print("=" * 70)

# ============================================================================
# STEP 1: LOAD NMR METABOLOMICS FOR FH CARRIERS
# ============================================================================
print("\n[1/5] Loading NMR metabolomics...")

# Key NMR fields from Nightingale platform
# Using calon_batch_nmr2d which has GlycA (p23480)
# and calon_extra_nmr files which have lipoprotein subclasses

# GlycA
glyca = pd.read_csv(f"{DATA}/ukb_wide/calon_batch_nmr2d.csv",
                     usecols=['participant.eid', 'participant.p23480_i0'])
glyca.columns = ['eid', 'glyca']
glyca['eid'] = glyca['eid'].astype(str)
glyca['glyca'] = pd.to_numeric(glyca['glyca'], errors='coerce')

# Load a broader set of NMR metabolites
# calon_batch_nmr1a has lipoprotein subclasses
nmr1a = pd.read_csv(f"{DATA}/ukb_wide/calon_batch_nmr1a.csv")
nmr1a.columns = [c.replace('participant.', '') for c in nmr1a.columns]
nmr1a['eid'] = nmr1a['eid'].astype(str)

# Get column names for NMR metabolites
nmr_cols = [c for c in nmr1a.columns if c.startswith('p23') and c != 'eid']
print(f"  NMR metabolites in batch 1a: {len(nmr_cols)}")

# Convert to numeric
for c in nmr_cols:
    nmr1a[c] = pd.to_numeric(nmr1a[c], errors='coerce')

# Load Lp(a)
lpa = pd.read_csv(f"{DATA}/ukb_wide/calon_lpa_CORRECT.csv")
lpa.columns = [c.replace('participant.', '') for c in lpa.columns]
lpa['eid'] = lpa['eid'].astype(str)
lpa['lpa'] = pd.to_numeric(lpa['p30790_i0'], errors='coerce')

# Load LDL
ldl = pd.read_csv(f"{DATA}/ukb_wide/ukb_reviewer_longitudinal_lipids.csv",
                   usecols=['participant.eid', 'participant.p30780_i0'])
ldl.columns = ['eid', 'ldl']
ldl['eid'] = ldl['eid'].astype(str)
ldl['ldl'] = pd.to_numeric(ldl['ldl'], errors='coerce')

# Load FH carriers
carriers = pd.read_csv(os.path.join(ANALYSIS, "ukb_carriers_annotated.csv"),
                        dtype={'eid': str})
patients = carriers.sort_values('sss', ascending=False).drop_duplicates('eid')
fh_eids = set(patients['eid'])

# ICD-10 ASCVD
icd = pd.read_csv(f"{DATA}/ukb_wide/calon_batch_icd10.csv",
                   usecols=['participant.eid', 'participant.p41270'], dtype=str)
icd.columns = ['eid', 'icd10_raw']
icd['ascvd'] = icd['icd10_raw'].str.contains('"I2[0-5]', na=False, regex=True).astype(int)

print(f"  GlycA: {glyca['glyca'].notna().sum()}")
print(f"  Lp(a): {lpa['lpa'].notna().sum()}")

# ============================================================================
# STEP 2: DEFINE PSEUDO-FH vs TRUE FH
# ============================================================================
print("\n[2/5] Defining pseudo-FH groups...")

# Merge all for FH carriers
fh = patients[['eid', 'ldlr_domain', 'gene']].copy()
fh = fh.merge(lpa[['eid', 'lpa']], on='eid', how='left')
fh = fh.merge(ldl[['eid', 'ldl']], on='eid', how='left')
fh = fh.merge(glyca[['eid', 'glyca']], on='eid', how='left')
fh = fh.merge(nmr1a, on='eid', how='left')
fh = fh.merge(icd[['eid', 'ascvd']], on='eid', how='left')

# Lp(a)-cholesterol correction
fh['lpa_chol'] = (fh['lpa'] / 2.5) * 0.30 / 38.7
fh['ldl_corrected'] = fh['ldl'] - fh['lpa_chol']

# Define groups:
# Pseudo-FH: high Lp(a) AND Lp(a)-cholesterol contributes >5% of LDL
# True FH: low-moderate Lp(a) (Lp(a) not driving the phenotype)
fh['lpa_pct_of_ldl'] = 100 * fh['lpa_chol'] / fh['ldl']

has_all = fh[fh['lpa'].notna() & fh['ldl'].notna()].copy()

# Pseudo-FH: Lp(a) > 75 nmol/L AND contributes >5% of LDL
has_all['pseudo_fh'] = ((has_all['lpa'] > 75) & (has_all['lpa_pct_of_ldl'] > 5)).astype(int)
# True FH: Lp(a) < 30 nmol/L (clearly not Lp(a)-driven)
has_all['true_fh'] = (has_all['lpa'] < 30).astype(int)
# Intermediate: between 30-75
has_all['group'] = 'Intermediate'
has_all.loc[has_all['pseudo_fh'] == 1, 'group'] = 'Pseudo-FH'
has_all.loc[has_all['true_fh'] == 1, 'group'] = 'True FH'

print(f"  Total FH with Lp(a) + LDL + NMR: {len(has_all)}")
for g in ['True FH', 'Intermediate', 'Pseudo-FH']:
    sub = has_all[has_all['group'] == g]
    print(f"    {g:>15s}: n={len(sub)}, Lp(a)={sub['lpa'].median():.0f}, "
          f"LDL={sub['ldl'].median():.2f}, ASCVD={100*sub['ascvd'].mean():.1f}%")

# ============================================================================
# STEP 3: COMPARE NMR PROFILES
# ============================================================================
print(f"\n{'='*70}")
print("[3/5] NMR METABOLOMIC COMPARISON: True FH vs Pseudo-FH")
print("=" * 70)

# Find metabolites that differ between True FH and Pseudo-FH
available_nmr = [c for c in nmr_cols if has_all[c].notna().sum() > 100]
print(f"  Available NMR metabolites with >100 values: {len(available_nmr)}")

results = []
for col in available_nmr:
    true_vals = has_all[has_all['group'] == 'True FH'][col].dropna()
    pseudo_vals = has_all[has_all['group'] == 'Pseudo-FH'][col].dropna()

    if len(true_vals) >= 20 and len(pseudo_vals) >= 20:
        u, p = stats.mannwhitneyu(true_vals, pseudo_vals, alternative='two-sided')
        d = (true_vals.mean() - pseudo_vals.mean()) / np.sqrt((true_vals.std()**2 + pseudo_vals.std()**2) / 2)
        results.append({
            'metabolite': col,
            'true_fh_median': true_vals.median(),
            'pseudo_fh_median': pseudo_vals.median(),
            'p_value': p,
            'cohen_d': d,
            'pct_diff': 100 * (pseudo_vals.median() - true_vals.median()) / max(abs(true_vals.median()), 0.001),
        })

results_df = pd.DataFrame(results).sort_values('p_value')

# Show top discriminating metabolites
print(f"\n  Top discriminating metabolites (True FH vs Pseudo-FH):")
print(f"  {'Metabolite':>15s}  {'True FH':>9s}  {'Pseudo':>9s}  {'% Diff':>8s}  {'P-value':>10s}  {'d':>6s}")
print("  " + "-" * 65)

# UKB NMR field name mapping (subset)
field_names = {
    'p23400_i0': 'Total-C', 'p23401_i0': 'VLDL-C', 'p23402_i0': 'Remnant-C',
    'p23403_i0': 'LDL-C(NMR)', 'p23404_i0': 'HDL-C(NMR)', 'p23405_i0': 'VLDL-TG',
    'p23406_i0': 'LDL-TG', 'p23407_i0': 'HDL-TG', 'p23408_i0': 'Total-TG',
    'p23409_i0': 'Total-PL', 'p23410_i0': 'VLDL-PL', 'p23411_i0': 'LDL-PL',
    'p23412_i0': 'HDL-PL', 'p23413_i0': 'ApoB(NMR)', 'p23414_i0': 'ApoA1(NMR)',
}

for _, r in results_df.head(20).iterrows():
    name = field_names.get(r['metabolite'], r['metabolite'])
    sig = '***' if r['p_value'] < 0.001 else '**' if r['p_value'] < 0.01 else '*' if r['p_value'] < 0.05 else ''
    print(f"  {name:>15s}  {r['true_fh_median']:>9.3f}  {r['pseudo_fh_median']:>9.3f}  "
          f"{r['pct_diff']:>+7.1f}%  {r['p_value']:>10.4f}  {r['cohen_d']:>6.3f} {sig}")

n_sig = (results_df['p_value'] < 0.05).sum()
n_bonf = (results_df['p_value'] < 0.05 / len(results_df)).sum()
print(f"\n  Significant (P<0.05): {n_sig}/{len(results_df)}")
print(f"  Bonferroni-corrected: {n_bonf}/{len(results_df)}")

# ============================================================================
# STEP 4: BUILD CLASSIFIER
# ============================================================================
print(f"\n{'='*70}")
print("[4/5] BUILD NMR CLASSIFIER: True FH vs Pseudo-FH")
print("=" * 70)

# Use top metabolites + Lp(a) for classification
# Exclude True FH vs Pseudo-FH only (drop intermediate)
classify = has_all[has_all['group'].isin(['True FH', 'Pseudo-FH'])].copy()
classify['target'] = (classify['group'] == 'Pseudo-FH').astype(int)

# Select features: top significant metabolites + Lp(a)
sig_mets = results_df[results_df['p_value'] < 0.05]['metabolite'].tolist()[:10]
features = sig_mets + ['lpa', 'glyca']
features = [f for f in features if f in classify.columns]

# Drop rows with NaN in features
classify_clean = classify[features + ['target', 'ascvd']].dropna()
print(f"  Classifier data: {len(classify_clean)} (True FH: {(classify_clean['target']==0).sum()}, "
      f"Pseudo-FH: {(classify_clean['target']==1).sum()})")

if len(classify_clean) >= 50 and classify_clean['target'].sum() >= 10:
    X = classify_clean[features].values
    y = classify_clean['target'].values

    # Standardise
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Logistic regression
    lr = LogisticRegression(max_iter=1000, solver='lbfgs')
    lr.fit(X_scaled, y)
    pred = lr.predict_proba(X_scaled)[:, 1]
    auc = roc_auc_score(y, pred)

    print(f"\n  NMR Classifier AUC: {auc:.3f}")
    print(f"  Features used: {len(features)}")
    print(f"\n  Feature importance (OR per SD):")
    for feat, coef in sorted(zip(features, lr.coef_[0]), key=lambda x: abs(x[1]), reverse=True):
        name = field_names.get(feat, feat)
        print(f"    {name:>15s}: OR={np.exp(coef):.3f} (beta={coef:.3f})")

    # Lp(a) alone vs NMR classifier
    auc_lpa_alone = roc_auc_score(y, classify_clean['lpa'].values)
    print(f"\n  Comparison:")
    print(f"    Lp(a) alone AUC:     {auc_lpa_alone:.3f}")
    print(f"    NMR classifier AUC:  {auc:.3f}")
    print(f"    Improvement:         {auc - auc_lpa_alone:+.3f}")

    # Does pseudo-FH have different ASCVD?
    print(f"\n  ASCVD by group:")
    for g in ['True FH', 'Pseudo-FH']:
        sub = has_all[has_all['group'] == g]
        print(f"    {g}: {100*sub['ascvd'].mean():.1f}% (n={len(sub)})")

# ============================================================================
# STEP 5: NON-FH COMPARISON — NMR signature of high Lp(a)
# ============================================================================
print(f"\n{'='*70}")
print("[5/5] NON-FH: Does high Lp(a) create a distinct NMR signature?")
print("=" * 70)

# Compare NMR in high vs low Lp(a) non-FH participants
nonfh = nmr1a[~nmr1a['eid'].isin(fh_eids)].copy()
nonfh = nonfh.merge(lpa[['eid', 'lpa']], on='eid', how='inner')
nonfh = nonfh[nonfh['lpa'].notna()].copy()

nonfh['lpa_group'] = 'Low'
nonfh.loc[nonfh['lpa'] > 75, 'lpa_group'] = 'High'

print(f"  Non-FH with NMR + Lp(a): {len(nonfh)}")
print(f"  High Lp(a) (>75): {(nonfh['lpa_group']=='High').sum()}")

nonfh_results = []
for col in available_nmr[:20]:  # Top 20 metabolites
    hi = nonfh[nonfh['lpa_group'] == 'High'][col].dropna()
    lo = nonfh[nonfh['lpa_group'] == 'Low'][col].dropna()
    if len(hi) >= 100 and len(lo) >= 100:
        u, p = stats.mannwhitneyu(hi, lo)
        d = (hi.mean() - lo.mean()) / np.sqrt((hi.std()**2 + lo.std()**2) / 2)
        nonfh_results.append({
            'metabolite': col,
            'high_lpa': hi.median(),
            'low_lpa': lo.median(),
            'p_value': p,
            'cohen_d': d,
        })

nonfh_df = pd.DataFrame(nonfh_results).sort_values('p_value')
print(f"\n  Top metabolites differing by Lp(a) (non-FH):")
for _, r in nonfh_df.head(10).iterrows():
    name = field_names.get(r['metabolite'], r['metabolite'])
    print(f"    {name:>15s}: d={r['cohen_d']:.3f}, P={r['p_value']:.2e}")

# ============================================================================
# FIGURE
# ============================================================================
print(f"\n{'='*70}")
print("Generating figure...")

fig, axes = plt.subplots(2, 2, figsize=(14, 11))
RED = '#b2182b'
BLUE = '#2166ac'
ORANGE = '#ef8a62'

# Panel A: Lp(a) distribution by group
ax = axes[0, 0]
for g, col in [('True FH', BLUE), ('Intermediate', '#999999'), ('Pseudo-FH', RED)]:
    vals = has_all[has_all['group'] == g]['lpa'].dropna()
    if len(vals) >= 10:
        ax.hist(vals, bins=30, alpha=0.5, color=col, label=f'{g} (n={len(vals)})',
                density=True, edgecolor='white', linewidth=0.3)
ax.axvline(x=75, color='black', linestyle='--', linewidth=1)
ax.set_xlabel('Lp(a) (nmol/L)', fontsize=10)
ax.set_ylabel('Density', fontsize=10)
ax.set_title('A. Lp(a) Distribution by FH Group', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: Top discriminating metabolites (effect size)
ax = axes[0, 1]
if len(results_df) > 0:
    top10 = results_df.head(10).copy()
    top10['name'] = top10['metabolite'].map(field_names).fillna(top10['metabolite'])
    top10 = top10.sort_values('cohen_d')
    colours = [RED if d > 0 else BLUE for d in top10['cohen_d']]
    ax.barh(range(len(top10)), top10['cohen_d'], color=colours,
            edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(top10)))
    ax.set_yticklabels(top10['name'], fontsize=8)
    ax.axvline(x=0, color='black', linewidth=0.5)
    ax.set_xlabel("Cohen's d (Pseudo-FH vs True FH)", fontsize=10)
    ax.set_title('B. Top Discriminating Metabolites', fontsize=12, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: ROC curve (NMR classifier vs Lp(a) alone)
ax = axes[1, 0]
if len(classify_clean) >= 50 and classify_clean['target'].sum() >= 10:
    fpr_nmr, tpr_nmr, _ = roc_curve(y, pred)
    fpr_lpa, tpr_lpa, _ = roc_curve(y, classify_clean['lpa'].values)
    ax.plot(fpr_nmr, tpr_nmr, color=RED, linewidth=2, label=f'NMR Classifier (AUC={auc:.3f})')
    ax.plot(fpr_lpa, tpr_lpa, color=BLUE, linewidth=2, label=f'Lp(a) Alone (AUC={auc_lpa_alone:.3f})')
    ax.plot([0, 1], [0, 1], '--', color='grey', linewidth=0.5)
    ax.set_xlabel('1 - Specificity', fontsize=10)
    ax.set_ylabel('Sensitivity', fontsize=10)
    ax.set_title('C. Pseudo-FH Detection: NMR vs Lp(a)', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: ASCVD by group
ax = axes[1, 1]
groups_plot = []
for g, col in [('True FH', BLUE), ('Intermediate', '#999999'), ('Pseudo-FH', RED)]:
    sub = has_all[has_all['group'] == g]
    if len(sub) >= 10:
        groups_plot.append((g, 100 * sub['ascvd'].mean(), len(sub), col))

if groups_plot:
    bars = ax.bar(range(len(groups_plot)), [g[1] for g in groups_plot],
                  color=[g[3] for g in groups_plot], edgecolor='black', linewidth=0.5)
    ax.set_xticks(range(len(groups_plot)))
    ax.set_xticklabels([f"{g[0]}\n(n={g[2]})" for g in groups_plot], fontsize=9)
    ax.set_ylabel('ASCVD Rate (%)', fontsize=10)
    ax.set_title('D. ASCVD by FH Classification', fontsize=12, fontweight='bold')
    for bar, g in zip(bars, groups_plot):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                f'{g[1]:.1f}%', ha='center', fontsize=10, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout(pad=2.0)
fig_path = os.path.join(FIGURES, "Figure_PseudoFH_NMR.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {fig_path}")

# Save results
results_df.to_csv(os.path.join(ANALYSIS, "pseudofh_nmr_discriminating_metabolites.csv"), index=False)

print("\n" + "=" * 70)
print("COMPLETE")
print("=" * 70)

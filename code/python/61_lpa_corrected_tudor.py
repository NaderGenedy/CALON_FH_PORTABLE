#!/usr/bin/env python3
"""
61_lpa_corrected_tudor.py
==========================
Impact of Lp(a)-cholesterol correction on TUDOR diagnostic performance.

Key question: If ~10-25% of "clinical FH" is driven by Lp(a)-inflated LDL-C,
how does Lp(a) correction affect:
  1. TUDOR's AUC for FH diagnosis
  2. DLCN score reclassification
  3. Sensitivity/specificity at optimal threshold
  4. Net reclassification improvement (NRI)

Uses Wales FH cohort (1,308 with Lp(a)) and UKB lipid clinic cohort.

Author: Dr Nader Genedy
Date:   April 2026
"""

import pandas as pd
import numpy as np
from scipy import stats
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
print("Lp(a)-CORRECTED TUDOR DIAGNOSTIC PERFORMANCE")
print("=" * 70)

# ============================================================================
# STEP 1: Load Wales FH cohort with Lp(a)
# ============================================================================
print("\n[1/5] Loading Wales FH cohort...")

ar = pd.read_csv("D:/CALON_AF3_PROJECT/data/fh_cohort/calon_ukb_analysis_ready.csv")
ar['lpa'] = pd.to_numeric(ar['lpa'], errors='coerce')
ar['ldl'] = pd.to_numeric(ar['ldl'], errors='coerce')
ar['tc'] = pd.to_numeric(ar['tc'], errors='coerce')

has_lpa = ar[ar['lpa'].notna() & ar['ldl'].notna()].copy()
print(f"  Wales FH with Lp(a) + LDL: {len(has_lpa)}")
print(f"  Lp(a): median={has_lpa['lpa'].median():.1f} nmol/L")
print(f"  LDL-C: median={has_lpa['ldl'].median():.2f} mmol/L")

# ============================================================================
# STEP 2: Compute Lp(a)-cholesterol correction
# ============================================================================
print("\n[2/5] Computing Lp(a)-cholesterol correction...")

# Conversion: Lp(a) nmol/L -> mg/dL: divide by 2.5 (approximate)
# Lp(a)-cholesterol (mmol/L) = Lp(a)_mg/dL * 0.30 / 38.7
# 0.30 = fraction of Lp(a) mass that is cholesterol
# 38.7 = molecular weight conversion factor for cholesterol (mg/dL to mmol/L)
has_lpa['lpa_mgdl'] = has_lpa['lpa'] / 2.5
has_lpa['lpa_chol_mmol'] = has_lpa['lpa_mgdl'] * 0.30 / 38.7
has_lpa['ldl_corrected'] = has_lpa['ldl'] - has_lpa['lpa_chol_mmol']

# Ensure corrected LDL doesn't go below zero
has_lpa['ldl_corrected'] = has_lpa['ldl_corrected'].clip(lower=0.5)

print(f"  Lp(a)-cholesterol: mean={has_lpa['lpa_chol_mmol'].mean():.3f} mmol/L")
print(f"  As % of measured LDL: {100 * has_lpa['lpa_chol_mmol'].mean() / has_lpa['ldl'].mean():.1f}%")
print(f"  LDL raw: median={has_lpa['ldl'].median():.2f} -> corrected: {has_lpa['ldl_corrected'].median():.2f}")

# Distribution of correction magnitude
print(f"\n  Correction magnitude distribution:")
corr = has_lpa['lpa_chol_mmol']
for pct in [25, 50, 75, 90, 95, 99]:
    print(f"    p{pct}: {corr.quantile(pct/100):.3f} mmol/L")

# ============================================================================
# STEP 3: Impact on DLCN diagnostic thresholds
# ============================================================================
print(f"\n{'='*70}")
print("[3/5] DLCN THRESHOLD RECLASSIFICATION")
print("=" * 70)

# DLCN LDL-C scoring:
# >8.5 mmol/L (>330 mg/dL): 8 points
# 6.5-8.4 mmol/L (250-329): 5 points
# 5.0-6.4 mmol/L (190-249): 3 points
# 4.0-4.9 mmol/L (155-189): 1 point
thresholds = [(8.5, '8 points (definite FH)'),
              (6.5, '5 points (probable FH)'),
              (5.0, '3 points (possible FH)'),
              (4.0, '1 point (elevated)')]

print(f"\n  {'Threshold':>30s}  {'Raw':>6s}  {'Corrected':>10s}  {'Reclassified':>13s}  {'% Drop':>7s}")
print("  " + "-" * 75)

for thresh, label in thresholds:
    above_raw = (has_lpa['ldl'] > thresh).sum()
    above_corr = (has_lpa['ldl_corrected'] > thresh).sum()
    reclass = above_raw - above_corr
    pct_drop = 100 * reclass / max(above_raw, 1)
    print(f"  LDL>{thresh} ({label:>25s})  {above_raw:>6d}  {above_corr:>10d}  {reclass:>13d}  {pct_drop:>6.1f}%")

# DLCN score recalculation
def dlcn_ldl_points(ldl):
    if ldl > 8.5: return 8
    elif ldl > 6.5: return 5
    elif ldl > 5.0: return 3
    elif ldl > 4.0: return 1
    else: return 0

has_lpa['dlcn_raw'] = has_lpa['ldl'].apply(dlcn_ldl_points)
has_lpa['dlcn_corr'] = has_lpa['ldl_corrected'].apply(dlcn_ldl_points)
has_lpa['dlcn_change'] = has_lpa['dlcn_corr'] - has_lpa['dlcn_raw']

downgraded = (has_lpa['dlcn_change'] < 0).sum()
unchanged = (has_lpa['dlcn_change'] == 0).sum()
print(f"\n  DLCN LDL-C point changes:")
print(f"    Downgraded: {downgraded} ({100*downgraded/len(has_lpa):.1f}%)")
print(f"    Unchanged:  {unchanged} ({100*unchanged/len(has_lpa):.1f}%)")

# Who gets downgraded the most?
dg = has_lpa[has_lpa['dlcn_change'] < 0]
if len(dg) > 0:
    print(f"    Downgraded patients: mean Lp(a)={dg['lpa'].mean():.0f} nmol/L "
          f"(vs {has_lpa['lpa'].mean():.0f} overall)")
    print(f"    Mean correction: {dg['lpa_chol_mmol'].mean():.3f} mmol/L")

# ============================================================================
# STEP 4: Impact on TUDOR AUC
# ============================================================================
print(f"\n{'='*70}")
print("[4/5] IMPACT ON TUDOR-LIKE DIAGNOSTIC PERFORMANCE")
print("=" * 70)

# We can test: does Lp(a) correction improve or worsen FH classification?
# In the Wales cohort, all patients ARE FH (genetically confirmed).
# So the question becomes: among the subgroup with Lp(a) data,
# does high Lp(a) explain why some patients were diagnosed?

# We can test this differently:
# 1. Correlate Lp(a)-correction magnitude with ASCVD (if Lp(a) drives risk independent of FH)
# 2. Check if corrected LDL still separates FH severity
# 3. Check impact on TUDOR's engineered features

# Test: Lp(a) contribution to LDL-C by gene
print("\n  Lp(a)-cholesterol by gene:")
for gene in ['LDLR', 'APOB']:
    sub = has_lpa[has_lpa['gene'] == gene]
    if len(sub) >= 10:
        print(f"    {gene}: n={len(sub)}, mean Lp(a)-C = {sub['lpa_chol_mmol'].mean():.3f} mmol/L "
              f"({100*sub['lpa_chol_mmol'].mean()/sub['ldl'].mean():.1f}% of LDL)")

# Test: Lp(a) correction and ASCVD
if 'ascvd_combined' in has_lpa.columns:
    ascvd = pd.to_numeric(has_lpa['ascvd_combined'], errors='coerce')
    valid = has_lpa[ascvd.notna()].copy()
    valid['ascvd'] = ascvd[valid.index]

    try:
        from sklearn.metrics import roc_auc_score

        # AUC: raw LDL for ASCVD
        auc_raw = roc_auc_score(valid['ascvd'], valid['ldl'])
        # AUC: corrected LDL for ASCVD
        auc_corr = roc_auc_score(valid['ascvd'], valid['ldl_corrected'])
        # AUC: Lp(a) alone for ASCVD
        auc_lpa = roc_auc_score(valid['ascvd'], valid['lpa'])

        print(f"\n  AUC for ASCVD prediction:")
        print(f"    Raw LDL-C:       {auc_raw:.3f}")
        print(f"    Corrected LDL-C: {auc_corr:.3f} (delta {auc_corr - auc_raw:+.3f})")
        print(f"    Lp(a) alone:     {auc_lpa:.3f}")

        # Does correction improve ASCVD prediction?
        if auc_corr > auc_raw:
            print(f"    -> Correction IMPROVES ASCVD prediction by {auc_corr - auc_raw:.3f}")
        else:
            print(f"    -> Correction does not improve (or slightly worsens) prediction")
    except ImportError:
        pass

# ============================================================================
# STEP 5: Lp(a) quintile analysis
# ============================================================================
print(f"\n{'='*70}")
print("[5/5] Lp(a) QUINTILE ANALYSIS — WHO IS MOST AFFECTED?")
print("=" * 70)

has_lpa['lpa_quintile'] = pd.qcut(has_lpa['lpa'], 5, labels=['Q1', 'Q2', 'Q3', 'Q4', 'Q5'])

print(f"\n  {'Quintile':>10s}  {'Lp(a) range':>15s}  {'Mean LDL raw':>12s}  {'Mean LDL corr':>13s}  "
      f"{'Mean correction':>15s}  {'% correction':>13s}")
print("  " + "-" * 85)

for q in ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']:
    sub = has_lpa[has_lpa['lpa_quintile'] == q]
    lpa_range = f"{sub['lpa'].min():.0f}-{sub['lpa'].max():.0f}"
    mean_raw = sub['ldl'].mean()
    mean_corr = sub['ldl_corrected'].mean()
    mean_correction = sub['lpa_chol_mmol'].mean()
    pct_corr = 100 * mean_correction / mean_raw
    print(f"  {q:>10s}  {lpa_range:>15s}  {mean_raw:>12.2f}  {mean_corr:>13.2f}  "
          f"{mean_correction:>15.3f}  {pct_corr:>12.1f}%")

# High Lp(a) subgroup analysis
print(f"\n  Impact by Lp(a) category:")
for thresh, label in [(50, 'Elevated (>50)'), (75, 'High (>75)'), (125, 'Very high (>125)')]:
    sub = has_lpa[has_lpa['lpa'] > thresh]
    if len(sub) >= 5:
        mean_corr_mag = sub['lpa_chol_mmol'].mean()
        pct = 100 * mean_corr_mag / sub['ldl'].mean()
        n_reclass = (sub['dlcn_change'] < 0).sum()
        print(f"    {label}: n={len(sub)}, correction={mean_corr_mag:.3f} mmol/L ({pct:.1f}% of LDL), "
              f"DLCN downgraded={n_reclass}")

# ============================================================================
# FIGURE: 4-panel
# ============================================================================
print(f"\n{'='*70}")
print("Generating figure...")

fig, axes = plt.subplots(2, 2, figsize=(14, 11))
RED = '#b2182b'
BLUE = '#2166ac'

# Panel A: LDL raw vs corrected scatter
ax = axes[0, 0]
ax.scatter(has_lpa['ldl'], has_lpa['ldl_corrected'], c=has_lpa['lpa'],
           cmap='RdYlBu_r', s=10, alpha=0.6)
ax.plot([0, 12], [0, 12], '--', color='grey', linewidth=1)
ax.set_xlabel('Measured LDL-C (mmol/L)', fontsize=10)
ax.set_ylabel('Lp(a)-Corrected LDL-C (mmol/L)', fontsize=10)
ax.set_title('A. LDL-C Before vs After Lp(a) Correction', fontsize=12, fontweight='bold')
cb = plt.colorbar(ax.collections[0], ax=ax, shrink=0.8)
cb.set_label('Lp(a) (nmol/L)', fontsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: Correction magnitude by Lp(a) quintile
ax = axes[0, 1]
q_data = []
q_labels = []
for q in ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']:
    sub = has_lpa[has_lpa['lpa_quintile'] == q]
    q_data.append(sub['lpa_chol_mmol'].values)
    q_labels.append(f"{q}\n({sub['lpa'].min():.0f}-{sub['lpa'].max():.0f})")

bp = ax.boxplot(q_data, labels=q_labels, patch_artist=True,
                medianprops=dict(color='black', linewidth=2), showfliers=False)
colours = ['#2166ac', '#67a9cf', '#fddbc7', '#ef8a62', '#b2182b']
for patch, c in zip(bp['boxes'], colours):
    patch.set_facecolor(c)
    patch.set_alpha(0.7)
ax.set_ylabel('Lp(a)-Cholesterol (mmol/L)', fontsize=10)
ax.set_xlabel('Lp(a) Quintile', fontsize=10)
ax.set_title('B. Correction Magnitude by Quintile', fontsize=12, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: DLCN reclassification waterfall
ax = axes[1, 0]
changes = has_lpa['dlcn_change'].value_counts().sort_index()
colours_c = [RED if x < 0 else BLUE if x > 0 else '#999999' for x in changes.index]
ax.bar(changes.index, changes.values, color=colours_c, edgecolor='black', linewidth=0.5)
ax.set_xlabel('Change in DLCN LDL-C Points', fontsize=10)
ax.set_ylabel('Number of Patients', fontsize=10)
ax.set_title(f'C. DLCN Score Change\n({downgraded} downgraded, {100*downgraded/len(has_lpa):.1f}%)',
             fontsize=12, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: % LDL from Lp(a) by Lp(a) level
ax = axes[1, 1]
has_lpa['pct_from_lpa'] = 100 * has_lpa['lpa_chol_mmol'] / has_lpa['ldl']
ax.scatter(has_lpa['lpa'], has_lpa['pct_from_lpa'], c=BLUE, alpha=0.3, s=10)
ax.axhline(y=10, color=RED, linestyle='--', linewidth=1, label='10% threshold')
ax.set_xlabel('Lp(a) (nmol/L)', fontsize=10)
ax.set_ylabel('% of LDL-C from Lp(a)', fontsize=10)
ax.set_title('D. Lp(a) Contribution to Measured LDL-C', fontsize=12, fontweight='bold')
above_10 = (has_lpa['pct_from_lpa'] > 10).sum()
ax.legend(fontsize=8, title=f'{above_10} patients >10%')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout(pad=2.0)
fig_path = os.path.join(FIGURES, "Figure_Lpa_Corrected_TUDOR.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {fig_path}")

# Save results
has_lpa[['eid', 'lpa', 'ldl', 'ldl_corrected', 'lpa_chol_mmol', 'dlcn_raw',
          'dlcn_corr', 'dlcn_change', 'gene']].to_csv(
    os.path.join(ANALYSIS, "lpa_corrected_ldl.csv"), index=False)

print("\n" + "=" * 70)
print("COMPLETE")
print("=" * 70)

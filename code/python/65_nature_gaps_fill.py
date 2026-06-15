#!/usr/bin/env python3
"""
65_nature_gaps_fill.py
=======================
Fill remaining Nature-level research gaps with available data:

1. FH-specific Lp(a) risk thresholds (vs general population)
2. Lp(a) as continuous predictor: FH vs non-FH (AUC comparison)
3. Pseudo-FH quantification: how many FH diagnoses are Lp(a)-driven?
4. Combined SSS x Lp(a) risk model
5. Domain-specific ASCVD with Lp(a) interaction

Author: Dr Nader Genedy
Date:   April 2026
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LogisticRegression
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

ANALYSIS = r"C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis"
FIGURES = os.path.join(ANALYSIS, "figures")

print("=" * 70)
print("NATURE GAPS: FH-SPECIFIC Lp(a) THRESHOLDS + COMBINED RISK")
print("=" * 70)

# Load data
carriers = pd.read_csv(os.path.join(ANALYSIS, "ukb_carriers_annotated.csv"), dtype={'eid': str})
patients = carriers.sort_values('sss', ascending=False).drop_duplicates('eid')

lpa = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/calon_lpa_CORRECT.csv")
lpa.columns = [c.replace('participant.', '') for c in lpa.columns]
lpa['eid'] = lpa['eid'].astype(str)
lpa['lpa'] = pd.to_numeric(lpa['p30790_i0'], errors='coerce')

sex = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/calon_sex.csv")
sex.columns = [c.replace('participant.', '') for c in sex.columns]
sex.rename(columns={'p31': 'sex'}, inplace=True)
sex['eid'] = sex['eid'].astype(str)

demo = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/ukb_reviewer_demographics.csv")
demo.columns = [c.replace('participant.', '') for c in demo.columns]
demo['eid'] = demo['eid'].astype(str)
demo['age'] = 2010 - pd.to_numeric(demo['p34'], errors='coerce')

icd = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/calon_batch_icd10.csv",
                   usecols=['participant.eid', 'participant.p41270'], dtype=str)
icd.columns = ['eid', 'icd10_raw']
icd['ascvd'] = icd['icd10_raw'].str.contains('"I2[0-5]', na=False, regex=True).astype(int)

# Merge FH carriers
fh = patients[['eid', 'ldlr_domain', 'mech_simple', 'gene']].copy()
fh = fh.merge(lpa[['eid', 'lpa']], on='eid', how='left')
fh = fh.merge(sex[['eid', 'sex']], on='eid', how='left')
fh = fh.merge(demo[['eid', 'age']], on='eid', how='left')
fh = fh.merge(icd[['eid', 'ascvd']], on='eid', how='left')
fh_lpa = fh[fh['lpa'].notna()].copy()

# Non-FH population
nonfh = icd[~icd['eid'].isin(set(patients['eid']))].copy()
nonfh = nonfh.merge(lpa[['eid', 'lpa']], on='eid', how='inner')
nonfh = nonfh.merge(sex[['eid', 'sex']], on='eid', how='left')
nonfh = nonfh.merge(demo[['eid', 'age']], on='eid', how='left')
nonfh = nonfh[nonfh['lpa'].notna()].copy()

print(f"FH with Lp(a): {len(fh_lpa)}, ASCVD: {fh_lpa['ascvd'].sum():.0f}")
print(f"Non-FH with Lp(a): {len(nonfh)}, ASCVD: {nonfh['ascvd'].sum():.0f}")

# =====================================================================
# GAP 1: FH-SPECIFIC Lp(a) RISK THRESHOLDS
# =====================================================================
print(f"\n{'='*70}")
print("GAP 1: FH-SPECIFIC Lp(a) THRESHOLDS FOR ASCVD")
print("Are general population thresholds valid in FH?")
print("=" * 70)

# AUC comparison
auc_fh = roc_auc_score(fh_lpa['ascvd'], fh_lpa['lpa'])
auc_nonfh = roc_auc_score(nonfh['ascvd'], nonfh['lpa'])
print(f"\n  Lp(a) AUC for ASCVD:")
print(f"    FH carriers:     {auc_fh:.3f}")
print(f"    Non-FH:          {auc_nonfh:.3f}")
print(f"    Difference:      {auc_fh - auc_nonfh:+.3f}")
print(f"    Interpretation:  Lp(a) is {'USELESS' if auc_fh < 0.55 else 'weak'} for ASCVD prediction in FH")

# Optimal thresholds by Youden's index
fpr_fh, tpr_fh, thresh_fh = roc_curve(fh_lpa['ascvd'], fh_lpa['lpa'])
youden_fh = tpr_fh - fpr_fh
optimal_idx_fh = np.argmax(youden_fh)
optimal_thresh_fh = thresh_fh[optimal_idx_fh]

fpr_nf, tpr_nf, thresh_nf = roc_curve(nonfh['ascvd'], nonfh['lpa'])
youden_nf = tpr_nf - fpr_nf
optimal_idx_nf = np.argmax(youden_nf)
optimal_thresh_nf = thresh_nf[optimal_idx_nf]

print(f"\n  Optimal Youden threshold:")
print(f"    FH:     {optimal_thresh_fh:.1f} nmol/L (sens={tpr_fh[optimal_idx_fh]:.2f}, spec={1-fpr_fh[optimal_idx_fh]:.2f})")
print(f"    Non-FH: {optimal_thresh_nf:.1f} nmol/L (sens={tpr_nf[optimal_idx_nf]:.2f}, spec={1-fpr_nf[optimal_idx_nf]:.2f})")

# Compare standard thresholds
print(f"\n  Standard threshold performance:")
print(f"  {'Threshold':>12s}  {'FH Sens':>8s}  {'FH Spec':>8s}  {'FH PPV':>7s}  {'NonFH Sens':>10s}  {'NonFH Spec':>10s}")
print("  " + "-" * 65)

for thresh in [30, 50, 75, 100, 125, 150]:
    for label, df_a in [('FH', fh_lpa), ('NonFH', nonfh)]:
        hi = df_a[df_a['lpa'] > thresh]
        lo = df_a[df_a['lpa'] <= thresh]
        tp = hi['ascvd'].sum()
        fp = len(hi) - tp
        fn = lo['ascvd'].sum()
        tn = len(lo) - fn
        sens = tp / max(tp + fn, 1)
        spec = tn / max(tn + fp, 1)
        ppv = tp / max(tp + fp, 1)
        if label == 'FH':
            print(f"  >{thresh:>3d} nmol/L  {sens:>8.2f}  {spec:>8.2f}  {ppv:>7.2f}", end='')
        else:
            print(f"  {sens:>10.2f}  {spec:>10.2f}")

# =====================================================================
# GAP 2: COMBINED STRUCTURAL SEVERITY + Lp(a) RISK MODEL
# =====================================================================
print(f"\n{'='*70}")
print("GAP 2: COMBINED SSS x Lp(a) RISK MODEL")
print("=" * 70)

# Use analysis-ready for SSS + Lp(a) combination
ar = pd.read_csv("D:/CALON_AF3_PROJECT/data/fh_cohort/calon_ukb_analysis_ready.csv")
ar['lpa_val'] = pd.to_numeric(ar['lpa'], errors='coerce')
ar['ascvd'] = pd.to_numeric(ar['ascvd_combined'], errors='coerce')

# We don't have SSS in analysis_ready, but we have domain + mechanism
# Use domain-level penetrance as SSS proxy
atlas = pd.read_csv(os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas_v3.csv"))
pen_map = atlas.groupby('domain')['domain_penetrance_by_60_pct'].first().to_dict()

# For UKB carriers: domain + Lp(a) compound risk
fh_compound = fh_lpa[fh_lpa['ldlr_domain'].notna()].copy()
fh_compound['penetrance'] = fh_compound['ldlr_domain'].map(pen_map)
fh_compound['high_pen'] = (fh_compound['penetrance'] > 20).astype(int)
fh_compound['high_lpa'] = (fh_compound['lpa'] > fh_compound['lpa'].median()).astype(int)

if len(fh_compound) >= 50:
    print(f"\n  Compound risk (UKB carriers, n={len(fh_compound)}):")
    print(f"  {'Group':>30s}  {'N':>6s}  {'ASCVD':>6s}  {'Rate%':>7s}")
    print("  " + "-" * 55)

    for pen_lab, pen_val in [('Low penetrance', 0), ('High penetrance', 1)]:
        for lpa_lab, lpa_val in [('Low Lp(a)', 0), ('High Lp(a)', 1)]:
            sub = fh_compound[(fh_compound['high_pen'] == pen_val) &
                              (fh_compound['high_lpa'] == lpa_val)]
            if len(sub) >= 10:
                rate = 100 * sub['ascvd'].mean()
                print(f"  {pen_lab + ' + ' + lpa_lab:>30s}  {len(sub):>6d}  "
                      f"{sub['ascvd'].sum():>6.0f}  {rate:>7.1f}")

    # Gradient
    ll = fh_compound[(fh_compound['high_pen'] == 0) & (fh_compound['high_lpa'] == 0)]['ascvd'].mean()
    hh = fh_compound[(fh_compound['high_pen'] == 1) & (fh_compound['high_lpa'] == 1)]['ascvd'].mean()
    if ll > 0:
        print(f"\n  Gradient: {hh/ll:.1f}-fold (High pen + High Lp(a) vs Low + Low)")

    # Interaction test
    X = fh_compound[['high_pen', 'high_lpa']].copy()
    X['interaction'] = X['high_pen'] * X['high_lpa']
    y = fh_compound['ascvd']
    valid = X.notna().all(axis=1) & y.notna()
    if valid.sum() >= 50:
        lr = LogisticRegression(max_iter=1000)
        lr.fit(X[valid].values, y[valid].values)
        print(f"\n  Logistic regression:")
        for feat, coef in zip(['High penetrance', 'High Lp(a)', 'Pen x Lp(a)'], lr.coef_[0]):
            print(f"    {feat:>20s}: OR={np.exp(coef):.3f}")

# =====================================================================
# GAP 3: PSEUDO-FH QUANTIFICATION
# =====================================================================
print(f"\n{'='*70}")
print("GAP 3: PSEUDO-FH — How many FH diagnoses are Lp(a)-driven?")
print("=" * 70)

ar_lpa = ar[ar['lpa_val'].notna()].copy()
ar_lpa['ldl_val'] = pd.to_numeric(ar_lpa['ldl'], errors='coerce')

if len(ar_lpa) > 0 and 'ldl_val' in ar_lpa.columns:
    # Lp(a)-corrected LDL
    ar_lpa['lpa_chol'] = (ar_lpa['lpa_val'] / 2.5) * 0.30 / 38.7
    ar_lpa['ldl_corrected'] = ar_lpa['ldl_val'] - ar_lpa['lpa_chol']

    # How many have LDL > 4.9 only because of Lp(a)?
    above_raw = ar_lpa['ldl_val'] > 4.9
    above_corr = ar_lpa['ldl_corrected'] > 4.9
    pseudo_fh = above_raw & ~above_corr

    print(f"\n  Wales FH cohort (n={len(ar_lpa)}):")
    print(f"    LDL >4.9 raw: {above_raw.sum()} ({100*above_raw.mean():.1f}%)")
    print(f"    LDL >4.9 corrected: {above_corr.sum()} ({100*above_corr.mean():.1f}%)")
    print(f"    'Pseudo-FH' (LDL >4.9 only due to Lp(a)): {pseudo_fh.sum()} ({100*pseudo_fh.mean():.1f}%)")

    # Lp(a) in pseudo-FH vs true high LDL
    if pseudo_fh.sum() >= 5:
        pseudo = ar_lpa[pseudo_fh]
        true_hi = ar_lpa[above_raw & above_corr]
        print(f"\n    Pseudo-FH group: mean Lp(a)={pseudo['lpa_val'].mean():.0f} nmol/L")
        print(f"    True high LDL: mean Lp(a)={true_hi['lpa_val'].mean():.0f} nmol/L")

    # By gene: are APOB carriers more likely to be "pseudo-FH"?
    print(f"\n    Pseudo-FH by gene:")
    for gene in ['LDLR', 'APOB']:
        sub = ar_lpa[ar_lpa['gene'] == gene]
        n_pseudo = (sub['ldl_val'] > 4.9).sum() - (sub['ldl_corrected'] > 4.9).sum()
        n_above = (sub['ldl_val'] > 4.9).sum()
        pct = 100 * n_pseudo / max(n_above, 1)
        print(f"      {gene}: {n_pseudo}/{n_above} pseudo-FH ({pct:.1f}%)")

# =====================================================================
# GAP 4: DOMAIN-SPECIFIC ASCVD WITH Lp(a) INTERACTION
# =====================================================================
print(f"\n{'='*70}")
print("GAP 4: DOMAIN x Lp(a) INTERACTION FOR ASCVD")
print("Does Lp(a) amplify ASCVD risk differently by domain?")
print("=" * 70)

fh_dom = fh_lpa[fh_lpa['ldlr_domain'].notna()].copy()
fh_dom['lpa_high'] = (fh_dom['lpa'] > fh_dom['lpa'].median()).astype(int)

dom_interaction = fh_dom.groupby(['ldlr_domain', 'lpa_high']).agg(
    n=('ascvd', 'count'),
    ascvd_rate=('ascvd', 'mean'),
).reset_index()

# Show domains with enough data
print(f"\n  {'Domain':>25s}  {'Low Lp(a)':>12s}  {'High Lp(a)':>12s}  {'Diff':>8s}")
print("  " + "-" * 65)

for dom in fh_dom['ldlr_domain'].unique():
    lo = dom_interaction[(dom_interaction['ldlr_domain'] == dom) & (dom_interaction['lpa_high'] == 0)]
    hi = dom_interaction[(dom_interaction['ldlr_domain'] == dom) & (dom_interaction['lpa_high'] == 1)]
    if len(lo) > 0 and len(hi) > 0 and lo.iloc[0]['n'] >= 10 and hi.iloc[0]['n'] >= 10:
        r_lo = 100 * lo.iloc[0]['ascvd_rate']
        r_hi = 100 * hi.iloc[0]['ascvd_rate']
        diff = r_hi - r_lo
        short = dom.replace('Ligand-binding ', 'LB-').replace('EGF-like ', 'EGF-')
        print(f"  {short:>25s}  {r_lo:>10.1f}%  {r_hi:>10.1f}%  {diff:>+7.1f}")

# =====================================================================
# FIGURE: 4-panel Nature-quality
# =====================================================================
print(f"\n{'='*70}")
print("Generating Nature-quality figure...")

fig, axes = plt.subplots(2, 2, figsize=(14, 11))
RED = '#b2182b'
BLUE = '#2166ac'
ORANGE = '#ef8a62'

# Panel A: ROC curves FH vs non-FH
ax = axes[0, 0]
ax.plot(fpr_fh, tpr_fh, color=RED, linewidth=2, label=f'FH (AUC={auc_fh:.3f})')
ax.plot(fpr_nf, tpr_nf, color=BLUE, linewidth=2, label=f'Non-FH (AUC={auc_nonfh:.3f})')
ax.plot([0, 1], [0, 1], '--', color='grey', linewidth=0.5)
ax.set_xlabel('1 - Specificity', fontsize=10)
ax.set_ylabel('Sensitivity', fontsize=10)
ax.set_title('A. Lp(a) for ASCVD: FH vs Non-FH\n(Lp(a) is useless in FH)', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: Threshold comparison (sensitivity at each cutoff)
ax = axes[0, 1]
thresholds = [25, 50, 75, 100, 125, 150]
sens_fh_list = []
sens_nf_list = []
for t in thresholds:
    tp_fh = fh_lpa[(fh_lpa['lpa'] > t) & (fh_lpa['ascvd'] == 1)].shape[0]
    fn_fh = fh_lpa[(fh_lpa['lpa'] <= t) & (fh_lpa['ascvd'] == 1)].shape[0]
    sens_fh_list.append(tp_fh / max(tp_fh + fn_fh, 1))

    tp_nf = nonfh[(nonfh['lpa'] > t) & (nonfh['ascvd'] == 1)].shape[0]
    fn_nf = nonfh[(nonfh['lpa'] <= t) & (nonfh['ascvd'] == 1)].shape[0]
    sens_nf_list.append(tp_nf / max(tp_nf + fn_nf, 1))

ax.plot(thresholds, sens_fh_list, 'o-', color=RED, linewidth=2, label='FH')
ax.plot(thresholds, sens_nf_list, 's-', color=BLUE, linewidth=2, label='Non-FH')
ax.set_xlabel('Lp(a) Threshold (nmol/L)', fontsize=10)
ax.set_ylabel('Sensitivity for ASCVD', fontsize=10)
ax.set_title('B. Sensitivity at Standard Thresholds', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: Compound risk (penetrance x Lp(a))
ax = axes[1, 0]
if len(fh_compound) >= 50:
    groups = []
    labels = []
    colours = []
    for pen_lab, pen_val, col in [('Low Pen', 0, BLUE), ('High Pen', 1, RED)]:
        for lpa_lab, lpa_val in [('Low Lp(a)', 0), ('High Lp(a)', 1)]:
            sub = fh_compound[(fh_compound['high_pen'] == pen_val) &
                              (fh_compound['high_lpa'] == lpa_val)]
            if len(sub) >= 10:
                rate = 100 * sub['ascvd'].mean()
                groups.append(rate)
                labels.append(f"{pen_lab}\n{lpa_lab}\n(n={len(sub)})")
                colours.append(col if lpa_val == 0 else ORANGE if pen_val == 0 else '#67001f')

    if groups:
        bars = ax.bar(range(len(groups)), groups, color=colours, edgecolor='black', linewidth=0.5)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_ylabel('ASCVD Rate (%)', fontsize=10)
        ax.set_title('C. Domain Penetrance x Lp(a) Compound Risk', fontsize=11, fontweight='bold')
        for bar, val in zip(bars, groups):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                    f'{val:.1f}%', ha='center', fontsize=9, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Domain x Lp(a) interaction heatmap
ax = axes[1, 1]
pivot = dom_interaction[dom_interaction['n'] >= 10].pivot_table(
    index='ldlr_domain', columns='lpa_high', values='ascvd_rate', aggfunc='first')
if len(pivot) >= 3:
    pivot.columns = ['Low Lp(a)', 'High Lp(a)']
    pivot = pivot.dropna()
    pivot.index = [d.replace('Ligand-binding ', 'LB-').replace('EGF-like ', 'EGF-')
                   for d in pivot.index]
    pivot_pct = pivot * 100

    im = ax.imshow(pivot_pct.values, cmap='RdYlBu_r', aspect='auto')
    ax.set_xticks(range(len(pivot_pct.columns)))
    ax.set_xticklabels(pivot_pct.columns, fontsize=10)
    ax.set_yticks(range(len(pivot_pct.index)))
    ax.set_yticklabels(pivot_pct.index, fontsize=8)
    ax.set_title('D. ASCVD Rate: Domain x Lp(a)', fontsize=11, fontweight='bold')

    for i in range(len(pivot_pct.index)):
        for j in range(len(pivot_pct.columns)):
            val = pivot_pct.values[i, j]
            if not np.isnan(val):
                ax.text(j, i, f'{val:.0f}%', ha='center', va='center', fontsize=8,
                        color='white' if val > 15 else 'black')

    plt.colorbar(im, ax=ax, shrink=0.7, label='ASCVD Rate (%)')

plt.tight_layout(pad=2.0)
fig_path = os.path.join(FIGURES, "Figure_Nature_Gaps.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {fig_path}")

print("\n" + "=" * 70)
print("COMPLETE")
print("=" * 70)

#!/usr/bin/env python3
"""
47_split_cohort_and_clinvar_mega_validation.py
===============================================
The two analyses that make or break a Nature paper:

PART A: SPLIT-COHORT VALIDATION
  - Derive SSS on Wales+Dragon3 (n=3,841)
  - Validate on UKB FH carriers (n=337)
  - Reports AUC, calibration, NRI, IDI

PART B: COMPREHENSIVE CLINVAR VALIDATION
  - 1,859 ClinVar missense variants matched to 16,340 FoldX saturation mutations
  - Pathogenic vs Benign ddG discrimination
  - ROC curve for pathogenicity prediction
  - Dose-response: ddG gradient across ClinVar categories

Outputs:
  alphafold/analysis/split_cohort_validation.csv
  alphafold/analysis/clinvar_mega_validation.csv
  alphafold/analysis/figures/Figure_SplitCohort_Validation.png
  alphafold/analysis/figures/Figure_ClinVar_MegaValidation.png
"""

import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve, brier_score_loss
from sklearn.calibration import calibration_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")
FIGURES = os.path.join(ANALYSIS, "figures")
os.makedirs(FIGURES, exist_ok=True)

# ============================================================================
# ============================================================================
#  PART A: SPLIT-COHORT VALIDATION
# ============================================================================
# ============================================================================
print("=" * 80)
print("PART A: SPLIT-COHORT VALIDATION")
print("Derive on Wales+Dragon3 -> Validate on UKB")
print("=" * 80)

# Load data
dev_val = pd.read_csv(os.path.join(ANALYSIS, "sss_v2_ukb_develop_wales_validate.csv"))
ukb_ext = pd.read_csv(os.path.join(ANALYSIS, "ukb_external_validation_sss.csv"))
val_data = pd.read_csv(os.path.join(ANALYSIS, "sss_validation_data.csv"))

print(f"\n  Development-Validation file: {len(dev_val)} patients")
print(f"    UKB: {(dev_val['cohort']=='UKB').sum()}")
print(f"    Wales: {(dev_val['cohort']=='Wales').sum()}")
print(f"  UKB external validation: {len(ukb_ext)} patients")
print(f"  SSS validation data: {len(val_data)} patients")

# Combine: Dragon3 + Wales = Development, UKB = Validation
# From val_data: Dragon3 + Wales patients with SSS + ASCVD
dragon3 = val_data[val_data['source'] == 'DRAGON3'].copy()
wales = val_data[val_data['source'] == 'WALES_ALL'].copy()
dev_set = pd.concat([dragon3, wales], ignore_index=True)

# UKB validation from dev_val
ukb_val = dev_val[dev_val['cohort'] == 'UKB'].copy()
ukb_val = ukb_val.rename(columns={'sss_v2': 'sss'})

print(f"\n  DEVELOPMENT SET:")
print(f"    Dragon3: {len(dragon3)} ({dragon3['ascvd'].sum()} events, {dragon3['ascvd'].mean()*100:.1f}%)")
print(f"    Wales:   {len(wales)} ({wales['ascvd'].sum()} events, {wales['ascvd'].mean()*100:.1f}%)")
print(f"    Total:   {len(dev_set)} ({dev_set['ascvd'].sum()} events, {dev_set['ascvd'].mean()*100:.1f}%)")

print(f"\n  VALIDATION SET:")
print(f"    UKB:     {len(ukb_val)} ({ukb_val['ascvd'].sum()} events, {ukb_val['ascvd'].mean()*100:.1f}%)")

# Also use the richer UKB external validation file
ukb_rich = ukb_ext.copy()
print(f"    UKB (rich): {len(ukb_rich)} ({ukb_rich['ascvd_combined'].sum()} events, "
      f"{ukb_rich['ascvd_combined'].mean()*100:.1f}%)")

# ── A1: SSS alone - development AUC ──
print("\n  --- A1: SSS Discrimination ---")

# Development: SSS -> ASCVD
dev_complete = dev_set.dropna(subset=['sss', 'ascvd'])
dev_y = dev_complete['ascvd'].values
dev_sss = dev_complete['sss'].values

if dev_y.sum() > 5:
    auc_dev = roc_auc_score(dev_y, dev_sss)
    print(f"  Development AUC (SSS alone): {auc_dev:.4f} (n={len(dev_y)}, events={dev_y.sum()})")
else:
    auc_dev = np.nan
    print(f"  Development: insufficient events ({dev_y.sum()})")

# Validation: SSS -> ASCVD (UKB)
val_complete = ukb_val.dropna(subset=['sss', 'ascvd'])
val_y = val_complete['ascvd'].values
val_sss = val_complete['sss'].values

if val_y.sum() > 5:
    auc_val = roc_auc_score(val_y, val_sss)
    print(f"  Validation AUC (SSS alone):  {auc_val:.4f} (n={len(val_y)}, events={val_y.sum()})")
else:
    auc_val = np.nan
    print(f"  Validation: insufficient events ({val_y.sum()})")

# Also use the rich UKB set
ukb_rich_complete = ukb_rich.dropna(subset=['sss', 'ascvd_combined'])
ukb_rich_y = ukb_rich_complete['ascvd_combined'].values
ukb_rich_sss = ukb_rich_complete['sss'].values

if ukb_rich_y.sum() > 5:
    auc_rich = roc_auc_score(ukb_rich_y, ukb_rich_sss)
    print(f"  Validation AUC (rich UKB):   {auc_rich:.4f} (n={len(ukb_rich_y)}, events={ukb_rich_y.sum()})")

# ── A2: Multivariable model - develop on Wales+Dragon3, validate on UKB ──
print("\n  --- A2: Multivariable Model ---")

# Development model: SSS + age + sex + ldl + on_statin
dev_mv = dev_complete.dropna(subset=['age', 'sex', 'ldl']).copy()
dev_mv_vars = ['sss', 'age', 'sex']
# Add ldl and statin if available
if 'ldl' in dev_mv.columns and dev_mv['ldl'].notna().sum() > 50:
    dev_mv_vars.append('ldl')

dev_mv_clean = dev_mv.dropna(subset=dev_mv_vars + ['ascvd'])
print(f"\n  Development (multivariable): n={len(dev_mv_clean)}, events={dev_mv_clean['ascvd'].sum()}")

if dev_mv_clean['ascvd'].sum() > 10:
    X_dev = dev_mv_clean[dev_mv_vars].values.astype(float)
    y_dev = dev_mv_clean['ascvd'].values

    model = LogisticRegression(max_iter=5000, penalty=None, solver='lbfgs')
    model.fit(X_dev, y_dev)

    y_dev_pred = model.predict_proba(X_dev)[:, 1]
    auc_dev_mv = roc_auc_score(y_dev, y_dev_pred)
    print(f"  Development AUC (SSS+age+sex+LDL): {auc_dev_mv:.4f}")

    # Print coefficients
    print(f"\n  Model coefficients (development):")
    print(f"    Intercept: {model.intercept_[0]:.4f}")
    for v, c in zip(dev_mv_vars, model.coef_[0]):
        print(f"    {v:15s}: beta={c:7.4f}  OR={np.exp(c):.3f}")

    # Validate on UKB
    # Map UKB columns
    ukb_for_val = ukb_rich.copy()
    ukb_for_val['ascvd'] = ukb_for_val['ascvd_combined']

    val_mv_vars = list(dev_mv_vars)  # same variables as development

    val_mv_clean = ukb_for_val.dropna(subset=val_mv_vars + ['ascvd']).copy()
    print(f"\n  Validation (UKB): n={len(val_mv_clean)}, events={val_mv_clean['ascvd'].sum()}")

    if val_mv_clean['ascvd'].sum() > 5:
        X_val = val_mv_clean[val_mv_vars].values.astype(float)
        y_val = val_mv_clean['ascvd'].values

        y_val_pred = model.predict_proba(X_val)[:, 1]
        auc_val_mv = roc_auc_score(y_val, y_val_pred)
        print(f"  Validation AUC (SSS+age+sex+LDL): {auc_val_mv:.4f}")

        # Brier score
        brier_dev = brier_score_loss(y_dev, y_dev_pred)
        brier_val = brier_score_loss(y_val, y_val_pred)
        print(f"\n  Brier score (development): {brier_dev:.4f}")
        print(f"  Brier score (validation):  {brier_val:.4f}")

        # Model without SSS (age + sex only)
        model_nosss = LogisticRegression(max_iter=5000, penalty=None, solver='lbfgs')
        nosss_vars = [v for v in dev_mv_vars if v != 'sss']
        model_nosss.fit(dev_mv_clean[nosss_vars].values.astype(float), y_dev)

        y_dev_nosss = model_nosss.predict_proba(dev_mv_clean[nosss_vars].values.astype(float))[:, 1]
        y_val_nosss = model_nosss.predict_proba(val_mv_clean[nosss_vars].values.astype(float))[:, 1]

        auc_dev_nosss = roc_auc_score(y_dev, y_dev_nosss)
        auc_val_nosss = roc_auc_score(y_val, y_val_nosss)

        print(f"\n  --- Model comparison ---")
        print(f"  {'Model':>30s}  {'Dev AUC':>8s}  {'Val AUC':>8s}  {'Delta':>7s}")
        print(f"  {'Age + Sex':>30s}  {auc_dev_nosss:>8.4f}  {auc_val_nosss:>8.4f}  {'ref':>7s}")
        print(f"  {'Age + Sex + SSS':>30s}  {auc_dev_mv:>8.4f}  {auc_val_mv:>8.4f}  "
              f"{auc_val_mv - auc_val_nosss:>+7.4f}")

        # NRI (Net Reclassification Improvement)
        # Binary NRI with 10% threshold
        threshold = 0.10
        # Without SSS
        pred_nosss_pos = y_val_nosss >= threshold
        pred_sss_pos = y_val_pred >= threshold

        # Events
        events = y_val == 1
        nonevents = y_val == 0

        # Event NRI: proportion reclassified up minus down (among events)
        up_events = (pred_sss_pos & ~pred_nosss_pos & events).sum()
        down_events = (~pred_sss_pos & pred_nosss_pos & events).sum()
        event_nri = (up_events - down_events) / events.sum() if events.sum() > 0 else 0

        # Non-event NRI: proportion reclassified down minus up (among non-events)
        up_nonevents = (pred_sss_pos & ~pred_nosss_pos & nonevents).sum()
        down_nonevents = (~pred_sss_pos & pred_nosss_pos & nonevents).sum()
        nonevent_nri = (down_nonevents - up_nonevents) / nonevents.sum() if nonevents.sum() > 0 else 0

        nri = event_nri + nonevent_nri

        print(f"\n  --- NRI (threshold {threshold:.0%}) ---")
        print(f"  Event NRI:     {event_nri:+.4f}")
        print(f"  Non-event NRI: {nonevent_nri:+.4f}")
        print(f"  Total NRI:     {nri:+.4f}")

        # IDI (Integrated Discrimination Improvement)
        idi_events = y_val_pred[events].mean() - y_val_nosss[events].mean()
        idi_nonevents = y_val_pred[nonevents].mean() - y_val_nosss[nonevents].mean()
        idi = idi_events - idi_nonevents

        print(f"\n  --- IDI ---")
        print(f"  IDI events:     {idi_events:+.6f}")
        print(f"  IDI non-events: {idi_nonevents:+.6f}")
        print(f"  IDI total:      {idi:+.6f}")

# ── A3: SSS tertile analysis in validation set ──
print("\n  --- A3: SSS Tertile Analysis (Validation Set) ---")

if len(val_complete) > 30 and val_y.sum() > 3:
    val_complete_copy = val_complete.copy()
    val_complete_copy['sss_tertile'] = pd.qcut(
        val_complete_copy['sss'], 3, labels=['T1 (Low)', 'T2 (Mid)', 'T3 (High)']
    )

    print(f"\n  {'Tertile':>12s}  {'N':>5s}  {'Events':>7s}  {'ASCVD%':>7s}  {'Mean SSS':>9s}")
    for t in ['T1 (Low)', 'T2 (Mid)', 'T3 (High)']:
        sub = val_complete_copy[val_complete_copy['sss_tertile'] == t]
        print(f"  {t:>12s}  {len(sub):>5d}  {sub['ascvd'].sum():>7.0f}  "
              f"{sub['ascvd'].mean()*100:>6.1f}%  {sub['sss'].mean():>9.3f}")

    # Trend test
    val_complete_copy['sss_t_num'] = val_complete_copy['sss_tertile'].map(
        {'T1 (Low)': 1, 'T2 (Mid)': 2, 'T3 (High)': 3})
    r_trend, p_trend = stats.spearmanr(val_complete_copy['sss_t_num'], val_complete_copy['ascvd'])
    print(f"  Trend test: rho={r_trend:.3f}, P={p_trend:.4f}")

# ── A4: Domain-level validation ──
print("\n  --- A4: Domain-Level ASCVD in Validation Set ---")

if 'domain' in ukb_val.columns:
    domain_val = ukb_val.groupby('domain').agg(
        n=('ascvd', 'count'),
        events=('ascvd', 'sum'),
        rate=('ascvd', 'mean'),
        mean_sss=('sss', 'mean'),
    ).reset_index()
    domain_val = domain_val[domain_val['n'] >= 3].sort_values('rate', ascending=False)

    print(f"\n  {'Domain':>30s}  {'N':>5s}  {'Events':>7s}  {'ASCVD%':>7s}  {'Mean SSS':>9s}")
    for _, r in domain_val.iterrows():
        print(f"  {r['domain']:>30s}  {int(r['n']):>5d}  {int(r['events']):>7d}  "
              f"{r['rate']*100:>6.1f}%  {r['mean_sss']:>9.3f}")

# ── A5: Calibration ──
print("\n  --- A5: Calibration ---")
if 'y_val_pred' in dir() and len(y_val) > 30:
    # Calibration by decile
    try:
        prob_true, prob_pred = calibration_curve(y_val, y_val_pred, n_bins=5, strategy='quantile')
        print(f"  Calibration (5 bins):")
        print(f"  {'Predicted':>10s}  {'Observed':>10s}")
        for pt, pp in zip(prob_true, prob_pred):
            print(f"  {pp:>10.4f}  {pt:>10.4f}")

        # Hosmer-Lemeshow-like: observed/expected ratio
        oe_ratio = y_val.mean() / y_val_pred.mean()
        print(f"\n  O/E ratio: {oe_ratio:.3f} (1.0 = perfect calibration)")
    except Exception as e:
        print(f"  Calibration error: {e}")

# ── Figure: Split-Cohort Validation ──
print("\n  Generating split-cohort figure...")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: ROC curves (development vs validation)
ax = axes[0, 0]
if 'y_dev_pred' in dir() and not np.isnan(auc_dev_mv):
    fpr_d, tpr_d, _ = roc_curve(y_dev, y_dev_pred)
    ax.plot(fpr_d, tpr_d, '-', color='#2166ac', linewidth=2,
            label=f'Development (AUC={auc_dev_mv:.3f}, n={len(y_dev)})')
if 'y_val_pred' in dir() and not np.isnan(auc_val_mv):
    fpr_v, tpr_v, _ = roc_curve(y_val, y_val_pred)
    ax.plot(fpr_v, tpr_v, '-', color='#b2182b', linewidth=2,
            label=f'Validation (AUC={auc_val_mv:.3f}, n={len(y_val)})')
ax.plot([0, 1], [0, 1], 'k--', linewidth=0.5)
ax.set_xlabel('1 - Specificity', fontsize=11, fontfamily='Arial')
ax.set_ylabel('Sensitivity', fontsize=11, fontfamily='Arial')
ax.set_title('A. ROC: Development vs External Validation', fontsize=12,
             fontweight='bold', fontfamily='Arial')
ax.legend(fontsize=9, loc='lower right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: Calibration plot
ax = axes[0, 1]
if 'prob_true' in dir():
    ax.plot(prob_pred, prob_true, 'o-', color='#b2182b', markersize=8, linewidth=2,
            label='Observed')
    ax.plot([0, max(prob_pred.max(), prob_true.max())],
            [0, max(prob_pred.max(), prob_true.max())], 'k--', linewidth=0.5,
            label='Perfect calibration')
    ax.set_xlabel('Predicted Probability', fontsize=11, fontfamily='Arial')
    ax.set_ylabel('Observed Probability', fontsize=11, fontfamily='Arial')
    ax.set_title('B. Calibration (Validation Set)', fontsize=12,
                 fontweight='bold', fontfamily='Arial')
    ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: SSS tertile bar chart (validation)
ax = axes[1, 0]
if 'val_complete_copy' in dir() and val_y.sum() > 3:
    tertiles = ['T1 (Low)', 'T2 (Mid)', 'T3 (High)']
    rates = []
    ns = []
    for t in tertiles:
        sub = val_complete_copy[val_complete_copy['sss_tertile'] == t]
        rates.append(sub['ascvd'].mean() * 100)
        ns.append(len(sub))
    bars = ax.bar(range(3), rates, color=['#2166ac', '#f4a582', '#b2182b'], edgecolor='black')
    ax.set_xticks(range(3))
    ax.set_xticklabels([f'{t}\n(n={n})' for t, n in zip(tertiles, ns)], fontsize=9)
    ax.set_ylabel('ASCVD Rate (%)', fontsize=11, fontfamily='Arial')
    ax.set_title('C. ASCVD by SSS Tertile (Validation)', fontsize=12,
                 fontweight='bold', fontfamily='Arial')
    # Add percentage labels on bars
    for bar, rate in zip(bars, rates):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.5,
                f'{rate:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Forest plot - SSS effect across cohorts
ax = axes[1, 1]
cohort_results = []
for label, y_true, y_score in [
    ('Dragon3', dragon3.dropna(subset=['sss','ascvd'])['ascvd'].values,
     dragon3.dropna(subset=['sss','ascvd'])['sss'].values),
    ('Wales', wales.dropna(subset=['sss','ascvd'])['ascvd'].values,
     wales.dropna(subset=['sss','ascvd'])['sss'].values),
    ('UKB', val_y, val_sss),
]:
    if y_true.sum() > 3:
        auc = roc_auc_score(y_true, y_score)
        # Bootstrap CI
        np.random.seed(42)
        boot_aucs = []
        for _ in range(1000):
            idx = np.random.choice(len(y_true), len(y_true), replace=True)
            try:
                boot_aucs.append(roc_auc_score(y_true[idx], y_score[idx]))
            except:
                pass
        ci_lo = np.percentile(boot_aucs, 2.5)
        ci_hi = np.percentile(boot_aucs, 97.5)
        cohort_results.append({
            'cohort': label, 'auc': auc, 'ci_lo': ci_lo, 'ci_hi': ci_hi,
            'n': len(y_true), 'events': y_true.sum()
        })

if cohort_results:
    y_pos = list(range(len(cohort_results)))
    for i, r in enumerate(cohort_results):
        ax.errorbar(r['auc'], i, xerr=[[r['auc']-r['ci_lo']], [r['ci_hi']-r['auc']]],
                   fmt='o', color='#b2182b', markersize=8, capsize=5, linewidth=2)
        ax.text(r['ci_hi'] + 0.01, i,
                f"{r['auc']:.3f} ({r['ci_lo']:.3f}-{r['ci_hi']:.3f})\nn={r['n']}, events={int(r['events'])}",
                va='center', fontsize=8, fontfamily='Arial')

    ax.set_yticks(y_pos)
    ax.set_yticklabels([r['cohort'] for r in cohort_results], fontsize=10)
    ax.axvline(x=0.5, color='gray', linestyle='--', linewidth=0.5)
    ax.set_xlabel('AUC (95% CI)', fontsize=11, fontfamily='Arial')
    ax.set_title('D. SSS Discrimination Across Cohorts', fontsize=12,
                 fontweight='bold', fontfamily='Arial')
    ax.set_xlim(0.3, 1.0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig_path = os.path.join(FIGURES, "Figure_SplitCohort_Validation.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Figure saved: {fig_path}")

# Save split-cohort results
split_results = []
for r in cohort_results:
    split_results.append(r)
pd.DataFrame(split_results).to_csv(
    os.path.join(ANALYSIS, "split_cohort_validation.csv"), index=False)


# ============================================================================
# ============================================================================
#  PART B: COMPREHENSIVE CLINVAR VALIDATION
# ============================================================================
# ============================================================================
print("\n\n" + "=" * 80)
print("PART B: COMPREHENSIVE CLINVAR VALIDATION")
print(f"1,859 ClinVar missense variants vs 16,340 FoldX saturation mutations")
print("=" * 80)

# Load data
clinvar = pd.read_csv(os.path.join(ANALYSIS, "clinvar_ldlr_full_missense.csv"))
sat = pd.read_csv(os.path.join(ANALYSIS, "foldx_saturation_mutagenesis.csv"))
atlas = pd.read_csv(os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas_v3.csv"))

print(f"\n  ClinVar missense: {len(clinvar)}")
print(f"  Saturation mutations: {len(sat)}")

# Standardise ClinVar classifications
def classify_cv(sig):
    if pd.isna(sig):
        return 'Unknown'
    s = str(sig).lower()
    if 'pathogenic/likely' in s:
        return 'Pathogenic/LP'
    elif 'pathogenic' in s and 'likely' not in s and 'benign' not in s and 'conflicting' not in s:
        return 'Pathogenic'
    elif 'likely pathogenic' in s:
        return 'Likely Pathogenic'
    elif 'benign/likely' in s:
        return 'Benign/LB'
    elif 'benign' in s and 'likely' not in s and 'pathogenic' not in s and 'conflicting' not in s:
        return 'Benign'
    elif 'likely benign' in s:
        return 'Likely Benign'
    elif 'conflicting' in s:
        return 'Conflicting'
    elif 'uncertain' in s:
        return 'VUS'
    elif 'not provided' in s:
        return 'Not provided'
    else:
        return 'Unknown'

clinvar['cv_class'] = clinvar['clinical_significance'].apply(classify_cv)

print(f"\n  ClinVar classifications:")
for cls in ['Pathogenic', 'Pathogenic/LP', 'Likely Pathogenic', 'VUS',
            'Conflicting', 'Likely Benign', 'Benign/LB', 'Benign', 'Not provided', 'Unknown']:
    n = (clinvar['cv_class'] == cls).sum()
    if n > 0:
        print(f"    {cls:25s}: {n:>5d}")

# ── B1: EXACT VARIANT MATCH ──
print("\n  --- B1: Exact Variant Match (ClinVar -> Saturation) ---")

sat_valid = sat.dropna(subset=['ddG']).copy()
sat_valid['position'] = sat_valid['position'].astype(int)

clinvar['position'] = pd.to_numeric(clinvar['position'], errors='coerce')
clinvar = clinvar.dropna(subset=['position'])
clinvar['position'] = clinvar['position'].astype(int)

# Filter to LDLR range (1-860)
clinvar_ldlr = clinvar[(clinvar['position'] >= 1) & (clinvar['position'] <= 860)].copy()
print(f"  ClinVar in LDLR range (1-860): {len(clinvar_ldlr)}")

# Merge by position + wt_aa + mut_aa
merged = clinvar_ldlr.merge(
    sat_valid[['position', 'wt_aa', 'mut_aa', 'ddG', 'effect_class']],
    on=['position', 'wt_aa', 'mut_aa'],
    how='left'
)

matched = merged[merged['ddG'].notna()].copy()
print(f"  ClinVar matched to exact FoldX ddG: {len(matched)}/{len(clinvar_ldlr)} "
      f"({len(matched)/len(clinvar_ldlr)*100:.1f}%)")

# Also add position-level stats from atlas
matched = matched.merge(
    atlas[['position', 'wildtype_aa', 'domain', 'sat_max_ddG', 'sat_mean_ddG',
           'mutational_sensitivity', 'risk_classification']],
    on='position', how='left'
)

print(f"\n  Exact ddG by ClinVar class:")
print(f"  {'Class':>25s}  {'N':>5s}  {'Median ddG':>10s}  {'Mean ddG':>9s}  {'% >2':>6s}  {'% >4':>6s}")
for cls in ['Pathogenic', 'Pathogenic/LP', 'Likely Pathogenic', 'VUS',
            'Conflicting', 'Likely Benign', 'Benign/LB', 'Benign']:
    sub = matched[matched['cv_class'] == cls]
    if len(sub) >= 3:
        pct_gt2 = (sub['ddG'] > 2).mean() * 100
        pct_gt4 = (sub['ddG'] > 4).mean() * 100
        print(f"  {cls:>25s}  {len(sub):>5d}  {sub['ddG'].median():>10.2f}  "
              f"{sub['ddG'].mean():>9.2f}  {pct_gt2:>5.1f}%  {pct_gt4:>5.1f}%")

# ── B2: Binary discrimination ──
print("\n  --- B2: Pathogenic vs Benign Discrimination ---")

# Combine pathogenic categories
path_mask = matched['cv_class'].isin(['Pathogenic', 'Pathogenic/LP', 'Likely Pathogenic'])
ben_mask = matched['cv_class'].isin(['Benign', 'Benign/LB', 'Likely Benign'])

path_ddg = matched[path_mask]['ddG']
ben_ddg = matched[ben_mask]['ddG']

print(f"  Pathogenic+LP (n={len(path_ddg)}): median={path_ddg.median():.2f}, mean={path_ddg.mean():.2f}")
print(f"  Benign+LB (n={len(ben_ddg)}):      median={ben_ddg.median():.2f}, mean={ben_ddg.mean():.2f}")

if len(path_ddg) > 5 and len(ben_ddg) > 3:
    u, p_mw = stats.mannwhitneyu(path_ddg, ben_ddg, alternative='greater')
    auc_bin = roc_auc_score(
        [1]*len(path_ddg) + [0]*len(ben_ddg),
        list(path_ddg) + list(ben_ddg)
    )
    print(f"  Mann-Whitney P = {p_mw:.2e}")
    print(f"  AUC (ddG discriminates pathogenic): {auc_bin:.3f}")

    # Effect size (Cohen's d)
    pooled_std = np.sqrt((path_ddg.std()**2 * (len(path_ddg)-1) +
                          ben_ddg.std()**2 * (len(ben_ddg)-1)) /
                         (len(path_ddg) + len(ben_ddg) - 2))
    cohens_d = (path_ddg.mean() - ben_ddg.mean()) / pooled_std if pooled_std > 0 else 0
    print(f"  Cohen's d = {cohens_d:.3f}")
else:
    auc_bin = np.nan
    print(f"  Insufficient benign variants for comparison (n={len(ben_ddg)})")

# ── B3: ROC curve for pathogenicity ──
print("\n  --- B3: ddG Thresholds for Pathogenicity ---")

if len(path_ddg) > 5 and len(ben_ddg) > 3:
    all_binary = pd.concat([
        pd.DataFrame({'ddG': path_ddg, 'pathogenic': 1}),
        pd.DataFrame({'ddG': ben_ddg, 'pathogenic': 0}),
    ])

    print(f"\n  {'Threshold':>10s}  {'Sens':>6s}  {'Spec':>6s}  {'PPV':>6s}  {'NPV':>6s}  {'Youden':>7s}")
    best_youden = -1
    best_thresh = 0
    for thresh in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
        pred = all_binary['ddG'] >= thresh
        true = all_binary['pathogenic'] == 1
        tp = (pred & true).sum()
        fp = (pred & ~true).sum()
        fn = (~pred & true).sum()
        tn = (~pred & ~true).sum()
        sens = tp/(tp+fn) if (tp+fn) > 0 else 0
        spec = tn/(tn+fp) if (tn+fp) > 0 else 0
        ppv = tp/(tp+fp) if (tp+fp) > 0 else 0
        npv = tn/(tn+fn) if (tn+fn) > 0 else 0
        youden = sens + spec - 1
        if youden > best_youden:
            best_youden = youden
            best_thresh = thresh
        print(f"  {thresh:>10.1f}  {sens:>5.1%}  {spec:>5.1%}  {ppv:>5.1%}  {npv:>5.1%}  {youden:>7.3f}")

    print(f"\n  Optimal threshold (Youden): ddG >= {best_thresh:.1f} kcal/mol")

# ── B4: Position-level analysis ──
print("\n  --- B4: Position-Level Analysis ---")

# Positions with ONLY pathogenic variants vs ONLY benign
path_positions = set(clinvar_ldlr[clinvar_ldlr['cv_class'].isin(
    ['Pathogenic', 'Pathogenic/LP', 'Likely Pathogenic'])]['position'].unique())
ben_positions = set(clinvar_ldlr[clinvar_ldlr['cv_class'].isin(
    ['Benign', 'Benign/LB', 'Likely Benign'])]['position'].unique())

path_only_pos = path_positions - ben_positions
ben_only_pos = ben_positions - path_positions

print(f"  Positions with pathogenic variants only: {len(path_only_pos)}")
print(f"  Positions with benign variants only:     {len(ben_only_pos)}")
print(f"  Positions with both:                     {len(path_positions & ben_positions)}")

# Get atlas saturation stats for these positions
atlas_path_pos = atlas[atlas['position'].isin(path_only_pos)]
atlas_ben_pos = atlas[atlas['position'].isin(ben_only_pos)]

for metric, label in [('sat_max_ddG', 'Position Max ddG'),
                       ('sat_mean_ddG', 'Position Mean ddG'),
                       ('mutational_sensitivity', 'Mutational Sensitivity')]:
    p_vals = atlas_path_pos[metric].dropna()
    b_vals = atlas_ben_pos[metric].dropna()
    if len(p_vals) > 3 and len(b_vals) > 3:
        u, p = stats.mannwhitneyu(p_vals, b_vals, alternative='greater')
        try:
            auc_pos = roc_auc_score(
                [1]*len(p_vals) + [0]*len(b_vals),
                list(p_vals) + list(b_vals)
            )
        except:
            auc_pos = np.nan
        print(f"\n  {label}:")
        print(f"    Pathogenic (n={len(p_vals)}): median={p_vals.median():.2f}, mean={p_vals.mean():.2f}")
        print(f"    Benign (n={len(b_vals)}):     median={b_vals.median():.2f}, mean={b_vals.mean():.2f}")
        print(f"    P = {p:.2e}, AUC = {auc_pos:.3f}")

# ── B5: Domain enrichment ──
print("\n  --- B5: Domain Enrichment of Pathogenic Variants ---")

matched_with_domain = matched[matched['domain'].notna()].copy()
domain_cv = matched_with_domain.groupby('domain').agg(
    n_total=('ddG', 'count'),
    n_pathogenic=('cv_class', lambda x: x.isin(['Pathogenic', 'Pathogenic/LP', 'Likely Pathogenic']).sum()),
    n_benign=('cv_class', lambda x: x.isin(['Benign', 'Benign/LB', 'Likely Benign']).sum()),
    n_vus=('cv_class', lambda x: (x == 'VUS').sum()),
    median_ddg=('ddG', 'median'),
).reset_index()
domain_cv['pct_pathogenic'] = domain_cv['n_pathogenic'] / domain_cv['n_total'] * 100
domain_cv = domain_cv[domain_cv['n_total'] >= 5].sort_values('pct_pathogenic', ascending=False)

print(f"\n  {'Domain':>25s}  {'Total':>5s}  {'Path':>5s}  {'Ben':>4s}  {'VUS':>4s}  {'%Path':>6s}  {'Med ddG':>8s}")
for _, r in domain_cv.iterrows():
    print(f"  {r['domain']:>25s}  {int(r['n_total']):>5d}  {int(r['n_pathogenic']):>5d}  "
          f"{int(r['n_benign']):>4d}  {int(r['n_vus']):>4d}  {r['pct_pathogenic']:>5.1f}%  "
          f"{r['median_ddg']:>8.2f}")

# ── B6: Review status stratification ──
print("\n  --- B6: Expert-Reviewed vs Single Submitter ---")

matched['expert'] = matched['review_status'].str.contains('expert|multiple|practice', case=False, na=False)
matched['multi_sub'] = matched['review_status'].str.contains('multiple', case=False, na=False)

for label, mask in [('Expert/Multiple submitters', matched['multi_sub']),
                     ('Single submitter', ~matched['multi_sub'])]:
    sub = matched[mask]
    path_sub = sub[sub['cv_class'].isin(['Pathogenic', 'Pathogenic/LP', 'Likely Pathogenic'])]['ddG']
    ben_sub = sub[sub['cv_class'].isin(['Benign', 'Benign/LB', 'Likely Benign'])]['ddG']
    if len(path_sub) > 3:
        print(f"\n  {label}:")
        print(f"    Pathogenic: n={len(path_sub)}, median ddG={path_sub.median():.2f}")
        if len(ben_sub) > 0:
            print(f"    Benign:     n={len(ben_sub)}, median ddG={ben_sub.median():.2f}")

# ── Figure: ClinVar Mega-Validation ──
print("\n  Generating ClinVar mega-validation figure...")

fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# Panel A: Violin/Box by ClinVar class
ax = axes[0, 0]
plot_classes = ['Pathogenic', 'Pathogenic/LP', 'Likely Pathogenic', 'VUS',
                'Conflicting', 'Likely Benign', 'Benign']
plot_data = []
plot_labels = []
plot_colours = ['#67001f', '#b2182b', '#d6604d', '#999999', '#fddbc7', '#4393c3', '#053061']
actual_colours = []

for cls, col in zip(plot_classes, plot_colours):
    vals = matched[matched['cv_class'] == cls]['ddG']
    if len(vals) >= 3:
        plot_data.append(vals.values)
        plot_labels.append(f"{cls}\n(n={len(vals)})")
        actual_colours.append(col)

bp = ax.boxplot(plot_data, labels=plot_labels, patch_artist=True,
                medianprops=dict(color='black', linewidth=2),
                showfliers=True, flierprops=dict(markersize=3))
for patch, col in zip(bp['boxes'], actual_colours):
    patch.set_facecolor(col)
    patch.set_alpha(0.7)
ax.axhline(y=2, color='red', linestyle='--', alpha=0.5, linewidth=1)
ax.set_ylabel('FoldX ddG (kcal/mol)', fontsize=11, fontfamily='Arial')
ax.set_title('A. Structural Damage by ClinVar Classification', fontsize=12,
             fontweight='bold', fontfamily='Arial')
ax.tick_params(axis='x', labelsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: ROC curve
ax = axes[0, 1]
if len(path_ddg) > 5 and len(ben_ddg) > 3:
    y_true_roc = np.array([1]*len(path_ddg) + [0]*len(ben_ddg))
    y_score_roc = np.concatenate([path_ddg.values, ben_ddg.values])
    fpr_cv, tpr_cv, thresh_cv = roc_curve(y_true_roc, y_score_roc)
    ax.plot(fpr_cv, tpr_cv, '-', color='#b2182b', linewidth=2,
            label=f'AUC = {auc_bin:.3f}\nPath n={len(path_ddg)}, Ben n={len(ben_ddg)}')
    ax.plot([0, 1], [0, 1], 'k--', linewidth=0.5)
    ax.set_xlabel('1 - Specificity', fontsize=11, fontfamily='Arial')
    ax.set_ylabel('Sensitivity', fontsize=11, fontfamily='Arial')
    ax.set_title('B. ROC: ddG Discriminates Pathogenicity', fontsize=12,
                 fontweight='bold', fontfamily='Arial')
    ax.legend(fontsize=10, loc='lower right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: Dose-response — ddG quintile -> % pathogenic
ax = axes[1, 0]
matched_classified = matched[matched['cv_class'].isin(
    ['Pathogenic', 'Pathogenic/LP', 'Likely Pathogenic', 'Benign', 'Benign/LB', 'Likely Benign', 'VUS']
)].copy()
if len(matched_classified) > 20:
    matched_classified['ddg_quintile'] = pd.qcut(
        matched_classified['ddG'], 5,
        labels=['Q1\n(most stable)', 'Q2', 'Q3', 'Q4', 'Q5\n(most destab.)'],
        duplicates='drop'
    )
    matched_classified['is_path'] = matched_classified['cv_class'].isin(
        ['Pathogenic', 'Pathogenic/LP', 'Likely Pathogenic']).astype(int)

    quintile_data = matched_classified.groupby('ddg_quintile').agg(
        n=('is_path', 'count'),
        n_path=('is_path', 'sum'),
        pct_path=('is_path', 'mean'),
    ).reset_index()

    colours_q = ['#2166ac', '#67a9cf', '#d1e5f0', '#ef8a62', '#b2182b']
    bars = ax.bar(range(len(quintile_data)), quintile_data['pct_path'] * 100,
                  color=colours_q[:len(quintile_data)], edgecolor='black')
    ax.set_xticks(range(len(quintile_data)))
    ax.set_xticklabels([f"{q}\n(n={n})" for q, n in
                        zip(quintile_data['ddg_quintile'], quintile_data['n'])], fontsize=8)
    ax.set_ylabel('% Pathogenic/Likely Pathogenic', fontsize=11, fontfamily='Arial')
    ax.set_title('C. Dose-Response: ddG -> Pathogenicity Rate', fontsize=12,
                 fontweight='bold', fontfamily='Arial')
    for bar, pct in zip(bars, quintile_data['pct_path']):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1,
                f'{pct*100:.0f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Position landscape - pathogenic vs benign positions
ax = axes[1, 1]
p_max = atlas_path_pos['sat_max_ddG'].dropna()
b_max = atlas_ben_pos['sat_max_ddG'].dropna()
if len(p_max) > 3 and len(b_max) > 3:
    data_bp = [p_max.values, b_max.values]
    labels_bp = [f"Pathogenic\npositions\n(n={len(p_max)})",
                 f"Benign\npositions\n(n={len(b_max)})"]
    bp3 = ax.boxplot(data_bp, labels=labels_bp, patch_artist=True,
                     medianprops=dict(color='black', linewidth=2))
    bp3['boxes'][0].set_facecolor('#b2182b')
    bp3['boxes'][0].set_alpha(0.7)
    bp3['boxes'][1].set_facecolor('#2166ac')
    bp3['boxes'][1].set_alpha(0.7)
    for i, d in enumerate(data_bp):
        x = np.random.normal(i+1, 0.04, size=len(d))
        ax.scatter(x, d, alpha=0.2, s=8, color='black')
    ax.axhline(y=2, color='red', linestyle='--', alpha=0.5)
    ax.set_ylabel('Position Max ddG (kcal/mol)', fontsize=11, fontfamily='Arial')
    ax.set_title('D. Saturation Landscape: ClinVar Positions', fontsize=12,
                 fontweight='bold', fontfamily='Arial')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig_path = os.path.join(FIGURES, "Figure_ClinVar_MegaValidation.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Figure saved: {fig_path}")

# ── Save everything ──
matched.to_csv(os.path.join(ANALYSIS, "clinvar_mega_validation.csv"), index=False)

# ============================================================================
# FINAL SUMMARY TABLE
# ============================================================================
print("\n" + "=" * 80)
print("FINAL SUMMARY: VALIDATION METRICS")
print("=" * 80)

print(f"""
  SPLIT-COHORT VALIDATION (Development: Wales+Dragon3, External: UKB)
  -------------------------------------------------------------------
  Development n = {len(dev_complete)}, events = {dev_y.sum()}
  Validation n  = {len(val_complete)}, events = {val_y.sum()}
""")

if 'auc_dev_mv' in dir():
    print(f"  SSS alone:        Dev AUC = {auc_dev:.4f}, Val AUC = {auc_val:.4f}")
    print(f"  SSS+Age+Sex+LDL:  Dev AUC = {auc_dev_mv:.4f}, Val AUC = {auc_val_mv:.4f}")
    print(f"  Age+Sex (no SSS): Dev AUC = {auc_dev_nosss:.4f}, Val AUC = {auc_val_nosss:.4f}")
    print(f"  SSS incremental:  +{auc_val_mv - auc_val_nosss:.4f}")
    if 'nri' in dir():
        print(f"  NRI = {nri:+.4f}, IDI = {idi:+.6f}")

print(f"""
  CLINVAR VALIDATION (1,859 missense variants)
  -------------------------------------------------------------------
  Exact variant match: {len(matched)}/{len(clinvar_ldlr)} ({len(matched)/len(clinvar_ldlr)*100:.1f}%)
  Pathogenic+LP: n={len(path_ddg)}, median ddG = {path_ddg.median():.2f}
  Benign+LB:     n={len(ben_ddg)}, median ddG = {ben_ddg.median():.2f}
""")

if not np.isnan(auc_bin):
    print(f"  ddG AUC (pathogenic vs benign): {auc_bin:.3f}")
    print(f"  Mann-Whitney P = {p_mw:.2e}")
    print(f"  Cohen's d = {cohens_d:.3f}")

print(f"""
  Position-level: {len(path_only_pos)} pathogenic-only, {len(ben_only_pos)} benign-only positions
  Position max ddG: Pathogenic {atlas_path_pos['sat_max_ddG'].dropna().median():.2f} vs """
      f"""Benign {atlas_ben_pos['sat_max_ddG'].dropna().median():.2f}
""")

print("=" * 80)
print("COMPLETE")
print("=" * 80)

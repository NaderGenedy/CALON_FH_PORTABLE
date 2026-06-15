#!/usr/bin/env python3
"""
43_DCA_and_TUDOR_calculator.py
===============================
1. Decision Curve Analysis comparing models A/B/C
2. TUDOR clinical calculator with extracted coefficients
3. Publication-ready DCA figure

Outputs:
  alphafold/analysis/dca_results.csv
  alphafold/analysis/tudor_calculator_coefficients.csv
  alphafold/analysis/figures/Figure_DCA_v2.png
"""

import pandas as pd
import numpy as np
import os
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")
FIGURES = os.path.join(ANALYSIS, "figures")
os.makedirs(FIGURES, exist_ok=True)

INPUT_FILE = os.path.join(ANALYSIS, "comprehensive_sss_analysis.csv")

# ── Load data ────────────────────────────────────────────────────────────────
print("=" * 80)
print("DECISION CURVE ANALYSIS & TUDOR CALCULATOR")
print("=" * 80)

df = pd.read_csv(INPUT_FILE)
print(f"\nLoaded: {len(df)} patients, {df['ascvd'].sum()} ASCVD events ({df['ascvd'].mean()*100:.1f}%)")

# ── Define models ────────────────────────────────────────────────────────────
# Model A: Traditional risk factors only
# Model B: Traditional + SSS + xanthomata (the AUC 0.85 model)
# Model C: SSS alone

model_vars = {
    'A_traditional': ['age', 'sex', 'ldl1', 'on_statin', 'dm', 'smoking', 'hypertension'],
    'B_full_sss': ['age', 'sex', 'ldl1', 'on_statin', 'dm', 'smoking', 'hypertension', 'xanthomata', 'sss'],
    'C_sss_only': ['sss'],
}

# Complete cases for Model B (most restrictive)
complete_vars = ['ascvd', 'age', 'sex', 'ldl1', 'on_statin', 'dm', 'smoking',
                 'hypertension', 'xanthomata', 'sss']
cc = df[complete_vars].dropna()
print(f"Complete cases: {len(cc)} ({cc['ascvd'].sum()} events)")

y = cc['ascvd'].values

# ── Fit models and extract coefficients ──────────────────────────────────────
print("\n[1/3] Fitting logistic regression models...")

fitted_models = {}
predictions = {}
coefficients_all = []

for name, vars_list in model_vars.items():
    X = cc[vars_list].values.astype(float)

    model = LogisticRegression(max_iter=5000, penalty=None, solver='lbfgs')
    model.fit(X, y)

    y_pred = model.predict_proba(X)[:, 1]
    auc = roc_auc_score(y, y_pred)

    fitted_models[name] = model
    predictions[name] = y_pred

    print(f"\n  {name}:")
    print(f"    AUC: {auc:.4f}")
    print(f"    Intercept: {model.intercept_[0]:.4f}")

    for v, c in zip(vars_list, model.coef_[0]):
        or_val = np.exp(c)
        print(f"    {v:20s}: beta={c:7.4f}  OR={or_val:6.3f}")
        coefficients_all.append({
            'model': name,
            'variable': v,
            'coefficient': round(c, 6),
            'odds_ratio': round(or_val, 4),
            'intercept': round(model.intercept_[0], 6),
        })

# Save coefficients
coef_df = pd.DataFrame(coefficients_all)
coef_file = os.path.join(ANALYSIS, "tudor_calculator_coefficients.csv")
coef_df.to_csv(coef_file, index=False)
print(f"\n  Coefficients saved: {coef_file}")

# ── Decision Curve Analysis ──────────────────────────────────────────────────
print("\n[2/3] Running Decision Curve Analysis...")

def decision_curve(y_true, y_pred, thresholds):
    """Calculate net benefit at each threshold."""
    n = len(y_true)
    net_benefits = []
    for t in thresholds:
        # Predicted positive
        pred_pos = y_pred >= t
        tp = np.sum((pred_pos) & (y_true == 1))
        fp = np.sum((pred_pos) & (y_true == 0))

        # Net benefit
        nb = (tp / n) - (fp / n) * (t / (1 - t))
        net_benefits.append(nb)
    return np.array(net_benefits)

# Threshold range: 1% to 40% (clinically relevant for FH)
thresholds = np.arange(0.01, 0.41, 0.005)

# Treat all baseline
prevalence = y.mean()
nb_all = prevalence - (1 - prevalence) * (thresholds / (1 - thresholds))

# Treat none
nb_none = np.zeros_like(thresholds)

# Model net benefits
nb_models = {}
for name, y_pred in predictions.items():
    nb_models[name] = decision_curve(y, y_pred, thresholds)

# Save DCA results
dca_rows = []
for i, t in enumerate(thresholds):
    row = {
        'threshold': round(t, 4),
        'treat_all': round(nb_all[i], 6),
        'treat_none': 0,
    }
    for name in nb_models:
        row[f'nb_{name}'] = round(nb_models[name][i], 6)
    dca_rows.append(row)

dca_df = pd.DataFrame(dca_rows)
dca_file = os.path.join(ANALYSIS, "dca_results.csv")
dca_df.to_csv(dca_file, index=False)
print(f"  DCA results saved: {dca_file}")

# Key thresholds summary
print("\n  Net benefit at key thresholds:")
print(f"  {'Threshold':>10s}  {'Treat All':>10s}  {'Traditional':>12s}  {'Full+SSS':>10s}  {'SSS Only':>10s}")
for t_check in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]:
    idx = np.argmin(np.abs(thresholds - t_check))
    print(f"  {t_check:>10.0%}  {nb_all[idx]:>10.4f}  "
          f"{nb_models['A_traditional'][idx]:>12.4f}  "
          f"{nb_models['B_full_sss'][idx]:>10.4f}  "
          f"{nb_models['C_sss_only'][idx]:>10.4f}")

# ── Publication DCA Figure ───────────────────────────────────────────────────
print("\n[3/3] Generating DCA figure...")

fig, ax = plt.subplots(figsize=(8, 6))

# Style: Nature-compatible
ax.plot(thresholds * 100, nb_none, 'k-', linewidth=1, label='Treat none', alpha=0.5)
ax.plot(thresholds * 100, nb_all, 'k--', linewidth=1, label='Treat all', alpha=0.5)
ax.plot(thresholds * 100, nb_models['A_traditional'], '-',
        color='#2166ac', linewidth=2, label='Model A: Traditional (AUC={:.3f})'.format(
            roc_auc_score(y, predictions['A_traditional'])))
ax.plot(thresholds * 100, nb_models['B_full_sss'], '-',
        color='#b2182b', linewidth=2.5, label='Model B: Full+SSS (AUC={:.3f})'.format(
            roc_auc_score(y, predictions['B_full_sss'])))
ax.plot(thresholds * 100, nb_models['C_sss_only'], '-',
        color='#4dac26', linewidth=1.5, label='Model C: SSS alone (AUC={:.3f})'.format(
            roc_auc_score(y, predictions['C_sss_only'])))

ax.set_xlabel('Threshold Probability (%)', fontsize=12, fontfamily='Arial')
ax.set_ylabel('Net Benefit', fontsize=12, fontfamily='Arial')
ax.set_title('Decision Curve Analysis: ASCVD Prediction in FH', fontsize=13, fontfamily='Arial')
ax.legend(loc='upper right', fontsize=9, frameon=True, fancybox=False, edgecolor='black')
ax.set_xlim(1, 40)
ax.set_ylim(-0.05, max(prevalence * 1.1, 0.2))
ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig_path = os.path.join(FIGURES, "Figure_DCA_v2.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Figure saved: {fig_path}")

# ── TUDOR Calculator Formula ────────────────────────────────────────────────
print("\n" + "=" * 80)
print("TUDOR CLINICAL CALCULATOR")
print("=" * 80)

model_b = fitted_models['B_full_sss']
vars_b = model_vars['B_full_sss']

print("\nModel B (Full + SSS) — AUC {:.3f}".format(
    roc_auc_score(y, predictions['B_full_sss'])))
print("\nFormula:")
print("  P(ASCVD) = 1 / (1 + exp(-LP))")
print(f"\n  LP = {model_b.intercept_[0]:.4f}")
for v, c in zip(vars_b, model_b.coef_[0]):
    sign = "+" if c >= 0 else "-"
    print(f"       {sign} {abs(c):.4f} x {v}")

print("\n\nVariable encoding:")
print("  age:          years (continuous)")
print("  sex:          1=male, 0=female")
print("  ldl1:         LDL-C mmol/L at diagnosis (continuous)")
print("  on_statin:    1=yes, 0=no")
print("  dm:           1=diabetes, 0=no")
print("  smoking:      1=current/ex, 0=never")
print("  hypertension: 1=yes, 0=no")
print("  xanthomata:   1=present, 0=absent")
print("  sss:          Structural Severity Score (0-1, continuous)")

# Example predictions
print("\n\nExample predictions:")
examples = [
    {'desc': 'Low risk: 40yo F, LDL 5.0, no statin, no comorbidities, SSS=0.1',
     'vals': [40, 0, 5.0, 0, 0, 0, 0, 0, 0.1]},
    {'desc': 'Moderate: 55yo M, LDL 7.0, on statin, smoker, SSS=0.3',
     'vals': [55, 1, 7.0, 1, 0, 1, 0, 0, 0.3]},
    {'desc': 'High risk: 60yo M, LDL 9.0, statin, DM, HTN, xanthomata, SSS=0.6',
     'vals': [60, 1, 9.0, 1, 1, 1, 1, 1, 0.6]},
]

for ex in examples:
    X_ex = np.array(ex['vals']).reshape(1, -1)
    prob = model_b.predict_proba(X_ex)[0, 1]
    print(f"\n  {ex['desc']}")
    print(f"  -> P(ASCVD) = {prob:.1%}")

print("\n" + "=" * 80)
print("COMPLETE")
print("=" * 80)

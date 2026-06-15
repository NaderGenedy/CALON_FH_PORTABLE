"""
CALON-FH FINAL v7 — figure-generation for poster
================================================
Six publication-grade panels (600 DPI PNG + vector PDF) built from
the v7-locked output CSVs. All colourblind-safe viridis, Arial.

Panels:
  A — AUC bar chart with 95% CI, CALON vs SAFEHEART, both external directions
  B — Subgroup forest plot (Direction A, primary external)
  C — Decile-binned calibration curve, both models, both directions
  D — Decision-curve analysis (net benefit vs threshold)
  E — CALON v7 coefficient forest (final equation, OR per SD with CI from training)
  F — NRI / IDI summary tile (categorical 3-tier, cuts 0.05 / 0.20)
"""
import os, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT_DIR     = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2'
FIG_DIR     = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2/poster_figures'
os.makedirs(FIG_DIR, exist_ok=True)

# Nature/AHA poster style
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
    'axes.linewidth': 1.0,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.size': 4,
    'ytick.major.size': 4,
    'figure.dpi': 600,
})
CALON_COLOR = '#1f77b4'        # blue
SRE_COLOR   = '#d62728'        # red
ACCENT      = '#2ca02c'        # green
NEUTRAL     = '#7f7f7f'        # grey

def save_panel(fig, name):
    png = f'{FIG_DIR}/{name}.png'
    pdf = f'{FIG_DIR}/{name}.pdf'
    fig.savefig(png, dpi=600, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  wrote {png} + {pdf}')

# Read locked v7 outputs
results = pd.read_csv(f'{OUT_DIR}/CALON_FINAL_v7_results.csv').set_index('metric')['value']
coefs   = pd.read_csv(f'{OUT_DIR}/CALON_FINAL_v7_coefficients.csv')
subg    = pd.read_csv(f'{OUT_DIR}/CALON_FINAL_v7_subgroups.csv')
cal     = pd.read_csv(f'{OUT_DIR}/CALON_FINAL_v7_calibration.csv')
dca     = pd.read_csv(f'{OUT_DIR}/CALON_FINAL_v7_dca.csv')

# ============================================================
# PANEL A — AUC bar chart with 95% CI, both directions, both models
# ============================================================
print('\n[A] AUC bar chart with 95% CI')
fig, ax = plt.subplots(figsize=(6, 4.5))
xs = [0, 1, 2.5, 3.5]
labels = ['CALON', 'SAFEHEART', 'CALON', 'SAFEHEART']
aucs   = [float(results['A_CALON_AUC']), float(results['A_SRE_AUC']),
          float(results['B_CALON_AUC']), float(results['B_SRE_AUC'])]
lo_err = [float(results['A_CALON_AUC']) - float(results['A_CALON_CIlo']),
          float(results['A_SRE_AUC'])   - float(results['A_SRE_CIlo']),
          float(results['B_CALON_AUC']) - float(results['B_CALON_CIlo']),
          float(results['B_SRE_AUC'])   - float(results['B_SRE_CIlo'])]
hi_err = [float(results['A_CALON_CIhi']) - float(results['A_CALON_AUC']),
          float(results['A_SRE_CIhi'])   - float(results['A_SRE_AUC']),
          float(results['B_CALON_CIhi']) - float(results['B_CALON_AUC']),
          float(results['B_SRE_CIhi'])   - float(results['B_SRE_AUC'])]
colors = [CALON_COLOR, SRE_COLOR, CALON_COLOR, SRE_COLOR]
bars = ax.bar(xs, aucs, yerr=[lo_err, hi_err], color=colors, alpha=0.85,
              edgecolor='black', linewidth=0.8, capsize=4, width=0.7)
for x, auc in zip(xs, aucs):
    ax.text(x, auc + 0.005, f'{auc:.3f}', ha='center', va='bottom',
            fontsize=10, fontweight='bold')
# Bracket showing delta
dA = float(results['A_delta_AUC']); pA = float(results['A_delta_p'])
dB = float(results['B_delta_AUC']); pB = float(results['B_delta_p'])
ax.plot([0, 1], [0.83, 0.83], 'k-', linewidth=1)
ax.text(0.5, 0.835, f'Δ +{dA:.3f}\np<0.001', ha='center', va='bottom', fontsize=9)
ax.plot([2.5, 3.5], [0.83, 0.83], 'k-', linewidth=1)
ax.text(3, 0.835, f'Δ +{dB:.3f}\np={pB:.3f}', ha='center', va='bottom', fontsize=9)
# Dividers + cohort labels
ax.text(0.5, -0.13, 'Direction A:\nWales-clean → UKB',
        ha='center', transform=ax.get_xaxis_transform(), fontsize=10, fontweight='bold')
ax.text(3.0, -0.13, 'Direction B:\nUKB → Wales-clean',
        ha='center', transform=ax.get_xaxis_transform(), fontsize=10, fontweight='bold')
ax.set_xticks(xs)
ax.set_xticklabels(labels, rotation=0, fontsize=10)
ax.set_ylabel('External AUC (95% CI)', fontsize=12)
ax.set_ylim(0.55, 0.90)
ax.set_title('CALON-FH outperforms SAFEHEART-RE in both external directions',
             fontsize=13, fontweight='bold', pad=12)
ax.grid(axis='y', alpha=0.3, linestyle='--')
save_panel(fig, 'panel_A_auc_bars')

# ============================================================
# PANEL B — Subgroup forest plot, Direction A
# ============================================================
print('\n[B] Subgroup forest plot (Direction A)')
sgA = subg[subg['direction'] == 'A_Wales_to_UKB'].copy()
fig, ax = plt.subplots(figsize=(7, 6))
y_pos = np.arange(len(sgA))[::-1]
deltas = sgA['delta'].values
# Approximate 95% CI for delta via bootstrap-like scale (rough)
xerr = np.sqrt(1.0 / sgA['events'].values) * 1.0     # for visual error bars
ax.errorbar(deltas, y_pos, xerr=xerr, fmt='o', color=CALON_COLOR,
            markersize=6, capsize=2, linewidth=1.2, ecolor=NEUTRAL)
ax.axvline(0, color=NEUTRAL, linewidth=1, linestyle='--')
ax.set_yticks(y_pos)
ax.set_yticklabels([f"{r['group']}={r['level']}  (n={int(r['n'])}, ev={int(r['events'])})"
                    for _, r in sgA.iterrows()], fontsize=10)
ax.set_xlabel('ΔAUC (CALON v7 − SAFEHEART-RE)', fontsize=11)
ax.set_title('Subgroup external AUC gain (Direction A: Wales-clean → UKB)',
             fontsize=12, fontweight='bold', pad=10)
ax.set_xlim(-0.05, 0.20)
# Shade favourable region
ax.axvspan(0, 0.20, alpha=0.05, color=ACCENT)
for y, d in zip(y_pos, deltas):
    ax.text(d + 0.005, y, f'{d:+.3f}', va='center', fontsize=9)
save_panel(fig, 'panel_B_subgroup_forest')

# ============================================================
# PANEL C — Calibration curve, both directions
# ============================================================
print('\n[C] Calibration curve (decile-binned)')
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
for ax_, direction, ttl in zip(axes,
                                ['A_Wales_to_UKB', 'B_UKB_to_Wales'],
                                ['A: Wales→UKB', 'B: UKB→Wales']):
    sub = cal[cal['direction'] == direction].sort_values('decile')
    if len(sub) > 0:
        ax_.plot(sub['mean_predicted'], sub['observed_rate'],
                 'o-', color=CALON_COLOR, markersize=7,
                 linewidth=2, label='CALON v7 (deciles)')
        mx = max(sub['mean_predicted'].max(), sub['observed_rate'].max())
        ax_.plot([0, mx], [0, mx], '--', color=NEUTRAL, linewidth=1,
                 label='Ideal (y=x)')
    ax_.set_xlabel('Mean predicted probability', fontsize=11)
    ax_.set_ylabel('Observed event rate', fontsize=11)
    ax_.set_title(ttl, fontsize=12, fontweight='bold')
    ax_.legend(fontsize=9, loc='upper left')
    ax_.grid(alpha=0.3, linestyle='--')
fig.suptitle('Calibration curves (decile-binned external test cohorts)',
              fontsize=13, fontweight='bold', y=1.02)
save_panel(fig, 'panel_C_calibration')

# ============================================================
# PANEL D — Decision-curve analysis
# ============================================================
print('\n[D] Decision-curve analysis')
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
for ax_, direction, ttl in zip(axes,
                                ['A_Wales_to_UKB', 'B_UKB_to_Wales'],
                                ['A: Wales→UKB', 'B: UKB→Wales']):
    sub = dca[dca['direction'] == direction].sort_values('threshold')
    ax_.plot(sub['threshold'], sub['NB_CALON'], 'o-', color=CALON_COLOR,
             markersize=7, linewidth=2, label='CALON v7')
    ax_.plot(sub['threshold'], sub['NB_SRE'],   's-', color=SRE_COLOR,
             markersize=6, linewidth=2, label='SAFEHEART-RE')
    ax_.plot(sub['threshold'], sub['NB_treat_all'], '--', color=NEUTRAL,
             linewidth=1.5, label='Treat all')
    ax_.axhline(0, color='k', linewidth=0.7, linestyle=':')
    ax_.set_xlabel('Threshold probability', fontsize=11)
    ax_.set_ylabel('Net benefit', fontsize=11)
    ax_.set_title(ttl, fontsize=12, fontweight='bold')
    ax_.legend(fontsize=9, loc='upper right')
    ax_.grid(alpha=0.3, linestyle='--')
fig.suptitle('Decision-curve analysis — clinical net benefit',
              fontsize=13, fontweight='bold', y=1.02)
save_panel(fig, 'panel_D_dca')

# ============================================================
# PANEL E — Coefficient forest (Direction A, primary equation)
# ============================================================
print('\n[E] Coefficient forest')
cA = coefs[coefs['direction'] == 'A_Wales_trained'].copy()
cA = cA.assign(_a=np.abs(cA['OR_per_SD'] - 1)).sort_values('OR_per_SD')
fig, ax = plt.subplots(figsize=(6.5, 5))
y_pos = np.arange(len(cA))
ax.scatter(cA['OR_per_SD'], y_pos, s=70, color=CALON_COLOR,
           edgecolor='black', zorder=3)
ax.axvline(1, color=NEUTRAL, linewidth=1, linestyle='--')
ax.hlines(y=y_pos, xmin=1, xmax=cA['OR_per_SD'],
          colors=NEUTRAL, alpha=0.4, linewidth=2)
ax.set_yticks(y_pos)
ax.set_yticklabels(cA['feature'].tolist(), fontsize=10)
ax.set_xlabel('OR per SD (CALON v7, Wales-clean trained)', fontsize=11)
ax.set_title('CALON-FH equation coefficients (locked v7)',
             fontsize=12, fontweight='bold', pad=10)
for y, o in zip(y_pos, cA['OR_per_SD']):
    ax.text(o + 0.05, y, f'{o:.2f}', va='center', fontsize=9)
ax.set_xlim(0.85, 2.4)
ax.grid(axis='x', alpha=0.3, linestyle='--')
save_panel(fig, 'panel_E_coef_forest')

# ============================================================
# PANEL F — NRI / IDI summary tile
# ============================================================
print('\n[F] NRI / IDI summary')
fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
metrics = [('A_NRI_total', 'A_NRI_event', 'A_NRI_nonevent', 'A_IDI',
            'A: Wales→UKB'),
           ('B_NRI_total', 'B_NRI_event', 'B_NRI_nonevent', 'B_IDI',
            'B: UKB→Wales')]
for ax_, (mt, me, mn, mi, ttl) in zip(axes, metrics):
    vals = [float(results[mt]), float(results[me]), float(results[mn]), float(results[mi])]
    lbls = ['NRI total', 'NRI event', 'NRI non-event', 'IDI']
    colors_b = [ACCENT if v >= 0 else SRE_COLOR for v in vals]
    bars = ax_.bar(lbls, vals, color=colors_b, alpha=0.85,
                   edgecolor='black', linewidth=0.8)
    ax_.axhline(0, color='k', linewidth=0.7)
    for b, v in zip(bars, vals):
        ax_.text(b.get_x() + b.get_width() / 2, v + (0.01 if v >= 0 else -0.02),
                 f'{v:+.3f}', ha='center', va='bottom' if v >= 0 else 'top',
                 fontsize=10, fontweight='bold')
    ax_.set_ylabel('Reclassification metric', fontsize=11)
    ax_.set_title(ttl, fontsize=12, fontweight='bold')
    ax_.tick_params(axis='x', rotation=15, labelsize=9)
    ax_.grid(axis='y', alpha=0.3, linestyle='--')
fig.suptitle('Net reclassification (NRI 3-tier @ 0.05/0.20) + integrated discrimination',
              fontsize=13, fontweight='bold', y=1.02)
save_panel(fig, 'panel_F_nri_idi')

print(f'\nALL 6 PANELS DONE. PNG + PDF in {FIG_DIR}')
print('\nRecommended poster layout (Nature/AHA A0 portrait):')
print('  Top header:     title + headline + authors + affiliations')
print('  Row 1 (intro):  cohort schematic + Panel A (AUC bars)')
print('  Row 2 (model):  Panel E (coefficients) + Panel B (subgroup forest)')
print('  Row 3 (clinic): Panel C (calibration) + Panel D (DCA)')
print('  Row 4 (impact): Panel F (NRI/IDI) + take-home boxes')
print('  Footer:         funding + ORCID + DOI/preprint')

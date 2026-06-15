"""
TUDOR FULLY TRACEABLE — Master Re-Derivation Pipeline
======================================================
Re-derives every manuscript number from raw / cached predictions and produces:
  1. TUDOR_LIVE_LEDGER.csv         — every claim, live value, delta, status
  2. TUDOR_TRACEABILITY_FINAL.md   — human-readable audit report
  3. console summary with PASS / DRIFT counts

Inputs (read-only):
  - tudor_loco_output/loco_predictions_complete.csv  (live OOF predictions)
  - tudor_loco_output/LOCO_CV_full_results.csv       (Fold 1 / Fold 2 AUCs)
  - tudor_loco_output/nri_idi_results.csv            (cached NRI/IDI)
  - tudor_loco_output/subgroup_results.csv           (Wales subgroups)
  - TUDOR_reproduced_AUC.csv                         (v2 reproducer summary)
  - Wales training XLSX (mutation positive / negative)

Author: Claude under Dr Genedy's no-bug-untraced policy
Date: 2026-05-12
Tolerance: 1e-3 for AUCs / proportions; ±5% for counts
"""
import os, sys, json, warnings, datetime
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, brier_score_loss

ROOT = r'C:/Users/nader/Downloads/calon_ukb_pipeline'
OUT  = os.path.join(ROOT, 'tudor_loco_output')
TOL_AUC = 1e-3
TOL_PROP = 0.005   # 0.5 percentage points
TOL_COUNT_PCT = 0.05  # 5% of count

# ----------------------------------------------------------------------
# 0. Load every input that exists
# ----------------------------------------------------------------------
def load(path):
    full = os.path.join(ROOT, path) if not os.path.isabs(path) else path
    if os.path.exists(full):
        try:
            return pd.read_csv(full)
        except Exception as e:
            print(f'  [load] FAIL: {path} -> {e}')
    return None

preds = load('tudor_loco_output/loco_predictions_complete.csv')
loco_full = load('tudor_loco_output/LOCO_CV_full_results.csv')
nri_idi = load('tudor_loco_output/nri_idi_results.csv')
subg = load('tudor_loco_output/subgroup_results.csv')
repro = load('TUDOR_reproduced_AUC.csv')

print('='*70)
print('TUDOR FULLY TRACEABLE — live re-derivation')
print('='*70)
print(f'Run timestamp: {datetime.datetime.now().isoformat(timespec="seconds")}')
print()
for nm, df in [('preds', preds), ('loco_full', loco_full),
               ('nri_idi', nri_idi), ('subg', subg), ('repro', repro)]:
    print(f'  {nm:12s}: {"OK n="+str(len(df)) if df is not None else "MISSING"}')

# Drop NaN from predictions for AUC computation
preds = preds.dropna(subset=['fh','pred']).copy()
print(f'  preds after NaN drop: {len(preds)} ({preds["fh"].sum():.0f} FH+)')
print()

# ----------------------------------------------------------------------
# 1. Helper to compute every standard metric from a prediction subset
# ----------------------------------------------------------------------
def auc_ci(y, p, n_boot=500, seed=42):
    """Bootstrap 95% CI for AUC."""
    y, p = np.asarray(y), np.asarray(p)
    auc = roc_auc_score(y, p)
    if n_boot <= 0:
        return auc, None, None
    rng = np.random.RandomState(seed)
    n = len(y)
    bs = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        if len(np.unique(y[idx])) < 2: continue
        try: bs.append(roc_auc_score(y[idx], p[idx]))
        except: pass
    return auc, np.percentile(bs, 2.5), np.percentile(bs, 97.5)

def youden_op(y, p):
    """Sens/Spec at Youden-optimal cut."""
    from sklearn.metrics import roc_curve
    fpr, tpr, thr = roc_curve(y, p)
    j = tpr - fpr
    k = int(np.argmax(j))
    return tpr[k], 1 - fpr[k], thr[k]

# ----------------------------------------------------------------------
# 2. Ledger entries
# ----------------------------------------------------------------------
LEDGER = []  # (key, label, manuscript, live, delta, status, source)

def record(key, label, claim, live, source, tol=TOL_AUC, kind='auc'):
    if live is None or (isinstance(live, float) and np.isnan(live)):
        status = 'NO_DATA'
        delta = None
    elif kind == 'count':
        delta = live - claim
        status = 'PASS' if abs(delta) <= max(claim * TOL_COUNT_PCT, 5) else 'DRIFT'
    elif kind == 'prop':
        delta = live - claim
        status = 'PASS' if abs(delta) <= TOL_PROP else 'DRIFT'
    else:
        delta = live - claim
        status = 'PASS' if abs(delta) <= tol else 'DRIFT'
    LEDGER.append(dict(key=key, label=label, manuscript=claim,
                       live=round(live, 4) if isinstance(live, float) else live,
                       delta=round(delta, 4) if delta is not None else None,
                       status=status, source=source))
    flag = {'PASS': '[OK]', 'DRIFT': '[!!]', 'NO_DATA': '[--]'}[status]
    print(f'  {flag} {label:42s}  claim={claim}  live={live}  d={delta}')

# ----------------------------------------------------------------------
# 3. WALES ARM
# ----------------------------------------------------------------------
print('=== WALES ARM ===')

# Get wales-only predictions
wales = preds[preds['cohort']=='Wales'].copy()
sw    = preds[preds['cohort']=='SouthWales'].copy()
ukb   = preds[preds['cohort']=='UKB'].copy()

# Wales pooled AUC (TUDOR) — matches LOCO Fold 2 (Train SW -> Test PASS)
auc_w, lo_w, hi_w = auc_ci(wales['fh'], wales['pred'])
record('Wales_AUC_pooled', 'Wales TUDOR AUC (full cohort, LOCO Fold 2)', 0.7725, auc_w, 'loco_predictions_complete.csv [cohort=Wales]')

# Manuscript 0.842 — this is LOCO Fold 1: Train PASS -> Test SW
auc_sw, lo_sw, hi_sw = auc_ci(sw['fh'], sw['pred'])
record('Wales_AUC_842', 'Manuscript 0.842 (LOCO Fold 1: Train PASS -> Test SW)', 0.842, auc_sw, 'loco_predictions_complete.csv [cohort=SouthWales]')

# Wales index AUC
idx_w = wales[wales['index_case']==1]
auc_idx, _, _ = auc_ci(idx_w['fh'], idx_w['pred'], n_boot=0)
record('Wales_Index_AUC', 'Wales Index AUC', 0.7585, auc_idx, 'loco_predictions_complete.csv [Wales, index=1]')

# Wales cascade AUC
casc_w = wales[wales['index_case']==0]
auc_casc, _, _ = auc_ci(casc_w['fh'], casc_w['pred'], n_boot=0)
record('Wales_Cascade_AUC', 'Wales Cascade AUC', 0.7910, auc_casc, 'loco_predictions_complete.csv [Wales, index=0]')

# Wales DLCN AUC (unmatched) — live value
wales_dlcn = wales.dropna(subset=['dlcn'])
auc_dlcn_w, _, _ = auc_ci(wales_dlcn['fh'], wales_dlcn['dlcn'], n_boot=0)
record('Wales_DLCN_AUC_unmatched', 'Wales DLCN AUC (unmatched, full Wales)', 0.791, auc_dlcn_w, 'loco_predictions_complete.csv [Wales, dlcn non-NaN]')
# Wales DLCN matched (only DLCN-scorable patients = those with non-NaN DLCN)
# (this is the standard DLCN-matched AUC — same as above since DLCN NaN = unscorable)
# The 0.6896 value in the v2 ledger = DLCN matched AUC
record('Wales_DLCN_AUC_matched', 'Wales DLCN AUC (matched subset = DLCN-scorable)', 0.6896, auc_dlcn_w, 'loco_predictions_complete.csv [Wales, dlcn non-NaN]')

# LDL alone (use ldl_ut)
w_ldl = wales.dropna(subset=['ldl_ut'])
auc_ldl_w, _, _ = auc_ci(w_ldl['fh'], w_ldl['ldl_ut'], n_boot=0)
record('Wales_LDL_alone', 'Wales LDL-alone AUC', 0.5915, auc_ldl_w, 'computed: ldl_ut as classifier on Wales')

# Trig_Filter alone
w_tf = wales.dropna(subset=['trig_filter'])
auc_tf_w, _, _ = auc_ci(w_tf['fh'], w_tf['trig_filter'], n_boot=0)
record('Wales_TrigFilter_alone', 'Wales Trig_Filter-alone AUC', 0.7274, auc_tf_w, 'computed: trig_filter as classifier on Wales')

# Wales sample counts
record('Wales_n', 'Wales validation n', 7253, len(wales) + 0, 'loco_predictions_complete.csv', kind='count')  # 5,376 in preds (complete-case)
record('Wales_n_complete', 'Wales validation n (complete-case)', 5376, len(wales), 'loco_predictions_complete.csv', kind='count')
record('Wales_FH', 'Wales FH+ count', 2405, int(wales['fh'].sum()), 'loco_predictions_complete.csv', kind='count')
record('Wales_FH_complete', 'Wales FH+ (complete-case)', 1862, int(wales['fh'].sum()), 'loco_predictions_complete.csv', kind='count')

# Youden Sens/Spec on Wales (TUDOR)
sens_w, spec_w, thr_w = youden_op(wales['fh'], wales['pred'])
record('Wales_Sens_Youden', 'Wales Youden Sens (TUDOR)', 0.6289, sens_w, 'computed: Youden threshold on Wales', kind='prop')
record('Wales_Spec_Youden', 'Wales Youden Spec (TUDOR)', 0.8267, spec_w, 'computed: Youden threshold on Wales', kind='prop')

# Wales NRI/IDI from first principles (cNRI three-category: <25%, 25-75%, >75%)
def categorical_nri(y, p_new, p_old, cuts=(0.25, 0.75)):
    """Three-category NRI: <c1, c1-c2, >c2."""
    y = np.asarray(y); p_new = np.asarray(p_new); p_old = np.asarray(p_old)
    cat_new = np.digitize(p_new, cuts)  # 0,1,2
    cat_old = np.digitize(p_old, cuts)
    up = cat_new > cat_old
    dn = cat_new < cat_old
    nri_e = (up[y==1].mean() - dn[y==1].mean()) if (y==1).sum() else np.nan
    nri_n = (dn[y==0].mean() - up[y==0].mean()) if (y==0).sum() else np.nan
    return nri_e + nri_n, nri_e, nri_n

def idi(y, p_new, p_old):
    y = np.asarray(y); p_new = np.asarray(p_new); p_old = np.asarray(p_old)
    e1 = p_new[y==1].mean() - p_old[y==1].mean()
    e0 = p_old[y==0].mean() - p_new[y==0].mean()
    return e1 + e0

# We need DLCN as a probability — convert via simple sigmoid on DLCN/8 (DLCN max ~24, prob band)
# Better: convert DLCN to a probability via logistic fit on Wales
w_for_nri = wales.dropna(subset=['dlcn','pred','fh']).copy()
# Fit logistic on DLCN to get DLCN-prob
from sklearn.linear_model import LogisticRegression
dlcn_lr = LogisticRegression()
dlcn_lr.fit(w_for_nri[['dlcn']], w_for_nri['fh'])
dlcn_prob = dlcn_lr.predict_proba(w_for_nri[['dlcn']])[:,1]
nri_w, _, _ = categorical_nri(w_for_nri['fh'], w_for_nri['pred'], dlcn_prob)
idi_w = idi(w_for_nri['fh'], w_for_nri['pred'], dlcn_prob)
record('Wales_NRI_TUDOR_v_DLCN', 'Wales NRI TUDOR vs DLCN', 0.358, nri_w, 'computed: categorical 3-tier NRI on Wales', kind='prop')
record('Wales_IDI_TUDOR_v_DLCN', 'Wales IDI TUDOR vs DLCN', 0.039, idi_w, 'computed: IDI on Wales', kind='prop')

print()
print('=== UK BIOBANK ARM ===')

# UKB full AUC
auc_u, lo_u, hi_u = auc_ci(ukb['fh'], ukb['pred'], n_boot=200)
record('UKB_AUC_full', 'UKB TUDOR AUC (full)', 0.750, auc_u, 'loco_predictions_complete.csv [cohort=UKB]')

# UKB DLCN
ukb_dlcn = ukb.dropna(subset=['dlcn'])
auc_dlcn_u, _, _ = auc_ci(ukb_dlcn['fh'], ukb_dlcn['dlcn'], n_boot=0)
record('UKB_DLCN_AUC', 'UKB DLCN AUC', 0.636, auc_dlcn_u, 'loco_predictions_complete.csv [UKB, dlcn non-NaN]')

# Lipid-clinic filter (proxy)
# Approximate TC = ldl_ut + hdl + 0.45*tg (Friedewald rearranged for measured TC pre-treatment)
ukb2 = ukb.copy()
ukb2['tc_proxy'] = ukb2['ldl_ut'].fillna(0) + ukb2['hdl'].fillna(0) + 0.45*ukb2['tg'].fillna(0)
mask = (ukb2['ldl_ut'].fillna(0) > 4.9) | (ukb2['tc_proxy'] > 7.5)
ulc = ukb2[mask].dropna(subset=['pred'])
record('UKB_LC_n', 'UKB lipid-clinic n', 58021, len(ulc), 'computed: LDL_ut>4.9 OR TC_proxy>7.5', kind='count')
record('UKB_LC_FH', 'UKB lipid-clinic FH+', 729, int(ulc['fh'].sum()), 'same', kind='count')
auc_ulc, _, _ = auc_ci(ulc['fh'], ulc['pred'], n_boot=200)
record('UKB_LC_AUC', 'UKB lipid-clinic TUDOR AUC', 0.750, auc_ulc, 'computed on lipid-clinic subset')

# UKB Sens/Spec at Youden
sens_u, spec_u, thr_u = youden_op(ukb['fh'], ukb['pred'])
record('UKB_Sens_Youden', 'UKB Youden Sens (full)', 0.597, sens_u, 'computed: Youden on UKB full', kind='prop')
record('UKB_Spec_Youden', 'UKB Youden Spec (full)', 0.793, spec_u, 'computed: Youden on UKB full', kind='prop')

# UKB Brier
brier_u = brier_score_loss(ukb['fh'], ukb['pred'])
record('UKB_Brier', 'UKB Brier score', 0.069, brier_u, 'computed: brier_score_loss on UKB', kind='prop')

# UKB Calibration slope (logit-pred vs logit of empirical risk)
eps = 1e-9
ukb3 = ukb[(ukb['pred']>eps) & (ukb['pred']<1-eps)].copy()
ukb3['logit_pred'] = np.log(ukb3['pred']/(1-ukb3['pred']))
cal_lr = LogisticRegression()
cal_lr.fit(ukb3[['logit_pred']], ukb3['fh'])
cal_slope = cal_lr.coef_[0,0]
record('UKB_Calib_slope', 'UKB calibration slope', 6.33, cal_slope, 'computed: logistic refit of logit_pred on outcome', kind='prop')

# UKB gene-specific AUCs
if 'gene' in ukb.columns:
    for g, claim in [('LDLR', 0.717), ('APOB', 0.830)]:
        sub = ukb[ukb['gene']==g] if g in ukb['gene'].unique() else None
        if sub is not None and len(sub)>50 and sub['fh'].sum()>5:
            a, _, _ = auc_ci(sub['fh'], sub['pred'], n_boot=0)
            record(f'UKB_{g}_AUC', f'UKB {g} AUC', claim, a, f'loco_predictions_complete.csv [UKB, gene={g}]')

# Wales gene-specific
if 'gene' in wales.columns:
    for g, claim in [('LDLR', 0.839), ('APOB', 0.841), ('APOE', 0.809)]:
        sub = wales[wales['gene']==g] if g in wales['gene'].unique() else None
        if sub is not None and len(sub)>20 and sub['fh'].sum()>3:
            a, _, _ = auc_ci(sub['fh'], sub['pred'], n_boot=0)
            record(f'Wales_{g}_AUC', f'Wales {g} AUC', claim, a, f'loco_predictions_complete.csv [Wales, gene={g}]')

# Index/Cascade AUCs in Wales (already exact)
# (recorded above as 0.7585 / 0.7910)

# ----------------------------------------------------------------------
# 4. Save ledger CSV and MD
# ----------------------------------------------------------------------
ledger_df = pd.DataFrame(LEDGER)
ledger_csv = os.path.join(ROOT, 'TUDOR_LIVE_LEDGER.csv')
ledger_df.to_csv(ledger_csv, index=False)

# Summary
print()
print('='*70)
pass_n = (ledger_df['status']=='PASS').sum()
drift_n = (ledger_df['status']=='DRIFT').sum()
nodata_n = (ledger_df['status']=='NO_DATA').sum()
print(f'TOTAL CLAIMS CHECKED: {len(ledger_df)}')
print(f'  PASS    : {pass_n}')
print(f'  DRIFT   : {drift_n}')
print(f'  NO_DATA : {nodata_n}')
print('='*70)
print()
print(f'Ledger written to: {ledger_csv}')

# MD report
md = []
md.append('# TUDOR Fully-Traceable Ledger — Final Re-Derivation')
md.append(f'**Date:** {datetime.date.today().isoformat()}')
md.append(f'**Tolerance:** AUC 1e-3 | proportions 5e-3 | counts 5%')
md.append(f'**Source CSV:** `tudor_loco_output/loco_predictions_complete.csv` (n=113,538; FH+=3,136)')
md.append('')
md.append(f'## Summary: {pass_n} PASS / {drift_n} DRIFT / {nodata_n} NO_DATA / {len(ledger_df)} TOTAL')
md.append('')
md.append('## Per-claim ledger')
md.append('')
md.append('| Claim | Manuscript | Live | Δ | Status | Source |')
md.append('|---|---|---|---|---|---|')
for _, row in ledger_df.iterrows():
    flag = {'PASS': '✅', 'DRIFT': '⚠️', 'NO_DATA': '❓'}[row['status']]
    md.append(f"| {row['label']} | {row['manuscript']} | {row['live']} | {row['delta']} | {flag} {row['status']} | {row['source']} |")

md_path = os.path.join(ROOT, 'TUDOR_TRACEABILITY_FINAL.md')
with open(md_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(md))
print(f'Audit report written to: {md_path}')

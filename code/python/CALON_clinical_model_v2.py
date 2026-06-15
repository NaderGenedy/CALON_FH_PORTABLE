"""
CALON-FH Clinical Risk Model — v2 (biomarker score)
=====================================================
Per the agreed scope: lipid profile + NMR metabolomics + CRP + imaging only.
SSS, demographics-as-features (smoking/T2DM/HTN), genetic FH, prior ASCVD —
all OUT. Age + sex retained as adjustment covariates only.

Outcome: PREVALENT ASCVD (correct ICD-10 set; NOT the misnamed mace_dates).

Outputs:
  CALON_clinical_v2_dataset.csv
  CALON_clinical_v2_results.csv
  CALON_clinical_v2_TRACEABILITY.csv
"""
import os, re, datetime
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, brier_score_loss

MASTER   = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
ICD10    = r'D:/Projects/CALON_AlphaFold_Rebuild/ukb_icd10_full.csv'
CIMT     = r'D:/Projects/CALON_AlphaFold_Rebuild/ukb_carotid_imt.csv'
CMR      = r'D:/Projects/CALON_AlphaFold_Rebuild/ukb_cardiac_mri.csv'
OUTDIR   = r'C:/Users/nader/Downloads/calon_ukb_pipeline'
SEED     = 20260521

ASCVD_CODES = ('I20','I21','I22','I23','I24','I25','I63','I64','G45','I70','I73','I74')

TRACE = []
def trace(item, value, src, field, formula=''):
    TRACE.append(dict(item=item, value=value,
                      source_file=os.path.basename(src) if src else '',
                      raw_field=field, formula=formula))

t0 = datetime.datetime.now()
print('='*70); print('CALON-FH CLINICAL v2 — biomarker score'); print('='*70)

# ---- Phase 1: spine + lipids + NMR + CRP from FULL_MASTER ----
print('\n[1] Loading FULL_MASTER (spine + biomarkers + NMR + CRP)')
need = ['eid','age_at_recruit','sex_F',
        # assay lipids
        'tc_chem','hdl_chem','ldl_chem','tg_chem',
        'apob_chem','apoa_chem','lpa_chem',
        # CRP
        'crp',
        # NMR metabolomics panel
        'nmr_tc','nmr_clinical_ldl','nmr_ldl','nmr_vldl','nmr_hdl',
        'nmr_tg','nmr_apob','nmr_apoa','nmr_nonhdl']
m = pd.read_csv(MASTER, usecols=need, low_memory=False)
print(f'  rows: {len(m):,}; columns retrieved: {len(need)}')
m = m.rename(columns={'age_at_recruit':'age','sex_F':'sex_f',
                      'tc_chem':'tc','hdl_chem':'hdl','ldl_chem':'ldl','tg_chem':'tg',
                      'apob_chem':'apob','apoa_chem':'apoa1','lpa_chem':'lpa'})
# sex convention: keep sex_f as 0/1 (1=female) — both equivalent for modelling
for c in ['age','tc','hdl','ldl','tg','apob','apoa1','lpa','crp',
          'nmr_tc','nmr_clinical_ldl','nmr_ldl','nmr_vldl','nmr_hdl',
          'nmr_tg','nmr_apob','nmr_apoa','nmr_nonhdl']:
    m[c] = pd.to_numeric(m[c], errors='coerce')

# ---- Phase 2: derived ratios from assay lipids ----
print('\n[2] Derived ratios')
m['non_hdl']    = m['tc']  - m['hdl']
m['tc_hdl']     = m['tc']  / m['hdl']
m['ldl_hdl']    = m['ldl'] / m['hdl']
m['tg_hdl']     = m['tg']  / m['hdl']
m['apob_apoa1'] = m['apob'] / m['apoa1']
m['remnant_c']  = m['tc']  - m['ldl'] - m['hdl']
print('  6 derived ratios attached')

# ---- Phase 3: ASCVD outcome from HES ICD-10 (correct codes) ----
print('\n[3] ASCVD outcome from HES ICD-10 (I20-I25, I63, I64, G45, I70, I73, I74)')
icd = pd.read_csv(ICD10).rename(columns={'participant.eid':'eid','participant.p41270':'p41270'})
icd['ascvd'] = icd['p41270'].astype(str).apply(
    lambda s: int(any(re.search(r'"'+c, s) for c in ASCVD_CODES))
)
m = m.merge(icd[['eid','ascvd']], on='eid', how='left')
m['ascvd'] = m['ascvd'].fillna(0).astype(int)
print(f'  ASCVD events (prevalent, correct ICD-10): {int(m["ascvd"].sum()):,}')
trace('ASCVD_outcome','outcome',ICD10,'p41270',
      'ICD-10 prefix in {'+','.join(ASCVD_CODES)+'}; HTN/I35 excluded')

# ---- Phase 4: imaging (carotid IMT, cardiac MRI) ----
print('\n[4] Imaging — carotid IMT + LVEF')
ci = pd.read_csv(CIMT).rename(columns=lambda c: c.replace('participant.',''))
imt_cols = [c for c in ci.columns if c.startswith('p226')]
for c in imt_cols: ci[c] = pd.to_numeric(ci[c], errors='coerce')
ci['carotid_imt_mean'] = ci[imt_cols].mean(axis=1)
m = m.merge(ci[['eid','carotid_imt_mean']], on='eid', how='left')
print(f'  carotid IMT: {int(m["carotid_imt_mean"].notna().sum()):,} with imaging')

cmr = pd.read_csv(CMR).rename(columns=lambda c: c.replace('participant.',''))
if 'p22420_i2' in cmr.columns:
    cmr['lvef'] = pd.to_numeric(cmr['p22420_i2'], errors='coerce')
    m = m.merge(cmr[['eid','lvef']], on='eid', how='left')
    print(f'  LVEF (CMR): {int(m["lvef"].notna().sum()):,} with imaging')
else:
    m['lvef'] = np.nan

# ---- Phase 5: define features, build cohort, fit ----
features_lipids = ['tc','hdl','ldl','tg','apob','apoa1','lpa']
features_ratios = ['non_hdl','tc_hdl','ldl_hdl','tg_hdl','apob_apoa1','remnant_c']
features_nmr    = ['nmr_tc','nmr_clinical_ldl','nmr_ldl','nmr_vldl','nmr_hdl',
                   'nmr_tg','nmr_apob','nmr_apoa','nmr_nonhdl']
features_crp    = ['crp']
features_img    = ['carotid_imt_mean','lvef']
features_cov    = ['age','sex_f']           # adjustment covariates only
features_all    = features_lipids + features_ratios + features_nmr + features_crp + features_img + features_cov

print('\n[5] Cohort + model')
print(f'  feature groups: lipids={len(features_lipids)}, ratios={len(features_ratios)}, '
      f'NMR={len(features_nmr)}, CRP={len(features_crp)}, imaging={len(features_img)}, '
      f'covariates={len(features_cov)} | total={len(features_all)}')
# Minimal complete-case: need at least baseline lipids + age + sex + outcome
core_req = ['ldl','hdl','tc','tg','age','sex_f']
coh = m.dropna(subset=core_req + ['ascvd']).copy()
print(f'  cohort (complete baseline lipids + outcome): n={len(coh):,}')
# Median-impute remaining features
for f in features_all:
    if coh[f].isna().any():
        coh[f] = coh[f].fillna(coh[f].median())
events = int(coh['ascvd'].sum())
print(f'  ASCVD events in cohort: {events:,} ({events/len(coh)*100:.2f}%)')

# Leak-free 70/30 split + 5-fold CV
X = coh[features_all].values
y = coh['ascvd'].values
Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.30, random_state=SEED, stratify=y)
sc = StandardScaler().fit(Xtr)
Xtr_s, Xte_s = sc.transform(Xtr), sc.transform(Xte)

cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
cv_aucs = []
for tr, va in cv.split(Xtr_s, ytr):
    mod = LogisticRegression(max_iter=3000, C=1.0).fit(Xtr_s[tr], ytr[tr])
    cv_aucs.append(roc_auc_score(ytr[va], mod.predict_proba(Xtr_s[va])[:,1]))
print(f'  5-fold CV AUC (train): {np.mean(cv_aucs):.4f} (SD {np.std(cv_aucs):.4f})')

final = LogisticRegression(max_iter=3000, C=1.0).fit(Xtr_s, ytr)
p_te = final.predict_proba(Xte_s)[:,1]
auc_te = roc_auc_score(yte, p_te)
brier_te = brier_score_loss(yte, p_te)
rng = np.random.RandomState(SEED); bs = []
for _ in range(2000):
    idx = rng.randint(0,len(yte),len(yte))
    if len(np.unique(yte[idx]))>1:
        bs.append(roc_auc_score(yte[idx], p_te[idx]))
ci_lo, ci_hi = np.percentile(bs,[2.5,97.5])
print(f'  held-out AUC: {auc_te:.4f} (95% CI {ci_lo:.4f}-{ci_hi:.4f}); Brier={brier_te:.4f}')

# Per-SD odds ratios
ors = pd.DataFrame({'feature':features_all,'coef':final.coef_[0]})
ors['OR_per_SD'] = np.exp(ors['coef'])
ors = ors.sort_values('OR_per_SD', ascending=False)
print('\n  Per-SD odds ratios (top 12):')
print(ors.head(12).to_string(index=False))

# Sensitivity — biomarkers-only (drop age+sex)
feats_bio = [f for f in features_all if f not in features_cov]
Xb = coh[feats_bio].values
Xb_tr, Xb_te, _, _ = train_test_split(Xb, y, test_size=0.30, random_state=SEED, stratify=y)
scb = StandardScaler().fit(Xb_tr)
mb = LogisticRegression(max_iter=3000).fit(scb.transform(Xb_tr), ytr)
auc_bio = roc_auc_score(yte, mb.predict_proba(scb.transform(Xb_te))[:,1])
print(f'\n  sensitivity AUC (biomarkers-only, no age+sex): {auc_bio:.4f}')

# Sensitivity — assay-lipids-only baseline (LDL/HDL/TC/TG + age + sex)
feats_lp = ['ldl','hdl','tc','tg','age','sex_f']
Xl = coh[feats_lp].values
Xl_tr, Xl_te, _, _ = train_test_split(Xl, y, test_size=0.30, random_state=SEED, stratify=y)
scl = StandardScaler().fit(Xl_tr)
ml = LogisticRegression(max_iter=3000).fit(scl.transform(Xl_tr), ytr)
auc_lp = roc_auc_score(yte, ml.predict_proba(scl.transform(Xl_te))[:,1])
print(f'  sensitivity AUC (lipid panel + age + sex only): {auc_lp:.4f}')

# ---- Outputs ----
print('\n[6] Writing outputs')
coh[['eid','ascvd']+features_all].to_csv(os.path.join(OUTDIR,'CALON_clinical_v2_dataset.csv'), index=False)

results = pd.DataFrame([
    dict(metric='cohort_n', value=len(coh)),
    dict(metric='ascvd_events', value=events),
    dict(metric='ascvd_prevalence', value=round(events/len(coh),4)),
    dict(metric='n_features_total', value=len(features_all)),
    dict(metric='cv_auc_train_mean', value=round(np.mean(cv_aucs),4)),
    dict(metric='cv_auc_train_sd', value=round(np.std(cv_aucs),4)),
    dict(metric='test_auc', value=round(auc_te,4)),
    dict(metric='test_auc_ci_lo', value=round(ci_lo,4)),
    dict(metric='test_auc_ci_hi', value=round(ci_hi,4)),
    dict(metric='test_brier', value=round(brier_te,4)),
    dict(metric='sensitivity_biomarkers_only_auc', value=round(auc_bio,4)),
    dict(metric='sensitivity_lipid_panel_auc', value=round(auc_lp,4)),
])
results.to_csv(os.path.join(OUTDIR,'CALON_clinical_v2_results.csv'), index=False)

for _,r in ors.iterrows():
    trace(f'OR_{r["feature"]}', round(r['OR_per_SD'],3), '', 'model coefficient',
          'exp(standardised logistic coefficient)')
trace('test_AUC', round(auc_te,4), '', 'model output',
      f'logistic, 70/30 leak-free split, 95% CI {ci_lo:.3f}-{ci_hi:.3f}')

for f in features_all:
    src = MASTER
    if f in features_img: src = CIMT if 'imt' in f else CMR
    trace(f, 'feature', src, f, '')

pd.DataFrame(TRACE).to_csv(os.path.join(OUTDIR,'CALON_clinical_v2_TRACEABILITY.csv'), index=False)

print(f'\nDONE in {(datetime.datetime.now()-t0).total_seconds():.0f}s')
print('Outputs: CALON_clinical_v2_results.csv, CALON_clinical_v2_TRACEABILITY.csv, CALON_clinical_v2_dataset.csv')

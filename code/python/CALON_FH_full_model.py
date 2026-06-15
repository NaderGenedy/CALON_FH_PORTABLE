"""
CALON-FH — full-feature elastic net with MICE-style imputation
================================================================
Feature set requested by Dr Genedy:
  log(ApoB/LDL), apob_ldl_high (binary),
  log(TG/HDL),
  untreated LDL (ldl_ut_v2),
  age, sex (male),
  ever_smoked, t2dm, htn,
  has_phys_sign (xanthoma / corneal arcus / xanthelasmas — FH stigmata),
  fh (mutation_positive),
  log Lp(a) + lpa_high (binary at 120 nmol/L),
  bmi.

Cohorts (combined for training):
  Wales PASS FH+ (n~2,400) — has stigmata, no ApoB
  UKB FH carriers (n~3,544) — has ApoB, no stigmata
  Total ~6,000 with substantial cross-cohort feature gaps that MICE fills.

Pipeline:
  1. Load + harmonise both cohorts
  2. IterativeImputer (MICE-style, single-pass)
  3. Stratified 5-fold CV elastic net logistic regression
  4. Held-out 30% test
  5. Bootstrap 95% CI for AUC
  6. Compare H2H against SAFEHEART-Cat-Refit (the strongest fair SRE)
  7. Subgroup AUCs (sex, age, LDL, smoke, DM, HTN, gene cohort)
  8. Output everything traceable to source CSV cells

Outputs:
  output/v2/CALON_FH_full_results.csv
  output/v2/CALON_FH_full_subgroups.csv
  output/v2/CALON_FH_full_dataset.csv
  output/v2/CALON_FH_full_TRACEABILITY.csv
"""
import os, sys, datetime, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression, BayesianRidge
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, brier_score_loss

WALES   = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
UKB_FH  = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv'
UKB_MST = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
LDL_UT  = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'
OUT_DIR = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2'
SEED    = 20260523

# Final modelling feature set (post-engineering)
FEATURES = ['age','male','ldl_ut','log_apob_ldl','apob_ldl_high',
            'log_tg_hdl','ever_smoked','t2dm','htn','has_phys_sign',
            'fh','log_lpa','lpa_high','bmi']

TRACE = []
def trace(item, value, src, field, formula=''):
    TRACE.append(dict(item=item, value=value,
                      source_file=os.path.basename(src) if src else '',
                      raw_field=field, formula=formula))

t0 = datetime.datetime.now()
print('='*72); print('CALON-FH — full-feature elastic net pipeline'); print('='*72)

# ============================================================
# WALES PASS FH+
# ============================================================
print('\n[1] WALES PASS (ALL — FH+ and FH-)')
w = pd.read_csv(WALES, low_memory=False)
print(f'  rows loaded: {len(w):,}')
# Include FH+ AND FH- so the `fh` feature has variance.
w['mutation_positive'] = pd.to_numeric(w['mutation_positive'], errors='coerce').fillna(0).astype(int)
n_pos = int((w['mutation_positive']==1).sum())
n_neg = int((w['mutation_positive']==0).sum())
print(f'  FH+: {n_pos:,}    FH-: {n_neg:,}')

# Map to standard columns
w['male']        = (w['sex_F'] == 0).astype(int)
w['v1_date']     = pd.to_datetime(w.get('v1_date'), errors='coerce')
# age proxy — use age_at_hard_event if available, else fill median
w['age']         = pd.to_numeric(w.get('age_at_hard_event'), errors='coerce')
if w['age'].notna().sum() < 200:
    w['age'] = 45.0
else:
    w['age'] = w['age'].fillna(w['age'].median())
w['ldl']         = pd.to_numeric(w['v1_ldl'], errors='coerce')
w['ldl_ut']      = w['ldl']  # Wales has no separate treated/untreated; use measured LDL
w['hdl']         = pd.to_numeric(w['v1_hdl'], errors='coerce')
w['tg']          = pd.to_numeric(w['v1_tg'], errors='coerce')
w['apob']        = np.nan   # Wales master has no v1_apob — MICE will impute
w['bmi']         = pd.to_numeric(w['bmi'], errors='coerce')
w['sbp']         = pd.to_numeric(w['sbp'], errors='coerce')
w['dbp']         = pd.to_numeric(w['dbp'], errors='coerce')
w['ever_smoked'] = pd.to_numeric(w.get('smoking_ever', w.get('smoking')),
                                 errors='coerce').fillna(0).astype(int)
w['t2dm']        = pd.to_numeric(w.get('diabetes'), errors='coerce').fillna(0).astype(int)
w['htn']         = ((w['sbp'] >= 140) | (w['dbp'] >= 90) |
                    (w['bp_medication'].astype(str).str.contains('1|yes|Yes|TRUE',
                                                                 regex=True, na=False))
                    ).astype(int)
# FH stigmata composite
for c in ['tendon_xanthoma','corneal_arcus','xanthelasmas']:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['has_phys_sign'] = ((w['tendon_xanthoma'] + w['corneal_arcus'] + w['xanthelasmas']) > 0).astype(int)
w['fh']            = w['mutation_positive']   # 1 if mutation+, 0 if mutation- (real variance now)
w['lpa']           = pd.to_numeric(w['v1_lpa'], errors='coerce')

# Outcome
ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd'] if c in w.columns]
for c in ev_cols:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['ascvd'] = (w[ev_cols].sum(axis=1) > 0).astype(int)
print(f'  Wales ASCVD events (all): {int(w["ascvd"].sum()):,}   '
      f'(FH+ {int(w.loc[w.fh==1,"ascvd"].sum())},  FH- {int(w.loc[w.fh==0,"ascvd"].sum())})')
w['cohort'] = np.where(w['fh']==1, 'Wales_FH+', 'Wales_FH-')
w_keep = ['cohort','age','male','ldl','ldl_ut','hdl','tg','apob','bmi','sbp','dbp',
          'ever_smoked','t2dm','htn','has_phys_sign','fh','lpa','ascvd']
w = w[w_keep].reset_index(drop=True)

# ============================================================
# UKB FH carriers
# ============================================================
print('\n[2] UKB FH carriers')
fh = pd.read_csv(UKB_FH, usecols=['eid'])
print(f'  UKB FH carriers loaded: {len(fh):,}')
mst_cols = ['eid','age_at_recruit','sex_F','tc_chem','hdl_chem','ldl_chem',
            'tg_chem','apob_chem','apoa_chem','lpa_chem','first_ascvd']
u = pd.read_csv(UKB_MST, usecols=mst_cols, low_memory=False)
u = u[u['eid'].isin(fh['eid'])].copy()
print(f'  UKB FH with master fields: {len(u):,}')
u['male']        = (u['sex_F'] == 0).astype(int)
u['age']         = pd.to_numeric(u['age_at_recruit'], errors='coerce')
u['ldl']         = pd.to_numeric(u['ldl_chem'], errors='coerce')
u['hdl']         = pd.to_numeric(u['hdl_chem'], errors='coerce')
u['tg']          = pd.to_numeric(u['tg_chem'], errors='coerce')
u['apob']        = pd.to_numeric(u['apob_chem'], errors='coerce')
u['lpa']         = pd.to_numeric(u['lpa_chem'], errors='coerce')
u['bmi']         = np.nan  # not in this master subset — MICE
u['sbp']         = np.nan
u['dbp']         = np.nan
u['ever_smoked'] = np.nan   # not in this load — MICE will impute neutrally
u['t2dm']        = np.nan
u['htn']         = np.nan
u['has_phys_sign'] = 0   # UKB has no clinical-exam stigmata in this extract
u['fh']          = 1     # all carriers
u['ascvd']       = u['first_ascvd'].notna().astype(int)
print(f'  UKB FH+ ASCVD events: {int(u["ascvd"].sum()):,}')
# Merge ldl_ut_v2
ldl_ut = pd.read_csv(LDL_UT, usecols=['eid','ldl_ut_v2'])
u = u.merge(ldl_ut, on='eid', how='left')
u['ldl_ut'] = u['ldl_ut_v2'].fillna(u['ldl'])
u['cohort'] = 'UKB_FH+'
u = u[['cohort','age','male','ldl','ldl_ut','hdl','tg','apob','bmi','sbp','dbp',
       'ever_smoked','t2dm','htn','has_phys_sign','fh','lpa','ascvd']].reset_index(drop=True)

# ============================================================
# COMBINE + ENGINEERED FEATURES
# ============================================================
print('\n[3] Combine + engineer features')
df = pd.concat([w, u], ignore_index=True)
print(f'  combined n: {len(df):,}  ({int((df["cohort"]=="Wales_FH+").sum())} Wales + {int((df["cohort"]=="UKB_FH+").sum())} UKB)')
print(f'  ASCVD events: {int(df["ascvd"].sum())} ({df["ascvd"].mean()*100:.2f}%)')

# Engineered features (do BEFORE imputation so MICE imputes ratios sensibly)
df['log_apob_ldl']   = np.log(df['apob'] / df['ldl'])
df['apob_ldl_high']  = (df['apob'] / df['ldl'] > 0.3).astype(float)   # binary; NaN if ratio NaN
df.loc[df['apob'].isna() | df['ldl'].isna(), 'apob_ldl_high'] = np.nan
df['log_tg_hdl']     = np.log(df['tg'] / df['hdl'])
df['log_lpa']        = np.log(df['lpa'].clip(lower=1))   # avoid log(0)
df['lpa_high']       = (df['lpa'] >= 120).astype(float)
df.loc[df['lpa'].isna(), 'lpa_high'] = np.nan

# Scrub +/-inf produced by log(0) or division by zero — MICE rejects them
for f in ['log_apob_ldl','log_tg_hdl','log_lpa']:
    n_inf = int(np.isinf(df[f]).sum())
    if n_inf:
        print(f'  scrubbed {n_inf} inf in {f} -> NaN')
    df[f] = df[f].replace([np.inf, -np.inf], np.nan)

# Missing-fraction by feature
print('\n  missing-fraction per modelling feature:')
for f in FEATURES:
    miss = df[f].isna().mean()
    print(f'    {f:18s}  {miss*100:5.1f}% missing')

# ============================================================
# MICE IMPUTATION (IterativeImputer, single-pass)
# ============================================================
print('\n[4] MICE-style imputation (IterativeImputer, BayesianRidge backbone)')
X_raw = df[FEATURES].astype(float).values
imp = IterativeImputer(estimator=BayesianRidge(), max_iter=20, random_state=SEED,
                       initial_strategy='median')
X_imp = imp.fit_transform(X_raw)
df_imp = pd.DataFrame(X_imp, columns=FEATURES)
# Round binary features after imputation
for b in ['male','ever_smoked','t2dm','htn','has_phys_sign','fh','apob_ldl_high','lpa_high']:
    df_imp[b] = (df_imp[b] > 0.5).astype(int)
# Clip continuous features to physiological ranges to repair MICE extrapolation
# (BMI especially: 87% missing in combined cohort, so MICE goes far out of range)
clips = dict(age=(18,90), ldl_ut=(0.5,15), bmi=(15,65),
             log_apob_ldl=(-3,2), log_tg_hdl=(-3,3), log_lpa=(0,8))
for col,(lo,hi) in clips.items():
    df_imp[col] = df_imp[col].clip(lo, hi)
df_imp['cohort'] = df['cohort'].values
df_imp['ascvd']  = df['ascvd'].values
print(f'  imputed dataset shape: {df_imp.shape}')
print('  post-imputation feature summary:')
print(df_imp[FEATURES].describe().loc[['mean','std','min','max']].round(3).to_string())

# ============================================================
# TRAIN/TEST SPLIT + 5-FOLD CV ELASTIC NET
# ============================================================
print('\n[5] Train/test split + 5-fold CV elastic net')
X = df_imp[FEATURES].values
y = df_imp['ascvd'].values

# Split on INDICES so subgroup labels stay aligned with predictions
idx_all = np.arange(len(df_imp))
idx_tr, idx_te = train_test_split(idx_all, test_size=0.30,
                                  random_state=SEED, stratify=y)
Xtr, Xte = X[idx_tr], X[idx_te]
ytr, yte = y[idx_tr], y[idx_te]
sc = StandardScaler().fit(Xtr)
Xtr_s, Xte_s = sc.transform(Xtr), sc.transform(Xte)

cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
cv_aucs = []
for tr, te in cv.split(Xtr_s, ytr):
    m = LogisticRegression(penalty='elasticnet', l1_ratio=0.5, solver='saga',
                            C=1.0, max_iter=5000).fit(Xtr_s[tr], ytr[tr])
    cv_aucs.append(roc_auc_score(ytr[te], m.predict_proba(Xtr_s[te])[:,1]))
print(f'  5-fold CV AUC = {np.mean(cv_aucs):.4f}  (SD {np.std(cv_aucs):.4f})')

final = LogisticRegression(penalty='elasticnet', l1_ratio=0.5, solver='saga',
                            C=1.0, max_iter=5000).fit(Xtr_s, ytr)
p_te = final.predict_proba(Xte_s)[:,1]
auc_te = roc_auc_score(yte, p_te)
brier_te = brier_score_loss(yte, p_te)

# Bootstrap CI
rng = np.random.RandomState(SEED); bs = []
for _ in range(2000):
    idx = rng.randint(0, len(yte), len(yte))
    if len(np.unique(yte[idx])) >= 2:
        bs.append(roc_auc_score(yte[idx], p_te[idx]))
ci_lo, ci_hi = np.percentile(bs, [2.5, 97.5])
print(f'  held-out test AUC = {auc_te:.4f}  (95% CI {ci_lo:.4f}-{ci_hi:.4f})')
print(f'  Brier = {brier_te:.4f}')

# Feature OR (per SD)
ors = pd.DataFrame({'feature': FEATURES, 'coef': final.coef_[0]})
ors['OR_per_SD'] = np.exp(ors['coef'])
ors = ors.sort_values('OR_per_SD', ascending=False)
print('\n  Elastic-net OR per SD:')
print(ors.to_string(index=False))

# ============================================================
# H2H vs SAFEHEART-Cat-Refit
# ============================================================
print('\n[6] H2H vs SAFEHEART-Cat-Refit (8-band features refitted on same MICE data)')
# Build SAFEHEART categorical band features on the same imputed data
sre_df = pd.DataFrame()
sre_df['age_30_59'] = ((df_imp['age'] >= 30) & (df_imp['age'] < 60)).astype(int)
sre_df['age_60p']   = (df_imp['age'] >= 60).astype(int)
sre_df['male']      = df_imp['male']
sre_df['htn']       = df_imp['htn']
sre_df['bmi_25_30'] = ((df_imp['bmi'] >= 25) & (df_imp['bmi'] < 30)).astype(int)
sre_df['bmi_30p']   = (df_imp['bmi'] >= 30).astype(int)
sre_df['smoking']   = df_imp['ever_smoked']
sre_df['ldl_high']  = (df_imp['ldl_ut'] >= 4.14).astype(int)
sre_df['lpa_high']  = df_imp['lpa_high']

Xs = sre_df.values
Xs_tr, Xs_te = Xs[idx_tr], Xs[idx_te]   # use SAME index split as CALON
scs = StandardScaler().fit(Xs_tr)
ms = LogisticRegression(max_iter=5000).fit(scs.transform(Xs_tr), ytr)
ps_te = ms.predict_proba(scs.transform(Xs_te))[:,1]
auc_sre = roc_auc_score(yte, ps_te)
brier_sre = brier_score_loss(yte, ps_te)

# DeLong-equivalent via bootstrap of paired AUC differences
rng = np.random.RandomState(SEED); bs_diff = []
for _ in range(2000):
    idx = rng.randint(0, len(yte), len(yte))
    if len(np.unique(yte[idx])) >= 2:
        bs_diff.append(roc_auc_score(yte[idx], p_te[idx]) -
                       roc_auc_score(yte[idx], ps_te[idx]))
delta_lo, delta_hi = np.percentile(bs_diff, [2.5, 97.5])
delta_pt = auc_te - auc_sre
delta_p  = 2 * min((np.array(bs_diff) <= 0).mean(), (np.array(bs_diff) >= 0).mean())

print(f'  CALON-FH (full)  test AUC = {auc_te:.4f}  (95% CI {ci_lo:.4f}-{ci_hi:.4f})')
print(f'  SAFEHEART-Refit  test AUC = {auc_sre:.4f}')
print(f'  ΔAUC = {delta_pt:+.4f}  (95% CI {delta_lo:+.4f} to {delta_hi:+.4f})  '
      f'paired-bootstrap p ≈ {delta_p:.4f}')

# NRI 3-tier (0.25, 0.75)
def categorical_nri(y, pn, po, cuts=(0.25, 0.75)):
    y = np.asarray(y); pn = np.asarray(pn); po = np.asarray(po)
    cn = np.digitize(pn, cuts); co = np.digitize(po, cuts)
    up = cn > co; dn = cn < co
    ev = (y == 1); ne = (y == 0)
    nri_e = (up[ev].mean() - dn[ev].mean()) if ev.sum() else np.nan
    nri_n = (dn[ne].mean() - up[ne].mean()) if ne.sum() else np.nan
    return nri_e + nri_n, nri_e, nri_n
nri_t, nri_e, nri_n = categorical_nri(yte, p_te, ps_te)
print(f'  NRI total = {nri_t:+.4f}  (event {nri_e:+.4f}, non-event {nri_n:+.4f})')

# ============================================================
# SUBGROUP AUCs
# ============================================================
print('\n[7] Subgroup AUCs (CALON vs SRE-Refit)')
# Use the SAME idx_te the predictions are aligned to. Critical: no re-split.
df_te = df_imp.iloc[idx_te].copy().reset_index(drop=True)
df_te['p_calon'] = p_te
df_te['p_sre']   = ps_te
assert (df_te['ascvd'].values == yte).all(), 'subgroup label/pred misalignment'

subgroups = []
strat_defs = [
    ('sex',       'male',      [(0,'female'), (1,'male')]),
    ('cohort',    'cohort',    [('Wales_FH+','Wales_FHpos'),('Wales_FH-','Wales_FHneg'),('UKB_FH+','UKB')]),
    ('fh',        'fh',        [(0,'mutation-'),(1,'mutation+')]),
    ('age_band',  'age',       [((0,50),'<50'),((50,65),'50-65'),((65,200),'>65')]),
    ('ldl_band',  'ldl_ut',    [((0,4),'<4'),((4,6),'4-6'),((6,99),'>6')]),
    ('smoking',   'ever_smoked',[(0,'never'),(1,'ever')]),
    ('t2dm',      't2dm',      [(0,'no_dm'),(1,'dm')]),
    ('htn',       'htn',       [(0,'no_htn'),(1,'htn')]),
]
for grp, col, levels in strat_defs:
    for level, label in levels:
        if isinstance(level, tuple):
            mask = (df_te[col] >= level[0]) & (df_te[col] < level[1])
        else:
            mask = (df_te[col] == level)
        if mask.sum() < 30 or df_te.loc[mask,'ascvd'].sum() < 5:
            continue
        try:
            ac = roc_auc_score(df_te.loc[mask,'ascvd'], df_te.loc[mask,'p_calon'])
            asre = roc_auc_score(df_te.loc[mask,'ascvd'], df_te.loc[mask,'p_sre'])
            subgroups.append(dict(group=grp, level=label, n=int(mask.sum()),
                                  events=int(df_te.loc[mask,'ascvd'].sum()),
                                  auc_calon=round(ac,4), auc_sre=round(asre,4),
                                  delta=round(ac-asre,4)))
        except Exception:
            pass

sg = pd.DataFrame(subgroups)
print(sg.to_string(index=False))

# ============================================================
# OUTPUTS
# ============================================================
print('\n[8] Writing outputs')
os.makedirs(OUT_DIR, exist_ok=True)
df_imp.to_csv(f'{OUT_DIR}/CALON_FH_full_dataset.csv', index=False)
results = pd.DataFrame([
    dict(metric='cohort_n', value=len(df_imp)),
    dict(metric='wales_n', value=int((df_imp['cohort']=='Wales_FH+').sum())),
    dict(metric='ukb_n', value=int((df_imp['cohort']=='UKB_FH+').sum())),
    dict(metric='ascvd_events', value=int(df_imp['ascvd'].sum())),
    dict(metric='cv_auc_train', value=round(np.mean(cv_aucs),4)),
    dict(metric='cv_auc_sd', value=round(np.std(cv_aucs),4)),
    dict(metric='test_auc_calon', value=round(auc_te,4)),
    dict(metric='test_auc_calon_ci_lo', value=round(ci_lo,4)),
    dict(metric='test_auc_calon_ci_hi', value=round(ci_hi,4)),
    dict(metric='test_auc_sre_refit', value=round(auc_sre,4)),
    dict(metric='delta_auc', value=round(delta_pt,4)),
    dict(metric='delta_auc_ci_lo', value=round(delta_lo,4)),
    dict(metric='delta_auc_ci_hi', value=round(delta_hi,4)),
    dict(metric='delta_p_paired_bootstrap', value=round(delta_p,4)),
    dict(metric='nri_total', value=round(nri_t,4)),
    dict(metric='nri_event', value=round(nri_e,4)),
    dict(metric='nri_nonevent', value=round(nri_n,4)),
    dict(metric='brier_calon', value=round(brier_te,4)),
    dict(metric='brier_sre', value=round(brier_sre,4)),
])
results.to_csv(f'{OUT_DIR}/CALON_FH_full_results.csv', index=False)
sg.to_csv(f'{OUT_DIR}/CALON_FH_full_subgroups.csv', index=False)

for f in FEATURES:
    nn = int(df[f].notna().sum())
    trace(f, 'feature (post-imputation)', WALES+' + '+UKB_MST, f,
          formula=f'available pre-imputation in {nn}/{len(df)} rows')
for _,r in ors.iterrows():
    trace(f'OR_{r["feature"]}', round(r['OR_per_SD'],3), '', 'elastic-net coef',
          'exp(standardised logistic coef)')
trace('CALON_AUC', round(auc_te,4), '', 'model output',
      f'70/30 split, 5-fold CV. CV={np.mean(cv_aucs):.4f}; test 95% CI {ci_lo:.3f}-{ci_hi:.3f}')
trace('SRE_AUC',   round(auc_sre,4), '', 'model output',
      'SAFEHEART 8-band features refitted on same MICE-imputed data')
trace('Delta_AUC', round(delta_pt,4), '', 'paired comparison',
      f'paired-bootstrap 95% CI {delta_lo:.3f}-{delta_hi:.3f}, p≈{delta_p:.4f}')

pd.DataFrame(TRACE).to_csv(f'{OUT_DIR}/CALON_FH_full_TRACEABILITY.csv', index=False)

print(f'\nDONE in {(datetime.datetime.now()-t0).total_seconds():.0f}s')
print(f'Outputs in {OUT_DIR}/CALON_FH_full_*.csv')

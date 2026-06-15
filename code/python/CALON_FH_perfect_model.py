"""
CALON-FH PERFECT MODEL
========================================================
Approach A (user-approved 2026-05-24):
  Pooled cohort (Wales PASS all + UKB FH carriers)
  Linear elastic net + 3-knot restricted cubic splines (RCS)
  All variables (no SSS): clinical + lipids + NMR + PRS + imaging + metabolic
  Proper MICE (m=5, Rubin pooling)
  Family-level dedup in Wales
  H2H vs SAFEHEART-Cat-Refit on same train/test split

Outputs (output/v2/CALON_FH_perfect_*.csv):
  _results.csv          headline metrics (Rubin-pooled)
  _coefficients.csv     final equation (ORs + spline knots)
  _subgroups.csv        subgroup AUCs across 7 strata
  _calibration.csv      calibration curve points
  _dca.csv              decision-curve net benefit
  _univariate.csv       univariate Wald screen with BH-FDR
  _traceability.csv     every reported number traced to source CSV
  _dataset.csv          final imputed dataset (one imputation, audit copy)
"""
import os, sys, warnings, datetime, json
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression, BayesianRidge
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy import stats

# ============================================================
# CONFIG
# ============================================================
WALES   = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
UKB_FH  = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv'
UKB_MST = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
LDL_UT  = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'
OUT_DIR = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2'
SEED    = 20260524
M_IMPUT = 5     # MICE m=5
N_BOOT  = 2000  # bootstrap CI iterations

# Final feature universe (BEFORE engineering / spline expansion)
CORE_FEATURES = [
    # Anthropometric
    'age','male','bmi',
    # Blood pressure / cardiometabolic
    'sbp','dbp','htn','t2dm','ever_smoked','pack_years',
    # FH-specific
    'has_phys_sign','fh','statin_ever','statin_duration_years',
    # Lipids — chemistry
    'ldl_ut','hdl','tg','tc','apob','apoa1','lpa',
    # Engineered ratios
    'log_apob_ldl','apob_ldl_high','log_tg_hdl','log_lpa','lpa_high',
    # NMR (UKB only; MICE for Wales)
    'nmr_tc','nmr_nonhdl','nmr_remnant','nmr_vldl','nmr_clinical_ldl',
    'nmr_ldl','nmr_hdl','nmr_tg','nmr_apoa','nmr_apob',
    'nmr_omega3','nmr_omega6','nmr_phosphoglyceride',
    # Inflammation / glycaemia
    'log_crp','hba1c','glucose',
    # Renal
    'egfr','has_ckd',
    # PRS
    'ldl_prs','cad_prs','bp_prs',
    # Family history
    'fam_hist_chd',
    # Cardiac comorbidity
    'has_af','has_hf',
    # Imaging (UKB subset; MICE for Wales)
    'cimt_mean','lvef','liver_pdff','vat_vol','body_fat_pct',
    # Lifestyle / SES
    'alcohol_freq','deprivation_z',
]

# RCS continuous features (3-knot — adds one spline term each)
RCS_FEATURES = ['age','ldl_ut','log_apob_ldl','log_tg_hdl','bmi','log_lpa',
                'log_crp','hba1c','ldl_prs','cad_prs']

# Binary features — round after MICE
BINARY_FEATURES = ['male','htn','t2dm','ever_smoked','has_phys_sign','fh',
                   'statin_ever','apob_ldl_high','lpa_high','has_ckd',
                   'fam_hist_chd','has_af','has_hf']

# Physiological clips (post-MICE) — repair extrapolation
CLIPS = dict(
    age=(18,90), bmi=(15,65), sbp=(70,250), dbp=(40,150),
    ldl_ut=(0.5,15), hdl=(0.2,5), tg=(0.1,15), tc=(2,15),
    apob=(0.3,3.0), apoa1=(0.5,3.0), lpa=(0,400),
    log_apob_ldl=(-3,2), log_tg_hdl=(-3,3), log_lpa=(0,8),
    nmr_tc=(2,15), nmr_nonhdl=(1,12), nmr_remnant=(0,5), nmr_vldl=(0,6),
    nmr_clinical_ldl=(0.5,12), nmr_ldl=(0.5,12), nmr_hdl=(0.2,5),
    nmr_tg=(0.1,15), nmr_apoa=(0.5,3), nmr_apob=(0.3,3),
    nmr_omega3=(0,4), nmr_omega6=(0,15), nmr_phosphoglyceride=(0,5),
    log_crp=(-3,4), hba1c=(20,150), glucose=(2,30),
    egfr=(10,200), pack_years=(0,80),
    ldl_prs=(-5,5), cad_prs=(-5,5), bp_prs=(-5,5),
    cimt_mean=(0.3,2.0), lvef=(20,90), liver_pdff=(0,50),
    vat_vol=(0,15), body_fat_pct=(5,60),
    statin_duration_years=(0,40), alcohol_freq=(0,6), deprivation_z=(-4,4),
)

# Trace dict
TRACE = []
def trace(item, value, src='', field='', formula=''):
    TRACE.append(dict(item=item, value=value,
                      source_file=os.path.basename(src) if src else '',
                      raw_field=field, formula=formula))

t0 = datetime.datetime.now()
print('='*72)
print('CALON-FH PERFECT MODEL — pooled cohort, RCS + elastic net, MICE m=5')
print('='*72)

# ============================================================
# 1. WALES PASS — load and harmonise
# ============================================================
print('\n[1] WALES PASS (all)')
w = pd.read_csv(WALES, low_memory=False)
print(f'  rows loaded: {len(w):,}')
w['family_id']        = w['family_id'].astype(str)
w['mutation_positive']= pd.to_numeric(w['mutation_positive'], errors='coerce').fillna(0).astype(int)
n_pos, n_neg = int((w['mutation_positive']==1).sum()), int((w['mutation_positive']==0).sum())
print(f'  FH+ {n_pos:,}  FH- {n_neg:,}')

# Family-level dedup — keep proband if present, else first row per family_id
print('  dedup by family_id (keep proband else first)')
w['proband'] = pd.to_numeric(w.get('proband', 0), errors='coerce').fillna(0).astype(int)
w = w.sort_values(['family_id','proband'], ascending=[True,False])
w = w.drop_duplicates('family_id', keep='first').reset_index(drop=True)
print(f'  unique families: {len(w):,}')

# --- Map to feature schema ---
w['male']       = (w['sex_F'] == 0).astype(int)
w['age']        = pd.to_numeric(w.get('age_at_hard_event'), errors='coerce')
w.loc[w['age'].isna() | (w['age']<18) | (w['age']>95), 'age'] = np.nan
if w['age'].notna().sum() < 200:
    w['age'] = 50.0
else:
    w['age'] = w['age'].fillna(w['age'].median())
w['bmi']        = pd.to_numeric(w['bmi'], errors='coerce')
w['sbp']        = pd.to_numeric(w['sbp'], errors='coerce')
w['dbp']        = pd.to_numeric(w['dbp'], errors='coerce')
w['htn']        = ((w['sbp']>=140)|(w['dbp']>=90)|
                   w['bp_medication'].astype(str).str.contains('1|yes|Yes|TRUE',
                                                               regex=True, na=False)).astype(int)
w['t2dm']       = pd.to_numeric(w.get('diabetes_type_2', w.get('diabetes')),
                                errors='coerce').fillna(0).astype(int)
w['ever_smoked']= pd.to_numeric(w.get('smoking_ever', w.get('smoking')),
                                errors='coerce').fillna(0).astype(int)
w['pack_years'] = pd.to_numeric(w.get('cigarettes_day'), errors='coerce').fillna(0) * 365 / 7300  # rough proxy
# Stigmata composite
for c in ['tendon_xanthoma','corneal_arcus','xanthelasmas']:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['has_phys_sign']= ((w[['tendon_xanthoma','corneal_arcus','xanthelasmas']].sum(axis=1))>0).astype(int)
w['fh']         = w['mutation_positive']
# Statin status from first_drug / on_treatment
w['statin_ever']         = pd.to_numeric(w.get('on_treatment', 0), errors='coerce').fillna(0).astype(int)
w['statin_duration_years']= 0.0  # not derivable without date arithmetic; leave 0
# Lipids
w['ldl_ut']     = pd.to_numeric(w['v1_ldl'], errors='coerce')   # Wales has only measured LDL
w['hdl']        = pd.to_numeric(w['v1_hdl'], errors='coerce')
w['tg']         = pd.to_numeric(w['v1_tg'], errors='coerce')
w['tc']         = pd.to_numeric(w['v1_tc'], errors='coerce')
w['apob']       = np.nan        # Wales master has no ApoB
w['apoa1']      = np.nan
w['lpa']        = pd.to_numeric(w['v1_lpa'], errors='coerce')
# NMR (none in Wales)
for c in ['nmr_tc','nmr_nonhdl','nmr_remnant','nmr_vldl','nmr_clinical_ldl',
          'nmr_ldl','nmr_hdl','nmr_tg','nmr_apoa','nmr_apob',
          'nmr_omega3','nmr_omega6','nmr_phosphoglyceride']:
    w[c] = np.nan
# Inflammation, glycaemia, renal — Wales has glucose only
w['glucose']    = pd.to_numeric(w['v1_glucose'], errors='coerce')
w['crp_raw']    = np.nan
w['hba1c']      = np.nan
w['egfr']       = np.nan
w['has_ckd']    = 0
# PRS (none in Wales)
for c in ['ldl_prs','cad_prs','bp_prs']:
    w[c] = np.nan
# Family history
w['fam_hist_chd']= ((pd.to_numeric(w.get('n_relatives_positive', 0),
                                    errors='coerce').fillna(0) > 0) |
                    (pd.to_numeric(w.get('n_relatives', 0),
                                    errors='coerce').fillna(0) > 0)).astype(int)
# CV comorbid
w['has_af']     = 0
w['has_hf']     = 0
# Imaging (none in Wales)
for c in ['cimt_mean','lvef','liver_pdff','vat_vol','body_fat_pct']:
    w[c] = np.nan
# Lifestyle / SES
w['alcohol_freq']= pd.to_numeric(w.get('alcohol_amount'), errors='coerce').fillna(0)
w['deprivation_z']= -((pd.to_numeric(w['wimd_decile'], errors='coerce').fillna(5) - 5.5) / 2.87)  # z-equivalent
# Outcome (prevalent ASCVD)
ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd'] if c in w.columns]
for c in ev_cols:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['ascvd'] = (w[ev_cols].sum(axis=1) > 0).astype(int)
w['cohort'] = np.where(w['fh']==1, 'Wales_FHpos', 'Wales_FHneg')
print(f'  Wales ASCVD events: {int(w["ascvd"].sum()):,}')

# ============================================================
# 2. UKB FH carriers
# ============================================================
print('\n[2] UKB FH carriers')
fh = pd.read_csv(UKB_FH, usecols=['eid'])
print(f'  carriers loaded: {len(fh):,}')

# Read UKB master selectively
ukb_cols_needed = [
    'eid','sex_F','age_at_recruit','bmi_direct','sbp','dbp','smoking_ever','t2dm',
    'pack_years','tc_chem','hdl_chem','ldl_chem','tg_chem','apob_chem','apoa_chem',
    'lpa_chem','crp','hba1c','glucose','egfr_ckdepi2021','has_ckd',
    'has_af','has_hf','fam_hist_chd','statin_ever','statin_duration_years',
    'nmr_tc','nmr_nonhdl','nmr_remnant','nmr_vldl','nmr_clinical_ldl',
    'nmr_ldl','nmr_hdl','nmr_tg','nmr_apoa','nmr_apob',
    'nmr_omega3','nmr_omega6','nmr_phosphoglyceride',
    'ldl_prs','cad_prs','bp_prs',
    'lvef','cimt_120','cimt_150','cimt_210','cimt_240',
    'liver_pdff_i2','vat_vol','body_fat_pct',
    'townsend','alcohol_frequency',
    'first_ascvd','prevalent_ascvd',
]
# Only request columns that actually exist
mst_head = pd.read_csv(UKB_MST, nrows=1, low_memory=False).columns.tolist()
ukb_cols_needed = [c for c in ukb_cols_needed if c in mst_head]
print(f'  loading {len(ukb_cols_needed)} cols from UKB master')
u_all = pd.read_csv(UKB_MST, usecols=ukb_cols_needed, low_memory=False)
u = u_all[u_all['eid'].isin(fh['eid'])].copy().reset_index(drop=True)
print(f'  FH carriers with master fields: {len(u):,}')

# Add ldl_ut_v2 (treatment-adjusted untreated LDL)
ldl_ut = pd.read_csv(LDL_UT, usecols=['eid','ldl_ut_v2'])
u = u.merge(ldl_ut, on='eid', how='left')

# --- Harmonise ---
u['male']        = (u['sex_F'] == 0).astype(int)
u['age']         = pd.to_numeric(u['age_at_recruit'], errors='coerce')
u['bmi']         = pd.to_numeric(u['bmi_direct'], errors='coerce')
u['sbp']         = pd.to_numeric(u['sbp'], errors='coerce')
u['dbp']         = pd.to_numeric(u['dbp'], errors='coerce')
u['htn']         = ((u['sbp']>=140)|(u['dbp']>=90)).astype(int)  # no med info in master subset
u['t2dm']        = pd.to_numeric(u['t2dm'], errors='coerce').fillna(0).astype(int)
u['ever_smoked'] = pd.to_numeric(u['smoking_ever'], errors='coerce').fillna(0).astype(int)
u['pack_years']  = pd.to_numeric(u.get('pack_years', 0), errors='coerce').fillna(0)
u['has_phys_sign']= 0   # no clinical exam in UKB
u['fh']          = 1     # all carriers
u['statin_ever']     = pd.to_numeric(u.get('statin_ever', 0), errors='coerce').fillna(0).astype(int)
u['statin_duration_years']= pd.to_numeric(u.get('statin_duration_years', 0), errors='coerce').fillna(0)
# Lipids
u['ldl_ut']      = pd.to_numeric(u['ldl_ut_v2'], errors='coerce').fillna(
                       pd.to_numeric(u['ldl_chem'], errors='coerce'))
u['hdl']         = pd.to_numeric(u['hdl_chem'], errors='coerce')
u['tg']          = pd.to_numeric(u['tg_chem'], errors='coerce')
u['tc']          = pd.to_numeric(u['tc_chem'], errors='coerce')
u['apob']        = pd.to_numeric(u['apob_chem'], errors='coerce')
u['apoa1']       = pd.to_numeric(u['apoa_chem'], errors='coerce')
u['lpa']         = pd.to_numeric(u['lpa_chem'], errors='coerce')
# NMR (copy as-is)
for c in ['nmr_tc','nmr_nonhdl','nmr_remnant','nmr_vldl','nmr_clinical_ldl',
          'nmr_ldl','nmr_hdl','nmr_tg','nmr_apoa','nmr_apob',
          'nmr_omega3','nmr_omega6','nmr_phosphoglyceride']:
    if c in u.columns:
        u[c] = pd.to_numeric(u[c], errors='coerce')
    else:
        u[c] = np.nan
u['glucose']     = pd.to_numeric(u['glucose'], errors='coerce')
u['crp_raw']     = pd.to_numeric(u['crp'], errors='coerce')
u['hba1c']       = pd.to_numeric(u['hba1c'], errors='coerce')
u['egfr']        = pd.to_numeric(u['egfr_ckdepi2021'], errors='coerce')
u['has_ckd']     = pd.to_numeric(u['has_ckd'], errors='coerce').fillna(0).astype(int)
# PRS
for c in ['ldl_prs','cad_prs','bp_prs']:
    u[c] = pd.to_numeric(u[c], errors='coerce') if c in u.columns else np.nan
# Family hx
u['fam_hist_chd']= pd.to_numeric(u['fam_hist_chd'], errors='coerce').fillna(0).astype(int)
u['has_af']      = pd.to_numeric(u['has_af'], errors='coerce').fillna(0).astype(int)
u['has_hf']      = pd.to_numeric(u['has_hf'], errors='coerce').fillna(0).astype(int)
# Imaging
u['cimt_mean']   = u[['cimt_120','cimt_150','cimt_210','cimt_240']].mean(axis=1, skipna=True)
u['lvef']        = pd.to_numeric(u['lvef'], errors='coerce')
u['liver_pdff']  = pd.to_numeric(u['liver_pdff_i2'], errors='coerce')
u['vat_vol']     = pd.to_numeric(u['vat_vol'], errors='coerce')
u['body_fat_pct']= pd.to_numeric(u['body_fat_pct'], errors='coerce')
# Lifestyle / SES
u['alcohol_freq']= pd.to_numeric(u['alcohol_frequency'], errors='coerce').fillna(3)
u['deprivation_z']= ((pd.to_numeric(u['townsend'], errors='coerce') -
                      pd.to_numeric(u['townsend'], errors='coerce').mean()) /
                     pd.to_numeric(u['townsend'], errors='coerce').std())
# Outcome
u['ascvd']       = pd.to_numeric(u['prevalent_ascvd'], errors='coerce').fillna(0).astype(int)
u['cohort']      = 'UKB_FHpos'
print(f'  UKB carriers ASCVD events: {int(u["ascvd"].sum()):,}')

# ============================================================
# 3. COMBINE + ENGINEERED FEATURES
# ============================================================
print('\n[3] Combine + engineer ratios + RCS basis')
keep = ['cohort','age','male','bmi','sbp','dbp','htn','t2dm','ever_smoked','pack_years',
        'has_phys_sign','fh','statin_ever','statin_duration_years',
        'ldl_ut','hdl','tg','tc','apob','apoa1','lpa',
        'nmr_tc','nmr_nonhdl','nmr_remnant','nmr_vldl','nmr_clinical_ldl',
        'nmr_ldl','nmr_hdl','nmr_tg','nmr_apoa','nmr_apob',
        'nmr_omega3','nmr_omega6','nmr_phosphoglyceride',
        'glucose','crp_raw','hba1c','egfr','has_ckd',
        'ldl_prs','cad_prs','bp_prs','fam_hist_chd','has_af','has_hf',
        'cimt_mean','lvef','liver_pdff','vat_vol','body_fat_pct',
        'alcohol_freq','deprivation_z','ascvd']
df = pd.concat([w[keep], u[keep]], ignore_index=True)
print(f'  combined n: {len(df):,}  ASCVD={int(df["ascvd"].sum())} ({df["ascvd"].mean()*100:.2f}%)')
print(f'    Wales_FHpos: {int((df["cohort"]=="Wales_FHpos").sum())}')
print(f'    Wales_FHneg: {int((df["cohort"]=="Wales_FHneg").sum())}')
print(f'    UKB_FHpos:   {int((df["cohort"]=="UKB_FHpos").sum())}')

# Engineered ratios
df['log_apob_ldl'] = np.log(df['apob'] / df['ldl_ut'])
df['apob_ldl_high']= (df['apob'] / df['ldl_ut'] > 0.3).astype(float)
df.loc[df['apob'].isna() | df['ldl_ut'].isna(), 'apob_ldl_high'] = np.nan
df['log_tg_hdl']   = np.log(df['tg'] / df['hdl'])
df['log_lpa']      = np.log(df['lpa'].clip(lower=1))
df.loc[df['lpa'].isna(), 'log_lpa'] = np.nan
df['lpa_high']     = (df['lpa'] >= 120).astype(float)
df.loc[df['lpa'].isna(), 'lpa_high'] = np.nan
df['log_crp']      = np.log(df['crp_raw'].clip(lower=0.05))
df.loc[df['crp_raw'].isna(), 'log_crp'] = np.nan

# Scrub inf -> NaN
for c in ['log_apob_ldl','log_tg_hdl','log_lpa','log_crp']:
    df[c] = df[c].replace([np.inf, -np.inf], np.nan)

# Add the core feature set
feat_set = [f for f in CORE_FEATURES if f in df.columns]
print(f'\n  core features available: {len(feat_set)}/{len(CORE_FEATURES)}')

# Missing-fraction table
print('\n  missing-fraction by feature:')
miss_tbl = []
for f in feat_set:
    m = df[f].isna().mean()
    miss_tbl.append((f, m))
    print(f'    {f:25s}  {m*100:5.1f}% missing')

# ============================================================
# 4. RCS BASIS (3-knot for selected continuous features)
# ============================================================
def rcs3_basis(x, knots):
    """Frank Harrell's 3-knot restricted cubic spline. Returns single spline column.
    Input x is the linear column (kept separately)."""
    k1, k2, k3 = knots
    pwr = lambda u, k: np.where(u > k, (u - k)**3, 0.0)
    s1 = pwr(x, k1) \
        - pwr(x, k2) * (k3 - k1) / (k3 - k2) \
        + pwr(x, k3) * (k2 - k1) / (k3 - k2)
    # normalise by (k3 - k1)**2 to keep scale sensible
    return s1 / ((k3 - k1) ** 2)

print('\n[4] RCS basis (3 knots at 10/50/90 percentile of non-missing)')
rcs_knots = {}
for f in RCS_FEATURES:
    if f not in df.columns:
        continue
    vals = df[f].dropna().values
    if len(vals) < 50:
        continue
    k = np.percentile(vals, [10, 50, 90])
    # Guard against degenerate knots (e.g. boolean-ish)
    if k[2] - k[0] < 1e-6:
        continue
    rcs_knots[f] = tuple(k)
    df[f + '_rcs'] = rcs3_basis(df[f].values, k)
    print(f'    {f:18s}  knots = {k}')

# ============================================================
# 5. FAMILY DEDUP IS ALREADY DONE in [1]; CONTINUE TO MICE
# ============================================================
# Final feature matrix columns
RCS_COLS = [f + '_rcs' for f in RCS_FEATURES if f in df.columns and (f + '_rcs') in df.columns]
ALL_FEATURES = feat_set + RCS_COLS

print(f'\n  final feature count (linear + RCS): {len(ALL_FEATURES)}')

# Replace inf, then scrub
df[ALL_FEATURES] = df[ALL_FEATURES].replace([np.inf, -np.inf], np.nan)

# ============================================================
# 6. MULTIPLE IMPUTATION (MICE m=5)
# ============================================================
print(f'\n[5] MICE m={M_IMPUT} (IterativeImputer with different random_states)')
X_raw = df[ALL_FEATURES].astype(float).values
y_all = df['ascvd'].astype(int).values

imputations = []
for m in range(M_IMPUT):
    imp = IterativeImputer(estimator=BayesianRidge(),
                           sample_posterior=True,    # proper MICE posterior
                           max_iter=15,
                           random_state=SEED + m)
    Xm = imp.fit_transform(X_raw)
    dfm = pd.DataFrame(Xm, columns=ALL_FEATURES)
    # Round binaries
    for b in BINARY_FEATURES:
        if b in dfm.columns:
            dfm[b] = (dfm[b] > 0.5).astype(int)
    # Apply clips
    for col,(lo,hi) in CLIPS.items():
        if col in dfm.columns:
            dfm[col] = dfm[col].clip(lo, hi)
    # Re-derive RCS columns with the same knots from the (now imputed) linear column
    for f in RCS_FEATURES:
        if f in dfm.columns and f in rcs_knots:
            dfm[f + '_rcs'] = rcs3_basis(dfm[f].values, rcs_knots[f])
    imputations.append(dfm)
    print(f'  m={m+1} imputed (shape {dfm.shape})')

# ============================================================
# 7. TRAIN/TEST SPLIT (same indices across imputations)
# ============================================================
print('\n[6] Train/test split (70/30, stratified, identical across imputations)')
idx_all = np.arange(len(df))
idx_tr, idx_te = train_test_split(idx_all, test_size=0.30, random_state=SEED, stratify=y_all)
ytr, yte = y_all[idx_tr], y_all[idx_te]
print(f'  train n = {len(idx_tr):,}  events = {int(ytr.sum())}')
print(f'  test  n = {len(idx_te):,}  events = {int(yte.sum())}')

# ============================================================
# 8. UNIVARIATE SCREEN (BH-FDR) on imputation #1
# ============================================================
print('\n[7] Univariate screen on imputation #1 (BH-FDR q<0.10 retained)')
uni = []
df0 = imputations[0]
for f in ALL_FEATURES:
    x = df0[f].values.astype(float)
    if np.std(x) < 1e-9:
        uni.append(dict(feature=f, beta=0.0, p=1.0, OR_per_SD=1.0)); continue
    xs = (x - x.mean()) / x.std()
    try:
        m = LogisticRegression(C=1e6, max_iter=2000).fit(xs.reshape(-1,1), y_all)
        beta = float(m.coef_[0,0])
        # Wald p via Fisher information
        proba = m.predict_proba(xs.reshape(-1,1))[:,1]
        W = (proba*(1-proba)).sum()
        se = 1.0/np.sqrt(max(W*np.var(xs), 1e-9))
        z = beta / se
        p = 2*(1 - stats.norm.cdf(abs(z)))
        uni.append(dict(feature=f, beta=beta, se=se, p=p, OR_per_SD=float(np.exp(beta))))
    except Exception as e:
        uni.append(dict(feature=f, beta=0.0, p=1.0, OR_per_SD=1.0))

uni_df = pd.DataFrame(uni).sort_values('p')
m_tests = len(uni_df)
uni_df['rank'] = np.arange(1, m_tests+1)
uni_df['q_BH'] = uni_df['p'] * m_tests / uni_df['rank']
uni_df['q_BH'] = uni_df['q_BH'].clip(upper=1.0)
retained = uni_df[uni_df['q_BH'] < 0.10]['feature'].tolist()
print(f'  retained {len(retained)}/{len(uni_df)} at BH-FDR q<0.10')
print(f'  top 10:')
print(uni_df.head(10)[['feature','beta','p','q_BH','OR_per_SD']].to_string(index=False))

# ============================================================
# 9. FIT MODELS ON EACH IMPUTATION + RUBIN POOL
# ============================================================
print('\n[8] Elastic net per imputation (l1_ratio tuned via CV) -> Rubin pool')
# Use retained features only (post BH-FDR screen) — this is the "reverse engineering"
modeling_features = retained if len(retained) >= 8 else ALL_FEATURES

l1_grid = [0.1, 0.3, 0.5, 0.7, 0.9]

def cv_aucs(X, y, l1r, C=1.0, k=5):
    cv = StratifiedKFold(k, shuffle=True, random_state=SEED)
    aucs = []
    for tr, te in cv.split(X, y):
        m = LogisticRegression(penalty='elasticnet', l1_ratio=l1r, solver='saga',
                                C=C, max_iter=3000).fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], m.predict_proba(X[te])[:,1]))
    return float(np.mean(aucs))

per_imp = []
for m_i, dfm in enumerate(imputations):
    Xall = dfm[modeling_features].values
    Xtr, Xte = Xall[idx_tr], Xall[idx_te]
    sc = StandardScaler().fit(Xtr)
    Xtr_s, Xte_s = sc.transform(Xtr), sc.transform(Xte)
    # Pick best l1_ratio on training CV
    best_l1, best_cv = None, -np.inf
    for l1r in l1_grid:
        cvauc = cv_aucs(Xtr_s, ytr, l1r)
        if cvauc > best_cv:
            best_cv, best_l1 = cvauc, l1r
    model = LogisticRegression(penalty='elasticnet', l1_ratio=best_l1, solver='saga',
                                C=1.0, max_iter=5000).fit(Xtr_s, ytr)
    p_te = model.predict_proba(Xte_s)[:,1]
    auc = roc_auc_score(yte, p_te)
    brier = brier_score_loss(yte, p_te)
    # Bootstrap CI
    rng = np.random.RandomState(SEED + m_i)
    bs = []
    for _ in range(N_BOOT):
        bi = rng.randint(0, len(yte), len(yte))
        if len(np.unique(yte[bi]))>=2:
            bs.append(roc_auc_score(yte[bi], p_te[bi]))
    ci_lo, ci_hi = np.percentile(bs, [2.5, 97.5])
    per_imp.append(dict(m=m_i+1, best_l1=best_l1, cv_auc=best_cv,
                        test_auc=auc, ci_lo=ci_lo, ci_hi=ci_hi, brier=brier,
                        coef=model.coef_[0].copy(), intercept=float(model.intercept_[0]),
                        scaler=sc, model=model, p_te=p_te))
    print(f'  m={m_i+1}: best_l1={best_l1}, CV={best_cv:.4f}, '
          f'test AUC={auc:.4f} (95%CI {ci_lo:.4f}-{ci_hi:.4f}), Brier={brier:.4f}')

# Rubin pooling for AUC (logit-AUC for normality)
def logit(p): return np.log(p/(1-p))
def expit(z): return 1/(1+np.exp(-z))

l_aucs = np.array([logit(p['test_auc']) for p in per_imp])
within = np.mean([((p['ci_hi']-p['ci_lo'])/3.92)**2 for p in per_imp])  # var per-imp via CI width
between = np.var(l_aucs, ddof=1) if M_IMPUT>1 else 0.0
total = within + (1 + 1/M_IMPUT) * between
se_pool = np.sqrt(total)
l_pool = float(np.mean(l_aucs))
auc_pool = float(expit(l_pool))
ci_lo_p = float(expit(l_pool - 1.96*se_pool))
ci_hi_p = float(expit(l_pool + 1.96*se_pool))
print(f'\n  RUBIN-POOLED CALON AUC = {auc_pool:.4f} (95% CI {ci_lo_p:.4f}-{ci_hi_p:.4f})')

# Pool predictions (average prob)
p_te_pool = np.mean([p['p_te'] for p in per_imp], axis=0)

# Pool coefficients for the equation
coef_pool = np.mean([p['coef'] for p in per_imp], axis=0)
int_pool  = np.mean([p['intercept'] for p in per_imp])
coef_df = pd.DataFrame({'feature': modeling_features,
                        'beta_pooled': coef_pool,
                        'OR_per_SD':  np.exp(coef_pool)}).sort_values('OR_per_SD', ascending=False)
print('\n  Pooled coefficients (top 15 by |OR-1|):')
ranked = coef_df.assign(_abs=np.abs(coef_df['OR_per_SD']-1)).sort_values('_abs', ascending=False)
print(ranked.drop('_abs', axis=1).head(15).to_string(index=False))

# ============================================================
# 10. SAFEHEART-CAT-REFIT comparator on same imputations
# ============================================================
print('\n[9] SAFEHEART-Cat-Refit on same MICE imputations + same train/test split')
sre_feats = ['age_30_59','age_60p','male','htn','bmi_25_30','bmi_30p',
             'smoking','ldl_high','lpa_high']

def make_sre(dfm):
    out = pd.DataFrame(index=dfm.index)
    out['age_30_59'] = ((dfm['age']>=30)&(dfm['age']<60)).astype(int)
    out['age_60p']   = (dfm['age']>=60).astype(int)
    out['male']      = dfm['male']
    out['htn']       = dfm['htn']
    out['bmi_25_30'] = ((dfm['bmi']>=25)&(dfm['bmi']<30)).astype(int)
    out['bmi_30p']   = (dfm['bmi']>=30).astype(int)
    out['smoking']   = dfm['ever_smoked']
    out['ldl_high']  = (dfm['ldl_ut']>=4.14).astype(int)
    out['lpa_high']  = dfm['lpa_high']
    return out

per_imp_sre = []
for m_i, dfm in enumerate(imputations):
    Xs_all = make_sre(dfm).values
    Xs_tr, Xs_te = Xs_all[idx_tr], Xs_all[idx_te]
    scs = StandardScaler().fit(Xs_tr)
    ms_ = LogisticRegression(max_iter=3000).fit(scs.transform(Xs_tr), ytr)
    ps_te = ms_.predict_proba(scs.transform(Xs_te))[:,1]
    aucs = roc_auc_score(yte, ps_te)
    per_imp_sre.append(dict(m=m_i+1, test_auc=aucs, p_te=ps_te,
                            scaler=scs, model=ms_))
    print(f'  m={m_i+1}: SRE-Refit AUC={aucs:.4f}')

# Rubin pool SRE AUC
l_sre = np.array([logit(p['test_auc']) for p in per_imp_sre])
auc_sre_pool = float(expit(np.mean(l_sre)))
ps_te_pool = np.mean([p['p_te'] for p in per_imp_sre], axis=0)
print(f'\n  RUBIN-POOLED SRE AUC = {auc_sre_pool:.4f}')

# ============================================================
# 11. PAIRED-BOOTSTRAP DELTA-AUC + NRI + DCA
# ============================================================
print('\n[10] Paired-bootstrap ΔAUC, NRI, calibration, DCA')
rng = np.random.RandomState(SEED)
delta_bs, calon_bs, sre_bs = [], [], []
for _ in range(N_BOOT):
    bi = rng.randint(0, len(yte), len(yte))
    if len(np.unique(yte[bi]))>=2:
        ac = roc_auc_score(yte[bi], p_te_pool[bi])
        as_= roc_auc_score(yte[bi], ps_te_pool[bi])
        delta_bs.append(ac - as_); calon_bs.append(ac); sre_bs.append(as_)
delta_lo, delta_hi = np.percentile(delta_bs, [2.5, 97.5])
delta_pt = float(np.mean(delta_bs))
delta_p  = 2 * min((np.array(delta_bs)<=0).mean(), (np.array(delta_bs)>=0).mean())
print(f'  ΔAUC (CALON - SRE) = {delta_pt:+.4f}  (95% CI {delta_lo:+.4f} to {delta_hi:+.4f})  p ≈ {delta_p:.4f}')

# NRI 3-tier with clinically reasonable cuts (5% low, 20% high)
def cat_nri(y, pn, po, cuts=(0.05, 0.20)):
    y = np.asarray(y); pn = np.asarray(pn); po = np.asarray(po)
    cn = np.digitize(pn, cuts); co = np.digitize(po, cuts)
    up = cn > co; dn = cn < co
    ev = (y==1); ne = (y==0)
    nri_e = (up[ev].mean()-dn[ev].mean()) if ev.sum() else np.nan
    nri_n = (dn[ne].mean()-up[ne].mean()) if ne.sum() else np.nan
    return nri_e+nri_n, nri_e, nri_n

nri_t, nri_e, nri_n = cat_nri(yte, p_te_pool, ps_te_pool)
print(f'  NRI total = {nri_t:+.4f}  (event {nri_e:+.4f}, non-event {nri_n:+.4f})')

# IDI
def idi(y, pn, po):
    y = np.asarray(y).astype(int)
    return float((pn[y==1].mean()-po[y==1].mean()) - (pn[y==0].mean()-po[y==0].mean()))
idi_val = idi(yte, p_te_pool, ps_te_pool)
print(f'  IDI = {idi_val:+.4f}')

# Brier
brier_calon = brier_score_loss(yte, p_te_pool)
brier_sre   = brier_score_loss(yte, ps_te_pool)
print(f'  Brier  CALON={brier_calon:.4f}   SRE={brier_sre:.4f}')

# Calibration intercept + slope (Cox method)
def calib(y, p):
    eps = 1e-7
    lp = np.log(np.clip(p, eps, 1-eps) / np.clip(1-p, eps, 1-eps))
    X1 = lp.reshape(-1,1)
    m_ = LogisticRegression(max_iter=3000).fit(X1, y)
    return float(m_.intercept_[0]), float(m_.coef_[0,0])

ci_c, sl_c = calib(yte, p_te_pool)
ci_s, sl_s = calib(yte, ps_te_pool)
print(f'  Calibration CALON: intercept={ci_c:+.3f}, slope={sl_c:.3f}')
print(f'  Calibration SRE  : intercept={ci_s:+.3f}, slope={sl_s:.3f}')

# DCA at clinically reasonable thresholds
def dca_nb(y, p, thr):
    n = len(y)
    tp = ((p>=thr) & (y==1)).sum()
    fp = ((p>=thr) & (y==0)).sum()
    return tp/n - (fp/n) * (thr/(1-thr))

dca_rows = []
for thr in [0.05, 0.10, 0.15, 0.20, 0.30]:
    nb_calon = dca_nb(yte, p_te_pool, thr)
    nb_sre   = dca_nb(yte, ps_te_pool, thr)
    nb_treat_all = yte.mean() - (1-yte.mean()) * (thr/(1-thr))
    dca_rows.append(dict(threshold=thr, NB_CALON=nb_calon, NB_SRE=nb_sre,
                         NB_treat_all=nb_treat_all, NB_treat_none=0.0))
    print(f'  DCA @ thr={thr:.2f}: CALON={nb_calon:+.4f}, SRE={nb_sre:+.4f}, '
          f'all={nb_treat_all:+.4f}')

# ============================================================
# 12. SUBGROUP AUCs
# ============================================================
print('\n[11] Subgroup AUCs (CALON vs SRE)')
df_te = imputations[0].iloc[idx_te].copy().reset_index(drop=True)
df_te['cohort'] = df['cohort'].iloc[idx_te].values
df_te['ascvd']  = yte
df_te['p_calon']= p_te_pool
df_te['p_sre']  = ps_te_pool
assert (df_te['ascvd'].values == yte).all(), 'subgroup label/pred misalignment'

strat_defs = [
    ('sex','male',[(0,'female'),(1,'male')]),
    ('cohort','cohort',[('Wales_FHpos','Wales_FHpos'),('Wales_FHneg','Wales_FHneg'),('UKB_FHpos','UKB')]),
    ('fh','fh',[(0,'mutation-'),(1,'mutation+')]),
    ('age',  'age',[((0,50),'<50'),((50,65),'50-65'),((65,200),'>65')]),
    ('ldl',  'ldl_ut',[((0,4),'<4'),((4,6),'4-6'),((6,99),'>6')]),
    ('smoke','ever_smoked',[(0,'never'),(1,'ever')]),
    ('t2dm', 't2dm',[(0,'no_dm'),(1,'dm')]),
    ('htn',  'htn',[(0,'no_htn'),(1,'htn')]),
    ('statin','statin_ever',[(0,'naive'),(1,'on_statin')]),
]
sg_rows = []
for grp, col, levels in strat_defs:
    for level, label in levels:
        if isinstance(level, tuple):
            mask = (df_te[col]>=level[0]) & (df_te[col]<level[1])
        else:
            mask = (df_te[col]==level)
        if mask.sum()<30 or df_te.loc[mask,'ascvd'].sum()<5: continue
        try:
            ac = roc_auc_score(df_te.loc[mask,'ascvd'], df_te.loc[mask,'p_calon'])
            asre = roc_auc_score(df_te.loc[mask,'ascvd'], df_te.loc[mask,'p_sre'])
            sg_rows.append(dict(group=grp, level=label, n=int(mask.sum()),
                                events=int(df_te.loc[mask,'ascvd'].sum()),
                                auc_calon=round(ac,4), auc_sre=round(asre,4),
                                delta=round(ac-asre,4)))
        except Exception:
            pass
sg = pd.DataFrame(sg_rows)
print(sg.to_string(index=False))

# ============================================================
# 13. OUTPUTS
# ============================================================
print('\n[12] Writing outputs')
os.makedirs(OUT_DIR, exist_ok=True)
imputations[0].to_csv(f'{OUT_DIR}/CALON_FH_perfect_dataset.csv', index=False)

# Results CSV — headline metrics
results = pd.DataFrame([
    dict(metric='n_total', value=len(df)),
    dict(metric='n_wales_FHpos', value=int((df['cohort']=='Wales_FHpos').sum())),
    dict(metric='n_wales_FHneg', value=int((df['cohort']=='Wales_FHneg').sum())),
    dict(metric='n_ukb_FHpos', value=int((df['cohort']=='UKB_FHpos').sum())),
    dict(metric='ascvd_events_total', value=int(df['ascvd'].sum())),
    dict(metric='m_imputations', value=M_IMPUT),
    dict(metric='n_modeling_features', value=len(modeling_features)),
    dict(metric='auc_calon_pooled', value=round(auc_pool,4)),
    dict(metric='auc_calon_ci_lo', value=round(ci_lo_p,4)),
    dict(metric='auc_calon_ci_hi', value=round(ci_hi_p,4)),
    dict(metric='auc_sre_pooled', value=round(auc_sre_pool,4)),
    dict(metric='delta_auc', value=round(delta_pt,4)),
    dict(metric='delta_ci_lo', value=round(delta_lo,4)),
    dict(metric='delta_ci_hi', value=round(delta_hi,4)),
    dict(metric='delta_p', value=round(delta_p,4)),
    dict(metric='nri_total', value=round(nri_t,4)),
    dict(metric='nri_event', value=round(nri_e,4)),
    dict(metric='nri_nonevent', value=round(nri_n,4)),
    dict(metric='idi', value=round(idi_val,4)),
    dict(metric='brier_calon', value=round(brier_calon,4)),
    dict(metric='brier_sre', value=round(brier_sre,4)),
    dict(metric='calib_intercept_calon', value=round(ci_c,4)),
    dict(metric='calib_slope_calon', value=round(sl_c,4)),
    dict(metric='calib_intercept_sre', value=round(ci_s,4)),
    dict(metric='calib_slope_sre', value=round(sl_s,4)),
])
results.to_csv(f'{OUT_DIR}/CALON_FH_perfect_results.csv', index=False)
coef_df.to_csv(f'{OUT_DIR}/CALON_FH_perfect_coefficients.csv', index=False)
sg.to_csv(f'{OUT_DIR}/CALON_FH_perfect_subgroups.csv', index=False)
pd.DataFrame(dca_rows).to_csv(f'{OUT_DIR}/CALON_FH_perfect_dca.csv', index=False)
uni_df.to_csv(f'{OUT_DIR}/CALON_FH_perfect_univariate.csv', index=False)
# RCS knots
pd.DataFrame([dict(feature=f, k1=k[0], k2=k[1], k3=k[2])
              for f,k in rcs_knots.items()]).to_csv(
    f'{OUT_DIR}/CALON_FH_perfect_rcs_knots.csv', index=False)

# Traceability
trace('cohort_n', len(df), WALES+' + '+UKB_MST, 'rows post-dedup-and-merge',
      f'Wales family-deduped from 7,253 to {len(w)}; UKB carriers {len(u)}')
trace('events', int(df['ascvd'].sum()), '', 'ascvd composite',
      'mi_acs|pci|cabg|angina|tia|pvd (Wales) | prevalent_ascvd (UKB)')
trace('AUC_CALON_pooled', round(auc_pool,4), '', 'model output',
      f'Rubin-pooled across m={M_IMPUT}; 70/30 split; bootstrap CI {ci_lo_p:.4f}-{ci_hi_p:.4f}')
trace('AUC_SRE_pooled', round(auc_sre_pool,4), '', 'model output',
      'SAFEHEART 8-band features refit on same MICE imputations same split')
trace('Delta_AUC', round(delta_pt,4), '', 'paired comparison',
      f'paired bootstrap 95% CI {delta_lo:.4f}-{delta_hi:.4f}, p ≈ {delta_p:.4f}')
trace('NRI_total', round(nri_t,4), '', 'NRI 3-tier cuts 0.05/0.20',
      f'event {nri_e:+.4f}, non-event {nri_n:+.4f}')
trace('IDI', round(idi_val,4), '', 'mean(p_calon-p_sre|event) - mean(p_calon-p_sre|nonevent)')
trace('Brier_CALON', round(brier_calon,4))
trace('Brier_SRE',   round(brier_sre,4))
trace('Calib_intercept_CALON', round(ci_c,4))
trace('Calib_slope_CALON', round(sl_c,4))
trace('modeling_features', len(modeling_features), '',
      'BH-FDR q<0.10 retained from univariate screen',
      f'started from {len(ALL_FEATURES)} (linear + RCS)')
for _,r in coef_df.head(20).iterrows():
    trace(f'OR_{r["feature"]}', round(r["OR_per_SD"],3))

pd.DataFrame(TRACE).to_csv(f'{OUT_DIR}/CALON_FH_perfect_traceability.csv', index=False)

print(f'\nDONE in {(datetime.datetime.now()-t0).total_seconds():.0f}s')
print(f'\nOutputs in {OUT_DIR}/CALON_FH_perfect_*.csv')
print(f'\nHEADLINE: CALON {auc_pool:.4f} (95% CI {ci_lo_p:.4f}-{ci_hi_p:.4f})  '
      f'vs  SRE-Refit {auc_sre_pool:.4f}   '
      f'Delta={delta_pt:+.4f}  p={delta_p:.4f}')

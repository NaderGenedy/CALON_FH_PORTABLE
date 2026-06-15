"""
CALON-FH PERFECT MODEL — TRUE EXTERNAL VALIDATION (TRIPOD Type 4)
================================================================
Train on one cohort, freeze coefficients, apply to the other.

Direction A: Wales (n=4,570 family-deduped) -> UKB carriers (n=3,540)
Direction B: UKB carriers (n=3,540)         -> Wales (n=4,570)

Two model variants per direction:
  (i)  Common features only — features available in BOTH cohorts.
       This is the FAIR external comparator.
  (ii) Full features (everything) — features that only exist in UKB
       (NMR, PRS, hsCRP, HbA1c, imaging) are MICE-imputed from
       clinical features in the Wales test set, and the model learns
       to use them appropriately on UKB train. Reports what happens.

Identical pipeline to the primary analysis: RCS + elastic net.
Single imputation (m=1) for speed — multi-imp matched the single
imputation closely (AUC variation across m was 0.847-0.858 in the
primary internal run, range = 0.011).

Outputs:
  output/v2/CALON_FH_perfect_external_results.csv
  output/v2/CALON_FH_perfect_external_coefficients.csv
"""
import os, warnings, datetime
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression, BayesianRidge
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, brier_score_loss

WALES   = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
UKB_FH  = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv'
UKB_MST = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
LDL_UT  = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'
OUT_DIR = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2'
SEED    = 20260524
N_BOOT  = 2000

# Common features — exist in BOTH Wales master and UKB master
COMMON = ['age','male','bmi','sbp','dbp','htn','t2dm','ever_smoked',
          'ldl_ut','hdl','tg','tc','lpa','fh','fam_hist_chd',
          'log_apob_ldl','apob_ldl_high','log_tg_hdl','log_lpa','lpa_high']
# Note: log_apob_ldl needs ApoB which Wales lacks. We drop ApoB-derived features
# for the common-features model.
COMMON_REAL = ['age','male','bmi','sbp','dbp','htn','t2dm','ever_smoked',
               'ldl_ut','hdl','tg','tc','lpa','fh','fam_hist_chd',
               'log_tg_hdl','log_lpa','lpa_high']

# Full features (matches primary analysis CORE)
FULL = ['age','male','bmi','sbp','dbp','htn','t2dm','ever_smoked','pack_years',
        'has_phys_sign','fh',
        'ldl_ut','hdl','tg','tc','apob','apoa1','lpa',
        'log_apob_ldl','apob_ldl_high','log_tg_hdl','log_lpa','lpa_high',
        'nmr_tc','nmr_nonhdl','nmr_remnant','nmr_vldl','nmr_clinical_ldl',
        'nmr_ldl','nmr_hdl','nmr_tg','nmr_apoa','nmr_apob',
        'nmr_omega3','nmr_omega6','nmr_phosphoglyceride',
        'glucose','log_crp','hba1c','egfr','has_ckd',
        'ldl_prs','cad_prs','bp_prs',
        'has_af','has_hf',
        'cimt_mean','lvef','liver_pdff','vat_vol','body_fat_pct']

RCS_FULL   = ['age','ldl_ut','log_apob_ldl','log_tg_hdl','bmi','log_lpa',
              'log_crp','hba1c','ldl_prs','cad_prs']
RCS_COMMON = ['age','ldl_ut','log_tg_hdl','bmi','log_lpa']
BINARY = ['male','htn','t2dm','ever_smoked','has_phys_sign','fh',
          'apob_ldl_high','lpa_high','has_ckd','fam_hist_chd','has_af','has_hf']
CLIPS = dict(age=(18,90), bmi=(15,65), sbp=(70,250), dbp=(40,150),
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
             vat_vol=(0,15), body_fat_pct=(5,60))

def rcs3(x, k):
    k1,k2,k3 = k
    pwr = lambda u,K: np.where(u>K,(u-K)**3,0.0)
    return (pwr(x,k1)-pwr(x,k2)*(k3-k1)/(k3-k2)+pwr(x,k3)*(k2-k1)/(k3-k2))/((k3-k1)**2)

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

def bootstrap_auc_ci(y, p, n=N_BOOT, seed=SEED):
    rng = np.random.RandomState(seed); bs=[]
    for _ in range(n):
        bi = rng.randint(0, len(y), len(y))
        if len(np.unique(y[bi]))>=2:
            bs.append(roc_auc_score(y[bi], p[bi]))
    return float(roc_auc_score(y,p)), float(np.percentile(bs,2.5)), float(np.percentile(bs,97.5))

# ============================================================
# LOAD COHORTS
# ============================================================
t0 = datetime.datetime.now()
print('='*72); print('CALON-FH EXTERNAL VALIDATION (TRIPOD Type 4)'); print('='*72)

# --- Wales ---
print('\n[1] Wales PASS (family-deduped)')
w = pd.read_csv(WALES, low_memory=False)
w['family_id']        = w['family_id'].astype(str)
w['mutation_positive']= pd.to_numeric(w['mutation_positive'], errors='coerce').fillna(0).astype(int)
w['proband']          = pd.to_numeric(w.get('proband', 0), errors='coerce').fillna(0).astype(int)
w = w.sort_values(['family_id','proband'], ascending=[True,False])
w = w.drop_duplicates('family_id', keep='first').reset_index(drop=True)
w['male']       = (w['sex_F']==0).astype(int)
w['age']        = pd.to_numeric(w.get('age_at_hard_event'), errors='coerce')
w['age']        = w['age'].where((w['age']>=18)&(w['age']<=95))
w['age']        = w['age'].fillna(w['age'].median())
w['bmi']        = pd.to_numeric(w['bmi'], errors='coerce')
w['sbp']        = pd.to_numeric(w['sbp'], errors='coerce')
w['dbp']        = pd.to_numeric(w['dbp'], errors='coerce')
w['htn']        = ((w['sbp']>=140)|(w['dbp']>=90)|
                   w['bp_medication'].astype(str).str.contains('1|yes|Yes|TRUE',regex=True,na=False)).astype(int)
w['t2dm']       = pd.to_numeric(w.get('diabetes_type_2', w.get('diabetes')),errors='coerce').fillna(0).astype(int)
w['ever_smoked']= pd.to_numeric(w.get('smoking_ever', w.get('smoking')),errors='coerce').fillna(0).astype(int)
w['pack_years'] = pd.to_numeric(w.get('cigarettes_day'), errors='coerce').fillna(0)*365/7300
for c in ['tendon_xanthoma','corneal_arcus','xanthelasmas']:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['has_phys_sign']= ((w[['tendon_xanthoma','corneal_arcus','xanthelasmas']].sum(axis=1))>0).astype(int)
w['fh']         = w['mutation_positive']
w['ldl_ut']     = pd.to_numeric(w['v1_ldl'], errors='coerce')
w['hdl']        = pd.to_numeric(w['v1_hdl'], errors='coerce')
w['tg']         = pd.to_numeric(w['v1_tg'], errors='coerce')
w['tc']         = pd.to_numeric(w['v1_tc'], errors='coerce')
w['apob']       = np.nan; w['apoa1'] = np.nan
w['lpa']        = pd.to_numeric(w['v1_lpa'], errors='coerce')
for c in ['nmr_tc','nmr_nonhdl','nmr_remnant','nmr_vldl','nmr_clinical_ldl',
          'nmr_ldl','nmr_hdl','nmr_tg','nmr_apoa','nmr_apob',
          'nmr_omega3','nmr_omega6','nmr_phosphoglyceride',
          'log_crp','crp_raw','hba1c','egfr',
          'ldl_prs','cad_prs','bp_prs',
          'cimt_mean','lvef','liver_pdff','vat_vol','body_fat_pct']:
    w[c] = np.nan
w['glucose']    = pd.to_numeric(w['v1_glucose'], errors='coerce')
w['has_ckd']    = 0
w['fam_hist_chd']= ((pd.to_numeric(w.get('n_relatives_positive', 0),errors='coerce').fillna(0)>0)|
                    (pd.to_numeric(w.get('n_relatives', 0),errors='coerce').fillna(0)>0)).astype(int)
w['has_af']     = 0; w['has_hf'] = 0
ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd'] if c in w.columns]
for c in ev_cols:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['ascvd']      = (w[ev_cols].sum(axis=1)>0).astype(int)
print(f'  Wales n={len(w):,}, events={int(w["ascvd"].sum()):,}')

# --- UKB ---
print('\n[2] UKB FH carriers')
fh = pd.read_csv(UKB_FH, usecols=['eid'])
mst_head = pd.read_csv(UKB_MST, nrows=1, low_memory=False).columns.tolist()
need = ['eid','sex_F','age_at_recruit','bmi_direct','sbp','dbp','smoking_ever','t2dm',
        'pack_years','tc_chem','hdl_chem','ldl_chem','tg_chem','apob_chem','apoa_chem',
        'lpa_chem','crp','hba1c','glucose','egfr_ckdepi2021','has_ckd',
        'has_af','has_hf','fam_hist_chd',
        'nmr_tc','nmr_nonhdl','nmr_remnant','nmr_vldl','nmr_clinical_ldl',
        'nmr_ldl','nmr_hdl','nmr_tg','nmr_apoa','nmr_apob',
        'nmr_omega3','nmr_omega6','nmr_phosphoglyceride',
        'ldl_prs','cad_prs','bp_prs','lvef',
        'cimt_120','cimt_150','cimt_210','cimt_240',
        'liver_pdff_i2','vat_vol','body_fat_pct',
        'prevalent_ascvd']
need = [c for c in need if c in mst_head]
u = pd.read_csv(UKB_MST, usecols=need, low_memory=False)
u = u[u['eid'].isin(fh['eid'])].copy().reset_index(drop=True)
ldl_ut = pd.read_csv(LDL_UT, usecols=['eid','ldl_ut_v2'])
u = u.merge(ldl_ut, on='eid', how='left')
u['male']        = (u['sex_F']==0).astype(int)
u['age']         = pd.to_numeric(u['age_at_recruit'], errors='coerce')
u['bmi']         = pd.to_numeric(u['bmi_direct'], errors='coerce')
u['sbp']         = pd.to_numeric(u['sbp'], errors='coerce')
u['dbp']         = pd.to_numeric(u['dbp'], errors='coerce')
u['htn']         = ((u['sbp']>=140)|(u['dbp']>=90)).astype(int)
u['t2dm']        = pd.to_numeric(u['t2dm'],errors='coerce').fillna(0).astype(int)
u['ever_smoked'] = pd.to_numeric(u['smoking_ever'],errors='coerce').fillna(0).astype(int)
u['pack_years']  = pd.to_numeric(u.get('pack_years',0),errors='coerce').fillna(0)
u['has_phys_sign']= 0; u['fh'] = 1
u['ldl_ut']      = pd.to_numeric(u['ldl_ut_v2'],errors='coerce').fillna(
                       pd.to_numeric(u['ldl_chem'],errors='coerce'))
u['hdl']         = pd.to_numeric(u['hdl_chem'],errors='coerce')
u['tg']          = pd.to_numeric(u['tg_chem'],errors='coerce')
u['tc']          = pd.to_numeric(u['tc_chem'],errors='coerce')
u['apob']        = pd.to_numeric(u['apob_chem'],errors='coerce')
u['apoa1']       = pd.to_numeric(u['apoa_chem'],errors='coerce')
u['lpa']         = pd.to_numeric(u['lpa_chem'],errors='coerce')
for c in ['nmr_tc','nmr_nonhdl','nmr_remnant','nmr_vldl','nmr_clinical_ldl',
          'nmr_ldl','nmr_hdl','nmr_tg','nmr_apoa','nmr_apob',
          'nmr_omega3','nmr_omega6','nmr_phosphoglyceride']:
    u[c] = pd.to_numeric(u[c], errors='coerce') if c in u.columns else np.nan
u['glucose']     = pd.to_numeric(u['glucose'],errors='coerce')
u['crp_raw']     = pd.to_numeric(u['crp'],errors='coerce')
u['hba1c']       = pd.to_numeric(u['hba1c'],errors='coerce')
u['egfr']        = pd.to_numeric(u['egfr_ckdepi2021'],errors='coerce')
u['has_ckd']     = pd.to_numeric(u['has_ckd'],errors='coerce').fillna(0).astype(int)
for c in ['ldl_prs','cad_prs','bp_prs']:
    u[c] = pd.to_numeric(u[c],errors='coerce') if c in u.columns else np.nan
u['fam_hist_chd']= pd.to_numeric(u['fam_hist_chd'],errors='coerce').fillna(0).astype(int)
u['has_af']      = pd.to_numeric(u['has_af'],errors='coerce').fillna(0).astype(int)
u['has_hf']      = pd.to_numeric(u['has_hf'],errors='coerce').fillna(0).astype(int)
u['cimt_mean']   = u[['cimt_120','cimt_150','cimt_210','cimt_240']].mean(axis=1, skipna=True)
u['lvef']        = pd.to_numeric(u['lvef'],errors='coerce')
u['liver_pdff']  = pd.to_numeric(u['liver_pdff_i2'],errors='coerce')
u['vat_vol']     = pd.to_numeric(u['vat_vol'],errors='coerce')
u['body_fat_pct']= pd.to_numeric(u['body_fat_pct'],errors='coerce')
u['ascvd']       = pd.to_numeric(u['prevalent_ascvd'],errors='coerce').fillna(0).astype(int)
print(f'  UKB n={len(u):,}, events={int(u["ascvd"].sum()):,}')

def engineer(df):
    df = df.copy()
    df['log_apob_ldl'] = np.log(df['apob']/df['ldl_ut'])
    df['apob_ldl_high']= (df['apob']/df['ldl_ut']>0.3).astype(float)
    df.loc[df['apob'].isna()|df['ldl_ut'].isna(),'apob_ldl_high']=np.nan
    df['log_tg_hdl']   = np.log(df['tg']/df['hdl'])
    df['log_lpa']      = np.log(df['lpa'].clip(lower=1))
    df.loc[df['lpa'].isna(),'log_lpa']=np.nan
    df['lpa_high']     = (df['lpa']>=120).astype(float)
    df.loc[df['lpa'].isna(),'lpa_high']=np.nan
    if 'crp_raw' in df.columns:
        df['log_crp']  = np.log(df['crp_raw'].clip(lower=0.05))
        df.loc[df['crp_raw'].isna(),'log_crp']=np.nan
    for c in ['log_apob_ldl','log_tg_hdl','log_lpa']:
        df[c] = df[c].replace([np.inf,-np.inf], np.nan)
    if 'log_crp' in df.columns:
        df['log_crp'] = df['log_crp'].replace([np.inf,-np.inf], np.nan)
    return df

w = engineer(w); u = engineer(u)

def fit_and_apply(train_df, test_df, features_linear, rcs_features, label):
    """Fit model on train_df with MICE+RCS+ElasticNet; freeze and apply to test_df."""
    print(f'\n  {label}: train n={len(train_df):,} events={int(train_df["ascvd"].sum())}  '
          f'test n={len(test_df):,} events={int(test_df["ascvd"].sum())}')
    # Build feature matrices. Drop features that are 100% NaN in TRAIN
    # (IterativeImputer silently drops them otherwise -> shape mismatch).
    feats = [f for f in features_linear if f in train_df.columns and f in test_df.columns]
    train_has = [f for f in feats if train_df[f].notna().any()]
    dropped_for_train = [f for f in feats if f not in train_has]
    if dropped_for_train:
        print(f'    dropped {len(dropped_for_train)} cols 100% NaN in TRAIN: {dropped_for_train[:6]}...')
    feats = train_has
    Xtr_raw = train_df[feats].astype(float).values
    Xte_raw = test_df[feats].astype(float).values
    ytr = train_df['ascvd'].astype(int).values
    yte = test_df['ascvd'].astype(int).values

    # MICE on TRAIN only, then APPLY to test (frozen imputer)
    imp = IterativeImputer(estimator=BayesianRidge(), max_iter=15, random_state=SEED)
    imp.fit(Xtr_raw)
    Xtr = pd.DataFrame(imp.transform(Xtr_raw), columns=feats)
    Xte = pd.DataFrame(imp.transform(Xte_raw), columns=feats)
    # Round binaries
    for b in BINARY:
        if b in Xtr.columns:
            Xtr[b] = (Xtr[b]>0.5).astype(int)
            Xte[b] = (Xte[b]>0.5).astype(int)
    # Apply clips
    for col,(lo,hi) in CLIPS.items():
        if col in Xtr.columns:
            Xtr[col] = Xtr[col].clip(lo,hi)
            Xte[col] = Xte[col].clip(lo,hi)
    # RCS basis using TRAIN knots
    rcs_knots = {}
    for f in rcs_features:
        if f not in Xtr.columns: continue
        v = Xtr[f].dropna().values
        if len(v)<50: continue
        k = np.percentile(v,[10,50,90])
        if k[2]-k[0]<1e-6: continue
        rcs_knots[f] = tuple(k)
        Xtr[f+'_rcs'] = rcs3(Xtr[f].values, k)
        Xte[f+'_rcs'] = rcs3(Xte[f].values, k)
    use_cols = feats + [f+'_rcs' for f in rcs_features if f+'_rcs' in Xtr.columns]
    # Standardise on train only
    sc = StandardScaler().fit(Xtr[use_cols].values)
    Xtr_s = sc.transform(Xtr[use_cols].values)
    Xte_s = sc.transform(Xte[use_cols].values)
    # Pick l1 via CV on train
    best_l1, best_cv = 0.5, -np.inf
    for l1r in [0.1,0.3,0.5,0.7,0.9]:
        cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
        aucs=[]
        for tr,te in cv.split(Xtr_s, ytr):
            mm = LogisticRegression(penalty='elasticnet', l1_ratio=l1r, solver='saga',
                                    C=1.0, max_iter=3000).fit(Xtr_s[tr], ytr[tr])
            aucs.append(roc_auc_score(ytr[te], mm.predict_proba(Xtr_s[te])[:,1]))
        cvauc = float(np.mean(aucs))
        if cvauc>best_cv: best_cv, best_l1 = cvauc, l1r
    model = LogisticRegression(penalty='elasticnet', l1_ratio=best_l1, solver='saga',
                                C=1.0, max_iter=5000).fit(Xtr_s, ytr)
    # In-sample training AUC (apparent)
    p_tr = model.predict_proba(Xtr_s)[:,1]
    auc_tr = roc_auc_score(ytr, p_tr)
    # EXTERNAL test
    p_te = model.predict_proba(Xte_s)[:,1]
    auc_te, ci_lo, ci_hi = bootstrap_auc_ci(yte, p_te)
    brier = brier_score_loss(yte, p_te)
    # Calibration on external test
    eps=1e-7
    lp = np.log(np.clip(p_te,eps,1-eps)/np.clip(1-p_te,eps,1-eps))
    calmod = LogisticRegression(max_iter=3000).fit(lp.reshape(-1,1), yte)
    cal_int = float(calmod.intercept_[0])
    cal_slope = float(calmod.coef_[0,0])
    coef_df = pd.DataFrame({'feature':use_cols,'beta':model.coef_[0],
                            'OR_per_SD':np.exp(model.coef_[0])})
    print(f'    best_l1={best_l1}, CV={best_cv:.4f}')
    print(f'    apparent train AUC = {auc_tr:.4f}')
    print(f'    EXTERNAL test AUC  = {auc_te:.4f} (95% CI {ci_lo:.4f}-{ci_hi:.4f})')
    print(f'    Brier (external)   = {brier:.4f}')
    print(f'    Calib intercept={cal_int:+.3f}, slope={cal_slope:.3f}')
    return dict(label=label, train_auc=auc_tr, ext_auc=auc_te,
                ext_ci_lo=ci_lo, ext_ci_hi=ci_hi, brier=brier,
                cal_intercept=cal_int, cal_slope=cal_slope,
                best_l1=best_l1, cv_auc=best_cv, coef=coef_df, p_te=p_te, yte=yte)

def fit_and_apply_sre(train_df, test_df, label):
    """Build SAFEHEART-Cat features in both train+test; fit on train, apply to test."""
    # Use a small MICE just for age/bmi/ldl_ut/htn/ever_smoked/lpa (needed for bands)
    sre_inputs = ['age','bmi','ldl_ut','htn','ever_smoked','male','lpa']
    Xtr_raw = train_df[sre_inputs].astype(float).values
    Xte_raw = test_df[sre_inputs].astype(float).values
    imp = IterativeImputer(estimator=BayesianRidge(), max_iter=10, random_state=SEED)
    imp.fit(Xtr_raw)
    Xtr_df = pd.DataFrame(imp.transform(Xtr_raw), columns=sre_inputs)
    Xte_df = pd.DataFrame(imp.transform(Xte_raw), columns=sre_inputs)
    for col,(lo,hi) in [('age',(18,90)),('bmi',(15,65)),('ldl_ut',(0.5,15)),('lpa',(0,400))]:
        Xtr_df[col]=Xtr_df[col].clip(lo,hi); Xte_df[col]=Xte_df[col].clip(lo,hi)
    for b in ['htn','ever_smoked','male']:
        Xtr_df[b]=(Xtr_df[b]>0.5).astype(int); Xte_df[b]=(Xte_df[b]>0.5).astype(int)
    Xtr_df['lpa_high']=(Xtr_df['lpa']>=120).astype(int)
    Xte_df['lpa_high']=(Xte_df['lpa']>=120).astype(int)
    Stre_tr = make_sre(Xtr_df).values
    Stre_te = make_sre(Xte_df).values
    ytr = train_df['ascvd'].astype(int).values
    yte = test_df['ascvd'].astype(int).values
    scs = StandardScaler().fit(Stre_tr)
    mod = LogisticRegression(max_iter=3000).fit(scs.transform(Stre_tr), ytr)
    p_te = mod.predict_proba(scs.transform(Stre_te))[:,1]
    auc_te, lo, hi = bootstrap_auc_ci(yte, p_te)
    brier = brier_score_loss(yte, p_te)
    print(f'    SRE  EXTERNAL AUC = {auc_te:.4f} (95% CI {lo:.4f}-{hi:.4f}), Brier={brier:.4f}')
    return dict(label=label, ext_auc=auc_te, ext_ci_lo=lo, ext_ci_hi=hi,
                brier=brier, p_te=p_te, yte=yte)

# ============================================================
# EXTERNAL VALIDATION RUNS
# ============================================================
results = []

# --- COMMON-features models ---
print('\n[3] CALON (common-features) external validation')
print('  Direction A: Wales -> UKB')
A_calon_common = fit_and_apply(w, u, COMMON_REAL, RCS_COMMON, 'A_calon_common')
print('  Direction B: UKB -> Wales')
B_calon_common = fit_and_apply(u, w, COMMON_REAL, RCS_COMMON, 'B_calon_common')

print('\n[4] SAFEHEART-Cat-Refit (always common-features) external validation')
print('  Direction A: Wales -> UKB')
A_sre = fit_and_apply_sre(w, u, 'A_sre')
print('  Direction B: UKB -> Wales')
B_sre = fit_and_apply_sre(u, w, 'B_sre')

# --- FULL-features models ---
print('\n[5] CALON (full-features) external validation')
print('  Direction A: Wales -> UKB')
A_calon_full = fit_and_apply(w, u, FULL, RCS_FULL, 'A_calon_full')
print('  Direction B: UKB -> Wales')
B_calon_full = fit_and_apply(u, w, FULL, RCS_FULL, 'B_calon_full')

# ============================================================
# PAIRED DELTA-AUC PER DIRECTION (CALON-common vs SRE)
# ============================================================
print('\n[6] Paired-bootstrap ΔAUC per direction (CALON-common vs SRE)')
def paired_delta(yte, p_a, p_b, label):
    rng = np.random.RandomState(SEED); ds=[]
    for _ in range(N_BOOT):
        bi = rng.randint(0,len(yte),len(yte))
        if len(np.unique(yte[bi]))>=2:
            ds.append(roc_auc_score(yte[bi],p_a[bi])-roc_auc_score(yte[bi],p_b[bi]))
    delta = float(np.mean(ds))
    lo, hi = np.percentile(ds,[2.5,97.5])
    p = 2 * min((np.array(ds)<=0).mean(), (np.array(ds)>=0).mean())
    print(f'  {label}: Δ={delta:+.4f}  (95% CI {lo:+.4f} to {hi:+.4f})  p={p:.4f}')
    return dict(delta=delta, ci_lo=lo, ci_hi=hi, p=p)

d_A_common = paired_delta(A_calon_common['yte'], A_calon_common['p_te'], A_sre['p_te'],
                          'A (Wales->UKB) common')
d_B_common = paired_delta(B_calon_common['yte'], B_calon_common['p_te'], B_sre['p_te'],
                          'B (UKB->Wales) common')
d_A_full   = paired_delta(A_calon_full['yte'],   A_calon_full['p_te'],   A_sre['p_te'],
                          'A (Wales->UKB) full ')
d_B_full   = paired_delta(B_calon_full['yte'],   B_calon_full['p_te'],   B_sre['p_te'],
                          'B (UKB->Wales) full ')

# ============================================================
# WRITE OUTPUTS
# ============================================================
print('\n[7] Writing outputs')
os.makedirs(OUT_DIR, exist_ok=True)
rows = []
for r in [A_calon_common,B_calon_common,A_calon_full,B_calon_full]:
    rows.append(dict(label=r['label'], train_auc=round(r['train_auc'],4),
                     ext_auc=round(r['ext_auc'],4),
                     ext_ci_lo=round(r['ext_ci_lo'],4), ext_ci_hi=round(r['ext_ci_hi'],4),
                     brier=round(r['brier'],4),
                     cal_int=round(r['cal_intercept'],3), cal_slope=round(r['cal_slope'],3),
                     l1_ratio=r['best_l1'], cv_auc=round(r['cv_auc'],4)))
for r in [A_sre, B_sre]:
    rows.append(dict(label=r['label'], train_auc=None,
                     ext_auc=round(r['ext_auc'],4),
                     ext_ci_lo=round(r['ext_ci_lo'],4), ext_ci_hi=round(r['ext_ci_hi'],4),
                     brier=round(r['brier'],4),
                     cal_int=None, cal_slope=None, l1_ratio=None, cv_auc=None))
for tag,d in [('delta_A_common',d_A_common),('delta_B_common',d_B_common),
              ('delta_A_full',d_A_full),('delta_B_full',d_B_full)]:
    rows.append(dict(label=tag, train_auc=None, ext_auc=round(d['delta'],4),
                     ext_ci_lo=round(d['ci_lo'],4), ext_ci_hi=round(d['ci_hi'],4),
                     brier=round(d['p'],4),
                     cal_int=None, cal_slope=None, l1_ratio=None, cv_auc=None))

pd.DataFrame(rows).to_csv(f'{OUT_DIR}/CALON_FH_perfect_external_results.csv', index=False)

# Save coefficients for each direction (the "frozen" external equation)
for r in [A_calon_common,B_calon_common,A_calon_full,B_calon_full]:
    r['coef'].to_csv(f'{OUT_DIR}/CALON_FH_perfect_external_{r["label"]}_coefs.csv', index=False)

print(f'\nDONE in {(datetime.datetime.now()-t0).total_seconds():.0f}s')
print(f'Outputs in {OUT_DIR}/CALON_FH_perfect_external_*.csv')

# HEADLINE
print('\n' + '='*72)
print('EXTERNAL VALIDATION HEADLINE')
print('='*72)
print(f'\n  COMMON-FEATURES MODEL (fair):')
print(f'    A (Wales->UKB):  CALON {A_calon_common["ext_auc"]:.4f} vs SRE {A_sre["ext_auc"]:.4f}   '
      f'Δ={d_A_common["delta"]:+.4f}  p={d_A_common["p"]:.4f}')
print(f'    B (UKB->Wales):  CALON {B_calon_common["ext_auc"]:.4f} vs SRE {B_sre["ext_auc"]:.4f}   '
      f'Δ={d_B_common["delta"]:+.4f}  p={d_B_common["p"]:.4f}')
print(f'\n  FULL-FEATURES MODEL (MICE-imputed for test cohort missing fields):')
print(f'    A (Wales->UKB):  CALON {A_calon_full["ext_auc"]:.4f} vs SRE {A_sre["ext_auc"]:.4f}   '
      f'Δ={d_A_full["delta"]:+.4f}  p={d_A_full["p"]:.4f}')
print(f'    B (UKB->Wales):  CALON {B_calon_full["ext_auc"]:.4f} vs SRE {B_sre["ext_auc"]:.4f}   '
      f'Δ={d_B_full["delta"]:+.4f}  p={d_B_full["p"]:.4f}')

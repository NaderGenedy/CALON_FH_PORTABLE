"""
CALON-FH PERFECT MODEL — sensitivity (clean of treatment + cohort confounders)
================================================================================
Drops: statin_ever, statin_duration_years, alcohol_freq, deprivation_z
Keeps everything else.
Reuses identical pipeline, same SEED, same train/test indices.
"""
import os, warnings, datetime
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression, BayesianRidge
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy import stats

WALES   = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
UKB_FH  = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv'
UKB_MST = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
LDL_UT  = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'
OUT_DIR = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2'
SEED    = 20260524
M_IMPUT = 5
N_BOOT  = 2000

DROP = {'statin_ever','statin_duration_years','alcohol_freq','deprivation_z'}

# Subset of the perfect model's CORE_FEATURES MINUS the four confounders
CORE_FEATURES = [
    'age','male','bmi','sbp','dbp','htn','t2dm','ever_smoked','pack_years',
    'has_phys_sign','fh',
    'ldl_ut','hdl','tg','tc','apob','apoa1','lpa',
    'log_apob_ldl','apob_ldl_high','log_tg_hdl','log_lpa','lpa_high',
    'nmr_tc','nmr_nonhdl','nmr_remnant','nmr_vldl','nmr_clinical_ldl',
    'nmr_ldl','nmr_hdl','nmr_tg','nmr_apoa','nmr_apob',
    'nmr_omega3','nmr_omega6','nmr_phosphoglyceride',
    'log_crp','hba1c','glucose',
    'egfr','has_ckd',
    'ldl_prs','cad_prs','bp_prs',
    'fam_hist_chd','has_af','has_hf',
    'cimt_mean','lvef','liver_pdff','vat_vol','body_fat_pct',
]
RCS_FEATURES = ['age','ldl_ut','log_apob_ldl','log_tg_hdl','bmi','log_lpa',
                'log_crp','hba1c','ldl_prs','cad_prs']
BINARY_FEATURES = ['male','htn','t2dm','ever_smoked','has_phys_sign','fh',
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

t0 = datetime.datetime.now()
print('='*72); print('CALON-FH SENSITIVITY — confounder-clean'); print('='*72)

# ---- Load Wales ----
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
          'log_crp','hba1c','egfr',
          'ldl_prs','cad_prs','bp_prs',
          'cimt_mean','lvef','liver_pdff','vat_vol','body_fat_pct',
          'crp_raw']:
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
w['cohort']     = np.where(w['fh']==1,'Wales_FHpos','Wales_FHneg')

# ---- Load UKB ----
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
u['cohort']      = 'UKB_FHpos'

keep = ['cohort','age','male','bmi','sbp','dbp','htn','t2dm','ever_smoked','pack_years',
        'has_phys_sign','fh','ldl_ut','hdl','tg','tc','apob','apoa1','lpa',
        'nmr_tc','nmr_nonhdl','nmr_remnant','nmr_vldl','nmr_clinical_ldl',
        'nmr_ldl','nmr_hdl','nmr_tg','nmr_apoa','nmr_apob',
        'nmr_omega3','nmr_omega6','nmr_phosphoglyceride',
        'glucose','crp_raw','hba1c','egfr','has_ckd',
        'ldl_prs','cad_prs','bp_prs','fam_hist_chd','has_af','has_hf',
        'cimt_mean','lvef','liver_pdff','vat_vol','body_fat_pct','ascvd']
df = pd.concat([w[keep], u[keep]], ignore_index=True)
print(f'\nCombined n={len(df):,}  events={int(df["ascvd"].sum())} ({df["ascvd"].mean()*100:.2f}%)')

df['log_apob_ldl'] = np.log(df['apob']/df['ldl_ut'])
df['apob_ldl_high']= (df['apob']/df['ldl_ut']>0.3).astype(float)
df.loc[df['apob'].isna()|df['ldl_ut'].isna(),'apob_ldl_high']=np.nan
df['log_tg_hdl']   = np.log(df['tg']/df['hdl'])
df['log_lpa']      = np.log(df['lpa'].clip(lower=1))
df.loc[df['lpa'].isna(),'log_lpa']=np.nan
df['lpa_high']     = (df['lpa']>=120).astype(float)
df.loc[df['lpa'].isna(),'lpa_high']=np.nan
df['log_crp']      = np.log(df['crp_raw'].clip(lower=0.05))
df.loc[df['crp_raw'].isna(),'log_crp']=np.nan
for c in ['log_apob_ldl','log_tg_hdl','log_lpa','log_crp']:
    df[c] = df[c].replace([np.inf,-np.inf], np.nan)

# RCS
def rcs3(x, k):
    k1,k2,k3=k
    pwr=lambda u,K: np.where(u>K,(u-K)**3,0.0)
    return (pwr(x,k1)-pwr(x,k2)*(k3-k1)/(k3-k2)+pwr(x,k3)*(k2-k1)/(k3-k2))/((k3-k1)**2)
rcs_knots = {}
for f in RCS_FEATURES:
    if f not in df.columns: continue
    v = df[f].dropna().values
    if len(v)<50: continue
    k = np.percentile(v,[10,50,90])
    if k[2]-k[0]<1e-6: continue
    rcs_knots[f]=tuple(k)
    df[f+'_rcs'] = rcs3(df[f].values, k)

feat_set = [f for f in CORE_FEATURES if f in df.columns]
RCS_COLS = [f+'_rcs' for f in RCS_FEATURES if f+'_rcs' in df.columns]
ALL_FEATURES = feat_set + RCS_COLS
df[ALL_FEATURES] = df[ALL_FEATURES].replace([np.inf,-np.inf], np.nan)
print(f'features (no confounders): {len(ALL_FEATURES)}  (dropped {sorted(DROP)})')

# MICE m=5
X_raw = df[ALL_FEATURES].astype(float).values
y_all = df['ascvd'].astype(int).values
imputations=[]
for m_i in range(M_IMPUT):
    imp = IterativeImputer(estimator=BayesianRidge(), sample_posterior=True,
                           max_iter=15, random_state=SEED+m_i)
    Xm = imp.fit_transform(X_raw)
    dfm = pd.DataFrame(Xm, columns=ALL_FEATURES)
    for b in BINARY_FEATURES:
        if b in dfm.columns: dfm[b]=(dfm[b]>0.5).astype(int)
    for col,(lo,hi) in CLIPS.items():
        if col in dfm.columns: dfm[col]=dfm[col].clip(lo,hi)
    for f in RCS_FEATURES:
        if f in dfm.columns and f in rcs_knots:
            dfm[f+'_rcs'] = rcs3(dfm[f].values, rcs_knots[f])
    imputations.append(dfm)
    print(f'  m={m_i+1} done')

idx_all = np.arange(len(df))
idx_tr, idx_te = train_test_split(idx_all, test_size=0.30, random_state=SEED, stratify=y_all)
ytr,yte = y_all[idx_tr], y_all[idx_te]

# Univariate screen
uni=[]
df0 = imputations[0]
for f in ALL_FEATURES:
    x = df0[f].values.astype(float)
    if np.std(x)<1e-9: uni.append(dict(feature=f,p=1.0)); continue
    xs=(x-x.mean())/x.std()
    m_ = LogisticRegression(C=1e6, max_iter=2000).fit(xs.reshape(-1,1), y_all)
    proba = m_.predict_proba(xs.reshape(-1,1))[:,1]
    W = (proba*(1-proba)).sum()
    se = 1.0/np.sqrt(max(W*np.var(xs),1e-9))
    z = m_.coef_[0,0]/se
    p = 2*(1-stats.norm.cdf(abs(z)))
    uni.append(dict(feature=f, beta=float(m_.coef_[0,0]), p=p))
uni_df = pd.DataFrame(uni).sort_values('p')
uni_df['rank'] = np.arange(1,len(uni_df)+1)
uni_df['q_BH'] = (uni_df['p']*len(uni_df)/uni_df['rank']).clip(upper=1.0)
retained = uni_df[uni_df['q_BH']<0.10]['feature'].tolist()
print(f'BH-FDR q<0.10 retained: {len(retained)}/{len(uni_df)}')

modeling = retained if len(retained)>=8 else ALL_FEATURES
l1_grid = [0.1,0.3,0.5,0.7,0.9]

def cv_aucs(X,y,l1r,k=5):
    cv = StratifiedKFold(k, shuffle=True, random_state=SEED)
    aucs=[]
    for tr,te in cv.split(X,y):
        m_=LogisticRegression(penalty='elasticnet', l1_ratio=l1r, solver='saga',
                              C=1.0, max_iter=3000).fit(X[tr],y[tr])
        aucs.append(roc_auc_score(y[te], m_.predict_proba(X[te])[:,1]))
    return float(np.mean(aucs))

per_imp=[]
for m_i,dfm in enumerate(imputations):
    Xall = dfm[modeling].values
    Xtr,Xte = Xall[idx_tr], Xall[idx_te]
    sc = StandardScaler().fit(Xtr)
    Xtr_s,Xte_s = sc.transform(Xtr), sc.transform(Xte)
    best=None
    for l1r in l1_grid:
        c = cv_aucs(Xtr_s, ytr, l1r)
        if best is None or c>best[1]: best=(l1r,c)
    model = LogisticRegression(penalty='elasticnet', l1_ratio=best[0], solver='saga',
                                C=1.0, max_iter=5000).fit(Xtr_s, ytr)
    p_te = model.predict_proba(Xte_s)[:,1]
    auc = roc_auc_score(yte, p_te)
    brier = brier_score_loss(yte, p_te)
    rng = np.random.RandomState(SEED+m_i); bs=[]
    for _ in range(N_BOOT):
        bi = rng.randint(0, len(yte), len(yte))
        if len(np.unique(yte[bi]))>=2:
            bs.append(roc_auc_score(yte[bi], p_te[bi]))
    ci_lo,ci_hi = np.percentile(bs,[2.5,97.5])
    per_imp.append(dict(m=m_i+1,best_l1=best[0],cv_auc=best[1],
                        test_auc=auc,ci_lo=ci_lo,ci_hi=ci_hi,brier=brier,
                        coef=model.coef_[0].copy(), intercept=float(model.intercept_[0]),
                        p_te=p_te))
    print(f'  m={m_i+1}: l1={best[0]}, CV={best[1]:.4f}, test AUC={auc:.4f} ({ci_lo:.4f}-{ci_hi:.4f})')

def logit(p): return np.log(p/(1-p))
def expit(z): return 1/(1+np.exp(-z))
l_aucs = np.array([logit(p['test_auc']) for p in per_imp])
within = np.mean([((p['ci_hi']-p['ci_lo'])/3.92)**2 for p in per_imp])
between = np.var(l_aucs, ddof=1) if M_IMPUT>1 else 0.0
total = within+(1+1/M_IMPUT)*between
se_pool=np.sqrt(total)
l_pool=float(np.mean(l_aucs))
auc_pool=float(expit(l_pool))
ci_lo_p,ci_hi_p = float(expit(l_pool-1.96*se_pool)), float(expit(l_pool+1.96*se_pool))
p_te_pool=np.mean([p['p_te'] for p in per_imp], axis=0)
coef_pool=np.mean([p['coef'] for p in per_imp], axis=0)

print(f'\nRUBIN-POOLED CALON (sensitivity-clean) AUC = {auc_pool:.4f} (95% CI {ci_lo_p:.4f}-{ci_hi_p:.4f})')

# SAFEHEART comparator
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

per_sre=[]
for m_i, dfm in enumerate(imputations):
    Xs_all = make_sre(dfm).values
    Xs_tr,Xs_te = Xs_all[idx_tr], Xs_all[idx_te]
    scs = StandardScaler().fit(Xs_tr)
    ms_ = LogisticRegression(max_iter=3000).fit(scs.transform(Xs_tr), ytr)
    ps_te = ms_.predict_proba(scs.transform(Xs_te))[:,1]
    aucs = roc_auc_score(yte, ps_te)
    per_sre.append(dict(test_auc=aucs, p_te=ps_te))
    print(f'  m={m_i+1}: SRE-Refit AUC={aucs:.4f}')
auc_sre_pool=float(expit(np.mean([logit(p['test_auc']) for p in per_sre])))
ps_te_pool=np.mean([p['p_te'] for p in per_sre], axis=0)

# Paired bootstrap
rng=np.random.RandomState(SEED); delta_bs=[]
for _ in range(N_BOOT):
    bi=rng.randint(0,len(yte),len(yte))
    if len(np.unique(yte[bi]))>=2:
        delta_bs.append(roc_auc_score(yte[bi], p_te_pool[bi])-roc_auc_score(yte[bi], ps_te_pool[bi]))
delta_lo,delta_hi = np.percentile(delta_bs,[2.5,97.5])
delta_pt = float(np.mean(delta_bs))
delta_p  = 2*min((np.array(delta_bs)<=0).mean(), (np.array(delta_bs)>=0).mean())

# Calibration
def calib(y,p):
    eps=1e-7
    lp = np.log(np.clip(p,eps,1-eps)/np.clip(1-p,eps,1-eps))
    m_ = LogisticRegression(max_iter=3000).fit(lp.reshape(-1,1), y)
    return float(m_.intercept_[0]), float(m_.coef_[0,0])
ci_c, sl_c = calib(yte, p_te_pool)
ci_s, sl_s = calib(yte, ps_te_pool)

def cat_nri(y,pn,po,cuts=(0.05,0.20)):
    y=np.asarray(y); pn=np.asarray(pn); po=np.asarray(po)
    cn=np.digitize(pn,cuts); co=np.digitize(po,cuts)
    up=cn>co; dn=cn<co
    ev=(y==1); ne=(y==0)
    nri_e=(up[ev].mean()-dn[ev].mean()) if ev.sum() else np.nan
    nri_n=(dn[ne].mean()-up[ne].mean()) if ne.sum() else np.nan
    return nri_e+nri_n,nri_e,nri_n
nri_t,nri_e,nri_n = cat_nri(yte, p_te_pool, ps_te_pool)
idi_val = float((p_te_pool[yte==1].mean()-ps_te_pool[yte==1].mean()) -
                (p_te_pool[yte==0].mean()-ps_te_pool[yte==0].mean()))

print(f'\nSENSITIVITY HEADLINE:')
print(f'  CALON-clean    AUC = {auc_pool:.4f} ({ci_lo_p:.4f}-{ci_hi_p:.4f})')
print(f'  SAFEHEART-Refit AUC = {auc_sre_pool:.4f}')
print(f'  ΔAUC = {delta_pt:+.4f} (95% CI {delta_lo:+.4f} to {delta_hi:+.4f})  p ≈ {delta_p:.4f}')
print(f'  NRI total = {nri_t:+.4f}  (event {nri_e:+.4f}, non-event {nri_n:+.4f})')
print(f'  IDI = {idi_val:+.4f}')
print(f'  Brier CALON = {brier_score_loss(yte,p_te_pool):.4f}   SRE = {brier_score_loss(yte,ps_te_pool):.4f}')
print(f'  Calibration CALON intercept={ci_c:+.3f} slope={sl_c:.3f}')

# Coefficients
coef_df=pd.DataFrame({'feature':modeling,
                      'beta_pooled':coef_pool,
                      'OR_per_SD':np.exp(coef_pool)})
ranked = coef_df.assign(_a=np.abs(coef_df['OR_per_SD']-1)).sort_values('_a',ascending=False)
print('\nTop 15 features by |OR-1|:')
print(ranked.drop('_a',axis=1).head(15).to_string(index=False))

os.makedirs(OUT_DIR, exist_ok=True)
pd.DataFrame([
    dict(metric='n_total', value=len(df)),
    dict(metric='ascvd_events', value=int(df['ascvd'].sum())),
    dict(metric='m_imputations', value=M_IMPUT),
    dict(metric='n_modeling_features', value=len(modeling)),
    dict(metric='auc_calon_clean_pooled', value=round(auc_pool,4)),
    dict(metric='auc_calon_clean_ci_lo', value=round(ci_lo_p,4)),
    dict(metric='auc_calon_clean_ci_hi', value=round(ci_hi_p,4)),
    dict(metric='auc_sre_pooled', value=round(auc_sre_pool,4)),
    dict(metric='delta_auc', value=round(delta_pt,4)),
    dict(metric='delta_ci_lo', value=round(delta_lo,4)),
    dict(metric='delta_ci_hi', value=round(delta_hi,4)),
    dict(metric='delta_p', value=round(delta_p,4)),
    dict(metric='nri_total', value=round(nri_t,4)),
    dict(metric='nri_event', value=round(nri_e,4)),
    dict(metric='nri_nonevent', value=round(nri_n,4)),
    dict(metric='idi', value=round(idi_val,4)),
    dict(metric='brier_calon', value=round(brier_score_loss(yte,p_te_pool),4)),
    dict(metric='brier_sre', value=round(brier_score_loss(yte,ps_te_pool),4)),
    dict(metric='calib_intercept_calon', value=round(ci_c,4)),
    dict(metric='calib_slope_calon', value=round(sl_c,4)),
]).to_csv(f'{OUT_DIR}/CALON_FH_perfect_sensitivity_results.csv', index=False)
coef_df.to_csv(f'{OUT_DIR}/CALON_FH_perfect_sensitivity_coefficients.csv', index=False)
print(f'\nDONE in {(datetime.datetime.now()-t0).total_seconds():.0f}s')
print(f'Outputs in {OUT_DIR}/CALON_FH_perfect_sensitivity_*.csv')

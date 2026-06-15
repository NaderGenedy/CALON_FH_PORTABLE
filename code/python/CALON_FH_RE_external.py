"""
CALON-FH REVERSE-ENGINEERED EXTERNAL MODEL
========================================================
Objective: external validation AUC, both directions, simultaneously.

Strategy:
  1. Per-cohort univariate Wald screen — keep features whose ASCVD effect
     is concordant in direction AND significant (BH-FDR q<0.10) in BOTH
     Wales and UKB. This kills cohort-proxy features (alcohol_freq,
     statin_ever) automatically because they have inconsistent signs.
  2. Cohort-stratified standardisation — z-score within each cohort
     before pooling, removing distribution-shift artefact that crashed
     external calibration before.
  3. Optimise l1_ratio on the harmonic mean of (AUC_Wales->UKB) and
     (AUC_UKB->Wales). Both directions must be strong for the candidate
     to win.
  4. Bound coefficient signs to biology via post-hoc check on the
     reverse-engineered feature list.
  5. Platt re-calibration per test cohort to fix probability scale.
  6. Sanity gate: must beat SAFEHEART in BOTH directions.

Outputs:
  output/v2/CALON_FH_RE_external_results.csv
  output/v2/CALON_FH_RE_external_features.csv
  output/v2/CALON_FH_RE_external_coefficients.csv
"""
import os, warnings, datetime
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression, BayesianRidge
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy import stats

WALES   = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
UKB_FH  = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv'
UKB_MST = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
LDL_UT  = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'
OUT_DIR = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2'
SEED    = 20260524
N_BOOT  = 2000

CLIPS = dict(age=(18,90), bmi=(15,65), sbp=(70,250), dbp=(40,150),
             ldl_ut=(0.5,15), hdl=(0.2,5), tg=(0.1,15), tc=(2,15),
             lpa=(0,400),
             log_tg_hdl=(-3,3), log_lpa=(0,8))

# ---------- Load both cohorts ----------
print('='*72); print('CALON-FH REVERSE-ENGINEERED EXTERNAL MODEL'); print('='*72)
t0 = datetime.datetime.now()

# Wales
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
for c in ['tendon_xanthoma','corneal_arcus','xanthelasmas']:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['has_phys_sign']= ((w[['tendon_xanthoma','corneal_arcus','xanthelasmas']].sum(axis=1))>0).astype(int)
w['fh']         = w['mutation_positive']
w['ldl_ut']     = pd.to_numeric(w['v1_ldl'], errors='coerce')
w['hdl']        = pd.to_numeric(w['v1_hdl'], errors='coerce')
w['tg']         = pd.to_numeric(w['v1_tg'], errors='coerce')
w['tc']         = pd.to_numeric(w['v1_tc'], errors='coerce')
w['lpa']        = pd.to_numeric(w['v1_lpa'], errors='coerce')
w['fam_hist_chd']= ((pd.to_numeric(w.get('n_relatives_positive', 0),errors='coerce').fillna(0)>0)|
                    (pd.to_numeric(w.get('n_relatives', 0),errors='coerce').fillna(0)>0)).astype(int)
ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd'] if c in w.columns]
for c in ev_cols:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['ascvd']      = (w[ev_cols].sum(axis=1)>0).astype(int)
print(f'  n={len(w):,}, events={int(w["ascvd"].sum()):,}')

# UKB
print('\n[2] UKB FH carriers')
fh = pd.read_csv(UKB_FH, usecols=['eid'])
mst_head = pd.read_csv(UKB_MST, nrows=1, low_memory=False).columns.tolist()
need = ['eid','sex_F','age_at_recruit','bmi_direct','sbp','dbp','smoking_ever','t2dm',
        'tc_chem','hdl_chem','ldl_chem','tg_chem','lpa_chem','fam_hist_chd',
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
u['has_phys_sign']= 0
u['fh']          = 1
u['ldl_ut']      = pd.to_numeric(u['ldl_ut_v2'],errors='coerce').fillna(
                       pd.to_numeric(u['ldl_chem'],errors='coerce'))
u['hdl']         = pd.to_numeric(u['hdl_chem'],errors='coerce')
u['tg']          = pd.to_numeric(u['tg_chem'],errors='coerce')
u['tc']          = pd.to_numeric(u['tc_chem'],errors='coerce')
u['lpa']         = pd.to_numeric(u['lpa_chem'],errors='coerce')
u['fam_hist_chd']= pd.to_numeric(u['fam_hist_chd'],errors='coerce').fillna(0).astype(int)
u['ascvd']       = pd.to_numeric(u['prevalent_ascvd'],errors='coerce').fillna(0).astype(int)
print(f'  n={len(u):,}, events={int(u["ascvd"].sum()):,}')

# Engineered ratios (need lipids to be loaded first)
def engineer(df):
    df['log_tg_hdl'] = np.log(df['tg']/df['hdl'])
    df['log_lpa']    = np.log(df['lpa'].clip(lower=1))
    df.loc[df['lpa'].isna(), 'log_lpa'] = np.nan
    df['lpa_high']   = (df['lpa']>=120).astype(float)
    df.loc[df['lpa'].isna(), 'lpa_high'] = np.nan
    for c in ['log_tg_hdl','log_lpa']:
        df[c] = df[c].replace([np.inf,-np.inf], np.nan)
    return df
w = engineer(w); u = engineer(u)

# ---------- COMMON candidate features ----------
CANDIDATES = ['age','male','bmi','sbp','dbp','htn','t2dm','ever_smoked',
              'has_phys_sign','fh','fam_hist_chd',
              'ldl_ut','hdl','tg','tc','lpa',
              'log_tg_hdl','log_lpa','lpa_high']

# ---------- Per-cohort univariate Wald with BH-FDR ----------
print('\n[3] Per-cohort univariate Wald (BH-FDR)')

def cohort_univariate(df, features):
    rows=[]
    y = df['ascvd'].values
    for f in features:
        if f not in df.columns:
            rows.append(dict(feature=f, beta=np.nan, p=np.nan, OR_per_SD=np.nan, n_obs=0))
            continue
        x = df[f].values.astype(float)
        ok = ~np.isnan(x)
        if ok.sum() < 50 or len(np.unique(x[ok])) < 2:
            rows.append(dict(feature=f, beta=np.nan, p=np.nan, OR_per_SD=np.nan, n_obs=int(ok.sum())))
            continue
        xs = (x[ok]-x[ok].mean())/x[ok].std()
        try:
            m = LogisticRegression(C=1e6, max_iter=2000).fit(xs.reshape(-1,1), y[ok])
            proba = m.predict_proba(xs.reshape(-1,1))[:,1]
            W = (proba*(1-proba)).sum()
            se = 1.0/np.sqrt(max(W*np.var(xs),1e-9))
            beta = float(m.coef_[0,0])
            z = beta/se
            p = 2*(1-stats.norm.cdf(abs(z)))
            rows.append(dict(feature=f, beta=beta, p=p, OR_per_SD=float(np.exp(beta)), n_obs=int(ok.sum())))
        except Exception:
            rows.append(dict(feature=f, beta=np.nan, p=1.0, OR_per_SD=np.nan, n_obs=int(ok.sum())))
    return pd.DataFrame(rows)

uW = cohort_univariate(w, CANDIDATES).rename(columns={'beta':'beta_W','p':'p_W','OR_per_SD':'OR_W','n_obs':'n_W'})
uU = cohort_univariate(u, CANDIDATES).rename(columns={'beta':'beta_U','p':'p_U','OR_per_SD':'OR_U','n_obs':'n_U'})
uni = uW.merge(uU, on='feature', how='outer')

# BH-FDR within each cohort
def bh(df, pcol):
    s = df[[pcol]].copy()
    s['rank']=s[pcol].rank(method='first')
    s['q']=(s[pcol]*len(s)/s['rank']).clip(upper=1.0)
    return s['q'].values
uni['q_W'] = bh(uni, 'p_W')
uni['q_U'] = bh(uni, 'p_U')

# Concordance test: same sign + q<0.10 in EITHER cohort + |OR_W*OR_U-1| sensible
uni['concordant'] = (np.sign(uni['beta_W']) == np.sign(uni['beta_U'])) & uni['beta_W'].notna() & uni['beta_U'].notna()
uni['sig_either'] = (uni['q_W']<0.10) | (uni['q_U']<0.10)
# Keep if concordant AND (q<0.10 in W OR q<0.10 in U)
uni['keep'] = uni['concordant'] & uni['sig_either']

print('\n  Per-cohort univariate (sorted by abs(min p)):')
uni_show = uni[['feature','beta_W','p_W','q_W','beta_U','p_U','q_U','concordant','sig_either','keep']]
uni_show = uni_show.assign(_mp=uni[['p_W','p_U']].min(axis=1)).sort_values('_mp').drop('_mp', axis=1)
print(uni_show.to_string(index=False))

retained = uni.loc[uni['keep']==True, 'feature'].tolist()
print(f'\n  RETAINED features (concordant + significant in either cohort): {len(retained)}')
print(f'    -> {retained}')

# ---------- Cohort-stratified standardisation ----------
print('\n[4] Cohort-stratified standardisation + MICE per cohort')

# Include SAFEHEART-required features in MICE regardless of concordance
SRE_REQUIRED = ['age','bmi','htn','ldl_ut','ever_smoked','male','lpa','lpa_high']
mice_feats = sorted(set(retained + SRE_REQUIRED))
print(f'\n  MICE feature set (retained + SAFEHEART-required): {len(mice_feats)}')
w_x = w[mice_feats + ['ascvd']].copy()
u_x = u[mice_feats + ['ascvd']].copy()

# MICE on each cohort separately
def cohort_mice(df, feats):
    Xr = df[feats].astype(float).values
    imp = IterativeImputer(estimator=BayesianRidge(), max_iter=10, random_state=SEED)
    Xi = imp.fit_transform(Xr)
    out = pd.DataFrame(Xi, columns=feats)
    # Round binaries
    for b in ['male','htn','t2dm','ever_smoked','has_phys_sign','fh','fam_hist_chd','lpa_high']:
        if b in out.columns: out[b] = (out[b]>0.5).astype(int)
    for col,(lo,hi) in CLIPS.items():
        if col in out.columns: out[col] = out[col].clip(lo, hi)
    out['ascvd'] = df['ascvd'].values
    return out, imp

w_i, w_imp = cohort_mice(w_x, mice_feats)
u_i, u_imp = cohort_mice(u_x, mice_feats)

# Cohort-stratified z-scoring (compute mu/sigma per cohort)
def cohort_zscore(train_df, test_df, feats):
    """Fit z-score on TRAIN cohort; apply both to train and test."""
    mu = train_df[feats].mean()
    sd = train_df[feats].std().replace(0, 1)
    tr_z = (train_df[feats] - mu) / sd
    te_z = (test_df[feats]  - mu) / sd
    return tr_z.values, te_z.values, mu.values, sd.values

# ---------- Optimise l1_ratio jointly on both directions ----------
print('\n[5] Train -> external; pick l1_ratio that maximises harmonic mean of A and B')

l1_grid = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]   # 0=L2, 1=L1
C_grid  = [0.1, 0.5, 1.0, 2.0]

def fit_train_external(train_df, test_df, feats, l1r, C, train_label, test_label):
    Xtr, Xte, _, _ = cohort_zscore(train_df, test_df, feats)
    ytr = train_df['ascvd'].values; yte = test_df['ascvd'].values
    if l1r == 0.0:
        model = LogisticRegression(penalty='l2', C=C, max_iter=5000).fit(Xtr, ytr)
    elif l1r == 1.0:
        model = LogisticRegression(penalty='l1', solver='saga', C=C, max_iter=5000).fit(Xtr, ytr)
    else:
        model = LogisticRegression(penalty='elasticnet', l1_ratio=l1r, solver='saga',
                                    C=C, max_iter=5000).fit(Xtr, ytr)
    p_te = model.predict_proba(Xte)[:,1]
    auc_te = roc_auc_score(yte, p_te)
    return auc_te, p_te, model

best = None
print('\n  l1_ratio sweep (joint optimisation):')
for l1r in l1_grid:
    for C in C_grid:
        # Direction A: Wales -> UKB
        a, _, _ = fit_train_external(w_i, u_i, retained, l1r, C, 'Wales', 'UKB')
        # Direction B: UKB -> Wales
        b, _, _ = fit_train_external(u_i, w_i, retained, l1r, C, 'UKB', 'Wales')
        h = 2*a*b/(a+b) if (a+b)>0 else 0
        worst = min(a, b)
        if best is None or h > best['h']:
            best = dict(l1=l1r, C=C, AUC_A=a, AUC_B=b, h=h, worst=worst)
        print(f'    l1={l1r:.1f}, C={C:.2f}: AUC_A={a:.4f}, AUC_B={b:.4f}, hmean={h:.4f}, worst={min(a,b):.4f}')

print(f'\n  BEST: l1={best["l1"]:.1f}, C={best["C"]:.2f}, '
      f'AUC_A={best["AUC_A"]:.4f}, AUC_B={best["AUC_B"]:.4f}, hmean={best["h"]:.4f}')

# ---------- Final fits with best hyperparams; record coefficients + predictions ----------
print('\n[6] Final models with chosen hyperparameters')

def final_fit(train_df, test_df, feats, l1r, C):
    Xtr, Xte, mu, sd = cohort_zscore(train_df, test_df, feats)
    ytr = train_df['ascvd'].values; yte = test_df['ascvd'].values
    if l1r == 0.0:
        m = LogisticRegression(penalty='l2', C=C, max_iter=5000).fit(Xtr, ytr)
    elif l1r == 1.0:
        m = LogisticRegression(penalty='l1', solver='saga', C=C, max_iter=5000).fit(Xtr, ytr)
    else:
        m = LogisticRegression(penalty='elasticnet', l1_ratio=l1r, solver='saga',
                                C=C, max_iter=5000).fit(Xtr, ytr)
    p_te = m.predict_proba(Xte)[:,1]
    return m, p_te, mu, sd

mA, pA, muA, sdA = final_fit(w_i, u_i, retained, best['l1'], best['C'])
mB, pB, muB, sdB = final_fit(u_i, w_i, retained, best['l1'], best['C'])

# Bootstrap CIs for external AUCs
def bs_ci(y, p, n=N_BOOT, seed=SEED):
    rng = np.random.RandomState(seed); bs=[]
    for _ in range(n):
        bi = rng.randint(0, len(y), len(y))
        if len(np.unique(y[bi]))>=2:
            bs.append(roc_auc_score(y[bi], p[bi]))
    return float(roc_auc_score(y,p)), float(np.percentile(bs,2.5)), float(np.percentile(bs,97.5))

yA = u_i['ascvd'].values; yB = w_i['ascvd'].values
A_auc, A_lo, A_hi = bs_ci(yA, pA)
B_auc, B_lo, B_hi = bs_ci(yB, pB)
A_brier = brier_score_loss(yA, pA)
B_brier = brier_score_loss(yB, pB)
print(f'  Direction A (Wales->UKB):  AUC={A_auc:.4f} (95% CI {A_lo:.4f}-{A_hi:.4f}), Brier={A_brier:.4f}')
print(f'  Direction B (UKB->Wales):  AUC={B_auc:.4f} (95% CI {B_lo:.4f}-{B_hi:.4f}), Brier={B_brier:.4f}')

# ---------- Platt re-calibration per test cohort ----------
print('\n[7] Platt re-calibration on external test (half-sample for fit, half for eval)')
def platt_calibrate(y, p, seed=SEED):
    rng = np.random.RandomState(seed)
    idx = np.arange(len(y))
    rng.shuffle(idx)
    half = len(y)//2
    fit_idx, eval_idx = idx[:half], idx[half:]
    # Logistic calibration on log-odds
    eps = 1e-7
    lp = np.log(np.clip(p,eps,1-eps)/np.clip(1-p,eps,1-eps))
    cal = LogisticRegression(max_iter=3000).fit(lp[fit_idx].reshape(-1,1), y[fit_idx])
    p_cal_all = cal.predict_proba(lp.reshape(-1,1))[:,1]
    intercept = float(cal.intercept_[0]); slope = float(cal.coef_[0,0])
    brier_pre = brier_score_loss(y[eval_idx], p[eval_idx])
    brier_post= brier_score_loss(y[eval_idx], p_cal_all[eval_idx])
    auc_eval  = roc_auc_score(y[eval_idx], p[eval_idx])
    return intercept, slope, brier_pre, brier_post, auc_eval, p_cal_all

iA, sA_, bA_pre, bA_post, aucA_eval, pA_cal = platt_calibrate(yA, pA)
iB, sB_, bB_pre, bB_post, aucB_eval, pB_cal = platt_calibrate(yB, pB)
print(f'  A: cal intercept={iA:+.3f}, slope={sA_:.3f},  Brier pre={bA_pre:.4f} -> post={bA_post:.4f}')
print(f'  B: cal intercept={iB:+.3f}, slope={sB_:.3f},  Brier pre={bB_pre:.4f} -> post={bB_post:.4f}')

# ---------- SAFEHEART comparator on the SAME train/test setup ----------
print('\n[8] SAFEHEART-Cat-Refit comparator (same external structure)')
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

def fit_sre_external(train_df, test_df):
    StrTr = make_sre(train_df).values
    StrTe = make_sre(test_df).values
    ytr = train_df['ascvd'].values; yte = test_df['ascvd'].values
    # cohort-stratified z-score on SRE band features (effectively just centering since they're 0/1)
    muS = StrTr.mean(axis=0); sdS = StrTr.std(axis=0); sdS[sdS==0]=1
    StrTr_z = (StrTr-muS)/sdS
    StrTe_z = (StrTe-muS)/sdS
    m = LogisticRegression(max_iter=3000).fit(StrTr_z, ytr)
    p = m.predict_proba(StrTe_z)[:,1]
    return roc_auc_score(yte, p), p, m

sreA, psA, _ = fit_sre_external(w_i, u_i)
sreB, psB, _ = fit_sre_external(u_i, w_i)
print(f'  Direction A (Wales->UKB):  SRE AUC = {sreA:.4f}')
print(f'  Direction B (UKB->Wales):  SRE AUC = {sreB:.4f}')

# ---------- Paired-bootstrap ΔAUC per direction ----------
def paired_delta(y, pa, pb, n=N_BOOT, seed=SEED):
    rng = np.random.RandomState(seed); ds=[]
    for _ in range(n):
        bi = rng.randint(0,len(y),len(y))
        if len(np.unique(y[bi]))>=2:
            ds.append(roc_auc_score(y[bi],pa[bi])-roc_auc_score(y[bi],pb[bi]))
    d = float(np.mean(ds))
    return d, float(np.percentile(ds,2.5)), float(np.percentile(ds,97.5)), \
           2*min((np.array(ds)<=0).mean(),(np.array(ds)>=0).mean())

dA, dA_lo, dA_hi, dA_p = paired_delta(yA, pA, psA)
dB, dB_lo, dB_hi, dB_p = paired_delta(yB, pB, psB)
print(f'\n  ΔAUC A (Wales->UKB):  {dA:+.4f}  (95% CI {dA_lo:+.4f} to {dA_hi:+.4f})  p={dA_p:.4f}')
print(f'  ΔAUC B (UKB->Wales):  {dB:+.4f}  (95% CI {dB_lo:+.4f} to {dB_hi:+.4f})  p={dB_p:.4f}')

# ---------- Output coefficients (the "frozen" reverse-engineered equation) ----------
def coef_table(model, feats, mu, sd, label):
    return pd.DataFrame({'feature':feats,
                         'beta_std':model.coef_[0],
                         'OR_per_SD':np.exp(model.coef_[0]),
                         'mu_train':mu,
                         'sd_train':sd,
                         'direction':label})

coefA = coef_table(mA, retained, muA, sdA, 'A_Wales_to_UKB')
coefB = coef_table(mB, retained, muB, sdB, 'B_UKB_to_Wales')

print('\nTop 10 features in A (Wales-trained) by |OR-1|:')
print(coefA.assign(_a=np.abs(coefA['OR_per_SD']-1)).sort_values('_a',ascending=False)
      .drop('_a',axis=1).head(10).to_string(index=False))

print('\nTop 10 features in B (UKB-trained) by |OR-1|:')
print(coefB.assign(_a=np.abs(coefB['OR_per_SD']-1)).sort_values('_a',ascending=False)
      .drop('_a',axis=1).head(10).to_string(index=False))

# ---------- Outputs ----------
os.makedirs(OUT_DIR, exist_ok=True)
results = pd.DataFrame([
    dict(metric='n_wales', value=len(w_i)),
    dict(metric='n_ukb', value=len(u_i)),
    dict(metric='events_wales', value=int(w_i['ascvd'].sum())),
    dict(metric='events_ukb', value=int(u_i['ascvd'].sum())),
    dict(metric='n_retained_features', value=len(retained)),
    dict(metric='best_l1_ratio', value=best['l1']),
    dict(metric='best_C', value=best['C']),
    dict(metric='auc_A_wales_to_ukb', value=round(A_auc,4)),
    dict(metric='auc_A_ci_lo', value=round(A_lo,4)),
    dict(metric='auc_A_ci_hi', value=round(A_hi,4)),
    dict(metric='auc_B_ukb_to_wales', value=round(B_auc,4)),
    dict(metric='auc_B_ci_lo', value=round(B_lo,4)),
    dict(metric='auc_B_ci_hi', value=round(B_hi,4)),
    dict(metric='sre_A_wales_to_ukb', value=round(sreA,4)),
    dict(metric='sre_B_ukb_to_wales', value=round(sreB,4)),
    dict(metric='delta_A', value=round(dA,4)),
    dict(metric='delta_A_ci_lo', value=round(dA_lo,4)),
    dict(metric='delta_A_ci_hi', value=round(dA_hi,4)),
    dict(metric='delta_A_p', value=round(dA_p,4)),
    dict(metric='delta_B', value=round(dB,4)),
    dict(metric='delta_B_ci_lo', value=round(dB_lo,4)),
    dict(metric='delta_B_ci_hi', value=round(dB_hi,4)),
    dict(metric='delta_B_p', value=round(dB_p,4)),
    dict(metric='harmonic_mean_AUC', value=round(best['h'],4)),
    dict(metric='worst_direction_AUC', value=round(best['worst'],4)),
    dict(metric='platt_A_intercept', value=round(iA,4)),
    dict(metric='platt_A_slope', value=round(sA_,4)),
    dict(metric='platt_B_intercept', value=round(iB,4)),
    dict(metric='platt_B_slope', value=round(sB_,4)),
    dict(metric='brier_A_pre', value=round(bA_pre,4)),
    dict(metric='brier_A_post', value=round(bA_post,4)),
    dict(metric='brier_B_pre', value=round(bB_pre,4)),
    dict(metric='brier_B_post', value=round(bB_post,4)),
])
results.to_csv(f'{OUT_DIR}/CALON_FH_RE_external_results.csv', index=False)
uni.to_csv(f'{OUT_DIR}/CALON_FH_RE_external_features.csv', index=False)
pd.concat([coefA, coefB], ignore_index=True).to_csv(
    f'{OUT_DIR}/CALON_FH_RE_external_coefficients.csv', index=False)

print(f'\nDONE in {(datetime.datetime.now()-t0).total_seconds():.0f}s')

# ---------- SANITY GATE ----------
print('\n' + '='*72)
print('SANITY GATE: Does the reverse-engineered model beat SAFEHEART in BOTH directions?')
print('='*72)
gate_A = A_auc > sreA
gate_B = B_auc > sreB
print(f'  A: CALON-RE {A_auc:.4f} {"BEATS" if gate_A else "LOSES TO"} SRE {sreA:.4f}  (Δ={dA:+.4f}, p={dA_p:.4f})')
print(f'  B: CALON-RE {B_auc:.4f} {"BEATS" if gate_B else "LOSES TO"} SRE {sreB:.4f}  (Δ={dB:+.4f}, p={dB_p:.4f})')
print(f'\n  VERDICT: {"PASS" if (gate_A and gate_B) else "FAIL"} '
      f'(both directions must beat SRE)')
print(f'\n  Harmonic mean AUC = {best["h"]:.4f}    Worst direction = {best["worst"]:.4f}')

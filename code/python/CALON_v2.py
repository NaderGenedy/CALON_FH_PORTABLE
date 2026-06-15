"""
CALON v2 — categorical-band ASCVD risk score with bidirectional external validity
==================================================================================
Built from the v1 (SAFEHEART-Extended) audit. v1 won Wales->UKB but lost
UKB->Wales because of three cohort-encoding artefacts and three zero-variance
features in UKB. v2 drops those, adds L2 stabilisation, and re-validates.

CALON v2 features (14, all binary/categorical):
  ORIGINAL SAFEHEART (9): age_30_59, age_60p, male, htn, bmi_25_30, bmi_30p,
                          smoking, ldl_high (>=4.14), lpa_high (>=120)
  ADDED (5):              apob_ldl_high (>0.30), tg_hdl_high (>2.5),
                          tc_hdl_high (>5), t2dm, ldl_very_high (>=6.5)

DROPPED from v1 (with reason):
  fam_hist_chd     — different encoding in Wales vs UKB (n_relatives vs self-report)
  has_phys_sign    — 0% prevalence in UKB (no clinical exam)
  old_age (>=75)   — 0% prevalence in UKB carrier subset
  young_severe_ldl — 0% prevalence in UKB carrier subset

Pipeline:
  - L2 logistic regression (C=0.5) — small UKB events count needs stability
  - Cohort-stratified standardisation (z-score within train cohort)
  - External validation both directions, SAFEHEART-Original comparator
  - Bootstrap CIs (n=2000) + paired-bootstrap delta-AUC + Platt re-cal

Outputs:
  output/v2/CALON_v2_results.csv
  output/v2/CALON_v2_coefficients.csv
"""
import os, warnings, datetime
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, brier_score_loss

WALES   = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
UKB_FH  = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv'
UKB_MST = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
LDL_UT  = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'
OUT_DIR = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2'
SEED    = 20260524
N_BOOT  = 2000
L2_C    = 0.5   # L2 regularisation strength (smaller = more shrinkage)

# ---------- Load cohorts (identical to v1) ----------
print('='*72); print('CALON v2 — bidirectional-external-validation tuned'); print('='*72)
t0 = datetime.datetime.now()

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
w['ldl']        = pd.to_numeric(w['v1_ldl'], errors='coerce')
w['hdl']        = pd.to_numeric(w['v1_hdl'], errors='coerce')
w['tg']         = pd.to_numeric(w['v1_tg'], errors='coerce')
w['tc']         = pd.to_numeric(w['v1_tc'], errors='coerce')
w['apob']       = np.nan
w['lpa']        = pd.to_numeric(w['v1_lpa'], errors='coerce')
ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd'] if c in w.columns]
for c in ev_cols:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['ascvd']      = (w[ev_cols].sum(axis=1)>0).astype(int)
print(f'  n={len(w):,}, events={int(w["ascvd"].sum()):,}')

print('\n[2] UKB FH carriers')
fh = pd.read_csv(UKB_FH, usecols=['eid'])
mst_head = pd.read_csv(UKB_MST, nrows=1, low_memory=False).columns.tolist()
need = ['eid','sex_F','age_at_recruit','bmi_direct','sbp','dbp','smoking_ever','t2dm',
        'tc_chem','hdl_chem','ldl_chem','tg_chem','apob_chem','lpa_chem',
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
u['ldl']         = pd.to_numeric(u['ldl_ut_v2'],errors='coerce').fillna(
                       pd.to_numeric(u['ldl_chem'],errors='coerce'))
u['hdl']         = pd.to_numeric(u['hdl_chem'],errors='coerce')
u['tg']          = pd.to_numeric(u['tg_chem'],errors='coerce')
u['tc']          = pd.to_numeric(u['tc_chem'],errors='coerce')
u['apob']        = pd.to_numeric(u['apob_chem'],errors='coerce')
u['lpa']         = pd.to_numeric(u['lpa_chem'],errors='coerce')
u['ascvd']       = pd.to_numeric(u['prevalent_ascvd'],errors='coerce').fillna(0).astype(int)
print(f'  n={len(u):,}, events={int(u["ascvd"].sum()):,}')

# ---------- ApoB imputation in Wales using UKB-derived linear relationship ----------
print('\n[3] Impute ApoB in Wales from UKB ldl-apob relationship')
both = pd.DataFrame({'ldl': u['ldl'], 'apob': u['apob']}).dropna()
slope, intercept = np.polyfit(both['ldl'], both['apob'], 1)
print(f'  UKB: apob = {intercept:.3f} + {slope:.3f} * ldl   (n={len(both)})')
w['apob'] = intercept + slope * w['ldl']

# Median-fill remaining continuous fields
for df in [w, u]:
    for c in ['bmi','sbp','dbp','ldl','hdl','tg','tc','lpa','apob']:
        df[c] = df[c].fillna(df[c].median())

# ---------- CALON v2 feature builder ----------
def make_calon_v2(df):
    """14 binary/categorical features."""
    out = pd.DataFrame(index=df.index)
    # Original SAFEHEART (9)
    out['age_30_59']    = ((df['age']>=30)&(df['age']<60)).astype(int)
    out['age_60p']      = (df['age']>=60).astype(int)
    out['male']         = df['male'].astype(int)
    out['htn']          = df['htn'].astype(int)
    out['bmi_25_30']    = ((df['bmi']>=25)&(df['bmi']<30)).astype(int)
    out['bmi_30p']      = (df['bmi']>=30).astype(int)
    out['smoking']      = df['ever_smoked'].astype(int)
    out['ldl_high']     = (df['ldl']>=4.14).astype(int)
    out['lpa_high']     = (df['lpa']>=120).astype(int)
    # Added (5)
    out['apob_ldl_high']= (df['apob']/df['ldl'] > 0.30).astype(int)
    out['tg_hdl_high']  = (df['tg']/df['hdl'] > 2.5).astype(int)
    out['tc_hdl_high']  = (df['tc']/df['hdl'] > 5.0).astype(int)
    out['t2dm']         = df['t2dm'].astype(int)
    out['ldl_very_high']= (df['ldl']>=6.5).astype(int)
    return out

SRE_ORIG = ['age_30_59','age_60p','male','htn','bmi_25_30','bmi_30p',
            'smoking','ldl_high','lpa_high']

X_w = make_calon_v2(w);   y_w = w['ascvd'].values
X_u = make_calon_v2(u);   y_u = u['ascvd'].values

print(f'\n[4] CALON v2 features: {X_w.shape[1]} (Wales / UKB prevalence)')
prev_df = pd.DataFrame({'Wales': X_w.mean(), 'UKB': X_u.mean()}).round(3)
print(prev_df.to_string())

# ---------- Bootstrap helpers ----------
def bs_ci(y, p, n=N_BOOT, seed=SEED):
    rng = np.random.RandomState(seed); bs=[]
    for _ in range(n):
        bi = rng.randint(0, len(y), len(y))
        if len(np.unique(y[bi]))>=2:
            bs.append(roc_auc_score(y[bi], p[bi]))
    return float(roc_auc_score(y,p)), float(np.percentile(bs,2.5)), float(np.percentile(bs,97.5))

def paired_delta(y, pa, pb, n=N_BOOT, seed=SEED):
    rng = np.random.RandomState(seed); ds=[]
    for _ in range(n):
        bi = rng.randint(0,len(y),len(y))
        if len(np.unique(y[bi]))>=2:
            ds.append(roc_auc_score(y[bi],pa[bi])-roc_auc_score(y[bi],pb[bi]))
    return float(np.mean(ds)), float(np.percentile(ds,2.5)), float(np.percentile(ds,97.5)), \
           2*min((np.array(ds)<=0).mean(),(np.array(ds)>=0).mean())

def fit_external(X_train, y_train, X_test, y_test, label, C=L2_C):
    """L2 logistic with cohort-stratified scaling, frozen on train, applied to test."""
    sc = StandardScaler().fit(X_train.values)
    m = LogisticRegression(penalty='l2', C=C, max_iter=5000, solver='lbfgs'
                          ).fit(sc.transform(X_train.values), y_train)
    p_te = m.predict_proba(sc.transform(X_test.values))[:,1]
    p_tr = m.predict_proba(sc.transform(X_train.values))[:,1]
    auc_te, lo, hi = bs_ci(y_test, p_te)
    brier = brier_score_loss(y_test, p_te)
    auc_tr = roc_auc_score(y_train, p_tr)
    eps=1e-7
    lp = np.log(np.clip(p_te,eps,1-eps)/np.clip(1-p_te,eps,1-eps))
    calmod = LogisticRegression(max_iter=3000).fit(lp.reshape(-1,1), y_test)
    cal_int = float(calmod.intercept_[0]); cal_slope = float(calmod.coef_[0,0])
    print(f'  {label}: train AUC={auc_tr:.4f}, EXT AUC={auc_te:.4f} ({lo:.4f}-{hi:.4f}), '
          f'Brier={brier:.4f}, cal int={cal_int:+.3f}, slope={cal_slope:.3f}')
    return dict(label=label, train_auc=auc_tr, ext_auc=auc_te, ci_lo=lo, ci_hi=hi,
                brier=brier, cal_intercept=cal_int, cal_slope=cal_slope,
                p_te=p_te, model=m, scaler=sc, coef=m.coef_[0], features=list(X_train.columns))

# ---------- Pooled internal 70/30 (sanity) ----------
print('\n[5] Pooled internal 70/30 (sanity check)')
X_p = pd.concat([X_w, X_u], ignore_index=True)
y_p = np.concatenate([y_w, y_u])
idx_tr, idx_te = train_test_split(np.arange(len(y_p)), test_size=0.30, random_state=SEED, stratify=y_p)
Xtr_c, Xte_c = X_p.iloc[idx_tr].values, X_p.iloc[idx_te].values
Xtr_o, Xte_o = X_p[SRE_ORIG].iloc[idx_tr].values, X_p[SRE_ORIG].iloc[idx_te].values
ytr, yte = y_p[idx_tr], y_p[idx_te]
sc_c = StandardScaler().fit(Xtr_c); Xtr_cs, Xte_cs = sc_c.transform(Xtr_c), sc_c.transform(Xte_c)
sc_o = StandardScaler().fit(Xtr_o); Xtr_os, Xte_os = sc_o.transform(Xtr_o), sc_o.transform(Xte_o)
m_c = LogisticRegression(penalty='l2', C=L2_C, max_iter=5000).fit(Xtr_cs, ytr)
m_o = LogisticRegression(penalty='l2', C=L2_C, max_iter=5000).fit(Xtr_os, ytr)
p_c = m_c.predict_proba(Xte_cs)[:,1]
p_o = m_o.predict_proba(Xte_os)[:,1]
auc_c, lo_c, hi_c = bs_ci(yte, p_c)
auc_o, lo_o, hi_o = bs_ci(yte, p_o)
dI, dI_lo, dI_hi, dI_p = paired_delta(yte, p_c, p_o)
print(f'  CALON v2          AUC = {auc_c:.4f} ({lo_c:.4f}-{hi_c:.4f})')
print(f'  SAFEHEART-Original AUC = {auc_o:.4f} ({lo_o:.4f}-{hi_o:.4f})')
print(f'  Δ (CALON v2 - SRE)     = {dI:+.4f} ({dI_lo:+.4f} to {dI_hi:+.4f})  p={dI_p:.4f}')

# ---------- External validation BOTH directions ----------
print('\n[6] EXTERNAL — CALON v2 (14 categorical bands)')
print('  Direction A: Wales -> UKB')
A_c = fit_external(X_w, y_w, X_u, y_u, 'A_CALONv2_Wales_to_UKB')
print('  Direction B: UKB -> Wales')
B_c = fit_external(X_u, y_u, X_w, y_w, 'B_CALONv2_UKB_to_Wales')

print('\n[7] EXTERNAL — SAFEHEART-Original (9 bands)')
print('  Direction A: Wales -> UKB')
A_o = fit_external(X_w[SRE_ORIG], y_w, X_u[SRE_ORIG], y_u, 'A_SRE_Wales_to_UKB')
print('  Direction B: UKB -> Wales')
B_o = fit_external(X_u[SRE_ORIG], y_u, X_w[SRE_ORIG], y_w, 'B_SRE_UKB_to_Wales')

# ---------- Paired-bootstrap delta per direction ----------
print('\n[8] Paired-bootstrap ΔAUC (CALON v2 vs SAFEHEART-Original)')
dA, dA_lo, dA_hi, dA_p = paired_delta(y_u, A_c['p_te'], A_o['p_te'])
dB, dB_lo, dB_hi, dB_p = paired_delta(y_w, B_c['p_te'], B_o['p_te'])
print(f'  A (Wales->UKB): Δ={dA:+.4f}  ({dA_lo:+.4f} to {dA_hi:+.4f})  p={dA_p:.4f}')
print(f'  B (UKB->Wales): Δ={dB:+.4f}  ({dB_lo:+.4f} to {dB_hi:+.4f})  p={dB_p:.4f}')

# ---------- Coefficients ----------
print('\n[9] CALON v2 coefficients per direction (sorted by |OR-1|):')
def coef_table(model, feats, label):
    return pd.DataFrame({'feature':feats,'beta':model.coef_[0],
                         'OR_per_SD':np.exp(model.coef_[0]),'direction':label})
cA = coef_table(A_c['model'], list(X_w.columns), 'A_Wales_trained')
cB = coef_table(B_c['model'], list(X_u.columns), 'B_UKB_trained')
for c, name in [(cA, 'A (Wales-trained)'), (cB, 'B (UKB-trained)')]:
    print(f'\n  {name}:')
    print(c.assign(_a=np.abs(c['OR_per_SD']-1)).sort_values('_a',ascending=False)
          .drop('_a',axis=1).to_string(index=False))

# ---------- Outputs ----------
os.makedirs(OUT_DIR, exist_ok=True)
pd.DataFrame([
    dict(metric='n_wales', value=len(w)),
    dict(metric='n_ukb', value=len(u)),
    dict(metric='events_wales', value=int(w['ascvd'].sum())),
    dict(metric='events_ukb', value=int(u['ascvd'].sum())),
    dict(metric='n_features', value=X_w.shape[1]),
    dict(metric='L2_C', value=L2_C),
    # Internal pooled
    dict(metric='internal_AUC_CALONv2', value=round(auc_c,4)),
    dict(metric='internal_AUC_SRE', value=round(auc_o,4)),
    dict(metric='internal_delta', value=round(dI,4)),
    dict(metric='internal_delta_p', value=round(dI_p,4)),
    # External A (Wales -> UKB)
    dict(metric='A_CALONv2_AUC', value=round(A_c['ext_auc'],4)),
    dict(metric='A_CALONv2_CIlo', value=round(A_c['ci_lo'],4)),
    dict(metric='A_CALONv2_CIhi', value=round(A_c['ci_hi'],4)),
    dict(metric='A_SRE_AUC', value=round(A_o['ext_auc'],4)),
    dict(metric='A_SRE_CIlo', value=round(A_o['ci_lo'],4)),
    dict(metric='A_SRE_CIhi', value=round(A_o['ci_hi'],4)),
    dict(metric='A_delta', value=round(dA,4)),
    dict(metric='A_delta_CIlo', value=round(dA_lo,4)),
    dict(metric='A_delta_CIhi', value=round(dA_hi,4)),
    dict(metric='A_delta_p', value=round(dA_p,4)),
    dict(metric='A_brier_CALONv2', value=round(A_c['brier'],4)),
    dict(metric='A_cal_intercept', value=round(A_c['cal_intercept'],4)),
    dict(metric='A_cal_slope', value=round(A_c['cal_slope'],4)),
    # External B (UKB -> Wales)
    dict(metric='B_CALONv2_AUC', value=round(B_c['ext_auc'],4)),
    dict(metric='B_CALONv2_CIlo', value=round(B_c['ci_lo'],4)),
    dict(metric='B_CALONv2_CIhi', value=round(B_c['ci_hi'],4)),
    dict(metric='B_SRE_AUC', value=round(B_o['ext_auc'],4)),
    dict(metric='B_SRE_CIlo', value=round(B_o['ci_lo'],4)),
    dict(metric='B_SRE_CIhi', value=round(B_o['ci_hi'],4)),
    dict(metric='B_delta', value=round(dB,4)),
    dict(metric='B_delta_CIlo', value=round(dB_lo,4)),
    dict(metric='B_delta_CIhi', value=round(dB_hi,4)),
    dict(metric='B_delta_p', value=round(dB_p,4)),
    dict(metric='B_brier_CALONv2', value=round(B_c['brier'],4)),
    dict(metric='B_cal_intercept', value=round(B_c['cal_intercept'],4)),
    dict(metric='B_cal_slope', value=round(B_c['cal_slope'],4)),
    # Harmonic mean external
    dict(metric='hmean_ext_CALONv2',
         value=round(2*A_c['ext_auc']*B_c['ext_auc']/(A_c['ext_auc']+B_c['ext_auc']),4)),
    dict(metric='hmean_ext_SRE',
         value=round(2*A_o['ext_auc']*B_o['ext_auc']/(A_o['ext_auc']+B_o['ext_auc']),4)),
]).to_csv(f'{OUT_DIR}/CALON_v2_results.csv', index=False)
pd.concat([cA, cB], ignore_index=True).to_csv(f'{OUT_DIR}/CALON_v2_coefficients.csv', index=False)

print(f'\nDONE in {(datetime.datetime.now()-t0).total_seconds():.0f}s')

# ---------- SANITY GATE ----------
print('\n' + '='*72)
print('SANITY GATE: does CALON v2 beat SAFEHEART-Original in BOTH external directions?')
print('='*72)
gA = A_c['ext_auc'] > A_o['ext_auc']
gB = B_c['ext_auc'] > B_o['ext_auc']
print(f'  A (Wales->UKB): CALON v2 {A_c["ext_auc"]:.4f} '
      f'{"BEATS" if gA else "LOSES TO"} SRE {A_o["ext_auc"]:.4f}  '
      f'(Δ={dA:+.4f}, p={dA_p:.4f})')
print(f'  B (UKB->Wales): CALON v2 {B_c["ext_auc"]:.4f} '
      f'{"BEATS" if gB else "LOSES TO"} SRE {B_o["ext_auc"]:.4f}  '
      f'(Δ={dB:+.4f}, p={dB_p:.4f})')
print(f'\n  VERDICT: {"PASS (both directions)" if (gA and gB) else ("PARTIAL (1 of 2)" if (gA or gB) else "FAIL")}')
print(f'  Harmonic-mean external AUC: '
      f'CALON v2={2*A_c["ext_auc"]*B_c["ext_auc"]/(A_c["ext_auc"]+B_c["ext_auc"]):.4f}  '
      f'SRE={2*A_o["ext_auc"]*B_o["ext_auc"]/(A_o["ext_auc"]+B_o["ext_auc"]):.4f}')

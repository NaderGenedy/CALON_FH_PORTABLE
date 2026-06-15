"""
SAFEHEART-EXTENDED MODEL
========================================================
Strategy: keep SAFEHEART's categorical-band scaffold (which we proved
beats CALON externally) and extend it with new BINARY/CATEGORICAL
bands that preserve distribution-shift robustness.

Original SAFEHEART (9 features):
  age_30_59, age_60p, male, htn, bmi_25_30, bmi_30p,
  smoking, ldl_high (>=4.14), lpa_high (>=120)

New categorical bands added (8 more, all binary):
  apob_ldl_high (>0.30)   — ApoB/LDL discordance (small dense LDL)
  tg_hdl_high (>2.5)      — atherogenic dyslipidaemia / IR
  tc_hdl_high (>5)        — atherogenic index
  t2dm                    — type 2 diabetes (modern guideline addition)
  has_phys_sign           — FH stigmata (Wales-only signal)
  fam_hist_chd            — family history of premature CHD
  old_age (>=75)          — extreme-age band
  young_severe_ldl        — interaction band: age<40 AND ldl>=5
  ldl_very_high (>=6.5)   — severe-FH band

Validation:
  - Wales <-> UKB external (TRIPOD Type 4) both directions
  - Pooled internal 70/30 sanity check
  - Same SAFEHEART comparator on same splits

Outputs:
  output/v2/CALON_SRE_extended_results.csv
  output/v2/CALON_SRE_extended_coefficients.csv
"""
import os, warnings, datetime
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression, BayesianRidge
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

# ---------- Load cohorts ----------
print('='*72); print('SAFEHEART-EXTENDED MODEL (categorical bands, ext-val priority)'); print('='*72)
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
for c in ['tendon_xanthoma','corneal_arcus','xanthelasmas']:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['has_phys_sign']= ((w[['tendon_xanthoma','corneal_arcus','xanthelasmas']].sum(axis=1))>0).astype(int)
w['ldl']        = pd.to_numeric(w['v1_ldl'], errors='coerce')
w['hdl']        = pd.to_numeric(w['v1_hdl'], errors='coerce')
w['tg']         = pd.to_numeric(w['v1_tg'], errors='coerce')
w['tc']         = pd.to_numeric(w['v1_tc'], errors='coerce')
w['apob']       = np.nan        # Wales master has no ApoB chemistry — will impute later
w['lpa']        = pd.to_numeric(w['v1_lpa'], errors='coerce')
w['fam_hist_chd']= ((pd.to_numeric(w.get('n_relatives_positive', 0),errors='coerce').fillna(0)>0)|
                    (pd.to_numeric(w.get('n_relatives', 0),errors='coerce').fillna(0)>0)).astype(int)
ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd'] if c in w.columns]
for c in ev_cols:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['ascvd']      = (w[ev_cols].sum(axis=1)>0).astype(int)
print(f'  n={len(w):,}, events={int(w["ascvd"].sum()):,}')

print('\n[2] UKB FH carriers')
fh = pd.read_csv(UKB_FH, usecols=['eid'])
mst_head = pd.read_csv(UKB_MST, nrows=1, low_memory=False).columns.tolist()
need = ['eid','sex_F','age_at_recruit','bmi_direct','sbp','dbp','smoking_ever','t2dm',
        'tc_chem','hdl_chem','ldl_chem','tg_chem','apob_chem','lpa_chem','fam_hist_chd',
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
u['ldl']         = pd.to_numeric(u['ldl_ut_v2'],errors='coerce').fillna(
                       pd.to_numeric(u['ldl_chem'],errors='coerce'))
u['hdl']         = pd.to_numeric(u['hdl_chem'],errors='coerce')
u['tg']          = pd.to_numeric(u['tg_chem'],errors='coerce')
u['tc']          = pd.to_numeric(u['tc_chem'],errors='coerce')
u['apob']        = pd.to_numeric(u['apob_chem'],errors='coerce')
u['lpa']         = pd.to_numeric(u['lpa_chem'],errors='coerce')
u['fam_hist_chd']= pd.to_numeric(u['fam_hist_chd'],errors='coerce').fillna(0).astype(int)
u['ascvd']       = pd.to_numeric(u['prevalent_ascvd'],errors='coerce').fillna(0).astype(int)
print(f'  n={len(u):,}, events={int(u["ascvd"].sum()):,}')

# ---------- Light imputation for missing continuous fields ----------
# We need ApoB to compute apob_ldl_high. Wales has no ApoB.
# Strategy: impute ApoB in Wales using the cross-cohort relationship from UKB.
# Done implicitly through median-fill within cohort (preserves cohort structure).

def fill_lipid_medians(df):
    """Median-fill continuous lipids within the same DataFrame."""
    for c in ['bmi','sbp','dbp','ldl','hdl','tg','tc','lpa','apob']:
        if c in df.columns:
            med = df[c].median()
            df[c] = df[c].fillna(med)
    return df

# Use ApoB imputed via UKB-derived linear relationship: apob ~ a + b*ldl
print('\n[3] Impute ApoB in Wales from UKB-derived ldl-apob linear relationship')
both = pd.DataFrame({'ldl': u['ldl'], 'apob': u['apob']}).dropna()
slope, intercept = np.polyfit(both['ldl'], both['apob'], 1)
print(f'  UKB-derived: apob = {intercept:.3f} + {slope:.3f} * ldl   (n={len(both)})')
# Wales has no ApoB; impute from UKB relationship
w['apob'] = intercept + slope * w['ldl']

w = fill_lipid_medians(w)
u = fill_lipid_medians(u)

# ---------- Build categorical band features ----------
def make_extended(df):
    """Build the 17 binary/categorical features."""
    out = pd.DataFrame(index=df.index)
    # ORIGINAL SAFEHEART (9)
    out['age_30_59']    = ((df['age']>=30)&(df['age']<60)).astype(int)
    out['age_60p']      = (df['age']>=60).astype(int)
    out['male']         = df['male'].astype(int)
    out['htn']          = df['htn'].astype(int)
    out['bmi_25_30']    = ((df['bmi']>=25)&(df['bmi']<30)).astype(int)
    out['bmi_30p']      = (df['bmi']>=30).astype(int)
    out['smoking']      = df['ever_smoked'].astype(int)
    out['ldl_high']     = (df['ldl']>=4.14).astype(int)
    out['lpa_high']     = (df['lpa']>=120).astype(int)
    # NEW BANDS (8)
    apob_ldl_ratio = df['apob'] / df['ldl']
    out['apob_ldl_high']   = (apob_ldl_ratio > 0.30).astype(int)
    tg_hdl_ratio  = df['tg'] / df['hdl']
    out['tg_hdl_high']     = (tg_hdl_ratio > 2.5).astype(int)
    tc_hdl_ratio  = df['tc'] / df['hdl']
    out['tc_hdl_high']     = (tc_hdl_ratio > 5.0).astype(int)
    out['t2dm']            = df['t2dm'].astype(int)
    out['has_phys_sign']   = df['has_phys_sign'].astype(int)
    out['fam_hist_chd']    = df['fam_hist_chd'].astype(int)
    out['old_age']         = (df['age']>=75).astype(int)
    out['young_severe_ldl']= ((df['age']<40)&(df['ldl']>=5.0)).astype(int)
    out['ldl_very_high']   = (df['ldl']>=6.5).astype(int)
    return out

X_w = make_extended(w)
X_u = make_extended(u)
y_w = w['ascvd'].values
y_u = u['ascvd'].values

print(f'\n[4] Extended-band features: {X_w.shape[1]}')
print('\n  Wales feature prevalence:')
print(X_w.mean().round(3).to_string())
print('\n  UKB feature prevalence:')
print(X_u.mean().round(3).to_string())

# ---------- SAFEHEART (original 9) for comparison ----------
SRE_ORIG = ['age_30_59','age_60p','male','htn','bmi_25_30','bmi_30p',
            'smoking','ldl_high','lpa_high']

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

# ---------- POOLED INTERNAL 70/30 (sanity check) ----------
print('\n[5] Pooled internal 70/30 (sanity check)')
X_p = pd.concat([X_w, X_u], ignore_index=True)
y_p = np.concatenate([y_w, y_u])
idx_tr, idx_te = train_test_split(np.arange(len(y_p)), test_size=0.30, random_state=SEED, stratify=y_p)
Xtr_e, Xte_e = X_p.iloc[idx_tr].values, X_p.iloc[idx_te].values
Xtr_o, Xte_o = X_p[SRE_ORIG].iloc[idx_tr].values, X_p[SRE_ORIG].iloc[idx_te].values
ytr, yte = y_p[idx_tr], y_p[idx_te]
sc_e = StandardScaler().fit(Xtr_e); Xtr_e_s, Xte_e_s = sc_e.transform(Xtr_e), sc_e.transform(Xte_e)
sc_o = StandardScaler().fit(Xtr_o); Xtr_o_s, Xte_o_s = sc_o.transform(Xtr_o), sc_o.transform(Xte_o)
m_e = LogisticRegression(max_iter=5000).fit(Xtr_e_s, ytr)
m_o = LogisticRegression(max_iter=5000).fit(Xtr_o_s, ytr)
p_e = m_e.predict_proba(Xte_e_s)[:,1]
p_o = m_o.predict_proba(Xte_o_s)[:,1]
auc_e, lo_e, hi_e = bs_ci(yte, p_e)
auc_o, lo_o, hi_o = bs_ci(yte, p_o)
dE, dE_lo, dE_hi, dE_p = paired_delta(yte, p_e, p_o)
print(f'  SAFEHEART-Extended AUC = {auc_e:.4f} ({lo_e:.4f}-{hi_e:.4f})')
print(f'  SAFEHEART-Original AUC = {auc_o:.4f} ({lo_o:.4f}-{hi_o:.4f})')
print(f'  Δ (extended - original) = {dE:+.4f} ({dE_lo:+.4f} to {dE_hi:+.4f})  p={dE_p:.4f}')

# ---------- EXTERNAL VALIDATION ----------
def external_run(X_train, y_train, X_test, y_test, label):
    """Fit on train, freeze, apply to test. Cohort-stratified standardization."""
    sc = StandardScaler().fit(X_train.values)
    m = LogisticRegression(max_iter=5000).fit(sc.transform(X_train.values), y_train)
    p_te = m.predict_proba(sc.transform(X_test.values))[:,1]
    p_tr = m.predict_proba(sc.transform(X_train.values))[:,1]
    auc_te, lo, hi = bs_ci(y_test, p_te)
    brier = brier_score_loss(y_test, p_te)
    auc_tr = roc_auc_score(y_train, p_tr)
    # Calibration
    eps=1e-7
    lp = np.log(np.clip(p_te,eps,1-eps)/np.clip(1-p_te,eps,1-eps))
    calmod = LogisticRegression(max_iter=3000).fit(lp.reshape(-1,1), y_test)
    cal_int = float(calmod.intercept_[0]); cal_slope = float(calmod.coef_[0,0])
    print(f'  {label}: train AUC={auc_tr:.4f}, EXT AUC={auc_te:.4f} ({lo:.4f}-{hi:.4f}), '
          f'Brier={brier:.4f}, cal int={cal_int:+.3f}, slope={cal_slope:.3f}')
    return dict(label=label, train_auc=auc_tr, ext_auc=auc_te, ci_lo=lo, ci_hi=hi,
                brier=brier, cal_intercept=cal_int, cal_slope=cal_slope,
                p_te=p_te, model=m, scaler=sc, coef=m.coef_[0], features=list(X_train.columns))

print('\n[6] EXTERNAL — SAFEHEART-Extended (17 bands)')
print('  Direction A (Wales -> UKB)')
A_ext = external_run(X_w, y_w, X_u, y_u, 'A_extended_Wales_to_UKB')
print('  Direction B (UKB -> Wales)')
B_ext = external_run(X_u, y_u, X_w, y_w, 'B_extended_UKB_to_Wales')

print('\n[7] EXTERNAL — SAFEHEART-Original (9 bands)')
print('  Direction A (Wales -> UKB)')
A_orig = external_run(X_w[SRE_ORIG], y_w, X_u[SRE_ORIG], y_u, 'A_original_Wales_to_UKB')
print('  Direction B (UKB -> Wales)')
B_orig = external_run(X_u[SRE_ORIG], y_u, X_w[SRE_ORIG], y_w, 'B_original_UKB_to_Wales')

# ---------- Paired-bootstrap delta per direction ----------
print('\n[8] Paired-bootstrap ΔAUC (Extended vs Original) per direction')
dA, dA_lo, dA_hi, dA_p = paired_delta(y_u, A_ext['p_te'], A_orig['p_te'])
dB, dB_lo, dB_hi, dB_p = paired_delta(y_w, B_ext['p_te'], B_orig['p_te'])
print(f'  A (Wales->UKB): Δ={dA:+.4f}  ({dA_lo:+.4f} to {dA_hi:+.4f})  p={dA_p:.4f}')
print(f'  B (UKB->Wales): Δ={dB:+.4f}  ({dB_lo:+.4f} to {dB_hi:+.4f})  p={dB_p:.4f}')

# ---------- Coefficient inspection ----------
print('\n[9] Per-direction coefficients (top |OR-1|):')
def coef_table(model, feats, label):
    return pd.DataFrame({'feature':feats,'beta':model.coef_[0],
                         'OR_per_SD':np.exp(model.coef_[0]),'direction':label})
cA = coef_table(A_ext['model'], list(X_w.columns), 'A_Wales_trained')
cB = coef_table(B_ext['model'], list(X_u.columns), 'B_UKB_trained')
print('\n  Direction A (Wales-trained, applied to UKB):')
print(cA.assign(_a=np.abs(cA['OR_per_SD']-1)).sort_values('_a',ascending=False).drop('_a',axis=1).to_string(index=False))
print('\n  Direction B (UKB-trained, applied to Wales):')
print(cB.assign(_a=np.abs(cB['OR_per_SD']-1)).sort_values('_a',ascending=False).drop('_a',axis=1).to_string(index=False))

# ---------- Outputs ----------
os.makedirs(OUT_DIR, exist_ok=True)
results = pd.DataFrame([
    dict(metric='n_wales', value=len(w)),
    dict(metric='n_ukb', value=len(u)),
    dict(metric='events_wales', value=int(w['ascvd'].sum())),
    dict(metric='events_ukb', value=int(u['ascvd'].sum())),
    dict(metric='n_features_extended', value=X_w.shape[1]),
    dict(metric='n_features_original', value=len(SRE_ORIG)),
    # Pooled internal
    dict(metric='internal_pooled_AUC_extended', value=round(auc_e,4)),
    dict(metric='internal_pooled_AUC_original', value=round(auc_o,4)),
    dict(metric='internal_delta', value=round(dE,4)),
    dict(metric='internal_delta_p', value=round(dE_p,4)),
    # External A
    dict(metric='A_extended_AUC', value=round(A_ext['ext_auc'],4)),
    dict(metric='A_extended_CI_lo', value=round(A_ext['ci_lo'],4)),
    dict(metric='A_extended_CI_hi', value=round(A_ext['ci_hi'],4)),
    dict(metric='A_original_AUC', value=round(A_orig['ext_auc'],4)),
    dict(metric='A_original_CI_lo', value=round(A_orig['ci_lo'],4)),
    dict(metric='A_original_CI_hi', value=round(A_orig['ci_hi'],4)),
    dict(metric='A_delta', value=round(dA,4)),
    dict(metric='A_delta_CI_lo', value=round(dA_lo,4)),
    dict(metric='A_delta_CI_hi', value=round(dA_hi,4)),
    dict(metric='A_delta_p', value=round(dA_p,4)),
    dict(metric='A_brier_extended', value=round(A_ext['brier'],4)),
    dict(metric='A_calib_intercept_extended', value=round(A_ext['cal_intercept'],4)),
    dict(metric='A_calib_slope_extended', value=round(A_ext['cal_slope'],4)),
    # External B
    dict(metric='B_extended_AUC', value=round(B_ext['ext_auc'],4)),
    dict(metric='B_extended_CI_lo', value=round(B_ext['ci_lo'],4)),
    dict(metric='B_extended_CI_hi', value=round(B_ext['ci_hi'],4)),
    dict(metric='B_original_AUC', value=round(B_orig['ext_auc'],4)),
    dict(metric='B_original_CI_lo', value=round(B_orig['ci_lo'],4)),
    dict(metric='B_original_CI_hi', value=round(B_orig['ci_hi'],4)),
    dict(metric='B_delta', value=round(dB,4)),
    dict(metric='B_delta_CI_lo', value=round(dB_lo,4)),
    dict(metric='B_delta_CI_hi', value=round(dB_hi,4)),
    dict(metric='B_delta_p', value=round(dB_p,4)),
    dict(metric='B_brier_extended', value=round(B_ext['brier'],4)),
    dict(metric='B_calib_intercept_extended', value=round(B_ext['cal_intercept'],4)),
    dict(metric='B_calib_slope_extended', value=round(B_ext['cal_slope'],4)),
])
results.to_csv(f'{OUT_DIR}/CALON_SRE_extended_results.csv', index=False)
pd.concat([cA, cB], ignore_index=True).to_csv(
    f'{OUT_DIR}/CALON_SRE_extended_coefficients.csv', index=False)

print(f'\nDONE in {(datetime.datetime.now()-t0).total_seconds():.0f}s')

# ---------- SANITY GATE ----------
print('\n' + '='*72)
print('SANITY GATE: SAFEHEART-Extended vs SAFEHEART-Original on EXTERNAL data')
print('='*72)
gateA = A_ext['ext_auc'] > A_orig['ext_auc']
gateB = B_ext['ext_auc'] > B_orig['ext_auc']
print(f'  A (Wales->UKB): Extended {A_ext["ext_auc"]:.4f} '
      f'{"BEATS" if gateA else "LOSES TO"} Original {A_orig["ext_auc"]:.4f}  '
      f'(Δ={dA:+.4f}, p={dA_p:.4f})')
print(f'  B (UKB->Wales): Extended {B_ext["ext_auc"]:.4f} '
      f'{"BEATS" if gateB else "LOSES TO"} Original {B_orig["ext_auc"]:.4f}  '
      f'(Δ={dB:+.4f}, p={dB_p:.4f})')
print(f'\n  VERDICT: {"PASS" if (gateA and gateB) else ("PARTIAL" if (gateA or gateB) else "FAIL")} '
      f'(must improve external in both directions for full pass)')
print(f'\n  Harmonic mean external AUC: '
      f'extended={2*A_ext["ext_auc"]*B_ext["ext_auc"]/(A_ext["ext_auc"]+B_ext["ext_auc"]):.4f}  '
      f'original={2*A_orig["ext_auc"]*B_orig["ext_auc"]/(A_orig["ext_auc"]+B_orig["ext_auc"]):.4f}')

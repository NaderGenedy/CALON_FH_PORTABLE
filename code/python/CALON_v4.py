"""
CALON v4 — adds TUDOR treatment-adjusted untreated LDL to BOTH cohorts
======================================================================
Change from v3:
  Wales LDL: v1_ldl (TREATED, measured) ->  ldl_ut_v2 (UNTREATED, TUDOR back-calc)
  UKB LDL:   ldl_ut_v2 (UNTREATED, TUDOR back-calc) — unchanged

This removes the treatment-confounding artefact that forced us to drop
ldl_high and ldl_very_high in v3. Untreated LDL is the true biological
driver; the sign-constrained model should now KEEP these bands.

Same pipeline as v3 sign-constrained:
  - 14 candidate categorical features
  - Iterative sign-violator drop (biology-anchored)
  - Wales-trained -> UKB-external (primary, TRIPOD Type 4)
  - UKB-trained -> Wales-external (secondary, methodological check)
  - SAFEHEART-Original baseline (same external structure)

Outputs:
  output/v2/CALON_v4_results.csv
  output/v2/CALON_v4_coefficients.csv
"""
import os, warnings, datetime
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, brier_score_loss

WALES       = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
WALES_LDL_UT= r'D:/Projects/CALON_AlphaFold_Rebuild/data/wales_ldl_ut_v2.csv'
UKB_FH      = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv'
UKB_MST     = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
UKB_LDL_UT  = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'
OUT_DIR     = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2'
SEED        = 20260524
N_BOOT      = 2000

EXPECTED_SIGNS = dict(
    age_30_59     =  0,
    age_60p       = +1,
    male          = +1,
    htn           = +1,
    bmi_25_30     = +1,
    bmi_30p       = +1,
    smoking       = +1,
    ldl_high      = +1,
    lpa_high      = +1,
    apob_ldl_high = +1,
    tg_hdl_high   = +1,
    tc_hdl_high   = +1,
    t2dm          = +1,
    ldl_very_high = +1,
)
SRE_ORIG = ['age_30_59','age_60p','male','htn','bmi_25_30','bmi_30p',
            'smoking','ldl_high','lpa_high']

# ---------- Load Wales ----------
print('='*72); print('CALON v4 — sign-constrained + TUDOR untreated LDL'); print('='*72)
t0 = datetime.datetime.now()

print('\n[1] Wales PASS + TUDOR ldl_ut_v2 (row-aligned stitch)')
w_master = pd.read_csv(WALES, low_memory=False)
w_ldl    = pd.read_csv(WALES_LDL_UT)
assert len(w_master) == len(w_ldl), \
    f'Wales row mismatch: master {len(w_master)} vs ldl_ut {len(w_ldl)}'
print(f'  master rows={len(w_master):,}, ldl_ut rows={len(w_ldl):,} — row-aligned')

# Stitch (row-order alignment, since wales_ldl_ut_v2 was generated FROM pass_FULL_MASTER
# in the same row order)
w_master['ldl_ut_v2_tudor'] = w_ldl['ldl_ut_v2'].values
w_master['ldl_measured_tudor'] = w_ldl['ldl_measured'].values
w_master['reduction_applied'] = w_ldl['reduction_applied'].values
w_master['tudor_drug'] = w_ldl['dominant_drug'].values
print(f'  treated rows (reduction>0): {int((w_master["reduction_applied"]>0).sum()):,}')
print(f'  measured LDL median  (treated): {w_master.loc[w_master["reduction_applied"]>0,"ldl_measured_tudor"].median():.2f}')
print(f'  ldl_ut_v2 median     (treated): {w_master.loc[w_master["reduction_applied"]>0,"ldl_ut_v2_tudor"].median():.2f}')

w = w_master.copy()
w['family_id']        = w['family_id'].astype(str)
w['mutation_positive']= pd.to_numeric(w['mutation_positive'], errors='coerce').fillna(0).astype(int)
w['proband']          = pd.to_numeric(w.get('proband', 0), errors='coerce').fillna(0).astype(int)
w = w.sort_values(['family_id','proband'], ascending=[True,False])
w = w.drop_duplicates('family_id', keep='first').reset_index(drop=True)
print(f'  Wales family-deduped: n={len(w):,}')

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
# *** KEY CHANGE: Use TUDOR untreated LDL (back-calculated for treated patients) ***
w['ldl']        = pd.to_numeric(w['ldl_ut_v2_tudor'], errors='coerce').fillna(
                       pd.to_numeric(w['v1_ldl'], errors='coerce'))
w['hdl']        = pd.to_numeric(w['v1_hdl'], errors='coerce')
w['tg']         = pd.to_numeric(w['v1_tg'], errors='coerce')
w['tc']         = pd.to_numeric(w['v1_tc'], errors='coerce')
w['apob']       = np.nan
w['lpa']        = pd.to_numeric(w['v1_lpa'], errors='coerce')
ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd'] if c in w.columns]
for c in ev_cols:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['ascvd']      = (w[ev_cols].sum(axis=1)>0).astype(int)
print(f'  Wales events={int(w["ascvd"].sum()):,}, '
      f'using LDL: ldl_ut_v2 (treated patients back-calculated to untreated)')

# ---------- Load UKB ----------
print('\n[2] UKB FH carriers')
fh = pd.read_csv(UKB_FH, usecols=['eid'])
mst_head = pd.read_csv(UKB_MST, nrows=1, low_memory=False).columns.tolist()
need = ['eid','sex_F','age_at_recruit','bmi_direct','sbp','dbp','smoking_ever','t2dm',
        'tc_chem','hdl_chem','ldl_chem','tg_chem','apob_chem','lpa_chem',
        'prevalent_ascvd']
need = [c for c in need if c in mst_head]
u = pd.read_csv(UKB_MST, usecols=need, low_memory=False)
u = u[u['eid'].isin(fh['eid'])].copy().reset_index(drop=True)
ldl_ut = pd.read_csv(UKB_LDL_UT, usecols=['eid','ldl_ut_v2'])
u = u.merge(ldl_ut, on='eid', how='left')
u['male']        = (u['sex_F']==0).astype(int)
u['age']         = pd.to_numeric(u['age_at_recruit'], errors='coerce')
u['bmi']         = pd.to_numeric(u['bmi_direct'], errors='coerce')
u['sbp']         = pd.to_numeric(u['sbp'], errors='coerce')
u['dbp']         = pd.to_numeric(u['dbp'], errors='coerce')
u['htn']         = ((u['sbp']>=140)|(u['dbp']>=90)).astype(int)
u['t2dm']        = pd.to_numeric(u['t2dm'],errors='coerce').fillna(0).astype(int)
u['ever_smoked'] = pd.to_numeric(u['smoking_ever'],errors='coerce').fillna(0).astype(int)
# UKB already uses ldl_ut_v2 — same as before
u['ldl']         = pd.to_numeric(u['ldl_ut_v2'],errors='coerce').fillna(
                       pd.to_numeric(u['ldl_chem'],errors='coerce'))
u['hdl']         = pd.to_numeric(u['hdl_chem'],errors='coerce')
u['tg']          = pd.to_numeric(u['tg_chem'],errors='coerce')
u['tc']          = pd.to_numeric(u['tc_chem'],errors='coerce')
u['apob']        = pd.to_numeric(u['apob_chem'],errors='coerce')
u['lpa']         = pd.to_numeric(u['lpa_chem'],errors='coerce')
u['ascvd']       = pd.to_numeric(u['prevalent_ascvd'],errors='coerce').fillna(0).astype(int)
print(f'  UKB events={int(u["ascvd"].sum()):,}')

# ---------- Impute Wales ApoB from UKB ldl-apob relationship ----------
print('\n[3] Impute Wales ApoB from UKB untreated-ldl - apob relationship')
both = pd.DataFrame({'ldl': u['ldl'], 'apob': u['apob']}).dropna()
slope, intercept = np.polyfit(both['ldl'], both['apob'], 1)
print(f'  UKB: apob = {intercept:.3f} + {slope:.3f} * ldl_ut   (n={len(both)})')
w['apob'] = intercept + slope * w['ldl']
for df in [w, u]:
    for c in ['bmi','sbp','dbp','ldl','hdl','tg','tc','lpa','apob']:
        df[c] = df[c].fillna(df[c].median())

# ---------- Diagnostic: LDL distribution after TUDOR adjustment ----------
print('\n[4] LDL distribution comparison')
print(f'  Wales ldl (TUDOR untreated): median={w["ldl"].median():.2f}, '
      f'mean={w["ldl"].mean():.2f}, IQR={w["ldl"].quantile(0.25):.2f}-{w["ldl"].quantile(0.75):.2f}')
print(f'  UKB   ldl (TUDOR untreated): median={u["ldl"].median():.2f}, '
      f'mean={u["ldl"].mean():.2f}, IQR={u["ldl"].quantile(0.25):.2f}-{u["ldl"].quantile(0.75):.2f}')
print(f'  Wales ldl_high (>=4.14):  {(w["ldl"]>=4.14).mean()*100:.1f}%')
print(f'  UKB   ldl_high (>=4.14):  {(u["ldl"]>=4.14).mean()*100:.1f}%')
print(f'  Wales ldl_very_high (>=6.5): {(w["ldl"]>=6.5).mean()*100:.1f}%')
print(f'  UKB   ldl_very_high (>=6.5): {(u["ldl"]>=6.5).mean()*100:.1f}%')

# ---------- CALON v4 features ----------
def make_calon(df):
    out = pd.DataFrame(index=df.index)
    out['age_30_59']    = ((df['age']>=30)&(df['age']<60)).astype(int)
    out['age_60p']      = (df['age']>=60).astype(int)
    out['male']         = df['male'].astype(int)
    out['htn']          = df['htn'].astype(int)
    out['bmi_25_30']    = ((df['bmi']>=25)&(df['bmi']<30)).astype(int)
    out['bmi_30p']      = (df['bmi']>=30).astype(int)
    out['smoking']      = df['ever_smoked'].astype(int)
    out['ldl_high']     = (df['ldl']>=4.14).astype(int)
    out['lpa_high']     = (df['lpa']>=120).astype(int)
    out['apob_ldl_high']= (df['apob']/df['ldl'] > 0.30).astype(int)
    out['tg_hdl_high']  = (df['tg']/df['hdl'] > 2.5).astype(int)
    out['tc_hdl_high']  = (df['tc']/df['hdl'] > 5.0).astype(int)
    out['t2dm']         = df['t2dm'].astype(int)
    out['ldl_very_high']= (df['ldl']>=6.5).astype(int)
    return out

X_w = make_calon(w);   y_w = w['ascvd'].values
X_u = make_calon(u);   y_u = u['ascvd'].values
FEATURES = list(X_w.columns)
print(f'\n[5] CALON v4 candidate features: {len(FEATURES)}')
print('  Wales / UKB band prevalence:')
print(pd.DataFrame({'Wales': X_w.mean(), 'UKB': X_u.mean()}).round(3).to_string())

# ---------- Helpers ----------
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

def sign_constrained_fit(X_train, y_train, C=0.5, expected=EXPECTED_SIGNS, max_iter=20):
    feats = list(X_train.columns)
    dropped = []
    for _ in range(max_iter):
        if not feats: break
        Xtr = X_train[feats].values
        sc = StandardScaler().fit(Xtr)
        m = LogisticRegression(penalty='l2', C=C, max_iter=5000, solver='lbfgs'
                              ).fit(sc.transform(Xtr), y_train)
        violators = []
        for i, f in enumerate(feats):
            exp = expected.get(f, 0)
            if exp == 0: continue
            if (exp > 0 and m.coef_[0,i] < 0) or (exp < 0 and m.coef_[0,i] > 0):
                violators.append((f, abs(m.coef_[0,i])))
        if not violators:
            return m, sc, feats, dropped
        violators.sort(key=lambda x: -x[1])
        feats = [x for x in feats if x != violators[0][0]]
        dropped.append(violators[0][0])
    Xtr = X_train[feats].values
    sc = StandardScaler().fit(Xtr)
    m = LogisticRegression(penalty='l2', C=C, max_iter=5000, solver='lbfgs'
                          ).fit(sc.transform(Xtr), y_train)
    return m, sc, feats, dropped

def fit_and_apply(X_train, y_train, X_test, y_test, label, sign_constrained=True, C=0.5):
    if sign_constrained:
        m, sc, feats, dropped = sign_constrained_fit(X_train, y_train, C=C)
    else:
        sc = StandardScaler().fit(X_train.values)
        m = LogisticRegression(penalty='l2', C=C, max_iter=5000, solver='lbfgs'
                              ).fit(sc.transform(X_train.values), y_train)
        feats = list(X_train.columns); dropped = []
    p_te = m.predict_proba(sc.transform(X_test[feats].values))[:,1]
    p_tr = m.predict_proba(sc.transform(X_train[feats].values))[:,1]
    auc_te, lo, hi = bs_ci(y_test, p_te)
    brier = brier_score_loss(y_test, p_te)
    auc_tr = roc_auc_score(y_train, p_tr)
    eps=1e-7
    lp = np.log(np.clip(p_te,eps,1-eps)/np.clip(1-p_te,eps,1-eps))
    calmod = LogisticRegression(max_iter=3000).fit(lp.reshape(-1,1), y_test)
    cal_int = float(calmod.intercept_[0]); cal_slope = float(calmod.coef_[0,0])
    print(f'  {label}: train AUC={auc_tr:.4f}, EXT AUC={auc_te:.4f} ({lo:.4f}-{hi:.4f}), '
          f'Brier={brier:.4f}, cal int={cal_int:+.3f}, slope={cal_slope:.3f}')
    if dropped: print(f'    dropped (sign-violators): {dropped}')
    print(f'    kept ({len(feats)}): {feats}')
    return dict(label=label, train_auc=auc_tr, ext_auc=auc_te, ci_lo=lo, ci_hi=hi,
                brier=brier, cal_intercept=cal_int, cal_slope=cal_slope,
                p_te=p_te, model=m, scaler=sc, features=feats, dropped=dropped)

# ---------- External validations ----------
print('\n[6] EXTERNAL — CALON v4 (sign-constrained, TUDOR untreated LDL)')
print('  Direction A: Wales -> UKB (PRIMARY, Wales is the large training cohort)')
A_c = fit_and_apply(X_w, y_w, X_u, y_u, 'A_CALONv4_Wales_to_UKB', sign_constrained=True, C=0.5)
print('  Direction B: UKB -> Wales (secondary, methodological check)')
B_c = fit_and_apply(X_u, y_u, X_w, y_w, 'B_CALONv4_UKB_to_Wales', sign_constrained=True, C=0.5)

print('\n[7] EXTERNAL — SAFEHEART-Original baseline')
print('  Direction A:')
A_s = fit_and_apply(X_w[SRE_ORIG], y_w, X_u[SRE_ORIG], y_u, 'A_SRE_Wales_to_UKB', sign_constrained=False)
print('  Direction B:')
B_s = fit_and_apply(X_u[SRE_ORIG], y_u, X_w[SRE_ORIG], y_w, 'B_SRE_UKB_to_Wales', sign_constrained=False)

# ---------- Paired ΔAUC ----------
print('\n[8] Paired-bootstrap ΔAUC vs SAFEHEART-Original')
dA, dA_lo, dA_hi, dA_p = paired_delta(y_u, A_c['p_te'], A_s['p_te'])
dB, dB_lo, dB_hi, dB_p = paired_delta(y_w, B_c['p_te'], B_s['p_te'])
print(f'  A (Wales->UKB): Δ={dA:+.4f} ({dA_lo:+.4f} to {dA_hi:+.4f})  p={dA_p:.4f}')
print(f'  B (UKB->Wales): Δ={dB:+.4f} ({dB_lo:+.4f} to {dB_hi:+.4f})  p={dB_p:.4f}')

# ---------- Coefficients ----------
def coef_table(model, feats, label):
    return pd.DataFrame({'feature':feats,'beta':model.coef_[0],
                         'OR_per_SD':np.exp(model.coef_[0]),'direction':label})
cA = coef_table(A_c['model'], A_c['features'], 'A_Wales_trained')
cB = coef_table(B_c['model'], B_c['features'], 'B_UKB_trained')
print('\n[9] CALON v4 coefficients (sorted by |OR-1|):')
for c, name in [(cA, f'A (Wales-trained, {len(cA)} kept)'),(cB, f'B (UKB-trained, {len(cB)} kept)')]:
    print(f'\n  {name}:')
    print(c.assign(_a=np.abs(c['OR_per_SD']-1)).sort_values('_a',ascending=False).drop('_a',axis=1).to_string(index=False))

# ---------- Outputs ----------
os.makedirs(OUT_DIR, exist_ok=True)
hmean = lambda a,b: 2*a*b/(a+b)
results = pd.DataFrame([
    dict(metric='n_wales', value=len(w)),
    dict(metric='n_ukb', value=len(u)),
    dict(metric='events_wales', value=int(w['ascvd'].sum())),
    dict(metric='events_ukb', value=int(u['ascvd'].sum())),
    dict(metric='wales_treated_n', value=int((w_master['reduction_applied']>0).sum())),
    # External A
    dict(metric='A_CALONv4_AUC', value=round(A_c['ext_auc'],4)),
    dict(metric='A_CALONv4_CIlo', value=round(A_c['ci_lo'],4)),
    dict(metric='A_CALONv4_CIhi', value=round(A_c['ci_hi'],4)),
    dict(metric='A_SRE_AUC', value=round(A_s['ext_auc'],4)),
    dict(metric='A_SRE_CIlo', value=round(A_s['ci_lo'],4)),
    dict(metric='A_SRE_CIhi', value=round(A_s['ci_hi'],4)),
    dict(metric='A_delta', value=round(dA,4)),
    dict(metric='A_delta_CIlo', value=round(dA_lo,4)),
    dict(metric='A_delta_CIhi', value=round(dA_hi,4)),
    dict(metric='A_delta_p', value=round(dA_p,4)),
    dict(metric='A_brier', value=round(A_c['brier'],4)),
    dict(metric='A_cal_intercept', value=round(A_c['cal_intercept'],4)),
    dict(metric='A_cal_slope', value=round(A_c['cal_slope'],4)),
    dict(metric='A_n_features_kept', value=len(A_c['features'])),
    dict(metric='A_dropped_violators', value=','.join(A_c['dropped'])),
    # External B
    dict(metric='B_CALONv4_AUC', value=round(B_c['ext_auc'],4)),
    dict(metric='B_CALONv4_CIlo', value=round(B_c['ci_lo'],4)),
    dict(metric='B_CALONv4_CIhi', value=round(B_c['ci_hi'],4)),
    dict(metric='B_SRE_AUC', value=round(B_s['ext_auc'],4)),
    dict(metric='B_SRE_CIlo', value=round(B_s['ci_lo'],4)),
    dict(metric='B_SRE_CIhi', value=round(B_s['ci_hi'],4)),
    dict(metric='B_delta', value=round(dB,4)),
    dict(metric='B_delta_CIlo', value=round(dB_lo,4)),
    dict(metric='B_delta_CIhi', value=round(dB_hi,4)),
    dict(metric='B_delta_p', value=round(dB_p,4)),
    dict(metric='B_brier', value=round(B_c['brier'],4)),
    dict(metric='B_cal_intercept', value=round(B_c['cal_intercept'],4)),
    dict(metric='B_cal_slope', value=round(B_c['cal_slope'],4)),
    dict(metric='B_n_features_kept', value=len(B_c['features'])),
    dict(metric='B_dropped_violators', value=','.join(B_c['dropped'])),
    # Harmonic means
    dict(metric='hmean_CALONv4', value=round(hmean(A_c['ext_auc'], B_c['ext_auc']),4)),
    dict(metric='hmean_SRE', value=round(hmean(A_s['ext_auc'], B_s['ext_auc']),4)),
])
results.to_csv(f'{OUT_DIR}/CALON_v4_results.csv', index=False)
pd.concat([cA, cB], ignore_index=True).to_csv(f'{OUT_DIR}/CALON_v4_coefficients.csv', index=False)

print(f'\nDONE in {(datetime.datetime.now()-t0).total_seconds():.0f}s')

# ---------- Sanity gate vs SAFEHEART ----------
print('\n' + '='*72)
print('SANITY GATE: CALON v4 beats SAFEHEART-Original in BOTH external directions?')
print('='*72)
gA = A_c['ext_auc'] > A_s['ext_auc']
gB = B_c['ext_auc'] > B_s['ext_auc']
print(f'  A (Wales->UKB): CALON v4 {A_c["ext_auc"]:.4f} '
      f'{"BEATS" if gA else "LOSES TO"} SRE {A_s["ext_auc"]:.4f}  (Δ={dA:+.4f}, p={dA_p:.4f})')
print(f'  B (UKB->Wales): CALON v4 {B_c["ext_auc"]:.4f} '
      f'{"BEATS" if gB else "LOSES TO"} SRE {B_s["ext_auc"]:.4f}  (Δ={dB:+.4f}, p={dB_p:.4f})')
print(f'\n  VERDICT: {"PASS BOTH" if (gA and gB) else ("PARTIAL" if (gA or gB) else "FAIL")}')
print(f'  Harmonic mean external AUC: CALON v4 = {hmean(A_c["ext_auc"],B_c["ext_auc"]):.4f}  '
      f'vs SRE = {hmean(A_s["ext_auc"],B_s["ext_auc"]):.4f}')

# ---------- v3 vs v4 comparison ----------
print('\n  v3 sign-constrained reference (from prior run):')
print('    A: 0.7361  B: 0.7432  hmean: 0.7397')
print(f'  v4 (this run with TUDOR untreated LDL):')
print(f'    A: {A_c["ext_auc"]:.4f}  B: {B_c["ext_auc"]:.4f}  hmean: {hmean(A_c["ext_auc"],B_c["ext_auc"]):.4f}')

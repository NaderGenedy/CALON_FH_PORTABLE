"""
CALON v4b — proper Wales untreated-LDL via on_treatment + first_drug parsing
=============================================================================
Fix from v4: TUDOR's drug regex missed 4,512 of 4,871 treated Wales patients.
Wales master has explicit `on_treatment` (67.2% treated) and `first_drug` with
drug + dose (Atorvastatin 20 mg, Rosuvastatin 40 mg, etc.). Parse those.

Back-calculation:
  if on_treatment == 1:
    detect drug from first_drug text
    apply midpoint reduction (atorva 0.365, rosuva 0.45, simva 0.31,
                              prava 0.22, fluva 0.185, unknown 0.30)
    if ezetimibe in first_drug or second_drug: add 0.20
    cap at 0.85
    ldl_ut = ldl_measured / (1 - reduction)

UKB: unchanged (uses TUDOR's ldl_ut_v2 which works for UKB GP records).

Then: sign-constrained CALON, external validate both directions.
"""
import os, warnings, datetime, re
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, brier_score_loss

WALES     = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
UKB_FH    = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv'
UKB_MST   = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
UKB_LDL_UT= r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'
OUT_DIR   = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2'
SEED      = 20260524
N_BOOT    = 2000

STATIN_REDUCTION = dict(
    atorvastatin = 0.365, rosuvastatin = 0.45, simvastatin = 0.31,
    pravastatin  = 0.22,  fluvastatin  = 0.185, pitavastatin = 0.35,
    statin_unknown = 0.30,
)
EZE_ADD = 0.20
CAP     = 0.85

def parse_drug(txt):
    if pd.isna(txt) or txt == '' or txt == 'nan': return 'none'
    t = str(txt).lower()
    if 'atorva' in t: return 'atorvastatin'
    if 'rosuva' in t: return 'rosuvastatin'
    if 'simva'  in t: return 'simvastatin'
    if 'prava'  in t: return 'pravastatin'
    if 'fluva'  in t: return 'fluvastatin'
    if 'pitava' in t: return 'pitavastatin'
    if 'statin' in t: return 'statin_unknown'
    return 'none'

def has_eze(first, second):
    for txt in [first, second]:
        if pd.isna(txt): continue
        if 'ezetimibe' in str(txt).lower(): return True
    return False

EXPECTED_SIGNS = dict(
    age_30_59=0, age_60p=+1, male=+1, htn=+1,
    bmi_25_30=+1, bmi_30p=+1, smoking=+1,
    ldl_high=+1, lpa_high=+1, apob_ldl_high=+1,
    tg_hdl_high=+1, tc_hdl_high=+1, t2dm=+1, ldl_very_high=+1,
)
SRE_ORIG = ['age_30_59','age_60p','male','htn','bmi_25_30','bmi_30p',
            'smoking','ldl_high','lpa_high']

print('='*72); print('CALON v4b — proper Wales untreated-LDL via on_treatment'); print('='*72)
t0 = datetime.datetime.now()

# ---------- Wales ----------
print('\n[1] Wales PASS + on_treatment + first_drug parsing')
w = pd.read_csv(WALES, low_memory=False)
w['family_id']        = w['family_id'].astype(str)
w['on_treatment']     = pd.to_numeric(w['on_treatment'], errors='coerce').fillna(0).astype(int)
w['first_drug']       = w['first_drug'].fillna('').astype(str)
w['second_drug']      = w['second_drug'].fillna('').astype(str)
w['drug']             = w['first_drug'].apply(parse_drug)
# If on_treatment=1 but no specific drug detected, mark statin_unknown
w.loc[(w['on_treatment']==1) & (w['drug']=='none'), 'drug'] = 'statin_unknown'
# If on_treatment=0, force drug=none even if a drug name appears
w.loc[w['on_treatment']==0, 'drug'] = 'none'
w['has_ezetimibe']    = w.apply(lambda r: has_eze(r['first_drug'], r['second_drug']), axis=1)
print(f'  on_treatment=1: {int((w["on_treatment"]==1).sum()):,} '
      f'({(w["on_treatment"]==1).mean()*100:.1f}%)')
print(f'  drug breakdown:')
print(w['drug'].value_counts().to_string())
print(f'  ezetimibe add-on: {int(w["has_ezetimibe"].sum()):,}')

# Back-calculate untreated LDL
def back_calc(row):
    if pd.isna(row['v1_ldl']): return np.nan
    if row['on_treatment'] == 0: return row['v1_ldl']
    red = STATIN_REDUCTION.get(row['drug'], 0.0)
    if row['has_ezetimibe']: red += EZE_ADD
    red = min(red, CAP)
    if red <= 0: return row['v1_ldl']
    return row['v1_ldl'] / (1 - red)

w['v1_ldl'] = pd.to_numeric(w['v1_ldl'], errors='coerce')
w['ldl_ut_wales'] = w.apply(back_calc, axis=1)
print(f'  LDL stats:')
print(f'    measured median (treated): '
      f'{w.loc[w["on_treatment"]==1,"v1_ldl"].median():.2f}')
print(f'    ldl_ut median (treated):   '
      f'{w.loc[w["on_treatment"]==1,"ldl_ut_wales"].median():.2f}')
print(f'    measured median (naive):   '
      f'{w.loc[w["on_treatment"]==0,"v1_ldl"].median():.2f}')

# Family dedup (proband-or-first)
w['proband']          = pd.to_numeric(w.get('proband', 0), errors='coerce').fillna(0).astype(int)
w = w.sort_values(['family_id','proband'], ascending=[True,False])
w = w.drop_duplicates('family_id', keep='first').reset_index(drop=True)
w['mutation_positive']= pd.to_numeric(w['mutation_positive'], errors='coerce').fillna(0).astype(int)

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
# *** KEY: Wales LDL is the back-calculated untreated value ***
w['ldl']        = w['ldl_ut_wales']
w['hdl']        = pd.to_numeric(w['v1_hdl'], errors='coerce')
w['tg']         = pd.to_numeric(w['v1_tg'], errors='coerce')
w['tc']         = pd.to_numeric(w['v1_tc'], errors='coerce')
w['apob']       = np.nan
w['lpa']        = pd.to_numeric(w['v1_lpa'], errors='coerce')
ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd'] if c in w.columns]
for c in ev_cols:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['ascvd']      = (w[ev_cols].sum(axis=1)>0).astype(int)
print(f'  Wales family-deduped: n={len(w):,}, events={int(w["ascvd"].sum()):,}')

# ---------- UKB ----------
print('\n[2] UKB FH carriers (using TUDOR ldl_ut_v2 already)')
fh = pd.read_csv(UKB_FH, usecols=['eid'])
mst_head = pd.read_csv(UKB_MST, nrows=1, low_memory=False).columns.tolist()
need = ['eid','sex_F','age_at_recruit','bmi_direct','sbp','dbp','smoking_ever','t2dm',
        'tc_chem','hdl_chem','ldl_chem','tg_chem','apob_chem','lpa_chem','prevalent_ascvd']
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
u['ldl']         = pd.to_numeric(u['ldl_ut_v2'],errors='coerce').fillna(
                       pd.to_numeric(u['ldl_chem'],errors='coerce'))
u['hdl']         = pd.to_numeric(u['hdl_chem'],errors='coerce')
u['tg']          = pd.to_numeric(u['tg_chem'],errors='coerce')
u['tc']          = pd.to_numeric(u['tc_chem'],errors='coerce')
u['apob']        = pd.to_numeric(u['apob_chem'],errors='coerce')
u['lpa']         = pd.to_numeric(u['lpa_chem'],errors='coerce')
u['ascvd']       = pd.to_numeric(u['prevalent_ascvd'],errors='coerce').fillna(0).astype(int)
print(f'  UKB events={int(u["ascvd"].sum()):,}')

# Impute Wales ApoB
both = pd.DataFrame({'ldl': u['ldl'], 'apob': u['apob']}).dropna()
slope, intercept = np.polyfit(both['ldl'], both['apob'], 1)
w['apob'] = intercept + slope * w['ldl']
for df in [w, u]:
    for c in ['bmi','sbp','dbp','ldl','hdl','tg','tc','lpa','apob']:
        df[c] = df[c].fillna(df[c].median())

print('\n[3] LDL distribution comparison (post-untreatment-back-calc):')
print(f'  Wales ldl: median={w["ldl"].median():.2f}, IQR={w["ldl"].quantile(0.25):.2f}-{w["ldl"].quantile(0.75):.2f}')
print(f'  UKB   ldl: median={u["ldl"].median():.2f}, IQR={u["ldl"].quantile(0.25):.2f}-{u["ldl"].quantile(0.75):.2f}')
print(f'  Wales ldl_high (>=4.14):    {(w["ldl"]>=4.14).mean()*100:.1f}%')
print(f'  UKB   ldl_high (>=4.14):    {(u["ldl"]>=4.14).mean()*100:.1f}%')
print(f'  Wales ldl_very_high (>=6.5): {(w["ldl"]>=6.5).mean()*100:.1f}%')
print(f'  UKB   ldl_very_high (>=6.5): {(u["ldl"]>=6.5).mean()*100:.1f}%')

# ---------- Features ----------
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

X_w = make_calon(w); y_w = w['ascvd'].values
X_u = make_calon(u); y_u = u['ascvd'].values
print(f'\n[4] Band prevalences (Wales / UKB):')
print(pd.DataFrame({'Wales': X_w.mean(), 'UKB': X_u.mean()}).round(3).to_string())

# ---------- Sign-constrained fit + bootstrap ----------
def bs_ci(y, p, n=N_BOOT, seed=SEED):
    rng = np.random.RandomState(seed); bs=[]
    for _ in range(n):
        bi = rng.randint(0, len(y), len(y))
        if len(np.unique(y[bi]))>=2: bs.append(roc_auc_score(y[bi], p[bi]))
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
    feats = list(X_train.columns); dropped = []
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
        if not violators: return m, sc, feats, dropped
        violators.sort(key=lambda x: -x[1])
        feats = [x for x in feats if x != violators[0][0]]
        dropped.append(violators[0][0])
    Xtr = X_train[feats].values
    sc = StandardScaler().fit(Xtr)
    m = LogisticRegression(penalty='l2', C=C, max_iter=5000, solver='lbfgs'
                          ).fit(sc.transform(Xtr), y_train)
    return m, sc, feats, dropped

def external_run(X_tr, y_tr, X_te, y_te, label, sign_constrained=True, C=0.5):
    if sign_constrained:
        m, sc, feats, dropped = sign_constrained_fit(X_tr, y_tr, C=C)
    else:
        sc = StandardScaler().fit(X_tr.values)
        m = LogisticRegression(penalty='l2', C=C, max_iter=5000, solver='lbfgs'
                              ).fit(sc.transform(X_tr.values), y_tr)
        feats = list(X_tr.columns); dropped = []
    p_te = m.predict_proba(sc.transform(X_te[feats].values))[:,1]
    auc_te, lo, hi = bs_ci(y_te, p_te)
    brier = brier_score_loss(y_te, p_te)
    print(f'  {label}: EXT AUC={auc_te:.4f} ({lo:.4f}-{hi:.4f}), Brier={brier:.4f}')
    if dropped: print(f'    dropped: {dropped}')
    print(f'    kept ({len(feats)}): {feats}')
    return dict(label=label, ext_auc=auc_te, ci_lo=lo, ci_hi=hi, brier=brier,
                p_te=p_te, model=m, scaler=sc, features=feats, dropped=dropped)

print('\n[5] EXTERNAL — CALON v4b vs SAFEHEART-Original')
print('  Direction A: Wales -> UKB')
A_c = external_run(X_w, y_w, X_u, y_u, 'A_CALONv4b', sign_constrained=True)
A_s = external_run(X_w[SRE_ORIG], y_w, X_u[SRE_ORIG], y_u, 'A_SRE', sign_constrained=False)
print('  Direction B: UKB -> Wales')
B_c = external_run(X_u, y_u, X_w, y_w, 'B_CALONv4b', sign_constrained=True)
B_s = external_run(X_u[SRE_ORIG], y_u, X_w[SRE_ORIG], y_w, 'B_SRE', sign_constrained=False)

print('\n[6] Paired-bootstrap ΔAUC vs SAFEHEART-Original')
dA, dA_lo, dA_hi, dA_p = paired_delta(y_u, A_c['p_te'], A_s['p_te'])
dB, dB_lo, dB_hi, dB_p = paired_delta(y_w, B_c['p_te'], B_s['p_te'])
print(f'  A (Wales->UKB): Δ={dA:+.4f} ({dA_lo:+.4f} to {dA_hi:+.4f}) p={dA_p:.4f}')
print(f'  B (UKB->Wales): Δ={dB:+.4f} ({dB_lo:+.4f} to {dB_hi:+.4f}) p={dB_p:.4f}')

# Coefficients
def coef_table(model, feats, label):
    return pd.DataFrame({'feature':feats, 'beta':model.coef_[0],
                         'OR_per_SD':np.exp(model.coef_[0]), 'direction':label})
cA = coef_table(A_c['model'], A_c['features'], 'A_Wales_trained')
cB = coef_table(B_c['model'], B_c['features'], 'B_UKB_trained')
print('\n[7] CALON v4b coefficients (sorted by |OR-1|):')
for c, name in [(cA, f'A (Wales-trained, {len(cA)} kept)'),(cB, f'B (UKB-trained, {len(cB)} kept)')]:
    print(f'\n  {name}:')
    print(c.assign(_a=np.abs(c['OR_per_SD']-1)).sort_values('_a',ascending=False).drop('_a',axis=1).to_string(index=False))

os.makedirs(OUT_DIR, exist_ok=True)
hmean = lambda a,b: 2*a*b/(a+b)
pd.DataFrame([
    dict(metric='n_wales', value=len(w)),
    dict(metric='n_ukb', value=len(u)),
    dict(metric='wales_treated_n', value=int((w_master := pd.read_csv(WALES, usecols=['on_treatment']))['on_treatment'].fillna(0).astype(int).sum())),
    dict(metric='A_CALONv4b_AUC', value=round(A_c['ext_auc'],4)),
    dict(metric='A_CALONv4b_CIlo', value=round(A_c['ci_lo'],4)),
    dict(metric='A_CALONv4b_CIhi', value=round(A_c['ci_hi'],4)),
    dict(metric='A_SRE_AUC', value=round(A_s['ext_auc'],4)),
    dict(metric='A_delta', value=round(dA,4)),
    dict(metric='A_delta_p', value=round(dA_p,4)),
    dict(metric='B_CALONv4b_AUC', value=round(B_c['ext_auc'],4)),
    dict(metric='B_CALONv4b_CIlo', value=round(B_c['ci_lo'],4)),
    dict(metric='B_CALONv4b_CIhi', value=round(B_c['ci_hi'],4)),
    dict(metric='B_SRE_AUC', value=round(B_s['ext_auc'],4)),
    dict(metric='B_delta', value=round(dB,4)),
    dict(metric='B_delta_p', value=round(dB_p,4)),
    dict(metric='hmean_CALONv4b', value=round(hmean(A_c['ext_auc'], B_c['ext_auc']),4)),
    dict(metric='hmean_SRE', value=round(hmean(A_s['ext_auc'], B_s['ext_auc']),4)),
    dict(metric='A_dropped', value=','.join(A_c['dropped'])),
    dict(metric='B_dropped', value=','.join(B_c['dropped'])),
]).to_csv(f'{OUT_DIR}/CALON_v4b_results.csv', index=False)
pd.concat([cA, cB], ignore_index=True).to_csv(f'{OUT_DIR}/CALON_v4b_coefficients.csv', index=False)

print(f'\nDONE in {(datetime.datetime.now()-t0).total_seconds():.0f}s')

# Sanity gate
print('\n' + '='*72)
print('SANITY GATE: CALON v4b beats SAFEHEART-Original in BOTH directions?')
print('='*72)
gA = A_c['ext_auc'] > A_s['ext_auc']; gB = B_c['ext_auc'] > B_s['ext_auc']
print(f'  A: CALON v4b {A_c["ext_auc"]:.4f} {"BEATS" if gA else "LOSES TO"} SRE {A_s["ext_auc"]:.4f} (Δ={dA:+.4f}, p={dA_p:.4f})')
print(f'  B: CALON v4b {B_c["ext_auc"]:.4f} {"BEATS" if gB else "LOSES TO"} SRE {B_s["ext_auc"]:.4f} (Δ={dB:+.4f}, p={dB_p:.4f})')
print(f'\n  VERDICT: {"PASS BOTH" if (gA and gB) else ("PARTIAL" if (gA or gB) else "FAIL")}')
print(f'  hmean: CALON v4b={hmean(A_c["ext_auc"],B_c["ext_auc"]):.4f}  vs SRE={hmean(A_s["ext_auc"],B_s["ext_auc"]):.4f}')
print('\nReference (prior runs):')
print('  v3 sign-constrained:  A=0.7361  B=0.7432  hmean=0.7397')
print('  v4 (TUDOR LDL):       A=0.7361  B=0.7430  hmean=0.7396')
print(f'  v4b (on_treatment):   A={A_c["ext_auc"]:.4f}  B={B_c["ext_auc"]:.4f}  hmean={hmean(A_c["ext_auc"],B_c["ext_auc"]):.4f}')

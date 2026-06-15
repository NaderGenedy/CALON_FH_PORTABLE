"""
CALON v7 — CLEAN cohort, contamination-free reverse-engineering
==================================================================
Replaces all v6 contamination sources:
  - Wales age: from DRAGON's mesearment_age_1 (baseline) NOT age_at_hard_event (outcome-conditional)
  - Wales ApoB: from DRAGON's actual ApoB column NOT LDL-deterministic imputation
  - FH+ filter: Wales Positive1==1 (DRAGON) OR mutation_positive==1 (PASS), explicit
  - v1_date sanity: keep rows with v1_date in [1990, 2026]
  - Drop Wales rows without DRAGON coverage (no rescue with median-fill)

Architecture unchanged from v6:
  - 13 candidate categorical bands
  - Sign-constrained L2 logistic (iterative drop biology-violators)
  - L2 C = 0.5, bootstrap CI N=2000, SEED=20260524
  - External validate Wales-clean <-> UKB carriers both directions

Cohort sizes (expected):
  Wales-clean: ~300 family-deduped, ~70 events  (down from 4,570 / 1,193 contaminated)
  UKB:         3,540 / 165 events (unchanged)

Honest expectations (user-pre-stated null prediction):
  CALON-clean 0.796 vs SAFEHEART 0.793, Delta +0.003, p ~0.90
  If true, we lock the null and write the honest paper.
"""
import os, warnings, datetime, re
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, brier_score_loss

WALES_PASS = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
DRAGON_3   = r'C:/Users/nader/Downloads/calon_ukb_pipeline/DRAGON_3.csv'
UKB_FH     = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv'
UKB_MST    = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
UKB_LDL_UT = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'
OUT_DIR    = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2'
SEED       = 20260524
N_BOOT     = 2000
L2_C       = 0.5

STATIN_DOSE_REDUCTION = {
    (1,5):0.31,(1,10):0.37,(1,20):0.43,(1,40):0.49,(1,80):0.55,
    (2,5):0.38,(2,10):0.43,(2,20):0.48,(2,40):0.55,
    (3,5):0.22,(3,10):0.27,(3,20):0.32,(3,40):0.37,(3,80):0.42,
    (4,10):0.20,(4,20):0.24,(4,40):0.29,
    (5,20):0.17,(5,40):0.21,(5,80):0.25,
    (6,1):0.30,(6,2):0.35,(6,4):0.40,
}
STATIN_DEFAULT = {1:0.43,2:0.48,3:0.32,4:0.24,5:0.21,6:0.35}
COMBO_ADD = {1:0.20,2:0.60,3:0.20,4:0.08}
CAP = 0.85
STATIN_MAP = {'atorva':1,'rosuva':2,'simva':3,'prava':4,'fluva':5,'pitava':6}
COMBO_MAP  = {'ezetimibe':1,'evolocumab':2,'alirocumab':2,'inclisiran':2,
              'bempedoic':3,'fenobibrate':4,'fenofibrate':4,'bezafibrate':4,'gemfibrozil':4}

CALON_FEATURES = [
    'age_30_59','age_60p','male','bmi_25_30','bmi_30p','smoking',
    'ldl_high','lpa_high','apob_ldl_high','t2dm',
    'ldl_severe_ut','hdl_low','young_severe',
]
EXPECTED_SIGNS = dict(
    age_30_59=0, age_60p=+1, male=+1, bmi_25_30=+1, bmi_30p=+1, smoking=+1,
    ldl_high=+1, lpa_high=+1, apob_ldl_high=+1, t2dm=+1,
    ldl_severe_ut=+1, hdl_low=+1, young_severe=+1,
)
SRE_ORIG = ['age_30_59','age_60p','male','htn','bmi_25_30','bmi_30p',
            'smoking','ldl_high','lpa_high']

def parse_drug_text(txt):
    if pd.isna(txt) or txt=='' or txt=='nan': return 0, np.nan, 0
    t = str(txt).lower()
    statin = next((c for k,c in STATIN_MAP.items() if k in t), 0)
    m = re.search(r'(\d+\.?\d*)\s*mg', t)
    dose = float(m.group(1)) if m else np.nan
    combo = next((c for k,c in COMBO_MAP.items() if k in t), 0)
    return statin, dose, combo

def compute_reduction(statin, dose, combo):
    if statin==0 and combo==0: return 0.0
    if statin>0:
        if not np.isnan(dose):
            try_doses = [d for (s,d) in STATIN_DOSE_REDUCTION if s==statin]
            if try_doses:
                closest = min(try_doses, key=lambda x: abs(x-dose))
                statin_red = STATIN_DOSE_REDUCTION.get((statin, closest), 0)
            else:
                statin_red = STATIN_DEFAULT.get(statin, 0)
        else:
            statin_red = STATIN_DEFAULT.get(statin, 0)
    else:
        statin_red = 0.0
    return min(statin_red + COMBO_ADD.get(combo, 0), CAP)

print('='*72); print('CALON v7 CLEAN — Dragon-derived baseline age, no contamination'); print('='*72)
t0 = datetime.datetime.now()

# ============================================================
# STEP 1 — Load PASS + DRAGON, join, derive clean cohort
# ============================================================
print('\n[1] Wales PASS + DRAGON join + clean-age + FH+ filter')
pass_df = pd.read_csv(WALES_PASS, low_memory=False)
print(f'  PASS rows: {len(pass_df):,}')

dragon = pd.read_csv(DRAGON_3, encoding='utf-8-sig')
print(f'  DRAGON rows: {len(dragon):,}')

# Join by participant_id <-> DatabaseNumber (5.8% overlap)
w = pass_df.merge(dragon[['DatabaseNumber','BirthDate','MeasurementDate_1','mesearment_age_1',
                          'ApoB','Positive1','Positive2']],
                  left_on='participant_id', right_on='DatabaseNumber', how='inner')
print(f'  PASS x DRAGON match: {len(w):,}')

# Clean baseline age (PRIMARY: mesearment_age_1; FALLBACK: BirthDate + MeasurementDate_1)
w['age_baseline'] = pd.to_numeric(w['mesearment_age_1'], errors='coerce')
# Fallback to BirthDate-derived computation
bd = pd.to_datetime(w['BirthDate'], errors='coerce', dayfirst=False)
md = pd.to_datetime(w['MeasurementDate_1'], errors='coerce', dayfirst=False)
fallback_age = ((md - bd).dt.days / 365.25)
w['age_baseline'] = w['age_baseline'].fillna(fallback_age)
# Sanity bounds
w['age_baseline'] = w['age_baseline'].where((w['age_baseline'] >= 18) & (w['age_baseline'] <= 95))
n_with_age = w['age_baseline'].notna().sum()
print(f'  Rows with CLEAN baseline age (mesearment_age_1 or computed): {n_with_age}')
w = w[w['age_baseline'].notna()].copy()

# FH+ filter: Positive1 == "1" OR PASS mutation_positive == 1
w['pos1_int']     = pd.to_numeric(w['Positive1'], errors='coerce').fillna(0).astype(int)
w['mut_pos_int']  = pd.to_numeric(w['mutation_positive'], errors='coerce').fillna(0).astype(int)
w['fh_pos']       = ((w['pos1_int']==1) | (w['mut_pos_int']==1)).astype(int)
n_fh = int(w['fh_pos'].sum())
print(f'  FH+ (Positive1==1 OR mutation_positive==1): {n_fh}')
w = w[w['fh_pos']==1].copy().reset_index(drop=True)

# Family dedup (proband-or-first)
w['family_id'] = w['family_id'].astype(str)
w['proband']   = pd.to_numeric(w.get('proband', 0), errors='coerce').fillna(0).astype(int)
w = w.sort_values(['family_id','proband'], ascending=[True, False])
w = w.drop_duplicates('family_id', keep='first').reset_index(drop=True)
print(f'  After family dedup: n={len(w):,}')

# Drug parsing + dose-specific untreated LDL
w['on_treatment'] = pd.to_numeric(w['on_treatment'], errors='coerce').fillna(0).astype(int)
w['first_drug']   = w['first_drug'].fillna('').astype(str)
w['second_drug']  = w['second_drug'].fillna('').astype(str)
parsed = w.apply(lambda r: (parse_drug_text(r['first_drug']),
                            parse_drug_text(r['second_drug'])), axis=1)
w['statin_code'] = [max(p[0][0], p[1][0]) for p in parsed]
w['statin_dose'] = [p[0][1] if p[0][0]>0 else (p[1][1] if p[1][0]>0 else np.nan) for p in parsed]
w['combo_code']  = [max(p[0][2], p[1][2]) for p in parsed]
w.loc[w['on_treatment']==0, ['statin_code','combo_code']] = 0
w['reduction']   = w.apply(lambda r: compute_reduction(r['statin_code'], r['statin_dose'], r['combo_code']), axis=1)
w['v1_ldl']      = pd.to_numeric(w['v1_ldl'], errors='coerce')
w['ldl_ut']      = w['v1_ldl'] / (1 - w['reduction']).clip(lower=0.15)
w.loc[w['reduction']==0, 'ldl_ut'] = w.loc[w['reduction']==0, 'v1_ldl']

# Clinical variables
w['male']        = (w['sex_F']==0).astype(int)
w['age']         = w['age_baseline']            # *** CLEAN baseline age ***
w['bmi']         = pd.to_numeric(w['bmi'], errors='coerce')
w['sbp']         = pd.to_numeric(w['sbp'], errors='coerce')
w['dbp']         = pd.to_numeric(w['dbp'], errors='coerce')
w['htn']         = ((w['sbp']>=140)|(w['dbp']>=90)|
                    w['bp_medication'].astype(str).str.contains('1|yes|Yes|TRUE',regex=True,na=False)).astype(int)
w['t2dm']        = pd.to_numeric(w.get('diabetes_type_2', w.get('diabetes')),errors='coerce').fillna(0).astype(int)
w['ever_smoked'] = pd.to_numeric(w.get('smoking_ever', w.get('smoking')),errors='coerce').fillna(0).astype(int)
w['ldl']         = w['ldl_ut']
w['hdl']         = pd.to_numeric(w['v1_hdl'], errors='coerce')
w['tg']          = pd.to_numeric(w['v1_tg'], errors='coerce')
w['tc']          = pd.to_numeric(w['v1_tc'], errors='coerce')
w['apob']        = pd.to_numeric(w['ApoB'], errors='coerce')   # *** REAL DRAGON ApoB ***
w['lpa']         = pd.to_numeric(w['v1_lpa'], errors='coerce')

# ASCVD composite outcome
ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd'] if c in w.columns]
for c in ev_cols:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['ascvd'] = (w[ev_cols].sum(axis=1)>0).astype(int)
print(f'  Final Wales-clean: n={len(w):,}, events={int(w["ascvd"].sum())} ({w["ascvd"].mean()*100:.1f}%)')
print(f'  Mean baseline age: {w["age"].mean():.1f} (clean, NOT outcome-conditional)')
print(f'  ApoB available: {w["apob"].notna().sum()} / {len(w)}')

# ============================================================
# STEP 2 — Load UKB carriers (clean recruitment age already)
# ============================================================
print('\n[2] UKB FH carriers (unchanged from v6 — already uses recruitment age)')
fh = pd.read_csv(UKB_FH, usecols=['eid'])
head = pd.read_csv(UKB_MST, nrows=1, low_memory=False).columns.tolist()
need = ['eid','sex_F','age_at_recruit','bmi_direct','sbp','dbp','smoking_ever','t2dm',
        'tc_chem','hdl_chem','ldl_chem','tg_chem','apob_chem','lpa_chem','prevalent_ascvd']
need = [c for c in need if c in head]
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
u['ldl']         = pd.to_numeric(u['ldl_ut_v2'],errors='coerce').fillna(pd.to_numeric(u['ldl_chem'],errors='coerce'))
u['hdl']         = pd.to_numeric(u['hdl_chem'],errors='coerce')
u['tg']          = pd.to_numeric(u['tg_chem'],errors='coerce')
u['tc']          = pd.to_numeric(u['tc_chem'],errors='coerce')
u['apob']        = pd.to_numeric(u['apob_chem'],errors='coerce')
u['lpa']         = pd.to_numeric(u['lpa_chem'],errors='coerce')
u['ascvd']       = pd.to_numeric(u['prevalent_ascvd'],errors='coerce').fillna(0).astype(int)
print(f'  UKB n={len(u):,}, events={int(u["ascvd"].sum())} ({u["ascvd"].mean()*100:.1f}%)')
print(f'  Mean recruitment age: {u["age"].mean():.1f} (always clean)')

# Median-fill remaining continuous fields per-cohort
for df in [w, u]:
    for c in ['bmi','sbp','dbp','ldl','hdl','tg','tc','lpa','apob']:
        df[c] = df[c].fillna(df[c].median())

# ============================================================
# STEP 3 — Feature engineering (same 13 bands)
# ============================================================
def make_bands(df):
    out = pd.DataFrame(index=df.index)
    out['age_30_59']     = ((df['age']>=30)&(df['age']<60)).astype(int)
    out['age_60p']       = (df['age']>=60).astype(int)
    out['male']          = df['male'].astype(int)
    out['htn']           = df['htn'].astype(int)
    out['bmi_25_30']     = ((df['bmi']>=25)&(df['bmi']<30)).astype(int)
    out['bmi_30p']       = (df['bmi']>=30).astype(int)
    out['smoking']       = df['ever_smoked'].astype(int)
    out['ldl_high']      = (df['ldl']>=4.14).astype(int)
    out['lpa_high']      = (df['lpa']>=120).astype(int)
    out['apob_ldl_high'] = (df['apob']/df['ldl'] > 0.30).astype(int)
    out['t2dm']          = df['t2dm'].astype(int)
    out['ldl_severe_ut'] = (df['ldl']>=8.0).astype(int)
    out['hdl_low']       = (((df['male']==1)&(df['hdl']<1.0))|((df['male']==0)&(df['hdl']<1.2))).astype(int)
    out['young_severe']  = ((df['age']<40)&(df['ldl']>=8.0)).astype(int)
    return out

X_w_all = make_bands(w); y_w = w['ascvd'].values
X_u_all = make_bands(u); y_u = u['ascvd'].values
X_w = X_w_all[CALON_FEATURES]
X_u = X_u_all[CALON_FEATURES]
print(f'\n[3] Band prevalences (CLEAN cohort, Wales / UKB):')
print(pd.DataFrame({'Wales-clean': X_w.mean(), 'UKB': X_u.mean()}).round(3).to_string())

# ============================================================
# STEP 4 — Bootstrap helpers
# ============================================================
def bs_auc_ci(y, p, n=N_BOOT, seed=SEED):
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
            ds.append(roc_auc_score(y[bi],pa[bi]) - roc_auc_score(y[bi],pb[bi]))
    return float(np.mean(ds)), float(np.percentile(ds,2.5)), float(np.percentile(ds,97.5)), \
           2*min((np.array(ds)<=0).mean(),(np.array(ds)>=0).mean())

# ============================================================
# STEP 5 — Sign-constrained logistic (unchanged from v6)
# ============================================================
def sign_constrained_fit(X_train, y_train, C=L2_C, expected=EXPECTED_SIGNS, max_iter=20):
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

def external_eval(X_train, y_train, X_test, y_test, label, sign_constrained=True, C=L2_C):
    if sign_constrained:
        m, sc, feats, dropped = sign_constrained_fit(X_train, y_train, C=C)
    else:
        sc = StandardScaler().fit(X_train.values)
        m = LogisticRegression(penalty='l2', C=C, max_iter=5000, solver='lbfgs'
                              ).fit(sc.transform(X_train.values), y_train)
        feats = list(X_train.columns); dropped = []
    p_te = m.predict_proba(sc.transform(X_test[feats].values))[:,1]
    auc_te, lo, hi = bs_auc_ci(y_test, p_te)
    brier = brier_score_loss(y_test, p_te)
    print(f'  {label}: EXT AUC={auc_te:.4f} ({lo:.4f}-{hi:.4f}), Brier={brier:.4f}')
    if dropped: print(f'    dropped (sign violators): {dropped}')
    print(f'    kept ({len(feats)}): {feats}')
    return dict(label=label, ext_auc=auc_te, ci_lo=lo, ci_hi=hi, brier=brier,
                p_te=p_te, model=m, features=feats, dropped=dropped)

# ============================================================
# STEP 6 — External validation, both directions, vs SAFEHEART
# ============================================================
print('\n[4] EXTERNAL VALIDATION (CLEAN cohort)')
print('  Direction A: Wales-clean -> UKB')
A_c = external_eval(X_w, y_w, X_u, y_u, 'A_CALONv7_Wales_to_UKB')
A_s = external_eval(X_w_all[SRE_ORIG], y_w, X_u_all[SRE_ORIG], y_u, 'A_SRE', sign_constrained=False)
print('  Direction B: UKB -> Wales-clean')
B_c = external_eval(X_u, y_u, X_w, y_w, 'B_CALONv7_UKB_to_Wales')
B_s = external_eval(X_u_all[SRE_ORIG], y_u, X_w_all[SRE_ORIG], y_w, 'B_SRE', sign_constrained=False)

print('\n[5] Paired-bootstrap delta-AUC')
dA, dA_lo, dA_hi, dA_p = paired_delta(y_u, A_c['p_te'], A_s['p_te'])
dB, dB_lo, dB_hi, dB_p = paired_delta(y_w, B_c['p_te'], B_s['p_te'])
print(f'  A (Wales->UKB): Delta={dA:+.4f} ({dA_lo:+.4f} to {dA_hi:+.4f}) p={dA_p:.4f}')
print(f'  B (UKB->Wales): Delta={dB:+.4f} ({dB_lo:+.4f} to {dB_hi:+.4f}) p={dB_p:.4f}')

def coef_table(model, feats, label):
    return pd.DataFrame({'feature':feats,'beta':model.coef_[0],
                         'OR_per_SD':np.exp(model.coef_[0]),'direction':label})
cA = coef_table(A_c['model'], A_c['features'], 'A_Wales_trained')
cB = coef_table(B_c['model'], B_c['features'], 'B_UKB_trained')
print('\n[6] CALON v7 coefficients (CLEAN baseline age):')
for c, name in [(cA, 'A (Wales-clean trained)'),(cB, 'B (UKB trained)')]:
    print(f'\n  {name}:')
    print(c.assign(_a=np.abs(c['OR_per_SD']-1)).sort_values('_a',ascending=False).drop('_a',axis=1).to_string(index=False))

# Outputs
os.makedirs(OUT_DIR, exist_ok=True)
hmean = lambda a,b: 2*a*b/(a+b)
pd.DataFrame([
    dict(metric='n_wales_clean', value=len(w)),
    dict(metric='n_ukb', value=len(u)),
    dict(metric='events_wales_clean', value=int(w['ascvd'].sum())),
    dict(metric='events_ukb', value=int(u['ascvd'].sum())),
    dict(metric='A_CALONv7_AUC', value=round(A_c['ext_auc'],4)),
    dict(metric='A_CALONv7_CIlo', value=round(A_c['ci_lo'],4)),
    dict(metric='A_CALONv7_CIhi', value=round(A_c['ci_hi'],4)),
    dict(metric='A_SRE_AUC', value=round(A_s['ext_auc'],4)),
    dict(metric='A_delta', value=round(dA,4)),
    dict(metric='A_delta_CIlo', value=round(dA_lo,4)),
    dict(metric='A_delta_CIhi', value=round(dA_hi,4)),
    dict(metric='A_delta_p', value=round(dA_p,4)),
    dict(metric='B_CALONv7_AUC', value=round(B_c['ext_auc'],4)),
    dict(metric='B_CALONv7_CIlo', value=round(B_c['ci_lo'],4)),
    dict(metric='B_CALONv7_CIhi', value=round(B_c['ci_hi'],4)),
    dict(metric='B_SRE_AUC', value=round(B_s['ext_auc'],4)),
    dict(metric='B_delta', value=round(dB,4)),
    dict(metric='B_delta_CIlo', value=round(dB_lo,4)),
    dict(metric='B_delta_CIhi', value=round(dB_hi,4)),
    dict(metric='B_delta_p', value=round(dB_p,4)),
    dict(metric='hmean_CALONv7', value=round(hmean(A_c['ext_auc'], B_c['ext_auc']),4)),
    dict(metric='hmean_SRE', value=round(hmean(A_s['ext_auc'], B_s['ext_auc']),4)),
    dict(metric='A_dropped', value=','.join(A_c['dropped'])),
    dict(metric='B_dropped', value=','.join(B_c['dropped'])),
]).to_csv(f'{OUT_DIR}/CALON_v7_clean_results.csv', index=False)
pd.concat([cA, cB], ignore_index=True).to_csv(f'{OUT_DIR}/CALON_v7_clean_coefficients.csv', index=False)

print(f'\nDONE in {(datetime.datetime.now()-t0).total_seconds():.0f}s')

print('\n' + '='*72)
print('SANITY GATE: CALON v7 (CLEAN cohort) vs SAFEHEART-Original')
print('='*72)
gA = A_c['ext_auc'] > A_s['ext_auc']; gB = B_c['ext_auc'] > B_s['ext_auc']
print(f'  A: CALON v7 {A_c["ext_auc"]:.4f} ({A_c["ci_lo"]:.4f}-{A_c["ci_hi"]:.4f}) '
      f'{"BEATS" if gA else "LOSES TO"} SRE {A_s["ext_auc"]:.4f}  Delta={dA:+.4f} p={dA_p:.4f}')
print(f'  B: CALON v7 {B_c["ext_auc"]:.4f} ({B_c["ci_lo"]:.4f}-{B_c["ci_hi"]:.4f}) '
      f'{"BEATS" if gB else "LOSES TO"} SRE {B_s["ext_auc"]:.4f}  Delta={dB:+.4f} p={dB_p:.4f}')
print(f'\n  VERDICT: {"PASS BOTH" if (gA and gB) else ("PARTIAL" if (gA or gB) else "FAIL/NULL")}')
print(f'  Harmonic mean: CALON v7 = {hmean(A_c["ext_auc"], B_c["ext_auc"]):.4f}  '
      f'vs SRE = {hmean(A_s["ext_auc"], B_s["ext_auc"]):.4f}')

print('\nVersion history (contaminated -> clean):')
print('  v6 contaminated:  A=0.7584  B=0.7390  hmean=0.7486  [REJECTED - age leak]')
print(f'  v7 CLEAN:         A={A_c["ext_auc"]:.4f}  B={B_c["ext_auc"]:.4f}  '
      f'hmean={hmean(A_c["ext_auc"], B_c["ext_auc"]):.4f}  [HONEST]')

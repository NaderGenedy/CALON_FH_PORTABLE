"""
CALON v4c — proper Wales untreated-LDL via dose-specific drug parsing
======================================================================
Drug encoding (per Wales clinical convention):
  X.Y where X = statin type, Y = combination LLT
    X: 0=none, 1=Atorva, 2=Rosuva, 3=Simva, 4=Prava, 5=Fluva, 6=Pitava
    Y: 0=none, 1=Ezetimibe, 2=PCSK9i, 3=Bempedoic, 4=Fibrate
  Examples: 1.0 = Atorvastatin alone
            1.1 = Atorvastatin + Ezetimibe
            2.2 = Rosuvastatin + PCSK9i (evolocumab)
            0.4 = Fibrate alone (no statin)

Dose-specific reduction factors (per drug Methods 2.3 + dose-response literature):
                10mg   20mg   40mg   80mg
  Atorvastatin: 0.37   0.43   0.49   0.55
  Rosuvastatin: 0.43   0.48   0.55   0.63       (5mg=0.38)
  Simvastatin:  0.27   0.32   0.37   0.42
  Pravastatin:  0.20   0.24   0.29   ---
  Fluvastatin:  0.17   0.21   ---    ---
  Pitavastatin: 0.30   0.35   0.40   ---
Combination add-ons:
  Ezetimibe:   +0.20
  PCSK9i:      +0.60
  Bempedoic:   +0.20
  Fibrate:     +0.08 (modest LDL effect, mostly TG-lowering)
Cap at 0.85.

ldl_ut = ldl_measured / (1 - total_reduction)   for treated patients.

Outputs:
  output/v2/CALON_v4c_results.csv
  output/v2/CALON_v4c_coefficients.csv
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

# Dose-specific LDL reduction (fraction) — keyed by (statin_code, dose_mg)
STATIN_DOSE_REDUCTION = {
    # Atorvastatin
    (1, 5):  0.31, (1, 10): 0.37, (1, 20): 0.43, (1, 40): 0.49, (1, 80): 0.55,
    # Rosuvastatin
    (2, 5):  0.38, (2, 10): 0.43, (2, 20): 0.48, (2, 40): 0.55,
    # Simvastatin
    (3, 5):  0.22, (3, 10): 0.27, (3, 20): 0.32, (3, 40): 0.37, (3, 80): 0.42,
    # Pravastatin
    (4, 10): 0.20, (4, 20): 0.24, (4, 40): 0.29,
    # Fluvastatin
    (5, 20): 0.17, (5, 40): 0.21, (5, 80): 0.25,
    # Pitavastatin
    (6, 1):  0.30, (6, 2):  0.35, (6, 4):  0.40,
}
STATIN_DEFAULT_REDUCTION = {1:0.43, 2:0.48, 3:0.32, 4:0.24, 5:0.21, 6:0.35}  # if dose unknown

COMBINATION_ADD = {1:0.20, 2:0.60, 3:0.20, 4:0.08}  # Eze, PCSK9i, Bempedoic, Fibrate
CAP = 0.85

STATIN_MAP = {
    'atorva':1, 'rosuva':2, 'simva':3, 'prava':4, 'fluva':5, 'pitava':6,
}
COMBO_MAP = {
    'ezetimibe':1, 'evolocumab':2, 'alirocumab':2, 'inclisiran':2,
    'bempedoic':3, 'fenobibrate':4, 'fenofibrate':4, 'bezafibrate':4,
    'gemfibrozil':4,
}

def parse_drug_text(txt):
    """Return (statin_code, dose_mg, combo_code) from a drug text."""
    if pd.isna(txt) or txt == '' or txt == 'nan':
        return 0, np.nan, 0
    t = str(txt).lower()
    # Statin detection
    statin_code = 0
    for key, code in STATIN_MAP.items():
        if key in t: statin_code = code; break
    # Dose detection (number followed by mg)
    m = re.search(r'(\d+\.?\d*)\s*mg', t)
    dose_mg = float(m.group(1)) if m else np.nan
    # Combination detection
    combo_code = 0
    for key, code in COMBO_MAP.items():
        if key in t: combo_code = code; break
    return statin_code, dose_mg, combo_code

def parse_row(first, second):
    """Returns (statin_code, dose_mg, combo_code, encoding_label).
    Combines first_drug + second_drug to detect combination therapy."""
    s1, d1, c1 = parse_drug_text(first)
    s2, d2, c2 = parse_drug_text(second)
    # Pick the statin (likely in first slot)
    statin = s1 if s1 > 0 else s2
    dose = d1 if s1 > 0 and not np.isnan(d1) else (d2 if s2 > 0 and not np.isnan(d2) else np.nan)
    # Pick the combination (could be in either slot)
    combo = max(c1, c2)
    label = f'{statin}.{combo}'
    return statin, dose, combo, label

def compute_reduction(statin_code, dose_mg, combo_code):
    """Total LDL reduction fraction. 0 = treatment-naive."""
    if statin_code == 0 and combo_code == 0:
        return 0.0
    # Statin component
    if statin_code > 0:
        if not np.isnan(dose_mg):
            # Round dose to nearest standard value in table
            try_doses = [d for (s, d) in STATIN_DOSE_REDUCTION if s == statin_code]
            if try_doses:
                closest = min(try_doses, key=lambda x: abs(x - dose_mg))
                statin_red = STATIN_DOSE_REDUCTION.get((statin_code, closest), 0)
            else:
                statin_red = STATIN_DEFAULT_REDUCTION.get(statin_code, 0)
        else:
            statin_red = STATIN_DEFAULT_REDUCTION.get(statin_code, 0)
    else:
        statin_red = 0.0
    # Combination add-on (can apply without a statin too, for combo-only)
    combo_red = COMBINATION_ADD.get(combo_code, 0)
    total = statin_red + combo_red
    return min(total, CAP)

EXPECTED_SIGNS = dict(
    age_30_59=0, age_60p=+1, male=+1, htn=+1,
    bmi_25_30=+1, bmi_30p=+1, smoking=+1,
    ldl_high=+1, lpa_high=+1, apob_ldl_high=+1,
    tg_hdl_high=+1, tc_hdl_high=+1, t2dm=+1, ldl_very_high=+1,
)
SRE_ORIG = ['age_30_59','age_60p','male','htn','bmi_25_30','bmi_30p',
            'smoking','ldl_high','lpa_high']

print('='*72); print('CALON v4c — dose-specific Wales untreated-LDL'); print('='*72)
t0 = datetime.datetime.now()

# ---------- Wales ----------
print('\n[1] Wales PASS + dose-specific drug parsing')
w = pd.read_csv(WALES, low_memory=False)
w['family_id']        = w['family_id'].astype(str)
w['on_treatment']     = pd.to_numeric(w['on_treatment'], errors='coerce').fillna(0).astype(int)
w['first_drug']       = w['first_drug'].fillna('').astype(str)
w['second_drug']      = w['second_drug'].fillna('').astype(str)

# Parse every row
parsed = w.apply(lambda r: parse_row(r['first_drug'], r['second_drug']), axis=1)
w['statin_code'] = [p[0] for p in parsed]
w['statin_dose'] = [p[1] for p in parsed]
w['combo_code']  = [p[2] for p in parsed]
w['drug_encoding']= [p[3] for p in parsed]
# Force no-treatment if on_treatment=0
w.loc[w['on_treatment']==0, ['statin_code','combo_code']] = 0
w.loc[w['on_treatment']==0, 'statin_dose'] = np.nan
w.loc[w['on_treatment']==0, 'drug_encoding'] = '0.0'

print(f'  on_treatment=1: {int((w["on_treatment"]==1).sum()):,}')
print(f'\n  Drug encoding distribution (X.Y):')
enc = w['drug_encoding'].value_counts().head(20)
print('  ' + enc.to_string().replace('\n', '\n  '))

# Compute LDL reduction
w['reduction_calc'] = w.apply(lambda r: compute_reduction(
    r['statin_code'], r['statin_dose'], r['combo_code']), axis=1)
w['v1_ldl']         = pd.to_numeric(w['v1_ldl'], errors='coerce')

# Back-calc untreated LDL
w['ldl_ut_calon']   = w['v1_ldl'] / (1 - w['reduction_calc']).clip(lower=0.15)
w.loc[w['reduction_calc'] == 0, 'ldl_ut_calon'] = w.loc[w['reduction_calc'] == 0, 'v1_ldl']

print(f'\n  LDL back-calc summary:')
trt = w[w['reduction_calc']>0]
print(f'    treated n: {len(trt):,}')
print(f'    measured LDL median (treated): {trt["v1_ldl"].median():.2f}')
print(f'    reduction median (treated):    {trt["reduction_calc"].median()*100:.1f}%')
print(f'    ldl_ut median (treated):       {trt["ldl_ut_calon"].median():.2f}')

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
w['ldl']        = w['ldl_ut_calon']
w['hdl']        = pd.to_numeric(w['v1_hdl'], errors='coerce')
w['tg']         = pd.to_numeric(w['v1_tg'], errors='coerce')
w['tc']         = pd.to_numeric(w['v1_tc'], errors='coerce')
w['apob']       = np.nan
w['lpa']        = pd.to_numeric(w['v1_lpa'], errors='coerce')
ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd'] if c in w.columns]
for c in ev_cols:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['ascvd']      = (w[ev_cols].sum(axis=1)>0).astype(int)
print(f'\n  Wales family-deduped: n={len(w):,}, events={int(w["ascvd"].sum()):,}')

# ---------- UKB ----------
print('\n[2] UKB FH carriers')
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

# Impute Wales ApoB from UKB ldl_ut <-> ApoB regression
both = pd.DataFrame({'ldl': u['ldl'], 'apob': u['apob']}).dropna()
slope, intercept = np.polyfit(both['ldl'], both['apob'], 1)
print(f'  Wales ApoB imputation: apob = {intercept:.3f} + {slope:.3f} * ldl_ut')
w['apob'] = intercept + slope * w['ldl']
for df in [w, u]:
    for c in ['bmi','sbp','dbp','ldl','hdl','tg','tc','lpa','apob']:
        df[c] = df[c].fillna(df[c].median())

print('\n[3] LDL distribution (post dose-specific back-calc):')
print(f'  Wales ldl median={w["ldl"].median():.2f}, IQR={w["ldl"].quantile(0.25):.2f}-{w["ldl"].quantile(0.75):.2f}')
print(f'  UKB   ldl median={u["ldl"].median():.2f}, IQR={u["ldl"].quantile(0.25):.2f}-{u["ldl"].quantile(0.75):.2f}')
print(f'  Wales ldl_high (>=4.14):     {(w["ldl"]>=4.14).mean()*100:.1f}%')
print(f'  UKB   ldl_high (>=4.14):     {(u["ldl"]>=4.14).mean()*100:.1f}%')
print(f'  Wales ldl_very_high (>=6.5): {(w["ldl"]>=6.5).mean()*100:.1f}%')
print(f'  UKB   ldl_very_high (>=6.5): {(u["ldl"]>=6.5).mean()*100:.1f}%')

# ---------- Features + sign-constrained ----------
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
                p_te=p_te, model=m, features=feats, dropped=dropped)

print('\n[5] EXTERNAL — CALON v4c vs SAFEHEART-Original')
print('  Direction A: Wales -> UKB')
A_c = external_run(X_w, y_w, X_u, y_u, 'A_CALONv4c')
A_s = external_run(X_w[SRE_ORIG], y_w, X_u[SRE_ORIG], y_u, 'A_SRE', sign_constrained=False)
print('  Direction B: UKB -> Wales')
B_c = external_run(X_u, y_u, X_w, y_w, 'B_CALONv4c')
B_s = external_run(X_u[SRE_ORIG], y_u, X_w[SRE_ORIG], y_w, 'B_SRE', sign_constrained=False)

print('\n[6] Paired-bootstrap ΔAUC')
dA, dA_lo, dA_hi, dA_p = paired_delta(y_u, A_c['p_te'], A_s['p_te'])
dB, dB_lo, dB_hi, dB_p = paired_delta(y_w, B_c['p_te'], B_s['p_te'])
print(f'  A (Wales->UKB): Δ={dA:+.4f} ({dA_lo:+.4f} to {dA_hi:+.4f}) p={dA_p:.4f}')
print(f'  B (UKB->Wales): Δ={dB:+.4f} ({dB_lo:+.4f} to {dB_hi:+.4f}) p={dB_p:.4f}')

def coef_table(model, feats, label):
    return pd.DataFrame({'feature':feats,'beta':model.coef_[0],
                         'OR_per_SD':np.exp(model.coef_[0]),'direction':label})
cA = coef_table(A_c['model'], A_c['features'], 'A_Wales_trained')
cB = coef_table(B_c['model'], B_c['features'], 'B_UKB_trained')
print('\n[7] CALON v4c coefficients:')
for c, name in [(cA, f'A (Wales-trained, {len(cA)} kept)'),(cB, f'B (UKB-trained, {len(cB)} kept)')]:
    print(f'\n  {name}:')
    print(c.assign(_a=np.abs(c['OR_per_SD']-1)).sort_values('_a',ascending=False).drop('_a',axis=1).to_string(index=False))

os.makedirs(OUT_DIR, exist_ok=True)
hmean = lambda a,b: 2*a*b/(a+b)
pd.DataFrame([
    dict(metric='n_wales', value=len(w)),
    dict(metric='n_ukb', value=len(u)),
    dict(metric='A_CALONv4c_AUC', value=round(A_c['ext_auc'],4)),
    dict(metric='A_CALONv4c_CIlo', value=round(A_c['ci_lo'],4)),
    dict(metric='A_CALONv4c_CIhi', value=round(A_c['ci_hi'],4)),
    dict(metric='A_SRE_AUC', value=round(A_s['ext_auc'],4)),
    dict(metric='A_delta', value=round(dA,4)),
    dict(metric='A_delta_p', value=round(dA_p,4)),
    dict(metric='B_CALONv4c_AUC', value=round(B_c['ext_auc'],4)),
    dict(metric='B_CALONv4c_CIlo', value=round(B_c['ci_lo'],4)),
    dict(metric='B_CALONv4c_CIhi', value=round(B_c['ci_hi'],4)),
    dict(metric='B_SRE_AUC', value=round(B_s['ext_auc'],4)),
    dict(metric='B_delta', value=round(dB,4)),
    dict(metric='B_delta_p', value=round(dB_p,4)),
    dict(metric='hmean_CALONv4c', value=round(hmean(A_c['ext_auc'], B_c['ext_auc']),4)),
    dict(metric='hmean_SRE', value=round(hmean(A_s['ext_auc'], B_s['ext_auc']),4)),
    dict(metric='A_dropped', value=','.join(A_c['dropped'])),
    dict(metric='B_dropped', value=','.join(B_c['dropped'])),
]).to_csv(f'{OUT_DIR}/CALON_v4c_results.csv', index=False)
pd.concat([cA, cB], ignore_index=True).to_csv(f'{OUT_DIR}/CALON_v4c_coefficients.csv', index=False)

print(f'\nDONE in {(datetime.datetime.now()-t0).total_seconds():.0f}s')

print('\n' + '='*72)
print('SANITY GATE: CALON v4c (dose-specific Wales LDL back-calc) vs SAFEHEART-Original')
print('='*72)
gA = A_c['ext_auc'] > A_s['ext_auc']; gB = B_c['ext_auc'] > B_s['ext_auc']
print(f'  A: CALON v4c {A_c["ext_auc"]:.4f} {"BEATS" if gA else "LOSES TO"} SRE {A_s["ext_auc"]:.4f} (Δ={dA:+.4f}, p={dA_p:.4f})')
print(f'  B: CALON v4c {B_c["ext_auc"]:.4f} {"BEATS" if gB else "LOSES TO"} SRE {B_s["ext_auc"]:.4f} (Δ={dB:+.4f}, p={dB_p:.4f})')
print(f'\n  VERDICT: {"PASS BOTH" if (gA and gB) else ("PARTIAL" if (gA or gB) else "FAIL")}')
print(f'  hmean: CALON v4c={hmean(A_c["ext_auc"],B_c["ext_auc"]):.4f}  vs SRE={hmean(A_s["ext_auc"],B_s["ext_auc"]):.4f}')
print('\nReference:')
print('  v3 (no LDL back-calc):       A=0.7361  B=0.7432  hmean=0.7397')
print('  v4 (TUDOR LDL):              A=0.7361  B=0.7430  hmean=0.7396')
print(f'  v4c (dose-specific Wales):   A={A_c["ext_auc"]:.4f}  B={B_c["ext_auc"]:.4f}  hmean={hmean(A_c["ext_auc"],B_c["ext_auc"]):.4f}')

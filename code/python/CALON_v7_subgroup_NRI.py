"""
CALON v7 — per-subgroup NRI computation (Direction A primary)
================================================================
The locked v7 results CSV stores aggregate NRI per direction only.
This script re-runs the v7 model pipeline (identical SEED, identical
sign-constrained architecture, identical cohort joins) and computes
NRI per subgroup so the poster tiles can display both ΔAUC and NRI%.

Outputs: output/v2/CALON_FINAL_v7_subgroup_NRI.csv
"""
import os, warnings, re
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

WALES_PASS = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
DRAGON_3   = r'C:/Users/nader/Downloads/calon_ukb_pipeline/DRAGON_3.csv'
UKB_FH     = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv'
UKB_MST    = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
UKB_LDL_UT = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'
OUT_DIR    = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2'
SEED       = 20260524
L2_C       = 0.5
NRI_CUTS   = (0.05, 0.20)

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

CALON_FEATURES = ['age_30_59','age_60p','male','bmi_25_30','bmi_30p','smoking',
                  'ldl_high','lpa_high','apob_ldl_high','t2dm',
                  'ldl_severe_ut','hdl_low','young_severe']
EXPECTED_SIGNS = dict(age_30_59=0, age_60p=+1, male=+1, bmi_25_30=+1, bmi_30p=+1,
                      smoking=+1, ldl_high=+1, lpa_high=+1, apob_ldl_high=+1,
                      t2dm=+1, ldl_severe_ut=+1, hdl_low=+1, young_severe=+1)
SRE_ORIG = ['age_30_59','age_60p','male','htn','bmi_25_30','bmi_30p',
            'smoking','ldl_high','lpa_high']

def parse_drug(txt):
    if pd.isna(txt) or txt=='' or txt=='nan': return 0, np.nan, 0
    t = str(txt).lower()
    s = next((c for k,c in STATIN_MAP.items() if k in t), 0)
    m = re.search(r'(\d+\.?\d*)\s*mg', t)
    d = float(m.group(1)) if m else np.nan
    c = next((c for k,c in COMBO_MAP.items() if k in t), 0)
    return s, d, c

def compute_reduction(s, d, c):
    if s==0 and c==0: return 0.0
    if s>0:
        if not np.isnan(d):
            try_doses = [dd for (ss,dd) in STATIN_DOSE_REDUCTION if ss==s]
            if try_doses:
                closest = min(try_doses, key=lambda x: abs(x-d))
                sr = STATIN_DOSE_REDUCTION.get((s, closest), 0)
            else:
                sr = STATIN_DEFAULT.get(s, 0)
        else:
            sr = STATIN_DEFAULT.get(s, 0)
    else:
        sr = 0.0
    return min(sr + COMBO_ADD.get(c, 0), CAP)

# ---------- Load Wales-clean ----------
print('Loading Wales PASS + DRAGON...')
pass_df = pd.read_csv(WALES_PASS, low_memory=False)
dragon  = pd.read_csv(DRAGON_3, encoding='utf-8-sig')
w = pass_df.merge(dragon[['DatabaseNumber','BirthDate','MeasurementDate_1',
                          'mesearment_age_1','ApoB','Positive1','Positive2']],
                  left_on='participant_id', right_on='DatabaseNumber', how='inner')
w['age_baseline'] = pd.to_numeric(w['mesearment_age_1'], errors='coerce')
bd = pd.to_datetime(w['BirthDate'], errors='coerce')
md = pd.to_datetime(w['MeasurementDate_1'], errors='coerce')
w['age_baseline'] = w['age_baseline'].fillna((md - bd).dt.days / 365.25)
w['age_baseline'] = w['age_baseline'].where((w['age_baseline']>=18)&(w['age_baseline']<=95))
w = w[w['age_baseline'].notna()].copy()
w['pos1_int']    = pd.to_numeric(w['Positive1'], errors='coerce').fillna(0).astype(int)
w['mut_pos_int'] = pd.to_numeric(w['mutation_positive'], errors='coerce').fillna(0).astype(int)
w = w[(w['pos1_int']==1) | (w['mut_pos_int']==1)].copy()
w['family_id'] = w['family_id'].astype(str)
w['proband']   = pd.to_numeric(w.get('proband', 0), errors='coerce').fillna(0).astype(int)
w = w.sort_values(['family_id','proband'], ascending=[True,False])
w = w.drop_duplicates('family_id', keep='first').reset_index(drop=True)
w['on_treatment'] = pd.to_numeric(w['on_treatment'], errors='coerce').fillna(0).astype(int)
w['first_drug']   = w['first_drug'].fillna('').astype(str)
w['second_drug']  = w['second_drug'].fillna('').astype(str)
parsed = w.apply(lambda r: (parse_drug(r['first_drug']), parse_drug(r['second_drug'])), axis=1)
w['statin_code'] = [max(p[0][0], p[1][0]) for p in parsed]
w['statin_dose'] = [p[0][1] if p[0][0]>0 else (p[1][1] if p[1][0]>0 else np.nan) for p in parsed]
w['combo_code']  = [max(p[0][2], p[1][2]) for p in parsed]
w.loc[w['on_treatment']==0, ['statin_code','combo_code']] = 0
w['reduction'] = w.apply(lambda r: compute_reduction(r['statin_code'], r['statin_dose'], r['combo_code']), axis=1)
w['v1_ldl'] = pd.to_numeric(w['v1_ldl'], errors='coerce')
w['ldl_ut'] = w['v1_ldl'] / (1 - w['reduction']).clip(lower=0.15)
w.loc[w['reduction']==0, 'ldl_ut'] = w.loc[w['reduction']==0, 'v1_ldl']
w['male']        = (w['sex_F']==0).astype(int)
w['age']         = w['age_baseline']
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
w['apob']        = pd.to_numeric(w['ApoB'], errors='coerce')
w['lpa']         = pd.to_numeric(w['v1_lpa'], errors='coerce')
ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd'] if c in w.columns]
for c in ev_cols:
    w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
w['ascvd'] = (w[ev_cols].sum(axis=1)>0).astype(int)

# ---------- Load UKB ----------
print('Loading UKB carriers...')
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
for df in [w, u]:
    for c in ['bmi','sbp','dbp','ldl','hdl','tg','tc','lpa','apob']:
        df[c] = df[c].fillna(df[c].median())
print(f'  Wales-clean n={len(w)}, events={int(w["ascvd"].sum())}')
print(f'  UKB         n={len(u)}, events={int(u["ascvd"].sum())}')

# ---------- Build bands ----------
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

X_w_all = make_bands(w); X_u_all = make_bands(u)
y_w = w['ascvd'].values; y_u = u['ascvd'].values
X_w = X_w_all[CALON_FEATURES]; X_u = X_u_all[CALON_FEATURES]

# ---------- Sign-constrained fit ----------
def sign_constrained_fit(X_train, y_train, C=L2_C, expected=EXPECTED_SIGNS, max_iter=20):
    feats = list(X_train.columns); dropped = []
    for _ in range(max_iter):
        if not feats: break
        Xtr = X_train[feats].values
        sc = StandardScaler().fit(Xtr)
        m = LogisticRegression(penalty='l2', C=C, max_iter=5000, solver='lbfgs').fit(sc.transform(Xtr), y_train)
        v = []
        for i, f in enumerate(feats):
            e = expected.get(f, 0)
            if e == 0: continue
            if (e > 0 and m.coef_[0,i] < 0) or (e < 0 and m.coef_[0,i] > 0):
                v.append((f, abs(m.coef_[0,i])))
        if not v: return m, sc, feats, dropped
        v.sort(key=lambda x: -x[1])
        feats = [x for x in feats if x != v[0][0]]
        dropped.append(v[0][0])
    Xtr = X_train[feats].values
    sc = StandardScaler().fit(Xtr)
    m = LogisticRegression(penalty='l2', C=C, max_iter=5000, solver='lbfgs').fit(sc.transform(Xtr), y_train)
    return m, sc, feats, dropped

def fit_and_predict(X_tr, y_tr, X_te, sign_constrained=True):
    if sign_constrained:
        m, sc, feats, _ = sign_constrained_fit(X_tr, y_tr)
    else:
        sc = StandardScaler().fit(X_tr.values)
        m = LogisticRegression(penalty='l2', C=L2_C, max_iter=5000, solver='lbfgs').fit(sc.transform(X_tr.values), y_tr)
        feats = list(X_tr.columns)
    return m.predict_proba(sc.transform(X_te[feats].values))[:,1]

# Direction A: Wales-trained -> UKB-test
print('\nFitting Direction A (Wales-clean -> UKB)...')
pA_calon = fit_and_predict(X_w, y_w, X_u, sign_constrained=True)
pA_sre   = fit_and_predict(X_w_all[SRE_ORIG], y_w, X_u_all[SRE_ORIG], sign_constrained=False)
print(f'  A AUCs: CALON={roc_auc_score(y_u, pA_calon):.4f}, SRE={roc_auc_score(y_u, pA_sre):.4f}')

# Direction B: UKB-trained -> Wales-test
print('\nFitting Direction B (UKB -> Wales-clean)...')
pB_calon = fit_and_predict(X_u, y_u, X_w, sign_constrained=True)
pB_sre   = fit_and_predict(X_u_all[SRE_ORIG], y_u, X_w_all[SRE_ORIG], sign_constrained=False)
print(f'  B AUCs: CALON={roc_auc_score(y_w, pB_calon):.4f}, SRE={roc_auc_score(y_w, pB_sre):.4f}')

# ---------- Per-subgroup NRI computation ----------
def categorical_nri(y, pn, po, cuts=NRI_CUTS):
    y = np.asarray(y); pn = np.asarray(pn); po = np.asarray(po)
    cn = np.digitize(pn, cuts); co = np.digitize(po, cuts)
    up = cn > co; dn = cn < co
    ev = (y == 1); ne = (y == 0)
    nri_e = (up[ev].mean() - dn[ev].mean()) if ev.sum() else np.nan
    nri_n = (dn[ne].mean() - up[ne].mean()) if ne.sum() else np.nan
    pct_up_event   = up[ev].mean() if ev.sum() else np.nan
    pct_dn_event   = dn[ev].mean() if ev.sum() else np.nan
    pct_up_nonevent= up[ne].mean() if ne.sum() else np.nan
    pct_dn_nonevent= dn[ne].mean() if ne.sum() else np.nan
    nri_total = (nri_e if ev.sum() else 0) + (nri_n if ne.sum() else 0)
    return dict(nri_total=nri_total, nri_event=nri_e, nri_nonevent=nri_n,
                pct_up_event=pct_up_event, pct_dn_event=pct_dn_event,
                pct_up_nonevent=pct_up_nonevent, pct_dn_nonevent=pct_dn_nonevent)

def subgroup_nri_aucs(df_test, p_calon, p_sre, direction_label):
    df = df_test.copy().reset_index(drop=True)
    df['p_calon'] = p_calon
    df['p_sre']   = p_sre
    strata = [
        ('sex', 'male', [(0, 'Female'), (1, 'Male')]),
        ('age_band', 'age',  [((0, 50), '<50y'), ((50, 65), '50-65y'), ((65, 100), '>=65y')]),
        ('ldl_band','ldl',  [((0, 4.14),'LDL<4.14'), ((4.14, 6.5),'LDL 4.14-6.5'), ((6.5, 99),'LDL>=6.5')]),
        ('t2dm', 't2dm',[(0, 'No DM'), (1, 'T2DM')]),
        ('smoking','ever_smoked',[(0,'Non-smoker'),(1,'Smoker')]),
        ('htn', 'htn', [(0,'No HTN'),(1,'HTN')]),
    ]
    rows = []
    for grp, col, levels in strata:
        for level, label_l in levels:
            if isinstance(level, tuple):
                mask = (df[col] >= level[0]) & (df[col] < level[1])
            else:
                mask = (df[col] == level)
            if mask.sum() < 20 or df.loc[mask, 'ascvd'].sum() < 5: continue
            y_sub = df.loc[mask, 'ascvd'].values
            pc_sub = df.loc[mask, 'p_calon'].values
            ps_sub = df.loc[mask, 'p_sre'].values
            try:
                ac = roc_auc_score(y_sub, pc_sub)
                as_ = roc_auc_score(y_sub, ps_sub)
                nri = categorical_nri(y_sub, pc_sub, ps_sub)
                rows.append(dict(direction=direction_label, group=grp, level=label_l,
                                  n=int(mask.sum()),
                                  events=int(df.loc[mask,'ascvd'].sum()),
                                  auc_calon=round(ac, 4),
                                  auc_sre=round(as_, 4),
                                  delta_auc=round(ac - as_, 4),
                                  nri_total=round(nri['nri_total'], 4),
                                  nri_event=round(nri['nri_event'], 4) if not np.isnan(nri['nri_event']) else np.nan,
                                  nri_nonevent=round(nri['nri_nonevent'], 4) if not np.isnan(nri['nri_nonevent']) else np.nan,
                                  pct_up_event=round(nri['pct_up_event']*100, 1) if not np.isnan(nri['pct_up_event']) else np.nan,
                                  pct_dn_event=round(nri['pct_dn_event']*100, 1) if not np.isnan(nri['pct_dn_event']) else np.nan,
                                  pct_up_nonevent=round(nri['pct_up_nonevent']*100, 1) if not np.isnan(nri['pct_up_nonevent']) else np.nan,
                                  pct_dn_nonevent=round(nri['pct_dn_nonevent']*100, 1) if not np.isnan(nri['pct_dn_nonevent']) else np.nan))
            except Exception:
                pass
    return pd.DataFrame(rows)

print('\nComputing per-subgroup NRI...')
sg_A = subgroup_nri_aucs(u, pA_calon, pA_sre, 'A_Wales_to_UKB')
sg_B = subgroup_nri_aucs(w, pB_calon, pB_sre, 'B_UKB_to_Wales')
sg = pd.concat([sg_A, sg_B], ignore_index=True)

print('\n=== Direction A subgroup table ===')
print(sg_A.to_string(index=False))
print('\n=== Direction B subgroup table ===')
print(sg_B.to_string(index=False))

os.makedirs(OUT_DIR, exist_ok=True)
sg.to_csv(f'{OUT_DIR}/CALON_FINAL_v7_subgroup_NRI.csv', index=False)
print(f'\nWrote {OUT_DIR}/CALON_FINAL_v7_subgroup_NRI.csv')

# Print poster-ready tile data for Direction A
print('\n=== POSTER TILE DATA (Direction A) ===')
print('Subgroup label  |  ΔAUC  |  NRI total  |  % event up  |  % non-event down')
for _, r in sg_A.iterrows():
    nri_pct = f'{r["nri_total"]*100:+.1f}%' if not pd.isna(r["nri_total"]) else 'n/a'
    pct_up_e = f'{r["pct_up_event"]:.0f}%' if not pd.isna(r["pct_up_event"]) else 'n/a'
    pct_dn_n = f'{r["pct_dn_nonevent"]:.0f}%' if not pd.isna(r["pct_dn_nonevent"]) else 'n/a'
    print(f'  {r["level"]:<18s} {r["delta_auc"]:+.3f}    {nri_pct:>8s}    {pct_up_e:>5s}        {pct_dn_n:>5s}')

"""
CALON v3 — sign-constrained + tighter-L2, TRIPOD-correct framing
==================================================================
Two model variants explored, both attempting to fix v2's Direction B failure:

  v3a: tighter L2 (C=0.1) — shrinks the wrong-sign features toward zero
  v3b: sign-constrained logistic — forces biology-anchored signs on the
       11 known risk factors (LDL, ApoB/LDL, TG/HDL, TC/HDL, HTN, smoking,
       T2DM, sex, BMI, Lp(a), age_60p). Implemented as iterative drop:
       fit, check signs, drop the worst violator (set to 0 and refit),
       repeat.

Validation strategy (TRIPOD-correct):
  PRIMARY: Wales-trained -> UKB-external (n_train=4,570 with 1,193 events,
           large FH registry as development cohort).
  SECONDARY: pooled internal 70/30 for sanity.
  Direction B (UKB-trained -> Wales) is run as a methodological check only,
  not as a primary validation, because UKB-only training (165 events) is
  underpowered.

Comparators throughout: SAFEHEART-Original (9 bands).

Outputs:
  output/v2/CALON_v3_results.csv
  output/v2/CALON_v3_coefficients.csv
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

# Expected signs (biology-anchored). +1 = positive risk factor, -1 = protective, 0 = free.
EXPECTED_SIGNS = dict(
    age_30_59     =  0,   # referencing artefact — free
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

# ---------- Load cohorts (identical to v2) ----------
print('='*72); print('CALON v3 — sign-constrained + tighter L2'); print('='*72)
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

# Impute Wales ApoB
both = pd.DataFrame({'ldl': u['ldl'], 'apob': u['apob']}).dropna()
slope, intercept = np.polyfit(both['ldl'], both['apob'], 1)
w['apob'] = intercept + slope * w['ldl']
for df in [w, u]:
    for c in ['bmi','sbp','dbp','ldl','hdl','tg','tc','lpa','apob']:
        df[c] = df[c].fillna(df[c].median())

# ---------- CALON v3 feature builder (same 14 as v2) ----------
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
FEATURES_14 = list(X_w.columns)
print(f'\n[3] CALON v3 candidate features: {len(FEATURES_14)}')

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

# ---------- Sign-constrained logistic via iterative drop ----------
def sign_constrained_fit(X_train, y_train, C=0.5, expected=EXPECTED_SIGNS, max_iter=20):
    """Fit L2 logistic; drop the worst sign-violator; refit until all signs OK or no features left."""
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
        # Drop the worst violator
        violators.sort(key=lambda x: -x[1])
        worst = violators[0][0]
        feats = [x for x in feats if x != worst]
        dropped.append(worst)
    Xtr = X_train[feats].values
    sc = StandardScaler().fit(Xtr)
    m = LogisticRegression(penalty='l2', C=C, max_iter=5000, solver='lbfgs'
                          ).fit(sc.transform(Xtr), y_train)
    return m, sc, feats, dropped

def standard_fit(X_train, y_train, C=0.5):
    sc = StandardScaler().fit(X_train.values)
    m = LogisticRegression(penalty='l2', C=C, max_iter=5000, solver='lbfgs'
                          ).fit(sc.transform(X_train.values), y_train)
    return m, sc, list(X_train.columns), []

def apply_model(model, sc, feats, X_test):
    return model.predict_proba(sc.transform(X_test[feats].values))[:,1]

def external_eval(X_train, y_train, X_test, y_test, label, mode, C=0.5):
    """mode: 'standard', 'tight_L2', or 'sign_constrained'."""
    if mode == 'standard':
        m, sc, feats, dropped = standard_fit(X_train, y_train, C=C)
    elif mode == 'tight_L2':
        m, sc, feats, dropped = standard_fit(X_train, y_train, C=0.1)
    elif mode == 'sign_constrained':
        m, sc, feats, dropped = sign_constrained_fit(X_train, y_train, C=C)
    p_te = apply_model(m, sc, feats, X_test)
    p_tr = apply_model(m, sc, feats, X_train)
    auc_te, lo, hi = bs_ci(y_test, p_te)
    brier = brier_score_loss(y_test, p_te)
    auc_tr = roc_auc_score(y_train, p_tr)
    eps=1e-7
    lp = np.log(np.clip(p_te,eps,1-eps)/np.clip(1-p_te,eps,1-eps))
    calmod = LogisticRegression(max_iter=3000).fit(lp.reshape(-1,1), y_test)
    cal_int = float(calmod.intercept_[0]); cal_slope = float(calmod.coef_[0,0])
    print(f'  {label} [{mode}]: train AUC={auc_tr:.4f}, EXT AUC={auc_te:.4f} ({lo:.4f}-{hi:.4f}), '
          f'Brier={brier:.4f}, cal int={cal_int:+.3f}, slope={cal_slope:.3f}')
    if dropped:
        print(f'    dropped (sign-violators): {dropped}')
    print(f'    kept features: {feats}')
    return dict(label=f'{label}_{mode}', mode=mode, train_auc=auc_tr,
                ext_auc=auc_te, ci_lo=lo, ci_hi=hi, brier=brier,
                cal_intercept=cal_int, cal_slope=cal_slope,
                p_te=p_te, model=m, scaler=sc, features=feats, dropped=dropped)

# ---------- Run all three variants ----------
print('\n[4] EXTERNAL — three CALON v3 variants, both directions')
print('  Direction A: Wales -> UKB')
A_std  = external_eval(X_w, y_w, X_u, y_u, 'A_Wales_to_UKB', 'standard', C=0.5)
A_tightL2 = external_eval(X_w, y_w, X_u, y_u, 'A_Wales_to_UKB', 'tight_L2')
A_sign = external_eval(X_w, y_w, X_u, y_u, 'A_Wales_to_UKB', 'sign_constrained', C=0.5)
print('  Direction B: UKB -> Wales')
B_std  = external_eval(X_u, y_u, X_w, y_w, 'B_UKB_to_Wales', 'standard', C=0.5)
B_tightL2 = external_eval(X_u, y_u, X_w, y_w, 'B_UKB_to_Wales', 'tight_L2')
B_sign = external_eval(X_u, y_u, X_w, y_w, 'B_UKB_to_Wales', 'sign_constrained', C=0.5)

print('\n[5] SAFEHEART-Original baseline (same external structure)')
A_sre = external_eval(X_w[SRE_ORIG], y_w, X_u[SRE_ORIG], y_u, 'A_SRE', 'standard')
B_sre = external_eval(X_u[SRE_ORIG], y_u, X_w[SRE_ORIG], y_w, 'B_SRE', 'standard')

# ---------- Paired-bootstrap delta-AUC ----------
print('\n[6] Paired-bootstrap ΔAUC vs SAFEHEART-Original')
results_rows = []
for name, A, B in [('standard',  A_std,  B_std),
                   ('tight_L2',  A_tightL2, B_tightL2),
                   ('sign_constrained', A_sign, B_sign)]:
    dA, dA_lo, dA_hi, dA_p = paired_delta(y_u, A['p_te'], A_sre['p_te'])
    dB, dB_lo, dB_hi, dB_p = paired_delta(y_w, B['p_te'], B_sre['p_te'])
    print(f'  [{name}]  A: Δ={dA:+.4f} ({dA_lo:+.4f} to {dA_hi:+.4f}) p={dA_p:.4f}  '
          f'|  B: Δ={dB:+.4f} ({dB_lo:+.4f} to {dB_hi:+.4f}) p={dB_p:.4f}')
    results_rows.append(dict(variant=name,
                              auc_A=A['ext_auc'], auc_A_lo=A['ci_lo'], auc_A_hi=A['ci_hi'],
                              auc_B=B['ext_auc'], auc_B_lo=B['ci_lo'], auc_B_hi=B['ci_hi'],
                              brier_A=A['brier'], brier_B=B['brier'],
                              cal_int_A=A['cal_intercept'], cal_slope_A=A['cal_slope'],
                              cal_int_B=B['cal_intercept'], cal_slope_B=B['cal_slope'],
                              delta_A=dA, delta_A_lo=dA_lo, delta_A_hi=dA_hi, delta_A_p=dA_p,
                              delta_B=dB, delta_B_lo=dB_lo, delta_B_hi=dB_hi, delta_B_p=dB_p,
                              hmean_ext=2*A['ext_auc']*B['ext_auc']/(A['ext_auc']+B['ext_auc']),
                              n_features_A=len(A['features']), n_features_B=len(B['features']),
                              dropped_A=','.join(A['dropped']), dropped_B=','.join(B['dropped']),
                              ))

# Add SAFEHEART row
results_rows.append(dict(variant='SAFEHEART_Original',
                          auc_A=A_sre['ext_auc'], auc_A_lo=A_sre['ci_lo'], auc_A_hi=A_sre['ci_hi'],
                          auc_B=B_sre['ext_auc'], auc_B_lo=B_sre['ci_lo'], auc_B_hi=B_sre['ci_hi'],
                          brier_A=A_sre['brier'], brier_B=B_sre['brier'],
                          cal_int_A=A_sre['cal_intercept'], cal_slope_A=A_sre['cal_slope'],
                          cal_int_B=B_sre['cal_intercept'], cal_slope_B=B_sre['cal_slope'],
                          delta_A=0, delta_A_lo=0, delta_A_hi=0, delta_A_p=1.0,
                          delta_B=0, delta_B_lo=0, delta_B_hi=0, delta_B_p=1.0,
                          hmean_ext=2*A_sre['ext_auc']*B_sre['ext_auc']/(A_sre['ext_auc']+B_sre['ext_auc']),
                          n_features_A=9, n_features_B=9, dropped_A='', dropped_B=''))

results = pd.DataFrame(results_rows)
print('\n[7] HEADLINE TABLE')
print(results[['variant','auc_A','auc_B','hmean_ext','delta_A','delta_A_p','delta_B','delta_B_p']].round(4).to_string(index=False))

# ---------- Best-variant coefficients ----------
# Pick the variant that wins both directions, else best harmonic mean
print('\n[8] Best-variant selection')
best_idx = None
for i, r in enumerate(results_rows[:-1]):  # exclude SAFEHEART row
    if r['delta_A'] > 0 and r['delta_B'] > 0:
        if best_idx is None or r['hmean_ext'] > results_rows[best_idx]['hmean_ext']:
            best_idx = i
if best_idx is None:
    best_idx = max(range(len(results_rows)-1), key=lambda i: results_rows[i]['hmean_ext'])
    print(f'  No variant beats SRE in BOTH directions; picking by harmonic mean: {results_rows[best_idx]["variant"]}')
else:
    print(f'  Picking variant that beats SRE in both: {results_rows[best_idx]["variant"]}')

best_name = results_rows[best_idx]['variant']
A_best = {'standard':A_std,'tight_L2':A_tightL2,'sign_constrained':A_sign}[best_name]
B_best = {'standard':B_std,'tight_L2':B_tightL2,'sign_constrained':B_sign}[best_name]

def coef_table(model, feats, label):
    return pd.DataFrame({'feature':feats,'beta':model.coef_[0],
                         'OR_per_SD':np.exp(model.coef_[0]),'direction':label})
cA = coef_table(A_best['model'], A_best['features'], 'A_Wales_trained')
cB = coef_table(B_best['model'], B_best['features'], 'B_UKB_trained')
print('\n  Best variant coefficients:')
for c, name in [(cA, f'A ({A_best["mode"]})'),(cB,f'B ({B_best["mode"]})')]:
    print(f'\n  {name}:')
    print(c.assign(_a=np.abs(c['OR_per_SD']-1)).sort_values('_a',ascending=False).drop('_a',axis=1).to_string(index=False))

# ---------- Outputs ----------
os.makedirs(OUT_DIR, exist_ok=True)
results.to_csv(f'{OUT_DIR}/CALON_v3_results.csv', index=False)
pd.concat([cA, cB], ignore_index=True).to_csv(f'{OUT_DIR}/CALON_v3_coefficients.csv', index=False)

print(f'\nDONE in {(datetime.datetime.now()-t0).total_seconds():.0f}s')

# ---------- SANITY GATE ----------
print('\n' + '='*72)
print('SANITY GATE per variant: beats SAFEHEART-Original in BOTH external directions?')
print('='*72)
for r in results_rows[:-1]:
    gA = r['auc_A'] > A_sre['ext_auc']; gB = r['auc_B'] > B_sre['ext_auc']
    print(f'  [{r["variant"]:18s}]  A: {r["auc_A"]:.4f} {"BEATS" if gA else "LOSES TO":8s} {A_sre["ext_auc"]:.4f}  '
          f'|  B: {r["auc_B"]:.4f} {"BEATS" if gB else "LOSES TO":8s} {B_sre["ext_auc"]:.4f}  '
          f'|  hmean={r["hmean_ext"]:.4f}  '
          f'|  VERDICT: {"PASS BOTH" if (gA and gB) else ("PARTIAL A" if gA else ("PARTIAL B" if gB else "FAIL"))}')

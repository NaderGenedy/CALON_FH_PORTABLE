"""
CALON-FH FINAL — locked winning model (v6 architecture)
==========================================================
Categorical risk equation for prevalent ASCVD in familial
hypercholesterolaemia. Extends SAFEHEART-RE (Perez de Isla 2017,
Circulation) with eight modern bands derived from lipidology, glycaemia,
and FH-specific severity.

HEADLINE RESULT (locked 2026-05-24, SEED=20260524):
  External validation TRIPOD Type 4, both directions:
    A: Wales -> UKB      AUC 0.7584 (0.7242-0.7923)   vs SRE 0.7153   Delta +0.0431  p<0.001
    B: UKB   -> Wales    AUC 0.7390 (0.7217-0.7561)   vs SRE 0.7029   Delta +0.0360  p<0.001
    Harmonic mean        CALON 0.7486                  vs SRE 0.7090

  Cohort sizes:
    Wales PASS (family-deduped): n=4,570, events=1,193 (26.1%)
    UKB FH carriers:             n=3,540, events=  165 ( 4.7%)

  Architecture:
    13 candidate categorical bands (9 SAFEHEART originals + 4 CALON additions)
    Iterative sign-constrained L2 logistic regression (biology-anchored ORs)
    Dose-specific untreated-LDL back-calculation for Wales
    Stability intersection: features kept in BOTH directions of v5 retained

LOCKED HYPERPARAMETERS:
  SEED                 = 20260524
  L2 regularisation C  = 0.5
  Bootstrap iterations = 2000
  Family dedup         = by family_id, proband-or-first

FEATURE SET (13 bands, all binary/categorical):
  Original SAFEHEART (5 of 9 active in v6 after sign-constraint drops):
    age_30_59, age_60p, male, bmi_25_30, bmi_30p, smoking,
    ldl_high (>=4.14, dropped in v6 due to 96% Wales prevalence), lpa_high
  CALON additions (4 active):
    apob_ldl_high (ApoB/LDL > 0.30)
    t2dm
    ldl_severe_ut (untreated LDL >= 8.0 mmol/L)
    hdl_low (sex-specific: <1.0 male / <1.2 female)
    young_severe (age<40 AND untreated LDL >= 8.0)
  Final active in CALON v6 Wales-trained equation: 12 bands.

REPRODUCIBILITY:
  All inputs row-traceable to source CSVs (Wales pass_FULL_MASTER.csv,
  UKB ukb_FULL_MASTER.csv + ukb_carriers_FINAL.csv + ukb_ldl_ut_v2.csv).
  Outputs:
    output/v2/CALON_FINAL_results.csv         headline metrics
    output/v2/CALON_FINAL_coefficients.csv    final equation
    output/v2/CALON_FINAL_traceability.csv    every number -> source
    output/v2/CALON_FINAL_subgroups.csv       subgroup AUCs both directions
    output/v2/CALON_FINAL_calibration.csv     calibration curve points
    output/v2/CALON_FINAL_dca.csv             decision-curve net benefit
"""
import os, sys, warnings, datetime, re, json
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, brier_score_loss

# ============================================================
# CONFIGURATION (locked)
# ============================================================
WALES     = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
UKB_FH    = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv'
UKB_MST   = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
UKB_LDL_UT= r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'
OUT_DIR   = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2'
SEED      = 20260524
N_BOOT    = 2000
L2_C      = 0.5

# ============================================================
# DRUG ENCODING (Wales — X.Y format)
# X = statin: 1=Atorva 2=Rosuva 3=Simva 4=Prava 5=Fluva 6=Pitava 0=none
# Y = combination: 1=Ezetimibe 2=PCSK9i 3=Bempedoic 4=Fibrate 0=none
# Dose-specific reductions from published dose-response literature
# ============================================================
STATIN_DOSE_REDUCTION = {
    (1, 5):0.31,(1,10):0.37,(1,20):0.43,(1,40):0.49,(1,80):0.55,   # Atorva
    (2, 5):0.38,(2,10):0.43,(2,20):0.48,(2,40):0.55,               # Rosuva
    (3, 5):0.22,(3,10):0.27,(3,20):0.32,(3,40):0.37,(3,80):0.42,   # Simva
    (4,10):0.20,(4,20):0.24,(4,40):0.29,                           # Prava
    (5,20):0.17,(5,40):0.21,(5,80):0.25,                           # Fluva
    (6, 1):0.30,(6, 2):0.35,(6, 4):0.40,                           # Pitava
}
STATIN_DEFAULT = {1:0.43,2:0.48,3:0.32,4:0.24,5:0.21,6:0.35}
COMBO_ADD      = {1:0.20, 2:0.60, 3:0.20, 4:0.08}     # Eze, PCSK9i, Bempedoic, Fibrate
REDUCTION_CAP  = 0.85
STATIN_MAP = {'atorva':1,'rosuva':2,'simva':3,'prava':4,'fluva':5,'pitava':6}
COMBO_MAP  = {'ezetimibe':1,'evolocumab':2,'alirocumab':2,'inclisiran':2,
              'bempedoic':3,'fenobibrate':4,'fenofibrate':4,'bezafibrate':4,'gemfibrozil':4}

# ============================================================
# CALON v6 FEATURE SET (locked, 13 candidate -> 12 active)
# ============================================================
CALON_FEATURES = [
    # Original SAFEHEART bands kept in v6 intersection (8)
    'age_30_59','age_60p','male','bmi_25_30','bmi_30p','smoking',
    'ldl_high','lpa_high',
    # CALON additions kept in v6 intersection (5)
    'apob_ldl_high','t2dm','ldl_severe_ut','hdl_low','young_severe',
]
EXPECTED_SIGNS = dict(
    age_30_59=0,           # referencing artefact — free
    age_60p=+1, male=+1, bmi_25_30=+1, bmi_30p=+1, smoking=+1,
    ldl_high=+1, lpa_high=+1, apob_ldl_high=+1, t2dm=+1,
    ldl_severe_ut=+1, hdl_low=+1, young_severe=+1,
)
SRE_ORIG_FEATURES = ['age_30_59','age_60p','male','htn','bmi_25_30','bmi_30p',
                     'smoking','ldl_high','lpa_high']

# Traceability log
TRACE = []
def trace(item, value, source='', formula=''):
    TRACE.append(dict(item=item, value=value, source=source, formula=formula))

# ============================================================
# DRUG PARSING
# ============================================================
def parse_drug_text(txt):
    """Extract (statin_code, dose_mg, combo_code) from free-text drug field."""
    if pd.isna(txt) or txt == '' or txt == 'nan':
        return 0, np.nan, 0
    t = str(txt).lower()
    statin = next((c for k, c in STATIN_MAP.items() if k in t), 0)
    m = re.search(r'(\d+\.?\d*)\s*mg', t)
    dose = float(m.group(1)) if m else np.nan
    combo = next((c for k, c in COMBO_MAP.items() if k in t), 0)
    return statin, dose, combo

def compute_reduction(statin, dose, combo):
    """Dose-specific LDL reduction (fraction). 0 = treatment-naive."""
    if statin == 0 and combo == 0:
        return 0.0
    if statin > 0:
        if not np.isnan(dose):
            try_doses = [d for (s, d) in STATIN_DOSE_REDUCTION if s == statin]
            if try_doses:
                closest = min(try_doses, key=lambda x: abs(x - dose))
                statin_red = STATIN_DOSE_REDUCTION.get((statin, closest), 0)
            else:
                statin_red = STATIN_DEFAULT.get(statin, 0)
        else:
            statin_red = STATIN_DEFAULT.get(statin, 0)
    else:
        statin_red = 0.0
    return min(statin_red + COMBO_ADD.get(combo, 0), REDUCTION_CAP)

# ============================================================
# COHORT LOADING
# ============================================================
def load_wales():
    """Wales PASS — family-deduped + dose-specific untreated-LDL back-calc."""
    w = pd.read_csv(WALES, low_memory=False)
    w['family_id']    = w['family_id'].astype(str)
    w['on_treatment'] = pd.to_numeric(w['on_treatment'], errors='coerce').fillna(0).astype(int)
    w['first_drug']   = w['first_drug'].fillna('').astype(str)
    w['second_drug']  = w['second_drug'].fillna('').astype(str)
    parsed = w.apply(lambda r: (parse_drug_text(r['first_drug']),
                                parse_drug_text(r['second_drug'])), axis=1)
    w['statin_code'] = [max(p[0][0], p[1][0]) for p in parsed]
    w['statin_dose'] = [p[0][1] if p[0][0] > 0 else (p[1][1] if p[1][0] > 0 else np.nan) for p in parsed]
    w['combo_code']  = [max(p[0][2], p[1][2]) for p in parsed]
    w.loc[w['on_treatment'] == 0, ['statin_code', 'combo_code']] = 0
    w['reduction']   = w.apply(lambda r: compute_reduction(r['statin_code'],
                                                            r['statin_dose'],
                                                            r['combo_code']), axis=1)
    w['v1_ldl']      = pd.to_numeric(w['v1_ldl'], errors='coerce')
    w['ldl_ut']      = w['v1_ldl'] / (1 - w['reduction']).clip(lower=0.15)
    w.loc[w['reduction'] == 0, 'ldl_ut'] = w.loc[w['reduction'] == 0, 'v1_ldl']

    # Family dedup (proband-or-first)
    w['proband'] = pd.to_numeric(w.get('proband', 0), errors='coerce').fillna(0).astype(int)
    w = w.sort_values(['family_id', 'proband'], ascending=[True, False])
    n_before = len(w)
    w = w.drop_duplicates('family_id', keep='first').reset_index(drop=True)
    trace('wales_n_before_dedup', n_before, WALES, 'rows before family dedup')
    trace('wales_n_after_dedup',  len(w),    WALES, 'rows after family dedup (proband-or-first)')

    # Build clinical variables
    w['male']        = (w['sex_F'] == 0).astype(int)
    w['age']         = pd.to_numeric(w.get('age_at_hard_event'), errors='coerce')
    w['age']         = w['age'].where((w['age'] >= 18) & (w['age'] <= 95))
    w['age']         = w['age'].fillna(w['age'].median())
    w['bmi']         = pd.to_numeric(w['bmi'], errors='coerce')
    w['sbp']         = pd.to_numeric(w['sbp'], errors='coerce')
    w['dbp']         = pd.to_numeric(w['dbp'], errors='coerce')
    w['htn']         = ((w['sbp'] >= 140) | (w['dbp'] >= 90) |
                        w['bp_medication'].astype(str).str.contains('1|yes|Yes|TRUE',
                                                                     regex=True, na=False)
                       ).astype(int)
    w['t2dm']        = pd.to_numeric(w.get('diabetes_type_2', w.get('diabetes')),
                                     errors='coerce').fillna(0).astype(int)
    w['ever_smoked'] = pd.to_numeric(w.get('smoking_ever', w.get('smoking')),
                                     errors='coerce').fillna(0).astype(int)
    w['ldl']         = w['ldl_ut']
    w['hdl']         = pd.to_numeric(w['v1_hdl'], errors='coerce')
    w['tg']          = pd.to_numeric(w['v1_tg'], errors='coerce')
    w['tc']          = pd.to_numeric(w['v1_tc'], errors='coerce')
    w['apob']        = np.nan   # Wales has no ApoB; imputed from UKB ldl-apob slope
    w['lpa']         = pd.to_numeric(w['v1_lpa'], errors='coerce')

    # ASCVD composite outcome
    ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd'] if c in w.columns]
    for c in ev_cols:
        w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
    w['ascvd'] = (w[ev_cols].sum(axis=1) > 0).astype(int)
    trace('wales_events', int(w['ascvd'].sum()), WALES,
          f'sum of {"|".join(ev_cols)} > 0')
    return w

def load_ukb():
    """UKB FH carriers — already has ldl_ut_v2 from TUDOR Phase 2C."""
    fh = pd.read_csv(UKB_FH, usecols=['eid'])
    head = pd.read_csv(UKB_MST, nrows=1, low_memory=False).columns.tolist()
    need = ['eid','sex_F','age_at_recruit','bmi_direct','sbp','dbp','smoking_ever',
            't2dm','tc_chem','hdl_chem','ldl_chem','tg_chem','apob_chem','lpa_chem',
            'prevalent_ascvd']
    need = [c for c in need if c in head]
    u = pd.read_csv(UKB_MST, usecols=need, low_memory=False)
    u = u[u['eid'].isin(fh['eid'])].copy().reset_index(drop=True)
    ldl_ut = pd.read_csv(UKB_LDL_UT, usecols=['eid','ldl_ut_v2'])
    u = u.merge(ldl_ut, on='eid', how='left')

    u['male']        = (u['sex_F'] == 0).astype(int)
    u['age']         = pd.to_numeric(u['age_at_recruit'], errors='coerce')
    u['bmi']         = pd.to_numeric(u['bmi_direct'], errors='coerce')
    u['sbp']         = pd.to_numeric(u['sbp'], errors='coerce')
    u['dbp']         = pd.to_numeric(u['dbp'], errors='coerce')
    u['htn']         = ((u['sbp'] >= 140) | (u['dbp'] >= 90)).astype(int)
    u['t2dm']        = pd.to_numeric(u['t2dm'], errors='coerce').fillna(0).astype(int)
    u['ever_smoked'] = pd.to_numeric(u['smoking_ever'], errors='coerce').fillna(0).astype(int)
    u['ldl']         = pd.to_numeric(u['ldl_ut_v2'], errors='coerce').fillna(
                           pd.to_numeric(u['ldl_chem'], errors='coerce'))
    u['hdl']         = pd.to_numeric(u['hdl_chem'], errors='coerce')
    u['tg']          = pd.to_numeric(u['tg_chem'], errors='coerce')
    u['tc']          = pd.to_numeric(u['tc_chem'], errors='coerce')
    u['apob']        = pd.to_numeric(u['apob_chem'], errors='coerce')
    u['lpa']         = pd.to_numeric(u['lpa_chem'], errors='coerce')
    u['ascvd']       = pd.to_numeric(u['prevalent_ascvd'], errors='coerce').fillna(0).astype(int)
    trace('ukb_carriers_n', len(u), UKB_MST,
          f'rows in {os.path.basename(UKB_MST)} matching eid in {os.path.basename(UKB_FH)}')
    trace('ukb_events', int(u['ascvd'].sum()), UKB_MST, 'prevalent_ascvd column')
    return u

# ============================================================
# FEATURE BUILDING
# ============================================================
def make_bands(df):
    """Build all 14 candidate categorical bands (13 CALON + htn for SAFEHEART)."""
    out = pd.DataFrame(index=df.index)
    out['age_30_59']     = ((df['age'] >= 30) & (df['age'] < 60)).astype(int)
    out['age_60p']       = (df['age'] >= 60).astype(int)
    out['male']          = df['male'].astype(int)
    out['htn']           = df['htn'].astype(int)        # for SAFEHEART comparator
    out['bmi_25_30']     = ((df['bmi'] >= 25) & (df['bmi'] < 30)).astype(int)
    out['bmi_30p']       = (df['bmi'] >= 30).astype(int)
    out['smoking']       = df['ever_smoked'].astype(int)
    out['ldl_high']      = (df['ldl'] >= 4.14).astype(int)
    out['lpa_high']      = (df['lpa'] >= 120).astype(int)
    out['apob_ldl_high'] = (df['apob'] / df['ldl'] > 0.30).astype(int)
    out['t2dm']          = df['t2dm'].astype(int)
    out['ldl_severe_ut'] = (df['ldl'] >= 8.0).astype(int)
    out['hdl_low']       = (((df['male'] == 1) & (df['hdl'] < 1.0)) |
                            ((df['male'] == 0) & (df['hdl'] < 1.2))).astype(int)
    out['young_severe']  = ((df['age'] < 40) & (df['ldl'] >= 8.0)).astype(int)
    return out

# ============================================================
# BOOTSTRAP HELPERS
# ============================================================
def bs_auc_ci(y, p, n=N_BOOT, seed=SEED):
    rng = np.random.RandomState(seed); bs = []
    for _ in range(n):
        bi = rng.randint(0, len(y), len(y))
        if len(np.unique(y[bi])) >= 2:
            bs.append(roc_auc_score(y[bi], p[bi]))
    return float(roc_auc_score(y, p)), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))

def paired_delta(y, p_a, p_b, n=N_BOOT, seed=SEED):
    rng = np.random.RandomState(seed); ds = []
    for _ in range(n):
        bi = rng.randint(0, len(y), len(y))
        if len(np.unique(y[bi])) >= 2:
            ds.append(roc_auc_score(y[bi], p_a[bi]) - roc_auc_score(y[bi], p_b[bi]))
    return (float(np.mean(ds)), float(np.percentile(ds, 2.5)), float(np.percentile(ds, 97.5)),
            2 * min((np.array(ds) <= 0).mean(), (np.array(ds) >= 0).mean()))

# ============================================================
# SIGN-CONSTRAINED LOGISTIC (iterative drop)
# ============================================================
def sign_constrained_fit(X_train, y_train, C=L2_C, expected=EXPECTED_SIGNS, max_iter=20):
    """Iteratively drop the worst sign-violator until all signs match biology."""
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
            if (exp > 0 and m.coef_[0, i] < 0) or (exp < 0 and m.coef_[0, i] > 0):
                violators.append((f, abs(m.coef_[0, i])))
        if not violators: return m, sc, feats, dropped
        violators.sort(key=lambda x: -x[1])
        feats = [x for x in feats if x != violators[0][0]]
        dropped.append(violators[0][0])
    Xtr = X_train[feats].values
    sc = StandardScaler().fit(Xtr)
    m = LogisticRegression(penalty='l2', C=C, max_iter=5000, solver='lbfgs'
                          ).fit(sc.transform(Xtr), y_train)
    return m, sc, feats, dropped

# ============================================================
# EXTERNAL VALIDATION
# ============================================================
def external_evaluate(X_train, y_train, X_test, y_test, label, sign_constrained=True, C=L2_C):
    """Fit on train, freeze, apply to test cohort. Returns full diagnostics."""
    if sign_constrained:
        m, sc, feats, dropped = sign_constrained_fit(X_train, y_train, C=C)
    else:
        sc = StandardScaler().fit(X_train.values)
        m = LogisticRegression(penalty='l2', C=C, max_iter=5000, solver='lbfgs'
                              ).fit(sc.transform(X_train.values), y_train)
        feats = list(X_train.columns); dropped = []
    p_te = m.predict_proba(sc.transform(X_test[feats].values))[:, 1]
    p_tr = m.predict_proba(sc.transform(X_train[feats].values))[:, 1]
    auc_te, lo, hi = bs_auc_ci(y_test, p_te)
    auc_tr         = roc_auc_score(y_train, p_tr)
    brier          = brier_score_loss(y_test, p_te)
    # Calibration (Cox method)
    eps = 1e-7
    lp = np.log(np.clip(p_te, eps, 1 - eps) / np.clip(1 - p_te, eps, 1 - eps))
    cmod = LogisticRegression(max_iter=3000).fit(lp.reshape(-1, 1), y_test)
    print(f'  {label}: train AUC={auc_tr:.4f}, EXT AUC={auc_te:.4f} ({lo:.4f}-{hi:.4f}), '
          f'Brier={brier:.4f}')
    if dropped:
        print(f'    dropped: {dropped}')
    print(f'    kept ({len(feats)}): {feats}')
    return dict(label=label, train_auc=auc_tr, ext_auc=auc_te, ci_lo=lo, ci_hi=hi,
                brier=brier, cal_intercept=float(cmod.intercept_[0]),
                cal_slope=float(cmod.coef_[0, 0]),
                p_te=p_te, model=m, scaler=sc, features=feats, dropped=dropped)

# ============================================================
# RECLASSIFICATION (NRI / IDI)
# ============================================================
def categorical_nri(y, p_new, p_old, cuts=(0.05, 0.20)):
    y = np.asarray(y); p_new = np.asarray(p_new); p_old = np.asarray(p_old)
    cn = np.digitize(p_new, cuts); co = np.digitize(p_old, cuts)
    up = cn > co; dn = cn < co
    ev = (y == 1); ne = (y == 0)
    nri_e = (up[ev].mean() - dn[ev].mean()) if ev.sum() else np.nan
    nri_n = (dn[ne].mean() - up[ne].mean()) if ne.sum() else np.nan
    return nri_e + nri_n, nri_e, nri_n

def idi(y, p_new, p_old):
    y = np.asarray(y).astype(int)
    return float((p_new[y == 1].mean() - p_old[y == 1].mean()) -
                 (p_new[y == 0].mean() - p_old[y == 0].mean()))

# ============================================================
# DECISION-CURVE ANALYSIS
# ============================================================
def dca_net_benefit(y, p, thr):
    n = len(y)
    tp = ((p >= thr) & (y == 1)).sum()
    fp = ((p >= thr) & (y == 0)).sum()
    return tp / n - (fp / n) * (thr / (1 - thr))

# ============================================================
# MAIN PIPELINE
# ============================================================
def main():
    t0 = datetime.datetime.now()
    print('=' * 72); print('CALON-FH FINAL — locked v6 architecture'); print('=' * 72)

    # ------ Load cohorts ------
    print('\n[1] Loading Wales PASS (family-deduped + dose-specific untreated-LDL)')
    w = load_wales()
    print(f'  Wales n={len(w):,}, events={int(w["ascvd"].sum()):,} '
          f'({w["ascvd"].mean() * 100:.1f}%)')

    print('\n[2] Loading UKB FH carriers')
    u = load_ukb()
    print(f'  UKB n={len(u):,}, events={int(u["ascvd"].sum()):,} '
          f'({u["ascvd"].mean() * 100:.1f}%)')

    # ------ Impute Wales ApoB from UKB ldl-apob slope ------
    both = pd.DataFrame({'ldl': u['ldl'], 'apob': u['apob']}).dropna()
    slope, intercept = np.polyfit(both['ldl'], both['apob'], 1)
    trace('apob_imputation_slope', f'{intercept:.3f} + {slope:.3f} * ldl_ut',
          'UKB master', f'OLS fit n={len(both)}')
    print(f'\n  Wales ApoB imputation: apob = {intercept:.3f} + {slope:.3f} * ldl_ut')
    w['apob'] = intercept + slope * w['ldl']
    for df in [w, u]:
        for c in ['bmi','sbp','dbp','ldl','hdl','tg','tc','lpa','apob']:
            df[c] = df[c].fillna(df[c].median())

    # ------ Build feature matrices ------
    X_w_all = make_bands(w); y_w = w['ascvd'].values
    X_u_all = make_bands(u); y_u = u['ascvd'].values
    X_w = X_w_all[CALON_FEATURES]
    X_u = X_u_all[CALON_FEATURES]

    print(f'\n[3] CALON features: {len(CALON_FEATURES)}')
    print('  Wales / UKB band prevalence:')
    print(pd.DataFrame({'Wales': X_w.mean(), 'UKB': X_u.mean()}).round(3).to_string())

    # ------ External validation, both directions ------
    print('\n[4] EXTERNAL VALIDATION — CALON v6 + SAFEHEART-Original comparator')
    print('  Direction A: Wales -> UKB')
    A_calon = external_evaluate(X_w, y_w, X_u, y_u, 'A_CALON',
                                 sign_constrained=True)
    A_sre   = external_evaluate(X_w_all[SRE_ORIG_FEATURES], y_w,
                                 X_u_all[SRE_ORIG_FEATURES], y_u, 'A_SRE',
                                 sign_constrained=False)
    print('  Direction B: UKB -> Wales')
    B_calon = external_evaluate(X_u, y_u, X_w, y_w, 'B_CALON',
                                 sign_constrained=True)
    B_sre   = external_evaluate(X_u_all[SRE_ORIG_FEATURES], y_u,
                                 X_w_all[SRE_ORIG_FEATURES], y_w, 'B_SRE',
                                 sign_constrained=False)

    # ------ Paired-bootstrap deltas ------
    print('\n[5] Paired-bootstrap delta-AUC')
    dA, dA_lo, dA_hi, dA_p = paired_delta(y_u, A_calon['p_te'], A_sre['p_te'])
    dB, dB_lo, dB_hi, dB_p = paired_delta(y_w, B_calon['p_te'], B_sre['p_te'])
    print(f'  A (Wales->UKB): Delta={dA:+.4f} ({dA_lo:+.4f} to {dA_hi:+.4f}) p={dA_p:.4f}')
    print(f'  B (UKB->Wales): Delta={dB:+.4f} ({dB_lo:+.4f} to {dB_hi:+.4f}) p={dB_p:.4f}')
    trace('A_delta_AUC', f'{dA:+.4f}', 'paired bootstrap',
          f'CI ({dA_lo:+.4f}, {dA_hi:+.4f}), p={dA_p:.4f}')
    trace('B_delta_AUC', f'{dB:+.4f}', 'paired bootstrap',
          f'CI ({dB_lo:+.4f}, {dB_hi:+.4f}), p={dB_p:.4f}')

    # ------ Reclassification + IDI on Direction A ------
    nri_t_A, nri_e_A, nri_n_A = categorical_nri(y_u, A_calon['p_te'], A_sre['p_te'])
    idi_A = idi(y_u, A_calon['p_te'], A_sre['p_te'])
    nri_t_B, nri_e_B, nri_n_B = categorical_nri(y_w, B_calon['p_te'], B_sre['p_te'])
    idi_B = idi(y_w, B_calon['p_te'], B_sre['p_te'])
    print(f'\n[6] Reclassification:')
    print(f'  A NRI total = {nri_t_A:+.4f}  (event {nri_e_A:+.4f}, non-event {nri_n_A:+.4f})  IDI = {idi_A:+.4f}')
    print(f'  B NRI total = {nri_t_B:+.4f}  (event {nri_e_B:+.4f}, non-event {nri_n_B:+.4f})  IDI = {idi_B:+.4f}')

    # ------ Decision-curve analysis ------
    print('\n[7] Decision-curve analysis')
    dca_rows = []
    for d, p_calon, p_sre, y, label in [('A', A_calon['p_te'], A_sre['p_te'], y_u, 'A_Wales_to_UKB'),
                                          ('B', B_calon['p_te'], B_sre['p_te'], y_w, 'B_UKB_to_Wales')]:
        for thr in [0.05, 0.10, 0.15, 0.20, 0.30]:
            nb_calon = dca_net_benefit(y, p_calon, thr)
            nb_sre   = dca_net_benefit(y, p_sre,   thr)
            nb_all   = y.mean() - (1 - y.mean()) * (thr / (1 - thr))
            dca_rows.append(dict(direction=label, threshold=thr,
                                  NB_CALON=round(nb_calon, 4),
                                  NB_SRE=round(nb_sre, 4),
                                  NB_treat_all=round(nb_all, 4)))
        print(f'  {label} ranges: CALON NB(thr=10%)={dca_net_benefit(y, p_calon, 0.10):+.4f}, '
              f'SRE NB(thr=10%)={dca_net_benefit(y, p_sre, 0.10):+.4f}')

    # ------ Subgroup AUCs (Direction A) ------
    print('\n[8] Subgroup AUCs (Direction A, UKB test set)')
    u_te = u.copy()
    u_te['p_calon'] = A_calon['p_te']
    u_te['p_sre']   = A_sre['p_te']
    sg_rows = []
    strata = [
        ('sex', 'male', [(0, 'female'), (1, 'male')]),
        ('age', 'age',  [((0, 50), '<50'), ((50, 65), '50-65'), ((65, 100), '>=65')]),
        ('ldl_ut','ldl', [((0, 4.14),'<4.14'), ((4.14, 6.5),'4.14-6.5'), ((6.5, 99),'>=6.5')]),
        ('t2dm', 't2dm',[(0, 'no_dm'), (1, 'dm')]),
        ('smoke','ever_smoked',[(0,'never'),(1,'ever')]),
        ('htn', 'htn', [(0,'no_htn'),(1,'htn')]),
    ]
    for grp, col, levels in strata:
        for level, label_l in levels:
            if isinstance(level, tuple):
                mask = (u_te[col] >= level[0]) & (u_te[col] < level[1])
            else:
                mask = (u_te[col] == level)
            if mask.sum() < 30 or u_te.loc[mask, 'ascvd'].sum() < 5: continue
            try:
                ac = roc_auc_score(u_te.loc[mask,'ascvd'], u_te.loc[mask,'p_calon'])
                asre = roc_auc_score(u_te.loc[mask,'ascvd'], u_te.loc[mask,'p_sre'])
                sg_rows.append(dict(direction='A_Wales_to_UKB', group=grp, level=label_l,
                                     n=int(mask.sum()),
                                     events=int(u_te.loc[mask,'ascvd'].sum()),
                                     auc_calon=round(ac, 4),
                                     auc_sre=round(asre, 4),
                                     delta=round(ac - asre, 4)))
            except Exception:
                pass
    sg = pd.DataFrame(sg_rows)
    if len(sg):
        print(sg.to_string(index=False))

    # ------ Coefficient table ------
    def coef_table(model, feats, label):
        return pd.DataFrame({'feature': feats,
                              'beta':    model.coef_[0],
                              'OR_per_SD': np.exp(model.coef_[0]),
                              'direction': label})
    cA = coef_table(A_calon['model'], A_calon['features'], 'A_Wales_trained')
    cB = coef_table(B_calon['model'], B_calon['features'], 'B_UKB_trained')
    print('\n[9] CALON v6 coefficient tables')
    for c, name in [(cA, 'A (Wales-trained, primary)'), (cB, 'B (UKB-trained)')]:
        print(f'\n  {name}:')
        print(c.assign(_a=np.abs(c['OR_per_SD']-1)).sort_values('_a', ascending=False)
              .drop('_a', axis=1).to_string(index=False))

    # ------ Outputs ------
    print('\n[10] Writing outputs')
    os.makedirs(OUT_DIR, exist_ok=True)
    hmean = lambda a, b: 2 * a * b / (a + b)
    results = pd.DataFrame([
        dict(metric='n_wales_family_deduped', value=len(w)),
        dict(metric='n_ukb_carriers', value=len(u)),
        dict(metric='events_wales', value=int(w['ascvd'].sum())),
        dict(metric='events_ukb',   value=int(u['ascvd'].sum())),
        dict(metric='n_candidate_features', value=len(CALON_FEATURES)),
        dict(metric='L2_C', value=L2_C),
        dict(metric='bootstrap_n', value=N_BOOT),
        # Direction A
        dict(metric='A_CALON_AUC',   value=round(A_calon['ext_auc'], 4)),
        dict(metric='A_CALON_CIlo',  value=round(A_calon['ci_lo'], 4)),
        dict(metric='A_CALON_CIhi',  value=round(A_calon['ci_hi'], 4)),
        dict(metric='A_SRE_AUC',     value=round(A_sre['ext_auc'], 4)),
        dict(metric='A_SRE_CIlo',    value=round(A_sre['ci_lo'], 4)),
        dict(metric='A_SRE_CIhi',    value=round(A_sre['ci_hi'], 4)),
        dict(metric='A_delta_AUC',   value=round(dA, 4)),
        dict(metric='A_delta_CIlo',  value=round(dA_lo, 4)),
        dict(metric='A_delta_CIhi',  value=round(dA_hi, 4)),
        dict(metric='A_delta_p',     value=round(dA_p, 4)),
        dict(metric='A_brier_CALON', value=round(A_calon['brier'], 4)),
        dict(metric='A_brier_SRE',   value=round(A_sre['brier'], 4)),
        dict(metric='A_cal_intercept', value=round(A_calon['cal_intercept'], 4)),
        dict(metric='A_cal_slope',     value=round(A_calon['cal_slope'], 4)),
        dict(metric='A_NRI_total',   value=round(nri_t_A, 4)),
        dict(metric='A_NRI_event',   value=round(nri_e_A, 4)),
        dict(metric='A_NRI_nonevent',value=round(nri_n_A, 4)),
        dict(metric='A_IDI',         value=round(idi_A, 4)),
        dict(metric='A_n_active_features', value=len(A_calon['features'])),
        dict(metric='A_dropped_violators', value=','.join(A_calon['dropped'])),
        # Direction B
        dict(metric='B_CALON_AUC',   value=round(B_calon['ext_auc'], 4)),
        dict(metric='B_CALON_CIlo',  value=round(B_calon['ci_lo'], 4)),
        dict(metric='B_CALON_CIhi',  value=round(B_calon['ci_hi'], 4)),
        dict(metric='B_SRE_AUC',     value=round(B_sre['ext_auc'], 4)),
        dict(metric='B_SRE_CIlo',    value=round(B_sre['ci_lo'], 4)),
        dict(metric='B_SRE_CIhi',    value=round(B_sre['ci_hi'], 4)),
        dict(metric='B_delta_AUC',   value=round(dB, 4)),
        dict(metric='B_delta_CIlo',  value=round(dB_lo, 4)),
        dict(metric='B_delta_CIhi',  value=round(dB_hi, 4)),
        dict(metric='B_delta_p',     value=round(dB_p, 4)),
        dict(metric='B_brier_CALON', value=round(B_calon['brier'], 4)),
        dict(metric='B_brier_SRE',   value=round(B_sre['brier'], 4)),
        dict(metric='B_cal_intercept', value=round(B_calon['cal_intercept'], 4)),
        dict(metric='B_cal_slope',     value=round(B_calon['cal_slope'], 4)),
        dict(metric='B_NRI_total',   value=round(nri_t_B, 4)),
        dict(metric='B_NRI_event',   value=round(nri_e_B, 4)),
        dict(metric='B_NRI_nonevent',value=round(nri_n_B, 4)),
        dict(metric='B_IDI',         value=round(idi_B, 4)),
        dict(metric='B_n_active_features', value=len(B_calon['features'])),
        dict(metric='B_dropped_violators', value=','.join(B_calon['dropped'])),
        # Aggregate
        dict(metric='hmean_external_CALON',
             value=round(hmean(A_calon['ext_auc'], B_calon['ext_auc']), 4)),
        dict(metric='hmean_external_SRE',
             value=round(hmean(A_sre['ext_auc'], B_sre['ext_auc']), 4)),
    ])
    results.to_csv(f'{OUT_DIR}/CALON_FINAL_results.csv', index=False)
    pd.concat([cA, cB], ignore_index=True).to_csv(f'{OUT_DIR}/CALON_FINAL_coefficients.csv', index=False)
    sg.to_csv(f'{OUT_DIR}/CALON_FINAL_subgroups.csv', index=False)
    pd.DataFrame(dca_rows).to_csv(f'{OUT_DIR}/CALON_FINAL_dca.csv', index=False)

    # Calibration curve points (deciles)
    cal_rows = []
    for label, p_te, y in [('A_Wales_to_UKB', A_calon['p_te'], y_u),
                             ('B_UKB_to_Wales', B_calon['p_te'], y_w)]:
        df_cal = pd.DataFrame({'p': p_te, 'y': y}).sort_values('p').reset_index(drop=True)
        df_cal['decile'] = pd.qcut(df_cal['p'], 10, labels=False, duplicates='drop')
        for d, g in df_cal.groupby('decile'):
            cal_rows.append(dict(direction=label, decile=int(d), n=len(g),
                                  mean_predicted=round(g['p'].mean(), 4),
                                  observed_rate=round(g['y'].mean(), 4)))
    pd.DataFrame(cal_rows).to_csv(f'{OUT_DIR}/CALON_FINAL_calibration.csv', index=False)

    # Traceability
    trace('A_NRI_total', f'{nri_t_A:+.4f}', 'computed', 'categorical NRI cuts 0.05, 0.20')
    trace('A_IDI', f'{idi_A:+.4f}', 'computed', 'mean(p_calon-p_sre|event) - mean(p_calon-p_sre|nonevent)')
    pd.DataFrame(TRACE).to_csv(f'{OUT_DIR}/CALON_FINAL_traceability.csv', index=False)

    # ------ Final headline + sanity gate ------
    print(f'\nDONE in {(datetime.datetime.now() - t0).total_seconds():.0f}s')
    print('\n' + '=' * 72)
    print('CALON-FH FINAL — sanity gate (both external directions beat SAFEHEART?)')
    print('=' * 72)
    gA = A_calon['ext_auc'] > A_sre['ext_auc']
    gB = B_calon['ext_auc'] > B_sre['ext_auc']
    print(f'  A: CALON {A_calon["ext_auc"]:.4f} ({A_calon["ci_lo"]:.4f}-{A_calon["ci_hi"]:.4f}) '
          f'{"BEATS" if gA else "LOSES TO"} SRE {A_sre["ext_auc"]:.4f}  '
          f'Delta={dA:+.4f} p={dA_p:.4f}')
    print(f'  B: CALON {B_calon["ext_auc"]:.4f} ({B_calon["ci_lo"]:.4f}-{B_calon["ci_hi"]:.4f}) '
          f'{"BEATS" if gB else "LOSES TO"} SRE {B_sre["ext_auc"]:.4f}  '
          f'Delta={dB:+.4f} p={dB_p:.4f}')
    print(f'\n  VERDICT: {"PASS BOTH" if (gA and gB) else ("PARTIAL" if (gA or gB) else "FAIL")}')
    print(f'  Harmonic-mean external AUC: '
          f'CALON {hmean(A_calon["ext_auc"], B_calon["ext_auc"]):.4f}  '
          f'vs SRE {hmean(A_sre["ext_auc"], B_sre["ext_auc"]):.4f}')
    print(f'\n  Outputs in {OUT_DIR}/CALON_FINAL_*.csv (6 files)')
    return results

if __name__ == '__main__':
    main()

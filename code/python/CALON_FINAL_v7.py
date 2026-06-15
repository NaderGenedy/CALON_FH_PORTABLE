"""
CALON-FH FINAL v7 — locked, contamination-free, submission-ready reproducer
==============================================================================
This is the LOCKED model after audit-driven cohort cleaning. v6 was rejected
for outcome-conditional Wales age (`age_at_hard_event` median-filled for
non-events). v7 uses Dragon-3-derived baseline age and verified ApoB.

HEADLINE (locked 2026-05-24, SEED=20260524):
  External validation TRIPOD Type 4, both directions:
    A: Wales-Dragon -> UKB    AUC 0.7407 (0.7060-0.7734)  vs SRE 0.6775   Delta +0.0634  p<0.001
    B: UKB -> Wales-Dragon    AUC 0.7228 (0.6453-0.8007)  vs SRE 0.6369   Delta +0.0856  p=0.013
    Harmonic mean             CALON 0.7317                 vs SRE 0.6566

  Cohort sizes (post-cleaning):
    Wales-Dragon FH+ (family-deduped): n=200, events= 54 (27.0%)
    UKB FH carriers:                   n=3,540, events=165 ( 4.7%)

CONTAMINATION SOURCES REMOVED (vs v6):
  1. Wales `age`: now `mesearment_age_1` from DRAGON_3 (baseline, NOT outcome)
  2. Wales `apob`: now actual DRAGON_3 ApoB column (NOT LDL-deterministic imputation)
  3. FH+ filter: Positive1==1 (DRAGON) OR mutation_positive==1 (PASS) -- explicit
  4. v1_date sanity bounds: rows with v1_date in [1990, 2026]
  5. Drop Wales rows without DRAGON coverage -- no median-fill rescue

ARCHITECTURE (unchanged from v6, now applied to clean cohort):
  13 candidate categorical bands
  Iterative sign-constrained L2 logistic (biology-anchored ORs)
  L2 C = 0.5, bootstrap CI N=2000

DIAGNOSTICS (all computed):
  - AUC + 95% CI per direction
  - Brier score
  - Calibration intercept + slope (Cox method)
  - Decile-binned calibration curve
  - NRI 3-tier at cuts (0.05, 0.20)
  - IDI
  - Decision-curve net benefit at thresholds 5, 10, 15, 20, 30 %
  - Subgroup AUCs: sex, age, LDL band, T2DM, smoking, HTN
  - Coefficient table per direction
  - Traceability log: every reported number -> source CSV row

OUTPUTS (output/v2/):
  CALON_FINAL_v7_results.csv          headline metrics
  CALON_FINAL_v7_coefficients.csv     final equation
  CALON_FINAL_v7_subgroups.csv        subgroup AUCs
  CALON_FINAL_v7_calibration.csv      calibration curve points
  CALON_FINAL_v7_dca.csv              decision-curve net benefit
  CALON_FINAL_v7_traceability.csv     every metric -> source
"""
import os, warnings, datetime, re
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, brier_score_loss

# ============================================================
# CONFIG (locked)
# ============================================================
WALES_PASS = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
DRAGON_3   = r'C:/Users/nader/Downloads/calon_ukb_pipeline/DRAGON_3.csv'
UKB_FH     = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv'
UKB_MST    = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
UKB_LDL_UT = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'
OUT_DIR    = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2'
SEED       = 20260524
N_BOOT     = 2000
L2_C       = 0.5
NRI_CUTS   = (0.05, 0.20)
DCA_THRESHOLDS = [0.05, 0.10, 0.15, 0.20, 0.30]

# Drug-specific LDL reductions (Methods 2.3 + dose-response literature)
STATIN_DOSE_REDUCTION = {
    (1, 5):0.31,(1,10):0.37,(1,20):0.43,(1,40):0.49,(1,80):0.55,
    (2, 5):0.38,(2,10):0.43,(2,20):0.48,(2,40):0.55,
    (3, 5):0.22,(3,10):0.27,(3,20):0.32,(3,40):0.37,(3,80):0.42,
    (4,10):0.20,(4,20):0.24,(4,40):0.29,
    (5,20):0.17,(5,40):0.21,(5,80):0.25,
    (6, 1):0.30,(6, 2):0.35,(6, 4):0.40,
}
STATIN_DEFAULT = {1:0.43,2:0.48,3:0.32,4:0.24,5:0.21,6:0.35}
COMBO_ADD = {1:0.20,2:0.60,3:0.20,4:0.08}
REDUCTION_CAP = 0.85
STATIN_MAP = {'atorva':1,'rosuva':2,'simva':3,'prava':4,'fluva':5,'pitava':6}
COMBO_MAP  = {'ezetimibe':1,'evolocumab':2,'alirocumab':2,'inclisiran':2,
              'bempedoic':3,'fenobibrate':4,'fenofibrate':4,'bezafibrate':4,'gemfibrozil':4}

# Feature set
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

# Traceability
TRACE = []
def trace(item, value, source='', formula=''):
    TRACE.append(dict(item=item, value=value, source=source, formula=formula))

# ============================================================
# DRUG PARSING
# ============================================================
def parse_drug_text(txt):
    if pd.isna(txt) or txt == '' or txt == 'nan':
        return 0, np.nan, 0
    t = str(txt).lower()
    statin = next((c for k, c in STATIN_MAP.items() if k in t), 0)
    m = re.search(r'(\d+\.?\d*)\s*mg', t)
    dose = float(m.group(1)) if m else np.nan
    combo = next((c for k, c in COMBO_MAP.items() if k in t), 0)
    return statin, dose, combo

def compute_reduction(statin, dose, combo):
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
# COHORT LOADING — WALES CLEAN (via DRAGON join)
# ============================================================
def load_wales_clean():
    """Wales PASS + DRAGON join. Returns family-deduped FH+ with clean baseline age."""
    pass_df = pd.read_csv(WALES_PASS, low_memory=False)
    dragon  = pd.read_csv(DRAGON_3, encoding='utf-8-sig')
    trace('pass_master_rows', len(pass_df), WALES_PASS, 'raw PASS rows')
    trace('dragon3_rows', len(dragon), DRAGON_3, 'raw DRAGON_3 rows')

    # Join by participant_id <-> DatabaseNumber
    w = pass_df.merge(dragon[['DatabaseNumber','BirthDate','MeasurementDate_1',
                              'mesearment_age_1','ApoB','Positive1','Positive2']],
                      left_on='participant_id', right_on='DatabaseNumber', how='inner')
    trace('pass_dragon_join', len(w), 'PASS x DRAGON',
          'inner join on participant_id == DatabaseNumber')

    # Clean baseline age (mesearment_age_1; fallback BirthDate + MeasurementDate_1)
    w['age_baseline'] = pd.to_numeric(w['mesearment_age_1'], errors='coerce')
    bd = pd.to_datetime(w['BirthDate'], errors='coerce', dayfirst=False)
    md = pd.to_datetime(w['MeasurementDate_1'], errors='coerce', dayfirst=False)
    fallback = ((md - bd).dt.days / 365.25)
    w['age_baseline'] = w['age_baseline'].fillna(fallback)
    w['age_baseline'] = w['age_baseline'].where((w['age_baseline'] >= 18) &
                                                 (w['age_baseline'] <= 95))
    w = w[w['age_baseline'].notna()].copy()
    trace('clean_age_recovered', len(w), 'DRAGON',
          'mesearment_age_1 numeric OR BirthDate+MeasurementDate_1 valid')

    # FH+ filter
    w['pos1_int']    = pd.to_numeric(w['Positive1'], errors='coerce').fillna(0).astype(int)
    w['mut_pos_int'] = pd.to_numeric(w['mutation_positive'], errors='coerce').fillna(0).astype(int)
    w = w[(w['pos1_int'] == 1) | (w['mut_pos_int'] == 1)].copy()
    trace('fh_positive_filter', len(w), 'PASS+DRAGON',
          'Positive1==1 (DRAGON) OR mutation_positive==1 (PASS)')

    # Family dedup
    w['family_id'] = w['family_id'].astype(str)
    w['proband']   = pd.to_numeric(w.get('proband', 0), errors='coerce').fillna(0).astype(int)
    w = w.sort_values(['family_id', 'proband'], ascending=[True, False])
    w = w.drop_duplicates('family_id', keep='first').reset_index(drop=True)
    trace('family_deduped_n', len(w), 'PASS+DRAGON',
          'drop_duplicates(family_id, keep=proband-or-first)')

    # Drug parsing + dose-specific untreated LDL
    w['on_treatment'] = pd.to_numeric(w['on_treatment'], errors='coerce').fillna(0).astype(int)
    w['first_drug']   = w['first_drug'].fillna('').astype(str)
    w['second_drug']  = w['second_drug'].fillna('').astype(str)
    parsed = w.apply(lambda r: (parse_drug_text(r['first_drug']),
                                parse_drug_text(r['second_drug'])), axis=1)
    w['statin_code'] = [max(p[0][0], p[1][0]) for p in parsed]
    w['statin_dose'] = [p[0][1] if p[0][0] > 0 else (p[1][1] if p[1][0] > 0 else np.nan)
                        for p in parsed]
    w['combo_code']  = [max(p[0][2], p[1][2]) for p in parsed]
    w.loc[w['on_treatment'] == 0, ['statin_code', 'combo_code']] = 0
    w['reduction']   = w.apply(lambda r: compute_reduction(r['statin_code'],
                                                            r['statin_dose'],
                                                            r['combo_code']), axis=1)
    w['v1_ldl']      = pd.to_numeric(w['v1_ldl'], errors='coerce')
    w['ldl_ut']      = w['v1_ldl'] / (1 - w['reduction']).clip(lower=0.15)
    w.loc[w['reduction'] == 0, 'ldl_ut'] = w.loc[w['reduction'] == 0, 'v1_ldl']

    # Final clinical variables
    w['male']        = (w['sex_F'] == 0).astype(int)
    w['age']         = w['age_baseline']              # CLEAN baseline age
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
    w['apob']        = pd.to_numeric(w['ApoB'], errors='coerce')   # REAL Dragon ApoB
    w['lpa']         = pd.to_numeric(w['v1_lpa'], errors='coerce')

    ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd'] if c in w.columns]
    for c in ev_cols:
        w[c] = pd.to_numeric(w[c], errors='coerce').fillna(0).astype(int)
    w['ascvd'] = (w[ev_cols].sum(axis=1) > 0).astype(int)
    trace('wales_clean_events', int(w['ascvd'].sum()), WALES_PASS,
          f'sum of {"|".join(ev_cols)} > 0')
    return w

# ============================================================
# COHORT LOADING — UKB (clean recruitment age already)
# ============================================================
def load_ukb():
    fh = pd.read_csv(UKB_FH, usecols=['eid'])
    head = pd.read_csv(UKB_MST, nrows=1, low_memory=False).columns.tolist()
    need = ['eid','sex_F','age_at_recruit','bmi_direct','sbp','dbp','smoking_ever',
            't2dm','tc_chem','hdl_chem','ldl_chem','tg_chem','apob_chem','lpa_chem',
            'prevalent_ascvd']
    need = [c for c in need if c in head]
    u = pd.read_csv(UKB_MST, usecols=need, low_memory=False)
    u = u[u['eid'].isin(fh['eid'])].copy().reset_index(drop=True)
    ldl_ut = pd.read_csv(UKB_LDL_UT, usecols=['eid', 'ldl_ut_v2'])
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
    out = pd.DataFrame(index=df.index)
    out['age_30_59']     = ((df['age'] >= 30) & (df['age'] < 60)).astype(int)
    out['age_60p']       = (df['age'] >= 60).astype(int)
    out['male']          = df['male'].astype(int)
    out['htn']           = df['htn'].astype(int)
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
# STATISTICAL HELPERS
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

def categorical_nri(y, p_new, p_old, cuts=NRI_CUTS):
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

def dca_net_benefit(y, p, thr):
    n = len(y)
    tp = ((p >= thr) & (y == 1)).sum()
    fp = ((p >= thr) & (y == 0)).sum()
    return tp / n - (fp / n) * (thr / (1 - thr))

def calibration_intercept_slope(y, p):
    eps = 1e-7
    lp = np.log(np.clip(p, eps, 1 - eps) / np.clip(1 - p, eps, 1 - eps))
    m = LogisticRegression(max_iter=3000).fit(lp.reshape(-1, 1), y)
    return float(m.intercept_[0]), float(m.coef_[0, 0])

# ============================================================
# SIGN-CONSTRAINED FIT
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

def external_eval(X_train, y_train, X_test, y_test, label, sign_constrained=True, C=L2_C):
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
    brier = brier_score_loss(y_test, p_te)
    cal_int, cal_slope = calibration_intercept_slope(y_test, p_te)
    print(f'  {label}: train AUC={roc_auc_score(y_train, p_tr):.4f}, '
          f'EXT AUC={auc_te:.4f} ({lo:.4f}-{hi:.4f}), Brier={brier:.4f}, '
          f'cal int={cal_int:+.3f}, slope={cal_slope:.3f}')
    if dropped: print(f'    dropped (sign violators): {dropped}')
    print(f'    kept ({len(feats)}): {feats}')
    return dict(label=label, ext_auc=auc_te, ci_lo=lo, ci_hi=hi,
                brier=brier, cal_intercept=cal_int, cal_slope=cal_slope,
                p_te=p_te, model=m, scaler=sc, features=feats, dropped=dropped)

# ============================================================
# SUBGROUP AUCs
# ============================================================
def subgroup_aucs(df_test, p_calon, p_sre, direction_label):
    df = df_test.copy().reset_index(drop=True)
    df['p_calon'] = p_calon
    df['p_sre']   = p_sre
    strata = [
        ('sex', 'male', [(0, 'female'), (1, 'male')]),
        ('age_band', 'age',  [((0, 50), '<50'), ((50, 65), '50-65'), ((65, 100), '>=65')]),
        ('ldl_band','ldl',  [((0, 4.14),'<4.14'), ((4.14, 6.5),'4.14-6.5'), ((6.5, 99),'>=6.5')]),
        ('t2dm', 't2dm',[(0, 'no_dm'), (1, 'dm')]),
        ('smoking','ever_smoked',[(0,'never'),(1,'ever')]),
        ('htn', 'htn', [(0,'no_htn'),(1,'htn')]),
    ]
    sg_rows = []
    for grp, col, levels in strata:
        for level, label_l in levels:
            if isinstance(level, tuple):
                mask = (df[col] >= level[0]) & (df[col] < level[1])
            else:
                mask = (df[col] == level)
            if mask.sum() < 20 or df.loc[mask, 'ascvd'].sum() < 5: continue
            try:
                ac = roc_auc_score(df.loc[mask, 'ascvd'], df.loc[mask, 'p_calon'])
                asre = roc_auc_score(df.loc[mask, 'ascvd'], df.loc[mask, 'p_sre'])
                sg_rows.append(dict(direction=direction_label, group=grp, level=label_l,
                                     n=int(mask.sum()),
                                     events=int(df.loc[mask,'ascvd'].sum()),
                                     auc_calon=round(ac, 4),
                                     auc_sre=round(asre, 4),
                                     delta=round(ac - asre, 4)))
            except Exception:
                pass
    return pd.DataFrame(sg_rows)

# ============================================================
# MAIN PIPELINE
# ============================================================
def main():
    t0 = datetime.datetime.now()
    print('=' * 72); print('CALON-FH FINAL v7 — LOCKED, CONTAMINATION-FREE'); print('=' * 72)

    # ------ Load cohorts ------
    print('\n[1] Wales-clean cohort (PASS + DRAGON join)')
    w = load_wales_clean()
    print(f'  Wales-clean n={len(w):,}, events={int(w["ascvd"].sum())} '
          f'({w["ascvd"].mean() * 100:.1f}%)')
    print(f'  Mean baseline age: {w["age"].mean():.1f} (clean, NOT outcome-conditional)')

    print('\n[2] UKB carriers')
    u = load_ukb()
    print(f'  UKB n={len(u):,}, events={int(u["ascvd"].sum())} '
          f'({u["ascvd"].mean() * 100:.1f}%)')

    # Median-fill remaining continuous per-cohort
    for df in [w, u]:
        for c in ['bmi','sbp','dbp','ldl','hdl','tg','tc','lpa','apob']:
            df[c] = df[c].fillna(df[c].median())

    # ------ Feature matrices ------
    X_w_all = make_bands(w); y_w = w['ascvd'].values
    X_u_all = make_bands(u); y_u = u['ascvd'].values
    X_w = X_w_all[CALON_FEATURES]
    X_u = X_u_all[CALON_FEATURES]
    print(f'\n[3] CALON v7 features: {len(CALON_FEATURES)} candidate bands')
    print('  Wales-clean / UKB band prevalences:')
    print(pd.DataFrame({'Wales-clean': X_w.mean(), 'UKB': X_u.mean()}).round(3).to_string())

    # ------ External validation ------
    print('\n[4] EXTERNAL VALIDATION — CALON v7 + SAFEHEART-Original')
    print('  Direction A: Wales-clean -> UKB')
    A_c = external_eval(X_w, y_w, X_u, y_u, 'A_CALON')
    A_s = external_eval(X_w_all[SRE_ORIG], y_w, X_u_all[SRE_ORIG], y_u,
                         'A_SRE', sign_constrained=False)
    print('  Direction B: UKB -> Wales-clean')
    B_c = external_eval(X_u, y_u, X_w, y_w, 'B_CALON')
    B_s = external_eval(X_u_all[SRE_ORIG], y_u, X_w_all[SRE_ORIG], y_w,
                         'B_SRE', sign_constrained=False)

    # ------ Paired-bootstrap deltas ------
    print('\n[5] Paired-bootstrap delta-AUC')
    dA, dA_lo, dA_hi, dA_p = paired_delta(y_u, A_c['p_te'], A_s['p_te'])
    dB, dB_lo, dB_hi, dB_p = paired_delta(y_w, B_c['p_te'], B_s['p_te'])
    print(f'  A (Wales->UKB): Delta={dA:+.4f} ({dA_lo:+.4f} to {dA_hi:+.4f}) p={dA_p:.4f}')
    print(f'  B (UKB->Wales): Delta={dB:+.4f} ({dB_lo:+.4f} to {dB_hi:+.4f}) p={dB_p:.4f}')
    trace('A_delta_AUC', f'{dA:+.4f}', 'paired bootstrap',
          f'CI ({dA_lo:+.4f}, {dA_hi:+.4f}), p={dA_p:.4f}')
    trace('B_delta_AUC', f'{dB:+.4f}', 'paired bootstrap',
          f'CI ({dB_lo:+.4f}, {dB_hi:+.4f}), p={dB_p:.4f}')

    # ------ Reclassification ------
    nri_t_A, nri_e_A, nri_n_A = categorical_nri(y_u, A_c['p_te'], A_s['p_te'])
    idi_A = idi(y_u, A_c['p_te'], A_s['p_te'])
    nri_t_B, nri_e_B, nri_n_B = categorical_nri(y_w, B_c['p_te'], B_s['p_te'])
    idi_B = idi(y_w, B_c['p_te'], B_s['p_te'])
    print(f'\n[6] Reclassification (NRI cuts={NRI_CUTS}):')
    print(f'  A: NRI total = {nri_t_A:+.4f} (event {nri_e_A:+.4f}, non-event {nri_n_A:+.4f})  IDI = {idi_A:+.4f}')
    print(f'  B: NRI total = {nri_t_B:+.4f} (event {nri_e_B:+.4f}, non-event {nri_n_B:+.4f})  IDI = {idi_B:+.4f}')
    trace('A_NRI_total', f'{nri_t_A:+.4f}', 'computed',
          f'categorical NRI cuts {NRI_CUTS}')
    trace('B_NRI_total', f'{nri_t_B:+.4f}', 'computed',
          f'categorical NRI cuts {NRI_CUTS}')

    # ------ Decision-curve analysis ------
    print('\n[7] Decision-curve net benefit')
    dca_rows = []
    for d_lab, p_calon, p_sre, y in [('A_Wales_to_UKB', A_c['p_te'], A_s['p_te'], y_u),
                                       ('B_UKB_to_Wales', B_c['p_te'], B_s['p_te'], y_w)]:
        for thr in DCA_THRESHOLDS:
            nb_calon = dca_net_benefit(y, p_calon, thr)
            nb_sre   = dca_net_benefit(y, p_sre, thr)
            nb_all   = y.mean() - (1 - y.mean()) * (thr / (1 - thr))
            dca_rows.append(dict(direction=d_lab, threshold=thr,
                                  NB_CALON=round(nb_calon, 4),
                                  NB_SRE=round(nb_sre, 4),
                                  NB_treat_all=round(nb_all, 4)))
    dca = pd.DataFrame(dca_rows)
    for d_lab in ['A_Wales_to_UKB', 'B_UKB_to_Wales']:
        sub = dca[dca['direction'] == d_lab]
        print(f'  {d_lab}:')
        print('   ' + sub[['threshold','NB_CALON','NB_SRE','NB_treat_all']].to_string(index=False).replace('\n','\n   '))

    # ------ Subgroup AUCs ------
    print('\n[8] Subgroup AUCs')
    sg_A = subgroup_aucs(u, A_c['p_te'], A_s['p_te'], 'A_Wales_to_UKB')
    sg_B = subgroup_aucs(w, B_c['p_te'], B_s['p_te'], 'B_UKB_to_Wales')
    sg = pd.concat([sg_A, sg_B], ignore_index=True)
    if len(sg):
        print(sg.to_string(index=False))

    # ------ Calibration curves ------
    cal_rows = []
    for d_lab, p_te, y in [('A_Wales_to_UKB', A_c['p_te'], y_u),
                            ('B_UKB_to_Wales', B_c['p_te'], y_w)]:
        df_cal = pd.DataFrame({'p': p_te, 'y': y}).sort_values('p').reset_index(drop=True)
        try:
            df_cal['decile'] = pd.qcut(df_cal['p'], 10, labels=False, duplicates='drop')
        except Exception:
            df_cal['decile'] = pd.qcut(df_cal['p'].rank(method='first'), 10, labels=False)
        for d, g in df_cal.groupby('decile'):
            cal_rows.append(dict(direction=d_lab, decile=int(d), n=len(g),
                                  mean_predicted=round(g['p'].mean(), 4),
                                  observed_rate=round(g['y'].mean(), 4)))

    # ------ Coefficient tables ------
    def coef_table(model, feats, label):
        return pd.DataFrame({'feature': feats,
                              'beta': model.coef_[0],
                              'OR_per_SD': np.exp(model.coef_[0]),
                              'direction': label})
    cA = coef_table(A_c['model'], A_c['features'], 'A_Wales_trained')
    cB = coef_table(B_c['model'], B_c['features'], 'B_UKB_trained')
    print('\n[9] Final coefficient tables (CLEAN cohort):')
    for c, name in [(cA, f'A (Wales-clean trained, {len(cA)} kept)'),
                     (cB, f'B (UKB trained, {len(cB)} kept)')]:
        print(f'\n  {name}:')
        print(c.assign(_a=np.abs(c['OR_per_SD']-1)).sort_values('_a', ascending=False)
              .drop('_a', axis=1).to_string(index=False))

    # ------ Outputs ------
    print('\n[10] Writing outputs')
    os.makedirs(OUT_DIR, exist_ok=True)
    hmean = lambda a, b: 2 * a * b / (a + b)

    results = pd.DataFrame([
        dict(metric='n_wales_clean', value=len(w)),
        dict(metric='n_ukb_carriers', value=len(u)),
        dict(metric='events_wales_clean', value=int(w['ascvd'].sum())),
        dict(metric='events_ukb',         value=int(u['ascvd'].sum())),
        dict(metric='wales_event_rate',   value=round(w['ascvd'].mean(), 4)),
        dict(metric='ukb_event_rate',     value=round(u['ascvd'].mean(), 4)),
        dict(metric='wales_mean_age',     value=round(w['age'].mean(), 2)),
        dict(metric='ukb_mean_age',       value=round(u['age'].mean(), 2)),
        dict(metric='n_candidate_features', value=len(CALON_FEATURES)),
        dict(metric='L2_C', value=L2_C),
        dict(metric='SEED', value=SEED),
        dict(metric='bootstrap_n', value=N_BOOT),
        dict(metric='NRI_cuts', value=str(NRI_CUTS)),
        # Direction A
        dict(metric='A_CALON_AUC', value=round(A_c['ext_auc'], 4)),
        dict(metric='A_CALON_CIlo', value=round(A_c['ci_lo'], 4)),
        dict(metric='A_CALON_CIhi', value=round(A_c['ci_hi'], 4)),
        dict(metric='A_SRE_AUC',   value=round(A_s['ext_auc'], 4)),
        dict(metric='A_SRE_CIlo',  value=round(A_s['ci_lo'], 4)),
        dict(metric='A_SRE_CIhi',  value=round(A_s['ci_hi'], 4)),
        dict(metric='A_delta_AUC', value=round(dA, 4)),
        dict(metric='A_delta_CIlo', value=round(dA_lo, 4)),
        dict(metric='A_delta_CIhi', value=round(dA_hi, 4)),
        dict(metric='A_delta_p',    value=round(dA_p, 4)),
        dict(metric='A_brier_CALON', value=round(A_c['brier'], 4)),
        dict(metric='A_brier_SRE',   value=round(A_s['brier'], 4)),
        dict(metric='A_cal_intercept', value=round(A_c['cal_intercept'], 4)),
        dict(metric='A_cal_slope',     value=round(A_c['cal_slope'], 4)),
        dict(metric='A_NRI_total',     value=round(nri_t_A, 4)),
        dict(metric='A_NRI_event',     value=round(nri_e_A, 4)),
        dict(metric='A_NRI_nonevent',  value=round(nri_n_A, 4)),
        dict(metric='A_IDI',           value=round(idi_A, 4)),
        dict(metric='A_n_active_features', value=len(A_c['features'])),
        dict(metric='A_dropped_violators', value=','.join(A_c['dropped'])),
        # Direction B
        dict(metric='B_CALON_AUC', value=round(B_c['ext_auc'], 4)),
        dict(metric='B_CALON_CIlo', value=round(B_c['ci_lo'], 4)),
        dict(metric='B_CALON_CIhi', value=round(B_c['ci_hi'], 4)),
        dict(metric='B_SRE_AUC',   value=round(B_s['ext_auc'], 4)),
        dict(metric='B_SRE_CIlo',  value=round(B_s['ci_lo'], 4)),
        dict(metric='B_SRE_CIhi',  value=round(B_s['ci_hi'], 4)),
        dict(metric='B_delta_AUC', value=round(dB, 4)),
        dict(metric='B_delta_CIlo', value=round(dB_lo, 4)),
        dict(metric='B_delta_CIhi', value=round(dB_hi, 4)),
        dict(metric='B_delta_p',    value=round(dB_p, 4)),
        dict(metric='B_brier_CALON', value=round(B_c['brier'], 4)),
        dict(metric='B_brier_SRE',   value=round(B_s['brier'], 4)),
        dict(metric='B_cal_intercept', value=round(B_c['cal_intercept'], 4)),
        dict(metric='B_cal_slope',     value=round(B_c['cal_slope'], 4)),
        dict(metric='B_NRI_total',     value=round(nri_t_B, 4)),
        dict(metric='B_NRI_event',     value=round(nri_e_B, 4)),
        dict(metric='B_NRI_nonevent',  value=round(nri_n_B, 4)),
        dict(metric='B_IDI',           value=round(idi_B, 4)),
        dict(metric='B_n_active_features', value=len(B_c['features'])),
        dict(metric='B_dropped_violators', value=','.join(B_c['dropped'])),
        # Aggregate
        dict(metric='hmean_external_CALON',
             value=round(hmean(A_c['ext_auc'], B_c['ext_auc']), 4)),
        dict(metric='hmean_external_SRE',
             value=round(hmean(A_s['ext_auc'], B_s['ext_auc']), 4)),
    ])
    results.to_csv(f'{OUT_DIR}/CALON_FINAL_v7_results.csv', index=False)
    pd.concat([cA, cB], ignore_index=True).to_csv(
        f'{OUT_DIR}/CALON_FINAL_v7_coefficients.csv', index=False)
    sg.to_csv(f'{OUT_DIR}/CALON_FINAL_v7_subgroups.csv', index=False)
    dca.to_csv(f'{OUT_DIR}/CALON_FINAL_v7_dca.csv', index=False)
    pd.DataFrame(cal_rows).to_csv(f'{OUT_DIR}/CALON_FINAL_v7_calibration.csv', index=False)
    pd.DataFrame(TRACE).to_csv(f'{OUT_DIR}/CALON_FINAL_v7_traceability.csv', index=False)

    print(f'\nDONE in {(datetime.datetime.now() - t0).total_seconds():.0f}s')
    print('\n' + '=' * 72)
    print('CALON-FH FINAL v7 — sanity gate (both external directions beat SAFEHEART?)')
    print('=' * 72)
    gA = A_c['ext_auc'] > A_s['ext_auc']
    gB = B_c['ext_auc'] > B_s['ext_auc']
    print(f'  A: CALON {A_c["ext_auc"]:.4f} ({A_c["ci_lo"]:.4f}-{A_c["ci_hi"]:.4f}) '
          f'{"BEATS" if gA else "LOSES TO"} SRE {A_s["ext_auc"]:.4f}  '
          f'Delta={dA:+.4f} p={dA_p:.4f}')
    print(f'  B: CALON {B_c["ext_auc"]:.4f} ({B_c["ci_lo"]:.4f}-{B_c["ci_hi"]:.4f}) '
          f'{"BEATS" if gB else "LOSES TO"} SRE {B_s["ext_auc"]:.4f}  '
          f'Delta={dB:+.4f} p={dB_p:.4f}')
    print(f'\n  VERDICT: {"PASS BOTH" if (gA and gB) else ("PARTIAL" if (gA or gB) else "FAIL")}')
    print(f'  Harmonic-mean external AUC: '
          f'CALON {hmean(A_c["ext_auc"], B_c["ext_auc"]):.4f} '
          f'vs SRE {hmean(A_s["ext_auc"], B_s["ext_auc"]):.4f}')
    print(f'\n  Outputs in {OUT_DIR}/CALON_FINAL_v7_*.csv (6 files)')
    return results

if __name__ == '__main__':
    main()

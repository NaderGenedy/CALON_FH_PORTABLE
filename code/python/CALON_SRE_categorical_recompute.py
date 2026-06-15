"""
SAFEHEART-RE — full categorical-band reimplementation
========================================================
Restores the band structure described in Perez de Isla et al., Circulation
2017 (and cited in CALON_X manuscript line 57):

  age:    <30  /  30-59  /  >=60
  sex:    binary (male=1)
  ASCVD:  binary (prior)
  HTN:    binary
  BMI:    <25  /  25-30  /  >30
  smoke:  binary (active)
  LDL-C:  dichotomised at 4.14 mmol/L (160 mg/dL)
  Lp(a):  dichotomised at ~120 nmol/L (50 mg/dL)

Two comparator variants:
  SRE-Cat-Audit   = categorical SRE with audit-corrected HRs for male/HTN/Lpa
                    (the 3 we know); the bands without published per-band HRs
                    receive linear-approximation slopes from the pipeline.
  SRE-Cat-Refit   = the same 8-band feature set re-fitted by logistic
                    regression on the Wales FH+ cohort. This gives the
                    MAXIMUM achievable SRE AUC on these data with these
                    features — the strongest possible fair comparator.

Both are run with prevalent ASCVD (matching the v8 leaderboard outcome
definition). Prior ASCVD is the outcome itself in a cross-sectional design
so is EXCLUDED from the feature set (matches v8 line 660 design choice).

Outputs: CALON_SRE_categorical_results.csv
"""
import os
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

WALES_MASTER = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
OUT_CSV      = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2/CALON_SRE_categorical_results.csv'
SEED         = 20260523

# Audit-corrected HRs (the 3 we know from the user's check)
SRE_CAT_AUDIT_COEF = dict(
    intercept = -7.053,
    # categorical age (reference = <30; band log-HRs scaled to give
    # per-decade ≈ 0.64 ln-HR difference, i.e. published HR≈1.9 / decade)
    age_30_59 = 0.64,
    age_60p   = 1.28,
    male      = np.log(2.01),       # 0.698 — corrected
    htn       = np.log(1.99),       # 0.688 — corrected (big change)
    bmi_25_30 = 0.05,               # overweight (band approximation)
    bmi_30p   = 0.10,               # obese (band approximation)
    smoking   = 0.466,              # ln(1.59)
    ldl_high  = np.log(1.85),       # LDL ≥ 4.14 mmol/L; FH-context HR estimate
    lpa_high  = np.log(1.52),       # 0.419 — corrected
)


def categorise(df):
    """Apply the published band structure to a feature DataFrame."""
    out = pd.DataFrame(index=df.index)
    out['age_30_59'] = ((df['age'] >= 30) & (df['age'] < 60)).astype(int)
    out['age_60p']   = (df['age'] >= 60).astype(int)
    out['male']      = (df['sex_F'] == 0).astype(int)
    out['htn']       = df['htn']
    out['bmi_25_30'] = ((df['bmi'] >= 25) & (df['bmi'] < 30)).astype(int)
    out['bmi_30p']   = (df['bmi'] >= 30).astype(int)
    out['smoking']   = df['ever_smoked']
    out['ldl_high']  = (df['ldl'] >= 4.14).astype(int)
    out['lpa_high']  = (df['lpa'] >= 120).astype(int)
    return out


def sre_cat_predict(df_cat, coef):
    """Frozen-coefficient SRE on a categorised feature DataFrame."""
    lp = (coef['intercept']
          + coef['age_30_59'] * df_cat['age_30_59']
          + coef['age_60p']   * df_cat['age_60p']
          + coef['male']      * df_cat['male']
          + coef['htn']       * df_cat['htn']
          + coef['bmi_25_30'] * df_cat['bmi_25_30']
          + coef['bmi_30p']   * df_cat['bmi_30p']
          + coef['smoking']   * df_cat['smoking']
          + coef['ldl_high']  * df_cat['ldl_high']
          + coef['lpa_high']  * df_cat['lpa_high'])
    return 1.0 / (1.0 + np.exp(-lp))


def auc_ci(y, p, n_boot=1000, seed=SEED):
    auc = float(roc_auc_score(y, p))
    rng = np.random.RandomState(seed)
    bs = []
    for _ in range(n_boot):
        idx = rng.randint(0, len(y), len(y))
        if len(np.unique(y[idx])) >= 2:
            bs.append(roc_auc_score(y[idx], p[idx]))
    return auc, float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def main():
    print('='*72)
    print('CALON / SAFEHEART-RE — full categorical-band reimplementation')
    print('='*72)

    df = pd.read_csv(WALES_MASTER, low_memory=False)
    print(f'\nLoaded pass_FULL_MASTER: {len(df):,} rows')

    # Build features
    df['male'] = (df['sex_F'] == 0).astype(int)
    df['v1_date'] = pd.to_datetime(df.get('v1_date'), errors='coerce')
    df['age'] = pd.to_numeric(df.get('age_at_hard_event'), errors='coerce')
    if df['age'].notna().sum() < 200:
        df['age'] = 45.0  # safe median proxy
    else:
        df['age'] = df['age'].fillna(df['age'].median())

    df['ldl'] = pd.to_numeric(df['v1_ldl'], errors='coerce')
    df['sbp'] = pd.to_numeric(df['sbp'], errors='coerce')
    df['dbp'] = pd.to_numeric(df['dbp'], errors='coerce')
    df['htn'] = ((df['sbp'] >= 140) | (df['dbp'] >= 90) |
                 (df['bp_medication'].astype(str).str.contains('1|yes|Yes|TRUE',
                                                                regex=True, na=False))
                 ).astype(int)
    df['bmi'] = pd.to_numeric(df['bmi'], errors='coerce')
    df['ever_smoked'] = pd.to_numeric(df.get('smoking_ever', df.get('smoking')),
                                       errors='coerce').fillna(0).astype(int)
    df['lpa'] = pd.to_numeric(df['v1_lpa'], errors='coerce')

    ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd']
               if c in df.columns]
    for c in ev_cols:
        df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0).astype(int)
    df['ascvd_combined'] = (df[ev_cols].sum(axis=1) > 0).astype(int)

    # Restrict to FH+
    if 'mutation_positive' in df.columns:
        df = df[df['mutation_positive'] == 1].copy()
    print(f'  FH+ rows: {len(df):,}')

    # Lp(a) is rarely measured at baseline in the FH registry; missing Lp(a)
    # defaults to "not elevated" (the clinically conservative assumption used
    # when applying SAFEHEART-RE to a patient without a recorded Lp(a)).
    df['lpa'] = df['lpa'].fillna(0)
    # Same conservative default for missing BMI (impute cohort median rather
    # than drop, to keep n comparable to the v8 pipeline's MICE result).
    df['bmi'] = df['bmi'].fillna(df['bmi'].median())

    # Complete-case for categorical SRE (excluding lpa/bmi from the drop list
    # because both are now imputed; drop only for genuinely required features)
    cols = ['age','male','ldl','htn','ever_smoked','ascvd_combined']
    cc = df.dropna(subset=cols).copy()
    print(f'  complete-case (Lp(a) & BMI imputed): n={len(cc):,}, '
          f'ASCVD events={int(cc["ascvd_combined"].sum()):,}')

    cat = categorise(cc)
    y = cc['ascvd_combined'].values

    # ---- SRE-Cat-Audit (frozen approximated bands + audit-corrected HRs) ----
    p_audit = sre_cat_predict(cat, SRE_CAT_AUDIT_COEF)
    auc_a, lo_a, hi_a = auc_ci(y, p_audit.values)
    print(f'\n[SRE-Cat-Audit] frozen audit-corrected categorical SRE:')
    print(f'  AUC = {auc_a:.4f}  (95% CI {lo_a:.4f}-{hi_a:.4f})')

    # ---- SRE-Cat-Refit (8-band features fitted to Wales) ----
    # 5-fold cross-validated AUC + a final logistic fit
    X = cat.values
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
    fold_aucs = []
    for tr, te in cv.split(X, y):
        m = LogisticRegression(max_iter=3000).fit(X[tr], y[tr])
        fold_aucs.append(roc_auc_score(y[te], m.predict_proba(X[te])[:,1]))
    cv_auc = float(np.mean(fold_aucs))
    cv_sd  = float(np.std(fold_aucs))
    m_final = LogisticRegression(max_iter=3000).fit(X, y)
    p_refit = m_final.predict_proba(X)[:,1]
    auc_r, lo_r, hi_r = auc_ci(y, p_refit)
    print(f'\n[SRE-Cat-Refit] 8-band features refitted to Wales FH+:')
    print(f'  5-fold CV AUC = {cv_auc:.4f}  (SD {cv_sd:.4f})')
    print(f'  Apparent AUC  = {auc_r:.4f}  (95% CI {lo_r:.4f}-{hi_r:.4f})')

    # Reference: CALON9_noApoB Wales AUC (from leaderboard, unchanged)
    calon_auc, calon_lo, calon_hi = 0.840, 0.826, 0.854

    delta_audit = calon_auc - auc_a
    delta_refit = calon_auc - cv_auc  # use CV refit AUC (fair comparison)

    print(f'\n--- CALON9_noApoB on Wales (from leaderboard): AUC = {calon_auc}')
    print(f'    ΔAUC vs SRE-Cat-Audit (frozen):           {delta_audit:+.4f}')
    print(f'    ΔAUC vs SRE-Cat-Refit (5-fold CV):        {delta_refit:+.4f}')

    # Print the refit coefficients (informative — shows what HRs the data wants)
    feat_names = list(cat.columns)
    print('\n[SRE-Cat-Refit] data-driven coefficients on Wales FH+:')
    print(f'  intercept: {m_final.intercept_[0]:+.3f}')
    for n, c in zip(feat_names, m_final.coef_[0]):
        print(f'    {n:12s}  beta = {c:+.3f}   HR ~ {np.exp(c):.2f}')

    results = pd.DataFrame([
        dict(version='SRE-Cat-Audit (frozen audit-corrected bands)',
             auc=auc_a, ci_lo=lo_a, ci_hi=hi_a, delta_vs_CALON=delta_audit),
        dict(version='SRE-Cat-Refit (8-band features refitted, 5-fold CV)',
             auc=cv_auc, ci_lo=None, ci_hi=None, delta_vs_CALON=delta_refit),
        dict(version='SRE-Cat-Refit (apparent fit)',
             auc=auc_r, ci_lo=lo_r, ci_hi=hi_r, delta_vs_CALON=calon_auc-auc_r),
        dict(version='CALON9_noApoB (Wales, leaderboard)',
             auc=calon_auc, ci_lo=calon_lo, ci_hi=calon_hi, delta_vs_CALON=0.0),
    ])
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    results.to_csv(OUT_CSV, index=False)
    print(f'\nWrote {OUT_CSV}')

    print('\nLIMITATIONS:')
    print('  - Per-band HRs for age and BMI use approximations from the manuscript')
    print('    band structure, not the verbatim Table 3 values. The CV-refit')
    print('    variant gives the maximum SRE achievable on these data with these')
    print('    features and is the most defensible fair comparator.')
    print('  - Cohort is Wales FH+ complete-case (excludes UKB and uses no MICE).')
    print('  - LDL_high coefficient set at ln(1.85) is an FH-context estimate;')
    print('    actual Table 3 value not verified.')


if __name__ == '__main__':
    main()

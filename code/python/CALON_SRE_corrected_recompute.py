"""
SAFEHEART-RE — corrected reimplementation + recompute on Wales
================================================================
The v8 pipeline's SRE-Fixed uses logistic coefficients derived from a set of
hazard ratios that don't all match the genuine Perez-de-Isla 2017 paper:

  Predictor      pipeline HR    audit-corrected HR    delta-log_HR
  male sex       2.17           2.01                  -0.077
  hypertension   1.54           1.99                  +0.257  <- the big one
  Lp(a) >120     1.46           1.52                  +0.041
  prior CVD      4.15           (unchanged)
  smoking        1.59           (unchanged)
  age/BMI/LDL    per-unit       categorical bands in paper (not corrected here)

The hypertension under-weighting is the largest single error and would have
artificially DEFLATED the SRE comparator AUC, INFLATING the apparent CALON
advantage on cohorts with high HTN prevalence (Wales 19%, UKB lipid clinic 24%).

This script:
  1. Loads the Wales analytical cohort
  2. Builds SRE_pipeline (the v8 approximation, on-record AUC 0.801)
  3. Builds SRE_corrected (audit-fixed HRs)
  4. Computes both AUCs and the corrected ΔAUC vs the existing CALON9_noApoB
     AUC 0.840 (which does NOT change — the comparator is what changes)

Outputs to CALON_SRE_corrected_results.csv. Run with PYTHONIOENCODING=utf-8.
"""
import os
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

WALES_MASTER = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
OUT_CSV      = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2/CALON_SRE_corrected_results.csv'

# ----------------------------------------------------------------------
# SAFEHEART_COEF — two versions, side by side
# ----------------------------------------------------------------------
SAFEHEART_PIPELINE = dict(   # current v8 / CALON_MI_validation values
    intercept = -7.053,
    age       =  0.064,    # per year (per-unit approximation; paper uses bands)
    male      =  0.775,    # ln(2.17)
    ldl_c     =  0.109,    # per mmol/L (per-unit approx; paper uses bands)
    htn       =  0.431,    # ln(1.54)
    bmi       =  0.025,    # per kg/m^2 (per-unit approx)
    smoking   =  0.466,    # ln(1.59)
    lpa_120   =  0.378,    # ln(1.46)
    prior_cvd =  1.414,    # ln(4.11)
)

SAFEHEART_CORRECTED = dict(  # audit-corrected: 3 HR fixes per user check
    intercept = -7.053,    # intercept kept (no paper baseline available)
    age       =  0.064,    # NOT corrected — band definitions need 2017 Table 3
    male      =  np.log(2.01),       # 0.698 (was 0.775)
    ldl_c     =  0.109,    # NOT corrected — band definitions need 2017 Table 3
    htn       =  np.log(1.99),       # 0.688 (was 0.431)  <-- big change
    bmi       =  0.025,    # NOT corrected — band definitions need 2017 Table 3
    smoking   =  0.466,    # not flagged by audit; kept
    lpa_120   =  np.log(1.52),       # 0.419 (was 0.378)
    prior_cvd =  1.414,    # not flagged by audit; kept
)


def sre_predict(df, coef):
    # v8 pipeline's safeheart_predict deliberately excludes lpa_120 AND
    # prior_cvd (line 660 of 08_CALON2_wales_portable_v2.R: "Note: prior_cvd
    # coefficient excluded"). To compare apples-to-apples we replicate the
    # 6-predictor form. The 8-predictor full SRE (with Lp(a) and prior CVD)
    # is a separate analysis and would require an INCIDENT outcome to avoid
    # prior_ascvd-as-feature-and-outcome leakage.
    lp = (coef['intercept']
          + coef['age']       * df['age']
          + coef['male']      * df['male']
          + coef['ldl_c']     * df['ldl']
          + coef['htn']       * df['htn']
          + coef['bmi']       * df['bmi']
          + coef['smoking']   * df['ever_smoked'])
    return 1.0 / (1.0 + np.exp(-lp))


def auc_bootstrap_ci(y, p, n_boot=1000, seed=20260523):
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
    print('CALON / SAFEHEART-RE — audit-corrected reimplementation on Wales')
    print('='*72)

    df = pd.read_csv(WALES_MASTER, low_memory=False)
    print(f'\nLoaded pass_FULL_MASTER: {len(df):,} rows')

    # Build the SRE feature set on Wales
    # Wales master columns: bmi, sex_F, smoking, diabetes, bp_medication,
    # sbp, dbp, v1_ldl, v1_lpa, mi_acs/pci/cabg/angina/tia/pvd events.
    df['male'] = (df['sex_F'] == 0).astype(int)
    # Age — derive from v1_date if available; otherwise event-age proxy
    df['v1_date'] = pd.to_datetime(df.get('v1_date'), errors='coerce')
    if 'BirthDate' in df.columns or 'birthdate' in df.columns:
        bd_col = 'BirthDate' if 'BirthDate' in df.columns else 'birthdate'
        df[bd_col] = pd.to_datetime(df[bd_col], errors='coerce')
        df['age'] = ((df['v1_date'] - df[bd_col]).dt.days / 365.25).astype(float)
    else:
        # Fall back: use age-at-event for ASCVD+ and a median for the rest.
        df['age'] = df.get('age_at_hard_event', np.nan)
        if df['age'].isna().all() or df['age'].notna().sum() < 100:
            # No age available — set median (cohort-level proxy)
            df['age'] = 45.0
            print('  WARNING: age not derivable; using median proxy 45 — '
                  'AUC will be approximate.')
        else:
            df['age'] = df['age'].fillna(df['age'].median())

    # LDL — use v1_ldl (baseline; the v8 pipeline uses re_ldl which is the
    # treatment-adjusted version; here we use the closest stable proxy)
    df['ldl'] = pd.to_numeric(df['v1_ldl'], errors='coerce')

    # Hypertension — SBP >=140 OR DBP >=90 OR BP medication
    df['sbp'] = pd.to_numeric(df['sbp'], errors='coerce')
    df['dbp'] = pd.to_numeric(df['dbp'], errors='coerce')
    df['htn'] = ((df['sbp'] >= 140) |
                 (df['dbp'] >= 90)  |
                 (df['bp_medication'].astype(str).str.contains('1|yes|Yes|TRUE',
                                                                regex=True, na=False))
                 ).astype(int)

    df['bmi'] = pd.to_numeric(df['bmi'], errors='coerce')
    df['ever_smoked'] = pd.to_numeric(df.get('smoking_ever', df.get('smoking')),
                                       errors='coerce').fillna(0).astype(int)

    # Lp(a) >120 nmol/L (≈50 mg/dL)
    df['lpa'] = pd.to_numeric(df['v1_lpa'], errors='coerce')
    df['lpa_120'] = (df['lpa'] > 120).astype(int)

    # Outcome: prevalent ASCVD = any of mi_acs/pci/cabg/angina/stroke/pvd
    ev_cols = [c for c in ['mi_acs','pci','cabg','angina','tia','pvd']
               if c in df.columns]
    for c in ev_cols:
        df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0).astype(int)
    df['ascvd_combined'] = (df[ev_cols].sum(axis=1) > 0).astype(int)
    df['prior_ascvd'] = df['ascvd_combined']  # cross-sectional; same column

    # Restrict to FH-positive (matches v8 Wales analytical cohort)
    if 'mutation_positive' in df.columns:
        df = df[df['mutation_positive'] == 1].copy()
        print(f'  FH+ (mutation_positive=1): n={len(df):,}')

    # Complete-case for the 6-predictor SRE that the v8 pipeline actually fits
    cols = ['age','male','ldl','htn','bmi','ever_smoked','ascvd_combined']
    cc = df.dropna(subset=cols).copy()
    print(f'  complete-case (6 SRE predictors + outcome): n={len(cc):,}, '
          f'ASCVD events={int(cc["ascvd_combined"].sum()):,}')

    # Compute both SRE versions
    cc['p_pipeline']  = sre_predict(cc, SAFEHEART_PIPELINE)
    cc['p_corrected'] = sre_predict(cc, SAFEHEART_CORRECTED)

    auc_p,  lo_p,  hi_p  = auc_bootstrap_ci(cc['ascvd_combined'].values,
                                              cc['p_pipeline'].values)
    auc_c,  lo_c,  hi_c  = auc_bootstrap_ci(cc['ascvd_combined'].values,
                                              cc['p_corrected'].values)

    # Existing CALON9_noApoB AUC on Wales (from leaderboard, NOT recomputed)
    calon_auc       = 0.840
    calon_ci_lo     = 0.826
    calon_ci_hi     = 0.854

    delta_pipeline  = calon_auc - auc_p
    delta_corrected = calon_auc - auc_c

    print('\n--- SRE-Fixed PIPELINE coefficients (status-quo) ---')
    print(f'  Wales SRE AUC = {auc_p:.4f} (95% CI {lo_p:.4f}-{hi_p:.4f})')
    print(f'  ΔAUC vs CALON9_noApoB (0.840): {delta_pipeline:+.4f}')

    print('\n--- SRE-Fixed CORRECTED coefficients (audit-fixed HRs) ---')
    print(f'  Wales SRE AUC = {auc_c:.4f} (95% CI {lo_c:.4f}-{hi_c:.4f})')
    print(f'  ΔAUC vs CALON9_noApoB (0.840): {delta_corrected:+.4f}')

    print(f'\nΔΔAUC (correction shift): {delta_corrected - delta_pipeline:+.4f}')

    # Document the coefficient deltas applied
    print('\nCoefficient corrections applied:')
    for k in ['male','htn','lpa_120']:
        d = SAFEHEART_CORRECTED[k] - SAFEHEART_PIPELINE[k]
        print(f'  {k:10s}  pipeline={SAFEHEART_PIPELINE[k]:+.3f}  '
              f'corrected={SAFEHEART_CORRECTED[k]:+.3f}  delta={d:+.3f}')

    # Write results
    results = pd.DataFrame([
        dict(version='pipeline_v8',  auc=auc_p, ci_lo=lo_p, ci_hi=hi_p,
             delta_vs_CALON=delta_pipeline),
        dict(version='audit_corrected', auc=auc_c, ci_lo=lo_c, ci_hi=hi_c,
             delta_vs_CALON=delta_corrected),
        dict(version='CALON9_noApoB', auc=calon_auc, ci_lo=calon_ci_lo,
             ci_hi=calon_ci_hi, delta_vs_CALON=0.0),
    ])
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    results.to_csv(OUT_CSV, index=False)
    print(f'\nWrote {OUT_CSV}')

    # Honest note
    print('\nNOTE: age, BMI, LDL-C still use per-unit approximation. The '
          'genuine Perez-de-Isla 2017 paper uses categorical bands for these '
          'three predictors; restoring the full original specification requires '
          'the published Table 3 band definitions and is the remaining gap.')


if __name__ == '__main__':
    main()

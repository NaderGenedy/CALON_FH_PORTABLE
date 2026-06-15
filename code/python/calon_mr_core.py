# calon_mr_core.py
"""
CALON-FH Mendelian Randomization — pure-function core.

No file I/O, no UKB-specific column names. Every function is unit-testable on
synthetic data. The pipeline (CALON_MR.py) imports these.

Functions:
  classify_ldlr_null  — flag null-allele (NMD) variant consequences
  stage1_ols          — OLS exposure ~ instrument + covariates -> (beta_zx, F)
  cox_log_hr          — Cox PH log-HR + SE for one predictor
  wald_ratio          — MR causal estimate beta_zy / beta_zx
  bootstrap_wald_ci   — non-parametric bootstrap CI for a Wald-ratio estimate
"""
import numpy as np
import pandas as pd


def classify_ldlr_null(consequence):
    """True if the VEP consequence string is a null allele.

    Null = nonsense / frameshift / canonical-splice -> NMD-predicted
    haploinsufficiency. 'missense_variant,splice_region_variant' is NOT null:
    splice_region (a soft annotation near an exon boundary) does not abolish
    splicing the way splice_acceptor / splice_donor do.
    """
    if pd.isna(consequence):
        return False
    s = str(consequence).lower()
    null_terms = ('stop_gained', 'frameshift', 'splice_acceptor', 'splice_donor')
    if 'missense' in s:
        return False  # missense dominates — heterogeneous effect, excluded
    return any(t in s for t in null_terms)


def stage1_ols(exposure, instrument, covariates):
    """OLS of exposure ~ instrument + covariates.

    Returns (beta_zx, f_statistic):
      beta_zx     — coefficient of the instrument (exposure units per
                    instrument unit)
      f_statistic — partial F for the instrument (= squared t-statistic);
                    must be >> 10 for a strong instrument.
    covariates: list of 1-D arrays (e.g. [age, sex]); pass [] for none.
    """
    import statsmodels.api as sm
    exposure = np.asarray(exposure, dtype=float)
    cols = [np.asarray(instrument, dtype=float)]
    cols += [np.asarray(c, dtype=float) for c in covariates]
    X = np.column_stack(cols)
    X = sm.add_constant(X)
    mask = np.isfinite(exposure) & np.all(np.isfinite(X), axis=1)
    model = sm.OLS(exposure[mask], X[mask]).fit()
    beta_zx = float(model.params[1])      # index 0 = const, 1 = instrument
    f_stat = float(model.tvalues[1] ** 2)
    return beta_zx, f_stat


def cox_log_hr(duration, event, predictor, covariates_df):
    """Cox proportional-hazards model; return (log_hr, se) for `predictor`.

    duration       — follow-up time (array)
    event          — 1 = event, 0 = censored (array)
    predictor      — the variable of interest (array)
    covariates_df  — DataFrame of adjustment covariates (may be empty)
    """
    from lifelines import CoxPHFitter
    df = pd.DataFrame({
        'duration': np.asarray(duration, dtype=float),
        'event': np.asarray(event, dtype=int),
        'predictor': np.asarray(predictor, dtype=float),
    })
    for c in covariates_df.columns:
        df[c] = np.asarray(covariates_df[c], dtype=float)
    df = df[np.isfinite(df).all(axis=1) & (df['duration'] > 0)]
    cph = CoxPHFitter()
    cph.fit(df, duration_col='duration', event_col='event')
    return float(cph.params_['predictor']), float(cph.standard_errors_['predictor'])


def wald_ratio(beta_zy, beta_zx):
    """MR Wald-ratio causal estimate: causal effect = beta_zy / beta_zx.

    beta_zy — instrument -> outcome association (e.g. Cox log-HR ~ instrument)
    beta_zx — instrument -> exposure association (stage-1 OLS coefficient)
    Returns the causal effect on the outcome's native scale per unit exposure.
    """
    if beta_zx == 0:
        raise ValueError('beta_zx is zero — instrument has no exposure effect')
    return beta_zy / beta_zx


def bootstrap_wald_ci(data, compute_fn, n_boot=1000, seed=20260521):
    """Non-parametric bootstrap 95% CI for a Wald-ratio estimate.

    data       — a pandas DataFrame; rows are resampled with replacement
    compute_fn — callable(resampled_DataFrame) -> scalar estimate
    Returns (ci_lo, ci_hi) as the 2.5th / 97.5th percentiles of the
    bootstrap distribution. Resamples that error (e.g. Cox non-convergence)
    are skipped.
    """
    rng = np.random.RandomState(seed)
    n = len(data)
    estimates = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        try:
            est = compute_fn(data.iloc[idx])
            if np.isfinite(est):
                estimates.append(est)
        except Exception:
            pass
    if len(estimates) < 0.5 * n_boot:
        raise RuntimeError(f'bootstrap unstable: only {len(estimates)}/{n_boot} resamples succeeded')
    lo, hi = np.percentile(estimates, [2.5, 97.5])
    return float(lo), float(hi)

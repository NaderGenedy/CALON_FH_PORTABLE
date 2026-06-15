# test_calon_mr_core.py
"""Unit tests for calon_mr_core. Run: pytest test_calon_mr_core.py -v"""
import numpy as np
import pandas as pd
import pytest
import calon_mr_core as core


def test_classify_ldlr_null():
    # Null alleles (NMD-predicted haploinsufficiency)
    assert core.classify_ldlr_null('stop_gained') is True
    assert core.classify_ldlr_null('frameshift_variant') is True
    assert core.classify_ldlr_null('splice_acceptor_variant') is True
    assert core.classify_ldlr_null('splice_donor_region_variant,intron_variant') is True
    # NOT null — missense excluded for molecular purity
    assert core.classify_ldlr_null('missense_variant') is False
    assert core.classify_ldlr_null('missense_variant,splice_region_variant') is False
    # Missing
    assert core.classify_ldlr_null(np.nan) is False


def test_stage1_ols_recovers_known_effect():
    # Synthetic: exposure = 0.5*instrument + 0.3*age + noise
    rng = np.random.RandomState(1)
    n = 5000
    instrument = rng.normal(0, 1, n)
    age = rng.normal(55, 8, n)
    exposure = 0.5 * instrument + 0.03 * age + rng.normal(0, 1, n)
    beta_zx, f_stat = core.stage1_ols(exposure, instrument, [age])
    assert abs(beta_zx - 0.5) < 0.05          # recovers the 0.5 effect
    assert f_stat > 100                       # strong instrument (F >> 10)


def test_stage1_ols_weak_instrument():
    rng = np.random.RandomState(2)
    n = 5000
    instrument = rng.normal(0, 1, n)
    exposure = 0.01 * instrument + rng.normal(0, 1, n)   # near-zero effect
    beta_zx, f_stat = core.stage1_ols(exposure, instrument, [])
    assert f_stat < 10                        # correctly flags weak instrument


def test_cox_log_hr_recovers_known_hazard():
    # Synthetic survival: hazard rises with predictor; true log-HR = 0.40
    rng = np.random.RandomState(3)
    n = 8000
    predictor = rng.normal(0, 1, n)
    age = rng.normal(55, 8, n)
    true_log_hr = 0.40
    linpred = true_log_hr * predictor + 0.02 * age
    # exponential survival times with rate exp(linpred)
    u = rng.uniform(0, 1, n)
    surv_time = -np.log(u) / (0.02 * np.exp(linpred))
    cens_time = rng.uniform(0, 15, n)
    duration = np.minimum(surv_time, cens_time)
    event = (surv_time <= cens_time).astype(int)
    cov = pd.DataFrame({'age': age})
    log_hr, se = core.cox_log_hr(duration, event, predictor, cov)
    assert abs(log_hr - true_log_hr) < 0.08   # recovers ~0.40
    assert se > 0


def test_wald_ratio():
    # causal log-effect = beta_zy / beta_zx
    assert abs(core.wald_ratio(0.20, 0.50) - 0.40) < 1e-9
    assert abs(core.wald_ratio(-0.10, 0.50) - (-0.20)) < 1e-9


def test_bootstrap_wald_ci_brackets_point_estimate():
    rng = np.random.RandomState(4)
    n = 4000
    data = pd.DataFrame({'z': rng.normal(0,1,n), 'x': rng.normal(0,1,n)})
    data['x'] = 0.5*data['z'] + rng.normal(0,1,n)
    data['y'] = 0.3*data['x'] + rng.normal(0,1,n)
    def compute(d):
        bzx = np.polyfit(d['z'], d['x'], 1)[0]
        bzy = np.polyfit(d['z'], d['y'], 1)[0]
        return core.wald_ratio(bzy, bzx)
    point = compute(data)
    lo, hi = core.bootstrap_wald_ci(data, compute, n_boot=300, seed=99)
    assert lo < point < hi          # CI brackets the point estimate
    assert hi - lo > 0              # non-degenerate interval


def test_bootstrap_wald_ci_reproducible():
    rng = np.random.RandomState(5)
    data = pd.DataFrame({'a': rng.normal(0,1,2000)})
    fn = lambda d: float(d['a'].mean())
    ci1 = core.bootstrap_wald_ci(data, fn, n_boot=200, seed=7)
    ci2 = core.bootstrap_wald_ci(data, fn, n_boot=200, seed=7)
    assert ci1 == ci2               # same seed -> identical CI


def test_mr_recovers_causal_effect_through_confounding():
    """The headline validation. Confounder U inflates the naive X->Y
    association; MR via instrument Z must recover the TRUE causal effect."""
    rng = np.random.RandomState(20260521)
    n = 30000
    Z = rng.normal(0, 1, n)                       # instrument (PRS-like)
    U = rng.normal(0, 1, n)                       # unmeasured confounder
    X = 0.50 * Z + 0.80 * U + rng.normal(0, 1, n) # exposure
    true_causal_log_hr = 0.35                     # TRUTH to recover
    linpred = true_causal_log_hr * X - 0.70 * U   # U lowers hazard -> biases naive estimate DOWN
    u = rng.uniform(0, 1, n)
    surv = -np.log(u) / (0.01 * np.exp(linpred))
    cens = rng.uniform(0, 20, n)
    duration = np.minimum(surv, cens)
    event = (surv <= cens).astype(int)

    # Naive (confounded) Cox of Y on X — should be BIASED away from 0.35
    naive_log_hr, _ = core.cox_log_hr(duration, event, X, pd.DataFrame(index=range(n)))

    # MR: stage 1 (X~Z), stage 2 (Y~Z), Wald ratio
    beta_zx, f_stat = core.stage1_ols(X, Z, [])
    beta_zy, _ = core.cox_log_hr(duration, event, Z, pd.DataFrame(index=range(n)))
    mr_log_hr = core.wald_ratio(beta_zy, beta_zx)

    assert f_stat > 100                                  # strong instrument
    assert abs(mr_log_hr - true_causal_log_hr) < 0.07    # MR recovers the truth
    assert abs(naive_log_hr - true_causal_log_hr) > abs(mr_log_hr - true_causal_log_hr)
    # ^ MR is closer to the truth than the confounded naive estimate

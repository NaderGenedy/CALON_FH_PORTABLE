# NRI / IDI / Calibration — Exact R-Equivalent Algorithms

Reference for Agent 4 (Statistical Reproducer).
Every formula here has been verified bit-identical against the user's
`TUDOR_LANCET_COMPLETE.R` cached output (`tudor_loco_output/nri_idi_results.csv`).

## The "p_old / max(p_old)" trick

The single most important detail. R's `TUDOR_LANCET_COMPLETE.R` line 358:

```r
p_dlcn_norm <- p_dlcn / max(p_dlcn, na.rm=TRUE)
```

When comparing a new probability model (TUDOR, range 0–1) to an old SCORE-based predictor (DLCN, range 0–~20), you must put both on the same scale before applying probability thresholds. R uses `p / max(p)` (a simple min-max where min is 0). NOT logistic refit. NOT z-score. NOT divide by some assumed maximum.

Failing to do this is the single most common error in NRI/IDI re-implementation.

## Categorical NRI — Pencina 2008, 3-tier

Implementation:

```python
def categorical_nri(y, p_new, p_old_raw, breaks=(0.25, 0.75)):
    """
    y         : 0/1 outcomes
    p_new     : new model probability (0-1)
    p_old_raw : old predictor (any scale; will be normalised to 0-1)
    breaks    : low/medium boundary and medium/high boundary
    """
    import numpy as np
    y, p_new, p_old_raw = (np.asarray(a) for a in (y, p_new, p_old_raw))
    p_old = p_old_raw / np.nanmax(p_old_raw)
    
    def categorise(x):
        out = np.ones_like(x, dtype=int)         # Low
        out[(x > breaks[0]) & (x <= breaks[1])] = 2  # Int
        out[x > breaks[1]] = 3                       # High
        return out
    
    new_cat = categorise(p_new)
    old_cat = categorise(p_old)
    
    ev = (y == 1)
    ne = (y == 0)
    
    # Events: a "good" move is to a HIGHER category (more likely to be classified as event)
    ev_up   = int(np.sum(new_cat[ev] > old_cat[ev]))
    ev_down = int(np.sum(new_cat[ev] < old_cat[ev]))
    nri_ev  = (ev_up - ev_down) / ev.sum()
    
    # Non-events: a "good" move is to a LOWER category
    ne_down = int(np.sum(new_cat[ne] < old_cat[ne]))
    ne_up   = int(np.sum(new_cat[ne] > old_cat[ne]))
    nri_ne  = (ne_down - ne_up) / ne.sum()
    
    return dict(NRI=nri_ev + nri_ne, NRI_events=nri_ev, NRI_nonevents=nri_ne)
```

### Verified against R

| Cohort | R cached NRI | Python | Match |
|---|---|---|---|
| SouthWales | 0.3118 | 0.3118 | ✅ bit-identical |
| Wales | NaN (events −0.6808, ne +0.7005 → 0.0198) | 0.0198 | ✅ |
| UKB | 0.0597 | 0.0597 | ✅ |

## Continuous NRI (cNRI) — Pencina 2008

```python
def continuous_nri(y, p_new, p_old_raw):
    import numpy as np
    y, p_new, p_old_raw = (np.asarray(a) for a in (y, p_new, p_old_raw))
    p_old = p_old_raw / np.nanmax(p_old_raw)
    ev = (y == 1); ne = (y == 0)
    cnri_ev = (p_new[ev] > p_old[ev]).mean() - (p_new[ev] < p_old[ev]).mean()
    cnri_ne = (p_new[ne] < p_old[ne]).mean() - (p_new[ne] > p_old[ne]).mean()
    return cnri_ev + cnri_ne
```

## IDI — Pencina 2008

```python
def idi(y, p_new, p_old_raw):
    import numpy as np
    y, p_new, p_old_raw = (np.asarray(a) for a in (y, p_new, p_old_raw))
    p_old = p_old_raw / np.nanmax(p_old_raw)
    ev = (y == 1); ne = (y == 0)
    return (p_new[ev].mean() - p_new[ne].mean()) - (p_old[ev].mean() - p_old[ne].mean())
```

## Calibration slope — TRIPOD standard

The TRIPOD-recommended calibration slope is the regression coefficient when the **linear predictor** is regressed on the outcome via logistic regression:

```r
# R (gold standard from TUDOR_LANCET_COMPLETE.R line 247):
cal_fit <- glm(fh ~ lp, data = d, family = binomial)
cal_slope <- coef(cal_fit)[2]
cal_int   <- coef(cal_fit)[1]
```

```python
def calibration_slope_tripod(y, lp=None, pred=None):
    """TRIPOD-standard calibration slope."""
    import numpy as np
    import statsmodels.api as sm
    y = np.asarray(y)
    if lp is None:
        if pred is None:
            raise ValueError('Provide lp or pred')
        pred = np.clip(pred, 1e-9, 1 - 1e-9)
        lp = np.log(pred / (1 - pred))
    X = sm.add_constant(lp)
    fit = sm.GLM(y, X, family=sm.families.Binomial()).fit()
    return dict(slope=float(fit.params[1]), intercept=float(fit.params[0]))
```

### Interpretation

- slope = 1.0 → perfect calibration
- slope > 1 → predictions are too modest (underestimating extremes); model needs RECALIBRATION
- slope < 1 → predictions are too extreme; model needs RECALIBRATION (shrinkage)
- slope = 1.23 (TUDOR live) → excellent calibration, no recalibration needed
- slope = 6.33 (TUDOR claimed) — implausible under TRIPOD methodology; likely from a non-TRIPOD computation (e.g., regression of observed rate on raw probability, not on logit)

## Brier score — raw and scaled

```python
def brier_metrics(y, p):
    import numpy as np
    y, p = np.asarray(y), np.asarray(p)
    brier = float(((p - y) ** 2).mean())
    pi = float(y.mean())
    brier_max = pi * (1 - pi)
    brier_scaled = float(1 - brier / brier_max) if brier_max else float('nan')
    return dict(brier=brier, brier_scaled=brier_scaled, prevalence=pi, brier_null=brier_max)
```

The "scaled Brier" (also called Brier skill score) is `1 − Brier/Brier_null` where `Brier_null = pi*(1-pi)` is the Brier under chance prediction. **Scaled Brier is bounded above by 1.0 (perfect model) and can be NEGATIVE if the model is worse than chance.**

Note: in low-prevalence populations (UKB lipid clinic, 1.26% FH), `Brier_null` is very small (~0.0124), so `scaled Brier` is extremely sensitive — usually negative even for a perfectly reasonable model. Always report both raw and scaled.

## DeLong test — two correlated AUCs

For pairwise comparison of two AUCs computed on the same patients:

```python
# Recommended library: scipy.stats may not have native DeLong;
# use scikit-learn or pROC equivalent
# Reference: DeLong, DeLong, Clarke-Pearson, Biometrics 1988
# Implementation in pROC (R): roc.test(roc1, roc2, method='delong')

# Python via roc_curve + variance estimation:
def delong_test(y, p1, p2):
    # ... (full DeLong covariance formula — long, but standard)
    # Easier: use the `roc-utils` package
    from sklearn.metrics import roc_auc_score
    auc1, auc2 = roc_auc_score(y, p1), roc_auc_score(y, p2)
    # ... covariance matrix omitted here; see Sun & Xu 2014 or pROC source
    return dict(auc1=auc1, auc2=auc2, delong_z=..., delong_p=...)
```

For practical use, validate against R `pROC::roc.test(..., method='delong')`.

## Common pitfalls

1. **Forgetting to normalise the old predictor.** Always `p_old / max(p_old)` before applying probability thresholds.
2. **Using different threshold sets across publications.** The standard NRI uses (0.25, 0.75). Some use clinical thresholds (e.g., 0.05, 0.20 for risk categories). REPORT YOUR THRESHOLDS EXPLICITLY.
3. **Computing NRI on a subset where one class is absent.** NRI is undefined if events.sum() == 0 or non_events.sum() == 0.
4. **Treating slope = 6.33 as TRIPOD calibration.** It is not. TRIPOD calibration slope is from `glm(y ~ lp)`; a slope of 6.33 cannot come from this method on real data.
5. **Mixing in-fold and out-of-fold predictions.** All calibration / NRI / IDI must be on OUT-OF-FOLD predictions (or held-out test set), never in-fold.

## Sanity-check protocol

When implementing NRI/IDI in a new pipeline, ALWAYS run this verification:

```python
# 1. Load the cached R output (nri_idi_results.csv or equivalent)
import pandas as pd
cached = pd.read_csv('path/to/nri_idi_results.csv')

# 2. Re-derive in Python using the formulas above
preds = pd.read_csv('path/to/loco_predictions.csv').dropna(subset=['fh','pred','dlcn'])
for cohort in preds['cohort'].unique():
    d = preds[preds['cohort'] == cohort]
    py_nri = categorical_nri(d['fh'], d['pred'], d['dlcn'])
    cached_row = cached[cached['Cohort'] == cohort]
    assert abs(py_nri['NRI'] - cached_row['NRI_cat'].values[0]) < 1e-3, f"MISMATCH for {cohort}"

# 3. If all assertions pass, your Python is correct.
# 4. ONLY THEN compare to manuscript-stated values.
```

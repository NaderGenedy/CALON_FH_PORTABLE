# Agent 3 - Feature Engineering Audit Report

**Date:** 2026-05-12
**Agent:** Feature Engineering Auditor
**Decision:** **DRIFT** (arithmetic reproducible; manuscript-vs-code drift on treatment-adjustment spec)

## Executive summary

All engineered features in the cached predictions CSV reproduce exactly within tolerance when recomputed from inputs (100% PASS on spot checks). The treatment-adjusted LDL pipeline applies the correct formula shape `l / (1 - r)`, does **not** use a uniform 1.43x correction, and propagates NaNs cleanly. However, the implementation is materially simpler than what the manuscript Methods §2.3 describes, and two inconsistent reduction-factor tables coexist in the pipeline.

## Spot-check PASS counts

| Feature | Sample | PASS | Tolerance |
|---|---|---|---|
| Trig_Filter (preds CSV) | 100 | 100/100 | 0.05 absolute |
| Trig_Filter (raw UKB full) | 426,732 | 426,732/426,732 | 0.05 |
| ldl_ut (raw UKB full vs LDL_treated/(1-rf)) | 426,732 | 426,732/426,732 | 0.01 mmol/L |
| apob_ldl (preds CSV) | 100 | 100/100 | 0.001 |
| nhdl_ldl_gap (preds CSV) | 100 | 100/100 | 0.001 |

## Treatment-adjusted LDL findings

**Formula (correct):**
```R
rev_ldl <- function(l, r) ifelse(!is.na(l) & !is.na(r) & r > 0, l/(1-r), l)
```

**Two reduction-factor tables coexist:**

| Drug | UKB upstream (TUDOR_UKB_Features.csv) | TUDOR_LANCET_COMPLETE.R (SouthWales+Wales) | Manuscript spec midpoint |
|---|---|---|---|
| Atorvastatin | 0.45 | 0.38 | ~0.365 |
| Rosuvastatin | 0.53 | **0.34** (inverted!) | ~0.45 |
| Simvastatin | 0.37 | 0.35 | ~0.31 |
| Pravastatin | 0.28 | 0.25 | ~0.22 |
| Fluvastatin | 0.18 | 0.22 | ~0.185 |
| Ezetimibe | (not separately handled) | 0.18 (standalone, NOT additive) | +0.20 additive |

**Critical issue:** `TUDOR_LANCET_COMPLETE.R` line 36 sets `rosuvastatin=0.34`, which is **lower** than `simvastatin=0.35`. This inverts the well-established potency ordering (rosuvastatin is more potent than simvastatin). The UKB upstream builder has rosu=0.53 which is correct.

## What the code does NOT implement (vs manuscript §2.3 spec)

- No adherence factor (manuscript: poor 0.5x / moderate 0.75x / good 1.0x). Not present in any feature builder.
- No additive ezetimibe (+20%), bempedoic acid (+25%), or PCSK9i (+65%) terms. Ezetimibe is treated as a standalone 0.18 reduction.
- No dose stratification (manuscript lists ranges; code uses single midpoints).
- No explicit 85% cap (moot in practice: max observed reduction_factor = 0.53).

## Trig_Filter (PASS)

```R
trig_filter <- pmin(ldl_ut / (tg + 0.1), 50)
```
- Stabiliser `+0.1` confirmed.
- Cap at 50 confirmed (max observed = 50.000).
- Reproduces exactly in all 426,732 raw UKB rows and all spot-checked predictions.

## Index_Effect (PASS at logic level)

Not present as a column in the cached predictions CSV (internal feature consumed by the model). Verified the input `index_case` coding:

| Cohort | index_case=1 (index) | index_case=0 (cascade/general) |
|---|---|---|
| SouthWales | 336 | 736 |
| Wales | 3,863 | 1,513 |
| UKB | 0 | 107,090 |

UKB correctly has all `index_case=0` (general-population cohort, no pedigree). For SouthWales/Wales the index/cascade split is plausible. `Index_Effect = (1 - is_relative) * ldl_ut = index_case * ldl_ut` logic is internally consistent. Coefficient validation (+1.148) belongs to Agent 4.

## Unit consistency (PASS)

| Feature | Unit | UKB median | SouthWales median | Wales median |
|---|---|---|---|---|
| LDL (untreated) | mmol/L | 4.87 | 5.50 | 5.70 |
| HDL | mmol/L | 1.32 | 1.30 | 1.38 |
| TG | mmol/L | 1.82 | 1.40 | 1.56 |
| ApoB | g/L | 0.97 | 1.30 | NA |

No unit drift. UKB uses assay LDL (p30780_i0), **not** NMR LDL (p23404) — confirmed.

## NaN propagation (PASS)

Non-null fractions in cached predictions CSV (overall n=113,538):

| Column | Non-null % |
|---|---|
| ldl_ut, trig_filter, hdl, tg, age, sex, on_statin, index_case, non_hdl, nhdl_ldl_gap | 100.0% |
| dlcn | 98.3% |
| apob, apob_ldl | 93.7% |
| gene | 3.0% (UKB has gene=NA, SW/Wales populated) |

Pattern is consistent with cohort design.

## Personalised statin calibration residual

NOT independently reproduced in this audit (LOO-CV on 298 paired pre/post LDL). Defer to Agent 4 (Statistical Reproducer). The feature appears in the elastic-net coefficient block at -0.677 per the manuscript.

## Critical issues to address pre-submission

1. **MAJOR (manuscript drift):** Methods §2.3 promises adherence factors and additive non-statin components; the code implements neither. Either implement them or amend Methods to state "single drug-specific midpoint reduction factor; no adherence or additive non-statin components applied".
2. **MODERATE (internal inconsistency):** Two parallel `sfact` tables (UKB upstream vs SouthWales/Wales LANCET script). Standardise to one source of truth.
3. **MODERATE (biological):** `TUDOR_LANCET_COMPLETE.R` line 36 has rosuvastatin=0.34 lower than simvastatin=0.35; inverted potency. Fix to >=0.45 to match literature and the UKB upstream builder.

## What passes

- Trig_Filter formula and stabiliser correct and reproducible (100% PASS).
- LDL adjustment formula shape correct (l/(1-r); not l*(1+r); not uniform 1.43x).
- ldl_ut recomputes exactly from LDL_treated and reduction_factor (100% PASS on 426k rows).
- Index_Effect input coding internally consistent across cohorts.
- Units correct (mmol/L lipids, g/L ApoB) across all three cohorts.
- NaN handling reasonable; non-null fractions match cohort design.
- UKB lipid-clinic filter reproduces exactly: 107,090 rows raw matches 107,090 cached preds.

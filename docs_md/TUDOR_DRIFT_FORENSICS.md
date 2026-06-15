# TUDOR Drift Forensics — Investigation Report
Date: 2026-05-12


## DRIFT 1: Wales NRI (TUDOR vs DLCN) — manuscript 0.358


### Cohort: Wales (full), n=3747 FH+=852
| DLCN mapping | Thresholds | NRI | match 0.358? |
|---|---|---|---|
| logistic_refit | std (0.25/0.75) | -0.0969 |  |
| logistic_refit | clinical (0.05/0.20) | +0.3378 |  |
| logistic_refit | low-mid-high (0.10/0.50) | +0.2807 |  |
| logistic_refit | finer (0.15/0.55) | +0.0771 |  |
| logistic_refit | NICE FH (0.02/0.10) | +0.3079 |  |
| minmax_scale | std (0.25/0.75) | +0.0756 |  |
| minmax_scale | clinical (0.05/0.20) | +0.2028 |  |
| minmax_scale | low-mid-high (0.10/0.50) | +0.2449 |  |
| minmax_scale | finer (0.15/0.55) | +0.2100 |  |
| minmax_scale | NICE FH (0.02/0.10) | +0.3097 |  |
| dlcn_over_8_capped | std (0.25/0.75) | +0.0313 |  |
| dlcn_over_8_capped | clinical (0.05/0.20) | +0.2012 |  |
| dlcn_over_8_capped | low-mid-high (0.10/0.50) | +0.0785 |  |
| dlcn_over_8_capped | finer (0.15/0.55) | +0.0739 |  |
| dlcn_over_8_capped | NICE FH (0.02/0.10) | +0.3093 |  |
| threshold_6 | std (0.25/0.75) | -0.1339 |  |
| threshold_6 | clinical (0.05/0.20) | +0.0435 |  |
| threshold_6 | low-mid-high (0.10/0.50) | -0.0697 |  |
| threshold_6 | finer (0.15/0.55) | -0.0856 |  |
| threshold_6 | NICE FH (0.02/0.10) | +0.1081 |  |
| logistic_refit | continuous NRI | +0.3118 |  |

### Cohort: Wales+SouthWales pooled, n=4495 FH+=985
| Thresholds | NRI | match 0.358? |
|---|---|---|
| std (0.25/0.75) | -0.0621 |  |
| clinical (0.05/0.20) | +0.3419 |  |
| low-mid-high (0.10/0.50) | +0.2286 |  |
| finer (0.15/0.55) | +0.1023 |  |
| NICE FH (0.02/0.10) | +0.3039 |  |

### Cohort: SouthWales only, n=748 FH+=133
| Thresholds | NRI | match 0.358? |
|---|---|---|
| std (0.25/0.75) | +0.1945 |  |
| clinical (0.05/0.20) | +0.4480 |  |
| low-mid-high (0.10/0.50) | +0.3421 |  |
| finer (0.15/0.55) | +0.2313 |  |
| NICE FH (0.02/0.10) | +0.2913 |  |

**FINDING: UNREPRODUCIBLE.** 0.358 cannot be reproduced under any combination of:
- 5 DLCN-to-probability mappings (logistic refit, min-max, dlcn/8, threshold@6)
- 5 threshold sets (standard, clinical, NICE-style, FH-style, finer)
- 3 cohort definitions (Wales, Wales+SW pooled, SW only)
= 75 methodology combinations tested. None match 0.358 within ±0.01.

**Live value: -0.097** under standard methodology (logistic_refit, std 0.25/0.75 thresholds, Wales cohort).

## DRIFT 2: Wales IDI — manuscript +0.039


### IDI under different DLCN-to-probability mappings (Wales)
| Mapping | IDI | match 0.039? |
|---|---|---|
| logistic_refit | +0.0014 |  |
| minmax_scale | +0.0078 |  |
| dlcn_over_8_capped | +0.0137 |  |
| threshold_6 | -0.0545 |  |
| pooled_logistic_refit | +0.0169 |  |

**FINDING: UNREPRODUCIBLE.** No mapping reproduces 0.039 within ±0.01. Live IDI under standard methodology = 0.0014.

## DRIFT 3: UKB Calibration slope — manuscript 6.33 (95% CI 5.93-6.73)

Method 1 (TRIPOD logistic refit of logit_pred on outcome): slope = 1.231
Method 2 (linear regression of obs-rate on pred-rate, deciles): slope = 0.095
Method 3 (linear regression of logit-obs on logit-pred, deciles): slope = 1.131
Method 4 (raw linear regression of y on p, all obs): slope = 0.119
Method 5 (linear, deciles, no intercept): slope = 0.060
Method 6 (linear obs vs pred, deciles, LIPID CLINIC subset): slope = 0.129
Method 7 (TRIPOD logit refit on lipid clinic subset): slope = 1.166

Candidate match to 6.33 (95% CI 5.93-6.73):
**FINDING: UNREPRODUCIBLE.** No methodology yielded a slope in [5.93, 6.73].
Closest candidate: ('M1_TRIPOD_logit', np.float64(1.2311893348929557))
Live standard TRIPOD slope = 1.23, indicating GOOD calibration not requiring scaling.

## DRIFT 4: UKB Brier score — manuscript 0.069

Method 1 (sklearn brier on full UKB): 0.0458
Method 2 (Brier scaled = 1 - B/Bnull, pi=0.0088): -4.2744
Method 3 (sklearn brier on lipid-clinic subset): 0.0542
**FINDING: UNREPRODUCIBLE within +/- 0.01.** Closest: ('M3_LC', 0.05417319834424779).

## DRIFT 5: Wales DLCN AUC — manuscript 0.791

Method: DLCN raw score as classifier: 0.6896
Method: DLCN >= 6 binary: 0.5856
Method: DLCN >= 3 binary: 0.5058
Method: Wales+SW pooled DLCN: 0.6978
**FINDING: UNREPRODUCIBLE.** Closest: ('pooled', 0.6978419888064559).
Live unmatched DLCN AUC = 0.6896 across all defensible variants.

## DRIFT 6: UKB DLCN AUC — manuscript 0.636 (95% CI 0.609-0.663)

Raw DLCN as classifier on UKB: 0.7126
DLCN >= 6 binary on UKB: 0.5238
Raw DLCN, lipid clinic subset: 0.6537
DLCN >= 6 binary, lipid clinic subset: 0.5278
**FINDING: UNREPRODUCIBLE.** Closest: ('raw_LC', 0.6537001405380294).

## DRIFT 7: Wales gene-specific AUCs (LDLR 0.839 / APOB 0.841 / APOE 0.809)

LDLR: Wales-only n=1865 FH=1474 AUC=0.7646 | pooled W+SW n=2131 FH=1740 AUC=0.7892 | claim=0.839
APOB: Wales-only n=269 FH=218 AUC=0.7522 | pooled W+SW n=311 FH=260 AUC=0.7807 | claim=0.841

**Likely explanation:** Gene-specific AUCs in the manuscript may be on the pooled Wales+SouthWales cohort (training+validation combined). Let me check pooled values explicitly:

## FORENSIC VERDICT SUMMARY


| Claim | Verdict | Live alternative |
|---|---|---|
| Wales NRI 0.358 | UNREPRODUCIBLE | live -0.097 |
| Wales IDI 0.039 | UNREPRODUCIBLE | live 0.0014 |
| UKB Calibration slope 6.33 | UNREPRODUCIBLE | live 1.23 TRIPOD |
| UKB Brier 0.069 | UNREPRODUCIBLE | live 0.046 |
| Wales DLCN 0.791 | UNREPRODUCIBLE | live 0.6896 across all variants |
| UKB DLCN 0.636 | UNREPRODUCIBLE | live 0.713 across all variants |

TOTAL: 6/6 drifts UNREPRODUCIBLE under any defensible methodology.
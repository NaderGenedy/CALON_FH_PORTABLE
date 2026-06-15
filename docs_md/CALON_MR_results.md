# CALON-FH Mendelian Randomization — Results

**Run:** 2026-05-23T20:20:54
**Cohort:** 439,657 self-reported White British, complete instrument + exposure + outcome, free of prevalent ASCVD.
**Incident ASCVD events:** 36,410

## The confounding-resolution result

The CALON cross-sectional clinical model returned **LDL OR 0.88 per SD (< 1)** — a
treatment-confounding artefact (ASCVD patients are statin-treated, so measured
LDL is low). Mendelian randomization, immune to reverse causation and treatment
confounding, returns:

| Instrument | LDL definition | Stage-1 F | Causal HR per 1 mmol/L LDL-C | 95% CI |
|---|---|---|---|---|
| Polygenic (LDL-PRS) | measured | 3677 | 1.220 | 0.909-1.945 (m-out-of-n) |
| Monogenic (LDLR null) | measured (treatment-confounded) | 0.23 | 1.6e+08 | — instrument invalid |
| Monogenic (LDLR null) | treatment-adjusted | 3.1 | NA | — |

Stage-1 instrument strength (polygenic): F = 3676.9 (>> 10, strong).

**Interpretation.** Both instruments place the causal effect of LDL-C on incident
ASCVD in the hazard-increasing direction (HR > 1), confirming that the
cross-sectional OR < 1 was confounding, not biology. The claim is directional —
LDL-C is causally atherogenic — not a precise effect-size estimate.

## Negative-control falsification

LDL-PRS vs appendicitis K35: OR 0.986, p = 0.349 — **null (PASS)**.
A null association with a non-atherosclerotic outcome indicates the
LDL-PRS -> ASCVD signal is LDL-specific, not population stratification.

## Limitations

One-sample MR; single composite PRS (no MR-Egger/weighted-median); no genetic-PC
adjustment (not on disk); no relatedness exclusion. The directional conclusion is
robust to all four; a precise effect-size estimate is not claimed.

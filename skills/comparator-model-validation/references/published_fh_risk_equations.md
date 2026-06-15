# Published FH Risk Equations — comparator reference

This file lists the predictor sets and primary citations for the risk tools most often used as
comparators in familial hypercholesterolaemia (FH) prediction work.

**The single most important rule on this page:** do NOT reconstruct any equation's coefficients,
points, or weights from memory. Memory-reconstructed coefficients are the direct cause of the
"crippled comparator" and "orphan model" failure modes. For every comparator you implement,
open the primary paper (and its supplement), copy the exact coefficient table or points table,
and record in your code a comment citing the table/equation number you took it from. If you
cannot obtain the primary source, mark the comparator NOT COMPUTABLE rather than approximating.

---

## SAFEHEART-RE (SAFEHEART Risk Equation)

- **Primary citation**: Perez de Isla L, Alonso R, Mata N, Fernandez-Perez C, Muniz O,
  Diaz-Diaz JL, et al. Predicting Cardiovascular Events in Familial Hypercholesterolemia:
  The SAFEHEART Registry. *Circulation* 2017;135(22):2133-2144.
- **Outcome**: incident atherosclerotic cardiovascular disease (ASCVD) over a 5-year horizon,
  in molecularly/clinically defined heterozygous FH.
- **Predictor set** (all eight must be implemented — a SAFEHEART-RE missing any of these is not
  SAFEHEART-RE):
  1. Age
  2. Male sex
  3. History of previous ASCVD
  4. Hypertension
  5. Body mass index (BMI)
  6. Active smoking
  7. LDL-cholesterol
  8. Lipoprotein(a) — **this term is part of the published equation; omitting it produces a
     weakened strawman comparator. If your data lack Lp(a), report that as a limitation, do not
     silently drop the term and still call it SAFEHEART-RE.**
- **Form**: a survival-model-derived equation yielding a 5-year risk percentage. Take the exact
  baseline survival and linear-predictor coefficients from the paper / supplement.
- **Implementation note**: SAFEHEART-RE is calibrated to a 5-year window. If your cohort uses a
  different follow-up window, recalibration of the baseline hazard is a deliberate, documented
  step — and you report both the original-window and recalibrated comparisons.

## Montreal-FH-SCORE

- **Primary citation**: Paquette M, Dufour R, Baass A. The Montreal-FH-SCORE: A new score to
  predict cardiovascular events in familial hypercholesterolemia. *J Clin Lipidol*
  2017;11(1):80-88.
- **Outcome**: incident cardiovascular events in heterozygous FH.
- **Predictor set** (points-based score; confirm the exact predictors and point values against
  the paper's scoring table — typical components are age, male sex, smoking, hypertension and
  diabetes; verify whether the published version includes or excludes LDL-C / Lp(a) before
  implementing).
- **Form**: an integer points score that stratifies patients into risk categories. Take the
  exact point allocations and category cutpoints from the paper's scoring table.
- **Implementation note**: a points score and a probability model are not directly comparable on
  every metric. Discrimination (C-statistic / AUC) is comparable; calibration and NRI require the
  points score to be mapped to an observed-risk estimate per category first. Document that mapping.

## DLCN (Dutch Lipid Clinic Network) criteria

- **Use**: a diagnostic-likelihood score for FH, not an event-prediction model. If it appears as
  a comparator for *event* prediction, that is a category error — flag it. DLCN predicts
  "is this FH", not "will this FH patient have an event".

## FAMCAT

- **Primary citation**: Weng SF, Kai J, Andrew Neil H, Humphries SE, Qureshi N. Improving
  identification of familial hypercholesterolaemia in primary care: derivation and validation
  of the familial hypercholesterolaemia case ascertainment tool (FAMCAT). *Atherosclerosis*
  2015;238(2):336-343.
- **Use**: case-ascertainment / diagnostic tool for FH in primary care — again a diagnostic
  model, not an event-prediction model. Same category-error caution as DLCN.

---

## Comparator-implementation checklist

For each comparator used in a head-to-head:

- [ ] Primary paper obtained; coefficient/points table located and cited by table number.
- [ ] Full published predictor set implemented (no silent omissions).
- [ ] Coefficients frozen to a committed file; comparator coded as a pure no-fitting function.
- [ ] Outcome definition and time horizon matched to the comparator's original, or recalibration
      documented.
- [ ] Comparator is an event-prediction model if the task is event prediction (not a diagnostic
      score used out of role).
- [ ] If any required predictor is unavailable in the data, the comparator is reported as a
      named approximation ("SAFEHEART-RE without Lp(a)") or NOT COMPUTABLE — never passed off as
      the full published score.

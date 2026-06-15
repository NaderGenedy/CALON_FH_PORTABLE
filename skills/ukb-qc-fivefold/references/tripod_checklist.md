# TRIPOD Checklist — for Prediction Model Manuscripts

A reduced checklist used by Agent 4 (Statistical Reproducer) and Agent 5 (Provenance Detector).
Reference: Collins et al. *Ann Intern Med* 2015;162:55-63.

## TRIPOD-required reporting items (Type 1–4)

### Title and abstract
- [ ] Identify study as prediction model development or validation (or both)
- [ ] Specify target population and outcome
- [ ] State key performance measures (AUC/c-statistic, calibration)

### Methods

**Source of data:**
- [ ] Describe study design (development cohort, validation cohort)
- [ ] Specify date of recruitment and outcome assessment

**Participants:**
- [ ] State setting (e.g., specialist clinic, primary care, biobank)
- [ ] Inclusion / exclusion criteria explicit
- [ ] Sample size at each stage

**Outcome:**
- [ ] Define outcome explicitly (ICD-10 codes, definitions)
- [ ] Time horizon if applicable
- [ ] Outcome adjudication described

**Predictors:**
- [ ] List every candidate predictor
- [ ] Define measurement units, transformations
- [ ] Handle missing data (multiple imputation / complete-case / etc.)

**Sample size:**
- [ ] Justify n per Riley et al. 2020 (≥10 events per parameter as a minimum)

**Missing data:**
- [ ] Approach for missing data described and justified

**Statistical analysis:**
- [ ] Model type (Cox / logistic / Elastic Net / XGBoost) specified
- [ ] Variable selection method described
- [ ] How model presented (equation, nomogram, calculator)
- [ ] Performance measures + how computed
- [ ] Calibration measures + how computed
- [ ] Internal validation method (split-sample, bootstrap, cross-validation)

### Results

**Participants:**
- [ ] Flow of participants (CONSORT-style figure)
- [ ] Baseline characteristics in development and validation cohorts

**Model development:**
- [ ] Final model parameters (coefficients with 95% CIs)
- [ ] Methods to identify functional form / non-linearity / interactions

**Model performance:**
- [ ] Discrimination (C-statistic / AUC with 95% CI)
- [ ] Calibration (intercept + slope + visual)
- [ ] Brier score (raw)
- [ ] Calibration plot

**Model updating:**
- [ ] Recalibration / re-estimation if applied
- [ ] Reasons for update

**Comparison:**
- [ ] Comparison with existing models or simple predictors
- [ ] DeLong p-value or equivalent for paired AUCs
- [ ] NRI / IDI if reclassification claimed

### Discussion

- [ ] Strengths and limitations
- [ ] Comparison with previous research
- [ ] Generalisability / external validity
- [ ] Implications for clinical practice

### Other

- [ ] Funding sources
- [ ] Conflicts of interest
- [ ] Ethical approval reference (REC number)

## Reproducibility-specific items

The user's CLAUDE.md emphasis on full reproducibility expands TRIPOD with:

- [ ] Every numerical claim traces to a CSV row (no hand-typed literals)
- [ ] Reproducer script that asserts manuscript value == CSV value within tolerance
- [ ] Provenance table: claim → file → row → cohort
- [ ] Multiple-testing correction for >10 simultaneous tests (Bonferroni or BH-FDR)
- [ ] E-value sensitivity for headline ORs/HRs (VanderWeele & Ding 2017; threshold > 2.0)
- [ ] Locked rerun discipline: methodological bugs trigger end-to-end re-run, not hot-patch

## Common TRIPOD failures

1. **Calibration slope reported without intercept** — both are needed
2. **AUC reported without 95% CI** — invalid
3. **No internal validation** — model performance optimistic
4. **NRI without specifying threshold set** — value is uninterpretable
5. **Sample size justification missing** — limits inference
6. **No declaration of model updating** — if recalibration applied, must say so explicitly

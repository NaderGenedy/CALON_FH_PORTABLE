# CALON-FH Clinical Risk Model — Build Specification & Traceability

**Author:** Dr Nader Genedy, University Hospital of Wales, Cardiff
**Build:** from scratch, 2026-05-13
**Runner:** `CALON_clinical_model.py`
**UKB Application:** 1002450

---

## 1. Purpose

A clinical cardiovascular risk model for the CALON-FH programme that predicts
**incident atherosclerotic cardiovascular disease (ASCVD)** from routinely
available clinical data: lipid profile, lipid ratios, prior ASCVD, familial
hypercholesterolaemia status, sex, smoking, type-2 diabetes, hypertension, and
vascular/cardiac imaging.

Every numerical claim this pipeline produces is traceable to a named raw UK
Biobank CSV and field ID. No hard-coded literals; no derived master files used
as primary input.

---

## 2. THE MACE-vs-HYPERTENSION CORRECTION (read first)

A previous version of this analysis drew its ASCVD outcome from
`ukb_dates_mace.csv`. **That file is misnamed.** It contains the UK Biobank
first-occurrence series **p131286, p131288, p131290, p131292, p131294,
p131296**, which code **ICD-10 I10–I15 — the hypertension family — not ASCVD.**

Proof (5,000-row sample of `ukb_dates_mace.csv`):

| Field | ICD-10 | Populated | Interpretation |
|---|---|---|---|
| p131286 | I10 essential hypertension | **41.6 %** | hypertension prevalence — NOT ASCVD |
| p131288 | I11 | 0.3 % | hypertensive heart disease |
| p131290 | I12 | 0.4 % | hypertensive renal disease |
| p131292 | I13 | 0.1 % | hypertensive heart+renal |
| p131294 | I15 | 0.0 % | secondary hypertension |
| p131296 | I20 angina | 7.4 % | the only ASCVD-relevant field here |

A model that labelled "p131286 has a date" as an ASCVD event would classify
**~42 % of the population as having ASCVD** — a catastrophic outcome-definition
error. `ukb_dates_mace.csv` is therefore **not used** by this pipeline except
for the death-date fields (p40000/p40001).

**Correct ASCVD outcome** is built from the HES ICD-10 diagnosis array
(`ukb_icd10_full.csv`, field p41270) using the code set below.

### ASCVD code set (used)

| Endpoint | ICD-10 prefixes |
|---|---|
| Ischaemic heart disease | I20, I21, I22, I23, I24, I25 |
| Ischaemic / unspecified stroke | I63, I64 |
| Transient ischaemic attack | G45 |
| Peripheral arterial disease | I70, I73, I74 |

### Codes explicitly EXCLUDED

| Excluded | Reason |
|---|---|
| I10, I11, I12, I13, I15 | hypertension — not atherosclerotic disease |
| I35 | aortic stenosis — structural valve disease, not ASCVD |
| I60, I61, I62 | haemorrhagic stroke — different mechanism (LDL J-curve) |
| I50 | heart failure — separate endpoint |

---

## 3. Raw data sources (inventory-verified 2026-05-13)

All twelve sources confirmed on disk by `CALON_CLINICAL_inventory.py`.

| Domain | File | Fields used |
|---|---|---|
| Assay lipids | `ukb_reviewer_longitudinal_lipids.csv` | p30780 LDL, p30690 TC, p30760 HDL, p30870 TG (i0 baseline) |
| ApoB / ApoA1 | `ukb_reviewer_apob_lpa.csv` | p30890 ApoB, p30900 ApoA1 (i0) |
| Lp(a) | `ukb_FULL_MASTER.csv` | lpa_chem (assay Lp(a)) |
| HES diagnoses | `ukb_icd10_full.csv` | p41270 ICD-10 array, p41271 dates |
| FH (genetic) | `ukb_carriers_FINAL.csv` | eid list of LDLR pathogenic-variant carriers |
| Sex | `calon_sex.csv` | p31 (0 = female, 1 = male) |
| Smoking + BP | `paper3_smoking_bp.csv` | p20116 smoking status, p4080 SBP, p4079 DBP |
| HbA1c / T2DM | `ukb_reviewer_hba1c.csv` | p30750 HbA1c, p2976 age-diabetes-diagnosed |
| Medications | `ukb_reviewer_medications.csv` | p6153/p6177 (BP & cholesterol meds) |
| Carotid IMT | `ukb_carotid_imt.csv` | p22671/p22674/p22677/p22680 (mean/max L/R, instance 2) |
| Cardiac MRI | `ukb_cardiac_mri.csv` | p22420–p22425 (LVEF, LV mass, EDV, ESV, …) |
| Recruitment | `ukb_recruitment_dates.csv` | p53 assessment date, p21022 age, p34 year of birth |

All UKB array columns carry a `participant.` prefix in the raw files; the
pipeline strips it on load.

---

## 4. Cohort definition

**Default cohort:** all UK Biobank participants with a complete baseline assay
lipid profile (LDL, HDL, TC, TG all non-missing at instance 0) and a valid
assessment-centre date. This is configurable in the script's `CONFIG` block.

**Why not carriers-only:** "FH" is used here as a binary *feature* (genetic
LDLR pathogenic-variant carrier — `ukb_carriers_FINAL.csv`), not a cohort
filter, so the cohort must be broader than carriers. Approximately 0.6 % of the
default cohort are FH carriers.

> Note: if "FH" was intended as *family history of premature CVD* (the QRISK /
> SCORE2 convention) rather than genetic FH, the UKB family-history fields
> p20107 / p20110 / p20111 would need extraction — they are not in the current
> inventory. The script uses genetic-carrier FH and flags this choice.

---

## 5. Outcome definition

**Primary outcome:** incident ASCVD = first HES ICD-10 record of any ASCVD code
(§2) occurring **after** the baseline assessment date.

- **Prior ASCVD** (any ASCVD code dated *before* baseline) becomes a model
  *feature*, not the outcome — this is the "previous ASCVD" predictor.
- Date source: p41271 (HES diagnosis dates). The script reports date coverage
  at runtime. If date coverage is inadequate for survival analysis, it falls
  back to a **prevalent-ASCVD logistic model** and flags the change explicitly.
- Death (p40000) is handled as a competing event / censoring point.

---

## 6. Feature specification

Every feature, its raw field, and its source file:

| # | Feature | Definition | Raw field | Source file |
|---|---|---|---|---|
| 1 | age | age at assessment | p21022 | ukb_recruitment_dates |
| 2 | sex | 1 = male | p31 | calon_sex |
| 3 | LDL-C | direct LDL, mmol/L | p30780_i0 | longitudinal_lipids |
| 4 | HDL-C | mmol/L | p30760_i0 | longitudinal_lipids |
| 5 | total cholesterol | mmol/L | p30690_i0 | longitudinal_lipids |
| 6 | triglycerides | mmol/L | p30870_i0 | longitudinal_lipids |
| 7 | ApoB | g/L | p30890_i0 | apob_lpa |
| 8 | ApoA1 | g/L | p30900_i0 | apob_lpa |
| 9 | Lp(a) | nmol/L | lpa_chem | ukb_FULL_MASTER |
| 10 | non-HDL-C | TC − HDL | derived | — |
| 11 | TC/HDL ratio | derived | derived | — |
| 12 | LDL/HDL ratio | derived | derived | — |
| 13 | TG/HDL ratio | derived (insulin-resistance proxy) | derived | — |
| 14 | ApoB/ApoA1 ratio | derived | derived | — |
| 15 | remnant cholesterol | TC − LDL − HDL | derived | — |
| 16 | prior ASCVD | ASCVD ICD-10 before baseline | p41270/p41271 | ukb_icd10_full |
| 17 | FH (genetic) | LDLR pathogenic-variant carrier | eid match | ukb_carriers_FINAL |
| 18 | current smoking | p20116 == 2 | p20116_i0 | smoking_bp |
| 19 | T2DM | ICD-10 E11 OR HbA1c ≥ 48 mmol/mol OR age-diabetes non-missing | p41270 / p30750 / p2976 | icd10 / hba1c |
| 20 | hypertension | ICD-10 I10–I15 OR SBP ≥ 140 OR DBP ≥ 90 OR BP medication | p41270 / p4080 / p4079 / p6153-p6177 | icd10 / smoking_bp / medications |
| 21 | carotid IMT (mean) | mean of 4 angles | p22671/4/7/80_i2 | carotid_imt |
| 22 | LVEF | left-ventricular ejection fraction | p22420_i2 | cardiac_mri |

Lipid ratios are computed in the script; each derived feature's formula is
recorded in the traceability ledger.

---

## 7. Model

- **Primary:** Cox proportional-hazards model for incident ASCVD, follow-up
  time as timescale, with a leak-free 70/30 train/test split and 5-fold
  cross-validation on the training set (per the `cox-analysis` methodology
  standard).
- **Companion:** logistic regression on the same features, reporting AUC, for
  discrimination comparison and because not all readers think in hazard ratios.
- **Reporting (TRIPOD):** C-statistic + 95 % CI, calibration slope + intercept,
  Brier score, hazard ratios with 95 % CI, proportional-hazards assumption
  check.
- **Pre-specified sensitivity variant:** NoAgeLDL — model refitted without age
  and LDL-C, to test whether the remaining features carry independent signal.

---

## 8. Traceability framework

The script writes `CALON_clinical_TRACEABILITY.csv` — one row per feature and
per result — recording: value, raw source file, raw field ID, derivation
formula (if derived), n non-missing, and a re-derivation check. Goal: every
number in any resulting manuscript can be traced file → field → row in one
lookup.

---

## 9. How to run

```
python CALON_clinical_model.py
```

Python 3.12. Requires: pandas, numpy, scipy, scikit-learn, lifelines,
statsmodels. Runtime ≈ 3–6 min (the 72 MB ICD-10 file is the bottleneck).

Outputs (all in the project directory):

| File | Contents |
|---|---|
| `CALON_clinical_analysis_dataset.csv` | the assembled per-participant feature + outcome matrix |
| `CALON_clinical_results.csv` | C-statistic, AUC, calibration, hazard ratios |
| `CALON_clinical_TRACEABILITY.csv` | per-feature / per-result provenance ledger |
| `CALON_clinical_RUN_LOG.txt` | full console log of the run |

---

## 10. Known limitations (honest)

1. HES diagnosis dates (p41271) may be incompletely populated; the script
   reports coverage and degrades gracefully to a prevalent-ASCVD model if
   survival analysis is not supportable.
2. "FH" is genetic LDLR-carrier status; if family-history-of-CVD was intended,
   fields p20107/p20110/p20111 must be extracted.
3. Imaging features (carotid IMT, cardiac MRI) are available only for the
   imaging sub-study (~50,000 participants); the script reports imaging
   coverage and offers a complete-case-with-imaging sensitivity model.
4. Lipid values are baseline and may be on-treatment; a treatment-adjustment
   step is out of scope for this clinical model and noted as future work.

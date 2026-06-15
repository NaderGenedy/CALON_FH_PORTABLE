# TUDOR Analysis Run Report

Date: 2026-04-05

## Scripts Run

1. `REPRODUCE_TUDOR_DIAGNOSTIC.R`
2. `tudor_loco_output/TUDOR_flaws_analysis.R`
3. `TUDOR_UKB_LIPID_CLINIC.R`

R runtime used:

- `C:/Program Files/R/R-4.5.2/bin/Rscript.exe`

## Key Results

### 1. Main diagnostic reproduction

Source file:

- `TUDOR_reproduced_AUC.csv`
- `TUDOR_reproduced_coefficients.csv`

Results:

- PASS AUC: `0.7725` (95% CI `0.7594-0.7863`)
- LOCO fold 1 AUC: `0.8412`
- LOCO fold 2 AUC: `0.7725`
- LOCO average AUC: `0.8068`
- DLCN AUC: `0.6896`
- LDL alone AUC: `0.5915`
- Trig_Filter alone AUC: `0.7274`
- Youden sensitivity: `62.9%`
- Youden specificity: `82.7%`

Reproduced coefficient file contains:

- `(Intercept) 1.3932`
- `tendon_xanth 1.0352`
- `hdl -0.9628`
- `trig_filter 0.4832`
- `on_statin 0.4392`
- `corneal_arcus 0.3951`
- `sex -0.3933`
- `age -0.0713`
- `ldl_ut 0.0528`
- `tg -0.0509`

### 2. Flaw remediation analysis

Source files:

- `tudor_loco_output/tudor_flaw_results/flaw_remediation_status.csv`
- `tudor_loco_output/tudor_flaw_results/index_vs_cascade_auc.csv`
- `tudor_loco_output/tudor_flaw_results/ukb_recalibration_summary.csv`
- `tudor_loco_output/tudor_flaw_results/roc_operating_points.csv`

Results:

- UKB recalibration: done
- NRI/IDI bootstrap: done
- Missing data audit: done
- Index vs cascade AUC split: done
- ApoB coverage and cost-effectiveness: done
- ROC operating thresholds: done

Still open per remediation status:

- F5 ethnicity stratification
- F6 full TRG multivariable adjustment with T2DM/BMI/alcohol
- F7 premature ASCVD with event dates
- F10 cumulative/untreated LDL revalidation

### 3. UKB lipid clinic validation script

Source files:

- `tudor_loco_output/ukb_lipid_clinic_discrimination.csv`
- `tudor_loco_output/ukb_lipid_clinic_delong.csv`

Results:

- UKB lipid clinic cohort size: `271819`
- FH positive in that cohort: `1334`
- FH prevalence: `0.49%`
- TUDOR AUC: `0.6881` (95% CI `0.6731-0.7037`)
- eDLCN AUC: `0.7259`
- LDL-C alone AUC: `0.7461`
- Trig_Filter AUC: `0.6785`
- TUDOR vs eDLCN DeLong p: `7.71e-06`
- Calibration intercept: `-6.014`
- Calibration slope: `7.074`
- Brier score: `0.0127`

## Critical Interpretation

The project currently contains multiple incompatible analysis states:

1. The manuscript/rebuttal state uses UKB figures such as `729 FH cases`, UKB AUC `0.750`, and calibration slope `6.33`.
2. `lancet_statistics_summary.csv` reports UKB AUC `0.7873` and calibration slope `12.395`.
3. The freshly rerun `TUDOR_UKB_LIPID_CLINIC.R` reports a much broader UKB lipid clinic cohort (`271819`) with TUDOR AUC `0.6881`, worse than both eDLCN and LDL-C alone.

The package is therefore not on one frozen submission-ready analysis state.

## Recommendation

Before any submission update:

1. Freeze the intended UKB cohort definition.
2. Freeze the intended TUDOR coefficient set.
3. Regenerate manuscript, response letter, tables, and summary CSVs from that one state only.

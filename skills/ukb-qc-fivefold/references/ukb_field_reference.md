# UKB Field Reference — Cardiology FH Work

Canonical map of UKB fields used in Cardiff CALON / TUDOR / NB01 / NB02 pipelines.
UKB Application ID: 1002450.

## Lipids — assay-based
- p30780_i0/i1 — LDL direct (mmol/L, 0.5-15)
- p30690_i0/i1 — Total cholesterol (1.5-20)
- p30760_i0/i1 — HDL (0.3-5)
- p30870_i0/i1 — Triglycerides (0.2-30)
- p30890_i0/i1 — ApoB (g/L, 0.3-3.0)
- p30900_i0/i1 — ApoA1
- p30790_i0/i1 — Lp(a) (nmol/L, 0-400)

## Lipids — NMR (Nightingale)
- p23400_i0 — Total cholesterol (NMR)
- p23404/p23405_i0 — Clinical LDL / LDL (NMR)
- p23406_i0 — HDL (NMR)
- p23410_i0 — ApoB (NMR)
- p23450-p23457_i0 — HDL sub-fractions 1-4

**TRAP:** NMR LDL (p23404/p23405) != assay LDL (p30780). Do not mix.

## Outcomes — first-occurrence dates
- p131296 — I20 unstable angina (CORRECT ASCVD)
- p131298 — I21 MI (CORRECT ASCVD)
- p131306 — I25 chronic IHD (CORRECT ASCVD)
- p131346 — I63 cerebral infarction (stroke)
- p131308 — I50 heart failure

**TRAP:** p131286-p131294 code I10-I15 (hypertension) — NOT ASCVD.

## Death
- p40000_i0 — Date of death
- p40001_i0 — Underlying cause
- p40002_i0_a0-a9 — Contributory causes (10 slots)
- p40007 — Age at death

## Diagnoses
- p41270 — ICD-10 diagnoses array
- p41271 — ICD-10 diagnosis dates
- p41272/p41273 — OPCS-4 codes/dates

## Medications
- p6153_i0/i1/i2 — Cholesterol/BP/DM (female)
- p6177_i0/i1/i2 — Same (male)
- p42039 — GP prescription linkage flag

For statin TYPE and DOSE: tier-2 GP prescriptions (BNF/drug name/quantity) at D:/Projects/CALON_AlphaFold_Rebuild/Paper3_ASCVD_Prediction/data.csv.

## Demographics
p21022 (age) | p31 (sex 0F 1M) | p34 (YOB) | p52 (MOB) | p53_i0/i1 (visit dates) | p21000_i0/i1 (ethnicity) | p21001_i0/i1 (BMI)

## Vitals
p4080_i0 (SBP) | p4079_i0 (DBP)

## Lifestyle
p1558_i0 (alcohol) | p20116_i0/i1/i2 (smoking 3 instances)

## Inflammation / metabolic
p30710 (CRP) | p30750 (HbA1c) | p2976_i0 (age DM dx) | p2443_i0 (self-report DM)

## Genetics
p26206 (LDL-PRS) | p26228 (CAD-PRS) | p30105 (CHIP) | p30106-p30107 (CHIP variants)

## Imaging
p22420-p22425_i2 — Cardiac MRI (LVEF, mass, EDV, ESV, LA, dist.)
p22671/22674/22677/22680_i2 — Carotid IMT
p21088_i2/i3 — Liver PDFF
p21085_i2 — Trunk fat
p21086_i2 — Visceral adipose

## Reproductive
p3581_i0 (menarche) | p2814_i0 (hyster.) | p2724_i0 (menopause)

## Deprivation
p22189 (Townsend) | p26410 (IMD England) | p26426 (IMD Wales)

## Files on disk (D: drive)
```
D:/Projects/CALON_AlphaFold_Rebuild/data/
├── ukb_lpa.csv
├── ukb_nmr_batch1.csv / batch2 / batch3
├── ukb_nmr_hdl_subfractions.csv
├── ukb_prs.csv
├── ukb_dates_mace.csv
├── ukb_icd10_full.csv
├── ukb_gp_prescriptions.csv
├── ukb_cardiac_mri.csv
├── ukb_carotid_imt.csv
├── ukb_chip.csv
├── ukb_liver_pdff.csv
├── ukb_visceral_fat.csv
├── ukb_menopause.csv
├── ukb_recruitment_dates.csv
└── ukb_reviewer_*.csv (lipids, ApoB/Lp(a), CRP, HbA1c, demographics, meds, smoking, death, deprivation, ICD-10)
```

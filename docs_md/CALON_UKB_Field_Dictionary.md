# CALON-FH: UK Biobank Field Dictionary
## Complete Variable Mapping for External Validation

**Project:** CALON-FH External Validation in UK Biobank
**Purpose:** Predict ASCVD in genetically confirmed FH (≠ TUDOR which diagnoses FH)
**Comparator:** SAFEHEART-RE / SAFEHEART-UK
**Date:** February 2026

---

## 1. CORE CALON MODEL VARIABLES

| CALON Variable | UKB Field | UKB Description | Units | Derivation |
|---|---|---|---|---|
| **Age** | p21022 | Age at recruitment | Years | Direct |
| **Reverse-engineered LDL-C** | p30780 + p20003 | LDL direct ÷ residual factor | mmol/L | LDL / RF (from statin codes) |
| **Lipid-year exposure** | Derived | RE_LDL × (age_tx - 18) / 10 | mmol-years/10 | Cholesterol pack-years |
| **Log(ApoB/LDL-C)** | p30640 / p30780 | log(ApoB ÷ LDL) | dimensionless | Particle discordance |
| **Lp(a) elevated** | p30790 | Lipoprotein(a) > 125 nmol/L | Binary | Threshold at 125 nmol/L |
| **Inverse HDL-C** | p30760 | 1 / HDL-C | 1/(mmol/L) | Reciprocal transform |

## 2. SAFEHEART-RE BENCHMARK VARIABLES

| SAFEHEART Variable | UKB Field | UKB Description | Units | Notes |
|---|---|---|---|---|
| Age | p21022 | Age at recruitment | Years | Direct |
| Male sex | p31 | Sex | Binary | 0=F, 1=M |
| LDL-C | p30780 | LDL direct | mmol/L | ON-TREATMENT value |
| Hypertension | p20002 + p4080 | Self-report + SBP≥140 | Binary | Composite definition |
| BMI | p21001_i0 | Body mass index | kg/m² | Visit 0 |
| Current smoking | p20116_i0 | Smoking status | Binary | Code 2 = current |
| Prior CVD | p131298 + p6150 | MI/Stroke before baseline | Binary | Prevalent ASCVD |

## 3. COMPLETE FIELD LIST BY BATCH

### Batch 1: Demographics & Anthropometrics
| Field | Description | Category |
|---|---|---|
| p31 | Sex | Demographics |
| p21022 | Age at recruitment | Demographics |
| p34 | Year of birth | Demographics |
| p52 | Month of birth | Demographics |
| p53_i0 | Date of assessment (Visit 0) | Date |
| p53_i1 | Date of assessment (Visit 1) | Date |
| p21001_i0 | BMI Visit 0 | Anthropometric |
| p21001_i1 | BMI Visit 1 | Anthropometric |
| p48_i0 | Waist circumference Visit 0 | Anthropometric |
| p48_i1 | Waist circumference Visit 1 | Anthropometric |
| p49_i0 | Hip circumference Visit 0 | Anthropometric |
| p4080_i0 | Systolic BP (2 readings) | Blood pressure |
| p4079_i0 | Diastolic BP (2 readings) | Blood pressure |
| p22006 | Genetic ethnic grouping | Genetics |
| p22001 | Genetic sex | Genetics |

### Batch 2: Lipid Biomarkers
| Field | Description | Units | CALON Role |
|---|---|---|---|
| p30690_i0 | Total cholesterol | mmol/L | Descriptive |
| p30760_i0 | HDL cholesterol | mmol/L | **inv_hdl = 1/HDL** |
| p30780_i0 | LDL direct | mmol/L | **RE_LDL denominator** |
| p30870_i0 | Triglycerides | mmol/L | Descriptive |
| p30640_i0 | **Apolipoprotein B** | g/L | **ApoB/LDL ratio** |
| p30630_i0 | Apolipoprotein A | g/L | inv_ApoA (alternative) |
| p30790_i0 | **Lipoprotein(a)** | nmol/L | **Lp(a) >125 binary** |

### Batch 3: Other Biomarkers
| Field | Description | Units | Purpose |
|---|---|---|---|
| p30750_i0 | HbA1c | mmol/mol | Diabetes phenotyping |
| p30740_i0 | Glucose | mmol/L | Metabolic status |
| p30710_i0 | C-reactive protein | mg/L | Inflammation |
| p30620_i0 | ALT | U/L | Liver function |
| p30650_i0 | AST | U/L | Liver function |
| p30680_i0 | Calcium | mmol/L | Metabolic |
| p30700_i0 | Creatinine | umol/L | eGFR calculation |
| p30720_i0 | Cystatin C | mg/L | eGFR (CKD-EPI) |

### Batch 4: Medications (Field 20003)
| Code | Medication | Statin Intensity | Residual Factor |
|---|---|---|---|
| 1140888648 | **Atorvastatin** | High | 0.50 |
| 1140910632 | **Rosuvastatin** | High | 0.50 |
| 1140861958 | Simvastatin | Moderate | 0.65 |
| 1141146234 | Pravastatin | Moderate | 0.65 |
| 1141192414 | Fluvastatin | Low-Moderate | 0.70 |
| 1140861922 | Lovastatin | Moderate | 0.65 |
| 1141146138 | **Ezetimibe** | Add-on | ×0.80 |

**CRITICAL LIMITATION:** UKB does NOT record statin **dose**. Assignment of high vs moderate intensity is by drug name only. Atorvastatin/Rosuvastatin are assumed high-intensity; others moderate. Sensitivity analysis with ±10% RF adjustment addresses this.

### Batch 5: Smoking & Alcohol
| Field | Description | Codes |
|---|---|---|
| p20116_i0 | Smoking status | 0=Never, 1=Previous, 2=Current |
| p20161 | Pack years of smoking | Continuous |
| p2867_i0 | Age started smoking | Years |
| p2897_i0 | Age stopped smoking | Years |
| p20117_i0 | Alcohol drinker status | Categorical |
| p1558_i0 | Alcohol intake frequency | Categorical |

### Batch 6: Comorbidities
| Field | Description | Codes |
|---|---|---|
| p2443_i0 | Diabetes (doctor-diagnosed) | 0=No, 1=Yes |
| p6177_i0 | Medication for BP/cholesterol (Male) | 1=Cholesterol, 2=BP, 3=Insulin |
| p6153_i0 | Medication for BP/cholesterol (Female) | Same codes |
| p2966_i0 | Age hypertension diagnosed | Years |

### Batch 7: ASCVD First-Occurrence Dates (GOLD STANDARD)
| Field | ICD-10 | Description | Event Type |
|---|---|---|---|
| p131296 | I20 | Angina pectoris | Coronary |
| p131298 | **I21** | **Acute MI** | Coronary |
| p131300 | I22 | Subsequent MI | Coronary |
| p131306 | I25 | Chronic IHD | Coronary |
| p131364 | **I63** | **Cerebral infarction** | Stroke |
| p131366 | I64 | Stroke NOS | Stroke |
| p131380 | I70 | Atherosclerosis | PVD |
| p131386 | I73 | Other PVD | PVD |
| p131388 | I74 | Arterial embolism | PVD |
| p40000_i0 | — | Date of death | Death |
| p40001_i0 | — | Cause of death (ICD-10) | Death |

### Batch 8: Self-Reported ASCVD
| Field | Description | Codes |
|---|---|---|
| p6150_i0 | Vascular/heart problems | 1=MI, 2=Angina, 3=Stroke, 4=HBP |
| p3894_i0 | Age MI diagnosed | Years |
| p3627_i0 | Age angina diagnosed | Years |
| p4056_i0 | Age stroke diagnosed | Years |

### Batch 9: Carotid IMT (Imaging Subset ~100k)
| Field | Description | Notes |
|---|---|---|
| p22671-p22678_i2 | Mean carotid IMT at various angles | Imaging visit only |

**NOTE:** Available for ~100,000 imaging participants. Expect ~200-500 FH patients with CIMT data.

### Batch 10: Coronary Artery Calcium
**IMPORTANT:** Standard UKB imaging is cardiac MRI, NOT CT. CAC scores are NOT routinely available. If derived CAC variables exist, they are in Category 157 derived imaging fields. This will be listed as a limitation.

### Batch 11: Self-Reported Illness (Field 20002)
| Code | Condition | CALON Relevance |
|---|---|---|
| 1065 | Hypertension | HTN phenotyping |
| 1075 | Heart attack / MI | Prevalent ASCVD |
| 1066 | Heart failure | Comorbidity |
| 1067 | Angina | Prevalent ASCVD |
| 1081 | Stroke | Prevalent ASCVD |
| 1082 | TIA | Cerebrovascular |
| 1083 | PVD | Prevalent ASCVD |
| 1223 | Type 2 diabetes | DM phenotyping |
| 1473 | High cholesterol | Treatment proxy |

### Batch 12-13: HES ICD-10 & OPCS-4
| Source | Fields | Content |
|---|---|---|
| HES ICD-10 codes | p41270 (arrays 0-149) | Hospital diagnoses |
| HES ICD-10 dates | p41280 (arrays 0-149) | Dates of diagnoses |
| OPCS-4 codes | p41272 (arrays 0-59) | Operative procedures |
| OPCS-4 dates | p41282 (arrays 0-59) | Dates of procedures |

**ASCVD Composite (OPCS-4 Revascularisation):**
- PCI: K49, K50, K75 series
- CABG: K40-K46 series

---

## 4. KNOWN LIMITATIONS

| Limitation | Impact | Mitigation |
|---|---|---|
| No statin **dose** in UKB | RF assignment by drug name only | Sensitivity ±10% RF |
| No age at treatment **start** | Cannot calculate exact lipid-years | 3 scenarios (A/B/C) |
| No CAC score (routine) | Cannot validate CAC substudy | Listed as limitation |
| CTCA not available | No coronary anatomy data | Listed as limitation |
| No index/cascade distinction | All UKB FH are probands | Acknowledged; compare with Wales where index/cascade known |
| PCSK9i not available at baseline | Cannot adjust for PCSK9i | <1% of FH on PCSK9i at baseline (pre-2015) |
| Carotid IMT in subset only | ~3-5% of FH will have CIMT | Exploratory analysis |
| FH ascertainment differs from Wales | UKB = population screening; Wales = clinical service | Discuss transportability |

---

## 5. GENETIC FH IDENTIFICATION

FH patients are identified from Whole Exome Sequencing (WES) data using the same pipeline as TUDOR:

1. **Variant calling:** UK Biobank WES processed through GATK pipeline
2. **Gene targets:** LDLR, APOB, PCSK9
3. **Classification:** ClinVar pathogenic / likely pathogenic variants
4. **Output:** `is_fh_genetic` flag in TUDOR_UKB_Features.csv
5. **Expected yield:** ~1,500-2,500 genetically confirmed FH from ~470,000 WES participants

---

## 6. ASCVD OUTCOME DEFINITION

**Primary composite endpoint:** First occurrence of ANY of:
- Myocardial infarction (I21, I22)
- Ischaemic stroke (I63)
- Coronary revascularisation (PCI: K49/K50/K75; CABG: K40-K46)
- Peripheral vascular disease (I70, I73, I74)

**Time zero:** Date of UKB assessment (p53_i0)
**Censoring:** Min(death date, latest HES linkage, loss to follow-up)

**Prevalent events:** ASCVD date ≤ assessment date
**Incident events:** ASCVD date > assessment date

---

*This dictionary accompanies the CALON-FH UKB extraction pipeline (00_CALON_extract_ukbrap.sh)*

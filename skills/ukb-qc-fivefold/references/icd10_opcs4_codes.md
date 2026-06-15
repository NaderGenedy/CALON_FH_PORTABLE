# ICD-10 / OPCS-4 Code Reference — Cardiology Endpoints

Complete code lists for the user's cardiology research domain.
Reference for Agent 1 (Raw Data Integrity) audit checks.

## ASCVD composite (CALON / FH / TUDOR / NB02 standard)

**ICD-10 codes:**
- I20 Angina pectoris
- I21 Acute myocardial infarction
- I22 Subsequent MI
- I23 Complications post-MI
- I24 Other acute IHD
- I25 Chronic IHD
- I63 Cerebral infarction
- G45 Transient ischaemic attack
- I70 Atherosclerosis (peripheral arterial disease)
- I73 Other peripheral vascular disease
- I74 Arterial embolism and thrombosis

**Must NOT be included in ASCVD:**
- I35 — aortic stenosis (this is structural valve disease, not atherosclerotic)
- I50 — heart failure (separate endpoint)
- All I60–I62 — haemorrhagic stroke (different mechanism)

**UKB first-occurrence date fields (correct ones for ASCVD):**
- p131296 — I20 (unstable angina)
- p131298 — I21 (MI)
- p131306 — I25 (chronic IHD)

**UKB first-occurrence fields that should NOT be used for ASCVD:**
- p131286–p131290 — these code I10/I11/I12 (essential, secondary hypertension)
- p131292/p131294 — code I13/I15 (hypertensive heart/kidney disease)

## Coronary revascularisation — OPCS-4

**CABG (Coronary Artery Bypass Graft):**
- K40 Saphenous vein graft replacement
- K41 Saphenous vein graft replacement
- K42 Other autologous graft
- K43 Other allograft
- K44 Other prosthesis
- K45 Combination CABG
- K46 Other unspecified CABG

**PCI (Percutaneous Coronary Intervention):**
- K49 Transluminal balloon angioplasty (PTCA)
- K50 Coronary intervention with stent
- K75 Other percutaneous intervention

## Severe aortic stenosis intervention

- K261, K262, K263 — Surgical AVR
- **K611 — Balloon valvuloplasty (TRAP: commonly omitted from severe-AS lists, loses ~3,000 cases in UKB)**
- K268, K269 — Other valve replacements

## Heart failure

- I50.0 Congestive heart failure
- I50.1 Left ventricular failure
- I50.9 Heart failure, unspecified
- I11.0 Hypertensive heart disease with HF
- I13.0 Hypertensive heart and renal with HF
- UKB first-occurrence: p131308

## Stroke specific

**Ischaemic stroke:**
- I63 Cerebral infarction
- I65 Pre-cerebral artery occlusion (no infarct)
- I66 Cerebral artery occlusion (no infarct)
- UKB first-occurrence: p131346 / p131350

**Haemorrhagic stroke:**
- I60 Subarachnoid haemorrhage
- I61 Intracerebral haemorrhage
- I62 Other non-traumatic intracranial haemorrhage

**TIA:**
- G45 (all subcategories)

## Atrial fibrillation

- I48 Atrial fibrillation and flutter
- I48.0–I48.9 (all subcategories)

## Familial Hypercholesterolaemia (rare ICD-10 use; mostly via genetic confirmation)

- E78.0 Pure hypercholesterolaemia
- E78.2 Mixed hyperlipidaemia (some misclassified FH)
- Z83.4 Family history of other endocrine, nutritional and metabolic diseases

## Diabetes

**Type 2 diabetes mellitus:**
- E11 (all subcategories)
- UKB self-report: p2443_i0 == 1
- UKB HbA1c-based: p30750 ≥ 48 mmol/mol

## Trap summary — recurring failure modes

| Trap | What happens | How to avoid |
|---|---|---|
| K611 missing from severe-AS | Loses 30%+ of cases (~3,000 in UKB) | Always include in valve-procedure lists |
| I35 in ASCVD composite | Inflates events with non-atherosclerotic disease | Exclude I35 from ASCVD |
| p131286/8/90/2/4 used as ASCVD | Codes hypertension, not ASCVD | Use only p131296, p131298, p131306 |
| Haemorrhagic + ischaemic stroke combined | Confounds direction of LDL effect (J-curve) | Separate I63/G45 from I60-I62 |
| I50 in ASCVD | Heart failure has different mechanism | Keep HF as separate endpoint |

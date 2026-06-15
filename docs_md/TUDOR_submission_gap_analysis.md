# TUDOR Manuscript — Submission Gap Analysis
# Date: 2026-03-29 | Target: Expert Review of Clinical Lipidology
# Status: Pre-submission review

---

## ✅ COMPLETED (Evidence in hand, results computed)

| # | Item | Evidence | Location |
|---|------|----------|----------|
| F2 | UKB recalibration | Slope 1.000 post-correction; original slope 12.4 | ukb_recalibration_summary.csv |
| F4 | NRI/IDI bootstrap | SW NRI_events +0.173 (0.014–0.338) p=0.033; IDI +0.119 | nri_idi_bootstrap.csv |
| F6 | TRG shield adjusted | T2DM β=−0.477 p=8.2×10⁻⁷⁰; OR 1.82 (age+sex adj) | trgshield_*.csv |
| F9 | Missing data audit | ApoB 6.3% overall; 62.7% in FH+; DLCN 1.7% | missing_data_audit.csv |
| F11 | Index vs cascade AUC | Index 0.807 vs cascade 0.708 (Wales) | index_vs_cascade_auc.csv |
| F12 | ApoB coverage | 50% LC cohort; Friedewald RMSE 2.86 too imprecise | apob_coverage.csv |
| F13 | Operating thresholds | 3 thresholds, PPV heatmap, LR+ table | roc_operating_points.csv |
| F5 | Ethnicity distribution | 96.9% White; South Asian T2DM 21.4%, ASCVD 19% | ethnicity_lc_cohort_summary.csv |
| — | 5 R figures | All rendered as PDF + PNG (Nature specs) | Figure1–5 in flaw_results/ |
| — | 4 manuscript tables | All 4 tables in TUDOR_Tables.docx | TUDOR_Tables.docx |
| — | Canva Pro prompts | 5 detailed design prompts | TUDOR_canva_pro_prompts.md |

---

## 🔴 CRITICAL GAPS — Cannot submit without these

### GAP 1: Premature ASCVD with event DATES (Flaw 7 — incomplete)
**What's missing:** ICD-10 field `p41280` (date of diagnoses) was not extracted from UKB.
**Why it matters:** "Premature ASCVD in FH" requires age-at-event < 55, not age-at-recruitment.
Currently we have 3.0% premature ASCVD proxy (ASCVD present AND age <55 at recruitment), which is CONSERVATIVE — misses patients who were older at recruitment but had early events.
**Action required:**
1. Re-run `tudor_ukb_extract_supplement.sh` on RAP adding `participant.p41280` to FIELDS
2. Run `tudor_supplement_process.R` with event dates to compute true age-at-first-ASCVD
3. Report: % with event < 55, median age at first event, Kaplan-Meier by TUDOR tertile

**Estimated time:** 1–2 hours (RAP extraction + local processing)

---

### GAP 2: Ethnicity-stratified AUC (Flaw 5 — incomplete)
**What's missing:** No EID column in `loco_predictions_complete.csv` → cannot link TUDOR scores to ethnicity from UKB supplement.
**Why it matters:** Reviewer will ask: "Does TUDOR work in South Asian patients?" We only have ethnicity distribution, not performance metrics by ethnicity.
**Action required — Option A (preferred):**
1. Re-run `TUDOR_UKB_LIPID_CLINIC.R` on RAP, modify to output `eid` column in predictions
2. Merge with supplement ethnicity data
3. Compute AUC stratified by ethnicity (expect lower in South Asian due to higher T2DM confounding)

**Action required — Option B (fast):**
Write in limitations: "Ethnicity-stratified AUC validation is not yet available as EID linkage to TUDOR predictions requires re-extraction; this is a priority for the next analysis cycle."

---

### GAP 3: Flaw 1 — CALON Predictive Model (external validation not complete)
**What's missing:** The CALON structural biology pipeline (AlphaFold3 → FoldX → Rosetta) as TUDOR feature is described but the model performance with structural features incorporated is not validated.
**Why it matters:** The paper claims structural biology features improve TUDOR; this claim needs AUC comparison.
**Action required:**
- Either report TUDOR WITHOUT structural features (current AUC 0.834) and describe CALON-structural as future work
- OR complete the CALON FoldX/Rosetta pipeline and rerun TUDOR including ddG features

---

### GAP 4: TUDOR prediction equation / calculator
**What's missing:** No published formula or online calculator.
**Why it matters:** Expert Review reviewers will ask "how would a clinician compute this?"
**Action required:**
- Extract TUDOR model coefficients (gradient boosting feature importances + logistic wrapper)
- OR provide simplified logistic regression approximation with coefficients
- Ideally: simple online calculator at a URL

---

## 🟡 MODERATE GAPS — Should address but can submit with caveats

### GAP 5: Flaw 10 — Cholesterol-years underpowered
**Status:** Only 3.3% (n=1,205/37,050) of LC cohort have 2 LDL timepoints.
**Write in paper:** "Cholesterol-years analysis was limited to 3.3% of the lipid-clinic cohort with repeat measures (median gap 6 years). Definitive cumulative LDL burden analysis requires primary care record linkage (GP lipid data in the CPRD/SAIL Databank)."
**No further analysis needed.**

### GAP 6: Flaw 3 — Non-FH comparator group definition
**Status:** Non-FH group includes all UKB without confirmed variant — may include undiagnosed FH.
**Write in paper:** "The non-FH group is defined by absence of confirmed pathogenic variant; some participants may carry variants below current detection thresholds. This biases TUDOR performance metrics toward underestimation of true discriminatory ability."

### GAP 7: ApoB Imputation — informative missingness
**Status:** ApoB missing in 62.7% of FH+ — almost certainly informative (FH+clinics where ApoB wasn't measured vs specialist centres).
**Write in paper:** "Informative missing ApoB (62.7% in FH+ vs 4.7% in FH−) precludes multiple imputation. TUDOR performs equivalently with and without ApoB substitution, suggesting LDL_untreated and the TRG Filter capture sufficient variance."

---

## 🟢 OPTIONAL ENHANCEMENTS — Would strengthen but not required

### OPT 1: Decision Curve Analysis (DCA)
A net benefit plot comparing TUDOR vs DLCN vs "treat all" vs "treat none" at threshold probabilities 0.01–0.30 would satisfy TRIPOD reporting standard and any methodology-focused reviewer.
**Code:** Use `rmda` package in R. All data available.

### OPT 2: External validation in a non-UK population
TUDOR is trained on a UK cohort. Reviewers from European or international lipid societies may ask about transferability. If you have access to any Dutch/Spanish/Israeli FH registry data, even a small external validation (n>200) would substantially strengthen the paper.

### OPT 3: TRIPOD-AI / PROBAST checklist
Expert Review of Clinical Lipidology expects TRIPOD checklist. A completed TRIPOD table (as supplementary) would pre-empt reviewer comment.

### OPT 4: Sensitivity analysis — treating the UKB FH+ group
The UKB FH+ (n=938) are likely on statin (64% overall). Check if excluding UKB FH+ cases (leaving only Wales for outcome) changes AUC materially.

### OPT 5: Lp(a) as additional feature
Lp(a) was measured in 11% of LC cohort (n≈4,000). Could add Lp(a) as a covariate to see if it improves discrimination in that subgroup.

---

## SUBMISSION CHECKLIST

| Item | Status |
|------|--------|
| Abstract ≤250 words | ✅ Written |
| Introduction | ✅ Written |
| Methods (TRIPOD-compliant) | ✅ Written |
| Results with all tables/figures | ✅ Written |
| Discussion + limitations | ✅ Written |
| 15 Vancouver references | ✅ Complete |
| Table 1 — Baseline characteristics | ✅ TUDOR_Tables.docx |
| Table 2 — Performance metrics | ✅ TUDOR_Tables.docx |
| Table 3 — TRG Shield | ✅ TUDOR_Tables.docx |
| Table 4 — Operating characteristics | ✅ TUDOR_Tables.docx |
| Figure 1 — ROC + Calibration | ✅ Figure1_ROC_Calibration.pdf/png |
| Figure 2 — TRG Shield | ✅ Figure2_TRGShield.pdf/png |
| Figure 3 — Clinical operating | ✅ Figure3_ClinicalOperating.pdf/png |
| Figure 4 — ASCVD + NRI | ✅ Figure4_ASCVD_NRI.pdf/png |
| Figure 5 — Penetrance + architecture | ✅ Figure5_Penetrance_Architecture.pdf/png |
| 5 Canva Pro design prompts | ✅ TUDOR_canva_pro_prompts.md |
| Graphical abstract | 🔴 Pending (use Canva Prompt 1) |
| EID merge for ethnicity AUC | 🔴 GAP 2 above |
| Premature ASCVD with dates | 🔴 GAP 1 above |
| TUDOR score calculator | 🔴 GAP 4 above |
| TRIPOD checklist | 🟡 Recommended |
| Cover letter | 🟡 Write last |
| Ethics statement (UKB application no.) | 🟡 Add application number |

---

## PRIORITY NEXT STEPS (ordered by impact on acceptance)

1. **[HIGH]** Re-extract `p41280` on RAP → compute true premature ASCVD age-at-event
2. **[HIGH]** Re-run TUDOR pipeline on RAP saving `eid` → ethnicity-stratified AUC
3. **[HIGH]** Extract/describe TUDOR prediction formula for clinical use
4. **[MEDIUM]** Add DCA (decision curve analysis) figure — 2 hours with existing data
5. **[MEDIUM]** Complete TRIPOD checklist
6. **[LOW]** Write cover letter targeting Expert Review clinical scope

---

## FILES READY FOR SUBMISSION

```
C:/Users/nader/Downloads/calon_ukb_pipeline/
├── TUDOR_Manuscript.docx           ← Main manuscript (Word)
├── TUDOR_Tables.docx               ← 4 tables (standalone)
├── TUDOR_canva_pro_prompts.md      ← 5 Canva Pro design specs
├── TUDOR_submission_gap_analysis.md ← This file
├── tudor_loco_output/tudor_flaw_results/
│   ├── Figure1_ROC_Calibration.pdf/png
│   ├── Figure2_TRGShield.pdf/png
│   ├── Figure3_ClinicalOperating.pdf/png
│   ├── Figure4_ASCVD_NRI.pdf/png
│   ├── Figure5_Penetrance_Architecture.pdf/png
│   ├── Table1_Baseline.csv
│   ├── Table2_Performance.csv
│   ├── Table3_TRGShield.csv
│   └── Table4_OperatingCharacteristics.csv
├── TUDOR_figures_nature.R          ← Reproducible figure code
├── TUDOR_tables_nature.R           ← Reproducible table code
└── TUDOR_write_manuscript.R        ← Manuscript generation code
```

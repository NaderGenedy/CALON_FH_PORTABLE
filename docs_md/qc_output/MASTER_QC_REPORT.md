# Master QC Report -- five-agent audit
**Project:** `C:\Users\nader\Downloads\calon_ukb_pipeline`
**Generated:** 2026-05-12T12:39:25

## Status summary

| # | Agent | Status | Pass | Drift | Fail |
|---|---|---|---|---|---|
| 1. Raw Data | Raw Data Integrity | [!!] DRIFT | 0 | 0 | 0 |
| 2. Cohort | Cohort Definition Validator | ? ? | 0 | 0 | 0 |
| 3. Features | Feature Engineering Auditor | ? ? | 0 | 0 | 0 |
| 4. Statistics | Statistical Reproducer | ? ? | 0 | 0 | 0 |
| 5. Provenance | Provenance Detector | ? ? | 0 | 0 | 0 |

## Overall: **DRIFT**

> **REVIEW BEFORE SUBMISSION.** Drift findings require user adjudication.

## 1. Raw Data. Raw Data Integrity
**Status:** [!!] DRIFT
**Summary:** PASS on data-layer structural integrity (parsing, eid uniqueness, biomarker ranges, recruitment dates, death sanity). DRIFT on code-list completeness: G45 (TIA) missing from pipeline ASCVD code-list, and the downloaded UKB MACE-date extract contains hypertension first-occurrence fields (p131286-p131294) rather than the correct ASCVD fields (p131298 I21, p131306 I25). The single correct ASCVD-adjacent field present is p131296 (I20 unstable angina, 7.84% prevalence — plausible). Recommend Agent 3 (Feature Engineering) verifies which field IDs the TUDOR pipeline actually maps to ASCVD events, and that any I35 reference in code is NOT in an ASCVD composite. Severe-AS profile not applicable to this FH-lipid paper.


## 2. Cohort. Cohort Definition Validator
**Status:** ? ?
**Summary:** (none)


## 3. Features. Feature Engineering Auditor
**Status:** ? ?
**Summary:** Engineered features arithmetically reproduce cleanly within tolerance for all spot-checks (Trig_Filter, ldl_ut, apob_ldl, nhdl_ldl_gap recomputed 100/100 PASS within 0.001). However, the LDL treatment-adjustment implementation deviates substantially from the manuscript spec.


## 4. Statistics. Statistical Reproducer
**Status:** ? ?
**Summary:** (none)


## 5. Provenance. Provenance Detector
**Status:** ? ?
**Summary:** {'total_claims_extracted': 190, 'n_traced': 0, 'n_hardcoded_drift': 10, 'n_orphan': 0, 'n_label_drift': 2, 'decision_rationale': 'All ten primary metrics in the manuscript prose match hardcoded string literals in TUDOR_write_manuscript.R and TUDOR_marked_manuscript.R that DO NOT match the current live pipeline values cached in tudor_loco_output/*.csv. The TUDOR_LIVE_LEDGER.csv project ledger itself flags every one of these as DRIFT. This is the canonical TUDOR-style failure mode: numbers correctly computed at one historical pipeline state, hand-transferred into the manuscript-generator script as cat()/rd() string literals, never refreshed when the pipeline evolved.'}

- **NRI of TUDOR over DLCN was 0.358 (Wales)** -- 
  - Implication: Headline reclassification claim does not reproduce. The literal 0.358 has no live source in any current CSV; the closest match (0.3118) belongs to SouthWales cohort, not Wales. Sign reversal between literal (+0.358) and live (-0.0969) means the conclusion 'more than one-third reclassified correctly' is contradicted by current pipeline.
  - Proposed fix: Replace hardcoded literal in both R generators with live read from tudor_loco_output/nri_idi_results.csv. Decide which cohort label is intended and rerun NRI on that cohort; update prose accordingly.
- **IDI = +0.039 (Wales)** -- 
  - Implication: An order-of-magnitude overstatement of integrated discrimination improvement. The literal 0.039 sits between Wales (0.0024) and SouthWales (0.124) but matches neither.
  - Proposed fix: Read IDI live from CSV; pick the cohort intended by the manuscript and refresh.
- **Calibration slope was 6.33 (95% CI 5.93–6.73) in UK Biobank** -- 
  - Implication: Headline calibration claim of 6.33 (rationalised in prose as 'expected consequence of training under 33.2% prevalence and validating at 1.26%') matches no CSV. Live ledger value 1.23 is close to ideal (1.0), which would INVERT the manuscript's narrative — TUDOR is actually well-calibrated, not over-predicting. The 95% CI 5.93–6.73 is also fully hardcoded (no CSV trace).
  - Proposed fix: Decide one canonical UKB lipid-clinic subset definition; recompute calibration slope + CI; refresh prose and Figure 3 caption.
- **Brier score was 0.069 in UK Biobank** -- 
  - Implication: Literal 0.069 is ~50% higher than the live value (0.046). Internal CSV inconsistency (0.046 vs 0.013 vs 0.0085) indicates the manuscript value was computed against a different (now superseded) UKB subset definition.
  - Proposed fix: Read Brier live from a single canonical CSV; refresh literal.
- **DLCN scoring AUC 0.791 in Wales** -- 
  - Implication: Headline TUDOR-vs-DLCN comparison overstates the DLCN baseline by ~10 AUC points. The literal 0.791 came from an old computation on a different sub-population (likely cascade-relatives-with-DLCN-imputed-zero) and was never refreshed when the analysis was restricted to DLCN-scorable patients only.
  - Proposed fix: Decide canonical denominator for DLCN AUC (matched-subset is the methodologically defensible choice). Update prose and Figure 1 caption.
- **eDLCN AUC 0.636 in UK Biobank (with DeLong Z = 10.08, p = 6.73 × 10⁻²⁴ vs TUDOR)** -- 
  - Implication: Current pipeline gives eDLCN AUC 0.71 vs manuscript 0.64. Combined with the TUDOR drift below (0.75 → 0.7532), the manuscript's headline comparison 'TUDOR 0.750 vs eDLCN 0.636 (Δ=0.114)' shrinks under the live pipeline to 0.753 vs 0.713 (Δ=0.04). The DeLong Z=10.08 and p=6.73×10⁻²⁴ are therefore also stale literals — current Δ would not produce that magnitude of significance.
  - Proposed fix: Recompute DeLong on live cohorts; refresh all 6 abstract/results/caption mentions; the conclusion 'TUDOR substantially outperforms eDLCN in UKB' may need to be softened to 'modestly outperforms'.
- **TUDOR AUC 0.842 (95% CI 0.822–0.863) in Wales** -- 
  - Implication: Headline figure for the paper's flagship validation cohort is either an outdated literal (delta only -0.008) or a label-cohort mismatch (Wales-labelled-as-SouthWales, delta -0.05).
  - Proposed fix: Clarify which cohort is meant. Either: (a) use Wales = full registry and quote 0.78, or (b) keep 0.84 but relabel as 'South Wales PASS' throughout.
- **TUDOR AUC 0.750 (95% CI 0.731–0.770) in UK Biobank** -- 
  - Implication: Small drift, but the 95% CI 0.731–0.770 is also a hardcoded literal. Pipeline now gives lipid-clinic-only AUC closer to 0.77, which would actually STRENGTHEN the paper if the canonical subset were resolved.
  - Proposed fix: Fix the subset definition once; refresh AUC + CI from that single source.
- **LDLR AUC 0.839 (Wales)** -- 
  - Implication: Gene-specific Wales LDLR AUC stated as 0.839 but live = 0.765. A ~7-AUC-point drop weakens the 'consistent performance across LDLR, APOB, APOE' subclaim.
  - Proposed fix: Refresh gene-stratified table from current loco_predictions_complete.csv.
- **APOB AUC 0.841 (Wales)** -- 
  - Implication: Gene-specific Wales APOB AUC stated as 0.841 but live = 0.752. Largest gene-specific drift (~9 AUC points). Combined with F9, the manuscript's 'consistent across genes' framing is weakened.
  - Proposed fix: Refresh from current gene-stratified pipeline output.


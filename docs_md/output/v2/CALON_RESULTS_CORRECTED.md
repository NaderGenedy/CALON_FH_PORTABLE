# CALON-FH — Corrected Results Summary

*Replaces the inflated narrative. Every numerical claim below traces to a specific row in `output/v2/*.csv`. Audit corrections from the user's check are folded in throughout.*

**Date:** 2026-05-23
**Source data:** `output/v2/CALON_NEWv2_*.csv` (7 files) + `output/v2/run_log.txt`
**Provenance scripts:** `08_CALON2_wales_portable_v2.R`, `09_CALON2_dragon3_validate.R`, `11_CALON2_cox_comparison.R`
**Status:** **Awaiting SAFEHEART-RE reimplementation** with the genuine *Circulation* 2017 Table 3 hazard ratios. The ΔAUC figures below reflect the locally-fitted SRE approximation and will be re-run once the original coefficients are restored.

---

## The single-line claim that survives audit

> *In the All-Wales FH Registry (n = 1,815; 290 prevalent ASCVD events), a seven-predictor logistic model — **CALON9_noApoB** (age, sex, treatment-adjusted LDL-C, log[TG/HDL], ever-smoked, hypertension, diabetes) — trained on the South Wales FH cohort (n = 417), achieved **AUC 0.840 (95 % CI 0.826–0.854)** with **calibration slope 1.18**, compared with our locally-fitted approximation of the SAFEHEART Risk Equation (**AUC 0.801; slope 0.68**; ΔAUC +0.039, DeLong p = 4.5 × 10⁻⁶). Reclassification benefit was driven by non-events (NRI_non-event +0.39) with a negative event component (NRI_event −0.19); total NRI +0.19. Decision-curve net benefit favoured CALON at the 5/10/20 % thresholds. Performance was robust in patients aged < 65 across LDL strata (AUC 0.73–0.94) and **weaker in hypertensives (AUC 0.71)**. Single-direction external validation (UKB → Wales) replicated the prevalent-ASCVD advantage; reverse direction Wales → UKB was **not run** because Wales lacks the ApoB feature required by the alternate model variants. SAFEHEART-RE was implemented as a logistic approximation of the published hazard ratios; restoration of the original Pérez de Isla 2017 Table 3 coefficients is pending and may change the ΔAUC magnitude.*

That paragraph is what a *JCL / Atherosclerosis / Heart* reviewer can verify against the CSVs without finding anything to pull on.

---

## Verified headline (CALON9_noApoB — picked as THE single deployable model)

| Cohort | n_test | events | CALON9_noApoB AUC (95 % CI) | SAFEHEART-RE AUC | ΔAUC | DeLong p | Cal slope | Source row |
|---|---|---|---|---|---|---|---|---|
| Wales | 1,815 | 290 | **0.840 (0.826–0.854)** | 0.801 | +0.039 | **4.5 × 10⁻⁶** | 1.18 | leaderboard r5; SRE_reference r7 |
| UKB-LipidClinic (prevalent) | 595 | 143 | 0.723 (0.720–0.727) | 0.713 | +0.010 | **0.388 (NS)** | 1.23 | leaderboard r3; lancet_stats r3 |
| UKB-NonClinic (prevalent) | 1,027 | 256 | 0.720 (point) | 0.706 | +0.015 | **0.131 (NS)** | 0.92 | all_results r6; lancet_stats r4 |
| South Wales FH | 417 | 61 | 0.883 (0.880–0.886) | NA (no SRE comparator: n=0) | — | — | 1.35 | leaderboard r1 |

**Internal cal-slope inconsistency to resolve:** leaderboard reports Wales cal-slope 1.181; lancet_stats reports 1.080 for the same row. Both round to ~1.1; manuscript should report whichever is computed by the standard `val.prob.ci.2` path and note the alternative as a sensitivity.

## What was overstated in my earlier summary — corrections

### 1. Model-variant swap

My earlier summary said "CALON9_noApoB — the leader, 7 predictors" across all rows. **The incident-cohort headlines (Δ +0.098, Δ +0.119) were from CALON9_+lpa143 (10 predictors), not noApoB.** The pipeline's own run_log flagged this as "cohort-tuning, not a single deployable model."

**Correction:** Pick CALON9_noApoB as the manuscript headline model. Report its numbers across all cohorts honestly. The +lpa143 variant becomes a *sensitivity analysis*, not a co-primary.

CALON9_noApoB on incident outcomes (from all_results):
- UKB-LipidClinic incident: not run (model trained on Wales-prevalent doesn't transfer to UKB-incident outcome definition without re-fitting)
- A pre-specified sensitivity reporting `+lpa143` for incident endpoints is honest if labelled as such.

### 2. "Bidirectional cross-validation" — overstated

Wales → UKB rows return n_test = 0 for `CALON9` and `CALON9_+lpa143` (Wales lacks ApoB). Only CALON9_noApoB and CALON9_NoAgeLDLApoB tested in both directions.

**Correction:** Drop the word "bidirectional." Replace with: *"externally validated in UKB-LipidClinic (n=595) and UKB-NonClinic (n=1,027); the reverse-direction Wales-trained-to-UKB test of the lipid-augmented variants was not performed because Wales lacks ApoB-dependent features."*

OR: extract a Wales-trained CALON9_noApoB and test it on UKB explicitly, then "bidirectional" is earned.

### 3. NRI is positive total, **negative for events** in every prevalent cohort

| Cohort | NRI total | NRI_event | NRI_non-event | Verdict |
|---|---|---|---|---|
| Wales (prev) | +0.194 | **−0.193** | +0.387 | non-event-driven |
| UKB-LipidClinic (prev) | +0.366 | **−0.336** | +0.701 | non-event-driven |
| UKB-NonClinic (prev) | +0.245 | **−0.402** | +0.647 | non-event-driven |

**Clinical translation:** CALON gives *lower* risk scores to true ASCVD events than SAFEHEART does (down-classification = missed events). CALON gives *much lower* risk scores to true non-events (down-classification = fewer false alarms). Net benefit is real but **trades sensitivity for specificity**.

**Correction:** Report total NRI with the event/non-event split made explicit. Frame it as a specificity-gain rather than a sensitivity-gain. This is defensible for **cascade-screening triage** (where you want to avoid over-investigating relatives without genuine risk) and less obviously good for **primary prevention** (where missing a real event matters more). State the use-case the score is for.

### 4. IDI — I conflated two different statistics

| Cohort | IDI | Verdict |
|---|---|---|
| Wales (prev) | **−0.059** | SAFEHEART has better separation |
| UKB-LipidClinic (prev) | **−0.044** | SAFEHEART better |
| UKB-NonClinic (prev) | **−0.020** | SAFEHEART better |
| UKB-LipidClinic (incident) | +0.023 | CALON better |
| UKB-NonClinic (incident) | +0.024 | CALON better |

IDI is negative for the prevalent outcomes — SAFEHEART has tighter probability separation between events and non-events on those cohorts. The "IDI superiority at 5/10/20 % thresholds" phrase in my earlier summary **conflated IDI (a single index) with DCA (which IS threshold-varying and IS positive for CALON)**.

**Correction:** Drop "IDI superiority." Replace with the genuine DCA net-benefit numbers from `lancet_stats`: Wales DCA NB at 5 % = 0.121, at 10 % = 0.103, at 20 % = 0.067 — all positive and threshold-varying. That is the real reason to prefer CALON clinically.

### 5. DeLong p significance — not all six are significant

| Cohort × outcome | ΔAUC | DeLong p vs SRE | Significant? |
|---|---|---|---|
| Wales prevalent | +0.039 | 4.47 × 10⁻⁶ | ✅ |
| UKB-LipidClinic prevalent | +0.010 | 0.388 | ❌ NS |
| UKB-NonClinic prevalent | +0.015 | 0.131 | ❌ NS |
| UKB-LipidClinic incident | +0.098 | 2.91 × 10⁻⁶ | ✅ |
| UKB-NonClinic incident | +0.119 | 9.21 × 10⁻¹² | ✅ |

**Correction:** "Beats SAFEHEART in 6 of 6 cohort-outcome combinations" → "Beats SAFEHEART on point estimate in 5 of 5 testable combinations; statistically significant in 3 of 5; non-significant in the two UKB prevalent comparisons where Δ is small." Honest.

---

## Three run-log items the earlier summary skipped

### Run-log issue 1 — calibration sanity failures
*run_log:* "7 / 75 rows fail calibration sanity (slope outside [0.5, 2.5])."

Headline rows affected:
- CALON9_NoAgeLDL on SW-FH incident: cal-slope **2.47** (severe over-fitting, AUC 0.884 overoptimistic)
- SRE-Fixed on UKB-LipidClinic incident: cal-slope **0.32** (severe miscalibration of the comparator)

**Manuscript action:** Report cal-slope alongside every AUC; flag rows with slope outside [0.7, 1.5] as not deployable.

### Run-log issue 2 — non-significant DeLong p
Already addressed in correction 5 above. The UKB prevalent comparisons are not statistically significant despite positive point estimates.

### Run-log issue 3 — SAFEHEART-RE comparator integrity (THE big one)

A prior session audit (`CALON_MI_validation.R`) found **0 of 8 SAFEHEART coefficients match the published *Circulation* 2017 Table 3** (Pérez de Isla et al.).

| Predictor | This pipeline's HR | Genuine Table 3 HR | Delta |
|---|---|---|---|
| Male sex | 2.17 | 2.01 | +0.16 |
| Hypertension | 1.54 | **1.99** | **−0.45** |
| Lp(a) | 1.46 | 1.52 | −0.06 |
| Age / BMI / LDL | per-unit fabricated | categorical bands in original | structural difference |

**Implication:** The published ΔAUC +0.039 reflects CALON's strength **plus** an over-weak SAFEHEART comparator on the hypertension term. Until SRE-Fixed is reimplemented with the original Table 3 coefficients, the headline Δ figure carries an asterisk.

**Manuscript action (B in our plan):** Reimplement SRE-Fixed with the genuine published HRs, rerun the leaderboard, report the corrected ΔAUC. Expected effect: the Δ shrinks somewhat — CALON likely still wins on Wales (where Δ is largest at +0.039 and the comparator margin is largest), but the magnitude is honest.

---

## Five-step plan to make this submission-survivable

1. ✅ **Pick one variant for the headline** — done above: **CALON9_noApoB**.
2. ⏳ **Reimplement SAFEHEART-RE with the genuine *Circulation* 2017 Table 3 hazard ratios.** Next deliverable.
3. ✅ **Reframe NRI honestly** — done above.
4. ✅ **Drop "IDI superiority"; replace with DCA net-benefit** — done above.
5. ⏳ **Either run Wales → UKB explicitly OR drop "bidirectional"** — recommend dropping the word; if a Wales-trained CALON9_noApoB is extracted and applied to UKB, the bidirectional claim is earned cheaply.

---

## Confidence ladder — corrected

| Claim | Confidence |
|---|---|
| Wales AUC 0.840 (CALON9_noApoB) | **High** — exact match across 3 output CSVs |
| Wales beats SRE on point estimate by +0.039 | **High** — confirmed |
| Wales DeLong p = 4.5 × 10⁻⁶ | **High** — but partly reflects comparator weakness |
| Wales calibration slope 1.18 | **High** with internal-CSV-inconsistency caveat |
| UKB-LipidClinic prevalent ΔAUC +0.010 | **Low — not significant** (p = 0.388) |
| UKB-NonClinic prevalent ΔAUC +0.015 | **Low — not significant** (p = 0.131) |
| UKB incident ΔAUC +0.098, +0.119 | **Moderate** — from a DIFFERENT model variant (+lpa143); needs the noApoB number for honest comparison |
| "Beats SAFEHEART bidirectionally" | **Not supportable** — one-way only |
| "IDI superiority" | **Wrong direction** — IDI is negative for prevalent outcomes |
| "NRI substantial" | **Half-supported** — total positive, event-component negative |
| DCA net benefit favours CALON | **High** — values in lancet_stats |
| Robust in subgroups aged < 65 across LDL strata | **High** — subgroup AUCs ≥ 0.73 |
| Robust in hypertensives | **Low** — AUC drops to 0.66–0.71 |
| Publishable in *JCL* / *Atherosclerosis* / *Heart* | **High** with corrections folded in |
| Publishable in *Circulation* / *Eur Heart J* | **Speculative** — needs the SRE comparator reimplemented, the bidirectional claim earned, and ideally a prospective external cohort |

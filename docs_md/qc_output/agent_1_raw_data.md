# Agent 1 — Raw UKB Data Integrity Report

**Project**: TUDOR (FH lipid-clinic diagnostic algorithm)
**Date**: 2026-05-12
**Status**: DRIFT (data layer clean; pipeline-side code-list incomplete)

## Traffic light summary

| Check | Status | Note |
|---|---|---|
| File parseability | GREEN | All 6 sampled files parse cleanly |
| EID uniqueness | GREEN | 501,936 unique eids in every file |
| Biomarker ranges | GREEN | 10 / ~1.85M outside-range (<0.001%) |
| Recruitment-date window | GREEN | All 2006-03-13 to 2010-10-01 |
| Death-before-recruitment | GREEN | 0 violations |
| UKB sentinel dates | AMBER | 3 rows of 1902-02-02; must filter |
| ASCVD ICD-10 completeness | AMBER | G45 (TIA) missing from pipeline |
| ASCVD first-occurrence fields | RED | Extract contains p131286-94 (hypertension I10-I15), missing the correct ASCVD fields p131298 (I21) and p131306 (I25) |
| Forbidden code I35 in code | AMBER | I35 string appears in pipeline — needs Agent 3 disambiguation |
| Severe-AS profile (K611 etc.) | N/A | TUDOR is a lipid paper, not a valve paper |

## Files sampled

| File | Rows | Cols | Eid unique | Notes |
|---|---|---|---|---|
| ukb_lpa.csv | 501,936 | 3 | YES | Lp(a) baseline 75% coverage; instance-1 3% |
| ukb_nmr_batch1.csv | 501,936 | 7 | YES | NMR core metabolites |
| ukb_dates_mace.csv | 501,936 | 10 | YES | Contains p131286-p131296 + death date |
| ukb_recruitment_dates.csv | 501,936 | 5 | YES | Age 37-73, median 58 |
| ukb_reviewer_longitudinal_lipids.csv | 501,936 | 9 | YES | LDL/TC/HDL/TG i0 + i1 |
| ukb_reviewer_demographics.csv | 501,936 | 10 | YES | BMI / ethnicity / alcohol |

All files use UKB RAP export naming: `participant.eid`, `participant.p<FIELD>_i<INSTANCE>`. Analysis scripts must alias or strip the `participant.` prefix.

## Biomarker range checks

| Field | Biology | Range | Min | Max | Median | Non-null | Outside | %Outside |
|---|---|---|---|---|---|---|---|---|
| p30780_i0 | LDL baseline | 0.5–15 mmol/L | 0.266 | 9.797 | 3.516 | 468,180 | 2 | 0.0004% |
| p30690_i0 | TC baseline | 1.5–20 mmol/L | 0.601 | 15.46 | 5.65 | 469,062 | 3 | 0.0006% |
| p30760_i0 | HDL baseline | 0.3–5 mmol/L | 0.219 | 4.401 | 1.398 | 429,393 | 5 | 0.001% |
| p30870_i0 | TG baseline | 0.2–30 mmol/L | 0.231 | 11.278 | 1.483 | 468,688 | 0 | 0% |
| p30790_i0 | Lp(a) baseline | 0–400 nmol/L | 3.8 | 189.0 | 21.1 | 375,200 | 0 | 0% |
| p21022 | Age at recruit | 37–73 | 37 | 73 | 58 | 501,936 | 0 | 0% |
| p21001_i0 | BMI baseline | 12–75 | 12.12 | 74.68 | 26.75 | 498,833 | 0 | 0% |

10 outside-range values across ~1.85M measurements — assay-floor / extreme-tail artefact, not calibration drift.

## Code-list audit (ASCVD composite profile)

Pipeline scripts in `C:/Users/nader/Downloads/calon_ukb_pipeline/` searched for ICD-10 string literals.

**Required (FOUND, 10/11)**: I20, I21, I22, I23, I24, I25, I63, I70, I73, I74
**MISSING**: **G45** (TIA) — modest case-count impact but TRIPOD-incomplete if ASCVD is used as TUDOR secondary outcome.

**Forbidden but FOUND**: **I35** (aortic stenosis). Agent 3 must verify whether I35 is referenced inside an ASCVD code list (would inflate events with non-atherosclerotic disease) or only as an unrelated valve-disease comparator.

**First-occurrence ASCVD fields** (per QC reference: should be `p131296` I20, `p131298` I21, `p131306` I25):
- Present in `ukb_dates_mace.csv`: `p131286`, `p131288`, `p131290`, `p131292`, `p131294`, `p131296`
- Of these, only `p131296` is genuinely an ASCVD field (unstable angina I20).
- `p131286–p131294` are flagged in CLAUDE.md as coding I10–I15 (hypertension family), NOT ASCVD.
- Prevalence evidence: `p131286` non-null = 209,310 (41.7%) — matches UKB hypertension prevalence; clearly not MI/ASCVD.
- The correct MI field `p131298` (I21) and chronic-IHD field `p131306` (I25) are **NOT downloaded**.

**Implication for TUDOR**: if any analysis script uses `p131286` (or its siblings) as an ASCVD event date, the event count will be inflated ~50× and the diagnostic algorithm performance will be invalid. Agent 3 must verify the mapping.

## Date sanity

- Recruitment window (p53_i0): all values 2006-03-13 to 2010-10-01 — PASS
- Death dates (p40000_i0): 2006-05-10 to 2024-12-02; n=56,961; zero deaths before recruitment — PASS
- UKB sentinel `1902-02-02`: 2 rows in p131286, 1 row in p131296 — trivial but downstream code must filter explicitly
- No future dates beyond today (2026-05-12) — PASS

## Severe aortic stenosis profile

**NOT APPLICABLE** to TUDOR (FH lipid-clinic algorithm, no valve-procedure component). K611 / K261–K263 audit not performed and not needed. Recorded here only because the QC five-fold harness has a `severe_AS_intervention` profile and the user asked it be acknowledged as not applicable.

## Decision

**DRIFT** — the raw data layer is structurally clean (parse, eid, ranges, dates all pass). The drift is at the **field-to-endpoint mapping** layer: the downloaded MACE-date file contains the wrong UKB first-occurrence fields for ASCVD, and the pipeline code-list misses G45 while ambiguously referencing I35. These are pipeline-design issues (Agent 3 domain) that surface here because they affect interpretation of the raw extract.

## Handoffs to other agents

- **Agent 3 (Feature Engineering)**: verify which `p131xxx` fields the TUDOR scripts map to ASCVD events; confirm I35 is not inside an ASCVD code list; consider re-extracting p131298 and p131306 from UKB RAP.
- **Agent 4 (Statistical Reproducer)**: if the wrong first-occurrence fields were used, all secondary-outcome HRs / event counts must be recomputed.
- **Agent 5 (Manuscript)**: verify TUDOR text states the exact ICD-10 codes and UKB field IDs used for any ASCVD secondary outcome.

# CALON-FH Mendelian Randomization — Design Specification

**Date:** 2026-05-21
**Author:** Dr Nader Genedy, University Hospital of Wales, Cardiff
**Project:** CALON-FH clinical risk model
**Status:** approved by user 2026-05-21; ready for implementation plan

---

## 1. Purpose

A Mendelian randomization (MR) analysis serving as a **confounding-resolution
section** within the CALON-FH clinical-model paper.

The cross-sectional clinical model (`CALON_clinical_model.py`) returned a
counter-intuitive odds ratio for LDL-C of **0.88 per SD (< 1)**. This is a
known cross-sectional artefact: in a prevalent-ASCVD design, participants with
ASCVD are heavily statin-treated, so their *measured* LDL-C is low, making LDL
appear spuriously protective. MR — which instruments LDL-C with germline
genetic variants and is therefore immune to reverse causation and treatment
confounding — is used to demonstrate the **true causal direction**: genetically
instrumented LDL-C raises incident ASCVD risk (HR > 1).

The deliverable is a manuscript section plus its reproducible code, **not** a
standalone MR paper. LDL→ASCVD causality is already among the most replicated
findings in genetic epidemiology (Ference et al. 2017); this analysis is
explicitly **confirmatory** and is framed as such.

## 2. Scope decisions (locked)

| Decision | Choice | Rationale |
|---|---|---|
| Purpose | Confounding-resolution section in the clinical paper | Does real scientific work without over-claiming |
| Instruments | **Both** polygenic (LDL-PRS) + monogenic (LDLR null) | Triangulation across non-shared failure modes |
| Ancestry | Self-reported White British; **no** PC adjustment | Only self-reported ethnicity on disk; no genetic PCs |
| Sensitivity | Core + one negative-control outcome | Cheapest defence of the no-PC-adjustment limitation |

## 3. Data sources (confirmed on disk)

| Element | Field / file |
|---|---|
| Polygenic instrument | `ldl_prs` — `ukb_FULL_MASTER.csv` |
| Exposure | `ldl_chem` (measured LDL-C, mmol/L) — `ukb_FULL_MASTER.csv` |
| Monogenic instrument | `ukb_carriers_FINAL.csv` — LDLR carriers; `consequence` column for null-allele subset |
| Outcome (incident ASCVD) | `first_ascvd`, `t_event_years` — `ukb_FULL_MASTER.csv` |
| Prevalent-ASCVD flag | derived from HES ICD-10 array (`ukb_icd10_full.csv`, p41270) — correct ASCVD codes |
| Ancestry restriction | `ethnicity_code` — `ukb_FULL_MASTER.csv` |
| Covariates | age, sex — `ukb_FULL_MASTER.csv` / `calon_sex.csv` |
| Negative-control outcome | non-atherosclerotic ICD-10 (appendicitis K35 or fracture S–T) — `ukb_icd10_full.csv` |

The `first_ascvd` / `t_event_years` fields are derived from the correctly-coded
first-occurrence components (acute MI, stroke, other IHD). The misnamed
`ukb_dates_mace.csv` (p131286–p131296 = hypertension I10–I15) is **not used**.

## 4. Cohort definition

1. Base: UKB participants with a complete baseline lipid profile (the
   ~428,000-participant clinical-model cohort).
2. Restrict: self-reported White British (`ethnicity_code`).
3. Exclude: prevalent ASCVD at baseline — an incident analysis cannot count an
   event in a participant who already has ASCVD.
4. Follow-up time: `t_event_years`; event: incident ASCVD (`first_ascvd`).

## 5. Statistical method — Wald-ratio MR

For both instruments, the causal estimate is the Wald ratio:

```
causal log-HR per mmol/L LDL-C  =  beta_ZY  /  beta_ZX
```

### Instrument 1 — polygenic (LDL-PRS)
- `beta_ZX` = coefficient of `ldl_chem ~ ldl_prs + age + sex` (OLS) — mmol/L LDL per unit PRS.
- `beta_ZY` = log-HR from `Cox(incident ASCVD) ~ ldl_prs + age + sex`.
- Causal HR = `exp(beta_ZY / beta_ZX)` per 1 mmol/L genetically-predicted LDL-C.
- **Instrument strength:** report the stage-1 F-statistic for `ldl_prs`; require F ≫ 10.
- **95% CI:** non-parametric bootstrap, 1,000 resamples (primary — handles the
  ratio-of-coefficients and the Cox numerator correctly). The delta-method
  CI is computed as a secondary cross-check only.

### Instrument 2 — monogenic (LDLR null alleles)
- Null-allele subset = `consequence` in {stop_gained, frameshift_variant,
  splice_acceptor_variant, splice_donor_variant, canonical splice}. Missense
  carriers are **excluded** — heterogeneous molecular effect would dilute the
  instrument. Approximately 459 carriers expected.
- `beta_ZX` = mean LDL difference (null-carrier − non-carrier), age/sex-adjusted.
- `beta_ZY` = log-HR from `Cox(incident ASCVD) ~ null_carrier + age + sex`.
- Causal HR = Wald ratio, scaled per mmol/L. Wide CI expected (small n) — this
  is accepted; the monogenic arm is the orthogonal confirmation, not the
  precision estimate.

### Triangulation
The two instruments have **non-shared failure modes** — the PRS risks
pleiotropy; the monogenic null alleles risk small-n weak-instrument bias and
selection. Concordant HR > 1 from both is genuine triangulation.

## 6. Negative-control falsification

Test `LDL-PRS → negative-control outcome`. The negative control is
**appendicitis (ICD-10 K35)** — an acute surgical condition with no plausible
LDL-cholesterol pathway, defined by a single clean code. If incident-K35 events
are too few for a stable Cox fit (< ~200), the pipeline falls back to fracture
(ICD-10 S–T codes) and records the substitution. Method:
`Cox(negative outcome) ~ ldl_prs + age + sex`.

- **Expected (pass):** null association (HR ≈ 1, CI spanning 1).
- **Interpretation:** a null here shows the LDL-PRS→ASCVD signal is
  LDL-specific, not population-structure masquerading as causation — a direct
  defence of the absent PC adjustment.
- A non-null negative control would flag residual stratification and would be
  reported honestly as a caveat.

## 7. Outputs

| File | Contents |
|---|---|
| `CALON_MR.py` | the runnable pipeline (single script, phased) |
| `CALON_MR_results.csv` | causal HRs + 95% CI for both instruments; stage-1 F-statistics; negative-control HR |
| `CALON_MR_TRACEABILITY.csv` | per-quantity provenance — value → raw file → field → formula |
| `CALON_MR_results.md` | the confounding-resolution narrative: clinical OR 0.88 vs MR HR > 1, triangulation, honest limitations |

## 8. Reporting requirements

- Headline comparison table: clinical-model LDL OR 0.88 (cross-sectional) vs MR
  causal HR per mmol/L (polygenic) vs MR causal HR (monogenic).
- The narrative claim is strictly **directional**: LDL-C is causally
  atherogenic; the cross-sectional OR < 1 was treatment confounding. No precise
  effect-size claim is made.
- TRIPOD/STROBE-MR-style transparency: instrument, assumptions, F-statistic,
  cohort flow.

## 9. Limitations (pre-committed to the write-up)

1. **One-sample MR** — exposure and outcome from the same cohort; weak-instrument
   bias is minimal with a strong PRS but the design is less robust than
   two-sample MR.
2. **Single composite PRS** — precludes MR-Egger and weighted-median pleiotropy
   sensitivity analyses (these require per-variant data). The monogenic
   instrument partially mitigates.
3. **No genetic principal-component adjustment** — genetic PCs are not on disk;
   only self-reported ethnicity. Residual population stratification is possible
   but is not expected to be directional enough to invert HR > 1. The
   negative-control outcome tests this.
4. **No relatedness exclusion** — kinship data not on disk; related pairs cannot
   be removed. Minor precision inflation; negligible point-estimate bias.

The conclusion claimed (LDL is causal; the OR < 1 was confounding) is robust to
all four limitations. A precise effect-size estimate is **not** claimed.

## 10. Explicitly out of scope (YAGNI)

- Per-component MI vs ischaemic-stroke MR breakdown.
- Multivariable lipid MR (no HDL or TG PRS on disk).
- Standalone-paper framing or Nobel-tier causal claims.
- Two-sample MR using external GWAS summary statistics.

## 11. Success criteria

1. `CALON_MR.py` runs end-to-end from raw CSVs and writes all four outputs.
2. Stage-1 F-statistic for the LDL-PRS instrument is ≫ 10 (strong instrument).
3. Both instruments return a causal HR point estimate > 1 for incident ASCVD.
4. The negative-control outcome returns a null association (CI spans 1).
5. Every number in `CALON_MR_results.csv` traces to a raw file and field in
   `CALON_MR_TRACEABILITY.csv`.

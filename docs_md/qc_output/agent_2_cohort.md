# Agent 2 — Cohort Definition Validator: TUDOR Audit

**Decision: FAIL**

Two large cohort-count drifts and an absent family-level deduplication step, plus a recoverable label-ambiguity issue between "Wales" prose and "SouthWales/Wales" code.

## Cohort counts: manuscript vs live pipeline

| Cohort | Manuscript claim | Cached predictions | Rebuilt from raw | Verdict |
|---|---|---|---|---|
| Wales All-PASS | 7,253 / 2,405 (33.2%) | 5,376 / 1,862 (34.6%) | 7,253 / 2,405 raw, drops to 5,376 after complete-cases on features | DRIFT — manuscript reports raw N, live analysis uses smaller filtered n |
| SouthWales / CAVUHB development | (not named) | 1,072 / 336 (31.3%) | 1,366 raw, 1,072 after filter | Hidden by ambiguous "Wales" label in prose |
| UKB lipid clinic | 58,021 / 729 (1.26%) | 107,090 / 938 (0.88%) | 60,028 / 777 (1.29%) when strict §2.2 criteria applied | CRITICAL DRIFT — 1.85× n discrepancy |
| Total genetically confirmed FH | 4,028 across 6 genes | 3,136 in predictions | — | Off by 892 |

## Why the UKB cohort is wrong

Manuscript §2.2 prose states:
> TC > 7.5 mmol/L OR statin-corrected LDL > 4.9 mmol/L OR premature ASCVD before 55M/60F

`TUDOR_LANCET_COMPLETE.R` lines 170-174 actually computes:
```
ukb_raw$lc <- (tc > 7.5) | (ldl_ut > 4.9) | (non_hdl > 5.9) | (on_statin == 1)
```

The `on_statin == 1` clause alone admits ~50,000 statin users regardless of any lipid threshold, doubling the cohort. Independent rebuild:

- Lancet-script filter: **n = 107,090, FH+ = 938, prev 0.88%** — matches cached predictions exactly
- Manuscript §2.2 strict (3-criterion, sex-stratified ASCVD): **n = 60,028, FH+ = 777, prev 1.29%** — within 3.5% of the manuscript's 58,021/729/1.26%

Conclusion: the manuscript headline 58,021/729/1.26% reflects the §2.2 prose criteria, but the **entire downstream pipeline** (LOCO-CV, AUC, NRI, calibration, mortality gradients) used the broader 107,090-row dataset. Every UKB statistic in the manuscript is computed on a cohort almost twice the stated size and with materially lower FH prevalence.

## Family-level deduplication — ABSENT

Grep across all `TUDOR_*.R` scripts: **zero matches for FamilyNumber, Pedigree, family_id**.

The raw `PASS_wrong_dob.sav` contains a `FamilyNumber` column (confirmed via `pyreadstat`). It is loaded but never used. Cascade-screened relatives — who share variants and produce correlated predictions by construction — are treated as independent observations in:
- AUC and DeLong tests
- NRI / IDI / Brier
- Calibration slope/intercept
- Subgroup analyses

Per user CLAUDE.md this is a **non-negotiable** requirement for any Wales-cohort FH analysis. Expected effect of fixing: Wales AUC will drop because intra-family correlation no longer inflates effective sample size.

## Cohort label cross-check

| Source | Cohort tokens used |
|---|---|
| Cached predictions CSV `cohort` column | `SouthWales`, `Wales`, `UKB` |
| `TUDOR_LANCET_COMPLETE.R` | `sw` (= "SouthWales"), `wales` (= PASS All-Wales) |
| Manuscript prose | "CAVUHB" (development), "Wales" / "All Wales FH Registry" / "PASS" (validation), "UK Biobank" |

The TUDOR scripts cleanly separate `SouthWales` (CAVUHB development cohort, n=1,072) from `Wales` (All-Wales PASS validation, n=5,376). The manuscript prose, however, says "TUDOR was developed at CAVUHB" and later "first external validation in the All Wales FH Registry" without ever telling the reader that the development cohort is itself a subset of Welsh patients (the same patient pool, in CAVUHB-only form, before being merged into PASS).

This is *not* a code bug, but it IS a reader-facing ambiguity. Per user CLAUDE.md ("South Wales / Dragon-3 is a subset of All-Wales PASS by FamilyNumber"), there is a real risk that some CAVUHB development patients appear ALSO in the All-Wales PASS validation set — and the absence of family-level dedup means this overlap is not removed.

## FH ascertainment — CORRECT

| Cohort | Definition | Script reference |
|---|---|---|
| Wales/SouthWales | `Positive1 == "Yes"` (genetic test) | TUDOR_LANCET_COMPLETE.R line 131 |
| UKB | `is_fh_genetic` (WES rare pathogenic variants in LDLR/APOB/PCSK9) | TUDOR_UKB_LIPID_CLINIC.R line 72 |

Both are gene-based; neither uses lipid-threshold-only ascertainment. CONSISTENT with the manuscript's claim of "genetically confirmed FH".

## Summary

| Check | Result |
|---|---|
| Cohort logic reproduced | FAIL (UKB filter in code is 4-criterion, prose says 3-criterion) |
| Family-level deduplication | FAIL (FamilyNumber present in raw, never used) |
| Cohort labels consistent | DRIFT (script SouthWales/Wales not reflected in prose) |
| FH ascertainment definition | PASS |
| Prevalence sanity (Wales 33.2%) | PASS for raw, off by +1.4% in predictions |
| Prevalence sanity (UKB 1.26%) | FAIL — live pipeline prevalence is 0.88% |

## Required actions before resubmission

1. **Pick one UKB filter and apply it consistently** — either rewrite §2.2 prose to match the 4-criterion code, or rerun the entire pipeline with the 3-criterion strict filter. Update n, FH+, prevalence, AUC, calibration, NRI, IDI, DCA, mortality gradient throughout.
2. **Add family-level deduplication** for the Wales analyses. Recompute Wales AUC (95% CI), DeLong vs DLCN, NRI, IDI on the dedup'd set.
3. **Disambiguate cohort labels** — first use of "Wales" in prose should say "Cardiff and Vale University Health Board (CAVUHB) FH clinic, n=1,072 — the South Wales / development cohort within the wider All-Wales PASS Registry". Second use should say "All-Wales PASS Registry (n=7,253), of which the CAVUHB cohort is a subset (removed from validation set)".
4. **Resolve total FH count** — manuscript 4,028, cached predictions 3,136. Likely the 4,028 includes patients dropped during complete-cases filtering; if so, prose should report the figure that matches the analytic sample.

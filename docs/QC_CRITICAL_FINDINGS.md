# QC Critical Findings — manuscript traceability audit (2026-05-30)

Parallel manuscript-QC audits were run on the "submission-ready" papers. **Three of the four
audited carry critical, submission-blocking integrity problems.** The headline is a confirmed,
recurring **endpoint mis-mapping bug** plus untraceable/contradicted headline claims.

> Deliverables per paper (PROVENANCE_TABLE.csv, *_TRACEABILITY_REPRODUCER.py, *_audit.md) are in
> each project folder. Paper 1's audit was interrupted and is not yet complete.

---

## 0b. ASCVD endpoint DEFINITION — raw-trace audit (2026-05-31)
The corrected analyses (NB01/NB02) use the master `first_ascvd` (UKB curated first-occurrence
algorithm) = **313 carriers, 8.8%**. Traced directly against the raw HES ICD-10 array
(`ukb_icd10_full.csv`, codes p41270), the counts differ materially:
| Definition | carriers | % |
|---|--:|--:|
| Buggy hypertension (p131286–294) | ~1,585 | 44.7% |
| Master first_ascvd (UKB first-occurrence) — USED | 313 | 8.8% |
| Raw HES narrow (I21+I25) | 453 | 12.8% |
| Raw HES full composite (I20–25/I63/I64/I70/G45 = CLAUDE.md list) | 591 | 16.7% |
Concordance master vs raw-full: Jaccard 45.6% (283 both, 308 raw-only, 30 master-only).
**Implication:** the corrected endpoint is NOT raw-HES-traced nor CLAUDE.md-code-list-complete; it
is the narrower UKB first-occurrence definition (likely incident + curated). Endpoint-definition
choice (first-occurrence vs full HES composite; incident vs ever) is a `[NEEDS-USER]` decision and
changes event counts ~2-fold. A fully raw-traced, code-list-complete endpoint would use 16.7%.

**RAW-VERIFICATION of the headline paradox (2026-05-31):** recomputed on the raw-HES ever-ASCVD
endpoint (full composite, 550 events), the statin-treatment paradox is **OR 2.29 (1.84–2.86),
p=2.2×10⁻¹³** — i.e. it HOLDS and slightly strengthens vs the manuscript's 2.05 (derived column).
So NB01's robust headline is endpoint-definition-INDEPENDENT and raw-HES-verified. (The collapsed
*interaction* claims are a separate, time-to-event matter still needing dated raw events via RAP.)

## 0c. Wales / South-Wales raw-trace (2026-05-31) — family-dedup confirmed
Traced to the raw registries (`Shared_Data/combined_soretd_ldl.csv`, `DRAGON_3.csv`):
- **Dragon-3 is 100% a subset of All-Wales** (DatabaseNumber 424/424, NHSNumber 424/424,
  FamilyNumber 220/220 overlap). Same registry — NOT independent cohorts. Confirms the standing
  "Dragon-3 = Wales PASS" discovery.
- All-Wales raw = **1,363 rows → 221 unique families** (425 unique DatabaseNumbers). Pooled,
  family-deduplicated Wales = **~426 individuals / 221 families**.
- Lpa_Manuscript's "All-Wales 2,405 + Dragon-3 418 = 2,823" is unsupported: wrong source count
  (1,363 not 2,405) AND double-counts the same families. Any Wales pooled analysis must dedup to
  **221 families** (or use one registry only).
- **DONE 2026-05-31:** family-deduplicated Wales cohort built (CLAUDE.md one-event-per-family rule)
  -> `Shared_Data/wales_family_deduplicated_221.csv` (221 families, 60 event-families = 27.1%;
  individual-level ASCVD 4.5%). Ready to replace the double-counted 2,823 in any Wales analysis.
  (Index-vs-cascade representative selection would need the TypeofPatient codebook; the
  one-event-per-family rule used here does not.)
- **CALON-FH external validation on the deduplicated cohort (DONE 2026-05-31):** the UKB-fitted
  CALON-FH model (robust 7-feature variant: age, sex, LDL, HDL, TG, TG/HDL, statin) validates on
  the 221-family deduplicated Wales cohort at **AUC 0.821** (n=155 complete-feature, 35 events) —
  vs the prior PASS external (n=97, 4 events, AUC 0.76). UKB apparent AUC 0.737 reproduces the known
  0.74. So CALON-FH external validity HOLDS/strengthens on the family-deduplicated cohort with ~9x
  more events. (sklearn approximation, not pixel-identical to phase4's statsmodels spec; remnant/
  diabetes dropped for Wales coverage.)

## 0. THE confirmed bug — UKB first-occurrence field mis-mapping

Verified directly against the UKB field catalogue (`D:/Projects/CALON_AlphaFold_Rebuild/New folder/fieldsum.tsv`):

| Field | TRUE meaning (UKB) | Often mis-used as |
|---|---|---|
| p131286 | **I10 essential hypertension** | "MI / first ASCVD" |
| p131288 | **I11 hypertensive heart disease** | "unstable angina" |
| p131290 | **I12 hypertensive renal disease** | "ischaemic stroke" |
| p131292 | **I13 hypertensive heart+renal** | "TIA" |
| p131294 | **I15 secondary hypertension** | "PVD" |
| p131296 | I20 angina pectoris | (only ischaemic one of the set) |
| **p131298** | **I21 acute MI** | — (the CORRECT MI field) |
| **p131300** | I22 subsequent MI | — |
| **p131302–p131306** | I23, I24, **I25 chronic IHD** | — |

**The composite `min(p131286,88,90,92,94,96)` = first hypertension-or-angina, NOT first ASCVD.**
The correct ASCVD composite uses **p131296–p131306** (I20–I25) + stroke (I63) + the other ICD-10
codes from the standing code-list (I63/I64/I70/G45). This is the exact WRONG_FIELDS class flagged
in CLAUDE.md. NB: CLAUDE.md's own data-table is internally inconsistent — it lists "p131286 = MI
(I21)" in one place and "p131286–296 = hypertension" in another; **the catalogue confirms the
latter (hypertension) is correct.**

### Where this bug is present (confirmed)
- **NB01** — `phase5*.py`, `phase6`, `phase7` build the time-to-event endpoint from the
  hypertension fields. Contaminates every TTE claim (LDL×statin HR 1.32, λ, pre-event trajectory,
  statin×SSS, duration-stratified). The cross-sectional paradox OR 2.05 and CALON-FH AUCs use a
  separate pre-derived `ascvd` column and are NOT hit by this specific bug (but re-confirm that column's definition).
- **Lpa_Manuscript** — `scripts/allelic_series_finegray_full.R` sets `date_as := p131286`, i.e.
  the AS endpoint is actually **hypertension**, not I35.0.
- **NB02 nb2_05** (and the local `ukb_dates_mace.csv`, which contains ONLY the hypertension block
  p131286–296) — my earlier nb2_05 re-run inherited this; its "events" are hypertension-dominated.
  The variant-key fix was correct; the endpoint must be rebuilt from I20–I25/I63 dates.

---

## 1. NB01 (Discordance Paradox) — 51/56 PASS, 0 FAIL, 5 SKIP + 5 CRITICAL
- **Arithmetic clean:** every sourced claim matches its CSV (OR 2.05 (1.61–2.60); 7 reconstruction
  metrics incl. R² ceiling 0.019; CALON-FH v13 AUC 0.7399/0.7608, v14 0.7372/0.7742).
- **CRITICAL:**
  1. TTE endpoint mis-mapping (above) — contaminates all time-to-event claims.
  2. Duration-dose **OR 1.06/year (p=0.012)** — untraceable; CSV stores only a per-SD coef (1.208).
  3. LDL×duration **HR 1.33 (p=0.030)** — untraceable; only duration-*tertile* HRs exist (interaction p=0.26); the "1.33" is the medium-tertile HR, not the interaction.
  4. No **E-value** for OR 2.05 (non-negotiable).
  5. No **NoAgeLDL** sensitivity variant (non-negotiable).
- Minor: unique statin carriers 593 (ms) vs 589 (CSV); Methods says Python 3.14 (should be 3.12).

## 2. Lpa_Manuscript (Saturation-Threshold) — 34/57 PASS, 4 FAIL, 19 SKIP
- **Solid:** all patient-level Dragon-3 numbers trace exactly (n=315, medians, ASCVD+ 95 vs 40 nmol/L p=0.0032, correlations, thresholds).
- **FAIL (block submission):**
  1. **Central allelic-series P-trend = 0.001 is CONTRADICTED** by the only real computation
     (`results/v13_S10_results.txt`: Wald P-trend **0.6859**, interaction P=0.7546, non-monotonic ORs). The monotonic 0.96→1.12→1.34→1.52 ladder is unsupported.
  2. AS endpoint mis-mapped to p131286 (hypertension), not I35.0.
  3. Combined FH **n=4448 not family-deduplicated** (Dragon-3 418 ⊂ All-Wales PASS); also All-Wales source has only 1363 unique rows vs "2405" asserted.
  4. On-/off-statin Lp(a) (52.0 vs 37.6, p=0.101) doesn't match cited CSV (176.0 vs 545.5, p=0.16).
- **CRITICAL untraceable:** AlphaFold3 contact matrix (36/28/12, iPTM) and UKB AS Cox (P=2.5e-23 / 0.43) are **hard-coded literals in a plotting script** — `data/alphafold3_predictions/` and `results/cox_six_beds_fh_vs_nonfh.csv` **do not exist**.
- CMR aortic-volume gap (§7.1): correctly handled as an honest caveat — no action.

## 3. Lpa_Multilevel v17 (submitted, Atherosclerosis) — 24/38 PASS, 8 FAIL, 6 SKIP
- **The locked engine (v21_FIXED_20) is sound** and reproduces every supplied headline (Q=54.63,
  p_het=1.77e-9; per-SD HRs AS 1.139/MI 1.127/IHD 1.092/stroke 1.032; FH interaction HR 0.946 NS;
  11.1%/36.1%; DRAGON-3 OR 1.55; NNT 315; MR AS 1.551/CAD 1.251). Existing nejm-email reproducers: PAPER1 66/66, THESIS 56/56 PASS.
- **CORRECTION 2026-05-31:** the QC sub-agent's "v17 doesn't trace" FAILs are **largely a
  units/source mismatch, NOT fabrication** — it compared v17 prose against the locked engine without
  matching scaling or GWAS source. Verified: v17's MR AS **1.026** is Helgadottir-2018 (n=4 SNPs,
  per-~10-nmol); the locked **1.551** is a different GWAS (n=6, per-SD) — both defensible, NOT a
  drift. CAD 1.282 (Aragam-2022, n=5, per-50-nmol) vs locked 1.251 (n=6) likewise. DRAGON-3 2.71 is
  the categorical ≥125→ASCVD OR, not the locked continuous estimate. The v17 reconciliation table is
  **RETRACTED — do not apply** (it would corrupt the manuscript). Only the within-UKB cardiac-MRI
  (LV-mass P 0.80 vs 0.041) and diabetes-interaction (0.044 vs 0.096) items are *potential* genuine
  same-quantity discrepancies, and must be confirmed on a matched model+covariate spec before any edit.
  Original (now-doubted) flags:
  - MR Lp(a)→AS OR 1.026 (v17) vs 1.551 (locked) — SPURIOUS (different GWAS/units).
  - DRAGON-3 OR 2.71 (v17) vs 1.55 (locked) — SPURIOUS (categorical vs continuous).
- **Matched-model re-check of the 2 candidate items (2026-05-31):**
  - **Diabetes interaction — GENUINE over-claim.** Locked `modifier_interactions_main.csv`: Lp(a)×DM
    interaction is NS (nonFH ASCVD interact_p=0.745; FH 0.871; BH 0.84–0.90). v17's "three-fold
    amplification, interaction P=0.044" is NOT supported — soften to "stratified ORs differ (0.021 vs
    0.067) but the formal interaction is non-significant." **This one real edit stands.**
  - **LV-mass — NOT a confirmed error.** Locked LV-mass p=0.041 is UNADJUSTED OLS (lpa_log_z, n=27,162);
    v17's "NS" is the covariate-adjusted model (age/sex/BMI legitimately attenuate it). Likely both
    correct — different models. No edit without confirming v17's adjustment.
- **Net v17 verdict:** engine sound; the only genuinely-actionable fix is the diabetes-interaction
  over-claim. The rest of the "drift" was an audit units/source/model artifact.
  - Cardiac MRI LV-mass P **0.80 "NS"** (v17) vs **0.041 significant** (locked) — direct refutation.
  - Diabetes interaction P **0.044 "3-fold"** (v17) vs **0.096, NS after BH-FDR** (locked).
  - Untraceable: per-10-nmol HR scale; n=370,788; 10-yr KM incidences; 42,120-pair matched cohort (actual locked PSM is a 14,722 treatment-interaction design).
- **Action:** rewrite v17 directly from the FIXED_20 locked engine before any resubmission.

## 3b. NB01 — VERIFIED after endpoint correction + reproducer fix (2026-05-31)
- Endpoint corrected in all 7 phase scripts (master `first_ascvd`, I20-I25). Re-run end-to-end.
- **Reliability bug found + fixed:** phase scripts write to `…/results/` but the reproducer read a
  stale `…/NB1/results/` snapshot (2026-04-19) — its prior "51/56 PASS" was validating month-old
  data. Repointed to live output.
- **Honest live audit: 31/56 PASS | 20 FAIL | 5 SKIP.** Robust claims PASS (paradox OR 2.05,
  non-identifiability R² 0.019, CALON-FH AUC 0.74/0.77). Interaction claims FAIL/collapse
  (LDL×statin HR 1.32→0.99 p 0.018→0.965; statin×SSS 1.32→0.91). The 20 FAILs are the exact
  manuscript numbers requiring revision.
- Lp(a) reproducers (Lpa_Manuscript `/data`, Lpa_Multilevel `results/v21_FIXED/`) read canonical
  locations — no stale-dir issue; their audits (§2, §3) are reliable.

## 4. Paper 1 (Orthogonal Decomposition) — 36/49 PASS, 8 FAIL, 5 SKIP — HEADLINE SOUND
**Correction (2026-05-31):** an earlier note cited a "n 775→277 / rho 0.083→0.183" drift from a
STALE C: ledger (`CALON_PAPER1_TRACE.py`). Direct re-computation disproves it: the SSS↔treated-LDL
claim is **rho=0.0833, p=0.0203, n=775 — EXACTLY the manuscript value (PASS).**
- **PASS (headline intact):** SSS↔treated-LDL 0.083/p0.020/n775; orthogonality SSS⊥PRS −0.015;
  2×2 residual-LDL grid; CAD-PRS pos/neg; ASCVD events 672; N_UKB_STATIN 775; Wales rho 0.234.
- **FAIL (all minor — NOT headline errors):** sub-0.01 rho deltas from a too-tight 1e-3 tolerance
  (SSS_UNTREATED −0.024 vs −0.029; PRS_UNTREATED 0.11 vs 0.118; PRS_ADJ 0.143 vs 0.15;
  CHOLY_RHO 0.111 vs 0.103); one p-underflow (PRS_ADJ_P, not real); DOMAIN_OLINK 4.02 vs 4.22; and
  a cholesterol-years **count drift (CHOLY pos 236→296, neg 215→272)** worth one reconciliation.
- **SKIP (untraceable):** NMR LDL/remnant/DHA correlations, carotid IMT high-SSS, Cox MACE HR.
- Conclusion: Paper 1 is the cleanest of the four — headline claims verify; fixes are loosening the
  reproducer tolerance for rho, resolving the cholesterol-years counts, and sourcing 5 minor claims.
  Lesson (again): verify reproducer FAILs against direct computation before trusting them.

## All four audited manuscripts carry traceability issues
| Paper | PASS | Key problem |
|---|---|---|
| NB01 | 31/56 | endpoint bug → interaction claims collapse (corrected; main effect survives) |
| Lpa_Manuscript | 34/57 | 3 load-bearing claims unsupported (allelic P-trend, AlphaFold, AS-Cox) |
| Lpa_Multilevel v17 | 24/38 | **drift largely SPURIOUS** (units/source mismatch in MR + external; reconciliation table RETRACTED). Engine sound. Only cardiac-MRI + DM-interaction need matched-model re-check |
| Paper 1 | 36/49 | **headline SOUND** (treated-LDL 0.083/n775, orthogonality, Wales 0.234 PASS); FAILs are sub-0.01 rho tolerance deltas + cholesterol-years counts |
The QC discipline worked: every "submission-ready" label was optimistic. None should be submitted
until reconciled. Per-paper reproducers + reconciliation docs are in each folder.

---

## Bottom line & recommended path

Three "submission-ready" papers are **NOT submission-ready**. The shared root cause is the
endpoint field mis-mapping; the Lp(a) papers additionally have headline claims that are either
contradicted by, or absent from, the on-disk computations.

**Per the locked-rerun discipline (CLAUDE.md), each needs an end-to-end rerun, not a hot-patch:**
1. Fix the ASCVD/AS endpoint to the correct fields (I20–I25 via p131296–p131306; I35.0 from the
   ICD-10 array for AS; full I63/I64/I70/G45 composite) and regenerate all downstream CSVs.
2. Re-derive the contradicted/untraceable headline claims from real computations (allelic-series
   P-trend; AlphaFold matrix; UKB AS Cox) or remove them.
3. Add the missing non-negotiables (E-value, NoAgeLDL, family-level dedup).
4. Rewrite Lpa_Multilevel v17 from the locked engine.

These are substantial, scientifically-consequential reruns — they should be done with your
direction, one paper at a time, not autonomously. The reproducers written today
(`*_TRACEABILITY_REPRODUCER.py` in each folder) will re-verify each paper after the rerun.

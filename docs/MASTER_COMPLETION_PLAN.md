# Master Completion Plan — finishing all projects

**Compiled:** 2026-05-30 (autonomous assessment of all project folders, verified against source files).

This is the single decision-ready roadmap for finishing the whole programme. Each project has a
status, the work I can safely complete autonomously, and the items that **need your decision**
(journal choice, ORCID/funding, scientific framing, external data, submission). Numerical/agent
claims have been verified against the actual files where feasible — corrections are flagged.

> **Honesty note on autonomous scope.** I will not silently rewrite manuscript claims, reframe
> features, or fabricate model parameters — those are in your documented red-flag zone. Such items
> are listed as `[NEEDS-USER]` even where technically "doable". I only auto-complete work whose
> correctness is verifiable (clear bugs, code that runs against an existing reference, QC).

---

## Status at a glance

| Project | Status | Headline outstanding |
|---|---|---|
> ⚠ **2026-05-30 QC UPDATE:** manuscript-QC audits (see `QC_CRITICAL_FINDINGS.md`) found that
> NB01, Lpa_Manuscript and Lpa_Multilevel v17 are **NOT submission-ready** — a confirmed UKB
> endpoint mis-mapping bug (p131286–294 = hypertension, not ASCVD) plus untraceable/contradicted
> headline claims. Status rows below revised accordingly.

| **CALON Paper 1** (Drug Response) | QC INTERRUPTED | audit not completed; re-run QC `[AUTO]` |
| **NB01** (Discordance Paradox) | ✅✅✅ CORRECTED MANUSCRIPT COMPLETE | endpoint corrected (7 phases) + reproducer fixed (31/56) + E-value + NoAgeLDL + duration-OR verified + ALL descriptives recomputed (LDL-excess, NMR/TG/HDL). **Full corrected manuscript draft — zero placeholders**, every number verified on corrected endpoint. Thesis collapse comprehensive (interactions, descriptives all NS); robust core (OR 2.05, non-id, CALON-FH) intact. **Reproducer now certifies the corrected draft: 51/56 PASS, 0 FAIL, 0 CRITICAL, exit 0** (stale flags resolved). Only thesis-framing decision + journal submission left `[your call]` |
| **NB02** (Penetrance 3D) | ✅ ANALYSIS CORRECTED + VERIFIED | endpoint fixed (nb2_04/05) + variant-key fix. **Verified: overall CI@65 36.1%→14.0% (correct this number); but Law-3 SSS×age pattern HOLDS** (age-45 SSS-high 1.99% vs low 1.16%, OR~1.7, attenuates by 65). See NB2_ENDPOINT_CORRECTED_FINDINGS.md. Remaining: nb2_06 calibration `[NEEDS-USER]`; doc-number update |
| **NB03–NB10** | SCAFFOLD / NOT-STARTED | NB07–09 ready to code `[AUTO]`; NB05/06 blocked on external data `[NEEDS-USER]` |
| **Paper 2** (VUS Reclass.) | SCAFFOLD-ONLY | initialise codebase (overlaps NB08) `[AUTO]` |
| **Paper 3** (ASCVD) | DRAFT (7 versions) | choose version + 7-flaw corrections `[NEEDS-USER]` |
| **TUDOR** | R3 in prep | reproducer extension `[AUTO]`; 16 DRIFT reconciliation + manuscript fixes `[NEEDS-USER]` |
| **Lpa_Manuscript** | ⚠ headline collapsed — honest reframe drafted | 3 load-bearing claims unsupported (allelic P-trend triple-verified NS 0.25–0.69 vs claimed 0.001; AlphaFold + AS-Cox = literals, no source). **Honest reframe draft written** (Lpa_Saturation_REFRAMED_DRAFT.md): verified Dragon-3 layer retained, 3 pillars withdrawn, 3 options laid out. Paper viability/direction = your scientific decision `[NEEDS-USER]` |
| **Lpa_Multilevel** | ✓✓ v18 CORRECTED DRAFT | engine sound; v17 "drift" was audit artifact (table retracted). The one real edit (diabetes-interaction over-claim) **applied to v18 draft** (Lpa_Multilevel_Manuscript_v18_DM_CORRECTED_DRAFT.md): all 3 locations softened to NS-verified, three-fold claim withdrawn. Remaining: your review + journal resubmission |
| **Lpa** (AS Fine-Gray) | DATA-READY | **MR-pathway decision C-1/C-2/C-3** `[NEEDS-USER]`, then run `[AUTO]` |
| **CALON-FH Atlas** | DEPLOYED | — |

---

## Work completed in this session

- **NB02 / nb2_05_per_variant.py — FIXED.** Root cause: the carrier extract
  (`ukb_carriers_FINAL.csv`) has **empty** `variant_id`/`cdna`/`protein` columns, so both the
  original and the so-called `*_CORRECTED` reference produced no per-variant result (the
  "corrected" reference was a meaningless 2-row by-`consequence` table). Fix: build the unique
  variant key from the fully-populated genomic coordinates `chrom:pos:ref:alt`, annotate with
  `protein_pos`, and guard the empty-result crash. Re-ran: **161 variants, 60 with ≥5 carriers**;
  honest result — no variant survives BH-FDR (all p_BH = 1.0), empirical-Bayes shrinkage factor
  0.08 pulls all toward cohort-mean OR 0.885 (consistent with SSS as a continuous severity axis).
  Output: `NB2/results/nb2_05_per_variant.csv`.
  - *Correction to assessment:* the earlier claim that the file "only has cdna, not variant_id"
    was wrong — it has `variant_id`, but it (and `cdna`, `protein`) are entirely null; coordinates
    are the right key.

---

- **Lpa AS Fine-Gray (track authorised C-1) — partially advanced.**
  - `[DONE]` Cross-sectional AS association ran (honest, underpowered): prevalence 1.84% (65/3,540);
    Lp(a) OR 1.26/SD [0.92–1.74] p=0.154 (NS); tier p-trend 0.636 (NS). Files in `D:\Projects\Lpa\`.
  - `[DONE]` Fine-Gray prep pipeline built; RAP extraction script **bug-fixed** (GRCh38→GRCh37,
    rsID-based) + AS-date block added.
  - `[NEEDS-USER / RAP]` Time-to-event Fine-Gray blocked on AS diagnosis **dates** (absent locally,
    `p41271` empty for all 65 flagged carriers). Run corrected `EXTRACT_LPA_SNPS_ON_RAP.sh` on RAP
    for `lpa_snps.csv` + `as_event_dates.csv`, then re-run prep. See `D:\Projects\Lpa\STATUS_2026-05-30_finegray.md`.

## Per-project remaining work

### CALON Paper 1 — SUBMISSION-READY
- `[NEEDS-USER]` Select target journal; fill journal/editor names in `manuscript/cover_letter.txt`.
- `[NEEDS-USER]` Format supplementary materials to journal style.

### NB01 — JCL-READY
- `[AUTO]` Age-strata (<65 vs ≥65) sensitivity to confirm the LDL×statin inversion holds.
- `[NEEDS-USER]` Choose journal (*Circulation* / *JAMA Cardiology*) and submit.

### NB02 — IN-PROGRESS
- `[DONE]` nb2_05 per-variant (above).
- `[AUTO]` Re-run nb2_02→nb2_04 on the CORRECTED cohort if their outputs predate the cohort fix
  (verify timestamps first; do not overwrite manuscript-feeding CSVs without a diff check).
- `[NEEDS-USER]` nb2_06 ascertainment simulation — the log-onset coefficients (−0.35·SSS,
  −0.05·PRS) are ad-hoc. Calibrating them to published FH onset ages (Mundal 2020, Besseling 2016,
  Pérez de Isla) is a **modelling choice you should own** — I can implement once you confirm the
  target onset distribution.
- `[NEEDS-USER]` Reframe the "Law 3 Biobank Event Horizon" narrative for the corrected event
  counts (22.2% vs 36.1% cumulative incidence at 65).

### NB03–NB10
- `[AUTO]` NB07 (pre-event NMR), NB08 (VUS reclassification — inputs on disk), NB09 (sex-specific)
  are codeable now from their design specs.
- `[NEEDS-USER]` NB05/NB06 blocked on FinnGen / All of Us access (institutional).
- `[NEEDS-USER]` NB03 needs GTEx download; NB10 lowest priority.

### Paper 2 (VUS) — SCAFFOLD
- `[AUTO]` Initialise codebase reusing SSS v3 + ClinVar (shares logic with NB08).

### Paper 3 (ASCVD) — DRAFT
- `[NEEDS-USER]` Choose the final manuscript version (v5/v6/v7); confirm the 7-flaw corrections
  (prevalent-vs-incident contamination, treatment-variable handling, LDL-direction "Raal
  correction") — these are scientific-framing decisions.

### TUDOR — R3 in prep
- `[AUTO]` Extend `TUDOR_TRACEABILITY_REPRODUCER.py` to assert all locked coefficients and
  calibration intercept/slope (currently asserts AUC/NRI/Se/Sp/n only) → target ≥56/56 PASS.
- `[NEEDS-USER]` The **16 DRIFT ledger items** are documented as *unreproducible* under 75
  methodology combinations (e.g. Wales NRI 0.358→−0.097, UKB calib slope 6.33→1.23). These must be
  replaced with live values or cited as superseded — a manuscript-integrity decision for you.
- `[NEEDS-USER]` Manuscript fixes flagged by QC (coefficient list §2.5, Design-B "frozen" framing,
  Trig_Filter narrative, WES methods detail) — all involve scientific framing/claims.
- `[NEEDS-USER]` 45-item author checklist (ORCID, CRediT, ethics, cover letter, figure DPI).

### Lpa_Manuscript — SUBMISSION-READY
- `[AUTO]` Resolve cardiac-MRI aortic-volume audit (re-derive from raw CMR p22420–p22425 or retain caveat).
- `[NEEDS-USER]` Funding + competing-interests + co-author affiliations (currently "[to be confirmed]").

### Lpa_Multilevel — UNDER REVIEW (v17, *Atherosclerosis* ATH-D-26-00700)
- `[AUTO]` Add 9 new 2025–26 citations; rewrite incomplete `v21_FIXED_23` trial-enrichment script;
  add MVMR cross-effects caveat.
- `[NEEDS-USER]` Confirm which version was actually submitted (editor disclosure needed only if ≤v15).

### Lpa (AS Fine-Gray) — DATA-READY
- `[NEEDS-USER]` **MR-pathway decision:** C-1 (extract rs10455872 + rs3798220 from RAP), C-2
  (two-sample summary-stat MR), or C-3 (skip MR, observational Fine-Gray only).
- `[AUTO]` Once chosen: run the AS Fine-Gray regression + sensitivity ASCVD composite.

---

## Autonomous work — VERIFIED status (corrected 2026-05-30 after checking each)
1. ✅ **NB02 nb2_05 fix + full pipeline re-run** (nb2_01→06 all execute on corrected cohort). DONE.
2. ✅ **Lpa cross-sectional AS** + Fine-Gray prep + RAP-script bug-fix. DONE to local-data limit.
3. ❌ **TUDOR reproducer extension** — BLOCKED: the traceability reproducer + locked JSON live in
   the `.claude/worktrees/fervent-mayer-52a6f3` worktree, which is an **actively running session**.
   Editing it risks clobbering live work. Do this from inside that session, not here.
4. ❌ **NB01 age-strata sensitivity** — no existing script; building it = inventing the analysis.
   NB01 is otherwise JCL-ready (gated on journal choice).
5. ❌ **NB07/NB08/NB09 + Paper 2** — empty scaffolds, **no methodology spec written**; the project
   notes are explicitly skeptical of autonomously expanding the Nobel series. Needs your design.
6. ⚠ **Lpa_Manuscript CMR re-derivation** — UKB has no "aortic root volume" field; available aortic
   fields (p24118–24123) need UKB-dictionary verification before use (field-mislabel risk). Needs
   your confirmation of which measure the manuscript intends.

> Honest conclusion: the safe, established autonomous work is essentially complete. Everything
> remaining requires a decision, a written spec, an external action (RAP / journal), or access to
> the running TUDOR session. Listed below.

## What only you can close (the true "finish" gate)
- All journal selections + submissions, ORCID/CRediT/funding/affiliations.
- TUDOR DRIFT reconciliation + manuscript framing decisions.
- NB02 nb2_06 calibration target + narrative; Paper 3 version + flaw-correction framing.
- Lpa MR-pathway choice; external-data access for NB05/NB06.

# Skills — the CALON-FH / TUDOR working toolkit

These are the reusable **skills** (procedural playbooks) that drive this project's analysis and
QC discipline. On Dr Genedy's machine they live in `C:\Users\nader\.claude\skills\`. They are
documented here from scratch so the methodology is portable even without the skill files.

Each skill is invoked when its trigger phrases appear. They encode the project's
**non-negotiable methodology** so that quality is enforced by process, not memory.

---

## 1. Data & cohort skills

### `ukb-inventory`
**When:** locating any UK Biobank field, cohort size, NMR/PRS/imaging measure, or result CSV.
**Does:** single source of truth for the three-tier storage layout — the 501,936-row master,
the 32 result CSVs, headline cohort numbers, and outstanding RAP extractions. Cross-checks
candidate fields against the UKB field catalogue (`fieldsum.tsv`) before declaring "not available".

### `ukb-preflight`
**When:** BEFORE launching any pipeline that loads CSVs, fits Cox models, or writes results.
**Does:** verifies environment (Python 3.12 + required packages + cp1252 encoding guard) AND
inventories on-disk data to avoid duplicate RAP extractions. Saves the retry cycle when pyarrow
is missing or the master CSV has not been re-extracted.

### `stage-files`
**When:** "collect", "stage", "gather", "bundle", "prepare files for transfer".
**Does:** **collect first, inventory second.** Hard-link → reflink → copy fallback chain;
detects locked Word/Excel files (`~$name`); dedups by 1 MB sha256; emits `staging_manifest.csv`.
*(This bundle was produced with this skill.)*

---

## 2. Analysis skills

### `cox-analysis`
**When:** any Cox proportional-hazards / survival / hazard-ratio / incident-ASCVD work.
**Does:** enforces the methodology TDD checklist — **leak-free splits, family-level
deduplication, NoAgeLDL sensitivity variant, OPCS-4/ICD-10 completeness (incl. K611),
TRIPOD-compliant reporting** (C-stat + CI, calibration, Brier, NRI, IDI, DCA). Checks the PH
assumption, flags immortal-time bias, collider bias, treatment paradox, Table-2 fallacy, VIF>5.

### `comparator-model-validation`
**When:** any claim of the form "our model beats / outperforms X" (SAFEHEART, Montreal-FH-SCORE,
DLCN, FAMCAT, PCE, SCORE2). Invoke **before** writing the comparison.
**Does:** enforces frozen published comparator coefficients, like-for-like out-of-sample
comparison, a pre-registered sign convention, full (uncherry-picked) subgroup tables, real
NRI/IDI computation, and traceability of every reported metric to a live results file.

---

## 3. QC & reproducibility skills

### `manuscript-qc`
**When:** "QC the manuscript", "audit numbers", "are all results traceable", before submission.
**Does:** validates every numerical claim against source CSVs; produces an automated reproducer
script + a per-paper provenance table. Goal: **100% PASS / 0 FAIL**.

### `verify-manuscript` / `autonomous-qc` / `ukb-qc-fivefold`
**When:** end-to-end manuscript QC with **parallel sub-agents** (one per verification domain).
**Does:** dispatches independent agents — numerical-claims, code-completeness, figure-reproducer,
cohort-validator, provenance-detector — each returning PASS/DRIFT/FAIL, then reconciles. Catches
the failure mode where a number was correct once, hand-copied into prose, and the pipeline then
drifted. `ukb-qc-fivefold` is the five-agent UKB-specific variant.

### `tudor-qc`
**When:** auditing the TUDOR pipeline end-to-end, or testing whether a candidate feature (T2DM,
age-per-decade, LDLR genotype tier, ApoB-continuous) improves discrimination/calibration/NRI.
**Does:** dispatches five agents — Raw Data Integrity, Statistical Methodology, Clinical
Plausibility, Molecular/Genetic, Feature-Augmentation — then a unified report. Designed to catch
the WRONG_FIELDS class of bug (UKB date fields p131286–p131294 silently coding I10–I20
hypertension, not ASCVD), plus treatment-adjustment, family-leakage, calibration-misweighting.

### `ukb-data-audit`
**When:** auditing UKB extracts, hunting code-list bugs, tracing manuscript claims to CSV rows,
verifying cohort filters before a portal amendment.
**Does:** five parallel agents — inventory, field verification, code-list completeness,
manuscript provenance, cohort consistency. Catches "ledger says missing but file is in another
folder", field-mislabel (p24100 CMR-not-CT), omitted K611, and HR-not-matching-any-CSV-row.

### `reproduce-paper` / `reproducibility-tests`
**When:** build/run a frozen ground-truth reproducer for a paper, or encode invariants as pytest.
**Does:** `reproduce-paper` re-runs the analysis end-to-end and asserts every published number
against fresh computation (the "110/110 PASS" pattern). `reproducibility-tests` converts those
invariants into a continuous pytest suite that fails CI when one breaks.

---

## 4. Output & figure skills

### `md2docx`
**When:** any Word version of a manuscript / reply / cover letter / feedback.
**Does:** dual-format default (md + docx side-by-side), Cambria/Arial fonts, robust find-and-
replace across **fragmented Word runs** (verify substitution count > 0), preserves embedded images.

### `figure-audit`
**When:** auditing AI-generated schematic figures (Gemini/DALL-E/Midjourney) before submission.
**Does:** catches text artefacts, garbled labels (the "VLDL lLopiite" / "abundatl" hallucination
class), impossible anatomy, unit-formatting errors.

### `student-feedback`
**When:** academic feedback for Cardiff medical students.
**Does:** Cardiff academic register, 3-paragraph structure (strength → development → action),
British English, per-student md+docx plus a cohort summary file.

---

## 5. Domain expert skills (always-active perspectives)

`academic-writing`, `biostatistics`, `cardiology`, `metabolic-medicine`, `genomics`,
`molecular-genomics`, `physiology`, `graphic-design`, `nature-reviewer`, `nobel-medicine`,
`data-analysis` — applied simultaneously as analytical lenses (survival stats, lipidology,
LDLR molecular pathway, ACMG/AMP classification, AlphaFold confidence bands, publication figures).

---

## How the skills enforce quality

The skills exist because of specific past failures:
- A UKB endpoint-mapping bug (date fields coding hypertension not ASCVD) → `tudor-qc`, `ukb-data-audit`.
- Manuscript numbers drifting from the live pipeline → `manuscript-qc`, `reproduce-paper`, traceability ledgers.
- Garbled AI figure text reaching a *Heart* submission → `figure-audit`.
- Inventorying before staging when asked to collect → `stage-files` (collect-first rule).
- Comparator "we beat X" claims without frozen coefficients → `comparator-model-validation`.

See `docs/INSTRUCTIONS.md` for the underlying statistical non-negotiables these skills enforce.

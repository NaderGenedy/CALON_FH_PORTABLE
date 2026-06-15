# Instructions & Conventions — how this project is run

These are the standing rules for the CALON-FH / TUDOR programme. They are applied **by default**,
not on request. They exist to keep the analysis reproducible and the manuscripts defensible.

---

## A. Statistical non-negotiables (every FH / UKB analysis)

1. **Family-level deduplication.** South Wales (Dragon-3) is a subset of All-Wales PASS by
   `FamilyNumber`. Always remove FamilyNumber overlap, not merely DatabaseNumber overlap.
   One event per family; never per-person without a sibling check.

2. **NoAgeLDL sensitivity variant.** For any FH risk-prediction or discordance model, run a
   variant with age + LDL removed, to test whether novel predictors carry signal beyond the
   dominant traditional backbone.

3. **Code-list completeness audit before any cohort extraction.**
   - **ASCVD composite (ICD-10):** I20, I21, I22, I23, I24, I25, I63, I64, I70, G45 — plus the
     matching first-occurrence fields p131286 (I21), p131288 (I20), p131290 (I63), p131292 (G45),
     p131294 (I70), p131296 (I25). **Never use I35 (aortic stenosis) for an ASCVD endpoint.**
   - **OPCS-4 valve:** must include **K611 (balloon valvuloplasty)** — omitting it loses ~3,000
     severe-AS cases. K40–K46 CABG; K49/K50/K75 PCI.

4. **Endpoint provenance check.** In UKB the date fields p131286–p131296 code specific I-codes;
   confirm each maps to the intended endpoint (the historic bug had them silently coding
   I10–I20 hypertension rather than ASCVD).

5. **TRIPOD compliance for prediction models:** C-statistic + 95% CI (DeLong or bootstrap),
   calibration intercept + slope, Brier (raw + scaled), NRI (bilateral bootstrap B≥100), IDI,
   DCA net benefit at 5/10/20% thresholds, DeLong p vs comparator. Prefer a frozen-coefficient
   external comparator (TRIPOD Type 4).

6. **Multiple-testing correction.** Benjamini–Hochberg FDR for any panel of >10 simultaneous
   tests (NMR field correlations, PRS panels, subgroup forests). Bonferroni when independent.

7. **E-value sensitivity** for every headline OR/HR (VanderWeele & Ding 2017). E-value > 2.0 =
   "robust to plausible unmeasured confounding".

8. **Leak-free / bi-external validation.** Train and validation cohorts share zero family IDs.
   Report internal and external discrimination separately.

9. **Validate every numerical claim against source CSV** before "submission-ready". Build a
   per-paper provenance table (one row per claim → file:row:n:events:date) and an automated
   reproducer that asserts manuscript value == CSV value within tolerance. Aim for 100% PASS.

10. **Locked rerun discipline.** When a methodological bug is found, do an end-to-end locked
    rerun (not a hot-patch), refresh all downstream CSVs, append to the PROVENANCE ledger,
    and update the reproducer.

---

## B. Analytical discipline (interpretation)

- **Report numbers first, propose mechanisms after.** Do not jump from one regression
  coefficient to a mechanism.
- **Label every claim:** "supported by data" / "consistent with data" / "speculative".
- **QC your own interpretation** — if the framing changes mid-analysis, say so and why.
- **Diminishing returns** — after ~20 sensitivity analyses that converge, stop. Further tests
  are noise-farming.
- **Red flags (stop):** claiming an untested mechanism; ignoring multiple-testing correction;
  writing "proof" when you mean "consistent with"; framing a null as "borderline" or
  "approaching significance"; inflating a paper's claims to hit a higher journal.
- Realistic targets: **JCL, Circulation, JAMA Cardiology**. Nature Med / NEJM need multi-cohort
  replication or mechanistic data not yet available.

---

## C. Environment & encoding (Windows)

- **Python 3.12 only** — `C:/Users/nader/AppData/Local/Programs/Python/Python312/python.exe`.
  3.14 lacks pyarrow/fastparquet wheels and breaks UKB workflows. Do not silently fall back.
- **R 4.5.x** — `/c/Program Files/R/R-4.5.2/bin/Rscript.exe` (not on PATH; call by absolute path).
- **Console is cp1252.** Never print raw Unicode arrows (→), em-dashes (—), or emoji. Use ASCII
  (`->`, `--`, `[OK]`, `[FAIL]`). For UTF-8 stdout: `sys.stdout.reconfigure(encoding='utf-8')`.
- **File I/O:** always `open(..., encoding='utf-8')` — the default codec silently corrupts
  Greek (β, ρ, χ²) and superscripts.
- **Subprocess:** pass `encoding='utf-8', errors='replace', env={'PYTHONIOENCODING':'utf-8', ...}`.

---

## D. Workflow conventions

- **"Collect / stage" means COPY first, inventory after** — never start with `ls`/`find`.
- **Before editing .docx:** check the file is not open in Word (locked → `python-docx` fails
  silently). `python-docx` find-replace returns 0 on fragmented runs — merge runs or edit
  `<w:t>` XML directly; verify substitution count > 0.
- **Long outputs:** write to files in `./out/` or the manuscript dir; print only 5–10-line
  summaries. Process in batches of 5, pause for OK between batches.
- **Disk-first, RAP last.** Exhaustively inventory local files before any UKB RAP extraction;
  most planned pulls turn out unnecessary.
- **Naming:** `<topic>_<version>_<purpose>.py` — e.g. `LOCK_RERUN_composite_ASCVD.py`,
  `paper4_ukb_apob_reproducer.py`. Lock dates in the docstring.

---

## E. Deliverables

- **Always dual-format:** every manuscript / reply / cover letter / feedback ships in BOTH `.md`
  AND `.docx`, same stem.
- **Manuscript QC is a hard gate** — every HR, CI, P-value, n, % traced to source CSV and
  recomputed before "submission-ready". Output: `audit_report.md` with N/N certified count.
- **Figures:** 300 DPI PNG (inline) + vector PDF (journal upload) — both files, both paths printed.
- **Reproducer scripts** must hit a literal pass-count target (e.g. "110/110 PASS"). The script
  is the deliverable, not the number in chat.

---

## F. Style

British English (colour, analyse, organisation). Nature/Lancet/NEJM register: active voice,
precision over verbosity, one claim per sentence, one hedge word maximum, every claim carries
a number + CI + P-value. IMRAD. Direct and clinical — no filler, no apologies.

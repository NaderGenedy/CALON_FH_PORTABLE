# CALON-FH / TUDOR — Portable Project Bundle

**Dr Nader Genedy** — Cardiology & Clinical Research, University Hospital of Wales, Cardiff
**Bundle generated:** 2026-05-30
**Source repo:** `github.com/NaderGenedy/calon-ukb-pipeline` (local: `C:\Users\nader\Downloads\calon_ukb_pipeline`)

This is a **self-contained, portable snapshot** of the CALON-FH familial-hypercholesterolaemia
research programme and the TUDOR diagnostic-algorithm pipeline. Everything needed to read,
understand, re-run, or hand to another device — or to an online Claude session — is in this
one folder. No external repository checkout is required.

---

## What this bundle contains

| Folder | Contents |
|---|---|
| `data/` | **All 266 CSVs used by the project** (315 MB), repo-relative paths preserved. See `docs/DATA_DICTIONARY.md`. |
| `code/python/` | All 145 Python analysis scripts (`17_*.py` … `73_*.py`, `CALON_*.py`, `TUDOR_*.py`) |
| `code/R/` | All 76 R scripts (`01_*.R` … `16_*.R`, TUDOR modelling, figures) |
| `code/shell/` | All 18 shell scripts (UKB-RAP extraction, helpers) |
| `manuscripts/` | All `.docx` / `.tex` / `.pdf` manuscripts, cover letters, tables, response-to-reviewers |
| `docs_md/` | All original `.md` documents from the repo (onboarding, traceability, field dictionaries) |
| `assets/` | Figures (PNG/PDF/SVG), JSON, AlphaFold structural outputs, miscellaneous artefacts |
| `website/` | The CALON-FH Atlas interactive site (`index.html`, data JSON) |
| `projects/` | **All other projects** mirrored from `D:\Projects` — NOBEL series NB01–NB10, Lp(a) sub-programme, MD thesis, risk models, AlphaFold SSS, shared data, teaching. Each preserves its own structure. See `docs/PROJECTS_OVERVIEW.md`. |
| `docs/` | **Read these first** — authored documentation (below) |
| `staging_manifest.csv` | Per-file lineage: source → destination, copy method, size, hash |
| `docs/LARGE_DATA_INDEX.csv` | Giant raw CSVs (> 300 MB) not copied — original D: paths, so nothing is lost-to-knowledge |

### The four documents to read first (`docs/`)

0. **`docs/PROJECTS_OVERVIEW.md`** — map of every project in this bundle (root pipeline + the
   25 folders under `projects/`). Read this to see the whole programme at a glance.
1. **`docs/WORK_DONE.md`** — the full narrative of the root pipeline (CALON-FH / TUDOR / Atlas):
   what was built, why, the validated findings, and how each piece fits together.
2. **`docs/INSTRUCTIONS.md`** — the working conventions and statistical non-negotiables
   (family-level dedup, NoAgeLDL sensitivity, ICD-10/OPCS code lists, TRIPOD, E-values).
3. **`docs/SKILLS.md`** — the 16 project skills (cox-analysis, tudor-qc, manuscript-qc,
   ukb-inventory, …): what each does and when to invoke it.
4. **`docs/DATA_DICTIONARY.md`** — every CSV profiled (rows × cols, size, headers), grouped.

5. **`docs/SESSION_HISTORY.md`** — the work catalogued discussion-by-discussion (26 Claude Code
   sessions, chronological), showing when and in what order each piece was built.

Plus **`CLAUDE.md`** at the root — paste-ready context to bootstrap an online Claude session.

---

## How to use this on another device

1. **Copy the whole folder** (`D:\CALON_FH_PORTABLE`) to the target machine / drive / USB.
2. Install Python **3.12** (not 3.14 — pyarrow/fastparquet need 3.12 wheels) and R 4.5.x.
3. Install packages: `pip install pandas numpy scipy lifelines statsmodels scikit-learn pyarrow python-docx`
   and in R: `install.packages(c("survival","glmnet","ggplot2","patchwork","mice"))`.
4. Open the folder in your editor / Claude Code. The `code/` scripts reference data by path —
   see the **Path note** below.

## How to use this with online Claude (claude.ai)

1. Open `CLAUDE.md` and paste its contents into a new conversation — that gives Claude the full
   project context, validated findings, and conventions.
2. Upload the specific CSV(s) and script(s) you want to work on (online Claude has per-file
   upload limits — upload the relevant slice, not the whole `data/` tree).
3. For large analyses, the small/medium derived CSVs (Table*, *_results, *_odds_ratios,
   coefficients, traceability ledgers) are the ones you'll usually need — the 70–92 MB
   full-cohort files (`CALON_clinical_analysis_dataset.csv`, the `tudor_loco_output/` and
   `tudor_flaw_results/` supplements) are too large for most online uploads but are present
   here for local re-runs.

---

## ⚠ Path note (important for "run straight away")

Many scripts in `code/` were written against **absolute Windows paths** on Dr Genedy's machine
(e.g. `D:/Projects/CALON_AlphaFold_Rebuild/data/...`, `C:/Users/nader/...`). On a fresh device
those paths will not resolve unmodified. To re-run:

- The **data each script needs is in `data/`** in this bundle (same filenames, structure preserved).
- Update the input path at the top of the script you are running to point at this bundle's
  `data/` folder, or set a working directory and use relative paths.
- The `docs/DATA_DICTIONARY.md` maps every filename to its location here.

This bundle is optimised for **completeness and portability of the material** (all code + all data
+ full documentation in one place), not for zero-edit execution on an arbitrary machine.

---

## The three projects in this bundle

- **CALON-FH** — structural-genomic FH severity scoring (SSS v3, 9-layer composite over 16,340
  LDLR variants) and its clinical validation in UK Biobank (3,544 carriers) and Wales cohorts.
- **TUDOR** — a UK Biobank lipid-threshold FH **diagnostic** risk algorithm (Elastic Net,
  leave-one-cohort-out validation), submitted to *Journal of Clinical Lipidology*.
- **CALON-FH Atlas** — the deployed interactive website (`website/`).

Full detail in `docs/WORK_DONE.md`.

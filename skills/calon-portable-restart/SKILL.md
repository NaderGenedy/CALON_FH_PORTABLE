---
name: calon-portable-restart
description: Use when rehydrating or restarting the CALON-FH / TUDOR programme from the CALON_FH_PORTABLE bundle on a new machine or an online Claude session. Triggers on "restart the project", "rehydrate the bundle", "set up CALON on this machine", "restore the data", "run it on Claude online", or when handed the CALON_FH_PORTABLE repo. Enforces the DUA rule (participant data stays on D:, restored locally; never pushed), points to the data inventory map, the code entry points, and the bootstrap docs.
---

# CALON portable restart

Rehydrate the CALON-FH / TUDOR programme from `CALON_FH_PORTABLE`.

## What is and isn't in the bundle
- **In repo:** `code/` (R + Python + shell), `skills/`, `docs/`, `docs_md/`,
  `manuscripts/`, `website/`, `results_public/` (aggregate tables), `CLAUDE.md`,
  and `DATA_INVENTORY_ON_D.csv` (the map of every data file's D: location).
- **NOT in repo (UKB DUA):** `data/`, `projects/`, `assets/` — participant-level
  data and heavy trees. They live on the D: drive only.

## Restart procedure
1. `git clone` the repo; read `RESTART_GUIDE.md` and `docs/PROJECTS_OVERVIEW.md`.
2. If the task needs participant data, restore `./data/` from the D: drive (or the
   user's secure storage) using the relative paths in `DATA_INVENTORY_ON_D.csv`.
   **Never upload participant data to GitHub or any uncovered cloud.**
3. Environment: Python **3.12** (pandas/numpy/scipy/lifelines/statsmodels/sklearn/
   pyarrow/matplotlib/python-docx), R **4.5.2** (survival/glmnet/mice/ggplot2/rms).
   Keep stdout ASCII on Windows; open files `encoding="utf-8"`.
4. Re-run via the entry points in `RESTART_GUIDE.md` §4 (cohort build → validate →
   develop v5-v7 → head-to-head `11_CALON2_cox_comparison.R` → IECV → TUDOR).
5. Online Claude: select the repo at claude.ai/code, paste `CLAUDE.md`, point at
   `skills/` and `docs/`.

## Non-negotiables (carry over from the programme)
Family-level dedup; NoAgeLDL sensitivity; full ICD-10/OPCS code lists (incl. K611);
TRIPOD reporting; frozen-coefficient comparators for any head-to-head claim;
every manuscript number traceable to a `results_public/` table.

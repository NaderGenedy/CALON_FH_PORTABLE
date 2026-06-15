# How to download & restart the whole project

Everything needed to rehydrate the CALON-FH / TUDOR programme on a new machine or
an online Claude session. **Code, skills, docs and manuscripts are in this repo;
the data stays on the D: drive** (UKB DUA) and is mapped by `DATA_INVENTORY_ON_D.csv`.

---

## 1. Download the code

```bash
git clone https://github.com/NaderGenedy/CALON_FH_PORTABLE.git
cd CALON_FH_PORTABLE
```

## 2. Restore the data (not in the repo — by design)

The repo ships only **aggregate results** (`results_public/`) and a **map** of every
data file (`DATA_INVENTORY_ON_D.csv`, column `abs_path_on_D`). To run analyses that
need participant data, restore it **within your UKB DUA** to a `data/` folder at the
repo root, mirroring the relative paths in the inventory:

- from your D: drive: copy `D:\CALON_FH_PORTABLE\data\` → `./data/`, or
- from your own secure storage: `aws s3 sync s3://<you>/calon-data ./data`.

Large raw extracts not staged locally are listed in `docs/LARGE_DATA_INDEX.csv`
with their original D: paths.

## 3. Set up the environment

| Tool | Version | Notes |
|---|---|---|
| Python | **3.12** | `pip install pandas numpy scipy lifelines statsmodels scikit-learn pyarrow matplotlib python-docx openpyxl` |
| R | **4.5.2** | `install.packages(c("survival","glmnet","mice","ggplot2","patchwork","pROC","rms"))` |

Console: keep stdout ASCII on Windows (cp1252); open files with `encoding="utf-8"`.

## 4. Re-run the work

| Piece | Entry point (in `code/`) |
|---|---|
| CALON-FH cohort build → validate → figures | `R/01_CALON_build_cohort.R` → `02_…` → `03_…` |
| CALON v5–v7 development | `R/05_…` … `R/07_CALON2_develop_v7.R` |
| Head-to-head vs SAFEHEART / Montreal-FH-SCORE | `R/11_CALON2_cox_comparison.R` + `R/99_extract_frozen_coefficients.R` |
| Internal-external cross-validation | `R/15_CALON2_IECV.R` |
| TUDOR model + reproducers | `code/python/TUDOR_*.py` |
| Per-paper traceability reproducers | `python/*_TRACEABILITY_REPRODUCER.py`, `*/run_all_*.py` |

Order and rationale are documented in `docs/WORK_DONE.md` and
`docs/SESSION_HISTORY.md`. Statistical non-negotiables (family-level dedup,
NoAgeLDL, ICD-10/OPCS lists, TRIPOD, E-values) are in `docs/INSTRUCTIONS.md`.

## 5. Restart an online Claude session

1. Open https://claude.ai/code and select this repo (`CALON_FH_PORTABLE`, private —
   your connected GitHub can access it).
2. Paste the contents of **`CLAUDE.md`** (root) — it bootstraps the full programme
   context.
3. Point Claude at **`docs/PROJECTS_OVERVIEW.md`** (the whole map) and
   **`skills/`** (the 16 + restart skill). Then ask for the piece you need; restore
   `data/` first if the task needs participant data.

## 6. What's verified / outperforms

The head-to-head results (CALON / TUDOR vs published FH risk tools) and their
frozen-comparator coefficients are in `results_public/` (e.g. `*_results.csv`,
`*reproduced_AUC*`, `LOCO_CV_results.csv`, `Subgroup_analysis_results.csv`) and the
manuscripts in `manuscripts/` (`TUDOR_Manuscript_v6_TRACEABLE.docx`,
`TUDOR_ResponseToReviewers.docx`). Every headline number traces to a results CSV
via the `*_TRACEABILITY*` / `*_LEDGER` tables.

---

*Data never leaves the secure environment. The repo is the code + map + aggregate
results; participant records are restored locally under your DUA.*

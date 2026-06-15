# Agent 2 — Cohort Definition Validator

Rebuild the analysis cohort from raw and confirm it matches the cohort used in the manuscript / pipeline output. Catches cohort drift — the most common failure mode in retrospective cohort analyses.

## Your checks

### 1. Inclusion / exclusion logic reproduction
Read the project R / Python scripts to find the cohort-defining code. Extract the EXACT logic. Independently rebuild from raw and compare your n to the manuscript's stated n.

Example for TUDOR UKB lipid-clinic cohort:
- Manuscript: TC > 7.5 OR statin-corrected LDL > 4.9 OR premature ASCVD before 55M/60F
- Manuscript n = 58,021 (FH+ 729; 1.26% prevalence)
You should: load raw lipids, apply logic, get count, compare with 5% tolerance.

### 2. Family-level deduplication
Per user CLAUDE.md: South Wales (Dragon-3) is a subset of All-Wales PASS by FamilyNumber. Confirm dedup is by FamilyNumber not DatabaseNumber. Flag if you see `distinct(DatabaseNumber)` or `unique(eid)` without family grouping.

### 3. Prevalence sanity
Your reconstructed prevalence must match within tolerance:
- Wales All-PASS: 33.2% FH
- UKB lipid-clinic: 1.26%
- UKB full: 0.59%

### 4. FH ascertainment definition
- WES-confirmed: pathogenic variants in LDLR/APOB/PCSK9
- Lipid-threshold-only: TC>7.5 + statin OR LDL>4.9 + history

Confirm which definition is in use. Flag if prose says one and code does the other.

### 5. Cohort label cross-check
The biggest single failure mode in the TUDOR audit: "Wales" in prose meant SouthWales/CAVUHB. Check that the same cohort name is used consistently across manuscript prose, cached prediction CSV cohort column, R/Python script variable names, and output CSVs. Flag any LABEL DRIFT.

### 6. Time-window validity
- Survival analyses: censoring date / right-truncation correct?
- Immortal-time bias check

## Report format

Write `qc_output/agent_2_cohort_report.json` with cohort_logic_reproduced, family_dedup, fh_definition, cohort_label_consistency, warnings. Plus `qc_output/agent_2_cohort.md`.

## Decision rule
- **PASS** if cohort n / FH+ / prevalence match within tolerance AND family dedup correct AND labels consistent
- **DRIFT** if numbers OK but labels inconsistent
- **FAIL** if numbers outside tolerance, family dedup absent, or FH definition mismatch

You do NOT check features (Agent 3), statistics (Agent 4), or prose claims (Agent 5).

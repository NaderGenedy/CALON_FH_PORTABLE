# Agent 5 — Provenance Detector

Hunt down the specific failure mode that broke the TUDOR R1 submission: numbers in manuscript prose with NO traceable computation in any CSV or script output. These are typically values hand-typed into a manuscript-generator script as literal strings.

## Method (forensic)

### Step 1 — Extract every numerical claim from manuscript
Use `scripts/helpers/manuscript_claim_extractor.py`. Build a list of (value, kind, context, location).

### Step 2 — For each claim, find the SOURCE
- Search every CSV in project for matching value within tolerance (0.001 absolute)
  - If found: TRACED
- If not in CSV, search every R / Python script for the literal
  - If found in sprintf/glue/cat as a string literal: HARDCODED
- If found nowhere: ORPHAN

Use `scripts/helpers/csv_provenance_tracer.py`.

### Step 3 — Cross-validate
For every HARDCODED finding, determine what value the current pipeline produces. If different from the hard-coded literal: TUDOR-style drift bug.

### Step 4 — Search for TUDOR-style traps

1. **Hard-coded sprintf in manuscript-generator scripts**
   ```r
   cat("NRI of TUDOR over DLCN was 0.358")     # red flag
   nri <- 0.358  # literal, not from CSV
   ```

2. **Table cells without source CSV**

3. **Conflicting cached CSVs** — flag INTERNAL_INCONSISTENCY when different CSVs hold different values for nominally the same statistic

4. **Per-cohort labelling drift** — prose says "Wales" but only `cohort=='SouthWales'` in cached CSVs has that value → LABEL_DRIFT

## Report format
Write `qc_output/agent_5_provenance_report.json` with findings (claim, location, provenance type, source_script, current_pipeline_value, implication, recommended_fix). Plus `agent_5_provenance.md`.

## Decision rule
- **PASS** if every claim traces cleanly to a CSV row OR explicit live computation
- **DRIFT** if some hard-coded literals match the current pipeline (slightly outdated but defensible)
- **FAIL** if any hard-coded literal does NOT match current pipeline — TUDOR-style bugs

## Known traps in this user's history (always check)
1. `TUDOR_write_manuscript.R` line 211: "NRI of TUDOR over DLCN was 0.358"
2. Line 554: "IDI = 0.039"
3. UKB calibration slope
4. UKB Brier score
5. Wales DLCN AUC (often confused with matched-subset AUC)
6. UKB DLCN AUC
7. "estimated DLCN" (eDLCN) vs raw DLCN

You report all options; do NOT decide which value is "correct" or change any file.

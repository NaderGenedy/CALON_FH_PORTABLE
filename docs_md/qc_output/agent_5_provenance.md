# Agent 5 — Provenance Detector — TUDOR Manuscript v5 (2026-05-12)

**Decision: FAIL** — every primary metric in the manuscript is a hardcoded string literal that does not match the current live pipeline.

## Headline numbers (independently re-traced)

| # | Claim (manuscript prose) | Manuscript value | Live pipeline value | Δ | Status |
|---|---|---|---|---|---|
| F1 | NRI of TUDOR over DLCN (Wales) | 0.358 | -0.0969 | -0.4549 | HARDCODED_DRIFT (sign flip) |
| F2 | IDI Wales | +0.039 | 0.0014 | -0.0376 | HARDCODED_DRIFT (28x overstatement) |
| F3 | Calibration slope (UKB) | 6.33 (5.93–6.73) | 1.2312 | -5.099 | HARDCODED_DRIFT (narrative inversion) |
| F4 | Brier score (UKB) | 0.069 | 0.0458 | -0.0232 | HARDCODED_DRIFT |
| F5 | DLCN scoring AUC (Wales) | 0.791 | 0.6896 | -0.1014 | HARDCODED_DRIFT + subset confusion |
| F6 | eDLCN AUC (UKB) | 0.636 | 0.7126 | +0.0766 | HARDCODED_DRIFT |
| F7 | TUDOR AUC (Wales) | 0.842 | 0.8335 (SouthWales) / 0.78 (Wales) | -0.0085 / -0.05 | HARDCODED_DRIFT + LABEL_DRIFT |
| F8 | TUDOR AUC (UKB) | 0.750 | 0.7532 (full) / 0.7684 (LC) / 0.7873 (lancet) | +0.003 | HARDCODED_DRIFT + INTERNAL_INCONSISTENCY |
| F9 | LDLR AUC (Wales) | 0.839 | 0.7646 | -0.0744 | HARDCODED_DRIFT |
| F10 | APOB AUC (Wales) | 0.841 | 0.7522 | -0.0888 | HARDCODED_DRIFT |

## How the drift entered the manuscript

Both `TUDOR_write_manuscript.R` (line 211) and `TUDOR_marked_manuscript.R` (line 201) contain the literal string

```
"Net Reclassification Improvement of TUDOR over DLCN was 0.358, indicating that more than one-third of patients were correctly reclassified."
```

The number 0.358 is hand-typed inside the `cat(...)` / `rd(doc, ...)` call. The script does **not** read any CSV to populate it. The same pattern repeats for every metric in the table above — see the corresponding line numbers in `agent_5_provenance_report.json`.

`TUDOR_LIVE_LEDGER.csv` in the project root is the internal drift ledger: every one of these ten findings is **already flagged** there as `DRIFT`, with the live pipeline value computed and stored. The ledger is being maintained in parallel to the manuscript rather than as a hard pre-submission gate, which is how the drift survived to v5.

## Label drift (prose says "Wales", number is "SouthWales")

- **L1 — TUDOR AUC 0.842 "in Wales"**: The literal 0.842 traces (per `TUDOR_LIVE_LEDGER.csv` row 3) to `loco_predictions_complete.csv [cohort=SouthWales]`. Wales = the larger all-Wales FH Registry. Wales-as-test AUC in `head_to_head_discrimination.csv` is 0.7814 (ENET) or 0.7882 (RF), not 0.842.
- **L2 — Wales DLCN AUC 0.791 vs 0.6896**: Same cohort label, two different denominators. 0.791 is the unmatched/full-Wales AUC (with DLCN imputed to zero for unscorable patients); 0.6896 is the matched/scorable-only AUC. The manuscript quotes 0.791 (the optimistic one) without disclosing this.

## Internal CSV inconsistencies (same statistic, different CSVs disagree)

- **UKB calibration slope**: 1.23 (ledger) vs 7.07 (lipid-clinic CSV) vs 12.40 (lancet summary) vs 6.33 (manuscript).
- **UKB Brier**: 0.046 vs 0.013 vs 0.0085 vs 0.069 (manuscript).
- **UKB TUDOR AUC**: 0.7524 vs 0.7532 vs 0.7684 vs 0.7873 vs 0.750 (manuscript).

Each spread reflects an unresolved "which UKB subset is canonical?" question (full UKB vs lipid-clinic-mimicking subset vs index-only). The manuscript prose talks about "the UK Biobank lipid clinic-mimicking cohort (n = 58,021)", but the canonical lipid-clinic computation in the ledger gives n = 57,965 and AUC = 0.7684, neither matching the manuscript.

## Summary

- **Total claims extracted**: 190 (manuscript_claim_extractor.py)
- **Traced cleanly to a CSV**: 0 of the 10 audited primary metrics
- **Hardcoded literals (DRIFT)**: 10
- **Orphan (nowhere in the project)**: 0
- **Label drift**: 2
- **Internal CSV inconsistency**: 3 metrics

## Recommended fix

1. Refactor `TUDOR_write_manuscript.R` and `TUDOR_marked_manuscript.R` so every numerical literal in `cat()/rd()/sprintf()` calls reads from `tudor_loco_output/*.csv` at runtime. Add `stopifnot(!is.na(x))` after each read.
2. Lock cohort-subset definitions (Wales vs SouthWales vs PASS; UKB-full vs UKB-lipid-clinic). Pick one of each, label consistently in prose.
3. Promote `TUDOR_LIVE_LEDGER.csv` from passive ledger to active pre-submission gate: if any row has status=DRIFT, block `make manuscript`.
4. Re-export the `.docx` from refactored scripts. Re-run all five agents.

## Files written

- `qc_output/agent_5_provenance_report.json` (structured, 10 findings)
- `qc_output/agent_5_provenance.md` (this file)
- `qc_output/agent_5_claims.json` (raw 190 claims from manuscript_claim_extractor.py)
- `qc_output/agent_5_manuscript_text.txt` (extracted manuscript prose)

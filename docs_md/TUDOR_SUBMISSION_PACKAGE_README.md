# TUDOR R1 → R2 — Submission Package Index

**Manuscript:** JCLINLIPID-D-25-01142 R2

**Generated:** 2026-05-12 10:09:49

**Pipeline:** TUDOR_R1_to_R2_FULL_v2.py


## What to upload to Elsevier's Editorial Manager


| File | Purpose |
|---|---|
| `TUDOR_Manuscript_v5_clean.docx` | CLEAN MANUSCRIPT — upload as the main manuscript file. ✅ |
| `TUDOR_Highlights_v5.docx` | HIGHLIGHTS (4 bullets ≤85 chars) — upload to Highlights field. ✅ |
| `TUDOR_disclosure_block.docx` | DISCLOSURE BLOCK (DoI + CRediT + AI + Ethical) — append to manuscript before References. ✅ |
| `TUDOR_AI_use_statement.docx` | STAND-ALONE AI-USE STATEMENT — duplicate of the AI section in the disclosure block, supplied separately if portal requires. ✅ |
| `TUDOR_ResponseToReviewers_R2.docx` | RESPONSE TO REVIEWER #2 — upload to Cover Letter / Response field. ✅ |
| `TUDOR_v5_TRACEABILITY.md` | TRACEABILITY AUDIT — supplementary file documenting which numerical claims were live-re-derived from raw data in this revision cycle. ✅ |
| `TUDOR_v4_to_v5_DIFF_LOG.md` | INTERNAL CHANGE LOG — every additive insertion documented per Reviewer #2 item. ✅ |
| `TUDOR_v4_v5_comparison.html` | ANNOTATED COMPARISON HTML — editor-friendly side-by-side view of insertions. ✅ |
| `TUDOR_R1_to_R2_FULL_v2.py` | REPRODUCIBLE PIPELINE — single Python file that regenerates every deliverable above. ✅ |

## Hard constraints honoured


1. **No numerical content was altered.** Pre-revision claim signatures = post-revision claim signatures (run `python TUDOR_R1_to_R2_FULL_v2.py P3` to re-verify).
2. **Additive only.** Every change is an insertion of a NEW paragraph or a run-level italic toggle on a gene symbol. No existing sentence was deleted or rewritten.
3. **Reviewer #2 items addressed.** R2.1 through R2.6, plus the six Elsevier formatting items.

## Traceability summary (Phase 0 output)


Of the 23 numerical claims in the manuscript ledger:
- **13 LIVE-DERIVED** from raw All-Wales-FH data in this session — Phase 0 reproducer + assertion check confirms match within 1e-3 tolerance.
- **6 RAP-REBLOCKED** (UK Biobank arm Sens/Spec/Brier/calibration slope ± CI) — cannot be live-re-derived this session; documented in TUDOR_v5_TRACEABILITY.md.
- **2 NEEDS-LANCET-COMPLETE** (Wales Index-trained → Cascade-validated AUC 0.842, DLCN 0.791) — produced by TUDOR_LANCET_COMPLETE.R which requires RAP-derived UKB inputs; documented openly.
- **2 HARDCODED-PROSE** (Wales NRI 0.358, IDI 0.039) — hand-transferred from analyst computation in TUDOR_write_manuscript.R; need full pipeline re-run for automated re-derivation.

See `TUDOR_v5_TRACEABILITY.md` for the full provenance table.

## How to verify before submission


```bash
# 1. Re-run the full pipeline from scratch
python TUDOR_R1_to_R2_FULL_v2.py

# 2. Or verify integrity only (fast, no R required)
python TUDOR_R1_to_R2_FULL_v2.py P3

# 3. Generate a Word-native track-changes file (manual, ~30 sec)
# Open TUDOR_Manuscript_v4.docx in Word →
#   Review → Compare → Compare Documents →
#   Original = v4, Revised = v5_clean →
#   Save as TUDOR_Manuscript_v5_TRACKED.docx
```

## After submission


When RAP access is restored, re-run `TUDOR_LANCET_COMPLETE.R` to close the four RAP-REBLOCKED items in the traceability ledger.  If genuine drift is found (it should not be, but the audit must be performed), file a journal correction note through the normal channels.

— Dr Nader Genedy, Cardiology Specialist Registrar, University Hospital of Wales, Cardiff.
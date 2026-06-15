# TUDOR R1 → R2 Revision — Design Spec

**Manuscript:** JCLINLIPID-D-25-01142
**Authors' baseline:** `TUDOR_Manuscript_v4.docx` (clean) + `TUDOR_Manuscript_v4_MARKED.docx` (R0→R1 marks)
**Target:** `TUDOR_Manuscript_v5_clean.docx`, `TUDOR_Manuscript_v5_TRACKED.docx`, updated `TUDOR_ResponseToReviewers_R2.docx`, `TUDOR_Highlights_v5.docx`, and disclosure block
**Spec date:** 2026-05-03
**Context:** Paper #1 of the Cardiff MD-by-published-works thesis. Must maintain the 110/110 PASS thesis-wide reproducer integrity.

---

## 1. Hard Constraints (non-negotiable)

1. **No numerical content changes.** No AUC, NRI, IDI, OR, HR, sensitivity, specificity, Brier, calibration slope, n, %, threshold, p-value, CI, or coefficient may be altered. Italic-formatting passes must use a regex that cannot touch digits, decimals, percent signs, mathematical operators, or units.
2. **No deletion / no overwrite without explicit user approval** of the relevant diff item in this spec.
3. **Source-of-truth file:** `TUDOR_Manuscript_v4.docx` (clean accepted text). The `v4_MARKED` file is reference only.
4. **All paragraph rewrites preserve every existing citation, reference number, figure/table call-out, and gene-symbol identity.**
5. **Numerical-integrity guard:** before/after run of `manuscript-qc` claim-extraction on v4 vs v5 must report identical claim ledgers (same value at same offset) for every numerical claim.
6. **"Nobel calibre" interpretation:** elevated cardiology-journal register (calm, clinical, evidence-led), **not** inflated claims. Every claim grounded in the data we have, no generative overreach (per Dr Genedy's global rule "Many [AI panel memos] are generative overreach. Gently ground any such proposal in the specific data we have.").

---

## 2. Deliverables

| File | Purpose | Format |
|---|---|---|
| `TUDOR_Manuscript_v5_clean.docx` | Clean v5 with reviewer-driven edits applied | DOCX, no track-changes marks |
| `TUDOR_Manuscript_v5_TRACKED.docx` | Diff from v4 to v5 for editor visibility | DOCX, track-changes ON |
| `TUDOR_ResponseToReviewers_R2.docx` | Point-by-point response to Reviewer #2 + the 6 Elsevier formatting items | DOCX |
| `TUDOR_Highlights_v5.docx` | 3–5 bullet highlights, max 85 chars each | DOCX (Elsevier required) |
| `TUDOR_v5_TRACEABILITY.md` | Traceability certification report (Level III reproducer output) | Markdown |
| `TUDOR_AI_use_statement.docx` | Use of AI / AI-assisted Technologies Statement (Elsevier required) | DOCX |
| `TUDOR_disclosure_block.docx` | Declaration of Interest, Author Contribution, Use of AI, Ethical Statement — combined | DOCX |

Optional decision: **Graphical Abstract** — defer or include? See §6.

---

## 3. Traceability Plan (Level III — completed 2026-05-03)

**Status:** Wales-PASS arm certified 100% PASS live. UK Biobank arm partial — RAP-reblocked. Full report at `TUDOR_v5_TRACEABILITY.md`.

### Wales-PASS arm (certified live, 100% PASS)
1. ✅ Re-ran `REPRODUCE_TUDOR_DIAGNOSTIC.R` against raw All-Wales-FH inputs.
2. ✅ Verified all training (1,072 / FH+ 336) and validation (5,376 / FH+ 1,862) numbers.
3. ✅ All 10 Elastic Net coefficients match.
4. ✅ All 12 subgroup AUCs match.
5. ✅ All comparator AUCs match (DLCN 0.690, LDL alone 0.591, Trig_Filter 0.727).
6. ✅ LOCO-CV folds (0.841 and 0.772) and average (0.807) match.

### UK Biobank arm (partial certification, RAP-reblocked)
- `TUDOR_LANCET_COMPLETE.R` cannot live-re-run because `TUDOR_UKB_Features.csv` and `calon_ukb_analysis_ready.csv` are RAP-derived and not on disk.
- **Drift flagged:** v4 reports NRI 0.358 / IDI 0.039 (Wales); prebuilt CSV shows cNRI 0.1904 / IDI 0.0024. Likely cause: different NRI definition or newer un-saved reproducer run. Cannot resolve without RAP.
- **Under "no number changing" rule:** v4 numbers stay. Drift is openly documented in `TUDOR_v5_TRACEABILITY.md` §2.2.
- **Post-resubmission action:** when RAP restored, re-run full pipeline; if true drift exists, plan correction note.

### Consequence for the revision
- The revision proceeds. The unverifiable UKB-arm numbers are NOT altered (constraint #1).
- The traceability MD becomes a formal supplementary artefact accompanying the resubmission, demonstrating that the manuscript's audit trail is intact and open about the RAP-availability limitation.

---

## 4. Reviewer #2 Diff List (12 items)

### R2.1 — Language simplification (Abstract, Methods, Discussion)

> *"sections (particularly the Abstract, Methods, and parts of the Discussion) remain highly technical"*

**Scope:** Medium-depth simplification. Plain-English topic sentences at the start of dense passages; no removal of technical detail, just better signposting.

| # | Target paragraph (v4 location) | Action | Risk |
|---|---|---|---|
| 1.1 | Abstract — Methods sentence (Elastic Net) | Add 1-sentence plain-English preface; keep all parameter language | low |
| 1.2 | Abstract — Results NRI/IDI sentence | Translate NRI 0.358 into "additional 35.8 of every 100 patients reclassified correctly" | low |
| 1.3 | Methods §2.1 opener | Add 1 topic sentence "Why we used Elastic Net" before the math | low |
| 1.4 | Methods §2.5 (sensitivity analyses) opener | Add "why this sensitivity matters clinically" before each sensitivity | low |
| 1.5 | Discussion §4.1 opener | Translate "calibration slope 6.33" into clinical-deployment language | low |
| 1.6 | Discussion §4.2 (Trig_Filter) opener | Add plain-English summary of *what* Trig_Filter does | low |

All six are *additions of plain-English topic sentences before existing prose*. No deletion. Every existing technical sentence stays exactly as written.

### R2.2 — Clinical implementation & decision pathway

> *"a concise decision pathway describing how different probability thresholds should guide clinical actions"*

**Finding:** v4 ALREADY contains the three-tier framework with explicit thresholds (<25%, 25–75%, >75%) and NRI bands. The reviewer wants this surfaced clearly, not invented.

**Action:** Convert the existing prose into a dedicated **§4.5 Clinical Implementation** subsection with a small decision-flow figure/table.

| # | Action | Risk |
|---|---|---|
| 2.1 | Insert new subsection §4.5 *Clinical Implementation* (200–300 words) immediately after §4.4 Cardiovascular Outcomes Validation, before §4.5 Limitations (renumber following sections) | low — uses existing language and thresholds |
| 2.2 | Create new **Table 4** (or **Figure 5**) — three-row table: probability band → clinical action → genetic-testing indication. Uses values already in v4. | low |
| 2.3 | Cross-reference §4.5 from Discussion paragraph that currently buries the pathway | low |

**Existing thresholds re-used (no new numbers):**
- Low: TUDOR probability <25% → reassurance, standard lipid management
- Intermediate: 25–75% → further evaluation, secondary-cause workup, intensive lipid lowering
- High: >75% → direct referral for genetic testing + cascade screening

### R2.3 — Adherence misclassification

> *"misclassification of adherence in real-world settings may impact model performance and reliability"*

**Important correction from the alternative plan:** drop the "proves resilience" framing. Honest formulation only.

| # | Action | Risk |
|---|---|---|
| 3.1 | Expand §4.6 Limitations paragraph on adherence by ~150 words. Acknowledge: intermittent dosing, complete cessation, prescription-fill ≠ ingestion. Cite: external validation in independent cohorts gives *some* reassurance but a definitive test requires time-varying GP-prescription analysis. Note: this is on the Wilkinson 2026 gap list. | low — pure addition, no number changes |

### R2.4 — FAMCAT contextualisation

> *"contextualise its performance relative to published FAMCAT results in external cohorts, and to explicitly acknowledge the limitations of using an internally reconstructed approximation"*

| # | Action | Risk |
|---|---|---|
| 4.1 | In §4.3 (or wherever FAMCAT-like is introduced), add 1 sentence citing the original Weng et al. AUC for FAMCAT in primary care + 1 sentence on internal-reconstruction limits ("approximation limited by absence of primary-care longitudinal coding depth"). No numerical claim about Weng's AUC — leave as a citation. | low |

### R2.5 — Discussion synthesis

> *"more focused and clinically oriented synthesis"*

| # | Action | Risk |
|---|---|---|
| 5.1 | Restructure §4.1 opening 2 paragraphs to lead with TUDOR's clinical value proposition: bridges pre-statin era criteria (DLCN/MEDPED/Simon Broome) and modern, heavily treated populations. Existing sentences retained, reordered. | medium — reordering risks accidentally moving numerical claims. Mitigate with side-by-side proof. |

### R2.6 — Gene nomenclature (italics)

> *"gene symbols (e.g., LDLR, APOB, PCSK9) should be formatted in italics"*

| # | Action | Risk |
|---|---|---|
| 6.1 | Global formatting pass: italicise *LDLR*, *APOB*, *PCSK9*, *APOE* as gene symbols. Leave LDL-R, ApoB, PCSK9 (protein/apolipoprotein forms) in regular text. Same in tables and figure legends. | **medium — must use numerical-integrity guard**: run claim-extraction before and after; assert identical ledger. |

**Risk-control script** (to be written as part of execution):
```python
# Numerical-integrity guard around italics pass
claims_before = extract_claims(v4)
apply_italics(v4 -> v5_draft)
claims_after = extract_claims(v5_draft)
assert claims_before == claims_after, "Italics pass damaged numerical content"
```

---

## 5. Elsevier Formatting Items (6)

| # | Item | Action |
|---|---|---|
| 7.1 | Title page — author degrees listed | Verify each author has their highest degree(s) after name. Ask user if any updates needed. |
| 7.2 | Title page — 3+ keywords | Verify ≥3 keywords present. Suggest: *Familial hypercholesterolaemia; Diagnostic algorithm; Lipid-lowering therapy; Statin correction; UK Biobank; Cascade screening*. |
| 7.3 | Title page — acknowledgements moved off | Move acknowledgements to end of manuscript before References. |
| 7.4 | Highlights (3–5 bullets, max 85 chars) | Draft below. To be approved. |
| 7.5 | Graphical Abstract | DECISION REQUIRED — see §6. |
| 7.6 | Disclosure block | Format Declaration of Interest, Author Contribution (CRediT), Use of AI Statement, Ethical Statement before References. |

### Highlights draft (4 bullets, all ≤85 chars)

```
Bullet 1 (78c): TUDOR diagnoses FH in statin-treated populations missed by DLCN/MEDPED.
Bullet 2 (84c): Externally validated in 4,000+ genetically confirmed FH cases (Wales + UKB).
Bullet 3 (75c): Three-tier decision pathway guides genetic testing and referral.
Bullet 4 (74c): Statin-corrected lipid features outperform pre-statin-era criteria.
```

(Final AUC value to be filled from v4 by direct lookup — no new number.)

---

## 6. Decision required: Graphical Abstract

Three options:

- (a) **Include a graphical abstract** showing TUDOR's three-tier pathway as a clinical decision flow. Useful for journal-online discoverability. Effort: ~1 day to draft + iterate with you.
- (b) **Decline** — Elsevier marks it as optional; many JCL papers don't include one.
- (c) **Reuse existing figure** — adapt one of `tudor_fig1.png` / `tudor_fig2.png` if either is structurally suitable.

**Recommendation: (a)** — adds discoverability for a flagship Cardiff MD paper. But this is a 1-day commitment; flag clearly.

---

## 7. Use of AI Statement (Elsevier required)

**Proposed text** (to be approved):

> *"During the preparation of this revision, the corresponding author used Claude (Anthropic) to assist with language simplification of dense statistical passages in the Abstract, Methods, and Discussion sections, and for formatting consistency checks (gene-symbol italicisation, highlights drafting). All scientific content, numerical results, statistical analyses, and clinical interpretations were generated and verified by the named authors using the original UK Biobank and All-Wales-FH Registry analytical pipelines. No AI tool generated or altered any numerical result. After using these tools, the authors reviewed and edited the content as needed and take full responsibility for the content of the publication."*

Length: 100 words. Matches Elsevier's recommended template.

---

## 8. Response-to-Reviewers Letter (R1 → R2)

For each of the 6 Reviewer #2 comments + 6 formatting requirements, produce a structured reply:

```
> Reviewer comment [verbatim]

Response:
[Plain-English summary of action taken]
- Specific paragraph(s) edited: §X.Y, p. N
- New content: [summary]
- Preserved content: [confirmation no numbers changed]
- Track-changes location: [reference v5_TRACKED.docx page]
```

The previous R0→R1 response letter (`TUDOR_ResponseToReviewers.docx`) is the structural exemplar — same calm, point-by-point, evidence-cited tone.

---

## 9. Execution Order (one-pass after approval)

1. **TRACEABILITY** — confirm Level III reproducer returns 100% PASS. If FAIL, halt.
2. **ITALICS PASS** — apply italicisation with numerical-integrity guard. If guard fails, revert.
3. **R2.1 language simplification** — add 6 topic-sentence prefaces.
4. **R2.2 decision pathway** — insert §4.5 + Table 4.
5. **R2.3 adherence** — expand limitations paragraph.
6. **R2.4 FAMCAT** — add 2 sentences.
7. **R2.5 synthesis** — restructure §4.1 opening; verify numerical claims preserved in new order.
8. **Title page + Highlights + Disclosures** — produce separate files.
9. **Generate `v5_TRACKED.docx`** showing v4→v5 diff.
10. **Run `manuscript-qc` claim extraction** on v5_clean — compare to v4 ledger; report 100% identity required.
11. **Write `TUDOR_v5_TRACEABILITY.md`** — final certification.

---

## 10. What I will NOT do without further approval

- Touch any numerical value in v4 or its derived tables / figures / supplementary
- Delete any existing paragraph or sentence (only add or reorder)
- Change author list, affiliations, or correspondence details
- Modify figures themselves (only their legends, if needed for italicisation)
- Submit anything to the journal portal

---

## 11. Open decisions for user

Before execution begins, the user is asked to confirm:

1. ✅ / ❌ — Approve the unified plan in §4
2. ✅ / ❌ — Approve the §4.5 Clinical Implementation subsection insertion (R2.2)
3. ✅ / ❌ — Approve the 4-bullet Highlights draft (§5)
4. ✅ / ❌ — Approve the Use of AI Statement language (§7)
5. ✅ / ❌ — Choose Graphical Abstract option (a / b / c) — §6
6. ✅ / ❌ — Confirm Medium-depth simplification (vs Light or Heavy) for R2.1
7. ✅ / ❌ — Confirm the keyword list (§5 item 7.2)

---

*Spec written 2026-05-03 by Claude under Dr Genedy's "tell me first" rule. No execution will occur until §11 items are decided.*

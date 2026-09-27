# Examiner's Report on an MD Thesis

**Thesis examined:** *Atherogenic Particle-Burden Phenotyping in Familial Hypercholesterolaemia: Measurement, Ascertainment, Identification and Prediction of Atherosclerotic Cardiovascular Risk.* File MD_SUBMISSION_v3_FINAL.docx, eight chapters, Appendices A to G, 100 references, about 67,700 words.

**Candidate:** Dr Nader Genedy. **Degree:** Doctor of Medicine, Cardiff University.

**Examiner's role assumed:** external examiner in lipid medicine, biostatistics and genetic medicine, applying the doctoral standard: a substantial and original contribution to knowledge, sound in method, defensible in interpretation, and reported with integrity.

**Basis of this report.** The full text of the thesis was read paragraph by paragraph. All 100 references were checked for existence against PubMed (91 resolved; 9 are the candidate's own abstracts and manuscripts, NICE CG71 and one paper that did not resolve in the tools used). Eighteen external numbers quoted in the thesis were checked against the source abstracts or full texts. Internal consistency was checked between chapter text, tables and the claim ledger. The underlying data, analysis scripts, and the two source manuscripts (Paper 14 R4 and CALON r7) were not available to the examiner; judgements about them are text-based and are labelled as such. Date of examination: 27 September 2026.

---

## 1. Recommendation

**Pass, subject to corrections, without re-examination.** The corrections are listed in Section 8. They are numerous, and several go beyond typography, but none requires new data collection, none alters the thesis's principal conclusions, and all can be completed and checked by the internal examiner within the University's corrections period. Should the University's regulations classify a list of this length and character as major rather than minor corrections, my recommendation is unchanged in substance: the thesis should be approved once the listed corrections are made, and a second oral examination is not needed.

**The single reason that drives the recommendation.** The thesis asks one clear question, when the same LDL-C value means the same thing in familial hypercholesterolaemia, and answers it with five measured conditions in two genotyped cohorts under an explicit evidence rule. It is unusually candid about its own limits. Its most direct finding, that three in four registry carriers with a severe untreated LDL-C would not have been graded severe on the value in their record (Section 7.7), is clinically legible and needs no model to stand. What prevents an unconditional pass is not the science but the housekeeping of a large assembled programme: internal inconsistencies between chapters and ledgers, a published description of TUDOR that the evaluated object does not match and that has not yet been corrected in the literature, some inferential wording that overstates what was tested, and provenance gaps that the thesis discloses but has not closed.

---

## 2. Marking

The University awards the MD without a numerical mark. The scores below are the examiner's working assessment on a ten-point scale for each criterion, given so the candidate can see where the strength and the exposure lie. Weights reflect the doctoral standard.

| Criterion | Weight | Score /10 | Basis |
|---|---|---|---|
| Originality and contribution | 20% | 7 | Integrative and quantitative contribution, not conceptual priority; the thesis says so itself (Section 1.6). The under-grading result (7.7) and the bidirectional transport design (Chapter 7) are the strongest increments. Each component principle is established elsewhere. |
| Command of the literature | 15% | 7 | Wide, current and critically read; correctly identifies the nearest prior work for each contribution (Table 2.4). Marked down for an incomplete search (no Embase, AI-assisted unduplicated screening) and four citation-wording inaccuracies found on checking. |
| Research design | 15% | 7 | Estimands declared per task (Table 3.1); frozen application separated from adaptation; competing risks handled. Marked down because neither prognostic model has an untouched evaluation cohort and because the framework and falsifiers were written after the results. |
| Statistical analysis | 15% | 6 | Appropriate methods throughout; family-aware inference in Chapters 5 and 7; whole-procedure bootstrap; paired comparisons on common rows. Marked down for the 28-event registry fit, absent family-aware inference for TUDOR, an incorrect missing-at-random inference, no interval for the key increment over age and sex, and an unverified competing-risk implementation. |
| Interpretation and claims discipline | 15% | 8 | The thesis's best quality. Non-claims are stated in every chapter; retired estimates are kept visible; "no statistically detectable difference" is never read as equivalence. Marked down for a handful of sentences that outrun the test ("relative effects travel", "errors predictable in direction", "the laboratory value was accurate", falsifier status in Table 8.1). |
| Presentation and structure | 5% | 7 | Clear architecture, consistent terminology, good tables. Long; Section 8.2 is mechanical; navigation lists must be rebuilt after corrections; the title over-weights one condition. |
| Research integrity, governance and reproducibility | 10% | 6 | Exemplary disclosure: evidence tiers, analysis-package identifiers, version histories, a zero-silent-correction ledger. Marked down because disclosure has not yet been followed by correction of the published record (TUDOR), because several values rest on retained aggregates that cannot be re-verified, and because one thesis-original result has no located producing output. |
| Publications and impact | 5% | 8 | Two version-of-record papers, one in press, two under revision or review, all in the field's journals; erratum handled correctly. |
| **Weighted total** | | **6.9** | **Pass with corrections. A thesis of clear doctoral standard whose written record needs to be brought to the level of its own stated rules.** |

---

## 3. Strengths

The strengths are real and should be recorded as fully as the weaknesses.

1. **One question, honestly pursued.** The thesis is organised around a single interpretive question and five conditions (Table 1.1), and it does not pretend that the conditions form a causal chain (Section 1.4). Few assembled programmes achieve this coherence.

2. **Claims discipline of a high order.** Every chapter states what it does not claim (Table 1.1, Sections 3.2, 4.10, 5.7, 6.9, 7.15). The retired hazard ratio of 38.55 is kept visible rather than deleted (Section 4.6). The imprecise within-pedigree rung is called imprecise, with both variance conventions shown (Section 5.3). TUDOR's frozen fall to 0.669 is reported ahead of its refitted 0.756 (Section 6.5). CALON-5's calibration failure is a headline, not a footnote (Section 7.6).

3. **Evidence governance.** Tiers T1 to T6, thirteen analysis-package identifiers, a claim ledger (Appendix A), a lineage appendix (B), a version history for the ascertainment manuscript (B.10) and a denominator ontology (C) allow an examiner to trace almost every number. This is rare in a clinical MD and it is what made the present examination possible.

4. **Correct handling of the hard statistical problems.** Competing death by Aalen-Johansen rather than one minus Kaplan-Meier (Sections 3.5, 4.5, 7.1); cause-specific and subdistribution estimands kept apart; frozen application distinguished from recalibration, updating and refitting (Section 3.8); pedigree-clustered errors, wild cluster bootstrap and pedigree fixed effects in Chapter 5; family-grouped paired resampling in Chapter 7; whole-procedure bootstrap optimism for CALON-5; a twenty-event floor and a precision statement for every comparator comparison (Sections 7.8, 7.11).

5. **The under-grading result (Section 7.7).** In 261 genotype-positive patients with measured untreated and on-treatment LDL-C, 145 of 193 above 4.9 mmol/L untreated were not above it on the recorded value, and 66.3 per cent fell a stratum. It needs no reconstruction and no biobank. It is the finding a lipid clinic will remember.

6. **The bidirectional transport design (Chapter 7).** Applying two models in both directions between differently ascertained cohorts, on common rows, with ranking and calibration reported separately and calculability reported as a result (Table 7.3), is a genuine methodological step beyond the one-cohort, one-direction validations in the field (McKay 2022; Mansilla-Rodriguez 2025; Tamehri Zadeh 2025, 2026).

7. **The attenuation ladder (Chapter 5).** Reading the sequence +1.09, +0.81, +0.65, +0.39 as three different operations, adjustment, restriction and change of comparison, rather than as a single decay, is exactly right, and the robustness panel (Table 5.2) is thorough.

8. **Fairness to competitors.** The thesis volunteers the counter-evidence: Akyea 2026, where fully scored DLCN beat an electronic algorithm (Section 6.4); Trinder 2024's polygenic alternative to the route contrast (Section 5.7); the published scores' own miscalibration reported alongside CALON-5's (Section 7.8).

9. **Reference quality.** Of 100 references, 91 resolved in PubMed and all 18 spot-checked numbers were correct or correct with a wording nuance. No fabricated or misattributed reference was found.

---

## 4. Weaknesses, ranked, each with a fix

Severity: **Critical** would, uncorrected, be incompatible with award. **Major** must be corrected before award. **Minor** should be corrected and will be checked.

### Critical

None. No finding in this examination is incompatible with award once corrected.

### Major

**M1. The published TUDOR description does not match the evaluated object, and the record is uncorrected** (Sections 6.2, 6.4, Box 6.1). The paper's title calls the model ascertainment-aware and its abstract says it encodes index versus cascade ascertainment; the recovered eleven-input object has no such term. The thesis discloses this and evaluates only the object, which is the right scientific choice. But disclosure in a thesis does not correct the literature, and the thesis states that no corrigendum has been requested "by the candidate's decision". *Fix:* request a corrigendum from the journal before the corrections are signed off, and cite the request (or the published correction) in Box 6.1. Quote the abstract's actual wording, which does not contain the phrase "Index Effect".

**M2. Internal inconsistencies between chapters and the claim ledger.** Decision-curve counts (Section 7.12: 13 of 30 and 10 of 30; Table A.2: 12 and 9). Comparator tally (Section 7.8: 3 of 4, one resolved; Table A.2: 4 of 4, none resolved). Two-term model lock (Sections 3.4.1, 7.2 and Table 7.1 say released; Table B.1 says not located). Falsifier status for Objective 4 in Table 8.1 reads the opposite of what Table 1.2 defines. A thesis whose governing rule is zero silent correction cannot carry these. *Fix:* reconcile each against manuscript r7 and the model lock, correct the stale side, and add the correction to the change ledger.

**M3. An incorrect missing-data inference, repeated** (Sections 5.2, 7.4, 8.7, Table 8.3). That participants without a recorded diabetes status had about half the hazard (0.40, 0.21 to 0.77) shows that missingness is outcome-associated and that the complete-case set is selected. It does not establish that the data are missing not at random, and "reconstructed rows are not missing at random" misuses the term for a selection pattern. The decision not to impute is defensible; the justification as written is not. *Fix:* reword to "missingness was associated with the outcome; complete-case analysis is therefore a selected analysis; no imputation was performed", and remove "not missing at random" from the limitation register unless a sensitivity analysis supports it.

**M4. Sentences that outrun the test.** "Relative effects travelled between them; absolute risk did not" (Sections 7.6, 8.2.4) rests on directionally agreeing hazard ratios in all-carrier fits, not on a test of coefficient transport. "The error is systematic and predictable in direction" (Section 8.3) is a group-level observation stated as if it held for individuals. "The laboratory value was accurate" (Section 1.1) asserts analytical accuracy that was never tested. *Fix:* restrict each sentence to what was observed; the argument survives the restriction.

**M5. The two prognostic models have no untouched evaluation cohort, and the text sometimes forgets it.** The thesis labels Chapter 7 "cross-cohort evaluation" (correct) but the reverse-direction P/LP fit rests on 28 events (5.6 per parameter, below the prespecified floor), the two-term model was selected after CALON-5's results were known, and no interval exists for CALON-5's increment over age and sex, so the claim that ranking beat demographics cannot be judged. *Fix:* compute the paired interval for CALON-5 over age and sex on registry rows (no UK Biobank access needed); state the two-term model's label as "development-only selected" wherever it appears; keep the 28-event fit explicitly exploratory in every table that uses it.

**M6. Provenance gaps that remain open.** The 71-variant synthesis (THESIS-AN-02, +0.91 mmol/L) has no located producing output yet appears in Table 5.2 and Appendix A. Several UK Biobank values are retained aggregates that cannot be re-verified. The reason for closure of the Research Analysis Platform is not stated. *Fix:* either locate the output or move the 71-variant value to a footnote labelled unverified; state the reason and date of platform closure in the governance statement; keep the retained-aggregate label, which is correctly applied.

**M7. Family dependence in TUDOR inference.** DeLong's test treats relatives as independent; the Welsh head-to-head differences carry no paired interval; no family-aware optimism estimate exists (Sections 6.2, 6.4). *Fix:* a family-grouped bootstrap of the Welsh paired AUC differences is computable from the registry alone and should be added; if it cannot be, state so and downgrade the head-to-head claim to "point estimates without paired inference".

**M8. Citation wording inaccuracies found on checking.** Tybjaerg-Hansen 2005 is quoted for LDL-C increments; the gradient is total cholesterol (Sections 1.3, 5.7, Table 2.4). Jin 2023 adjusted for VLDL, LDL and HDL particle concentrations, not LDL particles alone (Sections 2.3, 4.7). Marston 2022's abstract names non-HDL-C and triglycerides, not LDL-C, in the mutual-adjustment sentence (Section 2.3). Bogsrud 2025's null versus non-null contrast is within 61 FH newborns, not 113 (Section 5.2). The candidate's own Atherosclerosis Plus abstract (reference 82) claims "dual external validation in 4,028 genetically confirmed cases", which the thesis rightly disowns but does not explain. *Fix:* correct each sentence; for reference 82 add one sentence in Section 6.2 explaining what the 4,028 denominator was and why it is not independent validation.

### Minor

**m1.** The Summary omits the confidence interval for the primary hazard ratio; restore 1.070 to 1.718.
**m2.** Reference 52 is listed as "submitted" while the text uses revision R4; update.
**m3.** Objective 2's falsifier is worded differently in Table 1.2 and Section 5.3; harmonise.
**m4.** Chapter 4 does not state the prevalent-ASCVD exclusion rule and count for the 1,461 frame; add one sentence.
**m5.** Chapter 7 does not enumerate the "dominant FH genes" or give per-gene counts; add them.
**m6.** The UK Biobank ancestry restriction is mentioned once (Section 5.7) and never defined; define it and confirm it is the same across chapters.
**m7.** Reference 77 supports two different statements (Sections 2.8 and 8.8); check or split.
**m8.** Three reconstruction-agreement figures (r 0.09; MAE 1.20 and r 0.32; MAE 1.25 and R-squared 0.0975) and a stated 4.5 mmol/L bias should be brought into one table with one sentence on why they differ.
**m9.** TUDOR has narrative risk-of-bias considerations (Table 6.5) but no PROBAST+AI appraisal, unlike the Chapter 7 models (Appendix E.4); complete it.
**m10.** The 362-event figure for a decisive study is an extrapolation from observed differences; label it illustrative and state its assumptions (horizon, incidence, censoring).
**m11.** The undefined decision-curve cells should be inspected to distinguish an estimator failure from mathematically undefined net benefit, and preserved as missing.
**m12.** Section 8.2's seven-step template should be written as connected argument.
**m13.** The glossary conflates model state, development independence, calibration summary and publication status; separate them.
**m14.** The title leads with particle burden, one of five conditions; consider a subtitle naming the question.
**m15.** Navigation lists and cross-references must be rebuilt after corrections.

---

## 5. Chapter-by-chapter assessment

**Front matter.** Summary is accurate and admirably restrained ("No model is ready for clinical use"). Declarations are adequate; the AI statement is brief and should say that title screening in the search used a model first pass. The governance statement should give the reason the Research Analysis Platform closed.

**Chapter 1.** A model introduction: a clinical paradox, one question, five conditions, explicit non-claims, and objectives with falsifiers. The falsifiers' retrospective status is disclosed three times; that is the right handling, but the text should call them interpretive stress tests rather than protection against selection. Score 8.

**Chapter 2.** A claim-led critical review that correctly separates cholesterol mass from particle number, ratio from residual, recognition from exposure, and discrimination from calibration. Table 2.4 is the most useful table in the thesis for an examiner. Weak points: the search is incomplete and single-reviewer; four cited numbers need rewording (M8). Score 7.

**Chapter 3.** The strongest methods chapter I have examined in an MD: estimand registry, identity dictionary, cohort-role matrix, reconstruction conventions, five model states, evidence tiers, an adoption gate. Its honesty about what was not done (diagnostics not recovered, TUDOR family optimism absent) is a strength. Its weakness is that three reconstruction conventions coexist without a unifying rationale, which the chapter concedes. Score 8.

**Chapter 4.** A correct and courageous re-analysis: the published ratio threshold is retired on stated criteria, a continuous residual replaces it, and the frozen Welsh equation transports directionally. The chapter is candid that the increment in discrimination beyond LDL-C is about 0.011 with no interval, that four events are unconfirmed, and that the competing-risk sensitivity is unverified. The frame is operational rare LDLR/APOB variants, not confirmed HeFH, and the chapter never forgets it. Missing: the prevalent-disease exclusion statement (m4), a gene-stratified sensitivity, and a comparison with risk-weighted apoB. Score 7.

**Chapter 5.** A careful ascertainment study with the right estimand language ("selection into observation, not a biological effect"). The attenuation ladder, robustness panel, provenance restrictions and leave-one-out analyses are thorough. The within-pedigree result is honestly imprecise and convention-dependent, and the version history in Appendix B.10 shows it has moved. Two things weaken the chapter: the 71-variant synthesis with no located output (M6), and the MAR wording (M3). Chapter 5's values are those of a manuscript under revision and were not re-run; that dependency is disclosed and acceptable. Score 7.

**Chapter 6.** The comparison is well framed as an electronic, common-data problem and the transport ladder is complete. The chapter is frank that the comparators were input-starved and that Akyea 2026 points the other way. The description-object discrepancy (M1) and the absent family-aware inference (M7) are the chapter's exposures; the reconstruction error (MAE 1.20 mmol/L) is well handled. Calibration is correctly described as unavailable for the transport cohort. Score 6.

**Chapter 7.** Methodologically the most ambitious chapter and the one with the most to correct. The bidirectional design, the calculability result, the paired comparator comparisons with an events-to-resolve column, and the under-grading analysis are all strong. The exposures are the 28-event fit, the post-results selection of the two-term model, the missing interval over age and sex, the R-squared of 0.0975 for the reconstruction equation, the MAR wording, and the ledger inconsistencies. The chapter labels itself developmental and high risk of bias, which is correct. Score 6.

**Chapter 8.** The answer is given part by part with falsifier status, which is good practice, but Table 8.1 misstates Objective 4's status (M2). The sex section is properly bounded. The limitations table is ranked, which examiners appreciate. The resolving studies are concrete. Section 8.2 reads as a template. Score 7.

**Appendices.** Appendix A and B are exemplary in intent; A carries two of the inconsistencies (M2). Appendix D is honest about the three conventions. Appendix F is a model of how to report an incomplete search. Score 7.

---

## 6. Verification log

What was checked, and how.

| Check | Method | Result |
|---|---|---|
| Full text read | Every paragraph of Chapters 1 to 8 and Appendices A, B, D, F, G | Done |
| Reference existence | PubMed ID converter, citation matcher, search, metadata for all 100 references | 91 resolved; 9 not (refs 17, 22, 50, 52, 66, 73, 82, 83, 86) |
| External numbers | 18 quoted figures compared with source abstracts or PMC full text | 13 exact; 5 exact with a wording nuance (M8); 0 discrepancies; 1 unverifiable (ref 82) |
| Internal consistency | Chapter text against Tables A.2, B.1, 1.2, 8.1 | Four inconsistencies (M2) |
| Arithmetic spot checks | 145/193 = 75.1%; 173/261 = 66.3%; 74/5 = 14.8; 28/5 = 5.6; 1,009/3,209 = 31.4%; 945/1,639 = 57.7%; 715/2,398 = 29.8%; 1,904 = 815 + 1,089; 235 = 1,904 − 780 − 889 | All correct |
| Design-claim language | Search for causal and temporal verbs in the Summary and Discussion | Three sentences flagged (M4); otherwise disciplined |
| Data, scripts, locks | Not available to the examiner | Not verified; judgements about analyses are text-based |
| Source manuscripts (Paper 14 R4; CALON r7) | Not available | Not verified; M2 must be resolved against them |

---

## 7. Oral examination

The candidate should expect the oral to concentrate on: the TUDOR description-object discrepancy and why the literature was not corrected; the independence of Chapter 7's cohorts and the 28-event fit; the within-pedigree result and its variance conventions; the missing-data reasoning; the retained aggregates and the platform closure; and what, concretely, a clinic should do differently. A candidate who concedes each point in the first sentence and then states the bounded claim will pass the oral comfortably; the thesis already contains those bounded claims.

---

## 8. Required corrections (to be checked by the internal examiner)

1. Request a corrigendum for the TUDOR paper and record it in Box 6.1 (M1).
2. Reconcile decision-curve counts, comparator tally, two-term lock status and Table 8.1 falsifier status; log each in the change ledger (M2).
3. Reword the missing-data inferences in Sections 5.2, 7.4, 8.7 and Table 8.3 (M3).
4. Restrict the three over-reaching sentences (M4).
5. Add the paired interval for CALON-5 over age and sex on registry rows, or state it cannot be computed and downgrade the claim; label the two-term model "development-only selected" throughout; keep the 28-event fit exploratory in every table (M5).
6. Locate the THESIS-AN-02 output or demote the 71-variant value; state the reason and date of platform closure (M6).
7. Add a family-grouped bootstrap of the Welsh TUDOR paired differences, or downgrade the head-to-head inference (M7).
8. Correct the five citation wordings and explain the 4,028 denominator (M8).
9. Minor items m1 to m15.

Corrections 1 to 8 should be accompanied by a one-page note to the internal examiner listing where each change was made.

---

## 9. Recommended, not required

- Complete the Embase search and duplicate screening, then re-state the novelty sentences.
- Build the participant-overlap and definition matrix across chapters (the kit's K012-013) and replace "common definitions" with the accurate description.
- Add a gene-stratified sensitivity for Chapter 4 (LDLR versus APOB) even if imprecise.
- Consider moving Appendices B, C, E and F to supplementary files if regulations permit, shortening the bound thesis.
- Consider a subtitle that names the question.

---

## 10. Bottom line

This is a thesis of doctoral standard written by a candidate who understands the difference between what data show and what one would like them to show. Its contribution is integrative and quantitative, exactly as it claims, and its most durable finding needs no model at all. The corrections are about bringing the written record up to the thesis's own rules: correct the published description of TUDOR, remove the inconsistencies between chapters and ledgers, say only what was tested, and close or label the provenance gaps. Once that is done, I would be content to see the degree awarded.

*Examiner's report prepared 27 September 2026. Text-based examination; data, scripts and source manuscripts were not available.*

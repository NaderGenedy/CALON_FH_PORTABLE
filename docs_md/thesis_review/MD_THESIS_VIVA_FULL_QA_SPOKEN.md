# MD Thesis Viva: Every Expected Question, Section by Section, with Spoken Answers

**Thesis:** *Atherogenic Particle-Burden Phenotyping in Familial Hypercholesterolaemia: Measurement, Ascertainment, Identification and Prediction of Atherosclerotic Cardiovascular Risk* (MD_SUBMISSION_v3_FINAL.docx). **Candidate:** Dr Nader Genedy, Cardiff University.

**Purpose.** A complete question bank for the annual review or viva, ordered by thesis section, with each answer written as it should be spoken. The examiner is assumed to hold three hats: lipid medicine, biostatistics and genetic medicine. Questions are grouped under the section they arise from; cross-cutting questions follow.

**Sources.** Every number comes from the thesis and carries its section. External papers are cited by the thesis reference number and, where the DOI resolved in PubMed during preparation, the DOI is given as "verified". Nine references did not resolve (the candidate's own abstracts and manuscripts, NICE CG71 and reference 17) and carry no DOI here; quote them only as "the thesis cites". The full register is Appendix 2 of this document.

**Delivery cues used in the answers.**

| Cue | Meaning |
|---|---|
| [pause] | a beat of one second; let the point land |
| [long pause] | two to three seconds; used before a concession or a key number |
| [slow] | slow down; used on numbers and on the sentence that carries the claim |
| [firm] | steady voice, no hedge words, eye contact |
| [concede] | lower the volume slightly; accept the point without defensiveness |
| [lift] | lighter tone; used when moving from limitation to what was done about it |
| [stop] | end the answer; do not add a further sentence |

Speak the numbers slowly. Never fill a pause with "sort of" or "I think". When a concession is due, make it in the first sentence, not the last.

---

## Contents

- Front matter (Summary, Statements)
- Chapter 1 (Sections 1.1 to 1.8)
- Chapter 2 (Sections 2.1 to 2.12)
- Chapter 3 (Sections 3.1 to 3.10)
- Chapter 4 (Sections 4.1 to 4.10)
- Chapter 5 (Sections 5.1 to 5.7)
- Chapter 6 (Sections 6.1 to 6.9)
- Chapter 7 (Sections 7.1 to 7.15)
- Chapter 8 (Sections 8.1 to 8.9)
- Appendices A to G
- Cross-cutting: integrity, causality, translation, the candidate
- Appendix 1: corrections to make before the meeting
- Appendix 2: DOI register (PubMed-verified)

---

## Front matter: Summary and Statements

**F1. Summarise the thesis in two minutes.**

[firm] "It asks one question. [pause] When does the same LDL-C value mean the same thing in familial hypercholesterolaemia, and what must be added when it does not? [pause] I tested five conditions in two linked genotyped cohorts, the All-Wales FH registry and UK Biobank. [slow] Particles: in 1,461 carriers of rare LDLR or APOB variants with 75 events, apoB above the value expected for LDL-C carried a hazard ratio of 1.356 per standard deviation, and a Welsh equation applied unchanged kept the association. [pause] Route into care: at matched LDLR residues, probands had untreated LDL-C 1.09 millimoles per litre higher than cascade-detected relatives, narrowing to 0.39 within families. [pause] Treatment: a treatment-aware triage score, TUDOR, ranked carriers with an AUC of 0.760 against 0.652 for reconstructed criteria, but fell to 0.669 when frozen in UK Biobank; and in the registry, three in four carriers with a severe untreated LDL-C were graded milder on the recorded value. [pause] Cohort: a five-term routine-data risk model kept its ranking in both directions, concordance about 0.70, but its absolute risks failed in opposite directions. [pause] And outside FH, the same apoB cut-point carried different absolute risk in women and men. [long pause] [firm] LDL-C stays indispensable. It should not travel alone. No model here is ready for clinical use. [stop]"

Sections: Summary; 1.1; 1.5; 8.1.

**F2. The Statement of original work says you did the statistics. What did Dr Aubin do?**

"[firm] Statistical supervision and methodological advice. [pause] I designed, ran and verified the analyses in R and Python. Where a chapter's values come from a manuscript rather than a thesis rerun, the thesis says so: Chapter 5 reports Paper 14 R4 values that were not independently re-run, and Chapter 7 reports manuscript r7. [pause] I take responsibility for every number and for the evidence grading. [stop]" (Statement of original work; Sections 5.3, 7.1.)

**F3. Your AI statement is three sentences. Is that enough in 2026?**

"[concede] It is the minimum, and I would expand it if the committee wishes. [pause] What it establishes is the boundary: AI tools helped with literature searching, drafting, editing and consistency checking; they had no access to participant-level data and ran no analysis; I reviewed every passage. [pause] Appendix F adds that title screening used a language-model first pass with single-reviewer adjudication and was not duplicated, and lists that as a limit. [lift] The evidence rule in Section 3.1 is the real protection: model-generated prose is never numerical or bibliographic evidence, and every number must trace to a publication, manuscript, package or script. [stop]" (AI statement; Sections 3.1; Appendix F.8.)

**F4. Why did UK Biobank access cease, and what did you lose?**

"[firm] Access to the Research Analysis Platform closed during the final phase of the programme. [pause] I lost the ability to re-run participant-level UK Biobank analyses. What I kept are frozen aggregate packages; the Chapter 4 primary estimates sit in THESIS-AN-01, frozen on 12 July 2026. [pause] Any UK Biobank value that cannot be traced to a located producing output is labelled 'retained aggregate; not re-verifiable while the platform is closed' and is used only for sensitivity or description. [concede] The reason for closure should be stated in the thesis in one sentence, and I will add it. [stop]" (Data governance statement; Sections 3.3.4; Appendix A.1.)

---

## Chapter 1: Introduction

### Section 1.1 The clinical paradox

**1.1a. Why is FH a "stringent" setting for this question?**

"[firm] Because LDL-C should be at its most informative in a disorder of the LDL receptor. [pause] If particles, route into care, treatment and cohort still change what the number means here, those dependencies cannot be dismissed as peculiarities of mixed dyslipidaemia in the general population. [stop]" (Section 1.1.)

**1.1b. Give me the single example you open with.**

"[slow] The guideline severity cut-point is an untreated LDL-C of 4.9 millimoles per litre. In the registry, among genotyped carriers with a documented untreated value at or above it, applying the same cut-point to the value in their current record would have missed three in four. [pause] The laboratory value was accurate. What it represented had changed. [stop]" (Sections 1.1, 7.7: 145 of 193, 75.1 per cent.)

**1.1c. Isn't this just "treated LDL-C is lower than untreated"? Where is the science?**

"[lift] The science is in measuring how often, by how much, and what else does the same thing. [pause] Treatment is one of five conditions. The route into care changes the phenotype seen by 1.09 millimoles per litre at the same residue. The cohort of development changes a model's absolute risks by a factor of two in opposite directions. Particle concentration changes the event rate at the same LDL-C. [pause] Each is quantified in genotyped cohorts with common definitions, which had not been done together. [stop]" (Sections 1.1, 1.6.)

### Section 1.2 Scientific basis

**1.2a. Why is a current LDL-C not a measure of exposure?**

"[firm] Because it records one point on a lifelong trajectory. [pause] It does not contain the preceding concentrations, treatment interruptions or accumulated arterial exposure. Similar present values can follow early diagnosis with sustained therapy, decades of untreated elevation, or intermittent adherence. [pause] The cumulative-exposure evidence I cite is the EAS consensus and the CARDIA-based analysis of cumulative LDL-C in young adulthood. [stop]" (Section 1.2; refs 1 and 9, both DOI-verified.)

**1.2b. Why is apoB "not an LDL-specific particle assay"?**

"[slow] Because every apoB-containing particle carries one apoB molecule: VLDL, IDL, LDL and lipoprotein(a) carry apoB-100, chylomicrons and their remnants apoB-48. [pause] Routine plasma apoB therefore counts particles across those classes. It approximates particle number; it is not an LDL count. [stop]" (Sections 1.2, 2.3; refs 10 to 12.)

**1.2c. What is "complementarity, not superiority", and why did you settle on it?**

"[firm] The premise that apoB refines LDL-C rather than replacing it. [pause] Population studies show appreciable variation in apoB at the same LDL-C, and risk in discordant groups aligns more closely with apoB. Other studies find similar associations for apoB and non-HDL cholesterol, or conclusions that depend on how correlated lipids are modelled. [pause] Given that split, complementarity is the claim the evidence supports. [stop]" (Section 1.2; refs 13 to 18.)

**1.2d. Why does the thesis prefer a residual to a ratio, in one sentence each?**

"[slow] A ratio is mechanically coupled to its denominator and assumes proportionality across the whole range. [pause] A residual asks how far observed apoB lies from the value expected at that LDL-C, and it declares the equation it depends on, which can then be frozen and tested. [stop]" (Section 1.2; refs 19, 20.)

**1.2e. What does back-calculating an untreated LDL-C actually give you?**

"[firm] A model-based estimate, not a recovered measurement. [pause] Its accuracy depends on drug, dose, combination therapy, adherence, individual response and timing. It can preserve cohort means or rank order while being wrong for an individual. [pause] It may support triage or model development. It does not recover a person's true untreated value. [stop]" (Section 1.2; refs 21 to 23.)

### Section 1.3 Selection into observation

**1.3a. Define the ascertainment-conditioned phenotype without jargon.**

"[slow] The FH phenotype as it appears in people selected by a particular route into care. [pause] It is a property of the selected group, not of the variant. It describes how the data were generated. It does not say that referral changes anyone's biology. [stop]" (Sections 1.2, 1.3; Glossary.)

**1.3b. Who showed this before you?**

"[firm] Tybjaerg-Hansen and colleagues in 2005. Carriers of the same LDLR mutations showed progressively larger cholesterol increments when identified in the general population, among patients with ischaemic heart disease, and in FH clinics. They concluded a mutation's phenotype should be estimated in unselected carriers. [concede] I should note that their gradient is reported in total cholesterol; my text says LDL-C increments and I will correct the wording. [stop]" (Section 1.3; ref 91, DOI verified 10.1161/01.ATV.0000149380.94984.f0.)

**1.3c. Why does selection matter for a diagnostic test's AUC?**

"[firm] Because sensitivity, specificity and AUC are properties of a test in a stated population, not of the score alone. [pause] Ransohoff and Feinstein called it spectrum bias in 1978. A score that separates conspicuous clinic cases from unaffected controls faces a harder problem among treated relatives or genotype-first carriers. [stop]" (Section 1.3; ref 27, DOI verified 10.1056/NEJM197810262991705.)

**1.3d. Are UK Biobank and the Welsh registry interchangeable halves of one population?**

"[firm] No. [pause] UK Biobank offers genotype-first carriers with linked outcomes in a volunteer population that is itself selected. The Welsh service offers specialist, cascade and pedigree information within clinical care. [pause] Disagreement between them can expose spectrum, measurement or transport limits, if the model, outcome, time origin and receiving cohort are declared. Agreement does not establish universal validity. [stop]" (Section 1.3; refs 6 to 8, DOI verified.)

### Section 1.4 The conceptual framework

**1.4a. Is Measure to Select to Identify to Predict a causal pathway?**

"[firm] No. It is the order in which the thesis reads the evidence, and Figure 1.1 says so. [pause] It is an inferential dependency: measurement defines the phenotype; selection determines which part of its distribution reaches analysis; identification asks whether the disorder can be recognised after treatment; prognostic modelling follows once population, starting point and outcome are defined. [pause] Uncertainty introduced early can propagate downstream if it is hidden in a reconstructed value or a model score. [stop]" (Section 1.4.)

**1.4b. When was the framework written?**

"[concede] After the programme, as a thesis-level synthesis. It was not prespecified as one biological theory. [pause] Section 1.4 says that in the first paragraph, and Section 8.1 repeats it when the falsifiers are assessed. [stop]"

**1.4c. Name the five operational checks.**

"[slow] Estimand. Target population. Sampling frame. Measurement model. Evidence for transport. [pause] Every result is meant to carry all five. [stop]" (Section 1.4.)

**1.4d. Why does the fifth condition, sex, sit in Chapter 8 rather than in an empirical chapter?**

"[firm] Because its evidence comes from a population cohort, not an FH cohort. [pause] It is used as a limit on reading LDL-C as a single marker, examined in a bounded section, and explicitly not presented as an FH-specific or causal result. [stop]" (Sections 1.2, 1.4, 8.4.)

### Section 1.5 Question, aim, objectives and falsifiers

**1.5a. Your falsifiers were written after the results. Is that not hypothesising after the results are known?**

"[concede] The falsifiers were written across a completed programme, and the thesis says so in Sections 1.4, 1.5 and 8.1. [pause] They do not alter the prespecification or evidence status of the source studies. [lift] What they do is force my reading of my own results to be checkable, and they set the terms of the resolving studies in Section 8.8, which are prospective. [pause] Where a source had a real prospective register, I report that instead: Paper 14 registered 120 analyses in fifteen families with false-discovery control. [stop]" (Sections 1.5, 5.2.)

**1.5b. Objective 4 was "restated after results for both models were known". Explain.**

"[firm] Objective 4 originally named CALON-C. When CALON-5 replaced it, the objective was rewritten to describe two models, CALON-5 and a two-term model, with neither designated primary, and the thesis states that this happened after results for both were known. [pause] I disclose it rather than backdate it. [stop]" (Section 1.5.)

**1.5c. State each objective's estimand in one phrase.**

"[slow] One: the association between a continuous LDL-conditioned apoB residual and first ASCVD. Two: the difference in untreated LDL-C between probands and cascade-detected relatives at matched residues. Three: discrimination between genotype-confirmed carriers and non-carriers by TUDOR against reconstructed criteria. Four: the cumulative incidence of first MACE with competing death, ranking and calibration separately. Five: absolute risk at the same value and cut-point by sex and stage. [stop]" (Section 1.5, Table 1.2.)

**1.5d. Which falsifier came closest to being met?**

"[firm] Objective 2, within families. [slow] The within-pedigree difference was 0.39 millimoles per litre with an interval from 0.05 to 0.73, but under a dummy-variable degrees-of-freedom correction the interval spans zero, minus 0.06 to 0.85. [pause] Table 8.1 records it as 'not met, narrowly'. [stop]" (Sections 5.3, 8.1.)

**1.5e. What is the "object evaluated" in Objective 3, and why does that phrase appear?**

"[firm] The recovered executable TUDOR score. [pause] It differs from the published description in having no ascertainment-route term. The objective names the object so that no result is attributed to a model that does not exist as an executable. [stop]" (Sections 1.5, 6.2.)

### Section 1.6 Contribution in relation to prior work

**1.6a. What is new, given every component principle is established?**

"[firm] The measurement of the five conditions together in linked genotyped cohorts with common definitions, so that their direction and relative size can be compared and the repair for each named. [pause] Section 1.6 gives the nearest prior work for each and states the increment. [concede] The claim is integrative and quantitative, not one of conceptual priority, and it is conditional on a search that did not include Embase. [stop]" (Section 1.6; Appendix F.)

**1.6b. Nearest prior work for the particle condition?**

"[slow] Johannesen 2024, excess apoB in the Copenhagen General Population Study; Pan 2026, a frozen excess-apoB equation in statin-treated coronary disease validated externally in UK Biobank; and my own discordance paper in Welsh FH. [pause] The increment is a continuous residual in rare-variant carriers, a Welsh equation applied unchanged in UK Biobank, and the retirement of my own ratio threshold. [stop]" (Section 1.6; refs 38, 39, 51, 68, all DOI-verified.)

**1.6c. Nearest prior work for the treatment condition?**

"[slow] Besseling 2017, where statin use predicts an FH-causing variant; FAMCAT, where drug class and potency are inputs; the Welsh service's own scoring criteria, Haralambos 2015; and Akyea 2026, where fully scored DLCN beat FAMCAT in a tertiary clinic. [pause] Chapter 6 adds pretreatment reconstruction with a shared-data comparison and frozen application; Chapter 7 adds the direct measurement of under-grading. [stop]" (Section 1.6; refs 42, 43, 44, 45, 92, 93, all DOI-verified.)

**1.6d. Nearest prior work for the cohort condition?**

"[slow] McKay 2022 transported SAFEHEART-RE into English primary care by imputation; Mansilla-Rodriguez 2025 validated it in Australia; Tamehri Zadeh 2025 and 2026 compared FH scores on common data. [pause] Each ran in one cohort, one direction. Chapter 7 runs in both directions between two differently ascertained UK cohorts, separates ranking from calibration, models competing death, and reports calculability. [stop]" (Section 1.6; refs 28, 29, 94, 95, all DOI-verified.)

### Section 1.7 Scope

**1.7a. What is deliberately out of scope?**

"[firm] Homozygous FH, structural prediction and AlphaFold-derived severity scoring, population-wide apoB discordance as a contribution, and Lp(a) and polygenic burden as anything more than modifiers. [pause] Incidental and privately obtained genetic findings are considered briefly in Section 8.6. [stop]" (Section 1.7.)

**1.7b. Three different genetic definitions are used. Why not one?**

"[concede] Because the programme assembled them at different times under different data extracts. [pause] Chapter 4 uses an operational rare LDLR/APOB variant frame; Chapter 7 uses ClinVar pathogenic or likely-pathogenic single-variant carriers; the CALON-C lineage uses an operational LDLR-carrier flag. [firm] None is relabelled genetically confirmed HeFH, and Box 3.1 keeps them apart. Uniform adjudication is listed as what would remove the limitation. [stop]" (Sections 1.7, 3.2; Table 8.3, L09.)

### Section 1.8 Structure

**1.8a. Why do the chapters "inherit an open problem" from one another?**

"[slow] Chapter 4 finds that particle concentration carries information at the same LDL-C, but measurement cannot be separated from who was measured; so Chapter 5 asks how the route into care shapes what is seen. Route and treatment together obscure the classical phenotype; so Chapter 6 asks whether a treatment-aware score can still find carriers. Finding a carrier does not determine prognosis; so Chapter 7 asks whether a risk model's ranking and absolute risks survive a change of cohort. [stop]" (Section 1.8.)

---

## Chapter 2: Literature review

### Section 2.1 Purpose and stance

**2.1a. This is not a systematic review. Why not, and what did you do instead?**

"[concede] It is a critical, claim-led review, and the run state is recorded as 'search incomplete'. [pause] PubMed, Crossref including EHJ congress abstracts, the Europe PMC preprint index and OpenAlex were searched from 23 September 2016 to 23 September 2026, with one adversarial query per claim of absence. PubMed and Crossref were searched again on 27 September: 27 queries, 339 unique records, 23 retained. [pause] Embase was not searched; Elicit, Consensus, Scite and SciSpace returned no records; screening stopped at the top 20 to 25 records per query and was single-reviewer after a labelled language-model first pass. [firm] Every 'no study identified' sentence is conditional on that. [stop]" (Section 2.1; Appendix F.)

**2.1b. What are the four load-bearing propositions?**

"[slow] One: LDL-C and apoB are correlated but not interchangeable. Two: clinical criteria were built around a visible, commonly untreated phenotype, whereas detection now often occurs after treatment, through cascade or genotype-first routes. Three: treatment and ascertainment determine which version of the phenotype reaches the analyst. Four: prediction models require a declared estimand and evaluation across discrimination, calibration, overall error and utility. [stop]" (Section 2.1, Table 2.1.)

### Section 2.2 FH as exposure and recognition

**2.2a. What are the three category errors you want to prevent?**

"[slow] A low on-treatment LDL-C does not negate previous exposure. Failing a clinical score does not prove absence of a disease-causing variant. A pathogenic variant does not specify present absolute risk. [stop]" (Section 2.2.)

**2.2b. Why call FH both an exposure and a recognition event?**

"[firm] Because the biological exposure begins before recognition and may never be recognised. [pause] Consensus statements place cumulative LDL exposure, not any single result, at the centre of risk. The recognition event, referral, testing, registration, is what makes the person visible to a dataset. [stop]" (Section 2.2; refs 1 to 3, 57, DOI-verified.)

### Section 2.3 LDL-C, apoB and atherogenic burden

**2.3a. Give me the strongest evidence that apoB carries information beyond LDL-C.**

"[slow] In NHANES, at LDL-C of 100 milligrams per decilitre the middle 95 per cent of apoB ran from 66 to 99. [pause] In statin-treated cohorts, apoB retained 1.27 per standard deviation for myocardial infarction when the other measures did not after mutual adjustment. [pause] In Copenhagen, discordantly high apoB with low LDL-C carried a hazard ratio of 1.49 for infarction. [stop]" (Section 2.3; refs 15, 13, 14, DOIs verified: 10.1001/jamacardio.2024.1310; 10.1001/jamacardio.2021.5083; 10.1016/j.jacc.2021.01.027.)

**2.3b. And the strongest evidence against apoB superiority?**

"[slow] The Emerging Risk Factors Collaboration, 302,430 participants: broadly similar association shapes for apoB and non-HDL cholesterol. [pause] Sniderman's own INTERHEART re-analysis: adjustment among strongly correlated lipids can produce unstable claims of independence. [pause] And the Mendelian randomisation studies disagree by specification, Richardson favouring apoB, Helgadottir favouring cholesterol content. [stop]" (Section 2.3; refs 16, 17, 11, 18. DOIs verified for 16, 11 and 18; ref 17 is listed in the thesis with a PMID but did not resolve in the tools used here.)

**2.3c. Does NMR settle whether it is particle number or cholesterol per particle?**

"[firm] It moves the argument toward particle number without closing it. [slow] In 89,422 statin-free UK Biobank participants with 3,821 coronary events, cholesterol molecules per LDL particle were not associated with disease after adjustment for particle concentrations, hazard ratio 1.03, and apoB correlated 0.99 with LDL-particle concentration. [concede] The adjustment there was for VLDL, LDL and HDL particle concentrations together; my text says LDL-particle concentration and I will make it exact. [stop]" (Section 2.3; ref 61, DOI verified 10.1161/JAHA.123.029552.)

**2.3d. Where does lipoprotein(a) sit?**

"[firm] As a particle-composition and prognostic modifier. [pause] Each Lp(a) particle carries one apoB-100 and contributes variably to measured LDL-C, but its risk is not reducible to LDL-particle concentration. Extreme Lp(a) has been described as a risk equivalent in HeFH. It is in SAFEHEART-RE and the FH-Risk-Score. [stop]" (Section 2.3; refs 53, 56, 64, DOI-verified.)

### Section 2.4 Why discordance in FH remains uncertain

**2.4a. Why is equipoise appropriate rather than an expectation of a positive result?**

"[slow] Two priors oppose each other. Receptor impairment could tighten the mass-to-particle relation so that LDL-C is an unusually faithful proxy and the residual adds nothing. Or metabolic heterogeneity coexists with the variant and leaves clinically relevant variation around the relation. [pause] Selection and measurement can obscure either. So the chapter was designed to be able to return null. [stop]" (Section 2.4.)

**2.4b. What is the nearest construction to yours, and what is the gap?**

"[firm] Johannesen 2024: expected apoB from regressing apoB on LDL-C among people with triglycerides at or below 1 millimole per litre, in people not taking statins. [pause] Competitors are Kim's residual discordance against coronary calcium progression and Rehman's risk-weighted apoB. [slow] The gap is five-fold: population, FH not general; equation, derived in one FH setting and applied frozen in another; endpoint, incident ASCVD; treatment state, mixed; and validation design, an independent evaluation, which Chapter 4 does not supply. [stop]" (Section 2.4; refs 38, 67, 68, DOIs verified.)

**2.4c. Pan 2026 already validated a frozen excess-apoB equation in UK Biobank. What is left?**

"[firm] Pan's cohort was statin-treated coronary disease and the outcome was mortality; excess apoB carried adjusted hazard ratios of 1.12 for all-cause and 1.24 for cardiovascular death. [pause] That establishes the method outside FH. What was left was FH, incident ASCVD, and an equation derived in an FH registry. [stop]" (Section 2.4; ref 39, DOI verified 10.1186/s12944-026-02928-z.)

### Section 2.5 Ratio versus residual

**2.5a. Why can adjusting a ratio for its denominator not fix the coupling?**

"[slow] Because if the ratio is X over Y, and Y enters the model too, you are comparing a transformed variable with one of its own components. Archie described the error in 1981; Lolli showed it recently in a workload ratio. [pause] A quotient also assumes proportionality across the range, and in FH the range spans treated and untreated states. [stop]" (Section 2.5; refs 19, 20, DOIs verified.)

**2.5b. Write the residual on the board.**

"[slow] Residual apoB for person i equals observed apoB minus expected apoB given that person's LDL-C, under a prespecified equation. Positive means more particles than the cholesterol mass implies. [pause] It is not a count per unit cholesterol, not a category, not a particle-size class. [stop]" (Section 2.5.)

**2.5c. Is a threshold ever legitimate?**

"[firm] Yes, when its assumptions are appropriate, its components reliably measured, and its portability shown rather than assumed. [pause] An outcome-selected threshold in a small cohort captures sample-specific variation and invites regression to the mean; Barnett's paper and its correction are the reference. [stop]" (Section 2.5; ref 71, DOI verified 10.1093/ije/dyh299.)

**2.5d. What did the erratum to your discordance paper change?**

"[firm] Typographical errors in the published highlights. The third highlight had been truncated; the four bullets were replaced. [pause] No estimate, table, method or conclusion changed. The thesis's retirement of the ratio estimate is the thesis's re-analysis, not the erratum. [stop]" (Sections 2.5, 4.3; refs 51, 72, DOIs verified 10.1016/j.jacl.2025.11.008 and 10.1016/j.jacl.2026.03.024.)

### Section 2.6 Treatment and reconstruction

**2.6a. Why is medication status "a care-process variable"?**

"[firm] Because in observational data it encodes an intervention and the severity that prompted it, plus access and adherence. [pause] Its coefficient reflects confounding by indication as much as pharmacology. [stop]" (Section 2.6.)

**2.6b. What does the Simon Broome Register tell us about reconstruction?**

"[slow] That a stable mean relation between lipid measures can coexist with wide individual limits of agreement. Soran and colleagues showed it for non-HDL versus LDL-C. [pause] So a reconstruction can be approximately unbiased at cohort level and unreliable for one person. Correlation does not demonstrate agreement. [stop]" (Section 2.6; ref 23, DOI verified 10.1097/MOL.0000000000000692.)

### Section 2.7 Ascertainment and spectrum effects

**2.7a. What did Mourre 2025 show and why is it not your analysis?**

"[slow] Cascade-screened individuals started statins about 14 years earlier and had fewer cardiovascular events than opportunistically screened ones, and crude route differences attenuated after accounting for case mix. [pause] It shows route conditions the treatment and outcome spectrum. It is not residue-matched or within-pedigree, and it cannot show that route changes biology. [stop]" (Section 2.7; ref 40, DOI verified 10.1093/eurjpc/zwaf234.)

**2.7b. What is the polygenic alternative explanation?**

"[firm] Trinder 2024: clinically diagnosed FH in a Canadian registry had higher LDL-C and more ASCVD than genetically identified FH in UK Biobank, with a different polygenic background. [pause] Residue matching cannot exclude that. Polygenic scores were not available in the Welsh frame. [stop]" (Section 2.7; ref 41, DOI verified 10.1161/ATVBAHA.123.320287.)

**2.7c. UK Biobank is 5 per cent of invitees. Does reweighting fix it?**

"[concede] Reweighting on measured determinants can reduce bias in exposure-outcome associations, as van Alten and Schoeler showed, but it cannot recover unmeasured determinants or people outside the sampling frame. [pause] Large sample size reduces random error. It does not remove selection. [stop]" (Section 2.7; refs 6 to 8, DOIs verified.)

### Section 2.8 Clinical criteria and case-finding

**2.8a. Why is "missing" not "absent" for a tendon examination?**

"[firm] Because an unrecorded examination is not a documented negative one, and an unknown pedigree is not a negative family history. [pause] Assigning zero to unavailable information lowers sensitivity; optimistic imputation no longer reproduces the published criterion. [stop]" (Section 2.8; refs 74, 75, DOIs verified.)

**2.8b. What did Mohammadnia show?**

"[slow] In genetically confirmed patients, algorithm sensitivity differed when assessment relied on coded records rather than all available information. [pause] Criterion performance and extraction performance are different evidence objects. [stop]" (Section 2.8; ref 76, DOI verified 10.1093/ehjdh/ztac059.)

**2.8c. Why is FAMCAT's PPV low, and is that FAMCAT's fault?**

"[firm] Because FH is uncommon in primary care; external validation gave a positive predictive value under 1 per cent. [pause] That is prevalence, not the algorithm. It makes confirmatory workload integral to implementation. [stop]" (Section 2.8; ref 43, DOI verified 10.1016/S2468-2667(19)30061-1.)

### Section 2.9 From identification to prognosis

**2.9a. Why can't you rank SAFEHEART-RE, Montreal and the FH-Risk-Score from their published statistics?**

"[slow] Because they target different populations and outcomes. SAFEHEART-RE: incident events in a registry including people with prior disease. Montreal: principally prevalent disease. FH-Risk-Score: incident events in primary prevention. [pause] A prevalent-disease AUC, an incident concordance and a primary-prevention estimate are not interchangeable. [stop]" (Section 2.9; refs 53 to 56, DOIs verified.)

**2.9b. What does the polygenic-score example teach?**

"[slow] In an FH cohort a polygenic score had a hazard ratio of 1.77 yet moved the C-statistic only from 0.746 to 0.750, p 0.60. [pause] Association without incremental prediction. [stop]" (Section 2.9; ref 79, DOI verified 10.1093/eurjpc/zwag203.)

### Section 2.10 Evaluating prediction

**2.10a. Name the four questions of model evaluation and one measure for each.**

"[slow] Ranking: concordance. Agreement: calibration intercept, slope and curve. Overall error: Brier. Decisions: net benefit. [pause] None substitutes for another. [stop]" (Section 2.10, Table 2.3; refs 30 to 33, DOIs verified.)

**2.10b. Why must competing death be handled?**

"[firm] Because non-cardiovascular death prevents a later first ASCVD event. Treating it as censoring assumes the person remained observable and event-free. [pause] One minus Kaplan-Meier overestimates absolute risk; Aalen-Johansen does not. Cause-specific and subdistribution hazards answer different questions. [stop]" (Section 2.10; refs 34, 35, DOIs verified.)

**2.10c. Why is checklist adherence not validation?**

"[firm] Because STROBE, RECORD, TRIPOD+AI and PROBAST+AI expose omissions and risks. They cannot correct biased sampling, temporal leakage, incomplete outcomes, overfitting, absent calibration or inappropriate transport. [stop]" (Section 2.10; refs 30, 31, 36, 37, 49, DOIs verified.)

### Section 2.11 Life course, sex and modifiers

**2.11a. What is the FH-specific evidence on menopause?**

"[slow] Da Roza 2026: within-person lipid change across the menopausal transition in women with FH, larger in monogenic FH. [pause] Johansen: over 12 years young women with FH accumulated a higher LDL-C burden than men. Klevmoen: loss of statin-treatment years in pregnancy and breastfeeding. [pause] None establishes a uniform effect across genotypes or treatment. [stop]" (Section 2.11; refs 24, 25, 26, DOIs verified.)

### Section 2.12 Synthesis

**2.12a. What single rule does the review converge on?**

"[firm] Each claim must remain the size of the design that supports it. [stop]" (Section 2.12.)

---

## Chapter 3: Shared methods and governance

### Section 3.1 Stance

**3.1a. What are the six components of your "inferential contract"?**

"[slow] Distinct cohort roles. A declared eligibility, unit, time zero, endpoint and estimand for every analysis. Observed and reconstructed values kept as different objects. Family dependence, missingness, competing events and transport treated as science, not afterthought. Every number mapped to a frozen cohort, producing analysis and tier. And a bidirectional lock between methods and results. [stop]" (Section 3.1.)

**3.1b. What are the THESIS-AN packages?**

"[firm] Thirteen thesis-reanalysis identifiers, THESIS-AN-01 to 13, each recording its frame, producing artefact or its unavailability, and freeze date. [pause] A reanalysis is cited to its identifier, never to a publication that does not contain it. [stop]" (Section 3.1; Appendix A.)

### Section 3.2 Estimands

**3.2a. Why can one statistical model not answer the thesis question?**

"[firm] Because the four tasks have different populations, outcomes, temporal structures and estimands: a cause-specific hazard, a mean difference, a binary AUC, and a fixed-horizon cumulative incidence with competing death. [stop]" (Section 3.2, Table 3.1.)

**3.2b. What is the causal boundary of the thesis?**

"[firm] Chapter 5 does not estimate an intervention on ascertainment. Chapters 4 and 7 do not estimate effects of altering residuals or scores. Chapter 6 does not estimate outcomes caused by testing. [pause] The thesis does not acquire those estimands through stronger wording. [stop]" (Section 3.2.)

**3.2c. Two digits, 0.660 and 0.725, appear for different quantities. Explain.**

"[slow] In the CALON-C lineage, 0.7252 is the frozen nine-term concordance in the Welsh primary frame, and 0.6600 is the reverse frozen application of a seven-term Welsh refit to UK Biobank; a different 0.660 is the endpoint-date stress test in UK Biobank. [pause] That is why every value is written with its model, frame and horizon. [stop]" (Sections 3.2, 7.14; Appendix B.)

### Section 3.3 Data sources

**3.3a. What is the frozen registry, and is 4,570 an analytical denominator?**

"[firm] 4,570 deduplicated individuals including 887 LDLR carriers, the TUDOR evaluation resource. It describes the data resource, not any chapter's denominator. [stop]" (Section 3.3.1.)

**3.3b. What is the UK Biobank extract?**

"[slow] 501,936 participants; 3,540 with the operational LDLR-carrier flag; a lipid-clinic-eligible subset of 49,427 with 921 flagged carriers. [pause] The flag is broader than ClinVar P/LP and is not labelled HeFH. [stop]" (Section 3.3.2.)

**3.3c. Why is the flag frequency, one in 142, a problem?**

"[concede] Because it exceeds published estimates of pathogenic FH-variant prevalence, and the extract lacked identifiers for complete re-adjudication. [pause] That is why Chapter 7 restricts to ClinVar P/LP, and why the CALON-C development population is described as an operational-flag cohort. [stop]" (Appendix B.2; refs 4, 5, DOIs verified.)

**3.3d. What is disclosure control?**

"[firm] Aggregate results only: no participant or family identifiers, no exact dates, no row-level predictions, no rare combinations. [pause] Statistical verification and governance approval are separate release gates. [stop]" (Section 3.3.4.)

### Section 3.4 Data engineering

**3.4a. What is the row-count ledger?**

"[firm] Pre- and post-join row counts and unique-participant counts, with unmatched keys and duplicates recorded. [pause] Unexplained multiplicative expansion blocks the result; it is not repaired by silent deduplication. [stop]" (Section 3.4.1.)

**3.4b. Why are outliers not excluded at 1.5 IQR?**

"[firm] Because extreme FH measurements can be genuine. Values are checked against units, assay and coding provenance; influence is assessed by diagnostics; named sensitivities examine influential points without altering the source record. [stop]" (Section 3.4.2.)

**3.4c. How was LDL-C derived in Wales?**

"[slow] Friedewald below the documented triglyceride threshold, 4.5 millimoles per litre in Chapter 5, and direct assay above it. [pause] Registry-wide assay homogeneity is not assumed; measurements span calendar periods and platforms. [stop]" (Sections 3.4.3, 5.2.)

**3.4d. What is your missing-data policy?**

"[firm] No programme-wide mandatory method. Complete-case analysis is acceptable only when its assumptions and denominator are stated. Multiple imputation, where used, must respect outcome, interactions, survival structure and clustering, and happen inside resampling splits. [pause] Missing, unrecorded, structurally unavailable and examined-but-absent stay distinct. [stop]" (Section 3.4.4.)

### Section 3.5 Outcomes, time zero and competing events

**3.5a. Why not pool the event sets across chapters to gain power?**

"[firm] Because they are incompatible endpoint lineages. Chapter 4 retains four confirmation-dependent events; the CALON-C lineage has uncertain component dates for 142 of 289 events. [pause] A larger event count does not supersede a narrower source-locked endpoint. [stop]" (Section 3.5.)

**3.5b. What is time zero in each chapter?**

"[slow] Chapter 4: the baseline measurement defining the apoB-LDL-C phenotype. Chapter 5: none, it is a phenotype contrast. Chapter 6: assessment, a classification. Chapter 7: baseline free of atherosclerotic disease, MACE at five and ten years with competing death. [stop]" (Section 3.5.)

**3.5c. How many competing deaths, and why did the count change?**

"[slow] 224 in UK Biobank and 35 in each Welsh frame after a transparent death-first reconstruction; legacy lineages had reported 193 and one. [pause] Reconciling the counts did not complete CALON-C absolute-risk validation, because the historical fitted object and imputation pathway are incomplete. [stop]" (Section 3.5; Appendix B.4.)

### Section 3.6 Treatment reconstruction

**3.6a. State the three conventions and what each can support.**

"[slow] A: intensity-specific residual fractions, for discordance and TUDOR; supports agreement and error analysis. B: population-average scalars, 0.70 and 0.80, for CALON-C; supports rank stability across a factor sweep. C: a regression on 684 Welsh pairs, untreated non-HDL-C equals 5.2299 plus 0.3345 times treated, for CALON-5; supports group-level targeting of the untreated quantity. [pause] None supports an individual untreated value. [stop]" (Section 3.6, Table 3.4.)

**3.6b. Why do you have three?**

"[concede] Because three studies were built at different times for different modelling needs, and the record does not establish one rationale unifying them. Appendix D says so. [firm] They are kept separate and never used to validate one another. [stop]" (Section 3.6; Appendix D.)

### Section 3.7 Statistical principles

**3.7a. What is your position on p-values?**

"[firm] Effect estimates and 95 per cent intervals lead. P-values are compatibility summaries. An interval including zero is 'no statistically detectable difference in that sample', never equivalence, which needs a prespecified margin. [stop]" (Section 3.7.1; ref 80, DOI verified 10.1001/jama.2012.87802.)

**3.7b. Which diagnostics were completed and which were not?**

"[concede] Reference R diagnostics exist for reconstructed CALON-C states. For CALON-5 the record documents the proportional-hazards check and whole-procedure bootstrap optimism, but not influence or functional-form diagnostics for the log-ratio term. The exact Chapter 4 proportional-hazards object was not recovered. [firm] Producing a coefficient never establishes that assumptions passed. [stop]" (Section 3.7.1.)

**3.7c. Family dependence, chapter by chapter.**

"[slow] Chapter 5: pedigree-clustered errors, wild cluster bootstrap, pedigree fixed effects. Chapter 7: registry resampling grouped by family, both models on the same resamples. CALON-C high-bar rerun: FamilyNumber for folds and bootstrap. [concede] TUDOR: DeLong treats participants as independent, and no family-aware optimism estimate exists. UK Biobank: singleton clusters because the extract had no kinship field; independence is not assumed. [stop]" (Section 3.7.1.)

**3.7d. What is your multiplicity position?**

"[firm] No programme-wide correction. Estimates are read by interval and by prespecification. Where a source applied a procedure it is reported with its family: BH-FDR across the Paper 14 register; Holm in the CALON-C lineage. [stop]" (Section 3.8.)

### Section 3.8 Prediction development and transport

**3.8a. Name the five model states and the only one that is external validation.**

"[slow] Frozen application. Target-adapted preprocessing. Recalibration. Updating. Target refitting. [pause] Only the first evaluates the source model unchanged. [stop]" (Section 3.8; refs 28, 29, 32, DOIs verified.)

**3.8b. What is the sequential interval rule?**

"[slow] Baseline and slope re-estimated where the slope interval excludes one; baseline alone where only the observed-to-expected interval excludes one; otherwise no update. [pause] It is not a closed testing procedure, and the thesis says so. [stop]" (Section 3.8.)

**3.8c. What is the twenty-event floor?**

"[firm] Comparisons below 20 events at a declared horizon are not interpreted. [pause] Precision for concordance is governed by events and, for paired comparisons, by the correlation between scores, not by cohort size. [stop]" (Section 3.8.)

**3.8d. What does the Brier score not tell you?**

"[firm] It is not a calibration measure, and when events are rare a small Brier can coexist with poor ranking. It is stated once in Section 3.8 and referred to elsewhere. [stop]"

### Section 3.9 Reporting and evidence maturity

**3.9a. Recite the evidence tiers.**

"[slow] Published or corrected. Locked internal. Externally evaluated. Exploratory. Developmental. Null or inconclusive. Gated or excluded. [pause] And the tier codes: T1 version of record, T2 accepted or in press, T3 submitted, T4 frozen thesis analysis, T5 developmental, T6 not admissible. [stop]" (Section 3.9, Table 3.6; Appendix A.1.)

**3.9b. What is currently gated?**

"[slow] CALON-C absolute risk, calibration, overall error and utility. No interval for CALON-5's increment over age and sex. TUDOR's prospective calibration, utility and impact. The reconstruction-uncertainty propagation, frozen as a protocol but not run. [stop]" (Section 3.9.)

**3.9c. What is the no-do-harm adoption gate?**

"[firm] Before any construct moves toward use: no material loss in calibration, stability, subgroup performance or utility on transport; no reduction in testing or treatment where a known variant, pedigree or physical sign already indicates it; no substitution for sequencing or guideline-directed lipid lowering; then prospective assessment. [pause] Discrimination alone and retrospective decision curves do not pass it. [stop]" (Section 3.9.)

### Section 3.10 Reproducibility

**3.10a. What does "zero silent correction" mean?**

"[firm] Every change to a denominator, term, value or source status, and every deletion or move, is recorded in the change ledger. [pause] The audit package does not retrospectively supply missing participant-level provenance; alternative builds stay quarantined until adjudicated. [stop]" (Section 3.10.)

**3.10b. Which results are independently reproduced and which are source-locked?**

"[concede] The final CALON-C high-bar rerun met the manifest-level target, and the CALON-5 lock records its object under one hash. Several inherited Chapter 4 to 6 results did not retain the whole chain and remain source-locked rather than independently reproduced. [stop]" (Section 3.10.)

---

## Chapter 4: Measure. Do particles carry information at the same LDL-C?

### Section 4.1 The measurement question

**4.1a. State the chapter's question and its non-claims.**

"[firm] At the same LDL-C, does apoB-containing particle concentration carry additional information about first ASCVD events in rare LDLR or APOB variant carriers? [pause] It does not estimate the causal effect of lowering apoB at fixed LDL-C, it does not compare apoB-guided with LDL-C-guided treatment, and it does not validate a bedside calculator. [stop]" (Sections 4.1, 4.2.)

**4.1b. Why a "staged correction" rather than simply presenting the residual?**

"[firm] Because the published paper used a ratio and a threshold, and the residual is post hoc relative to it. [pause] The chapter preserves the published analysis, diagnoses its sparse-data, temporal and measurement limits, then re-specifies the exposure. The three stages have dates: the paper was submitted in August 2025 and published in March 2026; the residual reanalysis was frozen on 12 July 2026 as THESIS-AN-01. [stop]" (Section 4.1.)

### Section 4.2 Cohorts, measurements and estimands

**4.2a. Describe the two frames precisely.**

"[slow] The published Welsh study: 424 genotype-confirmed adults, median follow-up 9.1 years, 61 with ASCVD of whom 41 were prevalent and 20 incident. [pause] The thesis primary: 1,461 UK Biobank participants selected by an operational rare LDLR/APOB variant rule, 75 incident events, none of which appears in the published paper. [stop]" (Section 4.2; ref 51, DOI verified 10.1016/j.jacl.2025.11.008.)

**4.2b. Were participants with prevalent ASCVD excluded from the 1,461?**

"[firm] Follow-up began after the baseline measurement and the outcome is first incident ASCVD. [concede] The exclusion rule and the number excluded are not stated explicitly in Chapter 4 as they are in Chapter 7, and I will add that sentence. [stop]" (Sections 3.5, 4.2; see Appendix 1, item P7.)

**4.2c. Four of 75 events are unconfirmed. Why keep them?**

"[firm] Seventy-one of the 75 are supported by a core ASCVD code, none arises solely from aortic-stenosis coding, and four await procedure or cause-of-death confirmation. [pause] They stay in the source-locked primary because the influence and deletion sensitivities show the association does not depend on the most influential events. [concede] Their status remains a load-bearing outcome limitation and is stated as such. [stop]" (Section 4.2.)

**4.2d. What is the difference between the local and the transported residual?**

"[slow] The local residual is fitted within UK Biobank and standardised by the UK Biobank SD. The transported residual uses the apoB-on-LDL-C equation fitted in Wales and the Welsh residual SD, applied unchanged. [pause] They answer different questions and are not interchangeable. [stop]" (Section 4.2.)

**4.2e. Three estimands were kept separate. Name them.**

"[slow] Association: the cause-specific rate of first incident ASCVD, by Cox. Absolute incidence: ten-year cumulative incidence with competing death, by Aalen-Johansen. Incremental performance: whether the residual changed concordance beyond an existing predictor set. [pause] Evidence for one does not establish the others. [stop]" (Section 4.2.)

### Section 4.3 The published ratio analysis

**4.3a. What did the published paper actually find?**

"[slow] In 424 carriers, an apoB to LDL-C ratio of 0.31 grams per millimole or above gave an AUC of 0.726 for combined prevalent and incident ASCVD, sensitivity 41.2 per cent, specificity 86.7 per cent. Age alone gave 0.845. [pause] In 14 treatment-naive pairs the ratio ran from 0.18 to 0.38, a 2.1-fold spread. [pause] The ridge-penalised threshold hazard ratio was 38.55 with an interval from 3.72 to 399.36, and the paper called itself hypothesis-generating. [stop]" (Section 4.3.)

**4.3b. Your version of record contains three different proportions for the same comparison. Explain.**

"[concede] Yes. The prose reports 27.6 versus 11.8 per cent; Table 3 reports 25 of 121, 20.7 per cent, versus 36 of 303, 11.9 per cent; and a third rendering uses a 0.3 boundary with 46.9 versus 22.4 per cent. [firm] The thesis preserves all three as version-of-record inconsistencies rather than choosing one, and none supplies the thesis's active estimate. The erratum does not address them. [stop]" (Section 4.3.)

**4.3c. What were the paper's limitations, in your own words?**

"[slow] ApoB missing in 40 per cent, handled by 20 imputations. Only 20 incident events. A dichotomised, noisy ratio. And a headline comparison that mixed disease before lipid measurement with disease during follow-up. [stop]" (Section 4.3.)

**4.3d. You have two congress abstracts on the same cohort. Are they replication?**

"[firm] No. One study, two abstract publications with near-identical text, 425 patients against the paper's 424, using an LDL-conditioned definition rather than the ratio. [pause] They are earlier, differently specified outputs from a near-identical cohort, disclosed as such, and no estimate from them is used. [stop]" (Section 4.3; ref 66, not indexed in PubMed.)

### Section 4.4 Why the historical hazard ratio was unstable

**4.4a. Walk me through the sparse-data geometry.**

"[slow] Sixty-one disease observations, 20 incident. Cross-classify a binary exposure with several covariates and many cells are near-empty. That permits near-separation and highly mobile coefficients. [pause] Penalisation prevents numerical divergence; it cannot manufacture information absent from the risk sets. [stop]" (Section 4.4; ref 84, DOI verified 10.1111/j.0006-341X.2001.00114.x.)

**4.4b. Why does temporal ordering matter here?**

"[firm] Because 41 of 61 disease observations were prevalent. Prior disease influences later treatment, lipids and surveillance, so comparing a contemporary ratio with earlier disease risks reverse ordering and confounding by indication. [pause] Restricting to incident events restores temporality and leaves 20 events, exactly where the model became unstable. [stop]" (Section 4.4.)

**4.4c. Is 0.31 a biological constant?**

"[firm] No. It is a study-specific cut-point. Published discordance studies define excess apoB by a residual from a reference equation, not a fixed ratio, and the bounded search found no independent validation of this orientation or boundary. [stop]" (Section 4.4; refs 38, 67, DOIs verified.)

**4.4d. What happened to the Firth estimate and the negative controls?**

"[firm] The Firth estimate mentioned in programme documentation could not be traced to a verified source value or producing object and is withheld. The negative-control analyses were non-estimable from the aggregate inputs. [pause] The only active adjusted threshold re-analysis is 1.24, p 0.44, in the UK Biobank frame. [stop]" (Section 4.4.)

### Section 4.5 Residual re-specification and robustness

**4.5a. Give me the primary estimate and its adjustment ladder.**

"[slow] Hazard ratio 1.356 per standard deviation, 95 per cent interval 1.070 to 1.718, in 1,461 participants with 75 events. [pause] With age and sex: 1.556. Adding log-triglycerides: 1.469. Adding HDL-C and diabetes: 1.356, the primary. Age as the timescale with delayed entry: 1.467. [pause] The movement is consistent with shared information between residual apoB and metabolic risk factors; it does not identify a causal specification. [stop]" (Section 4.5, Table 4.2.)

**4.5b. Why is the competing-risk estimate "non-load-bearing"?**

"[firm] Because it was an inverse-probability-of-censoring-weighted Fine-Gray-style approximation, not an exact fit, and its implementation was not validated under realistic censoring. [pause] It gave 1.334, directionally similar, and is reported as unverified. The estimands also differ: rate among the event-free versus the cumulative-incidence process. [stop]" (Section 4.5; ref 34, DOI verified 10.1161/CIRCULATIONAHA.115.017719.)

**4.5c. What is the absolute-incidence result?**

"[slow] Aalen-Johansen ten-year cumulative incidence 7.38 per cent in local residual quartile 4, 366 participants, against 2.37 per cent in quartiles 1 to 3, 1,095 participants. [firm] A locally defined descriptive grouping. Not a threshold, not calibrated risk. [stop]" (Section 4.5.)

**4.5d. Does the residual improve discrimination?**

"[concede] Barely, at this precision. Overall concordance 0.697. Adding the residual to age and sex: plus 0.050, interval 0.003 to 0.108, but that base had no LDL-C. Beyond LDL-C and conventional predictors: about plus 0.011, no frozen interval. [pause] A marker can keep a stable hazard coefficient and add little ranking because age already orders many people, the effect is moderate, and 75 events limit precision. [stop]" (Section 4.5.)

**4.5e. You tested apoB information in a prediction model and it failed. Why is that in this chapter?**

"[firm] Because it is the programme's only direct test of whether the measurement signal adds ranking to a model. In the CALON-F grey zone of 5 to 20 per cent predicted risk, 1,685 of 3,333 scored with 218 events, adding log apoB over LDL-C changed concordance by plus 0.0146, interval minus 0.0050 to 0.0330; Lp(a) minus 0.0038; both plus 0.0118. [pause] It used the ratio, an unfrozen model and a different frame, so it does not test the residual. The link from measure to predict is carried forward as a hypothesis. [stop]" (Section 4.5; ref 22, developmental manuscript.)

**4.5f. What is the whole-population result for?**

"[firm] Context only. 406,902 participants, 17,680 events, hazard ratio 1.20 per SD, interval 1.18 to 1.22, retained aggregate. [pause] It is consistent with the population literature and does not enlarge the FH claim. [stop]" (Section 4.5; refs 13 to 15, DOIs verified.)

### Section 4.6 Retirement of the ratio threshold

**4.6a. On what criteria did you retire the published hazard ratio?**

"[slow] Sparse incident events. An interval spanning two orders of magnitude. Mathematical coupling. Dichotomisation. No independent validation. And failure of the adjusted threshold association to transport: 1.24, interval 0.72 to 2.16. [pause] Retirement is not repudiation: the paper identified the measurement problem and made its uncertainty visible. [stop]" (Section 4.6.)

**4.6b. Could a clinician use 0.31 tomorrow?**

"[firm] No. Action would need reproducible classification across assays, an absolute-risk model with calibration, a justified action threshold, decision-curve evidence and prospective assessment. None exists. It must not be described as a treatment trigger. [stop]" (Section 4.6.)

### Section 4.7 Biological interpretation

**4.7a. What is the residual, biologically?**

"[slow] ApoB unexplained by LDL-C under the fitted relation. Positive values may arise from more cholesterol-depleted LDL particles, triglyceride-rich remnants, lipoprotein(a), treatment-related compositional change, assay effects, or combinations. [firm] The coefficient cannot allocate among them. [stop]" (Section 4.7.)

**4.7b. Did remnant cholesterol explain it?**

"[slow] Calculated remnant cholesterol accounted for an estimated 16.0 per cent of the association, interval minus 0.9 to 44.6. NMR remnant cholesterol 6.68 per cent, interval minus 4.34 to 39.54, in 1,450 participants with 75 events. [pause] Both cross zero; compatible with nothing and with a relevant contribution. [stop]" (Section 4.7.)

**4.7c. In the joint model LDL-C was protective. Is it?**

"[firm] No. Conditional on both, apoB was 4.94 and LDL-C 0.26. With predictors this correlated each coefficient is the residual variation after conditioning on the other, and reversal reflects unstable partitioning. [pause] That model is shown to explain why the residual is the interpretable representation. [stop]" (Section 4.7.)

**4.7d. What did the polygenic score add?**

"[slow] Association at 1.43 per SD, interval 1.12 to 1.84, in 1,442 participants with 61 events, with concordance 0.569 and 64.1 per cent missing. [pause] It cannot explain the residual or support an incremental claim. [stop]" (Section 4.7.)

**4.7e. Does the Mendelian randomisation literature convert your association into a causal target?**

"[firm] No. Richardson's multivariable analysis gave apoB an odds ratio of 1.92 while LDL-C attenuated; Helgadottir emphasised cholesterol content. Neither converts an observational residual into a treatment contrast, and no intervention selectively altered residual apoB. [stop]" (Section 4.7; refs 11, 18, DOIs verified 10.1371/journal.pmed.1003062 and 10.1093/eurjpc/zwac219.)

### Section 4.8 Frozen-equation application

**4.8a. Why freeze the Welsh equation at all?**

"[firm] Because a locally fitted residual may exploit the receiving cohort's own apoB-LDL-C distribution. Freezing the derivation relation and the Welsh scaling is the stronger test. [pause] The frozen residual gave 2.22, interval 1.40 to 3.53, in UK Biobank. [stop]" (Section 4.8.)

**4.8b. Is the effect twice as large in UK Biobank?**

"[firm] That inference is not permitted. A Welsh-scaled SD is a different absolute apoB departure; assays, treatment and event structure differ. [pause] What transported is direction and an interval excluding unity. Not magnitude. [stop]" (Section 4.8.)

**4.8c. Is this external validation?**

"[firm] It is external evaluation of a measurement construct across an ascertainment boundary. Not clinical validation of a risk model: no baseline hazard, probability, threshold or decision curve was transported. [stop]" (Section 4.8.)

### Section 4.9 The never-treated boundary

**4.9a. Why not analyse only untreated carriers?**

"[concede] I did, as a boundary: 1.248 per SD, interval 0.885 to 1.760, on 30 events, retained aggregate. It neither confirms nor refutes. [pause] Never-treated participants are selected, younger and less severe, so removing treatment adds another selection mechanism. The primary claim stays in the mixed frame. [stop]" (Section 4.9.)

### Section 4.10 Contribution and limits

**4.10a. State the answer to the particles part in three sentences.**

"[firm] At the same LDL-C, particle concentration carried additional information about first ASCVD events. A residual defined in one setting kept its association when applied unchanged in another; the ratio threshold did not. Whether that information improves risk ranking in a model was not shown. [stop]" (Section 4.10.)

**4.10b. What did you not compare the residual with?**

"[concede] Risk-weighted apoB, Rehman's weighted sum of apoB, triglycerides and lipoprotein(a), which improved Harrell's C over apoB alone in general populations. [pause] So I cannot say which representation carries more prognostic information; that is limitation L10 and a resolving-study item. [stop]" (Section 4.10; ref 68, DOI verified 10.1093/eurheartj/ehaf1124.)

**4.10c. Perimeter of the claim, in one breath.**

"[slow] Observational; 75 events, four unconfirmed; the Welsh source had 40 per cent apoB missingness and 20 events; residualisation depends on its derivation relation and assay; treatment history incompletely characterised; never-treated analysis underpowered; one receiving cohort; family and referral structure in Wales, volunteer selection in UK Biobank. [stop]" (Section 4.10.)

---

## Chapter 5: Select. Does the route into care change the visible phenotype?

### Section 5.1 Selection as part of the phenotype

**5.1a. What is the estimand, and why is it not causal?**

"[firm] The observed difference in assigned untreated LDL-C between carriers entering the same service by different routes, after controlling approximately for the affected LDLR residue and then reducing age, pedigree and provenance differences. [pause] Ascertainment is not a biological intervention, and proband status is partly defined by the phenotype under study. The group contrast describes selected distributions, not an individual counterfactual. [stop]" (Section 5.1.)

**5.1b. Explain Figure 5.1.**

"[slow] Underlying LDL-C, variant severity and family history influence the route by which a carrier is identified. Route influences age at encounter and treatment history, which influence the observed untreated LDL-C. Inclusion in the sample is the selection node the analysis conditions on. [firm] There is no arrow from route to underlying biology. It is conceptual, not an analysis. [stop]" (Section 5.1.)

**5.1c. Distinguish expressivity, recognition and prognosis.**

"[slow] Expressivity: the phenotype observed among carriers. Recognition: whether that phenotype prompts detection. Prognosis: outcomes under a defined time origin. [pause] This chapter examines expressivity as observed through different recognition pathways; its event analyses are secondary. [stop]" (Section 5.1.)

### Section 5.2 Design, cohorts and provenance

**5.2a. Give me the denominators.**

"[slow] 1,904 eligible Welsh LDLR carriers: 815 probands and 1,089 cascade-detected relatives; computable untreated LDL-C in 780 and 889; 235 with no computable value contribute nothing. The residue-matched primary frame: 1,100 carriers, 464 probands and 636 relatives, at 50 residues in 488 pedigrees, median two carriers per pedigree, 48 per cent of pedigrees contributing one carrier. [pause] UK Biobank: 3,540 flagged, 2,398 untreated with a measured LDL-C. [stop]" (Section 5.2, Table 5.1.)

**5.2b. What is the provenance hierarchy for the untreated value?**

"[slow] First, a documented pretreatment value. Second, a measured value in a carrier with no recorded treatment. Third and only then, back-calculation from a treated value by a drug-specific rule. [pause] Directly observed in 88.7 per cent of probands and 95.3 per cent of relatives; reconstructed in 11.3 and 4.7. [stop]" (Section 5.2.)

**5.2c. Are the reconstructed rows a random subset?**

"[firm] No, and that matters. Carriers with reconstructed values had far more prior hard events: 78.4 versus 13.2 per cent in probands, 61.9 versus 4.4 in relatives. Individual reconstruction agreed poorly with documented values, r 0.09. The association of untreated LDL-C with prevalent disease differed by source, odds ratio 1.20 in observed rows against 1.00 in reconstructed. [pause] So reconstruction is a cohort-level device, and the observed-only and documented-only restrictions are the checks. [stop]" (Section 5.2.)

**5.2d. Chapter 7 says the registry's drug-specific back-calculation was biased by 4.5 millimoles per litre. Chapter 5 uses a drug-specific rule. Reconcile.**

"[concede] The two chapters use different rules on different pairs for different targets, and the thesis should bring the agreement figures into one table. [firm] What protects Chapter 5 is that its primary contrast does not need the rule: documented pretreatment values only give plus 1.16; directly observed values only give plus 1.09; nine alternative rules give plus 1.00 to plus 1.09. [stop]" (Sections 5.2, 5.4, 7.2; Appendix 1, item P11.)

**5.2e. What does residue matching control, and what does it not?**

"[slow] It blocks the explanation that groups carried variants in different receptor regions. It does not control the exact substitution, functional severity, polygenic background, lipoprotein(a), metabolic context, calendar era, treatment history or referral practice. [pause] Exact-allele matching controls the substitution; within-pedigree comparison holds the familial allele and much of the shared background, but only 236 of 488 pedigrees contain both routes. [stop]" (Section 5.2.)

**5.2f. How much LDL-C variance does the variant explain?**

"[slow] Little. In the TUDOR analyses the exact variant had an intraclass correlation of about 0.145, leave-one-variant-out 0.127 to 0.149. UK Biobank carriers of the same LDLR variants as Welsh clinic carriers had LDL-C 1.36 millimoles per litre lower, interval minus 2.03 to minus 0.76. [pause] And a predicted measure of variant function tracked untreated LDL-C in relatives, rho 0.26, but not in probands, 0.02, a post hoc pattern consistent with selection on the phenotype. [stop]" (Section 5.2; ref 21, DOI verified 10.1016/j.jacl.2026.06.030.)

**5.2g. What does the newborn evidence show?**

"[slow] In the Norwegian cascade programme, LDL-C did not differ between newborns with null and non-null variants, 1.26 versus 1.18 millimoles per litre, p 0.552. [concede] The contrast is within the 61 FH newborns of a 113-newborn sample; my text says 113 and I will correct it. It speaks to screening at birth, not adult phenotype by route. [stop]" (Section 5.2; ref 85, DOI verified 10.1093/eurheartj/ehaf815.)

### Section 5.3 Primary result and attenuation

**5.3a. The primary estimate, with its inference.**

"[slow] Plus 1.09 millimoles per litre, 95 per cent interval 0.83 to 1.36, p under 0.001, with residue fixed effects and pedigree-clustered standard errors. HC1-robust and a Rademacher wild cluster bootstrap with 1,999 pedigree resamples gave the same interval. [pause] These are Paper 14 R4 estimates, not independently re-run for the thesis. [stop]" (Section 5.3.)

**5.3b. Age was the largest threat. What did you do?**

"[slow] Probands were median 60.2 years, relatives 38.6. Linear age adjustment: plus 0.81, interval 0.52 to 1.09. Natural spline with an age-by-sex interaction: plus 0.78. Adding sex to linear age: plus 0.82. Restricting to ages 40 to 65, prespecified, medians 55 and 52: plus 0.65, interval 0.24 to 1.07, n 457. [stop]" (Section 5.3.)

**5.3c. The within-pedigree result. Is it real?**

"[long pause] Detectable, but imprecise, and the conclusion depends on the variance convention. [slow] Plus 0.39, p 0.023, pedigree bootstrap 0.04 to 0.74, directly observed only plus 0.41. Under a degrees-of-freedom correction that counts every pedigree effect, minus 0.06 to 0.85. [pause] Appendix B.10 records that the rung has moved across four manuscript versions and that revision R3 crossed zero. [firm] What it supports is a residual within-family association the design cannot attribute to referral selection or to biology. I would not claim more. [stop]" (Section 5.3; Appendix B.10.)

**5.3d. Why is the ladder "more informative than any coefficient"?**

"[slow] Because its rungs are different operations. 1.09 to 0.81 is adjustment in the same carriers. 0.81 to 0.65 is restriction to a common age band, a change of population. 0.65 to 0.39 is a change of comparison to within-family contrasts, informed by 236 pedigrees. [pause] On fixed membership, pedigree effects without age gave plus 0.59 and with age and sex plus 0.40, so the fall is model, not sample. [stop]" (Section 5.3.)

**5.3e. Where is the sex-stratified within-pedigree estimate from the earlier draft?**

"[firm] Withdrawn in the revision because it could not be refitted on the corrected frame, and the study was not designed or powered for interaction by sex. It was not retained with a stale value. [stop]" (Section 5.3.)

### Section 5.4 Robustness and the three strata

**5.4a. List the robustness rows and their purpose.**

"[slow] Directly observed only, plus 1.09. Documented pretreatment only, plus 1.16. Inverse-probability-of-observation weighting, plus 1.07, because 4.3 per cent of probands and 18.4 per cent of relatives lacked a computable value. Nine back-calculation rules, plus 1.00 to plus 1.09. Non-HDL cholesterol, plus 1.67 crude, plus 1.44 adjusted. Excluding probands with prior ASCVD, plus 0.96. Age, sex and treatment adjusted, plus 0.86, an over-adjustment. Adding BMI in 380 carriers, plus 0.65. [pause] Leave-one-out moved the estimate by at most 0.03 for pedigrees and 0.13 for residues or variants. These re-examine one frame; they are not replication. [stop]" (Section 5.4, Table 5.2.)

**5.4b. Why is the non-HDL difference larger than the LDL-C difference?**

"[slow] Mostly timing. On the same 974 carriers the LDL-C difference was 1.14; at the same first visit the two differences were 1.64 and 1.49. [pause] A small additional remnant component in probands is consistent with, but not shown by, the data. [stop]" (Section 5.4.)

**5.4c. What is the variant-level synthesis and why is its status weaker?**

"[slow] A random-effects meta-analysis across 71 LDLR/APOB variants with at least three carriers per route: plus 0.91, interval 0.75 to 1.07, I-squared 57 per cent. [concede] Its producing output was not located, so it rests on a secondary record, and a later 53-variant rerun gave a different pooled value and is treated as non-equivalent. [stop]" (Section 5.4; THESIS-AN-02.)

**5.4d. Is 5.80, 4.76, 3.87 a matched comparison?**

"[firm] No. It compares differing allele spectra. Its use is descriptive: the phenotype encountered as detection moves from referral to cascade to sequencing. [stop]" (Section 5.4.)

**5.4e. The first-versus-later reading result: regression to the mean?**

"[slow] A proband's first pre-treatment reading averaged 2.20 millimoles per litre above the next, interval 1.62 to 2.76, 95 pairs; a relative's 0.91, interval 0.43 to 1.42, 84 pairs; drift with interval not detectable. The empirical reference-change value was 93 per cent against a biological CV of about 7.9 per cent for calculated LDL-C. [concede] Compatible with referral-associated measurement selection and regression to the mean, inseparable from biological change or unrecorded treatment. Descriptive, post hoc, and no share of the between-group contrast is attributed to it. [stop]" (Section 5.4; ref 71, DOI verified.)

### Section 5.5 Cumulative exposure and event associations

**5.5a. Did earlier detection mean lower exposure?**

"[firm] The data cannot say. The manuscript approximated exposure as baseline LDL-C times attained age; in 3,204 UK Biobank carriers with 137 events it added 0.002 to Harrell's C over age, sex and LDL-C. [pause] A single-measurement product is not a life-course integral. [stop]" (Section 5.5.)

**5.5b. Were the per-unit LDL-C associations similar by setting?**

"[slow] Wales, cross-sectional: odds ratio 1.15 per millimole for prior hard ASCVD, 223 events in 1,669; relatives alone 1.14, 63 events. UK Biobank: odds ratio 1.19 for any ASCVD; hazard ratio 1.21 for incident, 137 events; exact-allele P/LP 1.10, interval 0.93 to 1.29, 49 events. [firm] No equivalence was tested, the magnitudes are not comparable, and proband analyses are collider-prone. [stop]" (Section 5.5.)

**5.5c. Can I add 1.09 to a cascade relative's value?**

"[firm] No. It would ignore age, exact allele, treatment provenance and between-variant heterogeneity. The clinical response is to preserve the molecular diagnosis, establish provenance, assess full risk context and follow guidance. [stop]" (Section 5.5.)

### Section 5.6 Consequences for identification and prediction

**5.6a. What is the reporting unit for identification research?**

"[firm] Model by ascertainment route by comparator spectrum by data availability. [pause] Discrimination with calibration and threshold yield, route-stratified where information permits. [stop]" (Section 5.6.)

**5.6b. What is "selection leakage"?**

"[slow] A predictor that partly reconstructs the case-definition process: proband status, referral source, pedigree completeness. Whether it is legitimate depends on use; it may prioritise testing within a cascade programme and be tautological in primary care. [pause] Gains attributed to route need a locked ablation, route-stratified evaluation and frozen application. [stop]" (Section 5.6.)

**5.6c. What proportion of P/LP carriers meet the Simon Broome LDL-C criterion?**

"[slow] Among 715 UK Biobank P/LP carriers with no recorded treatment, 84.5 per cent, interval 81.6 to 86.9, were below 4.9. In Wales, 38.1 per cent of probands, 56.9 per cent of relatives and 40.5 per cent of genotype-negative referrals were below it. [pause] Untreated LDL-C separated carriers from genotype-negatives with AUC 0.58 on the proband route and 0.76 on the cascade route. [stop]" (Section 5.6.)

**5.6d. How does ascertainment reach Chapter 7?**

"[slow] The P/LP fraction was 31.4 per cent in UK Biobank and 57.7 per cent in the registry, a gradient in the genetic label itself. [pause] A model developed in one and applied to the other inherits that. [stop]" (Sections 5.6, 7.4.)

### Section 5.7 Contribution and limits

**5.7a. State the answer to the route part.**

"[firm] The same affected LDLR residue carried a higher untreated LDL-C in carriers found through their own phenotype than in relatives found by cascade testing, plus 1.09. The difference narrowed to plus 0.39 within families. [pause] The route into care changes which phenotype is seen, not the biology. [stop]" (Section 5.7.)

**5.7b. The limitations that prevent a larger claim.**

"[slow] Residue matching is positional. The within-pedigree rung rests on 236 pedigrees and its interval depends on convention. Lipoprotein(a) was recorded for 143 of 1,904; apoB was not in this extract; polygenic scores were unavailable. Era and assay heterogeneity. Complete-case covariate models. BMI in 380. UK Biobank volunteer-selected and ancestry-restricted. And the Welsh carriers overlap with the Chapter 4 registry sample by an unestablished number. [stop]" (Section 5.7.)

**5.7c. Which ancestry restriction?**

"[concede] The thesis names an ancestry restriction without defining it. I need to state the group used and confirm consistency across Chapters 4 and 7. [stop]" (Section 5.7; Appendix 1, item P9.)

---

## Chapter 6: Identify. Can carriers be found after treatment?

### Section 6.1 Identification after treatment and selection

**6.1a. What is the contemporary identification question, and how does it differ from Simon Broome's?**

"[firm] Not whether an untreated proband satisfies a score. [pause] Among people already in a lipid, cardiology or genetic-testing pathway, whom should a service prioritise for confirmatory testing when the untreated phenotype is no longer observable? [stop]" (Section 6.1.)

**6.1b. Quantify the failure of a single LDL-C gate.**

"[slow] Among UK Biobank ClinVar P/LP carriers with an untreated LDL-C, 80.2 per cent were below the 4.9 Simon Broome trigger, and a 4.9 gate would have required 404 people sequenced per carrier found. [pause] That is the failure mode a multivariable triage score is meant to reduce. [stop]" (Section 6.1; ref 21, DOI verified.)

**6.1c. What would weaken the chapter's hypothesis?**

"[slow] Absence of a paired discrimination advantage; marked family leakage; poor frozen application; failure to reproduce the recovered score. [pause] Calibration failure would not erase discrimination but would prevent reading the score as a probability. [stop]" (Section 6.1.)

### Section 6.2 Model specification and evaluation states

**6.2a. Specify TUDOR exactly.**

"[slow] Elastic-net logistic regression, mixing parameter 0.5, C equal to 1.0, saga solver, seed 20260518. Eleven declared inputs in fixed order: reconstructed pretreatment LDL-C, HDL-C, triglycerides, total cholesterol, non-HDL-C, age, sex, treatment status, family history of cardiovascular disease, a diabetes-by-LDL interaction, and premature cardiovascular events. Ten contribute; family history is zero-variance with coefficient zero. [firm] No proband-route or ascertainment term. [stop]" (Section 6.2.)

**6.2b. The published paper says ascertainment-aware. The object has no route term. Explain, and tell me why there is no corrigendum.**

"[long pause] [firm] The publication's title says ascertainment-aware and its abstract says the model encodes index-referral versus cascade-relative ascertainment. When the frozen executable pipeline was recovered and reproduced by score readback, its eleven inputs contained no such term. [pause] Every TUDOR result in the thesis refers to the recovered object. No Index-Effect increment is reported, and a number that circulated in earlier programme material is omitted because no producing object supports it. [long pause] [concede] The discrepancy is disclosed in the thesis only. No corrigendum has been requested, and that was my decision. On reflection the order is wrong: a reader of the paper attributes to the model a feature the evaluated object does not contain. I will write to the journal before final submission. [stop]"

If pressed, "Does it invalidate the AUC?": "[firm] No. Readback reproduces the stored predictions within tolerance, so 0.760 is the performance of the object that exists. It invalidates the route-aware claim, not the ranking claim. [stop]" (Section 6.2, Box 6.1; ref 21, DOI verified 10.1016/j.jacl.2026.06.030; see Appendix 1, item P18 on the exact abstract wording.)

**6.2c. Your Atherosclerosis Plus abstract claims "dual external validation in 4,028 genetically confirmed cases". Is that consistent with the thesis?**

"[concede] No, and the thesis corrects it. The abstract draws on the same two cohorts as the chapter and is cited as related work, not independent validation. [pause] The 4,028 figure appears nowhere else in the thesis and I need to reconcile it against the abstract before quoting it. [stop]" (Section 6.2; ref 82, not indexed in PubMed; Appendix 1, item P19.)

**6.2d. Name the evaluation frames and say why they cannot be merged.**

"[slow] Welsh complete-case head-to-head, 1,274 with 311 carriers. Full deduplicated registry, 4,570 with 887. Reciprocal geographic splits, 3,099 with 554 and 1,471 with 333. UK Biobank lipid-clinic-eligible, 49,427 with 921. UK Biobank local sensitivity, 43,594 with 652. Whole UK Biobank, 501,936 with 3,540. And a 649-pair reconstruction resource. [pause] Different eligibility, comparator-availability and carrier-definition rules; merging them would create a result no analysis produced. [stop]" (Section 6.2, Table 6.1.)

**6.2e. Which paired test, and what is wrong with it here?**

"[firm] DeLong's test for correlated AUCs, used in the local UK Biobank paired comparisons and the reconstruction sensitivity. [concede] It treats participants as independent; in a family-structured registry that may understate uncertainty. The Welsh head-to-head differences were reported without a named paired test or interval, which is a gap in my own source report. [stop]" (Section 6.2; ref 87, DOI verified 10.2307/2531595.)

### Section 6.3 Welsh head-to-head discrimination

**6.3a. The result.**

"[slow] In 1,274 participants, 311 carriers: TUDOR AUC 0.760, interval 0.727 to 0.790. Reconstructed eDLCN 0.652, FAMCAT 0.600, Simon Broome 0.569, MEDPED 0.524. Net reclassification improvement against FAMCAT plus 0.204, interval 0.141 to 0.266. [stop]" (Section 6.3, Table 6.2.)

**6.3b. What does that ordering answer, and what does it not?**

"[firm] Given the information represented for every tool under the same complete-case restriction, which score ranked carriers above non-carriers more consistently? [pause] It does not show TUDOR beats the criteria applied with observed untreated lipids, full pedigrees and examination. And 0.760 is not a clinical threshold. [stop]" (Section 6.3.)

**6.3c. Why is NRI secondary?**

"[firm] Because it is sensitive to prevalence, score distributions, category definitions and miscalibration, and a large NRI against a weak comparator need not mean a better decision. [stop]" (Section 6.3; ref 88, DOI verified 10.1093/aje/kwx374.)

**6.3d. Does TUDOR beat pretreatment LDL-C alone?**

"[concede] In one defined local pipeline: 0.747 against 0.705, plus 0.043 by DeLong. An earlier analysis with a different 450-carrier definition gave plus 0.017 without statistical separation. [firm] The increment beyond LDL-C is carrier-definition-sensitive, and I do not claim a stable increment. [stop]" (Section 6.3.)

### Section 6.4 Comparator fairness and the ascertainment boundary

**6.4a. Make the case against your own comparison.**

"[slow] eDLCN ordinarily combines LDL-C with vascular history, tendon xanthomata, corneal arcus, pedigree and molecular evidence; the electronic reconstruction lacked several, and unrecorded findings scored as absent. FAMCAT was built for primary-care records and ran on a specialist-registry subset. Simon Broome and MEDPED depend on untreated cholesterol and physical signs. [firm] The comparison is asymmetrical. It is still clinically relevant because the electronic-data problem is real; it prohibits the claim that TUDOR beats a fully informed specialist. [stop]" (Section 6.4; refs 43, 74 to 76, DOIs verified.)

**6.4b. What is the counterweight?**

"[slow] Akyea 2026. In an Australian tertiary clinic, 885 referred for genetic testing, 267 with an FH-causing variant, DLCN scored with examination, pedigree and untreated lipids had an AUROC of 0.816 against FAMCAT's 0.748. [pause] Scored with that information, the criteria beat an electronic algorithm. TUDOR has not been compared with fully scored DLCN. [stop]" (Section 6.4; ref 45, DOI verified 10.1016/j.jacl.2026.07.009.)

**6.4c. Was TUDOR compared with the Welsh service's own criteria?**

"[concede] No. The Haralambos 2015 criteria, which modify DLCN, are TUDOR's local predecessor, and the increment over them is unknown. [stop]" (Section 6.4; ref 44, DOI verified 10.1016/j.atherosclerosis.2015.03.003.)

**6.4d. Can you partition TUDOR's advantage into treatment-awareness and route?**

"[firm] No. There is no executable route term and no full-versus-minus-route ablation. Descriptive registry AUCs of 0.760 in probands and 0.794 in relatives show no obvious loss in relatives, but case-mix differences prevent attribution. [stop]" (Section 6.4.)

**6.4e. Family leakage?**

"[concede] Unquantified. No family-aware optimism estimate with a located producing record exists, and no fold manifest has been reconstructed. Expected direction: optimism. [stop]" (Section 6.4.)

### Section 6.5 Frozen application and the intended-use boundary

**6.5a. Transport ladder, in numbers.**

"[slow] Head-to-head 0.760. South Wales to rest of Wales 0.732, interval 0.707 to 0.757. Rest of Wales to South Wales 0.770, interval 0.739 to 0.797. Pooled frozen Wales 0.746. Frozen UK Biobank lipid-clinic-eligible 0.669, interval 0.650 to 0.687. UK Biobank target refit 0.756 apparent. Whole UK Biobank 0.631. [stop]" (Section 6.5, Table 6.3.)

**6.5b. Which is the most consequential number?**

"[firm] 0.669. The frozen Welsh score did not preserve its source discrimination in a cohort that differed in age, treatment, prevalence, measurement and volunteer selection. [pause] It does not say which difference caused the loss. And it must not be displaced by the refitted 0.756. [stop]" (Section 6.5.)

**6.5c. Why is 0.756 not external validation?**

"[firm] Coefficients were re-estimated in UK Biobank; relative predictor weights changed; target outcomes informed the model. That is local development. It shows the feature architecture contained reweightable information. It does not say a new service should expect 0.756. [stop]" (Section 6.5.)

**6.5d. What does 0.631 mean?**

"[firm] The intended-use boundary. TUDOR was not built for indiscriminate population screening, and in an unselected local analysis it did not exceed LDL-C alone. A pathway-aware model can be useful after the pathway exists and add nothing before it. [stop]" (Section 6.5.)

**6.5e. Subgroups?**

"[slow] Five-gene-plus-CNV outcome 0.779. Type 2 diabetes 0.648 in Wales and 0.642, interval 0.559 to 0.724, in UK Biobank: a metabolically difficult setting, not a subgroup licence. [concede] A statin-naive 0.801 was reported without a verified denominator or interval and carries no strong conclusion. [stop]" (Section 6.5.)

### Section 6.6 Reconstruction as a measurement model

**6.6a. How large is the reconstruction error?**

"[slow] Mean absolute error 1.20 millimoles per litre in 649 paired participants, correlation 0.32. [pause] Clinically material near thresholds. It reflects adherence, dose, combinations, response and timing that a deterministic factor cannot recover. [stop]" (Section 6.6; Appendix D.)

**6.6b. Why did the AUC not move across reconstruction methods?**

"[slow] Across dose-specific, class-level and fixed-factor approaches the AUC range was under 0.02, paired p above 0.3. [pause] Ranking tolerates some systematic error because other predictors contribute and a shared shift preserves order. [firm] That is a statement about ranking, not agreement. [stop]" (Section 6.6.)

**6.6c. The clinical boundary of a TUDOR score?**

"[firm] A high score may justify priority for confirmatory testing; it does not establish a genotype. A low score should not overrule a known familial variant, tendon xanthomata, a compelling pedigree or expert assessment. [stop]" (Section 6.6.)

### Section 6.7 Calibration, workload and utility

**6.7a. Is TUDOR calibrated?**

"[concede] Not shown to be. In a local UK Biobank stream, slope 1.257, intercept minus 3.03, Brier 0.076, retained aggregate: substantial overprediction in a lower-prevalence setting. [firm] Those values must not be transferred to the locked 49,427 transport cohort, for which external calibration was not available. TUDOR's output is a ranking score. [stop]" (Section 6.7; THESIS-AN-03b.)

**6.7b. The operating points.**

"[slow] In the local 43,594 with 652 carriers: Youden point refers 10,433, captures 409, sensitivity 0.627, PPV 3.9 per cent, 25.5 sequenced per carrier. Top 10 per cent refers 4,360, captures 266, PPV 6.1, 16.4 per carrier. A high-sensitivity point refers 36,417, captures 619, sensitivity 0.949, PPV 1.7. [pause] A statistical optimum encodes none of the harms. [stop]" (Section 6.7, Table 6.4.)

**6.7c. Why do numbers needed to screen differ tenfold?**

"[firm] 2.0 in Wales, 15.1 in UK Biobank: prevalence and pathway, not the model. A Welsh threshold cannot be transferred mechanically. [stop]" (Section 6.7; THESIS-AN-03a.)

**6.7d. Why no decision curve for TUDOR?**

"[firm] Because a decision curve needs valid target probabilities, and TUDOR's are not calibrated for the transport cohort. It is the appropriate next method, after calibration, and even then it would not prove improved outcomes. [stop]" (Section 6.7; ref 33, DOI verified.)

### Section 6.8 Reporting appraisal and gates

**6.8a. Why is there no PROBAST+AI rating for TUDOR when you gave one for CALON-5?**

"[concede] Because the signalling-question assessment was not completed for TUDOR; Table 6.5 gives narrative considerations only. That is inconsistent with Appendix E.4 and I will complete it. [stop]" (Section 6.8; Appendix 1, item P12.)

**6.8b. What blocks independent deployment?**

"[slow] Incomplete target calibration, comparator asymmetry, absent prospective threshold evaluation, absent impact evidence, and the description-object discrepancy. [pause] Reproducibility is no longer the block: the object is recovered and reads back within tolerance. [stop]" (Section 6.8.)

**6.8c. The next validation stage, in order.**

"[slow] Freeze the pipeline with its software environment and readback tests. Fix carrier definition, predictor timing, missing-data policy and intended pathway, with preprocessing inside family-safe resampling. Frozen evaluation in an independently assembled specialist-triage cohort, reporting calibration-in-the-large, slope and curves separately from discrimination. Prespecified thresholds with workload, decision curves and missed-carrier consequences. Then prospective impact. [stop]" (Section 6.8.)

### Section 6.9 Contribution and transition

**6.9a. State the answer to the treatment part.**

"[firm] Reconstructing the untreated phenotype improved the ranking of carriers against reconstructed criteria in the development setting, and less so when the score was applied unchanged elsewhere. [pause] Identification is not prognosis: a high probability of carrying a variant does not say who has an event first. [stop]" (Section 6.9.)

---

## Chapter 7: Predict. Does risk ranking survive a change of cohort?

### Section 7.1 The question

**7.1a. Two meanings of "the same thing" for a risk score?**

"[firm] The order in which the model places patients, discrimination, and the absolute risk it assigns them, calibration. The chapter keeps them apart throughout. [stop]" (Section 7.1.)

**7.1b. Why can a within-FH model never inform whether to treat?**

"[firm] Because lipid lowering is already indicated by the diagnosis. It could inform intensity, timing in young relatives, and allocation of add-on therapy or imaging. [stop]" (Section 7.1; refs 1 to 3, DOIs verified.)

**7.1c. You call both models high risk of bias on your own appraisal. Why, in four points?**

"[slow] CALON-5's terms were chosen with knowledge of its transported registry performance. The two-term model was selected after those results were known. The registry fit rests on few events. And neither cohort is untouched by the programme's development. [stop]" (Section 7.1; Appendix E.4.)

### Section 7.2 Model identities

**7.2a. Specify CALON-5.**

"[slow] Cox proportional hazards, no stated penalty. Five terms: age, male sex, hypertension defined as recorded antihypertensive medication at baseline, diabetes, and the log of untreated non-HDL-C over HDL-C. [pause] The hypertension definition replaced a blood-pressure-or-medication definition on 24 September 2026, before any refit result was seen; the earlier definition is a sensitivity. [stop]" (Section 7.2.)

**7.2b. Where does the untreated lipid value come from?**

"[slow] The registry supplies a recorded pre-treatment value. UK Biobank holds none, so for the 24.7 per cent of biobank carriers on therapy, untreated non-HDL-C is back-calculated as 5.2299 plus 0.3345 times the recorded value, from 684 within-person registry pairs, cross-validated mean absolute error 1.25, R-squared 0.0975. [pause] A drug- and dose-specific equation was not used: dated, dosed prescriptions existed for 27.6 per cent of treated carriers, and the registry's drug-specific rule was biased by 4.5. [stop]" (Section 7.2.)

**7.2c. R-squared of 0.0975 and a slope of 0.33. Defend using it.**

"[concede] It explains almost no individual variation and compresses the treated between-person spread; it was fitted mostly on non-carrier pairs. [firm] It is used for one purpose: group-level targeting of the untreated quantity guidelines specify, so that the two cohorts model the same nominal exposure. It is not accurate for an individual, the error is plausibly differential, and the registry uses recorded values, so the cohorts measure exposure differently. All of that is stated. [stop]" (Sections 7.2, 7.10; Appendix D.3.)

**7.2d. What is the two-term model and why does it exist?**

"[slow] Age and the same lipid ratio. Because CALON-5's terms were chosen knowing their registry performance, selection was repeated on UK Biobank data alone, blind to the registry. It chose these two, the smallest of 25 indistinguishable subsets. In 200 bootstraps age was retained 83 per cent of the time and the lipid ratio 32 per cent, so its identity is a convention. Optimism-corrected concordance 0.612. [stop]" (Sections 7.2, 7.3.)

**7.2e. "Prespecified after the results were known" is a contradiction.**

"[concede] Both words are true of different things, and I should not hide behind the phrase. The selection rule was written into the revision protocol of 23 September 2026 before it was run, but after CALON-5's transported results were known. [firm] The honest label is 'development-only selected', which is why neither model is primary. [stop]" (Section 7.2, Table 7.1.)

**7.2f. Why exclude Lp(a), apoB and imaging?**

"[firm] Not for irrelevance. The question is whether routinely recorded variables support ranking where specialised inputs are unavailable. A model with fewer inputs computes for more patients while omitting information; a richer one performs where its inputs exist and cannot be computed where they were never collected. [stop]" (Section 7.2.)

### Section 7.3 Development and internal validation

**7.3a. Events per parameter, both directions.**

"[slow] Primary: UK Biobank P/LP, 888 complete cases, 74 events, 14.8 per parameter. Reverse: registry P/LP, 469 complete cases, 28 events, 5.6 per parameter, below the prespecified floor of ten; exploratory. The all-carrier reverse fit, with 52 events, is the one interpreted. [stop]" (Section 7.3.)

**7.3b. Why is this "cross-cohort evaluation" and not external validation?**

"[firm] The registry is not untouched: the untreated-exposure equation used its lipid pairs, and CALON-5's terms were chosen with knowledge of its transported registry performance. [stop]" (Section 7.3.)

**7.3c. Internal validation numbers.**

"[slow] Whole-procedure bootstrap: apparent concordance 0.6972, optimism 0.0209, corrected 0.6763; calibration slope corrected to 0.901. No term showed evidence against proportional hazards, 0 of 10 tests, a check with little power on 52 registry events. [stop]" (Section 7.3.)

**7.3d. Did the hazard ratios agree across cohorts?**

"[slow] In the all-carrier fits all five agreed in direction. The lipid term: 2.33 per log unit in UK Biobank, interval 1.62 to 3.35, and 2.40 in the registry, interval 1.34 to 4.31. At the same covariates registry carriers had 1.33 times the biobank hazard, interval 0.83 to 2.13. [stop]" (Section 7.3.)

### Section 7.4 Cohorts, flow and missingness

**7.4a. What did standardising to ClinVar P/LP cost?**

"[slow] UK Biobank fell from 3,209 carriers and 289 events to 1,009 and 87. The registry from 1,639 genotyped attendees and 101 events to 945 and 56. The asymmetry is the pathogenic fraction, 31.4 against 57.7 per cent, itself an ascertainment signal. [pause] At two or more ClinVar review stars, 69.0 per cent of biobank P/LP variants qualified against 97.7 per cent of registry variants; zygosity is recorded in neither. [stop]" (Section 7.4.)

**7.4b. Why is the registry complete-case set not missing at random, and what did you do?**

"[slow] Requiring all five inputs kept 800 of 1,639 all-carriers and 469 of 945 P/LP carriers, driven by unrecorded diabetes: recorded for 54 per cent, hypertension for 99. Those without a recorded diabetes status had about half the hazard, age- and sex-adjusted 0.40, interval 0.21 to 0.77. [firm] Missingness was associated with the outcome, so MAR is doubtful, and I did not impute. The set is selected toward longer-observed, clinically assessed patients; that is limitation L15. [stop]" (Section 7.4.)

**7.4c. Endpoint composition and censoring in the registry.**

"[slow] Of 101 events: 52 infarction or acute coronary syndrome, 16 percutaneous and 19 surgical revascularisations, 8 angina, 6 TIA or stroke. Event times from recorded ages. Censoring survival 0.83 at five years and 0.52 at ten, so ten-year estimates lean on the weights. [stop]" (Section 7.4.)

### Section 7.5 Discrimination

**7.5a. The headline.**

"[slow] Developed in UK Biobank and applied unchanged to the registry, ten-year concordance 0.699, interval 0.611 to 0.790, time-dependent AUC 0.718. Reverse: 0.698, interval 0.609 to 0.778, AUC 0.713. Across seven estimable cells, 0.698 to 0.722, every interval overlapping. [stop]" (Section 7.5.)

**7.5b. How much is age and sex?**

"[firm] Most of it. Age and sex alone reached 0.640 to 0.669 across the four cells. [pause] On identical registry rows CALON-5 beat the two-term model by 0.043, interval 0.005 to 0.083. [concede] No interval exists for CALON-5 over age and sex, so whether that increment is demonstrated cannot be judged. [stop]" (Section 7.5; limitation L14.)

**7.5c. The 4.9 cut-point as a ranker?**

"[slow] Close to chance, concordance 0.552 to 0.593; CALON-5 exceeded it by 0.111 to 0.170 with intervals excluding zero in all four cells. [pause] A severity marker is not a risk ranker. [stop]" (Section 7.5.)

### Section 7.6 Calibration

**7.6a. The headline.**

"[slow] Absolute risk did not transport and failed in opposite directions. Into the registry: 3.05 per cent predicted against 5.98 observed at ten years, observed-to-expected 1.96. Into UK Biobank: 6.76 against 4.73, ratio 0.70. None of seven ratios lay near one. [pause] Slopes 0.83, interval 0.47 to 1.22, into the registry; 0.55, interval 0.36 to 0.72, out of it. [stop]" (Section 7.6.)

**7.6b. What fixed it, and is the fixed model validated?**

"[firm] Re-estimating the baseline fixed the level in all seven cells; the slope needed a second correction only out of the registry. The two-term model needed no update out of the registry and both corrections into it. [pause] Updated results describe adaptation to the receiving cohort. They are not external validation. [stop]" (Section 7.6.)

**7.6c. Why is the gap unsurprising?**

"[slow] By attained age 70, cumulative MACE incidence was 9.0 per cent in UK Biobank carriers, interval 6.4 to 11.9, and 24.4 per cent in registry carriers, interval 19.3 to 30.6. A model inherits its development cohort's baseline hazard. Relative effects travelled; absolute risk did not. [stop]" (Section 7.6.)

### Section 7.7 What the recorded value costs guideline grading

**7.7a. The result, exactly.**

"[slow] In 261 genotype-positive registry patients with both a recorded pre-treatment and an on-treatment LDL-C, all measured, none reconstructed: median 6.1 untreated, 3.6 on treatment. 173 of 261, 66.3 per cent, sat in a lower ESC/EAS stratum on the recorded value; 29.5 unchanged; 4.2 higher. 193 were at or above 4.9 untreated and 57 on treatment; 145 of the 193, 75.1 per cent, would not have been identified on the recorded value. [stop]" (Section 7.7.)

**7.7b. Its weakness?**

"[concede] It is estimated only in patients who had a documented untreated value, and who has one is not random; they may be earlier-identified or more severe. So it is not a service-wide rate. [lift] Its strength is that it needs no reconstruction. [stop]" (Section 7.7; limitation L16.)

**7.7c. What is new relative to Singh 2026 and Mancini 2026?**

"[firm] Singh applied correction factors to trial participants; Mancini imputed untreated LDL-C in homozygous FH. The increment here is a direct measurement, in genotyped heterozygous carriers, of how often the recorded value under-grades severity. [stop]" (Section 7.7; refs 92, 93, DOIs verified 10.1016/j.jacl.2026.03.018 and 10.1016/j.atherosclerosis.2025.120590.)

**7.7d. Goal attainment.**

"[slow] Below 1.8: 7 of 594 treated UK Biobank participants, 1.2 per cent, and 5 of 1,194 treated registry patients, 0.4 per cent. Of the 261 with a pair, 26.8 per cent met the NICE 50 per cent reduction, 5.0 reached below 1.8, 3.1 met both; 88.6 per cent of those achieving the reduction still missed the goal. Attainment depended on the starting value, 39.4 above the median baseline against 12.9 below, as regression to the mean predicts. [firm] No goal contrast was estimable for events. [stop]" (Section 7.7; ref 73, NICE CG71, not DOI-indexed.)

**7.7e. Two measurement properties?**

"[slow] The triglyceride term's median share of Friedewald LDL-C rose from 9 per cent at the severity cut-point to 47 per cent below 1.4, so the lowest goals are judged on the least reliable number. And of 51 registry patients at the 2.6 goal with a pre-treatment value, 54.9 per cent had been at or above 4.9 untreated. [stop]" (Section 7.7.)

**7.7f. Does the cut-point separate risk?**

"[slow] Yes, more sharply in the clinic: event-rate ratio above versus below 4.9 was 2.97 in the registry and 1.61 in UK Biobank. That is the attenuation ascertainment predicts. [stop]" (Section 7.7.)

### Section 7.8 Comparison with published scores

**7.8a. Rules of the comparison.**

"[slow] Each score implemented as published with its own variables, coding and missing-input rules, Lp(a) converted to the units it takes. Compared with CALON-5 on rows where both are computable, one set of family-grouped resamples so every difference is paired. No non-inferiority verdicts, because no margin exists. A comparison resolves only when its interval excludes zero. [stop]" (Section 7.8.)

**7.8b. The results.**

"[slow] Standardised P/LP, ten years: FH-Risk-Score into the registry, 402 rows, 16 events, below the floor, point estimate minus 0.039 favouring the comparator. FH-Risk-Score out, 759 rows, 31 events, plus 0.029, interval minus 0.065 to 0.120. Montreal into the registry, 449 rows, 20 events, plus 0.052, interval plus 0.0003 to 0.104, resolves by the rule and did not under the earlier hypertension definition. Montreal out, 888 rows, 42 events, minus 0.006. SAFEHEART-RE into the registry, 26 rows, one event, not estimable. SAFEHEART-RE out, 725 rows, 31 events, plus 0.076, interval minus 0.041 to 0.189. [firm] CALON-5 higher in 3 of 4 estimable comparisons; one resolved, fragile. [stop]" (Section 7.8, Table 7.2.)

**7.8c. Did the published scores calibrate?**

"[slow] No. On UK Biobank P/LP rows SAFEHEART-RE under-predicted, observed-to-expected 1.72, interval 1.22 to 2.33, slope 0.52; the FH-Risk-Score over-predicted, 0.34, interval 0.23 to 0.47, slope 0.78. [stop]" (Section 7.8; refs 53, 56, DOIs verified 10.1161/CIRCULATIONAHA.116.024541 and 10.1161/ATVBAHA.121.316106.)

**7.8d. How does this sit against McKay and the Australian work?**

"[slow] McKay transported SAFEHEART-RE into English primary care by imputing unrecorded predictors: Harrell's c 0.67, interval 0.61 to 0.72, with substantial miscalibration. In Australian cohorts with measured inputs, SAFEHEART-RE validated for incident events and the Montreal, Combined and FH-Risk scores discriminated prevalent disease. [pause] One cohort, one direction each. Mine is paired, common rows, both directions, and finds most comparisons too small to resolve. [stop]" (Section 7.8; refs 28, 29, 94, 95, DOIs verified.)

### Section 7.9 Risk enhancers

**7.9a. Does apoB add to CALON-5? Reconcile with Chapter 4.**

"[slow] In 888 P/LP carriers with 74 events, no enhancer's interval excluded one: apoB/LDL-C discordance 0.79, interval 0.38 to 1.62; Lp(a) at or above 125 nanomoles per litre 0.88, interval 0.40 to 1.94. In all carriers the polygenic score was 1.37 per SD, interval 1.12 to 1.68, with concordance plus 0.012, interval minus 0.000 to 0.040. [firm] Chapter 4 asks about association at fixed LDL-C; this asks about increment over a model already using the untreated lipid ratio, age, sex and comorbidity, in a smaller frame. Different questions; moderate effects are not excluded. [stop]" (Section 7.9.)

### Section 7.10 Treatment reconstruction and measurement error

**7.10a. Can stable concordance validate the reconstruction?**

"[firm] No. If one correction moves most treated participants similarly, order changes little even when each value is wrong. Ranking stability cannot validate a threshold or a probability. [pause] Cumulative-exposure constructs did not improve on a single concentration: concordance difference 0.007, point estimate only. [stop]" (Section 7.10.)

**7.10b. What would be the stronger design?**

"[slow] Documented pre-treatment samples, serial lipids, time-linked prescribing, dose and adherence; and where untreated values remain unobserved, a distribution of plausible values propagated through fitting and calibration. [stop]" (Section 7.10; Appendix D.4.)

### Section 7.11 Calculability and the decisive study

**7.11a. Why is calculability "a result in its own right"?**

"[slow] SAFEHEART-RE requires Lp(a) and body-mass index and publishes no missing-input rule; after merging the laboratory export Lp(a) was available for 185 standardised registry patients but SAFEHEART-RE was computable for 26, with one event, because BMI was the binding input. The FH-Risk-Score publishes a zero-when-unavailable rule for Lp(a). CALON-5 itself was computable for 469 of 945 registry carriers because diabetes was so often unrecorded. [firm] The thesis reports what a service can compute with the data it holds. [stop]" (Section 7.11, Table 7.3.)

**7.11b. Why not impute SAFEHEART-RE's inputs, as McKay did?**

"[firm] Because the score publishes no missing-input rule, so imputing from derivation means would create an adapted score and label it as the published one. Possible, deliberately not done. [stop]" (Section 7.11; ref 28, DOI verified 10.1016/j.atherosclerosis.2022.07.011.)

**7.11c. The size of a decisive study.**

"[slow] Held at observed differences, comparisons out of the registry would resolve at about 71 events against SAFEHEART-RE, 312 against the FH-Risk-Score, and 6,656 against Montreal. For calibration-in-the-large within plus or minus 0.1, about 362 ten-year events; for the slope, 248 to 1,285. [pause] A genotype-confirmed cohort untouched by either model's development. [stop]" (Section 7.11.)

### Section 7.12 Decision curves

**7.12a. Do they show utility?**

"[slow] Frozen CALON-5 exceeded both treat-all and treat-none at 13 of 30 thresholds into the registry and 10 of 30 out of it; 18 and 7 after recalibration; the cut-point alone at 6 and 4. Net benefit undefined at six thresholds into the registry. [firm] Fewer than half the thresholds: net benefit does not follow from discrimination. [concede] Appendix A records 12 and 9; one is stale and I will correct it. [stop]" (Section 7.12; Appendix 1, item P1; ref 33, DOI verified.)

**7.12b. What is the credible action in FH?**

"[firm] Not a treatment threshold, because lipid lowering is already indicated. Added review, imaging or phenotyping. A low score must never withdraw indicated therapy. [stop]" (Section 7.12.)

**7.12c. Reclassification against the cut-point?**

"[slow] 47.5 per cent of registry and 35.0 per cent of UK Biobank carriers changed category. Frozen in the registry, the model placed 50.1 per cent lower than the cut-point and 0.6 higher; patients who had an event were moved down more often than up, event reclassification minus 0.61, not prespecified. [pause] That is the reclassification face of the under-prediction. [stop]" (Section 7.12.)

### Section 7.13 Effect modification

**7.13a. Sex and age?**

"[slow] Interactions with the lipid term: none of nine estimable excluded one; sex by lipid 0.958, interval 0.716 to 1.282. For the cut-points the ratio of hazard ratios per decade was below one at 12 of 14 thresholds, so cut-points discriminate hardest in the young. Sex modified 2 of 10, consistent with chance. No multiplicity adjustment. [stop]" (Section 7.13.)

### Section 7.14 The CALON-C lineage

**7.14a. What does CALON-C add that CALON-5 does not?**

"[slow] A sharper version of the same warning. Its frozen nine-term Welsh concordance of 0.7252 was matched by age and sex; the other seven terms added plus 0.0008. Its Welsh calibration failed, observed-to-expected 4.500 at five years and 4.374 at ten, repaired only by re-estimating the baseline. [stop]" (Section 7.14; Appendix B.)

**7.14b. Why is it lineage rather than evidence?**

"[firm] The exact historical fitted object and the missing-data pathway are unavailable; endpoint dates were unambiguous for 147 of 289 events and a stress test dropped concordance from 0.708 to 0.660; competing deaths had to be reconciled from 193 to 224. Its absolute-risk, Brier and decision-curve outputs are withheld. [stop]" (Appendix B.4.)

### Section 7.15 Answer

**7.15a. State the answer to the cohort part.**

"[firm] Ranking survived: concordance 0.699 and 0.698 in the two directions. Absolute risk did not: observed-to-expected 1.96 and 0.70, and the published scores miscalibrated the same way. [pause] What restores it is local re-estimation of baseline risk; what would settle it is a cohort with several hundred events. [stop]" (Section 7.15.)

**7.15b. In one sentence, why is neither model clinical?**

"[firm] Because most of the discrimination was demographic, most comparisons were too small to resolve, and neither has calibration, utility or impact evidence in a new setting. [stop]" (Section 7.15.)

---

## Chapter 8: Discussion

### Section 8.1 The answer

**8.1a. Give the answer to the central question in one breath.**

"[firm] An LDL-C value means the same thing only when five conditions also match: the particles carrying it, the route by which the person was found, the treatment already given, the cohort a model was built in and, in women, reproductive stage. [pause] Where any differs, the value needs something added before it can support the decision. [stop]" (Section 8.1.)

**8.1b. Where do the guidelines stand, exactly?**

"[slow] The 2026 ACC/AHA multisociety guideline rates apoB measurement COR 2a, LOE B-NR, in adults on lipid-lowering therapy, and COR 2b in adults not on therapy including to characterise inherited disorders. The 2025 ESC/EAS focused update makes no apoB recommendation and treats LDL-C as a direct cause of ASCVD. [firm] Nothing here challenges that. [stop]" (Section 8.1; refs 2, 3, DOIs verified 10.1161/CIR.0000000000001423 and 10.1093/eurheartj/ehaf190.)

**8.1c. Table 8.1 says the Objective 4 absolute-risk falsifier was "met". Was it?**

"[concede] No, the wording is wrong. The falsifier was calibration that transports without updating. Calibration did not transport, so the falsifier is not met. [pause] And the ranking falsifier, no better than age and sex, cannot be cleanly judged because no interval exists for that increment. I will correct the row. [stop]" (Table 8.1; Appendix 1, item P4.)

**8.1d. Is the answer interventional?**

"[firm] No. No chapter tested a change in care. It establishes where a single recorded value misleads, by how much in these cohorts, and what information would correct each error. [stop]" (Section 8.1.)

### Section 8.2 Principal findings by task

**8.2a. Why the seven-step template in every subsection?**

"[firm] So that finding, literature, explanation with its tier, competing explanation, boundary, implication and resolving study cannot be quietly omitted for any task. [concede] The cost is a mechanical read; the kit review suggests connecting them into argued paragraphs, and I agree. [stop]" (Section 8.2; K012-049.)

**8.2b. For each task, the strongest competing explanation?**

"[slow] Particles: the residual summarises metabolic comorbidity rather than particle number. Route: polygenic enrichment among probands, and era. Treatment: TUDOR's advantage reflects how weak the electronic reconstructions were. Cohort: differences in treatment era, access and event capture between a volunteer cohort and a service, nothing FH-specific. [firm] I accept each as plausible. [stop]" (Sections 8.2.1 to 8.2.4.)

**8.2c. "Relative effects travel, baseline risk does not." Did you test that?**

"[concede] Partially. In the all-carrier fits the five hazard ratios agreed in direction and the lipid term was 2.33 and 2.40 in the two cohorts; but coefficient or interaction transport was not formally tested, and the kit flags the sentence as an inference beyond the test. [firm] The supported statement is that ranking was comparable while baseline incidence differed by attained age 70, 9.0 against 24.4 per cent. [stop]" (Sections 7.3, 7.6, 8.2.4; K012-034.)

### Section 8.3 One value, three uses

**8.3a. The three uses and the condition that breaks each.**

"[slow] Grading severity: broken by treatment; restored by a documented untreated value. Finding carriers: weakened by treatment and by route; partly restored by treatment-aware reconstruction with genetic testing, locally until tested prospectively. Ranking risk: ordering survived a change of cohort, absolute risk did not; restored by local baseline re-estimation and, in women, sex-specific calibration. [stop]" (Section 8.3.)

**8.3b. You say the errors are "systematic and predictable in direction". Is that an overclaim?**

"[concede] For individuals, yes, and the kit flags it. At group level the directions observed were consistent: under-grading after treatment, milder phenotypes in cascade and genotype-first carriers, and miscalibration on transfer. [firm] Random error and unvalidated individual remedies remain; I will restrict the sentence to group-level patterns. [stop]" (Section 8.3; K012-051.)

**8.3c. Is 4.9 mmol/L an "FH severity" cut-point?**

"[concede] The guidelines use 4.9 as a severe-hypercholesterolaemia and diagnostic-suspicion threshold and, in the ESC/EAS bands, as a stratum boundary; whether it is properly called FH severity grading needs the exact guideline tables. [firm] The empirical result stands regardless: 75.1 per cent of registry carriers above it untreated were not above it on the record. I will reword the category. [stop]" (Sections 1.1, 7.7, 8.3; K012-007.)

**8.3d. "The laboratory value was accurate." How do you know?**

"[concede] I do not know analytically; the sentence asserts more than the evidence. The claim I can support is that the value changed its meaning, not that its assay accuracy was established. I will remove the accuracy assertion. [stop]" (Section 1.1; K012-008.)

**8.3e. "Linked genotyped cohorts with common definitions." Are the definitions common?**

"[concede] Only partly, and the kit is right to press. The two data sources are common; the genetic definitions differ by chapter, and participant overlap between Chapters 4 and 5, and between UK Biobank frames, is not established. [firm] The honest description is linked heterogeneous studies in two cohorts with a common evidence rule, and I will build the overlap and definition matrix. [stop]" (Sections 1.6, 8.3; K012-013.)

### Section 8.4 Beyond FH: sex and reproductive stage

**8.4a. State the finding and its design.**

"[slow] 273,036 UK Biobank women across the menopausal transition, age-matched male reference. At age 58 an apoB of at least 1.2 grams per litre corresponded to a ten-year MACE risk of 10.88 per cent in men and 4.34 in postmenopausal women. A sex-neutral internal model over-predicted in women, observed-to-expected 0.65 to 0.78, and under-predicted in men, 1.19. [pause] Population cohort, self-reported status, cross-sectional stages, internal model. [stop]" (Section 8.4; ref 50, in press, not indexed.)

**8.4b. Isn't that just lower baseline hazard in women at 58?**

"[firm] Yes, and that is the point being made: a marker's meaning for absolute risk depends on the baseline it is read against. [concede] It is a familiar age-sex baseline difference, not evidence of a sex-specific lipid effect or a menopause mechanism, and the section must not be read as either. [stop]" (Section 8.4; K012-005, 053.)

**8.4c. Da Roza found a larger rise in monogenic FH; you found a smaller percentage rise in LDLR carriers. Contradiction?**

"[slow] Different metrics from different starting levels: an absolute rise in a small FH cohort against a percentage rise, plus 11.5 versus plus 15.3, from a higher base in a population cohort. [firm] Not directly comparable, potentially compatible. I will not claim more than that. [stop]" (Section 8.4; ref 24, DOI verified 10.1016/j.atherosclerosis.2025.120587; K012-052.)

**8.4d. Has sex-specific calibration been validated for FH?**

"[firm] No. The recommendation to use sex-specific absolute risk in interpretation comes from a population cohort. The FH resolving study would need within-person, final-menstrual-period-anchored lipids in molecularly characterised women. [stop]" (Section 8.4; K012-048.)

### Section 8.5 What the thesis adds

**8.5a. Two layers of contribution?**

"[slow] What the studies established individually: apoB carries information at fixed LDL-C in genotyped carriers; treatment-aware triage ranks better than electronic reconstructions in one Welsh comparison; the route contrast at matched residues; ranking transports while calibration does not; the same cut-point differs by sex. [pause] What the thesis adds: the five measured together with one evidence rule, the identification of the one condition under which the value does transport, and a named repair for each. [stop]" (Section 8.5.)

**8.5b. Is claim-to-evidence governance a contribution?**

"[concede] It is good practice, not a new method, and the kit review says so. I present it as a strength that supports trust, not as doctoral novelty. [stop]" (Section 8.5; K012-015.)

**8.5c. What would refute the framework?**

"[slow] A frozen residual failing under harmonised assays in an independent FH cohort. No residue-matched route difference in an independently designed cohort. Frozen TUDOR failing against fully scored DLCN on complete data. A frozen prognostic model calibrating across independently ascertained cohorts without updating. [stop]" (Section 8.5.)

### Section 8.6 Implications

**8.6a. What changes tomorrow?**

"[firm] One recording practice: record the first untreated LDL-C with date and method as a fixed field and grade on it; where none exists, say so. [pause] Everything else is better understanding, more accurate communication of scope, or service evaluation to be done. No change in care has been tested. [stop]" (Section 8.6.)

**8.6b. Should risk scores ever reduce FH treatment?**

"[firm] No. Lipid lowering is already indicated. UK guidance advises against general-population tools in FH, and the 2026 ACC/AHA guideline gives FH-specific scores only COR 2b for short-term risk. [pause] Until a model passes calibration, competing-risk, utility, implementation and impact tests, these models support research and benchmarking. [stop]" (Section 8.6; refs 2, 73.)

**8.6c. What could a within-FH model legitimately inform?**

"[concede] The thesis says intensity, timing in young relatives and allocation of add-on therapy or imaging. The kit asks me to check exact guideline support for each example and to state plainly that a low developmental score cannot delay recommended care. I accept both. [stop]" (Section 7.1; K012-024.)

**8.6d. Minimum data set?**

"[slow] Untreated LDL-C, lipoprotein(a), body-mass index, diabetes and hypertension status, and dated events. Without them, SAFEHEART-RE was computable for 26 standardised registry patients. [stop]" (Section 8.6.)

### Section 8.7 Limitations

**8.7a. Rank the top five.**

"[slow] One: no untouched evaluation cohort, both settings in the UK. Two: few events and informative pedigrees, a registry fit on 28 events. Three: treatment reconstruction as a model-based surrogate with non-random reconstructed rows. Four: the grading result estimated only in patients with a documented untreated value. Five: registry event times from recorded ages with heavy censoring. [stop]" (Table 8.3.)

**8.7b. You write that reconstructed rows are "not missing at random" and that missing diabetes makes MAR "doubtful". Is that the right inference?**

"[concede] Not as written. Observing that missingness is associated with the outcome shows selection; it does not by itself establish MNAR, and it is not a proof that imputation would have failed. [firm] The defensible statement is: missingness was outcome-associated, I chose not to impute, and the complete-case set is selected. I will correct the conditional inference in Sections 5.2, 7.4 and 8.7. [stop]" (K012-029.)

**8.7c. Which limitations are resolved by objects that exist?**

"[slow] L13 is described as resolved by the released two-term lock, and L05 and L06 as superseded by L15. [concede] Table B.1 says no two-term frozen object was located, so the lock's existence and usability must be verified before L13 is called resolved. [stop]" (Section 8.7; Appendix 1, item P3; K012-054.)

### Section 8.8 Resolving studies

**8.8a. The principal study and its size.**

"[slow] Frozen validation of CALON-5 and the two-term model in an independently assembled genotype-confirmed cohort, for example an EAS FHSC registry, with harmonised MACE components and death as a competing event; about 362 ten-year events for calibration-in-the-large within plus or minus 0.1, 248 to 1,285 for the slope. [concede] Those numbers are extrapolations from observed differences and should be labelled illustrative; a future design should be built afresh with target incidence, censoring and desired precision. [stop]" (Section 8.8; K012-044.)

**8.8b. The four secondary studies.**

"[slow] SAIL linkage of the registry to primary-care lipids and prescribing, to measure under-grading service-wide. A cohort recording route at eligibility with exact-allele matching and family-aware estimation. Silent prospective evaluation of frozen TUDOR against fully scored DLCN and the local criteria. And a prospectively locked apoB-residual evaluation in adjudicated FH, head to head with risk-weighted apoB. [stop]" (Section 8.8, Table 8.4; refs 77, 90, DOIs verified.)

**8.8c. Three rules for all of them?**

"[firm] Ascertainment recorded at eligibility as a design variable. Frozen performance characterised before any updating, each updated model with a new identifier. And no biomarker, score or model adopted because it raises an AUC. [stop]" (Section 8.8.)

### Section 8.9 Conclusion

**8.9a. Your closing sentence, if asked for one.**

"[slow] LDL-C remains clinically indispensable. [pause] What this thesis asks is that it never travel alone: every inference drawn from it, or from a model built on it, should carry the conditions that give the number its meaning. [stop]" (Section 8.9.)

---

## Appendices A to G

**A1. What is the admission rule?**

"[firm] A numerical or methodological statement carries weight only when traceable to a publication, manuscript, aggregate package, producing script or explicit author decision; otherwise it is held or omitted. [pause] Tiers T1 to T6; retained aggregates support only sensitivity or descriptive statements. [stop]" (Appendix A.1.)

**A2. Which Appendix A rows are T5 or T6?**

"[slow] T5: the competing-risk sensitivity; the CALON-F grey-zone apoB test; the 71-variant synthesis; the CALON-5 increment over age and sex; the comparator tally; the decision curves; the CALON-C lineage rows. T6: TUDOR family-aware optimism, deleted because unresolvable. [stop]" (Table A.2.)

**B1. Why is Appendix B so long?**

"[firm] Because three models and four manuscript versions share parts of the programme, and no value may be read against the wrong object. It records CALON-C's specification, internal validation, endpoint-date sensitivity, competing-death reconciliation, Welsh transport, comparator benchmarking, reconstruction convention and risk of bias; the Paper 14 version history; and the Welsh denominator reconciliation. [stop]" (Appendix B.)

**B2. What did the CALON-C manuscript claim that the thesis withholds?**

"[slow] Holm-adjusted full-follow-up comparator wins: plus 0.070 against SAFEHEART-RE, plus 0.032 against Montreal, plus 0.015 against the FH-Risk-Score described as a tie; and a post hoc non-inferiority margin. [firm] The horizon-specific thesis analyses supersede them, and the margin analysis is held. [stop]" (Appendix B.6.)

**B3. The Paper 14 version history: read it to me.**

"[slow] As submitted: primary plus 1.11, within-pedigree plus 0.34, interval 0.01 to 0.67. Intermediate: plus 1.23; plus 0.34, interval minus 0.07 to 0.76. R3 v5: plus 1.10; plus 0.40, interval minus 0.003 to 0.81, crossing zero. R4: plus 1.09; plus 0.39, bootstrap 0.04 to 0.74, dummy-variable correction minus 0.06 to 0.85. [pause] Had R3 stood, Objective 2's falsifier would have been met at that rung. [stop]" (Appendix B.10.)

**C1. What is the denominator ontology?**

"[firm] A classification that keeps data-resource totals, chapter analysis sets, complete-case subsets and comparator-evaluable subsets as different classes, so that a comparator subset is never presented as a cohort total. [stop]" (Appendix C.)

**D1. Summarise the reconstruction agreement evidence.**

"[slow] Convention A, 133 pairs: at a residual fraction of 0.65, bias plus 0.09, limits of agreement minus 3.00 to plus 3.19, MAE 1.24, Lin's concordance 0.593; at 0.50, bias plus 1.79. Person-specific fractions had median 0.69 with a CV of about 25 per cent. TUDOR's 649-pair resource: MAE 1.20, correlation 0.32. Convention B: a scalar sweep 0.65 to 0.75 moved Harrell's C by 0.0022 and 0.0019. Convention C: MAE 1.25, R-squared 0.0975. [firm] None validates an individual value. [stop]" (Appendix D.)

**D2. The claim that error worsened in the high-risk tail?**

"[firm] Retired. Against true untreated LDL-C the slopes were negative; the earlier positive slope had placed the reconstruction itself on the horizontal axis. [stop]" (Appendix D.1.)

**E1. Which reporting standards, and what does "addressed" mean?**

"[slow] STROBE and RECORD for Chapters 4 and 5; TRIPOD+AI and PROBAST+AI for 6 and 7; TRIPOD-Cluster for family-structured data. [pause] Addressed means information or an explicit limitation is present; design-level bias remains. [stop]" (Appendix E; refs 30, 31, 36, 37, 49, DOIs verified.)

**F1. What would complete the search?**

"[slow] Embase through Cardiff University Library; a refresh of Elicit, Consensus and Scite after 1 October 2026 with the logged strings verbatim; the eight unread priority full texts; and independent duplicate screening, after which each novelty sentence is re-stated. [concede] The kit adds that historical access statuses must not be overwritten and that Crossref verification is not full-text appraisal. [stop]" (Appendix F.9; K012-055.)

**G1. Your glossary defines calibration as O/E and slope. Is that complete?**

"[concede] It is the summary used, not the concept. Calibration also includes flexible curves and local agreement, and the glossary should separate model state, development independence, O/E versus broader calibration, and publication state versus risk of bias. [stop]" (Appendix G; K012-056.)

---

## Cross-cutting questions

### Integrity

**X1. If you had to disclose one thing to the examiners today that they might otherwise miss, what is it?**

"[long pause] [firm] That the published description of TUDOR includes an ascertainment term the evaluated object does not contain, that no corrigendum has yet been requested, and that I will request one before final submission. [stop]" (Box 6.1.)

**X2. Your own abstracts describe "dual external validation". Your thesis says cross-cohort evaluation. Which is right?**

"[firm] The thesis. The abstracts were written under a looser wording rule than the thesis applies, and they will be corrected to match, not the other way round. [stop]" (Section 6.2; ref 82.)

**X3. How much of this thesis was written by a language model?**

"[firm] Drafting, editing and consistency checking were assisted; every passage was reviewed by me; no analysis was run by a model and no model saw participant data. [pause] Title screening in the evidence search used a model first pass, single-reviewer, not duplicated, and that is listed as a limit with its fix. [stop]" (AI statement; Appendix F.8.)

**X4. Is the published TUDOR AUC 0.760 or 0.7585?**

"[firm] The published paper reports 0.760 in 1,274 participants with 311 carriers, and the thesis uses 0.760. An earlier internal lineage carried 0.7585; it is not the published figure and I do not quote it. [stop]" (Section 6.3; ref 21, DOI verified.)

### Causality and design

**X5. Is anything in this thesis causal?**

"[firm] No. Chapter 4 is association at fixed LDL-C; Chapter 5 is selection into observation; Chapters 6 and 7 are model tasks. The only causal premise, that LDL causes ASCVD, is cited to the EAS consensus and the guidelines. [stop]" (Section 3.2; ref 1, DOI verified 10.1093/eurheartj/ehx144.)

**X6. Name the collider in Chapter 5.**

"[slow] Inclusion in the analysed sample, conditioned on referral. Prior disease can prompt referral, so proband outcome analyses are collider-prone, which is why the event analyses are secondary and why probands with prior ASCVD were excluded in a sensitivity: plus 0.96. [stop]" (Sections 5.1, 5.4, 5.5.)

**X7. What is the target trial Chapter 7 would emulate, and why can't it?**

"[firm] Eligibility: adults with a P/LP variant, free of ASCVD. Strategy: none, it is prognostic, not interventional. Time zero: baseline. Outcome: first MACE with competing death. Estimand: ten-year cumulative incidence. [pause] It cannot emulate a trial because there is no treatment contrast, and it says so: the models rank, they do not estimate effects. [stop]" (Section 7.1.)

**X8. Where is immortal time a risk?**

"[slow] In the registry, where event times derive from recorded ages and baseline is a dated first visit; participants recorded as treated without a start date were excluded from the paired first-versus-later analysis for that reason, and predictors after the outcome are barred by the temporal rule in Section 3.4.3. [stop]"

### Translation

**X9. If a lipid clinic asked what to do with this thesis on Monday, what would you say?**

"[firm] Record the first untreated LDL-C with its date and method, and grade on it. Do not adopt the 0.31 ratio. Do not run TUDOR or CALON-5 on patients. Measure apoB where the question is particle burden, which guidance already supports in treated adults. And collect the minimum data set so that any future FH score can be computed. [stop]" (Section 8.6.)

**X10. A cascade-detected 30-year-old with LDL-C 3.8 and a confirmed familial variant asks if she really has FH.**

"[slow] Yes. Her value is the phenotype observed through a different route, at a younger age, and 84.5 per cent of untreated P/LP carriers in a population cohort were below 4.9. [pause] It does not overturn the molecular diagnosis or imply low cumulative burden. Earlier detection is the opportunity to prevent cholesterol-years. [stop]" (Sections 5.6, 5.7.)

### The candidate

**X11. What would you do differently?**

"[slow] Register the framework and falsifiers prospectively. Fix one reconstruction convention with uncertainty propagation. Record route at eligibility. Set aside a validation cohort before any model choice. Request the corrigendum the day the discrepancy was found. Complete the Embase search. [stop]"

**X12. What are you least confident about?**

"[concede] The within-pedigree result, because it depends on the variance convention and has moved across versions; and the CALON-5 increment over age and sex, because there is no interval. [stop]"

**X13. What are you most confident about?**

"[firm] The grading result in Section 7.7, because it needs no model, no reconstruction and no UK Biobank data: 145 of 193 severe carriers were not severe on the value in the record. [stop]"

---

## Kit 012 flagged issues as viva questions

The paragraph-audit kit (PARAGRAPH_TASKS.csv, 3,146 rows) flags 19 P0, 60 P1 and 9 P2 paragraphs under issue codes K012-001 to K012-056. Each code is an examiner objection waiting to be asked. The table gives the code, the sections it touches, the question it generates and the answer's core. Where the answer already appears above, the question number is given.

| Code | Sections | Examiner question | Answer core |
|---|---|---|---|
| K012-001 | Summary, 1.5, 8.3 | Where is the crosswalk from each condition to population, lipid quantity, estimand and clinical use, with untested cells shown? | Concede it is implicit in Tables 1.1 and 3.1; commit to an explicit crosswalk with untested cells. |
| K012-002 | Summary | Why does the Summary drop the interval and call them "carriers"? | Restore 1.070 to 1.718 and keep "rare LDLR/APOB variant" wording. |
| K012-003 | Summary, 5.3, 8.2.2 | Which interval convention is preferred for 0.39, and what is the +0.47 family-model contrast the kit mentions? | State the preferred estimate and both conventions; separate the route association from any within-family model contrast (5.3c). |
| K012-004 | Summary | Is TUDOR's whole gain due to reconstruction? | No; architecture-level performance, no ablation (6.4d; 8.2.3). |
| K012-005, 052, 053 | Summary, 8.4 | Baseline hazard or marker effect? | Baseline hazard; not an FH result (8.4b, 8.4c). |
| K012-006 | Summary, 8.9 | Which of the five "additions" are demonstrated, which are tools, which are proposals? | Demonstrated interpretation: particle information, route contrast, under-grading. Partially supported tool: TUDOR, local recalibration. Proposed: sex-specific calibration in FH. |
| K012-007, 008 | 1.1, 7.7, 8.3 | Is 4.9 an FH-severity grade, and was the assay "accurate"? | Reword category; remove accuracy assertion (8.3c, 8.3d). |
| K012-009 | 1.1, 1.5 | Universal interchangeability or limits of interpretation? | The aim is the limits; reword "establish" to "evaluate". |
| K012-010 | 1.2, 8.2.1 | Is the particle mechanism measured or assumed? | Plausible background, not measured mediation (4.7a, 4.7b). |
| K012-011 | 1.2, 1.5 | Does the rare-variant frame answer the HeFH question? | Part of it; no population upgrade at the objective (4.2e, 4.9). |
| K012-012 | 1.5 | Are falsifiers protection or stress tests? | Interpretive stress tests written retrospectively (1.5a). |
| K012-013 | 1.6, 8.3 | Common definitions or linked heterogeneous studies? | The latter; build the overlap matrix (8.3e). |
| K012-014 | 1.6, F.10 | Was the joint-novelty claim tested adversarially for combinations? | Only in PubMed and Crossref; narrow to the supported combination (2.1a). |
| K012-015 | 1.6 | Is governance a method? | Good practice, not novelty (8.5b). |
| K012-016 | 1.7 | Define each genetic label positively. | Align with Box 3.1 (1.7b). |
| K012-017 | 2.5 | Ratios, components and residuals compared under explicit assumptions? | Retain the qualification; compare all three where data allow (A04). |
| K012-018 | 2.8 | Is Akyea 2026 the "only" specialist comparison? | Only one identified; search further before "only" (6.4b). |
| K012-019 | 5.3 | Lead with magnitude and precision. | Done in 5.3a to 5.3d. |
| K012-020 | 5.4 | Is the non-HDL excess a remnant mechanism? | Timing explains most; mechanism not shown (5.4b). |
| K012-021 | 5.4 | THESIS-AN-02's output was not located. Why is the number in the thesis? | Demote or hold; do not substitute the 53-variant rerun (5.4c). |
| K012-022 | 5.4 | Is the first-versus-later reading a mechanism? | A paired phenotype observation, not a mechanism (5.4e). |
| K012-023 | 6.2 | Is the object locked with metadata and readback? | Yes within tolerance; no Index Effect added from text (6.2a, 6.2b). |
| K012-024 | 7.1 | Guideline support for intensity, timing, add-on examples? | Check each; state a low score cannot delay care (8.6c). |
| K012-025 | 7.1 | Is r7 "submitted" evidenced? | Use "manuscript revision r7" until receipt is evidenced. |
| K012-026 | 7.2 | Blank medication field read as "no" for 436 of 945. Justified? | Recording convention with dated lipid visit; provenance retained (7.2a). |
| K012-027 | 7.2 | Was the reconstruction cross-validation person-grouped and carrier-trained? | Verify; do not conflate CV with external agreement (7.2c). |
| K012-028 | 7.3 | What was nested in the bootstrap? | The five-term refit; not the term choice; hold independent-validation language (7.3b, 7.3c). |
| K012-029 | 5.2, 7.4, 8.7 | Is outcome-associated missingness proof of MNAR? | No; correct the inference (8.7b). |
| K012-030 | 7.4 | Endpoint rounding, ties, censoring positivity? | Qualify unstable ten-year estimates (7.4c). |
| K012-031 | 7.5 | "Every interval overlapping" is not equivalence. | Agreed; report range and uncertainty only (7.5a). |
| K012-032 | 7.5 | Is beating the 4.9 cut-point a fair comparison? | It tests a misapplication, not a legitimate model; the fair baseline is age and sex (7.5b, 7.5c). |
| K012-033 | 7.6 | Was calibration "restored" by a ratio rule? | Level fixed on held-out rows; spread and local calibration separate (7.6b). |
| K012-034 | 7.6, 8.2.4 | Did relative effects travel? | Partially shown; not formally tested (8.2c). |
| K012-035, 050 | 7.7, 8.2.3 | Does any reconstruction restore individual severity? | No; recorded pretreatment retrieval is the simpler remedy (7.7b). |
| K012-036 | 7.7 | Goal eligibility, timing, mathematical coupling? | Verify denominators; discuss coupling and room for reduction (7.7d). |
| K012-037 | 7.7 | "Less penetrant" by what definition? | Descriptive; follow-up, treatment, age and endpoint differ (7.7f). |
| K012-038 | 7.7 | Is the TG-term share an error metric? | No; needs reference-method comparison (7.7e). |
| K012-039, 040 | 7.8 | Comparator equations, units, exact lower limit; naive pooling? | Retain +0.0003 and the hypertension sensitivity; pooling exploratory (7.8b). |
| K012-041 | 7.8 | Are the Australian validations being minimised? | State their layers exactly (7.8d). |
| K012-042 | 7.9 | Different marker definitions across chapters? | Ratio here, residual there; moderate effects not excluded (7.9a). |
| K012-043 | 7.11 | Is the BMI bottleneck verified, and is research imputation forbidden? | Verify; not forbidden, deliberately not done (7.11a, 7.11b). |
| K012-044 | 7.11, 8.2.4, 8.8 | What does 362 assume? | Observed-effect extrapolation; label illustrative (8.8a). |
| K012-045 | 7.12 | Undefined DCA cells: estimator or mathematics? | Inspect; preserve as missing (7.12a). |
| K012-046 | 7.12 | Reclassification categories and action equivalence? | Exploratory NRI; no causal harm claim (7.12c). |
| K012-047 | 7.13 | Expected nominal positives from the actual test count? | Compute from tests and alpha (7.13a). |
| K012-048 | 7.15, 8.3 | Sex-specific calibration validated for FH? | No (8.4d). |
| K012-049 | 8.2 | Connect the seven steps into argument. | Agreed (8.2a). |
| K012-051 | 8.3 | Errors "predictable in direction"? | Group-level only (8.3b). |
| K012-054 | 8.7 | Is the two-term lock located and usable? | Verify against Table B.1 (8.7c). |
| K012-055 | F.10 | Five-platform session logged? | Not yet; do not overwrite historical statuses (F1). |
| K012-056 | Glossary | Separate the four concepts. | Agreed (G1). |

**Kit novelty plan, N01 to N08.** The sceptical readings in the kit's novelty plan are the same objections in condensed form: familiar principles assembled after results (N01); tautological treatment effect (N02); expected proband selection (N03); established residual and broad genetic frame (N04); input-deprived comparators (N05); non-independent development and sparse events (N06); familiar age-sex baseline (N07); reporting discipline is not a method (N08). The defensible increments listed against each are the answers given in Sections 1.6, 4.10, 5.7, 6.4, 7.15, 8.4 and 8.5 above.

**Kit analysis cards, A01 to A15.** If asked "what would you do next with the data you have", the essential cards are A01 (reproduce the locked objects by readback), A02 (fixed-membership attenuation with family uncertainty), A03 (missingness wording and mapping), A05 (paired classification consequences with uncertainty), A09 (endpoint and censoring sensitivity), A10 (calibration and decision-curve estimator checks) and A11 (a target-specific sample-size card). A13 to A15 are new-data studies and should be named as such.

---

## Appendix 1: corrections to make before the meeting

Items P1 to P13 are from the paragraph-by-paragraph read; P14 to P19 from the citation check; K-items from the kit.

1. **P1.** Decision-curve counts: Section 7.12 says 13 of 30 and 10 of 30; Table A.2 says 12 and 9. Align to r7.
2. **P2.** Comparator tally: Section 7.8 says 3 of 4 and one resolved; Table A.2 says 4 of 4 and none resolved. Table 7.2 supports the chapter. Correct the ledger.
3. **P3.** Two-term lock: Sections 3.4.1, 7.2, Table 7.1 say released; Table B.1 says not located. Reconcile (K012-054).
4. **P4.** Table 8.1 Objective 4 falsifier status reads backwards. Correct; note the ranking falsifier cannot be judged without an interval (L14).
5. **P5.** Reference 52 status: update to R4 under revision.
6. **P6.** Objective 2 falsifier wording differs between Table 1.2 and Section 5.3.
7. **P7.** Chapter 4: state the prevalent-ASCVD exclusion rule and count for the 1,461 frame.
8. **P8.** Chapter 7: enumerate genes and per-gene counts.
9. **P9.** Define the UK Biobank ancestry restriction and confirm consistency across chapters.
10. **P10.** Reference 77 cited for two different things (Sections 2.8 and 8.8); check or split.
11. **P11.** Bring the three reconstruction-agreement figures (r 0.09; MAE 1.20 and r 0.32; MAE 1.25 and R-squared 0.0975) into one table with one sentence on why they differ; reconcile with the 4.5 mmol/L bias statement in Section 7.2.
12. **P12.** Complete TUDOR PROBAST+AI or say why not (Table 6.5 versus Appendix E.4).
13. **P13.** Title leads with particle burden; prepare a defence or a subtitle.
14. **P14.** Tybjaerg-Hansen 2005: total cholesterol, not LDL-C increments (Sections 1.3, 5.7, Table 2.4).
15. **P15.** Jin 2023: adjustment for VLDL, LDL and HDL particle concentrations (Sections 2.3, 4.7).
16. **P16.** Marston 2022: the abstract names non-HDL-C and triglycerides, not LDL-C (Section 2.3).
17. **P17.** Bogsrud 2025: contrast within 61 FH newborns, not 113 (Section 5.2).
18. **P18.** "Index Effect" does not appear in the TUDOR title or abstract; quote the abstract's wording (Section 6.2).
19. **P19.** Reference 82 not indexed; verify title and the 4,028 figure against the abstract PDF.
20. **K012-007, 008.** Reword the 4.9 category and remove "the laboratory value was accurate" (Sections 1.1, 7.7, 8.3).
21. **K012-013.** Replace "common definitions" with an overlap and definition matrix (Sections 1.6, 8.3).
22. **K012-021.** Demote or hold the 71-variant synthesis until its output is located (Section 5.4).
23. **K012-029.** Correct the MAR/MNAR inference (Sections 5.2, 7.4, 8.7, Table 8.3).
24. **K012-038.** Do not present the TG-term share as an error metric (Section 7.7).
25. **K012-044.** Label the 362-event figure illustrative and specify its assumptions (Sections 7.11, 8.2.4, 8.8).
26. **K012-045.** Inspect the undefined decision-curve cells and preserve them as missing (Section 7.12).
27. **K012-034, 051.** Restrict "relative effects travel" and "errors predictable in direction" to what was tested (Sections 7.6, 8.2.4, 8.3).
28. **K012-025.** Use "manuscript revision r7" until submission is evidenced (Section 7.1).
29. **K012-002, 003.** Restore the interval in the Summary; state the preferred within-pedigree convention (Summary, Section 5.3).
30. **K012-049, 056.** Connect the Section 8.2 steps into argument; expand the glossary definitions.

---

## Appendix 2: DOI register (PubMed-verified during preparation)

Verified means the DOI resolved to a PubMed record with the cited authors, journal and year during preparation on 27 September 2026. Not resolved means no PubMed record was found by DOI, citation lookup or search; quote these only as "the thesis cites".

| Ref | First author, year | DOI | Status |
|---|---|---|---|
| 1 | Ference 2017 | 10.1093/eurheartj/ehx144 | Verified |
| 2 | ACC/AHA 2026 | 10.1161/CIR.0000000000001423 | Verified |
| 3 | Mach 2025 | 10.1093/eurheartj/ehaf190 | Verified |
| 4 | Gratton 2023 | 10.1161/ATVBAHA.123.319438 | Verified |
| 5 | Hu 2020 | 10.1161/CIRCULATIONAHA.119.044795 | Verified |
| 6 | Fry 2017 | 10.1093/aje/kwx246 | Verified |
| 7 | van Alten 2024 | 10.1093/ije/dyae054 | Verified |
| 8 | Schoeler 2023 | 10.1038/s41562-023-01579-9 | Verified |
| 9 | Zhang 2021 | 10.1001/jamacardio.2021.3508 | Verified |
| 10 | Sniderman 2011 | 10.1161/CIRCOUTCOMES.110.959247 | Verified |
| 11 | Richardson 2020 | 10.1371/journal.pmed.1003062 | Verified |
| 12 | Sniderman 2024 | 10.1093/eurheartj/ehae258 | Verified |
| 13 | Marston 2022 | 10.1001/jamacardio.2021.5083 | Verified |
| 14 | Johannesen 2021 | 10.1016/j.jacc.2021.01.027 | Verified |
| 15 | Sayed 2024 | 10.1001/jamacardio.2024.1310 | Verified |
| 16 | ERFC 2009 | 10.1001/jama.2009.1619 | Verified |
| 17 | Sniderman 2013 | 10.1016/j.jacl.2013.08.004 | Not resolved by the tools used; thesis lists PMID 24314360 |
| 18 | Helgadottir 2022 | 10.1093/eurjpc/zwac219 | Verified |
| 19 | Archie 1981 | 10.1097/00000658-198103000-00008 | Verified |
| 20 | Lolli 2019 | 10.1136/bjsports-2017-098110 | Verified |
| 21 | Genedy 2026 (TUDOR) | 10.1016/j.jacl.2026.06.030 | Verified (PMID 42601321) |
| 22 | Genedy (CALON-C manuscript) | none | Not a publication |
| 23 | Soran 2020 | 10.1097/MOL.0000000000000692 | Verified |
| 24 | Da Roza 2026 | 10.1016/j.atherosclerosis.2025.120587 | Verified |
| 25 | Johansen 2023 | 10.1016/j.athplu.2023.01.001 | Verified |
| 26 | Klevmoen 2021 | 10.1016/j.atherosclerosis.2021.09.003 | Verified |
| 27 | Ransohoff 1978 | 10.1056/NEJM197810262991705 | Verified |
| 28 | McKay 2022 | 10.1016/j.atherosclerosis.2022.07.011 | Verified |
| 29 | Tamehri Zadeh 2026 | 10.1016/j.atherosclerosis.2026.120799 | Verified |
| 30 | Collins 2024 | 10.1136/bmj-2023-078378 | Verified |
| 31 | Moons 2025 | 10.1136/bmj-2024-082505 | Verified |
| 32 | Riley 2024 | 10.1136/bmj-2023-074820 | Verified |
| 33 | Vickers 2019 | 10.1186/s41512-019-0064-7 | Verified |
| 34 | Austin 2016 | 10.1161/CIRCULATIONAHA.115.017719 | Verified |
| 35 | van Geloven 2022 | 10.1136/bmj-2021-069249 | Verified |
| 36 | Benchimol 2015 | 10.1371/journal.pmed.1001885 | Verified |
| 37 | von Elm 2008 | 10.1016/j.jclinepi.2007.11.008 | Verified |
| 38 | Johannesen 2024 | 10.1016/j.jacc.2024.03.423 | Verified |
| 39 | Pan 2026 | 10.1186/s12944-026-02928-z | Verified |
| 40 | Mourre 2025 | 10.1093/eurjpc/zwaf234 | Verified |
| 41 | Trinder 2024 | 10.1161/ATVBAHA.123.320287 | Verified |
| 42 | Besseling 2017 | 10.1093/eurheartj/ehw135 | Verified |
| 43 | Weng 2019 | 10.1016/S2468-2667(19)30061-1 | Verified |
| 44 | Haralambos 2015 | 10.1016/j.atherosclerosis.2015.03.003 | Verified |
| 45 | Akyea 2026 | 10.1016/j.jacl.2026.07.009 | Verified |
| 46 | Ramos 2020 | 10.1016/j.atherosclerosis.2019.10.016 | Verified |
| 47 | Trinder 2020 | 10.1016/j.jacc.2020.03.065 | Verified |
| 48 | Van Calster 2019 | 10.1186/s12916-019-1466-7 | Verified |
| 49 | Debray 2023 | 10.1136/bmj-2022-071018 | Verified |
| 50 | Genedy 2026 (midlife) | 10.1016/j.jacl.2026.07.022 | Not resolved (in press) |
| 51 | Genedy 2026 (discordance) | 10.1016/j.jacl.2025.11.008 | Verified (PMID 41617625) |
| 52 | Genedy (Paper 14) | none | Manuscript under revision |
| 53 | Perez de Isla 2017 | 10.1161/CIRCULATIONAHA.116.024541 | Verified |
| 54 | Paquette 2017a | 10.1016/j.jacl.2016.10.004 | Verified |
| 55 | Paquette 2017b | 10.1016/j.jacl.2017.07.008 | Verified |
| 56 | Paquette 2021 | 10.1161/ATVBAHA.121.316106 | Verified |
| 57 | Ahmad 2026 | 10.1016/j.jacl.2026.01.011 | Verified |
| 58 | Gallo 2020 | 10.1016/j.atherosclerosis.2020.06.011 | Verified |
| 59 | Paquette 2025a | 10.1016/j.jacl.2025.01.004 | Verified |
| 60 | Trinder 2020 (JAMA Cardiol) | 10.1001/jamacardio.2019.5954 | Verified |
| 61 | Jin 2023 | 10.1161/JAHA.123.029552 | Verified |
| 62 | Morze 2025 | 10.1093/eurheartj/ehaf207 | Verified |
| 63 | Du 2025 | 10.1093/eurjpc/zwaf750 | Verified |
| 64 | Paquette 2025b | 10.1016/j.jacl.2025.08.009 | Verified |
| 65 | Lyu 2026 | 10.1186/s12944-025-02852-8 | Verified |
| 66 | Genedy 2025 (Atheroscler Plus abstract) | 10.1016/j.athplu.2025.10.006 | Not resolved |
| 67 | Kim 2026 | 10.1016/j.jacl.2026.07.015 | Verified |
| 68 | Rehman 2026 | 10.1093/eurheartj/ehaf1124 | Verified |
| 69 | Drexel 2021 | 10.1016/j.atherosclerosis.2021.05.010 | Verified |
| 70 | Xiao 2023 | 10.1186/s12944-023-01869-1 | Verified |
| 71 | Barnett 2005 | 10.1093/ije/dyh299 | Verified |
| 72 | Genedy 2026 (erratum) | 10.1016/j.jacl.2026.03.024 | Verified (PMID 42203541) |
| 73 | NICE CG71 | none | Guideline, not DOI-indexed |
| 74 | Ferch 2025 | 10.1093/eurjpc/zwae331 | Verified |
| 75 | Chan 2018 | 10.1210/jc.2017-02622 | Verified |
| 76 | Mohammadnia 2023 | 10.1093/ehjdh/ztac059 | Verified |
| 77 | Stevens 2024 | 10.1161/JAHA.123.034434 | Verified |
| 78 | Gallo 2021 | 10.1016/j.jcmg.2021.06.011 | Verified |
| 79 | Paquette 2026 | 10.1093/eurjpc/zwag203 | Verified |
| 80 | Piaggio 2012 | 10.1001/jama.2012.87802 | Verified |
| 81 | Maia 2025 | 10.1016/j.jacl.2025.05.006 | Verified |
| 82 | Genedy 2026 (Atheroscler Plus TUDOR abstract) | 10.1016/j.athplu.2026.100590 | Not resolved |
| 83 | Genedy (CALON r7 manuscript) | none | Manuscript |
| 84 | Heinze 2001 | 10.1111/j.0006-341X.2001.00114.x | Verified |
| 85 | Bogsrud 2025 | 10.1093/eurheartj/ehaf815 | Verified |
| 86 | Genedy 2025 (EHJ abstract) | 10.1093/eurheartj/ehaf784.3670 | Not resolved |
| 87 | DeLong 1988 | 10.2307/2531595 | Verified (PMID 3203132) |
| 88 | McKearnan 2018 | 10.1093/aje/kwx374 | Verified |
| 89 | Tandy-Connor 2018 | 10.1038/gim.2018.38 | Verified |
| 90 | Ford 2009 | 10.1186/1472-6963-9-157 | Verified |
| 91 | Tybjaerg-Hansen 2005 | 10.1161/01.ATV.0000149380.94984.f0 | Verified |
| 92 | Singh 2026 | 10.1016/j.jacl.2026.03.018 | Verified |
| 93 | Mancini 2026 | 10.1016/j.atherosclerosis.2025.120590 | Verified |
| 94 | Mansilla-Rodriguez 2025 | 10.1093/eurjpc/zwaf631 | Verified |
| 95 | Tamehri Zadeh 2025 | 10.1016/j.cjca.2025.07.042 | Verified |
| 96 | Qureshi 2021 | 10.1136/openhrt-2021-001752 | Verified |
| 97 | Perez de Isla 2024 | 10.1016/S2213-8587(24)00192-X | Verified |
| 98 | Klevmoen 2023 | 10.1007/s11883-023-01155-6 | Verified |
| 99 | Van Calster 2023 | 10.1186/s12916-023-02779-w | Verified |
| 100 | Andaur Navarro 2024 | 10.1016/j.jclinepi.2024.111364 | Verified |

Tools used for verification: PubMed (ID converter, citation matcher, search and metadata). Scite's monthly quota was exhausted on the day of preparation; Consensus was not needed for DOI resolution. No DOI in this register was taken from memory.

*Prepared 27 September 2026 as a preparation document for the candidate. Every thesis number is quoted from MD_SUBMISSION_v3_FINAL.docx at the section cited.*

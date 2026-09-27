# MD Thesis Viva: Plain Questions and Plain Answers

Dr Nader Genedy. Thesis: Atherogenic Particle-Burden Phenotyping in Familial Hypercholesterolaemia. Training document, 27 September 2026.

454 questions. Read the question, answer aloud in your own words, then read the answer and correct yourself. Answers are in the first person and follow the thesis section named in brackets. Where the honest answer is a concession, the concession comes first.

**Q1. Summarise the thesis in two minutes.**

It asks one question. When does the same LDL-C value mean the same thing in familial hypercholesterolaemia, and what must be added when it does not? I tested five conditions in two linked genotyped cohorts, the All-Wales FH registry and UK Biobank. Particles: in 1,461 carriers of rare LDLR or APOB variants with 75 events, apoB above the value expected for LDL-C carried a hazard ratio of 1.356 per standard deviation, and a Welsh equation applied unchanged kept the association. Route into care: at matched LDLR residues, probands had untreated LDL-C 1.09 millimoles per litre higher than cascade-detected relatives, narrowing to 0.39 within families. Treatment: a treatment-aware triage score, TUDOR, ranked carriers with an AUC of 0.760 against 0.652 for reconstructed criteria, but fell to 0.669 when frozen in UK Biobank; and in the registry, three in four carriers with a severe untreated LDL-C were graded milder on the recorded value. Cohort: a five-term routine-data risk model kept its ranking in both directions, concordance about 0.70, but its absolute risks failed in opposite directions. And outside FH, the same apoB cut-point carried different absolute risk in women and men. LDL-C stays indispensable. It should not travel alone. No model here is ready for clinical use. Sections: Summary; 1.1; 1.5; 8.1.

**Q2. The Statement of original work says you did the statistics. What did Dr Aubin do?**

Statistical supervision and methodological advice. I designed, ran and verified the analyses in R and Python. Where a chapter's values come from a manuscript rather than a thesis rerun, the thesis says so: Chapter 5 reports Paper 14 R4 values that were not independently re-run, and Chapter 7 reports manuscript r7. I take responsibility for every number and for the evidence grading. (Statement of original work; Sections 5.3, 7.1.)

**Q3. Your AI statement is three sentences. Is that enough in 2026?**

It is the minimum, and I would expand it if the committee wishes. What it establishes is the boundary: AI tools helped with literature searching, drafting, editing and consistency checking; they had no access to participant-level data and ran no analysis; I reviewed every passage. Appendix F adds that title screening used a language-model first pass with single-reviewer adjudication and was not duplicated, and lists that as a limit. The evidence rule in Section 3.1 is the real protection: model-generated prose is never numerical or bibliographic evidence, and every number must trace to a publication, manuscript, package or script. (AI statement; Sections 3.1; Appendix F.8.)

**Q4. Why did UK Biobank access cease, and what did you lose?**

Access to the Research Analysis Platform closed during the final phase of the programme. I lost the ability to re-run participant-level UK Biobank analyses. What I kept are frozen aggregate packages; the Chapter 4 primary estimates sit in THESIS-AN-01, frozen on 12 July 2026. Any UK Biobank value that cannot be traced to a located producing output is labelled 'retained aggregate; not re-verifiable while the platform is closed' and is used only for sensitivity or description. The reason for closure should be stated in the thesis in one sentence, and I will add it. (Data governance statement; Sections 3.3.4; Appendix A.1.) ---


## Chapter 1: Introduction


### Section 1.1 The clinical paradox

**Q5. Why is FH a stringent setting for this question?**

Because LDL-C should be at its most informative in a disorder of the LDL receptor. If particles, route into care, treatment and cohort still change what the number means here, those dependencies cannot be dismissed as peculiarities of mixed dyslipidaemia in the general population. (Section 1.1.)

**Q6. Give me the single example you open with.**

The guideline severity cut-point is an untreated LDL-C of 4.9 millimoles per litre. In the registry, among genotyped carriers with a documented untreated value at or above it, applying the same cut-point to the value in their current record would have missed three in four. The value in the record was a valid on-treatment measurement; what it represented had changed. I will remove the thesis sentence that asserts analytical accuracy, because that was not tested. (Sections 1.1, 7.7: 145 of 193, 75.1 per cent.)

**Q7. Isn't this just treated LDL-C is lower than untreated? Where is the science?**

The science is in measuring how often, by how much, and what else does the same thing. Treatment is one of five conditions. The route into care changes the phenotype seen by 1.09 millimoles per litre at the same residue. The cohort of development changes a model's absolute risks by a factor of two in opposite directions. Particle concentration changes the event rate at the same LDL-C. Each is quantified in genotyped cohorts with common definitions, which had not been done together. (Sections 1.1, 1.6.)


### Section 1.2 Scientific basis

**Q8. Why is a current LDL-C not a measure of exposure?**

Because it records one point on a lifelong trajectory. It does not contain the preceding concentrations, treatment interruptions or accumulated arterial exposure. Similar present values can follow early diagnosis with sustained therapy, decades of untreated elevation, or intermittent adherence. The cumulative-exposure evidence I cite is the EAS consensus and the CARDIA-based analysis of cumulative LDL-C in young adulthood. (Section 1.2; refs 1 and 9, both DOI-verified.)

**Q9. Why is apoB not an LDL-specific particle assay?**

Because every apoB-containing particle carries one apoB molecule: VLDL, IDL, LDL and lipoprotein(a) carry apoB-100, chylomicrons and their remnants apoB-48. Routine plasma apoB therefore counts particles across those classes. It approximates particle number; it is not an LDL count. (Sections 1.2, 2.3; refs 10 to 12.)

**Q10. What is complementarity, not superiority, and why did you settle on it?**

The premise that apoB refines LDL-C rather than replacing it. Population studies show appreciable variation in apoB at the same LDL-C, and risk in discordant groups aligns more closely with apoB. Other studies find similar associations for apoB and non-HDL cholesterol, or conclusions that depend on how correlated lipids are modelled. Given that split, complementarity is the claim the evidence supports. (Section 1.2; refs 13 to 18.)

**Q11. Why does the thesis prefer a residual to a ratio, in one sentence each?**

A ratio is mechanically coupled to its denominator and assumes proportionality across the whole range. A residual asks how far observed apoB lies from the value expected at that LDL-C, and it declares the equation it depends on, which can then be frozen and tested. (Section 1.2; refs 19, 20.)

**Q12. What does back-calculating an untreated LDL-C actually give you?**

A model-based estimate, not a recovered measurement. Its accuracy depends on drug, dose, combination therapy, adherence, individual response and timing. It can preserve cohort means or rank order while being wrong for an individual. It may support triage or model development. It does not recover a person's true untreated value. (Section 1.2; refs 21 to 23.)


### Section 1.3 Selection into observation

**Q13. Define the ascertainment-conditioned phenotype without jargon.**

The FH phenotype as it appears in people selected by a particular route into care. It is a property of the selected group, not of the variant. It describes how the data were generated. It does not say that referral changes anyone's biology. (Sections 1.2, 1.3; Glossary.)

**Q14. Who showed this before you?**

Tybjaerg-Hansen and colleagues in 2005. Carriers of the same LDLR mutations showed progressively larger cholesterol increments when identified in the general population, among patients with ischaemic heart disease, and in FH clinics. They concluded a mutation's phenotype should be estimated in unselected carriers. I should note that their gradient is reported in total cholesterol; my text says LDL-C increments and I will correct the wording. (Section 1.3; ref 91, DOI verified 10.1161/01.ATV.0000149380.94984.f0.)

**Q15. Why does selection matter for a diagnostic test's AUC?**

Because sensitivity, specificity and AUC are properties of a test in a stated population, not of the score alone. Ransohoff and Feinstein called it spectrum bias in 1978. A score that separates conspicuous clinic cases from unaffected controls faces a harder problem among treated relatives or genotype-first carriers. (Section 1.3; ref 27, DOI verified 10.1056/NEJM197810262991705.)

**Q16. Are UK Biobank and the Welsh registry interchangeable halves of one population?**

No. UK Biobank offers genotype-first carriers with linked outcomes in a volunteer population that is itself selected. The Welsh service offers specialist, cascade and pedigree information within clinical care. Disagreement between them can expose spectrum, measurement or transport limits, if the model, outcome, time origin and receiving cohort are declared. Agreement does not establish universal validity. (Section 1.3; refs 6 to 8, DOI verified.)


### Section 1.4 The conceptual framework

**Q17. Is Measure to Select to Identify to Predict a causal pathway?**

No. It is the order in which the thesis reads the evidence, and Figure 1.1 says so. It is an inferential dependency: measurement defines the phenotype; selection determines which part of its distribution reaches analysis; identification asks whether the disorder can be recognised after treatment; prognostic modelling follows once population, starting point and outcome are defined. Uncertainty introduced early can propagate downstream if it is hidden in a reconstructed value or a model score. (Section 1.4.)

**Q18. When was the framework written?**

After the programme, as a thesis-level synthesis. It was not prespecified as one biological theory. Section 1.4 says that in the first paragraph, and Section 8.1 repeats it when the falsifiers are assessed.

**Q19. Name the five operational checks.**

Estimand. Target population. Sampling frame. Measurement model. Evidence for transport. Every result is meant to carry all five. (Section 1.4.)

**Q20. Why does the fifth condition, sex, sit in Chapter 8 rather than in an empirical chapter?**

Because its evidence comes from a population cohort, not an FH cohort. It is used as a limit on reading LDL-C as a single marker, examined in a bounded section, and explicitly not presented as an FH-specific or causal result. (Sections 1.2, 1.4, 8.4.)


### Section 1.5 Question, aim, objectives and falsifiers

**Q21. Your falsifiers were written after the results. Is that not hypothesising after the results are known?**

The falsifiers were written across a completed programme, and the thesis says so in Sections 1.4, 1.5 and 8.1. They do not alter the prespecification or evidence status of the source studies. What they do is force my reading of my own results to be checkable, and they set the terms of the resolving studies in Section 8.8, which are prospective. Where a source had a real prospective register, I report that instead: Paper 14 registered 120 analyses in fifteen families with false-discovery control. (Sections 1.5, 5.2.)

**Q22. Objective 4 was restated after results for both models were known. Explain.**

Objective 4 originally named CALON-C. When CALON-5 replaced it, the objective was rewritten to describe two models, CALON-5 and a two-term model, with neither designated primary, and the thesis states that this happened after results for both were known. I disclose it rather than backdate it. (Section 1.5.)

**Q23. State each objective's estimand in one phrase.**

One: the association between a continuous LDL-conditioned apoB residual and first ASCVD. Two: the difference in untreated LDL-C between probands and cascade-detected relatives at matched residues. Three: discrimination between genotype-confirmed carriers and non-carriers by TUDOR against reconstructed criteria. Four: the cumulative incidence of first MACE with competing death, ranking and calibration separately. Five: absolute risk at the same value and cut-point by sex and stage. (Section 1.5, Table 1.2.)

**Q24. Which falsifier came closest to being met?**

Objective 2, within families. The within-pedigree difference was 0.39 millimoles per litre with an interval from 0.05 to 0.73, but under a dummy-variable degrees-of-freedom correction the interval spans zero, minus 0.06 to 0.85. Table 8.1 records it as 'not met, narrowly'. (Sections 5.3, 8.1.)

**Q25. What is the object evaluated in Objective 3, and why does that phrase appear?**

The recovered executable TUDOR score. It differs from the published description in having no ascertainment-route term. The objective names the object so that no result is attributed to a model that does not exist as an executable. (Sections 1.5, 6.2.)


### Section 1.6 Contribution in relation to prior work

**Q26. What is new, given every component principle is established?**

The measurement of the five conditions together in linked genotyped cohorts with common definitions, so that their direction and relative size can be compared and the repair for each named. Section 1.6 gives the nearest prior work for each and states the increment. The claim is integrative and quantitative, not one of conceptual priority, and it is conditional on a search that did not include Embase. (Section 1.6; Appendix F.)

**Q27. Nearest prior work for the particle condition?**

Johannesen 2024, excess apoB in the Copenhagen General Population Study; Pan 2026, a frozen excess-apoB equation in statin-treated coronary disease validated externally in UK Biobank; and my own discordance paper in Welsh FH. The increment is a continuous residual in rare-variant carriers, a Welsh equation applied unchanged in UK Biobank, and the retirement of my own ratio threshold. (Section 1.6; refs 38, 39, 51, 68, all DOI-verified.)

**Q28. Nearest prior work for the treatment condition?**

Besseling 2017, where statin use predicts an FH-causing variant; FAMCAT, where drug class and potency are inputs; the Welsh service's own scoring criteria, Haralambos 2015; and Akyea 2026, where fully scored DLCN beat FAMCAT in a tertiary clinic. Chapter 6 adds pretreatment reconstruction with a shared-data comparison and frozen application; Chapter 7 adds the direct measurement of under-grading. (Section 1.6; refs 42, 43, 44, 45, 92, 93, all DOI-verified.)

**Q29. Nearest prior work for the cohort condition?**

McKay 2022 transported SAFEHEART-RE into English primary care by imputation; Mansilla-Rodriguez 2025 validated it in Australia; Tamehri Zadeh 2025 and 2026 compared FH scores on common data. Each ran in one cohort, one direction. Chapter 7 runs in both directions between two differently ascertained UK cohorts, separates ranking from calibration, models competing death, and reports calculability. (Section 1.6; refs 28, 29, 94, 95, all DOI-verified.)


### Section 1.7 Scope

**Q30. What is deliberately out of scope?**

Homozygous FH, structural prediction and AlphaFold-derived severity scoring, population-wide apoB discordance as a contribution, and Lp(a) and polygenic burden as anything more than modifiers. Incidental and privately obtained genetic findings are considered briefly in Section 8.6. (Section 1.7.)

**Q31. Three different genetic definitions are used. Why not one?**

Because the programme assembled them at different times under different data extracts. Chapter 4 uses an operational rare LDLR/APOB variant frame; Chapter 7 uses ClinVar pathogenic or likely-pathogenic single-variant carriers; the CALON-C lineage uses an operational LDLR-carrier flag. None is relabelled genetically confirmed HeFH, and Box 3.1 keeps them apart. Uniform adjudication is listed as what would remove the limitation. (Sections 1.7, 3.2; Table 8.3, L09.)


### Section 1.8 Structure

**Q32. Why do the chapters inherit an open problem from one another?**

Chapter 4 finds that particle concentration carries information at the same LDL-C, but measurement cannot be separated from who was measured; so Chapter 5 asks how the route into care shapes what is seen. Route and treatment together obscure the classical phenotype; so Chapter 6 asks whether a treatment-aware score can still find carriers. Finding a carrier does not determine prognosis; so Chapter 7 asks whether a risk model's ranking and absolute risks survive a change of cohort. (Section 1.8.) ---


## Chapter 2: Literature review


### Section 2.1 Purpose and stance

**Q33. This is not a systematic review. Why not, and what did you do instead?**

It is a critical, claim-led review, and the run state is recorded as 'search incomplete'. PubMed, Crossref including EHJ congress abstracts, the Europe PMC preprint index and OpenAlex were searched from 23 September 2016 to 23 September 2026, with one adversarial query per claim of absence. PubMed and Crossref were searched again on 27 September: 27 queries, 339 unique records, 23 retained. Embase was not searched; Elicit, Consensus, Scite and SciSpace returned no records; screening stopped at the top 20 to 25 records per query and was single-reviewer after a labelled language-model first pass. Every 'no study identified' sentence is conditional on that. (Section 2.1; Appendix F.)

**Q34. What are the four load-bearing propositions?**

One: LDL-C and apoB are correlated but not interchangeable. Two: clinical criteria were built around a visible, commonly untreated phenotype, whereas detection now often occurs after treatment, through cascade or genotype-first routes. Three: treatment and ascertainment determine which version of the phenotype reaches the analyst. Four: prediction models require a declared estimand and evaluation across discrimination, calibration, overall error and utility. (Section 2.1, Table 2.1.)


### Section 2.2 FH as exposure and recognition

**Q35. What are the three category errors you want to prevent?**

A low on-treatment LDL-C does not negate previous exposure. Failing a clinical score does not prove absence of a disease-causing variant. A pathogenic variant does not specify present absolute risk. (Section 2.2.)

**Q36. Why call FH both an exposure and a recognition event?**

Because the biological exposure begins before recognition and may never be recognised. Consensus statements place cumulative LDL exposure, not any single result, at the centre of risk. The recognition event, referral, testing, registration, is what makes the person visible to a dataset. (Section 2.2; refs 1 to 3, 57, DOI-verified.)


### Section 2.3 LDL-C, apoB and atherogenic burden

**Q37. Give me the strongest evidence that apoB carries information beyond LDL-C.**

In NHANES, at LDL-C of 100 milligrams per decilitre the middle 95 per cent of apoB ran from 66 to 99. In statin-treated cohorts, apoB retained 1.27 per standard deviation for myocardial infarction when the other measures did not after mutual adjustment. In Copenhagen, discordantly high apoB with low LDL-C carried a hazard ratio of 1.49 for infarction. (Section 2.3; refs 15, 13, 14, DOIs verified: 10.1001/jamacardio.2024.1310; 10.1001/jamacardio.2021.5083; 10.1016/j.jacc.2021.01.027.)

**Q38. And the strongest evidence against apoB superiority?**

The Emerging Risk Factors Collaboration, 302,430 participants: broadly similar association shapes for apoB and non-HDL cholesterol. Sniderman's own INTERHEART re-analysis: adjustment among strongly correlated lipids can produce unstable claims of independence. And the Mendelian randomisation studies disagree by specification, Richardson favouring apoB, Helgadottir favouring cholesterol content. (Section 2.3; refs 16, 17, 11, 18. DOIs verified for 16, 11 and 18; ref 17 is listed in the thesis with a PMID but did not resolve in the tools used here.)

**Q39. Does NMR settle whether it is particle number or cholesterol per particle?**

It moves the argument toward particle number without closing it. In 89,422 statin-free UK Biobank participants with 3,821 coronary events, cholesterol molecules per LDL particle were not associated with disease after adjustment for particle concentrations, hazard ratio 1.03, and apoB correlated 0.99 with LDL-particle concentration. The adjustment there was for VLDL, LDL and HDL particle concentrations together; my text says LDL-particle concentration and I will make it exact. (Section 2.3; ref 61, DOI verified 10.1161/JAHA.123.029552.)

**Q40. Where does lipoprotein(a) sit?**

As a particle-composition and prognostic modifier. Each Lp(a) particle carries one apoB-100 and contributes variably to measured LDL-C, but its risk is not reducible to LDL-particle concentration. Extreme Lp(a) has been described as a risk equivalent in HeFH. It is in SAFEHEART-RE and the FH-Risk-Score. (Section 2.3; refs 53, 56, 64, DOI-verified.)


### Section 2.4 Why discordance in FH remains uncertain

**Q41. Why is equipoise appropriate rather than an expectation of a positive result?**

Two priors oppose each other. Receptor impairment could tighten the mass-to-particle relation so that LDL-C is an unusually faithful proxy and the residual adds nothing. Or metabolic heterogeneity coexists with the variant and leaves clinically relevant variation around the relation. Selection and measurement can obscure either. So the chapter was designed to be able to return null. (Section 2.4.)

**Q42. What is the nearest construction to yours, and what is the gap?**

Johannesen 2024: expected apoB from regressing apoB on LDL-C among people with triglycerides at or below 1 millimole per litre, in people not taking statins. Competitors are Kim's residual discordance against coronary calcium progression and Rehman's risk-weighted apoB. The gap is five-fold: population, FH not general; equation, derived in one FH setting and applied frozen in another; endpoint, incident ASCVD; treatment state, mixed; and validation design, an independent evaluation, which Chapter 4 does not supply. (Section 2.4; refs 38, 67, 68, DOIs verified.)

**Q43. Pan 2026 already validated a frozen excess-apoB equation in UK Biobank. What is left?**

Pan's cohort was statin-treated coronary disease and the outcome was mortality; excess apoB carried adjusted hazard ratios of 1.12 for all-cause and 1.24 for cardiovascular death. That establishes the method outside FH. What was left was FH, incident ASCVD, and an equation derived in an FH registry. (Section 2.4; ref 39, DOI verified 10.1186/s12944-026-02928-z.)


### Section 2.5 Ratio versus residual

**Q44. Why can adjusting a ratio for its denominator not fix the coupling?**

Because if the ratio is X over Y, and Y enters the model too, you are comparing a transformed variable with one of its own components. Archie described the error in 1981; Lolli showed it recently in a workload ratio. A quotient also assumes proportionality across the range, and in FH the range spans treated and untreated states. (Section 2.5; refs 19, 20, DOIs verified.)

**Q45. Write the residual on the board.**

Residual apoB for person i equals observed apoB minus expected apoB given that person's LDL-C, under a prespecified equation. Positive means more particles than the cholesterol mass implies. It is not a count per unit cholesterol, not a category, not a particle-size class. (Section 2.5.)

**Q46. Is a threshold ever legitimate?**

Yes, when its assumptions are appropriate, its components reliably measured, and its portability shown rather than assumed. An outcome-selected threshold in a small cohort captures sample-specific variation and invites regression to the mean; Barnett's paper and its correction are the reference. (Section 2.5; ref 71, DOI verified 10.1093/ije/dyh299.)

**Q47. What did the erratum to your discordance paper change?**

Typographical errors in the published highlights. The third highlight had been truncated; the four bullets were replaced. No estimate, table, method or conclusion changed. The thesis's retirement of the ratio estimate is the thesis's re-analysis, not the erratum. (Sections 2.5, 4.3; refs 51, 72, DOIs verified 10.1016/j.jacl.2025.11.008 and 10.1016/j.jacl.2026.03.024.)


### Section 2.6 Treatment and reconstruction

**Q48. Why is medication status a care-process variable?**

Because in observational data it encodes an intervention and the severity that prompted it, plus access and adherence. Its coefficient reflects confounding by indication as much as pharmacology. (Section 2.6.)

**Q49. What does the Simon Broome Register tell us about reconstruction?**

That a stable mean relation between lipid measures can coexist with wide individual limits of agreement. Soran and colleagues showed it for non-HDL versus LDL-C. So a reconstruction can be approximately unbiased at cohort level and unreliable for one person. Correlation does not demonstrate agreement. (Section 2.6; ref 23, DOI verified 10.1097/MOL.0000000000000692.)


### Section 2.7 Ascertainment and spectrum effects

**Q50. What did Mourre 2025 show and why is it not your analysis?**

Cascade-screened individuals started statins about 14 years earlier and had fewer cardiovascular events than opportunistically screened ones, and crude route differences attenuated after accounting for case mix. It shows route conditions the treatment and outcome spectrum. It is not residue-matched or within-pedigree, and it cannot show that route changes biology. (Section 2.7; ref 40, DOI verified 10.1093/eurjpc/zwaf234.)

**Q51. What is the polygenic alternative explanation?**

Trinder 2024: clinically diagnosed FH in a Canadian registry had higher LDL-C and more ASCVD than genetically identified FH in UK Biobank, with a different polygenic background. Residue matching cannot exclude that. Polygenic scores were not available in the Welsh frame. (Section 2.7; ref 41, DOI verified 10.1161/ATVBAHA.123.320287.)

**Q52. UK Biobank is 5 per cent of invitees. Does reweighting fix it?**

Reweighting on measured determinants can reduce bias in exposure-outcome associations, as van Alten and Schoeler showed, but it cannot recover unmeasured determinants or people outside the sampling frame. Large sample size reduces random error. It does not remove selection. (Section 2.7; refs 6 to 8, DOIs verified.)


### Section 2.8 Clinical criteria and case-finding

**Q53. Why is missing not absent for a tendon examination?**

Because an unrecorded examination is not a documented negative one, and an unknown pedigree is not a negative family history. Assigning zero to unavailable information lowers sensitivity; optimistic imputation no longer reproduces the published criterion. (Section 2.8; refs 74, 75, DOIs verified.)

**Q54. What did Mohammadnia show?**

In genetically confirmed patients, algorithm sensitivity differed when assessment relied on coded records rather than all available information. Criterion performance and extraction performance are different evidence objects. (Section 2.8; ref 76, DOI verified 10.1093/ehjdh/ztac059.)

**Q55. Why is FAMCAT's PPV low, and is that FAMCAT's fault?**

Because FH is uncommon in primary care; external validation gave a positive predictive value under 1 per cent. That is prevalence, not the algorithm. It makes confirmatory workload integral to implementation. (Section 2.8; ref 43, DOI verified 10.1016/S2468-2667(19)30061-1.)


### Section 2.9 From identification to prognosis

**Q56. Why can't you rank SAFEHEART-RE, Montreal and the FH-Risk-Score from their published statistics?**

Because they target different populations and outcomes. SAFEHEART-RE: incident events in a registry including people with prior disease. Montreal: principally prevalent disease. FH-Risk-Score: incident events in primary prevention. A prevalent-disease AUC, an incident concordance and a primary-prevention estimate are not interchangeable. (Section 2.9; refs 53 to 56, DOIs verified.)

**Q57. What does the polygenic-score example teach?**

In an FH cohort a polygenic score had a hazard ratio of 1.77 yet moved the C-statistic only from 0.746 to 0.750, p 0.60. Association without incremental prediction. (Section 2.9; ref 79, DOI verified 10.1093/eurjpc/zwag203.)


### Section 2.10 Evaluating prediction

**Q58. Name the four questions of model evaluation and one measure for each.**

Ranking: concordance. Agreement: calibration intercept, slope and curve. Overall error: Brier. Decisions: net benefit. None substitutes for another. (Section 2.10, Table 2.3; refs 30 to 33, DOIs verified.)

**Q59. Why must competing death be handled?**

Because non-cardiovascular death prevents a later first ASCVD event. Treating it as censoring assumes the person remained observable and event-free. One minus Kaplan-Meier overestimates absolute risk; Aalen-Johansen does not. Cause-specific and subdistribution hazards answer different questions. (Section 2.10; refs 34, 35, DOIs verified.)

**Q60. Why is checklist adherence not validation?**

Because STROBE, RECORD, TRIPOD+AI and PROBAST+AI expose omissions and risks. They cannot correct biased sampling, temporal leakage, incomplete outcomes, overfitting, absent calibration or inappropriate transport. (Section 2.10; refs 30, 31, 36, 37, 49, DOIs verified.)


### Section 2.11 Life course, sex and modifiers

**Q61. What is the FH-specific evidence on menopause?**

Da Roza 2026: within-person lipid change across the menopausal transition in women with FH, larger in monogenic FH. Johansen: over 12 years young women with FH accumulated a higher LDL-C burden than men. Klevmoen: loss of statin-treatment years in pregnancy and breastfeeding. None establishes a uniform effect across genotypes or treatment. (Section 2.11; refs 24, 25, 26, DOIs verified.)


### Section 2.12 Synthesis

**Q62. What single rule does the review converge on?**

Each claim must remain the size of the design that supports it. (Section 2.12.) ---


## Chapter 3: Shared methods and governance


### Section 3.1 Stance

**Q63. What are the six components of your inferential contract?**

Distinct cohort roles. A declared eligibility, unit, time zero, endpoint and estimand for every analysis. Observed and reconstructed values kept as different objects. Family dependence, missingness, competing events and transport treated as science, not afterthought. Every number mapped to a frozen cohort, producing analysis and tier. And a bidirectional lock between methods and results. (Section 3.1.)

**Q64. What are the THESIS-AN packages?**

Thirteen thesis-reanalysis identifiers, THESIS-AN-01 to 13, each recording its frame, producing artefact or its unavailability, and freeze date. A reanalysis is cited to its identifier, never to a publication that does not contain it. (Section 3.1; Appendix A.)


### Section 3.2 Estimands

**Q65. Why can one statistical model not answer the thesis question?**

Because the four tasks have different populations, outcomes, temporal structures and estimands: a cause-specific hazard, a mean difference, a binary AUC, and a fixed-horizon cumulative incidence with competing death. (Section 3.2, Table 3.1.)

**Q66. What is the causal boundary of the thesis?**

Chapter 5 does not estimate an intervention on ascertainment. Chapters 4 and 7 do not estimate effects of altering residuals or scores. Chapter 6 does not estimate outcomes caused by testing. The thesis does not acquire those estimands through stronger wording. (Section 3.2.)

**Q67. Two digits, 0.660 and 0.725, appear for different quantities. Explain.**

In the CALON-C lineage, 0.7252 is the frozen nine-term concordance in the Welsh primary frame, and 0.6600 is the reverse frozen application of a seven-term Welsh refit to UK Biobank; a different 0.660 is the endpoint-date stress test in UK Biobank. That is why every value is written with its model, frame and horizon. (Sections 3.2, 7.14; Appendix B.)


### Section 3.3 Data sources

**Q68. What is the frozen registry, and is 4,570 an analytical denominator?**

4,570 deduplicated individuals including 887 LDLR carriers, the TUDOR evaluation resource. It describes the data resource, not any chapter's denominator. (Section 3.3.1.)

**Q69. What is the UK Biobank extract?**

501,936 participants; 3,540 with the operational LDLR-carrier flag; a lipid-clinic-eligible subset of 49,427 with 921 flagged carriers. The flag is broader than ClinVar P/LP and is not labelled HeFH. (Section 3.3.2.)

**Q70. Why is the flag frequency, one in 142, a problem?**

Because it exceeds published estimates of pathogenic FH-variant prevalence, and the extract lacked identifiers for complete re-adjudication. That is why Chapter 7 restricts to ClinVar P/LP, and why the CALON-C development population is described as an operational-flag cohort. (Appendix B.2; refs 4, 5, DOIs verified.)

**Q71. What is disclosure control?**

Aggregate results only: no participant or family identifiers, no exact dates, no row-level predictions, no rare combinations. Statistical verification and governance approval are separate release gates. (Section 3.3.4.)


### Section 3.4 Data engineering

**Q72. What is the row-count ledger?**

Pre- and post-join row counts and unique-participant counts, with unmatched keys and duplicates recorded. Unexplained multiplicative expansion blocks the result; it is not repaired by silent deduplication. (Section 3.4.1.)

**Q73. Why are outliers not excluded at 1.5 IQR?**

Because extreme FH measurements can be genuine. Values are checked against units, assay and coding provenance; influence is assessed by diagnostics; named sensitivities examine influential points without altering the source record. (Section 3.4.2.)

**Q74. How was LDL-C derived in Wales?**

Friedewald below the documented triglyceride threshold, 4.5 millimoles per litre in Chapter 5, and direct assay above it. Registry-wide assay homogeneity is not assumed; measurements span calendar periods and platforms. (Sections 3.4.3, 5.2.)

**Q75. What is your missing-data policy?**

No programme-wide mandatory method. Complete-case analysis is acceptable only when its assumptions and denominator are stated. Multiple imputation, where used, must respect outcome, interactions, survival structure and clustering, and happen inside resampling splits. Missing, unrecorded, structurally unavailable and examined-but-absent stay distinct. (Section 3.4.4.)


### Section 3.5 Outcomes, time zero and competing events

**Q76. Why not pool the event sets across chapters to gain power?**

Because they are incompatible endpoint lineages. Chapter 4 retains four confirmation-dependent events; the CALON-C lineage has uncertain component dates for 142 of 289 events. A larger event count does not supersede a narrower source-locked endpoint. (Section 3.5.)

**Q77. What is time zero in each chapter?**

Chapter 4: the baseline measurement defining the apoB-LDL-C phenotype. Chapter 5: none, it is a phenotype contrast. Chapter 6: assessment, a classification. Chapter 7: baseline free of atherosclerotic disease, MACE at five and ten years with competing death. (Section 3.5.)

**Q78. How many competing deaths, and why did the count change?**

224 in UK Biobank and 35 in each Welsh frame after a transparent death-first reconstruction; legacy lineages had reported 193 and one. Reconciling the counts did not complete CALON-C absolute-risk validation, because the historical fitted object and imputation pathway are incomplete. (Section 3.5; Appendix B.4.)


### Section 3.6 Treatment reconstruction

**Q79. State the three conventions and what each can support.**

A: intensity-specific residual fractions, for discordance and TUDOR; supports agreement and error analysis. B: population-average scalars, 0.70 and 0.80, for CALON-C; supports rank stability across a factor sweep. C: a regression on 684 Welsh pairs, untreated non-HDL-C equals 5.2299 plus 0.3345 times treated, for CALON-5; supports group-level targeting of the untreated quantity. None supports an individual untreated value. (Section 3.6, Table 3.4.)

**Q80. Why do you have three?**

Because three studies were built at different times for different modelling needs, and the record does not establish one rationale unifying them. Appendix D says so. They are kept separate and never used to validate one another. (Section 3.6; Appendix D.)


### Section 3.7 Statistical principles

**Q81. What is your position on p-values?**

Effect estimates and 95 per cent intervals lead. P-values are compatibility summaries. An interval including zero is 'no statistically detectable difference in that sample', never equivalence, which needs a prespecified margin. (Section 3.7.1; ref 80, DOI verified 10.1001/jama.2012.87802.)

**Q82. Which diagnostics were completed and which were not?**

Reference R diagnostics exist for reconstructed CALON-C states. For CALON-5 the record documents the proportional-hazards check and whole-procedure bootstrap optimism, but not influence or functional-form diagnostics for the log-ratio term. The exact Chapter 4 proportional-hazards object was not recovered. Producing a coefficient never establishes that assumptions passed. (Section 3.7.1.)

**Q83. Family dependence, chapter by chapter.**

Chapter 5: pedigree-clustered errors, wild cluster bootstrap, pedigree fixed effects. Chapter 7: registry resampling grouped by family, both models on the same resamples. CALON-C high-bar rerun: FamilyNumber for folds and bootstrap. TUDOR: DeLong treats participants as independent, and no family-aware optimism estimate exists. UK Biobank: singleton clusters because the extract had no kinship field; independence is not assumed. (Section 3.7.1.)

**Q84. What is your multiplicity position?**

No programme-wide correction. Estimates are read by interval and by prespecification. Where a source applied a procedure it is reported with its family: BH-FDR across the Paper 14 register; Holm in the CALON-C lineage. (Section 3.8.)


### Section 3.8 Prediction development and transport

**Q85. Name the five model states and the only one that is external validation.**

Frozen application. Target-adapted preprocessing. Recalibration. Updating. Target refitting. Only the first evaluates the source model unchanged. (Section 3.8; refs 28, 29, 32, DOIs verified.)

**Q86. What is the sequential interval rule?**

Baseline and slope re-estimated where the slope interval excludes one; baseline alone where only the observed-to-expected interval excludes one; otherwise no update. It is not a closed testing procedure, and the thesis says so. (Section 3.8.)

**Q87. What is the twenty-event floor?**

Comparisons below 20 events at a declared horizon are not interpreted. Precision for concordance is governed by events and, for paired comparisons, by the correlation between scores, not by cohort size. (Section 3.8.)

**Q88. What does the Brier score not tell you?**

It is not a calibration measure, and when events are rare a small Brier can coexist with poor ranking. It is stated once in Section 3.8 and referred to elsewhere.


### Section 3.9 Reporting and evidence maturity

**Q89. Recite the evidence tiers.**

Published or corrected. Locked internal. Externally evaluated. Exploratory. Developmental. Null or inconclusive. Gated or excluded. And the tier codes: T1 version of record, T2 accepted or in press, T3 submitted, T4 frozen thesis analysis, T5 developmental, T6 not admissible. (Section 3.9, Table 3.6; Appendix A.1.)

**Q90. What is currently gated?**

CALON-C absolute risk, calibration, overall error and utility. No interval for CALON-5's increment over age and sex. TUDOR's prospective calibration, utility and impact. The reconstruction-uncertainty propagation, frozen as a protocol but not run. (Section 3.9.)

**Q91. What is the no-do-harm adoption gate?**

Before any construct moves toward use: no material loss in calibration, stability, subgroup performance or utility on transport; no reduction in testing or treatment where a known variant, pedigree or physical sign already indicates it; no substitution for sequencing or guideline-directed lipid lowering; then prospective assessment. Discrimination alone and retrospective decision curves do not pass it. (Section 3.9.)


### Section 3.10 Reproducibility

**Q92. What does zero silent correction mean?**

Every change to a denominator, term, value or source status, and every deletion or move, is recorded in the change ledger. The audit package does not retrospectively supply missing participant-level provenance; alternative builds stay quarantined until adjudicated. (Section 3.10.)

**Q93. Which results are independently reproduced and which are source-locked?**

The final CALON-C high-bar rerun met the manifest-level target, and the CALON-5 lock records its object under one hash. Several inherited Chapter 4 to 6 results did not retain the whole chain and remain source-locked rather than independently reproduced. (Section 3.10.) ---


## Chapter 4: Measure. Do particles carry information at the same LDL-C?


### Section 4.1 The measurement question

**Q94. State the chapter's question and its non-claims.**

At the same LDL-C, does apoB-containing particle concentration carry additional information about first ASCVD events in rare LDLR or APOB variant carriers? It does not estimate the causal effect of lowering apoB at fixed LDL-C, it does not compare apoB-guided with LDL-C-guided treatment, and it does not validate a bedside calculator. (Sections 4.1, 4.2.)

**Q95. Why a staged correction rather than simply presenting the residual?**

Because the published paper used a ratio and a threshold, and the residual is post hoc relative to it. The chapter preserves the published analysis, diagnoses its sparse-data, temporal and measurement limits, then re-specifies the exposure. The three stages have dates: the paper was submitted in August 2025 and published in March 2026; the residual reanalysis was frozen on 12 July 2026 as THESIS-AN-01. (Section 4.1.)


### Section 4.2 Cohorts, measurements and estimands

**Q96. Describe the two frames precisely.**

The published Welsh study: 424 genotype-confirmed adults, median follow-up 9.1 years, 61 with ASCVD of whom 41 were prevalent and 20 incident. The thesis primary: 1,461 UK Biobank participants selected by an operational rare LDLR/APOB variant rule, 75 incident events, none of which appears in the published paper. (Section 4.2; ref 51, DOI verified 10.1016/j.jacl.2025.11.008.)

**Q97. Were participants with prevalent ASCVD excluded from the 1,461?**

Follow-up began after the baseline measurement and the outcome is first incident ASCVD. The exclusion rule and the number excluded are not stated explicitly in Chapter 4 as they are in Chapter 7, and I will add that sentence. (Sections 3.5, 4.2; see Appendix 1, item P7.)

**Q98. Four of 75 events are unconfirmed. Why keep them?**

Seventy-one of the 75 are supported by a core ASCVD code, none arises solely from aortic-stenosis coding, and four await procedure or cause-of-death confirmation. They stay in the source-locked primary because the influence and deletion sensitivities show the association does not depend on the most influential events. Their status remains a load-bearing outcome limitation and is stated as such. (Section 4.2.)

**Q99. What is the difference between the local and the transported residual?**

The local residual is fitted within UK Biobank and standardised by the UK Biobank SD. The transported residual uses the apoB-on-LDL-C equation fitted in Wales and the Welsh residual SD, applied unchanged. They answer different questions and are not interchangeable. (Section 4.2.)

**Q100. Three estimands were kept separate. Name them.**

Association: the cause-specific rate of first incident ASCVD, by Cox. Absolute incidence: ten-year cumulative incidence with competing death, by Aalen-Johansen. Incremental performance: whether the residual changed concordance beyond an existing predictor set. Evidence for one does not establish the others. (Section 4.2.)


### Section 4.3 The published ratio analysis

**Q101. What did the published paper actually find?**

In 424 carriers, an apoB to LDL-C ratio of 0.31 grams per millimole or above gave an AUC of 0.726 for combined prevalent and incident ASCVD, sensitivity 41.2 per cent, specificity 86.7 per cent. Age alone gave 0.845. In 14 treatment-naive pairs the ratio ran from 0.18 to 0.38, a 2.1-fold spread. The ridge-penalised threshold hazard ratio was 38.55 with an interval from 3.72 to 399.36, and the paper called itself hypothesis-generating. (Section 4.3.)

**Q102. Your version of record contains three different proportions for the same comparison. Explain.**

Yes. The prose reports 27.6 versus 11.8 per cent; Table 3 reports 25 of 121, 20.7 per cent, versus 36 of 303, 11.9 per cent; and a third rendering uses a 0.3 boundary with 46.9 versus 22.4 per cent. The thesis preserves all three as version-of-record inconsistencies rather than choosing one, and none supplies the thesis's active estimate. The erratum does not address them. (Section 4.3.)

**Q103. What were the paper's limitations, in your own words?**

ApoB missing in 40 per cent, handled by 20 imputations. Only 20 incident events. A dichotomised, noisy ratio. And a headline comparison that mixed disease before lipid measurement with disease during follow-up. (Section 4.3.)

**Q104. You have two congress abstracts on the same cohort. Are they replication?**

No. One study, two abstract publications with near-identical text, 425 patients against the paper's 424, using an LDL-conditioned definition rather than the ratio. They are earlier, differently specified outputs from a near-identical cohort, disclosed as such, and no estimate from them is used. (Section 4.3; ref 66, not indexed in PubMed.)


### Section 4.4 Why the historical hazard ratio was unstable

**Q105. Walk me through the sparse-data geometry.**

Sixty-one disease observations, 20 incident. Cross-classify a binary exposure with several covariates and many cells are near-empty. That permits near-separation and highly mobile coefficients. Penalisation prevents numerical divergence; it cannot manufacture information absent from the risk sets. (Section 4.4; ref 84, DOI verified 10.1111/j.0006-341X.2001.00114.x.)

**Q106. Why does temporal ordering matter here?**

Because 41 of 61 disease observations were prevalent. Prior disease influences later treatment, lipids and surveillance, so comparing a contemporary ratio with earlier disease risks reverse ordering and confounding by indication. Restricting to incident events restores temporality and leaves 20 events, exactly where the model became unstable. (Section 4.4.)

**Q107. Is 0.31 a biological constant?**

No. It is a study-specific cut-point. Published discordance studies define excess apoB by a residual from a reference equation, not a fixed ratio, and the bounded search found no independent validation of this orientation or boundary. (Section 4.4; refs 38, 67, DOIs verified.)

**Q108. What happened to the Firth estimate and the negative controls?**

The Firth estimate mentioned in programme documentation could not be traced to a verified source value or producing object and is withheld. The negative-control analyses were non-estimable from the aggregate inputs. The only active adjusted threshold re-analysis is 1.24, p 0.44, in the UK Biobank frame. (Section 4.4.)


### Section 4.5 Residual re-specification and robustness

**Q109. Give me the primary estimate and its adjustment ladder.**

Hazard ratio 1.356 per standard deviation, 95 per cent interval 1.070 to 1.718, in 1,461 participants with 75 events. With age and sex: 1.556. Adding log-triglycerides: 1.469. Adding HDL-C and diabetes: 1.356, the primary. Age as the timescale with delayed entry: 1.467. The movement is consistent with shared information between residual apoB and metabolic risk factors; it does not identify a causal specification. (Section 4.5, Table 4.2.)

**Q110. Why is the competing-risk estimate non-load-bearing?**

Because it was an inverse-probability-of-censoring-weighted Fine-Gray-style approximation, not an exact fit, and its implementation was not validated under realistic censoring. It gave 1.334, directionally similar, and is reported as unverified. The estimands also differ: rate among the event-free versus the cumulative-incidence process. (Section 4.5; ref 34, DOI verified 10.1161/CIRCULATIONAHA.115.017719.)

**Q111. What is the absolute-incidence result?**

Aalen-Johansen ten-year cumulative incidence 7.38 per cent in local residual quartile 4, 366 participants, against 2.37 per cent in quartiles 1 to 3, 1,095 participants. A locally defined descriptive grouping. Not a threshold, not calibrated risk. (Section 4.5.)

**Q112. Does the residual improve discrimination?**

Barely, at this precision. Overall concordance 0.697. Adding the residual to age and sex: plus 0.050, interval 0.003 to 0.108, but that base had no LDL-C. Beyond LDL-C and conventional predictors: about plus 0.011, no frozen interval. A marker can keep a stable hazard coefficient and add little ranking because age already orders many people, the effect is moderate, and 75 events limit precision. (Section 4.5.)

**Q113. You tested apoB information in a prediction model and it failed. Why is that in this chapter?**

Because it is the programme's only direct test of whether the measurement signal adds ranking to a model. In the CALON-F grey zone of 5 to 20 per cent predicted risk, 1,685 of 3,333 scored with 218 events, adding log apoB over LDL-C changed concordance by plus 0.0146, interval minus 0.0050 to 0.0330; Lp(a) minus 0.0038; both plus 0.0118. It used the ratio, an unfrozen model and a different frame, so it does not test the residual. The link from measure to predict is carried forward as a hypothesis. (Section 4.5; ref 22, developmental manuscript.)

**Q114. What is the whole-population result for?**

Context only. 406,902 participants, 17,680 events, hazard ratio 1.20 per SD, interval 1.18 to 1.22, retained aggregate. It is consistent with the population literature and does not enlarge the FH claim. (Section 4.5; refs 13 to 15, DOIs verified.)


### Section 4.6 Retirement of the ratio threshold

**Q115. On what criteria did you retire the published hazard ratio?**

Sparse incident events. An interval spanning two orders of magnitude. Mathematical coupling. Dichotomisation. No independent validation. And failure of the adjusted threshold association to transport: 1.24, interval 0.72 to 2.16. Retirement is not repudiation: the paper identified the measurement problem and made its uncertainty visible. (Section 4.6.)

**Q116. Could a clinician use 0.31 tomorrow?**

No. Action would need reproducible classification across assays, an absolute-risk model with calibration, a justified action threshold, decision-curve evidence and prospective assessment. None exists. It must not be described as a treatment trigger. (Section 4.6.)


### Section 4.7 Biological interpretation

**Q117. What is the residual, biologically?**

ApoB unexplained by LDL-C under the fitted relation. Positive values may arise from more cholesterol-depleted LDL particles, triglyceride-rich remnants, lipoprotein(a), treatment-related compositional change, assay effects, or combinations. The coefficient cannot allocate among them. (Section 4.7.)

**Q118. Did remnant cholesterol explain it?**

Calculated remnant cholesterol accounted for an estimated 16.0 per cent of the association, interval minus 0.9 to 44.6. NMR remnant cholesterol 6.68 per cent, interval minus 4.34 to 39.54, in 1,450 participants with 75 events. Both cross zero; compatible with nothing and with a relevant contribution. (Section 4.7.)

**Q119. In the joint model LDL-C was protective. Is it?**

No. Conditional on both, apoB was 4.94 and LDL-C 0.26. With predictors this correlated each coefficient is the residual variation after conditioning on the other, and reversal reflects unstable partitioning. That model is shown to explain why the residual is the interpretable representation. (Section 4.7.)

**Q120. What did the polygenic score add?**

Association at 1.43 per SD, interval 1.12 to 1.84, in 1,442 participants with 61 events, with concordance 0.569 and 64.1 per cent missing. It cannot explain the residual or support an incremental claim. (Section 4.7.)

**Q121. Does the Mendelian randomisation literature convert your association into a causal target?**

No. Richardson's multivariable analysis gave apoB an odds ratio of 1.92 while LDL-C attenuated; Helgadottir emphasised cholesterol content. Neither converts an observational residual into a treatment contrast, and no intervention selectively altered residual apoB. (Section 4.7; refs 11, 18, DOIs verified 10.1371/journal.pmed.1003062 and 10.1093/eurjpc/zwac219.)


### Section 4.8 Frozen-equation application

**Q122. Why freeze the Welsh equation at all?**

Because a locally fitted residual may exploit the receiving cohort's own apoB-LDL-C distribution. Freezing the derivation relation and the Welsh scaling is the stronger test. The frozen residual gave 2.22, interval 1.40 to 3.53, in UK Biobank. (Section 4.8.)

**Q123. Is the effect twice as large in UK Biobank?**

That inference is not permitted. A Welsh-scaled SD is a different absolute apoB departure; assays, treatment and event structure differ. What transported is direction and an interval excluding unity. Not magnitude. (Section 4.8.)

**Q124. Is this external validation?**

It is external evaluation of a measurement construct across an ascertainment boundary. Not clinical validation of a risk model: no baseline hazard, probability, threshold or decision curve was transported. (Section 4.8.)


### Section 4.9 The never-treated boundary

**Q125. Why not analyse only untreated carriers?**

I did, as a boundary: 1.248 per SD, interval 0.885 to 1.760, on 30 events, retained aggregate. It neither confirms nor refutes. Never-treated participants are selected, younger and less severe, so removing treatment adds another selection mechanism. The primary claim stays in the mixed frame. (Section 4.9.)


### Section 4.10 Contribution and limits

**Q126. State the answer to the particles part in three sentences.**

At the same LDL-C, particle concentration carried additional information about first ASCVD events. A residual defined in one setting kept its association when applied unchanged in another; the ratio threshold did not. Whether that information improves risk ranking in a model was not shown. (Section 4.10.)

**Q127. What did you not compare the residual with?**

Risk-weighted apoB, Rehman's weighted sum of apoB, triglycerides and lipoprotein(a), which improved Harrell's C over apoB alone in general populations. So I cannot say which representation carries more prognostic information; that is limitation L10 and a resolving-study item. (Section 4.10; ref 68, DOI verified 10.1093/eurheartj/ehaf1124.)

**Q128. Perimeter of the claim, in one breath.**

Observational; 75 events, four unconfirmed; the Welsh source had 40 per cent apoB missingness and 20 events; residualisation depends on its derivation relation and assay; treatment history incompletely characterised; never-treated analysis underpowered; one receiving cohort; family and referral structure in Wales, volunteer selection in UK Biobank. (Section 4.10.) ---


## Chapter 5: Select. Does the route into care change the visible phenotype?


### Section 5.1 Selection as part of the phenotype

**Q129. What is the estimand, and why is it not causal?**

The observed difference in assigned untreated LDL-C between carriers entering the same service by different routes, after controlling approximately for the affected LDLR residue and then reducing age, pedigree and provenance differences. Ascertainment is not a biological intervention, and proband status is partly defined by the phenotype under study. The group contrast describes selected distributions, not an individual counterfactual. (Section 5.1.)

**Q130. Explain Figure 5.1.**

Underlying LDL-C, variant severity and family history influence the route by which a carrier is identified. Route influences age at encounter and treatment history, which influence the observed untreated LDL-C. Inclusion in the sample is the selection node the analysis conditions on. There is no arrow from route to underlying biology. It is conceptual, not an analysis. (Section 5.1.)

**Q131. Distinguish expressivity, recognition and prognosis.**

Expressivity: the phenotype observed among carriers. Recognition: whether that phenotype prompts detection. Prognosis: outcomes under a defined time origin. This chapter examines expressivity as observed through different recognition pathways; its event analyses are secondary. (Section 5.1.)


### Section 5.2 Design, cohorts and provenance

**Q132. Give me the denominators.**

1,904 eligible Welsh LDLR carriers: 815 probands and 1,089 cascade-detected relatives; computable untreated LDL-C in 780 and 889; 235 with no computable value contribute nothing. The residue-matched primary frame: 1,100 carriers, 464 probands and 636 relatives, at 50 residues in 488 pedigrees, median two carriers per pedigree, 48 per cent of pedigrees contributing one carrier. UK Biobank: 3,540 flagged, 2,398 untreated with a measured LDL-C. (Section 5.2, Table 5.1.)

**Q133. What is the provenance hierarchy for the untreated value?**

First, a documented pretreatment value. Second, a measured value in a carrier with no recorded treatment. Third and only then, back-calculation from a treated value by a drug-specific rule. Directly observed in 88.7 per cent of probands and 95.3 per cent of relatives; reconstructed in 11.3 and 4.7. (Section 5.2.)

**Q134. Are the reconstructed rows a random subset?**

No, and that matters. Carriers with reconstructed values had far more prior hard events: 78.4 versus 13.2 per cent in probands, 61.9 versus 4.4 in relatives. Individual reconstruction agreed poorly with documented values, r 0.09. The association of untreated LDL-C with prevalent disease differed by source, odds ratio 1.20 in observed rows against 1.00 in reconstructed. So reconstruction is a cohort-level device, and the observed-only and documented-only restrictions are the checks. (Section 5.2.)

**Q135. Chapter 7 says the registry's drug-specific back-calculation was biased by 4.5 millimoles per litre. Chapter 5 uses a drug-specific rule. Reconcile.**

The two chapters use different rules on different pairs for different targets, and the thesis should bring the agreement figures into one table. What protects Chapter 5 is that its primary contrast does not need the rule: documented pretreatment values only give plus 1.16; directly observed values only give plus 1.09; nine alternative rules give plus 1.00 to plus 1.09. (Sections 5.2, 5.4, 7.2; Appendix 1, item P11.)

**Q136. What does residue matching control, and what does it not?**

It blocks the explanation that groups carried variants in different receptor regions. It does not control the exact substitution, functional severity, polygenic background, lipoprotein(a), metabolic context, calendar era, treatment history or referral practice. Exact-allele matching controls the substitution; within-pedigree comparison holds the familial allele and much of the shared background, but only 236 of 488 pedigrees contain both routes. (Section 5.2.)

**Q137. How much LDL-C variance does the variant explain?**

Little. In the TUDOR analyses the exact variant had an intraclass correlation of about 0.145, leave-one-variant-out 0.127 to 0.149. UK Biobank carriers of the same LDLR variants as Welsh clinic carriers had LDL-C 1.36 millimoles per litre lower, interval minus 2.03 to minus 0.76. And a predicted measure of variant function tracked untreated LDL-C in relatives, rho 0.26, but not in probands, 0.02, a post hoc pattern consistent with selection on the phenotype. (Section 5.2; ref 21, DOI verified 10.1016/j.jacl.2026.06.030.)

**Q138. What does the newborn evidence show?**

In the Norwegian cascade programme, LDL-C did not differ between newborns with null and non-null variants, 1.26 versus 1.18 millimoles per litre, p 0.552. The contrast is within the 61 FH newborns of a 113-newborn sample; my text says 113 and I will correct it. It speaks to screening at birth, not adult phenotype by route. (Section 5.2; ref 85, DOI verified 10.1093/eurheartj/ehaf815.)


### Section 5.3 Primary result and attenuation

**Q139. The primary estimate, with its inference.**

Plus 1.09 millimoles per litre, 95 per cent interval 0.83 to 1.36, p under 0.001, with residue fixed effects and pedigree-clustered standard errors. HC1-robust and a Rademacher wild cluster bootstrap with 1,999 pedigree resamples gave the same interval. These are Paper 14 R4 estimates, not independently re-run for the thesis. (Section 5.3.)

**Q140. Age was the largest threat. What did you do?**

Probands were median 60.2 years, relatives 38.6. Linear age adjustment: plus 0.81, interval 0.52 to 1.09. Natural spline with an age-by-sex interaction: plus 0.78. Adding sex to linear age: plus 0.82. Restricting to ages 40 to 65, prespecified, medians 55 and 52: plus 0.65, interval 0.24 to 1.07, n 457. (Section 5.3.)

**Q141. The within-pedigree result. Is it real?**

Detectable, but imprecise, and the conclusion depends on the variance convention. Plus 0.39, p 0.023, pedigree bootstrap 0.04 to 0.74, directly observed only plus 0.41. Under a degrees-of-freedom correction that counts every pedigree effect, minus 0.06 to 0.85. Appendix B.10 records that the rung has moved across four manuscript versions and that revision R3 crossed zero. What it supports is a residual within-family association the design cannot attribute to referral selection or to biology. I would not claim more. (Section 5.3; Appendix B.10.)

**Q142. Why is the ladder more informative than any coefficient?**

Because its rungs are different operations. 1.09 to 0.81 is adjustment in the same carriers. 0.81 to 0.65 is restriction to a common age band, a change of population. 0.65 to 0.39 is a change of comparison to within-family contrasts, informed by 236 pedigrees. On fixed membership, pedigree effects without age gave plus 0.59 and with age and sex plus 0.40, so the fall is model, not sample. (Section 5.3.)

**Q143. Where is the sex-stratified within-pedigree estimate from the earlier draft?**

Withdrawn in the revision because it could not be refitted on the corrected frame, and the study was not designed or powered for interaction by sex. It was not retained with a stale value. (Section 5.3.)


### Section 5.4 Robustness and the three strata

**Q144. List the robustness rows and their purpose.**

Directly observed only, plus 1.09. Documented pretreatment only, plus 1.16. Inverse-probability-of-observation weighting, plus 1.07, because 4.3 per cent of probands and 18.4 per cent of relatives lacked a computable value. Nine back-calculation rules, plus 1.00 to plus 1.09. Non-HDL cholesterol, plus 1.67 crude, plus 1.44 adjusted. Excluding probands with prior ASCVD, plus 0.96. Age, sex and treatment adjusted, plus 0.86, an over-adjustment. Adding BMI in 380 carriers, plus 0.65. Leave-one-out moved the estimate by at most 0.03 for pedigrees and 0.13 for residues or variants. These re-examine one frame; they are not replication. (Section 5.4, Table 5.2.)

**Q145. Why is the non-HDL difference larger than the LDL-C difference?**

Mostly timing. On the same 974 carriers the LDL-C difference was 1.14; at the same first visit the two differences were 1.64 and 1.49. A small additional remnant component in probands is consistent with, but not shown by, the data. (Section 5.4.)

**Q146. What is the variant-level synthesis and why is its status weaker?**

A random-effects meta-analysis across 71 LDLR/APOB variants with at least three carriers per route: plus 0.91, interval 0.75 to 1.07, I-squared 57 per cent. Its producing output was not located, so it rests on a secondary record, and a later 53-variant rerun gave a different pooled value and is treated as non-equivalent. (Section 5.4; THESIS-AN-02.)

**Q147. Is 5.80, 4.76, 3.87 a matched comparison?**

No. It compares differing allele spectra. Its use is descriptive: the phenotype encountered as detection moves from referral to cascade to sequencing. (Section 5.4.)

**Q148. The first-versus-later reading result: regression to the mean?**

A proband's first pre-treatment reading averaged 2.20 millimoles per litre above the next, interval 1.62 to 2.76, 95 pairs; a relative's 0.91, interval 0.43 to 1.42, 84 pairs; drift with interval not detectable. The empirical reference-change value was 93 per cent against a biological CV of about 7.9 per cent for calculated LDL-C. Compatible with referral-associated measurement selection and regression to the mean, inseparable from biological change or unrecorded treatment. Descriptive, post hoc, and no share of the between-group contrast is attributed to it. (Section 5.4; ref 71, DOI verified.)


### Section 5.5 Cumulative exposure and event associations

**Q149. Did earlier detection mean lower exposure?**

The data cannot say. The manuscript approximated exposure as baseline LDL-C times attained age; in 3,204 UK Biobank carriers with 137 events it added 0.002 to Harrell's C over age, sex and LDL-C. A single-measurement product is not a life-course integral. (Section 5.5.)

**Q150. Were the per-unit LDL-C associations similar by setting?**

Wales, cross-sectional: odds ratio 1.15 per millimole for prior hard ASCVD, 223 events in 1,669; relatives alone 1.14, 63 events. UK Biobank: odds ratio 1.19 for any ASCVD; hazard ratio 1.21 for incident, 137 events; exact-allele P/LP 1.10, interval 0.93 to 1.29, 49 events. No equivalence was tested, the magnitudes are not comparable, and proband analyses are collider-prone. (Section 5.5.)

**Q151. Can I add 1.09 to a cascade relative's value?**

No. It would ignore age, exact allele, treatment provenance and between-variant heterogeneity. The clinical response is to preserve the molecular diagnosis, establish provenance, assess full risk context and follow guidance. (Section 5.5.)


### Section 5.6 Consequences for identification and prediction

**Q152. What is the reporting unit for identification research?**

Model by ascertainment route by comparator spectrum by data availability. Discrimination with calibration and threshold yield, route-stratified where information permits. (Section 5.6.)

**Q153. What is selection leakage?**

A predictor that partly reconstructs the case-definition process: proband status, referral source, pedigree completeness. Whether it is legitimate depends on use; it may prioritise testing within a cascade programme and be tautological in primary care. Gains attributed to route need a locked ablation, route-stratified evaluation and frozen application. (Section 5.6.)

**Q154. What proportion of P/LP carriers meet the Simon Broome LDL-C criterion?**

Among 715 UK Biobank P/LP carriers with no recorded treatment, 84.5 per cent, interval 81.6 to 86.9, were below 4.9. In Wales, 38.1 per cent of probands, 56.9 per cent of relatives and 40.5 per cent of genotype-negative referrals were below it. Untreated LDL-C separated carriers from genotype-negatives with AUC 0.58 on the proband route and 0.76 on the cascade route. (Section 5.6.)

**Q155. How does ascertainment reach Chapter 7?**

The P/LP fraction was 31.4 per cent in UK Biobank and 57.7 per cent in the registry, a gradient in the genetic label itself. A model developed in one and applied to the other inherits that. (Sections 5.6, 7.4.)


### Section 5.7 Contribution and limits

**Q156. State the answer to the route part.**

The same affected LDLR residue carried a higher untreated LDL-C in carriers found through their own phenotype than in relatives found by cascade testing, plus 1.09. The difference narrowed to plus 0.39 within families. The route into care changes which phenotype is seen, not the biology. (Section 5.7.)

**Q157. The limitations that prevent a larger claim.**

Residue matching is positional. The within-pedigree rung rests on 236 pedigrees and its interval depends on convention. Lipoprotein(a) was recorded for 143 of 1,904; apoB was not in this extract; polygenic scores were unavailable. Era and assay heterogeneity. Complete-case covariate models. BMI in 380. UK Biobank volunteer-selected and ancestry-restricted. And the Welsh carriers overlap with the Chapter 4 registry sample by an unestablished number. (Section 5.7.)

**Q158. Which ancestry restriction?**

The thesis names an ancestry restriction without defining it. I need to state the group used and confirm consistency across Chapters 4 and 7. (Section 5.7; Appendix 1, item P9.) ---


## Chapter 6: Identify. Can carriers be found after treatment?


### Section 6.1 Identification after treatment and selection

**Q159. What is the contemporary identification question, and how does it differ from Simon Broome's?**

Not whether an untreated proband satisfies a score. Among people already in a lipid, cardiology or genetic-testing pathway, whom should a service prioritise for confirmatory testing when the untreated phenotype is no longer observable? (Section 6.1.)

**Q160. Quantify the failure of a single LDL-C gate.**

Among UK Biobank ClinVar P/LP carriers with an untreated LDL-C, 80.2 per cent were below the 4.9 Simon Broome trigger, and a 4.9 gate would have required 404 people sequenced per carrier found. That is the failure mode a multivariable triage score is meant to reduce. (Section 6.1; ref 21, DOI verified.)

**Q161. What would weaken the chapter's hypothesis?**

Absence of a paired discrimination advantage; marked family leakage; poor frozen application; failure to reproduce the recovered score. Calibration failure would not erase discrimination but would prevent reading the score as a probability. (Section 6.1.)


### Section 6.2 Model specification and evaluation states

**Q162. Specify TUDOR exactly.**

Elastic-net logistic regression, mixing parameter 0.5, C equal to 1.0, saga solver, seed 20260518. Eleven declared inputs in fixed order: reconstructed pretreatment LDL-C, HDL-C, triglycerides, total cholesterol, non-HDL-C, age, sex, treatment status, family history of cardiovascular disease, a diabetes-by-LDL interaction, and premature cardiovascular events. Ten contribute; family history is zero-variance with coefficient zero. No proband-route or ascertainment term. (Section 6.2.)

**Q163. The published paper says ascertainment-aware. The object has no route term. Explain, and tell me why there is no corrigendum.**

The publication's title says ascertainment-aware and its abstract says the model encodes index-referral versus cascade-relative ascertainment. When the frozen executable pipeline was recovered and reproduced by score readback, its eleven inputs contained no such term. Every TUDOR result in the thesis refers to the recovered object. No Index-Effect increment is reported, and a number that circulated in earlier programme material is omitted because no producing object supports it. The discrepancy is disclosed in the thesis only. No corrigendum has been requested, and that was my decision. On reflection the order is wrong: a reader of the paper attributes to the model a feature the evaluated object does not contain. I will write to the journal before final submission. If pressed, Does it invalidate the AUC?: No. Readback reproduces the stored predictions within tolerance, so 0.760 is the performance of the object that exists. It invalidates the route-aware claim, not the ranking claim. (Section 6.2, Box 6.1; ref 21, DOI verified 10.1016/j.jacl.2026.06.030; see Appendix 1, item P18 on the exact abstract wording.)

**Q164. Your Atherosclerosis Plus abstract claims dual external validation in 4,028 genetically confirmed cases. Is that consistent with the thesis?**

No, and the thesis corrects it. The abstract draws on the same two cohorts as the chapter and is cited as related work, not independent validation. The 4,028 figure appears nowhere else in the thesis and I need to reconcile it against the abstract before quoting it. (Section 6.2; ref 82, not indexed in PubMed; Appendix 1, item P19.)

**Q165. Name the evaluation frames and say why they cannot be merged.**

Welsh complete-case head-to-head, 1,274 with 311 carriers. Full deduplicated registry, 4,570 with 887. Reciprocal geographic splits, 3,099 with 554 and 1,471 with 333. UK Biobank lipid-clinic-eligible, 49,427 with 921. UK Biobank local sensitivity, 43,594 with 652. Whole UK Biobank, 501,936 with 3,540. And a 649-pair reconstruction resource. Different eligibility, comparator-availability and carrier-definition rules; merging them would create a result no analysis produced. (Section 6.2, Table 6.1.)

**Q166. Which paired test, and what is wrong with it here?**

DeLong's test for correlated AUCs, used in the local UK Biobank paired comparisons and the reconstruction sensitivity. It treats participants as independent; in a family-structured registry that may understate uncertainty. The Welsh head-to-head differences were reported without a named paired test or interval, which is a gap in my own source report. (Section 6.2; ref 87, DOI verified 10.2307/2531595.)


### Section 6.3 Welsh head-to-head discrimination

**Q167. The result.**

In 1,274 participants, 311 carriers: TUDOR AUC 0.760, interval 0.727 to 0.790. Reconstructed eDLCN 0.652, FAMCAT 0.600, Simon Broome 0.569, MEDPED 0.524. Net reclassification improvement against FAMCAT plus 0.204, interval 0.141 to 0.266. (Section 6.3, Table 6.2.)

**Q168. What does that ordering answer, and what does it not?**

Given the information represented for every tool under the same complete-case restriction, which score ranked carriers above non-carriers more consistently? It does not show TUDOR beats the criteria applied with observed untreated lipids, full pedigrees and examination. And 0.760 is not a clinical threshold. (Section 6.3.)

**Q169. Why is NRI secondary?**

Because it is sensitive to prevalence, score distributions, category definitions and miscalibration, and a large NRI against a weak comparator need not mean a better decision. (Section 6.3; ref 88, DOI verified 10.1093/aje/kwx374.)

**Q170. Does TUDOR beat pretreatment LDL-C alone?**

In one defined local pipeline: 0.747 against 0.705, plus 0.043 by DeLong. An earlier analysis with a different 450-carrier definition gave plus 0.017 without statistical separation. The increment beyond LDL-C is carrier-definition-sensitive, and I do not claim a stable increment. (Section 6.3.)


### Section 6.4 Comparator fairness and the ascertainment boundary

**Q171. Make the case against your own comparison.**

eDLCN ordinarily combines LDL-C with vascular history, tendon xanthomata, corneal arcus, pedigree and molecular evidence; the electronic reconstruction lacked several, and unrecorded findings scored as absent. FAMCAT was built for primary-care records and ran on a specialist-registry subset. Simon Broome and MEDPED depend on untreated cholesterol and physical signs. The comparison is asymmetrical. It is still clinically relevant because the electronic-data problem is real; it prohibits the claim that TUDOR beats a fully informed specialist. (Section 6.4; refs 43, 74 to 76, DOIs verified.)

**Q172. What is the counterweight?**

Akyea 2026. In an Australian tertiary clinic, 885 referred for genetic testing, 267 with an FH-causing variant, DLCN scored with examination, pedigree and untreated lipids had an AUROC of 0.816 against FAMCAT's 0.748. Scored with that information, the criteria beat an electronic algorithm. TUDOR has not been compared with fully scored DLCN. (Section 6.4; ref 45, DOI verified 10.1016/j.jacl.2026.07.009.)

**Q173. Was TUDOR compared with the Welsh service's own criteria?**

No. The Haralambos 2015 criteria, which modify DLCN, are TUDOR's local predecessor, and the increment over them is unknown. (Section 6.4; ref 44, DOI verified 10.1016/j.atherosclerosis.2015.03.003.)

**Q174. Can you partition TUDOR's advantage into treatment-awareness and route?**

No. There is no executable route term and no full-versus-minus-route ablation. Descriptive registry AUCs of 0.760 in probands and 0.794 in relatives show no obvious loss in relatives, but case-mix differences prevent attribution. (Section 6.4.)

**Q175. Family leakage?**

Unquantified. No family-aware optimism estimate with a located producing record exists, and no fold manifest has been reconstructed. Expected direction: optimism. (Section 6.4.)


### Section 6.5 Frozen application and the intended-use boundary

**Q176. Transport ladder, in numbers.**

Head-to-head 0.760. South Wales to rest of Wales 0.732, interval 0.707 to 0.757. Rest of Wales to South Wales 0.770, interval 0.739 to 0.797. Pooled frozen Wales 0.746. Frozen UK Biobank lipid-clinic-eligible 0.669, interval 0.650 to 0.687. UK Biobank target refit 0.756 apparent. Whole UK Biobank 0.631. (Section 6.5, Table 6.3.)

**Q177. Which is the most consequential number?**

0.669. The frozen Welsh score did not preserve its source discrimination in a cohort that differed in age, treatment, prevalence, measurement and volunteer selection. It does not say which difference caused the loss. And it must not be displaced by the refitted 0.756. (Section 6.5.)

**Q178. Why is 0.756 not external validation?**

Coefficients were re-estimated in UK Biobank; relative predictor weights changed; target outcomes informed the model. That is local development. It shows the feature architecture contained reweightable information. It does not say a new service should expect 0.756. (Section 6.5.)

**Q179. What does 0.631 mean?**

The intended-use boundary. TUDOR was not built for indiscriminate population screening, and in an unselected local analysis it did not exceed LDL-C alone. A pathway-aware model can be useful after the pathway exists and add nothing before it. (Section 6.5.)

**Q180. Subgroups?**

Five-gene-plus-CNV outcome 0.779. Type 2 diabetes 0.648 in Wales and 0.642, interval 0.559 to 0.724, in UK Biobank: a metabolically difficult setting, not a subgroup licence. A statin-naive 0.801 was reported without a verified denominator or interval and carries no strong conclusion. (Section 6.5.)


### Section 6.6 Reconstruction as a measurement model

**Q181. How large is the reconstruction error?**

Mean absolute error 1.20 millimoles per litre in 649 paired participants, correlation 0.32. Clinically material near thresholds. It reflects adherence, dose, combinations, response and timing that a deterministic factor cannot recover. (Section 6.6; Appendix D.)

**Q182. Why did the AUC not move across reconstruction methods?**

Across dose-specific, class-level and fixed-factor approaches the AUC range was under 0.02, paired p above 0.3. Ranking tolerates some systematic error because other predictors contribute and a shared shift preserves order. That is a statement about ranking, not agreement. (Section 6.6.)

**Q183. The clinical boundary of a TUDOR score?**

A high score may justify priority for confirmatory testing; it does not establish a genotype. A low score should not overrule a known familial variant, tendon xanthomata, a compelling pedigree or expert assessment. (Section 6.6.)


### Section 6.7 Calibration, workload and utility

**Q184. Is TUDOR calibrated?**

Not shown to be. In a local UK Biobank stream, slope 1.257, intercept minus 3.03, Brier 0.076, retained aggregate: substantial overprediction in a lower-prevalence setting. Those values must not be transferred to the locked 49,427 transport cohort, for which external calibration was not available. TUDOR's output is a ranking score. (Section 6.7; THESIS-AN-03b.)

**Q185. The operating points.**

In the local 43,594 with 652 carriers: Youden point refers 10,433, captures 409, sensitivity 0.627, PPV 3.9 per cent, 25.5 sequenced per carrier. Top 10 per cent refers 4,360, captures 266, PPV 6.1, 16.4 per carrier. A high-sensitivity point refers 36,417, captures 619, sensitivity 0.949, PPV 1.7. A statistical optimum encodes none of the harms. (Section 6.7, Table 6.4.)

**Q186. Why do numbers needed to screen differ tenfold?**

2.0 in Wales, 15.1 in UK Biobank: prevalence and pathway, not the model. A Welsh threshold cannot be transferred mechanically. (Section 6.7; THESIS-AN-03a.)

**Q187. Why no decision curve for TUDOR?**

Because a decision curve needs valid target probabilities, and TUDOR's are not calibrated for the transport cohort. It is the appropriate next method, after calibration, and even then it would not prove improved outcomes. (Section 6.7; ref 33, DOI verified.)


### Section 6.8 Reporting appraisal and gates

**Q188. Why is there no PROBAST+AI rating for TUDOR when you gave one for CALON-5?**

Because the signalling-question assessment was not completed for TUDOR; Table 6.5 gives narrative considerations only. That is inconsistent with Appendix E.4 and I will complete it. (Section 6.8; Appendix 1, item P12.)

**Q189. What blocks independent deployment?**

Incomplete target calibration, comparator asymmetry, absent prospective threshold evaluation, absent impact evidence, and the description-object discrepancy. Reproducibility is no longer the block: the object is recovered and reads back within tolerance. (Section 6.8.)

**Q190. The next validation stage, in order.**

Freeze the pipeline with its software environment and readback tests. Fix carrier definition, predictor timing, missing-data policy and intended pathway, with preprocessing inside family-safe resampling. Frozen evaluation in an independently assembled specialist-triage cohort, reporting calibration-in-the-large, slope and curves separately from discrimination. Prespecified thresholds with workload, decision curves and missed-carrier consequences. Then prospective impact. (Section 6.8.)


### Section 6.9 Contribution and transition

**Q191. State the answer to the treatment part.**

Reconstructing the untreated phenotype improved the ranking of carriers against reconstructed criteria in the development setting, and less so when the score was applied unchanged elsewhere. Identification is not prognosis: a high probability of carrying a variant does not say who has an event first. (Section 6.9.) ---


## Chapter 7: Predict. Does risk ranking survive a change of cohort?


### Section 7.1 The question

**Q192. Two meanings of the same thing for a risk score?**

The order in which the model places patients, discrimination, and the absolute risk it assigns them, calibration. The chapter keeps them apart throughout. (Section 7.1.)

**Q193. Why can a within-FH model never inform whether to treat?**

Because lipid lowering is already indicated by the diagnosis. It could inform intensity, timing in young relatives, and allocation of add-on therapy or imaging. (Section 7.1; refs 1 to 3, DOIs verified.)

**Q194. You call both models high risk of bias on your own appraisal. Why, in four points?**

CALON-5's terms were chosen with knowledge of its transported registry performance. The two-term model was selected after those results were known. The registry fit rests on few events. And neither cohort is untouched by the programme's development. (Section 7.1; Appendix E.4.)


### Section 7.2 Model identities

**Q195. Specify CALON-5.**

Cox proportional hazards, no stated penalty. Five terms: age, male sex, hypertension defined as recorded antihypertensive medication at baseline, diabetes, and the log of untreated non-HDL-C over HDL-C. The hypertension definition replaced a blood-pressure-or-medication definition on 24 September 2026, before any refit result was seen; the earlier definition is a sensitivity. (Section 7.2.)

**Q196. Where does the untreated lipid value come from?**

The registry supplies a recorded pre-treatment value. UK Biobank holds none, so for the 24.7 per cent of biobank carriers on therapy, untreated non-HDL-C is back-calculated as 5.2299 plus 0.3345 times the recorded value, from 684 within-person registry pairs, cross-validated mean absolute error 1.25, R-squared 0.0975. A drug- and dose-specific equation was not used: dated, dosed prescriptions existed for 27.6 per cent of treated carriers, and the registry's drug-specific rule was biased by 4.5. (Section 7.2.)

**Q197. R-squared of 0.0975 and a slope of 0.33. Defend using it.**

It explains almost no individual variation and compresses the treated between-person spread; it was fitted mostly on non-carrier pairs. It is used for one purpose: group-level targeting of the untreated quantity guidelines specify, so that the two cohorts model the same nominal exposure. It is not accurate for an individual, the error is plausibly differential, and the registry uses recorded values, so the cohorts measure exposure differently. All of that is stated. (Sections 7.2, 7.10; Appendix D.3.)

**Q198. What is the two-term model and why does it exist?**

Age and the same lipid ratio. Because CALON-5's terms were chosen knowing their registry performance, selection was repeated on UK Biobank data alone, blind to the registry. It chose these two, the smallest of 25 indistinguishable subsets. In 200 bootstraps age was retained 83 per cent of the time and the lipid ratio 32 per cent, so its identity is a convention. Optimism-corrected concordance 0.612. (Sections 7.2, 7.3.)

**Q199. Prespecified after the results were known is a contradiction.**

Both words are true of different things, and I should not hide behind the phrase. The selection rule was written into the revision protocol of 23 September 2026 before it was run, but after CALON-5's transported results were known. The honest label is 'development-only selected', which is why neither model is primary. (Section 7.2, Table 7.1.)

**Q200. Why exclude Lp(a), apoB and imaging?**

Not for irrelevance. The question is whether routinely recorded variables support ranking where specialised inputs are unavailable. A model with fewer inputs computes for more patients while omitting information; a richer one performs where its inputs exist and cannot be computed where they were never collected. (Section 7.2.)


### Section 7.3 Development and internal validation

**Q201. Events per parameter, both directions.**

Primary: UK Biobank P/LP, 888 complete cases, 74 events, 14.8 per parameter. Reverse: registry P/LP, 469 complete cases, 28 events, 5.6 per parameter, below the prespecified floor of ten; exploratory. The all-carrier reverse fit, with 52 events, is the one interpreted. (Section 7.3.)

**Q202. Why is this cross-cohort evaluation and not external validation?**

The registry is not untouched: the untreated-exposure equation used its lipid pairs, and CALON-5's terms were chosen with knowledge of its transported registry performance. (Section 7.3.)

**Q203. Internal validation numbers.**

Whole-procedure bootstrap: apparent concordance 0.6972, optimism 0.0209, corrected 0.6763; calibration slope corrected to 0.901. No term showed evidence against proportional hazards, 0 of 10 tests, a check with little power on 52 registry events. (Section 7.3.)

**Q204. Did the hazard ratios agree across cohorts?**

In the all-carrier fits all five agreed in direction. The lipid term: 2.33 per log unit in UK Biobank, interval 1.62 to 3.35, and 2.40 in the registry, interval 1.34 to 4.31. At the same covariates registry carriers had 1.33 times the biobank hazard, interval 0.83 to 2.13. (Section 7.3.)


### Section 7.4 Cohorts, flow and missingness

**Q205. What did standardising to ClinVar P/LP cost?**

UK Biobank fell from 3,209 carriers and 289 events to 1,009 and 87. The registry from 1,639 genotyped attendees and 101 events to 945 and 56. The asymmetry is the pathogenic fraction, 31.4 against 57.7 per cent, itself an ascertainment signal. At two or more ClinVar review stars, 69.0 per cent of biobank P/LP variants qualified against 97.7 per cent of registry variants; zygosity is recorded in neither. (Section 7.4.)

**Q206. Why is the registry complete-case set not missing at random, and what did you do?**

Requiring all five inputs kept 800 of 1,639 all-carriers and 469 of 945 P/LP carriers, driven by unrecorded diabetes: recorded for 54 per cent, hypertension for 99. Those without a recorded diabetes status had about half the hazard, age- and sex-adjusted 0.40, interval 0.21 to 0.77. Missingness was associated with the outcome, so MAR is doubtful, and I did not impute. The set is selected toward longer-observed, clinically assessed patients; that is limitation L15. (Section 7.4.)

**Q207. Endpoint composition and censoring in the registry.**

Of 101 events: 52 infarction or acute coronary syndrome, 16 percutaneous and 19 surgical revascularisations, 8 angina, 6 TIA or stroke. Event times from recorded ages. Censoring survival 0.83 at five years and 0.52 at ten, so ten-year estimates lean on the weights. (Section 7.4.)


### Section 7.5 Discrimination

**Q208. The headline.**

Developed in UK Biobank and applied unchanged to the registry, ten-year concordance 0.699, interval 0.611 to 0.790, time-dependent AUC 0.718. Reverse: 0.698, interval 0.609 to 0.778, AUC 0.713. Across seven estimable cells, 0.698 to 0.722, every interval overlapping. (Section 7.5.)

**Q209. How much is age and sex?**

Most of it. Age and sex alone reached 0.640 to 0.669 across the four cells. On identical registry rows CALON-5 beat the two-term model by 0.043, interval 0.005 to 0.083. No interval exists for CALON-5 over age and sex, so whether that increment is demonstrated cannot be judged. (Section 7.5; limitation L14.)

**Q210. The 4.9 cut-point as a ranker?**

Close to chance, concordance 0.552 to 0.593; CALON-5 exceeded it by 0.111 to 0.170 with intervals excluding zero in all four cells. A severity marker is not a risk ranker. (Section 7.5.)


### Section 7.6 Calibration

**Q211. The headline.**

Absolute risk did not transport and failed in opposite directions. Into the registry: 3.05 per cent predicted against 5.98 observed at ten years, observed-to-expected 1.96. Into UK Biobank: 6.76 against 4.73, ratio 0.70. None of seven ratios lay near one. Slopes 0.83, interval 0.47 to 1.22, into the registry; 0.55, interval 0.36 to 0.72, out of it. (Section 7.6.)

**Q212. What fixed it, and is the fixed model validated?**

Re-estimating the baseline fixed the level in all seven cells; the slope needed a second correction only out of the registry. The two-term model needed no update out of the registry and both corrections into it. Updated results describe adaptation to the receiving cohort. They are not external validation. (Section 7.6.)

**Q213. Why is the gap unsurprising?**

By attained age 70, cumulative MACE incidence was 9.0 per cent in UK Biobank carriers, interval 6.4 to 11.9, and 24.4 per cent in registry carriers, interval 19.3 to 30.6. A model inherits its development cohort's baseline hazard. Relative effects travelled; absolute risk did not. (Section 7.6.)


### Section 7.7 What the recorded value costs guideline grading

**Q214. The result, exactly.**

In 261 genotype-positive registry patients with both a recorded pre-treatment and an on-treatment LDL-C, all measured, none reconstructed: median 6.1 untreated, 3.6 on treatment. 173 of 261, 66.3 per cent, sat in a lower ESC/EAS stratum on the recorded value; 29.5 unchanged; 4.2 higher. 193 were at or above 4.9 untreated and 57 on treatment; 145 of the 193, 75.1 per cent, would not have been identified on the recorded value. (Section 7.7.)

**Q215. Its weakness?**

It is estimated only in patients who had a documented untreated value, and who has one is not random; they may be earlier-identified or more severe. So it is not a service-wide rate. Its strength is that it needs no reconstruction. (Section 7.7; limitation L16.)

**Q216. What is new relative to Singh 2026 and Mancini 2026?**

Singh applied correction factors to trial participants; Mancini imputed untreated LDL-C in homozygous FH. The increment here is a direct measurement, in genotyped heterozygous carriers, of how often the recorded value under-grades severity. (Section 7.7; refs 92, 93, DOIs verified 10.1016/j.jacl.2026.03.018 and 10.1016/j.atherosclerosis.2025.120590.)

**Q217. Goal attainment.**

Below 1.8: 7 of 594 treated UK Biobank participants, 1.2 per cent, and 5 of 1,194 treated registry patients, 0.4 per cent. Of the 261 with a pair, 26.8 per cent met the NICE 50 per cent reduction, 5.0 reached below 1.8, 3.1 met both; 88.6 per cent of those achieving the reduction still missed the goal. Attainment depended on the starting value, 39.4 above the median baseline against 12.9 below, as regression to the mean predicts. No goal contrast was estimable for events. (Section 7.7; ref 73, NICE CG71, not DOI-indexed.)

**Q218. Two measurement properties?**

The triglyceride term's median share of Friedewald LDL-C rose from 9 per cent at the severity cut-point to 47 per cent below 1.4, so the lowest goals are judged on the least reliable number. And of 51 registry patients at the 2.6 goal with a pre-treatment value, 54.9 per cent had been at or above 4.9 untreated. (Section 7.7.)

**Q219. Does the cut-point separate risk?**

Yes, more sharply in the clinic: event-rate ratio above versus below 4.9 was 2.97 in the registry and 1.61 in UK Biobank. That is the attenuation ascertainment predicts. (Section 7.7.)


### Section 7.8 Comparison with published scores

**Q220. Rules of the comparison.**

Each score implemented as published with its own variables, coding and missing-input rules, Lp(a) converted to the units it takes. Compared with CALON-5 on rows where both are computable, one set of family-grouped resamples so every difference is paired. No non-inferiority verdicts, because no margin exists. A comparison resolves only when its interval excludes zero. (Section 7.8.)

**Q221. The results.**

Standardised P/LP, ten years: FH-Risk-Score into the registry, 402 rows, 16 events, below the floor, point estimate minus 0.039 favouring the comparator. FH-Risk-Score out, 759 rows, 31 events, plus 0.029, interval minus 0.065 to 0.120. Montreal into the registry, 449 rows, 20 events, plus 0.052, interval plus 0.0003 to 0.104, resolves by the rule and did not under the earlier hypertension definition. Montreal out, 888 rows, 42 events, minus 0.006. SAFEHEART-RE into the registry, 26 rows, one event, not estimable. SAFEHEART-RE out, 725 rows, 31 events, plus 0.076, interval minus 0.041 to 0.189. CALON-5 higher in 3 of 4 estimable comparisons; one resolved, fragile. (Section 7.8, Table 7.2.)

**Q222. Did the published scores calibrate?**

No. On UK Biobank P/LP rows SAFEHEART-RE under-predicted, observed-to-expected 1.72, interval 1.22 to 2.33, slope 0.52; the FH-Risk-Score over-predicted, 0.34, interval 0.23 to 0.47, slope 0.78. (Section 7.8; refs 53, 56, DOIs verified 10.1161/CIRCULATIONAHA.116.024541 and 10.1161/ATVBAHA.121.316106.)

**Q223. How does this sit against McKay and the Australian work?**

McKay transported SAFEHEART-RE into English primary care by imputing unrecorded predictors: Harrell's c 0.67, interval 0.61 to 0.72, with substantial miscalibration. In Australian cohorts with measured inputs, SAFEHEART-RE validated for incident events and the Montreal, Combined and FH-Risk scores discriminated prevalent disease. One cohort, one direction each. Mine is paired, common rows, both directions, and finds most comparisons too small to resolve. (Section 7.8; refs 28, 29, 94, 95, DOIs verified.)


### Section 7.9 Risk enhancers

**Q224. Does apoB add to CALON-5? Reconcile with Chapter 4.**

In 888 P/LP carriers with 74 events, no enhancer's interval excluded one: apoB/LDL-C discordance 0.79, interval 0.38 to 1.62; Lp(a) at or above 125 nanomoles per litre 0.88, interval 0.40 to 1.94. In all carriers the polygenic score was 1.37 per SD, interval 1.12 to 1.68, with concordance plus 0.012, interval minus 0.000 to 0.040. Chapter 4 asks about association at fixed LDL-C; this asks about increment over a model already using the untreated lipid ratio, age, sex and comorbidity, in a smaller frame. Different questions; moderate effects are not excluded. (Section 7.9.)


### Section 7.10 Treatment reconstruction and measurement error

**Q225. Can stable concordance validate the reconstruction?**

No. If one correction moves most treated participants similarly, order changes little even when each value is wrong. Ranking stability cannot validate a threshold or a probability. Cumulative-exposure constructs did not improve on a single concentration: concordance difference 0.007, point estimate only. (Section 7.10.)

**Q226. What would be the stronger design?**

Documented pre-treatment samples, serial lipids, time-linked prescribing, dose and adherence; and where untreated values remain unobserved, a distribution of plausible values propagated through fitting and calibration. (Section 7.10; Appendix D.4.)


### Section 7.11 Calculability and the decisive study

**Q227. Why is calculability a result in its own right?**

SAFEHEART-RE requires Lp(a) and body-mass index and publishes no missing-input rule; after merging the laboratory export Lp(a) was available for 185 standardised registry patients but SAFEHEART-RE was computable for 26, with one event, because BMI was the binding input. The FH-Risk-Score publishes a zero-when-unavailable rule for Lp(a). CALON-5 itself was computable for 469 of 945 registry carriers because diabetes was so often unrecorded. The thesis reports what a service can compute with the data it holds. (Section 7.11, Table 7.3.)

**Q228. Why not impute SAFEHEART-RE's inputs, as McKay did?**

Because the score publishes no missing-input rule, so imputing from derivation means would create an adapted score and label it as the published one. Possible, deliberately not done. (Section 7.11; ref 28, DOI verified 10.1016/j.atherosclerosis.2022.07.011.)

**Q229. The size of a decisive study.**

Held at observed differences, comparisons out of the registry would resolve at about 71 events against SAFEHEART-RE, 312 against the FH-Risk-Score, and 6,656 against Montreal. For calibration-in-the-large within plus or minus 0.1, about 362 ten-year events; for the slope, 248 to 1,285. A genotype-confirmed cohort untouched by either model's development. (Section 7.11.)


### Section 7.12 Decision curves

**Q230. Do they show utility?**

Frozen CALON-5 exceeded both treat-all and treat-none at 13 of 30 thresholds into the registry and 10 of 30 out of it; 18 and 7 after recalibration; the cut-point alone at 6 and 4. Net benefit undefined at six thresholds into the registry. Fewer than half the thresholds: net benefit does not follow from discrimination. Appendix A records 12 and 9; one is stale and I will correct it. (Section 7.12; Appendix 1, item P1; ref 33, DOI verified.)

**Q231. What is the credible action in FH?**

Not a treatment threshold, because lipid lowering is already indicated. Added review, imaging or phenotyping. A low score must never withdraw indicated therapy. (Section 7.12.)

**Q232. Reclassification against the cut-point?**

47.5 per cent of registry and 35.0 per cent of UK Biobank carriers changed category. Frozen in the registry, the model placed 50.1 per cent lower than the cut-point and 0.6 higher; patients who had an event were moved down more often than up, event reclassification minus 0.61, not prespecified. That is the reclassification face of the under-prediction. (Section 7.12.)


### Section 7.13 Effect modification

**Q233. Sex and age?**

Interactions with the lipid term: none of nine estimable excluded one; sex by lipid 0.958, interval 0.716 to 1.282. For the cut-points the ratio of hazard ratios per decade was below one at 12 of 14 thresholds, so cut-points discriminate hardest in the young. Sex modified 2 of 10, consistent with chance. No multiplicity adjustment. (Section 7.13.)


### Section 7.14 The CALON-C lineage

**Q234. What does CALON-C add that CALON-5 does not?**

A sharper version of the same warning. Its frozen nine-term Welsh concordance of 0.7252 was matched by age and sex; the other seven terms added plus 0.0008. Its Welsh calibration failed, observed-to-expected 4.500 at five years and 4.374 at ten, repaired only by re-estimating the baseline. (Section 7.14; Appendix B.)

**Q235. Why is it lineage rather than evidence?**

The exact historical fitted object and the missing-data pathway are unavailable; endpoint dates were unambiguous for 147 of 289 events and a stress test dropped concordance from 0.708 to 0.660; competing deaths had to be reconciled from 193 to 224. Its absolute-risk, Brier and decision-curve outputs are withheld. (Appendix B.4.)


### Section 7.15 Answer

**Q236. State the answer to the cohort part.**

Ranking survived: concordance 0.699 and 0.698 in the two directions. Absolute risk did not: observed-to-expected 1.96 and 0.70, and the published scores miscalibrated the same way. What restores it is local re-estimation of baseline risk; what would settle it is a cohort with several hundred events. (Section 7.15.)

**Q237. In one sentence, why is neither model clinical?**

Because most of the discrimination was demographic, most comparisons were too small to resolve, and neither has calibration, utility or impact evidence in a new setting. (Section 7.15.) ---


## Chapter 8: Discussion


### Section 8.1 The answer

**Q238. Give the answer to the central question in one breath.**

An LDL-C value means the same thing only when five conditions also match: the particles carrying it, the route by which the person was found, the treatment already given, the cohort a model was built in and, in women, reproductive stage. Where any differs, the value needs something added before it can support the decision. (Section 8.1.)

**Q239. Where do the guidelines stand, exactly?**

The 2026 ACC/AHA multisociety guideline rates apoB measurement COR 2a, LOE B-NR, in adults on lipid-lowering therapy, and COR 2b in adults not on therapy including to characterise inherited disorders. The 2025 ESC/EAS focused update makes no apoB recommendation and treats LDL-C as a direct cause of ASCVD. Nothing here challenges that. (Section 8.1; refs 2, 3, DOIs verified 10.1161/CIR.0000000000001423 and 10.1093/eurheartj/ehaf190.)

**Q240. Table 8.1 says the Objective 4 absolute-risk falsifier was met. Was it?**

No, the wording is wrong. The falsifier was calibration that transports without updating. Calibration did not transport, so the falsifier is not met. And the ranking falsifier, no better than age and sex, cannot be cleanly judged because no interval exists for that increment. I will correct the row. (Table 8.1; Appendix 1, item P4.)

**Q241. Is the answer interventional?**

No. No chapter tested a change in care. It establishes where a single recorded value misleads, by how much in these cohorts, and what information would correct each error. (Section 8.1.)


### Section 8.2 Principal findings by task

**Q242. Why the seven-step template in every subsection?**

So that finding, literature, explanation with its tier, competing explanation, boundary, implication and resolving study cannot be quietly omitted for any task. The cost is a mechanical read; the kit review suggests connecting them into argued paragraphs, and I agree. (Section 8.2; K012-049.)

**Q243. For each task, the strongest competing explanation?**

Particles: the residual summarises metabolic comorbidity rather than particle number. Route: polygenic enrichment among probands, and era. Treatment: TUDOR's advantage reflects how weak the electronic reconstructions were. Cohort: differences in treatment era, access and event capture between a volunteer cohort and a service, nothing FH-specific. I accept each as plausible. (Sections 8.2.1 to 8.2.4.)

**Q244. Relative effects travel, baseline risk does not. Did you test that?**

Partially. In the all-carrier fits the five hazard ratios agreed in direction and the lipid term was 2.33 and 2.40 in the two cohorts; but coefficient or interaction transport was not formally tested, and the kit flags the sentence as an inference beyond the test. The supported statement is that ranking was comparable while baseline incidence differed by attained age 70, 9.0 against 24.4 per cent. (Sections 7.3, 7.6, 8.2.4; K012-034.)


### Section 8.3 One value, three uses

**Q245. The three uses and the condition that breaks each.**

Grading severity: broken by treatment; restored by a documented untreated value. Finding carriers: weakened by treatment and by route; partly restored by treatment-aware reconstruction with genetic testing, locally until tested prospectively. Ranking risk: ordering survived a change of cohort, absolute risk did not; restored by local baseline re-estimation and, in women, sex-specific calibration. (Section 8.3.)

**Q246. You say the errors are systematic and predictable in direction. Is that an overclaim?**

For individuals, yes, and the kit flags it. At group level the directions observed were consistent: under-grading after treatment, milder phenotypes in cascade and genotype-first carriers, and miscalibration on transfer. Random error and unvalidated individual remedies remain; I will restrict the sentence to group-level patterns. (Section 8.3; K012-051.)

**Q247. Is 4.9 mmol/L an FH severity cut-point?**

The guidelines use 4.9 as a severe-hypercholesterolaemia and diagnostic-suspicion threshold and, in the ESC/EAS bands, as a stratum boundary; whether it is properly called FH severity grading needs the exact guideline tables. The empirical result stands regardless: 75.1 per cent of registry carriers above it untreated were not above it on the record. I will reword the category. (Sections 1.1, 7.7, 8.3; K012-007.)

**Q248. The laboratory value was accurate. How do you know?**

I do not know analytically; the sentence asserts more than the evidence. The claim I can support is that the value changed its meaning, not that its assay accuracy was established. I will remove the accuracy assertion. (Section 1.1; K012-008.)

**Q249. Linked genotyped cohorts with common definitions. Are the definitions common?**

Only partly, and the kit is right to press. The two data sources are common; the genetic definitions differ by chapter, and participant overlap between Chapters 4 and 5, and between UK Biobank frames, is not established. The honest description is linked heterogeneous studies in two cohorts with a common evidence rule, and I will build the overlap and definition matrix. (Sections 1.6, 8.3; K012-013.)


### Section 8.4 Beyond FH: sex and reproductive stage

**Q250. State the finding and its design.**

273,036 UK Biobank women across the menopausal transition, age-matched male reference. At age 58 an apoB of at least 1.2 grams per litre corresponded to a ten-year MACE risk of 10.88 per cent in men and 4.34 in postmenopausal women. A sex-neutral internal model over-predicted in women, observed-to-expected 0.65 to 0.78, and under-predicted in men, 1.19. Population cohort, self-reported status, cross-sectional stages, internal model. (Section 8.4; ref 50, in press, not indexed.)

**Q251. Isn't that just lower baseline hazard in women at 58?**

Yes, and that is the point being made: a marker's meaning for absolute risk depends on the baseline it is read against. It is a familiar age-sex baseline difference, not evidence of a sex-specific lipid effect or a menopause mechanism, and the section must not be read as either. (Section 8.4; K012-005, 053.)

**Q252. Da Roza found a larger rise in monogenic FH; you found a smaller percentage rise in LDLR carriers. Contradiction?**

Different metrics from different starting levels: an absolute rise in a small FH cohort against a percentage rise, plus 11.5 versus plus 15.3, from a higher base in a population cohort. Not directly comparable, potentially compatible. I will not claim more than that. (Section 8.4; ref 24, DOI verified 10.1016/j.atherosclerosis.2025.120587; K012-052.)

**Q253. Has sex-specific calibration been validated for FH?**

No. The recommendation to use sex-specific absolute risk in interpretation comes from a population cohort. The FH resolving study would need within-person, final-menstrual-period-anchored lipids in molecularly characterised women. (Section 8.4; K012-048.)


### Section 8.5 What the thesis adds

**Q254. Two layers of contribution?**

What the studies established individually: apoB carries information at fixed LDL-C in genotyped carriers; treatment-aware triage ranks better than electronic reconstructions in one Welsh comparison; the route contrast at matched residues; ranking transports while calibration does not; the same cut-point differs by sex. What the thesis adds: the five measured together with one evidence rule, the identification of the one condition under which the value does transport, and a named repair for each. (Section 8.5.)

**Q255. Is claim-to-evidence governance a contribution?**

It is good practice, not a new method, and the kit review says so. I present it as a strength that supports trust, not as doctoral novelty. (Section 8.5; K012-015.)

**Q256. What would refute the framework?**

A frozen residual failing under harmonised assays in an independent FH cohort. No residue-matched route difference in an independently designed cohort. Frozen TUDOR failing against fully scored DLCN on complete data. A frozen prognostic model calibrating across independently ascertained cohorts without updating. (Section 8.5.)


### Section 8.6 Implications

**Q257. What changes tomorrow?**

One recording practice: record the first untreated LDL-C with date and method as a fixed field and grade on it; where none exists, say so. Everything else is better understanding, more accurate communication of scope, or service evaluation to be done. No change in care has been tested. (Section 8.6.)

**Q258. Should risk scores ever reduce FH treatment?**

No. Lipid lowering is already indicated. UK guidance advises against general-population tools in FH, and the 2026 ACC/AHA guideline gives FH-specific scores only COR 2b for short-term risk. Until a model passes calibration, competing-risk, utility, implementation and impact tests, these models support research and benchmarking. (Section 8.6; refs 2, 73.)

**Q259. What could a within-FH model legitimately inform?**

The thesis says intensity, timing in young relatives and allocation of add-on therapy or imaging. The kit asks me to check exact guideline support for each example and to state plainly that a low developmental score cannot delay recommended care. I accept both. (Section 7.1; K012-024.)

**Q260. Minimum data set?**

Untreated LDL-C, lipoprotein(a), body-mass index, diabetes and hypertension status, and dated events. Without them, SAFEHEART-RE was computable for 26 standardised registry patients. (Section 8.6.)


### Section 8.7 Limitations

**Q261. Rank the top five.**

One: no untouched evaluation cohort, both settings in the UK. Two: few events and informative pedigrees, a registry fit on 28 events. Three: treatment reconstruction as a model-based surrogate with non-random reconstructed rows. Four: the grading result estimated only in patients with a documented untreated value. Five: registry event times from recorded ages with heavy censoring. (Table 8.3.)

**Q262. You write that reconstructed rows are not missing at random and that missing diabetes makes MAR doubtful. Is that the right inference?**

Not as written. Observing that missingness is associated with the outcome shows selection; it does not by itself establish MNAR, and it is not a proof that imputation would have failed. The defensible statement is: missingness was outcome-associated, I chose not to impute, and the complete-case set is selected. I will correct the conditional inference in Sections 5.2, 7.4 and 8.7. (K012-029.)

**Q263. Which limitations are resolved by objects that exist?**

L13 is described as resolved by the released two-term lock, and L05 and L06 as superseded by L15. Table B.1 says no two-term frozen object was located, so the lock's existence and usability must be verified before L13 is called resolved. (Section 8.7; Appendix 1, item P3; K012-054.)


### Section 8.8 Resolving studies

**Q264. The principal study and its size.**

Frozen validation of CALON-5 and the two-term model in an independently assembled genotype-confirmed cohort, for example an EAS FHSC registry, with harmonised MACE components and death as a competing event; about 362 ten-year events for calibration-in-the-large within plus or minus 0.1, 248 to 1,285 for the slope. Those numbers are extrapolations from observed differences and should be labelled illustrative; a future design should be built afresh with target incidence, censoring and desired precision. (Section 8.8; K012-044.)

**Q265. The four secondary studies.**

SAIL linkage of the registry to primary-care lipids and prescribing, to measure under-grading service-wide. A cohort recording route at eligibility with exact-allele matching and family-aware estimation. Silent prospective evaluation of frozen TUDOR against fully scored DLCN and the local criteria. And a prospectively locked apoB-residual evaluation in adjudicated FH, head to head with risk-weighted apoB. (Section 8.8, Table 8.4; refs 77, 90, DOIs verified.)

**Q266. Three rules for all of them?**

Ascertainment recorded at eligibility as a design variable. Frozen performance characterised before any updating, each updated model with a new identifier. And no biomarker, score or model adopted because it raises an AUC. (Section 8.8.)


### Section 8.9 Conclusion

**Q267. Your closing sentence, if asked for one.**

LDL-C remains clinically indispensable. What this thesis asks is that it never travel alone: every inference drawn from it, or from a model built on it, should carry the conditions that give the number its meaning. (Section 8.9.) ---


## Appendices A to G

**Q268. What is the admission rule?**

A numerical or methodological statement carries weight only when traceable to a publication, manuscript, aggregate package, producing script or explicit author decision; otherwise it is held or omitted. Tiers T1 to T6; retained aggregates support only sensitivity or descriptive statements. (Appendix A.1.)

**Q269. Which Appendix A rows are T5 or T6?**

T5: the competing-risk sensitivity; the CALON-F grey-zone apoB test; the 71-variant synthesis; the CALON-5 increment over age and sex; the comparator tally; the decision curves; the CALON-C lineage rows. T6: TUDOR family-aware optimism, deleted because unresolvable. (Table A.2.)

**Q270. Why is Appendix B so long?**

Because three models and four manuscript versions share parts of the programme, and no value may be read against the wrong object. It records CALON-C's specification, internal validation, endpoint-date sensitivity, competing-death reconciliation, Welsh transport, comparator benchmarking, reconstruction convention and risk of bias; the Paper 14 version history; and the Welsh denominator reconciliation. (Appendix B.)

**Q271. What did the CALON-C manuscript claim that the thesis withholds?**

Holm-adjusted full-follow-up comparator wins: plus 0.070 against SAFEHEART-RE, plus 0.032 against Montreal, plus 0.015 against the FH-Risk-Score described as a tie; and a post hoc non-inferiority margin. The horizon-specific thesis analyses supersede them, and the margin analysis is held. (Appendix B.6.)

**Q272. The Paper 14 version history: read it to me.**

As submitted: primary plus 1.11, within-pedigree plus 0.34, interval 0.01 to 0.67. Intermediate: plus 1.23; plus 0.34, interval minus 0.07 to 0.76. R3 v5: plus 1.10; plus 0.40, interval minus 0.003 to 0.81, crossing zero. R4: plus 1.09; plus 0.39, bootstrap 0.04 to 0.74, dummy-variable correction minus 0.06 to 0.85. Had R3 stood, Objective 2's falsifier would have been met at that rung. (Appendix B.10.)

**Q273. What is the denominator ontology?**

A classification that keeps data-resource totals, chapter analysis sets, complete-case subsets and comparator-evaluable subsets as different classes, so that a comparator subset is never presented as a cohort total. (Appendix C.)

**Q274. Summarise the reconstruction agreement evidence.**

Convention A, 133 pairs: at a residual fraction of 0.65, bias plus 0.09, limits of agreement minus 3.00 to plus 3.19, MAE 1.24, Lin's concordance 0.593; at 0.50, bias plus 1.79. Person-specific fractions had median 0.69 with a CV of about 25 per cent. TUDOR's 649-pair resource: MAE 1.20, correlation 0.32. Convention B: a scalar sweep 0.65 to 0.75 moved Harrell's C by 0.0022 and 0.0019. Convention C: MAE 1.25, R-squared 0.0975. None validates an individual value. (Appendix D.)

**Q275. The claim that error worsened in the high-risk tail?**

Retired. Against true untreated LDL-C the slopes were negative; the earlier positive slope had placed the reconstruction itself on the horizontal axis. (Appendix D.1.)

**Q276. Which reporting standards, and what does addressed mean?**

STROBE and RECORD for Chapters 4 and 5; TRIPOD+AI and PROBAST+AI for 6 and 7; TRIPOD-Cluster for family-structured data. Addressed means information or an explicit limitation is present; design-level bias remains. (Appendix E; refs 30, 31, 36, 37, 49, DOIs verified.)

**Q277. What would complete the search?**

Embase through Cardiff University Library; a refresh of Elicit, Consensus and Scite after 1 October 2026 with the logged strings verbatim; the eight unread priority full texts; and independent duplicate screening, after which each novelty sentence is re-stated. The kit adds that historical access statuses must not be overwritten and that Crossref verification is not full-text appraisal. (Appendix F.9; K012-055.)

**Q278. Your glossary defines calibration as O/E and slope. Is that complete?**

It is the summary used, not the concept. Calibration also includes flexible curves and local agreement, and the glossary should separate model state, development independence, O/E versus broader calibration, and publication state versus risk of bias. (Appendix G; K012-056.) ---


## Cross-cutting questions


### Integrity

**Q279. If you had to disclose one thing to the examiners today that they might otherwise miss, what is it?**

That the published description of TUDOR includes an ascertainment term the evaluated object does not contain, that no corrigendum has yet been requested, and that I will request one before final submission. (Box 6.1.)

**Q280. Your own abstracts describe dual external validation. Your thesis says cross-cohort evaluation. Which is right?**

The thesis. The abstracts were written under a looser wording rule than the thesis applies, and they will be corrected to match, not the other way round. (Section 6.2; ref 82.)

**Q281. How much of this thesis was written by a language model?**

Drafting, editing and consistency checking were assisted; every passage was reviewed by me; no analysis was run by a model and no model saw participant data. Title screening in the evidence search used a model first pass, single-reviewer, not duplicated, and that is listed as a limit with its fix. (AI statement; Appendix F.8.)

**Q282. Is the published TUDOR AUC 0.760 or 0.7585?**

The published paper reports 0.760 in 1,274 participants with 311 carriers, and the thesis uses 0.760. An earlier internal lineage carried 0.7585; it is not the published figure and I do not quote it. (Section 6.3; ref 21, DOI verified.)


### Causality and design

**Q283. Is anything in this thesis causal?**

No. Chapter 4 is association at fixed LDL-C; Chapter 5 is selection into observation; Chapters 6 and 7 are model tasks. The only causal premise, that LDL causes ASCVD, is cited to the EAS consensus and the guidelines. (Section 3.2; ref 1, DOI verified 10.1093/eurheartj/ehx144.)

**Q284. Name the collider in Chapter 5.**

Inclusion in the analysed sample, conditioned on referral. Prior disease can prompt referral, so proband outcome analyses are collider-prone, which is why the event analyses are secondary and why probands with prior ASCVD were excluded in a sensitivity: plus 0.96. (Sections 5.1, 5.4, 5.5.)

**Q285. What is the target trial Chapter 7 would emulate, and why can't it?**

Eligibility: adults with a P/LP variant, free of ASCVD. Strategy: none, it is prognostic, not interventional. Time zero: baseline. Outcome: first MACE with competing death. Estimand: ten-year cumulative incidence. It cannot emulate a trial because there is no treatment contrast, and it says so: the models rank, they do not estimate effects. (Section 7.1.)

**Q286. Where is immortal time a risk?**

In the registry, where event times derive from recorded ages and baseline is a dated first visit; participants recorded as treated without a start date were excluded from the paired first-versus-later analysis for that reason, and predictors after the outcome are barred by the temporal rule in Section 3.4.3.


### Translation

**Q287. If a lipid clinic asked what to do with this thesis on Monday, what would you say?**

Record the first untreated LDL-C with its date and method, and grade on it. Do not adopt the 0.31 ratio. Do not run TUDOR or CALON-5 on patients. Measure apoB where the question is particle burden, which guidance already supports in treated adults. And collect the minimum data set so that any future FH score can be computed. (Section 8.6.)

**Q288. A cascade-detected 30-year-old with LDL-C 3.8 and a confirmed familial variant asks if she really has FH.**

Yes. Her value is the phenotype observed through a different route, at a younger age, and 84.5 per cent of untreated P/LP carriers in a population cohort were below 4.9. It does not overturn the molecular diagnosis or imply low cumulative burden. Earlier detection is the opportunity to prevent cholesterol-years. (Sections 5.6, 5.7.)


### The candidate

**Q289. What would you do differently?**

Register the framework and falsifiers prospectively. Fix one reconstruction convention with uncertainty propagation. Record route at eligibility. Set aside a validation cohort before any model choice. Request the corrigendum the day the discrepancy was found. Complete the Embase search.

**Q290. What are you least confident about?**

The within-pedigree result, because it depends on the variance convention and has moved across versions; and the CALON-5 increment over age and sex, because there is no interval.

**Q291. What are you most confident about?**

The grading result in Section 7.7, because it needs no model, no reconstruction and no UK Biobank data: 145 of 193 severe carriers were not severe on the value in the record. ---


## Additional hard questions (from the annual-review bank)

**Q292. Is any of this causal?**

No chapter estimates a causal effect (Section 3.2). Chapter 4 is an association at fixed LDL-C; Chapter 5 is selection into observation; Chapter 6 and 7 are model tasks. The verbs are held to that: 'associated with', 'ranked', 'carried information'. The one place where causal knowledge is invoked, that LDL causes ASCVD, is cited to the EAS consensus and the guidelines, not to my data.

**Q293. If a sharp reviewer had one objection per chapter, what would it be, and what disarms it?**

Chapter 4: residual confounding by metabolic state; disarmed only partly by the adjustment ladder and the remnant decomposition, fully only by the resolving study. Chapter 5: polygenic enrichment in probands; not disarmed, disclosed. Chapter 6: comparator starvation; disarmed by naming it and by Akyea 2026 as counterweight, resolved only by a fully scored DLCN comparison. Chapter 7: leakage of receiving-cohort performance into term selection and 28 events; disarmed only by an untouched cohort with about 362 events.

**Q294. What would you do differently if you started again?**

Prospectively register the framework and falsifiers; fix one treatment-reconstruction convention with uncertainty propagation; record ascertainment route at eligibility; separate a validation cohort before any model choice; request the TUDOR corrigendum the day the object discrepancy was found; and complete the Embase search.

**Q295. Why should an MD be awarded for a thesis whose every model is 'not ready for clinical use'?**

Because the thesis's contribution is measurement of the conditions under which a familiar number misleads, with the direction and size of each error in linked genotyped cohorts, and the named repair for each (Section 8.3). Negative and bounded results are the evidence: the retired ratio, the imprecise within-family rung, TUDOR's frozen fall, CALON-5's calibration failure. Declaring a model ready on this evidence would have been the failure.

**Q296. Sample overlap: are Chapters 4 and 5 independent?**

No, and the number of shared Welsh participants was not established (Section 5.7, limitation L03). The UK Biobank frames of Chapters 4, 5 and 7 come from one application and their overlap has not been established either. I treat them as one programme, not as replications.

**Q297. Spin: apply the SPIN-PM framework to your own summary.**

The summary reports the frozen TUDOR fall to 0.669, the calibration failure of CALON-5 and 'no model is ready for clinical use' in the same paragraph as the favourable findings; comparator wins are reported as unresolved; refits are labelled apparent. I cite SPIN-PM and Van Calster 2023 as the standards (references 99, 100). If the committee finds spin, I would like it pointed to. ---

**Q298. What is the operational LDLR-carrier flag, and why does it give 1 in 142 when FH prevalence is nearer 1 in 250 to 300?**

It is a broad variant flag in the locked extract that selected 3,540 of 501,936 participants, and Appendix B.2 says explicitly that this frequency exceeds published estimates of pathogenic FH-variant prevalence and lacked identifiers for complete re-adjudication of pathogenicity, zygosity or functional class. That is why it is never labelled HeFH and why Chapter 7 restricts to ClinVar P/LP single-variant carriers, 1,009 of the 3,209 flagged participants free of ASCVD (Section 7.4). Among the 2,398 untreated flagged carriers in Chapter 5, 29.8% were P/LP and 46.3% VUS (Section 5.6).

**Q299. ClinVar P/LP by database, not expert adjudication. How different are the two cohorts' calls?**

Materially. At two or more review stars, 69.0% of the 1,009 biobank P/LP variants qualified against 97.7% of 954 registry variants; zygosity is recorded in neither cohort (Section 7.4, limitation L09). The P/LP fraction was 31.4% in the biobank and 57.7% in the registry, which is itself an ascertainment gradient in the genetic label. Uniform adjudication is listed as what would remove the limitation (Table 8.3).

**Q300. Why combine LDLR and APOB variants in Chapter 4 when APOB p.(Arg3527Gln) typically gives a milder phenotype?**

The Chapter 4 frame is an operational rare LDLR/APOB variant rule inherited from the source-locked package and its results are frame-specific (Section 4.2). The thesis does not report gene-stratified estimates for Chapter 4 and, at 75 events, could not do so with precision. I accept that APOB and LDLR carriers may differ in particle composition and that a gene-stratified sensitivity should be listed among the resolving-study requirements. Chapter 7 tested LDLR versus other gene as an interaction with the lipid term and the registry cell had too few events (Section 7.13).

**Q301. What does residue matching actually control, and what does variant identity explain?**

Position on the receptor, not the substitution or its function (Section 5.2). Variant identity explains a small share of LDL-C variance: intraclass correlation about 0.145, leave-one-variant-out 0.127 to 0.149 in the TUDOR analyses, and carriers of the same LDLR variants had LDL-C 1.36 mmol/L lower in UK Biobank than in the Welsh clinic (minus 2.03 to minus 0.76). A predicted measure of variant function tracked the untreated LDL-C stratum in relatives, rho 0.26, but not in probands, 0.02, a post hoc pattern consistent with selection acting on the phenotype (Section 5.2).

**Q302. Penetrance: does a cascade-detected relative with LDL-C 3.8 mmol/L have FH?**

Yes, in the molecular sense, and the thesis's interpretive rule is that a lower LDL-C observed through a different route should not overturn an expert-adjudicated familial LDLR diagnosis or imply low cumulative burden (Section 5.7). Expressivity, recognition and prognosis are separated deliberately (Section 5.1). Management follows molecular and clinical context and guidance, not a manufactured proband value.

**Q303. Polygenic risk: what did you find and what does it mean?**

Association without incremental ranking. Chapter 4's exploratory coronary PRS: 1.43 per SD (1.12 to 1.84) with concordance 0.569 and 64.1% missing. Chapter 7 all-carrier: 1.37 (1.12 to 1.68) with a concordance gain of +0.012 (minus 0.000 to +0.040). Consistent with Paquette 2026's FH result of hazard ratio 1.77 and a C change from 0.746 to 0.750 (Section 2.9). Polygenic enrichment among probands remains the unexcluded competing explanation in Chapter 5.

**Q304. Ancestry?**

Section 5.7 states an ancestry restriction limits transportability of the UK Biobank comparison. I need to state the restriction explicitly and confirm it was the same across Chapters 4 and 7 (Part 1, P9).

**Q305. Incidental and direct-to-consumer findings?**

Outside the data; considered in Section 8.6. Their phenotype may resemble genotype-first carriers, raw direct-to-consumer genotypes can be false positive on confirmation (reference 89), and a score developed in a referred population should not be assumed to rank or calibrate for them. A prospective audit is deferred work (Section 8.8).

**Q306. Homozygous FH and PCSK9?**

Homozygous FH is outside the empirical architecture (Section 1.7). Chapter 7's gene set is described as dominant FH genes and I should enumerate it with per-gene counts (Part 1, P8). ---


## Part 7. Biostatistics rapid-fire

**Q307. Events per parameter for CALON-5?**

14.8 in UK Biobank, 5.6 in the registry P/LP fit; floor was 10 (Section 7.3).

**Q308. Optimism?**

Whole-procedure bootstrap: apparent 0.6972, optimism 0.0209, corrected 0.6763; slope corrected 0.901 (Section 7.3).

**Q309. Proportional hazards?**

0 of 10 tests at p under 0.05 for CALON-5, little power on 52 registry events; Chapter 4's diagnostic object not recovered (Sections 7.3, 4.5).

**Q310. Why cause-specific Cox rather than Fine-Gray for the association?**

The estimand was the rate among the event-free; Fine-Gray answers the cumulative-incidence question and is reported only as a sensitivity (Table 3.3).

**Q311. Why Aalen-Johansen?**

One minus Kaplan-Meier overestimates absolute risk when competing death is censored (Section 3.5).

**Q312. Why no imputation in Chapter 7?**

Missingness associated with the outcome, hazard ratio 0.40, so MAR doubtful (Section 7.4).

**Q313. Why 20 imputations in the discordance paper?**

Source study choice for 40% apoB missingness; imputed apoB is not an assay (Section 4.4).

**Q314. What determines the width of a concordance interval?**

Number of events and, for paired comparisons, correlation between scores, not cohort size (Section 3.8).

**Q315. Twenty-event floor?**

Comparisons below 20 events at a declared horizon are not interpreted (Section 3.8).

**Q316. Sequential interval rule for updating?**

Baseline and slope where the slope interval excludes one; baseline alone where only the O/E interval does; otherwise none; not a closed testing procedure (Section 3.8).

**Q317. Why no equivalence claims?**

No prespecified margin; an interval including zero is 'undemonstrated, not demonstrated absent' (Section 3.7.1, reference 80).

**Q318. NRI weaknesses?**

Sensitive to prevalence, category definitions and miscalibration; secondary only (Sections 6.3, 6.7, reference 88).

**Q319. Brier score?**

Overall error at a horizon combining discrimination and calibration; not a calibration metric; small when events are rare (Section 3.8).

**Q320. Wild cluster bootstrap?**

Rademacher, 1,999 pedigree resamples, same interval 0.83 to 1.36 as the clustered SE (Section 5.3).

**Q321. Why is the dummy-variable df correction relevant?**

With 488 pedigree fixed effects it counts each as a parameter and widens the within-pedigree interval to minus 0.06 to 0.85 (Section 5.3).

**Q322. Why is the +0.052 Montreal win fragile?**

Lower limit +0.0003 on 20 events, and it did not resolve under the earlier hypertension definition (Section 7.8).

**Q323. Why did you change the hypertension definition?**

To recorded medication only, on 24 September 2026, before any refit result was seen; the earlier definition is a sensitivity (Section 7.2).

**Q324. Time-dependent AUC versus Harrell C?**

Both reported: 0.718 and 0.713 at ten years alongside 0.699 and 0.698 (Section 7.5).

**Q325. Why is target refit AUC 0.756 not comparable with frozen 0.669?**

Different model state; outcomes informed the coefficients (Section 6.5).

**Q326. What is 'calculability'?**

Rows on which an instrument can be computed as published; SAFEHEART-RE 26 registry rows with one event (Table 7.3).


## Part 8. Lipid-medicine rapid-fire

**Q327. Why does LDL-C under-grade after treatment?**

Median 6.1 untreated versus 3.6 on treatment in 261 registry pairs; 66.3% drop a stratum (Section 7.7).

**Q328. Why not just back-calculate?**

Because every convention has MAE around 1.2 mmol/L and correlation 0.09 to 0.32 with documented values (Part 1, P11).

**Q329. Where is Friedewald least reliable?**

Below 1.4 mmol/L the triglyceride term is 47% of the calculated value (Section 7.7).

**Q330. Non-HDL-C alternative in Chapter 5?**

+1.67 mmol/L crude, +1.44 age- and sex-adjusted; same-visit differences 1.64 and 1.49 (Section 5.4).

**Q331. Does a severity cut-point rank risk?**

Poorly: concordance 0.552 to 0.593 for 4.9 mmol/L alone; but the event-rate ratio above versus below was 2.97 in the registry and 1.61 in the biobank (Sections 7.5, 7.7).

**Q332. Reaching goal means low lifetime exposure?**

No: of 51 registry patients at the 2.6 goal with a recorded pre-treatment value, 54.9% had been at or above 4.9 (Section 7.7).

**Q333. Lp(a)?**

Modifier, not subtractable component; measured in only 143 of 1,904 Welsh carriers in Chapter 5; Lp(a) of at least 125 nmol/L over CALON-5 gave 0.88 (0.40 to 1.94) (Sections 5.7, 7.9).

**Q334. Should clinics measure apoB?**

Where the question concerns particle burden; guidance already supports it in treated adults (COR 2a); the 0.31 ratio threshold should not be adopted (Section 8.6).

**Q335. Which ratio direction did the paper use?**

apoB/LDL-C in g/mmol, threshold 0.31; other literature uses LDL-C/apoB as a size proxy (Section 2.5).

**Q336. Diabetes and TUDOR?**

AUC 0.648 Wales, 0.642 UK Biobank: a metabolically challenging setting, not a subgroup licence (Section 6.5).

**Q337. Cascade versus proband age?**

Median 60.2 versus 38.6 years; treated 89.0% versus 67.5%; prior ASCVD 21.3% versus 6.9% (Table 5.1).

**Q338. Minimum data set a service needs for any FH risk score?**

Untreated LDL-C, Lp(a), BMI, diabetes and hypertension status, dated events (Section 8.6).

**Q339. Menopause and FH LDL-C?**

LDLR carriers +11.5% versus +15.3% percentage rise from a higher base; Da Roza 2026 found a larger absolute rise in monogenic FH; compatible (Section 8.4).

**Q340. Pregnancy?**

Klevmoen 2021 documented loss of statin-treatment years; adds untreated exposure current status cannot summarise (Section 2.11).


## Annual-review questions: progress, risks and plan

**Q341. What is the publication status of each component?**

ApoB/LDL-C discordance paper: Published J Clin Lipidol 2026 with erratum (highlights only) (reference 51, 72). TUDOR: Published online J Clin Lipidol 2026; description-object discrepancy uncorrected (reference 21). TUDOR congress abstract: Atheroscler Plus 2026; related work, not independent (reference 82). Ascertainment (Paper 14): Under revision, R4, J Clin Lipidol (reference 52). CALON-5 cross-cohort evaluation: Submitted manuscript, revision r7 (reference 83). CALON-C: Developmental manuscript, lineage only (reference 22). Menopausal transition study: In press (reference 50).

**Q342. What is on your critical path to submission?**

Recommended answer, ordered by dependency and impact (executive-planner view): **Correct the internal inconsistencies in Part 1** (P1 to P6, P12): one week; no data access needed. **TUDOR corrigendum request** to J Clin Lipidol for the Index Effect descriptor: draft this month; integrity item that examiners weigh most heavily. **Registry-side analyses that need no UK Biobank access:** family-grouped bootstrap of the Welsh TUDOR paired differences (closes L11); interval for CALON-5 over age and sex on registry rows (closes half of L14); enumerate Chapter 7 genes and ancestry restriction (P8, P9); complete TUDOR PROBAST+AI (P12). **UK Biobank RAP:** either restore access under Application 1002450 or record formally why not; if restored, re-verify the retained aggregates flagged in Table A.2 and run the frozen reconstruction-uncertainty protocol (Appendix D.4). **Paper 14 R4** through review; if the within-pedigree rung changes again, Chapter 5 and Table B.2 change with it. **CALON-5 r7** through review; align Chapter 7 and Appendix A to the accepted version. **Complete the search:** Embase via Cardiff University Library, refresh Elicit, Consensus and Scite after 1 October 2026, read the eight unread priority full texts, duplicate screening (F.9). Only then: final formatting, word count against Cardiff regulations, and viva preparation.

**Q343. What could stop you submitting on time?**

Three things: the RAP closure if any examiner requires re-verification of a UK Biobank primary value; a reversal of the Paper 14 within-pedigree result on review; and the TUDOR corrigendum process if the journal asks for more than a wording correction. None affects the Chapter 7 grading result, which needs no reconstruction and no UK Biobank data.

**Q344. Training needs?**

Two the thesis exposes: formal multiple-imputation and MNAR sensitivity methods for the registry (Chapter 7 avoided imputation rather than modelling it); and systematic-search methodology so that the novelty statements can be made unconditional.

**Q345. The single next action after this meeting?**

Correct P1 to P6 and draft the corrigendum letter; both are within my control and neither needs data access. ---


## Questions arising from the Kit 012 paragraph audit

**Q346. Where is the crosswalk from each condition to population, lipid quantity, estimand and clinical use, with untested cells shown? (Sections Summary, 1.5, 8.3)**

Concede it is implicit in Tables 1.1 and 3.1; commit to an explicit crosswalk with untested cells. Audit code K012-001.

**Q347. Why does the Summary drop the interval and call them "carriers"? (Sections Summary)**

Restore 1.070 to 1.718 and keep "rare LDLR/APOB variant" wording. Audit code K012-002.

**Q348. Which interval convention is preferred for 0.39, and what is the +0.47 family-model contrast the kit mentions? (Sections Summary, 5.3, 8.2.2)**

State the preferred estimate and both conventions; separate the route association from any within-family model contrast (5.3c). Audit code K012-003.

**Q349. Is TUDOR's whole gain due to reconstruction? (Sections Summary)**

No; architecture-level performance, no ablation (6.4d; 8.2.3). Audit code K012-004.

**Q350. Baseline hazard or marker effect? (Sections Summary, 8.4)**

Baseline hazard; not an FH result (8.4b, 8.4c). Audit code K012-005, 052, 053.

**Q351. Which of the five "additions" are demonstrated, which are tools, which are proposals? (Sections Summary, 8.9)**

Demonstrated interpretation: particle information, route contrast, under-grading. Partially supported tool: TUDOR, local recalibration. Proposed: sex-specific calibration in FH. Audit code K012-006.

**Q352. Is 4.9 an FH-severity grade, and was the assay "accurate"? (Sections 1.1, 7.7, 8.3)**

Reword category; remove accuracy assertion (8.3c, 8.3d). Audit code K012-007, 008.

**Q353. Universal interchangeability or limits of interpretation? (Sections 1.1, 1.5)**

The aim is the limits; reword "establish" to "evaluate". Audit code K012-009.

**Q354. Is the particle mechanism measured or assumed? (Sections 1.2, 8.2.1)**

Plausible background, not measured mediation (4.7a, 4.7b). Audit code K012-010.

**Q355. Does the rare-variant frame answer the HeFH question? (Sections 1.2, 1.5)**

Part of it; no population upgrade at the objective (4.2e, 4.9). Audit code K012-011.

**Q356. Are falsifiers protection or stress tests? (Sections 1.5)**

Interpretive stress tests written retrospectively (1.5a). Audit code K012-012.

**Q357. Common definitions or linked heterogeneous studies? (Sections 1.6, 8.3)**

The latter; build the overlap matrix (8.3e). Audit code K012-013.

**Q358. Was the joint-novelty claim tested adversarially for combinations? (Sections 1.6, F.10)**

Only in PubMed and Crossref; narrow to the supported combination (2.1a). Audit code K012-014.

**Q359. Is governance a method? (Sections 1.6)**

Good practice, not novelty (8.5b). Audit code K012-015.

**Q360. Define each genetic label positively. (Sections 1.7)**

Align with Box 3.1 (1.7b). Audit code K012-016.

**Q361. Ratios, components and residuals compared under explicit assumptions? (Sections 2.5)**

Retain the qualification; compare all three where data allow (A04). Audit code K012-017.

**Q362. Is Akyea 2026 the "only" specialist comparison? (Sections 2.8)**

Only one identified; search further before "only" (6.4b). Audit code K012-018.

**Q363. Lead with magnitude and precision. (Sections 5.3)**

Done in 5.3a to 5.3d. Audit code K012-019.

**Q364. Is the non-HDL excess a remnant mechanism? (Sections 5.4)**

Timing explains most; mechanism not shown (5.4b). Audit code K012-020.

**Q365. THESIS-AN-02's output was not located. Why is the number in the thesis? (Sections 5.4)**

Demote or hold; do not substitute the 53-variant rerun (5.4c). Audit code K012-021.

**Q366. Is the first-versus-later reading a mechanism? (Sections 5.4)**

A paired phenotype observation, not a mechanism (5.4e). Audit code K012-022.

**Q367. Is the object locked with metadata and readback? (Sections 6.2)**

Yes within tolerance; no Index Effect added from text (6.2a, 6.2b). Audit code K012-023.

**Q368. Guideline support for intensity, timing, add-on examples? (Sections 7.1)**

Check each; state a low score cannot delay care (8.6c). Audit code K012-024.

**Q369. Is r7 "submitted" evidenced? (Sections 7.1)**

Use "manuscript revision r7" until receipt is evidenced. Audit code K012-025.

**Q370. Blank medication field read as "no" for 436 of 945. Justified? (Sections 7.2)**

Recording convention with dated lipid visit; provenance retained (7.2a). Audit code K012-026.

**Q371. Was the reconstruction cross-validation person-grouped and carrier-trained? (Sections 7.2)**

Verify; do not conflate CV with external agreement (7.2c). Audit code K012-027.

**Q372. What was nested in the bootstrap? (Sections 7.3)**

The five-term refit; not the term choice; hold independent-validation language (7.3b, 7.3c). Audit code K012-028.

**Q373. Is outcome-associated missingness proof of MNAR? (Sections 5.2, 7.4, 8.7)**

No; correct the inference (8.7b). Audit code K012-029.

**Q374. Endpoint rounding, ties, censoring positivity? (Sections 7.4)**

Qualify unstable ten-year estimates (7.4c). Audit code K012-030.

**Q375. "Every interval overlapping" is not equivalence. (Sections 7.5)**

Agreed; report range and uncertainty only (7.5a). Audit code K012-031.

**Q376. Is beating the 4.9 cut-point a fair comparison? (Sections 7.5)**

It tests a misapplication, not a legitimate model; the fair baseline is age and sex (7.5b, 7.5c). Audit code K012-032.

**Q377. Was calibration "restored" by a ratio rule? (Sections 7.6)**

Level fixed on held-out rows; spread and local calibration separate (7.6b). Audit code K012-033.

**Q378. Did relative effects travel? (Sections 7.6, 8.2.4)**

Partially shown; not formally tested (8.2c). Audit code K012-034.

**Q379. Does any reconstruction restore individual severity? (Sections 7.7, 8.2.3)**

No; recorded pretreatment retrieval is the simpler remedy (7.7b). Audit code K012-035, 050.

**Q380. Goal eligibility, timing, mathematical coupling? (Sections 7.7)**

Verify denominators; discuss coupling and room for reduction (7.7d). Audit code K012-036.

**Q381. "Less penetrant" by what definition? (Sections 7.7)**

Descriptive; follow-up, treatment, age and endpoint differ (7.7f). Audit code K012-037.

**Q382. Is the TG-term share an error metric? (Sections 7.7)**

No; needs reference-method comparison (7.7e). Audit code K012-038.

**Q383. Comparator equations, units, exact lower limit; naive pooling? (Sections 7.8)**

Retain +0.0003 and the hypertension sensitivity; pooling exploratory (7.8b). Audit code K012-039, 040.

**Q384. Are the Australian validations being minimised? (Sections 7.8)**

State their layers exactly (7.8d). Audit code K012-041.

**Q385. Different marker definitions across chapters? (Sections 7.9)**

Ratio here, residual there; moderate effects not excluded (7.9a). Audit code K012-042.

**Q386. Is the BMI bottleneck verified, and is research imputation forbidden? (Sections 7.11)**

Verify; not forbidden, deliberately not done (7.11a, 7.11b). Audit code K012-043.

**Q387. What does 362 assume? (Sections 7.11, 8.2.4, 8.8)**

Observed-effect extrapolation; label illustrative (8.8a). Audit code K012-044.

**Q388. Undefined DCA cells: estimator or mathematics? (Sections 7.12)**

Inspect; preserve as missing (7.12a). Audit code K012-045.

**Q389. Reclassification categories and action equivalence? (Sections 7.12)**

Exploratory NRI; no causal harm claim (7.12c). Audit code K012-046.

**Q390. Expected nominal positives from the actual test count? (Sections 7.13)**

Compute from tests and alpha (7.13a). Audit code K012-047.

**Q391. Sex-specific calibration validated for FH? (Sections 7.15, 8.3)**

No (8.4d). Audit code K012-048.

**Q392. Connect the seven steps into argument. (Sections 8.2)**

Agreed (8.2a). Audit code K012-049.

**Q393. Errors "predictable in direction"? (Sections 8.3)**

Group-level only (8.3b). Audit code K012-051.

**Q394. Is the two-term lock located and usable? (Sections 8.7)**

Verify against Table B.1 (8.7c). Audit code K012-054.

**Q395. Five-platform session logged? (Sections F.10)**

Not yet; do not overwrite historical statuses (F1). Audit code K012-055.

**Q396. Separate the four concepts. (Sections Glossary)**

Agreed (G1). Audit code K012-056.


## General viva questions

**Q397. Why did you choose this topic, and why an MD rather than a PhD?**

The topic came from clinical practice: a treated patient whose recorded LDL-C would have graded them as mild when their untreated value was severe. The MD suited work that is clinical-epidemiological and embedded in a service, using registry and biobank data rather than laboratory science.

**Q398. What was the original plan and how did it change?**

The programme began with the discordance paper and the CALON-C model. CALON-5 replaced CALON-C when the exact CALON-C object could not be recovered. The framework and falsifiers were written across the completed programme, and the thesis says so. Key dates: the discordance paper was submitted in August 2025 and published in March 2026; the residual reanalysis was frozen on 12 July 2026; the two-term selection protocol is dated 23 September 2026 and the hypertension definition changed on 24 September 2026.

**Q399. What research training did you complete, and what gaps did the thesis expose?**

Two gaps: multiple imputation with sensitivity analysis for data that may not be missing at random, because Chapter 7 avoided imputation rather than modelling it; and systematic-search methodology, because the novelty statements rest on an incomplete search. I have planned supervision or a course for each.

**Q400. What approvals cover the Welsh data and UK Biobank?**

UK Biobank analyses were under Application 1002450. Welsh analyses use the approvals, consent and governance arrangements recorded in the controlling source manuscripts, and the identifiers and wording are reproduced from the UK governance records in the declarations.

**Q401. Could any table re-identify a participant?**

No. The thesis reports aggregates only: no participant or family identifiers, no exact dates, no row-level predictions and no rare combinations. Small cells were checked before release, and statistical verification and governance approval are separate gates.

**Q402. Who are the co-authors on each incorporated paper and what did each do?**

TUDOR and the ascertainment paper: Genedy, Haralambos, Yousef, Cole, Datta and Zouwail. The discordance paper and the midlife study: Genedy and Zouwail. The CALON manuscripts are sole-authored. I led conception, design, analysis, interpretation and writing on each; co-authors provided data access, clinical adjudication, biochemistry and supervision.

**Q403. Your supervisors are co-authors. How is the thesis your own work?**

The statement of original work says it: scientific conception, study design, data curation, the producing analyses, interpretation and final editorial decisions are mine. Statistical supervision came from Dr Aubin, clinical and lipidological supervision from Professors Cole and Yousef. I take responsibility for every number.

**Q404. What is the status of each paper today?**

The discordance paper is published with an erratum that changed only the highlights. TUDOR is published online. The ascertainment paper is at revision four, under revision. The CALON-5 manuscript is at revision seven. The midlife study is in press. The congress abstracts are related work, not independent evidence.

**Q405. Have any of your outputs been corrected or retracted?**

One erratum, correcting typographical errors in the highlights of the discordance paper; no estimate, table, method or conclusion changed. A corrigendum for the TUDOR paper is to be requested because the published description includes an ascertainment term that the evaluated object does not contain.

**Q406. Which parts of the text did a language model draft?**

Generative AI assisted with literature searching and with drafting, editing and consistency checking throughout. Every passage was reviewed by me. No model ran an analysis or had access to participant-level data. Title screening in the evidence search used a language-model first pass with single-reviewer adjudication, which is listed as a limit.

**Q407. Can the examiners see the code?**

Producing scripts and locks exist for CALON-5 and for the CALON-C high-bar rerun, with hashes. Several inherited analyses are source-locked rather than independently reproduced. UK Biobank scripts cannot be re-run while the Research Analysis Platform is closed to me. I can provide an index of scripts, locks and hashes.

**Q408. Why is the UK Biobank platform closed to you?**

Access ceased during the final phase of the programme. The thesis states the consequence, retained aggregates that cannot be re-verified, but not the reason, and I will add one sentence giving it.

**Q409. What happens if a reviewer of the ascertainment paper changes the within-pedigree result again?**

Chapter 5 and the version history in Appendix B change with it, and the Objective 2 falsifier status in Table 8.1 is reassessed. No other chapter depends on that value.

**Q410. What is the single most original thing in the thesis?**

Measuring five conditions that change what an LDL-C value means, together, in two genotyped cohorts under one evidence rule; and, as the most direct finding, that three in four registry carriers with a severe untreated LDL-C would not have been graded severe on the value in their record.

**Q411. What is the weakest chapter and why?**

Chapter 7, because neither cohort is untouched by model development, the registry fit rests on 28 events, and most of the discrimination is demographic. I say that before the examiner does.

**Q412. If you removed one chapter, which and what would be lost?**

Without Chapter 4 the particle condition would rest on general-population evidence. Without Chapter 5 the route condition would be asserted, not measured. Without Chapter 6 the treatment condition would have no identification test. Without Chapter 7 the separation of ranking from absolute risk would have no FH demonstration.

**Q413. Is this one thesis or four papers stapled together?**

One question, five conditions, one evidence rule and the same two cohorts. The chapters are linked by inferential dependency, each inheriting an open problem from the one before, not by narrative.

**Q414. Why Cox and not a flexible parametric model?**

The estimands were a cause-specific hazard and a ranking, with proportional hazards checked where the record allows. The thesis states that flexible parametric models would be preferable where smooth absolute risk or time-varying effects are primary.

**Q415. Why elastic net for TUDOR?**

The inputs are correlated lipid measures. Penalisation stabilises the fit without making coefficients biological effects, and the thesis does not interpret them as such.

**Q416. Why no penalty for CALON-5?**

Five terms at 14.8 events per parameter in the development cohort, and the record states no penalty. CALON-C used ridge with nine correlated terms; CALON-5 does not need it.

**Q417. Why five and ten years?**

They are the guideline and comparator horizons; ten years is primary and five years is reported alongside.

**Q418. How were ties and event ages handled in the registry?**

Event times derive from recorded ages, which is an approximation, and censoring is heavy. Ties and rounding are stated as limitations and the ten-year registry estimates are qualified accordingly.

**Q419. What would change your mind about the residual?**

A frozen residual failing to reproduce under harmonised measurement in an independent FH cohort.

**Q420. What would change your mind about ascertainment?**

No residue-matched proband to cascade difference in an independently designed cohort with route recorded at eligibility.

**Q421. What would change your mind about TUDOR?**

Frozen TUDOR failing against fully scored DLCN on complete data.

**Q422. What would change your mind about CALON-5?**

A frozen prognostic model calibrating across independently ascertained cohorts without updating.

**Q423. Would you order apoB in every FH patient?**

No. Where the clinical question concerns particle burden not represented by LDL-C. Guidance already supports apoB in treated adults, and that recommendation does not rest on this thesis.

**Q424. Would you use TUDOR in your clinic tomorrow?**

No. It is a research triage score, not calibrated for a new setting, and not yet compared with fully scored DLCN or with the Welsh service criteria.

**Q425. Would you tell a patient their CALON-5 risk?**

No. Its absolute risks miscalibrated in both directions on transfer. It is a research model.

**Q426. What is the harm if someone misuses your findings?**

Withholding or reducing indicated treatment on a low score, or dismissing a cascade relative’s molecular diagnosis because of a lower LDL-C. Both are explicitly prohibited in the thesis.

**Q427. What are you doing next, in order?**

Correct the internal inconsistencies; request the TUDOR corrigendum; run the registry-side analyses that need no UK Biobank access; restore or formally close UK Biobank access; take the ascertainment paper and the CALON-5 manuscript through review; complete the Embase search; then final formatting.

**Q428. How long to final submission?**

I will give a date agreed with my supervisors, with buffers, and name the critical-path item, which is either UK Biobank access or the ascertainment paper’s review.

**Q429. Which resolving study would you fund first?**

Frozen validation of the prognostic models in an independent genotype-confirmed cohort with several hundred ten-year events, because it decides whether any routine-data model can prioritise review in FH services.

**Q430. What did you learn about yourself as a researcher?**

To request corrections the day a discrepancy is found, and to prespecify before results are known.

**Q431. What advice would you give a new MD student in this field?**

Register the framework first; record ascertainment route at eligibility; set aside a validation cohort before choosing any model; keep one treatment-reconstruction convention.

**Q432. Why is the thesis so long, and could it be shorter?**

About 67,700 words, of which the appendices carry provenance, lineage and the search record. The appendices could move to supplementary files if the regulations allow.

**Q433. Why are some numbers given to three or four decimals?**

To match the source manuscripts and to avoid rounding a lower confidence limit through zero, as with the Montreal comparison at plus 0.0003.


## Questions about specific citations

**Q434. You cite Akyea 2026 (reference 45). What exactly does it report, and is your wording accurate?**

The paper reports 885 patients referred for FH genetic testing, 267 with a confirmed variant; DLCN AUROC 0.816 (0.784 to 0.847) against FAMCAT 0.748 (0.712 to 0.784). Verified in PubMed. DOI 10.1016/j.jacl.2026.07.009.

**Q435. You cite Johannesen 2024 (reference 38). What exactly does it report, and is your wording accurate?**

Expected apoB was defined by regressing apoB on LDL-C among people with triglycerides at or below 1 mmol/L, in statin-free participants of the Copenhagen General Population Study. Verified. DOI 10.1016/j.jacc.2024.03.423.

**Q436. You cite Sayed 2024 (reference 15). What exactly does it report, and is your wording accurate?**

12,688 NHANES adults not on statins; at LDL-C 100 mg/dL the 95 per cent range of apoB was 66 to 99 mg/dL. Verified. DOI 10.1001/jamacardio.2024.1310.

**Q437. You cite Marston 2022 (reference 13). What exactly does it report, and is your wording accurate?**

ApoB hazard ratio 1.27 (1.15 to 1.40) per SD for myocardial infarction is exact. The abstract says non-HDL-C and triglycerides were not associated after adjustment for apoB; it does not name LDL-C in that sentence, so I will reword. DOI 10.1001/jamacardio.2021.5083.

**Q438. You cite Johannesen 2021 (reference 14). What exactly does it report, and is your wording accurate?**

Discordantly high apoB with low LDL-C in statin-treated patients: hazard ratio 1.49 (1.15 to 1.92) for myocardial infarction. Verified. DOI 10.1016/j.jacc.2021.01.027.

**Q439. You cite Jin 2023 (reference 61). What exactly does it report, and is your wording accurate?**

89,422 statin-free UK Biobank participants, 3,821 coronary events; cholesterol per LDL particle hazard ratio 1.03 (99 per cent CI 0.94 to 1.14); apoB correlated 0.99 with LDL particle concentration. The adjustment was for VLDL, LDL and HDL particle concentrations together, and I will make that exact. DOI 10.1161/JAHA.123.029552.

**Q440. You cite Morze 2025 (reference 62). What exactly does it report, and is your wording accurate?**

207,368 participants; adding lipoprotein(a) increased the AUC from 0.769 to 0.774. Verified. DOI 10.1093/eurheartj/ehaf207.

**Q441. You cite Paquette 2026 (reference 79). What exactly does it report, and is your wording accurate?**

Polygenic score hazard ratio 1.77 (1.20 to 2.61) in FH; C-statistic 0.746 to 0.750, p 0.60. Verified. DOI 10.1093/eurjpc/zwag203.

**Q442. You cite Tybjaerg-Hansen 2005 (reference 91). What exactly does it report, and is your wording accurate?**

Carriers of the same LDLR mutations had progressively higher cholesterol when identified in the general population, among patients with ischaemic heart disease and in an FH clinic. The gradient reported is total cholesterol, 2.9, 4.1 and 4.9 mmol/L; LDL-C is given as a single difference of about 1.6 mmol/L. My text says LDL-C increments and I will correct it. DOI 10.1161/01.ATV.0000149380.94984.f0.

**Q443. You cite Richardson 2020 (reference 11). What exactly does it report, and is your wording accurate?**

Multivariable Mendelian randomisation: apoB odds ratio 1.92 (1.31 to 2.81) with LDL-C attenuating to 0.85 (0.57 to 1.27). Verified. DOI 10.1371/journal.pmed.1003062.

**Q444. You cite Sniderman 2011 (reference 10). What exactly does it report, and is your wording accurate?**

Relative risk ratios 1.43 for apoB, 1.34 for non-HDL-C and 1.25 for LDL-C. Verified. DOI 10.1161/CIRCOUTCOMES.110.959247.

**Q445. You cite Pan 2026 (reference 39). What exactly does it report, and is your wording accurate?**

Excess apoB in statin-treated coronary disease: adjusted hazard ratios 1.12 (1.06 to 1.18) for all-cause and 1.24 (1.15 to 1.34) for cardiovascular mortality, with external validation in 13,702 UK Biobank participants. Verified. DOI 10.1186/s12944-026-02928-z.

**Q446. You cite McKay 2022 (reference 28). What exactly does it report, and is your wording accurate?**

SAFEHEART-RE in English routine care: Harrell’s c 0.67 (0.61 to 0.72) with substantial miscalibration. Verified. DOI 10.1016/j.atherosclerosis.2022.07.011.

**Q447. You cite Mourre 2025 (reference 40). What exactly does it report, and is your wording accurate?**

Cascade-screened individuals started statins about 14 years earlier and had 51 per cent fewer ASCVD events; the association disappeared after age and sex matching. Verified. DOI 10.1093/eurjpc/zwaf234.

**Q448. You cite Bogsrud 2025 (reference 85). What exactly does it report, and is your wording accurate?**

Newborn LDL-C 1.26 versus 1.18 mmol/L for null versus non-null variants, p 0.552. The contrast is within the 61 FH newborns of a 113-newborn sample; my text says 113 and I will correct it. DOI 10.1093/eurheartj/ehaf815.

**Q449. You cite Trinder 2024 (reference 41). What exactly does it report, and is your wording accurate?**

Clinically diagnosed FH in a Canadian registry (1,123) had higher LDL-C and more ASCVD than genetically identified FH in UK Biobank (723), with polygenic enrichment. Verified. DOI 10.1161/ATVBAHA.123.320287.

**Q450. You cite Weng 2019 (reference 43). What exactly does it report, and is your wording accurate?**

FAMCAT external validation in primary care: positive predictive value 0.84 per cent. Verified. DOI 10.1016/S2468-2667(19)30061-1.

**Q451. You cite Your TUDOR paper (reference 21). What exactly does it report, and is your wording accurate?**

Title: an ascertainment-aware model for LDLR genetic-testing triage. The abstract says the model encodes index-referral versus cascade-relative ascertainment; the phrase Index Effect does not appear in the title or abstract. Head-to-head AUC 0.760 in 1,274 with 311 carriers against eDLCN 0.652, FAMCAT 0.600, Simon Broome 0.569 and MEDPED 0.524. DOI 10.1016/j.jacl.2026.06.030.

**Q452. You cite Your Atherosclerosis Plus abstract (reference 82). What exactly does it report, and is your wording accurate?**

It claims dual external validation in 4,028 genetically confirmed cases. The thesis treats it as related work on the same two cohorts, not independent validation. The record could not be located in PubMed and the 4,028 figure appears nowhere else in the thesis; I will verify it against the abstract PDF before quoting it.

**Q453. You cite Sniderman 2013 INTERHEART (reference 17). What exactly does it report, and is your wording accurate?**

Cited for the warning that apparent superiority of apoB can arise from confounding by related variables. The thesis lists a PMID; the record did not resolve in the tools used for checking, so I quote it as the thesis cites it.

**Q454. You cite Your midlife study (reference 50). What exactly does it report, and is your wording accurate?**

273,036 UK Biobank women; at age 58 an apoB of at least 1.2 g/L corresponded to ten-year MACE risk of 10.88 per cent in men and 4.34 per cent in postmenopausal women; a sex-neutral model over-predicted in women and under-predicted in men. In press, not yet indexed.


## Corrections to make before the meeting (concede these if not yet fixed)

- P1. Decision-curve counts: Section 7.12 says 13 of 30 and 10 of 30; Table A.2 says 12 and 9. Align to r7.
- P2. Comparator tally: Section 7.8 says 3 of 4 and one resolved; Table A.2 says 4 of 4 and none resolved. Table 7.2 supports the chapter. Correct the ledger.
- P3. Two-term lock: Sections 3.4.1, 7.2, Table 7.1 say released; Table B.1 says not located. Reconcile (K012-054).
- P4. Table 8.1 Objective 4 falsifier status reads backwards. Correct; note the ranking falsifier cannot be judged without an interval (L14).
- P5. Reference 52 status: update to R4 under revision.
- P6. Objective 2 falsifier wording differs between Table 1.2 and Section 5.3.
- P7. Chapter 4: state the prevalent-ASCVD exclusion rule and count for the 1,461 frame.
- P8. Chapter 7: enumerate genes and per-gene counts.
- P9. Define the UK Biobank ancestry restriction and confirm consistency across chapters.
- P10. Reference 77 cited for two different things (Sections 2.8 and 8.8); check or split.
- P11. Bring the three reconstruction-agreement figures (r 0.09; MAE 1.20 and r 0.32; MAE 1.25 and R-squared 0.0975) into one table with one sentence on why they differ; reconcile with the 4.5 mmol/L bias statement in Section 7.2.
- P12. Complete TUDOR PROBAST+AI or say why not (Table 6.5 versus Appendix E.4).
- P13. Title leads with particle burden; prepare a defence or a subtitle.
- P14. Tybjaerg-Hansen 2005: total cholesterol, not LDL-C increments (Sections 1.3, 5.7, Table 2.4).
- P15. Jin 2023: adjustment for VLDL, LDL and HDL particle concentrations (Sections 2.3, 4.7).
- P16. Marston 2022: the abstract names non-HDL-C and triglycerides, not LDL-C (Section 2.3).
- P17. Bogsrud 2025: contrast within 61 FH newborns, not 113 (Section 5.2).
- P18. "Index Effect" does not appear in the TUDOR title or abstract; quote the abstract's wording (Section 6.2).
- P19. Reference 82 not indexed; verify title and the 4,028 figure against the abstract PDF.
- K012-007, 008. Reword the 4.9 category and remove "the laboratory value was accurate" (Sections 1.1, 7.7, 8.3).
- K012-013. Replace "common definitions" with an overlap and definition matrix (Sections 1.6, 8.3).
- K012-021. Demote or hold the 71-variant synthesis until its output is located (Section 5.4).
- K012-029. Correct the MAR/MNAR inference (Sections 5.2, 7.4, 8.7, Table 8.3).
- K012-038. Do not present the TG-term share as an error metric (Section 7.7).
- K012-044. Label the 362-event figure illustrative and specify its assumptions (Sections 7.11, 8.2.4, 8.8).
- K012-045. Inspect the undefined decision-curve cells and preserve them as missing (Section 7.12).
- K012-034, 051. Restrict "relative effects travel" and "errors predictable in direction" to what was tested (Sections 7.6, 8.2.4, 8.3).
- K012-025. Use "manuscript revision r7" until submission is evidenced (Section 7.1).
- K012-002, 003. Restore the interval in the Summary; state the preferred within-pedigree convention (Summary, Section 5.3).
- K012-049, 056. Connect the Section 8.2 steps into argument; expand the glossary definitions.

*All numbers are quoted from MD_SUBMISSION_v3_FINAL.docx at the section named. DOIs are given only where they resolved in PubMed on 27 September 2026.*

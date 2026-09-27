# Evidence and novelty search specification

## Scope and dates

Prioritise 2021–27 September 2026, cover 2016–2026 and include named older landmarks when they bound priority or establish a method. Update the end date when the executor actually searches. Search genetic-confirmation-restricted and broader clinical-FH literature in separate strata; broader literature can refute a conceptual novelty claim even when it does not answer the exact empirical question.

All five required platforms must receive a seed query and at least one adversarial query for the relevant claim families: PubMed, Elicit, SciSpace, Consensus and Scite. Use explicit fields/MeSH/Boolean syntax where supported and natural-language questions where required. Do not paste one platform's syntax into another without checking semantics. Metadata, snippets, AI summaries and full text are different evidence levels.

## Current planning pilot

Read CONNECTOR_PILOT.md, BROWSER_RESEARCH_ADDENDUM.md and the saved raw/downloaded records. SciSpace, Scite and a PubMed route returned discovery records. Elicit reported API-plan access denial; Consensus reported exhausted monthly quota. Historical claims that SciSpace is disconnected or Scite unavailable must not be carried forward as current facts. After the user authorised Chrome access, Elicit returned a completed browser response and SciSpace produced search exports. The browser route therefore resolves Elicit access for discovery; the API restriction remains route-specific. Pilot/browser records are UNADMITTED, and this limited work is not a completed novelty review.

On execution, retry only when the route or state can reasonably have changed; inspect authorised browser sessions as a permitted fallback without bypassing access controls. Do not create paid reports, purchase credits or upload the thesis. If blocked, log AUTHENTICATION_REQUIRED, CONNECTOR_UNAVAILABLE or QUERY_FAILED with the actual error, route and time. A quota/API-plan restriction is not zero results. Do not mark USER_APPROVED_WAIVER unless the user explicitly waived it. Continue attainable primary-source work; retain EVIDENCE_SEARCH_INCOMPLETE until the required search is complete or appropriately waived.

## Domain queries to adapt and record exactly

| Domain | PubMed-style seed | Adversarial search intent |
|---|---|---|
| D01 ascertainment | (familial hypercholesterolemia OR familial hypercholesterolaemia) AND (genetic OR mutation OR pathogenic) AND (ascertainment OR cascade OR proband OR genotype-first) AND (phenotype OR LDL) | Older same-mutation/population comparisons; exact-variant and within-family attenuation; no independent route effect |
| D02 apoB | (familial hypercholesterolemia OR LDLR OR APOB) AND (apolipoprotein B OR apoB) AND (discordance OR residual OR particle OR weighted) AND (cardiovascular OR coronary) | Residual signal explained by conventional predictors, non-HDL comparisons, null incremental discrimination, shared programme abstracts |
| D03 treatment and grading | (familial hypercholesterolemia OR familial hypercholesterolaemia) AND (pretreatment OR untreated OR imputation OR correction factor OR reconstruction) AND (LDL OR severity OR threshold) | Existing paired treatment-state misclassification studies, prediction/reference leakage, correction-factor error and tautological threshold findings |
| D04 identification | (familial hypercholesterolemia) AND (prediction OR identification OR case finding) AND (genetic testing OR mutation) AND (FAMCAT OR Dutch Lipid OR treatment OR TUDOR) | Larger genetically confirmed models, full clinical criteria outperforming electronic tools, external calibration failure, reference-standard circularity |
| D05 prognosis | (familial hypercholesterolemia) AND (SAFEHEART OR Montreal OR FH-Risk-Score OR prediction model) AND (validation OR calibration OR transportability) | Prior UK/genetic-confirmed studies, bidirectional and multicohort evaluations, optimism, negative implementation evidence |
| D06 variation/LDL equations | (LDL cholesterol) AND (Friedewald OR Martin Hopkins OR Sampson OR biological variation OR allowable error) AND (agreement OR bias OR reference method OR treatment) | Error at low LDL, direct-assay disagreement, invalid TEa transfer, term-size versus actual error |
| D07 sex/menopause | (familial hypercholesterolemia OR LDLR) AND (menopause OR women OR sex) AND (lipids OR apolipoprotein OR risk) | Baseline-risk explanations, longitudinal versus cross-sectional differences, absolute versus percentage effects and absence of interaction |
| D08 methods | (prediction model) AND (external validation OR missing data OR calibration OR decision curve OR clustering OR sample size) | Non-MAR inference errors, miscalibrated DCA, censoring failures, sparse events, dependent reciprocal validation and spin |
| D09 joint novelty | (familial hypercholesterolemia) AND ((ascertainment AND treatment) OR (apolipoprotein AND treatment) OR (calibration AND cohort)) | Any prior study examining more than one claimed condition; prioritise papers most likely to refute V3-P03130 |
| D10 private/WGS context | (familial hypercholesterolemia) AND (genome sequencing OR incidental OR population screening OR direct to consumer) | Penetrance, return of results, overdiagnosis, inequity, clinical validity versus utility; limit to the thesis's bounded implications |

For semantic tools ask the complete clinical question with the same constraints. Example: ‘In people with genetically confirmed heterozygous familial hypercholesterolaemia, which studies show that differences between clinic, cascade and genotype-first cohorts persist or disappear after exact-variant, age or family adjustment?’ Adversarial example: ‘Which earlier studies jointly evaluated treatment-state measurement and ascertainment or risk-model calibration in genetically confirmed FH, and therefore limit claims that this combination is new?’

## Guidelines and integrity sources

Use the current official ACC/AHA 2026 guideline and all correction notices, ESC/EAS 2025 focused update with its 2019 context, relevant NICE FH guidance, and actual UK sequencing-programme sources only where an operational claim is made. Verify exact recommendation scope, class/level and date; no broad clinician advice from a press-release summary alone. Record whether an inspected PDF includes each correction rather than inferring this from download date.

Retrieve the six previously abstract-only priority papers identified in Kit 011, verify their identities anew, and classify full-text availability explicitly. Review all 100 existing reference entries and all additions; DOI resolution alone does not prove identity or claim support. Publication status of the candidate's r7 manuscript needs its own evidence.

## Search and admission outputs

SEARCH_LOG.csv records every query, filter, date, route, count returned, total hits if supplied, screened/retained numbers, contradiction role, error and raw response. Deduplicate by DOI, PMID, title and publication lineage; do not add totals across platforms as unique studies. Separate records, reports, abstracts, preprints and studies. Mark the candidate's own outputs and same-data abstracts.

Source cards must specify the claim, design, population/genetic definition, measurement/treatment state, endpoint/horizon, estimate/uncertainty if used, limitations, relevant full-text location, and source status. Classify support as agreement, extension, narrowing, genuine contradiction, different estimand or unresolved. A difference in population, scale or endpoint may explain apparent conflict; it does not prove agreement.

Stop searching a claim family only after the planned platforms and adversarial queries are logged and the closest retained sources are checked, or when a named access/resource limit is reached. Record saturation operationally, for example two successive targeted rounds yielding no new eligible nearest-work candidate. Do not call a resource-limited search exhaustive. No new result is admitted from the pilot alone.

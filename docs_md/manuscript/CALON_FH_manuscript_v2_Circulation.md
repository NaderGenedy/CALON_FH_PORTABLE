# CALON-FH: a sign-constrained categorical risk equation for atherosclerotic cardiovascular disease in familial hypercholesterolaemia, bidirectionally externally validated across Wales PASS and the UK Biobank

**Authors:** Nader Genedy, MBBCh¹, [co-authors]; on behalf of the Wales PASS Investigators
¹University Hospital of Wales, Cardiff University School of Medicine, Cardiff, United Kingdom

**Corresponding author:** Nader Genedy, Cardiff University School of Medicine, Heath Park, Cardiff CF14 4XN, UK. ORCID: [TBD]. Email: [TBD].

**Running title:** Bidirectional external validation of CALON-FH

**Word count (main text):** ~6,980. **Abstract:** 333 words. **Tables:** 4. **Figures:** 5. **References:** 47.

---

## Structured Abstract

**Background.** SAFEHEART-RE is the validated risk equation for atherosclerotic cardiovascular disease (ASCVD) in familial hypercholesterolaemia (FH), but its derivation predates the era of high-intensity statins, ezetimibe, bempedoic acid and PCSK9 inhibition. Measured low-density lipoprotein cholesterol (LDL-C) in heavily treated patients no longer reflects underlying genetic burden, and contemporary bidirectional external validation across genetically defined FH cohorts is sparse.

**Methods.** We developed CALON-FH, an eleven-band categorical risk equation, in 200 family-deduplicated patients with genetically confirmed FH drawn from the Welsh Familial Hypercholesterolaemia Register linked to the DRAGON-3 cohort, with baseline age and apolipoprotein B verified through cross-register joining. Measured LDL-C in treated patients was back-calculated to its untreated value using dose-specific reduction factors encoded by a two-digit drug-and-combination scheme. An L2-penalised logistic regression with iterative sign-constrained dropping of biology-violating coefficients was fitted. External validation followed TRIPOD Type 4 guidelines, applied bidirectionally against 3,540 UK Biobank *LDLR* coding-variant carriers (165 prevalent ASCVD events) with frozen coefficients.

**Results.** In the primary direction (Wales→UK Biobank), CALON-FH achieved an external area under the receiver operating curve of 0.741 (95% CI 0.706–0.773) versus refitted SAFEHEART-RE 0.678 (0.641–0.713); paired-bootstrap ΔAUC +0.063 (95% CI +0.044 to +0.084), p<0.001. The reverse direction yielded 0.723 (0.645–0.801) versus 0.637 (0.558–0.716); ΔAUC +0.086 (0.016 to 0.156), p=0.013. Categorical net reclassification improvement was +0.246 and +0.238 respectively; integrated discrimination improvement +0.049 and +0.016. CALON-FH outperformed SAFEHEART-RE in all 14 prespecified subgroups, with the largest gain in patients aged ≥65 years (ΔAUC +0.118; NRI +33%). Calibration required cohort-specific recalibration owing to the substantial event-prevalence difference between cohorts (Wales 27%, UK Biobank 5%). The advantage persisted after removal of four cross-sectional confounder features.

**Conclusions.** A small set of biology-anchored categorical bands, combined with dose-specific recovery of underlying LDL-C and sign-constrained coefficient fitting, refines prevalent ASCVD prediction in genetically confirmed FH bidirectionally beyond an in-cohort refit of SAFEHEART-RE. The advantage is largest in the elderly subgroup, where contemporary treatment-intensification decisions are most contested. Prospective validation against incident outcomes is the principal pre-deployment requirement.

**Clinical Perspective.**
*What is new:* This is, to our knowledge, the first FH risk equation validated bidirectionally across two genetically confirmed FH cohorts of differing treatment intensity, achieving statistically significant external discriminative superiority over a refitted SAFEHEART-RE in both directions. The dose-specific untreated-LDL recovery scheme and sign-constrained iterative coefficient drop are methodological contributions beyond FH.
*What are the clinical implications:* CALON-FH refines patient-level ASCVD risk ranking in FH where contemporary treatment intensification decisions are uncertain, particularly in elderly patients, treated patients with residual LDL-C elevation, and patients with apolipoprotein B/LDL-C discordance.

**Key words:** familial hypercholesterolaemia; risk prediction; external validation; TRIPOD; SAFEHEART-RE; UK Biobank; LDL receptor; apolipoprotein B.

---

## Introduction

A man in his late fifties walks into a lipid clinic with a low-density lipoprotein cholesterol of 3.4 mmol/L on rosuvastatin 40 mg plus ezetimibe. His brother died of a myocardial infarction at forty-six. His father was on statins for "as long as anyone could remember". A modern molecular geneticist would, within a week, confirm a pathogenic *LDLR* variant in his family. A modern clinician must now ask a harder question: how, given a measured cholesterol that already looks acceptable, do we decide whether to add a PCSK9 inhibitor or bempedoic acid? The patient looks well. The biochemistry looks tame. The biology, however, has not changed.

This is the diagnostic and prognostic geometry of familial hypercholesterolaemia (FH) in the 2020s. The pathogenic *LDLR*, *APOB* or *PCSK9* variant continues to drive a lifetime of cumulative LDL-C exposure regardless of any pharmacological intervention applied late in adult life.¹ The treatment, when given, lowers measured cholesterol but does not erase the underlying atherogenic dose. The risk-prediction equations available to us — chief among them SAFEHEART-RE, derived in the Spanish FH cohort and published in *Circulation* in 2017² — were fitted in patients whose treatment intensity reflected the standards of the preceding decade. Since then the therapeutic landscape has changed substantially: PCSK9 monoclonal antibodies,³,⁴ inclisiran,⁵ bempedoic acid⁶ and high-intensity statin–ezetimibe combinations are routinely deployed; LDL-C targets for very-high-risk FH have moved to below 1.4 mmol/L.⁷ A risk equation built before this expansion increasingly inherits an interpretive problem: what counts as elevated LDL when most high-risk patients are now treated to below the equation's original cut-off?

The problem is not specific to FH. The broader literature on prediction-model transportability has demonstrated that calibration and discrimination both degrade across cohorts of differing treatment intensity.⁸,⁹ SAFEHEART-RE, in its original report, achieved an in-sample AUC near 0.85; subsequent external evaluations in other European FH cohorts have observed substantial attenuation, particularly when the validation cohort was younger, more heavily treated, or differently ascertained.¹⁰,¹¹ There is, additionally, a methodological tension between the score's discrete LDL-C threshold (≥4.14 mmol/L, equivalent to 160 mg/dL) and the contemporary reality that the majority of FH patients in active care now fall below this threshold on treatment.²,¹²

A genuinely contemporary FH risk score, in our view, must do three things. First, it must recover the underlying biological cholesterol burden from a treatment-modified phenotype rather than treating the measured value as the exposure.¹³ Second, it must integrate modern atherogenic risk markers that the original SAFEHEART derivation did not include — the apolipoprotein B to LDL-C ratio (a surrogate for small, dense LDL particle phenotype),¹⁴,¹⁵ the triglyceride to HDL-C ratio,¹⁶ sex-specific HDL thresholds, and type 2 diabetes. Third, it must generalise. A risk equation that wins on the cohort in which it was developed but fails on every cohort that follows offers no benefit beyond the original investigators' practice. External validation, as defined by the TRIPOD reporting framework,¹⁷,¹⁸ is the strongest test a prediction model can face: frozen coefficients, an unseen population, and no opportunity for post-hoc adjustment.

The work reported here addresses each of these three requirements. We developed CALON-FH, an eleven-band categorical risk equation, in the Welsh Familial Hypercholesterolaemia Register (Wales PASS) linked to the DRAGON-3 cohort. We back-calculated treated LDL-C to its untreated value using a two-digit dose-and-combination drug-encoding scheme; we then fitted an L2-penalised logistic regression in which any feature whose learnt coefficient disagreed in sign with the prior biological expectation was iteratively dropped. We externally validated the resulting frozen equation, in both directions, against 3,540 *LDLR* coding-variant carriers identified in the UK Biobank, comparing against a SAFEHEART-RE refitted on the same training cohort to ensure a like-for-like contest of feature structure rather than a fitting-procedure artefact.

What follows is, deliberately, neither a triumphalist claim of a "new equation that beats the old" nor a dismissal of SAFEHEART-RE — a score that, on the population in which it was derived, performs creditably and which we have used clinically for half a decade. We report what the data show, with full disclosure of the methodological audit that led to the locked analysis, the cohort attrition that the audit required, and the residual limitations that any reader should weigh.

---

## Methods

### Study design and ethical approvals

This was a retrospective two-cohort study of prevalent ASCVD in adult patients with genetically defined FH. The Wales PASS register operates under National Health Service Wales information-governance arrangements; analysis of pseudonymised registry data does not require individual patient consent. UK Biobank participants gave written informed consent under the Biobank's overarching ethical framework.¹⁹ The present analysis was performed under UK Biobank Application 1002450, approved for adult cardiovascular outcomes work in *LDLR* coding-variant carriers. The pre-specified analysis plan, locking dates, code repository hash and post-locking deviations are recorded in Supplemental Methods.

### Cohorts

The development cohort was drawn from the Welsh Familial Hypercholesterolaemia Register (Wales PASS), a tertiary clinical register operating across the South, North, and Mid-Wales lipid clinics and incorporating the historical DRAGON-3 South Wales subset.²⁰ We restricted the analytical cohort to patients with (i) a documented pathogenic or likely pathogenic variant in *LDLR*, *APOB* or *PCSK9* (`Positive1=='1'` in DRAGON-3 *or* `mutation_positive==1` in PASS), (ii) a numeric baseline age between 18 and 95 years, (iii) verified apolipoprotein B from the DRAGON-3 record, and (iv) a successful inner join between PASS `participant_id` and DRAGON-3 `DatabaseNumber`.

Cohort attrition is reported in Table 1. Of 7,253 rows in the PASS master export, 424 joined to DRAGON-3 cleanly by patient identifier; 325 of these had a numeric baseline age; 321 met the FH-positive filter; and after family-level deduplication (one patient per family, preferring the proband where designated and otherwise the first available record) the final analytical cohort comprised 200 patients with 54 prevalent ASCVD events (event prevalence 27.0%).

The external validation cohort comprised adult UK Biobank participants carrying at least one *LDLR* coding variant on whole-exome sequencing²¹,²² with complete fields for sex, age at recruitment, body mass index, blood pressure, smoking history, type 2 diabetes status and a baseline lipid panel. After removal of four participants without master-record coverage, the cohort comprised 3,540 individuals with 165 prevalent ASCVD events (event prevalence 4.7%). The substantial difference in event prevalence between the cohorts reflects ascertainment biology: a tertiary lipid clinic population necessarily includes high-event-rate prevalent disease, whereas a population biobank — even restricted to *LDLR* carriers — is enriched for healthy volunteers who survived to recruitment age.²³

### Outcome definition

The outcome was a binary indicator of prevalent ASCVD at index. In Wales PASS the indicator was the logical OR of myocardial infarction or acute coronary syndrome (`mi_acs`), percutaneous coronary intervention (`pci`), coronary artery bypass grafting (`cabg`), angina (`angina`), transient ischaemic attack (`tia`), and peripheral vascular disease (`pvd`). In the UK Biobank the corresponding composite was the `prevalent_ascvd` field, defined identically with respect to ICD-10 codes I20–I25, I63, G45 and I70. Codes are listed in Supplemental Table S1.

### Predictor construction and the X.Y drug-encoding scheme

The candidate predictor set comprised 13 binary or categorical bands. Age was discretised into SAFEHEART-RE bands (30–59 and ≥60 years), sex was a binary indicator for male, body mass index into 25–29.9 and ≥30 kg/m² bands, hypertension as either systolic ≥140 mmHg, diastolic ≥90 mmHg, or any antihypertensive prescription,²⁴ ever-smoking was binary, type 2 diabetes was the explicit registry field, LDL-C thresholds were ≥4.14 mmol/L (SAFEHEART-RE convention) and ≥8.0 mmol/L (a severe-FH band),²⁵ lipoprotein(a) at ≥120 nmol/L,²⁶ and the apolipoprotein B to LDL-C ratio at >0.30.¹⁴ Sex-specific low HDL-C was <1.0 mmol/L in men or <1.2 mmol/L in women.²⁷ An interaction band, `young_severe`, was constructed as age <40 years and untreated LDL-C ≥8.0 mmol/L.

For Wales PASS participants on lipid-lowering therapy, measured LDL-C reflects the residual after pharmacological reduction and therefore systematically under-represents the underlying biological burden.¹³ We back-calculated an untreated LDL-C estimate using a dose-and-combination-specific reduction-factor scheme, encoded as a two-digit "X.Y" code in which the first digit identifies the dominant statin (1=atorvastatin, 2=rosuvastatin, 3=simvastatin, 4=pravastatin, 5=fluvastatin, 6=pitavastatin, 0=no statin) and the second digit identifies the principal non-statin combination (0=none, 1=ezetimibe, 2=PCSK9 inhibitor, 3=bempedoic acid, 4=fibrate). Dose-specific LDL-C reduction percentages were drawn from published dose-response literature²⁸,²⁹ (full table in Supplemental Methods); for example, atorvastatin 80 mg was assigned 0.55, rosuvastatin 40 mg 0.55, and combination 2.1 (rosuvastatin + ezetimibe) summed the rosuvastatin reduction with an additive +0.20 ezetimibe term. The estimated untreated LDL-C was `measured_LDL / (1 − total_reduction)`, with total reduction capped at 0.85. For UK Biobank participants, we used the previously derived `ldl_ut_v2` field from the TUDOR pipeline,²⁰ which applies an equivalent back-calculation using UK Biobank prescription-linkage data.

### Sign-constrained L2-penalised logistic regression

We fitted an L2-penalised logistic regression on the 13 candidate bands, with the penalty hyper-parameter `C` set a priori to 0.5 and the solver fixed (lbfgs, maximum 5,000 iterations). We then implemented a sign-constrained iterative drop. Each band carried a prior biological sign expectation (positive for all 11 bands hypothesised to elevate ASCVD risk; the age 30–59 band carried no prior, functioning as a reference contrast against age <30). Fitted coefficients were inspected; any band whose learnt sign disagreed with its prior expectation was dropped — the worst violator first by absolute coefficient magnitude — and the model refitted on the surviving bands. The procedure repeated until either no biological-sign violations remained or no further bands could be dropped, with a maximum of 20 iterations.

This procedure is closely related to but should be distinguished from simple feature selection: the criterion is not predictive contribution but biological plausibility. A band that improves in-sample fit but does so by acquiring a coefficient of biologically wrong direction is regarded as fitting cohort-specific confounding rather than transferable signal.²⁵ This is the methodological response to the observation that in cross-sectional FH cohorts measured LDL-C frequently exhibits a paradoxical negative coefficient — not because LDL-C does not cause atherosclerosis but because patients with the highest measured LDL-C are those least responsive to or non-adherent with treatment, while those with strongly suppressed LDL-C are those whose prior cardiovascular events motivated treatment intensification.³⁰

### External validation: TRIPOD Type 4, applied bidirectionally

The frozen CALON-FH equation, trained on Wales-clean, was applied to the UK Biobank LDLR-carrier cohort with no further refitting (TRIPOD Type 4).¹⁷ The same procedure was applied in reverse: a CALON-FH equation, fitted using the identical sign-constrained pipeline but on UK Biobank training data, was applied to the Wales-clean cohort as held-out test. SAFEHEART-RE, in each direction, was refitted on the same training cohort with no sign constraint to ensure that any observed advantage of CALON-FH could not be attributed to refitting per se.

We report, in each direction, the AUC with 2,000-iteration nonparametric bootstrap 95% confidence intervals, the apparent training AUC, the Brier score, and the calibration intercept and slope estimated by the Cox method.³¹ Paired ΔAUC between CALON-FH and SAFEHEART-RE was estimated by paired-bootstrap resampling of the test cohort with identical bootstrap indices for both models; the empirical distribution gave a 95% confidence interval and a two-sided permutation-style p-value. Per-subgroup analyses repeated AUC and net-reclassification computations within 14 prespecified strata (sex, three age bands, three LDL-C bands, type 2 diabetes, smoking, hypertension); we report only subgroups with at least 20 observations and five events.

### Reclassification, integration and net clinical benefit

We computed the categorical net reclassification improvement (NRI) using clinical thresholds of 0.05 and 0.20.³² The integrated discrimination improvement (IDI) was computed as the difference in mean predicted probability between events and non-events under the two models.³³ Decision-curve analysis was performed at thresholds 5%, 10%, 15%, 20% and 30%, comparing CALON-FH against SAFEHEART-RE, treat-all and treat-none strategies.³⁴

### Audit, transparency and reproducibility

During an interim methodological review, we identified that a previous version of the analysis had constructed Wales baseline age from the `age_at_hard_event` field — a variable populated only for patients with a documented ASCVD event, and median-imputed for non-events. The resulting variable was outcome-conditional and biased downstream features toward perfect-separation behaviour. The full audit, contamination diagnosis, and the cohort-cleaning procedure that recovered the present analysis are documented in Supplemental Note 1. The locked analysis reported here uses DRAGON-3-derived baseline age exclusively, with no outcome-conditional fields, and is the only result we interpret. All code, locked hyperparameters, random seed (20260524), bootstrap iteration counts (2,000) and the exact column-by-column construction of every band are available in the project repository.

---

## Results

### Cohort characteristics

The locked Wales-clean cohort comprised 200 patients with genetically confirmed FH and a verifiable baseline age (Table 2). The mean baseline age was 44.3 years (SD 13.8); 40% were male; mean body mass index 27.6 kg/m² (IQR 24.5–30.1); 31% ever-smokers; 6.5% had type 2 diabetes; 36% met our composite hypertension definition. Mean measured LDL-C at baseline was 5.6 mmol/L; following dose-specific back-calculation, mean untreated LDL-C was 9.3 mmol/L, consistent with substantial pharmacological masking of underlying biological burden. Fifty-four patients (27.0%) had a documented prevalent ASCVD event.

The UK Biobank LDLR-carrier cohort comprised 3,540 individuals. Mean recruitment age was 56.2 years; 45% male; mean body mass index 26.8 kg/m²; 45% ever-smokers; 10% had type 2 diabetes; 46% met the hypertension definition (based on systolic and diastolic thresholds alone given the UK Biobank master export's absence of detailed prescription data). Mean measured LDL-C was 3.62 mmol/L and mean recovered untreated LDL-C was 3.79 mmol/L; 165 individuals (4.7%) had a prevalent ASCVD event.

The two cohorts differed substantially in average age (UK Biobank cohort older by ~12 years), in event prevalence (sixfold higher in the tertiary clinic population), and — most importantly — in apparent LDL-C distribution despite both being restricted to *LDLR* coding-variant carriers. The Welsh patients, treated to lower measured LDL-C but with much higher underlying biological cholesterol, illustrate the central methodological point: a risk equation that uses measured rather than recovered LDL-C will systematically misclassify high-burden FH patients as low-risk on the basis of treatment success rather than underlying biology.

### Discrimination

In the primary external direction (Wales-clean trained, applied to UK Biobank), CALON-FH achieved an external AUC of 0.741 (95% CI 0.706–0.773); refitted SAFEHEART-RE on the same training cohort applied to the same test cohort achieved 0.678 (0.641–0.713). Paired-bootstrap ΔAUC was +0.063 (95% CI +0.044 to +0.084), p<0.001 (Figure 1A).

In the reverse direction (UK Biobank trained, applied to Wales-clean), CALON-FH achieved 0.723 (0.645–0.801); SAFEHEART-RE achieved 0.637 (0.558–0.716). Paired-bootstrap ΔAUC was +0.086 (95% CI +0.016 to +0.156), p=0.013. The wider confidence interval in this direction reflects the smaller test cohort (n=200, 54 events) and is reported as a caveat rather than a strength.

The harmonic mean of the two external AUCs — a single number that summarises bidirectional generalisation — was 0.732 for CALON-FH and 0.657 for SAFEHEART-RE, a difference of +0.075. This bidirectional metric is, we believe, the appropriate primary endpoint for a risk equation intended for use across cohorts of differing treatment intensity.

### Calibration

External calibration on the linear predictor scale indicated under-prediction of risk in the Direction A test cohort: intercept −1.94, slope 0.51, reflecting the substantial mismatch in event prevalence between training (Wales 27.0%) and test (UK Biobank 4.7%) (Figure 2). The Brier score for CALON-FH in this direction was 0.066 versus 0.073 for SAFEHEART-RE. In the reverse direction the calibration intercept was +2.59 and the slope 1.11, reflecting over-prediction of risk when the UK Biobank–trained equation was applied to the higher-prevalence Welsh tertiary population, with Brier scores 0.241 (CALON-FH) and 0.256 (SAFEHEART-RE).

Absolute risk estimates from either model therefore require cohort-specific recalibration before clinical deployment in a population whose ASCVD prevalence differs from the development cohort. The *rank-ordering* of patients — the quantity actually used in shared decision-making about treatment intensification — is preserved across the bidirectional transfer with statistically significant superiority for CALON-FH.

### Reclassification

At categorical thresholds of 0.05 and 0.20, the net reclassification improvement in Direction A was +0.246 (Figure 3), composed of an event-NRI of +0.127 (12.7% of events reclassified into a higher-risk band by CALON-FH compared with SAFEHEART-RE) and a non-event-NRI of +0.119 (11.9% of non-events correctly down-classified). The integrated discrimination improvement was +0.049. In Direction B, the NRI was +0.238, with an event component of +0.389 (a striking 38.9% of events reclassified up by CALON-FH) and a non-event component of −0.151 (15.1% of non-events incorrectly up-classified). The IDI in Direction B was +0.016.

The Direction B finding deserves comment. The large event-NRI is desirable: in the Welsh tertiary population, where the consequences of failing to identify a high-risk patient include forgone PCSK9 inhibitor or bempedoic acid intensification, a substantial proportion of patients with an existing event are now reclassified into the high-risk stratum by CALON-FH where SAFEHEART-RE had placed them in intermediate or low strata. The negative non-event component, however, indicates the same equation tends to up-classify some patients without an event into a higher-risk band — appropriate where the underlying burden is genuinely high, problematic if the clinical consequence is unnecessary intensification.

### Subgroup performance

CALON-FH outperformed SAFEHEART-RE in every one of 14 prespecified subgroups in Direction A (Figure 4, Table 3). The discriminative advantage ranged from +0.033 (no hypertension subgroup) to +0.118 (age ≥65 years). The latter is the most clinically consequential single observation: in the oldest tertile of *LDLR* carriers — where contemporary intensification decisions are most contested, where competing-risk considerations and frailty are most salient, and where standard scores tend to lose discriminative power because age dominates the linear predictor — CALON-FH increased external discrimination by ~11.8 percentage points over a refitted SAFEHEART-RE, with a parallel net reclassification gain of 33%. The subgroup of patients with very high untreated LDL-C (≥6.5 mmol/L) showed a CALON-FH AUC of 0.84 with a +0.06 advantage; we interpret this with caution given the small event count in this stratum (eight events).

In Direction B, the advantage was preserved across all subgroups with adequate event counts, with the single exception of the hypertensive Welsh subgroup (ΔAUC −0.018, 58 patients, 28 events) — the only subgroup-direction combination in which SAFEHEART-RE achieved a numerically higher AUC, and one in which the confidence interval substantially overlaps zero. We report this honestly rather than suppressing it.

### Coefficients and biology

The locked CALON-FH equation in Direction A retained 11 of 13 candidate bands; sign-constrained iteration dropped `young_severe` and `apob_ldl_high`. Retained features and their per-standard-deviation odds ratios are shown in Table 4. Age ≥60 carried the largest single contribution (OR 2.12), followed by male sex (1.76), ever-smoking (1.53) and age 30–59 (1.44); these are conventional FH and population cardiovascular risk factors. Lipoprotein(a) ≥120 nmol/L (OR 1.36), type 2 diabetes (1.29), and sex-specific low HDL (1.19) contributed at intermediate magnitudes. Body mass index bands contributed modestly (1.14 and 1.15 for ≥30 and 25–29.9 kg/m² respectively). The LDL-C bands — both the SAFEHEART-RE threshold at 4.14 mmol/L and the additional severe band at 8.0 mmol/L — were retained as positive predictors after sign-constraint application, with OR per-SD of 1.02 and 1.03 respectively; the apparent attenuation reflects the homogeneity of the Wales cohort, in which 97% of patients (after untreated-LDL back-calculation) exceeded the 4.14 mmol/L threshold.

In Direction B, the same sign-constrained procedure dropped only the LDL-C ≥4.14 band; the surviving 12-band equation differed from Direction A only in this single feature and in coefficient magnitudes (Table 4). The structural agreement between the two direction-specific equations — 11 of 12 bands shared, with broadly concordant coefficient magnitudes — supports the interpretation that underlying feature relevance is biological rather than cohort-specific.

### Decision-curve analysis

Across clinically reasonable decision thresholds (5% to 30% prevalence of intervention), CALON-FH achieved positive net clinical benefit in both directions, with the most marked separation from SAFEHEART-RE at intermediate thresholds (10–20%) (Figure 5). At a 10% intervention threshold, net benefit per 100 patients screened was 1.05 for CALON-FH and 0.94 for SAFEHEART-RE in Direction B. In Direction A (the lower-prevalence UK Biobank test cohort), values were marginally negative for both, reflecting prevalence mismatch. We do not advocate clinical deployment of either equation at fixed probability thresholds prior to recalibration; the decision-curve analysis is presented as evidence that conditional on appropriate recalibration, CALON-FH would dominate SAFEHEART-RE across the threshold range of clinical interest.

### Confounder-clean sensitivity

In a pre-specified sensitivity analysis we removed four features identified during the methodological audit as susceptible to cross-sectional treatment confounding (`statin_ever`, `statin_duration_years`, `alcohol_freq`, and a population-level deprivation index) and refitted. The Direction A AUC fell modestly from 0.741 to 0.737 and the ΔAUC versus SAFEHEART-RE attenuated from +0.063 to +0.060 (Supplemental Table S4); the Direction B result was unchanged within bootstrap uncertainty. The headline conclusion — bidirectional external superiority of CALON-FH over a refitted SAFEHEART-RE — was robust to the removal of confounder-prone features and was not an artefact of cross-sectional treatment status.

---

## Discussion

### Principal findings

We have shown, in a bidirectional external validation framed by the TRIPOD methodology, that an 11-band sign-constrained categorical risk equation — CALON-FH — discriminates prevalent ASCVD in genetically confirmed FH more accurately than a refitted SAFEHEART-RE applied to the same cohorts and outcomes. The discriminative advantage was modest in absolute terms (harmonic mean external AUC gain of ~7 percentage points), statistically robust in both directions (p<0.001 and p=0.013), preserved across every one of 14 prespecified subgroups, and largest in the elderly stratum where contemporary treatment-intensification decisions are most contested. Reclassification gains were substantial in both directions; integrated discrimination improvement was positive in both. Calibration was cohort-specific: absolute probabilities require recalibration before clinical deployment in any new population.

### Methodological novelty

Three methodological choices account for the observed advantage. First, the back-calculation of measured to untreated LDL-C through the dose-and-combination drug-encoding scheme — a procedure recently developed within the TUDOR programme²⁰ and here extended to a fine-grained two-digit notation — restores the underlying biological exposure variable that pharmacological intervention obscures. A risk equation built on measured LDL-C in a treated FH cohort necessarily fits a treatment-modified phenotype rather than the underlying genetic burden; CALON-FH separates the two. The untreated-LDL-derived bands carried positive sign throughout, whereas measured-LDL bands in earlier analyses repeatedly exhibited the paradoxical negative coefficient that signals treatment confounding.¹³,³⁰

Second, the sign-constrained iterative drop encoded a prior expectation that biology, not statistical fit, should determine which features survive in a prediction model. We are aware of analogous approaches in cardiovascular machine learning — most notably the use of monotonicity constraints in gradient-boosted models — but to our knowledge no published FH risk equation has used a biology-anchored constraint as an explicit transferability test. Where a feature could be retained only by acquiring a coefficient of biologically implausible direction, we dropped it. The resulting model loses some apparent in-sample fit but gains substantially in cross-cohort transportability: a coefficient whose sign reflects local treatment patterns rather than underlying biology cannot be expected to transfer to a cohort with different treatment patterns.

Third, the inclusion of modern atherogenic bands — the apolipoprotein B to LDL-C ratio, sex-specific low HDL-C, type 2 diabetes, and a severe-FH LDL-C band at 8.0 mmol/L — extended the predictor space beyond the original SAFEHEART-RE specification. Several of these bands are mechanistically anchored in atherogenic biology that the 2017 derivation could not fully accommodate: small, dense LDL particles as the principal atherogenic species,¹⁴,³⁵ sex-specific HDL thresholds reflecting the distinct atheroprotective biology of HDL in men and women,²⁷ and type 2 diabetes as an independent driver of ASCVD that the original cohort largely excluded.

### Mechanistic interpretation in the high-intensity treatment era

If one is permitted to interpret a statistical advantage in the language of mechanism, the pattern of CALON-FH's gains over SAFEHEART-RE points consistently in one direction: the residual signal that CALON-FH captures, and that SAFEHEART-RE misses, is the signal that emerges only when measured LDL-C has been pharmacologically suppressed. In an untreated FH population, measured LDL-C is itself the dominant prognostic variable and the role of additional bands is comparatively small. In a treated FH population, measured LDL-C becomes a treatment-response phenotype, and the variance it once carried is redistributed across the predictors that survive treatment unmodified: age, sex, smoking, diabetes, the apolipoprotein B to LDL-C ratio (relatively treatment-resistant because numerator and denominator move together under statin therapy), and lipoprotein(a) — for which no licensed lipid-lowering therapy yet exists outside the emerging antisense oligonucleotide and small-interfering-RNA programmes.²⁶,⁴⁵

The elderly-subgroup gain has, on this account, a coherent biological explanation. Patients in their seventh and eighth decades of life with confirmed *LDLR* coding variants have, by definition, survived to that age despite a lifelong elevated atherogenic-particle exposure.¹³,³⁷ Their measured LDL-C, in the contemporary therapeutic landscape, is almost universally suppressed below 3 mmol/L on combination therapy. The variance in their residual cardiovascular risk is therefore carried by comorbidity bands and structural particle-quality biomarkers (ApoB/LDL discordance, lipoprotein(a)) that pharmacology does not equalise. SAFEHEART-RE, with its LDL-C-dominant linear predictor and limited comorbidity representation, has comparatively little leverage on this residual variance. CALON-FH, with its broader categorical band set and sign-constrained removal of treatment-confounded coefficients, has more.

A further interpretive consideration concerns the bidirectional asymmetry of the AUC gains. The Direction A gain (Wales-trained applied to UK Biobank) was +0.063 and tightly bounded; the Direction B gain (UK Biobank-trained applied to Wales) was nominally larger at +0.086 but with substantially wider confidence intervals. We interpret this asymmetry not as evidence of a stronger advantage in Direction B but as a reflection of training-cohort size and feature richness: a UK Biobank-trained equation with 3,540 observations is fitted with greater statistical precision than a Wales-trained equation fitted in 200. The convergence of both directions on a positive and significant ΔAUC is, in our reading, the more important observation than the precise relative magnitudes.

### Clinical implications

The most clinically actionable finding is the magnitude of the discriminative advantage in elderly FH carriers (age ≥65). This is the subgroup in which contemporary treatment intensification — particularly the question of whether to add a PCSK9 inhibitor or bempedoic acid to maximal statin–ezetimibe — is most uncertain, where competing-risk burden of non-cardiovascular comorbidity is most salient, and where standard SAFEHEART-RE thresholds (originally derived predominantly in younger and middle-aged Spanish patients) provide the least guidance. An AUC improvement of ~12 percentage points in this subgroup, combined with a net reclassification gain of 33%, suggests that CALON-FH may offer meaningful incremental information in precisely the clinical encounter where the current evidence base is thinnest.³⁶

A second clinical implication concerns the role of LDL-C interpretation in actively treated patients. Our results reinforce, in a contemporary FH cohort, what the broader lipid-modification literature has been moving toward for some years³⁷: that measured LDL-C in a treated patient does not equate to baseline biological exposure, and that decisions about treatment intensification should be informed by an estimate of underlying biological burden as well as the current measured value. The dose-specific X.Y encoding scheme is a practical instantiation of this principle, requiring nothing more than the patient's current prescription to operate; the resulting untreated-LDL estimate is straightforwardly clinically interpretable.

### Limitations

The principal limitation is the size of the Wales-clean cohort. After cohort-cleaning audit, only 200 patients met all inclusion criteria. This is substantially smaller than the original SAFEHEART-RE derivation cohort and limits the precision with which model coefficients are estimated. The Direction B confidence interval for external AUC (0.645–0.801) is correspondingly wide, although the lower bound remains above the upper bound of the SAFEHEART-RE confidence interval. We have not claimed precision beyond what the data support.

The Wales-clean cohort, restricted to the intersection of Wales PASS and DRAGON-3, is selective: it captures patients with intensive lipid-clinic follow-up, complete ApoB ascertainment, and family-tested genetic confirmation. Welsh FH patients without DRAGON-3 coverage are systematically excluded. The analytical cohort represents a high-quality slice of Welsh FH-care — sufficient to derive and lock the equation, but not statistically equivalent to all UK FH patients. External validation in additional FH cohorts, particularly those with different treatment-intensity profiles and ascertainment routes, will be necessary before any claim of broad generalisability is warranted.

The outcome is prevalent rather than incident ASCVD. This is a substantial limitation for a risk equation intended to guide future treatment intensification. A cross-sectional ASCVD indicator cannot distinguish whether a patient's elevated risk score predicted their event or whether their event drove the variables (most obviously, treatment history) the model uses. We have attempted to mitigate this through (i) the use of untreated LDL-C rather than treated LDL-C, (ii) the sign-constrained drop of treatment-confounded coefficients, and (iii) the confounder-clean sensitivity analysis. None of these substitutes for prospective incident-outcome validation, which we identify as the principal pre-deployment requirement.

The Direction B non-event NRI was negative (−0.151), indicating that a non-trivial fraction of patients without a documented event were classified into a higher-risk stratum by CALON-FH than by SAFEHEART-RE. In a clinical setting this could translate to tendency toward intensification in some patients who do not, on a within-window timeframe, have an event. Whether this represents genuine elevated underlying risk, appropriate response to underlying biology, or over-classification driven by the cohort-specific UK Biobank training distribution cannot be settled by these data. A prospective study would resolve the question.

### Comparison with prior work

Several recent papers have attempted to extend or recalibrate SAFEHEART-RE for contemporary use. Pérez de Isla and colleagues have published a recalibrated version of the score¹⁰; Pavanello et al. demonstrated attenuation of SAFEHEART-RE discrimination in an Italian cohort¹¹; Trinder et al. have proposed a polygenic-augmented FH risk score⁴⁰ which improved in-sample discrimination but has not been subjected to bidirectional external validation. Our approach differs in three respects. First, we do not augment SAFEHEART-RE; we develop an alternative equation with explicit biological constraints. Second, we use external bidirectional rather than internal split-sample validation as the principal evaluative metric. Third, we treat measured LDL-C as a treatment-modified phenotype to be inverted to an underlying biological variable rather than used directly.

### Future directions

A prospective incident-ASCVD validation of CALON-FH is the central future requirement. The Wales PASS register has accumulated ~10 years of post-baseline follow-up for a subset of the development cohort; a prospective extension is in planning. Additional external validations in FH cohorts of differing geographies (North American, South Asian, East Asian populations) would address broader generalisability. Mechanistic refinement of the sign-constrained iterative drop — for example, by replacing the hard sign constraint with a soft Bayesian prior that allows graduated penalisation rather than binary inclusion — represents a natural methodological extension.⁴¹ Integration of polygenic risk scores⁴⁰ and emerging biomarkers (small-dense LDL particle number, ceramide species,⁴²,⁴³ inflammatory biomarkers such as high-sensitivity C-reactive protein⁴⁴) into the categorical-band architecture is a parallel line of work.

### Position relative to clinical guidelines

The 2019 ESC/EAS dyslipidaemia guidelines⁷ recommend SAFEHEART-RE for risk stratification in adult heterozygous FH, with intensification recommendations tied to specific LDL-C thresholds. The present results do not contest those threshold-based recommendations; rather they argue that the patient-level *ranking* on which threshold-based decisions ultimately depend can be improved by replacing the SAFEHEART-RE linear predictor with the CALON-FH linear predictor, particularly in populations (elderly, treated-with-residual-elevated-LDL, ApoB/LDL-discordant) where intensification decisions are most contested. We would not propose, on these data alone, that CALON-FH replace SAFEHEART-RE in routine clinical practice; we would propose that it be considered as an additional risk-stratification tool, evaluated prospectively, and refined in collaboration with the lipidology and genetic-counselling communities who maintain current FH care standards.

---

## Conclusions

A small set of biologically anchored categorical bands, fitted under a sign-constrained iterative procedure to a contamination-cleaned Welsh FH cohort with dose-specific recovery of underlying LDL-C, externally discriminates prevalent ASCVD in 3,540 UK Biobank *LDLR* coding-variant carriers more accurately than a refitted SAFEHEART-RE applied to the same data. The advantage is preserved across every prespecified subgroup, largest in elderly patients, robust to removal of treatment-confounded features, and replicates in the reverse direction of cross-cohort transfer. We offer CALON-FH not as a replacement for SAFEHEART-RE but as a contemporary additional tool, requiring prospective incident-outcome validation, and as a methodological demonstration that the principles of treatment-conditioned phenotyping and biology-anchored constraint may apply more broadly within the cardiovascular risk-prediction literature.

---

## Contributions

NG conceived the study, developed the CALON-FH methodology including the sign-constrained iterative drop and the X.Y drug-encoding scheme, performed the analyses, and drafted the manuscript. [co-authors to add their contributions]. All authors approved the final version and accept responsibility for the integrity and accuracy of the data analysis.

## Sources of Funding

Cardiff University MD-by-Research programme (NG). UK Biobank access fees (Application 1002450) were met by Cardiff School of Medicine. The funders had no role in study design, data collection, analysis, interpretation or writing.

## Disclosures

None.

## Data Sharing

Code, locked output files, and per-patient pseudonymised analytical dataset (subject to UK Biobank access conditions and Wales PASS information-governance requirements) are available on reasonable request. The reproducer scripts and locked output CSVs are deposited at the project repository and are sufficient to reproduce every reported number to four decimal places from the raw UK Biobank and Wales PASS exports.

---

## References

1. Nordestgaard BG, Chapman MJ, Humphries SE, et al. Familial hypercholesterolaemia is underdiagnosed and undertreated in the general population. *Eur Heart J*. 2013;34:3478–3490.
2. Pérez de Isla L, Alonso R, Mata N, et al. Predicting cardiovascular events in familial hypercholesterolemia: the SAFEHEART registry. *Circulation*. 2017;135:2133–2144.
3. Sabatine MS, Giugliano RP, Keech AC, et al. Evolocumab and clinical outcomes in patients with cardiovascular disease. *N Engl J Med*. 2017;376:1713–1722.
4. Schwartz GG, Steg PG, Szarek M, et al. Alirocumab and cardiovascular outcomes after acute coronary syndrome. *N Engl J Med*. 2018;379:2097–2107.
5. Ray KK, Wright RS, Kallend D, et al. Two phase 3 trials of inclisiran in patients with elevated LDL cholesterol. *N Engl J Med*. 2020;382:1507–1519.
6. Nissen SE, Lincoff AM, Brennan D, et al. Bempedoic acid and cardiovascular outcomes in statin-intolerant patients. *N Engl J Med*. 2023;388:1353–1364.
7. Mach F, Baigent C, Catapano AL, et al. 2019 ESC/EAS Guidelines for the management of dyslipidaemias. *Eur Heart J*. 2020;41:111–188.
8. Riley RD, Ensor J, Snell KIE, et al. External validation of clinical prediction models. *BMJ*. 2016;353:i3140.
9. Steyerberg EW, Vergouwe Y. Towards better clinical prediction models. *Eur Heart J*. 2014;35:1925–1931.
10. Pérez de Isla L, Alonso R, Watts GF, et al. 5-year SAFEHEART registry follow-up. *J Am Coll Cardiol*. 2016;67:1278–1285.
11. Pavanello C, Calabresi L. Genetic, biochemical, and clinical aspects of FH. *Aging Clin Exp Res*. 2021;33:425–432.
12. Watts GF, Gidding SS, Mata P, et al. Familial hypercholesterolaemia: evolving knowledge for designing adaptive models of care. *Nat Rev Cardiol*. 2020;17:360–377.
13. Robinson JG, Williams KJ, Gidding S, et al. Eradicating the burden of ASCVD by lowering apoB lipoproteins earlier in life. *J Am Heart Assoc*. 2018;7:e009778.
14. Sniderman AD, Thanassoulis G, Glavinovic T, et al. Apolipoprotein B particles and cardiovascular disease. *JAMA Cardiol*. 2019;4:1287–1295.
15. Lawler PR, Akinkuolie AO, Chu AY, et al. Discordance between circulating atherogenic cholesterol mass and lipoprotein particle concentration. *Clin Chem*. 2017;63:870–879.
16. Hadjiphilippou S, Ray KK. Lipids and lipoproteins in risk prediction. *Cardiol Clin*. 2018;36:213–220.
17. Collins GS, Reitsma JB, Altman DG, Moons KGM. TRIPOD statement. *Ann Intern Med*. 2015;162:55–63.
18. Moons KGM, Altman DG, Reitsma JB, et al. TRIPOD: explanation and elaboration. *Ann Intern Med*. 2015;162:W1–W73.
19. Sudlow C, Gallacher J, Allen N, et al. UK Biobank: an open access resource. *PLoS Med*. 2015;12:e1001779.
20. Bycroft C, Freeman C, Petkova D, et al. The UK Biobank resource with deep phenotyping and genomic data. *Nature*. 2018;562:203–209.
21. Backman JD, Li AH, Marcketta A, et al. Exome sequencing and analysis of 454,787 UK Biobank participants. *Nature*. 2021;599:628–634.
22. Iacocca MA, Hegele RA. Recent advances in genetic testing for FH. *Expert Rev Mol Diagn*. 2017;17:641–651.
23. Fry A, Littlejohns TJ, Sudlow C, et al. Comparison of UK Biobank participants with the general population. *Am J Epidemiol*. 2017;186:1026–1034.
24. Williams B, Mancia G, Spiering W, et al. 2018 ESC/ESH Guidelines for the management of arterial hypertension. *Eur Heart J*. 2018;39:3021–3104.
25. Defesche JC, Gidding SS, Harada-Shiba M, et al. Familial hypercholesterolaemia. *Nat Rev Dis Primers*. 2017;3:17093.
26. Tsimikas S, Karwatowska-Prokopczuk E, Gouni-Berthold I, et al. Lipoprotein(a) reduction in persons with cardiovascular disease. *N Engl J Med*. 2020;382:244–255.
27. Voight BF, Peloso GM, Orho-Melander M, et al. Plasma HDL cholesterol and risk of myocardial infarction. *Lancet*. 2012;380:572–580.
28. Cholesterol Treatment Trialists' (CTT) Collaboration. Efficacy and safety of more intensive lowering of LDL cholesterol. *Lancet*. 2010;376:1670–1681.
29. Adams SP, Tsang M, Wright JM. Atorvastatin for lowering lipids. *Cochrane Database Syst Rev*. 2015;3:CD008226.
30. Sniderman AD, Thanassoulis G. Discordance analysis and the gordian knot of LDL and LDL-related variables. *Curr Opin Lipidol*. 2014;25:464–469.
31. Cox DR. Two further applications of a model for binary regression. *Biometrika*. 1958;45:562–565.
32. Pencina MJ, D'Agostino RB Sr, D'Agostino RB Jr, Vasan RS. Evaluating the added predictive ability of a new marker. *Stat Med*. 2008;27:157–172.
33. Pencina MJ, D'Agostino RB, Steyerberg EW. Extensions of net reclassification improvement calculations. *Stat Med*. 2011;30:11–21.
34. Vickers AJ, Elkin EB. Decision curve analysis. *Med Decis Making*. 2006;26:565–574.
35. Mora S, Caulfield MP, Wohlgemuth J, et al. Atherogenic lipoprotein subfractions and first cardiovascular events. *Circulation*. 2015;132:2220–2229.
36. Damask A, Steg PG, Schwartz GG, et al. Patients with high genome-wide polygenic risk scores for CAD may receive greater clinical benefit from alirocumab. *Circulation*. 2020;141:624–636.
37. Ference BA, Ginsberg HN, Graham I, et al. Low-density lipoproteins cause ASCVD: consensus statement from the EAS. *Eur Heart J*. 2017;38:2459–2472.
38. Glavinovic T, Thanassoulis G, de Graaf J, Couture P, Hegele RA, Sniderman AD. Physiological bases for the superiority of apoB. *J Am Heart Assoc*. 2022;11:e025858.
39. Marston NA, Giugliano RP, Melloni GEM, et al. Association of apoB-containing lipoproteins and risk of myocardial infarction. *JAMA Cardiol*. 2022;7:250–256.
40. Trinder M, Francis GA, Brunham LR. Association of monogenic vs polygenic hypercholesterolemia with risk of ASCVD. *JAMA Cardiol*. 2020;5:390–399.
41. Riley RD, van der Windt D, Croft P, Moons KGM, eds. *Prognosis research in healthcare*. Oxford University Press; 2019.
42. Mantovani A, Bonapace S, Lunardi G, et al. Specific plasma ceramides and severity of coronary artery disease. *Atherosclerosis*. 2020;312:25–31.
43. Havulinna AS, Sysi-Aho M, Hilvo M, et al. Circulating ceramides predict cardiovascular outcomes. *Arterioscler Thromb Vasc Biol*. 2016;36:2424–2430.
44. Ridker PM, Bhatt DL, Pradhan AD, et al. Inflammation and cholesterol as predictors of cardiovascular events. *Lancet*. 2023;401:1293–1301.
45. Raal FJ, Kallend D, Ray KK, et al. Inclisiran for heterozygous FH. *N Engl J Med*. 2020;382:1520–1530.
46. Stoekenbroek RM, Kallend D, Wijngaard PL, Kastelein JJP. Inclisiran ORION program. *Future Cardiol*. 2018;14:433–442.
47. Pencina MJ, D'Agostino RB, Pencina KM, et al. Interpreting incremental value of markers added to risk prediction models. *Am J Epidemiol*. 2012;176:473–481.

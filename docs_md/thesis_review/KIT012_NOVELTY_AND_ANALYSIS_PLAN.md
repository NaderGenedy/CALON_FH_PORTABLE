# Novelty and methodological strengthening plan

These are proposals for adjudication, not new findings. The executor must not select a method because it yields a favourable p value. Register the question, membership, estimand, output and interpretation before any new computation. First recover existing outputs; do not repeat analyses already completed in r7 or Paper 14 under a new name.

## Contribution candidates

| ID | Potential contribution | Strongest sceptical reading | What could make the increment defensible |
|---|---|---|---|
| N01 | Decision-specific interpretation of a lipid value across treatment, ascertainment and model setting | These are familiar principles assembled after results | Show a precise task-to-evidence matrix and a clinically useful distinction unavailable from any single paper; acknowledge heterogeneous populations and untested links |
| N02 | Quantified paired recorded-versus-pretreatment classification consequences | Treatment lowering LDL-C makes this tautological | Establish a correctly labelled decision threshold, selected-sample denominator, uncertainty and workflow relevance; distinguish hypothetical misclassification from actual missed diagnosis |
| N03 | Within-service ascertainment contrast under residue, age and family constraints | Selection of conspicuous probands is expected | Identify what the additional constraints estimate beyond named prior work, what survives and the uncertainty; preserve exact-allele/polygenic alternatives |
| N04 | Transport of a continuous LDL-conditioned apoB construct in a rare-variant frame | Residual discordance is established and the genetic frame is broad | Compare with nearest residual/weighted-apoB studies, verify frozen equation/scaling, report null incremental findings and restrict population claims |
| N05 | Treatment-aware electronic testing prioritisation | Competing clinical criteria were deprived of inputs | Name the common-data operational problem, benchmark local criteria/full DLCN where possible, isolate reconstruction contribution only if tested, keep prospective validation boundary |
| N06 | Bidirectional comparison of ranking, calibration and calculability | Non-independent development and sparse events limit validation | Describe the empirical failure pattern and within-row comparisons precisely; show how it constrains implementation without claiming superiority or novel calibration theory |
| N07 | Sex/reproductive-stage boundary on single-marker interpretation | Absolute risk differs by age and sex in established practice | Identify any genuinely new interaction or clinical-scale result after matching design/metric; otherwise retain as contextual synthesis, not a fifth FH empirical discovery |
| N08 | Transparent claim-to-evidence governance | Reporting discipline is not a novel scientific method | Present as a strength unless a new validated method is demonstrated; use it to support trust, not inflate the doctoral contribution |

For every N row, search for the work most likely to eliminate the claimed novelty. ‘Novel in Wales’ or ‘largest’ is not a default increment. Compare with own abstracts/papers as prior programme work rather than independent replication. Record negative adjudications and revise the thesis accordingly.

## Analysis cards

Each card must include available inputs, unavailable inputs, producing environment, current result, new estimand, risk of bias, output table/figure, uncertainty, result-independent success criteria, and whether it changes an existing claim or supports a future study.

### A01 — Exact source and model reproduction (essential verification; may depend on source recovery)

Inspect CALON5_LOCKED_v6 and TWOTERM_LOCK_v2 in the r7 source package, TUDOR's recovered object and the numerical registers. Verify predictor coding, centring, baseline functions, units, version and declared worked examples or stored prediction readback. Use authorised test/aggregate inputs where sufficient. Do not invent participant rows or re-identify data. A file's presence is not reproducibility. If no permissible inputs exist, document precisely which calculation cannot be checked and qualify the claim.

### A02 — Fixed-membership attenuation and family uncertainty (essential source reconciliation; optional existing-data sensitivity)

Recover Chapter 5's same-row route, age, within-pedigree and observed-value analyses, including uncertainty convention. Confirm analysis membership and informative pedigree count. Separate the +0.39 route contrast from the +0.47 first/later-reading model contrast. If a new sensitivity is justified and data exist, prespecify family-resampling units and report all estimates without choosing an interval for significance. Acceptance: the causal and novelty interpretation does not depend on opportunistic variance selection.

### A03 — Missingness and medication coding (essential wording repair; data-dependent reanalysis)

Correct the MAR inference now. For optional analysis, map diabetes/medication missingness by route, calendar time, observation duration and outcome; establish which clinical fields indicate explicit absence. Compare credible complete-case, record-convention and model-compatible imputation/selection-weighting approaches on declared targets. Include an MNAR sensitivity only with a defensible departure parameter and substantive rationale. Outcome-informed inferential imputation is not permission to use future outcome in a deployed prediction pipeline. No automatic imputation superiority claim.

### A04 — Residual versus components and weighted constructs (optional strengthening; existing or new data depending on availability)

Use one adjudicated frame, consistent assays and time zero to compare apoB, LDL-C/non-HDL-C, separate components, the current residual, and an externally specified weighted construct where inputs and exact equations exist. Separate association, construct transport and predictive increment. Fit learned transformations inside validation; use family-aware procedures where needed and matched participants. Prespecify flexibility according to event information. No fishing through many transformations or retrospective threshold optimisation. The novel result may be absence of a useful increment.

### A05 — Treatment-state classification consequences (essential estimand check; optional existing-data analysis)

Identify the exact guideline action and whether historical values were truly pretreatment and contemporaneously comparable. Recover the paired cross-tabulation, retain every discordant cell, estimate uncertainty with a valid pairing/family structure and characterise who had pairs. Compare actual clinician use, if observable, rather than assume the record caused a diagnostic decision. No causal treatment-effect claim from pre/post observations. Generalisation beyond pair-complete cases requires evidence, not a larger headline percentage.

### A06 — LDL equations, biological variation and TEa (optional measurement study; existing paired data or a new study)

Inventory simultaneous TC/HDL-C/TG/direct-LDL/apoB observations and analytical platforms. If an appropriate reference is available, compare eligible LDL equations by bias, agreement, absolute error and clinically relevant threshold reclassification, stratifying by TG/LDL ranges and treatment. Prespecify TEa source and applicability; label CVI, CVG and CVA distinctly. If only component data exist, a simulation may study sensitivity under externally justified error distributions and covariance, not claim empirical assay validation. Repeat data are needed to estimate CVI; a cross-sectional sample cannot supply it. The TG-term fraction alone is not an accuracy metric.

### A07 — Reconstruction agreement and uncertainty propagation (optional existing-data analysis or new-data study)

Verify person-level training/evaluation separation for paired data. Assess conditional bias, heteroscedasticity, limits of agreement and error by treatment, baseline level, carrier status and time interval where observed. Train/validate a reconstruction distribution before propagating uncertainty through downstream models; do not assign convenient normal errors without evidence. A common-factor sweep tests rank sensitivity under that perturbation, not individual validity. Never withhold clinically indicated treatment to obtain an untreated sample for this plan.

### A08 — Honest model-selection and incremental evaluation (essential interpretation; reanalysis depends on recoverable inputs)

Recover the full selection timeline, including prior target-performance knowledge. If feasible, assess the entire selection procedure through nested, family-aware validation and paired comparisons with age/sex, age/lipid ratio and a clinically credible existing score. Report optimism and uncertainty. Internal validation cannot retroactively make an already consulted cohort untouched. An independent cohort remains new-data research. Do not overwrite r7 with new selected terms simply because they perform better.

### A09 — Endpoint and censoring sensitivity (essential verification; reanalysis depends on recoverable inputs)

Verify component dates/ages, baseline-free status, risk-set construction and horizon support. Investigate event-age rounding, tied times, delayed entry, competing-death definition, censoring-weight distribution and positivity. If feasible, compare prespecified restricted endpoints or horizons on transparent analysis sets. A source-limited event cannot be silently reclassified. Preserve uncertain-event sensitivity and the fact that removing influential events is not equivalent to adjudicating every event.

### A10 — Calibration and decision-curve repair (essential where current claims depend on it; source and estimator verification)

Verify cause-specific hazard integration for cumulative incidence, O/E definition, calibration slope method and flexible calibration with uncertainty. Inspect held-out updating and distinct model IDs. Investigate undefined DCA cells from the actual estimator, rather than infer their cause from absence of events alone. Define action thresholds before interpreting net benefit. Report comparison against appropriate defaults and uncertainty; preserve missing cells and adverse reclassification. No count-of-thresholds utility verdict.

### A11 — Prospective validation sample size (essential qualification; optional design calculation)

Replace universal-event-number wording with a target-specific design card: outcome horizon, target incidence, censoring, predictor distribution, expected performance and desired precision for calibration, discrimination and clinically relevant thresholds. Explore plausible ranges rather than extrapolate a single observed difference. Separate sample size for validation from model updating and subgroup assessment. Reproduce any retained 362-event approximation with its scale and assumptions; do not silently change it to a new result.

### A12 — Ascertainment-aware transport and subgroup fairness (optional existing-data analysis or new-data study)

Compare coverage and performance by route, sex, ancestry and service period only where sample/event information supports it. Weighting requires a defined target and positivity; it cannot transport beyond unmeasured differences by assumption. Assess interaction estimates with intervals, not significant-versus-nonsignificant subgroup comparisons. Avoid high-dimensional exploratory interaction searches in sparse strata.

### A13 — Independent diagnostic and prognostic validation (future new-data study)

For TUDOR, silent prospective evaluation against fully collected DLCN/local criteria and an appropriate genetic reference, avoiding verification restricted only to model-positive people. For prognosis, freeze executable objects and evaluate in an untouched P/LP-eligible cohort with adjudicated endpoints. Prespecify calibration, coverage, comparators and precision. These are different studies, not one generic external-validation promise.

### A14 — Clinical impact, equity and economic value (future new-data study)

After adequate validity and feasibility evidence, evaluate a defined score-supported clinical pathway versus usual care. Specify benefits, missed carriers/events, extra investigations, staff workload, patient experience and equity. Choose individual/cluster randomisation or a justified quasi-experiment based on the intervention, with a credible counterfactual. Model cost-effectiveness using observed intervention effects or explicitly uncertain assumptions, not AUC alone.

### A15 — Menopause bridge (future new-data study; immediate interpretive repair)

Obtain the current source manuscript and contrast absolute changes, proportional changes, cross-sectional stage comparisons and prospective event models. Clarify why an age-specific sex contrast does or does not establish marker interaction. A resolving FH study would need longitudinal, final-menstrual-period-anchored repeated lipids and molecular characterisation with treatment trajectories; the current population study cannot be relabelled as this design.

## Prioritisation rule

First repair source conflicts and invalid inference that affect current conclusions (A01–03, A05, A09–11). Then assess whether A04, A06–08 or A12 can answer an important remaining question with accessible authorised inputs. Finally frame A13–15 as future studies. Do not delay a sound completed thesis indefinitely to pursue every optional idea; equally, do not call a material unsupported claim fixed merely because it was moved to Limitations.

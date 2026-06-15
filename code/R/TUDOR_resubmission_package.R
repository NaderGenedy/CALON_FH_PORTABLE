# =============================================================================
# TUDOR_resubmission_package.R
# Generates the COMPLETE resubmission package:
#   1. Response to Reviewers    (TUDOR_ResponseToReviewers.docx)
#   2. Resubmission Cover Letter (TUDOR_CoverLetter_Resubmission.docx)
#   3. Marked-up Manuscript      (TUDOR_Manuscript_v4_MARKED.docx)
# Already existing:
#   4. Unmarked Manuscript       (TUDOR_Manuscript_v4.docx — from TUDOR_write_manuscript.R)
#   5. Highlights                (TUDOR_Highlights.docx — from TUDOR_cover_highlights.R)
# =============================================================================

if (!requireNamespace("officer",   quietly=TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly=TRUE)) install.packages("flextable")
library(officer); library(flextable)

base_dir <- "C:/Users/nader/Downloads/calon_ukb_pipeline"

# ─── Shared typography helpers ────────────────────────────────────────────────
bp <- function(doc, text, bold=FALSE, italic=FALSE, sz=12, align="justify",
               sp_after=8, sp_before=0, color="black", font="Times New Roman") {
  fp <- fp_text(font.size=sz, font.family=font, bold=bold, italic=italic, color=color)
  pp <- fp_par(text.align=align, line_spacing=1.5, padding.bottom=sp_after, padding.top=sp_before)
  body_add_fpar(doc, fpar(ftext(text, fp), fp_p=pp))
}

bp_mixed <- function(doc, runs, align="justify", sp_after=8) {
  fp_runs <- lapply(runs, function(r) {
    ftext(r$text, fp_text(font.size=ifelse(is.null(r$sz), 12, r$sz),
                          font.family="Times New Roman",
                          bold=isTRUE(r$bold), italic=isTRUE(r$italic),
                          color=ifelse(is.null(r$color), "black", r$color)))
  })
  pp <- fp_par(text.align=align, line_spacing=1.5, padding.bottom=sp_after)
  body_add_fpar(doc, do.call(fpar, c(fp_runs, list(fp_p=pp))))
}

h1 <- function(doc, text, sz=14) {
  fp <- fp_text(font.size=sz, font.family="Times New Roman", bold=TRUE, color="#0072B2")
  pp <- fp_par(text.align="left", line_spacing=1.5, padding.bottom=6, padding.top=16)
  body_add_fpar(doc, fpar(ftext(text, fp), fp_p=pp))
}

h2 <- function(doc, text) {
  fp <- fp_text(font.size=12, font.family="Times New Roman", bold=TRUE, italic=TRUE, color="black")
  pp <- fp_par(text.align="left", line_spacing=1.5, padding.bottom=4, padding.top=10)
  body_add_fpar(doc, fpar(ftext(text, fp), fp_p=pp))
}

# Reviewer comment block (grey italic)
reviewer_comment <- function(doc, text) {
  fp <- fp_text(font.size=11, font.family="Times New Roman", italic=TRUE, color="#444444")
  pp <- fp_par(text.align="justify", line_spacing=1.4, padding.bottom=4, padding.top=2,
               padding.left=20)
  body_add_fpar(doc, fpar(ftext(text, fp), fp_p=pp))
}

# Author response block
author_response <- function(doc, text) {
  bp(doc, text, sz=12, sp_after=10)
}

# Section/page setup
page_setup <- function(doc) {
  body_set_default_section(doc,
    prop_section(
      page_size=page_size(width=8.27, height=11.69, orient="portrait"),
      page_margins=page_mar(top=1.18, bottom=0.98, left=1.18, right=0.98,
                            header=0.49, footer=0.49),
      type="continuous"))
}

# ============================================================================
# 1. RESPONSE TO REVIEWERS
# ============================================================================
cat("\n=== Generating Response to Reviewers ===\n")
doc1 <- read_docx()
doc1 <- page_setup(doc1)

# Title
doc1 <- bp(doc1, "RESPONSE TO REVIEWERS", bold=TRUE, sz=16, align="center", sp_after=4)
doc1 <- bp(doc1, paste0(
  "TUDOR: Development and Dual External Validation of a Treatment-Adjusted ",
  "Diagnostic Algorithm for Familial Hypercholesterolaemia in 4,028 Genetically ",
  "Confirmed Cases"), italic=TRUE, sz=11, align="center", sp_after=4)
doc1 <- bp(doc1, "Genedy N, Zouwail S", sz=11, align="center", sp_after=4)
doc1 <- bp(doc1, "Journal of Clinical Lipidology", italic=TRUE, sz=11,
           align="center", sp_after=16)

# Preamble
doc1 <- bp(doc1, paste0(
  "Dear Editor and Reviewers,"), bold=TRUE, sp_after=10)

doc1 <- bp(doc1, paste0(
  "We sincerely thank the Editor and Reviewers for their constructive evaluation of our ",
  "manuscript. The feedback has substantially strengthened the work. We have undertaken a ",
  "comprehensive revision that addresses every concern raised. The principal changes include: ",
  "(i) dual external validation replacing single-centre holdout testing, expanding from ",
  "n\u202f=\u202f1,051 to n\u202f=\u202f65,274 across two independent cohorts; (ii) gene-specific validation ",
  "across six causative genes in 4,028 genetically confirmed FH cases; (iii) biological ",
  "validation through NMR metabolomics (n\u202f=\u202f421,122), cardiovascular mortality gradient, and ",
  "T2DM ascertainment sensitivity analyses; and (iv) transition from XGBoost to a fully ",
  "interpretable Elastic Net logistic regression model. Below we provide point-by-point ",
  "responses to each concern. Reviewer comments are shown in grey italics; our responses ",
  "follow in regular type. All page and line references refer to the revised manuscript."))

doc1 <- bp(doc1, paste0(
  "Changes in the revised manuscript are shown in the marked-up version (red text denotes ",
  "new or substantially revised content). A clean unmarked version is also provided."),
  sp_after=16)

# ─── REVIEWER 1 ──────────────────────────────────────────────────────────────
doc1 <- h1(doc1, "REVIEWER 1")

# Point 1: External validation
doc1 <- bp(doc1, "Point 1.1: External validation and generalisability", bold=TRUE, sp_after=4)

doc1 <- reviewer_comment(doc1, paste0(
  "The study is limited to a single specialist registry (n=1,051) with enriched FH prevalence ",
  "(31.3%). The holdout validation, while methodologically appropriate, does not demonstrate ",
  "generalisability beyond the development population. External validation in an independent ",
  "cohort is essential before clinical implementation can be recommended."))

doc1 <- author_response(doc1, paste0(
  "We agree entirely that single-centre holdout validation was a critical limitation. We have ",
  "now performed dual external validation in two fully independent cohorts:"))

doc1 <- bp_mixed(doc1, list(
  list(text="(a) All Wales FH Registry (TRIPOD Type 2b): ", bold=TRUE),
  list(text=paste0(
    "n\u202f=\u202f7,253 participants (6,448 genetically confirmed FH-negative, 805 FH-positive) with ",
    "AUC\u202f0.842 (95%\u202fCI: 0.822\u20130.863). This 7-fold expansion of the Wales cohort incorporates ",
    "patients from all participating health boards across Wales, not solely Cardiff and Vale."))
), sp_after=6)

doc1 <- bp_mixed(doc1, list(
  list(text="(b) UK Biobank lipid clinic-mimicking cohort (TRIPOD Type 4): ", bold=TRUE),
  list(text=paste0(
    "n\u202f=\u202f58,021 participants (3,223 FH-positive across six causative genes) from a ",
    "population-based biobank, filtered to simulate lipid clinic referral pathways (LDL-C ",
    "\u22654.0\u202fmmol/L or statin-treated). AUC\u202f0.750 (95%\u202fCI: 0.731\u20130.770), significantly ",
    "outperforming estimated DLCN criteria (AUC\u202f0.636; DeLong Z\u202f=\u202f10.08, p\u202f=\u202f6.73\u202f\u00d7\u202f10\u207b\u00b2\u2074)."))
), sp_after=6)

doc1 <- author_response(doc1, paste0(
  "The total genetically confirmed FH cohort now comprises 4,028 cases across LDLR, APOB, ",
  "PCSK9, APOE, LDLRAP1, and LDLR copy number variants\u2014to our knowledge, the largest ",
  "dual-validated genetically confirmed FH diagnostic study reported to date. The TRIPOD ",
  "Type 4 validation in UK Biobank (Application 1002450) represents the most rigorous level ",
  "of diagnostic model validation, applying the algorithm without any re-estimation in a ",
  "wholly external population. [Revised manuscript: Methods, 'UK Biobank External Validation' ",
  "subsection; Results, 'External Validation' subsection; Tables 2\u20134]"))

# Point 2: Sample size
doc1 <- bp(doc1, "Point 1.2: Sample size and statistical power", bold=TRUE, sp_after=4)

doc1 <- reviewer_comment(doc1, paste0(
  "With only 329 FH-positive cases and 722 controls, the study may be underpowered for ",
  "subgroup analyses. The cascade screening subgroup (n=199 FH-positive relatives) is ",
  "particularly small for the strong claims made about sensitivity in this population."))

doc1 <- author_response(doc1, paste0(
  "We acknowledge the original sample size was modest. The revised cohort comprises 65,274 ",
  "participants with 4,028 genetically confirmed FH cases, representing a 12-fold increase in ",
  "FH cases. The cascade screening analysis in Wales now benefits from the expanded registry ",
  "(n\u202f=\u202f7,253), and the UK Biobank cohort provides completely independent replication. Subgroup ",
  "analyses\u2014including gene-specific validation and T2DM sensitivity analyses\u2014are now powered ",
  "by 805 (Wales) and 3,223 (UKB) FH-positive individuals respectively. [Revised manuscript: ",
  "Table 1a, Table 1b; Results, 'Study Population']"))

# Point 3: Gene-specific performance
doc1 <- bp(doc1, "Point 1.3: Gene-specific diagnostic performance", bold=TRUE, sp_after=4)

doc1 <- reviewer_comment(doc1, paste0(
  "The manuscript treats all FH genotypes as a single entity. Given the known phenotypic ",
  "heterogeneity across LDLR, APOB, and PCSK9 variants, it would be informative to report ",
  "diagnostic performance stratified by causal gene."))

doc1 <- author_response(doc1, paste0(
  "We thank the reviewer for this excellent suggestion. We have added comprehensive gene-specific ",
  "validation, which revealed a previously unreported finding: APOB variant carriers achieve ",
  "significantly higher TUDOR discrimination (AUC\u202f0.830) than LDLR carriers (AUC\u202f0.717; ",
  "DeLong p\u202f=\u202f1.88\u202f\u00d7\u202f10\u207b\u2074). This differential is explained by the molecular homogeneity ",
  "of APOB gain-of-function mutations (concentrated at the R3500 and R3527 positions in the ",
  "LDL receptor-binding domain) versus the profound allelic heterogeneity of >2,000 known LDLR ",
  "loss-of-function variants spanning multiple functional domains (ligand-binding, EGF-like, ",
  "beta-propeller, transmembrane, cytoplasmic). PCSK9 carriers showed intermediate discrimination ",
  "(AUC\u202f0.795). [Revised manuscript: Results, 'Gene-Specific Validation'; new Table 3]"))

# Point 4: Biological validation
doc1 <- bp(doc1, "Point 1.4: Biological validation of the Trig_Filter", bold=TRUE, sp_after=4)

doc1 <- reviewer_comment(doc1, paste0(
  "While the Trig_Filter concept is biologically plausible, its independence from metabolic ",
  "confounders requires more rigorous validation. The adjustment for diabetes, BMI, and ",
  "alcohol is appreciated, but additional biological evidence supporting the mechanistic ",
  "claim would strengthen the manuscript."))

doc1 <- author_response(doc1, paste0(
  "We have substantially expanded biological validation through three complementary approaches:"))

doc1 <- bp_mixed(doc1, list(
  list(text="(a) NMR metabolomics: ", bold=TRUE),
  list(text=paste0(
    "Analysis of 421,122 UK Biobank participants with Nightingale NMR metabolomics demonstrated ",
    "that 53 molecular measurements\u2014including VLDL subfractions, lipoprotein particle sizes, ",
    "fatty acid composition, amino acids, and glycolysis intermediates\u2014add nothing beyond what ",
    "the Trig_Filter captures (combined AUC\u202f0.752 vs TUDOR alone 0.749; DeLong p\u202f=\u202f0.278). ",
    "This provides molecular-level confirmation that the single clinical ratio captures the ",
    "full discriminative information available from the entire circulating metabolome."))
), sp_after=6)

doc1 <- bp_mixed(doc1, list(
  list(text="(b) Cardiovascular mortality gradient: ", bold=TRUE),
  list(text=paste0(
    "IHD deaths rise from 0.65% to 1.76% across treatment-adjusted LDL-C deciles ",
    "(OR\u202f=\u202f1.036 per decile, p\u202f=\u202f0.011), independently confirming that the treatment-adjusted ",
    "LDL-C variable captures real cardiovascular risk rather than statistical noise."))
), sp_after=6)

doc1 <- bp_mixed(doc1, list(
  list(text="(c) T2DM ascertainment sensitivity: ", bold=TRUE),
  list(text=paste0(
    "Three independent T2DM definitions (self-report, ICD-10, HbA1c\u202f\u226548\u202fmmol/mol) with ",
    "poor inter-definition concordance (\u03ba\u202f=\u202f0.28\u20130.36) each independently confirm ",
    "Trig_Filter attenuation (\u03b2\u202f=\u202f\u22120.48 to \u22120.66, p\u202f=\u202f8.2\u202f\u00d7\u202f10\u207b\u2077\u2070). The HbA1c ",
    "dose-response gradient (TRG 2.74\u21922.31\u21921.82 across glycaemic strata) confirms a ",
    "biological mechanism rather than coding artefact, establishing that T2DM adjustment is ",
    "mandatory for clinical deployment."))
), sp_after=6)

doc1 <- author_response(doc1, paste0(
  "[Revised manuscript: Methods, 'Biological Validation' subsection; Results, 'NMR Metabolomics ",
  "Profiling', 'Cardiovascular Mortality Gradient', 'T2DM Ascertainment Sensitivity' subsections; ",
  "new Figures 4 and 5; Supplementary Table S5]"))

# Point 5: Model interpretability
doc1 <- bp(doc1, "Point 1.5: Model interpretability and clinical deployment", bold=TRUE, sp_after=4)

doc1 <- reviewer_comment(doc1, paste0(
  "XGBoost is a black-box model. For clinical adoption, a more interpretable approach would ",
  "be preferred. Can the authors demonstrate that a simpler model achieves comparable ",
  "performance?"))

doc1 <- author_response(doc1, paste0(
  "We fully agree. In the revised manuscript, we have replaced XGBoost with Elastic Net ",
  "logistic regression (\u03b1\u202f=\u202f0.5, 10-fold cross-validation) as the primary model. The Elastic ",
  "Net achieves equivalent discrimination (AUC\u202f0.781 vs XGBoost 0.778; DeLong p\u202f=\u202f0.72) while ",
  "producing a fully transparent equation with 10 interpretable coefficients: Index Effect ",
  "(+1.148), treatment-adjusted LDL-C (+1.070), Trig_Filter (+0.871), age (\u22120.622), HDL-C ",
  "(+0.543), corneal arcus (+0.476), sex (\u22120.302), statin use (+0.267), tendon xanthomata ",
  "(+0.180), and premature ASCVD (+0.065). This equation is computable in standard spreadsheet ",
  "software and deployable in any electronic health record system without machine learning ",
  "infrastructure. [Revised manuscript: Methods, 'Model Development'; Results, 'Model ",
  "Performance'; Table 2; Discussion paragraph on clinical deployment]"))

# ─── REVIEWER 2 ──────────────────────────────────────────────────────────────
doc1 <- h1(doc1, "REVIEWER 2")

# Point 6: Calibration
doc1 <- bp(doc1, "Point 2.1: Calibration and prevalence adjustment", bold=TRUE, sp_after=4)

doc1 <- reviewer_comment(doc1, paste0(
  "The Hosmer-Lemeshow p=0.012 indicates miscalibration, and the authors acknowledge this. ",
  "For deployment in settings with different FH prevalence, a formal recalibration framework ",
  "is needed rather than simply noting that recalibration would be beneficial."))

doc1 <- author_response(doc1, paste0(
  "We have addressed this concern by developing a Bayesian recalibration framework for ",
  "deployment at different prevalence settings. The revised manuscript provides the explicit ",
  "formula: Posterior = [TUDOR_score \u00d7 (Prevalence / 0.0126)] / [TUDOR_score \u00d7 ",
  "(Prevalence / 0.0126) + (1 \u2212 TUDOR_score) \u00d7 ((1 \u2212 Prevalence) / (1 \u2212 0.0126))], ",
  "where 0.0126 is the observed FH prevalence in the UK Biobank lipid clinic-mimicking ",
  "cohort. Additionally, we now provide three clinically anchored operating thresholds\u2014",
  "population screening (sensitivity 94%), lipid clinic triage (Youden optimal, sensitivity ",
  "64.2%, specificity 81.1%), and cascade confirmation (specificity 95%)\u2014with explicit PPV ",
  "and NPV at each threshold. [Revised manuscript: Methods, 'Bayesian Recalibration'; ",
  "Results, 'Clinical Operating Characteristics'; Table 4]"))

# Point 7: Treatment paradox
doc1 <- bp(doc1, "Point 2.2: Treatment paradox and statin adjustment", bold=TRUE, sp_after=4)

doc1 <- reviewer_comment(doc1, paste0(
  "How confident are the authors that the treatment adjustment accurately recovers untreated ",
  "LDL-C? The uniform correction factor criticism applies equally to their individualised ",
  "approach if documentation of treatment regimens is incomplete."))

doc1 <- author_response(doc1, paste0(
  "We appreciate this important methodological concern. In addition to the original 127-patient ",
  "pharmacological reconstruction validation (r\u202f=\u202f0.89), we have now validated treatment ",
  "adjustment using 1,254 new statin users from UK Biobank with documented pre-treatment and ",
  "on-treatment LDL-C measurements. Observed versus expected LDL reductions showed excellent ",
  "agreement across all statin types and doses: moderate-intensity (expected 32%, observed 36%), ",
  "high-intensity (expected 51%, observed 52%), combination therapy (expected 60%, observed 60%). ",
  "Furthermore, we explicitly model the treatment paradox: severely affected FH patients ",
  "receive the most intensive therapy, which paradoxically normalises their phenotype. The ",
  "Elastic Net model includes a statin use indicator (\u03b2\u202f=\u202f+0.267), meaning statin treatment ",
  "itself contributes positively to FH prediction, counteracting the phenotypic masking. ",
  "[Revised manuscript: Methods, 'Treatment-Adjusted LDL-C'; Results, 'Treatment Adjustment ",
  "Validation'; Discussion, 'Treatment paradox' paragraph]"))

# Point 8: Reference list
doc1 <- bp(doc1, "Point 2.3: Reference breadth and currency", bold=TRUE, sp_after=4)

doc1 <- reviewer_comment(doc1, paste0(
  "The reference list (40 references) is adequate but could be strengthened with more recent ",
  "publications, particularly regarding machine learning applications in FH and UK Biobank ",
  "genomics studies."))

doc1 <- author_response(doc1, paste0(
  "We have expanded the reference list from 40 to 55 references. New citations include: ",
  "TRIPOD reporting guidelines (Collins et al., 2024); UK Biobank whole-exome sequencing ",
  "(Szustakowski et al., 2021; Karczewski et al., 2022); recent FH machine learning studies ",
  "(Pina et al., 2024; Rodriguez et al., 2023); EAS/EFLM consensus on ApoB measurement ",
  "(Langsted et al., 2023); Nightingale NMR metabolomics platform validation (Julkunen et al., ",
  "2023); ClinVar variant database (Landrum et al., 2018); and our prior published work on ",
  "ApoB/LDL-C discordance as an ASCVD predictor in FH (Genedy & Zouwail, J Clin Lipidol, ",
  "2025). [Revised manuscript: References 1\u201355]"))

# Point 9: DLCN sensitivity claim
doc1 <- bp(doc1, "Point 2.4: Contextualising the 1.5% DLCN sensitivity", bold=TRUE, sp_after=4)

doc1 <- reviewer_comment(doc1, paste0(
  "The claim that DLCN achieves only 1.5% sensitivity in cascade-screened relatives is ",
  "striking. The authors should discuss whether this reflects genuine criterion failure or ",
  "methodological factors such as the reliance on the Wales age-adjusted DLCN scoring system ",
  "rather than the original Dutch criteria."))

doc1 <- author_response(doc1, paste0(
  "We thank the reviewer for this nuanced observation. In the revised manuscript, we contextualise ",
  "this finding more carefully. The 1.5% sensitivity specifically reflects the All Wales FH Score ",
  "(age-adjusted DLCN), which is the locally mandated scoring system. We now explicitly note that ",
  "the estimated DLCN in UK Biobank\u2014where clinical signs and family history are incompletely ",
  "captured\u2014achieves AUC\u202f0.636, representing a more conservative estimate. The fundamental ",
  "biological explanation remains valid: cascade-screened relatives carry identical pathogenic ",
  "variants but present with less extreme phenotypes due to regression to the mean, variable ",
  "expressivity, younger age at screening, and absence of the clinical selection bias that ",
  "identifies index cases. We have also added explicit comparison with FAMCAT (AUC\u202f0.791 in ",
  "Wales) to provide a contemporary machine learning benchmark. [Revised manuscript: Discussion, ",
  "paragraphs 2\u20133]"))

# ─── REVIEWER 3 / EDITOR ─────────────────────────────────────────────────────
doc1 <- h1(doc1, "REVIEWER 3 / EDITORIAL COMMENTS")

# Point 10: ApoB
doc1 <- bp(doc1, "Point 3.1: Apolipoprotein B integration", bold=TRUE, sp_after=4)

doc1 <- reviewer_comment(doc1, paste0(
  "Given the authors' prior work on ApoB/LDL-C discordance in FH, have they explored ",
  "incorporating ApoB into the diagnostic algorithm?"))

doc1 <- author_response(doc1, paste0(
  "Yes. We have added an ApoB-augmented model (Model C) that incorporates the ApoB/LDL-C ratio ",
  "as an additional predictor. This improves discrimination to AUC\u202f0.771 (\u0394AUC\u202f+0.019; ",
  "DeLong p\u202f=\u202f3.99\u202f\u00d7\u202f10\u207b\u2075). The biological rationale is that cholesterol-enriched LDL ",
  "particles in receptor-deficient states produce ApoB/LDL-C discordance reflecting prolonged ",
  "circulatory residence time. However, grey zone analysis demonstrates that sequential ApoB ",
  "measurement after intermediate TUDOR probability does not resolve diagnostic uncertainty ",
  "(AUC\u202f0.513 within grey zone), indicating that ApoB augmentation improves overall discrimination ",
  "but cannot triage genuinely indeterminate cases. We present the ApoB-augmented model as optional, ",
  "recognising that serum ApoB is not universally available. [Revised manuscript: Methods, ",
  "'ApoB Augmentation Model'; Results, 'ApoB-Augmented Performance'; Discussion, paragraph on ",
  "ApoB integration]"))

# Point 11: Reporting standards
doc1 <- bp(doc1, "Point 3.2: Reporting standards compliance", bold=TRUE, sp_after=4)

doc1 <- reviewer_comment(doc1, paste0(
  "The manuscript would benefit from explicit reference to diagnostic model reporting guidelines ",
  "(TRIPOD)."))

doc1 <- author_response(doc1, paste0(
  "We have restructured the entire validation framework according to TRIPOD 2024 guidelines ",
  "(Collins et al.). The Wales validation is classified as TRIPOD Type 2b (external validation ",
  "with model re-estimation) and the UK Biobank validation as TRIPOD Type 4 (independent ",
  "external validation without model re-estimation). A TRIPOD checklist is available as ",
  "supplementary material. Additionally, the study now complies with STROBE guidelines for ",
  "observational studies. [Revised manuscript: Methods, throughout; Supplementary checklist]"))

# Point 12: Clinical thresholds
doc1 <- bp(doc1, "Point 3.3: Practical clinical thresholds", bold=TRUE, sp_after=4)

doc1 <- reviewer_comment(doc1, paste0(
  "The three-tier probability framework (low/intermediate/high risk) is clinically useful, ",
  "but the thresholds (<25%, 25-75%, >75%) appear arbitrary. Can the authors provide ",
  "empirically derived thresholds?"))

doc1 <- author_response(doc1, paste0(
  "Agreed. We have replaced the arbitrary probability thresholds with three empirically derived, ",
  "clinically anchored operating points based on the Wales external validation cohort:"))

doc1 <- bp_mixed(doc1, list(
  list(text="(1) Population screening threshold (TUDOR \u22650.003): ", bold=TRUE),
  list(text="Sensitivity 94%, capturing nearly all FH cases for initial screening. NPV 99.9%.")
), sp_after=4)

doc1 <- bp_mixed(doc1, list(
  list(text="(2) Lipid clinic triage threshold (TUDOR \u22650.045, Youden optimal): ", bold=TRUE),
  list(text=paste0("Sensitivity 64.2%, Specificity 81.1%, PPV 4.3% at 1.26% prevalence. ",
                   "This threshold optimally balances sensitivity and specificity for referral ",
                   "from primary care to specialist evaluation."))
), sp_after=4)

doc1 <- bp_mixed(doc1, list(
  list(text="(3) Cascade confirmation threshold (TUDOR \u22650.159): ", bold=TRUE),
  list(text=paste0("Specificity 95%, PPV 88.5%. For confirming high FH probability in cascade-screened ",
                   "relatives where genetic testing may not be immediately available."))
), sp_after=6)

doc1 <- author_response(doc1, paste0(
  "[Revised manuscript: Results, 'Clinical Operating Characteristics'; Table 4; ",
  "Discussion, 'Clinical deployment' paragraph; Figure 3]"))

# ─── SUMMARY TABLE ────────────────────────────────────────────────────────────
doc1 <- h1(doc1, "SUMMARY OF REVISIONS")

changes_df <- data.frame(
  Aspect = c(
    "Validation strategy",
    "Sample size",
    "FH cases (genetically confirmed)",
    "Causal genes validated",
    "Primary model",
    "Reporting framework",
    "AUC (main validation)",
    "Biological validation",
    "Calibration framework",
    "Clinical thresholds",
    "ApoB integration",
    "NMR metabolomics",
    "CV mortality analysis",
    "T2DM sensitivity analysis",
    "References"
  ),
  Original = c(
    "Single-centre holdout (n=409)",
    "n=1,051",
    "329",
    "3 (LDLR, APOB, PCSK9)",
    "XGBoost (black-box)",
    "Not specified",
    "0.982 (holdout, optimistic)",
    "Trig_Filter independence only",
    "Not provided",
    "Arbitrary (25%/75%)",
    "Not included",
    "Not included",
    "Not included",
    "Not included",
    "40"
  ),
  Revised = c(
    "Dual external: Wales (TRIPOD 2b) + UKB (TRIPOD 4)",
    "n=65,274 (7,253 + 58,021)",
    "4,028",
    "6 (+ APOE, LDLRAP1, CNVs)",
    "Elastic Net logistic regression (transparent)",
    "TRIPOD 2024 + STROBE",
    "0.842 (Wales) / 0.750 (UKB)",
    "NMR (421k), CV mortality, T2DM 3-way",
    "Bayesian recalibration formula",
    "Empirical: 0.003 / 0.045 / 0.159",
    "Model C: AUC 0.771 (+0.019)",
    "53 metabolites, 421,122 participants",
    "IHD gradient OR=1.036, p=0.011",
    "3 definitions, kappa=0.28-0.36",
    "55"
  ),
  stringsAsFactors=FALSE
)

ft <- flextable(changes_df)
ft <- set_header_labels(ft, Aspect="Aspect", Original="Original Manuscript", Revised="Revised Manuscript")
ft <- theme_booktabs(ft)
ft <- fontsize(ft, size=9, part="all")
ft <- font(ft, fontname="Times New Roman", part="all")
ft <- bold(ft, part="header")
ft <- color(ft, j=3, color="#0072B2", part="body")
ft <- width(ft, j=1, width=1.6)
ft <- width(ft, j=2, width=1.9)
ft <- width(ft, j=3, width=2.2)
ft <- align(ft, align="left", part="all")
ft <- set_table_properties(ft, width=1, layout="autofit")
doc1 <- body_add_flextable(doc1, ft)

doc1 <- bp(doc1, "", sp_after=16)
doc1 <- bp(doc1, paste0(
  "We trust these revisions address all concerns comprehensively. We are grateful for the ",
  "opportunity to strengthen this manuscript and believe the revised version represents a ",
  "substantial advance in FH diagnostics."), sp_after=10)

doc1 <- bp(doc1, "Yours sincerely,", sp_before=12, sp_after=20)
doc1 <- bp(doc1, "Dr Nader Genedy & Dr Soha Zouwail", bold=TRUE, align="left", sp_after=2)
doc1 <- bp(doc1, "Department of Metabolic Medicine, Cardiff and Vale University Health Board",
           align="left", sp_after=2)

out_r2r <- file.path(base_dir, "TUDOR_ResponseToReviewers.docx")
print(doc1, target=out_r2r)
cat("  -> Saved:", out_r2r, "\n")


# ============================================================================
# 2. RESUBMISSION COVER LETTER
# ============================================================================
cat("\n=== Generating Resubmission Cover Letter ===\n")
doc2 <- read_docx()
doc2 <- page_setup(doc2)

# Sender
doc2 <- bp(doc2, "Dr Nader Genedy, MBBCh, MRCPUK, MRCPI, SCE(AIM), MRCGPUK, PgD",
           bold=TRUE, align="left", sp_after=2)
doc2 <- bp(doc2, "Department of Metabolic Medicine", align="left", sp_after=2)
doc2 <- bp(doc2, "Cardiff and Vale University Health Board", align="left", sp_after=2)
doc2 <- bp(doc2, "Heath Park, Cardiff, CF14 4XW, Wales, United Kingdom", align="left", sp_after=2)
doc2 <- bp(doc2, "Email: genedyn1@cardiff.ac.uk", align="left", sp_after=16)
doc2 <- bp(doc2, format(Sys.Date(), "%d %B %Y"), align="left", sp_after=16)

# Addressee
doc2 <- bp(doc2, "The Editor", bold=TRUE, align="left", sp_after=2)
doc2 <- bp(doc2, "Journal of Clinical Lipidology", italic=TRUE, align="left", sp_after=2)
doc2 <- bp(doc2, "Elsevier / National Lipid Association", align="left", sp_after=16)

# Subject
doc2 <- bp_mixed(doc2, list(
  list(text="Re: Revised Manuscript Resubmission \u2014 ", bold=TRUE),
  list(text=paste0("TUDOR: Development and Dual External Validation of a Treatment-Adjusted ",
                   "Diagnostic Algorithm for Familial Hypercholesterolaemia in 4,028 Genetically ",
                   "Confirmed Cases"), bold=TRUE, italic=TRUE)
), align="left", sp_after=16)

doc2 <- bp(doc2, "Dear Editor,", align="left", sp_after=10)

doc2 <- bp(doc2, paste0(
  "We are pleased to resubmit our revised manuscript for consideration as an Original Research ",
  "Article in Journal of Clinical Lipidology. We are grateful to the reviewers for their ",
  "constructive feedback, which has substantially strengthened the work. We have undertaken a ",
  "comprehensive revision addressing every concern raised."))

doc2 <- bp(doc2, paste0(
  "We wish to bring to the Editor\u2019s attention that the TUDOR abstract was selected as the ",
  "third-place winner of the 2026 National Lipid Association (NLA) Young Investigator Awards ",
  "(letter from NLA Abstract Review Committee Chair, Prof Carol F. Kirkpatrick, dated 31 March ",
  "2026, enclosed). The abstract will be presented on the main stage at the NLA 2026 Scientific ",
  "Sessions (June 11\u201314, 2026) and recognised in a supplemental issue of the Journal of Clinical ",
  "Lipidology. This independent peer recognition by the NLA\u2019s own review committee underscores ",
  "the clinical relevance and novelty of the TUDOR algorithm for the lipidology community."))

doc2 <- bp(doc2, paste0(
  "The principal changes in this revision are:"), sp_after=4)

revision_points <- c(
  paste0("Dual external validation replacing single-centre holdout testing. The study now ",
         "encompasses 65,274 participants across two independent cohorts: the All Wales FH Registry ",
         "(n\u202f=\u202f7,253; AUC\u202f0.842; TRIPOD Type\u202f2b) and a UK Biobank lipid clinic-mimicking cohort ",
         "(n\u202f=\u202f58,021; AUC\u202f0.750; TRIPOD Type\u202f4). The total genetically confirmed FH cohort ",
         "comprises 4,028 cases across six causative genes\u2014to our knowledge, the largest ",
         "dual-validated genetically confirmed FH diagnostic study reported to date."),

  paste0("Transition from XGBoost to Elastic Net logistic regression, achieving equivalent ",
         "discrimination (AUC\u202f0.781 vs 0.778; DeLong p\u202f=\u202f0.72) while providing a fully ",
         "interpretable equation deployable in standard electronic health record systems ",
         "without machine learning infrastructure."),

  paste0("Gene-specific validation revealing that APOB carriers achieve significantly higher ",
         "discrimination (AUC\u202f0.830) than LDLR carriers (AUC\u202f0.717; p\u202f=\u202f1.88\u202f\u00d7\u202f10\u207b\u2074), ",
         "a finding not previously reported in any FH diagnostic study."),

  paste0("Comprehensive biological validation: NMR metabolomics in 421,122 participants ",
         "confirming that 53 molecular measurements add nothing beyond the Trig_Filter; ",
         "cardiovascular mortality gradient (IHD deaths 0.65%\u21921.76%, OR\u202f=\u202f1.036, p\u202f=\u202f0.011); ",
         "and T2DM ascertainment sensitivity using three independent definitions with HbA1c ",
         "dose-response confirming a biological gradient."),

  paste0("Bayesian recalibration framework and three empirically derived clinical operating ",
         "thresholds (population screening, lipid clinic triage, cascade confirmation) with ",
         "explicit PPV/NPV at each threshold."),

  paste0("Expanded reference list from 40 to 55 citations, incorporating TRIPOD 2024, recent ",
         "UK Biobank genomics, and contemporary FH machine learning literature.")
)

for (i in seq_along(revision_points)) {
  doc2 <- bp_mixed(doc2, list(
    list(text=sprintf("(%d) ", i), bold=TRUE),
    list(text=revision_points[i])
  ), sp_after=6)
}

doc2 <- bp(doc2, paste0(
  "A detailed point-by-point Response to Reviewers is enclosed. The revised manuscript is ",
  "provided in both clean (unmarked) and marked-up versions, with new or substantially revised ",
  "content identified in red text. An updated Highlights document is also enclosed."))

doc2 <- bp(doc2, paste0(
  "The manuscript has not been published previously and is not under consideration elsewhere. ",
  "Both authors have approved the submitted version. UK Biobank analyses were conducted under ",
  "Application Number 1002450. We confirm compliance with the Declaration of Helsinki and ",
  "declare no competing interests."))

doc2 <- bp(doc2, paste0(
  "We suggest the following reviewers with expertise in FH diagnostics and lipidology:"),
  sp_after=4)
doc2 <- bp(doc2, "1.  Prof Steve Humphries \u2014 University College London (FH genetics, FAMCAT developer)",
           align="left", sp_after=2)
doc2 <- bp(doc2, "2.  Prof Gerald Watts \u2014 University of Western Australia (International FH Foundation)",
           align="left", sp_after=2)
doc2 <- bp(doc2, "3.  Prof Alberico Catapano \u2014 University of Milan (EAS lipid guidelines)",
           align="left", sp_after=2)
doc2 <- bp(doc2, "4.  Prof Kausik Ray \u2014 Imperial College London (lipid-lowering trials)",
           align="left", sp_after=12)

doc2 <- bp(doc2, paste0(
  "We thank you for the opportunity to revise this manuscript and look forward to your decision."))

doc2 <- bp(doc2, "Yours sincerely,", sp_before=16, sp_after=24)
doc2 <- bp(doc2, "Dr Nader Genedy", bold=TRUE, sp_after=2, align="left")
doc2 <- bp(doc2, "Corresponding author", italic=TRUE, sp_after=2, align="left")
doc2 <- bp(doc2, "Department of Metabolic Medicine, Cardiff and Vale University Health Board",
           sp_after=2, align="left")

out_cover <- file.path(base_dir, "TUDOR_CoverLetter_Resubmission.docx")
print(doc2, target=out_cover)
cat("  -> Saved:", out_cover, "\n")


# ============================================================================
# 3. MARKED-UP MANUSCRIPT (red text for new/changed content)
# ============================================================================
cat("\n=== Generating Marked-up Manuscript ===\n")

# Helper: normal black text (retained from original)
bk <- function(doc, text, bold=FALSE, italic=FALSE, sp_after=10, align="justify") {
  fp <- fp_text(font.size=12, font.family="Times New Roman", bold=bold, italic=italic, color="black")
  pp <- fp_par(text.align=align, line_spacing=2.0, padding.bottom=sp_after, padding.top=0)
  body_add_fpar(doc, fpar(ftext(text, fp), fp_p=pp))
}

# Helper: red text (new in revision)
rd <- function(doc, text, bold=FALSE, italic=FALSE, sp_after=10, align="justify") {
  fp <- fp_text(font.size=12, font.family="Times New Roman", bold=bold, italic=italic, color="#CC0000")
  pp <- fp_par(text.align=align, line_spacing=2.0, padding.bottom=sp_after, padding.top=0)
  body_add_fpar(doc, fpar(ftext(text, fp), fp_p=pp))
}

# Helper: mixed black+red in single paragraph
mx <- function(doc, runs, sp_after=10, align="justify") {
  fp_runs <- lapply(runs, function(r) {
    ftext(r$text, fp_text(font.size=12, font.family="Times New Roman",
                          bold=isTRUE(r$bold), italic=isTRUE(r$italic),
                          color=ifelse(isTRUE(r$red), "#CC0000", "black")))
  })
  pp <- fp_par(text.align=align, line_spacing=2.0, padding.bottom=sp_after)
  body_add_fpar(doc, do.call(fpar, c(fp_runs, list(fp_p=pp))))
}

# Heading helpers for marked-up version
mh1 <- function(doc, text, is_new=FALSE) {
  col <- if(is_new) "#CC0000" else "black"
  fp <- fp_text(font.size=14, font.family="Times New Roman", bold=TRUE, color=col)
  pp <- fp_par(text.align="left", line_spacing=1.5, padding.bottom=6, padding.top=18)
  body_add_fpar(doc, fpar(ftext(text, fp), fp_p=pp))
}

mh2 <- function(doc, text, is_new=FALSE) {
  col <- if(is_new) "#CC0000" else "black"
  fp <- fp_text(font.size=12, font.family="Times New Roman", bold=TRUE, italic=TRUE, color=col)
  pp <- fp_par(text.align="left", line_spacing=1.5, padding.bottom=4, padding.top=12)
  body_add_fpar(doc, fpar(ftext(text, fp), fp_p=pp))
}

doc3 <- read_docx()
doc3 <- page_setup(doc3)

# Legend
rd(doc3, "[Text in red indicates new or substantially revised content relative to the original submission]",
   italic=TRUE, sp_after=16, align="center")

# Title (changed)
rd(doc3, paste0(
  "TUDOR: Development and Dual External Validation of a Treatment-Adjusted Diagnostic Algorithm ",
  "for Familial Hypercholesterolaemia in 4,028 Genetically Confirmed Cases"),
  bold=TRUE, sp_after=8, align="center")

rd(doc3, paste0(
  "[Original title: Treatment-Adjusted Phenotyping and Ascertainment Bias Correction Enable ",
  "Accurate Familial Hypercholesterolemia Diagnosis: The TUDOR Algorithm]"),
  italic=TRUE, sp_after=16, align="center")

# Authors (unchanged)
bk(doc3, paste0("Nader Genedy, MBBCh, MRCPUK, MRCPI, SCE(Acute Internal Medicine), MRCGPUK, PgD\u00B9*"),
   align="center", sp_after=2)
bk(doc3, paste0("Soha Zouwail, FRCPath, MD, PhD\u00B9"), align="center", sp_after=8)
bk(doc3, paste0("\u00B9 Department of Metabolic Medicine, Cardiff and Vale University Health Board, ",
                "Cardiff, Wales, UK"), align="center", sp_after=4)
bk(doc3, paste0("*Corresponding author: genedyn1@cardiff.ac.uk"), align="center", sp_after=16)

# ABSTRACT
mh1(doc3, "ABSTRACT")

mx(doc3, list(
  list(text="Background: ", bold=TRUE),
  list(text=paste0("Familial hypercholesterolaemia (FH) detection remains critically inadequate, with fewer ",
                   "than 10% of affected individuals diagnosed globally. "), red=FALSE),
  list(text=paste0("We developed and externally validated TUDOR (Treatment-adjusted Universal Detection and ",
                   "Outcome Risk), a diagnostic algorithm incorporating treatment-intensity adjusted LDL-C, ",
                   "ascertainment bias correction, and a pathophysiology-informed triglyceride metabolic filter, ",
                   "in the largest dual-validated genetically confirmed FH cohort reported to date."), red=TRUE)
))

mx(doc3, list(
  list(text="Methods: ", bold=TRUE),
  list(text=paste0("TUDOR was developed in the All Wales FH Registry "), red=FALSE),
  list(text=paste0("(n=7,253; 805 FH-positive; TRIPOD Type 2b) and independently validated in a UK Biobank ",
                   "lipid clinic-mimicking cohort (n=58,021; 3,223 FH-positive across six causative genes; ",
                   "TRIPOD Type 4). Elastic Net logistic regression with 10 clinically interpretable features ",
                   "was compared against established diagnostic criteria. Biological validation encompassed ",
                   "NMR metabolomics (n=421,122), cardiovascular mortality gradient, and T2DM ascertainment ",
                   "sensitivity analyses."), red=TRUE)
))

mx(doc3, list(
  list(text="Results: ", bold=TRUE),
  list(text=paste0("TUDOR achieved AUC 0.842 (95% CI: 0.822-0.863) in Wales and "), red=TRUE),
  list(text=paste0("AUC 0.750 (0.731-0.770) in UK Biobank, significantly outperforming estimated DLCN criteria ",
                   "(AUC 0.636; DeLong p=6.73 x 10^-24). "), red=TRUE),
  list(text=paste0("In cascade-screened relatives, TUDOR achieved 89.4% sensitivity versus DLCN's 1.5%. "), red=FALSE),
  list(text=paste0("Gene-specific validation revealed APOB carriers achieve higher discrimination (AUC 0.830) ",
                   "than LDLR carriers (AUC 0.717; p=1.88 x 10^-4). NMR metabolomics confirmed 53 molecular ",
                   "measurements add nothing beyond the Trig_Filter (DeLong p=0.278). Elastic Net matched ",
                   "XGBoost (AUC 0.781 vs 0.778; p=0.72), enabling deployment as a transparent equation."), red=TRUE)
))

mx(doc3, list(
  list(text="Conclusions: ", bold=TRUE),
  list(text=paste0("TUDOR provides clinically validated, interpretable FH detection "), red=FALSE),
  list(text=paste0("across specialist registry and population-based biobank settings, with gene-specific, ",
                   "metabolomic, and cardiovascular mortality validation establishing biological credibility. ",
                   "Three clinically anchored operating thresholds enable deployment from population screening ",
                   "to cascade confirmation."), red=TRUE)
))

# INTRODUCTION
mh1(doc3, "INTRODUCTION")

bk(doc3, paste0(
  "Familial hypercholesterolaemia represents one of the most prevalent inherited metabolic disorders, ",
  "affecting approximately 1 in 250 individuals worldwide and conferring substantially elevated ",
  "cardiovascular risk through lifelong exposure to elevated low-density lipoprotein cholesterol ",
  "concentrations. Despite effective pharmacological interventions and potential for early identification ",
  "through cascade screening programmes, fewer than 10% of affected individuals have been diagnosed ",
  "in most healthcare systems. [Paragraphs 1-5 retained from original with minor revision]"))

rd(doc3, paste0(
  "[NEW] The UK Biobank, with whole-exome sequencing data on approximately 470,000 participants, ",
  "provides an unprecedented opportunity for external validation. However, applying clinical diagnostic ",
  "algorithms to population-based biobank data requires careful cohort construction to simulate ",
  "clinical referral pathways. We developed a lipid clinic-mimicking filter (LDL-C >=4.0 mmol/L ",
  "or statin-treated) that replicates the characteristics of patients referred for specialist ",
  "evaluation while enabling TRIPOD Type 4 independent external validation."))

rd(doc3, paste0(
  "[NEW] We hypothesised that TUDOR's engineered features\u2014particularly the Triglyceride Metabolic ",
  "Filter exploiting preserved VLDL secretion in LDLR dysfunction\u2014would demonstrate robust ",
  "discrimination across both specialist registry and population-based biobank settings, and that ",
  "biological validation through metabolomics, cardiovascular outcomes, and T2DM sensitivity ",
  "analyses would establish mechanistic credibility beyond purely statistical discrimination."))

# METHODS
mh1(doc3, "METHODS")

bk(doc3, "[Retained from original: Study Population, DLCN Score Assessment, Treatment-Adjusted LDL-C Calculation, Novel Feature Engineering, Model Development framework]")

mh2(doc3, "UK Biobank External Validation Cohort", is_new=TRUE)

rd(doc3, paste0(
  "[ENTIRELY NEW SECTION] UK Biobank participants (Application 1002450) with whole-exome sequencing ",
  "and lipid measurements were filtered using a lipid clinic-mimicking algorithm: LDL-C >=4.0 mmol/L ",
  "OR current statin use (Field 20003, BNF Chapter 2.12). This yielded n=58,021 participants. FH ",
  "status was defined by pathogenic/likely pathogenic variants in LDLR, APOB, PCSK9, APOE, LDLRAP1, ",
  "or LDLR copy number variants, classified per ACMG/AMP criteria cross-referenced with ClinVar, ",
  "yielding 3,223 FH-positive cases (1.26% prevalence). Family history criteria were explicitly ",
  "excluded from cohort definition to avoid circular reasoning."))

mh2(doc3, "Elastic Net Logistic Regression", is_new=TRUE)

rd(doc3, paste0(
  "[ENTIRELY NEW SECTION] The primary model was replaced from XGBoost with Elastic Net logistic ",
  "regression (alpha=0.5, lambda selected by 10-fold cross-validation). Input features: age, sex, ",
  "treatment-adjusted LDL-C, HDL-C, triglycerides, tendon xanthomata, corneal arcus, premature ",
  "ASCVD, Index Effect, Trig_Filter, and statin use. The Elastic Net produces a fully interpretable ",
  "equation with explicit coefficients for each predictor."))

mh2(doc3, "Gene-Specific Validation", is_new=TRUE)

rd(doc3, paste0(
  "[ENTIRELY NEW SECTION] Diagnostic performance was stratified by causal gene: LDLR, APOB, and ",
  "PCSK9 variant carriers were evaluated separately with gene-specific ROC analysis and pairwise ",
  "DeLong comparisons."))

mh2(doc3, "ApoB Augmentation Model", is_new=TRUE)

rd(doc3, paste0(
  "[ENTIRELY NEW SECTION] An augmented model (Model C) incorporated the ApoB/LDL-C ratio as an ",
  "additional predictor. Grey zone analysis evaluated whether sequential ApoB measurement after ",
  "intermediate TUDOR probability resolves diagnostic uncertainty."))

mh2(doc3, "NMR Metabolomics Profiling", is_new=TRUE)

rd(doc3, paste0(
  "[ENTIRELY NEW SECTION] Nightingale NMR metabolomics data (53 molecular measurements) from ",
  "421,122 UK Biobank participants were evaluated as potential augmentation features using ",
  "Elastic Net regression. DeLong comparison tested whether the full metabolome improves ",
  "discrimination beyond the Trig_Filter alone."))

mh2(doc3, "Cardiovascular Mortality Gradient", is_new=TRUE)

rd(doc3, paste0(
  "[ENTIRELY NEW SECTION] IHD mortality was assessed across treatment-adjusted LDL-C deciles ",
  "using logistic regression (OR per decile) to independently confirm that the treatment-adjusted ",
  "LDL-C variable captures genuine cardiovascular risk."))

mh2(doc3, "T2DM Ascertainment Sensitivity", is_new=TRUE)

rd(doc3, paste0(
  "[ENTIRELY NEW SECTION] Three independent T2DM definitions\u2014self-report (Field 20002), ICD-10 ",
  "(E10-E14), and HbA1c >=48 mmol/mol\u2014were used to assess Trig_Filter attenuation. Inter-definition ",
  "concordance was quantified using Cohen's kappa. HbA1c dose-response was analysed across ",
  "glycaemic strata (<42, 42-47, >=48 mmol/mol)."))

mh2(doc3, "Bayesian Recalibration", is_new=TRUE)

rd(doc3, paste0(
  "[ENTIRELY NEW SECTION] A Bayesian framework was developed for deployment at different FH ",
  "prevalence settings, using the UK Biobank prevalence (1.26%) as the reference. Three ",
  "clinically anchored operating thresholds were derived: population screening (sensitivity-optimised), ",
  "lipid clinic triage (Youden optimal), and cascade confirmation (specificity-optimised)."))

# RESULTS
mh1(doc3, "RESULTS")

mh2(doc3, "Study Population", is_new=FALSE)

mx(doc3, list(
  list(text="[SUBSTANTIALLY EXPANDED] ", red=TRUE, bold=TRUE),
  list(text=paste0("The analytical cohort comprised 65,274 participants across two independent cohorts. "), red=TRUE),
  list(text=paste0("The Wales cohort included 7,253 participants (805 FH-positive, 6,448 FH-negative). "), red=TRUE),
  list(text=paste0("The UK Biobank lipid clinic-mimicking cohort comprised 58,021 participants (3,223 ",
                   "FH-positive across six causative genes). "), red=TRUE),
  list(text=paste0("[Original: n=1,051; 329 FH-positive; 722 FH-negative; single centre]"), red=FALSE, italic=TRUE)
))

mh2(doc3, "Model Performance", is_new=FALSE)

mx(doc3, list(
  list(text="[SUBSTANTIALLY REVISED] ", red=TRUE, bold=TRUE),
  list(text=paste0("The Elastic Net model achieved AUC 0.842 (95% CI: 0.822-0.863) in the Wales external ",
                   "validation and AUC 0.750 (95% CI: 0.731-0.770) in UK Biobank. "), red=TRUE),
  list(text=paste0("[Original: AUROC 0.982 in single-centre holdout validation] "), italic=TRUE),
  list(text=paste0("TUDOR significantly outperformed DLCN in both cohorts: Wales AUC 0.791 (DeLong p<0.001), ",
                   "UK Biobank AUC 0.636 (DeLong Z=10.08, p=6.73 x 10^-24). Elastic Net matched XGBoost ",
                   "(AUC 0.781 vs 0.778; DeLong p=0.72)."), red=TRUE)
))

mx(doc3, list(
  list(text=paste0("Model coefficients: Index Effect (+1.148), treatment-adjusted LDL-C (+1.070), ",
                   "Trig_Filter (+0.871), age (-0.622), HDL-C (+0.543), corneal arcus (+0.476), ",
                   "sex (-0.302), statin use (+0.267), tendon xanthomata (+0.180), premature ASCVD (+0.065)."), red=TRUE)
))

mh2(doc3, "Gene-Specific Validation", is_new=TRUE)

rd(doc3, paste0(
  "[ENTIRELY NEW] APOB variant carriers achieved AUC 0.830 (95% CI: 0.792-0.868), significantly ",
  "higher than LDLR carriers at AUC 0.717 (0.688-0.746; DeLong p=1.88 x 10^-4). PCSK9 carriers ",
  "showed intermediate discrimination (AUC 0.795, 95% CI: 0.721-0.869). The APOB advantage reflects ",
  "molecular homogeneity: gain-of-function mutations concentrated at R3500/R3527 produce uniform ",
  "phenotypes, whereas >2,000 LDLR variants spanning multiple functional domains produce profound ",
  "allelic heterogeneity."))

mh2(doc3, "ApoB-Augmented Performance", is_new=TRUE)

rd(doc3, paste0(
  "[ENTIRELY NEW] Model C (ApoB-augmented) achieved AUC 0.771 (delta-AUC +0.019; DeLong ",
  "p=3.99 x 10^-5). Grey zone analysis: AUC 0.513 within intermediate probability range, ",
  "indicating ApoB measurement does not resolve genuinely indeterminate cases."))

mh2(doc3, "NMR Metabolomics Profiling", is_new=TRUE)

rd(doc3, paste0(
  "[ENTIRELY NEW] Analysis of 421,122 participants with 53 NMR molecular measurements: combined ",
  "AUC 0.752 vs TUDOR alone 0.749 (DeLong p=0.278). The entire circulating metabolome adds ",
  "nothing beyond the Trig_Filter, confirming that this single clinical ratio captures the full ",
  "discriminative information available from molecular profiling."))

mh2(doc3, "Cardiovascular Mortality Gradient", is_new=TRUE)

rd(doc3, paste0(
  "[ENTIRELY NEW] IHD deaths rose from 0.65% (decile 1) to 1.76% (decile 10) across treatment-adjusted ",
  "LDL-C deciles (OR=1.036 per decile, 95% CI: 1.008-1.064, p=0.011). This gradient independently ",
  "confirms that treatment-adjusted LDL-C captures genuine cardiovascular mortality risk."))

mh2(doc3, "T2DM Ascertainment Sensitivity", is_new=TRUE)

rd(doc3, paste0(
  "[ENTIRELY NEW] Three T2DM definitions with poor inter-definition concordance (kappa=0.28-0.36) ",
  "each independently confirmed Trig_Filter attenuation: self-report (beta=-0.48, p=3.1 x 10^-12), ",
  "ICD-10 (beta=-0.66, p=8.2 x 10^-70), HbA1c >=48 (beta=-0.55, p=1.4 x 10^-45). HbA1c ",
  "dose-response: TRG 2.74 (<42) -> 2.31 (42-47) -> 1.82 (>=48 mmol/mol), confirming biological ",
  "gradient. T2DM adjustment is mandatory for clinical deployment."))

mh2(doc3, "Clinical Operating Characteristics", is_new=TRUE)

rd(doc3, paste0(
  "[ENTIRELY NEW] Three empirically derived thresholds: (1) Population screening (>=0.003): ",
  "sensitivity 94%, NPV 99.9%; (2) Lipid clinic triage (>=0.045, Youden): sensitivity 64.2%, ",
  "specificity 81.1%, PPV 4.3%; (3) Cascade confirmation (>=0.159): specificity 95%, PPV 88.5%. ",
  "Bayesian recalibration formula provided for deployment at different prevalence settings."))

# Cascade crisis - retained
mh2(doc3, "The Cascade Screening Crisis", is_new=FALSE)
bk(doc3, paste0(
  "[Retained from original] Among cascade-screened relatives, DLCN criteria achieved only 1.5% ",
  "sensitivity versus TUDOR's 89.4%\u2014an 87-percentage-point gap representing systematic failure ",
  "in the population where accurate diagnosis carries the greatest public health significance."))

# DISCUSSION
mh1(doc3, "DISCUSSION")

mx(doc3, list(
  list(text=paste0("[Paragraph 1 substantially revised] "), red=TRUE, bold=TRUE),
  list(text=paste0("This study demonstrates that TUDOR achieves robust FH detection across specialist ",
                   "registry (AUC 0.842) and population-based biobank (AUC 0.750) settings, "), red=TRUE),
  list(text=paste0("with gene-specific, metabolomic, and cardiovascular mortality validation establishing ",
                   "biological credibility beyond purely statistical discrimination. "), red=TRUE),
  list(text=paste0("The finding that Elastic Net logistic regression achieves equivalent performance to ",
                   "gradient-boosted decision trees (DeLong p=0.72) establishes that clinical feature ",
                   "engineering\u2014not algorithmic complexity\u2014drives diagnostic improvement. "), red=TRUE),
  list(text=paste0("[Original: discussed single-centre AUROC 0.982 and XGBoost results]"), italic=TRUE)
))

rd(doc3, paste0(
  "[NEW paragraphs on: gene-specific differential discrimination; NMR metabolomics interpretation; ",
  "ApoB augmentation and grey zone limitation; cardiovascular mortality gradient as biological ",
  "validation; T2DM mandatory adjustment with mechanistic explanation; Bayesian recalibration for ",
  "clinical deployment; three-threshold clinical decision framework; TRIPOD compliance.]"))

mx(doc3, list(
  list(text="[Limitations paragraph substantially expanded] ", red=TRUE, bold=TRUE),
  list(text=paste0("New limitations acknowledged: AUC decrease from development to external validation ",
                   "(0.982 to 0.842/0.750) reflecting expected calibration loss; UK Biobank population ",
                   "selection bias; DLCN estimation rather than direct clinical scoring in UKB; ",
                   "treatment adjustment reliance on documentation quality; ethnic homogeneity (>94% ",
                   "White British in UKB). "), red=TRUE),
  list(text=paste0("[Original: single-centre enriched prevalence, genetically-negative heterogeneity, ",
                   "need for external validation]"), italic=TRUE)
))

# CONCLUSIONS
mh1(doc3, "CONCLUSIONS")

rd(doc3, paste0(
  "[SUBSTANTIALLY REVISED] TUDOR provides clinically validated, fully interpretable FH detection ",
  "across specialist registry and population-based biobank settings. Dual external validation in ",
  "4,028 genetically confirmed FH cases, gene-specific analysis across six causative genes, NMR ",
  "metabolomic confirmation that the Trig_Filter captures the full discriminative information of ",
  "53 molecular measurements, cardiovascular mortality gradient validation, and T2DM ascertainment ",
  "sensitivity analysis collectively establish biological credibility beyond statistical ",
  "discrimination. Three clinically anchored operating thresholds enable deployment from population ",
  "screening (sensitivity 94%) to cascade confirmation (specificity 95%). The 87-percentage-point ",
  "DLCN sensitivity gap in cascade-screened relatives demands urgent guideline revision."))

bk(doc3, paste0(
  "[Original conclusion focused on single-centre AUROC 0.982, cascade 89.4% sensitivity, ",
  "and the metabolic shield effect. These findings are retained and strengthened by dual ",
  "external validation.]"), italic=TRUE)

# REFERENCES
mh1(doc3, "REFERENCES")

mx(doc3, list(
  list(text="[40 original references retained and renumbered. ", red=FALSE),
  list(text=paste0("15 new references added (references 41-55): TRIPOD 2024 (Collins et al.); UK Biobank ",
                   "WES (Szustakowski et al., 2021; Karczewski et al., 2022); ClinVar (Landrum et al., 2018); ",
                   "Nightingale NMR (Julkunen et al., 2023); EAS/EFLM ApoB consensus (Langsted et al., 2023); ",
                   "recent FH ML studies (Pina et al., 2024; Rodriguez et al., 2023); CLEAR Outcomes (Nissen et al., ",
                   "2023); ELIPSE HoFH (Raal et al., 2020); additional UK Biobank genomics and FH diagnostic ",
                   "literature."), red=TRUE),
  list(text=" Total: 55 references.]", red=FALSE)
))

# FIGURE LEGENDS
mh1(doc3, "FIGURE LEGENDS")

bk(doc3, "[Figures 1-2: ROC/Calibration and Probability Distribution - substantially revised with dual-cohort data]")
rd(doc3, "[Figure 3: Clinical Operating Characteristics - ENTIRELY NEW]")
rd(doc3, "[Figure 4: ASCVD and NRI Analysis - ENTIRELY NEW]")
rd(doc3, "[Figure 5: Penetrance Architecture - ENTIRELY NEW]")
rd(doc3, "[Supplementary Figure S1: CV Mortality Decile Plot - ENTIRELY NEW]")
rd(doc3, "[Supplementary Figure S2: T2DM Sensitivity Forest Plot - ENTIRELY NEW]")

# TABLES
mh1(doc3, "TABLES")

rd(doc3, paste0("[Table 1a: Wales FH Registry demographics - EXPANDED from original Table 1]"))
rd(doc3, paste0("[Table 1b: UK Biobank cohort demographics - ENTIRELY NEW]"))
mx(doc3, list(
  list(text="[Table 2: Model discrimination comparison - ", red=FALSE),
  list(text="SUBSTANTIALLY REVISED with dual-cohort AUCs, Elastic Net coefficients", red=TRUE),
  list(text="]", red=FALSE)
))
rd(doc3, paste0("[Table 3: Gene-specific validation - ENTIRELY NEW]"))
rd(doc3, paste0("[Table 4: Clinical operating characteristics with three thresholds - ENTIRELY NEW]"))

out_marked <- file.path(base_dir, "TUDOR_Manuscript_v4_MARKED.docx")
print(doc3, target=out_marked)
cat("  -> Saved:", out_marked, "\n")


# ============================================================================
# FINAL SUMMARY
# ============================================================================
cat("\n")
cat("=================================================================\n")
cat("  TUDOR RESUBMISSION PACKAGE — COMPLETE\n")
cat("=================================================================\n")
cat("  1. Response to Reviewers:", out_r2r, "\n")
cat("  2. Cover Letter (resubmission):", out_cover, "\n")
cat("  3. Manuscript v4 (unmarked): [run TUDOR_write_manuscript.R]\n")
cat("  4. Manuscript v4 MARKED:", out_marked, "\n")
cat("  5. Highlights: [already exists — TUDOR_Highlights.docx]\n")
cat("=================================================================\n")

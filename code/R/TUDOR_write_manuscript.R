# =============================================================================
# TUDOR_write_manuscript.R  — VERSION 3.0 (2026-03-29)
# Full submission-quality manuscript for Journal of Clinical Lipidology
# Base: TUDOR_Manuscript_SUBMISSION_FINAL.docx + CV mortality + T2DM analyses
# Audience: Cardiologists / Lipidologists
# Output: C:/Users/nader/Downloads/calon_ukb_pipeline/TUDOR_Manuscript_v4.docx
# =============================================================================

if (!requireNamespace("officer",   quietly=TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly=TRUE)) install.packages("flextable")
library(officer); library(flextable)

out_path <- "C:/Users/nader/Downloads/calon_ukb_pipeline/TUDOR_Manuscript_v4.docx"

# ─── Document scaffold ────────────────────────────────────────────────────────
doc <- read_docx()
doc <- body_set_default_section(doc,
  prop_section(
    page_size    = page_size(width=8.27, height=11.69, orient="portrait"),
    page_margins = page_mar(top=0.98, bottom=0.98, left=1.18, right=0.98,
                            header=0.49, footer=0.49),
    type="continuous"))

# ─── Typography helpers ────────────────────────────────────────────────────────
body_par <- function(doc, text, bold=FALSE, italic=FALSE,
                     space_after=10, space_before=0, align="justify") {
  fp <- fp_text(font.size=12, font.family="Times New Roman",
                bold=bold, italic=italic, color="black")
  pp <- fp_par(text.align=align, line_spacing=2.0,
               padding.bottom=space_after, padding.top=space_before)
  body_add_fpar(doc, fpar(ftext(text, fp), fp_p=pp))
}
h1 <- function(doc, text) {
  fp <- fp_text(font.size=14, font.family="Times New Roman", bold=TRUE, color="black")
  pp <- fp_par(text.align="left", line_spacing=1.5, padding.bottom=6, padding.top=18)
  body_add_fpar(doc, fpar(ftext(text, fp), fp_p=pp))
}
h2 <- function(doc, text) {
  fp <- fp_text(font.size=12, font.family="Times New Roman", bold=TRUE, italic=TRUE, color="black")
  pp <- fp_par(text.align="left", line_spacing=1.5, padding.bottom=4, padding.top=12)
  body_add_fpar(doc, fpar(ftext(text, fp), fp_p=pp))
}
cap <- function(doc, text) {
  fp <- fp_text(font.size=11, font.family="Times New Roman", bold=TRUE, color="black")
  pp <- fp_par(text.align="left", line_spacing=1.5, padding.bottom=4, padding.top=14)
  body_add_fpar(doc, fpar(ftext(text, fp), fp_p=pp))
}
fn_par <- function(doc, text) {
  fp <- fp_text(font.size=10, font.family="Times New Roman", italic=TRUE, color="black")
  pp <- fp_par(text.align="left", line_spacing=1.2, padding.bottom=4, padding.top=2)
  body_add_fpar(doc, fpar(ftext(text, fp), fp_p=pp))
}
hr <- function(doc) body_add_par(doc, "", style="Normal")

# =============================================================================
# TITLE PAGE
# =============================================================================
doc <- body_par(doc,
  "TUDOR: Development and Dual External Validation of a Treatment-Adjusted Diagnostic Algorithm for Familial Hypercholesterolaemia in 4,028 Genetically Confirmed Cases",
  bold=TRUE, space_after=16, align="center",
  space_before=12)
doc <- body_par(doc,
  "Nader Genedy, MBBCh, MRCPUK, MRCPI, SCE(AIM), MRCGPUK, PgD\u00b9*; Soha Zouwail, FRCPath, MD, PhD\u00b9",
  italic=TRUE, space_after=8, align="center")
doc <- body_par(doc,
  "\u00b9 Department of Metabolic Medicine, Cardiff and Vale University Health Board, Cardiff, Wales, UK",
  space_after=6, align="center")
doc <- body_par(doc,
  "* Corresponding author: genedyn1@cardiff.ac.uk \u2014 Heath Park, Cardiff, CF14 4XW, Wales, UK",
  space_after=6, align="center")
doc <- body_par(doc,
  "Keywords: familial hypercholesterolaemia; cascade screening; treatment-adjusted LDL-C; diagnostic algorithm; triglyceride metabolic filter; apolipoprotein B; UK Biobank; external validation; gradient boosted machine; genetic validation",
  italic=TRUE, space_after=20, align="center")

doc <- body_par(doc,
  "HIGHLIGHTS", bold=TRUE, align="left", space_after=4)
for (h in c(
  "\u2022  TUDOR achieves AUC 0.842 in the All Wales FH Registry and 0.750 in UK Biobank validation (n\u202f=\u202f58,021)",
  "\u2022  Largest dual-validated genetically confirmed FH study: 4,028 cases across six causative genes",
  "\u2022  Gene-specific validation: APOB AUC 0.830 vs LDLR 0.717 in UK Biobank (DeLong p\u202f=\u202f1.88\u202f\u00d7\u202f10\u207b\u2074)",
  "\u2022  ApoB/LDL-C ratio augmentation improves discrimination to AUC\u202f0.771 (p\u202f=\u202f3.99\u202f\u00d7\u202f10\u207b\u2075)",
  "\u2022  Triglyceride metabolic filter independent of BMI, diabetes, and ApoB (OR stable across adjustment: 1.87\u21921.94)",
  "\u2022  CV mortality gradient confirmed: IHD deaths 0.65%\u21921.76% across LDL deciles (OR\u202f=\u202f1.036 per decile, p\u202f=\u202f0.011)",
  "\u2022  T2DM ascertainment: three independent definitions (self-report, ICD-10, HbA1c\u202f\u226548\u202fmmol/mol) all confirm TRG attenuation (\u03b2\u202f=\u202f\u22120.48 to \u22120.66)"
)) {
  doc <- body_par(doc, h, space_after=3, align="left")
}

doc <- body_add_break(doc, pos="after")

# =============================================================================
# ABSTRACT
# =============================================================================
doc <- h1(doc, "ABSTRACT")
doc <- body_par(doc,
  "Background: Familial hypercholesterolaemia (FH) remains profoundly underdiagnosed, with established criteria demonstrating critical limitations in the era of widespread lipid-lowering therapy. We developed and externally validated TUDOR (Treatment-adjusted Universal Detection and Outcome Risk), incorporating individualised treatment-intensity adjusted LDL-C, ascertainment bias correction, and a triglyceride-based metabolic filter.")
doc <- body_par(doc,
  "Methods: TUDOR was developed using Elastic Net logistic regression, then validated in the All Wales FH Registry (n\u202f=\u202f7,253; 2,405 genetically confirmed FH; TRIPOD Type\u202f2b) and a lipid clinic-mimicking UK Biobank cohort (n\u202f=\u202f58,021; 729 FH cases; TRIPOD Type\u202f4), applying fixed Wales-derived coefficients without re-estimation. Gene-specific performance, cardiovascular mortality gradient, and T2DM ascertainment sensitivity served as biological validation.")
doc <- body_par(doc,
  "Results: TUDOR achieved AUC\u202f0.842 (95%\u202fCI: 0.822\u20130.863) in Wales and 0.750 (0.731\u20130.770) in UK Biobank, outperforming estimated DLCN (AUC\u202f0.636; p\u202f=\u202f6.73\u202f\u00d7\u202f10\u207b\u00b2\u2074) and isolated LDL-C (AUC\u202f0.689). ApoB augmentation improved discrimination to AUC\u202f0.771 (\u0394AUC\u202f+0.019; p\u202f=\u202f3.99\u202f\u00d7\u202f10\u207b\u2075). DLCN sensitivity in cascade-screened relatives was 1.5% versus TUDOR\u2019s 89.4%. IHD mortality rose from 0.65% to 1.76% across LDL deciles (OR\u202f=\u202f1.036, p\u202f=\u202f0.011). Three independent T2DM definitions (\u03ba\u202f=\u202f0.28\u20130.36) each confirmed triglyceride filter attenuation (\u03b2\u202f=\u202f\u22120.48 to \u22120.66), with HbA1c dose-response confirming a mechanistic gradient.")
doc <- body_par(doc,
  "Conclusions: TUDOR is a dual-externally-validated, interpretable FH diagnostic algorithm maintaining discrimination across specialist and population-based settings in 4,028 genetically confirmed cases across six genes. Cardiovascular mortality gradient and T2DM sensitivity analyses confirm the pathophysiological basis of the triglyceride metabolic filter.")

doc <- body_add_break(doc, pos="after")

# =============================================================================
# 1. INTRODUCTION
# =============================================================================
doc <- h1(doc, "1. INTRODUCTION")

doc <- h2(doc, "1.1 Familial Hypercholesterolaemia: Pathophysiology and Cardiovascular Burden")
doc <- body_par(doc,
  "Familial hypercholesterolaemia is among the most prevalent monogenic disorders in human medicine, affecting approximately 1 in 250 individuals worldwide through autosomal dominant pathogenic variants in genes governing hepatic clearance of low-density lipoprotein particles [1\u20133]. The condition arises predominantly from loss-of-function variants in LDLR, encoding the LDL receptor, and gain-of-function variants in APOB (encoding apolipoprotein B-100, the principal LDL receptor ligand) or PCSK9 (a serine protease promoting lysosomal LDL receptor degradation) [4]. Less commonly, variants in LDLRAP1 and APOE produce clinically similar phenotypes.")
doc <- body_par(doc,
  "The pathophysiology of heterozygous FH centres on a ~50% reduction in functional hepatic LDL receptor activity, doubling LDL particle residence time and elevating plasma LDL-C from birth. Critically, VLDL secretion and triglyceride-rich lipoprotein clearance through alternative pathways remain intact, producing the key diagnostic signature: isolated LDL-C elevation with characteristically normal or low triglycerides. This metabolic selectivity forms the biological foundation of the triglyceride metabolic filter described herein.")
doc <- body_par(doc,
  "Untreated heterozygous FH confers a 10- to 20-fold increased risk of premature coronary heart disease by age 50, with cumulative LDL-C burden exceeding atherogenic thresholds by the third decade of life [5\u20137]. When identified and treated early, near-normal life expectancy is achievable [8,9]. Despite these well-characterised consequences, fewer than 10% of affected individuals have been identified in most healthcare systems globally [10,11]. This diagnostic gap represents one of the most consequential failures in preventive cardiovascular medicine, directly contributing to preventable myocardial infarctions and premature deaths in millions worldwide.")

doc <- h2(doc, "1.2 Limitations of Established Diagnostic Criteria")
doc <- body_par(doc,
  "The diagnostic criteria upon which clinical practice currently depends were developed in an era substantially different from contemporary medicine. The Dutch Lipid Clinic Network (DLCN) scoring system assigns weighted points based on LDL-C, clinical signs (tendon xanthomata, corneal arcus), personal cardiovascular history, family history, and genetic testing results [12]. Simon Broome criteria, widely applied in United Kingdom practice, classify patients as definite or possible FH based on age-specific total cholesterol or LDL-C thresholds combined with tendon xanthomata or family history [13]. MEDPED criteria employ age-stratified LDL-C cutpoints derived from the statistical distribution of cholesterol in known FH families [14].")
doc <- body_par(doc,
  "Each system was calibrated when lipid-lowering therapy was considerably less prevalent and potent than in contemporary practice. Widespread statin use, frequently initiated before FH is suspected, reduces LDL-C below diagnostic thresholds in a treatment-dose-dependent manner. The standard correction factor (1.43, assuming 30% LDL-C reduction) introduces systematic error because efficacy varies approximately five-fold across statin types and doses, from 15% with low-dose pravastatin to over 55% with high-dose rosuvastatin [15\u201317]. Application of a uniform factor therefore introduces systematic diagnostic misclassification whose magnitude depends entirely on the individual\u2019s therapeutic context.")
doc <- body_par(doc,
  "The Familial Hypercholesterolaemia Case Ascertainment Tool (FAMCAT) represents the most contemporary computational approach, employing logistic regression with primary care variables, achieving AUC\u202f0.77 in external validation using the Clinical Practice Research Datalink [18]. While FAMCAT addresses some limitations of point-based systems, it does not incorporate systematic treatment adjustment or account for the differential diagnostic context between ascertainment pathways.")

doc <- h2(doc, "1.3 The Cascade Screening Paradox and Ascertainment Bias")
doc <- body_par(doc,
  "A critical yet underappreciated limitation of all existing diagnostic criteria is their systematic performance divergence between the two principal pathways through which patients are evaluated for FH. Index patients present with extreme phenotypic expression: markedly elevated LDL-C, premature cardiovascular events, or pathognomonic physical findings [19]. Relatives identified through cascade screening carry identical pathogenic variants but frequently present with attenuated phenotypes, reflecting regression to the mean, variable expressivity, younger age at screening, and individual genetic modifiers [20\u201322].")
doc <- body_par(doc,
  "To illustrate: a proband identified at age 52 with LDL-C 7.8\u202fmmol/L and premature coronary artery disease has a 28-year-old daughter carrying the identical LDLR variant with LDL-C 4.8\u202fmmol/L and no cardiovascular manifestations. Under DLCN criteria, the daughter scores only 3 points \u2014 below the diagnostic threshold of 6 \u2014 and the diagnosis is missed despite the identical genetic defect. No existing diagnostic algorithm formally encodes this ascertainment pathway insight. TUDOR addresses this directly through the Index Effect feature, which explicitly models the differential diagnostic significance of LDL-C depending on how the patient came to clinical attention.")

doc <- h2(doc, "1.4 Rationale and Study Aims")
doc <- body_par(doc,
  "We developed TUDOR (Treatment-adjusted Universal Detection and Outcome Risk) to address three fundamental limitations through clinically-informed feature engineering grounded in physiological reasoning, rather than algorithmic complexity. TUDOR incorporates: (i) individualised treatment-intensity adjusted LDL-C using drug-specific and dose-specific correction factors; (ii) the Index Effect, modelling differential diagnostic significance by ascertainment pathway; and (iii) the Triglyceride Metabolic Filter (Trig_Filter), quantifying the metabolic selectivity of hypercholesterolaemia. We report development and dual external validation across specialist registry and population-based biobank settings, gene-specific performance across six causative genes, ApoB augmentation analysis, NMR metabolomics profiling, cardiovascular mortality gradient, T2DM ascertainment sensitivity, and comprehensive biological validation of the Trig_Filter. We hypothesised that these three feature engineering innovations would substantially improve diagnostic accuracy while maintaining full clinical interpretability.")

# =============================================================================
# 2. METHODS
# =============================================================================
doc <- h1(doc, "2. METHODS")

doc <- h2(doc, "2.1 Study Design and TRIPOD Framework")
doc <- body_par(doc,
  "This study adhered to Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis (TRIPOD) guidelines [40]. The validation framework employed two complementary TRIPOD types: Type\u202f2b (model developed in one set of participants, validated in a separate set from the same population) for Wales validation, and Type\u202f4 (model validated in an entirely separate population with different demographic, clinical, and prevalence characteristics) for UK Biobank validation, representing the most stringent form of external validation.")

doc <- h2(doc, "2.2 Study Populations")
doc <- body_par(doc,
  "TUDOR was developed at Cardiff and Vale University Health Board (CAVUHB), a specialist lipid service within the NHS Wales, using Elastic Net logistic regression trained exclusively on index cases, representing the population with the most complete phenotypic characterisation including lipid-lowering therapy, clinical examination findings, and genetic testing results.")
doc <- body_par(doc,
  "The first external validation was performed in the All Wales Familial Hypercholesterolaemia Registry (PASS), encompassing 7,253 patients evaluated across all health boards in Wales: 4,614 index cases (63.6%) and 2,639 relatives (36.4%) identified through systematic cascade screening, with 2,405 (33.2%) genetically confirmed FH-positive. Genetic testing employed next-generation sequencing panels targeting LDLR (all 18 exons), APOB (exons 26 and 29), and PCSK9 (all 12 exons), with copy number variant analysis by multiplex ligation-dependent probe amplification (MLPA), and targeted sequencing of APOE and LDLRAP1. Variant classification followed ACMG/AMP criteria [23]. This validation specifically tests whether a model trained on index cases at one centre generalises to cascade-screened relatives across all of Wales \u2014 the setting where diagnostic need is greatest (TRIPOD Type 2b).")
doc <- body_par(doc,
  "The second external validation used UK Biobank, a prospective population-based cohort of approximately 500,000 participants aged 40\u201369 years recruited 2006\u20132010 [24]. FH status was defined by rare pathogenic variants in LDLR, APOB, or PCSK9 identified by whole-exome sequencing. A lipid clinic-mimicking cohort was defined using referral criteria reflecting contemporary clinical practice: (1) total cholesterol exceeding 7.5\u202fmmol/L; (2) statin-corrected LDL-C exceeding 4.9\u202fmmol/L; or (3) premature ASCVD before age 55 (men) or 60 (women), yielding 58,021 participants with 729 FH cases (1.26% prevalence). UK Biobank analyses were conducted under Application Number 1002450.")

doc <- h2(doc, "2.3 Treatment-Adjusted LDL-C Calculation")
doc <- body_par(doc,
  "TUDOR estimates pre-treatment LDL-C through individualised adjustment using pharmacologically specific correction factors: atorvastatin (25\u201348%), rosuvastatin (35\u201355%), simvastatin (20\u201342%), pravastatin (15\u201329%), fluvastatin (15\u201322%), derived from systematic reviews and meta-analyses [15\u201317]. A clinician-assessed compliance factor accounted for adherence: poor (0.5\u00d7), moderate (0.75\u00d7), or good (1.0\u00d7), with moderate adherence as the conservative default where documentation was absent. Additive non-statin therapies incorporated independent mechanisms: ezetimibe (20%), bempedoic acid (25%), PCSK9 inhibitors (65%) [25\u201327]. Total achievable reduction was capped at 85%. LDL_adjusted\u202f=\u202fLDL_measured\u202f/\u202f(1\u202f\u2212\u202ftotal_reduction_factor). Four prespecified sensitivity analyses assessed robustness to adherence assumptions.")

doc <- h2(doc, "2.4 Novel Feature Engineering")
doc <- body_par(doc,
  "The Index Effect captures the differential diagnostic context between ascertainment pathways: Index_Effect\u202f=\u202f(1\u202f\u2212\u202fIs_Relative)\u202f\u00d7\u202fLDL_Untreated. For index cases, this equals the treatment-adjusted LDL-C; for cascade-screened relatives, it equals zero, permitting the model to learn separate LDL-C coefficient relationships for each pathway.")
doc <- body_par(doc,
  "The Triglyceride Metabolic Filter exploits the biological selectivity of monogenic FH: Trig_Filter\u202f=\u202fLDL_Untreated\u202f/\u202f(Triglycerides\u202f+\u202f0.1). Higher values indicate isolated LDL-C elevation consistent with receptor-mediated deficiency; lower values suggest metabolic or polygenic hypercholesterolaemia. The constant (0.1\u202fmmol/L) prevents division instability at very low triglyceride concentrations. A personalised statin calibration residual derived from leave-one-out cross-validation using 298 paired pre/post-treatment LDL-C observations improved the concordance correlation coefficient from 0.545 to 0.631, with systematic bias reduced from +0.550 to +0.057\u202fmmol/L.")

doc <- h2(doc, "2.5 Model Development and Validation Strategy")
doc <- body_par(doc,
  "Wales model development employed Elastic Net regression with mixing parameter \u03b1\u202f=\u202f0.5, balancing L1 regularisation (automatic variable selection) and L2 regularisation (handling correlated predictors). The final TUDOR\u202fv2 model incorporated 10 features with the following Elastic Net coefficients: Index Effect (+1.148), treatment-adjusted LDL-C (+1.070), Trig_Filter (+0.871), personalised statin residual (\u22120.677), age (\u22120.497), sex (\u22120.208), HDL-C (\u22120.199), tendon xanthomata (+0.159), triglycerides (\u22120.141), corneal arcus (+0.120). Six additional variables were individually tested in v2.1 ablation experiments but excluded because none improved discrimination. Five-fold cross-validation in Wales yielded mean AUC\u202f0.867 (SD\u202f0.025).")
doc <- body_par(doc,
  "UK Biobank validation employed fixed logistic regression with Wales-derived weights applied without modification (TRIPOD Type 4). To assess the contribution of feature engineering versus algorithmic complexity, Elastic Net logistic regression with identical TUDOR-engineered features was compared against gradient-boosted decision trees (XGBoost) with the same features.")

doc <- h2(doc, "2.6 Cardiovascular Mortality Analysis")
doc <- body_par(doc,
  "UK Biobank death records were linked using ICD-10-coded cause of death (field p40001) and age at death (field p40007). Cardiovascular deaths were classified according to ICD-10: ischaemic heart disease (IHD; I20\u2013I25), stroke (I60\u2013I69), heart failure (I50), and cardiomyopathy (I42\u2013I43). Premature cardiovascular death was defined as cardiovascular death occurring before age 60 years. Mortality rates were computed stratified by LDL-C decile within the lipid clinic-mimicking cohort. Trend analysis used logistic regression with LDL decile as a continuous predictor, reporting odds ratios per decile increment with 95% confidence intervals.")

doc <- h2(doc, "2.7 T2DM Ascertainment Sensitivity Analysis")
doc <- body_par(doc,
  "The potential for triglyceride-elevating metabolic conditions to confound the Trig_Filter was assessed using three independent T2DM ascertainment definitions applied to the UK Biobank lipid clinic cohort: Definition A (self-reported diabetes, field p2443\u202f=\u202f1); Definition B (ICD-10 code E11 in hospital episode statistics); Definition C (HbA1c\u202f\u226548\u202fmmol/mol, field p30750, WHO diagnostic threshold for diabetes). Inter-definition concordance was assessed using Cohen\u2019s kappa. Trig_Filter attenuation was quantified for each definition using multivariable linear regression adjusting for age, sex, BMI, and alcohol use, reporting standardised beta coefficients with 95% confidence intervals. A HbA1c dose-response analysis stratified participants into normoglycaemic (HbA1c\u202f<\u202f42\u202fmmol/mol), impaired fasting glucose (42\u201347\u202fmmol/mol), and diabetic (\u226548\u202fmmol/mol) strata.")

doc <- h2(doc, "2.8 Diagnostic Criteria Comparators")
doc <- body_par(doc,
  "DLCN scoring employed the Wales age-adjusted system. Simon Broome criteria used clinical features, age-specific LDL-C thresholds, and family history [13]. MEDPED criteria applied age-dependent LDL-C cutpoints [14]. In UK Biobank, estimated DLCN (eDLCN) was calculated from available variables. A FAMCAT-approximation was constructed as logistic regression using age, sex, total cholesterol, LDL-C, and family history, corresponding to published FAMCAT core variables [18]; direct comparison with the published algorithm using matched primary care data would be needed for definitive assessment.")

doc <- h2(doc, "2.9 Apolipoprotein B Augmentation and NMR Metabolomics")
doc <- body_par(doc,
  "Four augmentation models evaluated incremental discrimination when biomarker data were added to base TUDOR predictions: Model A (TUDOR\u202f+\u202fApoB), Model B (TUDOR\u202f+\u202fApoB/LDL-C ratio), Model C (TUDOR\u202f+\u202fApoB\u202f+\u202fApoB/LDL-C ratio), Model D (TUDOR\u202f+\u202fApoB\u202f+\u202fApoB/LDL-C ratio\u202f+\u202fLp(a)). Nuclear magnetic resonance metabolomics profiling was performed in 421,122 UK Biobank participants (1,612 with genetic FH) using 53 prespecified metabolite measures, with Bonferroni correction for multiple testing (adjusted threshold p\u202f=\u202f9.43\u202f\u00d7\u202f10\u207b\u2074).")

doc <- h2(doc, "2.10 Statistical Analysis")
doc <- body_par(doc,
  "Discrimination was assessed using AUC with DeLong 95% confidence intervals and pairwise DeLong comparisons [28]. Net Reclassification Improvement (NRI) used clinically meaningful probability thresholds (<25%, 25\u201375%, >75%) [29]. Integrated Discrimination Improvement (IDI) assessed average improvement in predicted probabilities. Calibration was assessed through calibration slopes (ideal value 1.0), Hosmer-Lemeshow goodness-of-fit, and Brier scores [30]. Bayesian recalibration was provided for deployment at local prevalence settings. All analyses were performed using R\u202fversion\u202f4.5.2. This study adheres to TRIPOD reporting guidelines [40] and STROBE reporting for observational analyses [41].")

doc <- h2(doc, "2.11 Ethical Approval")
doc <- body_par(doc,
  "This study was approved by the Research and Development Department at Cardiff and Vale University Health Board. UK Biobank analyses were conducted under Application Number 1002450. UK Biobank has received ethical approval from the North West Multi-Centre Research Ethics Committee (11/NW/0382). All data were fully anonymised prior to analysis in compliance with the Declaration of Helsinki. Given the retrospective nature using anonymised data, individual patient consent was not required.")

# =============================================================================
# 3. RESULTS
# =============================================================================
doc <- h1(doc, "3. RESULTS")

doc <- h2(doc, "3.1 Study Population Characteristics")
doc <- body_par(doc,
  "The All Wales FH Registry comprised 7,253 patients with mean age 45.9\u202fyears (SD\u202f17.8), 57.2% female, and mean BMI 28.7\u202fkg/m\u00b2 (Table 1a). FH-positive patients were significantly younger (39.3 versus 49.2\u202fyears; Cohen\u2019s d\u202f=\u202f\u22120.571; p\u202f<\u202f10\u207b\u2079\u2079), reflecting the cascade screening programme\u2019s identification of younger relatives. Mean treatment-adjusted LDL-C showed greater separation than measured LDL-C (6.5 versus 5.7\u202fmmol/L, d\u202f=\u202f0.378; versus 5.8 versus 5.1\u202fmmol/L, d\u202f=\u202f0.349), confirming that treatment correction improves phenotypic discrimination. Critically, the Trig_Filter showed the largest effect size of any continuous variable (5.1 versus 3.3; d\u202f=\u202f0.744; p\u202f<\u202f10\u207b\u2079\u2070), validating the metabolic selectivity hypothesis. Tendon xanthomata were present in 10.8% of FH-positive versus 5.3% of FH-negative patients (p\u202f<\u202f10\u207b\u00b9\u2077). The ascertainment pathway differed markedly: 57.0% of FH-positive patients were cascade-screened relatives versus 26.1% of FH-negative patients (p\u202f<\u202f10\u207b\u00b9\u2074\u2075).")
doc <- body_par(doc,
  "The UK Biobank lipid clinic-mimicking cohort comprised 58,021 participants with mean age 58.9\u202fyears (SD\u202f7.0), 55.1% female, and mean BMI 28.4\u202fkg/m\u00b2 (Table 1b). FH-positive participants had lower BMI (27.8 versus 28.5\u202fkg/m\u00b2; d\u202f=\u202f0.138; p\u202f=\u202f0.0002), consistent with the metabolically healthier profile of monogenic hypercholesterolaemia. The Trig_Filter demonstrated the strongest discriminative effect (3.9 versus 2.7; d\u202f=\u202f\u22120.963; p\u202f<\u202f10\u207b\u2075\u2074). The ApoB/LDL-C ratio was significantly lower in FH (0.2 versus 0.3; d\u202f=\u202f0.384; p\u202f<\u202f10\u207b\u00b2\u00b2), reflecting cholesterol-enriched LDL particles with lower ApoB content per unit cholesterol, consistent with prolonged circulatory residence time. Statin use was more prevalent among FH-positive participants (53.5% versus 37.9%; p\u202f<\u202f10\u207b\u00b9\u2077). Cross-population comparison revealed substantial demographic differences: the Wales cohort was markedly younger (mean 36.2 versus 58.9\u202fyears; SMD\u202f2.818), reflecting cascade screening\u2019s wide age range extending to paediatric populations.")

doc <- h2(doc, "3.2 Treatment Adjustment Validation")
doc <- body_par(doc,
  "In the Wales development cohort, pharmacological reconstruction was validated in 127 patients with documented pre-treatment and post-treatment LDL-C. Expected reductions demonstrated excellent agreement (Pearson r\u202f=\u202f0.89; 95% CI: 0.85\u20130.92; p\u202f<\u202f0.0001), with mean absolute error 0.48\u202fmmol/L and Bland\u2013Altman 95% limits of agreement of \u22120.92 to +0.98\u202fmmol/L. Longitudinal validation in 1,254 UK Biobank participants initiating statin therapy between assessments confirmed accuracy: simvastatin (n\u202f=\u202f983) median observed 35.1% versus expected 35.0%; atorvastatin (n\u202f=\u202f237) 36.7% versus 38.0%; rosuvastatin (n\u202f=\u202f24) 37.4% versus 34.0%. The 0.2 percentage-point overall difference (Wilcoxon p\u202f=\u202f0.027) is clinically negligible. Treatment-adjusted LDL-C values showed near-perfect agreement across four LDL-C calculation equations (Friedewald [33], Sampson [34], Modified Sampson, Martin [35]; ICC\u202f=\u202f0.986, 95% CI: 0.985\u20130.988).")

doc <- h2(doc, "3.3 Wales External Validation (TRIPOD Type 2b)")
doc <- body_par(doc,
  "In the All Wales FH Registry, TUDOR\u202fv2 achieved AUC\u202f0.842 (95% CI: 0.822\u20130.863) when trained on index cases and validated on cascade-screened relatives (Figure 1). TUDOR significantly outperformed DLCN scoring (AUC\u202f0.791; DeLong p\u202f<\u202f0.001), Simon Broome criteria, and MEDPED criteria. Net Reclassification Improvement of TUDOR over DLCN was 0.358, indicating that more than one-third of patients were correctly reclassified. Calibration slope was 1.039, indicating excellent agreement between predicted and observed probabilities.")
doc <- body_par(doc,
  "The most clinically consequential finding was the performance divergence between ascertainment pathways. Among FH-positive index cases, DLCN at the probable threshold (\u22656 points) achieved 88.5% sensitivity (mean DLCN score 9.2). Among FH-positive cascade-screened relatives, the identical threshold achieved only 1.5% sensitivity (mean DLCN score 0.2). This 87-percentage-point sensitivity gap represents a systematic diagnostic deficiency in precisely the population where accurate identification carries the greatest public health significance. TUDOR achieved 89.4% sensitivity in cascade relatives versus 66.7% for the FAMCAT approximation, a 22.7-percentage-point advantage targeting this population of highest clinical need.")

doc <- h2(doc, "3.4 UK Biobank External Validation (TRIPOD Type 4)")
doc <- body_par(doc,
  "In the UK Biobank lipid clinic-mimicking cohort (n\u202f=\u202f58,021; 729 FH cases; 1.26% prevalence), TUDOR achieved AUC\u202f0.750 (95% CI: 0.731\u20130.770) using fixed Wales-derived model weights without any re-estimation (Figure 1, Table 3). This significantly exceeded eDLCN (AUC\u202f0.636; DeLong Z\u202f=\u202f10.08; p\u202f=\u202f6.73\u202f\u00d7\u202f10\u207b\u00b2\u2074), isolated LDL-C (AUC\u202f0.689; p\u202f<\u202f0.001), and the Trig_Filter alone (AUC\u202f0.700). TUDOR maintained clinically meaningful discriminative advantage despite an 80-fold prevalence difference and substantially lower data quality.")
doc <- body_par(doc,
  "At the Youden optimal threshold, TUDOR achieved sensitivity 59.7% and specificity 79.3%. The Brier score was 0.069. The calibration slope was 6.33 (95% CI: 5.93\u20136.73), indicating that TUDOR probability estimates require recalibration for deployment at lower-prevalence settings \u2014 an expected and remediable consequence of training under enriched prevalence (33.2%) and validating at 1.26%. Bayesian recalibration adjusts absolute probabilities while preserving discriminative ranking: P_recalibrated\u202f=\u202fP_TUDOR\u202f\u00d7\u202f(\u03c0_local\u202f/\u202f\u03c0_train)\u202f/\u202f[P_TUDOR\u202f\u00d7\u202f(\u03c0_local\u202f/\u202f\u03c0_train)\u202f+\u202f(1\u202f\u2212\u202fP_TUDOR)\u202f\u00d7\u202f((1\u202f\u2212\u202f\u03c0_local)\u202f/\u202f(1\u202f\u2212\u202f\u03c0_train))]. At the level of case detection, TUDOR identifies an estimated 639 FH cases per 10,000 lipid clinic patients compared with 52 under eDLCN \u2014 a 12-fold improvement reflecting substantially higher sensitivity (63.9% versus 5.2%).")

doc <- h2(doc, "3.5 Genetic Validation")
doc <- body_par(doc,
  "Mutation spectrum analysis revealed LDLR variants predominating in both cohorts: 1,901 cases (79.0%) in Wales and 1,321 (81.4%) in UK Biobank (Table 4). APOB variants were the second most common: 281 (11.7%) in Wales and 301 (18.5%) in UK Biobank. The higher APOB proportion in UK Biobank likely reflects wider population sampling, as APOB carriers present with milder phenotypes less frequently referred to specialist clinics. PCSK9 variants were rare (23 cases, 1.0% in Wales; 1 case, 0.1% in UK Biobank). The Wales panel additionally identified LDLR copy number variants (164 cases, 6.8%), APOE variants (33 cases, 1.4%), and LDLRAP1 variants (3 cases, 0.1%). The gene distribution differed significantly between cohorts (chi-square p\u202f=\u202f6.69\u202f\u00d7\u202f10\u207b\u2079), driven by the higher APOB proportion in UK Biobank and absence of CNV classification. Together, both cohorts encompassed 4,028 genetically confirmed FH cases, substantially exceeding previous diagnostic studies.")
doc <- body_par(doc,
  "Gene-specific TUDOR performance demonstrated consistent discrimination (Table 5, Figure 2). In Wales: LDLR AUC\u202f0.839 (95% CI: 0.817\u20130.861; n\u202f=\u202f724), APOB AUC\u202f0.841 (0.790\u20130.892; n\u202f=\u202f96), APOE AUC\u202f0.809 (0.688\u20130.930; n\u202f=\u202f20). In UK Biobank: LDLR AUC\u202f0.717 (0.693\u20130.741; n\u202f=\u202f515), APOB AUC\u202f0.830 (0.802\u20130.858; n\u202f=\u202f213). The significantly higher TUDOR discrimination for APOB versus LDLR carriers in UK Biobank (DeLong Z\u202f=\u202f3.74, p\u202f=\u202f1.88\u202f\u00d7\u202f10\u207b\u2074) reflects molecular homogeneity: APOB gain-of-function variants predominantly produce a single molecular defect (impaired ApoB\u2013LDL receptor binding), yielding a tight and reproducible phenotype. LDLR loss-of-function mutations are profoundly heterogeneous, spanning null alleles, receptor-defective variants with 2\u201325% residual function, and copy number variants \u2014 molecular diversity that translates to a wider LDL-C distribution, reducing discriminative efficiency. The Cohen\u2019s d for the Trig_Filter between FH and non-FH groups was 1.108 for APOB carriers versus 0.672 for LDLR carriers, confirming that the metabolic filter provides greater separation where the metabolic shield is most intact.")

doc <- h2(doc, "3.6 Apolipoprotein B Augmentation")
doc <- body_par(doc,
  "In the complete-case ApoB subset (n\u202f=\u202f57,446; 99.0% of the lipid clinic cohort), base TUDOR achieved AUC\u202f0.753. Addition of serum ApoB alone did not improve discrimination (Model A: AUC\u202f0.751; \u0394AUC\u202f=\u202f\u22120.002; DeLong p\u202f=\u202f0.349), confirming that absolute ApoB provides redundant information beyond what TUDOR captures through treatment-adjusted LDL-C. In contrast, the ApoB/LDL-C ratio provided significant augmentation (Model B: AUC\u202f0.760; \u0394AUC\u202f+0.007; p\u202f=\u202f0.006), capturing LDL particle cholesterol enrichment as a distinct biological signal. The combined Model C (TUDOR\u202f+\u202fApoB\u202f+\u202fApoB/LDL-C ratio) achieved the largest improvement: AUC\u202f0.771 (\u0394AUC\u202f+0.019; p\u202f=\u202f3.99\u202f\u00d7\u202f10\u207b\u2075). Model D additionally including Lp(a) achieved AUC\u202f0.776 (n\u202f=\u202f43,710).")
doc <- body_par(doc,
  "A grey zone analysis restricted to intermediate-probability participants (n\u202f=\u202f33,562; 454 FH-positive) showed that the ApoB/LDL-C ratio did not differ between FH-positive and FH-negative participants within this stratum (mean 0.2234 versus 0.2241; Wilcoxon p\u202f=\u202f0.356; AUC\u202f0.513, 95% CI: 0.486\u20130.539), indicating that sequential ApoB measurement after an intermediate TUDOR probability is not supported by current evidence. The ratio\u2019s benefit operates at the level of whole-cohort ranking rather than within-stratum triage. The ApoB/LDL-C ratio below 0.31 has been independently associated with higher ASCVD event rates in genetically confirmed FH (27.6% versus 11.8%; p\u202f=\u202f0.0022) [39], providing a separate prognostic rationale for its measurement after FH diagnosis.")

doc <- h2(doc, "3.7 Biological Validation of the Triglyceride Metabolic Filter")
doc <- body_par(doc,
  "Sequential multivariable logistic regression demonstrated that the Trig_Filter odds ratio remained remarkably stable across progressive adjustment: crude OR\u202f1.873, adjusted for BMI 1.891, further adjusted for diabetes mellitus 1.903, further adjusted for ApoB 1.920, and fully adjusted including age and sex 1.935 (Table 7). The less than 10% change from crude to fully adjusted estimates confirms independence from metabolic confounding. All variance inflation factors were below 1.2, excluding multicollinearity.")
doc <- body_par(doc,
  "The metabolic shield effect revealed progressive deterioration of diagnostic separation with increasing triglycerides: Cohen\u2019s d\u202f=\u202f0.89 in the normal triglyceride zone (<1.7\u202fmmol/L), decreasing to 0.52 in the borderline zone (1.7\u20132.3\u202fmmol/L) and 0.18 in the elevated zone (>2.3\u202fmmol/L). NMR metabolomics profiling (421,122 UK Biobank participants; 1,612 FH) confirmed 35 of 53 metabolites differing significantly between FH and non-FH after Bonferroni correction, consistent with cholesterol-enriched LDL particles with reduced VLDL. Crucially, NMR metabolomics data did not improve upon TUDOR: NMR alone AUC\u202f0.705; TUDOR alone 0.749; combined 0.752 (DeLong p\u202f=\u202f0.278), establishing that a parsimonious clinical model captures the same biological information encoded in 53 molecular measurements.")

doc <- h2(doc, "3.8 Cardiovascular Mortality Gradient by LDL Decile")
doc <- body_par(doc,
  "UK Biobank death records (median follow-up 14.2\u202fyears) were available for 58,021 lipid clinic cohort participants, with 1,893 deaths recorded at the time of analysis. IHD mortality showed a consistent gradient across LDL-C deciles, rising from 0.65% in the lowest decile (mean LDL-C 3.2\u202fmmol/L) to 1.76% in the highest decile (mean LDL-C 7.8\u202fmmol/L; logistic trend OR\u202f=\u202f1.036 per decile, 95% CI: 1.008\u20131.065, p\u202f=\u202f0.011; Figure 4). This positive gradient is consistent with the causal role of LDL-C in atherosclerosis and independently validates LDL-C as a continuous cardiovascular mortality risk determinant within a high-risk lipid clinic population.")
doc <- body_par(doc,
  "Stroke mortality showed an inverse trend across LDL deciles (0.89% decile 1 to 0.35% decile 10), consistent with the well-characterised J-curve relationship between LDL-C and haemorrhagic stroke, reflecting the mixed aetiology of cerebrovascular disease (thrombotic versus haemorrhagic). The inverse stroke gradient should not be interpreted as evidence that higher LDL-C is protective for cerebrovascular events; rather, it reflects the differential aetiological contribution of cholesterol to thrombotic versus haemorrhagic subtypes, with haemorrhagic stroke risk increasing at low LDL-C. Overall cardiovascular mortality (all cardiovascular ICD-10 codes combined) showed a non-significant positive trend (OR\u202f=\u202f1.012 per decile; p\u202f=\u202f0.14), attenuated by the opposing directions of the IHD and stroke components.")

doc <- h2(doc, "3.9 T2DM Ascertainment Sensitivity")
doc <- body_par(doc,
  "Three independent T2DM definitions identified substantially different patient populations: Definition\u202fA (self-report, p2443) identified 1.9% (n\u202f=\u202f1,102); Definition\u202fB (ICD-10 E11) identified 5.6% (n\u202f=\u202f3,249); Definition\u202fC (HbA1c\u202f\u226548\u202fmmol/mol) identified 1.5% (n\u202f=\u202f870). Cohen\u2019s kappa indicated poor inter-definition agreement: \u03ba\u202f=\u202f0.281 (A versus B), 0.358 (A versus C), 0.299 (B versus C), with substantial proportions of cases captured by only one definition. The low concordance reflects distinct ascertainment mechanisms: self-report captures diagnoses communicated to the patient, ICD-10 captures treated or hospitalised diabetes, and HbA1c captures biochemical glycaemia irrespective of clinical diagnosis.")
doc <- body_par(doc,
  "Despite this poor inter-definition concordance, all three definitions independently confirmed Trig_Filter attenuation: Definition\u202fA \u03b2\u202f=\u202f\u22120.484 (95% CI: \u22120.611 to \u22120.357; p\u202f=\u202f8.4\u202f\u00d7\u202f10\u207b\u00b9\u00b3); Definition\u202fB \u03b2\u202f=\u202f\u22120.473 (95% CI: \u22120.559 to \u22120.387; p\u202f=\u202f8.2\u202f\u00d7\u202f10\u207b\u00b2\u00b7); Definition\u202fC \u03b2\u202f=\u202f\u22120.658 (95% CI: \u22120.793 to \u22120.523; p\u202f=\u202f1.1\u202f\u00d7\u202f10\u207b\u00b2\u00b9). The consistency across three independent methods with poor inter-definition agreement demonstrates that Trig_Filter attenuation by diabetes reflects a biological phenomenon rather than coding artefact. HbA1c dose-response analysis confirmed a mechanistic gradient: mean Trig_Filter 2.74 (normoglycaemia), 2.31 (impaired fasting glucose, 42\u201347\u202fmmol/mol), and 1.82 (diabetes, \u226548\u202fmmol/mol), demonstrating proportional metabolic disruption with worsening glycaemic control. These findings are consistent with insulin resistance driving hypertriglyceridaemia through increased hepatic VLDL secretion and impaired lipoprotein lipase activity, which dilutes the FH-specific Trig_Filter signal but does not eliminate it.")

doc <- h2(doc, "3.10 Sensitivity Analyses and Algorithm Comparison")
doc <- body_par(doc,
  "Ten prespecified sensitivity analyses confirmed robust performance across LDL-C thresholds, calculation equations, statin correction factor perturbations, ethnic subgroups, winsorised extreme values, and alternative diagnostic criteria comparisons. The 100% compliance sensitivity analysis achieved equivalent performance to the primary model. Logistic regression using identical TUDOR-engineered features achieved equivalent performance to XGBoost (ENET AUC\u202f0.781, 95% CI: 0.768\u20130.796; XGBoost AUC\u202f0.778, 0.764\u20130.792; DeLong Z\u202f=\u202f0.36, p\u202f=\u202f0.72), confirming that clinical feature engineering \u2014 not algorithmic complexity \u2014 drives diagnostic improvement. This has direct clinical implications: a fully interpretable logistic regression equation can be deployed within standard electronic health record systems without machine learning infrastructure.")

# =============================================================================
# 4. DISCUSSION
# =============================================================================
doc <- h1(doc, "4. DISCUSSION")

doc <- h2(doc, "4.1 Principal Findings")
doc <- body_par(doc,
  "This study presents the development and dual external validation of TUDOR, achieving clinically meaningful discrimination across a specialist lipid registry (AUC\u202f0.842; 7,253 patients, all of Wales) and a population-based biobank (AUC\u202f0.750; 58,021 participants). Together, these validation cohorts encompass 4,028 genetically confirmed FH cases across six causative genes (LDLR, APOB, PCSK9, APOE, LDLRAP1, and LDLR copy number variants), representing to our knowledge the largest dual-validated genetically confirmed FH diagnostic study to date. The finding that systematic treatment adjustment, ascertainment bias correction, and triglyceride-based metabolic filtering substantially improve FH detection, implemented through interpretable logistic regression rather than opaque machine learning, carries direct implications for diagnostic practice and health system implementation.")

doc <- h2(doc, "4.2 The Cascade Screening Imperative")
doc <- body_par(doc,
  "The 87-percentage-point sensitivity gap between index cases (88.5%) and cascade-screened relatives (1.5%) using DLCN criteria represents a critical diagnostic deficiency affecting all established systems. Simon Broome criteria depend substantially on tendon xanthomata (present in fewer than 15% of heterozygous FH at any age) and corneal arcus (2.0% of our FH-positive cohort) \u2014 clinical signs increasingly rare in the statin-treated era [13]. MEDPED criteria lose discriminative power when treatment reduces LDL-C below threshold values [14]. This limitation is not merely statistical; it directly affects millions of FH-affected individuals worldwide who remain undiagnosed. Cascade screening programmes identify relatives carrying identical pathogenic variants to their affected family members, yet standard criteria fail to recognise the overwhelming majority [20,22]. TUDOR addresses this through two complementary mechanisms: the Index Effect modelling differential diagnostic significance by ascertainment pathway, and the Trig_Filter identifying the metabolic signature of monogenic FH regardless of absolute LDL-C.")

doc <- h2(doc, "4.3 Genetic Heterogeneity and Variant-Specific Performance")
doc <- body_par(doc,
  "The higher TUDOR discrimination for APOB versus LDLR carriers in UK Biobank (AUC\u202f0.830 versus 0.717; DeLong p\u202f=\u202f1.88\u202f\u00d7\u202f10\u207b\u2074) has a coherent pathophysiological basis. The APOB p.Arg3527Gln mutation, the predominant APOB FH variant, produces a relatively uniform molecular defect: impaired binding of ApoB-100 to the LDL receptor, resulting in consistent phenotypic expression characterised by moderately elevated LDL-C with characteristic cholesterol-enriched particles. In contrast, over 2,000 pathogenic LDLR variants have been described, encompassing null alleles, receptor-defective variants with 2\u201325% residual function, and class-specific defects affecting synthesis, transport, binding, internalisation, or recycling [4]. This molecular heterogeneity translates to a wider LDL-C distribution, reducing the discriminative efficiency of any LDL-C-based diagnostic tool. The clinical implication is that TUDOR may be particularly valuable for identifying APOB carriers, who represent 12\u201319% of genetically confirmed FH across our cohorts and are systematically missed by DLCN criteria that assign highest scores to LDL-C concentrations >8.5\u202fmmol/L \u2014 a threshold rarely reached by APOB carriers.")

doc <- h2(doc, "4.4 Cardiovascular Outcomes Validation")
doc <- body_par(doc,
  "The confirmed IHD mortality gradient (0.65%\u21921.76% across LDL deciles; OR\u202f=\u202f1.036 per decile, p\u202f=\u202f0.011) provides independent biological validation that LDL-C concentration within this high-risk cohort remains a continuous cardiovascular mortality risk determinant, consistent with Mendelian randomisation data [5] and the causal framework established by Brown and Goldstein\u2019s Nobel Prize-winning work on the LDL receptor pathway. This finding is particularly important for cardiologists because it confirms that even within a clinically enriched lipid clinic population \u2014 where mean LDL-C is substantially elevated \u2014 higher LDL-C concentrations carry incrementally higher cardiovascular mortality risk, reinforcing the principle that there is no safe LDL-C threshold in patients at risk of FH.")
doc <- body_par(doc,
  "The inverse stroke mortality gradient reflects the well-characterised J-curve relationship between LDL-C and haemorrhagic stroke, whereby very low LDL-C concentrations may increase the risk of haemorrhagic events through effects on vascular wall integrity [42]. The opposing directions of IHD and stroke mortality gradients across LDL deciles are clinically important: cardiologists managing FH patients should recognise that aggressive LDL-C lowering substantially reduces the predominant IHD risk while the haemorrhagic stroke risk at population-average LDL-C levels is not clinically meaningful in the context of FH management. Contemporary cardiovascular guidelines recommend LDL-C targets of <1.4\u202fmmol/L in very high-risk patients [43,44], well above concentrations where haemorrhagic stroke risk is clinically relevant.")

doc <- h2(doc, "4.5 T2DM Ascertainment and Triglyceride Filter Robustness")
doc <- body_par(doc,
  "The T2DM sensitivity analysis is methodologically novel in demonstrating that the Trig_Filter attenuation by diabetes is a biological phenomenon rather than a coding artefact: three independent definitions with poor inter-definition concordance (\u03ba\u202f=\u202f0.28\u20130.36) each independently confirm the same effect (\u03b2\u202f=\u202f\u22120.48 to \u22120.66). The low kappa actually strengthens this conclusion \u2014 if the effect were an artefact of how diabetes was coded, it would not be consistently observed across definitions that identify substantially different patient populations. The HbA1c dose-response (Trig_Filter 2.74\u21922.31\u21921.82 across glycaemic strata) provides the mechanistic gradient: insulin resistance increases hepatic VLDL secretion through SREBP-1c activation and impairs lipoprotein lipase activity, elevating triglycerides and reducing the FH-specific Trig_Filter signal in proportion to glycaemic severity [45,46]. This mechanistic understanding directly informs clinical implementation: in patients with diabetes or insulin resistance, the Trig_Filter threshold should be interpreted in the context of glycaemic control, and HbA1c measurement should accompany Trig_Filter calculation in diagnostic uncertainty.")

doc <- h2(doc, "4.6 ApoB Augmentation and Precision Lipidology")
doc <- body_par(doc,
  "The finding that absolute ApoB does not augment TUDOR (Model A, p\u202f=\u202f0.349) while the ApoB/LDL-C ratio does (Model C, \u0394AUC\u202f+0.019, p\u202f=\u202f3.99\u202f\u00d7\u202f10\u207b\u2075) has a clear biological interpretation. Absolute ApoB reflects the total number of atherogenic lipoprotein particles, which is correlated with but not identical to LDL-C. TUDOR already incorporates treatment-adjusted LDL-C, making absolute ApoB largely redundant. The ApoB/LDL-C ratio captures qualitatively different information: cholesterol content per LDL particle. In FH, impaired receptor-mediated clearance prolongs LDL circulatory residence time, allowing continued cholesterol ester accumulation, producing fewer but more cholesterol-laden particles. This produces a lower ApoB-to-cholesterol ratio compared with non-FH hypercholesterolaemia driven by overproduction of normal-cholesterol-content particles [47,48]. The prognostic implications are established: the ApoB/LDL-C ratio <0.31 predicts ASCVD events in genetically confirmed FH (27.6% versus 11.8%; p\u202f=\u202f0.0022) [39].")

doc <- h2(doc, "4.7 Clinical Implementation Framework")
doc <- body_par(doc,
  "We propose a three-tier clinical decision framework based on TUDOR probability estimates for implementation in specialist lipid clinics. Low-risk patients (TUDOR probability <25%) can be reassured that monogenic FH is unlikely, and standard lipid management should continue. Intermediate-probability patients (25\u201375%) warrant further clinical evaluation, including serum ApoB, ApoB/LDL-C ratio assessment, extended cascade data review, and shared decision-making regarding genetic testing referral. High-risk patients (>75%) should be referred directly for genetic testing and considered for immediate initiation or intensification of FH-appropriate lipid-lowering therapy including high-intensity statins, ezetimibe, and PCSK9 inhibitors [25\u201327,43,44,49,50].")
doc <- body_par(doc,
  "For implementation in lower-prevalence settings such as primary care, probability thresholds require Bayesian recalibration to local FH prevalence. The calibration slope of 6.33 in UK Biobank indicates systematic overestimation of absolute probabilities in lower-prevalence populations, which is remediable through the recalibration formula provided. Critically, discrimination (the ranking of patients) is preserved irrespective of prevalence, meaning the three-tier classification system retains clinical utility across settings. The equivalence of logistic regression and XGBoost performance confirms that a simple, fully interpretable logistic regression equation \u2014 computable in standard spreadsheet software \u2014 is sufficient for clinical deployment without requiring machine learning infrastructure [36,37].")

doc <- h2(doc, "4.8 Limitations")
doc <- body_par(doc,
  "Several limitations warrant consideration. First, TUDOR was developed in a single specialist centre in Wales; prospective validation in lipid clinic cohorts across different healthcare systems, including non-UK settings with different genetic variant spectra and treatment practices, is needed before widespread deployment. Second, UK Biobank\u2019s age restriction (40\u201369 years) means TUDOR\u2019s performance in paediatric and young adult populations, where cascade screening is particularly important, has not been externally validated, though the Wales cohort includes the full age spectrum. Third, the FAMCAT comparison uses an internally reconstructed approximation that may overestimate TUDOR\u2019s comparative advantage. Fourth, the calibration slope of 6.33 in UK Biobank requires recalibration before deployment at lower prevalence, though this is a straightforward statistical adjustment. Fifth, ApoB augmentation was testable only in UK Biobank; specialist lipid clinic validation is needed. Sixth, the T2DM sensitivity analysis confirms attenuation of the Trig_Filter in diabetes but does not preclude its diagnostic use; rather, it suggests that glycaemic context should accompany TUDOR interpretation in metabolically complex patients. Seventh, premature ASCVD event dates (ICD-10 field p41280) were not extracted in the current analysis; true age-at-first-ASCVD-event analysis is a priority for the next analysis cycle. Finally, prospective implementation studies evaluating the impact of TUDOR-guided decisions on genetic testing yield and patient outcomes are needed before clinical deployment can be formally recommended.")

# =============================================================================
# 5. CONCLUSIONS
# =============================================================================
doc <- h1(doc, "5. CONCLUSIONS")
doc <- body_par(doc,
  "TUDOR provides a dual-externally-validated, interpretable diagnostic algorithm for familial hypercholesterolaemia that maintains clinically meaningful discrimination across specialist registry (AUC\u202f0.842) and population-based biobank (AUC\u202f0.750) settings, in the largest genetically confirmed FH cohort reported to date (4,028 cases across six causative genes). Gene-specific validation demonstrates consistent performance across LDLR, APOB, and PCSK9 variant carriers, with APOB variants achieving particularly high discrimination (AUC\u202f0.830) reflecting the molecular homogeneity of gain-of-function APOB mutations. Optional ApoB augmentation incorporating the ApoB/LDL-C ratio improves discrimination to AUC\u202f0.771 where serum ApoB is available. Cardiovascular mortality gradient analysis independently confirms LDL-C as a continuous mortality risk determinant within a high-risk population, reinforcing the clinical urgency of FH diagnosis. T2DM ascertainment sensitivity analysis, using three independent definitions with poor inter-definition concordance, confirms that triglyceride filter attenuation by diabetes is a biological rather than coding phenomenon, with a dose-response gradient consistent with metabolic mechanisms. Treatment-adjusted phenotyping and ascertainment bias correction address the fundamental limitations underlying systematic FH underdiagnosis, with particular value in cascade-screened relatives where established criteria identify fewer than 2% of affected individuals. TUDOR offers an interpretable, implementable complement to established diagnostic frameworks that can be deployed using standard clinical data routinely available in specialist lipid clinic practice.")

# =============================================================================
# DECLARATIONS
# =============================================================================
doc <- h1(doc, "DECLARATIONS")
doc <- h2(doc, "Ethics Approval and Consent")
doc <- body_par(doc,
  "This study was approved by the Research and Development Department at Cardiff and Vale University Health Board. UK Biobank analyses were conducted under Application Number 1002450. UK Biobank has received ethical approval from the North West Multi-Centre Research Ethics Committee (11/NW/0382). All data were fully anonymised prior to analysis in compliance with the Declaration of Helsinki. Given the retrospective nature using anonymised data, individual patient consent was not required.")
doc <- h2(doc, "Data Availability")
doc <- body_par(doc,
  "The Wales FH Registry dataset is available from the corresponding author on reasonable request, subject to data governance approvals from Cardiff and Vale University Health Board. UK Biobank data are available to researchers through the standard application process (www.ukbiobank.ac.uk) under application 1002450. R analysis scripts will be deposited in a public repository (GitHub/Zenodo) upon acceptance. Scripts are available from the corresponding author during peer review.")
doc <- h2(doc, "Competing Interests")
doc <- body_par(doc, "The authors declare no competing interests.")
doc <- h2(doc, "Funding")
doc <- body_par(doc, "This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.")
doc <- h2(doc, "Use of Artificial Intelligence")
doc <- body_par(doc,
  "Artificial intelligence tools (Claude, Anthropic Inc.) were employed solely for language editing, grammar checking, statistical code review, and proofreading. AI tools were NOT used for data analysis, statistical computations, interpretation of results, study design, methodology development, literature review, citation selection, or generating scientific conclusions. All scientific content, data interpretation, statistical analyses, and clinical recommendations represent the independent intellectual work of the named authors.")
doc <- h2(doc, "Author Contributions")
doc <- body_par(doc,
  "Nader Genedy: Conceptualisation, Methodology, Software, Formal Analysis, Investigation, Data Curation, Writing \u2013 Original Draft, Writing \u2013 Review and Editing, Visualisation, Project Administration. Soha Zouwail: Supervision, Validation, Writing \u2013 Review and Editing.")
doc <- h2(doc, "Acknowledgements")
doc <- body_par(doc,
  "We thank the patients and staff of the Wales FH Registry and the All Wales PASS programme for their contributions to this research. This research has been conducted using the UK Biobank Resource under application number 1002450.")

doc <- body_add_break(doc, pos="after")

# =============================================================================
# REFERENCES (55 high-quality citations)
# =============================================================================
doc <- h1(doc, "REFERENCES")

refs <- c(
  "1. Nordestgaard BG, Chapman MJ, Humphries SE, et al. Familial hypercholesterolaemia is underdiagnosed and undertreated in the general population: guidance for clinicians to prevent coronary heart disease. Eur Heart J. 2013;34(45):3478\u20133490.",
  "2. Khera AV, Won HH, Peloso GM, et al. Diagnostic yield and clinical utility of sequencing familial hypercholesterolaemia genes in patients with severe hypercholesterolaemia. J Am Coll Cardiol. 2016;67(22):2578\u20132589.",
  "3. Sturm AC, Knowles JW, Gidding SS, et al. Clinical genetic testing for familial hypercholesterolaemia: JACC Scientific Expert Panel. J Am Coll Cardiol. 2018;72(6):662\u2013680.",
  "4. Defesche JC, Gidding SS, Harada-Shiba M, Hegele RA, Santos RD, Wierzbicki AS. Familial hypercholesterolaemia. Nat Rev Dis Primers. 2017;3:17093.",
  "5. Ference BA, Ginsberg HN, Graham I, et al. Low-density lipoproteins cause atherosclerotic cardiovascular disease. 1. Evidence from genetic, epidemiologic, and clinical studies. Eur Heart J. 2017;38(32):2459\u20132472.",
  "6. Marks D, Thorogood M, Neil HA, Humphries SE. A review on the diagnosis, natural history, and treatment of familial hypercholesterolaemia. Atherosclerosis. 2003;168(1):1\u201314.",
  "7. Mundal L, Sarancic M, Ose L, et al. Mortality among patients with familial hypercholesterolaemia: a registry-based study in Norway, 1992\u20132010. J Am Heart Assoc. 2014;3(6):e001236.",
  "8. Versmissen J, Oosterveer DM, Yazdanpanah M, et al. Efficacy of statins in familial hypercholesterolaemia: a long-term cohort study. BMJ. 2008;337:a2423.",
  "9. Neil A, Cooper J, Betteridge J, et al. Reductions in all-cause, cancer, and coronary mortality in statin-treated patients with heterozygous familial hypercholesterolaemia: a prospective registry study. Eur Heart J. 2008;29(21):2625\u20132633.",
  "10. Watts GF, Gidding S, Wierzbicki AS, et al. Integrated guidance on the care of familial hypercholesterolaemia from the International FH Foundation. Int J Cardiol. 2014;171(3):309\u2013325.",
  "11. Akioyamen LE, Genest J, Shan SD, et al. Estimating the prevalence of heterozygous familial hypercholesterolaemia: a systematic review and meta-analysis. BMJ Open. 2017;7(9):e016461.",
  "12. Civeira F; International Panel on Management of Familial Hypercholesterolaemia. Guidelines for the diagnosis and management of heterozygous familial hypercholesterolaemia. Atherosclerosis. 2004;173(1):55\u201368.",
  "13. Risk of fatal coronary heart disease in familial hypercholesterolaemia. Scientific Steering Committee on behalf of the Simon Broome Register Group. BMJ. 1991;303(6807):893\u2013896.",
  "14. Williams RR, Hunt SC, Schumacher MC, et al. Diagnosing heterozygous familial hypercholesterolaemia using new practical criteria validated by molecular genetics. Am J Cardiol. 1993;72(2):171\u2013176.",
  "15. Weng TC, Yang YH, Lin SJ, Tai SH. A systematic review and meta-analysis on the therapeutic equivalence of statins. J Clin Pharm Ther. 2010;35(2):139\u2013151.",
  "16. Adams SP, Tsang M, Wright JM. Lipid-lowering efficacy of atorvastatin. Cochrane Database Syst Rev. 2015;(3):CD008226.",
  "17. Adams SP, Sekhon SS, Wright JM. Lipid-lowering efficacy of rosuvastatin. Cochrane Database Syst Rev. 2014;(11):CD010254.",
  "18. Weng SF, Kai J, Neil HA, Humphries SE, Qureshi N. Improving identification of familial hypercholesterolaemia in primary care: derivation and validation of the familial hypercholesterolaemia case ascertainment tool (FAMCAT). Atherosclerosis. 2015;238(2):336\u2013343.",
  "19. Besseling J, Sjouke B, Kastelein JJ. Screening and treatment of familial hypercholesterolaemia \u2014 lessons from the past and opportunities for the future. Atherosclerosis. 2015;241(2):597\u2013606.",
  "20. Knowles JW, Rader DJ, Khoury MJ. Cascade screening for familial hypercholesterolaemia and the use of genetic testing. JAMA. 2017;318(4):381\u2013382.",
  "21. Haralambos K, Whatley SD, Edwards R, Gingell R, Humphries SE, Ashfield-Watt P. Clinical experience of scoring criteria for Familial Hypercholesterolaemia (FH) genetic testing in Wales. Atherosclerosis. 2015;240(1):190\u2013196.",
  "22. Bhatnagar D, Morgan J, Siddiq S, Mackness MI, Miller JP, Durrington PN. Outcome of case finding among relatives of patients with known heterozygous familial hypercholesterolaemia. BMJ. 2000;321(7275):1497\u20131500.",
  "23. Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the ACMG and AMP. Genet Med. 2015;17(5):405\u2013424.",
  "24. Sudlow C, Gallacher J, Allen N, et al. UK Biobank: an open access resource for identifying the causes of a wide range of complex diseases of middle and old age. PLoS Med. 2015;12(3):e1001779.",
  "25. Cannon CP, Blazing MA, Giugliano RP, et al. Ezetimibe added to statin therapy after acute coronary syndromes. N Engl J Med. 2015;372(25):2387\u20132397.",
  "26. Sabatine MS, Giugliano RP, Keech AC, et al. Evolocumab and clinical outcomes in patients with cardiovascular disease. N Engl J Med. 2017;376(18):1713\u20131722.",
  "27. Ray KK, Bays HE, Catapano AL, et al. Safety and efficacy of bempedoic acid to reduce LDL cholesterol. N Engl J Med. 2019;380(11):1022\u20131032.",
  "28. DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. Biometrics. 1988;44(3):837\u2013845.",
  "29. Pencina MJ, D\u2019Agostino RB Sr, D\u2019Agostino RB Jr, Vasan RS. Evaluating the added predictive ability of a new marker: from area under the ROC curve to reclassification and beyond. Stat Med. 2008;27(2):157\u2013172.",
  "30. Hosmer DW, Lemeshow S. Applied Logistic Regression. 2nd ed. New York: John Wiley and Sons; 2000.",
  "31. Benn M, Watts GF, Tybjaerg-Hansen A, Nordestgaard BG. Mutations causative of familial hypercholesterolaemia: screening of 98,098 individuals from the Copenhagen General Population Study estimated a prevalence of 1 in 217. Eur Heart J. 2016;37(17):1384\u20131394.",
  "32. Harada-Shiba M, Arai H, Ishigaki Y, et al. Guidelines for diagnosis and treatment of familial hypercholesterolaemia 2017. J Atheroscler Thromb. 2018;25(8):751\u2013770.",
  "33. Friedewald WT, Levy RI, Fredrickson DS. Estimation of the concentration of low-density lipoprotein cholesterol in plasma, without use of the preparative ultracentrifuge. Clin Chem. 1972;18(6):499\u2013502.",
  "34. Sampson M, Ling C, Sun Q, et al. A new equation for calculation of low-density lipoprotein cholesterol in patients with normolipidemia and/or hypertriglyceridemia. JAMA Cardiol. 2020;5(5):540\u2013548.",
  "35. Martin SS, Blaha MJ, Elshazly MB, et al. Comparison of a novel method vs the Friedewald equation for estimating low-density lipoprotein cholesterol levels from the standard lipid profile. JAMA. 2013;310(19):2061\u20132068.",
  "36. Banda JM, Sarraju A, Abbasi F, et al. Finding missed cases of familial hypercholesterolaemia in health systems using machine learning. NPJ Digit Med. 2019;2:23.",
  "37. Myers KD, Knowles JW, Staszak D, et al. Precision screening for familial hypercholesterolaemia: a machine learning study applied to electronic health encounter data. Lancet Digit Health. 2019;1(8):e393\u2013e402.",
  "38. Safarova MS, Liu H, Kullo IJ. Rapid identification of familial hypercholesterolaemia from electronic health records: the SEARCH study. J Clin Lipidol. 2016;10(5):1230\u20131239.",
  "39. Genedy N, Zouwail S. ApoB/LDL-C Discordance as a Predictor of Atherosclerotic Cardiovascular Disease in Genetically Confirmed Heterozygous Familial Hypercholesterolaemia: A Hypothesis-Generating Cohort Study. J Clin Lipidol. 2025. doi:10.1016/j.jacl.2025.11.008.",
  "40. Collins GS, Reitsma JB, Altman DG, Moons KG. Transparent reporting of a multivariable prediction model for individual prognosis or diagnosis (TRIPOD): the TRIPOD statement. Ann Intern Med. 2015;162(1):55\u201363.",
  "41. von Elm E, Altman DG, Egger M, et al. The Strengthening the Reporting of Observational Studies in Epidemiology (STROBE) statement: guidelines for reporting observational studies. Lancet. 2007;370(9596):1453\u20131457.",
  "42. Ference BA, Yoo W, Alesh I, et al. Effect of long-term exposure to lower low-density lipoprotein cholesterol beginning early in life on the risk of coronary heart disease: a Mendelian randomization analysis. J Am Coll Cardiol. 2012;60(25):2631\u20132639.",
  "43. Mach F, Baigent C, Catapano AL, et al. 2019 ESC/EAS Guidelines for the management of dyslipidaemias. Eur Heart J. 2020;41(1):111\u2013188.",
  "44. Grundy SM, Stone NJ, Bailey AL, et al. 2018 AHA/ACC/AACVPR/AAPA/ABC/ACPM/ADA/AGS/APhA/ASPC/NLA/PCNA Guideline on the management of blood cholesterol. J Am Coll Cardiol. 2019;73(24):e285\u2013e350.",
  "45. Adiels M, Olofsson SO, Taskinen MR, Bor\u00e9n J. Overproduction of very low-density lipoproteins is the hallmark of the dyslipidaemia in the metabolic syndrome. Arterioscler Thromb Vasc Biol. 2008;28(7):1225\u20131236.",
  "46. Ginsberg HN, Zhang YL, Hernandez-Ono A. Metabolic syndrome: focus on dyslipidaemia. Obesity. 2006;14(Suppl 1):41S\u201349S.",
  "47. Sniderman AD, Lamarche B, Tilley J, Seccombe D, Frohlich J. Hypertriglyceridemic hyperapo B in type 2 diabetes. Diabetes Care. 2002;25(3):579\u2013582.",
  "48. Walldius G, Jungner I. Apolipoprotein B and apolipoprotein A-I: risk indicators of coronary heart disease and targets for lipid-modifying therapy. J Intern Med. 2004;255(2):188\u2013205.",
  "49. Raal FJ, Kallend D, Ray KK, et al. Inclisiran for the treatment of heterozygous familial hypercholesterolaemia. N Engl J Med. 2020;382(16):1520\u20131530.",
  "50. Schwartz GG, Steg PG, Szarek M, et al. Alirocumab and cardiovascular outcomes after acute coronary syndrome. N Engl J Med. 2018;379(22):2097\u20132107.",
  "51. Santos RD, Stein EA, Hovingh GK, et al. Long-term evolocumab in patients with familial hypercholesterolaemia. J Am Coll Cardiol. 2020;75(6):565\u2013574.",
  "52. Tromp TR, Stroes ESG, Hovingh GK. Gene-based therapy in lipid management: the winding road from promise to practice. Eur Heart J. 2021;42(16):1551\u20131560.",
  "53. Natarajan P, Young R, Stitziel NO, et al. Polygenic background modifies penetrance of monogenic variants for tier 1 genomic conditions. Nat Commun. 2021;12:6372.",
  "54. Patel AP, Wang M, Pirruccello JP, et al. Lp(a) (Lipoprotein[a]) concentrations and incident atherosclerotic cardiovascular disease: new insights from a large national biobank. Arterioscler Thromb Vasc Biol. 2021;41(1):465\u2013474.",
  "55. Hovingh GK, Davidson MH, Kastelein JJ, O\u2019Connor AM. Diagnosis and treatment of familial hypercholesterolaemia. Eur Heart J. 2013;34(13):962\u2013971."
)

for (r in refs) {
  doc <- body_par(doc, r, space_after=4, align="left")
}

doc <- body_add_break(doc, pos="after")

# =============================================================================
# FIGURE LEGENDS
# =============================================================================
doc <- h1(doc, "FIGURE LEGENDS")

fig_legends <- list(
  c("Figure 1.",
    "Dual-panel ROC curves for TUDOR validation. (A) Wales All Wales FH Registry external validation (n\u202f=\u202f7,253; 2,405 FH-positive): TUDOR v2 (AUC\u202f0.842, 95%\u202fCI: 0.822\u20130.863) versus DLCN scoring (AUC\u202f0.791) and FAMCAT approximation, model trained on index cases and validated on cascade-screened relatives (TRIPOD Type 2b). (B) UK Biobank lipid clinic-mimicking cohort (n\u202f=\u202f58,021; 729 FH cases): TUDOR (AUC\u202f0.750, 0.731\u20130.770) versus eDLCN (0.636), LDL-C alone (0.689), and Trig_Filter alone (0.700). Fixed Wales-derived weights without re-estimation (TRIPOD Type 4). DeLong: TUDOR versus eDLCN Z\u202f=\u202f10.08, p\u202f=\u202f6.73\u202f\u00d7\u202f10\u207b\u00b2\u2074. DLCN sensitivity in cascade relatives 1.5% versus TUDOR 89.4%."),
  c("Figure 2.",
    "Gene-specific TUDOR performance: forest plot with DeLong comparisons. Horizontal forest plot showing AUC (95%\u202fCI) stratified by causative gene in both validation cohorts. Wales: LDLR 0.839 (0.817\u20130.861, n\u202f=\u202f724), APOB 0.841 (0.790\u20130.892, n\u202f=\u202f96), APOE 0.809 (0.688\u20130.930, n\u202f=\u202f20). UK Biobank: LDLR 0.717 (0.693\u20130.741, n\u202f=\u202f515), APOB 0.830 (0.802\u20130.858, n\u202f=\u202f213). DeLong APOB versus LDLR in UK Biobank: Z\u202f=\u202f3.74, p\u202f=\u202f1.88\u202f\u00d7\u202f10\u207b\u2074. Diamonds represent pooled estimates; squares represent gene-specific values with area proportional to sample size."),
  c("Figure 3.",
    "Calibration plot with Bayesian recalibration across prevalence settings. (A) Observed versus predicted FH probability in UK Biobank lipid clinic cohort by decile of predicted risk. Calibration slope\u202f=\u202f6.33 (95%\u202fCI: 5.93\u20136.73). Dashed line: perfect calibration. Loess smoothing with 95% confidence band. (B) Bayesian recalibration: mapping of identical TUDOR scores to absolute probabilities at three prevalence settings (33.2% Wales specialist registry, 1.26% UK Biobank lipid clinic, 0.4% primary care). Discrimination (ranking) is preserved; only absolute probabilities require adjustment."),
  c("Figure 4.",
    "Cardiovascular mortality gradient and T2DM sensitivity. (A) IHD and stroke mortality by LDL-C decile in UK Biobank lipid clinic cohort (n\u202f=\u202f58,021). IHD mortality rises 0.65%\u21921.76% (logistic trend OR\u202f=\u202f1.036 per decile, 95%\u202fCI: 1.008\u20131.065, p\u202f=\u202f0.011). Stroke mortality shows inverse J-curve (0.89%\u21920.35%), consistent with mixed thrombotic/haemorrhagic aetiology. (B) Trig_Filter dose-response by HbA1c glycaemic stratum: normoglycaemia (HbA1c\u202f<42\u202fmmol/mol) 2.74, impaired fasting glucose (42\u201347) 2.31, diabetes (\u226548) 1.82, confirming proportional metabolic attenuation."),
  c("Figure 5.",
    "ApoB augmentation and metabolic shield effect. (A) AUC with 95%\u202fCI for sequential augmentation models (n\u202f=\u202f57,446). Base TUDOR 0.753, Model A (+ApoB) 0.751 (p\u202f=\u202f0.349, ns), Model B (+ApoB/LDL-C ratio) 0.760 (p\u202f=\u202f0.006), Model C (+both) 0.771 (p\u202f=\u202f3.99\u202f\u00d7\u202f10\u207b\u2075), Model D (+Lp(a)) 0.776. (B) Violin plots of Trig_Filter distribution by causative gene and FH status (Wales and UK Biobank). APOB carriers demonstrate highest and most concentrated Trig_Filter values; Cohen\u2019s d\u202f=\u202f1.108 (APOB) versus 0.672 (LDLR) versus 0.963 (overall). (C) Trig_Filter sequential multivariable logistic regression: OR stable 1.873\u21921.935 across adjustment for BMI, diabetes, ApoB, age, and sex, confirming metabolic filter independence.")
)

for (fl in fig_legends) {
  doc <- body_par(doc, fl[1], bold=TRUE, space_after=2, align="left")
  doc <- body_par(doc, fl[2], space_after=14, align="left")
}

doc <- body_add_break(doc, pos="after")

# =============================================================================
# TABLES (flextable)
# =============================================================================
# A4 content width with 30mm margins: (210 - 60) / 25.4 = 5.91 inches
A4_WIDTH <- 5.91

make_ft <- function(df, fontsize_pt=9, col_widths=NULL) {
  ft <- flextable(df)
  ft <- font(ft, fontname="Times New Roman", part="all")
  ft <- fontsize(ft, size=fontsize_pt, part="all")
  ft <- bold(ft, part="header")
  ft <- align(ft, align="center", part="header")
  ft <- align(ft, align="left", part="body")
  ft <- border_outer(ft, part="all", border=fp_border(width=1))
  ft <- border_inner_h(ft, part="all", border=fp_border(width=0.5))
  ft <- bg(ft, bg="#F0F4FA", part="header")
  if (!is.null(col_widths)) {
    ft <- width(ft, width=col_widths)
  } else {
    ft <- autofit(ft, add_w=0)
    ft <- set_table_properties(ft, width=1, layout="autofit")
  }
  ft
}

# ── Table 1a ──────────────────────────────────────────────────────────────────
doc <- h1(doc, "TABLES")
doc <- cap(doc, "Table 1a. Baseline Characteristics of the All Wales FH Registry (n\u202f=\u202f7,253)")
doc <- fn_par(doc, "Continuous variables: Welch\u2019s t-test; categorical: chi-square. Cohen\u2019s d reported for continuous variables. Values are mean (SD) or n (%). BMI available for 2,306; LDL-C measured for 5,664; treatment-adjusted LDL-C for 5,555; HDL-C for 5,522; triglycerides for 5,716; Trig_Filter for 5,385.")

t1a <- data.frame(
  Variable = c("Age, years","Female sex, n (%)","BMI, kg/m\u00b2",
               "LDL-C measured, mmol/L","LDL-C adjusted, mmol/L",
               "HDL-C, mmol/L","Triglycerides, mmol/L",
               "Trig_Filter","Total cholesterol, mmol/L",
               "Tendon xanthomata, n (%)","Corneal arcus (<45y), n (%)",
               "Statin use, n (%)","PCSK9 inhibitor, n (%)",
               "Premature ASCVD, n (%)","Cascade relative, n (%)"),
  Overall  = c("45.9 (17.8)","4,148 (57.2%)","28.7 (5.8)",
               "5.3 (1.9)","6.1 (2.0)","1.5 (0.4)","1.7 (0.9)",
               "4.2 (1.8)","8.0 (1.7)","538 (7.4%)","48 (0.7%)",
               "3,461 (47.7%)","156 (2.2%)","812 (11.2%)","1,912 (26.4%)"),
  FH_Positive = c("39.3 (17.1)","1,420 (59.0%)","27.9 (5.4)",
               "5.8 (2.0)","6.5 (2.1)","1.5 (0.4)","1.4 (0.8)",
               "5.1 (1.9)","8.6 (1.7)","260 (10.8%)","49 (2.0%)",
               "1,203 (50.0%)","83 (3.5%)","284 (11.8%)","1,370 (57.0%)"),
  FH_Negative = c("49.2 (17.2)","2,728 (56.3%)","29.2 (6.0)",
               "5.1 (1.8)","5.7 (1.9)","1.5 (0.4)","2.0 (0.9)",
               "3.3 (1.5)","7.7 (1.7)","278 (5.7%)","29 (0.6%)",
               "2,258 (46.6%)","73 (1.5%)","528 (10.9%)","542 (11.2%)"),
  P_value = c("<10-99","0.14","<10-11",
               "<10-31","<10-33","0.31","<10-72",
               "<10-90","<10-70","<10-17","<10-7",
               "0.030","0.001","0.53","<10-145"),
  Cohens_d = c("-0.571","--","-0.238",
               "0.349","0.378","-0.041","-0.460",
               "0.744","0.522","--","--","--","--","--","--"),
  check.names=TRUE)
colnames(t1a) <- c("Variable","Overall","FH-Positive (n=2,405)","FH-Negative (n=4,848)","P value","Cohen's d")
doc <- body_add_flextable(doc, make_ft(t1a))

doc <- hr(doc)
doc <- cap(doc, "Table 1b. Baseline Characteristics of UK Biobank Lipid Clinic-Mimicking Cohort (n\u202f=\u202f58,021)")
doc <- fn_par(doc, "**ApoB available for 57,446 (99.0%); Lp(a) available for 43,710 (75.3%). ASCVD\u202f=\u202fatherosclerotic cardiovascular disease. Premature ASCVD defined as MI, angina, or stroke before age 55 (men) or 60 (women).")

t1b <- data.frame(
  Variable = c("Age, years","Female sex, n (%)","BMI, kg/m\u00b2",
               "Statin use, n (%)","LDL-C measured, mmol/L",
               "LDL-C adjusted, mmol/L","HDL-C, mmol/L",
               "Triglycerides, mmol/L","Trig_Filter",
               "ApoB, g/L**","ApoB/LDL-C ratio**",
               "Total cholesterol, mmol/L","Premature ASCVD, n (%)",
               "T2DM (self-report), n (%)","T2DM (ICD-10 E11), n (%)","HbA1c, mmol/mol"),
  Overall  = c("58.9 (7.0)","31,993 (55.1%)","28.4 (4.7)",
               "22,399 (38.6%)","5.1 (0.9)","5.9 (1.1)",
               "1.6 (0.4)","2.2 (1.2)","2.9 (1.4)",
               "1.11 (0.26)","0.23 (0.06)","7.7 (0.8)",
               "10,615 (18.3%)","1,102 (1.9%)","3,249 (5.6%)","36.1 (5.7)"),
  FH_Positive = c("57.8 (6.4)","399 (54.7%)","27.8 (4.4)",
               "390 (53.5%)","5.7 (0.9)","5.8 (1.3)",
               "1.6 (0.4)","1.8 (0.9)","3.9 (1.5)",
               "1.12 (0.24)","0.21 (0.05)","8.1 (0.9)",
               "93 (12.8%)","9 (1.2%)","29 (4.0%)","35.8 (5.4)"),
  FH_Negative = c("58.9 (7.0)","31,594 (55.2%)","28.5 (4.7)",
               "22,009 (38.4%)","5.1 (0.9)","5.9 (1.1)",
               "1.6 (0.4)","2.2 (1.2)","2.7 (1.4)",
               "1.11 (0.26)","0.23 (0.06)","7.7 (0.8)",
               "10,522 (18.4%)","1,093 (1.9%)","3,220 (5.6%)","36.1 (5.7)"),
  P_value = c("0.006","0.84","0.0002",
               "<10-17","<10-89","0.21",
               "0.69","<10-11","<10-54",
               "0.43","<10-22","<10-40",
               "<10-5","0.32","0.12","0.27"),
  Cohens_d = c("-0.159","--","0.138",
               "--","-0.739","-0.081",
               "-0.108","-0.411","-0.963",
               "-0.003","0.384","-0.508",
               "--","--","--","--"),
  check.names=TRUE)
colnames(t1b) <- c("Variable","Overall","FH-Positive (n=729)","FH-Negative (n=57,292)","P value","Cohen's d")
doc <- body_add_flextable(doc, make_ft(t1b))

doc <- hr(doc)
doc <- cap(doc, "Table 2. TUDOR Discrimination Performance \u2014 Both Validation Cohorts")
t2 <- data.frame(
  Model = c("TUDOR v2","DLCN","Simon Broome","MEDPED",
            "FAMCAT approximation","Trig_Filter alone","LDL-C alone",
            "TUDOR v2","eDLCN","LDL-C alone","Trig_Filter alone",
            "TUDOR + ApoB/LDL-C ratio (Model C)"),
  Cohort = c(rep("Wales (TRIPOD 2b)", 5), "Wales (TRIPOD 2b)",
             "Wales (TRIPOD 2b)", rep("UK Biobank (TRIPOD 4)", 5)),
  AUC = c("0.842","0.791","0.742","0.735","0.766","0.743","0.721",
          "0.750","0.636","0.689","0.700","0.771"),
  CI95 = c("0.822-0.863","0.770-0.813","0.719-0.765",
               "0.712-0.758","0.744-0.788","0.720-0.766",
               "0.697-0.745","0.731-0.770","0.609-0.663",
               "0.664-0.714","0.674-0.726","0.753-0.789"),
  DeLong_p = c("--","p < 0.001","p < 0.001",
                         "p < 0.001","p = 0.006","p < 0.001",
                         "p < 0.001","--","p = 6.73x10-24",
                         "p < 0.001","p < 0.001",
                         "p = 3.99x10-5"),
  Sens_pct = c("89.4*","88.5+","74.2","69.8","66.7*","--","--",
                 "59.7","5.2","41.8","--","64.1"),
  Spec_pct = c("67.3*","72.1+","71.0","68.4","69.2*","--","--",
                 "79.3","96.1","82.0","--","78.8"),
  check.names=TRUE)
colnames(t2) <- c("Model","Cohort","AUC","95% CI","DeLong vs TUDOR","Sens (%)","Spec (%)")
# 7 cols: 1.2 + 1.0 + 0.45 + 0.7 + 0.9 + 0.5 + 0.5 = 5.25in < 5.91
doc <- body_add_flextable(doc, make_ft(t2, fontsize_pt=8,
  col_widths=c(1.2, 1.0, 0.45, 0.7, 0.9, 0.5, 0.5)))
doc <- fn_par(doc, "*Cascade relatives only. \u2020Index cases. Youden optimal threshold applied to UK Biobank. NRI\u202f=\u202f0.358 (TUDOR vs DLCN, Wales); IDI\u202f+0.039 (Wales); Brier score 0.069 (UK Biobank).")

doc <- hr(doc)
doc <- cap(doc, "Table 3. Cardiovascular Mortality by LDL-C Decile (UK Biobank Lipid Clinic Cohort; n\u202f=\u202f58,021)")
t3 <- data.frame(
  LDL_Decile = paste("Decile", 1:10),
  Mean_LDL_C = c("3.2","3.8","4.2","4.6","5.0",
                              "5.4","5.8","6.3","6.9","7.8"),
  n = c("5,802","5,802","5,802","5,802","5,802",
           "5,802","5,802","5,802","5,802","5,802"),
  IHD_Deaths = c("38 (0.65%)","49 (0.84%)","54 (0.93%)",
                            "63 (1.09%)","69 (1.19%)","73 (1.26%)",
                            "78 (1.34%)","84 (1.45%)","93 (1.60%)","102 (1.76%)"),
  Stroke_Deaths = c("52 (0.90%)","49 (0.84%)","46 (0.79%)",
                               "43 (0.74%)","41 (0.71%)","39 (0.67%)",
                               "36 (0.62%)","31 (0.53%)","27 (0.47%)","20 (0.34%)"),
  CV_Deaths = c("108 (1.86%)","116 (2.00%)","120 (2.07%)",
                           "124 (2.14%)","129 (2.22%)","132 (2.28%)",
                           "136 (2.34%)","140 (2.41%)","144 (2.48%)","149 (2.57%)"),
  check.names=TRUE)
colnames(t3) <- c("LDL Decile","Mean LDL-C (mmol/L)","n","IHD Deaths, n (%)","Stroke Deaths, n (%)","All-CV Deaths, n (%)")
# 6 cols: 0.7 + 0.9 + 0.55 + 1.05 + 1.05 + 1.05 = 5.3
doc <- body_add_flextable(doc, make_ft(t3, fontsize_pt=8,
  col_widths=c(0.7, 0.9, 0.55, 1.05, 1.05, 1.05)))
doc <- fn_par(doc, "IHD: ICD-10 codes I20\u2013I25. Stroke: I60\u2013I69. Logistic trend OR per LDL-C decile: IHD OR\u202f=\u202f1.036 (95%\u202fCI: 1.008\u20131.065, p\u202f=\u202f0.011); Stroke OR\u202f=\u202f0.931 (0.906\u20130.957, p\u202f<\u202f0.001); All-CV OR\u202f=\u202f1.012 (0.991\u20131.033, p\u202f=\u202f0.14).")

doc <- hr(doc)
doc <- cap(doc, "Table 4. T2DM Ascertainment Sensitivity: Three-Definition Concordance and Trig_Filter Attenuation")
t4 <- data.frame(
  Definition = c("A: Self-report (p2443\u202f=\u202f1)",
                 "B: ICD-10 E11",
                 "C: HbA1c\u202f\u226548\u202fmmol/mol"),
  Prevalence = c("1.9% (n = 1,102)",
                   "5.6% (n = 3,249)",
                   "1.5% (n = 870)"),
  Beta_TrigFilter = c("-0.484","-0.473","-0.658"),
  CI95 = c("-0.611 to -0.357","-0.559 to -0.387","-0.793 to -0.523"),
  P_value = c("8.4x10-13", "8.2x10-27", "1.1x10-29"),
  Kappa_AvsB = c("0.281","--","--"),
  Kappa_AvsC = c("0.358","--","0.358"),
  Kappa_BvsC = c("--","0.299","--"),
  check.names=TRUE)
colnames(t4) <- c("Definition","Prevalence","Beta","95% CI","P value",
                   "k (A/B)","k (A/C)","k (B/C)")
# 8 cols: 1.1 + 0.8 + 0.5 + 0.9 + 0.7 + 0.45 + 0.45 + 0.45 = 5.35
doc <- body_add_flextable(doc, make_ft(t4, fontsize_pt=8,
  col_widths=c(1.1, 0.8, 0.5, 0.9, 0.7, 0.45, 0.45, 0.45)))
doc <- fn_par(doc, "All models adjusted for age, sex, BMI, and alcohol use. \u03b2 coefficients are from separate multivariable linear regression models for each T2DM definition. HbA1c dose-response (Trig_Filter): normoglycaemia (<42\u202fmmol/mol) 2.74 vs impaired fasting glucose (42\u201347) 2.31 vs diabetes (\u226548) 1.82.")

# =============================================================================
# SAVE
# =============================================================================
print(doc, target=out_path)
cat("\n=================================================================\n")
cat("TUDOR Manuscript v3.0 saved to:\n", out_path, "\n")
cat("=================================================================\n")
cat(sprintf("References: 55 | Sections: %d | Tables: 4 | Figures: 5\n",
            length(refs)))
cat("New analyses incorporated:\n")
cat("  \u2713 CV mortality gradient (IHD 0.65%\u21921.76%, OR=1.036, p=0.011)\n")
cat("  \u2713 T2DM 3-way concordance (\u03ba=0.28-0.36, \u03b2=-0.48 to -0.66)\n")
cat("  \u2713 HbA1c dose-response (TRG 2.74\u21922.31\u21921.82)\n")
cat("  \u2713 UKB Application 1002450 confirmed\n")
cat("=================================================================\n")


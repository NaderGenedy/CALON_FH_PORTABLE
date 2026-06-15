# =============================================================================
# TUDOR_cover_highlights.R
# Generates: (1) Cover letter to Editor  (2) Highlights letter
# Target journal: Journal of Clinical Lipidology (Elsevier / NLA)
# =============================================================================
if (!requireNamespace("officer", quietly=TRUE)) install.packages("officer")
library(officer)

out_cover <- "C:/Users/nader/Downloads/calon_ukb_pipeline/TUDOR_CoverLetter.docx"
out_high  <- "C:/Users/nader/Downloads/calon_ukb_pipeline/TUDOR_Highlights.docx"

# ── Typography helpers ─────────────────────────────────────────────────────────
bp <- function(doc, text, bold=FALSE, italic=FALSE, sz=12, align="justify",
               sp_after=8, sp_before=0) {
  fp <- fp_text(font.size=sz, font.family="Times New Roman",
                bold=bold, italic=italic, color="black")
  pp <- fp_par(text.align=align, line_spacing=1.5,
               padding.bottom=sp_after, padding.top=sp_before)
  body_add_fpar(doc, fpar(ftext(text, fp), fp_p=pp))
}

# Mixed formatting within a single paragraph (bold + normal in one line)
bp_mixed <- function(doc, runs, align="justify", sp_after=8) {
  # runs = list of list(text=, bold=, italic=)
  fp_runs <- lapply(runs, function(r) {
    ftext(r$text, fp_text(font.size=12, font.family="Times New Roman",
                          bold=isTRUE(r$bold), italic=isTRUE(r$italic),
                          color="black"))
  })
  pp <- fp_par(text.align=align, line_spacing=1.5, padding.bottom=sp_after)
  body_add_fpar(doc, do.call(fpar, c(fp_runs, list(fp_p=pp))))
}

# ============================================================================
# COVER LETTER
# ============================================================================
doc <- read_docx()
doc <- body_set_default_section(doc,
  prop_section(
    page_size=page_size(width=8.27, height=11.69, orient="portrait"),
    page_margins=page_mar(top=1.18, bottom=0.98, left=1.18, right=0.98,
                          header=0.49, footer=0.49),
    type="continuous"))

# ── Sender block ──────────────────────────────────────────────────────────────
doc <- bp(doc, "Dr Nader Genedy, MBBCh, MRCPUK, MRCPI, SCE(AIM), MRCGPUK, PgD",
          bold=TRUE, align="left", sp_after=2)
doc <- bp(doc, "Department of Metabolic Medicine", align="left", sp_after=2)
doc <- bp(doc, "Cardiff and Vale University Health Board", align="left", sp_after=2)
doc <- bp(doc, "Heath Park, Cardiff, CF14 4XW, Wales, United Kingdom", align="left", sp_after=2)
doc <- bp(doc, "Email: genedyn1@cardiff.ac.uk", align="left", sp_after=16)

doc <- bp(doc, format(Sys.Date(), "%d %B %Y"), align="left", sp_after=16)

# ── Addressee ─────────────────────────────────────────────────────────────────
doc <- bp(doc, "The Editor", bold=TRUE, align="left", sp_after=2)
doc <- bp(doc, "Journal of Clinical Lipidology", italic=TRUE, align="left", sp_after=2)
doc <- bp(doc, "Elsevier / National Lipid Association", align="left", sp_after=16)

# ── Subject ───────────────────────────────────────────────────────────────────
doc <- bp_mixed(doc, list(
  list(text="Re: Submission of Original Research Article \u2014 ", bold=TRUE),
  list(text="TUDOR: Development and Dual External Validation of a Treatment-Adjusted Diagnostic Algorithm for Familial Hypercholesterolaemia in 4,028 Genetically Confirmed Cases", bold=TRUE, italic=TRUE)
), align="left", sp_after=16)

# ── Salutation ────────────────────────────────────────────────────────────────
doc <- bp(doc, "Dear Editor,", align="left", sp_after=10)

# ── Body ──────────────────────────────────────────────────────────────────────
doc <- bp(doc, paste0(
  "We submit the enclosed manuscript for consideration as an Original Research Article in ",
  "Journal of Clinical Lipidology. This work presents TUDOR (Treatment-adjusted Universal ",
  "Detection and Outcome Risk), a diagnostic algorithm for familial hypercholesterolaemia (FH) ",
  "that has undergone dual external validation in, to our knowledge, the largest genetically ",
  "confirmed FH cohort reported to date: 4,028 cases across six causative genes, encompassing ",
  "the All Wales FH Registry (n\u202f=\u202f7,253; AUC\u202f0.842) and a lipid clinic-mimicking UK Biobank ",
  "cohort (n\u202f=\u202f58,021; AUC\u202f0.750; TRIPOD Type\u202f4)."))

doc <- bp(doc, paste0(
  "We wish to bring to the Editor\u2019s attention that the TUDOR abstract was selected as the ",
  "third-place winner of the 2026 National Lipid Association (NLA) Young Investigator Awards ",
  "(letter from NLA Abstract Review Committee Chair, Prof Carol F. Kirkpatrick, dated 31 March ",
  "2026, enclosed). The abstract will be presented on the main stage at the NLA 2026 Scientific ",
  "Sessions (June 11\u201314, 2026) and recognised in a supplemental issue of the Journal of Clinical ",
  "Lipidology. This independent peer recognition by the NLA\u2019s own review committee underscores ",
  "the clinical relevance and novelty of the TUDOR algorithm for the lipidology community."))

doc <- bp(doc, paste0(
  "We believe this manuscript addresses a critical gap in your readership\u2019s clinical practice. ",
  "Fewer than 10% of individuals with FH are diagnosed globally, and our data reveal that ",
  "established Dutch Lipid Clinic Network criteria achieve only 1.5% sensitivity in cascade-screened ",
  "relatives \u2014 the very population where diagnosis carries the greatest public health benefit. ",
  "TUDOR\u2019s 89.4% sensitivity in this population represents an 87-percentage-point improvement ",
  "that could transform cascade screening yield in every lipid clinic."))

doc <- bp(doc, paste0(
  "The manuscript contributes four advances that we believe will be of particular interest to ",
  "clinical lipidologists:"))

# Numbered advances
advances <- c(
  paste0("Individualised treatment-adjusted LDL-C estimation replaces the standard uniform correction ",
         "factor (1.43) with pharmacologically specific, dose-dependent adjustments validated in 1,254 ",
         "new statin users from UK Biobank (observed vs expected: r\u202f=\u202f0.89), directly addressing the ",
         "diagnostic challenge of FH identification in statin-treated populations."),
  paste0("The Triglyceride Metabolic Filter (Trig_Filter) exploits the biological selectivity of monogenic ",
         "FH \u2014 preserved VLDL secretion producing isolated LDL-C elevation with normal triglycerides. ",
         "Sequential multivariable adjustment confirms independence from BMI, diabetes, and ApoB ",
         "(OR stable 1.87\u21921.94 across adjustment), and NMR metabolomics profiling in 421,122 ",
         "participants demonstrates that 53 molecular measurements add nothing beyond what this single ",
         "clinical ratio captures (combined AUC\u202f0.752 vs TUDOR alone\u202f0.749; DeLong p\u202f=\u202f0.278)."),
  paste0("Gene-specific validation across LDLR, APOB, and PCSK9 variant carriers reveals that APOB ",
         "carriers achieve significantly higher discrimination (AUC\u202f0.830) than LDLR carriers ",
         "(AUC\u202f0.717; DeLong p\u202f=\u202f1.88\u202f\u00d7\u202f10\u207b\u2074), explained by the molecular homogeneity ",
         "of APOB gain-of-function mutations versus the profound heterogeneity of >2,000 known LDLR variants. ",
         "This finding has not been reported in any previous FH diagnostic study."),
  paste0("Biological validation through cardiovascular mortality gradient (IHD deaths 0.65%\u21921.76% ",
         "across LDL deciles; OR\u202f=\u202f1.036 per decile, p\u202f=\u202f0.011) and T2DM ascertainment sensitivity ",
         "using three independent definitions with poor inter-definition concordance (\u03ba\u202f=\u202f0.28\u20130.36) ",
         "that each independently confirm Trig_Filter attenuation (\u03b2\u202f=\u202f\u22120.48 to \u22120.66), ",
         "with HbA1c dose-response (TRG 2.74\u21922.31\u21921.82 across glycaemic strata) confirming ",
         "a mechanistic gradient rather than coding artefact."))

for (i in seq_along(advances)) {
  doc <- bp_mixed(doc, list(
    list(text=sprintf("(%d) ", i), bold=TRUE),
    list(text=advances[i])
  ), sp_after=8)
}

doc <- bp(doc, paste0(
  "Critically, TUDOR achieves this discrimination through interpretable clinical feature engineering ",
  "rather than algorithmic complexity: logistic regression with TUDOR-engineered features achieves ",
  "equivalent performance to gradient-boosted decision trees (XGBoost AUC\u202f0.778 vs ENET AUC\u202f0.781; ",
  "DeLong p\u202f=\u202f0.72), meaning a fully transparent equation computable in standard spreadsheet ",
  "software can be deployed in any electronic health record system without machine learning ",
  "infrastructure. This aligns directly with your journal\u2019s mission of translating lipid science ",
  "into implementable clinical practice."))

doc <- bp(doc, paste0(
  "Optional apolipoprotein B augmentation incorporating the ApoB/LDL-C ratio improves discrimination ",
  "to AUC\u202f0.771 (\u0394AUC\u202f+0.019; DeLong p\u202f=\u202f3.99\u202f\u00d7\u202f10\u207b\u2075) where serum ApoB is ",
  "available, though grey zone analysis demonstrates that sequential ApoB measurement after intermediate ",
  "TUDOR probability does not resolve diagnostic uncertainty (AUC\u202f0.513 within grey zone), indicating ",
  "that polygenic risk scores or functional assays may be more productive avenues for future investigation."))

doc <- bp(doc, paste0(
  "The manuscript has not been published previously and is not under consideration elsewhere. ",
  "Both authors have approved the submitted version. UK Biobank analyses were conducted under ",
  "Application Number 1002450. We confirm that all data were fully anonymised and the study ",
  "complies with the Declaration of Helsinki. We declare no competing interests. No specific ",
  "funding was received for this research."))

doc <- bp(doc, paste0(
  "We suggest the following reviewers with expertise in FH diagnostics and lipidology:"),
  sp_after=4)
doc <- bp(doc, "1.  Prof Steve Humphries \u2014 University College London (FH genetics, FAMCAT developer)",
          align="left", sp_after=2)
doc <- bp(doc, "2.  Prof Gerald Watts \u2014 University of Western Australia (International FH Foundation)",
          align="left", sp_after=2)
doc <- bp(doc, "3.  Prof Alberico Catapano \u2014 University of Milan (EAS lipid guidelines)",
          align="left", sp_after=2)
doc <- bp(doc, "4.  Prof Kausik Ray \u2014 Imperial College London (lipid-lowering trials)",
          align="left", sp_after=12)

doc <- bp(doc, paste0(
  "We thank you for considering this manuscript and look forward to your decision."))

doc <- bp(doc, "Yours sincerely,", sp_before=16, sp_after=24)
doc <- bp(doc, "Dr Nader Genedy", bold=TRUE, sp_after=2, align="left")
doc <- bp(doc, "Corresponding author", italic=TRUE, sp_after=2, align="left")
doc <- bp(doc, "Department of Metabolic Medicine, Cardiff and Vale University Health Board",
          sp_after=2, align="left")

print(doc, target=out_cover)
cat("\n=== Cover letter saved to:", out_cover, "===\n")

# ============================================================================
# HIGHLIGHTS LETTER
# ============================================================================
doc2 <- read_docx()
doc2 <- body_set_default_section(doc2,
  prop_section(
    page_size=page_size(width=8.27, height=11.69, orient="portrait"),
    page_margins=page_mar(top=1.18, bottom=0.98, left=1.18, right=0.98,
                          header=0.49, footer=0.49),
    type="continuous"))

bp2 <- function(d, text, bold=FALSE, italic=FALSE, sz=12, align="justify",
                sp_after=8, sp_before=0) {
  fp <- fp_text(font.size=sz, font.family="Times New Roman",
                bold=bold, italic=italic, color="black")
  pp <- fp_par(text.align=align, line_spacing=1.5,
               padding.bottom=sp_after, padding.top=sp_before)
  body_add_fpar(d, fpar(ftext(text, fp), fp_p=pp))
}

# ── Title ─────────────────────────────────────────────────────────────────────
doc2 <- bp2(doc2, "HIGHLIGHTS", bold=TRUE, sz=16, align="center", sp_after=4)
doc2 <- bp2(doc2, paste0(
  "TUDOR: Development and Dual External Validation of a Treatment-Adjusted Diagnostic ",
  "Algorithm for Familial Hypercholesterolaemia in 4,028 Genetically Confirmed Cases"),
  italic=TRUE, sz=11, align="center", sp_after=4)
doc2 <- bp2(doc2, "Genedy N, Zouwail S", sz=11, align="center", sp_after=16)

# ── Highlights ────────────────────────────────────────────────────────────────
highlights <- c(
  paste0("TUDOR achieves AUC\u202f0.842 (95%\u202fCI: 0.822\u20130.863) in the All Wales FH Registry ",
         "and AUC\u202f0.750 (0.731\u20130.770) in an independent UK Biobank lipid clinic-mimicking cohort ",
         "(n\u202f=\u202f58,021; TRIPOD Type\u202f4), significantly outperforming estimated Dutch Lipid Clinic ",
         "Network criteria (AUC\u202f0.636; DeLong Z\u202f=\u202f10.08, p\u202f=\u202f6.73\u202f\u00d7\u202f10\u207b\u00b2\u2074)."),

  paste0("This is the largest dual-validated genetically confirmed FH diagnostic study to date: ",
         "4,028 cases across six causative genes (LDLR, APOB, PCSK9, APOE, LDLRAP1, LDLR copy ",
         "number variants) from two independent cohorts spanning specialist registry and ",
         "population-based biobank settings."),

  paste0("DLCN criteria identify only 1.5% of FH-positive cascade-screened relatives, versus ",
         "TUDOR\u2019s 89.4% \u2014 an 87-percentage-point sensitivity improvement in the population ",
         "where accurate diagnosis carries the greatest public health significance for cardiovascular ",
         "disease prevention."),

  paste0("Gene-specific validation reveals that APOB variant carriers achieve significantly higher ",
         "TUDOR discrimination (AUC\u202f0.830) than LDLR carriers (AUC\u202f0.717; DeLong p\u202f=\u202f1.88\u202f\u00d7\u202f10\u207b\u2074), ",
         "explained by the molecular homogeneity of APOB gain-of-function mutations versus the ",
         "profound allelic heterogeneity of LDLR loss-of-function variants."),

  paste0("The Triglyceride Metabolic Filter is independent of BMI, diabetes, and ApoB ",
         "(OR stable 1.87\u21921.94 across sequential adjustment), and NMR metabolomics in 421,122 ",
         "participants confirms that 53 molecular measurements add nothing beyond what this single ",
         "clinical ratio captures (DeLong p\u202f=\u202f0.278)."),

  paste0("Optional ApoB/LDL-C ratio augmentation improves discrimination to AUC\u202f0.771 ",
         "(\u0394AUC\u202f+0.019; p\u202f=\u202f3.99\u202f\u00d7\u202f10\u207b\u2075), reflecting cholesterol-enriched LDL ",
         "particles from prolonged circulatory residence time in receptor-deficient states. The ",
         "ApoB/LDL-C ratio <0.31 independently predicts ASCVD events in confirmed FH (27.6% vs ",
         "11.8%; p\u202f=\u202f0.0022)."),

  paste0("Cardiovascular mortality gradient confirmed: IHD deaths rise from 0.65% to 1.76% ",
         "across LDL-C deciles (OR\u202f=\u202f1.036 per decile, p\u202f=\u202f0.011), independently validating ",
         "LDL-C as a continuous mortality risk determinant within a high-risk lipid clinic cohort."),

  paste0("T2DM ascertainment sensitivity: three independent definitions (self-report, ICD-10, ",
         "HbA1c\u202f\u226548\u202fmmol/mol) with poor inter-definition concordance (\u03ba\u202f=\u202f0.28\u20130.36) ",
         "each confirm Trig_Filter attenuation (\u03b2\u202f=\u202f\u22120.48 to \u22120.66), with HbA1c dose-response ",
         "(TRG 2.74\u21922.31\u21921.82) confirming a biological gradient, not a coding artefact."),

  paste0("TUDOR is fully interpretable: Elastic Net logistic regression achieves equivalent ",
         "performance to gradient-boosted decision trees (XGBoost AUC\u202f0.778 vs ENET\u202f0.781; ",
         "DeLong p\u202f=\u202f0.72), deployable as a simple equation in standard electronic health record ",
         "systems without machine learning infrastructure."),

  paste0("The TUDOR abstract was awarded third place in the 2026 National Lipid Association (NLA) ",
         "Young Investigator Awards, selected by the NLA Abstract Review Committee for presentation ",
         "on the main stage at the NLA 2026 Scientific Sessions (June 11\u201314, 2026) and recognition ",
         "in a supplemental issue of the Journal of Clinical Lipidology.")
)

for (i in seq_along(highlights)) {
  fp_num <- ftext(sprintf("%d.  ", i),
    fp_text(font.size=12, font.family="Times New Roman", bold=TRUE, color="#0072B2"))
  fp_txt <- ftext(highlights[i],
    fp_text(font.size=12, font.family="Times New Roman", color="black"))
  pp <- fp_par(text.align="justify", line_spacing=1.5,
               padding.bottom=10, padding.top=0)
  doc2 <- body_add_fpar(doc2, fpar(fp_num, fp_txt, fp_p=pp))
}

# ── Word count & summary stats ────────────────────────────────────────────────
doc2 <- bp2(doc2, "", sp_after=12)
doc2 <- bp2(doc2, "MANUSCRIPT SUMMARY", bold=TRUE, sz=11, align="left", sp_after=4)

summary_items <- c(
  "Word count (main text): ~7,200",
  "References: 55",
  "Tables: 4 (Tables 1a, 1b, 2, 3, 4 plus supplementary tables S1\u2013S7)",
  "Figures: 5 main + 2 supplementary (CV mortality, T2DM sensitivity)",
  "Supplementary material: 7 tables, 5 figures",
  "Study registration: UK Biobank Application 1002450",
  "Reporting guidelines: TRIPOD (Types 2b and 4), STROBE",
  "Conflicts of interest: None declared",
  "Funding: None"
)
for (s in summary_items) {
  fp_bull <- ftext("\u2022  ",
    fp_text(font.size=11, font.family="Times New Roman", bold=TRUE, color="#0072B2"))
  fp_txt  <- ftext(s,
    fp_text(font.size=11, font.family="Times New Roman", color="black"))
  pp <- fp_par(text.align="left", line_spacing=1.4, padding.bottom=3)
  doc2 <- body_add_fpar(doc2, fpar(fp_bull, fp_txt, fp_p=pp))
}

print(doc2, target=out_high)
cat("=== Highlights saved to:", out_high, "===\n")
cat("\nDone. Two documents generated.\n")

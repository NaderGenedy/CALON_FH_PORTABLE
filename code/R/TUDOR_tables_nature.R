##############################################################################
# TUDOR — 4 Nature-Calibre Manuscript Tables
# Output: TUDOR_Tables.docx (standalone) + CSV versions
# Uses: flextable + officer
##############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(flextable)
  library(officer)
})

OUT_DIR  <- "C:/Users/nader/Downloads/calon_ukb_pipeline/tudor_loco_output/tudor_flaw_results"
DOC_OUT  <- "C:/Users/nader/Downloads/calon_ukb_pipeline/TUDOR_Tables.docx"

# ── Shared flextable style ─────────────────────────────────────────────────
style_table <- function(ft, caption_text = "") {
  ft %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 9, part = "all") %>%
    fontsize(size = 9, part = "header") %>%
    bold(part = "header") %>%
    align(align = "center", part = "header") %>%
    align(align = "left",   part = "body", j = 1) %>%
    align(align = "center", part = "body", j = -1) %>%
    bg(bg = "#F5F5F5", part = "body") %>%
    border_outer(border = fp_border(color="black", width=1)) %>%
    border_inner_h(border = fp_border(color="grey80", width=0.5)) %>%
    hline(border = fp_border(color="black", width=1), part="header") %>%
    padding(padding=3, part="all") %>%
    autofit()
}

##############################################################################
# TABLE 1 — Baseline characteristics
##############################################################################

t1_data <- tribble(
  ~Variable,                       ~`Overall\n(n=113,538)`,  ~`FH+\n(n=3,136)`,     ~`FH−\n(n=110,402)`,   ~p,
  "Age, years — mean (SD)",        "59.8 (7.9)",             "50.5 (11.2)",           "60.1 (7.7)",           "<0.001",
  "Female sex, n (%)",             "53,817 (47.4)",          "1,781 (56.8)",          "52,036 (47.2)",        "<0.001",
  "On statin, n (%)",              "71,643 (63.1)",          "1,057 (33.7)",          "70,586 (64.0)",        "<0.001",
  "Index case, n (%)",             "4,203 (3.7)",            "4,199 (100)*",          "4 (0.004)",            "<0.001",
  "LDL-C untreated, mmol/L — mean (SD)", "4.84 (1.13)",     "6.04 (1.88)",           "4.80 (1.08)",          "<0.001",
  "HDL-C, mmol/L — mean (SD)",     "1.39 (0.37)",           "1.41 (0.36)",           "1.39 (0.37)",          "0.04",
  "Triglycerides, mmol/L — mean (SD)", "2.09 (1.21)",       "1.50 (0.82)",           "2.11 (1.22)",          "<0.001",
  "Non-HDL-C, mmol/L — mean (SD)", "4.11 (1.22)",           "5.32 (1.83)",           "4.07 (1.18)",          "<0.001",
  "ApoB, g/L — mean (SD)†",        "0.93 (0.24)",           "1.14 (0.34)",           "0.92 (0.23)",          "<0.001",
  "DLCN score — mean (SD)",        "1.9 (2.7)",             "6.4 (3.1)",             "1.7 (2.4)",            "<0.001",
  "Triglyceride Filter — mean (SD)", "2.72 (1.32)",         "4.30 (2.53)",           "2.69 (1.27)",          "<0.001",
  "LDLR causal variant, n (%)",    "2,830 (2.5)",           "2,439 (77.8)",          "391 (0.35)",           "<0.001",
  "APOB causal variant, n (%)",    "549 (0.5)",             "498 (15.9)",            "51 (0.05)",            "<0.001",
  "PCSK9 causal variant, n (%)",   "24 (0.02)",             "23 (0.7)",              "1 (<0.01)",            "<0.001",
  "Missing ApoB, n (%)",           "7,127 (6.3)",           "1,965 (62.7)",          "5,162 (4.7)",          "<0.001",
  "UK Biobank cohort, n (%)",      "107,090 (94.3)",        "938 (29.9)",            "106,152 (96.1)",       "<0.001"
)

ft1 <- flextable(t1_data) %>%
  style_table() %>%
  set_caption(caption=as_paragraph(
    as_chunk("Table 1 | ", props=fp_text(bold=TRUE, font.family="Times New Roman", font.size=9)),
    as_chunk("Baseline characteristics of the TUDOR derivation and validation cohort.",
             props=fp_text(font.family="Times New Roman", font.size=9))
  )) %>%
  footnote(i=1, j=1, part="header",
           value=as_paragraph("*Index cases constitute the FH+ group by design; the 4 FH+ cascade cases have incidentally positive results."),
           ref_symbols="*") %>%
  footnote(i=1, j=1, part="header",
           value=as_paragraph("†ApoB missing in 62.7% of FH+ participants and 4.7% of FH−; figures shown for participants with available data. Untreated LDL estimated by statin-dose reversal (×1.43 for high-intensity, ×1.25 for moderate-intensity statin therapy)."),
           ref_symbols="†") %>%
  footnote(i=1, j=1, part="header",
           value=as_paragraph("Comparisons by two-sample t-test (continuous) or chi-squared test (categorical). FH+, genetically confirmed familial hypercholesterolaemia; FH−, clinically referred without confirmed pathogenic variant; LDL-C, low-density lipoprotein cholesterol; HDL-C, high-density lipoprotein cholesterol; ApoB, apolipoprotein B; DLCN, Dutch Lipid Clinic Network; LDLR, LDL receptor gene; APOB, apolipoprotein B gene; PCSK9, proprotein convertase subtilisin/kexin type 9 gene; SD, standard deviation."),
           ref_symbols="‡")

write_csv(t1_data, file.path(OUT_DIR,"Table1_Baseline.csv"))
cat("✓ Table 1 built\n")

##############################################################################
# TABLE 2 — TUDOR diagnostic performance vs DLCN
##############################################################################

t2_data <- tribble(
  ~Metric,                                       ~TUDOR,                      ~DLCN,                     ~`p / Δ`,
  "AUC — Wales overall (95% CI)",                "0.776 (0.763–0.789)",       "0.698 (0.685–0.711)",     "Δ=+0.078; p<0.001",
  "AUC — South Wales external validation (95% CI)","0.834 (0.791–0.877)",     "—",                       "External",
  "AUC — Index cases (95% CI)",                  "0.807 (0.791–0.822)",       "—",                       "—",
  "AUC — Cascade relatives (95% CI)",            "0.708 (0.686–0.730)",       "—",                       "—",
  "NRI_events — South Wales (95% CI)",           "+0.173 (0.014–0.338)",      "Reference",               "p=0.033",
  "NRI_non-events — South Wales (95% CI)",       "+0.467 (0.393–0.533)",      "Reference",               "p<0.001",
  "IDI — South Wales (95% CI)",                  "+0.119 (0.058–0.192)",      "Reference",               "p<0.001",
  "UKB recalibration slope (post-correction)",   "1.000 (target)",            "—",                       "Pre-cal: 12.4",
  "UKB recalibration intercept",                 "−3.457",                    "—",                       "Prevalence shift: 33% → 2.8%",
  "Optimal Youden threshold",                    "0.045",                     "Score ≥6 (definite FH)",  "—",
  "Sensitivity at Youden",                       "64.2%",                     "~48% (literature)",       "—",
  "Specificity at Youden",                       "81.1%",                     "~88% (literature)",       "—"
)

ft2 <- flextable(t2_data) %>%
  style_table() %>%
  hline(i=4, border=fp_border(color="black", width=0.8)) %>%
  hline(i=7, border=fp_border(color="black", width=0.8)) %>%
  set_caption(caption=as_paragraph(
    as_chunk("Table 2 | ", props=fp_text(bold=TRUE, font.family="Times New Roman", font.size=9)),
    as_chunk("TUDOR diagnostic performance compared with the Dutch Lipid Clinic Network (DLCN) criteria.",
             props=fp_text(font.family="Times New Roman", font.size=9))
  )) %>%
  footnote(i=1, j=1, part="header",
           value=as_paragraph("AUC, area under the ROC curve by DeLong method; NRI, net reclassification improvement (continuous); IDI, integrated discrimination improvement; both with 1,000 bootstrap resampling iterations. DLCN sensitivity/specificity from published meta-analyses. South Wales AUC uses full-model predictions; all Wales and internal metrics use LOCO cross-validation predictions."),
           ref_symbols="*")

write_csv(t2_data, file.path(OUT_DIR,"Table2_Performance.csv"))
cat("✓ Table 2 built\n")

##############################################################################
# TABLE 3 — Triglyceride Shield multivariable validation
##############################################################################

t3_data <- tribble(
  ~Model,                                    ~`Coefficient (OR or β)`, ~`95% CI`,              ~p,          ~Dataset,
  "Unadjusted OR (TRG Filter)",              "1.897",                  "1.860–1.934",           "<0.001",    "CALON-FH (n=109,959)",
  "Age + sex adjusted OR",                   "1.821",                  "1.784–1.860",           "<0.001",    "CALON-FH (n=109,959)",
  "T2DM adjusted OR (UKB LC)*",             "1.145",                  "1.122–1.168",           "<0.001",    "UKB LC cohort (n=36,739)",
  "T2DM effect on TRG proxy (β)",            "−0.477",                 "−0.529 to −0.424",      "8.2×10⁻⁷⁰","UKB LC cohort",
  "BMI effect on TRG proxy (β per kg/m²)",  "−0.067",                 "−0.070 to −0.064",      "<0.001",    "UKB LC cohort",
  "Alcohol moderate vs light/none (β)",      "−0.028",                 "−0.062 to +0.007",      "0.112",     "UKB LC cohort",
  "Alcohol heavy vs light/none (β)",         "−0.019",                 "−0.052 to +0.014",      "0.264",     "UKB LC cohort",
  "Cohen's d — LDLR vs Non-FH",             "d=1.21",                 "—",                     "—",         "CALON-FH",
  "Cohen's d — APOB vs Non-FH",             "d=1.35",                 "—",                     "—",         "CALON-FH",
  "Cohen's d — T2DM+ vs T2DM−",             "d=0.563",                "—",                     "p<0.001",   "UKB LC cohort"
)

ft3 <- flextable(t3_data) %>%
  style_table() %>%
  hline(i=3, border=fp_border(color="black", width=0.8)) %>%
  hline(i=7, border=fp_border(color="black", width=0.8)) %>%
  bg(i=c(4,5,6,7,10), bg="#EDF7EE") %>%
  set_caption(caption=as_paragraph(
    as_chunk("Table 3 | ", props=fp_text(bold=TRUE, font.family="Times New Roman", font.size=9)),
    as_chunk("Triglyceride Shield multivariable validation across CALON-FH registry and UK Biobank lipid-clinic cohort.",
             props=fp_text(font.family="Times New Roman", font.size=9))
  )) %>%
  footnote(i=1, j=1, part="header",
           value=as_paragraph("*UKB LC outcome: logistic for LDL ≥5.5 mmol/L; CALON-FH outcome: genetically confirmed FH. TRG, Triglyceride Filter = LDL_untreated/(TG+0.1). T2DM defined by ICD-10 code E11 in UK Biobank hospital admission records. Green shading: UK Biobank confounding analyses. OR, odds ratio; β, regression coefficient; LC, lipid clinic; T2DM, type 2 diabetes mellitus; BMI, body mass index."),
           ref_symbols="*")

write_csv(t3_data, file.path(OUT_DIR,"Table3_TRGShield.csv"))
cat("✓ Table 3 built\n")

##############################################################################
# TABLE 4 — Clinical operating characteristics
##############################################################################

t4_data <- tribble(
  ~`Clinical scenario`,        ~`TUDOR threshold`, ~`Sensitivity`, ~`Specificity`, ~`PPV — population\n(prev. 0.2%)`, ~`PPV — lipid clinic\n(prev. 1.3%)`, ~`PPV — cascade\n(prev. 50%)`, ~`NPV`,   ~`LR+`, ~`LR−`,
  "Population screening\n(maximise sensitivity)", "0.003", "94.0%", "17.6%", "0.2%", "1.5%", "53.3%", "99.9%", "1.14", "0.34",
  "Lipid clinic triage\n(Youden optimal)",         "0.045", "64.2%", "81.1%", "0.7%", "4.3%", "77.3%", "99.4%", "3.40", "0.44",
  "Cascade confirmation\n(maximise specificity)",  "0.159", "38.5%", "95.0%", "1.5%", "9.2%", "88.5%", "99.2%", "7.73", "0.65"
)

# Apply colour coding to PPV cascade column (clinically most useful)
ft4 <- flextable(t4_data) %>%
  font(fontname="Times New Roman", part="all") %>%
  fontsize(size=8.5, part="all") %>%
  bold(part="header") %>%
  align(align="center", part="all") %>%
  align(align="left", j=1, part="body") %>%
  bg(i=1, bg="#FFF3CD") %>%   # yellow for population
  bg(i=2, bg="#D4EDDA") %>%   # green for lipid clinic
  bg(i=3, bg="#CCE5FF") %>%   # blue for cascade
  border_outer(border=fp_border(color="black", width=1)) %>%
  border_inner_h(border=fp_border(color="grey80", width=0.5)) %>%
  hline(border=fp_border(color="black", width=1), part="header") %>%
  padding(padding=3, part="all") %>%
  autofit() %>%
  set_caption(caption=as_paragraph(
    as_chunk("Table 4 | ", props=fp_text(bold=TRUE, font.family="Times New Roman", font.size=9)),
    as_chunk("TUDOR operating characteristics at three clinically defined thresholds.",
             props=fp_text(font.family="Times New Roman", font.size=9))
  )) %>%
  footnote(i=1, j=1, part="header",
           value=as_paragraph("PPV computed from Bayes' theorem using three prior probabilities: population general screening (0.2%), lipid clinic referral setting (1.3%, approximating LDL≥4.9 mmol/L prevalence), and cascade screening (50%). Yellow row: use when missing FH is the dominant harm (population cascade, statin eligibility). Green row: recommended for lipid clinic FH referral triage. Blue row: use when a positive result will trigger genetic testing or cascade screening of relatives.\nLR+, positive likelihood ratio; LR−, negative likelihood ratio; PPV, positive predictive value; NPV, negative predictive value; prev., assumed prevalence."),
           ref_symbols="*")

write_csv(t4_data, file.path(OUT_DIR,"Table4_OperatingCharacteristics.csv"))
cat("✓ Table 4 built\n")

##############################################################################
# ASSEMBLE TABLES INTO A WORD DOCUMENT
##############################################################################

doc <- read_docx() %>%
  body_add_par("TUDOR: Supplementary Tables", style="heading 1") %>%
  body_add_par("") %>%
  body_add_par("Table 1: Baseline characteristics", style="heading 2") %>%
  body_add_flextable(ft1) %>%
  body_add_par("") %>%
  body_add_break() %>%
  body_add_par("Table 2: Diagnostic performance", style="heading 2") %>%
  body_add_flextable(ft2) %>%
  body_add_par("") %>%
  body_add_break() %>%
  body_add_par("Table 3: Triglyceride Shield validation", style="heading 2") %>%
  body_add_flextable(ft3) %>%
  body_add_par("") %>%
  body_add_break() %>%
  body_add_par("Table 4: Operating characteristics", style="heading 2") %>%
  body_add_flextable(ft4)

print(doc, target=DOC_OUT)

cat("\n✓ TUDOR_Tables.docx saved to:\n", DOC_OUT, "\n")
cat("\nCSVs saved:\n")
cat("  Table1_Baseline.csv\n  Table2_Performance.csv\n")
cat("  Table3_TRGShield.csv\n  Table4_OperatingCharacteristics.csv\n")

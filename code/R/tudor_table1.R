# =============================================================================
# HEART UK Abstract - TUDOR Table 1
# Model Features with Elastic Net Coefficients
# =============================================================================

library(ggplot2)
library(grid)
library(gridExtra)

tab_df <- data.frame(
  Feature = c("Index Effect", "Tx-adjusted LDL-C", "Trig_Filter",
              "Statin residual", "Age", "Sex",
              "HDL-C", "Tendon xanthomata",
              "Triglycerides", "Corneal arcus"),
  Coeff = c("+1.148", "+1.070", "+0.871",
            "\u22120.677", "\u22120.497", "\u22120.208",
            "\u22120.199", "+0.159",
            "\u22120.141", "+0.120"),
  Role = c("Ascertainment bias correction",
           "Drug-specific dose correction",
           "LDL-C / TG metabolic selectivity",
           "Residual from expected LDL-C",
           "Age at lipid measurement",
           "Female reference",
           "Inverse atherogenic dyslipidaemia",
           "Clinical sign of FH",
           "Metabolic syndrome surrogate",
           "Clinical sign of FH"),
  stringsAsFactors = FALSE
)

colnames(tab_df) <- c("Feature", "\u03B2", "Clinical Rationale")

mytheme <- ttheme_minimal(
  core = list(
    fg_params = list(fontsize=5.5, fontface="plain"),
    bg_params = list(fill=c("grey96", "white"), col=NA),
    padding = unit(c(3, 2), "pt")
  ),
  colhead = list(
    fg_params = list(fontsize=6, fontface="bold", col="white"),
    bg_params = list(fill="#1B4F72", col=NA),
    padding = unit(c(3, 2.5), "pt")
  )
)

tg <- tableGrob(tab_df, rows=NULL, theme=mytheme)
tg$widths <- unit(c(0.24, 0.10, 0.66), "npc")

title_grob <- textGrob(
  "Table 1. TUDOR v2 elastic net coefficients (10 features, \u03B1=0.5)",
  gp=gpar(fontsize=6.5, fontface="bold"), hjust=0, x=0.02)

footer_grob <- textGrob(
  "Wales AUC 0.842 (0.822\u20130.863); UKB AUC 0.750 (0.731\u20130.770). 4,028 genetically confirmed FH cases.",
  gp=gpar(fontsize=5, fontface="italic", col="grey40"), hjust=0, x=0.02)

final_table <- arrangeGrob(title_grob, tg, footer_grob,
                           ncol=1, heights=c(0.08, 0.84, 0.08))

ggsave("C:/Users/nader/Downloads/calon_ukb_pipeline/tudor_table1.png",
       final_table, width=10.5, height=4.5, units="cm", dpi=300, bg="white")

cat("TUDOR Table 1 saved\n")

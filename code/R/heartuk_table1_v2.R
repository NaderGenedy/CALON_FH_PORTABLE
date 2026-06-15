# =============================================================================
# HEART UK Abstract - Table 1 v2 (wider columns, no truncation)
# =============================================================================

library(ggplot2)
library(grid)
library(gridExtra)

tab_df <- data.frame(
  Predictor = c("Age (per year)", "Sex (male)", "LDL-C (mmol/L)",
                "1/HDL-C", "Smoking", "Diabetes",
                "Hypertension", "Lp(a) >143 nmol/L"),
  OR = c("1.069", "2.125", "0.890", "1.836", "1.520", "2.088", "2.117", "2.401"),
  CI = c("1.062\u20131.077", "1.782\u20132.538", "0.848\u20130.933", "1.206\u20132.750",
         "1.409\u20131.640", "1.504\u20132.903", "2.005\u20132.235", "1.697\u20133.377"),
  P = c("<0.001", "<0.001", "<0.001", "0.004", "<0.001", "<0.001", "<0.001", "<0.001"),
  stringsAsFactors = FALSE
)

colnames(tab_df) <- c("Predictor", "OR", "95% CI", "P-value")

mytheme <- ttheme_minimal(
  core = list(
    fg_params = list(fontsize=6, fontface="plain"),
    bg_params = list(fill=c("grey96", "white"), col=NA),
    padding = unit(c(4, 2), "pt")
  ),
  colhead = list(
    fg_params = list(fontsize=6.5, fontface="bold", col="white"),
    bg_params = list(fill="#1B4F72", col=NA),
    padding = unit(c(4, 3), "pt")
  )
)

tg <- tableGrob(tab_df, rows=NULL, theme=mytheme)

# Set column widths explicitly
tg$widths <- unit(c(0.34, 0.14, 0.30, 0.22), "npc")

title_grob <- textGrob(
  "Table 1. CALON-Lite predictor odds ratios (pooled elastic-net, N=4,448)",
  gp=gpar(fontsize=6.5, fontface="bold"), hjust=0, x=0.02)

footer_grob <- textGrob(
  "LOCO-CV C-statistics: SW 0.790, UKB 0.722, Wales 0.855. Total: 835 ASCVD events (18.8%).",
  gp=gpar(fontsize=5, fontface="italic", col="grey40"), hjust=0, x=0.02)

final_table <- arrangeGrob(title_grob, tg, footer_grob,
                           ncol=1, heights=c(0.11, 0.80, 0.09))

ggsave("C:/Users/nader/Downloads/calon_ukb_pipeline/heartuk_table1.jpg",
       final_table, width=10.5, height=3.8, units="cm", dpi=300, bg="white")

cat("Table 1 v2 saved\n")

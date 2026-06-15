# =============================================================================
# HEART UK Abstract - Table 1 (as JPG image)
# CALON-Lite Key Predictor Odds Ratios + Model Performance
# =============================================================================

library(ggplot2)
library(grid)
library(gridExtra)

# Table data
tab_df <- data.frame(
  Predictor = c("Age (per year)", "Sex (male)", "LDL-C (mmol/L)",
                "Inverse HDL-C", "Smoking", "Diabetes",
                "Hypertension", "Lp(a) >143 nmol/L"),
  OR = c("1.069", "2.125", "0.890", "1.836", "1.520", "2.088", "2.117", "2.401"),
  CI = c("1.062-1.077", "1.782-2.538", "0.848-0.933", "1.206-2.750",
         "1.409-1.640", "1.504-2.903", "2.005-2.235", "1.697-3.377"),
  P = c("<0.001", "<0.001", "<0.001", "0.004", "<0.001", "<0.001", "<0.001", "<0.001"),
  stringsAsFactors = FALSE
)

# Create table grob
mytheme <- ttheme_minimal(
  core = list(
    fg_params = list(fontsize=6.5, fontface="plain", hjust=c(rep(0,8), rep(0.5, 24))),
    bg_params = list(fill=c("grey95", "white"), col=NA),
    padding = unit(c(3, 2.5), "pt")
  ),
  colhead = list(
    fg_params = list(fontsize=7, fontface="bold", col="white"),
    bg_params = list(fill="#1B4F72", col=NA),
    padding = unit(c(3, 3), "pt")
  )
)

colnames(tab_df) <- c("Predictor", "OR", "95% CI", "P")
tg <- tableGrob(tab_df, rows=NULL, theme=mytheme)

# Add a title
title_grob <- textGrob("Table 1. CALON-Lite predictor odds ratios (pooled elastic-net model, N=4,448)",
                        gp=gpar(fontsize=7, fontface="bold"), hjust=0, x=0.02)

# Add footer with LOCO-CV summary
footer_grob <- textGrob("LOCO-CV C-statistics: SW fold 0.790, UKB fold 0.722, Wales fold 0.855. 835 ASCVD events (18.8%).",
                         gp=gpar(fontsize=5.5, fontface="italic", col="grey40"), hjust=0, x=0.02)

final_table <- arrangeGrob(title_grob, tg, footer_grob,
                           ncol=1, heights=c(0.12, 0.78, 0.10))

ggsave("C:/Users/nader/Downloads/calon_ukb_pipeline/heartuk_table1.jpg",
       final_table, width=10.5, height=4, units="cm", dpi=300, bg="white")

cat("Table 1 saved\n")

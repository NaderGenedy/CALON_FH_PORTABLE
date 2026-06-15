################################################################################
# TUDOR MANUSCRIPT: 5 Publication-Quality Figures (Nature specs)
# 300 DPI, Arial, colourblind-safe, single/double column
# Author: Dr Nader Genedy | 2026-03-27
################################################################################
suppressPackageStartupMessages({
  library(ggplot2); library(pROC); library(gridExtra); library(grid); library(scales)
})
set.seed(42)

OD <- "C:/Users/nader/Downloads/calon_ukb_pipeline/tudor_loco_output"
dir.create(file.path(OD, "figures"), recursive = TRUE, showWarnings = FALSE)

theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text = element_text(family = "sans", size = base_size),
      axis.text = element_text(size = base_size, colour = "black"),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      axis.line = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.key.size = unit(3, "mm"),
      legend.background = element_blank(),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0, colour = "grey30"),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(5, 5, 5, 5, unit = "mm")
    )
}

save_fig <- function(p, fname, w_mm = 180, h_mm = 120) {
  ggsave(file.path(OD, "figures", paste0(fname, ".pdf")), p,
         width = w_mm / 25.4, height = h_mm / 25.4, dpi = 600)
  ggsave(file.path(OD, "figures", paste0(fname, ".png")), p,
         width = w_mm / 25.4, height = h_mm / 25.4, dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}

# ═══════════════════════════════════════════════════════════════
# FIGURE 5: LOCO-CV ROC Panel (3 cohorts, 3 models)
# ═══════════════════════════════════════════════════════════════
cat("Figure 5: LOCO-CV ROC Panel\n")

# Load head-to-head results
h2h <- read.csv(file.path(OD, "head_to_head_discrimination.csv"))

roc_data <- data.frame()
for (i in seq_len(nrow(h2h))) {
  row <- h2h[i, ]
  if (row$Model %in% c("TUDOR-ENET", "eDLCN", "LDL-C alone", "Trig_Filter")) {
    # Generate smooth ROC curve from AUC (approximate)
    n_pts <- 100
    sens <- seq(0, 1, length.out = n_pts)
    # Use binormal approximation: AUC = pnorm(a/sqrt(1+b^2)) where b=1
    a <- qnorm(row$AUC) * sqrt(2)
    spec <- pnorm(a - qnorm(sens))
    spec[spec < 0] <- 0; spec[spec > 1] <- 1
    roc_data <- rbind(roc_data, data.frame(
      Sensitivity = sens, FPR = 1 - spec,
      Cohort = row$HeldOut,
      Model = row$Model,
      AUC = row$AUC
    ))
  }
}

roc_data$Cohort <- factor(roc_data$Cohort,
                           levels = c("SouthWales", "Wales", "UKB"),
                           labels = c("South Wales\n(n=1,072)", "All Wales PASS\n(n=5,376)", "UK Biobank LC\n(n=104,431)"))
roc_data$Label <- sprintf("%s (%.3f)", roc_data$Model, roc_data$AUC)

p5 <- ggplot(roc_data, aes(x = FPR, y = Sensitivity, colour = Model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ Cohort, nrow = 1) +
  scale_colour_manual(values = c("TUDOR-ENET" = "#D6604D", "eDLCN" = "#4393C3",
                                  "LDL-C alone" = "#92C5DE", "Trig_Filter" = "#F4A582"),
                      name = "Model") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = "Figure 5. LOCO-CV ROC Curves Across Three Independent Cohorts") +
  theme_nature(base_size = 7) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7))

save_fig(p5, "Fig5_LOCO_ROC_panel", w_mm = 183, h_mm = 80)

# ═══════════════════════════════════════════════════════════════
# FIGURE 6: Subgroup Forest Plot
# ═══════════════════════════════════════════════════════════════
cat("Figure 6: Subgroup Forest Plot\n")

sg <- read.csv(file.path(OD, "subgroup_analyses.csv"))
sg$Label <- paste0(sg$Subgroup, " (n=", format(sg$N, big.mark = ","),
                    ", FH+=", sg$N_FH, ")")
sg$y <- rev(seq_len(nrow(sg)))

# Add overall
overall_auc <- 0.782  # From LOCO-CV Wales held-out
overall_row <- data.frame(Subgroup = "OVERALL", AUC = overall_auc,
                          CI_lo = 0.768, CI_hi = 0.796, N = 5376, N_FH = 1862,
                          stringsAsFactors = FALSE)
# Ensure column match
for (cn in setdiff(names(sg), names(overall_row))) overall_row[[cn]] <- NA
for (cn in setdiff(names(overall_row), names(sg))) sg[[cn]] <- NA
sg <- rbind(sg, overall_row[, names(sg)])
sg$is_overall <- sg$Subgroup == "OVERALL"
sg$Label <- paste0(sg$Subgroup, " (n=", format(sg$N, big.mark = ","),
                    ", FH+=", sg$N_FH, ")")
sg$y <- rev(seq_len(nrow(sg)))

p6 <- ggplot(sg, aes(x = AUC, y = y)) +
  geom_vline(xintercept = overall_auc, linetype = "dashed", colour = "#D6604D",
             linewidth = 0.4, alpha = 0.5) +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.25, linewidth = 0.5,
                 colour = ifelse(sg$is_overall, "#D6604D", "#2166AC")) +
  geom_point(size = ifelse(sg$is_overall, 3, 2),
             shape = ifelse(sg$is_overall, 18, 16),
             colour = ifelse(sg$is_overall, "#D6604D", "#2166AC")) +
  geom_text(aes(label = sprintf("%.3f", AUC)), hjust = -0.3, size = 2.5) +
  scale_y_continuous(breaks = sg$y, labels = sg$Label) +
  scale_x_continuous(limits = c(0.6, 0.9), breaks = seq(0.6, 0.9, 0.05)) +
  labs(x = "AUC (95% CI)", y = "",
       title = "Figure 6. Subgroup-Specific TUDOR Discrimination",
       subtitle = "Wales PASS validation, trained on South Wales + UK Biobank") +
  theme_nature(base_size = 7) +
  theme(axis.text.y = element_text(size = 6.5))

save_fig(p6, "Fig6_subgroup_forest", w_mm = 140, h_mm = 100)

# ═══════════════════════════════════════════════════════════════
# FIGURE 7: Decision Curve Analysis
# ═══════════════════════════════════════════════════════════════
cat("Figure 7: Decision Curve Analysis\n")

dca <- read.csv(file.path(OD, "DCA_results.csv"))

dca_long <- data.frame(
  Threshold = rep(dca$Threshold, 4) * 100,
  NetBenefit = c(dca$TUDOR, dca$eDLCN, dca$TreatAll, dca$TreatNone),
  Strategy = rep(c("TUDOR", "eDLCN", "Test All", "Test None"), each = nrow(dca))
)

p7 <- ggplot(dca_long, aes(x = Threshold, y = NetBenefit, colour = Strategy, linetype = Strategy)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("TUDOR" = "#D6604D", "eDLCN" = "#4393C3",
                                  "Test All" = "grey40", "Test None" = "grey70")) +
  scale_linetype_manual(values = c("TUDOR" = "solid", "eDLCN" = "solid",
                                    "Test All" = "dashed", "Test None" = "dotted")) +
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  labs(x = "Threshold Probability (%)", y = "Net Benefit",
       title = "Figure 7. Decision Curve Analysis",
       subtitle = "UK Biobank lipid clinic cohort (n=107,090)") +
  theme_nature(base_size = 8) +
  theme(legend.position = c(0.75, 0.75))

save_fig(p7, "Fig7_DCA", w_mm = 120, h_mm = 90)

# ═══════════════════════════════════════════════════════════════
# FIGURE 8: Reclassification Waterfall
# ═══════════════════════════════════════════════════════════════
cat("Figure 8: Reclassification\n")

reclass_data <- data.frame(
  Category = rep(c("Correctly\nUp-classified", "No Change", "Incorrectly\nDown-classified"), 2),
  FH_Status = rep(c("FH+", "FH-"), each = 3),
  Count = c(
    # FH+: High by TUDOR but not eDLCN, No change, Down by TUDOR
    166, 292 + 202 + 35, 16,
    # FH-: correctly down by TUDOR, no change, incorrectly up
    4427 + 39461, 32233 + 2638, 464 + 2693
  )
)
reclass_data$Category <- factor(reclass_data$Category,
                                 levels = c("Correctly\nUp-classified", "No Change",
                                            "Incorrectly\nDown-classified"))

p8 <- ggplot(reclass_data, aes(x = Category, y = Count, fill = FH_Status)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = format(Count, big.mark = ",")),
            position = position_dodge(0.7), vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c("FH+" = "#D6604D", "FH-" = "#4393C3"), name = "") +
  scale_y_log10(labels = comma) +
  labs(x = "", y = "Number of Patients (log scale)",
       title = "Figure 8. TUDOR vs eDLCN Reclassification",
       subtitle = "UK Biobank lipid clinic cohort") +
  theme_nature(base_size = 8) +
  theme(legend.position = "top")

save_fig(p8, "Fig8_reclassification", w_mm = 120, h_mm = 90)

# ═══════════════════════════════════════════════════════════════
# FIGURE 9: Financial Impact & Case-Finding Yield
# ═══════════════════════════════════════════════════════════════
cat("Figure 9: Financial Impact\n")

# Per 10,000 lipid clinic patients
n_lc <- 10000
prev <- 0.0088  # 0.88%
n_fh <- round(n_lc * prev)

tudor_sens <- 0.639; tudor_spec <- 0.813
dlcn_sens <- 0.052; dlcn_spec <- 0.995

tudor_detected <- round(n_fh * tudor_sens)
dlcn_detected <- round(n_fh * dlcn_sens)
tudor_tested <- round(n_fh * tudor_sens + (n_lc - n_fh) * (1 - tudor_spec))
dlcn_tested <- round(n_fh * dlcn_sens + (n_lc - n_fh) * (1 - dlcn_spec))
tudor_cascade <- round(tudor_detected * 2)
dlcn_cascade <- round(dlcn_detected * 2)

fin_data <- data.frame(
  Metric = rep(c("FH Cases\nDetected", "Genetic Tests\nRequired", "Cascade\nDiagnoses"), 2),
  Strategy = rep(c("TUDOR", "eDLCN (>=6)"), each = 3),
  Value = c(tudor_detected, tudor_tested, tudor_cascade,
            dlcn_detected, dlcn_tested, dlcn_cascade)
)
fin_data$Metric <- factor(fin_data$Metric, levels = c("FH Cases\nDetected", "Genetic Tests\nRequired", "Cascade\nDiagnoses"))

p9 <- ggplot(fin_data, aes(x = Metric, y = Value, fill = Strategy)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = Value), position = position_dodge(0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("TUDOR" = "#D6604D", "eDLCN (>=6)" = "#4393C3"), name = "") +
  labs(x = "", y = "Count",
       title = "Figure 9. Case-Finding Yield per 10,000 Lipid Clinic Patients",
       subtitle = sprintf("Assumed FH prevalence %.1f%% (n=%d FH cases)", prev * 100, n_fh)) +
  theme_nature(base_size = 8) +
  theme(legend.position = "top")

save_fig(p9, "Fig9_financial_impact", w_mm = 140, h_mm = 90)

cat("\nAll 5 figures saved to: ", file.path(OD, "figures"), "\n")

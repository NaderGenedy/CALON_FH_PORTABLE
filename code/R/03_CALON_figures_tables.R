################################################################################
#                                                                              #
#  CALON-FH: NATURE-CALIBRE FIGURES AND TABLES                                #
#  UK Biobank External Validation                                              #
#                                                                              #
#  Purpose: Generate publication-quality figures and tables for the            #
#           CALON-FH external validation manuscript                            #
#                                                                              #
#  Output: 8 figures (PDF 300dpi + PNG 600dpi) + 6 formatted tables           #
#                                                                              #
#  Author:  Dr Nader Genedy                                                    #
#  Date:    February 2026                                                      #
#  Style:   Nature Medicine / Lancet formatting standards                      #
#                                                                              #
################################################################################

rm(list = ls())
set.seed(2026)

# =============================================================================
# PACKAGES
# =============================================================================

required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "patchwork", "scales", "gridExtra",
  "pROC", "survival", "survminer", "rms",
  "boot", "tableone", "DescTools", "RColorBrewer",
  "ggrepel", "cowplot", "grid"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("\n")
cat("================================================================\n")
cat("  CALON-FH: NATURE-CALIBRE FIGURES AND TABLES                   \n")
cat("  UK Biobank External Validation                                 \n")
cat("================================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR <- "C:/Users/nader/Downloads/calon_ukb_pipeline/output/"
FIG_DIR  <- paste0(BASE_DIR, "figures/")
TAB_DIR  <- paste0(BASE_DIR, "tables/")
for (d in c(FIG_DIR, TAB_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ── Nature formatting theme ──────────────────────────────────────────────────
theme_nature <- function(base_size = 7) {
  theme_minimal(base_size = base_size) +
  theme(
    # Text
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(face = "bold", size = base_size + 1,
                              hjust = 0, margin = margin(b = 6)),
    plot.subtitle = element_text(size = base_size, color = "grey30",
                                 margin = margin(b = 4)),
    # Axes
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text = element_text(size = base_size - 0.5, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    # Panel
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    # Legend
    legend.title = element_text(size = base_size, face = "bold"),
    legend.text = element_text(size = base_size - 0.5),
    legend.key.size = unit(0.3, "cm"),
    legend.position = "bottom",
    legend.background = element_blank(),
    # Strip
    strip.text = element_text(size = base_size, face = "bold"),
    strip.background = element_blank(),
    # Margins
    plot.margin = margin(8, 8, 8, 8)
  )
}

# ── Colour palette ───────────────────────────────────────────────────────────
COL <- list(
  calon      = "#2166AC",   # Deep blue (CALON model)
  safeheart  = "#B2182B",   # Deep red (SAFEHEART-RE)
  core       = "#4393C3",   # Medium blue (CALON-CORE)
  event      = "#D6604D",   # Coral red (ASCVD events)
  noevent    = "#92C5DE",   # Light blue (no events)
  male       = "#4393C3",   # Blue
  female     = "#D6604D",   # Red
  ldlr       = "#2166AC",   # Blue (LDLR gene)
  apob       = "#B2182B",   # Red (APOB gene)
  pcsk9      = "#35978F",   # Teal (PCSK9 gene)
  grey_dark  = "#636363",
  grey_light = "#BDBDBD",
  tertile1   = "#92C5DE",   # Low risk
  tertile2   = "#FDDBC7",   # Medium risk
  tertile3   = "#D6604D"    # High risk
)

# ── Figure dimensions (Nature: single column 88mm, double column 180mm) ─────
WIDTH_SINGLE <- 88 / 25.4   # inches
WIDTH_DOUBLE <- 180 / 25.4  # inches
DPI_PDF <- 300
DPI_PNG <- 600

# ── Load analysis-ready dataset ──────────────────────────────────────────────
DATA_FILE <- paste0(BASE_DIR, "calon_ukb_analysis_ready.csv")
if (file.exists(DATA_FILE)) {
  df <- read.csv(DATA_FILE, stringsAsFactors = FALSE)
  cat(sprintf("Loaded: %d patients from analysis-ready dataset\n", nrow(df)))
} else {
  cat("!! Analysis-ready dataset not found. Generating figures with simulated data.\n")
  cat("   Run 01_CALON_build_cohort.R first to create the real dataset.\n")
  cat("   Figures will be generated with placeholder structure for template review.\n\n")
  # Generate minimal structure for template purposes
  df <- NULL
}

# ── Load validation results ──────────────────────────────────────────────────
RESULTS_FILE <- paste0(TAB_DIR, "calon_validation_results.csv")
if (file.exists(RESULTS_FILE)) {
  results <- read.csv(RESULTS_FILE, stringsAsFactors = FALSE)
  cat(sprintf("Loaded validation results: %d rows\n\n", nrow(results)))
} else {
  results <- NULL
  cat("!! Validation results not found. Run 02_CALON_validate.R first.\n\n")
}

# ── Helper: save figure in both formats ──────────────────────────────────────
save_figure <- function(plot, name, width = WIDTH_DOUBLE, height = 5) {
  ggsave(paste0(FIG_DIR, name, ".pdf"), plot, width = width, height = height,
         dpi = DPI_PDF, device = cairo_pdf)
  ggsave(paste0(FIG_DIR, name, ".png"), plot, width = width, height = height,
         dpi = DPI_PNG, bg = "white")
  cat(sprintf("  Saved: %s.pdf + .png\n", name))
}

# #############################################################################
#                                                                             #
#  FIGURE 1: DUAL-PANEL ROC CURVES                                           #
#  Panel A: Development (Wales) ROC                                           #
#  Panel B: Validation (UKB) ROC — CALON vs SAFEHEART-RE                     #
#                                                                             #
# #############################################################################

cat("\n--- Figure 1: Dual-Panel ROC Curves ---\n")

if (!is.null(df)) {

  # Compute ROC curves (requires predicted probabilities from 02 script)
  if (all(c("calon_prob", "sh_prob", "ascvd_combined") %in% names(df))) {
    roc_calon <- roc(df$ascvd_combined, df$calon_prob, quiet = TRUE)
    roc_sh    <- roc(df$ascvd_combined, df$sh_prob, quiet = TRUE)

    auc_calon <- sprintf("%.3f", as.numeric(roc_calon$auc))
    auc_sh    <- sprintf("%.3f", as.numeric(roc_sh$auc))

    # DeLong test
    delong <- roc.test(roc_calon, roc_sh, method = "delong")

    # Build ROC data
    roc_data <- rbind(
      data.frame(
        model = paste0("CALON (AUC = ", auc_calon, ")"),
        fpr = 1 - roc_calon$specificities,
        tpr = roc_calon$sensitivities
      ),
      data.frame(
        model = paste0("SAFEHEART-RE (AUC = ", auc_sh, ")"),
        fpr = 1 - roc_sh$specificities,
        tpr = roc_sh$sensitivities
      )
    )

    # Panel B: UKB Validation ROC
    fig1b <- ggplot(roc_data, aes(x = fpr, y = tpr, color = model)) +
      geom_line(linewidth = 0.6) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                  color = COL$grey_light, linewidth = 0.3) +
      scale_color_manual(values = c(COL$calon, COL$safeheart)) +
      labs(
        title = "B  UK Biobank External Validation",
        x = "1 - Specificity (False Positive Rate)",
        y = "Sensitivity (True Positive Rate)",
        color = NULL
      ) +
      coord_equal() +
      annotate("text", x = 0.65, y = 0.15,
               label = sprintf("DeLong p = %.4f", delong$p.value),
               size = 2.2, fontface = "italic", color = COL$grey_dark) +
      theme_nature() +
      theme(legend.position = c(0.65, 0.25))

    # If CALON-CORE available
    if ("calon_core_prob" %in% names(df)) {
      roc_core <- roc(df$ascvd_combined, df$calon_core_prob, quiet = TRUE)
      auc_core <- sprintf("%.3f", as.numeric(roc_core$auc))
      roc_core_data <- data.frame(
        model = paste0("CALON-CORE (AUC = ", auc_core, ")"),
        fpr = 1 - roc_core$specificities,
        tpr = roc_core$sensitivities
      )
      roc_data <- rbind(roc_data, roc_core_data)

      fig1b <- ggplot(roc_data, aes(x = fpr, y = tpr, color = model)) +
        geom_line(linewidth = 0.6) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                    color = COL$grey_light, linewidth = 0.3) +
        scale_color_manual(values = c(COL$calon, COL$core, COL$safeheart)) +
        labs(
          title = "B  UK Biobank External Validation",
          x = "1 - Specificity",
          y = "Sensitivity",
          color = NULL
        ) +
        coord_equal() +
        theme_nature() +
        theme(legend.position = c(0.65, 0.25))
    }

    save_figure(fig1b, "Figure1_ROC_Validation", width = WIDTH_SINGLE, height = WIDTH_SINGLE)
  } else {
    cat("  Skipped: predicted probabilities not found in dataset\n")
  }
} else {
  cat("  Skipped: dataset not loaded\n")
}

# #############################################################################
#                                                                             #
#  FIGURE 2: CALIBRATION PLOT (Observed vs Expected)                          #
#                                                                             #
# #############################################################################

cat("\n--- Figure 2: Calibration Plot ---\n")

if (!is.null(df) && "calon_prob" %in% names(df)) {

  # Decile calibration
  df_cal <- df[!is.na(df$calon_prob) & !is.na(df$ascvd_combined), ]
  df_cal$decile <- ntile(df_cal$calon_prob, 10)

  cal_data <- df_cal %>%
    group_by(decile) %>%
    summarise(
      n = n(),
      observed = mean(ascvd_combined),
      expected = mean(calon_prob),
      obs_se = sqrt(observed * (1 - observed) / n),
      .groups = "drop"
    )

  fig2 <- ggplot(cal_data, aes(x = expected, y = observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = COL$grey_light, linewidth = 0.3) +
    geom_point(aes(size = n), color = COL$calon, alpha = 0.8) +
    geom_errorbar(aes(ymin = pmax(0, observed - 1.96 * obs_se),
                      ymax = pmin(1, observed + 1.96 * obs_se)),
                  width = 0, color = COL$calon, linewidth = 0.3) +
    geom_smooth(method = "loess", se = FALSE, color = COL$calon,
                linewidth = 0.5, span = 1.2) +
    scale_size_continuous(range = c(1.5, 5), guide = "none") +
    labs(
      title = "Calibration: Observed vs Expected ASCVD Risk",
      subtitle = "CALON model in UK Biobank FH cohort (by predicted risk decile)",
      x = "Mean Predicted Probability (CALON)",
      y = "Observed Proportion (ASCVD events)"
    ) +
    coord_equal(xlim = c(0, max(cal_data$expected, cal_data$observed) * 1.1),
                ylim = c(0, max(cal_data$expected, cal_data$observed) * 1.1)) +
    theme_nature()

  save_figure(fig2, "Figure2_Calibration", width = WIDTH_SINGLE, height = WIDTH_SINGLE)
} else {
  cat("  Skipped\n")
}

# #############################################################################
#                                                                             #
#  FIGURE 3: GENE-SPECIFIC FOREST PLOT                                        #
#  AUC by gene (LDLR, APOB, PCSK9) with 95% CI                              #
#                                                                             #
# #############################################################################

cat("\n--- Figure 3: Gene-Specific Forest Plot ---\n")

if (!is.null(df) && all(c("calon_prob", "gene", "ascvd_combined") %in% names(df))) {

  gene_results <- list()
  genes <- c("LDLR", "APOB", "PCSK9", "Overall")

  for (g in genes) {
    if (g == "Overall") {
      sub <- df
    } else {
      sub <- df[df$gene == g, ]
    }
    sub <- sub[!is.na(sub$calon_prob) & !is.na(sub$ascvd_combined), ]

    if (nrow(sub) >= 20 && sum(sub$ascvd_combined) >= 3) {
      r <- roc(sub$ascvd_combined, sub$calon_prob, quiet = TRUE)
      ci <- ci.auc(r)
      gene_results[[g]] <- data.frame(
        gene = g, n = nrow(sub), events = sum(sub$ascvd_combined),
        auc = as.numeric(r$auc), lower = ci[1], upper = ci[3],
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(gene_results) > 0) {
    gene_df <- do.call(rbind, gene_results)
    gene_df$gene <- factor(gene_df$gene, levels = rev(c("Overall", "LDLR", "APOB", "PCSK9")))
    gene_df$label <- sprintf("%.3f (%.3f-%.3f)", gene_df$auc, gene_df$lower, gene_df$upper)
    gene_df$is_overall <- gene_df$gene == "Overall"

    fig3 <- ggplot(gene_df, aes(x = auc, y = gene)) +
      geom_vline(xintercept = 0.5, linetype = "dotted", color = COL$grey_light) +
      geom_errorbarh(aes(xmin = lower, xmax = upper),
                     height = 0.15, linewidth = 0.4, color = COL$grey_dark) +
      geom_point(aes(shape = is_overall, color = gene),
                 size = 3) +
      geom_text(aes(label = label), hjust = -0.15, size = 2, color = COL$grey_dark) +
      geom_text(aes(label = sprintf("N=%d (%d events)", n, events)),
                x = 0.45, size = 2, color = COL$grey_dark, hjust = 1) +
      scale_color_manual(values = c("Overall" = "black", "LDLR" = COL$ldlr,
                                    "APOB" = COL$apob, "PCSK9" = COL$pcsk9),
                         guide = "none") +
      scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16), guide = "none") +
      labs(
        title = "CALON Discrimination by Gene",
        subtitle = "AUC (95% CI) for ASCVD prediction in genetically confirmed FH",
        x = "Area Under the ROC Curve (AUC)"
      ) +
      scale_x_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.1)) +
      theme_nature() +
      theme(axis.title.y = element_blank())

    save_figure(fig3, "Figure3_Gene_Forest", width = WIDTH_DOUBLE, height = 3)
  }
} else {
  cat("  Skipped\n")
}

# #############################################################################
#                                                                             #
#  FIGURE 4: KAPLAN-MEIER CURVES BY RISK TERTILE                             #
#                                                                             #
# #############################################################################

cat("\n--- Figure 4: Kaplan-Meier by Risk Tertile ---\n")

if (!is.null(df) && all(c("calon_prob", "time_to_event", "ascvd_incident") %in% names(df))) {

  df_km <- df[!is.na(df$calon_prob) & !is.na(df$time_to_event) & df$time_to_event > 0, ]

  if (nrow(df_km) >= 30) {
    df_km$risk_group <- cut(df_km$calon_prob,
                            breaks = quantile(df_km$calon_prob, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                            labels = c("Low Risk (Tertile 1)", "Medium Risk (Tertile 2)",
                                       "High Risk (Tertile 3)"),
                            include.lowest = TRUE)

    km_fit <- survfit(Surv(time_to_event, ascvd_incident) ~ risk_group, data = df_km)

    fig4 <- ggsurvplot(
      km_fit,
      data = df_km,
      palette = c(COL$tertile1, COL$tertile2, COL$tertile3),
      risk.table = TRUE,
      risk.table.height = 0.25,
      pval = TRUE,
      pval.method = TRUE,
      conf.int = TRUE,
      conf.int.alpha = 0.15,
      xlab = "Follow-up Time (years)",
      ylab = "ASCVD-Free Survival",
      title = "ASCVD-Free Survival by CALON Risk Tertile",
      subtitle = "UK Biobank genetically confirmed FH cohort",
      legend.title = "",
      legend.labs = c("Low Risk", "Medium Risk", "High Risk"),
      font.main = c(9, "bold"),
      font.x = c(7, "bold"),
      font.y = c(7, "bold"),
      font.tickslab = c(6),
      font.legend = c(6),
      ggtheme = theme_nature(base_size = 7)
    )

    pdf(paste0(FIG_DIR, "Figure4_KM_Risk_Tertile.pdf"), width = WIDTH_DOUBLE, height = 5)
    print(fig4)
    dev.off()

    png(paste0(FIG_DIR, "Figure4_KM_Risk_Tertile.png"), width = WIDTH_DOUBLE,
        height = 5, units = "in", res = DPI_PNG)
    print(fig4)
    dev.off()

    cat("  Saved: Figure4_KM_Risk_Tertile.pdf + .png\n")
  }
} else {
  cat("  Skipped\n")
}

# #############################################################################
#                                                                             #
#  FIGURE 5: VARIABLE IMPORTANCE / CONTRIBUTION WATERFALL                     #
#                                                                             #
# #############################################################################

cat("\n--- Figure 5: Variable Contribution Waterfall ---\n")

# From development analysis (these should come from 99_extract script)
# Using manuscript values as defaults
contrib_data <- data.frame(
  variable = c("Age", "RE-LDL-C\n(Phantom Effect)",
               "Lipid-Year\nExposure", "Log(ApoB/LDL)\n(Particle Quality)",
               "Lp(a) >125\nnmol/L", "Inverse\nHDL-C"),
  contribution = c(30.8, 22.8, 19.4, 18.9, 7.2, 1.0),
  group = c("Standard", "Hidden Burden", "Hidden Burden",
            "Hidden Burden", "Standard", "Standard"),
  stringsAsFactors = FALSE
)

contrib_data$variable <- factor(contrib_data$variable,
                                levels = contrib_data$variable[order(-contrib_data$contribution)])

fig5 <- ggplot(contrib_data, aes(x = variable, y = contribution, fill = group)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.1f%%", contribution)),
            vjust = -0.5, size = 2.2, fontface = "bold") +
  scale_fill_manual(
    values = c("Hidden Burden" = COL$calon, "Standard" = COL$grey_light),
    name = NULL
  ) +
  labs(
    title = "CALON Model Variable Contributions",
    subtitle = "Relative importance (|standardised coefficient|)",
    x = NULL,
    y = "Contribution (%)"
  ) +
  annotate("segment", x = 0.5, xend = 4.5, y = 95, yend = 95,
           color = COL$calon, linewidth = 0.8) +
  annotate("text", x = 2.5, y = 98,
           label = "Three Hidden Burdens: 61.1% combined",
           size = 2.2, fontface = "bold.italic", color = COL$calon) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20)) +
  theme_nature() +
  theme(
    axis.text.x = element_text(size = 5.5, lineheight = 0.9),
    legend.position = c(0.85, 0.85)
  )

save_figure(fig5, "Figure5_Variable_Contributions", width = WIDTH_DOUBLE, height = 4)

# #############################################################################
#                                                                             #
#  FIGURE 6: RECLASSIFICATION COMPARISON (NRI Visualization)                  #
#                                                                             #
# #############################################################################

cat("\n--- Figure 6: Reclassification Comparison ---\n")

if (!is.null(df) && all(c("calon_prob", "sh_prob", "ascvd_combined") %in% names(df))) {

  df_nri <- df[!is.na(df$calon_prob) & !is.na(df$sh_prob) & !is.na(df$ascvd_combined), ]

  # Categorise risk
  risk_cat <- function(p) cut(p, breaks = c(0, 0.10, 0.20, 1),
                              labels = c("Low (<10%)", "Moderate (10-20%)", "High (>20%)"),
                              include.lowest = TRUE)

  df_nri$calon_cat <- risk_cat(df_nri$calon_prob)
  df_nri$sh_cat    <- risk_cat(df_nri$sh_prob)
  df_nri$event     <- factor(df_nri$ascvd_combined, levels = c(0, 1),
                             labels = c("No ASCVD", "ASCVD Event"))

  # Build reclassification table
  reclass_tab <- table(SAFEHEART = df_nri$sh_cat, CALON = df_nri$calon_cat)

  # Reclassification heatmap
  reclass_long <- as.data.frame(reclass_tab)
  names(reclass_long) <- c("SAFEHEART", "CALON", "Count")

  fig6 <- ggplot(reclass_long, aes(x = CALON, y = SAFEHEART, fill = Count)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = Count), size = 3, fontface = "bold") +
    scale_fill_gradient(low = "white", high = COL$calon, name = "N patients") +
    labs(
      title = "Risk Reclassification: CALON vs SAFEHEART-RE",
      subtitle = "UK Biobank FH cohort",
      x = "CALON Risk Category",
      y = "SAFEHEART-RE Risk Category"
    ) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))

  save_figure(fig6, "Figure6_Reclassification", width = WIDTH_SINGLE, height = WIDTH_SINGLE)
} else {
  cat("  Skipped\n")
}

# #############################################################################
#                                                                             #
#  FIGURE 7: DECISION CURVE ANALYSIS                                          #
#                                                                             #
# #############################################################################

cat("\n--- Figure 7: Decision Curve Analysis ---\n")

if (!is.null(df) && all(c("calon_prob", "sh_prob", "ascvd_combined") %in% names(df))) {

  # Compute net benefit manually
  thresholds <- seq(0.01, 0.50, 0.01)

  dca_data <- data.frame()
  for (t in thresholds) {
    # CALON
    tp_c <- sum(df$calon_prob >= t & df$ascvd_combined == 1, na.rm = TRUE)
    fp_c <- sum(df$calon_prob >= t & df$ascvd_combined == 0, na.rm = TRUE)
    n <- sum(!is.na(df$ascvd_combined))
    nb_c <- tp_c / n - fp_c / n * (t / (1 - t))

    # SAFEHEART
    tp_s <- sum(df$sh_prob >= t & df$ascvd_combined == 1, na.rm = TRUE)
    fp_s <- sum(df$sh_prob >= t & df$ascvd_combined == 0, na.rm = TRUE)
    nb_s <- tp_s / n - fp_s / n * (t / (1 - t))

    # Treat All
    prevalence <- mean(df$ascvd_combined, na.rm = TRUE)
    nb_all <- prevalence - (1 - prevalence) * (t / (1 - t))

    dca_data <- rbind(dca_data, data.frame(
      threshold = t,
      CALON = nb_c,
      SAFEHEART_RE = nb_s,
      Treat_All = nb_all,
      Treat_None = 0
    ))
  }

  dca_long <- pivot_longer(dca_data, -threshold, names_to = "Strategy", values_to = "Net_Benefit")

  fig7 <- ggplot(dca_long, aes(x = threshold, y = Net_Benefit, color = Strategy)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = c(
      "CALON" = COL$calon,
      "SAFEHEART_RE" = COL$safeheart,
      "Treat_All" = COL$grey_dark,
      "Treat_None" = "black"
    )) +
    labs(
      title = "Decision Curve Analysis",
      subtitle = "Net benefit at varying threshold probabilities",
      x = "Threshold Probability",
      y = "Net Benefit",
      color = NULL
    ) +
    scale_x_continuous(labels = percent_format(), limits = c(0, 0.50)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.2) +
    theme_nature()

  save_figure(fig7, "Figure7_DCA", width = WIDTH_DOUBLE, height = 4)
} else {
  cat("  Skipped\n")
}

# #############################################################################
#                                                                             #
#  FIGURE 8: CLINICAL DECISION PATHWAY SCHEMATIC                              #
#                                                                             #
# #############################################################################

cat("\n--- Figure 8: Clinical Decision Pathway ---\n")

# This figure is a flowchart — built with geom_rect + geom_segment

boxes <- data.frame(
  label = c(
    "Patient with\nGenetically Confirmed FH",
    "Collect:\nAge, LDL-C on treatment,\nStatin type, HDL-C,\nApoB, Lp(a)",
    "Calculate CALON Score\n(6-variable model)",
    "LOW RISK\n(< 10%)",
    "INTERMEDIATE\n(10-20%)",
    "HIGH RISK\n(> 20%)",
    "Standard Lipid\nManagement",
    "ApoB/LDL-C\nRatio Assessment",
    "Intensive Treatment\n+ Specialist Review",
    "Ratio < 0.31:\nStandard Care",
    "Ratio >= 0.31:\nUpgrade to High"
  ),
  x = c(5, 5, 5,    2, 5, 8,    2, 5, 8,    4, 6),
  y = c(10, 8, 6,    4, 4, 4,   2, 2, 2,    0, 0),
  w = c(3, 3.5, 3,   2.5, 2.5, 2.5,  2.5, 2.5, 2.5,  2, 2),
  h = c(1, 1.2, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.8, 0.8),
  fill = c("grey90", "grey90", COL$calon,
           COL$tertile1, COL$tertile2, COL$tertile3,
           COL$tertile1, COL$tertile2, COL$tertile3,
           COL$tertile1, COL$tertile3),
  text_col = c("black", "black", "white",
               "black", "black", "white",
               "black", "black", "white",
               "black", "white"),
  stringsAsFactors = FALSE
)

arrows <- data.frame(
  x = c(5, 5,    5, 5, 5,     2, 5, 8,    5, 5),
  y = c(9.5, 7.4,  5.55, 5.55, 5.55,  3.55, 3.55, 3.55,  1.55, 1.55),
  xend = c(5, 5,    2, 5, 8,    2, 5, 8,   4, 6),
  yend = c(8.6, 6.45, 4.45, 4.45, 4.45,  2.45, 2.45, 2.45,  0.4, 0.4)
)

fig8 <- ggplot() +
  geom_rect(data = boxes,
            aes(xmin = x - w/2, xmax = x + w/2,
                ymin = y - h/2, ymax = y + h/2, fill = fill),
            color = "grey50", linewidth = 0.3) +
  geom_text(data = boxes,
            aes(x = x, y = y, label = label, color = text_col),
            size = 1.8, lineheight = 0.85) +
  geom_segment(data = arrows,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
               linewidth = 0.3, color = "grey40") +
  scale_fill_identity() +
  scale_color_identity() +
  coord_cartesian(xlim = c(0, 10), ylim = c(-0.8, 11)) +
  labs(title = "Clinical Decision Pathway for ASCVD Risk Stratification in FH",
       subtitle = "Using CALON score with ApoB/LDL-C ratio augmentation for intermediate risk") +
  theme_void() +
  theme(
    plot.title = element_text(family = "Arial", face = "bold", size = 8, hjust = 0.5),
    plot.subtitle = element_text(family = "Arial", size = 6, hjust = 0.5, color = "grey40"),
    plot.margin = margin(5, 5, 5, 5)
  )

save_figure(fig8, "Figure8_Clinical_Pathway", width = WIDTH_DOUBLE, height = 5.5)

# #############################################################################
#                                                                             #
#  TABLE 1: BASELINE CHARACTERISTICS (UKB FH Cohort)                          #
#                                                                             #
# #############################################################################

cat("\n--- Table 1: Baseline Characteristics ---\n")

if (!is.null(df)) {

  # Define variables for Table 1
  tab1_vars <- c("age", "sex", "re_ldl", "lipid_years",
                 "log_apob_ldl", "lpa_binary", "inv_hdl",
                 "hdl", "apob", "triglycerides", "total_cholesterol",
                 "hypertension", "diabetes", "current_smoking", "bmi",
                 "on_statin", "ascvd_combined")

  # Keep only variables that exist
  tab1_vars_exist <- tab1_vars[tab1_vars %in% names(df)]

  if (length(tab1_vars_exist) >= 3) {
    # Identify categorical variables
    cat_vars <- intersect(tab1_vars_exist,
                          c("sex", "lpa_binary", "hypertension", "diabetes",
                            "current_smoking", "on_statin", "ascvd_combined"))

    # Stratify by gene if available
    if ("gene" %in% names(df)) {
      tab1 <- CreateTableOne(vars = tab1_vars_exist,
                             strata = "gene",
                             data = df,
                             factorVars = cat_vars,
                             addOverall = TRUE)
    } else {
      tab1 <- CreateTableOne(vars = tab1_vars_exist,
                             data = df,
                             factorVars = cat_vars)
    }

    # Save
    tab1_print <- print(tab1, printToggle = FALSE, showAllLevels = TRUE,
                        formatOptions = list(big.mark = ","))
    write.csv(tab1_print, paste0(TAB_DIR, "Table1_Baseline_Characteristics.csv"))
    cat("  Saved: Table1_Baseline_Characteristics.csv\n")
  }
} else {
  cat("  Skipped\n")
}

# #############################################################################
#                                                                             #
#  TABLE 2: MODEL DISCRIMINATION COMPARISON                                   #
#                                                                             #
# #############################################################################

cat("\n--- Table 2: Model Discrimination ---\n")

if (!is.null(results)) {
  write.csv(results, paste0(TAB_DIR, "Table2_Discrimination_Comparison.csv"), row.names = FALSE)
  cat("  Saved: Table2_Discrimination_Comparison.csv\n")
}

# #############################################################################
#                                                                             #
#  TABLE 3: CALIBRATION METRICS                                               #
#                                                                             #
# #############################################################################

cat("\n--- Table 3: Calibration Metrics ---\n")

if (!is.null(df) && "calon_prob" %in% names(df)) {
  df_cal2 <- df[!is.na(df$calon_prob) & !is.na(df$ascvd_combined), ]

  # Overall calibration
  oe <- sum(df_cal2$ascvd_combined) / sum(df_cal2$calon_prob)
  brier <- mean((df_cal2$calon_prob - df_cal2$ascvd_combined)^2)
  brier_null <- mean((mean(df_cal2$ascvd_combined) - df_cal2$ascvd_combined)^2)
  brier_scaled <- 1 - brier / brier_null

  # Calibration slope/intercept
  lp <- log(df_cal2$calon_prob / (1 - df_cal2$calon_prob + 1e-10))
  cal_mod <- glm(ascvd_combined ~ lp, data = df_cal2, family = binomial)
  cal_slope <- coef(cal_mod)[2]
  cal_int <- coef(cal_mod)[1]

  cal_table <- data.frame(
    Metric = c("O/E Ratio", "Calibration Slope", "Calibration Intercept",
               "Brier Score", "Scaled Brier Score"),
    CALON = c(sprintf("%.3f", oe), sprintf("%.3f", cal_slope),
              sprintf("%.3f", cal_int), sprintf("%.4f", brier),
              sprintf("%.4f", brier_scaled)),
    stringsAsFactors = FALSE
  )

  write.csv(cal_table, paste0(TAB_DIR, "Table3_Calibration_Metrics.csv"), row.names = FALSE)
  cat("  Saved: Table3_Calibration_Metrics.csv\n")
}

# #############################################################################
#                                                                             #
#  TABLE 4: GENE-STRATIFIED PERFORMANCE                                       #
#                                                                             #
# #############################################################################

cat("\n--- Table 4: Gene-Stratified Performance ---\n")

if (!is.null(df) && all(c("calon_prob", "gene", "ascvd_combined") %in% names(df))) {

  gene_table <- df %>%
    filter(!is.na(gene) & !is.na(calon_prob) & !is.na(ascvd_combined)) %>%
    group_by(gene) %>%
    summarise(
      N = n(),
      Events = sum(ascvd_combined),
      Event_Rate = sprintf("%.1f%%", 100 * mean(ascvd_combined)),
      .groups = "drop"
    )

  # Add AUCs
  for (i in 1:nrow(gene_table)) {
    g <- gene_table$gene[i]
    sub <- df[df$gene == g & !is.na(df$calon_prob) & !is.na(df$ascvd_combined), ]
    if (sum(sub$ascvd_combined) >= 3 && sum(sub$ascvd_combined) < nrow(sub)) {
      r <- roc(sub$ascvd_combined, sub$calon_prob, quiet = TRUE)
      ci <- ci.auc(r)
      gene_table$AUC[i] <- sprintf("%.3f (%.3f-%.3f)", as.numeric(r$auc), ci[1], ci[3])
    } else {
      gene_table$AUC[i] <- "N/A (insufficient events)"
    }
  }

  write.csv(gene_table, paste0(TAB_DIR, "Table4_Gene_Stratified.csv"), row.names = FALSE)
  cat("  Saved: Table4_Gene_Stratified.csv\n")
}

# #############################################################################
#                                                                             #
#  TABLE 5: SENSITIVITY ANALYSES SUMMARY                                      #
#                                                                             #
# #############################################################################

cat("\n--- Table 5: Sensitivity Analyses ---\n")

sensitivity_file <- paste0(TAB_DIR, "calon_sensitivity_results.csv")
if (file.exists(sensitivity_file)) {
  sens <- read.csv(sensitivity_file, stringsAsFactors = FALSE)
  write.csv(sens, paste0(TAB_DIR, "Table5_Sensitivity_Analyses.csv"), row.names = FALSE)
  cat("  Saved: Table5_Sensitivity_Analyses.csv\n")
} else {
  cat("  Skipped: sensitivity results not found\n")
}

# #############################################################################
#                                                                             #
#  TABLE 6: SUBGROUP ANALYSES                                                 #
#                                                                             #
# #############################################################################

cat("\n--- Table 6: Subgroup Analyses ---\n")

if (!is.null(df) && all(c("calon_prob", "ascvd_combined") %in% names(df))) {

  subgroups <- list()

  # By sex
  for (s in unique(df$sex[!is.na(df$sex)])) {
    sub <- df[df$sex == s & !is.na(df$calon_prob) & !is.na(df$ascvd_combined), ]
    if (sum(sub$ascvd_combined) >= 3) {
      r <- roc(sub$ascvd_combined, sub$calon_prob, quiet = TRUE)
      ci <- ci.auc(r)
      subgroups[[length(subgroups) + 1]] <- data.frame(
        Subgroup = paste("Sex:", s), N = nrow(sub), Events = sum(sub$ascvd_combined),
        AUC = sprintf("%.3f (%.3f-%.3f)", as.numeric(r$auc), ci[1], ci[3])
      )
    }
  }

  # By age (median split)
  if ("age" %in% names(df)) {
    age_med <- median(df$age, na.rm = TRUE)
    for (grp in c("Young", "Old")) {
      sub <- if (grp == "Young") {
        df[df$age < age_med & !is.na(df$calon_prob) & !is.na(df$ascvd_combined), ]
      } else {
        df[df$age >= age_med & !is.na(df$calon_prob) & !is.na(df$ascvd_combined), ]
      }
      if (sum(sub$ascvd_combined) >= 3) {
        r <- roc(sub$ascvd_combined, sub$calon_prob, quiet = TRUE)
        ci <- ci.auc(r)
        subgroups[[length(subgroups) + 1]] <- data.frame(
          Subgroup = paste0("Age ", grp, " (", ifelse(grp == "Young", "<", ">="),
                            sprintf("%.0f)", age_med)),
          N = nrow(sub), Events = sum(sub$ascvd_combined),
          AUC = sprintf("%.3f (%.3f-%.3f)", as.numeric(r$auc), ci[1], ci[3])
        )
      }
    }
  }

  # By statin status
  if ("on_statin" %in% names(df)) {
    for (st in c(0, 1)) {
      sub <- df[df$on_statin == st & !is.na(df$calon_prob) & !is.na(df$ascvd_combined), ]
      if (sum(sub$ascvd_combined) >= 3) {
        r <- roc(sub$ascvd_combined, sub$calon_prob, quiet = TRUE)
        ci <- ci.auc(r)
        subgroups[[length(subgroups) + 1]] <- data.frame(
          Subgroup = ifelse(st == 1, "On Statin", "Not on Statin"),
          N = nrow(sub), Events = sum(sub$ascvd_combined),
          AUC = sprintf("%.3f (%.3f-%.3f)", as.numeric(r$auc), ci[1], ci[3])
        )
      }
    }
  }

  if (length(subgroups) > 0) {
    subgroup_df <- do.call(rbind, subgroups)
    write.csv(subgroup_df, paste0(TAB_DIR, "Table6_Subgroup_Analyses.csv"), row.names = FALSE)
    cat("  Saved: Table6_Subgroup_Analyses.csv\n")
  }
}

# #############################################################################
#                                                                             #
#  SUPPLEMENTARY: CONSORT FLOW DIAGRAM DATA                                   #
#                                                                             #
# #############################################################################

cat("\n--- Supplementary: CONSORT Flow Data ---\n")

consort_file <- paste0(TAB_DIR, "consort_flowchart.csv")
if (file.exists(consort_file)) {
  consort <- read.csv(consort_file, stringsAsFactors = FALSE)
  cat("  CONSORT data available (", nrow(consort), " rows)\n")
}

# =============================================================================
# COMPLETION
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  FIGURE AND TABLE GENERATION COMPLETE                           \n")
cat("================================================================\n\n")

cat("  Figures saved to: ", FIG_DIR, "\n")
cat("  Tables saved to:  ", TAB_DIR, "\n\n")

cat("  Figure 1: ROC Curves (CALON vs SAFEHEART-RE)\n")
cat("  Figure 2: Calibration Plot (observed vs expected)\n")
cat("  Figure 3: Gene-Specific Forest Plot (AUC by gene)\n")
cat("  Figure 4: Kaplan-Meier by Risk Tertile\n")
cat("  Figure 5: Variable Contribution Waterfall\n")
cat("  Figure 6: Reclassification Comparison\n")
cat("  Figure 7: Decision Curve Analysis\n")
cat("  Figure 8: Clinical Decision Pathway Schematic\n")
cat("\n")
cat("  Table 1: Baseline Characteristics (by gene)\n")
cat("  Table 2: Model Discrimination Comparison\n")
cat("  Table 3: Calibration Metrics\n")
cat("  Table 4: Gene-Stratified Performance\n")
cat("  Table 5: Sensitivity Analyses\n")
cat("  Table 6: Subgroup Analyses\n")
cat("\n")
cat("================================================================\n")

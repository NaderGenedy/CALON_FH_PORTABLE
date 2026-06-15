#!/usr/bin/env Rscript
# =============================================================================
# AF3_R_validation.R
# Comprehensive Statistical Validation of the Structural Severity Score (SSS)
# Nature-calibre formal analysis with 10 sections
# =============================================================================

cat("=============================================================================\n")
cat("SSS FORMAL STATISTICAL VALIDATION\n")
cat("=============================================================================\n\n")

# --- Package Management -------------------------------------------------------
required <- c("dplyr", "ggplot2", "pROC", "boot", "patchwork", "scales",
              "ggrepel", "cowplot", "grid")
for (pkg in required) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# --- Nature Theme -------------------------------------------------------------
theme_nature <- function(base_size = 7) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "Arial", color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 0.5, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_blank(), panel.background = element_blank(),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 0.5),
      legend.position = "bottom", legend.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
    )
}

FIG_DIR <- "C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis/figures/"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
WIDTH_DOUBLE <- 180 / 25.4

# --- Load Data ----------------------------------------------------------------
cat("Loading data...\n")
df <- read.csv("C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis/comprehensive_sss_analysis.csv",
               stringsAsFactors = FALSE)
ukb <- read.csv("C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis/ukb_external_validation_sss.csv",
                stringsAsFactors = FALSE)
nobel <- read.csv("C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis/nobel_multimodal_integration.csv",
                  stringsAsFactors = FALSE)

cat(sprintf("Development cohort: %d rows, %d columns\n", nrow(df), ncol(df)))
cat(sprintf("UKB external validation: %d rows, %d columns\n", nrow(ukb), ncol(ukb)))
cat(sprintf("Nobel integration: %d rows, %d columns\n", nrow(nobel), ncol(nobel)))

# Ensure binary variables
df$ascvd <- as.integer(df$ascvd)
df$sex <- as.integer(df$sex)
df$on_statin <- as.integer(df$on_statin)
df$xanthomata <- as.integer(df$xanthomata)
df$smoking <- as.integer(df$smoking)

ukb$ascvd_combined <- as.integer(ukb$ascvd_combined)
ukb$sex <- as.integer(ukb$sex)
ukb$on_statin <- as.integer(ukb$on_statin)

cat(sprintf("\nDevelopment ASCVD prevalence: %d/%d (%.1f%%)\n",
            sum(df$ascvd, na.rm = TRUE), nrow(df),
            100 * mean(df$ascvd, na.rm = TRUE)))
cat(sprintf("UKB ASCVD prevalence: %d/%d (%.1f%%)\n",
            sum(ukb$ascvd_combined, na.rm = TRUE), nrow(ukb),
            100 * mean(ukb$ascvd_combined, na.rm = TRUE)))

# Master summary collector
summary_rows <- list()

# =============================================================================
# SECTION 1: Logistic Regression with Full Reporting
# =============================================================================
tryCatch({
  cat("\n=============================================================================\n")
  cat("SECTION 1: LOGISTIC REGRESSION - FULL MODEL COMPARISONS\n")
  cat("=============================================================================\n\n")

  models <- list(
    list(name = "Age+Sex",             formula = ascvd ~ age + sex),
    list(name = "Age+Sex+SSS",         formula = ascvd ~ age + sex + sss),
    list(name = "Age+Sex+LDL",         formula = ascvd ~ age + sex + ldl1),
    list(name = "Age+Sex+LDL+SSS",     formula = ascvd ~ age + sex + ldl1 + sss),
    list(name = "Age+Sex+Statin",      formula = ascvd ~ age + sex + on_statin),
    list(name = "Age+Sex+Statin+SSS",  formula = ascvd ~ age + sex + on_statin + sss),
    list(name = "Age+Sex+LDL+Statin",  formula = ascvd ~ age + sex + ldl1 + on_statin),
    list(name = "Age+Sex+LDL+Statin+SSS", formula = ascvd ~ age + sex + ldl1 + on_statin + sss),
    list(name = "Full",                formula = ascvd ~ age + sex + ldl1 + on_statin + smoking + xanthomata),
    list(name = "Full+SSS",            formula = ascvd ~ age + sex + ldl1 + on_statin + smoking + xanthomata + sss)
  )

  model_results <- data.frame(
    Model = character(), N = integer(), AIC = numeric(), BIC = numeric(),
    AUC = numeric(), AUC_lower = numeric(), AUC_upper = numeric(),
    SSS_OR = numeric(), SSS_OR_lower = numeric(), SSS_OR_upper = numeric(),
    SSS_P = numeric(),
    stringsAsFactors = FALSE
  )

  fitted_models <- list()
  roc_objects <- list()

  for (m in models) {
    cat(sprintf("\n--- Model: %s ---\n", m$name))
    vars_needed <- all.vars(m$formula)
    sub <- df[complete.cases(df[, vars_needed, drop = FALSE]), ]
    cat(sprintf("  N = %d (complete cases)\n", nrow(sub)))

    fit <- glm(m$formula, data = sub, family = binomial)
    fitted_models[[m$name]] <- fit

    cat("  Coefficients:\n")
    s <- summary(fit)
    coef_table <- s$coefficients
    ci_raw <- tryCatch(confint(fit), warning = function(w) confint.default(fit))

    for (v in rownames(coef_table)) {
      or_val <- exp(coef_table[v, "Estimate"])
      or_lo  <- exp(ci_raw[v, 1])
      or_hi  <- exp(ci_raw[v, 2])
      pv     <- coef_table[v, "Pr(>|z|)"]
      cat(sprintf("    %-15s  OR = %6.3f (95%% CI: %6.3f - %6.3f)  P = %.2e\n",
                  v, or_val, or_lo, or_hi, pv))
    }

    pred <- predict(fit, type = "response")
    outcome_vec <- sub$ascvd
    roc_obj <- roc(outcome_vec, pred, quiet = TRUE)
    roc_objects[[m$name]] <- roc_obj
    auc_ci <- ci.auc(roc_obj, method = "delong")

    cat(sprintf("  AUC = %.4f (95%% CI: %.4f - %.4f)\n", auc(roc_obj), auc_ci[1], auc_ci[3]))
    cat(sprintf("  AIC = %.1f  BIC = %.1f\n", AIC(fit), BIC(fit)))

    sss_or <- NA; sss_lo <- NA; sss_hi <- NA; sss_p <- NA
    if ("sss" %in% rownames(coef_table)) {
      sss_or <- exp(coef_table["sss", "Estimate"])
      sss_lo <- exp(ci_raw["sss", 1])
      sss_hi <- exp(ci_raw["sss", 2])
      sss_p  <- coef_table["sss", "Pr(>|z|)"]
    }

    model_results <- rbind(model_results, data.frame(
      Model = m$name, N = nrow(sub), AIC = AIC(fit), BIC = BIC(fit),
      AUC = as.numeric(auc(roc_obj)),
      AUC_lower = as.numeric(auc_ci[1]), AUC_upper = as.numeric(auc_ci[3]),
      SSS_OR = sss_or, SSS_OR_lower = sss_lo, SSS_OR_upper = sss_hi,
      SSS_P = sss_p, stringsAsFactors = FALSE
    ))
  }

  cat("\n\n=== MODEL COMPARISON SUMMARY TABLE ===\n")
  print(model_results, row.names = FALSE, digits = 4)
  summary_rows[["Section1"]] <- model_results

}, error = function(e) cat(sprintf("\n[ERROR] Section 1 failed: %s\n", e$message)))


# =============================================================================
# SECTION 2: ROC Curves with DeLong Tests
# =============================================================================
tryCatch({
  cat("\n=============================================================================\n")
  cat("SECTION 2: ROC CURVES WITH DELONG TESTS\n")
  cat("=============================================================================\n\n")

  # Key comparisons
  comparisons <- list(
    list(base = "Age+Sex",         aug = "Age+Sex+SSS"),
    list(base = "Age+Sex+LDL",     aug = "Age+Sex+LDL+SSS"),
    list(base = "Age+Sex+LDL+Statin", aug = "Age+Sex+LDL+Statin+SSS"),
    list(base = "Full",            aug = "Full+SSS")
  )

  for (cmp in comparisons) {
    if (!is.null(roc_objects[[cmp$base]]) && !is.null(roc_objects[[cmp$aug]])) {
      # Refit on same complete cases for valid DeLong
      vars_base <- all.vars(models[[which(sapply(models, function(x) x$name) == cmp$base)]]$formula)
      vars_aug  <- all.vars(models[[which(sapply(models, function(x) x$name) == cmp$aug)]]$formula)
      all_vars  <- unique(c(vars_base, vars_aug))
      sub_common <- df[complete.cases(df[, all_vars, drop = FALSE]), ]

      fit_b <- glm(models[[which(sapply(models, function(x) x$name) == cmp$base)]]$formula,
                    data = sub_common, family = binomial)
      fit_a <- glm(models[[which(sapply(models, function(x) x$name) == cmp$aug)]]$formula,
                    data = sub_common, family = binomial)

      roc_b <- roc(sub_common$ascvd, predict(fit_b, type = "response"), quiet = TRUE)
      roc_a <- roc(sub_common$ascvd, predict(fit_a, type = "response"), quiet = TRUE)
      dt <- roc.test(roc_b, roc_a, method = "delong")

      cat(sprintf("DeLong: %s (AUC=%.4f) vs %s (AUC=%.4f) -- Z=%.3f, P=%.4e\n",
                  cmp$base, auc(roc_b), cmp$aug, auc(roc_a), dt$statistic, dt$p.value))
    }
  }

  # Plot ROC curves for selected models
  plot_names <- c("Age+Sex", "Age+Sex+SSS", "Age+Sex+LDL+SSS",
                  "Age+Sex+LDL+Statin+SSS", "Full+SSS")
  colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#D55E00")

  roc_plot_data <- data.frame()
  for (i in seq_along(plot_names)) {
    nm <- plot_names[i]
    if (!is.null(roc_objects[[nm]])) {
      r <- roc_objects[[nm]]
      ci_val <- ci.auc(r)
      roc_df <- data.frame(
        sensitivity = r$sensitivities, specificity = 1 - r$specificities,
        model = sprintf("%s (AUC=%.3f)", nm, auc(r)),
        stringsAsFactors = FALSE
      )
      roc_plot_data <- rbind(roc_plot_data, roc_df)
    }
  }

  p_roc <- ggplot(roc_plot_data, aes(x = specificity, y = sensitivity, color = model)) +
    geom_line(linewidth = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    scale_color_manual(values = setNames(colors[seq_along(unique(roc_plot_data$model))],
                                         unique(roc_plot_data$model))) +
    labs(x = "1 - Specificity", y = "Sensitivity", title = "ROC Curves - SSS Model Comparisons",
         color = "Model") +
    theme_nature(base_size = 7) +
    theme(legend.position = "right", legend.text = element_text(size = 5))

  ggsave(paste0(FIG_DIR, "Figure_ROC_validation.png"), p_roc,
         width = WIDTH_DOUBLE, height = WIDTH_DOUBLE * 0.65, dpi = 600)
  cat("\nSaved: Figure_ROC_validation.png\n")

}, error = function(e) cat(sprintf("\n[ERROR] Section 2 failed: %s\n", e$message)))


# =============================================================================
# SECTION 3: Calibration Plot
# =============================================================================
tryCatch({
  cat("\n=============================================================================\n")
  cat("SECTION 3: CALIBRATION PLOT - Age+Sex+SSS MODEL\n")
  cat("=============================================================================\n\n")

  vars_cal <- c("ascvd", "age", "sex", "sss")
  sub_cal <- df[complete.cases(df[, vars_cal]), ]
  fit_cal <- glm(ascvd ~ age + sex + sss, data = sub_cal, family = binomial)
  sub_cal$pred <- predict(fit_cal, type = "response")

  # Decile calibration
  sub_cal$decile <- ntile(sub_cal$pred, 10)
  cal_df <- sub_cal %>%
    group_by(decile) %>%
    summarise(
      n = n(),
      observed = mean(ascvd),
      predicted = mean(pred),
      se = sqrt(observed * (1 - observed) / n),
      .groups = "drop"
    )

  # Hosmer-Lemeshow test (manual)
  hl_chi2 <- sum((cal_df$n * (cal_df$observed - cal_df$predicted)^2) /
                   (cal_df$predicted * (1 - cal_df$predicted) + 1e-10))
  hl_df <- nrow(cal_df) - 2
  hl_p <- 1 - pchisq(hl_chi2, df = hl_df)
  cat(sprintf("Hosmer-Lemeshow: chi2 = %.3f, df = %d, P = %.4f\n", hl_chi2, hl_df, hl_p))

  print(cal_df, n = 10)

  p_cal <- ggplot(cal_df, aes(x = predicted, y = observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_point(size = 2.5, color = "#D55E00") +
    geom_errorbar(aes(ymin = pmax(observed - 1.96 * se, 0),
                      ymax = pmin(observed + 1.96 * se, 1)),
                  width = 0.01, color = "#D55E00", linewidth = 0.4) +
    geom_smooth(method = "lm", se = FALSE, color = "#0072B2", linewidth = 0.5) +
    labs(x = "Predicted probability", y = "Observed proportion",
         title = "Calibration Plot - Age+Sex+SSS Model",
         subtitle = sprintf("Hosmer-Lemeshow: chi2=%.2f, P=%.3f", hl_chi2, hl_p)) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_nature(base_size = 7)

  ggsave(paste0(FIG_DIR, "Figure_calibration.png"), p_cal,
         width = WIDTH_DOUBLE * 0.55, height = WIDTH_DOUBLE * 0.55, dpi = 600)
  cat("Saved: Figure_calibration.png\n")

}, error = function(e) cat(sprintf("\n[ERROR] Section 3 failed: %s\n", e$message)))


# =============================================================================
# SECTION 4: Decision Curve Analysis
# =============================================================================
tryCatch({
  cat("\n=============================================================================\n")
  cat("SECTION 4: DECISION CURVE ANALYSIS\n")
  cat("=============================================================================\n\n")

  decision_curve <- function(outcome, pred, thresholds = seq(0.01, 0.50, 0.01)) {
    n <- length(outcome)
    nb <- data.frame(threshold = thresholds, model = NA, treat_all = NA, treat_none = 0)
    for (i in seq_along(thresholds)) {
      pt <- thresholds[i]
      pos <- pred >= pt
      tp <- sum(pos & outcome == 1)
      fp <- sum(pos & outcome == 0)
      nb$model[i] <- tp / n - fp / n * (pt / (1 - pt))
      nb$treat_all[i] <- sum(outcome == 1) / n - sum(outcome == 0) / n * (pt / (1 - pt))
    }
    nb
  }

  vars_dca <- c("ascvd", "age", "sex", "sss")
  sub_dca <- df[complete.cases(df[, vars_dca]), ]

  fit_base_dca <- glm(ascvd ~ age + sex, data = sub_dca, family = binomial)
  fit_sss_dca  <- glm(ascvd ~ age + sex + sss, data = sub_dca, family = binomial)

  pred_base <- predict(fit_base_dca, type = "response")
  pred_sss  <- predict(fit_sss_dca, type = "response")

  dc_base <- decision_curve(sub_dca$ascvd, pred_base)
  dc_sss  <- decision_curve(sub_dca$ascvd, pred_sss)

  dca_plot_data <- rbind(
    data.frame(threshold = dc_base$threshold, nb = dc_base$model, Strategy = "Age+Sex"),
    data.frame(threshold = dc_sss$threshold, nb = dc_sss$model, Strategy = "Age+Sex+SSS"),
    data.frame(threshold = dc_base$threshold, nb = dc_base$treat_all, Strategy = "Treat All"),
    data.frame(threshold = dc_base$threshold, nb = dc_base$treat_none, Strategy = "Treat None")
  )

  p_dca <- ggplot(dca_plot_data, aes(x = threshold, y = nb, color = Strategy, linetype = Strategy)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = c("Age+Sex" = "#999999", "Age+Sex+SSS" = "#D55E00",
                                  "Treat All" = "#0072B2", "Treat None" = "black")) +
    scale_linetype_manual(values = c("Age+Sex" = "solid", "Age+Sex+SSS" = "solid",
                                     "Treat All" = "dashed", "Treat None" = "dotted")) +
    labs(x = "Threshold probability", y = "Net benefit",
         title = "Decision Curve Analysis", color = "Strategy", linetype = "Strategy") +
    ylim(-0.05, NA) +
    theme_nature(base_size = 7)

  ggsave(paste0(FIG_DIR, "Figure_DCA.png"), p_dca,
         width = WIDTH_DOUBLE * 0.65, height = WIDTH_DOUBLE * 0.5, dpi = 600)
  cat("Saved: Figure_DCA.png\n")

}, error = function(e) cat(sprintf("\n[ERROR] Section 4 failed: %s\n", e$message)))


# =============================================================================
# SECTION 5: Forest Plot - SSS Effect Across Subgroups
# =============================================================================
tryCatch({
  cat("\n=============================================================================\n")
  cat("SECTION 5: FOREST PLOT - SSS EFFECT ACROSS SUBGROUPS\n")
  cat("=============================================================================\n\n")

  run_subgroup <- function(label, subset_df) {
    sub <- subset_df[complete.cases(subset_df[, c("ascvd", "age", "sex", "sss")]), ]
    if (nrow(sub) < 20 || length(unique(sub$ascvd)) < 2) {
      return(data.frame(Subgroup = label, N = nrow(sub), Events = sum(sub$ascvd),
                        OR = NA, OR_lower = NA, OR_upper = NA, P = NA))
    }
    fit <- glm(ascvd ~ age + sex + sss, data = sub, family = binomial)
    ci <- tryCatch(confint(fit), warning = function(w) confint.default(fit))
    s <- summary(fit)$coefficients
    data.frame(
      Subgroup = label, N = nrow(sub), Events = sum(sub$ascvd),
      OR = exp(s["sss", "Estimate"]),
      OR_lower = exp(ci["sss", 1]),
      OR_upper = exp(ci["sss", 2]),
      P = s["sss", "Pr(>|z|)"]
    )
  }

  forest_data <- rbind(
    run_subgroup("All patients", df),
    run_subgroup("Age >= 50", df[!is.na(df$age) & df$age >= 50, ]),
    run_subgroup("Age < 50", df[!is.na(df$age) & df$age < 50, ]),
    run_subgroup("Male", df[!is.na(df$sex) & df$sex == 1, ]),
    run_subgroup("Female", df[!is.na(df$sex) & df$sex == 0, ]),
    run_subgroup("LDLR only", df[!is.na(df$gene) & df$gene == "LDLR", ]),
    run_subgroup("On statin", df[!is.na(df$on_statin) & df$on_statin == 1, ]),
    run_subgroup("No statin", df[!is.na(df$on_statin) & df$on_statin == 0, ]),
    run_subgroup("With xanthomata", df[!is.na(df$xanthomata) & df$xanthomata == 1, ]),
    run_subgroup("Without xanthomata", df[!is.na(df$xanthomata) & df$xanthomata == 0, ]),
    run_subgroup("Simon Broome definite", df[!is.na(df$simon_broome) & df$simon_broome == 2, ])
  )

  forest_data <- forest_data[!is.na(forest_data$OR), ]
  forest_data$Subgroup <- factor(forest_data$Subgroup, levels = rev(forest_data$Subgroup))

  cat("Forest plot data:\n")
  print(forest_data, row.names = FALSE, digits = 3)

  p_forest <- ggplot(forest_data, aes(x = OR, y = Subgroup)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_point(size = 2, color = "#D55E00") +
    geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), height = 0.2,
                   color = "#D55E00", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f (%.2f-%.2f)", OR, OR_lower, OR_upper)),
              hjust = -0.1, size = 2, nudge_x = 0.05) +
    scale_x_log10() +
    labs(x = "Odds Ratio (log scale)", y = "",
         title = "SSS Effect on ASCVD Across Subgroups",
         subtitle = "Model: ascvd ~ age + sex + sss") +
    theme_nature(base_size = 7) +
    theme(axis.text.y = element_text(size = 6))

  ggsave(paste0(FIG_DIR, "Figure_forest_SSS.png"), p_forest,
         width = WIDTH_DOUBLE, height = WIDTH_DOUBLE * 0.6, dpi = 600)
  cat("Saved: Figure_forest_SSS.png\n")

}, error = function(e) cat(sprintf("\n[ERROR] Section 5 failed: %s\n", e$message)))


# =============================================================================
# SECTION 6: SSS Interaction Terms
# =============================================================================
tryCatch({
  cat("\n=============================================================================\n")
  cat("SECTION 6: SSS INTERACTION TERMS\n")
  cat("=============================================================================\n\n")

  # SSS x Statin
  vars_int <- c("ascvd", "age", "sex", "sss", "on_statin")
  sub_int <- df[complete.cases(df[, vars_int]), ]

  m_main   <- glm(ascvd ~ age + sex + sss + on_statin, data = sub_int, family = binomial)
  m_int_st <- glm(ascvd ~ age + sex + sss + on_statin + sss:on_statin, data = sub_int, family = binomial)

  cat("--- SSS x Statin Interaction ---\n")
  s_int <- summary(m_int_st)$coefficients
  ci_int <- tryCatch(confint(m_int_st), warning = function(w) confint.default(m_int_st))
  for (v in rownames(s_int)) {
    cat(sprintf("  %-25s  OR = %6.3f (95%% CI: %6.3f - %6.3f)  P = %.2e\n",
                v, exp(s_int[v, "Estimate"]), exp(ci_int[v, 1]), exp(ci_int[v, 2]),
                s_int[v, "Pr(>|z|)"]))
  }

  lr_statin <- anova(m_main, m_int_st, test = "Chisq")
  cat(sprintf("\nLR test (interaction term): chi2 = %.3f, P = %.4f\n",
              lr_statin$Deviance[2], lr_statin$`Pr(>Chi)`[2]))

  # SSS x Age
  m_main_age <- glm(ascvd ~ sex + sss + age, data = sub_int, family = binomial)
  m_int_age  <- glm(ascvd ~ sex + sss + age + sss:age, data = sub_int, family = binomial)

  cat("\n--- SSS x Age Interaction ---\n")
  s_age <- summary(m_int_age)$coefficients
  ci_age <- tryCatch(confint(m_int_age), warning = function(w) confint.default(m_int_age))
  for (v in rownames(s_age)) {
    cat(sprintf("  %-25s  OR = %6.3f (95%% CI: %6.3f - %6.3f)  P = %.2e\n",
                v, exp(s_age[v, "Estimate"]), exp(ci_age[v, 1]), exp(ci_age[v, 2]),
                s_age[v, "Pr(>|z|)"]))
  }

  lr_age <- anova(m_main_age, m_int_age, test = "Chisq")
  cat(sprintf("\nLR test (age interaction): chi2 = %.3f, P = %.4f\n",
              lr_age$Deviance[2], lr_age$`Pr(>Chi)`[2]))

}, error = function(e) cat(sprintf("\n[ERROR] Section 6 failed: %s\n", e$message)))


# =============================================================================
# SECTION 7: External Validation in UKB
# =============================================================================
tryCatch({
  cat("\n=============================================================================\n")
  cat("SECTION 7: EXTERNAL VALIDATION IN UK BIOBANK\n")
  cat("=============================================================================\n\n")

  cat(sprintf("UKB cohort: N=%d, Events=%d (%.1f%%)\n",
              nrow(ukb), sum(ukb$ascvd_combined, na.rm = TRUE),
              100 * mean(ukb$ascvd_combined, na.rm = TRUE)))

  ukb_models <- list(
    list(name = "UKB: Age+Sex",          formula = ascvd_combined ~ age + sex),
    list(name = "UKB: Age+Sex+SSS",      formula = ascvd_combined ~ age + sex + sss),
    list(name = "UKB: Age+Sex+LDL",      formula = ascvd_combined ~ age + sex + ldl),
    list(name = "UKB: Age+Sex+LDL+SSS",  formula = ascvd_combined ~ age + sex + ldl + sss),
    list(name = "UKB: Age+Sex+LDL+Statin",    formula = ascvd_combined ~ age + sex + ldl + on_statin),
    list(name = "UKB: Age+Sex+LDL+Statin+SSS", formula = ascvd_combined ~ age + sex + ldl + on_statin + sss)
  )

  ukb_results <- data.frame(
    Model = character(), N = integer(), AUC = numeric(),
    AUC_lower = numeric(), AUC_upper = numeric(),
    SSS_OR = numeric(), SSS_OR_lower = numeric(), SSS_OR_upper = numeric(),
    SSS_P = numeric(), stringsAsFactors = FALSE
  )

  for (m in ukb_models) {
    vars_needed <- all.vars(m$formula)
    sub <- ukb[complete.cases(ukb[, vars_needed, drop = FALSE]), ]
    if (nrow(sub) < 20 || length(unique(sub$ascvd_combined)) < 2) {
      cat(sprintf("  Skipping %s (insufficient data)\n", m$name))
      next
    }

    fit <- glm(m$formula, data = sub, family = binomial)
    s <- summary(fit)$coefficients
    ci <- tryCatch(confint(fit), warning = function(w) confint.default(fit))

    pred <- predict(fit, type = "response")
    roc_obj <- roc(sub$ascvd_combined, pred, quiet = TRUE)
    auc_ci <- ci.auc(roc_obj, method = "delong")

    cat(sprintf("\n--- %s (N=%d) ---\n", m$name, nrow(sub)))
    for (v in rownames(s)) {
      cat(sprintf("  %-15s  OR = %6.3f (95%% CI: %6.3f - %6.3f)  P = %.2e\n",
                  v, exp(s[v, "Estimate"]), exp(ci[v, 1]), exp(ci[v, 2]),
                  s[v, "Pr(>|z|)"]))
    }
    cat(sprintf("  AUC = %.4f (95%% CI: %.4f - %.4f)\n",
                auc(roc_obj), auc_ci[1], auc_ci[3]))

    sss_or <- NA; sss_lo <- NA; sss_hi <- NA; sss_p <- NA
    if ("sss" %in% rownames(s)) {
      sss_or <- exp(s["sss", "Estimate"])
      sss_lo <- exp(ci["sss", 1])
      sss_hi <- exp(ci["sss", 2])
      sss_p  <- s["sss", "Pr(>|z|)"]
    }

    ukb_results <- rbind(ukb_results, data.frame(
      Model = m$name, N = nrow(sub),
      AUC = as.numeric(auc(roc_obj)),
      AUC_lower = as.numeric(auc_ci[1]),
      AUC_upper = as.numeric(auc_ci[3]),
      SSS_OR = sss_or, SSS_OR_lower = sss_lo, SSS_OR_upper = sss_hi,
      SSS_P = sss_p, stringsAsFactors = FALSE
    ))
  }

  cat("\n=== UKB EXTERNAL VALIDATION SUMMARY ===\n")
  print(ukb_results, row.names = FALSE, digits = 4)

  # Calibration in UKB
  vars_ukb_cal <- c("ascvd_combined", "age", "sex", "sss")
  sub_ukb_cal <- ukb[complete.cases(ukb[, vars_ukb_cal]), ]
  fit_ukb_cal <- glm(ascvd_combined ~ age + sex + sss, data = sub_ukb_cal, family = binomial)
  sub_ukb_cal$pred <- predict(fit_ukb_cal, type = "response")
  sub_ukb_cal$decile <- ntile(sub_ukb_cal$pred, 5)  # quintiles for smaller N
  cal_ukb <- sub_ukb_cal %>%
    group_by(decile) %>%
    summarise(n = n(), observed = mean(ascvd_combined), predicted = mean(pred),
              .groups = "drop")
  cat("\nUKB Calibration (quintiles):\n")
  print(cal_ukb, n = 5)

  summary_rows[["Section7"]] <- ukb_results

}, error = function(e) cat(sprintf("\n[ERROR] Section 7 failed: %s\n", e$message)))


# =============================================================================
# SECTION 8: NRI and IDI (Manual Bootstrap)
# =============================================================================
tryCatch({
  cat("\n=============================================================================\n")
  cat("SECTION 8: NET RECLASSIFICATION IMPROVEMENT & IDI\n")
  cat("=============================================================================\n\n")

  vars_nri <- c("ascvd", "age", "sex", "sss")
  sub_nri <- df[complete.cases(df[, vars_nri]), ]

  base_model <- glm(ascvd ~ age + sex, data = sub_nri, family = binomial)
  aug_model  <- glm(ascvd ~ age + sex + sss, data = sub_nri, family = binomial)

  base_pred <- predict(base_model, type = "response")
  aug_pred  <- predict(aug_model, type = "response")
  events    <- sub_nri$ascvd == 1

  # Continuous NRI
  nri_events    <- mean(aug_pred[events] > base_pred[events]) - mean(aug_pred[events] < base_pred[events])
  nri_nonevents <- mean(aug_pred[!events] < base_pred[!events]) - mean(aug_pred[!events] > base_pred[!events])
  nri_total     <- nri_events + nri_nonevents

  # IDI
  idi <- mean(aug_pred[events]) - mean(base_pred[events]) -
    (mean(aug_pred[!events]) - mean(base_pred[!events]))

  cat(sprintf("Continuous NRI (events):     %.4f\n", nri_events))
  cat(sprintf("Continuous NRI (non-events): %.4f\n", nri_nonevents))
  cat(sprintf("Continuous NRI (total):      %.4f\n", nri_total))
  cat(sprintf("IDI:                         %.4f\n", idi))

  # Bootstrap 95% CI
  set.seed(42)
  n_boot <- 2000
  boot_nri <- numeric(n_boot)
  boot_idi <- numeric(n_boot)

  for (b in 1:n_boot) {
    idx <- sample(nrow(sub_nri), replace = TRUE)
    bdata <- sub_nri[idx, ]
    b_base <- glm(ascvd ~ age + sex, data = bdata, family = binomial)
    b_aug  <- glm(ascvd ~ age + sex + sss, data = bdata, family = binomial)
    bp <- predict(b_base, type = "response")
    ap <- predict(b_aug, type = "response")
    ev <- bdata$ascvd == 1

    if (sum(ev) > 0 && sum(!ev) > 0) {
      nri_ev  <- mean(ap[ev] > bp[ev]) - mean(ap[ev] < bp[ev])
      nri_nev <- mean(ap[!ev] < bp[!ev]) - mean(ap[!ev] > bp[!ev])
      boot_nri[b] <- nri_ev + nri_nev

      boot_idi[b] <- mean(ap[ev]) - mean(bp[ev]) - (mean(ap[!ev]) - mean(bp[!ev]))
    } else {
      boot_nri[b] <- NA
      boot_idi[b] <- NA
    }
  }

  nri_ci <- quantile(boot_nri, c(0.025, 0.975), na.rm = TRUE)
  idi_ci <- quantile(boot_idi, c(0.025, 0.975), na.rm = TRUE)
  nri_p <- 2 * min(mean(boot_nri > 0, na.rm = TRUE), mean(boot_nri < 0, na.rm = TRUE))
  idi_p <- 2 * min(mean(boot_idi > 0, na.rm = TRUE), mean(boot_idi < 0, na.rm = TRUE))

  cat(sprintf("\nNRI = %.4f (95%% CI: %.4f - %.4f), P = %.4f\n",
              nri_total, nri_ci[1], nri_ci[2], nri_p))
  cat(sprintf("IDI = %.4f (95%% CI: %.4f - %.4f), P = %.4f\n",
              idi, idi_ci[1], idi_ci[2], idi_p))

  # Category-based NRI (using 10% and 30% thresholds)
  cat("\n--- Category-based NRI (thresholds: 10%, 30%) ---\n")
  classify <- function(p) {
    ifelse(p < 0.10, "low", ifelse(p < 0.30, "intermediate", "high"))
  }
  base_cat <- classify(base_pred)
  aug_cat  <- classify(aug_pred)

  reclass_events <- table(base_cat[events], aug_cat[events])
  reclass_nonevents <- table(base_cat[!events], aug_cat[!events])

  cat("Reclassification among events:\n")
  print(reclass_events)
  cat("Reclassification among non-events:\n")
  print(reclass_nonevents)

}, error = function(e) cat(sprintf("\n[ERROR] Section 8 failed: %s\n", e$message)))


# =============================================================================
# SECTION 9: Treatment Response Analysis
# =============================================================================
tryCatch({
  cat("\n=============================================================================\n")
  cat("SECTION 9: TREATMENT RESPONSE ANALYSIS\n")
  cat("=============================================================================\n\n")

  treated <- df[!is.na(df$on_statin) & df$on_statin == 1 &
                  !is.na(df$ldl_pct_change) & !is.na(df$sss), ]
  cat(sprintf("Treated patients with LDL response: N = %d\n", nrow(treated)))

  if (nrow(treated) >= 10) {
    m_resp <- lm(ldl_pct_change ~ sss + age + sex, data = treated)
    cat("\nLinear regression: ldl_pct_change ~ sss + age + sex\n")
    s_resp <- summary(m_resp)
    print(s_resp)

    ci_resp <- confint(m_resp)
    cat("\n95% CIs for coefficients:\n")
    print(ci_resp)

    # Scatter plot
    p_resp <- ggplot(treated, aes(x = sss, y = ldl_pct_change)) +
      geom_point(alpha = 0.5, size = 1.5, color = "#0072B2") +
      geom_smooth(method = "lm", se = TRUE, color = "#D55E00", linewidth = 0.6, fill = "#D55E0033") +
      labs(x = "Structural Severity Score (SSS)",
           y = "LDL-C % Change on Statin",
           title = "SSS vs Treatment Response",
           subtitle = sprintf("beta=%.2f, P=%.3f, R2=%.3f",
                              coef(m_resp)["sss"],
                              s_resp$coefficients["sss", "Pr(>|t|)"],
                              s_resp$r.squared)) +
      theme_nature(base_size = 7)

    ggsave(paste0(FIG_DIR, "Figure_treatment_response.png"), p_resp,
           width = WIDTH_DOUBLE * 0.55, height = WIDTH_DOUBLE * 0.45, dpi = 600)
    cat("\nSaved: Figure_treatment_response.png\n")

    # Domain-specific ANOVA
    treated_domain <- treated[!is.na(treated$domain) & treated$domain != "" & treated$domain != "Unknown", ]
    if (nrow(treated_domain) >= 10 && length(unique(treated_domain$domain)) >= 2) {
      cat("\n--- Domain-specific Treatment Response (ANOVA) ---\n")
      aov_fit <- aov(ldl_pct_change ~ domain, data = treated_domain)
      cat("\nANOVA table:\n")
      print(summary(aov_fit))

      tukey_res <- TukeyHSD(aov_fit)
      cat("\nTukey HSD post-hoc:\n")
      print(tukey_res)
    }
  }

}, error = function(e) cat(sprintf("\n[ERROR] Section 9 failed: %s\n", e$message)))


# =============================================================================
# SECTION 10: ddG-Clinical Correlations
# =============================================================================
tryCatch({
  cat("\n=============================================================================\n")
  cat("SECTION 10: ddG-CLINICAL CORRELATIONS\n")
  cat("=============================================================================\n\n")

  # ddG vs last_ldl
  sub_ddg_ldl <- df[!is.na(df$ddG) & !is.na(df$last_ldl), ]
  cat(sprintf("ddG vs last_ldl: N = %d\n", nrow(sub_ddg_ldl)))
  if (nrow(sub_ddg_ldl) >= 5) {
    ct_ldl <- cor.test(sub_ddg_ldl$ddG, sub_ddg_ldl$last_ldl, method = "pearson")
    cat(sprintf("  Pearson r = %.4f (95%% CI: %.4f - %.4f), P = %.4e\n",
                ct_ldl$estimate, ct_ldl$conf.int[1], ct_ldl$conf.int[2], ct_ldl$p.value))

    sp_ldl <- cor.test(sub_ddg_ldl$ddG, sub_ddg_ldl$last_ldl, method = "spearman")
    cat(sprintf("  Spearman rho = %.4f, P = %.4e\n", sp_ldl$estimate, sp_ldl$p.value))
  }

  # ddG vs apob
  sub_ddg_apob <- df[!is.na(df$ddG) & !is.na(df$apob), ]
  cat(sprintf("\nddG vs ApoB: N = %d\n", nrow(sub_ddg_apob)))
  if (nrow(sub_ddg_apob) >= 5) {
    ct_apob <- cor.test(sub_ddg_apob$ddG, sub_ddg_apob$apob, method = "pearson")
    cat(sprintf("  Pearson r = %.4f (95%% CI: %.4f - %.4f), P = %.4e\n",
                ct_apob$estimate, ct_apob$conf.int[1], ct_apob$conf.int[2], ct_apob$p.value))

    sp_apob <- cor.test(sub_ddg_apob$ddG, sub_ddg_apob$apob, method = "spearman")
    cat(sprintf("  Spearman rho = %.4f, P = %.4e\n", sp_apob$estimate, sp_apob$p.value))
  }

  # Combined plot
  p_ddg1 <- ggplot(sub_ddg_ldl, aes(x = ddG, y = last_ldl)) +
    geom_point(alpha = 0.5, size = 1.5, color = "#0072B2") +
    geom_smooth(method = "lm", se = TRUE, color = "#D55E00", linewidth = 0.6, fill = "#D55E0033") +
    labs(x = expression(Delta * Delta * "G (kcal/mol)"), y = "Last LDL-C (mmol/L)",
         title = "ddG vs LDL-C",
         subtitle = sprintf("r=%.3f, P=%.2e", ct_ldl$estimate, ct_ldl$p.value)) +
    theme_nature(base_size = 7)

  p_ddg2 <- ggplot(sub_ddg_apob, aes(x = ddG, y = apob)) +
    geom_point(alpha = 0.5, size = 1.5, color = "#009E73") +
    geom_smooth(method = "lm", se = TRUE, color = "#D55E00", linewidth = 0.6, fill = "#D55E0033") +
    labs(x = expression(Delta * Delta * "G (kcal/mol)"), y = "ApoB (g/L)",
         title = "ddG vs ApoB",
         subtitle = sprintf("r=%.3f, P=%.2e", ct_apob$estimate, ct_apob$p.value)) +
    theme_nature(base_size = 7)

  p_ddg_combined <- p_ddg1 + p_ddg2 + plot_layout(ncol = 2)

  ggsave(paste0(FIG_DIR, "Figure_ddG_clinical.png"), p_ddg_combined,
         width = WIDTH_DOUBLE, height = WIDTH_DOUBLE * 0.45, dpi = 600)
  cat("\nSaved: Figure_ddG_clinical.png\n")

  # ddG vs SSS correlation
  sub_ddg_sss <- df[!is.na(df$ddG) & !is.na(df$sss), ]
  if (nrow(sub_ddg_sss) >= 5) {
    ct_sss_ddg <- cor.test(sub_ddg_sss$ddG, sub_ddg_sss$sss, method = "pearson")
    cat(sprintf("\nddG vs SSS: r = %.4f (95%% CI: %.4f - %.4f), P = %.4e, N = %d\n",
                ct_sss_ddg$estimate, ct_sss_ddg$conf.int[1], ct_sss_ddg$conf.int[2],
                ct_sss_ddg$p.value, nrow(sub_ddg_sss)))
  }

}, error = function(e) cat(sprintf("\n[ERROR] Section 10 failed: %s\n", e$message)))


# =============================================================================
# SAVE SUMMARY TABLE
# =============================================================================
tryCatch({
  cat("\n=============================================================================\n")
  cat("SAVING SUMMARY OUTPUT\n")
  cat("=============================================================================\n\n")

  if (!is.null(summary_rows[["Section1"]])) {
    write.csv(summary_rows[["Section1"]],
              paste0(FIG_DIR, "validation_model_comparison.csv"), row.names = FALSE)
    cat("Saved: validation_model_comparison.csv\n")
  }

  if (!is.null(summary_rows[["Section7"]])) {
    write.csv(summary_rows[["Section7"]],
              paste0(FIG_DIR, "ukb_external_validation_results.csv"), row.names = FALSE)
    cat("Saved: ukb_external_validation_results.csv\n")
  }

}, error = function(e) cat(sprintf("\n[ERROR] Summary save failed: %s\n", e$message)))


cat("\n=============================================================================\n")
cat("ALL SECTIONS COMPLETE\n")
cat("=============================================================================\n")

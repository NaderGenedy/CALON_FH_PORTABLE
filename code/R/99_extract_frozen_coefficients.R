################################################################################
#                                                                              #
#  99_EXTRACT_FROZEN_COEFFICIENTS.R                                            #
#  CALON-FH: Extract Exact Frozen Coefficients from Development Data           #
#                                                                              #
#  Purpose: Fit the EXACT 6-variable CALON model on REAL DRAGON2.csv data      #
#           Extract frozen coefficients for TRIPOD Type 4 external validation  #
#           Also fit CALON-CORE and SAFEHEART-RE for comparison                #
#                                                                              #
#  Author: CALON-FH Pipeline                                                   #
#  Date: February 2026                                                         #
#  Data: DRAGON2.csv (Wales FH Registry, N=424)                                #
#                                                                              #
#  CRITICAL: These coefficients are FROZEN and must be applied WITHOUT          #
#  re-estimation in the UKB validation cohort (TRIPOD Type 4)                  #
#                                                                              #
################################################################################

# =============================================================================
# SETUP
# =============================================================================

rm(list = ls())
set.seed(2026)

# Packages
library(dplyr)
library(pROC)
library(boot)

cat("\n")
cat("================================================================\n")
cat("  CALON-FH: FROZEN COEFFICIENT EXTRACTION                       \n")
cat("  From: Wales FH Registry (DRAGON2.csv)                         \n")
cat("  For:  UK Biobank External Validation (TRIPOD Type 4)          \n")
cat("================================================================\n\n")

# =============================================================================
# 1. LOAD DATA
# =============================================================================

DATA_FILE <- "C:/Users/nader/Downloads/DRAGON2.csv"
if (!file.exists(DATA_FILE)) {
  DATA_FILE <- "DRAGON2.csv"
}

df <- read.csv(DATA_FILE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
cat(sprintf("Loaded: %d patients from DRAGON2.csv\n", nrow(df)))

# =============================================================================
# 2. VARIABLE PREPARATION - Matching Manuscript Specification Exactly
# =============================================================================

cat("\n--- Variable Preparation ---\n\n")

# --- Outcome ---
df$ASCVD <- as.numeric(df$ASCVD_combined == 1)
cat(sprintf("  Outcome (ASCVD_combined): %d events (%.1f%%)\n",
            sum(df$ASCVD, na.rm = TRUE),
            100 * mean(df$ASCVD, na.rm = TRUE)))

# --- Age ---
df$age <- as.numeric(df$AGE_AT_TEST)
cat(sprintf("  Age at test: N=%d, mean=%.1f (SD=%.1f)\n",
            sum(!is.na(df$age)), mean(df$age, na.rm = TRUE), sd(df$age, na.rm = TRUE)))

# --- Reverse-Engineered LDL-C ---
# MtachedLDLC in DRAGON2.csv IS the reverse-engineered LDL (matched pre-treatment)
df$re_ldl <- as.numeric(df$MtachedLDLC)
cat(sprintf("  RE-LDL (MtachedLDLC): N=%d, mean=%.2f (SD=%.2f)\n",
            sum(!is.na(df$re_ldl)), mean(df$re_ldl, na.rm = TRUE), sd(df$re_ldl, na.rm = TRUE)))

# --- Age at Treatment Start ---
df$age_tx <- as.numeric(df$Age_Treatment1)
cat(sprintf("  Age at treatment: N=%d, mean=%.1f (SD=%.1f)\n",
            sum(!is.na(df$age_tx)), mean(df$age_tx, na.rm = TRUE), sd(df$age_tx, na.rm = TRUE)))

# --- Lipid-Year Exposure (Cholesterol Pack-Years) ---
# Formula: RE_LDL * (age_at_treatment - 18) / 10
# This captures cumulative untreated LDL-C exposure from age 18 to treatment start
df$lipid_years <- df$re_ldl * pmax(df$age_tx - 18, 0) / 10
cat(sprintf("  Lipid-years: N=%d, mean=%.1f (SD=%.1f)\n",
            sum(!is.na(df$lipid_years)), mean(df$lipid_years, na.rm = TRUE),
            sd(df$lipid_years, na.rm = TRUE)))

# --- Log(ApoB / LDL-C) ---
# Already computed in DRAGON2.csv as LOG_APOB_LDL
# Verify it matches our calculation
df$ApoB_raw <- as.numeric(df$ApoB)
df$LDL_1_raw <- as.numeric(df$LDL_1)
df$log_apob_ldl_existing <- as.numeric(df$LOG_APOB_LDL)

# Compute fresh from raw values for verification
# Note: ApoB is in g/L, LDL is in mmol/L → ratio is g/mmol
df$log_apob_ldl <- log(df$ApoB_raw / df$LDL_1_raw)

# Check correlation between existing and fresh computation
both_valid <- !is.na(df$log_apob_ldl) & !is.na(df$log_apob_ldl_existing)
if (sum(both_valid) > 10) {
  r <- cor(df$log_apob_ldl[both_valid], df$log_apob_ldl_existing[both_valid])
  cat(sprintf("  Log(ApoB/LDL): N=%d, mean=%.3f (SD=%.3f) [correlation with existing=%.4f]\n",
              sum(!is.na(df$log_apob_ldl)), mean(df$log_apob_ldl, na.rm = TRUE),
              sd(df$log_apob_ldl, na.rm = TRUE), r))
} else {
  cat(sprintf("  Log(ApoB/LDL): N=%d, mean=%.3f (SD=%.3f)\n",
              sum(!is.na(df$log_apob_ldl)), mean(df$log_apob_ldl, na.rm = TRUE),
              sd(df$log_apob_ldl, na.rm = TRUE)))
}

# If fresh computation has more missing, fall back to existing
if (sum(!is.na(df$log_apob_ldl_existing)) > sum(!is.na(df$log_apob_ldl))) {
  df$log_apob_ldl <- df$log_apob_ldl_existing
  cat("  --> Using existing LOG_APOB_LDL (better completeness)\n")
}

# --- Lp(a) Binary (>125 nmol/L) ---
df$lpa_raw <- as.numeric(df$Lpa)
df$lpa_binary <- as.numeric(df$lpa_raw > 125)
cat(sprintf("  Lp(a): N=%d, median=%.0f, elevated >125: %d (%.1f%%)\n",
            sum(!is.na(df$lpa_raw)), median(df$lpa_raw, na.rm = TRUE),
            sum(df$lpa_binary == 1, na.rm = TRUE),
            100 * mean(df$lpa_binary, na.rm = TRUE)))

# --- Inverse HDL-C ---
df$hdl_raw <- as.numeric(df$HDL_1)
df$inv_hdl <- 1 / df$hdl_raw
cat(sprintf("  HDL-C: N=%d, mean=%.2f (SD=%.2f) → inv_hdl mean=%.3f\n",
            sum(!is.na(df$hdl_raw)), mean(df$hdl_raw, na.rm = TRUE),
            sd(df$hdl_raw, na.rm = TRUE), mean(df$inv_hdl, na.rm = TRUE)))

# --- Additional variables for R-based Model 4 ---
df$IsIndex <- as.numeric(df$I_vs_R)
df$Xanthoma <- as.numeric(df$TendonXanthomata)
df$Xanthoma[is.na(df$Xanthoma)] <- 0
df$Arcus <- as.numeric(df$CornealArcus)
df$Arcus[is.na(df$Arcus)] <- 0
df$Xanthelasma <- as.numeric(df$Xanthelasmas)
df$Xanthelasma[is.na(df$Xanthelasma)] <- 0
df$N_Stigmata <- (df$Arcus > 0) + (df$Xanthoma > 0) + (df$Xanthelasma > 0)
df$Ratio_High <- as.numeric(as.numeric(df$APOB_LDL) >= 0.31)

# --- SAFEHEART-RE variables ---
df$male_sex <- as.numeric(df$GENDER_M)
df$htn <- as.numeric(df$HTN_M)
df$smoking <- as.numeric(df$Smoking_binary)
df$bmi <- 27.5  # Not directly available; will use population mean
df$prior_cvd <- as.numeric(df$ASCVD_combined == 1)  # prevalent ASCVD (for cross-ref)

# =============================================================================
# 3. COMPLETE-CASE ANALYSIS SETS
# =============================================================================

cat("\n--- Complete-Case Analysis Sets ---\n\n")

# CALON 6-variable model (manuscript specification)
calon_vars <- c("age", "re_ldl", "lipid_years", "log_apob_ldl", "lpa_binary", "inv_hdl", "ASCVD")
df_calon <- df[complete.cases(df[, calon_vars]), ]
cat(sprintf("  CALON 6-variable: N=%d (events=%d, EPV=%.1f)\n",
            nrow(df_calon), sum(df_calon$ASCVD), sum(df_calon$ASCVD) / 6))

# CALON-CORE (without ApoB/LDL ratio and Lp(a) — variables unavailable in many UKB FH)
core_vars <- c("age", "re_ldl", "lipid_years", "inv_hdl", "ASCVD")
df_core <- df[complete.cases(df[, core_vars]), ]
cat(sprintf("  CALON-CORE 4-variable: N=%d (events=%d, EPV=%.1f)\n",
            nrow(df_core), sum(df_core$ASCVD), sum(df_core$ASCVD) / 4))

# R-based Model 4 (from CALON_FH_ANALYSIS_VERIFIED.R)
model4_vars <- c("age_tx", "IsIndex", "N_Stigmata", "lpa_binary", "Ratio_High", "ASCVD")
df_model4 <- df[complete.cases(df[, model4_vars]), ]
cat(sprintf("  R-Model 4 (Age_Tx+Index+Stigmata+Lpa+Ratio): N=%d (events=%d)\n",
            nrow(df_model4), sum(df_model4$ASCVD)))

# SAFEHEART-RE variables
sh_vars <- c("age", "male_sex", "LDL_1_raw", "htn", "smoking", "ASCVD")
df_sh <- df[complete.cases(df[, sh_vars]), ]
cat(sprintf("  SAFEHEART-RE complete: N=%d (events=%d)\n",
            nrow(df_sh), sum(df_sh$ASCVD)))

# =============================================================================
# 4. FIT CALON 6-VARIABLE MODEL (Manuscript Specification)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  MODEL A: CALON 6-VARIABLE (Manuscript Specification)          \n")
cat("================================================================\n\n")

# Fit logistic regression on REAL data
calon_model <- glm(ASCVD ~ age + re_ldl + lipid_years + log_apob_ldl + lpa_binary + inv_hdl,
                   data = df_calon, family = binomial)

# Print summary
cat("--- Model Summary ---\n")
print(summary(calon_model))

# Extract coefficients
calon_coefs <- coef(calon_model)
calon_se <- summary(calon_model)$coefficients[, 2]
calon_p <- summary(calon_model)$coefficients[, 4]
calon_or <- exp(calon_coefs)

cat("\n--- CALON Frozen Coefficients ---\n")
cat(sprintf("%-20s %12s %8s %10s %10s\n", "Variable", "Coefficient", "SE", "OR", "P-value"))
cat(paste(rep("-", 65), collapse = ""), "\n")
for (i in 1:length(calon_coefs)) {
  cat(sprintf("%-20s %12.6f %8.4f %10.4f %10.6f\n",
              names(calon_coefs)[i], calon_coefs[i], calon_se[i], calon_or[i], calon_p[i]))
}

# AUC
roc_calon <- roc(df_calon$ASCVD, predict(calon_model, type = "response"), quiet = TRUE)
auc_calon <- as.numeric(roc_calon$auc)
auc_ci <- ci.auc(roc_calon)

cat(sprintf("\nAUC: %.4f (95%% CI: %.4f - %.4f)\n", auc_calon, auc_ci[1], auc_ci[3]))

# Cross-validated AUC (5-fold)
cv_aucs <- numeric(5)
folds <- sample(rep(1:5, length.out = nrow(df_calon)))
for (k in 1:5) {
  train <- df_calon[folds != k, ]
  test <- df_calon[folds == k, ]
  if (sum(test$ASCVD) < 2 || sum(test$ASCVD) == nrow(test)) next
  m <- glm(ASCVD ~ age + re_ldl + lipid_years + log_apob_ldl + lpa_binary + inv_hdl,
           data = train, family = binomial)
  p <- predict(m, newdata = test, type = "response")
  cv_aucs[k] <- as.numeric(roc(test$ASCVD, p, quiet = TRUE)$auc)
}
cv_auc <- mean(cv_aucs[cv_aucs > 0])
optimism <- auc_calon - cv_auc
cat(sprintf("CV-AUC: %.4f (optimism: %.4f)\n", cv_auc, optimism))

# Feature importance (absolute standardised coefficients)
cat("\n--- Feature Importance (Relative Contribution) ---\n")
X_calon <- as.matrix(df_calon[, c("age", "re_ldl", "lipid_years", "log_apob_ldl", "lpa_binary", "inv_hdl")])
sds <- apply(X_calon, 2, sd)
abs_scaled <- abs(calon_coefs[-1]) * sds  # exclude intercept
contributions <- abs_scaled / sum(abs_scaled) * 100

for (i in order(contributions, decreasing = TRUE)) {
  cat(sprintf("  %-20s: %5.1f%%\n", names(contributions)[i], contributions[i]))
}
cat(sprintf("\n  Hidden Burdens (RE-LDL + Lipid-years + Log(ApoB/LDL)): %.1f%%\n",
            contributions["re_ldl"] + contributions["lipid_years"] + contributions["log_apob_ldl"]))

# =============================================================================
# 5. FIT CALON-CORE 4-VARIABLE MODEL (Without ApoB Ratio and Lp(a))
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  MODEL B: CALON-CORE 4-VARIABLE (For UKB Tier 2)              \n")
cat("================================================================\n\n")

core_model <- glm(ASCVD ~ age + re_ldl + lipid_years + inv_hdl,
                  data = df_core, family = binomial)

print(summary(core_model))

core_coefs <- coef(core_model)
roc_core <- roc(df_core$ASCVD, predict(core_model, type = "response"), quiet = TRUE)
auc_core <- as.numeric(roc_core$auc)
auc_core_ci <- ci.auc(roc_core)

cat(sprintf("\nCALON-CORE AUC: %.4f (95%% CI: %.4f - %.4f)\n",
            auc_core, auc_core_ci[1], auc_core_ci[3]))

# =============================================================================
# 6. FIT R-BASED MODEL 4 (From CALON_FH_ANALYSIS_VERIFIED.R)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  MODEL C: R-Based Model 4 (Age_Tx + Index + Stigmata + Lpa +  \n")
cat("           Ratio_High) - Alternative Specification              \n")
cat("================================================================\n\n")

model4 <- glm(ASCVD ~ age_tx + IsIndex + N_Stigmata + lpa_binary + Ratio_High,
              data = df_model4, family = binomial)

print(summary(model4))

model4_coefs <- coef(model4)
roc_model4 <- roc(df_model4$ASCVD, predict(model4, type = "response"), quiet = TRUE)
auc_model4 <- as.numeric(roc_model4$auc)

cat(sprintf("\nModel 4 AUC: %.4f\n", auc_model4))

# =============================================================================
# 7. SAFEHEART-RE (Published Frozen Coefficients - Applied to Wales Data)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  MODEL D: SAFEHEART-RE (Published Coefficients on Wales Data)  \n")
cat("================================================================\n\n")

# Published SAFEHEART-RE coefficients
sh_intercept <- -7.053
sh_age <- 0.064
sh_male <- 0.775
sh_ldl <- 0.109
sh_htn <- 0.431
sh_bmi <- 0.025
sh_smoking <- 0.466
sh_prior_cvd <- 1.414

# Apply to Wales data (note: BMI not in DRAGON2, use 27.5; prior_cvd set to 0 for incident prediction)
df_sh$sh_lp <- sh_intercept + sh_age * df_sh$age + sh_male * df_sh$male_sex +
               sh_ldl * df_sh$LDL_1_raw + sh_htn * df_sh$htn + sh_bmi * 27.5 +
               sh_smoking * df_sh$smoking
df_sh$sh_prob <- 1 / (1 + exp(-df_sh$sh_lp))

roc_sh <- roc(df_sh$ASCVD, df_sh$sh_prob, quiet = TRUE)
auc_sh <- as.numeric(roc_sh$auc)
auc_sh_ci <- ci.auc(roc_sh)

cat(sprintf("SAFEHEART-RE AUC on Wales data: %.4f (95%% CI: %.4f - %.4f)\n",
            auc_sh, auc_sh_ci[1], auc_sh_ci[3]))

# =============================================================================
# 8. HEAD-TO-HEAD COMPARISON (DeLong Test)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  HEAD-TO-HEAD COMPARISON                                        \n")
cat("================================================================\n\n")

# Align datasets - need common patients with all models
common_ids <- intersect(rownames(df_calon), rownames(df_sh))
if (length(common_ids) > 10) {
  # Use patient IDs to align
  df_compare <- df_calon[!is.na(df_calon$male_sex) & !is.na(df_calon$LDL_1_raw) &
                         !is.na(df_calon$htn) & !is.na(df_calon$smoking), ]

  # CALON predictions
  compare_calon_prob <- predict(calon_model, newdata = df_compare, type = "response")

  # SAFEHEART-RE predictions
  compare_sh_lp <- sh_intercept + sh_age * df_compare$age + sh_male * df_compare$male_sex +
                   sh_ldl * df_compare$LDL_1_raw + sh_htn * df_compare$htn + sh_bmi * 27.5 +
                   sh_smoking * df_compare$smoking
  compare_sh_prob <- 1 / (1 + exp(-compare_sh_lp))

  roc1 <- roc(df_compare$ASCVD, compare_calon_prob, quiet = TRUE)
  roc2 <- roc(df_compare$ASCVD, compare_sh_prob, quiet = TRUE)

  delong <- roc.test(roc1, roc2, method = "delong")

  cat(sprintf("On matched subset (N=%d, events=%d):\n", nrow(df_compare), sum(df_compare$ASCVD)))
  cat(sprintf("  CALON AUC:       %.4f\n", as.numeric(roc1$auc)))
  cat(sprintf("  SAFEHEART-RE AUC: %.4f\n", as.numeric(roc2$auc)))
  cat(sprintf("  DELTAUC:          %.4f\n", as.numeric(roc1$auc) - as.numeric(roc2$auc)))
  cat(sprintf("  DeLong Z:        %.3f\n", delong$statistic))
  cat(sprintf("  DeLong p:        %.6f\n", delong$p.value))
} else {
  cat("  (Insufficient overlapping patients for DeLong test)\n")
}

# =============================================================================
# 9. CALIBRATION ASSESSMENT
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  CALIBRATION ASSESSMENT (CALON Development)                     \n")
cat("================================================================\n\n")

calon_pred <- predict(calon_model, type = "response")

# O/E ratio
oe_ratio <- sum(df_calon$ASCVD) / sum(calon_pred)
cat(sprintf("O/E ratio: %.3f\n", oe_ratio))

# Calibration slope and intercept
lp <- predict(calon_model, type = "link")
cal_model <- glm(ASCVD ~ lp, data = df_calon, family = binomial)
cal_slope <- coef(cal_model)[2]
cal_intercept <- coef(cal_model)[1]
cat(sprintf("Calibration slope: %.3f (ideal=1)\n", cal_slope))
cat(sprintf("Calibration intercept: %.3f (ideal=0)\n", cal_intercept))

# Hosmer-Lemeshow
n_groups <- 10
sorted_idx <- order(calon_pred)
groups <- split(sorted_idx, cut(seq_along(sorted_idx), n_groups, labels = FALSE))

hl_chi2 <- 0
cat(sprintf("\n%-8s %6s %10s %10s %8s\n", "Decile", "N", "Expected", "Observed", "O/E"))
cat(paste(rep("-", 46), collapse = ""), "\n")
for (i in 1:n_groups) {
  g <- groups[[i]]
  obs <- sum(df_calon$ASCVD[g])
  exp <- sum(calon_pred[g])
  n_g <- length(g)
  oe_g <- if (exp > 0) obs / exp else NA
  if (exp > 0 && exp < n_g) {
    hl_chi2 <- hl_chi2 + (obs - exp)^2 / (exp * (1 - exp / n_g))
  }
  cat(sprintf("%-8d %6d %10.2f %10d %8.2f\n", i, n_g, exp, obs, oe_g))
}
hl_p <- 1 - pchisq(hl_chi2, n_groups - 2)
cat(sprintf("\nHosmer-Lemeshow: chi2=%.2f, df=%d, p=%.4f\n", hl_chi2, n_groups - 2, hl_p))

# Brier score
brier <- mean((calon_pred - df_calon$ASCVD)^2)
brier_null <- mean((mean(df_calon$ASCVD) - df_calon$ASCVD)^2)
brier_scaled <- 1 - brier / brier_null
cat(sprintf("Brier score: %.4f (scaled: %.4f)\n", brier, brier_scaled))

# =============================================================================
# 10. SAVE FROZEN COEFFICIENTS AND STANDARDISATION PARAMETERS
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SAVING FROZEN COEFFICIENTS                                     \n")
cat("================================================================\n\n")

OUTPUT_DIR <- "C:/Users/nader/Downloads/calon_ukb_pipeline/"

# --- A: CALON 6-Variable Coefficients ---
calon_coef_df <- data.frame(
  variable = names(calon_coefs),
  coefficient = as.numeric(calon_coefs),
  se = as.numeric(calon_se),
  odds_ratio = as.numeric(calon_or),
  p_value = as.numeric(calon_p),
  stringsAsFactors = FALSE
)
write.csv(calon_coef_df, paste0(OUTPUT_DIR, "CALON_frozen_coefficients.csv"), row.names = FALSE)
cat("  Saved: CALON_frozen_coefficients.csv\n")

# --- B: CALON-CORE 4-Variable Coefficients ---
core_coef_df <- data.frame(
  variable = names(core_coefs),
  coefficient = as.numeric(core_coefs),
  se = as.numeric(summary(core_model)$coefficients[, 2]),
  odds_ratio = as.numeric(exp(core_coefs)),
  p_value = as.numeric(summary(core_model)$coefficients[, 4]),
  stringsAsFactors = FALSE
)
write.csv(core_coef_df, paste0(OUTPUT_DIR, "CALON_CORE_frozen_coefficients.csv"), row.names = FALSE)
cat("  Saved: CALON_CORE_frozen_coefficients.csv\n")

# --- C: R-Based Model 4 Coefficients ---
model4_coef_df <- data.frame(
  variable = names(model4_coefs),
  coefficient = as.numeric(model4_coefs),
  se = as.numeric(summary(model4)$coefficients[, 2]),
  odds_ratio = as.numeric(exp(model4_coefs)),
  p_value = as.numeric(summary(model4)$coefficients[, 4]),
  stringsAsFactors = FALSE
)
write.csv(model4_coef_df, paste0(OUTPUT_DIR, "CALON_Model4_frozen_coefficients.csv"), row.names = FALSE)
cat("  Saved: CALON_Model4_frozen_coefficients.csv\n")

# --- D: Variable Means and SDs (for standardisation if needed) ---
means_sds <- data.frame(
  variable = c("age", "re_ldl", "lipid_years", "log_apob_ldl", "lpa_binary", "inv_hdl"),
  mean = c(mean(df_calon$age), mean(df_calon$re_ldl), mean(df_calon$lipid_years),
           mean(df_calon$log_apob_ldl), mean(df_calon$lpa_binary), mean(df_calon$inv_hdl)),
  sd = c(sd(df_calon$age), sd(df_calon$re_ldl), sd(df_calon$lipid_years),
         sd(df_calon$log_apob_ldl), sd(df_calon$lpa_binary), sd(df_calon$inv_hdl)),
  median = c(median(df_calon$age), median(df_calon$re_ldl), median(df_calon$lipid_years),
             median(df_calon$log_apob_ldl), median(df_calon$lpa_binary), median(df_calon$inv_hdl)),
  min = c(min(df_calon$age), min(df_calon$re_ldl), min(df_calon$lipid_years),
          min(df_calon$log_apob_ldl), min(df_calon$lpa_binary), min(df_calon$inv_hdl)),
  max = c(max(df_calon$age), max(df_calon$re_ldl), max(df_calon$lipid_years),
          max(df_calon$log_apob_ldl), max(df_calon$lpa_binary), max(df_calon$inv_hdl)),
  n_complete = c(sum(!is.na(df_calon$age)), sum(!is.na(df_calon$re_ldl)),
                 sum(!is.na(df_calon$lipid_years)), sum(!is.na(df_calon$log_apob_ldl)),
                 sum(!is.na(df_calon$lpa_binary)), sum(!is.na(df_calon$inv_hdl))),
  stringsAsFactors = FALSE
)
write.csv(means_sds, paste0(OUTPUT_DIR, "CALON_development_variable_summary.csv"), row.names = FALSE)
cat("  Saved: CALON_development_variable_summary.csv\n")

# --- E: Model Performance Summary ---
performance_df <- data.frame(
  model = c("CALON_6var", "CALON_CORE_4var", "R_Model4_5var", "SAFEHEART_RE"),
  n_patients = c(nrow(df_calon), nrow(df_core), nrow(df_model4), nrow(df_sh)),
  n_events = c(sum(df_calon$ASCVD), sum(df_core$ASCVD), sum(df_model4$ASCVD), sum(df_sh$ASCVD)),
  auc = c(auc_calon, auc_core, auc_model4, auc_sh),
  auc_lower = c(auc_ci[1], auc_core_ci[1], NA, auc_sh_ci[1]),
  auc_upper = c(auc_ci[3], auc_core_ci[3], NA, auc_sh_ci[3]),
  n_variables = c(6, 4, 5, 7),
  stringsAsFactors = FALSE
)
write.csv(performance_df, paste0(OUTPUT_DIR, "CALON_model_performance_development.csv"), row.names = FALSE)
cat("  Saved: CALON_model_performance_development.csv\n")

# =============================================================================
# 11. GENERATE R CODE SNIPPET FOR 02_CALON_validate.R
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  R CODE SNIPPET FOR 02_CALON_validate.R                        \n")
cat("  Copy-paste these EXACT values to replace PLACEHOLDERs         \n")
cat("================================================================\n\n")

cat("# ====== FROZEN CALON COEFFICIENTS (From 99_extract_frozen_coefficients.R) ======\n")
cat("# Model: glm(ASCVD ~ age + re_ldl + lipid_years + log_apob_ldl + lpa_binary + inv_hdl)\n")
cat(sprintf("# Fitted on: Wales FH Registry (DRAGON2.csv), N=%d, events=%d\n", nrow(df_calon), sum(df_calon$ASCVD)))
cat(sprintf("# Development AUC: %.4f (95%% CI: %.4f - %.4f)\n", auc_calon, auc_ci[1], auc_ci[3]))
cat(sprintf("# CV-AUC: %.4f\n", cv_auc))
cat("\n")
cat("calon_coefficients <- list(\n")
cat(sprintf("  intercept     = %.10f,\n", calon_coefs["(Intercept)"]))
cat(sprintf("  age           = %.10f,\n", calon_coefs["age"]))
cat(sprintf("  re_ldl        = %.10f,\n", calon_coefs["re_ldl"]))
cat(sprintf("  lipid_years   = %.10f,\n", calon_coefs["lipid_years"]))
cat(sprintf("  log_apob_ldl  = %.10f,\n", calon_coefs["log_apob_ldl"]))
cat(sprintf("  lpa_binary    = %.10f,\n", calon_coefs["lpa_binary"]))
cat(sprintf("  inv_hdl       = %.10f\n",  calon_coefs["inv_hdl"]))
cat(")\n\n")

cat("# ====== FROZEN CALON-CORE COEFFICIENTS ======\n")
cat(sprintf("# Fitted on: Wales FH Registry, N=%d, events=%d\n", nrow(df_core), sum(df_core$ASCVD)))
cat(sprintf("# Development AUC: %.4f\n", auc_core))
cat("\n")
cat("calon_core_coefficients <- list(\n")
cat(sprintf("  intercept     = %.10f,\n", core_coefs["(Intercept)"]))
cat(sprintf("  age           = %.10f,\n", core_coefs["age"]))
cat(sprintf("  re_ldl        = %.10f,\n", core_coefs["re_ldl"]))
cat(sprintf("  lipid_years   = %.10f,\n", core_coefs["lipid_years"]))
cat(sprintf("  inv_hdl       = %.10f\n",  core_coefs["inv_hdl"]))
cat(")\n\n")

# =============================================================================
# 12. VARIABLE COMPLETENESS MATRIX (For Planning)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  VARIABLE COMPLETENESS IN DEVELOPMENT DATA                      \n")
cat("================================================================\n\n")

all_vars <- c("age", "re_ldl", "age_tx", "lipid_years", "log_apob_ldl",
              "lpa_binary", "inv_hdl", "IsIndex", "N_Stigmata", "Ratio_High",
              "male_sex", "htn", "smoking", "ASCVD")

for (v in all_vars) {
  n_total <- nrow(df)
  n_complete <- sum(!is.na(df[[v]]))
  pct <- 100 * n_complete / n_total
  cat(sprintf("  %-20s: %d/%d (%.1f%%)\n", v, n_complete, n_total, pct))
}

# =============================================================================
# 13. EXPECTED UKB VARIABLE AVAILABILITY
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  EXPECTED UKB VARIABLE AVAILABILITY                             \n")
cat("================================================================\n\n")

cat("  Variable          | UKB Field       | Available? | Notes\n")
cat("  ------------------|-----------------|------------|------\n")
cat("  age               | p21022          | YES        | Direct\n")
cat("  re_ldl            | p30780/p20003   | YES        | LDL/RF from statin codes\n")
cat("  age_tx            | NOT RECORDED    | PROXY      | 3 scenarios (A/B/C)\n")
cat("  lipid_years       | Derived         | YES*       | Uses proxy age_tx\n")
cat("  log_apob_ldl      | p30640/p30780   | YES        | log(ApoB/LDL)\n")
cat("  lpa_binary        | p30790          | YES        | >125 nmol/L\n")
cat("  inv_hdl           | p30760          | YES        | 1/HDL\n")
cat("  IsIndex           | NOT AVAILABLE   | NO         | All UKB are probands\n")
cat("  N_Stigmata        | NOT AVAILABLE   | NO         | Not recorded in UKB\n")
cat("  Ratio_High        | p30640/p30780   | YES        | ApoB/LDL >= 0.31\n")
cat("  male_sex          | p31             | YES        | Direct\n")
cat("  htn               | p20002/p4080    | YES        | Composite\n")
cat("  smoking           | p20116          | YES        | Code 2 = current\n")
cat("  ASCVD             | p131298/etc     | YES        | First-occurrence dates\n")
cat("\n")
cat("  * Lipid-years requires age_tx proxy; 3 sensitivity scenarios used\n")
cat("  CRITICAL: IsIndex and N_Stigmata are NOT available in UKB.\n")
cat("  This means Model 4 (R-based) CANNOT be validated in UKB.\n")
cat("  Only the 6-variable CALON model and CALON-CORE can be validated.\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n\n")
cat("================================================================\n")
cat("  SUMMARY: FROZEN COEFFICIENTS FOR UKB VALIDATION               \n")
cat("================================================================\n\n")

cat(sprintf("  CALON 6-variable:  AUC=%.4f (N=%d, %d events)\n", auc_calon, nrow(df_calon), sum(df_calon$ASCVD)))
cat(sprintf("  CALON-CORE 4-var:  AUC=%.4f (N=%d, %d events)\n", auc_core, nrow(df_core), sum(df_core$ASCVD)))
cat(sprintf("  R-Model 4 5-var:   AUC=%.4f (N=%d, %d events) [NOT validatable in UKB]\n",
            auc_model4, nrow(df_model4), sum(df_model4$ASCVD)))
cat(sprintf("  SAFEHEART-RE:      AUC=%.4f (N=%d, %d events)\n\n", auc_sh, nrow(df_sh), sum(df_sh$ASCVD)))

cat("  Output files saved to: ", OUTPUT_DIR, "\n")
cat("  1. CALON_frozen_coefficients.csv\n")
cat("  2. CALON_CORE_frozen_coefficients.csv\n")
cat("  3. CALON_Model4_frozen_coefficients.csv\n")
cat("  4. CALON_development_variable_summary.csv\n")
cat("  5. CALON_model_performance_development.csv\n")
cat("\n  >> Next step: Copy coefficients to 02_CALON_validate.R\n")
cat("  >> Then run pipeline on UKB RAP\n")

cat("\n================================================================\n")
cat("  EXTRACTION COMPLETE                                            \n")
cat("================================================================\n")

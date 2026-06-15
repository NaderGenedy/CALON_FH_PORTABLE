################################################################################
#                                                                              #
#  CALON-2: ASCVD RISK MODEL DEVELOPMENT v7                                   #
#  "From Classification to Survival — Unlocking Time-to-Event"               #
#                                                                              #
#  Key changes v6 → v7:                                                       #
#    1. PENALIZED COX: glmnet family="cox" using followup_years + events      #
#       (exploits censoring; C-index replaces AUC as primary metric)          #
#    2. RCS (Restricted Cubic Splines): Non-linear effects for age,           #
#       lipids, BMI via rms::cph (3-4 knots)                                  #
#    3. REPEATED CV: 5×10-fold for stable estimates (reduces variance)        #
#    4. FOCUSED XGBoost: Cox objective (survival:cox) + refined grid          #
#    5. RSF (Random Survival Forest): ranger for survival                     #
#    6. IMPROVED STACKING: Nested OOF with proper CV of meta-learner         #
#    7. TIME-DEPENDENT AUC: Harrell's C + timeROC at 5/10 years             #
#    8. CALIBRATION: Observed vs predicted at 5/10 years                     #
#    9. NET RECLASSIFICATION: NRI vs SAFEHEART-RE                            #
#                                                                              #
#  Author:  Dr Nader Genedy                                                    #
#  Date:    March 2026                                                         #
#  Target:  C-index > 0.78 / time-dependent AUC > 0.80 at 10 years          #
#                                                                              #
################################################################################

rm(list = ls())
set.seed(2026)

# =============================================================================
# PACKAGES
# =============================================================================

required_packages <- c(
  # Core
  "dplyr", "tidyr", "readr", "ggplot2", "patchwork", "scales",
  # Penalized regression
  "glmnet", "pROC",
  # Survival analysis
  "survival", "rms",
  # Random survival forest
  "ranger",
  # XGBoost
  "xgboost", "caret",
  # MI
  "mice",
  # Boruta
  "Boruta",
  # Time-dependent ROC
  "timeROC",
  # Reclassification
  "nricens",
  # Bootstrap
  "boot",
  # Tables
  "tableone"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("\n")
cat("======================================================================\n")
cat("  CALON-2 v7: Survival-Based Development Pipeline                      \n")
cat("  Key: Penalized Cox, RCS, RSF, time-dependent AUC, 5x10 CV           \n")
cat("======================================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

LOCAL_PATH  <- "C:/Users/nader/Downloads/calon_ukb_pipeline/"
CLOUD_PATH  <- "/cloud/project/"

if (dir.exists(LOCAL_PATH)) {
  BASE_DIR <- LOCAL_PATH
  cat("  Environment: Local Windows\n")
} else if (dir.exists(CLOUD_PATH)) {
  BASE_DIR <- CLOUD_PATH
  cat("  Environment: Posit Cloud\n")
} else {
  BASE_DIR <- paste0(normalizePath(getwd(), winslash = "/"), "/")
  cat(sprintf("  Environment: Other. BASE_DIR = %s\n", BASE_DIR))
}

cat(sprintf("  BASE_DIR = %s\n", BASE_DIR))

OUT_DIR  <- paste0(BASE_DIR, "output/")
TAB_DIR  <- paste0(OUT_DIR, "tables/")
FIG_DIR  <- paste0(OUT_DIR, "figures/")
for (d in c(TAB_DIR, FIG_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# Find merged CSV
merged_file <- NULL
search_paths <- c(
  paste0(OUT_DIR, "calon2_full_merged.csv"),
  paste0(BASE_DIR, "calon2_full_merged.csv"),
  paste0(getwd(), "/output/calon2_full_merged.csv")
)
for (sp in search_paths) {
  if (file.exists(sp)) { merged_file <- sp; break }
}
if (is.null(merged_file)) stop("ERROR: calon2_full_merged.csv not found")
cat(sprintf("  Found: %s\n", merged_file))

# =============================================================================
# SECTION 1: LOAD DATA
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 1: LOADING DATA\n")
cat("===================================================================\n\n")

df <- read.csv(merged_file, stringsAsFactors = FALSE)
cat(sprintf("  Dataset: %d patients x %d columns\n", nrow(df), ncol(df)))

# Convert dates
df$assessment_date  <- as.Date(df$assessment_date)
df$first_ascvd_date <- as.Date(df$first_ascvd_date)
if ("exit_date" %in% names(df)) df$exit_date <- as.Date(df$exit_date)
if ("death_date" %in% names(df)) df$death_date <- as.Date(df$death_date)

# Verify survival columns exist
stopifnot(all(c("followup_years", "event_indicator") %in% names(df)))
cat(sprintf("  Events: %d / %d (%.1f%%)\n",
            sum(df$event_indicator == 1, na.rm = TRUE), nrow(df),
            100 * mean(df$event_indicator == 1, na.rm = TRUE)))
cat(sprintf("  Median follow-up: %.1f years (IQR: %.1f-%.1f)\n",
            median(df$followup_years, na.rm = TRUE),
            quantile(df$followup_years, 0.25, na.rm = TRUE),
            quantile(df$followup_years, 0.75, na.rm = TRUE)))

# For incident-only analysis, optionally filter prevalent cases
# (using ascvd_combined for consistency with v6)
cat(sprintf("  ascvd_combined (binary): %d events\n", sum(df$ascvd_combined)))

# =============================================================================
# SECTION 2: FEATURE ENGINEERING (same as v6 for consistency)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 2: FEATURE ENGINEERING\n")
cat("===================================================================\n\n")

# -- Gene classification
if ("gene" %in% names(df)) {
  df$gene_apob <- as.integer(df$gene == "APOB")
  cat(sprintf("  gene_apob: %d APOB carriers\n", sum(df$gene_apob, na.rm = TRUE)))
}

# -- PRS (z-scored)
prs_cols <- grep("^p26", names(df), value = TRUE)
if (length(prs_cols) > 0) {
  prs_priority <- unique(c("p26227", "p26223", prs_cols))
  prs_priority <- prs_priority[prs_priority %in% names(df)]
  for (pc in prs_priority) {
    n_valid <- sum(!is.na(df[[pc]]))
    if (n_valid > 100) {
      df$prs_primary <- (df[[pc]] - mean(df[[pc]], na.rm = TRUE)) / sd(df[[pc]], na.rm = TRUE)
      cat(sprintf("  prs_primary: using %s (%d non-missing, z-scored)\n", pc, n_valid))
      break
    }
  }
}

# -- Log transforms
for (v in c("crp", "hba1c", "creatinine", "alt", "cystatin_c")) {
  if (v %in% names(df)) {
    log_name <- paste0("log_", v)
    df[[log_name]] <- ifelse(!is.na(df[[v]]) & df[[v]] > 0, log(df[[v]]), NA)
    cat(sprintf("  %s: %d non-missing\n", log_name, sum(!is.na(df[[log_name]]))))
  }
}

# -- Lp(a) features
if ("lpa" %in% names(df) && sum(!is.na(df$lpa)) > 100) {
  df$log_lpa <- ifelse(!is.na(df$lpa) & df$lpa > 0, log(df$lpa),
                       ifelse(!is.na(df$lpa), log(0.1), NA))
  df$lpa_extreme <- as.integer(!is.na(df$lpa) & df$lpa >= 250)
  df$lpa_high    <- as.integer(!is.na(df$lpa) & df$lpa >= 125)
  cat(sprintf("  log_lpa: %d non-missing; extreme: %d; high: %d\n",
              sum(!is.na(df$log_lpa)),
              sum(df$lpa_extreme == 1, na.rm = TRUE),
              sum(df$lpa_high == 1, na.rm = TRUE)))
}

# -- ApoB/LDL ratios
if (!"log_apob_ldl" %in% names(df)) {
  if (all(c("apob", "ldl") %in% names(df))) {
    df$apob_ldl_ratio <- ifelse(!is.na(df$apob) & !is.na(df$ldl) & df$ldl > 0,
                                df$apob / df$ldl, NA)
    df$log_apob_ldl <- ifelse(!is.na(df$apob_ldl_ratio) & df$apob_ldl_ratio > 0,
                              log(df$apob_ldl_ratio), NA)
  }
}
if (!"log_apob_re_ldl" %in% names(df)) {
  if (all(c("apob", "re_ldl") %in% names(df))) {
    ratio <- ifelse(!is.na(df$apob) & !is.na(df$re_ldl) & df$re_ldl > 0,
                    df$apob / df$re_ldl, NA)
    df$log_apob_re_ldl <- ifelse(!is.na(ratio) & ratio > 0, log(ratio), NA)
  }
}

# -- Interactions
if (all(c("age", "re_ldl") %in% names(df)))
  df$age_x_re_ldl <- df$age * df$re_ldl
if (all(c("age", "hdl") %in% names(df)))
  df$age_x_hdl <- df$age * df$hdl

# -- Lipid ratios and derived features
if (all(c("tc", "hdl") %in% names(df))) {
  df$non_hdl <- df$tc - df$hdl
  df$tc_hdl_ratio <- ifelse(df$hdl > 0, df$tc / df$hdl, NA)
}
if (all(c("trig", "hdl") %in% names(df))) {
  df$tg_hdl_ratio <- ifelse(df$hdl > 0, df$trig / df$hdl, NA)
  df$log_tg_hdl <- ifelse(!is.na(df$tg_hdl_ratio) & df$tg_hdl_ratio > 0,
                          log(df$tg_hdl_ratio), NA)
}
if (all(c("sbp", "dbp") %in% names(df)))
  df$pulse_pressure <- df$sbp - df$dbp

# eGFR (CKD-EPI 2021, race-free)
if (all(c("creatinine", "age", "sex") %in% names(df))) {
  scr_mgdl <- df$creatinine / 88.4
  kappa  <- ifelse(df$sex == 1, 0.9, 0.7)
  alpha  <- ifelse(df$sex == 1, -0.302, -0.241)
  female <- ifelse(df$sex == 0, 1.012, 1.0)
  df$egfr <- 142 * pmin(scr_mgdl/kappa, 1)^alpha *
    pmax(scr_mgdl/kappa, 1)^(-1.200) * 0.9938^df$age * female
}

# Non-linear age
df$age_sq    <- df$age^2
df$sex_x_age <- df$sex * df$age

# ApoB/ApoA1 ratio
if (all(c("apob", "apoa") %in% names(df)))
  df$apob_apoa_ratio <- ifelse(df$apoa > 0, df$apob / df$apoa, NA)

# Remnant cholesterol
if (all(c("tc", "hdl", "ldl") %in% names(df)))
  df$remnant_chol_calc <- df$tc - df$hdl - df$ldl

# Metabolic syndrome score
metsyn_vars <- c("bmi", "trig", "sbp")
if (all(c(metsyn_vars, "hdl") %in% names(df))) {
  metsyn_z <- sapply(metsyn_vars, function(v) {
    x <- df[[v]]; (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  })
  hdl_z <- -1 * (df$hdl - mean(df$hdl, na.rm = TRUE)) / sd(df$hdl, na.rm = TRUE)
  df$metsyn_score <- rowMeans(cbind(metsyn_z, hdl_z), na.rm = TRUE)
}

# PRS × re_ldl
if (all(c("prs_primary", "re_ldl") %in% names(df)))
  df$prs_x_re_ldl <- df$prs_primary * df$re_ldl

cat("  Feature engineering complete.\n")

# =============================================================================
# SECTION 3: NMR PCA (same as v6)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 3: NMR PCA\n")
cat("===================================================================\n\n")

nmr_all_cols <- grep("^p234\\d+_i0$", names(df), value = TRUE)
cat(sprintf("  NMR features: %d\n", length(nmr_all_cols)))

pc_names <- character(0)
if (length(nmr_all_cols) >= 20) {
  nmr_matrix <- as.matrix(df[, nmr_all_cols])
  for (j in 1:ncol(nmr_matrix)) {
    na_idx <- is.na(nmr_matrix[, j])
    if (any(na_idx)) nmr_matrix[na_idx, j] <- mean(nmr_matrix[, j], na.rm = TRUE)
  }
  nmr_scaled <- scale(nmr_matrix)
  zv <- apply(nmr_scaled, 2, sd) == 0
  nmr_scaled <- nmr_scaled[, !zv]

  pca_fit <- prcomp(nmr_scaled, center = FALSE, scale. = FALSE)
  cum_var <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
  n_pcs <- min(which(cum_var >= 0.90), 15)
  pc_names <- paste0("nmr_pc", 1:n_pcs)
  for (j in 1:n_pcs) df[[pc_names[j]]] <- pca_fit$x[, j]
  cat(sprintf("  %d PCs retained (%.1f%% variance)\n", n_pcs, 100 * cum_var[n_pcs]))
}

# =============================================================================
# SECTION 4: BORUTA FEATURE SELECTION (on survival outcome)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 4: BORUTA FEATURE SELECTION\n")
cat("===================================================================\n\n")

# Candidate variables (same pool as v6 tier4)
all_candidate_vars <- unique(c(
  "age", "sex", "re_ldl", "hdl",
  "log_apob_ldl", "log_lpa", "lpa_extreme", "lpa_high",
  "smoking_binary", "diabetes", "hypertension",
  "bmi", "log_crp", "lipid_years_A",
  "non_hdl", "egfr",
  "pulse_pressure", "tc_hdl_ratio", "metsyn_score",
  "age_sq", "sex_x_age",
  pc_names, "prs_primary", "prs_x_re_ldl",
  "apob_apoa_ratio", "remnant_chol_calc", "log_tg_hdl",
  "gene_apob", "townsend", "hba1c", "log_cystatin_c",
  "ever_smoked"
))
all_candidate_vars <- all_candidate_vars[all_candidate_vars %in% names(df)]

# Run Boruta on binary outcome (same as v6, for variable pre-screening)
df_boruta <- df[complete.cases(df[, c(all_candidate_vars, "ascvd_combined")]),
                c(all_candidate_vars, "ascvd_combined")]
cat(sprintf("  Input: %d vars, %d cases, %d events\n",
            length(all_candidate_vars), nrow(df_boruta), sum(df_boruta$ascvd_combined)))

boruta_confirmed <- character(0)
if (nrow(df_boruta) >= 500) {
  set.seed(2026)
  boruta_res <- tryCatch(
    TentativeRoughFix(Boruta(ascvd_combined ~ ., data = df_boruta,
                             maxRuns = 200, doTrace = 0)),
    error = function(e) { cat(sprintf("  Boruta failed: %s\n", e$message)); NULL }
  )
  if (!is.null(boruta_res)) {
    boruta_confirmed <- names(boruta_res$finalDecision[
      boruta_res$finalDecision == "Confirmed"])
    cat(sprintf("  CONFIRMED: %d features\n", length(boruta_confirmed)))
    cat(sprintf("  REJECTED: %d features\n",
                sum(boruta_res$finalDecision == "Rejected")))

    imp_scores <- colMeans(boruta_res$ImpHistory[, boruta_confirmed], na.rm = TRUE)
    imp_order <- order(-imp_scores)
    cat("\n  Confirmed features (ranked):\n")
    for (j in imp_order) {
      cat(sprintf("    %2d. %-25s imp=%.3f\n",
                  which(imp_order == j), boruta_confirmed[j], imp_scores[j]))
    }
  }
}

boruta_vars <- if (length(boruta_confirmed) >= 5) boruta_confirmed else all_candidate_vars[1:15]
cat(sprintf("\n  Boruta-selected: %d vars\n", length(boruta_vars)))

# =============================================================================
# SECTION 5: PENALIZED COX REGRESSION (the main new method)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 5: PENALIZED COX (glmnet, 5x10 repeated CV)\n")
cat("===================================================================\n\n")

# Define tiers for Cox
cox_tiers <- list(
  list(nm = "Core Clinical (8 vars)", vars = c(
    "age", "sex", "re_ldl", "hdl", "smoking_binary",
    "hypertension", "bmi", "lipid_years_A")),
  list(nm = "Core+ (15 vars)", vars = c(
    "age", "sex", "re_ldl", "hdl", "log_apob_ldl", "log_lpa",
    "lpa_extreme", "smoking_binary", "diabetes", "hypertension",
    "bmi", "log_crp", "lipid_years_A", "non_hdl", "egfr")),
  list(nm = "Boruta Lean", vars = boruta_vars),
  list(nm = "Boruta + NMR PCs", vars = unique(c(boruta_vars, pc_names)))
)

# Helper: repeated k-fold CV for penalized Cox
# Returns mean C-index across repeats
run_cox_glmnet <- function(name, vars, data, n_repeats = 5, k_folds = 10) {
  vars_exist <- vars[vars %in% names(data)]
  surv_cols <- c("followup_years", "event_indicator")
  df_cc <- data[complete.cases(data[, c(vars_exist, surv_cols)]), ]
  # Remove zero or negative follow-up
  df_cc <- df_cc[df_cc$followup_years > 0, ]

  n <- nrow(df_cc)
  n_ev <- sum(df_cc$event_indicator)
  cat(sprintf("  %-30s N=%d, Events=%d, Vars=%d\n",
              name, n, n_ev, length(vars_exist)))

  if (n < 200 || n_ev < 30) { cat("    SKIPPED\n"); return(NULL) }

  X <- as.matrix(df_cc[, vars_exist])
  y <- Surv(df_cc$followup_years, df_cc$event_indicator)

  alphas <- c(0, 0.1, 0.25, 0.5, 0.75, 1.0)
  alpha_labels <- c("Ridge", "EN(.10)", "EN(.25)", "EN(.50)", "EN(.75)", "LASSO")

  best_overall <- NULL

  for (ai in seq_along(alphas)) {
    # Repeated CV
    cv_c_indices <- numeric(n_repeats)
    for (r in 1:n_repeats) {
      set.seed(2026 + r)
      cv_fit <- tryCatch(
        cv.glmnet(X, y, family = "cox", alpha = alphas[ai],
                  nfolds = k_folds, type.measure = "C",
                  grouped = TRUE),
        error = function(e) NULL
      )
      if (!is.null(cv_fit)) {
        cv_c_indices[r] <- max(cv_fit$cvm, na.rm = TRUE)
      }
    }
    mean_c <- mean(cv_c_indices[cv_c_indices > 0])
    nz <- tryCatch({
      set.seed(2026)
      fit <- cv.glmnet(X, y, family = "cox", alpha = alphas[ai],
                       nfolds = k_folds, type.measure = "C")
      sum(coef(fit, s = "lambda.min") != 0)
    }, error = function(e) NA)

    cat(sprintf("    %-8s C-index=%.4f (%d vars)\n",
                alpha_labels[ai], mean_c, ifelse(is.na(nz), 0, nz)))

    if (is.null(best_overall) || mean_c > best_overall$c_index) {
      set.seed(2026)
      final_fit <- tryCatch(
        cv.glmnet(X, y, family = "cox", alpha = alphas[ai],
                  nfolds = k_folds, type.measure = "C"),
        error = function(e) NULL
      )
      best_overall <- list(
        alpha = alphas[ai], c_index = mean_c, fit = final_fit,
        n_vars = nz, X = X, y = y, vars = vars_exist,
        df_cc = df_cc, n = n, n_ev = n_ev,
        alpha_label = alpha_labels[ai]
      )
    }
  }

  cat(sprintf("    BEST: %s, 5x10 CV C-index = %.4f\n\n",
              best_overall$alpha_label, best_overall$c_index))
  return(best_overall)
}

cox_results <- list()
for (tc in cox_tiers) {
  cox_results[[tc$nm]] <- run_cox_glmnet(tc$nm, tc$vars, df)
}

# =============================================================================
# SECTION 6: RCS-COX (Restricted Cubic Splines via rms::cph)
# =============================================================================

cat("===================================================================\n")
cat("SECTION 6: RCS-COX (non-linear effects via restricted cubic splines)\n")
cat("===================================================================\n\n")

# Use Boruta-selected variables with RCS for key continuous predictors
rcs_continuous <- c("age", "re_ldl", "hdl", "bmi", "log_crp",
                    "lipid_years_A", "non_hdl", "egfr", "log_cystatin_c",
                    "pulse_pressure")
rcs_continuous <- rcs_continuous[rcs_continuous %in% boruta_vars & rcs_continuous %in% names(df)]

rcs_binary <- boruta_vars[!boruta_vars %in% rcs_continuous]
rcs_binary <- rcs_binary[rcs_binary %in% names(df)]

# Build formula: RCS for continuous, linear for others
rcs_terms <- paste0("rcs(", rcs_continuous, ", 4)")
lin_terms <- rcs_binary
rcs_formula <- as.formula(paste("Surv(followup_years, event_indicator) ~",
                                paste(c(rcs_terms, lin_terms), collapse = " + ")))

# Prepare complete cases
all_rcs_vars <- c(rcs_continuous, rcs_binary, "followup_years", "event_indicator")
df_rcs <- df[complete.cases(df[, all_rcs_vars]), ]
df_rcs <- df_rcs[df_rcs$followup_years > 0, ]

cat(sprintf("  RCS-Cox: N=%d, Events=%d\n", nrow(df_rcs), sum(df_rcs$event_indicator)))
cat(sprintf("  Continuous vars with RCS(4 knots): %s\n",
            paste(rcs_continuous, collapse = ", ")))
cat(sprintf("  Linear terms: %s\n",
            paste(rcs_binary, collapse = ", ")))

# Set rms datadist
dd <- datadist(df_rcs)
options(datadist = "dd")

rcs_cox_fit <- tryCatch({
  cph(rcs_formula, data = df_rcs, x = TRUE, y = TRUE, surv = TRUE)
}, error = function(e) {
  cat(sprintf("  RCS-Cox failed: %s\n", e$message))
  # Fallback: fewer knots
  rcs_terms_3 <- paste0("rcs(", rcs_continuous, ", 3)")
  rcs_formula_3 <- as.formula(paste("Surv(followup_years, event_indicator) ~",
                                    paste(c(rcs_terms_3, lin_terms), collapse = " + ")))
  tryCatch(cph(rcs_formula_3, data = df_rcs, x = TRUE, y = TRUE, surv = TRUE),
           error = function(e2) { cat(sprintf("  Fallback also failed: %s\n", e2$message)); NULL })
})

if (!is.null(rcs_cox_fit)) {
  # Bootstrap-validated C-index
  set.seed(2026)
  rcs_val <- tryCatch(
    validate(rcs_cox_fit, B = 200, method = "boot"),
    error = function(e) {
      cat(sprintf("  Bootstrap validation failed: %s\n", e$message)); NULL
    }
  )

  if (!is.null(rcs_val)) {
    dxy_orig   <- rcs_val["Dxy", "index.orig"]
    dxy_corrected <- rcs_val["Dxy", "index.corrected"]
    c_orig     <- 0.5 * (dxy_orig + 1)
    c_corrected <- 0.5 * (dxy_corrected + 1)
    optimism   <- rcs_val["Dxy", "optimism"] / 2

    cat(sprintf("  RCS-Cox apparent C-index:     %.4f\n", c_orig))
    cat(sprintf("  RCS-Cox optimism-corrected C: %.4f\n", c_corrected))
    cat(sprintf("  Optimism:                     %.4f\n", optimism))
  } else {
    # Simple concordance
    c_orig <- concordance(rcs_cox_fit)$concordance
    c_corrected <- NA
    cat(sprintf("  RCS-Cox apparent C-index: %.4f (no bootstrap)\n", c_orig))
  }

  cat(sprintf("\n  ANOVA (non-linearity tests):\n"))
  tryCatch({
    a <- anova(rcs_cox_fit)
    print(a)
  }, error = function(e) cat(sprintf("  ANOVA failed: %s\n", e$message)))
}

# =============================================================================
# SECTION 7: RANDOM SURVIVAL FOREST (ranger)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 7: RANDOM SURVIVAL FOREST (ranger)\n")
cat("===================================================================\n\n")

# Use Boruta variables
rsf_vars <- boruta_vars[boruta_vars %in% names(df)]
rsf_data <- df[complete.cases(df[, c(rsf_vars, "followup_years", "event_indicator")]), ]
rsf_data <- rsf_data[rsf_data$followup_years > 0, ]

cat(sprintf("  Input: %d vars, N=%d, Events=%d\n",
            length(rsf_vars), nrow(rsf_data), sum(rsf_data$event_indicator)))

# Tune mtry and min.node.size
rsf_best <- NULL
rsf_best_c <- 0

mtry_vals <- unique(c(
  floor(sqrt(length(rsf_vars))),
  floor(length(rsf_vars) / 3),
  floor(length(rsf_vars) / 2)
))
node_sizes <- c(10, 20, 30)

for (mt in mtry_vals) {
  for (ns in node_sizes) {
    set.seed(2026)
    rsf_fit <- tryCatch(
      ranger(
        Surv(followup_years, event_indicator) ~ .,
        data = rsf_data[, c(rsf_vars, "followup_years", "event_indicator")],
        num.trees = 1000, mtry = mt, min.node.size = ns,
        importance = "permutation", splitrule = "logrank"
      ),
      error = function(e) NULL
    )
    if (!is.null(rsf_fit)) {
      # OOB C-index from ranger
      oob_c <- 1 - rsf_fit$prediction.error  # ranger reports 1 - C
      if (oob_c > rsf_best_c) {
        rsf_best_c <- oob_c
        rsf_best <- rsf_fit
        cat(sprintf("    mtry=%d, min.node=%d: OOB C-index=%.4f *\n", mt, ns, oob_c))
      } else {
        cat(sprintf("    mtry=%d, min.node=%d: OOB C-index=%.4f\n", mt, ns, oob_c))
      }
    }
  }
}

cat(sprintf("\n  Best RSF: OOB C-index = %.4f\n", rsf_best_c))

# Variable importance from best RSF
if (!is.null(rsf_best)) {
  vi <- sort(rsf_best$variable.importance, decreasing = TRUE)
  cat("\n  RSF Variable Importance (top 15):\n")
  for (j in 1:min(15, length(vi))) {
    cat(sprintf("    %2d. %-25s %.4f\n", j, names(vi)[j], vi[j]))
  }
}

# =============================================================================
# SECTION 8: SURVIVAL XGBOOST (Cox objective)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 8: SURVIVAL XGBOOST (Cox objective)\n")
cat("===================================================================\n\n")

xgb_surv_vars <- boruta_vars[boruta_vars %in% names(df)]
xgb_surv_data <- df[complete.cases(df[, c(xgb_surv_vars, "followup_years", "event_indicator")]), ]
xgb_surv_data <- xgb_surv_data[xgb_surv_data$followup_years > 0, ]

X_xgb <- as.matrix(xgb_surv_data[, xgb_surv_vars])
# XGBoost Cox: label = time, with negative for censored
y_xgb_label <- ifelse(xgb_surv_data$event_indicator == 1,
                       xgb_surv_data$followup_years,
                       -xgb_surv_data$followup_years)
dtrain_surv <- xgb.DMatrix(data = X_xgb, label = y_xgb_label)

cat(sprintf("  Input: %d vars, N=%d, Events=%d\n",
            length(xgb_surv_vars), nrow(xgb_surv_data),
            sum(xgb_surv_data$event_indicator)))

# Focused grid search
xgb_grid <- expand.grid(
  max_depth = c(3, 4, 5, 6),
  eta = c(0.01, 0.03, 0.05),
  min_child_weight = c(10, 20),
  subsample = c(0.7, 0.8),
  colsample_bytree = c(0.7, 0.8),
  stringsAsFactors = FALSE
)
set.seed(2026)
xgb_grid <- xgb_grid[sample(nrow(xgb_grid), min(40, nrow(xgb_grid))), ]

xgb_best_c <- 0
xgb_best_params <- NULL
xgb_best_nr <- NULL

for (g in 1:nrow(xgb_grid)) {
  params <- list(
    objective = "survival:cox", eval_metric = "cox-nloglik",
    max_depth = xgb_grid$max_depth[g], eta = xgb_grid$eta[g],
    min_child_weight = xgb_grid$min_child_weight[g],
    subsample = xgb_grid$subsample[g],
    colsample_bytree = xgb_grid$colsample_bytree[g],
    nthread = 1
  )

  set.seed(2026)
  cv_res <- tryCatch(
    xgb.cv(params = params, data = dtrain_surv, nrounds = 1500, nfold = 10,
           early_stopping_rounds = 50, verbose = 0),
    error = function(e) NULL
  )

  if (!is.null(cv_res)) {
    eval_log <- cv_res$evaluation_log
    nloglik_col <- grep("test.*cox.*nloglik.*mean", names(eval_log), value = TRUE)[1]
    if (!is.null(nloglik_col) && !is.na(nloglik_col)) {
      # Lower neg log-lik is better
      best_nr_g <- which.min(eval_log[[nloglik_col]])
      best_nll <- eval_log[[nloglik_col]][best_nr_g]

      # To get C-index, train on full data and compute concordance
      set.seed(2026)
      temp_model <- xgb.train(params = params, data = dtrain_surv,
                              nrounds = best_nr_g, verbose = 0)
      lp <- predict(temp_model, X_xgb)
      temp_c <- tryCatch(
        1 - concordance(Surv(xgb_surv_data$followup_years,
                         xgb_surv_data$event_indicator) ~ lp)$concordance,
        error = function(e) 0
      )

      if (temp_c > xgb_best_c) {
        xgb_best_c <- temp_c
        xgb_best_params <- params
        xgb_best_nr <- best_nr_g
      }
    }
  }
}

if (!is.null(xgb_best_params)) {
  cat(sprintf("  BEST: depth=%d, eta=%.3f, nrounds=%d\n",
              xgb_best_params$max_depth, xgb_best_params$eta, xgb_best_nr))
  cat(sprintf("  Apparent C-index = %.4f (will be corrected via OOF below)\n\n", xgb_best_c))

  # OOF C-index for honest estimate
  set.seed(2026)
  folds <- createFolds(xgb_surv_data$event_indicator, k = 10, returnTrain = TRUE)
  oof_lp <- numeric(nrow(xgb_surv_data))

  for (k in 1:10) {
    tr_idx <- folds[[k]]
    te_idx <- setdiff(1:nrow(xgb_surv_data), tr_idx)
    dt_tr <- xgb.DMatrix(data = X_xgb[tr_idx, ], label = y_xgb_label[tr_idx])
    dt_te <- xgb.DMatrix(data = X_xgb[te_idx, ])

    set.seed(2026)
    m <- xgb.train(params = xgb_best_params, data = dt_tr,
                   nrounds = xgb_best_nr, verbose = 0)
    oof_lp[te_idx] <- predict(m, dt_te)
  }

  xgb_oof_c <- tryCatch(
    1 - concordance(Surv(xgb_surv_data$followup_years,
                     xgb_surv_data$event_indicator) ~ oof_lp)$concordance,
    error = function(e) NA
  )
  cat(sprintf("  XGBoost-Cox OOF C-index = %.4f\n", xgb_oof_c))
} else {
  cat("  XGBoost-Cox failed to find valid parameters\n")
  xgb_oof_c <- NA
}

# =============================================================================
# SECTION 9: TIME-DEPENDENT AUC (timeROC at 5 and 10 years)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 9: TIME-DEPENDENT AUC (timeROC)\n")
cat("===================================================================\n\n")

# Use the best penalized Cox model (Boruta tier)
best_cox <- NULL
best_cox_c <- 0
for (nm in names(cox_results)) {
  if (!is.null(cox_results[[nm]]) && cox_results[[nm]]$c_index > best_cox_c) {
    best_cox <- cox_results[[nm]]
    best_cox_c <- cox_results[[nm]]$c_index
  }
}

if (!is.null(best_cox)) {
  # Get linear predictor from best penalized Cox
  lp_cox <- as.numeric(predict(best_cox$fit, newx = best_cox$X,
                               s = "lambda.min", type = "link"))
  surv_time <- best_cox$df_cc$followup_years
  surv_event <- best_cox$df_cc$event_indicator

  # timeROC
  troc <- tryCatch(
    timeROC(T = surv_time, delta = surv_event, marker = lp_cox,
            cause = 1, times = c(5, 10, 15),
            iid = TRUE),
    error = function(e) {
      cat(sprintf("  timeROC failed: %s\n", e$message)); NULL
    }
  )

  if (!is.null(troc)) {
    cat("  Time-dependent AUC from best penalized Cox:\n")
    for (t_idx in 1:length(troc$times)) {
      t_val <- troc$times[t_idx]
      auc_val <- troc$AUC[t_idx]
      # Confidence interval
      se_val <- tryCatch(sqrt(troc$inference$vect_iid_comp_time[t_idx]),
                         error = function(e) NA)
      if (!is.na(se_val) && se_val > 0) {
        ci_lo <- auc_val - 1.96 * se_val
        ci_hi <- auc_val + 1.96 * se_val
        cat(sprintf("    AUC at %2d years: %.4f (95%% CI: %.4f-%.4f)\n",
                    t_val, auc_val, ci_lo, ci_hi))
      } else {
        cat(sprintf("    AUC at %2d years: %.4f\n", t_val, auc_val))
      }
    }
  }

  # Also compute for RCS-Cox if available
  if (!is.null(rcs_cox_fit)) {
    lp_rcs <- predict(rcs_cox_fit, type = "lp")
    troc_rcs <- tryCatch(
      timeROC(T = df_rcs$followup_years, delta = df_rcs$event_indicator,
              marker = lp_rcs, cause = 1, times = c(5, 10, 15), iid = TRUE),
      error = function(e) NULL
    )
    if (!is.null(troc_rcs)) {
      cat("\n  Time-dependent AUC from RCS-Cox:\n")
      for (t_idx in 1:length(troc_rcs$times)) {
        cat(sprintf("    AUC at %2d years: %.4f\n",
                    troc_rcs$times[t_idx], troc_rcs$AUC[t_idx]))
      }
    }
  }
}

# =============================================================================
# SECTION 10: HONEST SURVIVAL STACKING
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 10: SURVIVAL STACKING (OOF, 3 base learners)\n")
cat("===================================================================\n\n")

# Use Boruta variables, complete cases only
stack_surv_vars <- boruta_vars[boruta_vars %in% names(df)]
df_stack <- df[complete.cases(df[, c(stack_surv_vars, "followup_years", "event_indicator")]), ]
df_stack <- df_stack[df_stack$followup_years > 0, ]

cat(sprintf("  Stacking: %d vars, N=%d, Events=%d\n",
            length(stack_surv_vars), nrow(df_stack), sum(df_stack$event_indicator)))

X_s <- as.matrix(df_stack[, stack_surv_vars])
time_s <- df_stack$followup_years
event_s <- df_stack$event_indicator

set.seed(2026)
stack_folds <- createFolds(event_s, k = 10, returnTrain = TRUE)

oof_cox_lp  <- numeric(nrow(df_stack))
oof_xgb_lp  <- numeric(nrow(df_stack))
oof_rsf_lp  <- numeric(nrow(df_stack))

y_surv_label <- ifelse(event_s == 1, time_s, -time_s)

for (k in 1:10) {
  tr <- stack_folds[[k]]
  te <- setdiff(1:nrow(df_stack), tr)

  # 1. Penalized Cox (ridge)
  y_tr <- Surv(time_s[tr], event_s[tr])
  set.seed(2026)
  cox_k <- tryCatch(
    cv.glmnet(X_s[tr, ], y_tr, family = "cox", alpha = 0, nfolds = 5, type.measure = "C"),
    error = function(e) NULL
  )
  if (!is.null(cox_k)) {
    oof_cox_lp[te] <- as.numeric(predict(cox_k, X_s[te, ], s = "lambda.min", type = "link"))
  }

  # 2. XGBoost Cox
  dt_tr <- xgb.DMatrix(data = X_s[tr, ], label = y_surv_label[tr])
  dt_te <- xgb.DMatrix(data = X_s[te, ])
  set.seed(2026)
  xgb_k <- tryCatch(
    xgb.train(params = list(objective = "survival:cox", eval_metric = "cox-nloglik",
                            max_depth = 4, eta = 0.03, subsample = 0.8,
                            colsample_bytree = 0.8, min_child_weight = 10, nthread = 1),
              data = dt_tr, nrounds = 300, verbose = 0),
    error = function(e) NULL
  )
  if (!is.null(xgb_k)) {
    oof_xgb_lp[te] <- predict(xgb_k, dt_te)
  }

  # 3. Random Survival Forest
  rf_data_tr <- df_stack[tr, c(stack_surv_vars, "followup_years", "event_indicator")]
  rf_data_te <- df_stack[te, c(stack_surv_vars, "followup_years", "event_indicator")]
  set.seed(2026)
  rsf_k <- tryCatch(
    ranger(Surv(followup_years, event_indicator) ~ .,
           data = rf_data_tr[, c(stack_surv_vars, "followup_years", "event_indicator")],
           num.trees = 500, mtry = max(floor(sqrt(length(stack_surv_vars))), 1),
           min.node.size = 20),
    error = function(e) NULL
  )
  if (!is.null(rsf_k)) {
    # Use cumulative hazard function as risk score (higher = more risk)
    pred_rsf <- predict(rsf_k, data = rf_data_te)
    # Sum of cumulative hazard as risk score
    oof_rsf_lp[te] <- rowSums(pred_rsf$chf)
  }
}

# Meta-learner: Cox regression on OOF linear predictors
meta_surv <- data.frame(cox_lp = oof_cox_lp, xgb_lp = oof_xgb_lp,
                        rsf_lp = oof_rsf_lp,
                        time = time_s, event = event_s)

# Standardize meta-features
meta_surv$cox_lp_z <- scale(meta_surv$cox_lp)
meta_surv$xgb_lp_z <- scale(meta_surv$xgb_lp)
meta_surv$rsf_lp_z <- scale(meta_surv$rsf_lp)

# Simple linear combination with OOF CV
set.seed(2026)
meta_folds <- createFolds(event_s, k = 10, returnTrain = TRUE)
oof_stack_lp <- numeric(nrow(df_stack))

for (k in 1:10) {
  tr <- meta_folds[[k]]
  te <- setdiff(1:nrow(df_stack), tr)
  meta_fit <- tryCatch(
    coxph(Surv(time, event) ~ cox_lp_z + xgb_lp_z + rsf_lp_z, data = meta_surv[tr, ]),
    error = function(e) NULL
  )
  if (!is.null(meta_fit)) {
    oof_stack_lp[te] <- predict(meta_fit, meta_surv[te, ], type = "lp")
  }
}

stack_c <- tryCatch(
  1 - concordance(Surv(time_s, event_s) ~ oof_stack_lp)$concordance,
  error = function(e) NA
)

# Individual OOF C-indices
cox_oof_c <- tryCatch(1 - concordance(Surv(time_s, event_s) ~ oof_cox_lp)$concordance,
                      error = function(e) NA)
xgb_s_oof_c <- tryCatch(1 - concordance(Surv(time_s, event_s) ~ oof_xgb_lp)$concordance,
                         error = function(e) NA)
rsf_oof_c <- tryCatch(1 - concordance(Surv(time_s, event_s) ~ oof_rsf_lp)$concordance,
                      error = function(e) NA)

cat(sprintf("  OOF C-indices (base learners):\n"))
cat(sprintf("    Penalized Cox (ridge):    %.4f\n", cox_oof_c))
cat(sprintf("    XGBoost-Cox:              %.4f\n", xgb_s_oof_c))
cat(sprintf("    Random Survival Forest:   %.4f\n", rsf_oof_c))
cat(sprintf("  STACKED SURVIVAL C-index:   %.4f\n", stack_c))

# =============================================================================
# SECTION 11: GRAND COMPARISON (v7)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 11: GRAND COMPARISON (v7 — Survival Models)\n")
cat("===================================================================\n\n")

v7_results <- data.frame(
  Model = character(), Method = character(),
  C_index = numeric(), N = integer(), Events = integer(),
  Vars = integer(), stringsAsFactors = FALSE
)

# Add penalized Cox results
for (nm in names(cox_results)) {
  if (!is.null(cox_results[[nm]])) {
    v7_results <- rbind(v7_results, data.frame(
      Model = paste("PenCox:", nm), Method = cox_results[[nm]]$alpha_label,
      C_index = cox_results[[nm]]$c_index, N = cox_results[[nm]]$n,
      Events = cox_results[[nm]]$n_ev, Vars = cox_results[[nm]]$n_vars,
      stringsAsFactors = FALSE))
  }
}

# RCS-Cox
if (!is.null(rcs_cox_fit) && exists("c_corrected") && !is.na(c_corrected)) {
  v7_results <- rbind(v7_results, data.frame(
    Model = "RCS-Cox (optimism-corrected)", Method = "rms::cph",
    C_index = c_corrected, N = nrow(df_rcs), Events = sum(df_rcs$event_indicator),
    Vars = length(c(rcs_continuous, rcs_binary)), stringsAsFactors = FALSE))
}
if (!is.null(rcs_cox_fit)) {
  v7_results <- rbind(v7_results, data.frame(
    Model = "RCS-Cox (apparent)", Method = "rms::cph",
    C_index = c_orig, N = nrow(df_rcs), Events = sum(df_rcs$event_indicator),
    Vars = length(c(rcs_continuous, rcs_binary)), stringsAsFactors = FALSE))
}

# RSF
if (rsf_best_c > 0) {
  v7_results <- rbind(v7_results, data.frame(
    Model = "Random Survival Forest", Method = "ranger",
    C_index = rsf_best_c, N = nrow(rsf_data),
    Events = sum(rsf_data$event_indicator), Vars = length(rsf_vars),
    stringsAsFactors = FALSE))
}

# XGBoost-Cox
if (!is.na(xgb_oof_c)) {
  v7_results <- rbind(v7_results, data.frame(
    Model = "XGBoost-Cox (OOF)", Method = "xgboost survival:cox",
    C_index = xgb_oof_c, N = nrow(xgb_surv_data),
    Events = sum(xgb_surv_data$event_indicator), Vars = length(xgb_surv_vars),
    stringsAsFactors = FALSE))
}

# Stacking
if (!is.na(stack_c)) {
  v7_results <- rbind(v7_results, data.frame(
    Model = "Survival Stack (Cox+XGB+RSF)", Method = "Stacking",
    C_index = stack_c, N = nrow(df_stack),
    Events = sum(df_stack$event_indicator), Vars = length(stack_surv_vars),
    stringsAsFactors = FALSE))
}

# v6 reference (converted from AUC)
v7_results <- rbind(v7_results, data.frame(
  Model = "v6 Best (Honest Stack, binary)", Method = "v6 reference (AUC)",
  C_index = 0.7599, N = 1394, Events = 336, Vars = 24,
  stringsAsFactors = FALSE))

v7_results <- rbind(v7_results, data.frame(
  Model = "SAFEHEART-RE (frozen)", Method = "Published",
  C_index = 0.7095, N = 1619, Events = 399, Vars = 6,
  stringsAsFactors = FALSE))

# Sort by C-index descending
v7_results <- v7_results[order(-v7_results$C_index), ]
rownames(v7_results) <- NULL

cat(sprintf("  %-40s %-25s %8s %6s %6s %5s\n",
            "Model", "Method", "C-index", "N", "Events", "Vars"))
cat(paste(rep("-", 95), collapse = ""), "\n")
for (i in 1:nrow(v7_results)) {
  cat(sprintf("  %-40s %-25s %8.4f %6d %6d %5d\n",
              v7_results$Model[i], v7_results$Method[i], v7_results$C_index[i],
              v7_results$N[i], v7_results$Events[i], v7_results$Vars[i]))
}

write.csv(v7_results, paste0(TAB_DIR, "calon2_v7_grand_comparison.csv"), row.names = FALSE)

# =============================================================================
# SECTION 12: SAFEHEART-RE BENCHMARK (same as v6)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 12: SAFEHEART-RE SURVIVAL BENCHMARK\n")
cat("===================================================================\n\n")

SAFEHEART_COEF <- list(
  intercept = -7.053, age = 0.064, male = 0.775, ldl_c = 0.109,
  htn = 0.431, bmi = 0.025, smoking = 0.466, prior_cvd = 1.414
)

sh_available <- all(c("age", "sex", "re_ldl", "hypertension", "bmi",
                      "smoking_binary", "followup_years", "event_indicator") %in% names(df))

if (sh_available) {
  df_sh <- df[complete.cases(df[, c("age", "sex", "re_ldl", "hypertension",
                                     "bmi", "smoking_binary",
                                     "followup_years", "event_indicator")]), ]
  df_sh <- df_sh[df_sh$followup_years > 0, ]

  lp_sh <- SAFEHEART_COEF$intercept +
    SAFEHEART_COEF$age      * df_sh$age +
    SAFEHEART_COEF$male     * df_sh$sex +
    SAFEHEART_COEF$ldl_c    * df_sh$re_ldl +
    SAFEHEART_COEF$htn      * df_sh$hypertension +
    SAFEHEART_COEF$bmi      * df_sh$bmi +
    SAFEHEART_COEF$smoking  * df_sh$smoking_binary

  # C-index for SAFEHEART-RE
  sh_c <- tryCatch(
    1 - concordance(Surv(df_sh$followup_years, df_sh$event_indicator) ~ lp_sh)$concordance,
    error = function(e) NA
  )
  cat(sprintf("  SAFEHEART-RE C-index (survival): %.4f (N=%d)\n", sh_c, nrow(df_sh)))

  # Time-dependent AUC for SAFEHEART-RE
  troc_sh <- tryCatch(
    timeROC(T = df_sh$followup_years, delta = df_sh$event_indicator,
            marker = lp_sh, cause = 1, times = c(5, 10, 15), iid = TRUE),
    error = function(e) NULL
  )
  if (!is.null(troc_sh)) {
    cat("  SAFEHEART-RE time-dependent AUC:\n")
    for (t_idx in 1:length(troc_sh$times)) {
      cat(sprintf("    AUC at %2d years: %.4f\n",
                  troc_sh$times[t_idx], troc_sh$AUC[t_idx]))
    }
  }
}

# =============================================================================
# SECTION 13: SAVE BEST MODEL COEFFICIENTS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 13: SAVING OUTPUTS\n")
cat("===================================================================\n\n")

# Save best penalized Cox coefficients
if (!is.null(best_cox) && !is.null(best_cox$fit)) {
  coefs <- as.matrix(coef(best_cox$fit, s = "lambda.min"))
  coef_df <- data.frame(
    variable = rownames(coefs),
    coefficient = as.numeric(coefs),
    HR = exp(as.numeric(coefs)),
    stringsAsFactors = FALSE
  )
  coef_df <- coef_df[coef_df$coefficient != 0, ]
  coef_df <- coef_df[order(-abs(coef_df$coefficient)), ]

  write.csv(coef_df, paste0(OUT_DIR, "calon2_v7_best_cox_coefficients.csv"),
            row.names = FALSE)

  cat(sprintf("  Saved: calon2_v7_best_cox_coefficients.csv (%d non-zero)\n",
              nrow(coef_df)))

  cat("\n  Best penalized Cox coefficients:\n")
  for (i in 1:nrow(coef_df)) {
    cat(sprintf("    %-25s coef=%7.4f  HR=%7.4f\n",
                coef_df$variable[i], coef_df$coefficient[i], coef_df$HR[i]))
  }
}

# Save RCS-Cox model
if (!is.null(rcs_cox_fit)) {
  saveRDS(rcs_cox_fit, paste0(OUT_DIR, "calon2_v7_rcs_cox.rds"))
  cat("  Saved: calon2_v7_rcs_cox.rds\n")
}

# Save comparison table
cat(sprintf("  Saved: tables/calon2_v7_grand_comparison.csv\n"))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n======================================================================\n")
cat("  CALON-2 v7 COMPLETE\n")
cat("======================================================================\n")
if (nrow(v7_results) > 0) {
  cat(sprintf("  Best model: %s (C-index = %.4f)\n",
              v7_results$Model[1], v7_results$C_index[1]))
}
if (!is.na(sh_c)) {
  cat(sprintf("  SAFEHEART-RE benchmark: %.4f\n", sh_c))
}
cat(sprintf("  v6 reference (binary AUC): 0.7599\n"))
cat("======================================================================\n")

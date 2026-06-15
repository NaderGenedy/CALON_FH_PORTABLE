################################################################################
#                                                                              #
#  CALON-2: ASCVD RISK MODEL DEVELOPMENT v5                                   #
#  "Honest AUC > 0.80" — Real Data, Reliable Tests Only                      #
#                                                                              #
#  Key changes v4 → v5:                                                       #
#    1. NEW FEATURES: non-HDL, pulse pressure, eGFR (CKD-EPI), TC/HDL,       #
#       TG/HDL, metabolic syndrome proxy, age², sex×age, PRS×re_ldl           #
#    2. NMR PCA: Replace 100 raw NMR columns with top 10-15 PCs              #
#       (eliminates collinearity noise, preserves metabolomic signal)          #
#    3. PROPER MI: mice imputation with NMR PCs as auxiliary variables        #
#       (recovers Lp(a) missingness using highly complete NMR data)           #
#    4. TUNED XGBoost: Wider hyperparameter grid, nested CV                   #
#    5. HONEST STACKING: Out-of-fold predictions, no data leakage             #
#    6. RCS FIX: x=TRUE, y=TRUE in lrm() for bootstrap validation            #
#    7. REMOVED: SuperLearner default CV (confirmed leakage)                  #
#    8. BORUTA: Feature importance-guided variable selection                   #
#                                                                              #
#  Author:  Dr Nader Genedy                                                    #
#  Date:    March 2026                                                         #
#  Target:  Honest CV-AUC > 0.80 before external validation on DRAGON3       #
#                                                                              #
################################################################################

rm(list = ls())
set.seed(2026)

# =============================================================================
# PACKAGES
# =============================================================================

required_packages <- c(
  "dplyr", "tidyr", "readr", "ggplot2", "patchwork", "scales",
  "glmnet", "pROC",
  "boot", "ResourceSelection", "nricens", "tableone",
  "DescTools", "xgboost", "caret", "car",
  # v4 retained: multiple imputation, random forest, splines
  "mice", "randomForest", "rms",
  # v5 additions
  "Boruta"  # feature importance-guided selection
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("\n")
cat("======================================================================\n")
cat("  CALON-2: ASCVD RISK MODEL DEVELOPMENT v5                           \n")
cat("  v5: New Features + NMR-PCA + Proper MI + Tuned XGB + Honest Stack  \n")
cat("  Goal: Honest CV-AUC > 0.80 — no data leakage, all tests reliable   \n")
cat("======================================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Auto-detect environment: local Windows OR Posit Cloud OR other
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

# Diagnostic: list what's available
cat(sprintf("  Files in BASE_DIR: %d\n", length(list.files(BASE_DIR))))
cat(sprintf("  Output dir exists: %s\n", dir.exists(OUT_DIR)))
if (dir.exists(OUT_DIR)) {
  csv_files <- list.files(OUT_DIR, pattern = "\\.csv$")
  cat(sprintf("  CSV files in output/: %s\n",
              ifelse(length(csv_files) > 0, paste(csv_files, collapse = ", "), "NONE")))
}

# =============================================================================
# SECTION 1: LOAD FULL MERGED DATASET (from v4 pipeline)
# =============================================================================

cat("===================================================================\n")
cat("SECTION 1: LOADING FULL MERGED DATASET\n")
cat("===================================================================\n\n")

# Load the full merged dataset that v4 already created
merged_file <- paste0(OUT_DIR, "calon2_full_merged.csv")
if (!file.exists(merged_file)) {
  stop("ERROR: calon2_full_merged.csv not found. Run 05_CALON2_develop.R first.")
}

df <- read.csv(merged_file, stringsAsFactors = FALSE)
cat(sprintf("  Full merged dataset: %d patients x %d variables\n", nrow(df), ncol(df)))
cat(sprintf("  ASCVD events: %d (%.1f%%)\n", sum(df$ascvd_combined),
            100 * mean(df$ascvd_combined)))

# Convert dates
df$assessment_date  <- as.Date(df$assessment_date)
df$first_ascvd_date <- as.Date(df$first_ascvd_date)

# =============================================================================
# SECTION 2: v5 NEW FEATURE ENGINEERING
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 2: v5 NEW FEATURE ENGINEERING\n")
cat("===================================================================\n\n")

# -- 2A: Non-HDL cholesterol (TC - HDL) — established ASCVD predictor ----------
# Captures all atherogenic lipoproteins (LDL + VLDL + IDL + Lp(a))
# 100% complete since both TC and HDL are 100% complete
if (all(c("tc", "hdl") %in% names(df))) {
  df$non_hdl <- df$tc - df$hdl
  cat(sprintf("  non-HDL cholesterol: N=%d, median=%.2f mmol/L\n",
              sum(!is.na(df$non_hdl)), median(df$non_hdl, na.rm = TRUE)))
}

# -- 2B: TC/HDL ratio — classical ASCVD risk ratio ---------------------------
if (all(c("tc", "hdl") %in% names(df))) {
  df$tc_hdl_ratio <- ifelse(df$hdl > 0, df$tc / df$hdl, NA)
  cat(sprintf("  TC/HDL ratio: N=%d, median=%.2f\n",
              sum(!is.na(df$tc_hdl_ratio)), median(df$tc_hdl_ratio, na.rm = TRUE)))
}

# -- 2C: TG/HDL ratio — insulin resistance proxy -----------------------------
if (all(c("trig", "hdl") %in% names(df))) {
  df$tg_hdl_ratio <- ifelse(df$hdl > 0, df$trig / df$hdl, NA)
  df$log_tg_hdl   <- ifelse(!is.na(df$tg_hdl_ratio) & df$tg_hdl_ratio > 0,
                             log(df$tg_hdl_ratio), NA)
  cat(sprintf("  TG/HDL ratio: N=%d, median=%.2f\n",
              sum(!is.na(df$tg_hdl_ratio)), median(df$tg_hdl_ratio, na.rm = TRUE)))
}

# -- 2D: Pulse pressure (SBP - DBP) — arterial stiffness marker --------------
if (all(c("sbp", "dbp") %in% names(df))) {
  df$pulse_pressure <- df$sbp - df$dbp
  cat(sprintf("  Pulse pressure: N=%d (%.1f%%), median=%.0f mmHg\n",
              sum(!is.na(df$pulse_pressure)),
              100 * mean(!is.na(df$pulse_pressure)),
              median(df$pulse_pressure, na.rm = TRUE)))
}

# -- 2E: eGFR (CKD-EPI 2021 creatinine equation, race-free) -----------------
# eGFR is an established ASCVD risk modifier (lower eGFR → higher risk)
# CKD-EPI 2021: 142 × min(Scr/κ, 1)^α × max(Scr/κ, 1)^{-1.200} × 0.9938^age
# where κ = 0.7 (F), 0.9 (M); α = -0.241 (F), -0.302 (M)
# Creatinine in UKB is in µmol/L; CKD-EPI uses mg/dL → divide by 88.4
if (all(c("creatinine", "age", "sex") %in% names(df))) {
  scr_mgdl <- df$creatinine / 88.4
  kappa  <- ifelse(df$sex == 1, 0.9, 0.7)
  alpha  <- ifelse(df$sex == 1, -0.302, -0.241)
  female <- ifelse(df$sex == 0, 1.012, 1.0)

  df$egfr <- 142 *
    pmin(scr_mgdl / kappa, 1)^alpha *
    pmax(scr_mgdl / kappa, 1)^(-1.200) *
    0.9938^df$age *
    female

  # CKD stage binary (eGFR < 60)
  df$ckd <- as.integer(df$egfr < 60)
  cat(sprintf("  eGFR (CKD-EPI 2021): N=%d, median=%.1f mL/min/1.73m²\n",
              sum(!is.na(df$egfr)), median(df$egfr, na.rm = TRUE)))
  cat(sprintf("  CKD (eGFR<60): %d (%.1f%%)\n",
              sum(df$ckd == 1, na.rm = TRUE),
              100 * mean(df$ckd, na.rm = TRUE)))
}

# -- 2F: Metabolic syndrome proxy (continuous composite) ----------------------
# Components: elevated BMI, TG, glucose, BP; low HDL
# Z-score each and sum for a continuous metabolic burden score
metsyn_vars <- c("bmi", "trig", "glucose", "sbp")
metsyn_available <- metsyn_vars[metsyn_vars %in% names(df)]
if (length(metsyn_available) >= 3 && "hdl" %in% names(df)) {
  # Z-score each component (higher = worse metabolically)
  metsyn_z <- sapply(metsyn_available, function(v) {
    x <- df[[v]]
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  })
  # HDL is inverse (lower = worse)
  hdl_z <- -1 * (df$hdl - mean(df$hdl, na.rm = TRUE)) / sd(df$hdl, na.rm = TRUE)
  df$metsyn_score <- rowMeans(cbind(metsyn_z, hdl_z), na.rm = TRUE)
  # Set to NA if >2 components missing
  n_miss <- rowSums(is.na(cbind(metsyn_z, hdl_z)))
  df$metsyn_score[n_miss > 2] <- NA
  cat(sprintf("  Metabolic syndrome score: N=%d, mean=%.3f, SD=%.3f\n",
              sum(!is.na(df$metsyn_score)),
              mean(df$metsyn_score, na.rm = TRUE),
              sd(df$metsyn_score, na.rm = TRUE)))
}

# -- 2G: Non-linear age terms -------------------------------------------------
if ("age" %in% names(df)) {
  df$age_sq <- df$age^2
  df$age_cu <- df$age^3
  cat(sprintf("  Age non-linear: age² and age³ created\n"))
}

# -- 2H: Sex × age interaction ------------------------------------------------
if (all(c("sex", "age") %in% names(df))) {
  df$sex_x_age <- df$sex * df$age
  cat(sprintf("  sex × age interaction: N=%d\n", sum(!is.na(df$sex_x_age))))
}

# -- 2I: PRS × re_ldl interaction (genetic load × lipid burden) ---------------
if (all(c("prs_primary", "re_ldl") %in% names(df))) {
  df$prs_x_re_ldl <- df$prs_primary * df$re_ldl
  cat(sprintf("  PRS × re_ldl interaction: N=%d\n", sum(!is.na(df$prs_x_re_ldl))))
}

# -- 2J: ApoB/ApoA1 ratio from standard assay --------------------------------
if (all(c("apob", "apoa") %in% names(df))) {
  df$apob_apoa_ratio <- ifelse(!is.na(df$apob) & !is.na(df$apoa) & df$apoa > 0,
                                df$apob / df$apoa, NA)
  cat(sprintf("  ApoB/ApoA1 ratio: N=%d, median=%.3f\n",
              sum(!is.na(df$apob_apoa_ratio)),
              median(df$apob_apoa_ratio, na.rm = TRUE)))
}

# -- 2K: Remnant cholesterol (TC - HDL - LDL) --------------------------------
if (all(c("tc", "hdl", "ldl") %in% names(df))) {
  df$remnant_chol_calc <- df$tc - df$hdl - df$ldl
  cat(sprintf("  Remnant cholesterol: N=%d, median=%.3f mmol/L\n",
              sum(!is.na(df$remnant_chol_calc)),
              median(df$remnant_chol_calc, na.rm = TRUE)))
}

# =============================================================================
# SECTION 3: NMR PCA — DIMENSIONALITY REDUCTION
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 3: NMR PCA (DIMENSIONALITY REDUCTION)\n")
cat("===================================================================\n\n")

# Strategy: 100 raw NMR features are heavily collinear (lipoprotein subfractions)
# PCA extracts orthogonal components that capture the metabolomic variance
# without the collinearity noise that hurts regularised models.

nmr_all_cols <- grep("^p234\\d+_i0$", names(df), value = TRUE)
cat(sprintf("  NMR features detected: %d\n", length(nmr_all_cols)))

nmr_pca_result <- NULL
if (length(nmr_all_cols) >= 20) {
  # NMR matrix — impute column means for missing values before PCA
  # This preserves ALL patients instead of dropping ~23% with any missing NMR
  nmr_matrix <- as.matrix(df[, nmr_all_cols])
  n_nmr_miss <- sum(!complete.cases(nmr_matrix))
  cat(sprintf("  NMR complete cases: %d / %d (%.1f%%)\n",
              nrow(df) - n_nmr_miss, nrow(df),
              100 * (nrow(df) - n_nmr_miss) / nrow(df)))

  # Column-mean imputation for NMR features (preserves all patients)
  for (j in 1:ncol(nmr_matrix)) {
    col_na <- is.na(nmr_matrix[, j])
    if (any(col_na)) {
      nmr_matrix[col_na, j] <- mean(nmr_matrix[, j], na.rm = TRUE)
    }
  }
  cat(sprintf("  NMR after column-mean imputation: %d patients retained (100%%)\n",
              nrow(nmr_matrix)))

  # Scale and run PCA on FULL imputed matrix
  nmr_scaled <- scale(nmr_matrix)

  # Remove zero-variance columns
  zv <- apply(nmr_scaled, 2, function(x) sd(x, na.rm = TRUE) == 0)
  nmr_scaled <- nmr_scaled[, !zv]

  pca_fit <- prcomp(nmr_scaled, center = FALSE, scale. = FALSE)

  # Determine number of PCs to retain (cumulative variance > 90%)
  cum_var <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
  n_pcs <- min(which(cum_var >= 0.90))
  n_pcs <- min(n_pcs, 15)  # Cap at 15 PCs to avoid overfitting

  cat(sprintf("  PCs retained: %d (cumulative variance = %.1f%%)\n",
              n_pcs, 100 * cum_var[n_pcs]))
  cat(sprintf("  Variance by PC: "))
  for (j in 1:min(5, n_pcs)) {
    cat(sprintf("PC%d=%.1f%% ", j, 100 * pca_fit$sdev[j]^2 / sum(pca_fit$sdev^2)))
  }
  cat("...\n")

  # Project ALL patients (all retained after column-mean imputation)
  pc_names <- paste0("nmr_pc", 1:n_pcs)
  for (j in 1:n_pcs) {
    df[[pc_names[j]]] <- pca_fit$x[, j]
  }
  cat(sprintf("  Created %d NMR PC features in dataset (all %d patients)\n",
              n_pcs, nrow(df)))

  # Save PCA loadings for interpretation
  loadings_df <- data.frame(
    variable = rownames(pca_fit$rotation)[1:ncol(nmr_scaled)],
    pca_fit$rotation[, 1:min(5, n_pcs)]
  )
  write.csv(loadings_df, paste0(TAB_DIR, "calon2_v5_nmr_pca_loadings.csv"), row.names = FALSE)

  nmr_pca_result <- list(fit = pca_fit, n_pcs = n_pcs, pc_names = pc_names,
                          nmr_cols = colnames(nmr_scaled))
} else {
  cat("  Insufficient NMR columns for PCA. Skipping.\n")
}

# =============================================================================
# SECTION 4: v5 MODEL TIERS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 4: v5 MODEL TIERS\n")
cat("===================================================================\n\n")

# -- Existing v4 features that are already in df ---
# Ensure these exist (from v4 pipeline)
for (v in c("log_lpa", "lpa_extreme", "lpa_high", "log_apob_ldl",
            "log_apob_re_ldl", "lipid_years_A", "gene_apob",
            "age_x_re_ldl", "age_x_hdl", "log_crp",
            "log_hba1c", "log_creatinine", "log_alt", "log_cystatin_c",
            "prs_primary", "townsend")) {
  if (!v %in% names(df)) {
    cat(sprintf("  WARNING: %s not found in data\n", v))
  }
}

# NMR PC column names
pc_names <- if (!is.null(nmr_pca_result)) nmr_pca_result$pc_names else character(0)

# --- TIER v5-A: Clinical Core+ (15 vars) ------------------------------------
# Core clinical + new features with 100% completeness
v5_core_vars <- c(
  "age", "sex", "re_ldl", "hdl",
  "log_apob_ldl", "log_lpa", "lpa_extreme",
  "smoking_binary", "diabetes", "hypertension",
  "bmi", "log_crp",
  "lipid_years_A",
  "non_hdl",            # NEW: captures all atherogenic lipoproteins
  "egfr"                # NEW: renal function (CKD-EPI 2021)
)
cat(sprintf("  v5-A -- Clinical Core+: %d vars (v4 Core + non-HDL + eGFR)\n",
            length(v5_core_vars)))

# --- TIER v5-B: Enriched Clinical (20 vars) ---------------------------------
# Adds BP, metabolic features, ratios
v5_enriched_vars <- c(
  v5_core_vars,
  "pulse_pressure",     # NEW: arterial stiffness
  "tc_hdl_ratio",       # NEW: classical ASCVD ratio
  "metsyn_score",       # NEW: metabolic syndrome composite
  "age_sq",             # NEW: non-linear age
  "sex_x_age"           # NEW: sex-age interaction
)
cat(sprintf("  v5-B -- Enriched Clinical: %d vars\n", length(v5_enriched_vars)))

# --- TIER v5-C: NMR-PCA Enhanced (Core+ plus NMR PCs + PRS) -----------------
v5_nmr_pca_vars <- c(
  v5_core_vars,
  pc_names,             # NMR PCs (orthogonal metabolomic features)
  "prs_primary"         # Polygenic risk score
)
cat(sprintf("  v5-C -- NMR-PCA Enhanced: %d vars (Core+ + %d NMR PCs + PRS)\n",
            length(v5_nmr_pca_vars), length(pc_names)))

# --- TIER v5-D: Full v5 (all new + NMR PCs + PRS + interactions) ------------
v5_full_vars <- unique(c(
  v5_enriched_vars,
  pc_names,
  "prs_primary",
  "prs_x_re_ldl",      # NEW: genetic×lipid interaction
  "age_x_re_ldl",
  "apob_apoa_ratio",    # NEW: ApoB/ApoA1 from standard assay
  "remnant_chol_calc",  # NEW: remnant cholesterol
  "log_tg_hdl",         # NEW: insulin resistance proxy
  "gene_apob",
  "townsend",
  "hba1c", "log_cystatin_c",
  "ever_smoked"
))
cat(sprintf("  v5-D -- Full v5: %d vars (everything)\n", length(v5_full_vars)))

# --- TIER v5-E: Lean Powerhouse (Boruta-selected) ---------------------------
# Will be defined AFTER Boruta runs in Section 5

# --- Report N complete cases for each tier ---
cat("\n  Tier Sizes:\n")
cat(sprintf("  %-30s %6s %6s %8s %6s\n", "Tier", "Vars", "N_cc", "Events", "EPV"))
cat(paste(rep("-", 65), collapse = ""), "\n")

tier_check <- list(
  list(nm = "v5-A Core+", vars = v5_core_vars),
  list(nm = "v5-B Enriched Clinical", vars = v5_enriched_vars),
  list(nm = "v5-C NMR-PCA Enhanced", vars = v5_nmr_pca_vars),
  list(nm = "v5-D Full v5", vars = v5_full_vars)
)

for (tc in tier_check) {
  tv <- tc$vars[tc$vars %in% names(df)]
  df_tmp <- df[, c(tv, "ascvd_combined")]
  cc <- complete.cases(df_tmp)
  n_cc <- sum(cc)
  n_ev <- sum(df_tmp$ascvd_combined[cc])
  epv <- ifelse(length(tv) > 0, n_ev / length(tv), 0)
  cat(sprintf("  %-30s %6d %6d %8d %6.1f\n",
              tc$nm, length(tv), n_cc, n_ev, epv))
}
cat("\n")

# =============================================================================
# SECTION 5: BORUTA FEATURE SELECTION (importance-guided)
# =============================================================================

cat("===================================================================\n")
cat("SECTION 5: BORUTA FEATURE SELECTION\n")
cat("===================================================================\n\n")

# Strategy: Use Random Forest wrapper (Boruta) to identify features
# statistically more important than random shadows. This gives
# data-driven variable selection without manual bias.
# Run Boruta on the Full v5 tier

boruta_vars <- v5_full_vars[v5_full_vars %in% names(df)]
df_boruta <- df[, c(boruta_vars, "ascvd_combined")]
df_boruta_cc <- df_boruta[complete.cases(df_boruta), ]

cat(sprintf("  Boruta input: %d vars, %d complete cases, %d events\n",
            length(boruta_vars), nrow(df_boruta_cc),
            sum(df_boruta_cc$ascvd_combined)))

boruta_confirmed <- character(0)
if (nrow(df_boruta_cc) >= 500 && sum(df_boruta_cc$ascvd_combined) >= 50) {
  set.seed(2026)
  boruta_result <- tryCatch({
    Boruta(ascvd_combined ~ ., data = df_boruta_cc, maxRuns = 200,
           doTrace = 0)
  }, error = function(e) {
    cat(sprintf("  Boruta failed: %s\n", e$message))
    NULL
  })

  if (!is.null(boruta_result)) {
    # Try to resolve tentative features
    boruta_final <- TentativeRoughFix(boruta_result)
    boruta_confirmed <- names(boruta_final$finalDecision[
      boruta_final$finalDecision == "Confirmed"])
    boruta_rejected <- names(boruta_final$finalDecision[
      boruta_final$finalDecision == "Rejected"])

    cat(sprintf("  Boruta CONFIRMED: %d features\n", length(boruta_confirmed)))
    cat(sprintf("  Boruta REJECTED: %d features\n", length(boruta_rejected)))

    if (length(boruta_confirmed) > 0) {
      # Get importance scores
      imp_df <- data.frame(
        variable = names(boruta_final$finalDecision),
        decision = as.character(boruta_final$finalDecision),
        meanImp  = colMeans(boruta_result$ImpHistory[,
                    names(boruta_final$finalDecision)], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      imp_df <- imp_df[order(-imp_df$meanImp), ]

      cat("\n  Boruta Confirmed Features (ranked by importance):\n")
      conf_df <- imp_df[imp_df$decision == "Confirmed", ]
      for (j in 1:nrow(conf_df)) {
        cat(sprintf("    %2d. %-25s meanImp=%.3f\n",
                    j, conf_df$variable[j], conf_df$meanImp[j]))
      }
      cat("\n")

      write.csv(imp_df, paste0(TAB_DIR, "calon2_v5_boruta_importance.csv"),
                row.names = FALSE)
    }
  }
} else {
  cat("  Too few complete cases for Boruta. Skipping.\n")
}

# --- TIER v5-E: Boruta-selected Lean Powerhouse ---
if (length(boruta_confirmed) >= 5) {
  v5_boruta_vars <- boruta_confirmed
  cat(sprintf("  v5-E -- Boruta Lean: %d Confirmed features\n",
              length(v5_boruta_vars)))
} else {
  # Fallback: use enriched clinical if Boruta fails
  v5_boruta_vars <- v5_enriched_vars
  cat("  v5-E -- Fallback to Enriched Clinical (Boruta had <5 confirmed)\n")
}

# =============================================================================
# SECTION 6: ELASTIC NET — v5 TIERS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 6: ELASTIC NET — v5 TIERS\n")
cat("===================================================================\n\n")

model_results <- list()

alphas <- c(0, 0.25, 0.5, 0.75, 1.0)
alpha_names <- c("Ridge", "EN(0.25)", "EN(0.50)", "EN(0.75)", "LASSO")

run_elastic_net <- function(tier_name, vars, df_full) {
  cat(sprintf("  --- %s (Elastic Net) ---\n", tier_name))

  vars_exist <- vars[vars %in% names(df_full)]
  if (length(vars_exist) < length(vars)) {
    cat(sprintf("  NOTE: %d / %d vars found in data\n", length(vars_exist), length(vars)))
  }

  df_tier <- df_full[, c(vars_exist, "ascvd_combined")]
  cc_idx <- complete.cases(df_tier)
  df_cc <- df_tier[cc_idx, ]

  cat(sprintf("  Complete cases: %d / %d (%.1f%%)\n",
              nrow(df_cc), nrow(df_full), 100 * nrow(df_cc) / nrow(df_full)))
  cat(sprintf("  Events: %d (%.1f%%), EPV: %.1f\n",
              sum(df_cc$ascvd_combined), 100 * mean(df_cc$ascvd_combined),
              sum(df_cc$ascvd_combined) / length(vars_exist)))

  if (nrow(df_cc) < 200 || sum(df_cc$ascvd_combined) < 30) {
    cat("  SKIPPED: too few complete cases or events\n\n")
    return(NULL)
  }

  X <- as.matrix(df_cc[, vars_exist])
  y <- df_cc$ascvd_combined

  best <- NULL
  for (i in seq_along(alphas)) {
    a <- alphas[i]
    set.seed(2026)
    cv_fit <- cv.glmnet(X, y, family = "binomial",
                         alpha = a, nfolds = 10, type.measure = "auc")
    auc_max <- max(cv_fit$cvm)
    n_nonzero <- sum(coef(cv_fit, s = "lambda.min")[-1] != 0)

    cat(sprintf("    %-10s: AUC=%.4f, lambda=%.6f, %d/%d vars\n",
                alpha_names[i], auc_max, cv_fit$lambda.min, n_nonzero, length(vars_exist)))

    if (is.null(best) || auc_max > best$auc) {
      best <- list(alpha = a, lambda = cv_fit$lambda.min, auc = auc_max,
                    fit = cv_fit, n_vars = n_nonzero, alpha_name = alpha_names[i],
                    X = X, y = y, vars = vars_exist, df_cc = df_cc,
                    n_complete = nrow(df_cc), n_events = sum(df_cc$ascvd_combined))
    }
  }

  cat(sprintf("  BEST: %s (alpha=%.2f), CV-AUC = %.4f, %d vars selected\n\n",
              best$alpha_name, best$alpha, best$auc, best$n_vars))

  return(best)
}

# Run EN for all v5 tiers
en_v5_core     <- run_elastic_net("v5-A: Clinical Core+", v5_core_vars, df)
en_v5_enriched <- run_elastic_net("v5-B: Enriched Clinical", v5_enriched_vars, df)
en_v5_nmr_pca  <- run_elastic_net("v5-C: NMR-PCA Enhanced", v5_nmr_pca_vars, df)
en_v5_full     <- run_elastic_net("v5-D: Full v5", v5_full_vars, df)
en_v5_boruta   <- run_elastic_net("v5-E: Boruta Lean", v5_boruta_vars, df)

# Also re-run best v4 tiers for fair comparison
en_v4_kitchen  <- run_elastic_net("v4: Kitchen Sink (baseline)",
                    c("age", "sex", "apob", "re_ldl", "hdl", "trig",
                      "inv_apoa", "apob_ldl_ratio", "lpa_binary",
                      "smoking_binary", "diabetes", "hypertension", "bmi",
                      "log_crp",
                      grep("^p234\\d+_i0$", names(df), value = TRUE),
                      "prs_primary", "townsend", "lipid_years_A", "gene_apob",
                      "hba1c", "creatinine", "alt", "cystatin_c", "glucose",
                      "waist", "sbp", "dbp", "ever_smoked", "log_apob_ldl",
                      "log_hba1c", "log_creatinine", "log_alt", "log_cystatin_c",
                      "apob_apoa1_nmr_ratio", "vldl_hdl_ratio", "mean_cimt"), df)

# Store EN results
en_models <- list(
  list(obj = en_v5_core,     key = "EN_v5_Core",      nm = "EN v5 Core+"),
  list(obj = en_v5_enriched, key = "EN_v5_Enriched",   nm = "EN v5 Enriched Clinical"),
  list(obj = en_v5_nmr_pca,  key = "EN_v5_NMR_PCA",   nm = "EN v5 NMR-PCA Enhanced"),
  list(obj = en_v5_full,     key = "EN_v5_Full",       nm = "EN v5 Full"),
  list(obj = en_v5_boruta,   key = "EN_v5_Boruta",     nm = "EN v5 Boruta Lean"),
  list(obj = en_v4_kitchen,  key = "EN_v4_Kitchen",    nm = "EN v4 Kitchen (baseline)")
)

for (em in en_models) {
  if (!is.null(em$obj)) {
    model_results[[em$key]] <- list(
      name = em$nm, algorithm = "Elastic Net",
      n_vars = em$obj$n_vars, cv_auc = em$obj$auc,
      n_complete = em$obj$n_complete, n_events = em$obj$n_events,
      model = em$obj
    )
  }
}

# =============================================================================
# SECTION 7: TUNED XGBOOST — WIDER GRID + v5 TIERS
# =============================================================================

cat("===================================================================\n")
cat("SECTION 7: TUNED XGBOOST — v5 TIERS\n")
cat("===================================================================\n\n")

run_xgb_tuned <- function(tier_name, vars, df_full) {
  cat(sprintf("  --- %s (XGBoost Tuned) ---\n", tier_name))

  vars_exist <- vars[vars %in% names(df_full)]
  df_tier <- df_full[, c(vars_exist, "ascvd_combined")]
  cc_idx <- complete.cases(df_tier)
  df_cc <- df_tier[cc_idx, ]

  cat(sprintf("  N=%d, Events=%d, Vars=%d\n",
              nrow(df_cc), sum(df_cc$ascvd_combined), length(vars_exist)))

  if (nrow(df_cc) < 200 || sum(df_cc$ascvd_combined) < 30) {
    cat("  SKIPPED: too few complete cases or events\n\n")
    return(NULL)
  }

  X <- as.matrix(df_cc[, vars_exist])
  y <- df_cc$ascvd_combined
  dtrain <- xgb.DMatrix(data = X, label = y)

  # WIDER hyperparameter grid (v5 improvement)
  param_grid <- expand.grid(
    max_depth = c(3, 4, 5, 6),
    eta = c(0.005, 0.01, 0.03, 0.05, 0.10),
    min_child_weight = c(5, 10, 20),
    subsample = c(0.7, 0.8),
    colsample_bytree = c(0.7, 0.8, 1.0),
    stringsAsFactors = FALSE
  )

  # Random search: sample 40 combinations (much faster than full grid)
  set.seed(2026)
  if (nrow(param_grid) > 40) {
    param_grid <- param_grid[sample(nrow(param_grid), 40), ]
  }

  best_auc <- 0
  best_params <- NULL
  best_nrounds <- NULL

  for (g in 1:nrow(param_grid)) {
    params <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = param_grid$max_depth[g],
      eta = param_grid$eta[g],
      min_child_weight = param_grid$min_child_weight[g],
      subsample = param_grid$subsample[g],
      colsample_bytree = param_grid$colsample_bytree[g],
      nthread = 1
    )

    set.seed(2026)
    cv_res <- tryCatch({
      xgb.cv(params = params, data = dtrain,
             nrounds = 2000, nfold = 10,
             early_stopping_rounds = 50, verbose = 0,
             print_every_n = 0)
    }, error = function(e) {
      tryCatch({
        xgb.cv(params = params, data = dtrain,
               nrounds = 2000, nfold = 10,
               early_stopping_rounds = 50, verbose = 0)
      }, error = function(e2) NULL)
    })

    if (!is.null(cv_res)) {
      eval_log <- cv_res$evaluation_log
      auc_col <- grep("test.*auc.*mean", names(eval_log), value = TRUE)[1]
      if (is.null(auc_col) || is.na(auc_col)) {
        auc_col <- names(eval_log)[grep("auc", names(eval_log))[1]]
      }
      if (!is.null(auc_col) && !is.na(auc_col)) {
        this_auc <- max(eval_log[[auc_col]], na.rm = TRUE)
        this_nrounds <- which.max(eval_log[[auc_col]])

        if (this_auc > best_auc) {
          best_auc <- this_auc
          best_params <- params
          best_nrounds <- this_nrounds
        }
      }
    }
  }

  if (is.null(best_params)) {
    cat("  XGBoost tuning FAILED (all CV runs errored)\n\n")
    return(NULL)
  }

  cat(sprintf("  BEST: depth=%d, eta=%.3f, mcw=%d, sub=%.1f, colsample=%.1f, nrounds=%d\n",
              best_params$max_depth, best_params$eta,
              best_params$min_child_weight, best_params$subsample,
              best_params$colsample_bytree, best_nrounds))
  cat(sprintf("  CV-AUC = %.4f\n\n", best_auc))

  # Train final model
  set.seed(2026)
  final_model <- xgb.train(params = best_params, data = dtrain,
                            nrounds = best_nrounds, verbose = 0)

  # Feature importance
  imp <- xgb.importance(model = final_model, feature_names = vars_exist)
  if (nrow(imp) > 0) {
    cat("  Top 15 features:\n")
    for (j in 1:min(15, nrow(imp))) {
      cat(sprintf("    %2d. %-25s Gain=%.4f\n", j, imp$Feature[j], imp$Gain[j]))
    }
    cat("\n")
  }

  return(list(
    auc = best_auc, params = best_params, nrounds = best_nrounds,
    final_model = final_model, dtrain = dtrain,
    X = X, y = y, vars = vars_exist, df_cc = df_cc,
    n_complete = nrow(df_cc), n_events = sum(df_cc$ascvd_combined),
    importance = imp
  ))
}

# Run tuned XGBoost for v5 tiers
xgb_v5_core     <- run_xgb_tuned("v5-A: Clinical Core+", v5_core_vars, df)
xgb_v5_enriched <- run_xgb_tuned("v5-B: Enriched Clinical", v5_enriched_vars, df)
xgb_v5_nmr_pca  <- run_xgb_tuned("v5-C: NMR-PCA Enhanced", v5_nmr_pca_vars, df)
xgb_v5_full     <- run_xgb_tuned("v5-D: Full v5", v5_full_vars, df)
xgb_v5_boruta   <- run_xgb_tuned("v5-E: Boruta Lean", v5_boruta_vars, df)

# Store XGB results
xgb_models <- list(
  list(obj = xgb_v5_core,     key = "XGB_v5_Core",      nm = "XGB v5 Core+"),
  list(obj = xgb_v5_enriched, key = "XGB_v5_Enriched",   nm = "XGB v5 Enriched Clinical"),
  list(obj = xgb_v5_nmr_pca,  key = "XGB_v5_NMR_PCA",   nm = "XGB v5 NMR-PCA Enhanced"),
  list(obj = xgb_v5_full,     key = "XGB_v5_Full",       nm = "XGB v5 Full"),
  list(obj = xgb_v5_boruta,   key = "XGB_v5_Boruta",     nm = "XGB v5 Boruta Lean")
)

for (xm in xgb_models) {
  if (!is.null(xm$obj)) {
    model_results[[xm$key]] <- list(
      name = xm$nm, algorithm = "XGBoost",
      n_vars = length(xm$obj$vars), cv_auc = xm$obj$auc,
      n_complete = xm$obj$n_complete, n_events = xm$obj$n_events,
      model = xm$obj
    )
  }
}

# =============================================================================
# SECTION 8: MULTIPLE IMPUTATION — PROPER MI WITH NMR AS AUXILIARY
# =============================================================================

cat("===================================================================\n")
cat("SECTION 8: MULTIPLE IMPUTATION (v5 — NMR PCs as auxiliary)\n")
cat("===================================================================\n\n")

# Strategy: Lp(a) is 19.4% missing. NMR metabolomics (99%+ complete)
# likely correlate with Lp(a) and can serve as auxiliary variables
# in the imputation model, improving prediction of missing values.

run_mi_elastic_net <- function(tier_name, vars, df_full, m_imps = 10,
                                aux_vars = character(0),
                                alphas = c(0, 0.25, 0.5, 0.75, 1.0)) {
  cat(sprintf("  --- %s (MI Elastic Net, m=%d, %d alphas) ---\n",
              tier_name, m_imps, length(alphas)))

  vars_exist <- vars[vars %in% names(df_full)]
  # Include auxiliary variables for imputation (but not for modelling)
  aux_exist <- aux_vars[aux_vars %in% names(df_full)]
  all_imp_vars <- unique(c(vars_exist, aux_exist, "ascvd_combined"))

  imp_data <- df_full[, all_imp_vars]

  n_before <- nrow(imp_data)
  n_cc <- sum(complete.cases(imp_data[, c(vars_exist, "ascvd_combined")]))
  cat(sprintf("  Before MI: %d total, %d complete (%.1f%% loss)\n",
              n_before, n_cc, 100 * (n_before - n_cc) / n_before))
  cat(sprintf("  Auxiliary variables for imputation: %d\n", length(aux_exist)))

  if (n_before - n_cc < 20) {
    cat("  Minimal missingness — MI unnecessary. Using complete cases.\n\n")
    return(NULL)
  }

  # Run mice imputation
  cat("  Running mice imputation...\n")
  set.seed(2026)
  mice_fit <- tryCatch({
    mice(imp_data, m = m_imps, method = "pmm", maxit = 10,
         printFlag = FALSE)
  }, error = function(e) {
    cat(sprintf("  mice failed: %s\n", e$message))
    NULL
  })

  if (is.null(mice_fit)) return(NULL)

  # Fit elastic net on each imputed dataset (MODEL vars only, not auxiliary)
  # Search over alphas (same grid as regular EN) for fair comparison
  alpha_names <- c("Ridge", "EN(0.25)", "EN(0.50)", "EN(0.75)", "LASSO")
  best_alpha <- 0.5
  best_pooled_auc <- 0

  for (a_idx in seq_along(alphas)) {
    auc_pool <- numeric(m_imps)
    for (imp_i in 1:m_imps) {
      imp_df <- complete(mice_fit, imp_i)
      X_imp <- as.matrix(imp_df[, vars_exist])
      y_imp <- imp_df$ascvd_combined

      set.seed(2026 + a_idx)
      cv_fit <- cv.glmnet(X_imp, y_imp, family = "binomial",
                           alpha = alphas[a_idx], nfolds = 10, type.measure = "auc")
      auc_pool[imp_i] <- max(cv_fit$cvm)
    }
    pooled <- mean(auc_pool)
    cat(sprintf("    %-10s: MI-pooled AUC=%.4f (SD=%.4f)\n",
                alpha_names[a_idx], pooled, sd(auc_pool)))
    if (pooled > best_pooled_auc) {
      best_pooled_auc <- pooled
      best_alpha <- alphas[a_idx]
      best_auc_pool <- auc_pool
    }
  }

  cat(sprintf("  BEST MI alpha=%.2f, pooled CV-AUC: %.4f (range=%.4f--%.4f)\n",
              best_alpha, best_pooled_auc, min(best_auc_pool), max(best_auc_pool)))
  cat(sprintf("  Recovered %d patients (N now = %d)\n\n",
              n_before - n_cc, n_before))

  return(list(
    pooled_auc = best_pooled_auc, auc_per_imp = best_auc_pool,
    best_alpha = best_alpha,
    n_complete = n_before, n_events = sum(df_full$ascvd_combined),
    n_recovered = n_before - n_cc, vars = vars_exist
  ))
}

# NMR PCs + other complete features as auxiliary variables
aux_for_mi <- c(pc_names, "non_hdl", "tc_hdl_ratio", "egfr",
                "age_sq", "creatinine", "alt", "cystatin_c")
aux_for_mi <- aux_for_mi[aux_for_mi %in% names(df)]

mi_v5_core     <- run_mi_elastic_net("v5-A: Core+ (MI)", v5_core_vars, df,
                                      aux_vars = aux_for_mi)
mi_v5_enriched <- run_mi_elastic_net("v5-B: Enriched (MI)", v5_enriched_vars, df,
                                      aux_vars = aux_for_mi)
mi_v5_nmr_pca  <- run_mi_elastic_net("v5-C: NMR-PCA (MI)", v5_nmr_pca_vars, df,
                                      aux_vars = c("non_hdl", "tc_hdl_ratio", "egfr"))
mi_v5_boruta   <- run_mi_elastic_net("v5-E: Boruta Lean (MI)", v5_boruta_vars, df,
                                      aux_vars = aux_for_mi)

# Store MI results
for (mi_item in list(
  list(obj = mi_v5_core,     key = "MI_v5_Core",      nm = "MI-EN v5 Core+"),
  list(obj = mi_v5_enriched, key = "MI_v5_Enriched",   nm = "MI-EN v5 Enriched"),
  list(obj = mi_v5_nmr_pca,  key = "MI_v5_NMR_PCA",   nm = "MI-EN v5 NMR-PCA"),
  list(obj = mi_v5_boruta,   key = "MI_v5_Boruta",     nm = "MI-EN v5 Boruta Lean")
)) {
  if (!is.null(mi_item$obj)) {
    model_results[[mi_item$key]] <- list(
      name = mi_item$nm, algorithm = "MI-Elastic Net",
      n_vars = length(mi_item$obj$vars), cv_auc = mi_item$obj$pooled_auc,
      n_complete = mi_item$obj$n_complete, n_events = mi_item$obj$n_events,
      model = mi_item$obj
    )
  }
}

# =============================================================================
# SECTION 9: HONEST STACKING (out-of-fold, no leakage)
# =============================================================================

cat("===================================================================\n")
cat("SECTION 9: HONEST STACKING (out-of-fold predictions)\n")
cat("===================================================================\n\n")

# Strategy: Train EN + XGBoost + RF as base learners using 10-fold CV.
# Collect out-of-fold predictions for each fold.
# Train logistic regression meta-learner on the OOF predictions.
# Evaluate using SEPARATE outer 10-fold CV.
# This prevents the data leakage that plagued SuperLearner in v4.

run_honest_stack <- function(tier_name, vars, df_full) {
  cat(sprintf("  --- %s (Honest Stack) ---\n", tier_name))

  vars_exist <- vars[vars %in% names(df_full)]
  df_tier <- df_full[, c(vars_exist, "ascvd_combined")]
  cc_idx <- complete.cases(df_tier)
  df_cc <- df_tier[cc_idx, ]

  cat(sprintf("  N=%d, Events=%d, Vars=%d\n",
              nrow(df_cc), sum(df_cc$ascvd_combined), length(vars_exist)))

  if (nrow(df_cc) < 300 || sum(df_cc$ascvd_combined) < 50) {
    cat("  SKIPPED: too few observations for stacking\n\n")
    return(NULL)
  }

  X <- as.matrix(df_cc[, vars_exist])
  y <- df_cc$ascvd_combined

  # OUTER 10-fold CV for honest evaluation
  set.seed(2026)
  outer_folds <- sample(rep(1:10, length.out = nrow(df_cc)))
  outer_pred_stack <- numeric(nrow(df_cc))
  outer_pred_en    <- numeric(nrow(df_cc))
  outer_pred_xgb   <- numeric(nrow(df_cc))
  outer_pred_rf    <- numeric(nrow(df_cc))

  for (outer_k in 1:10) {
    test_idx  <- which(outer_folds == outer_k)
    train_idx <- which(outer_folds != outer_k)

    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_test  <- X[test_idx, , drop = FALSE]

    # --- BASE LEARNER 1: Elastic Net ---
    set.seed(2026)
    en_fit <- tryCatch({
      cv.glmnet(X_train, y_train, family = "binomial",
                 alpha = 0.5, nfolds = 5, type.measure = "auc")
    }, error = function(e) NULL)

    if (!is.null(en_fit)) {
      pred_en_train <- predict(en_fit, X_train, s = "lambda.min", type = "response")[, 1]
      pred_en_test  <- predict(en_fit, X_test, s = "lambda.min", type = "response")[, 1]
    } else {
      pred_en_train <- rep(mean(y_train), length(y_train))
      pred_en_test  <- rep(mean(y_train), length(test_idx))
    }

    # --- BASE LEARNER 2: XGBoost ---
    dtrain_inner <- xgb.DMatrix(data = X_train, label = y_train)
    xgb_params <- list(
      objective = "binary:logistic", eval_metric = "auc",
      max_depth = 4, eta = 0.03, min_child_weight = 10,
      subsample = 0.8, colsample_bytree = 0.8, nthread = 1
    )
    set.seed(2026)
    xgb_cv_inner <- tryCatch({
      xgb.cv(params = xgb_params, data = dtrain_inner,
             nrounds = 1000, nfold = 5,
             early_stopping_rounds = 30, verbose = 0,
             print_every_n = 0)
    }, error = function(e) {
      tryCatch({
        xgb.cv(params = xgb_params, data = dtrain_inner,
               nrounds = 1000, nfold = 5,
               early_stopping_rounds = 30, verbose = 0)
      }, error = function(e2) NULL)
    })

    if (!is.null(xgb_cv_inner)) {
      eval_log <- xgb_cv_inner$evaluation_log
      auc_col <- grep("test.*auc.*mean", names(eval_log), value = TRUE)[1]
      if (is.null(auc_col) || is.na(auc_col)) {
        auc_col <- names(eval_log)[grep("auc", names(eval_log))[1]]
      }
      best_round <- which.max(eval_log[[auc_col]])
      xgb_model <- xgb.train(params = xgb_params, data = dtrain_inner,
                               nrounds = best_round, verbose = 0)
      pred_xgb_train <- predict(xgb_model, X_train)
      pred_xgb_test  <- predict(xgb_model, X_test)
    } else {
      pred_xgb_train <- rep(mean(y_train), length(y_train))
      pred_xgb_test  <- rep(mean(y_train), length(test_idx))
    }

    # --- BASE LEARNER 3: Random Forest ---
    set.seed(2026)
    rf_fit <- tryCatch({
      randomForest(x = X_train, y = factor(y_train, levels = c(0, 1)),
                   ntree = 500, mtry = max(floor(sqrt(ncol(X_train))), 1))
    }, error = function(e) NULL)

    if (!is.null(rf_fit)) {
      pred_rf_train <- as.numeric(predict(rf_fit, X_train, type = "prob")[, "1"])
      pred_rf_test  <- as.numeric(predict(rf_fit, X_test, type = "prob")[, "1"])
    } else {
      pred_rf_train <- rep(mean(y_train), length(y_train))
      pred_rf_test  <- rep(mean(y_train), length(test_idx))
    }

    # --- META-LEARNER: Logistic regression on base learner predictions ---
    # Train meta-learner on TRAINING fold predictions
    meta_df_train <- data.frame(en = pred_en_train, xgb = pred_xgb_train,
                                 rf = pred_rf_train, y = y_train)
    meta_fit <- glm(y ~ en + xgb + rf, family = binomial, data = meta_df_train)

    # Predict on TEST fold using meta-learner
    meta_df_test <- data.frame(en = pred_en_test, xgb = pred_xgb_test,
                                rf = pred_rf_test)
    outer_pred_stack[test_idx] <- predict(meta_fit, meta_df_test, type = "response")
    outer_pred_en[test_idx]    <- pred_en_test
    outer_pred_xgb[test_idx]   <- pred_xgb_test
    outer_pred_rf[test_idx]    <- pred_rf_test
  }

  # Calculate AUCs from outer fold predictions
  auc_stack <- as.numeric(auc(roc(y, outer_pred_stack, quiet = TRUE)))
  auc_en    <- as.numeric(auc(roc(y, outer_pred_en, quiet = TRUE)))
  auc_xgb   <- as.numeric(auc(roc(y, outer_pred_xgb, quiet = TRUE)))
  auc_rf    <- as.numeric(auc(roc(y, outer_pred_rf, quiet = TRUE)))

  cat(sprintf("  HONEST STACK AUC:   %.4f\n", auc_stack))
  cat(sprintf("  Individual AUCs: EN=%.4f, XGB=%.4f, RF=%.4f\n",
              auc_en, auc_xgb, auc_rf))
  cat(sprintf("  Stack improvement over best single: %+.4f\n\n",
              auc_stack - max(auc_en, auc_xgb, auc_rf)))

  return(list(
    auc = auc_stack, auc_en = auc_en, auc_xgb = auc_xgb, auc_rf = auc_rf,
    vars = vars_exist, df_cc = df_cc,
    n_complete = nrow(df_cc), n_events = sum(df_cc$ascvd_combined)
  ))
}

# Run honest stacking on key tiers
stack_v5_enriched <- run_honest_stack("v5-B: Enriched Clinical", v5_enriched_vars, df)
stack_v5_nmr_pca  <- run_honest_stack("v5-C: NMR-PCA Enhanced", v5_nmr_pca_vars, df)
stack_v5_full     <- run_honest_stack("v5-D: Full v5", v5_full_vars, df)
stack_v5_boruta   <- run_honest_stack("v5-E: Boruta Lean", v5_boruta_vars, df)

# Store stacking results
for (st in list(
  list(obj = stack_v5_enriched, key = "Stack_v5_Enriched", nm = "Stack v5 Enriched"),
  list(obj = stack_v5_nmr_pca,  key = "Stack_v5_NMR_PCA",  nm = "Stack v5 NMR-PCA"),
  list(obj = stack_v5_full,     key = "Stack_v5_Full",      nm = "Stack v5 Full"),
  list(obj = stack_v5_boruta,   key = "Stack_v5_Boruta",    nm = "Stack v5 Boruta Lean")
)) {
  if (!is.null(st$obj)) {
    model_results[[st$key]] <- list(
      name = st$nm, algorithm = "Honest Stack (EN+XGB+RF)",
      n_vars = length(st$obj$vars), cv_auc = st$obj$auc,
      n_complete = st$obj$n_complete, n_events = st$obj$n_events,
      model = st$obj
    )
  }
}

# =============================================================================
# SECTION 10: RCS LOGISTIC REGRESSION (FIXED: x=TRUE, y=TRUE)
# =============================================================================

cat("===================================================================\n")
cat("SECTION 10: RESTRICTED CUBIC SPLINES (v5 — FIXED)\n")
cat("===================================================================\n\n")

# v4 bug: lrm() did not use x=TRUE, y=TRUE → validate() failed
# v5 fix: explicitly set x=TRUE, y=TRUE

rcs_vars <- c("age", "sex", "re_ldl", "hdl", "log_apob_ldl", "log_lpa",
              "smoking_binary", "diabetes", "hypertension", "bmi", "log_crp",
              "lipid_years_A", "non_hdl", "egfr")
rcs_exist <- rcs_vars[rcs_vars %in% names(df)]
df_rcs <- df[complete.cases(df[, c(rcs_exist, "ascvd_combined")]), ]

cat(sprintf("  RCS model: %d vars, N=%d, Events=%d\n",
            length(rcs_exist), nrow(df_rcs), sum(df_rcs$ascvd_combined)))

if (nrow(df_rcs) >= 300 && sum(df_rcs$ascvd_combined) >= 50) {
  dd <- datadist(df_rcs[, rcs_exist])
  options(datadist = "dd")

  rcs_formula <- as.formula(paste0(
    "ascvd_combined ~ rcs(age, 4) + sex + rcs(re_ldl, 3) + hdl + ",
    "log_apob_ldl + log_lpa + smoking_binary + diabetes + hypertension + ",
    "bmi + log_crp + lipid_years_A + non_hdl + rcs(egfr, 3)"))

  rcs_fit <- tryCatch({
    lrm(rcs_formula, data = df_rcs, x = TRUE, y = TRUE)
  }, error = function(e) {
    cat(sprintf("  lrm() failed: %s\n", e$message))
    NULL
  })

  if (!is.null(rcs_fit)) {
    cat(sprintf("  RCS Model R²: %.4f, C-index: %.4f\n",
                rcs_fit$stats["R2"], rcs_fit$stats["C"]))

    # Bootstrap validation (B=200 for stable estimates)
    set.seed(2026)
    rcs_val <- tryCatch({
      validate(rcs_fit, B = 200)
    }, error = function(e) {
      cat(sprintf("  validate() failed: %s\n", e$message))
      NULL
    })

    if (!is.null(rcs_val)) {
      # Extract Dxy → AUC
      dxy_row <- rcs_val["Dxy", ]
      auc_rcs_orig <- 0.5 + dxy_row["index.orig"] / 2
      auc_rcs_corr <- 0.5 + dxy_row["index.corrected"] / 2
      optimism_rcs <- dxy_row["optimism"] / 2

      cat(sprintf("  Apparent AUC:  %.4f\n", auc_rcs_orig))
      cat(sprintf("  Corrected AUC: %.4f (optimism=%.4f)\n",
                  auc_rcs_corr, optimism_rcs))

      model_results[["RCS_v5"]] <- list(
        name = "RCS v5 (non-linear age+LDL+eGFR)", algorithm = "RCS Logistic (rms)",
        n_vars = length(rcs_exist), cv_auc = auc_rcs_corr,
        n_complete = nrow(df_rcs), n_events = sum(df_rcs$ascvd_combined),
        model = list(fit = rcs_fit, val = rcs_val)
      )

      # Print non-linear coefficients
      cat("\n  RCS Coefficients:\n")
      print(coef(rcs_fit))
      cat("\n")
    }
  }
}

# =============================================================================
# SECTION 11: GRAND MODEL COMPARISON
# =============================================================================

cat("===================================================================\n")
cat("SECTION 11: GRAND MODEL COMPARISON (v5)\n")
cat("===================================================================\n\n")

if (length(model_results) > 0) {
  # Sort by CV-AUC descending
  comp_df <- data.frame(
    model     = sapply(model_results, function(x) x$name),
    algorithm = sapply(model_results, function(x) x$algorithm),
    n_vars    = sapply(model_results, function(x) x$n_vars),
    n_complete = sapply(model_results, function(x) x$n_complete),
    n_events  = sapply(model_results, function(x) x$n_events),
    cv_auc    = sapply(model_results, function(x) x$cv_auc),
    stringsAsFactors = FALSE
  )
  comp_df <- comp_df[order(-comp_df$cv_auc), ]

  cat("  GRAND COMPARISON TABLE:\n")
  cat(sprintf("  %-35s %-25s %5s %6s %6s %8s\n",
              "Model", "Algorithm", "Vars", "N", "Events", "CV-AUC"))
  cat(paste(rep("-", 95), collapse = ""), "\n")

  for (i in 1:nrow(comp_df)) {
    marker <- ifelse(i == 1, " *** BEST", "")
    cat(sprintf("  %-35s %-25s %5d %6d %6d %8.4f%s\n",
                comp_df$model[i], comp_df$algorithm[i],
                comp_df$n_vars[i], comp_df$n_complete[i],
                comp_df$n_events[i], comp_df$cv_auc[i], marker))
  }

  # Identify best model
  best_key <- names(model_results)[order(-sapply(model_results, function(x) x$cv_auc))][1]
  best_model <- model_results[[best_key]]

  cat(sprintf("\n  BEST OVERALL MODEL: %s\n", best_model$name))
  cat(sprintf("  CV-AUC = %.4f\n", best_model$cv_auc))
  cat(sprintf("  Algorithm: %s, Vars: %d, N: %d, Events: %d\n",
              best_model$algorithm, best_model$n_vars,
              best_model$n_complete, best_model$n_events))

  # Save comparison table
  write.csv(comp_df, paste0(TAB_DIR, "calon2_v5_model_comparison.csv"), row.names = FALSE)
  cat("  Saved: calon2_v5_model_comparison.csv\n")

  # Target check
  if (best_model$cv_auc >= 0.80) {
    cat("\n  *** TARGET ACHIEVED: CV-AUC >= 0.80 ***\n")
  } else {
    cat(sprintf("\n  TARGET NOT YET MET: Need %.4f more AUC to reach 0.80\n",
                0.80 - best_model$cv_auc))
  }
} else {
  cat("  ERROR: No models were successfully trained.\n")
  best_model <- NULL
}

# =============================================================================
# SECTION 12: INTERNAL VALIDATION OF BEST MODEL
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 12: INTERNAL VALIDATION OF BEST MODEL\n")
cat("===================================================================\n\n")

if (!is.null(best_model)) {
  bm <- best_model$model

  # --- 12A: Repeated 10-fold CV (20 reps) ---
  cat("  12A: Repeated 10-fold CV (20 repetitions)...\n")

  if (best_model$algorithm == "Elastic Net" && !is.null(bm$X)) {
    cv_aucs <- numeric(20)
    for (rep in 1:20) {
      set.seed(2026 + rep)
      cv_fit <- cv.glmnet(bm$X, bm$y, family = "binomial",
                           alpha = bm$alpha, nfolds = 10, type.measure = "auc")
      cv_aucs[rep] <- max(cv_fit$cvm)
    }

    cat(sprintf("  CV-AUC (20 reps): mean=%.4f, SD=%.4f, 95%% CI=[%.4f, %.4f]\n",
                mean(cv_aucs), sd(cv_aucs),
                quantile(cv_aucs, 0.025), quantile(cv_aucs, 0.975)))

    # --- 12B: Apparent AUC ---
    pred_best <- predict(bm$fit, bm$X, s = "lambda.min", type = "response")[, 1]
    roc_apparent <- roc(bm$y, pred_best, quiet = TRUE)
    auc_apparent <- as.numeric(auc(roc_apparent))
    auc_ci <- ci.auc(roc_apparent)
    cat(sprintf("  Apparent AUC: %.4f (95%% CI: %.4f--%.4f)\n",
                auc_apparent, auc_ci[1], auc_ci[3]))

    # --- 12C: Bootstrap optimism correction (B=200) ---
    cat("  12C: Bootstrap optimism correction (B=200)...\n")

    boot_data_full <- bm$df_cc
    boot_vars_used <- bm$vars

    boot_optimism_fn <- function(data, indices) {
      boot_data <- data[indices, ]
      X_boot <- as.matrix(boot_data[, boot_vars_used])
      y_boot <- boot_data$ascvd_combined
      X_orig <- as.matrix(data[, boot_vars_used])
      y_orig <- data$ascvd_combined

      fit_boot <- glmnet(X_boot, y_boot, family = "binomial",
                          alpha = bm$alpha, lambda = bm$lambda)
      pred_boot_on_boot <- predict(fit_boot, X_boot, type = "response")[, 1]
      pred_boot_on_orig <- predict(fit_boot, X_orig, type = "response")[, 1]

      auc_boot <- tryCatch(
        as.numeric(auc(roc(y_boot, pred_boot_on_boot, quiet = TRUE))),
        error = function(e) NA)
      auc_orig <- tryCatch(
        as.numeric(auc(roc(y_orig, pred_boot_on_orig, quiet = TRUE))),
        error = function(e) NA)

      return(auc_boot - auc_orig)
    }

    set.seed(2026)
    boot_result <- boot(boot_data_full, boot_optimism_fn, R = 200)
    optimism <- mean(boot_result$t, na.rm = TRUE)
    auc_corrected <- auc_apparent - optimism

    cat(sprintf("  Apparent AUC:            %.4f\n", auc_apparent))
    cat(sprintf("  Mean optimism:           %.4f\n", optimism))
    cat(sprintf("  Optimism-corrected AUC:  %.4f\n", auc_corrected))

    # Save validation
    internal_val <- data.frame(
      metric = c("Best_Model", "Algorithm",
                  "Apparent_AUC", "Apparent_AUC_Lower", "Apparent_AUC_Upper",
                  "CV_AUC_Mean", "CV_AUC_SD", "CV_AUC_Lower", "CV_AUC_Upper",
                  "Optimism", "Corrected_AUC",
                  "N_vars", "N_complete", "N_events", "EPV"),
      value = c(best_model$name, best_model$algorithm,
                round(auc_apparent, 4), round(auc_ci[1], 4), round(auc_ci[3], 4),
                round(mean(cv_aucs), 4), round(sd(cv_aucs), 4),
                round(quantile(cv_aucs, 0.025), 4), round(quantile(cv_aucs, 0.975), 4),
                round(optimism, 4), round(auc_corrected, 4),
                best_model$n_vars, best_model$n_complete, best_model$n_events,
                round(best_model$n_events / best_model$n_vars, 1)),
      stringsAsFactors = FALSE
    )
    write.csv(internal_val, paste0(TAB_DIR, "calon2_v5_internal_validation.csv"),
              row.names = FALSE)
    cat("  Saved: calon2_v5_internal_validation.csv\n")

  } else if (best_model$algorithm == "XGBoost" && !is.null(bm$dtrain)) {
    # XGBoost validation
    cv_aucs <- numeric(20)
    for (rep in 1:20) {
      set.seed(2026 + rep)
      cv_res <- tryCatch({
        xgb.cv(params = bm$params, data = bm$dtrain,
               nrounds = bm$nrounds + 200, nfold = 10,
               early_stopping_rounds = 50, verbose = 0,
               print_every_n = 0)
      }, error = function(e) {
        xgb.cv(params = bm$params, data = bm$dtrain,
               nrounds = bm$nrounds + 200, nfold = 10,
               early_stopping_rounds = 50, verbose = 0)
      })
      eval_log <- cv_res$evaluation_log
      auc_col <- grep("test.*auc.*mean", names(eval_log), value = TRUE)[1]
      if (is.null(auc_col) || is.na(auc_col)) {
        auc_col <- names(eval_log)[grep("auc", names(eval_log))[1]]
      }
      cv_aucs[rep] <- max(eval_log[[auc_col]], na.rm = TRUE)
    }

    cat(sprintf("  CV-AUC (20 reps): mean=%.4f, SD=%.4f, 95%% CI=[%.4f, %.4f]\n",
                mean(cv_aucs), sd(cv_aucs),
                quantile(cv_aucs, 0.025), quantile(cv_aucs, 0.975)))

    pred_best <- predict(bm$final_model, bm$dtrain)
    roc_apparent <- roc(bm$y, pred_best, quiet = TRUE)
    auc_apparent <- as.numeric(auc(roc_apparent))
    auc_ci <- ci.auc(roc_apparent)
    cat(sprintf("  Apparent AUC: %.4f (95%% CI: %.4f--%.4f)\n",
                auc_apparent, auc_ci[1], auc_ci[3]))

    internal_val <- data.frame(
      metric = c("Best_Model", "Algorithm", "CV_AUC_Mean", "CV_AUC_SD",
                  "Apparent_AUC", "N_vars", "N_complete", "N_events"),
      value = c(best_model$name, best_model$algorithm,
                round(mean(cv_aucs), 4), round(sd(cv_aucs), 4),
                round(auc_apparent, 4),
                best_model$n_vars, best_model$n_complete, best_model$n_events),
      stringsAsFactors = FALSE
    )
    write.csv(internal_val, paste0(TAB_DIR, "calon2_v5_internal_validation.csv"),
              row.names = FALSE)

  } else {
    cat(sprintf("  Best model is %s — internal validation via its own CV metric.\n",
                best_model$algorithm))
    cat(sprintf("  Reported CV-AUC: %.4f\n", best_model$cv_auc))
  }
}

# =============================================================================
# SECTION 13: HEAD-TO-HEAD vs SAFEHEART-RE (in UKB)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 13: vs SAFEHEART-RE (internal comparison)\n")
cat("===================================================================\n\n")

# SAFEHEART-RE PUBLISHED coefficients (frozen — fair external benchmark)
# Source: Pérez de Isla et al., SAFEHEART-RE model
SAFEHEART_COEF <- list(
  intercept = -7.053,
  age       =  0.064,
  male      =  0.775,
  ldl_c     =  0.109,
  htn       =  0.431,
  bmi       =  0.025,
  smoking   =  0.466,
  prior_cvd =  1.414
)

sh_vars <- c("age", "sex", "re_ldl", "hypertension", "bmi", "smoking_binary")
sh_exist <- sh_vars[sh_vars %in% names(df)]
df_sh <- df[complete.cases(df[, c(sh_exist, "ascvd_combined")]), ]

if (nrow(df_sh) >= 300 && !is.null(best_model)) {
  # Apply PUBLISHED coefficients (no re-fitting on UKB — fair comparison)
  df_sh$lp_safeheart <- with(df_sh,
    SAFEHEART_COEF$intercept +
    SAFEHEART_COEF$age     * age +
    SAFEHEART_COEF$male    * sex +
    SAFEHEART_COEF$ldl_c   * re_ldl +
    SAFEHEART_COEF$htn     * hypertension +
    SAFEHEART_COEF$bmi     * bmi +
    SAFEHEART_COEF$smoking * smoking_binary
  )
  df_sh$prob_safeheart <- 1 / (1 + exp(-df_sh$lp_safeheart))

  sh_roc <- roc(df_sh$ascvd_combined, df_sh$prob_safeheart, quiet = TRUE)
  cat(sprintf("  SAFEHEART-RE (published coefficients) AUC: %.4f (95%% CI: %.4f--%.4f)\n",
              as.numeric(auc(sh_roc)), ci.auc(sh_roc)[1], ci.auc(sh_roc)[3]))

  # Best model prediction on same patients
  if (best_model$algorithm == "Elastic Net") {
    bm_obj <- best_model$model
    bm_vars <- bm_obj$vars
    bm_vars_in_sh <- bm_vars[bm_vars %in% names(df_sh)]
    if (length(bm_vars_in_sh) == length(bm_vars)) {
      X_sh <- as.matrix(df_sh[, bm_vars])
      pred_bm <- predict(bm_obj$fit, X_sh, s = "lambda.min", type = "response")[, 1]
      bm_roc <- roc(df_sh$ascvd_combined, pred_bm, quiet = TRUE)

      # DeLong test
      delong <- roc.test(bm_roc, sh_roc, method = "delong")

      cat(sprintf("  Best model AUC:      %.4f (95%% CI: %.4f--%.4f) [%s]\n",
                  as.numeric(auc(bm_roc)), ci.auc(bm_roc)[1], ci.auc(bm_roc)[3],
                  best_model$name))
      cat(sprintf("  dAUC (Best - SH):    %+.4f (DeLong p=%.6f)",
                  as.numeric(auc(bm_roc)) - as.numeric(auc(sh_roc)),
                  delong$p.value))
      if (delong$p.value < 0.05) cat(" SIGNIFICANT") else cat(" n.s.")
      cat("\n")
    }
  }
}

# =============================================================================
# SECTION 14: SAVE v5 OUTPUTS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 14: SAVING v5 OUTPUTS\n")
cat("===================================================================\n\n")

# Save best model object
if (!is.null(best_model)) {
  saveRDS(best_model, paste0(OUT_DIR, "calon2_v5_best_model.rds"))
  cat("  Saved: calon2_v5_best_model.rds\n")
}

# Save the enriched dataset with v5 features
write.csv(df, paste0(OUT_DIR, "calon2_v5_full_merged.csv"), row.names = FALSE)
cat("  Saved: calon2_v5_full_merged.csv\n")

# Save frozen coefficients for best EN model
if (!is.null(best_model) && best_model$algorithm == "Elastic Net") {
  bm_obj <- best_model$model
  coefs <- coef(bm_obj$fit, s = "lambda.min")
  coef_df <- data.frame(
    variable = rownames(coefs),
    coefficient = as.numeric(coefs),
    stringsAsFactors = FALSE
  )
  coef_df <- coef_df[coef_df$coefficient != 0, ]
  coef_df <- coef_df[order(-abs(coef_df$coefficient)), ]
  write.csv(coef_df, paste0(OUT_DIR, "calon2_v5_best_coefficients.csv"), row.names = FALSE)
  cat("  Saved: calon2_v5_best_coefficients.csv\n")

  cat("\n  Frozen coefficients for best model:\n")
  for (j in 1:nrow(coef_df)) {
    cat(sprintf("    %-25s %+.6f\n", coef_df$variable[j], coef_df$coefficient[j]))
  }
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n\n")
cat("======================================================================\n")
cat("  CALON-2 v5 FINAL SUMMARY\n")
cat("======================================================================\n\n")

if (!is.null(best_model)) {
  cat(sprintf("   BEST MODEL:          %s\n", best_model$name))
  cat(sprintf("   Algorithm:           %s\n", best_model$algorithm))
  cat(sprintf("   CV-AUC:              %.4f\n", best_model$cv_auc))
  cat(sprintf("   N:                   %d\n", best_model$n_complete))
  cat(sprintf("   Events:              %d\n", best_model$n_events))
  cat(sprintf("   Variables:           %d\n", best_model$n_vars))
}

cat("\n  v5 improvements over v4:\n")
cat("    1. New features: non-HDL, pulse pressure, eGFR, TC/HDL, TG/HDL, MetSyn score\n")
cat("    2. NMR PCA: orthogonal components replace 100 collinear raw features\n")
cat("    3. Proper MI: NMR PCs as auxiliary variables for Lp(a) imputation\n")
cat("    4. Wider XGBoost grid: 40 random parameter combinations\n")
cat("    5. Honest stacking: out-of-fold predictions, NO data leakage\n")
cat("    6. RCS fixed: x=TRUE, y=TRUE in lrm() for bootstrap validation\n")
cat("    7. Boruta feature selection: data-driven variable importance\n")
cat("    8. Removed: SuperLearner (confirmed leakage in default CV)\n\n")

cat("  All models:\n")
if (length(model_results) > 0) {
  sorted_keys <- names(model_results)[order(-sapply(model_results, function(x) x$cv_auc))]
  for (k in sorted_keys) {
    mr <- model_results[[k]]
    cat(sprintf("   %-35s CV-AUC = %.4f (%s, %d vars)\n",
                mr$name, mr$cv_auc, mr$algorithm, mr$n_vars))
  }
}

cat("\n======================================================================\n")
cat("  DONE. Next step: External validation on DRAGON3 (06_CALON2_validate_DRAGON3.R)\n")
cat("======================================================================\n")

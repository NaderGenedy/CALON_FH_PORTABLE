################################################################################
#                                                                              #
#  CALON-2 Script 15: INTERNAL-EXTERNAL CROSS-VALIDATION (IECV)              #
#                                                                              #
#  Methodology (Steyerberg & Harrell, Debray et al. 2015):                   #
#    - Leave-one-dataset-out: develop on 2 datasets, validate on the 3rd     #
#    - Rotate: each dataset serves as external validation exactly once        #
#    - Pool external AUCs across folds for overall performance estimate       #
#    - This uses ALL data efficiently without data leakage                    #
#                                                                              #
#  Three folds:                                                               #
#    Fold 1: Train on UKB + Wales,    Validate on Dragon 3                   #
#    Fold 2: Train on UKB + Dragon 3, Validate on Wales                      #
#    Fold 3: Train on Wales + Dragon 3, Validate on UKB                      #
#                                                                              #
#  Only penalised logistic regression (ridge/elastic net).                    #
#  Variables restricted to those available in ALL 3 datasets.                 #
#                                                                              #
#  Author:  Dr Nader Genedy                                                    #
#  Date:    March 2026                                                         #
#                                                                              #
################################################################################

rm(list = ls())
set.seed(2026)

# =============================================================================
# PACKAGES
# =============================================================================

required_packages <- c("glmnet", "pROC", "dplyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("\n")
cat("======================================================================\n")
cat("  CALON-2 Script 15: INTERNAL-EXTERNAL CROSS-VALIDATION (IECV)       \n")
cat("  Leave-one-dataset-out with 3 cohorts                                \n")
cat("======================================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

CLOUD_PATH <- "/cloud/project/"
LOCAL_PATH <- "C:/Users/nader/Downloads/calon_ukb_pipeline/"

if (dir.exists(CLOUD_PATH) && file.exists(paste0(CLOUD_PATH, "project.Rproj"))) {
  BASE_DIR <- CLOUD_PATH
  cat("  Environment: Posit Cloud\n")
} else if (.Platform$OS.type == "windows" && dir.exists(LOCAL_PATH)) {
  BASE_DIR <- LOCAL_PATH
  cat("  Environment: Local Windows\n")
} else {
  BASE_DIR <- paste0(normalizePath(getwd(), winslash = "/"), "/")
  cat(sprintf("  Environment: Auto-detected (%s)\n", BASE_DIR))
}
cat(sprintf("  BASE_DIR = %s\n\n", BASE_DIR))

OUT_DIR <- paste0(BASE_DIR, "output/")
TAB_DIR <- paste0(OUT_DIR, "tables/")
for (d in c(TAB_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# =============================================================================
# SAFEHEART-RE COEFFICIENTS
# =============================================================================

SAFEHEART_COEF <- list(
  intercept = -7.053,
  betas = c(age = 0.064, sex = 0.775, re_ldl = 0.109,
            hypertension = 0.431, bmi = 0.025, smoking_binary = 0.466)
)

# =============================================================================
# HELPERS
# =============================================================================

safe_col <- function(df, col, as_num = FALSE) {
  if (col %in% names(df)) {
    vals <- df[[col]]
    if (as_num) vals <- suppressWarnings(as.numeric(vals))
    return(vals)
  }
  rep(NA_real_, nrow(df))
}

safe_auc <- function(y, pred) {
  valid <- !is.na(y) & !is.na(pred)
  if (sum(valid) < 20 || length(unique(y[valid])) < 2)
    return(list(auc = NA, lo = NA, hi = NA, n = sum(valid), events = sum(y[valid] == 1)))
  roc_obj <- pROC::roc(y[valid], pred[valid], quiet = TRUE)
  ci <- pROC::ci.auc(roc_obj)
  list(auc = as.numeric(pROC::auc(roc_obj)), lo = ci[1], hi = ci[3],
       n = sum(valid), events = sum(y[valid] == 1), roc = roc_obj)
}

apply_safeheart <- function(d, impute_bmi = 27.0) {
  lp <- rep(SAFEHEART_COEF$intercept, nrow(d))
  for (v in names(SAFEHEART_COEF$betas)) {
    vals <- d[[v]]
    if (v == "bmi") vals[is.na(vals)] <- impute_bmi
    vals[is.na(vals)] <- median(vals, na.rm = TRUE)
    lp <- lp + SAFEHEART_COEF$betas[v] * vals
  }
  1 / (1 + exp(-lp))
}

apply_frozen <- function(d, coefs, vars) {
  avail <- vars[vars %in% names(d)]
  lp <- rep(coefs["intercept"], nrow(d))
  valid <- rep(TRUE, nrow(d))
  for (v in avail) {
    if (v %in% names(coefs)) {
      vals <- d[[v]]
      valid <- valid & !is.na(vals)
      vals[is.na(vals)] <- 0
      lp <- lp + coefs[v] * vals
    }
  }
  valid <- valid & !is.na(d$ascvd_combined)
  pred <- 1 / (1 + exp(-lp))
  pred[!valid] <- NA  # ensure invalid rows return NA
  list(pred = pred, valid = valid)
}

# Simple median/mode imputation for development set
impute_simple <- function(df, vars) {
  for (v in vars) {
    if (!v %in% names(df) || sum(is.na(df[[v]])) == 0) next
    if (v %in% c("sex","smoking_binary","diabetes","hypertension","gene_apob")) {
      tab <- table(df[[v]], useNA = "no")
      if (length(tab) > 0) df[[v]][is.na(df[[v]])] <- as.numeric(names(tab)[which.max(tab)])
    } else {
      med <- median(df[[v]], na.rm = TRUE)
      if (!is.na(med)) df[[v]][is.na(df[[v]])] <- med
    }
  }
  df
}

# =============================================================================
# SECTION 1: LOAD & HARMONISE
# =============================================================================

cat("===================================================================\n")
cat("SECTION 1: LOAD & HARMONISE ALL DATASETS\n")
cat("===================================================================\n\n")

# --- UKB ---
merged_file <- NULL
search_paths <- c(paste0(OUT_DIR, "calon2_full_merged.csv"),
                  paste0(BASE_DIR, "calon2_full_merged.csv"),
                  "/cloud/project/output/calon2_full_merged.csv",
                  "/cloud/project/calon2_full_merged.csv")
for (sp in search_paths) if (file.exists(sp)) { merged_file <- sp; break }
if (is.null(merged_file)) stop("ERROR: calon2_full_merged.csv not found")
ukb_raw <- read.csv(merged_file, stringsAsFactors = FALSE)

ukb_h <- data.frame(
  age = ukb_raw$age, sex = ukb_raw$sex,
  re_ldl = ukb_raw$re_ldl, hdl = ukb_raw$hdl,
  tc = ukb_raw$tc, trig = ukb_raw$trig,
  apob = safe_col(ukb_raw, "apob", TRUE),
  lpa  = safe_col(ukb_raw, "lpa", TRUE),
  smoking_binary = ukb_raw$ever_smoked,
  diabetes = ukb_raw$diabetes,
  hypertension = ukb_raw$hypertension,
  bmi = ukb_raw$bmi,
  gene_apob = as.integer(ukb_raw$gene == "APOB"),
  ascvd_combined = ukb_raw$ascvd_combined,
  dataset = "UKB",
  stringsAsFactors = FALSE
)
ukb_h$non_hdl <- ukb_h$tc - ukb_h$hdl
ukb_h$tc_hdl_ratio <- ifelse(!is.na(ukb_h$hdl) & ukb_h$hdl > 0, ukb_h$tc / ukb_h$hdl, NA)
ukb_h$log_lpa <- ifelse(!is.na(ukb_h$lpa) & ukb_h$lpa > 0, log(ukb_h$lpa + 1), NA)

cat(sprintf("  UKB: N=%d, Events=%d (%.1f%%)\n",
            nrow(ukb_h), sum(ukb_h$ascvd_combined), 100*mean(ukb_h$ascvd_combined)))

# --- Dragon 3 ---
search_dirs <- unique(c(BASE_DIR, getwd(), "/cloud/project/"))
dragon_file <- NULL
for (sd in search_dirs) {
  if (!dir.exists(sd)) next
  found <- list.files(sd, pattern = "DRAGON.*3\\.csv$", full.names = TRUE, ignore.case = TRUE)
  if (length(found) > 0) { dragon_file <- found[1]; break }
}
if (is.null(dragon_file)) stop("ERROR: DRAGON_3.csv not found")
dragon_raw <- read.csv(dragon_file, stringsAsFactors = FALSE)
has_data <- !is.na(dragon_raw$ASCVD_combined) & dragon_raw$ASCVD_combined != "" &
  !is.na(dragon_raw$Currentage) & dragon_raw$Currentage != ""
dragon_raw <- dragon_raw[has_data, ]

d3 <- data.frame(row.names = 1:nrow(dragon_raw))
d3$age <- as.numeric(dragon_raw$Currentage)
d3$sex <- ifelse(dragon_raw$Gender == "M", 1, ifelse(dragon_raw$Gender == "F", 0, NA))
d3$re_ldl <- as.numeric(dragon_raw$LDL_1)
d3$hdl <- as.numeric(dragon_raw$HDL_1)
d3$tc <- as.numeric(dragon_raw$TC_1)
d3$trig <- as.numeric(dragon_raw$TRG_1)
for (suffix in c("_2","_3","_4")) {
  for (lip in list(c("LDL","re_ldl"),c("HDL","hdl"),c("TC","tc"),c("TRG","trig"))) {
    col <- paste0(lip[1], suffix)
    if (col %in% names(dragon_raw)) {
      fill <- is.na(d3[[lip[2]]]) & !is.na(suppressWarnings(as.numeric(dragon_raw[[col]])))
      d3[[lip[2]]][fill] <- as.numeric(dragon_raw[[col]][fill])
    }
  }
}
for (pair in list(c("LastLDL","re_ldl"),c("LastHDL","hdl"),c("LastTC","tc"),c("LastTrigs","trig"))) {
  if (pair[1] %in% names(dragon_raw)) {
    fill <- is.na(d3[[pair[2]]]) & !is.na(suppressWarnings(as.numeric(dragon_raw[[pair[1]]])))
    d3[[pair[2]]][fill] <- as.numeric(dragon_raw[[pair[1]]][fill])
  }
}
d3$apob <- suppressWarnings(as.numeric(dragon_raw$ApoB))
d3$lpa <- suppressWarnings(as.numeric(dragon_raw$Lpa))
if (sum(!is.na(d3$lpa)) < 20) d3$lpa <- suppressWarnings(as.numeric(dragon_raw$Lpa_1))
if (sum(!is.na(d3$lpa)) < 20) d3$lpa <- suppressWarnings(as.numeric(dragon_raw$Lpaunitsmgl))
smk <- dragon_raw$Smoking_binary
d3$smoking_binary <- ifelse(smk %in% c("1","2"), 1, ifelse(smk %in% c("0"), 0, NA))
dm <- dragon_raw$Diabetes_binary
d3$diabetes <- ifelse(dm %in% c(1,"1","Yes",TRUE), 1, ifelse(dm %in% c(0,"0","No",FALSE,""), 0, NA))
bp_treat <- dragon_raw$onBPtreat
d3$hypertension <- ifelse(bp_treat %in% c("Y","y","Yes","1",1), 1,
                         ifelse(bp_treat %in% c("N","n","No","0",0), 0, NA))
# Augment hypertension from BP readings
d3$sbp <- NA_real_; d3$dbp <- NA_real_
for (i in seq_len(nrow(dragon_raw))) {
  bp_str <- dragon_raw$BP[i]
  if (!is.na(bp_str) && grepl("/", bp_str)) {
    parts <- strsplit(bp_str, "/")[[1]]
    if (length(parts) == 2) {
      s <- suppressWarnings(as.numeric(trimws(parts[1])))
      dd <- suppressWarnings(as.numeric(trimws(parts[2])))
      if (!is.na(s) && s > 50 && s < 300) d3$sbp[i] <- s
      if (!is.na(dd) && dd > 20 && dd < 200) d3$dbp[i] <- dd
    }
  }
}
bp_sys <- suppressWarnings(as.numeric(dragon_raw$BloodPressureSystolic))
bp_dia <- suppressWarnings(as.numeric(dragon_raw$BloodPressureDiastolic))
d3$sbp[is.na(d3$sbp) & !is.na(bp_sys)] <- bp_sys[is.na(d3$sbp) & !is.na(bp_sys)]
d3$dbp[is.na(d3$dbp) & !is.na(bp_dia)] <- bp_dia[is.na(d3$dbp) & !is.na(bp_dia)]
high_bp <- (!is.na(d3$sbp) & d3$sbp >= 140) | (!is.na(d3$dbp) & d3$dbp >= 90)
d3$hypertension[high_bp & (is.na(d3$hypertension) | d3$hypertension == 0)] <- 1
d3$bmi <- suppressWarnings(as.numeric(dragon_raw$BMI))
d3$bmi[!is.na(d3$bmi) & (d3$bmi < 10 | d3$bmi > 80)] <- NA
d3$gene_apob <- as.integer(grepl("APOB", dragon_raw$Mutation1, ignore.case = TRUE))
d3$ascvd_combined <- as.integer(dragon_raw$ASCVD_combined %in% c(1, "1", "1.0"))
d3$non_hdl <- d3$tc - d3$hdl
d3$tc_hdl_ratio <- ifelse(!is.na(d3$hdl) & d3$hdl > 0, d3$tc / d3$hdl, NA)
d3$log_lpa <- ifelse(!is.na(d3$lpa) & d3$lpa > 0, log(d3$lpa + 1), NA)
d3$dataset <- "Dragon3"

cat(sprintf("  Dragon 3: N=%d, Events=%d (%.1f%%)\n",
            nrow(d3), sum(d3$ascvd_combined), 100*mean(d3$ascvd_combined)))

# --- Wales ---
wales_file <- NULL
for (sd in search_dirs) {
  if (!dir.exists(sd)) next
  found <- list.files(sd, pattern = "WALES.*FH.*CLEAN", full.names = TRUE, ignore.case = TRUE)
  if (length(found) > 0) { wales_file <- found[1]; break }
}
if (is.null(wales_file)) {
  wales_names <- c("WALES_FH_CLEANED (1) - Copy.csv","WALES_FH_CLEANED (1).csv","WALES_FH_CLEANED.csv")
  for (sd in search_dirs) {
    for (wn in wales_names) {
      wf <- paste0(sd, wn)
      if (file.exists(wf)) { wales_file <- wf; break }
    }
    if (!is.null(wales_file)) break
  }
}
if (is.null(wales_file)) stop("ERROR: Wales FH cleaned CSV not found")
wales_raw <- read.csv(wales_file, stringsAsFactors = FALSE)
fh_positive <- !is.na(wales_raw$Mutation1) & trimws(wales_raw$Mutation1) != ""
wales_raw <- wales_raw[fh_positive, ]

wales_h <- data.frame(row.names = 1:nrow(wales_raw))
WALES_REF_DATE <- as.Date("2025-01-01")
if ("DOB" %in% names(wales_raw)) {
  dob <- as.Date(wales_raw$DOB, format = "%Y-%m-%d")
  if (sum(!is.na(dob)) < 100) dob <- as.Date(wales_raw$DOB, format = "%d/%m/%Y")
  wales_h$age <- as.numeric(difftime(WALES_REF_DATE, dob, units = "days")) / 365.25
} else { wales_h$age <- NA }
wales_h$sex <- ifelse(wales_raw$Gender == "M", 1, ifelse(wales_raw$Gender == "F", 0, NA))
wales_h$re_ldl <- suppressWarnings(as.numeric(wales_raw$LDL.1))
wales_h$hdl <- suppressWarnings(as.numeric(wales_raw$HDL.1))
wales_h$tc <- suppressWarnings(as.numeric(wales_raw$TC.1))
wales_h$trig <- suppressWarnings(as.numeric(wales_raw$TRG.1))
for (suffix in c(".2",".3",".4")) {
  for (lip in list(c("LDL","re_ldl"),c("HDL","hdl"),c("TC","tc"),c("TRG","trig"))) {
    col <- paste0(lip[1], suffix)
    if (col %in% names(wales_raw)) {
      fill <- is.na(wales_h[[lip[2]]]) & !is.na(suppressWarnings(as.numeric(wales_raw[[col]])))
      wales_h[[lip[2]]][fill] <- suppressWarnings(as.numeric(wales_raw[[col]][fill]))
    }
  }
}
smk_w <- wales_raw$Smoking
wales_h$smoking_binary <- ifelse(smk_w %in% c(1,"1","Yes","yes","Current","current","Ex","ex",
                                    "Former","former","Ex-smoker","Current smoker","Previous"), 1,
                            ifelse(smk_w %in% c(0,"0","No","no","Never","never","Non-smoker"), 0, NA))
dm_w <- wales_raw$Diabetes
wales_h$diabetes <- ifelse(dm_w %in% c(1,"1","Yes","yes",TRUE), 1,
                     ifelse(dm_w %in% c(0,"0","No","no",FALSE,""), 0, NA))
bp_med <- wales_raw$BloodPressureMedication
wales_h$hypertension <- ifelse(bp_med %in% c(1,"1","Yes","yes",TRUE), 1,
                         ifelse(bp_med %in% c(0,"0","No","no",FALSE,""), 0, NA))
sbp_w <- suppressWarnings(as.numeric(wales_raw$BloodPressureSystolic))
dbp_w <- suppressWarnings(as.numeric(wales_raw$BloodPressureDiastolic))
high_bp_w <- (!is.na(sbp_w) & sbp_w >= 140) | (!is.na(dbp_w) & dbp_w >= 90)
wales_h$hypertension[high_bp_w & (is.na(wales_h$hypertension) | wales_h$hypertension == 0)] <- 1
wales_h$bmi <- suppressWarnings(as.numeric(wales_raw$BMI))
wales_h$gene_apob <- as.integer(grepl("APOB", wales_raw$Mutation1, ignore.case = TRUE))
wales_h$ascvd_combined <- as.integer(wales_raw$ascvd_combine %in% c(1, "1", "1.0"))
wales_h$non_hdl <- wales_h$tc - wales_h$hdl
wales_h$tc_hdl_ratio <- ifelse(!is.na(wales_h$hdl) & wales_h$hdl > 0, wales_h$tc / wales_h$hdl, NA)
wales_h$dataset <- "Wales"

cat(sprintf("  Wales: N=%d, Events=%d (%.1f%%)\n",
            nrow(wales_h), sum(wales_h$ascvd_combined), 100*mean(wales_h$ascvd_combined)))

# =============================================================================
# SECTION 2: DEFINE COMMON VARIABLE TIERS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 2: VARIABLE TIERS (common across all 3 datasets)\n")
cat("===================================================================\n\n")

# IECV tiers must use variables available in ALL datasets
# Wales has NO ApoB, ~2% Lpa — so ApoB/Lpa tiers only use UKB + D3

# Universal tiers (all 3 datasets have these)
universal_tiers <- list(
  "T1_safeheart" = c("age", "sex", "re_ldl", "hypertension", "bmi", "smoking_binary"),
  "T2_core" = c("age", "sex", "re_ldl", "hdl", "smoking_binary", "diabetes", "hypertension"),
  "T3_clinical" = c("age", "sex", "re_ldl", "hdl", "trig", "bmi",
                     "smoking_binary", "diabetes", "hypertension"),
  "T4_expanded" = c("age", "sex", "re_ldl", "hdl", "trig", "non_hdl",
                     "bmi", "smoking_binary", "diabetes", "hypertension")
)

# ApoB tiers (UKB + D3 only — Wales excluded from these folds)
apob_tiers <- list(
  "T5_apob" = c("age", "sex", "re_ldl", "hdl", "apob",
                  "smoking_binary", "diabetes", "hypertension"),
  "T6_apob_lpa" = c("age", "sex", "re_ldl", "hdl", "apob", "log_lpa",
                      "smoking_binary", "diabetes", "hypertension")
)

cat("  UNIVERSAL TIERS (3-fold IECV across all datasets):\n")
for (tn in names(universal_tiers)) {
  cat(sprintf("    %s: %d vars — %s\n", tn, length(universal_tiers[[tn]]),
              paste(universal_tiers[[tn]], collapse = ", ")))
}
cat("\n  ApoB TIERS (2-fold IECV: UKB + Dragon 3 only):\n")
for (tn in names(apob_tiers)) {
  cat(sprintf("    %s: %d vars — %s\n", tn, length(apob_tiers[[tn]]),
              paste(apob_tiers[[tn]], collapse = ", ")))
}

# =============================================================================
# SECTION 3: SAFEHEART-RE BASELINES
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 3: SAFEHEART-RE BASELINES\n")
cat("===================================================================\n\n")

sh_ukb   <- safe_auc(ukb_h$ascvd_combined, apply_safeheart(ukb_h))
sh_d3    <- safe_auc(d3$ascvd_combined, apply_safeheart(d3))
sh_wales <- safe_auc(wales_h$ascvd_combined, apply_safeheart(wales_h))

cat(sprintf("  SAFEHEART [UKB]:     AUC=%.4f (%.4f-%.4f), N=%d, Ev=%d\n",
            sh_ukb$auc, sh_ukb$lo, sh_ukb$hi, sh_ukb$n, sh_ukb$events))
cat(sprintf("  SAFEHEART [Dragon3]: AUC=%.4f (%.4f-%.4f), N=%d, Ev=%d\n",
            sh_d3$auc, sh_d3$lo, sh_d3$hi, sh_d3$n, sh_d3$events))
cat(sprintf("  SAFEHEART [Wales]:   AUC=%.4f (%.4f-%.4f), N=%d, Ev=%d\n",
            sh_wales$auc, sh_wales$lo, sh_wales$hi, sh_wales$n, sh_wales$events))

sh_aucs <- c(UKB = sh_ukb$auc, Dragon3 = sh_d3$auc, Wales = sh_wales$auc)

# =============================================================================
# SECTION 4: 3-FOLD IECV (UNIVERSAL TIERS)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 4: 3-FOLD IECV (Universal tiers — all 3 datasets)\n")
cat("===================================================================\n\n")

# Define the 3 datasets as a list
datasets <- list(UKB = ukb_h, Dragon3 = d3, Wales = wales_h)
ds_names <- names(datasets)

all_iecv_results <- data.frame()

# Store per-fold validation predictions for DeLong tests later
# Key: "tier|alpha" -> list of per-fold results with new_pred, sh_pred, outcome, dataset name
iecv_fold_preds <- list()

for (tier_name in names(universal_tiers)) {
  tier_vars <- universal_tiers[[tier_name]]
  cat(sprintf("  === %s (%d vars) ===\n", tier_name, length(tier_vars)))

  for (alpha_val in c(0, 0.10, 0.25, 0.50, 1.0)) {

    fold_aucs <- numeric(0)
    fold_deltas <- numeric(0)
    fold_details <- list()

    for (val_idx in 1:3) {
      val_name <- ds_names[val_idx]
      train_names <- ds_names[-val_idx]

      # Combine training datasets
      train_dfs <- list()
      for (tn in train_names) {
        # Select only the columns we need + outcome
        cols_needed <- c(tier_vars, "ascvd_combined")
        cols_avail <- cols_needed[cols_needed %in% names(datasets[[tn]])]
        train_dfs[[tn]] <- datasets[[tn]][, cols_avail, drop = FALSE]
      }

      # Ensure same columns before rbind
      common_cols <- Reduce(intersect, lapply(train_dfs, names))
      for (tn in names(train_dfs)) train_dfs[[tn]] <- train_dfs[[tn]][, common_cols, drop = FALSE]
      train_df <- do.call(rbind, train_dfs)

      # Impute training set (median/mode)
      train_df <- impute_simple(train_df, tier_vars)

      # Complete cases
      avail <- tier_vars[tier_vars %in% names(train_df)]
      cc <- complete.cases(train_df[, c(avail, "ascvd_combined")])
      train_cc <- train_df[cc, ]

      if (nrow(train_cc) < 100 || sum(train_cc$ascvd_combined) < 20) {
        cat(sprintf("    Fold %d (val=%s): SKIP — insufficient training data\n", val_idx, val_name))
        next
      }

      X_train <- as.matrix(train_cc[, avail])
      y_train <- train_cc$ascvd_combined

      # Fit penalised logistic
      set.seed(2026 + val_idx)
      cv_fit <- tryCatch(
        cv.glmnet(X_train, y_train, family = "binomial", alpha = alpha_val,
                  nfolds = 10, type.measure = "auc"),
        error = function(e) NULL
      )
      if (is.null(cv_fit)) {
        cat(sprintf("    Fold %d (val=%s): SKIP — glmnet failed\n", val_idx, val_name))
        next
      }

      coefs <- as.numeric(coef(cv_fit, s = "lambda.min"))
      names(coefs) <- c("intercept", avail)

      # Predict on validation dataset
      val_df <- datasets[[val_name]]
      val_res <- apply_frozen(val_df, coefs, avail)

      if (sum(val_res$valid) < 30) {
        cat(sprintf("    Fold %d (val=%s): SKIP — <30 valid in validation\n", val_idx, val_name))
        next
      }

      val_auc_obj <- safe_auc(val_df$ascvd_combined[val_res$valid], val_res$pred[val_res$valid])

      if (is.na(val_auc_obj$auc)) next

      delta <- val_auc_obj$auc - sh_aucs[val_name]
      fold_aucs <- c(fold_aucs, val_auc_obj$auc)
      fold_deltas <- c(fold_deltas, delta)

      cat(sprintf("    Fold %d: Train=%s (N=%d), Val=%s: AUC=%.4f (%.4f-%.4f) [SH=%.4f, d=%+.4f]\n",
                  val_idx, paste(train_names, collapse="+"), nrow(train_cc),
                  val_name, val_auc_obj$auc, val_auc_obj$lo, val_auc_obj$hi,
                  sh_aucs[val_name], delta))

      # Store per-fold predictions for DeLong test
      sh_pred_val <- apply_safeheart(val_df)
      both_ok <- val_res$valid & !is.na(sh_pred_val)
      fold_key <- paste0(tier_name, "|", alpha_val)
      if (is.null(iecv_fold_preds[[fold_key]])) iecv_fold_preds[[fold_key]] <- list()
      iecv_fold_preds[[fold_key]][[val_name]] <- list(
        new_pred = val_res$pred[both_ok],
        sh_pred  = sh_pred_val[both_ok],
        outcome  = val_df$ascvd_combined[both_ok],
        n        = sum(both_ok)
      )

      fold_details[[val_name]] <- list(auc = val_auc_obj$auc, lo = val_auc_obj$lo,
                                        hi = val_auc_obj$hi, delta = delta,
                                        n = val_auc_obj$n, events = val_auc_obj$events)
    }

    if (length(fold_aucs) == 0) next

    # Summary across IECV folds
    # Weighted average by number of events in each validation set
    weights <- sapply(fold_details, function(x) x$events)
    pooled_auc <- sum(sapply(fold_details, function(x) x$auc) * weights) / sum(weights)
    pooled_delta <- sum(sapply(fold_details, function(x) x$delta) * weights) / sum(weights)
    mean_auc <- mean(fold_aucs)
    mean_delta <- mean(fold_deltas)

    cat(sprintf("    >> Pooled IECV AUC=%.4f (weighted), Mean delta=%+.4f, %d/%d folds\n",
                pooled_auc, pooled_delta, length(fold_aucs), 3))

    # Store per-fold AUCs
    ukb_auc <- if ("UKB" %in% names(fold_details)) fold_details$UKB$auc else NA
    d3_auc <- if ("Dragon3" %in% names(fold_details)) fold_details$Dragon3$auc else NA
    wales_auc <- if ("Wales" %in% names(fold_details)) fold_details$Wales$auc else NA

    row <- data.frame(
      tier = tier_name, alpha = alpha_val, n_vars = length(tier_vars),
      iecv_type = "3fold_universal",
      ukb_auc = ukb_auc, d3_auc = d3_auc, wales_auc = wales_auc,
      pooled_auc = pooled_auc, mean_auc = mean_auc,
      mean_delta = mean_delta, pooled_delta = pooled_delta,
      n_folds = length(fold_aucs),
      stringsAsFactors = FALSE
    )
    all_iecv_results <- rbind(all_iecv_results, row)
  }
}

# =============================================================================
# SECTION 5: 2-FOLD IECV (ApoB TIERS — UKB + Dragon 3 only)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 5: 2-FOLD IECV (ApoB tiers — UKB + Dragon 3)\n")
cat("===================================================================\n\n")

apob_datasets <- list(UKB = ukb_h, Dragon3 = d3)
apob_ds_names <- names(apob_datasets)
apob_sh_aucs <- c(UKB = sh_ukb$auc, Dragon3 = sh_d3$auc)

for (tier_name in names(apob_tiers)) {
  tier_vars <- apob_tiers[[tier_name]]
  cat(sprintf("  === %s (%d vars) ===\n", tier_name, length(tier_vars)))

  for (alpha_val in c(0, 0.10, 0.25, 0.50, 1.0)) {

    fold_aucs <- numeric(0)
    fold_deltas <- numeric(0)
    fold_details <- list()

    for (val_idx in 1:2) {
      val_name <- apob_ds_names[val_idx]
      train_name <- apob_ds_names[-val_idx]

      train_df <- apob_datasets[[train_name]]

      # Impute
      train_df <- impute_simple(train_df, tier_vars)

      avail <- tier_vars[tier_vars %in% names(train_df)]
      cc <- complete.cases(train_df[, c(avail, "ascvd_combined")])
      train_cc <- train_df[cc, ]

      if (nrow(train_cc) < 50 || sum(train_cc$ascvd_combined) < 15) {
        cat(sprintf("    Fold %d (val=%s): SKIP\n", val_idx, val_name))
        next
      }

      X_train <- as.matrix(train_cc[, avail])
      y_train <- train_cc$ascvd_combined

      set.seed(2026 + val_idx + 100)
      cv_fit <- tryCatch(
        cv.glmnet(X_train, y_train, family = "binomial", alpha = alpha_val,
                  nfolds = min(10, floor(nrow(train_cc) / 10)), type.measure = "auc"),
        error = function(e) NULL
      )
      if (is.null(cv_fit)) next

      coefs <- as.numeric(coef(cv_fit, s = "lambda.min"))
      names(coefs) <- c("intercept", avail)

      val_df <- apob_datasets[[val_name]]
      val_res <- apply_frozen(val_df, coefs, avail)

      if (sum(val_res$valid) < 20) next

      val_auc_obj <- safe_auc(val_df$ascvd_combined[val_res$valid], val_res$pred[val_res$valid])
      if (is.na(val_auc_obj$auc)) next

      delta <- val_auc_obj$auc - apob_sh_aucs[val_name]
      fold_aucs <- c(fold_aucs, val_auc_obj$auc)
      fold_deltas <- c(fold_deltas, delta)

      cat(sprintf("    Fold %d: Train=%s (N=%d), Val=%s: AUC=%.4f (%.4f-%.4f) [d=%+.4f]\n",
                  val_idx, train_name, nrow(train_cc),
                  val_name, val_auc_obj$auc, val_auc_obj$lo, val_auc_obj$hi, delta))

      # Store per-fold predictions for DeLong test
      sh_pred_val2 <- apply_safeheart(val_df)
      both_ok2 <- val_res$valid & !is.na(sh_pred_val2)
      fold_key2 <- paste0(tier_name, "|", alpha_val)
      if (is.null(iecv_fold_preds[[fold_key2]])) iecv_fold_preds[[fold_key2]] <- list()
      iecv_fold_preds[[fold_key2]][[val_name]] <- list(
        new_pred = val_res$pred[both_ok2],
        sh_pred  = sh_pred_val2[both_ok2],
        outcome  = val_df$ascvd_combined[both_ok2],
        n        = sum(both_ok2)
      )

      fold_details[[val_name]] <- list(auc = val_auc_obj$auc, lo = val_auc_obj$lo,
                                        hi = val_auc_obj$hi, delta = delta,
                                        n = val_auc_obj$n, events = val_auc_obj$events)
    }

    if (length(fold_aucs) == 0) next

    weights <- sapply(fold_details, function(x) x$events)
    pooled_auc <- sum(sapply(fold_details, function(x) x$auc) * weights) / sum(weights)
    pooled_delta <- sum(sapply(fold_details, function(x) x$delta) * weights) / sum(weights)

    cat(sprintf("    >> Pooled IECV AUC=%.4f, Mean delta=%+.4f\n", pooled_auc, mean(fold_deltas)))

    row <- data.frame(
      tier = tier_name, alpha = alpha_val, n_vars = length(tier_vars),
      iecv_type = "2fold_apob",
      ukb_auc = if ("UKB" %in% names(fold_details)) fold_details$UKB$auc else NA,
      d3_auc = if ("Dragon3" %in% names(fold_details)) fold_details$Dragon3$auc else NA,
      wales_auc = NA,
      pooled_auc = pooled_auc, mean_auc = mean(fold_aucs),
      mean_delta = mean(fold_deltas), pooled_delta = pooled_delta,
      n_folds = length(fold_aucs),
      stringsAsFactors = FALSE
    )
    all_iecv_results <- rbind(all_iecv_results, row)
  }
}

# =============================================================================
# SECTION 6: DeLong TESTS — PER-FOLD (TRUE EXTERNAL VALIDATION)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 6: DeLong TESTS (per-fold IECV validation vs SAFEHEART)\n")
cat("===================================================================\n\n")

# Per-fold DeLong: uses stored IECV validation predictions (true external validation)
# Each fold's validation predictions come from a model that NEVER saw the validation data
if (nrow(all_iecv_results) > 0) {
  best_row <- all_iecv_results[which.max(all_iecv_results$pooled_delta), ]
  cat(sprintf("  Best IECV model: %s / alpha=%.2f, pooled delta=%+.4f\n",
              best_row$tier, best_row$alpha, best_row$pooled_delta))

  best_key <- paste0(best_row$tier, "|", best_row$alpha)

  if (best_key %in% names(iecv_fold_preds)) {
    fold_preds <- iecv_fold_preds[[best_key]]

    cat("\n  Per-fold DeLong tests (true external validation):\n")
    for (ds_name in names(fold_preds)) {
      fp <- fold_preds[[ds_name]]
      if (fp$n < 30) {
        cat(sprintf("    %s: SKIP — insufficient valid cases (%d)\n", ds_name, fp$n))
        next
      }

      roc_new <- tryCatch(pROC::roc(fp$outcome, fp$new_pred, quiet = TRUE), error = function(e) NULL)
      roc_sh  <- tryCatch(pROC::roc(fp$outcome, fp$sh_pred, quiet = TRUE), error = function(e) NULL)

      if (is.null(roc_new) || is.null(roc_sh)) {
        cat(sprintf("    %s: ROC computation failed (N=%d)\n", ds_name, fp$n))
        next
      }

      dl <- tryCatch(pROC::roc.test(roc_new, roc_sh, method = "delong"), error = function(e) NULL)
      new_auc <- as.numeric(pROC::auc(roc_new))
      sh_auc  <- as.numeric(pROC::auc(roc_sh))

      if (!is.null(dl)) {
        cat(sprintf("    %s: New=%.4f vs SH=%.4f, DeLong p=%.4f %s (N=%d)\n",
                    ds_name, new_auc, sh_auc, dl$p.value,
                    if (dl$p.value < 0.05) "*SIG*" else "(ns)", fp$n))
      } else {
        cat(sprintf("    %s: New=%.4f vs SH=%.4f, DeLong FAILED (N=%d)\n",
                    ds_name, new_auc, sh_auc, fp$n))
      }
    }
  } else {
    cat("  WARNING: No stored fold predictions for best model — skipping DeLong\n")
  }

  # -----------------------------------------------------------------------
  # Final model: fit on ALL datasets combined for coefficient export
  # NOTE: This is for reporting final coefficients only, NOT for validation
  # -----------------------------------------------------------------------
  if (best_row$tier %in% names(universal_tiers)) {
    best_vars <- universal_tiers[[best_row$tier]]
  } else {
    best_vars <- apob_tiers[[best_row$tier]]
  }

  cat("\n  Fitting final model on all data combined (for coefficient export)...\n")

  all_data_list <- list()
  for (ds_name in ds_names) {
    cols <- c(best_vars, "ascvd_combined")
    cols_avail <- cols[cols %in% names(datasets[[ds_name]])]
    all_data_list[[ds_name]] <- datasets[[ds_name]][, cols_avail, drop = FALSE]
  }
  common <- Reduce(intersect, lapply(all_data_list, names))
  for (n in names(all_data_list)) all_data_list[[n]] <- all_data_list[[n]][, common, drop = FALSE]
  all_combined <- do.call(rbind, all_data_list)
  all_combined <- impute_simple(all_combined, best_vars)

  avail <- best_vars[best_vars %in% names(all_combined)]
  cc <- complete.cases(all_combined[, c(avail, "ascvd_combined")])
  all_cc <- all_combined[cc, ]

  X_all <- as.matrix(all_cc[, avail])
  y_all <- all_cc$ascvd_combined

  set.seed(2026)
  final_fit <- cv.glmnet(X_all, y_all, family = "binomial", alpha = best_row$alpha,
                          nfolds = 10, type.measure = "auc")
  final_coefs <- as.numeric(coef(final_fit, s = "lambda.min"))
  names(final_coefs) <- c("intercept", avail)

  cat("\n  FINAL POOLED COEFFICIENTS (all data — for deployment):\n")
  for (v in names(final_coefs)) {
    if (v == "intercept") {
      cat(sprintf("    intercept = %.6f\n", final_coefs[v]))
    } else {
      cat(sprintf("    %-18s coef=% .6f  OR=%.4f\n", v, final_coefs[v], exp(final_coefs[v])))
    }
  }

  # Save final coefficients
  coef_df <- data.frame(variable = names(final_coefs), coefficient = final_coefs,
                         OR = exp(final_coefs), stringsAsFactors = FALSE)
  write.csv(coef_df, paste0(TAB_DIR, "calon2_s15_iecv_final_coefs.csv"), row.names = FALSE)
  cat("\n  Saved: calon2_s15_iecv_final_coefs.csv\n")
}

# =============================================================================
# SECTION 7: GRAND COMPARISON
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 7: GRAND COMPARISON TABLE\n")
cat("===================================================================\n\n")

# Add SAFEHEART baseline
sh_row <- data.frame(
  tier = "SAFEHEART-RE", alpha = NA, n_vars = 6,
  iecv_type = "baseline",
  ukb_auc = sh_ukb$auc, d3_auc = sh_d3$auc, wales_auc = sh_wales$auc,
  pooled_auc = NA, mean_auc = NA, mean_delta = 0, pooled_delta = 0, n_folds = NA,
  stringsAsFactors = FALSE
)

comparison <- rbind(sh_row, all_iecv_results)
comparison <- comparison[order(-comparison$pooled_delta), ]

cat(sprintf("  %-14s %5s %4s %-16s %7s %7s %7s %9s %+9s\n",
            "Tier", "Alpha", "Vars", "IECV_Type", "UKB", "D3", "Wales", "Pooled", "Delta"))
cat(paste(rep("-", 100), collapse = ""), "\n")

for (i in 1:nrow(comparison)) {
  r <- comparison[i, ]
  alpha_str <- if (is.na(r$alpha)) "  -" else sprintf("%.2f", r$alpha)
  ukb_str <- if (is.na(r$ukb_auc)) "  N/A" else sprintf("%.4f", r$ukb_auc)
  d3_str <- if (is.na(r$d3_auc)) "  N/A" else sprintf("%.4f", r$d3_auc)
  w_str <- if (is.na(r$wales_auc)) "  N/A" else sprintf("%.4f", r$wales_auc)
  pool_str <- if (is.na(r$pooled_auc)) "   N/A" else sprintf("%.4f", r$pooled_auc)
  delta_str <- if (is.na(r$pooled_delta)) "     N/A" else sprintf("%+.4f", r$pooled_delta)

  cat(sprintf("  %-14s %5s %4d %-16s %7s %7s %7s %9s %9s\n",
              r$tier, alpha_str, r$n_vars, r$iecv_type,
              ukb_str, d3_str, w_str, pool_str, delta_str))
}

# =============================================================================
# SECTION 8: SAVE RESULTS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 8: SAVE RESULTS\n")
cat("===================================================================\n\n")

write.csv(comparison, paste0(TAB_DIR, "calon2_s15_iecv_results.csv"), row.names = FALSE)
cat(sprintf("  Saved: calon2_s15_iecv_results.csv (%d rows)\n", nrow(comparison)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n======================================================================\n")
cat("  SCRIPT 15 COMPLETE: INTERNAL-EXTERNAL CROSS-VALIDATION             \n")
cat("======================================================================\n")
cat(sprintf("  UKB:     N=%d, Events=%d\n", nrow(ukb_h), sum(ukb_h$ascvd_combined)))
cat(sprintf("  Dragon3: N=%d, Events=%d\n", nrow(d3), sum(d3$ascvd_combined)))
cat(sprintf("  Wales:   N=%d, Events=%d\n", nrow(wales_h), sum(wales_h$ascvd_combined)))
cat(sprintf("  Total IECV configs: %d\n", nrow(all_iecv_results)))
cat(sprintf("  SAFEHEART baselines: UKB=%.4f, D3=%.4f, W=%.4f\n",
            sh_ukb$auc, sh_d3$auc, sh_wales$auc))
if (nrow(all_iecv_results) > 0) {
  best <- all_iecv_results[which.max(all_iecv_results$pooled_delta), ]
  cat(sprintf("  Best IECV: %s / a=%.2f, pooled_delta=%+.4f\n",
              best$tier, best$alpha, best$pooled_delta))
}
cat("======================================================================\n")

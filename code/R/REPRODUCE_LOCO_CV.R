#!/usr/bin/env Rscript
################################################################################
#
#  CALON-FH: STANDALONE LOCO-CV REPRODUCTION SCRIPT
#
#  Purpose: Reproduce the Leave-One-Cohort-Out cross-validation AUCs
#           for the CALON ASCVD prognostic model in genetically confirmed FH
#
#  Outcome: ASCVD (prevalent + incident cardiovascular disease)
#  Models:  CALON-Lite (8 features), CALON-Lite-noLpa (5 features)
#  Method:  Elastic Net (alpha=0.3), MICE (m=5), Rubin's rules pooling
#
#  Cohorts:
#    1. South Wales FH — mutation positive (DRAGON3, n~332-418)
#    2. UK Biobank — WES-confirmed FH (n~1,623)
#    3. Wales — PASS All Wales Registry, Positive1==1 (n~2,405, deduplicated)
#
#  Author: Dr Nader Genedy
#  Date: 2026-03-26
#  Reviewer: Claude (biostatistics, lipid medicine, molecular genetics)
#
################################################################################

cat("\n")
cat("================================================================\n")
cat("  CALON-FH: LOCO-CV REPRODUCTION\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

set.seed(42)

# ── 0. PACKAGES ──────────────────────────────────────────────────────────────

pkgs <- c("tidyverse", "readxl", "haven", "pROC", "mice", "glmnet", "boot")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(haven)
  library(pROC)
  library(mice)
  library(glmnet)
  library(boot)
})

# ── 1. PATHS ─────────────────────────────────────────────────────────────────

data_dir <- "C:/Users/nader/Downloads/calon_ukb_pipeline"

# Data files
DRAGON3_FILE  <- file.path(data_dir, "DRAGON_3.csv")
WALES_FILE    <- file.path(data_dir, "WALES_FH_CLEANED (1) - Copy.csv")
UKB_FEAT_FILE <- "D:/calon_backup_data/TUDOR_UKB_Features (1).csv"
UKB_PROC_FILE <- file.path(data_dir, "calon_ukb_analysis_ready.csv")
# Fallback: check output/ subdir
if (!file.exists(UKB_PROC_FILE))
  UKB_PROC_FILE <- file.path(data_dir, "output", "calon_ukb_analysis_ready.csv")

# Verify all exist
for (f in c(DRAGON3_FILE, WALES_FILE, UKB_FEAT_FILE, UKB_PROC_FILE)) {
  if (!file.exists(f)) stop(sprintf("MISSING: %s", f))
  cat(sprintf("  OK: %s (%.1f MB)\n", basename(f), file.info(f)$size / 1e6))
}

# ── 2. HELPER FUNCTIONS ─────────────────────────────────────────────────────

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

parse_gene <- function(mut) {
  ifelse(is.na(mut), NA_character_,
         ifelse(grepl("^LDLR", mut), "LDLR",
                ifelse(grepl("^APOB", mut), "APOB",
                       ifelse(grepl("^PCSK9", mut), "PCSK9", "Other"))))
}

pool_rubins <- function(estimates, ses) {
  m <- length(estimates)
  valid <- !is.na(estimates)
  if (sum(valid) == 0) return(list(est = NA, se = NA, ci_low = NA, ci_high = NA))
  Q_bar <- mean(estimates[valid])
  var_m <- ses[valid]^2
  W <- mean(var_m)
  B <- var(estimates[valid])
  T_var <- W + (1 + 1/sum(valid)) * B
  pooled_se <- sqrt(max(T_var, 0))
  list(est = Q_bar, se = pooled_se,
       ci_low = Q_bar - 1.96 * pooled_se,
       ci_high = Q_bar + 1.96 * pooled_se)
}

get_auc_ci <- function(train_df, test_df, features, outcome = "ascvd", n_boot = 500) {
  cols <- c(features, outcome)
  train <- train_df %>% select(all_of(cols)) %>% filter(complete.cases(.))
  test  <- test_df  %>% select(all_of(cols)) %>% filter(complete.cases(.))
  if (nrow(train) < 30 || nrow(test) < 30 || length(unique(test[[outcome]])) < 2)
    return(list(auc = NA, ci_low = NA, ci_high = NA, n_train = nrow(train), n_test = nrow(test)))
  X_train <- as.matrix(train[, features])
  y_train <- train[[outcome]]
  X_test  <- as.matrix(test[, features])
  y_test  <- test[[outcome]]
  cv_fit <- cv.glmnet(X_train, y_train, family = "binomial",
                       alpha = 0.3, nfolds = min(10, nrow(train)),
                       type.measure = "auc", standardize = TRUE)
  pred    <- as.numeric(predict(cv_fit, newx = X_test, s = "lambda.min", type = "response"))
  roc_obj <- roc(y_test, pred, quiet = TRUE)
  ci_obj  <- ci.auc(roc_obj, method = "bootstrap", boot.n = n_boot, quiet = TRUE)
  list(auc = as.numeric(auc(roc_obj)), ci_low = ci_obj[1], ci_high = ci_obj[3],
       n_train = nrow(train), n_test = nrow(test))
}

# ── 3. LOAD & HARMONISE: SOUTH WALES FH ─────────────────────────────────────

cat("\n=== Loading South Wales FH (DRAGON3) ===\n")
dragon_raw <- read.csv(DRAGON3_FILE, stringsAsFactors = FALSE, check.names = FALSE)
names(dragon_raw) <- trimws(names(dragon_raw))
cat(sprintf("  Raw: %d rows x %d cols\n", nrow(dragon_raw), ncol(dragon_raw)))

# Filter to genetically confirmed (Positive1 == 1)
pos1 <- safe_num(dragon_raw$Positive1)
dragon <- dragon_raw[which(pos1 == 1), ]
cat(sprintf("  After Positive1 filter: %d rows\n", nrow(dragon)))

# Harmonise columns
d <- data.frame(
  cohort       = "South Wales FH",
  id           = dragon$DatabaseNumber,
  age          = safe_num(dragon$Ageattest),
  sex          = ifelse(dragon$Gender == "M", 1, ifelse(dragon$Gender == "F", 0, NA)),
  bmi          = safe_num(dragon$BMI),
  ldl          = safe_num(dragon$LDL_1),
  hdl          = coalesce(safe_num(dragon$HDL_1), safe_num(dragon$LastHDL)),
  tg           = coalesce(safe_num(dragon$TRG_1), safe_num(dragon$LastTrigs)),
  apob         = safe_num(dragon$ApoB),
  lpa_nmol     = safe_num(dragon$Lpa),
  smoking_raw  = safe_num(dragon$Smoking_binary),
  diabetes     = safe_num(dragon$Diabetes_binary),
  bp_med       = safe_num(dragon$BloodPressureMedication),
  bp_sys       = safe_num(dragon$BloodPressureSystolic),
  bp_dia       = safe_num(dragon$BloodPressureDiastolic),
  ascvd        = safe_num(dragon$ASCVD_combined),
  gene         = parse_gene(dragon$Mutation1),
  stringsAsFactors = FALSE
)

d$inv_hdl        <- ifelse(!is.na(d$hdl) & d$hdl > 0, 1 / d$hdl, NA)
d$smoking_binary <- ifelse(d$smoking_raw %in% c(1, 2), 1, ifelse(d$smoking_raw == 0, 0, NA))
d$hypertension   <- ifelse(d$bp_med == 1 | d$bp_sys > 140 | d$bp_dia > 90, 1, 0)
d$lpa_143        <- ifelse(!is.na(d$lpa_nmol) & d$lpa_nmol > 143, 1, 0)
d$lipid_years    <- ifelse(!is.na(d$ldl) & !is.na(d$age) & d$age > 0, d$ldl * d$age, NA)
d$ever_smoked    <- d$smoking_binary

cat(sprintf("  South Wales FH: %d patients, ASCVD+ = %d (%.1f%%)\n",
            nrow(d), sum(d$ascvd == 1, na.rm = TRUE),
            100 * mean(d$ascvd, na.rm = TRUE)))

# ── 4. LOAD & HARMONISE: WALES (PASS) ───────────────────────────────────────

cat("\n=== Loading Wales FH (PASS Registry) ===\n")
wales_raw <- read.csv(WALES_FILE, stringsAsFactors = FALSE, check.names = FALSE)
names(wales_raw) <- trimws(names(wales_raw))
cat(sprintf("  Raw: %d rows x %d cols\n", nrow(wales_raw), ncol(wales_raw)))

# Filter to genetically confirmed
pos1_w <- safe_num(wales_raw$Positive1)
wales <- wales_raw[which(pos1_w == 1), ]
cat(sprintf("  After Positive1 filter: %d rows\n", nrow(wales)))

w <- data.frame(
  cohort       = "Wales",
  id           = wales$DatabaseNumber,
  age          = safe_num(wales$BMI_AGE),
  sex          = ifelse(wales$Gender == "M", 1, ifelse(wales$Gender == "F", 0, NA)),
  bmi          = safe_num(wales$BMI),
  ldl          = safe_num(wales$LDL.1),
  hdl          = coalesce(safe_num(wales$HDL.1), safe_num(wales$HDL.2)),
  tg           = coalesce(safe_num(wales$TRG.1), safe_num(wales$TRG.2)),
  apob         = NA_real_,
  lpa_nmol     = safe_num(wales$Lpa.1),
  smoking_raw  = wales$Smoking,
  diabetes_raw = wales$Diabetes,
  bp_med_raw   = wales$BloodPressureMedication,
  bp_sys       = safe_num(wales$BloodPressureSystolic),
  bp_dia       = safe_num(wales$BloodPressureDiastolic),
  on_treat     = safe_num(wales$OnTreatment),
  ascvd        = safe_num(wales$ascvd_combine),
  gene         = parse_gene(wales$Mutation1),
  stringsAsFactors = FALSE
)

w$inv_hdl        <- ifelse(!is.na(w$hdl) & w$hdl > 0, 1 / w$hdl, NA)
w$smoking_binary <- ifelse(w$smoking_raw %in% c("1", "True", "TRUE"), 1,
                           ifelse(w$smoking_raw %in% c("0", "False", "FALSE"), 0, NA))
w$ever_smoked    <- w$smoking_binary
w$diabetes       <- ifelse(w$diabetes_raw %in% c("1", "True"), 1,
                           ifelse(w$diabetes_raw %in% c("0", "False"), 0, NA))
w$hypertension   <- ifelse(w$bp_med_raw %in% c("1", "True") |
                             (!is.na(w$bp_sys) & w$bp_sys > 140) |
                             (!is.na(w$bp_dia) & w$bp_dia > 90), 1, 0)
w$lpa_143        <- ifelse(!is.na(w$lpa_nmol) & w$lpa_nmol > 143, 1, 0)
w$lipid_years    <- ifelse(!is.na(w$ldl) & !is.na(w$age) & w$age > 0, w$ldl * w$age, NA)

cat(sprintf("  Wales FH: %d patients, ASCVD+ = %d (%.1f%%)\n",
            nrow(w), sum(w$ascvd == 1, na.rm = TRUE),
            100 * mean(w$ascvd, na.rm = TRUE)))

# ── 5. LOAD & HARMONISE: UK BIOBANK ─────────────────────────────────────────

cat("\n=== Loading UK Biobank ===\n")
ukb_raw <- read.csv(UKB_FEAT_FILE, stringsAsFactors = FALSE, check.names = FALSE)
names(ukb_raw) <- make.names(names(ukb_raw), unique = TRUE)
cat(sprintf("  Raw: %d rows x %d cols\n", nrow(ukb_raw), ncol(ukb_raw)))

# Filter to genetically confirmed FH
ukb_fh <- ukb_raw[which(ukb_raw$is_fh_genetic == 1), ]
cat(sprintf("  Genetically confirmed FH: %d\n", nrow(ukb_fh)))

# Load processed UKB for Lp(a), ASCVD, smoking
ukb_proc <- read.csv(UKB_PROC_FILE, stringsAsFactors = FALSE, check.names = FALSE)
names(ukb_proc) <- make.names(names(ukb_proc), unique = TRUE)
cat(sprintf("  Processed UKB: %d rows x %d cols\n", nrow(ukb_proc), ncol(ukb_proc)))

# Merge processed fields into ukb_fh
proc_fields <- c("eid", "lpa", "ever_smoked", "smoking_binary",
                  "ascvd_combined", "log_apob_ldl", "hypertension",
                  "diabetes", "re_ldl")
proc_fields <- proc_fields[proc_fields %in% names(ukb_proc)]
ukb_proc_sub <- ukb_proc[, proc_fields, drop = FALSE]
ukb_proc_sub <- ukb_proc_sub[!duplicated(ukb_proc_sub$eid), ]

ukb_fh <- merge(ukb_fh, ukb_proc_sub, by.x = "participant.eid", by.y = "eid",
                 all.x = TRUE, sort = FALSE)

# Statin detection
statin_codes <- c(1140888648, 1140888594, 1140888690, 1141146234,
                  1141192414, 1140861970, 1140881748, 1141146138)
med_cols <- grep("^participant\\.p20003_i0", names(ukb_fh), value = TRUE)
ukb_fh$on_statin <- apply(ukb_fh[, med_cols], 1, function(row) {
  as.integer(any(row %in% statin_codes, na.rm = TRUE))
})

# Self-reported conditions
ill_cols <- grep("^participant\\.p20002_i0", names(ukb_fh), value = TRUE)
if (!"diabetes" %in% names(ukb_fh) || all(is.na(ukb_fh$diabetes))) {
  ukb_fh$diabetes <- apply(ukb_fh[, ill_cols], 1, function(row) {
    as.integer(any(row %in% c(1220, 1222, 1223), na.rm = TRUE))
  })
}
if (!"hypertension" %in% names(ukb_fh) || all(is.na(ukb_fh$hypertension))) {
  ukb_fh$hypertension <- apply(ukb_fh[, ill_cols], 1, function(row) {
    as.integer(any(row %in% c(1065, 1072), na.rm = TRUE))
  })
}

# Harmonise
u <- data.frame(
  cohort       = "UKB",
  id           = as.character(ukb_fh$participant.eid),
  age          = safe_num(ukb_fh$participant.p21022),
  sex          = safe_num(ukb_fh$participant.p31),
  bmi          = safe_num(ukb_fh$participant.p21001_i0),
  ldl          = safe_num(ukb_fh$participant.p30780_i0),
  hdl          = safe_num(ukb_fh$participant.p30760_i0),
  tg           = safe_num(ukb_fh$participant.p30870_i0),
  apob         = safe_num(ukb_fh$participant.p30640_i0),
  lpa_nmol     = safe_num(ukb_fh$lpa),
  smoking_binary = safe_num(ukb_fh$smoking_binary),
  ever_smoked  = safe_num(ukb_fh$ever_smoked),
  diabetes     = safe_num(ukb_fh$diabetes),
  hypertension = safe_num(ukb_fh$hypertension),
  ascvd        = safe_num(ukb_fh$ascvd_combined),
  gene         = ukb_fh$gene,
  stringsAsFactors = FALSE
)

u$inv_hdl     <- ifelse(!is.na(u$hdl) & u$hdl > 0, 1 / u$hdl, NA)
u$lpa_143     <- ifelse(!is.na(u$lpa_nmol) & u$lpa_nmol > 143, 1, 0)
u$lipid_years <- ifelse(!is.na(u$ldl) & !is.na(u$age) & u$age > 0, u$ldl * u$age, NA)

# Deduplicate
u <- u[!duplicated(u$id), ]

cat(sprintf("  UKB FH: %d patients, ASCVD+ = %d (%.1f%%)\n",
            nrow(u), sum(u$ascvd == 1, na.rm = TRUE),
            100 * mean(u$ascvd, na.rm = TRUE)))
cat(sprintf("  UKB Lp(a) available: %d (%.1f%%)\n",
            sum(!is.na(u$lpa_nmol)), 100 * mean(!is.na(u$lpa_nmol))))

# ── 6. COMBINE & DEDUPLICATE ────────────────────────────────────────────────

cat("\n=== COMBINING COHORTS ===\n")

common <- c("cohort", "id", "age", "sex", "bmi", "ldl", "hdl", "tg", "apob",
            "lpa_nmol", "inv_hdl", "smoking_binary", "ever_smoked",
            "diabetes", "hypertension", "lpa_143", "ascvd", "gene", "lipid_years")

for (cc in common) {
  if (!cc %in% names(d)) d[[cc]] <- NA
  if (!cc %in% names(w)) w[[cc]] <- NA
  if (!cc %in% names(u)) u[[cc]] <- NA
}

fh <- bind_rows(d[, common], w[, common], u[, common])
fh$cohort <- factor(fh$cohort, levels = c("South Wales FH", "UKB", "Wales"))

cat(sprintf("  Before dedup: %d (SW=%d, UKB=%d, Wales=%d)\n",
            nrow(fh),
            sum(fh$cohort == "South Wales FH"),
            sum(fh$cohort == "UKB"),
            sum(fh$cohort == "Wales")))

# CRITICAL: South Wales is a SUBSET of Wales — remove duplicates
sw_ids <- fh$id[fh$cohort == "South Wales FH"]
dup_mask <- fh$cohort == "Wales" & fh$id %in% sw_ids
n_dups <- sum(dup_mask)
fh <- fh[!dup_mask, ]

cat(sprintf("  Removed %d SW patients from Wales (prevent data leakage)\n", n_dups))
cat(sprintf("  After dedup: %d (SW=%d, UKB=%d, Wales=%d)\n",
            nrow(fh),
            sum(fh$cohort == "South Wales FH"),
            sum(fh$cohort == "UKB"),
            sum(fh$cohort == "Wales")))

# Clean impossible values
fh$tg[!is.na(fh$tg) & fh$tg <= 0]   <- NA
fh$hdl[!is.na(fh$hdl) & fh$hdl <= 0] <- NA

# ── 7. MISSINGNESS DIAGNOSTICS ──────────────────────────────────────────────

cat("\n=== MISSINGNESS BY COHORT ===\n")
lite_vars <- c("age", "sex", "ldl", "inv_hdl", "smoking_binary",
               "diabetes", "hypertension", "lpa_143")
nolpa_vars <- c("lipid_years", "sex", "diabetes", "ever_smoked", "inv_hdl")

for (coh in levels(fh$cohort)) {
  sub <- fh[fh$cohort == coh, ]
  cat(sprintf("\n  %s (N=%d, ASCVD+ = %d [%.1f%%]):\n",
              coh, nrow(sub), sum(sub$ascvd == 1, na.rm = TRUE),
              100 * mean(sub$ascvd, na.rm = TRUE)))
  for (v in lite_vars) {
    n_miss <- sum(is.na(sub[[v]]))
    if (n_miss > 0)
      cat(sprintf("    %-16s: %d/%d missing (%.1f%%)\n", v, n_miss, nrow(sub), 100*n_miss/nrow(sub)))
  }
  cat(sprintf("    Complete for CALON-Lite: %d (%.1f%%)\n",
              sum(complete.cases(sub[, lite_vars])),
              100 * mean(complete.cases(sub[, lite_vars]))))
}

# ── 8. MULTIPLE IMPUTATION (MICE) ───────────────────────────────────────────

cat("\n=== MULTIPLE IMPUTATION (m=5, maxit=10, by cohort) ===\n")

m_imp <- 5
imp_vars <- c("age", "sex", "ldl", "hdl", "inv_hdl", "tg",
              "smoking_binary", "ever_smoked", "diabetes", "hypertension",
              "lpa_143", "bmi", "ascvd", "lipid_years")

fh_imputed <- vector("list", m_imp)

for (i in 1:m_imp) {
  cat(sprintf("  Imputation %d/%d...\n", i, m_imp))
  fh_imp_i <- fh

  for (coh in levels(fh$cohort)) {
    idx <- which(fh_imp_i$cohort == coh)
    sub <- fh_imp_i[idx, imp_vars]

    # Skip if nothing to impute
    if (sum(!complete.cases(sub)) == 0) next
    if (nrow(sub) < 50) next

    # Set method: logistic for binary, pmm for continuous
    methods <- rep("pmm", ncol(sub))
    names(methods) <- names(sub)
    for (v in c("sex", "smoking_binary", "ever_smoked", "diabetes",
                "hypertension", "lpa_143", "ascvd")) {
      if (v %in% names(sub)) methods[v] <- "logreg"
    }
    # Don't impute variables with 0% available
    for (v in names(sub)) {
      if (all(is.na(sub[[v]]))) methods[v] <- ""
    }

    tryCatch({
      mice_obj <- mice(sub, m = 1, maxit = 10, method = methods,
                        seed = 42 + i * 100, printFlag = FALSE)
      imp_data <- complete(mice_obj, 1)
      fh_imp_i[idx, imp_vars] <- imp_data
    }, error = function(e) {
      cat(sprintf("    MICE error in %s: %s\n", coh, e$message))
    })
  }

  # Recompute derived variables
  fh_imp_i$inv_hdl <- ifelse(!is.na(fh_imp_i$hdl) & fh_imp_i$hdl > 0,
                              1 / fh_imp_i$hdl, fh_imp_i$inv_hdl)
  fh_imp_i$lipid_years <- ifelse(!is.na(fh_imp_i$ldl) & !is.na(fh_imp_i$age),
                                  fh_imp_i$ldl * fh_imp_i$age, fh_imp_i$lipid_years)

  fh_imputed[[i]] <- fh_imp_i
}

cat("  MICE complete.\n")

# ── 9. LOCO-CV ──────────────────────────────────────────────────────────────

cat("\n=== LOCO-CV (Leave-One-Cohort-Out) ===\n")
cat("  Training on 2 cohorts, testing on holdout\n")
cat("  Pooling across m=5 imputations via Rubin's rules\n\n")

LITE_features  <- c("age", "sex", "ldl", "inv_hdl", "smoking_binary",
                     "diabetes", "hypertension", "lpa_143")
NOLPA_features <- c("lipid_years", "sex", "diabetes", "ever_smoked", "inv_hdl")

model_defs <- list(
  "CALON-Lite-noLpa" = NOLPA_features,
  "CALON-Lite"       = LITE_features
)

loco_results <- data.frame()

for (holdout in c("South Wales FH", "UKB", "Wales")) {
  for (mname in names(model_defs)) {
    feats <- model_defs[[mname]]
    aucs <- ses <- n_tests <- numeric(m_imp)
    any_valid <- FALSE

    for (i in 1:m_imp) {
      train_loco <- fh_imputed[[i]] %>% filter(cohort != holdout)
      test_loco  <- fh_imputed[[i]] %>% filter(cohort == holdout)

      res <- get_auc_ci(train_loco, test_loco, feats, outcome = "ascvd", n_boot = 500)
      aucs[i]    <- res$auc
      n_tests[i] <- res$n_test

      if (!is.na(res$auc)) {
        ses[i] <- (res$ci_high - res$ci_low) / (2 * 1.96)
        any_valid <- TRUE
      } else {
        ses[i] <- NA
      }
    }

    if (any_valid) {
      pooled <- pool_rubins(aucs, ses)
      cat(sprintf("  LOCO holdout=%-15s %s: AUC=%.4f [%.4f-%.4f] (N_test~%d)\n",
                  holdout, mname, pooled$est, pooled$ci_low, pooled$ci_high,
                  round(mean(n_tests))))
    } else {
      pooled <- list(est = NA, ci_low = NA, ci_high = NA)
      cat(sprintf("  LOCO holdout=%-15s %s: AUC=NA\n", holdout, mname))
    }

    loco_results <- rbind(loco_results, data.frame(
      holdout = holdout, model = mname,
      auc = pooled$est, ci_low = pooled$ci_low, ci_high = pooled$ci_high,
      n_test = round(mean(n_tests)), stringsAsFactors = FALSE
    ))
  }
}

# ── 10. COMPARE WITH EXISTING RESULTS ───────────────────────────────────────

cat("\n=== COMPARISON WITH EXISTING RESULTS ===\n")

existing_file <- file.path(data_dir, "LOCO_CV_results.csv")
if (file.exists(existing_file)) {
  existing <- read.csv(existing_file, stringsAsFactors = FALSE)
  cat("\n  EXISTING results:\n")
  for (i in seq_len(nrow(existing))) {
    cat(sprintf("    %-15s %-20s AUC=%.4f [%.4f-%.4f]\n",
                existing$holdout[i], existing$model[i],
                existing$auc[i], existing$ci_low[i], existing$ci_high[i]))
  }
}

cat("\n  REPRODUCED results:\n")
for (i in seq_len(nrow(loco_results))) {
  cat(sprintf("    %-15s %-20s AUC=%.4f [%.4f-%.4f]\n",
              loco_results$holdout[i], loco_results$model[i],
              loco_results$auc[i], loco_results$ci_low[i], loco_results$ci_high[i]))
}

# Compute deltas
if (file.exists(existing_file)) {
  cat("\n  DELTAS (reproduced - existing):\n")
  merged <- merge(loco_results, existing,
                   by = c("holdout", "model"), suffixes = c("_new", "_old"))
  for (i in seq_len(nrow(merged))) {
    delta <- merged$auc_new[i] - merged$auc_old[i]
    cat(sprintf("    %-15s %-20s delta=%.4f (%s)\n",
                merged$holdout[i], merged$model[i], delta,
                ifelse(abs(delta) < 0.01, "CONCORDANT", "DIFFERS")))
  }
}

# ── 11. SAVE ─────────────────────────────────────────────────────────────────

out_file <- file.path(data_dir, "LOCO_CV_reproduced.csv")
write.csv(loco_results, out_file, row.names = FALSE)
cat(sprintf("\n  Saved: %s\n", out_file))

cat("\n================================================================\n")
cat("  LOCO-CV REPRODUCTION COMPLETE\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")

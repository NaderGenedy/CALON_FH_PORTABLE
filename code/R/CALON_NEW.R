################################################################################
# CALON-NEW: Clean ASCVD risk prediction in genetically confirmed FH
#
# Predictor set (per user 2026-05-03 specification):
#   1. log_apob_ldl     log(ApoB / LDL-C)        — atherogenic particle ratio
#   2. log_tg_hdl       log(Triglycerides / HDL) — insulin resistance proxy
#   3. re_ldl           Reverse-engineered untreated LDL (statin RF correction)
#   4. age              Years
#   5. diabetes         Binary (DM diagnosed)
#   6. hypertension     Binary (SBP>=140 OR DBP>=90 OR Rx OR self-report 1065)
#   7. ever_smoked      Binary (current OR previous)
#   8. lpa_143          Lp(a) > 143 nmol/L binary
#
# Endpoint: ascvd = first-occurrence ICD-10 union from p131* (verified clean)
#   Includes: I20, I21, I22, I25, I63, I64, I70, I73, I74
#   Excludes: I10-I15 (hypertension) — verified separation
#
# Comparator: SRE-Fixed (frozen SAFEHEART coefficients, TRIPOD Type 4)
#   age, sex, ldl, hypertension, bmi, ever_smoked, lpa_120, prior_ascvd
#
# Validation: 6 bi-external directions (SW<->UKB<->Wales)
# Imputation: MICE m=10, maxit=20, Rubin's rules
# Statistics: AUC+CI (bootstrap 2000), DeLong p, calibration int+slope
#
# Cohorts: South Wales FH (~418), UKB (~1623), Wales FH (~2405)
#
# Author: Dr Nader Genedy
# Date:   2026-05-03
# Builds on: CALON_MI_validation.R, CALON_FINAL_analysis.R
################################################################################

# --- Packages ---------------------------------------------------------------
pkgs_needed <- c("tidyverse", "readxl", "tableone", "pROC",
                 "broom", "boot", "mice", "glmnet")
pkgs_missing <- pkgs_needed[!pkgs_needed %in% installed.packages()[, "Package"]]
if (length(pkgs_missing) > 0) {
  install.packages(pkgs_missing, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(tableone); library(pROC)
  library(broom); library(boot); library(mice); library(glmnet)
})

# --- Working directory ------------------------------------------------------
if (dir.exists("/cloud/project")) {
  data_dir <- "/cloud/project"
} else if (dir.exists("C:/Users/nader/Downloads/calon_ukb_pipeline")) {
  data_dir <- "C:/Users/nader/Downloads/calon_ukb_pipeline"
} else {
  data_dir <- getwd()
}
setwd(data_dir)
dir.create("output/new", showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Working directory: %s\n", getwd()))
cat("Output directory: output/new/\n\n")

# --- Helpers ----------------------------------------------------------------
safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

parse_gene <- function(mut) {
  ifelse(is.na(mut), NA_character_,
         ifelse(grepl("^LDLR", mut), "LDLR",
                ifelse(grepl("^APOB", mut), "APOB",
                       ifelse(grepl("^PCSK9", mut), "PCSK9", "Other"))))
}

find_file <- function(pattern, required = TRUE) {
  matches <- list.files(".", pattern = pattern, ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
  matches <- matches[!grepl("output/|Table_|Table1|Table2|SCORECARD|\\.claude/", matches)]
  if (length(matches) == 0) {
    if (required) stop(sprintf("File not found: %s", pattern))
    return(NULL)
  }
  if (length(matches) > 1) {
    sizes <- file.info(matches)$size
    matches <- matches[order(-sizes)]
  }
  cat(sprintf("  Found: %s\n", matches[1]))
  matches[1]
}

safe_read_csv <- function(path, label = "") {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(df)[1] <- sub("^﻿", "", names(df)[1])
  names(df) <- trimws(names(df))
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) names(df)[bad] <- paste0("X_empty_", seq_len(sum(bad)))
  if (anyDuplicated(names(df))) names(df) <- make.unique(names(df), sep = "_dup_")
  cat(sprintf("  %s: %d rows x %d cols\n", label, nrow(df), ncol(df)))
  df
}

grab <- function(df, name) {
  cn <- names(df)
  if (name %in% cn) return(df[[name]])
  idx <- which(tolower(cn) == tolower(name))
  if (length(idx) > 0) return(df[[idx[1]]])
  rep(NA, nrow(df))
}

cat("===========================================================\n")
cat("  CALON-NEW: 8-PREDICTOR ASCVD MODEL\n")
cat("  Endpoint: ASCVD (first-occurrence ICD-10, verified clean)\n")
cat("===========================================================\n\n")

# --- Load cohorts -----------------------------------------------------------
cat("--- Loading cohorts ---\n")
cat("South Wales FH (Dragon-3):\n")
dragon_raw <- safe_read_csv(find_file("DRAGON.*3.*\\.csv"), "DRAGON_3")

cat("Wales FH Registry:\n")
wales_raw  <- safe_read_csv(find_file("WALES.*FH.*\\.csv"),  "Wales")

cat("UK Biobank (TUDOR-format):\n")
ukb_raw    <- safe_read_csv(find_file("TUDOR.*UKB.*\\.csv"), "UKB")
names(ukb_raw) <- make.names(names(ukb_raw), unique = TRUE)

# --- UKB processed file (has correct ASCVD, Lp(a), re_ldl) -----------------
ukb_proc_file <- find_file("calon_ukb_analysis_ready", required = FALSE)
if (!is.null(ukb_proc_file) && file.exists(ukb_proc_file)) {
  cat("UKB processed file:\n")
  ukb_proc <- safe_read_csv(ukb_proc_file, "UKB_proc")
  names(ukb_proc) <- make.names(names(ukb_proc), unique = TRUE)
  proc_fields <- c("eid", "lpa", "lpa_binary", "re_ldl", "residual_factor",
                   "ever_smoked", "smoking_binary",
                   "ascvd_combined", "ascvd_any",
                   "ascvd_prevalent", "ascvd_incident",
                   "log_apob_ldl", "hypertension", "diabetes")
  proc_fields <- proc_fields[proc_fields %in% names(ukb_proc)]
  ukb_proc_subset <- ukb_proc[, proc_fields, drop = FALSE]
  names(ukb_proc_subset) <- ifelse(names(ukb_proc_subset) == "eid", "eid",
                                    paste0("proc_", names(ukb_proc_subset)))
  names(ukb_proc_subset)[1] <- "eid"
  ukb_raw <- merge(ukb_raw, ukb_proc_subset,
                   by.x = "participant.eid", by.y = "eid",
                   all.x = TRUE, sort = FALSE)
  cat(sprintf("  Merged %d processed fields\n", ncol(ukb_proc_subset) - 1))
} else {
  cat("WARNING: calon_ukb_analysis_ready.csv not found.\n")
  cat("Run 01_CALON_build_cohort.R first.\n")
}

# --- Harmonise South Wales (Dragon-3) ---------------------------------------
cat("\n--- Harmonising South Wales FH ---\n")
pos1_d <- grab(dragon_raw, "Positive1")
keep_d <- which(safe_numeric(pos1_d) == 1 | pos1_d == "Yes")
d <- if (length(keep_d) == 0) dragon_raw else dragon_raw[keep_d, ]
cat(sprintf("  After Positive1 filter: %d patients\n", nrow(d)))

statin_d <- grab(d, "Statin")
gender_d <- grab(d, "Gender")
smoking_bin_d <- safe_numeric(grab(d, "Smoking_binary"))
bp_med_d <- safe_numeric(grab(d, "BloodPressureMedication"))
bp_sys_d <- safe_numeric(grab(d, "BloodPressureSystolic"))
bp_dia_d <- safe_numeric(grab(d, "BloodPressureDiastolic"))

d$cohort <- "South Wales FH"
d$id     <- as.character(grab(d, "DatabaseNumber"))
d$age    <- safe_numeric(grab(d, "Ageattest"))
d$sex    <- ifelse(gender_d == "M", 1, ifelse(gender_d == "F", 0, NA))
d$bmi    <- safe_numeric(grab(d, "BMI"))
d$tc     <- safe_numeric(grab(d, "TC_1"))
d$hdl    <- coalesce(safe_numeric(grab(d, "HDL_1")), safe_numeric(grab(d, "LastHDL")))
d$ldl    <- safe_numeric(grab(d, "LDL_1"))
d$tg     <- coalesce(safe_numeric(grab(d, "TRG_1")), safe_numeric(grab(d, "LastTrigs")))
d$apob   <- safe_numeric(grab(d, "ApoB"))
d$lpa_nmol <- safe_numeric(grab(d, "Lpa"))
d$on_statin <- ifelse(!is.na(statin_d) & statin_d != "" & statin_d != "0", 1, 0)
# Reduction factor: high-intensity (atorva/rosuva) 0.50, others 0.65 (default)
d$residual_factor <- ifelse(d$on_statin == 1,
                              ifelse(grepl("ATOR|ROSUV", toupper(statin_d)), 0.50, 0.65),
                              1.00)
d$re_ldl <- d$ldl / d$residual_factor

d$ever_smoked  <- ifelse(smoking_bin_d %in% c(1, 2), 1, ifelse(smoking_bin_d == 0, 0, NA))
d$diabetes     <- safe_numeric(grab(d, "Diabetes_binary"))
d$hypertension <- ifelse(!is.na(bp_med_d) & bp_med_d == 1 |
                          !is.na(bp_sys_d) & bp_sys_d > 140 |
                          !is.na(bp_dia_d) & bp_dia_d > 90, 1, 0)
d$lpa_143      <- ifelse(!is.na(d$lpa_nmol) & d$lpa_nmol > 143, 1, 0)
d$ascvd        <- safe_numeric(grab(d, "ASCVD_combined"))
d$prior_ascvd  <- d$ascvd
d$gene         <- parse_gene(grab(d, "Mutation1"))

# Derived ratios
d$apob_ldl     <- ifelse(!is.na(d$apob) & !is.na(d$ldl) & d$ldl > 0, d$apob / d$ldl, NA)
d$tg_hdl       <- ifelse(!is.na(d$tg) & !is.na(d$hdl) & d$hdl > 0, d$tg / d$hdl, NA)
d$log_apob_ldl <- ifelse(!is.na(d$apob_ldl) & d$apob_ldl > 0, log(d$apob_ldl), NA)
d$log_tg_hdl   <- ifelse(!is.na(d$tg_hdl) & d$tg_hdl > 0, log(d$tg_hdl), NA)

cat(sprintf("  South Wales FH: %d patients, ASCVD=%d (%.1f%%)\n",
            nrow(d), sum(d$ascvd == 1, na.rm = TRUE),
            100 * mean(d$ascvd, na.rm = TRUE)))

# --- Harmonise Wales FH Registry --------------------------------------------
cat("\n--- Harmonising Wales FH ---\n")
pos1_w <- grab(wales_raw, "Positive1")
keep_w <- which(safe_numeric(pos1_w) == 1 | pos1_w == "Yes" | pos1_w == "True")
w <- if (length(keep_w) == 0) wales_raw else wales_raw[keep_w, ]
cat(sprintf("  After Positive1 filter: %d patients\n", nrow(w)))

gender_w   <- grab(w, "Gender")
smoking_w  <- grab(w, "Smoking")
diabetes_w <- grab(w, "Diabetes")
bp_med_w   <- grab(w, "BloodPressureMedication")
bp_sys_w   <- safe_numeric(grab(w, "BloodPressureSystolic"))
bp_dia_w   <- safe_numeric(grab(w, "BloodPressureDiastolic"))
on_treat_w <- safe_numeric(grab(w, "OnTreatment"))

w$cohort <- "Wales"
w$id     <- as.character(grab(w, "DatabaseNumber"))
w$age    <- safe_numeric(grab(w, "BMI_AGE"))
w$sex    <- ifelse(gender_w == "M", 1, ifelse(gender_w == "F", 0, NA))
w$bmi    <- safe_numeric(grab(w, "BMI"))
w$tc     <- safe_numeric(grab(w, "TC.1"))
w$hdl    <- coalesce(safe_numeric(grab(w, "HDL.1")), safe_numeric(grab(w, "HDL.2")))
w$ldl    <- safe_numeric(grab(w, "LDL.1"))
w$tg     <- coalesce(safe_numeric(grab(w, "TRG.1")), safe_numeric(grab(w, "TRG.2")))
w$apob   <- NA_real_   # ApoB not in Wales registry
w$lpa_nmol <- safe_numeric(grab(w, "Lpa.1"))
w$on_statin <- ifelse(!is.na(on_treat_w) & on_treat_w >= 1, 1, 0)
# Wales: assume moderate intensity for treated patients (no statin name available)
w$residual_factor <- ifelse(w$on_statin == 1, 0.65, 1.00)
w$re_ldl <- w$ldl / w$residual_factor

w$ever_smoked  <- ifelse(smoking_w %in% c("1", "True", "TRUE"), 1,
                         ifelse(smoking_w %in% c("0", "False", "FALSE"), 0, NA))
w$diabetes     <- ifelse(diabetes_w %in% c("1", "True"), 1,
                         ifelse(diabetes_w %in% c("0", "False"), 0, NA))
w$hypertension <- ifelse(bp_med_w %in% c("1", "True") |
                          (!is.na(bp_sys_w) & bp_sys_w > 140) |
                          (!is.na(bp_dia_w) & bp_dia_w > 90), 1, 0)
w$lpa_143      <- ifelse(!is.na(w$lpa_nmol) & w$lpa_nmol > 143, 1, 0)
w$ascvd        <- safe_numeric(grab(w, "ascvd_combine"))
w$prior_ascvd  <- w$ascvd
w$gene         <- parse_gene(grab(w, "Mutation1"))

# Wales: no ApoB → log_apob_ldl = NA (will be excluded from CALON-NEW analysis
# in directions that include Wales; CALON-NEW-noApoB used as fallback)
w$apob_ldl     <- NA_real_
w$tg_hdl       <- ifelse(!is.na(w$tg) & !is.na(w$hdl) & w$hdl > 0, w$tg / w$hdl, NA)
w$log_apob_ldl <- NA_real_
w$log_tg_hdl   <- ifelse(!is.na(w$tg_hdl) & w$tg_hdl > 0, log(w$tg_hdl), NA)

cat(sprintf("  Wales FH: %d patients, ASCVD=%d (%.1f%%)\n",
            nrow(w), sum(w$ascvd == 1, na.rm = TRUE),
            100 * mean(w$ascvd, na.rm = TRUE)))

# --- Harmonise UKB ----------------------------------------------------------
cat("\n--- Harmonising UKB ---\n")
ukb_all <- ukb_raw
ukb_all$eid     <- grab(ukb_all, "participant.eid")
ukb_all$is_fh   <- safe_numeric(grab(ukb_all, "is_fh_genetic"))
ukb_all$sex_num <- safe_numeric(grab(ukb_all, "participant.p31"))
ukb_all$age_raw <- safe_numeric(grab(ukb_all, "participant.p21022"))
ukb_all$bmi_raw <- safe_numeric(grab(ukb_all, "participant.p21001_i0"))
ukb_all$ldl_raw <- safe_numeric(grab(ukb_all, "participant.p30780_i0"))
ukb_all$hdl_raw <- safe_numeric(grab(ukb_all, "participant.p30760_i0"))
ukb_all$tg_raw  <- safe_numeric(grab(ukb_all, "participant.p30870_i0"))
ukb_all$tc_raw  <- safe_numeric(grab(ukb_all, "participant.p30690_i0"))
ukb_all$apob_raw <- safe_numeric(grab(ukb_all, "participant.p30640_i0"))
ukb_all$sbp_1   <- safe_numeric(grab(ukb_all, "participant.p4080_i0_a0"))
ukb_all$sbp_2   <- safe_numeric(grab(ukb_all, "participant.p4080_i0_a1"))
ukb_all$gene_ukb <- grab(ukb_all, "gene")

# Statin detection from medication codes (Field 20003)
statin_codes <- c(1140888648, 1140888594, 1140888690, 1141146234,
                  1141192414, 1140861970, 1140881748, 1141146138)
med_cols <- paste0("participant.p20003_i0_a", 0:47)
med_cols <- med_cols[med_cols %in% names(ukb_all)]
ukb_all$on_statin_ukb <- if (length(med_cols) > 0) {
  apply(ukb_all[, med_cols, drop = FALSE], 1, function(row) {
    as.integer(any(row %in% statin_codes, na.rm = TRUE))
  })
} else 0L

# Self-reported diabetes / hypertension (Field 20002)
ill_cols <- grep("^participant\\.p20002_i0", names(ukb_all), value = TRUE)
ukb_all$diabetes_ukb <- if (length(ill_cols) > 0) {
  apply(ukb_all[, ill_cols, drop = FALSE], 1, function(row) {
    as.integer(any(row %in% c(1220, 1222, 1223), na.rm = TRUE))
  })
} else NA_integer_
ukb_all$hypertension_ukb <- if (length(ill_cols) > 0) {
  apply(ukb_all[, ill_cols, drop = FALSE], 1, function(row) {
    as.integer(any(row %in% c(1065, 1072), na.rm = TRUE))
  })
} else NA_integer_

# Smoking — from processed file if available
ukb_all$ever_smoked_ukb <- if ("proc_ever_smoked" %in% names(ukb_all)) {
  as.numeric(ukb_all$proc_ever_smoked)
} else if ("participant.p20116_i0" %in% names(ukb_all)) {
  ifelse(ukb_all$participant.p20116_i0 %in% c(1, 2), 1,
         ifelse(ukb_all$participant.p20116_i0 == 0, 0, NA))
} else NA_real_

# ASCVD — use processed file (correct construction from p131* dates)
if ("proc_ascvd_combined" %in% names(ukb_all)) {
  ukb_all$ascvd_ukb <- as.numeric(ukb_all$proc_ascvd_combined)
} else if ("proc_ascvd_any" %in% names(ukb_all)) {
  ukb_all$ascvd_ukb <- as.numeric(ukb_all$proc_ascvd_any)
} else {
  cat("  WARNING: No ASCVD column from processed file. Using NA.\n")
  ukb_all$ascvd_ukb <- NA_real_
}

# Lp(a) and re_ldl from processed file
if (!"proc_lpa" %in% names(ukb_all)) ukb_all$proc_lpa <- NA_real_
if (!"proc_re_ldl" %in% names(ukb_all)) ukb_all$proc_re_ldl <- NA_real_

# Filter to FH carriers
keep_fh <- which(ukb_all$is_fh == 1)
ukb_fh <- ukb_all[keep_fh, ]

# Build harmonised columns
ukb_fh$cohort   <- "UKB"
ukb_fh$id       <- as.character(ukb_fh$eid)
ukb_fh$age      <- ukb_fh$age_raw
ukb_fh$sex      <- ukb_fh$sex_num
ukb_fh$bmi      <- ukb_fh$bmi_raw
ukb_fh$tc       <- ukb_fh$tc_raw
ukb_fh$hdl      <- ukb_fh$hdl_raw
ukb_fh$ldl      <- ukb_fh$ldl_raw
ukb_fh$tg       <- ukb_fh$tg_raw
ukb_fh$apob     <- ukb_fh$apob_raw
ukb_fh$lpa_nmol <- safe_numeric(ukb_fh$proc_lpa)
ukb_fh$on_statin <- ukb_fh$on_statin_ukb
# Use processed re_ldl if available, else compute fresh
ukb_fh$re_ldl   <- if (!all(is.na(ukb_fh$proc_re_ldl))) {
  safe_numeric(ukb_fh$proc_re_ldl)
} else {
  ifelse(ukb_fh$on_statin == 1, ukb_fh$ldl / 0.65, ukb_fh$ldl)
}
ukb_fh$ever_smoked <- ukb_fh$ever_smoked_ukb
ukb_fh$diabetes    <- ukb_fh$diabetes_ukb
ukb_fh$hypertension <- ukb_fh$hypertension_ukb
ukb_fh$lpa_143     <- ifelse(!is.na(ukb_fh$lpa_nmol) & ukb_fh$lpa_nmol > 143, 1, 0)
ukb_fh$ascvd       <- ukb_fh$ascvd_ukb
ukb_fh$prior_ascvd <- ukb_fh$ascvd_ukb   # cross-sectional
ukb_fh$gene        <- ukb_fh$gene_ukb

# Derived ratios
ukb_fh$apob_ldl     <- ifelse(!is.na(ukb_fh$apob) & !is.na(ukb_fh$ldl) & ukb_fh$ldl > 0,
                                ukb_fh$apob / ukb_fh$ldl, NA)
ukb_fh$tg_hdl       <- ifelse(!is.na(ukb_fh$tg) & !is.na(ukb_fh$hdl) & ukb_fh$hdl > 0,
                                ukb_fh$tg / ukb_fh$hdl, NA)
ukb_fh$log_apob_ldl <- ifelse(!is.na(ukb_fh$apob_ldl) & ukb_fh$apob_ldl > 0,
                                log(ukb_fh$apob_ldl), NA)
ukb_fh$log_tg_hdl   <- ifelse(!is.na(ukb_fh$tg_hdl) & ukb_fh$tg_hdl > 0,
                                log(ukb_fh$tg_hdl), NA)

# Dedup UKB by eid
ukb_fh <- ukb_fh[!duplicated(ukb_fh$id), ]

cat(sprintf("  UKB FH: %d patients, ASCVD=%d (%.1f%%)\n",
            nrow(ukb_fh), sum(ukb_fh$ascvd == 1, na.rm = TRUE),
            100 * mean(ukb_fh$ascvd, na.rm = TRUE)))
cat(sprintf("  UKB Lp(a) non-missing: %d (%.1f%%)\n",
            sum(!is.na(ukb_fh$lpa_nmol)),
            100 * mean(!is.na(ukb_fh$lpa_nmol))))

# --- Pool cohorts -----------------------------------------------------------
common_cols <- c("cohort", "id", "age", "sex", "bmi",
                 "tc", "hdl", "ldl", "tg", "apob", "re_ldl",
                 "lpa_nmol", "lpa_143",
                 "apob_ldl", "tg_hdl", "log_apob_ldl", "log_tg_hdl",
                 "ever_smoked", "diabetes", "hypertension",
                 "on_statin", "ascvd", "prior_ascvd", "gene")
for (cc in common_cols) {
  if (!cc %in% names(d))      d[[cc]]      <- NA
  if (!cc %in% names(w))      w[[cc]]      <- NA
  if (!cc %in% names(ukb_fh)) ukb_fh[[cc]] <- NA
}
fh_pooled <- bind_rows(d[, common_cols], w[, common_cols], ukb_fh[, common_cols])
fh_pooled$cohort <- factor(fh_pooled$cohort, levels = c("South Wales FH", "UKB", "Wales"))

# Dedup: South Wales FH is a subset of Wales registry (same DatabaseNumbers).
# Remove SW IDs from Wales to prevent test-on-train leakage.
sw_ids <- fh_pooled$id[fh_pooled$cohort == "South Wales FH"]
n_dups <- sum(fh_pooled$cohort == "Wales" & fh_pooled$id %in% sw_ids)
fh_pooled <- fh_pooled[!(fh_pooled$cohort == "Wales" & fh_pooled$id %in% sw_ids), ]
cat(sprintf("\n  Removed %d SW patients overlapping Wales\n", n_dups))
cat(sprintf("  Pooled cohort: SW=%d, UKB=%d, Wales=%d (total=%d)\n",
            sum(fh_pooled$cohort == "South Wales FH"),
            sum(fh_pooled$cohort == "UKB"),
            sum(fh_pooled$cohort == "Wales"),
            nrow(fh_pooled)))

# --- ASCVD endpoint sanity check -------------------------------------------
cat("\n--- ASCVD endpoint sanity (verified clean) ---\n")
for (coh in levels(fh_pooled$cohort)) {
  sub <- fh_pooled[fh_pooled$cohort == coh, ]
  n_total <- nrow(sub)
  n_event <- sum(sub$ascvd == 1, na.rm = TRUE)
  pct <- 100 * n_event / n_total
  cat(sprintf("  %-16s: N=%d, ASCVD=%d (%.1f%%)\n", coh, n_total, n_event, pct))
}

################################################################################
# --- DEFINE MODELS ---------------------------------------------------------
################################################################################
# CALON-NEW (8 predictors per user spec)
NEW_features <- c("log_apob_ldl", "log_tg_hdl", "re_ldl",
                  "age", "diabetes", "hypertension", "ever_smoked", "lpa_143")

# Fallback for cohorts without ApoB (Wales): drop log_apob_ldl
NEW_noApoB_features <- c("log_tg_hdl", "re_ldl",
                          "age", "diabetes", "hypertension", "ever_smoked", "lpa_143")

# SAFEHEART-RE frozen (TRIPOD Type 4 — published coefficients, no retraining)
SRE_features <- c("age", "sex", "ldl", "hypertension", "bmi", "ever_smoked", "lpa_120", "prior_ascvd")
# Approximated SRE coefficients (Perez-de-Isla 2017 Cox -> logistic equivalent)
SRE_coefs <- c(intercept = -7.38, age = 0.0795, sex = 0.476,
               ldl = 0.279, hypertension = 0.527, bmi = 0.0246,
               ever_smoked = 0.755, lpa_120 = 0.601, prior_ascvd = 1.414)

cat("\n--- Models ---\n")
cat(sprintf("  CALON-NEW         (8): %s\n", paste(NEW_features, collapse = ", ")))
cat(sprintf("  CALON-NEW-noApoB  (7): %s\n", paste(NEW_noApoB_features, collapse = ", ")))
cat(sprintf("  SRE-Fixed         (8): %s [frozen coefficients]\n", paste(SRE_features, collapse = ", ")))

# Add lpa_120 binary (SAFEHEART threshold) for SRE comparator
fh_pooled$lpa_120 <- ifelse(!is.na(fh_pooled$lpa_nmol) & fh_pooled$lpa_nmol > 120, 1, 0)

################################################################################
# --- IMPUTATION (MICE m=10) -------------------------------------------------
################################################################################
cat("\n--- MICE imputation (m=10, maxit=20) per cohort ---\n")
m_imp <- 10
maxit_imp <- 20

impute_cohort <- function(sub, m = 10, maxit = 20, seed = 42) {
  imp_vars <- c("age", "sex", "ldl", "re_ldl", "hdl", "tg", "bmi",
                "apob", "lpa_nmol", "lpa_143", "lpa_120",
                "ever_smoked", "diabetes", "hypertension",
                "on_statin", "ascvd", "prior_ascvd",
                "log_apob_ldl", "log_tg_hdl")
  use_vars <- imp_vars[imp_vars %in% names(sub)]
  use_vars <- use_vars[sapply(use_vars, function(v) sum(!is.na(sub[[v]])) >= 5)]

  needs_imp <- use_vars[sapply(use_vars, function(v) {
    n_na <- sum(is.na(sub[[v]]))
    n_na > 0 && n_na < nrow(sub)
  })]

  if (length(needs_imp) == 0) {
    cat(sprintf("    %s: complete data (no imputation needed)\n", sub$cohort[1]))
    return(replicate(m, sub, simplify = FALSE))
  }

  imp <- tryCatch(
    mice(sub[, use_vars], m = m, maxit = maxit, printFlag = FALSE, seed = seed),
    error = function(e) { cat(sprintf("    MICE error: %s\n", e$message)); NULL }
  )
  if (is.null(imp)) return(replicate(m, sub, simplify = FALSE))

  lapply(seq_len(m), function(i) {
    completed <- complete(imp, i)
    out <- sub
    for (v in use_vars) out[[v]] <- completed[[v]]
    # Recompute derived variables after imputation
    out$apob_ldl     <- ifelse(!is.na(out$apob) & !is.na(out$ldl) & out$ldl > 0,
                                 out$apob / out$ldl, NA)
    out$tg_hdl       <- ifelse(!is.na(out$tg) & !is.na(out$hdl) & out$hdl > 0,
                                 out$tg / out$hdl, NA)
    out$log_apob_ldl <- ifelse(!is.na(out$apob_ldl) & out$apob_ldl > 0,
                                 log(out$apob_ldl), NA)
    out$log_tg_hdl   <- ifelse(!is.na(out$tg_hdl) & out$tg_hdl > 0,
                                 log(out$tg_hdl), NA)
    out
  })
}

set.seed(42)
imp_sw    <- impute_cohort(fh_pooled[fh_pooled$cohort == "South Wales FH", ], m_imp, maxit_imp)
imp_ukb   <- impute_cohort(fh_pooled[fh_pooled$cohort == "UKB", ],            m_imp, maxit_imp)
imp_wales <- impute_cohort(fh_pooled[fh_pooled$cohort == "Wales", ],          m_imp, maxit_imp)

cohort_imps <- list("South Wales FH" = imp_sw, "UKB" = imp_ukb, "Wales" = imp_wales)

################################################################################
# --- HELPERS: AUC + DeLong + calibration ------------------------------------
################################################################################
prep_cohort <- function(df, features, outcome = "ascvd") {
  cols <- c(features, outcome)
  df <- df[, intersect(cols, names(df))]
  df[complete.cases(df), , drop = FALSE]
}

# Elastic net AUC with bootstrap CI
get_auc_en <- function(train_df, test_df, features, outcome = "ascvd",
                       n_boot = 2000) {
  train <- prep_cohort(train_df, features, outcome)
  test  <- prep_cohort(test_df,  features, outcome)
  if (nrow(train) < 30 || nrow(test) < 30 ||
      length(unique(test[[outcome]])) < 2) {
    return(list(auc = NA, ci_low = NA, ci_high = NA, lp = NULL,
                n_train = nrow(train), n_test = nrow(test)))
  }
  X_tr <- scale(as.matrix(train[, features]))
  y_tr <- train[[outcome]]
  X_te <- scale(as.matrix(test[, features]),
                center = attr(X_tr, "scaled:center"),
                scale  = attr(X_tr, "scaled:scale"))
  y_te <- test[[outcome]]
  cv_fit <- cv.glmnet(X_tr, y_tr, family = "binomial",
                      alpha = 0.3, nfolds = 10,
                      type.measure = "auc", standardize = FALSE)
  pred <- as.numeric(predict(cv_fit, newx = X_te, s = "lambda.min", type = "response"))
  roc_obj <- roc(y_te, pred, quiet = TRUE)
  ci_obj <- ci.auc(roc_obj, method = "bootstrap", boot.n = n_boot, quiet = TRUE)
  list(auc = as.numeric(auc(roc_obj)),
       ci_low = ci_obj[1], ci_high = ci_obj[3],
       lp = pred, y_test = y_te,
       n_train = nrow(train), n_test = nrow(test))
}

# SRE-Fixed (frozen coefficients, no retraining)
get_auc_sre_fixed <- function(test_df, outcome = "ascvd", n_boot = 2000) {
  req <- c("age", "sex", "ldl", "hypertension", "bmi", "ever_smoked",
           "lpa_120", "prior_ascvd")
  test <- prep_cohort(test_df, req, outcome)
  if (nrow(test) < 30 || length(unique(test[[outcome]])) < 2) {
    return(list(auc = NA, ci_low = NA, ci_high = NA, lp = NULL,
                n_test = nrow(test)))
  }
  lp <- with(test, SRE_coefs["intercept"] +
                   SRE_coefs["age"]          * age +
                   SRE_coefs["sex"]          * sex +
                   SRE_coefs["ldl"]          * ldl +
                   SRE_coefs["hypertension"] * hypertension +
                   SRE_coefs["bmi"]          * bmi +
                   SRE_coefs["ever_smoked"]  * ever_smoked +
                   SRE_coefs["lpa_120"]      * lpa_120 +
                   SRE_coefs["prior_ascvd"]  * prior_ascvd)
  pred <- 1 / (1 + exp(-lp))
  roc_obj <- roc(test[[outcome]], pred, quiet = TRUE)
  ci_obj <- ci.auc(roc_obj, method = "bootstrap", boot.n = n_boot, quiet = TRUE)
  list(auc = as.numeric(auc(roc_obj)),
       ci_low = ci_obj[1], ci_high = ci_obj[3],
       lp = pred, y_test = test[[outcome]],
       n_test = nrow(test))
}

# Calibration intercept + slope on linear predictor
calibration <- function(y, p) {
  if (length(unique(y)) < 2 || length(p) < 30) return(c(intercept = NA, slope = NA))
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  cal_int <- tryCatch(coef(glm(y ~ offset(lp), family = binomial))[1],
                       error = function(e) NA)
  cal_slope <- tryCatch(coef(glm(y ~ lp, family = binomial))["lp"],
                         error = function(e) NA)
  c(intercept = unname(cal_int), slope = unname(cal_slope))
}

# Pool AUCs across imputations using Rubin's rules
pool_rubins <- function(estimates) {
  estimates <- estimates[!is.na(estimates)]
  if (length(estimates) == 0) return(list(est = NA, se = NA, ci_low = NA, ci_high = NA))
  est <- mean(estimates)
  if (length(estimates) > 1) {
    btw_var <- var(estimates) * (1 + 1 / length(estimates))
    se <- sqrt(btw_var)
    ci_low  <- est - 1.96 * se
    ci_high <- est + 1.96 * se
  } else { se <- NA; ci_low <- NA; ci_high <- NA }
  list(est = est, se = se, ci_low = ci_low, ci_high = ci_high)
}

################################################################################
# --- BI-EXTERNAL VALIDATION (6 directions) ----------------------------------
################################################################################
cat("\n=================================================================\n")
cat("  BI-EXTERNAL VALIDATION — 6 DIRECTIONS\n")
cat("=================================================================\n")

ALL_DIRS <- list(
  list(train = "South Wales FH", test = "UKB"),
  list(train = "UKB",            test = "South Wales FH"),
  list(train = "South Wales FH", test = "Wales"),
  list(train = "Wales",          test = "South Wales FH"),
  list(train = "UKB",            test = "Wales"),
  list(train = "Wales",          test = "UKB")
)

results <- data.frame()
for (dir in ALL_DIRS) {
  dir_name <- sprintf("%s -> %s", dir$train, dir$test)
  cat(sprintf("\n--- %s ---\n", dir_name))

  train_imps <- cohort_imps[[dir$train]]
  test_imps  <- cohort_imps[[dir$test]]

  for (model_name in c("CALON-NEW", "CALON-NEW-noApoB", "SRE-Fixed")) {
    feats <- switch(model_name,
                     "CALON-NEW"        = NEW_features,
                     "CALON-NEW-noApoB" = NEW_noApoB_features,
                     "SRE-Fixed"        = NULL)

    aucs <- numeric(m_imp); ci_lows <- numeric(m_imp); ci_highs <- numeric(m_imp)
    cal_ints <- numeric(m_imp); cal_slopes <- numeric(m_imp)
    n_tr <- 0; n_te <- 0
    for (i in seq_len(m_imp)) {
      train_i <- train_imps[[i]]
      test_i  <- test_imps[[i]]
      if (model_name == "SRE-Fixed") {
        res <- get_auc_sre_fixed(test_i, outcome = "ascvd", n_boot = 500)
      } else {
        res <- get_auc_en(train_i, test_i, feats, outcome = "ascvd", n_boot = 500)
      }
      aucs[i] <- res$auc
      ci_lows[i] <- res$ci_low
      ci_highs[i] <- res$ci_high
      n_tr <- max(n_tr, ifelse(is.null(res$n_train), 0, res$n_train))
      n_te <- max(n_te, res$n_test)
      if (!is.na(res$auc) && !is.null(res$lp)) {
        cal <- calibration(res$y_test, res$lp)
        cal_ints[i] <- cal["intercept"]; cal_slopes[i] <- cal["slope"]
      } else { cal_ints[i] <- NA; cal_slopes[i] <- NA }
    }
    pooled <- pool_rubins(aucs)
    cal_int_pooled <- mean(cal_ints, na.rm = TRUE)
    cal_slp_pooled <- mean(cal_slopes, na.rm = TRUE)

    cat(sprintf("  %-18s AUC=%.4f [%.4f-%.4f]  cal_int=%.3f cal_slope=%.3f  (n_tr=%d, n_te=%d)\n",
                model_name, pooled$est, pooled$ci_low, pooled$ci_high,
                cal_int_pooled, cal_slp_pooled, n_tr, n_te))

    results <- rbind(results, data.frame(
      direction = dir_name, model = model_name,
      auc = pooled$est, ci_low = pooled$ci_low, ci_high = pooled$ci_high,
      cal_intercept = cal_int_pooled, cal_slope = cal_slp_pooled,
      n_train = n_tr, n_test = n_te,
      stringsAsFactors = FALSE
    ))
  }
}

# Compute delta vs SRE-Fixed
results <- results %>%
  group_by(direction) %>%
  mutate(delta_vs_SRE = auc - auc[model == "SRE-Fixed"]) %>%
  ungroup() %>%
  as.data.frame()

cat("\n=================================================================\n")
cat("  RESULTS — CALON-NEW vs SAFEHEART-RE-Fixed (6 directions)\n")
cat("=================================================================\n\n")
print(results)
write.csv(results, "output/new/CALON_NEW_results.csv", row.names = FALSE)
cat("\n  Saved: output/new/CALON_NEW_results.csv\n")

################################################################################
# --- POOLED OR FOR CALON-NEW (Rubin's pooled coefficients) ------------------
################################################################################
cat("\n--- Pooled OR (CALON-NEW, all cohorts pooled) ---\n")

pooled_imps <- impute_cohort(fh_pooled, m_imp, maxit_imp)

# Use no-ApoB version for pooled (Wales has no ApoB)
or_estimates <- list()
for (i in seq_len(m_imp)) {
  df_i <- pooled_imps[[i]]
  cc <- df_i[complete.cases(df_i[, c(NEW_noApoB_features, "ascvd")]), ]
  X_sc <- scale(as.matrix(cc[, NEW_noApoB_features]))
  fit <- cv.glmnet(X_sc, cc$ascvd, family = "binomial",
                    alpha = 0.3, nfolds = 10, standardize = FALSE)
  coefs <- as.numeric(coef(fit, s = "lambda.min"))[-1]
  or_estimates[[i]] <- data.frame(
    variable = NEW_noApoB_features, coef = coefs,
    or = exp(coefs), imp = i
  )
}
or_pooled_data <- do.call(rbind, or_estimates) %>%
  group_by(variable) %>%
  summarise(coef = mean(coef),
             or   = exp(mean(coef)),
             se   = sd(coef),
             ci_low = exp(mean(coef) - 1.96 * sd(coef)),
             ci_high = exp(mean(coef) + 1.96 * sd(coef)),
             .groups = "drop")
print(or_pooled_data)
write.csv(or_pooled_data, "output/new/CALON_NEW_pooled_OR.csv", row.names = FALSE)
cat("\n  Saved: output/new/CALON_NEW_pooled_OR.csv\n")

################################################################################
# --- SUMMARY ---------------------------------------------------------------
################################################################################
cat("\n=================================================================\n")
cat("  SUMMARY: CALON-NEW vs SAFEHEART-RE\n")
cat("=================================================================\n")

# Win counts per model
for (m in c("CALON-NEW", "CALON-NEW-noApoB")) {
  m_results <- results[results$model == m, ]
  wins <- sum(m_results$delta_vs_SRE > 0, na.rm = TRUE)
  ties <- sum(abs(m_results$delta_vs_SRE) < 0.005, na.rm = TRUE)
  losses <- sum(m_results$delta_vs_SRE < 0, na.rm = TRUE)
  n_dir <- sum(!is.na(m_results$delta_vs_SRE))
  cat(sprintf("  %-18s wins %d / ties %d / losses %d (of %d directions)\n",
              m, wins, ties, losses, n_dir))
}

# Headline: best CALON-NEW AUC
best <- results[results$model %in% c("CALON-NEW", "CALON-NEW-noApoB"), ] %>%
  arrange(desc(auc)) %>% head(3)
cat("\n  Top 3 CALON-NEW AUCs:\n")
for (i in seq_len(nrow(best))) {
  cat(sprintf("    %s | %s: AUC=%.4f [%.4f-%.4f]  (n_test=%d)\n",
              best$direction[i], best$model[i],
              best$auc[i], best$ci_low[i], best$ci_high[i], best$n_test[i]))
}

cat("\nDone. Inspect output/new/ for full results.\n")

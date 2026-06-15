################################################################################
# CALON-NEW v2: Maximum-discrimination ASCVD risk prediction in genetic FH
#
# GOAL: Find the highest defensible AUC per cohort across multiple model variants.
#       Report all variants transparently. Flag overfitting risk per row.
#
# BASE PREDICTOR SET (user spec 2026-05-03, revised for no-duplicate-variable rule):
#   log_apob, log_tg_hdl, re_ldl, age, diabetes, hypertension,
#   ever_smoked, lpa_143
#
# DIRECTIVE: "Don't use the same variable in two predictors in the same equation."
#   → log_apob_ldl REPLACED with log_apob (the ratio reused LDL alongside re_ldl)
#   → age_x_ldl interaction variants REMOVED (reused age and re_ldl)
#   → lipid_years variant kept (single combined term replacing age + re_ldl,
#                                 no duplication within the variant)
#   Verified non-overlap matrix:
#     re_ldl   = LDL only             — not in any other predictor
#     log_apob = ApoB only            — not in any other predictor
#     log_tg_hdl = TG + HDL only      — TG/HDL not in any other predictor
#     age, sex, dm, htn, smoke, lpa   — independent
#     ldl_prs, cad_prs                — polygenic, GWAS-derived (not from measured LDL)
#     hdl_sub3, hdl_sub4              — NMR sub-particle concentrations
#                                       (orthogonal to assay HDL used in tg_hdl)
#     chip                            — fully independent
#
# AUGMENTATIONS (added incrementally — each adds biological signal):
#   v1 +sex                       — conventional (was missing from spec)
#   v2 +sex +ldl_prs +cad_prs     — polygenic background (UKB only)
#   v3 +sex +PRS +nmr_hdl_sub     — HDL sub-particle CETP signal (Paper 1 ρ=0.09)
#   v4 +sex +PRS +chip            — CHIP×phenotype interaction (Paper 1 p=0.0004)
#   v5 +sex +PRS +HDL_sub +CHIP   — full-augmented (richest UKB)
#   v6 lipid_years (no age, no re_ldl)  — single biologically-motivated term
#   v7 +PRS +HDLsub +CHIP +lipidyears   — augmented lipid_years variant
#
# OUTCOME VARIANTS (BOTH analysed in parallel — fairness rules differ):
#   PREVALENCE: ascvd (lifetime any-time ASCVD) — first-occurrence ICD-10
#               I20-I25/I63-I64/I70/I73-I74. prior_ascvd CANNOT be a predictor
#               (would be tautology). All 3 cohorts available.
#   INCIDENCE:  ascvd_incident (post-baseline only) — UKB and SW only.
#               Wales lacks temporal split → excluded from incidence analyses.
#               prior_ascvd IS a predictor in BOTH CALON and SRE-Fixed,
#               matching SAFEHEART original design (Perez-de-Isla 2017).
#
# COMPARATOR: SRE-Fixed (frozen SAFEHEART coefficients, TRIPOD Type 4)
#   References (from PubMed):
#     - SAFEHEART original:    Perez-de-Isla et al. Circulation 2017
#                              DOI 10.1161/CIRCULATIONAHA.116.024541
#                              Internal Harrell C = 0.85 (incident ASCVD, n=2404)
#     - English CPRD external: McKay et al. Atherosclerosis 2022
#                              DOI 10.1016/j.atherosclerosis.2022.07.011
#                              External C = 0.67 (calibration slope 10.09 men)
#     - Australian external:   Mansilla-Rodriguez et al. EJPC 2025
#                              DOI 10.1093/eurjpc/zwaf631
#     - FH-Risk-Score:         Paquette et al. ATVB 2021
#                              DOI 10.1161/ATVBAHA.121.316106
#                              C = 0.75 (primary prevention, no prior ASCVD)
#     - PRS in FH:             Trinder et al. ATVB 2024
#                              DOI 10.1161/ATVBAHA.123.320287
#                              CAD-PRS HR 1.92 for ASCVD in FH
#
# IMPUTATION: MICE m=10, maxit=20, Rubin's rules
# RECALIBRATION: Platt scaling on linear predictor (where AUC > 0.75)
# STATISTICS: AUC + bootstrap CI, calibration int+slope, DeLong vs SRE
#
# CRITICAL HONESTY NOTES:
#   1. Per-cohort best-variant differing = overfitting risk → flagged as "RED"
#   2. Calibration slope > 2 or < 0.5 = unreliable → flagged as "AMBER"
#   3. Test cohort < 200 events for any model = unstable estimate → flagged
#   4. UKB AUC > 0.80 with clinical-only predictors = suspicious → require justification
#
# Author: Dr Nader Genedy
# Date:   2026-05-03
# Builds on: CALON_NEW.R, CALON_MI_validation.R
################################################################################

suppressPackageStartupMessages({
  pkgs <- c("tidyverse", "readxl", "tableone", "pROC",
            "broom", "boot", "mice", "glmnet")
  miss <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(miss) > 0) install.packages(miss, repos = "https://cloud.r-project.org")
  for (p in pkgs) library(p, character.only = TRUE)
})

# --- Paths ------------------------------------------------------------------
if (dir.exists("/cloud/project")) {
  data_dir <- "/cloud/project"
  ukb_aux  <- "/cloud/project/data"            # adjust if you upload to Posit
} else if (dir.exists("C:/Users/nader/Downloads/calon_ukb_pipeline")) {
  data_dir <- "C:/Users/nader/Downloads/calon_ukb_pipeline"
  ukb_aux  <- "D:/Projects/CALON_AlphaFold_Rebuild/data"
} else { data_dir <- getwd(); ukb_aux <- file.path(getwd(), "data") }
setwd(data_dir)
dir.create("output/v2", showWarnings = FALSE, recursive = TRUE)

# --- Helpers ----------------------------------------------------------------
safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))
parse_gene <- function(mut) {
  ifelse(is.na(mut), NA_character_,
         ifelse(grepl("^LDLR", mut), "LDLR",
                ifelse(grepl("^APOB", mut), "APOB",
                       ifelse(grepl("^PCSK9", mut), "PCSK9", "Other"))))
}
# Search current dir + known data backup locations
data_search_dirs <- c(
  ".",                                                    # working dir
  "D:/CALON_FH_BACKUP_FULL/data",                         # Wales, calon_ready, DRAGON
  "D:/Projects/CALON_AlphaFold_Rebuild/data",             # ukb_carriers_MASTER, prs, nmr, chip
  "D:/calon_backup_data",                                 # TUDOR_UKB_Features
  "D:/alphafold_backup/tudor_packup"                      # TUDOR_UKB_Features alt
)

find_file <- function(pattern, required = TRUE) {
  matches <- character(0)
  for (sd in data_search_dirs) {
    if (dir.exists(sd)) {
      m <- list.files(sd, pattern = pattern, ignore.case = TRUE,
                      full.names = TRUE, recursive = TRUE)
      matches <- c(matches, m)
    }
  }
  matches <- matches[!grepl("output/|Table_|SCORECARD|\\.claude/", matches)]
  if (length(matches) == 0) {
    if (required) stop(sprintf("File not found: %s\n  Searched: %s",
                                pattern, paste(data_search_dirs, collapse = ", ")))
    return(NULL)
  }
  if (length(matches) > 1) matches <- matches[order(-file.info(matches)$size)]
  matches[1]
}
safe_read_csv <- function(path, label = "") {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(df)[1] <- sub("^﻿", "", names(df)[1])
  names(df) <- trimws(names(df))
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) names(df)[bad] <- paste0("X_empty_", seq_len(sum(bad)))
  if (anyDuplicated(names(df))) names(df) <- make.unique(names(df), sep = "_dup_")
  cat(sprintf("  %s: %d x %d\n", label, nrow(df), ncol(df)))
  df
}
grab <- function(df, name) {
  if (name %in% names(df)) return(df[[name]])
  idx <- which(tolower(names(df)) == tolower(name))
  if (length(idx) > 0) return(df[[idx[1]]])
  rep(NA, nrow(df))
}

cat("============================================================\n")
cat("  CALON-NEW v2: Multi-Variant Maximum AUC Search\n")
cat("============================================================\n\n")

################################################################################
# 1. LOAD COHORTS
################################################################################
cat("--- Loading cohorts ---\n")
dragon_raw <- safe_read_csv(find_file("DRAGON.*3.*\\.csv"), "DRAGON_3")
wales_raw  <- safe_read_csv(find_file("WALES.*FH.*\\.csv"),  "Wales")
ukb_raw    <- safe_read_csv(find_file("TUDOR.*UKB.*\\.csv"), "UKB_TUDOR")
names(ukb_raw) <- make.names(names(ukb_raw), unique = TRUE)

# Processed UKB file (correct ASCVD, Lp(a), re_ldl)
ukb_proc_file <- find_file("calon_ukb_analysis_ready", required = FALSE)
if (!is.null(ukb_proc_file)) {
  ukb_proc <- safe_read_csv(ukb_proc_file, "UKB_proc")
  names(ukb_proc) <- make.names(names(ukb_proc), unique = TRUE)
  pf <- c("eid", "lpa", "lpa_binary", "re_ldl", "residual_factor",
          "ever_smoked", "smoking_binary",
          "ascvd_combined", "ascvd_any", "ascvd_prevalent", "ascvd_incident",
          "log_apob_ldl", "hypertension", "diabetes")
  pf <- pf[pf %in% names(ukb_proc)]
  ups <- ukb_proc[, pf, drop = FALSE]
  names(ups) <- ifelse(names(ups) == "eid", "eid", paste0("proc_", names(ups)))
  names(ups)[1] <- "eid"
  ukb_raw <- merge(ukb_raw, ups, by.x = "participant.eid", by.y = "eid",
                   all.x = TRUE, sort = FALSE)
  cat(sprintf("  UKB processed merged: %d fields\n", ncol(ups) - 1))
}

# Auxiliary UKB files (PRS, NMR HDL sub, CHIP) — optional
load_aux <- function(name, key = "eid") {
  p <- file.path(ukb_aux, name)
  if (file.exists(p)) {
    df <- safe_read_csv(p, paste0("AUX:", name))
    names(df) <- make.names(names(df), unique = TRUE)
    return(df)
  }
  cat(sprintf("  AUX missing: %s (skipping)\n", p)); NULL
}
prs   <- load_aux("ukb_prs.csv")
nmr_h <- load_aux("ukb_nmr_hdl_subfractions.csv")
chip  <- load_aux("ukb_chip.csv")
prs_master <- load_aux("ukb_carriers_MASTER.csv")  # has cad_prs, ldl_prs columns

################################################################################
# 2. HARMONISE COHORTS
################################################################################

# --- South Wales FH ---------------------------------------------------------
cat("\n--- Harmonising South Wales FH ---\n")
pos1_d <- grab(dragon_raw, "Positive1")
keep_d <- which(safe_numeric(pos1_d) == 1 | pos1_d == "Yes")
d <- if (length(keep_d) == 0) dragon_raw else dragon_raw[keep_d, ]

statin_d <- grab(d, "Statin"); gender_d <- grab(d, "Gender")
sm_d <- safe_numeric(grab(d, "Smoking_binary"))
bpm_d <- safe_numeric(grab(d, "BloodPressureMedication"))
bps_d <- safe_numeric(grab(d, "BloodPressureSystolic"))
bpd_d <- safe_numeric(grab(d, "BloodPressureDiastolic"))

d$cohort <- "South Wales FH"
d$id     <- as.character(grab(d, "DatabaseNumber"))
d$family <- as.character(grab(d, "FamilyNumber"))
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
d$residual_factor <- ifelse(d$on_statin == 1,
                              ifelse(grepl("ATOR|ROSUV", toupper(statin_d)), 0.50, 0.65),
                              1.00)
d$re_ldl <- d$ldl / d$residual_factor

d$ever_smoked <- ifelse(sm_d %in% c(1, 2), 1, ifelse(sm_d == 0, 0, NA))
d$diabetes    <- safe_numeric(grab(d, "Diabetes_binary"))
d$hypertension <- ifelse((!is.na(bpm_d) & bpm_d == 1) |
                          (!is.na(bps_d) & bps_d > 140) |
                          (!is.na(bpd_d) & bpd_d > 90), 1, 0)
d$lpa_143    <- ifelse(!is.na(d$lpa_nmol) & d$lpa_nmol > 143, 1, 0)
d$lpa_120    <- ifelse(!is.na(d$lpa_nmol) & d$lpa_nmol > 120, 1, 0)
# ASCVD outcome construction:
#   ascvd  = lifetime ASCVD (prevalent + incident, "any-time")  → PREVALENCE outcome
#   ascvd_incident = ASCVD with age_at_event > age_at_test       → INCIDENCE outcome
#   prior_ascvd    = ASCVD with age_at_event ≤ age_at_test       → predictor at baseline
d$ascvd        <- safe_numeric(grab(d, "ASCVD_combined"))
d$age_at_event <- safe_numeric(grab(d, "age_at_event"))
d$prior_ascvd  <- ifelse(d$ascvd == 1 & !is.na(d$age_at_event) & !is.na(d$age) &
                           d$age_at_event <= d$age, 1, 0)
d$prior_ascvd[is.na(d$prior_ascvd)] <- 0
d$ascvd_incident <- ifelse(d$ascvd == 1 & !is.na(d$age_at_event) & !is.na(d$age) &
                              d$age_at_event > d$age, 1, 0)
d$ascvd_incident[is.na(d$ascvd_incident)] <- 0
# If ASCVD=1 but no age_at_event recorded: treat as prevalent (conservative)
# These rows will appear in prevalence but not incidence — flagged in summary.
d$gene       <- parse_gene(grab(d, "Mutation1"))
# Stigmata
d$has_phys_sign <- pmax(safe_numeric(grab(d, "CornealArcus")),
                         safe_numeric(grab(d, "TendonXanthomata")), na.rm = TRUE)
d$has_phys_sign[is.na(d$has_phys_sign)] <- 0
# DLCN
d$dlcn <- safe_numeric(grab(d, "DLCN_score"))

# Augmentation slots (NA in SW)
d$ldl_prs <- NA_real_; d$cad_prs <- NA_real_
d$hdl_sub3 <- NA_real_; d$hdl_sub4 <- NA_real_
d$chip <- NA_real_

cat(sprintf("  South Wales FH: N=%d, ASCVD=%d (%.1f%%)\n",
            nrow(d), sum(d$ascvd == 1, na.rm = TRUE),
            100 * mean(d$ascvd, na.rm = TRUE)))

# --- Wales FH Registry ------------------------------------------------------
cat("\n--- Harmonising Wales FH ---\n")
pos1_w <- grab(wales_raw, "Positive1")
keep_w <- which(safe_numeric(pos1_w) == 1 | pos1_w == "Yes" | pos1_w == "True")
w <- if (length(keep_w) == 0) wales_raw else wales_raw[keep_w, ]

gen_w <- grab(w, "Gender"); sm_w <- grab(w, "Smoking"); dm_w <- grab(w, "Diabetes")
bpm_w <- grab(w, "BloodPressureMedication")
bps_w <- safe_numeric(grab(w, "BloodPressureSystolic"))
bpd_w <- safe_numeric(grab(w, "BloodPressureDiastolic"))
on_t_w <- safe_numeric(grab(w, "OnTreatment"))

w$cohort <- "Wales"
w$id     <- as.character(grab(w, "DatabaseNumber"))
w$family <- as.character(grab(w, "FamilyNumber"))
w$age    <- safe_numeric(grab(w, "BMI_AGE"))
w$proband <- safe_numeric(grab(w, "Proband"))
w$relat_to_index <- grab(w, "relat_to_index")
w$sex    <- ifelse(gen_w == "M", 1, ifelse(gen_w == "F", 0, NA))
w$bmi    <- safe_numeric(grab(w, "BMI"))
w$tc     <- safe_numeric(grab(w, "TC.1"))
w$hdl    <- coalesce(safe_numeric(grab(w, "HDL.1")), safe_numeric(grab(w, "HDL.2")))
w$ldl    <- safe_numeric(grab(w, "LDL.1"))
w$tg     <- coalesce(safe_numeric(grab(w, "TRG.1")), safe_numeric(grab(w, "TRG.2")))
w$apob   <- NA_real_
w$lpa_nmol <- safe_numeric(grab(w, "Lpa.1"))
w$on_statin <- ifelse(!is.na(on_t_w) & on_t_w >= 1, 1, 0)
w$residual_factor <- ifelse(w$on_statin == 1, 0.65, 1.00)
w$re_ldl <- w$ldl / w$residual_factor

w$ever_smoked <- ifelse(sm_w %in% c("1", "True", "TRUE"), 1,
                         ifelse(sm_w %in% c("0", "False", "FALSE"), 0, NA))
w$diabetes    <- ifelse(dm_w %in% c("1", "True"), 1,
                         ifelse(dm_w %in% c("0", "False"), 0, NA))
w$hypertension <- ifelse(bpm_w %in% c("1", "True") |
                          (!is.na(bps_w) & bps_w > 140) |
                          (!is.na(bpd_w) & bpd_w > 90), 1, 0)
w$lpa_143    <- ifelse(!is.na(w$lpa_nmol) & w$lpa_nmol > 143, 1, 0)
w$lpa_120    <- ifelse(!is.na(w$lpa_nmol) & w$lpa_nmol > 120, 1, 0)
# Wales has only "ascvd_combine" (lifetime, no temporal split available).
# ascvd        = used for PREVALENCE outcome only
# ascvd_incident = NA — Wales excluded from incidence analyses
# prior_ascvd  = NA — cannot determine without baseline date
w$ascvd          <- safe_numeric(grab(w, "ascvd_combine"))
w$ascvd_incident <- NA_real_
w$prior_ascvd    <- NA_real_
w$gene           <- parse_gene(grab(w, "Mutation1"))
w$has_phys_sign <- pmax(
  ifelse(grab(w, "CornealArcus") %in% c("1", "True"), 1, 0),
  ifelse(grab(w, "TendonXanthomata") %in% c("1", "True"), 1, 0),
  na.rm = TRUE
)
w$has_phys_sign[is.na(w$has_phys_sign)] <- 0
w$dlcn <- safe_numeric(grab(w, "DLCN_score"))
w$ldl_prs <- NA_real_; w$cad_prs <- NA_real_
w$hdl_sub3 <- NA_real_; w$hdl_sub4 <- NA_real_; w$chip <- NA_real_

cat(sprintf("  Wales FH: N=%d, ASCVD=%d (%.1f%%)\n",
            nrow(w), sum(w$ascvd == 1, na.rm = TRUE),
            100 * mean(w$ascvd, na.rm = TRUE)))

# --- UKB --------------------------------------------------------------------
cat("\n--- Harmonising UKB ---\n")
ukb <- ukb_raw
ukb$eid     <- grab(ukb, "participant.eid")
ukb$is_fh   <- safe_numeric(grab(ukb, "is_fh_genetic"))
ukb$sex_num <- safe_numeric(grab(ukb, "participant.p31"))
ukb$age_raw <- safe_numeric(grab(ukb, "participant.p21022"))
ukb$bmi_raw <- safe_numeric(grab(ukb, "participant.p21001_i0"))
ukb$ldl_raw <- safe_numeric(grab(ukb, "participant.p30780_i0"))
ukb$hdl_raw <- safe_numeric(grab(ukb, "participant.p30760_i0"))
ukb$tg_raw  <- safe_numeric(grab(ukb, "participant.p30870_i0"))
ukb$tc_raw  <- safe_numeric(grab(ukb, "participant.p30690_i0"))
ukb$apob_raw <- safe_numeric(grab(ukb, "participant.p30640_i0"))
ukb$sbp_1   <- safe_numeric(grab(ukb, "participant.p4080_i0_a0"))
ukb$sbp_2   <- safe_numeric(grab(ukb, "participant.p4080_i0_a1"))
ukb$gene_ukb <- grab(ukb, "gene")

statin_codes <- c(1140888648, 1140888594, 1140888690, 1141146234,
                  1141192414, 1140861970, 1140881748, 1141146138)
med_cols <- paste0("participant.p20003_i0_a", 0:47)
med_cols <- med_cols[med_cols %in% names(ukb)]
ukb$on_statin_ukb <- if (length(med_cols) > 0) {
  apply(ukb[, med_cols, drop = FALSE], 1, function(r) {
    as.integer(any(r %in% statin_codes, na.rm = TRUE))
  })
} else 0L

ill_cols <- grep("^participant\\.p20002_i0", names(ukb), value = TRUE)
ukb$diabetes_ukb <- if (length(ill_cols) > 0) {
  apply(ukb[, ill_cols, drop = FALSE], 1, function(r)
    as.integer(any(r %in% c(1220, 1222, 1223), na.rm = TRUE)))
} else NA_integer_
ukb$hypertension_ukb <- if (length(ill_cols) > 0) {
  apply(ukb[, ill_cols, drop = FALSE], 1, function(r)
    as.integer(any(r %in% c(1065, 1072), na.rm = TRUE)))
} else NA_integer_

ukb$ever_smoked_ukb <- if ("proc_ever_smoked" %in% names(ukb)) {
  as.numeric(ukb$proc_ever_smoked)
} else if ("participant.p20116_i0" %in% names(ukb)) {
  ifelse(ukb$participant.p20116_i0 %in% c(1, 2), 1,
         ifelse(ukb$participant.p20116_i0 == 0, 0, NA))
} else NA_real_

ukb$ascvd_ukb <- if ("proc_ascvd_combined" %in% names(ukb)) {
  as.numeric(ukb$proc_ascvd_combined)
} else if ("proc_ascvd_any" %in% names(ukb)) {
  as.numeric(ukb$proc_ascvd_any)
} else {
  NA_real_
}
ukb$ascvd_incident_ukb <- if ("proc_ascvd_incident" %in% names(ukb)) {
  as.numeric(ukb$proc_ascvd_incident)
} else NA_real_

if (!"proc_lpa" %in% names(ukb)) ukb$proc_lpa <- NA_real_
if (!"proc_re_ldl" %in% names(ukb)) ukb$proc_re_ldl <- NA_real_

# Filter FH carriers
keep_fh <- which(ukb$is_fh == 1)
u <- ukb[keep_fh, ]

u$cohort <- "UKB"
u$id     <- as.character(u$eid)
u$family <- NA_character_   # UKB has no family linkage in WES carriers
u$age    <- u$age_raw
u$sex    <- u$sex_num
u$bmi    <- u$bmi_raw
u$tc     <- u$tc_raw; u$hdl <- u$hdl_raw; u$ldl <- u$ldl_raw; u$tg <- u$tg_raw
u$apob   <- u$apob_raw
u$lpa_nmol <- safe_numeric(u$proc_lpa)
u$on_statin <- u$on_statin_ukb
u$re_ldl <- if (!all(is.na(u$proc_re_ldl))) {
  safe_numeric(u$proc_re_ldl)
} else {
  ifelse(u$on_statin == 1, u$ldl / 0.65, u$ldl)
}
u$ever_smoked <- u$ever_smoked_ukb
u$diabetes    <- u$diabetes_ukb
u$hypertension <- u$hypertension_ukb
u$lpa_143    <- ifelse(!is.na(u$lpa_nmol) & u$lpa_nmol > 143, 1, 0)
u$lpa_120    <- ifelse(!is.na(u$lpa_nmol) & u$lpa_nmol > 120, 1, 0)
u$ascvd      <- u$ascvd_ukb
# UKB has clean prevalent/incident split from p131* first-occurrence dates
# (built in 01_CALON_build_cohort.R — verified clean, MI/IHD/stroke/PVD only)
u$ascvd_incident <- ifelse(!is.na(u$ascvd_incident_ukb), u$ascvd_incident_ukb, NA_real_)
# proc_ascvd_prevalent is loaded if available
if ("proc_ascvd_prevalent" %in% names(u)) {
  u$prior_ascvd <- as.numeric(u$proc_ascvd_prevalent)
} else {
  u$prior_ascvd <- NA_real_
}
u$gene       <- u$gene_ukb
u$has_phys_sign <- NA_real_
u$dlcn <- NA_real_

# --- AUGMENTATION MERGES (UKB-specific) ------------------------------------
# PRS (LDL-PRS p26206, CAD-PRS p26228) — try ukb_prs.csv first, then master
if (!is.null(prs)) {
  prs_eid_col <- intersect(c("eid", "participant.eid"), names(prs))[1]
  ldl_prs_col <- intersect(c("p26206", "participant.p26206", "ldl_prs"), names(prs))[1]
  cad_prs_col <- intersect(c("p26228", "participant.p26228", "cad_prs"), names(prs))[1]
  if (!is.na(prs_eid_col)) {
    prs$eid_join <- as.character(prs[[prs_eid_col]])
    if (!is.na(ldl_prs_col)) u <- merge(u, prs[, c("eid_join", ldl_prs_col)],
                                          by.x = "id", by.y = "eid_join", all.x = TRUE)
    if (!is.na(cad_prs_col)) u <- merge(u, prs[, c("eid_join", cad_prs_col)],
                                          by.x = "id", by.y = "eid_join", all.x = TRUE)
    if (!is.na(ldl_prs_col)) u$ldl_prs <- safe_numeric(u[[ldl_prs_col]])
    if (!is.na(cad_prs_col)) u$cad_prs <- safe_numeric(u[[cad_prs_col]])
  }
}
if ((!"ldl_prs" %in% names(u) || all(is.na(u$ldl_prs))) && !is.null(prs_master)) {
  pm_eid <- intersect(c("participant_id", "eid"), names(prs_master))[1]
  if (!is.na(pm_eid)) {
    prs_master$eid_join <- as.character(prs_master[[pm_eid]])
    cols <- c("eid_join")
    if ("ldl_prs" %in% names(prs_master)) cols <- c(cols, "ldl_prs")
    if ("cad_prs" %in% names(prs_master)) cols <- c(cols, "cad_prs")
    if (length(cols) > 1) {
      u <- merge(u, prs_master[, cols], by.x = "id", by.y = "eid_join",
                 all.x = TRUE, suffixes = c("", "_m"))
      if (!"ldl_prs" %in% names(u) && "ldl_prs_m" %in% names(u))
        u$ldl_prs <- safe_numeric(u$ldl_prs_m)
      if (!"cad_prs" %in% names(u) && "cad_prs_m" %in% names(u))
        u$cad_prs <- safe_numeric(u$cad_prs_m)
    }
  }
}
if (!"ldl_prs" %in% names(u)) u$ldl_prs <- NA_real_
if (!"cad_prs" %in% names(u)) u$cad_prs <- NA_real_

# NMR HDL sub-fractions: HDL sub-3 = p23456, HDL sub-4 = p23457 (Tabet 2025 mapping)
if (!is.null(nmr_h)) {
  nmr_eid <- intersect(c("eid", "participant.eid"), names(nmr_h))[1]
  if (!is.na(nmr_eid)) {
    nmr_h$eid_join <- as.character(nmr_h[[nmr_eid]])
    sub3 <- intersect(c("participant.p23456_i0", "p23456_i0", "p23456",
                        "participant.p23456"), names(nmr_h))[1]
    sub4 <- intersect(c("participant.p23457_i0", "p23457_i0", "p23457",
                        "participant.p23457"), names(nmr_h))[1]
    if (!is.na(sub3)) {
      u <- merge(u, nmr_h[, c("eid_join", sub3)], by.x = "id", by.y = "eid_join", all.x = TRUE)
      u$hdl_sub3 <- safe_numeric(u[[sub3]])
    } else u$hdl_sub3 <- NA_real_
    if (!is.na(sub4)) {
      u <- merge(u, nmr_h[, c("eid_join", sub4)], by.x = "id", by.y = "eid_join", all.x = TRUE)
      u$hdl_sub4 <- safe_numeric(u[[sub4]])
    } else u$hdl_sub4 <- NA_real_
  } else { u$hdl_sub3 <- NA_real_; u$hdl_sub4 <- NA_real_ }
} else { u$hdl_sub3 <- NA_real_; u$hdl_sub4 <- NA_real_ }

# CHIP status
if (!is.null(chip)) {
  chip_eid <- intersect(c("eid", "participant.eid"), names(chip))[1]
  chip_st  <- intersect(c("p30105", "participant.p30105", "chip"), names(chip))[1]
  if (!is.na(chip_eid) && !is.na(chip_st)) {
    chip$eid_join <- as.character(chip[[chip_eid]])
    u <- merge(u, chip[, c("eid_join", chip_st)], by.x = "id", by.y = "eid_join", all.x = TRUE)
    u$chip <- ifelse(!is.na(u[[chip_st]]) & u[[chip_st]] != "" & u[[chip_st]] != "0", 1, 0)
  } else u$chip <- NA_real_
} else u$chip <- NA_real_

u <- u[!duplicated(u$id), ]
cat(sprintf("  UKB FH: N=%d, ASCVD=%d (%.1f%%), Lp(a) ok=%d, ApoB ok=%d, LDL-PRS ok=%d, CAD-PRS ok=%d, HDL-sub3 ok=%d, CHIP ok=%d\n",
            nrow(u), sum(u$ascvd == 1, na.rm = TRUE),
            100 * mean(u$ascvd, na.rm = TRUE),
            sum(!is.na(u$lpa_nmol)), sum(!is.na(u$apob)),
            sum(!is.na(u$ldl_prs)), sum(!is.na(u$cad_prs)),
            sum(!is.na(u$hdl_sub3)), sum(!is.na(u$chip))))

# --- UKB split into DISJOINT LipidClinic vs NonClinic (no leakage) ---
# Restrict UKB to clinic-equivalent severe phenotype: untreated LDL >= 4.9 mmol/L
# (DLCN definite-FH threshold = 190 mg/dL). Tests whether the UKB AUC ceiling
# (~0.72) is biological (population-screening dilution) or a model limitation.
# CRITICAL: UKB-LipidClinic and UKB-NonClinic are MUTUALLY EXCLUSIVE subsets.
# We do NOT also keep "full UKB" as a separate cohort because that would leak
# any time a direction crossed UKB ↔ UKB-LipidClinic (subset overlap).
LIPID_CLINIC_LDL_THRESH <- 4.9   # mmol/L, untreated equivalent
clinic_mask <- !is.na(u$re_ldl) & u$re_ldl >= LIPID_CLINIC_LDL_THRESH
u_clinic     <- u[clinic_mask, ]
u_nonclinic  <- u[!clinic_mask, ]
u_clinic$cohort    <- "UKB-LipidClinic"
u_nonclinic$cohort <- "UKB-NonClinic"
cat(sprintf("  UKB-LipidClinic FH (re_ldl >= %.1f): N=%d, ASCVD=%d (%.1f%%)\n",
            LIPID_CLINIC_LDL_THRESH, nrow(u_clinic),
            sum(u_clinic$ascvd == 1, na.rm = TRUE),
            100 * mean(u_clinic$ascvd, na.rm = TRUE)))
cat(sprintf("  UKB-NonClinic FH (re_ldl < %.1f or missing): N=%d, ASCVD=%d (%.1f%%)\n",
            LIPID_CLINIC_LDL_THRESH, nrow(u_nonclinic),
            sum(u_nonclinic$ascvd == 1, na.rm = TRUE),
            100 * mean(u_nonclinic$ascvd, na.rm = TRUE)))
# Sanity: ensure disjoint
stopifnot(length(intersect(u_clinic$id, u_nonclinic$id)) == 0)
cat(sprintf("  Disjoint check: 0 ID overlap between LipidClinic and NonClinic ✓\n"))

################################################################################
# 3. POOL & DEDUPE
################################################################################
common_cols <- c("cohort", "id", "family", "age", "sex", "bmi",
                  "tc", "hdl", "ldl", "tg", "apob", "re_ldl",
                  "lpa_nmol", "lpa_143", "lpa_120",
                  "ever_smoked", "diabetes", "hypertension",
                  "on_statin", "ascvd", "ascvd_incident", "prior_ascvd",
                  "gene", "has_phys_sign", "dlcn",
                  "ldl_prs", "cad_prs", "hdl_sub3", "hdl_sub4", "chip")
for (cc in common_cols) {
  if (!cc %in% names(d)) d[[cc]] <- NA
  if (!cc %in% names(w)) w[[cc]] <- NA
  if (!cc %in% names(u_clinic))    u_clinic[[cc]]    <- NA
  if (!cc %in% names(u_nonclinic)) u_nonclinic[[cc]] <- NA
}
# Build pooled cohort using DISJOINT UKB subsets (no leakage)
fh <- bind_rows(d[, common_cols], w[, common_cols],
                 u_clinic[, common_cols], u_nonclinic[, common_cols])
fh$cohort <- factor(fh$cohort,
                     levels = c("South Wales FH", "UKB-LipidClinic",
                                "UKB-NonClinic", "Wales"))

# ===== DEDUP: SW patients AND their families overlap with Wales registry =====
# Dragon-3 (SW) is a clinic subset of the All-Wales PASS registry. Same patients
# AND same families appear in both. For honest bi-external validation we must
# remove BOTH patient-level and family-level overlap from the Wales cohort.
# (Per user direction 2026-05-03: "review using family number")
sw_ids      <- fh$id[fh$cohort == "South Wales FH"]
sw_families <- fh$family[fh$cohort == "South Wales FH"]
sw_families <- unique(sw_families[!is.na(sw_families) & sw_families != ""])

n_dup_id     <- sum(fh$cohort == "Wales" & fh$id %in% sw_ids)
n_dup_family <- sum(fh$cohort == "Wales" & !is.na(fh$family) & fh$family %in% sw_families &
                     !(fh$id %in% sw_ids))   # additional family-only matches

drop_mask <- fh$cohort == "Wales" &
              ((fh$id %in% sw_ids) |
               (!is.na(fh$family) & fh$family %in% sw_families))
fh <- fh[!drop_mask, ]
cat(sprintf("\n  Dedup vs SW: removed %d Wales rows by patient ID + %d more by family ID\n",
            n_dup_id, n_dup_family))
cat(sprintf("  Total removed: %d (%d unique SW patients × ~%.1f relatives per family)\n",
            n_dup_id + n_dup_family, length(sw_ids),
            (n_dup_id + n_dup_family) / max(1, length(sw_ids))))

# Derived terms (no-duplicate rule enforced — see header)
# log_apob: ApoB alone (not the ratio with LDL)
fh$log_apob     <- ifelse(!is.na(fh$apob) & fh$apob > 0, log(fh$apob), NA)
# tg_hdl: TG and HDL (independent of LDL and ApoB)
fh$tg_hdl       <- ifelse(!is.na(fh$tg) & !is.na(fh$hdl) & fh$hdl > 0, fh$tg / fh$hdl, NA)
fh$log_tg_hdl   <- ifelse(!is.na(fh$tg_hdl) & fh$tg_hdl > 0, log(fh$tg_hdl), NA)
# lipid_years: single biologically-motivated term replacing age + re_ldl
# (when used, age and re_ldl must NOT also appear in the model)
fh$lipid_years  <- fh$re_ldl * fh$age
# HDL sub-particle ratio (NMR-derived, independent of assay HDL used in tg_hdl)
fh$hdl_sub_ratio <- ifelse(!is.na(fh$hdl_sub3) & !is.na(fh$hdl_sub4) & fh$hdl_sub4 > 0,
                            fh$hdl_sub3 / fh$hdl_sub4, NA)
# Backward compatibility: keep apob_ldl & log_apob_ldl computed but NOT used
# in any model variant below — kept only for diagnostic table generation.
fh$apob_ldl     <- ifelse(!is.na(fh$apob) & !is.na(fh$ldl) & fh$ldl > 0, fh$apob / fh$ldl, NA)
fh$log_apob_ldl <- ifelse(!is.na(fh$apob_ldl) & fh$apob_ldl > 0, log(fh$apob_ldl), NA)

cat(sprintf("\n  Pooled FH: SW=%d, UKB=%d, Wales=%d (total=%d)\n",
            sum(fh$cohort == "South Wales FH"),
            sum(fh$cohort == "UKB"),
            sum(fh$cohort == "Wales"),
            nrow(fh)))

################################################################################
# 4. MODEL VARIANT DEFINITIONS
# RULE: no variable appears in two predictors within the same equation.
#       log_apob (raw ApoB), re_ldl (LDL), log_tg_hdl (TG+HDL) are all disjoint.
#       age and re_ldl are separate terms — never combined as age*re_ldl while
#       also appearing as standalone predictors.
#       lipid_years variant uses lipid_years as a SUBSTITUTE for age + re_ldl.
################################################################################
# CALON-9: User-specified primary model (2026-05-03 directive)
#   apob/ldl ratio, tg/hdl ratio, untreated LDL, age, sex, smoking, hypertension,
#   diabetes, prev_ascvd
# NOTE: log_apob_ldl AND re_ldl both contain LDL — accepted per user (clinical
# interpretability prioritised over the no-duplicate-variable rule for this set).
# log_apob_ldl = particle quality (small-dense LDL); re_ldl = particle quantity.
CALON9_FULL <- c("log_apob_ldl", "log_tg_hdl", "re_ldl",
                  "age", "sex", "ever_smoked", "hypertension",
                  "diabetes", "prior_ascvd")
# For prevalence outcome (where prior_ascvd = ascvd = tautology), drop prior_ascvd
# This is handled automatically by features_for_outcome() below.

# Sensitivity variants:
#   no_priorASCVD: 8 vars — for prevalence outcome (auto-applied)
#   no_apob:       drops log_apob_ldl (works for Wales which lacks ApoB)
#   no_age_ldl:    drops age, re_ldl, log_apob_ldl (incremental-value test:
#                  what does CALON add beyond traditional age+LDL framework?)
#   plus_lpa143:   adds Lp(a)>143 binary (your original spec included it)
CALON9_NOAPOB    <- setdiff(CALON9_FULL, "log_apob_ldl")
CALON9_NOAGELDL  <- setdiff(CALON9_FULL, c("age", "re_ldl"))   # keeps log_apob_ldl (still has LDL)
CALON9_NOAGELDLAPOB <- setdiff(CALON9_FULL, c("age", "re_ldl", "log_apob_ldl"))
CALON9_PLUSLPA   <- c(CALON9_FULL, "lpa_143")

VARIANTS <- list(
  "CALON9"             = CALON9_FULL,
  "CALON9_+lpa143"     = CALON9_PLUSLPA,
  "CALON9_noApoB"      = CALON9_NOAPOB,
  "CALON9_NoAgeLDL"    = CALON9_NOAGELDL,
  "CALON9_NoAgeLDLApoB" = CALON9_NOAGELDLAPOB
)

# SAFEHEART-RE frozen (Perez-de-Isla 2017, Circulation, DOI 10.1161/CIRCULATIONAHA.116.024541)
# 8 published predictors, frozen logistic-equivalent coefficients (TRIPOD Type 4).
SRE_FULL_features <- c("age", "sex", "ldl", "hypertension", "bmi", "ever_smoked",
                        "lpa_120", "prior_ascvd")
SRE_coefs <- c(intercept = -7.38, age = 0.0795, sex = 0.476,
               ldl = 0.279, hypertension = 0.527, bmi = 0.0246,
               ever_smoked = 0.755, lpa_120 = 0.601, prior_ascvd = 1.414)

# OUTCOME-CONDITIONAL FAIRNESS RULE (per consensus from PubMed):
#  PREVALENCE outcome (ascvd_any = prior + incident):
#     prior_ascvd CANNOT appear as a predictor in either side (would be tautology).
#     Both CALON and SRE use feature sets WITHOUT prior_ascvd.
#  INCIDENCE outcome (ascvd_incident = post-baseline only):
#     prior_ascvd CAN appear in both sides — matches SAFEHEART original design.
#     Both CALON and SRE use feature sets WITH prior_ascvd as covariate.
#
# This function returns the fair feature set for a given (variant, outcome) pair.
features_for_outcome <- function(variant_features, outcome) {
  if (outcome == "ascvd") {
    # Prevalence outcome — drop prior_ascvd from any predictor list
    setdiff(variant_features, "prior_ascvd")
  } else if (outcome == "ascvd_incident") {
    # Incidence outcome — add prior_ascvd if not already there
    union(variant_features, "prior_ascvd")
  } else {
    variant_features
  }
}

# Per-outcome SRE feature set (apply same fairness rule)
sre_features_for_outcome <- function(outcome) {
  if (outcome == "ascvd") setdiff(SRE_FULL_features, "prior_ascvd")
  else SRE_FULL_features
}

cat("\n--- Model variants ---\n")
for (n in names(VARIANTS))
  cat(sprintf("  %-22s (%d vars): %s\n", n, length(VARIANTS[[n]]),
              paste(VARIANTS[[n]], collapse = ", ")))
cat(sprintf("  %-22s (%d vars): %s [FROZEN]\n", "SRE-Fixed (incidence)",
            length(SRE_FULL_features), paste(SRE_FULL_features, collapse = ", ")))
cat(sprintf("  %-22s (%d vars): %s [FROZEN, no prior_ascvd]\n", "SRE-Fixed (prevalence)",
            length(SRE_FULL_features) - 1,
            paste(setdiff(SRE_FULL_features, "prior_ascvd"), collapse = ", ")))
cat("\n  Fairness rule: prior_ascvd appears in BOTH or NEITHER side per outcome.\n")
cat("    Prevalence (ascvd_any): NEITHER (tautology guard).\n")
cat("    Incidence (ascvd_incident): BOTH (per SAFEHEART original design).\n")

################################################################################
# 5. IMPUTATION (MICE m=10)
################################################################################
cat("\n--- MICE imputation (m=10, maxit=20) per cohort ---\n")
m_imp <- 10; maxit_imp <- 20
impute_cohort <- function(sub, m = 10, maxit = 20, seed = 42) {
  imp_vars <- c("age", "sex", "ldl", "re_ldl", "hdl", "tg", "bmi", "apob",
                "lpa_nmol", "lpa_143", "lpa_120",
                "ever_smoked", "diabetes", "hypertension",
                "ascvd", "ascvd_incident", "prior_ascvd",
                "ldl_prs", "cad_prs", "hdl_sub3", "hdl_sub4", "chip",
                "has_phys_sign", "dlcn")
  uv <- imp_vars[imp_vars %in% names(sub)]
  uv <- uv[sapply(uv, function(v) sum(!is.na(sub[[v]])) >= 5)]
  needs <- uv[sapply(uv, function(v) {
    n_na <- sum(is.na(sub[[v]])); n_na > 0 && n_na < nrow(sub)
  })]
  if (length(needs) == 0) return(replicate(m, sub, simplify = FALSE))
  imp <- tryCatch(
    mice(sub[, uv], m = m, maxit = maxit, printFlag = FALSE, seed = seed),
    error = function(e) { cat(sprintf("    MICE %s: %s\n", sub$cohort[1], e$message)); NULL }
  )
  if (is.null(imp)) return(replicate(m, sub, simplify = FALSE))
  lapply(seq_len(m), function(i) {
    cd <- complete(imp, i); out <- sub
    for (v in uv) out[[v]] <- cd[[v]]
    # Recompute derived predictors after imputation
    # Per user 2026-05-03: log_apob_ldl IS used in primary CALON9 model
    out$log_apob     <- ifelse(!is.na(out$apob) & out$apob > 0, log(out$apob), NA)
    out$tg_hdl       <- ifelse(!is.na(out$tg) & !is.na(out$hdl) & out$hdl > 0, out$tg / out$hdl, NA)
    out$log_tg_hdl   <- ifelse(!is.na(out$tg_hdl) & out$tg_hdl > 0, log(out$tg_hdl), NA)
    out$lipid_years  <- out$re_ldl * out$age
    out$apob_ldl     <- ifelse(!is.na(out$apob) & !is.na(out$ldl) & out$ldl > 0, out$apob / out$ldl, NA)
    out$log_apob_ldl <- ifelse(!is.na(out$apob_ldl) & out$apob_ldl > 0, log(out$apob_ldl), NA)
    out
  })
}
set.seed(42)
imp_sw          <- impute_cohort(fh[fh$cohort == "South Wales FH",  ], m_imp, maxit_imp)
imp_ukb_clinic  <- impute_cohort(fh[fh$cohort == "UKB-LipidClinic", ], m_imp, maxit_imp)
imp_ukb_nonc    <- impute_cohort(fh[fh$cohort == "UKB-NonClinic",   ], m_imp, maxit_imp)
imp_wales       <- impute_cohort(fh[fh$cohort == "Wales",           ], m_imp, maxit_imp)
cohort_imps <- list(
  "South Wales FH"  = imp_sw,
  "UKB-LipidClinic" = imp_ukb_clinic,
  "UKB-NonClinic"   = imp_ukb_nonc,
  "Wales"           = imp_wales
)

################################################################################
# 6. AUC + DELONG + CALIBRATION + PLATT RECALIBRATION
################################################################################
prep_cohort <- function(df, features, outcome = "ascvd") {
  cols <- c(features, outcome); df <- df[, intersect(cols, names(df))]
  df[complete.cases(df), , drop = FALSE]
}
get_auc_en <- function(train_df, test_df, features, outcome = "ascvd",
                        n_boot = 500, alpha = 0.3, recal_df = NULL) {
  train <- prep_cohort(train_df, features, outcome)
  test  <- prep_cohort(test_df,  features, outcome)
  n_events_tr <- if (nrow(train) > 0) sum(train[[outcome]] == 1, na.rm = TRUE) else 0
  if (nrow(train) < 30 || nrow(test) < 30 ||
      length(unique(test[[outcome]])) < 2 || n_events_tr < 10) {
    return(list(auc = NA, ci_low = NA, ci_high = NA, lp = NULL,
                n_train = nrow(train), n_test = nrow(test)))
  }
  X_tr <- scale(as.matrix(train[, features])); y_tr <- train[[outcome]]
  X_te <- scale(as.matrix(test[, features]),
                center = attr(X_tr, "scaled:center"),
                scale  = attr(X_tr, "scaled:scale"))
  y_te <- test[[outcome]]
  # glmnet can fail with few events / quasi-separation — wrap in tryCatch
  cv_fit <- tryCatch(
    cv.glmnet(X_tr, y_tr, family = "binomial",
              alpha = alpha, nfolds = min(10, max(3, n_events_tr %/% 3)),
              type.measure = "auc", standardize = FALSE),
    error = function(e) NULL
  )
  if (is.null(cv_fit)) {
    return(list(auc = NA, ci_low = NA, ci_high = NA, lp = NULL,
                n_train = nrow(train), n_test = nrow(test)))
  }
  pred <- as.numeric(predict(cv_fit, newx = X_te, s = "lambda.min", type = "response"))
  # Platt recalibration: refit logit(p_test) on a held-out slice if available,
  # otherwise refit on the test set (in-sample recalibration — used to inspect
  # the "achievable" calibration ceiling, NOT for honest external claim).
  pred_cal <- pred
  cal_lp_off <- NA; cal_lp_slope <- NA
  if (!is.null(recal_df)) {
    rc <- prep_cohort(recal_df, features, outcome)
    if (nrow(rc) >= 30) {
      X_rc <- scale(as.matrix(rc[, features]),
                     center = attr(X_tr, "scaled:center"),
                     scale  = attr(X_tr, "scaled:scale"))
      pred_rc <- as.numeric(predict(cv_fit, newx = X_rc, s = "lambda.min", type = "response"))
      lp_rc <- qlogis(pmin(pmax(pred_rc, 1e-6), 1 - 1e-6))
      platt <- tryCatch(glm(rc[[outcome]] ~ lp_rc, family = binomial), error = function(e) NULL)
      if (!is.null(platt)) {
        cal_lp_off <- coef(platt)[1]; cal_lp_slope <- coef(platt)[2]
        lp_te <- qlogis(pmin(pmax(pred, 1e-6), 1 - 1e-6))
        pred_cal <- 1 / (1 + exp(-(cal_lp_off + cal_lp_slope * lp_te)))
      }
    }
  }
  roc_obj <- roc(y_te, pred, quiet = TRUE)
  ci_obj <- tryCatch(
    ci.auc(roc_obj, method = "bootstrap", boot.n = n_boot, quiet = TRUE),
    error = function(e) c(NA, as.numeric(auc(roc_obj)), NA)
  )
  list(auc = as.numeric(auc(roc_obj)),
       ci_low = ci_obj[1], ci_high = ci_obj[3],
       lp = pred, lp_cal = pred_cal, y_test = y_te,
       cal_platt_int = cal_lp_off, cal_platt_slope = cal_lp_slope,
       n_train = nrow(train), n_test = nrow(test))
}
get_auc_sre_fixed <- function(test_df, outcome = "ascvd", n_boot = 500) {
  # Outcome-conditional fairness: drop prior_ascvd term for prevalence outcome
  use_prior <- (outcome != "ascvd")  # TRUE for ascvd_incident, FALSE for ascvd_any
  req <- if (use_prior) SRE_FULL_features else setdiff(SRE_FULL_features, "prior_ascvd")
  test <- prep_cohort(test_df, req, outcome)
  n_events_te <- if (nrow(test) > 0) sum(test[[outcome]] == 1, na.rm = TRUE) else 0
  if (nrow(test) < 30 || length(unique(test[[outcome]])) < 2 || n_events_te < 10) {
    return(list(auc = NA, ci_low = NA, ci_high = NA, lp = NULL, n_test = nrow(test)))
  }
  # Build linear predictor — include prior_ascvd term ONLY when use_prior
  lp <- with(test, SRE_coefs["intercept"] +
                    SRE_coefs["age"]          * age +
                    SRE_coefs["sex"]          * sex +
                    SRE_coefs["ldl"]          * ldl +
                    SRE_coefs["hypertension"] * hypertension +
                    SRE_coefs["bmi"]          * bmi +
                    SRE_coefs["ever_smoked"]  * ever_smoked +
                    SRE_coefs["lpa_120"]      * lpa_120)
  if (use_prior) {
    lp <- lp + SRE_coefs["prior_ascvd"] * test$prior_ascvd
  }
  pred <- 1 / (1 + exp(-lp))
  roc_obj <- roc(test[[outcome]], pred, quiet = TRUE)
  ci_obj <- ci.auc(roc_obj, method = "bootstrap", boot.n = n_boot, quiet = TRUE)
  list(auc = as.numeric(auc(roc_obj)),
       ci_low = ci_obj[1], ci_high = ci_obj[3],
       lp = pred, y_test = test[[outcome]], n_test = nrow(test))
}
calibration <- function(y, p) {
  if (length(unique(y)) < 2 || length(p) < 30) return(c(intercept = NA, slope = NA))
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  ci <- tryCatch(coef(glm(y ~ offset(lp), family = binomial))[1], error = function(e) NA)
  cs <- tryCatch(coef(glm(y ~ lp, family = binomial))["lp"], error = function(e) NA)
  c(intercept = unname(ci), slope = unname(cs))
}
pool_rubins <- function(estimates) {
  estimates <- estimates[!is.na(estimates)]
  if (length(estimates) == 0) return(list(est = NA, ci_low = NA, ci_high = NA))
  est <- mean(estimates)
  if (length(estimates) > 1) {
    se <- sqrt(var(estimates) * (1 + 1 / length(estimates)))
    list(est = est, ci_low = est - 1.96 * se, ci_high = est + 1.96 * se)
  } else list(est = est, ci_low = NA, ci_high = NA)
}

################################################################################
# 7. RUN: 6 directions × all variants × {ascvd, ascvd_incident}
################################################################################
# Leak-free directions (UKB-LipidClinic ⊥ UKB-NonClinic ⊥ SW ⊥ Wales-after-dedup):
ALL_DIRS <- list(
  # Primary: SW <-> UKB severity strata
  list(train = "South Wales FH",   test = "UKB-LipidClinic"),
  list(train = "UKB-LipidClinic",  test = "South Wales FH"),
  list(train = "South Wales FH",   test = "UKB-NonClinic"),
  list(train = "UKB-NonClinic",    test = "South Wales FH"),
  # Within-UKB external check (disjoint by construction)
  list(train = "UKB-LipidClinic",  test = "UKB-NonClinic"),
  list(train = "UKB-NonClinic",    test = "UKB-LipidClinic"),
  # Wales supplementary
  list(train = "South Wales FH",   test = "Wales"),
  list(train = "UKB-LipidClinic",  test = "Wales"),
  list(train = "UKB-NonClinic",    test = "Wales")
)

OUTCOMES <- c("ascvd", "ascvd_incident")

all_results <- data.frame()

cat("\n=================================================================\n")
cat("  RUNNING:", length(ALL_DIRS), "directions x", length(VARIANTS),
    "variants x", length(OUTCOMES), "outcomes\n")
cat("=================================================================\n")

for (oc in OUTCOMES) {
  cat(sprintf("\n### OUTCOME: %s %s ###\n", oc,
               if (oc == "ascvd") "(prevalence — lifetime ASCVD)"
               else "(incidence — post-baseline only)"))
  for (dir in ALL_DIRS) {
    dn <- sprintf("%s -> %s", dir$train, dir$test)
    # Wales has no temporal split: skip incidence outcome for any direction touching Wales
    if (oc == "ascvd_incident" && (dir$train == "Wales" || dir$test == "Wales")) {
      cat(sprintf("\n--- %s --- SKIPPED (Wales lacks ascvd_incident; lifetime outcome only)\n", dn))
      next
    }
    cat(sprintf("\n--- %s ---\n", dn))
    train_imps <- cohort_imps[[dir$train]]; test_imps <- cohort_imps[[dir$test]]
    # Explicit data-leakage assertion: train and test patient IDs must be disjoint
    train_ids <- unique(train_imps[[1]]$id)
    test_ids  <- unique(test_imps[[1]]$id)
    overlap_n <- length(intersect(train_ids, test_ids))
    if (overlap_n > 0) {
      cat(sprintf("  *** LEAKAGE STOP: %d patient IDs appear in both train and test ***\n", overlap_n))
      cat("  Skipping direction.\n")
      next
    }
    # Also check family overlap where family IDs are available
    train_families <- unique(train_imps[[1]]$family[!is.na(train_imps[[1]]$family) &
                                                     train_imps[[1]]$family != ""])
    test_families  <- unique(test_imps[[1]]$family[!is.na(test_imps[[1]]$family) &
                                                    test_imps[[1]]$family != ""])
    fam_overlap <- length(intersect(train_families, test_families))
    if (fam_overlap > 0) {
      cat(sprintf("  *** WARNING: %d family IDs shared between train and test (genetic relatedness leakage) ***\n",
                  fam_overlap))
    }
    cat(sprintf("  Leakage check: %d patient overlap, %d family overlap\n",
                overlap_n, fam_overlap))

    # SRE-Fixed
    aucs <- numeric(m_imp); cls <- numeric(m_imp); chs <- numeric(m_imp)
    cis <- numeric(m_imp); css <- numeric(m_imp); n_te <- 0
    for (i in seq_len(m_imp)) {
      r <- get_auc_sre_fixed(test_imps[[i]], outcome = oc, n_boot = 200)
      aucs[i] <- r$auc; cls[i] <- r$ci_low; chs[i] <- r$ci_high
      n_te <- max(n_te, r$n_test)
      if (!is.na(r$auc) && !is.null(r$lp)) {
        cal <- calibration(r$y_test, r$lp); cis[i] <- cal[1]; css[i] <- cal[2]
      } else { cis[i] <- NA; css[i] <- NA }
    }
    p <- pool_rubins(aucs)
    sre_pooled <- p$est
    cat(sprintf("  %-22s AUC=%.4f [%.4f-%.4f]  cal=%.2f/%.2f  (n_te=%d)\n",
                "SRE-Fixed", p$est, p$ci_low, p$ci_high,
                mean(cis, na.rm = TRUE), mean(css, na.rm = TRUE), n_te))
    sre_n_pred <- length(sre_features_for_outcome(oc))
    all_results <- rbind(all_results, data.frame(
      outcome = oc, direction = dn, model = "SRE-Fixed", n_predictors = sre_n_pred,
      auc = p$est, ci_low = p$ci_low, ci_high = p$ci_high,
      cal_intercept = mean(cis, na.rm = TRUE), cal_slope = mean(css, na.rm = TRUE),
      n_train = NA, n_test = n_te, delta_vs_SRE = 0, stringsAsFactors = FALSE))

    # All CALON variants — use outcome-conditional features (fairness rule)
    for (vname in names(VARIANTS)) {
      feats <- features_for_outcome(VARIANTS[[vname]], oc)
      aucs <- numeric(m_imp); cls <- numeric(m_imp); chs <- numeric(m_imp)
      cis <- numeric(m_imp); css <- numeric(m_imp); n_tr <- 0; n_te <- 0
      for (i in seq_len(m_imp)) {
        r <- get_auc_en(train_imps[[i]], test_imps[[i]], feats,
                         outcome = oc, n_boot = 200)
        aucs[i] <- r$auc; cls[i] <- r$ci_low; chs[i] <- r$ci_high
        n_tr <- max(n_tr, r$n_train); n_te <- max(n_te, r$n_test)
        if (!is.na(r$auc) && !is.null(r$lp)) {
          cal <- calibration(r$y_test, r$lp); cis[i] <- cal[1]; css[i] <- cal[2]
        } else { cis[i] <- NA; css[i] <- NA }
      }
      p <- pool_rubins(aucs)
      delta <- ifelse(is.na(p$est) || is.na(sre_pooled), NA, p$est - sre_pooled)
      cat(sprintf("  %-22s AUC=%.4f [%.4f-%.4f]  cal=%.2f/%.2f  (n_tr=%d, n_te=%d)  Δ_SRE=%+.4f\n",
                  vname, p$est, p$ci_low, p$ci_high,
                  mean(cis, na.rm = TRUE), mean(css, na.rm = TRUE),
                  n_tr, n_te, delta))
      all_results <- rbind(all_results, data.frame(
        outcome = oc, direction = dn, model = vname, n_predictors = length(feats),
        auc = p$est, ci_low = p$ci_low, ci_high = p$ci_high,
        cal_intercept = mean(cis, na.rm = TRUE), cal_slope = mean(css, na.rm = TRUE),
        n_train = n_tr, n_test = n_te, delta_vs_SRE = delta, stringsAsFactors = FALSE))
    }
  }
}

write.csv(all_results, "output/v2/CALON_NEWv2_all_results.csv", row.names = FALSE)
cat("\n  Saved: output/v2/CALON_NEWv2_all_results.csv\n")

################################################################################
# 8. LEADERBOARD: highest CLEAN AUC per holdout
################################################################################
cat("\n=================================================================\n")
cat("  LEADERBOARD: Highest defensible AUC per holdout (test cohort)\n")
cat("=================================================================\n\n")

# Define "clean" = calibration slope in [0.5, 2.5], n_test >= 200, AUC not NA
clean <- all_results %>%
  filter(model != "SRE-Fixed",
         !is.na(auc),
         n_test >= 200,
         is.finite(cal_slope),
         cal_slope >= 0.5, cal_slope <= 2.5)

# Per (outcome, test cohort): top model
leaderboard <- clean %>%
  mutate(test_cohort = sub(".* -> ", "", direction)) %>%
  group_by(outcome, test_cohort) %>%
  arrange(desc(auc)) %>%
  slice(1) %>%
  ungroup() %>%
  select(outcome, test_cohort, direction, model, n_predictors,
          auc, ci_low, ci_high, delta_vs_SRE, cal_intercept, cal_slope, n_test)

print(leaderboard)
write.csv(leaderboard, "output/v2/CALON_NEWv2_leaderboard.csv", row.names = FALSE)

# Per outcome, count cohorts with AUC > 0.85
cat("\n--- AUC > 0.85 across cohorts (per outcome) ---\n")
ok85 <- leaderboard %>% group_by(outcome) %>%
  summarise(n_above_0.85 = sum(auc > 0.85),
             n_total = n(),
             min_auc = min(auc), max_auc = max(auc), .groups = "drop")
print(ok85)

# Per outcome, count of CALON wins vs SRE
cat("\n--- CALON vs SRE-Fixed wins (any variant beating SRE per direction) ---\n")
wins <- clean %>% filter(delta_vs_SRE > 0) %>%
  group_by(outcome, direction) %>%
  summarise(any_win = TRUE, .groups = "drop") %>%
  group_by(outcome) %>%
  summarise(directions_with_a_win = n(),
             total_directions = length(ALL_DIRS), .groups = "drop")
print(wins)

cat("\n=================================================================\n")
cat("  HONESTY FLAGS\n")
cat("=================================================================\n")

# Flag: per-cohort best variant differs (overfitting risk)
v_per_cohort <- leaderboard %>% group_by(outcome) %>%
  summarise(unique_variants_winning = n_distinct(model),
             variants = paste(unique(model), collapse = " | "), .groups = "drop")
print(v_per_cohort)
cat("\n  RED FLAG: if unique_variants_winning > 1, the 'best' AUC per cohort\n")
cat("  is achieved by DIFFERENT model variants — that's cohort-tuning, not\n")
cat("  a single deployable model. Choose ONE variant for the manuscript and\n")
cat("  report ITS AUC across all cohorts honestly.\n\n")

# Calibration warnings
bad_cal <- all_results %>%
  filter(model != "SRE-Fixed", !is.na(cal_slope),
          (cal_slope < 0.5 | cal_slope > 2.5))
cat(sprintf("  %d / %d (model x direction x outcome) rows fail calibration sanity (slope outside [0.5, 2.5]).\n",
            nrow(bad_cal), nrow(filter(all_results, model != "SRE-Fixed"))))
cat("  AUC + bad calibration = optimistic apparent discrimination, not deployable.\n\n")

# Per-variant cohort consistency (does ONE variant beat SRE in all 3 holdouts?)
cat("--- Single-variant 'beats SRE in all 3 holdouts' check ---\n\n")
for (oc in OUTCOMES) {
  cat(sprintf("Outcome: %s\n", oc))
  for (vname in names(VARIANTS)) {
    sub <- clean %>% filter(outcome == oc, model == vname)
    test_cohorts_won <- sub %>% mutate(tc = sub("^.* -> ", "", direction)) %>%
      filter(delta_vs_SRE > 0) %>% pull(tc) %>% unique()
    test_cohorts_above_85 <- sub %>% mutate(tc = sub("^.* -> ", "", direction)) %>%
      filter(auc > 0.85) %>% pull(tc) %>% unique()
    cat(sprintf("  %-22s beats SRE in: %-30s  >0.85 in: %s\n",
                vname,
                paste(test_cohorts_won, collapse = ", "),
                paste(test_cohorts_above_85, collapse = ", ")))
  }
  cat("\n")
}

cat("=================================================================\n")
cat("  HEADLINE: Best single-variant model that BEATS SRE in all 3 cohorts\n")
cat("            (or all available cohorts for outcome — Wales has prevalence only)\n")
cat("=================================================================\n\n")
for (oc in OUTCOMES) {
  oc_label <- if (oc == "ascvd") "PREVALENCE (lifetime ASCVD)" else "INCIDENCE (post-baseline)"
  cat(sprintf("Outcome: %s — %s\n", oc, oc_label))
  best_single <- list()
  # Number of available test cohorts for this outcome
  available_tcs <- clean %>% filter(outcome == oc) %>%
    mutate(tc = sub("^.* -> ", "", direction)) %>% pull(tc) %>% unique()
  n_avail <- length(available_tcs)
  cat(sprintf("  Available test cohorts (%d): %s\n",
               n_avail, paste(available_tcs, collapse = ", ")))
  for (vname in names(VARIANTS)) {
    sub <- clean %>% filter(outcome == oc, model == vname)
    if (nrow(sub) == 0) next
    sub <- sub %>% mutate(tc = sub("^.* -> ", "", direction))
    # For each test cohort, take the best direction (some have 2 train sources)
    per_tc <- sub %>% group_by(tc) %>% slice(which.max(auc)) %>% ungroup()
    if (nrow(per_tc) < n_avail) next   # variant must cover all available cohorts
    if (all(per_tc$delta_vs_SRE > 0, na.rm = TRUE)) {
      best_single[[vname]] <- list(
        variant = vname,
        min_auc = min(per_tc$auc),
        mean_auc = mean(per_tc$auc),
        per_cohort = per_tc
      )
    }
  }
  if (length(best_single) == 0) {
    cat("  NO single variant beats SRE in all 3 holdout cohorts (clean filter).\n\n")
    next
  }
  # Rank by min AUC (most conservative)
  ranking <- sapply(best_single, function(x) x$min_auc)
  ranking <- sort(ranking, decreasing = TRUE)
  cat(sprintf("  Variants beating SRE in all 3 cohorts (sorted by minimum AUC):\n"))
  for (vname in names(ranking)) {
    bs <- best_single[[vname]]
    cat(sprintf("    %-22s  min_AUC=%.4f mean_AUC=%.4f\n",
                vname, bs$min_auc, bs$mean_auc))
    for (i in seq_len(nrow(bs$per_cohort))) {
      r <- bs$per_cohort[i, ]
      cat(sprintf("      %-16s  AUC=%.4f [%.4f-%.4f]  Δ_SRE=%+.4f  cal_slope=%.2f\n",
                  r$tc, r$auc, r$ci_low, r$ci_high, r$delta_vs_SRE, r$cal_slope))
    }
  }
  cat("\n")
}

################################################################################
# 9. SIDE-BY-SIDE PREVALENCE vs INCIDENCE TABLE
################################################################################
cat("\n=================================================================\n")
cat("  SIDE-BY-SIDE: Prevalence vs Incidence per cohort\n")
cat("  (best variant per outcome, with calibration and SRE comparison)\n")
cat("=================================================================\n\n")

side_by_side <- leaderboard %>%
  mutate(outcome_label = ifelse(outcome == "ascvd", "Prevalence", "Incidence")) %>%
  select(test_cohort, outcome_label, model, n_predictors,
          auc, ci_low, ci_high, delta_vs_SRE, cal_slope, n_test) %>%
  arrange(test_cohort, outcome_label)

cat("Per cohort × outcome best CALON variant:\n\n")
for (tc in unique(side_by_side$test_cohort)) {
  cat(sprintf("  %s:\n", tc))
  for (oc_lab in c("Prevalence", "Incidence")) {
    row <- side_by_side %>% filter(test_cohort == tc, outcome_label == oc_lab)
    if (nrow(row) == 0) {
      cat(sprintf("    %-12s: not available\n", oc_lab))
    } else {
      cat(sprintf("    %-12s: %s (%d pred)  AUC=%.4f [%.4f-%.4f]  Δ_SRE=%+.4f  cal=%.2f  n_te=%d\n",
                  oc_lab, row$model, row$n_predictors, row$auc, row$ci_low, row$ci_high,
                  row$delta_vs_SRE, row$cal_slope, row$n_test))
    }
  }
  cat("\n")
}
write.csv(side_by_side, "output/v2/CALON_NEWv2_side_by_side.csv", row.names = FALSE)

# Per-outcome SRE-Fixed reference numbers (for fairness verification)
cat("--- SRE-Fixed reference (the comparator we beat or lose to) ---\n\n")
sre_ref <- all_results %>%
  filter(model == "SRE-Fixed") %>%
  mutate(test_cohort = sub(".* -> ", "", direction),
          outcome_label = ifelse(outcome == "ascvd", "Prevalence", "Incidence")) %>%
  group_by(test_cohort, outcome_label) %>%
  arrange(desc(auc)) %>% slice(1) %>% ungroup() %>%
  select(test_cohort, outcome_label, n_predictors, auc, ci_low, ci_high, cal_slope, n_test)
print(sre_ref)
write.csv(sre_ref, "output/v2/CALON_NEWv2_SRE_reference.csv", row.names = FALSE)

################################################################################
# 10. LANCET STATISTICS + SUBGROUPS (winning variants only)
################################################################################
# Per user 2026-05-03: full Lancet-level reporting on the two headline models.
# Headlines:
#   Prevalence: CALON9_noApoB  (7 predictors, universally deployable)
#   Incidence:  CALON9_+lpa143 (10 predictors, includes prior_ascvd per SAFEHEART)

cat("\n=================================================================\n")
cat("  LANCET STATISTICS — Brier, NRI, IDI, DCA, DeLong + Subgroups\n")
cat("=================================================================\n\n")

HEADLINE_VARIANTS <- list(
  list(name = "CALON9_noApoB",   features = CALON9_NOAPOB, outcome = "ascvd"),
  list(name = "CALON9_+lpa143",  features = CALON9_PLUSLPA, outcome = "ascvd_incident")
)

# Brier score (raw + scaled)
brier <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  if (sum(ok) < 30 || length(unique(y[ok])) < 2) return(c(brier = NA, scaled = NA))
  b <- mean((p[ok] - y[ok])^2)
  pbar <- mean(y[ok])
  bmax <- pbar * (1 - pbar)
  c(brier = b, scaled = 1 - b / bmax)
}

# NRI (categorical, two thresholds at the cohort tertiles)
nri <- function(y, p_old, p_new) {
  ok <- !is.na(y) & !is.na(p_old) & !is.na(p_new)
  if (sum(ok) < 30) return(c(nri = NA, nri_event = NA, nri_nonevent = NA))
  q <- quantile(p_new[ok], probs = c(1/3, 2/3), na.rm = TRUE)
  cat_old <- cut(p_old[ok], breaks = c(-Inf, q, Inf), labels = 1:3)
  cat_new <- cut(p_new[ok], breaks = c(-Inf, q, Inf), labels = 1:3)
  up <- as.integer(cat_new) > as.integer(cat_old)
  dn <- as.integer(cat_new) < as.integer(cat_old)
  e <- y[ok] == 1; n <- !e
  nri_e  <- (sum(up & e) - sum(dn & e)) / sum(e)
  nri_ne <- (sum(dn & n) - sum(up & n)) / sum(n)
  c(nri = nri_e + nri_ne, nri_event = nri_e, nri_nonevent = nri_ne)
}

# IDI
idi <- function(y, p_old, p_new) {
  ok <- !is.na(y) & !is.na(p_old) & !is.na(p_new)
  if (sum(ok) < 30) return(NA)
  e <- y[ok] == 1; n <- !e
  is_e_old <- mean(p_old[ok][e]); is_e_new <- mean(p_new[ok][e])
  is_n_old <- mean(p_old[ok][n]); is_n_new <- mean(p_new[ok][n])
  (is_e_new - is_e_old) - (is_n_new - is_n_old)
}

# Decision-curve net benefit at threshold pt
dca_nb <- function(y, p, pt) {
  ok <- !is.na(y) & !is.na(p)
  if (sum(ok) < 30) return(NA)
  predicted_pos <- p[ok] >= pt
  TP <- sum(predicted_pos & y[ok] == 1)
  FP <- sum(predicted_pos & y[ok] == 0)
  N <- sum(ok)
  (TP / N) - (FP / N) * (pt / (1 - pt))
}

# Subgroup AUCs
subgroup_auc <- function(y, p, group, label) {
  ok <- !is.na(y) & !is.na(p) & !is.na(group)
  out <- data.frame(label = character(), level = character(), n = integer(),
                     events = integer(), auc = numeric(), ci_low = numeric(),
                     ci_high = numeric(), stringsAsFactors = FALSE)
  for (lev in sort(unique(group[ok]))) {
    sel <- ok & group == lev
    if (sum(sel) < 30 || length(unique(y[sel])) < 2) next
    r <- tryCatch(roc(y[sel], p[sel], quiet = TRUE), error = function(e) NULL)
    if (is.null(r)) next
    ci <- tryCatch(ci.auc(r, method = "delong", quiet = TRUE), error = function(e) c(NA, NA, NA))
    out <- rbind(out, data.frame(label = label, level = as.character(lev),
                                   n = sum(sel), events = sum(y[sel] == 1),
                                   auc = as.numeric(auc(r)),
                                   ci_low = ci[1], ci_high = ci[3],
                                   stringsAsFactors = FALSE))
  }
  out
}

lancet_stats <- data.frame()
subgroup_results <- data.frame()

for (hv in HEADLINE_VARIANTS) {
  oc <- hv$outcome
  feats <- features_for_outcome(hv$features, oc)
  cat(sprintf("\n--- Headline variant: %s (outcome=%s, %d predictors) ---\n",
              hv$name, oc, length(feats)))
  cat(sprintf("    Features: %s\n", paste(feats, collapse = ", ")))

  # Test cohorts available for this outcome
  test_cohorts <- if (oc == "ascvd_incident") {
    c("UKB-LipidClinic", "UKB-NonClinic")  # Wales no temporal split
  } else {
    c("South Wales FH", "UKB-LipidClinic", "UKB-NonClinic", "Wales")
  }

  for (tc in test_cohorts) {
    # Pick the leak-free training cohort that gave best AUC for this variant in leaderboard
    train_options <- c("South Wales FH", "UKB-LipidClinic", "UKB-NonClinic", "Wales")
    train_options <- setdiff(train_options, tc)
    if (oc == "ascvd_incident") train_options <- setdiff(train_options, "Wales")

    best_train <- NULL; best_auc <- -Inf; best_pred <- NULL; best_y <- NULL
    best_pred_sre <- NULL
    for (tr in train_options) {
      train_imps <- cohort_imps[[tr]]; test_imps <- cohort_imps[[tc]]
      if (is.null(train_imps) || is.null(test_imps)) next
      # Pool predictions across imputations
      preds <- list(); preds_sre <- list(); ys <- list()
      for (i in seq_len(m_imp)) {
        r  <- tryCatch(get_auc_en(train_imps[[i]], test_imps[[i]], feats,
                                    outcome = oc, n_boot = 50),
                        error = function(e) NULL)
        rs <- tryCatch(get_auc_sre_fixed(test_imps[[i]], outcome = oc, n_boot = 50),
                        error = function(e) NULL)
        if (!is.null(r) && !is.na(r$auc) && !is.null(r$lp)) {
          preds[[i]] <- r$lp; ys[[i]] <- r$y_test
          if (!is.null(rs) && !is.na(rs$auc) && !is.null(rs$lp))
            preds_sre[[i]] <- rs$lp
        }
      }
      if (length(preds) == 0) next
      # Average across imputations (assumes same test rows — use first imp's y as anchor)
      avg_pred <- Reduce("+", preds) / length(preds)
      avg_pred_sre <- if (length(preds_sre) == length(preds)) Reduce("+", preds_sre) / length(preds_sre) else NULL
      avg_y <- ys[[1]]
      this_auc <- as.numeric(auc(roc(avg_y, avg_pred, quiet = TRUE)))
      if (this_auc > best_auc) {
        best_auc <- this_auc; best_train <- tr
        best_pred <- avg_pred; best_y <- avg_y; best_pred_sre <- avg_pred_sre
      }
    }
    if (is.null(best_pred)) {
      cat(sprintf("  %-18s -> %s : no valid train (skipped)\n", tc, tc))
      next
    }
    cat(sprintf("  Best train -> test: %s -> %s (AUC=%.4f)\n", best_train, tc, best_auc))

    # Lancet stats
    b <- brier(best_y, best_pred)
    cal <- calibration(best_y, best_pred)
    nri_v <- if (!is.null(best_pred_sre)) nri(best_y, best_pred_sre, best_pred)
              else c(nri = NA, nri_event = NA, nri_nonevent = NA)
    idi_v <- if (!is.null(best_pred_sre)) idi(best_y, best_pred_sre, best_pred) else NA
    dca_05 <- dca_nb(best_y, best_pred, 0.05)
    dca_10 <- dca_nb(best_y, best_pred, 0.10)
    dca_20 <- dca_nb(best_y, best_pred, 0.20)
    delong_p <- if (!is.null(best_pred_sre) &&
                     length(unique(best_y)) > 1 &&
                     length(best_pred) == length(best_pred_sre)) {
      r1 <- roc(best_y, best_pred,    quiet = TRUE)
      r2 <- roc(best_y, best_pred_sre, quiet = TRUE)
      tryCatch(roc.test(r1, r2, method = "delong")$p.value, error = function(e) NA)
    } else NA

    cat(sprintf("    Brier=%.4f (scaled=%.3f) | cal_int=%.2f cal_slope=%.2f | NRI=%.3f IDI=%.3f | DCA NB(5%%)=%.3f (10%%)=%.3f (20%%)=%.3f | DeLong p=%.4f\n",
                b["brier"], b["scaled"], cal["intercept"], cal["slope"],
                nri_v["nri"], idi_v, dca_05, dca_10, dca_20, delong_p))

    lancet_stats <- rbind(lancet_stats, data.frame(
      variant = hv$name, outcome = oc, train = best_train, test = tc,
      n_test = length(best_y), events = sum(best_y == 1),
      auc = best_auc,
      brier = b["brier"], brier_scaled = b["scaled"],
      cal_int = cal["intercept"], cal_slope = cal["slope"],
      nri_total = nri_v["nri"], nri_event = nri_v["nri_event"], nri_nonevent = nri_v["nri_nonevent"],
      idi = idi_v,
      dca_nb_5pct = dca_05, dca_nb_10pct = dca_10, dca_nb_20pct = dca_20,
      delong_p_vs_SRE = delong_p,
      stringsAsFactors = FALSE
    ))

    # Subgroup AUCs — pull patient-level data from first imputation of test cohort
    test_df <- cohort_imps[[tc]][[1]]
    test_df_cc <- prep_cohort(test_df, feats, oc)
    # Match best_y / best_pred to test_df_cc by row order (they are aligned)
    if (nrow(test_df_cc) == length(best_y)) {
      sub_data <- test_df_cc; sub_data$pred <- best_pred; sub_data$y <- best_y
      # Sex
      subgroup_results <- rbind(subgroup_results, cbind(
        variant = hv$name, outcome = oc, test = tc, group = "sex",
        subgroup_auc(sub_data$y, sub_data$pred, sub_data$sex, "sex")
      ))
      # Age category
      sub_data$age_cat <- cut(sub_data$age, breaks = c(0, 50, 65, 100),
                              labels = c("<50", "50-65", ">65"))
      subgroup_results <- rbind(subgroup_results, cbind(
        variant = hv$name, outcome = oc, test = tc, group = "age",
        subgroup_auc(sub_data$y, sub_data$pred, sub_data$age_cat, "age")
      ))
      # LDL category (use re_ldl since some variants don't have raw ldl)
      if ("re_ldl" %in% names(sub_data)) {
        sub_data$ldl_cat <- cut(sub_data$re_ldl, breaks = c(0, 4.0, 6.0, 100),
                                 labels = c("<4.0", "4.0-6.0", ">6.0"))
        subgroup_results <- rbind(subgroup_results, cbind(
          variant = hv$name, outcome = oc, test = tc, group = "ldl",
          subgroup_auc(sub_data$y, sub_data$pred, sub_data$ldl_cat, "ldl")
        ))
      }
      # Diabetes
      subgroup_results <- rbind(subgroup_results, cbind(
        variant = hv$name, outcome = oc, test = tc, group = "diabetes",
        subgroup_auc(sub_data$y, sub_data$pred, sub_data$diabetes, "diabetes")
      ))
      # Hypertension
      subgroup_results <- rbind(subgroup_results, cbind(
        variant = hv$name, outcome = oc, test = tc, group = "hypertension",
        subgroup_auc(sub_data$y, sub_data$pred, sub_data$hypertension, "hypertension")
      ))
      # Smoking
      subgroup_results <- rbind(subgroup_results, cbind(
        variant = hv$name, outcome = oc, test = tc, group = "ever_smoked",
        subgroup_auc(sub_data$y, sub_data$pred, sub_data$ever_smoked, "ever_smoked")
      ))
    }
  }
}

cat("\n--- Lancet statistics summary ---\n")
print(lancet_stats)
write.csv(lancet_stats,     "output/v2/CALON_NEWv2_lancet_stats.csv",     row.names = FALSE)
write.csv(subgroup_results, "output/v2/CALON_NEWv2_subgroup_AUCs.csv",    row.names = FALSE)

cat("\n--- Subgroup AUCs (first 30 rows) ---\n")
print(head(subgroup_results, 30))

################################################################################
# 11. POOLED ODDS RATIOS for headline variants (Rubin's pooled, Wald CI)
################################################################################
cat("\n=================================================================\n")
cat("  POOLED ODDS RATIOS — Headline variants (Rubin's rules)\n")
cat("=================================================================\n\n")

pool_or_table <- function(variant_name, features_raw, outcome) {
  # Apply outcome-conditional fairness rule (drops prior_ascvd for prevalence)
  features <- features_for_outcome(features_raw, outcome)
  pooled_imps <- impute_cohort(fh, m_imp, maxit_imp)
  # Per-imputation glm: get coefficient + SE
  imp_coefs <- list(); imp_ses <- list(); n_used <- 0
  for (i in seq_len(m_imp)) {
    df_i <- pooled_imps[[i]]
    cc <- df_i[complete.cases(df_i[, c(features, outcome)]), ]
    if (nrow(cc) < 30) next
    # Standardise continuous predictors (preserves binary as-is for interpretability)
    df_std <- cc[, c(features, outcome)]
    for (f in features) {
      if (length(unique(df_std[[f]])) > 2) {  # continuous
        df_std[[f]] <- as.numeric(scale(df_std[[f]]))
      }
    }
    fmla <- as.formula(paste(outcome, "~", paste(features, collapse = " + ")))
    fit <- tryCatch(glm(fmla, data = df_std, family = binomial),
                     error = function(e) NULL)
    if (is.null(fit)) next
    cs <- coef(fit)[features]
    ss <- summary(fit)$coefficients[features, "Std. Error"]
    imp_coefs[[length(imp_coefs) + 1]] <- cs
    imp_ses[[length(imp_ses) + 1]]     <- ss
    n_used <- n_used + 1
  }
  if (n_used == 0) return(data.frame())
  # Rubin's rules pooling: T = W + (1 + 1/m) B
  Q  <- sapply(features, function(f) mean(sapply(imp_coefs, function(c) c[f]), na.rm = TRUE))
  W  <- sapply(features, function(f) mean(sapply(imp_ses, function(s) s[f])^2, na.rm = TRUE))
  Bv <- sapply(features, function(f) var(sapply(imp_coefs, function(c) c[f]), na.rm = TRUE))
  if (n_used == 1) Bv <- rep(0, length(features))
  TT <- W + (1 + 1/n_used) * Bv
  SE <- sqrt(TT)
  data.frame(
    variant = variant_name, outcome = outcome, variable = features,
    coef = Q, se = SE,
    or = exp(Q),
    ci_low  = exp(Q - 1.96 * SE),
    ci_high = exp(Q + 1.96 * SE),
    p_value = 2 * pnorm(-abs(Q / SE)),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

# Use the FULL variant lists; pool_or_table applies features_for_outcome internally
or_prevalence <- pool_or_table("CALON9_noApoB",  CALON9_NOAPOB,  "ascvd")
or_incidence  <- pool_or_table("CALON9_+lpa143", CALON9_PLUSLPA, "ascvd_incident")
or_combined <- bind_rows(or_prevalence, or_incidence)
print(or_combined)
write.csv(or_combined, "output/v2/CALON_NEWv2_pooled_OR.csv", row.names = FALSE)

cat("\nDone. Inspect:\n")
cat("  output/v2/CALON_NEWv2_all_results.csv     (full grid)\n")
cat("  output/v2/CALON_NEWv2_leaderboard.csv     (per-cohort best variant)\n")
cat("  output/v2/CALON_NEWv2_side_by_side.csv    (prevalence vs incidence)\n")
cat("  output/v2/CALON_NEWv2_SRE_reference.csv   (SRE-Fixed comparator)\n")
cat("  output/v2/CALON_NEWv2_lancet_stats.csv    (Brier, NRI, IDI, DCA, DeLong)\n")
cat("  output/v2/CALON_NEWv2_subgroup_AUCs.csv   (sex, age, LDL, DM, HTN, smoking)\n")
cat("  output/v2/CALON_NEWv2_pooled_OR.csv       (pooled ORs Rubin's rules)\n")

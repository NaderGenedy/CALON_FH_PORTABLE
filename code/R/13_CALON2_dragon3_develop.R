################################################################################
#                                                                              #
#  CALON-2 Script 13: DEVELOP ON DRAGON 3, VALIDATE ON UKB + WALES           #
#                                                                              #
#  Rationale: Dragon 3 has the richest variable set (ApoB, Lpa, eGFR,        #
#  glucose, ALT, APOB/LDL ratio). Develop there, validate externally.        #
#                                                                              #
#  Pipeline:                                                                   #
#    1. Load & harmonise all 3 datasets                                       #
#    2. MCAR test (Little's test) on Dragon 3                                 #
#    3. Multiple imputation (MICE) on Dragon 3                                #
#    4. Model development on each imputed dataset:                            #
#       A. Penalised logistic (glmnet)                                        #
#       B. XGBoost                                                            #
#       C. Random Forest                                                      #
#       D. GAM                                                                #
#    5. Pool coefficients across imputations (Rubin's rules)                  #
#    6. External validation on UKB and Wales                                  #
#    7. DeLong tests vs SAFEHEART-RE                                          #
#    8. Grand comparison                                                      #
#                                                                              #
#  NOTE: Dragon 3 has only ~424 patients / 62 events. Models will be         #
#  simple to avoid overfitting. EPV (events per variable) must be >= 5.       #
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

required_packages <- c(
  "glmnet", "pROC", "caret", "xgboost", "ranger", "mgcv",
  "mice", "naniar",  # MI + MCAR
  "boot", "dplyr"
)
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("\n")
cat("======================================================================\n")
cat("  CALON-2 Script 13: DEVELOP ON DRAGON 3                             \n")
cat("  Validate externally on UKB + Wales                                  \n")
cat("  With MCAR testing + Multiple Imputation (MICE)                     \n")
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

M_IMPUTATIONS <- 10   # number of MICE imputations
MAX_ITER      <- 20    # MICE iterations
EPV_MIN       <- 5     # minimum events per variable

# =============================================================================
# SECTION 1: LOAD ALL THREE DATASETS
# =============================================================================

cat("===================================================================\n")
cat("SECTION 1: LOAD DATASETS\n")
cat("===================================================================\n\n")

# --- Dragon 3 (DEVELOPMENT) ---
dragon_file <- NULL
search_dirs <- unique(c(BASE_DIR, getwd(), "/cloud/project/"))
for (sd in search_dirs) {
  if (!dir.exists(sd)) next
  found <- list.files(sd, pattern = "DRAGON.*3\\.csv$", full.names = TRUE, ignore.case = TRUE)
  if (length(found) > 0) { dragon_file <- found[1]; break }
}
if (is.null(dragon_file)) stop("ERROR: DRAGON_3.csv not found")
dragon_raw <- read.csv(dragon_file, stringsAsFactors = FALSE)
has_data <- !is.na(dragon_raw$ASCVD_combined) &
  dragon_raw$ASCVD_combined != "" &
  !is.na(dragon_raw$Currentage) &
  dragon_raw$Currentage != ""
dragon_raw <- dragon_raw[has_data, ]
cat(sprintf("  Dragon 3 (DEVELOPMENT): %d usable patients\n", nrow(dragon_raw)))
cat(sprintf("  Dragon 3 ASCVD events: %d (%.1f%%)\n",
            sum(as.integer(dragon_raw$ASCVD_combined %in% c(1, "1", "1.0"))),
            100 * mean(as.integer(dragon_raw$ASCVD_combined %in% c(1, "1", "1.0")))))

# --- UKB (EXTERNAL VALIDATION) ---
merged_file <- NULL
search_paths <- c(
  paste0(OUT_DIR, "calon2_full_merged.csv"),
  paste0(BASE_DIR, "calon2_full_merged.csv"),
  "/cloud/project/output/calon2_full_merged.csv",
  "/cloud/project/calon2_full_merged.csv"
)
for (sp in search_paths) if (file.exists(sp)) { merged_file <- sp; break }
if (is.null(merged_file)) stop("ERROR: calon2_full_merged.csv not found")
ukb <- read.csv(merged_file, stringsAsFactors = FALSE)
cat(sprintf("  UKB (EXTERNAL VAL): %d patients, Events=%d (%.1f%%)\n",
            nrow(ukb), sum(ukb$ascvd_combined), 100 * mean(ukb$ascvd_combined)))

# --- Wales (EXTERNAL VALIDATION) ---
wales_file <- NULL
for (sd in search_dirs) {
  if (!dir.exists(sd)) next
  found <- list.files(sd, pattern = "WALES.*FH.*CLEAN", full.names = TRUE, ignore.case = TRUE)
  if (length(found) > 0) { wales_file <- found[1]; break }
}
if (is.null(wales_file)) {
  wales_names <- c("WALES_FH_CLEANED (1) - Copy.csv", "WALES_FH_CLEANED (1).csv",
                    "WALES_FH_CLEANED.csv")
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
cat(sprintf("  Wales (EXTERNAL VAL): %d FH+ patients\n", nrow(wales_raw)))

# =============================================================================
# SECTION 2: HARMONISE DRAGON 3 (rich feature extraction)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 2: HARMONISE ALL DATASETS\n")
cat("===================================================================\n\n")

# --- Dragon 3: extract ALL available clinical variables ---
d3 <- data.frame(row.names = 1:nrow(dragon_raw))
d3$age <- as.numeric(dragon_raw$Currentage)
d3$sex <- ifelse(dragon_raw$Gender == "M", 1, ifelse(dragon_raw$Gender == "F", 0, NA))

# Lipids (fill from multiple timepoints)
d3$re_ldl <- as.numeric(dragon_raw$LDL_1)
d3$hdl    <- as.numeric(dragon_raw$HDL_1)
d3$tc     <- as.numeric(dragon_raw$TC_1)
d3$trig   <- as.numeric(dragon_raw$TRG_1)
for (suffix in c("_2", "_3", "_4")) {
  for (lip in list(c("LDL", "re_ldl"), c("HDL", "hdl"), c("TC", "tc"), c("TRG", "trig"))) {
    col <- paste0(lip[1], suffix)
    if (col %in% names(dragon_raw)) {
      fill <- is.na(d3[[lip[2]]]) & !is.na(suppressWarnings(as.numeric(dragon_raw[[col]])))
      d3[[lip[2]]][fill] <- as.numeric(dragon_raw[[col]][fill])
    }
  }
}
for (pair in list(c("LastLDL","re_ldl"), c("LastHDL","hdl"), c("LastTC","tc"), c("LastTrigs","trig"))) {
  if (pair[1] %in% names(dragon_raw)) {
    fill <- is.na(d3[[pair[2]]]) & !is.na(suppressWarnings(as.numeric(dragon_raw[[pair[1]]])))
    d3[[pair[2]]][fill] <- as.numeric(dragon_raw[[pair[1]]][fill])
  }
}

# ApoB (Dragon 3 unique strength)
d3$apob <- suppressWarnings(as.numeric(dragon_raw$ApoB))

# ApoA1
d3$apoa <- suppressWarnings(as.numeric(dragon_raw$ApoA1))

# Lpa (try multiple columns)
d3$lpa <- suppressWarnings(as.numeric(dragon_raw$Lpa))
if (sum(!is.na(d3$lpa)) < 20) {
  d3$lpa <- suppressWarnings(as.numeric(dragon_raw$Lpa_1))
}
if (sum(!is.na(d3$lpa)) < 20) {
  d3$lpa <- suppressWarnings(as.numeric(dragon_raw$Lpaunitsmgl))
}

# eGFR
d3$egfr <- suppressWarnings(as.numeric(dragon_raw$eGFR))

# ALT
d3$alt <- suppressWarnings(as.numeric(dragon_raw$ALT))

# Glucose
d3$glucose <- suppressWarnings(as.numeric(dragon_raw$Glucose_1))

# Smoking
smk <- dragon_raw$Smoking_binary
d3$smoking_binary <- ifelse(smk %in% c("1", "2"), 1, ifelse(smk %in% c("0"), 0, NA))

# Diabetes
dm <- dragon_raw$Diabetes_binary
d3$diabetes <- ifelse(dm %in% c(1, "1", "Yes", TRUE), 1,
                     ifelse(dm %in% c(0, "0", "No", FALSE, ""), 0, NA))

# Hypertension
bp_treat <- dragon_raw$onBPtreat
d3$hypertension <- ifelse(bp_treat %in% c("Y","y","Yes","1",1), 1,
                         ifelse(bp_treat %in% c("N","n","No","0",0), 0, NA))
# BP
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

# BMI
d3$bmi <- suppressWarnings(as.numeric(dragon_raw$BMI))
d3$bmi[!is.na(d3$bmi) & (d3$bmi < 10 | d3$bmi > 80)] <- NA

# Gene
d3$gene_apob <- as.integer(grepl("APOB", dragon_raw$Mutation1, ignore.case = TRUE))

# Pre-computed ratios from Dragon 3
d3$apob_ldl_ratio <- suppressWarnings(as.numeric(dragon_raw$APOB_LDL))

# Outcome
d3$ascvd_combined <- as.integer(dragon_raw$ASCVD_combined %in% c(1, "1", "1.0"))

# --- Derive additional features ---
d3$non_hdl <- d3$tc - d3$hdl
d3$tc_hdl_ratio <- ifelse(!is.na(d3$hdl) & d3$hdl > 0, d3$tc / d3$hdl, NA)
tg_hdl <- ifelse(!is.na(d3$hdl) & d3$hdl > 0 & !is.na(d3$trig) & d3$trig > 0,
                 d3$trig / d3$hdl, NA)
d3$log_tg_hdl <- ifelse(!is.na(tg_hdl) & tg_hdl > 0, log(tg_hdl), NA)
d3$pulse_pressure <- ifelse(!is.na(d3$sbp) & !is.na(d3$dbp), d3$sbp - d3$dbp, NA)
d3$log_lpa <- ifelse(!is.na(d3$lpa) & d3$lpa > 0, log(d3$lpa + 1), NA)
d3$log_apob_ldl <- ifelse(!is.na(d3$apob) & !is.na(d3$re_ldl) & d3$apob > 0 & d3$re_ldl > 0,
                          log(d3$apob / d3$re_ldl), NA)

cat(sprintf("  Dragon 3 harmonised: N=%d, Events=%d\n", nrow(d3), sum(d3$ascvd_combined)))

# Missingness summary for Dragon 3
cat("\n  Dragon 3 missingness:\n")
d3_analysis_vars <- c("age", "sex", "re_ldl", "hdl", "tc", "trig", "apob", "apoa",
                       "lpa", "egfr", "alt", "glucose", "smoking_binary", "diabetes",
                       "hypertension", "bmi", "sbp", "dbp", "gene_apob",
                       "non_hdl", "tc_hdl_ratio", "log_tg_hdl", "pulse_pressure",
                       "log_lpa", "log_apob_ldl", "apob_ldl_ratio")
for (v in d3_analysis_vars) {
  if (v %in% names(d3)) {
    n_miss <- sum(is.na(d3[[v]]))
    pct <- 100 * n_miss / nrow(d3)
    if (pct > 0) cat(sprintf("    %-18s: %d missing (%.1f%%)\n", v, n_miss, pct))
  }
}

# --- UKB harmonisation ---
# Safe column extraction: use NA vector if column missing
safe_col <- function(df, col, as_num = FALSE) {
  if (col %in% names(df)) {
    vals <- df[[col]]
    if (as_num) vals <- suppressWarnings(as.numeric(vals))
    return(vals)
  }
  rep(NA_real_, nrow(df))
}

ukb_h <- data.frame(
  age = ukb$age, sex = ukb$sex, re_ldl = ukb$re_ldl, hdl = ukb$hdl,
  tc = ukb$tc, trig = ukb$trig,
  apob = safe_col(ukb, "apob", TRUE),
  apoa = safe_col(ukb, "apoa", TRUE),
  lpa  = safe_col(ukb, "lpa", TRUE),
  egfr = safe_col(ukb, "egfr", TRUE),
  alt  = safe_col(ukb, "alt", TRUE),
  glucose = safe_col(ukb, "glucose", TRUE),
  smoking_binary = ukb$ever_smoked,
  diabetes = ukb$diabetes, hypertension = ukb$hypertension,
  bmi = ukb$bmi, sbp = ukb$sbp, dbp = ukb$dbp,
  gene_apob = as.integer(ukb$gene == "APOB"),
  ascvd_combined = ukb$ascvd_combined,
  stringsAsFactors = FALSE
)

# Report which UKB columns were found vs missing
ukb_biomarkers <- c("apob", "apoa", "lpa", "egfr", "alt", "glucose")
for (bm in ukb_biomarkers) {
  present <- bm %in% names(ukb) && sum(!is.na(ukb[[bm]])) > 0
  n_avail <- if (present) sum(!is.na(ukb[[bm]])) else 0
  cat(sprintf("  UKB %-10s: %s (n=%d)\n", bm,
              if (present) "FOUND" else "MISSING", n_avail))
}
# Derived
ukb_h$non_hdl <- ukb_h$tc - ukb_h$hdl
ukb_h$tc_hdl_ratio <- ifelse(!is.na(ukb_h$hdl) & ukb_h$hdl > 0, ukb_h$tc / ukb_h$hdl, NA)
tg_hdl_u <- ifelse(!is.na(ukb_h$hdl) & ukb_h$hdl > 0 & !is.na(ukb_h$trig) & ukb_h$trig > 0,
                   ukb_h$trig / ukb_h$hdl, NA)
ukb_h$log_tg_hdl <- ifelse(!is.na(tg_hdl_u) & tg_hdl_u > 0, log(tg_hdl_u), NA)
ukb_h$pulse_pressure <- ifelse(!is.na(ukb_h$sbp) & !is.na(ukb_h$dbp), ukb_h$sbp - ukb_h$dbp, NA)
ukb_h$log_lpa <- ifelse(!is.na(ukb_h$lpa) & ukb_h$lpa > 0, log(ukb_h$lpa + 1), NA)
ukb_h$log_apob_ldl <- ifelse(!is.na(ukb_h$apob) & !is.na(ukb_h$re_ldl) &
                               ukb_h$apob > 0 & ukb_h$re_ldl > 0,
                             log(ukb_h$apob / ukb_h$re_ldl), NA)
ukb_h$apob_ldl_ratio <- ifelse(!is.na(ukb_h$apob) & !is.na(ukb_h$re_ldl) & ukb_h$re_ldl > 0,
                                ukb_h$apob / ukb_h$re_ldl, NA)

cat(sprintf("  UKB harmonised: N=%d\n", nrow(ukb_h)))

# --- Wales harmonisation ---
wales_h <- data.frame(row.names = 1:nrow(wales_raw))
# Use fixed reference date for reproducibility (not Sys.Date())
WALES_REF_DATE <- as.Date("2025-01-01")
if ("DOB" %in% names(wales_raw)) {
  dob <- as.Date(wales_raw$DOB, format = "%Y-%m-%d")
  if (sum(!is.na(dob)) < 100) dob <- as.Date(wales_raw$DOB, format = "%d/%m/%Y")
  wales_h$age <- as.numeric(difftime(WALES_REF_DATE, dob, units = "days")) / 365.25
  cat(sprintf("  Wales age: computed from DOB with reference date %s\n", WALES_REF_DATE))
} else { wales_h$age <- NA }
wales_h$sex <- ifelse(wales_raw$Gender == "M", 1, ifelse(wales_raw$Gender == "F", 0, NA))
wales_h$re_ldl <- suppressWarnings(as.numeric(wales_raw$LDL.1))
wales_h$hdl    <- suppressWarnings(as.numeric(wales_raw$HDL.1))
wales_h$tc     <- suppressWarnings(as.numeric(wales_raw$TC.1))
wales_h$trig   <- suppressWarnings(as.numeric(wales_raw$TRG.1))
for (suffix in c(".2", ".3", ".4")) {
  for (lip in list(c("LDL", "re_ldl"), c("HDL", "hdl"), c("TC", "tc"), c("TRG", "trig"))) {
    col <- paste0(lip[1], suffix)
    if (col %in% names(wales_raw)) {
      fill <- is.na(wales_h[[lip[2]]]) & !is.na(suppressWarnings(as.numeric(wales_raw[[col]])))
      wales_h[[lip[2]]][fill] <- suppressWarnings(as.numeric(wales_raw[[col]][fill]))
    }
  }
}
smk <- wales_raw$Smoking
wales_h$smoking_binary <- ifelse(smk %in% c(1,"1","Yes","yes","Current","current","Ex","ex",
                                        "Former","former","Ex-smoker","Current smoker","Previous"), 1,
                            ifelse(smk %in% c(0,"0","No","no","Never","never","Non-smoker"), 0, NA))
dm <- wales_raw$Diabetes
wales_h$diabetes <- ifelse(dm %in% c(1,"1","Yes","yes",TRUE), 1,
                     ifelse(dm %in% c(0,"0","No","no",FALSE,""), 0, NA))
bp_med <- wales_raw$BloodPressureMedication
wales_h$hypertension <- ifelse(bp_med %in% c(1,"1","Yes","yes",TRUE), 1,
                         ifelse(bp_med %in% c(0,"0","No","no",FALSE,""), 0, NA))
sbp_w <- suppressWarnings(as.numeric(wales_raw$BloodPressureSystolic))
dbp_w <- suppressWarnings(as.numeric(wales_raw$BloodPressureDiastolic))
high_bp_w <- (!is.na(sbp_w) & sbp_w >= 140) | (!is.na(dbp_w) & dbp_w >= 90)
wales_h$hypertension[high_bp_w & (is.na(wales_h$hypertension) | wales_h$hypertension == 0)] <- 1
wales_h$sbp <- sbp_w; wales_h$dbp <- dbp_w
wales_h$bmi <- suppressWarnings(as.numeric(wales_raw$BMI))
wales_h$gene_apob <- as.integer(grepl("APOB", wales_raw$Mutation1, ignore.case = TRUE))
wales_h$ascvd_combined <- as.integer(wales_raw$ascvd_combine %in% c(1, "1", "1.0"))
# Derived
wales_h$non_hdl <- wales_h$tc - wales_h$hdl
wales_h$tc_hdl_ratio <- ifelse(!is.na(wales_h$hdl) & wales_h$hdl > 0, wales_h$tc / wales_h$hdl, NA)
tg_hdl_w <- ifelse(!is.na(wales_h$hdl) & wales_h$hdl > 0 & !is.na(wales_h$trig) & wales_h$trig > 0,
                   wales_h$trig / wales_h$hdl, NA)
wales_h$log_tg_hdl <- ifelse(!is.na(tg_hdl_w) & tg_hdl_w > 0, log(tg_hdl_w), NA)
wales_h$pulse_pressure <- ifelse(!is.na(wales_h$sbp) & !is.na(wales_h$dbp), wales_h$sbp - wales_h$dbp, NA)

cat(sprintf("  Wales harmonised: N=%d, Events=%d (%.1f%%)\n",
            nrow(wales_h), sum(wales_h$ascvd_combined), 100 * mean(wales_h$ascvd_combined)))

# =============================================================================
# SECTION 3: MCAR TEST ON DRAGON 3
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 3: MCAR TEST (Little's test on Dragon 3)\n")
cat("===================================================================\n\n")

# Select variables for MCAR test (only those with some missingness)
mcar_vars <- c("age", "sex", "re_ldl", "hdl", "tc", "trig", "apob", "lpa",
               "smoking_binary", "diabetes", "hypertension", "bmi",
               "sbp", "dbp", "egfr", "alt", "glucose")
mcar_vars <- mcar_vars[mcar_vars %in% names(d3)]

# Keep only vars with at least some missingness and some data
mcar_vars_use <- character(0)
for (v in mcar_vars) {
  n_miss <- sum(is.na(d3[[v]]))
  n_obs  <- sum(!is.na(d3[[v]]))
  if (n_miss > 0 && n_obs > 20) mcar_vars_use <- c(mcar_vars_use, v)
}

cat(sprintf("  Variables tested for MCAR: %d\n", length(mcar_vars_use)))
cat(sprintf("  Variables: %s\n", paste(mcar_vars_use, collapse = ", ")))

mcar_result <- tryCatch({
  naniar::mcar_test(d3[, mcar_vars_use])
}, error = function(e) {
  cat(sprintf("  Little's MCAR test failed: %s\n", e$message))
  cat("  Proceeding with MI regardless (conservative approach)\n")
  NULL
})

if (!is.null(mcar_result)) {
  cat(sprintf("\n  Little's MCAR test:\n"))
  cat(sprintf("    Chi-squared = %.2f\n", mcar_result$statistic))
  cat(sprintf("    df          = %d\n", mcar_result$df))
  cat(sprintf("    p-value     = %.6f\n", mcar_result$p.value))

  if (mcar_result$p.value < 0.05) {
    cat("    CONCLUSION: Data is NOT MCAR (p < 0.05). MI is appropriate.\n")
    cat("    Missing data mechanism is likely MAR or MNAR.\n")
  } else {
    cat("    CONCLUSION: Cannot reject MCAR (p >= 0.05).\n")
    cat("    MI is still appropriate and preferable to complete-case analysis.\n")
  }
}

# =============================================================================
# SECTION 4: MULTIPLE IMPUTATION (MICE) ON DRAGON 3
# =============================================================================

cat("\n===================================================================\n")
cat(sprintf("SECTION 4: MULTIPLE IMPUTATION (MICE, m=%d)\n", M_IMPUTATIONS))
cat("===================================================================\n\n")

# Define variables to impute — include outcome as predictor in imputation model
# but do NOT impute the outcome itself
imp_vars <- c("age", "sex", "re_ldl", "hdl", "tc", "trig", "apob", "lpa",
              "smoking_binary", "diabetes", "hypertension", "bmi",
              "sbp", "dbp", "egfr", "alt", "glucose", "gene_apob",
              "non_hdl", "tc_hdl_ratio", "log_tg_hdl",
              "ascvd_combined")
imp_vars <- imp_vars[imp_vars %in% names(d3)]

# Only keep vars with at least SOME data
imp_vars_final <- character(0)
for (v in imp_vars) {
  if (sum(!is.na(d3[[v]])) >= 20) imp_vars_final <- c(imp_vars_final, v)
}

d3_imp_data <- d3[, imp_vars_final]
cat(sprintf("  Imputing %d variables on %d Dragon 3 patients\n",
            length(imp_vars_final) - 1, nrow(d3_imp_data)))  # -1 for outcome

# Set up MICE: outcome should predict other variables but NOT be imputed
set.seed(2026)
ini <- mice(d3_imp_data, maxit = 0, printFlag = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

# Don't impute the outcome
meth["ascvd_combined"] <- ""

# Don't impute variables with 0 missing
for (v in names(meth)) {
  if (sum(is.na(d3_imp_data[[v]])) == 0) meth[v] <- ""
}

# Use PMM for continuous, logreg for binary
for (v in names(meth)) {
  if (meth[v] == "") next
  if (v %in% c("sex", "smoking_binary", "diabetes", "hypertension", "gene_apob")) {
    meth[v] <- "logreg"
  } else {
    meth[v] <- "pmm"
  }
}

cat("  Imputation methods:\n")
for (v in names(meth)) {
  if (meth[v] != "") {
    cat(sprintf("    %-18s: %s (%d missing)\n", v, meth[v], sum(is.na(d3_imp_data[[v]]))))
  }
}

cat(sprintf("\n  Running MICE with m=%d, maxit=%d...\n", M_IMPUTATIONS, MAX_ITER))

mice_fit <- tryCatch(
  mice(d3_imp_data, m = M_IMPUTATIONS, maxit = MAX_ITER,
       method = meth, predictorMatrix = pred,
       seed = 2026, printFlag = FALSE),
  error = function(e) {
    cat(sprintf("  MICE failed: %s\n", e$message))
    cat("  Falling back to single imputation with median/mode\n")
    NULL
  }
)

if (!is.null(mice_fit)) {
  cat("  MICE completed successfully.\n")
  cat(sprintf("  Convergence: check logged events = %d\n", nrow(mice_fit$loggedEvents)))
} else {
  cat("  WARNING: MICE failed. Using complete-case analysis.\n")
}

# =============================================================================
# SECTION 5: DEFINE FEATURE TIERS (EPV-aware)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 5: DEFINE FEATURE TIERS (EPV-aware)\n")
cat("===================================================================\n\n")

n_events_d3 <- sum(d3$ascvd_combined)
max_vars <- floor(n_events_d3 / EPV_MIN)
cat(sprintf("  Dragon 3 events: %d\n", n_events_d3))
cat(sprintf("  EPV minimum: %d\n", EPV_MIN))
cat(sprintf("  Maximum model variables: %d\n", max_vars))

# Tiers constrained by EPV
# Tier 1: SAFEHEART vars (6 vars, EPV = 10.3)
# Tier 2: Core + ApoB (8 vars, EPV = 7.75)
# Tier 3: Core + ApoB + Lpa (9 vars, EPV = 6.9)
# Tier 4: Extended (12 vars, EPV = 5.2) — at the limit

tiers <- list(
  "T1_safeheart" = c("age", "sex", "re_ldl", "hypertension", "bmi", "smoking_binary"),
  "T2_core" = c("age", "sex", "re_ldl", "hdl", "smoking_binary", "diabetes", "hypertension"),
  "T3_apob" = c("age", "sex", "re_ldl", "hdl", "apob", "smoking_binary", "diabetes", "hypertension"),
  "T4_apob_lpa" = c("age", "sex", "re_ldl", "hdl", "apob", "log_lpa",
                     "smoking_binary", "diabetes", "hypertension"),
  "T5_lipid_rich" = c("age", "sex", "re_ldl", "hdl", "trig", "apob", "log_apob_ldl",
                       "smoking_binary", "diabetes", "hypertension"),
  "T6_extended" = c("age", "sex", "re_ldl", "hdl", "trig", "apob", "log_lpa",
                    "non_hdl", "tc_hdl_ratio", "smoking_binary", "diabetes", "hypertension")
)

for (tn in names(tiers)) {
  nv <- length(tiers[[tn]])
  epv <- n_events_d3 / nv
  status <- if (epv >= EPV_MIN) "OK" else "RISKY"
  cat(sprintf("  %s: %d vars, EPV=%.1f [%s]\n", tn, nv, epv, status))
}

# =============================================================================
# SECTION 6: SAFEHEART-RE BASELINE ON ALL DATASETS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 6: SAFEHEART-RE BASELINE\n")
cat("===================================================================\n\n")

SAFEHEART_COEF <- list(
  intercept = -7.053,
  betas = c(age = 0.064, sex = 0.775, re_ldl = 0.109,
            hypertension = 0.431, bmi = 0.025, smoking_binary = 0.466)
)

apply_safeheart <- function(d, impute_bmi = 27.0) {
  lp <- rep(SAFEHEART_COEF$intercept, nrow(d))
  for (v in names(SAFEHEART_COEF$betas)) {
    vals <- d[[v]]
    if (v == "bmi") vals[is.na(vals)] <- impute_bmi
    lp <- lp + SAFEHEART_COEF$betas[v] * vals
  }
  1 / (1 + exp(-lp))
}

safe_auc <- function(y, pred) {
  valid <- !is.na(y) & !is.na(pred)
  if (sum(valid) < 20) return(list(auc = NA, lo = NA, hi = NA, n = sum(valid)))
  roc_obj <- pROC::roc(y[valid], pred[valid], quiet = TRUE)
  ci <- pROC::ci.auc(roc_obj)
  list(auc = as.numeric(pROC::auc(roc_obj)), lo = ci[1], hi = ci[3],
       n = sum(valid), events = sum(y[valid]), roc = roc_obj)
}

sh_d3    <- safe_auc(d3$ascvd_combined, apply_safeheart(d3))
sh_ukb   <- safe_auc(ukb_h$ascvd_combined, apply_safeheart(ukb_h))
sh_wales <- safe_auc(wales_h$ascvd_combined, apply_safeheart(wales_h))

cat(sprintf("  SAFEHEART [Dragon3-DEV]: AUC=%.4f (%.4f-%.4f), N=%d, Ev=%d\n",
            sh_d3$auc, sh_d3$lo, sh_d3$hi, sh_d3$n, sh_d3$events))
cat(sprintf("  SAFEHEART [UKB-EXT]:     AUC=%.4f (%.4f-%.4f), N=%d, Ev=%d\n",
            sh_ukb$auc, sh_ukb$lo, sh_ukb$hi, sh_ukb$n, sh_ukb$events))
cat(sprintf("  SAFEHEART [Wales-EXT]:   AUC=%.4f (%.4f-%.4f), N=%d, Ev=%d\n",
            sh_wales$auc, sh_wales$lo, sh_wales$hi, sh_wales$n, sh_wales$events))

# =============================================================================
# SECTION 7: MODEL DEVELOPMENT ON IMPUTED DATA (pooled across imputations)
# =============================================================================

cat("\n===================================================================\n")
cat(sprintf("SECTION 7: MODEL DEVELOPMENT (across %d imputations)\n", M_IMPUTATIONS))
cat("===================================================================\n\n")

# Helper: fit penalised logistic on one imputed dataset, return coefficients
fit_penlog <- function(imp_df, vars, alpha = 0) {
  # Recompute derived vars if needed
  if ("log_lpa" %in% vars && !"log_lpa" %in% names(imp_df)) {
    imp_df$log_lpa <- ifelse(!is.na(imp_df$lpa) & imp_df$lpa > 0, log(imp_df$lpa + 1), NA)
  }
  if ("log_apob_ldl" %in% vars && !"log_apob_ldl" %in% names(imp_df)) {
    imp_df$log_apob_ldl <- ifelse(!is.na(imp_df$apob) & !is.na(imp_df$re_ldl) &
                                    imp_df$apob > 0 & imp_df$re_ldl > 0,
                                  log(imp_df$apob / imp_df$re_ldl), NA)
  }
  if ("non_hdl" %in% vars && !"non_hdl" %in% names(imp_df)) {
    imp_df$non_hdl <- imp_df$tc - imp_df$hdl
  }
  if ("tc_hdl_ratio" %in% vars && !"tc_hdl_ratio" %in% names(imp_df)) {
    imp_df$tc_hdl_ratio <- ifelse(!is.na(imp_df$hdl) & imp_df$hdl > 0, imp_df$tc / imp_df$hdl, NA)
  }

  avail <- vars[vars %in% names(imp_df)]
  df_cc <- imp_df[complete.cases(imp_df[, c(avail, "ascvd_combined")]), ]
  if (nrow(df_cc) < 50) return(NULL)

  X <- as.matrix(df_cc[, avail])
  y <- df_cc$ascvd_combined

  set.seed(2026)
  cv_fit <- tryCatch(
    cv.glmnet(X, y, family = "binomial", alpha = alpha,
              nfolds = min(10, floor(nrow(df_cc) / 10)), type.measure = "auc"),
    error = function(e) NULL
  )
  if (is.null(cv_fit)) return(NULL)

  coefs <- as.numeric(coef(cv_fit, s = "lambda.min"))
  names(coefs) <- c("intercept", avail)
  cv_auc <- max(cv_fit$cvm, na.rm = TRUE)

  list(coefs = coefs, cv_auc = cv_auc, vars = avail, n = nrow(df_cc))
}

# Storage for all results
all_model_results <- data.frame()

for (tier_name in names(tiers)) {
  tier_vars <- tiers[[tier_name]]
  cat(sprintf("  --- %s (%d vars) ---\n", tier_name, length(tier_vars)))

  for (alpha in c(0, 0.1, 0.25, 0.5, 1.0)) {

    # Fit on each imputation and pool coefficients (Rubin's rules for point estimates)
    coef_matrix <- NULL
    cv_aucs <- numeric(0)

    if (!is.null(mice_fit)) {
      for (m_idx in 1:M_IMPUTATIONS) {
        imp_df <- complete(mice_fit, m_idx)
        result <- fit_penlog(imp_df, tier_vars, alpha)
        if (!is.null(result)) {
          if (is.null(coef_matrix)) {
            coef_matrix <- matrix(0, nrow = M_IMPUTATIONS, ncol = length(result$coefs))
            colnames(coef_matrix) <- names(result$coefs)
          }
          coef_matrix[m_idx, names(result$coefs)] <- result$coefs
          cv_aucs <- c(cv_aucs, result$cv_auc)
        }
      }
    } else {
      # Complete case fallback
      result <- fit_penlog(d3, tier_vars, alpha)
      if (!is.null(result)) {
        coef_matrix <- matrix(result$coefs, nrow = 1)
        colnames(coef_matrix) <- names(result$coefs)
        cv_aucs <- result$cv_auc
      }
    }

    if (is.null(coef_matrix) || length(cv_aucs) == 0) next

    # Pool: average coefficients across imputations (Rubin's rule — point estimate)
    pooled_coefs <- colMeans(coef_matrix, na.rm = TRUE)
    pooled_cv_auc <- mean(cv_aucs)

    # Apply pooled coefficients to external datasets
    apply_frozen <- function(d, coefs, vars) {
      avail <- vars[vars %in% names(d)]
      lp <- rep(coefs["intercept"], nrow(d))
      valid <- rep(TRUE, nrow(d))
      for (v in avail) {
        if (v %in% names(coefs)) {
          vals <- d[[v]]
          if (v == "bmi") vals[is.na(vals)] <- 27.0
          lp <- lp + coefs[v] * vals
          valid <- valid & !is.na(vals)
        }
      }
      valid <- valid & !is.na(d$ascvd_combined)
      pred <- 1 / (1 + exp(-lp))
      list(pred = pred, valid = valid)
    }

    vars_used <- setdiff(names(pooled_coefs), "intercept")

    # UKB external validation
    ukb_res <- apply_frozen(ukb_h, pooled_coefs, vars_used)
    if (sum(ukb_res$valid) < 30) next
    ukb_auc <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(ukb_h$ascvd_combined[ukb_res$valid],
                                       ukb_res$pred[ukb_res$valid], quiet = TRUE))),
      error = function(e) NA)

    # Wales external validation
    w_res <- apply_frozen(wales_h, pooled_coefs, vars_used)
    if (sum(w_res$valid) < 30) next
    w_auc <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(wales_h$ascvd_combined[w_res$valid],
                                       w_res$pred[w_res$valid], quiet = TRUE))),
      error = function(e) NA)

    if (is.na(ukb_auc) || is.na(w_auc)) next

    delta_ukb <- ukb_auc - sh_ukb$auc
    delta_w   <- w_auc - sh_wales$auc

    row <- data.frame(
      method = "PenLog_MI", tier = tier_name, alpha = alpha,
      n_vars = length(vars_used), d3_cv_auc = pooled_cv_auc,
      ukb_auc = ukb_auc, wales_auc = w_auc,
      delta_ukb = delta_ukb, delta_wales = delta_w,
      min_delta = min(delta_ukb, delta_w),
      stringsAsFactors = FALSE
    )
    all_model_results <- rbind(all_model_results, row)

    cat(sprintf("    alpha=%.2f: D3_CV=%.4f, UKB=%.4f(%+.4f), W=%.4f(%+.4f)\n",
                alpha, pooled_cv_auc, ukb_auc, delta_ukb, w_auc, delta_w))
  }
}

# =============================================================================
# SECTION 8: XGBOOST ON IMPUTED DATA (pooled predictions)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 8: XGBOOST ON IMPUTED DATA\n")
cat("===================================================================\n\n")

xgb_param_sets <- list(
  list(max_depth = 2, eta = 0.05, subsample = 0.8, colsample = 0.8, min_child = 10),
  list(max_depth = 3, eta = 0.03, subsample = 0.8, colsample = 0.8, min_child = 15),
  list(max_depth = 2, eta = 0.03, subsample = 0.7, colsample = 0.7, min_child = 20),
  list(max_depth = 3, eta = 0.05, subsample = 0.8, colsample = 0.6, min_child = 10),
  list(max_depth = 4, eta = 0.01, subsample = 0.8, colsample = 0.8, min_child = 20),
  list(max_depth = 2, eta = 0.1,  subsample = 0.8, colsample = 0.8, min_child = 15)
)

for (tier_name in names(tiers)) {
  tier_vars <- tiers[[tier_name]]
  cat(sprintf("  --- XGB %s ---\n", tier_name))

  for (pi in seq_along(xgb_param_sets)) {
    p <- xgb_param_sets[[pi]]

    # Pool predictions across imputations
    ukb_preds_all <- matrix(0, nrow = nrow(ukb_h), ncol = 0)
    w_preds_all   <- matrix(0, nrow = nrow(wales_h), ncol = 0)
    cv_aucs <- numeric(0)

    n_imp <- if (!is.null(mice_fit)) M_IMPUTATIONS else 1

    for (m_idx in 1:n_imp) {
      if (!is.null(mice_fit)) {
        imp_df <- complete(mice_fit, m_idx)
      } else {
        imp_df <- d3
      }

      # Recompute derived vars
      if ("log_lpa" %in% tier_vars && !"log_lpa" %in% names(imp_df)) {
        imp_df$log_lpa <- ifelse(!is.na(imp_df$lpa) & imp_df$lpa > 0, log(imp_df$lpa + 1), NA)
      }
      if ("log_apob_ldl" %in% tier_vars && !"log_apob_ldl" %in% names(imp_df)) {
        imp_df$log_apob_ldl <- ifelse(!is.na(imp_df$apob) & !is.na(imp_df$re_ldl) &
                                        imp_df$apob > 0 & imp_df$re_ldl > 0,
                                      log(imp_df$apob / imp_df$re_ldl), NA)
      }
      if ("non_hdl" %in% tier_vars && !"non_hdl" %in% names(imp_df)) {
        imp_df$non_hdl <- imp_df$tc - imp_df$hdl
      }
      if ("tc_hdl_ratio" %in% tier_vars && !"tc_hdl_ratio" %in% names(imp_df)) {
        imp_df$tc_hdl_ratio <- ifelse(!is.na(imp_df$hdl) & imp_df$hdl > 0, imp_df$tc / imp_df$hdl, NA)
      }

      avail <- tier_vars[tier_vars %in% names(imp_df)]
      df_cc <- imp_df[complete.cases(imp_df[, c(avail, "ascvd_combined")]), ]
      if (nrow(df_cc) < 50) next

      X_train <- as.matrix(df_cc[, avail])
      y_train <- df_cc$ascvd_combined

      dtrain <- xgb.DMatrix(data = X_train, label = y_train)

      params <- list(
        objective = "binary:logistic", eval_metric = "auc",
        max_depth = p$max_depth, eta = p$eta,
        subsample = p$subsample, colsample_bytree = p$colsample,
        min_child_weight = p$min_child, nthread = 1
      )

      set.seed(2026)
      cv_res <- tryCatch(
        xgb.cv(params = params, data = dtrain, nrounds = 300,
               nfold = 5, early_stopping_rounds = 20, verbose = 0, maximize = TRUE),
        error = function(e) NULL
      )
      if (is.null(cv_res)) next

      best_iter <- cv_res$best_iteration
      if (is.null(best_iter) || best_iter < 1) next
      cv_auc <- cv_res$evaluation_log$test_auc_mean[best_iter]
      cv_aucs <- c(cv_aucs, cv_auc)

      set.seed(2026)
      model <- xgb.train(params = params, data = dtrain, nrounds = best_iter, verbose = 0)

      # External predictions (need same vars)
      # UKB — XGBoost handles NAs natively, do NOT replace with 0
      ukb_avail <- avail[avail %in% names(ukb_h)]
      if (length(ukb_avail) == length(avail)) {
        ukb_X <- as.matrix(ukb_h[, avail])
        ukb_pred <- predict(model, xgb.DMatrix(data = ukb_X))
        ukb_preds_all <- cbind(ukb_preds_all, ukb_pred)
      }

      # Wales — XGBoost handles NAs natively
      w_avail <- avail[avail %in% names(wales_h)]
      if (length(w_avail) == length(avail)) {
        w_X <- as.matrix(wales_h[, avail])
        w_pred <- predict(model, xgb.DMatrix(data = w_X))
        w_preds_all <- cbind(w_preds_all, w_pred)
      }
    }

    if (ncol(ukb_preds_all) == 0 || ncol(w_preds_all) == 0) next

    # Pool predictions: average across imputations
    ukb_pred_pooled <- rowMeans(ukb_preds_all)
    w_pred_pooled   <- rowMeans(w_preds_all)

    ukb_auc <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(ukb_h$ascvd_combined, ukb_pred_pooled, quiet = TRUE))),
      error = function(e) NA)
    w_auc <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(wales_h$ascvd_combined, w_pred_pooled, quiet = TRUE))),
      error = function(e) NA)

    if (is.na(ukb_auc) || is.na(w_auc)) next

    delta_ukb <- ukb_auc - sh_ukb$auc
    delta_w   <- w_auc - sh_wales$auc

    row <- data.frame(
      method = "XGBoost_MI", tier = tier_name, alpha = NA,
      n_vars = length(tier_vars), d3_cv_auc = mean(cv_aucs),
      ukb_auc = ukb_auc, wales_auc = w_auc,
      delta_ukb = delta_ukb, delta_wales = delta_w,
      min_delta = min(delta_ukb, delta_w),
      stringsAsFactors = FALSE
    )
    all_model_results <- rbind(all_model_results, row)

    cat(sprintf("    p%d: D3_CV=%.4f, UKB=%.4f(%+.4f), W=%.4f(%+.4f)\n",
                pi, mean(cv_aucs), ukb_auc, delta_ukb, w_auc, delta_w))
  }
}

# =============================================================================
# SECTION 9: RANDOM FOREST ON IMPUTED DATA
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 9: RANDOM FOREST ON IMPUTED DATA\n")
cat("===================================================================\n\n")

rf_configs <- list(
  list(ntrees = 1000, mtry_frac = 0.33, min_node = 10),
  list(ntrees = 1000, mtry_frac = 0.50, min_node = 10),
  list(ntrees = 2000, mtry_frac = 0.33, min_node = 20),
  list(ntrees = 1000, mtry_frac = 0.67, min_node = 5),
  list(ntrees = 2000, mtry_frac = 0.50, min_node = 15)
)

for (tier_name in names(tiers)) {
  tier_vars <- tiers[[tier_name]]
  cat(sprintf("  --- RF %s ---\n", tier_name))

  for (ci in seq_along(rf_configs)) {
    cfg <- rf_configs[[ci]]
    mtry_val <- max(1, round(length(tier_vars) * cfg$mtry_frac))

    ukb_preds_all <- matrix(0, nrow = nrow(ukb_h), ncol = 0)
    w_preds_all   <- matrix(0, nrow = nrow(wales_h), ncol = 0)
    oob_aucs <- numeric(0)

    n_imp <- if (!is.null(mice_fit)) M_IMPUTATIONS else 1

    for (m_idx in 1:n_imp) {
      if (!is.null(mice_fit)) {
        imp_df <- complete(mice_fit, m_idx)
      } else {
        imp_df <- d3
      }

      # Recompute derived
      if ("log_lpa" %in% tier_vars && !"log_lpa" %in% names(imp_df))
        imp_df$log_lpa <- ifelse(!is.na(imp_df$lpa) & imp_df$lpa > 0, log(imp_df$lpa + 1), NA)
      if ("log_apob_ldl" %in% tier_vars && !"log_apob_ldl" %in% names(imp_df))
        imp_df$log_apob_ldl <- ifelse(!is.na(imp_df$apob) & !is.na(imp_df$re_ldl) &
                                        imp_df$apob > 0 & imp_df$re_ldl > 0,
                                      log(imp_df$apob / imp_df$re_ldl), NA)
      if ("non_hdl" %in% tier_vars && !"non_hdl" %in% names(imp_df))
        imp_df$non_hdl <- imp_df$tc - imp_df$hdl
      if ("tc_hdl_ratio" %in% tier_vars && !"tc_hdl_ratio" %in% names(imp_df))
        imp_df$tc_hdl_ratio <- ifelse(!is.na(imp_df$hdl) & imp_df$hdl > 0, imp_df$tc / imp_df$hdl, NA)

      avail <- tier_vars[tier_vars %in% names(imp_df)]
      df_cc <- imp_df[complete.cases(imp_df[, c(avail, "ascvd_combined")]), ]
      if (nrow(df_cc) < 50) next

      set.seed(2026)
      rf_model <- tryCatch(
        ranger(y = factor(df_cc$ascvd_combined), x = df_cc[, avail],
               num.trees = cfg$ntrees, mtry = mtry_val,
               min.node.size = cfg$min_node, probability = TRUE, seed = 2026),
        error = function(e) NULL
      )
      if (is.null(rf_model)) next

      oob_pred <- rf_model$predictions[, "1"]
      oob_auc <- tryCatch(
        as.numeric(pROC::auc(pROC::roc(df_cc$ascvd_combined, oob_pred, quiet = TRUE))),
        error = function(e) NA)
      if (!is.na(oob_auc)) oob_aucs <- c(oob_aucs, oob_auc)

      # External — ranger cannot handle NAs in prediction, use complete-case
      ukb_avail <- avail[avail %in% names(ukb_h)]
      if (length(ukb_avail) == length(avail)) {
        ukb_valid <- complete.cases(ukb_h[, avail])
        ukb_pred <- rep(NA_real_, nrow(ukb_h))
        if (sum(ukb_valid) > 10) {
          ukb_pred[ukb_valid] <- predict(rf_model, data = ukb_h[ukb_valid, avail])$predictions[, "1"]
        }
        ukb_preds_all <- cbind(ukb_preds_all, ukb_pred)
      }
      w_avail <- avail[avail %in% names(wales_h)]
      if (length(w_avail) == length(avail)) {
        w_valid <- complete.cases(wales_h[, avail])
        w_pred <- rep(NA_real_, nrow(wales_h))
        if (sum(w_valid) > 10) {
          w_pred[w_valid] <- predict(rf_model, data = wales_h[w_valid, avail])$predictions[, "1"]
        }
        w_preds_all <- cbind(w_preds_all, w_pred)
      }
    }

    if (ncol(ukb_preds_all) == 0 || ncol(w_preds_all) == 0) next

    # Pool predictions (na.rm=TRUE because RF uses complete-case per imputation)
    ukb_pred_pooled <- rowMeans(ukb_preds_all, na.rm = TRUE)
    w_pred_pooled   <- rowMeans(w_preds_all, na.rm = TRUE)

    # Evaluate only on rows that got predictions
    ukb_has_pred <- !is.nan(ukb_pred_pooled) & !is.na(ukb_pred_pooled) &
                    !is.na(ukb_h$ascvd_combined)
    w_has_pred   <- !is.nan(w_pred_pooled) & !is.na(w_pred_pooled) &
                    !is.na(wales_h$ascvd_combined)

    if (sum(ukb_has_pred) < 30 || sum(w_has_pred) < 30) next

    ukb_auc <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(ukb_h$ascvd_combined[ukb_has_pred],
                                       ukb_pred_pooled[ukb_has_pred], quiet = TRUE))),
      error = function(e) NA)
    w_auc <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(wales_h$ascvd_combined[w_has_pred],
                                       w_pred_pooled[w_has_pred], quiet = TRUE))),
      error = function(e) NA)

    if (is.na(ukb_auc) || is.na(w_auc)) next

    delta_ukb <- ukb_auc - sh_ukb$auc
    delta_w   <- w_auc - sh_wales$auc

    row <- data.frame(
      method = "RF_MI", tier = tier_name, alpha = NA,
      n_vars = length(tier_vars), d3_cv_auc = mean(oob_aucs, na.rm = TRUE),
      ukb_auc = ukb_auc, wales_auc = w_auc,
      delta_ukb = delta_ukb, delta_wales = delta_w,
      min_delta = min(delta_ukb, delta_w),
      stringsAsFactors = FALSE
    )
    all_model_results <- rbind(all_model_results, row)

    cat(sprintf("    c%d: OOB=%.4f, UKB=%.4f(%+.4f), W=%.4f(%+.4f)\n",
                ci, mean(oob_aucs, na.rm = TRUE), ukb_auc, delta_ukb, w_auc, delta_w))
  }
}

# =============================================================================
# SECTION 10: GRAND COMPARISON
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 10: GRAND COMPARISON\n")
cat("===================================================================\n\n")

# Add SAFEHEART baseline
safeheart_row <- data.frame(
  method = "SAFEHEART-RE", tier = "published", alpha = NA,
  n_vars = 6, d3_cv_auc = sh_d3$auc,
  ukb_auc = sh_ukb$auc, wales_auc = sh_wales$auc,
  delta_ukb = 0, delta_wales = 0, min_delta = 0,
  stringsAsFactors = FALSE
)
all_model_results <- rbind(safeheart_row, all_model_results)
all_model_results <- all_model_results[order(-all_model_results$min_delta), ]

# Print top 30
n_show <- min(30, nrow(all_model_results))
cat(sprintf("  %-14s %-14s %5s %5s %7s %7s %7s %8s %8s\n",
            "Method", "Tier", "Alpha", "Vars", "D3_CV", "UKB", "Wales", "dUKB", "dWales"))
cat(paste(rep("-", 100), collapse = ""), "\n")
for (i in 1:n_show) {
  r <- all_model_results[i, ]
  alpha_str <- if (is.na(r$alpha)) "  -" else sprintf("%.2f", r$alpha)
  cat(sprintf("  %-14s %-14s %5s %5d %7.4f %7.4f %7.4f %+8.4f %+8.4f\n",
              r$method, r$tier, alpha_str, r$n_vars, r$d3_cv_auc,
              r$ukb_auc, r$wales_auc, r$delta_ukb, r$delta_wales))
}

# =============================================================================
# SECTION 11: BEST MODEL — FROZEN COEFFICIENTS + DeLong
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 11: BEST MODEL DETAILS\n")
cat("===================================================================\n\n")

# Find best penalised logistic (has interpretable coefficients)
pen_rows <- all_model_results[all_model_results$method == "PenLog_MI", ]
if (nrow(pen_rows) > 0) {
  best_pen <- pen_rows[which.max(pen_rows$min_delta), ]
  cat(sprintf("  BEST PENALISED LOGISTIC:\n"))
  cat(sprintf("    Tier: %s, Alpha: %.2f\n", best_pen$tier, best_pen$alpha))
  cat(sprintf("    D3 CV-AUC: %.4f\n", best_pen$d3_cv_auc))
  cat(sprintf("    UKB AUC: %.4f (delta: %+.4f)\n", best_pen$ukb_auc, best_pen$delta_ukb))
  cat(sprintf("    Wales AUC: %.4f (delta: %+.4f)\n", best_pen$wales_auc, best_pen$delta_wales))

  # Refit to get pooled coefficients
  best_tier_vars <- tiers[[best_pen$tier]]
  coef_matrix <- NULL

  n_imp <- if (!is.null(mice_fit)) M_IMPUTATIONS else 1
  for (m_idx in 1:n_imp) {
    if (!is.null(mice_fit)) { imp_df <- complete(mice_fit, m_idx) } else { imp_df <- d3 }
    result <- fit_penlog(imp_df, best_tier_vars, best_pen$alpha)
    if (!is.null(result)) {
      if (is.null(coef_matrix)) {
        coef_matrix <- matrix(0, nrow = n_imp, ncol = length(result$coefs))
        colnames(coef_matrix) <- names(result$coefs)
      }
      coef_matrix[m_idx, names(result$coefs)] <- result$coefs
    }
  }

  if (!is.null(coef_matrix)) {
    pooled <- colMeans(coef_matrix, na.rm = TRUE)
    cat("\n  POOLED FROZEN COEFFICIENTS (across imputations):\n")
    for (v in names(pooled)) {
      if (v == "intercept") {
        cat(sprintf("    intercept = %.6f\n", pooled[v]))
      } else {
        or_val <- exp(pooled[v])
        cat(sprintf("    %-18s coef=% .6f  OR=%.4f\n", v, pooled[v], or_val))
      }
    }

    # Save coefficients
    coef_df <- data.frame(variable = names(pooled), coefficient = pooled,
                           OR = exp(pooled), stringsAsFactors = FALSE)
    write.csv(coef_df, paste0(TAB_DIR, "calon2_s13_best_coefs.csv"), row.names = FALSE)
    cat("\n  Saved: calon2_s13_best_coefs.csv\n")

    # --- DeLong test: best model vs SAFEHEART on each external dataset ---
    cat("\n  DeLong tests vs SAFEHEART-RE:\n")
    vars_used <- setdiff(names(pooled), "intercept")
    apply_pooled <- function(d) {
      avail <- vars_used[vars_used %in% names(d)]
      lp <- rep(pooled["intercept"], nrow(d))
      valid <- rep(TRUE, nrow(d))
      for (v in avail) {
        vals <- d[[v]]
        if (v == "bmi") vals[is.na(vals)] <- 27.0
        lp <- lp + pooled[v] * vals
        valid <- valid & !is.na(vals)
      }
      valid <- valid & !is.na(d$ascvd_combined)
      list(pred = 1 / (1 + exp(-lp)), valid = valid)
    }

    # UKB DeLong
    best_ukb <- apply_pooled(ukb_h)
    sh_ukb_pred <- apply_safeheart(ukb_h)
    ukb_both <- best_ukb$valid & !is.na(sh_ukb_pred) & !is.na(ukb_h$ascvd_combined)
    if (sum(ukb_both) > 30) {
      roc_best_ukb <- pROC::roc(ukb_h$ascvd_combined[ukb_both],
                                 best_ukb$pred[ukb_both], quiet = TRUE)
      roc_sh_ukb   <- pROC::roc(ukb_h$ascvd_combined[ukb_both],
                                 sh_ukb_pred[ukb_both], quiet = TRUE)
      dl_ukb <- tryCatch(pROC::roc.test(roc_best_ukb, roc_sh_ukb, method = "delong"),
                         error = function(e) NULL)
      if (!is.null(dl_ukb)) {
        cat(sprintf("    UKB: Best=%.4f vs SAFEHEART=%.4f, DeLong p=%.4f %s\n",
                    as.numeric(pROC::auc(roc_best_ukb)),
                    as.numeric(pROC::auc(roc_sh_ukb)),
                    dl_ukb$p.value,
                    if (dl_ukb$p.value < 0.05) "*SIGNIFICANT*" else "(ns)"))
      }
    }

    # Wales DeLong
    best_w <- apply_pooled(wales_h)
    sh_w_pred <- apply_safeheart(wales_h)
    w_both <- best_w$valid & !is.na(sh_w_pred) & !is.na(wales_h$ascvd_combined)
    if (sum(w_both) > 30) {
      roc_best_w <- pROC::roc(wales_h$ascvd_combined[w_both],
                                best_w$pred[w_both], quiet = TRUE)
      roc_sh_w   <- pROC::roc(wales_h$ascvd_combined[w_both],
                                sh_w_pred[w_both], quiet = TRUE)
      dl_w <- tryCatch(pROC::roc.test(roc_best_w, roc_sh_w, method = "delong"),
                       error = function(e) NULL)
      if (!is.null(dl_w)) {
        cat(sprintf("    Wales: Best=%.4f vs SAFEHEART=%.4f, DeLong p=%.4f %s\n",
                    as.numeric(pROC::auc(roc_best_w)),
                    as.numeric(pROC::auc(roc_sh_w)),
                    dl_w$p.value,
                    if (dl_w$p.value < 0.05) "*SIGNIFICANT*" else "(ns)"))
      }
    }
  }
}

# =============================================================================
# SECTION 12: SAVE ALL RESULTS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 12: SAVE RESULTS\n")
cat("===================================================================\n\n")

write.csv(all_model_results, paste0(TAB_DIR, "calon2_s13_all_results.csv"), row.names = FALSE)
cat(sprintf("  Saved: calon2_s13_all_results.csv (%d rows)\n", nrow(all_model_results)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n======================================================================\n")
cat("  SCRIPT 13 COMPLETE                                                  \n")
cat("======================================================================\n")
cat(sprintf("  Development: Dragon 3 (N=%d, Events=%d)\n", nrow(d3), sum(d3$ascvd_combined)))
cat(sprintf("  External Val: UKB (N=%d), Wales (N=%d)\n", nrow(ukb_h), nrow(wales_h)))
cat(sprintf("  MCAR test: %s\n",
            if (!is.null(mcar_result)) sprintf("p=%.4f", mcar_result$p.value) else "failed"))
cat(sprintf("  Imputation: MICE m=%d\n",
            if (!is.null(mice_fit)) M_IMPUTATIONS else 0))
cat(sprintf("  Total model configs: %d\n", nrow(all_model_results) - 1))  # -1 for SAFEHEART
cat(sprintf("  Best overall: %s / %s (min_delta=%+.4f)\n",
            all_model_results$method[1], all_model_results$tier[1],
            all_model_results$min_delta[1]))
cat(sprintf("  SAFEHEART baselines: UKB=%.4f, Wales=%.4f\n",
            sh_ukb$auc, sh_wales$auc))
cat("======================================================================\n")

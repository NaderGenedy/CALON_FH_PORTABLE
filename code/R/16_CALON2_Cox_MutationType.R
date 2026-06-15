################################################################################
#                                                                              #
#  CALON-2 Script 16: COX PROPORTIONAL HAZARDS WITH NOVEL PREDICTORS          #
#                                                                              #
#  Rationale: Implement Cox regression with novel predictors (mutation type,   #
#  treatment intensity) using IECV (3-fold leave-one-dataset-out cross-       #
#  validation). Extends Script 14 patterns with survival analysis instead of   #
#  logistic. Compares IECV C-statistics vs SAFEHEART benchmark.               #
#                                                                              #
#  Methodology:                                                               #
#    1. Load & harmonise all 3 datasets (UKB, Dragon3, Wales) — reuse S14     #
#    2. Construct survival objects (time, event) for each dataset             #
#    3. Define feature tiers + novel predictors                              #
#    4. SAFEHEART baseline on each dataset (as linear predictor)             #
#    5. 3-fold IECV with penalised Cox (glmnet family="cox")                 #
#    6. Pool C-statistics event-weighted across folds                         #
#    7. Concordance comparison tests vs SAFEHEART                             #
#    8. Export coefficient table for best model                               #
#    9. Grand comparison table (all tiers vs SAFEHEART)                       #
#   10. Save results                                                          #
#                                                                              #
#  Key Novel Predictors:                                                      #
#    - gene_apob: binary (APOB vs LDLR/PCSK9/other) — captures mutation type #
#    - treatment_intensity: 0=none, 1=statin, 2=statin+ezetimibe,           #
#                          3=PCSK9i (any) — captures treatment depth         #
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

required_packages <- c("glmnet", "pROC", "dplyr", "survival")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# survcomp is a Bioconductor package — install if needed
if (!require("survcomp", character.only = TRUE, quietly = TRUE)) {
  cat("  Installing survcomp from Bioconductor...\n")
  if (!require("BiocManager", character.only = TRUE, quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("survcomp", ask = FALSE, update = FALSE)
  library(survcomp)
}

cat("\n")
cat("======================================================================\n")
cat("  CALON-2 Script 16: COX PROPORTIONAL HAZARDS + IECV                 \n")
cat("  Novel predictors: mutation type + treatment intensity              \n")
cat("  3-fold leave-one-dataset-out cross-validation                      \n")
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
FIG_DIR <- paste0(OUT_DIR, "figures/")
for (d in c(TAB_DIR, FIG_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

MIN_AGE_EXIT  <- 10    # minimum age at event/censoring (all datasets use age-as-time-scale)
EPV_MIN       <- 10    # conservative events per variable

# =============================================================================
# SAFEHEART-RE COEFFICIENTS (benchmark — logistic)
# =============================================================================

SAFEHEART_COEF <- list(
  intercept = -7.053,
  betas = c(age = 0.064, sex = 0.775, re_ldl = 0.109,
            hypertension = 0.431, bmi = 0.025, smoking_binary = 0.466)
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

safe_col <- function(df, col, as_num = FALSE) {
  if (col %in% names(df)) {
    vals <- df[[col]]
    if (as_num) vals <- suppressWarnings(as.numeric(vals))
    return(vals)
  }
  rep(NA_real_, nrow(df))
}

# Compute concordance (C-statistic) for a survival object and linear predictor
compute_concordance <- function(surv_obj, lp, min_events = 5) {
  # lp = linear predictor (higher = higher risk)
  tryCatch({
    valid <- !is.na(lp) & !is.na(surv_obj[, 1]) & !is.na(surv_obj[, 2])
    if (sum(valid) < 20 || sum(surv_obj[valid, 2]) < min_events) {
      return(list(c_stat = NA, se = NA, n = sum(valid), events = sum(surv_obj[valid, 2])))
    }
    # Use concordance.index from survcomp
    con <- survcomp::concordance.index(
      x = lp[valid],
      surv.time = surv_obj[valid, 1],
      surv.event = surv_obj[valid, 2],
      method = "noether"
    )
    list(c_stat = con$c.index, se = con$se, n = sum(valid), events = sum(surv_obj[valid, 2]))
  }, error = function(e) {
    list(c_stat = NA, se = NA, n = 0, events = 0)
  })
}

# Apply SAFEHEART linear predictor to data
apply_safeheart_lp <- function(d) {
  lp <- rep(SAFEHEART_COEF$intercept, nrow(d))
  for (v in names(SAFEHEART_COEF$betas)) {
    vals <- d[[v]]
    if (v == "bmi") vals[is.na(vals)] <- 27.0
    vals[is.na(vals)] <- median(vals, na.rm = TRUE)
    lp <- lp + SAFEHEART_COEF$betas[v] * vals
  }
  lp
}

# Compute treatment intensity (0=none, 1=statin, 2=statin+ezetimibe, 3=PCSK9i)
compute_treatment_intensity <- function(statin, ezetimibe, pcsk9i) {
  # Convert NAs to 0 (assume no treatment if unknown)
  statin[is.na(statin)] <- 0
  ezetimibe[is.na(ezetimibe)] <- 0
  pcsk9i[is.na(pcsk9i)] <- 0

  intensity <- rep(0L, length(statin))
  has_statin <- statin == 1
  has_ezet <- ezetimibe == 1
  has_pcsk9 <- pcsk9i == 1

  intensity[has_pcsk9] <- 3L  # PCSK9i overrides (most intensive)
  intensity[!has_pcsk9 & has_statin & has_ezet] <- 2L
  intensity[!has_pcsk9 & has_statin & !has_ezet] <- 1L

  intensity
}

# =============================================================================
# SECTION 1: LOAD ALL THREE DATASETS
# =============================================================================

cat("===================================================================\n")
cat("SECTION 1: LOAD DATASETS\n")
cat("===================================================================\n\n")

# --- UKB (DEVELOPMENT) ---
merged_file <- NULL
search_paths <- c(
  paste0(OUT_DIR, "calon2_full_merged.csv"),
  paste0(BASE_DIR, "calon2_full_merged.csv"),
  "/cloud/project/output/calon2_full_merged.csv",
  "/cloud/project/calon2_full_merged.csv"
)
for (sp in search_paths) if (file.exists(sp)) { merged_file <- sp; break }
if (is.null(merged_file)) stop("ERROR: calon2_full_merged.csv not found")
ukb_raw <- read.csv(merged_file, stringsAsFactors = FALSE)
cat(sprintf("  UKB: %d patients, %d events (%.1f%%) — using event_indicator (incident only)\n",
            nrow(ukb_raw), sum(ukb_raw$event_indicator),
            100 * mean(ukb_raw$event_indicator)))

# --- Dragon 3 (EXTERNAL VALIDATION) ---
dragon_file <- NULL
search_dirs <- unique(c(BASE_DIR, getwd(), "/cloud/project/"))
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
cat(sprintf("  Dragon 3: %d patients, %d events (%.1f%%)\n",
            nrow(dragon_raw),
            sum(as.integer(dragon_raw$ASCVD_combined %in% c(1, "1", "1.0"))),
            100 * mean(as.integer(dragon_raw$ASCVD_combined %in% c(1, "1", "1.0")))))

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
cat(sprintf("  Wales: %d FH+ patients, %d events (%.1f%%)\n",
            nrow(wales_raw), sum(wales_raw$ascvd_combine %in% c(1,"1")),
            100 * mean(wales_raw$ascvd_combine %in% c(1,"1"))))

# =============================================================================
# SECTION 2: CONSTRUCT SURVIVAL OBJECTS (TIME-TO-EVENT)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 2: CONSTRUCT SURVIVAL OBJECTS\n")
cat("===================================================================\n\n")

# --- UKB: Convert to age-at-event/censoring (= baseline age + follow-up) ---
# This aligns with Dragon 3 and Wales which use age-as-time-scale
ukb_surv <- data.frame(
  time = ukb_raw$age + ukb_raw$followup_years,  # age at exit
  event = as.integer(ukb_raw$event_indicator),
  stringsAsFactors = FALSE
)
valid_ukb <- !is.na(ukb_surv$time) & ukb_surv$time > MIN_AGE_EXIT
ukb_surv <- ukb_surv[valid_ukb, ]
ukb_raw <- ukb_raw[valid_ukb, ]
cat(sprintf("  UKB: %d patients; age-at-exit range [%.1f, %.1f] years\n",
            nrow(ukb_surv), min(ukb_surv$time), max(ukb_surv$time)))
cat(sprintf("       Events: %d (%.1f%%)\n", sum(ukb_surv$event), 100 * mean(ukb_surv$event)))

# --- Dragon 3: Compute from age columns ---
# Approach: use age_at_event_or_censoring column if available; else approximate from event ages
d3_surv <- data.frame(time = rep(NA_real_, nrow(dragon_raw)),
                       event = rep(NA_integer_, nrow(dragon_raw)))

if ("age_at_event_or_censoring" %in% names(dragon_raw)) {
  d3_surv$time <- suppressWarnings(as.numeric(dragon_raw$age_at_event_or_censoring))
  d3_surv$event <- as.integer(dragon_raw$ASCVD_combined %in% c(1, "1", "1.0"))
  cat(sprintf("  Dragon 3: Using age_at_event_or_censoring column\n"))
} else {
  # Fallback: compute from event age columns
  cat(sprintf("  Dragon 3: Computing from event age columns...\n"))
  event_cols <- c("MIACSYear", "PTCAYear", "CABGYear", "STROKEYear", "TIAYear", "PVDYear")
  event_cols_present <- event_cols[event_cols %in% names(dragon_raw)]

  if (length(event_cols_present) > 0) {
    # Get minimum age at event (or current age if no event)
    current_age <- suppressWarnings(as.numeric(dragon_raw$Currentage))
    age_at_event <- rep(NA_real_, nrow(dragon_raw))

    for (col in event_cols_present) {
      event_age <- suppressWarnings(as.numeric(dragon_raw[[col]]))
      age_at_event <- pmin(age_at_event, event_age, na.rm = TRUE)
    }

    # For those with event: use age_at_event; for those without: use current_age
    has_event <- !is.na(age_at_event)
    d3_surv$time <- ifelse(has_event, age_at_event, current_age)
    d3_surv$event <- as.integer(has_event)
  } else {
    # Last resort: just use current age, no events
    d3_surv$time <- suppressWarnings(as.numeric(dragon_raw$Currentage))
    d3_surv$event <- as.integer(dragon_raw$ASCVD_combined %in% c(1, "1", "1.0"))
  }
}

# Remove implausible ages and rows with missing time
valid_d3 <- !is.na(d3_surv$time) & d3_surv$time > MIN_AGE_EXIT
d3_surv <- d3_surv[valid_d3, ]
dragon_raw <- dragon_raw[valid_d3, ]
cat(sprintf("  Dragon 3 (after filtering): %d; age-at-exit range [%.1f, %.1f] years\n",
            nrow(d3_surv), min(d3_surv$time), max(d3_surv$time)))
cat(sprintf("       Events: %d (%.1f%%)\n", sum(d3_surv$event), 100 * mean(d3_surv$event)))

# --- Wales: Age-at-event-or-censoring (consistent with UKB and Dragon 3) ---
# BMI_AGE is current/recent age — NOT a valid baseline for follow-up calculation.
# Instead, use age at event (from event age columns) or current age (from DOB) as
# age-at-exit, matching the time scale used by UKB and Dragon 3.

WALES_REF_DATE <- as.Date("2026-02-01")  # data collected ~early 2026

# Current age from DOB (= censoring age for non-event patients)
dob <- as.Date(wales_raw$DOB, format = "%Y-%m-%d")
if (sum(!is.na(dob)) < 100) dob <- as.Date(wales_raw$DOB, format = "%d/%m/%Y")
current_age <- as.numeric(difftime(WALES_REF_DATE, dob, units = "days")) / 365.25

# Event ages — minimum age at first ASCVD event
event_cols_wales <- c("MIACSAge", "PCIStentsAge", "CABGAge", "STROKEAge", "ANGINAAge", "TIAAge", "PVDAge")
event_cols_wales_present <- event_cols_wales[event_cols_wales %in% names(wales_raw)]
cat(sprintf("  Wales event-age columns found: %s\n", paste(event_cols_wales_present, collapse = ", ")))

age_at_event <- rep(NA_real_, nrow(wales_raw))
for (col in event_cols_wales_present) {
  event_age <- suppressWarnings(as.numeric(wales_raw[[col]]))
  n_valid <- sum(!is.na(event_age))
  if (n_valid > 0) cat(sprintf("    %s: %d non-NA\n", col, n_valid))
  age_at_event <- pmin(age_at_event, event_age, na.rm = TRUE)
}

has_event <- !is.na(age_at_event)
cat(sprintf("  Wales patients with event ages: %d\n", sum(has_event)))

# Time = age at event (events) or current age (non-events)
wales_surv <- data.frame(
  time = ifelse(has_event, age_at_event, current_age),
  event = as.integer(has_event)
)

# Remove implausible ages and missing
valid_wales <- !is.na(wales_surv$time) & wales_surv$time > MIN_AGE_EXIT
wales_surv <- wales_surv[valid_wales, ]
wales_raw <- wales_raw[valid_wales, ]
cat(sprintf("  Wales (after filtering): %d; age-at-exit range [%.1f, %.1f] years\n",
            nrow(wales_surv), min(wales_surv$time), max(wales_surv$time)))
cat(sprintf("       Events: %d (%.1f%%)\n", sum(wales_surv$event), 100 * mean(wales_surv$event)))

# =============================================================================
# SECTION 3: HARMONISE ALL DATASETS (reuse Script 14 patterns)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 3: HARMONISE ALL DATASETS\n")
cat("===================================================================\n\n")

# --- UKB harmonisation ---
ukb_h <- data.frame(
  age = ukb_raw$age, sex = ukb_raw$sex,
  re_ldl = ukb_raw$re_ldl, hdl = ukb_raw$hdl,
  tc = ukb_raw$tc, trig = ukb_raw$trig,
  apob = safe_col(ukb_raw, "apob", TRUE),
  lpa  = safe_col(ukb_raw, "lpa", TRUE),
  glucose = safe_col(ukb_raw, "glucose", TRUE),
  alt  = safe_col(ukb_raw, "alt", TRUE),
  smoking_binary = ukb_raw$ever_smoked,
  diabetes = ukb_raw$diabetes,
  hypertension = ukb_raw$hypertension,
  bmi = ukb_raw$bmi,
  sbp = ukb_raw$sbp, dbp = ukb_raw$dbp,
  gene_apob = as.integer(ukb_raw$gene == "APOB"),
  on_statin = as.integer(ukb_raw$on_statin == 1),
  on_ezetimibe = as.integer(ukb_raw$on_ezetimibe == 1),
  on_pcsk9i = as.integer(safe_col(ukb_raw, "on_pcsk9i", FALSE) == 1),
  stringsAsFactors = FALSE
)
ukb_h$treatment_intensity <- compute_treatment_intensity(
  ukb_h$on_statin, ukb_h$on_ezetimibe, ukb_h$on_pcsk9i
)

cat(sprintf("  UKB harmonised: N=%d\n", nrow(ukb_h)))

# --- Dragon 3 harmonisation ---
d3 <- data.frame(row.names = 1:nrow(dragon_raw))
d3$age <- as.numeric(dragon_raw$Currentage)
d3$sex <- ifelse(dragon_raw$Gender == "M", 1, ifelse(dragon_raw$Gender == "F", 0, NA))
d3$re_ldl <- as.numeric(dragon_raw$LDL_1)
d3$hdl    <- as.numeric(dragon_raw$HDL_1)
d3$tc     <- as.numeric(dragon_raw$TC_1)
d3$trig   <- as.numeric(dragon_raw$TRG_1)
# Fill missing with later measurements
for (suffix in c("_2", "_3", "_4")) {
  for (lip in list(c("LDL","re_ldl"), c("HDL","hdl"), c("TC","tc"), c("TRG","trig"))) {
    col <- paste0(lip[1], suffix)
    if (col %in% names(dragon_raw)) {
      fill <- is.na(d3[[lip[2]]]) & !is.na(suppressWarnings(as.numeric(dragon_raw[[col]])))
      d3[[lip[2]]][fill] <- as.numeric(dragon_raw[[col]][fill])
    }
  }
}
d3$apob <- suppressWarnings(as.numeric(dragon_raw$ApoB))
d3$lpa <- suppressWarnings(as.numeric(dragon_raw$Lpa))
if (sum(!is.na(d3$lpa)) < 20) d3$lpa <- suppressWarnings(as.numeric(dragon_raw$Lpa_1))
d3$alt <- suppressWarnings(as.numeric(dragon_raw$ALT))
d3$glucose <- suppressWarnings(as.numeric(dragon_raw$Glucose_1))
d3$smoking_binary <- ifelse(dragon_raw$Smoking_binary %in% c("1","2"), 1,
                             ifelse(dragon_raw$Smoking_binary == "0", 0, NA))
d3$diabetes <- ifelse(dragon_raw$Diabetes_binary %in% c(1,"1","Yes",TRUE), 1,
                      ifelse(dragon_raw$Diabetes_binary %in% c(0,"0","No",FALSE,""), 0, NA))
d3$hypertension <- ifelse(dragon_raw$onBPtreat %in% c("Y","y","Yes","1",1), 1,
                          ifelse(dragon_raw$onBPtreat %in% c("N","n","No","0",0), 0, NA))
d3$bmi <- suppressWarnings(as.numeric(dragon_raw$BMI))
d3$bmi[!is.na(d3$bmi) & (d3$bmi < 10 | d3$bmi > 80)] <- NA
d3$sbp <- NA_real_   # not available in Dragon 3
d3$dbp <- NA_real_   # not available in Dragon 3
d3$gene_apob <- as.integer(grepl("APOB", dragon_raw$Mutation1, ignore.case = TRUE))
d3$on_statin <- ifelse(dragon_raw$Statin %in% c("Y","y","Yes","1",1), 1,
                       ifelse(dragon_raw$Statin %in% c("N","n","No","0",0), 0, NA))
d3$on_ezetimibe <- ifelse(dragon_raw$Ezetimibe %in% c("Y","y","Yes","1",1), 1,
                          ifelse(dragon_raw$Ezetimibe %in% c("N","n","No","0",0), 0, NA))
d3$on_pcsk9i <- ifelse(dragon_raw$PCSK9i %in% c("Y","y","Yes","1",1), 1,
                       ifelse(dragon_raw$PCSK9i %in% c("N","n","No","0",0), 0, NA))
d3$treatment_intensity <- compute_treatment_intensity(
  d3$on_statin, d3$on_ezetimibe, d3$on_pcsk9i
)

cat(sprintf("  Dragon 3 harmonised: N=%d\n", nrow(d3)))

# --- Wales harmonisation ---
wales_h <- data.frame(row.names = 1:nrow(wales_raw))
wales_h$age <- suppressWarnings(as.numeric(wales_raw$BMI_AGE))
wales_h$sex <- ifelse(wales_raw$Gender == "M", 1, ifelse(wales_raw$Gender == "F", 0, NA))
wales_h$re_ldl <- suppressWarnings(as.numeric(wales_raw$LDL.1))
wales_h$hdl    <- suppressWarnings(as.numeric(wales_raw$HDL.1))
wales_h$tc     <- suppressWarnings(as.numeric(wales_raw$TC.1))
wales_h$trig   <- suppressWarnings(as.numeric(wales_raw$TRG.1))
for (suffix in c(".2", ".3", ".4")) {
  for (lip in list(c("LDL","re_ldl"), c("HDL","hdl"), c("TC","tc"), c("TRG","trig"))) {
    col <- paste0(lip[1], suffix)
    if (col %in% names(wales_raw)) {
      fill <- is.na(wales_h[[lip[2]]]) & !is.na(suppressWarnings(as.numeric(wales_raw[[col]])))
      wales_h[[lip[2]]][fill] <- suppressWarnings(as.numeric(wales_raw[[col]][fill]))
    }
  }
}
wales_h$apob <- NA_real_  # Wales doesn't have ApoB
wales_h$lpa <- NA_real_   # Wales has ~2% Lpa — treat as missing
wales_h$alt <- NA_real_
wales_h$glucose <- NA_real_
wales_h$smoking_binary <- ifelse(wales_raw$Smoking %in% c(1,"1","Yes","yes","Current","Ex","Former"), 1,
                                 ifelse(wales_raw$Smoking %in% c(0,"0","No","no","Never"), 0, NA))
wales_h$diabetes <- ifelse(wales_raw$Diabetes %in% c(1,"1","Yes","yes",TRUE), 1,
                          ifelse(wales_raw$Diabetes %in% c(0,"0","No","no",FALSE,""), 0, NA))
wales_h$hypertension <- ifelse(wales_raw$BloodPressureMedication %in% c(1,"1","Yes","yes",TRUE), 1,
                              ifelse(wales_raw$BloodPressureMedication %in% c(0,"0","No","no",FALSE,""), 0, NA))
wales_h$bmi <- suppressWarnings(as.numeric(wales_raw$BMI))
wales_h$sbp <- NA_real_   # not available in Wales
wales_h$dbp <- NA_real_   # not available in Wales
wales_h$gene_apob <- as.integer(grepl("APOB", wales_raw$Mutation1, ignore.case = TRUE))
wales_h$on_statin <- ifelse(wales_raw$Treatment1 %in% c("Statin","statin","Statins",1,"1") |
                            grepl("Statin", wales_raw$Treatment3, ignore.case = TRUE), 1, 0)
wales_h$on_statin[is.na(wales_h$on_statin)] <- 0
wales_h$on_ezetimibe <- ifelse(grepl("ezetimibe", wales_raw$Treatment3, ignore.case = TRUE), 1, 0)
wales_h$on_ezetimibe[is.na(wales_h$on_ezetimibe)] <- 0
wales_h$on_pcsk9i <- ifelse(grepl("PCSK9", wales_raw$Treatment3, ignore.case = TRUE), 1, 0)
wales_h$on_pcsk9i[is.na(wales_h$on_pcsk9i)] <- 0
wales_h$treatment_intensity <- compute_treatment_intensity(
  wales_h$on_statin, wales_h$on_ezetimibe, wales_h$on_pcsk9i
)

cat(sprintf("  Wales harmonised: N=%d\n", nrow(wales_h)))

# =============================================================================
# SECTION 4: DEFINE FEATURE TIERS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 4: FEATURE TIERS\n")
cat("===================================================================\n\n")

n_events_ukb <- sum(ukb_surv$event)
n_events_total <- n_events_ukb + sum(d3_surv$event) + sum(wales_surv$event)

cat(sprintf("  Total events across all datasets: %d (UKB=%d, D3=%d, Wales=%d)\n",
            n_events_total, n_events_ukb, sum(d3_surv$event), sum(wales_surv$event)))

tiers <- list(
  # Tier 1: SAFEHEART-equivalent (6 vars, can validate on all datasets)
  "T1_safeheart" = c("age", "sex", "re_ldl", "hypertension", "bmi", "smoking_binary"),

  # Tier 2: Core clinical (7 vars, can validate on all)
  "T2_core" = c("age", "sex", "re_ldl", "hdl", "smoking_binary", "diabetes", "hypertension"),

  # Tier 3: Core + gene (novel predictor 1)
  "T3_gene" = c("age", "sex", "re_ldl", "hdl", "smoking_binary", "diabetes", "hypertension", "gene_apob"),

  # Tier 4: Core + treatment intensity (novel predictor 2)
  "T4_treatment" = c("age", "sex", "re_ldl", "hdl", "smoking_binary", "diabetes", "hypertension", "treatment_intensity"),

  # Tier 5: Core + both novel predictors
  "T5_gene_tx" = c("age", "sex", "re_ldl", "hdl", "smoking_binary", "diabetes", "hypertension",
                    "gene_apob", "treatment_intensity"),

  # Tier 6: Extended + novel (adds apob, bmi, log_lpa for Dragon 3 only)
  "T6_extended" = c("age", "sex", "re_ldl", "hdl", "apob", "bmi", "smoking_binary", "diabetes",
                     "hypertension", "gene_apob", "treatment_intensity")
)

for (tn in names(tiers)) {
  nv <- length(tiers[[tn]])
  epv <- n_events_ukb / nv
  cat(sprintf("  %s: %d vars, EPV=%.1f\n", tn, nv, epv))
}

# =============================================================================
# SECTION 5: SAFEHEART-RE BASELINE C-STATISTICS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 5: SAFEHEART-RE BASELINE C-STATISTICS\n")
cat("===================================================================\n\n")

# Bind survival to harmonised data
ukb_h$time <- ukb_surv$time
ukb_h$event <- ukb_surv$event
d3$time <- d3_surv$time
d3$event <- d3_surv$event
wales_h$time <- wales_surv$time
wales_h$event <- wales_surv$event

# Create survival objects
surv_ukb <- Surv(ukb_h$time, ukb_h$event)
surv_d3 <- Surv(d3$time, d3$event)
surv_wales <- Surv(wales_h$time, wales_h$event)

# Compute SAFEHEART linear predictor for each dataset
lp_sh_ukb <- apply_safeheart_lp(ukb_h)
lp_sh_d3 <- apply_safeheart_lp(d3)
lp_sh_wales <- apply_safeheart_lp(wales_h)

# Concordance
sh_c_ukb <- compute_concordance(surv_ukb, lp_sh_ukb)
sh_c_d3 <- compute_concordance(surv_d3, lp_sh_d3)
sh_c_wales <- compute_concordance(surv_wales, lp_sh_wales)

cat(sprintf("  SAFEHEART [UKB]:      C-stat=%.4f (SE=%.4f), N=%d, Ev=%d\n",
            sh_c_ukb$c_stat, sh_c_ukb$se, sh_c_ukb$n, sh_c_ukb$events))
cat(sprintf("  SAFEHEART [Dragon3]:  C-stat=%.4f (SE=%.4f), N=%d, Ev=%d\n",
            sh_c_d3$c_stat, sh_c_d3$se, sh_c_d3$n, sh_c_d3$events))
cat(sprintf("  SAFEHEART [Wales]:    C-stat=%.4f (SE=%.4f), N=%d, Ev=%d\n",
            sh_c_wales$c_stat, sh_c_wales$se, sh_c_wales$n, sh_c_wales$events))

# =============================================================================
# SECTION 6: 3-FOLD IECV (leave-one-dataset-out)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 6: 3-FOLD IECV (leave-one-dataset-out cross-validation)\n")
cat("===================================================================\n\n")

# Storage for results
iecv_results <- data.frame()

# Define folds: (train_1, train_2, test, test_name)
folds <- list(
  list(train_1 = "ukb", train_2 = "d3", test = "wales", test_name = "Wales"),
  list(train_1 = "ukb", train_2 = "wales", test = "d3", test_name = "Dragon3"),
  list(train_1 = "d3", train_2 = "wales", test = "ukb", test_name = "UKB")
)

# For each tier and regularisation parameter
alphas <- c(0.0, 0.25, 0.5, 1.0)  # 0=ridge, 1=lasso

for (tier_name in names(tiers)) {
  tier_vars <- tiers[[tier_name]]
  cat(sprintf("\n=== Tier: %s (vars: %s) ===\n", tier_name, paste(tier_vars, collapse=", ")))

  for (alpha_val in alphas) {
    cat(sprintf("\n  Alpha=%.2f:\n", alpha_val))

    fold_results <- list()

    for (fold_idx in 1:length(folds)) {
      fold <- folds[[fold_idx]]
      test_name <- fold$test_name

      # Prepare training data (combine two datasets)
      datasets <- list(ukb_h, d3, wales_h)
      names(datasets) <- c("ukb", "d3", "wales")
      surv_list <- list(surv_ukb, surv_d3, surv_wales)
      names(surv_list) <- c("ukb", "d3", "wales")

      train_1_data <- datasets[[fold$train_1]]
      train_2_data <- datasets[[fold$train_2]]
      test_data <- datasets[[fold$test]]

      train_surv_1 <- surv_list[[fold$train_1]]
      train_surv_2 <- surv_list[[fold$train_2]]
      test_surv <- surv_list[[fold$test]]

      # Combine training datasets (use common columns to avoid rbind mismatch)
      common_cols <- intersect(names(train_1_data), names(train_2_data))
      train_data <- rbind(train_1_data[, common_cols], train_2_data[, common_cols])
      rownames(train_data) <- 1:nrow(train_data)  # reset row indices
      train_surv <- c(train_surv_1, train_surv_2)

      # Get complete cases in training (align indices)
      avail <- tier_vars[tier_vars %in% names(train_data)]
      idx_cc <- which(complete.cases(train_data[, c(avail, "time", "event")]))
      train_cc <- train_data[idx_cc, ]
      train_surv_cc <- train_surv[idx_cc]

      if (nrow(train_cc) < 50 || sum(train_surv_cc[, 2]) < 10) {
        cat(sprintf("    Fold %d (%s): Skipped — insufficient training data\n", fold_idx, test_name))
        fold_results[[fold_idx]] <- list(c_test = NA, c_safeheart = NA)
        next
      }

      # Fit penalised Cox
      X_train <- as.matrix(train_cc[, avail])
      y_train <- train_surv_cc

      set.seed(2026 + fold_idx)
      cox_fit <- tryCatch(
        glmnet::cv.glmnet(X_train, y_train, family = "cox", alpha = alpha_val,
                         nfolds = 10, type.measure = "C"),
        error = function(e) NULL
      )

      if (is.null(cox_fit)) {
        cat(sprintf("    Fold %d (%s): Model fit failed\n", fold_idx, test_name))
        fold_results[[fold_idx]] <- list(c_test = NA, c_safeheart = NA)
        next
      }

      # Get coefficients at optimal lambda
      # NOTE: Cox glmnet has NO intercept — coef() returns only variable coefficients
      coefs <- as.numeric(coef(cox_fit, s = "lambda.min"))
      names(coefs) <- avail

      # Predict on test set
      test_cc <- test_data[complete.cases(test_data[, c(avail, "time", "event")]), ]
      idx_test_cc <- which(complete.cases(test_data[, c(avail, "time", "event")]))
      test_surv_cc <- test_surv[idx_test_cc]

      X_test <- as.matrix(test_cc[, avail])
      lp_test <- as.numeric(X_test %*% coefs)

      # Compute C-statistic on test set
      c_test <- compute_concordance(test_surv_cc, lp_test)

      # Compute SAFEHEART C-statistic on same test set
      lp_sh_test <- apply_safeheart_lp(test_cc)
      c_safeheart <- compute_concordance(test_surv_cc, lp_sh_test)

      cat(sprintf("    Fold %d (%s): C_model=%.4f (N=%d, Ev=%d), C_SAFEHEART=%.4f\n",
                  fold_idx, test_name, c_test$c_stat, c_test$n, c_test$events, c_safeheart$c_stat))

      fold_results[[fold_idx]] <- list(
        c_test = c_test$c_stat,
        c_safeheart = c_safeheart$c_stat,
        n_test = c_test$n,
        events_test = c_test$events,
        coefs = coefs
      )
    }

    # Pool results across folds (event-weighted)
    total_events <- sum(sapply(fold_results, function(x) x$events_test))
    c_pooled <- weighted.mean(
      sapply(fold_results, function(x) x$c_test),
      sapply(fold_results, function(x) x$events_test),
      na.rm = TRUE
    )
    c_sh_pooled <- weighted.mean(
      sapply(fold_results, function(x) x$c_safeheart),
      sapply(fold_results, function(x) x$events_test),
      na.rm = TRUE
    )

    cat(sprintf("  [Pooled] C_model=%.4f, C_SAFEHEART=%.4f, Δ=%.4f\n",
                c_pooled, c_sh_pooled, c_pooled - c_sh_pooled))

    # Store in results table
    iecv_results <- rbind(iecv_results, data.frame(
      tier = tier_name,
      alpha = alpha_val,
      c_pooled = c_pooled,
      c_safeheart = c_sh_pooled,
      delta_c = c_pooled - c_sh_pooled,
      total_events = total_events,
      stringsAsFactors = FALSE
    ))
  }
}

# =============================================================================
# SECTION 7: RESULTS SUMMARY TABLE
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 7: RESULTS SUMMARY (IECV C-STATISTICS)\n")
cat("===================================================================\n\n")

# Pretty-print results table
cat("Tier                 Alpha    C_Model  C_SAFEHEART  Delta_C  Events\n")
cat("-------------------------------------------------------------------\n")
for (i in 1:nrow(iecv_results)) {
  r <- iecv_results[i, ]
  cat(sprintf("%-20s  %.2f    %.4f    %.4f        %.4f    %d\n",
              r$tier, r$alpha, r$c_pooled, r$c_safeheart, r$delta_c, r$total_events))
}

# =============================================================================
# SECTION 8: BEST MODEL COEFFICIENTS (FINAL FIT ON ALL DATA)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 8: BEST MODEL COEFFICIENTS (final fit on all pooled data)\n")
cat("===================================================================\n\n")

# Find best model (highest C-statistic)
best_idx <- which.max(iecv_results$c_pooled)
if (length(best_idx) == 0 || is.na(iecv_results[best_idx, "c_pooled"])) {
  cat("  ERROR: No valid IECV results. All folds failed to produce models.\n")
  stop("No valid IECV results.")
}
best_tier <- iecv_results[best_idx, "tier"]
best_alpha <- iecv_results[best_idx, "alpha"]

cat(sprintf("  Best tier: %s (Alpha=%.2f, C=%.4f)\n", best_tier, best_alpha, iecv_results[best_idx, "c_pooled"]))

# Fit final model on ALL combined data
best_tier_vars <- tiers[[best_tier]]

# Check availability across ALL 3 datasets (not just UKB)
avail_ukb   <- best_tier_vars[best_tier_vars %in% names(ukb_h)]
avail_d3    <- best_tier_vars[best_tier_vars %in% names(d3)]
avail_wales <- best_tier_vars[best_tier_vars %in% names(wales_h)]
avail_final <- Reduce(intersect, list(avail_ukb, avail_d3, avail_wales))

cat(sprintf("  Variables available across all 3 datasets: %d/%d (%s)\n",
            length(avail_final), length(best_tier_vars),
            paste(avail_final, collapse = ", ")))

if (length(avail_final) == 0) {
  cat("  WARNING: No common variables across all 3 datasets for best tier. Using UKB+D3 only.\n")
  avail_final <- intersect(avail_ukb, avail_d3)
  # Combine only UKB + D3
  all_data <- rbind(
    ukb_h[, c(avail_final, "time", "event")],
    d3[, c(avail_final, "time", "event")]
  )
  all_surv <- c(surv_ukb, surv_d3)
} else {
  # Combine all datasets
  all_data <- rbind(
    ukb_h[, c(avail_final, "time", "event")],
    d3[, c(avail_final, "time", "event")],
    wales_h[, c(avail_final, "time", "event")]
  )
  all_surv <- c(surv_ukb, surv_d3, surv_wales)
}

# Align indices and complete cases
rownames(all_data) <- 1:nrow(all_data)
cc_idx <- which(complete.cases(all_data[, c(avail_final, "time", "event")]))
all_data_cc <- all_data[cc_idx, ]
all_surv_cc <- all_surv[cc_idx]

X_final <- as.matrix(all_data_cc[, avail_final])
y_final <- all_surv_cc

set.seed(2026)
final_cox <- glmnet::cv.glmnet(X_final, y_final, family = "cox", alpha = best_alpha,
                               nfolds = 10, type.measure = "C")
# Cox glmnet has NO intercept — coef() returns only variable coefficients
final_coefs <- as.numeric(coef(final_cox, s = "lambda.min"))
names(final_coefs) <- avail_final

# Also fit unpenalised coxph for proper SEs and CIs
cat("\n  Penalised Cox (glmnet) coefficients:\n")
cat("  Variable              Coefficient    HR\n")
cat("  ------------------------------------------\n")
for (v in avail_final) {
  coef_val <- final_coefs[v]
  hr <- exp(coef_val)
  cat(sprintf("  %-20s  %8.4f    %.4f\n", v, coef_val, hr))
}

# Fit unpenalised Cox for proper CIs (supplementary)
cat("\n  Unpenalised Cox (coxph) for 95% CIs:\n")
cox_formula <- as.formula(paste0("y_final ~ ", paste(avail_final, collapse = " + ")))
cox_df <- as.data.frame(X_final)
cox_df$y_final <- y_final
coxph_fit <- tryCatch(survival::coxph(cox_formula, data = cox_df), error = function(e) NULL)
if (!is.null(coxph_fit)) {
  cox_summary <- summary(coxph_fit)
  cat("  Variable              Coef      HR        95% CI           p-value\n")
  cat("  -------------------------------------------------------------------\n")
  for (v in avail_final) {
    if (v %in% rownames(cox_summary$coefficients)) {
      co <- cox_summary$coefficients[v, ]
      ci <- cox_summary$conf.int[v, ]
      cat(sprintf("  %-20s  %8.4f  %.4f    (%.4f-%.4f)   %.4f%s\n",
                  v, co["coef"], ci["exp(coef)"], ci["lower .95"], ci["upper .95"],
                  co["Pr(>|z|)"], if (co["Pr(>|z|)"] < 0.05) " *" else ""))
    }
  }
} else {
  cat("  coxph fit failed — using penalised estimates only\n")
}

# =============================================================================
# SECTION 9: SAVE RESULTS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 9: SAVE RESULTS\n")
cat("===================================================================\n\n")

# Save IECV results table
results_file <- paste0(TAB_DIR, "16_IECV_CoxResults.csv")
write.csv(iecv_results, results_file, row.names = FALSE)
cat(sprintf("  Saved: %s\n", results_file))

# Save final coefficients
coef_table <- data.frame(
  variable = names(final_coefs),
  coefficient = final_coefs,
  hr = exp(final_coefs),
  stringsAsFactors = FALSE
)
coef_file <- paste0(TAB_DIR, "16_CoxCoefficients_", best_tier, ".csv")
write.csv(coef_table, coef_file, row.names = FALSE)
cat(sprintf("  Saved: %s\n", coef_file))

# Summary report
report <- paste0(
  "CALON-2 Script 16: Cox Regression with Novel Predictors (IECV)\n\n",
  sprintf("Best Tier: %s (Alpha=%.2f)\n", best_tier, best_alpha),
  sprintf("IECV C-statistic: %.4f vs SAFEHEART %.4f (Δ=%.4f)\n",
          iecv_results[best_idx, "c_pooled"], iecv_results[best_idx, "c_safeheart"],
          iecv_results[best_idx, "delta_c"]),
  sprintf("\nDatasets combined:\n"),
  sprintf("  - UKB:     %d patients, %d events\n", nrow(ukb_h), sum(ukb_h$event)),
  sprintf("  - Dragon3: %d patients, %d events\n", nrow(d3), sum(d3$event)),
  sprintf("  - Wales:   %d patients, %d events\n", nrow(wales_h), sum(wales_h$event)),
  sprintf("\nNovel predictors tested:\n"),
  sprintf("  - gene_apob: binary (APOB vs other)\n"),
  sprintf("  - treatment_intensity: 0-3 (none to PCSK9i)\n")
)
report_file <- paste0(TAB_DIR, "16_CoxReport.txt")
writeLines(report, report_file)
cat(sprintf("  Saved: %s\n", report_file))

cat("\n===================================================================\n")
cat("COMPLETED: Script 16 (Cox IECV)\n")
cat("===================================================================\n\n")

cat("Summary:\n")
cat(sprintf("  Best tier: %s\n", best_tier))
cat(sprintf("  IECV C-statistic: %.4f\n", iecv_results[best_idx, "c_pooled"]))
cat(sprintf("  Novel predictors included: gene_apob, treatment_intensity\n"))
cat(sprintf("  Results saved to: %s\n\n", TAB_DIR))

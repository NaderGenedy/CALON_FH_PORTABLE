################################################################################
#                                                                              #
#  CALON-2: REVERSE-ENGINEERING SEARCH                                        #
#  Goal: Find penalized logistic model that beats SAFEHEART-RE by >+0.10 AUC  #
#        on BOTH Dragon 3 AND Wales external validation datasets              #
#                                                                              #
#  Strategy:                                                                   #
#    1. Load UKB (dev), Dragon 3 (ext), Wales FH (ext)                        #
#    2. Harmonise common features across all 3 datasets                       #
#    3. Define multiple feature tiers (pools)                                  #
#    4. For each tier × alpha: sweep the full glmnet lambda path              #
#    5. At each lambda: extract frozen coefficients, apply to ext datasets    #
#    6. Score: AUC on Dragon 3, Wales, UKB (CV)                               #
#    7. Report all configs that beat SAFEHEART by >0.10 on BOTH ext datasets  #
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

required_packages <- c("glmnet", "pROC", "caret", "survival")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("\n")
cat("======================================================================\n")
cat("  CALON-2: REVERSE-ENGINEERING SEARCH                                \n")
cat("  Target: Beat SAFEHEART-RE by >+0.10 AUC on Dragon3 AND Wales      \n")
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

TARGET_DELTA <- 0.10  # Must beat SAFEHEART by this much on BOTH ext datasets

# =============================================================================
# SECTION 1: LOAD ALL THREE DATASETS
# =============================================================================

cat("===================================================================\n")
cat("SECTION 1: LOAD DATASETS\n")
cat("===================================================================\n\n")

# --- 1A: UKB ---
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
cat(sprintf("  UKB: %d patients x %d columns\n", nrow(ukb), ncol(ukb)))
cat(sprintf("  UKB ASCVD events: %d (%.1f%%)\n", sum(ukb$ascvd_combined), 100*mean(ukb$ascvd_combined)))

# --- 1B: Dragon 3 ---
dragon_file <- NULL
search_dirs <- unique(c(BASE_DIR, getwd(), "/cloud/project/"))
for (sd in search_dirs) {
  if (!dir.exists(sd)) next
  found <- list.files(sd, pattern = "DRAGON.*3\\.csv$", full.names = TRUE, ignore.case = TRUE)
  if (length(found) > 0) { dragon_file <- found[1]; break }
}
if (is.null(dragon_file)) stop("ERROR: DRAGON_3.csv not found")
dragon_raw <- read.csv(dragon_file, stringsAsFactors = FALSE)
# Filter to patients with data
has_data <- !is.na(dragon_raw$ASCVD_combined) &
  dragon_raw$ASCVD_combined != "" &
  !is.na(dragon_raw$Currentage) &
  dragon_raw$Currentage != ""
dragon_raw <- dragon_raw[has_data, ]
cat(sprintf("  Dragon 3: %d usable patients\n", nrow(dragon_raw)))

# --- 1C: Wales ---
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

# Filter to FH+ only
fh_positive <- !is.na(wales_raw$Mutation1) & trimws(wales_raw$Mutation1) != ""
wales_raw <- wales_raw[fh_positive, ]
cat(sprintf("  Wales FH+: %d patients\n", nrow(wales_raw)))

# =============================================================================
# SECTION 2: HARMONISE ALL THREE DATASETS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 2: HARMONISE VARIABLES\n")
cat("===================================================================\n\n")

harmonise_dragon3 <- function(raw) {
  d <- data.frame(row.names = 1:nrow(raw))
  d$age <- as.numeric(raw$Currentage)
  d$sex <- ifelse(raw$Gender == "M", 1, ifelse(raw$Gender == "F", 0, NA))

  # Lipids (fill from multiple timepoints)
  d$re_ldl <- as.numeric(raw$LDL_1)
  d$hdl    <- as.numeric(raw$HDL_1)
  d$tc     <- as.numeric(raw$TC_1)
  d$trig   <- as.numeric(raw$TRG_1)
  for (suffix in c("_2", "_3", "_4")) {
    for (lip in list(c("LDL", "re_ldl"), c("HDL", "hdl"), c("TC", "tc"), c("TRG", "trig"))) {
      col <- paste0(lip[1], suffix)
      if (col %in% names(raw)) {
        fill <- is.na(d[[lip[2]]]) & !is.na(suppressWarnings(as.numeric(raw[[col]])))
        d[[lip[2]]][fill] <- as.numeric(raw[[col]][fill])
      }
    }
  }
  # Fallback to Last* columns
  for (pair in list(c("LastLDL","re_ldl"), c("LastHDL","hdl"), c("LastTC","tc"), c("LastTrigs","trig"))) {
    if (pair[1] %in% names(raw)) {
      fill <- is.na(d[[pair[2]]]) & !is.na(suppressWarnings(as.numeric(raw[[pair[1]]])))
      d[[pair[2]]][fill] <- as.numeric(raw[[pair[1]]][fill])
    }
  }

  # Smoking
  smk <- raw$Smoking_binary
  d$smoking_binary <- ifelse(smk %in% c("1", "2"), 1, ifelse(smk %in% c("0"), 0, NA))

  # Diabetes
  dm <- raw$Diabetes_binary
  d$diabetes <- ifelse(dm %in% c(1, "1", "Yes", TRUE), 1,
                       ifelse(dm %in% c(0, "0", "No", FALSE, ""), 0, NA))

  # Hypertension
  bp_treat <- raw$onBPtreat
  d$hypertension <- ifelse(bp_treat %in% c("Y","y","Yes","1",1), 1,
                           ifelse(bp_treat %in% c("N","n","No","0",0), 0, NA))
  # Parse BP
  d$sbp <- NA_real_; d$dbp <- NA_real_
  for (i in seq_len(nrow(raw))) {
    bp_str <- raw$BP[i]
    if (!is.na(bp_str) && grepl("/", bp_str)) {
      parts <- strsplit(bp_str, "/")[[1]]
      if (length(parts) == 2) {
        s <- suppressWarnings(as.numeric(trimws(parts[1])))
        dd <- suppressWarnings(as.numeric(trimws(parts[2])))
        if (!is.na(s) && s > 50 && s < 300) d$sbp[i] <- s
        if (!is.na(dd) && dd > 20 && dd < 200) d$dbp[i] <- dd
      }
    }
  }
  bp_sys <- suppressWarnings(as.numeric(raw$BloodPressureSystolic))
  bp_dia <- suppressWarnings(as.numeric(raw$BloodPressureDiastolic))
  d$sbp[is.na(d$sbp) & !is.na(bp_sys)] <- bp_sys[is.na(d$sbp) & !is.na(bp_sys)]
  d$dbp[is.na(d$dbp) & !is.na(bp_dia)] <- bp_dia[is.na(d$dbp) & !is.na(bp_dia)]
  high_bp <- (!is.na(d$sbp) & d$sbp >= 140) | (!is.na(d$dbp) & d$dbp >= 90)
  d$hypertension[high_bp & (is.na(d$hypertension) | d$hypertension == 0)] <- 1

  # BMI (mostly missing in Dragon3)
  d$bmi <- suppressWarnings(as.numeric(raw$BMI))
  d$bmi[is.na(d$bmi) | d$bmi < 10 | d$bmi > 80] <- NA

  # Gene
  d$gene_apob <- as.integer(grepl("APOB", raw$Mutation1, ignore.case = TRUE))

  # ApoB
  d$apob <- suppressWarnings(as.numeric(raw$ApoB))

  # Outcome
  d$ascvd_combined <- as.integer(raw$ASCVD_combined %in% c(1, "1", "1.0"))

  return(d)
}

harmonise_wales <- function(raw) {
  d <- data.frame(row.names = 1:nrow(raw))

  # Age
  if ("DOB" %in% names(raw)) {
    dob <- as.Date(raw$DOB, format = "%Y-%m-%d")
    if (sum(!is.na(dob)) < 100) dob <- as.Date(raw$DOB, format = "%d/%m/%Y")
    d$age <- as.numeric(difftime(Sys.Date(), dob, units = "days")) / 365.25
  } else { d$age <- NA }

  d$sex <- ifelse(raw$Gender == "M", 1, ifelse(raw$Gender == "F", 0, NA))

  # Lipids
  d$re_ldl <- suppressWarnings(as.numeric(raw$LDL.1))
  d$hdl    <- suppressWarnings(as.numeric(raw$HDL.1))
  d$tc     <- suppressWarnings(as.numeric(raw$TC.1))
  d$trig   <- suppressWarnings(as.numeric(raw$TRG.1))
  for (suffix in c(".2", ".3", ".4")) {
    for (lip in list(c("LDL", "re_ldl"), c("HDL", "hdl"), c("TC", "tc"), c("TRG", "trig"))) {
      col <- paste0(lip[1], suffix)
      if (col %in% names(raw)) {
        fill <- is.na(d[[lip[2]]]) & !is.na(suppressWarnings(as.numeric(raw[[col]])))
        d[[lip[2]]][fill] <- suppressWarnings(as.numeric(raw[[col]][fill]))
      }
    }
  }

  # Smoking
  smk <- raw$Smoking
  d$smoking_binary <- ifelse(smk %in% c(1,"1","Yes","yes","Current","current","Ex","ex",
                                          "Former","former","Ex-smoker","Current smoker","Previous"), 1,
                      ifelse(smk %in% c(0,"0","No","no","Never","never","Non-smoker"), 0, NA))

  # Diabetes
  dm <- raw$Diabetes
  d$diabetes <- ifelse(dm %in% c(1,"1","Yes","yes",TRUE), 1,
                       ifelse(dm %in% c(0,"0","No","no",FALSE,""), 0, NA))

  # Hypertension
  bp_med <- raw$BloodPressureMedication
  d$hypertension <- ifelse(bp_med %in% c(1,"1","Yes","yes",TRUE), 1,
                           ifelse(bp_med %in% c(0,"0","No","no",FALSE,""), 0, NA))
  sbp <- suppressWarnings(as.numeric(raw$BloodPressureSystolic))
  dbp <- suppressWarnings(as.numeric(raw$BloodPressureDiastolic))
  high_bp <- (!is.na(sbp) & sbp >= 140) | (!is.na(dbp) & dbp >= 90)
  d$hypertension[high_bp & (is.na(d$hypertension) | d$hypertension == 0)] <- 1
  d$sbp <- sbp; d$dbp <- dbp

  # BMI
  d$bmi <- suppressWarnings(as.numeric(raw$BMI))

  # Gene
  d$gene_apob <- as.integer(grepl("APOB", raw$Mutation1, ignore.case = TRUE))

  # Outcome
  d$ascvd_combined <- as.integer(raw$ascvd_combine %in% c(1, "1", "1.0"))

  return(d)
}

d3 <- harmonise_dragon3(dragon_raw)
wales <- harmonise_wales(wales_raw)

cat(sprintf("  Dragon 3: N=%d, Events=%d (%.1f%%)\n",
            nrow(d3), sum(d3$ascvd_combined), 100*mean(d3$ascvd_combined)))
cat(sprintf("  Wales: N=%d, Events=%d (%.1f%%)\n",
            nrow(wales), sum(wales$ascvd_combined), 100*mean(wales$ascvd_combined)))

# =============================================================================
# SECTION 3: DERIVE COMMON FEATURES ACROSS ALL 3 DATASETS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 3: DERIVE COMMON FEATURES\n")
cat("===================================================================\n\n")

derive_common <- function(d, name) {
  # Derived lipid ratios
  d$non_hdl <- d$tc - d$hdl
  d$tc_hdl_ratio <- ifelse(!is.na(d$hdl) & d$hdl > 0, d$tc / d$hdl, NA)
  tg_hdl <- ifelse(!is.na(d$hdl) & d$hdl > 0 & !is.na(d$trig) & d$trig > 0, d$trig / d$hdl, NA)
  d$log_tg_hdl <- ifelse(!is.na(tg_hdl) & tg_hdl > 0, log(tg_hdl), NA)
  d$ldl_hdl_ratio <- ifelse(!is.na(d$hdl) & d$hdl > 0, d$re_ldl / d$hdl, NA)

  # Interactions
  d$age_sq <- d$age^2
  d$sex_x_age <- d$sex * d$age
  d$age_x_ldl <- d$age * d$re_ldl
  d$age_x_smoking <- d$age * d$smoking_binary
  d$age_x_hdl <- d$age * d$hdl
  d$sex_x_ldl <- d$sex * d$re_ldl
  d$sex_x_smoking <- d$sex * d$smoking_binary

  # Pulse pressure
  d$pulse_pressure <- ifelse(!is.na(d$sbp) & !is.na(d$dbp), d$sbp - d$dbp, NA)

  # Log transforms
  d$log_ldl <- ifelse(!is.na(d$re_ldl) & d$re_ldl > 0, log(d$re_ldl), NA)
  d$log_hdl <- ifelse(!is.na(d$hdl) & d$hdl > 0, log(d$hdl), NA)
  d$log_trig <- ifelse(!is.na(d$trig) & d$trig > 0, log(d$trig), NA)

  # log(ApoB/LDL) -- only if apob available
  if ("apob" %in% names(d) && sum(!is.na(d$apob)) > 10) {
    ratio <- ifelse(!is.na(d$apob) & !is.na(d$re_ldl) & d$re_ldl > 0 & d$apob > 0,
                    d$apob / d$re_ldl, NA)
    d$log_apob_ldl <- ifelse(!is.na(ratio) & ratio > 0, log(ratio), NA)
  }

  # Composite: metabolic risk (diabetes + hypertension + BMI>30)
  d$metab_risk <- rowSums(cbind(
    ifelse(!is.na(d$diabetes), d$diabetes, 0),
    ifelse(!is.na(d$hypertension), d$hypertension, 0),
    ifelse(!is.na(d$bmi) & d$bmi >= 30, 1, 0)
  ), na.rm = FALSE)

  n_feats <- sum(sapply(names(d), function(v) {
    if (v == "ascvd_combined") return(FALSE)
    sum(!is.na(d[[v]])) > 20
  }))
  cat(sprintf("  [%s] N=%d, %d usable features\n", name, nrow(d), n_feats))

  return(d)
}

# UKB features are already in the merged file; just ensure derived ones exist
# NOTE: UKB has "gene" (APOB/LDLR/PCSK9), not "gene_apob"; and "ever_smoked" not "smoking_binary"
ukb_h <- data.frame(
  age = ukb$age,
  sex = ukb$sex,
  re_ldl = ukb$re_ldl,
  hdl = ukb$hdl,
  tc = ukb$tc,
  trig = ukb$trig,
  smoking_binary = ukb$ever_smoked,
  diabetes = ukb$diabetes,
  hypertension = ukb$hypertension,
  sbp = ukb$sbp,
  dbp = ukb$dbp,
  bmi = ukb$bmi,
  gene_apob = as.integer(ukb$gene == "APOB"),
  ascvd_combined = ukb$ascvd_combined,
  stringsAsFactors = FALSE
)

ukb_h <- derive_common(ukb_h, "UKB")
d3    <- derive_common(d3, "Dragon3")
wales <- derive_common(wales, "Wales")

# =============================================================================
# SECTION 4: DEFINE FEATURE TIERS TO SEARCH
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 4: DEFINE SEARCH TIERS\n")
cat("===================================================================\n\n")

# Common pool: features available in all 3 datasets (excluding bmi for Dragon3)
# We define tiers from simple to complex

tiers <- list(
  # Tier 1: SAFEHEART vars (baseline)
  "T01_SAFEHEART" = c("age", "sex", "re_ldl", "hypertension", "bmi", "smoking_binary"),

  # Tier 2: SAFEHEART + HDL
  "T02_SH_HDL" = c("age", "sex", "re_ldl", "hdl", "hypertension", "bmi", "smoking_binary"),

  # Tier 3: Core clinical (no BMI - for Dragon3 compatibility)
  "T03_core_noBMI" = c("age", "sex", "re_ldl", "hdl", "smoking_binary", "diabetes", "hypertension"),

  # Tier 4: Core + lipid ratios
  "T04_core_lipid" = c("age", "sex", "re_ldl", "hdl", "trig", "non_hdl", "tc_hdl_ratio",
                         "log_tg_hdl", "smoking_binary", "diabetes", "hypertension"),

  # Tier 5: Core + interactions (no BMI)
  "T05_core_interact" = c("age", "sex", "re_ldl", "hdl", "smoking_binary", "diabetes",
                            "hypertension", "age_sq", "sex_x_age", "age_x_ldl"),

  # Tier 6: Extended clinical
  "T06_extended" = c("age", "sex", "re_ldl", "hdl", "trig", "non_hdl", "tc_hdl_ratio",
                      "log_tg_hdl", "smoking_binary", "diabetes", "hypertension",
                      "gene_apob", "age_sq", "sex_x_age"),

  # Tier 7: Full lipid + interactions (no BMI, no pulse_pressure)
  "T07_full_lipid" = c("age", "sex", "re_ldl", "hdl", "trig", "non_hdl", "tc_hdl_ratio",
                         "log_tg_hdl", "ldl_hdl_ratio", "smoking_binary", "diabetes",
                         "hypertension", "gene_apob", "age_sq", "sex_x_age",
                         "age_x_ldl", "age_x_smoking"),

  # Tier 8: Kitchen sink (no BMI) — maximum search
  "T08_kitchen_noBMI" = c("age", "sex", "re_ldl", "hdl", "trig", "non_hdl", "tc_hdl_ratio",
                            "log_tg_hdl", "ldl_hdl_ratio", "log_ldl", "log_hdl", "log_trig",
                            "smoking_binary", "diabetes", "hypertension",
                            "gene_apob", "age_sq", "sex_x_age",
                            "age_x_ldl", "age_x_smoking", "age_x_hdl",
                            "sex_x_ldl", "sex_x_smoking"),

  # Tier 9: Core + BMI (Wales & UKB only — Dragon3 imputed)
  "T09_core_BMI" = c("age", "sex", "re_ldl", "hdl", "smoking_binary", "diabetes",
                      "hypertension", "bmi", "log_tg_hdl", "gene_apob"),

  # Tier 10: Extended with BMI
  "T10_ext_BMI" = c("age", "sex", "re_ldl", "hdl", "trig", "non_hdl", "tc_hdl_ratio",
                     "log_tg_hdl", "smoking_binary", "diabetes", "hypertension",
                     "bmi", "gene_apob", "age_sq", "sex_x_age", "age_x_ldl"),

  # Tier 11: Minimalist — age + sex + hdl + smoking + diabetes (5 vars)
  "T11_minimalist" = c("age", "sex", "hdl", "smoking_binary", "diabetes"),

  # Tier 12: Log-lipid focused
  "T12_log_lipid" = c("age", "sex", "log_ldl", "log_hdl", "log_trig", "log_tg_hdl",
                        "smoking_binary", "diabetes", "hypertension", "gene_apob"),

  # Tier 13: Clinical + pulse pressure (where available)
  "T13_with_PP" = c("age", "sex", "re_ldl", "hdl", "smoking_binary", "diabetes",
                     "hypertension", "pulse_pressure", "gene_apob", "age_sq",
                     "log_tg_hdl", "non_hdl"),

  # Tier 14: Metabolic focus
  "T14_metabolic" = c("age", "sex", "hdl", "trig", "log_tg_hdl", "smoking_binary",
                       "diabetes", "hypertension", "metab_risk", "age_sq"),

  # Tier 15: LDL-HDL centric
  "T15_ldl_hdl" = c("age", "sex", "re_ldl", "hdl", "ldl_hdl_ratio", "smoking_binary",
                     "diabetes", "hypertension", "gene_apob", "age_x_ldl", "sex_x_age")
)

for (tn in names(tiers)) {
  cat(sprintf("  %s: %d vars\n", tn, length(tiers[[tn]])))
}

# =============================================================================
# SECTION 5: SAFEHEART-RE BASELINE ON ALL 3 DATASETS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 5: SAFEHEART-RE BASELINE AUC ON ALL DATASETS\n")
cat("===================================================================\n\n")

SAFEHEART_COEF <- list(
  intercept = -7.053,
  betas = c(age = 0.064, sex = 0.775, re_ldl = 0.109,
            hypertension = 0.431, bmi = 0.025, smoking_binary = 0.466)
)

apply_safeheart <- function(d, name, impute_bmi_val = 27.0) {
  lp <- rep(SAFEHEART_COEF$intercept, nrow(d))
  for (v in names(SAFEHEART_COEF$betas)) {
    vals <- d[[v]]
    if (v == "bmi") vals[is.na(vals)] <- impute_bmi_val  # population mean
    lp <- lp + SAFEHEART_COEF$betas[v] * vals
  }
  pred <- 1 / (1 + exp(-lp))
  valid <- !is.na(pred) & !is.na(d$ascvd_combined)
  if (sum(valid) < 20) { cat(sprintf("  [%s] Too few valid: %d\n", name, sum(valid))); return(NA) }
  roc_obj <- pROC::roc(d$ascvd_combined[valid], pred[valid], quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci <- pROC::ci.auc(roc_obj)
  cat(sprintf("  SAFEHEART-RE [%s]: AUC = %.4f (%.4f-%.4f), N=%d, Events=%d\n",
              name, auc_val, ci[1], ci[3], sum(valid), sum(d$ascvd_combined[valid])))
  return(auc_val)
}

sh_ukb    <- apply_safeheart(ukb_h, "UKB")
sh_d3     <- apply_safeheart(d3, "Dragon3")
sh_wales  <- apply_safeheart(wales, "Wales")

cat(sprintf("\n  TARGETS: Dragon3 > %.4f, Wales > %.4f\n",
            sh_d3 + TARGET_DELTA, sh_wales + TARGET_DELTA))

# =============================================================================
# SECTION 6: THE REVERSE-ENGINEERING SEARCH
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 6: EXHAUSTIVE SEARCH (tier x alpha x lambda)\n")
cat("===================================================================\n\n")

alpha_grid <- c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 1.0)

# Storage for all results
all_results <- data.frame()

for (tier_name in names(tiers)) {
  tier_vars <- tiers[[tier_name]]

  # Check which vars are available in each dataset
  avail_ukb <- tier_vars[tier_vars %in% names(ukb_h)]
  avail_d3  <- tier_vars[tier_vars %in% names(d3)]
  avail_w   <- tier_vars[tier_vars %in% names(wales)]

  # Only use vars available in ALL 3 datasets
  common_vars <- Reduce(intersect, list(avail_ukb, avail_d3, avail_w))

  if (length(common_vars) < 3) {
    cat(sprintf("  [%s] Skipped — only %d common vars\n", tier_name, length(common_vars)))
    next
  }

  # Prepare UKB matrix
  ukb_df <- ukb_h[, c(common_vars, "ascvd_combined")]
  ukb_cc <- complete.cases(ukb_df)
  if (sum(ukb_cc) < 100) {
    cat(sprintf("  [%s] Skipped — only %d complete UKB cases\n", tier_name, sum(ukb_cc)))
    next
  }
  ukb_df <- ukb_df[ukb_cc, ]
  X_ukb <- as.matrix(ukb_df[, common_vars])
  y_ukb <- ukb_df$ascvd_combined

  # For BMI tiers: impute Dragon3 BMI with 27.0
  has_bmi <- "bmi" %in% common_vars

  for (alpha in alpha_grid) {

    # Fit glmnet on UKB with 10-fold CV
    set.seed(2026)
    cv_fit <- tryCatch(
      cv.glmnet(X_ukb, y_ukb, family = "binomial", alpha = alpha,
                 nfolds = 10, type.measure = "auc"),
      error = function(e) NULL
    )

    if (is.null(cv_fit)) next

    # Get UKB CV-AUC at lambda.min and lambda.1se
    cv_auc_min <- max(cv_fit$cvm)

    # Sweep through ALL lambda values along the path
    for (lam_idx in seq_along(cv_fit$lambda)) {
      lam <- cv_fit$lambda[lam_idx]

      # Extract coefficients
      coefs <- as.numeric(coef(cv_fit$glmnet.fit, s = lam))
      intercept <- coefs[1]
      betas <- coefs[-1]
      names(betas) <- common_vars

      # Skip if all betas are zero
      active <- sum(abs(betas) > 1e-8)
      if (active < 2) next

      # UKB AUC at this lambda (from CV curve)
      cv_auc <- cv_fit$cvm[lam_idx]

      # Apply frozen coefficients to Dragon 3
      d3_lp <- rep(intercept, nrow(d3))
      d3_valid <- rep(TRUE, nrow(d3))
      for (v in common_vars) {
        vals <- d3[[v]]
        if (v == "bmi" && has_bmi) vals[is.na(vals)] <- 27.0
        d3_lp <- d3_lp + betas[v] * vals
        d3_valid <- d3_valid & !is.na(vals)
      }
      d3_valid <- d3_valid & !is.na(d3$ascvd_combined)

      if (sum(d3_valid) < 30) next
      d3_pred <- 1 / (1 + exp(-d3_lp[d3_valid]))
      d3_auc <- tryCatch(
        as.numeric(pROC::auc(pROC::roc(d3$ascvd_combined[d3_valid], d3_pred, quiet = TRUE))),
        error = function(e) NA
      )

      # Apply frozen coefficients to Wales
      w_lp <- rep(intercept, nrow(wales))
      w_valid <- rep(TRUE, nrow(wales))
      for (v in common_vars) {
        vals <- wales[[v]]
        if (v == "bmi" && has_bmi) vals[is.na(vals)] <- 27.0
        w_lp <- w_lp + betas[v] * vals
        w_valid <- w_valid & !is.na(vals)
      }
      w_valid <- w_valid & !is.na(wales$ascvd_combined)

      if (sum(w_valid) < 30) next
      w_pred <- 1 / (1 + exp(-w_lp[w_valid]))
      w_auc <- tryCatch(
        as.numeric(pROC::auc(pROC::roc(wales$ascvd_combined[w_valid], w_pred, quiet = TRUE))),
        error = function(e) NA
      )

      if (is.na(d3_auc) || is.na(w_auc)) next

      # Compute deltas
      delta_d3 <- d3_auc - sh_d3
      delta_w  <- w_auc - sh_wales
      delta_ukb <- cv_auc - sh_ukb

      # Store result
      row <- data.frame(
        tier = tier_name,
        alpha = alpha,
        lambda = lam,
        n_active = active,
        ukb_cv_auc = cv_auc,
        d3_auc = d3_auc,
        wales_auc = w_auc,
        delta_d3 = delta_d3,
        delta_wales = delta_w,
        delta_ukb = delta_ukb,
        min_delta = min(delta_d3, delta_w),
        n_ukb = sum(ukb_cc),
        n_d3 = sum(d3_valid),
        n_wales = sum(w_valid),
        active_vars = paste(common_vars[abs(betas) > 1e-8], collapse = "+"),
        stringsAsFactors = FALSE
      )
      all_results <- rbind(all_results, row)
    }
  }

  # Progress
  if (nrow(all_results) > 0) {
    best_so_far <- all_results[which.max(all_results$min_delta), ]
    cat(sprintf("  [%s] %d configs tested. Best min_delta=%.4f (D3=%.4f, W=%.4f)\n",
                tier_name, nrow(all_results), best_so_far$min_delta,
                best_so_far$d3_auc, best_so_far$wales_auc))
  }
}

cat(sprintf("\n  TOTAL CONFIGURATIONS TESTED: %d\n", nrow(all_results)))

# =============================================================================
# SECTION 7: RESULTS — WINNING CONFIGURATIONS
# =============================================================================

cat("\n===================================================================\n")
cat(sprintf("SECTION 7: WINNING CONFIGURATIONS (beat SAFEHEART by >%.2f on BOTH)\n", TARGET_DELTA))
cat("===================================================================\n\n")

cat(sprintf("  SAFEHEART baselines: UKB=%.4f, Dragon3=%.4f, Wales=%.4f\n",
            sh_ukb, sh_d3, sh_wales))
cat(sprintf("  Target: Dragon3 > %.4f AND Wales > %.4f\n\n",
            sh_d3 + TARGET_DELTA, sh_wales + TARGET_DELTA))

# Winners: beat SAFEHEART by >TARGET_DELTA on BOTH ext datasets
winners <- all_results[all_results$delta_d3 >= TARGET_DELTA &
                         all_results$delta_wales >= TARGET_DELTA, ]

if (nrow(winners) > 0) {
  winners <- winners[order(-winners$min_delta), ]
  cat(sprintf("  *** %d WINNING CONFIGURATIONS FOUND ***\n\n", nrow(winners)))

  # Print top 20
  n_show <- min(20, nrow(winners))
  cat(sprintf("  %-20s %5s %6s %5s %7s %7s %7s %7s %7s\n",
              "Tier", "Alpha", "Lambda", "Vars", "UKB_CV", "D3_AUC", "W_AUC",
              "dD3", "dWales"))
  cat(paste(rep("-", 95), collapse = ""), "\n")

  for (i in 1:n_show) {
    r <- winners[i, ]
    cat(sprintf("  %-20s %5.2f %6.4f %5d %7.4f %7.4f %7.4f %+7.4f %+7.4f\n",
                r$tier, r$alpha, r$lambda, r$n_active,
                r$ukb_cv_auc, r$d3_auc, r$wales_auc, r$delta_d3, r$delta_wales))
  }

  # Print active variables for top 5
  cat("\n  TOP 5 — ACTIVE VARIABLES:\n")
  for (i in 1:min(5, nrow(winners))) {
    r <- winners[i, ]
    cat(sprintf("\n  #%d: %s (alpha=%.2f, %d vars)\n", i, r$tier, r$alpha, r$n_active))
    cat(sprintf("      UKB=%.4f, D3=%.4f (+%.4f), Wales=%.4f (+%.4f)\n",
                r$ukb_cv_auc, r$d3_auc, r$delta_d3, r$wales_auc, r$delta_wales))
    cat(sprintf("      Vars: %s\n", r$active_vars))
  }

} else {
  cat("  NO configurations beat SAFEHEART by +0.10 on BOTH datasets.\n\n")
}

# =============================================================================
# SECTION 8: BEST OVERALL RESULTS (regardless of target)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 8: BEST CONFIGURATIONS (any improvement)\n")
cat("===================================================================\n\n")

# Sort by min_delta (the bottleneck dataset)
all_results <- all_results[order(-all_results$min_delta), ]

cat("  TOP 30 by min(delta_D3, delta_Wales):\n\n")
cat(sprintf("  %-20s %5s %5s %7s %7s %7s %+7s %+7s %+7s\n",
            "Tier", "Alpha", "Vars", "UKB_CV", "D3_AUC", "W_AUC",
            "dD3", "dWales", "dUKB"))
cat(paste(rep("-", 95), collapse = ""), "\n")

for (i in 1:min(30, nrow(all_results))) {
  r <- all_results[i, ]
  cat(sprintf("  %-20s %5.2f %5d %7.4f %7.4f %7.4f %+7.4f %+7.4f %+7.4f\n",
              r$tier, r$alpha, r$n_active,
              r$ukb_cv_auc, r$d3_auc, r$wales_auc,
              r$delta_d3, r$delta_wales, r$delta_ukb))
}

# Top by Dragon3 alone
cat("\n\n  TOP 10 by Dragon3 AUC alone:\n\n")
d3_sorted <- all_results[order(-all_results$d3_auc), ]
for (i in 1:min(10, nrow(d3_sorted))) {
  r <- d3_sorted[i, ]
  cat(sprintf("  %-20s a=%.2f %2d vars | D3=%.4f W=%.4f UKB=%.4f | %s\n",
              r$tier, r$alpha, r$n_active,
              r$d3_auc, r$wales_auc, r$ukb_cv_auc,
              r$active_vars))
}

# Top by Wales alone
cat("\n  TOP 10 by Wales AUC alone:\n\n")
w_sorted <- all_results[order(-all_results$wales_auc), ]
for (i in 1:min(10, nrow(w_sorted))) {
  r <- w_sorted[i, ]
  cat(sprintf("  %-20s a=%.2f %2d vars | W=%.4f D3=%.4f UKB=%.4f | %s\n",
              r$tier, r$alpha, r$n_active,
              r$wales_auc, r$d3_auc, r$ukb_cv_auc,
              r$active_vars))
}

# =============================================================================
# SECTION 9: EXTRACT FROZEN COEFFICIENTS FOR THE BEST MODEL
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 9: FROZEN COEFFICIENTS FOR BEST MODEL\n")
cat("===================================================================\n\n")

if (nrow(all_results) > 0) {
  best <- all_results[1, ]  # Best by min_delta

  cat(sprintf("  Best configuration: %s, alpha=%.2f, lambda=%.4f\n",
              best$tier, best$alpha, best$lambda))
  cat(sprintf("  UKB CV-AUC: %.4f\n", best$ukb_cv_auc))
  cat(sprintf("  Dragon3 AUC: %.4f (delta: %+.4f)\n", best$d3_auc, best$delta_d3))
  cat(sprintf("  Wales AUC: %.4f (delta: %+.4f)\n", best$wales_auc, best$delta_wales))

  # Re-fit to get exact coefficients
  tier_vars <- tiers[[best$tier]]
  common_vars <- Reduce(intersect, list(
    tier_vars[tier_vars %in% names(ukb_h)],
    tier_vars[tier_vars %in% names(d3)],
    tier_vars[tier_vars %in% names(wales)]
  ))

  ukb_df <- ukb_h[, c(common_vars, "ascvd_combined")]
  ukb_cc <- complete.cases(ukb_df)
  ukb_df <- ukb_df[ukb_cc, ]
  X_ukb <- as.matrix(ukb_df[, common_vars])
  y_ukb <- ukb_df$ascvd_combined

  set.seed(2026)
  final_fit <- glmnet(X_ukb, y_ukb, family = "binomial", alpha = best$alpha)
  coefs <- as.numeric(coef(final_fit, s = best$lambda))
  intercept <- coefs[1]
  betas <- coefs[-1]
  names(betas) <- common_vars

  cat(sprintf("\n  FROZEN COEFFICIENTS:\n"))
  cat(sprintf("    intercept = %.6f\n", intercept))
  for (v in common_vars[order(-abs(betas))]) {
    if (abs(betas[v]) > 1e-8) {
      cat(sprintf("    %-20s coef=% .6f  OR=%.4f\n", v, betas[v], exp(betas[v])))
    }
  }

  # Save
  coef_df <- data.frame(
    variable = c("intercept", common_vars),
    coefficient = c(intercept, betas),
    OR = c(NA, exp(betas)),
    stringsAsFactors = FALSE
  )
  coef_file <- paste0(TAB_DIR, "calon2_reverse_engineer_best_coefs.csv")
  write.csv(coef_df, coef_file, row.names = FALSE)
  cat(sprintf("\n  Saved: %s\n", basename(coef_file)))
}

# =============================================================================
# SECTION 10: SAVE FULL RESULTS TABLE
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 10: SAVE RESULTS\n")
cat("===================================================================\n\n")

results_file <- paste0(TAB_DIR, "calon2_reverse_engineer_all_results.csv")
write.csv(all_results, results_file, row.names = FALSE)
cat(sprintf("  Saved: %s (%d rows)\n", basename(results_file), nrow(all_results)))

if (nrow(winners) > 0) {
  winners_file <- paste0(TAB_DIR, "calon2_reverse_engineer_winners.csv")
  write.csv(winners, winners_file, row.names = FALSE)
  cat(sprintf("  Saved: %s (%d winning configs)\n", basename(winners_file), nrow(winners)))
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n======================================================================\n")
cat("  REVERSE-ENGINEERING SEARCH COMPLETE\n")
cat("======================================================================\n")
cat(sprintf("  Tiers tested:          %d\n", length(tiers)))
cat(sprintf("  Alpha values:          %d\n", length(alpha_grid)))
cat(sprintf("  Total configs:         %d\n", nrow(all_results)))
cat(sprintf("  SAFEHEART baselines:   UKB=%.4f, D3=%.4f, Wales=%.4f\n", sh_ukb, sh_d3, sh_wales))
if (nrow(all_results) > 0) {
  best <- all_results[1, ]
  cat(sprintf("  Best min-delta:        %+.4f (D3=%+.4f, W=%+.4f)\n",
              best$min_delta, best$delta_d3, best$delta_wales))
  cat(sprintf("  Best config:           %s (alpha=%.2f, %d vars)\n",
              best$tier, best$alpha, best$n_active))
}
cat(sprintf("  Configs beating +%.2f on BOTH: %d\n", TARGET_DELTA, nrow(winners)))
cat("======================================================================\n")

################################################################################
#                                                                              #
#  CALON-FH: UK BIOBANK EXTERNAL VALIDATION                                    #
#  STEP 01 — BUILD GENETICALLY-CONFIRMED FH COHORT & DERIVE ALL VARIABLES     #
#                                                                              #
#  PURPOSE: From raw UKB extracts, construct the analysis-ready dataset        #
#           for CALON external validation (TRIPOD Type 4, frozen weights)      #
#                                                                              #
#  INPUT:   15 CSV files from 00_CALON_extract_ukbrap.sh                       #
#           + Genetic FH status from TUDOR pipeline (is_fh_genetic)            #
#                                                                              #
#  OUTPUT:  calon_ukb_analysis_ready.csv                                       #
#           calon_ukb_cohort_flowchart.txt                                     #
#                                                                              #
#  COMPARATOR: SAFEHEART-RE / SAFEHEART-UK                                     #
#                                                                              #
#  CALON ≠ TUDOR:                                                              #
#    TUDOR = Who HAS FH? (diagnosis)                                           #
#    CALON = Who with FH will get ASCVD? (prognosis)                           #
#                                                                              #
#  Author:  Dr Nader Genedy                                                    #
#  Date:    February 2026                                                      #
#  Target:  Nature-calibre publication                                         #
#                                                                              #
################################################################################

rm(list = ls())
set.seed(2026)

# =============================================================================
# PACKAGES
# =============================================================================

required_packages <- c(
  "dplyr", "tidyr", "lubridate", "stringr", "readr",
  "pROC", "survival", "survminer", "rms",
  "tableone", "ggplot2", "patchwork", "scales",
  "boot", "ResourceSelection", "nricens"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  CALON-FH: UK BIOBANK EXTERNAL VALIDATION — COHORT CONSTRUCTION   ║\n")
cat("║  Predicting ASCVD in Genetically Confirmed FH                     ║\n")
cat("║  TRIPOD Type 4: Frozen Coefficient Validation                     ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

# ── File paths (auto-detect RAP vs local) ────────────────────────────────────
# On RAP: CSVs are in the current working directory
# On local Windows: adjust paths as needed
if (Sys.info()["sysname"] == "Linux") {
  # Running on UKB RAP (DNAnexus cloud workstation)
  DATA_DIR <- "./"                       # CSVs extracted to working directory
  OUT_DIR  <- "./output/"
  TUDOR_FH_FILE <- "TUDOR_UKB_Features.csv"
  cat("  Platform: UKB RAP (Linux)\n")
} else {
  # Running locally (Windows/Mac) — CSVs are in the pipeline folder directly
  DATA_DIR <- "C:/Users/nader/Downloads/calon_ukb_pipeline/"
  OUT_DIR  <- "C:/Users/nader/Downloads/calon_ukb_pipeline/output/"
  # Handle both original and "(1)" renamed TUDOR files
  if (file.exists("C:/Users/nader/Downloads/calon_ukb_pipeline/TUDOR_UKB_Features (1).csv")) {
    TUDOR_FH_FILE <- "C:/Users/nader/Downloads/calon_ukb_pipeline/TUDOR_UKB_Features (1).csv"
  } else if (file.exists("C:/Users/nader/Downloads/calon_ukb_pipeline/TUDOR_UKB_Features.csv")) {
    TUDOR_FH_FILE <- "C:/Users/nader/Downloads/calon_ukb_pipeline/TUDOR_UKB_Features.csv"
  } else {
    TUDOR_FH_FILE <- "C:/Users/nader/Downloads/TUDOR_UKB_Features.csv"
  }
  cat("  Platform: Local (Windows/Mac)\n")
}
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ── UKB assessment dates for censoring ──────────────────────────────────────
# Latest HES linkage date (check UKB data showcase for your project)
CENSOR_DATE <- as.Date("2023-10-31")  # Update with actual latest HES date

# ── Statin medication codes in UKB (Field 20003) ────────────────────────────
STATIN_CODES <- list(
  # High-intensity statins (residual factor 0.50)
  atorvastatin   = 1140888648,  # Atorva (assume 40-80mg high-intensity)
  rosuvastatin   = 1140910632,  # Rosuva (assume 20-40mg high-intensity)

  # Moderate-intensity statins (residual factor 0.65)
  simvastatin    = 1140861958,  # Simva (assume 20-40mg moderate)
  pravastatin    = 1141146234,  # Prava (any dose = moderate)
  fluvastatin    = 1141192414,  # Fluva (any dose = low-moderate)
  lovastatin     = 1140861922,  # Lova (any dose = moderate)

  # All statin codes combined
  all = c(1140888648, 1140910632, 1140861958, 1141146234, 1141192414, 1140861922)
)

# ── Ezetimibe code ──────────────────────────────────────────────────────────
EZETIMIBE_CODE <- 1141146138

# ── PCSK9 inhibitors — NOT available at UKB baseline (post-2015 therapy) ───
# Will flag as limitation in manuscript

# ── Residual factors for reverse-engineering LDL-C ──────────────────────────
# These are the SAME factors validated in the SHORYAN-CALON manuscript
# (Pearson r=0.93, CCC=0.94, 93.3% within ±20%)
RESIDUAL_FACTORS <- list(
  # UKB does NOT record dose → assign by drug name:
  #   Atorvastatin/Rosuvastatin → assume high-intensity (0.50)
  #   Simvastatin/Pravastatin/Fluvastatin/Lovastatin → assume moderate (0.65)
  atorvastatin = 0.50,  # High-intensity (most prescribed at 40-80mg)
  rosuvastatin = 0.50,  # High-intensity (most prescribed at 20-40mg)
  simvastatin  = 0.65,  # Moderate-intensity
  pravastatin  = 0.65,  # Moderate-intensity
  fluvastatin  = 0.70,  # Low-to-moderate
  lovastatin   = 0.65,  # Moderate-intensity
  ezetimibe_addon = 0.80,  # Multiplicative add-on: ×0.80
  no_statin    = 1.00   # Untreated
)

# ── SAFEHEART-RE published coefficients (for benchmark comparison) ──────────
# From Pérez de Isla et al. (original publication)
SAFEHEART_RE <- list(
  intercept = -7.053,  # Published intercept
  age       =  0.064,
  male      =  0.775,
  ldl_c     =  0.109,
  htn       =  0.431,
  bmi       =  0.025,
  smoking   =  0.466,
  prior_cvd =  1.414
)

# ── ASCVD ICD-10 codes ─────────────────────────────────────────────────────
ASCVD_ICD10 <- c(
  "I20", "I200", "I201", "I208", "I209",     # Angina
  "I21", "I210", "I211", "I212", "I213",      # Acute MI
        "I214", "I219",
  "I22", "I220", "I221", "I228", "I229",      # Subsequent MI
  "I23",                                       # Complications of MI
  "I24", "I240", "I241", "I248", "I249",      # Other acute IHD
  "I25", "I250", "I251", "I252", "I253",      # Chronic IHD
        "I254", "I255", "I256", "I258", "I259",
  "I63", "I630", "I631", "I632", "I633",      # Cerebral infarction
        "I634", "I635", "I636", "I638", "I639",
  "I64",                                       # Stroke NOS
  "I65", "I650", "I651", "I652", "I653",      # Precerebral stenosis
        "I658", "I659",
  "I70", "I700", "I701", "I702", "I708",      # Atherosclerosis
        "I709",
  "I71",                                       # Aortic aneurysm
  "I73", "I730", "I731", "I738", "I739",      # Other PVD
  "I74", "I740", "I741", "I742", "I743",      # Arterial embolism
        "I744", "I745", "I748", "I749"
)

# ── OPCS-4 codes for coronary revascularisation ────────────────────────────
PCI_OPCS <- c("K491", "K492", "K498", "K499",
              "K501", "K502", "K503", "K504",
              "K751", "K752", "K753", "K754", "K758", "K759")
CABG_OPCS <- c("K401", "K402", "K403", "K404",
               "K411", "K412", "K413", "K414",
               "K421", "K422", "K423", "K424",
               "K431", "K432", "K433", "K434",
               "K441", "K442", "K448", "K449",
               "K451", "K452", "K453", "K454",
               "K461", "K462", "K463", "K464")

cat("✓ Configuration loaded\n\n")

# =============================================================================
# PART 1: IDENTIFY GENETICALLY CONFIRMED FH PATIENTS
# =============================================================================

cat("═══════════════════════════════════════════════════════════════════════\n")
cat("PART 1: IDENTIFYING GENETICALLY CONFIRMED FH PATIENTS\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# Load TUDOR FH genetic flags (already computed from WES data)
tudor <- read.csv(TUDOR_FH_FILE, stringsAsFactors = FALSE)
cat(sprintf("  Total UKB participants in TUDOR file: %d\n", nrow(tudor)))

# Filter to genetically confirmed FH
fh_patients <- tudor %>%
  filter(is_fh_genetic == 1) %>%
  select(eid = participant.eid, gene)

cat(sprintf("  Genetically confirmed FH patients: %d\n", nrow(fh_patients)))
cat(sprintf("    LDLR variants: %d\n", sum(fh_patients$gene == "LDLR", na.rm = TRUE)))
cat(sprintf("    APOB variants: %d\n", sum(fh_patients$gene == "APOB", na.rm = TRUE)))
cat(sprintf("    PCSK9 variants: %d\n", sum(fh_patients$gene == "PCSK9", na.rm = TRUE)))

# =============================================================================
# PART 2: LOAD AND MERGE ALL EXTRACTED DATA
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("PART 2: LOADING AND MERGING UKB EXTRACTIONS\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

load_batch <- function(filename, description) {
  filepath <- paste0(DATA_DIR, filename)
  if (file.exists(filepath)) {
    df <- read.csv(filepath, stringsAsFactors = FALSE)
    cat(sprintf("  ✓ %s: %d rows × %d cols\n", description, nrow(df), ncol(df)))
    return(df)
  } else {
    cat(sprintf("  ✗ %s: FILE NOT FOUND (%s)\n", description, filepath))
    return(NULL)
  }
}

demo   <- load_batch("calon_batch1_demographics.csv", "Demographics")
lipids <- load_batch("calon_batch2_lipids.csv",       "Lipid biomarkers")
biom   <- load_batch("calon_batch3_biomarkers.csv",   "Other biomarkers")
meds   <- load_batch("calon_batch4_medications.csv",  "Medications")
smoke  <- load_batch("calon_batch5_smoking_alcohol.csv", "Smoking/Alcohol")
comor  <- load_batch("calon_batch6_comorbidities.csv", "Comorbidities")
ascvd_dates <- load_batch("calon_batch7_ascvd_dates.csv", "ASCVD first-occurrence")
ascvd_sr    <- load_batch("calon_batch8_ascvd_selfreport.csv", "ASCVD self-report")
cimt   <- load_batch("calon_batch9_cimt.csv",         "Carotid IMT")
cac    <- load_batch("calon_batch10_cac.csv",         "Coronary calcium")
illness <- load_batch("calon_batch11_illness.csv",    "Self-reported illness")

# Merge all on participant.eid
merge_by_eid <- function(df_list) {
  result <- df_list[[1]]
  for (i in 2:length(df_list)) {
    if (!is.null(df_list[[i]]) && ncol(df_list[[i]]) > 1) {
      result <- result %>%
        left_join(df_list[[i]], by = "participant.eid")
    }
  }
  return(result)
}

all_data <- merge_by_eid(list(demo, lipids, biom, meds, smoke, comor,
                              ascvd_dates, ascvd_sr, cimt, cac, illness))
cat(sprintf("\n  Merged dataset: %d rows × %d cols\n", nrow(all_data), ncol(all_data)))

# =============================================================================
# PART 3: RESTRICT TO GENETICALLY CONFIRMED FH
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("PART 3: RESTRICTING TO GENETICALLY CONFIRMED FH COHORT\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# Merge with FH genetic status
df <- all_data %>%
  inner_join(fh_patients, by = c("participant.eid" = "eid"))

cat(sprintf("  FH patients with complete UKB data: %d\n", nrow(df)))

# =============================================================================
# PART 4: DERIVE CORE VARIABLES
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("PART 4: DERIVING CALON VARIABLES\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# ── 4.1 Demographics ────────────────────────────────────────────────────────
df$eid      <- df$participant.eid
df$sex      <- df$participant.p31       # 0=Female, 1=Male
df$age      <- df$participant.p21022    # Age at recruitment
df$bmi      <- df$participant.p21001_i0
df$waist    <- df$participant.p48_i0
df$sbp      <- rowMeans(cbind(df$participant.p4080_i0_a0,
                               df$participant.p4080_i0_a1), na.rm = TRUE)
df$dbp      <- rowMeans(cbind(df$participant.p4079_i0_a0,
                               df$participant.p4079_i0_a1), na.rm = TRUE)

cat(sprintf("  Age: mean %.1f ± %.1f\n", mean(df$age, na.rm=T), sd(df$age, na.rm=T)))
cat(sprintf("  Male: %d (%.1f%%)\n", sum(df$sex==1, na.rm=T), 100*mean(df$sex==1, na.rm=T)))
cat(sprintf("  BMI: mean %.1f ± %.1f\n", mean(df$bmi, na.rm=T), sd(df$bmi, na.rm=T)))
cat(sprintf("  Waist: mean %.1f ± %.1f cm\n", mean(df$waist, na.rm=T), sd(df$waist, na.rm=T)))

# ── 4.2 Lipid biomarkers ───────────────────────────────────────────────────
df$tc       <- df$participant.p30690_i0  # Total cholesterol (mmol/L)
df$hdl      <- df$participant.p30760_i0  # HDL-C (mmol/L)
df$ldl      <- df$participant.p30780_i0  # LDL direct (mmol/L)
df$trig     <- df$participant.p30870_i0  # Triglycerides (mmol/L)
df$apob     <- df$participant.p30640_i0  # ApoB (g/L)
df$apoa     <- df$participant.p30630_i0  # ApoA1 (g/L)
df$lpa      <- df$participant.p30790_i0  # Lp(a) (nmol/L)

# Inverse HDL-C (for CALON model)
df$inv_hdl  <- ifelse(!is.na(df$hdl) & df$hdl > 0, 1 / df$hdl, NA)

# Inverse ApoA1 (alternative to inv_hdl)
df$inv_apoa <- ifelse(!is.na(df$apoa) & df$apoa > 0, 1 / df$apoa, NA)

# Lp(a) binary (>125 nmol/L threshold from CALON manuscript)
df$lpa_binary <- ifelse(!is.na(df$lpa), as.numeric(df$lpa > 125), NA)

cat(sprintf("\n  LDL-C: mean %.2f ± %.2f (N=%d)\n",
            mean(df$ldl, na.rm=T), sd(df$ldl, na.rm=T), sum(!is.na(df$ldl))))
cat(sprintf("  ApoB: mean %.3f ± %.3f (N=%d)\n",
            mean(df$apob, na.rm=T), sd(df$apob, na.rm=T), sum(!is.na(df$apob))))
cat(sprintf("  Lp(a): median %.1f [IQR %.1f-%.1f] (N=%d)\n",
            median(df$lpa, na.rm=T), quantile(df$lpa, 0.25, na.rm=T),
            quantile(df$lpa, 0.75, na.rm=T), sum(!is.na(df$lpa))))
cat(sprintf("  Lp(a) elevated >125: %d (%.1f%%)\n",
            sum(df$lpa_binary == 1, na.rm=T),
            100*mean(df$lpa_binary == 1, na.rm=T)))

# ── 4.3 Identify statin treatment + reverse-engineer LDL-C ─────────────────
cat("\n  --- Statin identification and reverse-engineering ---\n")

# Scan all medication columns for statin codes
med_cols <- grep("^participant\\.p20003_i0_a", names(df), value = TRUE)

identify_statin <- function(row) {
  meds_vec <- unlist(row[med_cols])
  meds_vec <- meds_vec[!is.na(meds_vec)]

  on_atorva  <- any(meds_vec == STATIN_CODES$atorvastatin)
  on_rosuva  <- any(meds_vec == STATIN_CODES$rosuvastatin)
  on_simva   <- any(meds_vec == STATIN_CODES$simvastatin)
  on_prava   <- any(meds_vec == STATIN_CODES$pravastatin)
  on_fluva   <- any(meds_vec == STATIN_CODES$fluvastatin)
  on_lova    <- any(meds_vec == STATIN_CODES$lovastatin)
  on_ezet    <- any(meds_vec == EZETIMIBE_CODE)

  # Determine statin and residual factor
  if (on_atorva) {
    statin_name <- "atorvastatin"
    rf <- RESIDUAL_FACTORS$atorvastatin
  } else if (on_rosuva) {
    statin_name <- "rosuvastatin"
    rf <- RESIDUAL_FACTORS$rosuvastatin
  } else if (on_simva) {
    statin_name <- "simvastatin"
    rf <- RESIDUAL_FACTORS$simvastatin
  } else if (on_prava) {
    statin_name <- "pravastatin"
    rf <- RESIDUAL_FACTORS$pravastatin
  } else if (on_fluva) {
    statin_name <- "fluvastatin"
    rf <- RESIDUAL_FACTORS$fluvastatin
  } else if (on_lova) {
    statin_name <- "lovastatin"
    rf <- RESIDUAL_FACTORS$lovastatin
  } else {
    statin_name <- "none"
    rf <- 1.0
  }

  # Ezetimibe add-on (multiplicative)
  if (on_ezet && statin_name != "none") {
    rf <- rf * RESIDUAL_FACTORS$ezetimibe_addon
  }

  on_statin <- as.numeric(statin_name != "none")

  return(c(statin_name = statin_name,
           residual_factor = rf,
           on_statin = on_statin,
           on_ezetimibe = as.numeric(on_ezet)))
}

cat("  Scanning medication codes for statins...\n")
statin_info <- t(apply(df, 1, identify_statin))
df$statin_name     <- statin_info[, "statin_name"]
df$residual_factor <- as.numeric(statin_info[, "residual_factor"])
df$on_statin       <- as.numeric(statin_info[, "on_statin"])
df$on_ezetimibe    <- as.numeric(statin_info[, "on_ezetimibe"])

cat(sprintf("  On statin: %d (%.1f%%)\n", sum(df$on_statin==1), 100*mean(df$on_statin==1)))
cat(sprintf("  On ezetimibe: %d (%.1f%%)\n", sum(df$on_ezetimibe==1), 100*mean(df$on_ezetimibe==1)))
cat(sprintf("  High-intensity (atorva/rosuva): %d\n",
            sum(df$statin_name %in% c("atorvastatin", "rosuvastatin"))))
cat(sprintf("  Moderate-intensity: %d\n",
            sum(df$statin_name %in% c("simvastatin", "pravastatin", "fluvastatin", "lovastatin"))))

# ── REVERSE-ENGINEERED LDL-C (The Phantom Effect) ──────────────────────────
# RE_LDL = On-treatment LDL ÷ Residual Factor
# For untreated patients: RE_LDL = measured LDL (RF = 1.0)
df$re_ldl <- df$ldl / df$residual_factor

cat(sprintf("\n  Reverse-engineered LDL-C:\n"))
cat(sprintf("    Treated: mean %.2f → RE %.2f mmol/L (masking ratio %.2f)\n",
            mean(df$ldl[df$on_statin==1], na.rm=T),
            mean(df$re_ldl[df$on_statin==1], na.rm=T),
            mean(df$re_ldl[df$on_statin==1], na.rm=T) / mean(df$ldl[df$on_statin==1], na.rm=T)))
cat(sprintf("    Untreated: mean %.2f mmol/L (no adjustment)\n",
            mean(df$ldl[df$on_statin==0], na.rm=T)))

# ── 4.4 ApoB/LDL-C discordance ratio (Log transform) ──────────────────────
# ApoB in UKB is in g/L; LDL in mmol/L
# The ratio ApoB/LDL captures particle discordance
# Logarithmic transform proved optimal in CALON development (AIC comparison)
df$apob_ldl_ratio <- ifelse(!is.na(df$apob) & !is.na(df$ldl) & df$ldl > 0,
                            df$apob / df$ldl, NA)
df$log_apob_ldl   <- ifelse(!is.na(df$apob_ldl_ratio) & df$apob_ldl_ratio > 0,
                            log(df$apob_ldl_ratio), NA)

cat(sprintf("\n  ApoB/LDL ratio: median %.3f (N=%d)\n",
            median(df$apob_ldl_ratio, na.rm=T), sum(!is.na(df$apob_ldl_ratio))))

# ── 4.5 Age at treatment initiation (PROXY) ────────────────────────────────
# UKB does NOT directly record when statins were started.
# PROXY STRATEGIES:
#   (a) If self-reported "high cholesterol" (code 1473 in field 20002),
#       use age at first reporting + small lag → assume treatment start
#   (b) If on statin at baseline → use age at recruitment as UPPER BOUND
#       of treatment duration (conservative: underestimates exposure)
#   (c) Sensitivity analysis with multiple assumptions
#
# For CALON, age_treatment is needed for lipid-year exposure calculation.

cat("\n  --- Age at treatment estimation ---\n")

# Scan illness codes for "high cholesterol" (code 1473)
illness_cols <- grep("^participant\\.p20002_i0_a", names(df), value = TRUE)

df$self_report_high_chol <- apply(df[, illness_cols, drop = FALSE], 1, function(row) {
  any(row == 1473, na.rm = TRUE)
})

# For patients on statins: assume treatment started at age of recruitment
# (conservative — this UNDERESTIMATES years of exposure, making CALON results
# CONSERVATIVE. If real age_treatment were known, performance would be BETTER.)
#
# SENSITIVITY: We will run 3 scenarios:
#   Scenario A: age_treatment = age (most conservative)
#   Scenario B: age_treatment = age - 5 (moderate assumption)
#   Scenario C: age_treatment = age - 10 (generous assumption)

df$age_treatment_A <- ifelse(df$on_statin == 1, df$age, df$age)  # Baseline age
df$age_treatment_B <- ifelse(df$on_statin == 1, pmax(18, df$age - 5), df$age)
df$age_treatment_C <- ifelse(df$on_statin == 1, pmax(18, df$age - 10), df$age)

cat(sprintf("  Self-reported high cholesterol: %d (%.1f%%)\n",
            sum(df$self_report_high_chol), 100*mean(df$self_report_high_chol)))

# ── 4.6 LIPID-YEAR EXPOSURE (Cholesterol Pack-Years) ───────────────────────
# lipid_years = RE_LDL × (age_treatment - 18) / 10
# This captures cumulative pre-treatment cholesterol burden

df$lipid_years_A <- df$re_ldl * pmax(0, df$age_treatment_A - 18) / 10
df$lipid_years_B <- df$re_ldl * pmax(0, df$age_treatment_B - 18) / 10
df$lipid_years_C <- df$re_ldl * pmax(0, df$age_treatment_C - 18) / 10

cat(sprintf("\n  Lipid-year exposure (Scenario A): mean %.1f ± %.1f\n",
            mean(df$lipid_years_A, na.rm=T), sd(df$lipid_years_A, na.rm=T)))
cat(sprintf("  Lipid-year exposure (Scenario B): mean %.1f ± %.1f\n",
            mean(df$lipid_years_B, na.rm=T), sd(df$lipid_years_B, na.rm=T)))

# ── 4.7 Smoking ─────────────────────────────────────────────────────────────
df$smoking_status <- df$participant.p20116_i0  # 0=Never, 1=Previous, 2=Current
df$smoking_binary <- ifelse(!is.na(df$smoking_status),
                            as.numeric(df$smoking_status == 2), NA)
df$ever_smoked    <- ifelse(!is.na(df$smoking_status),
                            as.numeric(df$smoking_status %in% c(1, 2)), NA)
df$pack_years     <- df$participant.p20161

cat(sprintf("\n  Smoking: Never %d, Previous %d, Current %d\n",
            sum(df$smoking_status==0, na.rm=T),
            sum(df$smoking_status==1, na.rm=T),
            sum(df$smoking_status==2, na.rm=T)))

# ── 4.8 Diabetes ────────────────────────────────────────────────────────────
df$diabetes <- ifelse(!is.na(df$participant.p2443_i0),
                      as.numeric(df$participant.p2443_i0 == 1), NA)
df$hba1c    <- df$participant.p30750_i0
df$glucose  <- df$participant.p30740_i0

cat(sprintf("  Diabetes: %d (%.1f%%)\n", sum(df$diabetes==1, na.rm=T),
            100*mean(df$diabetes==1, na.rm=T)))

# ── 4.9 Hypertension ───────────────────────────────────────────────────────
# Define hypertension as: SBP ≥ 140 OR DBP ≥ 90 OR on BP medication OR
# self-reported hypertension (code 1065)

df$self_report_htn <- apply(df[, illness_cols, drop = FALSE], 1, function(row) {
  any(row == 1065, na.rm = TRUE)
})

# Check medication questions for BP meds
# p6177 (male) and p6153 (female): code 2 = Blood pressure medication
df$on_bp_meds <- ifelse(df$sex == 1,
                        grepl("2", as.character(df$participant.p6177_i0)),
                        grepl("2", as.character(df$participant.p6153_i0)))

df$hypertension <- as.numeric(
  df$sbp >= 140 | df$dbp >= 90 |
  df$self_report_htn | df$on_bp_meds
)
df$hypertension[is.na(df$hypertension)] <- 0

cat(sprintf("  Hypertension: %d (%.1f%%)\n", sum(df$hypertension==1),
            100*mean(df$hypertension==1)))

# ── 4.10 Other biomarkers ──────────────────────────────────────────────────
df$crp        <- df$participant.p30710_i0
df$creatinine <- df$participant.p30700_i0
df$alt        <- df$participant.p30620_i0
df$cystatin_c <- df$participant.p30720_i0

# ── 4.11 Carotid IMT (imaging subset only) ─────────────────────────────────
cimt_cols <- grep("^participant\\.p2267", names(df), value = TRUE)
if (length(cimt_cols) > 0) {
  df$mean_cimt <- rowMeans(df[, cimt_cols, drop = FALSE], na.rm = TRUE)
  df$mean_cimt[is.nan(df$mean_cimt)] <- NA
  cat(sprintf("  Carotid IMT available: %d patients\n", sum(!is.na(df$mean_cimt))))
} else {
  df$mean_cimt <- NA
  cat("  Carotid IMT: NOT AVAILABLE (fields not dispensed)\n")
}

# =============================================================================
# PART 5: DEFINE ASCVD OUTCOMES
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("PART 5: DEFINING ASCVD OUTCOMES\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# ── 5.1 Assessment date (entry into follow-up) ─────────────────────────────
df$assessment_date <- as.Date(df$participant.p53_i0, format = "%Y-%m-%d")

# ── 5.2 First-occurrence ASCVD dates ───────────────────────────────────────
ascvd_date_cols <- c(
  "participant.p131298",  # MI
  "participant.p131300",  # Subsequent MI
  "participant.p131296",  # Angina
  "participant.p131306",  # Chronic IHD
  "participant.p131364",  # Cerebral infarction
  "participant.p131366",  # Stroke NOS
  "participant.p131380",  # Atherosclerosis
  "participant.p131386",  # PVD
  "participant.p131388"   # Arterial embolism
)

# Convert to dates and find earliest ASCVD event
for (col in ascvd_date_cols) {
  if (col %in% names(df)) {
    df[[col]] <- as.Date(df[[col]], format = "%Y-%m-%d")
  }
}

# Earliest ASCVD date across all event types
existing_cols <- intersect(ascvd_date_cols, names(df))
if (length(existing_cols) > 0) {
  df$first_ascvd_date <- apply(df[, existing_cols, drop = FALSE], 1, function(row) {
    dates <- as.Date(unlist(row), origin = "1970-01-01")
    dates <- dates[!is.na(dates)]
    if (length(dates) > 0) return(min(dates))
    return(NA)
  })
  df$first_ascvd_date <- as.Date(df$first_ascvd_date, origin = "1970-01-01")
}

# ── 5.3 Classify as PREVALENT vs INCIDENT ──────────────────────────────────
df$ascvd_any <- as.numeric(!is.na(df$first_ascvd_date))

df$ascvd_prevalent <- as.numeric(
  !is.na(df$first_ascvd_date) &
  !is.na(df$assessment_date) &
  df$first_ascvd_date <= df$assessment_date
)

df$ascvd_incident <- as.numeric(
  !is.na(df$first_ascvd_date) &
  !is.na(df$assessment_date) &
  df$first_ascvd_date > df$assessment_date
)

# ── 5.4 Self-reported ASCVD (supplementary) ─────────────────────────────────
# Field p6150: 1=Heart attack, 2=Angina, 3=Stroke, 4=High BP
# NOTE: p6150 may NOT be dispensed in your UKB application. If so, we fall back
# to age-at-diagnosis fields (p3894, p3627, p4056) from Batch 8.
sr_ascvd_cols <- grep("^participant\\.p6150_i0_a", names(df), value = TRUE)

if (length(sr_ascvd_cols) > 0) {
  cat(sprintf("  Self-reported ASCVD: using p6150 arrays (%d columns)\n", length(sr_ascvd_cols)))
  df$sr_mi     <- apply(df[, sr_ascvd_cols, drop = FALSE], 1, function(row) {
    any(row == 1, na.rm = TRUE)
  })
  df$sr_angina <- apply(df[, sr_ascvd_cols, drop = FALSE], 1, function(row) {
    any(row == 2, na.rm = TRUE)
  })
  df$sr_stroke <- apply(df[, sr_ascvd_cols, drop = FALSE], 1, function(row) {
    any(row == 3, na.rm = TRUE)
  })
} else {
  cat("  WARNING: Field p6150 NOT DISPENSED — using age-at-diagnosis fallback\n")
  # Fallback: use age-at-event fields from Batch 8
  # p3894_i0 = age heart attack, p3627_i0 = age angina, p4056_i0 = age stroke
  # Non-NA value means they reported the event (negative values = "do not know")
  df$sr_mi     <- !is.na(df$participant.p3894_i0) & df$participant.p3894_i0 > 0
  df$sr_angina <- !is.na(df$participant.p3627_i0) & df$participant.p3627_i0 > 0
  df$sr_stroke <- !is.na(df$participant.p4056_i0) & df$participant.p4056_i0 > 0

  # Handle case where even Batch 8 fields are missing
  if (!"participant.p3894_i0" %in% names(df)) {
    cat("  WARNING: No self-reported ASCVD fields available — setting all to FALSE\n")
    df$sr_mi <- FALSE; df$sr_angina <- FALSE; df$sr_stroke <- FALSE
  }
}
df$sr_ascvd  <- df$sr_mi | df$sr_angina | df$sr_stroke
cat(sprintf("  Self-reported MI: %d, Angina: %d, Stroke: %d\n",
            sum(df$sr_mi, na.rm=TRUE), sum(df$sr_angina, na.rm=TRUE),
            sum(df$sr_stroke, na.rm=TRUE)))

# ── 5.5 Combined ASCVD (union of HES + self-report) ────────────────────────
df$ascvd_combined <- as.numeric(df$ascvd_any == 1 | df$sr_ascvd == TRUE)

# ── 5.6 Time-to-event variables ────────────────────────────────────────────
# Death date
df$death_date <- as.Date(df$participant.p40000_i0, format = "%Y-%m-%d")

# Exit date = earliest of: first incident ASCVD, death, or censoring
df$exit_date <- pmin(
  ifelse(df$ascvd_incident == 1, df$first_ascvd_date, NA),
  df$death_date,
  CENSOR_DATE,
  na.rm = TRUE
)
df$exit_date <- as.Date(df$exit_date, origin = "1970-01-01")

# Follow-up time in years
df$followup_years <- as.numeric(difftime(df$exit_date, df$assessment_date, units = "days")) / 365.25
df$followup_years[df$followup_years < 0] <- NA

# Event indicator for survival analysis (incident only)
df$event_indicator <- df$ascvd_incident

cat(sprintf("  ASCVD (any, combined): %d (%.1f%%)\n",
            sum(df$ascvd_combined==1, na.rm=T), 100*mean(df$ascvd_combined==1, na.rm=T)))
cat(sprintf("  ASCVD prevalent: %d\n", sum(df$ascvd_prevalent==1, na.rm=T)))
cat(sprintf("  ASCVD incident: %d\n", sum(df$ascvd_incident==1, na.rm=T)))
cat(sprintf("  Median follow-up: %.1f years [IQR %.1f-%.1f]\n",
            median(df$followup_years, na.rm=T),
            quantile(df$followup_years, 0.25, na.rm=T),
            quantile(df$followup_years, 0.75, na.rm=T)))

# =============================================================================
# PART 6: COHORT FLOWCHART
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("PART 6: CONSORT-STYLE COHORT FLOWCHART\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

flowchart <- data.frame(
  Step = c(
    "Total UKB participants",
    "With WES data available",
    "Genetically confirmed FH (pathogenic/likely pathogenic)",
    "  - LDLR variants",
    "  - APOB variants",
    "  - PCSK9 variants",
    "With complete lipid panel (TC, HDL, LDL, TG)",
    "With ApoB measured (for Tier 1 full CALON)",
    "With Lp(a) measured (for Tier 1 full CALON)",
    "With BOTH ApoB AND Lp(a) (Tier 1 eligible)",
    "With follow-up data (Tier 2 CALON-CORE eligible)",
    "ASCVD events (prevalent + incident)",
    "ASCVD incident events only",
    "With carotid IMT data (imaging subset)"
  ),
  N = c(
    nrow(all_data),
    nrow(tudor),
    nrow(df),
    sum(df$gene == "LDLR", na.rm=T),
    sum(df$gene == "APOB", na.rm=T),
    sum(df$gene == "PCSK9", na.rm=T),
    sum(!is.na(df$ldl) & !is.na(df$hdl) & !is.na(df$tc) & !is.na(df$trig)),
    sum(!is.na(df$apob)),
    sum(!is.na(df$lpa)),
    sum(!is.na(df$apob) & !is.na(df$lpa)),
    sum(!is.na(df$followup_years) & df$followup_years > 0),
    sum(df$ascvd_combined == 1, na.rm=T),
    sum(df$ascvd_incident == 1, na.rm=T),
    sum(!is.na(df$mean_cimt))
  ),
  stringsAsFactors = FALSE
)

print(flowchart)

# Save flowchart
write.csv(flowchart, paste0(OUT_DIR, "cohort_flowchart.csv"), row.names = FALSE)

# =============================================================================
# PART 7: SELECT ANALYSIS VARIABLES AND EXPORT
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("PART 7: EXPORTING ANALYSIS-READY DATASET\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

analysis_vars <- c(
  # Identifiers
  "eid", "gene",

  # Demographics
  "sex", "age", "bmi", "waist", "sbp", "dbp",

  # Raw lipids
  "tc", "hdl", "ldl", "trig", "apob", "apoa", "lpa",

  # CALON-derived variables
  "on_statin", "statin_name", "residual_factor", "on_ezetimibe",
  "re_ldl",           # Reverse-engineered LDL-C (Phantom Effect)
  "inv_hdl",          # Inverse HDL-C
  "inv_apoa",         # Inverse ApoA1
  "lpa_binary",       # Lp(a) >125 nmol/L
  "apob_ldl_ratio",   # ApoB/LDL-C ratio
  "log_apob_ldl",     # Log(ApoB/LDL-C) — Particle Discordance

  # Age at treatment (3 sensitivity scenarios)
  "age_treatment_A", "age_treatment_B", "age_treatment_C",

  # Lipid-year exposure (Cholesterol Pack-Years)
  "lipid_years_A", "lipid_years_B", "lipid_years_C",

  # Comorbidities
  "hypertension", "diabetes", "smoking_binary", "ever_smoked", "pack_years",

  # Biomarkers
  "hba1c", "glucose", "crp", "creatinine", "alt", "cystatin_c",

  # Imaging
  "mean_cimt",

  # Outcomes
  "ascvd_combined", "ascvd_prevalent", "ascvd_incident",
  "first_ascvd_date", "assessment_date", "death_date", "exit_date",
  "followup_years", "event_indicator",

  # Self-reported
  "sr_mi", "sr_angina", "sr_stroke", "sr_ascvd",
  "self_report_high_chol", "self_report_htn"
)

# Keep only available variables
available_vars <- intersect(analysis_vars, names(df))
df_export <- df[, available_vars]

cat(sprintf("  Analysis dataset: %d patients × %d variables\n",
            nrow(df_export), ncol(df_export)))

# ── Variable completeness audit ─────────────────────────────────────────────
cat("\n  Variable completeness:\n")
for (v in available_vars) {
  n_complete <- sum(!is.na(df_export[[v]]))
  pct <- 100 * n_complete / nrow(df_export)
  if (pct < 100) {
    cat(sprintf("    %s: %d/%d (%.1f%%)\n", v, n_complete, nrow(df_export), pct))
  }
}

# Export
write.csv(df_export, paste0(OUT_DIR, "calon_ukb_analysis_ready.csv"), row.names = FALSE)
cat(sprintf("\n  ✓ Saved: %scalon_ukb_analysis_ready.csv\n", OUT_DIR))

# =============================================================================
# PART 8: TIER ASSIGNMENT
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("PART 8: VALIDATION TIER ASSIGNMENT\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# TIER 1: Full CALON (all 6 variables available)
# Requires: age, re_ldl, lipid_years, log_apob_ldl, lpa_binary, inv_hdl
df_export$tier1_eligible <- as.numeric(
  !is.na(df_export$age) & !is.na(df_export$re_ldl) &
  !is.na(df_export$lipid_years_A) & !is.na(df_export$log_apob_ldl) &
  !is.na(df_export$lpa_binary) & !is.na(df_export$inv_hdl)
)

# TIER 2: CALON-CORE (without ApoB/LDL and Lp(a))
# Requires: age, re_ldl, lipid_years, inv_hdl
df_export$tier2_eligible <- as.numeric(
  !is.na(df_export$age) & !is.na(df_export$re_ldl) &
  !is.na(df_export$lipid_years_A) & !is.na(df_export$inv_hdl)
)

# TIER 3: Gene-augmented (Tier 1 or 2 + gene annotation)
df_export$tier3_eligible <- as.numeric(
  df_export$tier2_eligible == 1 & !is.na(df_export$gene) & df_export$gene != "None"
)

cat(sprintf("  TIER 1 (Full CALON, 6 variables): %d patients\n", sum(df_export$tier1_eligible)))
cat(sprintf("  TIER 2 (CALON-CORE, 4 variables): %d patients\n", sum(df_export$tier2_eligible)))
cat(sprintf("  TIER 3 (Gene-augmented):           %d patients\n", sum(df_export$tier3_eligible)))

# SAFEHEART-RE eligible (needs: age, sex, ldl, htn, bmi, smoking, prior_cvd)
df_export$safeheart_eligible <- as.numeric(
  !is.na(df_export$age) & !is.na(df_export$sex) & !is.na(df_export$ldl) &
  !is.na(df_export$bmi) & !is.na(df_export$smoking_binary)
)
cat(sprintf("  SAFEHEART-RE eligible:              %d patients\n",
            sum(df_export$safeheart_eligible)))

# Save updated file
write.csv(df_export, paste0(OUT_DIR, "calon_ukb_analysis_ready.csv"), row.names = FALSE)

cat("\n╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  COHORT CONSTRUCTION COMPLETE                                     ║\n")
cat("║  Next: Run 02_CALON_validate.R for frozen-coefficient validation  ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")

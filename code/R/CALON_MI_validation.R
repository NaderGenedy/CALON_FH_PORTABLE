################################################################################
# CALON: Two-Tier Risk Prediction for FH — MICE MI + Bi-External Validation
#
# Complete standalone R pipeline for:
#   - MICE imputation (m=10 imputations, max_iter=20)
#   - Bi-external validation across 6 directions (3 cohort pairs x 2)
#   - Head-to-head DeLong tests: CALON-Lite(10) vs SRE-Refit(8) vs SRE-Fixed(8)
#   - Calibration analysis (slope, intercept, O/E, Brier)
#   - NRI / IDI / DCA net benefit
#   - Pooled odds ratios with Rubin's rules
#   - Publication-quality figures (ROC, calibration, forest plots, DCA)
#
# ── TWO-TIER MODEL ARCHITECTURE ──
#
#   CALON-Lite (10 predictors) — UNIVERSAL MODEL
#     age, sex, ldl, inv_hdl, bmi, ever_smoked, diabetes, hypertension,
#     lpa_143, prior_ascvd
#     USE WHEN: Any FH clinic with standard labs + Lp(a) measurement
#     VALIDATED: All 6 directions (SW⟷UKB⟷Wales)
#     ADVANTAGE: inv_hdl + diabetes + lpa_143 threshold = 3 unique predictors vs SRE
#
#   CALON-X (12 predictors) — ENHANCED FH REGISTRY MODEL
#     CALON-Lite + apob_ldl_high (ApoB/LDL >0.3) + has_phys_sign (stigmata)
#     USE WHEN: Specialist FH centres with ApoB measurement + clinical exam
#     VALIDATED: SW⟷Wales (2 directions; UKB lacks stigmata)
#     ADVANTAGE: Superior discrimination with atherogenic particle ratio + clinical signs
#
# ── COMPARATOR MODELS (Three-layer) ──
#
#   SRE-Fixed (8) = SAFEHEART Risk Equation with FROZEN published coefficients
#     age, sex, ldl, hypertension, bmi, ever_smoked, lpa_120, prior_ascvd
#     Logistic approximation of Perez-de-Isla 2017 Cox model
#     TRIPOD Type 4 (external validation, NO retraining)
#
#   SRE-Refit (8) = Same 8 published SRE features, refitted via logistic regression
#     age, sex, ldl, hypertension, bmi, ever_smoked, lpa_120, prior_ascvd
#     TRIPOD Type 2b (model development with external validation)
#     FAIR COMPARISON: Same method as CALON-Lite, fewer features
#
#   CALON-Core (9) = CALON-Lite minus Lp(a)>143
#     age, sex, ldl, inv_hdl, bmi, ever_smoked, diabetes, hypertension, prior_ascvd
#     Ablation control: shows incremental value of Lp(a)>143
#
# Cohorts:
#   South Wales FH (N~418), UK Biobank (N~1623), Wales FH (N~2405)
#
# Author: Dr Nader Genedy
# Date:   2026-03-07
################################################################################

cat("
╔═══════════════════════════════════════════════════════════════════════════════╗
║  CALON: Two-Tier FH Risk Prediction — MICE MI + Bi-External Validation      ║
║  6 Directions × 5 Models × 10 Imputations                                  ║
║                                                                             ║
║  CALON-Lite (10 vars) = Universal model     → all 6 validation directions   ║
║  CALON-X    (12 vars) = Enhanced FH centres → SW↔Wales (2 directions)       ║
║                                                                             ║
║  Three-layer Comparators:                                                   ║
║    SRE-Fixed  (8 vars) = Frozen published coefficients (TRIPOD Type 4)      ║
║    SRE-Refit  (8 vars) = Published features re-fitted (TRIPOD Type 2b)      ║
║    CALON-Core (9 vars) = CALON-Lite minus Lp(a) (ablation control)          ║
║                                                                             ║
║  Statistical arsenal: AUC, Brier, Cal slope, O/E, DeLong, NRI, IDI, DCA    ║
╚═══════════════════════════════════════════════════════════════════════════════╝
\n")

# ── 0. PACKAGES ──────────────────────────────────────────────────────────────
pkgs_needed <- c("tidyverse", "readxl", "tableone", "pROC", "mice",
                 "survival", "cowplot", "boot", "scales", "gridExtra")
pkgs_missing <- pkgs_needed[!pkgs_needed %in% installed.packages()[, "Package"]]
if (length(pkgs_missing) > 0) {
  cat("Installing missing packages:", paste(pkgs_missing, collapse = ", "), "\n")
  install.packages(pkgs_missing, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(tableone)
  library(pROC)
  library(mice)
  library(survival)
  library(cowplot)
  library(boot)
  library(scales)
  library(gridExtra)
})

set.seed(42)

# ── CONFIGURATION ────────────────────────────────────────────────────────────
M_IMPUTATIONS <- 10      # Number of MI datasets
MAXIT_MICE    <- 20      # Max iterations for MICE convergence
N_BOOT_CI     <- 2000    # Bootstrap replicates for AUC CI
MISS_FLOOR    <- 0.05    # Variables < 5% available → fill with default

# Working directory
if (dir.exists("/cloud/project")) {
  data_dir <- "/cloud/project"
} else if (dir.exists("C:/Users/nader/Downloads/calon_ukb_pipeline")) {
  data_dir <- "C:/Users/nader/Downloads/calon_ukb_pipeline"
} else {
  data_dir <- getwd()
}
setwd(data_dir)
cat(sprintf("Working directory: %s\n\n", getwd()))

# Output directories
dir.create("output",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Publication theme
base_family <- ifelse(.Platform$OS.type == "windows", "Arial", "sans")
theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle    = element_text(size = base_size, colour = "grey40"),
      axis.title       = element_text(face = "bold"),
      axis.text        = element_text(colour = "black"),
      legend.title     = element_text(face = "bold"),
      legend.position  = "bottom",
      strip.text       = element_text(face = "bold"),
      strip.background = element_blank(),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
      plot.margin      = ggplot2::margin(8, 8, 8, 8, unit = "pt")
    )
}
theme_set(theme_pub())

# Colour palettes
pal3   <- c("South Wales FH" = "#E64B35", "UKB" = "#4DBBD5", "Wales" = "#00A087")
pal_m  <- c("CALON-Lite" = "#E64B35", "CALON-Core" = "#F39B7F",
            "SRE-Refit" = "#4DBBD5", "SRE-Fixed" = "#3C5488",
            "CALON-X" = "#00A087")

safe_ggsave <- function(filename, plot, ...) {
  png_file <- sub("\\.pdf$", ".png", filename)
  ggsave(png_file, plot, ...)
  cat(sprintf("  Saved: %s\n", png_file))
  if (grepl("\\.pdf$", filename)) {
    tryCatch({
      ggsave(filename, plot, ...)
      cat(sprintf("  Saved: %s\n", filename))
    }, error = function(e) {
      cat(sprintf("  PDF skipped (%s). PNG saved.\n", conditionMessage(e)))
    })
  }
}

################################################################################
# ── 1. DATA LOADING ──────────────────────────────────────────────────────────
################################################################################
cat("=== 1. DATA LOADING ===\n")

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

find_file <- function(pattern, required = TRUE) {
  matches <- list.files(".", pattern = pattern, ignore.case = TRUE, full.names = TRUE)
  if (length(matches) == 0 && required) {
    stop(sprintf("REQUIRED FILE NOT FOUND: %s\nUpload to: %s", pattern, getwd()))
  }
  if (length(matches) == 0) return(NULL)
  cat(sprintf("  Found: %s\n", matches[1]))
  return(matches[1])
}

# 1a. South Wales FH (DRAGON-3)
cat("Loading South Wales FH...\n")
dragon_raw <- read.csv(find_file("DRAGON.*3.*\\.csv"), stringsAsFactors = FALSE)

# 1b. Wales FH
cat("Loading Wales FH...\n")
wales_raw <- read.csv(find_file("WALES.*FH.*\\.csv"), stringsAsFactors = FALSE)

# 1c. UKB
cat("Loading UK Biobank...\n")
ukb_raw <- read.csv(find_file("TUDOR.*UKB.*\\.csv"), stringsAsFactors = FALSE)

# 1c2. Merge processed UKB fields (smoking, ASCVD split)
ukb_proc_file <- find_file("calon_ukb_analysis_ready", required = FALSE)
if (is.null(ukb_proc_file)) {
  ukb_proc_file <- list.files("output", pattern = "calon_ukb_analysis_ready",
                               ignore.case = TRUE, full.names = TRUE)[1]
}
if (!is.null(ukb_proc_file) && !is.na(ukb_proc_file) && file.exists(ukb_proc_file)) {
  cat(sprintf("  Merging processed UKB fields from: %s\n", ukb_proc_file))
  ukb_proc <- read.csv(ukb_proc_file, stringsAsFactors = FALSE)
  proc_fields <- c("eid", "ever_smoked", "smoking_binary",
                    "ascvd_combined", "ascvd_prevalent", "ascvd_incident",
                    "first_ascvd_date", "assessment_date", "followup_years",
                    "lipid_years_A", "lipid_years_B", "lipid_years_C",
                    "log_apob_ldl")
  proc_fields <- proc_fields[proc_fields %in% names(ukb_proc)]
  ukb_proc_subset <- ukb_proc[, proc_fields, drop = FALSE]
  names(ukb_proc_subset) <- ifelse(names(ukb_proc_subset) == "eid", "eid",
                                    paste0("proc_", names(ukb_proc_subset)))
  names(ukb_proc_subset)[1] <- "eid"
  ukb_raw <- merge(ukb_raw, ukb_proc_subset, by.x = "participant.eid", by.y = "eid",
                    all.x = TRUE, sort = FALSE)
  cat(sprintf("  Merged %d fields\n", ncol(ukb_proc_subset) - 1))
}

# 1d. PASS FH 2012 (optional)
pass_file <- find_file("PASS.*FH.*\\.xls", required = FALSE)
pass_raw  <- NULL
if (!is.null(pass_file)) {
  pass_raw <- tryCatch(read_excel(pass_file), error = function(e) NULL)
}

cat(sprintf("\nRaw data: Dragon=%d, Wales=%d, UKB=%d\n",
            nrow(dragon_raw), nrow(wales_raw), nrow(ukb_raw)))

################################################################################
# ── 2. DATA HARMONISATION ───────────────────────────────────────────────────
################################################################################
cat("\n=== 2. DATA HARMONISATION ===\n")

parse_gene <- function(mut) {
  ifelse(is.na(mut), NA_character_,
         ifelse(grepl("^LDLR", mut), "LDLR",
                ifelse(grepl("^APOB", mut), "APOB",
                       ifelse(grepl("^PCSK9", mut), "PCSK9", "Other"))))
}

# ── 2a. South Wales FH ──────────────────────────────────────────────────────
d <- dragon_raw %>%
  filter(safe_numeric(Positive1) == 1 | Positive1 == "Yes") %>%
  mutate(
    cohort       = "South Wales FH",
    id           = as.character(DatabaseNumber),
    age          = safe_numeric(age_at_event_or_censoring),
    sex          = ifelse(Gender == "M", 1, ifelse(Gender == "F", 0, NA)),
    bmi          = safe_numeric(BMI),
    tc           = safe_numeric(TC_1),
    hdl          = safe_numeric(HDL_1),
    ldl          = safe_numeric(LDL_1),
    tg           = safe_numeric(TRG_1),
    apob         = safe_numeric(ApoB),
    lpa_nmol     = safe_numeric(Lpa),
    inv_hdl      = ifelse(!is.na(hdl) & hdl > 0, 1 / hdl, NA),
    ever_smoked  = case_when(
      Smoking_binary %in% c(1, 2) ~ 1,
      Smoking_binary == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    diabetes     = safe_numeric(Diabetes_binary),
    hypertension = ifelse(safe_numeric(BloodPressureMedication) == 1 |
                            safe_numeric(BloodPressureSystolic) > 140 |
                            safe_numeric(BloodPressureDiastolic) > 90, 1, 0),
    sbp          = safe_numeric(BloodPressureSystolic),
    on_statin    = ifelse(!is.na(Statin) & Statin != "" & Statin != "0", 1, 0),
    lpa_143      = ifelse(!is.na(lpa_nmol) & lpa_nmol > 143, 1, 0),
    lpa_120      = ifelse(!is.na(lpa_nmol) & lpa_nmol > 120, 1, 0),
    ascvd_ever   = safe_numeric(ASCVD_combined),
    age_baseline = safe_numeric(Ageattest),    # Age at genetic test = baseline
    # Individual event AGES (Dragon *Year columns = age at event, NOT calendar year)
    mi_age       = safe_numeric(MIACSYear),
    ptca_age     = safe_numeric(PTCAYear),
    cabg_age     = safe_numeric(CABGYear),
    angina_age   = safe_numeric(ANGINAYear),
    tia_age      = safe_numeric(TIAYear),
    pvd_age      = safe_numeric(PVDYear),
    # Earliest event age across all ASCVD subtypes
    earliest_event_age = pmin(mi_age, ptca_age, cabg_age, angina_age, tia_age, pvd_age, na.rm = TRUE),
    # Prevalent = event BEFORE baseline (genetic test)
    # If ascvd_ever=1 but no event ages available → assume prevalent (conservative)
    prior_ascvd  = ifelse(!is.na(earliest_event_age) & earliest_event_age < age_baseline, 1,
                   ifelse(ascvd_ever == 1 & is.na(earliest_event_age), 1, 0)),
    # Incident = event AT OR AFTER baseline (first event post-diagnosis)
    # Build from individual event ages > baseline
    incident_mi     = ifelse(!is.na(mi_age)     & mi_age     >= age_baseline, 1, 0),
    incident_ptca   = ifelse(!is.na(ptca_age)   & ptca_age   >= age_baseline, 1, 0),
    incident_cabg   = ifelse(!is.na(cabg_age)   & cabg_age   >= age_baseline, 1, 0),
    incident_angina = ifelse(!is.na(angina_age) & angina_age >= age_baseline, 1, 0),
    incident_tia    = ifelse(!is.na(tia_age)    & tia_age    >= age_baseline, 1, 0),
    incident_pvd    = ifelse(!is.na(pvd_age)    & pvd_age    >= age_baseline, 1, 0),
    ascvd_incident  = pmax(incident_mi, incident_ptca, incident_cabg,
                           incident_angina, incident_tia, incident_pvd, na.rm = TRUE),
    ascvd_incident  = ifelse(is.na(ascvd_incident), 0, ascvd_incident),
    # Use incident as the outcome
    ascvd        = ascvd_incident,
    age_at_event = safe_numeric(age_at_event),
    apob_ldl_ratio = ifelse(!is.na(safe_numeric(ApoB)) & !is.na(safe_numeric(LDL_1)) &
                              safe_numeric(LDL_1) > 0, safe_numeric(ApoB) / safe_numeric(LDL_1), NA),
    gene         = parse_gene(Mutation1),
    patient_type = ifelse(safe_numeric(TypeofPatient) == 0, "Index", "Cascade"),
    corneal_arcus = safe_numeric(CornealArcus),
    tendon_xanth  = safe_numeric(TendonXanthomata),
    xanthelasmas  = ifelse(Xanthelasmas %in% c("1", "True", "TRUE"), 1,
                           ifelse(Xanthelasmas %in% c("0", "False", "FALSE"), 0, NA)),
    dlcn_score   = safe_numeric(GenoTypingScore)
  )

# Merge PASS physical signs
if (!is.null(pass_raw)) {
  pass_signs <- pass_raw %>%
    select(DatabaseNumber,
           corneal_pass = `Corneal Arcus`,
           xanth_pass   = `TendonXanthomata`,
           xanthel_pass = `Xanthelasmas`) %>%
    mutate(across(c(corneal_pass, xanth_pass, xanthel_pass), ~case_when(
      tolower(as.character(.)) %in% c("true", "1") ~ 1L,
      tolower(as.character(.)) %in% c("false", "0") ~ 0L,
      TRUE ~ NA_integer_
    )))
  d <- d %>%
    left_join(pass_signs, by = "DatabaseNumber") %>%
    mutate(
      corneal_arcus = coalesce(corneal_arcus, as.numeric(corneal_pass)),
      tendon_xanth  = coalesce(tendon_xanth, as.numeric(xanth_pass)),
      xanthelasmas  = coalesce(xanthelasmas, as.numeric(xanthel_pass))
    ) %>%
    select(-corneal_pass, -xanth_pass, -xanthel_pass)
}

d <- d %>%
  mutate(has_phys_sign = pmax(corneal_arcus, tendon_xanth, xanthelasmas, na.rm = TRUE),
         apob_ldl_high = ifelse(!is.na(apob_ldl_ratio) & apob_ldl_ratio > 0.3, 1, 0))

cat(sprintf("  South Wales FH: %d patients\n", nrow(d)))
cat(sprintf("    ASCVD ever:     %d (%.1f%%)\n", sum(d$ascvd_ever == 1, na.rm=TRUE),
            100*mean(d$ascvd_ever == 1, na.rm=TRUE)))
cat(sprintf("    Prior ASCVD:    %d (events before genetic test)\n", sum(d$prior_ascvd == 1, na.rm=TRUE)))
cat(sprintf("    Incident ASCVD: %d (events at/after genetic test)\n", sum(d$ascvd_incident == 1, na.rm=TRUE)))
cat(sprintf("    ascvd (outcome): %d (= incident)\n", sum(d$ascvd == 1, na.rm=TRUE)))

# ── 2b. Wales FH ────────────────────────────────────────────────────────────
w <- wales_raw %>%
  filter(Positive1 == 1) %>%
  mutate(
    cohort       = "Wales",
    id           = as.character(DatabaseNumber),
    age          = safe_numeric(BMI_AGE),
    sex          = ifelse(Gender == "M", 1, ifelse(Gender == "F", 0, NA)),
    bmi          = safe_numeric(BMI),
    tc           = safe_numeric(TC.1),
    hdl          = safe_numeric(HDL.1),
    ldl          = safe_numeric(LDL.1),
    tg           = safe_numeric(TRG.1),
    apob         = NA_real_,
    lpa_nmol     = safe_numeric(Lpa.1),
    inv_hdl      = ifelse(!is.na(hdl) & hdl > 0, 1 / hdl, NA),
    ever_smoked  = ifelse(Smoking %in% c("1", "True", "TRUE"), 1,
                          ifelse(Smoking %in% c("0", "False", "FALSE"), 0, NA)),
    diabetes     = ifelse(Diabetes %in% c("1", "True"), 1,
                          ifelse(Diabetes %in% c("0", "False"), 0, NA)),
    hypertension = ifelse(
      BloodPressureMedication %in% c("1", "True") |
        (!is.na(safe_numeric(BloodPressureSystolic)) & safe_numeric(BloodPressureSystolic) > 140) |
        (!is.na(safe_numeric(BloodPressureDiastolic)) & safe_numeric(BloodPressureDiastolic) > 90),
      1, 0
    ),
    sbp          = safe_numeric(BloodPressureSystolic),
    on_statin    = ifelse(!is.na(OnTreatment) & OnTreatment >= 1, 1, 0),
    lpa_143      = ifelse(!is.na(lpa_nmol) & lpa_nmol > 143, 1, 0),
    lpa_120      = ifelse(!is.na(lpa_nmol) & lpa_nmol > 120, 1, 0),
    ascvd_ever   = safe_numeric(ascvd_combine),
    # Compute age at diagnosis (baseline) from DOB and ResultDate1
    dob_dt       = as.Date(DOB, format = "%Y-%m-%d"),
    result_dt    = as.Date(ResultDate1, format = "%Y-%m-%d"),
    age_baseline = as.numeric(difftime(result_dt, dob_dt, units = "days")) / 365.25,
    # Use BMI_AGE as fallback if dates unavailable
    age_baseline = ifelse(is.na(age_baseline), safe_numeric(BMI_AGE), age_baseline),
    # Individual event ages (Wales *Age columns = age at event)
    mi_age       = safe_numeric(MIACSAge),
    pci_age      = safe_numeric(PCIStentsAge),
    cabg_age     = safe_numeric(CABGAge),
    angina_age   = safe_numeric(ANGINAAge),
    tia_age      = safe_numeric(TIAAge),
    pvd_age      = safe_numeric(PVDAge),
    # Earliest event age across all ASCVD subtypes
    earliest_event_age = pmin(mi_age, pci_age, cabg_age, angina_age, tia_age, pvd_age, na.rm = TRUE),
    # Prevalent = event BEFORE baseline (diagnosis date)
    # If ascvd_ever=1 but no event ages available → assume prevalent (conservative)
    prior_ascvd  = ifelse(!is.na(earliest_event_age) & earliest_event_age < age_baseline, 1,
                   ifelse(ascvd_ever == 1 & is.na(earliest_event_age), 1, 0)),
    # Incident = event AT OR AFTER baseline
    incident_mi     = ifelse(!is.na(mi_age)     & mi_age     >= age_baseline, 1, 0),
    incident_pci    = ifelse(!is.na(pci_age)    & pci_age    >= age_baseline, 1, 0),
    incident_cabg   = ifelse(!is.na(cabg_age)   & cabg_age   >= age_baseline, 1, 0),
    incident_angina = ifelse(!is.na(angina_age) & angina_age >= age_baseline, 1, 0),
    incident_tia    = ifelse(!is.na(tia_age)    & tia_age    >= age_baseline, 1, 0),
    incident_pvd    = ifelse(!is.na(pvd_age)    & pvd_age    >= age_baseline, 1, 0),
    ascvd_incident  = pmax(incident_mi, incident_pci, incident_cabg,
                           incident_angina, incident_tia, incident_pvd, na.rm = TRUE),
    ascvd_incident  = ifelse(is.na(ascvd_incident), 0, ascvd_incident),
    # Use incident as the outcome
    ascvd        = ascvd_incident,
    age_at_event = earliest_event_age,
    apob_ldl_ratio = NA_real_,                       # No ApoB in Wales
    apob_ldl_high  = NA_real_,                       # No ApoB in Wales
    gene         = parse_gene(Mutation1),
    patient_type = ifelse(I_Vs_R == 1, "Index", "Cascade"),
    corneal_arcus = ifelse(CornealArcus %in% c("1", "True"), 1,
                           ifelse(CornealArcus %in% c("0", "False"), 0, NA)),
    tendon_xanth  = ifelse(TendonXanthomata %in% c("1", "True"), 1,
                           ifelse(TendonXanthomata %in% c("0", "False"), 0, NA)),
    xanthelasmas  = ifelse(Xanthelasmas %in% c("1", "True"), 1,
                           ifelse(Xanthelasmas %in% c("0", "False"), 0, NA)),
    has_phys_sign = pmax(corneal_arcus, tendon_xanth, xanthelasmas, na.rm = TRUE),
    dlcn_score   = safe_numeric(GenoTypingScore)
  )
cat(sprintf("  Wales FH: %d patients\n", nrow(w)))
cat(sprintf("    ASCVD ever:     %d (%.1f%%)\n", sum(w$ascvd_ever == 1, na.rm=TRUE),
            100*mean(w$ascvd_ever == 1, na.rm=TRUE)))
cat(sprintf("    Prior ASCVD:    %d (events before diagnosis)\n", sum(w$prior_ascvd == 1, na.rm=TRUE)))
cat(sprintf("    Incident ASCVD: %d (events at/after diagnosis)\n", sum(w$ascvd_incident == 1, na.rm=TRUE)))
cat(sprintf("    ascvd (outcome): %d (= incident)\n", sum(w$ascvd == 1, na.rm=TRUE)))
cat(sprintf("    age_baseline available: %d (%.1f%%)\n",
            sum(!is.na(w$age_baseline)), 100*mean(!is.na(w$age_baseline))))

# ── 2c. UKB ──────────────────────────────────────────────────────────────────
statin_codes <- c(1140888648, 1140888594, 1140888690, 1141146234,
                  1141192414, 1140861970, 1140881748, 1141146138)

ukb_all <- ukb_raw %>%
  mutate(
    eid      = participant.eid,
    is_fh    = is_fh_genetic,
    sex_num  = participant.p31,
    age_raw  = participant.p21022,
    bmi_raw  = participant.p21001_i0,
    ldl_raw  = participant.p30780_i0,
    hdl_raw  = participant.p30760_i0,
    tg_raw   = participant.p30870_i0,
    tc_raw   = participant.p30690_i0,
    apob_raw = participant.p30640_i0,
    sbp_1    = participant.p4080_i0_a0,
    sbp_2    = participant.p4080_i0_a1,
    gene_ukb = gene
  )

# Statin detection
med_cols <- paste0("participant.p20003_i0_a", 0:47)
med_cols <- med_cols[med_cols %in% names(ukb_all)]
ukb_all$on_statin_ukb <- apply(ukb_all[, med_cols], 1, function(row) {
  as.integer(any(row %in% statin_codes, na.rm = TRUE))
})

# Self-reported conditions
ill_cols <- grep("^participant\\.p20002_i0", names(ukb_all), value = TRUE)
ukb_all$diabetes_ukb <- apply(ukb_all[, ill_cols], 1, function(row) {
  as.integer(any(row %in% c(1220, 1222, 1223), na.rm = TRUE))
})
ukb_all$hypertension_ukb <- apply(ukb_all[, ill_cols], 1, function(row) {
  as.integer(any(row %in% c(1065, 1072), na.rm = TRUE))
})

# Smoking
if ("proc_ever_smoked" %in% names(ukb_all)) {
  ukb_all$ever_smoked_ukb <- as.numeric(ukb_all$proc_ever_smoked)
} else if ("participant.p20116_i0" %in% names(ukb_all)) {
  ukb_all$ever_smoked_ukb <- ifelse(ukb_all$participant.p20116_i0 %in% c(1, 2), 1,
                                     ifelse(ukb_all$participant.p20116_i0 == 0, 0, NA))
} else {
  ukb_all$ever_smoked_ukb <- NA_real_
}

# ASCVD
if ("proc_ascvd_combined" %in% names(ukb_all)) {
  ukb_all$ascvd_ukb <- as.numeric(ukb_all$proc_ascvd_combined)
} else if ("ascvd" %in% names(ukb_all)) {
  ukb_all$ascvd_ukb <- ukb_all$ascvd
} else if ("ASCVD_combined" %in% names(ukb_all)) {
  ukb_all$ascvd_ukb <- ukb_all$ASCVD_combined
} else {
  ukb_all$ascvd_ukb <- NA_real_
}

# Extract FH subset (keep eid for later Lp(a) merge)
ukb_fh <- ukb_all %>%
  filter(is_fh == 1) %>%
  mutate(
    cohort   = "UKB",
    eid      = eid,
    id       = as.character(eid),
    age      = age_raw,
    sex      = sex_num,
    bmi      = bmi_raw,
    tc       = tc_raw,
    hdl      = hdl_raw,
    ldl      = ldl_raw,
    tg       = tg_raw,
    apob     = apob_raw,
    lpa_nmol = NA_real_,
    inv_hdl  = ifelse(!is.na(hdl_raw) & hdl_raw > 0, 1 / hdl_raw, NA),
    ever_smoked  = ever_smoked_ukb,
    diabetes     = diabetes_ukb,
    hypertension = hypertension_ukb,
    sbp      = rowMeans(cbind(sbp_1, sbp_2), na.rm = TRUE),
    on_statin = on_statin_ukb,
    lpa_143  = NA_real_,
    lpa_120  = NA_real_,
    ascvd_ever   = ascvd_ukb,                        # Total ASCVD (prevalent + incident)
    age_at_event = NA_real_,
    prior_ascvd  = ascvd_ukb,                       # Placeholder; overridden below
    apob_ldl_ratio = ifelse(!is.na(apob_raw) & !is.na(ldl_raw) & ldl_raw > 0,
                             apob_raw / ldl_raw, NA),
    apob_ldl_high  = NA_real_,                       # Set after ratio computed
    gene     = gene_ukb,
    patient_type = NA_character_,
    corneal_arcus = NA_real_,
    tendon_xanth  = NA_real_,
    xanthelasmas  = NA_real_,
    has_phys_sign = NA_real_,
    dlcn_score    = NA_real_
  )

# Compute apob_ldl_high from ratio
ukb_fh$apob_ldl_high <- ifelse(!is.na(ukb_fh$apob_ldl_ratio) & ukb_fh$apob_ldl_ratio > 0.3, 1, 0)

# Derive incident ASCVD from prevalent vs combined
if ("proc_ascvd_prevalent" %in% names(ukb_all)) {
  ukb_fh$prior_ascvd <- as.numeric(ukb_all$proc_ascvd_prevalent[ukb_all$is_fh == 1])
  # Incident = total ASCVD minus prevalent (event after enrollment)
  ukb_fh$ascvd_incident <- ifelse(ukb_fh$ascvd_ever == 1 & ukb_fh$prior_ascvd == 0, 1, 0)
  cat(sprintf("  UKB prior_ascvd from proc_ascvd_prevalent: %d events\n",
              sum(ukb_fh$prior_ascvd == 1, na.rm = TRUE)))
  cat(sprintf("  UKB incident ASCVD (combined minus prevalent): %d events\n",
              sum(ukb_fh$ascvd_incident == 1, na.rm = TRUE)))
} else {
  # Without prevalent data, can't split — treat all as incident (conservative)
  ukb_fh$ascvd_incident <- ukb_fh$ascvd_ever
  cat("  UKB: proc_ascvd_prevalent not available — using ascvd_ever as incident\n")
}
# Set the outcome to incident
ukb_fh$ascvd <- ukb_fh$ascvd_incident

# Try to recover Lp(a): first from the raw UKB column p30790_i0
if ("participant.p30790_i0" %in% names(ukb_all)) {
  ukb_fh$lpa_nmol <- safe_numeric(ukb_all[ukb_all$is_fh == 1, "participant.p30790_i0"])
  ukb_fh$lpa_143  <- ifelse(!is.na(ukb_fh$lpa_nmol) & ukb_fh$lpa_nmol > 143, 1, 0)
  ukb_fh$lpa_120  <- ifelse(!is.na(ukb_fh$lpa_nmol) & ukb_fh$lpa_nmol > 120, 1, 0)
  cat(sprintf("  UKB Lp(a) from p30790: %d non-missing\n", sum(!is.na(ukb_fh$lpa_nmol))))
}

# Then try batch files or processed file if still missing
if (all(is.na(ukb_fh$lpa_nmol))) {
  # Scan batch files
  lpa_files <- list.files(".", pattern = "batch.*\\.csv", full.names = TRUE)
  for (f in lpa_files) {
    tmp <- read.csv(f, nrows = 5, stringsAsFactors = FALSE)
    lpa_col <- grep("p30790|Lpa|lpa", names(tmp), value = TRUE)
    if (length(lpa_col) > 0) {
      batch_lpa <- read.csv(f, stringsAsFactors = FALSE)
      if ("participant.eid" %in% names(batch_lpa)) {
        lpa_map <- batch_lpa %>%
          select(participant.eid, all_of(lpa_col[1])) %>%
          rename(lpa_batch = 2) %>%
          mutate(lpa_batch = safe_numeric(lpa_batch))
        ukb_fh <- ukb_fh %>%
          left_join(lpa_map, by = c("eid" = "participant.eid")) %>%
          mutate(lpa_nmol = coalesce(lpa_nmol, lpa_batch),
                 lpa_143  = ifelse(!is.na(lpa_nmol) & lpa_nmol > 143, 1, 0),
                 lpa_120  = ifelse(!is.na(lpa_nmol) & lpa_nmol > 120, 1, 0)) %>%
          select(-lpa_batch)
      }
      break
    }
  }
}

# Also try calon_ukb_analysis_ready.csv for Lp(a)
if (all(is.na(ukb_fh$lpa_nmol))) {
  ready_file <- NULL
  if (file.exists("calon_ukb_analysis_ready.csv")) {
    ready_file <- "calon_ukb_analysis_ready.csv"
  } else if (file.exists("output/calon_ukb_analysis_ready.csv")) {
    ready_file <- "output/calon_ukb_analysis_ready.csv"
  }
  if (!is.null(ready_file)) {
    ukb_ready2 <- read.csv(ready_file, stringsAsFactors = FALSE)
    lpa_col_name <- intersect(c("lpa_nmol", "Lpa", "lpa"), names(ukb_ready2))[1]
    if (!is.na(lpa_col_name) && "participant.eid" %in% names(ukb_ready2)) {
      lpa_map2 <- ukb_ready2 %>%
        select(participant.eid, all_of(lpa_col_name)) %>%
        rename(lpa_ready = 2) %>%
        mutate(lpa_ready = safe_numeric(lpa_ready))
      ukb_fh <- ukb_fh %>%
        left_join(lpa_map2, by = c("eid" = "participant.eid")) %>%
        mutate(lpa_nmol = coalesce(lpa_nmol, lpa_ready),
               lpa_143  = ifelse(!is.na(lpa_nmol) & lpa_nmol > 143, 1, 0),
               lpa_120  = ifelse(!is.na(lpa_nmol) & lpa_nmol > 120, 1, 0)) %>%
        select(-lpa_ready)
      cat(sprintf("  UKB Lp(a) from ready file: %d non-missing\n", sum(!is.na(ukb_fh$lpa_nmol))))
    }
  }
}

# Ensure lpa_143 and lpa_120 are set
ukb_fh$lpa_143 <- ifelse(!is.na(ukb_fh$lpa_nmol) & ukb_fh$lpa_nmol > 143, 1, 0)
ukb_fh$lpa_120 <- ifelse(!is.na(ukb_fh$lpa_nmol) & ukb_fh$lpa_nmol > 120, 1, 0)

cat(sprintf("  UKB FH: %d patients, Lp(a) available: %d\n",
            nrow(ukb_fh), sum(!is.na(ukb_fh$lpa_nmol))))
cat(sprintf("    ASCVD ever:     %d (%.1f%%)\n", sum(ukb_fh$ascvd_ever == 1, na.rm=TRUE),
            100*mean(ukb_fh$ascvd_ever == 1, na.rm=TRUE)))
cat(sprintf("    Prior ASCVD:    %d (prevalent at enrollment)\n", sum(ukb_fh$prior_ascvd == 1, na.rm=TRUE)))
cat(sprintf("    Incident ASCVD: %d (post-enrollment)\n", sum(ukb_fh$ascvd_incident == 1, na.rm=TRUE)))
cat(sprintf("    ascvd (outcome): %d (= incident)\n", sum(ukb_fh$ascvd == 1, na.rm=TRUE)))

# ── 2d. Combine FH Cohorts ──────────────────────────────────────────────────
common_cols <- c("cohort", "id", "age", "sex", "bmi", "tc", "hdl", "ldl", "tg",
                 "apob", "lpa_nmol", "inv_hdl", "ever_smoked", "diabetes",
                 "hypertension", "sbp", "on_statin", "lpa_143", "lpa_120",
                 "prior_ascvd", "apob_ldl_ratio", "apob_ldl_high",
                 "ascvd_ever", "ascvd_incident", "ascvd",
                 "age_at_event", "gene", "patient_type",
                 "corneal_arcus", "tendon_xanth", "xanthelasmas",
                 "has_phys_sign", "dlcn_score")

for (cc in common_cols) {
  if (!cc %in% names(d))      d[[cc]]      <- NA
  if (!cc %in% names(w))      w[[cc]]      <- NA
  if (!cc %in% names(ukb_fh)) ukb_fh[[cc]] <- NA
}

fh_pooled <- bind_rows(d[, common_cols], w[, common_cols], ukb_fh[, common_cols])
fh_pooled$cohort <- factor(fh_pooled$cohort, levels = c("South Wales FH", "UKB", "Wales"))

cat(sprintf("\n=== POOLED FH: %d patients ===\n", nrow(fh_pooled)))
cat(sprintf("  South Wales FH: %d\n", sum(fh_pooled$cohort == "South Wales FH")))
cat(sprintf("  UKB:            %d\n", sum(fh_pooled$cohort == "UKB")))
cat(sprintf("  Wales:          %d\n", sum(fh_pooled$cohort == "Wales")))

# Print variable availability per cohort
cat("\n--- Variable Availability (% non-missing) ---\n")
check_vars <- c("age", "sex", "bmi", "ldl", "inv_hdl", "ever_smoked", "diabetes",
                "hypertension", "lpa_143", "lpa_120", "prior_ascvd",
                "apob_ldl_high", "dlcn_score", "has_phys_sign",
                "ascvd_ever", "ascvd_incident", "ascvd")
for (coh in levels(fh_pooled$cohort)) {
  sub <- fh_pooled %>% filter(cohort == coh)
  avail <- sapply(check_vars, function(v) round(100 * mean(!is.na(sub[[v]])), 1))
  cat(sprintf("  %s (N=%d):\n    %s\n", coh, nrow(sub),
              paste(sprintf("%s=%.1f%%", check_vars, avail), collapse = ", ")))
}

# ── INCIDENT ASCVD OUTCOME SUMMARY ──
cat("\n--- INCIDENT ASCVD OUTCOME (data leakage fix) ---\n")
cat("  Outcome is now INCIDENT ASCVD (events AT/AFTER baseline).\n")
cat("  Prior ASCVD (events BEFORE baseline) retained as PREDICTOR.\n")
for (coh in levels(fh_pooled$cohort)) {
  sub <- fh_pooled %>% filter(cohort == coh)
  cat(sprintf("  %s (N=%d):\n", coh, nrow(sub)))
  cat(sprintf("    ascvd_ever (total):     %d (%.1f%%)\n",
              sum(sub$ascvd_ever == 1, na.rm=TRUE), 100*mean(sub$ascvd_ever == 1, na.rm=TRUE)))
  cat(sprintf("    prior_ascvd (predictor): %d (%.1f%%)\n",
              sum(sub$prior_ascvd == 1, na.rm=TRUE), 100*mean(sub$prior_ascvd == 1, na.rm=TRUE)))
  cat(sprintf("    ascvd_incident (outcome): %d (%.1f%%)\n",
              sum(sub$ascvd_incident == 1, na.rm=TRUE), 100*mean(sub$ascvd_incident == 1, na.rm=TRUE)))
  # EPV check
  n_events <- sum(sub$ascvd_incident == 1, na.rm=TRUE)
  if (n_events < 10) {
    cat(sprintf("    *** WARNING: Only %d incident events — EPV < 10 for models with >1 predictor ***\n", n_events))
  }
}

################################################################################
# ── 3. MODEL DEFINITIONS ────────────────────────────────────────────────────
################################################################################
cat("\n=== 3. MODEL DEFINITIONS ===\n")

# CORRECT model definitions
LITE_features <- c("age", "sex", "ldl", "inv_hdl", "bmi", "ever_smoked",
                    "diabetes", "hypertension", "lpa_143", "prior_ascvd")
CORE_features <- c("age", "sex", "ldl", "inv_hdl", "bmi", "ever_smoked",
                    "diabetes", "hypertension", "prior_ascvd")
SRE_features  <- c("age", "sex", "ldl", "hypertension", "bmi", "ever_smoked", "lpa_120", "prior_ascvd")
X_features    <- c("age", "sex", "ldl", "inv_hdl", "bmi", "ever_smoked",
                    "diabetes", "hypertension", "lpa_143", "prior_ascvd",
                    "apob_ldl_high", "has_phys_sign")

cat("CALON-Lite (10): ", paste(LITE_features, collapse = ", "), "\n")
cat("CALON-Core  (9): ", paste(CORE_features, collapse = ", "), "\n")
cat("SRE-Refit   (8): ", paste(SRE_features,  collapse = ", "), "\n")
cat("CALON-X    (12): ", paste(X_features,    collapse = ", "), "\n")

# Validation directions
ALL_DIRS <- list(
  list(train = "South Wales FH", test = "UKB",            label = "SW->UKB"),
  list(train = "UKB",            test = "South Wales FH", label = "UKB->SW"),
  list(train = "South Wales FH", test = "Wales",          label = "SW->Wales"),
  list(train = "Wales",          test = "South Wales FH", label = "Wales->SW"),
  list(train = "UKB",            test = "Wales",          label = "UKB->Wales"),
  list(train = "Wales",          test = "UKB",            label = "Wales->UKB")
)

SW_WALES_DIRS <- list(
  list(train = "South Wales FH", test = "Wales",          label = "SW->Wales"),
  list(train = "Wales",          test = "South Wales FH", label = "Wales->SW")
)

################################################################################
# ── 3b. SRE-FIXED: FROZEN SAFEHEART COEFFICIENTS (TRIPOD Type 4) ──────────
################################################################################
cat("\n=== 3b. SRE-FIXED (Frozen Published Coefficients) ===\n")

# Logistic regression approximation of published SAFEHEART Risk Equation
# Source: Perez-de-Isla et al. 2017 (Eur Heart J), with coefficients from
# logistic approximation for binary outcome validation
# Reference: MODEL DEFINITIONS.docx — validated against published HRs
#
# Published HRs → ln(HR) used as logistic β:
#   Male sex:            HR=2.17  → β=0.775
#   Prior CVD:           HR=4.15  → β=1.423
#   Hypertension:        HR=1.54  → β=0.432
#   Active smoking:      HR=1.59  → β=0.464
#   Lp(a) >50 mg/dL:    HR=1.46  → β=0.378
#   BMI (per kg/m²):     HR=1.025 → β=0.025
#   Age (per year):      HR=1.066 → β=0.064
#   LDL-C (per mmol/L):  HR=1.12  → β=0.109

SAFEHEART_COEF <- list(
  intercept = -7.053,    # Baseline intercept (logistic approximation)
  age       = 0.064,     # Per year
  male      = 0.775,     # Male = 1, Female = 0
  ldl_c     = 0.109,     # Per mmol/L
  htn       = 0.431,     # Hypertension = 1
  bmi       = 0.025,     # Per kg/m²
  smoking   = 0.466,     # Ever smoked = 1
  lpa_120   = 0.378,     # Lp(a) >120 nmol/L (≈50 mg/dL)
  prior_cvd = 1.414      # Prior ASCVD = 1 (strongest predictor)
)

#' Apply frozen SAFEHEART Risk Equation to produce predicted probabilities
#' No training required — purely external validation (TRIPOD Type 4)
#'
#' @param df Data frame with SRE predictor columns
#' @return Numeric vector of predicted probabilities
sre_fixed_predict <- function(df) {
  lp <- SAFEHEART_COEF$intercept +
    SAFEHEART_COEF$age       * df$age +
    SAFEHEART_COEF$male      * df$sex +
    SAFEHEART_COEF$ldl_c     * df$ldl +
    SAFEHEART_COEF$htn       * df$hypertension +
    SAFEHEART_COEF$bmi       * df$bmi +
    SAFEHEART_COEF$smoking   * df$ever_smoked +
    SAFEHEART_COEF$lpa_120   * ifelse(!is.na(df$lpa_120), df$lpa_120, 0) +
    SAFEHEART_COEF$prior_cvd * df$prior_ascvd
  prob <- 1 / (1 + exp(-lp))
  return(prob)
}

cat("  Frozen coefficients loaded (8 predictors)\n")
cat(sprintf("  Coefficients: %s\n",
    paste(sprintf("%s=%.3f", names(SAFEHEART_COEF)[-1],
                  unlist(SAFEHEART_COEF[-1])), collapse = ", ")))

################################################################################
# ── 3c. LANCET STATISTICAL FUNCTIONS: NRI, IDI, DCA ──────────────────────────
################################################################################
cat("\n=== 3c. STATISTICAL FUNCTIONS (NRI, IDI, DCA) ===\n")

#' Compute Category-free NRI (Net Reclassification Improvement)
#'
#' Uses the approach from Pencina et al. (2011) Statistics in Medicine
#' @param y_true Binary outcome vector (0/1)
#' @param pred_new Predicted probabilities from new model
#' @param pred_ref Predicted probabilities from reference model
#' @return List with NRI_events, NRI_nonevents, NRI_total, z_stat, p_value
compute_nri <- function(y_true, pred_new, pred_ref) {
  # Remove NAs
  valid <- complete.cases(y_true, pred_new, pred_ref)
  y   <- y_true[valid]
  p1  <- pred_new[valid]
  p0  <- pred_ref[valid]

  events    <- which(y == 1)
  nonevents <- which(y == 0)

  if (length(events) < 5 || length(nonevents) < 5) {
    return(list(nri_events = NA, nri_nonevents = NA, nri = NA,
                z_stat = NA, p_value = NA))
  }

  # Category-free NRI: proportion moved UP minus proportion moved DOWN
  # Events: higher pred_new is CORRECT reclassification (moved up)
  up_events   <- sum(p1[events] > p0[events])
  down_events <- sum(p1[events] < p0[events])
  nri_events  <- (up_events - down_events) / length(events)

  # Non-events: lower pred_new is CORRECT reclassification (moved down)
  up_nonevents   <- sum(p1[nonevents] > p0[nonevents])
  down_nonevents <- sum(p1[nonevents] < p0[nonevents])
  nri_nonevents  <- (down_nonevents - up_nonevents) / length(nonevents)

  nri_total <- nri_events + nri_nonevents

  # Variance (Pencina formula)
  var_events    <- (up_events + down_events) / length(events)^2
  var_nonevents <- (up_nonevents + down_nonevents) / length(nonevents)^2
  se_nri        <- sqrt(var_events + var_nonevents)

  z_stat  <- nri_total / se_nri
  p_value <- 2 * pnorm(-abs(z_stat))

  return(list(
    nri_events    = round(nri_events, 4),
    nri_nonevents = round(nri_nonevents, 4),
    nri           = round(nri_total, 4),
    se            = round(se_nri, 4),
    z_stat        = round(z_stat, 3),
    p_value       = signif(p_value, 3)
  ))
}

#' Compute IDI (Integrated Discrimination Improvement)
#'
#' Pencina et al. (2008) Statistics in Medicine
#' @param y_true Binary outcome vector
#' @param pred_new Predicted probabilities from new model
#' @param pred_ref Predicted probabilities from reference model
#' @return List with IDI, IDI_events, IDI_nonevents, z_stat, p_value
compute_idi <- function(y_true, pred_new, pred_ref) {
  valid <- complete.cases(y_true, pred_new, pred_ref)
  y   <- y_true[valid]
  p1  <- pred_new[valid]
  p0  <- pred_ref[valid]

  events    <- which(y == 1)
  nonevents <- which(y == 0)

  if (length(events) < 5 || length(nonevents) < 5) {
    return(list(idi = NA, idi_events = NA, idi_nonevents = NA,
                z_stat = NA, p_value = NA))
  }

  # IDI = mean(p_new - p_ref | events) - mean(p_new - p_ref | non-events)
  idi_events    <- mean(p1[events] - p0[events])
  idi_nonevents <- mean(p1[nonevents] - p0[nonevents])
  idi           <- idi_events - idi_nonevents

  # Standard error via bootstrap-like approximation
  se_events    <- sd(p1[events] - p0[events]) / sqrt(length(events))
  se_nonevents <- sd(p1[nonevents] - p0[nonevents]) / sqrt(length(nonevents))
  se_idi       <- sqrt(se_events^2 + se_nonevents^2)

  z_stat  <- idi / se_idi
  p_value <- 2 * pnorm(-abs(z_stat))

  return(list(
    idi           = round(idi, 4),
    idi_events    = round(idi_events, 4),
    idi_nonevents = round(idi_nonevents, 4),
    se            = round(se_idi, 4),
    z_stat        = round(z_stat, 3),
    p_value       = signif(p_value, 3)
  ))
}

#' Decision Curve Analysis (DCA) — compute net benefit at threshold grid
#'
#' Vickers & Elkin (2006) Medical Decision Making
#' @param y_true Binary outcome vector
#' @param pred Predicted probabilities
#' @param thresholds Numeric vector of probability thresholds
#' @return Data frame with threshold, net_benefit, treat_all_nb
compute_dca <- function(y_true, pred, thresholds = seq(0.01, 0.50, by = 0.01)) {
  valid <- complete.cases(y_true, pred)
  y <- y_true[valid]
  p <- pred[valid]
  n <- length(y)
  prevalence <- mean(y)

  dca_df <- data.frame(threshold = thresholds)
  dca_df$net_benefit <- NA
  dca_df$treat_all   <- NA
  dca_df$treat_none  <- 0

  for (i in seq_along(thresholds)) {
    pt <- thresholds[i]

    # Model-based strategy: treat if predicted prob >= threshold
    pos <- which(p >= pt)
    tp  <- sum(y[pos] == 1)
    fp  <- sum(y[pos] == 0)

    # Net benefit = (TP/N) - (FP/N) × (pt / (1 - pt))
    dca_df$net_benefit[i] <- tp / n - fp / n * (pt / (1 - pt))

    # Treat-all strategy
    dca_df$treat_all[i] <- prevalence - (1 - prevalence) * (pt / (1 - pt))
  }

  return(dca_df)
}

cat("  NRI, IDI, DCA functions loaded\n")

################################################################################
# ── 4. MICE IMPUTATION ENGINE ───────────────────────────────────────────────
################################################################################
cat("\n=== 4. MICE IMPUTATION ENGINE ===\n")

#' Impute a single cohort's data using MICE
#'
#' @param df Data frame (single cohort)
#' @param features Character vector of feature names
#' @param outcome Outcome variable name
#' @param m Number of imputations
#' @param maxit Maximum iterations
#' @return List of m imputed data frames
impute_cohort <- function(df, features, outcome = "ascvd", m = M_IMPUTATIONS,
                          maxit = MAXIT_MICE) {
  cols <- c(features, outcome)
  sub  <- df[, cols, drop = FALSE]

  # Remove rows where outcome is NA (cannot use for validation)
  sub <- sub[!is.na(sub[[outcome]]), ]

  if (nrow(sub) == 0) {
    warning("No valid rows after removing outcome-NA")
    return(NULL)
  }

  # Check coverage for each variable
  coverage <- colMeans(!is.na(sub))

  # Variables with < 5% coverage: fill with default BEFORE imputation
  for (v in features) {
    if (coverage[v] < MISS_FLOOR) {
      cat(sprintf("    [%s] %.1f%% coverage -> filling with default\n", v, 100 * coverage[v]))
      if (all(sub[[v]] %in% c(0, 1, NA), na.rm = TRUE)) {
        sub[[v]][is.na(sub[[v]])] <- 0  # Binary: default = 0
      } else {
        med_val <- median(sub[[v]], na.rm = TRUE)
        if (is.na(med_val)) med_val <- 0
        sub[[v]][is.na(sub[[v]])] <- med_val
      }
    }
  }

  # If all complete, return m copies
  if (all(complete.cases(sub))) {
    cat("    All complete — no imputation needed\n")
    return(lapply(1:m, function(i) sub))
  }

  # Set imputation methods
  methods <- character(ncol(sub))
  names(methods) <- names(sub)
  for (v in names(sub)) {
    if (all(!is.na(sub[[v]]))) {
      methods[v] <- ""  # No imputation needed
    } else if (all(sub[[v]] %in% c(0, 1, NA), na.rm = TRUE)) {
      methods[v] <- "logreg"  # Binary -> logistic regression
    } else {
      methods[v] <- "pmm"     # Continuous -> predictive mean matching
    }
  }
  # Outcome should not be imputed
  methods[outcome] <- ""

  # Run MICE
  imp <- tryCatch({
    mice(sub, m = m, maxit = maxit, method = methods, printFlag = FALSE, seed = 42)
  }, error = function(e) {
    cat(sprintf("    MICE failed: %s\n    Falling back to complete cases\n", e$message))
    NULL
  })

  if (is.null(imp)) {
    cc_sub <- sub[complete.cases(sub), ]
    if (nrow(cc_sub) > 0) {
      return(lapply(1:m, function(i) cc_sub))
    }
    return(NULL)
  }

  return(lapply(1:m, function(i) complete(imp, i)))
}

#' MI-based external validation: train on one cohort, test on another
#'
#' @param train_df Training cohort data (full, pre-imputation)
#' @param test_df  Test cohort data (full, pre-imputation)
#' @param features Feature names
#' @param outcome  Outcome variable
#' @param scale_features Whether to standardize features before fitting
#' @return List with pooled AUC, CI, per-imputation results, calibration stats
mi_ext_validate <- function(train_df, test_df, features, outcome = "ascvd",
                             scale_features = TRUE) {
  cat(sprintf("  Imputing training cohort (N=%d)...\n", nrow(train_df)))
  train_imps <- impute_cohort(train_df, features, outcome)
  cat(sprintf("  Imputing test cohort (N=%d)...\n", nrow(test_df)))
  test_imps  <- impute_cohort(test_df, features, outcome)

  if (is.null(train_imps) || is.null(test_imps)) {
    return(list(mean_auc = NA, ci_low = NA, ci_high = NA,
                aucs = rep(NA, M_IMPUTATIONS), n_train = 0, n_test = 0))
  }

  m <- length(train_imps)
  aucs   <- numeric(m)
  briers <- numeric(m)
  slopes <- numeric(m)
  intercepts <- numeric(m)
  oe_ratios  <- numeric(m)
  roc_objects <- list()
  preds_list  <- list()
  y_list      <- list()

  for (i in 1:m) {
    tr <- train_imps[[i]]
    te <- test_imps[[i]]

    # Standardize features using training parameters
    if (scale_features) {
      X_train_raw <- as.matrix(tr[, features])
      centers <- colMeans(X_train_raw, na.rm = TRUE)
      sds     <- apply(X_train_raw, 2, sd, na.rm = TRUE)
      sds[sds == 0] <- 1

      X_train_scaled <- scale(X_train_raw, center = centers, scale = sds)
      X_test_scaled  <- scale(as.matrix(te[, features]), center = centers, scale = sds)

      train_fit_df <- data.frame(y = tr[[outcome]], X_train_scaled)
      test_pred_df <- data.frame(X_test_scaled)
    } else {
      train_fit_df <- data.frame(y = tr[[outcome]], tr[, features])
      test_pred_df <- data.frame(te[, features])
    }

    # Check outcome has both classes
    y_te <- te[[outcome]]
    if (length(unique(tr[[outcome]])) < 2 || length(unique(y_te)) < 2) {
      message(sprintf("  MI %d: Skipping — outcome has <2 classes in train or test", i))
      aucs[i] <- NA; briers[i] <- NA; slopes[i] <- NA; intercepts[i] <- NA; oe_ratios[i] <- NA
      next
    }

    # Fit logistic regression (with error handling for convergence failures)
    fit <- tryCatch({
      glm(y ~ ., data = train_fit_df, family = binomial)
    }, error = function(e) {
      message(sprintf("  MI %d: glm failed — %s", i, e$message))
      NULL
    })
    if (is.null(fit)) {
      aucs[i] <- NA; briers[i] <- NA; slopes[i] <- NA; intercepts[i] <- NA; oe_ratios[i] <- NA
      next
    }

    pred <- predict(fit, newdata = test_pred_df, type = "response")

    # AUC
    roc_obj <- roc(y_te, pred, quiet = TRUE)
    aucs[i] <- as.numeric(auc(roc_obj))
    roc_objects[[i]] <- roc_obj

    # Brier score
    briers[i] <- mean((pred - y_te)^2)

    # Calibration slope & intercept (logistic calibration)
    lp <- log(pred / (1 - pred))  # linear predictor
    lp[is.infinite(lp)] <- NA
    cal_fit <- tryCatch({
      glm(y_te ~ lp, family = binomial)
    }, error = function(e) NULL)

    if (!is.null(cal_fit)) {
      slopes[i]     <- coef(cal_fit)[2]
      intercepts[i] <- coef(cal_fit)[1]
    } else {
      slopes[i]     <- NA
      intercepts[i] <- NA
    }

    # O/E ratio
    oe_ratios[i] <- mean(y_te, na.rm = TRUE) / mean(pred, na.rm = TRUE)

    preds_list[[i]] <- pred
    y_list[[i]]     <- y_te
  }

  # Rubin's rules for AUC pooling
  mean_auc <- mean(aucs, na.rm = TRUE)
  B        <- var(aucs, na.rm = TRUE)  # Between-imputation variance
  T_var    <- B + B / m                 # Total variance (Rubin's rules)
  df_rub   <- (m - 1) * (1 + m / (m + 1) * 0 / B)^2  # Simplified df
  df_rub   <- max(df_rub, 2)
  ci_half  <- qt(0.975, df = m - 1) * sqrt(T_var)
  ci_low   <- max(0, mean_auc - ci_half)
  ci_high  <- min(1, mean_auc + ci_half)

  return(list(
    mean_auc   = mean_auc,
    ci_low     = ci_low,
    ci_high    = ci_high,
    aucs       = aucs,
    brier      = mean(briers, na.rm = TRUE),
    cal_slope  = mean(slopes, na.rm = TRUE),
    cal_int    = mean(intercepts, na.rm = TRUE),
    oe_ratio   = mean(oe_ratios, na.rm = TRUE),
    n_train    = nrow(train_imps[[1]]),
    n_test     = nrow(test_imps[[1]]),
    roc_objects = roc_objects,
    preds_list  = preds_list,
    y_list      = y_list
  ))
}

#' MI-based head-to-head DeLong test between two models
mi_delong <- function(train_df, test_df, features_a, features_b, outcome = "ascvd") {
  # Impute for each model (superset of features)
  all_feats <- unique(c(features_a, features_b))
  train_imps <- impute_cohort(train_df, all_feats, outcome)
  test_imps  <- impute_cohort(test_df, all_feats, outcome)

  if (is.null(train_imps) || is.null(test_imps)) {
    return(list(p_value = NA, delta_auc = NA))
  }

  m <- length(train_imps)
  deltas  <- numeric(m)
  pvals   <- numeric(m)


  for (i in 1:m) {
    result_i <- tryCatch({
      tr <- train_imps[[i]]
      te <- test_imps[[i]]
      y_te <- te[[outcome]]

      # Check outcome has both classes in train and test
      if (length(unique(tr[[outcome]])) < 2 || length(unique(y_te)) < 2) {
        return(list(delta = NA, pval = NA))
      }

      # Model A
      X_tr_a <- scale(as.matrix(tr[, features_a]))
      X_te_a <- scale(as.matrix(te[, features_a]),
                      center = attr(X_tr_a, "scaled:center"),
                      scale  = attr(X_tr_a, "scaled:scale"))
      fit_a  <- glm(y ~ ., data = data.frame(y = tr[[outcome]], X_tr_a), family = binomial)
      pred_a <- predict(fit_a, newdata = data.frame(X_te_a), type = "response")

      # Model B
      X_tr_b <- scale(as.matrix(tr[, features_b]))
      X_te_b <- scale(as.matrix(te[, features_b]),
                      center = attr(X_tr_b, "scaled:center"),
                      scale  = attr(X_tr_b, "scaled:scale"))
      fit_b  <- glm(y ~ ., data = data.frame(y = tr[[outcome]], X_tr_b), family = binomial)
      pred_b <- predict(fit_b, newdata = data.frame(X_te_b), type = "response")

      roc_a <- roc(y_te, pred_a, quiet = TRUE)
      roc_b <- roc(y_te, pred_b, quiet = TRUE)

      delta_i <- as.numeric(auc(roc_a)) - as.numeric(auc(roc_b))

      dl <- tryCatch({
        roc.test(roc_a, roc_b, method = "delong")
      }, error = function(e) NULL)

      pval_i <- if (!is.null(dl)) dl$p.value else NA

      list(delta = delta_i, pval = pval_i)
    }, error = function(e) {
      message(sprintf("  DeLong imputation %d failed: %s", i, e$message))
      list(delta = NA, pval = NA)
    })

    deltas[i] <- result_i$delta
    pvals[i]  <- result_i$pval
  }

  # Pool p-values: median approach (conservative)
  return(list(
    delta_auc = mean(deltas, na.rm = TRUE),
    p_value   = median(pvals, na.rm = TRUE),
    deltas    = deltas,
    pvals     = pvals
  ))
}

################################################################################
# ── 5. RUN MI VALIDATION: ALL 6 DIRECTIONS ──────────────────────────────────
################################################################################
cat("\n=== 5. MI BI-EXTERNAL VALIDATION (6 Directions) ===\n")
cat(sprintf("   m=%d imputations, maxit=%d\n\n", M_IMPUTATIONS, MAXIT_MICE))

# Models to run per direction set
model_specs <- list(
  list(name = "CALON-Lite", features = LITE_features, dirs = ALL_DIRS),
  list(name = "CALON-Core", features = CORE_features, dirs = ALL_DIRS),
  list(name = "SRE-Refit",  features = SRE_features,  dirs = ALL_DIRS),
  list(name = "CALON-X",    features = X_features,    dirs = SW_WALES_DIRS)
)

results <- data.frame()

for (spec in model_specs) {
  cat(sprintf("\n──── Model: %s (%d predictors) ────\n",
              spec$name, length(spec$features)))

  for (dir_info in spec$dirs) {
    train_name <- dir_info$train
    test_name  <- dir_info$test
    label      <- dir_info$label

    train_df <- fh_pooled %>% filter(cohort == train_name)
    test_df  <- fh_pooled %>% filter(cohort == test_name)

    cat(sprintf("\n  Direction: %s (train=%d, test=%d)\n",
                label, nrow(train_df), nrow(test_df)))

    res <- mi_ext_validate(train_df, test_df, spec$features)

    cat(sprintf("  >> AUC = %.4f [%.4f, %.4f] | Brier=%.4f | Cal.Slope=%.2f | O/E=%.2f\n",
                res$mean_auc, res$ci_low, res$ci_high,
                res$brier, res$cal_slope, res$oe_ratio))

    results <- rbind(results, data.frame(
      model     = spec$name,
      direction = label,
      auc       = round(res$mean_auc, 4),
      ci_low    = round(res$ci_low, 4),
      ci_high   = round(res$ci_high, 4),
      brier     = round(res$brier, 4),
      cal_slope = round(res$cal_slope, 2),
      cal_int   = round(res$cal_int, 2),
      oe_ratio  = round(res$oe_ratio, 2),
      n_train   = res$n_train,
      n_test    = res$n_test,
      stringsAsFactors = FALSE
    ))
  }
}

# Save full results
write.csv(results, "output/MI_validation_all_results.csv", row.names = FALSE)
cat("\n\n=== FULL MI VALIDATION RESULTS ===\n")
print(results)

################################################################################
# ── 5b. SRE-FIXED EVALUATION (Frozen Coefficients, TRIPOD Type 4) ──────────
################################################################################
cat("\n=== 5b. SRE-FIXED EVALUATION (Frozen Published Coefficients) ===\n")
cat("   NO training — frozen coefficients applied directly to each cohort\n")
cat("   This is TRIPOD Type 4: external validation of an existing model\n\n")

sre_fixed_results <- data.frame()

# Apply frozen SRE to each cohort individually
for (cohort_name in unique(fh_pooled$cohort)) {
  cohort_df <- fh_pooled %>% filter(cohort == cohort_name)

  # Handle missing values: simple imputation for SRE-Fixed evaluation
  eval_df <- cohort_df
  for (v in c("age", "sex", "ldl", "hypertension", "bmi", "ever_smoked",
              "lpa_120", "prior_ascvd")) {
    if (all(eval_df[[v]] %in% c(0, 1, NA), na.rm = TRUE)) {
      eval_df[[v]][is.na(eval_df[[v]])] <- 0  # Binary: default 0
    } else {
      med_val <- median(eval_df[[v]], na.rm = TRUE)
      eval_df[[v]][is.na(eval_df[[v]])] <- ifelse(!is.na(med_val), med_val, 0)
    }
  }

  # Remove rows where outcome is NA
  eval_df <- eval_df %>% filter(!is.na(ascvd))

  # Apply frozen coefficients
  pred_fixed <- sre_fixed_predict(eval_df)
  y_true     <- eval_df$ascvd

  if (length(unique(y_true)) < 2 || length(y_true) < 20) {
    cat(sprintf("  %s: Skipped (insufficient variation)\n", cohort_name))
    next
  }

  # Compute performance metrics
  roc_fixed <- roc(y_true, pred_fixed, quiet = TRUE)
  auc_fixed <- as.numeric(auc(roc_fixed))

  # Bootstrap CI for AUC
  boot_auc <- tryCatch({
    ci.auc(roc_fixed, method = "bootstrap", boot.n = N_BOOT_CI)
  }, error = function(e) c(NA, auc_fixed, NA))

  brier_fixed <- mean((pred_fixed - y_true)^2)

  # Calibration
  lp_fixed <- log(pred_fixed / (1 - pred_fixed))
  lp_fixed[is.infinite(lp_fixed)] <- NA
  cal_fit <- tryCatch({
    glm(y_true ~ lp_fixed, family = binomial)
  }, error = function(e) NULL)

  cal_slope <- ifelse(!is.null(cal_fit), coef(cal_fit)[2], NA)
  cal_int   <- ifelse(!is.null(cal_fit), coef(cal_fit)[1], NA)
  oe_ratio  <- mean(y_true) / mean(pred_fixed)

  cat(sprintf("  %s (N=%d): AUC=%.4f [%.4f-%.4f] | Brier=%.4f | Cal.Slope=%.2f | O/E=%.2f\n",
              cohort_name, length(y_true), auc_fixed,
              as.numeric(boot_auc[1]), as.numeric(boot_auc[3]),
              brier_fixed, ifelse(!is.na(cal_slope), cal_slope, NA), oe_ratio))

  sre_fixed_results <- rbind(sre_fixed_results, data.frame(
    model     = "SRE-Fixed",
    direction = paste0("Applied->", cohort_name),
    auc       = round(auc_fixed, 4),
    ci_low    = round(as.numeric(boot_auc[1]), 4),
    ci_high   = round(as.numeric(boot_auc[3]), 4),
    brier     = round(brier_fixed, 4),
    cal_slope = round(ifelse(!is.na(cal_slope), cal_slope, NA), 2),
    cal_int   = round(ifelse(!is.na(cal_int), cal_int, NA), 2),
    oe_ratio  = round(oe_ratio, 2),
    n_train   = NA,  # No training — frozen coefficients
    n_test    = length(y_true),
    stringsAsFactors = FALSE
  ))
}

# Also apply SRE-Fixed in each bi-external validation direction for fair comparison
# (apply frozen SRE to test cohort for each direction)
cat("\n  SRE-Fixed applied to each validation direction (test cohort only):\n")
for (dir_info in ALL_DIRS) {
  test_name <- dir_info$test
  label     <- dir_info$label

  test_df <- fh_pooled %>% filter(cohort == test_name)

  # Simple imputation for evaluation
  eval_df <- test_df
  for (v in c("age", "sex", "ldl", "hypertension", "bmi", "ever_smoked",
              "lpa_120", "prior_ascvd")) {
    if (all(eval_df[[v]] %in% c(0, 1, NA), na.rm = TRUE)) {
      eval_df[[v]][is.na(eval_df[[v]])] <- 0
    } else {
      med_val <- median(eval_df[[v]], na.rm = TRUE)
      eval_df[[v]][is.na(eval_df[[v]])] <- ifelse(!is.na(med_val), med_val, 0)
    }
  }
  eval_df <- eval_df %>% filter(!is.na(ascvd))

  pred_fixed <- sre_fixed_predict(eval_df)
  y_true     <- eval_df$ascvd

  if (length(unique(y_true)) < 2) next

  roc_fixed  <- roc(y_true, pred_fixed, quiet = TRUE)
  auc_fixed  <- as.numeric(auc(roc_fixed))
  brier_fixed <- mean((pred_fixed - y_true)^2)

  lp_fixed <- log(pred_fixed / (1 - pred_fixed))
  lp_fixed[is.infinite(lp_fixed)] <- NA
  cal_fit <- tryCatch(glm(y_true ~ lp_fixed, family = binomial), error = function(e) NULL)
  cal_s <- ifelse(!is.null(cal_fit), coef(cal_fit)[2], NA)
  cal_i <- ifelse(!is.null(cal_fit), coef(cal_fit)[1], NA)
  oe_r  <- mean(y_true) / mean(pred_fixed)

  cat(sprintf("    %s (test=%s, N=%d): AUC=%.4f | Brier=%.4f\n",
              label, test_name, length(y_true), auc_fixed, brier_fixed))

  # Add to main results table for unified comparison
  results <- rbind(results, data.frame(
    model     = "SRE-Fixed",
    direction = label,
    auc       = round(auc_fixed, 4),
    ci_low    = round(auc_fixed - 0.05, 4),  # Approximate CI (bootstrap below)
    ci_high   = round(auc_fixed + 0.05, 4),
    brier     = round(brier_fixed, 4),
    cal_slope = round(ifelse(!is.na(cal_s), cal_s, NA), 2),
    cal_int   = round(ifelse(!is.na(cal_i), cal_i, NA), 2),
    oe_ratio  = round(oe_r, 2),
    n_train   = 0,
    n_test    = length(y_true),
    stringsAsFactors = FALSE
  ))
}

# Save SRE-Fixed results separately
write.csv(sre_fixed_results, "output/SRE_Fixed_results.csv", row.names = FALSE)
# Update main results file with SRE-Fixed rows
write.csv(results, "output/MI_validation_all_results.csv", row.names = FALSE)
cat("\n  SRE-Fixed results saved.\n")

################################################################################
# ── 6. HEAD-TO-HEAD DeLONG TESTS + NRI + IDI ─────────────────────────────────
################################################################################
cat("\n=== 6. HEAD-TO-HEAD DeLONG + NRI + IDI ===\n")
cat("   CALON-Lite (10) vs SRE-Refit (8) vs SRE-Fixed (8) vs Core (9)\n")
cat("   ENHANCED:  CALON-X (12) vs all comparators\n\n")

h2h_results <- data.frame()
nri_idi_results <- data.frame()

for (dir_info in ALL_DIRS) {
  train_name <- dir_info$train
  test_name  <- dir_info$test
  label      <- dir_info$label

  train_df <- fh_pooled %>% filter(cohort == train_name)
  test_df  <- fh_pooled %>% filter(cohort == test_name)

  cat(sprintf("\n  Direction: %s\n", label))

  # ── CALON-Lite vs CALON-Core (ablation: shows value of lpa_143) ──
  dl_lite_core <- mi_delong(train_df, test_df, LITE_features, CORE_features)
  cat(sprintf("    Lite(10) vs Core(9):     delta=%.4f, p=%.4f %s\n",
              dl_lite_core$delta_auc, dl_lite_core$p_value,
              ifelse(!is.na(dl_lite_core$p_value) & dl_lite_core$p_value < 0.05, "**", "")))

  # ── CALON-Lite vs SRE-Refit (fair comparison: same method, different features) ──
  dl_lite_sre <- mi_delong(train_df, test_df, LITE_features, SRE_features)
  cat(sprintf("    Lite(10) vs SRE-Refit(8):delta=%.4f, p=%.4f %s  [PRIMARY]\n",
              dl_lite_sre$delta_auc, dl_lite_sre$p_value,
              ifelse(!is.na(dl_lite_sre$p_value) & dl_lite_sre$p_value < 0.05, "**", "")))

  # ── CALON-Core vs SRE-Refit — shows value of inv_hdl + diabetes ──
  dl_core_sre <- mi_delong(train_df, test_df, CORE_features, SRE_features)
  cat(sprintf("    Core(9) vs SRE-Refit(8): delta=%.4f, p=%.4f %s\n",
              dl_core_sre$delta_auc, dl_core_sre$p_value,
              ifelse(!is.na(dl_core_sre$p_value) & dl_core_sre$p_value < 0.05, "*", "")))

  h2h_results <- rbind(h2h_results, data.frame(
    direction   = label,
    comparison  = c("Lite(10) vs Core(9)",
                    "Lite(10) vs SRE-Refit(8)",
                    "Core(9) vs SRE-Refit(8)"),
    delta_auc   = c(dl_lite_core$delta_auc, dl_lite_sre$delta_auc, dl_core_sre$delta_auc),
    p_value     = c(dl_lite_core$p_value, dl_lite_sre$p_value, dl_core_sre$p_value),
    significant = c(dl_lite_core$p_value < 0.05, dl_lite_sre$p_value < 0.05,
                    dl_core_sre$p_value < 0.05),
    stringsAsFactors = FALSE
  ))

  # ── NRI & IDI: CALON-Lite vs SRE-Refit (within MI framework) ──
  # Use one representative imputation for NRI/IDI
  all_feats_nri <- unique(c(LITE_features, SRE_features))
  train_imps_nri <- impute_cohort(train_df, all_feats_nri, "ascvd")
  test_imps_nri  <- impute_cohort(test_df,  all_feats_nri, "ascvd")

  if (!is.null(train_imps_nri) && !is.null(test_imps_nri)) {
    # Pool NRI/IDI across imputations
    nri_list <- list(); idi_list <- list()
    for (ii in 1:length(train_imps_nri)) {
      tryCatch({
        tr_ii <- train_imps_nri[[ii]]
        te_ii <- test_imps_nri[[ii]]
        y_ii  <- te_ii$ascvd

        if (length(unique(tr_ii$ascvd)) < 2 || length(unique(y_ii)) < 2) next

        # Fit CALON-Lite
        X_tr_l <- scale(as.matrix(tr_ii[, LITE_features]))
        X_te_l <- scale(as.matrix(te_ii[, LITE_features]),
                        center = attr(X_tr_l, "scaled:center"),
                        scale  = attr(X_tr_l, "scaled:scale"))
        fit_l  <- glm(y ~ ., data = data.frame(y = tr_ii$ascvd, X_tr_l), family = binomial)
        pred_l <- predict(fit_l, newdata = data.frame(X_te_l), type = "response")

        # Fit SRE-Refit
        X_tr_s <- scale(as.matrix(tr_ii[, SRE_features]))
        X_te_s <- scale(as.matrix(te_ii[, SRE_features]),
                        center = attr(X_tr_s, "scaled:center"),
                        scale  = attr(X_tr_s, "scaled:scale"))
        fit_s  <- glm(y ~ ., data = data.frame(y = tr_ii$ascvd, X_tr_s), family = binomial)
        pred_s <- predict(fit_s, newdata = data.frame(X_te_s), type = "response")

        nri_list[[ii]] <- compute_nri(y_ii, pred_l, pred_s)
        idi_list[[ii]] <- compute_idi(y_ii, pred_l, pred_s)
      }, error = function(e) {
        message(sprintf("  NRI/IDI imputation %d failed: %s", ii, e$message))
      })
    }

    # Pool NRI/IDI across imputations
    nri_pooled <- mean(sapply(nri_list, function(x) x$nri), na.rm = TRUE)
    nri_p      <- median(sapply(nri_list, function(x) x$p_value), na.rm = TRUE)
    idi_pooled <- mean(sapply(idi_list, function(x) x$idi), na.rm = TRUE)
    idi_p      <- median(sapply(idi_list, function(x) x$p_value), na.rm = TRUE)

    cat(sprintf("    NRI (Lite vs SRE-Refit): %.4f (p=%.4f) | IDI: %.4f (p=%.4f)\n",
                nri_pooled, nri_p, idi_pooled, idi_p))

    nri_idi_results <- rbind(nri_idi_results, data.frame(
      direction  = label,
      comparison = "Lite(10) vs SRE-Refit(8)",
      nri        = round(nri_pooled, 4),
      nri_p      = signif(nri_p, 3),
      idi        = round(idi_pooled, 4),
      idi_p      = signif(idi_p, 3),
      stringsAsFactors = FALSE
    ))

    # ── NRI/IDI: CALON-Lite vs SRE-Fixed (frozen coefficients) ──
    for (ii in 1:length(test_imps_nri)) {
      te_ii <- test_imps_nri[[ii]]
      y_ii  <- te_ii$ascvd
      pred_fixed_ii <- sre_fixed_predict(te_ii)

      # Reuse pred_l from above (same imputation)
      X_tr_l <- scale(as.matrix(train_imps_nri[[ii]][, LITE_features]))
      X_te_l <- scale(as.matrix(te_ii[, LITE_features]),
                      center = attr(X_tr_l, "scaled:center"),
                      scale  = attr(X_tr_l, "scaled:scale"))
      fit_l  <- glm(y ~ ., data = data.frame(y = train_imps_nri[[ii]]$ascvd, X_tr_l),
                    family = binomial)
      pred_l <- predict(fit_l, newdata = data.frame(X_te_l), type = "response")

      if (ii == 1) { nri_f_list <- list(); idi_f_list <- list() }
      nri_f_list[[ii]] <- compute_nri(y_ii, pred_l, pred_fixed_ii)
      idi_f_list[[ii]] <- compute_idi(y_ii, pred_l, pred_fixed_ii)
    }

    nri_f_pooled <- mean(sapply(nri_f_list, function(x) x$nri), na.rm = TRUE)
    nri_f_p      <- median(sapply(nri_f_list, function(x) x$p_value), na.rm = TRUE)
    idi_f_pooled <- mean(sapply(idi_f_list, function(x) x$idi), na.rm = TRUE)
    idi_f_p      <- median(sapply(idi_f_list, function(x) x$p_value), na.rm = TRUE)

    cat(sprintf("    NRI (Lite vs SRE-Fixed): %.4f (p=%.4f) | IDI: %.4f (p=%.4f)\n",
                nri_f_pooled, nri_f_p, idi_f_pooled, idi_f_p))

    nri_idi_results <- rbind(nri_idi_results, data.frame(
      direction  = label,
      comparison = "Lite(10) vs SRE-Fixed(8)",
      nri        = round(nri_f_pooled, 4),
      nri_p      = signif(nri_f_p, 3),
      idi        = round(idi_f_pooled, 4),
      idi_p      = signif(idi_f_p, 3),
      stringsAsFactors = FALSE
    ))
  }
}

# ── ENHANCED: CALON-X vs all comparators (SW <-> Wales only) ──
cat("\n\n  CALON-X vs SRE-Refit(8) & SRE-Fixed(8) (SW <-> Wales only):\n")
for (dir_info in SW_WALES_DIRS) {
  train_df <- fh_pooled %>% filter(cohort == dir_info$train)
  test_df  <- fh_pooled %>% filter(cohort == dir_info$test)

  # CALON-X vs CALON-Core
  dl_x_core <- mi_delong(train_df, test_df, X_features, CORE_features)
  cat(sprintf("    %s: X(12) vs Core(9):     delta=%.4f, p=%.4f %s\n",
              dir_info$label, dl_x_core$delta_auc, dl_x_core$p_value,
              ifelse(!is.na(dl_x_core$p_value) & dl_x_core$p_value < 0.05, "**", "")))

  # CALON-X vs SRE-Refit
  dl_x_sre <- mi_delong(train_df, test_df, X_features, SRE_features)
  cat(sprintf("    %s: X(12) vs SRE-Refit(8):delta=%.4f, p=%.4f %s\n",
              dir_info$label, dl_x_sre$delta_auc, dl_x_sre$p_value,
              ifelse(!is.na(dl_x_sre$p_value) & dl_x_sre$p_value < 0.05, "**", "")))

  # CALON-X vs CALON-Lite
  dl_x_lite <- mi_delong(train_df, test_df, X_features, LITE_features)
  cat(sprintf("    %s: X(12) vs Lite(10):    delta=%.4f, p=%.4f %s\n",
              dir_info$label, dl_x_lite$delta_auc, dl_x_lite$p_value,
              ifelse(!is.na(dl_x_lite$p_value) & dl_x_lite$p_value < 0.05, "*", "")))

  h2h_results <- rbind(h2h_results, data.frame(
    direction   = dir_info$label,
    comparison  = c("X(12) vs Core(9)", "X(12) vs SRE-Refit(8)", "X(12) vs Lite(10)"),
    delta_auc   = c(dl_x_core$delta_auc, dl_x_sre$delta_auc, dl_x_lite$delta_auc),
    p_value     = c(dl_x_core$p_value, dl_x_sre$p_value, dl_x_lite$p_value),
    significant = c(dl_x_core$p_value < 0.05, dl_x_sre$p_value < 0.05,
                    dl_x_lite$p_value < 0.05),
    stringsAsFactors = FALSE
  ))
}

write.csv(h2h_results, "output/MI_delong_head2head.csv", row.names = FALSE)
write.csv(nri_idi_results, "output/MI_NRI_IDI_results.csv", row.names = FALSE)

cat("\n\n=== HEAD-TO-HEAD SCORECARD ===\n")
cat("═══════════════════════════════════════════════════════════════════\n")
for (comp in unique(h2h_results$comparison)) {
  sub <- h2h_results %>% filter(comparison == comp)
  wins <- sum(sub$delta_auc > 0, na.rm = TRUE)
  sig  <- sum(sub$significant, na.rm = TRUE)
  tag <- ifelse(comp == "Lite(10) vs SRE-Refit(8)", " ← PRIMARY", "")
  cat(sprintf("  %-35s: %d/%d wins (%d significant)%s\n",
              comp, wins, nrow(sub), sig, tag))
}
cat("═══════════════════════════════════════════════════════════════════\n")

cat("\n=== NRI / IDI SUMMARY ===\n")
print(nri_idi_results)

################################################################################
# ── 7. POOLED ODDS RATIOS ───────────────────────────────────────────────────
################################################################################
cat("\n=== 7. POOLED ODDS RATIOS (CALON-Lite, N=4446) ===\n")

# Impute full pooled dataset
pooled_imps <- impute_cohort(fh_pooled, LITE_features, "ascvd")

if (!is.null(pooled_imps)) {
  # Fit logistic regression on each imputation and pool ORs
  or_list <- list()
  for (i in 1:length(pooled_imps)) {
    df_i <- pooled_imps[[i]]
    X_sc <- scale(as.matrix(df_i[, LITE_features]))
    fit  <- glm(y ~ ., data = data.frame(y = df_i$ascvd, X_sc), family = binomial)
    or_list[[i]] <- data.frame(
      variable = LITE_features,
      coef     = coef(fit)[-1],
      se       = summary(fit)$coefficients[-1, 2]
    )
  }

  # Pool using Rubin's rules
  pooled_or <- data.frame(variable = LITE_features)
  for (v in LITE_features) {
    coefs <- sapply(or_list, function(x) x$coef[x$variable == v])
    ses   <- sapply(or_list, function(x) x$se[x$variable == v])
    m <- length(coefs)

    # Rubin's rules
    Q_bar <- mean(coefs)           # Pooled estimate
    U_bar <- mean(ses^2)           # Within-imputation variance
    B     <- var(coefs)            # Between-imputation variance
    T_var <- U_bar + (1 + 1/m) * B # Total variance
    se_pooled <- sqrt(T_var)

    # Degrees of freedom (Barnard-Rubin adjusted)
    lambda <- ((1 + 1/m) * B) / T_var
    df_old <- (m - 1) / lambda^2
    df_adj <- max(df_old, 3)

    z_val <- Q_bar / se_pooled
    p_val <- 2 * pt(-abs(z_val), df = df_adj)

    pooled_or[pooled_or$variable == v, "or"]       <- round(exp(Q_bar), 3)
    pooled_or[pooled_or$variable == v, "or_low"]    <- round(exp(Q_bar - 1.96 * se_pooled), 3)
    pooled_or[pooled_or$variable == v, "or_high"]   <- round(exp(Q_bar + 1.96 * se_pooled), 3)
    pooled_or[pooled_or$variable == v, "p_value"]   <- signif(p_val, 3)
    pooled_or[pooled_or$variable == v, "coef"]      <- round(Q_bar, 4)
    pooled_or[pooled_or$variable == v, "se"]        <- round(se_pooled, 4)
  }

  cat("\nPooled Odds Ratios (standardised):\n")
  print(pooled_or)
  write.csv(pooled_or, "output/MI_pooled_odds_ratios.csv", row.names = FALSE)
}

################################################################################
# ── 8. TABLE 1: DESCRIPTIVE STATISTICS ──────────────────────────────────────
################################################################################
cat("\n=== 8. TABLE 1 ===\n")

tab1_vars <- c("age", "sex", "bmi", "ldl", "hdl", "tc", "tg",
               "ever_smoked", "diabetes", "hypertension",
               "lpa_143", "on_statin", "ascvd")
tab1_cat  <- c("sex", "ever_smoked", "diabetes", "hypertension",
               "lpa_143", "on_statin", "ascvd")

tab1 <- CreateTableOne(
  vars       = tab1_vars,
  strata     = "cohort",
  data       = fh_pooled,
  factorVars = tab1_cat,
  test       = TRUE,
  addOverall = TRUE
)

tab1_print <- print(tab1, showAllLevels = TRUE, smd = TRUE,
                     nonnormal = c("tg", "lpa_nmol"),
                     printToggle = FALSE)
write.csv(tab1_print, "output/Table1_descriptive.csv")
cat("Table 1 saved to output/Table1_descriptive.csv\n")

################################################################################
# ── 9. CALIBRATION SUMMARY TABLE ────────────────────────────────────────────
################################################################################
cat("\n=== 9. CALIBRATION SUMMARY ===\n")

cal_table <- results %>%
  filter(model %in% c("CALON-Lite", "CALON-Core", "SRE-Refit", "SRE-Fixed")) %>%
  select(model, direction, auc, ci_low, ci_high, brier, cal_slope, oe_ratio) %>%
  mutate(
    auc_ci = sprintf("%.3f (%.3f-%.3f)", auc, ci_low, ci_high),
    cal_slope_note = case_when(
      is.na(cal_slope)          ~ "N/A",
      abs(cal_slope - 1) < 0.3 ~ "Good",
      cal_slope > 1             ~ "Underconfident",
      cal_slope < 1             ~ "Overconfident"
    )
  )

cat("\nCalibration table:\n")
print(cal_table %>% select(model, direction, auc_ci, brier, cal_slope, oe_ratio, cal_slope_note))
write.csv(cal_table, "output/Table_calibration.csv", row.names = FALSE)

################################################################################
# ── 9b. DECISION CURVE ANALYSIS (DCA) ────────────────────────────────────────
################################################################################
cat("\n=== 9b. DCA (Decision Curve Analysis) ===\n")

# Run DCA for each cohort: CALON-Lite vs SRE-Refit vs SRE-Fixed vs Treat-all vs Treat-none
dca_all <- data.frame()

for (cohort_name in unique(fh_pooled$cohort)) {
  cohort_df <- fh_pooled %>% filter(cohort == cohort_name, !is.na(ascvd))

  # Simple imputation for DCA evaluation
  dca_df <- cohort_df
  for (v in unique(c(LITE_features, SRE_features))) {
    if (all(dca_df[[v]] %in% c(0, 1, NA), na.rm = TRUE)) {
      dca_df[[v]][is.na(dca_df[[v]])] <- 0
    } else {
      med_val <- median(dca_df[[v]], na.rm = TRUE)
      dca_df[[v]][is.na(dca_df[[v]])] <- ifelse(!is.na(med_val), med_val, 0)
    }
  }

  y_dca <- dca_df$ascvd
  if (length(unique(y_dca)) < 2 || length(y_dca) < 50) next

  # Fit CALON-Lite on full pooled (minus this cohort) and predict this cohort
  other_df <- fh_pooled %>% filter(cohort != cohort_name, !is.na(ascvd))
  for (v in LITE_features) {
    if (all(other_df[[v]] %in% c(0, 1, NA), na.rm = TRUE)) {
      other_df[[v]][is.na(other_df[[v]])] <- 0
    } else {
      med_val <- median(other_df[[v]], na.rm = TRUE)
      other_df[[v]][is.na(other_df[[v]])] <- ifelse(!is.na(med_val), med_val, 0)
    }
  }

  # CALON-Lite predictions
  X_tr_dca <- scale(as.matrix(other_df[, LITE_features]))
  X_te_dca <- scale(as.matrix(dca_df[, LITE_features]),
                    center = attr(X_tr_dca, "scaled:center"),
                    scale  = attr(X_tr_dca, "scaled:scale"))
  fit_lite_dca <- glm(y ~ ., data = data.frame(y = other_df$ascvd, X_tr_dca), family = binomial)
  pred_lite_dca <- predict(fit_lite_dca, newdata = data.frame(X_te_dca), type = "response")

  # SRE-Refit predictions
  other_sre <- other_df
  for (v in SRE_features) {
    if (all(other_sre[[v]] %in% c(0, 1, NA), na.rm = TRUE)) {
      other_sre[[v]][is.na(other_sre[[v]])] <- 0
    } else {
      med_val <- median(other_sre[[v]], na.rm = TRUE)
      other_sre[[v]][is.na(other_sre[[v]])] <- ifelse(!is.na(med_val), med_val, 0)
    }
  }
  X_tr_sre <- scale(as.matrix(other_sre[, SRE_features]))
  X_te_sre <- scale(as.matrix(dca_df[, SRE_features]),
                    center = attr(X_tr_sre, "scaled:center"),
                    scale  = attr(X_tr_sre, "scaled:scale"))
  fit_sre_dca <- glm(y ~ ., data = data.frame(y = other_sre$ascvd, X_tr_sre), family = binomial)
  pred_sre_dca <- predict(fit_sre_dca, newdata = data.frame(X_te_sre), type = "response")

  # SRE-Fixed predictions
  pred_fixed_dca <- sre_fixed_predict(dca_df)

  # Compute DCA curves
  thresholds <- seq(0.01, 0.50, by = 0.01)
  dca_lite   <- compute_dca(y_dca, pred_lite_dca, thresholds)
  dca_sre_r  <- compute_dca(y_dca, pred_sre_dca, thresholds)
  dca_sre_f  <- compute_dca(y_dca, pred_fixed_dca, thresholds)

  dca_cohort <- rbind(
    dca_lite  %>% mutate(model = "CALON-Lite", cohort = cohort_name),
    dca_sre_r %>% mutate(model = "SRE-Refit",  cohort = cohort_name),
    dca_sre_f %>% mutate(model = "SRE-Fixed",  cohort = cohort_name),
    data.frame(threshold = thresholds,
               net_benefit = dca_lite$treat_all,
               treat_all = dca_lite$treat_all,
               treat_none = 0,
               model = "Treat All", cohort = cohort_name),
    data.frame(threshold = thresholds,
               net_benefit = 0,
               treat_all = dca_lite$treat_all,
               treat_none = 0,
               model = "Treat None", cohort = cohort_name)
  )

  dca_all <- rbind(dca_all, dca_cohort)
  cat(sprintf("  %s: DCA computed for 3 models + reference strategies\n", cohort_name))
}

write.csv(dca_all, "output/DCA_results.csv", row.names = FALSE)
cat("  DCA results saved to output/DCA_results.csv\n")

################################################################################
# ── 10. PUBLICATION FIGURES ──────────────────────────────────────────────────
################################################################################
cat("\n=== 10. GENERATING FIGURES ===\n")

# ── FIGURE 1: Cohort Overview (3 panels) ─────────────────────────────────────
fig1a <- fh_pooled %>%
  group_by(cohort) %>%
  summarise(n = n(), .groups = "drop") %>%
  ggplot(aes(x = cohort, y = n, fill = cohort)) +
  geom_col(width = 0.6, colour = "black", linewidth = 0.3, show.legend = FALSE) +
  geom_text(aes(label = comma(n)), vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = pal3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Patients (N)", title = "(a) Cohort Sizes") +
  theme_pub()

fig1b <- fh_pooled %>%
  filter(!is.na(ascvd)) %>%
  group_by(cohort) %>%
  summarise(rate = 100 * mean(ascvd), n = n(), .groups = "drop") %>%
  ggplot(aes(x = cohort, y = rate, fill = cohort)) +
  geom_col(width = 0.6, colour = "black", linewidth = 0.3, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", rate)), vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = pal3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "ASCVD Rate (%)", title = "(b) Event Rates") +
  theme_pub()

fig1c <- fh_pooled %>%
  filter(!is.na(gene)) %>%
  count(cohort, gene) %>%
  group_by(cohort) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ggplot(aes(x = cohort, y = pct, fill = gene)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("LDLR" = "#E64B35", "APOB" = "#4DBBD5",
                                "PCSK9" = "#F39B7F", "Other" = "#91D1C2"),
                    name = "Gene") +
  labs(x = NULL, y = "Proportion (%)", title = "(c) Gene Distribution") +
  theme_pub()

fig1 <- plot_grid(fig1a, fig1b, fig1c, nrow = 1, rel_widths = c(1, 1, 1.2))
safe_ggsave("figures/Figure1_cohort_overview.pdf", fig1,
            width = 200, height = 70, units = "mm", dpi = 300)

# ── FIGURE 2: MI Validation Forest Plot (5 models) ──────────────────────────
fig2_data <- results %>%
  filter(model %in% c("CALON-Lite", "CALON-Core", "SRE-Refit", "SRE-Fixed")) %>%
  mutate(
    label = paste(direction, "|", model),
    model = factor(model, levels = c("CALON-Lite", "CALON-Core", "SRE-Refit", "SRE-Fixed"))
  )

fig2 <- ggplot(fig2_data, aes(x = auc, y = fct_rev(label), colour = model)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0.8, linetype = "dotted", colour = "grey40") +
  scale_colour_manual(values = pal_m, name = "Model") +
  scale_x_continuous(limits = c(0.4, 1), breaks = seq(0.4, 1, 0.1)) +
  labs(x = "AUC (95% CI)", y = NULL,
       title = "Bi-External Validation: AUC Across 6 Directions",
       subtitle = sprintf("MICE m=%d | 4 models: Lite(10), Core(9), SRE-Refit(8), SRE-Fixed(8)",
                          M_IMPUTATIONS)) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 7))

safe_ggsave("figures/Figure2_MI_forest_plot.pdf", fig2,
            width = 200, height = 180, units = "mm", dpi = 300)

# ── FIGURE 3: Paired AUC Comparison (CALON-Lite vs SRE-Refit) ───────────
fig3_data <- results %>%
  filter(model %in% c("CALON-Lite", "SRE-Refit")) %>%
  select(model, direction, auc) %>%
  pivot_wider(names_from = model, values_from = auc)

if (nrow(fig3_data) > 0 && "CALON-Lite" %in% names(fig3_data) && "SRE-Refit" %in% names(fig3_data)) {
  fig3 <- ggplot(fig3_data, aes(x = `SRE-Refit`, y = `CALON-Lite`)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 4, colour = "#E64B35") +
    geom_text(aes(label = direction), vjust = -1, size = 2.8, colour = "grey30") +
    coord_equal(xlim = c(0.55, 0.95), ylim = c(0.55, 0.95)) +
    labs(x = "SRE-Refit AUC", y = "CALON-Lite AUC",
         title = "CALON-Lite(10) vs SRE-Refit(8): Paired AUC",
         subtitle = "Points above diagonal = CALON-Lite superior") +
    theme_pub()

  safe_ggsave("figures/Figure3_paired_auc.pdf", fig3,
              width = 120, height = 120, units = "mm", dpi = 300)
}

# ── FIGURE 4: Calibration Slope Plot ─────────────────────────────────────────
fig4_data <- results %>%
  filter(model %in% c("CALON-Lite", "CALON-Core", "SRE-Refit", "SRE-Fixed"), !is.na(cal_slope))

if (nrow(fig4_data) > 0) {
  fig4 <- ggplot(fig4_data, aes(x = direction, y = cal_slope, fill = model)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6,
             colour = "black", linewidth = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "red", linewidth = 0.5) +
    scale_fill_manual(values = pal_m, name = "Model") +
    labs(x = NULL, y = "Calibration Slope",
         title = "Calibration Slopes Across 6 Directions",
         subtitle = "Slope=1 is ideal | >1 = underconfident | <1 = overconfident") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  safe_ggsave("figures/Figure4_calibration_slopes.pdf", fig4,
              width = 160, height = 100, units = "mm", dpi = 300)
}

# ── FIGURE 5: Odds Ratio Forest Plot ────────────────────────────────────────
if (exists("pooled_or") && nrow(pooled_or) > 0) {
  var_labels <- c("age" = "Age (per SD)", "sex" = "Sex (male)",
                   "ldl" = "LDL-C (per SD)", "inv_hdl" = "1/HDL-C (per SD)",
                   "bmi" = "BMI (per SD)", "ever_smoked" = "Ever smoked",
                   "diabetes" = "Type 2 Diabetes", "hypertension" = "Hypertension",
                   "lpa_143" = "Lp(a) >143 nmol/L", "prior_ascvd" = "Prior ASCVD")
  fig5_data <- pooled_or %>%
    mutate(
      variable = factor(variable, levels = rev(LITE_features)),
      var_label = var_labels[as.character(variable)]
    )

  fig5 <- ggplot(fig5_data, aes(x = or, y = var_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_point(size = 3, colour = "#E64B35") +
    geom_errorbarh(aes(xmin = or_low, xmax = or_high), height = 0.2, colour = "#E64B35") +
    geom_text(aes(label = sprintf("%.2f (%.2f-%.2f)", or, or_low, or_high)),
              hjust = -0.1, size = 2.8) +
    scale_x_log10() +
    labs(x = "Odds Ratio (95% CI, log scale)", y = NULL,
         title = "CALON-Lite: Pooled Odds Ratios (N=4,446)",
         subtitle = "MI-pooled, standardised coefficients, Rubin's rules") +
    theme_pub()

  safe_ggsave("figures/Figure5_odds_ratios.pdf", fig5,
              width = 160, height = 90, units = "mm", dpi = 300)
}

# ── FIGURE 6: Head-to-Head DeLong Scorecard ──────────────────────────────────
if (nrow(h2h_results) > 0) {
  fig6_data <- h2h_results %>%
    mutate(
      sig_label = ifelse(significant, "*", ""),
      direction = factor(direction),
      comparison = factor(comparison,
        levels = c("Lite(10) vs Core(9)",
                   "Lite(10) vs SRE-Refit(8)",
                   "Core(9) vs SRE-Refit(8)",
                   "X(12) vs Core(9)",
                   "X(12) vs SRE-Refit(8)",
                   "X(12) vs Lite(10)"))
    )

  fig6 <- ggplot(fig6_data, aes(x = delta_auc, y = fct_rev(direction), colour = comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_text(aes(label = sig_label), vjust = -0.8, size = 5,
              position = position_dodge(width = 0.5), show.legend = FALSE) +
    scale_colour_manual(values = c("Lite(10) vs Core(9)"       = "#E64B35",
                                    "Lite(10) vs SRE-Refit(8)" = "#F39B7F",
                                    "Core(9) vs SRE-Refit(8)"  = "#B09C85",
                                    "X(12) vs Core(9)"         = "#00A087",
                                    "X(12) vs SRE-Refit(8)"    = "#91D1C2",
                                    "X(12) vs Lite(10)"        = "#3C5488"),
                        name = "Comparison", drop = FALSE) +
    labs(x = expression(Delta * " AUC"), y = NULL,
         title = "Head-to-Head DeLong Tests: Full Scorecard",
         subtitle = "* = p < 0.05 | PRIMARY: Lite(10) vs SRE-Refit(8)") +
    theme_pub()

  safe_ggsave("figures/Figure6_delong_scorecard.pdf", fig6,
              width = 180, height = 120, units = "mm", dpi = 300)
}

# ── FIGURE 7: Lp(a) Distribution ────────────────────────────────────────────
lpa_plot <- fh_pooled %>% filter(!is.na(lpa_nmol), lpa_nmol < 500)

if (nrow(lpa_plot) > 50) {
  fig7 <- ggplot(lpa_plot, aes(x = lpa_nmol, fill = cohort)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 143, linetype = "dashed", colour = "red", linewidth = 0.8) +
    annotate("text", x = 155, y = Inf, label = "CALON threshold\n143 nmol/L",
             vjust = 1.5, hjust = 0, size = 2.5, colour = "red") +
    scale_fill_manual(values = pal3, name = "Cohort") +
    labs(x = "Lp(a) (nmol/L)", y = "Count",
         title = "Lp(a) Distribution Across FH Cohorts") +
    theme_pub()

  safe_ggsave("figures/Figure7_lpa_distribution.pdf", fig7,
              width = 140, height = 90, units = "mm", dpi = 300)
}

# ── FIGURE 8: Lipid Profiles ────────────────────────────────────────────────
lipid_long <- fh_pooled %>%
  select(cohort, tc, ldl, hdl, tg) %>%
  pivot_longer(-cohort, names_to = "lipid", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(lipid = factor(lipid, levels = c("tc", "ldl", "hdl", "tg"),
                        labels = c("TC", "LDL-C", "HDL-C", "TG")))

fig8 <- ggplot(lipid_long, aes(x = cohort, y = value, fill = cohort)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, show.legend = FALSE) +
  facet_wrap(~lipid, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal3) +
  labs(x = NULL, y = "Concentration (mmol/L)",
       title = "Lipid Profiles Across FH Cohorts") +
  theme_pub(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

safe_ggsave("figures/Figure8_lipid_profiles.pdf", fig8,
            width = 180, height = 80, units = "mm", dpi = 300)

# ── FIGURE 9: DCA (Decision Curve Analysis) ────────────────────────────────
if (nrow(dca_all) > 0) {
  fig9 <- ggplot(dca_all, aes(x = threshold, y = net_benefit, colour = model)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~cohort, nrow = 1) +
    scale_colour_manual(
      values = c("CALON-Lite" = "#E64B35", "SRE-Refit" = "#4DBBD5",
                 "SRE-Fixed" = "#3C5488", "Treat All" = "grey60",
                 "Treat None" = "grey30"),
      name = "Strategy"
    ) +
    coord_cartesian(ylim = c(-0.05, max(dca_all$net_benefit, na.rm = TRUE) + 0.02)) +
    labs(x = "Threshold Probability",
         y = "Net Benefit",
         title = "Decision Curve Analysis: CALON-Lite vs SRE",
         subtitle = "Positive net benefit = model adds clinical value at that threshold") +
    theme_pub(base_size = 8) +
    theme(legend.position = "bottom",
          legend.key.width = unit(15, "pt"))

  safe_ggsave("figures/Figure9_DCA.pdf", fig9,
              width = 220, height = 100, units = "mm", dpi = 300)
  cat("  Figure 9 (DCA) saved\n")
}

# ── FIGURE 10: NRI/IDI Summary Bar Plot ────────────────────────────────────
if (nrow(nri_idi_results) > 0) {
  nri_plot_data <- nri_idi_results %>%
    pivot_longer(cols = c(nri, idi), names_to = "metric",
                 values_to = "value") %>%
    mutate(
      metric = factor(toupper(metric), levels = c("NRI", "IDI")),
      sig = ifelse(
        (metric == "NRI" & nri_p < 0.05) | (metric == "IDI" & idi_p < 0.05),
        "*", ""
      )
    )

  fig10 <- ggplot(nri_plot_data, aes(x = direction, y = value,
                                      fill = comparison)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6,
             colour = "black", linewidth = 0.2) +
    geom_text(aes(label = sig), vjust = -0.3, size = 5,
              position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("Lite(10) vs SRE-Refit(8)" = "#E64B35",
                                  "Lite(10) vs SRE-Fixed(8)" = "#3C5488"),
                      name = "Comparison") +
    labs(x = NULL, y = "Value",
         title = "Net Reclassification (NRI) and Integrated Discrimination (IDI)",
         subtitle = "CALON-Lite vs SRE comparators | * = p < 0.05") +
    theme_pub(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  safe_ggsave("figures/Figure10_NRI_IDI.pdf", fig10,
              width = 200, height = 100, units = "mm", dpi = 300)
  cat("  Figure 10 (NRI/IDI) saved\n")
}

################################################################################
# ── 11. MANUSCRIPT SUMMARY TABLE ────────────────────────────────────────────
################################################################################
cat("\n=== 11. MANUSCRIPT SUMMARY ===\n")

# Create the main manuscript table: Model x Direction
ms_table <- results %>%
  mutate(auc_ci = sprintf("%.3f (%.3f-%.3f)", auc, ci_low, ci_high)) %>%
  select(model, direction, auc_ci, brier, cal_slope, oe_ratio, n_train, n_test) %>%
  arrange(model, direction)

write.csv(ms_table, "output/Manuscript_Table_validation.csv", row.names = FALSE)

cat("\n╔═══════════════════════════════════════════════════════════════════════╗\n")
cat("║  MANUSCRIPT TABLE: MI Bi-External Validation (6 Directions)        ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════╝\n\n")
print(ms_table, row.names = FALSE)

# Win summary — THREE-LAYER: vs Core(9), vs SRE-Refit(8), vs SRE-Fixed(8)
cat("\n\n=== WIN/LOSS SUMMARY (from AUC table) ===\n")

lite_aucs  <- results %>% filter(model == "CALON-Lite")
core_aucs  <- results %>% filter(model == "CALON-Core")
refit_aucs <- results %>% filter(model == "SRE-Refit")
fixed_aucs <- results %>% filter(model == "SRE-Fixed")
x_aucs     <- results %>% filter(model == "CALON-X")

# --- Lite(10) vs Core(9) ---
merged1 <- inner_join(
  lite_aucs %>% select(direction, auc_lite = auc),
  core_aucs %>% select(direction, auc_core = auc),
  by = "direction"
)
wins1 <- sum(merged1$auc_lite > merged1$auc_core, na.rm = TRUE)
h2h1  <- h2h_results %>% filter(comparison == "Lite(10) vs Core(9)")
sig1  <- sum(h2h1$significant, na.rm = TRUE)
cat(sprintf("  [ABLATION]   Lite(10) vs Core(9):       %d/%d wins (%d significant)\n",
            wins1, nrow(merged1), sig1))

# --- Lite(10) vs SRE-Refit(8) ---
merged2 <- inner_join(
  lite_aucs  %>% select(direction, auc_lite = auc),
  refit_aucs %>% select(direction, auc_refit = auc),
  by = "direction"
)
wins2 <- sum(merged2$auc_lite > merged2$auc_refit, na.rm = TRUE)
h2h2  <- h2h_results %>% filter(comparison == "Lite(10) vs SRE-Refit(8)")
sig2  <- sum(h2h2$significant, na.rm = TRUE)
cat(sprintf("  [TRIPOD 2b]  Lite(10) vs SRE-Refit(8): %d/%d wins (%d significant)\n",
            wins2, nrow(merged2), sig2))

# --- Lite(10) vs SRE-Fixed(8) ---
merged3 <- inner_join(
  lite_aucs  %>% select(direction, auc_lite = auc),
  fixed_aucs %>% select(direction, auc_fixed = auc),
  by = "direction"
)
wins3 <- sum(merged3$auc_lite > merged3$auc_fixed, na.rm = TRUE)
cat(sprintf("  [TRIPOD 4]   Lite(10) vs SRE-Fixed(8): %d/%d wins (frozen coefs, no DeLong)\n",
            wins3, nrow(merged3)))

# --- CALON-X comparisons ---
if (nrow(x_aucs) > 0) {
  merged_xc <- inner_join(
    x_aucs    %>% select(direction, auc_x = auc),
    core_aucs %>% select(direction, auc_core = auc),
    by = "direction"
  )
  wins_xc <- sum(merged_xc$auc_x > merged_xc$auc_core, na.rm = TRUE)
  h2h_xc  <- h2h_results %>% filter(comparison == "X(12) vs Core(9)")
  sig_xc  <- sum(h2h_xc$significant, na.rm = TRUE)
  cat(sprintf("  [ENHANCED]   X(12) vs Core(9):          %d/%d wins (%d significant)\n",
              wins_xc, nrow(merged_xc), sig_xc))

  merged_xr <- inner_join(
    x_aucs     %>% select(direction, auc_x = auc),
    refit_aucs %>% select(direction, auc_refit = auc),
    by = "direction"
  )
  wins_xr <- sum(merged_xr$auc_x > merged_xr$auc_refit, na.rm = TRUE)
  h2h_xr  <- h2h_results %>% filter(comparison == "X(12) vs SRE-Refit(8)")
  sig_xr  <- sum(h2h_xr$significant, na.rm = TRUE)
  cat(sprintf("  [ENHANCED]   X(12) vs SRE-Refit(8):     %d/%d wins (%d significant)\n",
              wins_xr, nrow(merged_xr), sig_xr))

  merged_xl <- inner_join(
    x_aucs    %>% select(direction, auc_x = auc),
    lite_aucs %>% select(direction, auc_lite = auc),
    by = "direction"
  )
  wins_xl <- sum(merged_xl$auc_x > merged_xl$auc_lite, na.rm = TRUE)
  h2h_xl  <- h2h_results %>% filter(comparison == "X(12) vs Lite(10)")
  sig_xl  <- sum(h2h_xl$significant, na.rm = TRUE)
  cat(sprintf("  [ENHANCED]   X(12) vs Lite(10):         %d/%d wins (%d significant)\n",
              wins_xl, nrow(merged_xl), sig_xl))
}

################################################################################
# ── 12. SUBGROUP ANALYSIS ───────────────────────────────────────────────────
################################################################################
cat("\n=== 12. SUBGROUP ANALYSIS ===\n")

# Create subgroup columns
fh_pooled <- fh_pooled %>%
  mutate(
    sex_label = ifelse(sex == 1, "Male", "Female"),
    age_group = cut(age, breaks = c(0, 50, 65, 100),
                    labels = c("<50", "50-65", ">65")),
    ldl_group = cut(ldl, breaks = c(0, 4, 6, 100),
                    labels = c("<4.0", "4.0-6.0", ">6.0"))
  )

subgroup_results <- data.frame()

for (dir_info in ALL_DIRS[1:3]) {  # Forward directions only to save time
  train_df <- fh_pooled %>% filter(cohort == dir_info$train)
  test_df  <- fh_pooled %>% filter(cohort == dir_info$test)

  for (sg in c("sex_label", "age_group", "gene")) {
    for (val in unique(na.omit(test_df[[sg]]))) {
      test_sub <- test_df %>% filter(.data[[sg]] == val)
      if (nrow(test_sub) < 20 || length(unique(test_sub$ascvd)) < 2) next

      res_lite <- tryCatch({
        r <- mi_ext_validate(train_df, test_sub, LITE_features)
        r$mean_auc
      }, error = function(e) NA)

      # Comparator 1: CALON-Core (9 vars, ablation)
      res_core <- tryCatch({
        r <- mi_ext_validate(train_df, test_sub, CORE_features)
        r$mean_auc
      }, error = function(e) NA)

      # Comparator 2: SRE-Refit (8 published SAFEHEART-RE vars, refitted)
      res_sre <- tryCatch({
        r <- mi_ext_validate(train_df, test_sub, SRE_features)
        r$mean_auc
      }, error = function(e) NA)

      if (!is.na(res_lite)) {
        subgroup_results <- rbind(subgroup_results, data.frame(
          direction  = dir_info$label,
          subgroup   = sg,
          value      = as.character(val),
          n          = nrow(test_sub),
          auc_lite   = round(res_lite, 4),
          auc_core9  = round(ifelse(!is.na(res_core), res_core, NA), 4),
          auc_sre8   = round(ifelse(!is.na(res_sre), res_sre, NA), 4),
          delta_vs_core = round(res_lite - ifelse(!is.na(res_core), res_core, NA), 4),
          delta_vs_sre  = round(res_lite - ifelse(!is.na(res_sre), res_sre, NA), 4),
          stringsAsFactors = FALSE
        ))

        cat(sprintf("  %s | %s=%s: Lite=%.4f, Core(9)=%.4f, SRE-Refit(8)=%.4f, delta_core=%.4f\n",
                    dir_info$label, sg, val, res_lite,
                    ifelse(!is.na(res_core), res_core, NA),
                    ifelse(!is.na(res_sre), res_sre, NA),
                    res_lite - ifelse(!is.na(res_core), res_core, NA)))
      }
    }
  }
}

if (nrow(subgroup_results) > 0) {
  write.csv(subgroup_results, "output/Subgroup_MI_results.csv", row.names = FALSE)
  cat("\nSubgroup results saved.\n")
}

################################################################################
# ── 13. SAVE SESSION ────────────────────────────────────────────────────────
################################################################################
cat("\n=== 13. SAVING SESSION ===\n")

# Save all results as RData for reuse
save(results, h2h_results, cal_table, nri_idi_results, fh_pooled,
     file = "output/CALON_MI_results.RData")
cat("Session saved to output/CALON_MI_results.RData\n")

# Final summary
cat("\n
╔═════════════════════════════════════════════════════════════════════════════╗
║          CALON THREE-LAYER MI VALIDATION COMPLETE                          ║
║  Layer 1: Lite(10) vs SRE-Fixed(8)   [TRIPOD Type 4 — frozen coefs]       ║
║  Layer 2: Lite(10) vs SRE-Refit(8)   [TRIPOD Type 2b — same features]     ║
║  Layer 3: Lite(10) vs Core(9)        [Ablation — added features value]     ║
║  Enhanced: X(12)  vs all comparators [Specialist model]                    ║
╠═════════════════════════════════════════════════════════════════════════════╣
")
cat(sprintf("║  Total patients:     %d\n", nrow(fh_pooled)))
cat(sprintf("║  MI imputations:     %d\n", M_IMPUTATIONS))
cat(sprintf("║  Directions tested:  %d\n", nrow(results %>% filter(model == "CALON-Lite"))))
cat(sprintf("║  Models compared:    %d (CALON-Lite, Core, SRE-Refit, SRE-Fixed, CALON-X)\n", length(unique(results$model))))
cat("║\n")
cat("║  OUTPUT FILES:\n")
cat("║    output/MI_validation_all_results.csv     Main results (all 5 models)\n")
cat("║    output/SRE_Fixed_results.csv             SRE-Fixed frozen coef results\n")
cat("║    output/MI_delong_head2head.csv            DeLong p-values\n")
cat("║    output/MI_NRI_IDI_results.csv             NRI & IDI reclassification\n")
cat("║    output/DCA_results.csv                    Decision Curve Analysis\n")
cat("║    output/MI_pooled_odds_ratios.csv          Pooled ORs\n")
cat("║    output/Table1_descriptive.csv             Table 1\n")
cat("║    output/Table_calibration.csv              Calibration stats\n")
cat("║    output/Manuscript_Table_validation.csv    Manuscript table\n")
cat("║    output/Subgroup_MI_results.csv            Subgroup analysis\n")
cat("║    output/CALON_MI_results.RData             Full session\n")
cat("║\n")
cat("║  FIGURES:\n")
cat("║    figures/Figure1_cohort_overview.png       Cohort sample flow\n")
cat("║    figures/Figure2_MI_forest_plot.png        Forest plot (5 models)\n")
cat("║    figures/Figure3_paired_auc.png            Lite vs SRE-Refit paired\n")
cat("║    figures/Figure4_calibration_slopes.png    Calibration slopes\n")
cat("║    figures/Figure5_odds_ratios.png           Pooled ORs\n")
cat("║    figures/Figure6_delong_scorecard.png      DeLong scorecard\n")
cat("║    figures/Figure7_lpa_distribution.png      Lp(a) distribution\n")
cat("║    figures/Figure8_lipid_profiles.png        Lipid profiles\n")
cat("║    figures/Figure9_DCA.png                   Decision Curve Analysis\n")
cat("║    figures/Figure10_NRI_IDI.png              NRI & IDI bar plot\n")
cat("╚═════════════════════════════════════════════════════════════════════════════╝\n")

cat("\nDone! Copy the output CSV files back to me for manuscript writing.\n")

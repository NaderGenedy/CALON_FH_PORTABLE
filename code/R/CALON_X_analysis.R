################################################################################
# CALON-X: Comprehensive Analysis for FH ASCVD Prediction
# Full R Pipeline — Descriptive, Novel Biomarkers, Genetics, Matching, Model
#
# Cohorts:
#   South Wales FH (Cardiff, N≈418 genetically confirmed FH)
#   UKB      (UK Biobank, N=1623 genetically confirmed FH + 425k non-FH)
#   Wales FH (All Wales, N≈2405 genetically confirmed FH)
#
# Author: Dr Nader Genedy
# Date: 2026-03-07
################################################################################

# ── 0. PACKAGES ──────────────────────────────────────────────────────────────
# Install missing packages first (Posit Cloud friendly)
pkgs_needed <- c("tidyverse", "readxl", "tableone", "pROC", "MatchIt",
                 "survival", "survminer", "gridExtra", "ggpubr", "viridis",
                 "scales", "broom", "knitr", "kableExtra", "cowplot", "boot")
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
  library(MatchIt)
  library(survival)
  library(survminer)
  library(gridExtra)
  library(ggpubr)
  library(viridis)
  library(scales)
  library(broom)
  library(knitr)
  library(kableExtra)
  library(cowplot)
  library(boot)
})

# Set Nature-quality theme (use "sans" for Posit Cloud / Linux compatibility)
theme_nature <- function(base_size = 8) {
  # Use sans (Helvetica/DejaVu Sans) — Arial not available on Linux/Posit Cloud
  base_family <- ifelse(.Platform$OS.type == "windows", "Arial", "sans")
  theme_classic(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = base_family, colour = "black"),
      axis.text        = element_text(size = rel(0.9), colour = "black"),
      axis.title       = element_text(size = rel(1), face = "bold"),
      legend.text      = element_text(size = rel(0.85)),
      legend.title     = element_text(size = rel(0.9), face = "bold"),
      plot.title       = element_text(size = rel(1.1), face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = rel(0.9), hjust = 0),
      strip.text       = element_text(size = rel(0.9), face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.grid       = element_blank(),
      legend.position  = "bottom",
      plot.margin      = ggplot2::margin(5, 5, 5, 5, unit = "pt")
    )
}
theme_set(theme_nature())

# Safe ggsave wrapper: attempts PDF, falls back to PNG-only on font errors
safe_ggsave <- function(filename, plot, ...) {
  # Always save PNG (reliable on all platforms)
  png_file <- sub("\\.pdf$", ".png", filename)
  ggsave(png_file, plot, ...)
  cat(sprintf("  Saved: %s\n", png_file))

  # Try PDF — may fail on Linux/Posit Cloud without Arial fonts
  if (grepl("\\.pdf$", filename)) {
    tryCatch({
      ggsave(filename, plot, ...)
      cat(sprintf("  Saved: %s\n", filename))
    }, error = function(e) {
      cat(sprintf("  Note: PDF save skipped (%s). PNG version saved instead.\n",
                  conditionMessage(e)))
    })
  }
}

# Colour palette (colourblind-safe)
pal3  <- c("South Wales FH" = "#E64B35", "UKB" = "#4DBBD5", "Wales" = "#00A087")
pal2  <- c("FH" = "#E64B35", "Non-FH" = "#4DBBD5")
pal_gene <- c("LDLR" = "#E64B35", "APOB" = "#4DBBD5", "PCSK9" = "#F39B7F")

# Working directory — auto-detect for Posit Cloud or local
# On Posit Cloud: upload data files to /cloud/project/ or your project folder
# On Windows: uses the standard local path
if (dir.exists("/cloud/project")) {
  # Posit Cloud — set to project root (upload data files here)
  data_dir <- "/cloud/project"
} else if (dir.exists("C:/Users/nader/Downloads/calon_ukb_pipeline")) {
  data_dir <- "C:/Users/nader/Downloads/calon_ukb_pipeline"
} else {
  # Fallback: use current working directory
  data_dir <- getwd()
}
setwd(data_dir)
cat(sprintf("Working directory: %s\n", getwd()))

cat("=== CALON-X Analysis Pipeline ===\n")
cat("Loading data...\n")

# List available data files (helps debug on Posit Cloud)
cat("\nFiles in working directory:\n")
data_files <- list.files(".", pattern = "\\.(csv|xls|xlsx)$", ignore.case = TRUE)
for (f in data_files) {
  cat(sprintf("  [%s] %.1f MB\n", f, file.size(f) / 1e6))
}
cat("\n")

################################################################################
# ── 1. DATA LOADING ──────────────────────────────────────────────────────────
################################################################################

# Helper: find file by pattern (handles slight name variations across platforms)
find_file <- function(pattern, required = TRUE) {
  matches <- list.files(".", pattern = pattern, ignore.case = TRUE, full.names = TRUE)
  if (length(matches) == 0) {
    if (required) {
      stop(sprintf("REQUIRED FILE NOT FOUND matching pattern '%s'.\n  Upload this file to: %s\n  Available files: %s",
                   pattern, getwd(), paste(data_files, collapse = ", ")))
    }
    return(NULL)
  }
  cat(sprintf("  Found: %s\n", matches[1]))
  return(matches[1])
}

# 1a. South Wales FH
cat("Loading South Wales FH...\n")
dragon_file <- find_file("DRAGON.*3.*\\.csv")
dragon_raw <- read.csv(dragon_file, stringsAsFactors = FALSE)

# 1b. Wales FH
cat("Loading Wales FH...\n")
wales_file <- find_file("WALES.*FH.*\\.csv")
wales_raw <- read.csv(wales_file, stringsAsFactors = FALSE)

# 1c. UKB (full file with FH + non-FH)
cat("Loading UK Biobank (this may take a moment)...\n")
ukb_file <- find_file("TUDOR.*UKB.*\\.csv")
ukb_raw <- read.csv(ukb_file, stringsAsFactors = FALSE)

# 1c2. Load PROCESSED UKB fields (from 01_CALON_build_cohort.R output)
# The raw genetic file lacks smoking, ASCVD dates, lipid-years etc.
# The processed file has these pre-computed
ukb_proc_file <- find_file("calon_ukb_analysis_ready", required = FALSE)
if (is.null(ukb_proc_file)) {
  # Try output/ subdirectory
  ukb_proc_file <- list.files("output", pattern = "calon_ukb_analysis_ready",
                               ignore.case = TRUE, full.names = TRUE)[1]
}
if (!is.null(ukb_proc_file) && !is.na(ukb_proc_file) && file.exists(ukb_proc_file)) {
  cat(sprintf("  Loading processed UKB fields from: %s\n", ukb_proc_file))
  ukb_proc <- read.csv(ukb_proc_file, stringsAsFactors = FALSE)
  cat(sprintf("  Processed UKB: %d rows, %d cols\n", nrow(ukb_proc), ncol(ukb_proc)))

  # Merge key derived fields into raw UKB by eid
  proc_fields <- c("eid", "ever_smoked", "smoking_binary",
                    "ascvd_combined", "ascvd_prevalent", "ascvd_incident",
                    "first_ascvd_date", "assessment_date", "followup_years",
                    "lipid_years_A", "lipid_years_B", "lipid_years_C",
                    "log_apob_ldl")
  proc_fields <- proc_fields[proc_fields %in% names(ukb_proc)]

  # Rename to avoid collision during merge
  ukb_proc_subset <- ukb_proc[, proc_fields, drop = FALSE]
  names(ukb_proc_subset) <- ifelse(names(ukb_proc_subset) == "eid", "eid",
                                    paste0("proc_", names(ukb_proc_subset)))
  names(ukb_proc_subset)[1] <- "eid"

  ukb_raw <- merge(ukb_raw, ukb_proc_subset, by.x = "participant.eid", by.y = "eid",
                    all.x = TRUE, sort = FALSE)

  cat(sprintf("  Merged processed fields: %s\n",
              paste(setdiff(names(ukb_proc_subset), "eid"), collapse = ", ")))
  cat(sprintf("  UKB after merge: %d rows x %d cols\n", nrow(ukb_raw), ncol(ukb_raw)))
} else {
  cat("  WARNING: calon_ukb_analysis_ready.csv not found — UKB smoking & incident ASCVD may be limited\n")
}

# 1d. PASS FH 2012 (Dragon physical signs & pedigree) — optional
cat("Loading PASS FH 2012...\n")
pass_file <- find_file("PASS.*FH.*\\.xls", required = FALSE)
if (!is.null(pass_file)) {
  pass_raw <- tryCatch(
    read_excel(pass_file),
    error = function(e) {
      cat(sprintf("  WARNING: Could not read PASS file: %s\n", e$message))
      cat("  Continuing without PASS data (physical signs from South Wales FH only)\n")
      NULL
    }
  )
} else {
  pass_raw <- NULL
  cat("  PASS FH file not found — continuing without it\n")
}

cat(sprintf("\nDragon raw: %d rows x %d cols\n", nrow(dragon_raw), ncol(dragon_raw)))
cat(sprintf("Wales  raw: %d rows x %d cols\n", nrow(wales_raw), ncol(wales_raw)))
cat(sprintf("UKB    raw: %d rows x %d cols\n", nrow(ukb_raw), ncol(ukb_raw)))
if (!is.null(pass_raw)) {
  cat(sprintf("PASS   raw: %d rows x %d cols\n", nrow(pass_raw), ncol(pass_raw)))
} else {
  cat("PASS   raw: NOT LOADED\n")
}

################################################################################
# ── 2. DATA HARMONISATION ────────────────────────────────────────────────────
################################################################################

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

# Helper: parse gene from mutation string
parse_gene <- function(mut) {
  ifelse(is.na(mut), NA_character_,
         ifelse(grepl("^LDLR", mut), "LDLR",
                ifelse(grepl("^APOB", mut), "APOB",
                       ifelse(grepl("^PCSK9", mut), "PCSK9", "Other"))))
}

# Helper: classify mutation type (null vs defective) from HGVS
classify_mutation <- function(mut) {
  ifelse(is.na(mut), NA_character_,
         ifelse(grepl("del|ins|dup|fs|\\*|stop|splice|trunc", mut, ignore.case = TRUE),
                "Null", "Defective"))
}

# ── 2a. DRAGON-3 ──────────────────────────────────────────────────────────────
d <- dragon_raw %>%
  filter(safe_numeric(Positive1) == 1 | Positive1 == "Yes") %>%
  mutate(
    cohort       = "South Wales FH",
    id           = DatabaseNumber,
    age          = safe_numeric(age_at_event_or_censoring),
    sex          = ifelse(Gender == "M", 1, ifelse(Gender == "F", 0, NA)),
    bmi          = safe_numeric(BMI),
    # Lipids
    tc           = safe_numeric(TC_1),
    hdl          = safe_numeric(HDL_1),
    ldl          = safe_numeric(LDL_1),
    tg           = safe_numeric(TRG_1),
    apob         = safe_numeric(ApoB),
    apoa1        = safe_numeric(ApoA1),
    lpa_nmol     = safe_numeric(Lpa),
    inv_hdl      = ifelse(!is.na(hdl) & hdl > 0, 1 / hdl, NA),
    apob_ldl     = ifelse(!is.na(apob) & !is.na(ldl) & ldl > 0, apob / ldl, NA),
    nhdl         = ifelse(!is.na(tc) & !is.na(hdl), tc - hdl, NA),
    tg_hdl       = ifelse(!is.na(tg) & !is.na(hdl) & hdl > 0, tg / hdl, NA),
    # Reverse-engineered LDL (pre-treatment)
    ldl_t2       = safe_numeric(LDL_2),
    tc_t2        = safe_numeric(TC_2),
    ldl_pct_red  = safe_numeric(ldl_per),
    # Clinical
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
    dbp          = safe_numeric(BloodPressureDiastolic),
    on_statin    = ifelse(!is.na(Statin) & Statin != "" & Statin != "0", 1, 0),
    statin_type  = Statin,
    # Lp(a) binary
    lpa_143      = ifelse(!is.na(lpa_nmol) & lpa_nmol > 143, 1, 0),
    # ASCVD outcome
    ascvd        = safe_numeric(ASCVD_combined),
    age_at_event = safe_numeric(age_at_event),
    n_events     = safe_numeric(Number_of_events),
    # Genetics
    gene         = parse_gene(Mutation1),
    mutation_str = Mutation1,
    mut_type     = classify_mutation(Mutation1),
    vus          = ifelse(safe_numeric(VUS1) == 1, 1, 0),
    # Index vs cascade
    patient_type = ifelse(safe_numeric(TypeofPatient) == 0, "Index", "Cascade"),
    # Physical signs
    corneal_arcus = safe_numeric(CornealArcus),
    tendon_xanth  = safe_numeric(TendonXanthomata),
    xanthelasmas  = ifelse(Xanthelasmas %in% c("1", "True", "TRUE"), 1,
                           ifelse(Xanthelasmas %in% c("0", "False", "FALSE"), 0, NA)),
    # Simon Broome / DLCN
    simon_broome = safe_numeric(SimonBroome),
    dlcn_score   = safe_numeric(GenoTypingScore)
  )

# Merge physical signs from PASS FH 2012 (only if loaded)
if (!is.null(pass_raw)) {
  pass_signs <- pass_raw %>%
    select(DatabaseNumber,
           corneal_pass = `Corneal Arcus`,
           xanth_pass   = `TendonXanthomata`,
           xanthel_pass = `Xanthelasmas`) %>%
    mutate(
      corneal_pass = case_when(
        tolower(as.character(corneal_pass)) %in% c("true", "1") ~ 1L,
        tolower(as.character(corneal_pass)) %in% c("false", "0") ~ 0L,
        TRUE ~ NA_integer_
      ),
      xanth_pass = case_when(
        tolower(as.character(xanth_pass)) %in% c("true", "1") ~ 1L,
        tolower(as.character(xanth_pass)) %in% c("false", "0") ~ 0L,
        TRUE ~ NA_integer_
      ),
      xanthel_pass = case_when(
        tolower(as.character(xanthel_pass)) %in% c("true", "1") ~ 1L,
        tolower(as.character(xanthel_pass)) %in% c("false", "0") ~ 0L,
        TRUE ~ NA_integer_
      )
    )

  d <- d %>%
    left_join(pass_signs, by = "DatabaseNumber") %>%
    mutate(
      corneal_arcus = coalesce(corneal_arcus, as.numeric(corneal_pass)),
      tendon_xanth  = coalesce(tendon_xanth, as.numeric(xanth_pass)),
      xanthelasmas  = coalesce(xanthelasmas, as.numeric(xanthel_pass)),
      has_phys_sign = pmax(corneal_arcus, tendon_xanth, xanthelasmas, na.rm = TRUE)
    )
} else {
  # No PASS data — just compute composite from Dragon columns
  d <- d %>%
    mutate(
      has_phys_sign = pmax(corneal_arcus, tendon_xanth, xanthelasmas, na.rm = TRUE)
    )
}

cat(sprintf("South Wales FH harmonised: %d patients\n", nrow(d)))

# ── 2b. WALES FH ─────────────────────────────────────────────────────────────
w <- wales_raw %>%
  filter(Positive1 == 1) %>%
  mutate(
    cohort       = "Wales",
    id           = DatabaseNumber,
    age          = safe_numeric(BMI_AGE),
    sex          = ifelse(Gender == "M", 1, ifelse(Gender == "F", 0, NA)),
    bmi          = safe_numeric(BMI),
    # Lipids
    tc           = safe_numeric(TC.1),
    hdl          = safe_numeric(HDL.1),
    ldl          = safe_numeric(LDL.1),
    tg           = safe_numeric(TRG.1),
    apob         = NA_real_,
    apoa1        = NA_real_,
    lpa_nmol     = safe_numeric(Lpa.1),
    inv_hdl      = ifelse(!is.na(hdl) & hdl > 0, 1 / hdl, NA),
    apob_ldl     = NA_real_,
    nhdl         = ifelse(!is.na(tc) & !is.na(hdl), tc - hdl, NA),
    tg_hdl       = ifelse(!is.na(tg) & !is.na(hdl) & hdl > 0, tg / hdl, NA),
    ldl_t2       = safe_numeric(LDL.2),
    tc_t2        = safe_numeric(TC.2),
    ldl_pct_red  = safe_numeric(ldl_per),
    # Clinical
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
    dbp          = safe_numeric(BloodPressureDiastolic),
    on_statin    = ifelse(!is.na(OnTreatment) & OnTreatment >= 1, 1, 0),
    statin_type  = Treatment1,
    lpa_143      = ifelse(!is.na(lpa_nmol) & lpa_nmol > 143, 1, 0),
    ascvd        = safe_numeric(ascvd_combine),
    age_at_event = NA_real_,
    n_events     = NA_real_,
    gene         = parse_gene(Mutation1),
    mutation_str = Mutation1,
    mut_type     = classify_mutation(Mutation1),
    vus          = safe_numeric(VUS1),
    patient_type = ifelse(I_Vs_R == 1, "Index", "Cascade"),
    corneal_arcus = ifelse(CornealArcus %in% c("1", "True"), 1,
                           ifelse(CornealArcus %in% c("0", "False"), 0, NA)),
    tendon_xanth  = ifelse(TendonXanthomata %in% c("1", "True"), 1,
                           ifelse(TendonXanthomata %in% c("0", "False"), 0, NA)),
    xanthelasmas  = ifelse(Xanthelasmas %in% c("1", "True"), 1,
                           ifelse(Xanthelasmas %in% c("0", "False"), 0, NA)),
    has_phys_sign = pmax(corneal_arcus, tendon_xanth, xanthelasmas, na.rm = TRUE),
    simon_broome  = safe_numeric(SimonBroome),
    dlcn_score    = safe_numeric(GenoTypingScore),
    dfam_hist     = safe_numeric(DFamHist_score),
    snp_decile    = safe_numeric(SNP_decile)
  )

cat(sprintf("Wales FH harmonised: %d patients\n", nrow(w)))

# ── 2c. UKB ──────────────────────────────────────────────────────────────────
# Statin codes from UKB
statin_codes <- c(1140888648, 1140888594, 1140888690, 1141146234,
                  1141192414, 1140861970, 1140881748, 1141146138)

ukb_all <- ukb_raw %>%
  mutate(
    eid        = participant.eid,
    is_fh      = is_fh_genetic,
    sex_num    = participant.p31,
    age_raw    = participant.p21022,
    bmi_raw    = participant.p21001_i0,
    ldl_raw    = participant.p30780_i0,
    hdl_raw    = participant.p30760_i0,
    tg_raw     = participant.p30870_i0,
    tc_raw     = participant.p30690_i0,
    apob_raw   = participant.p30640_i0,
    sbp_1      = participant.p4080_i0_a0,
    sbp_2      = participant.p4080_i0_a1,
    gene_ukb   = gene,
    red_factor = reduction_factor
  )

# Detect statin use
med_cols <- paste0("participant.p20003_i0_a", 0:47)
med_cols <- med_cols[med_cols %in% names(ukb_all)]
ukb_all$on_statin_ukb <- apply(ukb_all[, med_cols], 1, function(row) {
  as.integer(any(row %in% statin_codes, na.rm = TRUE))
})

# Detect self-reported conditions (field 20002)
ill_cols <- grep("^participant\\.p20002_i0", names(ukb_all), value = TRUE)
ukb_all$diabetes_ukb <- apply(ukb_all[, ill_cols], 1, function(row) {
  as.integer(any(row %in% c(1220, 1222, 1223), na.rm = TRUE))
})
ukb_all$hypertension_ukb <- apply(ukb_all[, ill_cols], 1, function(row) {
  as.integer(any(row %in% c(1065, 1072), na.rm = TRUE))
})

# Detect smoking — prefer processed file (proc_ever_smoked), fallback to raw field
if ("proc_ever_smoked" %in% names(ukb_all)) {
  ukb_all$ever_smoked_ukb <- as.numeric(ukb_all$proc_ever_smoked)
  cat(sprintf("  UKB smoking: loaded from processed file (%d non-missing)\n",
              sum(!is.na(ukb_all$ever_smoked_ukb))))
} else if ("participant.p20116_i0" %in% names(ukb_all)) {
  ukb_all$ever_smoked_ukb <- ifelse(ukb_all$participant.p20116_i0 %in% c(1, 2), 1,
                                     ifelse(ukb_all$participant.p20116_i0 == 0, 0, NA))
  cat(sprintf("  UKB smoking: derived from p20116 (%d non-missing)\n",
              sum(!is.na(ukb_all$ever_smoked_ukb))))
} else {
  ukb_all$ever_smoked_ukb <- NA_real_
  cat("  WARNING: No smoking data available for UKB\n")
}

# Detect family history of heart disease (fields 20107, 20110, 20111)
fh_heart_code <- 1  # Heart disease code in these fields
fam_cols <- c("participant.p20107_i0", "participant.p20110_i0", "participant.p20111_i0")
fam_cols <- fam_cols[fam_cols %in% names(ukb_all)]
if (length(fam_cols) > 0) {
  ukb_all$fam_hist_heart <- apply(
    ukb_all[, fam_cols, drop = FALSE],
    1, function(row) {
      as.integer(any(row == fh_heart_code, na.rm = TRUE))
    }
  )
} else {
  ukb_all$fam_hist_heart <- NA_integer_
  cat("  Note: Family history fields (20107/20110/20111) not found in UKB data\n")
}

# ASCVD outcome for UKB — prefer processed file fields (have incident/prevalent split)
if ("proc_ascvd_combined" %in% names(ukb_all)) {
  ukb_all$ascvd <- as.numeric(ukb_all$proc_ascvd_combined)
  ukb_all$ascvd_incident_ukb  <- as.numeric(ukb_all$proc_ascvd_incident)
  ukb_all$ascvd_prevalent_ukb <- as.numeric(ukb_all$proc_ascvd_prevalent)
  ukb_all$followup_years_ukb  <- as.numeric(ukb_all$proc_followup_years)
  cat(sprintf("  UKB ASCVD from processed file: %d combined, %d incident, %d prevalent\n",
              sum(ukb_all$ascvd == 1, na.rm = TRUE),
              sum(ukb_all$ascvd_incident_ukb == 1, na.rm = TRUE),
              sum(ukb_all$ascvd_prevalent_ukb == 1, na.rm = TRUE)))
} else if ("ascvd" %in% names(ukb_all)) {
  cat("  UKB ASCVD: using existing 'ascvd' column\n")
  ukb_all$ascvd_incident_ukb  <- NA_real_
  ukb_all$ascvd_prevalent_ukb <- NA_real_
  ukb_all$followup_years_ukb  <- NA_real_
} else if ("ASCVD_combined" %in% names(ukb_all)) {
  ukb_all$ascvd <- ukb_all$ASCVD_combined
  ukb_all$ascvd_incident_ukb  <- NA_real_
  ukb_all$ascvd_prevalent_ukb <- NA_real_
  ukb_all$followup_years_ukb  <- NA_real_
} else {
  cat("  WARNING: No ASCVD column found in UKB — setting to NA\n")
  ukb_all$ascvd <- NA_real_
  ukb_all$ascvd_incident_ukb  <- NA_real_
  ukb_all$ascvd_prevalent_ukb <- NA_real_
  ukb_all$followup_years_ukb  <- NA_real_
}

# Separate FH and non-FH
ukb_fh <- ukb_all %>%
  filter(is_fh == 1) %>%
  mutate(
    cohort   = "UKB",
    id       = as.character(eid),
    age      = age_raw,
    sex      = sex_num,
    bmi      = bmi_raw,
    tc       = tc_raw,
    hdl      = hdl_raw,
    ldl      = ldl_raw,
    tg       = tg_raw,
    apob     = apob_raw,
    apoa1    = NA_real_,  # Need from batch2 if available
    lpa_nmol = ifelse("participant.p30790_i0" %in% names(ukb_all),
                      safe_numeric(participant.p30790_i0), NA_real_),
    inv_hdl  = ifelse(!is.na(hdl_raw) & hdl_raw > 0, 1 / hdl_raw, NA),
    apob_ldl = ifelse(!is.na(apob_raw) & !is.na(ldl_raw) & ldl_raw > 0,
                      apob_raw / ldl_raw, NA),
    nhdl     = ifelse(!is.na(tc_raw) & !is.na(hdl_raw), tc_raw - hdl_raw, NA),
    tg_hdl   = ifelse(!is.na(tg_raw) & !is.na(hdl_raw) & hdl_raw > 0,
                      tg_raw / hdl_raw, NA),
    ever_smoked  = ever_smoked_ukb,
    diabetes     = diabetes_ukb,
    hypertension = hypertension_ukb,
    sbp          = rowMeans(cbind(sbp_1, sbp_2), na.rm = TRUE),
    dbp          = NA_real_,
    on_statin    = on_statin_ukb,
    statin_type  = ifelse("statin_name" %in% names(ukb_all), statin_name, NA_character_),
    lpa_143      = ifelse(!is.na(lpa_nmol) & lpa_nmol > 143, 1, 0),
    ascvd        = ascvd,
    ascvd_incident_ukb  = ascvd_incident_ukb,
    ascvd_prevalent_ukb = ascvd_prevalent_ukb,
    followup_years_ukb  = followup_years_ukb,
    gene         = gene_ukb,
    mutation_str = NA_character_,
    mut_type     = NA_character_,
    vus          = NA_real_,
    patient_type = NA_character_,
    corneal_arcus = NA_real_,
    tendon_xanth  = NA_real_,
    xanthelasmas  = NA_real_,
    has_phys_sign = NA_real_,
    simon_broome  = NA_real_,
    dlcn_score    = NA_real_,
    ldl_t2        = NA_real_,
    tc_t2         = NA_real_,
    ldl_pct_red   = NA_real_,
    age_at_event  = NA_real_,
    n_events      = NA_real_
  )

ukb_nonfh <- ukb_all %>%
  filter(is_fh == 0 | is.na(is_fh)) %>%
  mutate(
    cohort = "UKB-nonFH",
    id     = as.character(eid),
    age    = age_raw, sex = sex_num, bmi = bmi_raw,
    tc = tc_raw, hdl = hdl_raw, ldl = ldl_raw, tg = tg_raw,
    apob = apob_raw,
    on_statin = on_statin_ukb,
    ever_smoked = ever_smoked_ukb,
    diabetes = diabetes_ukb,
    hypertension = hypertension_ukb,
    sbp = rowMeans(cbind(sbp_1, sbp_2), na.rm = TRUE)
  )

cat(sprintf("UKB FH harmonised: %d patients\n", nrow(ukb_fh)))
cat(sprintf("UKB non-FH available: %d patients\n", nrow(ukb_nonfh)))

# ── 2d. TRY TO LOAD Lp(a) FOR UKB ───────────────────────────────────────────
cat("\nSearching for Lp(a) data in batch files...\n")
# Check for batch files with Lp(a) data
lpa_files <- list.files(".", pattern = "batch.*\\.csv", full.names = TRUE)
cat(sprintf("  Found %d batch files to scan\n", length(lpa_files)))
if (length(lpa_files) > 0) {
  for (f in lpa_files) {
    tmp <- read.csv(f, nrows = 5, stringsAsFactors = FALSE)
    lpa_col <- grep("p30790|Lpa|lpa", names(tmp), value = TRUE)
    if (length(lpa_col) > 0) {
      cat(sprintf("  Found Lp(a) in %s, col: %s\n", basename(f), paste(lpa_col, collapse=",")))
      batch_lpa <- read.csv(f, stringsAsFactors = FALSE)
      if ("participant.eid" %in% names(batch_lpa)) {
        lpa_map <- batch_lpa %>%
          select(participant.eid, all_of(lpa_col[1])) %>%
          rename(lpa_batch = 2) %>%
          mutate(lpa_batch = safe_numeric(lpa_batch))
        ukb_fh <- ukb_fh %>%
          left_join(lpa_map, by = c("eid" = "participant.eid")) %>%
          mutate(
            lpa_nmol = coalesce(lpa_nmol, lpa_batch),
            lpa_143  = ifelse(!is.na(lpa_nmol) & lpa_nmol > 143, 1, 0)
          ) %>%
          select(-lpa_batch)
      }
      break
    }
  }
}

# Also try calon_ukb_analysis_ready.csv for Lp(a)
if (all(is.na(ukb_fh$lpa_nmol))) {
  if (file.exists("calon_ukb_analysis_ready.csv")) {
    ukb_ready2 <- read.csv("calon_ukb_analysis_ready.csv", stringsAsFactors = FALSE)
    if ("lpa_nmol" %in% names(ukb_ready2) || "Lpa" %in% names(ukb_ready2)) {
      lpa_col_name <- ifelse("lpa_nmol" %in% names(ukb_ready2), "lpa_nmol", "Lpa")
      lpa_map2 <- ukb_ready2 %>%
        select(participant.eid, all_of(lpa_col_name)) %>%
        rename(lpa_ready = 2) %>%
        mutate(lpa_ready = safe_numeric(lpa_ready))
      ukb_fh <- ukb_fh %>%
        left_join(lpa_map2, by = c("eid" = "participant.eid")) %>%
        mutate(
          lpa_nmol = coalesce(lpa_nmol, lpa_ready),
          lpa_143  = ifelse(!is.na(lpa_nmol) & lpa_nmol > 143, 1, 0)
        ) %>%
        select(-lpa_ready)
    }
  }
}

# ── 2e. COMBINE FH COHORTS ──────────────────────────────────────────────────
common_cols <- c("cohort", "id", "age", "sex", "bmi",
                 "tc", "hdl", "ldl", "tg", "apob", "apoa1",
                 "lpa_nmol", "inv_hdl", "apob_ldl", "nhdl", "tg_hdl",
                 "ever_smoked", "diabetes", "hypertension",
                 "sbp", "dbp", "on_statin", "statin_type",
                 "lpa_143", "ascvd",
                 "ascvd_incident_ukb", "ascvd_prevalent_ukb", "followup_years_ukb",
                 "gene", "mutation_str", "mut_type", "vus",
                 "patient_type",
                 "corneal_arcus", "tendon_xanth", "xanthelasmas", "has_phys_sign",
                 "simon_broome", "dlcn_score",
                 "ldl_t2", "tc_t2", "ldl_pct_red",
                 "age_at_event", "n_events")

# Ensure all columns exist
for (cc in common_cols) {
  if (!cc %in% names(d))      d[[cc]]      <- NA
  if (!cc %in% names(w))      w[[cc]]      <- NA
  if (!cc %in% names(ukb_fh)) ukb_fh[[cc]] <- NA
}

fh_pooled <- bind_rows(
  d[, common_cols],
  w[, common_cols],
  ukb_fh[, common_cols]
)

fh_pooled$cohort <- factor(fh_pooled$cohort, levels = c("South Wales FH", "UKB", "Wales"))

cat(sprintf("\n=== POOLED FH COHORT: %d patients ===\n", nrow(fh_pooled)))
cat(sprintf("  South Wales FH: %d\n", sum(fh_pooled$cohort == "South Wales FH")))
cat(sprintf("  UKB:      %d\n", sum(fh_pooled$cohort == "UKB")))
cat(sprintf("  Wales:    %d\n", sum(fh_pooled$cohort == "Wales")))

# ── 2f. DERIVE ASCVD SPLITS, LIPID-YEARS, LOG APOB/LDL ──────────────────────
cat("\n=== DERIVING NEW VARIABLES ===\n")

# --- PRIMARY OUTCOME: Total ASCVD (ascvd) ---
# We use TOTAL ASCVD as the primary outcome for all cohorts because:
#   - Wales FH has NO event dates → cannot split incident/prevalent
#   - Using total ASCVD is consistent across all 3 cohorts
#   - Maximises number of events for model fitting & validation

# --- SENSITIVITY: Incident/Prevalent ASCVD (where available) ---
# UKB: incident/prevalent already computed in 01_CALON_build_cohort.R (via ICD-10 dates)
# SW FH: has age_at_event → can split
# Wales: CANNOT split (no event dates)
fh_pooled <- fh_pooled %>%
  mutate(
    # For UKB: use pre-computed incident/prevalent from processed file
    # For SW FH: use age_at_event vs age logic
    # For Wales: cannot split — set to NA
    ascvd_incident = case_when(
      cohort == "UKB" & !is.na(ascvd_incident_ukb) ~ as.integer(ascvd_incident_ukb),
      cohort == "South Wales FH" & ascvd == 1 & !is.na(age_at_event) & age_at_event >= age ~ 1L,
      cohort == "South Wales FH" & (ascvd == 0 | is.na(ascvd)) ~ 0L,
      cohort == "South Wales FH" & ascvd == 1 & !is.na(age_at_event) & age_at_event < age ~ 0L,
      cohort == "Wales" ~ NA_integer_,  # Cannot determine for Wales
      TRUE ~ NA_integer_
    ),
    ascvd_prevalent = case_when(
      cohort == "UKB" & !is.na(ascvd_prevalent_ukb) ~ as.integer(ascvd_prevalent_ukb),
      cohort == "South Wales FH" & ascvd == 1 & !is.na(age_at_event) & age_at_event < age ~ 1L,
      cohort == "South Wales FH" & (ascvd == 0 | is.na(ascvd)) ~ 0L,
      cohort == "South Wales FH" & ascvd == 1 & !is.na(age_at_event) & age_at_event >= age ~ 0L,
      cohort == "Wales" ~ NA_integer_,
      TRUE ~ NA_integer_
    )
  )

cat(sprintf("  Total ASCVD (primary outcome): %d (%.1f%%)\n",
            sum(fh_pooled$ascvd, na.rm = TRUE),
            100 * mean(fh_pooled$ascvd, na.rm = TRUE)))
cat("  Incident/Prevalent split (sensitivity, excl Wales):\n")
for (coh in c("South Wales FH", "UKB")) {
  sub <- fh_pooled %>% filter(cohort == coh)
  cat(sprintf("    %s: incident=%d (%.1f%%), prevalent=%d (%.1f%%)\n",
              coh,
              sum(sub$ascvd_incident == 1, na.rm = TRUE),
              100 * mean(sub$ascvd_incident == 1, na.rm = TRUE),
              sum(sub$ascvd_prevalent == 1, na.rm = TRUE),
              100 * mean(sub$ascvd_prevalent == 1, na.rm = TRUE)))
}

# --- Cumulative Lipid-Year Exposure (from birth, age 0) ---
# Formula: lipid_years = LDL_pre_treatment × age_at_LLT_start + LDL_on_treatment × years_on_treatment
# For cross-sectional registry data:
#   - If on statin: ldl (pre-treatment) × age is the exposure (treatment starts ~at diagnosis)
#   - If ldl_t2 available: that captures on-treatment LDL, but for registries data is at diagnosis
# Since FH patients have elevated LDL from birth, this is: pre-treatment LDL × current age
# More precisely: we use ldl (baseline, pre-treatment for untreated or first recorded) × age
# For patients with on-treatment LDL (ldl_t2), we still use pre-treatment LDL × age because
# the cumulative exposure happened BEFORE treatment
fh_pooled <- fh_pooled %>%
  mutate(
    # Cumulative lipid-year exposure from birth (mmol/L · years)
    lipid_years = ifelse(!is.na(ldl) & !is.na(age) & ldl > 0 & age > 0,
                         ldl * age, NA_real_)
  )

cat(sprintf("  Lipid-years (N with data): %d\n", sum(!is.na(fh_pooled$lipid_years))))
cat(sprintf("  Lipid-years: mean=%.1f, median=%.1f, IQR=[%.1f, %.1f]\n",
            mean(fh_pooled$lipid_years, na.rm = TRUE),
            median(fh_pooled$lipid_years, na.rm = TRUE),
            quantile(fh_pooled$lipid_years, 0.25, na.rm = TRUE),
            quantile(fh_pooled$lipid_years, 0.75, na.rm = TRUE)))

# --- Log ApoB/LDL Ratio (for CALON-X) ---
# Available: SW FH (ApoB measured), UKB (ApoB from p30640), Wales (NO ApoB)
fh_pooled <- fh_pooled %>%
  mutate(
    log_apob_ldl = ifelse(!is.na(apob_ldl) & apob_ldl > 0,
                          log(apob_ldl), NA_real_)
  )

cat(sprintf("  log(ApoB/LDL) available: %d (%.1f%%)\n",
            sum(!is.na(fh_pooled$log_apob_ldl)),
            100 * mean(!is.na(fh_pooled$log_apob_ldl))))

# Print per-cohort availability
for (coh in levels(fh_pooled$cohort)) {
  sub <- fh_pooled %>% filter(cohort == coh)
  cat(sprintf("\n  %s (N=%d):\n", coh, nrow(sub)))
  cat(sprintf("    ASCVD (total): %.1f%% (N=%d)\n",
              100 * mean(sub$ascvd == 1, na.rm = TRUE),
              sum(sub$ascvd == 1, na.rm = TRUE)))
  cat(sprintf("    Lp(a):         %.1f%% non-missing\n", 100 * mean(!is.na(sub$lpa_143))))
  cat(sprintf("    ever_smoked:   %.1f%% non-missing\n", 100 * mean(!is.na(sub$ever_smoked))))
  cat(sprintf("    lipid_years:   %.1f%% non-missing\n", 100 * mean(!is.na(sub$lipid_years))))
  cat(sprintf("    log_apob_ldl:  %.1f%% non-missing\n", 100 * mean(!is.na(sub$log_apob_ldl))))
  cat(sprintf("    has_phys_sign: %.1f%% non-missing\n", 100 * mean(!is.na(sub$has_phys_sign))))
  if (coh != "Wales") {
    cat(sprintf("    ascvd_incident: %.1f%% (N=%d)\n",
                100 * mean(sub$ascvd_incident == 1, na.rm = TRUE),
                sum(sub$ascvd_incident == 1, na.rm = TRUE)))
    cat(sprintf("    ascvd_prevalent:%.1f%% (N=%d)\n",
                100 * mean(sub$ascvd_prevalent == 1, na.rm = TRUE),
                sum(sub$ascvd_prevalent == 1, na.rm = TRUE)))
  } else {
    cat("    ascvd_incident: N/A (no event dates in Wales FH)\n")
  }
}

################################################################################
# ── 3. TABLE 1: BASELINE CHARACTERISTICS (3-COHORT COMPARISON) ──────────────
################################################################################
cat("\n=== TABLE 1: Baseline Characteristics ===\n")

tab1_vars <- c("age", "sex", "bmi", "tc", "hdl", "ldl", "tg",
               "nhdl", "apob", "apoa1", "lpa_nmol", "tg_hdl",
               "ever_smoked", "diabetes", "hypertension",
               "sbp", "dbp", "on_statin",
               "lpa_143", "ascvd")

tab1_cat <- c("sex", "ever_smoked", "diabetes", "hypertension",
              "on_statin", "lpa_143", "ascvd")

tab1 <- CreateTableOne(
  vars    = tab1_vars,
  strata  = "cohort",
  data    = fh_pooled,
  factorVars = tab1_cat,
  test    = TRUE
)

tab1_print <- print(tab1, showAllLevels = TRUE, nonnormal = c("lpa_nmol", "tg"),
                    smd = TRUE, missing = TRUE, printToggle = FALSE)
write.csv(tab1_print, "Table1_baseline_characteristics.csv")
cat("Table 1 saved.\n")

# ── 3b. TABLE 1b: ASCVD+ vs ASCVD- PER SITE (DESCRIPTIVE STATISTICS) ────────
cat("\n=== TABLE 1b: ASCVD+ vs ASCVD- by Site ===\n")

# All variables for comparison
tab1b_vars <- c("age", "sex", "bmi",
                "tc", "ldl", "hdl", "tg", "nhdl",
                "apob", "apoa1", "apob_ldl",
                "lpa_nmol", "tg_hdl",
                "diabetes", "hypertension", "ever_smoked",
                "on_statin",
                "corneal_arcus", "tendon_xanth", "xanthelasmas", "has_phys_sign",
                "gene", "mut_type", "vus")

tab1b_cat <- c("sex", "diabetes", "hypertension", "ever_smoked",
               "on_statin", "corneal_arcus", "tendon_xanth", "xanthelasmas",
               "has_phys_sign", "gene", "mut_type", "vus")

# Non-normal variables (right-skewed): use median [IQR]
tab1b_nonnorm <- c("lpa_nmol", "tg", "tg_hdl", "apob_ldl")

# Function to generate Table 1 ASCVD+ vs ASCVD- for a given cohort
make_ascvd_table <- function(data, cohort_name) {
  sub <- data %>% filter(cohort == cohort_name, !is.na(ascvd))
  if (nrow(sub) < 10) {
    cat(sprintf("  %s: too few patients, skipping.\n", cohort_name))
    return(NULL)
  }
  sub$ascvd_group <- factor(
    ifelse(sub$ascvd == 1, "ASCVD+", "ASCVD-"),
    levels = c("ASCVD-", "ASCVD+")
  )

  # Only include variables that have at least some non-missing data
  vars_available <- tab1b_vars[sapply(tab1b_vars, function(v) {
    sum(!is.na(sub[[v]])) >= 5
  })]
  cat_available <- intersect(tab1b_cat, vars_available)
  nonnorm_available <- intersect(tab1b_nonnorm, vars_available)

  # Shapiro-Wilk normality test for continuous variables
  cont_vars <- setdiff(vars_available, cat_available)
  nonnorm_final <- nonnorm_available
  for (cv in cont_vars) {
    vals <- sub[[cv]][!is.na(sub[[cv]])]
    if (length(vals) >= 8 && length(vals) <= 5000) {
      sw <- tryCatch(shapiro.test(vals)$p.value, error = function(e) 1)
      if (sw < 0.05 && !(cv %in% nonnorm_final)) {
        nonnorm_final <- c(nonnorm_final, cv)
      }
    }
  }

  tab <- CreateTableOne(
    vars       = vars_available,
    strata     = "ascvd_group",
    data       = sub,
    factorVars = cat_available,
    test       = TRUE
  )
  tab_print <- print(tab, showAllLevels = TRUE, nonnormal = nonnorm_final,
                     exact = cat_available, missing = TRUE, printToggle = FALSE)

  cat(sprintf("\n--- %s: ASCVD+ vs ASCVD- (N=%d) ---\n", cohort_name, nrow(sub)))
  cat(sprintf("  ASCVD+: N=%d (%.1f%%)\n",
              sum(sub$ascvd == 1), 100 * mean(sub$ascvd == 1)))
  cat(sprintf("  ASCVD-: N=%d (%.1f%%)\n",
              sum(sub$ascvd == 0), 100 * mean(sub$ascvd == 0)))

  # Print the full table
  print(tab_print)

  return(tab_print)
}

# Generate per-site ASCVD comparison tables
tab1b_sw   <- make_ascvd_table(fh_pooled, "South Wales FH")
tab1b_ukb  <- make_ascvd_table(fh_pooled, "UKB")
tab1b_wales <- make_ascvd_table(fh_pooled, "Wales")

# Save individual tables
if (!is.null(tab1b_sw))    write.csv(tab1b_sw,    "Table1b_ASCVD_SouthWales.csv")
if (!is.null(tab1b_ukb))   write.csv(tab1b_ukb,   "Table1b_ASCVD_UKB.csv")
if (!is.null(tab1b_wales)) write.csv(tab1b_wales,  "Table1b_ASCVD_Wales.csv")

# Also generate pooled ASCVD+ vs ASCVD- across all sites
cat("\n--- POOLED: ASCVD+ vs ASCVD- ---\n")
pooled_ascvd <- fh_pooled %>% filter(!is.na(ascvd))
pooled_ascvd$ascvd_group <- factor(
  ifelse(pooled_ascvd$ascvd == 1, "ASCVD+", "ASCVD-"),
  levels = c("ASCVD-", "ASCVD+")
)

vars_pool <- tab1b_vars[sapply(tab1b_vars, function(v) {
  sum(!is.na(pooled_ascvd[[v]])) >= 5
})]
cat_pool <- intersect(tab1b_cat, vars_pool)

tab1b_pooled <- CreateTableOne(
  vars       = vars_pool,
  strata     = "ascvd_group",
  data       = pooled_ascvd,
  factorVars = cat_pool,
  test       = TRUE
)
tab1b_pooled_print <- print(tab1b_pooled, showAllLevels = TRUE,
                            nonnormal = intersect(tab1b_nonnorm, vars_pool),
                            missing = TRUE, printToggle = FALSE)
print(tab1b_pooled_print)
write.csv(tab1b_pooled_print, "Table1b_ASCVD_Pooled.csv")
cat("Table 1b (ASCVD comparison) saved.\n")

################################################################################
# ── 4. TABLE 2: GENETIC ARCHITECTURE ────────────────────────────────────────
################################################################################
cat("\n=== TABLE 2: Genetic Architecture ===\n")

# Gene distribution
gene_tab <- fh_pooled %>%
  filter(!is.na(gene)) %>%
  group_by(cohort, gene) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cohort) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  pivot_wider(names_from = cohort, values_from = c(n, pct))

# Mutation type distribution
mut_tab <- fh_pooled %>%
  filter(!is.na(mut_type), gene == "LDLR") %>%
  group_by(cohort, mut_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cohort) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  pivot_wider(names_from = cohort, values_from = c(n, pct))

# VUS
vus_tab <- fh_pooled %>%
  filter(!is.na(vus)) %>%
  group_by(cohort) %>%
  summarise(
    total = n(),
    n_vus = sum(vus == 1, na.rm = TRUE),
    pct_vus = round(100 * n_vus / total, 1),
    .groups = "drop"
  )

# ASCVD by gene
ascvd_by_gene <- fh_pooled %>%
  filter(!is.na(gene), !is.na(ascvd)) %>%
  group_by(cohort, gene) %>%
  summarise(
    n = n(),
    events = sum(ascvd == 1, na.rm = TRUE),
    rate   = round(100 * events / n, 1),
    .groups = "drop"
  )

cat("Gene distribution:\n")
print(gene_tab)
cat("\nMutation type (LDLR only):\n")
print(mut_tab)
cat("\nVUS:\n")
print(vus_tab)
cat("\nASCVD by gene:\n")
print(ascvd_by_gene)

# Save
write.csv(gene_tab, "Table2a_gene_distribution.csv", row.names = FALSE)
write.csv(ascvd_by_gene, "Table2b_ascvd_by_gene.csv", row.names = FALSE)

# ── 4b. GENETIC VARIANTS & VUS: ASCVD+ vs ASCVD- COMPARISON ─────────────────
cat("\n=== TABLE 2c: Genetic Variants ASCVD+ vs ASCVD- ===\n")

# Gene distribution ASCVD+ vs ASCVD- by cohort
gene_ascvd <- fh_pooled %>%
  filter(!is.na(gene), !is.na(ascvd)) %>%
  mutate(ascvd_group = ifelse(ascvd == 1, "ASCVD+", "ASCVD-")) %>%
  group_by(cohort, ascvd_group, gene) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cohort, ascvd_group) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

cat("\nGene distribution by ASCVD status:\n")
for (coh in c("South Wales FH", "UKB", "Wales")) {
  sub <- gene_ascvd %>% filter(cohort == coh)
  if (nrow(sub) > 0) {
    cat(sprintf("\n  -- %s --\n", coh))
    sub_wide <- sub %>%
      pivot_wider(names_from = ascvd_group, values_from = c(n, pct),
                  values_fill = 0)
    print(sub_wide)
  }
}

# Chi-square: gene vs ASCVD by cohort
cat("\nChi-square: gene vs ASCVD by cohort:\n")
for (coh in c("South Wales FH", "UKB", "Wales")) {
  sub <- fh_pooled %>%
    filter(cohort == coh, !is.na(gene), !is.na(ascvd))
  if (nrow(sub) > 10) {
    ct <- table(sub$gene, sub$ascvd)
    if (nrow(ct) > 1 && ncol(ct) > 1) {
      test <- tryCatch(fisher.test(ct, simulate.p.value = TRUE, B = 10000),
                       error = function(e) NULL)
      if (!is.null(test)) {
        cat(sprintf("  %s: Fisher's exact p = %.4e\n", coh, test$p.value))
      }
    }
  }
}

# ORs for each gene vs LDLR (reference) per cohort
cat("\nASCVD odds ratios by gene (ref: LDLR):\n")
for (coh in c("South Wales FH", "UKB", "Wales")) {
  sub <- fh_pooled %>%
    filter(cohort == coh, !is.na(gene), !is.na(ascvd), gene %in% c("LDLR", "APOB", "PCSK9", "Other"))
  if (nrow(sub) > 10 && sum(sub$ascvd == 1) > 5) {
    sub$gene_f <- relevel(factor(sub$gene), ref = "LDLR")
    mod <- tryCatch(glm(ascvd ~ gene_f, data = sub, family = binomial),
                    error = function(e) NULL)
    if (!is.null(mod)) {
      cat(sprintf("\n  -- %s (N=%d, events=%d) --\n", coh, nrow(sub), sum(sub$ascvd)))
      or_tab <- tidy(mod, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(term != "(Intercept)") %>%
        mutate(
          term = gsub("gene_f", "", term),
          sig = case_when(p.value < 0.001 ~ "***",
                          p.value < 0.01 ~ "**",
                          p.value < 0.05 ~ "*",
                          TRUE ~ "")
        )
      for (i in seq_len(nrow(or_tab))) {
        cat(sprintf("    %s vs LDLR: OR=%.2f (95%% CI: %.2f-%.2f), p=%.4e %s\n",
                    or_tab$term[i], or_tab$estimate[i],
                    or_tab$conf.low[i], or_tab$conf.high[i],
                    or_tab$p.value[i], or_tab$sig[i]))
      }
    }
  }
}

# Mutation type (null vs defective) and ASCVD
cat("\nMutation type (Null vs Defective LDLR) & ASCVD:\n")
for (coh in c("South Wales FH", "Wales")) {
  sub <- fh_pooled %>%
    filter(cohort == coh, gene == "LDLR", !is.na(mut_type), !is.na(ascvd))
  if (nrow(sub) > 10) {
    null_s <- sub %>% filter(mut_type == "Null")
    def_s  <- sub %>% filter(mut_type == "Defective")
    null_rate <- 100 * mean(null_s$ascvd == 1)
    def_rate  <- 100 * mean(def_s$ascvd == 1)

    # Fisher's exact test
    ct <- table(sub$mut_type, sub$ascvd)
    ft <- tryCatch(fisher.test(ct), error = function(e) NULL)
    p_val <- if (!is.null(ft)) ft$p.value else NA
    or_val <- if (!is.null(ft)) ft$estimate else NA

    cat(sprintf("\n  -- %s (LDLR only, N=%d) --\n", coh, nrow(sub)))
    cat(sprintf("    Null: N=%d, ASCVD rate=%.1f%%\n", nrow(null_s), null_rate))
    cat(sprintf("    Defective: N=%d, ASCVD rate=%.1f%%\n", nrow(def_s), def_rate))
    cat(sprintf("    OR (Null vs Defective): %.2f, Fisher p=%.4e\n", or_val, p_val))
  }
}

# VUS vs pathogenic: ASCVD comparison
cat("\nVUS vs Pathogenic variant & ASCVD:\n")
for (coh in c("South Wales FH", "Wales")) {
  sub <- fh_pooled %>%
    filter(cohort == coh, !is.na(vus), !is.na(ascvd))
  if (nrow(sub) > 10) {
    vus_pos <- sub %>% filter(vus == 1)
    vus_neg <- sub %>% filter(vus == 0)
    if (nrow(vus_pos) >= 3 && nrow(vus_neg) >= 3) {
      vus_rate <- 100 * mean(vus_pos$ascvd == 1)
      path_rate <- 100 * mean(vus_neg$ascvd == 1)

      ct <- table(sub$vus, sub$ascvd)
      ft <- tryCatch(fisher.test(ct), error = function(e) NULL)
      p_val <- if (!is.null(ft)) ft$p.value else NA
      or_val <- if (!is.null(ft)) ft$estimate else NA

      cat(sprintf("\n  -- %s --\n", coh))
      cat(sprintf("    VUS carriers: N=%d, ASCVD rate=%.1f%%\n", nrow(vus_pos), vus_rate))
      cat(sprintf("    Pathogenic:   N=%d, ASCVD rate=%.1f%%\n", nrow(vus_neg), path_rate))
      cat(sprintf("    OR (VUS vs Pathogenic): %.2f, Fisher p=%.4e\n", or_val, p_val))

      # Baseline comparison VUS vs pathogenic
      vus_comp_vars <- c("age", "sex", "ldl", "tc", "hdl", "tg",
                         "diabetes", "hypertension", "on_statin", "ascvd")
      vus_comp_cat <- c("sex", "diabetes", "hypertension", "on_statin", "ascvd")

      # Only include available variables
      vus_vars_avail <- vus_comp_vars[sapply(vus_comp_vars, function(v) {
        sum(!is.na(sub[[v]])) >= 5
      })]
      vus_cat_avail <- intersect(vus_comp_cat, vus_vars_avail)

      sub$vus_label <- factor(ifelse(sub$vus == 1, "VUS", "Pathogenic"),
                              levels = c("Pathogenic", "VUS"))
      tab_vus <- CreateTableOne(
        vars       = vus_vars_avail,
        strata     = "vus_label",
        data       = sub,
        factorVars = vus_cat_avail,
        test       = TRUE
      )
      tab_vus_print <- print(tab_vus, showAllLevels = TRUE, printToggle = FALSE)
      cat(sprintf("\n  Baseline comparison VUS vs Pathogenic (%s):\n", coh))
      print(tab_vus_print)
      write.csv(tab_vus_print,
                sprintf("Table2c_VUS_vs_Pathogenic_%s.csv", gsub(" ", "_", coh)))
    }
  }
}

# ASCVD-adjusted model: gene + age + sex
cat("\nMultivariable logistic regression: ASCVD ~ gene + age + sex\n")
for (coh in c("South Wales FH", "Wales")) {
  sub <- fh_pooled %>%
    filter(cohort == coh, !is.na(gene), !is.na(ascvd), !is.na(age), !is.na(sex))
  if (nrow(sub) > 30) {
    sub$gene_f <- relevel(factor(sub$gene), ref = "LDLR")
    mod <- tryCatch(glm(ascvd ~ gene_f + age + sex, data = sub, family = binomial),
                    error = function(e) NULL)
    if (!is.null(mod)) {
      cat(sprintf("\n  -- %s (N=%d) --\n", coh, nrow(sub)))
      or_tab <- tidy(mod, conf.int = TRUE, exponentiate = TRUE) %>%
        mutate(sig = case_when(p.value < 0.001 ~ "***",
                               p.value < 0.01 ~ "**",
                               p.value < 0.05 ~ "*",
                               TRUE ~ ""))
      for (i in seq_len(nrow(or_tab))) {
        cat(sprintf("    %s: OR=%.2f (95%% CI: %.2f-%.2f), p=%.4e %s\n",
                    gsub("gene_f", "", or_tab$term[i]),
                    or_tab$estimate[i], or_tab$conf.low[i], or_tab$conf.high[i],
                    or_tab$p.value[i], or_tab$sig[i]))
      }
    }
  }
}

write.csv(gene_ascvd, "Table2c_gene_by_ascvd.csv", row.names = FALSE)
cat("\nTable 2c (Genetic variants & ASCVD comparison) saved.\n")

################################################################################
# ── 5. INDEX vs CASCADE ANALYSIS ────────────────────────────────────────────
################################################################################
cat("\n=== INDEX vs CASCADE ANALYSIS ===\n")

idx_casc <- fh_pooled %>%
  filter(!is.na(patient_type))

# Baseline comparison
idx_casc_vars <- c("age", "sex", "ldl", "hdl", "tc", "tg", "bmi",
                   "lpa_nmol", "ever_smoked", "diabetes", "hypertension",
                   "on_statin", "ascvd")
idx_casc_cat <- c("sex", "ever_smoked", "diabetes", "hypertension",
                  "on_statin", "ascvd")

tab_idx <- CreateTableOne(
  vars       = idx_casc_vars,
  strata     = "patient_type",
  data       = idx_casc,
  factorVars = idx_casc_cat,
  test       = TRUE
)
tab_idx_print <- print(tab_idx, showAllLevels = TRUE, printToggle = FALSE)
write.csv(tab_idx_print, "Table3_index_vs_cascade.csv")

# By cohort
for (coh in c("South Wales FH", "Wales")) {
  sub <- idx_casc %>% filter(cohort == coh)
  if (nrow(sub) > 10) {
    cat(sprintf("\n-- %s: Index vs Cascade --\n", coh))
    cat(sprintf("  Index:   N=%d, ASCVD=%.1f%%\n",
                sum(sub$patient_type == "Index"),
                100 * mean(sub$ascvd[sub$patient_type == "Index"], na.rm = TRUE)))
    cat(sprintf("  Cascade: N=%d, ASCVD=%.1f%%\n",
                sum(sub$patient_type == "Cascade"),
                100 * mean(sub$ascvd[sub$patient_type == "Cascade"], na.rm = TRUE)))
  }
}

################################################################################
# ── 6. NOVEL BIOMARKER ANALYSES ─────────────────────────────────────────────
################################################################################
cat("\n=== NOVEL BIOMARKER ANALYSES ===\n")

# ── 6a. ApoB Distribution and ApoB/LDL Ratio in FH ─────────────────────────
apob_data <- fh_pooled %>% filter(!is.na(apob))
cat(sprintf("ApoB available: Dragon=%d, UKB=%d, Wales=%d\n",
            sum(apob_data$cohort == "South Wales FH"),
            sum(apob_data$cohort == "UKB"),
            sum(apob_data$cohort == "Wales")))

# ApoB by ASCVD status
apob_ascvd <- apob_data %>%
  filter(!is.na(ascvd)) %>%
  group_by(cohort, ascvd) %>%
  summarise(
    n      = n(),
    mean   = round(mean(apob, na.rm = TRUE), 3),
    sd     = round(sd(apob, na.rm = TRUE), 3),
    median = round(median(apob, na.rm = TRUE), 3),
    .groups = "drop"
  )
cat("\nApoB by ASCVD status:\n")
print(apob_ascvd)

# ApoB/LDL ratio
apob_ldl_data <- fh_pooled %>%
  filter(!is.na(apob_ldl), apob_ldl > 0, apob_ldl < 1)

if (nrow(apob_ldl_data) > 0) {
  cat(sprintf("\nApoB/LDL ratio: N=%d, Mean=%.3f, SD=%.3f\n",
              nrow(apob_ldl_data),
              mean(apob_ldl_data$apob_ldl, na.rm = TRUE),
              sd(apob_ldl_data$apob_ldl, na.rm = TRUE)))

  # ApoB/LDL by ASCVD
  apob_ldl_ascvd <- apob_ldl_data %>%
    filter(!is.na(ascvd)) %>%
    group_by(ascvd) %>%
    summarise(
      n      = n(),
      mean   = round(mean(apob_ldl, na.rm = TRUE), 4),
      sd     = round(sd(apob_ldl, na.rm = TRUE), 4),
      .groups = "drop"
    )
  cat("ApoB/LDL ratio by ASCVD:\n")
  print(apob_ldl_ascvd)
}

# ── 6b. LDL Regression: ApoB vs LDL in FH vs Non-FH ────────────────────────
cat("\n--- ApoB vs LDL Regression ---\n")

# FH UKB
fh_apob_ldl <- ukb_fh %>%
  filter(!is.na(apob), !is.na(ldl), ldl > 0, apob > 0)

# Non-FH UKB (sample for speed)
set.seed(42)
nonfh_apob_ldl <- ukb_nonfh %>%
  filter(!is.na(apob_raw), !is.na(ldl_raw), ldl_raw > 0, apob_raw > 0) %>%
  sample_n(min(10000, nrow(.)))

if (nrow(fh_apob_ldl) > 30 & nrow(nonfh_apob_ldl) > 30) {
  lm_fh    <- lm(apob ~ ldl, data = fh_apob_ldl)
  lm_nonfh <- lm(apob_raw ~ ldl_raw, data = nonfh_apob_ldl)

  cat(sprintf("FH:     ApoB = %.3f + %.3f * LDL, R²=%.3f, N=%d\n",
              coef(lm_fh)[1], coef(lm_fh)[2], summary(lm_fh)$r.squared, nrow(fh_apob_ldl)))
  cat(sprintf("Non-FH: ApoB = %.3f + %.3f * LDL, R²=%.3f, N=%d\n",
              coef(lm_nonfh)[1], coef(lm_nonfh)[2], summary(lm_nonfh)$r.squared, nrow(nonfh_apob_ldl)))
}

# ── 6c. TG/HDL Ratio ────────────────────────────────────────────────────────
tg_hdl_summary <- fh_pooled %>%
  filter(!is.na(tg_hdl)) %>%
  group_by(cohort) %>%
  summarise(
    n      = n(),
    mean   = round(mean(tg_hdl, na.rm = TRUE), 2),
    median = round(median(tg_hdl, na.rm = TRUE), 2),
    .groups = "drop"
  )
cat("\nTG/HDL ratio by cohort:\n")
print(tg_hdl_summary)

################################################################################
# ── 7. Lp(a) DISTRIBUTION & THRESHOLD OPTIMISATION ─────────────────────────
################################################################################
cat("\n=== Lp(a) ANALYSIS ===\n")

lpa_data <- fh_pooled %>% filter(!is.na(lpa_nmol))
cat(sprintf("Lp(a) available: Dragon=%d, UKB=%d, Wales=%d, Total=%d\n",
            sum(lpa_data$cohort == "South Wales FH"),
            sum(lpa_data$cohort == "UKB"),
            sum(lpa_data$cohort == "Wales"),
            nrow(lpa_data)))

# Distribution stats
lpa_stats <- lpa_data %>%
  group_by(cohort) %>%
  summarise(
    n      = n(),
    mean   = round(mean(lpa_nmol), 1),
    median = round(median(lpa_nmol), 1),
    p25    = round(quantile(lpa_nmol, 0.25), 1),
    p75    = round(quantile(lpa_nmol, 0.75), 1),
    p90    = round(quantile(lpa_nmol, 0.90), 1),
    pct_above_143 = round(100 * mean(lpa_nmol > 143), 1),
    .groups = "drop"
  )
cat("Lp(a) distribution by cohort:\n")
print(lpa_stats)

# Threshold scan (ASCVD discrimination)
if (nrow(lpa_data %>% filter(!is.na(ascvd))) > 100) {
  thresholds <- seq(50, 250, by = 10)
  thresh_results <- data.frame()

  for (th in thresholds) {
    for (coh in c("South Wales FH", "UKB", "Wales")) {
      sub <- lpa_data %>% filter(cohort == coh, !is.na(ascvd))
      if (nrow(sub) > 30 & length(unique(sub$ascvd)) == 2) {
        sub$lpa_bin <- ifelse(sub$lpa_nmol > th, 1, 0)
        if (length(unique(sub$lpa_bin)) == 2) {
          roc_obj <- roc(sub$ascvd, sub$lpa_bin, quiet = TRUE)
          thresh_results <- rbind(thresh_results, data.frame(
            cohort    = coh,
            threshold = th,
            auc       = as.numeric(auc(roc_obj)),
            n_above   = sum(sub$lpa_bin == 1),
            pct_above = round(100 * mean(sub$lpa_bin), 1)
          ))
        }
      }
    }
  }

  if (nrow(thresh_results) > 0) {
    cat("\nLp(a) threshold scan (AUC for ASCVD by binary Lp(a)):\n")
    best_thresh <- thresh_results %>%
      group_by(cohort) %>%
      slice_max(auc, n = 1)
    print(best_thresh)
  }
}

################################################################################
# ── 8. FH vs NON-FH MATCHING IN UKB ────────────────────────────────────────
################################################################################
cat("\n=== FH vs NON-FH PROPENSITY SCORE MATCHING ===\n")

# Prepare matching data
match_fh <- ukb_fh %>%
  filter(!is.na(age), !is.na(sex), !is.na(ldl)) %>%
  mutate(fh_status = 1) %>%
  select(id = eid, fh_status, age, sex, bmi, tc, hdl, ldl, tg,
         apob, on_statin, ever_smoked, diabetes, hypertension, sbp, ascvd)

match_nonfh <- ukb_nonfh %>%
  filter(!is.na(age_raw), !is.na(sex_num), !is.na(ldl_raw)) %>%
  mutate(fh_status = 0) %>%
  select(id = eid, fh_status, age = age_raw, sex = sex_num,
         bmi = bmi_raw, tc = tc_raw, hdl = hdl_raw, ldl = ldl_raw,
         tg = tg_raw, apob = apob_raw, on_statin = on_statin_ukb,
         ever_smoked = ever_smoked_ukb, diabetes = diabetes_ukb,
         hypertension = hypertension_ukb,
         sbp = sbp_1, ascvd)

match_data <- bind_rows(match_fh, match_nonfh) %>%
  filter(complete.cases(age, sex, bmi))

cat(sprintf("Matching pool: FH=%d, non-FH=%d\n",
            sum(match_data$fh_status == 1), sum(match_data$fh_status == 0)))

# 1:3 Propensity Score Matching on age, sex, BMI
if (nrow(match_data) > 1000) {
  m.out <- matchit(fh_status ~ age + sex + bmi,
                   data    = match_data,
                   method  = "nearest",
                   ratio   = 3,
                   caliper = 0.2)

  matched <- match.data(m.out)
  cat(sprintf("Matched: FH=%d, non-FH=%d\n",
              sum(matched$fh_status == 1), sum(matched$fh_status == 0)))

  # Compare lipids after matching
  matched_tab_vars <- c("ldl", "hdl", "tc", "tg", "apob", "on_statin",
                        "ever_smoked", "diabetes", "hypertension", "ascvd")
  matched_cat <- c("on_statin", "ever_smoked", "diabetes", "hypertension", "ascvd")

  tab_matched <- CreateTableOne(
    vars       = matched_tab_vars,
    strata     = "fh_status",
    data       = matched,
    factorVars = matched_cat,
    test       = TRUE
  )
  tab_matched_print <- print(tab_matched, showAllLevels = TRUE, smd = TRUE,
                             printToggle = FALSE)
  write.csv(tab_matched_print, "Table4_FH_vs_NonFH_matched.csv")
  cat("Matched comparison table saved.\n")

  # ApoB vs LDL in matched groups
  matched_fh    <- matched %>% filter(fh_status == 1, !is.na(apob), !is.na(ldl))
  matched_nonfh <- matched %>% filter(fh_status == 0, !is.na(apob), !is.na(ldl))
  if (nrow(matched_fh) > 20 & nrow(matched_nonfh) > 20) {
    cat(sprintf("\nMatched ApoB/LDL: FH mean=%.3f, non-FH mean=%.3f\n",
                mean(matched_fh$apob / matched_fh$ldl, na.rm = TRUE),
                mean(matched_nonfh$apob / matched_nonfh$ldl, na.rm = TRUE)))
  }
}

################################################################################
# ── 9. LDL RESPONSE TO LLT ─────────────────────────────────────────────────
################################################################################
cat("\n=== LDL RESPONSE TO LLT ===\n")

# South Wales FH: has serial lipid panels and LDL % reduction
llt_dragon <- d %>%
  filter(!is.na(ldl), !is.na(ldl_t2), ldl > 0, ldl_t2 > 0) %>%
  mutate(ldl_change_pct = 100 * (ldl_t2 - ldl) / ldl)

if (nrow(llt_dragon) > 10) {
  cat(sprintf("South Wales FH LDL change: N=%d, Mean=%.1f%%, Median=%.1f%%\n",
              nrow(llt_dragon),
              mean(llt_dragon$ldl_change_pct, na.rm = TRUE),
              median(llt_dragon$ldl_change_pct, na.rm = TRUE)))

  # By gene
  llt_by_gene <- llt_dragon %>%
    filter(!is.na(gene)) %>%
    group_by(gene) %>%
    summarise(
      n = n(),
      mean_pct = round(mean(ldl_change_pct, na.rm = TRUE), 1),
      .groups = "drop"
    )
  cat("LDL change by gene:\n")
  print(llt_by_gene)
}

# Wales: same
llt_wales <- w %>%
  filter(!is.na(ldl), !is.na(ldl_t2), ldl > 0, ldl_t2 > 0) %>%
  mutate(ldl_change_pct = 100 * (ldl_t2 - ldl) / ldl)

if (nrow(llt_wales) > 10) {
  cat(sprintf("\nWales LDL change: N=%d, Mean=%.1f%%, Median=%.1f%%\n",
              nrow(llt_wales),
              mean(llt_wales$ldl_change_pct, na.rm = TRUE),
              median(llt_wales$ldl_change_pct, na.rm = TRUE)))
}

################################################################################
# ── 10. BMI IN FH ───────────────────────────────────────────────────────────
################################################################################
cat("\n=== BMI IN FH ===\n")

bmi_data <- fh_pooled %>% filter(!is.na(bmi), bmi > 10, bmi < 60)
bmi_summary <- bmi_data %>%
  group_by(cohort) %>%
  summarise(
    n      = n(),
    mean   = round(mean(bmi), 1),
    sd     = round(sd(bmi), 1),
    median = round(median(bmi), 1),
    pct_overweight = round(100 * mean(bmi >= 25), 1),
    pct_obese      = round(100 * mean(bmi >= 30), 1),
    .groups = "drop"
  )
cat("BMI by cohort:\n")
print(bmi_summary)

# BMI vs ASCVD
bmi_ascvd <- bmi_data %>%
  filter(!is.na(ascvd)) %>%
  mutate(bmi_cat = cut(bmi, breaks = c(0, 25, 30, 60),
                       labels = c("<25", "25-30", ">30"))) %>%
  group_by(cohort, bmi_cat) %>%
  summarise(
    n = n(),
    events = sum(ascvd, na.rm = TRUE),
    rate   = round(100 * events / n, 1),
    .groups = "drop"
  )
cat("\nASCVD event rate by BMI category:\n")
print(bmi_ascvd)

################################################################################
# ── 11. TWO-TIER MODEL: CALON-LITE & CALON-X ────────────────────────────────
################################################################################
cat("\n=== TWO-TIER MODEL ARCHITECTURE ===\n")

# --- CALON-Lite: Universal 6-predictor model ---
# Uses cumulative lipid-years (replaces separate age + LDL — captures both in one
# biologically meaningful predictor), sex, diabetes, ever_smoked, lpa_143, inv_hdl
LITE_features <- c("lipid_years", "sex", "diabetes", "ever_smoked",
                    "lpa_143", "inv_hdl")

# --- CALON-Lite-noLpa: 5-predictor variant for clinics without Lp(a) measurement ---
# Identical to CALON-Lite minus lpa_143 — for settings where Lp(a) is not routinely measured
LITE_NOLPA_features <- c("lipid_years", "sex", "diabetes", "ever_smoked", "inv_hdl")

# --- CALON-X: Enhanced model (Lite + additional predictors where available) ---
# Adds: log_apob_ldl (ApoB/LDL ratio) + has_phys_sign (clinical stigmata)
# NOTE: prior_ascvd DROPPED — causes quasi-separation in SW FH (0.2% prevalence)
X_features    <- c("lipid_years", "sex", "diabetes", "ever_smoked",
                    "lpa_143", "inv_hdl",
                    "log_apob_ldl", "has_phys_sign")

# Legacy comparison: old CALON-X2 (8 predictors) and SAFEHEART baseline
X2_features   <- c("age", "sex", "ldl", "inv_hdl", "ever_smoked",
                    "diabetes", "hypertension", "lpa_143")
SAFE_features <- c("age", "sex", "ldl", "inv_hdl", "ever_smoked",
                    "diabetes", "hypertension")

cat("CALON-Lite (6 pred):      ", paste(LITE_features, collapse=", "), "\n")
cat("CALON-Lite-noLpa (5 pred):", paste(LITE_NOLPA_features, collapse=", "), "\n")
cat("CALON-X (8 pred):         ", paste(X_features, collapse=", "), "\n")
cat("Legacy CALON-X2 (8 pred): ", paste(X2_features, collapse=", "), "\n")
cat("SAFEHEART baseline (7 pred):", paste(SAFE_features, collapse=", "), "\n")

# Prepare cohort-specific data (outcome = ascvd — TOTAL ASCVD as primary)
prep_cohort <- function(df, features, outcome = "ascvd") {
  cols <- c(features, outcome)
  # Ensure all columns exist
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    cat(sprintf("  WARNING: Missing columns: %s\n", paste(missing_cols, collapse=", ")))
    return(data.frame())
  }
  out <- df %>%
    select(all_of(cols)) %>%
    filter(complete.cases(.))
  return(out)
}

# AUC with CI via bootstrap
get_auc_ci <- function(train_df, test_df, features, outcome = "ascvd", n_boot = 1000) {
  train <- prep_cohort(train_df, features, outcome)
  test  <- prep_cohort(test_df, features, outcome)

  if (nrow(train) < 30 || nrow(test) < 30 || length(unique(test[[outcome]])) < 2) {
    return(list(auc = NA, ci_low = NA, ci_high = NA, n_train = nrow(train), n_test = nrow(test)))
  }

  # Scale & fit
  X_train <- scale(as.matrix(train[, features]))
  y_train <- train[[outcome]]
  X_test  <- scale(as.matrix(test[, features]),
                   center = attr(X_train, "scaled:center"),
                   scale  = attr(X_train, "scaled:scale"))
  y_test  <- test[[outcome]]

  fit <- glm(y ~ ., data = data.frame(y = y_train, X_train), family = binomial)
  pred <- predict(fit, newdata = data.frame(X_test), type = "response")

  roc_obj <- roc(y_test, pred, quiet = TRUE)
  ci_obj  <- ci.auc(roc_obj, method = "bootstrap", boot.n = n_boot, quiet = TRUE)

  return(list(
    auc     = as.numeric(auc(roc_obj)),
    ci_low  = ci_obj[1],
    ci_high = ci_obj[3],
    n_train = nrow(train),
    n_test  = nrow(test),
    roc     = roc_obj,
    pred    = pred,
    y_test  = y_test
  ))
}

# Run bi-external validation (outcome = total ASCVD)
# NOTE: CALON-Lite-noLpa used for directions involving Wales (Lp(a) only 2.9% available)
directions <- list(
  c("South Wales FH", "UKB",   "SW \u2192 UKB"),
  c("UKB",      "Wales", "UKB \u2192 W"),
  c("South Wales FH", "Wales", "SW \u2192 W")
)

results <- data.frame()

for (dir in directions) {
  train_name <- dir[1]
  test_name  <- dir[2]
  label      <- dir[3]

  train_df <- fh_pooled %>% filter(cohort == train_name)
  test_df  <- fh_pooled %>% filter(cohort == test_name)

  cat(sprintf("\n--- Direction: %s (train N=%d, test N=%d) ---\n",
              label, nrow(train_df), nrow(test_df)))

  # Determine if Wales is involved (use noLpa variants for those directions)
  involves_wales <- test_name == "Wales" | train_name == "Wales"

  # 1. CALON-Lite (6 predictors) — may fail if Lp(a) coverage too low
  res_lite <- get_auc_ci(train_df, test_df, LITE_features, outcome = "ascvd", n_boot = 2000)
  cat(sprintf("  CALON-Lite:      AUC=%.4f, N_train=%d, N_test=%d\n",
              ifelse(is.na(res_lite$auc), -1, res_lite$auc), res_lite$n_train, res_lite$n_test))

  # 2. CALON-Lite-noLpa (5 predictors) — robust for all directions
  res_lite5 <- get_auc_ci(train_df, test_df, LITE_NOLPA_features, outcome = "ascvd", n_boot = 2000)
  cat(sprintf("  CALON-Lite-noLpa:AUC=%.4f, N_train=%d, N_test=%d\n",
              ifelse(is.na(res_lite5$auc), -1, res_lite5$auc), res_lite5$n_train, res_lite5$n_test))

  # 3. CALON-X (enhanced: Lite + ApoB/LDL + stigmata)
  res_x <- get_auc_ci(train_df, test_df, X_features, outcome = "ascvd", n_boot = 2000)
  cat(sprintf("  CALON-X:         AUC=%.4f, N_train=%d, N_test=%d\n",
              ifelse(is.na(res_x$auc), -1, res_x$auc), res_x$n_train, res_x$n_test))

  # 4. Legacy CALON-X2 (8 predictors)
  res_x2 <- get_auc_ci(train_df, test_df, X2_features, outcome = "ascvd", n_boot = 2000)
  cat(sprintf("  CALON-X2:        AUC=%.4f, N_train=%d, N_test=%d\n",
              ifelse(is.na(res_x2$auc), -1, res_x2$auc), res_x2$n_train, res_x2$n_test))

  # 5. SAFEHEART baseline (7 predictors — no Lp(a))
  res_safe <- get_auc_ci(train_df, test_df, SAFE_features, outcome = "ascvd", n_boot = 2000)
  cat(sprintf("  SAFEHEART:       AUC=%.4f, N_train=%d, N_test=%d\n",
              ifelse(is.na(res_safe$auc), -1, res_safe$auc), res_safe$n_train, res_safe$n_test))

  results <- rbind(results, data.frame(
    direction = label,
    model     = c("CALON-Lite", "CALON-Lite-noLpa", "CALON-X", "CALON-X2", "SAFEHEART"),
    auc       = c(res_lite$auc, res_lite5$auc, res_x$auc, res_x2$auc, res_safe$auc),
    ci_low    = c(res_lite$ci_low, res_lite5$ci_low, res_x$ci_low, res_x2$ci_low, res_safe$ci_low),
    ci_high   = c(res_lite$ci_high, res_lite5$ci_high, res_x$ci_high, res_x2$ci_high, res_safe$ci_high),
    n_train   = c(res_lite$n_train, res_lite5$n_train, res_x$n_train, res_x2$n_train, res_safe$n_train),
    n_test    = c(res_lite$n_test, res_lite5$n_test, res_x$n_test, res_x2$n_test, res_safe$n_test),
    delta     = c(NA, NA, NA, NA, NA)
  ))
}

# Compute delta vs SAFEHEART
results <- results %>%
  group_by(direction) %>%
  mutate(delta = auc - auc[model == "SAFEHEART"]) %>%
  ungroup()

cat("\n=== BI-EXTERNAL VALIDATION RESULTS (Outcome: Total ASCVD) ===\n")
print(results %>% select(direction, model, auc, ci_low, ci_high, delta, n_test))

write.csv(results, "Table5_biexternal_validation.csv", row.names = FALSE)

################################################################################
# ── 12. SUBGROUP ANALYSIS ──────────────────────────────────────────────────
################################################################################
cat("\n=== SUBGROUP ANALYSIS ===\n")

run_subgroup <- function(train_df, test_df, subgroup_col, subgroup_label,
                         direction_label) {
  sub_results <- data.frame()

  test_all <- test_df %>% filter(!is.na(.data[[subgroup_col]]))
  if (nrow(test_all) < 30) return(sub_results)

  for (val in sort(unique(test_all[[subgroup_col]]))) {
    test_sub <- test_all %>% filter(.data[[subgroup_col]] == val)

    if (nrow(test_sub) < 15 || length(unique(test_sub$ascvd)) < 2) next

    auc_lite <- tryCatch({
      r <- get_auc_ci(train_df, test_sub, LITE_features, n_boot = 500)
      r$auc
    }, error = function(e) NA)

    auc_safe <- tryCatch({
      r <- get_auc_ci(train_df, test_sub, SAFE_features, n_boot = 500)
      r$auc
    }, error = function(e) NA)

    if (!is.na(auc_lite) & !is.na(auc_safe)) {
      sub_results <- rbind(sub_results, data.frame(
        direction = direction_label,
        subgroup  = subgroup_label,
        value     = as.character(val),
        n         = nrow(test_sub),
        events    = sum(test_sub$ascvd, na.rm = TRUE),
        auc_lite  = round(auc_lite, 4),
        auc_safe  = round(auc_safe, 4),
        delta     = round(auc_lite - auc_safe, 4)
      ))
    }
  }
  return(sub_results)
}

# Create binary subgroup columns
fh_pooled <- fh_pooled %>%
  mutate(
    sex_label  = ifelse(sex == 1, "Male", "Female"),
    age_group  = cut(age, breaks = c(0, 50, 60, 70, 100),
                     labels = c("<50", "50-60", "60-70", ">70")),
    ldl_group  = cut(ldl, breaks = c(0, 3, 4.5, 100),
                     labels = c("<3.0", "3.0-4.5", ">4.5")),
    lpa_high   = ifelse(!is.na(lpa_143) & lpa_143 == 1, "High", "Normal")
  )

# Run subgroups for each direction
all_subgroups <- data.frame()
for (dir in directions) {
  train_df <- fh_pooled %>% filter(cohort == dir[1])
  test_df  <- fh_pooled %>% filter(cohort == dir[2])
  label    <- dir[3]

  all_subgroups <- rbind(all_subgroups,
    run_subgroup(train_df, test_df, "sex_label", "Sex", label),
    run_subgroup(train_df, test_df, "age_group", "Age Group", label),
    run_subgroup(train_df, test_df, "gene", "Gene", label),
    run_subgroup(train_df, test_df, "ldl_group", "LDL Group", label)
  )
}

if (nrow(all_subgroups) > 0) {
  cat("Subgroup results:\n")
  print(all_subgroups)
  write.csv(all_subgroups, "Subgroup_analysis_results.csv", row.names = FALSE)
}

################################################################################
# ── 13. PHYSICAL SIGNS & ASCVD ──────────────────────────────────────────────
################################################################################
cat("\n=== PHYSICAL SIGNS ANALYSIS ===\n")

phys_data <- fh_pooled %>%
  filter(!is.na(has_phys_sign), !is.na(ascvd))

phys_summary <- phys_data %>%
  group_by(cohort, has_phys_sign) %>%
  summarise(
    n      = n(),
    events = sum(ascvd, na.rm = TRUE),
    rate   = round(100 * events / n, 1),
    .groups = "drop"
  )
cat("ASCVD event rate by physical signs:\n")
print(phys_summary)

# Individual signs
for (sign in c("corneal_arcus", "tendon_xanth", "xanthelasmas")) {
  sign_data <- fh_pooled %>% filter(!is.na(.data[[sign]]), !is.na(ascvd))
  if (nrow(sign_data) > 20) {
    sign_tab <- sign_data %>%
      group_by(cohort, .data[[sign]]) %>%
      summarise(
        n      = n(),
        events = sum(ascvd, na.rm = TRUE),
        rate   = round(100 * events / n, 1),
        .groups = "drop"
      )
    cat(sprintf("\n%s:\n", sign))
    print(sign_tab)
  }
}

################################################################################
# ── 14. FIGURES ─────────────────────────────────────────────────────────────
################################################################################
cat("\n=== GENERATING FIGURES ===\n")

# ── FIGURE 1: Study Design & Cohort Overview ─────────────────────────────────
# Panel A: Cohort flow / N per cohort (bar chart)
fig1a <- fh_pooled %>%
  group_by(cohort) %>%
  summarise(n = n(), .groups = "drop") %>%
  ggplot(aes(x = cohort, y = n, fill = cohort)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0("N=", comma(n))), vjust = -0.5, size = 3) +
  scale_fill_manual(values = pal3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Number of Patients", title = "a") +
  theme_nature()

# Panel B: ASCVD event rate by cohort
fig1b <- fh_pooled %>%
  filter(!is.na(ascvd)) %>%
  group_by(cohort) %>%
  summarise(
    rate = 100 * mean(ascvd, na.rm = TRUE),
    n    = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = cohort, y = rate, fill = cohort)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", rate)), vjust = -0.5, size = 3) +
  scale_fill_manual(values = pal3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "ASCVD Event Rate (%)", title = "b") +
  theme_nature()

# Panel C: Gene distribution (stacked bar)
fig1c <- fh_pooled %>%
  filter(!is.na(gene)) %>%
  count(cohort, gene) %>%
  group_by(cohort) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(x = cohort, y = pct, fill = gene)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = pal_gene, name = "Gene") +
  labs(x = NULL, y = "Proportion (%)", title = "c") +
  theme_nature()

fig1 <- plot_grid(fig1a, fig1b, fig1c, nrow = 1, rel_widths = c(1, 1, 1.2))
safe_ggsave("Figure1_cohort_overview.pdf", fig1, width = 180, height = 60, units = "mm", dpi = 300)
ggsave("Figure1_cohort_overview.png", fig1, width = 180, height = 60, units = "mm", dpi = 300)
cat("Figure 1 saved.\n")

# ── FIGURE 2: Lipid Profiles Comparison ──────────────────────────────────────
lipid_long <- fh_pooled %>%
  select(cohort, tc, ldl, hdl, tg, nhdl) %>%
  pivot_longer(cols = c(tc, ldl, hdl, tg, nhdl), names_to = "lipid", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(lipid = factor(lipid, levels = c("tc", "ldl", "hdl", "tg", "nhdl"),
                        labels = c("TC", "LDL-C", "HDL-C", "TG", "Non-HDL-C")))

fig2 <- ggplot(lipid_long, aes(x = cohort, y = value, fill = cohort)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  facet_wrap(~lipid, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal3, guide = "none") +
  labs(x = NULL, y = "Concentration (mmol/L)") +
  theme_nature(base_size = 7)

safe_ggsave("Figure2_lipid_profiles.pdf", fig2, width = 180, height = 60, units = "mm", dpi = 300)
ggsave("Figure2_lipid_profiles.png", fig2, width = 180, height = 60, units = "mm", dpi = 300)
cat("Figure 2 saved.\n")

# ── FIGURE 3: Lp(a) Distribution ─────────────────────────────────────────────
lpa_plot <- fh_pooled %>% filter(!is.na(lpa_nmol), lpa_nmol < 500)

fig3 <- ggplot(lpa_plot, aes(x = lpa_nmol, fill = cohort)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 143, linetype = "dashed", colour = "red", size = 0.8) +
  annotate("text", x = 155, y = Inf, label = "CALON threshold\n143 nmol/L",
           vjust = 1.5, hjust = 0, size = 2.5, colour = "red") +
  geom_vline(xintercept = 125, linetype = "dotted", colour = "blue", size = 0.6) +
  annotate("text", x = 85, y = Inf, label = "SAFEHEART/Paquette\n~50 mg/dL",
           vjust = 1.5, hjust = 1, size = 2.5, colour = "blue") +
  scale_fill_manual(values = pal3, name = "Cohort") +
  labs(x = "Lp(a) (nmol/L)", y = "Count",
       title = "Lp(a) Distribution Across FH Cohorts") +
  theme_nature()

safe_ggsave("Figure3_lpa_distribution.pdf", fig3, width = 120, height = 80, units = "mm", dpi = 300)
ggsave("Figure3_lpa_distribution.png", fig3, width = 120, height = 80, units = "mm", dpi = 300)
cat("Figure 3 saved.\n")

# ── FIGURE 4: ApoB vs LDL Regression (FH vs Non-FH) ────────────────────────
if (exists("lm_fh") & exists("lm_nonfh")) {
  plot_apob <- bind_rows(
    fh_apob_ldl %>% mutate(group = "FH", apob_val = apob, ldl_val = ldl),
    nonfh_apob_ldl %>% mutate(group = "Non-FH", apob_val = apob_raw, ldl_val = ldl_raw)
  ) %>%
    filter(ldl_val < 10, apob_val < 3)

  fig4 <- ggplot(plot_apob, aes(x = ldl_val, y = apob_val, colour = group)) +
    geom_point(alpha = 0.15, size = 0.5) +
    geom_smooth(method = "lm", se = TRUE, size = 0.8) +
    scale_colour_manual(values = pal2, name = NULL) +
    labs(x = "LDL-C (mmol/L)", y = "ApoB (g/L)",
         title = "ApoB vs LDL-C: FH vs Non-FH (UKB)") +
    annotate("text", x = 1, y = 2.5,
             label = sprintf("FH: slope=%.3f, R\u00B2=%.3f\nNon-FH: slope=%.3f, R\u00B2=%.3f",
                             coef(lm_fh)[2], summary(lm_fh)$r.squared,
                             coef(lm_nonfh)[2], summary(lm_nonfh)$r.squared),
             hjust = 0, size = 2.5) +
    theme_nature()

  safe_ggsave("Figure4_apob_vs_ldl.pdf", fig4, width = 120, height = 90, units = "mm", dpi = 300)
  ggsave("Figure4_apob_vs_ldl.png", fig4, width = 120, height = 90, units = "mm", dpi = 300)
  cat("Figure 4 saved.\n")
}

# ── FIGURE 5: Bi-External Validation Forest Plot ────────────────────────────
if (nrow(results) > 0) {
  fig5_data <- results %>%
    mutate(label = paste(direction, model, sep = " | "))

  fig5 <- ggplot(fig5_data, aes(x = auc, y = fct_rev(label), colour = model)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c("CALON-X2" = "#E64B35", "SAFEHEART" = "#4DBBD5"),
                        name = "Model") +
    labs(x = "AUC (95% CI)", y = NULL,
         title = "Bi-External Validation: CALON-X2 vs SAFEHEART") +
    theme_nature() +
    theme(axis.text.y = element_text(size = 7))

  safe_ggsave("Figure5_forest_plot.pdf", fig5, width = 140, height = 80, units = "mm", dpi = 300)
  ggsave("Figure5_forest_plot.png", fig5, width = 140, height = 80, units = "mm", dpi = 300)
  cat("Figure 5 saved.\n")
}

# ── FIGURE 6: Subgroup Forest Plot ──────────────────────────────────────────
if (nrow(all_subgroups) > 0) {
  fig6 <- all_subgroups %>%
    mutate(label = paste(subgroup, value, sep = ": ")) %>%
    ggplot(aes(x = delta, y = fct_rev(label), colour = direction)) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_brewer(palette = "Set1", name = "Direction") +
    labs(x = "\u0394 AUC (CALON-X2 \u2212 SAFEHEART)", y = NULL,
         title = "Subgroup Analysis") +
    theme_nature() +
    theme(axis.text.y = element_text(size = 6))

  safe_ggsave("Figure6_subgroup_forest.pdf", fig6, width = 160, height = 120, units = "mm", dpi = 300)
  ggsave("Figure6_subgroup_forest.png", fig6, width = 160, height = 120, units = "mm", dpi = 300)
  cat("Figure 6 saved.\n")
}

# ── FIGURE 7: ASCVD Event Rates by Physical Signs, Gene, Lp(a) ─────────────
# Panel A: Physical signs
fig7a_data <- phys_data %>%
  group_by(cohort, has_phys_sign) %>%
  summarise(rate = 100 * mean(ascvd, na.rm = TRUE), n = n(), .groups = "drop") %>%
  mutate(sign = ifelse(has_phys_sign == 1, "Present", "Absent"))

fig7a <- ggplot(fig7a_data, aes(x = cohort, y = rate, fill = sign)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", rate)),
            position = position_dodge(width = 0.6), vjust = -0.5, size = 2.5) +
  scale_fill_manual(values = c("Present" = "#E64B35", "Absent" = "#4DBBD5"),
                    name = "Physical Signs") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "ASCVD Event Rate (%)", title = "a  Physical Signs") +
  theme_nature()

# Panel B: Gene
fig7b_data <- fh_pooled %>%
  filter(!is.na(gene), !is.na(ascvd)) %>%
  group_by(cohort, gene) %>%
  summarise(rate = 100 * mean(ascvd, na.rm = TRUE), n = n(), .groups = "drop")

fig7b <- ggplot(fig7b_data, aes(x = cohort, y = rate, fill = gene)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = pal_gene, name = "Gene") +
  labs(x = NULL, y = "ASCVD Event Rate (%)", title = "b  Gene") +
  theme_nature()

# Panel C: Lp(a) > 143
fig7c_data <- fh_pooled %>%
  filter(!is.na(lpa_143), !is.na(ascvd)) %>%
  group_by(cohort, lpa_143) %>%
  summarise(rate = 100 * mean(ascvd, na.rm = TRUE), n = n(), .groups = "drop") %>%
  mutate(lpa_label = ifelse(lpa_143 == 1, ">143 nmol/L", "\u2264143 nmol/L"))

fig7c <- ggplot(fig7c_data, aes(x = cohort, y = rate, fill = lpa_label)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c(">143 nmol/L" = "#E64B35", "\u2264143 nmol/L" = "#4DBBD5"),
                    name = "Lp(a)") +
  labs(x = NULL, y = "ASCVD Event Rate (%)", title = "c  Lp(a)") +
  theme_nature()

fig7 <- plot_grid(fig7a, fig7b, fig7c, nrow = 1, rel_widths = c(1.1, 1.1, 1))
safe_ggsave("Figure7_event_rates.pdf", fig7, width = 180, height = 70, units = "mm", dpi = 300)
ggsave("Figure7_event_rates.png", fig7, width = 180, height = 70, units = "mm", dpi = 300)
cat("Figure 7 saved.\n")

################################################################################
# ── 14B. LANCET-CALIBER STATISTICS ──────────────────────────────────────────
#     Comprehensive model evaluation beyond discrimination (AUC)
#     Includes: Calibration, Brier Score, NRI, IDI, DCA, E-values, NNS
################################################################################
cat("\n=== LANCET-CALIBER STATISTICAL EVALUATION ===\n")

# ── Helper: Brier Score and Decomposition ────────────────────────────────────
brier_score <- function(observed, predicted) {
  bs <- mean((predicted - observed)^2)
  # Decomposition (Murphy 1973)
  n <- length(observed)
  o_bar <- mean(observed)
  # Reliability (calibration component)
  groups <- cut(predicted, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
  grp_data <- data.frame(obs = observed, pred = predicted, grp = groups) %>%
    group_by(grp) %>%
    summarise(n_g = n(), o_g = mean(obs), p_g = mean(pred), .groups = "drop")
  reliability <- sum(grp_data$n_g * (grp_data$p_g - grp_data$o_g)^2) / n
  # Resolution (discrimination component)
  resolution <- sum(grp_data$n_g * (grp_data$o_g - o_bar)^2) / n
  # Uncertainty (irreducible)
  uncertainty <- o_bar * (1 - o_bar)
  list(brier = bs, reliability = reliability, resolution = resolution,
       uncertainty = uncertainty, brier_scaled = 1 - bs / uncertainty)
}

# ── Helper: Hosmer-Lemeshow Test ─────────────────────────────────────────────
hosmer_lemeshow <- function(observed, predicted, g = 10) {
  grp <- cut(predicted, breaks = quantile(predicted, probs = seq(0, 1, 1/g)),
             include.lowest = TRUE)
  tab <- data.frame(obs = observed, pred = predicted, grp = grp) %>%
    group_by(grp) %>%
    summarise(n = n(), obs_events = sum(obs), exp_events = sum(pred),
              obs_nonevents = n - obs_events,
              exp_nonevents = n - exp_events, .groups = "drop") %>%
    filter(!is.na(grp))
  chisq <- sum((tab$obs_events - tab$exp_events)^2 / tab$exp_events +
               (tab$obs_nonevents - tab$exp_nonevents)^2 / tab$exp_nonevents,
               na.rm = TRUE)
  df <- nrow(tab) - 2
  pval <- 1 - pchisq(chisq, df)
  list(chisq = chisq, df = df, pval = pval, table = tab)
}

# ── Helper: Net Reclassification Improvement (NRI) ───────────────────────────
compute_nri <- function(obs, p_new, p_old, threshold = NULL) {
  if (is.null(threshold)) {
    # Category-free NRI (Pencina 2011)
    events <- obs == 1
    nonevents <- obs == 0
    p_up_events <- mean(p_new[events] > p_old[events])
    p_down_events <- mean(p_new[events] < p_old[events])
    p_up_nonevents <- mean(p_new[nonevents] > p_old[nonevents])
    p_down_nonevents <- mean(p_new[nonevents] < p_old[nonevents])
    nri_events <- p_up_events - p_down_events
    nri_nonevents <- p_down_nonevents - p_up_nonevents
    nri_total <- nri_events + nri_nonevents
  } else {
    # Category-based NRI
    cat_new <- ifelse(p_new >= threshold, 1, 0)
    cat_old <- ifelse(p_old >= threshold, 1, 0)
    events <- obs == 1
    nonevents <- obs == 0
    nri_events <- mean(cat_new[events] > cat_old[events]) -
                  mean(cat_new[events] < cat_old[events])
    nri_nonevents <- mean(cat_new[nonevents] < cat_old[nonevents]) -
                     mean(cat_new[nonevents] > cat_old[nonevents])
    nri_total <- nri_events + nri_nonevents
  }
  list(nri = nri_total, nri_events = nri_events, nri_nonevents = nri_nonevents)
}

# ── Helper: Integrated Discrimination Improvement (IDI) ──────────────────────
compute_idi <- function(obs, p_new, p_old) {
  events <- obs == 1
  nonevents <- obs == 0
  is_events <- mean(p_new[events]) - mean(p_old[events])
  ip_nonevents <- mean(p_new[nonevents]) - mean(p_old[nonevents])
  idi <- is_events - ip_nonevents
  # Discrimination slopes
  ds_new <- mean(p_new[events]) - mean(p_new[nonevents])
  ds_old <- mean(p_old[events]) - mean(p_old[nonevents])
  list(idi = idi, is_events = is_events, ip_nonevents = ip_nonevents,
       discrimination_slope_new = ds_new, discrimination_slope_old = ds_old)
}

# ── Helper: Decision Curve Analysis (DCA) ────────────────────────────────────
dca_analysis <- function(obs, pred, thresholds = seq(0.01, 0.99, 0.01)) {
  n <- length(obs)
  prevalence <- mean(obs)
  results <- lapply(thresholds, function(pt) {
    # Model net benefit
    tp <- sum(pred >= pt & obs == 1)
    fp <- sum(pred >= pt & obs == 0)
    nb_model <- tp / n - fp / n * (pt / (1 - pt))
    # Treat-all net benefit
    nb_all <- prevalence - (1 - prevalence) * (pt / (1 - pt))
    data.frame(threshold = pt, nb_model = nb_model,
               nb_all = nb_all, nb_none = 0)
  })
  do.call(rbind, results)
}

# ── Helper: E-value Calculation (VanderWeele & Ding 2017) ────────────────────
compute_evalue <- function(or, ci_lower = NULL) {
  # Convert OR to risk ratio approximation (conservative)
  rr <- or
  evalue_point <- rr + sqrt(rr * (rr - 1))
  if (!is.null(ci_lower)) {
    if (ci_lower > 1) {
      evalue_ci <- ci_lower + sqrt(ci_lower * (ci_lower - 1))
    } else {
      evalue_ci <- 1
    }
  } else {
    evalue_ci <- NA
  }
  list(evalue_point = evalue_point, evalue_ci = evalue_ci)
}

# ── Helper: Number Needed to Screen (NNS) ───────────────────────────────────
compute_nns <- function(obs, pred, cutoff) {
  high_risk <- pred >= cutoff
  if (sum(high_risk) == 0) return(list(nns = Inf, ppv = 0, sensitivity = 0))
  tp <- sum(high_risk & obs == 1)
  fp <- sum(high_risk & obs == 0)
  fn <- sum(!high_risk & obs == 1)
  ppv <- tp / (tp + fp)
  sensitivity <- tp / (tp + fn)
  nns <- ifelse(ppv > 0, 1 / ppv, Inf)
  list(nns = round(nns, 1), ppv = round(ppv, 4),
       sensitivity = round(sensitivity, 4),
       specificity = round(1 - fp / sum(obs == 0), 4))
}

# ────────────────────────────────────────────────────────────────────────────
# APPLY LANCET STATISTICS TO EACH VALIDATION DIRECTION
# ────────────────────────────────────────────────────────────────────────────

# Bridge: create val_pairs data frame and cohort_list from `directions`
val_pairs <- data.frame(
  train = sapply(directions, `[`, 1),
  test  = sapply(directions, `[`, 2),
  label = sapply(directions, `[`, 3),
  stringsAsFactors = FALSE
)
cohort_list <- list(
  "South Wales FH" = fh_pooled %>% filter(cohort == "South Wales FH"),
  "UKB"      = fh_pooled %>% filter(cohort == "UKB"),
  "Wales"    = fh_pooled %>% filter(cohort == "Wales")
)
cat("  val_pairs created:", nrow(val_pairs), "validation directions\n")

lancet_results <- list()

for (i in 1:nrow(val_pairs)) {
  train_name <- val_pairs$train[i]
  test_name  <- val_pairs$test[i]
  dir_label  <- paste(train_name, "->", test_name)
  cat(sprintf("\n--- Lancet Statistics: %s ---\n", dir_label))

  # Use CALON-Lite-noLpa features for Lancet statistics (works for ALL directions)
  lite_vars <- c("lipid_years", "sex", "diabetes", "ever_smoked", "inv_hdl")
  train_data <- cohort_list[[train_name]] %>%
    filter(complete.cases(lipid_years, sex, diabetes, ever_smoked,
                          inv_hdl, ascvd))
  test_data  <- cohort_list[[test_name]] %>%
    filter(complete.cases(lipid_years, sex, diabetes, ever_smoked,
                          inv_hdl, ascvd))

  if (nrow(train_data) < 50 | nrow(test_data) < 50) {
    cat("  Insufficient data, skipping.\n")
    next
  }

  # Standardise continuous variables
  cont_vars <- c("lipid_years", "inv_hdl")
  means_tr <- colMeans(train_data[, cont_vars])
  sds_tr   <- apply(train_data[, cont_vars], 2, sd)

  train_std <- train_data
  test_std  <- test_data
  for (v in cont_vars) {
    train_std[[v]] <- (train_data[[v]] - means_tr[v]) / sds_tr[v]
    test_std[[v]]  <- (test_data[[v]] - means_tr[v]) / sds_tr[v]
  }

  # Fit CALON-Lite (primary) and SAFEHEART baseline
  fit_lite <- glm(ascvd ~ lipid_years + sex + diabetes + ever_smoked + inv_hdl,
                  data = train_std, family = binomial)
  # SAFEHEART-like baseline (uses age+ldl instead of lipid_years for comparison)
  train_data_safe <- cohort_list[[train_name]] %>%
    filter(complete.cases(age, sex, ldl, inv_hdl, ever_smoked,
                          diabetes, hypertension, ascvd))
  test_data_safe  <- cohort_list[[test_name]] %>%
    filter(complete.cases(age, sex, ldl, inv_hdl, ever_smoked,
                          diabetes, hypertension, ascvd))
  fit_sh <- if (nrow(train_data_safe) > 50 & nrow(test_data_safe) > 50) {
    glm(ascvd ~ age + sex + ldl + inv_hdl + ever_smoked +
          diabetes + hypertension,
        data = train_data_safe, family = binomial)
  } else NULL

  pred_x2 <- predict(fit_lite, newdata = test_std, type = "response")
  pred_sh <- if (!is.null(fit_sh) && nrow(test_data_safe) > 50) {
    predict(fit_sh, newdata = test_data_safe, type = "response")
  } else rep(NA, nrow(test_std))
  obs     <- test_std$ascvd

  # 1. Brier Score
  bs_x2 <- brier_score(obs, pred_x2)
  bs_sh <- brier_score(obs, pred_sh)
  cat(sprintf("  Brier Score CALON-Lite: %.4f  (scaled: %.4f)\n",
              bs_x2$brier, bs_x2$brier_scaled))
  cat(sprintf("  Brier Score SAFEHEART:  %.4f (scaled: %.4f)\n",
              bs_sh$brier, bs_sh$brier_scaled))
  cat(sprintf("  Delta Brier: %.4f (negative = CALON-Lite better)\n",
              bs_x2$brier - bs_sh$brier))

  # 2. Hosmer-Lemeshow Calibration Test
  hl_x2 <- tryCatch(hosmer_lemeshow(obs, pred_x2), error = function(e) NULL)
  if (!is.null(hl_x2)) {
    cat(sprintf("  Hosmer-Lemeshow: chi2=%.2f, df=%d, p=%.4f %s\n",
                hl_x2$chisq, hl_x2$df, hl_x2$pval,
                ifelse(hl_x2$pval > 0.05, "(good calibration)", "(poor calibration)")))
  }

  # 3. Calibration slope and intercept
  cal_model <- glm(obs ~ pred_x2, family = binomial)
  cal_intercept <- coef(cal_model)[1]
  cal_slope     <- coef(cal_model)[2]
  cat(sprintf("  Calibration-in-the-large (intercept): %.4f\n", cal_intercept))
  cat(sprintf("  Calibration slope: %.4f (ideal=1)\n", cal_slope))

  # 4. NRI (Category-free)
  nri <- compute_nri(obs, pred_x2, pred_sh)
  cat(sprintf("  Category-free NRI: %.4f (events=%.4f, non-events=%.4f)\n",
              nri$nri, nri$nri_events, nri$nri_nonevents))

  # Bootstrap NRI CI
  nri_boot <- replicate(500, {
    idx <- sample(length(obs), replace = TRUE)
    tryCatch({
      n <- compute_nri(obs[idx], pred_x2[idx], pred_sh[idx])
      n$nri
    }, error = function(e) NA)
  })
  nri_ci <- quantile(nri_boot, c(0.025, 0.975), na.rm = TRUE)
  cat(sprintf("  NRI 95%% CI: [%.4f, %.4f]\n", nri_ci[1], nri_ci[2]))

  # 5. IDI
  idi <- compute_idi(obs, pred_x2, pred_sh)
  cat(sprintf("  IDI: %.4f (discrimination slope CALON=%.4f, SAFE=%.4f)\n",
              idi$idi, idi$discrimination_slope_new, idi$discrimination_slope_old))

  # Bootstrap IDI CI
  idi_boot <- replicate(500, {
    idx <- sample(length(obs), replace = TRUE)
    tryCatch({
      id <- compute_idi(obs[idx], pred_x2[idx], pred_sh[idx])
      id$idi
    }, error = function(e) NA)
  })
  idi_ci <- quantile(idi_boot, c(0.025, 0.975), na.rm = TRUE)
  cat(sprintf("  IDI 95%% CI: [%.4f, %.4f]\n", idi_ci[1], idi_ci[2]))

  # 6. NNS at various thresholds
  cat("  Number Needed to Screen (NNS):\n")
  for (thresh in c(0.10, 0.15, 0.20, 0.30)) {
    nns <- compute_nns(obs, pred_x2, thresh)
    cat(sprintf("    Threshold %.0f%%: NNS=%.1f, PPV=%.1f%%, Sensitivity=%.1f%%, Specificity=%.1f%%\n",
                thresh * 100, nns$nns, nns$ppv * 100, nns$sensitivity * 100,
                nns$specificity * 100))
  }

  # 7. E-values for pooled OR
  cat("  E-values (pooled model predictors):\n")
  if (exists("fit_pooled")) {
    or_tidy <- tidy(fit_pooled, exponentiate = TRUE, conf.int = TRUE)
    for (j in 2:nrow(or_tidy)) {
      ev <- compute_evalue(or_tidy$estimate[j], or_tidy$conf.low[j])
      cat(sprintf("    %s: OR=%.2f, E-value=%.2f (CI bound=%.2f)\n",
                  or_tidy$term[j], or_tidy$estimate[j],
                  ev$evalue_point, ev$evalue_ci))
    }
  }

  # 8. Absolute Risk Difference (high vs low predicted risk)
  risk_q <- quantile(pred_x2, c(0.25, 0.75))
  low_risk  <- obs[pred_x2 <= risk_q[1]]
  high_risk <- obs[pred_x2 >= risk_q[2]]
  ard <- mean(high_risk) - mean(low_risk)
  cat(sprintf("  Absolute Risk Difference (Q4 vs Q1): %.1f%% (%.1f%% vs %.1f%%)\n",
              ard * 100, mean(high_risk) * 100, mean(low_risk) * 100))

  # Store results
  lancet_results[[dir_label]] <- list(
    brier_x2 = bs_x2, brier_sh = bs_sh,
    hl = hl_x2, cal_intercept = cal_intercept, cal_slope = cal_slope,
    nri = nri, nri_ci = nri_ci, idi = idi, idi_ci = idi_ci,
    pred_x2 = pred_x2, pred_sh = pred_sh, obs = obs,
    dca = dca_analysis(obs, pred_x2)
  )
}

# Write Lancet statistics summary
lancet_summary <- data.frame(
  Direction = names(lancet_results),
  Brier_X2 = sapply(lancet_results, function(x) round(x$brier_x2$brier, 4)),
  Brier_SH = sapply(lancet_results, function(x) round(x$brier_sh$brier, 4)),
  Brier_Scaled = sapply(lancet_results, function(x) round(x$brier_x2$brier_scaled, 4)),
  Cal_Intercept = sapply(lancet_results, function(x) round(x$cal_intercept, 4)),
  Cal_Slope = sapply(lancet_results, function(x) round(x$cal_slope, 4)),
  NRI = sapply(lancet_results, function(x) round(x$nri$nri, 4)),
  NRI_CI_lo = sapply(lancet_results, function(x) round(x$nri_ci[1], 4)),
  NRI_CI_hi = sapply(lancet_results, function(x) round(x$nri_ci[2], 4)),
  IDI = sapply(lancet_results, function(x) round(x$idi$idi, 4)),
  IDI_CI_lo = sapply(lancet_results, function(x) round(x$idi_ci[1], 4)),
  IDI_CI_hi = sapply(lancet_results, function(x) round(x$idi_ci[2], 4)),
  stringsAsFactors = FALSE
)
write.csv(lancet_summary, "Lancet_statistics_summary.csv", row.names = FALSE)
cat("\nLancet statistics saved to Lancet_statistics_summary.csv\n")


################################################################################
# ── FIGURE 8: CALON-X2 CLINICAL PERFORMANCE DASHBOARD ───────────────────────
#     Nature-caliber six-panel figure with Lancet statistics
#     (a) Calibration Plot  (b) Decision Curve Analysis
#     (c) Risk Distribution (d) NRI & IDI Forest
#     (e) Brier Decomposition (f) E-value Sensitivity
################################################################################
cat("\n=== FIGURE 8: CLINICAL PERFORMANCE DASHBOARD ===\n")

# Use the first available validation result for the main panels
lr_key <- names(lancet_results)[1]
lr <- lancet_results[[lr_key]]

# ── Panel A: CALIBRATION PLOT ────────────────────────────────────────────────
cal_data <- data.frame(obs = lr$obs, pred = lr$pred_x2) %>%
  mutate(decile = ntile(pred, 10)) %>%
  group_by(decile) %>%
  summarise(
    observed   = mean(obs),
    predicted  = mean(pred),
    n          = n(),
    se         = sqrt(observed * (1 - observed) / n),
    lower      = pmax(0, observed - 1.96 * se),
    upper      = pmin(1, observed + 1.96 * se),
    .groups = "drop"
  )

fig8a <- ggplot(cal_data, aes(x = predicted, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey50", linewidth = 0.5) +
  geom_point(aes(size = n), colour = "#E64B35", alpha = 0.9) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.008, colour = "#E64B35", linewidth = 0.4) +
  geom_smooth(method = "loess", se = FALSE, colour = "#3C5488",
              linewidth = 0.7, span = 1) +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  annotate("text", x = 0.05, y = 0.92,
           label = sprintf("Cal. slope = %.2f\nCal. intercept = %.3f",
                           lr$cal_slope, lr$cal_intercept),
           hjust = 0, size = 2.2, colour = "grey30") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted Probability", y = "Observed Proportion",
       title = "a  Calibration") +
  theme_nature(base_size = 7)

# ── Panel B: DECISION CURVE ANALYSIS ────────────────────────────────────────
dca_data <- lr$dca %>% filter(threshold <= 0.6)

fig8b <- ggplot(dca_data) +
  geom_line(aes(x = threshold, y = nb_model, colour = "CALON-X2"),
            linewidth = 0.8) +
  geom_line(aes(x = threshold, y = nb_all, colour = "Treat All"),
            linewidth = 0.6, linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4,
             linetype = "dotted") +
  scale_colour_manual(values = c("CALON-X2" = "#E64B35",
                                  "Treat All" = "#4DBBD5"),
                      name = NULL) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(x = "Threshold Probability", y = "Net Benefit",
       title = "b  Decision Curve Analysis") +
  theme_nature(base_size = 7) +
  theme(legend.position = c(0.7, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5))

# ── Panel C: RISK DISTRIBUTION (Density) ────────────────────────────────────
risk_df <- data.frame(
  pred = lr$pred_x2,
  outcome = factor(lr$obs, levels = c(0, 1),
                   labels = c("No ASCVD", "ASCVD"))
)

# Discrimination slope
ds_val <- mean(lr$pred_x2[lr$obs == 1]) - mean(lr$pred_x2[lr$obs == 0])

fig8c <- ggplot(risk_df, aes(x = pred, fill = outcome)) +
  geom_density(alpha = 0.55, colour = NA) +
  geom_vline(xintercept = mean(lr$pred_x2[lr$obs == 1]),
             linetype = "dashed", colour = "#E64B35", linewidth = 0.4) +
  geom_vline(xintercept = mean(lr$pred_x2[lr$obs == 0]),
             linetype = "dashed", colour = "#4DBBD5", linewidth = 0.4) +
  annotate("segment",
           x = mean(lr$pred_x2[lr$obs == 0]),
           xend = mean(lr$pred_x2[lr$obs == 1]),
           y = max(density(lr$pred_x2)$y) * 0.9,
           yend = max(density(lr$pred_x2)$y) * 0.9,
           arrow = arrow(length = unit(0.15, "cm"), ends = "both"),
           colour = "grey30", linewidth = 0.4) +
  annotate("text",
           x = mean(c(mean(lr$pred_x2[lr$obs == 0]),
                       mean(lr$pred_x2[lr$obs == 1]))),
           y = max(density(lr$pred_x2)$y) * 0.95,
           label = sprintf("DS = %.3f", ds_val),
           size = 2.2, colour = "grey20") +
  scale_fill_manual(values = c("No ASCVD" = "#4DBBD5", "ASCVD" = "#E64B35"),
                    name = NULL) +
  labs(x = "Predicted Probability", y = "Density",
       title = "c  Risk Distribution") +
  theme_nature(base_size = 7) +
  theme(legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5))

# ── Panel D: NRI & IDI FOREST PLOT ──────────────────────────────────────────
nri_idi_data <- data.frame(
  metric = c(rep("NRI", length(lancet_results)),
             rep("IDI", length(lancet_results))),
  direction = rep(names(lancet_results), 2),
  estimate = c(
    sapply(lancet_results, function(x) x$nri$nri),
    sapply(lancet_results, function(x) x$idi$idi)
  ),
  ci_lo = c(
    sapply(lancet_results, function(x) x$nri_ci[1]),
    sapply(lancet_results, function(x) x$idi_ci[1])
  ),
  ci_hi = c(
    sapply(lancet_results, function(x) x$nri_ci[2]),
    sapply(lancet_results, function(x) x$idi_ci[2])
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(label = paste(metric, direction, sep = "\n"))

fig8d <- ggplot(nri_idi_data, aes(x = estimate, y = fct_rev(label),
                                   colour = metric)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50",
             linewidth = 0.4) +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25,
                 linewidth = 0.5) +
  scale_colour_manual(values = c("NRI" = "#E64B35", "IDI" = "#3C5488"),
                      name = NULL) +
  labs(x = "Estimate (95% Bootstrap CI)", y = NULL,
       title = "d  NRI & IDI") +
  theme_nature(base_size = 7) +
  theme(axis.text.y = element_text(size = 5),
        legend.position = c(0.85, 0.15),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5))

# ── Panel E: BRIER SCORE DECOMPOSITION ──────────────────────────────────────
brier_data <- data.frame(
  direction = rep(names(lancet_results), each = 3),
  component = rep(c("Reliability\n(Calibration)", "Resolution\n(Discrimination)",
                     "Uncertainty\n(Irreducible)"), length(lancet_results)),
  value = unlist(lapply(lancet_results, function(x) {
    c(x$brier_x2$reliability, x$brier_x2$resolution, x$brier_x2$uncertainty)
  })),
  stringsAsFactors = FALSE
)

fig8e <- ggplot(brier_data, aes(x = direction, y = value, fill = component)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = c("Reliability\n(Calibration)" = "#E64B35",
                                "Resolution\n(Discrimination)" = "#4DBBD5",
                                "Uncertainty\n(Irreducible)" = "#F39B7F"),
                    name = NULL) +
  labs(x = NULL, y = "Brier Score Component",
       title = "e  Brier Decomposition") +
  theme_nature(base_size = 7) +
  theme(axis.text.x = element_text(size = 5, angle = 30, hjust = 1),
        legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 4.5))

# Add Brier annotation
brier_labels <- data.frame(
  direction = names(lancet_results),
  brier = sapply(lancet_results, function(x) x$brier_x2$brier),
  total = sapply(lancet_results, function(x) {
    x$brier_x2$reliability + x$brier_x2$resolution + x$brier_x2$uncertainty
  })
)
fig8e <- fig8e +
  geom_text(data = brier_labels,
            aes(x = direction, y = total + 0.01,
                label = sprintf("BS=%.3f", brier)),
            inherit.aes = FALSE, size = 2, colour = "grey30")

# ── Panel F: E-VALUE SENSITIVITY ANALYSIS ───────────────────────────────────
if (exists("fit_pooled")) {
  or_tidy <- tidy(fit_pooled, exponentiate = TRUE, conf.int = TRUE)
  evalue_data <- or_tidy %>%
    filter(term != "(Intercept)") %>%
    rowwise() %>%
    mutate(
      evalue_pt = ifelse(estimate > 1,
                         estimate + sqrt(estimate * (estimate - 1)),
                         1 / estimate + sqrt(1 / estimate * (1 / estimate - 1))),
      evalue_ci = ifelse(is.na(conf.low) | conf.low <= 1, 1,
                         conf.low + sqrt(conf.low * (conf.low - 1)))
    ) %>%
    ungroup() %>%
    mutate(
      term_clean = case_when(
        term == "age" ~ "Age",
        term == "sex" ~ "Sex (Male)",
        term == "ldl" ~ "LDL-C",
        term == "inv_hdl" ~ "1/HDL-C",
        term == "ever_smoked" ~ "Ever Smoked",
        term == "diabetes" ~ "Diabetes",
        term == "hypertension" ~ "Hypertension",
        term == "lpa_143" ~ "Lp(a)>143",
        TRUE ~ term
      )
    )

  fig8f <- ggplot(evalue_data, aes(x = evalue_pt, y = fct_rev(term_clean))) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50",
               linewidth = 0.4) +
    geom_segment(aes(x = evalue_ci, xend = evalue_pt,
                     y = fct_rev(term_clean), yend = fct_rev(term_clean)),
                 colour = "#3C5488", linewidth = 0.8) +
    geom_point(aes(size = estimate), colour = "#E64B35", alpha = 0.9) +
    geom_point(aes(x = evalue_ci), colour = "#3C5488", size = 1.8, shape = 18) +
    scale_size_continuous(range = c(2, 5), name = "OR") +
    labs(x = "E-value", y = NULL,
         title = "f  E-value Sensitivity") +
    theme_nature(base_size = 7) +
    theme(legend.position = c(0.85, 0.25),
          legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 5))
} else {
  fig8f <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "E-values\n(Pooled model\nrequired)",
             size = 3) +
    theme_void() +
    labs(title = "f  E-value Sensitivity")
}

# ── ASSEMBLE FIGURE 8 (6-panel layout) ──────────────────────────────────────
# Top row: Calibration | DCA | Risk Distribution
# Bottom row: NRI/IDI | Brier Decomposition | E-value
top_row    <- plot_grid(fig8a, fig8b, fig8c, nrow = 1,
                        rel_widths = c(1, 1, 1))
bottom_row <- plot_grid(fig8d, fig8e, fig8f, nrow = 1,
                        rel_widths = c(1.1, 1, 1))
fig8 <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 1))

# Add overall title
fig8_titled <- ggdraw() +
  draw_label("CALON-X2 Clinical Performance Dashboard",
             fontface = "bold", size = 9,
             x = 0.5, y = 0.98, hjust = 0.5, vjust = 1) +
  draw_plot(fig8, x = 0, y = 0, width = 1, height = 0.96)

safe_ggsave("Figure8_clinical_dashboard.pdf", fig8_titled,
       width = 190, height = 150, units = "mm", dpi = 300)
ggsave("Figure8_clinical_dashboard.png", fig8_titled,
       width = 190, height = 150, units = "mm", dpi = 300)
cat("Figure 8 saved.\n")


################################################################################
# ── 15. POOLED MODEL (ALL COHORTS) — CALON-Lite + CALON-X ─────────────────
################################################################################
cat("\n=== POOLED MODEL ANALYSIS ===\n")

# --- 15a. CALON-Lite Pooled (6 predictors) ---
pooled_lite <- fh_pooled %>%
  filter(complete.cases(lipid_years, sex, diabetes, ever_smoked,
                        lpa_143, inv_hdl, ascvd))

cat(sprintf("CALON-Lite pooled complete cases: %d\n", nrow(pooled_lite)))

if (nrow(pooled_lite) > 100) {
  fit_pooled_lite <- glm(ascvd ~ lipid_years + sex + diabetes +
                           ever_smoked + lpa_143 + inv_hdl,
                         data = pooled_lite, family = binomial)
  cat("\nPooled CALON-Lite logistic regression:\n")
  print(summary(fit_pooled_lite))

  or_lite <- tidy(fit_pooled_lite, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(across(where(is.numeric), ~round(.x, 3)))
  cat("\nCALON-Lite Odds Ratios:\n")
  print(or_lite)
  write.csv(or_lite, "CALON_Lite_odds_ratios.csv", row.names = FALSE)

  pred_pooled_lite <- predict(fit_pooled_lite, type = "response")
  roc_pooled_lite  <- roc(pooled_lite$ascvd, pred_pooled_lite, quiet = TRUE)
  cat(sprintf("\nPooled CALON-Lite AUC: %.4f (95%% CI: %.4f-%.4f)\n",
              auc(roc_pooled_lite), ci.auc(roc_pooled_lite)[1], ci.auc(roc_pooled_lite)[3]))
}

# --- 15a2. CALON-Lite-noLpa Pooled (5 predictors — for clinics without Lp(a)) ---
pooled_lite5 <- fh_pooled %>%
  filter(complete.cases(lipid_years, sex, diabetes, ever_smoked, inv_hdl, ascvd))

cat(sprintf("\nCALON-Lite-noLpa pooled complete cases: %d\n", nrow(pooled_lite5)))

if (nrow(pooled_lite5) > 100) {
  fit_pooled_lite5 <- glm(ascvd ~ lipid_years + sex + diabetes +
                            ever_smoked + inv_hdl,
                          data = pooled_lite5, family = binomial)
  cat("\nPooled CALON-Lite-noLpa logistic regression:\n")
  print(summary(fit_pooled_lite5))

  or_lite5 <- tidy(fit_pooled_lite5, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(across(where(is.numeric), ~round(.x, 3)))
  cat("\nCALON-Lite-noLpa Odds Ratios:\n")
  print(or_lite5)
  write.csv(or_lite5, "CALON_Lite_noLpa_odds_ratios.csv", row.names = FALSE)

  pred_pooled_lite5 <- predict(fit_pooled_lite5, type = "response")
  roc_pooled_lite5  <- roc(pooled_lite5$ascvd, pred_pooled_lite5, quiet = TRUE)
  cat(sprintf("\nPooled CALON-Lite-noLpa AUC: %.4f (95%% CI: %.4f-%.4f)\n",
              auc(roc_pooled_lite5), ci.auc(roc_pooled_lite5)[1], ci.auc(roc_pooled_lite5)[3]))
}

# --- 15b. CALON-X Pooled (enhanced: Lite + ApoB/LDL + stigmata) ---
pooled_x <- fh_pooled %>%
  filter(complete.cases(lipid_years, sex, diabetes, ever_smoked,
                        lpa_143, inv_hdl, log_apob_ldl, has_phys_sign, ascvd))

cat(sprintf("\nCALON-X pooled complete cases: %d\n", nrow(pooled_x)))

if (nrow(pooled_x) > 50) {
  fit_pooled_x <- glm(ascvd ~ lipid_years + sex + diabetes +
                        ever_smoked + lpa_143 + inv_hdl +
                        log_apob_ldl + has_phys_sign,
                      data = pooled_x, family = binomial)
  cat("\nPooled CALON-X logistic regression:\n")
  print(summary(fit_pooled_x))

  or_x <- tidy(fit_pooled_x, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(across(where(is.numeric), ~round(.x, 3)))
  cat("\nCALON-X Odds Ratios:\n")
  print(or_x)
  write.csv(or_x, "CALON_X_odds_ratios.csv", row.names = FALSE)

  pred_pooled_x <- predict(fit_pooled_x, type = "response")
  roc_pooled_x  <- roc(pooled_x$ascvd, pred_pooled_x, quiet = TRUE)
  cat(sprintf("\nPooled CALON-X AUC: %.4f (95%% CI: %.4f-%.4f)\n",
              auc(roc_pooled_x), ci.auc(roc_pooled_x)[1], ci.auc(roc_pooled_x)[3]))
} else {
  cat("  Insufficient CALON-X complete cases for pooled model.\n")
  fit_pooled_x <- NULL
}

# Legacy CALON-X2 pooled (for comparison)
pooled_complete <- fh_pooled %>%
  filter(complete.cases(age, sex, ldl, inv_hdl, ever_smoked,
                        diabetes, hypertension, lpa_143, ascvd))

if (nrow(pooled_complete) > 100) {
  fit_pooled <- glm(ascvd ~ age + sex + ldl + inv_hdl + ever_smoked +
                      diabetes + hypertension + lpa_143,
                    data = pooled_complete, family = binomial)
  or_table <- tidy(fit_pooled, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(across(where(is.numeric), ~round(.x, 3)))
  write.csv(or_table, "CALON_X2_odds_ratios.csv", row.names = FALSE)
  pred_pooled <- predict(fit_pooled, type = "response")
  roc_pooled  <- roc(pooled_complete$ascvd, pred_pooled, quiet = TRUE)
  cat(sprintf("\nLegacy CALON-X2 Pooled AUC: %.4f (95%% CI: %.4f-%.4f)\n",
              auc(roc_pooled), ci.auc(roc_pooled)[1], ci.auc(roc_pooled)[3]))
}

# --- 15c. LOCO-CV (Leave-One-Cohort-Out Cross-Validation) ---
cat("\n=== LOCO-CV (Outcome: Total ASCVD) ===\n")
loco_results <- data.frame()

for (holdout in c("South Wales FH", "UKB", "Wales")) {
  train_loco <- fh_pooled %>% filter(cohort != holdout)
  test_loco  <- fh_pooled %>% filter(cohort == holdout)

  # CALON-Lite-noLpa LOCO (works for all holdouts including Wales)
  res_loco_lite <- get_auc_ci(train_loco, test_loco, LITE_NOLPA_features,
                               outcome = "ascvd", n_boot = 2000)
  # CALON-Lite LOCO (may fail for Wales holdout due to Lp(a))
  res_loco_lite6 <- get_auc_ci(train_loco, test_loco, LITE_features,
                                outcome = "ascvd", n_boot = 2000)
  # CALON-X LOCO (may have limited data for some holdouts)
  res_loco_x <- get_auc_ci(train_loco, test_loco, X_features,
                             outcome = "ascvd", n_boot = 2000)

  cat(sprintf("  LOCO holdout=%s:\n    Lite-noLpa AUC=%.4f (N=%d)\n    Lite AUC=%.4f (N=%d)\n    X AUC=%.4f (N=%d)\n",
              holdout,
              ifelse(is.na(res_loco_lite$auc), -1, res_loco_lite$auc), res_loco_lite$n_test,
              ifelse(is.na(res_loco_lite6$auc), -1, res_loco_lite6$auc), res_loco_lite6$n_test,
              ifelse(is.na(res_loco_x$auc), -1, res_loco_x$auc), res_loco_x$n_test))

  loco_results <- rbind(loco_results, data.frame(
    holdout = holdout,
    model   = c("CALON-Lite-noLpa", "CALON-Lite", "CALON-X"),
    auc     = c(res_loco_lite$auc, res_loco_lite6$auc, res_loco_x$auc),
    ci_low  = c(res_loco_lite$ci_low, res_loco_lite6$ci_low, res_loco_x$ci_low),
    ci_high = c(res_loco_lite$ci_high, res_loco_lite6$ci_high, res_loco_x$ci_high),
    n_test  = c(res_loco_lite$n_test, res_loco_lite6$n_test, res_loco_x$n_test)
  ))
}
cat("\nLOCO-CV Results:\n")
print(loco_results)
write.csv(loco_results, "LOCO_CV_results.csv", row.names = FALSE)

# --- 15d. Bootstrap CI for Calibration Slope (Defense) ---
cat("\n=== CALIBRATION SLOPE DEFENSE ===\n")
# For each validation direction, compute bootstrap CI for calibration slope
for (dir_name in names(lancet_results)) {
  lr <- lancet_results[[dir_name]]
  if (is.null(lr)) next

  cal_boot <- replicate(1000, {
    idx <- sample(length(lr$obs), replace = TRUE)
    tryCatch({
      cm <- glm(lr$obs[idx] ~ lr$pred_x2[idx], family = binomial)
      coef(cm)[2]
    }, error = function(e) NA)
  })
  cal_ci <- quantile(cal_boot, c(0.025, 0.975), na.rm = TRUE)
  cat(sprintf("  %s: Cal slope = %.2f (95%% CI: %.2f - %.2f)\n",
              dir_name, lr$cal_slope, cal_ci[1], cal_ci[2]))

  # DCA net benefit at key thresholds
  if (!is.null(lr$dca)) {
    dca_nb <- lr$dca %>% filter(threshold %in% c(0.05, 0.10, 0.15, 0.20, 0.30))
    if (nrow(dca_nb) > 0) {
      cat("  DCA Net Benefit at clinical thresholds:\n")
      for (j in 1:nrow(dca_nb)) {
        cat(sprintf("    Threshold %.0f%%: Net Benefit = %.4f (vs treat-all=%.4f)\n",
                    dca_nb$threshold[j] * 100, dca_nb$net_benefit[j],
                    dca_nb$net_benefit_all[j]))
      }
    }
  }

  # Recalibration: apply slope to produce recalibrated predictions
  if (!is.na(lr$cal_slope) && lr$cal_slope != 0) {
    logit_pred <- qlogis(pmax(0.001, pmin(0.999, lr$pred_x2)))
    recal_pred <- plogis(lr$cal_intercept + lr$cal_slope * logit_pred)
    roc_recal <- roc(lr$obs, recal_pred, quiet = TRUE)
    cat(sprintf("  AUC after recalibration: %.4f (unchanged, as expected)\n",
                auc(roc_recal)))
    brier_recal <- mean((lr$obs - recal_pred)^2)
    cat(sprintf("  Brier after recalibration: %.4f (vs original: %.4f)\n",
                brier_recal, lr$brier_x2$brier))
  }
}

################################################################################
# ── 16. PRS INTEGRATION: CALON-X3 ────────────────────────────────────────────
################################################################################
cat("\n=== PRS INTEGRATION ===\n")

# Load UKB PRS data (field p26227 = CAD PRS)
prs_file <- list.files(data_dir, pattern = "calon_extra_prs\\.csv",
                       full.names = TRUE)

if (length(prs_file) > 0) {
  prs_raw <- read.csv(prs_file[1], stringsAsFactors = FALSE)
  prs_data <- prs_raw %>%
    mutate(prs_cad_raw = as.numeric(participant.p26227)) %>%
    filter(!is.na(prs_cad_raw)) %>%
    mutate(prs_cad = (prs_cad_raw - mean(prs_cad_raw)) / sd(prs_cad_raw)) %>%
    select(participant.eid, prs_cad_raw, prs_cad)

  cat(sprintf("CAD PRS loaded: %d participants\n", nrow(prs_data)))
  cat(sprintf("  Raw: Mean=%.4f, SD=%.4f\n",
              mean(prs_data$prs_cad_raw), sd(prs_data$prs_cad_raw)))

  # Merge PRS into UKB FH
  ukb_fh <- ukb_fh %>%
    left_join(prs_data, by = c("eid" = "participant.eid"))

  prs_n <- sum(!is.na(ukb_fh$prs_cad))
  cat(sprintf("PRS available for %d / %d UKB FH patients (%.1f%%)\n",
              prs_n, nrow(ukb_fh), 100 * prs_n / nrow(ukb_fh)))

  # PRS by ASCVD status
  prs_ascvd <- ukb_fh %>%
    filter(!is.na(prs_cad), !is.na(ascvd))

  if (nrow(prs_ascvd) > 50) {
    pos <- prs_ascvd %>% filter(ascvd == 1)
    neg <- prs_ascvd %>% filter(ascvd == 0)
    cat(sprintf("\n  ASCVD+ (N=%d): PRS mean=%.4f, SD=%.4f\n",
                nrow(pos), mean(pos$prs_cad), sd(pos$prs_cad)))
    cat(sprintf("  ASCVD- (N=%d): PRS mean=%.4f, SD=%.4f\n",
                nrow(neg), mean(neg$prs_cad), sd(neg$prs_cad)))

    # Welch t-test
    tt <- t.test(pos$prs_cad, neg$prs_cad)
    cat(sprintf("  t-test: t=%.3f, p=%.4e\n", tt$statistic, tt$p.value))

    # Quartile analysis
    prs_ascvd$prs_q <- ntile(prs_ascvd$prs_cad, 4)
    cat("\n  PRS Quartile ASCVD rates:\n")
    for (q in 1:4) {
      sub <- prs_ascvd %>% filter(prs_q == q)
      rate <- 100 * mean(sub$ascvd)
      cat(sprintf("    Q%d: %.1f%% (%d/%d)\n", q, rate, sum(sub$ascvd), nrow(sub)))
    }

    # Q4 vs Q1 OR
    q4 <- prs_ascvd %>% filter(prs_q == 4)
    q1 <- prs_ascvd %>% filter(prs_q == 1)
    a <- sum(q4$ascvd); b <- nrow(q4) - a
    cc <- sum(q1$ascvd); dd <- nrow(q1) - cc
    or_q4q1 <- (a * dd) / (b * cc)
    se_log <- sqrt(1/a + 1/b + 1/cc + 1/dd)
    cat(sprintf("\n  Q4 vs Q1 OR: %.2f (95%% CI: %.2f-%.2f)\n",
                or_q4q1, exp(log(or_q4q1) - 1.96*se_log),
                exp(log(or_q4q1) + 1.96*se_log)))

    # Logistic regression: ASCVD ~ PRS (unadjusted)
    fit_prs <- glm(ascvd ~ prs_cad, data = prs_ascvd, family = binomial)
    or_prs <- exp(coef(fit_prs)["prs_cad"])
    ci_prs <- exp(confint.default(fit_prs)["prs_cad", ])
    cat(sprintf("  PRS OR per SD: %.3f (95%% CI: %.3f-%.3f)\n",
                or_prs, ci_prs[1], ci_prs[2]))

    # CALON-X3: X2 + PRS
    cat("\n  --- CALON-X3 (CALON-X2 + PRS) ---\n")
    x3_data <- prs_ascvd %>%
      filter(!is.na(age), !is.na(sex), !is.na(ldl), !is.na(inv_hdl),
             !is.na(ever_smoked), !is.na(diabetes), !is.na(hypertension),
             !is.na(lpa_143), !is.na(prs_cad))

    if (nrow(x3_data) > 50 & length(unique(x3_data$ascvd)) == 2) {
      fit_x3 <- glm(ascvd ~ age + sex + ldl + inv_hdl + ever_smoked +
                       diabetes + hypertension + lpa_143 + prs_cad,
                     data = x3_data, family = binomial)
      cat("  CALON-X3 model coefficients:\n")
      x3_or <- data.frame(
        Variable = names(coef(fit_x3)),
        Coefficient = coef(fit_x3),
        OR = exp(coef(fit_x3)),
        CI_low = exp(confint.default(fit_x3)[, 1]),
        CI_high = exp(confint.default(fit_x3)[, 2]),
        p_value = summary(fit_x3)$coefficients[, 4]
      )
      print(x3_or)
      write.csv(x3_or, file.path(data_dir, "CALON_X3_PRS_coefficients.csv"),
                row.names = FALSE)
      cat("  Saved: CALON_X3_PRS_coefficients.csv\n")
    } else {
      cat("  Insufficient complete cases for CALON-X3\n")
    }
  }
} else {
  cat("PRS file not found — skipping PRS integration\n")
}

################################################################################
# ── 17. METABOLIC BURDEN IN FH ──────────────────────────────────────────────
################################################################################
cat("\n=== METABOLIC BURDEN ANALYSIS ===\n")

# Create metabolic burden score: obesity + diabetes + high TG
mb_data <- fh_pooled %>%
  filter(!is.na(ascvd)) %>%
  mutate(
    obese   = ifelse(!is.na(bmi), as.integer(bmi >= 30), NA_integer_),
    high_tg = ifelse(!is.na(tg),  as.integer(tg >= 2.3), NA_integer_)
  )

# Metabolic burden score (sum of available factors)
mb_data$metabolic_burden <- rowSums(
  mb_data[, c("obese", "diabetes", "high_tg")], na.rm = FALSE)

cat("--- DATA AVAILABILITY ---\n")
for (coh in c("South Wales FH", "UKB", "Wales")) {
  sub <- mb_data %>% filter(cohort == coh)
  cat(sprintf("  %s (N=%d): BMI=%d, TG=%d, Diabetes=%d\n",
              coh, nrow(sub),
              sum(!is.na(sub$bmi)), sum(!is.na(sub$tg)),
              sum(!is.na(sub$diabetes))))
}

# Metabolic burden score distribution
mb_valid <- mb_data %>% filter(!is.na(metabolic_burden))
cat(sprintf("\n--- METABOLIC BURDEN SCORE (N=%d) ---\n", nrow(mb_valid)))
for (s in 0:3) {
  sub <- mb_valid %>% filter(metabolic_burden == s)
  if (nrow(sub) > 0) {
    rate <- 100 * mean(sub$ascvd)
    cat(sprintf("  Score %d: N=%d, ASCVD events=%d, Rate=%.1f%%\n",
                s, nrow(sub), sum(sub$ascvd), rate))
  }
}

# Trend test (Cochran-Armitage)
cat("\n  Trend p-value: ")
trend_tab <- mb_valid %>%
  group_by(metabolic_burden) %>%
  summarise(events = sum(ascvd), n = n(), .groups = "drop")
trend_fit <- glm(ascvd ~ metabolic_burden, data = mb_valid, family = binomial)
trend_p <- summary(trend_fit)$coefficients["metabolic_burden", 4]
cat(sprintf("%.4e (OR per unit increase: %.2f)\n", trend_p,
            exp(coef(trend_fit)["metabolic_burden"])))

# Individual metabolic factor ORs
cat("\n--- INDIVIDUAL METABOLIC FACTOR ORs ---\n")
for (var in c("obese", "diabetes", "high_tg")) {
  sub <- mb_data %>% filter(!is.na(.data[[var]]), !is.na(ascvd))
  if (nrow(sub) > 30) {
    tab <- table(sub[[var]], sub$ascvd)
    if (nrow(tab) == 2 & ncol(tab) == 2) {
      a <- tab[2, 2]; b <- tab[2, 1]; cc <- tab[1, 2]; dd <- tab[1, 1]
      OR <- (a * dd) / (b * cc)
      se <- sqrt(1/a + 1/b + 1/cc + 1/dd)
      chi <- chisq.test(tab)
      cat(sprintf("  %s: OR=%.2f (95%% CI: %.2f-%.2f), p=%.4e\n",
                  var, OR, exp(log(OR)-1.96*se), exp(log(OR)+1.96*se),
                  chi$p.value))
      cat(sprintf("    Factor+: ASCVD %.1f%% (%d/%d)  Factor-: ASCVD %.1f%% (%d/%d)\n",
                  100*a/(a+b), a, a+b, 100*cc/(cc+dd), cc, cc+dd))
    }
  }
}

# BMI categories & ASCVD
cat("\n--- BMI CATEGORIES & ASCVD ---\n")
bmi_ascvd2 <- mb_data %>%
  filter(!is.na(bmi), bmi > 10, bmi < 60) %>%
  mutate(bmi_cat = cut(bmi, breaks = c(0, 25, 30, 60),
                       labels = c("<25", "25-30", ">30"))) %>%
  group_by(bmi_cat) %>%
  summarise(n = n(), events = sum(ascvd),
            rate = round(100 * events / n, 1), .groups = "drop")
print(bmi_ascvd2)

# TG tertiles & ASCVD
cat("\n--- TRIGLYCERIDE TERTILES & ASCVD ---\n")
tg_ascvd <- mb_data %>%
  filter(!is.na(tg)) %>%
  mutate(tg_tertile = ntile(tg, 3)) %>%
  group_by(tg_tertile) %>%
  summarise(
    tg_range = sprintf("[%.2f-%.2f]", min(tg), max(tg)),
    n = n(), events = sum(ascvd),
    rate = round(100 * events / n, 1), .groups = "drop")
print(tg_ascvd)

# Metabolic burden forest plot
cat("\n  Generating metabolic burden forest plot...\n")
mb_forest <- data.frame(
  Factor = c("Obesity (BMI>=30)", "Diabetes", "TG>=2.3 mmol/L"),
  stringsAsFactors = FALSE
)
mb_ors <- list()
for (i in seq_along(c("obese", "diabetes", "high_tg"))) {
  var <- c("obese", "diabetes", "high_tg")[i]
  sub <- mb_data %>% filter(!is.na(.data[[var]]), !is.na(ascvd))
  fit <- glm(ascvd ~ .data[[var]], data = sub, family = binomial)
  # Simple 2x2 OR
  tab <- table(sub[[var]], sub$ascvd)
  a <- tab[2,2]; b <- tab[2,1]; cc <- tab[1,2]; dd <- tab[1,1]
  mb_ors[[i]] <- c(OR = (a*dd)/(b*cc),
                   low = exp(log((a*dd)/(b*cc)) - 1.96*sqrt(1/a+1/b+1/cc+1/dd)),
                   high = exp(log((a*dd)/(b*cc)) + 1.96*sqrt(1/a+1/b+1/cc+1/dd)))
}
mb_forest$OR   <- sapply(mb_ors, function(x) x[1])
mb_forest$low  <- sapply(mb_ors, function(x) x[2])
mb_forest$high <- sapply(mb_ors, function(x) x[3])

p_mb <- ggplot(mb_forest, aes(x = OR, y = reorder(Factor, OR))) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = low, xmax = high), colour = "#E64B35", size = 0.8) +
  labs(x = "Odds Ratio (95% CI)", y = "",
       title = "Metabolic Burden Factors and ASCVD Risk in FH",
       subtitle = "Pooled across South Wales FH, UKB, and Wales FH cohorts") +
  theme_nature(base_size = 10) +
  theme(plot.margin = margin(10, 15, 10, 10))
safe_ggsave(file.path(data_dir, "Figure9_metabolic_burden_forest.pdf"),
            p_mb, width = 7, height = 4, dpi = 300)

################################################################################
# ── 18. ApoB/LDL DISCORDANCE ANALYSIS (0.3 CUTOFF) ─────────────────────────
################################################################################
cat("\n=== ApoB/LDL DISCORDANCE ANALYSIS ===\n")

# ApoB/LDL ratio available in South Wales FH + UKB only (Wales has no ApoB)
abl_data <- fh_pooled %>%
  filter(!is.na(apob_ldl), apob_ldl > 0, apob_ldl < 1, !is.na(ascvd))
cat(sprintf("ApoB/LDL analysis: N=%d (Dragon=%d, UKB=%d)\n",
            nrow(abl_data),
            sum(abl_data$cohort == "South Wales FH"),
            sum(abl_data$cohort == "UKB")))

# Distribution
cat("\n--- ApoB/LDL RATIO DISTRIBUTION ---\n")
cat(sprintf("  Overall: Mean=%.4f, SD=%.4f, Median=%.4f\n",
            mean(abl_data$apob_ldl), sd(abl_data$apob_ldl),
            median(abl_data$apob_ldl)))

pos <- abl_data %>% filter(ascvd == 1)
neg <- abl_data %>% filter(ascvd == 0)
cat(sprintf("  ASCVD+: Mean=%.4f, SD=%.4f (N=%d)\n",
            mean(pos$apob_ldl), sd(pos$apob_ldl), nrow(pos)))
cat(sprintf("  ASCVD-: Mean=%.4f, SD=%.4f (N=%d)\n",
            mean(neg$apob_ldl), sd(neg$apob_ldl), nrow(neg)))
tt_abl <- t.test(pos$apob_ldl, neg$apob_ldl)
cat(sprintf("  t-test: t=%.3f, p=%.4e\n", tt_abl$statistic, tt_abl$p.value))

# 0.3 CUTOFF
cat("\n--- 0.3 CUTOFF LINE ---\n")
abl_data$discordant <- ifelse(abl_data$apob_ldl >= 0.3, 1, 0)
below <- abl_data %>% filter(discordant == 0)
above <- abl_data %>% filter(discordant == 1)
cat(sprintf("  Below 0.3 (concordant): N=%d\n", nrow(below)))
cat(sprintf("  Above 0.3 (discordant): N=%d (%.1f%%)\n",
            nrow(above), 100 * nrow(above) / nrow(abl_data)))
cat(sprintf("  Below 0.3: ASCVD+ %d (%.1f%%), ASCVD- %d (%.1f%%)\n",
            sum(below$ascvd), 100*mean(below$ascvd),
            sum(below$ascvd == 0), 100*(1-mean(below$ascvd))))
cat(sprintf("  Above 0.3: ASCVD+ %d (%.1f%%), ASCVD- %d (%.1f%%)\n",
            sum(above$ascvd), 100*mean(above$ascvd),
            sum(above$ascvd == 0), 100*(1-mean(above$ascvd))))

# OR for 0.3 cutoff
tab03 <- table(abl_data$discordant, abl_data$ascvd)
a <- tab03[2,2]; b <- tab03[2,1]; cc <- tab03[1,2]; dd <- tab03[1,1]
or_03 <- (a * dd) / (b * cc)
se_03 <- sqrt(1/a + 1/b + 1/cc + 1/dd)
cat(sprintf("  OR (>=0.3 vs <0.3): %.2f (95%% CI: %.2f-%.2f)\n",
            or_03, exp(log(or_03) - 1.96*se_03),
            exp(log(or_03) + 1.96*se_03)))

# Detailed bins with ASCVD+ / ASCVD- percentages
cat("\n--- ApoB/LDL RATIO BINS ---\n")
abl_data$ratio_bin <- cut(abl_data$apob_ldl,
                           breaks = c(0, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 1.0))
bin_summary <- abl_data %>%
  group_by(ratio_bin) %>%
  summarise(
    N = n(),
    ASCVD_pos = sum(ascvd),
    ASCVD_neg = sum(ascvd == 0),
    pct_pos   = round(100 * ASCVD_pos / N, 1),
    pct_neg   = round(100 * ASCVD_neg / N, 1),
    .groups = "drop"
  )
print(bin_summary)

# Adjusted logistic regression models
cat("\n--- ADJUSTED ApoB/LDL MODELS ---\n")

# Model 1: Unadjusted
adj1 <- abl_data %>% filter(!is.na(tg), !is.na(diabetes))
fit1 <- glm(ascvd ~ apob_ldl, data = adj1, family = binomial)
cat(sprintf("Model 1 (unadjusted, N=%d): ApoB/LDL coef=%.4f, OR=%.3f, p=%.4e\n",
            nrow(adj1), coef(fit1)["apob_ldl"],
            exp(coef(fit1)["apob_ldl"]),
            summary(fit1)$coefficients["apob_ldl", 4]))

# Model 2: Adjusted for T2DM + TG
fit2 <- glm(ascvd ~ apob_ldl + diabetes + tg, data = adj1, family = binomial)
cat("Model 2 (+ T2DM + TG):\n")
for (v in c("apob_ldl", "diabetes", "tg")) {
  cat(sprintf("  %s: coef=%.4f, OR=%.3f (%.3f-%.3f), p=%.4e\n",
              v, coef(fit2)[v], exp(coef(fit2)[v]),
              exp(confint.default(fit2)[v, 1]),
              exp(confint.default(fit2)[v, 2]),
              summary(fit2)$coefficients[v, 4]))
}

# Model 3: + BMI (UKB only)
adj3 <- adj1 %>% filter(!is.na(bmi), bmi > 10, bmi < 60)
if (nrow(adj3) > 50) {
  fit3 <- glm(ascvd ~ apob_ldl + diabetes + tg + bmi, data = adj3, family = binomial)
  cat(sprintf("Model 3 (+ BMI, N=%d):\n", nrow(adj3)))
  for (v in c("apob_ldl", "diabetes", "tg", "bmi")) {
    cat(sprintf("  %s: coef=%.4f, OR=%.3f (%.3f-%.3f), p=%.4e\n",
                v, coef(fit3)[v], exp(coef(fit3)[v]),
                exp(confint.default(fit3)[v, 1]),
                exp(confint.default(fit3)[v, 2]),
                summary(fit3)$coefficients[v, 4]))
  }
}

# Effect of metabolic factors on ApoB/LDL ratio
cat("\n--- METABOLIC FACTORS AFFECTING ApoB/LDL RATIO ---\n")

# By diabetes
dm_abl <- abl_data %>% filter(!is.na(diabetes))
dm1 <- dm_abl %>% filter(diabetes == 1)
dm0 <- dm_abl %>% filter(diabetes == 0)
cat(sprintf("  DM+: ApoB/LDL mean=%.4f, %%disc=%.1f%% (N=%d)\n",
            mean(dm1$apob_ldl), 100*mean(dm1$apob_ldl >= 0.3), nrow(dm1)))
cat(sprintf("  DM-: ApoB/LDL mean=%.4f, %%disc=%.1f%% (N=%d)\n",
            mean(dm0$apob_ldl), 100*mean(dm0$apob_ldl >= 0.3), nrow(dm0)))
tt_dm <- t.test(dm1$apob_ldl, dm0$apob_ldl)
cat(sprintf("  t-test: p=%.4e\n", tt_dm$p.value))

# By TG tertile
tg_abl <- abl_data %>%
  filter(!is.na(tg)) %>%
  mutate(tg_tertile = ntile(tg, 3))
cat("\n  By TG tertile:\n")
for (tt in 1:3) {
  sub <- tg_abl %>% filter(tg_tertile == tt)
  cat(sprintf("    T%d: ApoB/LDL mean=%.4f, %%disc=%.1f%%, N=%d\n",
              tt, mean(sub$apob_ldl), 100*mean(sub$apob_ldl >= 0.3), nrow(sub)))
}

# By BMI category
bmi_abl <- abl_data %>%
  filter(!is.na(bmi), bmi > 10, bmi < 60) %>%
  mutate(bmi_cat = cut(bmi, breaks = c(0, 25, 30, 60),
                       labels = c("<25", "25-30", ">30")))
if (nrow(bmi_abl) > 30) {
  cat("\n  By BMI category:\n")
  for (cat_name in c("<25", "25-30", ">30")) {
    sub <- bmi_abl %>% filter(bmi_cat == cat_name)
    if (nrow(sub) > 0) {
      cat(sprintf("    BMI %s: ApoB/LDL mean=%.4f, %%disc=%.1f%%, N=%d\n",
                  cat_name, mean(sub$apob_ldl),
                  100*mean(sub$apob_ldl >= 0.3), nrow(sub)))
    }
  }
}

# Correlation ApoB/LDL vs TG
tg_cor <- cor.test(abl_data$apob_ldl[!is.na(abl_data$tg)],
                    abl_data$tg[!is.na(abl_data$tg)])
cat(sprintf("\n  Pearson r (ApoB/LDL vs TG): r=%.4f, p=%.4e\n",
            tg_cor$estimate, tg_cor$p.value))

# Discordance scatter plot with 0.3 cutoff line
cat("\n  Generating ApoB/LDL discordance plot...\n")
p_disc <- ggplot(abl_data, aes(x = apob_ldl, fill = factor(ascvd))) +
  geom_histogram(position = "dodge", bins = 40, alpha = 0.7) +
  geom_vline(xintercept = 0.3, linetype = "dashed", colour = "red", size = 1) +
  annotate("text", x = 0.30, y = Inf, label = "0.3 cutoff",
           vjust = 2, hjust = -0.1, colour = "red", fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c("0" = "#4DBBD5", "1" = "#E64B35"),
                    labels = c("ASCVD-", "ASCVD+")) +
  labs(x = "ApoB/LDL-C Ratio", y = "Count", fill = "ASCVD Status",
       title = "ApoB/LDL-C Discordance in Genetically Confirmed FH",
       subtitle = sprintf("N=%d | Discordant (>=0.3): %.1f%% | OR=%.2f",
                          nrow(abl_data), 100*mean(abl_data$discordant), or_03)) +
  theme_nature(base_size = 10)
safe_ggsave(file.path(data_dir, "Figure10_apob_ldl_discordance.pdf"),
            p_disc, width = 8, height = 5, dpi = 300)

# Save discordance summary table
disc_table <- data.frame(
  Metric = c("N_total", "N_concordant", "N_discordant", "pct_discordant",
             "ASCVD_rate_concordant", "ASCVD_rate_discordant",
             "OR_discordant_vs_concordant", "OR_CI_low", "OR_CI_high"),
  Value = c(nrow(abl_data), nrow(below), nrow(above),
            round(100 * nrow(above)/nrow(abl_data), 1),
            round(100 * mean(below$ascvd), 1),
            round(100 * mean(above$ascvd), 1),
            round(or_03, 2),
            round(exp(log(or_03) - 1.96*se_03), 2),
            round(exp(log(or_03) + 1.96*se_03), 2))
)
write.csv(disc_table, file.path(data_dir, "ApoB_LDL_discordance_summary.csv"),
          row.names = FALSE)
cat("  Saved: ApoB_LDL_discordance_summary.csv\n")

################################################################################
# ── 19. SUMMARY OUTPUT ─────────────────────────────────────────────────────
################################################################################

cat("\n")
cat("================================================================\n")
cat("                  CALON-X ANALYSIS COMPLETE                      \n")
cat("================================================================\n")
cat("\nFiles generated:\n")
cat("  Tables:\n")
cat("    - Table1_baseline_characteristics.csv\n")
cat("    - Table2a_gene_distribution.csv\n")
cat("    - Table2b_ascvd_by_gene.csv\n")
cat("    - Table3_index_vs_cascade.csv\n")
cat("    - Table4_FH_vs_NonFH_matched.csv\n")
cat("    - Table5_biexternal_validation.csv\n")
cat("    - CALON_X2_odds_ratios.csv\n")
cat("    - Subgroup_analysis_results.csv\n")
cat("    - Lancet_statistics_summary.csv\n")
cat("    - CALON_X3_PRS_coefficients.csv  [NEW: PRS integration]\n")
cat("    - ApoB_LDL_discordance_summary.csv  [NEW: 0.3 cutoff analysis]\n")
cat("  Figures:\n")
cat("    - Figure1_cohort_overview.pdf/png\n")
cat("    - Figure2_lipid_profiles.pdf/png\n")
cat("    - Figure3_lpa_distribution.pdf/png\n")
cat("    - Figure4_apob_vs_ldl.pdf/png\n")
cat("    - Figure5_forest_plot.pdf/png\n")
cat("    - Figure6_subgroup_forest.pdf/png\n")
cat("    - Figure7_event_rates.pdf/png\n")
cat("    - Figure8_clinical_dashboard.pdf/png  [Lancet-caliber 6-panel]\n")
cat("    - Figure9_metabolic_burden_forest.pdf/png  [NEW: metabolic burden ORs]\n")
cat("    - Figure10_apob_ldl_discordance.pdf/png  [NEW: 0.3 cutoff]\n")
cat("\nDone.\n")

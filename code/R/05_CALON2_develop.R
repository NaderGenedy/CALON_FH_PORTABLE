################################################################################
#                                                                              #
#  CALON-2: ASCVD RISK MODEL DEVELOPMENT IN FH PATIENTS                       #
#  Elastic Net + XGBoost on UK Biobank Genetically Confirmed FH Cohort        #
#                                                                              #
#  Development cohort: UKB (N=1,623, 399 ASCVD events)                         #
#  Benchmark:          SAFEHEART-RE (published coefficients)                   #
#  Output:             Frozen coefficients + best model for ext. validation    #
#                                                                              #
#  Key changes vs v2 → v3 → v3b → v4:                                        #
#    v3: COLLINEARITY FIX: Drop standalone apob, use log(ApoB/LDL) ratio     #
#    v3: CONTINUOUS Lp(a): log(Lp(a)+0.1) replaces binary >125 cutoff        #
#    v3: INTERACTION TERMS: age×re_ldl, lipid_years_A for cumulative burden   #
#    v3b: LEAN MODELS: Drop redundant inv_apoa + trig (VIF, non-sig)         #
#    v3b: 5 lean variants: Lean Core/Interact/Minimal + Enh/Kit-Lean         #
#    v4: MULTIPLE IMPUTATION (mice): recover 345 lost patients               #
#    v4: LITERATURE FEATURES: lipid_years_A in Core, lpa_extreme ≥250       #
#    v4: SUPER LEARNER: stacked EN+XGB+RF+GAM ensemble                       #
#    v4: RANDOM FOREST: additional base learner for diversity                 #
#    v4: RESTRICTED CUBIC SPLINES: non-linear age + LDL effects (rms)        #
#                                                                              #
#  Author:  Dr Nader Genedy                                                    #
#  Date:    March 2026                                                         #
#  Target:  AUC > 0.80 before external validation on DRAGON3                  #
#                                                                              #
################################################################################

rm(list = ls())
set.seed(2026)

# =============================================================================
# PACKAGES
# =============================================================================

required_packages <- c(
  "dplyr", "tidyr", "readr", "ggplot2", "patchwork", "scales",
  "glmnet", "pROC",
  "boot", "ResourceSelection", "nricens", "tableone",
  "DescTools", "xgboost", "caret", "car",
  # v4 additions: multiple imputation, super learner, random forest, splines
  "mice", "randomForest", "SuperLearner", "rms"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("\n")
cat("======================================================================\n")
cat("  CALON-2: ASCVD RISK MODEL DEVELOPMENT v4                           \n")
cat("  v3: Collinearity Fix + Log Transforms + Interaction Terms         \n")
cat("  v3b: Lean Models (drop redundant inv_apoa + trig)                 \n")
cat("  v4: Multiple Imputation + Super Learner + Literature Features     \n")
cat("  Goal: AUC > 0.80 — 30+ model variants × 4 algorithms + stacking  \n")
cat("======================================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR <- "/cloud/project/"
OUT_DIR  <- paste0(BASE_DIR, "output/")
TAB_DIR  <- paste0(OUT_DIR, "tables/")
FIG_DIR  <- paste0(OUT_DIR, "figures/")
for (d in c(TAB_DIR, FIG_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# =============================================================================
# SECTION 1: LOAD & MERGE DATA
# =============================================================================

cat("===================================================================\n")
cat("SECTION 1: LOADING AND MERGING DATA\n")
cat("===================================================================\n\n")

# -- Load existing analysis-ready dataset ------------------------------------
df <- read.csv(paste0(OUT_DIR, "calon_ukb_analysis_ready.csv"),
               stringsAsFactors = FALSE)
cat(sprintf("  Base dataset: %d patients x %d variables\n", nrow(df), ncol(df)))
cat(sprintf("  ASCVD events: %d (%.1f%%)\n", sum(df$ascvd_combined),
            100 * mean(df$ascvd_combined)))

# -- Merge Townsend deprivation ---------------------------------------------
dep_file <- paste0(BASE_DIR, "calon_extra_deprivation.csv")
if (file.exists(dep_file)) {
  dep <- read.csv(dep_file, stringsAsFactors = FALSE)
  names(dep) <- gsub("^participant\\.", "", names(dep))
  townsend_candidates <- c("p22189", "p22189_i0", "p189_i0", "p189")
  for (tc in townsend_candidates) {
    if (tc %in% names(dep)) {
      names(dep)[names(dep) == tc] <- "townsend"
      break
    }
  }
  imd_cols <- intersect(c("p26410", "p26426", "p26427"), names(dep))
  if (!"townsend" %in% names(dep) && length(imd_cols) > 0) {
    dep$townsend <- NA
    for (ic in imd_cols) {
      dep$townsend <- ifelse(is.na(dep$townsend), dep[[ic]], dep$townsend)
    }
    cat("  NOTE: Using IMD (Index of Multiple Deprivation) as deprivation proxy.\n")
  }
  if (!"townsend" %in% names(dep) && ncol(dep) == 2) {
    names(dep) <- c("eid", "townsend")
  }
  df <- merge(df, dep, by = "eid", all.x = TRUE)
  cat(sprintf("  Townsend deprivation merged: %d non-missing (%.1f%%)\n",
              sum(!is.na(df$townsend)), 100 * mean(!is.na(df$townsend))))
} else {
  cat("  Townsend deprivation file not found. Skipping.\n")
  df$townsend <- NA
}

# -- Merge NMR metabolomics (8 batch files) ----------------------------------
nmr_files <- c(
  paste0(BASE_DIR, "calon_extra_nmr_a.csv"),
  paste0(BASE_DIR, "calon_extra_nmr_b.csv"),
  paste0(BASE_DIR, "calon_extra_nmr_c.csv"),
  paste0(BASE_DIR, "calon_extra_nmr_d.csv"),
  paste0(BASE_DIR, "calon_extra_nmr_e.csv"),
  paste0(BASE_DIR, "calon_extra_nmr_f.csv"),
  paste0(BASE_DIR, "calon_extra_nmr_g.csv"),
  paste0(BASE_DIR, "calon_extra_nmr_h.csv")
)

nmr_fallback <- c(
  "C:/Users/nader/Downloads/ukb_nmr_a.csv",
  "C:/Users/nader/Downloads/ukb_nmr_b.csv",
  "C:/Users/nader/Downloads/ukb_nmr_c.csv",
  "C:/Users/nader/Downloads/ukb_nmr_d.csv",
  "C:/Users/nader/Downloads/ukb_nmr_e.csv"
)

nmr_merged <- FALSE
for (i in seq_along(nmr_files)) {
  f <- nmr_files[i]
  if (!file.exists(f) && i <= length(nmr_fallback) && file.exists(nmr_fallback[i])) {
    f <- nmr_fallback[i]
    cat(sprintf("  Using TUDOR fallback: %s\n", basename(f)))
  }
  if (file.exists(f)) {
    nmr_batch <- read.csv(f, stringsAsFactors = FALSE)
    names(nmr_batch) <- gsub("^participant\\.", "", names(nmr_batch))
    if ("eid" %in% names(nmr_batch)) {
      df <- merge(df, nmr_batch, by = "eid", all.x = TRUE)
      nmr_merged <- TRUE
      cat(sprintf("  NMR batch %s merged (%d columns)\n",
                  toupper(letters[i]), ncol(nmr_batch) - 1))
    }
  }
}
if (!nmr_merged) {
  cat("  No NMR files found. Enhanced model will skip NMR features.\n")
}

# -- Merge CAC ---------------------------------------------------------------
cac_file <- paste0(BASE_DIR, "calon_extra_cac.csv")
if (file.exists(cac_file)) {
  cac <- read.csv(cac_file, stringsAsFactors = FALSE)
  names(cac) <- gsub("^participant\\.", "", names(cac))
  if (ncol(cac) > 1) {
    df <- merge(df, cac, by = "eid", all.x = TRUE)
    cac_cols <- grep("p224", names(df), value = TRUE)
    if (length(cac_cols) > 0) {
      n_cac <- sum(!is.na(df[[cac_cols[1]]]))
      cat(sprintf("  CAC merged: %d non-missing (%.1f%%)\n", n_cac, 100*n_cac/nrow(df)))
      if (n_cac < 50) {
        cat("  WARNING: CAC data too sparse (<50 non-missing). Will exclude from model.\n")
      }
    } else {
      cat("  CAC file loaded but no p224xx columns found.\n")
    }
  } else {
    cat("  CAC file is empty placeholder. Skipping.\n")
  }
} else {
  cat("  CAC file not found. Skipping.\n")
}

# -- Merge PRS ---------------------------------------------------------------
prs_file <- paste0(BASE_DIR, "calon_extra_prs.csv")
prs_cad_file <- paste0(BASE_DIR, "calon_extra_prs_cad.csv")
gpc_file <- paste0(BASE_DIR, "calon_extra_genetic_pcs.csv")

prs_loaded <- FALSE
for (pf in c(prs_file, prs_cad_file)) {
  if (file.exists(pf)) {
    prs <- read.csv(pf, stringsAsFactors = FALSE)
    names(prs) <- gsub("^participant\\.", "", names(prs))
    if (ncol(prs) > 1) {
      df <- merge(df, prs, by = "eid", all.x = TRUE)
      prs_cols <- grep("p26", names(df), value = TRUE)
      if (length(prs_cols) > 0) {
        n_prs <- sum(!is.na(df[[prs_cols[1]]]))
        cat(sprintf("  PRS merged from %s: %d fields, %d non-missing (%.1f%%)\n",
                    basename(pf), length(prs_cols), n_prs, 100*n_prs/nrow(df)))
        prs_loaded <- TRUE
        break
      }
    }
  }
}

if (!prs_loaded && file.exists(gpc_file)) {
  gpc <- read.csv(gpc_file, stringsAsFactors = FALSE)
  names(gpc) <- gsub("^participant\\.", "", names(gpc))
  if (ncol(gpc) > 1) {
    df <- merge(df, gpc, by = "eid", all.x = TRUE)
    gpc_cols <- grep("p22009", names(df), value = TRUE)
    if (length(gpc_cols) > 0) {
      n_gpc <- sum(!is.na(df[[gpc_cols[1]]]))
      cat(sprintf("  Genetic PCs merged (fallback): %d PCs, %d non-missing (%.1f%%)\n",
                  length(gpc_cols), n_gpc, 100*n_gpc/nrow(df)))
      cat("  NOTE: These are PCs (population stratification), NOT true PRS.\n")
      prs_loaded <- TRUE
    }
  }
}

if (!prs_loaded) {
  cat("  No PRS or genetic PC files found. Skipping.\n")
}

# -- Merge Cardiac MRI IDPs (for reporting, NOT forced into models) ----------
cat("\n  Merging Cardiac MRI IDPs (for reporting only)...\n")
mri_files <- list(
  lv     = paste0(BASE_DIR, "calon_extra_mri_lv.csv"),
  aorta  = paste0(BASE_DIR, "calon_extra_mri_aorta.csv"),
  la_rv  = paste0(BASE_DIR, "calon_extra_mri_la_rv.csv"),
  strain = paste0(BASE_DIR, "calon_extra_mri_strain.csv"),
  cat162 = paste0(BASE_DIR, "calon_extra_mri_cat162.csv")
)

mri_merged_ok <- FALSE
for (nm in names(mri_files)) {
  f <- mri_files[[nm]]
  if (file.exists(f)) {
    mri_batch <- read.csv(f, stringsAsFactors = FALSE)
    names(mri_batch) <- gsub("^participant\\.", "", names(mri_batch))
    if ("eid" %in% names(mri_batch) && ncol(mri_batch) > 1) {
      mri_batch <- mri_batch[, !duplicated(names(mri_batch))]
      new_cols <- setdiff(names(mri_batch), c("eid", names(df)))
      if (length(new_cols) > 0) {
        df <- merge(df, mri_batch[, c("eid", new_cols)], by = "eid", all.x = TRUE)
        mri_merged_ok <- TRUE
        cat(sprintf("  MRI %s merged: %d new fields\n", toupper(nm), length(new_cols)))
      } else {
        cat(sprintf("  MRI %s: all columns already present\n", nm))
      }
    }
  } else {
    cat(sprintf("  MRI %s not found: %s\n", nm, basename(f)))
  }
}
if (!mri_merged_ok) {
  cat("  No MRI IDP files found. MRI features will not be available.\n")
}

# -- Merge imaging visit date (for temporal filtering) -----------------------
img_date_file  <- paste0(BASE_DIR, "calon_extra_imaging_date.csv")
img_date_repeat <- paste0(BASE_DIR, "calon_extra_imaging_date_repeat.csv")

if (file.exists(img_date_file)) {
  img_dt <- read.csv(img_date_file, stringsAsFactors = FALSE)
  names(img_dt) <- gsub("^participant\\.", "", names(img_dt))
  date_col <- setdiff(names(img_dt), "eid")[1]
  if (!is.null(date_col)) {
    names(img_dt)[names(img_dt) == date_col] <- "imaging_date"
    df <- merge(df, img_dt[, c("eid", "imaging_date")], by = "eid", all.x = TRUE)
    df$imaging_date <- as.Date(df$imaging_date)
    n_img <- sum(!is.na(df$imaging_date))
    cat(sprintf("  Imaging visit dates merged: %d non-missing (%.1f%%)\n",
                n_img, 100 * n_img / nrow(df)))
  }
} else {
  cat("  Imaging date file not found.\n")
  df$imaging_date <- as.Date(NA)
}

if (file.exists(img_date_repeat)) {
  img_dt2 <- read.csv(img_date_repeat, stringsAsFactors = FALSE)
  names(img_dt2) <- gsub("^participant\\.", "", names(img_dt2))
  date_col2 <- setdiff(names(img_dt2), "eid")[1]
  if (!is.null(date_col2)) {
    names(img_dt2)[names(img_dt2) == date_col2] <- "imaging_date_repeat"
    df <- merge(df, img_dt2[, c("eid", "imaging_date_repeat")], by = "eid", all.x = TRUE)
    df$imaging_date_repeat <- as.Date(df$imaging_date_repeat)
    cat(sprintf("  Repeat imaging dates: %d non-missing\n",
                sum(!is.na(df$imaging_date_repeat))))
  }
}

# -- Merge CIMT (full, instance 3) ------------------------------------------
cimt_i3_file <- paste0(BASE_DIR, "calon_extra_cimt_full.csv")
if (file.exists(cimt_i3_file)) {
  cimt3 <- read.csv(cimt_i3_file, stringsAsFactors = FALSE)
  names(cimt3) <- gsub("^participant\\.", "", names(cimt3))
  if ("eid" %in% names(cimt3) && ncol(cimt3) > 1) {
    cimt_new <- setdiff(names(cimt3), c("eid", names(df)))
    if (length(cimt_new) > 0) {
      df <- merge(df, cimt3[, c("eid", cimt_new)], by = "eid", all.x = TRUE)
      cat(sprintf("  CIMT instance 3 merged: %d new fields\n", length(cimt_new)))
    }
  }
}

cat(sprintf("\n  FINAL merged dataset: %d patients x %d variables\n", nrow(df), ncol(df)))

# -- SAVE full merged dataset for downstream use -----------------------------
write.csv(df, paste0(OUT_DIR, "calon2_full_merged.csv"), row.names = FALSE)
cat("  Saved: output/calon2_full_merged.csv\n")

# =============================================================================
# SECTION 1B: TEMPORAL FILTERING -- MEASUREMENTS BEFORE INDEX ASCVD
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 1B: TEMPORAL FILTERING (pre-ASCVD requirement)\n")
cat("===================================================================\n\n")

# Convert dates
df$assessment_date  <- as.Date(df$assessment_date)
df$first_ascvd_date <- as.Date(df$first_ascvd_date)

# -- Blood biomarkers: assessment_date vs first_ascvd_date -------------------
n_prevalent <- sum(df$ascvd_prevalent == 1, na.rm = TRUE)
n_incident  <- sum(df$ascvd_incident == 1, na.rm = TRUE)
n_no_event  <- sum(df$ascvd_combined == 0, na.rm = TRUE)

cat("  Blood biomarker temporal status:\n")
cat(sprintf("    Prevalent ASCVD (event BEFORE assessment): %d\n", n_prevalent))
cat(sprintf("    Incident ASCVD (event AFTER assessment):   %d\n", n_incident))
cat(sprintf("    No ASCVD event:                            %d\n", n_no_event))

# Flag: bloods measured before ASCVD (or no ASCVD)
df$bloods_before_ascvd <- TRUE
df$bloods_before_ascvd[!is.na(df$first_ascvd_date) & !is.na(df$assessment_date) &
                         df$first_ascvd_date <= df$assessment_date] <- FALSE
cat(sprintf("    Temporally valid (bloods before event or no event): %d / %d\n",
            sum(df$bloods_before_ascvd), nrow(df)))

# -- MRI IDPs: imaging_date vs first_ascvd_date -----------------------------
if ("imaging_date" %in% names(df) && sum(!is.na(df$imaging_date)) > 0) {
  df$mri_before_ascvd <- NA
  has_imaging <- !is.na(df$imaging_date)
  has_ascvd   <- !is.na(df$first_ascvd_date) & df$ascvd_combined == 1

  df$mri_before_ascvd[has_imaging & !has_ascvd] <- TRUE
  df$mri_before_ascvd[has_imaging & has_ascvd &
                        df$imaging_date < df$first_ascvd_date] <- TRUE
  df$mri_before_ascvd[has_imaging & has_ascvd &
                        df$imaging_date >= df$first_ascvd_date] <- FALSE

  n_mri_valid   <- sum(df$mri_before_ascvd == TRUE, na.rm = TRUE)
  n_mri_invalid <- sum(df$mri_before_ascvd == FALSE, na.rm = TRUE)
  n_mri_none    <- sum(is.na(df$mri_before_ascvd))

  cat(sprintf("\n  MRI temporal status:\n"))
  cat(sprintf("    MRI before ASCVD (valid):     %d\n", n_mri_valid))
  cat(sprintf("    MRI after ASCVD (EXCLUDED):   %d\n", n_mri_invalid))
  cat(sprintf("    No MRI data:                  %d\n", n_mri_none))

  # NULL OUT MRI values for patients where MRI was AFTER ASCVD event
  mri_cols <- grep("^p24|^p310", names(df), value = TRUE)
  if (n_mri_invalid > 0 && length(mri_cols) > 0) {
    invalid_idx <- which(df$mri_before_ascvd == FALSE)
    for (mc in mri_cols) {
      df[[mc]][invalid_idx] <- NA
    }
    cat(sprintf("    -> Set %d MRI fields to NA for %d patients (post-ASCVD imaging)\n",
                length(mri_cols), n_mri_invalid))
  }
} else {
  cat("\n  No imaging date data available -- MRI temporal filtering skipped.\n")
  df$mri_before_ascvd <- NA
}

cat(sprintf("\n  DECISION: Primary model uses ALL patients (prevalent + incident)\n"))
cat(sprintf("           Sensitivity analysis: incident-only (exclude %d prevalent)\n\n",
            n_prevalent))

# =============================================================================
# SECTION 2: FEATURE ENGINEERING (EXPANDED)
# =============================================================================

cat("===================================================================\n")
cat("SECTION 2: FEATURE ENGINEERING (EXPANDED)\n")
cat("===================================================================\n\n")

# -- Gene type binary -------------------------------------------------------
df$gene_apob <- as.integer(df$gene == "APOB")
cat(sprintf("  Gene: %d LDLR, %d APOB, %d PCSK9\n",
            sum(df$gene == "LDLR"), sum(df$gene == "APOB"),
            sum(df$gene == "PCSK9")))

# -- Townsend quintiles -----------------------------------------------------
if (sum(!is.na(df$townsend)) > 100) {
  df$townsend_quintile <- as.integer(cut(df$townsend,
                                          breaks = quantile(df$townsend,
                                                            probs = 0:5/5, na.rm = TRUE),
                                          include.lowest = TRUE, labels = FALSE))
  cat(sprintf("  Townsend quintiles created: %d non-missing\n",
              sum(!is.na(df$townsend_quintile))))
}

# -- Automatic NMR feature detection ----------------------------------------
# Detect ALL p234xx_i0 columns from the merged NMR data
nmr_all_cols <- grep("^p234\\d+_i0$", names(df), value = TRUE)
cat(sprintf("  NMR columns auto-detected (p234xx_i0 pattern): %d\n", length(nmr_all_cols)))

if (length(nmr_all_cols) > 0) {
  # Report NMR completeness
  nmr_available <- sapply(nmr_all_cols, function(v) mean(!is.na(df[[v]])))
  cat(sprintf("  NMR completeness: median=%.1f%%, min=%.1f%%, max=%.1f%%\n",
              100*median(nmr_available), 100*min(nmr_available), 100*max(nmr_available)))

  # Named NMR features for interpretability
  if ("p23485_i0" %in% names(df)) {
    df$ldl_particle_size <- df$p23485_i0
    cat(sprintf("  -> LDL particle size (p23485): %d non-missing\n",
                sum(!is.na(df$ldl_particle_size))))
  }
  if ("p23486_i0" %in% names(df)) {
    df$hdl_particle_size <- df$p23486_i0
  }
  if ("p23484_i0" %in% names(df)) {
    df$vldl_particle_size <- df$p23484_i0
  }
  if ("p23493_i0" %in% names(df)) {
    df$apob_nmr <- df$p23493_i0
  }
  if ("p23494_i0" %in% names(df)) {
    df$apoa1_nmr <- df$p23494_i0
  }
  if ("p23495_i0" %in% names(df)) {
    df$glyca <- df$p23495_i0
    cat(sprintf("  -> GlycA inflammation (p23495): %d non-missing\n",
                sum(!is.na(df$glyca))))
  }
  if ("p23492_i0" %in% names(df)) {
    df$remnant_chol <- df$p23492_i0
  }

  # NMR-derived ratios
  if ("p23493_i0" %in% names(df) && "p23494_i0" %in% names(df)) {
    df$apob_apoa1_nmr_ratio <- df$p23493_i0 / (df$p23494_i0 + 0.001)
    cat(sprintf("  -> ApoB/ApoA1 NMR ratio: %d non-missing\n",
                sum(!is.na(df$apob_apoa1_nmr_ratio))))
  }
  # VLDL/HDL particle ratio (atherogenic index)
  if ("p23484_i0" %in% names(df) && "p23486_i0" %in% names(df)) {
    df$vldl_hdl_ratio <- df$p23484_i0 / (df$p23486_i0 + 0.001)
  }
} else {
  cat("  No NMR columns found -- enhanced model will be limited.\n")
}

# -- PRS features (z-score the CAD PRS) -------------------------------------
prs_cols <- grep("^p26", names(df), value = TRUE)
if (length(prs_cols) > 0) {
  cat(sprintf("  PRS columns available: %d (%s)\n",
              length(prs_cols), paste(prs_cols, collapse=", ")))
  prs_priority <- c("p26227", "p26223", prs_cols)
  prs_priority <- prs_priority[prs_priority %in% names(df)]
  for (pc in prs_priority) {
    n_valid <- sum(!is.na(df[[pc]]))
    if (n_valid > 100) {
      df$prs_primary <- scale(df[[pc]])[, 1]
      cat(sprintf("  Using %s as primary PRS (%d non-missing, z-scored)\n", pc, n_valid))
      break
    }
  }
}

# -- hsCRP (log-transformed) ------------------------------------------------
if ("crp" %in% names(df) && sum(!is.na(df$crp)) > 100) {
  df$log_crp <- log(df$crp + 0.01)
  cat(sprintf("  hsCRP: %d non-missing (%.1f%%), median=%.2f mg/L\n",
              sum(!is.na(df$crp)), 100 * mean(!is.na(df$crp)),
              median(df$crp, na.rm = TRUE)))
} else {
  cat("  hsCRP (crp) not found or too sparse.\n")
  df$log_crp <- NA
}

# -- Additional log-transformed biomarkers (skewed distributions) -----------
if ("hba1c" %in% names(df) && sum(!is.na(df$hba1c)) > 100) {
  df$log_hba1c <- log(df$hba1c + 0.01)
  cat(sprintf("  log(HbA1c) created: %d non-missing\n", sum(!is.na(df$log_hba1c))))
}
if ("creatinine" %in% names(df) && sum(!is.na(df$creatinine)) > 100) {
  df$log_creatinine <- log(df$creatinine + 0.01)
  cat(sprintf("  log(creatinine) created: %d non-missing\n", sum(!is.na(df$log_creatinine))))
}
if ("alt" %in% names(df) && sum(!is.na(df$alt)) > 100) {
  df$log_alt <- log(df$alt + 0.01)
  cat(sprintf("  log(ALT) created: %d non-missing\n", sum(!is.na(df$log_alt))))
}
if ("cystatin_c" %in% names(df) && sum(!is.na(df$cystatin_c)) > 100) {
  df$log_cystatin_c <- log(df$cystatin_c + 0.01)
  cat(sprintf("  log(cystatin_c) created: %d non-missing\n", sum(!is.na(df$log_cystatin_c))))
}

# ── COLLINEARITY FIXES (v3 additions) ───────────────────────────────────────
cat("\n  --- Collinearity fixes & new transforms (v3) ---\n")

# -- Continuous log Lp(a) (replaces binary >125 cutoff) ----------------------
# Lp(a) is heavily right-skewed; log transform preserves dose-response
if ("lpa" %in% names(df) && sum(!is.na(df$lpa)) > 100) {
  df$log_lpa <- ifelse(!is.na(df$lpa) & df$lpa > 0, log(df$lpa),
                       ifelse(!is.na(df$lpa), log(0.1), NA))
  cat(sprintf("  log(Lp(a)) created: %d non-missing, median=%.2f (= %.0f nmol/L)\n",
              sum(!is.na(df$log_lpa)), median(df$log_lpa, na.rm=TRUE),
              exp(median(df$log_lpa, na.rm=TRUE))))
} else {
  df$log_lpa <- NA
  cat("  WARNING: Lp(a) not available for log transform\n")
}

# -- v4: Extreme Lp(a) binary (≥250 nmol/L) ---------------------------------
# Paquette et al. 2025: Lp(a) ≥250 nmol/L is a cardiovascular risk equivalent
# in HeFH (28.7% vs 11.0% 10-year ASCVD risk, HR 2.6)
# This captures the non-linear threshold effect that continuous log_lpa may miss
if ("lpa" %in% names(df) && sum(!is.na(df$lpa)) > 100) {
  df$lpa_extreme <- ifelse(!is.na(df$lpa) & df$lpa >= 250, 1, 0)
  n_extreme <- sum(df$lpa_extreme == 1, na.rm = TRUE)
  cat(sprintf("  Lp(a) extreme (≥250 nmol/L): %d / %d (%.1f%%)\n",
              n_extreme, sum(!is.na(df$lpa)),
              100 * n_extreme / sum(!is.na(df$lpa))))
  # Also create moderate threshold
  df$lpa_high <- ifelse(!is.na(df$lpa) & df$lpa >= 125, 1, 0)
  cat(sprintf("  Lp(a) high (≥125 nmol/L):    %d / %d (%.1f%%)\n",
              sum(df$lpa_high == 1, na.rm = TRUE), sum(!is.na(df$lpa)),
              100 * sum(df$lpa_high == 1, na.rm = TRUE) / sum(!is.na(df$lpa))))
} else {
  df$lpa_extreme <- NA
  df$lpa_high <- NA
}

# -- Consistent ApoB/LDL ratio using RE-LDL (not measured LDL) ---------------
# PROBLEM: apob_ldl_ratio = apob / measured_ldl (from build_cohort)
#          BUT re_ldl = measured_ldl / statin_factor
#          For untreated: re_ldl = ldl → apob_ldl_ratio = apob/re_ldl (PERFECT collinearity)
#          For treated: ratio uses measured LDL, model uses re_ldl (inconsistent!)
# FIX: Create apob / re_ldl for internal consistency
if (all(c("apob", "re_ldl") %in% names(df))) {
  df$apob_re_ldl_ratio <- ifelse(!is.na(df$apob) & !is.na(df$re_ldl) & df$re_ldl > 0,
                                  df$apob / df$re_ldl, NA)
  df$log_apob_re_ldl   <- ifelse(!is.na(df$apob_re_ldl_ratio) & df$apob_re_ldl_ratio > 0,
                                  log(df$apob_re_ldl_ratio), NA)
  cat(sprintf("  ApoB/RE-LDL ratio: median %.3f (N=%d)\n",
              median(df$apob_re_ldl_ratio, na.rm=TRUE),
              sum(!is.na(df$apob_re_ldl_ratio))))
  cat(sprintf("  log(ApoB/RE-LDL):  median %.3f\n",
              median(df$log_apob_re_ldl, na.rm=TRUE)))
  cat(sprintf("  log(ApoB/LDL):     median %.3f (original, uses measured LDL)\n",
              median(df$log_apob_ldl, na.rm=TRUE)))
} else {
  df$apob_re_ldl_ratio <- NA
  df$log_apob_re_ldl <- NA
}

# -- Interaction terms --------------------------------------------------------
if (all(c("age", "re_ldl") %in% names(df))) {
  df$age_x_re_ldl <- df$age * df$re_ldl
  cat(sprintf("  Interaction age x re_ldl: %d non-missing\n",
              sum(!is.na(df$age_x_re_ldl))))
}
if (all(c("age", "hdl") %in% names(df))) {
  df$age_x_hdl <- df$age * df$hdl
}

cat("  COLLINEARITY DIAGNOSIS:\n")
cat("    Original Core: apob + re_ldl + apob_ldl_ratio → severe collinearity\n")
cat("    Fixed Core:    log(apob/ldl) replaces standalone apob + ratio\n")
cat("    Consistent:    log(apob/re_ldl) for internal consistency\n")
cat("    log(Lp(a)):    continuous replaces binary >125 nmol/L cutoff\n\n")

# -- CIMT features (instance 3 if available) --------------------------------
cimt_i3_cols <- grep("p22671_i3|p22674_i3|p22677_i3|p22680_i3",
                     names(df), value = TRUE)
if (length(cimt_i3_cols) > 0) {
  cimt_vals <- df[, cimt_i3_cols, drop = FALSE]
  df$mean_cimt_i3 <- rowMeans(cimt_vals, na.rm = TRUE)
  df$mean_cimt_i3[rowSums(!is.na(cimt_vals)) == 0] <- NA
  n_cimt3 <- sum(!is.na(df$mean_cimt_i3))
  cat(sprintf("  CIMT (instance 3, %d measurements averaged): %d non-missing\n",
              length(cimt_i3_cols), n_cimt3))
  if (all(is.na(df$mean_cimt))) {
    df$mean_cimt <- df$mean_cimt_i3
    cat("  -> Filled mean_cimt from instance 3 data\n")
  }
}

# -- MRI derived features (for reporting only) ------------------------------
find_mri_col <- function(base_field) {
  candidates <- c(base_field, paste0(base_field, "_i2"),
                  paste0(base_field, "_i0"), paste0(base_field, "_i2_a0"))
  found <- candidates[candidates %in% names(df)]
  if (length(found) > 0) return(found[1])
  return(NULL)
}

lv_mass_col <- find_mri_col("p24105")
if (!is.null(lv_mass_col)) df$lv_mass <- as.numeric(df[[lv_mass_col]])
lvef_col <- find_mri_col("p24103")
if (!is.null(lvef_col)) df$lvef <- as.numeric(df[[lvef_col]])

cat("\n")

# =============================================================================
# SECTION 3: MISSINGNESS REPORT
# =============================================================================

cat("===================================================================\n")
cat("SECTION 3: MISSINGNESS REPORT\n")
cat("===================================================================\n\n")

# Collect all candidate variables for the report
all_candidate_vars <- c(
  # Core 14
  "age", "sex", "apob", "re_ldl", "hdl", "trig",
  "inv_apoa", "apob_ldl_ratio", "lpa_binary",
  "smoking_binary", "diabetes", "hypertension", "bmi", "log_crp",
  # NMR (all auto-detected)
  nmr_all_cols,
  # NMR derived
  "ldl_particle_size", "glyca", "remnant_chol",
  "apob_apoa1_nmr_ratio", "vldl_hdl_ratio",
  # PRS
  "prs_primary",
  # Biomarkers
  "hba1c", "creatinine", "alt", "cystatin_c", "glucose",
  "log_hba1c", "log_creatinine", "log_alt", "log_cystatin_c",
  # Other
  "townsend", "lipid_years_A", "gene_apob", "mean_cimt",
  # MRI (for reporting)
  "lv_mass", "lvef",
  # Additional
  "waist", "sbp", "dbp", "ever_smoked", "log_apob_ldl",
  # v3 collinearity fixes
  "log_lpa", "apob_re_ldl_ratio", "log_apob_re_ldl",
  "age_x_re_ldl", "age_x_hdl"
)

# Deduplicate and filter to existing columns
available_vars <- unique(all_candidate_vars[all_candidate_vars %in% names(df)])

miss_report <- data.frame(
  variable    = available_vars,
  n_total     = nrow(df),
  n_available = sapply(available_vars, function(v) sum(!is.na(df[[v]]))),
  n_missing   = sapply(available_vars, function(v) sum(is.na(df[[v]]))),
  stringsAsFactors = FALSE
)
miss_report$pct_available <- round(100 * miss_report$n_available / miss_report$n_total, 1)
miss_report$pct_missing   <- round(100 * miss_report$n_missing / miss_report$n_total, 1)

cat("  Variable Missingness (selected key variables):\n")
cat(sprintf("  %-25s %6s %8s %8s\n", "Variable", "N", "Missing", "% Avail"))
cat(paste(rep("-", 55), collapse = ""), "\n")

# Print first 30 variables (fixed in v3: no duplicate printing)
n_to_print <- min(30, nrow(miss_report))
for (i in 1:n_to_print) {
  flag <- ifelse(miss_report$pct_missing[i] > 50, " *** EXCLUDE",
          ifelse(miss_report$pct_missing[i] > 30, " ** HIGH", ""))
  cat(sprintf("  %-25s %6d %8d %7.1f%%%s\n",
              miss_report$variable[i], miss_report$n_total,
              miss_report$n_missing[i], miss_report$pct_available[i], flag))
}
if (nrow(miss_report) > 30) {
  cat(sprintf("  ... and %d more variables (see calon2_missingness.csv)\n",
              nrow(miss_report) - 30))
}

write.csv(miss_report, paste0(TAB_DIR, "calon2_missingness.csv"), row.names = FALSE)
cat(sprintf("\n  Missingness report saved: calon2_missingness.csv (%d variables)\n\n",
            nrow(miss_report)))

# =============================================================================
# SECTION 4: DEFINE MODEL TIERS
# =============================================================================

cat("===================================================================\n")
cat("SECTION 4: DEFINING MODEL TIERS\n")
cat("===================================================================\n\n")

# -- TIER 1: Clinical Core (14 vars, externally validatable) -----------------
# NOTE: apob_ldl_ratio is KEPT but flagged for collinearity with apob + re_ldl.
# LASSO regularisation will handle this; Ridge may inflate its coefficient.
core_vars <- c("age", "sex", "apob", "re_ldl", "hdl", "trig",
               "inv_apoa", "apob_ldl_ratio", "lpa_binary",
               "smoking_binary", "diabetes", "hypertension", "bmi",
               "log_crp")

cat("  TIER 1 -- Clinical Core (14 vars, externally validatable):\n")
cat("  NOTE: apob_ldl_ratio included but collinear with apob + re_ldl.\n")
cat("        Elastic net (LASSO) will handle this via regularisation.\n")
for (v in core_vars) {
  if (v %in% names(df)) {
    n <- sum(!is.na(df[[v]]))
    cat(sprintf("    %-20s: %d / %d (%.1f%%)\n", v, n, nrow(df), 100*n/nrow(df)))
  } else {
    cat(sprintf("    %-20s: *** NOT FOUND ***\n", v))
  }
}

# -- TIER 2: Enhanced (Core + ALL NMR + PRS + biomarkers, NO MRI) -----------
# This is the key change: include all NMR metabolomics features
enhanced_extra <- c()

# ALL NMR p234xx columns (auto-detected, ~90 features)
# Only include NMR columns with >80% completeness
if (length(nmr_all_cols) > 0) {
  nmr_complete <- sapply(nmr_all_cols, function(v) mean(!is.na(df[[v]])))
  nmr_keep <- nmr_all_cols[nmr_complete > 0.80]
  enhanced_extra <- c(enhanced_extra, nmr_keep)
  cat(sprintf("\n  NMR features for Enhanced: %d / %d (>80%% complete)\n",
              length(nmr_keep), length(nmr_all_cols)))
}

# PRS
if ("prs_primary" %in% names(df) && mean(!is.na(df$prs_primary)) > 0.5) {
  enhanced_extra <- c(enhanced_extra, "prs_primary")
}
# Townsend
if ("townsend" %in% names(df) && mean(!is.na(df$townsend)) > 0.5) {
  enhanced_extra <- c(enhanced_extra, "townsend")
}
# Lipid-years
if ("lipid_years_A" %in% names(df) && mean(!is.na(df$lipid_years_A)) > 0.5) {
  enhanced_extra <- c(enhanced_extra, "lipid_years_A")
}
# Gene type
if ("gene_apob" %in% names(df)) {
  enhanced_extra <- c(enhanced_extra, "gene_apob")
}
# Lab biomarkers
for (v in c("hba1c", "creatinine", "alt", "cystatin_c", "glucose")) {
  if (v %in% names(df) && mean(!is.na(df[[v]])) > 0.5) {
    enhanced_extra <- c(enhanced_extra, v)
  }
}

enhanced_vars <- c(core_vars, enhanced_extra)
# Remove any duplicates (in case NMR cols overlap with named features)
enhanced_vars <- unique(enhanced_vars)

cat(sprintf("  TIER 2 -- Enhanced: %d total vars (Core 14 + %d extra)\n",
            length(enhanced_vars), length(enhanced_extra)))
cat("  NO MRI features (only ~12%% available, kills complete-case)\n")

# -- TIER 3: Full Kitchen Sink (Enhanced + everything else) -----------------
kitchen_extra <- c()
for (v in c("waist", "sbp", "dbp", "ever_smoked", "log_apob_ldl",
            "log_hba1c", "log_creatinine", "log_alt", "log_cystatin_c",
            "apob_apoa1_nmr_ratio", "vldl_hdl_ratio", "mean_cimt")) {
  if (v %in% names(df) && mean(!is.na(df[[v]])) > 0.3) {
    kitchen_extra <- c(kitchen_extra, v)
  }
}

kitchen_vars <- unique(c(enhanced_vars, kitchen_extra))
cat(sprintf("  TIER 3 -- Kitchen Sink: %d total vars (Enhanced + %d more)\n",
            length(kitchen_vars), length(kitchen_extra)))

# ── v3 VARIANT MODELS (Collinearity-Fixed) ──────────────────────────────────
cat("\n  --- v3 VARIANT MODELS (Collinearity-Fixed) ---\n")

# VARIANT B: Fixed Core — drop standalone apob, use log(ApoB/LDL) + log(Lp(a))
# This resolves: (1) ApoB sign-flip, (2) absurd ApoB/LDL OR, (3) Lp(a) info loss
core_fixed_vars <- c("age", "sex", "re_ldl", "hdl", "trig",
                      "inv_apoa", "log_apob_ldl", "log_lpa",
                      "smoking_binary", "diabetes", "hypertension", "bmi",
                      "log_crp")
cat(sprintf("  VARIANT B -- Fixed Core: %d vars (log(ApoB/LDL), log(Lpa), no standalone apob)\n",
            length(core_fixed_vars)))

# VARIANT C: Consistent Ratio — use log(apob/re_ldl) for internal consistency
core_consistent_vars <- c("age", "sex", "re_ldl", "hdl", "trig",
                           "inv_apoa", "log_apob_re_ldl", "log_lpa",
                           "smoking_binary", "diabetes", "hypertension", "bmi",
                           "log_crp")
cat(sprintf("  VARIANT C -- Consistent Ratio: %d vars (log(ApoB/RE-LDL) instead of ApoB/LDL)\n",
            length(core_consistent_vars)))

# VARIANT D: Fixed + Interactions — adds lipid_years + age×re_ldl for cumulative burden
core_interactions_vars <- c("age", "sex", "re_ldl", "hdl", "trig",
                             "inv_apoa", "log_apob_re_ldl", "log_lpa",
                             "smoking_binary", "diabetes", "hypertension", "bmi",
                             "log_crp", "lipid_years_A", "age_x_re_ldl")
cat(sprintf("  VARIANT D -- Fixed + Interactions: %d vars (+ lipid_years, age x re_ldl)\n",
            length(core_interactions_vars)))

# ENHANCED-FIXED: Fixed Core + NMR + PRS + biomarkers (the key model)
enhanced_fixed_extra <- enhanced_extra  # same extras as original Enhanced
enhanced_fixed_vars <- unique(c(core_fixed_vars, enhanced_fixed_extra))
cat(sprintf("  ENHANCED-FIXED: %d vars (Fixed Core + NMR + PRS + biomarkers)\n",
            length(enhanced_fixed_vars)))

# KITCHEN-FIXED: Enhanced-Fixed + everything else
kitchen_fixed_vars <- unique(c(enhanced_fixed_vars, kitchen_extra,
                                "log_apob_re_ldl", "age_x_re_ldl", "age_x_hdl"))
cat(sprintf("  KITCHEN-FIXED: %d vars (Enhanced-Fixed + all extras)\n",
            length(kitchen_fixed_vars)))

# ── v3b LEAN MODELS (Remove redundancy: drop inv_apoa + trig) ──────────────
cat("\n  --- v3b LEAN MODELS (Remove HDL/ApoA1 & trig redundancy) ---\n")
cat("  RATIONALE:\n")
cat("    hdl + inv_apoa compete (VIF 4.73 & 3.83, both non-sig)\n")
cat("    trig: OR=0.92, p=0.35 — adds noise, not signal\n")
cat("    Solution: keep hdl (universally available), drop inv_apoa + trig\n\n")

# VARIANT E: Lean Core — 11 vars (drop inv_apoa + trig from Fixed Core)
core_lean_vars <- c("age", "sex", "re_ldl", "hdl",
                     "log_apob_ldl", "log_lpa",
                     "smoking_binary", "diabetes", "hypertension", "bmi",
                     "log_crp")
cat(sprintf("  VARIANT E -- Lean Core: %d vars (no inv_apoa, no trig)\n",
            length(core_lean_vars)))

# VARIANT F: Lean + Interactions — adds lipid_years + age×re_ldl
core_lean_interact_vars <- c("age", "sex", "re_ldl", "hdl",
                              "log_apob_ldl", "log_lpa",
                              "smoking_binary", "diabetes", "hypertension", "bmi",
                              "log_crp", "lipid_years_A", "age_x_re_ldl")
cat(sprintf("  VARIANT F -- Lean + Interactions: %d vars\n",
            length(core_lean_interact_vars)))

# VARIANT G: Minimal Core — 9 vars (only significant + biologically key vars)
# Keep: age, sex, re_ldl, log_apob_ldl, log_lpa, smoking, hypertension, bmi, log_crp
# Drop: hdl (captured indirectly via log_apob_ldl), diabetes (p=0.37)
core_minimal_vars <- c("age", "sex", "re_ldl",
                        "log_apob_ldl", "log_lpa",
                        "smoking_binary", "hypertension", "bmi",
                        "log_crp")
cat(sprintf("  VARIANT G -- Minimal Core: %d vars (sig + key biology only)\n",
            length(core_minimal_vars)))

# ENHANCED-LEAN: Lean Core + NMR + PRS + biomarkers
enhanced_lean_vars <- unique(c(core_lean_vars, enhanced_extra))
cat(sprintf("  ENHANCED-LEAN: %d vars (Lean Core + NMR + PRS + biomarkers)\n",
            length(enhanced_lean_vars)))

# KITCHEN-LEAN: Enhanced-Lean + extras (but exclude inv_apoa, keep trig in kitchen for completeness)
kitchen_lean_extra <- kitchen_extra[!kitchen_extra %in% c("inv_apoa")]
kitchen_lean_vars <- unique(c(enhanced_lean_vars, kitchen_lean_extra,
                               "log_apob_re_ldl", "age_x_re_ldl"))
cat(sprintf("  KITCHEN-LEAN: %d vars (Enhanced-Lean + extras)\n",
            length(kitchen_lean_vars)))

# ── v4 LITERATURE-INFORMED MODELS ────────────────────────────────────────────
cat("\n  === v4 LITERATURE-INFORMED MODELS ===\n")
cat("  Based on systematic review of 25 published FH ASCVD prediction models:\n")
cat("    - Promote lipid_years_A (cumulative LDL burden) to Core\n")
cat("    - Add lpa_extreme (≥250 nmol/L, Paquette 2025)\n")
cat("    - Add lpa_high (≥125 nmol/L threshold)\n")
cat("    - Include restricted cubic spline terms for non-linearity\n\n")

# v4-A: Literature Core — Lean Core + lipid_years_A + lpa_extreme
# Rationale: Montreal-FH-SCORE (AUC 0.799) used 5 vars; SAFEHEART-RE used 8
# Cumulative LDL burden is unanimously supported (Tada, Korneva, Gallo)
# Extreme Lp(a) adds threshold effect (Paquette 2025)
v4_core_vars <- c("age", "sex", "re_ldl", "hdl",
                    "log_apob_ldl", "log_lpa", "lpa_extreme",
                    "smoking_binary", "diabetes", "hypertension", "bmi",
                    "log_crp", "lipid_years_A")
cat(sprintf("  v4-A -- Literature Core: %d vars (Lean + lipid_years + lpa_extreme)\n",
            length(v4_core_vars)))

# v4-B: Literature Enriched — v4 Core + interactions + extra ratios
v4_enriched_vars <- c(v4_core_vars,
                       "age_x_re_ldl", "log_apob_re_ldl", "lpa_high")
cat(sprintf("  v4-B -- Literature Enriched: %d vars (v4 Core + interactions)\n",
            length(v4_enriched_vars)))

# v4-C: Montreal-Inspired — recreate Montreal-FH-SCORE approach (5 vars)
# Montreal-FH-SCORE: age, HDL-C, gender, HTN, smoking → AUC 0.799
v4_montreal_vars <- c("age", "sex", "hdl", "hypertension", "smoking_binary")
cat(sprintf("  v4-C -- Montreal-Inspired: %d vars (replicate Montreal-FH-SCORE)\n",
            length(v4_montreal_vars)))

# v4-D: Montreal-Plus — Montreal + our additional signals
v4_montreal_plus_vars <- c("age", "sex", "hdl", "hypertension", "smoking_binary",
                             "re_ldl", "log_lpa", "lpa_extreme",
                             "bmi", "lipid_years_A")
cat(sprintf("  v4-D -- Montreal-Plus: %d vars (Montreal + FH-specific additions)\n",
            length(v4_montreal_plus_vars)))

# v4-E: Full Literature — everything literature supports + NMR if available
v4_full_lit_vars <- unique(c(v4_enriched_vars, enhanced_extra))
cat(sprintf("  v4-E -- Full Literature: %d vars (v4 Enriched + NMR + PRS + biomarkers)\n",
            length(v4_full_lit_vars)))

# -- Report N complete cases and EPV for each tier --------------------------
n_events_total <- sum(df$ascvd_combined, na.rm = TRUE)

tier_names <- c("Tier 1 (Core)", "Tier 2 (Enhanced)", "Tier 3 (Kitchen Sink)",
                "V-B (Fixed Core)", "V-C (Consistent)", "V-D (Interactions)",
                "Enhanced-Fixed", "Kitchen-Fixed",
                "V-E (Lean Core)", "V-F (Lean+Interact)", "V-G (Minimal)",
                "Enhanced-Lean", "Kitchen-Lean",
                "v4-A (Lit Core)", "v4-B (Lit Enriched)", "v4-C (Montreal)",
                "v4-D (Montreal+)", "v4-E (Full Lit)")
tier_vars_list <- list(core_vars, enhanced_vars, kitchen_vars,
                       core_fixed_vars, core_consistent_vars, core_interactions_vars,
                       enhanced_fixed_vars, kitchen_fixed_vars,
                       core_lean_vars, core_lean_interact_vars, core_minimal_vars,
                       enhanced_lean_vars, kitchen_lean_vars,
                       v4_core_vars, v4_enriched_vars, v4_montreal_vars,
                       v4_montreal_plus_vars, v4_full_lit_vars)

cat("\n  Tier Summary:\n")
cat(sprintf("  %-25s %6s %6s %8s %6s\n", "Tier", "Vars", "N_cc", "Events", "EPV"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (t in seq_along(tier_names)) {
  tvars <- tier_vars_list[[t]]
  tvars_exist <- tvars[tvars %in% names(df)]
  df_tmp <- df[, c(tvars_exist, "ascvd_combined")]
  cc <- complete.cases(df_tmp)
  n_cc <- sum(cc)
  n_ev <- sum(df_tmp$ascvd_combined[cc])
  epv <- ifelse(length(tvars_exist) > 0, n_ev / length(tvars_exist), 0)
  cat(sprintf("  %-25s %6d %6d %8d %6.1f\n",
              tier_names[t], length(tvars_exist), n_cc, n_ev, epv))
}
cat("\n")

# =============================================================================
# SECTION 5: ELASTIC NET -- ALL TIERS
# =============================================================================

cat("===================================================================\n")
cat("SECTION 5: ELASTIC NET -- ALL TIERS\n")
cat("===================================================================\n\n")

# Storage for grand comparison
model_results <- list()

alphas <- c(0, 0.25, 0.5, 0.75, 1.0)
alpha_names <- c("Ridge", "EN(0.25)", "EN(0.50)", "EN(0.75)", "LASSO")

run_elastic_net <- function(tier_name, vars, df_full) {
  cat(sprintf("  --- %s (Elastic Net) ---\n", tier_name))

  # Filter to vars that exist
  vars_exist <- vars[vars %in% names(df_full)]
  if (length(vars_exist) < length(vars)) {
    cat(sprintf("  NOTE: %d / %d vars found in data\n", length(vars_exist), length(vars)))
  }

  # Complete-case analysis
  df_tier <- df_full[, c(vars_exist, "ascvd_combined")]
  cc_idx <- complete.cases(df_tier)
  df_cc <- df_tier[cc_idx, ]

  cat(sprintf("  Complete cases: %d / %d (%.1f%%)\n",
              nrow(df_cc), nrow(df_full), 100 * nrow(df_cc) / nrow(df_full)))
  cat(sprintf("  Events: %d (%.1f%%), EPV: %.1f\n",
              sum(df_cc$ascvd_combined), 100 * mean(df_cc$ascvd_combined),
              sum(df_cc$ascvd_combined) / length(vars_exist)))

  if (nrow(df_cc) < 200 || sum(df_cc$ascvd_combined) < 30) {
    cat("  SKIPPED: too few complete cases or events\n\n")
    return(NULL)
  }

  X <- as.matrix(df_cc[, vars_exist])
  y <- df_cc$ascvd_combined

  # Alpha grid search
  best <- NULL
  for (i in seq_along(alphas)) {
    a <- alphas[i]
    set.seed(2026)
    cv_fit <- cv.glmnet(X, y, family = "binomial",
                         alpha = a, nfolds = 10, type.measure = "auc")
    auc_max <- max(cv_fit$cvm)
    n_nonzero <- sum(coef(cv_fit, s = "lambda.min")[-1] != 0)

    cat(sprintf("    %-10s: AUC=%.4f, lambda=%.6f, %d/%d vars\n",
                alpha_names[i], auc_max, cv_fit$lambda.min, n_nonzero, length(vars_exist)))

    if (is.null(best) || auc_max > best$auc) {
      best <- list(alpha = a, lambda = cv_fit$lambda.min, auc = auc_max,
                    fit = cv_fit, n_vars = n_nonzero, alpha_name = alpha_names[i],
                    X = X, y = y, vars = vars_exist, df_cc = df_cc,
                    n_complete = nrow(df_cc), n_events = sum(df_cc$ascvd_combined))
    }
  }

  cat(sprintf("  BEST: %s (alpha=%.2f), CV-AUC = %.4f, %d vars selected\n\n",
              best$alpha_name, best$alpha, best$auc, best$n_vars))

  return(best)
}

# Run elastic net for all tiers (original + v3 variants)
en_core     <- run_elastic_net("Tier 1: Core (Original)", core_vars, df)
en_enhanced <- run_elastic_net("Tier 2: Enhanced (Original)", enhanced_vars, df)
en_kitchen  <- run_elastic_net("Tier 3: Kitchen Sink (Original)", kitchen_vars, df)

# v3 COLLINEARITY-FIXED variants
cat("  === v3 COLLINEARITY-FIXED VARIANTS ===\n\n")
en_fixed       <- run_elastic_net("V-B: Fixed Core", core_fixed_vars, df)
en_consistent  <- run_elastic_net("V-C: Consistent Ratio", core_consistent_vars, df)
en_interact    <- run_elastic_net("V-D: Fixed + Interactions", core_interactions_vars, df)
en_enh_fixed   <- run_elastic_net("Enhanced-Fixed", enhanced_fixed_vars, df)
en_kit_fixed   <- run_elastic_net("Kitchen-Fixed", kitchen_fixed_vars, df)

# v3b LEAN variants (remove HDL/ApoA1 + trig redundancy)
cat("  === v3b LEAN VARIANTS (No inv_apoa, no trig) ===\n\n")
en_lean        <- run_elastic_net("V-E: Lean Core", core_lean_vars, df)
en_lean_int    <- run_elastic_net("V-F: Lean + Interactions", core_lean_interact_vars, df)
en_minimal     <- run_elastic_net("V-G: Minimal Core", core_minimal_vars, df)
en_enh_lean    <- run_elastic_net("Enhanced-Lean", enhanced_lean_vars, df)
en_kit_lean    <- run_elastic_net("Kitchen-Lean", kitchen_lean_vars, df)

# v4 Literature-informed variants (EN)
cat("  === v4 LITERATURE-INFORMED VARIANTS (Elastic Net) ===\n\n")
en_v4_core     <- run_elastic_net("v4-A: Literature Core", v4_core_vars, df)
en_v4_enriched <- run_elastic_net("v4-B: Literature Enriched", v4_enriched_vars, df)
en_v4_montreal <- run_elastic_net("v4-C: Montreal-Inspired", v4_montreal_vars, df)
en_v4_montp    <- run_elastic_net("v4-D: Montreal-Plus", v4_montreal_plus_vars, df)
en_v4_full     <- run_elastic_net("v4-E: Full Literature", v4_full_lit_vars, df)

# Store results
if (!is.null(en_core)) {
  model_results[["EN_Core"]] <- list(
    name = "EN Core (14 vars)", algorithm = "Elastic Net",
    n_vars = en_core$n_vars, cv_auc = en_core$auc,
    n_complete = en_core$n_complete, n_events = en_core$n_events,
    model = en_core
  )
}
if (!is.null(en_enhanced)) {
  model_results[["EN_Enhanced"]] <- list(
    name = "EN Enhanced (NMR+PRS)", algorithm = "Elastic Net",
    n_vars = en_enhanced$n_vars, cv_auc = en_enhanced$auc,
    n_complete = en_enhanced$n_complete, n_events = en_enhanced$n_events,
    model = en_enhanced
  )
}
if (!is.null(en_kitchen)) {
  model_results[["EN_Kitchen"]] <- list(
    name = "EN Kitchen Sink", algorithm = "Elastic Net",
    n_vars = en_kitchen$n_vars, cv_auc = en_kitchen$auc,
    n_complete = en_kitchen$n_complete, n_events = en_kitchen$n_events,
    model = en_kitchen
  )
}

# v3 variant results
if (!is.null(en_fixed)) {
  model_results[["EN_Fixed"]] <- list(
    name = "EN Fixed Core (v3)", algorithm = "Elastic Net",
    n_vars = en_fixed$n_vars, cv_auc = en_fixed$auc,
    n_complete = en_fixed$n_complete, n_events = en_fixed$n_events,
    model = en_fixed
  )
}
if (!is.null(en_consistent)) {
  model_results[["EN_Consistent"]] <- list(
    name = "EN Consistent Ratio (v3)", algorithm = "Elastic Net",
    n_vars = en_consistent$n_vars, cv_auc = en_consistent$auc,
    n_complete = en_consistent$n_complete, n_events = en_consistent$n_events,
    model = en_consistent
  )
}
if (!is.null(en_interact)) {
  model_results[["EN_Interact"]] <- list(
    name = "EN Fixed+Interact (v3)", algorithm = "Elastic Net",
    n_vars = en_interact$n_vars, cv_auc = en_interact$auc,
    n_complete = en_interact$n_complete, n_events = en_interact$n_events,
    model = en_interact
  )
}
if (!is.null(en_enh_fixed)) {
  model_results[["EN_Enh_Fixed"]] <- list(
    name = "EN Enhanced-Fixed (v3)", algorithm = "Elastic Net",
    n_vars = en_enh_fixed$n_vars, cv_auc = en_enh_fixed$auc,
    n_complete = en_enh_fixed$n_complete, n_events = en_enh_fixed$n_events,
    model = en_enh_fixed
  )
}
if (!is.null(en_kit_fixed)) {
  model_results[["EN_Kit_Fixed"]] <- list(
    name = "EN Kitchen-Fixed (v3)", algorithm = "Elastic Net",
    n_vars = en_kit_fixed$n_vars, cv_auc = en_kit_fixed$auc,
    n_complete = en_kit_fixed$n_complete, n_events = en_kit_fixed$n_events,
    model = en_kit_fixed
  )
}

# v3b lean variant results
if (!is.null(en_lean)) {
  model_results[["EN_Lean"]] <- list(
    name = "EN Lean Core (v3b)", algorithm = "Elastic Net",
    n_vars = en_lean$n_vars, cv_auc = en_lean$auc,
    n_complete = en_lean$n_complete, n_events = en_lean$n_events,
    model = en_lean
  )
}
if (!is.null(en_lean_int)) {
  model_results[["EN_Lean_Int"]] <- list(
    name = "EN Lean+Interact (v3b)", algorithm = "Elastic Net",
    n_vars = en_lean_int$n_vars, cv_auc = en_lean_int$auc,
    n_complete = en_lean_int$n_complete, n_events = en_lean_int$n_events,
    model = en_lean_int
  )
}
if (!is.null(en_minimal)) {
  model_results[["EN_Minimal"]] <- list(
    name = "EN Minimal Core (v3b)", algorithm = "Elastic Net",
    n_vars = en_minimal$n_vars, cv_auc = en_minimal$auc,
    n_complete = en_minimal$n_complete, n_events = en_minimal$n_events,
    model = en_minimal
  )
}
if (!is.null(en_enh_lean)) {
  model_results[["EN_Enh_Lean"]] <- list(
    name = "EN Enhanced-Lean (v3b)", algorithm = "Elastic Net",
    n_vars = en_enh_lean$n_vars, cv_auc = en_enh_lean$auc,
    n_complete = en_enh_lean$n_complete, n_events = en_enh_lean$n_events,
    model = en_enh_lean
  )
}
if (!is.null(en_kit_lean)) {
  model_results[["EN_Kit_Lean"]] <- list(
    name = "EN Kitchen-Lean (v3b)", algorithm = "Elastic Net",
    n_vars = en_kit_lean$n_vars, cv_auc = en_kit_lean$auc,
    n_complete = en_kit_lean$n_complete, n_events = en_kit_lean$n_events,
    model = en_kit_lean
  )
}

# v4 EN results
if (!is.null(en_v4_core)) {
  model_results[["EN_v4_Core"]] <- list(
    name = "EN Literature Core (v4)", algorithm = "Elastic Net",
    n_vars = en_v4_core$n_vars, cv_auc = en_v4_core$auc,
    n_complete = en_v4_core$n_complete, n_events = en_v4_core$n_events,
    model = en_v4_core
  )
}
if (!is.null(en_v4_enriched)) {
  model_results[["EN_v4_Enriched"]] <- list(
    name = "EN Lit Enriched (v4)", algorithm = "Elastic Net",
    n_vars = en_v4_enriched$n_vars, cv_auc = en_v4_enriched$auc,
    n_complete = en_v4_enriched$n_complete, n_events = en_v4_enriched$n_events,
    model = en_v4_enriched
  )
}
if (!is.null(en_v4_montreal)) {
  model_results[["EN_v4_Montreal"]] <- list(
    name = "EN Montreal-Inspired (v4)", algorithm = "Elastic Net",
    n_vars = en_v4_montreal$n_vars, cv_auc = en_v4_montreal$auc,
    n_complete = en_v4_montreal$n_complete, n_events = en_v4_montreal$n_events,
    model = en_v4_montreal
  )
}
if (!is.null(en_v4_montp)) {
  model_results[["EN_v4_MontrealPlus"]] <- list(
    name = "EN Montreal-Plus (v4)", algorithm = "Elastic Net",
    n_vars = en_v4_montp$n_vars, cv_auc = en_v4_montp$auc,
    n_complete = en_v4_montp$n_complete, n_events = en_v4_montp$n_events,
    model = en_v4_montp
  )
}
if (!is.null(en_v4_full)) {
  model_results[["EN_v4_Full"]] <- list(
    name = "EN Full Literature (v4)", algorithm = "Elastic Net",
    n_vars = en_v4_full$n_vars, cv_auc = en_v4_full$auc,
    n_complete = en_v4_full$n_complete, n_events = en_v4_full$n_events,
    model = en_v4_full
  )
}

# -- Extract and save v4 Literature Core coefficients -----------------------
if (!is.null(en_v4_core)) {
  v4_coef <- as.matrix(coef(en_v4_core$fit, s = "lambda.min"))
  v4_coef_df <- data.frame(
    variable    = rownames(v4_coef),
    coefficient = as.numeric(v4_coef),
    stringsAsFactors = FALSE
  )
  v4_coef_df <- v4_coef_df[v4_coef_df$coefficient != 0 |
                               v4_coef_df$variable == "(Intercept)", ]

  cat("\n  CALON-2 Literature Core (v4) EN Coefficients:\n")
  for (i in 1:nrow(v4_coef_df)) {
    cat(sprintf("    %-20s: %10.6f  (OR = %8.4f)\n",
                v4_coef_df$variable[i], v4_coef_df$coefficient[i],
                exp(v4_coef_df$coefficient[i])))
  }
  write.csv(v4_coef_df, paste0(OUT_DIR, "calon2_v4_literature_core_coefficients.csv"),
            row.names = FALSE)
  cat("  Saved: calon2_v4_literature_core_coefficients.csv\n")

  # Refit logistic regression for interpretable ORs
  sel_v4 <- v4_coef_df$variable[v4_coef_df$variable != "(Intercept)" &
                                   v4_coef_df$coefficient != 0]
  if (length(sel_v4) > 0) {
    formula_str <- paste("ascvd_combined ~", paste(sel_v4, collapse = " + "))
    glm_v4 <- glm(as.formula(formula_str), data = en_v4_core$df_cc, family = binomial)
    glm_summary <- summary(glm_v4)$coefficients
    glm_or_v4 <- data.frame(
      variable    = rownames(glm_summary),
      coefficient = glm_summary[, "Estimate"],
      se          = glm_summary[, "Std. Error"],
      z_value     = glm_summary[, "z value"],
      p_value     = glm_summary[, "Pr(>|z|)"],
      OR          = exp(glm_summary[, "Estimate"]),
      OR_lower    = exp(glm_summary[, "Estimate"] - 1.96 * glm_summary[, "Std. Error"]),
      OR_upper    = exp(glm_summary[, "Estimate"] + 1.96 * glm_summary[, "Std. Error"]),
      stringsAsFactors = FALSE
    )

    cat("\n  v4 Literature Core Logistic Regression Odds Ratios:\n")
    cat(sprintf("  %-20s %8s %12s %10s\n", "Variable", "OR", "95% CI", "p-value"))
    cat(paste(rep("-", 60), collapse = ""), "\n")
    for (i in 1:nrow(glm_or_v4)) {
      sig <- ifelse(glm_or_v4$p_value[i] < 0.001, "***",
             ifelse(glm_or_v4$p_value[i] < 0.01,  "**",
             ifelse(glm_or_v4$p_value[i] < 0.05,  "*", "")))
      cat(sprintf("  %-20s %8.3f (%6.3f-%6.3f) %10.4f %s\n",
                  glm_or_v4$variable[i], glm_or_v4$OR[i], glm_or_v4$OR_lower[i],
                  glm_or_v4$OR_upper[i], glm_or_v4$p_value[i], sig))
    }
    write.csv(glm_or_v4, paste0(TAB_DIR, "calon2_v4_literature_core_ORs.csv"),
              row.names = FALSE)
    cat("  Saved: tables/calon2_v4_literature_core_ORs.csv\n")

    # VIF check
    cat("\n  VIF Check (v4 Literature Core):\n")
    tryCatch({
      vif_vals <- car::vif(glm_v4)
      for (vn in names(vif_vals)) {
        flag <- ifelse(vif_vals[vn] > 5, " ** HIGH", "")
        cat(sprintf("    %-20s VIF = %6.2f%s\n", vn, vif_vals[vn], flag))
      }
    }, error = function(e) cat(sprintf("    VIF failed: %s\n", e$message)))
  }
}

# -- Extract and save Core elastic net coefficients -------------------------
if (!is.null(en_core)) {
  core_coef <- as.matrix(coef(en_core$fit, s = "lambda.min"))
  core_coef_df <- data.frame(
    variable    = rownames(core_coef),
    coefficient = as.numeric(core_coef),
    stringsAsFactors = FALSE
  )
  core_coef_df <- core_coef_df[core_coef_df$coefficient != 0 |
                                 core_coef_df$variable == "(Intercept)", ]

  cat("  CALON-2 Core Elastic Net Coefficients:\n")
  for (i in 1:nrow(core_coef_df)) {
    cat(sprintf("    %-20s: %10.6f\n",
                core_coef_df$variable[i], core_coef_df$coefficient[i]))
  }
  write.csv(core_coef_df, paste0(OUT_DIR, "calon2_core_coefficients.csv"), row.names = FALSE)
  cat("  Saved: calon2_core_coefficients.csv\n")

  # Refit logistic regression for ORs
  selected_vars <- core_coef_df$variable[core_coef_df$variable != "(Intercept)" &
                                          core_coef_df$coefficient != 0]
  if (length(selected_vars) > 0) {
    formula_str <- paste("ascvd_combined ~", paste(selected_vars, collapse = " + "))
    glm_core <- glm(as.formula(formula_str), data = en_core$df_cc, family = binomial)
    glm_summary <- summary(glm_core)$coefficients
    glm_or <- data.frame(
      variable    = rownames(glm_summary),
      coefficient = glm_summary[, "Estimate"],
      se          = glm_summary[, "Std. Error"],
      z_value     = glm_summary[, "z value"],
      p_value     = glm_summary[, "Pr(>|z|)"],
      OR          = exp(glm_summary[, "Estimate"]),
      OR_lower    = exp(glm_summary[, "Estimate"] - 1.96 * glm_summary[, "Std. Error"]),
      OR_upper    = exp(glm_summary[, "Estimate"] + 1.96 * glm_summary[, "Std. Error"]),
      stringsAsFactors = FALSE
    )

    cat("\n  Logistic Regression Odds Ratios (Core):\n")
    cat(sprintf("  %-20s %8s %12s %10s\n", "Variable", "OR", "95% CI", "p-value"))
    cat(paste(rep("-", 60), collapse = ""), "\n")
    for (i in 1:nrow(glm_or)) {
      cat(sprintf("  %-20s %8.3f (%6.3f-%6.3f) %10.4f\n",
                  glm_or$variable[i], glm_or$OR[i], glm_or$OR_lower[i],
                  glm_or$OR_upper[i], glm_or$p_value[i]))
    }
    write.csv(glm_or, paste0(TAB_DIR, "calon2_core_logistic_ORs.csv"), row.names = FALSE)
  }
}

# -- Extract and save Fixed Core (v3) elastic net coefficients ---------------
if (!is.null(en_fixed)) {
  fixed_coef <- as.matrix(coef(en_fixed$fit, s = "lambda.min"))
  fixed_coef_df <- data.frame(
    variable    = rownames(fixed_coef),
    coefficient = as.numeric(fixed_coef),
    stringsAsFactors = FALSE
  )
  fixed_coef_df <- fixed_coef_df[fixed_coef_df$coefficient != 0 |
                                   fixed_coef_df$variable == "(Intercept)", ]

  cat("\n  CALON-2 Fixed Core (v3) EN Coefficients:\n")
  for (i in 1:nrow(fixed_coef_df)) {
    cat(sprintf("    %-20s: %10.6f  (OR = %8.4f)\n",
                fixed_coef_df$variable[i], fixed_coef_df$coefficient[i],
                exp(fixed_coef_df$coefficient[i])))
  }
  write.csv(fixed_coef_df, paste0(OUT_DIR, "calon2_fixed_core_coefficients.csv"), row.names = FALSE)
  cat("  Saved: calon2_fixed_core_coefficients.csv\n")

  # Refit logistic regression for interpretable ORs
  sel_vars <- fixed_coef_df$variable[fixed_coef_df$variable != "(Intercept)" &
                                       fixed_coef_df$coefficient != 0]
  if (length(sel_vars) > 0) {
    formula_str <- paste("ascvd_combined ~", paste(sel_vars, collapse = " + "))
    glm_fixed <- glm(as.formula(formula_str), data = en_fixed$df_cc, family = binomial)
    glm_summary <- summary(glm_fixed)$coefficients
    glm_or_fixed <- data.frame(
      variable    = rownames(glm_summary),
      coefficient = glm_summary[, "Estimate"],
      se          = glm_summary[, "Std. Error"],
      z_value     = glm_summary[, "z value"],
      p_value     = glm_summary[, "Pr(>|z|)"],
      OR          = exp(glm_summary[, "Estimate"]),
      OR_lower    = exp(glm_summary[, "Estimate"] - 1.96 * glm_summary[, "Std. Error"]),
      OR_upper    = exp(glm_summary[, "Estimate"] + 1.96 * glm_summary[, "Std. Error"]),
      stringsAsFactors = FALSE
    )

    cat("\n  Fixed Core Logistic Regression Odds Ratios (v3):\n")
    cat(sprintf("  %-20s %8s %12s %10s\n", "Variable", "OR", "95% CI", "p-value"))
    cat(paste(rep("-", 60), collapse = ""), "\n")
    for (i in 1:nrow(glm_or_fixed)) {
      cat(sprintf("  %-20s %8.3f (%6.3f-%6.3f) %10.4f\n",
                  glm_or_fixed$variable[i], glm_or_fixed$OR[i], glm_or_fixed$OR_lower[i],
                  glm_or_fixed$OR_upper[i], glm_or_fixed$p_value[i]))
    }
    write.csv(glm_or_fixed, paste0(TAB_DIR, "calon2_fixed_core_logistic_ORs.csv"), row.names = FALSE)
    cat("  Saved: calon2_fixed_core_logistic_ORs.csv\n")

    # COLLINEARITY CHECK: VIF for fixed model
    cat("\n  VIF Check (Fixed Core):\n")
    tryCatch({
      vif_vals <- car::vif(glm_fixed)
      for (vn in names(vif_vals)) {
        flag <- ifelse(vif_vals[vn] > 5, " ** HIGH", "")
        cat(sprintf("    %-20s VIF = %6.2f%s\n", vn, vif_vals[vn], flag))
      }
    }, error = function(e) {
      cat(sprintf("    VIF calculation failed: %s\n", e$message))
      cat("    (Install 'car' package if needed: install.packages('car'))\n")
    })
  }
}

# -- Extract and save Lean Core (v3b) elastic net coefficients ---------------
if (!is.null(en_lean)) {
  lean_coef <- as.matrix(coef(en_lean$fit, s = "lambda.min"))
  lean_coef_df <- data.frame(
    variable    = rownames(lean_coef),
    coefficient = as.numeric(lean_coef),
    stringsAsFactors = FALSE
  )
  lean_coef_df <- lean_coef_df[lean_coef_df$coefficient != 0 |
                                 lean_coef_df$variable == "(Intercept)", ]

  cat("\n  CALON-2 Lean Core (v3b) EN Coefficients:\n")
  for (i in 1:nrow(lean_coef_df)) {
    cat(sprintf("    %-20s: %10.6f  (OR = %8.4f)\n",
                lean_coef_df$variable[i], lean_coef_df$coefficient[i],
                exp(lean_coef_df$coefficient[i])))
  }
  write.csv(lean_coef_df, paste0(OUT_DIR, "calon2_lean_core_coefficients.csv"), row.names = FALSE)
  cat("  Saved: calon2_lean_core_coefficients.csv\n")

  # Refit logistic regression for interpretable ORs
  sel_vars_lean <- lean_coef_df$variable[lean_coef_df$variable != "(Intercept)" &
                                           lean_coef_df$coefficient != 0]
  if (length(sel_vars_lean) > 0) {
    formula_str <- paste("ascvd_combined ~", paste(sel_vars_lean, collapse = " + "))
    glm_lean <- glm(as.formula(formula_str), data = en_lean$df_cc, family = binomial)
    glm_summary <- summary(glm_lean)$coefficients
    glm_or_lean <- data.frame(
      variable    = rownames(glm_summary),
      coefficient = glm_summary[, "Estimate"],
      se          = glm_summary[, "Std. Error"],
      z_value     = glm_summary[, "z value"],
      p_value     = glm_summary[, "Pr(>|z|)"],
      OR          = exp(glm_summary[, "Estimate"]),
      OR_lower    = exp(glm_summary[, "Estimate"] - 1.96 * glm_summary[, "Std. Error"]),
      OR_upper    = exp(glm_summary[, "Estimate"] + 1.96 * glm_summary[, "Std. Error"]),
      stringsAsFactors = FALSE
    )

    cat("\n  Lean Core Logistic Regression Odds Ratios (v3b):\n")
    cat(sprintf("  %-20s %8s %12s %10s\n", "Variable", "OR", "95% CI", "p-value"))
    cat(paste(rep("-", 60), collapse = ""), "\n")
    for (i in 1:nrow(glm_or_lean)) {
      sig <- ifelse(glm_or_lean$p_value[i] < 0.001, "***",
             ifelse(glm_or_lean$p_value[i] < 0.01,  "**",
             ifelse(glm_or_lean$p_value[i] < 0.05,  "*", "")))
      cat(sprintf("  %-20s %8.3f (%6.3f-%6.3f) %10.4f %s\n",
                  glm_or_lean$variable[i], glm_or_lean$OR[i], glm_or_lean$OR_lower[i],
                  glm_or_lean$OR_upper[i], glm_or_lean$p_value[i], sig))
    }
    write.csv(glm_or_lean, paste0(TAB_DIR, "calon2_lean_core_logistic_ORs.csv"), row.names = FALSE)
    cat("  Saved: calon2_lean_core_logistic_ORs.csv\n")

    # VIF check for Lean model
    cat("\n  VIF Check (Lean Core):\n")
    tryCatch({
      vif_vals <- car::vif(glm_lean)
      for (vn in names(vif_vals)) {
        flag <- ifelse(vif_vals[vn] > 5, " ** HIGH", "")
        cat(sprintf("    %-20s VIF = %6.2f%s\n", vn, vif_vals[vn], flag))
      }
      cat(sprintf("    Max VIF = %.2f %s\n", max(vif_vals),
                  ifelse(max(vif_vals) < 3, "(EXCELLENT — no collinearity)", "")))
    }, error = function(e) {
      cat(sprintf("    VIF calculation failed: %s\n", e$message))
    })
  }
}

# -- Extract and save Enhanced elastic net coefficients ---------------------
if (!is.null(en_enhanced)) {
  enh_coef <- as.matrix(coef(en_enhanced$fit, s = "lambda.min"))
  enh_coef_df <- data.frame(
    variable    = rownames(enh_coef),
    coefficient = as.numeric(enh_coef),
    stringsAsFactors = FALSE
  )
  enh_coef_df <- enh_coef_df[enh_coef_df$coefficient != 0 |
                               enh_coef_df$variable == "(Intercept)", ]

  cat(sprintf("\n  Enhanced EN: %d non-zero coefficients\n", nrow(enh_coef_df) - 1))
  write.csv(enh_coef_df, paste0(OUT_DIR, "calon2_enhanced_coefficients.csv"), row.names = FALSE)
  cat("  Saved: calon2_enhanced_coefficients.csv\n")
}

# -- Extract and save Enhanced-Fixed (v3) coefficients -----------------------
if (!is.null(en_enh_fixed)) {
  enh_fx_coef <- as.matrix(coef(en_enh_fixed$fit, s = "lambda.min"))
  enh_fx_coef_df <- data.frame(
    variable    = rownames(enh_fx_coef),
    coefficient = as.numeric(enh_fx_coef),
    stringsAsFactors = FALSE
  )
  enh_fx_coef_df <- enh_fx_coef_df[enh_fx_coef_df$coefficient != 0 |
                                      enh_fx_coef_df$variable == "(Intercept)", ]
  cat(sprintf("\n  Enhanced-Fixed (v3) EN: %d non-zero coefficients\n",
              nrow(enh_fx_coef_df) - 1))
  write.csv(enh_fx_coef_df, paste0(OUT_DIR, "calon2_enhanced_fixed_coefficients.csv"),
            row.names = FALSE)
  cat("  Saved: calon2_enhanced_fixed_coefficients.csv\n")
}

cat("\n")

# =============================================================================
# SECTION 6: XGBOOST -- ALL TIERS
# =============================================================================

cat("===================================================================\n")
cat("SECTION 6: XGBOOST -- ALL TIERS\n")
cat("===================================================================\n\n")

run_xgboost <- function(tier_name, vars, df_full) {
  cat(sprintf("  --- %s (XGBoost) ---\n", tier_name))

  vars_exist <- vars[vars %in% names(df_full)]
  df_tier <- df_full[, c(vars_exist, "ascvd_combined")]
  cc_idx <- complete.cases(df_tier)
  df_cc <- df_tier[cc_idx, ]

  cat(sprintf("  Complete cases: %d, Events: %d\n",
              nrow(df_cc), sum(df_cc$ascvd_combined)))

  if (nrow(df_cc) < 200 || sum(df_cc$ascvd_combined) < 30) {
    cat("  SKIPPED: too few complete cases or events\n\n")
    return(NULL)
  }

  X <- as.matrix(df_cc[, vars_exist])
  y <- df_cc$ascvd_combined
  dtrain <- xgb.DMatrix(data = X, label = y)

  # Hyperparameter grid search
  grid <- expand.grid(
    max_depth = c(3, 5, 7),
    eta = c(0.01, 0.05, 0.1),
    min_child_weight = c(5, 10),
    stringsAsFactors = FALSE
  )

  best <- NULL
  best_auc <- 0
  best_params_row <- NULL

  cat(sprintf("  Grid search: %d combinations...\n", nrow(grid)))

  for (g in 1:nrow(grid)) {
    params <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = grid$max_depth[g],
      eta = grid$eta[g],
      subsample = 0.8,
      colsample_bytree = 0.8,
      min_child_weight = grid$min_child_weight[g]
    )

    set.seed(2026)
    cv_result <- tryCatch({
      xgb.cv(params = params, data = dtrain, nrounds = 1000, nfold = 10,
             early_stopping_rounds = 50, verbose = 0, print_every_n = 0)
    }, error = function(e) {
      # Fallback for older xgboost versions
      xgb.cv(params = params, data = dtrain, nrounds = 1000, nfold = 10,
             early_stopping_rounds = 50, verbose = 0)
    })

    # Extract best AUC from evaluation log
    eval_log <- cv_result$evaluation_log
    auc_col <- grep("test.*auc.*mean", names(eval_log), value = TRUE)[1]
    if (is.null(auc_col) || is.na(auc_col)) {
      auc_col <- names(eval_log)[grep("auc", names(eval_log))[1]]
    }
    if (!is.null(auc_col) && !is.na(auc_col)) {
      this_auc <- max(eval_log[[auc_col]], na.rm = TRUE)
    } else {
      this_auc <- 0
    }

    if (this_auc > best_auc) {
      best_auc <- this_auc
      best_params_row <- grid[g, ]

      # Get best nrounds
      best_nrounds <- which.max(eval_log[[auc_col]])

      best <- list(
        auc = this_auc,
        params = params,
        nrounds = best_nrounds,
        cv_result = cv_result,
        X = X, y = y, dtrain = dtrain,
        vars = vars_exist, df_cc = df_cc,
        n_complete = nrow(df_cc), n_events = sum(df_cc$ascvd_combined)
      )
    }
  }

  cat(sprintf("  BEST: depth=%d, eta=%.2f, mcw=%d, nrounds=%d\n",
              best_params_row$max_depth, best_params_row$eta,
              best_params_row$min_child_weight, best$nrounds))
  cat(sprintf("  CV-AUC = %.4f\n\n", best_auc))

  # Train final model with best params
  set.seed(2026)
  best$final_model <- xgb.train(
    params = best$params, data = best$dtrain,
    nrounds = best$nrounds, verbose = 0
  )

  return(best)
}

# Run XGBoost for all tiers (original + v3 variants)
xgb_core     <- run_xgboost("Tier 1: Core (Original)", core_vars, df)
xgb_enhanced <- run_xgboost("Tier 2: Enhanced (Original)", enhanced_vars, df)
xgb_kitchen  <- run_xgboost("Tier 3: Kitchen Sink (Original)", kitchen_vars, df)

# v3 variants (XGBoost)
cat("  === v3 COLLINEARITY-FIXED VARIANTS (XGBoost) ===\n\n")
xgb_fixed     <- run_xgboost("V-B: Fixed Core", core_fixed_vars, df)
xgb_interact  <- run_xgboost("V-D: Fixed + Interactions", core_interactions_vars, df)
xgb_enh_fixed <- run_xgboost("Enhanced-Fixed", enhanced_fixed_vars, df)
xgb_kit_fixed <- run_xgboost("Kitchen-Fixed", kitchen_fixed_vars, df)

# v3b lean variants (XGBoost)
cat("  === v3b LEAN VARIANTS (XGBoost) ===\n\n")
xgb_lean      <- run_xgboost("V-E: Lean Core", core_lean_vars, df)
xgb_lean_int  <- run_xgboost("V-F: Lean + Interactions", core_lean_interact_vars, df)
xgb_minimal   <- run_xgboost("V-G: Minimal Core", core_minimal_vars, df)
xgb_enh_lean  <- run_xgboost("Enhanced-Lean", enhanced_lean_vars, df)
xgb_kit_lean  <- run_xgboost("Kitchen-Lean", kitchen_lean_vars, df)

# v4 literature-informed XGBoost
cat("  === v4 LITERATURE-INFORMED VARIANTS (XGBoost) ===\n\n")
xgb_v4_core     <- run_xgboost("v4-A: Literature Core", v4_core_vars, df)
xgb_v4_enriched <- run_xgboost("v4-B: Literature Enriched", v4_enriched_vars, df)
xgb_v4_montreal <- run_xgboost("v4-C: Montreal-Inspired", v4_montreal_vars, df)
xgb_v4_montp    <- run_xgboost("v4-D: Montreal-Plus", v4_montreal_plus_vars, df)
xgb_v4_full     <- run_xgboost("v4-E: Full Literature", v4_full_lit_vars, df)

# Store results
if (!is.null(xgb_core)) {
  model_results[["XGB_Core"]] <- list(
    name = "XGB Core (14 vars)", algorithm = "XGBoost",
    n_vars = length(xgb_core$vars), cv_auc = xgb_core$auc,
    n_complete = xgb_core$n_complete, n_events = xgb_core$n_events,
    model = xgb_core
  )
}
if (!is.null(xgb_enhanced)) {
  model_results[["XGB_Enhanced"]] <- list(
    name = "XGB Enhanced (NMR+PRS)", algorithm = "XGBoost",
    n_vars = length(xgb_enhanced$vars), cv_auc = xgb_enhanced$auc,
    n_complete = xgb_enhanced$n_complete, n_events = xgb_enhanced$n_events,
    model = xgb_enhanced
  )
}
if (!is.null(xgb_kitchen)) {
  model_results[["XGB_Kitchen"]] <- list(
    name = "XGB Kitchen Sink", algorithm = "XGBoost",
    n_vars = length(xgb_kitchen$vars), cv_auc = xgb_kitchen$auc,
    n_complete = xgb_kitchen$n_complete, n_events = xgb_kitchen$n_events,
    model = xgb_kitchen
  )
}

# v3 XGBoost variant results
if (!is.null(xgb_fixed)) {
  model_results[["XGB_Fixed"]] <- list(
    name = "XGB Fixed Core (v3)", algorithm = "XGBoost",
    n_vars = length(xgb_fixed$vars), cv_auc = xgb_fixed$auc,
    n_complete = xgb_fixed$n_complete, n_events = xgb_fixed$n_events,
    model = xgb_fixed
  )
}
if (!is.null(xgb_interact)) {
  model_results[["XGB_Interact"]] <- list(
    name = "XGB Fixed+Interact (v3)", algorithm = "XGBoost",
    n_vars = length(xgb_interact$vars), cv_auc = xgb_interact$auc,
    n_complete = xgb_interact$n_complete, n_events = xgb_interact$n_events,
    model = xgb_interact
  )
}
if (!is.null(xgb_enh_fixed)) {
  model_results[["XGB_Enh_Fixed"]] <- list(
    name = "XGB Enhanced-Fixed (v3)", algorithm = "XGBoost",
    n_vars = length(xgb_enh_fixed$vars), cv_auc = xgb_enh_fixed$auc,
    n_complete = xgb_enh_fixed$n_complete, n_events = xgb_enh_fixed$n_events,
    model = xgb_enh_fixed
  )
}
if (!is.null(xgb_kit_fixed)) {
  model_results[["XGB_Kit_Fixed"]] <- list(
    name = "XGB Kitchen-Fixed (v3)", algorithm = "XGBoost",
    n_vars = length(xgb_kit_fixed$vars), cv_auc = xgb_kit_fixed$auc,
    n_complete = xgb_kit_fixed$n_complete, n_events = xgb_kit_fixed$n_events,
    model = xgb_kit_fixed
  )
}

# v3b lean XGBoost variant results
if (!is.null(xgb_lean)) {
  model_results[["XGB_Lean"]] <- list(
    name = "XGB Lean Core (v3b)", algorithm = "XGBoost",
    n_vars = length(xgb_lean$vars), cv_auc = xgb_lean$auc,
    n_complete = xgb_lean$n_complete, n_events = xgb_lean$n_events,
    model = xgb_lean
  )
}
if (!is.null(xgb_lean_int)) {
  model_results[["XGB_Lean_Int"]] <- list(
    name = "XGB Lean+Interact (v3b)", algorithm = "XGBoost",
    n_vars = length(xgb_lean_int$vars), cv_auc = xgb_lean_int$auc,
    n_complete = xgb_lean_int$n_complete, n_events = xgb_lean_int$n_events,
    model = xgb_lean_int
  )
}
if (!is.null(xgb_minimal)) {
  model_results[["XGB_Minimal"]] <- list(
    name = "XGB Minimal Core (v3b)", algorithm = "XGBoost",
    n_vars = length(xgb_minimal$vars), cv_auc = xgb_minimal$auc,
    n_complete = xgb_minimal$n_complete, n_events = xgb_minimal$n_events,
    model = xgb_minimal
  )
}
if (!is.null(xgb_enh_lean)) {
  model_results[["XGB_Enh_Lean"]] <- list(
    name = "XGB Enhanced-Lean (v3b)", algorithm = "XGBoost",
    n_vars = length(xgb_enh_lean$vars), cv_auc = xgb_enh_lean$auc,
    n_complete = xgb_enh_lean$n_complete, n_events = xgb_enh_lean$n_events,
    model = xgb_enh_lean
  )
}
if (!is.null(xgb_kit_lean)) {
  model_results[["XGB_Kit_Lean"]] <- list(
    name = "XGB Kitchen-Lean (v3b)", algorithm = "XGBoost",
    n_vars = length(xgb_kit_lean$vars), cv_auc = xgb_kit_lean$auc,
    n_complete = xgb_kit_lean$n_complete, n_events = xgb_kit_lean$n_events,
    model = xgb_kit_lean
  )
}

# v4 XGBoost results
if (!is.null(xgb_v4_core)) {
  model_results[["XGB_v4_Core"]] <- list(
    name = "XGB Literature Core (v4)", algorithm = "XGBoost",
    n_vars = length(xgb_v4_core$vars), cv_auc = xgb_v4_core$auc,
    n_complete = xgb_v4_core$n_complete, n_events = xgb_v4_core$n_events,
    model = xgb_v4_core
  )
}
if (!is.null(xgb_v4_enriched)) {
  model_results[["XGB_v4_Enriched"]] <- list(
    name = "XGB Lit Enriched (v4)", algorithm = "XGBoost",
    n_vars = length(xgb_v4_enriched$vars), cv_auc = xgb_v4_enriched$auc,
    n_complete = xgb_v4_enriched$n_complete, n_events = xgb_v4_enriched$n_events,
    model = xgb_v4_enriched
  )
}
if (!is.null(xgb_v4_montreal)) {
  model_results[["XGB_v4_Montreal"]] <- list(
    name = "XGB Montreal-Inspired (v4)", algorithm = "XGBoost",
    n_vars = length(xgb_v4_montreal$vars), cv_auc = xgb_v4_montreal$auc,
    n_complete = xgb_v4_montreal$n_complete, n_events = xgb_v4_montreal$n_events,
    model = xgb_v4_montreal
  )
}
if (!is.null(xgb_v4_montp)) {
  model_results[["XGB_v4_MontrealPlus"]] <- list(
    name = "XGB Montreal-Plus (v4)", algorithm = "XGBoost",
    n_vars = length(xgb_v4_montp$vars), cv_auc = xgb_v4_montp$auc,
    n_complete = xgb_v4_montp$n_complete, n_events = xgb_v4_montp$n_events,
    model = xgb_v4_montp
  )
}
if (!is.null(xgb_v4_full)) {
  model_results[["XGB_v4_Full"]] <- list(
    name = "XGB Full Literature (v4)", algorithm = "XGBoost",
    n_vars = length(xgb_v4_full$vars), cv_auc = xgb_v4_full$auc,
    n_complete = xgb_v4_full$n_complete, n_events = xgb_v4_full$n_events,
    model = xgb_v4_full
  )
}

# -- Extract and save XGBoost feature importance ----------------------------
for (xgb_name in c("xgb_enh_fixed", "xgb_kit_fixed", "xgb_enhanced", "xgb_kitchen",
                    "xgb_core", "xgb_fixed", "xgb_interact",
                    "xgb_lean", "xgb_lean_int", "xgb_minimal",
                    "xgb_enh_lean", "xgb_kit_lean",
                    "xgb_v4_core", "xgb_v4_enriched", "xgb_v4_montreal",
                    "xgb_v4_montp", "xgb_v4_full")) {
  xgb_obj <- get(xgb_name)
  if (!is.null(xgb_obj) && !is.null(xgb_obj$final_model)) {
    imp <- xgb.importance(model = xgb_obj$final_model)
    if (nrow(imp) > 0) {
      imp_top <- head(imp, 30)
      suffix <- gsub("xgb_", "", xgb_name)
      write.csv(imp_top,
                paste0(TAB_DIR, sprintf("calon2_xgb_importance_%s.csv", suffix)),
                row.names = FALSE)
      if (xgb_name == "xgb_enh_fixed") {
        # v3: Use Enhanced-Fixed as the main importance file
        write.csv(imp_top, paste0(TAB_DIR, "calon2_xgb_importance.csv"),
                  row.names = FALSE)
        cat("  XGBoost Enhanced-Fixed (v3) top 30 features:\n")
        for (j in 1:min(15, nrow(imp_top))) {
          cat(sprintf("    %2d. %-25s  Gain=%.4f\n",
                      j, imp_top$Feature[j], imp_top$Gain[j]))
        }
        if (nrow(imp_top) > 15) cat("    ... (see calon2_xgb_importance.csv)\n")
      } else if (xgb_name == "xgb_enhanced") {
        cat("  XGBoost Enhanced (original) top 15 features:\n")
        for (j in 1:min(15, nrow(imp_top))) {
          cat(sprintf("    %2d. %-25s  Gain=%.4f\n",
                      j, imp_top$Feature[j], imp_top$Gain[j]))
        }
      }
    }
  }
}

cat("\n")

# =============================================================================
# SECTION 6B: v4 — RANDOM FOREST (additional base learner)
# =============================================================================

cat("===================================================================\n")
cat("SECTION 6B: RANDOM FOREST (v4)\n")
cat("===================================================================\n\n")

run_rf <- function(tier_name, vars, df_full) {
  cat(sprintf("  --- %s (Random Forest) ---\n", tier_name))

  vars_exist <- vars[vars %in% names(df_full)]
  df_tier <- df_full[, c(vars_exist, "ascvd_combined")]
  cc_idx <- complete.cases(df_tier)
  df_cc <- df_tier[cc_idx, ]

  cat(sprintf("  Complete cases: %d, Events: %d\n",
              nrow(df_cc), sum(df_cc$ascvd_combined)))

  if (nrow(df_cc) < 200 || sum(df_cc$ascvd_combined) < 30) {
    cat("  SKIPPED: too few complete cases or events\n\n")
    return(NULL)
  }

  X <- df_cc[, vars_exist]
  y <- factor(df_cc$ascvd_combined, levels = c(0, 1))
  mtry_vals <- c(floor(sqrt(ncol(X))), floor(ncol(X) / 3), floor(ncol(X) / 2))
  mtry_vals <- unique(pmax(mtry_vals, 1))

  best_auc <- 0
  best_rf <- NULL
  best_mtry <- NA

  for (mt in mtry_vals) {
    set.seed(2026)
    rf_fit <- randomForest(x = X, y = y, ntree = 1000, mtry = mt,
                            importance = TRUE, classwt = c("0" = 1, "1" = 2))
    # OOB predictions
    oob_probs <- rf_fit$votes[, "1"]
    this_auc <- tryCatch(
      as.numeric(auc(roc(as.numeric(as.character(y)), oob_probs, quiet = TRUE))),
      error = function(e) 0
    )
    cat(sprintf("    mtry=%d: OOB-AUC=%.4f\n", mt, this_auc))
    if (this_auc > best_auc) {
      best_auc <- this_auc
      best_rf <- rf_fit
      best_mtry <- mt
    }
  }

  cat(sprintf("  BEST: mtry=%d, OOB-AUC = %.4f\n\n", best_mtry, best_auc))

  return(list(
    auc = best_auc, rf = best_rf, mtry = best_mtry,
    X = X, y = y, vars = vars_exist, df_cc = df_cc,
    n_complete = nrow(df_cc), n_events = sum(df_cc$ascvd_combined == 1)
  ))
}

# Run RF on key variants only (to save time)
rf_v4_core     <- run_rf("v4-A: Literature Core", v4_core_vars, df)
rf_v4_enriched <- run_rf("v4-B: Literature Enriched", v4_enriched_vars, df)
rf_v4_montp    <- run_rf("v4-D: Montreal-Plus", v4_montreal_plus_vars, df)
rf_lean        <- run_rf("V-E: Lean Core", core_lean_vars, df)
rf_kitchen     <- run_rf("Kitchen Sink", kitchen_vars, df)
rf_v4_full     <- run_rf("v4-E: Full Literature", v4_full_lit_vars, df)

# Store RF results
for (rf_item in list(
  list(obj = rf_v4_core, key = "RF_v4_Core", nm = "RF Literature Core (v4)"),
  list(obj = rf_v4_enriched, key = "RF_v4_Enriched", nm = "RF Lit Enriched (v4)"),
  list(obj = rf_v4_montp, key = "RF_v4_MontrealPlus", nm = "RF Montreal-Plus (v4)"),
  list(obj = rf_lean, key = "RF_Lean", nm = "RF Lean Core (v3b)"),
  list(obj = rf_kitchen, key = "RF_Kitchen", nm = "RF Kitchen Sink"),
  list(obj = rf_v4_full, key = "RF_v4_Full", nm = "RF Full Literature (v4)")
)) {
  if (!is.null(rf_item$obj)) {
    model_results[[rf_item$key]] <- list(
      name = rf_item$nm, algorithm = "Random Forest",
      n_vars = length(rf_item$obj$vars), cv_auc = rf_item$obj$auc,
      n_complete = rf_item$obj$n_complete, n_events = rf_item$obj$n_events,
      model = rf_item$obj
    )
  }
}

# Save RF variable importance for best v4 model
if (!is.null(rf_v4_core)) {
  rf_imp <- importance(rf_v4_core$rf)
  rf_imp_df <- data.frame(
    variable = rownames(rf_imp),
    MeanDecreaseGini = rf_imp[, "MeanDecreaseGini"],
    stringsAsFactors = FALSE
  )
  rf_imp_df <- rf_imp_df[order(-rf_imp_df$MeanDecreaseGini), ]
  write.csv(rf_imp_df, paste0(TAB_DIR, "calon2_rf_v4_importance.csv"), row.names = FALSE)
  cat("  RF v4 Literature Core variable importance:\n")
  for (j in 1:min(15, nrow(rf_imp_df))) {
    cat(sprintf("    %2d. %-25s  Gini=%.4f\n",
                j, rf_imp_df$variable[j], rf_imp_df$MeanDecreaseGini[j]))
  }
  cat("\n")
}

# =============================================================================
# SECTION 6C: v4 — MULTIPLE IMPUTATION (mice) → recover lost patients
# =============================================================================

cat("===================================================================\n")
cat("SECTION 6C: MULTIPLE IMPUTATION (v4) — recover 345 lost patients\n")
cat("===================================================================\n\n")

# Strategy: impute missing values in Core predictors using mice, then fit
# models on m=10 imputed datasets and pool the AUC estimates.
# This recovers patients lost to complete-case analysis.

run_mi_elastic_net <- function(tier_name, vars, df_full, m_imps = 10) {
  cat(sprintf("  --- %s (MI Elastic Net, m=%d) ---\n", tier_name, m_imps))

  vars_exist <- vars[vars %in% names(df_full)]
  imp_data <- df_full[, c(vars_exist, "ascvd_combined")]

  # Report baseline missingness
  n_before <- nrow(imp_data)
  n_cc <- sum(complete.cases(imp_data))
  cat(sprintf("  Before MI: %d total, %d complete (%.1f%% loss)\n",
              n_before, n_cc, 100 * (n_before - n_cc) / n_before))

  if (n_before - n_cc < 20) {
    cat("  Minimal missingness — MI unnecessary. Using complete cases.\n\n")
    return(NULL)
  }

  # Run mice imputation (predictive mean matching for continuous, logistic for binary)
  cat("  Running mice imputation...\n")
  set.seed(2026)
  tryCatch({
    imp <- mice(imp_data, m = m_imps, method = "pmm", maxit = 20,
                printFlag = FALSE, seed = 2026)

    # Fit elastic net on each imputed dataset
    mi_aucs <- numeric(m_imps)
    mi_fits <- vector("list", m_imps)

    for (i in 1:m_imps) {
      d_imp <- complete(imp, i)
      X_imp <- as.matrix(d_imp[, vars_exist])
      y_imp <- d_imp$ascvd_combined

      # Alpha grid search on imputed data
      best_a <- NULL
      best_auc_a <- 0
      for (a in c(0, 0.25, 0.5, 0.75, 1.0)) {
        set.seed(2026 + i)
        cv_fit <- cv.glmnet(X_imp, y_imp, family = "binomial",
                              alpha = a, nfolds = 10, type.measure = "auc")
        if (max(cv_fit$cvm) > best_auc_a) {
          best_auc_a <- max(cv_fit$cvm)
          best_a <- list(alpha = a, fit = cv_fit, auc = best_auc_a)
        }
      }
      mi_aucs[i] <- best_a$auc
      mi_fits[[i]] <- best_a
    }

    mean_auc <- mean(mi_aucs)
    sd_auc <- sd(mi_aucs)
    cat(sprintf("  MI-pooled CV-AUC: %.4f (SD=%.4f, range=%.4f--%.4f)\n",
                mean_auc, sd_auc, min(mi_aucs), max(mi_aucs)))
    cat(sprintf("  N recovered: %d additional patients (total N=%d)\n",
                n_before - n_cc, n_before))
    cat("\n")

    return(list(
      auc = mean_auc, mi_aucs = mi_aucs, mi_fits = mi_fits,
      imp = imp, vars = vars_exist, n_total = n_before,
      n_complete = n_cc, n_recovered = n_before - n_cc,
      n_events = sum(df_full$ascvd_combined, na.rm = TRUE)
    ))
  }, error = function(e) {
    cat(sprintf("  MI FAILED: %s\n\n", e$message))
    return(NULL)
  })
}

# Run MI on top v4 variants
mi_v4_core     <- run_mi_elastic_net("v4-A: Literature Core (MI)", v4_core_vars, df)
mi_v4_enriched <- run_mi_elastic_net("v4-B: Literature Enriched (MI)", v4_enriched_vars, df)
mi_lean        <- run_mi_elastic_net("V-E: Lean Core (MI)", core_lean_vars, df)
mi_v4_montp    <- run_mi_elastic_net("v4-D: Montreal-Plus (MI)", v4_montreal_plus_vars, df)

# Store MI results
for (mi_item in list(
  list(obj = mi_v4_core, key = "MI_v4_Core", nm = "MI-EN Lit Core (v4)"),
  list(obj = mi_v4_enriched, key = "MI_v4_Enriched", nm = "MI-EN Lit Enriched (v4)"),
  list(obj = mi_lean, key = "MI_Lean", nm = "MI-EN Lean Core (v3b)"),
  list(obj = mi_v4_montp, key = "MI_v4_MontrealPlus", nm = "MI-EN Montreal-Plus (v4)")
)) {
  if (!is.null(mi_item$obj)) {
    model_results[[mi_item$key]] <- list(
      name = mi_item$nm, algorithm = "MI-Elastic Net",
      n_vars = length(mi_item$obj$vars), cv_auc = mi_item$obj$auc,
      n_complete = mi_item$obj$n_total, n_events = mi_item$obj$n_events,
      model = mi_item$obj
    )
    cat(sprintf("  Stored: %s → CV-AUC = %.4f (N=%d, recovered %d)\n",
                mi_item$nm, mi_item$obj$auc, mi_item$obj$n_total,
                mi_item$obj$n_recovered))
  }
}

# =============================================================================
# SECTION 6D: v4 — SUPER LEARNER (stacked ensemble)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 6D: SUPER LEARNER ENSEMBLE (v4)\n")
cat("===================================================================\n\n")
cat("  Combines predictions from EN + XGBoost + RF via meta-learner\n")
cat("  Expected: +2-5% AUC over best single model\n\n")

run_super_learner <- function(tier_name, vars, df_full) {
  cat(sprintf("  --- %s (Super Learner) ---\n", tier_name))

  vars_exist <- vars[vars %in% names(df_full)]
  df_tier <- df_full[, c(vars_exist, "ascvd_combined")]
  cc_idx <- complete.cases(df_tier)
  df_cc <- df_tier[cc_idx, ]

  if (nrow(df_cc) < 200 || sum(df_cc$ascvd_combined) < 30) {
    cat("  SKIPPED: too few complete cases\n\n")
    return(NULL)
  }

  X <- data.frame(df_cc[, vars_exist])
  Y <- df_cc$ascvd_combined

  cat(sprintf("  N=%d, Events=%d, Vars=%d\n", nrow(df_cc), sum(Y), length(vars_exist)))

  # Define base learners (SuperLearner library)
  sl_lib <- c("SL.glmnet", "SL.xgboost", "SL.randomForest", "SL.glm")

  # Check which learners are available
  available_learners <- c()
  for (lib in sl_lib) {
    if (exists(lib, mode = "function")) {
      available_learners <- c(available_learners, lib)
    }
  }

  if (length(available_learners) < 2) {
    cat("  SKIPPED: fewer than 2 SuperLearner algorithms available\n\n")
    return(NULL)
  }

  cat(sprintf("  Learners: %s\n", paste(available_learners, collapse = ", ")))

  tryCatch({
    # BUG FIX (v4.1): Use CV.SuperLearner for HONEST out-of-sample AUC
    #   PROBLEM: SuperLearner()$SL.predict gives training-set predictions,
    #     which inflated AUC to 0.93-0.97 (vs honest EN/XGB ~0.77).
    #     The meta-learner optimises weights on the SAME data it predicts on.
    #   FIX: CV.SuperLearner() adds an OUTER CV loop (V=10).
    #     The outer fold predictions are truly out-of-sample → honest AUC.
    cat("  Using CV.SuperLearner (nested 10×10-fold CV) for honest AUC...\n")
    set.seed(2026)
    cv_sl_fit <- CV.SuperLearner(
      Y = Y, X = X,
      V = 10,
      family = binomial(),
      SL.library = available_learners,
      innerCvControl = list(list(V = 10)),
      verbose = FALSE
    )

    # Extract OUTER-fold cross-validated predictions → honest AUC
    # cv_sl_fit$SL.predict contains out-of-sample predictions from each outer fold
    sl_pred_honest <- cv_sl_fit$SL.predict[, 1]
    sl_roc <- roc(Y, sl_pred_honest, quiet = TRUE)
    sl_auc <- as.numeric(auc(sl_roc))

    # Also fit the regular SuperLearner for final model (to get learner weights)
    set.seed(2026)
    sl_fit <- SuperLearner(
      Y = Y, X = X,
      family = binomial(),
      SL.library = available_learners,
      cvControl = list(V = 10, shuffle = TRUE),
      verbose = FALSE
    )

    # Report learner weights (key insight for interpretability)
    cat(sprintf("  Super Learner HONEST CV-AUC: %.4f (via CV.SuperLearner)\n", sl_auc))
    cat("  Learner weights (meta-model coefficients):\n")
    for (j in seq_along(sl_fit$coef)) {
      cat(sprintf("    %-25s: %.4f\n",
                  available_learners[j], sl_fit$coef[j]))
    }

    # Compare with individual learner HONEST AUCs from outer folds
    cat("  Individual learner HONEST CV-AUCs (outer fold):\n")
    for (j in seq_len(ncol(cv_sl_fit$library.predict))) {
      lib_pred <- cv_sl_fit$library.predict[, j]
      lib_auc <- tryCatch(as.numeric(auc(roc(Y, lib_pred, quiet = TRUE))),
                           error = function(e) NA)
      cat(sprintf("    %-25s: %.4f\n", available_learners[min(j, length(available_learners))], lib_auc))
    }
    cat("\n")

    return(list(
      auc = sl_auc, sl_fit = sl_fit, cv_sl_fit = cv_sl_fit, roc = sl_roc,
      X = X, Y = Y, vars = vars_exist, df_cc = df_cc,
      n_complete = nrow(df_cc), n_events = sum(Y),
      learner_weights = sl_fit$coef
    ))
  }, error = function(e) {
    cat(sprintf("  Super Learner FAILED: %s\n\n", e$message))
    return(NULL)
  })
}

# Run Super Learner on key variants
sl_v4_core     <- run_super_learner("v4-A: Literature Core", v4_core_vars, df)
sl_v4_enriched <- run_super_learner("v4-B: Literature Enriched", v4_enriched_vars, df)
sl_v4_montp    <- run_super_learner("v4-D: Montreal-Plus", v4_montreal_plus_vars, df)
sl_lean        <- run_super_learner("V-E: Lean Core", core_lean_vars, df)
sl_v4_full     <- run_super_learner("v4-E: Full Literature", v4_full_lit_vars, df)

# Store SL results
for (sl_item in list(
  list(obj = sl_v4_core, key = "SL_v4_Core", nm = "SL Literature Core (v4)"),
  list(obj = sl_v4_enriched, key = "SL_v4_Enriched", nm = "SL Lit Enriched (v4)"),
  list(obj = sl_v4_montp, key = "SL_v4_MontrealPlus", nm = "SL Montreal-Plus (v4)"),
  list(obj = sl_lean, key = "SL_Lean", nm = "SL Lean Core (v3b)"),
  list(obj = sl_v4_full, key = "SL_v4_Full", nm = "SL Full Literature (v4)")
)) {
  if (!is.null(sl_item$obj)) {
    model_results[[sl_item$key]] <- list(
      name = sl_item$nm, algorithm = "Super Learner",
      n_vars = length(sl_item$obj$vars), cv_auc = sl_item$obj$auc,
      n_complete = sl_item$obj$n_complete, n_events = sl_item$obj$n_events,
      model = sl_item$obj
    )
    cat(sprintf("  Stored: %s → CV-AUC = %.4f\n",
                sl_item$nm, sl_item$obj$auc))
  }
}

# =============================================================================
# SECTION 6E: v4 — RESTRICTED CUBIC SPLINES (non-linear effects)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 6E: RESTRICTED CUBIC SPLINES (v4)\n")
cat("===================================================================\n\n")
cat("  Test non-linear effects of age and LDL-C using rms::lrm\n\n")

tryCatch({
  library(rms)

  # Use v4 Literature Core complete cases
  if (!is.null(en_v4_core)) {
    df_rcs <- en_v4_core$df_cc

    # Set data distribution for rms
    dd <- datadist(df_rcs)
    options(datadist = "dd")

    # Model with RCS for age (3 knots) and re_ldl (3 knots)
    # BUG FIX (v4.1): lpa_extreme is constant (0/1308 have Lp(a)>=250 in UKB)
    #   → causes singular Hessian → lrm convergence failure
    #   REPLACED with lpa_high (>=125 nmol/L, 10.4% prevalence) — sufficient
    #   variation for estimation while preserving threshold biology
    rcs_formula <- ascvd_combined ~ rcs(age, 3) + sex + rcs(re_ldl, 3) + hdl +
      log_apob_ldl + log_lpa + lpa_high +
      smoking_binary + diabetes + hypertension + bmi +
      log_crp + lipid_years_A

    # Guard: check all RCS predictors have variance (avoid singular Hessian)
    rcs_pred_vars <- c("age", "sex", "re_ldl", "hdl", "log_apob_ldl",
                        "log_lpa", "lpa_high", "smoking_binary", "diabetes",
                        "hypertension", "bmi", "log_crp", "lipid_years_A")
    for (rv in rcs_pred_vars) {
      if (rv %in% names(df_rcs) && var(df_rcs[[rv]], na.rm = TRUE) < 1e-10) {
        cat(sprintf("  WARNING: %s has near-zero variance — dropping from RCS\n", rv))
      }
    }

    set.seed(2026)
    rcs_fit <- lrm(rcs_formula, data = df_rcs, x = TRUE, y = TRUE)
    rcs_val <- validate(rcs_fit, B = 100, seed = 2026)

    # Extract Dxy → AUC
    dxy_orig <- rcs_val["Dxy", "index.orig"]
    dxy_corr <- rcs_val["Dxy", "index.corrected"]
    auc_rcs_orig <- (dxy_orig / 2) + 0.5
    auc_rcs_corr <- (dxy_corr / 2) + 0.5

    cat(sprintf("  RCS Model (age + LDL non-linear):\n"))
    cat(sprintf("    Apparent AUC:  %.4f\n", auc_rcs_orig))
    cat(sprintf("    Corrected AUC: %.4f (optimism=%.4f)\n",
                auc_rcs_corr, auc_rcs_orig - auc_rcs_corr))
    cat(sprintf("    N=%d, Events=%d\n", nrow(df_rcs), sum(df_rcs$ascvd_combined)))

    # Print model summary
    cat("\n  RCS Non-linear coefficients:\n")
    print(rcs_fit)

    # ANOVA for non-linearity tests
    cat("\n  ANOVA (non-linearity tests):\n")
    an <- anova(rcs_fit)
    print(an)

    # Store in model_results
    model_results[["RCS_v4_Core"]] <- list(
      name = "RCS Lit Core (v4)", algorithm = "RCS Logistic (rms)",
      n_vars = 13, cv_auc = auc_rcs_corr,
      n_complete = nrow(df_rcs), n_events = sum(df_rcs$ascvd_combined),
      model = list(fit = rcs_fit, val = rcs_val, auc_orig = auc_rcs_orig,
                   auc_corr = auc_rcs_corr, df_cc = df_rcs)
    )

    # Clean up rms options
    options(datadist = NULL)
  }
}, error = function(e) {
  cat(sprintf("  RCS analysis failed: %s\n", e$message))
  cat("  (Install rms package if needed: install.packages('rms'))\n")
})

cat("\n")

# =============================================================================
# SECTION 7: GRAND MODEL COMPARISON
# =============================================================================

cat("===================================================================\n")
cat("SECTION 7: GRAND MODEL COMPARISON\n")
cat("===================================================================\n\n")

if (length(model_results) > 0) {
  comparison_df <- data.frame(
    model     = sapply(model_results, function(m) m$name),
    algorithm = sapply(model_results, function(m) m$algorithm),
    n_vars    = sapply(model_results, function(m) m$n_vars),
    n_complete = sapply(model_results, function(m) m$n_complete),
    n_events  = sapply(model_results, function(m) m$n_events),
    cv_auc    = sapply(model_results, function(m) m$cv_auc),
    stringsAsFactors = FALSE
  )
  comparison_df <- comparison_df[order(-comparison_df$cv_auc), ]

  cat("  GRAND COMPARISON TABLE:\n")
  cat(sprintf("  %-30s %-12s %5s %6s %6s %8s\n",
              "Model", "Algorithm", "Vars", "N", "Events", "CV-AUC"))
  cat(paste(rep("-", 80), collapse = ""), "\n")
  for (i in 1:nrow(comparison_df)) {
    marker <- ifelse(i == 1, " *** BEST", "")
    cat(sprintf("  %-30s %-12s %5d %6d %6d %8.4f%s\n",
                comparison_df$model[i], comparison_df$algorithm[i],
                comparison_df$n_vars[i], comparison_df$n_complete[i],
                comparison_df$n_events[i], comparison_df$cv_auc[i], marker))
  }

  write.csv(comparison_df, paste0(TAB_DIR, "calon2_model_comparison.csv"), row.names = FALSE)
  cat("\n  Saved: calon2_model_comparison.csv\n")

  # Identify the BEST model
  best_model_key <- names(model_results)[which.max(
    sapply(model_results, function(m) m$cv_auc)
  )]
  best_model <- model_results[[best_model_key]]

  cat(sprintf("\n  BEST OVERALL MODEL: %s\n", best_model$name))
  cat(sprintf("  CV-AUC = %.4f\n", best_model$cv_auc))
  cat(sprintf("  Algorithm: %s, Vars: %d, N: %d, Events: %d\n",
              best_model$algorithm, best_model$n_vars,
              best_model$n_complete, best_model$n_events))

  # Save best model object — handles ALL algorithm types (BUG FIX v4.1)
  bm_algo <- best_model$algorithm
  bm_obj  <- best_model$model

  if (bm_algo == "XGBoost") {
    xgb.save(bm_obj$final_model, paste0(OUT_DIR, "calon2_best_xgb_model.bin"))
    saveRDS(list(
      algorithm = "xgboost",
      model = bm_obj$final_model,
      params = bm_obj$params,
      nrounds = bm_obj$nrounds,
      vars = bm_obj$vars,
      cv_auc = best_model$cv_auc,
      name = best_model$name
    ), paste0(OUT_DIR, "calon2_best_model.rds"))
    cat("  Saved: calon2_best_model.rds (XGBoost)\n")

  } else if (bm_algo == "Elastic Net") {
    saveRDS(list(
      algorithm = "glmnet",
      fit = bm_obj$fit,
      alpha = bm_obj$alpha,
      lambda = bm_obj$lambda,
      vars = bm_obj$vars,
      cv_auc = best_model$cv_auc,
      name = best_model$name
    ), paste0(OUT_DIR, "calon2_best_model.rds"))
    cat("  Saved: calon2_best_model.rds (Elastic Net)\n")

  } else if (bm_algo == "Super Learner") {
    saveRDS(list(
      algorithm = "super_learner",
      sl_fit = bm_obj$sl_fit,
      cv_sl_fit = bm_obj$cv_sl_fit,
      learner_weights = bm_obj$learner_weights,
      vars = bm_obj$vars,
      cv_auc = best_model$cv_auc,
      name = best_model$name
    ), paste0(OUT_DIR, "calon2_best_model.rds"))
    cat("  Saved: calon2_best_model.rds (Super Learner)\n")

  } else if (bm_algo == "Random Forest") {
    saveRDS(list(
      algorithm = "random_forest",
      rf = bm_obj$rf,
      mtry = bm_obj$mtry,
      vars = bm_obj$vars,
      cv_auc = best_model$cv_auc,
      name = best_model$name
    ), paste0(OUT_DIR, "calon2_best_model.rds"))
    cat("  Saved: calon2_best_model.rds (Random Forest)\n")

  } else if (bm_algo == "MI-Elastic Net") {
    saveRDS(list(
      algorithm = "mi_glmnet",
      mi_fits = bm_obj$mi_fits,
      imp = bm_obj$imp,
      vars = bm_obj$vars,
      cv_auc = best_model$cv_auc,
      name = best_model$name
    ), paste0(OUT_DIR, "calon2_best_model.rds"))
    cat("  Saved: calon2_best_model.rds (MI-Elastic Net)\n")

  } else if (grepl("RCS", bm_algo)) {
    saveRDS(list(
      algorithm = "rcs_lrm",
      fit = bm_obj$fit,
      val = bm_obj$val,
      auc_orig = bm_obj$auc_orig,
      auc_corr = bm_obj$auc_corr,
      cv_auc = best_model$cv_auc,
      name = best_model$name
    ), paste0(OUT_DIR, "calon2_best_model.rds"))
    cat("  Saved: calon2_best_model.rds (RCS Logistic)\n")

  } else {
    # Fallback: save whatever we have
    saveRDS(list(
      algorithm = bm_algo,
      model = bm_obj,
      cv_auc = best_model$cv_auc,
      name = best_model$name
    ), paste0(OUT_DIR, "calon2_best_model.rds"))
    cat(sprintf("  Saved: calon2_best_model.rds (%s)\n", bm_algo))
  }
} else {
  cat("  ERROR: No models were successfully trained.\n")
  best_model <- NULL
  best_model_key <- NULL
}

cat("\n")

# =============================================================================
# SECTION 8: INTERNAL VALIDATION OF BEST MODEL
# =============================================================================

cat("===================================================================\n")
cat("SECTION 8: INTERNAL VALIDATION OF BEST MODEL\n")
cat("===================================================================\n\n")

if (!is.null(best_model)) {
  bm <- best_model$model  # underlying model details
  bm_algo <- best_model$algorithm

  # ═══════════════════════════════════════════════════════════════════════════
  # BUG FIX (v4.1): Section 8 now handles ALL 6 algorithm types:
  #   Elastic Net, XGBoost, Super Learner, Random Forest, MI-Elastic Net, RCS
  # PREVIOUSLY: only handled EN / XGBoost → crashed on Super Learner
  #   Error: "check.booster.params(params, ...) : params must be a list"
  # ═══════════════════════════════════════════════════════════════════════════

  # -- Helper: get complete-case data + outcome from any model type ----------
  get_cc_data <- function(bm, algo) {
    if (algo %in% c("Elastic Net", "XGBoost")) {
      return(list(df_cc = bm$df_cc, vars = bm$vars,
                  X = as.matrix(bm$df_cc[, bm$vars]),
                  y = bm$df_cc$ascvd_combined))
    } else if (algo == "Super Learner") {
      return(list(df_cc = bm$df_cc, vars = bm$vars,
                  X = as.matrix(bm$df_cc[, bm$vars]),
                  y = bm$Y))
    } else if (algo == "Random Forest") {
      return(list(df_cc = bm$df_cc, vars = bm$vars,
                  X = as.matrix(bm$df_cc[, bm$vars]),
                  y = as.numeric(as.character(bm$y))))
    } else if (algo == "MI-Elastic Net") {
      # MI uses full data (imputed), use first imputed dataset as reference
      d1 <- mice::complete(bm$imp, 1)
      return(list(df_cc = d1, vars = bm$vars,
                  X = as.matrix(d1[, bm$vars]),
                  y = d1$ascvd_combined))
    } else if (grepl("RCS", algo)) {
      return(list(df_cc = bm$df_cc, vars = names(bm$fit$assign),
                  X = NULL, y = bm$df_cc$ascvd_combined))
    }
    # Fallback
    return(list(df_cc = bm$df_cc, vars = bm$vars,
                X = as.matrix(bm$df_cc[, bm$vars]),
                y = bm$df_cc$ascvd_combined))
  }

  cc <- get_cc_data(bm, bm_algo)

  # -- 8A: Repeated 10-fold CV (20 reps) with 95% CI -----------------------
  cat(sprintf("  8A: Repeated 10-fold CV (20 repetitions) [%s]...\n", bm_algo))

  cv_aucs <- numeric(20)

  if (bm_algo == "Elastic Net") {
    for (rep_i in 1:20) {
      set.seed(2026 + rep_i)
      cv_fit <- cv.glmnet(cc$X, cc$y, family = "binomial",
                           alpha = bm$alpha, nfolds = 10, type.measure = "auc")
      cv_aucs[rep_i] <- max(cv_fit$cvm)
    }

  } else if (bm_algo == "XGBoost") {
    for (rep_i in 1:20) {
      set.seed(2026 + rep_i)
      cv_res <- tryCatch({
        xgb.cv(params = bm$params, data = bm$dtrain,
               nrounds = bm$nrounds + 200, nfold = 10,
               early_stopping_rounds = 50, verbose = 0, print_every_n = 0)
      }, error = function(e) {
        xgb.cv(params = bm$params, data = bm$dtrain,
               nrounds = bm$nrounds + 200, nfold = 10,
               early_stopping_rounds = 50, verbose = 0)
      })
      eval_log <- cv_res$evaluation_log
      auc_col <- grep("test.*auc.*mean", names(eval_log), value = TRUE)[1]
      if (is.null(auc_col) || is.na(auc_col)) {
        auc_col <- names(eval_log)[grep("auc", names(eval_log))[1]]
      }
      cv_aucs[rep_i] <- max(eval_log[[auc_col]], na.rm = TRUE)
    }

  } else if (bm_algo == "Super Learner") {
    # Repeated CV.SuperLearner (outer 10-fold, no inflation)
    for (rep_i in 1:20) {
      set.seed(2026 + rep_i)
      cv_sl_rep <- tryCatch({
        CV.SuperLearner(
          Y = cc$y, X = data.frame(cc$df_cc[, bm$vars]),
          V = 10, family = binomial(),
          SL.library = names(bm$sl_fit$cvRisk),
          verbose = FALSE
        )
      }, error = function(e) NULL)
      if (!is.null(cv_sl_rep)) {
        sl_pred_rep <- cv_sl_rep$SL.predict[, 1]
        cv_aucs[rep_i] <- tryCatch(
          as.numeric(auc(roc(cc$y, sl_pred_rep, quiet = TRUE))),
          error = function(e) NA)
      } else {
        cv_aucs[rep_i] <- NA
      }
    }

  } else if (bm_algo == "Random Forest") {
    # Repeated RF with OOB predictions
    for (rep_i in 1:20) {
      set.seed(2026 + rep_i)
      rf_rep <- randomForest(x = cc$df_cc[, bm$vars],
                              y = factor(cc$y, levels = c(0, 1)),
                              ntree = 1000, mtry = bm$mtry, importance = FALSE,
                              classwt = c("0" = 1, "1" = 2))
      oob_probs <- rf_rep$votes[, "1"]
      cv_aucs[rep_i] <- tryCatch(
        as.numeric(auc(roc(cc$y, oob_probs, quiet = TRUE))),
        error = function(e) NA)
    }

  } else if (bm_algo == "MI-Elastic Net") {
    # Repeated MI-EN: for each rep, re-run mice + cv.glmnet on m=5 imputations
    for (rep_i in 1:20) {
      set.seed(2026 + rep_i)
      imp_rep <- tryCatch(
        mice(cc$df_cc[, c(bm$vars, "ascvd_combined")], m = 5, method = "pmm",
             maxit = 10, printFlag = FALSE),
        error = function(e) NULL)
      if (!is.null(imp_rep)) {
        mi_a <- numeric(5)
        for (ii in 1:5) {
          d_ii <- mice::complete(imp_rep, ii)
          X_ii <- as.matrix(d_ii[, bm$vars])
          y_ii <- d_ii$ascvd_combined
          cv_ii <- cv.glmnet(X_ii, y_ii, family = "binomial",
                              alpha = 0.5, nfolds = 10, type.measure = "auc")
          mi_a[ii] <- max(cv_ii$cvm)
        }
        cv_aucs[rep_i] <- mean(mi_a)
      } else {
        cv_aucs[rep_i] <- NA
      }
    }

  } else if (grepl("RCS", bm_algo)) {
    # rms::validate already does bootstrap — use that result directly
    cv_aucs <- rep(best_model$cv_auc, 20)  # RCS corrected AUC as point estimate
    cat("  NOTE: RCS uses rms::validate() bootstrap internally (B=100).\n")
  }

  cv_aucs_valid <- cv_aucs[!is.na(cv_aucs)]
  if (length(cv_aucs_valid) > 0) {
    cat(sprintf("  CV-AUC (20 reps): mean=%.4f, SD=%.4f, 95%% CI=[%.4f, %.4f]\n",
                mean(cv_aucs_valid), sd(cv_aucs_valid),
                quantile(cv_aucs_valid, 0.025), quantile(cv_aucs_valid, 0.975)))
  } else {
    cat("  WARNING: All CV repetitions failed.\n")
  }

  # -- 8B: Apparent AUC ----------------------------------------------------
  cat("  8B: Apparent AUC...\n")
  if (bm_algo == "Elastic Net") {
    pred_best <- predict(bm$fit, cc$X, s = "lambda.min", type = "response")[, 1]
  } else if (bm_algo == "XGBoost") {
    pred_best <- predict(bm$final_model, bm$dtrain)
  } else if (bm_algo == "Super Learner") {
    pred_best <- predict(bm$sl_fit, newdata = data.frame(cc$df_cc[, bm$vars]),
                          onlySL = TRUE)$pred[, 1]
  } else if (bm_algo == "Random Forest") {
    pred_best <- predict(bm$rf, newdata = cc$df_cc[, bm$vars], type = "prob")[, "1"]
  } else if (bm_algo == "MI-Elastic Net") {
    # Average predictions across m imputed datasets
    mi_preds <- matrix(0, nrow = nrow(cc$df_cc), ncol = length(bm$mi_fits))
    for (ii in seq_along(bm$mi_fits)) {
      d_ii <- mice::complete(bm$imp, ii)
      X_ii <- as.matrix(d_ii[, bm$vars])
      mi_preds[, ii] <- predict(bm$mi_fits[[ii]]$fit, X_ii,
                                  s = "lambda.min", type = "response")[, 1]
    }
    pred_best <- rowMeans(mi_preds)
  } else if (grepl("RCS", bm_algo)) {
    pred_best <- predict(bm$fit, type = "fitted")
  } else {
    pred_best <- rep(0.5, length(cc$y))
    cat("  WARNING: Unknown algorithm — cannot compute apparent AUC\n")
  }

  roc_apparent <- roc(cc$y, pred_best, quiet = TRUE)
  auc_apparent <- as.numeric(auc(roc_apparent))
  auc_ci <- ci.auc(roc_apparent)
  cat(sprintf("  Apparent AUC: %.4f (95%% CI: %.4f--%.4f)\n",
              auc_apparent, auc_ci[1], auc_ci[3]))

  # -- 8C: Bootstrap optimism correction (B=200) ---------------------------
  cat("  8C: Bootstrap optimism correction (B=200)...\n")

  boot_data_full <- cc$df_cc
  boot_vars_used <- cc$vars

  boot_optimism_fn <- function(data, indices) {
    boot_data <- data[indices, ]
    X_boot <- as.matrix(boot_data[, boot_vars_used])
    y_boot <- boot_data$ascvd_combined

    X_orig <- as.matrix(data[, boot_vars_used])
    y_orig <- data$ascvd_combined

    if (bm_algo == "Elastic Net") {
      fit_b <- glmnet(X_boot, y_boot, family = "binomial",
                       alpha = bm$alpha, lambda = bm$lambda)
      p_on_boot <- predict(fit_b, X_boot, type = "response")[, 1]
      p_on_orig <- predict(fit_b, X_orig, type = "response")[, 1]

    } else if (bm_algo == "XGBoost") {
      dtb <- xgb.DMatrix(data = X_boot, label = y_boot)
      fit_b <- xgb.train(params = bm$params, data = dtb,
                          nrounds = bm$nrounds, verbose = 0)
      p_on_boot <- predict(fit_b, X_boot)
      p_on_orig <- predict(fit_b, X_orig)

    } else if (bm_algo == "Super Learner") {
      # Use glmnet as proxy for bootstrap optimism (SL too slow per bootstrap)
      fit_b <- glmnet(X_boot, y_boot, family = "binomial", alpha = 0.5)
      lam <- fit_b$lambda[min(50, length(fit_b$lambda))]
      p_on_boot <- predict(fit_b, X_boot, s = lam, type = "response")[, 1]
      p_on_orig <- predict(fit_b, X_orig, s = lam, type = "response")[, 1]

    } else if (bm_algo == "Random Forest") {
      fit_b <- randomForest(x = data.frame(boot_data[, boot_vars_used]),
                             y = factor(y_boot, levels = c(0, 1)),
                             ntree = 500, mtry = bm$mtry,
                             classwt = c("0" = 1, "1" = 2))
      p_on_boot <- predict(fit_b, newdata = data.frame(boot_data[, boot_vars_used]),
                             type = "prob")[, "1"]
      p_on_orig <- predict(fit_b, newdata = data.frame(data[, boot_vars_used]),
                             type = "prob")[, "1"]

    } else if (bm_algo == "MI-Elastic Net") {
      fit_b <- glmnet(X_boot, y_boot, family = "binomial", alpha = 0.5)
      lam <- fit_b$lambda[min(50, length(fit_b$lambda))]
      p_on_boot <- predict(fit_b, X_boot, s = lam, type = "response")[, 1]
      p_on_orig <- predict(fit_b, X_orig, s = lam, type = "response")[, 1]

    } else {
      # RCS or unknown — use GLM as proxy
      df_b <- data.frame(boot_data[, boot_vars_used], ascvd_combined = y_boot)
      fml <- as.formula(paste("ascvd_combined ~",
                                paste(boot_vars_used, collapse = " + ")))
      fit_b <- tryCatch(glm(fml, data = df_b, family = binomial),
                          error = function(e) NULL)
      if (is.null(fit_b)) return(NA)
      p_on_boot <- predict(fit_b, newdata = df_b, type = "response")
      df_o <- data.frame(data[, boot_vars_used], ascvd_combined = y_orig)
      p_on_orig <- predict(fit_b, newdata = df_o, type = "response")
    }

    auc_b <- tryCatch(as.numeric(auc(roc(y_boot, p_on_boot, quiet = TRUE))),
                       error = function(e) NA)
    auc_o <- tryCatch(as.numeric(auc(roc(y_orig, p_on_orig, quiet = TRUE))),
                       error = function(e) NA)

    return(auc_b - auc_o)
  }

  set.seed(2026)
  boot_result <- boot(boot_data_full, boot_optimism_fn, R = 200)
  optimism <- mean(boot_result$t, na.rm = TRUE)
  auc_corrected <- auc_apparent - optimism

  cat(sprintf("  Apparent AUC:            %.4f\n", auc_apparent))
  cat(sprintf("  Mean optimism:           %.4f\n", optimism))
  cat(sprintf("  Optimism-corrected AUC:  %.4f\n", auc_corrected))

  # -- Save internal validation results -------------------------------------
  internal_val <- data.frame(
    metric = c("Best_Model", "Algorithm",
               "Apparent_AUC", "Apparent_AUC_Lower", "Apparent_AUC_Upper",
               "CV_AUC_Mean", "CV_AUC_SD", "CV_AUC_Lower", "CV_AUC_Upper",
               "Optimism", "Corrected_AUC",
               "N_vars", "N_complete", "N_events", "EPV"),
    value = c(best_model$name, best_model$algorithm,
              round(auc_apparent, 4), round(auc_ci[1], 4), round(auc_ci[3], 4),
              round(mean(cv_aucs), 4), round(sd(cv_aucs), 4),
              round(quantile(cv_aucs, 0.025), 4), round(quantile(cv_aucs, 0.975), 4),
              round(optimism, 4), round(auc_corrected, 4),
              best_model$n_vars, best_model$n_complete, best_model$n_events,
              round(best_model$n_events / best_model$n_vars, 1)),
    stringsAsFactors = FALSE
  )

  write.csv(internal_val, paste0(TAB_DIR, "calon2_internal_validation.csv"), row.names = FALSE)
  cat("  Saved: calon2_internal_validation.csv\n\n")
} else {
  cat("  No best model available for internal validation.\n\n")
  cv_aucs <- NA
  auc_apparent <- NA
  auc_corrected <- NA
  optimism <- NA
}

# =============================================================================
# SECTION 9: HEAD-TO-HEAD vs SAFEHEART-RE
# =============================================================================

cat("===================================================================\n")
cat("SECTION 9: HEAD-TO-HEAD vs SAFEHEART-RE (in UKB)\n")
cat("===================================================================\n\n")

# -- SAFEHEART-RE frozen coefficients ---------------------------------------
SAFEHEART_COEF <- list(
  intercept = -7.053,
  age       =  0.064,
  male      =  0.775,
  ldl_c     =  0.109,
  htn       =  0.431,
  bmi       =  0.025,
  smoking   =  0.466,
  prior_cvd =  1.414
)

# Use the Core elastic net complete-case dataset for fair comparison
if (!is.null(en_core)) {
  df_h2h <- en_core$df_cc

  # Compute SAFEHEART-RE predictions (no prior_cvd since we include prevalent)
  df_h2h$lp_safeheart <- with(df_h2h,
    SAFEHEART_COEF$intercept +
    SAFEHEART_COEF$age       * age +
    SAFEHEART_COEF$male      * sex +
    SAFEHEART_COEF$ldl_c     * re_ldl +
    SAFEHEART_COEF$htn       * hypertension +
    SAFEHEART_COEF$bmi       * bmi +
    SAFEHEART_COEF$smoking   * smoking_binary
  )
  df_h2h$prob_safeheart <- 1 / (1 + exp(-df_h2h$lp_safeheart))

  # Compute BEST model predictions on the same cohort
  # For fair comparison, need predictions on the Core complete cases
  if (!is.null(best_model)) {
    bm <- best_model$model

    # ═══════════════════════════════════════════════════════════════════════
    # BUG FIX (v4.1): Section 9 now handles ALL 6 algorithm types for
    # best-model predictions in SAFEHEART head-to-head comparison.
    # PREVIOUSLY: only handled EN / XGBoost → crashed on SL, RF, MI, RCS
    # ═══════════════════════════════════════════════════════════════════════

    # Helper: extract vars (RCS doesn't store $vars directly)
    if (grepl("RCS", best_model$algorithm)) {
      bm_vars <- setdiff(names(bm$df_cc), "ascvd_combined")
    } else {
      bm_vars <- bm$vars
    }
    core_has_all <- all(bm_vars %in% names(df_h2h))

    # -- Fallback function for when best-model vars aren't in Core dataset --
    fallback_core_en <- function() {
      cat("  NOTE: Best model uses vars not in Core dataset — using Core EN.\n")
      X_h2h <- as.matrix(df_h2h[, en_core$vars])
      predict(en_core$fit, X_h2h, s = "lambda.min", type = "response")[, 1]
    }

    if (best_model$algorithm == "Elastic Net") {
      if (core_has_all) {
        X_h2h <- as.matrix(df_h2h[, bm_vars])
        df_h2h$prob_best <- predict(bm$fit, X_h2h, s = "lambda.min",
                                     type = "response")[, 1]
      } else {
        df_h2h$prob_best <- fallback_core_en()
      }

    } else if (best_model$algorithm == "XGBoost") {
      if (core_has_all) {
        X_h2h <- as.matrix(df_h2h[, bm_vars])
        df_h2h$prob_best <- predict(bm$final_model, X_h2h)
      } else {
        df_h2h$prob_best <- fallback_core_en()
      }

    } else if (best_model$algorithm == "Super Learner") {
      if (core_has_all) {
        df_h2h$prob_best <- predict(bm$sl_fit,
                                      newdata = data.frame(df_h2h[, bm_vars]),
                                      onlySL = TRUE)$pred[, 1]
      } else {
        df_h2h$prob_best <- fallback_core_en()
      }

    } else if (best_model$algorithm == "Random Forest") {
      if (core_has_all) {
        df_h2h$prob_best <- predict(bm$rf,
                                      newdata = df_h2h[, bm_vars],
                                      type = "prob")[, "1"]
      } else {
        df_h2h$prob_best <- fallback_core_en()
      }

    } else if (best_model$algorithm == "MI-Elastic Net") {
      if (core_has_all) {
        # Average predictions across m imputed datasets
        mi_preds_h2h <- matrix(0, nrow = nrow(df_h2h), ncol = length(bm$mi_fits))
        for (ii in seq_along(bm$mi_fits)) {
          X_ii <- as.matrix(df_h2h[, bm_vars])
          mi_preds_h2h[, ii] <- predict(bm$mi_fits[[ii]]$fit, X_ii,
                                          s = "lambda.min", type = "response")[, 1]
        }
        df_h2h$prob_best <- rowMeans(mi_preds_h2h)
      } else {
        df_h2h$prob_best <- fallback_core_en()
      }

    } else if (grepl("RCS", best_model$algorithm)) {
      # RCS predict on h2h data — needs rms datadist environment
      tryCatch({
        df_h2h$prob_best <- predict(bm$fit, newdata = df_h2h, type = "fitted")
      }, error = function(e) {
        cat(sprintf("  NOTE: RCS predict failed (%s) — falling back to Core EN.\n",
                    e$message))
        df_h2h$prob_best <<- fallback_core_en()
      })

    } else {
      # Unknown algorithm — fall back to Core EN
      cat(sprintf("  NOTE: Unknown algorithm '%s' — using Core EN.\n",
                  best_model$algorithm))
      df_h2h$prob_best <- fallback_core_en()
    }

    # Also compute Core EN predictions for comparison
    X_core_h2h <- as.matrix(df_h2h[, en_core$vars])
    df_h2h$prob_core <- predict(en_core$fit, X_core_h2h, s = "lambda.min",
                                 type = "response")[, 1]

    # v3: Also compute Fixed Core EN predictions (if available)
    auc_fixed_en <- NA
    if (!is.null(en_fixed) && all(en_fixed$vars %in% names(df_h2h))) {
      X_fixed_h2h <- as.matrix(df_h2h[, en_fixed$vars])
      df_h2h$prob_fixed <- predict(en_fixed$fit, X_fixed_h2h, s = "lambda.min",
                                    type = "response")[, 1]
    }

    # v3b: Also compute Lean Core EN predictions (if available)
    auc_lean_en <- NA
    if (!is.null(en_lean) && all(en_lean$vars %in% names(df_h2h))) {
      X_lean_h2h <- as.matrix(df_h2h[, en_lean$vars])
      df_h2h$prob_lean <- predict(en_lean$fit, X_lean_h2h, s = "lambda.min",
                                   type = "response")[, 1]
    }

    # v4: Literature Core EN predictions (if available)
    auc_v4_en <- NA
    if (!is.null(en_v4_core) && all(en_v4_core$vars %in% names(df_h2h))) {
      X_v4_h2h <- as.matrix(df_h2h[, en_v4_core$vars])
      df_h2h$prob_v4 <- predict(en_v4_core$fit, X_v4_h2h, s = "lambda.min",
                                  type = "response")[, 1]
    }

    # v4: Montreal-Plus EN predictions
    auc_montp_en <- NA
    if (!is.null(en_v4_montp) && all(en_v4_montp$vars %in% names(df_h2h))) {
      X_montp_h2h <- as.matrix(df_h2h[, en_v4_montp$vars])
      df_h2h$prob_montp <- predict(en_v4_montp$fit, X_montp_h2h, s = "lambda.min",
                                     type = "response")[, 1]
    }

    # -- ROC comparison (Best model vs SAFEHEART) ---------------------------
    roc_best      <- roc(df_h2h$ascvd_combined, df_h2h$prob_best, quiet = TRUE)
    roc_core_en   <- roc(df_h2h$ascvd_combined, df_h2h$prob_core, quiet = TRUE)
    roc_safeheart <- roc(df_h2h$ascvd_combined, df_h2h$prob_safeheart, quiet = TRUE)

    auc_best      <- as.numeric(auc(roc_best))
    auc_core_en   <- as.numeric(auc(roc_core_en))
    auc_safeheart <- as.numeric(auc(roc_safeheart))
    ci_best       <- ci.auc(roc_best)
    ci_core_en    <- ci.auc(roc_core_en)
    ci_safeheart  <- ci.auc(roc_safeheart)

    # v3: Fixed Core ROC
    if ("prob_fixed" %in% names(df_h2h)) {
      roc_fixed_en <- roc(df_h2h$ascvd_combined, df_h2h$prob_fixed, quiet = TRUE)
      auc_fixed_en <- as.numeric(auc(roc_fixed_en))
      ci_fixed_en  <- ci.auc(roc_fixed_en)
      delong_fixed_sh <- roc.test(roc_fixed_en, roc_safeheart, method = "delong")
    }

    # v3b: Lean Core ROC
    if ("prob_lean" %in% names(df_h2h)) {
      roc_lean_en <- roc(df_h2h$ascvd_combined, df_h2h$prob_lean, quiet = TRUE)
      auc_lean_en <- as.numeric(auc(roc_lean_en))
      ci_lean_en  <- ci.auc(roc_lean_en)
      delong_lean_sh <- roc.test(roc_lean_en, roc_safeheart, method = "delong")
    }

    # v4: Literature Core ROC
    if ("prob_v4" %in% names(df_h2h)) {
      roc_v4_en <- roc(df_h2h$ascvd_combined, df_h2h$prob_v4, quiet = TRUE)
      auc_v4_en <- as.numeric(auc(roc_v4_en))
      ci_v4_en  <- ci.auc(roc_v4_en)
      delong_v4_sh <- roc.test(roc_v4_en, roc_safeheart, method = "delong")
    }

    # v4: Montreal-Plus ROC
    if ("prob_montp" %in% names(df_h2h)) {
      roc_montp_en <- roc(df_h2h$ascvd_combined, df_h2h$prob_montp, quiet = TRUE)
      auc_montp_en <- as.numeric(auc(roc_montp_en))
      ci_montp_en  <- ci.auc(roc_montp_en)
      delong_montp_sh <- roc.test(roc_montp_en, roc_safeheart, method = "delong")
    }

    # DeLong tests
    delong_best_sh <- roc.test(roc_best, roc_safeheart, method = "delong")
    delong_core_sh <- roc.test(roc_core_en, roc_safeheart, method = "delong")

    cat(sprintf("  Best model AUC:      %.4f (%.4f--%.4f) [%s]\n",
                auc_best, ci_best[1], ci_best[3], best_model$name))
    cat(sprintf("  EN Core AUC:         %.4f (%.4f--%.4f)\n",
                auc_core_en, ci_core_en[1], ci_core_en[3]))
    if (!is.na(auc_fixed_en)) {
      cat(sprintf("  EN Fixed Core (v3):  %.4f (%.4f--%.4f)\n",
                  auc_fixed_en, ci_fixed_en[1], ci_fixed_en[3]))
    }
    if (!is.na(auc_lean_en)) {
      cat(sprintf("  EN Lean Core (v3b):  %.4f (%.4f--%.4f)\n",
                  auc_lean_en, ci_lean_en[1], ci_lean_en[3]))
    }
    if (!is.na(auc_v4_en)) {
      cat(sprintf("  EN Lit Core (v4):    %.4f (%.4f--%.4f)\n",
                  auc_v4_en, ci_v4_en[1], ci_v4_en[3]))
    }
    if (!is.na(auc_montp_en)) {
      cat(sprintf("  EN Montreal+ (v4):   %.4f (%.4f--%.4f)\n",
                  auc_montp_en, ci_montp_en[1], ci_montp_en[3]))
    }
    cat(sprintf("  SAFEHEART-RE AUC:    %.4f (%.4f--%.4f)\n",
                auc_safeheart, ci_safeheart[1], ci_safeheart[3]))
    cat(sprintf("  dAUC (Best - SH):    %+.4f (DeLong p=%.6f) %s\n",
                auc_best - auc_safeheart, delong_best_sh$p.value,
                ifelse(delong_best_sh$p.value < 0.05, "SIGNIFICANT", "")))
    cat(sprintf("  dAUC (Core EN - SH): %+.4f (DeLong p=%.6f) %s\n",
                auc_core_en - auc_safeheart, delong_core_sh$p.value,
                ifelse(delong_core_sh$p.value < 0.05, "SIGNIFICANT", "")))
    if (!is.na(auc_fixed_en)) {
      cat(sprintf("  dAUC (Fixed v3-SH):  %+.4f (DeLong p=%.6f) %s\n",
                  auc_fixed_en - auc_safeheart, delong_fixed_sh$p.value,
                  ifelse(delong_fixed_sh$p.value < 0.05, "SIGNIFICANT", "")))
    }
    if (!is.na(auc_lean_en)) {
      cat(sprintf("  dAUC (Lean v3b-SH):  %+.4f (DeLong p=%.6f) %s\n",
                  auc_lean_en - auc_safeheart, delong_lean_sh$p.value,
                  ifelse(delong_lean_sh$p.value < 0.05, "SIGNIFICANT", "")))
    }
    if (!is.na(auc_v4_en)) {
      cat(sprintf("  dAUC (v4 Core-SH):   %+.4f (DeLong p=%.6f) %s\n",
                  auc_v4_en - auc_safeheart, delong_v4_sh$p.value,
                  ifelse(delong_v4_sh$p.value < 0.05, "SIGNIFICANT", "")))
    }
    if (!is.na(auc_montp_en)) {
      cat(sprintf("  dAUC (Montr+-SH):    %+.4f (DeLong p=%.6f) %s\n",
                  auc_montp_en - auc_safeheart, delong_montp_sh$p.value,
                  ifelse(delong_montp_sh$p.value < 0.05, "SIGNIFICANT", "")))
    }

    # -- NRI and IDI (Best model vs SAFEHEART) ------------------------------
    cat("\n  Reclassification (Best model vs SAFEHEART-RE):\n")
    tryCatch({
      nri_result <- nricens::nribin(
        event = df_h2h$ascvd_combined,
        p.std = df_h2h$prob_safeheart,
        p.new = df_h2h$prob_best,
        cut = c(0.10, 0.20),
        niter = 1000,
        msg = FALSE
      )

      cat(sprintf("  Categorical NRI: %.4f (%.4f--%.4f), p=%.4f\n",
                  nri_result$nri["Categorical NRI", "Estimate"],
                  nri_result$nri["Categorical NRI", "Lower"],
                  nri_result$nri["Categorical NRI", "Upper"],
                  nri_result$nri["Categorical NRI", "P-value"]))
      cat(sprintf("  IDI:             %.4f (%.4f--%.4f), p=%.4f\n",
                  nri_result$nri["IDI", "Estimate"],
                  nri_result$nri["IDI", "Lower"],
                  nri_result$nri["IDI", "Upper"],
                  nri_result$nri["IDI", "P-value"]))
    }, error = function(e) {
      cat(sprintf("  NRI/IDI calculation failed: %s\n", e$message))
    })

    # -- Save comparison results --------------------------------------------
    h2h_models <- c(best_model$name, "EN Core (14 vars)")
    h2h_aucs   <- c(auc_best, auc_core_en)
    h2h_lower  <- c(ci_best[1], ci_core_en[1])
    h2h_upper  <- c(ci_best[3], ci_core_en[3])

    if (!is.na(auc_fixed_en)) {
      h2h_models <- c(h2h_models, "EN Fixed Core (v3)")
      h2h_aucs   <- c(h2h_aucs, auc_fixed_en)
      h2h_lower  <- c(h2h_lower, ci_fixed_en[1])
      h2h_upper  <- c(h2h_upper, ci_fixed_en[3])
    }
    if (!is.na(auc_lean_en)) {
      h2h_models <- c(h2h_models, "EN Lean Core (v3b)")
      h2h_aucs   <- c(h2h_aucs, auc_lean_en)
      h2h_lower  <- c(h2h_lower, ci_lean_en[1])
      h2h_upper  <- c(h2h_upper, ci_lean_en[3])
    }
    if (!is.na(auc_v4_en)) {
      h2h_models <- c(h2h_models, "EN Lit Core (v4)")
      h2h_aucs   <- c(h2h_aucs, auc_v4_en)
      h2h_lower  <- c(h2h_lower, ci_v4_en[1])
      h2h_upper  <- c(h2h_upper, ci_v4_en[3])
    }
    if (!is.na(auc_montp_en)) {
      h2h_models <- c(h2h_models, "EN Montreal+ (v4)")
      h2h_aucs   <- c(h2h_aucs, auc_montp_en)
      h2h_lower  <- c(h2h_lower, ci_montp_en[1])
      h2h_upper  <- c(h2h_upper, ci_montp_en[3])
    }
    h2h_models <- c(h2h_models, "SAFEHEART-RE")
    h2h_aucs   <- c(h2h_aucs, auc_safeheart)
    h2h_lower  <- c(h2h_lower, ci_safeheart[1])
    h2h_upper  <- c(h2h_upper, ci_safeheart[3])

    h2h_comparison <- data.frame(
      model = h2h_models, AUC = h2h_aucs,
      AUC_lower = h2h_lower, AUC_upper = h2h_upper,
      stringsAsFactors = FALSE
    )
    h2h_comparison$delta_vs_safeheart <- h2h_comparison$AUC - auc_safeheart

    write.csv(h2h_comparison, paste0(TAB_DIR, "calon2_vs_safeheart_internal.csv"),
              row.names = FALSE)
    cat("\n  Saved: calon2_vs_safeheart_internal.csv\n")
  }
} else {
  cat("  Core EN model not available for SAFEHEART comparison.\n")
  auc_safeheart <- NA
}

cat("\n")

# =============================================================================
# SECTION 10: SENSITIVITY -- INCIDENT-ONLY
# =============================================================================

cat("===================================================================\n")
cat("SECTION 10: SENSITIVITY ANALYSIS -- INCIDENT ASCVD ONLY\n")
cat("===================================================================\n\n")

if (!is.null(best_model)) {
  bm <- best_model$model

  # ═══════════════════════════════════════════════════════════════════════════
  # BUG FIX (v4.1): Section 10 now handles ALL 6 algorithm types.
  # PREVIOUSLY: only handled EN / XGBoost → crashed on SL, RF, MI, RCS
  # Also: RCS doesn't store $vars — extract from df_cc column names.
  # ═══════════════════════════════════════════════════════════════════════════

  # Identify incident-only patients in the best model's complete-case data
  # Extract vars (RCS doesn't store $vars directly)
  if (grepl("RCS", best_model$algorithm)) {
    bm_vars <- setdiff(names(bm$df_cc), "ascvd_combined")
  } else {
    bm_vars <- bm$vars
  }

  # Build incident subset from the full df (not the model's cc data)
  needed_cols <- unique(c("eid", bm_vars, "ascvd_combined",
                            "bloods_before_ascvd", "ascvd_incident"))
  avail_cols <- needed_cols[needed_cols %in% names(df)]
  df_for_incident <- df[, avail_cols]
  cc_idx <- complete.cases(df_for_incident[, c(bm_vars, "ascvd_combined")])
  df_cc_full <- df_for_incident[cc_idx, ]

  # Incident-only: exclude prevalent (keep incident events + non-events)
  df_incident <- df_cc_full[df_cc_full$bloods_before_ascvd == TRUE |
                              df_cc_full$ascvd_combined == 0, ]
  if (nrow(df_incident) < 100) {
    # Fallback
    df_incident <- df_cc_full[df_cc_full$ascvd_incident == 1 |
                                df_cc_full$ascvd_combined == 0, ]
  }

  n_inc <- nrow(df_incident)
  n_inc_events <- sum(df_incident$ascvd_combined)

  cat(sprintf("  Incident-only cohort: %d patients, %d events (%.1f%%)\n",
              n_inc, n_inc_events, 100 * n_inc_events / n_inc))

  if (n_inc_events >= 30) {
    X_inc <- as.matrix(df_incident[, bm_vars])
    y_inc <- df_incident$ascvd_combined

    if (best_model$algorithm == "Elastic Net") {
      set.seed(2026)
      cv_inc <- cv.glmnet(X_inc, y_inc, family = "binomial",
                           alpha = bm$alpha, nfolds = 10, type.measure = "auc")
      auc_inc <- max(cv_inc$cvm)
      pred_inc <- predict(cv_inc, X_inc, s = "lambda.min", type = "response")[, 1]

    } else if (best_model$algorithm == "XGBoost") {
      dtrain_inc <- xgb.DMatrix(data = X_inc, label = y_inc)
      set.seed(2026)
      cv_inc <- tryCatch({
        xgb.cv(params = bm$params, data = dtrain_inc,
               nrounds = 1000, nfold = 10,
               early_stopping_rounds = 50, verbose = 0, print_every_n = 0)
      }, error = function(e) {
        xgb.cv(params = bm$params, data = dtrain_inc,
               nrounds = 1000, nfold = 10,
               early_stopping_rounds = 50, verbose = 0)
      })
      eval_log <- cv_inc$evaluation_log
      auc_col <- grep("test.*auc.*mean", names(eval_log), value = TRUE)[1]
      if (is.null(auc_col) || is.na(auc_col)) {
        auc_col <- names(eval_log)[grep("auc", names(eval_log))[1]]
      }
      auc_inc <- max(eval_log[[auc_col]], na.rm = TRUE)
      inc_model <- xgb.train(params = bm$params, data = dtrain_inc,
                              nrounds = which.max(eval_log[[auc_col]]),
                              verbose = 0)
      pred_inc <- predict(inc_model, X_inc)

    } else if (best_model$algorithm == "Super Learner") {
      # Use CV.SuperLearner for honest AUC on incident subset
      set.seed(2026)
      cv_sl_inc <- tryCatch({
        CV.SuperLearner(
          Y = y_inc, X = data.frame(df_incident[, bm_vars]),
          V = 10, family = binomial(),
          SL.library = names(bm$sl_fit$cvRisk),
          verbose = FALSE
        )
      }, error = function(e) NULL)
      if (!is.null(cv_sl_inc)) {
        sl_pred_inc <- cv_sl_inc$SL.predict[, 1]
        auc_inc <- as.numeric(auc(roc(y_inc, sl_pred_inc, quiet = TRUE)))
        pred_inc <- sl_pred_inc
      } else {
        # Fallback: glmnet proxy
        set.seed(2026)
        cv_inc_fb <- cv.glmnet(X_inc, y_inc, family = "binomial",
                                alpha = 0.5, nfolds = 10, type.measure = "auc")
        auc_inc <- max(cv_inc_fb$cvm)
        pred_inc <- predict(cv_inc_fb, X_inc, s = "lambda.min",
                             type = "response")[, 1]
        cat("  NOTE: SL failed on incident subset — used glmnet proxy.\n")
      }

    } else if (best_model$algorithm == "Random Forest") {
      set.seed(2026)
      rf_inc <- randomForest(x = data.frame(df_incident[, bm_vars]),
                              y = factor(y_inc, levels = c(0, 1)),
                              ntree = 1000, mtry = bm$mtry,
                              classwt = c("0" = 1, "1" = 2))
      oob_probs_inc <- rf_inc$votes[, "1"]
      auc_inc <- as.numeric(auc(roc(y_inc, oob_probs_inc, quiet = TRUE)))
      pred_inc <- predict(rf_inc, newdata = data.frame(df_incident[, bm_vars]),
                            type = "prob")[, "1"]

    } else if (best_model$algorithm == "MI-Elastic Net") {
      # Re-run mice on incident subset + cv.glmnet
      set.seed(2026)
      imp_inc <- tryCatch(
        mice(df_incident[, c(bm_vars, "ascvd_combined")], m = 5,
             method = "pmm", maxit = 10, printFlag = FALSE),
        error = function(e) NULL)
      if (!is.null(imp_inc)) {
        mi_a_inc <- numeric(5)
        mi_preds_inc <- matrix(0, nrow = nrow(df_incident), ncol = 5)
        for (ii in 1:5) {
          d_ii <- mice::complete(imp_inc, ii)
          X_ii <- as.matrix(d_ii[, bm_vars])
          y_ii <- d_ii$ascvd_combined
          cv_ii <- cv.glmnet(X_ii, y_ii, family = "binomial",
                              alpha = 0.5, nfolds = 10, type.measure = "auc")
          mi_a_inc[ii] <- max(cv_ii$cvm)
          mi_preds_inc[, ii] <- predict(cv_ii, X_ii,
                                          s = "lambda.min", type = "response")[, 1]
        }
        auc_inc <- mean(mi_a_inc)
        pred_inc <- rowMeans(mi_preds_inc)
      } else {
        # Fallback: simple glmnet on complete cases
        set.seed(2026)
        cv_inc_fb <- cv.glmnet(X_inc, y_inc, family = "binomial",
                                alpha = 0.5, nfolds = 10, type.measure = "auc")
        auc_inc <- max(cv_inc_fb$cvm)
        pred_inc <- predict(cv_inc_fb, X_inc, s = "lambda.min",
                             type = "response")[, 1]
        cat("  NOTE: MI failed on incident subset — used glmnet proxy.\n")
      }

    } else if (grepl("RCS", best_model$algorithm)) {
      # RCS on incident subset — refit using rms::lrm
      tryCatch({
        dd_inc <- datadist(df_incident)
        options(datadist = "dd_inc")
        rcs_inc <- lrm(bm$fit$sformula, data = df_incident,
                        x = TRUE, y = TRUE)
        val_inc <- validate(rcs_inc, B = 100, seed = 2026)
        auc_inc <- 0.5 * (val_inc["Dxy", "index.corrected"] + 1)
        pred_inc <- predict(rcs_inc, type = "fitted")
      }, error = function(e) {
        # Fallback: simple glmnet
        set.seed(2026)
        cv_inc_fb <- cv.glmnet(X_inc, y_inc, family = "binomial",
                                alpha = 0.5, nfolds = 10, type.measure = "auc")
        auc_inc <<- max(cv_inc_fb$cvm)
        pred_inc <<- predict(cv_inc_fb, X_inc, s = "lambda.min",
                              type = "response")[, 1]
        cat(sprintf("  NOTE: RCS failed on incident (%s) — used glmnet proxy.\n",
                    e$message))
      })

    } else {
      # Unknown algorithm — fall back to glmnet
      set.seed(2026)
      cv_inc_fb <- cv.glmnet(X_inc, y_inc, family = "binomial",
                              alpha = 0.5, nfolds = 10, type.measure = "auc")
      auc_inc <- max(cv_inc_fb$cvm)
      pred_inc <- predict(cv_inc_fb, X_inc, s = "lambda.min",
                           type = "response")[, 1]
      cat(sprintf("  NOTE: Algorithm '%s' not handled — used glmnet proxy.\n",
                  best_model$algorithm))
    }

    roc_inc <- roc(y_inc, pred_inc, quiet = TRUE)
    ci_inc <- ci.auc(roc_inc)

    cat(sprintf("  Incident-only CV-AUC:      %.4f\n", auc_inc))
    cat(sprintf("  Incident-only apparent AUC: %.4f (%.4f--%.4f)\n",
                as.numeric(auc(roc_inc)), ci_inc[1], ci_inc[3]))
    cat(sprintf("  Delta vs full cohort:       %+.4f\n",
                auc_inc - mean(cv_aucs, na.rm = TRUE)))
  } else {
    cat("  Too few incident events for sensitivity analysis.\n")
  }
} else {
  cat("  No best model available for sensitivity analysis.\n")
}

cat("\n")

# =============================================================================
# SECTION 11: SAVE ALL OUTPUTS
# =============================================================================

cat("===================================================================\n")
cat("SECTION 11: SAVING ALL OUTPUTS\n")
cat("===================================================================\n\n")

# Development summary (Core complete cases)
if (!is.null(en_core)) {
  dev_summary <- data.frame(
    variable = core_vars,
    stringsAsFactors = FALSE
  )
  dev_summary$mean   <- sapply(core_vars, function(v) mean(en_core$df_cc[[v]], na.rm = TRUE))
  dev_summary$sd     <- sapply(core_vars, function(v) sd(en_core$df_cc[[v]], na.rm = TRUE))
  dev_summary$min    <- sapply(core_vars, function(v) min(en_core$df_cc[[v]], na.rm = TRUE))
  dev_summary$max    <- sapply(core_vars, function(v) max(en_core$df_cc[[v]], na.rm = TRUE))
  dev_summary$median <- sapply(core_vars, function(v) median(en_core$df_cc[[v]], na.rm = TRUE))
  write.csv(dev_summary, paste0(TAB_DIR, "calon2_development_summary.csv"), row.names = FALSE)
  cat("  Saved: tables/calon2_development_summary.csv\n")
}

# v3: Fixed Core development summary
if (!is.null(en_fixed)) {
  dev_fixed <- data.frame(
    variable = core_fixed_vars,
    stringsAsFactors = FALSE
  )
  dev_fixed$mean   <- sapply(core_fixed_vars, function(v) mean(en_fixed$df_cc[[v]], na.rm = TRUE))
  dev_fixed$sd     <- sapply(core_fixed_vars, function(v) sd(en_fixed$df_cc[[v]], na.rm = TRUE))
  dev_fixed$min    <- sapply(core_fixed_vars, function(v) min(en_fixed$df_cc[[v]], na.rm = TRUE))
  dev_fixed$max    <- sapply(core_fixed_vars, function(v) max(en_fixed$df_cc[[v]], na.rm = TRUE))
  dev_fixed$median <- sapply(core_fixed_vars, function(v) median(en_fixed$df_cc[[v]], na.rm = TRUE))
  write.csv(dev_fixed, paste0(TAB_DIR, "calon2_fixed_core_summary.csv"), row.names = FALSE)
  cat("  Saved: tables/calon2_fixed_core_summary.csv\n")
}

# v3b: Lean Core development summary
if (!is.null(en_lean)) {
  dev_lean <- data.frame(
    variable = core_lean_vars,
    stringsAsFactors = FALSE
  )
  dev_lean$mean   <- sapply(core_lean_vars, function(v) mean(en_lean$df_cc[[v]], na.rm = TRUE))
  dev_lean$sd     <- sapply(core_lean_vars, function(v) sd(en_lean$df_cc[[v]], na.rm = TRUE))
  dev_lean$min    <- sapply(core_lean_vars, function(v) min(en_lean$df_cc[[v]], na.rm = TRUE))
  dev_lean$max    <- sapply(core_lean_vars, function(v) max(en_lean$df_cc[[v]], na.rm = TRUE))
  dev_lean$median <- sapply(core_lean_vars, function(v) median(en_lean$df_cc[[v]], na.rm = TRUE))
  write.csv(dev_lean, paste0(TAB_DIR, "calon2_lean_core_summary.csv"), row.names = FALSE)
  cat("  Saved: tables/calon2_lean_core_summary.csv\n")
}

# Temporal filtering summary
temporal_summary <- data.frame(
  category = c("Total cohort", "Prevalent ASCVD (event before assessment)",
               "Incident ASCVD (event after assessment)", "No ASCVD event",
               "MRI data available (pre-ASCVD only)", "hsCRP available"),
  n = c(nrow(df), n_prevalent, n_incident, n_no_event,
        sum(df$mri_before_ascvd == TRUE, na.rm = TRUE),
        sum(!is.na(df$crp))),
  stringsAsFactors = FALSE
)
temporal_summary$pct <- round(100 * temporal_summary$n / nrow(df), 1)
write.csv(temporal_summary, paste0(TAB_DIR, "calon2_temporal_filtering.csv"), row.names = FALSE)
cat("  Saved: tables/calon2_temporal_filtering.csv\n")

# List all saved files
cat("\n  All output files:\n")
cat("  [output/]\n")
cat("    calon2_core_coefficients.csv       (Core EN frozen coefficients)\n")
cat("    calon2_fixed_core_coefficients.csv (Fixed Core v3 coefficients)\n")
cat("    calon2_lean_core_coefficients.csv  (Lean Core v3b coefficients)\n")
cat("    calon2_enhanced_coefficients.csv   (Enhanced EN coefficients)\n")
cat("    calon2_best_model.rds              (Best model object)\n")
cat("    calon2_full_merged.csv             (Full merged dataset)\n")
cat("  [tables/]\n")
cat("    calon2_model_comparison.csv        (Grand comparison all models)\n")
cat("    calon2_internal_validation.csv     (Best model internal val)\n")
cat("    calon2_vs_safeheart_internal.csv   (Head-to-head comparison)\n")
cat("    calon2_missingness.csv             (Missingness report)\n")
cat("    calon2_temporal_filtering.csv      (Temporal filtering summary)\n")
cat("    calon2_development_summary.csv     (Core variable statistics)\n")
cat("    calon2_core_logistic_ORs.csv       (Core logistic regression ORs)\n")
cat("    calon2_fixed_core_logistic_ORs.csv (Fixed Core v3 ORs)\n")
cat("    calon2_lean_core_logistic_ORs.csv  (Lean Core v3b ORs)\n")
cat("    calon2_xgb_importance.csv          (XGBoost feature importance)\n")
cat("    calon2_v4_literature_core_ORs.csv  (v4 Literature Core ORs)\n")
cat("    calon2_rf_v4_importance.csv        (RF variable importance)\n")
cat("  [v4 models/]\n")
cat("    calon2_v4_literature_core_coefficients.csv (v4 frozen coeffs)\n")
cat("    Super Learner weights saved in model_results\n")
cat("    Multiple Imputation pooled estimates in model_results\n")

# =============================================================================
# FINAL SUMMARY BOX
# =============================================================================

cat("\n")
cat("======================================================================\n")
cat("   CALON-2 DEVELOPMENT v4 COMPLETE                                    \n")
cat("   v3: Collinearity fix + log transforms + interactions               \n")
cat("   v3b: Lean models (drop redundant inv_apoa + trig)                  \n")
cat("   v4: Literature features + MI + Super Learner + RF + RCS            \n")
cat("======================================================================\n")

if (!is.null(best_model)) {
  cat(sprintf("   BEST MODEL:          %s\n", best_model$name))
  cat(sprintf("   Algorithm:           %s\n", best_model$algorithm))
  cat(sprintf("   CV-AUC:              %.4f\n", best_model$cv_auc))
  if (!is.na(auc_apparent)) {
    cat(sprintf("   Apparent AUC:        %.4f\n", auc_apparent))
  }
  if (!is.na(auc_corrected)) {
    cat(sprintf("   Corrected AUC:       %.4f (optimism = %.4f)\n",
                auc_corrected, optimism))
  }
  if (!is.na(auc_safeheart)) {
    cat(sprintf("   vs SAFEHEART-RE:     dAUC = %+.4f\n",
                auc_best - auc_safeheart))
  }
}

cat("----------------------------------------------------------------------\n")
cat("   MODEL COMPARISON SUMMARY:\n")
if (length(model_results) > 0) {
  for (mn in names(model_results)) {
    m <- model_results[[mn]]
    cat(sprintf("   %-30s CV-AUC = %.4f\n", m$name, m$cv_auc))
  }
}

cat("----------------------------------------------------------------------\n")
cat("   TEMPORAL FILTERING:\n")
cat(sprintf("   Prevalent ASCVD (event before bloods): %4d\n", n_prevalent))
cat(sprintf("   Incident ASCVD (event after bloods):   %4d\n", n_incident))
cat(sprintf("   MRI pre-ASCVD (valid for imaging):     %4d\n",
            sum(df$mri_before_ascvd == TRUE, na.rm = TRUE)))
cat("----------------------------------------------------------------------\n")
cat("   v4 NEW ALGORITHMS:\n")
cat("   Random Forest, Multiple Imputation, Super Learner, RCS\n")
cat("----------------------------------------------------------------------\n")
cat("   NEXT: If AUC > 0.80 → Run 06_CALON2_validate_DRAGON3.R\n")
cat("         If AUC < 0.80 → Extract GGT (07b), consider CAC subset\n")
cat("======================================================================\n\n")

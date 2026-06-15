################################################################################
#                                                                              #
#  CALON-2: EXTERNAL VALIDATION ON DRAGON3 + WALES                            #
#  Frozen-Coefficient Validation (TRIPOD Type 4)                              #
#                                                                              #
#  Development:  UKB (N=1,623, 399 events) — Script 05                        #
#  Validation 1: DRAGON3 (N=1,362)                                            #
#  Validation 2: Wales ApoB-linked (N≈429)                                    #
#  Benchmark:    SAFEHEART-RE (published coefficients)                         #
#                                                                              #
#  Output: 8 Nature-quality figures + 6 manuscript tables                     #
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
  "dplyr", "tidyr", "readr", "ggplot2", "patchwork", "scales", "gridExtra",
  "pROC", "survival", "survminer", "rms",
  "boot", "ResourceSelection", "nricens", "tableone",
  "DescTools", "dcurves", "ggrepel", "cowplot", "grid"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("\n")
cat("\u2554\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2557\n")
cat("\u2551  CALON-2: EXTERNAL VALIDATION ON DRAGON3 + WALES                 \u2551\n")
cat("\u2551  TRIPOD Type 4 \u2014 Frozen Coefficients from UKB                    \u2551\n")
cat("\u2551  Goal: Demonstrate superiority over SAFEHEART-RE                  \u2551\n")
cat("\u255a\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u255d\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR <- "C:/Users/nader/Downloads/calon_ukb_pipeline/"
OUT_DIR  <- paste0(BASE_DIR, "output/")
FIG_DIR  <- paste0(OUT_DIR, "figures/")
TAB_DIR  <- paste0(OUT_DIR, "tables/")
for (d in c(FIG_DIR, TAB_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# External data paths
DRAGON3_FILE <- "C:/Users/nader/Downloads/DRAGON_3.csv"
WALES_CLEAN_FILE <- "C:/Users/nader/Downloads/WALES_FH_CLEANED.csv"
WALES_APOB_FILE  <- "C:/Users/nader/Downloads/ApoB_gl_Lpa_nmol_NHDL_1_first_LDL_u_TG_1_APOA_gl.csv"

# ── Nature formatting theme ──────────────────────────────────────────────────
theme_nature <- function(base_size = 7) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "Arial", color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 1,
                                hjust = 0, margin = margin(b = 6)),
      plot.subtitle = element_text(size = base_size, color = "grey30",
                                   margin = margin(b = 4)),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 0.5, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.size = unit(0.3, "cm"),
      legend.position = "bottom",
      legend.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# ── Colour palette ───────────────────────────────────────────────────────────
COL <- list(
  calon2     = "#2166AC",   # Deep blue (CALON-2)
  safeheart  = "#B2182B",   # Deep red (SAFEHEART-RE)
  dragon3    = "#4393C3",   # Medium blue (DRAGON3 cohort)
  wales      = "#35978F",   # Teal (Wales cohort)
  ukb        = "#762A83",   # Purple (UKB cohort)
  event      = "#D6604D",   # Coral red (ASCVD events)
  noevent    = "#92C5DE",   # Light blue (no events)
  male       = "#4393C3",
  female     = "#D6604D",
  ldlr       = "#2166AC",
  apob       = "#B2182B",
  pcsk9      = "#35978F",
  unknown    = "#999999",
  grey_dark  = "#636363",
  grey_light = "#BDBDBD",
  risk_low   = "#92C5DE",
  risk_med   = "#FDDBC7",
  risk_high  = "#D6604D"
)

# ── Figure dimensions (Nature: single 88mm, double 180mm) ────────────────────
WIDTH_SINGLE <- 88 / 25.4
WIDTH_DOUBLE <- 180 / 25.4
DPI_PDF <- 300
DPI_PNG <- 600

save_figure <- function(plot, name, width = WIDTH_DOUBLE, height = 5) {
  ggsave(paste0(FIG_DIR, name, ".pdf"), plot, width = width, height = height,
         dpi = DPI_PDF, device = cairo_pdf)
  ggsave(paste0(FIG_DIR, name, ".png"), plot, width = width, height = height,
         dpi = DPI_PNG, bg = "white")
  cat(sprintf("  Saved: %s.pdf + .png\n", name))
}

# =============================================================================
# SAFEHEART-RE FROZEN COEFFICIENTS (benchmark)
# =============================================================================

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

# =============================================================================
# CALIBRATION METRICS FUNCTION (reused from 02_CALON_validate.R)
# =============================================================================

calibration_metrics <- function(observed, predicted, model_name) {
  valid <- !is.na(observed) & !is.na(predicted) & !is.nan(predicted)
  if (sum(valid) < 20) {
    cat(sprintf("  %s: INSUFFICIENT VALID PREDICTIONS (%d)\n", model_name, sum(valid)))
    return(NULL)
  }
  obs  <- observed[valid]
  pred <- predicted[valid]
  pred <- pmax(pmin(pred, 1 - 1e-10), 1e-10)

  O <- sum(obs)
  E <- sum(pred)
  OE <- O / E

  recal <- tryCatch(
    glm(obs ~ offset(qlogis(pred)), family = binomial),
    error = function(e) NULL
  )
  intercept_recal <- if (!is.null(recal)) coef(recal)[1] else NA

  slope_model <- tryCatch(
    glm(obs ~ qlogis(pred), family = binomial),
    error = function(e) NULL
  )
  slope <- if (!is.null(slope_model)) coef(slope_model)[2] else NA

  hl <- tryCatch(
    hoslem.test(obs, pred, g = 10),
    error = function(e) {
      tryCatch(hoslem.test(obs, pred, g = 5),
               error = function(e2) list(statistic = NA, p.value = NA))
    }
  )

  brier <- mean((pred - obs)^2)
  brier_max <- mean(obs) * (1 - mean(obs))
  brier_scaled <- 1 - brier / brier_max

  cat(sprintf("  %s:\n", model_name))
  cat(sprintf("    N valid:            %d\n", sum(valid)))
  cat(sprintf("    O/E ratio:          %.3f (O=%d, E=%.1f)\n", OE, O, E))
  cat(sprintf("    Calibration slope:  %.3f (ideal=1.0)\n", slope))
  cat(sprintf("    Calibration int:    %.3f (ideal=0.0)\n", intercept_recal))
  if (!is.na(hl$statistic)) {
    cat(sprintf("    Hosmer-Lemeshow:    X2=%.2f, p=%.4f\n", hl$statistic, hl$p.value))
  } else {
    cat("    Hosmer-Lemeshow:    Could not compute\n")
  }
  cat(sprintf("    Brier score:        %.4f (scaled: %.3f)\n", brier, brier_scaled))

  return(data.frame(
    Model = model_name, N = sum(valid), Events = O,
    OE_ratio = OE, Cal_slope = slope, Cal_intercept = intercept_recal,
    HL_chi2 = ifelse(is.na(hl$statistic), NA, hl$statistic),
    HL_p = ifelse(is.na(hl$p.value), NA, hl$p.value),
    Brier = brier, Brier_scaled = brier_scaled
  ))
}

# =============================================================================
# SECTION 1: LOAD CALON-2 FROZEN COEFFICIENTS
# =============================================================================

cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("SECTION 1: LOADING CALON-2 FROZEN COEFFICIENTS\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

coef_file <- paste0(OUT_DIR, "calon2_core_coefficients.csv")
if (file.exists(coef_file)) {
  calon2_coef_df <- read.csv(coef_file, stringsAsFactors = FALSE)
  calon2_coef <- setNames(calon2_coef_df$coefficient, calon2_coef_df$variable)
  cat(sprintf("  Loaded %d CALON-2 Core coefficients from CSV\n", length(calon2_coef)))
  cat("\n  CALON-2 Core frozen coefficients:\n")
  for (nm in names(calon2_coef)) {
    cat(sprintf("    %-20s: %10.6f\n", nm, calon2_coef[nm]))
  }
} else {
  stop("ERROR: calon2_core_coefficients.csv not found!\n  Run 05_CALON2_develop.R first.")
}

# Also load internal validation metrics for comparison
internal_val_file <- paste0(TAB_DIR, "calon2_internal_validation.csv")
if (file.exists(internal_val_file)) {
  internal_val <- read.csv(internal_val_file, stringsAsFactors = FALSE)
  int_val <- setNames(internal_val$value, internal_val$metric)
  cat(sprintf("\n  Internal (UKB) CV-AUC: %.4f, Corrected AUC: %.4f\n",
              int_val["CV_AUC_Mean"], int_val["Corrected_AUC"]))
} else {
  cat("\n  Internal validation file not found. Will compute UKB metrics from scratch.\n")
  int_val <- NULL
}

# =============================================================================
# SECTION 2: LOAD & HARMONISE DRAGON3
# =============================================================================

cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("SECTION 2: LOADING & HARMONISING DRAGON3\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

d3 <- read.csv(DRAGON3_FILE, stringsAsFactors = FALSE)
cat(sprintf("  DRAGON3 loaded: %d patients x %d variables\n", nrow(d3), ncol(d3)))

# ── Harmonise Core variables to match UKB naming ─────────────────────────────

# Sex: M/F -> 1/0
d3$sex <- ifelse(d3$Gender == "M", 1,
           ifelse(d3$Gender == "F", 0, NA))

# Age
d3$age <- as.numeric(d3$Age_Treatment1)

# LDL (treatment-adjusted) — use LDL_UT if available, else derive
if ("LDL_UT" %in% names(d3)) {
  d3$re_ldl <- as.numeric(d3$LDL_UT)
} else {
  # Fallback: apply statin correction to LDL_1
  d3$re_ldl <- ifelse(!is.na(as.numeric(d3$Statin)) & as.numeric(d3$Statin) == 1,
                       as.numeric(d3$LDL_1) / 0.7, as.numeric(d3$LDL_1))
}

# HDL
d3$hdl <- as.numeric(d3$HDL_1)

# Triglycerides
d3$trig <- as.numeric(d3$TRG_1)

# ApoB
d3$apob <- as.numeric(d3$ApoB)

# Inverse ApoA
if ("inverse_APOA" %in% names(d3)) {
  d3$inv_apoa <- as.numeric(d3$inverse_APOA)
} else if ("ApoA1" %in% names(d3)) {
  apoa_val <- as.numeric(d3$ApoA1)
  d3$inv_apoa <- ifelse(!is.na(apoa_val) & apoa_val > 0, 1 / apoa_val, NA)
} else {
  d3$inv_apoa <- NA
}

# ApoB/LDL ratio
if ("APOB_LDL" %in% names(d3)) {
  d3$apob_ldl_ratio <- as.numeric(d3$APOB_LDL)
} else {
  d3$apob_ldl_ratio <- ifelse(!is.na(d3$apob) & !is.na(d3$re_ldl) & d3$re_ldl > 0,
                               d3$apob / d3$re_ldl, NA)
}

# Lp(a) binary (>50 nmol/L or >30 mg/dL)
lpa_val <- as.numeric(d3$Lpa)
d3$lpa_binary <- ifelse(!is.na(lpa_val) & lpa_val > 50, 1, 0)
# Handle missing Lp(a): set to NA (conservative)
d3$lpa_binary[is.na(lpa_val)] <- NA

# Smoking
d3$smoking_binary <- as.numeric(d3$Smoking_binary)

# Diabetes
d3$diabetes <- as.numeric(d3$Diabetes_binary)

# Hypertension
d3$hypertension <- as.numeric(d3$BloodPressureMedication)

# BMI
d3$bmi <- as.numeric(d3$BMI)

# Outcome
d3$ascvd_combined <- as.numeric(d3$ASCVD_combined)

# ── Gene type from Mutation1 text ─────────────────────────────────────────────
d3$gene <- case_when(
  grepl("LDLR", d3$Mutation1, ignore.case = TRUE) & d3$Positive1 == 1 ~ "LDLR",
  grepl("APOB", d3$Mutation1, ignore.case = TRUE) & d3$Positive1 == 1 ~ "APOB",
  grepl("PCSK9", d3$Mutation1, ignore.case = TRUE) & d3$Positive1 == 1 ~ "PCSK9",
  TRUE ~ "Unknown"
)

# ── Statin status ────────────────────────────────────────────────────────────
if ("Statin" %in% names(d3)) {
  d3$on_statin <- as.numeric(d3$Statin)
} else {
  d3$on_statin <- NA
}

# ── ASCVD prevalent (for SAFEHEART-RE prior_cvd term) ─────────────────────────
# Use ASCVD_P (prevalent) if available, else set to 0
if ("ASCVD_P" %in% names(d3)) {
  d3$ascvd_prevalent <- as.numeric(d3$ASCVD_P)
} else if ("ASCVD_prevalent" %in% names(d3)) {
  d3$ascvd_prevalent <- as.numeric(d3$ASCVD_prevalent)
} else {
  d3$ascvd_prevalent <- 0
  cat("  Note: No prevalent ASCVD column found in DRAGON3. Set to 0 for SAFEHEART-RE.\n")
}

# ── Link Deprivation from Wales Clean via DatabaseNumber ──────────────────────
cat("\n  Linking deprivation from Wales Clean...\n")
wales <- read.csv(WALES_CLEAN_FILE, stringsAsFactors = FALSE)
cat(sprintf("  Wales Clean loaded: %d patients\n", nrow(wales)))

# Select WIMD columns
wimd_cols <- c("DatabaseNumber")
# Find WIMD columns dynamically
wimd_candidates <- grep("^WIMD", names(wales), value = TRUE)
if (length(wimd_candidates) > 0) {
  # Prefer WIMD2025OverallDecile or WIMD or WIMD2025LSOARank
  if ("WIMD2025OverallDecile" %in% wimd_candidates) {
    wimd_cols <- c(wimd_cols, "WIMD2025OverallDecile")
  }
  if ("WIMD" %in% wimd_candidates) {
    wimd_cols <- c(wimd_cols, "WIMD")
  }
  if ("WIMD2025LSOARank" %in% wimd_candidates) {
    wimd_cols <- c(wimd_cols, "WIMD2025LSOARank")
  }
}

wimd_df <- wales[, wimd_cols[wimd_cols %in% names(wales)], drop = FALSE]
wimd_df <- wimd_df[!duplicated(wimd_df$DatabaseNumber), ]

d3 <- merge(d3, wimd_df, by = "DatabaseNumber", all.x = TRUE)

# Convert WIMD rank to centile (1917 LSOAs in Wales)
if ("WIMD" %in% names(d3)) {
  d3$deprivation_centile <- (d3$WIMD / 1917) * 100
  cat(sprintf("  WIMD linked: %d non-missing (%.1f%%)\n",
              sum(!is.na(d3$deprivation_centile)),
              100 * mean(!is.na(d3$deprivation_centile))))
} else if ("WIMD2025LSOARank" %in% names(d3)) {
  d3$deprivation_centile <- (d3$WIMD2025LSOARank / 1917) * 100
  cat(sprintf("  WIMD2025 linked: %d non-missing (%.1f%%)\n",
              sum(!is.na(d3$deprivation_centile)),
              100 * mean(!is.na(d3$deprivation_centile))))
} else {
  d3$deprivation_centile <- NA
  cat("  Warning: No WIMD column found in Wales Clean\n")
}

# Deprivation decile if available
if ("WIMD2025OverallDecile" %in% names(d3)) {
  d3$deprivation_decile <- d3$WIMD2025OverallDecile
} else if (!is.na(d3$deprivation_centile[1])) {
  d3$deprivation_decile <- ceiling(d3$deprivation_centile / 10)
} else {
  d3$deprivation_decile <- NA
}

# ── hsCRP (log-transformed) — check if DRAGON3 has CRP data ────────────────
crp_candidates <- c("CRP", "crp", "hsCRP", "HSCRP", "CRP_1", "CRP1")
crp_found <- FALSE
for (cc in crp_candidates) {
  if (cc %in% names(d3) && sum(!is.na(d3[[cc]])) > 10) {
    d3$crp <- as.numeric(d3[[cc]])
    d3$log_crp <- log(d3$crp + 0.01)
    crp_found <- TRUE
    cat(sprintf("\n  hsCRP found in DRAGON3 as '%s': %d non-missing\n", cc,
                sum(!is.na(d3$crp))))
    break
  }
}
if (!crp_found) {
  d3$log_crp <- NA
  cat("\n  NOTE: hsCRP not found in DRAGON3.\n")
  cat("        CALON-2 Core model log_crp coefficient will use median imputation.\n")
  cat("        Sensitivity analysis: Core-13 (without hsCRP) will also be reported.\n")
}

# Print DRAGON3 harmonisation summary
cat("\n  DRAGON3 harmonisation summary:\n")
core_vars <- c("age", "sex", "apob", "re_ldl", "hdl", "trig",
               "inv_apoa", "apob_ldl_ratio", "lpa_binary",
               "smoking_binary", "diabetes", "hypertension", "bmi",
               "log_crp")
for (v in c(core_vars, "ascvd_combined", "deprivation_centile")) {
  n_valid <- sum(!is.na(d3[[v]]))
  cat(sprintf("    %-20s: %d / %d (%.1f%%)\n", v, n_valid, nrow(d3),
              100 * n_valid / nrow(d3)))
}

cat(sprintf("\n  ASCVD events: %d / %d (%.1f%%)\n",
            sum(d3$ascvd_combined, na.rm = TRUE), nrow(d3),
            100 * mean(d3$ascvd_combined, na.rm = TRUE)))

cat(sprintf("  Gene distribution: LDLR=%d, APOB=%d, PCSK9=%d, Unknown=%d\n",
            sum(d3$gene == "LDLR"), sum(d3$gene == "APOB"),
            sum(d3$gene == "PCSK9"), sum(d3$gene == "Unknown")))

# =============================================================================
# SECTION 3: LOAD & HARMONISE WALES ApoB COHORT
# =============================================================================

cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("SECTION 3: LOADING & HARMONISING WALES ApoB COHORT\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

wales_apob_valid <- FALSE

if (file.exists(WALES_APOB_FILE)) {
  apob_file <- read.csv(WALES_APOB_FILE, stringsAsFactors = FALSE)
  cat(sprintf("  Wales ApoB file loaded: %d patients x %d variables\n",
              nrow(apob_file), ncol(apob_file)))

  # Check for DatabaseNumber to link back to Wales Clean
  if ("DatabaseNumber" %in% names(apob_file)) {
    # Merge ApoB file with Wales Clean for full variable set
    wl <- merge(wales, apob_file, by = "DatabaseNumber", all.x = FALSE, all.y = TRUE)
    cat(sprintf("  Merged with Wales Clean: %d patients\n", nrow(wl)))
  } else {
    # If no DatabaseNumber, use the ApoB file standalone
    wl <- apob_file
    cat("  Note: No DatabaseNumber in ApoB file. Using standalone.\n")
  }

  # ── Harmonise Wales variables ──────────────────────────────────────────────

  # Sex
  if ("Gender" %in% names(wl)) {
    wl$sex <- ifelse(wl$Gender == "M", 1,
               ifelse(wl$Gender == "F", 0, NA))
  } else if ("gender" %in% names(wl)) {
    wl$sex <- ifelse(tolower(wl$gender) == "m", 1,
               ifelse(tolower(wl$gender) == "f", 0, NA))
  }

  # Age
  if ("Ageattest" %in% names(wl)) {
    wl$age <- as.numeric(wl$Ageattest)
  } else if ("Age_Treatment1" %in% names(wl)) {
    wl$age <- as.numeric(wl$Age_Treatment1)
  }

  # ApoB — try multiple column names
  apob_cols <- c("ApoB_gl", "ApoB", "apob", "APOB")
  for (ac in apob_cols) {
    if (ac %in% names(wl) && sum(!is.na(wl[[ac]])) > 0) {
      wl$apob <- as.numeric(wl[[ac]])
      break
    }
  }

  # ApoA1 / inverse ApoA
  apoa_cols <- c("APOA_gl", "ApoA1", "apoa", "APOA")
  for (ac in apoa_cols) {
    if (ac %in% names(wl) && sum(!is.na(wl[[ac]])) > 0) {
      wl$apoa <- as.numeric(wl[[ac]])
      wl$inv_apoa <- ifelse(!is.na(wl$apoa) & wl$apoa > 0, 1 / wl$apoa, NA)
      break
    }
  }

  # LDL (treatment-adjusted)
  if ("LDL_u" %in% names(wl)) {
    # LDL_u from ApoB file = untreated LDL (already corrected)
    wl$re_ldl <- as.numeric(wl$LDL_u)
  } else if ("first_LDL_u" %in% names(wl)) {
    wl$re_ldl <- as.numeric(wl$first_LDL_u)
  } else if ("LDL.1" %in% names(wl)) {
    # Apply statin correction
    on_statin <- !is.na(wl$Statin) & wl$Statin == 1
    wl$re_ldl <- ifelse(on_statin, wl$LDL.1 / 0.7, wl$LDL.1)
  }

  # HDL
  if ("HDL.1" %in% names(wl)) {
    wl$hdl <- wl$HDL.1
  } else if ("HDL_1" %in% names(wl)) {
    wl$hdl <- wl$HDL_1
  }

  # Triglycerides
  if ("TRG.1" %in% names(wl)) {
    wl$trig <- wl$TRG.1
  } else if ("TG_1" %in% names(wl)) {
    wl$trig <- as.numeric(wl$TG_1)
  } else if ("TRG_1" %in% names(wl)) {
    wl$trig <- wl$TRG_1
  }

  # ApoB/LDL ratio
  wl$apob_ldl_ratio <- ifelse(!is.na(wl$apob) & !is.na(wl$re_ldl) & wl$re_ldl > 0,
                               wl$apob / wl$re_ldl, NA)

  # Lp(a)
  lpa_cols <- c("Lpa_nmol", "Lpa.1", "Lpa", "lpa")
  for (lc in lpa_cols) {
    if (lc %in% names(wl) && sum(!is.na(wl[[lc]])) > 0) {
      wl$lpa_raw <- as.numeric(wl[[lc]])
      break
    }
  }
  if ("lpa_raw" %in% names(wl)) {
    wl$lpa_binary <- ifelse(!is.na(wl$lpa_raw) & wl$lpa_raw > 50, 1, 0)
    wl$lpa_binary[is.na(wl$lpa_raw)] <- NA
  } else {
    wl$lpa_binary <- NA
  }

  # Clinical variables
  if ("Smoking" %in% names(wl)) {
    wl$smoking_binary <- as.numeric(wl$Smoking)
  } else if ("Smoking_binary" %in% names(wl)) {
    wl$smoking_binary <- as.numeric(wl$Smoking_binary)
  }

  if ("Diabetes" %in% names(wl)) {
    wl$diabetes <- as.numeric(wl$Diabetes)
  } else if ("Diabetes_binary" %in% names(wl)) {
    wl$diabetes <- as.numeric(wl$Diabetes_binary)
  }

  if ("BloodPressureMedication" %in% names(wl)) {
    wl$hypertension <- as.numeric(wl$BloodPressureMedication)
  }

  if ("BMI" %in% names(wl)) {
    wl$bmi <- as.numeric(wl$BMI)
  }

  # Outcome
  if ("ASCVD_combined" %in% names(wl)) {
    wl$ascvd_combined <- as.numeric(wl$ASCVD_combined)
  } else if ("ASCVD_Combined" %in% names(wl)) {
    wl$ascvd_combined <- as.numeric(wl$ASCVD_Combined)
  } else {
    # Derive from component fields if present
    wl$ascvd_combined <- NA
    cat("  Warning: No ASCVD outcome in Wales ApoB file\n")
  }

  # Prevalent ASCVD for SAFEHEART-RE
  if ("ASCVD_P" %in% names(wl)) {
    wl$ascvd_prevalent <- as.numeric(wl$ASCVD_P)
  } else {
    wl$ascvd_prevalent <- 0
  }

  # WIMD deprivation from Wales Clean
  if ("WIMD" %in% names(wl)) {
    wl$deprivation_centile <- (wl$WIMD / 1917) * 100
  } else if ("WIMD2025LSOARank" %in% names(wl)) {
    wl$deprivation_centile <- (wl$WIMD2025LSOARank / 1917) * 100
  }

  # hsCRP — check if available in Wales data
  crp_wl_cols <- c("CRP", "crp", "hsCRP", "HSCRP", "CRP_1", "CRP.1")
  wl$log_crp <- NA
  for (cc in crp_wl_cols) {
    if (cc %in% names(wl) && sum(!is.na(wl[[cc]])) > 10) {
      wl$crp <- as.numeric(wl[[cc]])
      wl$log_crp <- log(wl$crp + 0.01)
      cat(sprintf("  Wales hsCRP found as '%s': %d non-missing\n",
                  cc, sum(!is.na(wl$crp))))
      break
    }
  }
  if (all(is.na(wl$log_crp))) {
    cat("  Wales hsCRP not found — will use UKB median imputation\n")
  }

  # Check data adequacy
  n_apob_valid <- sum(!is.na(wl$apob))
  n_outcome    <- sum(!is.na(wl$ascvd_combined))
  n_events_wl  <- sum(wl$ascvd_combined, na.rm = TRUE)

  cat(sprintf("\n  Wales ApoB cohort summary:\n"))
  cat(sprintf("    N total:       %d\n", nrow(wl)))
  cat(sprintf("    ApoB valid:    %d (%.1f%%)\n", n_apob_valid, 100 * n_apob_valid / nrow(wl)))
  cat(sprintf("    Outcome valid: %d (%.1f%%)\n", n_outcome, 100 * n_outcome / nrow(wl)))
  cat(sprintf("    ASCVD events:  %d\n", n_events_wl))

  for (v in core_vars) {
    if (v %in% names(wl)) {
      n_v <- sum(!is.na(wl[[v]]))
      cat(sprintf("    %-20s: %d / %d (%.1f%%)\n", v, n_v, nrow(wl), 100 * n_v / nrow(wl)))
    } else {
      cat(sprintf("    %-20s: NOT FOUND\n", v))
    }
  }

  if (n_apob_valid >= 50 && n_events_wl >= 10) {
    wales_apob_valid <- TRUE
    cat("\n  Wales ApoB cohort is ADEQUATE for secondary validation.\n")
  } else {
    cat("\n  Wales ApoB cohort has INSUFFICIENT data. Will report descriptive only.\n")
  }
} else {
  cat("  Wales ApoB file not found. Skipping secondary validation.\n")
  wl <- NULL
}

# =============================================================================
# SECTION 4: LOAD UKB DEVELOPMENT DATA (for Table 1 comparison)
# =============================================================================

cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("SECTION 4: LOADING UKB DEVELOPMENT DATA\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

ukb <- read.csv(paste0(OUT_DIR, "calon_ukb_analysis_ready.csv"), stringsAsFactors = FALSE)
cat(sprintf("  UKB loaded: %d patients\n", nrow(ukb)))

# Create log_crp in UKB (hsCRP is stored as 'crp' in base dataset)
if ("crp" %in% names(ukb) && sum(!is.na(ukb$crp)) > 100) {
  ukb$log_crp <- log(ukb$crp + 0.01)
  cat(sprintf("  UKB log(hsCRP) created: %d non-missing\n", sum(!is.na(ukb$log_crp))))
} else {
  ukb$log_crp <- NA
  cat("  UKB: hsCRP not available\n")
}

# =============================================================================
# SECTION 5: APPLY FROZEN MODELS — COMPUTE PREDICTED PROBABILITIES
# =============================================================================

cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("SECTION 5: APPLYING FROZEN MODELS\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

# ── Helper: compute linear predictor from frozen elastic net coefficients ─────
# impute_values: named list of median imputation values from development cohort
#                for variables not available in validation data (e.g., log_crp)
apply_calon2 <- function(data, coefs, impute_values = NULL) {
  # The coefficient file has "(Intercept)" and variable names
  intercept <- coefs["(Intercept)"]
  if (is.na(intercept)) intercept <- 0

  var_coefs <- coefs[names(coefs) != "(Intercept)"]

  lp <- rep(intercept, nrow(data))
  vars_used <- c()
  vars_imputed <- c()
  vars_missing <- c()

  for (v in names(var_coefs)) {
    if (var_coefs[v] == 0) next  # Skip zero coefficients

    if (v %in% names(data)) {
      val <- as.numeric(data[[v]])
      # For columns that exist but are mostly NA, use imputation if available
      if (all(is.na(val)) && !is.null(impute_values) && v %in% names(impute_values)) {
        lp <- lp + var_coefs[v] * impute_values[[v]]
        vars_imputed <- c(vars_imputed, v)
      } else {
        # Use median imputation for individual NAs within a partially-observed column
        median_val <- median(val, na.rm = TRUE)
        if (is.na(median_val)) median_val <- 0
        lp <- lp + var_coefs[v] * ifelse(is.na(val), median_val, val)
        vars_used <- c(vars_used, v)
      }
    } else if (!is.null(impute_values) && v %in% names(impute_values)) {
      # Variable not in data at all — use development median
      lp <- lp + var_coefs[v] * impute_values[[v]]
      vars_imputed <- c(vars_imputed, v)
    } else {
      cat(sprintf("    Warning: Variable '%s' not in data and no imputation, coef %.6f IGNORED\n",
                  v, var_coefs[v]))
      vars_missing <- c(vars_missing, v)
    }
  }

  if (length(vars_imputed) > 0) {
    cat(sprintf("    Variables imputed with development median: %s\n",
                paste(vars_imputed, collapse = ", ")))
  }
  if (length(vars_missing) > 0) {
    cat(sprintf("    Variables MISSING (no imputation): %s\n",
                paste(vars_missing, collapse = ", ")))
  }

  prob <- 1 / (1 + exp(-lp))
  return(list(lp = lp, prob = prob, vars_used = vars_used,
              vars_imputed = vars_imputed, vars_missing = vars_missing))
}

apply_safeheart <- function(data) {
  lp <- with(data,
    SAFEHEART_COEF$intercept +
    SAFEHEART_COEF$age       * age +
    SAFEHEART_COEF$male      * sex +
    SAFEHEART_COEF$ldl_c     * re_ldl +
    SAFEHEART_COEF$htn       * hypertension +
    SAFEHEART_COEF$bmi       * bmi +
    SAFEHEART_COEF$smoking   * smoking_binary +
    SAFEHEART_COEF$prior_cvd * ascvd_prevalent
  )
  prob <- 1 / (1 + exp(-lp))
  return(list(lp = lp, prob = prob))
}

# ── Compute development-data medians for imputation ─────────────────────────
# When a validation cohort is missing a variable (e.g., log_crp not in DRAGON3),
# impute with the UKB development median — effectively "centres" the missing
# variable's contribution, minimising bias.
cat("  Computing UKB development medians for imputation...\n")
ukb_impute <- list()
for (v in names(calon2_coef)) {
  if (v != "(Intercept)" && v %in% names(ukb) && calon2_coef[v] != 0) {
    med <- median(as.numeric(ukb[[v]]), na.rm = TRUE)
    if (!is.na(med)) {
      ukb_impute[[v]] <- med
      cat(sprintf("    %s: median = %.4f\n", v, med))
    }
  }
}

# ── Apply to DRAGON3 ─────────────────────────────────────────────────────────
cat("\n  Applying models to DRAGON3...\n")
d3_calon2 <- apply_calon2(d3, calon2_coef, impute_values = ukb_impute)
d3$prob_calon2 <- d3_calon2$prob
d3$lp_calon2   <- d3_calon2$lp

d3_sh <- apply_safeheart(d3)
d3$prob_safeheart <- d3_sh$prob
d3$lp_safeheart   <- d3_sh$lp

cat(sprintf("  DRAGON3: %d CALON-2 predictions, %d SAFEHEART-RE predictions\n",
            sum(!is.na(d3$prob_calon2)), sum(!is.na(d3$prob_safeheart))))

# Track which variables were imputed for reporting
d3_vars_imputed <- d3_calon2$vars_imputed
d3_vars_missing <- d3_calon2$vars_missing

# ── Apply to Wales ApoB ──────────────────────────────────────────────────────
if (wales_apob_valid && !is.null(wl)) {
  cat("\n  Applying models to Wales ApoB cohort...\n")
  wl_calon2 <- apply_calon2(wl, calon2_coef, impute_values = ukb_impute)
  wl$prob_calon2 <- wl_calon2$prob
  wl$lp_calon2   <- wl_calon2$lp

  wl_sh <- apply_safeheart(wl)
  wl$prob_safeheart <- wl_sh$prob
  wl$lp_safeheart   <- wl_sh$lp

  cat(sprintf("  Wales: %d CALON-2 predictions, %d SAFEHEART-RE predictions\n",
              sum(!is.na(wl$prob_calon2)), sum(!is.na(wl$prob_safeheart))))
}

# ── Apply to UKB (re-compute for consistency) ────────────────────────────────
cat("\n  Applying models to UKB (internal comparison)...\n")
ukb_calon2 <- apply_calon2(ukb, calon2_coef, impute_values = ukb_impute)
ukb$prob_calon2 <- ukb_calon2$prob
ukb$lp_calon2   <- ukb_calon2$lp

# UKB needs ascvd_prevalent for SAFEHEART-RE
if (!"ascvd_prevalent" %in% names(ukb)) {
  ukb$ascvd_prevalent <- 0
}
ukb_sh <- apply_safeheart(ukb)
ukb$prob_safeheart <- ukb_sh$prob
ukb$lp_safeheart   <- ukb_sh$lp

# =============================================================================
# SECTION 6: DISCRIMINATION — AUC + DeLong
# =============================================================================

cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("SECTION 6: DISCRIMINATION (AUC)\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

disc_results <- list()

# ── Helper: compute AUC + DeLong between two models ──────────────────────────
compute_discrimination <- function(outcome, prob_new, prob_ref, cohort_name) {
  valid <- !is.na(outcome) & !is.na(prob_new) & !is.na(prob_ref)
  y    <- outcome[valid]
  pnew <- prob_new[valid]
  pref <- prob_ref[valid]

  if (sum(valid) < 30 || sum(y) < 5) {
    cat(sprintf("  %s: INSUFFICIENT DATA (n=%d, events=%d)\n", cohort_name, sum(valid), sum(y)))
    return(NULL)
  }

  roc_new <- roc(y, pnew, quiet = TRUE)
  roc_ref <- roc(y, pref, quiet = TRUE)

  auc_new <- as.numeric(auc(roc_new))
  auc_ref <- as.numeric(auc(roc_ref))
  ci_new  <- ci.auc(roc_new)
  ci_ref  <- ci.auc(roc_ref)

  delong <- roc.test(roc_new, roc_ref, method = "delong")

  cat(sprintf("  %s (n=%d, events=%d):\n", cohort_name, sum(valid), sum(y)))
  cat(sprintf("    CALON-2 AUC:      %.4f (%.4f\u2013%.4f)\n", auc_new, ci_new[1], ci_new[3]))
  cat(sprintf("    SAFEHEART-RE AUC: %.4f (%.4f\u2013%.4f)\n", auc_ref, ci_ref[1], ci_ref[3]))
  cat(sprintf("    \u0394AUC:              %+.4f\n", auc_new - auc_ref))
  cat(sprintf("    DeLong p-value:   %.6f %s\n\n", delong$p.value,
              ifelse(delong$p.value < 0.05, "\u2605 SIGNIFICANT", "")))

  return(data.frame(
    Cohort = cohort_name, N = sum(valid), Events = sum(y),
    CALON2_AUC = auc_new, CALON2_Lower = ci_new[1], CALON2_Upper = ci_new[3],
    SAFEHEART_AUC = auc_ref, SAFEHEART_Lower = ci_ref[1], SAFEHEART_Upper = ci_ref[3],
    Delta_AUC = auc_new - auc_ref, DeLong_p = delong$p.value,
    stringsAsFactors = FALSE
  ))
}

# ── UKB (internal / apparent) ────────────────────────────────────────────────
disc_results[[1]] <- compute_discrimination(
  ukb$ascvd_combined, ukb$prob_calon2, ukb$prob_safeheart, "UKB (Development)")

# ── DRAGON3 (primary external validation) ────────────────────────────────────
disc_results[[2]] <- compute_discrimination(
  d3$ascvd_combined, d3$prob_calon2, d3$prob_safeheart, "DRAGON3 (External)")

# ── Wales ApoB (secondary external validation) ──────────────────────────────
if (wales_apob_valid && !is.null(wl)) {
  disc_results[[3]] <- compute_discrimination(
    wl$ascvd_combined, wl$prob_calon2, wl$prob_safeheart, "Wales ApoB (External)")
}

# Store ROC objects for plotting
roc_d3_calon2 <- roc(d3$ascvd_combined[!is.na(d3$prob_calon2) & !is.na(d3$ascvd_combined)],
                      d3$prob_calon2[!is.na(d3$prob_calon2) & !is.na(d3$ascvd_combined)], quiet = TRUE)
roc_d3_sh     <- roc(d3$ascvd_combined[!is.na(d3$prob_safeheart) & !is.na(d3$ascvd_combined)],
                      d3$prob_safeheart[!is.na(d3$prob_safeheart) & !is.na(d3$ascvd_combined)], quiet = TRUE)
roc_ukb_calon2 <- roc(ukb$ascvd_combined[!is.na(ukb$prob_calon2)],
                       ukb$prob_calon2[!is.na(ukb$prob_calon2)], quiet = TRUE)
roc_ukb_sh     <- roc(ukb$ascvd_combined[!is.na(ukb$prob_safeheart)],
                       ukb$prob_safeheart[!is.na(ukb$prob_safeheart)], quiet = TRUE)

disc_table <- bind_rows(disc_results)
write.csv(disc_table, paste0(TAB_DIR, "calon2_discrimination.csv"), row.names = FALSE)
cat("  Saved: calon2_discrimination.csv\n")

# =============================================================================
# SECTION 7: CALIBRATION
# =============================================================================

cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("SECTION 7: CALIBRATION\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

cal_results <- list()

cal_results[[1]] <- calibration_metrics(ukb$ascvd_combined, ukb$prob_calon2, "CALON-2 in UKB")
cal_results[[2]] <- calibration_metrics(ukb$ascvd_combined, ukb$prob_safeheart, "SAFEHEART-RE in UKB")
cal_results[[3]] <- calibration_metrics(d3$ascvd_combined, d3$prob_calon2, "CALON-2 in DRAGON3")
cal_results[[4]] <- calibration_metrics(d3$ascvd_combined, d3$prob_safeheart, "SAFEHEART-RE in DRAGON3")

if (wales_apob_valid && !is.null(wl)) {
  cal_results[[5]] <- calibration_metrics(wl$ascvd_combined, wl$prob_calon2, "CALON-2 in Wales")
  cal_results[[6]] <- calibration_metrics(wl$ascvd_combined, wl$prob_safeheart, "SAFEHEART-RE in Wales")
}

cal_table <- bind_rows(cal_results)
write.csv(cal_table, paste0(TAB_DIR, "calon2_calibration.csv"), row.names = FALSE)
cat("\n  Saved: calon2_calibration.csv\n")

# =============================================================================
# SECTION 8: RECLASSIFICATION (NRI / IDI)
# =============================================================================

cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("SECTION 8: RECLASSIFICATION (NRI / IDI)\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

risk_cuts <- c(0.10, 0.20)

compute_nri <- function(outcome, p_old, p_new, cohort_name) {
  valid <- !is.na(outcome) & !is.na(p_old) & !is.na(p_new)
  y    <- outcome[valid]
  pold <- p_old[valid]
  pnew <- p_new[valid]

  if (sum(valid) < 30 || sum(y) < 5) {
    cat(sprintf("  %s: INSUFFICIENT DATA\n", cohort_name))
    return(NULL)
  }

  tryCatch({
    nri_result <- nricens::nribin(
      event = y, p.std = pold, p.new = pnew,
      cut = risk_cuts, niter = 1000, msg = FALSE
    )

    cat(sprintf("  %s:\n", cohort_name))
    cat(sprintf("    Categorical NRI:  %.4f (%.4f\u2013%.4f), p=%.4f\n",
                nri_result$nri["Categorical NRI", "Estimate"],
                nri_result$nri["Categorical NRI", "Lower"],
                nri_result$nri["Categorical NRI", "Upper"],
                nri_result$nri["Categorical NRI", "P-value"]))
    cat(sprintf("    Continuous NRI:   %.4f (%.4f\u2013%.4f), p=%.4f\n",
                nri_result$nri["Continuous NRI", "Estimate"],
                nri_result$nri["Continuous NRI", "Lower"],
                nri_result$nri["Continuous NRI", "Upper"],
                nri_result$nri["Continuous NRI", "P-value"]))
    cat(sprintf("    IDI:              %.4f (%.4f\u2013%.4f), p=%.4f\n\n",
                nri_result$nri["IDI", "Estimate"],
                nri_result$nri["IDI", "Lower"],
                nri_result$nri["IDI", "Upper"],
                nri_result$nri["IDI", "P-value"]))

    return(data.frame(
      Cohort = cohort_name,
      Cat_NRI = nri_result$nri["Categorical NRI", "Estimate"],
      Cat_NRI_Lower = nri_result$nri["Categorical NRI", "Lower"],
      Cat_NRI_Upper = nri_result$nri["Categorical NRI", "Upper"],
      Cat_NRI_p = nri_result$nri["Categorical NRI", "P-value"],
      Cont_NRI = nri_result$nri["Continuous NRI", "Estimate"],
      Cont_NRI_Lower = nri_result$nri["Continuous NRI", "Lower"],
      Cont_NRI_Upper = nri_result$nri["Continuous NRI", "Upper"],
      Cont_NRI_p = nri_result$nri["Continuous NRI", "P-value"],
      IDI = nri_result$nri["IDI", "Estimate"],
      IDI_Lower = nri_result$nri["IDI", "Lower"],
      IDI_Upper = nri_result$nri["IDI", "Upper"],
      IDI_p = nri_result$nri["IDI", "P-value"],
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    cat(sprintf("  %s: NRI computation failed: %s\n", cohort_name, e$message))
    return(NULL)
  })
}

nri_results <- list()
nri_results[[1]] <- compute_nri(ukb$ascvd_combined, ukb$prob_safeheart, ukb$prob_calon2, "UKB")
nri_results[[2]] <- compute_nri(d3$ascvd_combined, d3$prob_safeheart, d3$prob_calon2, "DRAGON3")
if (wales_apob_valid && !is.null(wl)) {
  nri_results[[3]] <- compute_nri(wl$ascvd_combined, wl$prob_safeheart, wl$prob_calon2, "Wales ApoB")
}

nri_table <- bind_rows(nri_results)
write.csv(nri_table, paste0(TAB_DIR, "calon2_reclassification.csv"), row.names = FALSE)
cat("  Saved: calon2_reclassification.csv\n")

# =============================================================================
# SECTION 9: SUBGROUP ANALYSIS
# =============================================================================

cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("SECTION 9: SUBGROUP ANALYSIS\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

subgroup_auc <- function(data, subgroup_col, subgroup_label, prob_col = "prob_calon2",
                          outcome_col = "ascvd_combined") {
  results <- list()
  levels <- unique(data[[subgroup_col]])
  levels <- levels[!is.na(levels)]

  for (lev in levels) {
    sub <- data[data[[subgroup_col]] == lev & !is.na(data[[prob_col]]) &
                  !is.na(data[[outcome_col]]), ]
    n <- nrow(sub)
    events <- sum(sub[[outcome_col]])

    if (n >= 20 && events >= 3) {
      roc_obj <- tryCatch(roc(sub[[outcome_col]], sub[[prob_col]], quiet = TRUE),
                           error = function(e) NULL)
      if (!is.null(roc_obj)) {
        ci <- ci.auc(roc_obj)
        results[[length(results) + 1]] <- data.frame(
          Subgroup = subgroup_label, Level = as.character(lev),
          N = n, Events = events,
          AUC = as.numeric(auc(roc_obj)),
          AUC_Lower = ci[1], AUC_Upper = ci[3],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  return(bind_rows(results))
}

# DRAGON3 subgroups
sub_results <- list()

# By gene
sub_results[[1]] <- subgroup_auc(d3, "gene", "Gene (DRAGON3)")

# By sex
d3$sex_label <- ifelse(d3$sex == 1, "Male", "Female")
sub_results[[2]] <- subgroup_auc(d3, "sex_label", "Sex (DRAGON3)")

# By age group
d3$age_group <- ifelse(d3$age < 50, "<50", ">=50")
sub_results[[3]] <- subgroup_auc(d3, "age_group", "Age (DRAGON3)")

# By statin status
if (!all(is.na(d3$on_statin))) {
  d3$statin_label <- ifelse(d3$on_statin == 1, "On Statin", "No Statin")
  sub_results[[4]] <- subgroup_auc(d3, "statin_label", "Statin (DRAGON3)")
}

# UKB subgroups
ukb$sex_label <- ifelse(ukb$sex == 1, "Male", "Female")
sub_results[[5]] <- subgroup_auc(ukb, "sex_label", "Sex (UKB)")

ukb$age_group <- ifelse(ukb$age < 50, "<50", ">=50")
sub_results[[6]] <- subgroup_auc(ukb, "age_group", "Age (UKB)")

subgroup_table <- bind_rows(sub_results)
cat("  Subgroup AUCs:\n")
if (nrow(subgroup_table) > 0) {
  for (i in 1:nrow(subgroup_table)) {
    cat(sprintf("    %-25s %-10s: AUC=%.3f (%.3f\u2013%.3f), n=%d, events=%d\n",
                subgroup_table$Subgroup[i], subgroup_table$Level[i],
                subgroup_table$AUC[i], subgroup_table$AUC_Lower[i], subgroup_table$AUC_Upper[i],
                subgroup_table$N[i], subgroup_table$Events[i]))
  }
}

write.csv(subgroup_table, paste0(TAB_DIR, "calon2_subgroups.csv"), row.names = FALSE)
cat("\n  Saved: calon2_subgroups.csv\n")

# =============================================================================
# SECTION 10: DECISION CURVE ANALYSIS
# =============================================================================

cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("SECTION 10: DECISION CURVE ANALYSIS\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

# DCA on DRAGON3
d3_dca <- d3[!is.na(d3$ascvd_combined) & !is.na(d3$prob_calon2) & !is.na(d3$prob_safeheart), ]
d3_dca$ascvd_combined <- as.integer(d3_dca$ascvd_combined)

dca_dragon3 <- tryCatch({
  dcurves::dca(
    ascvd_combined ~ prob_calon2 + prob_safeheart,
    data = d3_dca,
    thresholds = seq(0.01, 0.60, by = 0.01),
    label = list(prob_calon2 = "CALON-2", prob_safeheart = "SAFEHEART-RE")
  )
}, error = function(e) {
  cat(sprintf("  DCA computation failed: %s\n", e$message))
  NULL
})

if (!is.null(dca_dragon3)) {
  cat("  DCA computed successfully for DRAGON3.\n")
}

# #############################################################################
#                                                                             #
#  FIGURES                                                                    #
#                                                                             #
# #############################################################################

cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("GENERATING FIGURES\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

# #############################################################################
# FIGURE 1: DUAL ROC — A) UKB internal, B) DRAGON3 external
# #############################################################################

cat("  Figure 1: Dual ROC curves...\n")

roc_data_ukb <- data.frame(
  sens = rev(roc_ukb_calon2$sensitivities),
  spec = rev(1 - roc_ukb_calon2$specificities),
  model = "CALON-2"
) %>% bind_rows(data.frame(
  sens = rev(roc_ukb_sh$sensitivities),
  spec = rev(1 - roc_ukb_sh$specificities),
  model = "SAFEHEART-RE"
))

roc_data_d3 <- data.frame(
  sens = rev(roc_d3_calon2$sensitivities),
  spec = rev(1 - roc_d3_calon2$specificities),
  model = "CALON-2"
) %>% bind_rows(data.frame(
  sens = rev(roc_d3_sh$sensitivities),
  spec = rev(1 - roc_d3_sh$specificities),
  model = "SAFEHEART-RE"
))

auc_ukb_c2 <- sprintf("%.3f", as.numeric(auc(roc_ukb_calon2)))
auc_ukb_sh <- sprintf("%.3f", as.numeric(auc(roc_ukb_sh)))
auc_d3_c2  <- sprintf("%.3f", as.numeric(auc(roc_d3_calon2)))
auc_d3_sh  <- sprintf("%.3f", as.numeric(auc(roc_d3_sh)))

p1a <- ggplot(roc_data_ukb, aes(x = spec, y = sens, color = model)) +
  geom_line(linewidth = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = c("CALON-2" = COL$calon2, "SAFEHEART-RE" = COL$safeheart),
                     labels = c(paste0("CALON-2 (AUC ", auc_ukb_c2, ")"),
                                paste0("SAFEHEART-RE (AUC ", auc_ukb_sh, ")"))) +
  labs(title = "A. UKB (Development)", x = "1 - Specificity", y = "Sensitivity", color = NULL) +
  coord_equal() + theme_nature() +
  theme(legend.position = c(0.65, 0.2))

p1b <- ggplot(roc_data_d3, aes(x = spec, y = sens, color = model)) +
  geom_line(linewidth = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = c("CALON-2" = COL$calon2, "SAFEHEART-RE" = COL$safeheart),
                     labels = c(paste0("CALON-2 (AUC ", auc_d3_c2, ")"),
                                paste0("SAFEHEART-RE (AUC ", auc_d3_sh, ")"))) +
  labs(title = "B. DRAGON3 (External Validation)", x = "1 - Specificity", y = "Sensitivity", color = NULL) +
  coord_equal() + theme_nature() +
  theme(legend.position = c(0.65, 0.2))

fig1 <- p1a + p1b
save_figure(fig1, "calon2_fig1_dual_roc", width = WIDTH_DOUBLE, height = WIDTH_SINGLE)

# #############################################################################
# FIGURE 2: CALIBRATION PLOTS — A) UKB, B) DRAGON3
# #############################################################################

cat("  Figure 2: Calibration plots...\n")

calibration_plot_data <- function(observed, predicted, n_groups = 10) {
  pred <- pmax(pmin(predicted, 1 - 1e-10), 1e-10)
  valid <- !is.na(observed) & !is.na(pred)
  obs <- observed[valid]; pred <- pred[valid]
  groups <- cut(pred, breaks = quantile(pred, probs = seq(0, 1, length.out = n_groups + 1)),
                include.lowest = TRUE)
  data.frame(
    predicted_mean = tapply(pred, groups, mean),
    observed_mean  = tapply(obs, groups, mean),
    n = tapply(obs, groups, length),
    stringsAsFactors = FALSE
  )
}

cal_ukb_c2 <- calibration_plot_data(ukb$ascvd_combined, ukb$prob_calon2)
cal_ukb_c2$model <- "CALON-2"
cal_ukb_sh <- calibration_plot_data(ukb$ascvd_combined, ukb$prob_safeheart)
cal_ukb_sh$model <- "SAFEHEART-RE"

cal_d3_c2 <- calibration_plot_data(d3$ascvd_combined, d3$prob_calon2)
cal_d3_c2$model <- "CALON-2"
cal_d3_sh <- calibration_plot_data(d3$ascvd_combined, d3$prob_safeheart)
cal_d3_sh$model <- "SAFEHEART-RE"

p2a <- ggplot(bind_rows(cal_ukb_c2, cal_ukb_sh),
              aes(x = predicted_mean, y = observed_mean, color = model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("CALON-2" = COL$calon2, "SAFEHEART-RE" = COL$safeheart)) +
  labs(title = "A. UKB (Development)",
       x = "Predicted probability", y = "Observed proportion", color = NULL) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_nature() + theme(legend.position = c(0.25, 0.85))

p2b <- ggplot(bind_rows(cal_d3_c2, cal_d3_sh),
              aes(x = predicted_mean, y = observed_mean, color = model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("CALON-2" = COL$calon2, "SAFEHEART-RE" = COL$safeheart)) +
  labs(title = "B. DRAGON3 (External Validation)",
       x = "Predicted probability", y = "Observed proportion", color = NULL) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_nature() + theme(legend.position = c(0.25, 0.85))

fig2 <- p2a + p2b
save_figure(fig2, "calon2_fig2_calibration", width = WIDTH_DOUBLE, height = WIDTH_SINGLE)

# #############################################################################
# FIGURE 3: FOREST PLOT — Subgroup AUCs
# #############################################################################

cat("  Figure 3: Forest plot of subgroup AUCs...\n")

if (nrow(subgroup_table) > 0) {
  subgroup_table$label <- paste0(subgroup_table$Subgroup, ": ", subgroup_table$Level)
  subgroup_table$label <- factor(subgroup_table$label,
                                  levels = rev(subgroup_table$label))

  fig3 <- ggplot(subgroup_table, aes(x = AUC, y = label)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey70") +
    geom_errorbarh(aes(xmin = AUC_Lower, xmax = AUC_Upper),
                   height = 0.2, linewidth = 0.4, color = COL$calon2) +
    geom_point(size = 2.5, color = COL$calon2) +
    geom_text(aes(label = sprintf("%.3f", AUC)), hjust = -0.3, size = 2.2) +
    labs(title = "CALON-2 Discrimination by Subgroup",
         x = "AUC (95% CI)", y = NULL) +
    coord_cartesian(xlim = c(0.4, 1.0)) +
    theme_nature() + theme(legend.position = "none")

  save_figure(fig3, "calon2_fig3_forest_subgroup", width = WIDTH_DOUBLE, height = 5)
} else {
  cat("  Skipping: no subgroup data\n")
}

# #############################################################################
# FIGURE 4: DECISION CURVE ANALYSIS
# #############################################################################

cat("  Figure 4: Decision curve...\n")

if (!is.null(dca_dragon3)) {
  dca_df <- as.data.frame(dca_dragon3)
  # Rename for plotting
  dca_df$label[dca_df$label == "prob_calon2"] <- "CALON-2"
  dca_df$label[dca_df$label == "prob_safeheart"] <- "SAFEHEART-RE"
  dca_df$label[dca_df$label == "Treat All"] <- "Treat All"
  dca_df$label[dca_df$label == "Treat None"] <- "Treat None"

  fig4 <- ggplot(dca_df, aes(x = threshold, y = net_benefit, color = label, linetype = label)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = c("CALON-2" = COL$calon2, "SAFEHEART-RE" = COL$safeheart,
                                   "Treat All" = COL$grey_dark, "Treat None" = COL$grey_light)) +
    scale_linetype_manual(values = c("CALON-2" = "solid", "SAFEHEART-RE" = "solid",
                                      "Treat All" = "dashed", "Treat None" = "dotted")) +
    labs(title = "Decision Curve Analysis (DRAGON3)",
         x = "Threshold Probability", y = "Net Benefit", color = NULL, linetype = NULL) +
    coord_cartesian(ylim = c(-0.05, max(dca_df$net_benefit, na.rm = TRUE) * 1.1)) +
    theme_nature()

  save_figure(fig4, "calon2_fig4_decision_curve", width = WIDTH_DOUBLE, height = 4)
} else {
  cat("  Skipping: DCA not available\n")
}

# #############################################################################
# FIGURE 5: VARIABLE IMPORTANCE (Elastic Net Coefficients)
# #############################################################################

cat("  Figure 5: Variable importance...\n")

coef_plot_df <- calon2_coef_df[calon2_coef_df$variable != "(Intercept)" &
                                 calon2_coef_df$coefficient != 0, ]

if (nrow(coef_plot_df) > 0) {
  coef_plot_df$abs_coef <- abs(coef_plot_df$coefficient)
  coef_plot_df$direction <- ifelse(coef_plot_df$coefficient > 0, "Risk Factor", "Protective")
  coef_plot_df$variable <- factor(coef_plot_df$variable,
                                   levels = coef_plot_df$variable[order(coef_plot_df$abs_coef)])

  fig5 <- ggplot(coef_plot_df, aes(x = variable, y = coefficient, fill = direction)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = c("Risk Factor" = COL$event, "Protective" = COL$noevent)) +
    labs(title = "CALON-2 Core Model Coefficients",
         subtitle = "Elastic net selected predictors",
         x = NULL, y = "Coefficient", fill = NULL) +
    coord_flip() +
    theme_nature() + theme(legend.position = "bottom")

  save_figure(fig5, "calon2_fig5_variable_importance", width = WIDTH_DOUBLE, height = 4)
}

# #############################################################################
# FIGURE 6: NRI WATERFALL — Reclassification vs SAFEHEART-RE
# #############################################################################

cat("  Figure 6: NRI reclassification...\n")

# Build reclassification table for DRAGON3
d3_reclass <- d3[!is.na(d3$ascvd_combined) & !is.na(d3$prob_calon2) & !is.na(d3$prob_safeheart), ]

cat_old <- cut(d3_reclass$prob_safeheart, breaks = c(0, 0.10, 0.20, 1),
               include.lowest = TRUE, labels = c("Low (<10%)", "Medium (10-20%)", "High (>20%)"))
cat_new <- cut(d3_reclass$prob_calon2, breaks = c(0, 0.10, 0.20, 1),
               include.lowest = TRUE, labels = c("Low (<10%)", "Medium (10-20%)", "High (>20%)"))

d3_reclass$move <- as.character(cat_new) != as.character(cat_old)
d3_reclass$direction <- ifelse(as.numeric(cat_new) > as.numeric(cat_old), "Reclassified Up",
                        ifelse(as.numeric(cat_new) < as.numeric(cat_old), "Reclassified Down", "Same"))

reclass_summary <- d3_reclass %>%
  group_by(ascvd_combined, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(event_label = ifelse(ascvd_combined == 1, "Events", "Non-events"))

if (nrow(reclass_summary) > 0) {
  fig6 <- ggplot(reclass_summary, aes(x = event_label, y = n, fill = direction)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_text(aes(label = n), position = position_dodge(width = 0.7), vjust = -0.5, size = 2.2) +
    scale_fill_manual(values = c("Reclassified Up" = COL$event,
                                  "Reclassified Down" = COL$noevent,
                                  "Same" = COL$grey_light)) +
    labs(title = "Reclassification: CALON-2 vs SAFEHEART-RE (DRAGON3)",
         subtitle = "Risk categories: Low (<10%), Medium (10-20%), High (>20%)",
         x = NULL, y = "Number of patients", fill = NULL) +
    theme_nature()

  save_figure(fig6, "calon2_fig6_reclassification", width = WIDTH_DOUBLE, height = 4)
}

# #############################################################################
# FIGURE 7: GENE-SPECIFIC ROC (DRAGON3)
# #############################################################################

cat("  Figure 7: Gene-specific ROC...\n")

gene_rocs <- list()
gene_labels <- c()

for (g in c("LDLR", "APOB", "PCSK9")) {
  sub <- d3[d3$gene == g & !is.na(d3$prob_calon2) & !is.na(d3$ascvd_combined), ]
  if (nrow(sub) >= 20 && sum(sub$ascvd_combined) >= 3) {
    roc_g <- roc(sub$ascvd_combined, sub$prob_calon2, quiet = TRUE)
    gene_rocs[[g]] <- data.frame(
      sens = rev(roc_g$sensitivities),
      spec = rev(1 - roc_g$specificities),
      gene = g
    )
    gene_labels <- c(gene_labels, paste0(g, " (AUC ", sprintf("%.3f", as.numeric(auc(roc_g))),
                                          ", n=", nrow(sub), ")"))
  }
}

if (length(gene_rocs) > 0) {
  gene_roc_df <- bind_rows(gene_rocs)

  fig7 <- ggplot(gene_roc_df, aes(x = spec, y = sens, color = gene)) +
    geom_line(linewidth = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    scale_color_manual(values = c("LDLR" = COL$ldlr, "APOB" = COL$apob, "PCSK9" = COL$pcsk9),
                       labels = gene_labels) +
    labs(title = "CALON-2 Discrimination by Gene (DRAGON3)",
         x = "1 - Specificity", y = "Sensitivity", color = NULL) +
    coord_equal() + theme_nature() +
    theme(legend.position = c(0.65, 0.2))

  save_figure(fig7, "calon2_fig7_gene_roc", width = WIDTH_SINGLE * 1.5, height = WIDTH_SINGLE * 1.5)
} else {
  cat("  Skipping: insufficient gene-specific data\n")
}

# #############################################################################
# FIGURE 8: RISK DISTRIBUTION — Predicted probability histograms
# #############################################################################

cat("  Figure 8: Risk distribution...\n")

risk_df_d3 <- d3[!is.na(d3$prob_calon2) & !is.na(d3$ascvd_combined),
                  c("prob_calon2", "ascvd_combined")]
risk_df_d3$event_label <- ifelse(risk_df_d3$ascvd_combined == 1, "ASCVD Event", "No Event")

fig8 <- ggplot(risk_df_d3, aes(x = prob_calon2, fill = event_label)) +
  geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = c(0.10, 0.20), linetype = "dashed", color = "grey40", linewidth = 0.3) +
  annotate("text", x = 0.05, y = Inf, label = "Low", vjust = 1.5, size = 2, color = "grey40") +
  annotate("text", x = 0.15, y = Inf, label = "Med", vjust = 1.5, size = 2, color = "grey40") +
  annotate("text", x = 0.35, y = Inf, label = "High", vjust = 1.5, size = 2, color = "grey40") +
  scale_fill_manual(values = c("ASCVD Event" = COL$event, "No Event" = COL$noevent)) +
  labs(title = "CALON-2 Predicted Risk Distribution (DRAGON3)",
       x = "Predicted ASCVD Probability", y = "Count", fill = NULL) +
  theme_nature()

save_figure(fig8, "calon2_fig8_risk_distribution", width = WIDTH_DOUBLE, height = 3.5)

# #############################################################################
#                                                                             #
#  TABLES                                                                     #
#                                                                             #
# #############################################################################

cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("GENERATING TABLES\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

# #############################################################################
# TABLE 1: BASELINE CHARACTERISTICS — UKB vs DRAGON3 vs Wales
# #############################################################################

cat("  Table 1: Baseline characteristics...\n")

# Prepare UKB subset
ukb_t1 <- ukb %>%
  select(any_of(c("age", "sex", "apob", "re_ldl", "hdl", "trig",
                    "inv_apoa", "apob_ldl_ratio", "lpa_binary",
                    "smoking_binary", "diabetes", "hypertension", "bmi",
                    "ascvd_combined"))) %>%
  mutate(cohort = "UKB (N=1,623)")

# Prepare DRAGON3 subset
d3_t1 <- d3 %>%
  select(any_of(c("age", "sex", "apob", "re_ldl", "hdl", "trig",
                    "inv_apoa", "apob_ldl_ratio", "lpa_binary",
                    "smoking_binary", "diabetes", "hypertension", "bmi",
                    "ascvd_combined"))) %>%
  mutate(cohort = sprintf("DRAGON3 (N=%d)", nrow(d3)))

# Combine
combined_t1 <- bind_rows(ukb_t1, d3_t1)

# Add Wales if valid
if (wales_apob_valid && !is.null(wl)) {
  wl_t1 <- wl %>%
    select(any_of(c("age", "sex", "apob", "re_ldl", "hdl", "trig",
                      "inv_apoa", "apob_ldl_ratio", "lpa_binary",
                      "smoking_binary", "diabetes", "hypertension", "bmi",
                      "ascvd_combined"))) %>%
    mutate(cohort = sprintf("Wales ApoB (N=%d)", nrow(wl)))
  combined_t1 <- bind_rows(combined_t1, wl_t1)
}

# Define continuous and categorical variables
cont_vars <- intersect(c("age", "apob", "re_ldl", "hdl", "trig", "bmi"),
                        names(combined_t1))
cat_vars  <- intersect(c("sex", "lpa_binary", "smoking_binary", "diabetes",
                           "hypertension", "ascvd_combined"),
                        names(combined_t1))

table1 <- tryCatch({
  CreateTableOne(vars = c(cont_vars, cat_vars),
                 strata = "cohort",
                 data = combined_t1,
                 factorVars = cat_vars,
                 test = TRUE)
}, error = function(e) {
  cat(sprintf("    TableOne error: %s\n", e$message))
  NULL
})

if (!is.null(table1)) {
  table1_print <- print(table1, smd = TRUE, printToggle = FALSE)
  write.csv(table1_print, paste0(TAB_DIR, "calon2_table1_baseline.csv"))
  cat("  Saved: calon2_table1_baseline.csv\n")
}

# #############################################################################
# TABLE 2: CALON-2 CORE COEFFICIENTS with OR (from development script)
# #############################################################################

cat("  Table 2: Model coefficients...\n")

# Load OR table if available from 05 script
or_file <- paste0(TAB_DIR, "calon2_core_logistic_ORs.csv")
if (file.exists(or_file)) {
  or_table <- read.csv(or_file, stringsAsFactors = FALSE)
  cat("  Loaded logistic OR table from development script.\n")
} else {
  # Create from elastic net coefficients
  or_table <- calon2_coef_df
  or_table$OR <- exp(or_table$coefficient)
  cat("  Created approximate OR table from elastic net coefficients.\n")
}

write.csv(or_table, paste0(TAB_DIR, "calon2_table2_coefficients.csv"), row.names = FALSE)
cat("  Saved: calon2_table2_coefficients.csv\n")

# #############################################################################
# TABLE 3: DISCRIMINATION — already saved as calon2_discrimination.csv
# #############################################################################

cat("  Table 3: Discrimination (already saved as calon2_discrimination.csv)\n")

# #############################################################################
# TABLE 4: CALIBRATION — already saved as calon2_calibration.csv
# #############################################################################

cat("  Table 4: Calibration (already saved as calon2_calibration.csv)\n")

# #############################################################################
# TABLE 5: RECLASSIFICATION — already saved as calon2_reclassification.csv
# #############################################################################

cat("  Table 5: Reclassification (already saved as calon2_reclassification.csv)\n")

# #############################################################################
# TABLE 6: SUBGROUP ANALYSIS — already saved as calon2_subgroups.csv
# #############################################################################

cat("  Table 6: Subgroups (already saved as calon2_subgroups.csv)\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat("\u2554\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2557\n")
cat("\u2551  CALON-2 EXTERNAL VALIDATION COMPLETE                            \u2551\n")
cat("\u255a\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u255d\n\n")

cat("  DISCRIMINATION (AUC):\n")
if (nrow(disc_table) > 0) {
  for (i in 1:nrow(disc_table)) {
    cat(sprintf("    %-25s: CALON-2=%.3f, SAFEHEART-RE=%.3f, \u0394=%+.3f, p=%s\n",
                disc_table$Cohort[i], disc_table$CALON2_AUC[i], disc_table$SAFEHEART_AUC[i],
                disc_table$Delta_AUC[i],
                ifelse(disc_table$DeLong_p[i] < 0.001, "<0.001",
                       sprintf("%.4f", disc_table$DeLong_p[i]))))
  }
}

cat("\n  FILES GENERATED:\n")
cat("    Figures (8): calon2_fig1-fig8 (.pdf + .png)\n")
cat("    Tables  (6): calon2_table1-table6 (.csv)\n")
cat("    Plus: calon2_discrimination.csv, calon2_calibration.csv,\n")
cat("          calon2_reclassification.csv, calon2_subgroups.csv\n\n")

cat("  Next steps:\n")
cat("    1. Review figures for Nature submission\n")
cat("    2. Check calibration metrics (slope near 1, intercept near 0)\n")
cat("    3. Verify DeLong p-values for CALON-2 vs SAFEHEART-RE\n")
cat("    4. If CALON-2 wins: write manuscript!\n")
cat("    5. If calibration off: consider recalibration in DRAGON3\n\n")

cat("  CALON-2 External Validation Pipeline Complete.\n")
cat("================================================================\n\n")

################################################################################
#                                                                              #
#  CALON-2: PORTABLE MODEL — DEVELOP ON UKB, VALIDATE ON ALL WALES           #
#  "Only variables available in BOTH datasets"                                #
#                                                                              #
#  Strategy:                                                                   #
#    1. Map Wales variables → UKB variable names                              #
#    2. Build Portable model on UKB (only Wales-available vars)               #
#    3. Freeze coefficients                                                    #
#    4. Apply frozen model to Wales → External AUC + Calibration              #
#    5. Compare vs SAFEHEART-RE                                                #
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

## NOTE: Only load packages actually used. rms/caret/randomForest were removed
## because their massive dependency trees cause namespace conflicts on Posit Cloud
## (e.g., Hmisc from rms can mask base functions, causing cryptic errors).
required_packages <- c(
    "glmnet", "pROC", "xgboost", "car"
)

for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
        library(pkg, character.only = TRUE)
    }
}

cat("\n")
cat("======================================================================\n")
cat("  CALON-2: PORTABLE MODEL — Develop UKB → Validate Wales             \n")
cat("  Only variables available in BOTH datasets                           \n")
cat("======================================================================\n\n")

# =============================================================================
# CONFIGURATION — auto-detect environment
# CHECK CLOUD FIRST (Posit Cloud can false-positive on Windows paths)
# =============================================================================

CLOUD_PATH  <- "/cloud/project/"
LOCAL_PATH  <- "C:/Users/nader/Downloads/calon_ukb_pipeline/"

# Check Cloud FIRST, then Local, then fallback to working directory
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

cat(sprintf("  BASE_DIR = %s\n", BASE_DIR))
cat(sprintf("  Working directory = %s\n", getwd()))

OUT_DIR <- paste0(BASE_DIR, "output/")
TAB_DIR <- paste0(OUT_DIR, "tables/")
FIG_DIR <- paste0(OUT_DIR, "figures/")
for (d in c(TAB_DIR, FIG_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# =============================================================================
# SECTION 1: LOAD UKB DEVELOPMENT DATA
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 1: LOAD UKB DATA\n")
cat("===================================================================\n\n")

# Search for the merged CSV in multiple locations
merged_file <- NULL
search_paths <- c(
    paste0(OUT_DIR, "calon2_full_merged.csv"),
    paste0(BASE_DIR, "calon2_full_merged.csv"),
    paste0(getwd(), "/calon2_full_merged.csv"),
    paste0(getwd(), "/output/calon2_full_merged.csv"),
    "/cloud/project/calon2_full_merged.csv",
    "/cloud/project/output/calon2_full_merged.csv"
)

for (sp in search_paths) {
    if (file.exists(sp)) { merged_file <- sp; break }
}

if (is.null(merged_file)) {
    cat("  SEARCHED IN:\n")
    for (sp in search_paths) cat(sprintf("    %s → %s\n", sp, file.exists(sp)))
    stop("ERROR: calon2_full_merged.csv not found in any expected location")
}
cat(sprintf("  Found: %s\n", merged_file))

ukb <- read.csv(merged_file, stringsAsFactors = FALSE)
cat(sprintf("  UKB: %d patients x %d columns\n", nrow(ukb), ncol(ukb)))
cat(sprintf("  ASCVD events: %d (%.1f%%)\n", sum(ukb$ascvd_combined), 100*mean(ukb$ascvd_combined)))

# =============================================================================
# SECTION 2: LOAD WALES DATA
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 2: LOAD WALES DATA\n")
cat("===================================================================\n\n")

# Search for Wales CSV in multiple locations and name variants
wales_file <- NULL
search_dirs <- unique(c(BASE_DIR, getwd(), paste0(getwd(), "/"),
                         "/cloud/project/", "/cloud/project/output/"))

for (sd in search_dirs) {
    if (!dir.exists(sd)) next
    found <- list.files(sd, pattern = "WALES.*FH.*CLEAN", full.names = TRUE,
                        ignore.case = TRUE)
    if (length(found) > 0) { wales_file <- found[1]; break }
}

# Also try exact names in each directory
if (is.null(wales_file)) {
    wales_names <- c("WALES_FH_CLEANED (1) - Copy.csv",
                     "WALES_FH_CLEANED (1).csv",
                     "WALES_FH_CLEANED.csv",
                     "wales_fh_cleaned.csv")
    for (sd in search_dirs) {
        for (wn in wales_names) {
            wf <- paste0(sd, wn)
            if (file.exists(wf)) { wales_file <- wf; break }
        }
        if (!is.null(wales_file)) break
    }
}

if (is.null(wales_file)) {
    cat("  Wales file NOT FOUND. Searched in:\n")
    for (sd in search_dirs) {
        cat(sprintf("    %s → files: %s\n", sd,
                    paste(list.files(sd, pattern = "WALES|wales", ignore.case = TRUE),
                          collapse = ", ")))
    }
    stop("ERROR: Wales FH cleaned CSV not found. Upload it to Posit Cloud project root.")
}

wales_raw <- read.csv(wales_file, stringsAsFactors = FALSE)
cat(sprintf("  Wales raw: %d patients x %d columns\n", nrow(wales_raw), ncol(wales_raw)))
cat(sprintf("  File: %s\n", basename(wales_file)))

# =============================================================================
# SECTION 3: HARMONISE WALES → UKB VARIABLE NAMES
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 3: HARMONISE WALES VARIABLES\n")
cat("===================================================================\n\n")

wales <- data.frame(row.names = 1:nrow(wales_raw))

# -- 3A: Age (from DOB) --
if ("DOB" %in% names(wales_raw)) {
    dob <- as.Date(wales_raw$DOB, format = "%Y-%m-%d")
    if (sum(!is.na(dob)) < 100) dob <- as.Date(wales_raw$DOB, format = "%d/%m/%Y")
    # Censor date: use ascvd event age or last measurement date or today
    wales$age <- as.numeric(difftime(Sys.Date(), dob, units = "days")) / 365.25
    # If BMI_AGE or corneal_age gives better context, use those
    if ("BMI_AGE" %in% names(wales_raw)) {
        has_bmi_age <- !is.na(wales_raw$BMI_AGE) & wales_raw$BMI_AGE > 10
        # Use BMI_AGE as a reasonable clinical age where available
    }
    cat(sprintf("  age: %d non-missing (median %.0f)\n",
                sum(!is.na(wales$age)), median(wales$age, na.rm = TRUE)))
} else {
    wales$age <- NA
    cat("  WARNING: DOB not found\n")
}

# -- 3B: Sex --
wales$sex <- ifelse(wales_raw$Gender == "M", 1,
             ifelse(wales_raw$Gender == "F", 0, NA))
cat(sprintf("  sex: %d non-missing (%d male, %d female)\n",
            sum(!is.na(wales$sex)), sum(wales$sex == 1, na.rm = TRUE),
            sum(wales$sex == 0, na.rm = TRUE)))

# -- 3C: Lipids (use first available measurement) --
wales$re_ldl <- as.numeric(wales_raw$LDL.1)
wales$hdl    <- as.numeric(wales_raw$HDL.1)
wales$tc     <- as.numeric(wales_raw$TC.1)
wales$trig   <- as.numeric(wales_raw$TRG.1)

# Fill from later measurements if .1 is missing
for (suffix in c(".2", ".3", ".4")) {
    ldl_col <- paste0("LDL", suffix)
    hdl_col <- paste0("HDL", suffix)
    tc_col  <- paste0("TC", suffix)
    tg_col  <- paste0("TRG", suffix)

    if (ldl_col %in% names(wales_raw)) {
        fill_idx <- is.na(wales$re_ldl) & !is.na(as.numeric(wales_raw[[ldl_col]]))
        wales$re_ldl[fill_idx] <- as.numeric(wales_raw[[ldl_col]][fill_idx])
    }
    if (hdl_col %in% names(wales_raw)) {
        fill_idx <- is.na(wales$hdl) & !is.na(as.numeric(wales_raw[[hdl_col]]))
        wales$hdl[fill_idx] <- as.numeric(wales_raw[[hdl_col]][fill_idx])
    }
    if (tc_col %in% names(wales_raw)) {
        fill_idx <- is.na(wales$tc) & !is.na(as.numeric(wales_raw[[tc_col]]))
        wales$tc[fill_idx] <- as.numeric(wales_raw[[tc_col]][fill_idx])
    }
    if (tg_col %in% names(wales_raw)) {
        fill_idx <- is.na(wales$trig) & !is.na(as.numeric(wales_raw[[tg_col]]))
        wales$trig[fill_idx] <- as.numeric(wales_raw[[tg_col]][fill_idx])
    }
}

cat(sprintf("  re_ldl: %d non-missing\n", sum(!is.na(wales$re_ldl))))
cat(sprintf("  hdl:    %d non-missing\n", sum(!is.na(wales$hdl))))
cat(sprintf("  tc:     %d non-missing\n", sum(!is.na(wales$tc))))
cat(sprintf("  trig:   %d non-missing\n", sum(!is.na(wales$trig))))

# -- 3D: Smoking --
# Wales "Smoking" may be text or numeric
smk <- wales_raw$Smoking
wales$smoking_binary <- ifelse(smk %in% c(1, "1", "Yes", "yes", "Current", "current",
                                           "Ex", "ex", "Former", "former", "Ex-smoker",
                                           "Current smoker", "Previous"), 1,
                        ifelse(smk %in% c(0, "0", "No", "no", "Never", "never",
                                           "Non-smoker"), 0, NA))
cat(sprintf("  smoking_binary: %d non-missing (%d ever-smokers)\n",
            sum(!is.na(wales$smoking_binary)), sum(wales$smoking_binary == 1, na.rm = TRUE)))

# -- 3E: Diabetes --
dm <- wales_raw$Diabetes
wales$diabetes <- ifelse(dm %in% c(1, "1", "Yes", "yes", TRUE), 1,
                  ifelse(dm %in% c(0, "0", "No", "no", FALSE, ""), 0, NA))
# Also check IFGIGT for pre-diabetes
if ("IFGIGT" %in% names(wales_raw)) {
    has_ifg <- wales_raw$IFGIGT %in% c(1, "1", "Yes", "yes")
    # Keep as separate or combine — for now, diabetes only means actual DM
}
cat(sprintf("  diabetes: %d non-missing (%d diabetic)\n",
            sum(!is.na(wales$diabetes)), sum(wales$diabetes == 1, na.rm = TRUE)))

# -- 3F: Hypertension (from BP medication — 100% available) --
bp_med <- wales_raw$BloodPressureMedication
wales$hypertension <- ifelse(bp_med %in% c(1, "1", "Yes", "yes", TRUE), 1,
                      ifelse(bp_med %in% c(0, "0", "No", "no", FALSE, ""), 0, NA))
# Also flag if SBP >= 140 or DBP >= 90
sbp <- as.numeric(wales_raw$BloodPressureSystolic)
dbp <- as.numeric(wales_raw$BloodPressureDiastolic)
high_bp <- (!is.na(sbp) & sbp >= 140) | (!is.na(dbp) & dbp >= 90)
wales$hypertension[high_bp & (is.na(wales$hypertension) | wales$hypertension == 0)] <- 1
cat(sprintf("  hypertension: %d non-missing (%d hypertensive)\n",
            sum(!is.na(wales$hypertension)), sum(wales$hypertension == 1, na.rm = TRUE)))

# -- 3G: BMI --
wales$bmi <- as.numeric(wales_raw$BMI)
cat(sprintf("  bmi: %d non-missing\n", sum(!is.na(wales$bmi))))

# -- 3H: Blood pressure --
wales$sbp <- sbp
wales$dbp <- dbp
cat(sprintf("  sbp: %d non-missing\n", sum(!is.na(wales$sbp))))
cat(sprintf("  dbp: %d non-missing\n", sum(!is.na(wales$dbp))))

# -- 3I: Gene (APOB vs LDLR/PCSK9) --
mut <- wales_raw$Mutation1
wales$gene_apob <- as.integer(grepl("APOB", mut, ignore.case = TRUE))
cat(sprintf("  gene_apob: %d APOB carriers (of %d with mutation data)\n",
            sum(wales$gene_apob == 1, na.rm = TRUE), sum(mut != "" & !is.na(mut))))

# -- 3J: Deprivation (WIMD = Welsh equivalent of Townsend) --
wales$townsend <- as.numeric(wales_raw$WIMD2025OverallDecile)
# WIMD decile: 1 = most deprived, 10 = least deprived → invert for Townsend-like
wales$townsend <- ifelse(!is.na(wales$townsend), 11 - wales$townsend, NA)
cat(sprintf("  townsend (from WIMD): %d non-missing\n", sum(!is.na(wales$townsend))))

# -- 3K: Outcome --
wales$ascvd_combined <- as.integer(wales_raw$ascvd_combine %in% c(1, "1", "1.0"))
cat(sprintf("\n  ASCVD events: %d / %d (%.1f%%)\n",
            sum(wales$ascvd_combined), nrow(wales),
            100 * mean(wales$ascvd_combined)))

# =============================================================================
# SECTION 4: CREATE DERIVED FEATURES (both UKB and Wales)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 4: DERIVE FEATURES (both datasets)\n")
cat("===================================================================\n\n")

derive_features <- function(d, name) {
    cat(sprintf("  [%s] Deriving features...\n", name))

    # Non-HDL cholesterol
    d$non_hdl <- d$tc - d$hdl

    # TC/HDL ratio
    d$tc_hdl_ratio <- ifelse(d$hdl > 0, d$tc / d$hdl, NA)

    # log(TG/HDL)
    tg_hdl <- ifelse(d$hdl > 0 & d$trig > 0, d$trig / d$hdl, NA)
    d$log_tg_hdl <- ifelse(!is.na(tg_hdl) & tg_hdl > 0, log(tg_hdl), NA)

    # Pulse pressure (if available)
    d$pulse_pressure <- ifelse(!is.na(d$sbp) & !is.na(d$dbp), d$sbp - d$dbp, NA)

    # Age squared (non-linear age effect)
    d$age_sq <- d$age^2

    # Sex × Age interaction
    d$sex_x_age <- d$sex * d$age

    n_derived <- sum(!is.na(d$non_hdl)) + sum(!is.na(d$tc_hdl_ratio)) +
                 sum(!is.na(d$log_tg_hdl))
    cat(sprintf("    non_hdl: %d, tc_hdl_ratio: %d, log_tg_hdl: %d\n",
                sum(!is.na(d$non_hdl)), sum(!is.na(d$tc_hdl_ratio)),
                sum(!is.na(d$log_tg_hdl))))
    cat(sprintf("    pulse_pressure: %d, age_sq: %d, sex_x_age: %d\n",
                sum(!is.na(d$pulse_pressure)), sum(!is.na(d$age_sq)),
                sum(!is.na(d$sex_x_age))))

    return(d)
}

ukb   <- derive_features(ukb, "UKB")
wales <- derive_features(wales, "Wales")

# -- Derive gene_apob for UKB (was missing) --
if ("gene" %in% names(ukb) && !("gene_apob" %in% names(ukb))) {
    ukb$gene_apob <- as.integer(ukb$gene == "APOB")
    cat(sprintf("  [UKB] gene_apob created: %d APOB carriers\n", sum(ukb$gene_apob, na.rm = TRUE)))
} else if (!("gene_apob" %in% names(ukb))) {
    ukb$gene_apob <- 0
    cat("  [UKB] gene_apob: gene column not found, set to 0\n")
}

# -- Force ALL model variables to numeric in both datasets --
all_model_vars <- c("age", "sex", "re_ldl", "hdl", "tc", "trig",
                     "non_hdl", "tc_hdl_ratio", "log_tg_hdl",
                     "smoking_binary", "diabetes", "hypertension",
                     "pulse_pressure", "bmi", "gene_apob",
                     "age_sq", "sex_x_age", "sbp", "dbp",
                     "ascvd_combined")

cat("\n  Forcing all model variables to numeric...\n")
for (v in all_model_vars) {
    if (v %in% names(ukb))   ukb[[v]]   <- suppressWarnings(as.numeric(ukb[[v]]))
    if (v %in% names(wales)) wales[[v]] <- suppressWarnings(as.numeric(wales[[v]]))
}

# Verify no character columns remain in model variables
cat("  UKB column types: ")
for (v in all_model_vars[all_model_vars %in% names(ukb)]) {
    if (!is.numeric(ukb[[v]])) cat(sprintf("%s=%s ", v, class(ukb[[v]])))
}
cat("OK\n")
cat("  Wales column types: ")
for (v in all_model_vars[all_model_vars %in% names(wales)]) {
    if (!is.numeric(wales[[v]])) cat(sprintf("%s=%s ", v, class(wales[[v]])))
}
cat("OK\n")

# =============================================================================
# SECTION 5: DEFINE PORTABLE MODEL TIERS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 5: PORTABLE MODEL TIERS (Wales-compatible only)\n")
cat("===================================================================\n\n")

# MODEL A: Lean Portable (maximise N — skip BMI, SBP/DBP)
# Uses BP_med as hypertension proxy (100% available in Wales)
model_a_vars <- c("age", "sex", "re_ldl", "hdl",
                   "non_hdl", "tc_hdl_ratio", "log_tg_hdl",
                   "smoking_binary", "diabetes", "hypertension")

# MODEL B: Enriched Portable (adds BP values, BMI, gene — loses ~60% patients)
model_b_vars <- c(model_a_vars, "pulse_pressure", "bmi", "gene_apob")

# MODEL C: Extended Portable (adds non-linear and interaction terms)
model_c_vars <- c(model_b_vars, "age_sq", "sex_x_age")

# Check availability
cat("  Variable availability check:\n")
cat(sprintf("  %-20s %-8s %-8s\n", "Variable", "UKB", "Wales"))
cat(paste(rep("-", 40), collapse = ""), "\n")
for (v in model_c_vars) {
    ukb_n <- sum(!is.na(ukb[[v]]))
    wales_n <- sum(!is.na(wales[[v]]))
    cat(sprintf("  %-20s %6d   %6d\n", v,
                ifelse(is.null(ukb_n), 0, ukb_n),
                ifelse(is.null(wales_n), 0, wales_n)))
}

# Complete cases for each model
for (model_name in c("A: Lean Portable", "B: Enriched Portable", "C: Extended Portable")) {
    mv <- switch(model_name,
                 "A: Lean Portable" = model_a_vars,
                 "B: Enriched Portable" = model_b_vars,
                 "C: Extended Portable" = model_c_vars)

    ukb_cc  <- sum(complete.cases(ukb[, c(mv[mv %in% names(ukb)], "ascvd_combined")]))
    wales_cc <- sum(complete.cases(wales[, c(mv[mv %in% names(wales)], "ascvd_combined")]))
    ukb_ev  <- sum(ukb$ascvd_combined[complete.cases(ukb[, c(mv[mv %in% names(ukb)], "ascvd_combined")])])
    wales_ev <- sum(wales$ascvd_combined[complete.cases(wales[, c(mv[mv %in% names(wales)], "ascvd_combined")])])

    cat(sprintf("\n  %-25s UKB: N=%d, Ev=%d | Wales: N=%d, Ev=%d\n",
                model_name, ukb_cc, ukb_ev, wales_cc, wales_ev))
}

# =============================================================================
# SECTION 6: DEVELOP PORTABLE MODELS ON UKB
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 6: DEVELOP ON UKB (Elastic Net)\n")
cat("===================================================================\n\n")

alphas <- c(0, 0.25, 0.5, 0.75, 1.0)
alpha_names <- c("Ridge", "EN(.25)", "EN(.50)", "EN(.75)", "LASSO")

develop_en <- function(name, vars, data) {
    vars_exist <- vars[vars %in% names(data)]
    df_cc <- data[complete.cases(data[, c(vars_exist, "ascvd_combined")]),
                  c(vars_exist, "ascvd_combined")]

    cat(sprintf("\n  %s: N=%d, Events=%d, Vars=%d\n",
                name, nrow(df_cc), sum(df_cc$ascvd_combined), length(vars_exist)))

    if (nrow(df_cc) < 200 || sum(df_cc$ascvd_combined) < 30) {
        cat("    SKIPPED (too few cases/events)\n")
        return(NULL)
    }

    # Build numeric matrix manually (no apply/data.matrix — avoids Posit Cloud conflicts)
    n_obs <- nrow(df_cc)
    n_var <- length(vars_exist)
    X <- matrix(NA_real_, nrow = n_obs, ncol = n_var)
    colnames(X) <- vars_exist
    for (j in seq_along(vars_exist)) {
        X[, j] <- suppressWarnings(as.double(df_cc[[ vars_exist[j] ]]))
    }
    y <- as.double(df_cc$ascvd_combined)

    # Remove any rows with NAs after conversion
    good <- complete.cases(X) & !is.na(y)
    if (sum(!good) > 0) {
        cat(sprintf("    Removed %d rows with NAs\n", sum(!good)))
        X <- X[good, , drop = FALSE]
        y <- y[good]
        df_cc <- df_cc[good, , drop = FALSE]
    }
    cat(sprintf("    Matrix: %d x %d\n", nrow(X), ncol(X)))

    # Save means and SDs (explicit loop — no apply/colMeans)
    X_means <- setNames(numeric(ncol(X)), colnames(X))
    X_sds   <- setNames(numeric(ncol(X)), colnames(X))
    for (j in seq_len(ncol(X))) {
        X_means[j] <- mean(X[, j], na.rm = TRUE)
        X_sds[j]   <- stats::sd(X[, j], na.rm = TRUE)
        if (is.na(X_sds[j]) || X_sds[j] == 0) X_sds[j] <- 1
    }

    best <- NULL
    for (i in seq_along(alphas)) {
        set.seed(2026)
        cv_fit <- tryCatch(
            glmnet::cv.glmnet(X, y, family = "binomial",
                      alpha = alphas[i], nfolds = 10, type.measure = "auc"),
            error = function(e) {
                cat(sprintf("    %-8s FAILED: %s\n", alpha_names[i], e$message))
                NULL
            })

        if (!is.null(cv_fit)) {
            auc_max <- max(cv_fit$cvm)
            nz <- sum(coef(cv_fit, s = "lambda.min")[-1] != 0)
            cat(sprintf("    %-8s CV-AUC=%.4f (%d non-zero vars)\n", alpha_names[i], auc_max, nz))

            if (is.null(best) || auc_max > best$auc) {
                best <- list(alpha = alphas[i], auc = auc_max, fit = cv_fit,
                             n_vars = nz, X = X, y = y, vars = vars_exist,
                             df_cc = df_cc, n = nrow(df_cc), n_ev = sum(y),
                             X_means = X_means, X_sds = X_sds)
            }
        }
    }

    if (is.null(best)) {
        cat("    ALL ALPHAS FAILED — no model produced\n")
        return(NULL)
    }

    # Extract frozen coefficients
    coefs <- as.matrix(coef(best$fit, s = "lambda.min"))
    best$frozen_intercept <- coefs[1, 1]
    best$frozen_betas <- coefs[-1, 1]
    names(best$frozen_betas) <- vars_exist

    cat(sprintf("    BEST: %s, CV-AUC = %.4f\n", alpha_names[which(alphas == best$alpha)], best$auc))
    cat(sprintf("    Frozen coefficients:\n"))
    cat(sprintf("      (Intercept) = %.6f\n", best$frozen_intercept))
    for (j in seq_along(best$frozen_betas)) {
        if (best$frozen_betas[j] != 0) {
            cat(sprintf("      %-20s = %+.6f\n", names(best$frozen_betas)[j], best$frozen_betas[j]))
        }
    }

    return(best)
}

results <- list()
results[["A"]] <- develop_en("Model A: Lean Portable (10 vars)", model_a_vars, ukb)
results[["B"]] <- develop_en("Model B: Enriched Portable (13 vars)", model_b_vars, ukb)
results[["C"]] <- develop_en("Model C: Extended Portable (15 vars)", model_c_vars, ukb)

# =============================================================================
# SECTION 7: XGBOOST — PORTABLE TIERS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 7: DEVELOP ON UKB (XGBoost)\n")
cat("===================================================================\n\n")

develop_xgb <- function(name, vars, data) {
    vars_exist <- vars[vars %in% names(data)]
    df_cc <- data[complete.cases(data[, c(vars_exist, "ascvd_combined")]),
                  c(vars_exist, "ascvd_combined")]

    cat(sprintf("\n  %s: N=%d, Events=%d\n", name, nrow(df_cc), sum(df_cc$ascvd_combined)))

    if (nrow(df_cc) < 200) { cat("    SKIPPED\n"); return(NULL) }

    # Build numeric matrix manually (no data.matrix)
    n_obs <- nrow(df_cc)
    X <- matrix(NA_real_, nrow = n_obs, ncol = length(vars_exist))
    colnames(X) <- vars_exist
    for (j in seq_along(vars_exist)) {
        X[, j] <- suppressWarnings(as.double(df_cc[[ vars_exist[j] ]]))
    }
    y <- as.double(df_cc$ascvd_combined)
    good <- complete.cases(X) & !is.na(y)
    if (sum(!good) > 0) {
        X <- X[good, , drop = FALSE]; y <- y[good]; df_cc <- df_cc[good, , drop = FALSE]
    }
    dtrain <- xgboost::xgb.DMatrix(data = X, label = y)

    grid <- expand.grid(
        max_depth = c(3, 4, 5), eta = c(0.01, 0.03, 0.05),
        min_child_weight = c(5, 10, 20), subsample = c(0.7, 0.8),
        colsample_bytree = c(0.7, 0.8, 1.0), stringsAsFactors = FALSE)
    set.seed(2026)
    grid <- grid[sample(nrow(grid), min(30, nrow(grid))), ]

    best_auc <- 0; best_params <- NULL; best_nr <- NULL

    for (g in 1:nrow(grid)) {
        params <- list(objective = "binary:logistic", eval_metric = "auc",
                       max_depth = grid$max_depth[g], eta = grid$eta[g],
                       min_child_weight = grid$min_child_weight[g],
                       subsample = grid$subsample[g],
                       colsample_bytree = grid$colsample_bytree[g], nthread = 1)

        set.seed(2026)
        cv_res <- tryCatch(
            xgboost::xgb.cv(params = params, data = dtrain, nrounds = 2000, nfold = 10,
                   early_stopping_rounds = 50, verbose = 0, print_every_n = 0),
            error = function(e) tryCatch(
                xgboost::xgb.cv(params = params, data = dtrain, nrounds = 2000, nfold = 10,
                       early_stopping_rounds = 50, verbose = 0),
                error = function(e2) NULL))

        if (!is.null(cv_res)) {
            eval_log <- cv_res$evaluation_log
            auc_col <- grep("test.*auc.*mean", names(eval_log), value = TRUE)[1]
            if (!is.null(auc_col) && !is.na(auc_col)) {
                this_auc <- max(eval_log[[auc_col]], na.rm = TRUE)
                if (this_auc > best_auc) {
                    best_auc <- this_auc; best_params <- params
                    best_nr <- which.max(eval_log[[auc_col]])
                }
            }
        }
    }

    if (is.null(best_params)) { cat("    XGB FAILED\n"); return(NULL) }

    cat(sprintf("    BEST: depth=%d, eta=%.3f, nrounds=%d → CV-AUC=%.4f\n",
                best_params$max_depth, best_params$eta, best_nr, best_auc))

    set.seed(2026)
    final <- xgboost::xgb.train(params = best_params, data = dtrain, nrounds = best_nr, verbose = 0)

    return(list(auc = best_auc, model = final, vars = vars_exist,
                X = X, y = y, n = nrow(df_cc), n_ev = sum(y)))
}

xgb_results <- list()
xgb_results[["A"]] <- develop_xgb("XGB Model A: Lean Portable", model_a_vars, ukb)
xgb_results[["B"]] <- develop_xgb("XGB Model B: Enriched Portable", model_b_vars, ukb)

# =============================================================================
# SECTION 8: SAFEHEART-RE BENCHMARK (both datasets)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 8: SAFEHEART-RE BENCHMARK\n")
cat("===================================================================\n\n")

safeheart_predict <- function(d, name) {
    lp <- -7.053 +
        0.064 * d$age +
        0.775 * d$sex +
        0.109 * d$re_ldl +
        0.431 * d$hypertension +
        0.025 * d$bmi +
        0.466 * d$smoking_binary
    # Note: prior_cvd coefficient (1.414) excluded

    pred <- 1 / (1 + exp(-lp))
    valid <- !is.na(pred) & !is.na(d$ascvd_combined)

    if (sum(valid) < 50) {
        cat(sprintf("  %s: SAFEHEART-RE — too few valid cases (%d)\n", name, sum(valid)))
        return(NULL)
    }

    roc_obj <- pROC::roc(d$ascvd_combined[valid], pred[valid], quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    cat(sprintf("  %s: SAFEHEART-RE AUC = %.4f (N=%d, Events=%d)\n",
                name, auc_val, sum(valid), sum(d$ascvd_combined[valid])))
    return(list(auc = auc_val, pred = pred, valid = valid, n = sum(valid),
                n_ev = sum(d$ascvd_combined[valid])))
}

sh_ukb   <- safeheart_predict(ukb, "UKB")
sh_wales <- safeheart_predict(wales, "Wales")

# =============================================================================
# SECTION 9: EXTERNAL VALIDATION ON WALES (frozen coefficients)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 9: EXTERNAL VALIDATION ON WALES\n")
cat("===================================================================\n\n")

validate_en <- function(model, wales_data, model_name) {
    if (is.null(model)) { cat(sprintf("  %s: no model to validate\n", model_name)); return(NULL) }

    vars <- model$vars
    vars_in_wales <- vars[vars %in% names(wales_data)]
    missing_vars <- setdiff(vars, names(wales_data))

    if (length(missing_vars) > 0) {
        cat(sprintf("  %s: MISSING in Wales: %s\n", model_name, paste(missing_vars, collapse = ", ")))
    }

    if (length(vars_in_wales) < length(vars)) {
        cat(sprintf("  %s: SKIPPED — %d/%d vars missing\n",
                    model_name, length(missing_vars), length(vars)))
        return(NULL)
    }

    # Complete cases in Wales
    wales_cc <- wales_data[complete.cases(wales_data[, c(vars_in_wales, "ascvd_combined")]),
                            c(vars_in_wales, "ascvd_combined")]

    if (nrow(wales_cc) < 100) {
        cat(sprintf("  %s: SKIPPED — only %d complete cases in Wales\n", model_name, nrow(wales_cc)))
        return(NULL)
    }

    # Build numeric matrix manually (same approach as develop_en)
    X_wales <- matrix(NA_real_, nrow = nrow(wales_cc), ncol = length(vars_in_wales))
    colnames(X_wales) <- vars_in_wales
    for (j in seq_along(vars_in_wales)) {
        X_wales[, j] <- suppressWarnings(as.double(wales_cc[[ vars_in_wales[j] ]]))
    }
    y_wales <- as.double(wales_cc$ascvd_combined)

    # Apply frozen coefficients
    pred_wales <- as.numeric(predict(model$fit, X_wales, s = "lambda.min", type = "response"))

    # AUC
    roc_obj <- pROC::roc(y_wales, pred_wales, quiet = TRUE)
    ext_auc <- as.numeric(pROC::auc(roc_obj))
    ci_obj  <- pROC::ci.auc(roc_obj)

    cat(sprintf("  %s → Wales external AUC = %.4f (95%% CI: %.4f–%.4f)\n",
                model_name, ext_auc, ci_obj[1], ci_obj[3]))
    cat(sprintf("    N = %d, Events = %d (%.1f%%)\n",
                nrow(wales_cc), sum(y_wales), 100 * mean(y_wales)))

    # CALIBRATION — Hosmer-Lemeshow
    cal_groups <- 10
    pred_cut <- cut(pred_wales, breaks = quantile(pred_wales, probs = seq(0, 1, 1/cal_groups)),
                    include.lowest = TRUE)
    cal_df <- data.frame(obs = y_wales, pred = pred_wales, grp = pred_cut)
    cal_summary <- aggregate(cbind(obs, pred) ~ grp, data = cal_df, FUN = mean)
    cal_summary$n <- as.numeric(table(pred_cut))
    cal_summary$obs_events <- aggregate(obs ~ grp, data = cal_df, FUN = sum)$obs

    cat(sprintf("\n    CALIBRATION (observed vs predicted by decile):\n"))
    cat(sprintf("    %-10s %8s %8s %8s %8s\n", "Decile", "N", "Obs%", "Pred%", "O/E"))
    cat(paste(rep("-", 50), collapse = ""), "\n")
    for (r in 1:nrow(cal_summary)) {
        oe <- ifelse(cal_summary$pred[r] > 0, cal_summary$obs[r] / cal_summary$pred[r], NA)
        cat(sprintf("    %-10d %8d %8.1f %8.1f %8.2f\n",
                    r, cal_summary$n[r],
                    100 * cal_summary$obs[r], 100 * cal_summary$pred[r], oe))
    }

    # Calibration slope and intercept
    cal_model <- glm(y_wales ~ log(pred_wales / (1 - pred_wales)), family = binomial)
    cal_intercept <- coef(cal_model)[1]
    cal_slope     <- coef(cal_model)[2]
    cat(sprintf("\n    Calibration intercept: %.4f (ideal = 0)\n", cal_intercept))
    cat(sprintf("    Calibration slope:     %.4f (ideal = 1)\n", cal_slope))

    # RECALIBRATION — update intercept only
    offset_lp <- log(pred_wales / (1 - pred_wales))
    recal_fit <- glm(y_wales ~ offset(offset_lp), family = binomial)
    recal_intercept <- coef(recal_fit)[1]
    pred_recal <- 1 / (1 + exp(-(offset_lp + recal_intercept)))

    roc_recal <- pROC::roc(y_wales, pred_recal, quiet = TRUE)
    cat(sprintf("\n    After recalibration (intercept update = %+.4f):\n", recal_intercept))
    cat(sprintf("    AUC unchanged = %.4f (recalibration only affects calibration, not discrimination)\n",
                as.numeric(pROC::auc(roc_recal))))

    return(list(ext_auc = ext_auc, ci = ci_obj, n = nrow(wales_cc), n_ev = sum(y_wales),
                cal_intercept = cal_intercept, cal_slope = cal_slope,
                recal_intercept = recal_intercept, pred = pred_wales, y = y_wales,
                pred_recal = pred_recal))
}

val_a <- validate_en(results[["A"]], wales, "Model A: Lean Portable")
val_b <- validate_en(results[["B"]], wales, "Model B: Enriched Portable")
val_c <- validate_en(results[["C"]], wales, "Model C: Extended Portable")

# =============================================================================
# SECTION 10: GRAND COMPARISON TABLE
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 10: GRAND COMPARISON\n")
cat("===================================================================\n\n")

cat(sprintf("  %-40s %8s %8s %8s %6s %6s\n",
            "Model", "UKB-AUC", "Wales-AUC", "CI-low", "N(W)", "Ev(W)"))
cat(paste(rep("-", 85), collapse = ""), "\n")

# EN models
for (nm in c("A", "B", "C")) {
    model <- results[[nm]]
    val   <- switch(nm, "A" = val_a, "B" = val_b, "C" = val_c)
    label <- switch(nm,
                    "A" = "EN Model A: Lean Portable",
                    "B" = "EN Model B: Enriched Portable",
                    "C" = "EN Model C: Extended Portable")

    if (!is.null(model)) {
        ukb_auc <- model$auc
        w_auc   <- ifelse(!is.null(val), sprintf("%.4f", val$ext_auc), "N/A")
        ci_low  <- ifelse(!is.null(val), sprintf("%.4f", val$ci[1]), "N/A")
        n_w     <- ifelse(!is.null(val), val$n, 0)
        ev_w    <- ifelse(!is.null(val), val$n_ev, 0)
        cat(sprintf("  %-40s %8.4f %8s %8s %6d %6d\n", label, ukb_auc, w_auc, ci_low, n_w, ev_w))
    }
}

# XGBoost models
for (nm in c("A", "B")) {
    xgb <- xgb_results[[nm]]
    label <- switch(nm,
                    "A" = "XGB Model A: Lean Portable",
                    "B" = "XGB Model B: Enriched Portable")
    if (!is.null(xgb)) {
        cat(sprintf("  %-40s %8.4f %8s %8s %6d %6d\n",
                    label, xgb$auc, "(N/A)", "(N/A)", 0, 0))
    }
}

# SAFEHEART-RE
if (!is.null(sh_ukb)) {
    w_auc <- ifelse(!is.null(sh_wales), sprintf("%.4f", sh_wales$auc), "N/A")
    cat(sprintf("  %-40s %8.4f %8s %8s %6d %6d\n",
                "SAFEHEART-RE (frozen)", sh_ukb$auc, w_auc, "—",
                ifelse(!is.null(sh_wales), sh_wales$n, 0),
                ifelse(!is.null(sh_wales), sh_wales$n_ev, 0)))
}

# =============================================================================
# SECTION 11: SAVE FROZEN COEFFICIENTS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 11: SAVE FROZEN COEFFICIENTS\n")
cat("===================================================================\n\n")

for (nm in c("A", "B", "C")) {
    model <- results[[nm]]
    if (!is.null(model)) {
        coef_df <- data.frame(
            variable = c("(Intercept)", model$vars),
            coefficient = c(model$frozen_intercept, model$frozen_betas),
            stringsAsFactors = FALSE
        )
        fname <- paste0(TAB_DIR, "CALON2_portable_model_", nm, "_coefficients.csv")
        write.csv(coef_df, fname, row.names = FALSE)
        cat(sprintf("  Saved: %s\n", basename(fname)))
    }
}

# Save summary table
summary_rows <- list()
for (nm in c("A", "B", "C")) {
    model <- results[[nm]]
    val <- switch(nm, "A" = val_a, "B" = val_b, "C" = val_c)
    if (!is.null(model)) {
        summary_rows[[nm]] <- data.frame(
            Model = paste0("Model_", nm),
            Variables = paste(model$vars, collapse = "; "),
            N_vars = length(model$vars),
            UKB_CV_AUC = model$auc,
            UKB_N = model$n,
            UKB_Events = model$n_ev,
            Wales_AUC = ifelse(!is.null(val), val$ext_auc, NA),
            Wales_N = ifelse(!is.null(val), val$n, NA),
            Wales_Events = ifelse(!is.null(val), val$n_ev, NA),
            Cal_Intercept = ifelse(!is.null(val), val$cal_intercept, NA),
            Cal_Slope = ifelse(!is.null(val), val$cal_slope, NA),
            stringsAsFactors = FALSE
        )
    }
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, paste0(TAB_DIR, "CALON2_portable_validation_summary.csv"), row.names = FALSE)
cat(sprintf("  Saved: CALON2_portable_validation_summary.csv\n"))

# =============================================================================
# SECTION 12: VIF CHECK (collinearity)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 12: VIF — COLLINEARITY CHECK\n")
cat("===================================================================\n\n")

for (nm in c("A", "B", "C")) {
    model <- results[[nm]]
    if (!is.null(model)) {
        label <- switch(nm, "A" = "Model A", "B" = "Model B", "C" = "Model C")
        tryCatch({
            vif_fit <- glm(ascvd_combined ~ ., data = model$df_cc, family = binomial)
            vif_vals <- car::vif(vif_fit)
            cat(sprintf("  %s VIF (max = %.2f):\n", label, max(vif_vals)))
            for (j in order(-vif_vals)) {
                flag <- ifelse(vif_vals[j] > 5, " *** HIGH ***", "")
                cat(sprintf("    %-20s VIF = %.2f%s\n", names(vif_vals)[j], vif_vals[j], flag))
            }
            cat("\n")
        }, error = function(e) cat(sprintf("  %s VIF failed: %s\n", label, e$message)))
    }
}

# =============================================================================
# DONE
# =============================================================================

cat("\n===================================================================\n")
cat("DONE — CALON-2 PORTABLE MODEL VALIDATION\n")
cat("===================================================================\n\n")

cat("  OUTPUT FILES:\n")
cat(sprintf("    %s\n", paste0(TAB_DIR, "CALON2_portable_model_A_coefficients.csv")))
cat(sprintf("    %s\n", paste0(TAB_DIR, "CALON2_portable_model_B_coefficients.csv")))
cat(sprintf("    %s\n", paste0(TAB_DIR, "CALON2_portable_model_C_coefficients.csv")))
cat(sprintf("    %s\n", paste0(TAB_DIR, "CALON2_portable_validation_summary.csv")))

cat("\n  NEXT STEPS:\n")
cat("    1. Review which model has best Wales AUC vs parsimony trade-off\n")
cat("    2. If Model A performs similarly to B/C → use A (more patients)\n")
cat("    3. Generate calibration plots for the chosen model\n")
cat("    4. Repeat on DRAGON3 with the chosen model\n\n")

cat(sprintf("  Completed at: %s\n\n", Sys.time()))

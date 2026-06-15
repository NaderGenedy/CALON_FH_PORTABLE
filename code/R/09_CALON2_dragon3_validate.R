################################################################################
#                                                                              #
#  CALON-2 v2: DRAGON3 EXTERNAL VALIDATION                                   #
#  Apply frozen v2-A coefficients (from UKB) to DRAGON3 cohort               #
#                                                                              #
#  DRAGON3 data notes:                                                        #
#    - N = ~424 usable (rest are empty/relatives)                             #
#    - 62 ASCVD events (14.6%)                                                #
#    - BMI essentially unavailable → SAFEHEART uses population mean           #
#    - BP stored as "129/70" format → parsed                                  #
#    - Smoking coded 0/1/2 → recoded to binary (0 vs 1+2)                    #
#    - Lipids use underscore notation: LDL_1, HDL_1 etc                       #
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

required_packages <- c("glmnet", "pROC", "car")

for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
        library(pkg, character.only = TRUE)
    }
}

cat("\n")
cat("======================================================================\n")
cat("  CALON-2 v2: DRAGON3 EXTERNAL VALIDATION                           \n")
cat("  Apply frozen v2-A coefficients -> third external validation        \n")
cat("======================================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

CLOUD_PATH  <- "/cloud/project/"
LOCAL_PATH  <- "C:/Users/nader/Downloads/calon_ukb_pipeline/"

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

OUT_DIR <- paste0(BASE_DIR, "output/")
TAB_DIR <- paste0(OUT_DIR, "tables/")
FIG_DIR <- paste0(OUT_DIR, "figures/")
for (d in c(TAB_DIR, FIG_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# =============================================================================
# SECTION 1: LOAD DRAGON3 DATA
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 1: LOAD DRAGON3 DATA\n")
cat("===================================================================\n\n")

dragon_file <- NULL
search_dirs <- unique(c(BASE_DIR, getwd(), paste0(getwd(), "/"),
                         "/cloud/project/", "/cloud/project/output/"))

for (sd in search_dirs) {
    if (!dir.exists(sd)) next
    found <- list.files(sd, pattern = "DRAGON.*3\\.csv$", full.names = TRUE, ignore.case = TRUE)
    if (length(found) > 0) { dragon_file <- found[1]; break }
}

if (is.null(dragon_file)) {
    stop("ERROR: DRAGON_3.csv not found. Upload to project root.")
}

dragon_raw <- read.csv(dragon_file, stringsAsFactors = FALSE)
cat(sprintf("  DRAGON3 raw: %d rows x %d columns\n", nrow(dragon_raw), ncol(dragon_raw)))
cat(sprintf("  File: %s\n", basename(dragon_file)))

# Filter to patients with data (ASCVD_combined is populated)
has_data <- !is.na(dragon_raw$ASCVD_combined) &
            dragon_raw$ASCVD_combined != "" &
            !is.na(dragon_raw$Currentage) &
            dragon_raw$Currentage != ""
dragon_raw <- dragon_raw[has_data, ]
cat(sprintf("  After filtering to patients with data: %d rows\n", nrow(dragon_raw)))

# =============================================================================
# SECTION 2: HARMONISE DRAGON3 -> UKB VARIABLE NAMES
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 2: HARMONISE DRAGON3 VARIABLES\n")
cat("===================================================================\n\n")

d3 <- data.frame(row.names = 1:nrow(dragon_raw))

# -- 2A: Age --
d3$age <- as.numeric(dragon_raw$Currentage)
cat(sprintf("  age: %d non-missing (median %.0f)\n",
            sum(!is.na(d3$age)), median(d3$age, na.rm = TRUE)))

# -- 2B: Sex --
d3$sex <- ifelse(dragon_raw$Gender == "M", 1,
          ifelse(dragon_raw$Gender == "F", 0, NA))
cat(sprintf("  sex: %d non-missing (%d male, %d female)\n",
            sum(!is.na(d3$sex)), sum(d3$sex == 1, na.rm = TRUE),
            sum(d3$sex == 0, na.rm = TRUE)))

# -- 2C: Lipids (DRAGON3 uses LDL_1, HDL_1, TC_1, TRG_1 with underscores) --
d3$re_ldl <- as.numeric(dragon_raw$LDL_1)
d3$hdl    <- as.numeric(dragon_raw$HDL_1)
d3$tc     <- as.numeric(dragon_raw$TC_1)
d3$trig   <- as.numeric(dragon_raw$TRG_1)

# Fill from later measurements if _1 is missing
for (suffix in c("_2", "_3", "_4")) {
    ldl_col <- paste0("LDL", suffix)
    hdl_col <- paste0("HDL", suffix)
    tc_col  <- paste0("TC", suffix)
    tg_col  <- paste0("TRG", suffix)

    if (ldl_col %in% names(dragon_raw)) {
        fill_idx <- is.na(d3$re_ldl) & !is.na(as.numeric(dragon_raw[[ldl_col]]))
        d3$re_ldl[fill_idx] <- as.numeric(dragon_raw[[ldl_col]][fill_idx])
    }
    if (hdl_col %in% names(dragon_raw)) {
        fill_idx <- is.na(d3$hdl) & !is.na(as.numeric(dragon_raw[[hdl_col]]))
        d3$hdl[fill_idx] <- as.numeric(dragon_raw[[hdl_col]][fill_idx])
    }
    if (tc_col %in% names(dragon_raw)) {
        fill_idx <- is.na(d3$tc) & !is.na(as.numeric(dragon_raw[[tc_col]]))
        d3$tc[fill_idx] <- as.numeric(dragon_raw[[tc_col]][fill_idx])
    }
    if (tg_col %in% names(dragon_raw)) {
        fill_idx <- is.na(d3$trig) & !is.na(as.numeric(dragon_raw[[tg_col]]))
        d3$trig[fill_idx] <- as.numeric(dragon_raw[[tg_col]][fill_idx])
    }
}

# Also try LastLDL/LastHDL/LastTC/LastTrigs as fallback
if ("LastLDL" %in% names(dragon_raw)) {
    fill <- is.na(d3$re_ldl) & !is.na(as.numeric(dragon_raw$LastLDL))
    d3$re_ldl[fill] <- as.numeric(dragon_raw$LastLDL[fill])
}
if ("LastHDL" %in% names(dragon_raw)) {
    fill <- is.na(d3$hdl) & !is.na(as.numeric(dragon_raw$LastHDL))
    d3$hdl[fill] <- as.numeric(dragon_raw$LastHDL[fill])
}
if ("LastTC" %in% names(dragon_raw)) {
    fill <- is.na(d3$tc) & !is.na(as.numeric(dragon_raw$LastTC))
    d3$tc[fill] <- as.numeric(dragon_raw$LastTC[fill])
}
if ("LastTrigs" %in% names(dragon_raw)) {
    fill <- is.na(d3$trig) & !is.na(as.numeric(dragon_raw$LastTrigs))
    d3$trig[fill] <- as.numeric(dragon_raw$LastTrigs[fill])
}

cat(sprintf("  re_ldl: %d non-missing\n", sum(!is.na(d3$re_ldl))))
cat(sprintf("  hdl:    %d non-missing\n", sum(!is.na(d3$hdl))))
cat(sprintf("  tc:     %d non-missing\n", sum(!is.na(d3$tc))))
cat(sprintf("  trig:   %d non-missing\n", sum(!is.na(d3$trig))))

# -- 2D: Smoking (DRAGON3: Smoking_binary 0=never, 1=current, 2=ex) --
smk <- dragon_raw$Smoking_binary
d3$smoking_binary <- ifelse(smk %in% c("1", "2"), 1,
                     ifelse(smk %in% c("0"), 0, NA))
cat(sprintf("  smoking_binary: %d non-missing (%d ever-smokers)\n",
            sum(!is.na(d3$smoking_binary)), sum(d3$smoking_binary == 1, na.rm = TRUE)))

# -- 2E: Diabetes (DRAGON3: Diabetes_binary) --
dm <- dragon_raw$Diabetes_binary
d3$diabetes <- ifelse(dm %in% c(1, "1", "Yes", "yes", TRUE), 1,
               ifelse(dm %in% c(0, "0", "No", "no", FALSE, ""), 0, NA))
cat(sprintf("  diabetes: %d non-missing (%d diabetic)\n",
            sum(!is.na(d3$diabetes)), sum(d3$diabetes == 1, na.rm = TRUE)))

# -- 2F: Hypertension (from onBPtreat + parsed BP) --
bp_treat <- dragon_raw$onBPtreat
d3$hypertension <- ifelse(bp_treat %in% c("Y", "y", "Yes", "yes", "1", 1), 1,
                   ifelse(bp_treat %in% c("N", "n", "No", "no", "0", 0), 0, NA))

# Parse BP from "129/70" format
bp_raw <- dragon_raw$BP
d3$sbp <- NA_real_
d3$dbp <- NA_real_
for (i in seq_len(nrow(dragon_raw))) {
    bp_str <- bp_raw[i]
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

# Also use BloodPressureSystolic/Diastolic columns
bp_sys <- as.numeric(dragon_raw$BloodPressureSystolic)
bp_dia <- as.numeric(dragon_raw$BloodPressureDiastolic)
fill_sbp <- is.na(d3$sbp) & !is.na(bp_sys)
fill_dbp <- is.na(d3$dbp) & !is.na(bp_dia)
d3$sbp[fill_sbp] <- bp_sys[fill_sbp]
d3$dbp[fill_dbp] <- bp_dia[fill_dbp]

# Flag high BP as hypertensive
high_bp <- (!is.na(d3$sbp) & d3$sbp >= 140) | (!is.na(d3$dbp) & d3$dbp >= 90)
d3$hypertension[high_bp & (is.na(d3$hypertension) | d3$hypertension == 0)] <- 1
cat(sprintf("  hypertension: %d non-missing (%d hypertensive)\n",
            sum(!is.na(d3$hypertension)), sum(d3$hypertension == 1, na.rm = TRUE)))
cat(sprintf("  sbp: %d, dbp: %d (parsed from BP + BloodPressure columns)\n",
            sum(!is.na(d3$sbp)), sum(!is.na(d3$dbp))))

# -- 2G: BMI — essentially unavailable in DRAGON3 --
d3$bmi <- suppressWarnings(as.numeric(dragon_raw$BMI))
cat(sprintf("  bmi: %d non-missing (NOTE: essentially unavailable)\n",
            sum(!is.na(d3$bmi) & d3$bmi > 10 & d3$bmi < 80)))
# Set invalid BMI to NA
d3$bmi[is.na(d3$bmi) | d3$bmi < 10 | d3$bmi > 80] <- NA

# -- 2H: Gene --
mut <- dragon_raw$Mutation1
d3$gene_apob <- as.integer(grepl("APOB", mut, ignore.case = TRUE))
cat(sprintf("  gene_apob: %d APOB carriers\n", sum(d3$gene_apob == 1, na.rm = TRUE)))

# -- 2I: Outcome --
d3$ascvd_combined <- as.integer(dragon_raw$ASCVD_combined %in% c(1, "1", "1.0"))
cat(sprintf("\n  ASCVD events: %d / %d (%.1f%%)\n",
            sum(d3$ascvd_combined), nrow(d3),
            100 * mean(d3$ascvd_combined)))

# =============================================================================
# SECTION 3: DERIVE FEATURES
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 3: DERIVE FEATURES\n")
cat("===================================================================\n\n")

# log(TG/HDL)
tg_hdl <- ifelse(d3$hdl > 0 & d3$trig > 0, d3$trig / d3$hdl, NA)
d3$log_tg_hdl <- ifelse(!is.na(tg_hdl) & tg_hdl > 0, log(tg_hdl), NA)

# Pulse pressure
d3$pulse_pressure <- ifelse(!is.na(d3$sbp) & !is.na(d3$dbp), d3$sbp - d3$dbp, NA)

# LDL/HDL ratio
d3$ldl_hdl_ratio <- ifelse(d3$hdl > 0, d3$re_ldl / d3$hdl, NA)

# Age interactions
d3$age_x_ldl <- d3$age * d3$re_ldl
d3$age_x_smoking <- d3$age * d3$smoking_binary
d3$age_sq <- d3$age^2

# log(ApoB/LDL) — atherogenic particle discordance (CALON-C feature)
apob_raw <- suppressWarnings(as.numeric(dragon_raw$ApoB))
d3$apob <- apob_raw
apob_ldl_ratio <- ifelse(!is.na(d3$apob) & !is.na(d3$re_ldl) & d3$re_ldl > 0 & d3$apob > 0,
                         d3$apob / d3$re_ldl, NA)
d3$log_apob_ldl <- ifelse(!is.na(apob_ldl_ratio) & apob_ldl_ratio > 0,
                          log(apob_ldl_ratio), NA)
cat(sprintf("  log_apob_ldl: %d non-missing (from ApoB column)\n",
            sum(!is.na(d3$log_apob_ldl))))

# v1 features (for backward compatibility)
d3$non_hdl <- d3$tc - d3$hdl
d3$tc_hdl_ratio <- ifelse(d3$hdl > 0, d3$tc / d3$hdl, NA)

cat(sprintf("  log_tg_hdl: %d, ldl_hdl_ratio: %d, pulse_pressure: %d\n",
            sum(!is.na(d3$log_tg_hdl)), sum(!is.na(d3$ldl_hdl_ratio)),
            sum(!is.na(d3$pulse_pressure))))

# Force numeric
all_vars <- c("age", "sex", "re_ldl", "hdl", "tc", "trig",
              "non_hdl", "tc_hdl_ratio", "log_tg_hdl",
              "ldl_hdl_ratio", "age_x_ldl", "age_x_smoking",
              "smoking_binary", "diabetes", "hypertension",
              "pulse_pressure", "bmi", "gene_apob",
              "age_sq", "log_apob_ldl", "apob", "ascvd_combined")

for (v in all_vars) {
    if (v %in% names(d3)) d3[[v]] <- suppressWarnings(as.numeric(d3[[v]]))
}

# =============================================================================
# SECTION 4: DEFINE FROZEN MODELS (from v2 UKB development)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 4: FROZEN MODEL COEFFICIENTS\n")
cat("===================================================================\n\n")

# v1-A (original, 10 vars) — from Posit Cloud v2 run
v1a_coefs <- list(
    intercept = -3.655475,
    betas = c(age = 0.056304, sex = 0.165622, re_ldl = 0.156757,
              hdl = -0.730921, non_hdl = -0.182754, tc_hdl_ratio = -0.018701,
              log_tg_hdl = 0.112695, smoking_binary = 0.435609,
              diabetes = 0.591802, hypertension = 0.572546)
)

# v2-A (lean clean, 8 vars — but log_tg_hdl dropped to 0)
v2a_coefs <- list(
    intercept = -4.201192,
    betas = c(age = 0.063350, sex = 0.151549, re_ldl = 0.031590,
              hdl = -0.905700, log_tg_hdl = 0.000000,
              smoking_binary = 0.413440, diabetes = 0.683894,
              hypertension = 0.622197)
)

# v2-A recalibrated (intercept-only recalibration from Wales FH+ validation)
# Section 18 output: alpha_adj = -0.5820, recalibrated intercept = -4.783143
# Brier improvement: -4.9%
v2a_recal_coefs <- list(
    intercept = -4.783143,  # actual from Section 18 recalibration (FH+ only)
    betas = c(age = 0.063350, sex = 0.151549, re_ldl = 0.031590,
              hdl = -0.905700, log_tg_hdl = 0.000000,
              smoking_binary = 0.413440, diabetes = 0.683894,
              hypertension = 0.622197)
)

# v2-B (enriched, 11 vars — re_ldl, log_tg_hdl, pulse_pressure, gene_apob dropped)
v2b_coefs <- list(
    intercept = -4.578069,
    betas = c(age = 0.058228, sex = 0.163033, re_ldl = 0.000000,
              hdl = -0.729041, log_tg_hdl = 0.000000,
              smoking_binary = 0.341169, diabetes = 0.471959,
              hypertension = 0.524514, pulse_pressure = 0.000000,
              bmi = 0.023854, gene_apob = 0.000000)
)

# CALON-C (v2-A + log_apob_ldl) — from script 08 EN(.25) output
# UKB CV-AUC = 0.7390, 9 active vars, max VIF = 2.73
calon_c_coefs <- list(
    intercept = 1.259681,
    betas = c(age = 0.066326, sex = 0.114764, re_ldl = 0.069103,
              hdl = -0.564067, log_tg_hdl = -0.047040,
              smoking_binary = 0.419634, diabetes = 0.463288,
              hypertension = 0.581718,
              log_apob_ldl = 5.131619)
)
# SAFEHEART-RE (published frozen)
safeheart_coefs <- list(
    intercept = -7.053,
    betas = c(age = 0.064, sex = 0.775, re_ldl = 0.109,
              hypertension = 0.431, bmi = 0.025, smoking_binary = 0.466)
)

all_models <- list(
    "v1-A (original)" = v1a_coefs,
    "v2-A (lean clean)" = v2a_coefs,
    "v2-A (recalibrated)" = v2a_recal_coefs,
    "v2-B (enriched)" = v2b_coefs,
    "CALON-C (ApoB)" = calon_c_coefs,
    "SAFEHEART-RE" = safeheart_coefs
)

for (mn in names(all_models)) {
    m <- all_models[[mn]]
    active <- sum(m$betas != 0)
    cat(sprintf("  %s: intercept=%.4f, %d active vars\n", mn, m$intercept, active))
}

# =============================================================================
# SECTION 5: CHECK VARIABLE AVAILABILITY
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 5: VARIABLE AVAILABILITY IN DRAGON3\n")
cat("===================================================================\n\n")

for (mn in names(all_models)) {
    m <- all_models[[mn]]
    vars_needed <- names(m$betas[m$betas != 0])
    available <- vars_needed[vars_needed %in% names(d3)]
    missing <- setdiff(vars_needed, names(d3))

    cat(sprintf("  %s:\n", mn))
    for (v in vars_needed) {
        if (v %in% names(d3)) {
            n_valid <- sum(!is.na(d3[[v]]))
            cat(sprintf("    %-20s %d non-missing (%.0f%%)\n", v, n_valid, 100*n_valid/nrow(d3)))
        } else {
            cat(sprintf("    %-20s MISSING\n", v))
        }
    }
    cat("\n")
}

# =============================================================================
# SECTION 6: APPLY FROZEN MODELS TO DRAGON3
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 6: APPLY FROZEN MODELS TO DRAGON3\n")
cat("===================================================================\n\n")

apply_frozen <- function(data, model_coefs, model_name, impute_bmi = FALSE) {
    vars_needed <- names(model_coefs$betas)
    active_vars <- names(model_coefs$betas[model_coefs$betas != 0])

    # Check all active vars exist
    missing <- setdiff(active_vars, names(data))
    if (length(missing) > 0) {
        cat(sprintf("  %s: MISSING active vars: %s\n", model_name, paste(missing, collapse=", ")))
        return(NULL)
    }

    # Compute linear predictor
    lp <- rep(model_coefs$intercept, nrow(data))
    for (v in active_vars) {
        vals <- data[[v]]
        # Impute BMI with population mean if requested
        if (v == "bmi" && impute_bmi) {
            vals[is.na(vals)] <- 27.0  # FH population mean BMI
            cat(sprintf("    Imputed %d missing BMI values with 27.0\n",
                        sum(is.na(data[[v]]))))
        }
        lp <- lp + model_coefs$betas[v] * vals
    }

    pred <- 1 / (1 + exp(-lp))
    y <- data$ascvd_combined

    # Need all active vars non-missing (except imputed BMI)
    valid <- !is.na(pred) & !is.na(y)
    if (sum(valid) < 30) {
        cat(sprintf("  %s: only %d valid cases — SKIPPED\n", model_name, sum(valid)))
        return(NULL)
    }

    roc_obj <- pROC::roc(y[valid], pred[valid], quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    ci_obj  <- pROC::ci.auc(roc_obj)

    cat(sprintf("  %s -> DRAGON3 AUC = %.4f (95%% CI: %.4f-%.4f)\n",
                model_name, auc_val, ci_obj[1], ci_obj[3]))
    cat(sprintf("    N = %d, Events = %d (%.1f%%)\n",
                sum(valid), sum(y[valid]), 100 * mean(y[valid])))

    # Calibration
    pred_v <- pred[valid]
    y_v <- y[valid]

    cal_groups <- min(10, floor(sum(valid) / 20))
    if (cal_groups >= 5) {
        pred_cut <- cut(pred_v, breaks = quantile(pred_v, probs = seq(0, 1, 1/cal_groups)),
                        include.lowest = TRUE)
        cal_df <- data.frame(obs = y_v, pred = pred_v, grp = pred_cut)
        cal_summary <- aggregate(cbind(obs, pred) ~ grp, data = cal_df, FUN = mean)
        cal_summary$n <- as.numeric(table(pred_cut))

        cat(sprintf("\n    CALIBRATION (%d groups):\n", cal_groups))
        cat(sprintf("    %-6s %6s %8s %8s %6s\n", "Grp", "N", "Obs%", "Pred%", "O/E"))
        for (r in 1:nrow(cal_summary)) {
            oe <- ifelse(cal_summary$pred[r] > 0, cal_summary$obs[r] / cal_summary$pred[r], NA)
            cat(sprintf("    %-6d %6d %8.1f %8.1f %6.2f\n",
                        r, cal_summary$n[r],
                        100 * cal_summary$obs[r], 100 * cal_summary$pred[r], oe))
        }
    }

    # Calibration slope/intercept
    lp_pred <- log(pred_v / (1 - pred_v))
    cal_fit <- tryCatch(glm(y_v ~ lp_pred, family = binomial), error = function(e) NULL)
    cal_int <- NA; cal_slp <- NA
    if (!is.null(cal_fit)) {
        cal_int <- coef(cal_fit)[1]
        cal_slp <- coef(cal_fit)[2]
        cat(sprintf("\n    Cal intercept: %.4f (ideal=0), Cal slope: %.4f (ideal=1)\n",
                    cal_int, cal_slp))
    }

    # Intercept recalibration
    recal_fit <- tryCatch(glm(y_v ~ offset(lp_pred), family = binomial), error = function(e) NULL)
    recal_int <- NA
    if (!is.null(recal_fit)) {
        recal_int <- coef(recal_fit)[1]
        pred_recal <- 1 / (1 + exp(-(lp_pred + recal_int)))
        cat(sprintf("    Recalibration intercept update: %+.4f\n", recal_int))
    }

    return(list(auc = auc_val, ci = ci_obj, n = sum(valid), n_ev = sum(y[valid]),
                cal_intercept = cal_int, cal_slope = cal_slp,
                recal_intercept = recal_int,
                pred = pred[valid], y = y[valid]))
}

results <- list()

# CALON-2 models (no BMI imputation needed)
results[["v1-A"]] <- apply_frozen(d3, v1a_coefs, "EN v1-A (original)")
results[["v2-A"]] <- apply_frozen(d3, v2a_coefs, "EN v2-A (lean clean)")
results[["v2-A-recal"]] <- apply_frozen(d3, v2a_recal_coefs, "EN v2-A (recalibrated)")
results[["v2-B"]] <- apply_frozen(d3, v2b_coefs, "EN v2-B (enriched)", impute_bmi = TRUE)
results[["CALON-C"]] <- apply_frozen(d3, calon_c_coefs, "CALON-C (ApoB)")

# SAFEHEART-RE (needs BMI — impute with 27.0)
results[["SAFEHEART"]] <- apply_frozen(d3, safeheart_coefs, "SAFEHEART-RE", impute_bmi = TRUE)

# =============================================================================
# SECTION 7: HEAD-TO-HEAD COMPARISON
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 7: HEAD-TO-HEAD COMPARISON ON DRAGON3\n")
cat("===================================================================\n\n")

cat(sprintf("  %-30s %8s %12s %6s %6s\n",
            "Model", "AUC", "95% CI", "N", "Ev"))
cat(paste(rep("-", 75), collapse = ""), "\n")

for (mn in c("v1-A", "v2-A", "v2-A-recal", "v2-B", "CALON-C", "SAFEHEART")) {
    r <- results[[mn]]
    if (!is.null(r)) {
        label <- switch(mn,
            "v1-A" = "EN v1-A (original)",
            "v2-A" = "EN v2-A (lean clean)",
            "v2-A-recal" = "EN v2-A (recalibrated)",
            "v2-B" = "EN v2-B (enriched+impBMI)",
            "CALON-C" = "CALON-C (ApoB, 9 vars)",
            "SAFEHEART" = "SAFEHEART-RE (impBMI)")
        cat(sprintf("  %-35s %8.4f %5.4f-%.4f %6d %6d\n",
                    label, r$auc, r$ci[1], r$ci[3], r$n, r$n_ev))
    }
}

# Delta table: every CALON model vs SAFEHEART
cat("\n  Delta AUC vs SAFEHEART-RE:\n")
cat(sprintf("  %-35s %8s %8s %8s\n", "Model", "AUC", "SH AUC", "Delta"))
cat(paste(rep("-", 65), collapse = ""), "\n")
sr <- results[["SAFEHEART"]]
if (!is.null(sr)) {
    for (mn in c("v1-A", "v2-A", "v2-A-recal", "v2-B", "CALON-C")) {
        r <- results[[mn]]
        if (!is.null(r)) {
            label <- switch(mn,
                "v1-A" = "v1-A (original)",
                "v2-A" = "v2-A (lean clean)",
                "v2-A-recal" = "v2-A (recalibrated)",
                "v2-B" = "v2-B (enriched)",
                "CALON-C" = "CALON-C (ApoB)")
            cat(sprintf("  %-35s %8.4f %8.4f %+8.4f%s\n",
                        label, r$auc, sr$auc, r$auc - sr$auc,
                        ifelse(r$auc - sr$auc >= 0.05, " ** TARGET", "")))
        }
    }
}

# =============================================================================
# SECTION 8: DELTA AUC + DeLong TEST
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 8: DELTA AUC (CALON-2 vs SAFEHEART)\n")
cat("===================================================================\n\n")

# Compare CALON models vs SAFEHEART on matched patients
for (calon_name in c("v2-A", "v2-A-recal", "CALON-C", "v1-A")) {
    cr <- results[[calon_name]]
    sr <- results[["SAFEHEART"]]
    if (is.null(cr) || is.null(sr)) next

    # Both predictions are on different subsets — need to match
    # Since apply_frozen uses all patients where vars are non-missing,
    # we need to align. Simplest: recompute on intersection.
    # For now, report both AUCs (different N) and note this.
    cat(sprintf("  %s AUC = %.4f (N=%d)\n", calon_name, cr$auc, cr$n))
    cat(sprintf("  SAFEHEART  AUC = %.4f (N=%d)\n", sr$auc, sr$n))
    cat(sprintf("  Raw delta = %+.4f\n\n", cr$auc - sr$auc))
}

# Matched comparison: recompute on patients valid for BOTH v2-A and SAFEHEART
cat("  Matched-patient comparison (v2-A vs SAFEHEART):\n")

# v2-A active vars: age, sex, re_ldl, hdl, smoking_binary, diabetes, hypertension
# SAFEHEART active vars: age, sex, re_ldl, hypertension, bmi, smoking_binary
# Union (with bmi imputed)

v2a_active <- c("age", "sex", "re_ldl", "hdl", "smoking_binary", "diabetes", "hypertension")
sh_active  <- c("age", "sex", "re_ldl", "hypertension", "smoking_binary")
# bmi is imputed for SAFEHEART

all_needed <- unique(c(v2a_active, sh_active))
matched_valid <- complete.cases(d3[, c(all_needed, "ascvd_combined")])

if (sum(matched_valid) >= 30) {
    d3m <- d3[matched_valid, ]
    cat(sprintf("    Matched N = %d, Events = %d\n", nrow(d3m), sum(d3m$ascvd_combined)))

    # v2-A predictions
    lp_v2a <- v2a_coefs$intercept
    for (v in names(v2a_coefs$betas)) {
        if (v %in% names(d3m) && v2a_coefs$betas[v] != 0) {
            lp_v2a <- lp_v2a + v2a_coefs$betas[v] * d3m[[v]]
        }
    }
    pred_v2a <- 1 / (1 + exp(-lp_v2a))

    # SAFEHEART predictions (BMI imputed)
    bmi_imp <- d3m$bmi
    bmi_imp[is.na(bmi_imp)] <- 27.0
    lp_sh <- -7.053 + 0.064*d3m$age + 0.775*d3m$sex + 0.109*d3m$re_ldl +
             0.431*d3m$hypertension + 0.025*bmi_imp + 0.466*d3m$smoking_binary
    pred_sh <- 1 / (1 + exp(-lp_sh))

    both_valid <- !is.na(pred_v2a) & !is.na(pred_sh) & !is.na(d3m$ascvd_combined)

    if (sum(both_valid) >= 30) {
        roc_v2a <- pROC::roc(d3m$ascvd_combined[both_valid], pred_v2a[both_valid], quiet = TRUE)
        roc_sh  <- pROC::roc(d3m$ascvd_combined[both_valid], pred_sh[both_valid], quiet = TRUE)

        auc_v2a <- as.numeric(pROC::auc(roc_v2a))
        auc_sh  <- as.numeric(pROC::auc(roc_sh))
        delta   <- auc_v2a - auc_sh

        delong <- tryCatch(pROC::roc.test(roc_v2a, roc_sh, method = "delong"),
                           error = function(e) NULL)

        cat(sprintf("    v2-A AUC:      %.4f\n", auc_v2a))
        cat(sprintf("    SAFEHEART AUC: %.4f\n", auc_sh))
        cat(sprintf("    Delta:         %+.4f\n", delta))
        if (!is.null(delong)) {
            cat(sprintf("    DeLong P:      %.4f\n", delong$p.value))
        }

        # Bootstrap CI for delta
        set.seed(2026)
        n_boot <- 2000
        delta_boot <- numeric(n_boot)
        yb_src <- d3m$ascvd_combined[both_valid]
        pv_src <- pred_v2a[both_valid]
        ps_src <- pred_sh[both_valid]

        for (b in 1:n_boot) {
            idx <- sample(length(yb_src), replace = TRUE)
            r1 <- tryCatch(as.numeric(pROC::auc(pROC::roc(yb_src[idx], pv_src[idx], quiet=TRUE))),
                           error = function(e) NA)
            r2 <- tryCatch(as.numeric(pROC::auc(pROC::roc(yb_src[idx], ps_src[idx], quiet=TRUE))),
                           error = function(e) NA)
            delta_boot[b] <- r1 - r2
        }
        delta_ci <- quantile(delta_boot, c(0.025, 0.975), na.rm = TRUE)
        cat(sprintf("    Bootstrap 95%% CI: (%.4f to %.4f)\n", delta_ci[1], delta_ci[2]))

        # --- NRI/IDI ---
        cat("\n  NRI/IDI (v2-A vs SAFEHEART on DRAGON3):\n")
        y_m  <- d3m$ascvd_combined[both_valid]
        p_new <- pred_v2a[both_valid]
        p_old <- pred_sh[both_valid]
        diff_p <- p_new - p_old

        events <- y_m == 1; non_events <- y_m == 0

        # Category-free NRI
        nri_ev <- (sum(diff_p[events] > 0) - sum(diff_p[events] < 0)) / sum(events)
        nri_ne <- (sum(diff_p[non_events] < 0) - sum(diff_p[non_events] > 0)) / sum(non_events)
        nri <- nri_ev + nri_ne
        cat(sprintf("    Category-free NRI: %+.4f (event=%+.4f, non-event=%+.4f)\n",
                    nri, nri_ev, nri_ne))

        # IDI
        idi_ev <- mean(p_new[events]) - mean(p_old[events])
        idi_ne <- mean(p_new[non_events]) - mean(p_old[non_events])
        idi <- idi_ev - idi_ne
        idi_se <- sqrt(var(p_new[events] - p_old[events])/sum(events) +
                       var(p_new[non_events] - p_old[non_events])/sum(non_events))
        idi_p <- 2 * (1 - pnorm(abs(idi / idi_se)))
        cat(sprintf("    IDI: %+.6f (SE=%.6f, P=%.4f)\n", idi, idi_se, idi_p))
    }
}

# =============================================================================
# SECTION 9: SUBGROUP ANALYSIS ON DRAGON3
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 9: SUBGROUP ANALYSIS ON DRAGON3\n")
cat("===================================================================\n\n")

if (sum(matched_valid) >= 30) {
    d3m <- d3[matched_valid, ]
    bmi_imp <- d3m$bmi; bmi_imp[is.na(bmi_imp)] <- 27.0

    lp_v2a <- v2a_coefs$intercept
    for (v in names(v2a_coefs$betas)) {
        if (v %in% names(d3m) && v2a_coefs$betas[v] != 0)
            lp_v2a <- lp_v2a + v2a_coefs$betas[v] * d3m[[v]]
    }
    pred_v2a <- 1 / (1 + exp(-lp_v2a))

    lp_sh <- -7.053 + 0.064*d3m$age + 0.775*d3m$sex + 0.109*d3m$re_ldl +
             0.431*d3m$hypertension + 0.025*bmi_imp + 0.466*d3m$smoking_binary
    pred_sh <- 1 / (1 + exp(-lp_sh))

    subgroups <- list(
        "Age < 50" = d3m$age < 50,
        "Age >= 50" = d3m$age >= 50,
        "Male" = d3m$sex == 1,
        "Female" = d3m$sex == 0,
        "LDL >= 4.0" = d3m$re_ldl >= 4.0,
        "LDL < 4.0" = d3m$re_ldl < 4.0,
        "Ever smoker" = d3m$smoking_binary == 1,
        "Never smoker" = d3m$smoking_binary == 0,
        "Diabetic" = d3m$diabetes == 1,
        "Hypertensive" = d3m$hypertension == 1,
        "HDL < 1.0" = d3m$hdl < 1.0
    )

    cat(sprintf("  %-20s %5s %5s %8s %8s %8s\n",
                "Subgroup", "N", "Ev", "CALON-2", "SAFEHRT", "Delta"))
    cat(paste(rep("-", 60), collapse = ""), "\n")

    sg_results <- list()
    for (sg_name in names(subgroups)) {
        sg <- subgroups[[sg_name]]
        sg[is.na(sg)] <- FALSE
        sg_v <- sg & !is.na(pred_v2a) & !is.na(pred_sh) & !is.na(d3m$ascvd_combined)

        n_sg <- sum(sg_v)
        n_ev <- sum(d3m$ascvd_combined[sg_v])

        if (n_sg >= 20 && n_ev >= 5) {
            roc_n <- tryCatch(pROC::roc(d3m$ascvd_combined[sg_v], pred_v2a[sg_v], quiet=TRUE),
                              error = function(e) NULL)
            roc_o <- tryCatch(pROC::roc(d3m$ascvd_combined[sg_v], pred_sh[sg_v], quiet=TRUE),
                              error = function(e) NULL)

            if (!is.null(roc_n) && !is.null(roc_o)) {
                a_n <- as.numeric(pROC::auc(roc_n))
                a_o <- as.numeric(pROC::auc(roc_o))
                flag <- ifelse(a_n - a_o >= 0.10, " ***", ifelse(a_n - a_o >= 0.05, " **", ""))
                cat(sprintf("  %-20s %5d %5d %8.4f %8.4f %+8.4f%s\n",
                            sg_name, n_sg, n_ev, a_n, a_o, a_n - a_o, flag))

                sg_results[[sg_name]] <- data.frame(
                    subgroup = sg_name, n = n_sg, events = n_ev,
                    calon2_auc = a_n, safeheart_auc = a_o,
                    delta = a_n - a_o, stringsAsFactors = FALSE)
            }
        } else {
            cat(sprintf("  %-20s %5d %5d %8s %8s %8s (too few)\n",
                        sg_name, n_sg, n_ev, "-", "-", "-"))
        }
    }

    if (length(sg_results) > 0) {
        sg_df <- do.call(rbind, sg_results)
        sg_df <- sg_df[order(-sg_df$delta), ]
        write.csv(sg_df, paste0(TAB_DIR, "CALON2_v2_dragon3_subgroups.csv"), row.names = FALSE)
        cat(sprintf("\n  Saved: CALON2_v2_dragon3_subgroups.csv\n"))
    }
}

# =============================================================================
# SECTION 10: THREE-COHORT SUMMARY
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 10: THREE-COHORT SUMMARY\n")
cat("===================================================================\n\n")

cat("  =====================================================================\n")
cat("  CROSS-COHORT COMPARISON: ALL CALON MODELS vs SAFEHEART-RE\n")
cat("  =====================================================================\n\n")

# --- Panel A: AUC by Model x Cohort ---
cat("  Panel A: AUC across cohorts\n\n")
cat(sprintf("  %-25s %10s %10s %10s\n", "Model", "UKB(dev)", "Wales", "DRAGON3"))
cat(paste(rep("-", 60), collapse = ""), "\n")

# Hardcoded from script 08 outputs
ukb_aucs <- c("v1-A" = 0.7250, "v2-A" = 0.7247, "v2-A-recal" = NA,
              "CALON-C" = 0.7390, "SAFEHEART" = 0.7095)
wales_aucs <- c("v1-A" = 0.8262, "v2-A" = 0.8491, "v2-A-recal" = NA,
                "CALON-C" = NA, "SAFEHEART" = 0.8239)

model_labels <- c("v1-A" = "v1-A (original)",
                   "v2-A" = "v2-A (lean clean)",
                   "v2-A-recal" = "v2-A (recalibrated)",
                   "CALON-C" = "CALON-C (ApoB)",
                   "SAFEHEART" = "SAFEHEART-RE")

for (mn in c("v1-A", "v2-A", "v2-A-recal", "CALON-C", "SAFEHEART")) {
    ukb_str <- ifelse(is.na(ukb_aucs[mn]), "   —    ", sprintf("  %.4f  ", ukb_aucs[mn]))
    wal_str <- ifelse(is.na(wales_aucs[mn]), "   —    ", sprintf("  %.4f  ", wales_aucs[mn]))
    d3_str <- "   —    "
    if (!is.null(results[[mn]])) d3_str <- sprintf("  %.4f  ", results[[mn]]$auc)
    cat(sprintf("  %-25s %10s %10s %10s\n", model_labels[mn], ukb_str, wal_str, d3_str))
}

# --- Panel B: Delta vs SAFEHEART ---
cat("\n  Panel B: Delta AUC vs SAFEHEART-RE (target: +0.05)\n\n")
cat(sprintf("  %-25s %10s %10s %10s\n", "Model", "UKB", "Wales", "DRAGON3"))
cat(paste(rep("-", 60), collapse = ""), "\n")

sh_d3 <- ifelse(!is.null(results[["SAFEHEART"]]), results[["SAFEHEART"]]$auc, NA)

for (mn in c("v1-A", "v2-A", "v2-A-recal", "CALON-C")) {
    ukb_d <- ifelse(is.na(ukb_aucs[mn]), "   —    ",
                    sprintf("  %+.4f ", ukb_aucs[mn] - ukb_aucs["SAFEHEART"]))
    wal_d <- ifelse(is.na(wales_aucs[mn]), "   —    ",
                    sprintf("  %+.4f ", wales_aucs[mn] - wales_aucs["SAFEHEART"]))
    d3_d <- "   —    "
    if (!is.null(results[[mn]]) && !is.na(sh_d3)) {
        delta_val <- results[[mn]]$auc - sh_d3
        flag <- ifelse(delta_val >= 0.05, " **", "")
        d3_d <- sprintf("  %+.4f%s", delta_val, flag)
    }
    cat(sprintf("  %-25s %10s %10s %10s\n", model_labels[mn], ukb_d, wal_d, d3_d))
}

cat("\n  ** = meets +0.05 target\n")

# --- Panel C: Best model per cohort ---
cat("\n  Panel C: Best CALON model per cohort\n\n")

best_ukb <- "CALON-C"
best_ukb_delta <- ukb_aucs["CALON-C"] - ukb_aucs["SAFEHEART"]
cat(sprintf("  UKB:     %s  AUC=%.4f  delta=%+.4f\n",
            model_labels[best_ukb], ukb_aucs[best_ukb], best_ukb_delta))

best_wales <- "v2-A"
best_wales_delta <- wales_aucs["v2-A"] - wales_aucs["SAFEHEART"]
cat(sprintf("  Wales:   %s  AUC=%.4f  delta=%+.4f\n",
            model_labels[best_wales], wales_aucs[best_wales], best_wales_delta))

# DRAGON3: find best CALON model
best_d3_mn <- NULL; best_d3_auc <- -Inf
for (mn in c("v1-A", "v2-A", "v2-A-recal", "CALON-C")) {
    if (!is.null(results[[mn]]) && results[[mn]]$auc > best_d3_auc) {
        best_d3_auc <- results[[mn]]$auc
        best_d3_mn <- mn
    }
}
if (!is.null(best_d3_mn) && !is.na(sh_d3)) {
    cat(sprintf("  DRAGON3: %s  AUC=%.4f  delta=%+.4f\n",
                model_labels[best_d3_mn], best_d3_auc, best_d3_auc - sh_d3))
}

# --- Panel D: Detailed per-model tables ---
cat("\n\n  CALON-2 v2-A (lean clean, 8 vars, 7 active) across all cohorts:\n\n")
cat(sprintf("  %-20s %8s %12s %6s %6s %8s\n",
            "Cohort", "AUC", "95% CI", "N", "Ev", "Prev"))
cat(paste(rep("-", 70), collapse = ""), "\n")

cat(sprintf("  %-20s %8s %12s %6s %6s %8s\n",
            "UKB (development)", "0.7247", "(CV 10-fold)", "1619", "397", "24.5%"))
cat(sprintf("  %-20s %8s %12s %6s %6s %8s\n",
            "Wales (external 1)", "0.8491", "0.8242-0.8739", "1142", "208", "18.2%"))

if (!is.null(results[["v2-A"]])) {
    r <- results[["v2-A"]]
    cat(sprintf("  %-20s %8.4f %5.4f-%.4f %6d %6d %7.1f%%\n",
                "DRAGON3 (external 2)", r$auc, r$ci[1], r$ci[3],
                r$n, r$n_ev, 100 * r$n_ev / r$n))
}

# CALON-C across cohorts
if (!is.null(results[["CALON-C"]])) {
    cat("\n  CALON-C (ApoB, 9 vars) across available cohorts:\n\n")
    cat(sprintf("  %-20s %8s %12s %6s %6s %8s\n",
                "Cohort", "AUC", "95% CI", "N", "Ev", "Prev"))
    cat(paste(rep("-", 70), collapse = ""), "\n")
    cat(sprintf("  %-20s %8s %12s %6s %6s %8s\n",
                "UKB (development)", "0.7390", "(CV 10-fold)", "1591", "~387", "~24%"))
    r <- results[["CALON-C"]]
    cat(sprintf("  %-20s %8.4f %5.4f-%.4f %6d %6d %7.1f%%\n",
                "DRAGON3 (external)", r$auc, r$ci[1], r$ci[3],
                r$n, r$n_ev, 100 * r$n_ev / r$n))
    cat(sprintf("  %-20s %8s %12s %6s %6s %8s\n",
                "Wales", "  —  ", "(no ApoB)", "  —  ", "  —  ", "  —  "))
}

cat("\n  SAFEHEART-RE across all cohorts:\n\n")
cat(sprintf("  %-20s %8s %12s %6s\n", "Cohort", "AUC", "95% CI", "N"))
cat(paste(rep("-", 55), collapse = ""), "\n")
cat(sprintf("  %-20s %8s %12s %6s\n",
            "UKB", "0.7095", "0.6809-0.7381", "1619"))
cat(sprintf("  %-20s %8s %12s %6s\n",
            "Wales", "0.8239", "0.7880-0.8598", "601"))

if (!is.null(results[["SAFEHEART"]])) {
    r <- results[["SAFEHEART"]]
    cat(sprintf("  %-20s %8.4f %5.4f-%.4f %6d\n",
                "DRAGON3 (imp BMI)", r$auc, r$ci[1], r$ci[3], r$n))
}

# =============================================================================
# SECTION 11: SAVE RESULTS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 11: SAVE RESULTS\n")
cat("===================================================================\n\n")

# Save validation summary
summary_rows <- list()
for (mn in names(results)) {
    r <- results[[mn]]
    if (!is.null(r)) {
        summary_rows[[mn]] <- data.frame(
            Model = mn,
            Cohort = "DRAGON3",
            AUC = round(r$auc, 4),
            CI_low = round(r$ci[1], 4),
            CI_high = round(r$ci[3], 4),
            N = r$n,
            Events = r$n_ev,
            Cal_Intercept = ifelse(!is.na(r$cal_intercept), round(r$cal_intercept, 4), NA),
            Cal_Slope = ifelse(!is.na(r$cal_slope), round(r$cal_slope, 4), NA),
            stringsAsFactors = FALSE)
    }
}

if (length(summary_rows) > 0) {
    summary_df <- do.call(rbind, summary_rows)
    write.csv(summary_df, paste0(TAB_DIR, "CALON2_v2_dragon3_validation.csv"), row.names = FALSE)
    cat(sprintf("  Saved: CALON2_v2_dragon3_validation.csv\n"))
}

# =============================================================================
# DONE
# =============================================================================

cat("\n===================================================================\n")
cat("DONE — CALON-2 v2 DRAGON3 VALIDATION\n")
cat("===================================================================\n\n")

cat("  OUTPUT FILES:\n")
cat(sprintf("    %sCALON2_v2_dragon3_validation.csv\n", TAB_DIR))
cat(sprintf("    %sCALON2_v2_dragon3_subgroups.csv\n", TAB_DIR))

cat(sprintf("\n  Completed at: %s\n\n", Sys.time()))

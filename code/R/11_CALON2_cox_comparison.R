##############################################################################
#  11_CALON2_cox_comparison.R
#  Cox PH time-to-event comparison: CALON-2 v2-A vs SAFEHEART-RE
#  Across UKB (development), Wales (external 1), DRAGON3 (external 2)
##############################################################################

cat("\n======================================================================\n")
cat("  CALON-2: COX PROPORTIONAL HAZARDS TIME-TO-EVENT COMPARISON\n")
cat("  v2-A vs SAFEHEART-RE across three FH cohorts\n")
cat("======================================================================\n\n")

# =============================================================================
# SECTION 0: ENVIRONMENT & PACKAGES
# =============================================================================

if (dir.exists("/cloud/project")) {
    BASE_DIR <- "/cloud/project/"
    cat("  Environment: Posit Cloud\n")
} else {
    BASE_DIR <- getwd()
    cat("  Environment: Local\n")
}
cat(sprintf("  BASE_DIR = %s\n\n", BASE_DIR))

OUT_DIR <- paste0(BASE_DIR, "output/tables/")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load packages
for (pkg in c("survival", "pROC")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
}

# =============================================================================
# FROZEN COEFFICIENTS (from script 08 elastic net output)
# =============================================================================

# v2-A (lean clean) — 8 vars, 7 active (log_tg_hdl zeroed)
v2a_coefs <- list(
    intercept = -4.201192,
    betas = c(age = 0.063350, sex = 0.151549, re_ldl = 0.031590,
              hdl = -0.905700, log_tg_hdl = 0.000000,
              smoking_binary = 0.413440, diabetes = 0.683894,
              hypertension = 0.622197)
)

# SAFEHEART-RE (published)
safeheart_coefs <- list(
    intercept = -7.053,
    betas = c(age = 0.064, sex = 0.775, re_ldl = 0.109,
              hypertension = 0.431, bmi = 0.025, smoking_binary = 0.466)
)

# =============================================================================
# HELPER: Compute linear predictor from frozen logistic coefficients
# =============================================================================

compute_lp <- function(data, coefs, impute_bmi = FALSE) {
    active_vars <- names(coefs$betas[coefs$betas != 0])
    missing_vars <- setdiff(active_vars, names(data))
    if (length(missing_vars) > 0) {
        cat(sprintf("    WARNING: Missing vars: %s\n", paste(missing_vars, collapse = ", ")))
        return(rep(NA, nrow(data)))
    }

    lp <- rep(coefs$intercept, nrow(data))
    for (v in active_vars) {
        vals <- data[[v]]
        if (v == "bmi" && impute_bmi) {
            vals[is.na(vals)] <- 27.0
        }
        lp <- lp + coefs$betas[v] * vals
    }
    return(lp)
}

# =============================================================================
# HELPER: Cox analysis using logistic LP as single prognostic index
# =============================================================================

cox_analysis <- function(time, event, lp_calon, lp_safeheart, cohort_name,
                         time_label = "years") {

    cat(sprintf("\n  --- %s ---\n", cohort_name))

    # Complete cases
    valid <- !is.na(time) & !is.na(event) & !is.na(lp_calon) & !is.na(lp_safeheart) &
             time > 0
    if (sum(valid) < 30) {
        cat(sprintf("    Only %d valid cases — SKIPPED\n", sum(valid)))
        return(NULL)
    }

    t <- time[valid]
    e <- event[valid]
    lp_c <- lp_calon[valid]
    lp_s <- lp_safeheart[valid]

    cat(sprintf("    N = %d, Events = %d (%.1f%%)\n", length(t), sum(e), 100 * mean(e)))
    cat(sprintf("    Median follow-up (%s): %.2f\n", time_label, median(t)))
    cat(sprintf("    Follow-up range: %.2f - %.2f %s\n", min(t), max(t), time_label))

    # ---- Cox with CALON-2 LP as single covariate ----
    surv_obj <- Surv(t, e)

    cox_calon <- coxph(surv_obj ~ lp_c)
    cox_safe  <- coxph(surv_obj ~ lp_s)

    # C-index (concordance = discrimination for survival)
    c_calon <- summary(cox_calon)$concordance
    c_safe  <- summary(cox_safe)$concordance

    cat(sprintf("\n    Cox C-index (Harrell's C):\n"))
    cat(sprintf("      CALON-2 v2-A:   %.4f (SE=%.4f)\n", c_calon[1], c_calon[2]))
    cat(sprintf("      SAFEHEART-RE:   %.4f (SE=%.4f)\n", c_safe[1], c_safe[2]))
    cat(sprintf("      Delta:          %+.4f%s\n",
                c_calon[1] - c_safe[1],
                ifelse(c_calon[1] - c_safe[1] >= 0.05, " ** TARGET MET", "")))

    # ---- Cox regression summary ----
    cat(sprintf("\n    Cox model (CALON-2 LP):\n"))
    cat(sprintf("      HR per unit LP: %.4f (95%% CI: %.4f-%.4f)\n",
                exp(coef(cox_calon)), exp(confint(cox_calon))[1], exp(confint(cox_calon))[2]))
    cat(sprintf("      Wald P: %.2e\n", summary(cox_calon)$coefficients[, "Pr(>|z|)"]))
    cat(sprintf("      Likelihood ratio P: %.2e\n", summary(cox_calon)$logtest["pvalue"]))

    cat(sprintf("\n    Cox model (SAFEHEART LP):\n"))
    cat(sprintf("      HR per unit LP: %.4f (95%% CI: %.4f-%.4f)\n",
                exp(coef(cox_safe)), exp(confint(cox_safe))[1], exp(confint(cox_safe))[2]))
    cat(sprintf("      Wald P: %.2e\n", summary(cox_safe)$coefficients[, "Pr(>|z|)"]))
    cat(sprintf("      Likelihood ratio P: %.2e\n", summary(cox_safe)$logtest["pvalue"]))

    # ---- Nested model: does CALON-2 add to SAFEHEART? ----
    cox_both <- coxph(surv_obj ~ lp_s + lp_c)
    lr_test <- anova(cox_safe, cox_both)

    cat(sprintf("\n    Added value of CALON-2 LP over SAFEHEART (nested LR test):\n"))
    cat(sprintf("      LR chi-sq: %.4f, df=%d, P=%.4e\n",
                lr_test[["Chisq"]][2], lr_test[["Df"]][2], lr_test[["Pr(>|Chi|)"]][2]))

    # ---- Risk group comparison (tertiles of LP) ----
    calon_tert <- cut(lp_c, breaks = quantile(lp_c, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                      include.lowest = TRUE, labels = c("Low", "Med", "High"))
    safe_tert <- cut(lp_s, breaks = quantile(lp_s, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                     include.lowest = TRUE, labels = c("Low", "Med", "High"))

    cat(sprintf("\n    Risk tertile event rates:\n"))
    cat(sprintf("    %-8s  CALON-2 (N / Ev / Rate)  SAFEHEART (N / Ev / Rate)\n", "Tertile"))
    cat(paste(rep("-", 68), collapse = ""), "\n")

    for (tert in c("Low", "Med", "High")) {
        c_idx <- calon_tert == tert & !is.na(calon_tert)
        s_idx <- safe_tert == tert & !is.na(safe_tert)

        c_n <- sum(c_idx); c_ev <- sum(e[c_idx]); c_rate <- ifelse(c_n > 0, 100 * c_ev / c_n, NA)
        s_n <- sum(s_idx); s_ev <- sum(e[s_idx]); s_rate <- ifelse(s_n > 0, 100 * s_ev / s_n, NA)

        cat(sprintf("    %-8s  %5d / %3d / %5.1f%%       %5d / %3d / %5.1f%%\n",
                    tert, c_n, c_ev, c_rate, s_n, s_ev, s_rate))
    }

    # ---- Time-dependent AUC at specific timepoints ----
    cat(sprintf("\n    Time-dependent AUC (incident/dynamic):\n"))
    cat(sprintf("    %-12s %8s %8s %8s\n", "Timepoint", "CALON-2", "SAFEHRT", "Delta"))
    cat(paste(rep("-", 45), collapse = ""), "\n")

    # Use concordance at different time horizons via survConcordance
    # Simple approach: calculate AUC among those with event <= t or censored > t
    max_t <- max(t)
    time_cuts <- sort(unique(c(5, 10, 15, 20)))
    time_cuts <- time_cuts[time_cuts < max_t * 0.9]  # only feasible timepoints

    for (tc in time_cuts) {
        # Patients still at risk or had event by tc
        at_risk <- t >= tc | (e == 1 & t <= tc)
        event_by_tc <- e == 1 & t <= tc
        censor_after_tc <- t >= tc

        usable <- at_risk & (event_by_tc | censor_after_tc)
        if (sum(usable) < 30 || sum(event_by_tc[usable]) < 5) next

        y_tc <- as.numeric(event_by_tc[usable])
        roc_c <- tryCatch(pROC::roc(y_tc, lp_c[usable], quiet = TRUE), error = function(e) NULL)
        roc_s <- tryCatch(pROC::roc(y_tc, lp_s[usable], quiet = TRUE), error = function(e) NULL)

        if (!is.null(roc_c) && !is.null(roc_s)) {
            auc_c <- as.numeric(pROC::auc(roc_c))
            auc_s <- as.numeric(pROC::auc(roc_s))
            flag <- ifelse(auc_c - auc_s >= 0.05, " **", "")
            cat(sprintf("    %d-%s %12.4f %8.4f %+8.4f%s\n",
                        tc, time_label, auc_c, auc_s, auc_c - auc_s, flag))
        }
    }

    return(list(
        cohort = cohort_name,
        n = length(t), events = sum(e),
        c_calon = c_calon[1], c_safe = c_safe[1],
        delta = c_calon[1] - c_safe[1],
        lr_chisq = lr_test[["Chisq"]][2],
        lr_p = lr_test[["Pr(>|Chi|)"]][2]
    ))
}


# =============================================================================
# SECTION 1: UKB COHORT (Cox with followup_years)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 1: UKB COHORT — Cox PH Analysis\n")
cat("===================================================================\n")

ukb_file <- NULL
for (sd in c(BASE_DIR, paste0(BASE_DIR, "output/"))) {
    f <- list.files(sd, pattern = "calon2_full_merged\\.csv$", full.names = TRUE, ignore.case = TRUE)
    if (length(f) > 0) { ukb_file <- f[1]; break }
}

ukb_result <- NULL
if (!is.null(ukb_file)) {
    ukb <- read.csv(ukb_file, stringsAsFactors = FALSE)
    cat(sprintf("  Loaded UKB: %d rows, %d cols\n", nrow(ukb), ncol(ukb)))

    # Force numeric
    for (v in c("age", "sex", "re_ldl", "hdl", "log_tg_hdl", "smoking_binary",
                "diabetes", "hypertension", "bmi", "followup_years", "event_indicator")) {
        if (v %in% names(ukb)) ukb[[v]] <- suppressWarnings(as.numeric(ukb[[v]]))
    }

    # Derive log_tg_hdl if missing
    if (!("log_tg_hdl" %in% names(ukb)) || all(is.na(ukb$log_tg_hdl))) {
        tg <- suppressWarnings(as.numeric(ukb$trig))
        hdl <- ukb$hdl
        ratio <- ifelse(!is.na(tg) & !is.na(hdl) & hdl > 0 & tg > 0, tg / hdl, NA)
        ukb$log_tg_hdl <- ifelse(!is.na(ratio) & ratio > 0, log(ratio), NA)
    }

    # Use INCIDENT events only for Cox (not prevalent)
    # event_indicator = ascvd_incident = 211 events
    ukb_time  <- ukb$followup_years
    ukb_event <- ukb$event_indicator

    cat(sprintf("  Incident events: %d / %d (%.1f%%)\n",
                sum(ukb_event == 1, na.rm = TRUE), nrow(ukb),
                100 * mean(ukb_event == 1, na.rm = TRUE)))
    cat(sprintf("  Median follow-up: %.1f years\n", median(ukb_time, na.rm = TRUE)))

    # Compute linear predictors
    lp_v2a <- compute_lp(ukb, v2a_coefs)
    lp_sh  <- compute_lp(ukb, safeheart_coefs, impute_bmi = TRUE)

    ukb_result <- cox_analysis(ukb_time, ukb_event, lp_v2a, lp_sh,
                               "UKB (development)", time_label = "years")
} else {
    cat("  WARNING: UKB file not found\n")
}


# =============================================================================
# SECTION 2: WALES COHORT (Cox with age-at-event)
# =============================================================================

cat("\n\n===================================================================\n")
cat("SECTION 2: WALES COHORT — Cox PH Analysis\n")
cat("===================================================================\n")

wales_file <- NULL
for (sd in c(BASE_DIR, paste0(BASE_DIR, "output/"))) {
    f <- list.files(sd, pattern = "WALES.*FH.*CLEANED.*\\.csv$", full.names = TRUE, ignore.case = TRUE)
    if (length(f) > 0) { wales_file <- f[1]; break }
}

wales_result <- NULL
if (!is.null(wales_file)) {
    wales <- read.csv(wales_file, stringsAsFactors = FALSE)
    cat(sprintf("  Loaded Wales: %d rows, %d cols\n", nrow(wales), ncol(wales)))

    # ── CRITICAL FILTER: Keep only genetically confirmed FH+ patients ──
    fh_positive <- !is.na(wales$Mutation1) & trimws(wales$Mutation1) != ""
    cat(sprintf("  FH+ (genetically confirmed): %d / %d (%.1f%%)\n",
                sum(fh_positive), nrow(wales), 100 * mean(fh_positive)))
    wales <- wales[fh_positive, ]
    cat(sprintf("  Wales after FH+ filter: %d patients\n", nrow(wales)))

    # Force numeric core vars
    num_vars <- c("Currentage", "LDLCholesterol", "HDLCholesterol",
                  "Triglycerides", "TotalCholesterol", "BMI",
                  "MIACSAge", "PCIStentsAge", "CABGAge", "ANGINAAge",
                  "TIAAge", "PVDAge", "OtherAge", "AGE_AT_DECEASED")
    for (v in num_vars) {
        if (v %in% names(wales)) wales[[v]] <- suppressWarnings(as.numeric(wales[[v]]))
    }

    # Harmonise variables
    wales$age <- wales$Currentage
    wales$sex <- ifelse(wales$Gender == "M", 1,
                 ifelse(wales$Gender == "F", 0, NA))
    wales$re_ldl <- wales$LDLCholesterol
    wales$hdl <- wales$HDLCholesterol

    # log_tg_hdl
    tg <- wales$Triglycerides
    hdl <- wales$hdl
    ratio <- ifelse(!is.na(tg) & !is.na(hdl) & hdl > 0 & tg > 0, tg / hdl, NA)
    wales$log_tg_hdl <- ifelse(!is.na(ratio) & ratio > 0, log(ratio), NA)

    # Smoking (column: Smoking, matching script 08 pattern)
    smk <- wales$Smoking
    wales$smoking_binary <- ifelse(smk %in% c(1, "1", "Yes", "yes", "Current", "current",
                                               "Ex", "ex", "Former", "former", "Ex-smoker",
                                               "Current smoker", "Previous"), 1,
                            ifelse(smk %in% c(0, "0", "No", "no", "Never", "never",
                                               "Non-smoker"), 0, NA))

    # Diabetes (column: Diabetes, matching script 08 pattern)
    dm <- wales$Diabetes
    wales$diabetes <- ifelse(dm %in% c(1, "1", "Yes", "yes", TRUE), 1,
                      ifelse(dm %in% c(0, "0", "No", "no", FALSE, ""), 0, NA))

    # Hypertension (column: BloodPressureMedication + BP thresholds, matching script 08)
    bp_med <- wales$BloodPressureMedication
    wales$hypertension <- ifelse(bp_med %in% c(1, "1", "Yes", "yes", TRUE), 1,
                          ifelse(bp_med %in% c(0, "0", "No", "no", FALSE, ""), 0, NA))
    sbp_w <- suppressWarnings(as.numeric(wales$BloodPressureSystolic))
    dbp_w <- suppressWarnings(as.numeric(wales$BloodPressureDiastolic))
    high_bp <- (!is.na(sbp_w) & sbp_w >= 140) | (!is.na(dbp_w) & dbp_w >= 90)
    wales$hypertension[high_bp & (is.na(wales$hypertension) | wales$hypertension == 0)] <- 1

    wales$bmi <- wales$BMI

    # Derive age at first ASCVD event
    event_age_cols <- c("MIACSAge", "PCIStentsAge", "CABGAge", "ANGINAAge", "TIAAge", "PVDAge")
    existing_cols <- event_age_cols[event_age_cols %in% names(wales)]

    if (length(existing_cols) > 0) {
        age_at_first_event <- apply(wales[, existing_cols, drop = FALSE], 1, function(x) {
            vals <- suppressWarnings(as.numeric(x))
            vals <- vals[!is.na(vals) & vals > 10 & vals < 120]  # sanity filter
            if (length(vals) > 0) min(vals) else NA
        })
    } else {
        age_at_first_event <- rep(NA, nrow(wales))
    }

    # Event indicator
    wales$ascvd <- ifelse(!is.na(wales$ascvd_combine) & wales$ascvd_combine == 1, 1, 0)

    # Time variable: age at event (for cases) or current age (for controls)
    # This uses age as the time scale (standard in genetic/FH studies)
    wales$time_age <- ifelse(wales$ascvd == 1 & !is.na(age_at_first_event),
                             age_at_first_event,
                             wales$age)

    # Filter valid
    valid_w <- !is.na(wales$time_age) & wales$time_age > 10 & !is.na(wales$ascvd)
    wales_sub <- wales[valid_w, ]

    cat(sprintf("  Valid for Cox: %d patients, %d events (%.1f%%)\n",
                nrow(wales_sub), sum(wales_sub$ascvd),
                100 * mean(wales_sub$ascvd)))
    cat(sprintf("  Time scale: age at event/censoring\n"))
    cat(sprintf("  Median age: %.1f years\n", median(wales_sub$time_age, na.rm = TRUE)))

    # Compute linear predictors
    lp_v2a_w <- compute_lp(wales_sub, v2a_coefs)
    lp_sh_w  <- compute_lp(wales_sub, safeheart_coefs, impute_bmi = FALSE)

    wales_result <- cox_analysis(wales_sub$time_age, wales_sub$ascvd,
                                 lp_v2a_w, lp_sh_w,
                                 "Wales (external 1)", time_label = "age-years")
} else {
    cat("  WARNING: Wales file not found\n")
}


# =============================================================================
# SECTION 3: DRAGON3 COHORT (Cox with age_at_event_or_censoring)
# =============================================================================

cat("\n\n===================================================================\n")
cat("SECTION 3: DRAGON3 COHORT — Cox PH Analysis\n")
cat("===================================================================\n")

dragon_file <- NULL
for (sd in c(BASE_DIR, paste0(BASE_DIR, "output/"))) {
    f <- list.files(sd, pattern = "DRAGON.*3\\.csv$", full.names = TRUE, ignore.case = TRUE)
    if (length(f) > 0) { dragon_file <- f[1]; break }
}

dragon_result <- NULL
if (!is.null(dragon_file)) {
    dragon_raw <- read.csv(dragon_file, stringsAsFactors = FALSE)
    cat(sprintf("  Loaded DRAGON3: %d rows\n", nrow(dragon_raw)))

    # Filter to patients with data
    has_data <- !is.na(dragon_raw$ASCVD_combined) & dragon_raw$ASCVD_combined != "" &
                !is.na(dragon_raw$Currentage) & dragon_raw$Currentage != ""
    d3 <- dragon_raw[has_data, ]
    cat(sprintf("  Usable patients: %d\n", nrow(d3)))

    # Harmonise (same as script 09)
    d3$age <- suppressWarnings(as.numeric(d3$Currentage))
    d3$sex <- ifelse(tolower(trimws(d3$Sex)) %in% c("male", "m"), 1,
              ifelse(tolower(trimws(d3$Sex)) %in% c("female", "f"), 0, NA))
    d3$re_ldl <- suppressWarnings(as.numeric(d3$LDLCholesterol))
    d3$hdl <- suppressWarnings(as.numeric(d3$HDLCholesterol))
    d3$trig <- suppressWarnings(as.numeric(d3$Triglycerides))
    d3$bmi <- suppressWarnings(as.numeric(d3$BMI))
    d3$bmi[is.na(d3$bmi) | d3$bmi < 10 | d3$bmi > 80] <- NA

    # log_tg_hdl
    ratio <- ifelse(!is.na(d3$trig) & !is.na(d3$hdl) & d3$hdl > 0 & d3$trig > 0,
                    d3$trig / d3$hdl, NA)
    d3$log_tg_hdl <- ifelse(!is.na(ratio) & ratio > 0, log(ratio), NA)

    # Smoking
    smoke_raw <- suppressWarnings(as.numeric(d3$SmokingStatus))
    d3$smoking_binary <- ifelse(!is.na(smoke_raw) & smoke_raw >= 1, 1,
                         ifelse(!is.na(smoke_raw) & smoke_raw == 0, 0, NA))

    # Diabetes
    d3$diabetes <- ifelse(!is.na(d3$DiabetesYear) & d3$DiabetesYear != "" &
                          suppressWarnings(as.numeric(d3$DiabetesYear)) > 0, 1, 0)

    # Hypertension — parse from BP columns
    bp_cols <- c("BP", "BloodPressure")
    d3$sbp <- NA; d3$dbp <- NA
    for (bp_col in bp_cols) {
        if (bp_col %in% names(d3)) {
            for (i in 1:nrow(d3)) {
                if (!is.na(d3[[bp_col]][i]) && grepl("/", d3[[bp_col]][i])) {
                    parts <- strsplit(d3[[bp_col]][i], "/")[[1]]
                    s <- suppressWarnings(as.numeric(trimws(parts[1])))
                    dd <- suppressWarnings(as.numeric(trimws(parts[2])))
                    if (!is.na(s) && s > 50 && s < 300 && is.na(d3$sbp[i])) d3$sbp[i] <- s
                    if (!is.na(dd) && dd > 20 && dd < 200 && is.na(d3$dbp[i])) d3$dbp[i] <- dd
                }
            }
        }
    }
    d3$hypertension <- ifelse(!is.na(d3$sbp) & d3$sbp >= 140, 1,
                       ifelse(!is.na(d3$sbp), 0, NA))

    # Event
    d3$ascvd_combined <- suppressWarnings(as.numeric(d3$ASCVD_combined))

    # Time: age_at_event_or_censoring
    d3$time_age <- suppressWarnings(as.numeric(d3$age_at_event_or_censoring))

    cat(sprintf("  Events: %d / %d (%.1f%%)\n",
                sum(d3$ascvd_combined == 1, na.rm = TRUE), nrow(d3),
                100 * mean(d3$ascvd_combined == 1, na.rm = TRUE)))
    cat(sprintf("  time_age (age_at_event_or_censoring) available: %d (%.0f%%)\n",
                sum(!is.na(d3$time_age)), 100 * mean(!is.na(d3$time_age))))

    # Compute linear predictors
    lp_v2a_d <- compute_lp(d3, v2a_coefs)
    lp_sh_d  <- compute_lp(d3, safeheart_coefs, impute_bmi = TRUE)

    dragon_result <- cox_analysis(d3$time_age, d3$ascvd_combined,
                                  lp_v2a_d, lp_sh_d,
                                  "DRAGON3 (external 2)", time_label = "age-years")
} else {
    cat("  WARNING: DRAGON3 file not found\n")
}


# =============================================================================
# SECTION 4: CROSS-COHORT SUMMARY TABLE
# =============================================================================

cat("\n\n===================================================================\n")
cat("SECTION 4: CROSS-COHORT COX C-INDEX COMPARISON\n")
cat("===================================================================\n\n")

cat("  =====================================================================\n")
cat("  Cox Harrell's C-index: CALON-2 v2-A vs SAFEHEART-RE\n")
cat("  =====================================================================\n\n")

cat(sprintf("  %-22s %6s %6s %8s %8s %8s %10s\n",
            "Cohort", "N", "Ev", "CALON-2", "SAFEHRT", "Delta", "LR P"))
cat(paste(rep("-", 78), collapse = ""), "\n")

all_results <- list(UKB = ukb_result, Wales = wales_result, DRAGON3 = dragon_result)

for (nm in names(all_results)) {
    r <- all_results[[nm]]
    if (!is.null(r)) {
        flag <- ifelse(r$delta >= 0.05, " **", "")
        cat(sprintf("  %-22s %6d %6d %8.4f %8.4f %+8.4f %10.2e%s\n",
                    r$cohort, r$n, r$events,
                    r$c_calon, r$c_safe, r$delta,
                    r$lr_p, flag))
    }
}

cat("\n  ** = delta >= +0.05 target\n")
cat("  LR P = likelihood ratio test P-value for added value of CALON-2 LP over SAFEHEART\n")

# Compare with logistic AUC
cat("\n\n  Comparison: Logistic AUC vs Cox C-index\n\n")
cat(sprintf("  %-22s %10s %10s %10s %10s\n",
            "Cohort", "Log.AUC-C", "Log.AUC-S", "Cox.C-C", "Cox.C-S"))
cat(paste(rep("-", 66), collapse = ""), "\n")

# Logistic AUCs (from script 08/09 outputs)
log_aucs <- data.frame(
    cohort = c("UKB", "Wales", "DRAGON3"),
    calon = c(0.7247, 0.7551, 0.9023),
    safe  = c(0.7095, 0.7211, 0.8952),
    stringsAsFactors = FALSE
)

for (i in 1:nrow(log_aucs)) {
    r <- all_results[[log_aucs$cohort[i]]]
    cox_c <- ifelse(!is.null(r), sprintf("%.4f", r$c_calon), "  —  ")
    cox_s <- ifelse(!is.null(r), sprintf("%.4f", r$c_safe), "  —  ")
    cat(sprintf("  %-22s %10.4f %10.4f %10s %10s\n",
                log_aucs$cohort[i], log_aucs$calon[i], log_aucs$safe[i],
                cox_c, cox_s))
}


# =============================================================================
# SECTION 5: SAVE RESULTS
# =============================================================================

cat("\n\n===================================================================\n")
cat("SECTION 5: SAVE RESULTS\n")
cat("===================================================================\n\n")

summary_rows <- list()
for (nm in names(all_results)) {
    r <- all_results[[nm]]
    if (!is.null(r)) {
        summary_rows[[nm]] <- data.frame(
            Cohort = r$cohort,
            N = r$n, Events = r$events,
            Cox_C_CALON2 = round(r$c_calon, 4),
            Cox_C_SAFEHEART = round(r$c_safe, 4),
            Delta = round(r$delta, 4),
            LR_ChiSq = round(r$lr_chisq, 4),
            LR_P = r$lr_p,
            stringsAsFactors = FALSE)
    }
}

if (length(summary_rows) > 0) {
    summary_df <- do.call(rbind, summary_rows)
    write.csv(summary_df, paste0(OUT_DIR, "CALON2_cox_comparison.csv"), row.names = FALSE)
    cat("  Saved: CALON2_cox_comparison.csv\n")
}

cat("\n===================================================================\n")
cat("DONE — CALON-2 Cox PH Comparison\n")
cat("===================================================================\n\n")
cat(sprintf("  Completed at: %s\n\n", Sys.time()))

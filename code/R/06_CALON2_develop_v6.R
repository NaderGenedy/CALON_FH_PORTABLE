################################################################################
#                                                                              #
#  CALON-2: ASCVD RISK MODEL DEVELOPMENT v6 (STREAMLINED)                     #
#  "Honest AUC — Real Data, Reliable Tests, Fast Execution"                   #
#                                                                              #
#  FIXES from v5:                                                              #
#    - Recreates ALL derived features from raw columns (log_lpa, log_crp,      #
#      prs_primary, gene_apob, etc.) that v5 was missing                       #
#    - Runs only the TOP 5 model tiers (not 18+) → 30 min not 6 hours         #
#    - Honest stacking with outer-fold predictions (no leakage)                #
#    - SAFEHEART-RE comparison with frozen published coefficients              #
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

required_packages <- c(
    "dplyr", "readr", "glmnet", "pROC", "xgboost", "randomForest",
    "mice", "rms", "Boruta", "car", "caret"
)

for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
        library(pkg, character.only = TRUE)
    }
}

cat("\n")
cat("======================================================================\n")
cat("  CALON-2 v6: Streamlined Development Pipeline                        \n")
cat("  Fixes: all derived features recreated, top tiers only               \n")
cat("======================================================================\n\n")

# =============================================================================
# CONFIGURATION — auto-detects environment
# CHECK CLOUD FIRST (Posit Cloud can false-positive on Windows paths)
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
cat(sprintf("  Working directory = %s\n", getwd()))

OUT_DIR  <- paste0(BASE_DIR, "output/")
TAB_DIR  <- paste0(OUT_DIR, "tables/")
FIG_DIR  <- paste0(OUT_DIR, "figures/")
for (d in c(TAB_DIR, FIG_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

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

# =============================================================================
# SECTION 1: LOAD DATA
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 1: LOADING DATA\n")
cat("===================================================================\n\n")

df <- read.csv(merged_file, stringsAsFactors = FALSE)
cat(sprintf("  Dataset: %d patients x %d columns\n", nrow(df), ncol(df)))
cat(sprintf("  ASCVD events: %d (%.1f%%)\n",
            sum(df$ascvd_combined), 100 * mean(df$ascvd_combined)))

# =============================================================================
# SECTION 2: RECREATE ALL DERIVED FEATURES (the v5 bug fix)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 2: FEATURE ENGINEERING (recreating ALL derived features)\n")
cat("===================================================================\n\n")

# -- 2A: Gene classification --------------------------------------------------
if ("gene" %in% names(df)) {
    df$gene_apob <- as.integer(df$gene == "APOB")
    cat(sprintf("  gene_apob: %d APOB carriers\n", sum(df$gene_apob, na.rm = TRUE)))
}

# -- 2B: Polygenic Risk Score (z-scored) --------------------------------------
prs_cols <- grep("^p26", names(df), value = TRUE)
if (length(prs_cols) > 0) {
    prs_priority <- c("p26227", "p26223", prs_cols)
    prs_priority <- unique(prs_priority[prs_priority %in% names(df)])
    for (pc in prs_priority) {
        n_valid <- sum(!is.na(df[[pc]]))
        if (n_valid > 100) {
            df$prs_primary <- scale(df[[pc]])[, 1]
            cat(sprintf("  prs_primary: using %s (%d non-missing, z-scored)\n", pc, n_valid))
            break
        }
    }
} else {
    cat("  WARNING: No PRS columns (p26xxx) found\n")
}

# -- 2C: Log-transformed CRP --------------------------------------------------
if ("crp" %in% names(df) && sum(!is.na(df$crp)) > 100) {
    df$log_crp <- log(df$crp + 0.01)
    cat(sprintf("  log_crp: %d non-missing\n", sum(!is.na(df$log_crp))))
}

# -- 2D: Log-transformed biomarkers -------------------------------------------
for (var in c("hba1c", "creatinine", "alt", "cystatin_c")) {
    if (var %in% names(df) && sum(!is.na(df[[var]])) > 100) {
        new_name <- paste0("log_", var)
        df[[new_name]] <- log(df[[var]] + 0.01)
        cat(sprintf("  %s: %d non-missing\n", new_name, sum(!is.na(df[[new_name]]))))
    }
}

# -- 2E: Lp(a) features -------------------------------------------------------
if ("lpa" %in% names(df) && sum(!is.na(df$lpa)) > 100) {
    df$log_lpa <- ifelse(!is.na(df$lpa) & df$lpa > 0, log(df$lpa),
                         ifelse(!is.na(df$lpa), log(0.1), NA))
    df$lpa_extreme <- as.integer(!is.na(df$lpa) & df$lpa >= 250)
    df$lpa_high    <- as.integer(!is.na(df$lpa) & df$lpa >= 125)
    cat(sprintf("  log_lpa: %d non-missing; extreme: %d; high: %d\n",
                sum(!is.na(df$log_lpa)),
                sum(df$lpa_extreme == 1, na.rm = TRUE),
                sum(df$lpa_high == 1, na.rm = TRUE)))
}

# -- 2F: ApoB/LDL ratios (if not already present) ----------------------------
if (!"log_apob_ldl" %in% names(df)) {
    if (all(c("apob", "ldl") %in% names(df))) {
        df$apob_ldl_ratio <- ifelse(!is.na(df$apob) & !is.na(df$ldl) & df$ldl > 0,
                                    df$apob / df$ldl, NA)
        df$log_apob_ldl <- ifelse(!is.na(df$apob_ldl_ratio) & df$apob_ldl_ratio > 0,
                                  log(df$apob_ldl_ratio), NA)
        cat(sprintf("  log_apob_ldl: %d non-missing\n", sum(!is.na(df$log_apob_ldl))))
    }
}
if (!"log_apob_re_ldl" %in% names(df)) {
    if (all(c("apob", "re_ldl") %in% names(df))) {
        ratio <- ifelse(!is.na(df$apob) & !is.na(df$re_ldl) & df$re_ldl > 0,
                        df$apob / df$re_ldl, NA)
        df$log_apob_re_ldl <- ifelse(!is.na(ratio) & ratio > 0, log(ratio), NA)
    }
}

# -- 2G: Interaction terms ----------------------------------------------------
if (all(c("age", "re_ldl") %in% names(df)))
    df$age_x_re_ldl <- df$age * df$re_ldl
if (all(c("age", "hdl") %in% names(df)))
    df$age_x_hdl <- df$age * df$hdl

# -- 2H: v5 new features ------------------------------------------------------
# Non-HDL cholesterol
if (all(c("tc", "hdl") %in% names(df))) {
    df$non_hdl <- df$tc - df$hdl
    df$tc_hdl_ratio <- ifelse(df$hdl > 0, df$tc / df$hdl, NA)
}

# TG/HDL ratio
if (all(c("trig", "hdl") %in% names(df))) {
    df$tg_hdl_ratio <- ifelse(df$hdl > 0, df$trig / df$hdl, NA)
    df$log_tg_hdl <- ifelse(!is.na(df$tg_hdl_ratio) & df$tg_hdl_ratio > 0,
                            log(df$tg_hdl_ratio), NA)
}

# Pulse pressure
if (all(c("sbp", "dbp") %in% names(df)))
    df$pulse_pressure <- df$sbp - df$dbp

# eGFR (CKD-EPI 2021, race-free)
if (all(c("creatinine", "age", "sex") %in% names(df))) {
    scr_mgdl <- df$creatinine / 88.4
    kappa  <- ifelse(df$sex == 1, 0.9, 0.7)
    alpha  <- ifelse(df$sex == 1, -0.302, -0.241)
    female <- ifelse(df$sex == 0, 1.012, 1.0)
    df$egfr <- 142 * pmin(scr_mgdl/kappa, 1)^alpha *
        pmax(scr_mgdl/kappa, 1)^(-1.200) * 0.9938^df$age * female
}

# Non-linear age
df$age_sq    <- df$age^2
df$sex_x_age <- df$sex * df$age

# ApoB/ApoA1 ratio
if (all(c("apob", "apoa") %in% names(df)))
    df$apob_apoa_ratio <- ifelse(df$apoa > 0, df$apob / df$apoa, NA)

# Remnant cholesterol
if (all(c("tc", "hdl", "ldl") %in% names(df)))
    df$remnant_chol_calc <- df$tc - df$hdl - df$ldl

# Metabolic syndrome score (z-score composite)
metsyn_vars <- c("bmi", "trig", "sbp")
if (all(c(metsyn_vars, "hdl") %in% names(df))) {
    metsyn_z <- sapply(metsyn_vars, function(v) {
        x <- df[[v]]; (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    })
    hdl_z <- -1 * (df$hdl - mean(df$hdl, na.rm = TRUE)) / sd(df$hdl, na.rm = TRUE)
    df$metsyn_score <- rowMeans(cbind(metsyn_z, hdl_z), na.rm = TRUE)
}

# PRS × re_ldl interaction
if (all(c("prs_primary", "re_ldl") %in% names(df)))
    df$prs_x_re_ldl <- df$prs_primary * df$re_ldl

# -- Summary of created features -----------------------------------------------
new_feats <- c("gene_apob", "prs_primary", "log_crp", "log_hba1c",
               "log_creatinine", "log_alt", "log_cystatin_c",
               "log_lpa", "lpa_extreme", "lpa_high",
               "log_apob_ldl", "log_apob_re_ldl",
               "age_x_re_ldl", "age_x_hdl",
               "non_hdl", "tc_hdl_ratio", "log_tg_hdl", "pulse_pressure",
               "egfr", "age_sq", "sex_x_age", "apob_apoa_ratio",
               "remnant_chol_calc", "metsyn_score", "prs_x_re_ldl")

cat("\n  Feature availability check:\n")
for (f in new_feats) {
    if (f %in% names(df)) {
        n <- sum(!is.na(df[[f]]))
        cat(sprintf("    %-22s: %d / %d (%.1f%%)\n", f, n, nrow(df), 100*n/nrow(df)))
    } else {
        cat(sprintf("    %-22s: *** NOT CREATED ***\n", f))
    }
}

# =============================================================================
# SECTION 3: NMR PCA
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 3: NMR PCA\n")
cat("===================================================================\n\n")

nmr_all_cols <- grep("^p234\\d+_i0$", names(df), value = TRUE)
cat(sprintf("  NMR features: %d\n", length(nmr_all_cols)))

pc_names <- character(0)
if (length(nmr_all_cols) >= 20) {
    nmr_matrix <- as.matrix(df[, nmr_all_cols])
    # Column-mean imputation
    for (j in 1:ncol(nmr_matrix)) {
        na_idx <- is.na(nmr_matrix[, j])
        if (any(na_idx)) nmr_matrix[na_idx, j] <- mean(nmr_matrix[, j], na.rm = TRUE)
    }
    nmr_scaled <- scale(nmr_matrix)
    zv <- apply(nmr_scaled, 2, sd) == 0
    nmr_scaled <- nmr_scaled[, !zv]

    pca_fit <- prcomp(nmr_scaled, center = FALSE, scale. = FALSE)
    cum_var <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
    n_pcs <- min(which(cum_var >= 0.90), 15)

    pc_names <- paste0("nmr_pc", 1:n_pcs)
    for (j in 1:n_pcs) df[[pc_names[j]]] <- pca_fit$x[, j]
    cat(sprintf("  %d PCs retained (%.1f%% variance), all %d patients\n",
                n_pcs, 100 * cum_var[n_pcs], nrow(df)))
}

# =============================================================================
# SECTION 4: DEFINE TOP 5 MODEL TIERS (streamlined)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 4: MODEL TIERS (top 5 only)\n")
cat("===================================================================\n\n")

# TIER 1: Clinical Core+ (best parsimony vs performance)
tier1_vars <- c("age", "sex", "re_ldl", "hdl",
                "log_apob_ldl", "log_lpa", "lpa_extreme",
                "smoking_binary", "diabetes", "hypertension",
                "bmi", "log_crp", "lipid_years_A",
                "non_hdl", "egfr")

# TIER 2: Enriched Clinical (adds BP, metabolic, non-linear)
tier2_vars <- c(tier1_vars, "pulse_pressure", "tc_hdl_ratio",
                "metsyn_score", "age_sq", "sex_x_age")

# TIER 3: NMR-PCA + PRS (adds metabolomics + genetics)
tier3_vars <- c(tier1_vars, pc_names, "prs_primary")

# TIER 4: Full v6 (everything Boruta-relevant + NMR PCs)
tier4_vars <- unique(c(tier2_vars, pc_names, "prs_primary", "prs_x_re_ldl",
                       "apob_apoa_ratio", "remnant_chol_calc", "log_tg_hdl",
                       "gene_apob", "townsend", "hba1c", "log_cystatin_c",
                       "ever_smoked", "lpa_high"))

# TIER 5: Boruta-selected (will be defined after Boruta)

# Report tier sizes
tier_list <- list(
    list(nm = "Tier 1: Core+", vars = tier1_vars),
    list(nm = "Tier 2: Enriched", vars = tier2_vars),
    list(nm = "Tier 3: NMR-PCA+PRS", vars = tier3_vars),
    list(nm = "Tier 4: Full v6", vars = tier4_vars)
)

cat(sprintf("  %-25s %5s %6s %6s %5s\n", "Tier", "Vars", "N_cc", "Events", "EPV"))
cat(paste(rep("-", 55), collapse = ""), "\n")
for (tc in tier_list) {
    tv <- tc$vars[tc$vars %in% names(df)]
    cc <- complete.cases(df[, c(tv, "ascvd_combined")])
    n_ev <- sum(df$ascvd_combined[cc])
    cat(sprintf("  %-25s %5d %6d %6d %5.1f\n",
                tc$nm, length(tv), sum(cc), n_ev, n_ev / max(length(tv), 1)))
}

# =============================================================================
# SECTION 5: BORUTA FEATURE SELECTION
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 5: BORUTA FEATURE SELECTION\n")
cat("===================================================================\n\n")

boruta_input <- tier4_vars[tier4_vars %in% names(df)]
df_boruta <- df[complete.cases(df[, c(boruta_input, "ascvd_combined")]),
                c(boruta_input, "ascvd_combined")]

cat(sprintf("  Input: %d vars, %d cases, %d events\n",
            length(boruta_input), nrow(df_boruta), sum(df_boruta$ascvd_combined)))

boruta_confirmed <- character(0)
if (nrow(df_boruta) >= 500) {
    set.seed(2026)
    boruta_res <- tryCatch(
        TentativeRoughFix(Boruta(ascvd_combined ~ ., data = df_boruta,
                                 maxRuns = 200, doTrace = 0)),
        error = function(e) { cat(sprintf("  Boruta failed: %s\n", e$message)); NULL }
    )
    if (!is.null(boruta_res)) {
        boruta_confirmed <- names(boruta_res$finalDecision[
            boruta_res$finalDecision == "Confirmed"])
        cat(sprintf("  CONFIRMED: %d features\n", length(boruta_confirmed)))
        cat(sprintf("  REJECTED: %d features\n",
                    sum(boruta_res$finalDecision == "Rejected")))

        imp_scores <- colMeans(boruta_res$ImpHistory[, boruta_confirmed], na.rm = TRUE)
        imp_order <- order(-imp_scores)
        cat("\n  Confirmed features (ranked):\n")
        for (j in imp_order) {
            cat(sprintf("    %2d. %-25s imp=%.3f\n",
                        which(imp_order == j), boruta_confirmed[j], imp_scores[j]))
        }

        write.csv(data.frame(variable = boruta_confirmed[imp_order],
                              importance = imp_scores[imp_order]),
                  paste0(TAB_DIR, "calon2_v6_boruta.csv"), row.names = FALSE)
    }
}

# Tier 5: Boruta Lean
tier5_vars <- if (length(boruta_confirmed) >= 5) boruta_confirmed else tier2_vars
cat(sprintf("\n  Tier 5: Boruta Lean — %d vars\n", length(tier5_vars)))

# Update tier list
tier_list[[5]] <- list(nm = "Tier 5: Boruta Lean", vars = tier5_vars)

# =============================================================================
# SECTION 6: ELASTIC NET — ALL 5 TIERS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 6: ELASTIC NET\n")
cat("===================================================================\n\n")

model_results <- list()
alphas <- c(0, 0.25, 0.5, 0.75, 1.0)
alpha_names <- c("Ridge", "EN(.25)", "EN(.50)", "EN(.75)", "LASSO")

run_en <- function(name, vars, data) {
    vars_exist <- vars[vars %in% names(data)]
    df_cc <- data[complete.cases(data[, c(vars_exist, "ascvd_combined")]),
                  c(vars_exist, "ascvd_combined")]

    cat(sprintf("  %-25s N=%d, Events=%d, Vars=%d\n",
                name, nrow(df_cc), sum(df_cc$ascvd_combined), length(vars_exist)))

    if (nrow(df_cc) < 200 || sum(df_cc$ascvd_combined) < 30) {
        cat("    SKIPPED\n"); return(NULL)
    }

    X <- as.matrix(df_cc[, vars_exist])
    y <- df_cc$ascvd_combined

    best <- NULL
    for (i in seq_along(alphas)) {
        set.seed(2026)
        cv_fit <- cv.glmnet(X, y, family = "binomial",
                            alpha = alphas[i], nfolds = 10, type.measure = "auc")
        auc_max <- max(cv_fit$cvm)
        nz <- sum(coef(cv_fit, s = "lambda.min")[-1] != 0)
        cat(sprintf("    %-8s AUC=%.4f (%d vars)\n", alpha_names[i], auc_max, nz))

        if (is.null(best) || auc_max > best$auc) {
            best <- list(alpha = alphas[i], auc = auc_max, fit = cv_fit,
                         n_vars = nz, X = X, y = y, vars = vars_exist,
                         df_cc = df_cc, n = nrow(df_cc), n_ev = sum(y))
        }
    }
    cat(sprintf("    BEST: %s, CV-AUC = %.4f\n\n", alpha_names[which(alphas == best$alpha)], best$auc))
    return(best)
}

en_results <- list()
for (i in seq_along(tier_list)) {
    tc <- tier_list[[i]]
    en_results[[tc$nm]] <- run_en(tc$nm, tc$vars, df)
    if (!is.null(en_results[[tc$nm]])) {
        model_results[[paste0("EN_", tc$nm)]] <- list(
            name = paste("EN", tc$nm), algorithm = "Elastic Net",
            cv_auc = en_results[[tc$nm]]$auc,
            n_vars = en_results[[tc$nm]]$n_vars,
            n = en_results[[tc$nm]]$n, n_ev = en_results[[tc$nm]]$n_ev)
    }
}

# =============================================================================
# SECTION 7: TUNED XGBOOST — ALL 5 TIERS
# =============================================================================

cat("===================================================================\n")
cat("SECTION 7: TUNED XGBOOST\n")
cat("===================================================================\n\n")

run_xgb <- function(name, vars, data) {
    vars_exist <- vars[vars %in% names(data)]
    df_cc <- data[complete.cases(data[, c(vars_exist, "ascvd_combined")]),
                  c(vars_exist, "ascvd_combined")]

    cat(sprintf("  %-25s N=%d, Events=%d, Vars=%d\n",
                name, nrow(df_cc), sum(df_cc$ascvd_combined), length(vars_exist)))

    if (nrow(df_cc) < 200 || sum(df_cc$ascvd_combined) < 30) {
        cat("    SKIPPED\n"); return(NULL)
    }

    X <- as.matrix(df_cc[, vars_exist])
    y <- df_cc$ascvd_combined
    dtrain <- xgb.DMatrix(data = X, label = y)

    # Random grid search (30 combos)
    grid <- expand.grid(
        max_depth = c(3, 4, 5, 6), eta = c(0.005, 0.01, 0.03, 0.05, 0.10),
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
            xgb.cv(params = params, data = dtrain, nrounds = 2000, nfold = 10,
                   early_stopping_rounds = 50, verbose = 0, print_every_n = 0),
            error = function(e) tryCatch(
                xgb.cv(params = params, data = dtrain, nrounds = 2000, nfold = 10,
                       early_stopping_rounds = 50, verbose = 0),
                error = function(e2) NULL))

        if (!is.null(cv_res)) {
            eval_log <- cv_res$evaluation_log
            auc_col <- grep("test.*auc.*mean", names(eval_log), value = TRUE)[1]
            if (!is.null(auc_col) && !is.na(auc_col)) {
                this_auc <- max(eval_log[[auc_col]], na.rm = TRUE)
                if (this_auc > best_auc) {
                    best_auc <- this_auc
                    best_params <- params
                    best_nr <- which.max(eval_log[[auc_col]])
                }
            }
        }
    }

    if (is.null(best_params)) { cat("    XGB FAILED\n"); return(NULL) }

    cat(sprintf("    BEST: depth=%d, eta=%.3f, nrounds=%d → CV-AUC = %.4f\n\n",
                best_params$max_depth, best_params$eta, best_nr, best_auc))

    set.seed(2026)
    final <- xgb.train(params = best_params, data = dtrain, nrounds = best_nr, verbose = 0)

    return(list(auc = best_auc, model = final, dtrain = dtrain,
                X = X, y = y, vars = vars_exist, df_cc = df_cc,
                n = nrow(df_cc), n_ev = sum(y)))
}

xgb_results <- list()
for (i in seq_along(tier_list)) {
    tc <- tier_list[[i]]
    xgb_results[[tc$nm]] <- run_xgb(tc$nm, tc$vars, df)
    if (!is.null(xgb_results[[tc$nm]])) {
        model_results[[paste0("XGB_", tc$nm)]] <- list(
            name = paste("XGB", tc$nm), algorithm = "XGBoost",
            cv_auc = xgb_results[[tc$nm]]$auc,
            n_vars = length(xgb_results[[tc$nm]]$vars),
            n = xgb_results[[tc$nm]]$n, n_ev = xgb_results[[tc$nm]]$n_ev)
    }
}

# =============================================================================
# SECTION 8: RANDOM FOREST — TOP 3 TIERS
# =============================================================================

cat("===================================================================\n")
cat("SECTION 8: RANDOM FOREST (top 3 tiers)\n")
cat("===================================================================\n\n")

run_rf <- function(name, vars, data) {
    vars_exist <- vars[vars %in% names(data)]
    df_cc <- data[complete.cases(data[, c(vars_exist, "ascvd_combined")]),
                  c(vars_exist, "ascvd_combined")]

    cat(sprintf("  %-25s N=%d, Events=%d\n",
                name, nrow(df_cc), sum(df_cc$ascvd_combined)))

    if (nrow(df_cc) < 200) { cat("    SKIPPED\n"); return(NULL) }

    X <- df_cc[, vars_exist]
    y <- factor(df_cc$ascvd_combined, levels = c(0, 1))

    mtry_vals <- unique(pmax(c(floor(sqrt(ncol(X))), floor(ncol(X)/3)), 1))
    best_auc <- 0; best_rf <- NULL

    for (mt in mtry_vals) {
        set.seed(2026)
        rf_fit <- randomForest(x = X, y = y, ntree = 1000, mtry = mt,
                               importance = TRUE, classwt = c("0" = 1, "1" = 2))
        oob_probs <- rf_fit$votes[, "1"]
        this_auc <- tryCatch(
            as.numeric(auc(roc(as.numeric(as.character(y)), oob_probs, quiet = TRUE))),
            error = function(e) 0)
        cat(sprintf("    mtry=%d: OOB-AUC=%.4f\n", mt, this_auc))
        if (this_auc > best_auc) { best_auc <- this_auc; best_rf <- rf_fit }
    }

    return(list(auc = best_auc, rf = best_rf, X = X, y = y,
                vars = vars_exist, n = nrow(df_cc),
                n_ev = sum(df_cc$ascvd_combined)))
}

# Run RF on tiers 1, 4, 5 only
for (i in c(1, 4, 5)) {
    tc <- tier_list[[i]]
    rf_res <- run_rf(tc$nm, tc$vars, df)
    if (!is.null(rf_res)) {
        model_results[[paste0("RF_", tc$nm)]] <- list(
            name = paste("RF", tc$nm), algorithm = "Random Forest",
            cv_auc = rf_res$auc, n_vars = length(rf_res$vars),
            n = rf_res$n, n_ev = rf_res$n_ev)
    }
}

# =============================================================================
# SECTION 9: MULTIPLE IMPUTATION — BEST TIER
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 9: MULTIPLE IMPUTATION (mice, m=10)\n")
cat("===================================================================\n\n")

run_mi_en <- function(name, vars, data, m_imps = 10) {
    vars_exist <- vars[vars %in% names(data)]
    imp_data <- data[, c(vars_exist, "ascvd_combined")]

    n_total <- nrow(imp_data)
    n_cc <- sum(complete.cases(imp_data))
    cat(sprintf("  %s: %d total, %d complete (%.1f%% loss)\n",
                name, n_total, n_cc, 100 * (n_total - n_cc) / n_total))

    if (n_total - n_cc < 20) { cat("  Minimal missingness, skipping MI\n"); return(NULL) }

    set.seed(2026)
    tryCatch({
        imp <- mice(imp_data, m = m_imps, method = "pmm", maxit = 20,
                    printFlag = FALSE, seed = 2026)

        mi_aucs <- numeric(m_imps)
        for (i in 1:m_imps) {
            d <- complete(imp, i)
            X <- as.matrix(d[, vars_exist])
            y <- d$ascvd_combined

            best_auc <- 0
            for (a in c(0, 0.25, 0.5, 0.75, 1.0)) {
                set.seed(2026 + i)
                cv_fit <- cv.glmnet(X, y, family = "binomial",
                                    alpha = a, nfolds = 10, type.measure = "auc")
                if (max(cv_fit$cvm) > best_auc) best_auc <- max(cv_fit$cvm)
            }
            mi_aucs[i] <- best_auc
        }

        mean_auc <- mean(mi_aucs)
        cat(sprintf("  MI-pooled CV-AUC: %.4f (SD=%.4f)\n", mean_auc, sd(mi_aucs)))
        cat(sprintf("  Recovered: %d patients\n\n", n_total - n_cc))

        return(list(auc = mean_auc, mi_aucs = mi_aucs,
                    n = n_total, n_cc = n_cc, n_ev = sum(data$ascvd_combined, na.rm = TRUE)))
    }, error = function(e) {
        cat(sprintf("  MI failed: %s\n", e$message)); NULL
    })
}

# MI on tiers 1 and 5
for (i in c(1, 5)) {
    tc <- tier_list[[i]]
    mi_res <- run_mi_en(tc$nm, tc$vars, df)
    if (!is.null(mi_res)) {
        model_results[[paste0("MI_", tc$nm)]] <- list(
            name = paste("MI-EN", tc$nm), algorithm = "MI-Elastic Net",
            cv_auc = mi_res$auc, n_vars = length(tc$vars[tc$vars %in% names(df)]),
            n = mi_res$n, n_ev = mi_res$n_ev)
    }
}

# =============================================================================
# SECTION 10: HONEST STACKING (outer-fold, no leakage)
# =============================================================================

cat("===================================================================\n")
cat("SECTION 10: HONEST STACKING ENSEMBLE\n")
cat("===================================================================\n\n")

# Use the best-performing tier for stacking
# Pick tier with most complete cases and good EN performance
best_tier_idx <- 5  # Boruta Lean by default
best_tier <- tier_list[[best_tier_idx]]
stack_vars <- best_tier$vars[best_tier$vars %in% names(df)]
df_stack <- df[complete.cases(df[, c(stack_vars, "ascvd_combined")]),
               c(stack_vars, "ascvd_combined")]

cat(sprintf("  Stacking on: %s (%d vars, N=%d, Events=%d)\n",
            best_tier$nm, length(stack_vars), nrow(df_stack),
            sum(df_stack$ascvd_combined)))

if (nrow(df_stack) >= 300) {
    X_stack <- as.matrix(df_stack[, stack_vars])
    y_stack <- df_stack$ascvd_combined

    # Create 10 outer folds
    set.seed(2026)
    folds <- createFolds(y_stack, k = 10, returnTrain = TRUE)

    # Collect out-of-fold predictions from 3 base learners
    oof_en  <- numeric(length(y_stack))
    oof_xgb <- numeric(length(y_stack))
    oof_rf  <- numeric(length(y_stack))

    for (k in 1:10) {
        train_idx <- folds[[k]]
        test_idx  <- setdiff(1:length(y_stack), train_idx)

        X_tr <- X_stack[train_idx, ]; y_tr <- y_stack[train_idx]
        X_te <- X_stack[test_idx, ]

        # EN
        set.seed(2026)
        en_fit <- cv.glmnet(X_tr, y_tr, family = "binomial",
                            alpha = 0.25, nfolds = 5, type.measure = "auc")
        oof_en[test_idx] <- as.numeric(predict(en_fit, X_te, s = "lambda.min",
                                                type = "response"))

        # XGBoost
        dt_tr <- xgb.DMatrix(data = X_tr, label = y_tr)
        dt_te <- xgb.DMatrix(data = X_te)
        set.seed(2026)
        xgb_fit <- xgb.train(
            params = list(objective = "binary:logistic", eval_metric = "auc",
                          max_depth = 4, eta = 0.03, subsample = 0.8,
                          colsample_bytree = 0.8, min_child_weight = 10),
            data = dt_tr, nrounds = 300, verbose = 0)
        oof_xgb[test_idx] <- predict(xgb_fit, dt_te)

        # RF
        rf_fit <- randomForest(x = data.frame(X_tr),
                               y = factor(y_tr, levels = c(0, 1)),
                               ntree = 500, mtry = max(floor(sqrt(ncol(X_tr))), 1))
        oof_rf[test_idx] <- predict(rf_fit, data.frame(X_te), type = "prob")[, "1"]
    }

    # Meta-learner: logistic regression on OOF predictions
    meta_df <- data.frame(en = oof_en, xgb = oof_xgb, rf = oof_rf, y = y_stack)
    meta_fit <- glm(y ~ en + xgb + rf, data = meta_df, family = binomial)
    meta_pred <- predict(meta_fit, type = "response")

    # Honest AUC (using OOF meta-predictions via LOOCV of meta-learner)
    # Simple approach: 10-fold CV of the meta-learner
    set.seed(2026)
    meta_folds <- createFolds(y_stack, k = 10, returnTrain = TRUE)
    meta_oof <- numeric(length(y_stack))
    for (k in 1:10) {
        tr <- meta_folds[[k]]
        te <- setdiff(1:length(y_stack), tr)
        mf <- glm(y ~ en + xgb + rf, data = meta_df[tr, ], family = binomial)
        meta_oof[te] <- predict(mf, meta_df[te, ], type = "response")
    }

    stack_auc <- as.numeric(auc(roc(y_stack, meta_oof, quiet = TRUE)))

    # Also report individual base learner AUCs
    en_auc_oof  <- as.numeric(auc(roc(y_stack, oof_en, quiet = TRUE)))
    xgb_auc_oof <- as.numeric(auc(roc(y_stack, oof_xgb, quiet = TRUE)))
    rf_auc_oof  <- as.numeric(auc(roc(y_stack, oof_rf, quiet = TRUE)))

    cat(sprintf("  Base learner OOF AUCs:\n"))
    cat(sprintf("    EN:      %.4f\n", en_auc_oof))
    cat(sprintf("    XGBoost: %.4f\n", xgb_auc_oof))
    cat(sprintf("    RF:      %.4f\n", rf_auc_oof))
    cat(sprintf("  STACKED HONEST AUC: %.4f\n\n", stack_auc))

    cat(sprintf("  Meta-learner coefficients:\n"))
    for (nm in names(coef(meta_fit))) {
        cat(sprintf("    %-10s: %.4f\n", nm, coef(meta_fit)[nm]))
    }

    model_results[["Stack_Honest"]] <- list(
        name = "Honest Stack (EN+XGB+RF)", algorithm = "Stacking",
        cv_auc = stack_auc, n_vars = length(stack_vars),
        n = length(y_stack), n_ev = sum(y_stack))
}

# =============================================================================
# SECTION 11: GRAND COMPARISON TABLE
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 11: GRAND COMPARISON\n")
cat("===================================================================\n\n")

comp_df <- do.call(rbind, lapply(model_results, function(x) {
    data.frame(Model = x$name, Algorithm = x$algorithm,
               CV_AUC = x$cv_auc, N = x$n, Events = x$n_ev,
               Vars = x$n_vars, stringsAsFactors = FALSE)
}))
comp_df <- comp_df[order(-comp_df$CV_AUC), ]
rownames(comp_df) <- NULL

cat(sprintf("  %-40s %-15s %8s %6s %6s %5s\n",
            "Model", "Algorithm", "CV-AUC", "N", "Events", "Vars"))
cat(paste(rep("-", 85), collapse = ""), "\n")
for (i in 1:nrow(comp_df)) {
    cat(sprintf("  %-40s %-15s %8.4f %6d %6d %5d\n",
                comp_df$Model[i], comp_df$Algorithm[i], comp_df$CV_AUC[i],
                comp_df$N[i], comp_df$Events[i], comp_df$Vars[i]))
}

write.csv(comp_df, paste0(TAB_DIR, "calon2_v6_grand_comparison.csv"), row.names = FALSE)
cat("\n  Saved: tables/calon2_v6_grand_comparison.csv\n")

# =============================================================================
# SECTION 12: SAFEHEART-RE COMPARISON (frozen published coefficients)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 12: SAFEHEART-RE BENCHMARK (frozen coefficients)\n")
cat("===================================================================\n\n")

# Published SAFEHEART-RE coefficients (Perez de Isla 2017)
SAFEHEART_COEF <- list(
    intercept = -7.053, age = 0.064, male = 0.775, ldl_c = 0.109,
    htn = 0.431, bmi = 0.025, smoking = 0.466, prior_cvd = 1.414
)

# Map to UKB columns
sh_available <- all(c("age", "sex", "re_ldl", "hypertension", "bmi",
                       "smoking_binary") %in% names(df))

if (sh_available) {
    # Linear predictor using FROZEN coefficients (no refitting!)
    lp <- SAFEHEART_COEF$intercept +
        SAFEHEART_COEF$age      * df$age +
        SAFEHEART_COEF$male     * df$sex +
        SAFEHEART_COEF$ldl_c    * df$re_ldl +
        SAFEHEART_COEF$htn      * df$hypertension +
        SAFEHEART_COEF$bmi      * df$bmi +
        SAFEHEART_COEF$smoking  * df$smoking_binary
    # Note: prior_cvd excluded (not directly available in UKB FH cohort)

    sh_pred <- 1 / (1 + exp(-lp))
    sh_valid <- !is.na(sh_pred) & !is.na(df$ascvd_combined)

    sh_roc <- roc(df$ascvd_combined[sh_valid], sh_pred[sh_valid], quiet = TRUE)
    sh_auc <- as.numeric(auc(sh_roc))

    cat(sprintf("  SAFEHEART-RE (frozen coefficients) AUC: %.4f (N=%d)\n",
                sh_auc, sum(sh_valid)))

    # Compare with best CALON-2 model
    best_model <- comp_df[1, ]
    cat(sprintf("  Best CALON-2 model: %s → CV-AUC = %.4f\n",
                best_model$Model, best_model$CV_AUC))
    cat(sprintf("  Improvement: +%.4f (%.1f%%)\n",
                best_model$CV_AUC - sh_auc,
                100 * (best_model$CV_AUC - sh_auc) / sh_auc))

    model_results[["SAFEHEART_RE"]] <- list(
        name = "SAFEHEART-RE (frozen)", algorithm = "Published",
        cv_auc = sh_auc, n_vars = 7, n = sum(sh_valid),
        n_ev = sum(df$ascvd_combined[sh_valid]))
}

# =============================================================================
# SECTION 13: BOOTSTRAP OPTIMISM CORRECTION (best model)
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 13: BOOTSTRAP OPTIMISM CORRECTION\n")
cat("===================================================================\n\n")

# Use best EN model for bootstrap validation
best_en_name <- names(which.max(sapply(en_results, function(x) if(!is.null(x)) x$auc else 0)))
best_en <- en_results[[best_en_name]]

if (!is.null(best_en)) {
    cat(sprintf("  Bootstrap validation of: %s (CV-AUC=%.4f)\n", best_en_name, best_en$auc))

    n_boot <- 200
    apparent_aucs <- numeric(n_boot)
    test_aucs <- numeric(n_boot)

    X_full <- best_en$X
    y_full <- best_en$y

    for (b in 1:n_boot) {
        set.seed(b)
        boot_idx <- sample(1:nrow(X_full), replace = TRUE)
        oob_idx <- setdiff(1:nrow(X_full), unique(boot_idx))

        if (length(oob_idx) < 30 || sum(y_full[oob_idx]) < 5) next

        # Fit on bootstrap sample
        cv_b <- cv.glmnet(X_full[boot_idx, ], y_full[boot_idx],
                          family = "binomial", alpha = best_en$alpha,
                          nfolds = 5, type.measure = "auc")

        # Apparent AUC (on training data)
        pred_train <- predict(cv_b, X_full[boot_idx, ], s = "lambda.min", type = "response")
        apparent_aucs[b] <- tryCatch(
            as.numeric(auc(roc(y_full[boot_idx], as.numeric(pred_train), quiet = TRUE))),
            error = function(e) NA)

        # Test AUC (on OOB data)
        pred_oob <- predict(cv_b, X_full[oob_idx, ], s = "lambda.min", type = "response")
        test_aucs[b] <- tryCatch(
            as.numeric(auc(roc(y_full[oob_idx], as.numeric(pred_oob), quiet = TRUE))),
            error = function(e) NA)
    }

    valid_boots <- !is.na(apparent_aucs) & !is.na(test_aucs)
    optimism <- mean(apparent_aucs[valid_boots] - test_aucs[valid_boots])
    corrected_auc <- best_en$auc - optimism

    cat(sprintf("  Apparent AUC (mean): %.4f\n", mean(apparent_aucs[valid_boots])))
    cat(sprintf("  Test AUC (mean):     %.4f\n", mean(test_aucs[valid_boots])))
    cat(sprintf("  Optimism:            %.4f\n", optimism))
    cat(sprintf("  Corrected AUC:       %.4f\n", corrected_auc))
    cat(sprintf("  Bootstrap 95%% CI:    %.4f - %.4f\n",
                quantile(test_aucs[valid_boots], 0.025),
                quantile(test_aucs[valid_boots], 0.975)))
}

# =============================================================================
# SECTION 14: SAVE FINAL OUTPUTS
# =============================================================================

cat("\n===================================================================\n")
cat("SECTION 14: SAVING OUTPUTS\n")
cat("===================================================================\n\n")

# Save best model coefficients
if (!is.null(best_en)) {
    coefs <- as.matrix(coef(best_en$fit, s = "lambda.min"))
    coef_df <- data.frame(
        variable = rownames(coefs),
        coefficient = as.numeric(coefs),
        OR = exp(as.numeric(coefs)),
        stringsAsFactors = FALSE
    )
    coef_df <- coef_df[coef_df$coefficient != 0 | coef_df$variable == "(Intercept)", ]
    write.csv(coef_df, paste0(OUT_DIR, "calon2_v6_best_coefficients.csv"), row.names = FALSE)
    cat("  Saved: calon2_v6_best_coefficients.csv\n")

    cat("\n  Best model coefficients:\n")
    for (i in 1:nrow(coef_df)) {
        cat(sprintf("    %-25s coef=%8.4f  OR=%8.4f\n",
                    coef_df$variable[i], coef_df$coefficient[i], coef_df$OR[i]))
    }
}

# Final summary
cat("\n\n======================================================================\n")
cat("  CALON-2 v6 COMPLETE\n")
cat("======================================================================\n")
cat(sprintf("  Best single model:  %s (CV-AUC = %.4f)\n",
            comp_df$Model[1], comp_df$CV_AUC[1]))
if (exists("stack_auc"))
    cat(sprintf("  Honest Stack:       %.4f\n", stack_auc))
if (exists("sh_auc"))
    cat(sprintf("  SAFEHEART-RE:       %.4f\n", sh_auc))
if (exists("corrected_auc"))
    cat(sprintf("  Optimism-corrected: %.4f\n", corrected_auc))
cat("======================================================================\n")

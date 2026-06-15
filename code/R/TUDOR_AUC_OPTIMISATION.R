################################################################################
#
#  TUDOR FH DIAGNOSTIC MODEL: AUC OPTIMISATION — TARGET > 0.82 ALL COHORTS
#
#  DIAGNOSTIC model: predict genetically confirmed FH from clinical phenotype
#  Outcome: fh (1 = mutation positive, 0 = mutation negative / control)
#
#  Three strategies tested systematically:
#    A. Feature engineering — expanded feature sets with domain-guided selection
#    B. Ensemble stacking — ENET + RF + XGB with calibrated meta-learner
#    C. Pooled backbone + cohort-specific recalibration
#
#  Cohorts:
#    1. South Wales — mutation positive (n~425) vs negative controls (n~941)
#    2. Wales PASS — Positive1=Yes vs No (all-Wales FH registry)
#    3. UKB Lipid Clinic — WES-confirmed FH vs lipid-clinic controls
#
#  Validation: Leave-One-Cohort-Out (LOCO-CV)
#  Reporting: TRIPOD Type 3
#
#  Author: Dr Nader Genedy
#  Date: 2026-03-27
#
################################################################################

cat("\n================================================================\n")
cat("  TUDOR DIAGNOSTIC: AUC OPTIMISATION\n")
cat("  Outcome: Genetically confirmed FH (diagnostic)\n")
cat("  Target: AUC > 0.82 in all 3 cohorts\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

set.seed(42)

# ── 0. PACKAGES ──────────────────────────────────────────────────────────────

pkgs <- c("readxl", "haven", "pROC", "glmnet", "randomForest",
          "xgboost", "data.table")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(readxl); library(haven); library(pROC); library(glmnet)
  library(randomForest); library(xgboost); library(data.table)
})

# ── 1. PATHS & HELPERS ──────────────────────────────────────────────────────

DD  <- "C:/Users/nader/Downloads/calon_ukb_pipeline"
TD  <- "D:/alphafold_backup/tudor_packup"
OUT <- file.path(DD, "tudor_loco_output")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

statin_factors <- c(atorvastatin = 0.38, simvastatin = 0.35, rosuvastatin = 0.34,
                    pravastatin = 0.25, fluvastatin = 0.22, lovastatin = 0.25,
                    ezetimibe = 0.18)

get_red <- function(tx) {
  if (is.na(tx) || tx == "" || tx == "NaN") return(0)
  tx_low <- tolower(as.character(tx))
  for (nm in names(statin_factors))
    if (grepl(nm, tx_low, fixed = TRUE)) return(statin_factors[[nm]])
  if (grepl("statin", tx_low)) return(0.30)
  return(0)
}

rev_ldl <- function(ldl, rf) ifelse(!is.na(ldl) & !is.na(rf) & rf > 0, ldl / (1 - rf), ldl)

parse_dt <- function(x) {
  x <- as.character(x)
  d <- suppressWarnings(as.Date(x, format = "%d-%m-%Y"))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, format = "%d/%m/%Y"))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  d
}

# ── 2. LOAD COHORT 1: SOUTH WALES (mutation+ vs mutation-) ──────────────────

cat("=== Loading cohorts ===\n")

pos <- as.data.frame(read_excel(file.path(DD, "mutation positive group.xlsx")))
pos$fh <- 1
pos$dob_d  <- as.Date(sapply(pos$BirthDate, parse_dt), origin = "1970-01-01")
pos$meas_d <- as.Date(sapply(pos$MeasurementDate.1, parse_dt), origin = "1970-01-01")
pos$age    <- as.numeric(difftime(pos$meas_d, pos$dob_d, units = "days")) / 365.25
pos$sex    <- ifelse(pos$Gender %in% c("M", "Male"), 1, 0)
pos$hdl    <- safe_num(pos$HDL.1); pos$tg <- safe_num(pos$TRG.1)
pos$ldl_m  <- safe_num(pos$LDL.1); pos$tc <- safe_num(pos$TC.1)
pos$rf     <- sapply(pos$Treatment1.1, get_red)
pos$ldl_ut <- ifelse(pos$rf > 0, rev_ldl(pos$ldl_m, pos$rf), pos$ldl_m)
pos$trig_filter <- ifelse(!is.na(pos$ldl_ut) & !is.na(pos$tg) & pos$tg >= 0,
                           pmin(pos$ldl_ut / (pos$tg + 0.1), 50), NA)
pos$on_statin <- ifelse(pos$rf > 0, 1, 0)
pos$tendon_xanth <- ifelse(pos$TendonXanthomata %in% c(1, "1", "Yes", "TRUE"), 1, 0)
pos$tendon_xanth[is.na(pos$tendon_xanth)] <- 0
ca1 <- if ("Corneal Arcus" %in% names(pos)) "Corneal Arcus" else "CornealArcus"
pos$corneal_arcus <- ifelse(pos[[ca1]] %in% c(1, "1", "Yes", "TRUE"), 1, 0)
pos$corneal_arcus[is.na(pos$corneal_arcus)] <- 0
pos$dlcn <- safe_num(pos$GenoTypingScore)
apob_col_pos <- intersect(c("Apo-B", "ApoB", "Apo.B", "apoB"), names(pos))[1]
pos$apob <- if (!is.na(apob_col_pos)) safe_num(pos[[apob_col_pos]]) else NA
pos$cohort <- "SouthWales"

neg <- as.data.frame(read_excel(file.path(DD, "Mutation negative control updated.xlsx")))
neg$fh  <- 0
neg$age <- safe_num(neg[["age at the result"]])
neg$sex <- ifelse(neg$Gender %in% c("M", "Male"), 1, 0)
neg$hdl <- safe_num(neg$HDL.1); neg$tg <- safe_num(neg$TRG.1)
neg$ldl_m <- safe_num(neg$LDL.1); neg$tc <- safe_num(neg$TC.1)
neg$rf  <- sapply(neg$Treatment1.1, get_red)
neg$ldl_ut <- ifelse(neg$rf > 0, rev_ldl(neg$ldl_m, neg$rf), neg$ldl_m)
neg$trig_filter <- ifelse(!is.na(neg$ldl_ut) & !is.na(neg$tg) & neg$tg >= 0,
                           pmin(neg$ldl_ut / (neg$tg + 0.1), 50), NA)
neg$on_statin <- ifelse(neg$rf > 0, 1, 0)
neg$tendon_xanth <- ifelse(neg$TendonXanthomata %in% c(1, "1", "Yes", "TRUE"), 1, 0)
neg$tendon_xanth[is.na(neg$tendon_xanth)] <- 0
ca2 <- if ("Corneal Arcus" %in% names(neg)) "Corneal Arcus" else "CornealArcus"
neg$corneal_arcus <- ifelse(neg[[ca2]] %in% c(1, "1", "Yes", "TRUE"), 1, 0)
neg$corneal_arcus[is.na(neg$corneal_arcus)] <- 0
neg$dlcn <- safe_num(neg$GenoTypingScore)
apob_col_neg <- intersect(c("apoB...2", "apoB", "ApoB"), names(neg))[1]
neg$apob <- if (!is.na(apob_col_neg)) safe_num(neg[[apob_col_neg]]) else NA
neg$cohort <- "SouthWales"

# Engineered diagnostic features
compute_diag_features <- function(df) {
  df$non_hdl <- ifelse(!is.na(df$tc) & !is.na(df$hdl), df$tc - df$hdl, NA)
  df$nhdl_ldl_gap <- ifelse(!is.na(df$non_hdl) & !is.na(df$ldl_ut),
                             df$non_hdl - df$ldl_ut, NA)
  df$cumulative_ldl <- ifelse(!is.na(df$ldl_ut) & !is.na(df$age),
                               df$ldl_ut * df$age, NA)
  df$ldl_hdl_ratio <- ifelse(!is.na(df$ldl_ut) & !is.na(df$hdl) & df$hdl > 0,
                              df$ldl_ut / df$hdl, NA)
  df$tg_hdl_ratio <- ifelse(!is.na(df$tg) & !is.na(df$hdl) & df$hdl > 0,
                             df$tg / df$hdl, NA)
  df$apob_ldl <- ifelse(!is.na(df$apob) & !is.na(df$ldl_ut) & df$ldl_ut > 0,
                         df$apob / df$ldl_ut, NA)
  df$inv_hdl <- ifelse(!is.na(df$hdl) & df$hdl > 0, 1 / df$hdl, NA)
  df$log_trig_filter <- ifelse(!is.na(df$trig_filter) & df$trig_filter > 0,
                                log(df$trig_filter), NA)
  df$age_sq <- df$age^2
  df$has_phys_sign <- ifelse(df$tendon_xanth == 1 | df$corneal_arcus == 1, 1, 0)
  # Interaction: high LDL + low TG = classic FH pattern
  df$ldl_tg_contrast <- ifelse(!is.na(df$ldl_ut) & !is.na(df$tg),
                                df$ldl_ut - df$tg, NA)
  # Young + high LDL = more diagnostic of FH
  df$young_high_ldl <- ifelse(!is.na(df$age) & df$age < 50 &
                               !is.na(df$ldl_ut) & df$ldl_ut > 5, 1, 0)
  df
}

pos <- compute_diag_features(pos)
neg <- compute_diag_features(neg)

cols <- c("fh", "age", "sex", "hdl", "tg", "ldl_m", "ldl_ut", "trig_filter",
          "on_statin", "tendon_xanth", "corneal_arcus", "dlcn", "apob",
          "non_hdl", "nhdl_ldl_gap", "cumulative_ldl", "ldl_hdl_ratio",
          "tg_hdl_ratio", "apob_ldl", "inv_hdl", "log_trig_filter",
          "age_sq", "has_phys_sign", "ldl_tg_contrast", "young_high_ldl",
          "cohort")

sw <- rbind(pos[, cols], neg[, cols])
cat(sprintf("  SW: %d (FH+=%d, FH-=%d)\n", nrow(sw), sum(sw$fh), sum(sw$fh == 0)))

# ── 3. LOAD COHORT 2: WALES PASS ────────────────────────────────────────────

pass <- as.data.frame(read_sav(file.path(DD, "PASS_wrong_dob.sav")))
pass$fh  <- ifelse(pass$Positive1 == "Yes", 1, ifelse(pass$Positive1 == "No", 0, NA))
pass$age <- safe_num(pass$BMI_AGE)
pass$sex <- ifelse(pass$Gender == "M", 1, 0)
pass$hdl <- safe_num(pass$HDL.1); pass$tg <- safe_num(pass$TRG.1)
pass$ldl_m <- safe_num(pass$LDL.1); pass$tc <- safe_num(pass$TC.1)
pass$rf  <- sapply(pass$Treatment1.1, get_red)
pass$ldl_ut <- ifelse(pass$rf > 0, rev_ldl(pass$ldl_m, pass$rf), pass$ldl_m)
pass$trig_filter <- ifelse(!is.na(pass$ldl_ut) & !is.na(pass$tg) & pass$tg >= 0,
                            pmin(pass$ldl_ut / (pass$tg + 0.1), 50), NA)
pass$on_statin <- ifelse(pass$rf > 0, 1, 0)
pass$tendon_xanth <- ifelse(safe_num(pass$TendonXanthomata) == 1, 1, 0)
pass$tendon_xanth[is.na(pass$tendon_xanth)] <- 0
pass$corneal_arcus <- ifelse(safe_num(pass$CornealArcus) == 1, 1, 0)
pass$corneal_arcus[is.na(pass$corneal_arcus)] <- 0
pass$dlcn <- safe_num(pass$GenoTypingScore)
pass$apob <- NA
pass$cohort <- "Wales"

pass <- compute_diag_features(pass)
wales <- pass[!is.na(pass$fh), cols]
cat(sprintf("  Wales: %d (FH+=%d, FH-=%d)\n", nrow(wales), sum(wales$fh), sum(wales$fh == 0)))

# ── 4. LOAD COHORT 3: UKB LIPID CLINIC ──────────────────────────────────────

ukb_raw <- as.data.frame(fread(file.path(TD, "TUDOR_UKB_Features.csv"), showProgress = FALSE))
ukb_raw$fh <- safe_num(ukb_raw$is_fh_genetic); ukb_raw$fh[is.na(ukb_raw$fh)] <- 0
ukb_raw$age <- safe_num(ukb_raw$Age_at_LDL1)
ukb_raw$sex <- safe_num(ukb_raw$Gender_num)
ukb_raw$hdl <- safe_num(ukb_raw$HDL.1); ukb_raw$tg <- safe_num(ukb_raw$TRG.1)
ukb_raw$ldl_m <- safe_num(ukb_raw$LDL_treated)
ukb_raw$ldl_ut <- safe_num(ukb_raw$LDL_untreated)
ukb_raw$trig_filter <- safe_num(ukb_raw$Trig_Filter)
ukb_raw$on_statin <- ifelse(safe_num(ukb_raw$reduction_factor) > 0, 1, 0)
ukb_raw$on_statin[is.na(ukb_raw$on_statin)] <- 0
ukb_raw$tc <- safe_num(ukb_raw$CHOL)
ukb_raw$apob <- safe_num(ukb_raw$APOB)
ukb_raw$non_hdl <- ifelse(!is.na(ukb_raw$tc) & !is.na(ukb_raw$hdl),
                           ukb_raw$tc - ukb_raw$hdl, NA)

# Lipid clinic filter (NICE CG71)
ukb_raw$lc <- (!is.na(ukb_raw$tc) & ukb_raw$tc > 7.9) |
              (!is.na(ukb_raw$ldl_ut) & ukb_raw$ldl_ut > 4.9) |
              (!is.na(ukb_raw$non_hdl) & ukb_raw$non_hdl > 5.9) |
              (ukb_raw$on_statin == 1)
ukb <- ukb_raw[ukb_raw$lc == TRUE, ]

# eDLCN for UKB
ukb$dlcn <- ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut >= 8.5, 8,
            ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut >= 6.5, 5,
            ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut >= 5.0, 3,
            ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut >= 4.0, 1, 0))))
ukb$tendon_xanth <- 0; ukb$corneal_arcus <- 0
ukb$ldl_m <- ukb$ldl_m; ukb$cohort <- "UKB"

ukb <- compute_diag_features(as.data.frame(ukb))
ukb_df <- ukb[, cols]
cat(sprintf("  UKB LC: %d (FH+=%d, FH-=%d)\n", nrow(ukb_df), sum(ukb_df$fh), sum(ukb_df$fh == 0)))

# ── 5. COMBINE ──────────────────────────────────────────────────────────────

all_data <- rbind(sw, wales, ukb_df)
cat(sprintf("\n  TOTAL: %d (SW=%d, Wales=%d, UKB=%d)\n\n",
            nrow(all_data), nrow(sw), nrow(wales), nrow(ukb_df)))

# ══════════════════════════════════════════════════════════════════════════════
#  FEATURE TIERS — DIAGNOSTIC-SPECIFIC
#
#  The key diagnostic features for FH are:
#   - LDL (untreated) and trig_filter: core lipid signature
#   - HDL, TG: discriminators (FH has high LDL but NORMAL TG)
#   - Clinical signs: tendon xanthomata, corneal arcus (absent in UKB)
#   - Age, sex: confounders that need adjustment
#   - LDL/TG contrast: FH = isolated LDL elevation, not mixed dyslipidaemia
#   - ApoB/LDL ratio: reflects particle composition
# ══════════════════════════════════════════════════════════════════════════════

# ── STRATEGY A: FEATURE ENGINEERING ──────────────────────────────────────────

cat("================================================================\n")
cat("  STRATEGY A: FEATURE ENGINEERING (LOCO-CV)\n")
cat("================================================================\n")

feature_tiers <- list(
  # Tier 1: Original TUDOR (6 features)
  "TUDOR-Base" = c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex"),

  # Tier 2: + clinical signs + statin status
  "TUDOR-Clinical" = c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex",
                        "tendon_xanth", "corneal_arcus", "on_statin"),

  # Tier 3: + lipid ratios (LDL/HDL captures atherogenic phenotype)
  "TUDOR-Ratios" = c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex",
                      "on_statin", "ldl_hdl_ratio", "tg_hdl_ratio"),

  # Tier 4: + FH-specific contrast features
  "TUDOR-FHSig" = c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex",
                     "on_statin", "ldl_tg_contrast", "log_trig_filter",
                     "inv_hdl"),

  # Tier 5: + nonlinear age + composite physical sign
  "TUDOR-Enhanced" = c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex",
                        "tendon_xanth", "corneal_arcus", "on_statin",
                        "ldl_hdl_ratio", "ldl_tg_contrast", "age_sq",
                        "has_phys_sign"),

  # Tier 6: + cumulative LDL + young+high LDL interaction
  "TUDOR-Extended" = c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex",
                        "tendon_xanth", "corneal_arcus", "on_statin",
                        "ldl_hdl_ratio", "ldl_tg_contrast", "cumulative_ldl",
                        "young_high_ldl", "inv_hdl"),

  # Tier 7: Full kitchen sink — let elastic net select
  "TUDOR-Full" = c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex",
                    "tendon_xanth", "corneal_arcus", "on_statin",
                    "ldl_hdl_ratio", "tg_hdl_ratio", "ldl_tg_contrast",
                    "cumulative_ldl", "log_trig_filter", "inv_hdl",
                    "age_sq", "has_phys_sign", "young_high_ldl",
                    "nhdl_ldl_gap")
)

cohorts <- c("SouthWales", "Wales", "UKB")
strat_a <- list()

for (ho in cohorts) {
  cat(sprintf("\n  HELD-OUT: %s\n", ho))

  train <- all_data[all_data$cohort != ho, ]
  test  <- all_data[all_data$cohort == ho, ]

  for (tier_name in names(feature_tiers)) {
    feats <- feature_tiers[[tier_name]]

    train_cc <- train[complete.cases(train[, c("fh", feats)]), ]
    test_cc  <- test[complete.cases(test[, c("fh", feats)]), ]

    if (nrow(test_cc) < 50 || sum(test_cc$fh == 1) < 10) {
      cat(sprintf("    %-18s SKIPPED\n", tier_name)); next
    }

    X_tr <- as.matrix(train_cc[, feats]); y_tr <- train_cc$fh
    X_te <- as.matrix(test_cc[, feats]);  y_te <- test_cc$fh

    cv <- cv.glmnet(X_tr, y_tr, family = "binomial", alpha = 0.5,
                    nfolds = 10, type.measure = "auc", standardize = TRUE)
    pred <- as.numeric(predict(cv, newx = X_te, s = "lambda.min", type = "response"))

    r  <- roc(y_te, pred, quiet = TRUE)
    ci <- ci.auc(r, quiet = TRUE)

    # Also compute eDLCN for comparison
    edlcn_ok <- !is.na(test_cc$dlcn)
    auc_dlcn <- NA
    if (sum(edlcn_ok) > 50 && length(unique(test_cc$fh[edlcn_ok])) == 2) {
      r_d <- roc(test_cc$fh[edlcn_ok], test_cc$dlcn[edlcn_ok], quiet = TRUE)
      auc_dlcn <- round(as.numeric(auc(r_d)), 4)
    }

    flag <- ifelse(as.numeric(auc(r)) >= 0.82, " <<<", "")
    cat(sprintf("    %-18s AUC=%.4f [%.4f-%.4f] n=%d FH+=%d eDLCN=%.4f%s\n",
                tier_name, as.numeric(auc(r)), ci[1], ci[3],
                nrow(test_cc), sum(y_te), ifelse(is.na(auc_dlcn), 0, auc_dlcn), flag))

    strat_a[[length(strat_a) + 1]] <- data.frame(
      Strategy = "A_Features", HeldOut = ho, Model = tier_name,
      AUC = round(as.numeric(auc(r)), 4),
      CI_lo = round(ci[1], 4), CI_hi = round(ci[3], 4),
      N = nrow(test_cc), N_FH = sum(y_te), eDLCN_AUC = auc_dlcn,
      stringsAsFactors = FALSE
    )
  }
}

# ── STRATEGY B: ENSEMBLE STACKING ───────────────────────────────────────────

cat("\n\n================================================================\n")
cat("  STRATEGY B: ENSEMBLE STACKING (LOCO-CV)\n")
cat("================================================================\n")

# Use two feature sets for stacking
stack_tiers <- list(
  "Stack-Enhanced" = c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex",
                        "tendon_xanth", "corneal_arcus", "on_statin",
                        "ldl_hdl_ratio", "ldl_tg_contrast", "age_sq",
                        "has_phys_sign"),
  "Stack-Full" = c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex",
                    "tendon_xanth", "corneal_arcus", "on_statin",
                    "ldl_hdl_ratio", "tg_hdl_ratio", "ldl_tg_contrast",
                    "cumulative_ldl", "log_trig_filter", "inv_hdl",
                    "age_sq", "has_phys_sign", "young_high_ldl",
                    "nhdl_ldl_gap")
)

strat_b <- list()

for (ho in cohorts) {
  cat(sprintf("\n  HELD-OUT: %s\n", ho))

  train <- all_data[all_data$cohort != ho, ]
  test  <- all_data[all_data$cohort == ho, ]

  for (stack_name in names(stack_tiers)) {
    feats <- stack_tiers[[stack_name]]

    train_cc <- train[complete.cases(train[, c("fh", feats)]), ]
    test_cc  <- test[complete.cases(test[, c("fh", feats)]), ]

    if (nrow(test_cc) < 50 || sum(test_cc$fh == 1) < 10) {
      cat(sprintf("    %-18s SKIPPED\n", stack_name)); next
    }

    X_tr <- as.matrix(train_cc[, feats]); y_tr <- train_cc$fh
    X_te <- as.matrix(test_cc[, feats]);  y_te <- test_cc$fh
    n_pos <- sum(y_tr == 1); n_neg <- sum(y_tr == 0)

    # --- BASE 1: Elastic Net ---
    cv_enet <- cv.glmnet(X_tr, y_tr, family = "binomial", alpha = 0.5,
                          nfolds = 10, type.measure = "auc", standardize = TRUE)

    # --- BASE 2: Random Forest (class-balanced) ---
    df_rf <- data.frame(X_tr, fh = factor(y_tr, levels = c(0, 1)))
    rf_fit <- tryCatch({
      randomForest(fh ~ ., data = df_rf, ntree = 300,
                   mtry = max(2, floor(sqrt(ncol(X_tr)))),
                   sampsize = c("0" = min(n_neg, 3 * n_pos), "1" = n_pos),
                   strata = df_rf$fh)
    }, error = function(e) NULL)

    # --- BASE 3: XGBoost ---
    dtrain <- xgb.DMatrix(data = X_tr, label = y_tr)
    params <- list(objective = "binary:logistic", eval_metric = "auc",
                   max_depth = 3, eta = 0.05, min_child_weight = 10,
                   subsample = 0.8, colsample_bytree = 0.8,
                   gamma = 1, lambda = 5, alpha = 1,
                   scale_pos_weight = n_neg / n_pos)
    cv_xgb <- xgb.cv(params = params, data = dtrain, nrounds = 500,
                       nfold = 5, early_stopping_rounds = 50,
                       verbose = 0, stratified = TRUE)
    best_n <- cv_xgb$best_iteration
    if (is.null(best_n) || is.na(best_n)) best_n <- 100
    xgb_fit <- xgb.train(params = params, data = dtrain, nrounds = best_n, verbose = 0)

    # --- OUT-OF-FOLD PREDICTIONS (5-fold on training) ---
    n_tr <- nrow(train_cc)
    folds_cv <- sample(rep(1:5, length.out = n_tr))
    oof_enet <- oof_rf <- oof_xgb <- numeric(n_tr)

    for (k in 1:5) {
      ti <- folds_cv != k; vi <- folds_cv == k
      Xk_tr <- X_tr[ti, ]; yk_tr <- y_tr[ti]
      Xk_va <- X_tr[vi, ]
      npk <- sum(yk_tr == 1); nnk <- sum(yk_tr == 0)

      # ENET
      cvk <- cv.glmnet(Xk_tr, yk_tr, family = "binomial", alpha = 0.5,
                        nfolds = 5, type.measure = "auc", standardize = TRUE)
      oof_enet[vi] <- as.numeric(predict(cvk, newx = Xk_va, s = "lambda.min", type = "response"))

      # RF
      dfk <- data.frame(Xk_tr, fh = factor(yk_tr, levels = c(0, 1)))
      rfk <- tryCatch({
        randomForest(fh ~ ., data = dfk, ntree = 200,
                     mtry = max(2, floor(sqrt(ncol(Xk_tr)))),
                     sampsize = c("0" = min(nnk, 3 * npk), "1" = npk),
                     strata = dfk$fh)
      }, error = function(e) NULL)
      if (!is.null(rfk)) {
        oof_rf[vi] <- predict(rfk, newdata = data.frame(Xk_va), type = "prob")[, "1"]
      } else { oof_rf[vi] <- 0.5 }

      # XGB
      dtk <- xgb.DMatrix(data = Xk_tr, label = yk_tr)
      pk <- list(objective = "binary:logistic", eval_metric = "auc",
                 max_depth = 3, eta = 0.05, min_child_weight = 10,
                 subsample = 0.8, colsample_bytree = 0.8,
                 gamma = 1, lambda = 5, alpha = 1,
                 scale_pos_weight = nnk / npk)
      xk <- xgb.train(params = pk, data = dtk, nrounds = best_n, verbose = 0)
      oof_xgb[vi] <- predict(xk, xgb.DMatrix(data = Xk_va))
    }

    # --- META-LEARNER (logistic) ---
    meta_df <- data.frame(enet = oof_enet, rf = oof_rf, xgb = oof_xgb)
    meta_fit <- glm(y_tr ~ enet + rf + xgb, data = meta_df, family = binomial)

    # --- TEST PREDICTIONS ---
    p_enet <- as.numeric(predict(cv_enet, newx = X_te, s = "lambda.min", type = "response"))
    p_rf   <- if (!is.null(rf_fit)) {
      predict(rf_fit, newdata = data.frame(X_te), type = "prob")[, "1"]
    } else { rep(0.5, nrow(X_te)) }
    p_xgb  <- predict(xgb_fit, xgb.DMatrix(data = X_te))

    meta_te <- data.frame(enet = p_enet, rf = p_rf, xgb = p_xgb)
    pred_stack <- predict(meta_fit, newdata = meta_te, type = "response")

    r  <- roc(y_te, pred_stack, quiet = TRUE)
    ci <- ci.auc(r, quiet = TRUE)

    # Also report individual base learner AUCs
    r_e <- roc(y_te, p_enet, quiet = TRUE)
    r_r <- roc(y_te, p_rf, quiet = TRUE)
    r_x <- roc(y_te, p_xgb, quiet = TRUE)

    flag <- ifelse(as.numeric(auc(r)) >= 0.82, " <<<", "")
    cat(sprintf("    %-18s Stack=%.4f [%.4f-%.4f] | ENET=%.3f RF=%.3f XGB=%.3f | n=%d%s\n",
                stack_name, as.numeric(auc(r)), ci[1], ci[3],
                as.numeric(auc(r_e)), as.numeric(auc(r_r)), as.numeric(auc(r_x)),
                nrow(test_cc), flag))

    strat_b[[length(strat_b) + 1]] <- data.frame(
      Strategy = "B_Ensemble", HeldOut = ho, Model = stack_name,
      AUC = round(as.numeric(auc(r)), 4),
      CI_lo = round(ci[1], 4), CI_hi = round(ci[3], 4),
      AUC_ENET = round(as.numeric(auc(r_e)), 4),
      AUC_RF = round(as.numeric(auc(r_r)), 4),
      AUC_XGB = round(as.numeric(auc(r_x)), 4),
      N = nrow(test_cc), N_FH = sum(y_te),
      stringsAsFactors = FALSE
    )
  }
}

# ── STRATEGY C: POOLED BACKBONE + RECALIBRATION ─────────────────────────────

cat("\n\n================================================================\n")
cat("  STRATEGY C: POOLED + RECALIBRATION (LOCO-CV)\n")
cat("================================================================\n")

recal_feats <- c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex",
                  "tendon_xanth", "corneal_arcus", "on_statin",
                  "ldl_hdl_ratio", "ldl_tg_contrast", "inv_hdl",
                  "age_sq", "has_phys_sign", "cumulative_ldl")

strat_c <- list()

for (ho in cohorts) {
  cat(sprintf("\n  HELD-OUT: %s\n", ho))

  train <- all_data[all_data$cohort != ho, ]
  test  <- all_data[all_data$cohort == ho, ]

  train_cc <- train[complete.cases(train[, c("fh", recal_feats)]), ]
  test_cc  <- test[complete.cases(test[, c("fh", recal_feats)]), ]

  if (nrow(test_cc) < 50 || sum(test_cc$fh == 1) < 10) {
    cat(sprintf("    SKIPPED\n")); next
  }

  X_tr <- as.matrix(train_cc[, recal_feats]); y_tr <- train_cc$fh
  X_te <- as.matrix(test_cc[, recal_feats]);  y_te <- test_cc$fh

  # Pooled backbone
  cv_pool <- cv.glmnet(X_tr, y_tr, family = "binomial", alpha = 0.5,
                        nfolds = 10, type.measure = "auc", standardize = TRUE)
  pred_raw <- as.numeric(predict(cv_pool, newx = X_te, s = "lambda.min", type = "response"))
  lp_raw   <- as.numeric(predict(cv_pool, newx = X_te, s = "lambda.min", type = "link"))

  r_raw  <- roc(y_te, pred_raw, quiet = TRUE)
  ci_raw <- ci.auc(r_raw, quiet = TRUE)

  # Recalibration: fit intercept+slope on cross-predicted training LP
  train_cohorts <- unique(as.character(train_cc$cohort))
  recal_lp <- numeric(nrow(train_cc))
  for (tc in train_cohorts) {
    idx_tc <- which(train_cc$cohort == tc)
    idx_ot <- which(train_cc$cohort != tc)
    if (length(idx_ot) < 30 || length(idx_tc) < 20) {
      recal_lp[idx_tc] <- predict(cv_pool, newx = X_tr[idx_tc, , drop = FALSE],
                                   s = "lambda.min", type = "link")
      next
    }
    cv_ot <- cv.glmnet(X_tr[idx_ot, ], y_tr[idx_ot], family = "binomial",
                        alpha = 0.5, nfolds = min(10, length(idx_ot)),
                        type.measure = "auc", standardize = TRUE)
    recal_lp[idx_tc] <- as.numeric(predict(cv_ot, newx = X_tr[idx_tc, , drop = FALSE],
                                            s = "lambda.min", type = "link"))
  }

  recal_fit <- glm(y_tr ~ recal_lp, family = binomial)
  pred_recal <- predict(recal_fit, newdata = data.frame(recal_lp = lp_raw), type = "response")

  r_recal  <- roc(y_te, pred_recal, quiet = TRUE)
  ci_recal <- ci.auc(r_recal, quiet = TRUE)

  flag_r <- ifelse(as.numeric(auc(r_raw)) >= 0.82, " <<<", "")
  flag_c <- ifelse(as.numeric(auc(r_recal)) >= 0.82, " <<<", "")

  cat(sprintf("    Pooled-raw       AUC=%.4f [%.4f-%.4f] n=%d%s\n",
              as.numeric(auc(r_raw)), ci_raw[1], ci_raw[3], nrow(test_cc), flag_r))
  cat(sprintf("    Pooled-recal     AUC=%.4f [%.4f-%.4f] n=%d%s\n",
              as.numeric(auc(r_recal)), ci_recal[1], ci_recal[3], nrow(test_cc), flag_c))

  strat_c[[length(strat_c) + 1]] <- data.frame(
    Strategy = "C_Recal", HeldOut = ho, Model = "Pooled-raw",
    AUC = round(as.numeric(auc(r_raw)), 4),
    CI_lo = round(ci_raw[1], 4), CI_hi = round(ci_raw[3], 4),
    N = nrow(test_cc), N_FH = sum(y_te), stringsAsFactors = FALSE
  )
  strat_c[[length(strat_c) + 1]] <- data.frame(
    Strategy = "C_Recal", HeldOut = ho, Model = "Pooled-recal",
    AUC = round(as.numeric(auc(r_recal)), 4),
    CI_lo = round(ci_recal[1], 4), CI_hi = round(ci_recal[3], 4),
    N = nrow(test_cc), N_FH = sum(y_te), stringsAsFactors = FALSE
  )
}

# ══════════════════════════════════════════════════════════════════════════════
#  GRAND SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

cat("\n\n================================================================\n")
cat("  GRAND SUMMARY: FH DIAGNOSTIC MODEL OPTIMISATION\n")
cat("================================================================\n\n")

# Combine all — handle differing columns
strat_a_df <- do.call(rbind, strat_a)
strat_b_df <- do.call(rbind, strat_b)
strat_c_df <- do.call(rbind, strat_c)

# Standardise columns
core_cols <- c("Strategy", "HeldOut", "Model", "AUC", "CI_lo", "CI_hi", "N", "N_FH")
strat_a_slim <- strat_a_df[, core_cols]
strat_b_slim <- strat_b_df[, core_cols]
strat_c_slim <- strat_c_df[, core_cols]

all_results <- rbind(strat_a_slim, strat_b_slim, strat_c_slim)
all_results$Meets_082 <- ifelse(all_results$AUC >= 0.82, "YES", "no")

# Print by strategy
for (s in unique(all_results$Strategy)) {
  cat(sprintf("\n--- %s ---\n", s))
  sub <- all_results[all_results$Strategy == s, ]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    cat(sprintf("  %-12s %-18s AUC=%.4f [%.4f-%.4f] n=%-5d FH+=%d %s\n",
                r$HeldOut, r$Model, r$AUC, r$CI_lo, r$CI_hi,
                r$N, r$N_FH, r$Meets_082))
  }
}

# Best per cohort
cat("\n\n=== BEST MODEL PER COHORT ===\n")
for (ho in cohorts) {
  sub <- all_results[all_results$HeldOut == ho & !is.na(all_results$AUC), ]
  if (nrow(sub) == 0) next
  best <- sub[which.max(sub$AUC), ]
  flag <- ifelse(best$AUC >= 0.82, "TARGET MET", "BELOW 0.82")
  cat(sprintf("  %-12s  %s/%s  AUC=%.4f [%.4f-%.4f]  %s\n",
              ho, best$Strategy, best$Model, best$AUC, best$CI_lo, best$CI_hi, flag))
}

# Universal winners
cat("\n=== MODELS >= 0.82 IN ALL 3 COHORTS ===\n")
model_keys <- paste(all_results$Strategy, all_results$Model, sep = "/")
any_winner <- FALSE
for (mk in unique(model_keys)) {
  sub <- all_results[model_keys == mk, ]
  if (nrow(sub) == 3 && all(sub$AUC >= 0.82, na.rm = TRUE)) {
    cat(sprintf("  %s: SW=%.3f UKB=%.3f Wales=%.3f\n", mk,
                sub$AUC[sub$HeldOut == "SouthWales"],
                sub$AUC[sub$HeldOut == "UKB"],
                sub$AUC[sub$HeldOut == "Wales"]))
    any_winner <- TRUE
  }
}
if (!any_winner) {
  cat("  NONE yet. Closest models:\n")
  for (mk in unique(model_keys)) {
    sub <- all_results[model_keys == mk, ]
    if (nrow(sub) == 3) {
      min_auc <- min(sub$AUC, na.rm = TRUE)
      if (min_auc >= 0.78) {
        cat(sprintf("  %s: min=%.3f | SW=%.3f UKB=%.3f Wales=%.3f\n", mk,
                    min_auc,
                    sub$AUC[sub$HeldOut == "SouthWales"],
                    sub$AUC[sub$HeldOut == "UKB"],
                    sub$AUC[sub$HeldOut == "Wales"]))
      }
    }
  }
}

# Save all
write.csv(all_results, file.path(OUT, "diagnostic_AUC_optimisation.csv"), row.names = FALSE)
# Save detailed ensemble results separately
if (nrow(strat_b_df) > 0)
  write.csv(strat_b_df, file.path(OUT, "diagnostic_ensemble_detail.csv"), row.names = FALSE)

cat(sprintf("\n  Saved: %s/diagnostic_AUC_optimisation.csv\n", OUT))
cat(sprintf("  Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("================================================================\n")

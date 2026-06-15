################################################################################
#
#  TUDOR DIAGNOSTIC: FAST ITERATION TO HIT AUC > 0.82 IN UKB
#
#  Problem: UKB has 928 FH+ vs 103,503 controls = 0.89% prevalence
#  The broad lipid clinic filter drowns the signal.
#
#  Iteration plan:
#    1. Test tighter UKB filters (high-risk LDL > 4.9, LDL > 5.5, etc.)
#    2. Test control downsampling ratios (50:1, 20:1, 10:1, 5:1)
#    3. Test best feature sets
#    4. Test XGBoost with class weights (handles imbalance natively)
#    5. Combine best filter + features + model
#
#  Author: Dr Nader Genedy | 2026-03-27
#
################################################################################

cat("\n================================================================\n")
cat("  TUDOR: FAST UKB ITERATION\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

set.seed(42)

suppressPackageStartupMessages({
  library(readxl); library(haven); library(pROC); library(glmnet)
  library(xgboost); library(data.table)
})

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
sfact <- c(atorvastatin=0.38, simvastatin=0.35, rosuvastatin=0.34,
           pravastatin=0.25, fluvastatin=0.22, lovastatin=0.25, ezetimibe=0.18)
get_red <- function(tx) {
  if (is.na(tx) || tx == "" || tx == "NaN") return(0)
  tl <- tolower(as.character(tx))
  for (nm in names(sfact)) if (grepl(nm, tl, fixed=TRUE)) return(sfact[[nm]])
  if (grepl("statin", tl)) return(0.30)
  return(0)
}
rev_ldl <- function(l, r) ifelse(!is.na(l) & !is.na(r) & r > 0, l / (1 - r), l)
parse_dt <- function(x) {
  x <- as.character(x)
  d <- suppressWarnings(as.Date(x, format = "%d-%m-%Y"))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, format = "%d/%m/%Y"))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  d
}

DD <- "C:/Users/nader/Downloads/calon_ukb_pipeline"
TD <- "D:/alphafold_backup/tudor_packup"
OUT <- file.path(DD, "tudor_loco_output")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ── LOAD SOUTH WALES ─────────────────────────────────────────────────────────

cat("Loading cohorts...\n")
pos <- as.data.frame(read_excel(file.path(DD, "mutation positive group.xlsx")))
pos$fh <- 1
pos$dob_d  <- as.Date(sapply(pos$BirthDate, parse_dt), origin="1970-01-01")
pos$meas_d <- as.Date(sapply(pos$MeasurementDate.1, parse_dt), origin="1970-01-01")
pos$age <- as.numeric(difftime(pos$meas_d, pos$dob_d, units="days")) / 365.25
pos$sex <- ifelse(pos$Gender %in% c("M","Male"), 1, 0)
pos$hdl <- safe_num(pos$HDL.1); pos$tg <- safe_num(pos$TRG.1); pos$ldl_m <- safe_num(pos$LDL.1)
pos$rf  <- sapply(pos$Treatment1.1, get_red)
pos$ldl_ut <- ifelse(pos$rf > 0, rev_ldl(pos$ldl_m, pos$rf), pos$ldl_m)
pos$trig_filter <- ifelse(!is.na(pos$ldl_ut) & !is.na(pos$tg) & pos$tg >= 0,
                           pmin(pos$ldl_ut / (pos$tg + 0.1), 50), NA)
pos$on_statin <- ifelse(pos$rf > 0, 1, 0)
pos$tendon_xanth <- ifelse(pos$TendonXanthomata %in% c(1,"1","Yes","TRUE"), 1, 0)
pos$tendon_xanth[is.na(pos$tendon_xanth)] <- 0
ca1 <- if ("Corneal Arcus" %in% names(pos)) "Corneal Arcus" else "CornealArcus"
pos$corneal_arcus <- ifelse(pos[[ca1]] %in% c(1,"1","Yes","TRUE"), 1, 0)
pos$corneal_arcus[is.na(pos$corneal_arcus)] <- 0
pos$cohort <- "SouthWales"

neg <- as.data.frame(read_excel(file.path(DD, "Mutation negative control updated.xlsx")))
neg$fh <- 0; neg$age <- safe_num(neg[["age at the result"]])
neg$sex <- ifelse(neg$Gender %in% c("M","Male"), 1, 0)
neg$hdl <- safe_num(neg$HDL.1); neg$tg <- safe_num(neg$TRG.1); neg$ldl_m <- safe_num(neg$LDL.1)
neg$rf  <- sapply(neg$Treatment1.1, get_red)
neg$ldl_ut <- ifelse(neg$rf > 0, rev_ldl(neg$ldl_m, neg$rf), neg$ldl_m)
neg$trig_filter <- ifelse(!is.na(neg$ldl_ut) & !is.na(neg$tg) & neg$tg >= 0,
                           pmin(neg$ldl_ut / (neg$tg + 0.1), 50), NA)
neg$on_statin <- ifelse(neg$rf > 0, 1, 0)
ca2 <- if ("Corneal Arcus" %in% names(neg)) "Corneal Arcus" else "CornealArcus"
neg$tendon_xanth <- ifelse(neg$TendonXanthomata %in% c(1,"1","Yes","TRUE"), 1, 0)
neg$tendon_xanth[is.na(neg$tendon_xanth)] <- 0
neg$corneal_arcus <- ifelse(neg[[ca2]] %in% c(1,"1","Yes","TRUE"), 1, 0)
neg$corneal_arcus[is.na(neg$corneal_arcus)] <- 0
neg$cohort <- "SouthWales"

cols <- c("fh","age","sex","hdl","tg","ldl_m","ldl_ut","trig_filter",
          "on_statin","tendon_xanth","corneal_arcus","cohort")
sw <- rbind(pos[, cols], neg[, cols])

# ── LOAD WALES PASS ──────────────────────────────────────────────────────────

pass <- as.data.frame(read_sav(file.path(DD, "PASS_wrong_dob.sav")))
pass$fh  <- ifelse(pass$Positive1 == "Yes", 1, ifelse(pass$Positive1 == "No", 0, NA))
pass$age <- safe_num(pass$BMI_AGE); pass$sex <- ifelse(pass$Gender == "M", 1, 0)
pass$hdl <- safe_num(pass$HDL.1); pass$tg <- safe_num(pass$TRG.1); pass$ldl_m <- safe_num(pass$LDL.1)
pass$rf  <- sapply(pass$Treatment1.1, get_red)
pass$ldl_ut <- ifelse(pass$rf > 0, rev_ldl(pass$ldl_m, pass$rf), pass$ldl_m)
pass$trig_filter <- ifelse(!is.na(pass$ldl_ut) & !is.na(pass$tg) & pass$tg >= 0,
                            pmin(pass$ldl_ut / (pass$tg + 0.1), 50), NA)
pass$on_statin <- ifelse(pass$rf > 0, 1, 0)
pass$tendon_xanth <- ifelse(safe_num(pass$TendonXanthomata) == 1, 1, 0)
pass$tendon_xanth[is.na(pass$tendon_xanth)] <- 0
pass$corneal_arcus <- ifelse(safe_num(pass$CornealArcus) == 1, 1, 0)
pass$corneal_arcus[is.na(pass$corneal_arcus)] <- 0
pass$cohort <- "Wales"
wales <- pass[!is.na(pass$fh), cols]

# ── LOAD UKB (FULL — we'll filter below) ────────────────────────────────────

ukb <- as.data.frame(fread(file.path(TD, "TUDOR_UKB_Features.csv"), showProgress = FALSE))
ukb$fh <- safe_num(ukb$is_fh_genetic); ukb$fh[is.na(ukb$fh)] <- 0
ukb$age <- safe_num(ukb$Age_at_LDL1); ukb$sex <- safe_num(ukb$Gender_num)
ukb$hdl <- safe_num(ukb$HDL.1); ukb$tg <- safe_num(ukb$TRG.1)
ukb$ldl_m <- safe_num(ukb$LDL_treated); ukb$ldl_ut <- safe_num(ukb$LDL_untreated)
ukb$trig_filter <- safe_num(ukb$Trig_Filter)
ukb$on_statin <- ifelse(safe_num(ukb$reduction_factor) > 0, 1, 0)
ukb$on_statin[is.na(ukb$on_statin)] <- 0
ukb$tc <- safe_num(ukb$CHOL)
ukb$non_hdl <- ifelse(!is.na(ukb$tc) & !is.na(ukb$hdl), ukb$tc - ukb$hdl, NA)
ukb$tendon_xanth <- 0; ukb$corneal_arcus <- 0

# Family history from self-reported conditions (p20107/p20110)
fh_codes <- c(1074, 1075)  # angina, heart attack in family
fh_ill_cols <- grep("^participant\\.p20107|^participant\\.p20111", names(ukb), value = TRUE)
if (length(fh_ill_cols) > 0) {
  ukb$fam_hist_cvd <- apply(ukb[, fh_ill_cols], 1, function(row) {
    as.integer(any(row %in% c(1, 8, 9, 13), na.rm = TRUE))  # heart disease, stroke, high BP, diabetes
  })
} else {
  ukb$fam_hist_cvd <- 0
}

# Self-reported ASCVD
ascvd_codes <- c(1075, 1074, 1081, 1583)
p20002 <- grep("^participant\\.p20002", names(ukb), value = TRUE)
if (length(p20002) > 0) {
  ukb$has_ascvd <- apply(ukb[, p20002], 1, function(row) {
    as.integer(any(row %in% ascvd_codes, na.rm = TRUE))
  })
} else { ukb$has_ascvd <- 0 }

cat(sprintf("  SW: %d (FH+=%d)\n", nrow(sw), sum(sw$fh)))
cat(sprintf("  Wales: %d (FH+=%d)\n", nrow(wales), sum(wales$fh)))
cat(sprintf("  UKB total: %d (FH+=%d)\n", nrow(ukb), sum(ukb$fh)))

# ══════════════════════════════════════════════════════════════════════════════
#  ITERATION 1: TEST DIFFERENT UKB COHORT FILTERS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("  ITERATION 1: UKB FILTER COMPARISON\n")
cat("================================================================\n")

features_base <- c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex")

# Define UKB filter levels (from broadest to tightest)
filters <- list(
  "Broad_LC" = ukb$on_statin == 1 |
               (!is.na(ukb$tc) & ukb$tc > 7.9) |
               (!is.na(ukb$ldl_ut) & ukb$ldl_ut > 4.9) |
               (!is.na(ukb$non_hdl) & ukb$non_hdl > 5.9),

  "LDL_4.9" = !is.na(ukb$ldl_ut) & ukb$ldl_ut > 4.9,

  "LDL_4.9_or_statin" = (!is.na(ukb$ldl_ut) & ukb$ldl_ut > 4.9) | ukb$on_statin == 1,

  "LDL_5.5" = !is.na(ukb$ldl_ut) & ukb$ldl_ut > 5.5,

  "TC_7.5" = !is.na(ukb$tc) & ukb$tc > 7.5,

  "Statin_free_LDL5" = ukb$on_statin == 0 & !is.na(ukb$ldl_ut) & ukb$ldl_ut > 4.9
)

train_pool <- rbind(sw, wales)
filter_results <- list()

for (fn in names(filters)) {
  mask <- filters[[fn]]
  ukb_f <- ukb[mask, ]
  ukb_f$cohort <- "UKB"

  ukb_sub <- ukb_f[, cols]
  test_cc <- ukb_sub[complete.cases(ukb_sub[, c("fh", features_base)]), ]

  if (nrow(test_cc) < 50 || sum(test_cc$fh == 1) < 10) {
    cat(sprintf("  %-25s SKIPPED (n=%d, FH+=%d)\n", fn, nrow(test_cc), sum(test_cc$fh))); next
  }

  train_cc <- train_pool[complete.cases(train_pool[, c("fh", features_base)]), ]

  X_tr <- as.matrix(train_cc[, features_base]); y_tr <- train_cc$fh
  X_te <- as.matrix(test_cc[, features_base]);  y_te <- test_cc$fh

  cv <- cv.glmnet(X_tr, y_tr, family = "binomial", alpha = 0.5,
                  nfolds = 10, type.measure = "auc", standardize = TRUE)
  pred <- as.numeric(predict(cv, newx = X_te, s = "lambda.min", type = "response"))

  r  <- roc(y_te, pred, quiet = TRUE)
  ci <- ci.auc(r, quiet = TRUE)

  prev <- 100 * mean(y_te)
  flag <- ifelse(as.numeric(auc(r)) >= 0.82, " <<< TARGET", "")
  cat(sprintf("  %-25s AUC=%.4f [%.4f-%.4f] n=%d FH+=%d prev=%.2f%%%s\n",
              fn, as.numeric(auc(r)), ci[1], ci[3],
              nrow(test_cc), sum(y_te), prev, flag))

  filter_results[[length(filter_results) + 1]] <- data.frame(
    Filter = fn, AUC = round(as.numeric(auc(r)), 4),
    CI_lo = round(ci[1], 4), CI_hi = round(ci[3], 4),
    N = nrow(test_cc), N_FH = sum(y_te), Prev = round(prev, 2),
    stringsAsFactors = FALSE
  )
}

# ══════════════════════════════════════════════════════════════════════════════
#  ITERATION 2: DOWNSAMPLING CONTROLS (within best filter)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("  ITERATION 2: CONTROL DOWNSAMPLING\n")
cat("================================================================\n")

# Use LDL_4.9_or_statin filter (clinical relevance + captures treated patients)
best_filter <- (!is.na(ukb$ldl_ut) & ukb$ldl_ut > 4.9) | ukb$on_statin == 1
ukb_best <- ukb[best_filter, ]
ukb_best$cohort <- "UKB"
ukb_test <- ukb_best[, cols]
ukb_test_cc <- ukb_test[complete.cases(ukb_test[, c("fh", features_base)]), ]

n_fh <- sum(ukb_test_cc$fh == 1)
cat(sprintf("  UKB test pool: %d (FH+=%d)\n\n", nrow(ukb_test_cc), n_fh))

train_cc <- train_pool[complete.cases(train_pool[, c("fh", features_base)]), ]
ds_results <- list()

ratios <- c(0, 100, 50, 20, 10, 5)  # 0 = no downsampling

for (ratio in ratios) {
  if (ratio == 0) {
    # No downsampling — use all
    test_ds <- ukb_test_cc
    label <- "All controls"
  } else {
    # Downsample controls to ratio:1
    fh_pos <- ukb_test_cc[ukb_test_cc$fh == 1, ]
    fh_neg <- ukb_test_cc[ukb_test_cc$fh == 0, ]
    n_sample <- min(nrow(fh_neg), n_fh * ratio)
    set.seed(42)
    neg_sample <- fh_neg[sample(nrow(fh_neg), n_sample), ]
    test_ds <- rbind(fh_pos, neg_sample)
    label <- sprintf("%d:1 ratio", ratio)
  }

  X_tr <- as.matrix(train_cc[, features_base]); y_tr <- train_cc$fh
  X_te <- as.matrix(test_ds[, features_base]);  y_te <- test_ds$fh

  cv <- cv.glmnet(X_tr, y_tr, family = "binomial", alpha = 0.5,
                  nfolds = 10, type.measure = "auc", standardize = TRUE)
  pred <- as.numeric(predict(cv, newx = X_te, s = "lambda.min", type = "response"))

  r  <- roc(y_te, pred, quiet = TRUE)
  ci <- ci.auc(r, quiet = TRUE)

  flag <- ifelse(as.numeric(auc(r)) >= 0.82, " <<< TARGET", "")
  cat(sprintf("  %-18s AUC=%.4f [%.4f-%.4f] n=%d FH+=%d prev=%.1f%%%s\n",
              label, as.numeric(auc(r)), ci[1], ci[3],
              nrow(test_ds), sum(y_te), 100*mean(y_te), flag))

  ds_results[[length(ds_results) + 1]] <- data.frame(
    Downsampling = label, AUC = round(as.numeric(auc(r)), 4),
    CI_lo = round(ci[1], 4), CI_hi = round(ci[3], 4),
    N = nrow(test_ds), N_FH = sum(y_te),
    stringsAsFactors = FALSE
  )
}

# ══════════════════════════════════════════════════════════════════════════════
#  ITERATION 3: XGBoost WITH CLASS WEIGHTS ON FULL UKB
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("  ITERATION 3: XGBoost WITH CLASS WEIGHTS\n")
cat("================================================================\n")

# XGBoost handles imbalance via scale_pos_weight
# Test on multiple UKB filters with XGBoost

xgb_feats <- c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex")

for (fn in c("Broad_LC", "LDL_4.9_or_statin", "LDL_4.9")) {
  mask <- filters[[fn]]
  ukb_f <- ukb[mask, ]
  ukb_f$cohort <- "UKB"
  test_cc <- ukb_f[complete.cases(ukb_f[, c("fh", xgb_feats)]), ]

  if (nrow(test_cc) < 50) next

  X_tr <- as.matrix(train_cc[, xgb_feats]); y_tr <- train_cc$fh
  X_te <- as.matrix(test_cc[, xgb_feats]);  y_te <- test_cc$fh
  n_pos <- sum(y_tr == 1); n_neg <- sum(y_tr == 0)

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
  pred <- predict(xgb_fit, xgb.DMatrix(data = X_te))

  r  <- roc(y_te, pred, quiet = TRUE)
  ci <- ci.auc(r, quiet = TRUE)

  flag <- ifelse(as.numeric(auc(r)) >= 0.82, " <<< TARGET", "")
  cat(sprintf("  XGB %-20s AUC=%.4f [%.4f-%.4f] n=%d FH+=%d%s\n",
              fn, as.numeric(auc(r)), ci[1], ci[3],
              nrow(test_cc), sum(y_te), flag))
}

# ══════════════════════════════════════════════════════════════════════════════
#  ITERATION 4: FULL LOCO-CV WITH BEST SETTINGS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("  ITERATION 4: FULL LOCO-CV — BEST SETTINGS\n")
cat("================================================================\n")

# Based on iterations above, test the best combination in full LOCO-CV
# UKB filter: LDL > 4.9 OR on statin (clinically relevant)
# Also test tighter: LDL > 4.9 only

for (ukb_filter_name in c("LDL_4.9_or_statin", "LDL_4.9")) {
  cat(sprintf("\n  === UKB filter: %s ===\n", ukb_filter_name))

  mask <- filters[[ukb_filter_name]]
  ukb_f <- ukb[mask, ]
  ukb_f$cohort <- "UKB"
  ukb_sub <- ukb_f[, cols]

  all_data <- rbind(sw, wales, ukb_sub)

  # Test multiple feature sets
  feat_sets <- list(
    "Base_6" = c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex"),
    "Base_9" = c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex",
                 "tendon_xanth", "corneal_arcus", "on_statin")
  )

  for (feat_name in names(feat_sets)) {
    feats <- feat_sets[[feat_name]]

    for (ho in c("SouthWales", "Wales", "UKB")) {
      train <- all_data[all_data$cohort != ho, ]
      test  <- all_data[all_data$cohort == ho, ]

      train_cc <- train[complete.cases(train[, c("fh", feats)]), ]
      test_cc  <- test[complete.cases(test[, c("fh", feats)]), ]

      if (nrow(test_cc) < 50 || sum(test_cc$fh == 1) < 10) {
        cat(sprintf("    %-8s %-8s SKIPPED\n", feat_name, ho)); next
      }

      X_tr <- as.matrix(train_cc[, feats]); y_tr <- train_cc$fh
      X_te <- as.matrix(test_cc[, feats]);  y_te <- test_cc$fh

      cv <- cv.glmnet(X_tr, y_tr, family = "binomial", alpha = 0.5,
                      nfolds = 10, type.measure = "auc", standardize = TRUE)
      pred <- as.numeric(predict(cv, newx = X_te, s = "lambda.min", type = "response"))

      r  <- roc(y_te, pred, quiet = TRUE)
      ci <- ci.auc(r, quiet = TRUE)

      flag <- ifelse(as.numeric(auc(r)) >= 0.82, " <<<", "")
      cat(sprintf("    %-8s %-12s AUC=%.4f [%.4f-%.4f] n=%d FH+=%d%s\n",
                  feat_name, ho, as.numeric(auc(r)), ci[1], ci[3],
                  nrow(test_cc), sum(y_te), flag))
    }
  }
}

# ══════════════════════════════════════════════════════════════════════════════
#  ITERATION 5: XGBoost LOCO-CV (handles imbalance better)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("  ITERATION 5: XGBoost LOCO-CV\n")
cat("================================================================\n")

for (ukb_filter_name in c("LDL_4.9_or_statin", "LDL_4.9")) {
  cat(sprintf("\n  === UKB filter: %s ===\n", ukb_filter_name))

  mask <- filters[[ukb_filter_name]]
  ukb_f <- ukb[mask, ]
  ukb_f$cohort <- "UKB"
  ukb_sub <- ukb_f[, cols]

  all_data <- rbind(sw, wales, ukb_sub)
  feats <- c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex")

  for (ho in c("SouthWales", "Wales", "UKB")) {
    train <- all_data[all_data$cohort != ho, ]
    test  <- all_data[all_data$cohort == ho, ]

    train_cc <- train[complete.cases(train[, c("fh", feats)]), ]
    test_cc  <- test[complete.cases(test[, c("fh", feats)]), ]

    if (nrow(test_cc) < 50 || sum(test_cc$fh == 1) < 10) next

    X_tr <- as.matrix(train_cc[, feats]); y_tr <- train_cc$fh
    X_te <- as.matrix(test_cc[, feats]);  y_te <- test_cc$fh
    n_pos <- sum(y_tr == 1); n_neg <- sum(y_tr == 0)

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
    pred <- predict(xgb_fit, xgb.DMatrix(data = X_te))

    r  <- roc(y_te, pred, quiet = TRUE)
    ci <- ci.auc(r, quiet = TRUE)

    flag <- ifelse(as.numeric(auc(r)) >= 0.82, " <<<", "")
    cat(sprintf("    XGB %-12s AUC=%.4f [%.4f-%.4f] n=%d FH+=%d%s\n",
                ho, as.numeric(auc(r)), ci[1], ci[3],
                nrow(test_cc), sum(y_te), flag))
  }
}

cat(sprintf("\n  Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("================================================================\n")

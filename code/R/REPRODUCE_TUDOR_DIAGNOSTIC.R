################################################################################
#
#  TUDOR v2: DIAGNOSTIC FH MODEL — FULL REPRODUCTION
#
#  TRAINING:  mutation positive group.xlsx (425 FH+, South Wales)
#           + Mutation negative control updated.xlsx (941 FH-, cascade negatives)
#
#  VALIDATION: PASS_wrong_dob.sav (All Wales, 7368 patients, 2405 FH+ / 4848 FH-)
#
#  Model: Elastic Net (alpha=0.5), 10-fold CV
#  Target AUC: ~0.85
#
#  Author: Dr Nader Genedy
#  Reproduction: 2026-03-26
#
################################################################################

cat("\n================================================================\n")
cat("  TUDOR v2: DIAGNOSTIC FH MODEL REPRODUCTION\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

set.seed(42)

# ── 0. PACKAGES ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(readxl); library(haven); library(pROC); library(glmnet)
})

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

data_dir <- "C:/Users/nader/Downloads/calon_ukb_pipeline"

# ── 1. STATIN CORRECTION FACTORS ────────────────────────────────────────────

statin_tbl <- c(
  simvastatin = 0.35, atorvastatin = 0.38, rosuvastatin = 0.34,
  pravastatin = 0.23, fluvastatin = 0.21, lovastatin = 0.25, ezetimibe = 0.18
)

get_reduction <- function(tx) {
  if (is.na(tx) || tx == "" || tx == "NaN") return(0)
  tx_low <- tolower(as.character(tx))
  for (nm in names(statin_tbl)) {
    if (grepl(nm, tx_low, fixed = TRUE)) return(statin_tbl[[nm]])
  }
  # If any statin text detected but no match, use generic 0.35
  if (grepl("statin", tx_low)) return(0.35)
  return(0)
}

reverse_ldl <- function(ldl, rf) {
  ifelse(!is.na(ldl) & !is.na(rf) & rf > 0, ldl / (1 - rf), ldl)
}

# ── 2. LOAD TRAINING DATA ───────────────────────────────────────────────────

cat("=== LOADING TRAINING DATA ===\n")

# --- Mutation POSITIVE (FH+) ---
pos <- as.data.frame(read_excel(file.path(data_dir, "mutation positive group.xlsx")))
cat(sprintf("  Mutation positive: %d rows x %d cols\n", nrow(pos), ncol(pos)))

pos$fh <- 1

# Age: compute from BirthDate and MeasurementDate.1
parse_date <- function(x) {
  x <- as.character(x)
  # Try multiple formats
  d <- suppressWarnings(as.Date(x, format = "%d-%m-%Y"))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, format = "%d/%m/%Y"))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  return(d)
}
pos$dob_parsed <- sapply(pos$BirthDate, parse_date)
pos$dob_parsed <- as.Date(pos$dob_parsed, origin = "1970-01-01")
pos$meas_parsed <- sapply(pos$MeasurementDate.1, parse_date)
pos$meas_parsed <- as.Date(pos$meas_parsed, origin = "1970-01-01")
pos$age <- as.numeric(difftime(pos$meas_parsed, pos$dob_parsed, units = "days")) / 365.25
cat(sprintf("  Pos age: %d non-NA, median=%.0f\n", sum(!is.na(pos$age)), median(pos$age, na.rm=TRUE)))
pos$sex <- ifelse(pos$Gender == "M" | pos$Gender == "Male", 1,
            ifelse(pos$Gender == "F" | pos$Gender == "Female", 0, NA))

# Lipids at Visit 1
pos$hdl   <- safe_num(pos$HDL.1)
pos$tg    <- safe_num(pos$TRG.1)
pos$ldl_m <- safe_num(pos$LDL.1)

# Statin correction
pos$red_factor <- sapply(pos$Treatment1.1, get_reduction)
pos$ldl_ut <- reverse_ldl(pos$ldl_m, pos$red_factor)

# Trig Filter
pos$trig_filter <- ifelse(!is.na(pos$ldl_ut) & !is.na(pos$tg) & pos$tg >= 0,
                           pmin(pos$ldl_ut / (pos$tg + 0.1), 50), NA)

# On statin flag
pos$on_statin <- ifelse(pos$red_factor > 0, 1, 0)

# Clinical signs
pos$tendon_xanth <- ifelse(pos$TendonXanthomata %in% c(1, "1", "TRUE", "True", "Yes", "YES"), 1, 0)
pos$tendon_xanth[is.na(pos$tendon_xanth)] <- 0

ca_col <- if ("Corneal Arcus" %in% names(pos)) "Corneal Arcus" else "CornealArcus"
pos$corneal_arcus <- ifelse(pos[[ca_col]] %in% c(1, "1", "TRUE", "True", "Yes", "YES"), 1, 0)
pos$corneal_arcus[is.na(pos$corneal_arcus)] <- 0

# DLCN score
pos$dlcn <- safe_num(pos$GenoTypingScore)

# Index effect — no explicit column in pos, but ALL are index patients (probands)
pos$index_effect <- 1

cat(sprintf("  Pos: LDL non-NA=%d, HDL non-NA=%d, TG non-NA=%d, Age non-NA=%d\n",
            sum(!is.na(pos$ldl_m)), sum(!is.na(pos$hdl)),
            sum(!is.na(pos$tg)), sum(!is.na(pos$age))))

# --- Mutation NEGATIVE (FH-) ---
neg <- as.data.frame(read_excel(file.path(data_dir, "Mutation negative control updated.xlsx")))
cat(sprintf("  Mutation negative: %d rows x %d cols\n", nrow(neg), ncol(neg)))

neg$fh <- 0
neg$age <- safe_num(neg$`age at the result`)
neg$sex <- ifelse(neg$Gender == "M" | neg$Gender == "Male", 1,
            ifelse(neg$Gender == "F" | neg$Gender == "Female", 0, NA))

# Lipids at Visit 1
neg$hdl   <- safe_num(neg$HDL.1)
neg$tg    <- safe_num(neg$TRG.1)
neg$ldl_m <- safe_num(neg$LDL.1)

# Statin correction
neg$red_factor <- sapply(neg$Treatment1.1, get_reduction)
neg$ldl_ut <- reverse_ldl(neg$ldl_m, neg$red_factor)

# Trig Filter
neg$trig_filter <- ifelse(!is.na(neg$ldl_ut) & !is.na(neg$tg) & neg$tg >= 0,
                           pmin(neg$ldl_ut / (neg$tg + 0.1), 50), NA)

# On statin flag
neg$on_statin <- ifelse(neg$red_factor > 0, 1, 0)

# Clinical signs
neg$tendon_xanth <- ifelse(neg$TendonXanthomata %in% c(1, "1", "TRUE", "True", "Yes", "YES"), 1, 0)
neg$tendon_xanth[is.na(neg$tendon_xanth)] <- 0

ca_col_neg <- if ("Corneal Arcus" %in% names(neg)) "Corneal Arcus" else "CornealArcus"
neg$corneal_arcus <- ifelse(neg[[ca_col_neg]] %in% c(1, "1", "TRUE", "True", "Yes", "YES"), 1, 0)
neg$corneal_arcus[is.na(neg$corneal_arcus)] <- 0

# DLCN score
neg$dlcn <- safe_num(neg$GenoTypingScore)

# Index effect — negatives are cascade relatives
neg$index_effect <- 0

cat(sprintf("  Neg: LDL non-NA=%d, HDL non-NA=%d, TG non-NA=%d, Age non-NA=%d\n",
            sum(!is.na(neg$ldl_m)), sum(!is.na(neg$hdl)),
            sum(!is.na(neg$tg)), sum(!is.na(neg$age))))

# ── 3. COMBINE TRAINING SET ─────────────────────────────────────────────────

# NOTE: index_effect EXCLUDED from training features because it perfectly
# separates training classes (all pos=index, all neg=cascade) causing data leakage.
# It can be added as a post-hoc adjustment in validation where both groups
# contain index + cascade patients.
# statin_residual is mathematically degenerate: ldl_m - ldl_ut*(1-rf) ≈ 0 always
# Replace with binary on_statin flag (clinically meaningful, no collinearity)
features <- c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex",
              "tendon_xanth", "corneal_arcus", "on_statin")

keep_cols <- c("fh", features, "dlcn")

train_df <- rbind(pos[, keep_cols], neg[, keep_cols])
train_complete <- train_df[complete.cases(train_df[, c("fh", features)]), ]

cat(sprintf("\n  Combined training: %d total, %d complete cases\n",
            nrow(train_df), nrow(train_complete)))
cat(sprintf("  FH+: %d (%.1f%%), FH-: %d (%.1f%%)\n",
            sum(train_complete$fh == 1), 100 * mean(train_complete$fh == 1),
            sum(train_complete$fh == 0), 100 * mean(train_complete$fh == 0)))

# ── 4. TRAIN TUDOR MODEL ────────────────────────────────────────────────────

cat("\n=== TRAINING TUDOR MODEL ===\n")

X_train <- as.matrix(train_complete[, features])
y_train <- train_complete$fh

cv_fit <- cv.glmnet(X_train, y_train, family = "binomial",
                     alpha = 0.5, nfolds = 10,
                     type.measure = "auc", standardize = TRUE)

cat(sprintf("  Lambda.min: %.6f  |  Lambda.1se: %.6f\n",
            cv_fit$lambda.min, cv_fit$lambda.1se))
cat(sprintf("  10-fold CV AUC: %.4f\n", max(cv_fit$cvm)))

# Coefficients
coefs <- coef(cv_fit, s = "lambda.min")
coef_df <- data.frame(
  Feature = c("(Intercept)", features),
  Coefficient = as.numeric(coefs)
)
coef_df <- coef_df[order(-abs(coef_df$Coefficient)), ]

cat("\n  TUDOR v2 Coefficients:\n")
for (i in seq_len(nrow(coef_df))) {
  cat(sprintf("    %-20s %+.4f\n", coef_df$Feature[i], coef_df$Coefficient[i]))
}

# Training apparent AUC
pred_train <- as.numeric(predict(cv_fit, newx = X_train, s = "lambda.min", type = "response"))
roc_train <- roc(y_train, pred_train, quiet = TRUE)
cat(sprintf("\n  Apparent AUC (training): %.4f\n", as.numeric(auc(roc_train))))

# ── 5. VALIDATION: ALL WALES PASS REGISTRY ──────────────────────────────────

cat("\n=== VALIDATION: All Wales PASS Registry ===\n")

pass <- as.data.frame(read_sav(file.path(data_dir, "PASS_wrong_dob.sav")))
cat(sprintf("  PASS raw: %d rows x %d cols\n", nrow(pass), ncol(pass)))

# Outcome: Positive1 (Yes/No/NA)
pass$fh <- ifelse(pass$Positive1 == "Yes", 1,
            ifelse(pass$Positive1 == "No", 0, NA))
cat(sprintf("  FH+: %d, FH-: %d, Unknown: %d\n",
            sum(pass$fh == 1, na.rm = TRUE),
            sum(pass$fh == 0, na.rm = TRUE),
            sum(is.na(pass$fh))))

# Index effect: I_Vs_R (1=Index, 2=Relative)
pass$index_effect <- ifelse(safe_num(pass$I_Vs_R) == 1, 1, 0)
pass$index_effect[is.na(pass$index_effect)] <- 0
cat(sprintf("  Index: %d, Relative: %d\n",
            sum(pass$index_effect == 1), sum(pass$index_effect == 0)))

# Demographics
pass$age <- safe_num(pass$BMI_AGE)
pass$sex <- ifelse(pass$Gender == "M", 1, ifelse(pass$Gender == "F", 0, NA))

# Lipids
pass$hdl   <- safe_num(pass$HDL.1)
pass$tg    <- safe_num(pass$TRG.1)
pass$ldl_m <- safe_num(pass$LDL.1)

# Statin correction
pass$red_factor <- sapply(pass$Treatment1.1, get_reduction)
pass$ldl_ut <- reverse_ldl(pass$ldl_m, pass$red_factor)

# Trig Filter
pass$trig_filter <- ifelse(!is.na(pass$ldl_ut) & !is.na(pass$tg) & pass$tg >= 0,
                            pmin(pass$ldl_ut / (pass$tg + 0.1), 50), NA)

# On statin flag
pass$on_statin <- ifelse(pass$red_factor > 0, 1, 0)

# Clinical signs
pass$tendon_xanth <- ifelse(safe_num(pass$TendonXanthomata) == 1, 1, 0)
pass$tendon_xanth[is.na(pass$tendon_xanth)] <- 0

pass$corneal_arcus <- ifelse(safe_num(pass$CornealArcus) == 1, 1, 0)
pass$corneal_arcus[is.na(pass$corneal_arcus)] <- 0

# DLCN score
pass$dlcn <- safe_num(pass$GenoTypingScore)

# Complete cases
pass_valid <- pass[!is.na(pass$fh) & complete.cases(pass[, features]), ]
cat(sprintf("  PASS validation: %d complete (FH+=%d, FH-=%d)\n",
            nrow(pass_valid), sum(pass_valid$fh == 1), sum(pass_valid$fh == 0)))

# ── 6. PREDICT & EVALUATE ───────────────────────────────────────────────────

cat("\n=== RESULTS ===\n")

X_pass <- as.matrix(pass_valid[, features])
pred_pass <- as.numeric(predict(cv_fit, newx = X_pass, s = "lambda.min", type = "response"))

roc_pass <- roc(pass_valid$fh, pred_pass, quiet = TRUE)
ci_pass  <- ci.auc(roc_pass, method = "bootstrap", boot.n = 2000, quiet = TRUE)

cat(sprintf("\n  *** TUDOR AUC (PASS): %.3f (95%% CI: %.3f-%.3f) ***\n",
            as.numeric(auc(roc_pass)), ci_pass[1], ci_pass[3]))

# DLCN comparison (on same patients)
pass_dlcn <- pass_valid[!is.na(pass_valid$dlcn), ]
if (nrow(pass_dlcn) > 100) {
  X_dlcn <- as.matrix(pass_dlcn[, features])
  pred_dlcn_tudor <- as.numeric(predict(cv_fit, newx = X_dlcn, s = "lambda.min", type = "response"))

  roc_tudor_matched <- roc(pass_dlcn$fh, pred_dlcn_tudor, quiet = TRUE)
  roc_dlcn_matched  <- roc(pass_dlcn$fh, pass_dlcn$dlcn, quiet = TRUE)

  cat(sprintf("  TUDOR AUC (matched): %.3f\n", as.numeric(auc(roc_tudor_matched))))
  cat(sprintf("  DLCN AUC (matched):  %.3f\n", as.numeric(auc(roc_dlcn_matched))))

  delong <- roc.test(roc_tudor_matched, roc_dlcn_matched, method = "delong")
  cat(sprintf("  DeLong p-value: %s\n", formatC(delong$p.value, format = "e", digits = 2)))
}

# LDL-C alone
roc_ldl <- roc(pass_valid$fh, pass_valid$ldl_ut, quiet = TRUE)
cat(sprintf("  LDL-C alone AUC: %.3f\n", as.numeric(auc(roc_ldl))))

# Trig Filter alone
roc_tf <- roc(pass_valid$fh, pass_valid$trig_filter, quiet = TRUE)
cat(sprintf("  Trig_Filter alone AUC: %.3f\n", as.numeric(auc(roc_tf))))

# ── 7. CLINICAL PERFORMANCE ─────────────────────────────────────────────────

cat("\n=== CLINICAL PERFORMANCE ===\n")

coords_best <- coords(roc_pass, "best", best.method = "youden")
cat(sprintf("  Youden threshold: %.3f\n", coords_best$threshold))
cat(sprintf("  Sensitivity: %.1f%%\n", coords_best$sensitivity * 100))
cat(sprintf("  Specificity: %.1f%%\n", coords_best$specificity * 100))

# By patient type
for (ie in c(1, 0)) {
  label <- ifelse(ie == 1, "Index", "Cascade")
  sub <- pass_valid[pass_valid$index_effect == ie, ]
  if (nrow(sub) > 50 && sum(sub$fh == 1) > 10) {
    X_sub <- as.matrix(sub[, features])
    pred_sub <- as.numeric(predict(cv_fit, newx = X_sub, s = "lambda.min", type = "response"))
    roc_sub <- roc(sub$fh, pred_sub, quiet = TRUE)
    ci_sub  <- ci.auc(roc_sub, quiet = TRUE)
    cat(sprintf("  %s patients: AUC=%.3f (%.3f-%.3f) [n=%d, FH+=%d]\n",
                label, as.numeric(auc(roc_sub)), ci_sub[1], ci_sub[3],
                nrow(sub), sum(sub$fh == 1)))
  }
}

# ── 8. LOCO-CV (Leave-One-Cohort-Out) ───────────────────────────────────────

cat("\n=== LOCO-CV: LEAVE-ONE-COHORT-OUT ===\n")
cat("  Fold 1: Train on PASS → Validate on South Wales (pos+neg)\n")
cat("  Fold 2: Train on South Wales → Validate on PASS\n\n")

# Fold 1: Train PASS, validate South Wales
pass_train <- pass[!is.na(pass$fh) & complete.cases(pass[, features]), ]
sw_test    <- train_complete  # South Wales = mutation pos + neg

X_pass_tr <- as.matrix(pass_train[, features])
y_pass_tr <- pass_train$fh

cv_fold1 <- cv.glmnet(X_pass_tr, y_pass_tr, family = "binomial",
                       alpha = 0.5, nfolds = 10, type.measure = "auc",
                       standardize = TRUE)

X_sw_te <- as.matrix(sw_test[, features])
pred_fold1 <- as.numeric(predict(cv_fold1, newx = X_sw_te, s = "lambda.min", type = "response"))
roc_fold1 <- roc(sw_test$fh, pred_fold1, quiet = TRUE)
ci_fold1  <- ci.auc(roc_fold1, quiet = TRUE)

cat(sprintf("  Fold 1 (Train=PASS → Test=SouthWales): AUC=%.3f (%.3f-%.3f)\n",
            as.numeric(auc(roc_fold1)), ci_fold1[1], ci_fold1[3]))

# Fold 2: Train South Wales, validate PASS (this is the main result)
cat(sprintf("  Fold 2 (Train=SouthWales → Test=PASS):  AUC=%.3f (%.3f-%.3f)\n",
            as.numeric(auc(roc_pass)), ci_pass[1], ci_pass[3]))

# Average
avg_auc <- (as.numeric(auc(roc_fold1)) + as.numeric(auc(roc_pass))) / 2
cat(sprintf("\n  LOCO-CV Average AUC: %.3f\n", avg_auc))

# ── 9. SUMMARY ──────────────────────────────────────────────────────────────

cat("\n================================================================\n")
cat("  TUDOR v2 DIAGNOSTIC MODEL — FINAL RESULTS\n")
cat("================================================================\n")
cat(sprintf("  Training:  South Wales (FH+=%d, FH-=%d, Total=%d)\n",
            sum(train_complete$fh == 1), sum(train_complete$fh == 0), nrow(train_complete)))
cat(sprintf("  Validation: All Wales PASS (FH+=%d, FH-=%d, Total=%d)\n",
            sum(pass_valid$fh == 1), sum(pass_valid$fh == 0), nrow(pass_valid)))
cat("  ────────────────────────────────────────────────\n")
cat(sprintf("  TUDOR AUC (PASS):   %.3f (95%% CI: %.3f-%.3f)\n",
            as.numeric(auc(roc_pass)), ci_pass[1], ci_pass[3]))
cat(sprintf("  LOCO-CV Avg AUC:    %.3f\n", avg_auc))
cat(sprintf("  DLCN AUC (matched): %.3f\n",
            ifelse(exists("roc_dlcn_matched"), as.numeric(auc(roc_dlcn_matched)), NA)))
cat("================================================================\n")

# Save
results_df <- data.frame(
  Metric = c("TUDOR_PASS_AUC", "TUDOR_PASS_CI_low", "TUDOR_PASS_CI_high",
             "LOCO_Fold1_AUC", "LOCO_Fold2_AUC", "LOCO_Average",
             "DLCN_AUC", "LDL_alone_AUC", "TrigFilter_alone_AUC",
             "Sensitivity_Youden", "Specificity_Youden",
             "N_train", "N_valid", "N_FH_train", "N_FH_valid"),
  Value = c(as.numeric(auc(roc_pass)), ci_pass[1], ci_pass[3],
            as.numeric(auc(roc_fold1)), as.numeric(auc(roc_pass)), avg_auc,
            ifelse(exists("roc_dlcn_matched"), as.numeric(auc(roc_dlcn_matched)), NA),
            as.numeric(auc(roc_ldl)), as.numeric(auc(roc_tf)),
            coords_best$sensitivity, coords_best$specificity,
            nrow(train_complete), nrow(pass_valid),
            sum(train_complete$fh == 1), sum(pass_valid$fh == 1))
)
write.csv(results_df, file.path(data_dir, "TUDOR_reproduced_AUC.csv"), row.names = FALSE)
write.csv(coef_df, file.path(data_dir, "TUDOR_reproduced_coefficients.csv"), row.names = FALSE)
saveRDS(cv_fit, file.path(data_dir, "TUDOR_v2_model.rds"))

cat(sprintf("\n  Saved: TUDOR_reproduced_AUC.csv\n"))
cat(sprintf("  Saved: TUDOR_reproduced_coefficients.csv\n"))
cat(sprintf("  Saved: TUDOR_v2_model.rds\n"))
cat(sprintf("  Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

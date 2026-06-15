################################################################################
#
#  TUDOR v2: UKB LIPID CLINIC VALIDATION + FULL LANCET STATISTICS
#
#  Lipid clinic referral criteria (NICE CG71 + clinical practice):
#    TC > 7.5 mmol/L                                         OR
#    LDL-C > 4.9 mmol/L (pre-treatment)                     OR
#    Non-HDL-C > 5.9 mmol/L                                  OR
#    Premature ASCVD in participant (MI/stroke/angina < 60)  OR
#    Family history of FH (relative with high cholesterol)   OR
#    Family history of premature CVD (relative with heart    OR
#      disease/stroke; age-of-onset not available in UKB)    OR
#    On lipid-lowering therapy (proxy for clinical concern)
#
#  Full Lancet statistics:
#    - ROC/AUC with bootstrap CI
#    - DeLong pairwise comparisons
#    - Calibration (slope, intercept, Hosmer-Lemeshow)
#    - NRI (categorical + continuous)
#    - IDI
#    - Decision Curve Analysis
#    - Brier score
#    - Subgroup analyses (sex, age, gene, statin status)
#    - Reclassification tables
#
#  Author: Dr Nader Genedy | 2026-03-27
#
################################################################################

cat("\n================================================================\n")
cat("  TUDOR v2: UKB LIPID CLINIC + LANCET STATISTICS\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

set.seed(42)

suppressPackageStartupMessages({
  library(data.table); library(pROC); library(glmnet); library(boot)
})

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

DATA_DIR <- "C:/Users/nader/Downloads/calon_ukb_pipeline"
OUT_DIR  <- file.path(DATA_DIR, "tudor_loco_output")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load trained model
cv_fit <- readRDS(file.path(DATA_DIR, "TUDOR_v2_model.rds"))
features <- c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex",
              "tendon_xanth", "corneal_arcus", "on_statin")

# ══════════════════════════════════════════════════════════════════════════════
# PART 1: LOAD UKB DATA + LIPID CLINIC FILTER
# ══════════════════════════════════════════════════════════════════════════════

cat("=== LOADING UKB DATA ===\n")

# Try multiple paths
ukb_paths <- c(
  "D:/alphafold_backup/tudor_packup/TUDOR_UKB_Features.csv",
  file.path(DATA_DIR, "TUDOR_UKB_Features (1).csv"),
  "D:/calon_backup_data/TUDOR_UKB_Features (1).csv"
)
ukb_file <- NULL
for (p in ukb_paths) if (file.exists(p)) { ukb_file <- p; break }
if (is.null(ukb_file)) stop("UKB features file not found")

ukb <- as.data.frame(fread(ukb_file, showProgress = FALSE))
cat(sprintf("  UKB raw: %d rows\n", nrow(ukb)))

# ── Extract key variables ──
ukb$fh     <- safe_num(ukb$is_fh_genetic); ukb$fh[is.na(ukb$fh)] <- 0
ukb$age    <- safe_num(ukb$Age_at_LDL1)
ukb$sex    <- safe_num(ukb$Gender_num)
ukb$hdl    <- safe_num(ukb$HDL.1)
ukb$tg     <- safe_num(ukb$TRG.1)
ukb$ldl_m  <- safe_num(ukb$LDL_treated)
ukb$ldl_ut <- safe_num(ukb$LDL_untreated)
ukb$tc     <- safe_num(ukb$CHOL)
ukb$rf     <- safe_num(ukb$reduction_factor); ukb$rf[is.na(ukb$rf)] <- 0
ukb$on_statin <- ifelse(ukb$rf > 0, 1, 0)
ukb$non_hdl <- ifelse(!is.na(ukb$tc) & !is.na(ukb$hdl), ukb$tc - ukb$hdl, NA)
ukb$trig_filter <- safe_num(ukb$Trig_Filter)

# No clinical signs in UKB
ukb$tendon_xanth  <- 0
ukb$corneal_arcus <- 0

# Gene
ukb$gene <- as.character(ukb$gene)

# ── ASCVD: self-reported illness codes (p20002) ──────────────────────────────
# Codes: 1075=heart attack, 1074=angina, 1081=stroke/TIA, 1583=peripheral vasc
ascvd_codes <- c(1075, 1074, 1081, 1583)
p20002_cols <- grep("^participant\\.p20002", names(ukb), value = TRUE)

if (length(p20002_cols) > 0) {
  ukb$has_ascvd_p20002 <- apply(ukb[, p20002_cols], 1, function(row) {
    as.integer(any(row %in% ascvd_codes, na.rm = TRUE))
  })
  cat(sprintf("  ASCVD from self-reported illness (p20002): %d\n",
              sum(ukb$has_ascvd_p20002)))
} else {
  ukb$has_ascvd_p20002 <- 0
  cat("  WARNING: No p20002 columns — ASCVD (illness) not available\n")
}

# ── ASCVD: vascular/heart problems diagnosed by doctor (p6150) ───────────────
# Codes: 1=heart attack, 2=angina, 3=stroke, 4=high blood pressure
p6150_cols <- grep("^participant\\.p6150", names(ukb), value = TRUE)
if (length(p6150_cols) > 0) {
  ukb$has_ascvd_p6150 <- apply(ukb[, p6150_cols, drop = FALSE], 1, function(row) {
    as.integer(any(row %in% c(1, 2, 3), na.rm = TRUE))
  })
  cat(sprintf("  ASCVD from diagnosed conditions (p6150):   %d\n",
              sum(ukb$has_ascvd_p6150)))
} else {
  ukb$has_ascvd_p6150 <- 0
  cat("  WARNING: No p6150 columns — ASCVD (diagnosed) not available\n")
}

# Combined ASCVD: either source
ukb$has_ascvd <- as.integer(ukb$has_ascvd_p20002 == 1 | ukb$has_ascvd_p6150 == 1)
cat(sprintf("  ASCVD combined (p20002 OR p6150):          %d\n", sum(ukb$has_ascvd)))

# ── FAMILY HISTORY parsing (p20107=father, p20110=mother, p20111=siblings) ───
# Fields stored as concatenated strings e.g. "[1,12,-27]"
# Code 1  = heart disease  (proxy for CVD family history)
# Code 12 = high cholesterol (proxy for FH family history)
# Code 13 = stroke
parse_fh_code <- function(df, cols, code) {
  pattern <- paste0("(^|,|\\[)", code, "($|,|\\])")
  if (length(cols) == 0) return(rep(0L, nrow(df)))
  apply(df[, cols, drop = FALSE], 1, function(row) {
    as.integer(any(grepl(pattern, as.character(row)), na.rm = TRUE))
  })
}

fh_cols     <- grep("^participant\\.p20107|^participant\\.p20110|^participant\\.p20111",
                    names(ukb), value = TRUE)
father_cols  <- grep("^participant\\.p20107", names(ukb), value = TRUE)
mother_cols  <- grep("^participant\\.p20110", names(ukb), value = TRUE)
sibling_cols <- grep("^participant\\.p20111", names(ukb), value = TRUE)

# Family history of high cholesterol (FH indicator — code 12)
ukb$fhx_high_chol <- as.integer(
  parse_fh_code(ukb, father_cols,  12) == 1 |
  parse_fh_code(ukb, mother_cols,  12) == 1 |
  parse_fh_code(ukb, sibling_cols, 12) == 1
)
# Family history of heart disease or stroke (premature CVD proxy — codes 1, 13)
ukb$fhx_cvd <- as.integer(
  parse_fh_code(ukb, father_cols,  1)  == 1 | parse_fh_code(ukb, father_cols,  13) == 1 |
  parse_fh_code(ukb, mother_cols,  1)  == 1 | parse_fh_code(ukb, mother_cols,  13) == 1 |
  parse_fh_code(ukb, sibling_cols, 1)  == 1 | parse_fh_code(ukb, sibling_cols, 13) == 1
)
cat(sprintf("  Family history of high cholesterol (code 12): %d\n",
            sum(ukb$fhx_high_chol)))
cat(sprintf("  Family history of heart disease/stroke (1,13):%d\n",
            sum(ukb$fhx_cvd)))

# ── LIPID CLINIC REFERRAL FILTER ──
# NICE CG71 + clinical practice criteria

ukb$lc_tc       <- !is.na(ukb$tc) & ukb$tc > 7.5
ukb$lc_ldl      <- !is.na(ukb$ldl_ut) & ukb$ldl_ut > 4.9
ukb$lc_nhdl     <- !is.na(ukb$non_hdl) & ukb$non_hdl > 5.9
ukb$lc_ascvd    <- ukb$has_ascvd == 1 & !is.na(ukb$age) & ukb$age < 60
ukb$lc_fhx_fh   <- ukb$fhx_high_chol == 1            # FH of hypercholesterolaemia
ukb$lc_fhx_cvd  <- ukb$fhx_cvd == 1                  # FH of heart disease / stroke
ukb$lc_statin   <- ukb$on_statin == 1                 # On lipid-lowering therapy

ukb$lipid_clinic <- ukb$lc_tc | ukb$lc_ldl | ukb$lc_nhdl |
                    ukb$lc_ascvd | ukb$lc_fhx_fh | ukb$lc_fhx_cvd | ukb$lc_statin

# Label each participant's primary referral reason (for subgroup analysis)
ukb$ref_reason <- ifelse(ukb$lc_ldl,     "LDL>4.9",
                  ifelse(ukb$lc_tc,      "TC>7.5",
                  ifelse(ukb$lc_nhdl,    "nonHDL>5.9",
                  ifelse(ukb$lc_ascvd,   "Premature ASCVD",
                  ifelse(ukb$lc_fhx_fh,  "FHx of FH",
                  ifelse(ukb$lc_fhx_cvd, "FHx of CVD",
                  ifelse(ukb$lc_statin,  "On statin", NA)))))))

# Apply filter
ukb_lc <- ukb[ukb$lipid_clinic == TRUE, ]

# ── Count how many each criterion contributes (any criterion, not exclusive) ──
n_total   <- nrow(ukb)
n_lc      <- nrow(ukb_lc)
fh_prev   <- 100 * mean(ukb_lc$fh == 1)

# Unique contributors: how many would be ADDED by each criterion over the others
n_lc_no_tc      <- sum(!ukb$lc_tc & (ukb$lc_ldl | ukb$lc_nhdl | ukb$lc_ascvd | ukb$lc_fhx_fh | ukb$lc_fhx_cvd | ukb$lc_statin))
added_by_tc     <- n_lc - sum(ukb$lc_ldl | ukb$lc_nhdl | ukb$lc_ascvd | ukb$lc_fhx_fh | ukb$lc_fhx_cvd | ukb$lc_statin)
added_by_ldl    <- n_lc - sum(ukb$lc_tc  | ukb$lc_nhdl | ukb$lc_ascvd | ukb$lc_fhx_fh | ukb$lc_fhx_cvd | ukb$lc_statin)
added_by_nhdl   <- n_lc - sum(ukb$lc_tc  | ukb$lc_ldl  | ukb$lc_ascvd | ukb$lc_fhx_fh | ukb$lc_fhx_cvd | ukb$lc_statin)
added_by_ascvd  <- n_lc - sum(ukb$lc_tc  | ukb$lc_ldl  | ukb$lc_nhdl  | ukb$lc_fhx_fh | ukb$lc_fhx_cvd | ukb$lc_statin)
added_by_fhx_fh <- n_lc - sum(ukb$lc_tc  | ukb$lc_ldl  | ukb$lc_nhdl  | ukb$lc_ascvd  | ukb$lc_fhx_cvd | ukb$lc_statin)
added_by_fhx_cv <- n_lc - sum(ukb$lc_tc  | ukb$lc_ldl  | ukb$lc_nhdl  | ukb$lc_ascvd  | ukb$lc_fhx_fh  | ukb$lc_statin)
added_by_statin <- n_lc - sum(ukb$lc_tc  | ukb$lc_ldl  | ukb$lc_nhdl  | ukb$lc_ascvd  | ukb$lc_fhx_fh  | ukb$lc_fhx_cvd)

fh_in_added <- function(flag_new, flag_base) {
  only_new <- ukb$fh == 1 & flag_new & !flag_base
  c(n = sum(flag_new & !flag_base), fh = sum(only_new))
}
base_no_tc     <- ukb$lc_ldl | ukb$lc_nhdl | ukb$lc_ascvd | ukb$lc_fhx_fh | ukb$lc_fhx_cvd | ukb$lc_statin
base_no_ldl    <- ukb$lc_tc  | ukb$lc_nhdl | ukb$lc_ascvd | ukb$lc_fhx_fh | ukb$lc_fhx_cvd | ukb$lc_statin
base_no_nhdl   <- ukb$lc_tc  | ukb$lc_ldl  | ukb$lc_ascvd | ukb$lc_fhx_fh | ukb$lc_fhx_cvd | ukb$lc_statin
base_no_ascvd  <- ukb$lc_tc  | ukb$lc_ldl  | ukb$lc_nhdl  | ukb$lc_fhx_fh | ukb$lc_fhx_cvd | ukb$lc_statin
base_no_fhxfh  <- ukb$lc_tc  | ukb$lc_ldl  | ukb$lc_nhdl  | ukb$lc_ascvd  | ukb$lc_fhx_cvd | ukb$lc_statin
base_no_fhxcv  <- ukb$lc_tc  | ukb$lc_ldl  | ukb$lc_nhdl  | ukb$lc_ascvd  | ukb$lc_fhx_fh  | ukb$lc_statin
base_no_stat   <- ukb$lc_tc  | ukb$lc_ldl  | ukb$lc_nhdl  | ukb$lc_ascvd  | ukb$lc_fhx_fh  | ukb$lc_fhx_cvd

cat(sprintf("\n  ╔══════════════════════════════════════════════════════╗\n"))
cat(sprintf("  ║  EXPANDED LIPID CLINIC REFERRAL FILTER SUMMARY      ║\n"))
cat(sprintf("  ╠══════════════════════════════════════════════════════╣\n"))
cat(sprintf("  ║  Criterion              N (any)  Unique add  FH+ add ║\n"))
cat(sprintf("  ╠══════════════════════════════════════════════════════╣\n"))

print_row <- function(label, flag_any, flag_base, n_any, n_fh_any) {
  unique_n  <- sum(flag_any & !flag_base, na.rm = TRUE)
  unique_fh <- sum(ukb$fh == 1 & flag_any & !flag_base, na.rm = TRUE)
  cat(sprintf("  ║  %-22s  %7d   %7d    %5d  ║\n",
              label, n_any, unique_n, unique_fh))
}

print_row("TC > 7.5",         ukb$lc_tc,      base_no_tc,    sum(ukb$lc_tc),     sum(ukb$fh==1 & ukb$lc_tc))
print_row("LDL-C > 4.9",      ukb$lc_ldl,     base_no_ldl,   sum(ukb$lc_ldl),    sum(ukb$fh==1 & ukb$lc_ldl))
print_row("non-HDL > 5.9",    ukb$lc_nhdl,    base_no_nhdl,  sum(ukb$lc_nhdl),   sum(ukb$fh==1 & ukb$lc_nhdl))
print_row("Premature ASCVD<60",ukb$lc_ascvd,  base_no_ascvd, sum(ukb$lc_ascvd),  sum(ukb$fh==1 & ukb$lc_ascvd))
print_row("FHx high chol (12)",ukb$lc_fhx_fh, base_no_fhxfh, sum(ukb$lc_fhx_fh), sum(ukb$fh==1 & ukb$lc_fhx_fh))
print_row("FHx CVD (1/13)",   ukb$lc_fhx_cvd, base_no_fhxcv, sum(ukb$lc_fhx_cvd),sum(ukb$fh==1 & ukb$lc_fhx_cvd))
print_row("On statin",        ukb$lc_statin,  base_no_stat,  sum(ukb$lc_statin), sum(ukb$fh==1 & ukb$lc_statin))

cat(sprintf("  ╠══════════════════════════════════════════════════════╣\n"))
cat(sprintf("  ║  TOTAL lipid clinic     %7d       —       %5d  ║\n",
            n_lc, sum(ukb_lc$fh == 1)))
cat(sprintf("  ║  pct of UKB             %6.1f%%                       ║\n",
            100 * n_lc / n_total))
cat(sprintf("  ║  FH prevalence in LC    %6.2f%%                       ║\n",
            fh_prev))
cat(sprintf("  ╚══════════════════════════════════════════════════════╝\n\n"))

# Complete cases (lc_tc, lc_ldl, has_ascvd already exist in ukb_lc from ukb)
ukb_valid <- ukb_lc[complete.cases(ukb_lc[, c("fh", features)]), ]
cat(sprintf("    Complete cases:        %d (FH+=%d, FH-=%d)\n",
            nrow(ukb_valid), sum(ukb_valid$fh == 1), sum(ukb_valid$fh == 0)))

# ── eDLCN score (truncated — no family history/clinical signs in UKB) ──
ukb_valid$dlcn <- ifelse(!is.na(ukb_valid$ldl_ut) & ukb_valid$ldl_ut >= 8.5, 8,
                  ifelse(!is.na(ukb_valid$ldl_ut) & ukb_valid$ldl_ut >= 6.5, 5,
                  ifelse(!is.na(ukb_valid$ldl_ut) & ukb_valid$ldl_ut >= 5.0, 3,
                  ifelse(!is.na(ukb_valid$ldl_ut) & ukb_valid$ldl_ut >= 4.0, 1, 0))))
# Add premature ASCVD points
ukb_valid$dlcn <- ukb_valid$dlcn +
  ifelse(ukb_valid$has_ascvd == 1 & ukb_valid$age < 55 & ukb_valid$sex == 1, 2,
  ifelse(ukb_valid$has_ascvd == 1 & ukb_valid$age < 60 & ukb_valid$sex == 0, 2, 0))

# ══════════════════════════════════════════════════════════════════════════════
# PART 2: TUDOR PREDICTIONS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== TUDOR PREDICTIONS ON LIPID CLINIC COHORT ===\n")

X_ukb <- as.matrix(ukb_valid[, features])
ukb_valid$tudor_prob <- as.numeric(predict(cv_fit, newx = X_ukb,
                                           s = "lambda.min", type = "response"))

# ══════════════════════════════════════════════════════════════════════════════
# PART 3: DISCRIMINATION
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== DISCRIMINATION ===\n")

# TUDOR
roc_tudor <- roc(ukb_valid$fh, ukb_valid$tudor_prob, quiet = TRUE)
ci_tudor  <- ci.auc(roc_tudor, method = "bootstrap", boot.n = 2000, quiet = TRUE)

# eDLCN
roc_dlcn <- roc(ukb_valid$fh, ukb_valid$dlcn, quiet = TRUE)
ci_dlcn  <- ci.auc(roc_dlcn, quiet = TRUE)

# LDL-C alone
roc_ldl <- roc(ukb_valid$fh, ukb_valid$ldl_ut, quiet = TRUE)
ci_ldl  <- ci.auc(roc_ldl, quiet = TRUE)

# Trig_Filter alone
roc_tf <- roc(ukb_valid$fh, ukb_valid$trig_filter, quiet = TRUE)
ci_tf  <- ci.auc(roc_tf, quiet = TRUE)

cat(sprintf("  TUDOR:       AUC = %.3f (%.3f-%.3f)\n",
            auc(roc_tudor), ci_tudor[1], ci_tudor[3]))
cat(sprintf("  eDLCN:       AUC = %.3f (%.3f-%.3f)\n",
            auc(roc_dlcn), ci_dlcn[1], ci_dlcn[3]))
cat(sprintf("  LDL-C alone: AUC = %.3f (%.3f-%.3f)\n",
            auc(roc_ldl), ci_ldl[1], ci_ldl[3]))
cat(sprintf("  Trig_Filter: AUC = %.3f (%.3f-%.3f)\n",
            auc(roc_tf), ci_tf[1], ci_tf[3]))

# DeLong pairwise
dt1 <- roc.test(roc_tudor, roc_dlcn, method = "delong")
dt2 <- roc.test(roc_tudor, roc_ldl, method = "delong")
dt3 <- roc.test(roc_tudor, roc_tf, method = "delong")
dt4 <- roc.test(roc_dlcn, roc_ldl, method = "delong")

cat(sprintf("\n  DeLong comparisons (paired, same patients):\n"))
cat(sprintf("    TUDOR vs eDLCN:       Z=%.2f, p=%s\n", dt1$statistic, formatC(dt1$p.value, format="e", digits=2)))
cat(sprintf("    TUDOR vs LDL-C:       Z=%.2f, p=%s\n", dt2$statistic, formatC(dt2$p.value, format="e", digits=2)))
cat(sprintf("    TUDOR vs Trig_Filter: Z=%.2f, p=%s\n", dt3$statistic, formatC(dt3$p.value, format="e", digits=2)))
cat(sprintf("    eDLCN vs LDL-C:       Z=%.2f, p=%s\n", dt4$statistic, formatC(dt4$p.value, format="e", digits=2)))

# ══════════════════════════════════════════════════════════════════════════════
# PART 4: CLINICAL PERFORMANCE AT YOUDEN THRESHOLD
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== CLINICAL PERFORMANCE ===\n")

coords_tudor <- coords(roc_tudor, "best", best.method = "youden")
prev <- mean(ukb_valid$fh)
sens <- coords_tudor$sensitivity; spec <- coords_tudor$specificity
ppv <- (sens * prev) / (sens * prev + (1 - spec) * (1 - prev))
npv <- (spec * (1 - prev)) / ((1 - sens) * prev + spec * (1 - prev))
lr_pos <- sens / (1 - spec)
lr_neg <- (1 - sens) / spec

cat(sprintf("  Youden threshold: %.3f\n", coords_tudor$threshold))
cat(sprintf("  Sensitivity: %.1f%%\n", sens * 100))
cat(sprintf("  Specificity: %.1f%%\n", spec * 100))
cat(sprintf("  PPV:         %.1f%%\n", ppv * 100))
cat(sprintf("  NPV:         %.1f%%\n", npv * 100))
cat(sprintf("  LR+:         %.2f\n", lr_pos))
cat(sprintf("  LR-:         %.4f\n", lr_neg))
cat(sprintf("  Prevalence:  %.2f%%\n", prev * 100))

# eDLCN at >=6 (probable/definite FH)
dlcn_pos <- ukb_valid$dlcn >= 6
dlcn_sens <- sum(dlcn_pos & ukb_valid$fh == 1) / sum(ukb_valid$fh == 1)
dlcn_spec <- sum(!dlcn_pos & ukb_valid$fh == 0) / sum(ukb_valid$fh == 0)
dlcn_ppv <- sum(dlcn_pos & ukb_valid$fh == 1) / sum(dlcn_pos)
dlcn_npv <- sum(!dlcn_pos & ukb_valid$fh == 0) / sum(!dlcn_pos)

cat(sprintf("\n  eDLCN at >=6 (probable/definite):\n"))
cat(sprintf("  Sensitivity: %.1f%%\n", dlcn_sens * 100))
cat(sprintf("  Specificity: %.1f%%\n", dlcn_spec * 100))
cat(sprintf("  PPV:         %.1f%%\n", dlcn_ppv * 100))
cat(sprintf("  NPV:         %.1f%%\n", dlcn_npv * 100))

# ══════════════════════════════════════════════════════════════════════════════
# PART 5: CALIBRATION
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== CALIBRATION ===\n")

# Calibration-in-the-large: intercept and slope
calib_glm <- glm(fh ~ tudor_prob, data = ukb_valid, family = binomial)
calib_intercept <- coef(calib_glm)[1]
calib_slope     <- coef(calib_glm)[2]

cat(sprintf("  Calibration intercept: %.3f (ideal=0)\n", calib_intercept))
cat(sprintf("  Calibration slope:     %.3f (ideal=1)\n", calib_slope))

# Hosmer-Lemeshow (10 groups)
ukb_valid$prob_decile <- cut(ukb_valid$tudor_prob,
                              breaks = quantile(ukb_valid$tudor_prob,
                                                probs = seq(0, 1, 0.1), na.rm = TRUE),
                              include.lowest = TRUE, labels = FALSE)

hl_table <- aggregate(cbind(observed = fh, expected = tudor_prob) ~ prob_decile,
                       data = ukb_valid, FUN = sum)
hl_n <- aggregate(fh ~ prob_decile, data = ukb_valid, FUN = length)
names(hl_n)[2] <- "n"
hl_table <- merge(hl_table, hl_n, by = "prob_decile")

hl_chi2 <- sum((hl_table$observed - hl_table$expected)^2 /
               (hl_table$expected * (1 - hl_table$expected / hl_table$n) + 0.001))
hl_p <- pchisq(hl_chi2, df = nrow(hl_table) - 2, lower.tail = FALSE)

cat(sprintf("  Hosmer-Lemeshow: chi2=%.2f, df=%d, p=%.4f\n",
            hl_chi2, nrow(hl_table) - 2, hl_p))

# Brier score
brier <- mean((ukb_valid$tudor_prob - ukb_valid$fh)^2)
brier_max <- prev * (1 - prev)
brier_scaled <- 1 - brier / brier_max

cat(sprintf("  Brier score: %.4f (scaled=%.4f, max=%.4f)\n",
            brier, brier_scaled, brier_max))

# ── Recalibrated PPV/NPV at realistic prevalence scenarios ───────────────────
# The Wales-trained model is miscalibrated in UKB (prevalence mismatch).
# PPV/NPV from Bayes' theorem using AUC-derived Sens/Spec are prevalence-independent.
cat("\n  Recalibrated performance at clinical prevalence scenarios:\n")
cat(sprintf("  %-25s  %6s  %6s  %6s  %6s\n",
            "Setting (prevalence)", "Sens%", "Spec%", "PPV%", "NPV%"))
cat(sprintf("  %s\n", strrep("-", 62)))
coords_t <- coords(roc_tudor, "best", best.method = "youden")
sens_y   <- coords_t$sensitivity
spec_y   <- coords_t$specificity
for (prev_scenario in list(
  list(label = "Lipid clinic (1.26%)", prev = 0.0126),
  list(label = "Specialist FH clinic (10%)", prev = 0.10),
  list(label = "Cascade screening (33%)", prev = 0.33),
  list(label = "High-risk family (50%)", prev = 0.50)
)) {
  p  <- prev_scenario$prev
  ppv_s <- (sens_y * p) / (sens_y * p + (1 - spec_y) * (1 - p))
  npv_s <- (spec_y * (1 - p)) / ((1 - sens_y) * p + spec_y * (1 - p))
  cat(sprintf("  %-25s  %5.1f   %5.1f   %5.1f   %5.1f\n",
              prev_scenario$label, sens_y*100, spec_y*100, ppv_s*100, npv_s*100))
}
cat(sprintf("  %s\n\n", strrep("-", 62)))

# Calibration table
cat("\n  Calibration by decile:\n")
cat("  Decile   N     Obs    Exp    O/E\n")
for (i in seq_len(nrow(hl_table))) {
  oe <- hl_table$observed[i] / max(hl_table$expected[i], 0.01)
  cat(sprintf("    %2d   %5d   %4d  %6.1f  %.2f\n",
              hl_table$prob_decile[i], hl_table$n[i],
              hl_table$observed[i], hl_table$expected[i], oe))
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 6: RECLASSIFICATION (NRI + IDI)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== RECLASSIFICATION (TUDOR vs eDLCN) ===\n")

# Normalise eDLCN to probability scale
dlcn_glm <- glm(fh ~ dlcn, data = ukb_valid, family = binomial)
ukb_valid$dlcn_prob <- predict(dlcn_glm, type = "response")

# Categorical NRI at clinical thresholds: <1%, 1-3%, >3%
tudor_cat <- ifelse(ukb_valid$tudor_prob < 0.01, "Low",
             ifelse(ukb_valid$tudor_prob <= 0.03, "Intermediate", "High"))
dlcn_cat  <- ifelse(ukb_valid$dlcn_prob < 0.01, "Low",
             ifelse(ukb_valid$dlcn_prob <= 0.03, "Intermediate", "High"))

# Events (FH+)
ev <- ukb_valid$fh == 1
ev_up   <- sum(tudor_cat[ev] == "High" & dlcn_cat[ev] != "High") +
           sum(tudor_cat[ev] == "Intermediate" & dlcn_cat[ev] == "Low")
ev_down <- sum(tudor_cat[ev] == "Low" & dlcn_cat[ev] != "Low") +
           sum(tudor_cat[ev] == "Intermediate" & dlcn_cat[ev] == "High")
nri_ev  <- (ev_up - ev_down) / sum(ev)

# Non-events (FH-)
ne <- ukb_valid$fh == 0
ne_down <- sum(tudor_cat[ne] == "Low" & dlcn_cat[ne] != "Low") +
           sum(tudor_cat[ne] == "Intermediate" & dlcn_cat[ne] == "High")
ne_up   <- sum(tudor_cat[ne] == "High" & dlcn_cat[ne] != "High") +
           sum(tudor_cat[ne] == "Intermediate" & dlcn_cat[ne] == "Low")
nri_ne  <- (ne_down - ne_up) / sum(ne)
nri_total <- nri_ev + nri_ne

cat(sprintf("  Categorical NRI (<1%%, 1-3%%, >3%%):\n"))
cat(sprintf("    NRI(events):     %+.3f (up=%d, down=%d)\n", nri_ev, ev_up, ev_down))
cat(sprintf("    NRI(non-events): %+.3f (down=%d, up=%d)\n", nri_ne, ne_down, ne_up))
cat(sprintf("    NRI(total):      %+.3f\n", nri_total))

# Continuous NRI
p_tudor <- ukb_valid$tudor_prob
p_dlcn  <- ukb_valid$dlcn_prob
cnri_ev <- mean(p_tudor[ev] > p_dlcn[ev]) - mean(p_tudor[ev] < p_dlcn[ev])
cnri_ne <- mean(p_tudor[ne] < p_dlcn[ne]) - mean(p_tudor[ne] > p_dlcn[ne])
cnri_total <- cnri_ev + cnri_ne

cat(sprintf("\n  Continuous NRI:\n"))
cat(sprintf("    cNRI(events):     %+.3f\n", cnri_ev))
cat(sprintf("    cNRI(non-events): %+.3f\n", cnri_ne))
cat(sprintf("    cNRI(total):      %+.3f\n", cnri_total))

# IDI
idi <- (mean(p_tudor[ev]) - mean(p_tudor[ne])) -
       (mean(p_dlcn[ev]) - mean(p_dlcn[ne]))
cat(sprintf("\n  IDI: %+.4f\n", idi))

# Reclassification table
cat("\n  Reclassification matrix (FH+ events):\n")
cat("                   TUDOR-Low  TUDOR-Int  TUDOR-High\n")
for (d_cat in c("Low", "Intermediate", "High")) {
  mask_ev <- ev & dlcn_cat == d_cat
  cat(sprintf("    eDLCN-%-12s  %5d      %5d       %5d\n", d_cat,
              sum(tudor_cat[mask_ev] == "Low"),
              sum(tudor_cat[mask_ev] == "Intermediate"),
              sum(tudor_cat[mask_ev] == "High")))
}

cat("\n  Reclassification matrix (FH- non-events):\n")
cat("                   TUDOR-Low  TUDOR-Int  TUDOR-High\n")
for (d_cat in c("Low", "Intermediate", "High")) {
  mask_ne <- ne & dlcn_cat == d_cat
  cat(sprintf("    eDLCN-%-12s  %5d      %5d       %5d\n", d_cat,
              sum(tudor_cat[mask_ne] == "Low"),
              sum(tudor_cat[mask_ne] == "Intermediate"),
              sum(tudor_cat[mask_ne] == "High")))
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 7: DECISION CURVE ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== DECISION CURVE ANALYSIS ===\n")

thresholds <- seq(0.001, 0.10, 0.001)
dca_results <- data.frame()

for (pt in thresholds) {
  # TUDOR
  tp_t <- sum(ukb_valid$tudor_prob >= pt & ukb_valid$fh == 1)
  fp_t <- sum(ukb_valid$tudor_prob >= pt & ukb_valid$fh == 0)
  nb_tudor <- (tp_t / nrow(ukb_valid)) - (fp_t / nrow(ukb_valid)) * (pt / (1 - pt))

  # eDLCN (probability scale)
  tp_d <- sum(ukb_valid$dlcn_prob >= pt & ukb_valid$fh == 1)
  fp_d <- sum(ukb_valid$dlcn_prob >= pt & ukb_valid$fh == 0)
  nb_dlcn <- (tp_d / nrow(ukb_valid)) - (fp_d / nrow(ukb_valid)) * (pt / (1 - pt))

  # Treat all
  nb_all <- prev - (1 - prev) * (pt / (1 - pt))

  dca_results <- rbind(dca_results, data.frame(
    Threshold = pt, TUDOR = nb_tudor, eDLCN = nb_dlcn, Treat_All = nb_all
  ))
}

# Summary at key thresholds
cat("  Net benefit at clinical thresholds:\n")
cat("  Threshold   TUDOR     eDLCN     Treat-All  TUDOR advantage\n")
for (pt in c(0.005, 0.01, 0.02, 0.03, 0.05, 0.10)) {
  row <- dca_results[which.min(abs(dca_results$Threshold - pt)), ]
  cat(sprintf("    %.1f%%     %+.4f   %+.4f    %+.4f     %+.4f\n",
              pt * 100, row$TUDOR, row$eDLCN, row$Treat_All,
              row$TUDOR - row$eDLCN))
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 8: SUBGROUP ANALYSES
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== SUBGROUP ANALYSES ===\n")

subgroup_results <- list()

# Helper
compute_subgroup <- function(data, label) {
  if (nrow(data) < 50 || sum(data$fh == 1) < 5) return(NULL)
  X_s <- as.matrix(data[, features])
  p_s <- as.numeric(predict(cv_fit, newx = X_s, s = "lambda.min", type = "response"))
  r_s <- roc(data$fh, p_s, quiet = TRUE)
  ci_s <- tryCatch(ci.auc(r_s, quiet = TRUE), error = function(e) c(NA,NA,NA))
  data.frame(Subgroup = label, N = nrow(data), N_FH = sum(data$fh == 1),
             Prevalence = round(100 * mean(data$fh), 2),
             AUC = round(as.numeric(auc(r_s)), 3),
             CI_lo = round(ci_s[1], 3), CI_hi = round(ci_s[3], 3))
}

# Sex
subgroup_results[[1]] <- compute_subgroup(ukb_valid[ukb_valid$sex == 1, ], "Male")
subgroup_results[[2]] <- compute_subgroup(ukb_valid[ukb_valid$sex == 0, ], "Female")

# Age groups
subgroup_results[[3]] <- compute_subgroup(ukb_valid[ukb_valid$age < 40, ], "Age < 40")
subgroup_results[[4]] <- compute_subgroup(ukb_valid[ukb_valid$age >= 40 & ukb_valid$age < 55, ], "Age 40-54")
subgroup_results[[5]] <- compute_subgroup(ukb_valid[ukb_valid$age >= 55, ], "Age >= 55")

# Statin status
subgroup_results[[6]] <- compute_subgroup(ukb_valid[ukb_valid$on_statin == 1, ], "On statin")
subgroup_results[[7]] <- compute_subgroup(ukb_valid[ukb_valid$on_statin == 0, ], "Statin-free")

# ASCVD status
if ("has_ascvd" %in% names(ukb_valid) && sum(ukb_valid$has_ascvd == 1) > 50) {
  subgroup_results[[8]] <- compute_subgroup(ukb_valid[ukb_valid$has_ascvd == 1, ], "Prior ASCVD")
  subgroup_results[[9]] <- compute_subgroup(ukb_valid[ukb_valid$has_ascvd == 0, ], "No ASCVD")
}

# LDL range
subgroup_results[[10]] <- compute_subgroup(ukb_valid[ukb_valid$ldl_ut >= 4.0 & ukb_valid$ldl_ut < 5.0, ], "LDL 4.0-4.9")
subgroup_results[[11]] <- compute_subgroup(ukb_valid[ukb_valid$ldl_ut >= 5.0 & ukb_valid$ldl_ut < 6.5, ], "LDL 5.0-6.4")
subgroup_results[[12]] <- compute_subgroup(ukb_valid[ukb_valid$ldl_ut >= 6.5, ], "LDL >= 6.5")

# Gene-specific (FH+ only vs all FH-)
neg_data <- ukb_valid[ukb_valid$fh == 0, ]
for (g in c("LDLR", "APOB")) {
  fh_gene <- ukb_valid[ukb_valid$fh == 1 & ukb_valid$gene == g, ]
  if (nrow(fh_gene) >= 5) {
    combined <- rbind(fh_gene[, c("fh", features)], neg_data[, c("fh", features)])
    subgroup_results[[length(subgroup_results) + 1]] <-
      compute_subgroup(combined, paste0("Gene: ", g))
  }
}

# Referral reason subgroups (all criteria)
ukb_valid$ref_tc      <- !is.na(ukb_valid$tc)     & ukb_valid$tc > 7.5
ukb_valid$ref_ldl     <- !is.na(ukb_valid$ldl_ut) & ukb_valid$ldl_ut > 4.9
ukb_valid$ref_nhdl    <- !is.na(ukb_valid$non_hdl) & ukb_valid$non_hdl > 5.9
ukb_valid$ref_ascvd   <- if ("lc_ascvd"   %in% names(ukb_valid)) ukb_valid$lc_ascvd   else FALSE
ukb_valid$ref_fhx_fh  <- if ("lc_fhx_fh"  %in% names(ukb_valid)) ukb_valid$lc_fhx_fh  else FALSE
ukb_valid$ref_fhx_cvd <- if ("lc_fhx_cvd" %in% names(ukb_valid)) ukb_valid$lc_fhx_cvd else FALSE
ukb_valid$ref_statin  <- ukb_valid$on_statin == 1

for (ref in list(
  list(flag = "ref_tc",      label = "Referral: TC>7.5"),
  list(flag = "ref_ldl",     label = "Referral: LDL>4.9"),
  list(flag = "ref_nhdl",    label = "Referral: nonHDL>5.9"),
  list(flag = "ref_ascvd",   label = "Referral: Premature ASCVD"),
  list(flag = "ref_fhx_fh",  label = "Referral: FHx high chol"),
  list(flag = "ref_fhx_cvd", label = "Referral: FHx CVD"),
  list(flag = "ref_statin",  label = "Referral: On statin")
)) {
  mask <- ukb_valid[[ref$flag]] == TRUE
  if (sum(mask, na.rm = TRUE) >= 50)
    subgroup_results[[length(subgroup_results) + 1]] <-
      compute_subgroup(ukb_valid[mask & !is.na(mask), ], ref$label)
}

# FHx vs non-FHx: does TUDOR retain discrimination in those referred ONLY by family history?
if ("ref_fhx_fh" %in% names(ukb_valid)) {
  only_fhx <- ukb_valid$ref_fhx_fh == TRUE &
              !ukb_valid$ref_tc & !ukb_valid$ref_ldl & !ukb_valid$ref_nhdl
  if (sum(only_fhx, na.rm = TRUE) >= 50)
    subgroup_results[[length(subgroup_results) + 1]] <-
      compute_subgroup(ukb_valid[only_fhx & !is.na(only_fhx), ],
                       "FHx only (no lipid criteria)")
}

sg_df <- do.call(rbind, Filter(Negate(is.null), subgroup_results))
cat("\n")
print(sg_df)

# ══════════════════════════════════════════════════════════════════════════════
# PART 9: SAVE ALL RESULTS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== SAVING RESULTS ===\n")

# Discrimination table
disc_df <- data.frame(
  Model = c("TUDOR", "eDLCN", "LDL-C alone", "Trig_Filter"),
  AUC = c(auc(roc_tudor), auc(roc_dlcn), auc(roc_ldl), auc(roc_tf)),
  CI_lo = c(ci_tudor[1], ci_dlcn[1], ci_ldl[1], ci_tf[1]),
  CI_hi = c(ci_tudor[3], ci_dlcn[3], ci_ldl[3], ci_tf[3])
)
write.csv(disc_df, file.path(OUT_DIR, "ukb_lipid_clinic_discrimination.csv"), row.names = FALSE)

# DeLong table
delong_df <- data.frame(
  Comparison = c("TUDOR vs eDLCN", "TUDOR vs LDL-C", "TUDOR vs TrigFilter", "eDLCN vs LDL-C"),
  Z = c(dt1$statistic, dt2$statistic, dt3$statistic, dt4$statistic),
  P_value = c(dt1$p.value, dt2$p.value, dt3$p.value, dt4$p.value)
)
write.csv(delong_df, file.path(OUT_DIR, "ukb_lipid_clinic_delong.csv"), row.names = FALSE)

# Calibration
calib_df <- data.frame(
  Metric = c("Intercept", "Slope", "HL_chi2", "HL_p", "Brier", "Brier_scaled"),
  Value = c(calib_intercept, calib_slope, hl_chi2, hl_p, brier, brier_scaled)
)
write.csv(calib_df, file.path(OUT_DIR, "ukb_lipid_clinic_calibration.csv"), row.names = FALSE)

# NRI/IDI
nri_df <- data.frame(
  Metric = c("NRI_events", "NRI_nonevents", "NRI_total",
             "cNRI_events", "cNRI_nonevents", "cNRI_total", "IDI"),
  Value = c(nri_ev, nri_ne, nri_total, cnri_ev, cnri_ne, cnri_total, idi)
)
write.csv(nri_df, file.path(OUT_DIR, "ukb_lipid_clinic_nri_idi.csv"), row.names = FALSE)

# DCA
write.csv(dca_results, file.path(OUT_DIR, "ukb_lipid_clinic_dca.csv"), row.names = FALSE)

# Subgroups
write.csv(sg_df, file.path(OUT_DIR, "ukb_lipid_clinic_subgroups.csv"), row.names = FALSE)

# Calibration table
write.csv(hl_table, file.path(OUT_DIR, "ukb_lipid_clinic_calibration_table.csv"), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# PART 10: FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("  TUDOR v2 — UKB LIPID CLINIC VALIDATION COMPLETE\n")
cat("════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Lipid clinic population: %d (FH prevalence %.2f%%)\n",
            nrow(ukb_valid), 100 * mean(ukb_valid$fh)))
cat("  ────────────────────────────────────────────────\n")
cat(sprintf("  TUDOR AUC:       %.3f (%.3f-%.3f)\n", auc(roc_tudor), ci_tudor[1], ci_tudor[3]))
cat(sprintf("  eDLCN AUC:       %.3f (%.3f-%.3f)\n", auc(roc_dlcn), ci_dlcn[1], ci_dlcn[3]))
cat(sprintf("  DeLong p:        %s\n", formatC(dt1$p.value, format="e", digits=2)))
cat("  ────────────────────────────────────────────────\n")
cat(sprintf("  Sensitivity:     %.1f%%  (eDLCN: %.1f%%)\n", sens*100, dlcn_sens*100))
cat(sprintf("  Specificity:     %.1f%%  (eDLCN: %.1f%%)\n", spec*100, dlcn_spec*100))
cat(sprintf("  PPV:             %.1f%%  (eDLCN: %.1f%%)\n", ppv*100, dlcn_ppv*100))
cat(sprintf("  NPV:             %.1f%%  (eDLCN: %.1f%%)\n", npv*100, dlcn_npv*100))
cat("  ────────────────────────────────────────────────\n")
cat(sprintf("  Calibration:     intercept=%.3f, slope=%.3f\n", calib_intercept, calib_slope))
cat(sprintf("  Brier:           %.4f (scaled=%.4f)\n", brier, brier_scaled))
cat(sprintf("  NRI (total):     %+.3f\n", nri_total))
cat(sprintf("  cNRI (total):    %+.3f\n", cnri_total))
cat(sprintf("  IDI:             %+.4f\n", idi))
cat("════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

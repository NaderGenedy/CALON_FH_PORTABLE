################################################################################
#
#  TUDOR v2: FIXED WALES WEIGHTS — EXACT REPRODUCTION
#
#  Model: 6-feature logistic regression with FIXED coefficients from Wales
#
#  Score = 0.7557 + 0.0579×LDL_RW + 0.4924×Trig_Filter_RW
#          - 1.1280×HDL - 0.0334×Age - 0.0886×Sex
#
#  Prob = 1 / (1 + exp(-Score))
#
#  This is the EXACT model from 01_data_merge.R (TRIPOD Type 4)
#
#  Author: Dr Nader Genedy | 2026-03-27
#
################################################################################

cat("\n================================================================\n")
cat("  TUDOR v2: FIXED WALES WEIGHTS REPRODUCTION\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

set.seed(42)

suppressPackageStartupMessages({
  library(data.table); library(pROC); library(boot)
})

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# ══════════════════════════════════════════════════════════════════════════════
# FIXED TUDOR WEIGHTS (from Wales training — NEVER re-estimate)
# ══════════════════════════════════════════════════════════════════════════════

TUDOR <- list(
  intercept = 0.755722,
  beta_LDL  = 0.057911,
  beta_Trig = 0.492412,
  beta_HDL  = -1.128045,
  beta_Age  = -0.033393,
  beta_Sex  = -0.088550
)

tudor_score <- function(ldl_rw, trig_filter_rw, hdl, age, sex) {
  TUDOR$intercept +
    TUDOR$beta_LDL  * ldl_rw +
    TUDOR$beta_Trig * trig_filter_rw +
    TUDOR$beta_HDL  * hdl +
    TUDOR$beta_Age  * age +
    TUDOR$beta_Sex  * sex
}

tudor_prob <- function(score) 1 / (1 + exp(-score))

cat("  TUDOR weights:\n")
for (nm in names(TUDOR)) {
  cat(sprintf("    %-12s = %+.6f\n", nm, TUDOR[[nm]]))
}

# ══════════════════════════════════════════════════════════════════════════════
# STATIN CORRECTION (from original pipeline)
# ══════════════════════════════════════════════════════════════════════════════

statin_codes <- list(
  `1140888648` = list(name = "Atorvastatin", factor = 0.38),
  `1141146234` = list(name = "Atorvastatin", factor = 0.38),
  `1140888594` = list(name = "Simvastatin",  factor = 0.35),
  `1141192414` = list(name = "Rosuvastatin", factor = 0.34),
  `1140861958` = list(name = "Simvastatin",  factor = 0.35),
  `1140861970` = list(name = "Pravastatin",  factor = 0.23),
  `1140881748` = list(name = "Fluvastatin",  factor = 0.21),
  `1140910632` = list(name = "Ezetimibe",    factor = 0.18),
  `1141146138` = list(name = "Rosuvastatin", factor = 0.34),
  `1140888690` = list(name = "Pravastatin",  factor = 0.23),
  `1140861922` = list(name = "Lovastatin",   factor = 0.25)
)

# ══════════════════════════════════════════════════════════════════════════════
# LOAD UKB DATA
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== LOADING UKB DATA ===\n")

ukb_file <- "D:/alphafold_backup/tudor_packup/TUDOR_UKB_Features.csv"
ukb <- fread(ukb_file, showProgress = FALSE)
cat(sprintf("  UKB raw: %d rows\n", nrow(ukb)))

# Core fields
ukb$fh      <- as.integer(ukb$is_fh_genetic); ukb$fh[is.na(ukb$fh)] <- 0L
ukb$age     <- as.numeric(ukb$Age_at_LDL1)
ukb$sex     <- as.numeric(ukb$Gender_num)  # 1=Male, 0=Female
ukb$hdl     <- as.numeric(ukb$HDL.1)
ukb$tg      <- as.numeric(ukb$TRG.1)
ukb$tc      <- as.numeric(ukb$CHOL)
ukb$apob    <- as.numeric(ukb$APOB)
ukb$ldl_treated   <- as.numeric(ukb$LDL_treated)
ukb$ldl_untreated <- as.numeric(ukb$LDL_untreated)
ukb$trig_filter   <- as.numeric(ukb$Trig_Filter)
ukb$rf      <- as.numeric(ukb$reduction_factor); ukb$rf[is.na(ukb$rf)] <- 0
ukb$on_statin <- as.integer(ukb$rf > 0)
ukb$gene    <- as.character(ukb$gene)

# Non-HDL
ukb$non_hdl <- ifelse(!is.na(ukb$tc) & !is.na(ukb$hdl), ukb$tc - ukb$hdl, NA)

# LDL_RW = LDL_untreated (already reverse-engineered in the features file)
ukb$ldl_rw <- ukb$ldl_untreated

# Trig_Filter_RW = Trig_Filter (already computed on untreated LDL)
ukb$trig_filter_rw <- ukb$trig_filter

# ASCVD from self-reported conditions
ascvd_codes <- c(1075, 1074, 1081, 1583)
p20002_cols <- grep("^participant\\.p20002", names(ukb), value = TRUE)
ukb$has_ascvd <- apply(ukb[, ..p20002_cols], 1, function(row) {
  as.integer(any(row %in% ascvd_codes, na.rm = TRUE))
})

cat(sprintf("  FH+: %d (%.3f%%)\n", sum(ukb$fh), 100 * mean(ukb$fh)))
cat(sprintf("  On statin: %d\n", sum(ukb$on_statin)))
cat(sprintf("  Has ASCVD: %d\n", sum(ukb$has_ascvd)))

# ══════════════════════════════════════════════════════════════════════════════
# COMPUTE TUDOR SCORE (fixed weights)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== COMPUTING TUDOR SCORES ===\n")

ukb$tudor_score <- tudor_score(ukb$ldl_rw, ukb$trig_filter_rw, ukb$hdl, ukb$age, ukb$sex)
ukb$tudor_prob  <- tudor_prob(ukb$tudor_score)

valid <- !is.na(ukb$tudor_prob) & !is.na(ukb$fh)
cat(sprintf("  Valid scores: %d / %d\n", sum(valid), nrow(ukb)))
cat(sprintf("  Score range: %.3f to %.3f\n",
            min(ukb$tudor_score[valid]), max(ukb$tudor_score[valid])))

# ══════════════════════════════════════════════════════════════════════════════
# LIPID CLINIC FILTER
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== LIPID CLINIC FILTER ===\n")

ukb$lc_tc     <- !is.na(ukb$tc) & ukb$tc > 7.5
ukb$lc_ldl    <- !is.na(ukb$ldl_rw) & ukb$ldl_rw > 4.9
ukb$lc_nhdl   <- !is.na(ukb$non_hdl) & ukb$non_hdl > 5.9
ukb$lc_ascvd  <- ukb$has_ascvd == 1 & !is.na(ukb$age) & ukb$age < 60
ukb$lc_statin <- ukb$on_statin == 1

ukb$lipid_clinic <- ukb$lc_tc | ukb$lc_ldl | ukb$lc_nhdl | ukb$lc_ascvd | ukb$lc_statin

# Also create original "high-risk" = LDL_RW > 4.9 (from 02_external_validation.R)
ukb$high_risk <- !is.na(ukb$ldl_rw) & ukb$ldl_rw > 4.9

cat(sprintf("  TC > 7.5:              %d\n", sum(ukb$lc_tc)))
cat(sprintf("  LDL > 4.9:             %d\n", sum(ukb$lc_ldl)))
cat(sprintf("  Non-HDL > 5.9:         %d\n", sum(ukb$lc_nhdl)))
cat(sprintf("  Premature ASCVD (<60): %d\n", sum(ukb$lc_ascvd)))
cat(sprintf("  On statin:             %d\n", sum(ukb$lc_statin)))
cat(sprintf("  ─────────────────────────────\n"))
cat(sprintf("  Lipid clinic total:    %d (%.1f%%)\n", sum(ukb$lipid_clinic), 100*mean(ukb$lipid_clinic)))
cat(sprintf("  High-risk (LDL>4.9):   %d (%.1f%%)\n", sum(ukb$high_risk), 100*mean(ukb$high_risk)))

# ── eDLCN ──
ukb$edlcn <- ifelse(ukb$ldl_rw >= 8.5, 8,
             ifelse(ukb$ldl_rw >= 6.5, 5,
             ifelse(ukb$ldl_rw >= 5.0, 3,
             ifelse(ukb$ldl_rw >= 4.0, 1, 0))))
ukb$edlcn <- ukb$edlcn +
  ifelse(ukb$has_ascvd == 1 & ukb$age < 55 & ukb$sex == 1, 2,
  ifelse(ukb$has_ascvd == 1 & ukb$age < 60 & ukb$sex == 0, 2, 0))

# ══════════════════════════════════════════════════════════════════════════════
# RESULTS FOR EACH COHORT DEFINITION
# ══════════════════════════════════════════════════════════════════════════════

run_validation <- function(data, label) {
  d <- data[!is.na(data$tudor_prob) & !is.na(data$fh), ]
  if (nrow(d) < 100 || sum(d$fh) < 5) {
    cat(sprintf("  %s: insufficient data\n", label)); return(NULL)
  }

  cat(sprintf("\n═══ %s (n=%s, FH+=%d, prev=%.2f%%) ═══\n",
              label, format(nrow(d), big.mark=","), sum(d$fh), 100*mean(d$fh)))

  # TUDOR
  roc_t <- roc(d$fh, d$tudor_prob, quiet=TRUE)
  ci_t  <- ci.auc(roc_t, method="bootstrap", boot.n=2000, quiet=TRUE)

  # eDLCN
  roc_d <- roc(d$fh, d$edlcn, quiet=TRUE)
  ci_d  <- ci.auc(roc_d, quiet=TRUE)

  # LDL alone
  roc_l <- roc(d$fh, d$ldl_rw, quiet=TRUE)

  # Trig_Filter alone
  roc_tf <- roc(d$fh, d$trig_filter_rw, quiet=TRUE)

  # DeLong
  dt <- roc.test(roc_t, roc_d, method="delong")

  cat(sprintf("  TUDOR:       AUC=%.3f (%.3f-%.3f)\n", auc(roc_t), ci_t[1], ci_t[3]))
  cat(sprintf("  eDLCN:       AUC=%.3f (%.3f-%.3f)\n", auc(roc_d), ci_d[1], ci_d[3]))
  cat(sprintf("  LDL-C alone: AUC=%.3f\n", auc(roc_l)))
  cat(sprintf("  Trig_Filter: AUC=%.3f\n", auc(roc_tf)))
  cat(sprintf("  DeLong TUDOR vs eDLCN: Z=%.2f, p=%s\n", dt$statistic,
              formatC(dt$p.value, format="e", digits=2)))

  # Clinical performance at Youden
  co <- coords(roc_t, "best", best.method="youden")
  prev <- mean(d$fh)
  sens <- co$sensitivity; spec <- co$specificity
  ppv <- (sens*prev)/(sens*prev+(1-spec)*(1-prev))
  npv <- (spec*(1-prev))/((1-sens)*prev+spec*(1-prev))

  cat(sprintf("  Threshold=%.4f: Sens=%.1f%%, Spec=%.1f%%, PPV=%.1f%%, NPV=%.1f%%\n",
              co$threshold, sens*100, spec*100, ppv*100, npv*100))

  # Brier
  brier <- mean((d$tudor_prob - d$fh)^2)
  cat(sprintf("  Brier: %.4f\n", brier))

  data.frame(Cohort=label, N=nrow(d), N_FH=sum(d$fh), Prev=round(100*mean(d$fh),2),
             TUDOR_AUC=round(auc(roc_t),3), TUDOR_CI_lo=round(ci_t[1],3), TUDOR_CI_hi=round(ci_t[3],3),
             eDLCN_AUC=round(auc(roc_d),3), LDL_AUC=round(auc(roc_l),3), TF_AUC=round(auc(roc_tf),3),
             DeLong_p=dt$p.value, Sens=round(sens,3), Spec=round(spec,3),
             PPV=round(ppv,3), NPV=round(npv,3), Brier=round(brier,4))
}

results <- list()

# 1. Full UKB population
results[[1]] <- run_validation(as.data.frame(ukb), "Full UKB")

# 2. High-risk (LDL > 4.9) — original 02_external_validation definition
results[[2]] <- run_validation(as.data.frame(ukb[ukb$high_risk == TRUE, ]), "High-risk (LDL>4.9)")

# 3. Lipid clinic (broad NICE criteria)
results[[3]] <- run_validation(as.data.frame(ukb[ukb$lipid_clinic == TRUE, ]), "Lipid clinic (broad)")

# 4. Statin-free only
results[[4]] <- run_validation(as.data.frame(ukb[ukb$on_statin == 0 & ukb$high_risk == TRUE, ]),
                                "Statin-free + LDL>4.9")

# 5. On statin only
results[[5]] <- run_validation(as.data.frame(ukb[ukb$on_statin == 1, ]), "On statin only")

# 6. Gene-specific
for (g in c("LDLR", "APOB")) {
  fh_g <- ukb[ukb$fh == 1 & ukb$gene == g, ]
  neg  <- ukb[ukb$fh == 0 & ukb$high_risk == TRUE, ]
  combined <- rbind(as.data.frame(fh_g), as.data.frame(neg))
  results[[length(results)+1]] <- run_validation(combined, paste0("Gene: ", g, " (vs HR neg)"))
}

# Compile
res_df <- do.call(rbind, Filter(Negate(is.null), results))
cat("\n\n=== COMPLETE RESULTS TABLE ===\n")
print(res_df)

OUT_DIR <- file.path("C:/Users/nader/Downloads/calon_ukb_pipeline/tudor_loco_output")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
write.csv(res_df, file.path(OUT_DIR, "ukb_fixed_weights_validation.csv"), row.names=FALSE)
cat(sprintf("\nSaved: %s\n", file.path(OUT_DIR, "ukb_fixed_weights_validation.csv")))
cat(sprintf("Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

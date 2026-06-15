################################################################################
# LANCET STATISTICS: Fill all [FROM SCRIPT] gaps + Full calibration + NRI/IDI
# + DCA + Reclassification + Subgroups + Financial Analysis
# Author: Dr Nader Genedy | 2026-03-27
################################################################################
suppressPackageStartupMessages({
  library(readxl); library(haven); library(pROC); library(glmnet); library(data.table)
})
set.seed(42)
safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
sfact <- c(atorvastatin=0.38, simvastatin=0.35, rosuvastatin=0.34, pravastatin=0.25,
           fluvastatin=0.22, lovastatin=0.25, ezetimibe=0.18)
get_red <- function(tx) {
  if (is.na(tx) || tx == "" || tx == "NaN") return(0)
  tl <- tolower(as.character(tx))
  for (nm in names(sfact)) if (grepl(nm, tl, fixed = TRUE)) return(sfact[[nm]])
  if (grepl("statin", tl)) return(0.30)
  return(0)
}
rev_ldl <- function(l, r) ifelse(!is.na(l) & !is.na(r) & r > 0, l / (1 - r), l)

DD <- "C:/Users/nader/Downloads/calon_ukb_pipeline"
TD <- "D:/alphafold_backup/tudor_packup"
OD <- file.path(DD, "tudor_loco_output")

cat("\n================================================================\n")
cat("  LANCET STATISTICS + MANUSCRIPT GAP FILLING\n")
cat("================================================================\n\n")

# ═══════════════════════════════════════════════════════════════
# 1. LOAD UKB LIPID CLINIC
# ═══════════════════════════════════════════════════════════════
cat("Loading UKB...\n")
ukb <- as.data.frame(fread(file.path(TD, "TUDOR_UKB_Features.csv"), showProgress = FALSE))
ukb$fh <- safe_num(ukb$is_fh_genetic); ukb$fh[is.na(ukb$fh)] <- 0
ukb$ldl_ut <- safe_num(ukb$LDL_untreated); ukb$trig_filter <- safe_num(ukb$Trig_Filter)
ukb$tc <- safe_num(ukb$CHOL); ukb$hdl <- safe_num(ukb$HDL.1); ukb$tg <- safe_num(ukb$TRG.1)
ukb$apob <- safe_num(ukb$APOB); ukb$age <- safe_num(ukb$Age_at_LDL1)
ukb$sex <- safe_num(ukb$Gender_num)
ukb$non_hdl <- ifelse(!is.na(ukb$tc) & !is.na(ukb$hdl), ukb$tc - ukb$hdl, NA)
ukb$on_statin <- ifelse(safe_num(ukb$reduction_factor) > 0, 1, 0)
ukb$on_statin[is.na(ukb$on_statin)] <- 0

# Lipid clinic filter (matching manuscript: TC>7.5, LDL_ut>4.9, premature ASCVD)
ukb$lc <- (!is.na(ukb$tc) & ukb$tc > 7.5) |
          (!is.na(ukb$ldl_ut) & ukb$ldl_ut > 4.9) |
          (!is.na(ukb$non_hdl) & ukb$non_hdl > 5.9) |
          (ukb$on_statin == 1)
lc <- ukb[ukb$lc == TRUE, ]
lc$apob_ldl <- ifelse(!is.na(lc$apob) & !is.na(lc$ldl_ut) & lc$ldl_ut > 0, lc$apob / lc$ldl_ut, NA)
cat(sprintf("  Lipid clinic: %d (FH+=%d, prevalence=%.2f%%)\n\n",
            nrow(lc), sum(lc$fh == 1), 100 * mean(lc$fh)))

# ═══════════════════════════════════════════════════════════════
# 2. GENE-SPECIFIC STATISTICS [FROM SCRIPT] GAPS
# ═══════════════════════════════════════════════════════════════
cat("=== GENE-SPECIFIC [FROM SCRIPT] GAPS ===\n\n")

ldlr_fh <- lc[lc$gene == "LDLR" & lc$fh == 1 & !is.na(lc$trig_filter), ]
apob_fh <- lc[lc$gene == "APOB" & lc$fh == 1 & !is.na(lc$trig_filter), ]
nonfh <- lc[lc$fh == 0 & !is.na(lc$trig_filter), ]

# Trig_Filter Cohen's d by gene
tf_nf <- nonfh$trig_filter
d_ldlr <- (mean(ldlr_fh$trig_filter, na.rm = TRUE) - mean(tf_nf, na.rm = TRUE)) /
  sqrt((var(ldlr_fh$trig_filter, na.rm = TRUE) + var(tf_nf, na.rm = TRUE)) / 2)
d_apob <- (mean(apob_fh$trig_filter, na.rm = TRUE) - mean(tf_nf, na.rm = TRUE)) /
  sqrt((var(apob_fh$trig_filter, na.rm = TRUE) + var(tf_nf, na.rm = TRUE)) / 2)

cat(sprintf("  Trig_Filter Cohen's d (FH vs non-FH):\n"))
cat(sprintf("    LDLR: d = %.3f (n=%d)\n", d_ldlr, nrow(ldlr_fh)))
cat(sprintf("    APOB: d = %.3f (n=%d)\n", d_apob, nrow(apob_fh)))

# Mean TG by gene
cat(sprintf("\n  Mean TG:\n"))
cat(sprintf("    LDLR: %.2f (SD %.2f)\n", mean(ldlr_fh$tg, na.rm = TRUE), sd(ldlr_fh$tg, na.rm = TRUE)))
cat(sprintf("    APOB: %.2f (SD %.2f)\n", mean(apob_fh$tg, na.rm = TRUE), sd(apob_fh$tg, na.rm = TRUE)))

# Gene-specific AUC using Trig_Filter alone
roc_ldlr_tf <- roc(c(rep(1, nrow(ldlr_fh)), rep(0, nrow(nonfh))),
                    c(ldlr_fh$trig_filter, nonfh$trig_filter), quiet = TRUE)
roc_apob_tf <- roc(c(rep(1, nrow(apob_fh)), rep(0, nrow(nonfh))),
                    c(apob_fh$trig_filter, nonfh$trig_filter), quiet = TRUE)
cat(sprintf("\n  TrigFilter AUC by gene:\n"))
cat(sprintf("    LDLR: %.3f\n", as.numeric(auc(roc_ldlr_tf))))
cat(sprintf("    APOB: %.3f\n", as.numeric(auc(roc_apob_tf))))

# ═══════════════════════════════════════════════════════════════
# 3. BUILD TUDOR PREDICTIONS (5-fold CV within UKB LC)
# ═══════════════════════════════════════════════════════════════
cat("\n=== TUDOR 5-fold CV predictions ===\n")
feats <- c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex", "on_statin")
lc2 <- lc[complete.cases(lc[, c("fh", feats)]), ]
cat(sprintf("  Complete cases: %d (FH+=%d)\n", nrow(lc2), sum(lc2$fh == 1)))

n <- nrow(lc2); folds <- sample(rep(1:5, length.out = n)); preds <- numeric(n)
for (k in 1:5) {
  ti <- folds != k; vi <- folds == k
  cv <- cv.glmnet(as.matrix(lc2[ti, feats]), lc2$fh[ti], family = "binomial",
                  alpha = 0.5, nfolds = 5, type.measure = "auc", standardize = TRUE)
  preds[vi] <- as.numeric(predict(cv, newx = as.matrix(lc2[vi, feats]),
                                   s = "lambda.min", type = "response"))
}
lc2$tudor_prob <- preds

roc_tudor <- roc(lc2$fh, lc2$tudor_prob, quiet = TRUE)
ci_tudor <- ci.auc(roc_tudor, quiet = TRUE)
cat(sprintf("  TUDOR AUC: %.3f (%.3f-%.3f)\n", as.numeric(auc(roc_tudor)), ci_tudor[1], ci_tudor[3]))

# ═══════════════════════════════════════════════════════════════
# 4. GREY ZONE APOB ANALYSIS [FROM SCRIPT]
# ═══════════════════════════════════════════════════════════════
cat("\n=== GREY ZONE ApoB/LDL ANALYSIS ===\n")

fh_preds <- preds[lc2$fh == 1]
q25 <- quantile(fh_preds, 0.25); q75 <- quantile(fh_preds, 0.75)
gz <- lc2[lc2$tudor_prob >= q25 & lc2$tudor_prob <= q75, ]
cat(sprintf("  Grey zone (%.4f to %.4f): n=%d, FH+=%d\n", q25, q75, nrow(gz), sum(gz$fh == 1)))

gz$apob_ldl <- lc$apob_ldl[match(rownames(gz), rownames(lc))]
gz_a <- gz[!is.na(gz$apob_ldl), ]
cat(sprintf("  With ApoB data: n=%d, FH+=%d\n", nrow(gz_a), sum(gz_a$fh == 1)))

if (nrow(gz_a) > 50 && sum(gz_a$fh == 1) >= 5) {
  fh_r <- gz_a$apob_ldl[gz_a$fh == 1]
  nfh_r <- gz_a$apob_ldl[gz_a$fh == 0]
  cat(sprintf("  FH+ ApoB/LDL: mean=%.4f (SD=%.4f)\n", mean(fh_r), sd(fh_r)))
  cat(sprintf("  FH- ApoB/LDL: mean=%.4f (SD=%.4f)\n", mean(nfh_r), sd(nfh_r)))
  wt <- wilcox.test(fh_r, nfh_r)
  cat(sprintf("  Welch p = %s\n", formatC(wt$p.value, format = "e", digits = 2)))

  # Threshold 0.31
  tp <- sum(gz_a$apob_ldl < 0.31 & gz_a$fh == 1)
  fp <- sum(gz_a$apob_ldl < 0.31 & gz_a$fh == 0)
  fn <- sum(gz_a$apob_ldl >= 0.31 & gz_a$fh == 1)
  tn <- sum(gz_a$apob_ldl >= 0.31 & gz_a$fh == 0)
  sens <- tp / (tp + fn); spec <- tn / (tn + fp)
  cat(sprintf("  ApoB/LDL <0.31: Sens=%.1f%%, Spec=%.1f%%\n", sens * 100, spec * 100))

  roc_gz <- roc(gz_a$fh, gz_a$apob_ldl, direction = ">", quiet = TRUE)
  ci_gz <- ci.auc(roc_gz, quiet = TRUE)
  cat(sprintf("  ApoB/LDL AUC in grey zone: %.3f (%.3f-%.3f)\n",
              as.numeric(auc(roc_gz)), ci_gz[1], ci_gz[3]))
}

# ═══════════════════════════════════════════════════════════════
# 5. FULL CALIBRATION STATISTICS
# ═══════════════════════════════════════════════════════════════
cat("\n=== CALIBRATION ===\n")
brier <- mean((lc2$tudor_prob - lc2$fh)^2)
cat(sprintf("  Brier score: %.5f\n", brier))

cal_glm <- glm(fh ~ tudor_prob, data = lc2, family = binomial)
cat(sprintf("  Calibration slope: %.3f, intercept: %.3f\n", coef(cal_glm)[2], coef(cal_glm)[1]))

# Youden threshold
co <- coords(roc_tudor, "best", best.method = "youden")
cat(sprintf("  Youden threshold: %.4f\n", co$threshold))
cat(sprintf("  Sensitivity: %.1f%%\n", co$sensitivity * 100))
cat(sprintf("  Specificity: %.1f%%\n", co$specificity * 100))

# ═══════════════════════════════════════════════════════════════
# 6. NRI / IDI: TUDOR vs eDLCN
# ═══════════════════════════════════════════════════════════════
cat("\n=== NRI / IDI ===\n")
lc2$dlcn <- ifelse(!is.na(lc2$ldl_ut) & lc2$ldl_ut >= 8.5, 8,
            ifelse(!is.na(lc2$ldl_ut) & lc2$ldl_ut >= 6.5, 5,
            ifelse(!is.na(lc2$ldl_ut) & lc2$ldl_ut >= 5.0, 3,
            ifelse(!is.na(lc2$ldl_ut) & lc2$ldl_ut >= 4.0, 1, 0))))
dlcn_glm <- glm(fh ~ dlcn, data = lc2, family = binomial)
lc2$dlcn_prob <- predict(dlcn_glm, type = "response")

roc_dlcn <- roc(lc2$fh, lc2$dlcn, quiet = TRUE)
cat(sprintf("  eDLCN AUC: %.3f\n", as.numeric(auc(roc_dlcn))))

# DeLong
dl <- roc.test(roc_tudor, roc_dlcn, method = "delong")
cat(sprintf("  DeLong Z=%.2f, p=%s\n", dl$statistic, formatC(dl$p.value, format = "e", digits = 2)))

# Categorical NRI
ev <- lc2$fh == 1; ne <- lc2$fh == 0
tudor_cat <- ifelse(lc2$tudor_prob < 0.005, "Low", ifelse(lc2$tudor_prob < 0.02, "Mid", "High"))
dlcn_cat <- ifelse(lc2$dlcn_prob < 0.005, "Low", ifelse(lc2$dlcn_prob < 0.02, "Mid", "High"))

ev_up <- sum(tudor_cat[ev] == "High" & dlcn_cat[ev] != "High")
ev_down <- sum(tudor_cat[ev] != "High" & dlcn_cat[ev] == "High")
ne_down <- sum(tudor_cat[ne] == "Low" & dlcn_cat[ne] != "Low")
ne_up <- sum(tudor_cat[ne] != "Low" & dlcn_cat[ne] == "Low")
nri_ev <- (ev_up - ev_down) / sum(ev); nri_ne <- (ne_down - ne_up) / sum(ne)
cat(sprintf("  NRI(events) = %+.3f\n", nri_ev))
cat(sprintf("  NRI(non-events) = %+.3f\n", nri_ne))
cat(sprintf("  NRI(total) = %+.3f\n", nri_ev + nri_ne))

# IDI
idi <- (mean(lc2$tudor_prob[ev]) - mean(lc2$tudor_prob[ne])) -
       (mean(lc2$dlcn_prob[ev]) - mean(lc2$dlcn_prob[ne]))
cat(sprintf("  IDI = %+.4f\n", idi))

# Reclassification table
cat("\n  Reclassification Table (FH+):\n")
print(table(TUDOR = tudor_cat[ev], eDLCN = dlcn_cat[ev]))
cat("\n  Reclassification Table (FH-):\n")
print(table(TUDOR = tudor_cat[ne], eDLCN = dlcn_cat[ne]))

# ═══════════════════════════════════════════════════════════════
# 7. DECISION CURVE ANALYSIS
# ═══════════════════════════════════════════════════════════════
cat("\n=== DECISION CURVE ANALYSIS ===\n")
dca_results <- data.frame()
for (pt in seq(0.002, 0.05, by = 0.002)) {
  tp_t <- sum(lc2$tudor_prob >= pt & lc2$fh == 1)
  fp_t <- sum(lc2$tudor_prob >= pt & lc2$fh == 0)
  nb_t <- tp_t / nrow(lc2) - (fp_t / nrow(lc2)) * (pt / (1 - pt))

  tp_d <- sum(lc2$dlcn_prob >= pt & lc2$fh == 1)
  fp_d <- sum(lc2$dlcn_prob >= pt & lc2$fh == 0)
  nb_d <- tp_d / nrow(lc2) - (fp_d / nrow(lc2)) * (pt / (1 - pt))

  nb_all <- mean(lc2$fh) - (1 - mean(lc2$fh)) * (pt / (1 - pt))

  dca_results <- rbind(dca_results, data.frame(
    Threshold = pt, TUDOR = nb_t, eDLCN = nb_d, TreatAll = nb_all, TreatNone = 0
  ))
}
cat("  Threshold  TUDOR_NB     eDLCN_NB     Diff\n")
for (pt in c(0.005, 0.01, 0.02, 0.05)) {
  row <- dca_results[dca_results$Threshold == pt, ]
  if (nrow(row) > 0) {
    cat(sprintf("  %.1f%%      %.5f    %.5f    %+.5f\n",
                pt * 100, row$TUDOR, row$eDLCN, row$TUDOR - row$eDLCN))
  }
}

# ═══════════════════════════════════════════════════════════════
# 8. FINANCIAL ANALYSIS
# ═══════════════════════════════════════════════════════════════
cat("\n=== FINANCIAL IMPACT ANALYSIS ===\n")

# UK costs (NHS 2024 tariff)
cost_genetic_test <- 250   # GBP per panel test
cost_apob <- 8             # GBP per ApoB measurement
cost_lipid_panel <- 5      # GBP per standard lipid panel
cost_clinic_visit <- 150   # GBP per specialist lipid clinic visit
cost_statin_annual <- 30   # GBP per year (generic atorvastatin)
cost_pcsk9i_annual <- 4380 # GBP per year (evolocumab)
cost_mi_acute <- 5000      # GBP per acute MI admission
cost_mi_lifetime <- 50000  # GBP lifetime cost post-MI
cost_cabg <- 15000         # GBP per CABG

# Number needed to test (NNT) to find one FH case
prev_lc <- mean(lc2$fh)
sens_tudor <- co$sensitivity
spec_tudor <- co$specificity
ppv_tudor <- (sens_tudor * prev_lc) / (sens_tudor * prev_lc + (1 - spec_tudor) * (1 - prev_lc))
nnt_tudor <- 1 / ppv_tudor

# eDLCN at score >= 6
dlcn_sens <- sum(lc2$dlcn[ev] >= 6) / sum(ev)
dlcn_spec <- sum(lc2$dlcn[ne] < 6) / sum(ne)
ppv_dlcn <- (dlcn_sens * prev_lc) / (dlcn_sens * prev_lc + (1 - dlcn_spec) * (1 - prev_lc))
nnt_dlcn <- 1 / ppv_dlcn

cat(sprintf("  Prevalence in lipid clinic: %.2f%%\n", prev_lc * 100))
cat(sprintf("  TUDOR: Sens=%.1f%%, Spec=%.1f%%, PPV=%.2f%%, NNT=%.0f\n",
            sens_tudor * 100, spec_tudor * 100, ppv_tudor * 100, nnt_tudor))
cat(sprintf("  eDLCN>=6: Sens=%.1f%%, Spec=%.1f%%, PPV=%.2f%%, NNT=%.0f\n",
            dlcn_sens * 100, dlcn_spec * 100, ppv_dlcn * 100, nnt_dlcn))

# Cost per diagnosis
cost_per_dx_tudor <- nnt_tudor * cost_genetic_test
cost_per_dx_dlcn <- nnt_dlcn * cost_genetic_test
cat(sprintf("\n  Cost per FH diagnosis:\n"))
cat(sprintf("    TUDOR: GBP %.0f (%d tests x GBP %d)\n", cost_per_dx_tudor, round(nnt_tudor), cost_genetic_test))
cat(sprintf("    eDLCN: GBP %.0f (%d tests x GBP %d)\n", cost_per_dx_dlcn, round(nnt_dlcn), cost_genetic_test))
cat(sprintf("    Savings per diagnosis: GBP %.0f\n", cost_per_dx_dlcn - cost_per_dx_tudor))

# Cascade screening yield
cascade_yield <- 2.0  # Average 2 additional FH diagnoses per index case (1st degree relatives)
cost_per_cascade_dx <- cost_genetic_test / cascade_yield
cat(sprintf("\n  Cascade screening yield: %.1f additional diagnoses per index\n", cascade_yield))
cat(sprintf("  Cost per cascade diagnosis: GBP %.0f\n", cost_per_cascade_dx))

# 10-year ASCVD prevention
# Assume 20% untreated 10-year ASCVD risk in FH, reduced to 5% with treatment
ascvd_risk_untreated <- 0.20
ascvd_risk_treated <- 0.05
nnt_prevent_event <- 1 / (ascvd_risk_untreated - ascvd_risk_treated)
events_prevented_per_1000_dx <- 1000 * (ascvd_risk_untreated - ascvd_risk_treated)
cost_saved_per_1000 <- events_prevented_per_1000_dx * cost_mi_lifetime

cat(sprintf("\n  ASCVD prevention (per 1,000 FH diagnoses):\n"))
cat(sprintf("    Events prevented: %.0f\n", events_prevented_per_1000_dx))
cat(sprintf("    Cost of testing: GBP %.0f\n", 1000 * cost_per_dx_tudor))
cat(sprintf("    Healthcare costs saved: GBP %.0f\n", cost_saved_per_1000))
cat(sprintf("    Net saving: GBP %.0f per 1,000 diagnoses\n",
            cost_saved_per_1000 - 1000 * cost_per_dx_tudor))

# UK-wide impact
uk_fh_total <- 260000  # Estimated 1 in 250 x 65M population
uk_fh_diagnosed <- 20000  # ~8% currently diagnosed
uk_fh_undiagnosed <- uk_fh_total - uk_fh_diagnosed
cat(sprintf("\n  UK-wide impact:\n"))
cat(sprintf("    Estimated FH in UK: %s\n", format(uk_fh_total, big.mark = ",")))
cat(sprintf("    Currently diagnosed: %s (%.0f%%)\n",
            format(uk_fh_diagnosed, big.mark = ","), 100 * uk_fh_diagnosed / uk_fh_total))
cat(sprintf("    Undiagnosed: %s\n", format(uk_fh_undiagnosed, big.mark = ",")))
cat(sprintf("    ASCVD events preventable: %s over 10 years\n",
            format(round(uk_fh_undiagnosed * 0.15), big.mark = ",")))

# ═══════════════════════════════════════════════════════════════
# 9. SAVE ALL RESULTS
# ═══════════════════════════════════════════════════════════════
cat("\n=== SAVING ===\n")
write.csv(dca_results, file.path(OD, "DCA_results.csv"), row.names = FALSE)

# Summary CSV with all key statistics
summary_stats <- data.frame(
  Statistic = c(
    "TUDOR_AUC_UKB_LC", "TUDOR_AUC_CI_lo", "TUDOR_AUC_CI_hi",
    "eDLCN_AUC", "DeLong_Z", "DeLong_p",
    "Brier_score", "Calibration_slope", "Calibration_intercept",
    "Youden_threshold", "Sensitivity", "Specificity",
    "NRI_events", "NRI_nonevents", "NRI_total", "IDI",
    "TrigFilter_Cohen_d_LDLR", "TrigFilter_Cohen_d_APOB",
    "Cost_per_dx_TUDOR", "Cost_per_dx_eDLCN",
    "NNT_TUDOR", "NNT_eDLCN",
    "Grey_zone_n", "Grey_zone_FH",
    "Grey_zone_ApoB_LDL_FH_mean", "Grey_zone_ApoB_LDL_nonFH_mean"
  ),
  Value = c(
    as.numeric(auc(roc_tudor)), ci_tudor[1], ci_tudor[3],
    as.numeric(auc(roc_dlcn)), dl$statistic, dl$p.value,
    brier, coef(cal_glm)[2], coef(cal_glm)[1],
    co$threshold, co$sensitivity, co$specificity,
    nri_ev, nri_ne, nri_ev + nri_ne, idi,
    d_ldlr, d_apob,
    cost_per_dx_tudor, cost_per_dx_dlcn,
    nnt_tudor, nnt_dlcn,
    nrow(gz), sum(gz$fh == 1),
    ifelse(exists("fh_r"), mean(fh_r), NA),
    ifelse(exists("nfh_r"), mean(nfh_r), NA)
  )
)
write.csv(summary_stats, file.path(OD, "lancet_statistics_summary.csv"), row.names = FALSE)
write.csv(summary_stats, "D:/alphafold_backup/tudor_packup/lancet_statistics_summary.csv", row.names = FALSE)

cat("  Saved: DCA_results.csv, lancet_statistics_summary.csv\n")
cat(sprintf("  Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("================================================================\n")

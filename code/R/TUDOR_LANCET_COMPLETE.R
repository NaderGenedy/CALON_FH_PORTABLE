################################################################################
#
#  TUDOR FH DIAGNOSTIC: COMPLETE LANCET-GRADE STATISTICAL ANALYSIS
#
#  Contents:
#    PART 1: Data loading & cohort assembly (SW, Wales, UKB lipid clinic)
#    PART 2: LOCO-CV with all models (ENET, XGB)
#    PART 3: Calibration (slope, intercept, Hosmer-Lemeshow, plots)
#    PART 4: Decision Curve Analysis
#    PART 5: NRI / IDI (TUDOR vs eDLCN, categorical + continuous)
#    PART 6: Subgroup analyses (sex, age, statin, gene, index/cascade)
#    PART 7: ApoB/LDL discordance analysis (FH+ vs FH-)
#    PART 8: NMR lipoprotein profiling (UKB nested CV)
#    PART 9: Lp(a), ApoA, nonHDL-LDL gap analysis
#    PART 10: Publication figures (Nature specs)
#    PART 11: Summary tables (3 tables, Nature format)
#
#  Outcome: fh (genetically confirmed FH, diagnostic model)
#  Author: Dr Nader Genedy | 2026-03-28
#
################################################################################

cat("\n================================================================\n")
cat("  TUDOR: COMPLETE LANCET STATISTICS\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

set.seed(42)

suppressPackageStartupMessages({
  library(readxl); library(haven); library(pROC); library(glmnet)
  library(xgboost); library(data.table); library(boot)
})

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
sfact <- c(atorvastatin=0.38, simvastatin=0.35, rosuvastatin=0.34,
           pravastatin=0.25, fluvastatin=0.22, lovastatin=0.25, ezetimibe=0.18)
get_red <- function(tx) {
  if (is.na(tx) || tx == "" || tx == "NaN") return(0)
  tl <- tolower(as.character(tx))
  for (nm in names(sfact)) if (grepl(nm, tl, fixed=TRUE)) return(sfact[[nm]])
  if (grepl("statin", tl)) return(0.30); return(0)
}
rev_ldl <- function(l, r) ifelse(!is.na(l) & !is.na(r) & r > 0, l/(1-r), l)
parse_dt <- function(x) {
  x <- as.character(x); d <- suppressWarnings(as.Date(x, format="%d-%m-%Y"))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, format="%d/%m/%Y"))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, format="%Y-%m-%d")); d
}

DD  <- "C:/Users/nader/Downloads/calon_ukb_pipeline"
TD  <- "D:/alphafold_backup/tudor_packup"
OUT <- file.path(DD, "tudor_loco_output")
FIG <- file.path(OUT, "figures")
dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

# ══════════════════════════════════════════════════════════════════════════════
#  PART 1: DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════

cat("=== PART 1: Loading cohorts ===\n")

# SOUTH WALES
pos <- as.data.frame(read_excel(file.path(DD, "mutation positive group.xlsx")))
pos$fh <- 1
pos$dob_d  <- as.Date(sapply(pos$BirthDate, parse_dt), origin="1970-01-01")
pos$meas_d <- as.Date(sapply(pos$MeasurementDate.1, parse_dt), origin="1970-01-01")
pos$age <- as.numeric(difftime(pos$meas_d, pos$dob_d, units="days"))/365.25
pos$sex <- ifelse(pos$Gender %in% c("M","Male"), 1, 0)
pos$hdl <- safe_num(pos$HDL.1); pos$tg <- safe_num(pos$TRG.1)
pos$ldl_m <- safe_num(pos$LDL.1); pos$tc <- safe_num(pos$TC.1)
pos$rf <- sapply(pos$Treatment1.1, get_red)
pos$ldl_ut <- ifelse(pos$rf > 0, rev_ldl(pos$ldl_m, pos$rf), pos$ldl_m)
pos$trig_filter <- ifelse(!is.na(pos$ldl_ut) & !is.na(pos$tg) & pos$tg >= 0,
                           pmin(pos$ldl_ut/(pos$tg+0.1), 50), NA)
pos$on_statin <- ifelse(pos$rf > 0, 1, 0)
pos$tendon_xanth <- ifelse(pos$TendonXanthomata %in% c(1,"1","Yes","TRUE"), 1, 0)
pos$tendon_xanth[is.na(pos$tendon_xanth)] <- 0
ca1 <- if ("Corneal Arcus" %in% names(pos)) "Corneal Arcus" else "CornealArcus"
pos$corneal_arcus <- ifelse(pos[[ca1]] %in% c(1,"1","Yes","TRUE"), 1, 0)
pos$corneal_arcus[is.na(pos$corneal_arcus)] <- 0
pos$dlcn <- safe_num(pos$GenoTypingScore)
apob_p <- intersect(c("Apo-B","ApoB","Apo.B","apoB"), names(pos))[1]
pos$apob <- if (!is.na(apob_p)) safe_num(pos[[apob_p]]) else NA
pos$non_hdl <- ifelse(!is.na(pos$tc) & !is.na(pos$hdl), pos$tc - pos$hdl, NA)
pos$apob_ldl <- ifelse(!is.na(pos$apob) & !is.na(pos$ldl_ut) & pos$ldl_ut > 0,
                        pos$apob / pos$ldl_ut, NA)
pos$nhdl_ldl_gap <- ifelse(!is.na(pos$non_hdl) & !is.na(pos$ldl_ut),
                            pos$non_hdl - pos$ldl_ut, NA)
mut_col <- intersect(c("Mutation (1)","Mutation1","Mutation.1"), names(pos))[1]
pos_mut <- if (!is.na(mut_col)) pos[[mut_col]] else NA
pos$gene <- ifelse(grepl("^LDLR", pos_mut), "LDLR",
             ifelse(grepl("^APOB", pos_mut), "APOB",
             ifelse(grepl("^PCSK9", pos_mut), "PCSK9", "Other")))
ix_col <- intersect(c("I_vs_R","I_Vs_R","IndexVsRelative"), names(pos))[1]
pos$index_case <- if (!is.na(ix_col)) ifelse(safe_num(pos[[ix_col]]) == 1, 1, 0) else 1
pos$cohort <- "SouthWales"

neg <- as.data.frame(read_excel(file.path(DD, "Mutation negative control updated.xlsx")))
neg$fh <- 0; neg$age <- safe_num(neg[["age at the result"]])
neg$sex <- ifelse(neg$Gender %in% c("M","Male"), 1, 0)
neg$hdl <- safe_num(neg$HDL.1); neg$tg <- safe_num(neg$TRG.1)
neg$ldl_m <- safe_num(neg$LDL.1); neg$tc <- safe_num(neg$TC.1)
neg$rf <- sapply(neg$Treatment1.1, get_red)
neg$ldl_ut <- ifelse(neg$rf > 0, rev_ldl(neg$ldl_m, neg$rf), neg$ldl_m)
neg$trig_filter <- ifelse(!is.na(neg$ldl_ut) & !is.na(neg$tg) & neg$tg >= 0,
                           pmin(neg$ldl_ut/(neg$tg+0.1), 50), NA)
neg$on_statin <- ifelse(neg$rf > 0, 1, 0)
ca2 <- if ("Corneal Arcus" %in% names(neg)) "Corneal Arcus" else "CornealArcus"
neg$tendon_xanth <- ifelse(neg$TendonXanthomata %in% c(1,"1","Yes","TRUE"), 1, 0)
neg$tendon_xanth[is.na(neg$tendon_xanth)] <- 0
neg$corneal_arcus <- ifelse(neg[[ca2]] %in% c(1,"1","Yes","TRUE"), 1, 0)
neg$corneal_arcus[is.na(neg$corneal_arcus)] <- 0
neg$dlcn <- safe_num(neg$GenoTypingScore)
apob_n <- intersect(c("apoB...2","apoB","ApoB"), names(neg))[1]
neg$apob <- if (!is.na(apob_n)) safe_num(neg[[apob_n]]) else NA
neg$non_hdl <- ifelse(!is.na(neg$tc) & !is.na(neg$hdl), neg$tc - neg$hdl, NA)
neg$apob_ldl <- ifelse(!is.na(neg$apob) & !is.na(neg$ldl_ut) & neg$ldl_ut > 0,
                        neg$apob / neg$ldl_ut, NA)
neg$nhdl_ldl_gap <- ifelse(!is.na(neg$non_hdl) & !is.na(neg$ldl_ut),
                            neg$non_hdl - neg$ldl_ut, NA)
neg$gene <- NA; neg$index_case <- 0; neg$cohort <- "SouthWales"

cols <- c("fh","age","sex","hdl","tg","ldl_m","ldl_ut","trig_filter","tc",
          "on_statin","tendon_xanth","corneal_arcus","dlcn","apob",
          "non_hdl","apob_ldl","nhdl_ldl_gap","gene","index_case","cohort")
sw <- rbind(pos[, cols], neg[, cols])

# WALES PASS
pass <- as.data.frame(read_sav(file.path(DD, "PASS_wrong_dob.sav")))
pass$fh  <- ifelse(pass$Positive1 == "Yes", 1, ifelse(pass$Positive1 == "No", 0, NA))
pass$age <- safe_num(pass$BMI_AGE); pass$sex <- ifelse(pass$Gender == "M", 1, 0)
pass$hdl <- safe_num(pass$HDL.1); pass$tg <- safe_num(pass$TRG.1)
pass$ldl_m <- safe_num(pass$LDL.1); pass$tc <- safe_num(pass$TC.1)
pass$rf <- sapply(pass$Treatment1.1, get_red)
pass$ldl_ut <- ifelse(pass$rf > 0, rev_ldl(pass$ldl_m, pass$rf), pass$ldl_m)
pass$trig_filter <- ifelse(!is.na(pass$ldl_ut) & !is.na(pass$tg) & pass$tg >= 0,
                            pmin(pass$ldl_ut/(pass$tg+0.1), 50), NA)
pass$on_statin <- ifelse(pass$rf > 0, 1, 0)
pass$tendon_xanth <- ifelse(safe_num(pass$TendonXanthomata) == 1, 1, 0)
pass$tendon_xanth[is.na(pass$tendon_xanth)] <- 0
pass$corneal_arcus <- ifelse(safe_num(pass$CornealArcus) == 1, 1, 0)
pass$corneal_arcus[is.na(pass$corneal_arcus)] <- 0
pass$dlcn <- safe_num(pass$GenoTypingScore)
pass$apob <- NA; pass$non_hdl <- ifelse(!is.na(pass$tc) & !is.na(pass$hdl), pass$tc - pass$hdl, NA)
pass$apob_ldl <- NA; pass$nhdl_ldl_gap <- ifelse(!is.na(pass$non_hdl) & !is.na(pass$ldl_ut), pass$non_hdl - pass$ldl_ut, NA)
pass$gene <- ifelse(grepl("^LDLR", pass$Mutation1), "LDLR",
              ifelse(grepl("^APOB", pass$Mutation1), "APOB",
              ifelse(grepl("^PCSK9", pass$Mutation1), "PCSK9", NA)))
pass$index_case <- ifelse(safe_num(pass$I_Vs_R) == 1, 1, 0)
pass$index_case[is.na(pass$index_case)] <- 0; pass$cohort <- "Wales"
wales <- pass[!is.na(pass$fh), cols]

# UKB LIPID CLINIC
ukb_raw <- as.data.frame(fread(file.path(TD, "TUDOR_UKB_Features.csv"), showProgress=FALSE))
ukb_raw$fh <- safe_num(ukb_raw$is_fh_genetic); ukb_raw$fh[is.na(ukb_raw$fh)] <- 0
ukb_raw$age <- safe_num(ukb_raw$Age_at_LDL1); ukb_raw$sex <- safe_num(ukb_raw$Gender_num)
ukb_raw$hdl <- safe_num(ukb_raw$HDL.1); ukb_raw$tg <- safe_num(ukb_raw$TRG.1)
ukb_raw$ldl_m <- safe_num(ukb_raw$LDL_treated); ukb_raw$ldl_ut <- safe_num(ukb_raw$LDL_untreated)
ukb_raw$trig_filter <- safe_num(ukb_raw$Trig_Filter)
ukb_raw$on_statin <- ifelse(safe_num(ukb_raw$reduction_factor) > 0, 1, 0)
ukb_raw$on_statin[is.na(ukb_raw$on_statin)] <- 0
ukb_raw$tc <- safe_num(ukb_raw$CHOL); ukb_raw$apob <- safe_num(ukb_raw$APOB)
ukb_raw$non_hdl <- ifelse(!is.na(ukb_raw$tc) & !is.na(ukb_raw$hdl), ukb_raw$tc - ukb_raw$hdl, NA)
ukb_raw$apob_ldl <- ifelse(!is.na(ukb_raw$apob) & !is.na(ukb_raw$ldl_ut) & ukb_raw$ldl_ut > 0,
                            ukb_raw$apob / ukb_raw$ldl_ut, NA)
ukb_raw$nhdl_ldl_gap <- ifelse(!is.na(ukb_raw$non_hdl) & !is.na(ukb_raw$ldl_ut),
                                ukb_raw$non_hdl - ukb_raw$ldl_ut, NA)
# Lipid clinic filter
ukb_raw$lc <- (!is.na(ukb_raw$tc) & ukb_raw$tc > 7.5) |
              (!is.na(ukb_raw$ldl_ut) & ukb_raw$ldl_ut > 4.9) |
              (!is.na(ukb_raw$non_hdl) & ukb_raw$non_hdl > 5.9) |
              (ukb_raw$on_statin == 1)
ukb <- ukb_raw[ukb_raw$lc == TRUE, ]
ukb$dlcn <- ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut >= 8.5, 8,
            ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut >= 6.5, 5,
            ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut >= 5.0, 3,
            ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut >= 4.0, 1, 0))))
ukb$tendon_xanth <- 0; ukb$corneal_arcus <- 0
ukb$gene <- as.character(ukb$gene); ukb$index_case <- 0; ukb$cohort <- "UKB"
ukb_df <- ukb[, cols]

all_data <- rbind(sw, wales, ukb_df)
cat(sprintf("  SW: %d (FH+=%d) | Wales: %d (FH+=%d) | UKB LC: %d (FH+=%d)\n",
            nrow(sw), sum(sw$fh), nrow(wales), sum(wales$fh),
            nrow(ukb_df), sum(ukb_df$fh)))
cat(sprintf("  TOTAL: %d\n\n", nrow(all_data)))

features <- c("ldl_ut","trig_filter","hdl","tg","age","sex")
cohorts  <- c("SouthWales","Wales","UKB")

# ══════════════════════════════════════════════════════════════════════════════
#  PART 2: LOCO-CV — TRAIN, PREDICT, STORE
# ══════════════════════════════════════════════════════════════════════════════

cat("=== PART 2: LOCO-CV ===\n")

loco_preds <- list()  # Store predictions for downstream analyses

for (ho in cohorts) {
  train <- all_data[all_data$cohort != ho, ]
  test  <- all_data[all_data$cohort == ho, ]

  train_cc <- train[complete.cases(train[, c("fh", features)]), ]
  test_cc  <- test[complete.cases(test[, c("fh", features)]), ]

  X_tr <- as.matrix(train_cc[, features]); y_tr <- train_cc$fh
  X_te <- as.matrix(test_cc[, features]);  y_te <- test_cc$fh

  cv <- cv.glmnet(X_tr, y_tr, family="binomial", alpha=0.5,
                  nfolds=10, type.measure="auc", standardize=TRUE)
  pred <- as.numeric(predict(cv, newx=X_te, s="lambda.min", type="response"))
  lp   <- as.numeric(predict(cv, newx=X_te, s="lambda.min", type="link"))

  r  <- roc(y_te, pred, quiet=TRUE)
  ci <- ci.auc(r, method="bootstrap", boot.n=2000, quiet=TRUE)
  co <- coords(r, "best", best.method="youden")

  cat(sprintf("  %-12s AUC=%.4f [%.4f-%.4f] Sens=%.1f%% Spec=%.1f%% n=%d FH+=%d\n",
              ho, as.numeric(auc(r)), ci[1], ci[3],
              co$sensitivity*100, co$specificity*100, nrow(test_cc), sum(y_te)))

  # Store for downstream
  loco_preds[[ho]] <- data.frame(
    cohort = ho, fh = y_te, pred = pred, lp = lp,
    dlcn = test_cc$dlcn, age = test_cc$age, sex = test_cc$sex,
    on_statin = test_cc$on_statin, gene = test_cc$gene,
    index_case = test_cc$index_case, ldl_ut = test_cc$ldl_ut,
    apob = test_cc$apob, non_hdl = test_cc$non_hdl,
    apob_ldl = test_cc$apob_ldl, nhdl_ldl_gap = test_cc$nhdl_ldl_gap,
    trig_filter = test_cc$trig_filter, hdl = test_cc$hdl, tg = test_cc$tg,
    stringsAsFactors = FALSE
  )
}

# ══════════════════════════════════════════════════════════════════════════════
#  PART 3: CALIBRATION
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== PART 3: Calibration ===\n")

for (ho in cohorts) {
  d <- loco_preds[[ho]]
  cat(sprintf("\n  --- %s ---\n", ho))

  # Calibration slope + intercept
  cal_fit <- glm(fh ~ lp, data = d, family = binomial)
  cal_int   <- coef(cal_fit)[1]
  cal_slope <- coef(cal_fit)[2]
  cat(sprintf("  Calibration-in-the-large: intercept=%.4f\n", cal_int))
  cat(sprintf("  Calibration slope: %.4f (ideal=1.0)\n", cal_slope))

  # Hosmer-Lemeshow (manual — 10 groups)
  d$decile <- cut(d$pred, breaks = quantile(d$pred, probs = seq(0, 1, 0.1), na.rm = TRUE),
                  include.lowest = TRUE, labels = FALSE)
  hl_tab <- aggregate(cbind(obs = fh, n = rep(1, nrow(d))) ~ decile, data = d,
                       FUN = function(x) c(sum = sum(x), n = length(x)))
  obs <- tapply(d$fh, d$decile, mean)
  pred_m <- tapply(d$pred, d$decile, mean)
  n_g <- tapply(d$fh, d$decile, length)

  # Chi-squared HL test
  hl_chi2 <- sum(n_g * (obs - pred_m)^2 / (pred_m * (1 - pred_m) + 1e-10))
  hl_df <- length(n_g) - 2
  hl_p  <- 1 - pchisq(hl_chi2, df = hl_df)
  cat(sprintf("  Hosmer-Lemeshow: chi2=%.2f, df=%d, p=%.4f (%s)\n",
              hl_chi2, hl_df, hl_p,
              ifelse(hl_p > 0.05, "GOOD calibration", "Poor calibration")))

  # Brier score
  brier <- mean((d$pred - d$fh)^2)
  brier_max <- mean(d$fh) * (1 - mean(d$fh))
  brier_scaled <- 1 - brier / brier_max
  cat(sprintf("  Brier score: %.4f (scaled=%.4f, max=%.4f)\n", brier, brier_scaled, brier_max))

  # Calibration plot
  pdf(file.path(FIG, sprintf("calibration_%s.pdf", ho)), width=4, height=4)
  par(mar=c(4,4,2,1), family="sans")
  plot(pred_m, obs, xlim=c(0,1), ylim=c(0,1), pch=19, cex=1.2,
       col="#2166AC", xlab="Predicted probability", ylab="Observed proportion",
       main=sprintf("Calibration — %s", ho), cex.main=0.9, cex.lab=0.8, cex.axis=0.7)
  abline(0, 1, lty=2, col="grey50")
  lines(lowess(pred_m, obs, f=0.8), col="#B2182B", lwd=2)
  legend("bottomright", legend=c(sprintf("Slope=%.2f", cal_slope),
                                  sprintf("HL p=%.3f", hl_p),
                                  sprintf("Brier=%.4f", brier)),
         cex=0.6, bty="n")
  dev.off()
}

# ══════════════════════════════════════════════════════════════════════════════
#  PART 4: DECISION CURVE ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== PART 4: Decision Curve Analysis ===\n")

dca_compute <- function(y, pred, thresholds) {
  n <- length(y)
  prev <- mean(y)
  nb <- sapply(thresholds, function(pt) {
    tp <- sum(pred >= pt & y == 1)
    fp <- sum(pred >= pt & y == 0)
    tp/n - fp/n * pt/(1 - pt)
  })
  data.frame(threshold = thresholds, net_benefit = nb,
             treat_all = prev - (1 - prev) * thresholds / (1 - thresholds),
             treat_none = 0)
}

pdf(file.path(FIG, "decision_curve_all.pdf"), width=8, height=4)
par(mfrow=c(1,3), mar=c(4,4,2,1), family="sans")

for (ho in cohorts) {
  d <- loco_preds[[ho]]
  thresholds <- seq(0.01, 0.60, 0.01)

  # TUDOR
  dca_t <- dca_compute(d$fh, d$pred, thresholds)

  # eDLCN (normalise to 0-1)
  dlcn_ok <- !is.na(d$dlcn)
  if (sum(dlcn_ok) > 50) {
    dlcn_norm <- d$dlcn[dlcn_ok] / max(d$dlcn[dlcn_ok], na.rm=TRUE)
    dca_d <- dca_compute(d$fh[dlcn_ok], dlcn_norm, thresholds)
  }

  plot(dca_t$threshold, dca_t$net_benefit, type="l", lwd=2, col="#2166AC",
       xlab="Threshold probability", ylab="Net benefit",
       main=ho, ylim=c(-0.05, max(dca_t$treat_all, na.rm=TRUE) * 1.1),
       cex.main=0.9, cex.lab=0.8, cex.axis=0.7)
  lines(dca_t$threshold, dca_t$treat_all, lty=2, col="grey40")
  abline(h=0, lty=3, col="grey70")
  if (sum(dlcn_ok) > 50)
    lines(dca_d$threshold, dca_d$net_benefit, lwd=2, col="#B2182B")
  legend("topright", legend=c("TUDOR","eDLCN","Test all","Test none"),
         col=c("#2166AC","#B2182B","grey40","grey70"),
         lty=c(1,1,2,3), lwd=c(2,2,1,1), cex=0.5, bty="n")
}
dev.off()
cat("  DCA plots saved\n")

# ══════════════════════════════════════════════════════════════════════════════
#  PART 5: NRI / IDI
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== PART 5: NRI / IDI ===\n")

nri_results <- list()

for (ho in cohorts) {
  d <- loco_preds[[ho]]
  dlcn_ok <- !is.na(d$dlcn)
  if (sum(dlcn_ok) < 100) next

  y <- d$fh[dlcn_ok]
  p_tudor <- d$pred[dlcn_ok]
  p_dlcn  <- d$dlcn[dlcn_ok]
  p_dlcn_norm <- p_dlcn / max(p_dlcn, na.rm=TRUE)

  # Categorical NRI (thresholds: 0.25, 0.75)
  t_cat <- cut(p_tudor, breaks=c(0, 0.25, 0.75, 1), labels=c("Low","Int","High"), include.lowest=TRUE)
  d_cat <- cut(p_dlcn_norm, breaks=c(0, 0.25, 0.75, 1), labels=c("Low","Int","High"), include.lowest=TRUE)

  ev <- y == 1; ne <- y == 0
  # Events: moved to higher category
  ev_up   <- sum(as.numeric(t_cat[ev]) > as.numeric(d_cat[ev]))
  ev_down <- sum(as.numeric(t_cat[ev]) < as.numeric(d_cat[ev]))
  nri_ev  <- (ev_up - ev_down) / sum(ev)
  # Non-events: moved to lower category
  ne_down <- sum(as.numeric(t_cat[ne]) < as.numeric(d_cat[ne]))
  ne_up   <- sum(as.numeric(t_cat[ne]) > as.numeric(d_cat[ne]))
  nri_ne  <- (ne_down - ne_up) / sum(ne)
  nri_total <- nri_ev + nri_ne

  # Continuous NRI
  cnri_ev <- mean(p_tudor[ev] > p_dlcn_norm[ev]) - mean(p_tudor[ev] < p_dlcn_norm[ev])
  cnri_ne <- mean(p_tudor[ne] < p_dlcn_norm[ne]) - mean(p_tudor[ne] > p_dlcn_norm[ne])
  cnri_total <- cnri_ev + cnri_ne

  # IDI
  idi <- (mean(p_tudor[ev]) - mean(p_tudor[ne])) - (mean(p_dlcn_norm[ev]) - mean(p_dlcn_norm[ne]))

  cat(sprintf("\n  %s:\n", ho))
  cat(sprintf("    Categorical NRI: events=%+.3f, non-events=%+.3f, total=%+.3f\n",
              nri_ev, nri_ne, nri_total))
  cat(sprintf("    Continuous NRI:  events=%+.3f, non-events=%+.3f, total=%+.3f\n",
              cnri_ev, cnri_ne, cnri_total))
  cat(sprintf("    IDI: %+.4f\n", idi))
  cat(sprintf("    Events reclassified: up=%d down=%d | Non-events: down=%d up=%d\n",
              ev_up, ev_down, ne_down, ne_up))

  nri_results[[ho]] <- data.frame(
    Cohort = ho, NRI_cat = round(nri_total, 4),
    NRI_events = round(nri_ev, 4), NRI_nonevents = round(nri_ne, 4),
    cNRI = round(cnri_total, 4), IDI = round(idi, 4),
    stringsAsFactors = FALSE
  )
}

# ══════════════════════════════════════════════════════════════════════════════
#  PART 6: SUBGROUP ANALYSES
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== PART 6: Subgroup Analyses ===\n")

# Use Wales as primary subgroup cohort (largest, most complete)
d_sg <- loco_preds[["Wales"]]
sg_results <- list()

subgroups <- list(
  "Male"       = d_sg$sex == 1,
  "Female"     = d_sg$sex == 0,
  "Age <40"    = d_sg$age < 40,
  "Age 40-60"  = d_sg$age >= 40 & d_sg$age < 60,
  "Age >60"    = d_sg$age >= 60,
  "On statin"  = d_sg$on_statin == 1,
  "Statin-free"= d_sg$on_statin == 0,
  "Index case"   = d_sg$index_case == 1,
  "Cascade"      = d_sg$index_case == 0
)

cat("\n  Wales subgroups:\n")
for (sg in names(subgroups)) {
  m <- subgroups[[sg]]
  if (sum(m) < 50 || length(unique(d_sg$fh[m])) < 2) next
  r <- roc(d_sg$fh[m], d_sg$pred[m], quiet=TRUE)
  ci <- ci.auc(r, quiet=TRUE)
  cat(sprintf("    %-15s AUC=%.3f [%.3f-%.3f] n=%d FH+=%d\n",
              sg, as.numeric(auc(r)), ci[1], ci[3], sum(m), sum(d_sg$fh[m])))
  sg_results[[length(sg_results)+1]] <- data.frame(
    Subgroup = sg, AUC = round(as.numeric(auc(r)),4),
    CI_lo = round(ci[1],4), CI_hi = round(ci[3],4),
    N = sum(m), N_FH = sum(d_sg$fh[m]))
}

# Gene-specific
cat("\n  Gene-specific (FH+ only vs all FH-):\n")
for (g in c("LDLR","APOB","PCSK9")) {
  fg <- d_sg[!is.na(d_sg$gene) & d_sg$gene == g & d_sg$fh == 1, ]
  ng <- d_sg[d_sg$fh == 0, ]
  if (nrow(fg) < 5) next
  cg <- rbind(fg, ng)
  rg <- roc(cg$fh, cg$pred, quiet=TRUE); cig <- ci.auc(rg, quiet=TRUE)
  cat(sprintf("    %-15s AUC=%.3f [%.3f-%.3f] FH+=%d\n",
              g, as.numeric(auc(rg)), cig[1], cig[3], nrow(fg)))
  sg_results[[length(sg_results)+1]] <- data.frame(
    Subgroup = g, AUC = round(as.numeric(auc(rg)),4),
    CI_lo = round(cig[1],4), CI_hi = round(cig[3],4),
    N = nrow(cg), N_FH = nrow(fg))
}

# Interaction P-values (sex × model, statin × model)
cat("\n  Interaction tests:\n")
int_sex <- glm(fh ~ pred * sex, data = d_sg, family = binomial)
p_int_sex <- summary(int_sex)$coefficients["pred:sex", "Pr(>|z|)"]
cat(sprintf("    Sex interaction:    P=%s\n", formatC(p_int_sex, format="e", digits=2)))

int_statin <- glm(fh ~ pred * on_statin, data = d_sg, family = binomial)
p_int_statin <- summary(int_statin)$coefficients["pred:on_statin", "Pr(>|z|)"]
cat(sprintf("    Statin interaction: P=%s\n", formatC(p_int_statin, format="e", digits=2)))

# ══════════════════════════════════════════════════════════════════════════════
#  PART 7: ApoB/LDL-C DISCORDANCE (FH+ vs FH-)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== PART 7: ApoB/LDL-C Discordance ===\n")

# South Wales has ApoB data
d_ab <- sw[!is.na(sw$apob) & !is.na(sw$ldl_ut), ]
cat(sprintf("  SW with ApoB: %d (FH+=%d, FH-=%d)\n", nrow(d_ab), sum(d_ab$fh), sum(d_ab$fh==0)))

if (nrow(d_ab) > 50) {
  cat("\n  ApoB/LDL ratio:\n")
  cat(sprintf("    FH+: mean=%.3f, sd=%.3f, median=%.3f\n",
              mean(d_ab$apob_ldl[d_ab$fh==1], na.rm=TRUE),
              sd(d_ab$apob_ldl[d_ab$fh==1], na.rm=TRUE),
              median(d_ab$apob_ldl[d_ab$fh==1], na.rm=TRUE)))
  cat(sprintf("    FH-: mean=%.3f, sd=%.3f, median=%.3f\n",
              mean(d_ab$apob_ldl[d_ab$fh==0], na.rm=TRUE),
              sd(d_ab$apob_ldl[d_ab$fh==0], na.rm=TRUE),
              median(d_ab$apob_ldl[d_ab$fh==0], na.rm=TRUE)))
  wt <- wilcox.test(apob_ldl ~ fh, data = d_ab)
  cat(sprintf("    Wilcoxon P = %s\n", formatC(wt$p.value, format="e", digits=2)))

  cat("\n  nonHDL - LDL gap:\n")
  cat(sprintf("    FH+: mean=%.3f, sd=%.3f\n",
              mean(d_ab$nhdl_ldl_gap[d_ab$fh==1], na.rm=TRUE),
              sd(d_ab$nhdl_ldl_gap[d_ab$fh==1], na.rm=TRUE)))
  cat(sprintf("    FH-: mean=%.3f, sd=%.3f\n",
              mean(d_ab$nhdl_ldl_gap[d_ab$fh==0], na.rm=TRUE),
              sd(d_ab$nhdl_ldl_gap[d_ab$fh==0], na.rm=TRUE)))
  wt2 <- wilcox.test(nhdl_ldl_gap ~ fh, data = d_ab)
  cat(sprintf("    Wilcoxon P = %s\n", formatC(wt2$p.value, format="e", digits=2)))

  # nonHDL-LDL > 0.8 enrichment
  d_ab$nhdl_high <- !is.na(d_ab$nhdl_ldl_gap) & d_ab$nhdl_ldl_gap > 0.8
  t_nhdl <- table(d_ab$fh, d_ab$nhdl_high)
  cat("\n  nonHDL-LDL > 0.8:\n")
  print(t_nhdl)
  if (all(dim(t_nhdl) == c(2,2))) {
    ft <- fisher.test(t_nhdl)
    cat(sprintf("    Fisher OR=%.2f, P=%s\n", ft$estimate, formatC(ft$p.value, format="e", digits=2)))
  }

  # ApoB discordance plot
  pdf(file.path(FIG, "apob_ldl_discordance.pdf"), width=5, height=4)
  par(mar=c(4,4,2,1), family="sans")
  plot(d_ab$ldl_ut[d_ab$fh==0], d_ab$apob[d_ab$fh==0],
       pch=16, cex=0.6, col=adjustcolor("#4393C3", 0.4),
       xlab="LDL-C untreated (mmol/L)", ylab="ApoB (g/L)",
       xlim=c(2, 15), ylim=c(0.5, 3.5), main="ApoB/LDL-C Discordance",
       cex.main=0.9, cex.lab=0.8, cex.axis=0.7)
  points(d_ab$ldl_ut[d_ab$fh==1], d_ab$apob[d_ab$fh==1],
         pch=17, cex=0.7, col=adjustcolor("#D6604D", 0.6))
  abline(lm(apob ~ ldl_ut, data=d_ab[d_ab$fh==0,]), col="#4393C3", lwd=2, lty=2)
  abline(lm(apob ~ ldl_ut, data=d_ab[d_ab$fh==1,]), col="#D6604D", lwd=2)
  legend("topleft", legend=c("FH+","FH-"), pch=c(17,16),
         col=c("#D6604D","#4393C3"), cex=0.7, bty="n")
  dev.off()
}

# UKB ApoB analysis
d_ukb <- loco_preds[["UKB"]]
ukb_ab <- d_ukb[!is.na(d_ukb$apob) & !is.na(d_ukb$ldl_ut), ]
if (nrow(ukb_ab) > 100) {
  cat("\n  UKB ApoB/LDL ratio:\n")
  cat(sprintf("    FH+: mean=%.3f (n=%d)\n",
              mean(ukb_ab$apob_ldl[ukb_ab$fh==1], na.rm=TRUE), sum(ukb_ab$fh==1 & !is.na(ukb_ab$apob_ldl))))
  cat(sprintf("    FH-: mean=%.3f (n=%d)\n",
              mean(ukb_ab$apob_ldl[ukb_ab$fh==0], na.rm=TRUE), sum(ukb_ab$fh==0 & !is.na(ukb_ab$apob_ldl))))
  wt3 <- wilcox.test(apob_ldl ~ fh, data = ukb_ab)
  cat(sprintf("    Wilcoxon P = %s\n", formatC(wt3$p.value, format="e", digits=2)))
}

# ══════════════════════════════════════════════════════════════════════════════
#  PART 8: NMR LIPOPROTEIN PROFILING (UKB 5-fold CV)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== PART 8: NMR Analysis ===\n")

nmr_files <- paste0(TD, "/ukb_nmr_", c("a","b","c","d","e"), ".csv")
nmr_all <- NULL
for (nf in nmr_files) {
  if (file.exists(nf)) {
    tmp <- fread(nf, showProgress=FALSE)
    if (is.null(nmr_all)) nmr_all <- tmp
    else nmr_all <- merge(nmr_all, tmp, by="participant.eid", all=TRUE)
  }
}

if (!is.null(nmr_all)) {
  ukb_nmr <- as.data.frame(merge(as.data.table(ukb), nmr_all,
                                  by.x="participant.eid", by.y="participant.eid", all.x=TRUE))

  nk <- grep("^participant\\.p23", names(ukb_nmr), value=TRUE)
  nk <- nk[nk %in% names(ukb_nmr)]
  cat(sprintf("  NMR features available: %d\n", length(nk)))

  ef <- c(features, nk)
  ukb_cc <- ukb_nmr[complete.cases(ukb_nmr[, c("fh", ef)]), ]
  cat(sprintf("  Complete NMR cases: %d (FH+=%d)\n", nrow(ukb_cc), sum(ukb_cc$fh)))

  if (nrow(ukb_cc) > 500 && sum(ukb_cc$fh) >= 20) {
    n <- nrow(ukb_cc); folds <- sample(rep(1:5, length.out=n))
    pn <- pe <- numeric(n)

    for (k in 1:5) {
      cat(sprintf("  NMR fold %d/5...\n", k))
      ti <- folds != k; vi <- folds == k
      cn <- cv.glmnet(as.matrix(ukb_cc[ti, ef]), ukb_cc$fh[ti],
                      family="binomial", alpha=0.5, nfolds=5,
                      type.measure="auc", standardize=TRUE)
      pn[vi] <- as.numeric(predict(cn, newx=as.matrix(ukb_cc[vi, ef]),
                                    s="lambda.min", type="response"))
      ce <- cv.glmnet(as.matrix(ukb_cc[ti, features]), ukb_cc$fh[ti],
                      family="binomial", alpha=0.5, nfolds=5,
                      type.measure="auc", standardize=TRUE)
      pe[vi] <- as.numeric(predict(ce, newx=as.matrix(ukb_cc[vi, features]),
                                    s="lambda.min", type="response"))
    }

    rn <- roc(ukb_cc$fh, pn, quiet=TRUE); cin <- ci.auc(rn, quiet=TRUE)
    re <- roc(ukb_cc$fh, pe, quiet=TRUE); cie <- ci.auc(re, quiet=TRUE)
    dl <- tryCatch(roc.test(re, rn, method="delong", paired=TRUE), error=function(e) NULL)

    cat(sprintf("\n  NMR-TUDOR:  AUC=%.4f [%.4f-%.4f]\n", as.numeric(auc(rn)), cin[1], cin[3]))
    cat(sprintf("  Base-TUDOR: AUC=%.4f [%.4f-%.4f]\n", as.numeric(auc(re)), cie[1], cie[3]))
    if (!is.null(dl))
      cat(sprintf("  DeLong: dAUC=%+.4f, p=%s\n",
                  as.numeric(auc(rn)) - as.numeric(auc(re)),
                  formatC(dl$p.value, format="e", digits=2)))
  }
} else { cat("  NMR files not found\n") }

# ══════════════════════════════════════════════════════════════════════════════
#  PART 9: Lp(a), ApoA, nonHDL-LDL Analysis (UKB)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== PART 9: Lp(a), ApoA, Biomarker Profiles ===\n")

# Load processed UKB for Lp(a), ApoA
ukb_proc_file <- file.path(DD, "calon_ukb_analysis_ready.csv")
if (!file.exists(ukb_proc_file))
  ukb_proc_file <- file.path(DD, "output", "calon_ukb_analysis_ready.csv")

if (file.exists(ukb_proc_file)) {
  ukb_proc <- read.csv(ukb_proc_file, stringsAsFactors=FALSE)
  names(ukb_proc) <- make.names(names(ukb_proc), unique=TRUE)

  # Merge with UKB FH status
  ukb_bio <- merge(ukb[, c("participant.eid","fh","ldl_ut","hdl","tg","apob","on_statin","gene")],
                   ukb_proc[, intersect(c("eid","lpa","apoa","inv_apoa","hba1c","glucose","crp",
                                           "creatinine","alt","cystatin_c","mean_cimt"),
                                         names(ukb_proc))],
                   by.x="participant.eid", by.y="eid", all.x=TRUE)
  ukb_bio <- ukb_bio[!duplicated(ukb_bio$participant.eid), ]

  biomarkers <- c("lpa","apoa","hba1c","glucose","crp","creatinine","alt","cystatin_c","mean_cimt")
  biomarkers <- biomarkers[biomarkers %in% names(ukb_bio)]

  cat("\n  Biomarker comparison (FH+ vs FH-):\n")
  cat(sprintf("  %-15s %8s %8s %8s %8s %10s\n", "Biomarker", "FH+_mean", "FH+_sd", "FH-_mean", "FH-_sd", "P_value"))

  for (bm in biomarkers) {
    v_pos <- ukb_bio[[bm]][ukb_bio$fh == 1]
    v_neg <- ukb_bio[[bm]][ukb_bio$fh == 0]
    v_pos <- v_pos[!is.na(v_pos)]; v_neg <- v_neg[!is.na(v_neg)]
    if (length(v_pos) < 10 || length(v_neg) < 10) next
    wt <- wilcox.test(v_pos, v_neg)
    cat(sprintf("  %-15s %8.3f %8.3f %8.3f %8.3f %10s\n",
                bm, mean(v_pos), sd(v_pos), mean(v_neg), sd(v_neg),
                formatC(wt$p.value, format="e", digits=2)))
  }

  # Lp(a) > 143 nmol/L enrichment
  if ("lpa" %in% names(ukb_bio)) {
    ukb_bio$lpa_high <- !is.na(ukb_bio$lpa) & ukb_bio$lpa > 143
    t_lpa <- table(ukb_bio$fh, ukb_bio$lpa_high)
    cat("\n  Lp(a) > 143 nmol/L:\n")
    if (all(dim(t_lpa) == c(2,2))) {
      print(t_lpa)
      ft_lpa <- fisher.test(t_lpa)
      cat(sprintf("    Fisher OR=%.2f, P=%s\n", ft_lpa$estimate, formatC(ft_lpa$p.value, format="e", digits=2)))
    }
  }
} else { cat("  UKB processed file not found\n") }

# ══════════════════════════════════════════════════════════════════════════════
#  PART 10: PUBLICATION FIGURES
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== PART 10: Publication Figures ===\n")

# Figure 1: ROC curves (3-panel, one per cohort)
pdf(file.path(FIG, "Figure1_ROC_LOCO.pdf"), width=10, height=3.5)
par(mfrow=c(1,3), mar=c(4,4,2,1), family="sans")

for (ho in cohorts) {
  d <- loco_preds[[ho]]
  r_tudor <- roc(d$fh, d$pred, quiet=TRUE)
  ci_t <- ci.auc(r_tudor, quiet=TRUE)

  # eDLCN
  dlcn_ok <- !is.na(d$dlcn)
  r_dlcn <- if (sum(dlcn_ok) > 50) roc(d$fh[dlcn_ok], d$dlcn[dlcn_ok], quiet=TRUE) else NULL

  # LDL alone
  r_ldl <- roc(d$fh, d$ldl_ut, quiet=TRUE)

  # Trig_filter alone
  r_tf <- roc(d$fh, d$trig_filter, quiet=TRUE)

  plot(r_tudor, col="#2166AC", lwd=2.5, main=ho,
       cex.main=0.9, cex.lab=0.8, cex.axis=0.7, legacy.axes=TRUE)
  if (!is.null(r_dlcn)) plot(r_dlcn, add=TRUE, col="#B2182B", lwd=2, lty=2)
  plot(r_ldl, add=TRUE, col="#66BD63", lwd=1.5, lty=3)
  plot(r_tf, add=TRUE, col="#F46D43", lwd=1.5, lty=4)

  legend("bottomright",
         legend=c(sprintf("TUDOR (%.3f)", auc(r_tudor)),
                  sprintf("eDLCN (%.3f)", ifelse(!is.null(r_dlcn), auc(r_dlcn), NA)),
                  sprintf("LDL-C (%.3f)", auc(r_ldl)),
                  sprintf("Trig_Filter (%.3f)", auc(r_tf))),
         col=c("#2166AC","#B2182B","#66BD63","#F46D43"),
         lty=c(1,2,3,4), lwd=c(2.5,2,1.5,1.5), cex=0.5, bty="n")
}
dev.off()

# Figure 3: Subgroup forest plot
if (length(sg_results) > 0) {
  sg_df <- do.call(rbind, sg_results)

  pdf(file.path(FIG, "Figure3_subgroup_forest.pdf"), width=6, height=5)
  par(mar=c(4,10,2,2), family="sans")
  n_sg <- nrow(sg_df)
  plot(sg_df$AUC, n_sg:1, xlim=c(0.55, 1.0), pch=19, cex=1.2,
       yaxt="n", xlab="AUC", ylab="", main="Subgroup Analysis — Wales",
       cex.main=0.9, cex.lab=0.8, cex.axis=0.7)
  segments(sg_df$CI_lo, n_sg:1, sg_df$CI_hi, n_sg:1, lwd=2)
  axis(2, at=n_sg:1, labels=sprintf("%s (n=%d)", sg_df$Subgroup, sg_df$N),
       las=1, cex.axis=0.6)
  abline(v=0.78, lty=2, col="grey60")
  text(sg_df$AUC + 0.01, n_sg:1, labels=sprintf("%.3f", sg_df$AUC), cex=0.55, pos=4)
  dev.off()
}

cat("  Figures saved to:", FIG, "\n")

# ══════════════════════════════════════════════════════════════════════════════
#  PART 11: SUMMARY TABLES
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== PART 11: Summary Tables ===\n")

# TABLE 1: Baseline characteristics
cat("\n  TABLE 1: Baseline Characteristics\n")
cat("  ──────────────────────────────────────────────────────────\n")
cat(sprintf("  %-20s %12s %12s %12s\n", "Variable", "South Wales", "Wales", "UKB"))

for (v in c("age","ldl_ut","hdl","tg","trig_filter")) {
  vals <- sapply(list(sw, wales, ukb_df), function(df) {
    x <- df[[v]]
    sprintf("%.1f (%.1f)", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
  })
  cat(sprintf("  %-20s %12s %12s %12s\n", v, vals[1], vals[2], vals[3]))
}
for (v in c("sex","on_statin")) {
  vals <- sapply(list(sw, wales, ukb_df), function(df) {
    sprintf("%d (%.1f%%)", sum(df[[v]] == 1, na.rm=TRUE), 100*mean(df[[v]], na.rm=TRUE))
  })
  cat(sprintf("  %-20s %12s %12s %12s\n", v, vals[1], vals[2], vals[3]))
}
cat(sprintf("  %-20s %12s %12s %12s\n", "FH prevalence",
            sprintf("%.1f%%", 100*mean(sw$fh)),
            sprintf("%.1f%%", 100*mean(wales$fh)),
            sprintf("%.1f%%", 100*mean(ukb_df$fh))))

# TABLE 2: Model performance
cat("\n\n  TABLE 2: Diagnostic Performance\n")
cat("  ──────────────────────────────────────────────────────────\n")
cat(sprintf("  %-12s %8s %15s %8s %8s %8s %8s %8s\n",
            "Cohort","Model","AUC (95% CI)","Sens","Spec","PPV","NPV","Brier"))

for (ho in cohorts) {
  d <- loco_preds[[ho]]
  for (model_name in c("TUDOR","eDLCN")) {
    if (model_name == "TUDOR") {
      p <- d$pred; y <- d$fh
    } else {
      ok <- !is.na(d$dlcn)
      if (sum(ok) < 50) next
      p <- d$dlcn[ok] / max(d$dlcn[ok]); y <- d$fh[ok]
    }
    r <- roc(y, p, quiet=TRUE); ci <- ci.auc(r, quiet=TRUE)
    co <- coords(r, "best", best.method="youden")
    prev <- mean(y)
    ppv <- (co$sensitivity * prev) / (co$sensitivity * prev + (1-co$specificity)*(1-prev))
    npv <- (co$specificity * (1-prev)) / ((1-co$sensitivity)*prev + co$specificity*(1-prev))
    brier <- mean((p - y)^2)

    cat(sprintf("  %-12s %8s %6.3f (%5.3f-%5.3f) %6.1f%% %6.1f%% %6.1f%% %6.1f%% %6.4f\n",
                ho, model_name, as.numeric(auc(r)), ci[1], ci[3],
                co$sensitivity*100, co$specificity*100, ppv*100, npv*100, brier))
  }
}

# TABLE 3: NRI / IDI
if (length(nri_results) > 0) {
  nri_df <- do.call(rbind, nri_results)
  cat("\n\n  TABLE 3: Reclassification (TUDOR vs eDLCN)\n")
  cat("  ──────────────────────────────────────────────────────────\n")
  print(nri_df)
}

# Save all results
all_preds <- do.call(rbind, loco_preds)
write.csv(all_preds, file.path(OUT, "loco_predictions_complete.csv"), row.names=FALSE)
if (length(sg_results) > 0) {
  sg_df <- do.call(rbind, sg_results)
  write.csv(sg_df, file.path(OUT, "subgroup_results.csv"), row.names=FALSE)
}
if (length(nri_results) > 0)
  write.csv(nri_df, file.path(OUT, "nri_idi_results.csv"), row.names=FALSE)

cat(sprintf("\n  All results saved to: %s/\n", OUT))
cat(sprintf("  Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("================================================================\n")

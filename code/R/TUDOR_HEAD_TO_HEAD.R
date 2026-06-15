################################################################################
#
#  TUDOR: HEAD-TO-HEAD COMPARATOR ANALYSIS
#
#  Paired DeLong comparisons on IDENTICAL patients:
#    - TUDOR-ENET vs eDLCN
#    - TUDOR-RF vs eDLCN
#    - TUDOR-XGB vs eDLCN
#    - TUDOR-ENET vs LDL-C alone
#    - TUDOR-ENET vs Trig_Filter alone
#    - TUDOR-RF vs TUDOR-ENET
#    - eDLCN vs LDL-C alone
#
#  Per cohort + pooled. Bootstrap 95% CIs. Youden thresholds.
#  Sens/Spec/PPV/NPV at matched specificity. NRI/IDI.
#
#  Author: Dr Nader Genedy | 2026-03-27
#
################################################################################

cat("\n================================================================\n")
cat("  TUDOR: HEAD-TO-HEAD COMPARATOR ANALYSIS\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

set.seed(42)

suppressPackageStartupMessages({
  library(readxl); library(haven); library(pROC)
  library(glmnet); library(randomForest); library(xgboost)
  library(data.table)
})

# ── HELPERS ──────────────────────────────────────────────────────────────────

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

statin_factors <- c(atorvastatin=0.38, simvastatin=0.35, rosuvastatin=0.34,
                    pravastatin=0.25, fluvastatin=0.22, lovastatin=0.25,
                    ezetimibe=0.18)

get_red <- function(tx) {
  if (is.na(tx) || tx == "" || tx == "NaN") return(0)
  tx_low <- tolower(as.character(tx))
  for (nm in names(statin_factors))
    if (grepl(nm, tx_low, fixed = TRUE)) return(statin_factors[[nm]])
  if (grepl("statin", tx_low)) return(0.30)
  return(0)
}

rev_ldl <- function(ldl, rf) ifelse(!is.na(ldl) & !is.na(rf) & rf > 0, ldl/(1-rf), ldl)

parse_dt <- function(x) {
  x <- as.character(x)
  d <- suppressWarnings(as.Date(x, format="%d-%m-%Y"))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, format="%d/%m/%Y"))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, format="%Y-%m-%d"))
  d
}

DATA_DIR  <- "C:/Users/nader/Downloads/calon_ukb_pipeline"
TUDOR_DIR <- "D:/alphafold_backup/tudor_packup"
OUT_DIR   <- file.path(DATA_DIR, "tudor_loco_output")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

features <- c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex")

# ── LOAD COHORT 1: SOUTH WALES ──────────────────────────────────────────────

cat("Loading cohorts...\n")

pos <- as.data.frame(read_excel(file.path(DATA_DIR, "mutation positive group.xlsx")))
pos$fh <- 1
pos$dob_d <- as.Date(sapply(pos$BirthDate, parse_dt), origin="1970-01-01")
pos$meas_d <- as.Date(sapply(pos$MeasurementDate.1, parse_dt), origin="1970-01-01")
pos$age <- as.numeric(difftime(pos$meas_d, pos$dob_d, units="days"))/365.25
pos$sex <- ifelse(pos$Gender %in% c("M","Male"), 1, 0)
pos$hdl <- safe_num(pos$HDL.1); pos$tg <- safe_num(pos$TRG.1); pos$ldl_m <- safe_num(pos$LDL.1)
pos$tx1 <- as.character(pos$Treatment1.1)
pos$is_tx <- !is.na(pos$tx1) & pos$tx1 != "" & pos$tx1 != "NA"
pos$rf <- ifelse(pos$is_tx, sapply(pos$tx1, get_red), 0)
pos$ldl_ut <- ifelse(pos$is_tx, rev_ldl(pos$ldl_m, pos$rf), pos$ldl_m)
pos$trig_filter <- ifelse(!is.na(pos$ldl_ut)&!is.na(pos$tg)&pos$tg>=0, pmin(pos$ldl_ut/(pos$tg+0.1),50), NA)
ca <- if("Corneal Arcus"%in%names(pos)) "Corneal Arcus" else "CornealArcus"
pos$tendon_xanth <- ifelse(pos$TendonXanthomata %in% c(1,"1","Yes","TRUE"), 1, 0)
pos$corneal_arcus <- ifelse(pos[[ca]] %in% c(1,"1","Yes","TRUE"), 1, 0)
pos$tendon_xanth[is.na(pos$tendon_xanth)] <- 0; pos$corneal_arcus[is.na(pos$corneal_arcus)] <- 0
pos$dlcn <- safe_num(pos$GenoTypingScore); pos$cohort <- "SouthWales"

neg <- as.data.frame(read_excel(file.path(DATA_DIR, "Mutation negative control updated.xlsx")))
neg$fh <- 0; neg$age <- safe_num(neg$`age at the result`)
neg$sex <- ifelse(neg$Gender %in% c("M","Male"), 1, 0)
neg$hdl <- safe_num(neg$HDL.1); neg$tg <- safe_num(neg$TRG.1); neg$ldl_m <- safe_num(neg$LDL.1)
neg$tx1 <- as.character(neg$Treatment1.1)
neg$is_tx <- !is.na(neg$tx1) & neg$tx1 != "" & neg$tx1 != "NA"
neg$rf <- ifelse(neg$is_tx, sapply(neg$tx1, get_red), 0)
neg$ldl_ut <- ifelse(neg$is_tx, rev_ldl(neg$ldl_m, neg$rf), neg$ldl_m)
neg$trig_filter <- ifelse(!is.na(neg$ldl_ut)&!is.na(neg$tg)&neg$tg>=0, pmin(neg$ldl_ut/(neg$tg+0.1),50), NA)
ca2 <- if("Corneal Arcus"%in%names(neg)) "Corneal Arcus" else "CornealArcus"
neg$tendon_xanth <- ifelse(neg$TendonXanthomata %in% c(1,"1","Yes","TRUE"), 1, 0)
neg$corneal_arcus <- ifelse(neg[[ca2]] %in% c(1,"1","Yes","TRUE"), 1, 0)
neg$tendon_xanth[is.na(neg$tendon_xanth)] <- 0; neg$corneal_arcus[is.na(neg$corneal_arcus)] <- 0
neg$dlcn <- safe_num(neg$GenoTypingScore); neg$cohort <- "SouthWales"

cols <- c("fh","age","sex","hdl","tg","ldl_m","ldl_ut","trig_filter",
          "tendon_xanth","corneal_arcus","dlcn","cohort")
sw <- rbind(pos[,cols], neg[,cols])

# ── LOAD COHORT 2: WALES PASS ───────────────────────────────────────────────

pass <- as.data.frame(read_sav(file.path(DATA_DIR, "PASS_wrong_dob.sav")))
pass$fh <- ifelse(pass$Positive1=="Yes",1, ifelse(pass$Positive1=="No",0,NA))
pass$age <- safe_num(pass$BMI_AGE); pass$sex <- ifelse(pass$Gender=="M",1,0)
pass$hdl <- safe_num(pass$HDL.1); pass$tg <- safe_num(pass$TRG.1); pass$ldl_m <- safe_num(pass$LDL.1)
pass$tx1 <- as.character(pass$Treatment1.1)
pass$is_tx <- !is.na(pass$tx1) & pass$tx1 != "" & pass$tx1 != "NA"
pass$rf2 <- ifelse(pass$is_tx, sapply(pass$tx1, get_red), 0)
pass$ldl_ut <- ifelse(pass$is_tx, rev_ldl(pass$ldl_m, pass$rf2), pass$ldl_m)
pass$trig_filter <- ifelse(!is.na(pass$ldl_ut)&!is.na(pass$tg)&pass$tg>=0, pmin(pass$ldl_ut/(pass$tg+0.1),50), NA)
pass$tendon_xanth <- ifelse(safe_num(pass$TendonXanthomata)==1,1,0)
pass$corneal_arcus <- ifelse(safe_num(pass$CornealArcus)==1,1,0)
pass$tendon_xanth[is.na(pass$tendon_xanth)] <- 0; pass$corneal_arcus[is.na(pass$corneal_arcus)] <- 0
pass$dlcn <- safe_num(pass$GenoTypingScore); pass$cohort <- "Wales"
wales <- pass[!is.na(pass$fh), cols]

# ── LOAD COHORT 3: UKB LIPID CLINIC ─────────────────────────────────────────

ukb <- as.data.frame(fread(file.path(TUDOR_DIR, "TUDOR_UKB_Features.csv"), showProgress=FALSE))
ukb$fh <- safe_num(ukb$is_fh_genetic); ukb$fh[is.na(ukb$fh)] <- 0
ukb$age <- safe_num(ukb$Age_at_LDL1); ukb$sex <- safe_num(ukb$Gender_num)
ukb$hdl <- safe_num(ukb$HDL.1); ukb$tg <- safe_num(ukb$TRG.1)
ukb$ldl_m <- safe_num(ukb$LDL_treated); ukb$ldl_ut <- safe_num(ukb$LDL_untreated)
ukb$trig_filter <- safe_num(ukb$Trig_Filter)
ukb$on_statin <- ifelse(safe_num(ukb$reduction_factor)>0,1,0); ukb$on_statin[is.na(ukb$on_statin)] <- 0
ukb$tc <- safe_num(ukb$CHOL); ukb$non_hdl <- ifelse(!is.na(ukb$tc)&!is.na(ukb$hdl), ukb$tc-ukb$hdl, NA)

# Lipid clinic filter
ukb$lc <- (!is.na(ukb$tc) & ukb$tc > 7.9) |
          (!is.na(ukb$ldl_ut) & ukb$ldl_ut > 4.9) |
          (!is.na(ukb$non_hdl) & ukb$non_hdl > 5.9) |
          (ukb$on_statin == 1)
ukb <- ukb[ukb$lc == TRUE, ]

# eDLCN
ukb$dlcn <- ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut>=8.5, 8,
            ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut>=6.5, 5,
            ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut>=5.0, 3,
            ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut>=4.0, 1, 0))))
ukb$tendon_xanth <- 0; ukb$corneal_arcus <- 0; ukb$cohort <- "UKB"
ukb_df <- ukb[, cols]

all_data <- rbind(sw, wales, ukb_df)
cat(sprintf("  Total: %d (SW=%d, Wales=%d, UKB=%d)\n",
            nrow(all_data), nrow(sw), nrow(wales), nrow(ukb_df)))

# ══════════════════════════════════════════════════════════════════════════════
# HEAD-TO-HEAD: For each held-out cohort, train on other 2 and compare
# ══════════════════════════════════════════════════════════════════════════════

cohorts <- c("SouthWales", "Wales", "UKB")
h2h_results <- list()
delong_results <- list()

for (held_out in cohorts) {
  cat(sprintf("\n════════════════════════════════════════════\n"))
  cat(sprintf("  HELD-OUT: %s\n", held_out))
  cat(sprintf("════════════════════════════════════════════\n"))

  train <- all_data[all_data$cohort != held_out, ]
  test  <- all_data[all_data$cohort == held_out, ]

  # Complete cases for ML features
  train_cc <- train[complete.cases(train[, c("fh", features)]), ]
  test_cc  <- test[complete.cases(test[, c("fh", features)]), ]

  if (nrow(test_cc) < 50 || sum(test_cc$fh == 1) < 10) {
    cat("  Insufficient test data — skipping\n"); next
  }

  X_tr <- as.matrix(train_cc[, features]); y_tr <- train_cc$fh
  X_te <- as.matrix(test_cc[, features]);  y_te <- test_cc$fh

  cat(sprintf("  Train: %d (FH+=%d) | Test: %d (FH+=%d)\n",
              nrow(train_cc), sum(y_tr), nrow(test_cc), sum(y_te)))

  # ── TRAIN ALL MODELS ──

  # 1. Elastic Net
  cv_enet <- cv.glmnet(X_tr, y_tr, family="binomial", alpha=0.5,
                        nfolds=10, type.measure="auc", standardize=TRUE)
  pred_enet <- as.numeric(predict(cv_enet, newx=X_te, s="lambda.min", type="response"))

  # 2. Random Forest (with downsampling)
  n_pos <- sum(y_tr==1); n_neg <- sum(y_tr==0)
  if (n_neg > 50000) {
    idx_ds <- sort(c(sample(which(y_tr==0), min(n_neg, 5*n_pos)), which(y_tr==1)))
    X_rf <- X_tr[idx_ds,]; y_rf <- y_tr[idx_ds]
  } else { X_rf <- X_tr; y_rf <- y_tr }
  df_rf <- data.frame(X_rf, fh=factor(y_rf, levels=c(0,1)))
  np <- sum(y_rf==1); nn <- sum(y_rf==0)
  rf_fit <- randomForest(fh~., data=df_rf, ntree=200, mtry=max(2,floor(sqrt(ncol(X_rf)))),
                          sampsize=c("0"=min(nn,3*np), "1"=np), strata=df_rf$fh)
  pred_rf <- predict(rf_fit, newdata=data.frame(X_te), type="prob")[,"1"]

  # 3. XGBoost
  dtrain <- xgb.DMatrix(data=X_tr, label=y_tr)
  params <- list(objective="binary:logistic", eval_metric="auc", max_depth=3,
                 eta=0.05, min_child_weight=10, subsample=0.8, colsample_bytree=0.8,
                 gamma=1, lambda=5, alpha=1, scale_pos_weight=n_neg/n_pos)
  cv_xgb <- xgb.cv(params=params, data=dtrain, nrounds=500, nfold=5,
                     early_stopping_rounds=50, verbose=0, stratified=TRUE)
  best_n <- cv_xgb$best_iteration; if(is.null(best_n)||is.na(best_n)) best_n <- 100
  xgb_fit <- xgb.train(params=params, data=dtrain, nrounds=best_n, verbose=0)
  pred_xgb <- predict(xgb_fit, xgb.DMatrix(data=X_te))

  # 4. Comparators (on same test patients)
  pred_ldl <- test_cc$ldl_ut        # LDL-C alone
  pred_tf  <- test_cc$trig_filter   # Trig_Filter alone
  pred_dlcn <- test_cc$dlcn         # eDLCN

  # ── COMPUTE ROC FOR ALL ──

  models <- list(
    "TUDOR-ENET"  = pred_enet,
    "TUDOR-RF"    = pred_rf,
    "TUDOR-XGB"   = pred_xgb,
    "eDLCN"       = pred_dlcn,
    "LDL-C alone" = pred_ldl,
    "Trig_Filter"  = pred_tf
  )

  rocs <- list()
  cat("\n  ── Discrimination ──\n")
  for (nm in names(models)) {
    p <- models[[nm]]
    valid <- !is.na(p) & !is.na(y_te)
    if (sum(valid) < 50 || length(unique(y_te[valid])) < 2) {
      cat(sprintf("    %s: insufficient data\n", nm)); next
    }
    r <- roc(y_te[valid], p[valid], quiet=TRUE)
    ci <- ci.auc(r, quiet=TRUE)
    rocs[[nm]] <- r

    # Youden threshold
    co <- coords(r, "best", best.method="youden")
    sens <- co$sensitivity; spec <- co$specificity
    # PPV / NPV at Youden
    prev <- mean(y_te[valid])
    ppv <- (sens * prev) / (sens * prev + (1-spec) * (1-prev))
    npv <- (spec * (1-prev)) / ((1-sens) * prev + spec * (1-prev))

    cat(sprintf("    %-14s AUC=%.3f (%.3f-%.3f) | Sens=%.1f%% Spec=%.1f%% PPV=%.1f%% NPV=%.1f%%\n",
                nm, as.numeric(auc(r)), ci[1], ci[3],
                sens*100, spec*100, ppv*100, npv*100))

    h2h_results[[length(h2h_results)+1]] <- data.frame(
      HeldOut=held_out, Model=nm,
      AUC=round(as.numeric(auc(r)),4), CI_lo=round(ci[1],4), CI_hi=round(ci[3],4),
      Sensitivity=round(sens,4), Specificity=round(spec,4),
      PPV=round(ppv,4), NPV=round(npv,4),
      N=sum(valid), N_FH=sum(y_te[valid]==1),
      Threshold=round(co$threshold,4)
    )
  }

  # ── DELONG PAIRWISE COMPARISONS ──

  cat("\n  ── DeLong Pairwise Tests ──\n")

  comparisons <- list(
    c("TUDOR-ENET", "eDLCN"),
    c("TUDOR-RF", "eDLCN"),
    c("TUDOR-XGB", "eDLCN"),
    c("TUDOR-RF", "TUDOR-ENET"),
    c("TUDOR-ENET", "LDL-C alone"),
    c("TUDOR-ENET", "Trig_Filter"),
    c("eDLCN", "LDL-C alone"),
    c("TUDOR-RF", "TUDOR-XGB")
  )

  for (comp in comparisons) {
    m1 <- comp[1]; m2 <- comp[2]
    if (!(m1 %in% names(rocs)) || !(m2 %in% names(rocs))) next

    # Need matched patients (same indices)
    p1 <- models[[m1]]; p2 <- models[[m2]]
    valid <- !is.na(p1) & !is.na(p2) & !is.na(y_te)
    if (sum(valid) < 50) next

    r1 <- roc(y_te[valid], p1[valid], quiet=TRUE)
    r2 <- roc(y_te[valid], p2[valid], quiet=TRUE)

    dt <- tryCatch(roc.test(r1, r2, method="delong", paired=TRUE), error=function(e) NULL)
    if (is.null(dt)) next

    auc_diff <- as.numeric(auc(r1)) - as.numeric(auc(r2))
    sig <- ifelse(dt$p.value < 0.001, "***",
           ifelse(dt$p.value < 0.01, "**",
           ifelse(dt$p.value < 0.05, "*", "ns")))

    cat(sprintf("    %s vs %s: ΔAUC=%+.3f, Z=%.2f, p=%s %s\n",
                m1, m2, auc_diff, dt$statistic, formatC(dt$p.value, format="e", digits=2), sig))

    delong_results[[length(delong_results)+1]] <- data.frame(
      HeldOut=held_out, Model1=m1, Model2=m2,
      AUC1=round(as.numeric(auc(r1)),4), AUC2=round(as.numeric(auc(r2)),4),
      Delta_AUC=round(auc_diff,4), Z=round(dt$statistic,3),
      P_value=dt$p.value,
      Significant=sig, N_matched=sum(valid)
    )
  }

  # ── NRI: TUDOR-ENET vs eDLCN ──

  cat("\n  ── NRI / IDI: TUDOR-ENET vs eDLCN ──\n")

  valid_nri <- !is.na(pred_enet) & !is.na(pred_dlcn) & !is.na(y_te)
  if (sum(valid_nri) > 100) {
    y_nri <- y_te[valid_nri]
    p_tudor <- pred_enet[valid_nri]
    p_dlcn  <- pred_dlcn[valid_nri]

    # Normalise eDLCN to probability scale (0-1) for fair comparison
    p_dlcn_norm <- (p_dlcn - min(p_dlcn, na.rm=TRUE)) /
                   (max(p_dlcn, na.rm=TRUE) - min(p_dlcn, na.rm=TRUE) + 0.001)

    # Categorical NRI at median threshold
    thresh <- median(p_tudor)
    tudor_pos <- p_tudor >= thresh
    dlcn_pos  <- p_dlcn_norm >= 0.5  # eDLCN: score >=6 is probable/definite

    # Events
    ev <- y_nri == 1
    ev_up   <- sum(tudor_pos[ev] & !dlcn_pos[ev])
    ev_down <- sum(!tudor_pos[ev] & dlcn_pos[ev])
    nri_ev  <- (ev_up - ev_down) / sum(ev)

    # Non-events
    ne <- y_nri == 0
    ne_down <- sum(!tudor_pos[ne] & dlcn_pos[ne])
    ne_up   <- sum(tudor_pos[ne] & !dlcn_pos[ne])
    nri_ne  <- (ne_down - ne_up) / sum(ne)

    nri_total <- nri_ev + nri_ne

    # IDI
    idi <- (mean(p_tudor[ev]) - mean(p_tudor[ne])) -
           (mean(p_dlcn_norm[ev]) - mean(p_dlcn_norm[ne]))

    cat(sprintf("    NRI(events)=%+.3f, NRI(non-events)=%+.3f, NRI(total)=%+.3f\n",
                nri_ev, nri_ne, nri_total))
    cat(sprintf("    IDI=%+.4f\n", idi))
    cat(sprintf("    Events reclassified: up=%d, down=%d\n", ev_up, ev_down))
    cat(sprintf("    Non-events reclassified: down=%d, up=%d\n", ne_down, ne_up))
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY TABLES
# ══════════════════════════════════════════════════════════════════════════════

cat("\n\n════════════════════════════════════════════════════════════════\n")
cat("  COMPLETE HEAD-TO-HEAD RESULTS\n")
cat("════════════════════════════════════════════════════════════════\n\n")

h2h_df <- do.call(rbind, h2h_results)
cat("=== Discrimination & Clinical Performance ===\n")
print(h2h_df[, c("HeldOut","Model","AUC","CI_lo","CI_hi","Sensitivity","Specificity","PPV","NPV")])

cat("\n=== DeLong Pairwise Comparisons ===\n")
dl_df <- do.call(rbind, delong_results)
print(dl_df[, c("HeldOut","Model1","Model2","AUC1","AUC2","Delta_AUC","P_value","Significant")])

# Save
write.csv(h2h_df, file.path(OUT_DIR, "head_to_head_discrimination.csv"), row.names=FALSE)
write.csv(dl_df, file.path(OUT_DIR, "head_to_head_delong.csv"), row.names=FALSE)

# ── POOLED AVERAGES ──

cat("\n=== Pooled Average AUC Across 3 Folds ===\n")
for (m in unique(h2h_df$Model)) {
  sub <- h2h_df[h2h_df$Model == m, ]
  if (nrow(sub) >= 2) {
    cat(sprintf("  %-14s  Mean=%.3f  Range=[%.3f, %.3f]  Folds=%d\n",
                m, mean(sub$AUC), min(sub$AUC), max(sub$AUC), nrow(sub)))
  }
}

cat(sprintf("\n  Saved: %s/head_to_head_*.csv\n", OUT_DIR))
cat(sprintf("  Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("════════════════════════════════════════════════════════════════\n")

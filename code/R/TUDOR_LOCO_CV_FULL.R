################################################################################
# TUDOR FULL LOCO-CV + AUGMENTATION + SUBGROUPS
# Author: Dr Nader Genedy | 2026-03-27
################################################################################
cat("\n================================================================\n")
cat("  TUDOR FULL LOCO-CV\n  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")
set.seed(42)
suppressPackageStartupMessages({
  library(readxl); library(haven); library(pROC); library(glmnet); library(data.table)
})
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
pd <- function(x) {
  x <- as.character(x)
  d <- suppressWarnings(as.Date(x, format = "%d-%m-%Y"))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, format = "%d/%m/%Y"))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  d
}
DD <- "C:/Users/nader/Downloads/calon_ukb_pipeline"
TD <- "D:/alphafold_backup/tudor_packup"
OD <- file.path(DD, "tudor_loco_output")
dir.create(OD, recursive = TRUE, showWarnings = FALSE)

cat("=== Loading cohorts ===\n")

# SOUTH WALES POSITIVE
pos <- as.data.frame(read_excel(file.path(DD, "mutation positive group.xlsx")))
pos$fh <- 1
pos$dob_d <- as.Date(sapply(pos$BirthDate, pd), origin = "1970-01-01")
pos$meas_d <- as.Date(sapply(pos$MeasurementDate.1, pd), origin = "1970-01-01")
pos$age <- as.numeric(difftime(pos$meas_d, pos$dob_d, units = "days")) / 365.25
pos$sex <- ifelse(pos$Gender %in% c("M", "Male"), 1, 0)
pos$hdl <- safe_num(pos$HDL.1); pos$tg <- safe_num(pos$TRG.1); pos$ldl_m <- safe_num(pos$LDL.1)
pos$tc <- safe_num(pos$TC.1)
pos$rf <- sapply(pos$Treatment1.1, get_red)
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
pos$apob_ldl <- ifelse(!is.na(pos$apob) & !is.na(pos$ldl_ut) & pos$ldl_ut > 0, pos$apob / pos$ldl_ut, NA)
pos$non_hdl <- ifelse(!is.na(pos$tc) & !is.na(pos$hdl), pos$tc - pos$hdl, NA)
pos$nhdl_ldl_gap <- ifelse(!is.na(pos$non_hdl) & !is.na(pos$ldl_ut), pos$non_hdl - pos$ldl_ut, NA)
pos$cumulative_ldl <- ifelse(!is.na(pos$ldl_ut) & !is.na(pos$age), pos$ldl_ut * pos$age, NA)
mut_col_pos <- intersect(c("Mutation (1)", "Mutation1", "Mutation.1"), names(pos))[1]
pos_mut <- if (!is.na(mut_col_pos)) pos[[mut_col_pos]] else NA
pos$gene <- ifelse(grepl("^LDLR", pos_mut), "LDLR",
             ifelse(grepl("^APOB", pos_mut), "APOB",
             ifelse(grepl("^PCSK9", pos_mut), "PCSK9", "Other")))
pos$index_effect <- 1; pos$cohort <- "SouthWales"

# SOUTH WALES NEGATIVE
neg <- as.data.frame(read_excel(file.path(DD, "Mutation negative control updated.xlsx")))
neg$fh <- 0; neg$age <- safe_num(neg[["age at the result"]])
neg$sex <- ifelse(neg$Gender %in% c("M", "Male"), 1, 0)
neg$hdl <- safe_num(neg$HDL.1); neg$tg <- safe_num(neg$TRG.1); neg$ldl_m <- safe_num(neg$LDL.1)
neg$tc <- safe_num(neg$TC.1)
neg$rf <- sapply(neg$Treatment1.1, get_red)
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
# ApoB column in neg is named apoB...2 (disambiguated by readxl)
apob_col_neg <- intersect(c("apoB...2", "apoB", "ApoB"), names(neg))[1]
neg$apob <- if (!is.na(apob_col_neg)) safe_num(neg[[apob_col_neg]]) else NA
neg$apob_ldl <- ifelse(!is.na(neg$apob) & !is.na(neg$ldl_ut) & neg$ldl_ut > 0, neg$apob / neg$ldl_ut, NA)
neg$non_hdl <- ifelse(!is.na(neg$tc) & !is.na(neg$hdl), neg$tc - neg$hdl, NA)
neg$nhdl_ldl_gap <- ifelse(!is.na(neg$non_hdl) & !is.na(neg$ldl_ut), neg$non_hdl - neg$ldl_ut, NA)
neg$cumulative_ldl <- ifelse(!is.na(neg$ldl_ut) & !is.na(neg$age), neg$ldl_ut * neg$age, NA)
neg$gene <- NA; neg$index_effect <- 0; neg$cohort <- "SouthWales"

cols <- c("fh", "age", "sex", "hdl", "tg", "ldl_m", "ldl_ut", "trig_filter", "on_statin",
          "tendon_xanth", "corneal_arcus", "dlcn", "apob", "apob_ldl", "non_hdl",
          "nhdl_ldl_gap", "cumulative_ldl", "gene", "index_effect", "cohort")
sw <- rbind(pos[, cols], neg[, cols])
cat(sprintf("  SW: %d (FH+=%d, FH-=%d)\n", nrow(sw), sum(sw$fh == 1), sum(sw$fh == 0)))

# WALES PASS
pass <- as.data.frame(read_sav(file.path(DD, "PASS_wrong_dob.sav")))
pass$fh <- ifelse(pass$Positive1 == "Yes", 1, ifelse(pass$Positive1 == "No", 0, NA))
pass$age <- safe_num(pass$BMI_AGE); pass$sex <- ifelse(pass$Gender == "M", 1, 0)
pass$hdl <- safe_num(pass$HDL.1); pass$tg <- safe_num(pass$TRG.1); pass$ldl_m <- safe_num(pass$LDL.1)
pass$tc <- safe_num(pass$TC.1)
pass$rf <- sapply(pass$Treatment1.1, get_red)
pass$ldl_ut <- ifelse(pass$rf > 0, rev_ldl(pass$ldl_m, pass$rf), pass$ldl_m)
pass$trig_filter <- ifelse(!is.na(pass$ldl_ut) & !is.na(pass$tg) & pass$tg >= 0,
                            pmin(pass$ldl_ut / (pass$tg + 0.1), 50), NA)
pass$on_statin <- ifelse(pass$rf > 0, 1, 0)
pass$tendon_xanth <- ifelse(safe_num(pass$TendonXanthomata) == 1, 1, 0)
pass$tendon_xanth[is.na(pass$tendon_xanth)] <- 0
pass$corneal_arcus <- ifelse(safe_num(pass$CornealArcus) == 1, 1, 0)
pass$corneal_arcus[is.na(pass$corneal_arcus)] <- 0
pass$dlcn <- safe_num(pass$GenoTypingScore)
pass$apob <- NA; pass$apob_ldl <- NA
pass$non_hdl <- ifelse(!is.na(pass$tc) & !is.na(pass$hdl), pass$tc - pass$hdl, NA)
pass$nhdl_ldl_gap <- ifelse(!is.na(pass$non_hdl) & !is.na(pass$ldl_ut), pass$non_hdl - pass$ldl_ut, NA)
pass$cumulative_ldl <- ifelse(!is.na(pass$ldl_ut) & !is.na(pass$age), pass$ldl_ut * pass$age, NA)
pass$gene <- ifelse(grepl("^LDLR", pass$Mutation1), "LDLR",
              ifelse(grepl("^APOB", pass$Mutation1), "APOB",
              ifelse(grepl("^PCSK9", pass$Mutation1), "PCSK9", NA)))
pass$index_effect <- ifelse(safe_num(pass$I_Vs_R) == 1, 1, 0)
pass$index_effect[is.na(pass$index_effect)] <- 0; pass$cohort <- "Wales"
wales <- pass[!is.na(pass$fh), cols]
cat(sprintf("  Wales: %d (FH+=%d, FH-=%d)\n", nrow(wales), sum(wales$fh == 1), sum(wales$fh == 0)))

# UKB LIPID CLINIC
ukb_raw <- as.data.frame(fread(file.path(TD, "TUDOR_UKB_Features.csv"), showProgress = FALSE))
ukb_raw$fh <- safe_num(ukb_raw$is_fh_genetic); ukb_raw$fh[is.na(ukb_raw$fh)] <- 0
ukb_raw$age <- safe_num(ukb_raw$Age_at_LDL1); ukb_raw$sex <- safe_num(ukb_raw$Gender_num)
ukb_raw$hdl <- safe_num(ukb_raw$HDL.1); ukb_raw$tg <- safe_num(ukb_raw$TRG.1)
ukb_raw$ldl_m <- safe_num(ukb_raw$LDL_treated); ukb_raw$ldl_ut <- safe_num(ukb_raw$LDL_untreated)
ukb_raw$trig_filter <- safe_num(ukb_raw$Trig_Filter)
ukb_raw$on_statin <- ifelse(safe_num(ukb_raw$reduction_factor) > 0, 1, 0)
ukb_raw$on_statin[is.na(ukb_raw$on_statin)] <- 0
ukb_raw$tc <- safe_num(ukb_raw$CHOL)
ukb_raw$non_hdl <- ifelse(!is.na(ukb_raw$tc) & !is.na(ukb_raw$hdl), ukb_raw$tc - ukb_raw$hdl, NA)
ukb_raw$apob <- safe_num(ukb_raw$APOB)
ukb_raw$apob_ldl <- ifelse(!is.na(ukb_raw$apob) & !is.na(ukb_raw$ldl_ut) & ukb_raw$ldl_ut > 0,
                            ukb_raw$apob / ukb_raw$ldl_ut, NA)
ukb_raw$nhdl_ldl_gap <- ifelse(!is.na(ukb_raw$non_hdl) & !is.na(ukb_raw$ldl_ut),
                                ukb_raw$non_hdl - ukb_raw$ldl_ut, NA)
ukb_raw$cumulative_ldl <- ifelse(!is.na(ukb_raw$ldl_ut) & !is.na(ukb_raw$age),
                                  ukb_raw$ldl_ut * ukb_raw$age, NA)
# Lipid clinic filter
ukb_raw$lc <- (!is.na(ukb_raw$tc) & ukb_raw$tc > 7.9) |
              (!is.na(ukb_raw$ldl_ut) & ukb_raw$ldl_ut > 4.9) |
              (!is.na(ukb_raw$non_hdl) & ukb_raw$non_hdl > 5.9) |
              (ukb_raw$on_statin == 1)
ukb <- ukb_raw[ukb_raw$lc == TRUE, ]
ukb$dlcn <- ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut >= 8.5, 8,
            ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut >= 6.5, 5,
            ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut >= 5.0, 3,
            ifelse(!is.na(ukb$ldl_ut) & ukb$ldl_ut >= 4.0, 1, 0))))
ukb$tendon_xanth <- 0; ukb$corneal_arcus <- 0
ukb$gene <- ukb$gene; ukb$index_effect <- 0; ukb$cohort <- "UKB"
ukb_df <- ukb[, cols]
cat(sprintf("  UKB LC: %d (FH+=%d, FH-=%d)\n", nrow(ukb_df), sum(ukb_df$fh == 1), sum(ukb_df$fh == 0)))

all_data <- rbind(sw, wales, ukb_df)
cat(sprintf("  TOTAL: %d\n\n", nrow(all_data)))

# Feature distributions
cat("  nonHDL-LDL gap: FH+ =", round(mean(sw$nhdl_ldl_gap[sw$fh == 1], na.rm = TRUE), 2),
    " FH- =", round(mean(sw$nhdl_ldl_gap[sw$fh == 0], na.rm = TRUE), 2), "\n")
cat("  Cumulative LDL: FH+ =", round(mean(sw$cumulative_ldl[sw$fh == 1], na.rm = TRUE), 0),
    " FH- =", round(mean(sw$cumulative_ldl[sw$fh == 0], na.rm = TRUE), 0), "\n\n")

# ═══════════════════════════════════════════════════════════════
# PART B: LOCO-CV
# ═══════════════════════════════════════════════════════════════
cat("=== PART B: LOCO-CV ===\n")
bf <- c("ldl_ut", "trig_filter", "hdl", "tg", "age", "sex", "tendon_xanth", "corneal_arcus", "on_statin")
ef <- c(bf, "nhdl_ldl_gap", "cumulative_ldl")
af <- c(ef, "apob_ldl")
res <- list()

for (ho in c("SouthWales", "Wales", "UKB")) {
  cat(sprintf("\n-- HELD-OUT: %s --\n", ho))
  tr <- all_data[all_data$cohort != ho, ]
  te <- all_data[all_data$cohort == ho, ]

  tiers <- list(
    list(n = "Base", f = bf),
    list(n = "Enhanced", f = ef),
    list(n = "Enhanced+ApoB", f = af)
  )

  for (tier in tiers) {
    tr_c <- tr[complete.cases(tr[, c("fh", tier$f)]), ]
    te_c <- te[complete.cases(te[, c("fh", tier$f)]), ]
    if (nrow(te_c) < 50 || sum(te_c$fh == 1) < 5 || nrow(tr_c) < 50) {
      cat(sprintf("  %-16s SKIPPED\n", tier$n)); next
    }
    Xtr <- as.matrix(tr_c[, tier$f]); ytr <- tr_c$fh
    Xte <- as.matrix(te_c[, tier$f]); yte <- te_c$fh
    cv <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 0.5,
                    nfolds = 10, type.measure = "auc", standardize = TRUE)
    p <- as.numeric(predict(cv, newx = Xte, s = "lambda.min", type = "response"))
    r <- roc(yte, p, quiet = TRUE); ci <- ci.auc(r, quiet = TRUE)
    co <- coords(r, "best", best.method = "youden")
    pv <- mean(yte)
    ppv <- (co$sensitivity * pv) / (co$sensitivity * pv + (1 - co$specificity) * (1 - pv))
    npv <- (co$specificity * (1 - pv)) / ((1 - co$sensitivity) * pv + co$specificity * (1 - pv))
    cat(sprintf("  %-16s AUC=%.3f (%.3f-%.3f) Sens=%.1f%% Spec=%.1f%% [n=%d, FH+=%d]\n",
                tier$n, as.numeric(auc(r)), ci[1], ci[3],
                co$sensitivity * 100, co$specificity * 100, nrow(te_c), sum(yte == 1)))
    res[[length(res) + 1]] <- data.frame(
      HeldOut = ho, Model = tier$n,
      AUC = round(as.numeric(auc(r)), 4),
      CI_lo = round(ci[1], 4), CI_hi = round(ci[3], 4),
      Sens = round(co$sensitivity, 4), Spec = round(co$specificity, 4),
      PPV = round(ppv, 4), NPV = round(npv, 4),
      N = nrow(te_c), N_FH = sum(yte == 1)
    )
  }

  # eDLCN (raw score — no ML)
  te_b <- te[complete.cases(te[, c("fh", bf)]), ]
  dlcn_ok <- !is.na(te_b$dlcn)
  if (sum(dlcn_ok) > 50) {
    r_d <- roc(te_b$fh[dlcn_ok], te_b$dlcn[dlcn_ok], quiet = TRUE)
    cat(sprintf("  eDLCN            AUC=%.3f [n=%d]\n", as.numeric(auc(r_d)), sum(dlcn_ok)))
  }
}

# ═══════════════════════════════════════════════════════════════
# PART C: NMR (UKB nested CV)
# ═══════════════════════════════════════════════════════════════
cat("\n=== PART C: NMR ===\n")
nmr_files <- paste0(TD, "/ukb_nmr_", c("a", "b", "c", "d", "e"), ".csv")
nmr_all <- NULL
for (nf in nmr_files) {
  if (file.exists(nf)) {
    cat(sprintf("  %s...\n", basename(nf)))
    tmp <- fread(nf, showProgress = FALSE)
    if (is.null(nmr_all)) nmr_all <- tmp
    else nmr_all <- merge(nmr_all, tmp, by = "participant.eid", all = TRUE)
  }
}
if (!is.null(nmr_all)) {
  ukb_nmr <- as.data.frame(merge(as.data.table(ukb), nmr_all, by.x = "participant.eid", by.y = "participant.eid", all.x = TRUE))
  nk <- c("participant.p23407_i0", "participant.p23408_i0", "participant.p23409_i0",
           "participant.p23485_i0", "participant.p23484_i0", "participant.p23493_i0",
           "participant.p23491_i0", "participant.p23492_i0",
           "participant.p23477_i0", "participant.p23478_i0")
  nk <- nk[nk %in% names(ukb_nmr)]
  cat(sprintf("  NMR features: %d\n", length(nk)))
  anf <- c(ef, nk)
  ukb_cc <- ukb_nmr[complete.cases(ukb_nmr[, c("fh", anf)]), ]
  cat(sprintf("  Complete: %d (FH+=%d)\n", nrow(ukb_cc), sum(ukb_cc$fh == 1)))

  if (nrow(ukb_cc) > 1000 && sum(ukb_cc$fh == 1) >= 30) {
    n <- nrow(ukb_cc); folds <- sample(rep(1:5, length.out = n))
    pn <- pe <- numeric(n)
    for (k in 1:5) {
      cat(sprintf("  Fold %d/5...\n", k))
      ti <- folds != k; vi <- folds == k
      cn <- cv.glmnet(as.matrix(ukb_cc[ti, anf]), ukb_cc$fh[ti],
                      family = "binomial", alpha = 0.5, nfolds = 5,
                      type.measure = "auc", standardize = TRUE)
      pn[vi] <- as.numeric(predict(cn, newx = as.matrix(ukb_cc[vi, anf]),
                                    s = "lambda.min", type = "response"))
      ce <- cv.glmnet(as.matrix(ukb_cc[ti, ef]), ukb_cc$fh[ti],
                      family = "binomial", alpha = 0.5, nfolds = 5,
                      type.measure = "auc", standardize = TRUE)
      pe[vi] <- as.numeric(predict(ce, newx = as.matrix(ukb_cc[vi, ef]),
                                    s = "lambda.min", type = "response"))
    }
    rn <- roc(ukb_cc$fh, pn, quiet = TRUE); cin <- ci.auc(rn, quiet = TRUE)
    re <- roc(ukb_cc$fh, pe, quiet = TRUE); cie <- ci.auc(re, quiet = TRUE)
    dl <- tryCatch(roc.test(re, rn, method = "delong", paired = TRUE), error = function(e) NULL)
    cat(sprintf("\n  NMR-TUDOR:      AUC=%.3f (%.3f-%.3f)\n", as.numeric(auc(rn)), cin[1], cin[3]))
    cat(sprintf("  Enhanced-TUDOR: AUC=%.3f (%.3f-%.3f)\n", as.numeric(auc(re)), cie[1], cie[3]))
    if (!is.null(dl)) cat(sprintf("  DeLong: dAUC=%+.3f p=%s\n",
                                   as.numeric(auc(rn)) - as.numeric(auc(re)),
                                   formatC(dl$p.value, format = "e", digits = 2)))
    res[[length(res) + 1]] <- data.frame(
      HeldOut = "UKB_5fCV", Model = "NMR-TUDOR",
      AUC = round(as.numeric(auc(rn)), 4), CI_lo = round(cin[1], 4), CI_hi = round(cin[3], 4),
      Sens = NA, Spec = NA, PPV = NA, NPV = NA,
      N = nrow(ukb_cc), N_FH = sum(ukb_cc$fh == 1)
    )
  }
} else {
  cat("  NMR not found\n")
}

# ═══════════════════════════════════════════════════════════════
# PART D: SUBGROUPS
# ═══════════════════════════════════════════════════════════════
cat("\n=== PART D: Subgroups ===\n")
sgr <- list()
tr_s <- all_data[all_data$cohort %in% c("SouthWales", "UKB"), ]
tr_s <- tr_s[complete.cases(tr_s[, c("fh", bf)]), ]
cv_s <- cv.glmnet(as.matrix(tr_s[, bf]), tr_s$fh, family = "binomial",
                  alpha = 0.5, nfolds = 10, type.measure = "auc", standardize = TRUE)
te_s <- wales[complete.cases(wales[, c("fh", bf)]), ]
Xs <- as.matrix(te_s[, bf])
ps <- as.numeric(predict(cv_s, newx = Xs, s = "lambda.min", type = "response"))

sgs <- list(
  Male = te_s$sex == 1, Female = te_s$sex == 0,
  "Age<40" = te_s$age < 40, "Age40-60" = te_s$age >= 40 & te_s$age < 60,
  "Age>60" = te_s$age >= 60,
  On_statin = te_s$on_statin == 1, No_statin = te_s$on_statin == 0,
  Index = te_s$index_effect == 1, Cascade = te_s$index_effect == 0
)
for (sg in names(sgs)) {
  m <- sgs[[sg]]; ys <- te_s$fh[m]; pp <- ps[m]
  if (length(unique(ys)) < 2 || sum(ys == 1) < 5) next
  rs <- roc(ys, pp, quiet = TRUE); cs <- ci.auc(rs, quiet = TRUE)
  cat(sprintf("  %-12s AUC=%.3f (%.3f-%.3f) [n=%d, FH+=%d]\n",
              sg, as.numeric(auc(rs)), cs[1], cs[3], sum(m), sum(ys == 1)))
  sgr[[length(sgr) + 1]] <- data.frame(
    Subgroup = sg, AUC = round(as.numeric(auc(rs)), 4),
    CI_lo = round(cs[1], 4), CI_hi = round(cs[3], 4),
    N = sum(m), N_FH = sum(ys == 1)
  )
}
# Gene-specific
for (g in c("LDLR", "APOB")) {
  fg <- te_s[!is.na(te_s$gene) & te_s$gene == g & te_s$fh == 1, ]
  ng <- te_s[te_s$fh == 0, ]
  if (nrow(fg) < 5) next
  cg <- rbind(fg, ng)
  Xg <- as.matrix(cg[, bf])
  pg <- as.numeric(predict(cv_s, newx = Xg, s = "lambda.min", type = "response"))
  rg <- roc(cg$fh, pg, quiet = TRUE); cig <- ci.auc(rg, quiet = TRUE)
  cat(sprintf("  %-12s AUC=%.3f (%.3f-%.3f) [FH+=%d]\n",
              g, as.numeric(auc(rg)), cig[1], cig[3], nrow(fg)))
  sgr[[length(sgr) + 1]] <- data.frame(
    Subgroup = g, AUC = round(as.numeric(auc(rg)), 4),
    CI_lo = round(cig[1], 4), CI_hi = round(cig[3], 4),
    N = nrow(cg), N_FH = nrow(fg)
  )
}

# ═══════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════
cat("\n\n================================================================\n")
cat("  FULL RESULTS\n")
cat("================================================================\n\n")

rdf <- do.call(rbind, res)
print(rdf[, c("HeldOut", "Model", "AUC", "CI_lo", "CI_hi", "N", "N_FH")])

cat("\n=== Pooled Average ===\n")
for (m in unique(rdf$Model[rdf$HeldOut != "UKB_5fCV"])) {
  s <- rdf[rdf$Model == m & rdf$HeldOut != "UKB_5fCV", ]
  if (nrow(s) >= 2) cat(sprintf("  %-20s Mean=%.3f [%.3f, %.3f]\n",
                                  m, mean(s$AUC), min(s$AUC), max(s$AUC)))
}

if (length(sgr) > 0) {
  sdf <- do.call(rbind, sgr)
  cat("\n=== Subgroups ===\n")
  print(sdf)
}

# Save to D: drive as requested
write.csv(rdf, file.path(OD, "LOCO_CV_full_results.csv"), row.names = FALSE)
write.csv(rdf, "D:/alphafold_backup/tudor_packup/LOCO_CV_full_results.csv", row.names = FALSE)
if (length(sgr) > 0) {
  write.csv(sdf, file.path(OD, "subgroup_analyses.csv"), row.names = FALSE)
  write.csv(sdf, "D:/alphafold_backup/tudor_packup/subgroup_analyses.csv", row.names = FALSE)
}

cat(sprintf("\n  Saved to: %s/ AND D:/alphafold_backup/tudor_packup/\n", OD))
cat(sprintf("  Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("================================================================\n")

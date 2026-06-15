#!/usr/bin/env Rscript
# ===========================================================================
# ApoB/LDL-C Discordance: FH vs Non-FH in UKB
# Compares ASCVD+ vs ASCVD- within FH (pooled SW+UKB) and non-FH UKB
# Then: matched FH+ASCVD vs non-FH+ASCVD (age, BMI, TG, T2DM, statin)
# ===========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tableone)
  library(MatchIt)
})

cat("═══════════════════════════════════════════════════════════════════════\n")
cat("  ApoB/LDL-C DISCORDANCE: FH vs NON-FH COMPARISON\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# ── 1. LOAD FH DATA (UKB + South Wales) ─────────────────────────────────────
cat("Loading FH data...\n")

# UKB FH from analysis-ready file
fh_ukb <- read.csv("calon_ukb_analysis_ready.csv", stringsAsFactors = FALSE)
fh_ukb$cohort <- "UKB-FH"
fh_ukb$eid <- as.character(fh_ukb$eid)
cat(sprintf("  UKB FH: N=%d\n", nrow(fh_ukb)))

# South Wales FH from merged file
merged <- read.csv("calon2_full_merged.csv", stringsAsFactors = FALSE)
fh_sw <- merged %>% filter(cohort == "South Wales FH")
cat(sprintf("  South Wales FH: N=%d\n", nrow(fh_sw)))

# Harmonise FH UKB columns
fh_ukb_h <- fh_ukb %>%
  transmute(
    id = eid, cohort = "UKB-FH",
    age = age, sex = sex, bmi = bmi,
    ldl = ldl, hdl = hdl, tc = tc, tg = trig,
    apob = apob, inv_hdl = inv_hdl,
    apob_ldl = apob_ldl_ratio,
    log_apob_ldl = log_apob_ldl,
    diabetes = diabetes, hypertension = hypertension,
    smoking_binary = smoking_binary,
    on_statin = on_statin,
    lpa_nmol = lpa,
    ascvd = ascvd_combined,
    fh_status = 1
  )

# Harmonise South Wales columns
fh_sw_h <- fh_sw %>%
  transmute(
    id = as.character(id), cohort = "SW-FH",
    age = age, sex = sex, bmi = bmi,
    ldl = ldl, hdl = hdl, tc = tc, tg = tg,
    apob = apob, inv_hdl = inv_hdl,
    apob_ldl = apob_ldl,
    log_apob_ldl = log_apob_ldl,
    diabetes = diabetes, hypertension = hypertension,
    smoking_binary = smoking_binary,
    on_statin = on_statin,
    lpa_nmol = lpa_nmol,
    ascvd = ascvd,
    fh_status = 1
  )

fh_pooled <- bind_rows(fh_ukb_h, fh_sw_h)
cat(sprintf("  FH pooled (UKB + SW): N=%d (ASCVD+ = %d, %.1f%%)\n",
            nrow(fh_pooled), sum(fh_pooled$ascvd==1, na.rm=TRUE),
            100*mean(fh_pooled$ascvd==1, na.rm=TRUE)))

# ── 2. LOAD NON-FH UKB DATA ─────────────────────────────────────────────────
cat("\nLoading non-FH UKB data...\n")

# TUDOR file has is_fh_genetic flag + basic fields
tudor <- read.csv("TUDOR_UKB_Features (1).csv", stringsAsFactors = FALSE)
names(tudor) <- make.names(names(tudor), unique = TRUE)
cat(sprintf("  TUDOR total: N=%d\n", nrow(tudor)))

# Load batch2 lipids for ApoB
batch2 <- read.csv("calon_batch2_lipids.csv", stringsAsFactors = FALSE)
names(batch2) <- make.names(names(batch2), unique = TRUE)
cat(sprintf("  Batch2 lipids: N=%d\n", nrow(batch2)))

# Load batch1 demographics for BMI
batch1 <- read.csv("calon_batch1_demographics.csv", stringsAsFactors = FALSE)
names(batch1) <- make.names(names(batch1), unique = TRUE)

# Identify non-FH
tudor$is_fh <- as.numeric(tudor$is_fh_genetic)
nonfh <- tudor %>% filter(is_fh == 0 | is.na(is_fh))
cat(sprintf("  Non-FH in TUDOR: N=%d\n", nrow(nonfh)))

# Extract variables from TUDOR
# p31 = sex, p21022 = age, p21001_i0 = BMI, p30780 = LDL, p30760 = HDL
# p30870 = TG, p30690 = TC, p30640 = ApoB, p6150 = vascular diagnoses
nonfh_h <- nonfh %>%
  transmute(
    id = as.character(participant.eid),
    age = as.numeric(participant.p21022),
    sex = as.numeric(participant.p31),
    bmi = as.numeric(participant.p21001_i0),
    ldl = as.numeric(participant.p30780_i0),
    hdl = as.numeric(participant.p30760_i0),
    tg  = as.numeric(participant.p30870_i0),
    tc  = as.numeric(participant.p30690_i0),
    apob_tudor = as.numeric(participant.p30640_i0),
    sbp1 = as.numeric(participant.p4080_i0_a0),
    sbp2 = as.numeric(participant.p4080_i0_a1)
  )

# Merge batch2 ApoB (field p30640)
batch2_sub <- batch2 %>%
  transmute(
    id = as.character(participant.eid),
    apob_batch = as.numeric(participant.p30640_i0),
    lpa_nmol_batch = as.numeric(participant.p30790_i0)
  )
nonfh_h <- left_join(nonfh_h, batch2_sub, by = "id")

# Use TUDOR ApoB if available, else batch2
nonfh_h$apob <- ifelse(!is.na(nonfh_h$apob_tudor), nonfh_h$apob_tudor, nonfh_h$apob_batch)
nonfh_h$lpa_nmol <- nonfh_h$lpa_nmol_batch

# Compute ApoB/LDL-C ratio
nonfh_h <- nonfh_h %>%
  mutate(
    inv_hdl = ifelse(!is.na(hdl) & hdl > 0, 1/hdl, NA),
    apob_ldl = ifelse(!is.na(apob) & !is.na(ldl) & ldl > 0, apob/ldl, NA),
    log_apob_ldl = ifelse(!is.na(apob_ldl) & apob_ldl > 0, log(apob_ldl), NA)
  )

# ASCVD from self-report (p6150: 1=MI, 2=angina, 3=stroke)
# Need to go back to TUDOR for the p6150 field
sr_cols <- grep("^participant\\.p6150_i0", names(nonfh), value = TRUE)
cat(sprintf("  Self-report ASCVD columns (p6150): %d found\n", length(sr_cols)))

if (length(sr_cols) > 0) {
  sr_ascvd <- apply(nonfh[, sr_cols, drop = FALSE], 1, function(row) {
    vals <- as.numeric(row)
    any(vals %in% c(1, 2, 3), na.rm = TRUE)
  })
  nonfh_h$sr_ascvd <- as.numeric(sr_ascvd)
} else {
  cat("  WARNING: p6150 not found - using non-cancer illness codes\n")
  nonfh_h$sr_ascvd <- NA
}

# Also check non-cancer illness codes (p20002): 1065=MI, 1066=HF, 1081=stroke
illness_cols <- grep("^participant\\.p20002_i0", names(nonfh), value = TRUE)
if (length(illness_cols) > 0) {
  illness_ascvd <- apply(nonfh[, illness_cols, drop = FALSE], 1, function(row) {
    vals <- as.numeric(row)
    any(vals %in% c(1065, 1066, 1081, 1082, 1083), na.rm = TRUE)
  })
  nonfh_h$illness_ascvd <- as.numeric(illness_ascvd)
  nonfh_h$ascvd <- as.numeric(nonfh_h$sr_ascvd == 1 | nonfh_h$illness_ascvd == 1)
} else {
  nonfh_h$ascvd <- nonfh_h$sr_ascvd
}

# Diabetes from p20002 (code 1220 = diabetes, 1223 = T2DM)
if (length(illness_cols) > 0) {
  has_diabetes <- apply(nonfh[, illness_cols, drop = FALSE], 1, function(row) {
    vals <- as.numeric(row)
    any(vals %in% c(1220, 1223), na.rm = TRUE)
  })
  nonfh_h$diabetes <- as.numeric(has_diabetes)
} else {
  nonfh_h$diabetes <- NA
}

# Hypertension from p20002 (code 1065=MI already used; 1072 = essential HTN)
# Also p6150 code 4 = high blood pressure
if (length(sr_cols) > 0) {
  has_htn_sr <- apply(nonfh[, sr_cols, drop = FALSE], 1, function(row) {
    vals <- as.numeric(row)
    any(vals %in% c(4), na.rm = TRUE)
  })
} else {
  has_htn_sr <- rep(FALSE, nrow(nonfh_h))
}
if (length(illness_cols) > 0) {
  has_htn_ill <- apply(nonfh[, illness_cols, drop = FALSE], 1, function(row) {
    vals <- as.numeric(row)
    any(vals %in% c(1065, 1072), na.rm = TRUE)  # 1072 = essential HTN
  })
} else {
  has_htn_ill <- rep(FALSE, nrow(nonfh_h))
}
nonfh_h$hypertension <- as.numeric(has_htn_sr | has_htn_ill)

# Statin use from medication codes (p20003)
med_cols <- grep("^participant\\.p20003_i0", names(nonfh), value = TRUE)
if (length(med_cols) > 0) {
  # Statin codes: atorvastatin=1141146234, simvastatin=1140861958, rosuvastatin=1141192414,
  # pravastatin=1140888594, fluvastatin=1140888648
  statin_codes <- c(1141146234, 1140861958, 1141192414, 1140888594, 1140888648,
                    1141146138, 1141146188, 1141200040)
  has_statin <- apply(nonfh[, med_cols, drop = FALSE], 1, function(row) {
    vals <- as.numeric(row)
    any(vals %in% statin_codes, na.rm = TRUE)
  })
  nonfh_h$on_statin <- as.numeric(has_statin)
} else {
  nonfh_h$on_statin <- NA
}

# Smoking from p20116 (field not in TUDOR directly; check)
smoke_cols <- grep("p20116|p20107|p20111", names(nonfh), value = TRUE)
if (length(smoke_cols) > 0) {
  # p20107 = age started smoking, p20111 = pack years
  # If age started smoking is non-NA, they were ever-smokers
  nonfh_h$smoking_binary <- as.numeric(!is.na(as.numeric(nonfh$participant.p20107_i0)))
} else {
  nonfh_h$smoking_binary <- NA
}

nonfh_h$fh_status <- 0
nonfh_h$cohort <- "UKB-nonFH"

cat(sprintf("  Non-FH harmonised: N=%d\n", nrow(nonfh_h)))
cat(sprintf("  Non-FH with ApoB: %d (%.1f%%)\n",
            sum(!is.na(nonfh_h$apob)), 100*mean(!is.na(nonfh_h$apob))))
cat(sprintf("  Non-FH with ApoB/LDL: %d (%.1f%%)\n",
            sum(!is.na(nonfh_h$apob_ldl)), 100*mean(!is.na(nonfh_h$apob_ldl))))
cat(sprintf("  Non-FH ASCVD+: %d (%.1f%%)\n",
            sum(nonfh_h$ascvd==1, na.rm=TRUE), 100*mean(nonfh_h$ascvd==1, na.rm=TRUE)))

# ── 3. ANALYSIS 1: FH pooled ASCVD+ vs ASCVD- (ApoB/LDL discordance) ───────
cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("ANALYSIS 1: ApoB/LDL-C Discordance — ASCVD+ vs ASCVD- WITHIN Groups\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# FH pooled (SW + UKB) — only those with ApoB
fh_apob <- fh_pooled %>% filter(!is.na(apob_ldl))
cat(sprintf("FH pooled with ApoB/LDL: N=%d\n", nrow(fh_apob)))

fh_ascvd_pos <- fh_apob %>% filter(ascvd == 1)
fh_ascvd_neg <- fh_apob %>% filter(ascvd == 0)

cat(sprintf("  ASCVD+: N=%d, median ApoB/LDL = %.4f [%.4f-%.4f]\n",
            nrow(fh_ascvd_pos),
            median(fh_ascvd_pos$apob_ldl, na.rm=TRUE),
            quantile(fh_ascvd_pos$apob_ldl, 0.25, na.rm=TRUE),
            quantile(fh_ascvd_pos$apob_ldl, 0.75, na.rm=TRUE)))
cat(sprintf("  ASCVD-: N=%d, median ApoB/LDL = %.4f [%.4f-%.4f]\n",
            nrow(fh_ascvd_neg),
            median(fh_ascvd_neg$apob_ldl, na.rm=TRUE),
            quantile(fh_ascvd_neg$apob_ldl, 0.25, na.rm=TRUE),
            quantile(fh_ascvd_neg$apob_ldl, 0.75, na.rm=TRUE)))

wt_fh <- wilcox.test(apob_ldl ~ ascvd, data = fh_apob)
tt_fh <- t.test(log_apob_ldl ~ ascvd, data = fh_apob)
cat(sprintf("  Wilcoxon p = %.4e, t-test (log) p = %.4e\n", wt_fh$p.value, tt_fh$p.value))

# Effect sizes
fh_smd <- (mean(fh_ascvd_pos$log_apob_ldl, na.rm=TRUE) - mean(fh_ascvd_neg$log_apob_ldl, na.rm=TRUE)) /
  sqrt((var(fh_ascvd_pos$log_apob_ldl, na.rm=TRUE) + var(fh_ascvd_neg$log_apob_ldl, na.rm=TRUE)) / 2)
cat(sprintf("  SMD (log ApoB/LDL): %.3f\n", fh_smd))

# Non-FH UKB — only those with ApoB/LDL
nonfh_apob <- nonfh_h %>% filter(!is.na(apob_ldl))
cat(sprintf("\nNon-FH UKB with ApoB/LDL: N=%d\n", nrow(nonfh_apob)))

nonfh_ascvd_pos <- nonfh_apob %>% filter(ascvd == 1)
nonfh_ascvd_neg <- nonfh_apob %>% filter(ascvd == 0)

cat(sprintf("  ASCVD+: N=%d, median ApoB/LDL = %.4f [%.4f-%.4f]\n",
            nrow(nonfh_ascvd_pos),
            median(nonfh_ascvd_pos$apob_ldl, na.rm=TRUE),
            quantile(nonfh_ascvd_pos$apob_ldl, 0.25, na.rm=TRUE),
            quantile(nonfh_ascvd_pos$apob_ldl, 0.75, na.rm=TRUE)))
cat(sprintf("  ASCVD-: N=%d, median ApoB/LDL = %.4f [%.4f-%.4f]\n",
            nrow(nonfh_ascvd_neg),
            median(nonfh_ascvd_neg$apob_ldl, na.rm=TRUE),
            quantile(nonfh_ascvd_neg$apob_ldl, 0.25, na.rm=TRUE),
            quantile(nonfh_ascvd_neg$apob_ldl, 0.75, na.rm=TRUE)))

wt_nonfh <- wilcox.test(apob_ldl ~ ascvd, data = nonfh_apob)
tt_nonfh <- t.test(log_apob_ldl ~ ascvd, data = nonfh_apob)
cat(sprintf("  Wilcoxon p = %.4e, t-test (log) p = %.4e\n", wt_nonfh$p.value, tt_nonfh$p.value))

nonfh_smd <- (mean(nonfh_ascvd_pos$log_apob_ldl, na.rm=TRUE) - mean(nonfh_ascvd_neg$log_apob_ldl, na.rm=TRUE)) /
  sqrt((var(nonfh_ascvd_pos$log_apob_ldl, na.rm=TRUE) + var(nonfh_ascvd_neg$log_apob_ldl, na.rm=TRUE)) / 2)
cat(sprintf("  SMD (log ApoB/LDL): %.3f\n", nonfh_smd))

# ── 4. FH vs non-FH overall ApoB/LDL comparison ─────────────────────────────
cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("ANALYSIS 2: FH vs Non-FH — Overall ApoB/LDL-C Ratio\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

cat(sprintf("FH (pooled SW+UKB):  mean ApoB/LDL = %.4f +/- %.4f, N=%d\n",
            mean(fh_apob$apob_ldl, na.rm=TRUE), sd(fh_apob$apob_ldl, na.rm=TRUE), nrow(fh_apob)))
cat(sprintf("Non-FH (UKB):        mean ApoB/LDL = %.4f +/- %.4f, N=%d\n",
            mean(nonfh_apob$apob_ldl, na.rm=TRUE), sd(nonfh_apob$apob_ldl, na.rm=TRUE), nrow(nonfh_apob)))
fh_nonfh_test <- wilcox.test(fh_apob$apob_ldl, nonfh_apob$apob_ldl)
cat(sprintf("Wilcoxon p = %.4e\n", fh_nonfh_test$p.value))

# Log scale
cat(sprintf("\nFH:     mean log(ApoB/LDL) = %.4f +/- %.4f\n",
            mean(fh_apob$log_apob_ldl, na.rm=TRUE), sd(fh_apob$log_apob_ldl, na.rm=TRUE)))
cat(sprintf("Non-FH: mean log(ApoB/LDL) = %.4f +/- %.4f\n",
            mean(nonfh_apob$log_apob_ldl, na.rm=TRUE), sd(nonfh_apob$log_apob_ldl, na.rm=TRUE)))
tt_overall <- t.test(fh_apob$log_apob_ldl, nonfh_apob$log_apob_ldl)
cat(sprintf("t-test p = %.4e\n", tt_overall$p.value))

# ── 5. ANALYSIS 3: MATCHED FH+ASCVD vs Non-FH+ASCVD ────────────────────────
cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("ANALYSIS 3: MATCHED FH+ASCVD vs Non-FH+ASCVD\n")
cat("  Matching on: age, BMI, TG, diabetes, statin use\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# Prepare FH ASCVD+ with ApoB
fh_ascvd_match <- fh_pooled %>%
  filter(ascvd == 1, !is.na(apob_ldl), !is.na(age), !is.na(bmi),
         !is.na(tg), !is.na(diabetes), !is.na(on_statin)) %>%
  mutate(group = 1)

# Prepare non-FH ASCVD+ with ApoB
nonfh_ascvd_match <- nonfh_h %>%
  filter(ascvd == 1, !is.na(apob_ldl), !is.na(age), !is.na(bmi),
         !is.na(tg), !is.na(diabetes), !is.na(on_statin)) %>%
  mutate(group = 0)

cat(sprintf("FH ASCVD+ with complete data: N=%d\n", nrow(fh_ascvd_match)))
cat(sprintf("Non-FH ASCVD+ with complete data: N=%d\n", nrow(nonfh_ascvd_match)))

# Combine for matching
combined <- bind_rows(
  fh_ascvd_match %>% select(id, group, age, bmi, tg, diabetes, on_statin,
                             apob, ldl, apob_ldl, log_apob_ldl, hdl, inv_hdl,
                             smoking_binary, hypertension, sex),
  nonfh_ascvd_match %>% select(id, group, age, bmi, tg, diabetes, on_statin,
                                apob, ldl, apob_ldl, log_apob_ldl, hdl, inv_hdl,
                                smoking_binary, hypertension, sex)
)

cat(sprintf("Combined for matching: N=%d (FH=%d, non-FH=%d)\n",
            nrow(combined), sum(combined$group==1), sum(combined$group==0)))

# Propensity score matching 1:3 (FH:non-FH) on age, BMI, TG, diabetes, statin
if (nrow(fh_ascvd_match) >= 10 & nrow(nonfh_ascvd_match) >= 30) {
  m <- matchit(group ~ age + bmi + tg + diabetes + on_statin,
               data = combined, method = "nearest", ratio = 3,
               caliper = 0.2, std.caliper = TRUE)

  matched <- match.data(m)
  cat(sprintf("\nMatched dataset: N=%d (FH=%d, non-FH=%d)\n",
              nrow(matched), sum(matched$group==1), sum(matched$group==0)))

  # Check balance
  cat("\nBalance after matching:\n")
  match_vars <- c("age", "bmi", "tg", "diabetes", "on_statin")
  for (v in match_vars) {
    g1 <- matched[[v]][matched$group == 1]
    g0 <- matched[[v]][matched$group == 0]
    smd_v <- (mean(g1, na.rm=TRUE) - mean(g0, na.rm=TRUE)) /
      sqrt((var(g1, na.rm=TRUE) + var(g0, na.rm=TRUE)) / 2)
    cat(sprintf("  %15s: FH=%.2f, non-FH=%.2f, SMD=%.3f\n",
                v, mean(g1, na.rm=TRUE), mean(g0, na.rm=TRUE), abs(smd_v)))
  }

  # Compare ApoB/LDL-C in matched groups
  cat("\n--- ApoB/LDL-C in MATCHED ASCVD+ patients ---\n")
  fh_m <- matched %>% filter(group == 1)
  nonfh_m <- matched %>% filter(group == 0)

  cat(sprintf("FH+ASCVD:     median ApoB/LDL = %.4f [%.4f-%.4f], mean = %.4f +/- %.4f\n",
              median(fh_m$apob_ldl, na.rm=TRUE),
              quantile(fh_m$apob_ldl, 0.25, na.rm=TRUE),
              quantile(fh_m$apob_ldl, 0.75, na.rm=TRUE),
              mean(fh_m$apob_ldl, na.rm=TRUE),
              sd(fh_m$apob_ldl, na.rm=TRUE)))
  cat(sprintf("Non-FH+ASCVD: median ApoB/LDL = %.4f [%.4f-%.4f], mean = %.4f +/- %.4f\n",
              median(nonfh_m$apob_ldl, na.rm=TRUE),
              quantile(nonfh_m$apob_ldl, 0.25, na.rm=TRUE),
              quantile(nonfh_m$apob_ldl, 0.75, na.rm=TRUE),
              mean(nonfh_m$apob_ldl, na.rm=TRUE),
              sd(nonfh_m$apob_ldl, na.rm=TRUE)))

  # Log scale
  cat(sprintf("\nFH+ASCVD:     mean log(ApoB/LDL) = %.4f +/- %.4f\n",
              mean(fh_m$log_apob_ldl, na.rm=TRUE), sd(fh_m$log_apob_ldl, na.rm=TRUE)))
  cat(sprintf("Non-FH+ASCVD: mean log(ApoB/LDL) = %.4f +/- %.4f\n",
              mean(nonfh_m$log_apob_ldl, na.rm=TRUE), sd(nonfh_m$log_apob_ldl, na.rm=TRUE)))

  wt_matched <- wilcox.test(apob_ldl ~ group, data = matched)
  tt_matched <- t.test(log_apob_ldl ~ group, data = matched)
  cat(sprintf("Wilcoxon p = %.4e, t-test (log) p = %.4e\n",
              wt_matched$p.value, tt_matched$p.value))

  match_smd <- (mean(fh_m$log_apob_ldl, na.rm=TRUE) - mean(nonfh_m$log_apob_ldl, na.rm=TRUE)) /
    sqrt((var(fh_m$log_apob_ldl, na.rm=TRUE) + var(nonfh_m$log_apob_ldl, na.rm=TRUE)) / 2)
  cat(sprintf("SMD (log ApoB/LDL): %.3f\n", match_smd))

  # Also compare ApoB and LDL separately
  cat("\n--- ApoB (g/L) in matched groups ---\n")
  cat(sprintf("FH+ASCVD:     mean ApoB = %.3f +/- %.3f\n",
              mean(fh_m$apob, na.rm=TRUE), sd(fh_m$apob, na.rm=TRUE)))
  cat(sprintf("Non-FH+ASCVD: mean ApoB = %.3f +/- %.3f\n",
              mean(nonfh_m$apob, na.rm=TRUE), sd(nonfh_m$apob, na.rm=TRUE)))
  tt_apob <- t.test(apob ~ group, data = matched)
  cat(sprintf("t-test p = %.4e\n", tt_apob$p.value))

  cat("\n--- LDL-C (mmol/L) in matched groups ---\n")
  cat(sprintf("FH+ASCVD:     mean LDL = %.3f +/- %.3f\n",
              mean(fh_m$ldl, na.rm=TRUE), sd(fh_m$ldl, na.rm=TRUE)))
  cat(sprintf("Non-FH+ASCVD: mean LDL = %.3f +/- %.3f\n",
              mean(nonfh_m$ldl, na.rm=TRUE), sd(nonfh_m$ldl, na.rm=TRUE)))
  tt_ldl <- t.test(ldl ~ group, data = matched)
  cat(sprintf("t-test p = %.4e\n", tt_ldl$p.value))

  # Binary threshold >=0.31
  cat("\n--- Binary ApoB/LDL >= 0.31 threshold ---\n")
  fh_m$ratio_high <- as.numeric(fh_m$apob_ldl >= 0.31)
  nonfh_m$ratio_high <- as.numeric(nonfh_m$apob_ldl >= 0.31)
  cat(sprintf("FH+ASCVD:     %.1f%% >= 0.31\n", 100*mean(fh_m$ratio_high, na.rm=TRUE)))
  cat(sprintf("Non-FH+ASCVD: %.1f%% >= 0.31\n", 100*mean(nonfh_m$ratio_high, na.rm=TRUE)))
  ft <- fisher.test(table(matched$group, matched$apob_ldl >= 0.31))
  cat(sprintf("Fisher's exact p = %.4e, OR = %.2f\n", ft$p.value, ft$estimate))

  # Full tableone for matched groups
  cat("\n--- Full TableOne for matched groups ---\n")
  matched$fh_label <- ifelse(matched$group == 1, "FH+ASCVD", "Non-FH+ASCVD")
  tab_vars <- c("age", "sex", "bmi", "ldl", "hdl", "tg", "apob",
                "apob_ldl", "log_apob_ldl", "inv_hdl",
                "diabetes", "hypertension", "smoking_binary", "on_statin")
  cat_vars_m <- c("sex", "diabetes", "hypertension", "smoking_binary", "on_statin")
  nonnorm_m <- c("apob_ldl", "tg")

  tab_matched <- CreateTableOne(vars = tab_vars, strata = "fh_label",
                                 data = matched, factorVars = cat_vars_m)
  print(tab_matched, nonnormal = nonnorm_m, smd = TRUE, test = TRUE)

  # Save matched results
  write.csv(print(tab_matched, nonnormal = nonnorm_m, smd = TRUE, test = TRUE, printToggle = FALSE),
            "output/final/Table_FH_vs_NonFH_ASCVD_matched.csv")
  cat("\nSaved: output/final/Table_FH_vs_NonFH_ASCVD_matched.csv\n")

} else {
  cat("INSUFFICIENT DATA for matching\n")
}

# ── 6. UNMATCHED COMPARISON TABLE (for reference) ───────────────────────────
cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("ANALYSIS 4: Unmatched FH vs Non-FH — ASCVD+ only\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# Combine ASCVD+ from both groups (unmatched)
fh_ascvd_tab <- fh_pooled %>%
  filter(ascvd == 1, !is.na(apob_ldl)) %>%
  mutate(fh_label = "FH+ASCVD")

nonfh_ascvd_tab <- nonfh_h %>%
  filter(ascvd == 1, !is.na(apob_ldl)) %>%
  mutate(fh_label = "Non-FH+ASCVD")

unmatched <- bind_rows(
  fh_ascvd_tab %>% select(fh_label, age, sex, bmi, ldl, hdl, tg, apob,
                           apob_ldl, log_apob_ldl, inv_hdl,
                           diabetes, hypertension, smoking_binary, on_statin),
  nonfh_ascvd_tab %>% select(fh_label, age, sex, bmi, ldl, hdl, tg, apob,
                              apob_ldl, log_apob_ldl, inv_hdl,
                              diabetes, hypertension, smoking_binary, on_statin)
)

tab_vars <- c("age", "sex", "bmi", "ldl", "hdl", "tg", "apob",
              "apob_ldl", "log_apob_ldl", "inv_hdl",
              "diabetes", "hypertension", "smoking_binary", "on_statin")
cat_vars_u <- c("sex", "diabetes", "hypertension", "smoking_binary", "on_statin")
nonnorm_u <- c("apob_ldl", "tg")

tab_unmatched <- CreateTableOne(vars = tab_vars, strata = "fh_label",
                                 data = unmatched, factorVars = cat_vars_u)
print(tab_unmatched, nonnormal = nonnorm_u, smd = TRUE, test = TRUE)

write.csv(print(tab_unmatched, nonnormal = nonnorm_u, smd = TRUE, test = TRUE, printToggle = FALSE),
          "output/final/Table_FH_vs_NonFH_ASCVD_unmatched.csv")
cat("\nSaved: output/final/Table_FH_vs_NonFH_ASCVD_unmatched.csv\n")

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("  DONE\n")
cat("═══════════════════════════════════════════════════════════════════════\n")

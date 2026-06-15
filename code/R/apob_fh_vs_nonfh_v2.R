#!/usr/bin/env Rscript
# ===========================================================================
# ApoB/LDL-C Discordance: FH vs Non-FH in UKB
# ===========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tableone)
})

dir.create("output/final", showWarnings = FALSE, recursive = TRUE)

cat("===================================================================\n")
cat("  ApoB/LDL-C DISCORDANCE: FH vs NON-FH COMPARISON\n")
cat("===================================================================\n\n")

# ── 1. LOAD FH DATA ──────────────────────────────────────────────────────

cat("Loading FH UKB data...\n")
fh_ukb <- read.csv("calon_ukb_analysis_ready.csv", stringsAsFactors = FALSE)
fh_ukb <- fh_ukb %>%
  filter(is.finite(apob_ldl_ratio) & apob_ldl_ratio > 0) %>%
  transmute(
    id = as.character(eid), cohort = "UKB-FH",
    age = age, sex = sex, bmi = bmi,
    ldl = ldl, hdl = hdl, tc = tc, tg = trig,
    apob = apob, apob_ldl = apob_ldl_ratio,
    log_apob_ldl = log_apob_ldl,
    diabetes = diabetes, hypertension = hypertension,
    smoking = smoking_binary, on_statin = on_statin,
    ascvd = ascvd_combined, fh_status = 1L
  )
cat(sprintf("  UKB FH with ApoB/LDL: N=%d (ASCVD+ = %d)\n", nrow(fh_ukb), sum(fh_ukb$ascvd==1)))

cat("Loading FH South Wales data...\n")
sw_raw <- read.csv("DRAGON_3.csv", stringsAsFactors = FALSE)
# Clean BMI to numeric
sw_raw$BMI_clean <- suppressWarnings(as.numeric(sw_raw$BMI))
fh_sw <- sw_raw %>%
  filter(is.finite(APOB_LDL) & APOB_LDL > 0) %>%
  transmute(
    id = as.character(DatabaseNumber), cohort = "SW-FH",
    age = age_at_event_or_censoring,
    sex = ifelse(Gender == "M", 1, 0),
    bmi = BMI_clean,
    ldl = LastLDL, hdl = LastHDL, tc = LastTC, tg = LastTrigs,
    apob = ApoB, apob_ldl = APOB_LDL,
    log_apob_ldl = log(APOB_LDL),
    diabetes = ifelse(Diabetes_binary == 1, 1, 0),
    hypertension = ifelse(onBPtreat == 1 | BP == 1, 1, 0),
    smoking = ifelse(Smoking_binary == 1, 1, 0),
    on_statin = ifelse(Statin != "" & Statin != "NO" & Statin != " ", 1, 0),
    ascvd = ifelse(ASCVD_combined == 1, 1, 0),
    fh_status = 1L
  )
cat(sprintf("  SW FH with ApoB/LDL: N=%d (ASCVD+ = %d)\n", nrow(fh_sw), sum(fh_sw$ascvd==1, na.rm=TRUE)))

fh_pooled <- bind_rows(fh_ukb, fh_sw)
cat(sprintf("  FH pooled: N=%d (ASCVD+ = %d, %.1f%%)\n\n",
            nrow(fh_pooled), sum(fh_pooled$ascvd==1, na.rm=TRUE),
            100*mean(fh_pooled$ascvd==1, na.rm=TRUE)))

# ── 2. LOAD NON-FH UKB DATA ──────────────────────────────────────────────

cat("Loading non-FH UKB data (this may take a minute)...\n")
tudor <- read.csv("TUDOR_UKB_Features (1).csv", stringsAsFactors = FALSE)
names(tudor) <- make.names(names(tudor), unique = TRUE)
cat(sprintf("  TUDOR total: N=%d\n", nrow(tudor)))

# Load batch2 lipids for ApoB (p30640)
batch2 <- read.csv("calon_batch2_lipids.csv", stringsAsFactors = FALSE)
cat(sprintf("  Batch2 lipids: N=%d\n", nrow(batch2)))

# Non-FH: those NOT in the FH cohort
fh_eids <- as.character(fh_ukb$id)
tudor$eid_str <- as.character(tudor$eid)
nonfh <- tudor %>% filter(!(eid_str %in% fh_eids))
cat(sprintf("  Non-FH after excluding FH eids: N=%d\n", nrow(nonfh)))

# Merge ApoB from batch2
batch2$eid_str <- as.character(batch2$eid)
nonfh <- nonfh %>% left_join(batch2 %>% select(eid_str, p30640, p30780, p30690, p30760, p30870), by = "eid_str")
# p30640 = ApoB, p30780 = LDL, p30690 = total chol, p30760 = HDL, p30870 = TG

# Compute ApoB/LDL ratio
nonfh$apob_g <- nonfh$p30640 / 100  # mg/dL -> g/L
nonfh$ldl_mmol <- nonfh$p30780
nonfh$apob_ldl <- ifelse(is.finite(nonfh$apob_g) & is.finite(nonfh$ldl_mmol) & nonfh$ldl_mmol > 0,
                          nonfh$apob_g / nonfh$ldl_mmol, NA)

# Define ASCVD from p6150 (self-report: 1=MI, 2=angina, 3=stroke) + p20002
# p6150 columns
p6150_cols <- grep("^p6150", names(nonfh), value = TRUE)
nonfh$sr_mi <- apply(nonfh[p6150_cols], 1, function(x) any(x == 1, na.rm = TRUE))
nonfh$sr_angina <- apply(nonfh[p6150_cols], 1, function(x) any(x == 2, na.rm = TRUE))
nonfh$sr_stroke <- apply(nonfh[p6150_cols], 1, function(x) any(x == 3, na.rm = TRUE))

# p20002 illness codes (1065=MI, 1066=HF, 1081=stroke)
p20002_cols <- grep("^p20002", names(nonfh), value = TRUE)
nonfh$ill_mi <- apply(nonfh[p20002_cols], 1, function(x) any(x == 1065, na.rm = TRUE))
nonfh$ill_hf <- apply(nonfh[p20002_cols], 1, function(x) any(x == 1066, na.rm = TRUE))
nonfh$ill_stroke <- apply(nonfh[p20002_cols], 1, function(x) any(x == 1081, na.rm = TRUE))

nonfh$ascvd <- ifelse(nonfh$sr_mi | nonfh$sr_angina | nonfh$sr_stroke |
                       nonfh$ill_mi | nonfh$ill_hf | nonfh$ill_stroke, 1, 0)

# Diabetes from p20002 (1220=T2DM, 1223=T1DM)
nonfh$diabetes <- apply(nonfh[p20002_cols], 1, function(x) any(x %in% c(1220, 1223), na.rm = TRUE))
nonfh$diabetes <- as.integer(nonfh$diabetes)

# Statin from p20003 medication codes
statin_codes <- c(1141146234, 1140861958, 1141192414, 1140888594, 1140888648, 1141146138)
p20003_cols <- grep("^p20003", names(nonfh), value = TRUE)
if (length(p20003_cols) > 0) {
  nonfh$on_statin <- apply(nonfh[p20003_cols], 1, function(x) any(x %in% statin_codes, na.rm = TRUE))
  nonfh$on_statin <- as.integer(nonfh$on_statin)
} else {
  nonfh$on_statin <- NA_integer_
}

# Smoking from p20116 (current smoking status) or p20107
p20116_cols <- grep("^p20116", names(nonfh), value = TRUE)
if (length(p20116_cols) > 0) {
  nonfh$smoking <- ifelse(nonfh[[p20116_cols[1]]] %in% c(1, 2), 1, 0)
} else {
  p20107_cols <- grep("^p20107", names(nonfh), value = TRUE)
  if (length(p20107_cols) > 0) {
    nonfh$smoking <- ifelse(nonfh[[p20107_cols[1]]] %in% c(1, 2), 1, 0)
  } else {
    nonfh$smoking <- NA_integer_
  }
}

# Age, sex, BMI from TUDOR columns
age_col <- grep("^p21003", names(nonfh), value = TRUE)[1]  # age at recruitment
sex_col <- grep("^p31$|^p31\\.", names(nonfh), value = TRUE)[1]  # sex
bmi_col <- grep("^p21001", names(nonfh), value = TRUE)[1]  # BMI

cat(sprintf("  Age col: %s, Sex col: %s, BMI col: %s\n", age_col, sex_col, bmi_col))

# Build clean non-FH dataset
nonfh_clean <- nonfh %>%
  filter(is.finite(apob_ldl) & apob_ldl > 0) %>%
  transmute(
    id = eid_str, cohort = "UKB-nonFH",
    age = if (!is.null(age_col) && age_col %in% names(.)) .[[age_col]] else NA_real_,
    sex = if (!is.null(sex_col) && sex_col %in% names(.)) .[[sex_col]] else NA_real_,
    bmi = if (!is.null(bmi_col) && bmi_col %in% names(.)) .[[bmi_col]] else NA_real_,
    ldl = ldl_mmol, hdl = p30760, tc = p30690, tg = p30870,
    apob = apob_g, apob_ldl = apob_ldl,
    log_apob_ldl = log(apob_ldl),
    diabetes = diabetes, hypertension = NA_integer_,
    smoking = smoking, on_statin = on_statin,
    ascvd = ascvd, fh_status = 0L
  )
cat(sprintf("  Non-FH UKB with ApoB/LDL: N=%d (ASCVD+ = %d, %.1f%%)\n\n",
            nrow(nonfh_clean), sum(nonfh_clean$ascvd==1, na.rm=TRUE),
            100*mean(nonfh_clean$ascvd==1, na.rm=TRUE)))


# ══════════════════════════════════════════════════════════════════════
# ANALYSIS 1: ASCVD+ vs ASCVD- WITHIN FH (pooled)
# ══════════════════════════════════════════════════════════════════════
cat("===================================================================\n")
cat("  ANALYSIS 1: ApoB/LDL-C in ASCVD+ vs ASCVD- WITHIN FH (pooled)\n")
cat("===================================================================\n")

fh_ascvd_pos <- fh_pooled %>% filter(ascvd == 1, is.finite(apob_ldl))
fh_ascvd_neg <- fh_pooled %>% filter(ascvd == 0, is.finite(apob_ldl))

cat(sprintf("  FH ASCVD+: N=%d, mean ApoB/LDL=%.4f, median=%.4f\n",
            nrow(fh_ascvd_pos), mean(fh_ascvd_pos$apob_ldl), median(fh_ascvd_pos$apob_ldl)))
cat(sprintf("  FH ASCVD-: N=%d, mean ApoB/LDL=%.4f, median=%.4f\n",
            nrow(fh_ascvd_neg), mean(fh_ascvd_neg$apob_ldl), median(fh_ascvd_neg$apob_ldl)))

wt <- wilcox.test(fh_ascvd_pos$apob_ldl, fh_ascvd_neg$apob_ldl)
tt <- t.test(fh_ascvd_pos$apob_ldl, fh_ascvd_neg$apob_ldl)
cat(sprintf("  Wilcoxon p = %.6f\n", wt$p.value))
cat(sprintf("  t-test p = %.6f\n", tt$p.value))

# SMD
smd_fh <- (mean(fh_ascvd_pos$apob_ldl) - mean(fh_ascvd_neg$apob_ldl)) /
  sqrt((sd(fh_ascvd_pos$apob_ldl)^2 + sd(fh_ascvd_neg$apob_ldl)^2) / 2)
cat(sprintf("  SMD = %.3f\n", smd_fh))

# Binary threshold >=0.31
fh_pooled_apob <- fh_pooled %>% filter(is.finite(apob_ldl))
fh_pooled_apob$apob_high <- ifelse(fh_pooled_apob$apob_ldl >= 0.31, 1, 0)
tab_fh <- table(fh_pooled_apob$apob_high, fh_pooled_apob$ascvd)
cat("\n  Binary >=0.31 g/mmol:\n")
print(tab_fh)
ft_fh <- fisher.test(tab_fh)
cat(sprintf("  Fisher's exact OR = %.2f (%.2f-%.2f), p = %.6f\n",
            ft_fh$estimate, ft_fh$conf.int[1], ft_fh$conf.int[2], ft_fh$p.value))

# ASCVD rate by threshold
rate_high <- sum(fh_pooled_apob$ascvd[fh_pooled_apob$apob_high==1]==1) / sum(fh_pooled_apob$apob_high==1)
rate_low <- sum(fh_pooled_apob$ascvd[fh_pooled_apob$apob_high==0]==1) / sum(fh_pooled_apob$apob_high==0)
cat(sprintf("  ASCVD rate if >=0.31: %.1f%% vs <0.31: %.1f%%\n\n", 100*rate_high, 100*rate_low))


# ══════════════════════════════════════════════════════════════════════
# ANALYSIS 2: ASCVD+ vs ASCVD- WITHIN NON-FH
# ══════════════════════════════════════════════════════════════════════
cat("===================================================================\n")
cat("  ANALYSIS 2: ApoB/LDL-C in ASCVD+ vs ASCVD- WITHIN NON-FH\n")
cat("===================================================================\n")

nf_ascvd_pos <- nonfh_clean %>% filter(ascvd == 1, is.finite(apob_ldl))
nf_ascvd_neg <- nonfh_clean %>% filter(ascvd == 0, is.finite(apob_ldl))

cat(sprintf("  Non-FH ASCVD+: N=%d, mean ApoB/LDL=%.4f, median=%.4f\n",
            nrow(nf_ascvd_pos), mean(nf_ascvd_pos$apob_ldl), median(nf_ascvd_pos$apob_ldl)))
cat(sprintf("  Non-FH ASCVD-: N=%d, mean ApoB/LDL=%.4f, median=%.4f\n",
            nrow(nf_ascvd_neg), mean(nf_ascvd_neg$apob_ldl), median(nf_ascvd_neg$apob_ldl)))

wt2 <- wilcox.test(nf_ascvd_pos$apob_ldl, nf_ascvd_neg$apob_ldl)
tt2 <- t.test(nf_ascvd_pos$apob_ldl, nf_ascvd_neg$apob_ldl)
cat(sprintf("  Wilcoxon p = %.6f\n", wt2$p.value))
cat(sprintf("  t-test p = %.6f\n", tt2$p.value))

smd_nf <- (mean(nf_ascvd_pos$apob_ldl) - mean(nf_ascvd_neg$apob_ldl)) /
  sqrt((sd(nf_ascvd_pos$apob_ldl)^2 + sd(nf_ascvd_neg$apob_ldl)^2) / 2)
cat(sprintf("  SMD = %.3f\n", smd_nf))

nf_apob <- nonfh_clean %>% filter(is.finite(apob_ldl))
nf_apob$apob_high <- ifelse(nf_apob$apob_ldl >= 0.31, 1, 0)
tab_nf <- table(nf_apob$apob_high, nf_apob$ascvd)
cat("\n  Binary >=0.31 g/mmol:\n")
print(tab_nf)
ft_nf <- fisher.test(tab_nf)
cat(sprintf("  Fisher's exact OR = %.2f (%.2f-%.2f), p = %.6f\n",
            ft_nf$estimate, ft_nf$conf.int[1], ft_nf$conf.int[2], ft_nf$p.value))


# ══════════════════════════════════════════════════════════════════════
# ANALYSIS 3: FH vs NON-FH overall comparison
# ══════════════════════════════════════════════════════════════════════
cat("\n===================================================================\n")
cat("  ANALYSIS 3: FH vs NON-FH overall ApoB/LDL-C\n")
cat("===================================================================\n")

fh_apob <- fh_pooled %>% filter(is.finite(apob_ldl))
nf_apob_all <- nonfh_clean %>% filter(is.finite(apob_ldl))

cat(sprintf("  FH: N=%d, mean=%.4f, median=%.4f, SD=%.4f\n",
            nrow(fh_apob), mean(fh_apob$apob_ldl), median(fh_apob$apob_ldl), sd(fh_apob$apob_ldl)))
cat(sprintf("  Non-FH: N=%d, mean=%.4f, median=%.4f, SD=%.4f\n",
            nrow(nf_apob_all), mean(nf_apob_all$apob_ldl), median(nf_apob_all$apob_ldl), sd(nf_apob_all$apob_ldl)))

wt3 <- wilcox.test(fh_apob$apob_ldl, nf_apob_all$apob_ldl)
tt3 <- t.test(fh_apob$apob_ldl, nf_apob_all$apob_ldl)
cat(sprintf("  Wilcoxon p = %g\n", wt3$p.value))
cat(sprintf("  t-test p = %g, diff = %.4f\n", tt3$p.value, diff(tt3$estimate)))

smd_overall <- (mean(fh_apob$apob_ldl) - mean(nf_apob_all$apob_ldl)) /
  sqrt((sd(fh_apob$apob_ldl)^2 + sd(nf_apob_all$apob_ldl)^2) / 2)
cat(sprintf("  SMD = %.3f\n", smd_overall))


# ══════════════════════════════════════════════════════════════════════
# ANALYSIS 4: MATCHED FH+ASCVD vs NON-FH+ASCVD
# ══════════════════════════════════════════════════════════════════════
cat("\n===================================================================\n")
cat("  ANALYSIS 4: PROPENSITY-MATCHED FH+ASCVD vs NON-FH+ASCVD\n")
cat("===================================================================\n")

# Check if MatchIt is available
if (requireNamespace("MatchIt", quietly = TRUE)) {
  library(MatchIt)

  # Combine FH ASCVD+ and non-FH ASCVD+
  fh_ascvd <- fh_pooled %>%
    filter(ascvd == 1, is.finite(apob_ldl), is.finite(bmi), is.finite(tg), is.finite(age)) %>%
    mutate(fh_status = 1L, diabetes = ifelse(is.na(diabetes), 0, diabetes),
           on_statin = ifelse(is.na(on_statin), 0, on_statin))

  nf_ascvd <- nonfh_clean %>%
    filter(ascvd == 1, is.finite(apob_ldl), is.finite(bmi), is.finite(tg), is.finite(age)) %>%
    mutate(fh_status = 0L, diabetes = ifelse(is.na(diabetes), 0, diabetes),
           on_statin = ifelse(is.na(on_statin), 0, on_statin))

  cat(sprintf("  FH ASCVD+ (complete cases): N=%d\n", nrow(fh_ascvd)))
  cat(sprintf("  Non-FH ASCVD+ (complete cases): N=%d\n", nrow(nf_ascvd)))

  combined <- bind_rows(fh_ascvd, nf_ascvd)
  combined$fh_factor <- factor(combined$fh_status, levels = c(0, 1), labels = c("NonFH", "FH"))

  # Propensity score matching: 1:3 nearest neighbor
  m <- matchit(fh_status ~ age + bmi + tg + diabetes + on_statin,
               data = combined, method = "nearest", ratio = 3,
               caliper = 0.2)

  cat("\n  Match summary:\n")
  print(summary(m))

  matched <- match.data(m)
  cat(sprintf("\n  Matched dataset: N=%d (FH=%d, NonFH=%d)\n",
              nrow(matched), sum(matched$fh_status==1), sum(matched$fh_status==0)))

  # Compare ApoB/LDL in matched groups
  m_fh <- matched %>% filter(fh_status == 1)
  m_nf <- matched %>% filter(fh_status == 0)

  cat(sprintf("\n  Matched FH ASCVD+: mean ApoB/LDL = %.4f (SD %.4f), median = %.4f\n",
              mean(m_fh$apob_ldl), sd(m_fh$apob_ldl), median(m_fh$apob_ldl)))
  cat(sprintf("  Matched NonFH ASCVD+: mean ApoB/LDL = %.4f (SD %.4f), median = %.4f\n",
              mean(m_nf$apob_ldl), sd(m_nf$apob_ldl), median(m_nf$apob_ldl)))

  wt4 <- wilcox.test(m_fh$apob_ldl, m_nf$apob_ldl)
  tt4 <- t.test(m_fh$apob_ldl, m_nf$apob_ldl)
  cat(sprintf("  Wilcoxon p = %g\n", wt4$p.value))
  cat(sprintf("  t-test p = %g\n", tt4$p.value))

  smd_matched <- (mean(m_fh$apob_ldl) - mean(m_nf$apob_ldl)) /
    sqrt((sd(m_fh$apob_ldl)^2 + sd(m_nf$apob_ldl)^2) / 2)
  cat(sprintf("  SMD = %.3f\n", smd_matched))

  # Binary threshold
  m_fh_high <- sum(m_fh$apob_ldl >= 0.31) / nrow(m_fh)
  m_nf_high <- sum(m_nf$apob_ldl >= 0.31) / nrow(m_nf)
  cat(sprintf("\n  ApoB/LDL >=0.31: FH %.1f%% vs NonFH %.1f%%\n", 100*m_fh_high, 100*m_nf_high))

  tab_matched <- table(
    factor(ifelse(matched$apob_ldl >= 0.31, "High", "Low"), levels=c("Low","High")),
    matched$fh_factor
  )
  print(tab_matched)
  ft_matched <- fisher.test(tab_matched)
  cat(sprintf("  Fisher's exact OR = %.2f (%.2f-%.2f), p = %g\n",
              ft_matched$estimate, ft_matched$conf.int[1], ft_matched$conf.int[2], ft_matched$p.value))

  # Create comparison TableOne
  vars <- c("age", "sex", "bmi", "ldl", "hdl", "tc", "tg", "apob", "apob_ldl", "diabetes", "on_statin", "smoking")
  catvars <- c("sex", "diabetes", "on_statin", "smoking")
  tab_matched_tbl <- CreateTableOne(vars = vars, strata = "fh_factor", data = matched,
                                     factorVars = catvars, test = TRUE, smd = TRUE)
  cat("\n  Matched TableOne:\n")
  print(tab_matched_tbl, smd = TRUE, test = TRUE, showAllLevels = FALSE)

  # Save matched table
  matched_out <- print(tab_matched_tbl, smd = TRUE, test = TRUE, showAllLevels = FALSE, printToggle = FALSE)
  write.csv(matched_out, "output/final/Table_FH_vs_NonFH_ASCVD_matched.csv")
  cat("\n  Saved: output/final/Table_FH_vs_NonFH_ASCVD_matched.csv\n")

} else {
  cat("  MatchIt package not available. Skipping matched analysis.\n")
  cat("  Install with: install.packages('MatchIt')\n")
}


# ══════════════════════════════════════════════════════════════════════
# ANALYSIS 5: UNMATCHED COMPARISON TABLE
# ══════════════════════════════════════════════════════════════════════
cat("\n===================================================================\n")
cat("  ANALYSIS 5: UNMATCHED FH vs NON-FH (all with ApoB/LDL)\n")
cat("===================================================================\n")

all_combined <- bind_rows(
  fh_pooled %>% filter(is.finite(apob_ldl)) %>% mutate(group = "FH"),
  nonfh_clean %>% filter(is.finite(apob_ldl)) %>% mutate(group = "NonFH")
)

vars_um <- c("age", "sex", "bmi", "ldl", "hdl", "tc", "tg", "apob", "apob_ldl", "diabetes", "on_statin", "smoking", "ascvd")
catvars_um <- c("sex", "diabetes", "on_statin", "smoking", "ascvd")
tab_um <- CreateTableOne(vars = vars_um, strata = "group", data = all_combined,
                          factorVars = catvars_um, test = TRUE, smd = TRUE)
cat("\n  Unmatched TableOne:\n")
print(tab_um, smd = TRUE, test = TRUE, showAllLevels = FALSE)

um_out <- print(tab_um, smd = TRUE, test = TRUE, showAllLevels = FALSE, printToggle = FALSE)
write.csv(um_out, "output/final/Table_FH_vs_NonFH_unmatched.csv")
cat("\n  Saved: output/final/Table_FH_vs_NonFH_unmatched.csv\n")

cat("\n===================================================================\n")
cat("  DONE\n")
cat("===================================================================\n")

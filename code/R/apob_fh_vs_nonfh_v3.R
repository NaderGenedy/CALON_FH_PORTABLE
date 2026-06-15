#!/usr/bin/env Rscript
# ===========================================================================
# ApoB/LDL-C Discordance: FH vs Non-FH
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
sw_raw$BMI_clean <- suppressWarnings(as.numeric(sw_raw$BMI))
fh_sw <- sw_raw %>%
  filter(is.finite(APOB_LDL) & APOB_LDL > 0) %>%
  transmute(
    id = as.character(DatabaseNumber), cohort = "SW-FH",
    age = age_at_event_or_censoring,
    sex = ifelse(Gender == "M", 1, 0),
    bmi = BMI_clean,
    ldl = suppressWarnings(as.numeric(LastLDL)),
    hdl = LastHDL, tc = LastTC, tg = LastTrigs,
    apob = ApoB, apob_ldl = APOB_LDL,
    log_apob_ldl = log(APOB_LDL),
    diabetes = ifelse(Diabetes_binary == 1, 1, 0),
    hypertension = ifelse(onBPtreat == 1, 1, 0),
    smoking = ifelse(Smoking_binary == 1, 1, 0),
    on_statin = ifelse(Statin %in% c("", "NO", " "), 0, 1),
    ascvd = ifelse(ASCVD_combined == 1, 1, 0),
    fh_status = 1L
  )
cat(sprintf("  SW FH with ApoB/LDL: N=%d (ASCVD+ = %d)\n", nrow(fh_sw), sum(fh_sw$ascvd==1, na.rm=TRUE)))

fh_pooled <- bind_rows(fh_ukb, fh_sw)
cat(sprintf("  FH pooled: N=%d (ASCVD+ = %d, %.1f%%)\n\n",
            nrow(fh_pooled), sum(fh_pooled$ascvd==1, na.rm=TRUE),
            100*mean(fh_pooled$ascvd==1, na.rm=TRUE)))

# ── 2. LOAD NON-FH UKB DATA ──────────────────────────────────────────────

cat("Loading non-FH UKB data...\n")
tudor <- read.csv("TUDOR_UKB_Features (1).csv", stringsAsFactors = FALSE)
cat(sprintf("  TUDOR total: N=%d\n", nrow(tudor)))

# The TUDOR file has: participant.eid, participant.p30640_i0 (ApoB mg/dL),
# participant.p30780_i0 (LDL), APOB (derived g/L), LDL_untreated, etc.
# Also: is_fh_genetic, participant.p6150_i0, participant.p20002_i0_a*,
# participant.p20003_i0_a*, BMI_imputed, Age_at_LDL1, Gender_num

# Non-FH: exclude FH genetic
fh_eids <- as.character(fh_ukb$id)
tudor$eid_str <- as.character(tudor$participant.eid)
nonfh <- tudor %>% filter(is_fh_genetic == 0 | is.na(is_fh_genetic))
# Also exclude by eid match
nonfh <- nonfh %>% filter(!(eid_str %in% fh_eids))
cat(sprintf("  Non-FH after excluding FH: N=%d\n", nrow(nonfh)))

# Compute ApoB/LDL ratio
# participant.p30640_i0 is ApoB in g/L (UKB standard)
# participant.p30780_i0 is LDL in mmol/L
nonfh$apob_gL <- nonfh$participant.p30640_i0
nonfh$ldl_mmol <- nonfh$participant.p30780_i0
nonfh$apob_ldl <- ifelse(is.finite(nonfh$apob_gL) & is.finite(nonfh$ldl_mmol) & nonfh$ldl_mmol > 0,
                          nonfh$apob_gL / nonfh$ldl_mmol, NA)

# Define ASCVD: p6150_i0 (1=MI, 2=angina, 3=stroke)
nonfh$sr_ascvd <- nonfh$participant.p6150_i0 %in% c(1, 2, 3)

# p20002 illness codes: 1065=MI, 1066=HF, 1081=stroke
p20002_cols <- grep("^participant\\.p20002", names(nonfh), value = TRUE)
nonfh$ill_ascvd <- apply(nonfh[p20002_cols], 1, function(x) any(x %in% c(1065, 1066, 1081), na.rm = TRUE))

nonfh$ascvd <- ifelse(nonfh$sr_ascvd | nonfh$ill_ascvd, 1L, 0L)

# Diabetes from p20002 (1220=T2DM, 1223=T1DM)
nonfh$diabetes <- apply(nonfh[p20002_cols], 1, function(x) any(x %in% c(1220, 1223), na.rm = TRUE))
nonfh$diabetes <- as.integer(nonfh$diabetes)

# Statin from p20003
statin_codes <- c(1141146234, 1140861958, 1141192414, 1140888594, 1140888648, 1141146138)
p20003_cols <- grep("^participant\\.p20003", names(nonfh), value = TRUE)
nonfh$on_statin <- apply(nonfh[p20003_cols], 1, function(x) any(x %in% statin_codes, na.rm = TRUE))
nonfh$on_statin <- as.integer(nonfh$on_statin)

# Smoking from p20107 (0=never, 1=previous, 2=current)
nonfh$smoking <- ifelse(nonfh$participant.p20107_i0 %in% c(1, 2), 1L, 0L)

# Build clean non-FH dataset
nonfh_clean <- nonfh %>%
  filter(is.finite(apob_ldl) & apob_ldl > 0 & apob_ldl < 2) %>%
  transmute(
    id = eid_str, cohort = "UKB-nonFH",
    age = Age_at_LDL1,
    sex = Gender_num,
    bmi = BMI_imputed,
    ldl = ldl_mmol, hdl = participant.p30760_i0, tc = participant.p30690_i0,
    tg = participant.p30870_i0,
    apob = apob_gL, apob_ldl = apob_ldl,
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
cat("  ANALYSIS 1: ApoB/LDL-C in ASCVD+ vs ASCVD- WITHIN FH\n")
cat("===================================================================\n")

fh_pos <- fh_pooled %>% filter(ascvd == 1, is.finite(apob_ldl))
fh_neg <- fh_pooled %>% filter(ascvd == 0, is.finite(apob_ldl))

cat(sprintf("  FH ASCVD+: N=%d, mean=%.4f (SD %.4f), median=%.4f\n",
            nrow(fh_pos), mean(fh_pos$apob_ldl), sd(fh_pos$apob_ldl), median(fh_pos$apob_ldl)))
cat(sprintf("  FH ASCVD-: N=%d, mean=%.4f (SD %.4f), median=%.4f\n",
            nrow(fh_neg), mean(fh_neg$apob_ldl), sd(fh_neg$apob_ldl), median(fh_neg$apob_ldl)))

wt1 <- wilcox.test(fh_pos$apob_ldl, fh_neg$apob_ldl)
tt1 <- t.test(fh_pos$apob_ldl, fh_neg$apob_ldl)
smd1 <- (mean(fh_pos$apob_ldl) - mean(fh_neg$apob_ldl)) /
  sqrt((sd(fh_pos$apob_ldl)^2 + sd(fh_neg$apob_ldl)^2) / 2)
cat(sprintf("  Wilcoxon p = %g, t-test p = %g, SMD = %.3f\n", wt1$p.value, tt1$p.value, smd1))

# Binary >=0.31
fh_all <- fh_pooled %>% filter(is.finite(apob_ldl))
fh_all$apob_high <- ifelse(fh_all$apob_ldl >= 0.31, 1, 0)
tab1 <- table(fh_all$apob_high, fh_all$ascvd)
cat("\n  Binary >=0.31:\n"); print(tab1)
ft1 <- fisher.test(tab1)
cat(sprintf("  Fisher OR = %.2f (%.2f-%.2f), p = %g\n", ft1$estimate, ft1$conf.int[1], ft1$conf.int[2], ft1$p.value))
r_high <- sum(fh_all$ascvd[fh_all$apob_high==1]==1) / sum(fh_all$apob_high==1)
r_low <- sum(fh_all$ascvd[fh_all$apob_high==0]==1) / sum(fh_all$apob_high==0)
cat(sprintf("  ASCVD rate: >=0.31 = %.1f%%, <0.31 = %.1f%%\n\n", 100*r_high, 100*r_low))


# ══════════════════════════════════════════════════════════════════════
# ANALYSIS 2: ASCVD+ vs ASCVD- WITHIN NON-FH
# ══════════════════════════════════════════════════════════════════════
cat("===================================================================\n")
cat("  ANALYSIS 2: ApoB/LDL-C in ASCVD+ vs ASCVD- WITHIN NON-FH\n")
cat("===================================================================\n")

nf_pos <- nonfh_clean %>% filter(ascvd == 1, is.finite(apob_ldl))
nf_neg <- nonfh_clean %>% filter(ascvd == 0, is.finite(apob_ldl))

cat(sprintf("  NonFH ASCVD+: N=%d, mean=%.4f (SD %.4f), median=%.4f\n",
            nrow(nf_pos), mean(nf_pos$apob_ldl), sd(nf_pos$apob_ldl), median(nf_pos$apob_ldl)))
cat(sprintf("  NonFH ASCVD-: N=%d, mean=%.4f (SD %.4f), median=%.4f\n",
            nrow(nf_neg), mean(nf_neg$apob_ldl), sd(nf_neg$apob_ldl), median(nf_neg$apob_ldl)))

wt2 <- wilcox.test(nf_pos$apob_ldl, nf_neg$apob_ldl)
tt2 <- t.test(nf_pos$apob_ldl, nf_neg$apob_ldl)
smd2 <- (mean(nf_pos$apob_ldl) - mean(nf_neg$apob_ldl)) /
  sqrt((sd(nf_pos$apob_ldl)^2 + sd(nf_neg$apob_ldl)^2) / 2)
cat(sprintf("  Wilcoxon p = %g, t-test p = %g, SMD = %.3f\n", wt2$p.value, tt2$p.value, smd2))

nf_all <- nonfh_clean %>% filter(is.finite(apob_ldl))
nf_all$apob_high <- ifelse(nf_all$apob_ldl >= 0.31, 1, 0)
tab2 <- table(nf_all$apob_high, nf_all$ascvd)
cat("\n  Binary >=0.31:\n"); print(tab2)
ft2 <- fisher.test(tab2)
cat(sprintf("  Fisher OR = %.2f (%.2f-%.2f), p = %g\n\n", ft2$estimate, ft2$conf.int[1], ft2$conf.int[2], ft2$p.value))


# ══════════════════════════════════════════════════════════════════════
# ANALYSIS 3: FH vs NON-FH overall
# ══════════════════════════════════════════════════════════════════════
cat("===================================================================\n")
cat("  ANALYSIS 3: FH vs NON-FH overall ApoB/LDL-C\n")
cat("===================================================================\n")

fh_apob <- fh_pooled %>% filter(is.finite(apob_ldl))
nf_apob <- nonfh_clean %>% filter(is.finite(apob_ldl))

cat(sprintf("  FH: N=%d, mean=%.4f (SD %.4f), median=%.4f\n",
            nrow(fh_apob), mean(fh_apob$apob_ldl), sd(fh_apob$apob_ldl), median(fh_apob$apob_ldl)))
cat(sprintf("  NonFH: N=%d, mean=%.4f (SD %.4f), median=%.4f\n",
            nrow(nf_apob), mean(nf_apob$apob_ldl), sd(nf_apob$apob_ldl), median(nf_apob$apob_ldl)))

wt3 <- wilcox.test(fh_apob$apob_ldl, nf_apob$apob_ldl)
tt3 <- t.test(fh_apob$apob_ldl, nf_apob$apob_ldl)
smd3 <- (mean(fh_apob$apob_ldl) - mean(nf_apob$apob_ldl)) /
  sqrt((sd(fh_apob$apob_ldl)^2 + sd(nf_apob$apob_ldl)^2) / 2)
cat(sprintf("  Wilcoxon p = %g, t-test p = %g, SMD = %.3f\n\n", wt3$p.value, tt3$p.value, smd3))


# ══════════════════════════════════════════════════════════════════════
# ANALYSIS 4: MATCHED FH+ASCVD vs NON-FH+ASCVD
# ══════════════════════════════════════════════════════════════════════
cat("===================================================================\n")
cat("  ANALYSIS 4: PROPENSITY-MATCHED FH+ASCVD vs NON-FH+ASCVD\n")
cat("===================================================================\n")

has_matchit <- requireNamespace("MatchIt", quietly = TRUE)
if (!has_matchit) {
  cat("  Installing MatchIt...\n")
  install.packages("MatchIt", repos = "https://cloud.r-project.org", quiet = TRUE)
}
library(MatchIt)

fh_ascvd <- fh_pooled %>%
  filter(ascvd == 1, is.finite(apob_ldl), is.finite(bmi), is.finite(tg), is.finite(age)) %>%
  mutate(diabetes = ifelse(is.na(diabetes), 0, diabetes),
         on_statin = ifelse(is.na(on_statin), 0, on_statin))

nf_ascvd <- nonfh_clean %>%
  filter(ascvd == 1, is.finite(apob_ldl), is.finite(bmi), is.finite(tg), is.finite(age)) %>%
  mutate(diabetes = ifelse(is.na(diabetes), 0, diabetes),
         on_statin = ifelse(is.na(on_statin), 0, on_statin))

cat(sprintf("  FH ASCVD+ (complete): N=%d\n", nrow(fh_ascvd)))
cat(sprintf("  NonFH ASCVD+ (complete): N=%d\n", nrow(nf_ascvd)))

combined <- bind_rows(
  fh_ascvd %>% mutate(fh_status = 1L),
  nf_ascvd %>% mutate(fh_status = 0L)
)

m <- matchit(fh_status ~ age + bmi + tg + diabetes + on_statin,
             data = combined, method = "nearest", ratio = 3, caliper = 0.2)

cat("\n  Match summary:\n")
print(summary(m))

matched <- match.data(m)
cat(sprintf("\n  Matched: N=%d (FH=%d, NonFH=%d)\n",
            nrow(matched), sum(matched$fh_status==1), sum(matched$fh_status==0)))

m_fh <- matched %>% filter(fh_status == 1)
m_nf <- matched %>% filter(fh_status == 0)

cat(sprintf("\n  Matched FH ASCVD+: mean=%.4f (SD %.4f), median=%.4f\n",
            mean(m_fh$apob_ldl), sd(m_fh$apob_ldl), median(m_fh$apob_ldl)))
cat(sprintf("  Matched NonFH ASCVD+: mean=%.4f (SD %.4f), median=%.4f\n",
            mean(m_nf$apob_ldl), sd(m_nf$apob_ldl), median(m_nf$apob_ldl)))

wt4 <- wilcox.test(m_fh$apob_ldl, m_nf$apob_ldl)
tt4 <- t.test(m_fh$apob_ldl, m_nf$apob_ldl)
smd4 <- (mean(m_fh$apob_ldl) - mean(m_nf$apob_ldl)) /
  sqrt((sd(m_fh$apob_ldl)^2 + sd(m_nf$apob_ldl)^2) / 2)
cat(sprintf("  Wilcoxon p = %g, t-test p = %g, SMD = %.3f\n", wt4$p.value, tt4$p.value, smd4))

# Binary threshold
matched$fh_factor <- factor(matched$fh_status, levels=c(0,1), labels=c("NonFH","FH"))
m_fh_high <- 100 * sum(m_fh$apob_ldl >= 0.31) / nrow(m_fh)
m_nf_high <- 100 * sum(m_nf$apob_ldl >= 0.31) / nrow(m_nf)
cat(sprintf("\n  ApoB/LDL >=0.31: FH %.1f%% vs NonFH %.1f%%\n", m_fh_high, m_nf_high))

tab4 <- table(
  factor(ifelse(matched$apob_ldl >= 0.31, "High", "Low"), levels=c("Low","High")),
  matched$fh_factor
)
print(tab4)
ft4 <- fisher.test(tab4)
cat(sprintf("  Fisher OR = %.2f (%.2f-%.2f), p = %g\n", ft4$estimate, ft4$conf.int[1], ft4$conf.int[2], ft4$p.value))

# TableOne
vars <- c("age","sex","bmi","ldl","hdl","tc","tg","apob","apob_ldl","diabetes","on_statin","smoking")
catvars <- c("sex","diabetes","on_statin","smoking")
tab_m <- CreateTableOne(vars=vars, strata="fh_factor", data=matched, factorVars=catvars, test=TRUE, smd=TRUE)
cat("\n  Matched TableOne:\n")
print(tab_m, smd=TRUE, test=TRUE, showAllLevels=FALSE)

matched_out <- print(tab_m, smd=TRUE, test=TRUE, showAllLevels=FALSE, printToggle=FALSE)
write.csv(matched_out, "output/final/Table_FH_vs_NonFH_ASCVD_matched.csv")
cat("\n  Saved: output/final/Table_FH_vs_NonFH_ASCVD_matched.csv\n")


# ══════════════════════════════════════════════════════════════════════
# ANALYSIS 5: UNMATCHED TABLE
# ══════════════════════════════════════════════════════════════════════
cat("\n===================================================================\n")
cat("  ANALYSIS 5: UNMATCHED FH vs NON-FH\n")
cat("===================================================================\n")

all_comb <- bind_rows(
  fh_pooled %>% filter(is.finite(apob_ldl)) %>% mutate(group = "FH"),
  nonfh_clean %>% filter(is.finite(apob_ldl)) %>% mutate(group = "NonFH")
)
vars_um <- c("age","sex","bmi","ldl","hdl","tc","tg","apob","apob_ldl","diabetes","on_statin","smoking","ascvd")
catvars_um <- c("sex","diabetes","on_statin","smoking","ascvd")
tab_um <- CreateTableOne(vars=vars_um, strata="group", data=all_comb, factorVars=catvars_um, test=TRUE, smd=TRUE)
cat("\n  Unmatched TableOne:\n")
print(tab_um, smd=TRUE, test=TRUE, showAllLevels=FALSE)

um_out <- print(tab_um, smd=TRUE, test=TRUE, showAllLevels=FALSE, printToggle=FALSE)
write.csv(um_out, "output/final/Table_FH_vs_NonFH_unmatched.csv")
cat("\n  Saved: output/final/Table_FH_vs_NonFH_unmatched.csv\n")

cat("\n===================================================================\n")
cat("  DONE\n")
cat("===================================================================\n")

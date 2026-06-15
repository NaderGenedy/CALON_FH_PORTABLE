##############################################################################
# TUDOR — New UKB files: Death + Diabetes + Xanthomata
# Adds:
#   - CV mortality by TUDOR/LDL decile
#   - T2DM direct question (better than ICD-10 E11)
#   - HbA1c as continuous insulin-resistance measure
#   - Xanthomata prevalence by TUDOR decile
##############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
})

RAW  <- "E:/"
OUT  <- "C:/Users/nader/Downloads/calon_ukb_pipeline/tudor_loco_output/tudor_flaw_results"

# ── Load new files ──────────────────────────────────────────────────────────
cat("Loading new files...\n")

death  <- read_csv(file.path(RAW, "ukb_reviewer_death.csv"),      show_col_types=FALSE)
diab   <- read_csv(file.path(RAW, "ukb_reviewer_diabetes.csv"),   show_col_types=FALSE)
xanth  <- read_csv(file.path(RAW, "ukb_reviewer_xanthomata.csv"), show_col_types=FALSE)

# Rename EID
death  <- rename(death,  eid = participant.eid)
diab   <- rename(diab,   eid = participant.eid)
xanth  <- rename(xanth,  eid = participant.eid)

cat(sprintf("death: %d | diabetes: %d | xanthomata: %d\n",
            nrow(death), nrow(diab), nrow(xanth)))
cat("Death cols:",   paste(names(death),  collapse=", "), "\n")
cat("Diabetes cols:",paste(names(diab),   collapse=", "), "\n")
cat("Xanthomata cols:",paste(names(xanth),collapse=", "), "\n\n")

##############################################################################
# 1. DEATH — cause of death ICD-10
##############################################################################

cat("=== Death analysis ===\n")

# p40000 = date of death, p40001 = cause of death (ICD-10), p40007 = age at death
death_col_icd  <- grep("p40001", names(death), value=TRUE)[1]
death_col_date <- grep("p40000", names(death), value=TRUE)[1]
death_col_age  <- grep("p40007", names(death), value=TRUE)[1]

cat(sprintf("Using: ICD=%s | Date=%s | Age=%s\n",
            death_col_icd, death_col_date, death_col_age))

# CV death ICD-10 codes
cv_death_pat <- "^I[0-9]|^G45"

death <- death %>%
  mutate(
    death_icd      = .data[[death_col_icd]],
    age_at_death   = .data[[death_col_age]],
    any_death      = !is.na(death_icd) & death_icd != "",
    cv_death       = any_death & grepl(cv_death_pat, death_icd),
    premature_cv_death = cv_death & !is.na(age_at_death) & age_at_death < 60
  ) %>%
  select(eid, any_death, cv_death, premature_cv_death, age_at_death)

cat(sprintf("Any death:           n=%d (%.1f%%)\n",
            sum(death$any_death), 100*mean(death$any_death)))
cat(sprintf("CV death:            n=%d (%.1f%%)\n",
            sum(death$cv_death), 100*mean(death$cv_death)))
cat(sprintf("Premature CV death (<60): n=%d (%.1f%%)\n",
            sum(death$premature_cv_death), 100*mean(death$premature_cv_death)))

##############################################################################
# 2. DIABETES — direct question p2443 + HbA1c p30750
##############################################################################

cat("\n=== Diabetes — direct question ===\n")

t2dm_col  <- grep("p2443_i0", names(diab), value=TRUE)[1]
hba1c_col <- grep("p30750_i0", names(diab), value=TRUE)[1]

cat(sprintf("T2DM direct: %s | HbA1c: %s\n", t2dm_col, hba1c_col))

# p2443: 1=Yes (diabetes diagnosed), 0=No, -1/-3=prefer not/don't know
diab <- diab %>%
  mutate(
    t2dm_direct = case_when(
      .data[[t2dm_col]] == 1 ~ TRUE,
      .data[[t2dm_col]] == 0 ~ FALSE,
      TRUE ~ NA
    ),
    hba1c = .data[[hba1c_col]]
  ) %>%
  select(eid, t2dm_direct, hba1c)

cat(sprintf("T2DM (direct question): n=%d (%.1f%%)\n",
            sum(diab$t2dm_direct, na.rm=TRUE),
            100*mean(diab$t2dm_direct, na.rm=TRUE)))
cat(sprintf("HbA1c available: n=%d | mean=%.1f mmol/mol\n",
            sum(!is.na(diab$hba1c)),
            mean(diab$hba1c, na.rm=TRUE)))

##############################################################################
# 3. XANTHOMATA — p4631
##############################################################################

cat("\n=== Xanthomata ===\n")

xanth_col <- grep("p4631_i0", names(xanth), value=TRUE)[1]
cat(sprintf("Xanthomata column: %s\n", xanth_col))
cat("Values (sample):\n"); print(table(xanth[[xanth_col]], useNA="ifany")[1:10])

# p4631: typically 1=Yes, 0=No (self-reported)
xanth <- xanth %>%
  mutate(
    xanthomata = case_when(
      .data[[xanth_col]] == 1 ~ TRUE,
      .data[[xanth_col]] == 0 ~ FALSE,
      TRUE ~ NA
    )
  ) %>%
  select(eid, xanthomata)

cat(sprintf("Xanthomata present: n=%d (%.2f%%)\n",
            sum(xanth$xanthomata, na.rm=TRUE),
            100*mean(xanth$xanthomata, na.rm=TRUE)))

##############################################################################
# 4. LOAD EXISTING SUPPLEMENT + MERGE
##############################################################################

cat("\n=== Merging with existing supplement ===\n")

lc <- read_csv(file.path(OUT, "ukb_supplement_lc_cohort.csv"),
               show_col_types=FALSE)
cat(sprintf("Existing LC cohort: %d rows\n", nrow(lc)))

lc2 <- lc %>%
  left_join(death, by="eid") %>%
  left_join(diab,  by="eid") %>%
  left_join(xanth, by="eid")

cat(sprintf("After merge: %d rows | %d cols\n", nrow(lc2), ncol(lc2)))

##############################################################################
# 5. CV MORTALITY BY LDL DECILE
##############################################################################

cat("\n=== CV mortality by LDL decile ===\n")

cv_ldl <- lc2 %>%
  filter(!is.na(ldl_baseline), !is.na(cv_death)) %>%
  mutate(ldl_decile = ntile(ldl_baseline, 10)) %>%
  group_by(ldl_decile) %>%
  summarise(
    n             = n(),
    ldl_mean      = round(mean(ldl_baseline), 2),
    pct_any_death = round(100*mean(any_death,      na.rm=TRUE), 2),
    pct_cv_death  = round(100*mean(cv_death,       na.rm=TRUE), 2),
    pct_prem_cv   = round(100*mean(premature_cv_death, na.rm=TRUE), 2),
    .groups       = "drop"
  )
cat("\nCV mortality by LDL decile:\n"); print(cv_ldl)
write_csv(cv_ldl, file.path(OUT, "cv_mortality_by_ldl_decile.csv"))

##############################################################################
# 6. T2DM DIRECT QUESTION VS ICD-10 — COMPARISON
##############################################################################

cat("\n=== T2DM direct vs ICD-10 comparison ===\n")

t2dm_compare <- lc2 %>%
  filter(!is.na(t2dm), !is.na(t2dm_direct)) %>%
  summarise(
    n                    = n(),
    t2dm_icd10_pct       = round(100*mean(t2dm, na.rm=TRUE), 1),
    t2dm_direct_pct      = round(100*mean(t2dm_direct, na.rm=TRUE), 1),
    both_positive        = sum(t2dm & t2dm_direct, na.rm=TRUE),
    icd10_only           = sum(t2dm & !t2dm_direct, na.rm=TRUE),
    direct_only          = sum(!t2dm & t2dm_direct, na.rm=TRUE),
    neither              = sum(!t2dm & !t2dm_direct, na.rm=TRUE)
  )
cat("\nT2DM source comparison:\n"); print(t2dm_compare)
write_csv(t2dm_compare, file.path(OUT, "t2dm_source_comparison.csv"))

# Re-run TRG analysis with DIRECT T2DM (better measure)
trg_direct <- lc2 %>%
  filter(!is.na(trig_filter_proxy), !is.na(t2dm_direct),
         !is.na(bmi), !is.na(age_baseline), is.finite(trig_filter_proxy))

cat(sprintf("\nTRG analysis with direct T2DM: n=%d\n", nrow(trg_direct)))

if (nrow(trg_direct) > 1000) {
  lm_trg2 <- lm(trig_filter_proxy ~ t2dm_direct + bmi + alcohol_cat + age_baseline,
                data = trg_direct)
  t2dm_beta2 <- coef(lm_trg2)["t2dm_directTRUE"]
  t2dm_ci2   <- confint(lm_trg2)["t2dm_directTRUE",]
  cat(sprintf("T2DM (direct) β = %.3f (%.3f to %.3f)\n",
              t2dm_beta2, t2dm_ci2[1], t2dm_ci2[2]))
  cat(sprintf("Compare with ICD-10 T2DM β = -0.477\n"))

  # Group comparison
  tf_direct_t2dm    <- trg_direct$trig_filter_proxy[trg_direct$t2dm_direct == TRUE]
  tf_direct_no_t2dm <- trg_direct$trig_filter_proxy[trg_direct$t2dm_direct == FALSE]
  cat(sprintf("Mean TRG: Direct-T2DM+ = %.3f vs T2DM- = %.3f\n",
              mean(tf_direct_t2dm), mean(tf_direct_no_t2dm)))
}

##############################################################################
# 7. XANTHOMATA BY LDL DECILE
##############################################################################

cat("\n=== Xanthomata by LDL decile ===\n")

xanth_ldl <- lc2 %>%
  filter(!is.na(ldl_baseline), !is.na(xanthomata)) %>%
  mutate(ldl_decile = ntile(ldl_baseline, 10)) %>%
  group_by(ldl_decile) %>%
  summarise(
    n           = n(),
    ldl_mean    = round(mean(ldl_baseline), 2),
    pct_xanth   = round(100*mean(xanthomata, na.rm=TRUE), 2),
    .groups     = "drop"
  )
cat("\nXanthomata by LDL decile:\n"); print(xanth_ldl)
write_csv(xanth_ldl, file.path(OUT, "xanthomata_by_ldl_decile.csv"))

overall_xanth_pct <- 100*mean(lc2$xanthomata, na.rm=TRUE)
cat(sprintf("\nOverall xanthomata prevalence in LC cohort: %.2f%%\n", overall_xanth_pct))

##############################################################################
# 8. SAVE UPDATED COHORT
##############################################################################

write_csv(lc2, file.path(OUT, "ukb_supplement_lc_cohort_v2.csv"))

cat("\n========================================================\n")
cat("SUMMARY — NEW FILE ANALYSES COMPLETE\n")
cat("========================================================\n")
cat(sprintf("T2DM by ICD-10: %.1f%% | by direct question: %.1f%%\n",
            100*mean(lc2$t2dm, na.rm=TRUE),
            100*mean(lc2$t2dm_direct, na.rm=TRUE)))
cat(sprintf("CV deaths: %.1f%% | Premature CV death (<60): %.1f%%\n",
            100*mean(lc2$cv_death, na.rm=TRUE),
            100*mean(lc2$premature_cv_death, na.rm=TRUE)))
cat(sprintf("Xanthomata: %.2f%%\n", overall_xanth_pct))
cat("\nNew output files:\n")
cat("  cv_mortality_by_ldl_decile.csv\n")
cat("  t2dm_source_comparison.csv\n")
cat("  xanthomata_by_ldl_decile.csv\n")
cat("  ukb_supplement_lc_cohort_v2.csv\n")

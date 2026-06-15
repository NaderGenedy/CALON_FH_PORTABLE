################################################################################
# AF3 MANUSCRIPT TABLES - Nature/Lancet Calibre
# 7 Publication Tables for AlphaFold3 Familial Hypercholesterolaemia Manuscript
# Author: Dr Nader Genedy
# Date: March 2026
################################################################################

# ==============================================================================
# Section 0: Packages & Data Loading
# ==============================================================================

required_pkgs <- c("tableone", "pROC", "data.table")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tableone)
library(pROC)
library(data.table)

# --- Paths ---
base_dir   <- "C:/Users/nader/Downloads/calon_ukb_pipeline"
analy_dir  <- file.path(base_dir, "alphafold", "analysis")
nmr_dir    <- file.path(base_dir, "alphafold")
out_dir    <- file.path(analy_dir, "tables")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- Load core datasets ---
cat("Loading core datasets...\n")

comp_sss <- fread(file.path(analy_dir, "comprehensive_sss_analysis.csv"),
                  na.strings = c("", "NA", "NaN"))

nobel <- fread(file.path(analy_dir, "nobel_multimodal_integration.csv"),
               na.strings = c("", "NA", "NaN"))

sss_scores <- fread(file.path(analy_dir, "structural_severity_scores.csv"),
                    na.strings = c("", "NA", "NaN"))

foldx <- fread(file.path(analy_dir, "foldx_ddg_results.csv"),
               na.strings = c("", "NA", "NaN"))

llt_response <- fread(file.path(analy_dir, "variant_llt_response.csv"),
                      na.strings = c("", "NA", "NaN"))

loco_cv <- fread(file.path(analy_dir, "loco_cv_pooled_sss_results.csv"),
                 na.strings = c("", "NA", "NaN"))

cat("  comprehensive_sss_analysis:", nrow(comp_sss), "rows\n")
cat("  nobel_multimodal_integration:", nrow(nobel), "rows\n")
cat("  structural_severity_scores:", nrow(sss_scores), "rows\n")
cat("  foldx_ddg_results:", nrow(foldx), "rows\n")
cat("  variant_llt_response:", nrow(llt_response), "rows\n")
cat("  loco_cv_pooled_sss_results:", nrow(loco_cv), "rows\n")

# --- Helper functions ---
fmt2 <- function(x) formatC(x, format = "f", digits = 2)
fmt1 <- function(x) formatC(x, format = "f", digits = 1)
mean_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  paste0(fmt2(mean(x)), " (", fmt2(sd(x)), ")")
}
median_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  q <- quantile(x, c(0.25, 0.5, 0.75))
  paste0(fmt2(q[2]), " [", fmt2(q[1]), "-", fmt2(q[3]), "]")
}
pct <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  fmt1(100 * mean(x))
}


# ==============================================================================
# Section 1: Table 1 - Baseline Characteristics by Gene
# ==============================================================================

generate_table1 <- function() {
  cat("\n--- Generating Table 1: Baseline Characteristics by Gene ---\n")

  dt <- copy(comp_sss)

  # Ensure gene is clean factor

  dt[, gene := trimws(gene)]
  dt <- dt[gene %in% c("LDLR", "APOB", "PCSK9")]

  # Convert binary columns to factors for tableone
  binary_vars <- c("on_statin", "ezetimibe", "pcsk9i", "dm", "hypertension",
                   "smoking", "xanthomata", "ascvd", "miacs")
  for (v in binary_vars) {
    if (v %in% names(dt)) {
      dt[[v]] <- factor(ifelse(dt[[v]] == 1 | dt[[v]] == TRUE, "Yes", "No"),
                        levels = c("No", "Yes"))
    }
  }

  # Sex: 0=Male, 1=Female based on typical coding
  if ("sex" %in% names(dt)) {
    dt[, sex := factor(ifelse(sex == 1, "Female", "Male"),
                       levels = c("Male", "Female"))]
  }

  # tx_intensity: recode to high intensity flag
  if ("tx_intensity" %in% names(dt)) {
    dt[, high_intensity := factor(
      ifelse(!is.na(tx_intensity) & tx_intensity == 1, "Yes", "No"),
      levels = c("No", "Yes")
    )]
  }

  # Simon Broome: 1 = possible, 2 = definite
  if ("simon_broome" %in% names(dt)) {
    dt[, simon_broome_cat := factor(
      ifelse(is.na(simon_broome), "Not assessed",
             ifelse(simon_broome >= 2, "Definite",
                    ifelse(simon_broome >= 1, "Possible", "Not met"))),
      levels = c("Not met", "Possible", "Definite", "Not assessed")
    )]
  }

  # Define variables for tableone
  cont_vars <- c("age", "bmi", "tc1", "ldl1", "hdl1", "tg1",
                 "apob", "apoa1", "lpa", "last_ldl", "ldl_pct_change", "sss")
  cat_vars <- c("sex", "on_statin", "ezetimibe", "pcsk9i", "high_intensity",
                "dm", "hypertension", "smoking", "xanthomata",
                "simon_broome_cat", "ascvd", "miacs")

  # Keep only variables that exist
  cont_vars <- cont_vars[cont_vars %in% names(dt)]
  cat_vars  <- cat_vars[cat_vars %in% names(dt)]
  all_vars  <- c(cont_vars, cat_vars)

  # Create TableOne
  tab1 <- CreateTableOne(
    vars       = all_vars,
    strata     = "gene",
    data       = as.data.frame(dt),
    factorVars = cat_vars,
    addOverall = TRUE,
    test       = TRUE
  )

  # Print with non-normal continuous (Kruskal-Wallis)
  tab1_print <- print(tab1,
                      nonnormal   = cont_vars,
                      exact       = NULL,
                      smd         = FALSE,
                      showAllLevels = FALSE,
                      printToggle = FALSE,
                      noSpaces    = TRUE)

  # Convert to data.frame for export
  tab1_df <- as.data.frame(tab1_print)
  tab1_df <- cbind(Variable = rownames(tab1_df), tab1_df)
  rownames(tab1_df) <- NULL

  # Rename columns for clarity
  names(tab1_df)[names(tab1_df) == "p"] <- "P_value"
  names(tab1_df)[names(tab1_df) == "test"] <- "Test_type"

  out_path <- file.path(out_dir, "Table1_Baseline_Characteristics.csv")
  fwrite(tab1_df, out_path)
  cat("  Saved:", out_path, "\n")
  cat("  Dimensions:", nrow(tab1_df), "rows x", ncol(tab1_df), "cols\n")

  return(invisible(tab1_df))
}


# ==============================================================================
# Section 2: Table 2 - Top 30 Variant Catalogue
# ==============================================================================

generate_table2 <- function() {
  cat("\n--- Generating Table 2: Top 30 Variant Catalogue ---\n")

  # Start from nobel (44 variants with aggregated patient-level data)
  dt <- copy(nobel)

  # Merge with sss_scores for additional structural info
  sss_merge <- sss_scores[, .(variant_id, variant_type, plddt, interface_dist,
                               domain_severity)]
  # nobel uses 'variant', sss_scores uses 'variant_id'
  dt <- merge(dt, sss_merge, by.x = "variant", by.y = "variant_id",
              all.x = TRUE, suffixes = c("", ".sss"))

  # Use plddt from sss_scores if missing in nobel
  if ("plddt.sss" %in% names(dt)) {
    dt[is.na(plddt), plddt := plddt.sss]
    dt[, plddt.sss := NULL]
  }
  if ("interface_dist.sss" %in% names(dt)) {
    dt[is.na(interface_dist), interface_dist := interface_dist.sss]
    dt[, interface_dist.sss := NULL]
  }

  # Merge with foldx for ddG where available
  foldx_merge <- foldx[, .(variant_id, ddG_kcal_mol)]
  dt <- merge(dt, foldx_merge, by.x = "variant", by.y = "variant_id",
              all.x = TRUE)
  # Use foldx ddG if nobel ddG is missing
  dt[is.na(ddG), ddG := ddG_kcal_mol]
  dt[, ddG_kcal_mol := NULL]

  # Sort by n (patient count) descending, take top 30
  setorder(dt, -n)
  dt <- head(dt, 30)

  # Format output table
  out <- data.table(
    Variant             = dt$variant,
    Gene                = dt$gene,
    Domain              = dt$domain,
    SSS                 = fmt2(dt$sss),
    ddG_kcal_mol        = ifelse(is.na(dt$ddG), "NA", fmt2(dt$ddG)),
    pLDDT               = ifelse(is.na(dt$plddt), "NA", fmt2(dt$plddt)),
    N_patients          = dt$n,
    Mean_LDL_mmol_L     = ifelse(is.na(dt$mean_ldl), "NA", fmt2(dt$mean_ldl)),
    Mean_ApoB_g_L       = ifelse(is.na(dt$mean_apob), "NA", fmt2(dt$mean_apob)),
    ASCVD_pct           = ifelse(is.na(dt$ascvd_pct), "NA", fmt1(dt$ascvd_pct)),
    Mean_Treatment_Response = ifelse(is.na(dt$mean_response), "NA",
                                     fmt2(dt$mean_response)),
    Cluster             = dt$cluster
  )

  out_path <- file.path(out_dir, "Table2_Top30_Variant_Catalogue.csv")
  fwrite(out, out_path)
  cat("  Saved:", out_path, "\n")
  cat("  Dimensions:", nrow(out), "rows x", ncol(out), "cols\n")

  return(invisible(out))
}


# ==============================================================================
# Section 3: Table 3 - Domain-Specific Phenotypes
# ==============================================================================

generate_table3 <- function() {
  cat("\n--- Generating Table 3: Domain-Specific Phenotypes ---\n")

  dt <- copy(comp_sss)
  dt[, domain := trimws(domain)]

  # Aggregate by domain
  domain_stats <- dt[, .(
    N_variants      = uniqueN(mutation),
    N_patients      = .N,
    Mean_SSS        = mean_sd(sss),
    Mean_LDL        = mean_sd(ldl1),
    Mean_ApoB       = mean_sd(apob),
    ApoB_LDL_Ratio  = ifelse(
      sum(!is.na(apob) & !is.na(ldl1)) > 0,
      fmt2(mean(apob[!is.na(apob) & !is.na(ldl1)] /
                ldl1[!is.na(apob) & !is.na(ldl1)], na.rm = TRUE)),
      "NA"
    ),
    Mean_Lpa        = ifelse(sum(!is.na(lpa)) > 0,
                             fmt2(mean(lpa, na.rm = TRUE)), "NA"),
    ASCVD_pct       = pct(ascvd),
    Xanthomata_pct  = pct(xanthomata),
    Mean_LDL_Reduction_pct = ifelse(
      sum(!is.na(ldl_pct_change)) > 0,
      fmt2(mean(ldl_pct_change, na.rm = TRUE)),
      "NA"
    )
  ), by = domain]

  # Sort by N_patients descending
  setorder(domain_stats, -N_patients)

  # Rename for publication
  setnames(domain_stats, "domain", "Domain")

  out_path <- file.path(out_dir, "Table3_Domain_Specific_Phenotypes.csv")
  fwrite(domain_stats, out_path)
  cat("  Saved:", out_path, "\n")
  cat("  Dimensions:", nrow(domain_stats), "rows x", ncol(domain_stats), "cols\n")

  return(invisible(domain_stats))
}


# ==============================================================================
# Section 4: Table 4 - 16 Nested Logistic Regression Models
# ==============================================================================

generate_table4 <- function() {
  cat("\n--- Generating Table 4: 16 Nested Logistic Regression Models ---\n")

  dt <- copy(comp_sss)

  # Ensure ascvd is numeric binary
  dt[, ascvd := as.numeric(ascvd)]

  # Ensure gene is factor
  dt[, gene := as.factor(gene)]

  # Ensure all predictor columns are numeric where needed
  num_cols <- c("age", "sex", "sss", "ldl1", "on_statin", "dm",
                "smoking", "hypertension", "xanthomata")
  for (v in num_cols) {
    if (v %in% names(dt)) {
      dt[[v]] <- as.numeric(dt[[v]])
    }
  }

  # Define 16 model formulae (paired: without SSS, then with SSS)
  model_specs <- list(
    list(id = 1,  formula = ascvd ~ age + sex,                                                                         has_sss = FALSE),
    list(id = 2,  formula = ascvd ~ age + sex + sss,                                                                    has_sss = TRUE),
    list(id = 3,  formula = ascvd ~ age + sex + gene,                                                                   has_sss = FALSE),
    list(id = 4,  formula = ascvd ~ age + sex + gene + sss,                                                             has_sss = TRUE),
    list(id = 5,  formula = ascvd ~ age + sex + ldl1,                                                                   has_sss = FALSE),
    list(id = 6,  formula = ascvd ~ age + sex + ldl1 + sss,                                                             has_sss = TRUE),
    list(id = 7,  formula = ascvd ~ age + sex + on_statin,                                                              has_sss = FALSE),
    list(id = 8,  formula = ascvd ~ age + sex + on_statin + sss,                                                        has_sss = TRUE),
    list(id = 9,  formula = ascvd ~ age + sex + ldl1 + on_statin,                                                       has_sss = FALSE),
    list(id = 10, formula = ascvd ~ age + sex + ldl1 + on_statin + sss,                                                 has_sss = TRUE),
    list(id = 11, formula = ascvd ~ age + sex + ldl1 + on_statin + dm + smoking,                                        has_sss = FALSE),
    list(id = 12, formula = ascvd ~ age + sex + ldl1 + on_statin + dm + smoking + sss,                                  has_sss = TRUE),
    list(id = 13, formula = ascvd ~ age + sex + ldl1 + on_statin + dm + smoking + hypertension,                         has_sss = FALSE),
    list(id = 14, formula = ascvd ~ age + sex + ldl1 + on_statin + dm + smoking + hypertension + sss,                   has_sss = TRUE),
    list(id = 15, formula = ascvd ~ age + sex + ldl1 + on_statin + dm + smoking + hypertension + xanthomata,            has_sss = FALSE),
    list(id = 16, formula = ascvd ~ age + sex + ldl1 + on_statin + dm + smoking + hypertension + xanthomata + sss,      has_sss = TRUE)
  )

  results <- list()

  for (spec in model_specs) {
    # Get complete cases for this model
    vars_needed <- all.vars(spec$formula)
    cc <- dt[complete.cases(dt[, ..vars_needed])]

    n_obs    <- nrow(cc)
    n_events <- sum(cc$ascvd == 1, na.rm = TRUE)

    # Fit model
    fit <- tryCatch(
      glm(spec$formula, data = cc, family = binomial(link = "logit")),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      results[[spec$id]] <- data.table(
        Model     = spec$id,
        Formula   = deparse(spec$formula),
        N         = n_obs,
        Events    = n_events,
        AUC       = NA_real_,
        SSS_coef  = NA_real_,
        SSS_OR    = NA_real_,
        SSS_pval  = NA_character_,
        Delta_AUC = NA_real_
      )
      next
    }

    # Compute AUC
    pred <- predict(fit, type = "response")
    roc_obj <- tryCatch(
      roc(cc$ascvd, pred, quiet = TRUE),
      error = function(e) NULL
    )
    auc_val <- if (!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA_real_

    # SSS coefficient and p-value
    sss_coef <- NA_real_
    sss_or   <- NA_real_
    sss_pval <- NA_character_

    if (spec$has_sss) {
      coef_tab <- summary(fit)$coefficients
      if ("sss" %in% rownames(coef_tab)) {
        sss_coef <- coef_tab["sss", "Estimate"]
        sss_or   <- exp(sss_coef)
        raw_p    <- coef_tab["sss", "Pr(>|z|)"]
        sss_pval <- if (raw_p < 0.001) {
          formatC(raw_p, format = "e", digits = 2)
        } else {
          fmt2(raw_p) # use 2 decimal places that's enough for this context
        }
      }
    }

    results[[spec$id]] <- data.table(
      Model     = spec$id,
      Formula   = paste(deparse(spec$formula), collapse = ""),
      N         = n_obs,
      Events    = n_events,
      AUC       = auc_val,
      SSS_coef  = sss_coef,
      SSS_OR    = sss_or,
      SSS_pval  = sss_pval,
      Delta_AUC = NA_real_
    )
  }

  result_dt <- rbindlist(results)

  # Compute Delta-AUC: for each SSS model (even), compare to previous (odd)
  for (i in seq(2, 16, by = 2)) {
    base_auc <- result_dt[Model == (i - 1), AUC]
    sss_auc  <- result_dt[Model == i, AUC]
    if (!is.na(base_auc) && !is.na(sss_auc)) {
      result_dt[Model == i, Delta_AUC := sss_auc - base_auc]
    }
  }

  # Format numeric columns
  result_dt[, AUC       := fifelse(is.na(AUC),       NA_character_, fmt2(AUC))]
  result_dt[, SSS_coef  := fifelse(is.na(SSS_coef),  NA_character_, fmt2(SSS_coef))]
  result_dt[, SSS_OR    := fifelse(is.na(SSS_OR),    NA_character_, fmt2(SSS_OR))]
  result_dt[, Delta_AUC := fifelse(is.na(Delta_AUC), NA_character_, fmt2(Delta_AUC))]

  out_path <- file.path(out_dir, "Table4_Nested_Logistic_Regression.csv")
  fwrite(result_dt, out_path)
  cat("  Saved:", out_path, "\n")
  cat("  Dimensions:", nrow(result_dt), "rows x", ncol(result_dt), "cols\n")

  return(invisible(result_dt))
}


# ==============================================================================
# Section 5: Table 5 - LOCO-CV Results
# ==============================================================================

generate_table5 <- function() {
  cat("\n--- Generating Table 5: LOCO-CV Results ---\n")

  dt <- copy(loco_cv)

  # The loco_cv file has multiple SSS weighting schemes per patient

  # We use sss_original as the primary SSS and compute AUC per variant
  # using logistic regression: ascvd ~ age + sex (+/- sss)

  # Ensure numeric
  dt[, ascvd := as.numeric(ascvd)]
  dt[, age   := as.numeric(age)]
  dt[, sex   := as.numeric(sex)]
  dt[, ldl1  := as.numeric(ldl1)]
  dt[, sss_original := as.numeric(sss_original)]

  # Try to merge cluster info from nobel
  nobel_clust <- nobel[, .(variant, cluster)]
  dt <- merge(dt, nobel_clust, by = "variant", all.x = TRUE)

  # Group by variant - compute base AUC (age+sex) and SSS AUC (age+sex+sss)
  variants <- unique(dt$variant)

  results <- lapply(variants, function(v) {
    sub <- dt[variant == v]
    n   <- nrow(sub)

    # Need at least 10 observations and both outcomes
    if (n < 10 || length(unique(sub$ascvd)) < 2) {
      return(NULL)
    }

    cc <- sub[complete.cases(sub[, .(ascvd, age, sex, sss_original)])]
    if (nrow(cc) < 10 || length(unique(cc$ascvd)) < 2) return(NULL)

    # Base model: age + sex
    fit_base <- tryCatch(
      glm(ascvd ~ age + sex, data = cc, family = binomial),
      error = function(e) NULL
    )

    # SSS model: age + sex + sss_original
    fit_sss <- tryCatch(
      glm(ascvd ~ age + sex + sss_original, data = cc, family = binomial),
      error = function(e) NULL
    )

    auc_base <- NA_real_
    auc_sss  <- NA_real_

    if (!is.null(fit_base)) {
      roc_b <- tryCatch(roc(cc$ascvd, predict(fit_base, type = "response"),
                            quiet = TRUE), error = function(e) NULL)
      if (!is.null(roc_b)) auc_base <- as.numeric(auc(roc_b))
    }

    if (!is.null(fit_sss)) {
      roc_s <- tryCatch(roc(cc$ascvd, predict(fit_sss, type = "response"),
                            quiet = TRUE), error = function(e) NULL)
      if (!is.null(roc_s)) auc_sss <- as.numeric(auc(roc_s))
    }

    # Bootstrap 95% CI for delta AUC
    delta <- NA_real_
    ci_lo <- NA_real_
    ci_hi <- NA_real_

    if (!is.na(auc_base) && !is.na(auc_sss)) {
      delta <- auc_sss - auc_base

      # Simple bootstrap CI for delta AUC
      set.seed(42)
      n_boot  <- 500
      deltas  <- numeric(n_boot)
      for (b in seq_len(n_boot)) {
        idx <- sample(nrow(cc), replace = TRUE)
        boot_data <- cc[idx]
        if (length(unique(boot_data$ascvd)) < 2) {
          deltas[b] <- NA
          next
        }
        fb <- tryCatch(glm(ascvd ~ age + sex, data = boot_data,
                           family = binomial), error = function(e) NULL)
        fs <- tryCatch(glm(ascvd ~ age + sex + sss_original,
                           data = boot_data, family = binomial),
                       error = function(e) NULL)
        if (is.null(fb) || is.null(fs)) {
          deltas[b] <- NA
          next
        }
        rb <- tryCatch(roc(boot_data$ascvd,
                           predict(fb, type = "response"), quiet = TRUE),
                       error = function(e) NULL)
        rs <- tryCatch(roc(boot_data$ascvd,
                           predict(fs, type = "response"), quiet = TRUE),
                       error = function(e) NULL)
        if (is.null(rb) || is.null(rs)) {
          deltas[b] <- NA
          next
        }
        deltas[b] <- as.numeric(auc(rs)) - as.numeric(auc(rb))
      }
      valid_deltas <- deltas[!is.na(deltas)]
      if (length(valid_deltas) >= 50) {
        ci_lo <- quantile(valid_deltas, 0.025)
        ci_hi <- quantile(valid_deltas, 0.975)
      }
    }

    clust <- unique(cc$cluster)
    clust <- clust[!is.na(clust)]
    clust_str <- if (length(clust) > 0) clust[1] else NA_character_

    data.table(
      Variant       = v,
      Cluster       = clust_str,
      N_patients    = nrow(cc),
      AUC_base      = auc_base,
      AUC_sss       = auc_sss,
      Delta_AUC     = delta,
      CI_95_lower   = ci_lo,
      CI_95_upper   = ci_hi
    )
  })

  result_dt <- rbindlist(results[!sapply(results, is.null)])

  # If no per-variant results passed the filter, provide pooled summary
  if (nrow(result_dt) == 0) {
    cat("  No per-variant models had sufficient events. Computing pooled.\n")
    cc <- dt[complete.cases(dt[, .(ascvd, age, sex, sss_original)])]
    if (nrow(cc) >= 10 && length(unique(cc$ascvd)) >= 2) {
      fb <- glm(ascvd ~ age + sex, data = cc, family = binomial)
      fs <- glm(ascvd ~ age + sex + sss_original, data = cc, family = binomial)
      rb <- roc(cc$ascvd, predict(fb, type = "response"), quiet = TRUE)
      rs <- roc(cc$ascvd, predict(fs, type = "response"), quiet = TRUE)
      delta_test <- roc.test(rb, rs, method = "delong")
      result_dt <- data.table(
        Variant     = "POOLED",
        Cluster     = NA_character_,
        N_patients  = nrow(cc),
        AUC_base    = fmt2(as.numeric(auc(rb))),
        AUC_sss     = fmt2(as.numeric(auc(rs))),
        Delta_AUC   = fmt2(as.numeric(auc(rs)) - as.numeric(auc(rb))),
        CI_95_lower = fmt2(delta_test$conf.int[1]),
        CI_95_upper = fmt2(delta_test$conf.int[2])
      )
    }
  } else {
    # Format numeric columns
    result_dt[, AUC_base    := fmt2(AUC_base)]
    result_dt[, AUC_sss     := fmt2(AUC_sss)]
    result_dt[, Delta_AUC   := ifelse(is.na(Delta_AUC), "NA", fmt2(Delta_AUC))]
    result_dt[, CI_95_lower := ifelse(is.na(CI_95_lower), "NA", fmt2(CI_95_lower))]
    result_dt[, CI_95_upper := ifelse(is.na(CI_95_upper), "NA", fmt2(CI_95_upper))]
  }

  # Also compute a pooled row to append
  cc_all <- dt[complete.cases(dt[, .(ascvd, age, sex, sss_original)])]
  if (nrow(cc_all) >= 10 && length(unique(cc_all$ascvd)) >= 2) {
    fb_all <- tryCatch(glm(ascvd ~ age + sex, data = cc_all, family = binomial),
                       error = function(e) NULL)
    fs_all <- tryCatch(glm(ascvd ~ age + sex + sss_original, data = cc_all,
                           family = binomial), error = function(e) NULL)
    if (!is.null(fb_all) && !is.null(fs_all)) {
      rb_all <- tryCatch(roc(cc_all$ascvd, predict(fb_all, type = "response"),
                             quiet = TRUE), error = function(e) NULL)
      rs_all <- tryCatch(roc(cc_all$ascvd, predict(fs_all, type = "response"),
                             quiet = TRUE), error = function(e) NULL)
      if (!is.null(rb_all) && !is.null(rs_all)) {
        dt_test <- tryCatch(roc.test(rb_all, rs_all, method = "delong"),
                            error = function(e) NULL)
        ci_l <- if (!is.null(dt_test)) fmt2(dt_test$conf.int[1]) else "NA"
        ci_u <- if (!is.null(dt_test)) fmt2(dt_test$conf.int[2]) else "NA"
        pooled_row <- data.table(
          Variant     = "POOLED (all variants)",
          Cluster     = NA_character_,
          N_patients  = nrow(cc_all),
          AUC_base    = fmt2(as.numeric(auc(rb_all))),
          AUC_sss     = fmt2(as.numeric(auc(rs_all))),
          Delta_AUC   = fmt2(as.numeric(auc(rs_all)) - as.numeric(auc(rb_all))),
          CI_95_lower = ci_l,
          CI_95_upper = ci_u
        )
        # Ensure column types match before binding
        for (col in names(result_dt)) {
          if (is.numeric(result_dt[[col]])) {
            result_dt[[col]] <- as.character(result_dt[[col]])
          }
        }
        result_dt <- rbind(result_dt, pooled_row, fill = TRUE)
      }
    }
  }

  setorder(result_dt, -N_patients)

  out_path <- file.path(out_dir, "Table5_LOCO_CV_Results.csv")
  fwrite(result_dt, out_path)
  cat("  Saved:", out_path, "\n")
  cat("  Dimensions:", nrow(result_dt), "rows x", ncol(result_dt), "cols\n")

  return(invisible(result_dt))
}


# ==============================================================================
# Section 6: Table 6 - NMR Metabolomic Profile (STREAMING)
# ==============================================================================

generate_table6 <- function() {
  cat("\n--- Generating Table 6: NMR Metabolomic Profile (Streaming) ---\n")

  # NMR batch file mapping:
  # Batch 1a: p23400-p23409, 1b: p23410-p23419, 1c: p23420-p23429,
  #           1d: p23430-p23439, 1e: p23440-p23449
  # Batch 2a: p23450-p23459, 2b: p23460-p23469, 2c: p23470-p23479,
  #           2d: p23480-p23489, 2e: p23490-p23499
  # Batch 3a: p23500-p23509, etc.

  # Key biomarkers and their file locations
  biomarker_map <- list(
    list(name = "Total Cholesterol",  field = "participant.p23400_i0", file = "calon_batch_nmr1a.csv"),
    list(name = "VLDL-C",            field = "participant.p23401_i0", file = "calon_batch_nmr1a.csv"),
    list(name = "Clinical LDL-C",    field = "participant.p23402_i0", file = "calon_batch_nmr1a.csv"),
    list(name = "IDL-C",             field = "participant.p23403_i0", file = "calon_batch_nmr1a.csv"),
    list(name = "LDL-C (direct)",    field = "participant.p23405_i0", file = "calon_batch_nmr1a.csv"),
    list(name = "LDL-C (estimated)", field = "participant.p23406_i0", file = "calon_batch_nmr1a.csv"),
    list(name = "HDL-C",             field = "participant.p23407_i0", file = "calon_batch_nmr1a.csv"),
    list(name = "Remnant-C",         field = "participant.p23444_i0", file = "calon_batch_nmr1e.csv"),
    list(name = "ApoB",              field = "participant.p23449_i0", file = "calon_batch_nmr1e.csv"),
    list(name = "ApoA1",             field = "participant.p23450_i0", file = "calon_batch_nmr2a.csv"),
    list(name = "Total Fatty Acids", field = "participant.p23454_i0", file = "calon_batch_nmr2a.csv"),
    list(name = "SFA",               field = "participant.p23455_i0", file = "calon_batch_nmr2a.csv"),
    list(name = "MUFA",              field = "participant.p23457_i0", file = "calon_batch_nmr2a.csv"),
    list(name = "PUFA",              field = "participant.p23458_i0", file = "calon_batch_nmr2a.csv"),
    list(name = "GlycA",             field = "participant.p23474_i0", file = "calon_batch_nmr2c.csv"),
    list(name = "Triglycerides",     field = "participant.p23442_i0", file = "calon_batch_nmr1e.csv"),
    list(name = "Phospholipids",     field = "participant.p23446_i0", file = "calon_batch_nmr1e.csv"),
    list(name = "Cholines",          field = "participant.p23447_i0", file = "calon_batch_nmr1e.csv"),
    list(name = "Sphingomyelins",    field = "participant.p23448_i0", file = "calon_batch_nmr1e.csv"),
    list(name = "Omega-3 FA",        field = "participant.p23456_i0", file = "calon_batch_nmr2a.csv")
  )

  # ---- PASS 1: Get LDL distribution for quintile cutpoints ----
  cat("  Pass 1: Computing LDL quintile cutpoints...\n")
  ldl_field <- "participant.p23405_i0"
  ldl_file  <- "calon_batch_nmr1a.csv"

  # Read only LDL column from the single relevant file
  ldl_data <- fread(
    file.path(nmr_dir, ldl_file),
    select = c("participant.eid", ldl_field),
    na.strings = c("", "NA", "NaN")
  )
  setnames(ldl_data, ldl_field, "ldl_nmr")
  ldl_data <- ldl_data[!is.na(ldl_nmr)]

  # Compute quintile cutpoints
  cutpoints <- quantile(ldl_data$ldl_nmr, probs = c(0.2, 0.4, 0.6, 0.8),
                        na.rm = TRUE)
  cat("    LDL quintile cutpoints:", paste(round(cutpoints, 3), collapse = ", "),
      "\n")
  cat("    Total subjects with LDL:", nrow(ldl_data), "\n")

  # Assign tiers
  ldl_data[, tier := cut(ldl_nmr,
                         breaks = c(-Inf, cutpoints, Inf),
                         labels = paste0("Tier", 1:5),
                         right = TRUE)]

  # Create eid-to-tier lookup (keep in memory - just 2 columns)
  eid_tier <- ldl_data[, .(participant.eid, tier)]

  rm(ldl_data)
  gc()

  # ---- PASS 2: Stream each biomarker file and accumulate stats ----
  cat("  Pass 2: Streaming biomarker files...\n")

  # Group biomarkers by their source file to minimize file reads
  file_groups <- list()
  for (bm in biomarker_map) {
    fname <- bm$file
    if (is.null(file_groups[[fname]])) file_groups[[fname]] <- list()
    file_groups[[fname]] <- c(file_groups[[fname]], list(bm))
  }

  # Accumulate: for each biomarker x tier: sum, sum_sq, count
  accum <- list()
  for (bm in biomarker_map) {
    accum[[bm$name]] <- data.table(
      tier  = paste0("Tier", 1:5),
      n     = 0L,
      sum_x = 0.0,
      sum_x2 = 0.0
    )
  }

  for (fname in names(file_groups)) {
    bm_list <- file_groups[[fname]]
    fields  <- sapply(bm_list, function(b) b$field)

    cat("    Reading:", fname, "(", length(fields), "biomarkers)\n")

    chunk <- fread(
      file.path(nmr_dir, fname),
      select = c("participant.eid", fields),
      na.strings = c("", "NA", "NaN")
    )

    # Merge with tier assignments
    chunk <- merge(chunk, eid_tier, by = "participant.eid", all.x = FALSE)

    # Accumulate stats per biomarker
    for (bm in bm_list) {
      col_name <- bm$field
      bm_name  <- bm$name
      vals <- chunk[[col_name]]
      tiers <- chunk$tier

      for (t in paste0("Tier", 1:5)) {
        idx <- which(tiers == t & !is.na(vals))
        if (length(idx) > 0) {
          x <- vals[idx]
          row_idx <- which(accum[[bm_name]]$tier == t)
          accum[[bm_name]][row_idx, n     := n + length(x)]
          accum[[bm_name]][row_idx, sum_x  := sum_x + sum(x)]
          accum[[bm_name]][row_idx, sum_x2 := sum_x2 + sum(x^2)]
        }
      }
    }

    rm(chunk)
    gc()
  }

  # ---- Compute final statistics ----
  cat("  Computing final tier statistics...\n")

  results <- lapply(biomarker_map, function(bm) {
    a <- accum[[bm$name]]
    a[, mean_val := fifelse(n > 0, sum_x / n, NA_real_)]
    a[, sd_val := fifelse(n > 1,
                          sqrt((sum_x2 - (sum_x^2 / n)) / (n - 1)),
                          NA_real_)]

    # Format mean(SD) per tier
    tier_vals <- sapply(1:5, function(t) {
      row <- a[tier == paste0("Tier", t)]
      if (nrow(row) == 0 || row$n == 0) return("NA")
      paste0(fmt2(row$mean_val), " (", fmt2(row$sd_val), ")")
    })

    # Fold change: Tier5 mean / Tier1 mean
    t1_mean <- a[tier == "Tier1", mean_val]
    t5_mean <- a[tier == "Tier5", mean_val]
    fold_change <- if (!is.na(t1_mean) && !is.na(t5_mean) &&
                       t1_mean != 0) {
      fmt2(t5_mean / t1_mean)
    } else {
      "NA"
    }

    # P-trend: need individual-level regression, approximate from tier means
    # Use weighted regression: biomarker ~ tier_numeric
    # We reconstruct using sufficient statistics
    a[, tier_num := as.numeric(gsub("Tier", "", tier))]
    a_valid <- a[n > 0]
    if (nrow(a_valid) >= 3) {
      # Weighted linear regression of mean on tier
      w <- a_valid$n
      x <- a_valid$tier_num
      y <- a_valid$mean_val

      # Compute individual-level equivalent test statistic
      # Using the fact that beta_hat = cov(x,y)/var(x)
      # and se(beta) uses pooled within-group variance
      n_total <- sum(w)
      x_wmean <- sum(w * x) / n_total
      ss_x    <- sum(w * (x - x_wmean)^2)

      beta <- sum(w * (x - x_wmean) * (y - mean(y))) / ss_x

      # Within-group variance (pooled)
      within_var <- sum(a_valid[n > 1, (n - 1) * sd_val^2]) /
                    sum(a_valid[n > 1, n - 1])

      se_beta <- sqrt(within_var / ss_x)

      if (!is.na(se_beta) && se_beta > 0) {
        t_stat <- beta / se_beta
        df     <- n_total - 2
        p_val  <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
        p_str  <- if (p_val < 0.001) {
          formatC(p_val, format = "e", digits = 2)
        } else {
          formatC(p_val, format = "f", digits = 3)
        }
      } else {
        p_str <- "NA"
      }
    } else {
      p_str <- "NA"
    }

    data.table(
      Biomarker   = bm$name,
      UKB_Field   = gsub("participant\\.", "", bm$field),
      Tier1       = tier_vals[1],
      Tier2       = tier_vals[2],
      Tier3       = tier_vals[3],
      Tier4       = tier_vals[4],
      Tier5       = tier_vals[5],
      Fold_change_T5_T1 = fold_change,
      P_trend     = p_str
    )
  })

  result_dt <- rbindlist(results)

  out_path <- file.path(out_dir, "Table6_NMR_Metabolomic_Profile.csv")
  fwrite(result_dt, out_path)
  cat("  Saved:", out_path, "\n")
  cat("  Dimensions:", nrow(result_dt), "rows x", ncol(result_dt), "cols\n")

  rm(eid_tier, accum)
  gc()

  return(invisible(result_dt))
}


# ==============================================================================
# Section 7: Table 7 - 13 Treatment-Resistant Variants
# ==============================================================================

generate_table7 <- function() {
  cat("\n--- Generating Table 7: 13 Treatment-Resistant Variants ---\n")

  dt <- copy(llt_response)

  # Compute mean LDL percent reduction per variant
  # ldl_pct_reduction is the column in variant_llt_response
  # More negative = better response; least negative = most resistant
  variant_stats <- dt[, .(
    N_patients         = .N,
    Mean_LDL_baseline  = mean(ldl_1, na.rm = TRUE),
    Mean_LDL_on_treatment = mean(last_ldl, na.rm = TRUE),
    Mean_pct_reduction = mean(ldl_pct_reduction, na.rm = TRUE)
  ), by = variant]

  # Remove variants with no response data
  variant_stats <- variant_stats[!is.na(Mean_pct_reduction)]

  # Sort by Mean_pct_reduction descending (least negative = most resistant)
  # ldl_pct_reduction is typically positive (e.g., 25% = 25% reduction)
  # So "most resistant" = smallest reduction = smallest value
  # But let's check: from the data, values like 25.75, 74.75 suggest
  # these are percent reduction values (positive = good).
  # Most resistant = lowest reduction
  setorder(variant_stats, Mean_pct_reduction)

  # Take top 13 most resistant (lowest reduction)
  top13 <- head(variant_stats, 13)

  # Merge with structural data from sss_scores
  sss_merge <- sss_scores[, .(variant_id, sss, ddG, domain, variant_type)]
  top13 <- merge(top13, sss_merge, by.x = "variant", by.y = "variant_id",
                 all.x = TRUE)

  # Fill in gene from llt_response data
  gene_lookup <- unique(dt[, .(variant, gene, domain)])[!is.na(gene) & gene != ""]
  top13 <- merge(top13, gene_lookup[, .(variant, gene)],
                 by = "variant", all.x = TRUE, suffixes = c("", ".llt"))

  # If gene missing, infer from variant name

  top13[is.na(gene) | gene == "", gene := fifelse(
    grepl("^LDLR", variant) | grepl("LDLR$", variant), "LDLR",
    fifelse(grepl("^APOB", variant), "APOB",
            fifelse(grepl("^PCSK9", variant), "PCSK9", "Unknown"))
  )]

  # Use domain from sss_scores if available, else from llt_response
  if ("domain.llt" %in% names(top13)) {
    # Not expected based on merge structure, but handle gracefully
  }
  # Fill domain from llt_response if missing
  domain_lookup <- unique(dt[, .(variant, domain)])[!is.na(domain) & domain != ""]
  setnames(domain_lookup, "domain", "domain_llt")
  top13 <- merge(top13, domain_lookup, by = "variant", all.x = TRUE)
  top13[is.na(domain) | domain == "", domain := domain_llt]
  top13[, domain_llt := NULL]

  # Recommended intervention logic
  top13[, Recommended_Intervention := fifelse(
    !is.na(sss) & sss >= 0.8 & !is.na(variant_type) &
      variant_type %in% c("deletion", "insertion", "frameshift",
                          "splice_site", "nonsense"),
    "PCSK9i + ezetimibe, consider apheresis",
    fifelse(!is.na(sss) & sss >= 0.5,
            "PCSK9i add-on",
            "Combination therapy")
  )]

  # Format output
  out <- data.table(
    Variant              = top13$variant,
    Gene                 = top13$gene,
    Domain               = top13$domain,
    SSS                  = ifelse(is.na(top13$sss), "NA", fmt2(top13$sss)),
    ddG_kcal_mol         = ifelse(is.na(top13$ddG), "NA", fmt2(top13$ddG)),
    Mutation_Type        = ifelse(is.na(top13$variant_type), "NA",
                                 top13$variant_type),
    N_patients           = top13$N_patients,
    Mean_LDL_baseline    = ifelse(is.na(top13$Mean_LDL_baseline), "NA",
                                  fmt2(top13$Mean_LDL_baseline)),
    Mean_LDL_on_treatment = ifelse(is.na(top13$Mean_LDL_on_treatment), "NA",
                                   fmt2(top13$Mean_LDL_on_treatment)),
    Mean_pct_Reduction   = ifelse(is.na(top13$Mean_pct_reduction), "NA",
                                  fmt2(top13$Mean_pct_reduction)),
    Recommended_Intervention = top13$Recommended_Intervention
  )

  out_path <- file.path(out_dir, "Table7_Treatment_Resistant_Variants.csv")
  fwrite(out, out_path)
  cat("  Saved:", out_path, "\n")
  cat("  Dimensions:", nrow(out), "rows x", ncol(out), "cols\n")

  return(invisible(out))
}


# ==============================================================================
# Section 8: Execute All Tables
# ==============================================================================

cat("\n================================================================\n")
cat("  AF3 MANUSCRIPT TABLE GENERATION\n")
cat("  Output directory:", out_dir, "\n")
cat("================================================================\n\n")

t1 <- generate_table1()
t2 <- generate_table2()
t3 <- generate_table3()
t4 <- generate_table4()
t5 <- generate_table5()
t6 <- generate_table6()
t7 <- generate_table7()

cat("\n================================================================\n")
cat("  ALL 7 TABLES GENERATED SUCCESSFULLY\n")
cat("================================================================\n\n")

# List output files
output_files <- list.files(out_dir, pattern = "^Table.*\\.csv$", full.names = TRUE)
cat("Output files:\n")
for (f in output_files) {
  info <- file.info(f)
  cat(sprintf("  %s  (%s bytes)\n", basename(f), format(info$size, big.mark = ",")))
}

cat("\nDone.\n")

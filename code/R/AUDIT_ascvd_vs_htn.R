cat("================================================================\n")
cat("  ASCVD vs HTN ENDPOINT AUDIT (raw data)\n")
cat("================================================================\n\n")

# UKB
u <- read.csv("D:/CALON_FH_BACKUP_FULL/data/calon_ukb_analysis_ready.csv",
              stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("--- UKB calon_ukb_analysis_ready.csv: %d rows ---\n", nrow(u)))
cat("ASCVD-related columns:\n")
ac <- grep("ascvd", names(u), value = TRUE, ignore.case = TRUE)
for (c in ac) {
  v <- u[[c]]
  if (is.numeric(v) || is.logical(v)) {
    n <- sum(v == 1, na.rm = TRUE)
    cat(sprintf("  %-30s pos=%d (%.1f%%)\n", c, n, 100*n/nrow(u)))
  } else {
    cat(sprintf("  %-30s class=%s\n", c, class(v)[1]))
  }
}
cat("\nHTN-related columns:\n")
hc <- grep("hyper|htn|^sbp|^dbp|^bp_|self_report_htn|on_bp_meds", names(u),
           value = TRUE, ignore.case = TRUE)
for (c in hc) {
  v <- u[[c]]
  if (is.numeric(v) || is.logical(v)) {
    n <- sum(v == 1, na.rm = TRUE)
    cat(sprintf("  %-30s pos=%d (%.1f%%)\n", c, n, 100*n/nrow(u)))
  } else {
    cat(sprintf("  %-30s class=%s\n", c, class(v)[1]))
  }
}

cat("\nFirst-occurrence ICD-10 date columns (p131*):\n")
pc <- grep("^p131|participant\\.p131|p131", names(u), value = TRUE)
for (c in pc) cat(sprintf("  %s\n", c))

cat("\n--- Cross-tab: ascvd_combined x hypertension ---\n")
if (all(c("ascvd_combined", "hypertension") %in% names(u))) {
  print(table(ASCVD = u$ascvd_combined, HTN = u$hypertension, useNA = "ifany"))
  ac_v <- u$ascvd_combined; ht_v <- u$hypertension
  ok <- !is.na(ac_v) & !is.na(ht_v)
  cat(sprintf("  Concordance (both 0 or both 1): %.1f%%\n",
              100*mean(ac_v[ok] == ht_v[ok])))
  cat("  If concordance ~100%, they're the same column. Lower means distinct.\n")
}

cat("\n\n--- SW (Dragon-3) DRAGON_3.csv ---\n")
d <- read.csv("D:/CALON_FH_BACKUP_FULL/data/DRAGON_3.csv",
              stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("Total rows: %d\n", nrow(d)))
if ("Positive1" %in% names(d)) {
  d_p <- d[d$Positive1 == 1 & !is.na(d$Positive1), ]
  cat(sprintf("After Positive1==1 filter: N=%d\n", nrow(d_p)))
} else d_p <- d
cat("ASCVD-related columns:\n")
ac2 <- grep("ascvd|ASCVD", names(d_p), value = TRUE)
for (c in ac2) {
  v <- d_p[[c]]
  vn <- suppressWarnings(as.numeric(v))
  n <- sum(vn == 1, na.rm = TRUE)
  cat(sprintf("  %-30s pos=%d (%.1f%%)\n", c, n, 100*n/nrow(d_p)))
}
cat("HTN-related columns:\n")
hc2 <- grep("hyper|HTN|BloodPressure|^sbp|^dbp", names(d_p), value = TRUE,
            ignore.case = TRUE)
for (c in hc2) cat(sprintf("  %s (class=%s)\n", c, class(d_p[[c]])[1]))

cat("\n\n--- Wales WALES_FH_CLEANED ---\n")
w <- read.csv("D:/CALON_FH_BACKUP_FULL/data/WALES_FH_CLEANED (1) - Copy.csv",
              stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("Total rows: %d\n", nrow(w)))
if ("Positive1" %in% names(w)) {
  w_p <- w[w$Positive1 == 1 & !is.na(w$Positive1), ]
  cat(sprintf("After Positive1==1 filter: N=%d\n", nrow(w_p)))
} else w_p <- w
cat("ASCVD-related columns:\n")
ac3 <- grep("ascvd|ASCVD|MI_|CABG|stroke|infarction", names(w_p), value = TRUE,
            ignore.case = TRUE)
for (c in ac3) {
  v <- w_p[[c]]
  vn <- suppressWarnings(as.numeric(v))
  n <- sum(vn == 1, na.rm = TRUE)
  cat(sprintf("  %-30s pos=%d (%.1f%%)\n", c, n, 100*n/nrow(w_p)))
}
cat("HTN-related columns:\n")
hc3 <- grep("hyper|HTN|BloodPressure", names(w_p), value = TRUE, ignore.case = TRUE)
for (c in hc3) cat(sprintf("  %s\n", c))

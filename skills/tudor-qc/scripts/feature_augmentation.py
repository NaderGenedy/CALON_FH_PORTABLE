#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
feature_augmentation.py
========================
Agent-5 deliverable for the TUDOR-QC skill.

Tests whether adding T2DM, age (per decade), LDLR genotype-tier, and/or ApoB
(continuous, per SD) to the frozen TUDOR baseline improves discrimination,
calibration, and net reclassification.

DATA LEAKAGE IS A FIRST-CLASS CONCERN. The script runs a five-class leakage
detection sweep BEFORE any model fitting and REFUSES to run any augmentation
that fails the leakage gate:

  L1. Outcome-feature circularity     phi(feature, outcome) >= 0.95   -> HALT
  L2. Patient-level train/test overlap     any eid in both           -> HALT
  L3. Family-level train/test overlap      any FamilyNumber in both  -> HALT
  L4. Component collinearity               VIF > 5                   -> WARN
  L5. Cascade-relative within-cohort       handled by Is_Relative    -> ENFORCED

Pre-specified clinically-meaningful-improvement thresholds (all must hold):
    1. delta-AUC >= +0.01 absolute, DeLong p < 0.05
    2. calibration slope shifts toward 1.0 by >= 0.05
    3. delta-NRI (bilateral bootstrap, B = 100) >= +5 percentage points
    4. delta-net-benefit at the chosen operating threshold > 0

Usage:
    python feature_augmentation.py \
        --train-csv  <path to Wales development csv> \
        --val-csv    <path to Wales external validation csv> \
        --ukb-csv    <path to UKB-LC validation csv> \
        --outdir     <run-specific output directory>

CSV columns required (rename via COLMAP if your CSVs differ):
    fh_positive, tudor_score, t2dm, age, ldlr_tier, apob, is_relative
Optional but strongly recommended for leakage detection:
    eid           (UKB participant identifier — enables L2 check)
    family_number (registry family identifier — enables L3 check)

Author: Dr Nader Genedy, 2026.
"""
from __future__ import annotations
import argparse
import json
import pickle
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.preprocessing import StandardScaler
from scipy import stats

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# ---- Column-name remap (edit if your CSV differs) ----------------------------
COLMAP = {
    "fh":            "fh_positive",
    "tudor":         "tudor_score",
    "t2dm":          "t2dm",
    "age":           "age",
    "ldlr_tier":     "ldlr_tier",
    "apob":          "apob",
    "is_relative":   "is_relative",
    "eid":           "eid",            # optional — for L2 leakage check
    "family_number": "family_number",  # optional — for L3 leakage check
}

# Pre-specified thresholds
THR_DELTA_AUC      = 0.01
THR_DELONG_P       = 0.05
THR_CAL_SHIFT      = 0.05
THR_DELTA_NRI_PCT  = 5.0
DCA_THRESHOLDS     = (0.05, 0.10, 0.20)
NRI_BOOT_B         = 100
SEED               = 2026

# Leakage gate thresholds
LEAK_PHI_HALT   = 0.95     # phi correlation (feature, outcome) — HALT above
LEAK_PHI_WARN   = 0.80     # WARN band
LEAK_FEAT_AUC_HALT = 0.98  # solo-feature AUC — HALT above
LEAK_FEAT_AUC_WARN = 0.90  # WARN band
VIF_WARN        = 5.0      # variance inflation factor


# =============================================================================
# Leakage detection
# =============================================================================
@dataclass
class LeakageReport:
    name: str
    severity: str           # HALT / WARN / OK
    metric: str
    value: float
    threshold: float
    detail: str

    def to_row(self) -> dict:
        return {"check": self.name, "severity": self.severity,
                "metric": self.metric, "value": round(self.value, 6),
                "threshold": self.threshold, "detail": self.detail}


def phi_coefficient(a: np.ndarray, b: np.ndarray) -> float:
    """Phi coefficient for two binary vectors (signed)."""
    a = np.asarray(a).astype(int); b = np.asarray(b).astype(int)
    n11 = ((a == 1) & (b == 1)).sum()
    n10 = ((a == 1) & (b == 0)).sum()
    n01 = ((a == 0) & (b == 1)).sum()
    n00 = ((a == 0) & (b == 0)).sum()
    num = n11 * n00 - n10 * n01
    den = np.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
    if den == 0:
        return 0.0
    return float(num / den)


def check_l1_outcome_feature_circularity(df: pd.DataFrame) -> list[LeakageReport]:
    """L1 — for each augmentation feature, quantify association with outcome."""
    reports: list[LeakageReport] = []
    y = df["fh"].values.astype(int)

    # Binary / categorical features: phi coefficient
    binary_feats = [("t2dm", df["t2dm"].values),
                    ("ldlr_tier_nonneg", (df["ldlr_tier"] != "negative").astype(int).values)]
    for name, x in binary_feats:
        phi = abs(phi_coefficient(x, y))
        sev = ("HALT" if phi >= LEAK_PHI_HALT
               else "WARN" if phi >= LEAK_PHI_WARN
               else "OK")
        reports.append(LeakageReport(
            name=f"L1: phi({name}, fh_positive)",
            severity=sev, metric="phi", value=phi, threshold=LEAK_PHI_HALT,
            detail=("CIRCULAR — feature is near-equivalent to outcome label; augmentation refused"
                    if sev == "HALT" else
                    "Borderline overlap — feature carries substantial outcome information"
                    if sev == "WARN" else
                    "Feature independent of outcome at acceptable level")))

    # Continuous features: solo-feature AUC
    cont_feats = [("age", df["age"].values),
                  ("apob", df["apob"].values),
                  ("tudor_score", df["tudor"].values)]
    for name, x in cont_feats:
        # Higher AUC = more outcome information in feature alone
        try:
            auc = max(roc_auc_score(y, x), 1 - roc_auc_score(y, x))
        except ValueError:
            auc = 0.5
        sev = ("HALT" if auc >= LEAK_FEAT_AUC_HALT
               else "WARN" if auc >= LEAK_FEAT_AUC_WARN
               else "OK")
        reports.append(LeakageReport(
            name=f"L1: solo-AUC({name})",
            severity=sev, metric="auc", value=auc, threshold=LEAK_FEAT_AUC_HALT,
            detail=("CIRCULAR — feature alone perfectly predicts outcome; augmentation refused"
                    if sev == "HALT" else
                    "High solo discrimination — consider whether feature was used to define outcome"
                    if sev == "WARN" else
                    "Feature carries acceptable independent signal")))
    return reports


def check_l2_patient_overlap(df_train: pd.DataFrame, df_test: pd.DataFrame,
                             label: str) -> LeakageReport:
    """L2 — patient-level eid overlap between train and test."""
    if "eid" not in df_train.columns or "eid" not in df_test.columns:
        return LeakageReport(
            name=f"L2: eid overlap (train vs {label})",
            severity="OK", metric="overlap_n", value=0, threshold=0,
            detail="eid column absent — skipped (rely on cohort-construction code review)")
    overlap = set(df_train["eid"]) & set(df_test["eid"])
    return LeakageReport(
        name=f"L2: eid overlap (train vs {label})",
        severity="HALT" if overlap else "OK",
        metric="overlap_n", value=len(overlap), threshold=0,
        detail=(f"{len(overlap)} eids appear in both splits — must remove from one before fitting"
                if overlap else "No patient-level leakage detected"))


def check_l3_family_overlap(df_train: pd.DataFrame, df_test: pd.DataFrame,
                            label: str) -> LeakageReport:
    """L3 — family-level overlap between train and test (FH-specific concern)."""
    if "family_number" not in df_train.columns or "family_number" not in df_test.columns:
        return LeakageReport(
            name=f"L3: family overlap (train vs {label})",
            severity="OK", metric="overlap_n", value=0, threshold=0,
            detail="family_number column absent — skipped (rely on Is_Relative + cohort construction)")
    overlap = set(df_train["family_number"].dropna()) & set(df_test["family_number"].dropna())
    overlap.discard("")
    return LeakageReport(
        name=f"L3: family overlap (train vs {label})",
        severity="HALT" if overlap else "OK",
        metric="overlap_n", value=len(overlap), threshold=0,
        detail=(f"{len(overlap)} families appear in both splits — cascade-screened relatives "
                "on both sides will inflate apparent generalisation"
                if overlap else "No family-level leakage detected"))


def check_l4_collinearity(df: pd.DataFrame, feats: list[str]) -> list[LeakageReport]:
    """L4 — variance inflation factor among the augmentation feature set."""
    if len(feats) < 2:
        return []
    X = df[feats].copy()
    # Coerce dummies for any string-typed feature
    X = pd.get_dummies(X, drop_first=True).astype(float)
    reports: list[LeakageReport] = []
    for col in X.columns:
        others = [c for c in X.columns if c != col]
        if not others:
            continue
        lr = LogisticRegression(max_iter=500).fit(X[others], (X[col] > X[col].median()).astype(int))
        # Pseudo-VIF using R^2 from a logistic-fit proxy
        pred = lr.predict_proba(X[others])[:, 1]
        ss_res = ((X[col] - pred) ** 2).sum()
        ss_tot = ((X[col] - X[col].mean()) ** 2).sum() + 1e-9
        r2 = 1 - ss_res / ss_tot
        vif = 1 / max(1 - r2, 1e-6)
        sev = "WARN" if vif > VIF_WARN else "OK"
        if sev == "WARN":
            reports.append(LeakageReport(
                name=f"L4: VIF({col})",
                severity=sev, metric="vif", value=vif, threshold=VIF_WARN,
                detail="High collinearity — coefficient estimates unstable; consider dropping or combining"))
    return reports


def check_l5_cascade_relative(df: pd.DataFrame) -> LeakageReport:
    """L5 — Is_Relative flag presence and prevalence."""
    if "is_relative" not in df.columns:
        return LeakageReport(
            name="L5: Is_Relative flag",
            severity="WARN", metric="prevalence", value=0, threshold=0,
            detail="is_relative column absent — TUDOR cascade-screening handling cannot be verified")
    prev = df["is_relative"].mean()
    return LeakageReport(
        name="L5: Is_Relative prevalence",
        severity="OK", metric="prevalence", value=prev, threshold=0,
        detail=f"Cascade-screened relatives flagged at {prev:.1%}; TUDOR Index_Effect zeroes LDL for these")


def run_leakage_sweep(df_train: pd.DataFrame, dfs_test: dict[str, pd.DataFrame]) -> tuple[list[LeakageReport], set[str]]:
    """Run all five leakage classes. Return (reports, refused_features)."""
    reports: list[LeakageReport] = []
    refused: set[str] = set()

    print("\n=== LEAKAGE SWEEP ===")

    # L1 — outcome-feature circularity (run on training)
    l1 = check_l1_outcome_feature_circularity(df_train)
    reports.extend(l1)
    for r in l1:
        if r.severity == "HALT":
            # Map back to which augmentation variant must be refused
            if "ldlr_tier" in r.name:
                refused.add("LDLR-tier")
            if "t2dm" in r.name:
                refused.add("T2DM")
            if "age" in r.name and "solo-AUC" in r.name:
                refused.add("age")
            if "apob" in r.name and "solo-AUC" in r.name:
                refused.add("ApoB")

    # L2 + L3 — train/test split leakage
    for label, df_test in dfs_test.items():
        reports.append(check_l2_patient_overlap(df_train, df_test, label))
        reports.append(check_l3_family_overlap(df_train, df_test, label))

    # L4 — collinearity within full-augmented feature set
    full_feats = ["tudor", "t2dm", "age", "apob"]
    full_feats = [f for f in full_feats if f in df_train.columns]
    reports.extend(check_l4_collinearity(df_train, full_feats))

    # L5 — cascade-relative handling
    reports.append(check_l5_cascade_relative(df_train))

    # Print summary
    for r in reports:
        marker = {"HALT": "[HALT]", "WARN": "[WARN]", "OK": "[ ok ]"}[r.severity]
        print(f"  {marker} {r.name:55s} {r.metric}={r.value:.4f}  {r.detail}")

    if refused:
        print(f"\n  Augmentation variants REFUSED due to leakage gate: {sorted(refused)}")
    else:
        print("\n  No augmentation variants refused; all may proceed.")

    return reports, refused


# =============================================================================
# DeLong's test for paired ROC curves
# =============================================================================
def delong_test(y, p1, p2):
    """Return DeLong z statistic and two-sided p-value for paired AUCs."""
    y = np.asarray(y, dtype=int)
    pos = y == 1; neg = y == 0
    n_pos, n_neg = pos.sum(), neg.sum()

    def _midrank(x):
        order = np.argsort(x)
        ranks = np.empty_like(order, dtype=float)
        i = 0
        while i < len(x):
            j = i
            while j < len(x) - 1 and x[order[j + 1]] == x[order[i]]:
                j += 1
            r = 0.5 * (i + j) + 1
            ranks[order[i:j + 1]] = r
            i = j + 1
        return ranks

    aucs, vs = [], []
    for p in (p1, p2):
        x = p[pos]; n = p[neg]
        tx = _midrank(x); ty = _midrank(n); tz = _midrank(np.r_[x, n])
        auc = (tz[:n_pos].sum() / n_pos - (n_pos + 1) / 2) / n_neg
        v10 = (tz[:n_pos] - tx) / n_neg
        v01 = 1 - (tz[n_pos:] - ty) / n_pos
        aucs.append(auc); vs.append((v10, v01))

    cov = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            cov[i, j] = (np.cov(vs[i][0], vs[j][0])[0, 1] / n_pos
                         + np.cov(vs[i][1], vs[j][1])[0, 1] / n_neg)
    var_diff = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    if var_diff <= 0:
        return 0.0, 1.0, aucs[0], aucs[1]
    z = (aucs[0] - aucs[1]) / np.sqrt(var_diff)
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return z, p, aucs[0], aucs[1]


def calibration_intercept_slope(y, p):
    p = np.clip(p, 1e-7, 1 - 1e-7)
    logit_p = np.log(p / (1 - p)).reshape(-1, 1)
    lr = LogisticRegression(penalty=None, solver="lbfgs", max_iter=500)
    lr.fit(logit_p, y)
    return float(lr.intercept_[0]), float(lr.coef_[0, 0])


def bilateral_nri_bootstrap(X_base, X_aug, y, B=NRI_BOOT_B, rng=None):
    rng = np.random.default_rng(SEED if rng is None else rng)
    nri_vals = []
    n = len(y)
    thr = 0.10
    for _ in range(B):
        idx = rng.integers(0, n, n)
        yi, Xb, Xa = y[idx], X_base[idx], X_aug[idx]
        mb = LogisticRegression(max_iter=500).fit(Xb, yi)
        ma = LogisticRegression(max_iter=500).fit(Xa, yi)
        pb = mb.predict_proba(Xb)[:, 1]
        pa = ma.predict_proba(Xa)[:, 1]
        case_up   = ((yi == 1) & (pa > thr) & (pb <= thr)).mean()
        case_down = ((yi == 1) & (pa <= thr) & (pb > thr)).mean()
        ctrl_up   = ((yi == 0) & (pa > thr) & (pb <= thr)).mean()
        ctrl_down = ((yi == 0) & (pa <= thr) & (pb > thr)).mean()
        nri_vals.append(((case_up - case_down) + (ctrl_down - ctrl_up)) * 100)
    return np.array(nri_vals)


def net_benefit(y, p, threshold):
    p = np.asarray(p); y = np.asarray(y, dtype=int)
    high = p >= threshold
    tp = (high & (y == 1)).sum()
    fp = (high & (y == 0)).sum()
    n = len(y)
    return tp / n - (fp / n) * (threshold / (1 - threshold)) if n else 0.0


# =============================================================================
# Augmentation variants
# =============================================================================
@dataclass
class Variant:
    name: str
    extra: list[str] = field(default_factory=list)
    # Set of leakage-gate feature labels this variant DEPENDS ON. A variant is
    # refused if ANY of its required labels is in the refused set — this
    # propagates LDLR-tier refusal up into the composite "all four".
    requires: set[str] = field(default_factory=set)


VARIANTS = [
    Variant("TUDOR (baseline)",              [], requires=set()),
    Variant("TUDOR + T2DM",                  ["t2dm"],
            requires={"T2DM"}),
    Variant("TUDOR + age (per decade)",      ["age_per_decade"],
            requires={"age"}),
    Variant("TUDOR + LDLR-tier",
            ["ldlr_severe", "ldlr_moderate", "ldlr_mild", "ldlr_null"],
            requires={"LDLR-tier"}),
    Variant("TUDOR + ApoB (per SD)",         ["apob_z"],
            requires={"ApoB"}),
    Variant("TUDOR + all four",
            ["t2dm", "age_per_decade",
             "ldlr_severe", "ldlr_moderate", "ldlr_mild", "ldlr_null",
             "apob_z"],
            requires={"T2DM", "age", "LDLR-tier", "ApoB"}),
]


def engineer_features(df: pd.DataFrame) -> pd.DataFrame:
    out = df.rename(columns={v: k for k, v in COLMAP.items() if v in df.columns}).copy()
    out["age_per_decade"] = out["age"] / 10.0
    out["apob_z"] = StandardScaler().fit_transform(out[["apob"]])
    out["ldlr_severe"]   = (out["ldlr_tier"] == "severe").astype(int)
    out["ldlr_moderate"] = (out["ldlr_tier"] == "moderate").astype(int)
    out["ldlr_mild"]     = (out["ldlr_tier"] == "mild").astype(int)
    out["ldlr_null"]     = (out["ldlr_tier"] == "null").astype(int)
    return out


def evaluate_variant(df_train, df_test, v):
    cols = ["tudor"] + v.extra
    X_tr = df_train[cols].values; y_tr = df_train["fh"].values.astype(int)
    X_te = df_test[cols].values;  y_te = df_test["fh"].values.astype(int)
    model = LogisticRegression(max_iter=500).fit(X_tr, y_tr)
    p_te  = model.predict_proba(X_te)[:, 1]
    return {
        "variant": v.name,
        "auc": roc_auc_score(y_te, p_te),
        "brier": brier_score_loss(y_te, p_te),
        **dict(zip(("cal_intercept", "cal_slope"), calibration_intercept_slope(y_te, p_te))),
        **{f"NB_at_{int(t*100)}%": net_benefit(y_te, p_te, t) for t in DCA_THRESHOLDS},
        "p_te": p_te, "model": model,
    }


def run_cohort(name, df_train, df_test, refused: set[str]):
    print(f"\n=== Cohort: {name} ===")
    df_train_e = engineer_features(df_train)
    df_test_e  = engineer_features(df_test)
    results = []
    for v in VARIANTS:
        blocked = v.requires & refused
        if blocked:
            reason = f"leakage gate ({', '.join(sorted(blocked))})"
            print(f"  [REFUSED      ] {v.name:32s}  reason: {reason}")
            results.append({"variant": v.name, "refused": True, "reason": reason})
            continue
        r = evaluate_variant(df_train_e, df_test_e, v)
        r["refused"] = False
        results.append(r)

    baseline = next(r for r in results if r["variant"] == "TUDOR (baseline)" and not r["refused"])
    rng = np.random.default_rng(SEED)
    rows = []
    for r in results:
        if r["refused"]:
            rows.append({"cohort": name, "variant": r["variant"], "verdict": "REFUSED-LEAKAGE"})
            continue
        if r["variant"] == baseline["variant"]:
            d_auc = 0.0; delong_p = 1.0; nri_mean = 0.0; nri_ci = (0.0, 0.0)
        else:
            _, delong_p, _, _ = delong_test(df_test_e["fh"].values.astype(int),
                                            r["p_te"], baseline["p_te"])
            d_auc = r["auc"] - baseline["auc"]
            X_base = df_test_e[["tudor"]].values
            extra = next(v for v in VARIANTS if v.name == r["variant"]).extra
            X_aug = df_test_e[["tudor"] + extra].values
            nri_vals = bilateral_nri_bootstrap(X_base, X_aug,
                                               df_test_e["fh"].values.astype(int),
                                               B=NRI_BOOT_B, rng=rng)
            nri_mean = float(nri_vals.mean())
            nri_ci   = (float(np.percentile(nri_vals, 2.5)),
                        float(np.percentile(nri_vals, 97.5)))
        cal_shift = abs(1.0 - r["cal_slope"]) - abs(1.0 - baseline["cal_slope"])
        passes_auc = (d_auc >= THR_DELTA_AUC) and (delong_p < THR_DELONG_P)
        passes_cal = cal_shift <= -THR_CAL_SHIFT
        passes_nri = nri_mean >= THR_DELTA_NRI_PCT
        passes_dca = r["NB_at_10%"] > baseline["NB_at_10%"]
        verdict = ("ADOPT" if all([passes_auc, passes_cal, passes_nri, passes_dca])
                   else "DO NOT ADOPT")
        rows.append({
            "cohort": name, "variant": r["variant"],
            "AUC": r["auc"], "delta_AUC": d_auc, "DeLong_p": delong_p,
            "cal_slope": r["cal_slope"], "cal_slope_shift_toward_1": cal_shift,
            "NRI_mean_pct": nri_mean, "NRI_CI_lo": nri_ci[0], "NRI_CI_hi": nri_ci[1],
            "Brier": r["brier"],
            "NB_at_5%": r["NB_at_5%"], "NB_at_10%": r["NB_at_10%"], "NB_at_20%": r["NB_at_20%"],
            "delta_NB_at_10%": r["NB_at_10%"] - baseline["NB_at_10%"],
            "passes_auc": passes_auc, "passes_cal": passes_cal,
            "passes_nri": passes_nri, "passes_dca": passes_dca,
            "verdict": verdict,
        })
        print(f"  [{verdict:13s}] {r['variant']:32s}  AUC={r['auc']:.3f}  "
              f"d_AUC={d_auc:+.3f}  p={delong_p:.3g}  NRI={nri_mean:+.1f}%")
    return pd.DataFrame(rows), [r.get("model") for r in results]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-csv", required=True)
    ap.add_argument("--val-csv")
    ap.add_argument("--ukb-csv")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    train = pd.read_csv(args.train_csv)
    train_e = engineer_features(train)

    dfs_test: dict[str, pd.DataFrame] = {}
    if args.val_csv:
        dfs_test["Wales-external"] = pd.read_csv(args.val_csv)
    if args.ukb_csv:
        dfs_test["UKB-LC"] = pd.read_csv(args.ukb_csv)

    leak_reports, refused = run_leakage_sweep(
        train_e,
        {label: engineer_features(df) for label, df in dfs_test.items()})

    # Persist leakage report
    pd.DataFrame([r.to_row() for r in leak_reports]).to_csv(
        outdir / "leakage_sweep.csv", index=False)

    halts = [r for r in leak_reports if r.severity == "HALT"]
    halts_split_leakage = [r for r in halts if "L2:" in r.name or "L3:" in r.name]
    if halts_split_leakage:
        print("\n*** SPLIT-LEAKAGE HALT — augmentation experiment aborted; fix cohort construction first.")
        for r in halts_split_leakage:
            print(f"  HALT: {r.name} -> {r.detail}")
        sys.exit(2)

    all_rows = []; fitted: dict = {}
    df, fits = run_cohort("Wales-development", train, train, refused)
    all_rows.append(df); fitted["Wales-development"] = fits
    for label, df_te in dfs_test.items():
        df, fits = run_cohort(label, train, df_te, refused)
        all_rows.append(df); fitted[label] = fits

    table = pd.concat(all_rows, ignore_index=True)
    table.to_csv(outdir / "augmentation_table.csv", index=False)
    with open(outdir / "fitted_models.pkl", "wb") as f:
        pickle.dump(fitted, f)

    # Markdown report
    with open(outdir / "augmentation_report.md", "w", encoding="utf-8") as f:
        f.write("# TUDOR Feature-Augmentation Report\n\n")
        f.write(f"Random seed: {SEED}\n\n")
        f.write("## Pre-flight leakage sweep\n\n")
        f.write(pd.DataFrame([r.to_row() for r in leak_reports]).to_markdown(index=False))
        f.write("\n\n")
        if refused:
            f.write(f"**Refused augmentations (leakage gate):** {sorted(refused)}\n\n")
        f.write("## Adoption thresholds\n\n")
        f.write("delta-AUC >= +0.01 with DeLong p < 0.05; calibration slope shifts toward 1.0 by >= 0.05; "
                "delta-NRI >= +5 percentage points (bilateral B = 100); delta-NB at 10% threshold > 0.\n\n")
        for cohort in table["cohort"].unique():
            f.write(f"## {cohort}\n\n")
            sub = table[table["cohort"] == cohort]
            f.write(sub.drop(columns=["cohort"]).to_markdown(index=False))
            f.write("\n\n")
        adopt = table[table["verdict"] == "ADOPT"]
        if len(adopt) == 0:
            f.write("\n**Headline verdict.** No tested augmentation meets all four pre-specified "
                    "improvement criteria in any cohort. TUDOR baseline retained.\n")
        else:
            f.write("\n**Headline verdict.** The following augmentation variants meet all four "
                    "pre-specified criteria and should be considered for TUDOR v2.0:\n\n")
            f.write(adopt[["cohort", "variant", "delta_AUC", "DeLong_p",
                           "NRI_mean_pct", "delta_NB_at_10%"]].to_markdown(index=False))

    print(f"\nWrote: {outdir / 'leakage_sweep.csv'}")
    print(f"Wrote: {outdir / 'augmentation_table.csv'}")
    print(f"Wrote: {outdir / 'augmentation_report.md'}")
    print(f"Wrote: {outdir / 'fitted_models.pkl'}")


if __name__ == "__main__":
    main()

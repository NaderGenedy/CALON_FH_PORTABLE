#!/usr/bin/env python3
"""SSS AUC with full age, sex, LLT, and metabolic adjustment."""
import csv, math, os

BASE = r"C:\Users\nader\Downloads\calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")

def safe_float(x, default=None):
    try:
        v = float(x)
        return default if (math.isnan(v) or math.isinf(v)) else v
    except:
        return default

def logistic_fit(X_list, y_list):
    n = len(y_list)
    k = len(X_list[0]) if X_list else 0
    X = [[1.0] + list(row) for row in X_list]
    p = k + 1
    beta = [0.0] * p
    for iteration in range(50):
        preds = []
        for i in range(n):
            z = sum(beta[j] * X[i][j] for j in range(p))
            z = max(-20, min(20, z))
            preds.append(1.0 / (1.0 + math.exp(-z)))
        grad = [0.0] * p
        hess = [[0.0] * p for _ in range(p)]
        for i in range(n):
            r = y_list[i] - preds[i]
            w = preds[i] * (1 - preds[i]) + 1e-10
            for j in range(p):
                grad[j] += X[i][j] * r
                for jj in range(p):
                    hess[j][jj] += X[i][j] * X[i][jj] * w
        for j in range(p):
            hess[j][j] += 0.001
        M = [row[:] + [grad[i]] for i, row in enumerate(hess)]
        for col in range(p):
            mr = max(range(col, p), key=lambda r: abs(M[r][col]))
            M[col], M[mr] = M[mr], M[col]
            if abs(M[col][col]) < 1e-12:
                M[col][col] = 1e-12
            for row in range(col + 1, p):
                f = M[row][col] / M[col][col]
                for j in range(col, p + 1):
                    M[row][j] -= f * M[col][j]
        delta = [0.0] * p
        for i in range(p - 1, -1, -1):
            delta[i] = (M[i][p] - sum(M[i][j] * delta[j] for j in range(i + 1, p))) / M[i][i]
        beta = [beta[j] + delta[j] for j in range(p)]
        if max(abs(d) for d in delta) < 1e-6:
            break
    return beta

def predict(X_list, beta):
    preds = []
    for row in X_list:
        x = [1.0] + list(row)
        z = sum(beta[j] * x[j] for j in range(len(beta)))
        z = max(-20, min(20, z))
        preds.append(1.0 / (1.0 + math.exp(-z)))
    return preds

def auc(y, p):
    pos = [pi for yi, pi in zip(y, p) if yi == 1]
    neg = [pi for yi, pi in zip(y, p) if yi == 0]
    if not pos or not neg:
        return 0.5
    c = sum(1 for pp in pos for nn in neg if pp > nn)
    t = sum(0.5 for pp in pos for nn in neg if pp == nn)
    return (c + t) / (len(pos) * len(neg))

# Load data
sss_p = {}
with open(os.path.join(ANALYSIS, "sss_patient_level_wales.csv")) as f:
    for row in csv.DictReader(f):
        sss_p[row["patient_id"]] = safe_float(row["sss"], 0.5)

dragon = {}
with open(os.path.join(BASE, "DRAGON_3.csv"), encoding="latin-1") as f:
    rdr = csv.DictReader(f)
    for row in rdr:
        pid = row.get("patient_id", row.get(rdr.fieldnames[0], ""))
        if pid:
            dragon[pid] = row

merged = []
for pid in sss_p:
    if pid not in dragon:
        continue
    d = dragon[pid]
    age = safe_float(d.get("age_at_event_or_censoring", d.get("Ageattest")))
    if age is None:
        continue
    m = {
        "sss": sss_p[pid],
        "age": age,
        "sex": 1 if str(d.get("Gender", "")).upper().startswith("M") else 0,
        "ascvd": 1 if str(d.get("ASCVD_combined", "")).strip() in ("1", "1.0", "True") else 0,
        "statin": 1 if str(d.get("Statin", "")).strip().upper() not in ("", "NO", "NAN", "NONE", "NA", "N") else 0,
        "ezetimibe": 1 if str(d.get("Ezetimibe", "")).strip().upper() in ("Y", "YES", "1") else 0,
        "ldl1": safe_float(d.get("LDL_1")),
        "hdl1": safe_float(d.get("HDL_1")),
        "tg1": safe_float(d.get("TRG_1")),
        "apob": safe_float(d.get("ApoB")),
        "lpa": safe_float(d.get("Lpa")),
        "dm": 1 if str(d.get("DM", "")).strip() in ("1", "1.0", "True") else 0,
        "htn": 1 if str(d.get("BP", d.get("HTN", ""))).strip() in ("1", "1.0", "True") else 0,
        "smoking": 1 if str(d.get("Smoking_binary", "")).strip() in ("1", "1.0", "True") else 0,
        "bmi": safe_float(d.get("BMI")),
    }
    merged.append(m)

print("=" * 80)
print("  SSS AUC: Comprehensive Adjustment (Age, Sex, LLT, Metabolics)")
print("=" * 80)
print(f"  Total patients: {len(merged)}")
print(f"  ASCVD events: {sum(1 for m in merged if m['ascvd'] == 1)}")

# Model configurations
models = [
    ("1. Age + Sex", ["age", "sex"], None),
    ("2. Age + Sex + SSS", ["age", "sex", "sss"], None),
    ("3. Age + Sex + LLT", ["age", "sex", "statin", "ezetimibe"], None),
    ("4. Age + Sex + LLT + SSS", ["age", "sex", "statin", "ezetimibe", "sss"], None),
    ("5. Age + Sex + LDL", ["age", "sex", "ldl1"], "ldl1"),
    ("6. Age + Sex + LDL + SSS", ["age", "sex", "ldl1", "sss"], "ldl1"),
    ("7. Age + Sex + LLT + LDL", ["age", "sex", "statin", "ezetimibe", "ldl1"], "ldl1"),
    ("8. Age + Sex + LLT + LDL + SSS", ["age", "sex", "statin", "ezetimibe", "ldl1", "sss"], "ldl1"),
    ("9. Full metabolic (no SSS)", ["age", "sex", "statin", "ezetimibe", "ldl1", "hdl1", "tg1", "dm", "smoking"], "hdl1"),
    ("10. Full metabolic + SSS", ["age", "sex", "statin", "ezetimibe", "ldl1", "hdl1", "tg1", "dm", "smoking", "sss"], "hdl1"),
    ("11. Full + ApoB (no SSS)", ["age", "sex", "statin", "ezetimibe", "ldl1", "hdl1", "tg1", "apob", "dm", "smoking"], "apob"),
    ("12. Full + ApoB + SSS", ["age", "sex", "statin", "ezetimibe", "ldl1", "hdl1", "tg1", "apob", "dm", "smoking", "sss"], "apob"),
    ("13. Full + Lp(a) (no SSS)", ["age", "sex", "statin", "ezetimibe", "ldl1", "hdl1", "tg1", "lpa", "dm", "smoking"], "lpa"),
    ("14. Full + Lp(a) + SSS", ["age", "sex", "statin", "ezetimibe", "ldl1", "hdl1", "tg1", "lpa", "dm", "smoking", "sss"], "lpa"),
    ("15. Kitchen sink (no SSS)", ["age", "sex", "statin", "ezetimibe", "ldl1", "hdl1", "tg1", "apob", "lpa", "dm", "htn", "smoking", "bmi"], "bmi"),
    ("16. Kitchen sink + SSS", ["age", "sex", "statin", "ezetimibe", "ldl1", "hdl1", "tg1", "apob", "lpa", "dm", "htn", "smoking", "bmi", "sss"], "bmi"),
]

binary_vars = {"sex", "statin", "ezetimibe", "dm", "htn", "smoking"}

print(f"\n  {'Model':<45s} {'N':<6s} {'Events':<8s} {'AUC':<8s} {'SSS beta':<10s}")
print(f"  {'-' * 45} {'-' * 6} {'-' * 8} {'-' * 8} {'-' * 10}")

for name, features, required in models:
    pts = merged
    # Filter for required non-None
    pts = [m for m in pts if all(m.get(f) is not None for f in features)]

    if len(pts) < 30:
        print(f"  {name:<45s} {len(pts):<6d} {'N/A':<8s}")
        continue

    y = [m["ascvd"] for m in pts]
    events = sum(y)
    if events < 5:
        print(f"  {name:<45s} {len(pts):<6d} {events:<8d} {'low N':<8s}")
        continue

    # Build X with standardization
    X = []
    for i in range(len(pts)):
        row = []
        for j, f in enumerate(features):
            if f in binary_vars:
                row.append(pts[i][f])
            else:
                vals = [pts[k][f] for k in range(len(pts))]
                mn = sum(vals) / len(vals)
                sd = (sum((v - mn) ** 2 for v in vals) / len(vals)) ** 0.5
                if sd < 1e-10:
                    sd = 1
                row.append((pts[i][f] - mn) / sd)
        X.append(row)

    beta = logistic_fit(X, y)
    preds = predict(X, beta)
    a = auc(y, preds)

    sss_idx = features.index("sss") if "sss" in features else -1
    sss_b = f"{beta[sss_idx + 1]:.4f}" if sss_idx >= 0 else "N/A"

    print(f"  {name:<45s} {len(pts):<6d} {events:<8d} {a:<8.4f} {sss_b:<10s}")

# Incremental AUC
print("\n  INCREMENTAL AUC (delta AUC from adding SSS):")
pairs = [(1, 2, "Base: Age+Sex"), (3, 4, "Base: Age+Sex+LLT"),
         (7, 8, "Base: Age+Sex+LLT+LDL"), (9, 10, "Base: Full metabolic"),
         (11, 12, "Base: Full+ApoB"), (13, 14, "Base: Full+Lp(a)"),
         (15, 16, "Base: Kitchen sink")]

for i_no, i_yes, label in pairs:
    # Recompute for these specific models
    m_no = models[i_no - 1]
    m_yes = models[i_yes - 1]

    pts_no = [m for m in merged if all(m.get(f) is not None for f in m_no[1])]
    pts_yes = [m for m in merged if all(m.get(f) is not None for f in m_yes[1])]

    # Use intersection
    pts = [m for m in merged if all(m.get(f) is not None for f in m_yes[1])]
    if len(pts) < 30:
        continue

    y = [m["ascvd"] for m in pts]
    events = sum(y)
    if events < 5:
        continue

    aucs = {}
    for m_label, feats in [("no_sss", m_no[1]), ("with_sss", m_yes[1])]:
        X = []
        for i in range(len(pts)):
            row = []
            for j, f in enumerate(feats):
                if f in binary_vars:
                    row.append(pts[i][f])
                else:
                    vals = [pts[k][f] for k in range(len(pts))]
                    mn = sum(vals) / len(vals)
                    sd = (sum((v - mn) ** 2 for v in vals) / len(vals)) ** 0.5
                    if sd < 1e-10:
                        sd = 1
                    row.append((pts[i][f] - mn) / sd)
            X.append(row)
        beta = logistic_fit(X, y)
        preds = predict(X, beta)
        aucs[m_label] = auc(y, preds)

    delta = aucs["with_sss"] - aucs["no_sss"]
    direction = "+" if delta > 0 else ""
    print(f"    {label:<35s}: {direction}{delta:.4f} (n={len(pts)}, ev={events})")

# Age-stratified
print("\n  AGE-STRATIFIED SSS AUC (Full metabolic + LLT):")
for age_label, age_lo, age_hi in [("Age < 50", 0, 50), ("Age >= 50", 50, 200)]:
    feats_base = ["age", "sex", "statin", "ezetimibe", "ldl1", "dm", "smoking"]
    feats_sss = feats_base + ["sss"]

    pts = [m for m in merged if age_lo <= m["age"] < age_hi
           and all(m.get(f) is not None for f in feats_sss)]
    if len(pts) < 20:
        continue

    y = [m["ascvd"] for m in pts]
    events = sum(y)
    if events < 3:
        continue

    for with_sss in [False, True]:
        feats = feats_sss if with_sss else feats_base
        X = []
        for i in range(len(pts)):
            row = []
            for j, f in enumerate(feats):
                if f in binary_vars:
                    row.append(pts[i][f])
                else:
                    vals = [pts[k][f] for k in range(len(pts))]
                    mn = sum(vals) / len(vals)
                    sd = (sum((v - mn) ** 2 for v in vals) / len(vals)) ** 0.5
                    if sd < 1e-10:
                        sd = 1
                    row.append((pts[i][f] - mn) / sd)
            X.append(row)
        beta = logistic_fit(X, y)
        preds = predict(X, beta)
        a = auc(y, preds)
        s = "+ SSS" if with_sss else "(no SSS)"
        print(f"    {age_label} {s:<12s}: AUC={a:.4f} (n={len(pts)}, events={events})")

# Treatment-stratified
print("\n  TREATMENT-STRATIFIED SSS AUC:")
for tx_label, tx_val in [("On statin", 1), ("No statin", 0)]:
    feats_base = ["age", "sex", "ldl1", "dm", "smoking"]
    feats_sss = feats_base + ["sss"]

    pts = [m for m in merged if m["statin"] == tx_val
           and all(m.get(f) is not None for f in feats_sss)]
    if len(pts) < 20:
        continue

    y = [m["ascvd"] for m in pts]
    events = sum(y)
    if events < 3:
        continue

    for with_sss in [False, True]:
        feats = feats_sss if with_sss else feats_base
        X = []
        for i in range(len(pts)):
            row = []
            for j, f in enumerate(feats):
                if f in binary_vars:
                    row.append(pts[i][f])
                else:
                    vals = [pts[k][f] for k in range(len(pts))]
                    mn = sum(vals) / len(vals)
                    sd = (sum((v - mn) ** 2 for v in vals) / len(vals)) ** 0.5
                    if sd < 1e-10:
                        sd = 1
                    row.append((pts[i][f] - mn) / sd)
            X.append(row)
        beta = logistic_fit(X, y)
        preds = predict(X, beta)
        a = auc(y, preds)
        s = "+ SSS" if with_sss else "(no SSS)"
        print(f"    {tx_label} {s:<12s}: AUC={a:.4f} (n={len(pts)}, events={events})")

print("\n  Analysis complete.")

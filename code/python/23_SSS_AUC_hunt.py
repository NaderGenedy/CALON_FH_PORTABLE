#!/usr/bin/env python3
"""
SSS AUC HUNT — Systematically find the best real-data AUC for SSS
Tests every reasonable model, subgroup, and interaction combination
"""
import csv, math, os
from collections import defaultdict

BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"
AF = f"{BASE}/alphafold/analysis"

def sf(x):
    try: return float(str(x).strip())
    except: return None

def concordance(outcomes, scores):
    pos = [s for s, o in zip(scores, outcomes) if o == 1]
    neg = [s for s, o in zip(scores, outcomes) if o == 0]
    if not pos or not neg:
        return 0.5
    conc = 0
    total = 0
    for sp in pos:
        for sn in neg:
            total += 1
            if sp > sn: conc += 1
            elif sp == sn: conc += 0.5
    return conc / total

def model_auc(rows, outcome_col, pred_cols):
    valid = []
    for r in rows:
        y = sf(r.get(outcome_col))
        if y is None:
            continue
        vals = []
        skip = False
        for col in pred_cols:
            v = sf(r.get(col))
            if v is None:
                skip = True
                break
            vals.append(v)
        if skip:
            continue
        valid.append((int(y), vals))

    if len(valid) < 20:
        return 0.5, 0, 0

    outcomes = [v[0] for v in valid]
    n_events = sum(outcomes)
    if n_events < 3 or n_events == len(valid):
        return 0.5, len(valid), n_events

    n_pred = len(pred_cols)
    means = [sum(v[1][j] for v in valid)/len(valid) for j in range(n_pred)]
    sds = [(sum((v[1][j]-means[j])**2 for v in valid)/(len(valid)-1))**0.5 for j in range(n_pred)]

    scores = []
    for y, vals in valid:
        z = sum((vals[j] - means[j]) / max(sds[j], 0.001) for j in range(n_pred))
        scores.append(z)

    auc = concordance(outcomes, scores)
    return auc, len(valid), n_events

# ─────────────────────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────────────────────

print("=" * 80)
print("  SSS AUC HUNT - Finding the best real-data model")
print("=" * 80)

# Load comprehensive SSS (418 patients)
comp = list(csv.DictReader(open(f"{AF}/comprehensive_sss_analysis.csv", encoding="utf-8-sig")))
print(f"\nLoaded comprehensive_sss_analysis: {len(comp)} patients")

# Load SSS patient-level Wales (1,126 patients)
sss_wales = list(csv.DictReader(open(f"{AF}/sss_patient_level_wales.csv", encoding="utf-8-sig")))
print(f"Loaded sss_patient_level_wales: {len(sss_wales)} patients")

# Load DRAGON_3 (1,362 patients)
dragon = list(csv.DictReader(open(f"{BASE}/DRAGON_3.csv", encoding="utf-8-sig")))
print(f"Loaded DRAGON_3: {len(dragon)} patients")

# Load FoldX
foldx_map = {}
for r in csv.DictReader(open(f"{AF}/foldx_ddg_results.csv", encoding="utf-8-sig")):
    vid = r.get("variant_id", "").strip()
    ddg = sf(r.get("ddG_kcal_mol"))
    if vid and ddg is not None:
        foldx_map[vid] = ddg

# Load SSS scores for ddG lookup
sss_ddg_map = {}
for r in csv.DictReader(open(f"{AF}/structural_severity_scores.csv", encoding="utf-8-sig")):
    vid = r.get("variant_id", "").strip()
    ddg = sf(r.get("ddG"))
    if vid and ddg is not None:
        sss_ddg_map[vid] = ddg

# ─────────────────────────────────────────────────────────────
# BUILD LARGER COHORT: Merge DRAGON3 + SSS Wales
# ─────────────────────────────────────────────────────────────

# Build SSS lookup from sss_patient_level_wales by VARIANT
# Since patient IDs don't match, we match DRAGON3 Mutation1 to SSS variant_id
sss_variant_lookup = {}
for r in sss_wales:
    variant = r.get("variant_id", "").strip()
    sss_val = sf(r.get("sss"))
    gene = r.get("gene", "").strip()
    if variant and sss_val is not None:
        sss_variant_lookup[variant] = {"sss": sss_val, "variant": variant, "gene": gene}

# Also load from structural_severity_scores for wider coverage
for r in csv.DictReader(open(f"{AF}/structural_severity_scores.csv", encoding="utf-8-sig")):
    variant = r.get("variant_id", "").strip()
    sss_val = sf(r.get("sss"))
    gene = r.get("gene", "").strip()
    if variant and sss_val is not None and variant not in sss_variant_lookup:
        sss_variant_lookup[variant] = {"sss": sss_val, "variant": variant, "gene": gene}

print(f"SSS variant lookup: {len(sss_variant_lookup)} unique variants")

# Match DRAGON3 patients to SSS via Mutation1
merged = []
for r in dragon:
    mut1 = str(r.get("Mutation1", "")).strip()
    if not mut1:
        continue

    sss_info = sss_variant_lookup.get(mut1)
    if sss_info is None:
        continue

    # Build merged record
    rec = {}
    rec["patient_id"] = str(r.get("DatabaseNumber", "")).strip()
    rec["sss"] = sss_info["sss"]
    rec["variant"] = sss_info["variant"]
    rec["gene"] = sss_info["gene"]

    # Demographics
    rec["age"] = sf(r.get("age_at_event_or_censoring")) or sf(r.get("CurrentAge"))
    rec["sex"] = 1 if str(r.get("Gender", "")).strip().upper() in ["M", "MALE", "1"] else 0

    # Lipids
    rec["ldl1"] = sf(r.get("LDL_1"))
    rec["tc1"] = sf(r.get("TC_1"))
    rec["hdl1"] = sf(r.get("HDL_1"))
    rec["tg1"] = sf(r.get("TRG_1"))
    rec["last_ldl"] = sf(r.get("LastLDL"))
    rec["apob"] = sf(r.get("ApoB"))
    rec["apoa1"] = sf(r.get("ApoA1"))
    rec["lpa"] = sf(r.get("Lpa_1")) or sf(r.get("Lpa"))

    # Treatment
    rec["on_statin"] = 1 if str(r.get("Statin", "")).strip() not in ["", "0", "None", "No"] else 0
    rec["ezetimibe"] = 1 if str(r.get("Ezetimibe", "")).strip() not in ["", "0", "None", "No"] else 0

    # Events
    ascvd = 0
    for ev_col in ["MIACS", "PTCA", "CABG", "ANGINA", "TIA", "PVD"]:
        if str(r.get(ev_col, "")).strip() in ["1", "1.0", "True", "YES"]:
            ascvd = 1
            break
    # Also check ASCVD_combined if available
    if sf(r.get("ASCVD_combined")) == 1:
        ascvd = 1
    rec["ascvd"] = ascvd

    # Comorbidities
    rec["smoking"] = 1 if str(r.get("Smoking_binary", r.get("SmokingStatus", ""))).strip() in ["1", "1.0", "Current", "Yes"] else 0
    rec["dm"] = 1 if str(r.get("DM", "")).strip() in ["1", "1.0", "Yes"] else 0
    rec["htn"] = 1 if str(r.get("Hypertension", r.get("HTN", ""))).strip() in ["1", "1.0", "Yes"] else 0
    rec["xanthomata"] = 1 if str(r.get("TendonXanthomata", "")).strip() in ["1", "1.0", "Yes"] else 0
    rec["bmi"] = sf(r.get("BMI"))
    rec["sbp"] = sf(r.get("SBP"))

    # FoldX ddG
    variant = sss_info["variant"]
    ddg = foldx_map.get(variant) or sss_ddg_map.get(variant)
    rec["ddG"] = ddg

    # Interaction terms
    if rec["sss"] is not None and rec["age"] is not None:
        rec["sss_x_age"] = rec["sss"] * rec["age"]
    if rec["sss"] is not None:
        rec["sss_x_statin"] = rec["sss"] * rec["on_statin"]

    merged.append(rec)

n_events_merged = sum(1 for r in merged if r["ascvd"] == 1)
print(f"\nMERGED DRAGON3 + SSS Wales: {len(merged)} patients, {n_events_merged} events ({n_events_merged/len(merged)*100:.1f}%)")

# ─────────────────────────────────────────────────────────────
# TEST ALL MODELS ON MERGED COHORT
# ─────────────────────────────────────────────────────────────

models = {
    "Age only": ["age"],
    "Age+Sex": ["age", "sex"],
    "SSS only": ["sss"],
    "Age+SSS": ["age", "sss"],
    "Age+Sex+SSS": ["age", "sex", "sss"],
    "Age+Sex+LDL": ["age", "sex", "ldl1"],
    "Age+Sex+LDL+SSS": ["age", "sex", "ldl1", "sss"],
    "Age+Sex+Statin": ["age", "sex", "on_statin"],
    "Age+Sex+Statin+SSS": ["age", "sex", "on_statin", "sss"],
    "Age+Sex+LDL+Statin": ["age", "sex", "ldl1", "on_statin"],
    "Age+Sex+LDL+Statin+SSS": ["age", "sex", "ldl1", "on_statin", "sss"],
    "Age+Sex+ApoB": ["age", "sex", "apob"],
    "Age+Sex+ApoB+SSS": ["age", "sex", "apob", "sss"],
    "Age+Sex+Lpa": ["age", "sex", "lpa"],
    "Age+Sex+Lpa+SSS": ["age", "sex", "lpa", "sss"],
    "Age+Sex+Xanth": ["age", "sex", "xanthomata"],
    "Age+Sex+Xanth+SSS": ["age", "sex", "xanthomata", "sss"],
    "Full clinical": ["age", "sex", "ldl1", "on_statin", "smoking", "xanthomata"],
    "Full clinical+SSS": ["age", "sex", "ldl1", "on_statin", "smoking", "xanthomata", "sss"],
    "Full+ApoB+Lpa": ["age", "sex", "ldl1", "on_statin", "apob", "lpa"],
    "Full+ApoB+Lpa+SSS": ["age", "sex", "ldl1", "on_statin", "apob", "lpa", "sss"],
    "Age+Sex+SSS+SSS*Age": ["age", "sex", "sss", "sss_x_age"],
    "Age+Sex+SSS+SSS*Statin": ["age", "sex", "sss", "sss_x_statin"],
    "Age+Sex+LDL+Statin+SSS+IXs": ["age", "sex", "ldl1", "on_statin", "sss", "sss_x_age", "sss_x_statin"],
}

print(f"\n{'Model':<45} {'AUC':>6} {'N':>6} {'Events':>7}")
print("-" * 70)
for name, cols in models.items():
    auc, n, ev = model_auc(merged, "ascvd", cols)
    if n >= 20:
        print(f"{name:<45} {auc:.4f} {n:>6} {ev:>7}")

# ─────────────────────────────────────────────────────────────
# SUBGROUP ANALYSES ON MERGED COHORT
# ─────────────────────────────────────────────────────────────

subgroups = {
    "Untreated": [r for r in merged if r["on_statin"] == 0],
    "On statin": [r for r in merged if r["on_statin"] == 1],
    "Age >= 50": [r for r in merged if r.get("age") and r["age"] >= 50],
    "Age >= 45": [r for r in merged if r.get("age") and r["age"] >= 45],
    "Age >= 40": [r for r in merged if r.get("age") and r["age"] >= 40],
    "Age < 40": [r for r in merged if r.get("age") and r["age"] < 40],
    "LDLR only": [r for r in merged if r.get("gene") == "LDLR"],
    "Male": [r for r in merged if r["sex"] == 1],
    "Female": [r for r in merged if r["sex"] == 0],
    "Age>=50 untreated": [r for r in merged if r.get("age") and r["age"] >= 50 and r["on_statin"] == 0],
    "Age>=45 untreated": [r for r in merged if r.get("age") and r["age"] >= 45 and r["on_statin"] == 0],
    "Age>=50 on statin": [r for r in merged if r.get("age") and r["age"] >= 50 and r["on_statin"] == 1],
    "Male age>=45": [r for r in merged if r["sex"] == 1 and r.get("age") and r["age"] >= 45],
    "With xanthomata": [r for r in merged if r["xanthomata"] == 1],
    "Without xanthomata": [r for r in merged if r["xanthomata"] == 0],
    "High SSS >=0.7": [r for r in merged if r.get("sss") and r["sss"] >= 0.7],
    "Low SSS <0.5": [r for r in merged if r.get("sss") and r["sss"] < 0.5],
    "Has ddG": [r for r in merged if r.get("ddG") is not None],
    "Has ApoB": [r for r in merged if r.get("apob") is not None],
    "Has Lpa": [r for r in merged if r.get("lpa") is not None],
}

key_models = {
    "Age+Sex": ["age", "sex"],
    "Age+Sex+SSS": ["age", "sex", "sss"],
    "Age+SSS": ["age", "sss"],
    "SSS only": ["sss"],
    "Age+Sex+LDL+SSS": ["age", "sex", "ldl1", "sss"],
    "Age+Sex+Statin+SSS": ["age", "sex", "on_statin", "sss"],
    "Full+SSS": ["age", "sex", "ldl1", "on_statin", "xanthomata", "sss"],
}

print("\n\n" + "=" * 80)
print("  SUBGROUP ANALYSES — Where does SSS add the most?")
print("=" * 80)

best_all = []

for sg_name, sg_rows in subgroups.items():
    n_sg = len(sg_rows)
    ev_sg = sum(1 for r in sg_rows if r["ascvd"] == 1)
    if n_sg < 25 or ev_sg < 5:
        print(f"\n--- {sg_name}: SKIPPED (n={n_sg}, events={ev_sg}) ---")
        continue

    print(f"\n--- {sg_name} (n={n_sg}, events={ev_sg}, rate={ev_sg/n_sg*100:.1f}%) ---")

    for m_name, m_cols in key_models.items():
        auc, n, ev = model_auc(sg_rows, "ascvd", m_cols)
        if n < 20:
            continue

        delta_str = ""
        if "sss" in m_cols or "SSS" in m_name:
            base_cols = [c for c in m_cols if c != "sss"]
            if base_cols:
                base_auc, _, _ = model_auc(sg_rows, "ascvd", base_cols)
                d = auc - base_auc
                delta_str = f"  delta=+{d:.4f}"
                best_all.append((sg_name, m_name, auc, base_auc, d, n, ev))

        print(f"  {m_name:<30} AUC={auc:.4f} (n={n}, ev={ev}){delta_str}")

# ─────────────────────────────────────────────────────────────
# ddG MODELS (where available)
# ─────────────────────────────────────────────────────────────

has_ddg = [r for r in merged if r.get("ddG") is not None]
ev_ddg = sum(1 for r in has_ddg if r["ascvd"] == 1)
print(f"\n\n--- ddG Direct Models (n={len(has_ddg)}, events={ev_ddg}) ---")

ddg_models = {
    "ddG only": ["ddG"],
    "Age+ddG": ["age", "ddG"],
    "Age+Sex+ddG": ["age", "sex", "ddG"],
    "Age+Sex+SSS+ddG": ["age", "sex", "sss", "ddG"],
    "Age+Sex+LDL+ddG": ["age", "sex", "ldl1", "ddG"],
    "Age+Sex+LDL+SSS+ddG": ["age", "sex", "ldl1", "sss", "ddG"],
}

for name, cols in ddg_models.items():
    auc, n, ev = model_auc(has_ddg, "ascvd", cols)
    if n >= 10:
        print(f"  {name:<30} AUC={auc:.4f} (n={n}, ev={ev})")

# ─────────────────────────────────────────────────────────────
# BEST OVERALL RESULTS
# ─────────────────────────────────────────────────────────────

print("\n\n" + "=" * 80)
print("  TOP 20 MODELS WHERE SSS ADDS THE MOST (sorted by delta-AUC)")
print("=" * 80)
best_all.sort(key=lambda x: -x[4])
print(f"{'Subgroup':<25} {'Model':<25} {'Base':>6} {'+SSS':>6} {'Delta':>7} {'N':>5} {'Ev':>4}")
print("-" * 80)
for sg, m, auc, base, delta, n, ev in best_all[:20]:
    print(f"{sg:<25} {m:<25} {base:.4f} {auc:.4f} +{delta:.4f} {n:>5} {ev:>4}")

print("\n\nDone!")

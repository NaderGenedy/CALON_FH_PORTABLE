#!/usr/bin/env python3
"""
30_VUS_reclassification.py
VUS Reclassification Analysis using AlphaFold Structural Severity Scores (SSS)
Wales FH Cohort
"""

import csv
import math
import os

# ── Helper functions ──────────────────────────────────────────────────────────

def sf(x):
    try: return float(str(x).strip())
    except: return None

def odds_ratio_ci(a, b, c, d):
    if b * c == 0: return 0, 0, 0
    OR = (a * d) / (b * c)
    se = (1/max(a,1) + 1/max(b,1) + 1/max(c,1) + 1/max(d,1))**0.5
    ci_low = math.exp(math.log(max(OR, 0.001)) - 1.96 * se)
    ci_high = math.exp(math.log(max(OR, 0.001)) + 1.96 * se)
    return OR, ci_low, ci_high

def chi2_p(a, b, c, d):
    n = a + b + c + d
    if n == 0: return 0, 1.0
    ea = (a+b)*(a+c)/n; eb = (a+b)*(b+d)/n
    ec = (c+d)*(a+c)/n; ed = (c+d)*(b+d)/n
    chi2 = sum((o-e)**2/max(e, 0.001) for o, e in [(a,ea),(b,eb),(c,ec),(d,ed)])
    p = math.exp(-0.5 * chi2) if chi2 < 30 else 0.0001
    return chi2, p

def pct(num, denom):
    if denom == 0: return 0.0
    return 100.0 * num / denom

def mean_vals(vals):
    clean = [v for v in vals if v is not None]
    if not clean: return 0.0, 0
    return sum(clean)/len(clean), len(clean)

# ── File paths ────────────────────────────────────────────────────────────────

BASE = os.path.dirname(os.path.abspath(__file__))
CLINICAL_FILE = os.path.join(BASE, "WALES_FH_CLEANED (1) - Copy.csv")
SSS_FILE = os.path.join(BASE, "alphafold", "analysis", "structural_severity_scores.csv")
OUT_FILE = os.path.join(BASE, "alphafold", "analysis", "vus_reclassification_results.csv")

# ── Load SSS catalogue ───────────────────────────────────────────────────────

print("=" * 80)
print("VUS RECLASSIFICATION ANALYSIS")
print("Using AlphaFold Structural Severity Scores")
print("=" * 80)

sss_map = {}  # variant_id -> dict with sss, domain, gene, etc.
with open(SSS_FILE, "r", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    for row in reader:
        vid = row["variant_id"].strip()
        sss_val = sf(row.get("sss", ""))
        sss_map[vid] = {
            "sss": sss_val,
            "domain": row.get("domain", ""),
            "gene": row.get("gene", ""),
            "variant_type": row.get("variant_type", ""),
            "domain_severity": sf(row.get("domain_severity", "")),
        }

print(f"\nLoaded {len(sss_map)} variants from SSS catalogue")

# ── Load clinical data ────────────────────────────────────────────────────────

patients = []
with open(CLINICAL_FILE, "r", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    for row in reader:
        patients.append(row)

print(f"Loaded {len(patients)} patients from clinical data")

# ── Classify patients into groups ─────────────────────────────────────────────

pathogenic = []
vus = []
negative = []

for p in patients:
    vus1 = str(p.get("VUS1", "")).strip()
    pos1 = str(p.get("Positive1", "")).strip()

    if vus1 == "1.0":
        vus.append(p)
    elif pos1 == "1":
        pathogenic.append(p)
    elif pos1 == "0":
        negative.append(p)
    # else: unclassified (no Positive1 value), skip

print(f"\n--- Patient Classification ---")
print(f"Pathogenic (Positive1=1, not VUS):  {len(pathogenic)}")
print(f"VUS (VUS1=1.0):                     {len(vus)}")
print(f"Negative (Positive1=0, not VUS):    {len(negative)}")
print(f"Total classified:                   {len(pathogenic)+len(vus)+len(negative)}")

# ── Map SSS to patients via Mutation1 ─────────────────────────────────────────

def get_sss(patient):
    mut = str(patient.get("Mutation1", "")).strip()
    if mut in sss_map:
        return sss_map[mut]["sss"]
    return None

def get_domain(patient):
    mut = str(patient.get("Mutation1", "")).strip()
    if mut in sss_map:
        return sss_map[mut]["domain"]
    return "Unknown"

# Count how many patients in each group have SSS
for label, group in [("Pathogenic", pathogenic), ("VUS", vus), ("Negative", negative)]:
    matched = sum(1 for p in group if get_sss(p) is not None)
    print(f"  {label}: {matched}/{len(group)} have SSS mapping ({pct(matched, len(group)):.1f}%)")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 1: SSS Distribution by Group
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("ANALYSIS 1: SSS Distribution by Classification Group")
print("=" * 80)

results_rows = []

for label, group in [("Pathogenic", pathogenic), ("VUS", vus), ("Negative", negative)]:
    sss_vals = [get_sss(p) for p in group]
    sss_clean = [v for v in sss_vals if v is not None]
    if sss_clean:
        avg = sum(sss_clean) / len(sss_clean)
        sss_sorted = sorted(sss_clean)
        median = sss_sorted[len(sss_sorted)//2]
        mn = min(sss_clean)
        mx = max(sss_clean)
        high = sum(1 for v in sss_clean if v >= 0.5)
        low = sum(1 for v in sss_clean if v < 0.5)
    else:
        avg = median = mn = mx = 0
        high = low = 0

    print(f"\n  {label} (n={len(group)}, SSS mapped={len(sss_clean)}):")
    print(f"    Mean SSS:   {avg:.4f}")
    print(f"    Median SSS: {median:.4f}")
    print(f"    Range:      {mn:.2f} - {mx:.2f}")
    print(f"    High (>=0.5): {high} ({pct(high, len(sss_clean)):.1f}%)")
    print(f"    Low  (<0.5):  {low} ({pct(low, len(sss_clean)):.1f}%)")

    results_rows.append({
        "analysis": "SSS_distribution",
        "group": label,
        "n": len(group),
        "n_sss_mapped": len(sss_clean),
        "mean_sss": f"{avg:.4f}",
        "median_sss": f"{median:.4f}",
        "high_sss_count": high,
        "low_sss_count": low,
    })

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 2: Split VUS into HIGH vs LOW SSS
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("ANALYSIS 2: VUS Stratification by SSS Threshold (0.5)")
print("=" * 80)

vus_high = [p for p in vus if get_sss(p) is not None and get_sss(p) >= 0.5]
vus_low = [p for p in vus if get_sss(p) is not None and get_sss(p) < 0.5]
vus_unmapped = [p for p in vus if get_sss(p) is None]

print(f"\n  VUS HIGH SSS (>=0.5): {len(vus_high)}")
print(f"  VUS LOW SSS  (<0.5): {len(vus_low)}")
print(f"  VUS unmapped (no SSS): {len(vus_unmapped)}")

# List unique variants in each VUS subgroup
vus_high_variants = {}
for p in vus_high:
    mut = str(p.get("Mutation1", "")).strip()
    vus_high_variants[mut] = vus_high_variants.get(mut, 0) + 1
vus_low_variants = {}
for p in vus_low:
    mut = str(p.get("Mutation1", "")).strip()
    vus_low_variants[mut] = vus_low_variants.get(mut, 0) + 1

print(f"\n  VUS HIGH SSS - unique variants ({len(vus_high_variants)}):")
for v, c in sorted(vus_high_variants.items(), key=lambda x: -x[1])[:15]:
    sss_val = sss_map.get(v, {}).get("sss", "?")
    domain = sss_map.get(v, {}).get("domain", "?")
    print(f"    {v}: n={c}, SSS={sss_val}, domain={domain}")

print(f"\n  VUS LOW SSS - unique variants ({len(vus_low_variants)}):")
for v, c in sorted(vus_low_variants.items(), key=lambda x: -x[1])[:15]:
    sss_val = sss_map.get(v, {}).get("sss", "?")
    domain = sss_map.get(v, {}).get("domain", "?")
    print(f"    {v}: n={c}, SSS={sss_val}, domain={domain}")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 3: Clinical Comparisons across 4 groups
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("ANALYSIS 3: Clinical Phenotype Comparison")
print("=" * 80)

def ascvd_rate(group):
    """Return (count_ascvd, total, rate%)"""
    total = len(group)
    count = sum(1 for p in group if str(p.get("ascvd_combine", "")).strip() == "1.0")
    return count, total, pct(count, total)

def ldl_mean(group):
    """Return (mean_ldl, n_with_data)"""
    vals = [sf(p.get("LDL.1", "")) for p in group]
    return mean_vals(vals)

def tendon_rate(group):
    total = len(group)
    count = sum(1 for p in group if str(p.get("TendonXanthomata", "")).strip() not in ("", "0", "Unknown", "NoValue"))
    return count, total, pct(count, total)

def corneal_rate(group):
    total = len(group)
    count = sum(1 for p in group if str(p.get("CornealArcus", "")).strip() not in ("", "0", "Unknown", "NoValue"))
    return count, total, pct(count, total)

def simon_definite_rate(group):
    total = len(group)
    count = sum(1 for p in group if str(p.get("SimonBroome", "")).strip() == "1.0")
    return count, total, pct(count, total)

groups_4 = [
    ("Pathogenic", pathogenic),
    ("VUS_HIGH_SSS", vus_high),
    ("VUS_LOW_SSS", vus_low),
    ("Negative", negative),
]

print(f"\n{'Group':<20} {'N':>6} {'ASCVD':>10} {'LDL.1':>10} {'Tendon':>10} {'Corneal':>10} {'SB_Def':>10}")
print("-" * 80)

for label, group in groups_4:
    asc_n, asc_t, asc_r = ascvd_rate(group)
    ldl_m, ldl_n = ldl_mean(group)
    tx_n, tx_t, tx_r = tendon_rate(group)
    ca_n, ca_t, ca_r = corneal_rate(group)
    sb_n, sb_t, sb_r = simon_definite_rate(group)

    print(f"{label:<20} {len(group):>6} {asc_r:>8.1f}%  {ldl_m:>8.2f}  {tx_r:>8.1f}%  {ca_r:>8.1f}%  {sb_r:>8.1f}%")

    results_rows.append({
        "analysis": "clinical_comparison",
        "group": label,
        "n": len(group),
        "ascvd_n": asc_n,
        "ascvd_rate": f"{asc_r:.2f}",
        "ldl_mean": f"{ldl_m:.2f}",
        "ldl_n": ldl_n,
        "tendon_xanthomata_n": tx_n,
        "tendon_xanthomata_rate": f"{tx_r:.2f}",
        "corneal_arcus_n": ca_n,
        "corneal_arcus_rate": f"{ca_r:.2f}",
        "simon_broome_definite_n": sb_n,
        "simon_broome_definite_rate": f"{sb_r:.2f}",
    })

# Detailed printout
print("\n--- Detailed Clinical Data ---")
for label, group in groups_4:
    asc_n, asc_t, asc_r = ascvd_rate(group)
    ldl_m, ldl_n = ldl_mean(group)
    tx_n, tx_t, tx_r = tendon_rate(group)
    ca_n, ca_t, ca_r = corneal_rate(group)
    sb_n, sb_t, sb_r = simon_definite_rate(group)

    print(f"\n  {label} (n={len(group)}):")
    print(f"    ASCVD:                {asc_n}/{asc_t} ({asc_r:.1f}%)")
    print(f"    Mean LDL.1:           {ldl_m:.2f} mmol/L (n={ldl_n})")
    print(f"    Tendon Xanthomata:    {tx_n}/{tx_t} ({tx_r:.1f}%)")
    print(f"    Corneal Arcus:        {ca_n}/{ca_t} ({ca_r:.1f}%)")
    print(f"    Simon Broome Def:     {sb_n}/{sb_t} ({sb_r:.1f}%)")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 4: OR with 95% CI - HIGH-SSS VUS vs LOW-SSS VUS for ASCVD
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("ANALYSIS 4: Odds Ratios for ASCVD")
print("=" * 80)

def get_2x2_ascvd(group1, group2):
    """Returns a,b,c,d for 2x2 table: group1 exposed, group2 unexposed"""
    a = sum(1 for p in group1 if str(p.get("ascvd_combine", "")).strip() == "1.0")  # exposed + outcome
    b = len(group1) - a  # exposed + no outcome
    c = sum(1 for p in group2 if str(p.get("ascvd_combine", "")).strip() == "1.0")  # unexposed + outcome
    d = len(group2) - c  # unexposed + no outcome
    return a, b, c, d

# HIGH vs LOW VUS
a, b, c, d = get_2x2_ascvd(vus_high, vus_low)
OR, ci_lo, ci_hi = odds_ratio_ci(a, b, c, d)
chi2, p_val = chi2_p(a, b, c, d)

print(f"\n  HIGH-SSS VUS vs LOW-SSS VUS (ASCVD):")
print(f"    2x2 table: a={a}, b={b}, c={c}, d={d}")
print(f"    OR = {OR:.3f}  (95% CI: {ci_lo:.3f} - {ci_hi:.3f})")
print(f"    Chi-square = {chi2:.3f}, p = {p_val:.6f}")

results_rows.append({
    "analysis": "OR_ASCVD_high_vs_low_VUS",
    "group": "VUS_HIGH vs VUS_LOW",
    "OR": f"{OR:.3f}",
    "CI_low": f"{ci_lo:.3f}",
    "CI_high": f"{ci_hi:.3f}",
    "chi2": f"{chi2:.3f}",
    "p_value": f"{p_val:.6f}",
    "a": a, "b": b, "c": c, "d": d,
})

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 5: OR with 95% CI - HIGH-SSS VUS vs Negative for ASCVD
# ══════════════════════════════════════════════════════════════════════════════

a, b, c, d = get_2x2_ascvd(vus_high, negative)
OR, ci_lo, ci_hi = odds_ratio_ci(a, b, c, d)
chi2, p_val = chi2_p(a, b, c, d)

print(f"\n  HIGH-SSS VUS vs Negative (ASCVD):")
print(f"    2x2 table: a={a}, b={b}, c={c}, d={d}")
print(f"    OR = {OR:.3f}  (95% CI: {ci_lo:.3f} - {ci_hi:.3f})")
print(f"    Chi-square = {chi2:.3f}, p = {p_val:.6f}")

results_rows.append({
    "analysis": "OR_ASCVD_high_VUS_vs_negative",
    "group": "VUS_HIGH vs Negative",
    "OR": f"{OR:.3f}",
    "CI_low": f"{ci_lo:.3f}",
    "CI_high": f"{ci_hi:.3f}",
    "chi2": f"{chi2:.3f}",
    "p_value": f"{p_val:.6f}",
    "a": a, "b": b, "c": c, "d": d,
})

# Also: Pathogenic vs Negative
a, b, c, d = get_2x2_ascvd(pathogenic, negative)
OR, ci_lo, ci_hi = odds_ratio_ci(a, b, c, d)
chi2, p_val = chi2_p(a, b, c, d)

print(f"\n  Pathogenic vs Negative (ASCVD):")
print(f"    2x2 table: a={a}, b={b}, c={c}, d={d}")
print(f"    OR = {OR:.3f}  (95% CI: {ci_lo:.3f} - {ci_hi:.3f})")
print(f"    Chi-square = {chi2:.3f}, p = {p_val:.6f}")

results_rows.append({
    "analysis": "OR_ASCVD_pathogenic_vs_negative",
    "group": "Pathogenic vs Negative",
    "OR": f"{OR:.3f}",
    "CI_low": f"{ci_lo:.3f}",
    "CI_high": f"{ci_hi:.3f}",
    "chi2": f"{chi2:.3f}",
    "p_value": f"{p_val:.6f}",
    "a": a, "b": b, "c": c, "d": d,
})

# Also: LOW-SSS VUS vs Negative
a, b, c, d = get_2x2_ascvd(vus_low, negative)
OR, ci_lo, ci_hi = odds_ratio_ci(a, b, c, d)
chi2, p_val = chi2_p(a, b, c, d)

print(f"\n  LOW-SSS VUS vs Negative (ASCVD):")
print(f"    2x2 table: a={a}, b={b}, c={c}, d={d}")
print(f"    OR = {OR:.3f}  (95% CI: {ci_lo:.3f} - {ci_hi:.3f})")
print(f"    Chi-square = {chi2:.3f}, p = {p_val:.6f}")

results_rows.append({
    "analysis": "OR_ASCVD_low_VUS_vs_negative",
    "group": "VUS_LOW vs Negative",
    "OR": f"{OR:.3f}",
    "CI_low": f"{ci_lo:.3f}",
    "CI_high": f"{ci_hi:.3f}",
    "chi2": f"{chi2:.3f}",
    "p_value": f"{p_val:.6f}",
    "a": a, "b": b, "c": c, "d": d,
})

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 6: Reclassification column parsing
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("ANALYSIS 6: Reclassification Status (Reclass column)")
print("=" * 80)

# Collect all Reclass values from VUS patients
reclass_counts = {}
for p in vus:
    rc = str(p.get("Reclass", "")).strip()
    reclass_counts[rc] = reclass_counts.get(rc, 0) + 1

print(f"\n  Reclass column values among VUS patients (n={len(vus)}):")
for val, cnt in sorted(reclass_counts.items(), key=lambda x: -x[1]):
    display_val = val if val else "(empty)"
    print(f"    '{display_val}': {cnt}")

# Classify reclassification: look for pathogenic/benign keywords
reclass_path = []  # reclassified to pathogenic
reclass_benign = []  # reclassified to benign
reclass_other = []  # other reclassification
reclass_none = []  # no reclassification

pathogenic_keywords = ["pathogenic", "likely pathogenic", "lp", "p/lp", "path"]
benign_keywords = ["benign", "likely benign", "lb", "b/lb"]

for p in vus:
    rc = str(p.get("Reclass", "")).strip().lower()
    if not rc:
        reclass_none.append(p)
    elif any(kw in rc for kw in pathogenic_keywords):
        reclass_path.append(p)
    elif any(kw in rc for kw in benign_keywords):
        reclass_benign.append(p)
    else:
        reclass_other.append(p)

print(f"\n  Reclassification summary:")
print(f"    To Pathogenic/Likely Pathogenic: {len(reclass_path)}")
print(f"    To Benign/Likely Benign:         {len(reclass_benign)}")
print(f"    Other reclassification:          {len(reclass_other)}")
print(f"    No reclassification (empty):     {len(reclass_none)}")

# Show variants reclassified
if reclass_path:
    print(f"\n  Variants reclassified to PATHOGENIC:")
    for p in reclass_path:
        mut = str(p.get("Mutation1", "")).strip()
        rc = str(p.get("Reclass", "")).strip()
        sss = get_sss(p)
        sss_str = f"{sss:.2f}" if sss is not None else "N/A"
        print(f"    {mut} -> Reclass='{rc}', SSS={sss_str}")

if reclass_benign:
    print(f"\n  Variants reclassified to BENIGN:")
    for p in reclass_benign:
        mut = str(p.get("Mutation1", "")).strip()
        rc = str(p.get("Reclass", "")).strip()
        sss = get_sss(p)
        sss_str = f"{sss:.2f}" if sss is not None else "N/A"
        print(f"    {mut} -> Reclass='{rc}', SSS={sss_str}")

if reclass_other:
    print(f"\n  Variants with OTHER reclassification:")
    shown = set()
    for p in reclass_other:
        mut = str(p.get("Mutation1", "")).strip()
        rc = str(p.get("Reclass", "")).strip()
        sss = get_sss(p)
        sss_str = f"{sss:.2f}" if sss is not None else "N/A"
        key = f"{mut}|{rc}"
        if key not in shown:
            shown.add(key)
            print(f"    {mut} -> Reclass='{rc}', SSS={sss_str}")

results_rows.append({
    "analysis": "reclassification_summary",
    "reclass_to_pathogenic": len(reclass_path),
    "reclass_to_benign": len(reclass_benign),
    "reclass_other": len(reclass_other),
    "reclass_none": len(reclass_none),
})

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 7: SSS >= 0.5 as predictor of reclassification to pathogenic
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("ANALYSIS 7: SSS >= 0.5 as Predictor of Reclassification to Pathogenic")
print("=" * 80)

# Among VUS patients with reclassification data (pathogenic or benign)
reclassified = reclass_path + reclass_benign + reclass_other

if len(reclassified) > 0:
    # True positive: SSS >= 0.5 AND reclassified to pathogenic
    tp = sum(1 for p in reclass_path if get_sss(p) is not None and get_sss(p) >= 0.5)
    # False negative: SSS < 0.5 AND reclassified to pathogenic
    fn = sum(1 for p in reclass_path if get_sss(p) is not None and get_sss(p) < 0.5)
    fn_unmapped = sum(1 for p in reclass_path if get_sss(p) is None)
    # False positive: SSS >= 0.5 AND NOT reclassified to pathogenic
    fp = sum(1 for p in (reclass_benign + reclass_other) if get_sss(p) is not None and get_sss(p) >= 0.5)
    # True negative: SSS < 0.5 AND NOT reclassified to pathogenic
    tn = sum(1 for p in (reclass_benign + reclass_other) if get_sss(p) is not None and get_sss(p) < 0.5)

    sensitivity = pct(tp, tp + fn) if (tp + fn) > 0 else 0
    specificity = pct(tn, tn + fp) if (tn + fp) > 0 else 0
    ppv = pct(tp, tp + fp) if (tp + fp) > 0 else 0
    npv = pct(tn, tn + fn) if (tn + fn) > 0 else 0

    print(f"\n  Among VUS patients with ANY reclassification data (n={len(reclassified)}):")
    print(f"    To pathogenic: {len(reclass_path)}, To benign/other: {len(reclass_benign) + len(reclass_other)}")
    print(f"\n  SSS >= 0.5 prediction of pathogenic reclassification:")
    print(f"    True Positives  (SSS>=0.5 & pathogenic):  {tp}")
    print(f"    False Negatives (SSS<0.5 & pathogenic):   {fn}  (unmapped: {fn_unmapped})")
    print(f"    False Positives (SSS>=0.5 & not path):    {fp}")
    print(f"    True Negatives  (SSS<0.5 & not path):     {tn}")
    print(f"\n    Sensitivity: {sensitivity:.1f}%")
    print(f"    Specificity: {specificity:.1f}%")
    print(f"    PPV:         {ppv:.1f}%")
    print(f"    NPV:         {npv:.1f}%")

    results_rows.append({
        "analysis": "SSS_predicts_reclass",
        "tp": tp, "fn": fn, "fp": fp, "tn": tn,
        "sensitivity": f"{sensitivity:.1f}",
        "specificity": f"{specificity:.1f}",
        "ppv": f"{ppv:.1f}",
        "npv": f"{npv:.1f}",
    })
else:
    print("\n  No VUS patients with reclassification data found.")
    print("  Checking across ALL patients with Reclass column...")

    # Broader check: look at Reclass for all patients
    all_reclass = {}
    for p in patients:
        rc = str(p.get("Reclass", "")).strip()
        if rc:
            all_reclass[rc] = all_reclass.get(rc, 0) + 1
    print(f"\n  All non-empty Reclass values across entire cohort:")
    for val, cnt in sorted(all_reclass.items(), key=lambda x: -x[1]):
        print(f"    '{val}': {cnt}")

    results_rows.append({
        "analysis": "SSS_predicts_reclass",
        "note": "No reclassification data found among VUS patients",
    })

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 8: Domain-specific ASCVD rates among VUS patients
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("ANALYSIS 8: Domain-Specific ASCVD Rates Among VUS Patients")
print("=" * 80)

domain_groups = {}
for p in vus:
    dom = get_domain(p)
    if dom not in domain_groups:
        domain_groups[dom] = []
    domain_groups[dom].append(p)

print(f"\n  {'Domain':<30} {'N':>5} {'ASCVD':>8} {'Rate':>8} {'Mean SSS':>10} {'Mean LDL':>10}")
print("  " + "-" * 75)

for dom in sorted(domain_groups.keys()):
    group = domain_groups[dom]
    asc_n, asc_t, asc_r = ascvd_rate(group)
    ldl_m, ldl_n = ldl_mean(group)
    sss_vals = [get_sss(p) for p in group]
    sss_clean = [v for v in sss_vals if v is not None]
    sss_m = sum(sss_clean)/len(sss_clean) if sss_clean else 0

    print(f"  {dom:<30} {len(group):>5} {asc_n:>5}/{asc_t:<3} {asc_r:>6.1f}%  {sss_m:>8.3f}   {ldl_m:>8.2f}")

    results_rows.append({
        "analysis": "domain_ascvd_VUS",
        "domain": dom,
        "n": len(group),
        "ascvd_n": asc_n,
        "ascvd_rate": f"{asc_r:.2f}",
        "mean_sss": f"{sss_m:.3f}",
        "mean_ldl": f"{ldl_m:.2f}",
    })

# ══════════════════════════════════════════════════════════════════════════════
# Additional: LDL comparison stats
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("ADDITIONAL: LDL Distribution Details")
print("=" * 80)

for label, group in groups_4:
    vals = [sf(p.get("LDL.1", "")) for p in group]
    clean = [v for v in vals if v is not None]
    if clean:
        avg = sum(clean)/len(clean)
        srt = sorted(clean)
        med = srt[len(srt)//2]
        sd = (sum((v - avg)**2 for v in clean) / len(clean))**0.5
        print(f"\n  {label} (n={len(clean)}):")
        print(f"    Mean:   {avg:.2f} mmol/L")
        print(f"    Median: {med:.2f} mmol/L")
        print(f"    SD:     {sd:.2f}")
        print(f"    Range:  {min(clean):.2f} - {max(clean):.2f}")
    else:
        print(f"\n  {label}: No LDL data")

# ══════════════════════════════════════════════════════════════════════════════
# Additional: Chi-square tests for clinical features
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("ADDITIONAL: Chi-square Tests - VUS HIGH vs VUS LOW")
print("=" * 80)

for feature_name, feature_fn in [
    ("ASCVD", lambda p: str(p.get("ascvd_combine", "")).strip() == "1.0"),
    ("TendonXanthomata", lambda p: str(p.get("TendonXanthomata", "")).strip() not in ("", "0", "Unknown", "NoValue")),
    ("CornealArcus", lambda p: str(p.get("CornealArcus", "")).strip() not in ("", "0", "Unknown", "NoValue")),
    ("SimonBroome_Definite", lambda p: str(p.get("SimonBroome", "")).strip() == "1.0"),
]:
    a = sum(1 for p in vus_high if feature_fn(p))
    b = len(vus_high) - a
    c = sum(1 for p in vus_low if feature_fn(p))
    d = len(vus_low) - c
    OR_val, ci_lo, ci_hi = odds_ratio_ci(a, b, c, d)
    chi2_val, p_val = chi2_p(a, b, c, d)
    print(f"\n  {feature_name}:")
    print(f"    HIGH: {a}/{len(vus_high)} ({pct(a, len(vus_high)):.1f}%)  LOW: {c}/{len(vus_low)} ({pct(c, len(vus_low)):.1f}%)")
    print(f"    OR = {OR_val:.3f} (95% CI: {ci_lo:.3f} - {ci_hi:.3f})")
    print(f"    Chi2 = {chi2_val:.3f}, p = {p_val:.6f}")

# ══════════════════════════════════════════════════════════════════════════════
# Save results to CSV
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("SAVING RESULTS")
print("=" * 80)

# Collect all unique keys across result rows
all_keys = []
for row in results_rows:
    for k in row.keys():
        if k not in all_keys:
            all_keys.append(k)

os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)

with open(OUT_FILE, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=all_keys)
    writer.writeheader()
    for row in results_rows:
        writer.writerow(row)

print(f"\n  Results saved to: {OUT_FILE}")
print(f"  Total result rows: {len(results_rows)}")

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)

#!/usr/bin/env python3
"""Variant-level ApoB/LDL discordance and Lp(a) analysis."""
import csv, math, os
from collections import defaultdict

BASE = r"C:\Users\nader\Downloads\calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")

def safe_float(x, default=None):
    try:
        v = float(x)
        return default if (math.isnan(v) or math.isinf(v)) else v
    except:
        return default

# Load SSS
sss_data = {}
with open(os.path.join(ANALYSIS, "structural_severity_scores.csv")) as f:
    for row in csv.DictReader(f):
        sss_data[row["variant_id"]] = {
            "sss": safe_float(row["sss"], 0.5),
            "gene": row["gene"],
            "domain": row.get("domain", "Unknown"),
        }

sss_patients = {}
with open(os.path.join(ANALYSIS, "sss_patient_level_wales.csv")) as f:
    for row in csv.DictReader(f):
        sss_patients[row["patient_id"]] = row["variant_id"]

# Load DRAGON_3
patients = []
with open(os.path.join(BASE, "DRAGON_3.csv"), encoding="latin-1") as f:
    rdr = csv.DictReader(f)
    for row in rdr:
        pid = row.get("patient_id", row.get(rdr.fieldnames[0], ""))
        if pid not in sss_patients:
            continue
        variant = sss_patients[pid]
        vi = sss_data.get(variant, {})

        ldl1 = safe_float(row.get("LDL_1"))
        apob = safe_float(row.get("ApoB"))
        lpa = safe_float(row.get("Lpa"))
        apoa1 = safe_float(row.get("ApoA1"))
        hdl1 = safe_float(row.get("HDL_1"))
        tg1 = safe_float(row.get("TRG_1"))
        ascvd = 1 if str(row.get("ASCVD_combined", "")).strip() in ("1", "1.0", "True") else 0
        age = safe_float(row.get("age_at_event_or_censoring", row.get("Ageattest")))
        sex = "M" if str(row.get("Gender", "")).upper().startswith("M") else "F"
        statin = 1 if str(row.get("Statin", "")).strip().upper() not in ("", "NO", "NAN", "NONE", "NA", "N") else 0
        xanth = 1 if str(row.get("TendonXanthomata", "")).strip() in ("1", "1.0", "True") else 0

        apob_ldl = apob / ldl1 if apob and ldl1 and ldl1 > 0 else None

        patients.append({
            "pid": pid, "variant": variant,
            "gene": vi.get("gene", ""), "domain": vi.get("domain", ""),
            "sss": vi.get("sss", 0.5),
            "ldl1": ldl1, "apob": apob, "lpa": lpa, "apoa1": apoa1,
            "hdl1": hdl1, "tg1": tg1,
            "apob_ldl": apob_ldl,
            "ascvd": ascvd, "age": age, "sex": sex,
            "statin": statin, "xanth": xanth,
        })

# ========================================================================
print("=" * 100)
print("  VARIANT-LEVEL ApoB/LDL-C DISCORDANCE ANALYSIS")
print("=" * 100)

variant_apob = defaultdict(list)
for p in patients:
    if p["apob_ldl"] is not None:
        variant_apob[p["variant"]].append(p)

variant_stats = []
for variant, vpts in variant_apob.items():
    if len(vpts) >= 2:
        ratios = [p["apob_ldl"] for p in vpts]
        mean_r = sum(ratios) / len(ratios)
        apob_vals = [p["apob"] for p in vpts if p["apob"]]
        ldl_vals = [p["ldl1"] for p in vpts if p["ldl1"]]
        ev = sum(1 for p in vpts if p["ascvd"] == 1)
        sss = vpts[0]["sss"]
        gene = vpts[0]["gene"]
        domain = vpts[0]["domain"]

        variant_stats.append({
            "variant": variant, "n": len(vpts), "gene": gene,
            "domain": domain, "sss": sss,
            "mean_ratio": mean_r,
            "mean_apob": sum(apob_vals) / len(apob_vals) if apob_vals else 0,
            "mean_ldl": sum(ldl_vals) / len(ldl_vals) if ldl_vals else 0,
            "ascvd_n": ev, "ascvd_pct": 100 * ev / len(vpts),
        })

# Sort: HIGHEST ApoB/LDL
variant_stats.sort(key=lambda x: -x["mean_ratio"])

print("\n  HIGHEST ApoB/LDL-C RATIO (ApoB-predominant discordance)")
print("  More ApoB-containing particles per unit LDL-C -> small dense LDL pattern\n")
hdr = f"  {'Variant':<40s} {'N':<4s} {'Gene':<6s} {'Domain':<22s} {'SSS':<6s} {'ApoB/LDL':<9s} {'ApoB':<6s} {'LDL':<6s} {'ASCVD%':<7s}"
sep = f"  {'-'*40} {'-'*4} {'-'*6} {'-'*22} {'-'*6} {'-'*9} {'-'*6} {'-'*6} {'-'*7}"
print(hdr)
print(sep)
for vs in variant_stats[:20]:
    v = vs["variant"][:39]
    print(f"  {v:<40s} {vs['n']:<4d} {vs['gene']:<6s} {vs['domain']:<22s} {vs['sss']:<6.3f} {vs['mean_ratio']:<9.4f} {vs['mean_apob']:<6.2f} {vs['mean_ldl']:<6.2f} {vs['ascvd_pct']:<7.1f}")

print("\n  LOWEST ApoB/LDL-C RATIO (cholesterol-enriched large LDL particles)\n")
variant_stats.sort(key=lambda x: x["mean_ratio"])
print(hdr)
print(sep)
for vs in variant_stats[:20]:
    v = vs["variant"][:39]
    print(f"  {v:<40s} {vs['n']:<4d} {vs['gene']:<6s} {vs['domain']:<22s} {vs['sss']:<6.3f} {vs['mean_ratio']:<9.4f} {vs['mean_apob']:<6.2f} {vs['mean_ldl']:<6.2f} {vs['ascvd_pct']:<7.1f}")

# Domain-level discordance
print("\n  ApoB/LDL-C DISCORDANCE BY DOMAIN:")
domain_disc = defaultdict(list)
for p in patients:
    if p["apob_ldl"] is not None:
        domain_disc[p["domain"]].append(p["apob_ldl"])

print(f"  {'Domain':<25s} {'N':<6s} {'Mean ratio':<12s} {'SD':<8s}")
print(f"  {'-'*25} {'-'*6} {'-'*12} {'-'*8}")
for domain in sorted(domain_disc.keys(), key=lambda d: -sum(domain_disc[d]) / len(domain_disc[d])):
    vals = domain_disc[domain]
    if len(vals) >= 3:
        m = sum(vals) / len(vals)
        sd = (sum((v - m) ** 2 for v in vals) / len(vals)) ** 0.5
        print(f"  {domain:<25s} {len(vals):<6d} {m:<12.4f} {sd:<8.4f}")

# ========================================================================
print(f"\n{'=' * 100}")
print("  VARIANT-LEVEL Lp(a) ANALYSIS")
print("=" * 100)

variant_lpa = defaultdict(list)
for p in patients:
    if p["lpa"] is not None and p["lpa"] > 0:
        variant_lpa[p["variant"]].append(p)

lpa_stats = []
for variant, vpts in variant_lpa.items():
    if len(vpts) >= 2:
        lpa_vals = [p["lpa"] for p in vpts]
        mean_lpa = sum(lpa_vals) / len(lpa_vals)
        ev = sum(1 for p in vpts if p["ascvd"] == 1)
        high_lpa = sum(1 for p in vpts if p["lpa"] >= 143)
        sss = vpts[0]["sss"]
        gene = vpts[0]["gene"]
        domain = vpts[0]["domain"]
        ldl_vals = [p["ldl1"] for p in vpts if p["ldl1"]]

        lpa_stats.append({
            "variant": variant, "n": len(vpts), "gene": gene,
            "domain": domain, "sss": sss,
            "mean_lpa": mean_lpa,
            "median_lpa": sorted(lpa_vals)[len(lpa_vals) // 2],
            "high_lpa_pct": 100 * high_lpa / len(vpts),
            "mean_ldl": sum(ldl_vals) / len(ldl_vals) if ldl_vals else 0,
            "ascvd_n": ev, "ascvd_pct": 100 * ev / len(vpts),
        })

lpa_stats.sort(key=lambda x: -x["mean_lpa"])

print("\n  HIGHEST Lp(a) BY VARIANT")
print("  Lp(a) is genetically determined (LPA gene) - INDEPENDENT of LDLR/APOB.")
print("  High Lp(a) + FH = compound risk.\n")
hdr2 = f"  {'Variant':<40s} {'N':<4s} {'Gene':<6s} {'Domain':<22s} {'SSS':<6s} {'Lp(a)mean':<10s} {'Lp(a)med':<9s} {'High%':<7s} {'LDL':<6s} {'ASCVD%':<7s}"
sep2 = f"  {'-'*40} {'-'*4} {'-'*6} {'-'*22} {'-'*6} {'-'*10} {'-'*9} {'-'*7} {'-'*6} {'-'*7}"
print(hdr2)
print(sep2)
for ls in lpa_stats[:25]:
    v = ls["variant"][:39]
    print(f"  {v:<40s} {ls['n']:<4d} {ls['gene']:<6s} {ls['domain']:<22s} {ls['sss']:<6.3f} {ls['mean_lpa']:<10.1f} {ls['median_lpa']:<9.1f} {ls['high_lpa_pct']:<7.1f} {ls['mean_ldl']:<6.2f} {ls['ascvd_pct']:<7.1f}")

# Domain-level Lp(a)
print("\n  Lp(a) BY DOMAIN:")
domain_lpa = defaultdict(list)
for p in patients:
    if p["lpa"] is not None and p["lpa"] > 0:
        domain_lpa[p["domain"]].append(p["lpa"])

print(f"  {'Domain':<25s} {'N':<6s} {'Mean':<10s} {'Median':<10s} {'>=143 (%)':<10s}")
print(f"  {'-'*25} {'-'*6} {'-'*10} {'-'*10} {'-'*10}")
for domain in sorted(domain_lpa.keys(), key=lambda d: -sum(domain_lpa[d]) / len(domain_lpa[d])):
    vals = domain_lpa[domain]
    if len(vals) >= 3:
        m = sum(vals) / len(vals)
        s = sorted(vals)
        med = s[len(s) // 2]
        high = sum(1 for v in vals if v >= 143)
        print(f"  {domain:<25s} {len(vals):<6d} {m:<10.1f} {med:<10.1f} {100 * high / len(vals):<10.1f}")

# ========================================================================
print(f"\n{'=' * 100}")
print("  COMPOUND RISK: High Lp(a) + High ApoB/LDL Discordance + ASCVD")
print("=" * 100)

high_risk = [p for p in patients if p["lpa"] is not None and p["apob_ldl"] is not None]

if high_risk:
    groups = {
        "High Lp(a) + High ApoB/LDL": [p for p in high_risk if p["lpa"] >= 143 and p["apob_ldl"] >= 0.31],
        "High Lp(a) only": [p for p in high_risk if p["lpa"] >= 143 and p["apob_ldl"] < 0.31],
        "High ApoB/LDL only": [p for p in high_risk if p["lpa"] < 143 and p["apob_ldl"] >= 0.31],
        "Neither (low risk)": [p for p in high_risk if p["lpa"] < 143 and p["apob_ldl"] < 0.31],
    }

    print("\n  ASCVD RATES BY COMPOUND RISK GROUP:")
    for label, gpts in groups.items():
        if gpts:
            ev = sum(1 for p in gpts if p["ascvd"] == 1)
            mean_sss = sum(p["sss"] for p in gpts) / len(gpts)
            mean_lpa = sum(p["lpa"] for p in gpts) / len(gpts)
            mean_r = sum(p["apob_ldl"] for p in gpts) / len(gpts)
            print(f"    {label:<35s}: n={len(gpts):3d}, ASCVD={ev:2d} ({100 * ev / len(gpts):5.1f}%), SSS={mean_sss:.3f}, Lp(a)={mean_lpa:.0f}, ApoB/LDL={mean_r:.3f}")

    # Top individual patients
    lpa_sorted = sorted(p["lpa"] for p in high_risk)
    ratio_sorted = sorted(p["apob_ldl"] for p in high_risk)
    for p in high_risk:
        lpa_pctile = sum(1 for v in lpa_sorted if v <= p["lpa"]) / len(lpa_sorted)
        ratio_pctile = sum(1 for v in ratio_sorted if v <= p["apob_ldl"]) / len(ratio_sorted)
        p["risk_score"] = lpa_pctile + ratio_pctile

    high_risk.sort(key=lambda x: -x["risk_score"])

    print("\n  TOP 20 HIGHEST COMPOUND RISK PATIENTS:")
    print(f"  {'PID':<8s} {'Variant':<35s} {'SSS':<6s} {'Lp(a)':<8s} {'ApoB/LDL':<9s} {'ApoB':<6s} {'LDL':<6s} {'ASCVD':<6s} {'Domain':<20s}")
    print(f"  {'-'*8} {'-'*35} {'-'*6} {'-'*8} {'-'*9} {'-'*6} {'-'*6} {'-'*6} {'-'*20}")
    for p in high_risk[:20]:
        v = p["variant"][:34]
        asc = "YES" if p["ascvd"] else "no"
        ldl_v = f"{p['ldl1']:.2f}" if p["ldl1"] else "N/A"
        print(f"  {p['pid']:<8s} {v:<35s} {p['sss']:<6.3f} {p['lpa']:<8.1f} {p['apob_ldl']:<9.4f} {p['apob']:<6.2f} {ldl_v:<6s} {asc:<6s} {p['domain']:<20s}")

# Variant-level compound risk
print("\n  VARIANT-LEVEL COMPOUND RISK (sorted by ASCVD rate among those with both data):")
compound_v = defaultdict(list)
for p in high_risk:
    compound_v[p["variant"]].append(p)

cv_stats = []
for variant, vpts in compound_v.items():
    if len(vpts) >= 3:
        ev = sum(1 for p in vpts if p["ascvd"] == 1)
        mean_lpa = sum(p["lpa"] for p in vpts) / len(vpts)
        mean_ratio = sum(p["apob_ldl"] for p in vpts) / len(vpts)
        both = sum(1 for p in vpts if p["lpa"] >= 143 and p["apob_ldl"] >= 0.31)
        cv_stats.append({
            "variant": variant, "n": len(vpts),
            "gene": vpts[0]["gene"], "domain": vpts[0]["domain"],
            "sss": vpts[0]["sss"],
            "ascvd_pct": 100 * ev / len(vpts), "ascvd_n": ev,
            "mean_lpa": mean_lpa, "mean_ratio": mean_ratio,
            "both_high": both,
        })

cv_stats.sort(key=lambda x: -x["ascvd_pct"])
print(f"  {'Variant':<35s} {'N':<4s} {'ASCVD%':<8s} {'Lp(a)':<8s} {'ApoB/LDL':<9s} {'SSS':<6s} {'Both#':<6s} {'Domain':<20s}")
print(f"  {'-'*35} {'-'*4} {'-'*8} {'-'*8} {'-'*9} {'-'*6} {'-'*6} {'-'*20}")
for cs in cv_stats[:20]:
    v = cs["variant"][:34]
    print(f"  {v:<35s} {cs['n']:<4d} {cs['ascvd_pct']:<8.1f} {cs['mean_lpa']:<8.1f} {cs['mean_ratio']:<9.4f} {cs['sss']:<6.3f} {cs['both_high']:<6d} {cs['domain']:<20s}")

print("\n  Analysis complete.")

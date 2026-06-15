#!/usr/bin/env python3
"""
26_deep_signal_discovery.py
Nobel Prize Signal Discovery — Finding every clinically meaningful pattern
in the AlphaFold3/FoldX structural biology data
"""
import csv, math, sys
from collections import defaultdict

def sf(x):
    try: return float(str(x).strip())
    except: return None

def pearson(x, y):
    pairs = [(xi,yi) for xi,yi in zip(x,y) if xi is not None and yi is not None]
    n = len(pairs)
    if n < 5: return 0, 1.0, n
    x2, y2 = zip(*pairs)
    mx, my = sum(x2)/n, sum(y2)/n
    sxx = sum((xi-mx)**2 for xi in x2)
    syy = sum((yi-my)**2 for yi in y2)
    sxy = sum((xi-mx)*(yi-my) for xi,yi in zip(x2,y2))
    if sxx == 0 or syy == 0: return 0, 1.0, n
    r = sxy / (sxx*syy)**0.5
    return r, 0, n

def concordance(outcomes, scores):
    pos = [s for s, o in zip(scores, outcomes) if o == 1]
    neg = [s for s, o in zip(scores, outcomes) if o == 0]
    if not pos or not neg: return 0.5
    conc = sum(1 for sp in pos for sn in neg if sp > sn) + 0.5*sum(1 for sp in pos for sn in neg if sp == sn)
    return conc / (len(pos) * len(neg))

BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"
AF = f"{BASE}/alphafold/analysis"

# Load data
comp = list(csv.DictReader(open(f"{AF}/comprehensive_sss_analysis.csv", encoding="utf-8-sig")))
nobel = list(csv.DictReader(open(f"{AF}/nobel_multimodal_integration.csv", encoding="utf-8-sig")))
sss_scores = list(csv.DictReader(open(f"{AF}/structural_severity_scores.csv", encoding="utf-8-sig")))

# Try interface data
try:
    interface = list(csv.DictReader(open(f"{AF}/interface_distance_analysis.csv", encoding="utf-8-sig")))
    iface_lookup = {r["variant_id"]: sf(r["interface_dist"]) for r in interface}
except:
    iface_lookup = {}

# Build lookups
ddg_lookup = {}
plddt_lookup = {}
for r in sss_scores:
    v = r["variant_id"]
    d = sf(r.get("ddG"))
    p = sf(r.get("plddt"))
    if d is not None: ddg_lookup[v] = d
    if p is not None: plddt_lookup[v] = p

PCSK9_INTERFACE = set([25,35,50,53,58,59,60,61,78,79,81,96,316,318,319,320,321,322,325,326,328,329,330,331,332,339,341,342,344,351,352,372,386,387,390])

print("=" * 80)
print("  DEEP SIGNAL DISCOVERY - Nobel Prize Analysis")
print("=" * 80)

# SIGNAL 1: STRUCTURAL PENETRANCE
print("\n" + "=" * 80)
print("  SIGNAL 1: STRUCTURAL PENETRANCE")
print("  Does protein instability predict WHEN patients get heart attacks?")
print("=" * 80)

var_age_event = defaultdict(list)
var_ddg = {}
for r in comp:
    m = r.get("mutation", "").strip()
    ae = sf(r.get("age_event"))
    asc = sf(r.get("ascvd"))
    if asc == 1 and ae is not None:
        var_age_event[m].append(ae)
    d = ddg_lookup.get(m)
    if d is not None:
        var_ddg[m] = d

x_ddg, y_age = [], []
for m in var_age_event:
    if m in var_ddg and len(var_age_event[m]) >= 1:
        x_ddg.append(var_ddg[m])
        y_age.append(sum(var_age_event[m]) / len(var_age_event[m]))

r_val, _, n = pearson(x_ddg, y_age)
print(f"  ddG vs Mean Age at Event: r = {r_val:.3f} (n = {n} variants)")
if r_val > 0:
    print("  POSITIVE: Higher ddG -> LATER events (treatment paradox masking)")
else:
    print("  NEGATIVE: Higher ddG -> EARLIER events (structural damage accelerates)")

# Age at event by SSS tertile
events_with_age = [(sf(r["sss"]), sf(r["age_event"])) for r in comp
                   if sf(r.get("ascvd")) == 1 and sf(r.get("age_event")) is not None and sf(r.get("sss")) is not None]
events_with_age.sort(key=lambda x: x[0])
n3 = max(len(events_with_age) // 3, 1)
for tname, data in [("Low SSS", events_with_age[:n3]), ("Mid SSS", events_with_age[n3:2*n3]), ("High SSS", events_with_age[2*n3:])]:
    if data:
        ages = [d[1] for d in data]
        sss_v = [d[0] for d in data]
        print(f"  {tname} (SSS {min(sss_v):.2f}-{max(sss_v):.2f}): age at event = {sum(ages)/len(ages):.1f} +/- {(sum((a-sum(ages)/len(ages))**2 for a in ages)/max(len(ages)-1,1))**0.5:.1f} (n={len(data)})")

# SIGNAL 2: TREATMENT INTENSITY x SSS
print("\n" + "=" * 80)
print("  SIGNAL 2: TREATMENT PARADOX PROOF")
print("  Are severe variants treated more aggressively?")
print("=" * 80)

sss_treated = [sf(r["sss"]) for r in comp if sf(r.get("on_statin")) == 1 and sf(r["sss"])]
sss_untreated = [sf(r["sss"]) for r in comp if sf(r.get("on_statin")) != 1 and sf(r["sss"])]

if sss_treated and sss_untreated:
    print(f"  On statin: mean SSS = {sum(sss_treated)/len(sss_treated):.3f} (n={len(sss_treated)})")
    print(f"  No statin: mean SSS = {sum(sss_untreated)/len(sss_untreated):.3f} (n={len(sss_untreated)})")
    diff = sum(sss_treated)/len(sss_treated) - sum(sss_untreated)/len(sss_untreated)
    print(f"  Difference: {diff:+.3f}")

# Statin by event status
for label, subset in [("ASCVD event", [r for r in comp if sf(r.get("ascvd"))==1]),
                       ("No event", [r for r in comp if sf(r.get("ascvd"))!=1])]:
    statin_rate = sum(1 for r in subset if sf(r.get("on_statin"))==1) / max(len(subset), 1) * 100
    sss_v = [sf(r["sss"]) for r in subset if sf(r["sss"])]
    print(f"  {label}: n={len(subset)}, statin={statin_rate:.1f}%, mean SSS={sum(sss_v)/len(sss_v):.3f}")

# SIGNAL 3: LDL TRAJECTORY
print("\n" + "=" * 80)
print("  SIGNAL 3: LDL TRAJECTORY - Does SSS predict LDL change over time?")
print("=" * 80)

x_sss, y_ldl_change = [], []
for r in comp:
    s = sf(r.get("sss"))
    l1 = sf(r.get("ldl1"))
    ll = sf(r.get("last_ldl"))
    if s and l1 and ll and l1 > 0:
        pct_change = (ll - l1) / l1 * 100
        x_sss.append(s)
        y_ldl_change.append(pct_change)

r_val, _, n = pearson(x_sss, y_ldl_change)
print(f"  SSS vs LDL% change: r = {r_val:.3f} (n = {n})")

if y_ldl_change:
    sss_sorted = sorted(zip(x_sss, y_ldl_change))
    n3 = max(len(sss_sorted) // 3, 1)
    for tname, data in [("Low SSS", sss_sorted[:n3]), ("Mid SSS", sss_sorted[n3:2*n3]), ("High SSS", sss_sorted[2*n3:])]:
        changes = [d[1] for d in data]
        print(f"  {tname}: mean LDL change = {sum(changes)/len(changes):+.1f}% (n={len(changes)})")

# SIGNAL 4: XANTHOMATA AS STRUCTURAL BIOMARKER
print("\n" + "=" * 80)
print("  SIGNAL 4: XANTHOMATA - Physical manifestation of structural damage")
print("=" * 80)

dom_xan = defaultdict(lambda: [0, 0])
for r in comp:
    d = r.get("domain", "")
    if not d: continue
    dom_xan[d][1] += 1
    if sf(r.get("xanthomata")) == 1:
        dom_xan[d][0] += 1

print(f"  {'Domain':<30} {'Xanth%':>8} {'N':>5}")
for d in sorted(dom_xan.keys()):
    x, n = dom_xan[d]
    if n >= 3:
        print(f"  {d:<30} {x/n*100:>7.1f}% {n:>5}")

sss_xan = [sf(r["sss"]) for r in comp if sf(r.get("xanthomata")) == 1 and sf(r["sss"])]
sss_noxan = [sf(r["sss"]) for r in comp if sf(r.get("xanthomata")) != 1 and sf(r["sss"])]
if sss_xan and sss_noxan:
    print(f"\n  With xanthomata: SSS = {sum(sss_xan)/len(sss_xan):.3f} (n={len(sss_xan)})")
    print(f"  Without: SSS = {sum(sss_noxan)/len(sss_noxan):.3f} (n={len(sss_noxan)})")

# SIGNAL 5: VARIANT RARITY x SEVERITY
print("\n" + "=" * 80)
print("  SIGNAL 5: VARIANT RARITY - Rare = more pathogenic?")
print("=" * 80)

x_n = [sf(r["n"]) for r in nobel if sf(r.get("n")) and sf(r.get("ascvd_pct")) is not None]
y_asc = [sf(r["ascvd_pct"]) for r in nobel if sf(r.get("n")) and sf(r.get("ascvd_pct")) is not None]
r_val, _, n = pearson(x_n, y_asc)
print(f"  Variant frequency vs ASCVD%: r = {r_val:.3f} (n = {n})")

x_n2 = [sf(r["n"]) for r in nobel if sf(r.get("n")) and sf(r.get("sss")) is not None]
y_sss = [sf(r["sss"]) for r in nobel if sf(r.get("n")) and sf(r.get("sss")) is not None]
r_val2, _, n2 = pearson(x_n2, y_sss)
print(f"  Variant frequency vs SSS: r = {r_val2:.3f} (n = {n2})")

# SIGNAL 6: INTERFACE PROXIMITY
print("\n" + "=" * 80)
print("  SIGNAL 6: PCSK9 INTERFACE PROXIMITY - Distance to binding site")
print("=" * 80)

if iface_lookup:
    for label, lo, hi in [("AT interface (0)", -1, 0), ("NEAR (1-5)", 1, 5), ("MODERATE (6-20)", 6, 20), ("FAR (>20)", 21, 999)]:
        sub_events = []
        sub_disc = []
        sub_resp = []
        for r in comp:
            m = r.get("mutation", "").strip()
            idist = iface_lookup.get(m)
            if idist is None: continue
            if not (lo <= idist <= hi): continue
            asc = sf(r.get("ascvd"))
            apob = sf(r.get("apob"))
            ldl = sf(r.get("ldl1"))
            resp = sf(r.get("ldl_pct_change"))
            disc = apob / ldl if apob and ldl and ldl > 0 else None
            if asc is not None: sub_events.append(asc)
            if disc is not None: sub_disc.append(disc)
            if resp is not None: sub_resp.append(resp)
        n = len(sub_events)
        if n == 0: continue
        ev = sum(sub_events)
        line = f"  {label}: n={n}, ASCVD={ev:.0f}({ev/n*100:.1f}%)"
        if sub_disc: line += f", Disc={sum(sub_disc)/len(sub_disc):.3f}"
        if sub_resp: line += f", Resp={sum(sub_resp)/len(sub_resp):.1f}%"
        print(line)

# SIGNAL 7: COMPOUND STRUCTURAL RISK SCORE
print("\n" + "=" * 80)
print("  SIGNAL 7: COMPOUND STRUCTURAL RISK SCORE")
print("  SSS + Lp(a) + Discordance + Interface + Age")
print("=" * 80)

scores = []
outcomes = []
for r in comp:
    s = sf(r.get("sss"))
    age = sf(r.get("age"))
    asc = sf(r.get("ascvd"))
    lpa = sf(r.get("lpa"))
    apob = sf(r.get("apob"))
    ldl = sf(r.get("ldl1"))
    m = r.get("mutation", "").strip()
    idist = iface_lookup.get(m)

    if s is None or age is None or asc is None: continue

    s_norm = s
    age_norm = min(age / 80, 1.0)
    lpa_norm = min((lpa or 0) / 300, 1.0)
    disc_norm = min((apob / ldl if apob and ldl and ldl > 0 else 0) / 0.5, 1.0)
    iface_norm = 1.0 - min((idist or 100) / 100, 1.0)

    compound = s_norm * 0.25 + age_norm * 0.30 + lpa_norm * 0.15 + disc_norm * 0.15 + iface_norm * 0.15
    scores.append(compound)
    outcomes.append(int(asc))

auc_compound = concordance(outcomes, scores)

# Compare with age alone and SSS alone
age_only = [sf(r.get("age")) for r in comp if sf(r.get("sss")) is not None and sf(r.get("ascvd")) is not None and sf(r.get("age")) is not None]
out_only = [int(sf(r.get("ascvd"))) for r in comp if sf(r.get("sss")) is not None and sf(r.get("ascvd")) is not None and sf(r.get("age")) is not None]
auc_age = concordance(out_only, age_only)

sss_only = [sf(r.get("sss")) for r in comp if sf(r.get("sss")) is not None and sf(r.get("ascvd")) is not None]
out_sss = [int(sf(r.get("ascvd"))) for r in comp if sf(r.get("sss")) is not None and sf(r.get("ascvd")) is not None]
auc_sss = concordance(out_sss, sss_only)

print(f"  Compound Score AUC: {auc_compound:.4f} (n={len(outcomes)}, ev={sum(outcomes)})")
print(f"  Age alone AUC: {auc_age:.4f}")
print(f"  SSS alone AUC: {auc_sss:.4f}")
print(f"  Compound vs Age: +{(auc_compound-auc_age)*100:.1f}%")
print(f"  Compound vs SSS: +{(auc_compound-auc_sss)*100:.1f}%")

# Quintile analysis
score_out = sorted(zip(scores, outcomes))
n5 = max(len(score_out) // 5, 1)
print("\n  Compound Score Quintiles:")
for i, qname in enumerate(["Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"]):
    start = i * n5
    end = start + n5 if i < 4 else len(score_out)
    sub = score_out[start:end]
    ev = sum(o for _, o in sub)
    n = len(sub)
    print(f"  {qname}: n={n}, ASCVD={ev} ({ev/n*100:.1f}%)")

# SIGNAL 8: CLUSTER CHARACTERISTICS DEEP DIVE
print("\n" + "=" * 80)
print("  SIGNAL 8: PHENOTYPE CLUSTER DEEP CHARACTERISTICS")
print("=" * 80)

clusters = defaultdict(list)
for r in nobel:
    c = r.get("cluster", "")
    clusters[c].append(r)

for c in sorted(clusters.keys()):
    rows = clusters[c]
    total_pts = sum(int(sf(r.get("n", 0)) or 0) for r in rows)
    sss_v = [sf(r.get("sss", 0)) or 0 for r in rows]
    ldl_v = [sf(r.get("mean_ldl", 0)) or 0 for r in rows]
    apob_v = [sf(r.get("mean_apob", 0)) or 0 for r in rows]
    asc_v = [sf(r.get("ascvd_pct", 0)) or 0 for r in rows]
    disc_v = [sf(r.get("mean_discordance", 0)) or 0 for r in rows]
    resp_v = [sf(r.get("mean_response", 0)) or 0 for r in rows]
    lpa_v = [sf(r.get("mean_lpa", 0)) or 0 for r in rows]
    xan_v = [sf(r.get("xanth_pct", 0)) or 0 for r in rows]

    print(f"\n  {c} ({len(rows)} variants, {total_pts} patients):")
    print(f"    SSS: {sum(sss_v)/len(sss_v):.3f}")
    print(f"    LDL: {sum(ldl_v)/len(ldl_v):.2f} mmol/L")
    print(f"    ApoB: {sum(apob_v)/len(apob_v):.2f} g/L")
    print(f"    Lp(a): {sum(lpa_v)/len(lpa_v):.1f} nmol/L")
    print(f"    Discordance: {sum(disc_v)/len(disc_v):.3f}")
    print(f"    ASCVD: {sum(asc_v)/len(asc_v):.1f}%")
    print(f"    Treatment response: {sum(resp_v)/len(resp_v):.1f}%")
    print(f"    Xanthomata: {sum(xan_v)/len(xan_v):.1f}%")

    # List variants in cluster
    for r in sorted(rows, key=lambda x: sf(x.get("ascvd_pct", 0)) or 0, reverse=True)[:5]:
        print(f"      {r['variant']}: n={r['n']}, ASCVD={sf(r.get('ascvd_pct',0)):.1f}%, SSS={sf(r.get('sss',0)):.3f}")

# SIGNAL 9: DOMAIN-SPECIFIC TREATMENT RESPONSE
print("\n" + "=" * 80)
print("  SIGNAL 9: DOMAIN-SPECIFIC TREATMENT RESPONSE")
print("=" * 80)

dom_resp = defaultdict(list)
for r in comp:
    d = r.get("domain", "")
    resp = sf(r.get("ldl_pct_change"))
    if d and resp is not None:
        dom_resp[d].append(resp)

print(f"  {'Domain':<30} {'Mean Resp':>10} {'SD':>8} {'N':>5}")
for d in sorted(dom_resp.keys()):
    vals = dom_resp[d]
    if len(vals) >= 3:
        m = sum(vals) / len(vals)
        sd = (sum((v - m) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5
        print(f"  {d:<30} {m:>+9.1f}% {sd:>7.1f} {len(vals):>5}")

# SIGNAL 10: SIMON BROOME x SSS
print("\n" + "=" * 80)
print("  SIGNAL 10: SIMON BROOME RECLASSIFICATION")
print("=" * 80)

sb_groups = defaultdict(list)
for r in comp:
    sb = str(r.get("simon_broome", "")).strip()
    if sb in ("1", "1.0"):
        sb_groups["Definite"].append(r)
    elif sb in ("2", "2.0"):
        sb_groups["Possible"].append(r)
    else:
        sb_groups["Unknown/None"].append(r)

for sb, rows in sorted(sb_groups.items()):
    n = len(rows)
    ev = sum(1 for r in rows if sf(r.get("ascvd")) == 1)
    sss_v = [sf(r["sss"]) for r in rows if sf(r["sss"])]
    mean_sss = sum(sss_v) / len(sss_v) if sss_v else 0
    print(f"  {sb}: n={n}, ASCVD={ev} ({ev/n*100:.1f}%), SSS={mean_sss:.3f}")

# Within each SB category, does SSS further stratify?
for sb, rows in sorted(sb_groups.items()):
    if len(rows) < 10: continue
    sss_vals = sorted([sf(r["sss"]) for r in rows if sf(r["sss"])])
    if len(sss_vals) < 6: continue
    mid = sss_vals[len(sss_vals) // 2]
    for label, lo, hi in [("Low SSS", 0, mid), ("High SSS", mid, 1.01)]:
        sub = [r for r in rows if sf(r.get("sss")) is not None and lo <= sf(r["sss"]) < hi]
        ev = sum(1 for r in sub if sf(r.get("ascvd")) == 1)
        print(f"    {sb} + {label}: n={len(sub)}, ASCVD={ev} ({ev/max(len(sub),1)*100:.1f}%)")

print("\n" + "=" * 80)
print("  ALL SIGNALS DISCOVERED")
print("=" * 80)

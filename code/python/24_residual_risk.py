#!/usr/bin/env python3
"""
RESIDUAL RISK ANALYSIS
Can AlphaFold3/FoldX predict which patients have ASCVD
DESPITE reaching LDL-C treatment targets?
"""
import csv, math
from collections import defaultdict

def sf(x):
    try: return float(str(x).strip())
    except: return None

BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"
AF = f"{BASE}/alphafold/analysis"

comp = list(csv.DictReader(open(f"{AF}/comprehensive_sss_analysis.csv", encoding='utf-8-sig')))

print("=" * 80)
print("  RESIDUAL RISK ANALYSIS")
print("  Who gets ASCVD DESPITE reaching LDL targets?")
print("=" * 80)

# Patients on statin with last LDL data
treated = [r for r in comp if sf(r.get('on_statin'))==1 and sf(r.get('last_ldl')) is not None]
print(f"Treated patients with Last LDL: {len(treated)}")
ev_treated = sum(1 for r in treated if sf(r.get('ascvd'))==1)
print(f"Events in treated: {ev_treated} ({ev_treated/max(len(treated),1)*100:.1f}%)")

# LDL thresholds
for tname, threshold in [("LDL < 1.4 (very high risk)", 1.4),
                          ("LDL < 1.8 (high risk)", 1.8),
                          ("LDL < 2.6 (moderate risk)", 2.6),
                          ("LDL < 3.0", 3.0),
                          ("LDL < 4.0", 4.0),
                          ("LDL < 5.0", 5.0)]:
    at = [r for r in treated if sf(r['last_ldl']) < threshold]
    above = [r for r in treated if sf(r['last_ldl']) >= threshold]
    at_ev = sum(1 for r in at if sf(r.get('ascvd'))==1)
    ab_ev = sum(1 for r in above if sf(r.get('ascvd'))==1)
    at_sss = [sf(r['sss']) for r in at if sf(r['sss'])]
    ab_sss = [sf(r['sss']) for r in above if sf(r['sss'])]

    print(f"\n{tname}:")
    print(f"  At target: n={len(at)}, ASCVD={at_ev} ({at_ev/max(len(at),1)*100:.1f}%)", end="")
    if at_sss: print(f", mean SSS={sum(at_sss)/len(at_sss):.3f}", end="")
    print()
    print(f"  Above:     n={len(above)}, ASCVD={ab_ev} ({ab_ev/max(len(above),1)*100:.1f}%)", end="")
    if ab_sss: print(f", mean SSS={sum(ab_sss)/len(ab_sss):.3f}", end="")
    print()

# KEY: SSS predicts residual risk at each threshold
print("\n" + "=" * 80)
print("  SSS PREDICTS RESIDUAL RISK AT EACH THRESHOLD")
print("=" * 80)

for threshold in [2.6, 3.0, 4.0, 5.0]:
    at_target = [r for r in treated if sf(r['last_ldl']) < threshold]
    if len(at_target) < 10: continue

    events = [r for r in at_target if sf(r.get('ascvd'))==1]
    no_events = [r for r in at_target if sf(r.get('ascvd'))!=1]

    ev_sss = [sf(r['sss']) for r in events if sf(r['sss'])]
    noev_sss = [sf(r['sss']) for r in no_events if sf(r['sss'])]

    print(f"\nLDL < {threshold} (n={len(at_target)}, events={len(events)}):")
    if ev_sss and noev_sss:
        m1 = sum(ev_sss)/len(ev_sss)
        m2 = sum(noev_sss)/len(noev_sss)
        print(f"  Events SSS:    {m1:.3f} +/- {(sum((s-m1)**2 for s in ev_sss)/max(len(ev_sss)-1,1))**0.5:.3f} (n={len(ev_sss)})")
        print(f"  No-events SSS: {m2:.3f} +/- {(sum((s-m2)**2 for s in noev_sss)/max(len(noev_sss)-1,1))**0.5:.3f} (n={len(noev_sss)})")

        # Welch t-test approximation
        n1, n2 = len(ev_sss), len(noev_sss)
        s1 = (sum((s-m1)**2 for s in ev_sss)/max(n1-1,1))**0.5
        s2 = (sum((s-m2)**2 for s in noev_sss)/max(n2-1,1))**0.5
        se = (s1**2/n1 + s2**2/n2)**0.5
        if se > 0:
            t = (m1-m2)/se
            # Rough p approximation
            df = n1+n2-2
            p = 2 * (1 - 0.5*(1+math.erf(abs(t)/math.sqrt(2))))
            print(f"  Difference: {m1-m2:+.3f}, t={t:.2f}, P={p:.3f}")

        # SSS tertiles within at-target
        all_sss = sorted([sf(r['sss']) for r in at_target if sf(r['sss'])])
        if len(all_sss) >= 9:
            t1 = all_sss[len(all_sss)//3]
            t2 = all_sss[2*len(all_sss)//3]
            for tlabel, lo, hi in [("Low SSS", 0, t1), ("Mid SSS", t1, t2), ("High SSS", t2, 1.01)]:
                sub = [r for r in at_target if sf(r.get('sss')) is not None and lo <= sf(r['sss']) < hi]
                ev = sum(1 for r in sub if sf(r.get('ascvd'))==1)
                if sub:
                    print(f"  {tlabel} ({lo:.2f}-{hi:.2f}): n={len(sub)}, ASCVD={ev} ({ev/len(sub)*100:.1f}%)")

# Domain residual risk
print("\n" + "=" * 80)
print("  DOMAIN-SPECIFIC RESIDUAL RISK")
print("=" * 80)

dom_r = defaultdict(lambda: {'treated':0, 'events':0, 'at3':0, 'at3_ev':0})
for r in treated:
    d = r.get('domain','')
    if not d: continue
    dom_r[d]['treated'] += 1
    if sf(r.get('ascvd'))==1: dom_r[d]['events'] += 1
    if sf(r['last_ldl']) < 3.0:
        dom_r[d]['at3'] += 1
        if sf(r.get('ascvd'))==1: dom_r[d]['at3_ev'] += 1

print(f"{'Domain':<30} {'Treated':>8} {'Events':>7} {'Rate':>6} {'AtLDL<3':>8} {'ResidRisk':>10}")
print("-" * 75)
for d in sorted(dom_r.keys()):
    v = dom_r[d]
    if v['treated'] < 3: continue
    rate = v['events']/v['treated']*100
    resid = v['at3_ev']/v['at3']*100 if v['at3'] > 0 else 0
    print(f"{d:<30} {v['treated']:>8} {v['events']:>7} {rate:>5.1f}% {v['at3']:>8} {resid:>9.1f}%")

# FoldX ddG and residual risk
print("\n" + "=" * 80)
print("  ddG AND RESIDUAL RISK")
print("=" * 80)

foldx_map = {}
for r in csv.DictReader(open(f"{AF}/foldx_ddg_results.csv", encoding='utf-8-sig')):
    vid = r.get('variant_id','').strip()
    ddg = sf(r.get('ddG_kcal_mol'))
    if vid and ddg is not None: foldx_map[vid] = ddg
for r in csv.DictReader(open(f"{AF}/structural_severity_scores.csv", encoding='utf-8-sig')):
    vid = r.get('variant_id','').strip()
    ddg = sf(r.get('ddG'))
    if vid and ddg is not None and vid not in foldx_map: foldx_map[vid] = ddg

ddg_ev, ddg_noev = [], []
for r in treated:
    m = r.get('mutation','').strip()
    ddg = foldx_map.get(m)
    if ddg is None: continue
    if sf(r.get('ascvd'))==1: ddg_ev.append(ddg)
    else: ddg_noev.append(ddg)

if ddg_ev and ddg_noev:
    m1 = sum(ddg_ev)/len(ddg_ev)
    m2 = sum(ddg_noev)/len(ddg_noev)
    print(f"Events: mean ddG = {m1:.2f} (n={len(ddg_ev)})")
    print(f"No-events: mean ddG = {m2:.2f} (n={len(ddg_noev)})")
    print(f"Difference: {m1-m2:+.2f}")

# Variant-level residual risk
print("\n" + "=" * 80)
print("  VARIANTS WITH ASCVD DESPITE TREATMENT")
print("=" * 80)

var_r = defaultdict(lambda: {'treated':0, 'events':0, 'sss':[], 'last_ldl':[], 'domain':''})
for r in treated:
    m = r.get('mutation','')
    if not m: continue
    var_r[m]['treated'] += 1
    if sf(r.get('sss')): var_r[m]['sss'].append(sf(r['sss']))
    if sf(r.get('last_ldl')): var_r[m]['last_ldl'].append(sf(r['last_ldl']))
    var_r[m]['domain'] = r.get('domain','')
    if sf(r.get('ascvd'))==1: var_r[m]['events'] += 1

resistant = [(m, v) for m, v in var_r.items() if v['treated'] >= 2 and v['events'] >= 1]
resistant.sort(key=lambda x: -x[1]['events']/x[1]['treated'])

print(f"{'Variant':<35} {'Tx':>3} {'Ev':>3} {'Rate':>6} {'SSS':>6} {'LastLDL':>8} {'Domain':<25}")
print("-" * 95)
for m, v in resistant[:25]:
    rate = v['events']/v['treated']*100
    sss = sum(v['sss'])/len(v['sss']) if v['sss'] else 0
    ldl = sum(v['last_ldl'])/len(v['last_ldl']) if v['last_ldl'] else 0
    print(f"{m:<35} {v['treated']:>3} {v['events']:>3} {rate:>5.1f}% {sss:>5.3f} {ldl:>7.1f} {v['domain']:<25}")

# UKB residual risk
print("\n" + "=" * 80)
print("  UKB RESIDUAL RISK (EXTERNAL VALIDATION)")
print("=" * 80)

ukb = list(csv.DictReader(open(f"{BASE}/calon_ukb_analysis_ready.csv", encoding='utf-8-sig')))
ukb_statin = [r for r in ukb if sf(r.get('on_statin'))==1]
print(f"UKB on statin: {len(ukb_statin)}")

# Gene-level event rates in treated UKB
for gene in ['LDLR', 'APOB']:
    sub = [r for r in ukb_statin if r.get('gene')==gene]
    ev = sum(1 for r in sub if sf(r.get('ascvd_combined'))==1)
    ldl_vals = [sf(r.get('ldl')) for r in sub if sf(r.get('ldl'))]
    print(f"\n{gene} on statin: n={len(sub)}, ASCVD={ev} ({ev/max(len(sub),1)*100:.1f}%)")
    if ldl_vals: print(f"  Mean LDL: {sum(ldl_vals)/len(ldl_vals):.2f}")

    # At target
    for threshold in [1.8, 2.6, 3.0]:
        at = [r for r in sub if sf(r.get('ldl')) is not None and sf(r['ldl']) < threshold]
        at_ev = sum(1 for r in at if sf(r.get('ascvd_combined'))==1)
        if at:
            print(f"  LDL < {threshold}: n={len(at)}, ASCVD={at_ev} ({at_ev/len(at)*100:.1f}%)")

# LDLR vs APOB residual risk comparison
print("\n" + "=" * 80)
print("  GENE-LEVEL RESIDUAL RISK: LDLR vs APOB")
print("=" * 80)

for threshold in [2.6, 3.0]:
    print(f"\nPatients at LDL < {threshold}:")
    for gene in ['LDLR', 'APOB']:
        sub = [r for r in ukb_statin if r.get('gene')==gene and sf(r.get('ldl')) is not None and sf(r['ldl']) < threshold]
        ev = sum(1 for r in sub if sf(r.get('ascvd_combined'))==1)
        if sub:
            apob_vals = [sf(r.get('apob')) for r in sub if sf(r.get('apob'))]
            lpa_vals = [sf(r.get('lpa')) for r in sub if sf(r.get('lpa'))]
            print(f"  {gene}: n={len(sub)}, ASCVD={ev} ({ev/len(sub)*100:.1f}%)", end="")
            if apob_vals: print(f", ApoB={sum(apob_vals)/len(apob_vals):.2f}", end="")
            if lpa_vals: print(f", Lp(a)={sum(lpa_vals)/len(lpa_vals):.1f}", end="")
            print()

print("\n" + "=" * 80)
print("  ANALYSIS COMPLETE")
print("=" * 80)

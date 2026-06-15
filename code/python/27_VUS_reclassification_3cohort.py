#!/usr/bin/env python3
"""
27_VUS_reclassification_3cohort.py
VUS Reclassification using AlphaFold3/FoldX + 3-Cohort External Validation
Development: All Wales (7,253) | Validation 1: DRAGON3 (1,362) | Validation 2: UKB (1,623)

Author: Dr Nader Genedy
"""

import csv, math, os
from collections import Counter, defaultdict

BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"
AF = f"{BASE}/alphafold/analysis"

def sf(x):
    try: return float(str(x).strip())
    except: return None

def mean_sd(vals):
    vals = [v for v in vals if v is not None]
    if not vals: return 0, 0, 0
    n = len(vals)
    m = sum(vals)/n
    sd = (sum((v-m)**2 for v in vals)/(n-1))**0.5 if n > 1 else 0
    return m, sd, n

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
    if abs(r) >= 1.0: return r, 0, n
    t = r * ((n-2)/(1-r*r+1e-15))**0.5
    df = n - 2
    p = (df/(df+t*t))**(df/2)
    return r, p, n

def concordance(outcomes, scores):
    pos = [s for s, o in zip(scores, outcomes) if o == 1]
    neg = [s for s, o in zip(scores, outcomes) if o == 0]
    if not pos or not neg: return 0.5
    conc = sum(1 for sp in pos for sn in neg if sp > sn) + 0.5*sum(1 for sp in pos for sn in neg if sp == sn)
    return conc / (len(pos) * len(neg))

def odds_ratio_ci(a, b, c, d):
    if b == 0 or c == 0 or a == 0 or d == 0:
        return 0, "NA"
    OR = (a*d)/(b*c)
    se = (1/max(a,1) + 1/max(b,1) + 1/max(c,1) + 1/max(d,1))**0.5
    lo = math.exp(math.log(max(OR,0.001)) - 1.96*se)
    hi = math.exp(math.log(max(OR,0.001)) + 1.96*se)
    return OR, f"{lo:.2f}-{hi:.2f}"

def chi2_test(a, b, c, d):
    n = a+b+c+d
    if n == 0: return 0, 1.0
    ea = (a+b)*(a+c)/n; eb = (a+b)*(b+d)/n
    ec = (c+d)*(a+c)/n; ed = (c+d)*(b+d)/n
    chi2 = sum((o-e)**2/max(e,0.001) for o,e in [(a,ea),(b,eb),(c,ec),(d,ed)])
    p = math.exp(-chi2/2) if chi2 < 20 else 0.0001
    return chi2, p

def group_report(name, data, ascvd_col='ascvd_combine', ldl_col='LDL.1'):
    n = len(data)
    if n == 0:
        print(f"  {name}: n=0")
        return
    ev = sum(1 for r in data if sf(r.get(ascvd_col)) == 1)
    ldl = [sf(r.get(ldl_col)) for r in data if sf(r.get(ldl_col))]
    tc = [sf(r.get('TC.1', r.get('tc',''))) for r in data if sf(r.get('TC.1', r.get('tc','')))]
    xan_count = sum(1 for r in data if str(r.get('TendonXanthomata', r.get('xanthomata',''))).strip() not in ['','0','No','no','0.0'])
    ca_count = sum(1 for r in data if str(r.get('CornealArcus','')).strip() not in ['','0','No','no','0.0'])

    m_ldl, sd_ldl, n_ldl = mean_sd(ldl)
    m_tc, sd_tc, n_tc = mean_sd(tc)

    print(f"  {name}: n={n}")
    print(f"    ASCVD: {ev}/{n} ({ev/n*100:.1f}%)")
    if n_ldl > 0: print(f"    LDL: {m_ldl:.2f} +/- {sd_ldl:.2f} (n={n_ldl})")
    if n_tc > 0: print(f"    TC: {m_tc:.2f} +/- {sd_tc:.2f} (n={n_tc})")
    print(f"    Xanthomata: {xan_count} ({xan_count/n*100:.1f}%)")
    print(f"    Corneal arcus: {ca_count} ({ca_count/n*100:.1f}%)")

# ═══════════════════════════════════════════════════════════════
# LOAD DATA
# ═══════════════════════════════════════════════════════════════

print("Loading data...")

# Load SSS lookup (with expanded FoldX)
sss_lookup = {}
for fname in ['structural_severity_scores_full.csv', 'structural_severity_scores.csv']:
    fpath = f"{AF}/{fname}"
    if os.path.exists(fpath):
        for r in csv.DictReader(open(fpath, encoding='utf-8-sig')):
            vid = r.get('variant_id','').strip()
            s = sf(r.get('sss'))
            ddg = sf(r.get('ddG'))
            dom = r.get('domain','')
            plddt = sf(r.get('plddt'))
            if vid and s is not None:
                if vid not in sss_lookup or (ddg is not None and sss_lookup[vid].get('ddG') is None):
                    sss_lookup[vid] = {'sss': s, 'ddG': ddg, 'domain': dom, 'plddt': plddt}

print(f"  SSS lookup: {len(sss_lookup)} variants")

# Cohort 1: All Wales (DEVELOPMENT)
wales = list(csv.DictReader(open(f"{BASE}/WALES_FH_CLEANED (1) - Copy.csv", encoding='utf-8-sig')))
print(f"  All Wales: {len(wales)} patients")

# Cohort 2: DRAGON3 (South Wales, EXTERNAL VALIDATION 1)
dragon = list(csv.DictReader(open(f"{BASE}/DRAGON_3.csv", encoding='utf-8-sig')))
print(f"  DRAGON3: {len(dragon)} patients")

# Cohort 3: UKB (EXTERNAL VALIDATION 2)
ukb = list(csv.DictReader(open(f"{BASE}/calon_ukb_analysis_ready.csv", encoding='utf-8-sig')))
print(f"  UKB: {len(ukb)} patients")

# ═══════════════════════════════════════════════════════════════
# SECTION A: MAP SSS TO ALL WALES
# ═══════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("  SECTION A: MAPPING SSS TO ALL WALES COHORT (n=7,253)")
print("="*80)

# Classify patients
for r in wales:
    mut = r.get('Mutation1','').strip()
    is_vus = sf(r.get('VUS1')) == 1.0
    is_positive = r.get('Positive1','').strip() == '1'

    if is_vus:
        r['_class'] = 'VUS'
    elif is_positive:
        r['_class'] = 'Pathogenic'
    else:
        r['_class'] = 'Negative'

    # Map SSS
    info = sss_lookup.get(mut)
    if info:
        r['_sss'] = info['sss']
        r['_ddg'] = info['ddG']
        r['_domain'] = info['domain']
        r['_plddt'] = info['plddt']
    else:
        r['_sss'] = None
        r['_ddg'] = None
        r['_domain'] = None
        r['_plddt'] = None

# Report mapping
classes = Counter(r['_class'] for r in wales)
print(f"\nClassification: {dict(classes)}")

for cls in ['Pathogenic', 'VUS', 'Negative']:
    sub = [r for r in wales if r['_class'] == cls]
    mapped = sum(1 for r in sub if r['_sss'] is not None)
    print(f"  {cls}: {len(sub)} total, {mapped} with SSS ({mapped/len(sub)*100:.1f}%)")

# ═══════════════════════════════════════════════════════════════
# SECTION B: VUS RECLASSIFICATION — THE KEY ANALYSIS
# ═══════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("  SECTION B: VUS RECLASSIFICATION USING AlphaFold3/FoldX")
print("="*80)

vus = [r for r in wales if r['_class'] == 'VUS']
vus_with_sss = [r for r in vus if r['_sss'] is not None]
pathogenic = [r for r in wales if r['_class'] == 'Pathogenic']
path_with_sss = [r for r in pathogenic if r['_sss'] is not None]
negative = [r for r in wales if r['_class'] == 'Negative']

print(f"\nVUS with AlphaFold3/FoldX data: {len(vus_with_sss)} / {len(vus)}")

# Split VUS by SSS threshold
for threshold_name, threshold in [("SSS >= 0.6 (likely pathogenic)", 0.6),
                                    ("SSS >= 0.5 (probably pathogenic)", 0.5),
                                    ("SSS >= 0.4 (possibly pathogenic)", 0.4)]:
    high = [r for r in vus_with_sss if r['_sss'] >= threshold]
    low = [r for r in vus_with_sss if r['_sss'] < threshold]

    print(f"\n--- Threshold: {threshold_name} ---")
    group_report(f"VUS SSS >= {threshold} (RECLASSIFY?)", high)
    group_report(f"VUS SSS < {threshold} (LIKELY BENIGN?)", low)
    group_report("Pathogenic (reference)", path_with_sss)
    group_report("Negative (reference)", negative)

    # Statistical test: High-SSS VUS vs Low-SSS VUS
    a = sum(1 for r in high if sf(r.get('ascvd_combine')) == 1)
    b = len(high) - a
    c = sum(1 for r in low if sf(r.get('ascvd_combine')) == 1)
    d = len(low) - c

    if (a+b) > 0 and (c+d) > 0:
        OR, ci = odds_ratio_ci(a, b, c, d)
        chi2, p = chi2_test(a, b, c, d)
        print(f"\n  ASCVD: High-SSS VUS {a}/{a+b} ({a/(a+b)*100:.1f}%) vs Low-SSS VUS {c}/{c+d} ({c/(c+d)*100:.1f}%)")
        print(f"  OR = {OR:.2f} (95% CI: {ci})")
        print(f"  Chi-squared = {chi2:.2f}, P ~ {p:.4f}")

# ═══════════════════════════════════════════════════════════════
# SECTION C: VUS VARIANT-LEVEL ANALYSIS
# ═══════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("  SECTION C: VUS VARIANT-LEVEL DEEP DIVE")
print("="*80)

vus_var = defaultdict(lambda: {'n':0, 'ascvd':0, 'ldl':[], 'xanth':0, 'ca':0, 'sss':None, 'ddG':None, 'domain':'', 'plddt':None})
for r in vus:
    mut = r.get('Mutation1','').strip()
    if not mut: continue
    v = vus_var[mut]
    v['n'] += 1
    if sf(r.get('ascvd_combine')) == 1: v['ascvd'] += 1
    l = sf(r.get('LDL.1'))
    if l: v['ldl'].append(l)
    if str(r.get('TendonXanthomata','')).strip() not in ['','0','No']: v['xanth'] += 1
    if str(r.get('CornealArcus','')).strip() not in ['','0','No']: v['ca'] += 1
    if r['_sss'] is not None: v['sss'] = r['_sss']
    if r['_ddg'] is not None: v['ddG'] = r['_ddg']
    v['domain'] = r.get('_domain','') or ''
    if r.get('_plddt') is not None: v['plddt'] = r['_plddt']

# VUS with SSS data
has_sss_var = [(m,v) for m,v in vus_var.items() if v['sss'] is not None]
has_sss_var.sort(key=lambda x: -x[1]['sss'])

print(f"\nVUS variants with SSS data: {len(has_sss_var)}")
print(f"{'Variant':<35} {'N':>4} {'ASCVD%':>7} {'LDL':>6} {'Xanth%':>7} {'SSS':>6} {'ddG':>7} {'pLDDT':>6} {'Domain':<20} {'Reclass?':<15}")
print("-"*120)

reclassified_path = 0
reclassified_benign = 0
for m, v in has_sss_var:
    n = v['n']
    ascvd_pct = v['ascvd']/n*100 if n > 0 else 0
    ldl_mean = sum(v['ldl'])/len(v['ldl']) if v['ldl'] else 0
    xanth_pct = v['xanth']/n*100
    ddg_str = f"{v['ddG']:.1f}" if v['ddG'] is not None else "NA"
    plddt_str = f"{v['plddt']:.0f}" if v['plddt'] is not None else "NA"

    # Reclassification logic
    if v['sss'] >= 0.5 and (v['ddG'] is not None and v['ddG'] >= 2.0):
        reclass = "LIKELY PATH"
        reclassified_path += 1
    elif v['sss'] >= 0.5:
        reclass = "PROB PATH"
        reclassified_path += 1
    elif v['sss'] < 0.3:
        reclass = "LIKELY BENIGN"
        reclassified_benign += 1
    else:
        reclass = "UNCERTAIN"

    if n >= 2:  # Show variants with >=2 carriers
        print(f"{m:<35} {n:>4} {ascvd_pct:>6.1f}% {ldl_mean:>5.2f} {xanth_pct:>6.1f}% {v['sss']:>5.3f} {ddg_str:>7} {plddt_str:>6} {v['domain']:<20} {reclass:<15}")

print(f"\nReclassification summary:")
print(f"  Likely/Probably Pathogenic: {reclassified_path} VUS variants")
print(f"  Likely Benign: {reclassified_benign} VUS variants")
print(f"  Uncertain: {len(has_sss_var) - reclassified_path - reclassified_benign} VUS variants")

# ═══════════════════════════════════════════════════════════════
# SECTION D: ddG-CLINICAL CORRELATIONS IN VUS
# ═══════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("  SECTION D: FoldX ddG PREDICTS CLINICAL PHENOTYPE IN VUS")
print("="*80)

# Variant-level correlations
vus_ddg_ldl_x, vus_ddg_ldl_y = [], []
vus_ddg_ascvd_x, vus_ddg_ascvd_y = [], []
vus_sss_ldl_x, vus_sss_ldl_y = [], []
vus_sss_ascvd_x, vus_sss_ascvd_y = [], []

for m, v in has_sss_var:
    if v['n'] < 2: continue
    if v['ddG'] is not None and v['ldl']:
        vus_ddg_ldl_x.append(v['ddG'])
        vus_ddg_ldl_y.append(sum(v['ldl'])/len(v['ldl']))
    if v['ddG'] is not None:
        vus_ddg_ascvd_x.append(v['ddG'])
        vus_ddg_ascvd_y.append(v['ascvd']/v['n']*100)
    if v['ldl']:
        vus_sss_ldl_x.append(v['sss'])
        vus_sss_ldl_y.append(sum(v['ldl'])/len(v['ldl']))
    vus_sss_ascvd_x.append(v['sss'])
    vus_sss_ascvd_y.append(v['ascvd']/v['n']*100)

r1, p1, n1 = pearson(vus_ddg_ldl_x, vus_ddg_ldl_y)
print(f"ddG vs LDL (VUS variants, n>= 2 carriers): r = {r1:.3f}, P = {p1:.4f}, n = {n1}")

r2, p2, n2 = pearson(vus_ddg_ascvd_x, vus_ddg_ascvd_y)
print(f"ddG vs ASCVD% (VUS variants): r = {r2:.3f}, P = {p2:.4f}, n = {n2}")

r3, p3, n3 = pearson(vus_sss_ldl_x, vus_sss_ldl_y)
print(f"SSS vs LDL (VUS variants): r = {r3:.3f}, P = {p3:.4f}, n = {n3}")

r4, p4, n4 = pearson(vus_sss_ascvd_x, vus_sss_ascvd_y)
print(f"SSS vs ASCVD% (VUS variants): r = {r4:.3f}, P = {p4:.4f}, n = {n4}")

# ═══════════════════════════════════════════════════════════════
# SECTION E: SSS AS PREDICTOR — AUC IN ALL WALES
# ═══════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("  SECTION E: SSS PREDICTIVE PERFORMANCE IN ALL WALES (DEVELOPMENT)")
print("="*80)

# All patients with SSS and ASCVD data
wales_analytic = [r for r in wales if r['_sss'] is not None and sf(r.get('ascvd_combine')) is not None]
print(f"Analytic sample: {len(wales_analytic)} patients with SSS + ASCVD")
ev_total = sum(1 for r in wales_analytic if sf(r.get('ascvd_combine')) == 1)
print(f"Events: {ev_total}")

# AUC: SSS alone
outcomes = [int(sf(r.get('ascvd_combine'))) for r in wales_analytic]
sss_scores_list = [r['_sss'] for r in wales_analytic]
auc_sss = concordance(outcomes, sss_scores_list)
print(f"\nSSS alone AUC: {auc_sss:.4f}")

# AUC: SSS + clinical
# Need age - compute from DOB or use available age columns
# Check what age columns exist
age_vals = []
for r in wales_analytic:
    for col in ['BMI_AGE', 'Currentage', 'corneal_age']:
        a = sf(r.get(col))
        if a and 10 < a < 100:
            r['_age'] = a
            break
    else:
        r['_age'] = None

wales_with_age = [r for r in wales_analytic if r['_age'] is not None]
print(f"With age data: {len(wales_with_age)}")

if wales_with_age:
    outcomes_a = [int(sf(r.get('ascvd_combine'))) for r in wales_with_age]

    # AUC: Age alone
    age_scores = [r['_age'] for r in wales_with_age]
    auc_age = concordance(outcomes_a, age_scores)

    # AUC: Age + SSS (simple sum of z-scores)
    sss_v = [r['_sss'] for r in wales_with_age]
    m_age = sum(age_scores)/len(age_scores)
    sd_age = (sum((a-m_age)**2 for a in age_scores)/(len(age_scores)-1))**0.5
    m_sss = sum(sss_v)/len(sss_v)
    sd_sss = (sum((s-m_sss)**2 for s in sss_v)/(len(sss_v)-1))**0.5

    combined = [(a-m_age)/max(sd_age,0.001) + (s-m_sss)/max(sd_sss,0.001) for a, s in zip(age_scores, sss_v)]
    auc_combined = concordance(outcomes_a, combined)

    ev_a = sum(outcomes_a)
    print(f"\nAge alone AUC: {auc_age:.4f} (n={len(wales_with_age)}, events={ev_a})")
    print(f"SSS alone AUC: {concordance(outcomes_a, sss_v):.4f}")
    print(f"Age + SSS AUC: {auc_combined:.4f}")
    print(f"Delta (Age+SSS vs Age): +{(auc_combined-auc_age)*100:.1f}%")

# ═══════════════════════════════════════════════════════════════
# SECTION F: EXTERNAL VALIDATION 1 — DRAGON3
# ═══════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("  SECTION F: EXTERNAL VALIDATION — DRAGON3 (South Wales, n=1,362)")
print("="*80)

# Map SSS to DRAGON3
for r in dragon:
    mut = r.get('Mutation1','').strip()
    info = sss_lookup.get(mut)
    if info:
        r['_sss'] = info['sss']
    else:
        r['_sss'] = None

dragon_with_sss = [r for r in dragon if r['_sss'] is not None and sf(r.get('ASCVD_combined')) is not None]
print(f"DRAGON3 with SSS: {len(dragon_with_sss)}")
ev_d = sum(1 for r in dragon_with_sss if sf(r.get('ASCVD_combined')) == 1)
print(f"Events: {ev_d}")

if dragon_with_sss:
    outcomes_d = [int(sf(r.get('ASCVD_combined'))) for r in dragon_with_sss]
    sss_d = [r['_sss'] for r in dragon_with_sss]
    auc_d = concordance(outcomes_d, sss_d)
    print(f"SSS alone AUC: {auc_d:.4f}")

    # Age + SSS
    ages_d = [sf(r.get('age_at_event_or_censoring')) for r in dragon_with_sss]
    has_age_d = [(r, a) for r, a in zip(dragon_with_sss, ages_d) if a is not None]
    if has_age_d:
        out_d2 = [int(sf(r.get('ASCVD_combined'))) for r, a in has_age_d]
        age_d2 = [a for r, a in has_age_d]
        sss_d2 = [r['_sss'] for r, a in has_age_d]

        auc_age_d = concordance(out_d2, age_d2)

        ma = sum(age_d2)/len(age_d2)
        sa = (sum((a-ma)**2 for a in age_d2)/(len(age_d2)-1))**0.5
        ms = sum(sss_d2)/len(sss_d2)
        ss = (sum((s-ms)**2 for s in sss_d2)/(len(sss_d2)-1))**0.5
        comb_d = [(a-ma)/max(sa,0.001) + (s-ms)/max(ss,0.001) for a, s in zip(age_d2, sss_d2)]
        auc_comb_d = concordance(out_d2, comb_d)

        print(f"Age alone AUC: {auc_age_d:.4f} (n={len(has_age_d)}, events={sum(out_d2)})")
        print(f"Age + SSS AUC: {auc_comb_d:.4f}")
        print(f"Delta: +{(auc_comb_d-auc_age_d)*100:.1f}%")

# ═══════════════════════════════════════════════════════════════
# SECTION G: EXTERNAL VALIDATION 2 — UKB
# ═══════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("  SECTION G: EXTERNAL VALIDATION — UK BIOBANK (n=1,623)")
print("="*80)

# UKB has gene but not variant — use gene-level SSS
gene_sss = {'LDLR': 0.645, 'APOB': 0.473, 'PCSK9': 0.575}
for r in ukb:
    g = r.get('gene','')
    r['_sss'] = gene_sss.get(g, 0.5)

ukb_analytic = [r for r in ukb if sf(r.get('ascvd_combined')) is not None]
print(f"UKB analytic: {len(ukb_analytic)}")
ev_u = sum(1 for r in ukb_analytic if sf(r.get('ascvd_combined')) == 1)
print(f"Events: {ev_u}")

outcomes_u = [int(sf(r.get('ascvd_combined'))) for r in ukb_analytic]
sss_u = [r['_sss'] for r in ukb_analytic]
ages_u = [sf(r.get('age')) for r in ukb_analytic]

auc_sss_u = concordance(outcomes_u, sss_u)
print(f"Gene-SSS alone AUC: {auc_sss_u:.4f}")

has_age_u = [(r, a) for r, a in zip(ukb_analytic, ages_u) if a is not None]
if has_age_u:
    out_u2 = [int(sf(r.get('ascvd_combined'))) for r, a in has_age_u]
    age_u2 = [a for r, a in has_age_u]
    sss_u2 = [r['_sss'] for r, a in has_age_u]

    auc_age_u = concordance(out_u2, age_u2)

    ma = sum(age_u2)/len(age_u2)
    sa = (sum((a-ma)**2 for a in age_u2)/(len(age_u2)-1))**0.5
    ms = sum(sss_u2)/len(sss_u2)
    ss = (sum((s-ms)**2 for s in sss_u2)/(len(sss_u2)-1))**0.5
    comb_u = [(a-ma)/max(sa,0.001) + (s-ms)/max(ss,0.001) for a, s in zip(age_u2, sss_u2)]
    auc_comb_u = concordance(out_u2, comb_u)

    print(f"Age alone AUC: {auc_age_u:.4f} (n={len(has_age_u)}, events={sum(out_u2)})")
    print(f"Age + Gene-SSS AUC: {auc_comb_u:.4f}")
    print(f"Delta: +{(auc_comb_u-auc_age_u)*100:.1f}%")

# ═══════════════════════════════════════════════════════════════
# SECTION H: CASCADE vs INDEX — STRUCTURAL DIFFERENCES
# ═══════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("  SECTION H: CASCADE vs INDEX PATIENT STRUCTURAL ANALYSIS")
print("="*80)

# I_Vs_R: 1.0 = Index, 2.0 = Relative/Cascade
index_pts = [r for r in wales if sf(r.get('I_Vs_R')) == 1.0 and r['_sss'] is not None]
cascade_pts = [r for r in wales if sf(r.get('I_Vs_R')) == 2.0 and r['_sss'] is not None]

print(f"Index patients with SSS: {len(index_pts)}")
print(f"Cascade patients with SSS: {len(cascade_pts)}")

group_report("Index (proband)", index_pts)
group_report("Cascade (relative)", cascade_pts)

# SSS comparison
if index_pts and cascade_pts:
    idx_sss = [r['_sss'] for r in index_pts]
    cas_sss = [r['_sss'] for r in cascade_pts]
    idx_ev = sum(1 for r in index_pts if sf(r.get('ascvd_combine')) == 1)
    cas_ev = sum(1 for r in cascade_pts if sf(r.get('ascvd_combine')) == 1)

    print(f"\n  Index SSS: {sum(idx_sss)/len(idx_sss):.3f} +/- {(sum((s-sum(idx_sss)/len(idx_sss))**2 for s in idx_sss)/(len(idx_sss)-1))**0.5:.3f}")
    print(f"  Cascade SSS: {sum(cas_sss)/len(cas_sss):.3f} +/- {(sum((s-sum(cas_sss)/len(cas_sss))**2 for s in cas_sss)/(len(cas_sss)-1))**0.5:.3f}")
    print(f"  Index ASCVD: {idx_ev}/{len(index_pts)} ({idx_ev/len(index_pts)*100:.1f}%)")
    print(f"  Cascade ASCVD: {cas_ev}/{len(cascade_pts)} ({cas_ev/len(cascade_pts)*100:.1f}%)")

# ═══════════════════════════════════════════════════════════════
# SECTION I: COMPREHENSIVE SUMMARY
# ═══════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("  COMPREHENSIVE 3-COHORT VALIDATION SUMMARY")
print("="*80)

print(f"""
┌─────────────────────────────────────────────────────────────────┐
│              3-COHORT EXTERNAL VALIDATION                       │
├─────────────────┬──────────┬────────────┬──────────────────────┤
│ Cohort          │ N        │ Events     │ SSS AUC              │
├─────────────────┼──────────┼────────────┼──────────────────────┤
│ All Wales (Dev) │ {len(wales):>6}   │ {sum(1 for r in wales if sf(r.get('ascvd_combine'))==1):>6}     │ See Section E        │
│ DRAGON3 (Val1)  │ {len(dragon):>6}   │ {sum(1 for r in dragon if sf(r.get('ASCVD_combined'))==1):>6}     │ See Section F        │
│ UKB (Val2)      │ {len(ukb):>6}   │ {sum(1 for r in ukb if sf(r.get('ascvd_combined'))==1):>6}     │ See Section G        │
├─────────────────┼──────────┼────────────┼──────────────────────┤
│ TOTAL           │ {len(wales)+len(dragon)+len(ukb):>6}   │ {sum(1 for r in wales if sf(r.get('ascvd_combine'))==1)+sum(1 for r in dragon if sf(r.get('ASCVD_combined'))==1)+sum(1 for r in ukb if sf(r.get('ascvd_combined'))==1):>6}     │                      │
└─────────────────┴──────────┴────────────┴──────────────────────┘

VUS RECLASSIFICATION:
  Total VUS carriers: {len(vus)}
  VUS with AlphaFold3/FoldX data: {len(vus_with_sss)}
  VUS reclassified as likely pathogenic: {reclassified_path}
  VUS reclassified as likely benign: {reclassified_benign}
""")

print("DONE.")

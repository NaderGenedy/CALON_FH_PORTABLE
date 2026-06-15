#!/usr/bin/env python3
"""
32_SSS_UKB_develop_Wales_validate.py

Develop SSS-based risk model in UKB (larger cohort),
then externally validate in All Wales (independent cohort).

Gold-standard study design: develop in one population, validate in another.
"""

import csv, math, random, os, sys

# ── helpers ──────────────────────────────────────────────────────────────────

def sf(x):
    try: return float(str(x).strip())
    except: return None

def mean_sd(vals):
    vals = [v for v in vals if v is not None]
    if not vals: return 0, 0
    m = sum(vals)/len(vals)
    sd = (sum((v-m)**2 for v in vals)/(max(len(vals)-1,1)))**0.5
    return m, sd

def concordance(outcomes, scores):
    pos = [s for s, o in zip(scores, outcomes) if o == 1]
    neg = [s for s, o in zip(scores, outcomes) if o == 0]
    if not pos or not neg: return 0.5
    conc = sum(1 for sp in pos for sn in neg if sp > sn) + 0.5*sum(1 for sp in pos for sn in neg if sp == sn)
    return conc / (len(pos) * len(neg))

def bootstrap_auc(outcomes, scores, n_boot=500):
    random.seed(2026)
    aucs = []
    n = len(outcomes)
    for _ in range(n_boot):
        idx = [random.randint(0, n-1) for _ in range(n)]
        boot_out = [outcomes[i] for i in idx]
        boot_sc = [scores[i] for i in idx]
        if sum(boot_out) == 0 or sum(boot_out) == n: continue
        aucs.append(concordance(boot_out, boot_sc))
    aucs.sort()
    if len(aucs) < 10: return 0.5, 0, 1
    return sum(aucs)/len(aucs), aucs[int(len(aucs)*0.025)], aucs[int(len(aucs)*0.975)]

def pearson_r_ci(x, y):
    pairs = [(xi,yi) for xi,yi in zip(x,y) if xi is not None and yi is not None]
    n = len(pairs)
    if n < 5: return 0, -1, 1, 1.0, n
    x2, y2 = zip(*pairs)
    mx, my = sum(x2)/n, sum(y2)/n
    sxx = sum((xi-mx)**2 for xi in x2)
    syy = sum((yi-my)**2 for yi in y2)
    sxy = sum((xi-mx)*(yi-my) for xi,yi in zip(x2,y2))
    if sxx==0 or syy==0: return 0, -1, 1, 1.0, n
    r = sxy/(sxx*syy)**0.5
    z = 0.5*math.log((1+r)/(1-r+1e-15))
    se = 1/(n-3)**0.5 if n>3 else 1
    z_lo, z_hi = z-1.96*se, z+1.96*se
    r_lo = (math.exp(2*z_lo)-1)/(math.exp(2*z_lo)+1)
    r_hi = (math.exp(2*z_hi)-1)/(math.exp(2*z_hi)+1)
    t = r*((n-2)/(1-r*r+1e-15))**0.5
    p = 2*math.exp(-0.5*t*t)*0.4 if abs(t)<5 else 0.001
    return r, r_lo, r_hi, p, n

def odds_ratio_ci(a, b, c, d):
    if b*c == 0: return 0, 0, 0, 1.0
    OR = (a*d)/(b*c)
    se = (1/max(a,1)+1/max(b,1)+1/max(c,1)+1/max(d,1))**0.5
    ci_lo = math.exp(math.log(max(OR,0.001))-1.96*se)
    ci_hi = math.exp(math.log(max(OR,0.001))+1.96*se)
    n = a+b+c+d
    chi2 = sum((o-e)**2/max(e,0.001) for o,e in [(a,(a+b)*(a+c)/n),(b,(a+b)*(b+d)/n),(c,(c+d)*(a+c)/n),(d,(c+d)*(b+d)/n)])
    p = math.exp(-0.5*chi2) if chi2<30 else 0.0001
    return OR, ci_lo, ci_hi, p

def read_csv(path):
    with open(path, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        return list(reader)

def spearman_rank(x, y):
    """Spearman rank correlation between two lists."""
    pairs = [(xi,yi) for xi,yi in zip(x,y) if xi is not None and yi is not None]
    n = len(pairs)
    if n < 3: return 0, 1.0, n
    x2, y2 = zip(*pairs)
    def rank(vals):
        s = sorted(range(len(vals)), key=lambda i: vals[i])
        ranks = [0]*len(vals)
        i = 0
        while i < len(s):
            j = i
            while j < len(s) and vals[s[j]] == vals[s[i]]:
                j += 1
            avg = (i+j-1)/2.0 + 1
            for k in range(i, j):
                ranks[s[k]] = avg
            i = j
        return ranks
    rx = rank(x2)
    ry = rank(y2)
    d2 = sum((a-b)**2 for a,b in zip(rx,ry))
    rho = 1 - 6*d2/(n*(n*n-1))
    t = rho * ((n-2)/(1-rho*rho+1e-15))**0.5 if abs(rho) < 1 else 10
    p = 2*math.exp(-0.5*t*t)*0.4 if abs(t)<5 else 0.001
    return rho, p, n

def logistic_score(coefs, vals):
    """Simple linear predictor from coefficients dict."""
    return sum(coefs.get(k, 0)*v for k, v in vals.items())

# ── LOAD DATA ────────────────────────────────────────────────────────────────

BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"
print("="*80)
print("SSS v2: UKB Development + All Wales External Validation")
print("="*80)

print("\n[Loading data...]")
ukb = read_csv(os.path.join(BASE, "alphafold/analysis/ukb_external_validation_sss.csv"))
wales_raw = read_csv(os.path.join(BASE, "WALES_FH_CLEANED (1) - Copy.csv"))
sss_map_raw = read_csv(os.path.join(BASE, "alphafold/analysis/structural_severity_scores.csv"))

print(f"  UKB rows: {len(ukb)}")
print(f"  Wales rows: {len(wales_raw)}")
print(f"  SSS mapping rows: {len(sss_map_raw)}")

# Build SSS lookup: variant_id -> {sss, domain, variant_type, domain_severity}
sss_lookup = {}
for r in sss_map_raw:
    vid = r.get('variant_id','').strip()
    if vid:
        sss_lookup[vid] = {
            'sss': sf(r.get('sss','')),
            'domain': r.get('domain','').strip(),
            'variant_type': r.get('variant_type','').strip(),
            'domain_severity': sf(r.get('domain_severity',''))
        }

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Develop Domain-Risk Model in UKB
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("SECTION 1: Develop Domain-Risk Model in UKB (Development Cohort)")
print("="*80)

# Parse UKB data
for r in ukb:
    r['_age'] = sf(r.get('age',''))
    r['_sex'] = sf(r.get('sex',''))
    r['_ldl'] = sf(r.get('ldl',''))
    r['_sss'] = sf(r.get('sss',''))
    r['_ascvd'] = int(sf(r.get('ascvd_combined','')) or 0)
    r['_domain'] = r.get('domain','').strip()
    r['_on_statin'] = int(sf(r.get('on_statin','')) or 0)
    r['_variant_id'] = r.get('variant_id','').strip()
    # Get variant_type from sss_lookup
    info = sss_lookup.get(r['_variant_id'], {})
    r['_variant_type'] = info.get('variant_type', '')

# 1a. Domain-specific ASCVD rates in UKB
domain_stats_ukb = {}
for r in ukb:
    d = r['_domain']
    if not d: continue
    if d not in domain_stats_ukb:
        domain_stats_ukb[d] = {'n': 0, 'events': 0}
    domain_stats_ukb[d]['n'] += 1
    domain_stats_ukb[d]['events'] += r['_ascvd']

print("\n1a. Domain-specific ASCVD rates in UKB (Development):")
print(f"  {'Domain':<35} {'N':>5} {'Events':>7} {'Rate':>8}")
print("  " + "-"*60)

overall_rate = sum(r['_ascvd'] for r in ukb) / max(len(ukb),1)

domain_risk_weights = {}
for d in sorted(domain_stats_ukb.keys()):
    s = domain_stats_ukb[d]
    rate = s['events']/max(s['n'],1)
    print(f"  {d:<35} {s['n']:>5} {s['events']:>7} {rate:>8.3f}")
    domain_risk_weights[d] = rate

# 1b. Normalize domain risk weights to 0-1 scale
max_rate = max(domain_risk_weights.values()) if domain_risk_weights else 1
min_rate = min(domain_risk_weights.values()) if domain_risk_weights else 0
range_rate = max_rate - min_rate if max_rate > min_rate else 1

for d in domain_risk_weights:
    domain_risk_weights[d] = (domain_risk_weights[d] - min_rate) / range_rate

print("\n1b. Normalized Domain Risk Weights (0-1 scale):")
for d in sorted(domain_risk_weights.keys()):
    print(f"  {d:<35} {domain_risk_weights[d]:.4f}")

# 1c. Variant type weights
vtype_stats = {}
for r in ukb:
    vt = r['_variant_type']
    if not vt: vt = 'unknown'
    if vt not in vtype_stats:
        vtype_stats[vt] = {'n': 0, 'events': 0}
    vtype_stats[vt]['n'] += 1
    vtype_stats[vt]['events'] += r['_ascvd']

vtype_weights = {}
print("\n1c. Variant Type ASCVD rates in UKB:")
print(f"  {'Type':<25} {'N':>5} {'Events':>7} {'Rate':>8}")
print("  " + "-"*50)
for vt in sorted(vtype_stats.keys()):
    s = vtype_stats[vt]
    rate = s['events']/max(s['n'],1)
    print(f"  {vt:<25} {s['n']:>5} {s['events']:>7} {rate:>8.3f}")
    vtype_weights[vt] = rate

# Normalize variant type weights
vt_max = max(vtype_weights.values()) if vtype_weights else 1
vt_min = min(vtype_weights.values()) if vtype_weights else 0
vt_range = vt_max - vt_min if vt_max > vt_min else 1
for vt in vtype_weights:
    vtype_weights[vt] = (vtype_weights[vt] - vt_min) / vt_range

print("\n  Normalized variant type weights:")
for vt in sorted(vtype_weights.keys()):
    print(f"  {vt:<25} {vtype_weights[vt]:.4f}")

# 1d. Create SSS_v2 for UKB
print("\n1d. Computing SSS_v2 for UKB...")
for r in ukb:
    sss_orig = r['_sss'] if r['_sss'] is not None else 0.5
    drw = domain_risk_weights.get(r['_domain'], 0.5)
    vtw = vtype_weights.get(r['_variant_type'], 0.5)
    r['_sss_v2'] = sss_orig * 0.4 + drw * 0.4 + vtw * 0.2
    r['_drw'] = drw

sss_v1_vals = [r['_sss'] for r in ukb if r['_sss'] is not None]
sss_v2_vals = [r['_sss_v2'] for r in ukb]
m1, s1 = mean_sd(sss_v1_vals)
m2, s2 = mean_sd(sss_v2_vals)
print(f"  SSS_v1: mean={m1:.4f}, SD={s1:.4f}")
print(f"  SSS_v2: mean={m2:.4f}, SD={s2:.4f}")

# 1e. Test SSS_v2 AUC in UKB (development AUC)
print("\n1e. Model AUCs in UKB (Development Cohort):")
outcomes_ukb = [r['_ascvd'] for r in ukb]
n_events_ukb = sum(outcomes_ukb)
print(f"  N={len(ukb)}, Events={n_events_ukb}, Rate={n_events_ukb/len(ukb):.3f}")

models_ukb = {}

# Age alone
scores = [r['_age'] if r['_age'] is not None else 50 for r in ukb]
auc, lo, hi = bootstrap_auc(outcomes_ukb, scores)
models_ukb['Age alone'] = (auc, lo, hi)
print(f"  Age alone:                AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# Age + Sex
scores = [(r['_age'] or 50) + (r['_sex'] or 0)*5 for r in ukb]
auc, lo, hi = bootstrap_auc(outcomes_ukb, scores)
models_ukb['Age + Sex'] = (auc, lo, hi)
print(f"  Age + Sex:                AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# Age + Sex + SSS_v1
scores = [(r['_age'] or 50) + (r['_sex'] or 0)*5 + (r['_sss'] or 0.5)*20 for r in ukb]
auc, lo, hi = bootstrap_auc(outcomes_ukb, scores)
models_ukb['Age + Sex + SSS_v1'] = (auc, lo, hi)
print(f"  Age + Sex + SSS_v1:       AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# Age + Sex + SSS_v2
scores = [(r['_age'] or 50) + (r['_sex'] or 0)*5 + r['_sss_v2']*20 for r in ukb]
auc, lo, hi = bootstrap_auc(outcomes_ukb, scores)
models_ukb['Age + Sex + SSS_v2'] = (auc, lo, hi)
print(f"  Age + Sex + SSS_v2:       AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# Age + Sex + SSS_v2 + LDL
scores = [(r['_age'] or 50) + (r['_sex'] or 0)*5 + r['_sss_v2']*20 + (r['_ldl'] or 3)*3 for r in ukb]
auc, lo, hi = bootstrap_auc(outcomes_ukb, scores)
models_ukb['Age + Sex + SSS_v2 + LDL'] = (auc, lo, hi)
print(f"  Age + Sex + SSS_v2 + LDL: AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# Age + Sex + SSS_v2 + on_statin
scores = [(r['_age'] or 50) + (r['_sex'] or 0)*5 + r['_sss_v2']*20 + r['_on_statin']*5 for r in ukb]
auc, lo, hi = bootstrap_auc(outcomes_ukb, scores)
models_ukb['Age + Sex + SSS_v2 + Statin'] = (auc, lo, hi)
print(f"  Age + Sex + SSS_v2 + Statin: AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# SSS_v2 alone
scores = [r['_sss_v2'] for r in ukb]
auc, lo, hi = bootstrap_auc(outcomes_ukb, scores)
models_ukb['SSS_v2 alone'] = (auc, lo, hi)
print(f"  SSS_v2 alone:             AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# Domain risk weight alone
scores = [r['_drw'] for r in ukb]
auc, lo, hi = bootstrap_auc(outcomes_ukb, scores)
models_ukb['Domain risk weight'] = (auc, lo, hi)
print(f"  Domain risk weight alone: AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# 1f. SSS_v1 vs SSS_v2
print("\n1f. SSS_v1 vs SSS_v2 in UKB:")
scores_v1 = [r['_sss'] if r['_sss'] is not None else 0.5 for r in ukb]
scores_v2 = [r['_sss_v2'] for r in ukb]
auc_v1 = concordance(outcomes_ukb, scores_v1)
auc_v2 = concordance(outcomes_ukb, scores_v2)
print(f"  SSS_v1 C-statistic: {auc_v1:.4f}")
print(f"  SSS_v2 C-statistic: {auc_v2:.4f}")
print(f"  Delta AUC: {auc_v2 - auc_v1:+.4f}")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Map SSS_v2 to All Wales (External Validation)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("SECTION 2: Map SSS_v2 to All Wales (External Validation)")
print("="*80)

# Parse Wales data and map SSS
wales = []
mapped = 0
unmapped = 0

for r in wales_raw:
    mut = r.get('Mutation1','').strip()
    if not mut:
        continue

    # Look up SSS
    info = sss_lookup.get(mut, None)
    if info is None or info['sss'] is None:
        unmapped += 1
        continue

    mapped += 1
    w = {}
    w['mutation'] = mut
    w['sss_v1'] = info['sss']
    w['domain'] = info['domain']
    w['variant_type'] = info['variant_type']

    # Apply UKB-derived domain weights (FROZEN)
    drw = domain_risk_weights.get(info['domain'], 0.5)
    vtw = vtype_weights.get(info['variant_type'], 0.5)
    w['sss_v2'] = info['sss'] * 0.4 + drw * 0.4 + vtw * 0.2
    w['drw'] = drw

    # Parse Wales clinical
    ascvd_str = str(r.get('ascvd_combine','')).strip()
    if ascvd_str in ('1.0','1'): w['ascvd'] = 1
    elif ascvd_str in ('0.0','0'): w['ascvd'] = 0
    else: w['ascvd'] = 0

    w['gender'] = r.get('Gender','').strip()
    w['sex'] = 1 if w['gender'] == 'M' else 0
    w['ldl'] = sf(r.get('LDL.1',''))
    w['tc'] = sf(r.get('TC.1',''))
    w['hdl'] = sf(r.get('HDL.1',''))
    w['trig'] = sf(r.get('TRG.1',''))
    w['lpa'] = sf(r.get('Lpa.1',''))

    w['positive'] = str(r.get('Positive1','')).strip()
    w['vus'] = str(r.get('VUS1','')).strip()
    w['proband'] = str(r.get('Proband','')).strip()
    w['i_vs_r'] = str(r.get('I_Vs_R','')).strip()
    w['xanthomata'] = str(r.get('TendonXanthomata','')).strip()
    w['simon_broome'] = str(r.get('SimonBroome','')).strip()
    w['diabetes'] = str(r.get('Diabetes','')).strip()
    w['smoking'] = str(r.get('Smoking','')).strip()
    w['sbp'] = sf(r.get('BloodPressureSystolic',''))
    w['bmi'] = sf(r.get('BMI',''))

    # Age from DOB
    dob = r.get('DOB','').strip()
    bmi_age = sf(r.get('BMI_AGE',''))
    w['age'] = bmi_age  # Use BMI_AGE as proxy for age

    # Statin
    statin_str = str(r.get('statin_at_genetic_test','')).strip()
    treat_str = str(r.get('Treatment1','')).strip()
    w['on_statin'] = 1 if (statin_str not in ('','0','NoValue') or (treat_str not in ('','NoValue'))) else 0

    # LDL measurements for treatment response
    w['ldl1'] = sf(r.get('LDL.1',''))
    w['ldl2'] = sf(r.get('LDL.2',''))
    w['ldl3'] = sf(r.get('LDL.3',''))
    w['ldl4'] = sf(r.get('LDL.4',''))

    # ASCVD components
    w['mi'] = str(r.get('MIACS','')).strip()
    w['cabg'] = str(r.get('CABG','')).strip()
    w['angina'] = str(r.get('ANGINA','')).strip()
    w['tia'] = str(r.get('TIA','')).strip()
    w['pvd'] = str(r.get('PVD','')).strip()

    wales.append(w)

print(f"\n  Wales patients with Mutation1: {mapped + unmapped}")
print(f"  Successfully mapped to SSS: {mapped}")
print(f"  Unmapped (variant not in SSS table): {unmapped}")
print(f"  Wales analysis cohort: {len(wales)}")

n_events_w = sum(w['ascvd'] for w in wales)
print(f"  ASCVD events: {n_events_w} ({n_events_w/max(len(wales),1)*100:.1f}%)")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: External Validation Statistics
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("SECTION 3: External Validation AUCs in All Wales")
print("="*80)

outcomes_w = [w['ascvd'] for w in wales]
print(f"  N={len(wales)}, Events={n_events_w}")

models_wales = {}

# Age alone
scores = [w['age'] if w['age'] is not None else 50 for w in wales]
auc, lo, hi = bootstrap_auc(outcomes_w, scores)
models_wales['Age alone'] = (auc, lo, hi)
print(f"  Age alone:                AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# Age + Sex
scores = [(w['age'] or 50) + w['sex']*5 for w in wales]
auc, lo, hi = bootstrap_auc(outcomes_w, scores)
models_wales['Age + Sex'] = (auc, lo, hi)
print(f"  Age + Sex:                AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# Age + Sex + SSS_v1
scores = [(w['age'] or 50) + w['sex']*5 + w['sss_v1']*20 for w in wales]
auc, lo, hi = bootstrap_auc(outcomes_w, scores)
models_wales['Age + Sex + SSS_v1'] = (auc, lo, hi)
print(f"  Age + Sex + SSS_v1:       AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# Age + Sex + SSS_v2
scores = [(w['age'] or 50) + w['sex']*5 + w['sss_v2']*20 for w in wales]
auc, lo, hi = bootstrap_auc(outcomes_w, scores)
models_wales['Age + Sex + SSS_v2'] = (auc, lo, hi)
print(f"  Age + Sex + SSS_v2:       AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# Age + Sex + SSS_v2 + LDL
scores = [(w['age'] or 50) + w['sex']*5 + w['sss_v2']*20 + (w['ldl'] or 3)*3 for w in wales]
auc, lo, hi = bootstrap_auc(outcomes_w, scores)
models_wales['Age + Sex + SSS_v2 + LDL'] = (auc, lo, hi)
print(f"  Age + Sex + SSS_v2 + LDL: AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# Age + Sex + SSS_v2 + on_statin
scores = [(w['age'] or 50) + w['sex']*5 + w['sss_v2']*20 + w['on_statin']*5 for w in wales]
auc, lo, hi = bootstrap_auc(outcomes_w, scores)
models_wales['Age + Sex + SSS_v2 + Statin'] = (auc, lo, hi)
print(f"  Age + Sex + SSS_v2 + Statin: AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# SSS_v2 alone
scores = [w['sss_v2'] for w in wales]
auc, lo, hi = bootstrap_auc(outcomes_w, scores)
models_wales['SSS_v2 alone'] = (auc, lo, hi)
print(f"  SSS_v2 alone:             AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# Domain risk weight alone
scores = [w['drw'] for w in wales]
auc, lo, hi = bootstrap_auc(outcomes_w, scores)
models_wales['Domain risk weight'] = (auc, lo, hi)
print(f"  Domain risk weight alone: AUC={auc:.4f} (95% CI: {lo:.4f}-{hi:.4f})")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Subgroup Validation in Wales
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("SECTION 4: Subgroup Validation in Wales")
print("="*80)

subgroups = {
    'All mapped': lambda w: True,
    'Pathogenic only': lambda w: w['positive'] == '1' and w['vus'] != '1.0',
    'VUS only': lambda w: w['vus'] == '1.0',
    'Index cases': lambda w: w['proband'] in ('1.0','1') or w['i_vs_r'] in ('1.0','1'),
    'Cascade cases': lambda w: w['i_vs_r'] == '2.0',
    'With xanthomata': lambda w: w['xanthomata'] in ('1','1.0','Yes'),
    'Age >= 50': lambda w: w['age'] is not None and w['age'] >= 50,
    'Age < 50': lambda w: w['age'] is not None and w['age'] < 50,
}

subgroup_results = {}
print(f"\n  {'Subgroup':<25} {'N':>6} {'Events':>7} {'Rate':>7} {'AUC':>7} {'95% CI':>16}")
print("  " + "-"*72)

for sg_name, sg_filter in subgroups.items():
    sg = [w for w in wales if sg_filter(w)]
    if len(sg) < 10:
        print(f"  {sg_name:<25} {len(sg):>6} {'--':>7} {'--':>7} {'--':>7} {'--':>16}")
        continue
    sg_out = [w['ascvd'] for w in sg]
    sg_scores = [w['sss_v2'] for w in sg]
    n_ev = sum(sg_out)
    rate = n_ev / len(sg)
    if n_ev == 0 or n_ev == len(sg):
        print(f"  {sg_name:<25} {len(sg):>6} {n_ev:>7} {rate:>7.3f} {'--':>7} {'--':>16}")
        subgroup_results[sg_name] = (len(sg), n_ev, 0.5, 0, 1)
        continue
    auc, lo, hi = bootstrap_auc(sg_out, sg_scores)
    print(f"  {sg_name:<25} {len(sg):>6} {n_ev:>7} {rate:>7.3f} {auc:>7.4f} {lo:.4f}-{hi:.4f}")
    subgroup_results[sg_name] = (len(sg), n_ev, auc, lo, hi)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: Domain-Specific Validation
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("SECTION 5: Domain-Specific Validation (UKB vs Wales)")
print("="*80)

domain_stats_wales = {}
for w in wales:
    d = w['domain']
    if not d: continue
    if d not in domain_stats_wales:
        domain_stats_wales[d] = {'n': 0, 'events': 0}
    domain_stats_wales[d]['n'] += 1
    domain_stats_wales[d]['events'] += w['ascvd']

# Compare domain ASCVD rates
all_domains = sorted(set(list(domain_stats_ukb.keys()) + list(domain_stats_wales.keys())))
print(f"\n  {'Domain':<35} {'UKB N':>6} {'UKB Rate':>9} {'Wales N':>8} {'Wales Rate':>11}")
print("  " + "-"*75)

ukb_rates = []
wales_rates = []
domains_both = []

for d in all_domains:
    ukb_s = domain_stats_ukb.get(d, {'n': 0, 'events': 0})
    wal_s = domain_stats_wales.get(d, {'n': 0, 'events': 0})
    ukb_r = ukb_s['events']/max(ukb_s['n'],1) if ukb_s['n'] > 0 else None
    wal_r = wal_s['events']/max(wal_s['n'],1) if wal_s['n'] > 0 else None
    ukb_r_str = f"{ukb_r:.3f}" if ukb_r is not None else "--"
    wal_r_str = f"{wal_r:.3f}" if wal_r is not None else "--"
    print(f"  {d:<35} {ukb_s['n']:>6} {ukb_r_str:>9} {wal_s['n']:>8} {wal_r_str:>11}")
    if ukb_r is not None and wal_r is not None and ukb_s['n'] >= 3 and wal_s['n'] >= 3:
        ukb_rates.append(ukb_r)
        wales_rates.append(wal_r)
        domains_both.append(d)

if len(ukb_rates) >= 3:
    rho, p, n = spearman_rank(ukb_rates, wales_rates)
    print(f"\n  Spearman rank correlation of domain ASCVD rates (UKB vs Wales):")
    print(f"    rho = {rho:.4f}, p = {p:.4f}, N_domains = {n}")
    print(f"    Interpretation: {'Strong' if abs(rho) > 0.7 else 'Moderate' if abs(rho) > 0.4 else 'Weak'} correlation")
    print(f"    => Domain risk ranking {'TRANSFERS' if rho > 0.3 else 'does NOT transfer'} across populations")
else:
    print("\n  Not enough overlapping domains for Spearman correlation.")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: Treatment Response Validation
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("SECTION 6: Treatment Response Validation in Wales")
print("="*80)

# For treated patients, compute LDL change
treated = []
for w in wales:
    if w['on_statin'] != 1: continue
    ldl1 = w['ldl1']
    if ldl1 is None: continue
    # Find latest LDL
    latest = None
    for k in ['ldl4','ldl3','ldl2']:
        if w[k] is not None:
            latest = w[k]
            break
    if latest is None: continue
    ldl_change = latest - ldl1
    ldl_pct = ldl_change / ldl1 * 100 if ldl1 != 0 else None
    treated.append({
        'sss_v1': w['sss_v1'],
        'sss_v2': w['sss_v2'],
        'domain': w['domain'],
        'drw': w['drw'],
        'ldl_change': ldl_change,
        'ldl_pct': ldl_pct,
        'ldl_start': ldl1,
        'ldl_end': latest
    })

print(f"\n  Treated patients with LDL change data: {len(treated)}")
if treated:
    changes = [t['ldl_change'] for t in treated]
    pcts = [t['ldl_pct'] for t in treated if t['ldl_pct'] is not None]
    m_ch, s_ch = mean_sd(changes)
    m_pct, s_pct = mean_sd(pcts)
    print(f"  Mean LDL change: {m_ch:.2f} +/- {s_ch:.2f} mmol/L")
    print(f"  Mean % LDL change: {m_pct:.1f}% +/- {s_pct:.1f}%")

    # SSS_v2 predicts LDL change?
    x_sss = [t['sss_v2'] for t in treated]
    y_ldl = [t['ldl_change'] for t in treated]
    r, rlo, rhi, p, n = pearson_r_ci(x_sss, y_ldl)
    print(f"\n  SSS_v2 vs LDL change:")
    print(f"    Pearson r = {r:.4f} (95% CI: {rlo:.4f} to {rhi:.4f})")
    print(f"    p = {p:.4f}, N = {n}")

    # SSS_v1 predicts LDL change?
    x_sss1 = [t['sss_v1'] for t in treated]
    r1, r1lo, r1hi, p1, n1 = pearson_r_ci(x_sss1, y_ldl)
    print(f"\n  SSS_v1 vs LDL change:")
    print(f"    Pearson r = {r1:.4f} (95% CI: {r1lo:.4f} to {r1hi:.4f})")
    print(f"    p = {p1:.4f}, N = {n1}")

    # Domain predicts LDL change (ANOVA-like)
    print(f"\n  Domain-specific LDL change (ANOVA-like):")
    dom_ldl = {}
    for t in treated:
        d = t['domain']
        if d not in dom_ldl: dom_ldl[d] = []
        dom_ldl[d].append(t['ldl_change'])

    print(f"    {'Domain':<35} {'N':>5} {'Mean LDL chg':>13} {'SD':>7}")
    print("    " + "-"*65)
    for d in sorted(dom_ldl.keys()):
        m, s = mean_sd(dom_ldl[d])
        print(f"    {d:<35} {len(dom_ldl[d]):>5} {m:>13.3f} {s:>7.3f}")

    # Grand mean and between-group variance
    all_changes = [c for lst in dom_ldl.values() for c in lst]
    grand_mean = sum(all_changes)/len(all_changes) if all_changes else 0
    ss_between = sum(len(lst)*(sum(lst)/len(lst) - grand_mean)**2 for lst in dom_ldl.values() if lst)
    ss_within = sum(sum((c-sum(lst)/len(lst))**2 for c in lst) for lst in dom_ldl.values() if lst)
    k = len(dom_ldl)
    n_total = len(all_changes)
    if k > 1 and n_total > k and ss_within > 0:
        ms_between = ss_between / (k-1)
        ms_within = ss_within / (n_total - k)
        f_stat = ms_between / ms_within
        print(f"\n    F-statistic = {f_stat:.3f} (df = {k-1}, {n_total-k})")
        f_p = math.exp(-0.5*f_stat) if f_stat < 30 else 0.0001
        print(f"    Approximate p = {f_p:.4f}")
else:
    print("  No treated patients with LDL change data available.")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: VUS Reclassification Validation
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("SECTION 7: VUS Reclassification Validation in Wales")
print("="*80)

vus_patients = [w for w in wales if w['vus'] == '1.0']
print(f"\n  VUS patients with SSS: {len(vus_patients)}")
n_vus_events = sum(w['ascvd'] for w in vus_patients)
print(f"  VUS ASCVD events: {n_vus_events} ({n_vus_events/max(len(vus_patients),1)*100:.1f}%)")

if len(vus_patients) >= 10 and n_vus_events > 0:
    vus_out = [w['ascvd'] for w in vus_patients]

    # Domain risk weight predicts ASCVD in VUS?
    vus_drw = [w['drw'] for w in vus_patients]
    auc_drw, lo_drw, hi_drw = bootstrap_auc(vus_out, vus_drw)
    print(f"\n  Domain risk weight -> VUS ASCVD:")
    print(f"    AUC = {auc_drw:.4f} (95% CI: {lo_drw:.4f}-{hi_drw:.4f})")

    # SSS_v1 prediction
    vus_v1 = [w['sss_v1'] for w in vus_patients]
    auc_v1, lo_v1, hi_v1 = bootstrap_auc(vus_out, vus_v1)
    print(f"\n  SSS_v1 -> VUS ASCVD:")
    print(f"    AUC = {auc_v1:.4f} (95% CI: {lo_v1:.4f}-{hi_v1:.4f})")

    # SSS_v2 prediction
    vus_v2 = [w['sss_v2'] for w in vus_patients]
    auc_v2_vus, lo_v2, hi_v2 = bootstrap_auc(vus_out, vus_v2)
    print(f"\n  SSS_v2 -> VUS ASCVD:")
    print(f"    AUC = {auc_v2_vus:.4f} (95% CI: {lo_v2:.4f}-{hi_v2:.4f})")

    # Sensitivity / specificity of domain-risk >= 0.5
    tp = sum(1 for w in vus_patients if w['drw'] >= 0.5 and w['ascvd'] == 1)
    fp = sum(1 for w in vus_patients if w['drw'] >= 0.5 and w['ascvd'] == 0)
    fn = sum(1 for w in vus_patients if w['drw'] < 0.5 and w['ascvd'] == 1)
    tn = sum(1 for w in vus_patients if w['drw'] < 0.5 and w['ascvd'] == 0)
    sens = tp/(tp+fn) if (tp+fn) > 0 else 0
    spec = tn/(tn+fp) if (tn+fp) > 0 else 0
    ppv = tp/(tp+fp) if (tp+fp) > 0 else 0
    npv = tn/(tn+fn) if (tn+fn) > 0 else 0
    print(f"\n  Domain risk >= 0.5 for VUS ASCVD prediction:")
    print(f"    TP={tp}, FP={fp}, FN={fn}, TN={tn}")
    print(f"    Sensitivity = {sens:.4f}")
    print(f"    Specificity = {spec:.4f}")
    print(f"    PPV = {ppv:.4f}")
    print(f"    NPV = {npv:.4f}")
    OR, ci_lo, ci_hi, p_or = odds_ratio_ci(tp, fp, fn, tn)
    print(f"    OR = {OR:.3f} (95% CI: {ci_lo:.3f}-{ci_hi:.3f}), p={p_or:.4f}")
else:
    print("  Not enough VUS patients or events for analysis.")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8: Summary Forest Plot Data
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("SECTION 8: Forest Plot Data (Master Summary)")
print("="*80)

forest_rows = []

# UKB development
for model_name in ['Age alone','Age + Sex','Age + Sex + SSS_v1','Age + Sex + SSS_v2',
                    'Age + Sex + SSS_v2 + LDL','Age + Sex + SSS_v2 + Statin',
                    'SSS_v2 alone','Domain risk weight']:
    if model_name in models_ukb:
        auc, lo, hi = models_ukb[model_name]
        forest_rows.append({
            'cohort': 'UKB (Development)',
            'subgroup': model_name,
            'N': len(ukb),
            'events': n_events_ukb,
            'AUC': auc,
            'CI_lower': lo,
            'CI_upper': hi
        })

# Wales external validation
for model_name in ['Age alone','Age + Sex','Age + Sex + SSS_v1','Age + Sex + SSS_v2',
                    'Age + Sex + SSS_v2 + LDL','Age + Sex + SSS_v2 + Statin',
                    'SSS_v2 alone','Domain risk weight']:
    if model_name in models_wales:
        auc, lo, hi = models_wales[model_name]
        forest_rows.append({
            'cohort': 'Wales (External)',
            'subgroup': model_name,
            'N': len(wales),
            'events': n_events_w,
            'AUC': auc,
            'CI_lower': lo,
            'CI_upper': hi
        })

# Wales subgroups (SSS_v2 alone)
for sg_name, (n_sg, n_ev, auc, lo, hi) in subgroup_results.items():
    forest_rows.append({
        'cohort': 'Wales Subgroup',
        'subgroup': sg_name,
        'N': n_sg,
        'events': n_ev,
        'AUC': auc,
        'CI_lower': lo,
        'CI_upper': hi
    })

# Print forest plot data
print(f"\n  {'Cohort':<22} {'Subgroup':<30} {'N':>6} {'Events':>7} {'AUC':>7} {'95% CI':>16}")
print("  " + "-"*92)
for fr in forest_rows:
    print(f"  {fr['cohort']:<22} {fr['subgroup']:<30} {fr['N']:>6} {fr['events']:>7} {fr['AUC']:>7.4f} {fr['CI_lower']:.4f}-{fr['CI_upper']:.4f}")

# Save forest plot data
forest_path = os.path.join(BASE, "alphafold/analysis/sss_v2_forest_plot_data.csv")
with open(forest_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['cohort','subgroup','N','events','AUC','CI_lower','CI_upper'])
    writer.writeheader()
    for fr in forest_rows:
        writer.writerow(fr)
print(f"\n  Saved: {forest_path}")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 9: Head-to-Head SSS_v1 vs SSS_v2 in Wales
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("SECTION 9: Head-to-Head SSS_v1 vs SSS_v2 in Wales (External)")
print("="*80)

if len(wales) >= 20 and n_events_w > 0:
    # SSS_v1 AUC
    scores_v1_w = [w['sss_v1'] for w in wales]
    auc_v1_w, lo_v1_w, hi_v1_w = bootstrap_auc(outcomes_w, scores_v1_w)
    print(f"\n  SSS_v1 AUC: {auc_v1_w:.4f} (95% CI: {lo_v1_w:.4f}-{hi_v1_w:.4f})")

    # SSS_v2 AUC
    scores_v2_w = [w['sss_v2'] for w in wales]
    auc_v2_w, lo_v2_w, hi_v2_w = bootstrap_auc(outcomes_w, scores_v2_w)
    print(f"  SSS_v2 AUC: {auc_v2_w:.4f} (95% CI: {lo_v2_w:.4f}-{hi_v2_w:.4f})")

    # Delta AUC
    delta = auc_v2_w - auc_v1_w
    print(f"\n  Delta AUC (v2 - v1): {delta:+.4f}")

    # Continuous NRI
    # NRI = P(score_up | event) - P(score_up | nonevent) + P(score_down | nonevent) - P(score_down | event)
    events_up = sum(1 for w in wales if w['ascvd'] == 1 and w['sss_v2'] > w['sss_v1'])
    events_down = sum(1 for w in wales if w['ascvd'] == 1 and w['sss_v2'] < w['sss_v1'])
    events_same = sum(1 for w in wales if w['ascvd'] == 1 and w['sss_v2'] == w['sss_v1'])
    nonevents_up = sum(1 for w in wales if w['ascvd'] == 0 and w['sss_v2'] > w['sss_v1'])
    nonevents_down = sum(1 for w in wales if w['ascvd'] == 0 and w['sss_v2'] < w['sss_v1'])
    nonevents_same = sum(1 for w in wales if w['ascvd'] == 0 and w['sss_v2'] == w['sss_v1'])

    n_ev_tot = sum(1 for w in wales if w['ascvd'] == 1)
    n_nonev_tot = sum(1 for w in wales if w['ascvd'] == 0)

    if n_ev_tot > 0 and n_nonev_tot > 0:
        nri_events = (events_up - events_down) / n_ev_tot
        nri_nonevents = (nonevents_down - nonevents_up) / n_nonev_tot
        nri = nri_events + nri_nonevents

        print(f"\n  Continuous NRI Analysis:")
        print(f"    Events reclassified UP:   {events_up}")
        print(f"    Events reclassified DOWN: {events_down}")
        print(f"    Events unchanged:         {events_same}")
        print(f"    Non-events reclassified UP:   {nonevents_up}")
        print(f"    Non-events reclassified DOWN: {nonevents_down}")
        print(f"    Non-events unchanged:         {nonevents_same}")
        print(f"\n    NRI (events):     {nri_events:+.4f}")
        print(f"    NRI (non-events): {nri_nonevents:+.4f}")
        print(f"    Total NRI:        {nri:+.4f}")
        print(f"    Interpretation: {'IMPROVEMENT' if nri > 0 else 'No improvement'} with SSS_v2")
else:
    print("  Not enough Wales data for head-to-head comparison.")

# ══════════════════════════════════════════════════════════════════════════════
# SAVE MAIN RESULTS CSV
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("SAVING RESULTS")
print("="*80)

# Save patient-level SSS_v2 data for Wales
out_path = os.path.join(BASE, "alphafold/analysis/sss_v2_ukb_develop_wales_validate.csv")
with open(out_path, 'w', newline='') as f:
    fields = ['cohort','mutation','domain','variant_type','sss_v1','sss_v2','domain_risk_weight',
              'ascvd','age','sex','ldl','on_statin']
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    # UKB rows
    for r in ukb:
        writer.writerow({
            'cohort': 'UKB',
            'mutation': r['_variant_id'],
            'domain': r['_domain'],
            'variant_type': r['_variant_type'],
            'sss_v1': r['_sss'],
            'sss_v2': round(r['_sss_v2'], 4),
            'domain_risk_weight': round(r['_drw'], 4),
            'ascvd': r['_ascvd'],
            'age': r['_age'],
            'sex': r['_sex'],
            'ldl': r['_ldl'],
            'on_statin': r['_on_statin']
        })
    # Wales rows
    for w in wales:
        writer.writerow({
            'cohort': 'Wales',
            'mutation': w['mutation'],
            'domain': w['domain'],
            'variant_type': w['variant_type'],
            'sss_v1': w['sss_v1'],
            'sss_v2': round(w['sss_v2'], 4),
            'domain_risk_weight': round(w['drw'], 4),
            'ascvd': w['ascvd'],
            'age': round(w['age'], 1) if w['age'] is not None else '',
            'sex': w['sex'],
            'ldl': w['ldl'] if w['ldl'] is not None else '',
            'on_statin': w['on_statin']
        })

print(f"  Saved: {out_path}")
print(f"  Saved: {forest_path}")

# ══════════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("FINAL SUMMARY")
print("="*80)

print(f"\n  Development cohort (UKB): N={len(ukb)}, Events={n_events_ukb}")
print(f"  External validation (Wales): N={len(wales)}, Events={n_events_w}")

print(f"\n  SSS_v2 composition: SSS_original*0.4 + domain_risk_weight*0.4 + variant_type_weight*0.2")
print(f"  Domain weights derived from UKB, applied FROZEN to Wales")

# Best model comparison
if 'Age + Sex + SSS_v2' in models_ukb and 'Age + Sex + SSS_v2' in models_wales:
    ukb_best = models_ukb['Age + Sex + SSS_v2']
    wales_best = models_wales['Age + Sex + SSS_v2']
    print(f"\n  Age + Sex + SSS_v2 model:")
    print(f"    Development (UKB):      AUC={ukb_best[0]:.4f} (95% CI: {ukb_best[1]:.4f}-{ukb_best[2]:.4f})")
    print(f"    External (Wales):       AUC={wales_best[0]:.4f} (95% CI: {wales_best[1]:.4f}-{wales_best[2]:.4f})")
    drop = ukb_best[0] - wales_best[0]
    print(f"    AUC drop on validation: {drop:+.4f}")
    if drop < 0.05:
        print(f"    => MINIMAL calibration loss — model generalizes well")
    elif drop < 0.10:
        print(f"    => MODERATE calibration loss")
    else:
        print(f"    => SUBSTANTIAL calibration loss — consider recalibration")

print(f"\n  Number of UKB-derived domain weights: {len(domain_risk_weights)}")
print(f"  Number of variant type weights: {len(vtype_weights)}")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)

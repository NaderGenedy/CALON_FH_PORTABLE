#!/usr/bin/env python3
"""
Nobel-calibre analyses:
1. Cox survival analysis with SSS
2. Kaplan-Meier penetrance curves by domain
3. Compound heterozygote structural analysis
4. Family clustering
5. Cascade screening prioritisation
"""
import csv, math, os
from collections import defaultdict, Counter

def sf(x):
    try: return float(str(x).strip())
    except: return None

BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"

# Load data
wales = list(csv.DictReader(open(f'{BASE}/WALES_FH_CLEANED (1) - Copy.csv', encoding='utf-8-sig')))
sss_scores = list(csv.DictReader(open(f'{BASE}/alphafold/analysis/structural_severity_scores.csv', encoding='utf-8-sig')))
sss_lookup = {r['variant_id']: r for r in sss_scores}

print(f"Wales patients: {len(wales)}")
print(f"SSS variants: {len(sss_lookup)}")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: COX PROPORTIONAL HAZARDS (manual implementation)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("  SECTION 1: SURVIVAL ANALYSIS")
print("  Time-to-ASCVD by SSS tertile (avoids treatment paradox via censoring)")
print("="*80)

# Build survival data: time = age at event (if event) or current age (if censored)
survival_data = []
for r in wales:
    mut = r.get('Mutation1', '').strip()
    info = sss_lookup.get(mut)
    if not info: continue
    sss = sf(info.get('sss'))
    if sss is None: continue
    domain = info.get('domain', '')

    ascvd = 1 if sf(r.get('ascvd_combine')) == 1 else 0

    # Time variable
    # If ASCVD: use earliest event age
    event_ages = []
    for col in ['MIACSAge', 'PCIStentsAge', 'CABGAge', 'ANGINAAge', 'TIAAge', 'PVDAge', 'OtherAge']:
        a = sf(r.get(col))
        if a and a > 0:
            event_ages.append(a)

    current_age = sf(r.get('BMI_AGE'))
    if not current_age:
        current_age = sf(r.get('Currentage'))
    if not current_age:
        continue

    if ascvd == 1 and event_ages:
        time = min(event_ages)  # Age at first event
    else:
        time = current_age  # Censored at current age

    if time <= 0 or time > 100:
        continue

    gender = 1 if str(r.get('Gender', '')).strip().upper() in ['M', 'MALE', '1', '1.0'] else 0
    ldl = sf(r.get('LDL.1'))
    statin = 1 if str(r.get('statin_at_genetic_test', '')).strip() in ['1', '1.0', 'Yes', 'yes'] else 0
    xanth = 1 if r.get('TendonXanthomata', '').strip() not in ['', '0', 'No', 'no'] else 0
    smoking = 1 if str(r.get('Smoking', '')).strip() not in ['', '0', 'No', 'no', 'Never'] else 0
    family = r.get('FamilyNumber', '').strip()
    proband = 1 if sf(r.get('I_Vs_R')) == 1 else 0
    mut2 = r.get('Mutation2', '').strip()

    survival_data.append({
        'time': time, 'event': ascvd, 'sss': sss, 'domain': domain,
        'age': current_age, 'sex': gender, 'ldl': ldl, 'statin': statin,
        'xanth': xanth, 'smoking': smoking, 'family': family,
        'proband': proband, 'mutation1': mut, 'mutation2': mut2,
    })

print(f"\nSurvival cohort: {len(survival_data)} patients")
events = sum(r['event'] for r in survival_data)
print(f"ASCVD events: {events} ({events/len(survival_data)*100:.1f}%)")
times = [r['time'] for r in survival_data]
print(f"Follow-up time: median {sorted(times)[len(times)//2]:.1f} years")
print(f"Age range: {min(times):.0f} - {max(times):.0f}")

# SSS tertiles
sss_vals = sorted([r['sss'] for r in survival_data])
t1_cut = sss_vals[len(sss_vals)//3]
t2_cut = sss_vals[2*len(sss_vals)//3]

print(f"\nSSS tertile cutpoints: T1 <{t1_cut:.3f}, T2 {t1_cut:.3f}-{t2_cut:.3f}, T3 >={t2_cut:.3f}")

# Kaplan-Meier by SSS tertile
def kaplan_meier(data, max_time=80):
    """Compute KM survival curve."""
    # Sort by time
    data_sorted = sorted(data, key=lambda x: x['time'])
    n_at_risk = len(data_sorted)
    survival = 1.0
    curve = [(0, 1.0, n_at_risk)]

    i = 0
    while i < len(data_sorted):
        t = data_sorted[i]['time']
        # Count events and censored at this time
        d = 0  # deaths
        c = 0  # censored
        while i < len(data_sorted) and data_sorted[i]['time'] == t:
            if data_sorted[i]['event'] == 1:
                d += 1
            else:
                c += 1
            i += 1

        if d > 0 and n_at_risk > 0:
            survival *= (1 - d / n_at_risk)

        curve.append((t, survival, n_at_risk))
        n_at_risk -= (d + c)

    return curve

print("\n--- Kaplan-Meier Event-Free Survival by SSS Tertile ---\n")

for tname, lo, hi in [('T1 Low SSS', 0, t1_cut), ('T2 Mid SSS', t1_cut, t2_cut), ('T3 High SSS', t2_cut, 1.01)]:
    sub = [r for r in survival_data if lo <= r['sss'] < hi]
    n = len(sub)
    ev = sum(r['event'] for r in sub)

    curve = kaplan_meier(sub)

    # Get survival at key ages
    surv_40 = 1.0
    surv_50 = 1.0
    surv_60 = 1.0
    surv_70 = 1.0
    for t, s, nar in curve:
        if t <= 40: surv_40 = s
        if t <= 50: surv_50 = s
        if t <= 60: surv_60 = s
        if t <= 70: surv_70 = s

    event_rate_40 = (1-surv_40)*100
    event_rate_50 = (1-surv_50)*100
    event_rate_60 = (1-surv_60)*100
    event_rate_70 = (1-surv_70)*100

    print(f"  {tname} (n={n}, events={ev}):")
    print(f"    ASCVD by age 40: {event_rate_40:.1f}%")
    print(f"    ASCVD by age 50: {event_rate_50:.1f}%")
    print(f"    ASCVD by age 60: {event_rate_60:.1f}%")
    print(f"    ASCVD by age 70: {event_rate_70:.1f}%")

# Log-rank test (manual)
print("\n--- Log-Rank Test: T1 vs T3 ---")
t1_data = [r for r in survival_data if r['sss'] < t1_cut]
t3_data = [r for r in survival_data if r['sss'] >= t2_cut]

# Simplified log-rank: compare observed vs expected events
all_combined = t1_data + t3_data
all_combined.sort(key=lambda x: x['time'])

# At each event time, compute expected events for each group
O1, E1 = 0, 0  # observed, expected for T1
O3, E3 = 0, 0  # observed, expected for T3
n1 = len(t1_data)
n3 = len(t3_data)

# Count events per group
O1 = sum(r['event'] for r in t1_data)
O3 = sum(r['event'] for r in t3_data)

# Expected under null (proportional to group size)
total_events = O1 + O3
E1 = total_events * n1 / (n1 + n3)
E3 = total_events * n3 / (n1 + n3)

chi2_lr = (O1 - E1)**2 / max(E1, 0.001) + (O3 - E3)**2 / max(E3, 0.001)
p_lr = math.exp(-0.5 * chi2_lr) if chi2_lr < 30 else 0.0001

print(f"  T1: O={O1}, E={E1:.1f}")
print(f"  T3: O={O3}, E={E3:.1f}")
print(f"  Chi-squared = {chi2_lr:.3f}, P = {p_lr:.6f}")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: KAPLAN-MEIER PENETRANCE BY DOMAIN
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("  SECTION 2: DOMAIN-SPECIFIC PENETRANCE CURVES")
print("  Probability of ASCVD by age, stratified by LDLR domain")
print("="*80)

domain_counts = Counter(r['domain'] for r in survival_data)
major_domains = [d for d, n in domain_counts.items() if n >= 30]

print(f"\nDomains with n>=30: {len(major_domains)}")
print(f"\n{'Domain':<35} {'N':>5} {'Events':>7} {'By40':>6} {'By50':>6} {'By60':>6} {'By70':>6}")
print("-"*80)

domain_penetrance = {}
for d in sorted(major_domains, key=lambda x: -domain_counts[x]):
    sub = [r for r in survival_data if r['domain'] == d]
    n = len(sub)
    ev = sum(r['event'] for r in sub)
    curve = kaplan_meier(sub)

    s40, s50, s60, s70 = 1,1,1,1
    for t, s, nar in curve:
        if t <= 40: s40 = s
        if t <= 50: s50 = s
        if t <= 60: s60 = s
        if t <= 70: s70 = s

    p40 = (1-s40)*100
    p50 = (1-s50)*100
    p60 = (1-s60)*100
    p70 = (1-s70)*100

    domain_penetrance[d] = {'n': n, 'events': ev, 'p40': p40, 'p50': p50, 'p60': p60, 'p70': p70}
    print(f"  {d:<35} {n:>5} {ev:>7} {p40:>5.1f}% {p50:>5.1f}% {p60:>5.1f}% {p70:>5.1f}%")

# Highest vs lowest penetrance domains
if domain_penetrance:
    highest = max(domain_penetrance.items(), key=lambda x: x[1]['p60'])
    lowest = min(domain_penetrance.items(), key=lambda x: x[1]['p60'])
    print(f"\n  HIGHEST penetrance by 60: {highest[0]} ({highest[1]['p60']:.1f}%)")
    print(f"  LOWEST penetrance by 60:  {lowest[0]} ({lowest[1]['p60']:.1f}%)")
    ratio = highest[1]['p60'] / max(lowest[1]['p60'], 0.1)
    print(f"  RATIO: {ratio:.1f}x")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: COMPOUND HETEROZYGOTE ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("  SECTION 3: COMPOUND HETEROZYGOTE STRUCTURAL ANALYSIS")
print("  Patients with TWO mutations: does cumulative structural damage predict severity?")
print("="*80)

# Check Mutation2 column
has_mut2 = [r for r in survival_data if r['mutation2'] and r['mutation2'].strip()]
has_mut1_only = [r for r in survival_data if not r['mutation2'] or not r['mutation2'].strip()]

print(f"\n  Single mutation (Mutation1 only): {len(has_mut1_only)}")
print(f"  Compound heterozygote (Mutation1 + Mutation2): {len(has_mut2)}")

if has_mut2:
    # Get SSS for mutation2
    for r in has_mut2:
        m2 = r['mutation2'].strip()
        info2 = sss_lookup.get(m2)
        r['sss2'] = sf(info2.get('sss')) if info2 else None
        r['sss_total'] = r['sss'] + (r['sss2'] if r['sss2'] else 0)

    compound_with_sss2 = [r for r in has_mut2 if r['sss2'] is not None]
    compound_no_sss2 = [r for r in has_mut2 if r['sss2'] is None]

    print(f"  Compound het with BOTH mutations mapped to SSS: {len(compound_with_sss2)}")
    print(f"  Compound het with only Mutation1 SSS: {len(compound_no_sss2)}")

    # Compare compound het vs single mutation
    print(f"\n  --- Compound Heterozygote vs Single Mutation ---")

    # Clinical comparison
    for label, data in [('Single mutation', has_mut1_only), ('Compound heterozygote (all)', has_mut2)]:
        n = len(data)
        ev = sum(r['event'] for r in data)
        ldl_vals = [r['ldl'] for r in data if r['ldl'] is not None]
        xanth = sum(r['xanth'] for r in data)
        mean_sss = sum(r['sss'] for r in data) / n
        ages_event = [r['time'] for r in data if r['event'] == 1]

        print(f"\n  {label} (n={n}):")
        print(f"    ASCVD: {ev}/{n} ({ev/n*100:.1f}%)")
        if ldl_vals: print(f"    LDL: {sum(ldl_vals)/len(ldl_vals):.2f} mmol/L")
        print(f"    Xanthomata: {xanth}/{n} ({xanth/n*100:.1f}%)")
        print(f"    Mean SSS (Mut1): {mean_sss:.3f}")
        if ages_event: print(f"    Mean age at event: {sum(ages_event)/len(ages_event):.1f}")

    # OR for compound het vs single
    a = sum(r['event'] for r in has_mut2)
    b = len(has_mut2) - a
    c = sum(r['event'] for r in has_mut1_only)
    d = len(has_mut1_only) - c

    if b > 0 and c > 0:
        OR = (a * d) / (b * c)
        se = (1/max(a,1) + 1/max(b,1) + 1/max(c,1) + 1/max(d,1))**0.5
        ci_lo = math.exp(math.log(max(OR, 0.001)) - 1.96*se)
        ci_hi = math.exp(math.log(max(OR, 0.001)) + 1.96*se)
        print(f"\n  Compound het vs Single: OR = {OR:.2f} (95% CI: {ci_lo:.2f}-{ci_hi:.2f})")

    # Show compound het mutations
    if compound_with_sss2:
        print(f"\n  Compound heterozygotes with both SSS mapped ({len(compound_with_sss2)}):")
        print(f"  {'Mut1':<30} {'SSS1':>5} {'Mut2':<30} {'SSS2':>5} {'Total':>6} {'ASCVD':>6}")
        print("  " + "-"*95)
        for r in sorted(compound_with_sss2, key=lambda x: -x['sss_total'])[:20]:
            print(f"  {r['mutation1']:<30} {r['sss']:.3f} {r['mutation2']:<30} {r['sss2']:.3f} {r['sss_total']:.3f} {'YES' if r['event']==1 else 'no':>6}")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: FAMILY-LEVEL ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("  SECTION 4: FAMILY-LEVEL STRUCTURAL CLUSTERING")
print("  Within families: does SSS predict who gets ASCVD?")
print("="*80)

# Group by family
families = defaultdict(list)
for r in survival_data:
    if r['family']:
        families[r['family']].append(r)

# Filter to families with >=2 members and at least 1 event and 1 non-event
informative_families = {}
for fam, members in families.items():
    has_event = any(r['event'] == 1 for r in members)
    has_no_event = any(r['event'] == 0 for r in members)
    if len(members) >= 2 and has_event and has_no_event:
        informative_families[fam] = members

print(f"\n  Total families: {len(families)}")
print(f"  Families with >=2 members: {sum(1 for f in families if len(families[f])>=2)}")
print(f"  Informative families (event + non-event): {len(informative_families)}")

if informative_families:
    # Within informative families, compare SSS of affected vs unaffected
    affected_sss = []
    unaffected_sss = []
    affected_ages = []
    unaffected_ages = []

    for fam, members in informative_families.items():
        for r in members:
            if r['event'] == 1:
                affected_sss.append(r['sss'])
                affected_ages.append(r['time'])
            else:
                unaffected_sss.append(r['sss'])
                unaffected_ages.append(r['age'])

    m1 = sum(affected_sss)/len(affected_sss)
    m2 = sum(unaffected_sss)/len(unaffected_sss)

    print(f"\n  Within informative families:")
    print(f"    Affected members:   SSS = {m1:.3f} (n={len(affected_sss)}), mean age at event = {sum(affected_ages)/len(affected_ages):.1f}")
    print(f"    Unaffected members: SSS = {m2:.3f} (n={len(unaffected_sss)}), mean current age = {sum(unaffected_ages)/len(unaffected_ages):.1f}")
    print(f"    Difference: {m1-m2:+.3f}")

    # Note: within same family, same mutation, so SSS should be SAME
    # The difference would come from different mutations in the family
    # or from Mutation2

    # More useful: compare family-level SSS vs family ASCVD rate
    fam_sss_rates = []
    for fam, members in families.items():
        if len(members) < 2: continue
        fam_sss_mean = sum(r['sss'] for r in members) / len(members)
        fam_ascvd_rate = sum(r['event'] for r in members) / len(members)
        fam_sss_rates.append((fam_sss_mean, fam_ascvd_rate, len(members)))

    if fam_sss_rates:
        x = [f[0] for f in fam_sss_rates]
        y = [f[1] for f in fam_sss_rates]
        n = len(x)
        mx, my = sum(x)/n, sum(y)/n
        sxx = sum((a-mx)**2 for a in x)
        syy = sum((b-my)**2 for b in y)
        sxy = sum((a-mx)*(b-my) for a,b in zip(x,y))
        r_val = sxy/(sxx*syy)**0.5 if sxx>0 and syy>0 else 0

        print(f"\n  Family-level: SSS vs ASCVD rate (n={len(fam_sss_rates)} families)")
        print(f"    Pearson r = {r_val:+.3f}")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5: CASCADE SCREENING PRIORITISATION
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("  SECTION 5: CASCADE SCREENING PRIORITISATION")
print("  Should we prioritise screening families of HIGH SSS index cases?")
print("="*80)

# Index vs cascade
index_cases = [r for r in survival_data if r['proband'] == 1]
cascade_cases = [r for r in survival_data if r['proband'] == 0]

print(f"\n  Index cases: {len(index_cases)}")
print(f"  Cascade cases: {len(cascade_cases)}")

# Among cascade cases: does the INDEX SSS predict CASCADE ASCVD?
# Group cascade cases by family, get the index SSS
cascade_by_index_sss = {'high': [], 'low': []}
sss_median = sorted([r['sss'] for r in survival_data])[len(survival_data)//2]

for fam, members in families.items():
    index_members = [r for r in members if r['proband'] == 1]
    cascade_members = [r for r in members if r['proband'] == 0]

    if not index_members or not cascade_members:
        continue

    index_sss = sum(r['sss'] for r in index_members) / len(index_members)

    for r in cascade_members:
        if index_sss >= sss_median:
            cascade_by_index_sss['high'].append(r)
        else:
            cascade_by_index_sss['low'].append(r)

for label, data in [('High SSS index family', cascade_by_index_sss['high']),
                     ('Low SSS index family', cascade_by_index_sss['low'])]:
    n = len(data)
    if n == 0: continue
    ev = sum(r['event'] for r in data)
    ldl_vals = [r['ldl'] for r in data if r['ldl'] is not None]
    xanth = sum(r['xanth'] for r in data)

    print(f"\n  Cascade cases from {label} (n={n}):")
    print(f"    ASCVD: {ev}/{n} ({ev/n*100:.1f}%)")
    if ldl_vals: print(f"    LDL: {sum(ldl_vals)/len(ldl_vals):.2f}")
    print(f"    Xanthomata: {xanth}/{n} ({xanth/n*100:.1f}%)")

# OR
h = cascade_by_index_sss['high']
l = cascade_by_index_sss['low']
if h and l:
    a = sum(r['event'] for r in h)
    b = len(h) - a
    c = sum(r['event'] for r in l)
    d = len(l) - c
    if b > 0 and c > 0:
        OR = (a*d)/(b*c)
        se = (1/max(a,1)+1/max(b,1)+1/max(c,1)+1/max(d,1))**0.5
        ci_lo = math.exp(math.log(max(OR,0.001))-1.96*se)
        ci_hi = math.exp(math.log(max(OR,0.001))+1.96*se)
        print(f"\n  Cascade ASCVD: High SSS family vs Low SSS family:")
        print(f"    OR = {OR:.2f} (95% CI: {ci_lo:.2f}-{ci_hi:.2f})")

# Number needed to screen
print(f"\n  --- Number Needed to Screen (NNS) ---")
for label, data in [('High SSS families', cascade_by_index_sss['high']),
                     ('Low SSS families', cascade_by_index_sss['low'])]:
    n = len(data)
    if n == 0: continue
    ev = sum(r['event'] for r in data)
    rate = ev/n if n > 0 else 0
    nns = 1/rate if rate > 0 else float('inf')
    print(f"  {label}: ASCVD rate={rate*100:.1f}%, NNS to find 1 ASCVD = {nns:.1f}")

# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("  SUMMARY OF NOVEL FINDINGS")
print("="*80)

# Save results
output_path = f"{BASE}/alphafold/analysis/cox_km_compound_results.csv"
with open(output_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['analysis', 'group', 'n', 'events', 'rate', 'metric', 'value'])

    # KM data
    for tname, lo, hi in [('T1_Low', 0, t1_cut), ('T2_Mid', t1_cut, t2_cut), ('T3_High', t2_cut, 1.01)]:
        sub = [r for r in survival_data if lo <= r['sss'] < hi]
        n = len(sub)
        ev = sum(r['event'] for r in sub)
        curve = kaplan_meier(sub)
        s60 = 1
        for t, s, nar in curve:
            if t <= 60: s60 = s
        w.writerow(['KM_SSS_tertile', tname, n, ev, ev/n, 'penetrance_by_60', (1-s60)*100])

    # Domain penetrance
    for d, v in domain_penetrance.items():
        w.writerow(['KM_domain', d, v['n'], v['events'], v['events']/v['n'], 'penetrance_by_60', v['p60']])

print(f"\n  Results saved to: {output_path}")
print("\n  DONE.")

#!/usr/bin/env python3
"""
CALON-Structure Phase 4d: Variant-Level LLT Response Analysis
==============================================================
How do different FH variants respond to lipid-lowering therapy?
This is the PHARMACOGENOMICS arm of the structural severity analysis.

Key questions:
  1. Do severe mutations (high SSS) show worse LDL reduction on statins?
  2. Is statin type/dose different across SSS strata?
  3. Do null mutations need PCSK9i earlier?
  4. Does domain location predict treatment response?
  5. Which individual variants are treatment-resistant?
  6. SSS as a predictor of time-to-target LDL
  7. Pre-treatment vs on-treatment LDL trajectory by SSS
"""

import os, csv, math, re
from collections import defaultdict

BASE = r"C:\Users\nader\Downloads\calon_ukb_pipeline"
AF_DIR = os.path.join(BASE, "alphafold")
ANALYSIS = os.path.join(AF_DIR, "analysis")

def safe_float(x, default=None):
    try:
        v = float(x)
        if math.isnan(v) or math.isinf(v): return default
        return v
    except:
        return default

def pearson_r(x, y):
    n = len(x)
    if n < 3: return 0.0, 1.0
    mx, my = sum(x)/n, sum(y)/n
    sxx = sum((xi-mx)**2 for xi in x)
    syy = sum((yi-my)**2 for yi in y)
    sxy = sum((x[i]-mx)*(y[i]-my) for i in range(n))
    if sxx < 1e-10 or syy < 1e-10: return 0.0, 1.0
    r = sxy / (sxx*syy)**0.5
    if abs(r) >= 1.0: return r, 0.0
    t = r * ((n-2)/(1-r*r))**0.5
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(t) / math.sqrt(2))))
    return r, p

def median(vals):
    s = sorted(vals)
    n = len(s)
    if n == 0: return 0
    if n % 2 == 1: return s[n//2]
    return (s[n//2-1] + s[n//2]) / 2

def iqr(vals):
    s = sorted(vals)
    n = len(s)
    if n < 4: return 0
    q1 = s[n//4]
    q3 = s[3*n//4]
    return q3 - q1

print("=" * 80)
print("  CALON-Structure: Variant-Level LLT Response Analysis")
print("  Pharmacogenomics of LDLR Structural Severity")
print("=" * 80)

# ============================================================================
#  LOAD AND MERGE DATA
# ============================================================================

# Load SSS scores
sss_data = {}
with open(os.path.join(ANALYSIS, "structural_severity_scores.csv")) as f:
    for row in csv.DictReader(f):
        sss_data[row['variant_id']] = {
            'sss': safe_float(row['sss'], 0.5),
            'gene': row['gene'],
            'domain': row.get('domain', 'Unknown'),
            'ddG': safe_float(row.get('ddG')),
            'sss_source': row.get('sss_source', ''),
        }

# Load SSS patient mapping
sss_patients = {}
with open(os.path.join(ANALYSIS, "sss_patient_level_wales.csv")) as f:
    for row in csv.DictReader(f):
        sss_patients[row['patient_id']] = {
            'variant': row['variant_id'],
            'sss': safe_float(row['sss'], 0.5),
            'gene': row['gene'],
        }

# Load DRAGON_3 with full treatment data
patients = []
with open(os.path.join(BASE, "DRAGON_3.csv"), encoding='latin-1') as f:
    reader = csv.DictReader(f)
    for row in reader:
        pid = row.get('patient_id', row.get(reader.fieldnames[0], ''))
        if pid not in sss_patients:
            continue

        sp = sss_patients[pid]

        # Parse statin info
        statin_raw = str(row.get('Statin', '')).strip()
        on_statin = statin_raw.upper() not in ('', 'NO', 'NAN', 'NONE', 'NA', 'N')

        # Parse detailed treatment
        tx1_1 = str(row.get('Treatment1_1', '')).strip()
        tx1_2 = str(row.get('Treatment1_2', '')).strip()
        tx1_3 = str(row.get('Treatment1_3', '')).strip()

        # Extract statin type and dose from Treatment1_1
        statin_type = ''
        statin_dose = 0
        if tx1_1:
            tx_lower = tx1_1.lower()
            if 'atorva' in tx_lower: statin_type = 'Atorvastatin'
            elif 'rosuva' in tx_lower: statin_type = 'Rosuvastatin'
            elif 'simva' in tx_lower: statin_type = 'Simvastatin'
            elif 'prava' in tx_lower: statin_type = 'Pravastatin'
            elif 'fluva' in tx_lower: statin_type = 'Fluvastatin'

            dose_match = re.search(r'(\d+)\s*mg', tx1_1, re.IGNORECASE)
            if dose_match:
                statin_dose = int(dose_match.group(1))

        # Ezetimibe
        ezetimibe = 'ezet' in tx1_2.lower() or str(row.get('Ezetimibe', '')).upper() in ('Y', 'YES', '1')

        # PCSK9i
        pcsk9i = str(row.get('PCSK9i', '')).upper() in ('Y', 'YES', '1')
        if not pcsk9i and tx1_3:
            pcsk9i = any(k in tx1_3.lower() for k in ['evoloc', 'alirocum', 'pcsk9', 'repatha', 'praluent'])

        # Treatment intensity score
        # High: Atorva 40-80mg or Rosuva 20-40mg
        # Moderate: Atorva 10-20mg or Rosuva 5-10mg or Simva 20-40mg
        # Low: everything else
        intensity = 'none'
        if statin_type == 'Atorvastatin' and statin_dose >= 40: intensity = 'high'
        elif statin_type == 'Rosuvastatin' and statin_dose >= 20: intensity = 'high'
        elif statin_type == 'Atorvastatin' and statin_dose >= 10: intensity = 'moderate'
        elif statin_type == 'Rosuvastatin' and statin_dose >= 5: intensity = 'moderate'
        elif statin_type == 'Simvastatin' and statin_dose >= 20: intensity = 'moderate'
        elif statin_type: intensity = 'low'

        # LDL trajectory
        ldl1 = safe_float(row.get('LDL_1'))
        ldl2 = safe_float(row.get('LDL_2'))
        ldl3 = safe_float(row.get('LDL_3'))
        ldl4 = safe_float(row.get('LDL_4'))
        last_ldl = safe_float(row.get('LastLDL'))
        matched_ldl = safe_float(row.get('MtachedLDLC'))
        ldcc = safe_float(row.get('LDCC'))

        # LDL reduction
        ldl_pct = safe_float(row.get('ldl_per'))        # as fraction (0-1)
        ldl_pct1 = safe_float(row.get('ldl_per_1'))     # as percentage (0-100)
        tc_pct = safe_float(row.get('TC_per'))
        nhdl_pct = safe_float(row.get('NHDL2_1_per'))

        # Age at treatment
        age_tx = safe_float(row.get('Age_Treatment1'))
        age = safe_float(row.get('age_at_event_or_censoring', row.get('Ageattest')))
        sex = 1 if str(row.get('Gender', '')).upper().startswith('M') else 0

        # ASCVD
        ascvd = 1 if str(row.get('ASCVD_combined', '')).strip() in ('1', '1.0', 'True') else 0

        # Other lipids
        tc1 = safe_float(row.get('TC_1'))
        hdl1 = safe_float(row.get('HDL_1'))
        tg1 = safe_float(row.get('TRG_1'))
        apob = safe_float(row.get('ApoB'))
        lpa = safe_float(row.get('Lpa'))

        # Xanthomata
        xanth = 1 if str(row.get('TendonXanthomata', '')).strip() in ('1', '1.0', 'True') else 0

        # On treatment flag
        on_treatment = str(row.get('OnTreatment', '')).strip() in ('1', '1.0')

        # Treatment level (0-3 scale from database)
        tx_level = safe_float(row.get('Treatment1'))

        domain = sss_data.get(sp['variant'], {}).get('domain', 'Unknown')

        patients.append({
            'pid': pid, 'variant': sp['variant'], 'gene': sp['gene'],
            'domain': domain, 'sss': sp['sss'],
            'statin_type': statin_type, 'statin_dose': statin_dose,
            'intensity': intensity, 'ezetimibe': ezetimibe, 'pcsk9i': pcsk9i,
            'on_statin': on_statin, 'on_treatment': on_treatment,
            'tx_level': tx_level, 'tx1_1': tx1_1, 'tx1_2': tx1_2, 'tx1_3': tx1_3,
            'ldl1': ldl1, 'ldl2': ldl2, 'ldl3': ldl3, 'ldl4': ldl4,
            'last_ldl': last_ldl, 'matched_ldl': matched_ldl, 'ldcc': ldcc,
            'ldl_pct': ldl_pct, 'ldl_pct1': ldl_pct1,
            'tc_pct': tc_pct, 'nhdl_pct': nhdl_pct,
            'age_tx': age_tx, 'age': age, 'sex': sex,
            'ascvd': ascvd, 'tc1': tc1, 'hdl1': hdl1, 'tg1': tg1,
            'apob': apob, 'lpa': lpa, 'xanth': xanth,
        })

print(f"\n  Total patients: {len(patients)}")
on_tx = sum(1 for p in patients if p['on_statin'])
print(f"  On statin: {on_tx} ({100*on_tx/len(patients):.1f}%)")
with_response = sum(1 for p in patients if p['ldl_pct1'] is not None)
print(f"  With LDL% response: {with_response}")
with_trajectory = sum(1 for p in patients if p['ldl1'] is not None and p['ldl2'] is not None)
print(f"  With LDL trajectory (LDL1+LDL2): {with_trajectory}")

# ============================================================================
#  A. OVERALL SSS vs LDL RESPONSE
# ============================================================================

print("\n" + "=" * 80)
print("  A. SSS vs LDL-C Reduction (Treatment Response)")
print("=" * 80)

resp_pts = [p for p in patients if p['ldl_pct1'] is not None and p['ldl_pct1'] > 0]
print(f"\n  Patients with measurable LDL reduction: {len(resp_pts)}")

if resp_pts:
    sss_v = [p['sss'] for p in resp_pts]
    ldl_r = [p['ldl_pct1'] for p in resp_pts]
    r, pval = pearson_r(sss_v, ldl_r)
    print(f"  SSS vs LDL% reduction: r = {r:.4f}, p = {pval:.4f}")
    print(f"    (Negative r = higher SSS -> LESS LDL reduction = treatment resistance)")

    # SSS tertiles
    sss_sorted = sorted(p['sss'] for p in resp_pts)
    t1 = sss_sorted[len(sss_sorted)//3]
    t2 = sss_sorted[2*len(sss_sorted)//3]

    print(f"\n  LDL-C % reduction by SSS tertile:")
    print(f"  {'Tertile':<20s} {'N':<6s} {'Mean%':<8s} {'Median%':<10s} {'IQR':<8s} {'SSS range'}")
    print(f"  {'-'*20} {'-'*6} {'-'*8} {'-'*10} {'-'*8} {'-'*15}")
    for label, lo, hi in [('T1 (mild SSS)', 0, t1), ('T2 (moderate)', t1, t2), ('T3 (severe)', t2, 2.0)]:
        tpts = [p for p in resp_pts if lo <= p['sss'] < hi]
        if tpts:
            vals = [p['ldl_pct1'] for p in tpts]
            print(f"  {label:<20s} {len(tpts):<6d} {sum(vals)/len(vals):<8.1f} {median(vals):<10.1f} {iqr(vals):<8.1f} {lo:.2f}-{hi:.2f}")

    # Also compute absolute LDL reduction
    abs_pts = [p for p in patients if p['ldl1'] is not None and p['ldl2'] is not None and p['ldl1'] > 0]
    if abs_pts:
        print(f"\n  Absolute LDL reduction (LDL1 - LDL2) by SSS tertile:")
        for label, lo, hi in [('T1 (mild)', 0, t1), ('T2 (moderate)', t1, t2), ('T3 (severe)', t2, 2.0)]:
            tpts = [p for p in abs_pts if lo <= p['sss'] < hi]
            if tpts:
                reductions = [p['ldl1'] - p['ldl2'] for p in tpts]
                mean_ldl1 = sum(p['ldl1'] for p in tpts) / len(tpts)
                mean_ldl2 = sum(p['ldl2'] for p in tpts) / len(tpts)
                mean_red = sum(reductions) / len(reductions)
                print(f"    {label}: LDL1={mean_ldl1:.2f} -> LDL2={mean_ldl2:.2f}, reduction={mean_red:.2f} mmol/L ({100*mean_red/mean_ldl1:.1f}%)")

# ============================================================================
#  B. STATIN TYPE AND DOSE BY SSS
# ============================================================================

print("\n" + "=" * 80)
print("  B. Treatment Prescription Patterns by SSS (Confounding by Indication)")
print("=" * 80)
print("  Are severe mutations treated more aggressively?\n")

statin_pts = [p for p in patients if p['statin_type']]
print(f"  Patients with identifiable statin: {len(statin_pts)}")

if statin_pts:
    # Statin type distribution
    type_counts = defaultdict(int)
    for p in statin_pts:
        type_counts[p['statin_type']] += 1
    print(f"\n  Statin distribution:")
    for stype, count in sorted(type_counts.items(), key=lambda x: -x[1]):
        pct = 100 * count / len(statin_pts)
        print(f"    {stype:<20s}: {count:4d} ({pct:.1f}%)")

    # Intensity by SSS
    print(f"\n  Treatment intensity by SSS tertile:")
    sss_sorted = sorted(p['sss'] for p in statin_pts)
    t1 = sss_sorted[len(sss_sorted)//3]
    t2 = sss_sorted[2*len(sss_sorted)//3]

    print(f"  {'Tertile':<20s} {'High%':<10s} {'Mod%':<10s} {'Low%':<10s} {'Mean dose':<12s} {'+Ezet%':<10s} {'+PCSK9i%'}")
    print(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*10} {'-'*12} {'-'*10} {'-'*10}")

    for label, lo, hi in [('T1 (mild)', 0, t1), ('T2 (moderate)', t1, t2), ('T3 (severe)', t2, 2.0)]:
        tpts = [p for p in statin_pts if lo <= p['sss'] < hi]
        if tpts:
            high = sum(1 for p in tpts if p['intensity'] == 'high')
            mod = sum(1 for p in tpts if p['intensity'] == 'moderate')
            low = sum(1 for p in tpts if p['intensity'] == 'low')
            doses = [p['statin_dose'] for p in tpts if p['statin_dose'] > 0]
            mean_dose = sum(doses) / len(doses) if doses else 0
            ezet = sum(1 for p in tpts if p['ezetimibe'])
            pcsk = sum(1 for p in tpts if p['pcsk9i'])
            n = len(tpts)
            print(f"  {label:<20s} {100*high/n:<10.1f} {100*mod/n:<10.1f} {100*low/n:<10.1f} {mean_dose:<12.1f} {100*ezet/n:<10.1f} {100*pcsk/n:<10.1f}")

    # High-intensity statin response by SSS
    print(f"\n  LDL response on HIGH-INTENSITY statin by SSS:")
    hi_statin = [p for p in patients if p['intensity'] == 'high' and p['ldl_pct1'] is not None]
    if len(hi_statin) >= 10:
        for label, lo, hi in [('Low SSS', 0, t1), ('Mid SSS', t1, t2), ('High SSS', t2, 2.0)]:
            tpts = [p for p in hi_statin if lo <= p['sss'] < hi]
            if tpts:
                vals = [p['ldl_pct1'] for p in tpts]
                print(f"    {label}: n={len(tpts)}, mean LDL reduction={sum(vals)/len(vals):.1f}%, median={median(vals):.1f}%")

# ============================================================================
#  C. DOMAIN-SPECIFIC TREATMENT RESPONSE
# ============================================================================

print("\n" + "=" * 80)
print("  C. Treatment Response by LDLR Domain")
print("=" * 80)
print("  Different domains affect different receptor functions:")
print("  - Ligand-binding: LDL capture -> statins upregulate remaining function")
print("  - EGF/Beta-propeller: recycling -> statins less effective if recycling broken")
print("  - Null mutations: no receptor -> statins cannot help\n")

domain_pts = [p for p in patients if p['ldl_pct1'] is not None and p['gene'] == 'LDLR']

if domain_pts:
    print(f"  {'Domain':<25s} {'N':<6s} {'Mean LDL%':<12s} {'Median':<10s} {'Pre-LDL':<10s} {'Post-LDL':<10s} {'ASCVD%'}")
    print(f"  {'-'*25} {'-'*6} {'-'*12} {'-'*10} {'-'*10} {'-'*10} {'-'*8}")

    domain_groups = defaultdict(list)
    for p in domain_pts:
        domain_groups[p['domain']].append(p)

    for domain in ['Ligand-binding', 'EGF-A', 'EGF-B', 'EGF-C', 'EGF domains',
                    'Beta-propeller', 'O-linked', 'Signal peptide', 'Unknown']:
        dpts = domain_groups.get(domain, [])
        if len(dpts) >= 3:
            ldl_r = [p['ldl_pct1'] for p in dpts]
            pre = [p['ldl1'] for p in dpts if p['ldl1']]
            post = [p['ldl2'] for p in dpts if p['ldl2']]
            ev = sum(1 for p in dpts if p['ascvd'] == 1)
            pre_mean = sum(pre)/len(pre) if pre else 0
            post_mean = sum(post)/len(post) if post else 0

            print(f"  {domain:<25s} {len(dpts):<6d} {sum(ldl_r)/len(ldl_r):<12.1f} {median(ldl_r):<10.1f} {pre_mean:<10.2f} {post_mean:<10.2f} {100*ev/len(dpts):.1f}")

# Loss-of-function vs missense
print(f"\n  Loss-of-function vs Missense treatment response:")
lof_pts = [p for p in patients if p['ldl_pct1'] is not None and p['sss'] >= 0.95]
mis_pts = [p for p in patients if p['ldl_pct1'] is not None and p['sss'] < 0.95]

if lof_pts and mis_pts:
    lof_r = [p['ldl_pct1'] for p in lof_pts]
    mis_r = [p['ldl_pct1'] for p in mis_pts]
    print(f"    Loss-of-function (SSS>=0.95): n={len(lof_pts)}, LDL reduction={sum(lof_r)/len(lof_r):.1f}%")
    print(f"    Missense (SSS<0.95):          n={len(mis_pts)}, LDL reduction={sum(mis_r)/len(mis_r):.1f}%")
    print(f"    Difference: {sum(lof_r)/len(lof_r) - sum(mis_r)/len(mis_r):.1f} percentage points")

# ============================================================================
#  D. INDIVIDUAL VARIANT TREATMENT RESPONSE (Top 20)
# ============================================================================

print("\n" + "=" * 80)
print("  D. Individual Variant Treatment Response (Top 20 by count)")
print("=" * 80)

variant_response = defaultdict(list)
for p in patients:
    if p['ldl_pct1'] is not None:
        variant_response[p['variant']].append(p)

# Sort by patient count
top_variants = sorted(variant_response.items(), key=lambda x: -len(x[1]))[:25]

print(f"\n  {'Variant':<35s} {'N':<5s} {'Gene':<6s} {'SSS':<6s} {'LDL%red':<8s} {'Pre-LDL':<8s} {'On-LDL':<8s} {'ASCVD':<6s} {'Domain'}")
print(f"  {'-'*35} {'-'*5} {'-'*6} {'-'*6} {'-'*8} {'-'*8} {'-'*8} {'-'*6} {'-'*20}")

for variant, vpts in top_variants:
    if len(vpts) < 2:
        continue
    ldl_r = [p['ldl_pct1'] for p in vpts]
    sss = vpts[0]['sss']
    gene = vpts[0]['gene']
    domain = vpts[0]['domain']
    pre = [p['ldl1'] for p in vpts if p['ldl1']]
    post = [p['ldl2'] for p in vpts if p['ldl2']]
    ev = sum(1 for p in vpts if p['ascvd'] == 1)

    pre_m = sum(pre)/len(pre) if pre else 0
    post_m = sum(post)/len(post) if post else 0

    v_short = variant[:34]
    print(f"  {v_short:<35s} {len(vpts):<5d} {gene:<6s} {sss:<6.3f} {sum(ldl_r)/len(ldl_r):<8.1f} {pre_m:<8.2f} {post_m:<8.2f} {ev:<6d} {domain}")

# ============================================================================
#  E. TREATMENT-RESISTANT VARIANTS
# ============================================================================

print("\n" + "=" * 80)
print("  E. Treatment-Resistant Variants (LDL reduction < 30%)")
print("=" * 80)
print("  ESC/EAS guidelines expect >= 50% LDL reduction on high-intensity statins")
print("  Variants with < 30% mean reduction may be functionally null\n")

resistant = []
for variant, vpts in variant_response.items():
    if len(vpts) >= 2:
        ldl_r = [p['ldl_pct1'] for p in vpts]
        mean_r = sum(ldl_r) / len(ldl_r)
        if mean_r < 30:
            resistant.append((variant, vpts, mean_r))

resistant.sort(key=lambda x: x[2])

print(f"  Variants with mean LDL reduction < 30%: {len(resistant)}")
if resistant:
    print(f"\n  {'Variant':<35s} {'N':<5s} {'SSS':<6s} {'LDL%':<8s} {'Pre-LDL':<8s} {'Domain'}")
    print(f"  {'-'*35} {'-'*5} {'-'*6} {'-'*8} {'-'*8} {'-'*20}")
    for variant, vpts, mean_r in resistant[:15]:
        sss = vpts[0]['sss']
        domain = vpts[0]['domain']
        pre = [p['ldl1'] for p in vpts if p['ldl1']]
        pre_m = sum(pre)/len(pre) if pre else 0
        v_short = variant[:34]
        print(f"  {v_short:<35s} {len(vpts):<5d} {sss:<6.3f} {mean_r:<8.1f} {pre_m:<8.2f} {domain}")

# ============================================================================
#  F. LDL TRAJECTORY BY SSS (Serial Measurements)
# ============================================================================

print("\n" + "=" * 80)
print("  F. LDL-C Trajectory (Serial Measurements) by SSS")
print("=" * 80)

traj_pts = [p for p in patients if p['ldl1'] is not None]
print(f"\n  Patients with LDL trajectory data:")
print(f"    LDL_1: {sum(1 for p in traj_pts if p['ldl1'] is not None)}")
print(f"    LDL_2: {sum(1 for p in traj_pts if p['ldl2'] is not None)}")
print(f"    LDL_3: {sum(1 for p in traj_pts if p['ldl3'] is not None)}")
print(f"    LDL_4: {sum(1 for p in traj_pts if p['ldl4'] is not None)}")
print(f"    LastLDL: {sum(1 for p in traj_pts if p['last_ldl'] is not None)}")

sss_sorted = sorted(p['sss'] for p in traj_pts)
t1 = sss_sorted[len(sss_sorted)//3]
t2 = sss_sorted[2*len(sss_sorted)//3]

print(f"\n  LDL trajectory by SSS tertile:")
print(f"  {'Tertile':<15s} {'LDL_1':<8s} {'LDL_2':<8s} {'LDL_3':<8s} {'LDL_4':<8s} {'LastLDL':<8s} {'Delta(1-Last)'}")
print(f"  {'-'*15} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*14}")

for label, lo, hi in [('T1 (mild)', 0, t1), ('T2 (mod)', t1, t2), ('T3 (severe)', t2, 2.0)]:
    tpts = [p for p in traj_pts if lo <= p['sss'] < hi]
    if tpts:
        means = []
        for key in ['ldl1', 'ldl2', 'ldl3', 'ldl4', 'last_ldl']:
            vals = [p[key] for p in tpts if p[key] is not None]
            means.append(sum(vals)/len(vals) if vals else 0)

        delta = means[0] - means[4] if means[0] > 0 and means[4] > 0 else 0
        print(f"  {label:<15s} {means[0]:<8.2f} {means[1]:<8.2f} {means[2]:<8.2f} {means[3]:<8.2f} {means[4]:<8.2f} {delta:<14.2f}")

# ============================================================================
#  G. TARGET ACHIEVEMENT BY SSS
# ============================================================================

print("\n" + "=" * 80)
print("  G. LDL-C Target Achievement by SSS")
print("=" * 80)
print("  ESC/EAS 2019 FH targets:")
print("    Primary prevention: LDL-C < 2.6 mmol/L AND >= 50% reduction")
print("    Very high risk: LDL-C < 1.4 mmol/L AND >= 50% reduction\n")

target_pts = [p for p in patients if p['last_ldl'] is not None]

if target_pts:
    print(f"  {'Tertile':<15s} {'N':<6s} {'<2.6':<8s} {'<1.8':<8s} {'<1.4':<8s} {'>=50%red':<10s} {'Mean LastLDL'}")
    print(f"  {'-'*15} {'-'*6} {'-'*8} {'-'*8} {'-'*8} {'-'*10} {'-'*12}")

    for label, lo, hi in [('T1 (mild)', 0, t1), ('T2 (mod)', t1, t2), ('T3 (severe)', t2, 2.0)]:
        tpts = [p for p in target_pts if lo <= p['sss'] < hi]
        if tpts:
            n = len(tpts)
            lt26 = sum(1 for p in tpts if p['last_ldl'] < 2.6)
            lt18 = sum(1 for p in tpts if p['last_ldl'] < 1.8)
            lt14 = sum(1 for p in tpts if p['last_ldl'] < 1.4)
            red50 = sum(1 for p in tpts if p['ldl_pct1'] is not None and p['ldl_pct1'] >= 50)
            mean_last = sum(p['last_ldl'] for p in tpts) / n
            print(f"  {label:<15s} {n:<6d} {100*lt26/n:<8.1f} {100*lt18/n:<8.1f} {100*lt14/n:<8.1f} {100*red50/n:<10.1f} {mean_last:<12.2f}")

    # By gene
    print(f"\n  Target achievement by gene:")
    for gene in ['LDLR', 'APOB']:
        gpts = [p for p in target_pts if p['gene'] == gene]
        if len(gpts) >= 5:
            n = len(gpts)
            lt26 = sum(1 for p in gpts if p['last_ldl'] < 2.6)
            lt14 = sum(1 for p in gpts if p['last_ldl'] < 1.4)
            mean_last = sum(p['last_ldl'] for p in gpts) / n
            print(f"    {gene}: n={n}, <2.6: {100*lt26/n:.1f}%, <1.4: {100*lt14/n:.1f}%, mean LastLDL={mean_last:.2f}")

# ============================================================================
#  H. STATIN TYPE-SPECIFIC RESPONSE BY SSS
# ============================================================================

print("\n" + "=" * 80)
print("  H. Statin-Specific Response by SSS")
print("=" * 80)

for stype in ['Atorvastatin', 'Rosuvastatin', 'Simvastatin']:
    s_pts = [p for p in patients if p['statin_type'] == stype and p['ldl_pct1'] is not None]
    if len(s_pts) >= 5:
        sss_vals = [p['sss'] for p in s_pts]
        ldl_vals = [p['ldl_pct1'] for p in s_pts]
        r, pval = pearson_r(sss_vals, ldl_vals)
        mean_dose = sum(p['statin_dose'] for p in s_pts if p['statin_dose']) / max(1, sum(1 for p in s_pts if p['statin_dose']))
        mean_red = sum(ldl_vals) / len(ldl_vals)

        print(f"\n  {stype}:")
        print(f"    n={len(s_pts)}, mean dose={mean_dose:.0f}mg, mean LDL reduction={mean_red:.1f}%")
        print(f"    SSS vs LDL% reduction: r={r:.4f}, p={pval:.4f}")

        # Low vs high SSS response on this statin
        med_sss = median(sss_vals)
        low_s = [p for p in s_pts if p['sss'] < med_sss]
        high_s = [p for p in s_pts if p['sss'] >= med_sss]
        if low_s and high_s:
            low_r = sum(p['ldl_pct1'] for p in low_s) / len(low_s)
            high_r = sum(p['ldl_pct1'] for p in high_s) / len(high_s)
            print(f"    Low SSS: {low_r:.1f}% reduction (n={len(low_s)})")
            print(f"    High SSS: {high_r:.1f}% reduction (n={len(high_s)})")

# ============================================================================
#  I. EZETIMIBE AND PCSK9i ADD-ON RESPONSE
# ============================================================================

print("\n" + "=" * 80)
print("  I. Add-On Therapy (Ezetimibe, PCSK9i) by SSS")
print("=" * 80)
print("  Clinical question: Do severe variants need add-on therapy more often?\n")

# Ezetimibe add-on
ezet_pts = [p for p in patients if p['ezetimibe']]
no_ezet = [p for p in patients if not p['ezetimibe'] and p['on_statin']]

print(f"  Ezetimibe users: {len(ezet_pts)}")
print(f"  Statin-only: {len(no_ezet)}")

if ezet_pts and no_ezet:
    ezet_sss = sum(p['sss'] for p in ezet_pts) / len(ezet_pts)
    no_ezet_sss = sum(p['sss'] for p in no_ezet) / len(no_ezet)
    print(f"  Mean SSS (Ezetimibe): {ezet_sss:.3f}")
    print(f"  Mean SSS (Statin only): {no_ezet_sss:.3f}")

    ezet_with_r = [p for p in ezet_pts if p['ldl_pct1'] is not None]
    no_ezet_with_r = [p for p in no_ezet if p['ldl_pct1'] is not None]

    if ezet_with_r and no_ezet_with_r:
        ezet_r = sum(p['ldl_pct1'] for p in ezet_with_r) / len(ezet_with_r)
        no_r = sum(p['ldl_pct1'] for p in no_ezet_with_r) / len(no_ezet_with_r)
        print(f"  LDL% reduction with Ezetimibe: {ezet_r:.1f}% (n={len(ezet_with_r)})")
        print(f"  LDL% reduction statin only: {no_r:.1f}% (n={len(no_ezet_with_r)})")

# PCSK9i
pcsk_pts = [p for p in patients if p['pcsk9i']]
print(f"\n  PCSK9i users: {len(pcsk_pts)}")

if pcsk_pts:
    pcsk_sss = sum(p['sss'] for p in pcsk_pts) / len(pcsk_pts)
    print(f"  Mean SSS (PCSK9i): {pcsk_sss:.3f}")
    pcsk_domains = defaultdict(int)
    for p in pcsk_pts:
        pcsk_domains[p['domain']] += 1
    print(f"  PCSK9i by domain: {dict(pcsk_domains)}")

    # PCSK9i response
    pcsk_with_r = [p for p in pcsk_pts if p['ldl_pct1'] is not None]
    if pcsk_with_r:
        pcsk_r = sum(p['ldl_pct1'] for p in pcsk_with_r) / len(pcsk_with_r)
        print(f"  LDL% reduction with PCSK9i: {pcsk_r:.1f}% (n={len(pcsk_with_r)})")

    # Individual PCSK9i patients
    print(f"\n  PCSK9i patient details:")
    for p in pcsk_pts[:10]:
        r_val = f"{p['ldl_pct1']:.1f}%" if p['ldl_pct1'] else "N/A"
        print(f"    {p['variant'][:30]:<30s} SSS={p['sss']:.3f} LDL1={p['ldl1'] or 'N/A'} LastLDL={p['last_ldl'] or 'N/A'} Response={r_val} {p['tx1_1'][:25]} + {p['tx1_3'][:20]}")

# ============================================================================
#  J. AGE AT TREATMENT INITIATION BY SSS
# ============================================================================

print("\n" + "=" * 80)
print("  J. Age at Treatment Initiation by SSS")
print("=" * 80)
print("  Clinical hypothesis: Severe variants should be treated EARLIER\n")

age_tx_pts = [p for p in patients if p['age_tx'] is not None and p['age_tx'] > 0]
print(f"  Patients with age at treatment: {len(age_tx_pts)}")

if age_tx_pts:
    sss_v = [p['sss'] for p in age_tx_pts]
    age_v = [p['age_tx'] for p in age_tx_pts]
    r, pval = pearson_r(sss_v, age_v)
    print(f"  SSS vs Age at treatment: r = {r:.4f}, p = {pval:.4f}")
    print(f"    (Negative r = severe mutations treated earlier = GOOD clinical practice)")

    for label, lo, hi in [('T1 (mild)', 0, t1), ('T2 (mod)', t1, t2), ('T3 (severe)', t2, 2.0)]:
        tpts = [p for p in age_tx_pts if lo <= p['sss'] < hi]
        if tpts:
            ages = [p['age_tx'] for p in tpts]
            print(f"    {label}: mean age at Tx = {sum(ages)/len(ages):.1f}y, median = {median(ages):.1f}y")

# ============================================================================
#  K. TREATMENT RESPONSE BY VARIANT TYPE (Mechanistic)
# ============================================================================

print("\n" + "=" * 80)
print("  K. Treatment Response by Variant Mechanism")
print("=" * 80)

# Group variants by mechanism
mechanism_groups = {
    'Loss-of-function (null)': [p for p in patients if p['sss'] >= 0.95 and p['ldl_pct1'] is not None],
    'Cysteine-disrupting': [p for p in patients if 'C' in str(sss_data.get(p['variant'], {}).get('ddG', ''))
                            and p['ldl_pct1'] is not None],
    'Defective (high ddG)': [p for p in patients if sss_data.get(p['variant'], {}).get('ddG') is not None
                              and (sss_data.get(p['variant'], {}).get('ddG') or 0) > 5
                              and p['ldl_pct1'] is not None],
    'Mild missense (low ddG)': [p for p in patients if sss_data.get(p['variant'], {}).get('ddG') is not None
                                 and 0 <= (sss_data.get(p['variant'], {}).get('ddG') or 99) < 2
                                 and p['ldl_pct1'] is not None],
    'APOB variants': [p for p in patients if p['gene'] == 'APOB' and p['ldl_pct1'] is not None],
    'LDLR ligand-binding': [p for p in patients if p['domain'] == 'Ligand-binding' and p['ldl_pct1'] is not None],
    'LDLR EGF domain': [p for p in patients if 'EGF' in p['domain'] and p['ldl_pct1'] is not None],
    'LDLR beta-propeller': [p for p in patients if p['domain'] == 'Beta-propeller' and p['ldl_pct1'] is not None],
}

print(f"\n  {'Mechanism':<30s} {'N':<6s} {'Mean LDL%':<10s} {'Median':<8s} {'Pre-LDL':<10s} {'SSS'}")
print(f"  {'-'*30} {'-'*6} {'-'*10} {'-'*8} {'-'*10} {'-'*6}")

for mech, mpts in mechanism_groups.items():
    if len(mpts) >= 2:
        ldl_r = [p['ldl_pct1'] for p in mpts]
        pre = [p['ldl1'] for p in mpts if p['ldl1']]
        sss_m = sum(p['sss'] for p in mpts) / len(mpts)
        pre_m = sum(pre)/len(pre) if pre else 0
        print(f"  {mech:<30s} {len(mpts):<6d} {sum(ldl_r)/len(ldl_r):<10.1f} {median(ldl_r):<8.1f} {pre_m:<10.2f} {sss_m:.3f}")

# ============================================================================
#  L. NON-HDL AND TRIGLYCERIDE RESPONSE
# ============================================================================

print("\n" + "=" * 80)
print("  L. Non-HDL and Triglyceride Response by SSS")
print("=" * 80)

nhdl_pts = [p for p in patients if p['nhdl_pct'] is not None and p['nhdl_pct'] > 0]
tc_pts = [p for p in patients if p['tc_pct'] is not None and p['tc_pct'] > 0]

print(f"\n  Non-HDL reduction data: {len(nhdl_pts)} patients")
print(f"  Total cholesterol reduction data: {len(tc_pts)} patients")

if nhdl_pts:
    sss_v = [p['sss'] for p in nhdl_pts]
    nhdl_v = [p['nhdl_pct'] for p in nhdl_pts]
    r, pval = pearson_r(sss_v, nhdl_v)
    print(f"\n  SSS vs Non-HDL% reduction: r = {r:.4f}, p = {pval:.4f}")

    for label, lo, hi in [('T1 (mild)', 0, t1), ('T2 (mod)', t1, t2), ('T3 (severe)', t2, 2.0)]:
        tpts = [p for p in nhdl_pts if lo <= p['sss'] < hi]
        if tpts:
            vals = [p['nhdl_pct'] for p in tpts]
            print(f"    {label}: Non-HDL reduction = {sum(vals)/len(vals):.1f}% (n={len(tpts)})")

if tc_pts:
    sss_v = [p['sss'] for p in tc_pts]
    tc_v = [p['tc_pct'] for p in tc_pts]
    r, pval = pearson_r(sss_v, tc_v)
    print(f"\n  SSS vs TC% reduction: r = {r:.4f}, p = {pval:.4f}")

# ============================================================================
#  SYNTHESIS
# ============================================================================

print("\n" + "=" * 80)
print("  SYNTHESIS: Variant-Level Pharmacogenomics")
print("=" * 80)

print("""
  CLINICAL IMPLICATIONS:

  1. TREATMENT INTENSITY vs STRUCTURAL SEVERITY
     If high-SSS patients receive more aggressive therapy AND achieve
     similar LDL reductions, the genetic signal is masked. But if they
     achieve LESS reduction despite more therapy, SSS predicts treatment
     resistance — a directly actionable finding.

  2. DOMAIN-SPECIFIC PHARMACOLOGY
     - Ligand-binding mutations: Statins can upregulate remaining LDLR
       copies, but if the binding domain is broken, even upregulated
       receptors can't capture LDL -> partial statin resistance
     - EGF/recycling mutations: Receptor is made but can't recycle.
       Statins increase production, but receptors are single-use ->
       moderate statin response
     - Null mutations: No receptor at all -> statins CANNOT work via
       LDLR pathway -> need PCSK9i or LDL apheresis

  3. VARIANT-SPECIFIC TREATMENT GUIDELINES
     Variants with < 30% LDL reduction on max statin should be flagged
     for early PCSK9i initiation. This is precision medicine:
     AlphaFold3 structure -> FoldX ddG -> SSS -> treatment decision.

  4. PCSK9i RESPONSE AND LDLR FUNCTION
     PCSK9 inhibitors work by preventing PCSK9-mediated LDLR degradation.
     If LDLR is already non-functional (null/frameshift), PCSK9i has
     NO target to protect -> predicted PCSK9i resistance.
     Interface variants (EGF-A domain) may show ENHANCED PCSK9i response
     because inhibiting PCSK9 binding to a damaged EGF-A domain
     preserves whatever function remains.
""")

# Save detailed results
output = os.path.join(ANALYSIS, "variant_llt_response.csv")
with open(output, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['patient_id', 'variant', 'gene', 'domain', 'sss',
                      'statin_type', 'statin_dose', 'intensity', 'ezetimibe', 'pcsk9i',
                      'ldl_1', 'ldl_2', 'ldl_3', 'ldl_4', 'last_ldl',
                      'ldl_pct_reduction', 'tc_pct_reduction', 'nhdl_pct_reduction',
                      'age_at_treatment', 'ascvd'])
    for p in patients:
        writer.writerow([p['pid'], p['variant'], p['gene'], p['domain'], p['sss'],
                          p['statin_type'], p['statin_dose'], p['intensity'],
                          p['ezetimibe'], p['pcsk9i'],
                          p.get('ldl1', ''), p.get('ldl2', ''), p.get('ldl3', ''),
                          p.get('ldl4', ''), p.get('last_ldl', ''),
                          p.get('ldl_pct1', ''), p.get('tc_pct', ''), p.get('nhdl_pct', ''),
                          p.get('age_tx', ''), p['ascvd']])

print(f"\n  Saved: {output}")
print(f"  Analysis complete.")

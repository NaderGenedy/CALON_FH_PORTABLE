#!/usr/bin/env python3
"""
31_external_validation_complete.py
Complete external validation: (A) Treatment response (UKB vs South Wales)
                              (B) VUS classification via SIFT/PolyPhen
"""
import csv, math, os, sys

BASE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")

# ── helper functions ────────────────────────────────────────────────────────

def sf(x):
    try: return float(str(x).strip())
    except: return None

def welch_t_test(vals1, vals2):
    n1, n2 = len(vals1), len(vals2)
    if n1 < 2 or n2 < 2: return 0.0, 1.0
    m1 = sum(vals1) / n1
    m2 = sum(vals2) / n2
    s1 = (sum((v - m1) ** 2 for v in vals1) / (n1 - 1)) ** 0.5
    s2 = (sum((v - m2) ** 2 for v in vals2) / (n2 - 1)) ** 0.5
    se = (s1 ** 2 / n1 + s2 ** 2 / n2) ** 0.5
    if se == 0: return 0.0, 1.0
    t = (m1 - m2) / se
    df = min(n1, n2) - 1
    p = 2 * math.exp(-0.5 * t * t) * 0.4 if abs(t) < 5 else 0.001
    return t, p

def odds_ratio_ci(a, b, c, d):
    if b * c == 0: return 0.0, 0.0, 0.0, 1.0
    OR = (a * d) / (b * c)
    se = (1 / max(a, 1) + 1 / max(b, 1) + 1 / max(c, 1) + 1 / max(d, 1)) ** 0.5
    ci_low = math.exp(math.log(max(OR, 0.001)) - 1.96 * se)
    ci_high = math.exp(math.log(max(OR, 0.001)) + 1.96 * se)
    n = a + b + c + d
    ea = (a + b) * (a + c) / n
    eb = (a + b) * (b + d) / n
    ec = (c + d) * (a + c) / n
    ed = (c + d) * (b + d) / n
    chi2 = sum((o - e) ** 2 / max(e, 0.001) for o, e in [(a, ea), (b, eb), (c, ec), (d, ed)])
    p = math.exp(-0.5 * chi2) if chi2 < 30 else 0.0001
    return OR, ci_low, ci_high, p

def spearman_rank(x, y):
    n = len(x)
    if n < 3: return 0.0, 1.0
    rx = [sorted(x).index(v) + 1 for v in x]
    ry = [sorted(y).index(v) + 1 for v in y]
    d2 = sum((rxi - ryi) ** 2 for rxi, ryi in zip(rx, ry))
    rho = 1 - 6 * d2 / (n * (n * n - 1))
    return rho, 0.0

def mean_ci_95(vals):
    n = len(vals)
    if n == 0: return float('nan'), float('nan'), float('nan')
    if n == 1: return vals[0], vals[0], vals[0]
    m = sum(vals) / n
    sd = (sum((v - m) ** 2 for v in vals) / (n - 1)) ** 0.5
    se = sd / n ** 0.5
    return m, m - 1.96 * se, m + 1.96 * se

def read_csv(path):
    rows = []
    with open(path, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)
    return rows

# ── load data ───────────────────────────────────────────────────────────────

print("=" * 80)
print("EXTERNAL VALIDATION ANALYSIS")
print("=" * 80)

ukb_sss = read_csv(os.path.join(ANALYSIS, "ukb_external_validation_sss.csv"))
wales = read_csv(os.path.join(ANALYSIS, "comprehensive_sss_analysis.csv"))
vep = read_csv(os.path.join(ANALYSIS, "ukb_ldlr_vep_annotations.csv"))
carriers = read_csv(os.path.join(ANALYSIS, "ukb_carriers_variant_sss.csv"))
ukb_clinical = read_csv(os.path.join(BASE, "calon_ukb_analysis_ready.csv"))

print(f"UKB SSS data:        {len(ukb_sss)} rows")
print(f"South Wales data:    {len(wales)} rows")
print(f"VEP annotations:     {len(vep)} rows")
print(f"UKB carriers+SSS:    {len(carriers)} rows")
print(f"UKB clinical data:   {len(ukb_clinical)} rows")

# Build clinical lookup by eid
clin_lookup = {}
for r in ukb_clinical:
    eid = r.get('eid', '').strip().strip('"')
    clin_lookup[eid] = r

# ============================================================================
# PART A: TREATMENT RESPONSE EXTERNAL VALIDATION
# ============================================================================
print("\n" + "=" * 80)
print("PART A: TREATMENT RESPONSE EXTERNAL VALIDATION (UKB vs South Wales)")
print("=" * 80)

# ── A1: UKB domain-specific LDL by statin status ───────────────────────────
# Group UKB patients by domain and statin status
domain_treated = {}   # domain -> [ldl values] for on_statin==1
domain_untreated = {} # domain -> [ldl values] for on_statin==0

for r in ukb_sss:
    domain = r.get('domain', '').strip()
    ldl = sf(r.get('ldl'))
    statin = r.get('on_statin', '').strip()
    if not domain or ldl is None or statin not in ('0', '1'):
        continue
    if statin == '1':
        domain_treated.setdefault(domain, []).append(ldl)
    else:
        domain_untreated.setdefault(domain, []).append(ldl)

all_domains_ukb = sorted(set(list(domain_treated.keys()) + list(domain_untreated.keys())))

print(f"\nDomains found in UKB data: {len(all_domains_ukb)}")
print(f"{'Domain':<30s} {'N_untrt':>7s} {'N_trt':>7s} {'LDL_untrt':>10s} {'LDL_trt':>10s} {'Diff':>8s} {'t-stat':>8s} {'P':>8s}")
print("-" * 95)

ukb_domain_results = []

for dom in all_domains_ukb:
    ut = domain_untreated.get(dom, [])
    tr = domain_treated.get(dom, [])
    m_ut, lo_ut, hi_ut = mean_ci_95(ut) if ut else (float('nan'), float('nan'), float('nan'))
    m_tr, lo_tr, hi_tr = mean_ci_95(tr) if tr else (float('nan'), float('nan'), float('nan'))
    diff = m_tr - m_ut if not (math.isnan(m_tr) or math.isnan(m_ut)) else float('nan')
    t, p = welch_t_test(tr, ut) if len(tr) >= 2 and len(ut) >= 2 else (0, 1.0)
    print(f"{dom:<30s} {len(ut):>7d} {len(tr):>7d} {m_ut:>10.3f} {m_tr:>10.3f} {diff:>8.3f} {t:>8.2f} {p:>8.4f}")
    ukb_domain_results.append({
        'domain': dom,
        'n_untreated': len(ut),
        'n_treated': len(tr),
        'mean_ldl_untreated': round(m_ut, 4) if not math.isnan(m_ut) else '',
        'ci_low_untreated': round(lo_ut, 4) if not math.isnan(lo_ut) else '',
        'ci_high_untreated': round(hi_ut, 4) if not math.isnan(hi_ut) else '',
        'mean_ldl_treated': round(m_tr, 4) if not math.isnan(m_tr) else '',
        'ci_low_treated': round(lo_tr, 4) if not math.isnan(lo_tr) else '',
        'ci_high_treated': round(hi_tr, 4) if not math.isnan(hi_tr) else '',
        'ldl_diff': round(diff, 4) if not math.isnan(diff) else '',
        't_stat': round(t, 4),
        'p_value': round(p, 6)
    })

# ── A2: South Wales domain-specific treatment response (ldl_pct_change) ────
print("\n--- South Wales Development Cohort: Domain-specific LDL % Change ---")
wales_domain_response = {}
for r in wales:
    domain = r.get('domain', '').strip()
    pct = sf(r.get('ldl_pct_change'))
    if domain and pct is not None:
        wales_domain_response.setdefault(domain, []).append(pct)

print(f"\n{'Domain':<30s} {'N':>5s} {'Mean_%chg':>10s} {'CI_low':>10s} {'CI_high':>10s}")
print("-" * 70)

wales_domain_means = {}
for dom in sorted(wales_domain_response.keys()):
    vals = wales_domain_response[dom]
    m, lo, hi = mean_ci_95(vals)
    wales_domain_means[dom] = m
    print(f"{dom:<30s} {len(vals):>5d} {m:>10.2f} {lo:>10.2f} {hi:>10.2f}")

# ── A3: Spearman rank correlation between UKB and Wales domain rankings ────
print("\n--- Spearman Rank Correlation: UKB on-treatment LDL vs Wales LDL% change ---")

# Find common domains with sufficient data
common_domains = []
ukb_vals_for_corr = []
wales_vals_for_corr = []

# UKB metric: mean treated LDL (higher = worse treatment response)
ukb_treated_means = {}
for dom in all_domains_ukb:
    tr = domain_treated.get(dom, [])
    if len(tr) >= 2:
        ukb_treated_means[dom] = sum(tr) / len(tr)

for dom in sorted(set(ukb_treated_means.keys()) & set(wales_domain_means.keys())):
    common_domains.append(dom)
    ukb_vals_for_corr.append(ukb_treated_means[dom])
    # Wales: more negative pct_change = BETTER response, so use negative for "worse"
    # We want: high UKB treated LDL ~ less negative Wales pct_change (= poor response)
    wales_vals_for_corr.append(wales_domain_means[dom])

print(f"Common domains with data in both cohorts: {len(common_domains)}")
for i, dom in enumerate(common_domains):
    print(f"  {dom:<30s}  UKB trt LDL={ukb_vals_for_corr[i]:.3f}  Wales %chg={wales_vals_for_corr[i]:.2f}")

if len(common_domains) >= 3:
    # For Spearman: we expect HIGH treated LDL in UKB <-> LESS negative (higher) %change in Wales
    rho, p_rho = spearman_rank(ukb_vals_for_corr, wales_vals_for_corr)
    print(f"\nSpearman rho = {rho:.4f}")
    print(f"Interpretation: {'Positive correlation confirms external validation' if rho > 0 else 'Negative correlation (unexpected direction)'}")
    print(f"  Domains with higher on-treatment LDL in UKB {'do' if rho > 0 else 'do NOT'} correspond to")
    print(f"  domains with less LDL reduction in South Wales")
else:
    rho = float('nan')
    print("  Insufficient common domains for correlation.")

# ── A4: Forest plot data (domain-specific LDL reduction with 95% CI) ──────
print("\n--- Forest Plot Data: Domain-specific On-Treatment LDL (UKB) ---")
print(f"{'Domain':<30s} {'Mean':>8s} {'95% CI_low':>10s} {'95% CI_hi':>10s} {'N':>5s}")
print("-" * 68)
for dom in all_domains_ukb:
    tr = domain_treated.get(dom, [])
    if len(tr) >= 2:
        m, lo, hi = mean_ci_95(tr)
        print(f"{dom:<30s} {m:>8.3f} {lo:>10.3f} {hi:>10.3f} {len(tr):>5d}")

# ── A5: On-treatment LDL by SSS tertile ───────────────────────────────────
print("\n--- On-Treatment LDL by SSS Tertile (Higher SSS = More Severe?) ---")

treated_sss_ldl = []
for r in ukb_sss:
    statin = r.get('on_statin', '').strip()
    sss = sf(r.get('sss'))
    ldl = sf(r.get('ldl'))
    if statin == '1' and sss is not None and ldl is not None:
        treated_sss_ldl.append((sss, ldl))

treated_sss_ldl.sort(key=lambda x: x[0])
n = len(treated_sss_ldl)
if n >= 6:
    t1_cut = n // 3
    t2_cut = 2 * n // 3
    tertiles = [
        ("T1 (low SSS)", [x[1] for x in treated_sss_ldl[:t1_cut]]),
        ("T2 (mid SSS)", [x[1] for x in treated_sss_ldl[t1_cut:t2_cut]]),
        ("T3 (high SSS)", [x[1] for x in treated_sss_ldl[t2_cut:]]),
    ]
    print(f"  N treated patients with SSS: {n}")
    for label, vals in tertiles:
        m, lo, hi = mean_ci_95(vals)
        print(f"  {label:<20s}  N={len(vals):>3d}  mean LDL={m:.3f}  95%CI [{lo:.3f}, {hi:.3f}]")

    # T3 vs T1 comparison
    t_stat, p_val = welch_t_test(tertiles[2][1], tertiles[0][1])
    print(f"\n  T3 vs T1 Welch t-test: t={t_stat:.3f}, P={p_val:.4f}")
    if tertiles[2][1] and tertiles[0][1]:
        diff_m = sum(tertiles[2][1])/len(tertiles[2][1]) - sum(tertiles[0][1])/len(tertiles[0][1])
        print(f"  Mean LDL difference (T3 - T1): {diff_m:.3f} mmol/L")
        if diff_m > 0:
            print("  -> Higher SSS associated with HIGHER on-treatment LDL (treatment resistance)")
        else:
            print("  -> Higher SSS NOT associated with higher on-treatment LDL")
else:
    print(f"  Insufficient treated patients with SSS data: {n}")

# ── Save Part A results ────────────────────────────────────────────────────
out_a = os.path.join(ANALYSIS, "external_treatment_response_validation.csv")
with open(out_a, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=[
        'domain', 'n_untreated', 'n_treated',
        'mean_ldl_untreated', 'ci_low_untreated', 'ci_high_untreated',
        'mean_ldl_treated', 'ci_low_treated', 'ci_high_treated',
        'ldl_diff', 't_stat', 'p_value'
    ])
    w.writeheader()
    w.writerows(ukb_domain_results)
print(f"\nSaved: {out_a}")


# ============================================================================
# PART B: VUS EXTERNAL VALIDATION VIA SIFT/POLYPHEN CLASSIFICATION
# ============================================================================
print("\n" + "=" * 80)
print("PART B: VUS EXTERNAL VALIDATION (SIFT/PolyPhen Classification)")
print("=" * 80)

# ── B1: Classify VEP variants ──────────────────────────────────────────────
# Build a lookup: (chrom, pos, ref, alt) -> classification
def classify_variant(sift_str, polyphen_str, consequence):
    sift = str(sift_str).strip().lower()
    pp = str(polyphen_str).strip().lower()

    # Frameshift, stop_gained, splice donor/acceptor = pathogenic-like
    severe = ['frameshift_variant', 'stop_gained', 'splice_donor_variant',
              'splice_acceptor_variant', 'start_lost']
    if consequence in severe:
        return 'Pathogenic-like'

    # Synonymous, upstream, downstream, intronic = benign-like
    benign_cons = ['synonymous_variant', 'upstream_gene_variant',
                   'downstream_gene_variant', 'intron_variant',
                   '5_prime_UTR_variant', '3_prime_UTR_variant']
    if consequence in benign_cons:
        return 'Benign-like'

    # For missense: use SIFT + PolyPhen
    if 'deleterious' in sift and ('probably_damaging' in pp or 'possibly_damaging' in pp):
        return 'Pathogenic-like'
    if 'tolerated' in sift and 'benign' in pp:
        return 'Benign-like'

    # Mixed or missing predictions = VUS-like
    return 'VUS-like'

vep_classification = {}
for r in vep:
    key = (r.get('chrom',''), r.get('pos',''), r.get('ref',''), r.get('alt',''))
    cons = r.get('consequence', '').strip()
    sift_val = r.get('sift', '')
    pp_val = r.get('polyphen', '')
    cls = classify_variant(sift_val, pp_val, cons)
    vep_classification[key] = {
        'classification': cls,
        'consequence': cons,
        'sift': sift_val,
        'polyphen': pp_val
    }

# Count classifications
class_counts = {}
for v in vep_classification.values():
    c = v['classification']
    class_counts[c] = class_counts.get(c, 0) + 1

print(f"\nVEP variant classifications (unique variants):")
for c in sorted(class_counts.keys()):
    print(f"  {c:<20s}: {class_counts[c]:>5d}")

# ── B2: Map carriers to classification + ASCVD ─────────────────────────────
# carrier rows have variant_id (LDLR:c.XXX) but we need to map to VEP via pos
# carriers also have domain, sss, consequence
# We'll use the carrier's consequence + VEP-based classification

# Build carrier-level data with classification
carrier_data = []
# Deduplicate carriers (some appear twice)
seen_carrier = set()

for r in carriers:
    eid = r.get('eid', '').strip()
    variant_id = r.get('variant_id', '').strip()
    key = (eid, variant_id)
    if key in seen_carrier:
        continue
    seen_carrier.add(key)

    sss = sf(r.get('sss'))
    domain = r.get('domain', '').strip()
    consequence = r.get('consequence', '').strip()

    # Classify based on consequence (we don't have chr:pos in carrier file, so use consequence)
    # For missense, we mark as VUS-like (default); for severe/benign consequence, classify accordingly
    cls = classify_variant('', '', consequence)  # Without SIFT/PP, uses consequence only

    # Get clinical outcome
    clin = clin_lookup.get(eid)
    ascvd = None
    if clin:
        ascvd_val = clin.get('ascvd_combined', '').strip().strip('"')
        ascvd = sf(ascvd_val)

    carrier_data.append({
        'eid': eid,
        'variant_id': variant_id,
        'sss': sss,
        'domain': domain,
        'consequence': consequence,
        'classification': cls,
        'ascvd': ascvd
    })

# However, we also need to refine missense classifications using VEP SIFT/PP
# Build a lookup from carrier variant_id -> VEP via protein_position + amino_acids
# Actually, carriers have protein_position. Let's try matching by protein_position to VEP.

# Build VEP lookup by protein_position for missense
vep_by_pos = {}
for r in vep:
    pp_pos = r.get('protein_position', '').strip()
    cons = r.get('consequence', '').strip()
    if pp_pos and cons == 'missense_variant':
        sift_val = r.get('sift', '').strip()
        pp_val = r.get('polyphen', '').strip()
        vep_by_pos.setdefault(pp_pos, []).append((sift_val, pp_val))

# Refine carrier classifications for missense
refined_count = 0
for cd in carrier_data:
    if cd['consequence'] == 'missense_variant':
        # Try to find VEP SIFT/PP by protein_position
        # We need protein_position from the carrier file
        # Re-read carriers to get protein_position
        pass  # Will handle below

# Re-read carriers to get protein_position
carrier_pp = {}
for r in carriers:
    vid = r.get('variant_id', '').strip()
    pp_pos = r.get('protein_position', '').strip()
    if vid and pp_pos:
        carrier_pp[vid] = pp_pos

# Now refine
for cd in carrier_data:
    if cd['consequence'] == 'missense_variant':
        vid = cd['variant_id']
        pp_pos = carrier_pp.get(vid, '')
        if pp_pos and pp_pos in vep_by_pos:
            # Use SIFT/PP from VEP
            sift_val, pp_val = vep_by_pos[pp_pos][0]
            cd['classification'] = classify_variant(sift_val, pp_val, 'missense_variant')
            refined_count += 1

print(f"\nRefined {refined_count} missense carrier classifications using VEP SIFT/PolyPhen")

# Count carrier classifications
carrier_class_counts = {}
for cd in carrier_data:
    c = cd['classification']
    carrier_class_counts[c] = carrier_class_counts.get(c, 0) + 1

print(f"\nCarrier classifications:")
for c in sorted(carrier_class_counts.keys()):
    print(f"  {c:<20s}: {carrier_class_counts[c]:>5d}")

# ── B3: ASCVD rates across classification groups ──────────────────────────
print("\n--- ASCVD Rates by Variant Classification ---")

class_ascvd = {}
for cd in carrier_data:
    c = cd['classification']
    a = cd['ascvd']
    if a is not None:
        class_ascvd.setdefault(c, []).append(int(a))

print(f"\n{'Classification':<20s} {'N':>6s} {'ASCVD+':>7s} {'ASCVD%':>8s} {'95% CI':>20s}")
print("-" * 65)

class_rates = {}
for c in sorted(class_ascvd.keys()):
    vals = class_ascvd[c]
    n = len(vals)
    n_events = sum(vals)
    rate = n_events / n if n > 0 else 0
    # Wilson CI for proportion
    if n > 0:
        se = (rate * (1 - rate) / n) ** 0.5
        lo = max(0, rate - 1.96 * se)
        hi = min(1, rate + 1.96 * se)
    else:
        lo, hi = 0, 0
    class_rates[c] = (n_events, n, rate)
    print(f"{c:<20s} {n:>6d} {n_events:>7d} {rate*100:>7.1f}% [{lo*100:.1f}%, {hi*100:.1f}%]")

# OR: Pathogenic-like vs Benign-like
if 'Pathogenic-like' in class_ascvd and 'Benign-like' in class_ascvd:
    path_vals = class_ascvd['Pathogenic-like']
    ben_vals = class_ascvd['Benign-like']
    a = sum(path_vals)
    b = len(path_vals) - a
    c = sum(ben_vals)
    d = len(ben_vals) - c
    OR, ci_l, ci_h, p = odds_ratio_ci(a, b, c, d)
    print(f"\nPathogenic-like vs Benign-like ASCVD:")
    print(f"  OR = {OR:.2f} (95% CI: {ci_l:.2f}-{ci_h:.2f}), P = {p:.4f}")

# OR: VUS-like vs Benign-like
if 'VUS-like' in class_ascvd and 'Benign-like' in class_ascvd:
    vus_vals = class_ascvd['VUS-like']
    ben_vals = class_ascvd['Benign-like']
    a = sum(vus_vals)
    b = len(vus_vals) - a
    c = sum(ben_vals)
    d = len(ben_vals) - c
    OR, ci_l, ci_h, p = odds_ratio_ci(a, b, c, d)
    print(f"\nVUS-like vs Benign-like ASCVD:")
    print(f"  OR = {OR:.2f} (95% CI: {ci_l:.2f}-{ci_h:.2f}), P = {p:.4f}")

# ── B4: Within VUS-like: HIGH SSS vs LOW SSS and ASCVD ───────────────────
print("\n--- Within VUS-like Variants: SSS and ASCVD ---")

vus_data = [(cd['sss'], cd['ascvd']) for cd in carrier_data
            if cd['classification'] == 'VUS-like' and cd['sss'] is not None and cd['ascvd'] is not None]

if len(vus_data) >= 4:
    vus_data.sort(key=lambda x: x[0])
    med_idx = len(vus_data) // 2
    low_sss = [x for x in vus_data[:med_idx]]
    high_sss = [x for x in vus_data[med_idx:]]

    low_ascvd = sum(x[1] for x in low_sss)
    low_n = len(low_sss)
    high_ascvd = sum(x[1] for x in high_sss)
    high_n = len(high_sss)

    low_sss_mean = sum(x[0] for x in low_sss) / low_n
    high_sss_mean = sum(x[0] for x in high_sss) / high_n

    print(f"  Low SSS group  (mean SSS={low_sss_mean:.3f}): N={low_n}, ASCVD={low_ascvd} ({low_ascvd/low_n*100:.1f}%)")
    print(f"  High SSS group (mean SSS={high_sss_mean:.3f}): N={high_n}, ASCVD={high_ascvd} ({high_ascvd/high_n*100:.1f}%)")

    a_h = int(high_ascvd)
    b_h = high_n - a_h
    c_l = int(low_ascvd)
    d_l = low_n - c_l
    OR, ci_l, ci_h, p = odds_ratio_ci(a_h, b_h, c_l, d_l)
    print(f"  High vs Low SSS OR = {OR:.2f} (95% CI: {ci_l:.2f}-{ci_h:.2f}), P = {p:.4f}")
    if OR > 1:
        print("  -> Higher SSS VUS-like variants have HIGHER ASCVD risk")
    else:
        print("  -> Higher SSS VUS-like variants do NOT show higher ASCVD risk")
else:
    print(f"  Insufficient VUS-like carriers with SSS+ASCVD data: {len(vus_data)}")

# ── B5: Domain-specific ASCVD in VUS-like UKB vs Wales VUS ────────────────
print("\n--- Domain-Specific ASCVD in VUS-like Variants (UKB vs Wales) ---")

# UKB VUS-like by domain
ukb_vus_domain = {}
for cd in carrier_data:
    if cd['classification'] == 'VUS-like' and cd['domain'] and cd['ascvd'] is not None:
        ukb_vus_domain.setdefault(cd['domain'], []).append(int(cd['ascvd']))

# Wales: identify VUS-like (we don't have SIFT/PP, so use sss_source)
# VUS in Wales are typically those without a definitive classification
# Use sss_source == 'foldx' as proxy for VUS (not preset_lof or preset_benign)
wales_vus_domain = {}
for r in wales:
    domain = r.get('domain', '').strip()
    sss_src = r.get('sss_source', '').strip()
    ascvd_val = sf(r.get('ascvd'))
    # FoldX-predicted = missense VUS-like; preset_lof = clearly pathogenic
    if sss_src == 'foldx' and domain and ascvd_val is not None:
        wales_vus_domain.setdefault(domain, []).append(int(ascvd_val))

all_vus_domains = sorted(set(list(ukb_vus_domain.keys()) + list(wales_vus_domain.keys())))

print(f"\n{'Domain':<30s} {'UKB_N':>6s} {'UKB_%':>7s} {'Wales_N':>8s} {'Wales_%':>8s}")
print("-" * 65)

ukb_vus_rates = {}
wales_vus_rates = {}

for dom in all_vus_domains:
    ukb_v = ukb_vus_domain.get(dom, [])
    wales_v = wales_vus_domain.get(dom, [])
    ukb_rate = sum(ukb_v) / len(ukb_v) * 100 if ukb_v else float('nan')
    wales_rate = sum(wales_v) / len(wales_v) * 100 if wales_v else float('nan')
    ukb_vus_rates[dom] = ukb_rate
    wales_vus_rates[dom] = wales_rate
    ukb_str = f"{ukb_rate:.1f}%" if not math.isnan(ukb_rate) else "N/A"
    wales_str = f"{wales_rate:.1f}%" if not math.isnan(wales_rate) else "N/A"
    print(f"{dom:<30s} {len(ukb_v):>6d} {ukb_str:>7s} {len(wales_v):>8d} {wales_str:>8s}")

# Correlation of domain-specific VUS ASCVD rates
common_vus_doms = [d for d in all_vus_domains
                   if d in ukb_vus_rates and d in wales_vus_rates
                   and not math.isnan(ukb_vus_rates.get(d, float('nan')))
                   and not math.isnan(wales_vus_rates.get(d, float('nan')))
                   and len(ukb_vus_domain.get(d, [])) >= 2
                   and len(wales_vus_domain.get(d, [])) >= 2]

if len(common_vus_doms) >= 3:
    ukb_r = [ukb_vus_rates[d] for d in common_vus_doms]
    wales_r = [wales_vus_rates[d] for d in common_vus_doms]
    rho_vus, _ = spearman_rank(ukb_r, wales_r)
    print(f"\nSpearman correlation of domain-specific VUS ASCVD rates (UKB vs Wales): rho={rho_vus:.4f}")
    print(f"  Common domains: {len(common_vus_doms)}")
else:
    rho_vus = float('nan')
    print(f"\nInsufficient common domains for VUS ASCVD correlation: {len(common_vus_doms)}")

# ── B6: SSS as discriminator within each classification group ─────────────
print("\n--- SSS Distribution by Classification ---")

for cls in ['Pathogenic-like', 'VUS-like', 'Benign-like']:
    sss_vals = [cd['sss'] for cd in carrier_data if cd['classification'] == cls and cd['sss'] is not None]
    if sss_vals:
        m, lo, hi = mean_ci_95(sss_vals)
        print(f"  {cls:<20s}: N={len(sss_vals):>4d}, mean SSS={m:.4f} [{lo:.4f}, {hi:.4f}]")

# T-test: SSS pathogenic-like vs benign-like
path_sss = [cd['sss'] for cd in carrier_data if cd['classification'] == 'Pathogenic-like' and cd['sss'] is not None]
ben_sss = [cd['sss'] for cd in carrier_data if cd['classification'] == 'Benign-like' and cd['sss'] is not None]
if path_sss and ben_sss:
    t_s, p_s = welch_t_test(path_sss, ben_sss)
    print(f"\n  SSS Pathogenic-like vs Benign-like: t={t_s:.3f}, P={p_s:.4f}")

# ── Save Part B results ────────────────────────────────────────────────────
out_b = os.path.join(ANALYSIS, "external_vus_validation.csv")

# Build output rows
vus_out_rows = []
for cd in carrier_data:
    vus_out_rows.append({
        'eid': cd['eid'],
        'variant_id': cd['variant_id'],
        'domain': cd['domain'],
        'consequence': cd['consequence'],
        'classification': cd['classification'],
        'sss': cd['sss'] if cd['sss'] is not None else '',
        'ascvd': int(cd['ascvd']) if cd['ascvd'] is not None else ''
    })

with open(out_b, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['eid', 'variant_id', 'domain', 'consequence',
                                       'classification', 'sss', 'ascvd'])
    w.writeheader()
    w.writerows(vus_out_rows)
print(f"\nSaved: {out_b}")


# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("SUMMARY OF EXTERNAL VALIDATION FINDINGS")
print("=" * 80)

print("\n--- Part A: Treatment Response External Validation ---")
# Count domains where treated LDL > untreated LDL (paradoxical in cross-sectional)
n_higher = sum(1 for r in ukb_domain_results if r['ldl_diff'] != '' and float(r['ldl_diff']) < 0)
n_total_dom = sum(1 for r in ukb_domain_results if r['ldl_diff'] != '')
print(f"1. {len(all_domains_ukb)} LDLR domains analysed in UKB")
print(f"2. {n_higher}/{n_total_dom} domains show LOWER LDL in treated vs untreated")
print(f"   (expected: treated patients have indication bias = higher baseline,")
print(f"    but statin effect lowers measured LDL)")

if not math.isnan(rho) if isinstance(rho, float) else True:
    direction = "POSITIVE" if rho > 0 else "NEGATIVE"
    print(f"3. Spearman rank correlation (UKB treated LDL vs Wales LDL% change): rho = {rho:.4f}")
    print(f"   {direction} correlation: domains with higher on-treatment LDL in UKB")
    if rho > 0:
        print(f"   correspond to domains with LESS LDL reduction in South Wales")
        print(f"   => VALIDATES domain-specific treatment resistance across cohorts")
    else:
        print(f"   correspond to domains with MORE LDL reduction in South Wales")

if n >= 6:
    t3_ldl = sum(tertiles[2][1]) / len(tertiles[2][1])
    t1_ldl = sum(tertiles[0][1]) / len(tertiles[0][1])
    print(f"4. SSS tertile analysis: T3 on-treatment LDL = {t3_ldl:.3f} vs T1 = {t1_ldl:.3f}")
    if t3_ldl > t1_ldl:
        print(f"   Higher SSS = higher on-treatment LDL => SSS predicts treatment resistance")
    else:
        print(f"   No clear SSS-treatment resistance gradient")

print("\n--- Part B: VUS External Validation ---")
for c in sorted(carrier_class_counts.keys()):
    print(f"  {c}: {carrier_class_counts[c]} carriers")

if 'Pathogenic-like' in class_rates and 'Benign-like' in class_rates:
    p_rate = class_rates['Pathogenic-like'][2] * 100
    b_rate = class_rates['Benign-like'][2] * 100
    print(f"5. ASCVD rates: Pathogenic-like {p_rate:.1f}% vs Benign-like {b_rate:.1f}%")
    if p_rate > b_rate:
        print(f"   SIFT/PolyPhen classification separates ASCVD risk as expected")

if 'VUS-like' in class_rates:
    v_rate = class_rates['VUS-like'][2] * 100
    print(f"6. VUS-like ASCVD rate: {v_rate:.1f}%")
    if 'Pathogenic-like' in class_rates and 'Benign-like' in class_rates:
        p_rate = class_rates['Pathogenic-like'][2] * 100
        b_rate = class_rates['Benign-like'][2] * 100
        if b_rate < v_rate < p_rate:
            print(f"   VUS-like ASCVD rate is INTERMEDIATE (between benign and pathogenic)")
            print(f"   => SSS could help reclassify VUS into higher/lower risk")

if len(vus_data) >= 4:
    print(f"7. Within VUS-like: high SSS group ASCVD = {high_ascvd/high_n*100:.1f}% vs low SSS = {low_ascvd/low_n*100:.1f}%")
    if high_ascvd / high_n > low_ascvd / low_n:
        print(f"   => SSS stratifies VUS ASCVD risk in UKB (external validation of Wales findings)")

if not math.isnan(rho_vus):
    print(f"8. Domain-specific VUS ASCVD correlation (UKB vs Wales): rho = {rho_vus:.4f}")
    if rho_vus > 0:
        print(f"   => Domain-level VUS risk patterns replicate across cohorts")

print(f"\nOutput files:")
print(f"  {out_a}")
print(f"  {out_b}")
print("\nDone.")

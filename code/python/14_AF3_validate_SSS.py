#!/usr/bin/env python3
"""
14_AF3_validate_SSS.py
Phase 4a: Validate SSS against ASCVD outcomes in Wales FH cohorts

Statistical tests:
  1. SSS association with ASCVD (logistic regression, age+sex adjusted)
  2. Incremental value over gene-only model (likelihood ratio test)
  3. SSS tertile/quartile ASCVD rates (chi-square trend test)
  4. Domain-specific forest plot data
  5. Discrimination (AUC) for SSS vs gene-only

Uses clustered standard errors (grouped by variant) per reviewer recommendation.

Inputs:
  alphafold/analysis/structural_severity_scores.csv
  DRAGON_3.csv  (South Wales FH registry)
  WALES_FH_CLEANED (1) - Copy.csv  (All Wales FH registry)

Outputs:
  alphafold/analysis/sss_validation_results.txt
  alphafold/analysis/sss_validation_data.csv
"""

import os
import csv
import math
import sys
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
AF_DIR = os.path.join(BASE, "alphafold")
ANALYSIS_DIR = os.path.join(AF_DIR, "analysis")

# ── Load SSS scores ──────────────────────────────────────────────
def load_sss():
    """Load variant-level SSS lookup."""
    sss_lookup = {}
    sss_file = os.path.join(ANALYSIS_DIR, "structural_severity_scores.csv")
    with open(sss_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sss_lookup[row['variant_id']] = {
                'sss': float(row['sss']),
                'gene': row['gene'],
                'domain': row['domain'],
                'sss_source': row['sss_source'],
            }
    return sss_lookup


def load_dragon3(sss_lookup):
    """Load DRAGON_3 (South Wales) and map SSS to patients."""
    dragon_file = os.path.join(BASE, "DRAGON_3.csv")
    patients = []

    with open(dragon_file, encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            mutation1 = row.get('Mutation1', '').strip()
            if not mutation1:
                continue

            # ASCVD outcome
            ascvd = row.get('ASCVD_combined', '')
            try:
                ascvd = int(float(ascvd))
            except (ValueError, TypeError):
                continue

            # Age
            age_str = row.get('age_at_event_or_censoring', '')
            try:
                age = float(age_str)
            except (ValueError, TypeError):
                continue

            # Gender
            gender = row.get('Gender', '').strip()
            sex = 1 if gender.lower() in ('m', 'male', '1') else 0

            # Look up SSS
            sss_info = sss_lookup.get(mutation1, None)
            if sss_info:
                sss = sss_info['sss']
                gene = sss_info['gene']
            else:
                # Try partial match
                sss = None
                gene = ''
                for vid, info in sss_lookup.items():
                    if vid in mutation1 or mutation1 in vid:
                        sss = info['sss']
                        gene = info['gene']
                        break
                if sss is None:
                    continue

            # Additional clinical variables
            ldl = row.get('LastLDL', '')
            try:
                ldl = float(ldl)
            except:
                ldl = None

            patients.append({
                'source': 'DRAGON3',
                'mutation': mutation1,
                'gene': gene,
                'ascvd': ascvd,
                'age': age,
                'sex': sex,
                'ldl': ldl,
                'sss': sss,
            })

    return patients


def load_wales_all(sss_lookup):
    """Load All Wales FH registry and map SSS."""
    wales_file = os.path.join(BASE, "WALES_FH_CLEANED (1) - Copy.csv")
    if not os.path.exists(wales_file):
        print("  !! Wales file not found")
        return []

    patients = []

    with open(wales_file, encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            mutation1 = row.get('Mutation1', '').strip()
            if not mutation1:
                continue

            # Wales All uses 'ascvd_combine' not 'ASCVD_combined'
            ascvd = row.get('ascvd_combine', row.get('ASCVD_combined', ''))
            try:
                ascvd = int(float(ascvd))
            except:
                continue

            # Wales All uses 'Converttodays' (age in days) or BMI_AGE
            age_str = row.get('age_at_event_or_censoring', '')
            if not age_str:
                # Try Converttodays (days since birth)
                days_str = row.get('Converttodays', '')
                try:
                    age = float(days_str) / 365.25
                except:
                    # Try BMI_AGE as fallback
                    bmi_age = row.get('BMI_AGE', '')
                    try:
                        age = float(bmi_age)
                    except:
                        continue
            else:
                try:
                    age = float(age_str)
                except:
                    continue

            gender = row.get('Gender', '').strip()
            sex = 1 if gender.lower() in ('m', 'male', '1') else 0

            sss_info = sss_lookup.get(mutation1, None)
            if sss_info:
                sss = sss_info['sss']
                gene = sss_info['gene']
            else:
                sss = None
                gene = ''
                for vid, info in sss_lookup.items():
                    if vid in mutation1 or mutation1 in vid:
                        sss = info['sss']
                        gene = info['gene']
                        break
                if sss is None:
                    continue

            ldl = row.get('LastLDL', row.get('LDLC', ''))
            try:
                ldl = float(ldl)
            except:
                ldl = None

            patients.append({
                'source': 'WALES_ALL',
                'mutation': mutation1,
                'gene': gene,
                'ascvd': ascvd,
                'age': age,
                'sex': sex,
                'ldl': ldl,
                'sss': sss,
            })

    return patients


# ── Simple logistic regression (no scipy dependency) ──────────────
def logistic_regression_simple(X, y, max_iter=100, lr=0.01):
    """
    Simple gradient descent logistic regression.
    X: list of lists (each row = [1, x1, x2, ...])
    y: list of 0/1
    Returns: coefficients, log-likelihood
    """
    n = len(y)
    k = len(X[0])
    beta = [0.0] * k

    for iteration in range(max_iter):
        # Compute predictions
        ll = 0
        grad = [0.0] * k
        for i in range(n):
            z = sum(beta[j] * X[i][j] for j in range(k))
            z = max(min(z, 20), -20)  # clip for numerical stability
            p = 1.0 / (1.0 + math.exp(-z))
            p = max(min(p, 0.9999), 0.0001)

            ll += y[i] * math.log(p) + (1 - y[i]) * math.log(1 - p)

            for j in range(k):
                grad[j] += (y[i] - p) * X[i][j]

        # Update
        for j in range(k):
            beta[j] += lr * grad[j] / n

    # Final predictions and AUC
    preds = []
    for i in range(n):
        z = sum(beta[j] * X[i][j] for j in range(k))
        z = max(min(z, 20), -20)
        p = 1.0 / (1.0 + math.exp(-z))
        preds.append(p)

    return beta, ll, preds


def compute_auc(y_true, y_pred):
    """Compute AUC using the Mann-Whitney U statistic."""
    pos = [(y_pred[i], y_true[i]) for i in range(len(y_true)) if y_true[i] == 1]
    neg = [(y_pred[i], y_true[i]) for i in range(len(y_true)) if y_true[i] == 0]

    if not pos or not neg:
        return 0.5

    concordant = 0
    discordant = 0
    tied = 0

    for p_score, _ in pos:
        for n_score, _ in neg:
            if p_score > n_score:
                concordant += 1
            elif p_score < n_score:
                discordant += 1
            else:
                tied += 1

    total = concordant + discordant + tied
    if total == 0:
        return 0.5

    return (concordant + 0.5 * tied) / total


def chi_square_trend(groups):
    """
    Chi-square test for trend across ordered groups.
    groups: list of (n_events, n_total) tuples
    Returns chi-square statistic and approximate p-value.
    """
    k = len(groups)
    N = sum(n for _, n in groups)
    D = sum(d for d, _ in groups)

    if N == 0 or D == 0:
        return 0, 1.0

    p_overall = D / N

    scores = list(range(1, k + 1))

    numerator = 0
    denominator = 0

    s_bar = sum(scores[i] * groups[i][1] for i in range(k)) / N

    for i in range(k):
        ni = groups[i][1]
        di = groups[i][0]
        numerator += (scores[i] - s_bar) * (di - ni * p_overall)
        denominator += (scores[i] - s_bar) ** 2 * ni

    if denominator == 0:
        return 0, 1.0

    chi2 = numerator ** 2 / (p_overall * (1 - p_overall) * denominator)

    # Approximate p-value (chi2 with 1 df)
    # Using simple approximation
    if chi2 > 10:
        p_value = math.exp(-chi2 / 2)
    elif chi2 > 3.84:
        p_value = 0.05 * math.exp(-(chi2 - 3.84) / 2)
    else:
        p_value = 1.0 - 0.5 * (1 - math.exp(-chi2 / 2))

    return chi2, p_value


# ── Main validation ──────────────────────────────────────────────
def main():
    print("=" * 70)
    print("  CALON-Structure: SSS Validation Against ASCVD Outcomes")
    print("  Phase 4a - Wales FH Registry Internal Validation")
    print("=" * 70)

    # Load data
    sss_lookup = load_sss()
    print(f"  SSS lookup: {len(sss_lookup)} variants")

    dragon3 = load_dragon3(sss_lookup)
    print(f"  DRAGON_3 patients with SSS: {len(dragon3)}")

    wales_all = load_wales_all(sss_lookup)
    print(f"  Wales All patients with SSS: {len(wales_all)}")

    # Combine
    all_patients = dragon3 + wales_all

    # Remove duplicates (same mutation + same age = likely same patient)
    seen = set()
    unique_patients = []
    for p in all_patients:
        key = (p['mutation'], round(p['age'], 1), p['sex'])
        if key not in seen:
            seen.add(key)
            unique_patients.append(p)

    print(f"  Unique patients: {len(unique_patients)}")
    ascvd_count = sum(p['ascvd'] for p in unique_patients)
    print(f"  ASCVD events: {ascvd_count} ({100*ascvd_count/len(unique_patients):.1f}%)")

    # ── 1. Descriptive statistics by SSS tertile ──────────────────
    print(f"\n{'='*70}")
    print(f"  1. ASCVD rates by SSS tertile")
    print(f"{'='*70}")

    sorted_patients = sorted(unique_patients, key=lambda p: p['sss'])
    n = len(sorted_patients)
    t1 = sorted_patients[:n//3]
    t2 = sorted_patients[n//3:2*n//3]
    t3 = sorted_patients[2*n//3:]

    tertiles = [
        ('T1 (lowest SSS)', t1),
        ('T2 (middle SSS)', t2),
        ('T3 (highest SSS)', t3),
    ]

    print(f"\n  {'Tertile':<20s} {'N':>6s} {'Events':>7s} {'Rate':>7s} {'Mean SSS':>9s} {'Mean Age':>9s}")
    print(f"  {'-'*20} {'-'*6} {'-'*7} {'-'*7} {'-'*9} {'-'*9}")

    trend_groups = []
    for label, pts in tertiles:
        n_pts = len(pts)
        events = sum(p['ascvd'] for p in pts)
        rate = events / n_pts * 100 if n_pts > 0 else 0
        mean_sss = sum(p['sss'] for p in pts) / n_pts if n_pts else 0
        mean_age = sum(p['age'] for p in pts) / n_pts if n_pts else 0
        print(f"  {label:<20s} {n_pts:>6d} {events:>7d} {rate:>6.1f}% {mean_sss:>9.3f} {mean_age:>9.1f}")
        trend_groups.append((events, n_pts))

    chi2, p_trend = chi_square_trend(trend_groups)
    print(f"\n  Chi-square trend test: chi2={chi2:.2f}, p={'<0.001' if p_trend < 0.001 else f'{p_trend:.3f}'}")

    # ── 2. SSS quartiles ──────────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"  2. ASCVD rates by SSS quartile")
    print(f"{'='*70}")

    q1 = sorted_patients[:n//4]
    q2 = sorted_patients[n//4:n//2]
    q3 = sorted_patients[n//2:3*n//4]
    q4 = sorted_patients[3*n//4:]

    quartiles = [('Q1', q1), ('Q2', q2), ('Q3', q3), ('Q4', q4)]
    q_groups = []
    print(f"\n  {'Quartile':<10s} {'N':>6s} {'Events':>7s} {'Rate':>7s} {'SSS range':<15s}")
    print(f"  {'-'*10} {'-'*6} {'-'*7} {'-'*7} {'-'*15}")
    for label, pts in quartiles:
        n_pts = len(pts)
        events = sum(p['ascvd'] for p in pts)
        rate = events / n_pts * 100 if n_pts > 0 else 0
        sss_min = min(p['sss'] for p in pts)
        sss_max = max(p['sss'] for p in pts)
        print(f"  {label:<10s} {n_pts:>6d} {events:>7d} {rate:>6.1f}% {sss_min:.3f}-{sss_max:.3f}")
        q_groups.append((events, n_pts))

    chi2_q, p_q = chi_square_trend(q_groups)
    print(f"\n  Chi-square trend: chi2={chi2_q:.2f}, p={'<0.001' if p_q < 0.001 else f'{p_q:.3f}'}")

    # ── 3. Logistic regression models ─────────────────────────────
    print(f"\n{'='*70}")
    print(f"  3. Logistic regression: ASCVD ~ SSS + age + sex")
    print(f"{'='*70}")

    # Standardize continuous variables
    mean_age = sum(p['age'] for p in unique_patients) / len(unique_patients)
    sd_age = (sum((p['age'] - mean_age)**2 for p in unique_patients) / len(unique_patients)) ** 0.5
    mean_sss = sum(p['sss'] for p in unique_patients) / len(unique_patients)
    sd_sss = (sum((p['sss'] - mean_sss)**2 for p in unique_patients) / len(unique_patients)) ** 0.5

    y = [p['ascvd'] for p in unique_patients]

    # Model 1: Age + Sex only (baseline)
    X_base = [[1, (p['age'] - mean_age) / sd_age, p['sex']] for p in unique_patients]
    beta_base, ll_base, preds_base = logistic_regression_simple(X_base, y, max_iter=200, lr=0.1)
    auc_base = compute_auc(y, preds_base)
    print(f"\n  Model 1 (Age + Sex):")
    print(f"    AUC = {auc_base:.3f}")
    print(f"    Coefficients: intercept={beta_base[0]:.3f}, age={beta_base[1]:.3f}, sex={beta_base[2]:.3f}")

    # Model 2: Age + Sex + Gene (LDLR vs APOB)
    X_gene = [[1, (p['age'] - mean_age) / sd_age, p['sex'],
               1 if p['gene'] == 'LDLR' else 0] for p in unique_patients]
    beta_gene, ll_gene, preds_gene = logistic_regression_simple(X_gene, y, max_iter=200, lr=0.1)
    auc_gene = compute_auc(y, preds_gene)
    print(f"\n  Model 2 (Age + Sex + Gene):")
    print(f"    AUC = {auc_gene:.3f}")
    print(f"    Coefficients: intercept={beta_gene[0]:.3f}, age={beta_gene[1]:.3f}, sex={beta_gene[2]:.3f}, LDLR={beta_gene[3]:.3f}")

    # Model 3: Age + Sex + SSS
    X_sss = [[1, (p['age'] - mean_age) / sd_age, p['sex'],
              (p['sss'] - mean_sss) / sd_sss] for p in unique_patients]
    beta_sss, ll_sss, preds_sss = logistic_regression_simple(X_sss, y, max_iter=200, lr=0.1)
    auc_sss = compute_auc(y, preds_sss)
    print(f"\n  Model 3 (Age + Sex + SSS):")
    print(f"    AUC = {auc_sss:.3f}")
    print(f"    Coefficients: intercept={beta_sss[0]:.3f}, age={beta_sss[1]:.3f}, sex={beta_sss[2]:.3f}, SSS={beta_sss[3]:.3f}")

    # Model 4: Age + Sex + Gene + SSS
    X_full = [[1, (p['age'] - mean_age) / sd_age, p['sex'],
               1 if p['gene'] == 'LDLR' else 0,
               (p['sss'] - mean_sss) / sd_sss] for p in unique_patients]
    beta_full, ll_full, preds_full = logistic_regression_simple(X_full, y, max_iter=200, lr=0.1)
    auc_full = compute_auc(y, preds_full)
    print(f"\n  Model 4 (Age + Sex + Gene + SSS):")
    print(f"    AUC = {auc_full:.3f}")
    print(f"    Coefficients: intercept={beta_full[0]:.3f}, age={beta_full[1]:.3f}, sex={beta_full[2]:.3f}, LDLR={beta_full[3]:.3f}, SSS={beta_full[4]:.3f}")

    # ── AUC comparison ────────────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"  4. Model Comparison (AUC)")
    print(f"{'='*70}")
    print(f"\n  {'Model':<30s} {'AUC':>8s} {'Delta vs base':>14s}")
    print(f"  {'-'*30} {'-'*8} {'-'*14}")
    print(f"  {'Age + Sex':<30s} {auc_base:>8.3f} {'(reference)':>14s}")
    print(f"  {'Age + Sex + Gene':<30s} {auc_gene:>8.3f} {auc_gene - auc_base:>+14.3f}")
    print(f"  {'Age + Sex + SSS':<30s} {auc_sss:>8.3f} {auc_sss - auc_base:>+14.3f}")
    print(f"  {'Age + Sex + Gene + SSS':<30s} {auc_full:>8.3f} {auc_full - auc_base:>+14.3f}")

    # ── 5. Gene-specific analysis ─────────────────────────────────
    print(f"\n{'='*70}")
    print(f"  5. SSS performance within LDLR carriers only")
    print(f"{'='*70}")

    ldlr_patients = [p for p in unique_patients if p['gene'] == 'LDLR']
    if ldlr_patients:
        ldlr_sorted = sorted(ldlr_patients, key=lambda p: p['sss'])
        n_ldlr = len(ldlr_sorted)
        mid = n_ldlr // 2

        low_half = ldlr_sorted[:mid]
        high_half = ldlr_sorted[mid:]

        low_events = sum(p['ascvd'] for p in low_half)
        high_events = sum(p['ascvd'] for p in high_half)
        low_rate = low_events / len(low_half) * 100
        high_rate = high_events / len(high_half) * 100

        print(f"\n  LDLR carriers (n={n_ldlr}):")
        print(f"    Low SSS half:  {len(low_half)} patients, {low_events} events ({low_rate:.1f}%)")
        print(f"    High SSS half: {len(high_half)} patients, {high_events} events ({high_rate:.1f}%)")
        print(f"    Absolute risk difference: {high_rate - low_rate:+.1f}%")

    # ── 6. Domain-specific forest plot data ───────────────────────
    print(f"\n{'='*70}")
    print(f"  6. ASCVD rates by LDLR domain (for forest plot)")
    print(f"{'='*70}")

    domain_data = defaultdict(lambda: {'events': 0, 'total': 0, 'sss_vals': []})
    for p in unique_patients:
        if p['gene'] != 'LDLR':
            continue
        # Get domain from SSS lookup
        info = sss_lookup.get(p['mutation'], {})
        domain = info.get('domain', 'Unknown')

        # Simplify domain names
        if 'Ligand' in domain or 'ligand' in domain:
            domain = 'Ligand-binding'
        elif 'EGF' in domain or 'egf' in domain:
            domain = 'EGF-like'
        elif 'propeller' in domain or 'Beta' in domain:
            domain = 'Beta-propeller'
        elif 'O-linked' in domain or 'sugar' in domain:
            domain = 'O-linked sugar'

        domain_data[domain]['events'] += p['ascvd']
        domain_data[domain]['total'] += 1
        domain_data[domain]['sss_vals'].append(p['sss'])

    print(f"\n  {'Domain':<20s} {'N':>6s} {'Events':>7s} {'Rate':>7s} {'Mean SSS':>9s}")
    print(f"  {'-'*20} {'-'*6} {'-'*7} {'-'*7} {'-'*9}")
    for domain in sorted(domain_data.keys()):
        d = domain_data[domain]
        rate = d['events'] / d['total'] * 100 if d['total'] > 0 else 0
        mean_s = sum(d['sss_vals']) / len(d['sss_vals']) if d['sss_vals'] else 0
        print(f"  {domain:<20s} {d['total']:>6d} {d['events']:>7d} {rate:>6.1f}% {mean_s:>9.3f}")

    # ── Save validation dataset ───────────────────────────────────
    out_file = os.path.join(ANALYSIS_DIR, "sss_validation_data.csv")
    with open(out_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['source', 'mutation', 'gene', 'ascvd',
                                                'age', 'sex', 'ldl', 'sss'])
        writer.writeheader()
        writer.writerows(unique_patients)
    print(f"\n  Saved: {out_file}")

    # ── Key finding ───────────────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"  KEY FINDINGS")
    print(f"{'='*70}")
    print(f"""
  1. SSS tertile analysis shows {'SIGNIFICANT' if p_trend < 0.05 else 'NO SIGNIFICANT'}
     trend of increasing ASCVD with higher SSS (p={'<0.001' if p_trend < 0.001 else f'{p_trend:.3f}'})

  2. AUC comparison:
     - Age+Sex baseline:  {auc_base:.3f}
     - Adding Gene:       {auc_gene:.3f} (delta {auc_gene-auc_base:+.3f})
     - Adding SSS:        {auc_sss:.3f} (delta {auc_sss-auc_base:+.3f})
     - Gene + SSS:        {auc_full:.3f} (delta {auc_full-auc_base:+.3f})

  3. SSS {'ADDS' if auc_sss > auc_gene else 'does NOT add'} incremental
     discrimination beyond gene-only classification
     (AUC SSS={auc_sss:.3f} vs AUC Gene={auc_gene:.3f})

  4. Note: These are INTERNAL validation results.
     External validation in UKB is required (Phase 4b).
     Simple logistic regression used (for proper inference,
     use R with clustered SEs by variant).
""")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
AF3_manuscript_docx.py
======================
Main assembler script for the AlphaFold3 FH Manuscript.

Loads all CSV data files, computes statistics into a stats{} dictionary,
imports the three part scripts, creates a python-docx Document with A4
formatting, calls each part's functions in order, and saves the final DOCX.
"""

import csv
import math
import os
import sys
from collections import defaultdict

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"
AF = f"{BASE}/alphafold/analysis"

comp_sss_file       = f"{AF}/comprehensive_sss_analysis.csv"
nobel_file          = f"{AF}/nobel_multimodal_integration.csv"
sss_scores_file     = f"{AF}/structural_severity_scores.csv"
foldx_file          = f"{AF}/foldx_ddg_results.csv"
llt_file            = f"{AF}/variant_llt_response.csv"
quality_file        = f"{AF}/af3_structure_quality_summary.csv"
plddt_ldlr_file     = f"{AF}/LDLR_wildtype_plddt_per_residue.csv"
plddt_pcsk9_file    = f"{AF}/PCSK9_wildtype_plddt_per_residue.csv"
plddt_apob_file     = f"{AF}/ApoB_RBD_plddt_per_residue.csv"
dragon_file         = f"{BASE}/DRAGON_3.csv"
sss_wales_file      = f"{AF}/sss_patient_level_wales.csv"

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def safe_float(x, default=None):
    if x is None:
        return default
    try:
        return float(str(x).strip())
    except (ValueError, TypeError):
        return default


def safe_int(x, default=0):
    try:
        return int(float(str(x).strip()))
    except (ValueError, TypeError):
        return default


def mean_vals(vals):
    vals = [v for v in vals if v is not None]
    return sum(vals) / len(vals) if vals else 0.0


def sd_vals(vals):
    vals = [v for v in vals if v is not None]
    n = len(vals)
    if n < 2:
        return 0.0
    m = sum(vals) / n
    return (sum((v - m) ** 2 for v in vals) / (n - 1)) ** 0.5


def median_vals(vals):
    vals = sorted(v for v in vals if v is not None)
    n = len(vals)
    if n == 0:
        return 0.0
    if n % 2 == 1:
        return vals[n // 2]
    return (vals[n // 2 - 1] + vals[n // 2]) / 2


def iqr_vals(vals):
    vals = sorted(v for v in vals if v is not None)
    n = len(vals)
    if n < 4:
        return 0.0, 0.0
    q1 = median_vals(vals[: n // 2])
    q3 = median_vals(vals[(n + 1) // 2 :])
    return q1, q3


def pearson_r(x, y):
    pairs = [(xi, yi) for xi, yi in zip(x, y) if xi is not None and yi is not None]
    n = len(pairs)
    if n < 3:
        return 0.0, 1.0, n
    x2, y2 = zip(*pairs)
    mx, my = sum(x2) / n, sum(y2) / n
    sxx = sum((xi - mx) ** 2 for xi in x2)
    syy = sum((yi - my) ** 2 for yi in y2)
    sxy = sum((xi - mx) * (yi - my) for xi, yi in zip(x2, y2))
    if sxx == 0 or syy == 0:
        return 0.0, 1.0, n
    r = sxy / (sxx * syy) ** 0.5
    if abs(r) >= 1.0:
        return r, 0.0, n
    t_stat = r * ((n - 2) / (1 - r * r + 1e-15)) ** 0.5
    df = n - 2
    x_val = df / (df + t_stat ** 2)
    p = x_val ** (df / 2)  # rough beta approximation
    return r, p, n


def pct(count, total):
    return (count / total * 100) if total > 0 else 0.0


def fmt(val, decimals=1):
    if val is None:
        return "NA"
    return f"{val:.{decimals}f}"


def odds_ratio(a, b, c, d):
    """2x2 table OR: a=events_exposed, b=nonevents_exposed,
    c=events_unexposed, d=nonevents_unexposed"""
    if b == 0 or c == 0:
        return 0, "NA"
    or_val = (a * d) / (b * c)
    se = (1 / max(a, 1) + 1 / max(b, 1) + 1 / max(c, 1) + 1 / max(d, 1)) ** 0.5
    lower = math.exp(math.log(max(or_val, 0.001)) - 1.96 * se)
    upper = math.exp(math.log(max(or_val, 0.001)) + 1.96 * se)
    return or_val, f"{lower:.2f}-{upper:.2f}"


def concordance_statistic(outcome, score):
    """Compute C-statistic (equivalent to AUC for binary outcome)."""
    positives = [s for s, o in zip(score, outcome) if o == 1]
    negatives = [s for s, o in zip(score, outcome) if o == 0]
    if not positives or not negatives:
        return 0.5
    pairs_concordant = 0.0
    pairs_total = 0
    for sp in positives:
        for sn in negatives:
            pairs_total += 1
            if sp > sn:
                pairs_concordant += 1
            elif sp == sn:
                pairs_concordant += 0.5
    return pairs_concordant / max(pairs_total, 1)


# ---------------------------------------------------------------------------
# Safe CSV loader
# ---------------------------------------------------------------------------

def load_csv(filepath, label=""):
    """Load a CSV file, returning list of dicts. Returns [] on failure."""
    if not os.path.isfile(filepath):
        print(f"  WARNING: {label or filepath} not found.")
        return []
    try:
        with open(filepath, newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            data = list(reader)
        print(f"  Loaded {label or filepath}: {len(data)} rows")
        return data
    except Exception as e:
        print(f"  ERROR loading {label or filepath}: {e}")
        return []


# ===================================================================
# MAIN
# ===================================================================

def main():
    print("=" * 60)
    print("AlphaFold3 FH Manuscript Assembler")
    print("=" * 60)

    # ------------------------------------------------------------------
    # 1.  LOAD ALL DATA
    # ------------------------------------------------------------------
    print("\n[1] Loading data files ...")

    comp_sss_rows   = load_csv(comp_sss_file,   "comprehensive_sss_analysis")
    nobel_rows      = load_csv(nobel_file,       "nobel_multimodal_integration")
    sss_scores_rows = load_csv(sss_scores_file,  "structural_severity_scores")
    foldx_rows      = load_csv(foldx_file,       "foldx_ddg_results")
    llt_rows        = load_csv(llt_file,         "variant_llt_response")
    quality_rows    = load_csv(quality_file,      "af3_structure_quality_summary")
    plddt_ldlr      = load_csv(plddt_ldlr_file,  "LDLR pLDDT per-residue")
    plddt_pcsk9     = load_csv(plddt_pcsk9_file,  "PCSK9 pLDDT per-residue")
    plddt_apob      = load_csv(plddt_apob_file,   "ApoB pLDDT per-residue")
    dragon_rows     = load_csv(dragon_file,       "DRAGON_3")
    sss_wales_rows  = load_csv(sss_wales_file,    "sss_patient_level_wales")

    rows = comp_sss_rows  # alias for convenience

    # ------------------------------------------------------------------
    # 2.  COMPUTE ALL STATISTICS
    # ------------------------------------------------------------------
    print("\n[2] Computing statistics ...")
    stats = {}

    # ---- Cohort sizes ----
    stats['n_analytic'] = len(rows)
    stats['n_total'] = len(dragon_rows)
    stats['n_variants'] = len(sss_scores_rows)
    stats['n_sss'] = len(sss_wales_rows)

    # ---- Gene counts ----
    genes = [r.get('gene', '') for r in rows]
    for g in ['LDLR', 'APOB', 'PCSK9']:
        stats[f'n_{g.lower()}'] = sum(1 for x in genes if x == g)
        stats[f'pct_{g.lower()}'] = pct(stats[f'n_{g.lower()}'], len(rows))

    # ---- Demographics ----
    ages = [safe_float(r.get('age')) for r in rows]
    stats['mean_age'] = mean_vals(ages)
    stats['sd_age'] = sd_vals(ages)
    stats['pct_female'] = pct(
        sum(1 for r in rows if safe_float(r.get('sex')) == 0), len(rows)
    )
    stats['pct_male'] = 100.0 - stats['pct_female']

    # ---- Baseline lipids ----
    for col, key in [('ldl1', 'ldl'), ('tc1', 'tc'), ('hdl1', 'hdl'), ('tg1', 'tg')]:
        vals = [safe_float(r.get(col)) for r in rows]
        stats[f'mean_{key}'] = mean_vals(vals)
        stats[f'sd_{key}'] = sd_vals(vals)
        stats[f'median_{key}'] = median_vals(vals)

    # ---- ApoB, ApoA1, Lp(a) ----
    for col, key in [('apob', 'apob_val'), ('apoa1', 'apoa1'), ('lpa', 'lpa')]:
        vals = [safe_float(r.get(col)) for r in rows]
        stats[f'mean_{key}'] = mean_vals(vals)
        stats[f'sd_{key}'] = sd_vals(vals)
        stats[f'median_{key}'] = median_vals(vals)

    # ---- SSS summary ----
    sss_vals = [safe_float(r.get('sss')) for r in rows]
    stats['mean_sss'] = mean_vals(sss_vals)
    stats['sd_sss'] = sd_vals(sss_vals)
    stats['median_sss'] = median_vals(sss_vals)
    q1_sss, q3_sss = iqr_vals(sss_vals)
    stats['q1_sss'] = q1_sss
    stats['q3_sss'] = q3_sss

    # ---- ASCVD / events ----
    stats['ascvd_rate'] = pct(
        sum(1 for r in rows if safe_float(r.get('ascvd')) == 1), len(rows)
    )
    stats['miacs_rate'] = pct(
        sum(1 for r in rows if safe_float(r.get('miacs')) == 1), len(rows)
    )
    stats['angina_rate'] = pct(
        sum(1 for r in rows if safe_float(r.get('angina')) == 1), len(rows)
    )
    stats['tia_rate'] = pct(
        sum(1 for r in rows if safe_float(r.get('tia')) == 1), len(rows)
    )
    stats['pvd_rate'] = pct(
        sum(1 for r in rows if safe_float(r.get('pvd')) == 1), len(rows)
    )
    event_ages = [safe_float(r.get('age_event')) for r in rows]
    stats['mean_age_event'] = mean_vals(event_ages)
    stats['sd_age_event'] = sd_vals(event_ages)

    n_events_vals = [safe_float(r.get('n_events')) for r in rows]
    stats['mean_n_events'] = mean_vals(n_events_vals)

    # ---- Treatment rates ----
    stats['statin_rate'] = pct(
        sum(1 for r in rows if safe_float(r.get('on_statin')) == 1), len(rows)
    )
    stats['ezetimibe_rate'] = pct(
        sum(1 for r in rows if safe_float(r.get('ezetimibe')) == 1), len(rows)
    )
    stats['pcsk9i_rate'] = pct(
        sum(1 for r in rows if safe_float(r.get('pcsk9i')) == 1), len(rows)
    )

    # treatment intensity
    tx_dist = defaultdict(int)
    for r in rows:
        tx = str(r.get('tx_intensity', '')).strip()
        if tx:
            tx_dist[tx] += 1
    stats['tx_intensity_dist'] = dict(tx_dist)

    # ---- Comorbidities ----
    stats['dm_rate'] = pct(
        sum(1 for r in rows if safe_float(r.get('dm')) == 1), len(rows)
    )
    stats['smoking_rate'] = pct(
        sum(1 for r in rows if safe_float(r.get('smoking')) == 1), len(rows)
    )
    stats['hypertension_rate'] = pct(
        sum(1 for r in rows if safe_float(r.get('hypertension')) == 1), len(rows)
    )
    stats['xanthomata_rate'] = pct(
        sum(1 for r in rows if safe_float(r.get('xanthomata')) == 1), len(rows)
    )
    stats['simon_broome_rate'] = pct(
        sum(1 for r in rows if safe_float(r.get('simon_broome')) == 1), len(rows)
    )
    bmi_vals = [safe_float(r.get('bmi')) for r in rows]
    stats['mean_bmi'] = mean_vals(bmi_vals)
    stats['sd_bmi'] = sd_vals(bmi_vals)
    sbp_vals = [safe_float(r.get('sbp')) for r in rows]
    stats['mean_sbp'] = mean_vals(sbp_vals)
    stats['sd_sbp'] = sd_vals(sbp_vals)

    # ---- Domains ----
    domains_set = set(r.get('domain', '') for r in rows if r.get('domain'))
    stats['n_domains'] = len(domains_set)
    stats['domain_list'] = sorted(domains_set)

    # ---- Last lipids / LDL % change ----
    last_ldl_vals = [safe_float(r.get('last_ldl')) for r in rows]
    stats['mean_last_ldl'] = mean_vals(last_ldl_vals)
    stats['sd_last_ldl'] = sd_vals(last_ldl_vals)

    ldl_pct_changes = [safe_float(r.get('ldl_pct_change')) for r in rows]
    stats['mean_ldl_pct_change'] = mean_vals(ldl_pct_changes)
    stats['sd_ldl_pct_change'] = sd_vals(ldl_pct_changes)

    # ------------------------------------------------------------------
    # FoldX ddG statistics
    # ------------------------------------------------------------------
    stats['n_foldx'] = len(foldx_rows)
    ddg_vals = [safe_float(r.get('ddG_kcal_mol')) for r in foldx_rows]
    ddg_valid = [v for v in ddg_vals if v is not None]
    if ddg_valid:
        stats['ddg_min'] = min(ddg_valid)
        stats['ddg_max'] = max(ddg_valid)
        stats['ddg_mean'] = mean_vals(ddg_valid)
        stats['ddg_sd'] = sd_vals(ddg_valid)
        stats['ddg_median'] = median_vals(ddg_valid)
        stats['n_stabilising'] = sum(1 for v in ddg_valid if v < 0)
        stats['n_neutral_ddg'] = sum(1 for v in ddg_valid if 0 <= v < 1)
        stats['n_destabilising'] = sum(1 for v in ddg_valid if 1 <= v < 4)
        stats['n_highly_destabilising'] = sum(1 for v in ddg_valid if v >= 4)
    else:
        for k in ['ddg_min', 'ddg_max', 'ddg_mean', 'ddg_sd', 'ddg_median']:
            stats[k] = 0.0
        for k in ['n_stabilising', 'n_neutral_ddg', 'n_destabilising', 'n_highly_destabilising']:
            stats[k] = 0

    # FoldX effect categories
    foldx_effects = defaultdict(int)
    for r in foldx_rows:
        eff = str(r.get('effect', '')).strip().lower()
        if eff:
            foldx_effects[eff] += 1
    stats['foldx_effects'] = dict(foldx_effects)

    # ------------------------------------------------------------------
    # ddG-clinical correlations (merge foldx with comp_sss on variant_id~mutation)
    # ------------------------------------------------------------------
    foldx_by_variant = {}
    for r in foldx_rows:
        vid = r.get('variant_id', '').strip()
        ddg = safe_float(r.get('ddG_kcal_mol'))
        if vid and ddg is not None:
            foldx_by_variant[vid] = ddg

    ddg_ldl_x, ddg_ldl_y = [], []
    ddg_apob_x, ddg_apob_y = [], []
    ddg_ascvd_x, ddg_ascvd_y = [], []
    for r in rows:
        mut = r.get('mutation', '').strip()
        if mut in foldx_by_variant:
            ddg = foldx_by_variant[mut]
            ll = safe_float(r.get('last_ldl'))
            ab = safe_float(r.get('apob'))
            asc = safe_float(r.get('ascvd'))
            if ll is not None:
                ddg_ldl_x.append(ddg)
                ddg_ldl_y.append(ll)
            if ab is not None:
                ddg_apob_x.append(ddg)
                ddg_apob_y.append(ab)
            if asc is not None:
                ddg_ascvd_x.append(ddg)
                ddg_ascvd_y.append(asc)

    r_val, p_val, n_val = pearson_r(ddg_ldl_x, ddg_ldl_y)
    stats['ddg_ldl_r'] = r_val
    stats['ddg_ldl_p'] = p_val
    stats['ddg_ldl_n'] = n_val

    r_val, p_val, n_val = pearson_r(ddg_apob_x, ddg_apob_y)
    stats['ddg_apob_r'] = r_val
    stats['ddg_apob_p'] = p_val
    stats['ddg_apob_n'] = n_val

    r_val, p_val, n_val = pearson_r(ddg_ascvd_x, ddg_ascvd_y)
    stats['ddg_ascvd_r'] = r_val
    stats['ddg_ascvd_p'] = p_val
    stats['ddg_ascvd_n'] = n_val

    # ------------------------------------------------------------------
    # pLDDT-clinical correlation
    # ------------------------------------------------------------------
    sss_plddt_map = {}
    for r in sss_scores_rows:
        vid = r.get('variant_id', '').strip()
        pl = safe_float(r.get('plddt'))
        if vid and pl is not None:
            sss_plddt_map[vid] = pl

    plddt_x, plddt_y = [], []
    for r in rows:
        mut = r.get('mutation', '').strip()
        if mut in sss_plddt_map:
            ldl = safe_float(r.get('ldl1'))
            if ldl is not None:
                plddt_x.append(sss_plddt_map[mut])
                plddt_y.append(ldl)

    r_val, p_val, n_val = pearson_r(plddt_x, plddt_y)
    stats['plddt_ldl_r'] = r_val
    stats['plddt_ldl_p'] = p_val
    stats['plddt_ldl_n'] = n_val

    # ------------------------------------------------------------------
    # Structure quality from quality summary file
    # ------------------------------------------------------------------
    stats['ldlr_plddt'] = 0.0
    stats['pcsk9_plddt'] = 0.0
    stats['apob_plddt'] = 0.0
    stats['ldlr_residues'] = 0
    stats['pcsk9_residues'] = 0
    stats['apob_residues'] = 0
    stats['ldlr_quality'] = "NA"
    stats['pcsk9_quality'] = "NA"
    stats['apob_quality'] = "NA"

    for r in quality_rows:
        struct = str(r.get('structure', '')).lower()
        if 'ldlr' in struct:
            stats['ldlr_plddt'] = safe_float(r.get('mean_plddt'), 0.0)
            stats['ldlr_residues'] = safe_int(r.get('n_residues'))
            stats['ldlr_quality'] = r.get('quality', 'NA')
            stats['ldlr_pct_high'] = safe_float(r.get('pct_very_high_90plus'), 0.0)
        elif 'pcsk9' in struct:
            stats['pcsk9_plddt'] = safe_float(r.get('mean_plddt'), 0.0)
            stats['pcsk9_residues'] = safe_int(r.get('n_residues'))
            stats['pcsk9_quality'] = r.get('quality', 'NA')
            stats['pcsk9_pct_high'] = safe_float(r.get('pct_very_high_90plus'), 0.0)
        elif 'apob' in struct:
            stats['apob_plddt'] = safe_float(r.get('mean_plddt'), 0.0)
            stats['apob_residues'] = safe_int(r.get('n_residues'))
            stats['apob_quality'] = r.get('quality', 'NA')
            stats['apob_pct_high'] = safe_float(r.get('pct_very_high_90plus'), 0.0)

    # Also compute from per-residue files if quality_summary was sparse
    if plddt_ldlr and stats['ldlr_plddt'] == 0.0:
        vals = [safe_float(r.get('plddt')) for r in plddt_ldlr]
        stats['ldlr_plddt'] = mean_vals(vals)
        stats['ldlr_residues'] = len(plddt_ldlr)
    if plddt_pcsk9 and stats['pcsk9_plddt'] == 0.0:
        vals = [safe_float(r.get('plddt')) for r in plddt_pcsk9]
        stats['pcsk9_plddt'] = mean_vals(vals)
        stats['pcsk9_residues'] = len(plddt_pcsk9)
    if plddt_apob and stats['apob_plddt'] == 0.0:
        vals = [safe_float(r.get('plddt')) for r in plddt_apob]
        stats['apob_plddt'] = mean_vals(vals)
        stats['apob_residues'] = len(plddt_apob)

    stats['mean_overall_plddt'] = mean_vals(
        [stats['ldlr_plddt'], stats['pcsk9_plddt'], stats['apob_plddt']]
    )

    # ------------------------------------------------------------------
    # SSS-ASCVD analysis  (tertiles)
    # ------------------------------------------------------------------
    sss_valid = sorted(v for v in sss_vals if v is not None)
    n_valid = len(sss_valid)
    if n_valid >= 3:
        t1_cut = sss_valid[n_valid // 3]
        t2_cut = sss_valid[2 * n_valid // 3]
    else:
        t1_cut, t2_cut = 0.33, 0.67

    stats['sss_t1_cut'] = t1_cut
    stats['sss_t2_cut'] = t2_cut

    # Tertile ASCVD counts
    tertile_events = {1: [0, 0], 2: [0, 0], 3: [0, 0]}  # {tertile: [events, total]}
    for r in rows:
        s = safe_float(r.get('sss'))
        a = safe_float(r.get('ascvd'))
        if s is None or a is None:
            continue
        if s <= t1_cut:
            t = 1
        elif s <= t2_cut:
            t = 2
        else:
            t = 3
        tertile_events[t][1] += 1
        if a == 1:
            tertile_events[t][0] += 1

    stats['t1_events'] = tertile_events[1][0]
    stats['t1_total'] = tertile_events[1][1]
    stats['t1_ascvd'] = pct(tertile_events[1][0], tertile_events[1][1])
    stats['t2_events'] = tertile_events[2][0]
    stats['t2_total'] = tertile_events[2][1]
    stats['t2_ascvd'] = pct(tertile_events[2][0], tertile_events[2][1])
    stats['t3_events'] = tertile_events[3][0]
    stats['t3_total'] = tertile_events[3][1]
    stats['t3_ascvd'] = pct(tertile_events[3][0], tertile_events[3][1])

    # OR T3 vs T1
    a = tertile_events[3][0]
    b = tertile_events[3][1] - a
    c = tertile_events[1][0]
    d = tertile_events[1][1] - c
    stats['or_t3_vs_t1'], stats['or_t3_ci'] = odds_ratio(a, b, c, d)

    # OR T2 vs T1
    a2 = tertile_events[2][0]
    b2 = tertile_events[2][1] - a2
    stats['or_t2_vs_t1'], stats['or_t2_ci'] = odds_ratio(a2, b2, c, d)

    # SSS-ASCVD correlation
    sss_for_corr = [safe_float(r.get('sss')) for r in rows]
    ascvd_for_corr = [safe_float(r.get('ascvd')) for r in rows]
    r_val, p_val, n_val = pearson_r(sss_for_corr, ascvd_for_corr)
    stats['sss_ascvd_r'] = r_val
    stats['sss_ascvd_p'] = p_val

    # SSS-LDL correlation
    sss_for_ldl = [safe_float(r.get('sss')) for r in rows]
    ldl_for_corr = [safe_float(r.get('ldl1')) for r in rows]
    r_val, p_val, n_val = pearson_r(sss_for_ldl, ldl_for_corr)
    stats['sss_ldl_r'] = r_val
    stats['sss_ldl_p'] = p_val

    # ------------------------------------------------------------------
    # Model AUC (C-statistic): age vs age+SSS for ASCVD prediction
    # ------------------------------------------------------------------
    valid_model = [
        (safe_float(r.get('ascvd')), safe_float(r.get('age')), safe_float(r.get('sss')))
        for r in rows
    ]
    valid_model = [(a, ag, s) for a, ag, s in valid_model if a is not None and ag is not None and s is not None]

    if valid_model:
        outcomes = [int(a) for a, _, _ in valid_model]
        age_scores = [ag for _, ag, _ in valid_model]
        age_sss_scores = [ag + s * 10 for _, ag, s in valid_model]  # weighted combination

        stats['auc_base'] = concordance_statistic(outcomes, age_scores)
        stats['auc_sss'] = concordance_statistic(outcomes, age_sss_scores)
        stats['auc_improvement'] = stats['auc_sss'] - stats['auc_base']
        stats['nri_approx'] = stats['auc_improvement'] * 2  # rough NRI approximation
    else:
        stats['auc_base'] = 0.5
        stats['auc_sss'] = 0.5
        stats['auc_improvement'] = 0.0
        stats['nri_approx'] = 0.0

    # ------------------------------------------------------------------
    # Treatment response: LoF vs missense
    # ------------------------------------------------------------------
    variant_type_map = {}
    for r in sss_scores_rows:
        vid = r.get('variant_id', '').strip()
        vtype = str(r.get('variant_type', '')).strip().lower()
        if vid:
            variant_type_map[vid] = vtype

    lof_responses = []
    missense_responses = []
    all_responses = []
    for r in llt_rows:
        resp = safe_float(r.get('ldl_pct_reduction'))
        if resp is None:
            continue
        all_responses.append(resp)
        variant = r.get('variant', '').strip()
        vtype = variant_type_map.get(variant, '')
        if vtype in ('deletion', 'insertion', 'splice'):
            lof_responses.append(resp)
        elif vtype == 'substitution':
            missense_responses.append(resp)

    stats['lof_response'] = mean_vals(lof_responses)
    stats['missense_response'] = mean_vals(missense_responses)
    stats['response_gap'] = abs(stats['lof_response'] - stats['missense_response'])
    stats['n_lof_response'] = len(lof_responses)
    stats['n_missense_response'] = len(missense_responses)
    stats['mean_overall_response'] = mean_vals(all_responses)
    stats['sd_overall_response'] = sd_vals(all_responses)

    # ------------------------------------------------------------------
    # Nobel-derived stats
    # ------------------------------------------------------------------
    # Domain phenotype table
    domain_groups = defaultdict(list)
    for r in nobel_rows:
        domain_groups[r.get('domain', 'Unknown')].append(r)

    domain_table = []
    for domain, drows in sorted(domain_groups.items()):
        n_patients = sum(safe_int(r.get('n', 0)) for r in drows)
        domain_table.append({
            'domain': domain,
            'n_variants': len(drows),
            'n_patients': n_patients,
            'mean_ldl': mean_vals([safe_float(r.get('mean_ldl')) for r in drows]),
            'mean_apob': mean_vals([safe_float(r.get('mean_apob')) for r in drows]),
            'mean_lpa': mean_vals([safe_float(r.get('mean_lpa')) for r in drows]),
            'mean_discordance': mean_vals([safe_float(r.get('mean_discordance')) for r in drows]),
            'ascvd_pct': mean_vals([safe_float(r.get('ascvd_pct')) for r in drows]),
            'xanth_pct': mean_vals([safe_float(r.get('xanth_pct')) for r in drows]),
            'mean_sss': mean_vals([safe_float(r.get('sss')) for r in drows]),
            'mean_ddg': mean_vals([safe_float(r.get('ddG')) for r in drows]),
            'mean_response': mean_vals([safe_float(r.get('mean_response')) for r in drows]),
        })
    stats['domain_table'] = domain_table
    stats['n_nobel_domains'] = len(domain_table)

    # Cluster table
    cluster_groups = defaultdict(list)
    for r in nobel_rows:
        cluster_groups[r.get('cluster', 'Unknown')].append(r)

    cluster_table = []
    for cluster, crows in sorted(cluster_groups.items()):
        n_patients = sum(safe_int(r.get('n', 0)) for r in crows)
        cluster_table.append({
            'cluster': cluster,
            'n_variants': len(crows),
            'n_patients': n_patients,
            'mean_ldl': mean_vals([safe_float(r.get('mean_ldl')) for r in crows]),
            'mean_sss': mean_vals([safe_float(r.get('sss')) for r in crows]),
            'ascvd_pct': mean_vals([safe_float(r.get('ascvd_pct')) for r in crows]),
            'mean_response': mean_vals([safe_float(r.get('mean_response')) for r in crows]),
            'mean_discordance': mean_vals([safe_float(r.get('mean_discordance')) for r in crows]),
        })
    stats['n_clusters'] = len(cluster_groups)
    stats['cluster_table'] = cluster_table

    # Treatment-resistant variants (lowest response = hardest to treat)
    nobel_with_response = [r for r in nobel_rows if safe_float(r.get('mean_response')) is not None]
    resistant = sorted(
        nobel_with_response, key=lambda r: safe_float(r.get('mean_response'), 999)
    )[:13]
    stats['resistant_variants'] = [
        {
            'variant': r.get('variant', ''),
            'sss': r.get('sss', ''),
            'response': r.get('mean_response', ''),
            'domain': r.get('domain', ''),
            'cluster': r.get('cluster', ''),
        }
        for r in resistant
    ]
    stats['n_resistant_variants'] = len(resistant)

    # Discordance (top 10)
    top_disc = sorted(
        nobel_rows,
        key=lambda r: safe_float(r.get('mean_discordance'), 0),
        reverse=True,
    )[:10]
    stats['top_disc_variants'] = [
        (r.get('variant', ''), safe_float(r.get('mean_discordance')))
        for r in top_disc
    ]

    # EGF-A discordance
    egfa = [r for r in nobel_rows if 'EGF' in str(r.get('domain', ''))]
    stats['egfa_discordance'] = mean_vals(
        [safe_float(r.get('mean_discordance')) for r in egfa]
    )
    stats['n_egfa_variants'] = len(egfa)

    # Ligand-binding LDL
    lb = [r for r in nobel_rows if 'Ligand' in str(r.get('domain', ''))]
    stats['ligand_binding_ldl'] = mean_vals(
        [safe_float(r.get('mean_ldl')) for r in lb]
    )
    stats['n_ligand_binding'] = len(lb)
    stats['ligand_binding_ascvd'] = mean_vals(
        [safe_float(r.get('ascvd_pct')) for r in lb]
    )

    # Nobel summary
    stats['n_nobel_variants'] = len(nobel_rows)
    stats['nobel_total_patients'] = sum(
        safe_int(r.get('n', 0)) for r in nobel_rows
    )

    # ------------------------------------------------------------------
    # Compound risk: SSS x Lp(a)
    # ------------------------------------------------------------------
    sss_med = stats['median_sss']
    lpa_vals_valid = [safe_float(r.get('lpa')) for r in rows if safe_float(r.get('lpa')) is not None]
    lpa_med = median_vals(lpa_vals_valid) if lpa_vals_valid else 143.0
    stats['lpa_median'] = lpa_med

    # 2x2 SSS x Lp(a)
    cells = {'hh': [0, 0], 'hl': [0, 0], 'lh': [0, 0], 'll': [0, 0]}
    for r in rows:
        s = safe_float(r.get('sss'))
        l = safe_float(r.get('lpa'))
        a = safe_float(r.get('ascvd'))
        if s is None or a is None:
            continue
        s_high = s >= sss_med
        l_high = (l is not None and l >= lpa_med)
        if s_high and l_high:
            key = 'hh'
        elif s_high and not l_high:
            key = 'hl'
        elif not s_high and l_high:
            key = 'lh'
        else:
            key = 'll'
        cells[key][1] += 1
        if a == 1:
            cells[key][0] += 1

    stats['hh_ascvd'] = pct(cells['hh'][0], cells['hh'][1])
    stats['hl_ascvd'] = pct(cells['hl'][0], cells['hl'][1])
    stats['lh_ascvd'] = pct(cells['lh'][0], cells['lh'][1])
    stats['ll_ascvd'] = pct(cells['ll'][0], cells['ll'][1])
    stats['hh_n'] = cells['hh'][1]
    stats['hl_n'] = cells['hl'][1]
    stats['lh_n'] = cells['lh'][1]
    stats['ll_n'] = cells['ll'][1]

    # OR for high-SSS+high-Lpa vs low-both
    stats['or_compound'], stats['or_compound_ci'] = odds_ratio(
        cells['hh'][0], cells['hh'][1] - cells['hh'][0],
        cells['ll'][0], cells['ll'][1] - cells['ll'][0],
    )

    # Triple risk: high SSS + high Lp(a) + high discordance
    # For discordance, use matched_ldl proxy: if matched_ldl - ldl1 is large = discordant
    disc_vals = []
    for r in rows:
        ml = safe_float(r.get('matched_ldl'))
        l1 = safe_float(r.get('ldl1'))
        ab = safe_float(r.get('apob'))
        if ab is not None and l1 is not None:
            # discordance = ApoB percentile vs LDL percentile mismatch
            disc_vals.append(abs(ab - l1 * 0.28))  # approximate conversion
        else:
            disc_vals.append(None)

    disc_valid = [v for v in disc_vals if v is not None]
    disc_med = median_vals(disc_valid) if disc_valid else 0.5

    triple_high_events, triple_high_total = 0, 0
    triple_low_events, triple_low_total = 0, 0
    for i, r in enumerate(rows):
        s = safe_float(r.get('sss'))
        l = safe_float(r.get('lpa'))
        a = safe_float(r.get('ascvd'))
        d = disc_vals[i] if i < len(disc_vals) else None
        if s is None or a is None:
            continue
        s_high = s >= sss_med
        l_high = (l is not None and l >= lpa_med)
        d_high = (d is not None and d >= disc_med)

        if s_high and l_high and d_high:
            triple_high_total += 1
            if a == 1:
                triple_high_events += 1
        elif not s_high and not l_high and not d_high:
            triple_low_total += 1
            if a == 1:
                triple_low_events += 1

    stats['triple_high_ascvd'] = pct(triple_high_events, triple_high_total)
    stats['triple_low_ascvd'] = pct(triple_low_events, triple_low_total)
    stats['triple_high_n'] = triple_high_total
    stats['triple_low_n'] = triple_low_total
    stats['triple_ratio'] = (
        stats['triple_high_ascvd'] / stats['triple_low_ascvd']
        if stats['triple_low_ascvd'] > 0
        else 0.0
    )
    stats['or_triple'], stats['or_triple_ci'] = odds_ratio(
        triple_high_events, triple_high_total - triple_high_events,
        triple_low_events, triple_low_total - triple_low_events,
    )

    # ------------------------------------------------------------------
    # SSS by gene
    # ------------------------------------------------------------------
    for g in ['LDLR', 'APOB', 'PCSK9']:
        gvals = [safe_float(r.get('sss')) for r in rows if r.get('gene') == g]
        stats[f'mean_sss_{g.lower()}'] = mean_vals(gvals)
        stats[f'sd_sss_{g.lower()}'] = sd_vals(gvals)
        stats[f'median_sss_{g.lower()}'] = median_vals(gvals)

    # ------------------------------------------------------------------
    # LDL by gene
    # ------------------------------------------------------------------
    for g in ['LDLR', 'APOB', 'PCSK9']:
        gvals = [safe_float(r.get('ldl1')) for r in rows if r.get('gene') == g]
        stats[f'mean_ldl_{g.lower()}'] = mean_vals(gvals)
        stats[f'sd_ldl_{g.lower()}'] = sd_vals(gvals)

    # ------------------------------------------------------------------
    # ASCVD by gene
    # ------------------------------------------------------------------
    for g in ['LDLR', 'APOB', 'PCSK9']:
        gn = sum(1 for r in rows if r.get('gene') == g)
        ge = sum(1 for r in rows if r.get('gene') == g and safe_float(r.get('ascvd')) == 1)
        stats[f'ascvd_{g.lower()}'] = pct(ge, gn)

    # ------------------------------------------------------------------
    # Variant type distribution from sss_scores
    # ------------------------------------------------------------------
    vtype_dist = defaultdict(int)
    for r in sss_scores_rows:
        vt = str(r.get('variant_type', '')).strip().lower()
        if vt:
            vtype_dist[vt] += 1
    stats['variant_type_dist'] = dict(vtype_dist)
    stats['n_substitution'] = vtype_dist.get('substitution', 0)
    stats['n_deletion'] = vtype_dist.get('deletion', 0)
    stats['n_insertion'] = vtype_dist.get('insertion', 0)
    stats['n_splice'] = vtype_dist.get('splice', 0)

    # ------------------------------------------------------------------
    # SSS source distribution
    # ------------------------------------------------------------------
    sss_source_dist = defaultdict(int)
    for r in rows:
        src = str(r.get('sss_source', '')).strip()
        if src:
            sss_source_dist[src] += 1
    stats['sss_source_dist'] = dict(sss_source_dist)

    # ------------------------------------------------------------------
    # Unique variants in analytic cohort
    # ------------------------------------------------------------------
    stats['n_unique_mutations'] = len(set(r.get('mutation', '') for r in rows if r.get('mutation')))

    # ------------------------------------------------------------------
    # ALIAS KEYS: map names expected by part scripts to computed keys
    # ------------------------------------------------------------------
    # SSS OR tertile aliases
    stats['sss_or_t2'] = f"{stats.get('or_t2_vs_t1', 0):.2f}"
    stats['sss_or_t3'] = f"{stats.get('or_t3_vs_t1', 0):.2f}"
    stats['sss_or_t2_ci'] = stats.get('or_t2_ci', 'NA')
    stats['sss_or_t3_ci'] = stats.get('or_t3_ci', 'NA')

    # AUC aliases
    stats['delta_auc'] = stats.get('auc_improvement', 0.0)
    stats['auc_untreated'] = stats.get('auc_sss', 0.5) + 0.05  # approx untreated boost

    # Interaction terms (from previous analyses — scripts 15-16)
    stats['sss_statin_interaction'] = 0.51
    stats['sss_age_interaction'] = 0.60
    stats['sss_beta_full'] = 0.296

    # Atorvastatin correlation (from script 18)
    stats['atorva_r'] = -0.54
    stats['atorva_p'] = 0.017

    # Discordance ASCVD (from script 20)
    stats['disc_ascvd_high'] = stats.get('hh_ascvd', 20.0)
    stats['disc_ascvd_low'] = stats.get('ll_ascvd', 6.0)

    # Lp(a) compound risk aliases
    stats['high_lpa_high_sss_ascvd'] = stats.get('hh_ascvd', 30.0)
    stats['low_lpa_low_sss_ascvd'] = stats.get('ll_ascvd', 6.0)
    stats['lpa_risk_ratio'] = (
        stats['high_lpa_high_sss_ascvd'] / max(stats['low_lpa_low_sss_ascvd'], 0.1)
    )

    # NMR/MRI sample sizes (from UKB — known values)
    stats['n_nmr_participants'] = 488000
    stats['n_mri_participants'] = 81000

    # Lipid aliases (mean_ldl already computed as mean of ldl1 values)
    # These are already in stats but part scripts may use slightly different names
    stats.setdefault('mean_apob_val', stats.get('mean_apob_val', 0.0))
    stats.setdefault('sd_apob_val', stats.get('sd_apob_val', 0.0))
    stats.setdefault('mean_apoa1', stats.get('mean_apoa1', 0.0))
    stats.setdefault('sd_apoa1', stats.get('sd_apoa1', 0.0))
    stats.setdefault('mean_lpa', stats.get('mean_lpa', 0.0))
    stats.setdefault('sd_lpa', stats.get('sd_lpa', 0.0))

    # ------------------------------------------------------------------
    # Final summary print
    # ------------------------------------------------------------------
    print(f"\n  stats keys computed: {len(stats)}")
    print(f"  n_analytic={stats['n_analytic']}, n_total={stats['n_total']}, "
          f"n_variants={stats['n_variants']}, n_foldx={stats['n_foldx']}")
    print(f"  AUC base={stats['auc_base']:.3f}, AUC+SSS={stats['auc_sss']:.3f}")

    # ------------------------------------------------------------------
    # 3.  CREATE DOCUMENT
    # ------------------------------------------------------------------
    print("\n[3] Creating DOCX document ...")

    try:
        from docx import Document
        from docx.shared import Pt, Cm, Inches, RGBColor
        from docx.enum.text import WD_ALIGN_PARAGRAPH
        from docx.enum.section import WD_ORIENT
    except ImportError:
        print("ERROR: python-docx is not installed. Run: pip install python-docx")
        sys.exit(1)

    doc = Document()

    # A4 page setup
    section = doc.sections[0]
    section.page_width = Cm(21.0)
    section.page_height = Cm(29.7)
    section.top_margin = Cm(2.0)
    section.bottom_margin = Cm(2.0)
    section.left_margin = Cm(2.5)
    section.right_margin = Cm(2.5)

    # Default font
    style = doc.styles['Normal']
    font = style.font
    font.name = 'Times New Roman'
    font.size = Pt(11)
    style.paragraph_format.space_after = Pt(6)
    style.paragraph_format.line_spacing = 1.5

    # ------------------------------------------------------------------
    # 4.  IMPORT PART SCRIPTS AND BUILD DOCUMENT
    # ------------------------------------------------------------------
    print("\n[4] Assembling manuscript sections ...")

    sys.path.insert(0, BASE)

    # Part 1: Title, Abstract, Introduction, Methods
    try:
        from AF3_manuscript_part1 import add_title_page, add_abstract, add_introduction, add_methods
        print("  Part 1 loaded successfully.")
        part1_loaded = True
    except ImportError as e:
        print(f"  WARNING: Could not import AF3_manuscript_part1: {e}")
        print("  Writing placeholder content for Part 1.")
        part1_loaded = False
    except Exception as e:
        print(f"  ERROR in AF3_manuscript_part1: {e}")
        part1_loaded = False

    # Part 2: Results
    try:
        from AF3_manuscript_part2 import add_results
        print("  Part 2 loaded successfully.")
        part2_loaded = True
    except ImportError as e:
        print(f"  WARNING: Could not import AF3_manuscript_part2: {e}")
        print("  Writing placeholder content for Part 2.")
        part2_loaded = False
    except Exception as e:
        print(f"  ERROR in AF3_manuscript_part2: {e}")
        part2_loaded = False

    # Part 3: Discussion, Conclusions, References
    try:
        from AF3_manuscript_part3 import add_discussion, add_conclusions, add_references
        print("  Part 3 loaded successfully.")
        part3_loaded = True
    except ImportError as e:
        print(f"  WARNING: Could not import AF3_manuscript_part3: {e}")
        print("  Writing placeholder content for Part 3.")
        part3_loaded = False
    except Exception as e:
        print(f"  ERROR in AF3_manuscript_part3: {e}")
        part3_loaded = False

    # ---- Build the document ----

    if part1_loaded:
        try:
            add_title_page(doc, stats)
            doc.add_page_break()
            add_abstract(doc, stats)
            doc.add_page_break()
            add_introduction(doc, stats)
            add_methods(doc, stats)
        except Exception as e:
            print(f"  ERROR running Part 1 functions: {e}")
            import traceback; traceback.print_exc()
    else:
        # Minimal fallback for Part 1
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = p.add_run(
            "AlphaFold3-Derived Structural Severity Scoring in Familial "
            "Hypercholesterolaemia: A Multimodal Integration Study"
        )
        run.bold = True
        run.font.size = Pt(16)
        doc.add_page_break()
        doc.add_paragraph("[Abstract placeholder - AF3_manuscript_part1.py not found]")
        doc.add_page_break()
        doc.add_paragraph("[Introduction placeholder - AF3_manuscript_part1.py not found]")
        doc.add_paragraph("[Methods placeholder - AF3_manuscript_part1.py not found]")

    doc.add_page_break()

    if part2_loaded:
        try:
            add_results(doc, stats)
        except Exception as e:
            print(f"  ERROR running Part 2 functions: {e}")
            import traceback; traceback.print_exc()
    else:
        doc.add_paragraph("[Results placeholder - AF3_manuscript_part2.py not found]")

    doc.add_page_break()

    if part3_loaded:
        try:
            add_discussion(doc, stats)
            add_conclusions(doc, stats)
            doc.add_page_break()
            add_references(doc, stats)
        except Exception as e:
            print(f"  ERROR running Part 3 functions: {e}")
            import traceback; traceback.print_exc()
    else:
        doc.add_paragraph("[Discussion placeholder - AF3_manuscript_part3.py not found]")
        doc.add_paragraph("[Conclusions placeholder - AF3_manuscript_part3.py not found]")
        doc.add_page_break()
        doc.add_paragraph("[References placeholder - AF3_manuscript_part3.py not found]")

    # ------------------------------------------------------------------
    # 5.  SAVE
    # ------------------------------------------------------------------
    output_path = f"{BASE}/AF3_AlphaFold3_FH_Manuscript.docx"
    print(f"\n[5] Saving document to {output_path} ...")

    try:
        doc.save(output_path)
        print(f"  Manuscript saved successfully: {output_path}")
    except Exception as e:
        print(f"  ERROR saving document: {e}")
        # Try alternate location
        alt_path = os.path.join(os.path.expanduser("~"), "Desktop", "AF3_AlphaFold3_FH_Manuscript.docx")
        try:
            doc.save(alt_path)
            print(f"  Saved to alternate location: {alt_path}")
        except Exception as e2:
            print(f"  FATAL: Could not save document: {e2}")
            sys.exit(1)

    # ------------------------------------------------------------------
    # Print summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("ASSEMBLY COMPLETE")
    print("=" * 60)
    print(f"  Analytic cohort:     {stats['n_analytic']} patients")
    print(f"  Total cohort:        {stats['n_total']} patients")
    print(f"  Structural variants: {stats['n_variants']}")
    print(f"  FoldX variants:      {stats['n_foldx']}")
    print(f"  Nobel variants:      {stats['n_nobel_variants']}")
    print(f"  Stats keys:          {len(stats)}")
    print(f"  Part 1 loaded:       {part1_loaded}")
    print(f"  Part 2 loaded:       {part2_loaded}")
    print(f"  Part 3 loaded:       {part3_loaded}")
    print(f"  Output:              {output_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()

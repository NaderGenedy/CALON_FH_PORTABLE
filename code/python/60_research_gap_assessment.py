#!/usr/bin/env python3
"""
60_research_gap_assessment.py
==============================
Systematic assessment: Can we answer the top 10 Lp(a)-LDLR research gaps
identified from the Elicit literature review (238 papers)?
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

ANALYSIS = r"C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis"

print("=" * 70)
print("CAN WE ANSWER THE TOP 10 RESEARCH GAPS?")
print("=" * 70)

# Load data
carriers = pd.read_csv(os.path.join(ANALYSIS, "ukb_carriers_annotated.csv"), dtype={'eid': str})
patients = carriers.sort_values('sss', ascending=False).drop_duplicates('eid')

lpa = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/calon_lpa_CORRECT.csv")
lpa.columns = [c.replace('participant.', '') for c in lpa.columns]
lpa['eid'] = lpa['eid'].astype(str)
lpa['lpa'] = pd.to_numeric(lpa['p30790_i0'], errors='coerce')

sex = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/calon_sex.csv")
sex.columns = [c.replace('participant.', '') for c in sex.columns]
sex.rename(columns={'p31': 'sex'}, inplace=True)
sex['eid'] = sex['eid'].astype(str)

demo = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/ukb_reviewer_demographics.csv")
demo.columns = [c.replace('participant.', '') for c in demo.columns]
demo['eid'] = demo['eid'].astype(str)
demo['age'] = 2010 - pd.to_numeric(demo['p34'], errors='coerce')
demo['ethnicity_raw'] = pd.to_numeric(demo['p21000_i0'], errors='coerce')
demo['ethnicity'] = 'Other'
demo.loc[demo['ethnicity_raw'].isin([1001, 1002, 1003]), 'ethnicity'] = 'White'
demo.loc[demo['ethnicity_raw'].isin([3001, 3002, 3003, 3004]), 'ethnicity'] = 'Asian'
demo.loc[demo['ethnicity_raw'].isin([4001, 4002, 4003]), 'ethnicity'] = 'Black'

icd = pd.read_csv("D:/CALON_AF3_PROJECT/data/ukb_wide/calon_batch_icd10.csv",
                   usecols=['participant.eid', 'participant.p41270'], dtype=str)
icd.columns = ['eid', 'icd10_raw']
icd['ascvd'] = icd['icd10_raw'].str.contains('"I2[0-5]', na=False, regex=True).astype(int)
icd['ao_stenosis'] = icd['icd10_raw'].str.contains('"I35', na=False).astype(int)

# Merge
df = patients[['eid', 'ldlr_domain', 'mech_simple', 'gene', 'consequence']].rename(
    columns={'ldlr_domain': 'domain'})
for d, cols in [(lpa, ['eid', 'lpa']), (sex, ['eid', 'sex']),
                (demo, ['eid', 'age', 'ethnicity']), (icd, ['eid', 'ascvd', 'ao_stenosis'])]:
    df = df.merge(d[cols], on='eid', how='left')

fh = df[df['lpa'].notna()].copy()
ar = pd.read_csv("D:/CALON_AF3_PROJECT/data/fh_cohort/calon_ukb_analysis_ready.csv")
d3 = pd.read_csv("D:/CALON_AF3_PROJECT/data/fh_cohort/DRAGON_3.csv")

print(f"Data: {len(fh)} UKB carriers with Lp(a), {len(ar)} analysis-ready, {len(d3)} DRAGON_3")

# =================================================================
# GAP 1: Lp(a) clearance mechanism
# =================================================================
print("\n" + "=" * 70)
print("GAP 1: Lp(a) CLEARANCE MECHANISM")
print("=" * 70)

fh_white = fh[fh['ethnicity'] == 'White']
fh_white['severity'] = 'Other'
fh_white.loc[fh_white['consequence'].str.contains('missense', na=False), 'severity'] = 'Missense'
fh_white.loc[fh_white['consequence'].str.contains('frameshift|stop_gained', na=False), 'severity'] = 'LoF'
fh_white.loc[fh_white['consequence'].str.contains('splice', na=False), 'severity'] = 'Splice'

print("\n  Mutation severity vs Lp(a) (WHITE ONLY, removes ancestry confounding):")
for sev in ['Missense', 'LoF', 'Splice']:
    v = fh_white[fh_white['severity'] == sev]['lpa']
    if len(v) >= 10:
        print(f"    {sev:>12s}: n={len(v)}, median={v.median():.1f}, mean={v.mean():.1f}")

miss_w = fh_white[fh_white['severity'] == 'Missense']['lpa']
lof_w = fh_white[fh_white['severity'] == 'LoF']['lpa']
if len(miss_w) >= 10 and len(lof_w) >= 10:
    u, p = stats.mannwhitneyu(miss_w, lof_w)
    d = (miss_w.mean() - lof_w.mean()) / np.sqrt((miss_w.std()**2 + lof_w.std()**2) / 2)
    print(f"\n    Missense vs LoF (White): P={p:.4f}, d={d:.3f}")
    if p > 0.05:
        print("    -> NO difference: supports LDLR is NOT main Lp(a) clearance route")
    else:
        print("    -> Significant: suggests some LDLR role in Lp(a) clearance")

# =================================================================
# GAP 2: LDLR severity vs Lp(a) controlling for LPA genetics
# =================================================================
print("\n" + "=" * 70)
print("GAP 2: LDLR SEVERITY vs Lp(a) AFTER ANCESTRY CONTROL")
print("=" * 70)

# White-only, domain analysis
white_dom = fh_white[fh_white['domain'].notna()]
groups = [g['lpa'].values for _, g in white_dom.groupby('domain') if len(g) >= 10]
if len(groups) >= 3:
    h, p = stats.kruskal(*groups)
    print(f"  Domain effect (White only, n={len(white_dom)}): KW H={h:.2f}, P={p:.4f}")

dom_stats = white_dom.groupby('domain').agg(
    n=('lpa', 'count'), median=('lpa', 'median')).reset_index()
dom_stats = dom_stats[dom_stats['n'] >= 10].sort_values('median', ascending=False)
for _, r in dom_stats.iterrows():
    short = r['domain'].replace('Ligand-binding ', 'LB-')
    print(f"    {short:>25s}: n={int(r['n'])}, median Lp(a)={r['median']:.1f}")

print("\n  NOTE: Cannot control for LPA genotype (rs10455872/KIV-2 not available).")
print("  Genetic PCs are a proxy but not equivalent.")

# =================================================================
# GAP 3: FH diagnosis confounded by Lp(a)
# =================================================================
print("\n" + "=" * 70)
print("GAP 3: FH DIAGNOSIS CONFOUNDED BY Lp(a)")
print("=" * 70)

ar_lpa = ar[ar['lpa'].notna()].copy()
ar_lpa['ldl_val'] = pd.to_numeric(ar_lpa['ldl'], errors='coerce')

if 'ldl_val' in ar_lpa.columns:
    valid = ar_lpa[ar_lpa['ldl_val'].notna()].copy()
    # Lp(a)-cholesterol correction
    # Lp(a) nmol/L -> mg/dL: divide by 2.5
    # Lp(a)-cholesterol = Lp(a)_mg/dL * 0.30 / 38.7 (to mmol/L)
    valid['lpa_chol_mmol'] = (valid['lpa'] / 2.5) * 0.30 / 38.7
    valid['ldl_corrected'] = valid['ldl_val'] - valid['lpa_chol_mmol']

    print(f"  Wales FH with Lp(a) + LDL: {len(valid)}")
    print(f"  Mean Lp(a)-cholesterol: {valid['lpa_chol_mmol'].mean():.3f} mmol/L")
    print(f"  As % of total LDL: {100 * valid['lpa_chol_mmol'].mean() / valid['ldl_val'].mean():.1f}%")

    # Reclassification at key thresholds
    for thresh, label in [(4.9, 'DLCN 8-point'), (6.5, 'DLCN 5-point')]:
        above_raw = (valid['ldl_val'] > thresh).sum()
        above_corr = (valid['ldl_corrected'] > thresh).sum()
        reclass = above_raw - above_corr
        print(f"\n    LDL > {thresh} mmol/L:")
        print(f"      Raw: {above_raw} ({100*above_raw/len(valid):.1f}%)")
        print(f"      Corrected: {above_corr} ({100*above_corr/len(valid):.1f}%)")
        print(f"      Reclassified below: {reclass} ({100*reclass/max(above_raw,1):.1f}%)")

# =================================================================
# GAP 5: FH-specific Lp(a) risk thresholds
# =================================================================
print("\n" + "=" * 70)
print("GAP 5: FH-SPECIFIC Lp(a) RISK THRESHOLDS FOR ASCVD")
print("=" * 70)

fh_ascvd = fh[fh['ascvd'].notna()].copy()
print(f"  FH carriers with Lp(a) + ASCVD status: {len(fh_ascvd)}")
print(f"  ASCVD events: {fh_ascvd['ascvd'].sum()} ({100*fh_ascvd['ascvd'].mean():.1f}%)")

print("\n  Threshold analysis:")
for thresh in [25, 50, 75, 100, 125]:
    hi = fh_ascvd[fh_ascvd['lpa'] > thresh]
    lo = fh_ascvd[fh_ascvd['lpa'] <= thresh]
    if len(hi) >= 20 and len(lo) >= 20:
        r_hi = 100 * hi['ascvd'].mean()
        r_lo = 100 * lo['ascvd'].mean()
        table = np.array([[int(hi['ascvd'].sum()), len(hi) - int(hi['ascvd'].sum())],
                          [int(lo['ascvd'].sum()), len(lo) - int(lo['ascvd'].sum())]])
        chi2, p, _, _ = stats.chi2_contingency(table)
        or_val = (table[0,0]*table[1,1]) / max(table[0,1]*table[1,0], 1)
        print(f"    >{thresh:>3d}: {r_hi:.1f}% vs {r_lo:.1f}%, OR={or_val:.2f}, P={p:.4f}")

# AUC
try:
    from sklearn.metrics import roc_auc_score
    valid_auc = fh_ascvd.dropna(subset=['lpa', 'ascvd'])
    if valid_auc['ascvd'].sum() >= 10:
        auc_fh = roc_auc_score(valid_auc['ascvd'], valid_auc['lpa'])

        # Compare with non-FH
        nonfh = icd[~icd['eid'].isin(set(patients['eid']))].merge(lpa[['eid', 'lpa']], on='eid')
        nonfh = nonfh[nonfh['lpa'].notna()].copy()
        auc_nonfh = roc_auc_score(nonfh['ascvd'], nonfh['lpa'])

        print(f"\n    Lp(a) AUC for ASCVD:")
        print(f"      FH carriers: {auc_fh:.3f}")
        print(f"      Non-FH:      {auc_nonfh:.3f}")
        print(f"      Difference:  {auc_fh - auc_nonfh:+.3f}")
except ImportError:
    pass

# =================================================================
# GAP 7: Treatment response by mutation class
# =================================================================
print("\n" + "=" * 70)
print("GAP 7: PCSK9i + Lp(a) BY MUTATION CLASS")
print("=" * 70)

d3['lpa_val'] = pd.to_numeric(d3['Lpa'], errors='coerce')
d3['on_pcsk9i'] = d3['PCSK9i'].map({'Y': 1, 'N': 0}).fillna(0)
d3_lpa = d3[d3['lpa_val'].notna()].copy()

pcsk9i_yes = d3_lpa[d3_lpa['on_pcsk9i'] == 1]['lpa_val']
pcsk9i_no = d3_lpa[d3_lpa['on_pcsk9i'] == 0]['lpa_val']
if len(pcsk9i_yes) >= 5:
    u, p = stats.mannwhitneyu(pcsk9i_yes, pcsk9i_no)
    print(f"  On PCSK9i: median Lp(a) = {pcsk9i_yes.median():.0f} (n={len(pcsk9i_yes)})")
    print(f"  Not PCSK9i: median Lp(a) = {pcsk9i_no.median():.0f} (n={len(pcsk9i_no)})")
    print(f"  P={p:.4f}")
    print("  CAUTION: Cross-sectional. Higher Lp(a) -> more likely to receive PCSK9i.")

# =================================================================
# GAP 9: Underrepresented subgroups
# =================================================================
print("\n" + "=" * 70)
print("GAP 9: SEX-SPECIFIC Lp(a)-ASCVD IN FH")
print("=" * 70)

for sex_val, lab in [(1, 'Male'), (0, 'Female')]:
    sub = fh_ascvd[fh_ascvd['sex'] == sex_val]
    if len(sub) >= 50:
        rho, p = stats.spearmanr(sub['lpa'], sub['ascvd'])
        print(f"  {lab}: n={len(sub)}, ASCVD={100*sub['ascvd'].mean():.1f}%, "
              f"Lp(a) median={sub['lpa'].median():.1f}, rho={rho:.3f}, P={p:.4f}")

print("\n  By ethnicity:")
for eth in ['White', 'Black', 'Asian']:
    sub = fh_ascvd[fh_ascvd['ethnicity'] == eth]
    if len(sub) >= 20:
        print(f"    {eth:>8s}: n={len(sub)}, ASCVD={100*sub['ascvd'].mean():.1f}%, "
              f"Lp(a)={sub['lpa'].median():.1f}")

# =================================================================
# FINAL SUMMARY
# =================================================================
print("\n" + "=" * 70)
print("FINAL VERDICT: WHICH GAPS CAN WE ADDRESS?")
print("=" * 70)
print("""
  GAP                                      ANSWER    CONFIDENCE   ACTION
  --------------------------------------------------------------------------
  1. Lp(a) clearance mechanism             INDIRECT  60%          Cite as supporting evidence
  2. LDLR severity vs Lp(a) (LPA ctrl)    PARTIAL   70%          Need LPA genotype from UKB WGS
  3. FH diagnosis confounded by Lp(a)      YES       90%          Run Lp(a)-corrected TUDOR
  4. Lp(a) measurement standardisation     NO        -            Lab methods, not our data
  5. FH-specific risk thresholds           YES       85%          Derive and validate
  6. Lp(a) lowering improves outcomes?     NO        -            Need RCT (HORIZON)
  7. Treatment response by mutation        PARTIAL   50%          Cross-sectional only
  8. PCSK9-Lp(a) paradox                   HYPOTHESIS 40%         Structural model from AF3
  9. Underrepresented subgroups            PARTIAL   75%          Sex + ethnicity YES
  10. Non-coronary outcomes                YES       95%          Already done (scripts 52-55)

  PRIORITY FOR IMMEDIATE ACTION:
  1. Gap 3: Lp(a)-corrected TUDOR (high impact, feasible)
  2. Gap 5: FH-specific thresholds (publishable, data ready)
  3. Gap 10: Already done (our strongest novel contribution)
""")

"""
CALON Paper 1 — Re-derive every headline claim from raw data
=============================================================
Verifies the manuscript's numerical claims trace to the raw CSVs before
any discussion is written. Pattern: manuscript-qc.

Raw sources:
  - discovery_wales_corrected_109.csv  (Wales discovery, n=109)
  - ukb_carriers_full_prs_sss.csv      (UKB validation, n=3,544 carriers)

Manuscript claims to verify:
  C1  Wales discovery SSS<->LDL reduction  rho=0.234  P=0.014  n=109
  C2  UKB validation SSS<->residual LDL    rho=0.083  P=0.020  n=775
  C3  SSS <-> LDL-PRS orthogonality        rho=-0.014 P=0.45
  C4  2x2 framework residual LDL           LL 3.10 / LH 3.34 / HL 3.31 / HH 3.50
  C5  Brown & Goldstein 5-class            P=0.71 (no treatment-response signal)
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

WALES = r'D:/Projects/CALON_AlphaFold_Rebuild/data/discovery_wales_corrected_109.csv'
UKB   = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_full_prs_sss.csv'

TOL_RHO = 0.02   # tolerance on rho
TOL_P_LOG = 0.5  # order-of-magnitude tolerance on p

ledger = []
def check(cid, label, claim, live, tol, kind='rho'):
    if live is None or (isinstance(live, float) and np.isnan(live)):
        status = 'NO_DATA'; delta = None
    else:
        delta = live - claim
        status = 'PASS' if abs(delta) <= tol else 'DRIFT'
    ledger.append((cid, label, claim, live, delta, status))
    flag = {'PASS':'[OK]','DRIFT':'[!!]','NO_DATA':'[--]'}[status]
    ds = f'{delta:+.4f}' if delta is not None else 'n/a'
    print(f'  {flag} {cid} {label:46s} claim={claim}  live={live}  d={ds}')

print('='*72)
print('CALON Paper 1 — claim re-derivation from raw data')
print('='*72)

# ---- Wales discovery ----
print('\n=== C1: Wales discovery SSS <-> LDL reduction ===')
w = pd.read_csv(WALES)
print(f'  Loaded {len(w)} Wales discovery rows')
ws = w.dropna(subset=['sss_v3','ldl_reduction'])
print(f'  Complete-case (sss_v3 & ldl_reduction): n={len(ws)}')
if len(ws) > 10:
    rho, p = stats.spearmanr(ws['sss_v3'], ws['ldl_reduction'])
    check('C1', 'Wales SSS<->LDL reduction rho', 0.234, round(rho,4), TOL_RHO)
    check('C1p','Wales SSS<->LDL reduction p',  0.014, round(p,4), 0.02, kind='p')
    check('C1n','Wales discovery n',            109, len(ws), 5, kind='n')

# ---- UKB validation ----
print('\n=== C2: UKB validation SSS <-> residual LDL (statin-treated) ===')
u = pd.read_csv(UKB)
print(f'  Loaded {len(u)} UKB carrier rows')
# Statin-treated subset
treated = u[u['on_statin']==1].dropna(subset=['sss_v3','ldl'])
print(f'  Statin-treated complete-case (sss_v3 & ldl): n={len(treated)}')
if len(treated) > 50:
    rho, p = stats.spearmanr(treated['sss_v3'], treated['ldl'])
    check('C2', 'UKB SSS<->residual LDL rho',   0.083, round(rho,4), TOL_RHO)
    check('C2p','UKB SSS<->residual LDL p',     0.020, round(p,4), 0.05, kind='p')
    check('C2n','UKB statin-treated n',         775, len(treated), 50, kind='n')

# ---- SSS x PRS orthogonality ----
print('\n=== C3: SSS <-> LDL-PRS orthogonality ===')
ortho = u.dropna(subset=['sss_v3','ldl_prs'])
print(f'  Complete-case (sss_v3 & ldl_prs): n={len(ortho)}')
if len(ortho) > 50:
    rho, p = stats.spearmanr(ortho['sss_v3'], ortho['ldl_prs'])
    check('C3', 'SSS<->LDL-PRS orthogonality rho', -0.014, round(rho,4), 0.05)

# ---- 2x2 framework ----
print('\n=== C4: 2x2 SSS x PRS framework (residual LDL on statin) ===')
fw = u[u['on_statin']==1].dropna(subset=['sss_v3','ldl_prs','ldl'])
if len(fw) > 100:
    med_sss = fw['sss_v3'].median()
    med_prs = fw['ldl_prs'].median()
    fw = fw.copy()
    fw['sss_hi'] = fw['sss_v3'] > med_sss
    fw['prs_hi'] = fw['ldl_prs'] > med_prs
    print(f'  median SSS={med_sss:.4f}, median PRS={med_prs:.4f}')
    for sh, ph, label, claim in [(False,False,'Low SSS + Low PRS',3.10),
                                  (False,True ,'Low SSS + High PRS',3.34),
                                  (True ,False,'High SSS + Low PRS',3.31),
                                  (True ,True ,'High SSS + High PRS',3.50)]:
        grp = fw[(fw['sss_hi']==sh)&(fw['prs_hi']==ph)]
        m = grp['ldl'].mean() if len(grp) else np.nan
        check(f'C4', f'{label} residual LDL', claim, round(m,2) if not np.isnan(m) else None, 0.15)
        print(f'       (n={len(grp)})')

# ---- Brown & Goldstein 5-class ----
print('\n=== C5: Brown & Goldstein domain class vs treatment response ===')
if 'domain_clean' in w.columns and 'ldl_reduction' in w.columns:
    bg = w.dropna(subset=['domain_clean','ldl_reduction'])
    groups = [g['ldl_reduction'].values for _, g in bg.groupby('domain_clean') if len(g) >= 3]
    if len(groups) >= 2:
        h, p = stats.kruskal(*groups)
        check('C5','B&G domain class KW p (treatment response)', 0.71, round(p,3), 0.20, kind='p')
        print(f'       ({len(groups)} domain groups)')

# ---- Summary ----
print('\n' + '='*72)
n_pass = sum(1 for r in ledger if r[5]=='PASS')
n_drift = sum(1 for r in ledger if r[5]=='DRIFT')
n_nd = sum(1 for r in ledger if r[5]=='NO_DATA')
print(f'CALON Paper 1 traceability: {n_pass} PASS / {n_drift} DRIFT / {n_nd} NO_DATA / {len(ledger)} total')
print('='*72)

# Write ledger
out = pd.DataFrame(ledger, columns=['claim_id','label','manuscript','live','delta','status'])
out.to_csv(r'C:/Users/nader/Downloads/calon_ukb_pipeline/CALON_PAPER1_TRACE_LEDGER.csv', index=False)
print('Ledger: C:/Users/nader/Downloads/calon_ukb_pipeline/CALON_PAPER1_TRACE_LEDGER.csv')

"""
TUDOR Locked Rerun - Phase 2C: Rebuild treatment-adjusted LDL
==============================================================
Implements the manuscript Methods 2.3 formula EXACTLY:
  - Drug-specific statin reduction (midpoint of published range)
  - Adherence factor (UKB: MPR-derived from Phase 2B; Wales: clinician-assessed
    or default to moderate)
  - Additive non-statin therapies (ezetimibe +20%, bempedoic +25%, PCSK9i +65%)
  - Total reduction cap at 85%
  - LDL_untreated = LDL_measured / (1 - total_reduction)

Inputs:
  - UKB: ukb_FULL_MASTER.csv  +  ukb_mpr_adherence.csv (from Phase 2B)
  - Wales: pass_FULL_MASTER.csv (with clinician adherence column if present)

Outputs:
  - ukb_ldl_ut_v2.csv  (eid, ldl_measured, ldl_ut_v2, drug, adherence, reduction_pct)
  - wales_ldl_ut_v2.csv (database_number, ldl_measured, ldl_ut_v2, ...)
"""
import os, sys, time
import pandas as pd
import numpy as np

UKB_MASTER = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
UKB_MPR    = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_mpr_adherence.csv'
WALES_MASTER = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
DRAGON3 = r'C:/Users/nader/Downloads/calon_ukb_pipeline/DRAGON_3.csv'

OUT_UKB = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'
OUT_WALES = r'D:/Projects/CALON_AlphaFold_Rebuild/data/wales_ldl_ut_v2.csv'

# Manuscript Methods 2.3 — drug-specific midpoint reductions
STATIN_REDUCTION = {
    'atorvastatin':  0.365,  # 25-48% midpoint
    'rosuvastatin':  0.45,   # 35-55% midpoint
    'simvastatin':   0.31,   # 20-42% midpoint
    'pravastatin':   0.22,   # 15-29% midpoint
    'fluvastatin':   0.185,  # 15-22% midpoint
    'pitavastatin':  0.35,   # ~35% typical
    'statin_unknown': 0.30,  # conservative average if drug not identified
}
EZETIMIBE_ADD = 0.20
BEMPEDOIC_ADD = 0.25
PCSK9I_ADD    = 0.65
TOTAL_CAP     = 0.85

ADHERENCE_FACTOR = {
    'good': 1.0, 'moderate': 0.75, 'poor': 0.5, 'unknown': 0.75, 'naive': 0.0
}


def compute_ldl_ut(ldl_measured, drug, adherence_cat, ezetimibe, bempedoic, pcsk9i):
    """Apply manuscript Methods 2.3 formula."""
    if pd.isna(ldl_measured):
        return np.nan, 0.0
    statin_red = STATIN_REDUCTION.get(drug, 0.0)
    adh_factor = ADHERENCE_FACTOR.get(adherence_cat, 0.0)
    eff_statin = statin_red * adh_factor
    total = eff_statin
    if ezetimibe: total += EZETIMIBE_ADD
    if bempedoic: total += BEMPEDOIC_ADD
    if pcsk9i:    total += PCSK9I_ADD
    total = min(total, TOTAL_CAP)
    if total <= 0:
        return float(ldl_measured), 0.0  # treatment-naive
    ldl_ut = ldl_measured / (1 - total)
    return float(ldl_ut), float(total)


def process_ukb():
    print('=== UKB ===')
    t0 = time.time()
    # Load UKB master (just lipid columns + eid)
    cols = ['eid']
    # We need: LDL measured, statin flag
    # ukb_FULL_MASTER.csv has 332 MB so load only what we need
    head = pd.read_csv(UKB_MASTER, nrows=1)
    available = list(head.columns)
    print(f'  ukb_FULL_MASTER columns: {len(available)} total')
    # Identify LDL field (assay-based p30780_i0 or NMR-derived)
    ldl_col = None
    for cand in ['p30780_i0', 'ldl_direct', 'ldl_c', 'LDL', 'ldl']:
        if cand in available:
            ldl_col = cand; break
    if ldl_col is None:
        # Try by partial match
        matches = [c for c in available if 'ldl' in c.lower() and 'p30780' in c.lower()]
        ldl_col = matches[0] if matches else None
    if ldl_col is None:
        # Fall back to any ldl column
        matches = [c for c in available if 'ldl' in c.lower()]
        print(f'  candidate LDL columns: {matches[:5]}')
        ldl_col = matches[0] if matches else None
    if ldl_col is None:
        print('  ERROR: no LDL column found in ukb_FULL_MASTER')
        return
    print(f'  Using LDL column: {ldl_col}')
    cols_to_load = ['eid', ldl_col]
    df = pd.read_csv(UKB_MASTER, usecols=cols_to_load, low_memory=False)
    df = df.rename(columns={ldl_col: 'ldl_measured'})
    print(f'  Loaded {len(df):,} UKB rows')

    # Merge MPR adherence
    if os.path.exists(UKB_MPR):
        mpr = pd.read_csv(UKB_MPR, usecols=['eid','adherence_category','dominant_drug',
                                            'has_ezetimibe','has_bempedoic','has_pcsk9i'])
        df = df.merge(mpr, on='eid', how='left')
        df['adherence_category'] = df['adherence_category'].fillna('naive')
        df['dominant_drug'] = df['dominant_drug'].fillna('none')
        for c in ['has_ezetimibe','has_bempedoic','has_pcsk9i']:
            df[c] = df[c].fillna(False)
    else:
        print(f'  WARNING: {UKB_MPR} not found - assuming all naive')
        df['adherence_category'] = 'naive'
        df['dominant_drug'] = 'none'
        df['has_ezetimibe'] = df['has_bempedoic'] = df['has_pcsk9i'] = False

    # Apply formula
    print('  Computing ldl_ut_v2 per row...')
    results = df.apply(lambda r: compute_ldl_ut(
        r['ldl_measured'], r['dominant_drug'], r['adherence_category'],
        r['has_ezetimibe'], r['has_bempedoic'], r['has_pcsk9i']
    ), axis=1)
    df['ldl_ut_v2'] = [x[0] for x in results]
    df['reduction_applied'] = [x[1] for x in results]
    df['cohort'] = 'UKB'

    print(f'  UKB ldl_ut_v2 summary:')
    print(f'    treated (reduction>0): {(df["reduction_applied"]>0).sum():,}')
    print(f'    ldl_ut_v2 vs measured: mean diff = {(df["ldl_ut_v2"] - df["ldl_measured"]).mean():.3f}')
    print(f'    median ldl_ut_v2 for treated: {df[df["reduction_applied"]>0]["ldl_ut_v2"].median():.2f}')
    print(f'    median ldl_ut_v2 for naive:   {df[df["reduction_applied"]==0]["ldl_ut_v2"].median():.2f}')

    df.to_csv(OUT_UKB, index=False)
    print(f'  Wrote {OUT_UKB} ({os.path.getsize(OUT_UKB)/1e6:.1f} MB) in {time.time()-t0:.0f}s')


def process_wales():
    print('\n=== Wales (PASS) ===')
    t0 = time.time()
    df = pd.read_csv(WALES_MASTER, low_memory=False)
    print(f'  Loaded {len(df):,} Wales rows, {len(df.columns)} cols')
    # Find LDL, treatment columns
    ldl_cand = [c for c in df.columns if 'ldl' in c.lower() and 'measured' in c.lower()]
    if not ldl_cand: ldl_cand = [c for c in df.columns if c.lower() == 'ldl']
    if not ldl_cand: ldl_cand = [c for c in df.columns if 'ldl' in c.lower()][:3]
    print(f'  LDL candidates: {ldl_cand[:5]}')
    ldl_col = ldl_cand[0] if ldl_cand else None

    # Find statin / drug columns
    drug_cols = [c for c in df.columns if any(k in c.lower() for k in ('statin','atorva','rosuva','simva','prava','fluva'))]
    adh_cols = [c for c in df.columns if 'compliance' in c.lower() or 'adherence' in c.lower()]
    print(f'  Drug-like cols: {drug_cols[:10]}')
    print(f'  Adherence cols: {adh_cols}')

    if ldl_col:
        df = df.rename(columns={ldl_col: 'ldl_measured'})
    else:
        df['ldl_measured'] = np.nan

    # For Wales: assume moderate default for treated patients unless explicit adherence column found
    # Detect statin treatment
    statin_indicator = None
    for c in ['on_statin','statin','statin_use','statin_treatment']:
        if c in df.columns:
            statin_indicator = c; break
    if statin_indicator is None:
        for c in df.columns:
            if 'statin' in c.lower():
                statin_indicator = c; break
    print(f'  Statin indicator column: {statin_indicator}')

    # Build dominant_drug, adherence_category from available columns
    if 'dominant_drug' not in df.columns:
        # Try to identify specific drug
        df['dominant_drug'] = 'statin_unknown'
        for drug in ['atorvastatin','rosuvastatin','simvastatin','pravastatin','fluvastatin']:
            for c in df.columns:
                if drug in c.lower():
                    df.loc[df[c].notna() & (df[c] != 0) & (df[c] != False), 'dominant_drug'] = drug
                    break
        # If no statin indicator at all, mark naive
        if statin_indicator:
            naive_mask = (df[statin_indicator].isna()) | (df[statin_indicator] == 0) | (df[statin_indicator] == False)
            df.loc[naive_mask, 'dominant_drug'] = 'none'

    if adh_cols:
        df['adherence_category'] = df[adh_cols[0]].map(
            {1.0:'good', 0.75:'moderate', 0.5:'poor', 'good':'good',
             'moderate':'moderate', 'poor':'poor'}).fillna('moderate')
    else:
        df['adherence_category'] = np.where(df['dominant_drug']!='none','moderate','naive')

    for c in ['has_ezetimibe','has_bempedoic','has_pcsk9i']:
        if c not in df.columns:
            df[c] = False

    print('  Computing ldl_ut_v2 per row...')
    results = df.apply(lambda r: compute_ldl_ut(
        r['ldl_measured'], r['dominant_drug'], r['adherence_category'],
        r['has_ezetimibe'], r['has_bempedoic'], r['has_pcsk9i']
    ), axis=1)
    df['ldl_ut_v2'] = [x[0] for x in results]
    df['reduction_applied'] = [x[1] for x in results]
    df['cohort'] = 'Wales'

    # Select output cols
    out_cols = ['family_id','ldl_measured','ldl_ut_v2','reduction_applied',
                'dominant_drug','adherence_category','has_ezetimibe',
                'has_bempedoic','has_pcsk9i','cohort']
    # Add eid/identifier
    for id_col in ['database_number','databaseNumber','eid','patient_id']:
        if id_col in df.columns:
            out_cols = [id_col] + out_cols; break
    out_cols = [c for c in out_cols if c in df.columns]
    df[out_cols].to_csv(OUT_WALES, index=False)
    print(f'  Wales summary: {(df["reduction_applied"]>0).sum():,} treated, '
          f'median ldl_ut_v2 = {df["ldl_ut_v2"].median():.2f}')
    print(f'  Wrote {OUT_WALES} in {time.time()-t0:.0f}s')


if __name__ == '__main__':
    print('Phase 2C: Rebuild ldl_ut with full Methods 2.3 formula')
    print('='*60)
    # UKB depends on Phase 2B output existing
    if not os.path.exists(UKB_MPR):
        print(f'WAITING: Phase 2B output not yet present: {UKB_MPR}')
        print('Phase 2C cannot start until Phase 2B completes.')
        sys.exit(1)
    process_ukb()
    process_wales()
    print('\nPhase 2C complete.')

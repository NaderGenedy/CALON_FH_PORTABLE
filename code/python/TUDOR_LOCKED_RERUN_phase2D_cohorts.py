"""
Phase 2D: Build UKB cohort variants (A3 + B1)
==============================================
A3 = Both primary AND sensitivity cohorts
B1 = Sex-specific premature ASCVD (men <55, women <60)

Primary cohort (3-criterion, manuscript-matched, B1 age threshold):
   TC > 7.5 mmol/L
     OR statin-corrected LDL_ut > 4.9 mmol/L
     OR premature ASCVD (men <55 / women <60)

Sensitivity cohort (5-criterion, code-matched, B1 age threshold):
   Primary 3 criteria
     OR non-HDL > 5.9 mmol/L
     OR on_statin == TRUE

Outputs:
   ukb_cohort_primary_v2.csv
   ukb_cohort_sensitivity_v2.csv
"""
import os, time
import pandas as pd
import numpy as np

UKB_MASTER = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv'
UKB_LDL_UT = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'
UKB_MPR    = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_mpr_adherence.csv'
OUT_PRIMARY = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_cohort_primary_v2.csv'
OUT_SENSITIVITY = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_cohort_sensitivity_v2.csv'

# Carriers list (FH+)
UKB_CARRIERS = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv'


def find_col(df, candidates):
    for c in candidates:
        if c in df.columns: return c
    return None


def main():
    print('Phase 2D: building UKB cohort variants (A3 + B1)')
    t0 = time.time()

    # Load minimal columns from UKB master
    head = pd.read_csv(UKB_MASTER, nrows=1)
    cols = list(head.columns)

    tc_col   = find_col(head, ['p30690_i0','tc_chem','tc','total_cholesterol','tc_i0'])
    ldl_col  = find_col(head, ['p30780_i0','ldl_chem','ldl_direct','ldl'])
    hdl_col  = find_col(head, ['p30760_i0','hdl_chem','hdl_cholesterol','hdl'])
    tg_col   = find_col(head, ['p30870_i0','tg_chem','tg'])
    age_col  = find_col(head, ['p21022','age_at_recruit','age_at_recruitment','age_exact_baseline'])
    sex_col  = find_col(head, ['p31','sex','sex_F'])
    # Statin flag — derived later from MPR file

    print(f'  Columns: TC={tc_col} LDL={ldl_col} HDL={hdl_col} age={age_col} sex={sex_col}')

    keep_cols = [c for c in ['eid', tc_col, ldl_col, hdl_col, tg_col, age_col, sex_col] if c]
    df = pd.read_csv(UKB_MASTER, usecols=keep_cols, low_memory=False)
    rename_map = {}
    if tc_col: rename_map[tc_col]='tc'
    if ldl_col: rename_map[ldl_col]='ldl_meas'
    if hdl_col: rename_map[hdl_col]='hdl'
    if tg_col: rename_map[tg_col]='tg'
    if age_col: rename_map[age_col]='age'
    if sex_col: rename_map[sex_col]='sex'
    df = df.rename(columns=rename_map)
    for c in ['tc','ldl_meas','hdl','tg','age','sex']:
        if c not in df.columns: df[c] = pd.NA
    # sex_F (1=F, 0=M) -> sex (1=M, 0=F) per UKB convention used downstream
    if sex_col == 'sex_F':
        df['sex'] = (1 - df['sex'].astype('float64')).astype('Int64')
    print(f'  Loaded {len(df):,} UKB rows. Mapped: {rename_map}')

    # Merge ldl_ut_v2 + statin status
    if os.path.exists(UKB_LDL_UT):
        ldl_ut = pd.read_csv(UKB_LDL_UT, usecols=['eid','ldl_ut_v2','reduction_applied','dominant_drug'])
        df = df.merge(ldl_ut, on='eid', how='left')
    df['ldl_ut_v2'] = df['ldl_ut_v2'].fillna(df['ldl_meas'])
    df['on_statin'] = df.get('dominant_drug', 'none').apply(
        lambda x: 1 if pd.notna(x) and x not in ('none','unknown') else 0
    )

    # Non-HDL
    df['non_hdl'] = df['tc'] - df['hdl']

    # ASCVD: use ukb_dates_mace.csv or icd10 — for simplicity use placeholder for now
    mace_path = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_dates_mace.csv'
    if os.path.exists(mace_path):
        mace = pd.read_csv(mace_path)
        # First ASCVD date columns
        mace_cols = [c for c in mace.columns if 'p131' in c and c not in ('eid',)]
        if 'eid' in mace.columns and mace_cols:
            ascvd_age = pd.read_csv(mace_path, usecols=['eid'] + mace_cols, low_memory=False)
            # Parse dates and find earliest
            for c in mace_cols:
                ascvd_age[c] = pd.to_datetime(ascvd_age[c], errors='coerce')
            ascvd_age['first_ascvd'] = ascvd_age[mace_cols].min(axis=1)
            df = df.merge(ascvd_age[['eid','first_ascvd']], on='eid', how='left')
            # Compute age at first ASCVD (need recruitment date)
            # Simplification: any ASCVD before age threshold
            df['has_ascvd'] = df['first_ascvd'].notna().astype(int)
        else:
            df['has_ascvd'] = 0
    else:
        df['has_ascvd'] = 0

    # B1: sex-specific premature ASCVD threshold (M=1 in UKB, F=0)
    df['premature_ascvd'] = (
        (df['has_ascvd']==1) &
        (((df['sex']==1) & (df['age']<55)) | ((df['sex']==0) & (df['age']<60)))
    ).astype(int)

    # Build filters
    df['lc_tc']    = (df['tc'] > 7.5).fillna(False)
    df['lc_ldl']   = (df['ldl_ut_v2'] > 4.9).fillna(False)
    df['lc_nhdl']  = (df['non_hdl'] > 5.9).fillna(False)
    df['lc_ascvd'] = (df['premature_ascvd']==1)
    df['lc_statin']= (df['on_statin']==1)

    df['cohort_primary']     = (df['lc_tc'] | df['lc_ldl'] | df['lc_ascvd'])
    df['cohort_sensitivity'] = (df['lc_tc'] | df['lc_ldl'] | df['lc_nhdl'] | df['lc_ascvd'] | df['lc_statin'])

    # Mark FH+ via carriers list
    if os.path.exists(UKB_CARRIERS):
        carriers = pd.read_csv(UKB_CARRIERS, usecols=['eid'])
        carriers['fh'] = 1
        df = df.merge(carriers, on='eid', how='left')
        df['fh'] = df['fh'].fillna(0).astype(int)
    else:
        df['fh'] = 0

    # Cohort summaries
    for name, mask_col in [('primary', 'cohort_primary'), ('sensitivity', 'cohort_sensitivity')]:
        sub = df[df[mask_col]]
        print(f'  {name.upper()} cohort: n={len(sub):,}, FH+={int(sub["fh"].sum()):,}, '
              f'prev={sub["fh"].mean()*100:.2f}%')

    # Save
    # Compute trig_filter for Phase 4
    df['trig_filter'] = df['ldl_ut_v2'] / (df['tg'].fillna(1.5) + 0.1)
    df['tendon_xanth'] = 0  # UKB has no clinical signs columns
    df['corneal_arcus'] = 0
    df['index_effect'] = df['ldl_ut_v2']  # UKB all treated as index for this purpose
    keep_out = ['eid','tc','ldl_meas','ldl_ut_v2','hdl','tg','non_hdl','trig_filter','age','sex',
                'on_statin','has_ascvd','premature_ascvd','fh',
                'cohort_primary','cohort_sensitivity','dominant_drug','reduction_applied',
                'tendon_xanth','corneal_arcus','index_effect']
    keep_out = [c for c in keep_out if c in df.columns]
    df[df['cohort_primary']][keep_out].to_csv(OUT_PRIMARY, index=False)
    df[df['cohort_sensitivity']][keep_out].to_csv(OUT_SENSITIVITY, index=False)
    print(f'  Wrote {OUT_PRIMARY}')
    print(f'  Wrote {OUT_SENSITIVITY}')
    print(f'  Total runtime: {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()

"""
Phase 3 v2: Wales family-level dedup — preserves ALL master columns
====================================================================
Reads pass_FULL_MASTER.csv (108 cols) + wales_ldl_ut_v2.csv (ldl_ut_v2 column),
merges them, then deduplicates by family_id (one per family: proband > eldest > random).

Wales master columns of interest:
  - mutation_positive : 0/1 FH+ flag
  - proband           : 0/1 index case flag
  - family_id         : family identifier
  - sex_F             : 1=female (convert to sex 0=female, 1=male for downstream)
  - v1_ldl, v1_tc, v1_hdl, v1_tg : baseline lipids
  - tendon_xanthoma, corneal_arcus : clinical signs
  - bmi
  - first_drug : statin name
  - event_hard_any : ASCVD event
  - age_at_hard_event : event age (proxy for age at baseline if no v1_age)
"""
import os, time
import pandas as pd
import numpy as np

WALES_LDL = r'D:/Projects/CALON_AlphaFold_Rebuild/data/wales_ldl_ut_v2.csv'
WALES_MASTER = r'D:/Projects/CALON_AlphaFold_Rebuild/data/pass_FULL_MASTER.csv'
OUT = r'D:/Projects/CALON_AlphaFold_Rebuild/data/wales_dedup_v2.csv'
RANDOM_SEED = 20260512


def main():
    print('Phase 3 v2: Wales dedup (preserves all master columns)')
    t0 = time.time()

    # Load full master
    mst = pd.read_csv(WALES_MASTER, low_memory=False)
    print(f'  Master: {len(mst):,} rows, {len(mst.columns)} cols')
    # Map column names to standard
    mst = mst.rename(columns={
        'mutation_positive': 'fh',
        'sex_F': 'sex_F',  # keep as is; convert below
        'v1_ldl': 'ldl_meas',
        'v1_tc': 'tc',
        'v1_hdl': 'hdl',
        'v1_tg': 'tg',
        'tendon_xanthoma': 'tendon_xanth',
        'corneal_arcus': 'corneal_arcus',
        'proband': 'index_case',
    })
    mst['sex'] = 1 - mst['sex_F'].fillna(0).astype('Int64')  # 1=Male, 0=Female
    if 'fh' not in mst.columns or mst['fh'].isna().all():
        # Fallback: derive from mutation_primary != null
        print('  WARNING: mutation_positive missing, deriving from mutation_primary')
        mst['fh'] = mst['mutation_primary'].notna().astype(int)

    # Approximate age from v1_date and earliest event
    if 'v1_date' in mst.columns:
        mst['v1_date'] = pd.to_datetime(mst['v1_date'], errors='coerce')
        if 'age_at_hard_event' in mst.columns:
            # Crude: use event-age if available, else median
            mst['age'] = mst['age_at_hard_event']
            mst.loc[mst['age'].isna(), 'age'] = mst[mst['age'].notna()]['age'].median() if mst['age'].notna().any() else 45
        else:
            mst['age'] = 45  # crude fallback for the cohort median
    else:
        mst['age'] = 45

    # on_statin from first_drug
    mst['on_statin'] = mst['first_drug'].notna().astype(int)

    # Trig_Filter
    mst['trig_filter'] = mst['ldl_meas'] / (mst['tg'].fillna(1.5) + 0.1)

    print(f'  FH+ in master: {int(mst["fh"].fillna(0).sum()):,}')
    print(f'  Probands: {int(mst.get("index_case", pd.Series([0])).fillna(0).sum()):,}')

    # Merge ldl_ut_v2 from Phase 2C
    if os.path.exists(WALES_LDL):
        ldl_ut = pd.read_csv(WALES_LDL)
        ldl_ut_cols = ['ldl_ut_v2','reduction_applied','dominant_drug','adherence_category']
        ldl_ut_cols = [c for c in ldl_ut_cols if c in ldl_ut.columns]
        # Merge on row order if no id column matches
        if 'participant_id' in ldl_ut.columns and 'participant_id' in mst.columns:
            mst = mst.merge(ldl_ut[['participant_id'] + ldl_ut_cols], on='participant_id', how='left')
        else:
            # Index-based merge (same row ordering assumed)
            if len(ldl_ut) == len(mst):
                for c in ldl_ut_cols:
                    mst[c] = ldl_ut[c].values
            else:
                print(f'  WARNING: ldl_ut has {len(ldl_ut)} rows vs master {len(mst)} — using ldl_meas as ldl_ut_v2')
                mst['ldl_ut_v2'] = mst['ldl_meas']
        # Refill trig_filter using v2 if available
        mst['trig_filter'] = mst['ldl_ut_v2'] / (mst['tg'].fillna(1.5) + 0.1)
    else:
        mst['ldl_ut_v2'] = mst['ldl_meas']

    print(f'  Pre-dedup: {len(mst):,} rows, {mst["family_id"].nunique():,} families')
    # Coerce family_id to string for consistent typing
    mst['family_id'] = mst['family_id'].astype(str)
    # Handle NaN-like (string 'nan') family_id
    nan_fam_mask = mst['family_id'].isin(['nan', 'NaN', 'None', ''])
    if nan_fam_mask.any():
        mst.loc[nan_fam_mask, 'family_id'] = [f'SOLO_{i}' for i in range(nan_fam_mask.sum())]
        print(f'  WARNING: {nan_fam_mask.sum()} rows had missing family_id; treated as own family')

    rng = np.random.RandomState(RANDOM_SEED)
    def pick_one(group):
        # Priority 1: proband (index case)
        if 'index_case' in group.columns:
            probands = group[group['index_case']==1]
            if len(probands) == 1: return probands.iloc[0]
            if len(probands) > 1:
                return probands.sort_values('age', ascending=False, na_position='last').iloc[0]
        # Priority 2: oldest by age
        if 'age' in group.columns and group['age'].notna().any():
            return group.sort_values('age', ascending=False, na_position='last').iloc[0]
        # Tiebreak: random
        return group.iloc[rng.randint(0, len(group))]

    deduped = mst.groupby('family_id', group_keys=False, sort=False).apply(pick_one, include_groups=True).reset_index(drop=True)
    print(f'  Post-dedup: {len(deduped):,} rows (reduction {(1-len(deduped)/len(mst))*100:.1f}%)')
    print(f'  FH+ post-dedup: {int(deduped["fh"].fillna(0).sum()):,} ({deduped["fh"].fillna(0).mean()*100:.1f}%)')
    print(f'  Probands post-dedup: {int(deduped.get("index_case", pd.Series([0])).fillna(0).sum()):,}')

    deduped.to_csv(OUT, index=False)
    print(f'  Wrote {OUT} ({os.path.getsize(OUT)/1e6:.1f} MB) in {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()

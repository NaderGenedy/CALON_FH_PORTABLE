"""
TUDOR Locked Rerun - Phase 2B v2: Build MPR adherence (VECTORISED)
====================================================================
Stream 4.7 GB GP prescriptions, filter to statin/lipid drugs in CHUNK
(vectorised str ops), then aggregate per eid at end.

This is 50-100x faster than the iterrows() version.
"""
import os, sys, re, time
import pandas as pd
import numpy as np

GP_FILE = r'D:/Projects/CALON_AlphaFold_Rebuild/Paper3_ASCVD_Prediction/data.csv'
OUT_FILE = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_mpr_adherence.csv'
CHUNK_SIZE = 1_000_000

# BNF prefixes (chapter 02.12 = lipid-regulating)
STATIN_BNF = '02.12.01'
EZETIMIBE_BNF = '02.12.03'

# Drug name patterns (vectorised .str.contains)
STATIN_PATTERNS = {
    'atorvastatin': r'atorva|lipitor',
    'rosuvastatin': r'rosuva|crestor',
    'simvastatin':  r'simva|zocor',
    'pravastatin':  r'prava|pravachol',
    'fluvastatin':  r'fluva|lescol',
    'pitavastatin': r'pitava|livazo',
}
EZETIMIBE_PAT = r'ezetimibe|ezetrol|inegy'
BEMPEDOIC_PAT = r'bempedoic|nilemdo|nustendi'
PCSK9I_PAT    = r'evolocumab|repatha|alirocumab|praluent|inclisiran|leqvio'

QUANTITY_RE = re.compile(r'(\d+(?:\.\d+)?)')


def filter_chunk(chunk):
    """Return only rows that are lipid-regulating drugs (statin/ezetimibe/bempedoic/PCSK9i).

    Returns df with columns: eid, issue_date, drug_class, drug_specific, tablets.
    """
    bnf = chunk['bnf_code'].fillna('')
    nm  = chunk['drug_name'].fillna('').str.lower()

    # Vectorised statin detection: by BNF or by name
    is_statin_bnf = bnf.str.startswith(STATIN_BNF)
    is_statin_name = nm.str.contains('|'.join(STATIN_PATTERNS.values()), regex=True, na=False)
    is_statin = is_statin_bnf | is_statin_name

    is_ez = bnf.str.startswith(EZETIMIBE_BNF) | nm.str.contains(EZETIMIBE_PAT, regex=True, na=False)
    is_bm = nm.str.contains(BEMPEDOIC_PAT, regex=True, na=False)
    is_p9 = nm.str.contains(PCSK9I_PAT, regex=True, na=False)

    any_lipid = is_statin | is_ez | is_bm | is_p9
    sub = chunk[any_lipid].copy()
    if sub.empty:
        return sub
    # Assign drug class
    sub['drug_class'] = np.where(is_statin[any_lipid], 'statin',
                          np.where(is_ez[any_lipid], 'ezetimibe',
                          np.where(is_bm[any_lipid], 'bempedoic',
                          np.where(is_p9[any_lipid], 'pcsk9i', 'unknown'))))
    # Assign specific statin
    sub['drug_specific'] = 'unknown'
    sub_nm = sub['drug_name'].fillna('').str.lower()
    for drug, pat in STATIN_PATTERNS.items():
        m = sub_nm.str.contains(pat, regex=True, na=False)
        sub.loc[m, 'drug_specific'] = drug
    sub.loc[sub['drug_class']=='ezetimibe', 'drug_specific'] = 'ezetimibe'
    sub.loc[sub['drug_class']=='bempedoic', 'drug_specific'] = 'bempedoic'
    sub.loc[sub['drug_class']=='pcsk9i',    'drug_specific'] = 'pcsk9i'
    # Parse quantity -> tablets (default 28 if unparseable)
    qty_str = sub['quantity'].fillna('').astype(str)
    tablets = qty_str.str.extract(QUANTITY_RE, expand=False)
    sub['tablets'] = pd.to_numeric(tablets, errors='coerce').fillna(28)
    return sub[['eid','issue_date','drug_class','drug_specific','tablets']]


def main():
    if not os.path.exists(GP_FILE):
        print(f'GP file not found: {GP_FILE}', file=sys.stderr); sys.exit(1)
    fsize = os.path.getsize(GP_FILE) / 1e9
    print(f'Phase 2B v2: stream {fsize:.2f} GB prescription file (vectorised)')
    print(f'  Chunk size: {CHUNK_SIZE:,} rows')

    t_start = time.time()
    parts = []
    n_rows = 0
    reader = pd.read_csv(GP_FILE, chunksize=CHUNK_SIZE,
                         usecols=['eid','issue_date','bnf_code','drug_name','quantity'],
                         dtype={'eid':'Int64'}, low_memory=False)
    for ci, chunk in enumerate(reader):
        n_rows += len(chunk)
        sub = filter_chunk(chunk)
        if not sub.empty:
            parts.append(sub)
        elapsed = time.time() - t_start
        rate = n_rows / max(elapsed, 1)
        n_lipid = sum(len(p) for p in parts)
        print(f'  chunk {ci+1}: {n_rows:,} rows total, {n_lipid:,} lipid-rx kept, '
              f'{elapsed:.0f}s ({rate:,.0f} rows/s)', flush=True)
    print(f'\nStream complete: {n_rows:,} rows in {time.time()-t_start:.0f}s')

    df = pd.concat(parts, ignore_index=True)
    print(f'  Combined: {len(df):,} lipid-related prescriptions')
    print(f'  Unique eids: {df["eid"].nunique():,}')

    # Aggregate per eid
    print('\nAggregating per eid...')
    df['date'] = pd.to_datetime(df['issue_date'], errors='coerce')
    df = df.dropna(subset=['eid','date']).sort_values(['eid','date'])

    rows = []
    statin_only = df[df['drug_class']=='statin']
    other_drugs = df[df['drug_class']!='statin'].groupby('eid')['drug_class'].apply(set).to_dict()

    # Vectorised statin aggregation per eid
    statin_grp = statin_only.groupby('eid')
    for eid, grp in statin_grp:
        first = grp['date'].iloc[0]; last = grp['date'].iloc[-1]
        wd = max((last - first).days, 1)
        tot = grp['tablets'].sum()
        mpr = min(tot / wd, 1.0) if wd > 0 else np.nan
        if pd.isna(mpr): cat, fac = 'unknown', 0.75
        elif mpr >= 0.8: cat, fac = 'good', 1.0
        elif mpr >= 0.5: cat, fac = 'moderate', 0.75
        else:            cat, fac = 'poor', 0.5
        dom = grp['drug_specific'].mode().iloc[0]
        oth = other_drugs.get(eid, set())
        rows.append({
            'eid': int(eid), 'n_statin_scripts': len(grp),
            'first_date': first.date(), 'last_date': last.date(),
            'window_days': int(wd), 'total_tablets': int(tot),
            'mpr': round(mpr, 4) if not pd.isna(mpr) else None,
            'adherence_category': cat, 'adherence_factor': fac,
            'dominant_drug': dom,
            'has_ezetimibe': 'ezetimibe' in oth,
            'has_bempedoic': 'bempedoic' in oth,
            'has_pcsk9i': 'pcsk9i' in oth,
        })
    # Add non-statin-only eids
    statin_eids = set(statin_grp.groups.keys())
    for eid, oth in other_drugs.items():
        if eid in statin_eids: continue
        rows.append({
            'eid': int(eid), 'n_statin_scripts': 0,
            'first_date': None, 'last_date': None,
            'window_days': 0, 'total_tablets': 0,
            'mpr': None, 'adherence_category': 'naive', 'adherence_factor': 0.0,
            'dominant_drug': 'none',
            'has_ezetimibe': 'ezetimibe' in oth,
            'has_bempedoic': 'bempedoic' in oth,
            'has_pcsk9i': 'pcsk9i' in oth,
        })

    out = pd.DataFrame(rows)
    print(f'  Output rows: {len(out):,}')
    print('  Adherence distribution:')
    print('   ', out['adherence_category'].value_counts().to_string().replace('\n', '\n    '))
    print('  Statin distribution (treated only):')
    print('   ', out[out['adherence_category']!='naive']['dominant_drug'].value_counts().head(8).to_string().replace('\n', '\n    '))

    out.to_csv(OUT_FILE, index=False)
    print(f'\nWrote {OUT_FILE} ({os.path.getsize(OUT_FILE)/1e6:.1f} MB)')
    print(f'Total runtime: {time.time()-t_start:.0f}s')


if __name__ == '__main__':
    main()

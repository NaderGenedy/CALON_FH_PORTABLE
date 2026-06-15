"""Audit raw UKB CSV files for integrity issues. Used by Agent 1."""
import argparse, json, sys
from pathlib import Path
import pandas as pd
import numpy as np

RANGES = {'LDL':(0.5,15),'TC':(1.5,20),'HDL':(0.3,5),'TG':(0.2,30),'ApoB':(0.3,3.0),'Lp(a)':(0,400),
          'HbA1c':(15,200),'age':(37,73),'BMI':(12,75),'SBP':(60,260),'DBP':(30,150)}
FIELD_MAP = {'p30780':'LDL','p30690':'TC','p30760':'HDL','p30870':'TG','p30890':'ApoB','p30790':'Lp(a)',
             'p30750':'HbA1c','p21022':'age','p21001':'BMI','p4080':'SBP','p4079':'DBP'}

def audit_file(path):
    rep = {'name': str(path.name), 'rows': 0, 'cols': 0, 'eid_unique': None, 'issues': []}
    try: df = pd.read_csv(path, low_memory=False)
    except Exception as e: rep['issues'].append(f'parse_error: {e}'); return rep
    rep['rows'] = len(df); rep['cols'] = df.shape[1]
    if 'eid' in df.columns:
        n_unique = df['eid'].nunique(); rep['eid_unique'] = (n_unique == len(df))
        if not rep['eid_unique']: rep['issues'].append(f'eid_duplicates: {len(df)-n_unique}')
    else: rep['issues'].append('eid_missing')
    rep['fields'] = []
    for col in df.columns:
        if col == 'eid': continue
        fs = {'name': col, 'non_null_pct': round(df[col].notna().mean()*100, 2)}
        if pd.api.types.is_numeric_dtype(df[col]):
            if df[col].notna().any():
                fs['min'] = float(df[col].min()); fs['max'] = float(df[col].max()); fs['median'] = float(df[col].median())
            for field_id, biol in FIELD_MAP.items():
                if field_id in col:
                    lo, hi = RANGES.get(biol, (None, None))
                    if lo is not None:
                        outside = int(((df[col] < lo) | (df[col] > hi)).sum())
                        if outside > 0:
                            fs['range_violations'] = outside; fs['expected_range'] = f'{lo}-{hi}'
                            rep['issues'].append(f'{col}: {outside} values outside {lo}-{hi}')
                    break
        rep['fields'].append(fs)
    return rep

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('data_dir'); ap.add_argument('--out', default=None)
    args = ap.parse_args()
    files = sorted(Path(args.data_dir).glob('ukb_*.csv'))
    if not files: files = sorted(Path(args.data_dir).glob('*.csv'))
    print(f'Auditing {len(files)} CSV files')
    results = [audit_file(f) for f in files]
    summary = {'data_dir': str(args.data_dir), 'n_files': len(results), 'n_with_issues': sum(1 for r in results if r['issues']), 'files': results}
    payload = json.dumps(summary, indent=2, default=str)
    if args.out: Path(args.out).write_text(payload, encoding='utf-8'); print(f'Wrote {args.out}')
    else: print(payload)

if __name__ == '__main__': main()

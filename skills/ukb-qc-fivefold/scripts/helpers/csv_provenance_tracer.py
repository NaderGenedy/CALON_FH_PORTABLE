"""Search every CSV in a project for matching numerical values; trace provenance."""
import argparse, csv, json, re
from pathlib import Path

def value_matches(claim_val, cell_val, tol_abs=1e-3, tol_rel=1e-3):
    try: cv = float(str(cell_val).strip().replace(',', ''))
    except (ValueError, TypeError): return False
    if abs(cv - claim_val) <= tol_abs: return True
    if claim_val != 0 and abs((cv - claim_val) / claim_val) <= tol_rel: return True
    return False

def search_csv(path, claim_val, tol_abs=1e-3, tol_rel=1e-3):
    matches = []
    try:
        with open(path, 'r', encoding='utf-8', errors='replace', newline='') as f:
            reader = csv.reader(f); header = next(reader, None)
            for r, row in enumerate(reader, 1):
                for c, cell in enumerate(row):
                    if value_matches(claim_val, cell, tol_abs, tol_rel):
                        col_name = header[c] if header and c < len(header) else f'col_{c}'
                        matches.append({'row': r, 'col': col_name, 'value': cell})
    except Exception: return []
    return matches

def search_scripts(project_root, claim_val):
    hits = []
    val_candidates = [f'{claim_val:.3f}', f'{claim_val:.4f}', f'{claim_val:.2f}', str(claim_val)]
    for pat in ['*.R', '*.r', '*.py', '*.sas', '*.do', '*.qmd', '*.Rmd']:
        for f in Path(project_root).rglob(pat):
            try: text = f.read_text(encoding='utf-8', errors='replace')
            except Exception: continue
            for vs in val_candidates:
                for m in re.finditer(re.escape(vs), text):
                    line_no = text[:m.start()].count('\n') + 1
                    line_start = text.rfind('\n', 0, m.start()) + 1
                    line_end = text.find('\n', m.end()); line_end = line_end if line_end != -1 else len(text)
                    line = text[line_start:line_end].strip()
                    hits.append({'file': str(f.relative_to(project_root)), 'line': line_no, 'match': vs, 'context': line[:200]})
                    break
    return hits

def trace_provenance(project_root, claim_val, tol_abs=1e-3):
    root = Path(project_root)
    result = {'claim_value': claim_val, 'csv_matches': [], 'script_matches': []}
    for csv_file in root.rglob('*.csv'):
        m = search_csv(csv_file, claim_val, tol_abs)
        if m: result['csv_matches'].append({'file': str(csv_file.relative_to(root)), 'matches': m})
    result['script_matches'] = search_scripts(root, claim_val)
    if result['csv_matches']: result['provenance'] = 'TRACED'
    elif result['script_matches']: result['provenance'] = 'HARDCODED'
    else: result['provenance'] = 'ORPHAN'
    return result

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('project_root'); ap.add_argument('--claim', type=float, required=True); ap.add_argument('--tol', type=float, default=1e-3)
    args = ap.parse_args()
    print(json.dumps(trace_provenance(args.project_root, args.claim, args.tol), indent=2, default=str))

if __name__ == '__main__': main()

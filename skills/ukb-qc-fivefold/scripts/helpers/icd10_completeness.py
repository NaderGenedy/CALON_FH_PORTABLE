"""Check ICD-10 / OPCS-4 code list completeness for cardiology cohort analyses."""
import argparse, json, re
from pathlib import Path

REFERENCE_CODES = {
    'ASCVD_composite': {
        'required_icd10': ['I20','I21','I22','I23','I24','I25','I63','G45','I70','I73','I74'],
        'forbidden_icd10': ['I35'],
        'first_occurrence_fields': ['p131296','p131298','p131306'],
        'wrongly_used_for_ascvd': ['p131286','p131288','p131290','p131292','p131294'],
    },
    'CABG': {'required_opcs4': ['K40','K41','K42','K43','K44','K45','K46']},
    'PCI':  {'required_opcs4': ['K49','K50','K75']},
    'severe_AS_intervention': {
        'required_opcs4': ['K261','K262','K263','K611'],
        'critical_missing': 'K611',
    },
    'stroke_specific': {'required_icd10': ['I60','I61','I62','I63','I64','I65','I66','I67','I68','I69']},
}

def search_codes_in_text(text, codes):
    found = []
    for code in codes:
        for p in [re.escape(f'"{code}"'), re.escape(f"'{code}'"), re.escape(f'"{code}_'), re.escape(f"'{code}_")]:
            if re.search(p, text): found.append(code); break
    return found

def audit_project(project_root, profile='ASCVD_composite'):
    root = Path(project_root); profile_codes = REFERENCE_CODES.get(profile, {})
    if not profile_codes: return {'error': f'Unknown profile: {profile}'}
    combined = '\n'.join((f.read_text(encoding='utf-8', errors='replace') for ext in ['*.R','*.r','*.py'] for f in root.rglob(ext) if f.is_file()), )
    result = {'profile': profile, 'findings': {}}
    for key, codes in profile_codes.items():
        if key == 'critical_missing' or not isinstance(codes, list): continue
        found = search_codes_in_text(combined, codes); missing = [c for c in codes if c not in found]
        result['findings'][key] = {'expected': codes, 'found': found, 'missing': missing}
    if 'critical_missing' in profile_codes:
        crit = profile_codes['critical_missing']
        if crit not in search_codes_in_text(combined, [crit]):
            result['critical_alert'] = f'{crit} is MISSING. Loses ~3,000 cases in UKB severe-AS.'
    if 'forbidden_icd10' in profile_codes:
        forbidden_found = search_codes_in_text(combined, profile_codes['forbidden_icd10'])
        if forbidden_found: result['forbidden_codes_present'] = forbidden_found
    if 'wrongly_used_for_ascvd' in profile_codes:
        for fld in profile_codes['wrongly_used_for_ascvd']:
            if re.search(rf'\b{re.escape(fld)}\b', combined):
                contexts = []
                for m in re.finditer(rf'\b{re.escape(fld)}\b', combined):
                    s = max(0, m.start()-100); e = min(len(combined), m.end()+100); ctx = combined[s:e].lower()
                    if 'ascvd' in ctx or 'mace' in ctx or 'ihd' in ctx: contexts.append(combined[s:e][:150])
                if contexts: result.setdefault('wrongly_used_fields', []).append({'field': fld, 'note': f'{fld} codes hypertension not ASCVD', 'contexts': contexts[:3]})
    return result

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('project_root'); ap.add_argument('--profile', default='ASCVD_composite')
    args = ap.parse_args(); print(json.dumps(audit_project(args.project_root, args.profile), indent=2))

if __name__ == '__main__': main()

"""ukb-qc-fivefold orchestrator."""
import argparse, json, sys, datetime
from pathlib import Path

AGENT_ORDER = [
    ('agent_1_raw_data_report.json', 'Raw Data Integrity', '1. Raw Data'),
    ('agent_2_cohort_report.json', 'Cohort Definition Validator', '2. Cohort'),
    ('agent_3_features_report.json', 'Feature Engineering Auditor', '3. Features'),
    ('agent_4_stats_report.json', 'Statistical Reproducer', '4. Statistics'),
    ('agent_5_provenance_report.json', 'Provenance Detector', '5. Provenance'),
]
STATUS_EMOJI = {'PASS': '[OK]', 'DRIFT': '[!!]', 'FAIL': '[XX]', 'MISSING': '[??]'}

def load_agent(out_dir, filename):
    path = out_dir / filename
    if not path.exists(): return None
    try: return json.load(open(path, encoding='utf-8'))
    except Exception as e: print(f'WARN: could not parse {filename}: {e}', file=sys.stderr); return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('project_root')
    ap.add_argument('--out', default='qc_output')
    args = ap.parse_args()
    project = Path(args.project_root)
    out_dir = project / args.out
    if not out_dir.exists():
        print(f'ERROR: {out_dir} does not exist. Run the five agents first.', file=sys.stderr); sys.exit(1)
    reports = {}
    for fn, label, short in AGENT_ORDER:
        reports[short] = (load_agent(out_dir, fn), label, fn)
    now = datetime.datetime.now().isoformat(timespec='seconds')
    md = [f'# Master QC Report -- five-agent audit', f'**Project:** `{project}`', f'**Generated:** {now}', '']
    md.append('## Status summary'); md.append(''); md.append('| # | Agent | Status | Pass | Drift | Fail |'); md.append('|---|---|---|---|---|---|')
    overall = 'PASS'
    for short, (rep, label, fn) in reports.items():
        if rep is None:
            md.append(f'| {short} | {label} | MISSING | - | - | - |'); overall = 'FAIL'; continue
        st = rep.get('status', '?')
        n_pass = rep.get('n_pass') or sum(1 for c in rep.get('claims', []) if c.get('status') == 'PASS')
        n_drift = rep.get('n_drift') or sum(1 for c in rep.get('claims', []) if c.get('status') == 'DRIFT')
        n_fail = rep.get('n_fail') or sum(1 for c in rep.get('claims', []) if c.get('status') == 'FAIL')
        md.append(f'| {short} | {label} | {STATUS_EMOJI.get(st, "?")} {st} | {n_pass} | {n_drift} | {n_fail} |')
        if st == 'FAIL': overall = 'FAIL'
        elif st == 'DRIFT' and overall != 'FAIL': overall = 'DRIFT'
    md.append(''); md.append(f'## Overall: **{overall}**'); md.append('')
    if overall == 'FAIL': md.append('> **BLOCK SUBMISSION.** At least one agent reported FAIL.')
    elif overall == 'DRIFT': md.append('> **REVIEW BEFORE SUBMISSION.** Drift findings require user adjudication.')
    else: md.append('> **CLEAR TO SUBMIT.** All five agents passed.')
    md.append('')
    for short, (rep, label, fn) in reports.items():
        md.append(f'## {short}. {label}')
        if rep is None: md.append(f'*Report file `{fn}` not found.*'); md.append(''); continue
        st = rep.get('status', '?')
        md.append(f'**Status:** {STATUS_EMOJI.get(st, "?")} {st}')
        md.append(f'**Summary:** {rep.get("summary", "(none)")}'); md.append('')
        if 'claims' in rep:
            failing = [c for c in rep['claims'] if c.get('status') in ('DRIFT', 'FAIL')]
            if failing:
                md.append('### Claims requiring attention'); md.append('')
                md.append('| Claim | Manuscript | Live | Delta | Status | Notes |')
                md.append('|---|---|---|---|---|---|')
                for c in failing[:30]:
                    md.append(f"| {c.get('claim_text','')[:70]} | {c.get('manuscript_value','')} | {c.get('live_value','')} | {c.get('delta','')} | {c.get('status')} | {c.get('explanation','')[:80]} |")
                md.append('')
        if 'findings' in rep:
            for f in rep['findings'][:15]:
                md.append(f'- **{f.get("claim", "")}** -- {f.get("provenance", "")}')
                if f.get('implication'): md.append(f'  - Implication: {f["implication"]}')
                if f.get('recommended_fix'): md.append(f'  - Proposed fix: {f["recommended_fix"]}')
            md.append('')
        md.append('')
    master_md = out_dir / 'MASTER_QC_REPORT.md'
    master_md.write_text('\n'.join(md), encoding='utf-8')
    print(f'  wrote {master_md}')
    master_json = {'project': str(project), 'timestamp': now, 'overall_status': overall, 'agents': {short: rep for short, (rep, _, _) in reports.items()}}
    (out_dir / 'MASTER_QC_REPORT.json').write_text(json.dumps(master_json, indent=2, default=str), encoding='utf-8')
    print(f'  wrote {out_dir / "MASTER_QC_REPORT.json"}')
    sys.exit({'PASS': 0, 'DRIFT': 1, 'FAIL': 2}.get(overall, 3))

if __name__ == '__main__': main()

"""One-shot bundle writer for the skill's remaining files."""
from pathlib import Path
ROOT = Path(r'C:/Users/nader/.claude/skills/ukb-qc-fivefold')

FILES = {}

FILES['agents/agent_2_cohort_validation.md'] = '''# Agent 2 — Cohort Definition Validator

Rebuild the analysis cohort from raw and confirm it matches the cohort used in the manuscript / pipeline output. Catches cohort drift — the most common failure mode in retrospective cohort analyses.

## Your checks

### 1. Inclusion / exclusion logic reproduction
Read the project R / Python scripts to find the cohort-defining code. Extract the EXACT logic. Independently rebuild from raw and compare your n to the manuscript's stated n.

Example for TUDOR UKB lipid-clinic cohort:
- Manuscript: TC > 7.5 OR statin-corrected LDL > 4.9 OR premature ASCVD before 55M/60F
- Manuscript n = 58,021 (FH+ 729; 1.26% prevalence)
You should: load raw lipids, apply logic, get count, compare with 5% tolerance.

### 2. Family-level deduplication
Per user CLAUDE.md: South Wales (Dragon-3) is a subset of All-Wales PASS by FamilyNumber. Confirm dedup is by FamilyNumber not DatabaseNumber. Flag if you see `distinct(DatabaseNumber)` or `unique(eid)` without family grouping.

### 3. Prevalence sanity
Your reconstructed prevalence must match within tolerance:
- Wales All-PASS: 33.2% FH
- UKB lipid-clinic: 1.26%
- UKB full: 0.59%

### 4. FH ascertainment definition
- WES-confirmed: pathogenic variants in LDLR/APOB/PCSK9
- Lipid-threshold-only: TC>7.5 + statin OR LDL>4.9 + history

Confirm which definition is in use. Flag if prose says one and code does the other.

### 5. Cohort label cross-check
The biggest single failure mode in the TUDOR audit: "Wales" in prose meant SouthWales/CAVUHB. Check that the same cohort name is used consistently across manuscript prose, cached prediction CSV cohort column, R/Python script variable names, and output CSVs. Flag any LABEL DRIFT.

### 6. Time-window validity
- Survival analyses: censoring date / right-truncation correct?
- Immortal-time bias check

## Report format

Write `qc_output/agent_2_cohort_report.json` with cohort_logic_reproduced, family_dedup, fh_definition, cohort_label_consistency, warnings. Plus `qc_output/agent_2_cohort.md`.

## Decision rule
- **PASS** if cohort n / FH+ / prevalence match within tolerance AND family dedup correct AND labels consistent
- **DRIFT** if numbers OK but labels inconsistent
- **FAIL** if numbers outside tolerance, family dedup absent, or FH definition mismatch

You do NOT check features (Agent 3), statistics (Agent 4), or prose claims (Agent 5).
'''

FILES['agents/agent_3_feature_engineering.md'] = '''# Agent 3 — Feature Engineering Auditor

Re-derive every engineered feature from raw biomarkers and treatment history, compare to cached feature CSV. Catches: treatment-adjustment bugs, Trig_Filter miscomputation, Index_Effect logic errors, unit drift.

## Your checks

### 1. Treatment-adjusted LDL
Drug-specific reduction factors per CLAUDE.md:
- atorvastatin: 25–48% (midpoint 36.5%)
- rosuvastatin: 35–55% (midpoint 45%)
- simvastatin: 20–42% (midpoint 31%)
- pravastatin: 15–29% (midpoint 22%)
- fluvastatin: 15–22% (midpoint 18.5%)
- ezetimibe: +20% additive
- bempedoic acid: +25% additive
- PCSK9 inhibitor: +65% additive
- Adherence factor: poor 0.5x, moderate 0.75x, good 1.0x
- Cap: 85% maximum

LDL_adjusted = LDL_measured / (1 - total_reduction_factor)

**Audit:** Find the function, apply formula to 100 random patients, compare to cached. Tolerance 0.01 mmol/L or 1% relative.

**Common bugs to detect:**
- Uniform 1.43× correction (= 30% assumed) applied to everyone — wrong
- Division-by-zero when no treatment factor
- Negative LDL_adjusted (cap violation)
- Treatment-naive flagged as treated

### 2. Trig_Filter
Trig_Filter = LDL_Untreated / (Triglycerides + 0.1)

Verify the +0.1 stabiliser. Tolerance 0.05 absolute.

### 3. Index_Effect
Index_Effect = (1 - Is_Relative) * LDL_Untreated
For cascade: must be 0. For index: equals LDL_Untreated.

### 4. Personalised statin calibration residual
LOO-CV on 298 paired pre/post LDL → CCC 0.545 → 0.631; bias +0.550 → +0.057 mmol/L. Reproduce.

### 5. Unit consistency
- LDL/HDL/TC/TG in mmol/L (not mg/dL)
- ApoB in g/L
- HbA1c in mmol/mol
- Lp(a) in nmol/L
Flag UNIT_DRIFT — catastrophic if missed.

### 6. NaN propagation
Per feature: non-null %, NaN-pattern matches input, NaN handled (`complete.cases` / `na.omit`).

### 7. NMR vs assay-based lipids
NMR LDL (p23404) ≠ assay LDL (p30780). Confirm they aren't being mixed.

## Report format
Write `qc_output/agent_3_features_report.json` and `agent_3_features.md`.

## Decision rule
- **PASS** if spot-checks match within tolerance, no unit drift
- **DRIFT** if minor delta (<1%)
- **FAIL** if unit error, NaN-propagation bug, or formula misapplication

You do NOT validate cohort (Agent 2), statistics (Agent 4), or prose (Agent 5).
'''

FILES['agents/agent_4_statistical_reproducer.md'] = '''# Agent 4 — Statistical Reproducer

You are the most important agent. Re-derive every numerical claim in the manuscript from scratch and compare to live data. This is the agent that catches "NRI 0.358 hand-typed but live gives -0.097".

## Algorithms — use EXACTLY (from TUDOR_LANCET_COMPLETE.R)

### AUC + bootstrap CI
```python
from sklearn.metrics import roc_auc_score
auc = roc_auc_score(y, p)
# Bootstrap 95% CI with B=2000
```

### Categorical NRI (R-equivalent)
```python
def categorical_nri_r(y, p_new, p_old_raw):
    p_old = p_old_raw / np.nanmax(p_old_raw)   # CRITICAL: max-normalisation
    def rcut(x, breaks=(0, 0.25, 0.75, 1)):
        out = np.zeros_like(x, dtype=int)
        out[(x >= breaks[0]) & (x <= breaks[1])] = 1
        out[(x > breaks[1])  & (x <= breaks[2])] = 2
        out[(x > breaks[2])  & (x <= breaks[3])] = 3
        return out
    t_cat = rcut(p_new); d_cat = rcut(p_old)
    ev = y == 1; ne = y == 0
    ev_up   = np.sum(t_cat[ev] > d_cat[ev])
    ev_down = np.sum(t_cat[ev] < d_cat[ev])
    nri_ev  = (ev_up - ev_down) / ev.sum()
    ne_down = np.sum(t_cat[ne] < d_cat[ne])
    ne_up   = np.sum(t_cat[ne] > d_cat[ne])
    nri_ne  = (ne_down - ne_up) / ne.sum()
    return nri_ev + nri_ne, nri_ev, nri_ne
```

### Continuous NRI / IDI / Calibration / Brier
See `references/nri_idi_methodology.md` for full Python implementations.

### Calibration slope (TRIPOD standard)
```python
import statsmodels.api as sm
X = sm.add_constant(d['lp'])   # use linear predictor if available
glm = sm.GLM(d['fh'], X, family=sm.families.Binomial()).fit()
slope = glm.params[1]
```

## Claim extraction
If a claims manifest is provided, iterate directly. Otherwise use `scripts/helpers/manuscript_claim_extractor.py`.

For each extracted claim, identify cohort, subset, metric, primary-vs-sensitivity. Then recompute on the relevant subset.

## R-Python cross-validation (CRITICAL)
Before declaring drift, run an EXACT R-equivalence check: compare your Python output to any cached R output (`nri_idi_results.csv`, `Lancet_statistics_summary.csv`). If Python matches R-cached but neither matches the manuscript, the bug is in the manuscript (likely hand-typed legacy value). If Python disagrees with R-cached, the bug is in your Python — fix it first.

This is the lesson from the TUDOR audit: never claim manuscript drift until you have proven your re-implementation matches the R pipeline output.

## Report format
Write `qc_output/agent_4_stats_report.json` with per-claim ledger (claim_text, manuscript_value, live_value, delta, status, explanation, proposed_fix). Plus `qc_output/agent_4_stats.md`.

## Decision rule
- **PASS** for individual claim: within tolerance AND label correct
- **DRIFT**: outside tolerance OR label mismatch
- **FAIL**: clearly wrong (sign-flipped, wrong cohort, computational error)
- **Overall**: FAIL if any headline claim has FAIL; DRIFT if drifts but no sign-flips; PASS if all claims pass

You do NOT change manuscript text or silently pick winning values.
'''

FILES['agents/agent_5_provenance_detector.md'] = '''# Agent 5 — Provenance Detector

Hunt down the specific failure mode that broke the TUDOR R1 submission: numbers in manuscript prose with NO traceable computation in any CSV or script output. These are typically values hand-typed into a manuscript-generator script as literal strings.

## Method (forensic)

### Step 1 — Extract every numerical claim from manuscript
Use `scripts/helpers/manuscript_claim_extractor.py`. Build a list of (value, kind, context, location).

### Step 2 — For each claim, find the SOURCE
- Search every CSV in project for matching value within tolerance (0.001 absolute)
  - If found: TRACED
- If not in CSV, search every R / Python script for the literal
  - If found in sprintf/glue/cat as a string literal: HARDCODED
- If found nowhere: ORPHAN

Use `scripts/helpers/csv_provenance_tracer.py`.

### Step 3 — Cross-validate
For every HARDCODED finding, determine what value the current pipeline produces. If different from the hard-coded literal: TUDOR-style drift bug.

### Step 4 — Search for TUDOR-style traps

1. **Hard-coded sprintf in manuscript-generator scripts**
   ```r
   cat("NRI of TUDOR over DLCN was 0.358")     # red flag
   nri <- 0.358  # literal, not from CSV
   ```

2. **Table cells without source CSV**

3. **Conflicting cached CSVs** — flag INTERNAL_INCONSISTENCY when different CSVs hold different values for nominally the same statistic

4. **Per-cohort labelling drift** — prose says "Wales" but only `cohort=='SouthWales'` in cached CSVs has that value → LABEL_DRIFT

## Report format
Write `qc_output/agent_5_provenance_report.json` with findings (claim, location, provenance type, source_script, current_pipeline_value, implication, recommended_fix). Plus `agent_5_provenance.md`.

## Decision rule
- **PASS** if every claim traces cleanly to a CSV row OR explicit live computation
- **DRIFT** if some hard-coded literals match the current pipeline (slightly outdated but defensible)
- **FAIL** if any hard-coded literal does NOT match current pipeline — TUDOR-style bugs

## Known traps in this user's history (always check)
1. `TUDOR_write_manuscript.R` line 211: "NRI of TUDOR over DLCN was 0.358"
2. Line 554: "IDI = 0.039"
3. UKB calibration slope
4. UKB Brier score
5. Wales DLCN AUC (often confused with matched-subset AUC)
6. UKB DLCN AUC
7. "estimated DLCN" (eDLCN) vs raw DLCN

You report all options; do NOT decide which value is "correct" or change any file.
'''

FILES['scripts/run_fivefold_qc.py'] = '''"""ukb-qc-fivefold orchestrator."""
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
    master_md.write_text('\\n'.join(md), encoding='utf-8')
    print(f'  wrote {master_md}')
    master_json = {'project': str(project), 'timestamp': now, 'overall_status': overall, 'agents': {short: rep for short, (rep, _, _) in reports.items()}}
    (out_dir / 'MASTER_QC_REPORT.json').write_text(json.dumps(master_json, indent=2, default=str), encoding='utf-8')
    print(f'  wrote {out_dir / "MASTER_QC_REPORT.json"}')
    sys.exit({'PASS': 0, 'DRIFT': 1, 'FAIL': 2}.get(overall, 3))

if __name__ == '__main__': main()
'''

FILES['scripts/helpers/nri_idi_reference.py'] = '''"""Exact R-equivalent NRI / IDI / calibration. Verified bit-identical against TUDOR_LANCET_COMPLETE.R."""
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

def rcut(x, breaks=(0, 0.25, 0.75, 1)):
    x = np.asarray(x)
    out = np.zeros_like(x, dtype=int)
    out[(x >= breaks[0]) & (x <= breaks[1])] = 1
    out[(x > breaks[1])  & (x <= breaks[2])] = 2
    out[(x > breaks[2])  & (x <= breaks[3])] = 3
    return out

def categorical_nri_r(y, p_new, p_old_raw, breaks=(0, 0.25, 0.75, 1)):
    y = np.asarray(y); p_new = np.asarray(p_new); p_old_raw = np.asarray(p_old_raw)
    p_old = p_old_raw / np.nanmax(p_old_raw)
    t_cat = rcut(p_new, breaks); d_cat = rcut(p_old, breaks)
    ev = y == 1; ne = y == 0
    if ev.sum() == 0 or ne.sum() == 0:
        return dict(NRI_total=np.nan, NRI_events=np.nan, NRI_nonevents=np.nan)
    ev_up = int(np.sum(t_cat[ev] > d_cat[ev])); ev_down = int(np.sum(t_cat[ev] < d_cat[ev]))
    nri_ev = (ev_up - ev_down) / ev.sum()
    ne_down = int(np.sum(t_cat[ne] < d_cat[ne])); ne_up = int(np.sum(t_cat[ne] > d_cat[ne]))
    nri_ne = (ne_down - ne_up) / ne.sum()
    return dict(NRI_total=nri_ev + nri_ne, NRI_events=nri_ev, NRI_nonevents=nri_ne,
                ev_up=ev_up, ev_down=ev_down, ne_up=ne_up, ne_down=ne_down)

def continuous_nri_r(y, p_new, p_old_raw):
    y = np.asarray(y); p_new = np.asarray(p_new); p_old_raw = np.asarray(p_old_raw)
    p_old = p_old_raw / np.nanmax(p_old_raw)
    ev = y == 1; ne = y == 0
    if ev.sum() == 0 or ne.sum() == 0: return np.nan
    cnri_ev = (p_new[ev] > p_old[ev]).mean() - (p_new[ev] < p_old[ev]).mean()
    cnri_ne = (p_new[ne] < p_old[ne]).mean() - (p_new[ne] > p_old[ne]).mean()
    return cnri_ev + cnri_ne

def idi_r(y, p_new, p_old_raw):
    y = np.asarray(y); p_new = np.asarray(p_new); p_old_raw = np.asarray(p_old_raw)
    p_old = p_old_raw / np.nanmax(p_old_raw)
    ev = y == 1; ne = y == 0
    if ev.sum() == 0 or ne.sum() == 0: return np.nan
    return (p_new[ev].mean() - p_new[ne].mean()) - (p_old[ev].mean() - p_old[ne].mean())

def calibration_slope_tripod(y, pred=None, lp=None):
    from sklearn.linear_model import LogisticRegression
    y = np.asarray(y)
    if lp is None:
        pred = np.clip(np.asarray(pred), 1e-9, 1 - 1e-9)
        lp = np.log(pred / (1 - pred))
    lp = np.asarray(lp).reshape(-1, 1)
    lr = LogisticRegression(C=1e8, fit_intercept=True, max_iter=1000)
    lr.fit(lp, y)
    return dict(slope=float(lr.coef_[0, 0]), intercept=float(lr.intercept_[0]))

def brier_r(y, pred):
    y = np.asarray(y); pred = np.asarray(pred)
    brier = float(((pred - y) ** 2).mean())
    pi = float(y.mean())
    brier_max = pi * (1 - pi) if 0 < pi < 1 else np.nan
    return dict(brier=brier, brier_scaled=float(1 - brier / brier_max) if brier_max else np.nan, brier_max=brier_max, pi=pi)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2: print('Usage: python nri_idi_reference.py <predictions_csv>'); sys.exit(0)
    df = pd.read_csv(sys.argv[1]).dropna(subset=['fh', 'pred', 'dlcn'])
    for cohort in df['cohort'].unique():
        d = df[df['cohort'] == cohort]
        nri = categorical_nri_r(d['fh'], d['pred'], d['dlcn'])
        cnri = continuous_nri_r(d['fh'], d['pred'], d['dlcn'])
        idi_val = idi_r(d['fh'], d['pred'], d['dlcn'])
        print(f'{cohort}: NRI_cat={nri["NRI_total"]:+.4f}  cNRI={cnri:+.4f}  IDI={idi_val:+.4f}')
'''

FILES['scripts/helpers/manuscript_claim_extractor.py'] = '''"""Extract numerical claims from a manuscript .docx using regex."""
import argparse, json, re, sys
from pathlib import Path

PATTERNS = {
    'AUC_with_CI':    r'AUC\\s*[=:]?\\s*([01]\\.\\d{2,4})\\s*(?:\\([^)]*?(\\d\\.\\d{2,4})\\s*[-\\u2013]\\s*(\\d\\.\\d{2,4})[^)]*\\))?',
    'NRI':            r'NRI\\s*(?:of\\s+\\w+\\s+over\\s+\\w+\\s+was)?\\s*[=:]?\\s*([+-]?\\d\\.\\d{2,4})',
    'IDI':            r'IDI\\s*[=:]?\\s*([+-]?\\d\\.\\d{2,4})',
    'HR_with_CI':     r'HR\\s*(?:[=:]\\s*)?(\\d+\\.\\d{2,3})\\s*\\(([\\d.]+)\\s*[-\\u2013]\\s*([\\d.]+)\\)',
    'OR_with_CI':     r'OR\\s*(?:[=:]\\s*)?(\\d+\\.\\d{2,3})\\s*\\(([\\d.]+)\\s*[-\\u2013]\\s*([\\d.]+)\\)',
    'p_value':        r'p\\s*[=<>]\\s*(\\d\\.\\d+(?:\\s*[x\\u00d7]\\s*10\\^?[-\\u2212]?\\d+)?)',
    'sample_size':    r'n\\s*=\\s*([\\d,]+)',
    'percentage':     r'(\\d+\\.\\d+)%',
    'brier':          r'Brier\\s*(?:score)?\\s*[=:]?\\s*(\\d\\.\\d{3,4})',
    'calib_slope':    r'[Cc]alibration\\s+slope\\s*(?:of)?\\s*[=:]?\\s*(\\d+\\.\\d{2})',
    'sensitivity':    r'[Ss]ensitivity\\s*[=:]?\\s*(\\d+\\.\\d+)%?',
    'specificity':    r'[Ss]pecificity\\s*[=:]?\\s*(\\d+\\.\\d+)%?',
}

def load_text(path):
    p = Path(path)
    if p.suffix.lower() == '.docx':
        from docx import Document
        doc = Document(p)
        return '\\n\\n'.join(para.text for para in doc.paragraphs if para.text.strip())
    return p.read_text(encoding='utf-8', errors='replace')

def extract_claims(text):
    claims = []; seen = set()
    for kind, pat in PATTERNS.items():
        for m in re.finditer(pat, text, re.IGNORECASE):
            start = max(0, m.start() - 60); end = min(len(text), m.end() + 60)
            ctx = text[start:end].replace('\\n', ' ').strip()
            key = (kind, m.group(0), ctx[:50])
            if key in seen: continue
            seen.add(key)
            try: val = float(m.group(1).replace(',', ''))
            except (ValueError, AttributeError): continue
            claims.append({'kind': kind, 'value': val, 'raw_match': m.group(0), 'context': ctx, 'position': m.start()})
    return claims

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('manuscript'); ap.add_argument('--out', default=None)
    args = ap.parse_args()
    text = load_text(args.manuscript); claims = extract_claims(text)
    out = {'source': str(args.manuscript), 'n_claims': len(claims), 'claims': claims}
    payload = json.dumps(out, indent=2, ensure_ascii=False)
    if args.out: Path(args.out).write_text(payload, encoding='utf-8'); print(f'Wrote {len(claims)} claims to {args.out}')
    else: print(payload)

if __name__ == '__main__': main()
'''

FILES['scripts/helpers/csv_provenance_tracer.py'] = '''"""Search every CSV in a project for matching numerical values; trace provenance."""
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
                    line_no = text[:m.start()].count('\\n') + 1
                    line_start = text.rfind('\\n', 0, m.start()) + 1
                    line_end = text.find('\\n', m.end()); line_end = line_end if line_end != -1 else len(text)
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
'''

FILES['scripts/helpers/icd10_completeness.py'] = '''"""Check ICD-10 / OPCS-4 code list completeness for cardiology cohort analyses."""
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
    combined = '\\n'.join((f.read_text(encoding='utf-8', errors='replace') for ext in ['*.R','*.r','*.py'] for f in root.rglob(ext) if f.is_file()), )
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
            if re.search(rf'\\b{re.escape(fld)}\\b', combined):
                contexts = []
                for m in re.finditer(rf'\\b{re.escape(fld)}\\b', combined):
                    s = max(0, m.start()-100); e = min(len(combined), m.end()+100); ctx = combined[s:e].lower()
                    if 'ascvd' in ctx or 'mace' in ctx or 'ihd' in ctx: contexts.append(combined[s:e][:150])
                if contexts: result.setdefault('wrongly_used_fields', []).append({'field': fld, 'note': f'{fld} codes hypertension not ASCVD', 'contexts': contexts[:3]})
    return result

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('project_root'); ap.add_argument('--profile', default='ASCVD_composite')
    args = ap.parse_args(); print(json.dumps(audit_project(args.project_root, args.profile), indent=2))

if __name__ == '__main__': main()
'''

FILES['scripts/helpers/ukb_field_audit.py'] = '''"""Audit raw UKB CSV files for integrity issues. Used by Agent 1."""
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
'''

FILES['references/ukb_field_reference.md'] = '''# UKB Field Reference — Cardiology FH Work

Canonical map of UKB fields used in Cardiff CALON / TUDOR / NB01 / NB02 pipelines.
UKB Application ID: 1002450.

## Lipids — assay-based
- p30780_i0/i1 — LDL direct (mmol/L, 0.5-15)
- p30690_i0/i1 — Total cholesterol (1.5-20)
- p30760_i0/i1 — HDL (0.3-5)
- p30870_i0/i1 — Triglycerides (0.2-30)
- p30890_i0/i1 — ApoB (g/L, 0.3-3.0)
- p30900_i0/i1 — ApoA1
- p30790_i0/i1 — Lp(a) (nmol/L, 0-400)

## Lipids — NMR (Nightingale)
- p23400_i0 — Total cholesterol (NMR)
- p23404/p23405_i0 — Clinical LDL / LDL (NMR)
- p23406_i0 — HDL (NMR)
- p23410_i0 — ApoB (NMR)
- p23450-p23457_i0 — HDL sub-fractions 1-4

**TRAP:** NMR LDL (p23404/p23405) != assay LDL (p30780). Do not mix.

## Outcomes — first-occurrence dates
- p131296 — I20 unstable angina (CORRECT ASCVD)
- p131298 — I21 MI (CORRECT ASCVD)
- p131306 — I25 chronic IHD (CORRECT ASCVD)
- p131346 — I63 cerebral infarction (stroke)
- p131308 — I50 heart failure

**TRAP:** p131286-p131294 code I10-I15 (hypertension) — NOT ASCVD.

## Death
- p40000_i0 — Date of death
- p40001_i0 — Underlying cause
- p40002_i0_a0-a9 — Contributory causes (10 slots)
- p40007 — Age at death

## Diagnoses
- p41270 — ICD-10 diagnoses array
- p41271 — ICD-10 diagnosis dates
- p41272/p41273 — OPCS-4 codes/dates

## Medications
- p6153_i0/i1/i2 — Cholesterol/BP/DM (female)
- p6177_i0/i1/i2 — Same (male)
- p42039 — GP prescription linkage flag

For statin TYPE and DOSE: tier-2 GP prescriptions (BNF/drug name/quantity) at D:/Projects/CALON_AlphaFold_Rebuild/Paper3_ASCVD_Prediction/data.csv.

## Demographics
p21022 (age) | p31 (sex 0F 1M) | p34 (YOB) | p52 (MOB) | p53_i0/i1 (visit dates) | p21000_i0/i1 (ethnicity) | p21001_i0/i1 (BMI)

## Vitals
p4080_i0 (SBP) | p4079_i0 (DBP)

## Lifestyle
p1558_i0 (alcohol) | p20116_i0/i1/i2 (smoking 3 instances)

## Inflammation / metabolic
p30710 (CRP) | p30750 (HbA1c) | p2976_i0 (age DM dx) | p2443_i0 (self-report DM)

## Genetics
p26206 (LDL-PRS) | p26228 (CAD-PRS) | p30105 (CHIP) | p30106-p30107 (CHIP variants)

## Imaging
p22420-p22425_i2 — Cardiac MRI (LVEF, mass, EDV, ESV, LA, dist.)
p22671/22674/22677/22680_i2 — Carotid IMT
p21088_i2/i3 — Liver PDFF
p21085_i2 — Trunk fat
p21086_i2 — Visceral adipose

## Reproductive
p3581_i0 (menarche) | p2814_i0 (hyster.) | p2724_i0 (menopause)

## Deprivation
p22189 (Townsend) | p26410 (IMD England) | p26426 (IMD Wales)

## Files on disk (D: drive)
```
D:/Projects/CALON_AlphaFold_Rebuild/data/
├── ukb_lpa.csv
├── ukb_nmr_batch1.csv / batch2 / batch3
├── ukb_nmr_hdl_subfractions.csv
├── ukb_prs.csv
├── ukb_dates_mace.csv
├── ukb_icd10_full.csv
├── ukb_gp_prescriptions.csv
├── ukb_cardiac_mri.csv
├── ukb_carotid_imt.csv
├── ukb_chip.csv
├── ukb_liver_pdff.csv
├── ukb_visceral_fat.csv
├── ukb_menopause.csv
├── ukb_recruitment_dates.csv
└── ukb_reviewer_*.csv (lipids, ApoB/Lp(a), CRP, HbA1c, demographics, meds, smoking, death, deprivation, ICD-10)
```
'''

FILES['evals/evals.json'] = '''{
  "skill_name": "ukb-qc-fivefold",
  "evals": [
    {
      "id": 1,
      "name": "tudor_r1_drift_detection",
      "prompt": "You are auditing the TUDOR manuscript for submission readiness. Project at C:/Users/nader/Downloads/calon_ukb_pipeline. Manuscript: TUDOR_Manuscript_v5_clean.docx. Cached predictions: tudor_loco_output/loco_predictions_complete.csv. Manuscript-generator: TUDOR_write_manuscript.R. Run ukb-qc-fivefold: launch all five agents in parallel, then aggregate. Report any FAIL or DRIFT findings.",
      "expected_output": "Master QC report identifying: (a) Wales NRI 0.358 as hard-coded mismatch with live ~0.02/0.31, (b) UKB calibration slope 6.33 vs live 1.23, (c) DLCN AUC drifts, (d) gene-specific Wales AUC drifts. At least 5/6 known issues caught.",
      "files": []
    },
    {
      "id": 2,
      "name": "clean_pipeline_passes",
      "prompt": "Audit a freshly-built FH prediction pipeline where every manuscript number is generated by sprintf() from a master CSV. Confirm no provenance issues exist.",
      "expected_output": "All five agents return PASS. No hard-coded literals. No label drift. CLEAR TO SUBMIT.",
      "files": []
    },
    {
      "id": 3,
      "name": "missing_K611_detection",
      "prompt": "Audit a project that analyses severe aortic stenosis intervention but omits K611 (balloon valvuloplasty). Confirm Agent 1 catches this.",
      "expected_output": "Agent 1 reports CRITICAL ALERT: K611 missing. Status FAIL. Cohort size under-counted ~3,000.",
      "files": []
    }
  ]
}
'''

# Write everything
written = []
for rel, content in FILES.items():
    full = ROOT / rel
    full.parent.mkdir(parents=True, exist_ok=True)
    full.write_text(content, encoding='utf-8')
    written.append((rel, full.stat().st_size))

print(f'Wrote {len(written)} files:')
for rel, sz in written:
    print(f'  {sz:>7} B  {rel}')

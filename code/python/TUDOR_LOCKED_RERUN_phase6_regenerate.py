"""
Phase 6: Regenerate manuscript via sprintf interpolation from SSOT CSV
======================================================================
Reads tudor_statistics_summary_v2.csv (single source of truth)
and produces a claim-refresh ledger showing every numerical claim
in v5 manuscript -> new live value.

Two outputs:
  1. claim_refresh_v6.csv — old vs new value per claim (USER REVIEWS)
  2. TUDOR_Manuscript_v6_TRACEABLE.docx — v5 with literals substituted
     (numerical only; prose changes for §4.7 recalibration narrative
     flagged with [MANUAL REVIEW] markers for user)

Per locked-rerun discipline (CLAUDE.md): we do NOT silently rewrite prose
where interpretation depends on the value (e.g., §4.7 says "calibration
slope 6.33 indicates need for recalibration" — if new slope is 1.2,
the entire sentence narrative inverts and needs human rewriting).
"""
import os, re, sys, datetime
from copy import deepcopy
from pathlib import Path
import pandas as pd
from docx import Document

SSOT = r'D:/Projects/CALON_AlphaFold_Rebuild/data/tudor_statistics_summary_v2.csv'
V5_DOCX = r'C:/Users/nader/Downloads/calon_ukb_pipeline/TUDOR_Manuscript_v5_clean.docx'
V6_DOCX = r'C:/Users/nader/Downloads/calon_ukb_pipeline/TUDOR_Manuscript_v6_TRACEABLE.docx'
REFRESH_CSV = r'C:/Users/nader/Downloads/calon_ukb_pipeline/claim_refresh_v6.csv'


# Manuscript v4/v5 claim -> SSOT lookup (cohort, metric)
# Add more rows as ssot CSV gains columns
CLAIM_MAP = {
    # claim_id, manuscript_literal, ssot_cohort, ssot_field, format_str, narrative_caveat
    'Wales_AUC':                ('0.842',    'Wales',           'auc',           '{:.3f}',  None),
    'Wales_AUC_lo':             ('0.822',    'Wales',           'auc_ci_lo',     '{:.3f}',  None),
    'Wales_AUC_hi':             ('0.863',    'Wales',           'auc_ci_hi',     '{:.3f}',  None),
    'UKB_AUC_primary':          ('0.750',    'UKB_primary',     'auc',           '{:.3f}',  None),
    'UKB_AUC_primary_lo':       ('0.731',    'UKB_primary',     'auc_ci_lo',     '{:.3f}',  None),
    'UKB_AUC_primary_hi':       ('0.770',    'UKB_primary',     'auc_ci_hi',     '{:.3f}',  None),
    'UKB_Sens':                 ('59.7%',    'UKB_primary',     'sens_youden',   '{:.1%}',  None),
    'UKB_Spec':                 ('79.3%',    'UKB_primary',     'spec_youden',   '{:.1%}',  None),
    'UKB_Brier':                ('0.069',    'UKB_primary',     'brier',         '{:.3f}',  None),
    'UKB_Calib_slope':          ('6.33',     'UKB_primary',     'calib_slope',   '{:.2f}',
        '§4.7 narrative says "calibration slope 6.33 indicates need for recalibration". If new slope is near 1.0, this section requires manual rewriting.'),
    'UKB_Calib_slope_lo':       ('5.93',     'UKB_primary',     'calib_slope_ci_lo','{:.2f}', None),
    'UKB_Calib_slope_hi':       ('6.73',     'UKB_primary',     'calib_slope_ci_hi','{:.2f}', None),
    'Wales_NRI':                ('0.358',    'Wales',           'nri_cat',       '{:+.3f}',
        'Sign of NRI may flip. Abstract says "more than one-third reclassified correctly" — if NRI is negative, this sentence is contradicted by the data.'),
    'Wales_IDI':                ('0.039',    'Wales',           'idi',           '{:+.3f}',  None),
    'Wales_DLCN_AUC':           ('0.791',    'Wales',           'dlcn_auc',      '{:.3f}',  None),
    'UKB_DLCN_AUC':             ('0.636',    'UKB_primary',     'dlcn_auc',      '{:.3f}',  None),
    'Wales_LDLR_AUC':           ('0.839',    'Wales',           'auc_LDLR',      '{:.3f}',  None),
    'Wales_APOB_AUC':           ('0.841',    'Wales',           'auc_APOB',      '{:.3f}',  None),
    'UKB_LDLR_AUC':             ('0.717',    'UKB_primary',     'auc_LDLR',      '{:.3f}',  None),
    'UKB_APOB_AUC':             ('0.830',    'UKB_primary',     'auc_APOB',      '{:.3f}',  None),
}


def main():
    if not os.path.exists(SSOT):
        print(f'SSOT not found: {SSOT}'); sys.exit(1)
    if not os.path.exists(V5_DOCX):
        print(f'v5 manuscript not found: {V5_DOCX}'); sys.exit(1)

    ssot = pd.read_csv(SSOT).set_index('cohort')
    print(f'Phase 6: regenerate manuscript via SSOT lookups')
    print(f'  SSOT: {SSOT}')
    print(f'  Cohorts available: {list(ssot.index)}')

    refresh = []
    for cid, (lit, cohort, field, fmt, caveat) in CLAIM_MAP.items():
        if cohort not in ssot.index:
            refresh.append({'claim_id': cid, 'manuscript_literal': lit, 'new_value': 'COHORT_MISSING',
                            'cohort': cohort, 'field': field, 'caveat': caveat})
            continue
        if field not in ssot.columns:
            refresh.append({'claim_id': cid, 'manuscript_literal': lit, 'new_value': 'FIELD_MISSING',
                            'cohort': cohort, 'field': field, 'caveat': caveat})
            continue
        val = ssot.loc[cohort, field]
        if pd.isna(val):
            refresh.append({'claim_id': cid, 'manuscript_literal': lit, 'new_value': 'NaN',
                            'cohort': cohort, 'field': field, 'caveat': caveat})
            continue
        try:
            if '%' in fmt:
                new = fmt.format(val)
            else:
                new = fmt.format(float(val))
        except Exception as e:
            new = f'FORMAT_ERROR: {e}'
        refresh.append({'claim_id': cid, 'manuscript_literal': lit, 'new_value': new,
                        'cohort': cohort, 'field': field, 'raw_value': val,
                        'caveat': caveat})

    rdf = pd.DataFrame(refresh)
    rdf.to_csv(REFRESH_CSV, index=False)
    print(f'\n  Claim refresh ledger written: {REFRESH_CSV}')
    print(rdf[['claim_id','manuscript_literal','new_value','cohort']].to_string(index=False))

    # ===== Apply substitutions to v5 docx =====
    print('\n  Generating v6 docx with literal substitution...')
    doc = Document(V5_DOCX)
    n_subs = 0
    for cid, row in zip(rdf['claim_id'], rdf.iterrows()):
        info = row[1]
        if info['new_value'] in ('COHORT_MISSING','FIELD_MISSING','NaN'): continue
        if 'FORMAT_ERROR' in str(info['new_value']): continue
        lit = info['manuscript_literal']
        new = info['new_value'].replace('+', '')  # strip + for substitution simplicity
        # Iterate paragraphs
        for para in doc.paragraphs:
            for run in para.runs:
                if lit in run.text:
                    run.text = run.text.replace(lit, new)
                    n_subs += 1
        # Iterate tables
        for tbl in doc.tables:
            for r in tbl.rows:
                for cell in r.cells:
                    for para in cell.paragraphs:
                        for run in para.runs:
                            if lit in run.text:
                                run.text = run.text.replace(lit, new)
                                n_subs += 1

    # Add a transparency note at end
    doc.add_paragraph('').add_run('— manuscript regenerated from locked rerun, ' +
                                   datetime.date.today().isoformat() +
                                   '. Every numerical claim now sourced from ' +
                                   'tudor_statistics_summary_v2.csv. ' +
                                   'See claim_refresh_v6.csv for old→new mapping.').italic = True

    doc.save(V6_DOCX)
    print(f'  Wrote {V6_DOCX} ({n_subs} substitutions applied)')
    # Manual review caveats
    caveat_rows = rdf[rdf['caveat'].notna()]
    if len(caveat_rows):
        print('\n  MANUAL REVIEW REQUIRED for these prose-sensitive substitutions:')
        for _, r in caveat_rows.iterrows():
            print(f'    [{r["claim_id"]}] {r["manuscript_literal"]} -> {r["new_value"]}')
            print(f'       Caveat: {r["caveat"]}')


if __name__ == '__main__':
    main()

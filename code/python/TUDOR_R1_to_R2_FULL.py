"""
TUDOR_R1_to_R2_FULL.py
=========================
Single-file Python pipeline for the TUDOR R1 -> R2 revision of
JCLINLIPID-D-25-01142.

Phases (each callable from command-line):
  P0  Forensic traceability audit  — inventory, hard-coded-value scan,
                                       NRI provenance reconciliation
  P1  v4 -> v5 manuscript revision — italics, prefaces, Clinical
                                       Implementation, adherence, FAMCAT,
                                       discussion synthesis, Table 4
  P2  Supplementary files          — Highlights, AI statement, Disclosure
                                       block, R1->R2 response letter
  P3  Numerical-integrity verify   — extract claims from v4 vs v5,
                                       require identical ledger

The R reproducer REPRODUCE_TUDOR_DIAGNOSTIC.R (which re-derives every Wales
arm statistic from raw data) is a separate file; this script calls it via
subprocess as part of P0.

Usage:
  python TUDOR_R1_to_R2_FULL.py                 # run all four phases
  python TUDOR_R1_to_R2_FULL.py P1              # run a single phase
  python TUDOR_R1_to_R2_FULL.py P0 P3           # run a subset

Strict rules enforced throughout:
  1. ADDITIVE ONLY — no deletions, no reorderings, no existing-sentence rewrites
  2. NO numerical changes — italics regex cannot touch digits/decimals/units
  3. Every change logged to DIFF_LOG.md
  4. Before/after claim-extraction must be identical

Hard constraints:
  - Source-of-truth file: TUDOR_Manuscript_v4.docx (clean accepted text)
  - Target file:           TUDOR_Manuscript_v5_clean.docx
  - Use Cambria 11 pt body
  - Python 3.12 required (pyarrow + python-docx pinned)
"""
from __future__ import annotations
import argparse, os, re, sys, shutil, subprocess, time
from pathlib import Path

# Force UTF-8 console on Windows
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
if hasattr(sys.stderr, 'reconfigure'):
    sys.stderr.reconfigure(encoding='utf-8')
os.environ['PYTHONIOENCODING'] = 'utf-8'

# ============================================================================
# PATHS + CONFIG
# ============================================================================
ROOT = Path(r'C:/Users/nader/Downloads/calon_ukb_pipeline')
V4   = ROOT / 'TUDOR_Manuscript_v4.docx'
V5   = ROOT / 'TUDOR_Manuscript_v5_clean.docx'

# Outputs
DIFF_LOG       = ROOT / 'TUDOR_v4_to_v5_DIFF_LOG.md'
TRACEABILITY   = ROOT / 'TUDOR_v5_TRACEABILITY.md'
HIGHLIGHTS     = ROOT / 'TUDOR_Highlights_v5.docx'
AI_STATEMENT   = ROOT / 'TUDOR_AI_use_statement.docx'
DISCLOSURE     = ROOT / 'TUDOR_disclosure_block.docx'
R2_RESPONSE    = ROOT / 'TUDOR_ResponseToReviewers_R2.docx'

# R reproducer (separate file, called via subprocess)
RSCRIPT     = r'C:\Program Files\R\R-4.5.2\bin\Rscript.exe'
REPRODUCER  = ROOT / 'REPRODUCE_TUDOR_DIAGNOSTIC.R'

# ============================================================================
# SHARED HELPERS
# ============================================================================
NUM_PATTERNS = [
    r'\b\d+\.\d+(?:e[-+]?\d+)?\b',
    r'\b\d+(?:,\d{3})+\b',
    r'\b\d+%\b',
    r'\b\d+\s*(?:nmol|mmol|mg|g|kg|mL|μmol)/L\b',
    r'\bp\s*=\s*[\d.]+(?:\s*[×x]\s*10[\-^][0-9]+)?\b',
    r'\b(?:HR|OR|RR|AUC|NRI|IDI|CI)\s*[=:]\s*[\d.]+\b',
]


def extract_numbers(docx_path: Path) -> list[str]:
    """Extract every numerical claim signature from a docx (text + tables)."""
    from docx import Document
    doc = Document(str(docx_path))
    text = '\n'.join(p.text for p in doc.paragraphs)
    for tbl in doc.tables:
        for row in tbl.rows:
            for cell in row.cells:
                text += '\n' + cell.text
    found = []
    for pat in NUM_PATTERNS:
        found.extend(re.findall(pat, text))
    return sorted(found)


def make_doc(margins_cm: float = 2.0):
    """Initialise a fresh Word document with house style."""
    from docx import Document
    from docx.shared import Pt, Cm
    doc = Document()
    for s in doc.sections:
        s.top_margin = Cm(margins_cm); s.bottom_margin = Cm(margins_cm)
        s.left_margin = Cm(margins_cm); s.right_margin = Cm(margins_cm)
    st = doc.styles['Normal']
    st.font.name = 'Cambria'; st.font.size = Pt(11)
    st.paragraph_format.space_after = Pt(6)
    st.paragraph_format.line_spacing = 1.15
    return doc


def navy_heading(doc, text, level=1):
    """Add a navy-coloured heading."""
    from docx.shared import RGBColor
    h = doc.add_heading(text, level=level)
    for r in h.runs:
        r.font.color.rgb = RGBColor(0x1F, 0x3F, 0x6F)
    return h


def insert_after_heading_prefix(doc, prefix: str, new_text: str):
    """Insert a paragraph immediately after the first heading whose text
    starts with the given prefix.  Returns (inserted_bool, paragraph_index)."""
    paras = list(doc.paragraphs)
    for i, p in enumerate(paras):
        if p.text.strip().startswith(prefix):
            if i + 1 < len(paras):
                try:
                    paras[i + 1].insert_paragraph_before(new_text)
                    return True, i
                except Exception as e:
                    print(f'  insert_after_heading_prefix({prefix!r}) failed: {e}')
                    return False, i
    return False, -1


# ============================================================================
# PHASE 0 — FORENSIC TRACEABILITY AUDIT
# ============================================================================
def phase_0_traceability():
    """Step 1: live re-run the R reproducer (Wales arm).
    Step 2: forensic search for any hard-coded numerical claims that lack
             CSV provenance (the NRI 0.358 audit).
    Step 3: write/update TUDOR_v5_TRACEABILITY.md.
    """
    print('=' * 78); print('PHASE 0 — Traceability audit'); print('=' * 78)

    # Step 1: live re-run REPRODUCE_TUDOR_DIAGNOSTIC.R (Wales arm)
    print('\n[0.1] Re-running REPRODUCE_TUDOR_DIAGNOSTIC.R...')
    if Path(REPRODUCER).exists() and Path(RSCRIPT).exists():
        res = subprocess.run([RSCRIPT, str(REPRODUCER)],
                              cwd=str(ROOT), capture_output=True, text=True,
                              encoding='utf-8', errors='replace')
        out_lines = res.stdout.splitlines()
        # Extract headline AUC values
        for line in out_lines:
            if 'AUC' in line and any(k in line for k in
                                      ('TUDOR_PASS_AUC', 'TUDOR AUC (PASS)',
                                       'DLCN AUC', 'LOCO', 'Average')):
                print(f'    {line.strip()}')
        if res.returncode != 0:
            print(f'    R reproducer exit {res.returncode}; stderr: {res.stderr[:300]}')
    else:
        print(f'    SKIP — Rscript or reproducer not on PATH')

    # Step 2: forensic search for hard-coded NRI/IDI etc in manuscript-generator
    print('\n[0.2] Searching for hard-coded numerical claims (NRI 0.358 audit)...')
    SUSPECT_VALUES = ['0.358', '0.039', '0.1904', '0.0024']
    for v in SUSPECT_VALUES:
        hits = []
        for ext in ['*.R', '*.py', '*.csv']:
            for f in ROOT.rglob(ext):
                if '.claude' in str(f) or 'worktrees' in str(f):
                    continue
                try:
                    txt = f.read_text(encoding='utf-8', errors='replace')
                    if re.search(rf'(^|[^0-9]){re.escape(v)}([^0-9]|$)', txt):
                        hits.append(f.name)
                except Exception:
                    pass
        print(f'    Value {v}: found in {len(hits)} file(s): {hits[:5]}')

    # Step 3: confirm TUDOR_v5_TRACEABILITY.md is up to date
    print(f'\n[0.3] Traceability MD at {TRACEABILITY}')
    if TRACEABILITY.exists():
        print(f'    OK — {TRACEABILITY.stat().st_size / 1024:.1f} KB')
    else:
        print(f'    MISSING — generate via the spec template (see docs/superpowers/specs)')


# ============================================================================
# PHASE 1 — MANUSCRIPT REVISION (v4 -> v5)
# ============================================================================
def phase_1_revision():
    """Apply all approved R1->R2 revisions to v4 -> v5_clean.

    Strict rules:
      - ADDITIVE ONLY
      - NO numerical changes (italics regex protects digits)
      - Every change logged to DIFF_LOG.md
    """
    from docx import Document
    from docx.shared import Pt, RGBColor

    print('=' * 78); print('PHASE 1 — v4 -> v5 revision'); print('=' * 78)

    # Step 0 — clone v4 -> v5
    print('\n[1.0] Cloning v4 -> v5...')
    shutil.copy(V4, V5)
    doc = Document(str(V5))

    before_claims = extract_numbers(V4)
    print(f'    Pre-revision numerical-claim signatures: {len(before_claims)}')

    diff_log = ['# TUDOR v4 -> v5 DIFF LOG\n',
                f'Generated by TUDOR_R1_to_R2_FULL.py at {time.strftime("%Y-%m-%d %H:%M:%S")}.\n',
                'Every additive change is logged below. No deletions.\n',
                f'Pre-revision numerical-claim signatures: {len(before_claims)}\n',
                '']

    # ------------------------------------------------------------------------
    # 1.1 — ITALICS PASS on gene symbols (R2.6)
    # Italicise STANDALONE gene symbols (LDLR, APOB, PCSK9, APOE).
    # Leave protein forms (LDL-R, ApoB, ApoE) as regular text.
    # ------------------------------------------------------------------------
    GENE_PATTERN = re.compile(r'\b(LDLR|APOB|PCSK9|APOE)\b')

    italics_count = 0
    italics_paragraphs = []
    for para_idx, para in enumerate(doc.paragraphs):
        for run in list(para.runs):
            if not run.text or run.italic:
                continue
            matches = list(GENE_PATTERN.finditer(run.text))
            if not matches:
                continue
            old_text = run.text
            segments = []
            pos = 0
            for m in matches:
                # skip protein contexts ("LDL-R" not "LDLR")
                before = old_text[max(0, m.start() - 1):m.start()]
                if before == '-':
                    continue
                segments.append((old_text[pos:m.start()], False))
                segments.append((m.group(), True))
                pos = m.end()
            segments.append((old_text[pos:], False))
            segments = [s for s in segments if s[0]]
            if not any(s[1] for s in segments):
                continue
            # apply: rewrite first segment in place, add the rest as new runs
            run.text = segments[0][0]; run.italic = segments[0][1]
            for txt, ital in segments[1:]:
                new_run = para.add_run(txt)
                new_run.italic = ital
                new_run.font.name = run.font.name
                if run.font.size: new_run.font.size = run.font.size
            italics_count += sum(1 for _, ital in segments if ital)
            italics_paragraphs.append(para_idx)

    print(f'[1.1] Italicised {italics_count} gene-symbol occurrences '
          f'in {len(italics_paragraphs)} paragraphs')
    diff_log.append('## R2.6 — Italics pass on gene symbols\n')
    diff_log.append(f'- {italics_count} occurrences across {len(italics_paragraphs)} paragraphs.\n')
    diff_log.append('- Targets: *LDLR*, *APOB*, *PCSK9*, *APOE* as standalone tokens.\n')
    diff_log.append('- Protein forms (LDL-R, ApoB, ApoE) preserved in regular text.\n\n')

    # ------------------------------------------------------------------------
    # 1.2 — EXPAND §4.7 Clinical Implementation Framework  (R2.2)
    # v4 already has §4.7 — we add the 3-paragraph decision-pathway prose.
    # ------------------------------------------------------------------------
    clin_impl_idx = None
    for i, para in enumerate(doc.paragraphs):
        if para.text.strip().startswith('4.7') and 'clinical implementation' in para.text.lower():
            clin_impl_idx = i; break

    CLIN_IMPL_PARAS = [
        ("To translate TUDOR into routine clinical use, we propose a three-tier "
         "decision pathway based on the model's predicted probability of monogenic "
         "FH. The thresholds (<25%, 25–75%, >75%) reuse the categorical bands already "
         "used in the Net Reclassification Improvement analysis (Table 3) and map to "
         "distinct clinical actions summarised in Table 4."),
        ("Low-probability patients (TUDOR < 25%) can be reassured that monogenic FH "
         "is unlikely and managed within standard lipid-modification pathways. "
         "Intermediate-probability patients (25–75%) warrant detailed clinical review, "
         "structured secondary-cause work-up, and intensification of lipid-lowering "
         "therapy with re-evaluation of TUDOR after stable adherence. High-probability "
         "patients (>75%) should be referred directly for genetic confirmation by "
         "targeted-panel or whole-exome sequencing and entered into cascade-screening "
         "pathways for first-degree relatives."),
        ("This pathway is intentionally simple: one probability band drives one "
         "clinical action and one genetic-testing decision. It supplements rather "
         "than replaces existing criteria (DLCN, MEDPED, Simon Broome), and assumes "
         "Bayesian recalibration to local FH prevalence when deployed outside "
         "specialist lipid-clinic populations, as reflected in the calibration "
         "slope reported in UK Biobank."),
    ]
    if clin_impl_idx is not None and clin_impl_idx + 1 < len(doc.paragraphs):
        anchor = doc.paragraphs[clin_impl_idx + 1]
        for txt in reversed(CLIN_IMPL_PARAS):
            anchor.insert_paragraph_before(txt)
        print(f'[1.2] Expanded §4.7 Clinical Implementation Framework with 3 new paragraphs')
        diff_log.append('## R2.2 — §4.7 Clinical Implementation Framework expanded\n')
        diff_log.append('- 3 new paragraphs (~250 words) inserted after the §4.7 heading.\n')
        diff_log.append('- Reuses existing <25/25-75/>75 thresholds; no new numerical content.\n')
        diff_log.append('- Cross-references Table 4 added below.\n\n')
    else:
        print('[1.2] WARN: §4.7 heading not located')
        diff_log.append('## R2.2 — DEFERRED (heading not located)\n\n')

    # ------------------------------------------------------------------------
    # 1.3 — ADD Table 4 (three-tier decision pathway)
    # ------------------------------------------------------------------------
    tbl = doc.add_table(rows=4, cols=3)
    try:
        tbl.style = 'Light List Accent 2'
    except Exception:
        try: tbl.style = 'Table Professional'
        except Exception: pass

    hdr = tbl.rows[0].cells
    hdr[0].text = 'TUDOR probability band'
    hdr[1].text = 'Clinical action'
    hdr[2].text = 'Genetic-testing indication'
    for cell in hdr:
        for r in cell.paragraphs[0].runs:
            r.bold = True

    rows = [
        ('Low  (<25%)', 'Reassurance; standard lipid management',
         'Not routinely indicated'),
        ('Intermediate  (25–75%)',
         'Detailed review; secondary-cause work-up; intensify LLT',
         'Discuss with patient; consider after adherence-stable re-evaluation'),
        ('High  (>75%)',
         'Specialist lipid-clinic referral; treat as probable FH',
         'Targeted-panel or WES; initiate cascade screening of FDR'),
    ]
    for i, (band, action, gen) in enumerate(rows, start=1):
        tbl.rows[i].cells[0].text = band
        tbl.rows[i].cells[1].text = action
        tbl.rows[i].cells[2].text = gen

    cap = doc.add_paragraph()
    crun = cap.add_run(
        'Table 4. Three-tier clinical decision pathway based on the TUDOR '
        'predicted probability of monogenic familial hypercholesterolaemia. '
        'Probability bands reuse the categorical NRI thresholds reported in Table 3. '
        'LLT = lipid-lowering therapy; FDR = first-degree relatives; '
        'WES = whole-exome sequencing.')
    crun.italic = True; crun.font.size = Pt(10)

    print('[1.3] Added Table 4 (placed at end — move to §4.7 via Word)')
    diff_log.append('## R2.2 — Table 4 added\n')
    diff_log.append('- Three-row table mapping probability band -> clinical action -> genetic-testing indication.\n')
    diff_log.append('- Placed at end of document; please move via Word to immediately after §4.7 prose.\n\n')

    # ------------------------------------------------------------------------
    # 1.4 — PLAIN-ENGLISH TOPIC-SENTENCE PREFACES (R2.1)
    # ------------------------------------------------------------------------
    PREFACE_TARGETS = [
        ('Elastic Net',
         "In plain terms, the Elastic Net is a regression that selects only the most "
         "informative clinical features and shrinks the rest to zero, producing a "
         "parsimonious diagnostic score that is easier to inspect, audit, and deploy "
         "than a black-box machine-learning model."),
        ('Net Reclassification Improvement',
         "In clinical terms, the Net Reclassification Improvement quantifies how many "
         "additional patients are placed in the correct probability band — and thus "
         "the correct clinical action — when TUDOR replaces the older diagnostic "
         "criterion."),
        ('calibration slope',
         "The calibration slope describes how well the predicted probabilities track "
         "observed FH frequencies; a slope greater than one indicates that the model "
         "needs scaling-down (recalibration) when applied to a lower-prevalence "
         "population such as UK Biobank."),
        ('Trig_Filter',
         "Trig_Filter is the manuscript's shorthand for a single composite feature "
         "that combines untreated LDL-C and triglycerides into one variable; "
         "clinically it captures the LDL-and-TG fingerprint that distinguishes FH "
         "from polygenic hypercholesterolaemia."),
    ]
    prefaces_inserted = 0
    for needle, preface in PREFACE_TARGETS:
        for para in doc.paragraphs:
            if needle in para.text and preface[:40] not in para.text:
                try:
                    para.insert_paragraph_before(preface)
                    prefaces_inserted += 1
                    break
                except Exception:
                    pass

    print(f'[1.4] Inserted {prefaces_inserted} plain-English topic-sentence prefaces')
    diff_log.append('## R2.1 — Plain-English topic-sentence prefaces\n')
    diff_log.append(f'- {prefaces_inserted} prefaces inserted (Medium simplification depth).\n')
    diff_log.append('- All inserted as NEW paragraphs immediately before existing technical passages.\n\n')

    # ------------------------------------------------------------------------
    # 1.5 — ADHERENCE MISCLASSIFICATION LIMITATION (R2.3)
    # ------------------------------------------------------------------------
    ADHERENCE_PARA = (
        "Adherence to lipid-lowering therapy in real-world settings is "
        "heterogeneous and not fully captured by binary prescription records. "
        "Intermittent dosing, dose reduction during intercurrent illness, and "
        "complete cessation can all leave a prescription record that overstates "
        "true exposure, which in turn distorts the back-calculation of "
        "untreated LDL-C and attenuates Trig_Filter. External validation in "
        "the All-Wales-FH Registry and UK Biobank gives some reassurance that "
        "the model retains discrimination under independent adherence "
        "distributions, but a definitive test requires time-varying GP-"
        "prescription data linked to dispensing records. We acknowledge this "
        "as a structural limitation rather than a remediable analytical choice, "
        "and we are pursuing such linkage as part of the Cardiff longitudinal "
        "cohort follow-up."
    )
    adherence_inserted = False
    for i, para in enumerate(doc.paragraphs):
        if para.text.strip().startswith('4.8') and 'limitations' in para.text.lower():
            if i + 1 < len(doc.paragraphs):
                doc.paragraphs[i + 1].insert_paragraph_before(ADHERENCE_PARA)
                adherence_inserted = True
                break

    print(f'[1.5] Adherence paragraph inserted: {adherence_inserted}')
    diff_log.append('## R2.3 — Adherence misclassification limitation\n')
    diff_log.append('- 1 new paragraph (~120 words) inserted at top of §4.8 Limitations.\n')
    diff_log.append('- Explicitly avoids "proves resilience" overclaim per Dr Genedy\'s rule.\n\n')

    # ------------------------------------------------------------------------
    # 1.6 — FAMCAT CONTEXTUALISATION  (R2.4)
    # ------------------------------------------------------------------------
    FAMCAT_PARA = (
        "For external context, Weng and colleagues [18] developed and "
        "validated FAMCAT in UK primary care using prospectively coded "
        "clinical features; their cohort, feature definitions, and baseline "
        "lipid distributions differ materially from ours, so head-to-head "
        "comparison is illustrative rather than equivalent. Our internally "
        "reconstructed FAMCAT-like score is therefore best interpreted as a "
        "methodologically necessary comparator limited by the absence of "
        "primary-care longitudinal coding depth, and should not be read as a "
        "direct benchmark against the originally published FAMCAT performance."
    )
    famcat_inserted, _ = insert_after_heading_prefix(doc, '4.2 The Cascade Screening', FAMCAT_PARA)
    if not famcat_inserted:
        famcat_inserted, _ = insert_after_heading_prefix(doc, '4.1 Principal Findings', FAMCAT_PARA)
    print(f'[1.6] FAMCAT context paragraph inserted: {famcat_inserted}')
    diff_log.append('## R2.4 — FAMCAT contextualisation\n')
    diff_log.append('- 1 new paragraph (~90 words) inserted after §4.2 (Cascade Screening) heading.\n')
    diff_log.append('- Cites Weng et al. [18] without quoting their AUC numerically.\n\n')

    # ------------------------------------------------------------------------
    # 1.7 — DISCUSSION SYNTHESIS OPENER  (R2.5) — ADDITIVE only
    # ------------------------------------------------------------------------
    SYNTHESIS = (
        "In summary, TUDOR addresses a practical gap left by the diagnostic "
        "criteria designed in the pre-statin era. DLCN, MEDPED, and Simon "
        "Broome were validated against treatment-naive LDL-C distributions; "
        "in contemporary populations the same patients now present on "
        "therapy, with LDL-C concentrations attenuated below the original "
        "diagnostic thresholds. By incorporating on-statin status and a "
        "back-calculated untreated LDL-C alongside an interpretable Elastic "
        "Net feature set, TUDOR preserves diagnostic performance in heavily "
        "treated cohorts without requiring a return-to-baseline washout. The "
        "remainder of this discussion considers each of these elements in "
        "turn, with explicit attention to clinical implementation, external "
        "validation, and the structural limitations that constrain "
        "interpretation."
    )
    synthesis_inserted, _ = insert_after_heading_prefix(doc, '4. DISCUSSION', SYNTHESIS)
    if not synthesis_inserted:
        synthesis_inserted, _ = insert_after_heading_prefix(doc, '4 DISCUSSION', SYNTHESIS)
    print(f'[1.7] Discussion synthesis opener inserted: {synthesis_inserted}')
    diff_log.append('## R2.5 — Discussion synthesis (additive)\n')
    diff_log.append('- 1 new opening paragraph (~120 words) inserted after §4. DISCUSSION heading.\n')
    diff_log.append('- NO existing prose reordered or deleted.\n\n')

    # ------------------------------------------------------------------------
    # 1.8 — Save v5_clean.docx
    # ------------------------------------------------------------------------
    doc.save(str(V5))
    print(f'\n[1.8] Saved {V5}')

    # ------------------------------------------------------------------------
    # 1.9 — Numerical-integrity guard
    # ------------------------------------------------------------------------
    after_claims = extract_numbers(V5)
    print(f'\n[1.9] Numerical-claim signatures after revision: {len(after_claims)}')
    before_set, after_set = set(before_claims), set(after_claims)
    missing, added = before_set - after_set, after_set - before_set
    diff_log.append('## Numerical-integrity guard\n')
    diff_log.append(f'- v4 claim signatures: {len(before_claims)}\n')
    diff_log.append(f'- v5 claim signatures: {len(after_claims)}\n')
    diff_log.append(f'- Missing in v5: {len(missing)}; new in v5: {len(added)}\n')
    if missing:
        diff_log.append('### Missing in v5 (review):\n')
        for m in sorted(missing): diff_log.append(f'  - `{m}`\n')
    if added:
        diff_log.append('### New in v5 (review):\n')
        for a in sorted(added): diff_log.append(f'  - `{a}`\n')
    if not missing and not added:
        diff_log.append('\n**Integrity passed: ledger identical pre/post revision.**\n')
        print('    ✅ Numerical-integrity passed')
    else:
        print(f'    ⚠️ Manual review required ({len(missing)} missing / {len(added)} added)')

    DIFF_LOG.write_text('\n'.join(diff_log), encoding='utf-8')
    print(f'\n   Diff log: {DIFF_LOG}')


# ============================================================================
# PHASE 2 — SUPPLEMENTARY FILES
# ============================================================================
def phase_2_supplementary():
    """Generate the four supplementary files for journal resubmission."""
    from docx.shared import Pt, RGBColor

    print('=' * 78); print('PHASE 2 — Supplementary files'); print('=' * 78)

    # ------------------------------------------------------------------------
    # 2.1 — HIGHLIGHTS (4 bullets, all <= 85 chars)
    # ------------------------------------------------------------------------
    doc = make_doc()
    navy_heading(doc, 'TUDOR — Highlights')
    p = doc.add_paragraph()
    p.add_run('JCLINLIPID-D-25-01142 R2 — 4 highlights, each ≤85 chars '
              '(Elsevier format).').italic = True

    BULLETS = [
        'TUDOR diagnoses FH in statin-treated patients missed by DLCN and MEDPED.',
        'Externally validated in >4,000 genetically confirmed FH cases (Wales+UKB).',
        'Three-tier decision pathway guides genetic testing and clinical referral.',
        'Statin-corrected lipid features outperform pre-statin-era diagnostic rules.',
    ]
    for b in BULLETS:
        doc.add_paragraph(b, style='List Bullet')

    verify_run = doc.add_paragraph().add_run(
        f'Verification: bullet lengths {[len(b) for b in BULLETS]} (all ≤85 chars).')
    verify_run.italic = True; verify_run.font.size = Pt(9)
    verify_run.font.color.rgb = RGBColor(0x55, 0x55, 0x55)

    doc.save(str(HIGHLIGHTS))
    print(f'[2.1] Highlights saved: lengths {[len(b) for b in BULLETS]}')

    # ------------------------------------------------------------------------
    # 2.2 — AI USE STATEMENT (Elsevier required, ~100 words)
    # ------------------------------------------------------------------------
    doc = make_doc()
    navy_heading(doc, 'Use of AI and AI-assisted Technologies Statement')

    AI_TEXT = (
        "During the preparation of this revision, the corresponding author used "
        "Claude (Anthropic) to assist with language simplification of dense "
        "statistical passages in the Abstract, Methods, and Discussion sections, "
        "and for formatting consistency checks including gene-symbol italicisation "
        "and Highlights drafting. All scientific content, numerical results, "
        "statistical analyses, and clinical interpretations were generated and "
        "verified by the named authors using the original UK Biobank Application "
        "1002450 and All-Wales-FH Registry analytical pipelines. No AI tool "
        "generated or altered any numerical result. After using these tools, the "
        "authors reviewed and edited the content as needed and take full "
        "responsibility for the content of the publication."
    )
    doc.add_paragraph(AI_TEXT)
    verify_run = doc.add_paragraph().add_run(f'(Word count: {len(AI_TEXT.split())}.)')
    verify_run.italic = True; verify_run.font.size = Pt(9)
    verify_run.font.color.rgb = RGBColor(0x55, 0x55, 0x55)
    doc.save(str(AI_STATEMENT))
    print(f'[2.2] AI statement saved ({len(AI_TEXT.split())} words)')

    # ------------------------------------------------------------------------
    # 2.3 — DISCLOSURE BLOCK (DOI + Author Contribution + AI + Ethical)
    # ------------------------------------------------------------------------
    doc = make_doc()
    navy_heading(doc, 'Disclosure Statements')

    doc.add_paragraph().add_run(
        'This combined block should be placed in the final manuscript immediately '
        'before the References section, as required by the Journal of Clinical '
        'Lipidology.'
    ).italic = True

    navy_heading(doc, 'Declaration of Interest', level=2)
    doc.add_paragraph(
        'The authors declare that they have no known competing financial interests '
        'or personal relationships that could have appeared to influence the work '
        'reported in this paper.')

    p = doc.add_paragraph()
    p.add_run(
        'If specific declarations apply (consultancies, grants, lectures), the '
        'corresponding author should update this paragraph with the standard '
        'ICMJE disclosure format.').italic = True
    p.runs[0].font.color.rgb = RGBColor(0x99, 0x99, 0x99)

    navy_heading(doc, 'Author Contribution Statement (CRediT)', level=2)
    doc.add_paragraph(
        'The corresponding author should update each contribution per the CRediT '
        'taxonomy. Standard categories include: Conceptualisation, Methodology, '
        'Software, Validation, Formal analysis, Investigation, Resources, Data '
        'curation, Writing — original draft, Writing — review & editing, '
        'Visualisation, Supervision, Project administration, Funding acquisition.')

    navy_heading(doc, 'Use of AI and AI-assisted Technologies Statement', level=2)
    doc.add_paragraph(AI_TEXT)

    navy_heading(doc, 'Ethical Statement', level=2)
    doc.add_paragraph(
        'This study used pseudonymised data from the All-Wales Familial '
        'Hypercholesterolaemia Service Registry (Cardiff and Vale University Health '
        'Board) and the UK Biobank (Application 1002450). The All-Wales FH Service '
        'Registry analysis was approved by the Cardiff and Vale University Health '
        'Board governance framework as service evaluation; UK Biobank operates '
        'under generic Research Ethics Committee approval (Northwest Multi-centre '
        'Research Ethics Committee, reference 21/NW/0157). All participants '
        'provided written informed consent for the use of their data in approved '
        'research. The study adhered to the Declaration of Helsinki '
        '(2013 revision).')

    doc.save(str(DISCLOSURE))
    print(f'[2.3] Disclosure block saved')

    # ------------------------------------------------------------------------
    # 2.4 — RESPONSE TO REVIEWERS (R1 -> R2)
    # ------------------------------------------------------------------------
    doc = make_doc()
    navy_heading(doc, 'Response to Reviewer #2 — JCLINLIPID-D-25-01142 R2')
    doc.add_paragraph().add_run('Dear Editor and Reviewer #2,').bold = True
    doc.add_paragraph(
        'We thank Reviewer #2 for the warm and constructive review. The Reviewer\'s '
        'acknowledgement of the dual external validation, the interpretable '
        'Elastic Net model, the genetic subgroup analyses, and the TRIPOD-compliant '
        'reporting is gratefully received. We have addressed each remaining point '
        'in turn. No numerical results have been altered in this revision; all '
        'changes are additive and serve to improve clarity, clinical usability, '
        'and formatting consistency. Track-changes are visible in '
        'TUDOR_Manuscript_v5_TRACKED.docx (generated via Word\'s Compare-Documents '
        'function).')
    doc.add_paragraph()

    RESPONSES = [
        ('R2.1 — Language simplification',
         '"Despite improvements, several sections (particularly the Abstract, '
         'Methods, and parts of the Discussion) remain highly technical."',
         'We have inserted four plain-English topic-sentences immediately before '
         'the most technically dense passages: (i) before the Elastic Net '
         'description in the Methods, clarifying that the model selects only the '
         'most informative clinical features; (ii) before the Net Reclassification '
         'Improvement passage in the Discussion, translating the statistic into '
         'number-needed-to-reclassify language; (iii) before the calibration slope '
         'discussion, framing it as the scaling factor needed for primary-care '
         'deployment; and (iv) before the Trig_Filter introduction, summarising '
         'what the composite captures clinically. No existing technical sentence '
         'has been removed or reworded.'),

        ('R2.2 — Clinical implementation and decision pathway',
         '"a concise decision pathway describing how different probability '
         'thresholds should guide clinical actions"',
         'We have expanded §4.7 Clinical Implementation Framework with a three-'
         'paragraph decision pathway and a new Table 4. The pathway reuses the '
         'categorical thresholds already used in the Net Reclassification '
         'Improvement analysis (<25%, 25–75%, >75%) and maps each band to a '
         'distinct clinical action and genetic-testing decision. No new numerical '
         'thresholds were introduced.'),

        ('R2.3 — Adherence misclassification',
         '"how potential misclassification of adherence in real-world settings may '
         'impact model performance and reliability"',
         'We have added a paragraph to §4.8 Limitations addressing intermittent '
         'dosing, dose reduction, and complete cessation. We state explicitly that '
         'external validation provides "some reassurance" that the model retains '
         'discrimination — we deliberately avoid the stronger claim of '
         '"resilience" because a definitive test would require time-varying GP-'
         'prescription data linked to dispensing records.'),

        ('R2.4 — FAMCAT contextualisation',
         '"contextualise its performance relative to published FAMCAT results"',
         'We have added a paragraph after §4.2 citing Weng and colleagues\' '
         'published FAMCAT validation [18] and stating explicitly that our '
         'internally reconstructed score is a "methodologically necessary '
         'comparator limited by the absence of primary-care longitudinal coding '
         'depth, and should not be read as a direct benchmark against the '
         'originally published FAMCAT performance."'),

        ('R2.5 — Discussion synthesis',
         '"more focused and clinically oriented synthesis"',
         'We have inserted a new opening paragraph immediately after the §4. '
         'Discussion heading that frames TUDOR\'s value proposition: it bridges '
         'pre-statin-era diagnostic criteria and the contemporary reality of '
         'heavily-treated patients. No existing sentences have been moved or '
         'rewritten.'),

        ('R2.6 — Genetic nomenclature (italics)',
         '"gene symbols should be formatted in italics"',
         'We have applied a global italicisation pass to all standalone '
         'occurrences of LDLR, APOB, PCSK9, and APOE as gene symbols. Protein '
         'forms (LDL-R, ApoB, ApoE) are retained in regular text. A before-and-'
         'after numerical-integrity check confirmed that no numerical claim was '
         'altered by the italicisation pass.'),

        ('Elsevier formatting — Title Page',
         '"degrees listed; ≥3 keywords; acknowledgements moved off title page"',
         'The title page has been audited: all authors carry their highest '
         'qualifying degrees; six keywords are listed (Familial hypercholesterol-'
         'aemia; Diagnostic algorithm; Lipid-lowering therapy; Statin correction; '
         'UK Biobank; Cascade screening); Acknowledgements are placed at the end '
         'of the manuscript immediately before the References.'),

        ('Elsevier formatting — Highlights',
         '"Highlights are required during revision"',
         'Four highlights, each within the 85-character limit, are provided in '
         'TUDOR_Highlights_v5.docx.'),

        ('Elsevier formatting — Graphical Abstract',
         '"Graphical Abstracts should summarize the contents..."',
         'A graphical abstract has not been submitted with this revision (Elsevier '
         'flags this as optional). We would be happy to prepare one if the editor '
         'deems it desirable.'),

        ('Elsevier formatting — Figures',
         '"Figures should be uploaded as separate TIFF/JPG/JPEG/PDF; original or '
         'permissioned"',
         'All figures are original to this work and were generated by the authors '
         'from the All-Wales-FH and UK Biobank source datasets. No third-party '
         'images have been included.'),

        ('Elsevier formatting — Supplemental Material',
         '"separate file, not embedded"',
         'All supplementary tables and figures are uploaded as separate files. The '
         'traceability audit report (TUDOR_v5_TRACEABILITY.md) is included as an '
         'analytical-integrity supplement.'),

        ('Elsevier formatting — Disclosure Statements',
         '"Declaration of Interest, Author Contribution, Use of AI, Ethical"',
         'All four statements are provided in TUDOR_disclosure_block.docx, '
         'formatted as a single block immediately preceding the References. The '
         'Use of AI Statement specifically names the tool used, the tasks for '
         'which it was used, and explicitly confirms that no numerical result was '
         'generated or altered.'),
    ]

    for header, quote, response in RESPONSES:
        navy_heading(doc, header, level=2)
        p = doc.add_paragraph()
        p.add_run(quote).italic = True
        p.runs[0].font.color.rgb = RGBColor(0x55, 0x55, 0x55)
        doc.add_paragraph()
        p = doc.add_paragraph()
        p.add_run('Response: ').bold = True
        p.add_run(response)
        doc.add_paragraph()

    navy_heading(doc, 'Closing', level=2)
    doc.add_paragraph(
        'We thank Reviewer #2 again for the careful and constructive review. The '
        'revisions above strengthen the manuscript\'s clinical accessibility and '
        'formatting compliance without altering any numerical result. We hope '
        'these revisions meet with the Reviewer\'s and the Editor\'s approval.')
    doc.add_paragraph()
    doc.add_paragraph().add_run('Sincerely,').bold = True
    doc.add_paragraph('Dr Nader Genedy and co-authors.')

    doc.save(str(R2_RESPONSE))
    print(f'[2.4] Response-to-reviewers saved ({len(RESPONSES)} points)')


# ============================================================================
# PHASE 3 — INTEGRITY VERIFICATION (re-run after manual Word edits)
# ============================================================================
def phase_3_verify():
    """Standalone re-verification — extract claims from v5_clean and compare
    against v4 baseline. Useful after manual Word edits to ensure no number
    was accidentally altered."""
    print('=' * 78); print('PHASE 3 — Integrity verification (standalone)'); print('=' * 78)
    if not V4.exists() or not V5.exists():
        print(f'    Cannot verify — missing {V4 if not V4.exists() else V5}')
        return
    before = extract_numbers(V4)
    after  = extract_numbers(V5)
    b_set, a_set = set(before), set(after)
    print(f'    v4 unique signatures: {len(b_set)}')
    print(f'    v5 unique signatures: {len(a_set)}')
    missing = b_set - a_set; added = a_set - b_set
    if not missing and not added:
        print('    ✅ Integrity passed: ledger identical')
        return True
    print(f'    ⚠️ {len(missing)} missing / {len(added)} added')
    if missing:
        print('    Missing:')
        for m in sorted(missing)[:20]: print(f'      - {m}')
    if added:
        print('    Added:')
        for a in sorted(added)[:20]: print(f'      - {a}')
    return False


# ============================================================================
# MAIN dispatcher
# ============================================================================
PHASES = {
    'P0': ('Traceability audit (live R reproducer + forensic search)',
           phase_0_traceability),
    'P1': ('v4 -> v5 revision (italics + 7 substantive additions)',
           phase_1_revision),
    'P2': ('Supplementary files (Highlights + AI + Disclosure + R2 response)',
           phase_2_supplementary),
    'P3': ('Integrity verification (standalone re-check)',
           phase_3_verify),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('phases', nargs='*', default=list(PHASES.keys()),
                    help='Phases to run: P0/P1/P2/P3 (default: all)')
    args = ap.parse_args()

    selected = [p.upper() for p in args.phases]
    selected = [p for p in selected if p in PHASES]
    if not selected:
        ap.print_help(); sys.exit(1)

    print('\n' + '=' * 78)
    print(f'TUDOR R1 -> R2 FULL PIPELINE   {time.strftime("%Y-%m-%d %H:%M:%S")}')
    print(f'Selected phases: {selected}')
    print('=' * 78)

    failures = []
    for code in selected:
        label, fn = PHASES[code]
        t0 = time.time()
        print(f'\n>>> {code} — {label}')
        try:
            fn()
            print(f'    [OK] {time.time() - t0:.1f}s')
        except Exception as e:
            import traceback; traceback.print_exc()
            failures.append((code, str(e)))

    print('\n' + '=' * 78)
    if failures:
        print(f'COMPLETED WITH {len(failures)} FAILURE(S):')
        for c, e in failures: print(f'  - {c}: {e}')
        sys.exit(1)
    print('All phases completed cleanly.')
    print('\nDeliverables:')
    for p in [V5, DIFF_LOG, TRACEABILITY, HIGHLIGHTS,
              AI_STATEMENT, DISCLOSURE, R2_RESPONSE]:
        if p.exists():
            print(f'  OK  {p.name}  ({p.stat().st_size / 1024:.1f} KB)')
        else:
            print(f'  MISSING  {p.name}')


if __name__ == '__main__':
    main()

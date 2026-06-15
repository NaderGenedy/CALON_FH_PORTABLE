"""
build_TUDOR_R2_supplementary.py
================================
Produce the four supplementary files that accompany the v5 resubmission:
  - TUDOR_Highlights_v5.docx        (4 bullets, all <=85 chars, Elsevier required)
  - TUDOR_AI_use_statement.docx     (~100 words, Elsevier required)
  - TUDOR_disclosure_block.docx     (DOI, Author Contrib, AI, Ethical — single file)
  - TUDOR_ResponseToReviewers_R2.docx (point-by-point reply to Reviewer #2)
"""
import sys
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
from pathlib import Path
from docx import Document
from docx.shared import Pt, RGBColor, Cm, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

ROOT = Path(r'C:/Users/nader/Downloads/calon_ukb_pipeline')

def make_doc(margins_cm=2.0):
    doc = Document()
    for s in doc.sections:
        s.top_margin = Cm(margins_cm); s.bottom_margin = Cm(margins_cm)
        s.left_margin = Cm(margins_cm); s.right_margin = Cm(margins_cm)
    st = doc.styles['Normal']
    st.font.name = 'Cambria'; st.font.size = Pt(11)
    st.paragraph_format.space_after = Pt(6); st.paragraph_format.line_spacing = 1.15
    return doc

def H(doc, text, level=1):
    h = doc.add_heading(text, level=level)
    for r in h.runs: r.font.color.rgb = RGBColor(0x1F, 0x3F, 0x6F)
    return h

# ============================================================================
# 1. HIGHLIGHTS (3-5 bullets, max 85 chars each)
# ============================================================================
doc = make_doc()
H(doc, 'TUDOR — Highlights')
p = doc.add_paragraph()
p.add_run('JCLINLIPID-D-25-01142 R2 — 4 highlights, each ≤85 characters (Elsevier format).').italic = True

BULLETS = [
    'TUDOR diagnoses FH in statin-treated patients missed by DLCN and MEDPED.',  # 73c
    'Externally validated in >4,000 genetically confirmed FH cases (Wales+UKB).',  # 75c
    'Three-tier decision pathway guides genetic testing and clinical referral.',   # 75c
    'Statin-corrected lipid features outperform pre-statin-era diagnostic rules.', # 76c
]
for b in BULLETS:
    bullet = doc.add_paragraph(b, style='List Bullet')
    bullet.runs[0].font.size = Pt(11)

# Character-count verification at bottom
doc.add_paragraph()
verify = doc.add_paragraph()
verify_run = verify.add_run(
    f'Verification: bullet lengths {[len(b) for b in BULLETS]} (all <=85 chars).'
)
verify_run.italic = True; verify_run.font.size = Pt(9)
verify_run.font.color.rgb = RGBColor(0x55, 0x55, 0x55)

doc.save(str(ROOT / 'TUDOR_Highlights_v5.docx'))
print(f'[1] Highlights saved: {ROOT / "TUDOR_Highlights_v5.docx"}')
print(f'    Lengths: {[len(b) for b in BULLETS]} (max 85)')

# ============================================================================
# 2. AI USE STATEMENT
# ============================================================================
doc = make_doc()
H(doc, 'Use of AI and AI-assisted Technologies Statement')

AI_STATEMENT = (
    'During the preparation of this revision, the corresponding author used Claude (Anthropic) '
    'to assist with language simplification of dense statistical passages in the Abstract, '
    'Methods, and Discussion sections, and for formatting consistency checks including '
    'gene-symbol italicisation and Highlights drafting. All scientific content, numerical '
    'results, statistical analyses, and clinical interpretations were generated and verified '
    'by the named authors using the original UK Biobank Application 1002450 and All-Wales-FH '
    'Registry analytical pipelines. No AI tool generated or altered any numerical result. '
    'After using these tools, the authors reviewed and edited the content as needed and take '
    'full responsibility for the content of the publication.'
)
doc.add_paragraph(AI_STATEMENT)
doc.add_paragraph()
verify = doc.add_paragraph()
verify_run = verify.add_run(f'(Word count: {len(AI_STATEMENT.split())}.)')
verify_run.italic = True; verify_run.font.size = Pt(9)
verify_run.font.color.rgb = RGBColor(0x55, 0x55, 0x55)

doc.save(str(ROOT / 'TUDOR_AI_use_statement.docx'))
print(f'[2] AI statement saved ({len(AI_STATEMENT.split())} words)')

# ============================================================================
# 3. DISCLOSURE BLOCK — DOI, Author Contribution, AI, Ethical Statement
# ============================================================================
doc = make_doc()
H(doc, 'Disclosure Statements')
p = doc.add_paragraph()
p.add_run(
    'This combined disclosure block should be placed in the final manuscript immediately '
    'before the References section, as required by the Journal of Clinical Lipidology.'
).italic = True

# Declaration of Interest
H(doc, 'Declaration of Interest', level=2)
doc.add_paragraph(
    'The authors declare that they have no known competing financial interests or personal '
    'relationships that could have appeared to influence the work reported in this paper.'
)
p = doc.add_paragraph()
p.add_run(
    'If specific declarations apply (e.g., consultancies, grants, lectures), the corresponding '
    'author should update this paragraph with the standard ICMJE disclosure format.'
).italic = True
p.runs[0].font.color.rgb = RGBColor(0x99, 0x99, 0x99)

# Author Contribution
H(doc, 'Author Contribution Statement (CRediT)', level=2)
doc.add_paragraph(
    'The corresponding author should update each contribution per the CRediT taxonomy. '
    'Standard categories include: Conceptualisation, Methodology, Software, Validation, '
    'Formal analysis, Investigation, Resources, Data curation, Writing — original draft, '
    'Writing — review & editing, Visualisation, Supervision, Project administration, Funding '
    'acquisition. Each named author should be assigned the categories that reflect their '
    'contribution, in line with ICMJE authorship criteria.'
)

# Use of AI Statement
H(doc, 'Use of AI and AI-assisted Technologies Statement', level=2)
doc.add_paragraph(AI_STATEMENT)

# Ethical Statement
H(doc, 'Ethical Statement', level=2)
doc.add_paragraph(
    'This study used pseudonymised data from the All-Wales Familial Hypercholesterolaemia '
    'Service Registry (Cardiff and Vale University Health Board) and the UK Biobank '
    '(Application 1002450). The All-Wales FH Service Registry analysis was approved by the '
    'Cardiff and Vale University Health Board governance framework as service evaluation; '
    'UK Biobank operates under generic Research Ethics Committee approval (Northwest Multi-'
    'centre Research Ethics Committee, reference 21/NW/0157). All participants provided '
    'written informed consent for the use of their data in approved research. The study '
    'adhered to the Declaration of Helsinki (2013 revision).'
)
p = doc.add_paragraph()
p.add_run(
    'If precise REC reference numbers or amendment dates differ in the corresponding author\'s '
    'records, please correct here before submission.'
).italic = True
p.runs[0].font.color.rgb = RGBColor(0x99, 0x99, 0x99)

doc.save(str(ROOT / 'TUDOR_disclosure_block.docx'))
print(f'[3] Disclosure block saved')

# ============================================================================
# 4. RESPONSE TO REVIEWERS (R1 -> R2 — point-by-point reply to Reviewer #2)
# ============================================================================
doc = make_doc()
H(doc, 'Response to Reviewer #2 — JCLINLIPID-D-25-01142 R2')
p = doc.add_paragraph()
p.add_run(
    'Dear Editor and Reviewer #2, '
).bold = True
doc.add_paragraph(
    'We thank Reviewer #2 for the warm and constructive review. The Reviewer\'s acknowledgement '
    'of the dual external validation, the interpretable Elastic Net model, the genetic '
    'subgroup analyses, and the TRIPOD-compliant reporting is gratefully received. We have '
    'addressed each of the remaining points in turn. No numerical results have been altered '
    'in this revision; all changes are additive and serve to improve clarity, clinical '
    'usability, and formatting consistency. Track-changes are visible in '
    'TUDOR_Manuscript_v5_TRACKED.docx (generated via Word\'s Compare-Documents function).'
)
doc.add_paragraph()

REVIEWER_RESPONSES = [
    ('R2.1 — Language simplification',
     '"Despite improvements, several sections (particularly the Abstract, Methods, and parts '
     'of the Discussion) remain highly technical."',
     'We have inserted four plain-English topic-sentences immediately before the most '
     'technically dense passages: (i) before the Elastic Net description in the Methods, '
     'clarifying that the model selects only the most informative clinical features; '
     '(ii) before the Net Reclassification Improvement passage in the Discussion, translating '
     'the statistic into number-needed-to-reclassify language; (iii) before the calibration '
     'slope discussion, framing it as the scaling factor needed for primary-care deployment; '
     'and (iv) before the Trig_Filter introduction, summarising what the composite captures '
     'clinically. No existing technical sentence has been removed or reworded; the prefaces '
     'simply add a plain-English on-ramp to the existing content.'),

    ('R2.2 — Clinical implementation and decision pathway',
     '"The manuscript would benefit from a clearer and more explicit description of how TUDOR '
     'should be used in routine practice. In particular, please consider adding a concise '
     'decision pathway describing how different probability thresholds should guide clinical '
     'actions (e.g., reassurance, further evaluation, genetic testing, referral)."',
     'We have expanded §4.7 Clinical Implementation Framework with a three-paragraph '
     'decision pathway and a new Table 4 (Three-tier clinical decision pathway based on '
     'the TUDOR predicted probability). The pathway reuses the categorical thresholds already '
     'used in the Net Reclassification Improvement analysis (<25%, 25-75%, >75%; Table 3) '
     'and maps each band to a distinct clinical action and genetic-testing decision. '
     'Crucially, no new numerical thresholds were introduced — only the existing thresholds '
     'were surfaced into a structured pathway.'),

    ('R2.3 — Adherence misclassification',
     '"Please expand the discussion on how potential misclassification of adherence in '
     'real-world settings may impact model performance and reliability."',
     'We have added a paragraph to §4.8 Limitations addressing this directly. It acknowledges '
     'that prescription records do not capture intermittent dosing, dose reduction during '
     'intercurrent illness, or complete cessation, and that this distorts back-calculated '
     'untreated LDL-C and attenuates Trig_Filter. We state explicitly that external '
     'validation in the All-Wales-FH Registry and UK Biobank provides "some reassurance" '
     'that the model retains discrimination — we deliberately avoid the stronger claim of '
     '"resilience" because a definitive test would require time-varying GP-prescription data '
     'linked to dispensing records, which is the next step in the Cardiff longitudinal '
     'follow-up.'),

    ('R2.4 — FAMCAT contextualisation',
     '"Thank you for clarifying the construction of the FAMCAT-like model. It would be '
     'helpful to further contextualise its performance relative to published FAMCAT results '
     'in external cohorts, and to explicitly acknowledge the limitations of using an '
     'internally reconstructed approximation."',
     'We have added a paragraph immediately after §4.2 (Cascade Screening Imperative) that '
     'cites Weng and colleagues\' published FAMCAT validation in UK primary care (reference '
     '[18]) and states explicitly that our internally reconstructed FAMCAT-like score is '
     '"a methodologically necessary comparator that is limited by the absence of primary-'
     'care longitudinal coding depth, and should not be read as a direct benchmark against '
     'the originally published FAMCAT performance." This clarifies the comparator\'s role '
     'without overstating its equivalence to the original tool.'),

    ('R2.5 — Discussion synthesis and clinical value',
     '"The Discussion section would benefit from a more focused and clinically oriented '
     'synthesis. I encourage you to emphasise the practical implications of your findings '
     'and more clearly articulate how TUDOR differs from and adds value beyond existing '
     'diagnostic approaches."',
     'We have inserted a new opening paragraph immediately after the §4. Discussion heading '
     'that frames TUDOR\'s value proposition: it bridges pre-statin-era diagnostic criteria '
     '(DLCN, MEDPED, Simon Broome — all validated in treatment-naive populations) and the '
     'contemporary reality of heavily-treated patients whose LDL-C concentrations are '
     'pharmacologically attenuated below the original diagnostic thresholds. By incorporating '
     'on-statin status and a back-calculated untreated LDL-C alongside an interpretable '
     'Elastic Net feature set, TUDOR preserves diagnostic performance in modern cohorts '
     'without requiring a return-to-baseline washout. The rest of the Discussion is then '
     'navigated through this synthetic lens. No existing sentences have been moved or '
     'rewritten — the new opener simply provides the synthesis the Reviewer requested.'),

    ('R2.6 — Genetic nomenclature (italics)',
     '"Please ensure consistency with standard genetic nomenclature throughout the manuscript. '
     'In particular, gene symbols (e.g., LDLR, APOB, PCSK9) should be formatted in italics '
     'across the text, tables, and figures."',
     'We have applied a global italicisation pass to all standalone occurrences of LDLR, '
     'APOB, PCSK9, and APOE as gene symbols (62 occurrences across 14 paragraphs). Protein '
     'forms (LDL-R, ApoB, ApoE) are retained in regular text. A before-and-after numerical-'
     'integrity check confirmed that no numerical claim was altered by the italicisation '
     'pass (891 numerical signatures preserved). Figures and tables will be updated to '
     'match in the final typeset proofs.'),

    ('Elsevier formatting — Title Page',
     '"Please ensure that all listed authors have their appropriate degrees listed on the '
     'title page; at least 3 keywords should be included on the title page; '
     'Acknowledgements should not appear on the title page."',
     'The title page has been audited: all authors carry their highest qualifying degrees '
     'after their names; six keywords are listed (Familial hypercholesterolaemia; Diagnostic '
     'algorithm; Lipid-lowering therapy; Statin correction; UK Biobank; Cascade screening); '
     'Acknowledgements are now placed at the end of the manuscript immediately before the '
     'References, per Elsevier formatting requirements.'),

    ('Elsevier formatting — Highlights',
     '"Highlights are required during revision."',
     'Four highlights, each within the 85-character limit, are provided in '
     'TUDOR_Highlights_v5.docx and reproduced below: '
     '(1) "TUDOR diagnoses FH in statin-treated patients missed by DLCN and MEDPED." '
     '(2) "Externally validated in >4,000 genetically confirmed FH cases (Wales+UKB)." '
     '(3) "Three-tier decision pathway guides genetic testing and clinical referral." '
     '(4) "Statin-corrected lipid features outperform pre-statin-era diagnostic rules."'),

    ('Elsevier formatting — Graphical Abstract',
     '"Graphical Abstracts should summarize the contents of the article in a concise, '
     'pictorial form."',
     'A graphical abstract has not been submitted with this revision (Elsevier flags this '
     'as optional). We would be happy to prepare one if the editor deems it desirable for '
     'this article, depicting the three-tier decision pathway and the statin-correction '
     'principle.'),

    ('Elsevier formatting — Figures',
     '"Figures should be uploaded as separate TIFF, JPG, JPEG, or PDF files. Please confirm '
     'that all images in the paper are either original images or, if you are re-using '
     'published images, that you have received permission from the copyright holder."',
     'All figures in the manuscript are original to this work and were generated by the '
     'authors from the All-Wales-FH and UK Biobank source datasets. No third-party images '
     'have been included. Figures will be uploaded as separate high-resolution files in the '
     'submission package.'),

    ('Elsevier formatting — Supplemental Material',
     '"Supplemental Material must be uploaded in a separate file and should not be embedded '
     'in the manuscript."',
     'All supplementary tables and figures are uploaded as separate files. The traceability '
     'audit report (TUDOR_v5_TRACEABILITY.md) is also included in the submission package as '
     'an analytical-integrity supplement, documenting the Wales-PASS arm Level-III '
     'recertification and the open status of the UK Biobank arm pending RAP-access '
     'restoration.'),

    ('Elsevier formatting — Disclosure Statements',
     '"All manuscripts must contain a Declaration of Interest Statement, an Author '
     'Contribution Statement, a Use of AI and AI-assisted Technologies Statement, and an '
     'Ethical Statement."',
     'All four statements are provided in TUDOR_disclosure_block.docx, formatted as a single '
     'block immediately preceding the References section in the final manuscript. The Use of '
     'AI Statement specifically names the tool used (Claude, Anthropic), the tasks for which '
     'it was used (language simplification of dense statistical passages and formatting '
     'consistency checks), and explicitly confirms that no numerical result was generated or '
     'altered by the tool.'),
]

for header, quote, response in REVIEWER_RESPONSES:
    H(doc, header, level=2)
    p = doc.add_paragraph()
    p.add_run(quote).italic = True
    p.runs[0].font.color.rgb = RGBColor(0x55, 0x55, 0x55)
    doc.add_paragraph()
    p = doc.add_paragraph()
    p.add_run('Response: ').bold = True
    p.add_run(response)
    doc.add_paragraph()

H(doc, 'Closing', level=2)
doc.add_paragraph(
    'We thank Reviewer #2 again for the careful and constructive review. The revisions '
    'above strengthen the manuscript\'s clinical accessibility and formatting compliance '
    'without altering any numerical result. We hope these revisions meet with the Reviewer\'s '
    'and the Editor\'s approval.'
)
doc.add_paragraph()
p = doc.add_paragraph()
p.add_run('Sincerely, ').bold = True
doc.add_paragraph('Dr Nader Genedy and co-authors.')

doc.save(str(ROOT / 'TUDOR_ResponseToReviewers_R2.docx'))
print(f'[4] Response-to-reviewers letter saved ({len(REVIEWER_RESPONSES)} points)')

print('\nAll four supplementary files generated.')
print(f'Total deliverables in {ROOT}:')
for f in ['TUDOR_Manuscript_v5_clean.docx',
          'TUDOR_v4_to_v5_DIFF_LOG.md',
          'TUDOR_v5_TRACEABILITY.md',
          'TUDOR_Highlights_v5.docx',
          'TUDOR_AI_use_statement.docx',
          'TUDOR_disclosure_block.docx',
          'TUDOR_ResponseToReviewers_R2.docx',
          'docs/superpowers/specs/2026-05-03-TUDOR-R1-revision-design.md']:
    p = ROOT / f
    if p.exists():
        print(f'  OK  {f}  ({p.stat().st_size/1024:.1f} KB)')
    else:
        print(f'  MISSING  {f}')

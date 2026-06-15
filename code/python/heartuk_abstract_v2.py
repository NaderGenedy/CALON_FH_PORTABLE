#!/usr/bin/env python3
"""Generate Heart UK Abstract + Table as Word document — CALON-FH SSS."""

from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
import sys, os

doc = Document()

# ── Page setup ──
for section in doc.sections:
    section.top_margin = Cm(2.54)
    section.bottom_margin = Cm(2.54)
    section.left_margin = Cm(2.54)
    section.right_margin = Cm(2.54)

style = doc.styles['Normal']
font = style.font
font.name = 'Arial'
font.size = Pt(11)
font.color.rgb = RGBColor(0, 0, 0)

# ── Helpers ──
def add_heading_custom(text, size=14, color=(0, 51, 102)):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    p.space_before = Pt(12)
    p.space_after = Pt(6)
    run = p.add_run(text)
    run.bold = True
    run.font.size = Pt(size)
    run.font.name = 'Arial'
    run.font.color.rgb = RGBColor(*color)
    return p

def add_body(text, bold=False, italic=False, size=11, space_after=6):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.space_after = Pt(space_after)
    p.paragraph_format.line_spacing = 1.15
    run = p.add_run(text)
    run.font.size = Pt(size)
    run.font.name = 'Arial'
    run.bold = bold
    run.italic = italic
    return p

def add_mixed(parts, alignment=WD_ALIGN_PARAGRAPH.JUSTIFY, space_after=6):
    p = doc.add_paragraph()
    p.alignment = alignment
    p.space_after = Pt(space_after)
    p.paragraph_format.line_spacing = 1.15
    for text, bold, italic in parts:
        run = p.add_run(text)
        run.font.size = Pt(11)
        run.font.name = 'Arial'
        run.bold = bold
        run.italic = italic
    return p

def set_cell_shading(cell, color):
    shading = cell._tc.get_or_add_tcPr()
    elm = shading.makeelement(qn('w:shd'), {
        qn('w:fill'): color, qn('w:val'): 'clear', qn('w:color'): 'auto'
    })
    shading.append(elm)

def format_cell(cell, text, bold=False, size=9, alignment=WD_ALIGN_PARAGRAPH.CENTER):
    cell.text = ''
    p = cell.paragraphs[0]
    p.alignment = alignment
    p.space_before = Pt(2)
    p.space_after = Pt(2)
    run = p.add_run(text)
    run.font.size = Pt(size)
    run.font.name = 'Arial'
    run.bold = bold
    return run

# ═══════════════════════════════════════════════════════════
# TITLE
# ═══════════════════════════════════════════════════════════
add_heading_custom(
    'CALON-FH: AlphaFold3-Derived Structural Severity Score Enables '
    'Domain-Specific ASCVD Penetrance Prediction and Mutation-Guided '
    'Treatment Selection in Familial Hypercholesterolaemia',
    size=14, color=(0, 51, 102)
)

add_body(
    'Genedy N, [co-authors], Sherwood K, AbiFadel M, Sherwood R, Watts GF',
    italic=True, size=10, space_after=4
)

add_body(
    'Cardiff University School of Medicine; Wales FH Registry (DRAGON3); '
    'UK Biobank Research Analysis Platform',
    italic=True, size=9, space_after=12
)

# ═══════════════════════════════════════════════════════════
# ABSTRACT
# ═══════════════════════════════════════════════════════════
add_heading_custom('Abstract', size=13, color=(0, 51, 102))

# Background
add_mixed([
    ('Background: ', True, False),
    ('Current familial hypercholesterolaemia (FH) management applies uniform treatment '
     'algorithms regardless of mutation location, despite the LDL receptor (LDLR) comprising '
     'functionally distinct domains with different mechanisms of dysfunction. The ligand-binding '
     'repeats (LA1\u2013LA7) capture LDL particles, the EGF-A domain mediates PCSK9 binding, '
     'the beta-propeller enables receptor recycling at endosomal pH, and the cytoplasmic NPXY '
     'motif signals clathrin-mediated endocytosis. Disruption of each mechanism has distinct '
     'therapeutic implications, yet no framework exists to translate mutation location into '
     'treatment selection. We developed the Structural Severity Score (SSS) to classify LDLR '
     'mutations by mechanism of receptor dysfunction and predict domain-specific ASCVD '
     'penetrance and treatment response.', False, False)
])

# Methods
add_mixed([
    ('Methods: ', True, False),
    ('We constructed the SSS by integrating 14 independent data sources across all 860 LDLR '
     'residue positions: AlphaFold3 structure prediction (pLDDT confidence), FoldX 5.0 '
     'thermodynamic stability (\u0394\u0394G) for 359 clinical variants and saturation '
     'mutagenesis, Rosetta energy decomposition, AlphaFold3 mutant LDLR\u2013PCSK9 complex '
     'modelling, interface distance to the LDL binding surface derived from the Reimund ', False, False),
    ('et al.', False, True),
    (' (', False, False),
    ('Nature', False, True),
    (' 2025) cryo-EM structure, and Islam ', False, False),
    ('et al.', False, True),
    (' (2024) cell-surface LDLR functional assay (234 variants, 200 positions). ClinVar '
     'pathogenicity classifications were incorporated. Clinical validation used three independent '
     'cohorts: the Wales FH Registry (n=7,253; 2,405 genetically positive; 323 unique LDLR '
     'mutations; 1,515 ASCVD events), the DRAGON3 clinical cohort (n=1,362; 424 with identified '
     'mutations), and UK Biobank whole-exome sequencing (n=426,732; 1,623 FH carriers). '
     'Kaplan\u2013Meier survival analysis with time-to-first-ASCVD methodology was used for '
     'domain-level penetrance estimation. ASCVD prediction used logistic regression adjusted for '
     'age, sex, LDL-C, statin use, diabetes, smoking, hypertension, and tendon xanthomata. An '
     'interactive clinical atlas with 3D AlphaFold structure visualisation was deployed at '
     'https://nadergenedy.github.io/calon-fh-atlas/.', False, False)
])

# Results
add_mixed([
    ('Results: ', True, False),
    ('Domain-level SSS correlated with ASCVD penetrance across 19 LDLR domains (Spearman '
     '\u03c1=0.644, P=0.003), with a 5.7-fold range: ligand-binding R4 (21.4%, n=42) versus '
     'ligand-binding R3 (3.8%, n=26). SSS was significantly higher in patients with tendon '
     'xanthomata (0.343 vs 0.297, P=2.7\u00d710\u207b\u2074). In 1,095 patients with '
     'longitudinal LDL-C data, treatment response varied significantly by domain: ligand-binding '
     'R4 achieved 25.6% LDL-C reduction (n=21), EGF-A 21.9% (n=31), versus cytoplasmic domain '
     '4.1% (n=19) \u2014 consistent with the mechanistic prediction that NPXY internalisation '
     'defects are refractory to conventional lipid-lowering therapy. Twenty-eight patients on '
     'PCSK9 inhibitors spanning 14 LDLR domains and 23 unique variants achieved 52% mean LDL-C '
     'reduction, demonstrating efficacy across all mutation types through wild-type allele '
     'rescue in heterozygous FH. VUS carriers (n=232) demonstrated ASCVD rates comparable to '
     'pathogenic carriers (18.1% vs 15.6% vs 3.5% in mutation-negative). FoldX thermodynamic '
     'predictions were independently validated against the Islam ', False, False),
    ('et al.', False, True),
    (' functional assay: variants classified as decreased binding/uptake had mean '
     '\u0394\u0394G=3.37 kcal/mol versus 0.06 kcal/mol for undetermined variants (56-fold '
     'difference). A fully adjusted ASCVD prediction model incorporating SSS with conventional '
     'risk factors and xanthomata achieved AUC=0.85 (n=350, 39 events).', False, False)
])

# Conclusions
add_mixed([
    ('Conclusions: ', True, False),
    ('CALON-FH demonstrates that LDLR domain architecture determines both ASCVD penetrance '
     'and treatment response in heterozygous FH. The 5.7-fold variation in domain-specific '
     'ASCVD risk and the 6.2-fold variation in treatment response (25.6% vs 4.1%) support a '
     'paradigm shift from uniform genotype-aware management to structure-informed, '
     'mutation-guided precision treatment. PCSK9 inhibitor efficacy across all 14 tested '
     'domains provides mechanistic reassurance for early use through wild-type allele '
     'protection, while the identification of cytoplasmic NPXY mutations as treatment-refractory '
     'supports early referral for LDLR-independent therapies including evinacumab.', False, False)
])

# ═══════════════════════════════════════════════════════════
# PAGE BREAK — TABLE
# ═══════════════════════════════════════════════════════════
doc.add_page_break()

add_heading_custom(
    'Table 1. Domain-Specific ASCVD Penetrance, Treatment Response, and Structural '
    'Classification Across 3,064 LDLR Mutation Carriers',
    size=11, color=(0, 51, 102)
)

headers = [
    'LDLR Domain', 'n', 'ASCVD\nRate (%)', 'LDL-C\nReduction\non LLT (%)',
    'Mean\nSSS', 'Mechanistic\nClass', 'Recommended\nApproach'
]

rows_data = [
    ['Ligand-binding R4', '42', '21.4', '25.6', '0.303', 'Binding failure\n(BS1)', 'Early PCSK9i +\nezetimibe'],
    ['Transmembrane', '37', '21.6', '\u2014', '0.629', 'Trafficking\ndefect', 'High-intensity\nstatin + ezetimibe'],
    ['Cytoplasmic\n(NPXY)', '48', '18.8', '4.1', '0.447', 'Internalisation\ndefect', 'Evinacumab /\napheresis'],
    ['EGF-A', '102', '15.7', '21.9', '0.397', 'Dual: BS1 +\nPCSK9 binding', 'Statin + ezetimibe\n\u2192 monitor PCSK9i'],
    ['Ligand-binding\n(null)\u2020', '328', '14.3', '20.3', '0.539', 'Loss-of-\nfunction', 'PCSK9i\n(WT rescue)'],
    ['Linker\n(EGF-C/O-linked)', '161', '13.0', '12.2', '0.330', 'O-glycosylation\ndefect', 'High-intensity\nstatin'],
    ['EGF-B', '158', '12.7', '15.2', '0.274', 'Recycling\npivot', 'Standard\nescalation'],
    ['EGF-precursor\nhomology', '258', '10.9', '15.6', '0.409', 'Recycling /\nstructural', 'Standard\nescalation'],
    ['Beta-propeller', '894', '10.7', '11.6', '0.132', 'Recycling\n(BS2)', 'Standard\nescalation'],
    ['Ligand-binding R2', '214', '8.4', '16.9', '0.191', 'Binding\n(BS1)', 'Statin +\nezetimibe'],
    ['Ligand-binding R3', '26', '3.8', '25.3', '0.314', 'Binding\n(BS1)', 'Standard\ncare'],
]

col_widths = [Cm(2.8), Cm(1.0), Cm(1.3), Cm(1.5), Cm(1.0), Cm(2.2), Cm(2.5)]

table = doc.add_table(rows=1 + len(rows_data), cols=7)
table.alignment = WD_TABLE_ALIGNMENT.CENTER
table.style = 'Table Grid'

# Remove all borders
for row in table.rows:
    for cell in row.cells:
        tc = cell._tc
        tcPr = tc.get_or_add_tcPr()
        borders = tcPr.makeelement(qn('w:tcBorders'), {})
        for bn in ['top', 'left', 'bottom', 'right']:
            b = borders.makeelement(qn(f'w:{bn}'), {qn('w:val'): 'nil'})
            borders.append(b)
        tcPr.append(borders)

# Header row
header_row = table.rows[0]
for i, (h, w) in enumerate(zip(headers, col_widths)):
    cell = header_row.cells[i]
    cell.width = w
    format_cell(cell, h, bold=True, size=9)
    set_cell_shading(cell, '1F4E79')
    for p in cell.paragraphs:
        for run in p.runs:
            run.font.color.rgb = RGBColor(255, 255, 255)

# Header borders
for cell in header_row.cells:
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    borders = tcPr.makeelement(qn('w:tcBorders'), {})
    for bn in ['top', 'bottom']:
        b = borders.makeelement(qn(f'w:{bn}'), {
            qn('w:val'): 'single', qn('w:sz'): '12', qn('w:color'): '1F4E79'
        })
        borders.append(b)
    tcPr.append(borders)

# Data rows
for r_idx, rd in enumerate(rows_data):
    row = table.rows[r_idx + 1]
    for c_idx, (val, w) in enumerate(zip(rd, col_widths)):
        cell = row.cells[c_idx]
        cell.width = w
        align = WD_ALIGN_PARAGRAPH.LEFT if c_idx in [0, 5, 6] else WD_ALIGN_PARAGRAPH.CENTER
        format_cell(cell, val, bold=(c_idx == 0), size=9, alignment=align)

    if r_idx % 2 == 0:
        for cell in row.cells:
            set_cell_shading(cell, 'E8F0FE')

    # Traffic-light ASCVD column
    ascvd_val = rd[2]
    if ascvd_val != '\u2014':
        ascvd = float(ascvd_val)
        if ascvd >= 18:
            set_cell_shading(row.cells[2], 'FFCCCC')
        elif ascvd >= 12:
            set_cell_shading(row.cells[2], 'FFE5CC')
        else:
            set_cell_shading(row.cells[2], 'CCFFCC')

# Bottom border
for cell in table.rows[-1].cells:
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    borders = tcPr.makeelement(qn('w:tcBorders'), {})
    b = borders.makeelement(qn('w:bottom'), {
        qn('w:val'): 'single', qn('w:sz'): '12', qn('w:color'): '1F4E79'
    })
    borders.append(b)
    tcPr.append(borders)

# Table footnote
add_body(
    'Data from Wales FH Registry (n=3,064 LDLR carriers; 1,095 with longitudinal LDL-C). '
    'ASCVD = composite of MI, ACS, PCI, CABG, angina, TIA, PVD. LDL-C reduction = percentage '
    'change between first and second clinic visit on lipid-lowering therapy. SSS = Structural '
    'Severity Score integrating 14 data sources (AlphaFold3, FoldX 5.0, Rosetta, Islam et al. '
    'functional assay, UK Biobank WES outcomes, ClinVar). BS1 = primary LDL binding site '
    '(LA3\u2013LA7 + EGF-A, Reimund et al. Nature 2025); BS2 = secondary binding site '
    '(beta-propeller). WT = wild-type. ASCVD column colour-coded: red \u226518%, '
    'amber 12\u201317.9%, green <12%. '
    '\u2020Splice-site and frameshift mutations assigned SSS based on predicted '
    'loss-of-function. \u2014 = insufficient longitudinal data.',
    size=8, italic=True, space_after=12
)

# ═══════════════════════════════════════════════════════════
# SAVE
# ═══════════════════════════════════════════════════════════
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CALON_FH_HeartUK_Abstract.docx')
doc.save(out)
print(f'Saved: {out}')
print(f'Size: {os.path.getsize(out):,} bytes')

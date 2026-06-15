"""
Convert manuscript markdown files to Word .docx for submission.
Uses python-docx (no pandoc dependency).
Cambria 11pt body, 1.5 spacing, Cardiff-academic style.
"""
import re
import sys
from pathlib import Path
from docx import Document
from docx.shared import Pt, Cm, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

MANUSCRIPT_DIR = Path(r'C:/Users/nader/Downloads/calon_ukb_pipeline/manuscript')

def set_font(run, font='Cambria', size_pt=11, bold=False, italic=False):
    run.font.name = font
    rPr = run._element.get_or_add_rPr()
    rFonts = rPr.find(qn('w:rFonts'))
    if rFonts is None:
        rFonts = OxmlElement('w:rFonts')
        rPr.append(rFonts)
    rFonts.set(qn('w:ascii'), font)
    rFonts.set(qn('w:hAnsi'), font)
    run.font.size = Pt(size_pt)
    run.bold = bold
    run.italic = italic

def add_para(doc, text, style='Normal', bold=False, italic=False, size=11, align=None):
    p = doc.add_paragraph()
    if align is not None:
        p.alignment = align
    pf = p.paragraph_format
    pf.line_spacing = 1.5
    pf.space_after = Pt(6)
    run = p.add_run(text)
    set_font(run, size_pt=size, bold=bold, italic=italic)
    return p

def add_heading(doc, text, level=1):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.line_spacing = 1.5
    pf.space_before = Pt(12)
    pf.space_after = Pt(6)
    sizes = {1: 16, 2: 13, 3: 12}
    run = p.add_run(text)
    set_font(run, font='Cambria', size_pt=sizes.get(level, 12), bold=True)
    return p

def parse_inline(p, text):
    """Handle **bold**, *italic*, and `code` in inline text."""
    # Split on **bold** and *italic* markers
    parts = re.split(r'(\*\*[^*]+?\*\*|\*[^*]+?\*|`[^`]+?`|\^[^^]+?\^|_[^_]+?_)', text)
    for part in parts:
        if not part:
            continue
        if part.startswith('**') and part.endswith('**'):
            r = p.add_run(part[2:-2])
            set_font(r, bold=True)
        elif part.startswith('*') and part.endswith('*'):
            r = p.add_run(part[1:-1])
            set_font(r, italic=True)
        elif part.startswith('`') and part.endswith('`'):
            r = p.add_run(part[1:-1])
            set_font(r, font='Consolas', size_pt=10)
        elif part.startswith('^') and part.endswith('^'):
            r = p.add_run(part[1:-1])
            set_font(r)
            r.font.superscript = True
        elif part.startswith('_') and part.endswith('_'):
            r = p.add_run(part[1:-1])
            set_font(r, italic=True)
        else:
            r = p.add_run(part)
            set_font(r)

def md_to_docx(md_path: Path, docx_path: Path):
    doc = Document()
    # Page setup: A4, 1-inch margins
    for section in doc.sections:
        section.top_margin = Inches(1)
        section.bottom_margin = Inches(1)
        section.left_margin = Inches(1)
        section.right_margin = Inches(1)
    # Style defaults
    styles = doc.styles
    normal = styles['Normal']
    normal.font.name = 'Cambria'
    normal.font.size = Pt(11)

    lines = md_path.read_text(encoding='utf-8').splitlines()
    in_table = False
    table_rows = []

    i = 0
    while i < len(lines):
        line = lines[i].rstrip()

        if not line:
            i += 1
            continue

        # Headings
        if line.startswith('# '):
            add_heading(doc, line[2:].strip(), level=1)
        elif line.startswith('## '):
            add_heading(doc, line[3:].strip(), level=2)
        elif line.startswith('### '):
            add_heading(doc, line[4:].strip(), level=3)
        elif line.startswith('---'):
            # horizontal rule -> small spacer
            add_para(doc, '')
        elif line.startswith('| ') and '|' in line[2:]:
            # Table — collect rows
            cells = [c.strip() for c in line.strip('|').split('|')]
            # Skip the separator row (---)
            if not all(re.match(r'^-+$', c) for c in cells):
                table_rows.append(cells)
            # Look ahead — if next non-empty line is not a table row, render
            j = i + 1
            while j < len(lines) and lines[j].rstrip() == '':
                j += 1
            if j >= len(lines) or not lines[j].lstrip().startswith('|'):
                if table_rows:
                    t = doc.add_table(rows=len(table_rows), cols=len(table_rows[0]))
                    t.style = 'Light Grid Accent 1'
                    for r_idx, row in enumerate(table_rows):
                        for c_idx, cell_text in enumerate(row):
                            if c_idx < len(t.rows[r_idx].cells):
                                cell = t.rows[r_idx].cells[c_idx]
                                cell.text = ''
                                p = cell.paragraphs[0]
                                pf = p.paragraph_format
                                pf.space_after = Pt(0)
                                run = p.add_run(cell_text)
                                set_font(run, size_pt=10, bold=(r_idx == 0))
                    table_rows = []
                    add_para(doc, '')
        else:
            # Regular paragraph; handle inline markdown
            p = doc.add_paragraph()
            pf = p.paragraph_format
            pf.line_spacing = 1.5
            pf.space_after = Pt(6)
            pf.first_line_indent = Cm(0)
            parse_inline(p, line)
        i += 1

    doc.save(docx_path)
    print(f'  wrote {docx_path.name}')

if __name__ == '__main__':
    targets = [
        ('CALON_FH_manuscript_v2_Circulation.md', 'CALON_FH_manuscript_v2_Circulation.docx'),
        ('cover_letter_Circulation.md',           'cover_letter_Circulation.docx'),
        ('reviewer_anticipation_QA.md',           'reviewer_anticipation_QA.docx'),
    ]
    for md_name, docx_name in targets:
        md = MANUSCRIPT_DIR / md_name
        if not md.exists():
            print(f'  SKIP {md_name} (not found)')
            continue
        md_to_docx(md, MANUSCRIPT_DIR / docx_name)
    print('\nDONE')

"""
HEART UK Abstract - TUDOR Algorithm
Generates formatted Word document with embedded figures and table
"""

from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
import os

doc = Document()

# ---- Page Setup: A4 with narrow margins ----
section = doc.sections[0]
section.page_width = Cm(21)
section.page_height = Cm(29.7)
section.top_margin = Cm(1.0)
section.bottom_margin = Cm(0.8)
section.left_margin = Cm(1.5)
section.right_margin = Cm(1.5)

# ---- Helpers ----
def set_run(run, size=11, bold=False, italic=False, font_name="Calibri", superscript=False):
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.name = font_name
    if superscript:
        run.font.superscript = True
    rPr = run._element.get_or_add_rPr()
    rFonts = rPr.find(qn('w:rFonts'))
    if rFonts is None:
        rFonts = rPr.makeelement(qn('w:rFonts'), {})
        rPr.insert(0, rFonts)
    rFonts.set(qn('w:ascii'), font_name)
    rFonts.set(qn('w:hAnsi'), font_name)
    rFonts.set(qn('w:cs'), font_name)
    return run

def add_spacing(para, before=0, after=0, line=None):
    pPr = para._element.get_or_add_pPr()
    spacing = pPr.find(qn('w:spacing'))
    if spacing is None:
        spacing = pPr.makeelement(qn('w:spacing'), {})
        pPr.append(spacing)
    if before:
        spacing.set(qn('w:before'), str(before))
    if after:
        spacing.set(qn('w:after'), str(after))
    if line:
        spacing.set(qn('w:line'), str(line))
        spacing.set(qn('w:lineRule'), 'auto')

# ================================================================
# TITLE
# ================================================================
title_text = ("TUDOR: a treatment-adjusted diagnostic algorithm for familial hypercholesterolaemia "
              "with dual external validation in 4,028 genetically confirmed cases")
print(f"Title: {len(title_text)} chars (max 250)")

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.LEFT
add_spacing(p, after=60)
r = p.add_run(title_text)
set_run(r, size=11, bold=True)

# ================================================================
# AUTHORS
# ================================================================
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.LEFT
add_spacing(p, after=30)

r = p.add_run("N. Genedy*")
set_run(r, size=11)
r = p.add_run("1")
set_run(r, size=9, superscript=True)
r = p.add_run(", S. Zouwail")
set_run(r, size=11)
r = p.add_run("1")
set_run(r, size=9, superscript=True)

# ================================================================
# AFFILIATIONS
# ================================================================
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.LEFT
add_spacing(p, after=100)

r = p.add_run("1")
set_run(r, size=8, superscript=True, italic=True)
r = p.add_run("Department of Metabolic Medicine, Cardiff and Vale University Health Board, Cardiff, UK")
set_run(r, size=9, italic=True)

# ================================================================
# ABSTRACT BODY
# ================================================================

bg = ("Familial hypercholesterolaemia (FH) remains underdiagnosed worldwide, with fewer than 10% "
      "of cases identified. Established criteria (DLCN, Simon Broome) were calibrated before "
      "widespread statin use and apply uniform treatment correction, introducing systematic error. "
      "Cascade-screened relatives are missed because diagnostic thresholds were set on index cases "
      "with extreme phenotypes.")

methods = ("TUDOR (Treatment-adjusted Universal Detection and Outcome Risk) is a 10-feature elastic "
           "net logistic regression incorporating three innovations: (1) drug-specific dose-adjusted "
           "LDL-C, (2) an Index Effect encoding ascertainment pathway, and (3) a triglyceride "
           "metabolic filter (Trig_Filter = adjusted-LDL-C/TG) exploiting the metabolic selectivity "
           "of monogenic FH. Development used Cardiff index cases. First external validation: All "
           "Wales FH Registry (n=7,253; 2,405 FH; TRIPOD 2b). Second: UK Biobank lipid clinic cohort "
           "(n=58,021; 729 FH; TRIPOD Type 4). Gene-specific performance was assessed across LDLR "
           "and APOB carriers. ApoB augmentation models were evaluated.")

results = ("TUDOR achieved AUC 0.842 (95% CI 0.822\u20130.863) in Wales and 0.750 (0.731\u20130.770) in "
           "UK Biobank (Figure 1), significantly outperforming eDLCN (0.636; DeLong p=6.7\u00d710\u207b\u00b2\u2074) "
           "and isolated LDL-C (0.689). Gene-specific AUC was consistent: LDLR 0.839 and APOB 0.841 "
           "in Wales; LDLR 0.717 and APOB 0.830 in UKB (Figure 2). The ApoB/LDL-C ratio augmented "
           "TUDOR to AUC 0.771 (+0.019; p=4\u00d710\u207b\u2075). DLCN sensitivity in cascade relatives was "
           "1.5% versus TUDOR 89.4% (Table 1). Total: 4,028 genetically confirmed FH cases across "
           "both cohorts.")

conclusions = ("TUDOR is the first FH diagnostic algorithm incorporating individualised treatment "
               "adjustment, ascertainment bias correction, and metabolic filtering, with dual external "
               "validation across specialist and biobank settings. It addresses the 87-percentage-point "
               "sensitivity gap in cascade screening and supports optional ApoB augmentation.")

body_total = bg + " " + methods + " " + results + " " + conclusions
print(f"Body: {len(body_total)} chars (max 2500)")

# Background
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
add_spacing(p, after=40, line=264)
r = p.add_run("Background: ")
set_run(r, size=11, bold=True)
r = p.add_run(bg)
set_run(r, size=11)

# Methods
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
add_spacing(p, after=40, line=264)
r = p.add_run("Methods: ")
set_run(r, size=11, bold=True)
r = p.add_run(methods)
set_run(r, size=11)

# Results
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
add_spacing(p, after=40, line=264)
r = p.add_run("Results: ")
set_run(r, size=11, bold=True)
r = p.add_run(results)
set_run(r, size=11)

# Conclusions
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
add_spacing(p, after=80, line=264)
r = p.add_run("Conclusions: ")
set_run(r, size=11, bold=True)
r = p.add_run(conclusions)
set_run(r, size=11)

# ================================================================
# TABLE 1
# ================================================================
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
add_spacing(p, after=40)
r = p.add_run()
r.add_picture("C:/Users/nader/Downloads/calon_ukb_pipeline/tudor_table1.png",
              width=Cm(10.5))

# ================================================================
# FIGURE 1
# ================================================================
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
add_spacing(p, after=30)
r = p.add_run()
r.add_picture("C:/Users/nader/Downloads/calon_ukb_pipeline/tudor_fig1.png",
              width=Cm(10.5))

# ================================================================
# FIGURE 2
# ================================================================
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
add_spacing(p, after=10)
r = p.add_run()
r.add_picture("C:/Users/nader/Downloads/calon_ukb_pipeline/tudor_fig2.png",
              width=Cm(10.5))

# Save
output_path = "C:/Users/nader/Downloads/calon_ukb_pipeline/HEARTUK_Abstract_TUDOR.docx"
doc.save(output_path)
print(f"\nSaved: {output_path}")
print(f"File size: {os.path.getsize(output_path) / 1024:.1f} KB")

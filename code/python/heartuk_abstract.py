"""
HEART UK Abstract - CALON-Lite
Generates formatted Word document with embedded figures and table
"""

from docx import Document
from docx.shared import Pt, Cm, Inches, RGBColor, Emu
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

# ---- Helper to set font ----
def set_run(run, size=11, bold=False, italic=False, font_name="Calibri", color=None, superscript=False):
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.name = font_name
    if color:
        run.font.color.rgb = RGBColor(*color)
    if superscript:
        run.font.superscript = True
    # Set east asian and complex script fonts too
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
# TITLE (≤250 chars)
# ================================================================
title_text = ("CALON-Lite: a novel ASCVD risk score developed and validated exclusively in "
              "genetically confirmed familial hypercholesterolaemia using leave-one-cohort-out "
              "cross-validation")
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
r = p.add_run("2")
set_run(r, size=9, superscript=True)
r = p.add_run(", S. Sherif")
set_run(r, size=11)
r = p.add_run("3")
set_run(r, size=9, superscript=True)

# ================================================================
# AFFILIATIONS
# ================================================================
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.LEFT
add_spacing(p, after=100)

affil_parts = [
    ("1", "Department of Cardiovascular Sciences, Cardiff University, Cardiff, UK "),
    ("2", "Department of Chemical Pathology, Swansea Bay UHB, Swansea, UK "),
    ("3", "Department of Cardiology, Swansea Bay UHB, Swansea, UK"),
]
for num, text in affil_parts:
    r = p.add_run(num)
    set_run(r, size=8, superscript=True, italic=True)
    r = p.add_run(text)
    set_run(r, size=9, italic=True)

# ================================================================
# ABSTRACT BODY (≤2500 chars)
# ================================================================

bg = ("Familial hypercholesterolaemia (FH) confers lifelong elevated ASCVD risk, yet no risk score "
      "has been developed and externally validated exclusively in genetically confirmed FH. Existing "
      "tools (SAFEHEART-RE, MFHS, CFHS) show substantial performance heterogeneity across populations, "
      "and general-population models systematically underestimate risk by ignoring monogenic pathophysiology.")

methods = ("CALON-Lite is an 8-predictor penalised logistic regression model (elastic net, \u03B1=0.1) "
           "trained across three independent genetically confirmed FH cohorts: South Wales FH Registry "
           "(n=418), UK Biobank FH subset (n=1,625), and Wales FH Registry (n=2,405), totalling 4,448 "
           "patients with 835 ASCVD events (18.8%). External validation used TRIPOD Type 2b "
           "leave-one-cohort-out cross-validation (LOCO-CV). Predictors: age, sex, LDL-C, inverse HDL-C, "
           "smoking, diabetes, hypertension, and Lp(a) >143 nmol/L. ApoB/LDL-C discordance was analysed "
           "in 1,788 FH versus 423,002 non-FH UK Biobank participants.")

results = ("CALON-Lite achieved C-statistics of 0.790 (95% CI 0.741\u20130.839), 0.722 (0.694\u20130.748), "
           "and 0.855 (0.837\u20130.872) across three LOCO-CV folds (Figure 1), outperforming the refitted "
           "SRE in 2/3 folds (DeLong p=0.040 and p<0.001). The strongest predictors were Lp(a) >143 nmol/L "
           "(OR 2.40), hypertension (OR 2.12), sex (OR 2.13), and diabetes (OR 2.09; Table 1). LDL-C "
           "carried a negative coefficient (OR 0.89), reflecting statin-treatment paradox. ApoB/LDL-C "
           "discordance (\u22650.31 g/mmol) conferred an OR of 3.57 for ASCVD in FH versus 2.88 in non-FH "
           "(Figure 2), with event rates of 44.3% versus 18.2% in discordant FH\u2014a 2.4-fold enrichment.")

conclusions = ("CALON-Lite is the first ASCVD risk model developed and validated exclusively in genetically "
               "confirmed FH using multi-cohort cross-validation. It demonstrates robust discrimination "
               "using 8 routine clinical predictors and reveals that ApoB/LDL-C discordance has amplified "
               "prognostic value in the FH context, identifying a high-risk subgroup amenable to intensified "
               "therapy.")

body_total = bg + " " + methods + " " + results + " " + conclusions
print(f"Body: {len(body_total)} chars (max 2500)")

# Background
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
add_spacing(p, after=40, line=264)  # ~1.1 spacing
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
# TABLE 1 (as image)
# ================================================================
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
add_spacing(p, after=40)
r = p.add_run()
r.add_picture("C:/Users/nader/Downloads/calon_ukb_pipeline/heartuk_table1.png",
              width=Cm(10.5))

# ================================================================
# FIGURE 1
# ================================================================
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
add_spacing(p, after=30)
r = p.add_run()
r.add_picture("C:/Users/nader/Downloads/calon_ukb_pipeline/heartuk_fig1.png",
              width=Cm(10.5))

# ================================================================
# FIGURE 2
# ================================================================
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
add_spacing(p, after=10)
r = p.add_run()
r.add_picture("C:/Users/nader/Downloads/calon_ukb_pipeline/heartuk_fig2.png",
              width=Cm(10.5))

# Save
output_path = "C:/Users/nader/Downloads/calon_ukb_pipeline/HEARTUK_Abstract_CALON_Lite.docx"
doc.save(output_path)
print(f"\nSaved: {output_path}")
print(f"File size: {os.path.getsize(output_path) / 1024:.1f} KB")

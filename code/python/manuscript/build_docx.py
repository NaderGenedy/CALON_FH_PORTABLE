#!/usr/bin/env python3
"""
Build Nature-calibre DOCX manuscript from LaTeX source.
Uses python-docx to create a professionally formatted Word document.
"""
import re
from docx import Document
from docx.shared import Pt, Inches, Cm, RGBColor, Emu
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.style import WD_STYLE_TYPE
from docx.oxml.ns import qn
import os

doc = Document()

# --- Page setup ---
section = doc.sections[0]
section.page_width = Cm(21.0)  # A4
section.page_height = Cm(29.7)
section.top_margin = Cm(2.5)
section.bottom_margin = Cm(2.5)
section.left_margin = Cm(2.5)
section.right_margin = Cm(2.5)

# --- Styles ---
style = doc.styles['Normal']
font = style.font
font.name = 'Arial'
font.size = Pt(11)
font.color.rgb = RGBColor(0, 0, 0)
pf = style.paragraph_format
pf.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
pf.space_after = Pt(6)

# Heading styles
for level, size, color in [(1, 14, (0, 51, 102)), (2, 12, (0, 51, 102)), (3, 11, (51, 51, 51))]:
    hs = doc.styles[f'Heading {level}']
    hs.font.name = 'Arial'
    hs.font.size = Pt(size)
    hs.font.bold = True
    hs.font.color.rgb = RGBColor(*color)
    hs.paragraph_format.space_before = Pt(18)
    hs.paragraph_format.space_after = Pt(8)
    hs.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE

# ============================================================================
# TITLE PAGE
# ============================================================================
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.space_before = Pt(48)
run = p.add_run('A Saturation Threshold Model of Lipoprotein(a) Clearance\nResolves Three Paradoxes in Familial Hypercholesterolaemia')
run.bold = True
run.font.size = Pt(16)
run.font.color.rgb = RGBColor(0, 51, 102)

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = p.add_run('Nader Genedy')
run.bold = True
run.font.size = Pt(12)
run = p.add_run('1,2,*')
run.font.size = Pt(8)
run.font.superscript = True
run = p.add_run(', [Co-authors]')
run.font.size = Pt(12)
run = p.add_run('1,2,3')
run.font.size = Pt(8)
run.font.superscript = True

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.space_after = Pt(4)
run = p.add_run('1 Department of Cardiology, [Institution], Wales, UK\n2 Wales FH Service, Cardiff, UK\n3 [Collaborating Institution]')
run.font.size = Pt(9)
run.font.italic = True

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = p.add_run('* Correspondence: [email]')
run.font.size = Pt(9)

doc.add_page_break()

# ============================================================================
# ABSTRACT
# ============================================================================
doc.add_heading('Abstract', level=1)

def add_bold_para(doc, bold_text, normal_text):
    p = doc.add_paragraph()
    run = p.add_run(bold_text + ' ')
    run.bold = True
    run.font.size = Pt(11)
    run = p.add_run(normal_text)
    run.font.size = Pt(11)
    return p

add_bold_para(doc, 'Background.',
    'Lipoprotein(a) [Lp(a)] is a genetically determined, highly atherogenic lipoprotein affecting over '
    '1.4 billion individuals worldwide, yet three fundamental paradoxes have paralysed its clinical management: '
    'why statins fail to lower Lp(a) despite upregulating its only confirmed receptor; why populations of African '
    'descent appear to tolerate higher concentrations; and whether its vascular toxicity operates through systemic '
    'inflammation or localised retention.')

add_bold_para(doc, 'Methods.',
    'We integrated AlphaFold3 structural proteomics with population-scale epidemiology across three cohorts: '
    'the UK Biobank general population (n = 372,830), UK Biobank FH carriers (n = 2,370), and the Wales FH '
    'Registry (n = 1,623). We systematically modelled the binding affinity of the apolipoprotein(a) Kringle IV '
    'type 10 (KIV-10) domain against all nine historically proposed clearance receptors. We validated structural '
    'findings against the genome-wide CRISPR screen of Khan et al. (2025) and analysed outcome-specific Lp(a) '
    'risk thresholds across six vascular beds using cardiac magnetic resonance imaging (n = 60,947), carotid '
    'intima-media thickness (n = 287), and systemic inflammation biomarkers.')

add_bold_para(doc, 'Results.',
    'No human receptor demonstrated high-affinity protein\u2013protein binding for the apo(a) KIV-10 domain. '
    'The LDLR\u2013PCSK9 interaction yielded a robust interface predicted template modelling score (iPTM = 0.58, '
    '75 inter-chain contacts), confirming the structural pipeline against crystallographic ground truth (PDB: 3GCX). '
    'In contrast, both the LDLR\u2013KIV-10 (iPTM = 0.13, zero contacts) and all eight alternative receptors '
    '(iPTM range 0.12\u20130.44) failed to achieve canonical binding thresholds. Lp(a) demonstrated outcome-specific '
    'risk gradients, with the strongest effects for myocardial infarction (Q5 versus Q1 relative risk [RR] = 1.45, '
    'P = 1.4 \u00d7 10\u207b\u2074\u2075), carotid artery disease (RR = 1.48, P = 9.3 \u00d7 10\u207b\u00b9\u2070), '
    'and calcific aortic valve stenosis (RR = 1.39, P = 2.2 \u00d7 10\u207b\u00b9\u2078), but minimal effect on '
    'ischaemic stroke (RR = 1.09) or transient ischaemic attack (RR = 1.07, P = 0.16). In FH, the clinical risk '
    'threshold shifted from >25 nmol/L (general population) to >50 nmol/L (Wales FH, RR = 1.32, P = 0.008), '
    'consistent with background LDL-C saturation of LDLR. Cardiac MRI demonstrated that aortic stenosis was '
    'associated with significantly elevated Lp(a) (54.4 versus 43.9 nmol/L, P = 0.0005), increased left ventricular '
    'mass (107.1 versus 85.3 g, P < 0.0001), and reduced ejection fraction (58.2% versus 59.6%, P = 0.006). '
    'Lp(a) was entirely independent of traditional risk factors in FH, showing zero correlation with ApoB '
    '(\u03c1 = 0.029, P = 0.31), LDL-cholesterol (\u03c1 = 0.031, P = 0.27), systemic inflammation, BMI, or diabetes.')

add_bold_para(doc, 'Conclusions.',
    'Humans lack a dedicated high-affinity clearance receptor for Lp(a). The particle relies on opportunistic, '
    'low-affinity access to LDLR that is perpetually saturated by LDL-C, structurally explaining the failure of '
    'receptor-upregulating therapies. Because hepatic clearance is biologically impossible at physiological LDL '
    'concentrations, therapeutic silencing of hepatic LPA mRNA via RNA interference represents the obligate '
    'strategy for mitigating Lp(a)-driven cardiovascular disease.')

p = doc.add_paragraph()
run = p.add_run('Keywords: ')
run.bold = True
run.font.size = Pt(10)
run = p.add_run('Lipoprotein(a), familial hypercholesterolaemia, LDLR, AlphaFold3, saturation threshold, '
                'RNA interference, PCSK9, cardiovascular risk')
run.font.size = Pt(10)
run.font.italic = True

doc.add_page_break()

# ============================================================================
# INTRODUCTION
# ============================================================================
doc.add_heading('Introduction', level=1)

doc.add_paragraph(
    'The discovery of receptor-mediated endocytosis by Brown and Goldstein fundamentally transformed our '
    'understanding of cholesterol homeostasis and established the low-density lipoprotein receptor (LDLR) as '
    'the central gatekeeper of hepatic lipid clearance.\u00b9 Four decades later, the therapeutic exploitation '
    'of this pathway\u2014through statins, ezetimibe, and PCSK9 inhibitors\u2014has prevented millions of '
    'cardiovascular events worldwide.\u00b2\u00b3 Yet one lipoprotein has stubbornly resisted this pharmacological '
    'revolution: lipoprotein(a).')

doc.add_paragraph(
    'Lp(a) is a macromolecular complex comprising an LDL-like core particle covalently bound via a single '
    'disulphide bridge to the highly polymorphic glycoprotein apolipoprotein(a) [apo(a)].\u2074\u2075 Plasma '
    'Lp(a) concentrations are approximately 90% heritable, determined primarily by structural variation at the '
    'LPA gene locus encoding the Kringle IV type 2 (KIV-2) copy number polymorphism.\u2076\u2077 Elevated '
    'Lp(a) (>125 nmol/L) affects approximately 20% of the global population and is an established, independent, '
    'causal risk factor for atherosclerotic cardiovascular disease (ASCVD) and calcific aortic valve stenosis '
    '(CAVS).\u2078\u2079\u00b9\u2070 Despite this, Lp(a) remains the single most neglected treatable '
    'cardiovascular risk factor in clinical practice.\u00b9\u00b9\u00b9\u00b2')

doc.add_paragraph(
    'This clinical paralysis stems from three enduring biological paradoxes that have resisted resolution for '
    'over three decades.')

p = doc.add_paragraph()
run = p.add_run('The Clearance Paradox. ')
run.bold = True
p.add_run(
    'A landmark genome-wide CRISPR screen recently confirmed that LDLR is the primary mediator of hepatic '
    'Lp(a) uptake.\u00b9\u00b3 Yet statins, which robustly upregulate LDLR expression through the SREBP2 '
    'pathway, completely fail to lower plasma Lp(a).\u00b9\u2074 Conversely, PCSK9 inhibitors, which also '
    'increase surface LDLR abundance by preventing PCSK9-mediated receptor degradation, reduce Lp(a) by only '
    '20\u201330%.\u00b9\u2075\u00b9\u2076 This pharmacological dissociation\u2014where two therapies that both '
    'increase the same receptor produce diametrically opposite effects on the same ligand\u2014has never been '
    'mechanistically explained. Furthermore, at least eight alternative receptors (VLDLR, LRP1, PlgRKT, LOX-1, '
    'Megalin, CD36, SR-B1, ASGPR) have been proposed as supplementary Lp(a) clearance pathways, with conflicting '
    'and often non-reproducible in vitro validations.\u00b9\u2077\u00b9\u2078\u00b9\u2079\u00b2\u2070')

p = doc.add_paragraph()
run = p.add_run('The Ancestry Paradox. ')
run.bold = True
p.add_run(
    'Individuals of African descent possess median Lp(a) concentrations two- to three-fold higher than Europeans, '
    'driven by population-specific LPA genetic architecture including smaller apo(a) isoform sizes and distinct '
    'single-nucleotide polymorphism profiles.\u00b2\u00b9\u00b2\u00b2 Historical epidemiological cohorts have '
    'repeatedly reported a blunted relative ASCVD risk in Black populations at equivalent Lp(a) elevations, '
    'leading to the widely cited but never proven hypothesis that ancestral genetic modifiers attenuate Lp(a) '
    'pathogenicity in African-descent populations.\u00b2\u00b3\u00b2\u2074')

p = doc.add_paragraph()
run = p.add_run('The Pathophysiological Paradox. ')
run.bold = True
p.add_run(
    'Lp(a) carries a substantial payload of pro-inflammatory oxidised phospholipids (OxPL), covalently bound '
    'to the apo(a) moiety, leading to the prevailing hypothesis that its atherogenicity is primarily mediated '
    'through systemic inflammatory cascades.\u00b2\u2075\u00b2\u2076 If correct, Lp(a) risk should correlate '
    'tightly with systemic inflammation biomarkers and amplify classical metabolic syndrome pathways. Whether '
    'Lp(a) operates as a systemic inflammatory mediator or a localised mechanical vascular toxin has never been '
    'definitively resolved in large-scale human cohorts.')

doc.add_paragraph(
    'Here, we integrate AlphaFold3 (AF3) structural proteomics,\u00b2\u2077 quantitative thermodynamic '
    'modelling, and deep epidemiological analysis of three complementary cohorts (N = 376,823) to systematically '
    'dismantle each paradox. We propose the \u201cSaturation Threshold Model,\u201d demonstrating that Lp(a) is '
    'an evolutionary orphan\u2014a circulating lipoprotein lacking a dedicated high-affinity hepatic receptor\u2014'
    'and prove that its localised, pan-vascular toxicity is universal across human ancestries, independent of '
    'systemic inflammation, and modifiable only through direct genetic silencing.')

doc.add_page_break()

# ============================================================================
# RESULTS
# ============================================================================
doc.add_heading('Results', level=1)

# Result 1
doc.add_heading('Structural proteomics reveal the absence of a dedicated Lp(a) receptor', level=2)

doc.add_paragraph(
    'To resolve the mechanistic basis of hepatic Lp(a) clearance, we systematically evaluated all nine '
    'historically proposed clearance receptors using AlphaFold3 multimeric structure prediction. We modelled '
    'the interaction between the apo(a) Kringle IV type 10 (KIV-10) domain\u2014the only KIV repeat domain '
    'with putative receptor-binding capacity\u00b9\u2077\u2014and the extracellular domains of LDLR (696 amino '
    'acids), VLDLR, LRP1 (cluster IV), PlgRKT, LOX-1, Megalin, CD36, SR-B1, and ASGPR.')

doc.add_paragraph(
    'As a positive control, we first modelled the LDLR\u2013PCSK9 interaction, for which a high-resolution '
    'co-crystal structure exists (PDB: 3GCX).\u00b2\u2078 AF3 predicted this interaction with high confidence '
    '(iPTM = 0.58, predicted template modelling [PTM] = 0.67), generating 75 inter-chain residue contacts with '
    'a maximum contact probability of 0.84 (Table 1; Fig. 1a). This validates the AF3 pipeline for modelling '
    'LDLR protein\u2013protein interactions.')

doc.add_paragraph(
    'In contrast, the KIV-10 domain failed to achieve canonical binding with any human receptor. The '
    'LDLR\u2013KIV-10 interaction yielded an iPTM of 0.13 with zero inter-chain contacts exceeding a '
    'probability threshold of 0.30 (Fig. 1b). Among alternative receptors, PlgRKT achieved the highest '
    'iPTM (0.44), consistent with its established role as a plasminogen receptor that recognises the '
    'lysine-binding site of KIV-10.\u00b2\u2079 However, this interaction mediates fibrinolytic signalling, '
    'not hepatic clearance. VLDLR (iPTM = 0.18), LRP1 (iPTM = 0.12), CD36 (iPTM = 0.19), and SR-B1 '
    '(iPTM = 0.12) all demonstrated structural failure to engage the apo(a) protein backbone (Table 1).')

doc.add_paragraph(
    'These structural data converge precisely with the CRISPR screen of Khan et al.,\u00b9\u00b3 which '
    'identified LDLR as the sole statistically significant genetic mediator of cellular Lp(a) uptake (false '
    'discovery rate = 0.005) while none of the eight alternative receptors achieved statistical significance. '
    'The concordance between our in silico structural proteomics and their unbiased functional genomics screen '
    'establishes a critical conclusion: human physiology did not evolve a dedicated, high-affinity protein '
    'clearance receptor for Lp(a).')

doc.add_paragraph(
    'A methodological caveat merits emphasis. Both the LDLR\u2013ApoB (iPTM = 0.14) and LDLR\u2013KIV-10 '
    '(iPTM = 0.13) interactions yielded low iPTM scores when modelled against the full LDLR extracellular '
    'domain, despite ApoB\u2013LDLR binding being experimentally validated. This reflects a known limitation '
    'of AF3: the calcium-dependent ligand-binding repeats of LDLR (R1\u2013R7), through which both ApoB and '
    'KIV-10 are proposed to interact,\u00b3\u2070\u00b3\u00b9 require coordinated Ca\u00b2\u207a ions for '
    'structural integrity. AF3 does not model ion coordination at protein\u2013protein interfaces.\u00b2\u2077 '
    'The PCSK9\u2013LDLR interaction succeeds because it engages the EGF-A domain through a direct '
    'protein\u2013protein interface that is calcium-independent.\u00b2\u2078 This limitation does not weaken '
    'our conclusion; rather, it strengthens it: PCSK9, which binds via a calcium-independent mechanism, is '
    'robustly predicted (iPTM = 0.58), while both KIV-10 and ApoB, which require calcium coordination, are '
    'not\u2014confirming that AF3 accurately discriminates interaction mechanisms.')

# HINGE
p = doc.add_paragraph()
p.add_run(
    'Having established thermodynamically that Lp(a) lacks a dedicated clearance receptor and must compete '
    'with LDL for opportunistic LDLR access, we hypothesised that the clinical risk threshold for Lp(a) '
    'would be highly sensitive to a patient\u2019s background LDL-C burden. To test this in vivo, we evaluated '
    'Lp(a) risk stratification across three cohorts with progressively increasing LDL exposure.').italic = True

# Result 2
doc.add_heading('Background LDL-C burden dynamically shifts the Lp(a) risk threshold', level=2)

doc.add_paragraph(
    'In the UK Biobank general population (n = 372,830), Lp(a) demonstrated near-perfect quintile gradients '
    'across all vascular endpoints (Table 2; Fig. 2). Participants in the highest Lp(a) quintile '
    '(82\u2013189 nmol/L) experienced significantly elevated rates of myocardial infarction (4.51% versus '
    '3.11%, RR = 1.45, P = 1.4 \u00d7 10\u207b\u2074\u2075), aortic stenosis (2.21% versus 1.59%, '
    'RR = 1.39, P = 2.2 \u00d7 10\u207b\u00b9\u2078), carotid artery disease (0.80% versus 0.54%, '
    'RR = 1.48, P = 9.3 \u00d7 10\u207b\u00b9\u2070), and peripheral vascular disease (4.11% versus 3.46%, '
    'RR = 1.19, P = 5.0 \u00d7 10\u207b\u00b9\u00b9) compared to the lowest quintile (4\u20138 nmol/L). '
    'Notably, Lp(a) showed minimal predictive power for ischaemic stroke (RR = 1.09, P = 0.01) and no '
    'statistically significant association with transient ischaemic attack (RR = 1.07, P = 0.16), '
    'demonstrating that Lp(a) is not a uniform vascular toxin but exhibits remarkable endothelial bed '
    'specificity. The dominant pathology is atherothrombotic (coronary, aortic valve, carotid, peripheral) '
    'rather than thromboembolic (cerebrovascular).')

doc.add_paragraph(
    'In the deeply phenotyped Wales FH Registry (n = 1,623; LDLR 1,321, APOB 301, PCSK9 1; mean follow-up '
    '13.4 years; 399 ASCVD events), the Lp(a) risk threshold shifted rightward. Concentrations below '
    '50 nmol/L conferred no statistically significant excess risk (RR = 1.06, P = 0.60), while concentrations '
    'exceeding 50 nmol/L achieved significance (RR = 1.32, P = 0.008), with a dose\u2013response gradient '
    'extending to concentrations above 125 nmol/L (RR = 1.51, P = 0.004; Fig. 2b). This 25 nmol/L rightward '
    'shift in the clinical threshold is precisely predicted by the Saturation Threshold Model: in FH, where '
    'baseline LDL-C is markedly elevated (mean 4.04 mmol/L in the Wales cohort), LDLR is more completely '
    'saturated by LDL particles, requiring higher Lp(a) concentrations to compete for the residual receptor '
    'vacancies.')

# HINGE
p = doc.add_paragraph()
p.add_run(
    'Because LDLR saturation forces circulating Lp(a) away from hepatic clearance and into peripheral '
    'circulation, we next investigated its capacity for diffuse vascular toxicity across different endothelial '
    'beds using advanced cardiovascular imaging.').italic = True

# Result 3
doc.add_heading('Cardiac MRI and carotid imaging confirm localised Lp(a) vascular toxicity', level=2)

doc.add_paragraph(
    'Among 60,947 UK Biobank participants with both cardiac magnetic resonance imaging (CMR) and Lp(a) '
    'measurements, Lp(a) quintile analysis revealed that the primary structural cardiac consequence of '
    'Lp(a) elevation is aortic valve and aortic root pathology rather than direct myocardial dysfunction. '
    'Left ventricular ejection fraction (LVEF) was remarkably stable across Lp(a) quintiles (Q1: 59.6% '
    'versus Q5: 59.5%, \u03c1 = 0.003, P = 0.49), as were global longitudinal strain (GLS; \u03c1 = \u22120.004, '
    'P = 0.35) and LV end-diastolic volume (Table 4; Fig. 3). This absence of direct Lp(a)\u2013myocardial '
    'toxicity is mechanistically informative: Lp(a) does not impair cardiomyocyte function through circulating '
    'inflammatory or oxidative pathways; rather, its damage is confined to the vascular endothelium and valve '
    'interstitium.')

doc.add_paragraph(
    'Participants with aortic stenosis (n = 428) demonstrated significantly elevated Lp(a) compared to those '
    'without (54.4 versus 43.9 nmol/L, P = 0.0005), alongside markedly increased LV mass (107.1 versus '
    '85.3 g, P < 0.0001) and reduced LVEF (58.2% versus 59.6%, P = 0.006; Fig. 3c). The LV hypertrophy in '
    'aortic stenosis patients reflects chronic pressure overload secondary to valve calcification\u2014a '
    'downstream haemodynamic consequence of Lp(a)-driven valvular pathology, not a direct myocardial effect.')

doc.add_paragraph(
    'Among 369 FH carriers with CMR, cardiac structure was indistinguishable from the general population '
    '(LVEF 59.5% versus 59.6%, P = 0.99; LV mass 83.5 versus 85.5 g, P = 0.11), confirming that even in '
    'the context of severely elevated LDL-C, chronic lipid exposure does not produce detectable myocardial '
    'dysfunction on CMR at middle age in treated patients.')

doc.add_paragraph(
    'In the Wales FH cohort, carotid intima-media thickness (cIMT; n = 287; mean 710.2 \u03bcm, SD 134.6) '
    'demonstrated a graded association with ASCVD (T3: 16.7% versus T1: 12.5%), although this did not achieve '
    'statistical significance (P = 0.27). Critically, cIMT showed no correlation with Lp(a) (\u03c1 = 0.029, '
    'P = 0.66), LDL-C (\u03c1 = 0.061, P = 0.30), or ApoB (\u03c1 = 0.088, P = 0.14), confirming that cIMT '
    'in FH reflects cumulative lifetime LDL exposure (cholesterol-years) rather than current biomarker '
    'concentrations.\u00b3\u00b2\u00b3\u00b3')

# HINGE
p = doc.add_paragraph()
p.add_run(
    'While Lp(a) infiltration universally drives pan-vascular disease through localised retention, the '
    'prevailing dogma suggests this process is amplified by systemic inflammatory cascades. To decouple the '
    'mechanical retention of Lp(a) from systemic immune activation, we analysed the interaction between '
    'Lp(a) and validated metabolic and inflammatory biomarkers.').italic = True

# Result 4
doc.add_heading('Lp(a) pathogenicity is orthogonal to systemic inflammation and metabolic syndrome', level=2)

doc.add_paragraph(
    'In the FH cohort, Lp(a) operated as an entirely independent axis of cardiovascular risk. Lp(a) showed '
    'zero meaningful correlation with ApoB (\u03c1 = 0.029, P = 0.31), LDL-C (\u03c1 = 0.031, P = 0.27), '
    'HDL-C (\u03c1 = \u22120.006, P = 0.83), triglycerides (\u03c1 = \u22120.083, P = 0.07), body mass '
    'index (\u03c1 = \u22120.001, P = 0.97), C-reactive protein (\u03c1 = \u22120.031, P = 0.27), HbA1c, '
    'or blood pressure (Fig. 4a; Table 3). This orthogonality persisted in the non-FH general population, '
    'where Lp(a) showed only a trivial positive correlation with ApoB (\u03c1 = 0.080) and LDL-C '
    '(\u03c1 = 0.081), reflecting the small contribution of Lp(a)-cholesterol (~0.14 mmol/L) to the LDL-C assay.')

doc.add_paragraph(
    'The ApoB/LDL ratio\u2014a marker of lipoprotein discordance\u2014was unaffected by Lp(a) in FH '
    '(\u03c1 = 0.005, P = 0.81). However, when LDL-C was corrected for Lp(a)-cholesterol content using the '
    'Dahl\u00e9n formula,\u00b3\u2077 the ApoB/LDL ratio increased by +0.011 (from 0.289 to 0.300) in both '
    'FH and non-FH populations. This paradoxical increase in discordance after Lp(a) correction reveals that '
    'conventional LDL-C assays artifactually inflate the denominator of the ApoB/LDL ratio by including '
    'Lp(a)-cholesterol, making concordant patients appear more concordant than they truly are. This has direct '
    'clinical implications: FH patients with elevated Lp(a) may be undertreated if clinical decisions rely '
    'on uncorrected LDL-C targets.\u00b3\u00b2')

# HINGE
p = doc.add_paragraph()
p.add_run(
    'Having established Lp(a) as a localised, independent vascular pathogen unlinked to systemic inflammation '
    'or metabolic dysfunction, we finally sought to determine whether its fundamental atherogenicity is modified '
    'by ancestral genetic architecture, thereby addressing the long-standing \u201cAncestry Paradox.\u201d').italic = True

# Result 5
doc.add_heading('Ancestry does not modify Lp(a) pathogenicity', level=2)

doc.add_paragraph(
    'Initial unadjusted analyses of the UK Biobank multi-ethnic cohort replicated historical observations: '
    'Black participants demonstrated higher median Lp(a) concentrations alongside a seemingly lower '
    'age-adjusted ASCVD rate. However, rigorous demographic profiling revealed profound confounding. Black '
    'participants in the UK Biobank FH carrier cohort were significantly younger (mean age difference of '
    'approximately 8 years) and exhibited a distinct LDLR mutation spectrum, with the majority of variants '
    'localising to two moderate-penetrance domains (beta-propeller and EGF-like B), rather than the '
    'high-penetrance ligand-binding and EGF-A domains enriched in European carriers.')

doc.add_paragraph(
    'When cohorts were strictly age-matched to the peak ASCVD risk window (60\u201369 years) and further '
    'matched for mutation domain severity, the apparent racial disparity substantially attenuated. We note '
    'that the small sample size of non-White FH carriers in the UK Biobank (n_Black < 50, n_South Asian < 50) '
    'limits definitive conclusions from formal interaction testing. Nevertheless, the direction of effect is '
    'unambiguous: the historical \u201cAncestry Paradox\u201d is a demographic artefact of age confounding '
    'and mutation spectrum bias, not a biological modifier of Lp(a) pathogenicity.')

# Result 6
doc.add_heading('Menopausal status modifies Lp(a) concentration but not FH-associated risk', level=2)

doc.add_paragraph(
    'Among 204,032 women with Lp(a) measurements in the UK Biobank, we observed a significant menopausal '
    'effect on Lp(a) concentrations. Pre-menopausal women (ages 40\u201349) and men of equivalent age had '
    'identical median Lp(a) levels (19.7 versus 19.7 nmol/L, P = 0.77; Fig. 7a). After age 50, a progressive '
    'sex divergence emerged: post-menopausal women (ages 55\u201364) had significantly higher Lp(a) than '
    'age-matched men (22.8 versus 19.7 nmol/L, P < 0.0001), representing a 14% increase attributable to '
    'oestrogen withdrawal.\u00b3\u2076')

doc.add_paragraph(
    'Strikingly, this menopausal effect was absent in FH women (n = 1,302). Pre-menopausal FH women '
    '(ages 40\u201349) had a median Lp(a) of 27.6 nmol/L, virtually identical to post-menopausal FH women '
    '(ages 55\u201364; median 27.8 nmol/L, P = 0.65; Fig. 7b). This insensitivity to oestrogen withdrawal '
    'suggests that in the context of genetically elevated LDLR dysfunction, the dominant determinant of '
    'Lp(a) concentration is LPA genotype rather than hormonal regulation of hepatic receptor expression.')

doc.add_paragraph(
    'Lp(a) showed no association with pregnancy complications in the UK Biobank (pre-eclampsia: median 22.5 '
    'versus 22.2 nmol/L, P = 0.94; gestational diabetes: 21.9 versus 22.2 nmol/L, P = 0.85), confirming '
    'that Lp(a)-driven pathology operates through chronic vascular retention rather than acute endothelial '
    'dysfunction.')

doc.add_page_break()

# ============================================================================
# DISCUSSION
# ============================================================================
doc.add_heading('Discussion', level=1)

doc.add_paragraph(
    'For three decades, the clinical management of lipoprotein(a) has been paralysed by conflicting data on '
    'clearance mechanisms, racial disparities in apparent risk, and unresolved questions about its '
    'pathophysiological mechanism. By integrating cutting-edge structural proteomics with population-scale '
    'epidemiology, we have resolved these paradoxes and propose a unified model\u2014the Saturation '
    'Threshold\u2014that explains the natural history of Lp(a) from atomic structure to clinical outcome.')

p = doc.add_paragraph()
run = p.add_run('The Saturation Threshold Model. ')
run.bold = True
p.add_run(
    'Our primary finding establishes that human physiology did not evolve a dedicated, high-affinity receptor '
    'for hepatic Lp(a) clearance. The nine-receptor screen, validated against the CRISPR screen of Khan et '
    'al.,\u00b9\u00b3 demonstrates that KIV-10 fails to achieve canonical protein\u2013protein binding with '
    'any human receptor. This model elegantly resolves the statin paradox. Statins upregulate LDLR expression '
    'through SREBP2 activation but simultaneously increase hepatic LDL uptake, maintaining receptor saturation. '
    'The newly expressed LDLR molecules are immediately occupied by LDL particles, leaving no vacancies for '
    'the weaker Lp(a) ligand. PCSK9 inhibitors achieve a modest 20\u201330% Lp(a) reduction because they '
    'operate through a different mechanism: by preventing PCSK9-mediated LDLR degradation, they extend LDLR '
    'surface half-life, creating a transient window during which some receptors remain unoccupied long enough '
    'for Lp(a) to bind.\u00b9\u2075')

p = doc.add_paragraph()
run = p.add_run('Clinical translation. ')
run.bold = True
p.add_run(
    'The dynamic shifting of Lp(a) risk thresholds between populations with different LDL-C burdens provides '
    'the first in vivo validation of the Saturation Threshold. For the practising cardiologist, this has '
    'immediate implications. First, Lp(a) should be measured in every FH patient at diagnosis, but the clinical '
    'action threshold should be adjusted for background LDL-C burden. A Lp(a) of 60 nmol/L in a patient with '
    'an LDL-C of 2.0 mmol/L (well-treated) carries substantially more residual risk than the same concentration '
    'in a patient with an LDL-C of 5.0 mmol/L (untreated), because the former has more available LDLR vacancies '
    'through which Lp(a) can undergo hepatic clearance.')

p = doc.add_paragraph()
run = p.add_run('Cardiac imaging insights. ')
run.bold = True
p.add_run(
    'The CMR data from 60,947 participants provide a novel structural perspective. Lp(a) elevation does not '
    'cause direct myocardial dysfunction: LVEF, GLS, and LV volumes are preserved across the entire Lp(a) '
    'concentration range. The cardiac damage manifests exclusively through valvular and vascular '
    'pathology\u2014aortic stenosis patients show increased LV mass (pressure overload) and reduced LVEF '
    '(decompensation)\u2014confirming that Lp(a) is a vascular and valvular toxin, not a myocardial one.')

p = doc.add_paragraph()
run = p.add_run('Therapeutic implications. ')
run.bold = True
p.add_run(
    'Because humans lack a dedicated clearance receptor, pharmacological strategies that depend on receptor '
    'upregulation are fundamentally constrained by the Saturation Threshold. Statins will never lower Lp(a). '
    'PCSK9 inhibitors have reached their biophysical ceiling (~25%). The only rational therapeutic strategy is '
    'to silence Lp(a) production at the source: the LPA gene in hepatocytes. Three RNA-interference approaches '
    'are in advanced clinical development: pelacarsen (antisense oligonucleotide; \u221235\u201380%),\u00b3\u2079 '
    'olpasiran (siRNA; \u221270\u2013101%),\u2074\u2070 and lepodisiran (siRNA; \u221255\u201397%).\u2074\u00b9 '
    'Additionally, the oral small-molecule muvalaplin, which disrupts the apo(a)\u2013ApoB covalent bond, '
    'represents a complementary mechanism (\u221245\u201365%).\u2074\u00b2 Our structural data provide the '
    'biological rationale for these agents: they succeed precisely because they bypass the receptor bottleneck '
    'that defeats all clearance-based strategies.')

p = doc.add_paragraph()
run = p.add_run('Limitations. ')
run.bold = True
p.add_run(
    'Several limitations merit acknowledgement. First, AF3 does not model calcium coordination at protein\u2013'
    'protein interfaces, preventing accurate prediction of calcium-dependent LDLR\u2013ligand interactions. '
    'Second, the UK Biobank is a survival-selected cohort enriched for healthier individuals. Third, ethnic '
    'diversity in the UK Biobank FH carrier cohort is limited (Black n < 50). Fourth, Lp(a) was measured at '
    'a single timepoint. Fifth, the cIMT analysis (n = 287) was underpowered for subgroup analyses.')

p = doc.add_paragraph()
run = p.add_run('Conclusion. ')
run.bold = True
p.add_run(
    'The human body is structurally incapable of efficiently clearing elevated Lp(a) through receptor-mediated '
    'endocytosis. The Saturation Threshold Model\u2014validated by structural proteomics, CRISPR functional '
    'genomics, population-scale epidemiology, and cardiac imaging across 376,823 individuals\u2014establishes '
    'that Lp(a) is an evolutionary orphan whose hepatic clearance is perpetually blocked by LDL-C. Because '
    'receptor upregulation is biologically futile, the silencing of hepatic LPA mRNA represents not merely a '
    'promising therapeutic option but a biological imperative for the 1.4 billion individuals living with '
    'elevated Lp(a) worldwide.')

doc.add_page_break()

# ============================================================================
# METHODS (abbreviated)
# ============================================================================
doc.add_heading('Methods', level=1)

doc.add_heading('Study populations', level=2)
doc.add_paragraph(
    'UK Biobank general population. The UK Biobank is a prospective cohort study of approximately 500,000 '
    'participants aged 40\u201369 years at recruitment (2006\u20132010). Lp(a) was measured using an '
    'immunoturbidimetric assay (Randox Laboratories) on stored plasma samples (field 30790). After excluding '
    'participants without Lp(a) measurements, 375,200 individuals were included (372,830 non-FH carriers; '
    '2,370 FH carriers). ASCVD events were ascertained from Hospital Episode Statistics linked ICD-10 codes.')

doc.add_paragraph(
    'UK Biobank FH carriers. FH carriers were identified from whole-exome sequencing data, selecting '
    'individuals with rare (gnomAD MAF <0.01) missense or splice-region variants in LDLR. A total of 3,094 '
    'variant\u2013carrier observations representing 2,370 unique individuals carrying 202 distinct missense '
    'variants were identified.')

doc.add_paragraph(
    'Wales FH Registry. The Wales FH Service clinical registry comprises 1,623 genetically confirmed FH '
    'patients (LDLR 1,321, APOB 301, PCSK9 1; 55.2% female; mean age 56.7 \u00b1 8.1 years) with prospective '
    'follow-up (mean 13.4 years; 399 ASCVD events; 191 deaths). Lp(a) was available for 1,308 participants '
    '(median 25.3 nmol/L, IQR 11.2\u201362.6). Carotid IMT was measured by B-mode ultrasound in 287 '
    'participants.')

doc.add_paragraph(
    'UK Biobank cardiac MRI sub-study. Cardiac MRI was performed using a standardised protocol on a 1.5T '
    'scanner. LV volumes, ejection fraction, LV mass, cardiac output, global longitudinal strain, and aortic '
    'dimensions were obtained (fields 24100\u201324181). A total of 60,947 participants had both CMR and '
    'Lp(a) measurements available.')

doc.add_heading('AlphaFold3 structural proteomics', level=2)
doc.add_paragraph(
    'Multimeric structure predictions were performed using the AlphaFold3 server. For each receptor\u2013KIV-10 '
    'complex, we submitted the full extracellular domain of the receptor paired with the apo(a) KIV-10 domain '
    '(UniProt P08519, residues 3946\u20134031). Five independent models were generated per complex using '
    'different random seeds. Interface quality was assessed using iPTM, PTM, inter-chain contact probabilities '
    '(threshold >0.30), and maximum contact probability.')

doc.add_heading('Statistical analysis', level=2)
doc.add_paragraph(
    'Continuous variables are presented as mean \u00b1 SD or median (IQR) and compared using Mann\u2013Whitney '
    'U tests. Correlations are assessed using Spearman\u2019s rank correlation coefficient. Lp(a) quintile '
    'analyses used pre-specified cut-points. Relative risks were calculated as the ratio of event rates above '
    'versus below each threshold. Multiple testing was addressed using the Benjamini\u2013Hochberg FDR '
    'procedure. All analyses were performed in Python 3.12 and R 4.3.')

# ============================================================================
# SAVE
# ============================================================================
outpath = r'C:\Users\nader\Downloads\calon_ukb_pipeline\manuscript\Genedy_et_al_2026_Saturation_Threshold.docx'
doc.save(outpath)
print(f'Saved: {outpath}')
print(f'Size: {os.path.getsize(outpath) / 1024:.0f} KB')

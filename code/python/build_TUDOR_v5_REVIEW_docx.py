"""Build TUDOR_v5_REVIEW.docx -- reviewer-friendly walkthrough."""
from docx import Document
from docx.shared import Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH

doc = Document()

style = doc.styles['Normal']
style.font.name = 'Cambria'
style.font.size = Pt(11)
pf = style.paragraph_format
pf.line_spacing = 1.15
pf.space_after = Pt(4)

INS = RGBColor(0x10, 0x6E, 0x10)
TAG = RGBColor(0x1F, 0x77, 0xB4)
WARN = RGBColor(0x80, 0x40, 0x10)


def H(text, level=1):
    h = doc.add_heading(text, level=level)
    for r in h.runs:
        r.font.name = 'Cambria'
    return h


def P(text, bold=False, italic=False, colour=None):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.font.name = 'Cambria'
    r.font.size = Pt(11)
    if bold:
        r.bold = True
    if italic:
        r.italic = True
    if colour:
        r.font.color.rgb = colour
    return p


def P_mixed(parts):
    p = doc.add_paragraph()
    for text, bold, italic, colour in parts:
        r = p.add_run(text)
        r.font.name = 'Cambria'
        r.font.size = Pt(11)
        if bold:
            r.bold = True
        if italic:
            r.italic = True
        if colour:
            r.font.color.rgb = colour
    return p


# Title
title = doc.add_paragraph()
title.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = title.add_run('TUDOR -- Manuscript v5 (clean) -- Review Walkthrough')
r.font.name = 'Cambria'
r.font.size = Pt(16)
r.bold = True

sub = doc.add_paragraph()
sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = sub.add_run('JCLINLIPID-D-25-01142 R1 -> R2  |  Generated 2026-05-12')
r.font.name = 'Cambria'
r.font.size = Pt(11)
r.italic = True
doc.add_paragraph()

# Front matter
H('Title, authors, keywords', 2)
P('TUDOR: Development and Dual External Validation of a Treatment-Adjusted Diagnostic Algorithm for Familial Hypercholesterolaemia in 4,028 Genetically Confirmed Cases', bold=True)
P('Nader Genedy, MBBCh, MRCPUK, MRCPI, SCE(AIM), MRCGPUK, PgD (1)*; Soha Zouwail, FRCPath, MD, PhD (1)')
P('(1) Department of Metabolic Medicine, Cardiff and Vale University Health Board, Cardiff, Wales, UK')
P('* Corresponding author: genedyn1@cardiff.ac.uk')
P('Keywords: familial hypercholesterolaemia; cascade screening; treatment-adjusted LDL-C; diagnostic algorithm; triglyceride metabolic filter; apolipoprotein B; UK Biobank; external validation; gradient boosted machine; genetic validation', italic=True)

# Highlights
H('Highlights (Elsevier; 4 bullets, 72-75 chars, all <= 85)', 2)
for h in [
    'TUDOR achieves AUC 0.842 in the All Wales FH Registry and 0.750 in UK Biobank validation (n = 58,021)',
    'Largest dual-validated genetically confirmed FH study: 4,028 cases across six causative genes',
    'Gene-specific validation: APOB AUC 0.830 vs LDLR 0.717 in UK Biobank (DeLong p = 1.88 x 10^-4)',
    'ApoB/LDL-C ratio augmentation improves discrimination to AUC 0.771 (p = 3.99 x 10^-5)',
    'Triglyceride metabolic filter independent of BMI, diabetes, and ApoB (OR stable: 1.87 -> 1.94)',
    'CV mortality gradient: IHD deaths 0.65% -> 1.76% across LDL deciles (OR = 1.036/decile, p = 0.011)',
    'T2DM ascertainment: three independent definitions (self-report, ICD-10, HbA1c >= 48 mmol/mol) confirm TRG attenuation (beta = -0.48 to -0.66)',
]:
    p = doc.add_paragraph(h, style='List Bullet')
    for r in p.runs:
        r.font.name = 'Cambria'
        r.font.size = Pt(11)

# Abstract
H('Abstract', 2)
P_mixed([
    ('Background. ', True, False, None),
    ('Familial hypercholesterolaemia (FH) remains profoundly underdiagnosed, with established criteria demonstrating critical limitations in the era of widespread lipid-lowering therapy. We developed and externally validated TUDOR (Treatment-adjusted Universal Detection and Outcome Risk), incorporating individualised treatment-intensity adjusted LDL-C, ascertainment bias correction, and a triglyceride-based metabolic filter. ', False, False, None),
    ('In plain terms, the Elastic Net is a regression that selects only the most informative clinical features and shrinks the rest to zero, producing a parsimonious diagnostic score that is easier to inspect, audit, and deploy than a black-box machine-learning model. ', False, True, INS),
    ('[R2.1 preface -- inserted]', True, False, TAG),
])
P_mixed([
    ('Methods. ', True, False, None),
    ('TUDOR was developed using Elastic Net logistic regression, then validated in the All Wales FH Registry (n = 7,253; 2,405 genetically confirmed FH; TRIPOD Type 2b) and a lipid clinic-mimicking UK Biobank cohort (n = 58,021; 729 FH cases; TRIPOD Type 4), applying fixed Wales-derived coefficients without re-estimation. Gene-specific performance, cardiovascular mortality gradient, and T2DM ascertainment sensitivity served as biological validation.', False, False, None),
])
P_mixed([
    ('Results. ', True, False, None),
    ('TUDOR achieved AUC 0.842 (95% CI: 0.822-0.863) in Wales and 0.750 (0.731-0.770) in UK Biobank, outperforming estimated DLCN (AUC 0.636; p = 6.73 x 10^-24) and isolated LDL-C (AUC 0.689). ApoB augmentation improved discrimination to AUC 0.771 (delta-AUC +0.019; p = 3.99 x 10^-5). DLCN sensitivity in cascade-screened relatives was 1.5% versus TUDOR 89.4%. IHD mortality rose from 0.65% to 1.76% across LDL deciles (OR = 1.036, p = 0.011). Three independent T2DM definitions (kappa = 0.28-0.36) each confirmed triglyceride filter attenuation (beta = -0.48 to -0.66), with HbA1c dose-response confirming a mechanistic gradient.', False, False, None),
])
P_mixed([
    ('Conclusions. ', True, False, None),
    ('TUDOR is a dual-externally-validated, interpretable FH diagnostic algorithm maintaining discrimination across specialist and population-based settings in 4,028 genetically confirmed cases across six genes.', False, False, None),
])

# Section 1
H('1. Introduction', 1)
H('1.1 Familial hypercholesterolaemia: pathophysiology and CV burden', 3)
P('Autosomal dominant ~1/250; LDLR, APOB, PCSK9, less commonly LDLRAP1 and APOE. ~50% functional receptor reduction, doubled LDL particle residence time, intact VLDL/TG clearance -- biological foundation of the Trig_Filter. Untreated heterozygous FH confers 10-20x premature CHD risk by age 50; <10% diagnosed globally.')

H('1.2 Limitations of established criteria', 3)
P('DLCN, Simon Broome, MEDPED all calibrated before contemporary lipid-lowering practice. Uniform 1.43 correction factor introduces systematic error because efficacy spans 15-55%. FAMCAT (Weng 2015) external AUC 0.77 but no systematic treatment adjustment.')

H('1.3 Cascade screening paradox', 3)
P('A 28-year-old daughter carrying the identical LDLR variant as her proband scores 3 on DLCN (below threshold 6) while her father scores 12 -- the diagnosis is missed despite identical genetic defect. TUDOR addresses this through the Index Effect feature.')

H('1.4 Rationale and study aims', 3)
P_mixed([
    ('Trig_Filter is the manuscript shorthand for a single composite feature combining untreated LDL-C and triglycerides; clinically it captures the LDL-and-TG fingerprint that distinguishes FH from polygenic hypercholesterolaemia. ', False, True, INS),
    ('[R2.1 preface -- inserted]', True, False, TAG),
])
P('Three innovations: (i) individualised treatment-intensity adjusted LDL-C; (ii) Index Effect modelling differential diagnostic significance by ascertainment pathway; (iii) Triglyceride Metabolic Filter.')

# Section 2
H('2. Methods', 1)
P('TRIPOD Type 2b for Wales (model from one set, validated in same population); Type 4 for UK Biobank (entirely separate population -- most stringent external validation).', bold=True)
P('Wales training: CAVUHB specialist clinic, index cases only. External validation: All-Wales PASS Registry (n = 7,253; 2,405 genetic FH+; index 63.6% / cascade 36.4%).')
P('UK Biobank: lipid clinic-mimicking subset (TC > 7.5; statin-corrected LDL > 4.9; or premature ASCVD before 55M/60F) -> n = 58,021 with 729 FH (1.26%). Application 1002450.')
P('Treatment-adjusted LDL-C: drug-specific reductions (atorvastatin 25-48%, rosuvastatin 35-55%, simvastatin 20-42%, pravastatin 15-29%, fluvastatin 15-22%), adherence factor (0.5/0.75/1.0), non-statin add-ons (ezetimibe 20%, bempedoic 25%, PCSK9i 65%), cap 85%.')
P('Novel features: Index_Effect = (1 - Is_Relative) x LDL_Untreated; Trig_Filter = LDL_Untreated / (Triglycerides + 0.1); personalised statin calibration residual (LOO-CV on 298 paired pre/post LDL).')
P('Elastic Net alpha = 0.5. 10 retained features: Index_Effect (+1.148), LDL_adjusted (+1.070), Trig_Filter (+0.871), statin residual (-0.677), age (-0.497), sex (-0.208), HDL (-0.199), tendon xanthomata (+0.159), TG (-0.141), corneal arcus (+0.120). 5-fold CV AUC 0.867 (SD 0.025).')
P_mixed([
    ('In plain terms, the Elastic Net selects only the most informative clinical features and shrinks the rest to zero -- producing a parsimonious diagnostic score easier to inspect, audit, and deploy than a black-box ML model. ', False, True, INS),
    ('[R2.1 preface -- inserted in 2.5]', True, False, TAG),
])
P_mixed([
    ('In clinical terms, the NRI quantifies how many additional patients are placed in the correct probability band -- and thus the correct clinical action -- when TUDOR replaces the older diagnostic criterion. The calibration slope describes how well predicted probabilities track observed FH frequencies; slope > 1 indicates the model needs scaling-down (recalibration) for lower-prevalence populations. ', False, True, INS),
    ('[R2.1 preface -- inserted in 2.10]', True, False, TAG),
])
P('Ethical approval: CAVUHB R&D; UKB app 1002450; NW MREC 11/NW/0382.')

# Section 3
H('3. Results', 1)

H('3.1 Study populations', 3)
P('Wales: mean age 45.9, 57.2% F, BMI 28.7. FH+ younger (39.3 vs 49.2 y; d = -0.571). Trig_Filter d = 0.744 (largest of any continuous variable). Cascade 57.0% in FH+ vs 26.1% in FH-.')
P('UK Biobank: mean age 58.9, 55.1% F, BMI 28.4. FH+ lower BMI (27.8 vs 28.5). Trig_Filter d = -0.963. ApoB/LDL ratio d = 0.384.')

H('3.2 Treatment adjustment validation', 3)
P('Wales 127 paired pre/post LDL: r = 0.89, MAE 0.48 mmol/L, B-A LoA -0.92 to +0.98. UKB longitudinal (n = 1,254 statin initiators): simvastatin 35.1%/35.0%, atorvastatin 36.7%/38.0%, rosuvastatin 37.4%/34.0%. ICC across 4 LDL equations 0.986.')

H('3.3 Wales external validation (TRIPOD 2b)', 3)
P('TUDOR AUC 0.842 (0.822-0.863); DLCN 0.791 (p < 0.001); NRI 0.358; calibration slope 1.039. Index DLCN sensitivity 88.5%; cascade DLCN sensitivity 1.5% -- 87-pp gap. TUDOR cascade sensitivity 89.4% versus FAMCAT-approx 66.7%.', bold=True)

H('3.4 UK Biobank external validation (TRIPOD 4)', 3)
P('TUDOR AUC 0.750 (0.731-0.770) vs eDLCN 0.636 (DeLong Z = 10.08, p = 6.73 x 10^-24), LDL alone 0.689, Trig_Filter alone 0.700. Youden Sens 59.7%, Spec 79.3%, Brier 0.069, calibration slope 6.33 (5.93-6.73) -- recalibration formula provided. 12-fold detection improvement (639 vs 52 per 10,000).', bold=True)

H('3.5 Genetic validation', 3)
P('Wales: LDLR 1,901 (79.0%), APOB 281 (11.7%), PCSK9 23 (1.0%), LDLR CNV 164 (6.8%), APOE 33 (1.4%), LDLRAP1 3 (0.1%). UKB: LDLR 1,321 (81.4%), APOB 301 (18.5%), PCSK9 1 (0.1%). Total 4,028 genetic FH.')
P('Gene-specific AUCs (Wales/UKB): LDLR 0.839/0.717, APOB 0.841/0.830, APOE 0.809/-. APOB vs LDLR in UKB DeLong Z = 3.74, p = 1.88 x 10^-4. Cohen d for Trig_Filter 1.108 (APOB) vs 0.672 (LDLR).')

H('3.6 ApoB augmentation', 3)
P('Base TUDOR 0.753. Model A (+ApoB) 0.751 (p = 0.349, ns). Model B (+ratio) 0.760 (p = 0.006). Model C (+both) 0.771 (delta-AUC +0.019, p = 3.99 x 10^-5). Model D (+Lp(a)) 0.776.')
P('Grey-zone analysis (intermediate-probability stratum n = 33,562): ratio AUC 0.513 -- operates at whole-cohort ranking, not within-stratum triage.')

H('3.7 Trig_Filter biological validation', 3)
P('Sequential adjustment ORs 1.873 -> 1.891 (BMI) -> 1.903 (DM) -> 1.920 (ApoB) -> 1.935 (age+sex). VIF < 1.2. Metabolic shield Cohen d 0.89 (TG<1.7) -> 0.52 (1.7-2.3) -> 0.18 (>2.3). NMR (421,122 UKB; 1,612 genetic FH): 35/53 metabolites differ post-Bonferroni; TUDOR alone 0.749, NMR alone 0.705, combined 0.752 (p = 0.278) -- clinical model captures the molecular biology.')

H('3.8 CV mortality gradient', 3)
P('n = 58,021, 1,893 deaths, median FU 14.2 y. IHD 0.65% -> 1.76% across LDL deciles (OR 1.036/decile, p = 0.011). Stroke inverse gradient 0.89% -> 0.35% (haemorrhagic-stroke J-curve). All-CV OR 1.012 (p = 0.14).')

H('3.9 T2DM ascertainment', 3)
P('Self-report 1.9%, ICD-10 5.6%, HbA1c 1.5%. kappa = 0.28-0.36 (poor concordance). All three confirm attenuation: beta = -0.484 / -0.473 / -0.658. HbA1c dose-response: Trig_Filter 2.74 -> 2.31 -> 1.82 across normo/IFG/diabetes.')

H('3.10 Sensitivity analyses and algorithm comparison', 3)
P('10 prespecified sensitivity analyses pass. ENET (AUC 0.781) vs XGBoost (0.778) DeLong p = 0.72 -- clinical feature engineering, not algorithmic complexity, drives diagnostic improvement.')

# Section 4
H('4. Discussion', 1)
P_mixed([
    ('Synthesis opener (R2.5 insertion): ', True, False, TAG),
    ('TUDOR addresses a practical gap left by criteria designed in the pre-statin era. DLCN, MEDPED, and Simon Broome were validated against treatment-naive LDL-C distributions; in contemporary populations the same patients now present on therapy with LDL-C attenuated below the original thresholds. By incorporating on-statin status and back-calculated untreated LDL-C alongside an interpretable Elastic Net feature set, TUDOR preserves diagnostic performance in heavily treated cohorts without requiring a return-to-baseline washout.', False, True, INS),
])

H('4.1 Principal findings', 3)
P('Dual external validation, AUC 0.842 (Wales) / 0.750 (UKB), 4,028 genetic FH across six genes -- the largest dual-validated genetically confirmed FH diagnostic study to date.')

H('4.2 Cascade screening imperative', 3)
P_mixed([
    ('External context: Weng et al [18] developed FAMCAT in UK primary care using prospectively coded features; their cohort, definitions, and lipid distributions differ materially from ours, so head-to-head comparison is illustrative rather than equivalent. Our internally reconstructed FAMCAT-like score is therefore best interpreted as a methodologically necessary comparator limited by the absence of primary-care longitudinal coding depth, and should not be read as a direct benchmark. ', False, True, INS),
    ('[R2.3 FAMCAT contextualisation -- inserted]', True, False, TAG),
])
P('The 87-pp sensitivity gap (88.5% index -> 1.5% cascade) using DLCN represents the central diagnostic deficiency. Simon Broome depends on tendon xanthomata (<15% het FH) and corneal arcus (2.0% here) -- increasingly rare in the statin era.')

H('4.3 Genetic heterogeneity and variant-specific performance', 3)
P('APOB p.Arg3527Gln produces a uniform molecular defect (impaired LDL-receptor binding) -> reproducible phenotype -> higher TUDOR discrimination. >2,000 LDLR pathogenic variants span null, defective, CNV -> wider phenotypic distribution -> lower discriminative efficiency.')

H('4.4 Cardiovascular outcomes validation', 3)
P('IHD trend confirms continuous LDL-mortality relationship (Brown-Goldstein causal framework). Inverse stroke gradient reflects mixed thrombotic/haemorrhagic aetiology -- not protective LDL. Guidelines target <1.4 mmol/L in very-high-risk, well above haemorrhagic thresholds.')

H('4.5 T2DM ascertainment and Trig_Filter robustness', 3)
P('Low kappa (0.28-0.36) strengthens the conclusion: if the effect were a coding artefact, three discordant definitions would not converge. HbA1c dose-response is the mechanism (SREBP-1c VLDL overproduction, LPL impairment).')

H('4.6 ApoB augmentation and precision lipidology', 3)
P('Absolute ApoB redundant given treatment-adjusted LDL; the ApoB/LDL-C ratio captures cholesterol-content-per-particle (lower in FH due to prolonged residence). Ratio < 0.31 predicts ASCVD in genetic FH (27.6 vs 11.8%, p = 0.0022).')

H('4.7 Clinical Implementation Framework', 3)
P_mixed([
    ('This pathway is intentionally simple: one probability band drives one clinical action and one genetic-testing decision. It supplements rather than replaces existing criteria (DLCN, MEDPED, Simon Broome), and assumes Bayesian recalibration to local FH prevalence when deployed outside specialist lipid-clinic populations. ', False, True, INS),
    ('[R2.1 preface -- inserted]', True, False, TAG),
])
P('Low-probability (TUDOR < 25%): reassurance, standard lipid pathway.')
P('Intermediate (25-75%): detailed review, secondary-cause work-up, intensify LLT, re-evaluate after adherence-stable repeat.')
P('High (>75%): direct referral for genetic confirmation (panel/WES) and cascade.')
P_mixed([
    ('To translate TUDOR into routine clinical use we propose a three-tier decision pathway based on the model predicted probability of monogenic FH. The thresholds (<25%, 25-75%, >75%) reuse the categorical bands already used in the NRI analysis and map to distinct clinical actions summarised in Table 5. ', False, True, INS),
    ('[R2.2 expansion -- inserted]', True, False, TAG),
])
P('Note: the original v4 three-tier paragraph is retained immediately below per the no-deletion rule. Items #1-2 of the flag list at the end propose a 30-second manual trim.', italic=True, colour=WARN)

H('4.8 Limitations', 3)
P_mixed([
    ('Adherence to lipid-lowering therapy in real-world settings is heterogeneous and not fully captured by binary prescription records. Intermittent dosing, dose reduction during illness, and complete cessation can leave a prescription record that overstates true exposure, distorting back-calculation of untreated LDL-C and attenuating Trig_Filter. External validation in All-Wales-PASS and UK Biobank gives some reassurance, but a definitive test requires time-varying GP-prescription data linked to dispensing records. We acknowledge this as a structural limitation and are pursuing linkage as part of the Cardiff longitudinal cohort follow-up. ', False, True, INS),
    ('[R2.4 adherence limitation -- inserted]', True, False, TAG),
])
P('Other limitations: single-centre development; UKB age 40-69; FAMCAT comparison internal approximation; calibration slope 6.33 needs Bayesian recalibration; ApoB augmentation testable only in UKB; T2DM attenuates but does not preclude; ICD-10 p41280 ASCVD event dates not yet extracted; prospective implementation studies needed.')

# Section 5
H('5. Conclusions', 1)
P('TUDOR is a dual-externally-validated, interpretable FH diagnostic algorithm maintaining discrimination across specialist (AUC 0.842) and population (AUC 0.750) settings in 4,028 genetic FH across six genes. APOB AUC 0.830 reflects molecular homogeneity. Optional ApoB/LDL-C ratio augmentation lifts AUC to 0.771. IHD mortality gradient independently confirms LDL-C as continuous CV risk determinant. Three-definition T2DM concordance confirms the Trig_Filter attenuation is biological not artefactual. Treatment-adjusted phenotyping and ascertainment correction address the structural drivers of FH underdiagnosis.')

# Key Tables
H('Key Tables', 1)

H('Table 1a. Wales baseline (selected rows; n = 7,253)', 3)
t = doc.add_table(rows=1, cols=5)
t.style = 'Light List Accent 2'
hdr = t.rows[0].cells
for i, txt in enumerate(['Variable', 'FH+ (n=2,405)', 'FH- (n=4,848)', 'P', 'Cohen d']):
    hdr[i].text = txt
for r in [
    ('Age, y', '39.3 (17.1)', '49.2 (17.2)', '<10^-99', '-0.571'),
    ('LDL-C measured, mmol/L', '5.8 (2.0)', '5.1 (1.8)', '<10^-31', '0.349'),
    ('LDL-C adjusted, mmol/L', '6.5 (2.1)', '5.7 (1.9)', '<10^-33', '0.378'),
    ('Triglycerides, mmol/L', '1.4 (0.8)', '2.0 (0.9)', '<10^-72', '-0.460'),
    ('Trig_Filter', '5.1 (1.9)', '3.3 (1.5)', '<10^-90', '0.744'),
    ('Tendon xanthomata', '10.8%', '5.7%', '<10^-17', '-'),
    ('Cascade relative', '57.0%', '11.2%', '<10^-145', '-'),
]:
    row = t.add_row().cells
    for i, v in enumerate(r):
        row[i].text = v

H('Table 1b. UK Biobank baseline (selected rows; n = 58,021)', 3)
t = doc.add_table(rows=1, cols=5)
t.style = 'Light List Accent 2'
hdr = t.rows[0].cells
for i, txt in enumerate(['Variable', 'FH+ (n=729)', 'FH- (n=57,292)', 'P', 'Cohen d']):
    hdr[i].text = txt
for r in [
    ('LDL-C measured', '5.7 (0.9)', '5.1 (0.9)', '<10^-89', '-0.739'),
    ('Triglycerides', '1.8 (0.9)', '2.2 (1.2)', '<10^-11', '-0.411'),
    ('Trig_Filter', '3.9 (1.5)', '2.7 (1.4)', '<10^-54', '-0.963'),
    ('ApoB/LDL-C ratio', '0.21 (0.05)', '0.23 (0.06)', '<10^-22', '0.384'),
    ('Statin use', '53.5%', '38.4%', '<10^-17', '-'),
]:
    row = t.add_row().cells
    for i, v in enumerate(r):
        row[i].text = v

H('Table 2. Discrimination performance (both cohorts)', 3)
t = doc.add_table(rows=1, cols=6)
t.style = 'Light List Accent 2'
hdr = t.rows[0].cells
for i, txt in enumerate(['Model', 'Cohort', 'AUC', '95% CI', 'DeLong vs TUDOR', 'Sens / Spec %']):
    hdr[i].text = txt
for r in [
    ('TUDOR v2', 'Wales', '0.842', '0.822-0.863', '-', '89.4 / 67.3*'),
    ('DLCN', 'Wales', '0.791', '0.770-0.813', 'p<0.001', '88.5 / 72.1+'),
    ('FAMCAT-approx', 'Wales', '0.766', '0.744-0.788', 'p=0.006', '66.7 / 69.2*'),
    ('TUDOR v2', 'UK Biobank', '0.750', '0.731-0.770', '-', '59.7 / 79.3'),
    ('eDLCN', 'UK Biobank', '0.636', '0.609-0.663', 'p=6.73x10^-24', '5.2 / 96.1'),
    ('TUDOR + ApoB ratio', 'UK Biobank', '0.771', '0.753-0.789', 'p=3.99x10^-5', '64.1 / 78.8'),
]:
    row = t.add_row().cells
    for i, v in enumerate(r):
        row[i].text = v
P('* cascade relatives only ; + index cases', italic=True)

H('Table 3. CV mortality by LDL decile (UKB; n = 58,021)', 3)
t = doc.add_table(rows=1, cols=5)
t.style = 'Light List Accent 2'
hdr = t.rows[0].cells
for i, txt in enumerate(['Decile', 'Mean LDL', 'n', 'IHD deaths', 'Stroke deaths']):
    hdr[i].text = txt
for d, ldl, ihd, strk in [(1, 3.2, '38 (0.65%)', '52 (0.90%)'),
                          (5, 5.0, '69 (1.19%)', '41 (0.71%)'),
                          (10, 7.8, '102 (1.76%)', '20 (0.34%)')]:
    row = t.add_row().cells
    row[0].text = str(d)
    row[1].text = str(ldl)
    row[2].text = '5,802'
    row[3].text = ihd
    row[4].text = strk
P('IHD logistic OR/decile 1.036 (1.008-1.065), p = 0.011. Stroke OR 0.931 (p < 0.001). All-CV OR 1.012 (p = 0.14).', italic=True)

H('Table 4. T2DM ascertainment -- three-definition concordance and Trig_Filter attenuation', 3)
t = doc.add_table(rows=1, cols=5)
t.style = 'Light List Accent 2'
hdr = t.rows[0].cells
for i, txt in enumerate(['Definition', 'Prevalence', 'beta', '95% CI', 'P']):
    hdr[i].text = txt
for r in [
    ('A: Self-report (p2443)', '1.9% (n=1,102)', '-0.484', '-0.611 to -0.357', '8.4x10^-13'),
    ('B: ICD-10 E11', '5.6% (n=3,249)', '-0.473', '-0.559 to -0.387', '8.2x10^-27'),
    ('C: HbA1c >= 48 mmol/mol', '1.5% (n=870)', '-0.658', '-0.793 to -0.523', '1.1x10^-29'),
]:
    row = t.add_row().cells
    for i, v in enumerate(r):
        row[i].text = v
P('All adjusted age/sex/BMI/alcohol. kappa pairwise: A-B 0.281; A-C 0.358; B-C 0.299. HbA1c dose-response Trig_Filter: 2.74 (normo) -> 2.31 (IFG) -> 1.82 (diabetes).', italic=True)

H('Table 5. Three-tier clinical decision pathway (R2.2 addition -- proposed rename from Table 4)', 3)
t = doc.add_table(rows=1, cols=3)
t.style = 'Light List Accent 2'
hdr = t.rows[0].cells
for i, txt in enumerate(['TUDOR band', 'Clinical action', 'Genetic-testing indication']):
    hdr[i].text = txt
for r in [
    ('Low (<25%)', 'Reassurance; standard lipid management', 'Not routinely indicated'),
    ('Intermediate (25-75%)', 'Detailed review; secondary-cause work-up; intensify LLT', 'Discuss; consider after adherence-stable re-evaluation'),
    ('High (>75%)', 'Specialist lipid-clinic referral; treat as probable FH', 'Targeted-panel or WES; initiate cascade screening of FDR'),
]:
    row = t.add_row().cells
    for i, v in enumerate(r):
        row[i].text = v

# Items to fix
H('Items I would flag before you upload', 1)
items = [
    ('Table 4 numbering collision.', 'Two tables share the label "Table 4" -- the T2DM concordance table (original v4) and the new three-tier decision pathway table (R2.2 addition). Rename the new one to Table 5 in both caption and in the section 4.7 prose reference.'),
    ('Section 4.7 reads thrice.', 'The original v4 three-tier paragraph was retained per no-deletion rule, plus two new paragraphs inserted ahead. All three describe the same framework. Strictly rule-compliant, but reads redundantly. Suggest trimming the original v4 paragraph to a single sentence linking to Table 5.'),
    ('Table 4 (decision pathway) is appended at the end of the document.', 'Drag it into Section 4.7 in Word (10 seconds).'),
    ('Discussion synthesis opener sits between heading 4 and 4.1.', 'Intended placement per R2.5 -- confirm it reads as you expect.'),
    ('References 51-55 in bibliography.', 'Final cross-check that all are cited inline in the body before submission.'),
]
for ttl, body in items:
    p = doc.add_paragraph(style='List Number')
    r = p.add_run(ttl + ' ')
    r.bold = True
    r.font.name = 'Cambria'
    r.font.size = Pt(11)
    r2 = p.add_run(body)
    r2.font.name = 'Cambria'
    r2.font.size = Pt(11)
P('None of these touch any numerical claim. All are editorial fixes resolvable in Word in under five minutes.', italic=True, bold=True)

# Legend
H('Legend used in this review document', 2)
P_mixed([
    ('Italic green text = ', False, True, None),
    ('paragraph inserted in v5 (was not in v4).', False, True, INS),
])
P_mixed([
    ('Blue bold tag [R2.x] = ', False, True, None),
    ('which Reviewer #2 item this insertion addresses.', False, True, TAG),
])

out = 'TUDOR_v5_REVIEW.docx'
doc.save(out)
import os
print('wrote ' + out + ' (' + str(round(os.path.getsize(out) / 1024, 1)) + ' KB)')

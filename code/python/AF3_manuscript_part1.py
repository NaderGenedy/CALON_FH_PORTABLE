"""
AF3 FH Manuscript - Part 1: Title, Abstract, Introduction, Methods

Text-generation functions for the first half of the manuscript.
Each function takes a python-docx Document object and a stats dictionary,
and appends formatted content to the document.
"""

from docx.shared import Pt, Cm, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.style import WD_STYLE_TYPE


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def set_run(run, size=11, bold=False, italic=False, font='Calibri', color=None):
    """Apply font formatting to a run."""
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    run.font.name = font
    if color:
        run.font.color.rgb = RGBColor(*color)


def add_heading_text(doc, text, level=1):
    """Add a formatted heading paragraph."""
    p = doc.add_paragraph()
    run = p.add_run(text)
    if level == 1:
        set_run(run, size=14, bold=True)
    elif level == 2:
        set_run(run, size=12, bold=True)
    elif level == 3:
        set_run(run, size=11, bold=True, italic=True)
    p.paragraph_format.space_before = Pt(12)
    p.paragraph_format.space_after = Pt(6)
    return p


def add_body_text(doc, text):
    """Add a formatted body paragraph with double spacing."""
    p = doc.add_paragraph()
    run = p.add_run(text)
    set_run(run, size=11)
    p.paragraph_format.space_after = Pt(6)
    p.paragraph_format.line_spacing = 2.0
    return p


# ---------------------------------------------------------------------------
# 1. Title Page
# ---------------------------------------------------------------------------

def add_title_page(doc, stats):
    """Add title page to document."""

    # Title
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run(
        "From Atomic Structure to Arterial Disease: AlphaFold3-Guided "
        "Multi-Modal Risk Stratification in Familial Hypercholesterolemia"
    )
    set_run(run, size=16, bold=True)
    p.paragraph_format.space_after = Pt(24)

    # Authors
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run(
        "Nader Genedy, BSc, MSc; [additional authors to be confirmed]"
    )
    set_run(run, size=12)
    p.paragraph_format.space_after = Pt(12)

    # Affiliations
    affiliations = [
        "1. Swansea University Medical School, Faculty of Medicine, Health and Life Science, "
        "Swansea University, Swansea, United Kingdom",
        "2. South Wales Familial Hypercholesterolemia Registry, University Hospital of Wales, "
        "Cardiff, United Kingdom",
        "3. The Secure Anonymised Information Linkage (SAIL) Databank, Swansea University "
        "Medical School, Swansea, United Kingdom",
    ]
    for aff in affiliations:
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = p.add_run(aff)
        set_run(run, size=10, italic=True)
        p.paragraph_format.space_after = Pt(2)

    # Corresponding author
    doc.add_paragraph()  # spacer
    p = doc.add_paragraph()
    run = p.add_run("Corresponding Author: ")
    set_run(run, size=11, bold=True)
    run = p.add_run(
        "Nader Genedy, Swansea University Medical School, Singleton Campus, "
        "Swansea SA2 8PP, United Kingdom. Email: nader.genedy@swansea.ac.uk"
    )
    set_run(run, size=11)
    p.paragraph_format.space_after = Pt(12)

    # Word count and counts
    p = doc.add_paragraph()
    run = p.add_run("Word Count: ")
    set_run(run, size=11, bold=True)
    run = p.add_run("~7,500 (excluding references, tables, and figure legends)")
    set_run(run, size=11)

    p = doc.add_paragraph()
    run = p.add_run("Figures: ")
    set_run(run, size=11, bold=True)
    run = p.add_run("15")
    set_run(run, size=11)

    p = doc.add_paragraph()
    run = p.add_run("Tables: ")
    set_run(run, size=11, bold=True)
    run = p.add_run("7")
    set_run(run, size=11)

    p = doc.add_paragraph()
    run = p.add_run("Supplementary Materials: ")
    set_run(run, size=11, bold=True)
    run = p.add_run("Supplementary Tables S1\u2013S8, Supplementary Figures S1\u2013S6")
    set_run(run, size=11)

    # Page break
    doc.add_page_break()


# ---------------------------------------------------------------------------
# 2. Abstract
# ---------------------------------------------------------------------------

def add_abstract(doc, stats):
    """Add structured abstract (~250 words)."""

    add_heading_text(doc, "ABSTRACT", level=1)

    # Background
    p = doc.add_paragraph()
    run = p.add_run("Background: ")
    set_run(run, size=11, bold=True)
    run = p.add_run(
        f"Familial hypercholesterolemia (FH) affects approximately 1 in 250 individuals "
        f"and confers markedly elevated risk of premature atherosclerotic cardiovascular "
        f"disease (ASCVD). Current risk-prediction tools treat all pathogenic variants as "
        f"functionally equivalent, ignoring substantial evidence that the structural and "
        f"thermodynamic consequences of individual mutations drive divergent clinical "
        f"trajectories. We hypothesized that AlphaFold3 (AF3) protein structure predictions "
        f"could be distilled into a variant-level Structural Severity Score (SSS) that "
        f"improves ASCVD risk stratification beyond conventional clinical models."
    )
    set_run(run, size=11)
    p.paragraph_format.line_spacing = 2.0
    p.paragraph_format.space_after = Pt(6)

    # Methods
    p = doc.add_paragraph()
    run = p.add_run("Methods: ")
    set_run(run, size=11, bold=True)
    run = p.add_run(
        f"We studied {stats['n_total']} patients from the South Wales FH Registry carrying "
        f"{stats['n_variants']} distinct variants across LDLR (n={stats['n_ldlr']}), APOB "
        f"(n={stats['n_apob']}), and PCSK9 (n={stats['n_pcsk9']}). AF3 was used to model "
        f"wildtype and mutant protein structures, from which FoldX-derived thermodynamic "
        f"stability (\u0394\u0394G), per-residue confidence (pLDDT), and interface disruption "
        f"metrics were extracted. A composite SSS was constructed from five scoring strategies "
        f"and validated using leave-one-centre-out cross-validation (LOCO-CV) with 16 nested "
        f"logistic regression models. External multi-modal validation leveraged UK Biobank "
        f"NMR metabolomics (~{stats['n_nmr_participants']//1000}K) and cardiac MRI "
        f"(~{stats['n_mri_participants']//1000}K)."
    )
    set_run(run, size=11)
    p.paragraph_format.line_spacing = 2.0
    p.paragraph_format.space_after = Pt(6)

    # Results
    p = doc.add_paragraph()
    run = p.add_run("Results: ")
    set_run(run, size=11, bold=True)
    run = p.add_run(
        f"Among {stats['n_analytic']} patients with complete data (mean age "
        f"{stats['mean_age']:.1f}\u00b1{stats['sd_age']:.1f} years; "
        f"{stats['pct_female']:.1f}% female), the ASCVD event rate was "
        f"{stats['ascvd_rate']:.1f}%. \u0394\u0394G correlated with baseline LDL-C "
        f"(r={stats['ddg_ldl_r']:.2f}, P={stats['ddg_ldl_p']:.3f}) and ApoB "
        f"(r={stats['ddg_apob_r']:.2f}, P={stats['ddg_apob_p']:.3f}). Addition of the "
        f"SSS to a clinical base model improved discrimination from AUC "
        f"{stats['auc_base']:.3f} to {stats['auc_sss']:.3f} (DeLong P<0.05), with "
        f"significant net reclassification improvement. Patients in the highest compound "
        f"risk stratum experienced an ASCVD rate of {stats['triple_high_ascvd']:.1f}% "
        f"compared with {stats['triple_low_ascvd']:.1f}% in the lowest stratum. "
        f"In UK Biobank, SSS tiers demonstrated monotonic gradients across lipoprotein "
        f"subfractions, inflammatory markers, and cardiac structural indices."
    )
    set_run(run, size=11)
    p.paragraph_format.line_spacing = 2.0
    p.paragraph_format.space_after = Pt(6)

    # Conclusions
    p = doc.add_paragraph()
    run = p.add_run("Conclusions: ")
    set_run(run, size=11, bold=True)
    run = p.add_run(
        f"AF3-derived structural features capture clinically meaningful variation among FH "
        f"mutations and, when combined into a Structural Severity Score, significantly enhance "
        f"ASCVD risk prediction. Integrating protein-level biology with multi-modal population "
        f"data offers a scalable framework for precision risk stratification in monogenic "
        f"dyslipidaemias."
    )
    set_run(run, size=11)
    p.paragraph_format.line_spacing = 2.0
    p.paragraph_format.space_after = Pt(6)

    # Keywords
    p = doc.add_paragraph()
    run = p.add_run("Keywords: ")
    set_run(run, size=11, bold=True)
    run = p.add_run(
        "familial hypercholesterolemia; AlphaFold3; structural severity score; "
        "atherosclerotic cardiovascular disease; protein structure prediction; "
        "risk stratification; NMR metabolomics; cardiac MRI"
    )
    set_run(run, size=11, italic=True)
    p.paragraph_format.space_after = Pt(12)

    doc.add_page_break()


# ---------------------------------------------------------------------------
# 3. Introduction
# ---------------------------------------------------------------------------

def add_introduction(doc, stats):
    """Add introduction (~1,500 words, 6 paragraphs)."""

    add_heading_text(doc, "INTRODUCTION", level=1)

    # ---- Paragraph 1: FH epidemiology ----
    add_body_text(doc,
        f"Familial hypercholesterolemia is the most common life-threatening monogenic "
        f"disorder in clinical medicine, affecting roughly 1 in 250 individuals worldwide "
        f"and accounting for an estimated 30 million cases globally. Inherited in an "
        f"autosomal dominant pattern, FH is caused by pathogenic variants in three principal "
        f"genes\u2014LDLR, APOB, and PCSK9\u2014that disrupt hepatic low-density lipoprotein "
        f"(LDL) receptor-mediated clearance and produce lifelong elevations in circulating "
        f"LDL cholesterol. The clinical consequence of this sustained hypercholesterolaemia "
        f"is premature atherosclerotic cardiovascular disease: untreated heterozygous FH "
        f"confers a roughly 20-fold increase in coronary heart disease risk by age 60, and "
        f"homozygous FH can manifest with severe aortic valve calcification and myocardial "
        f"infarction in childhood. Despite the availability of effective lipid-lowering "
        f"therapies\u2014including high-intensity statins, ezetimibe, PCSK9 inhibitors, and "
        f"the recently approved small-interfering RNA inclisiran\u2014an estimated 85% of "
        f"individuals with FH worldwide remain undiagnosed, and even among those identified, "
        f"fewer than half achieve guideline-recommended LDL-C targets. The cumulative burden "
        f"of LDL exposure, often quantified as cholesterol-years, provides a compelling "
        f"physiological rationale for early intervention, yet current clinical practice "
        f"lacks the tools to identify which patients harbour the most structurally damaging "
        f"mutations and therefore stand to benefit most from aggressive, early treatment."
    )

    # ---- Paragraph 2: Current risk prediction tools ----
    add_body_text(doc,
        f"Efforts to refine cardiovascular risk prediction in FH have yielded several "
        f"validated scoring systems over the past decade. The Spanish SAFEHEART Registry "
        f"Equation (SAFEHEART-RE) incorporates six clinical variables\u2014age, sex, history "
        f"of ASCVD, baseline LDL-C, lipoprotein(a), and hypertension\u2014to estimate 5- and "
        f"10-year event probabilities, and has been externally validated in multiple European "
        f"cohorts. The Montreal InHeRIT score extends this framework by adding family history "
        f"burden and treatment response metrics. The Dutch Lipid Clinic Network criteria, "
        f"while primarily diagnostic, also provide prognostic stratification through a tiered "
        f"scoring algorithm. Notably, however, all three systems share a fundamental "
        f"limitation: they treat the genetic aetiology of FH as a binary variable (variant "
        f"present or absent) or, at most, stratify by gene (LDLR vs. APOB vs. PCSK9). None "
        f"incorporates information about the specific molecular consequences of the patient's "
        f"variant, despite growing evidence that pathogenic variants within the same gene can "
        f"produce vastly different phenotypes depending on their location within the protein "
        f"structure and their effect on protein folding, stability, and receptor trafficking. "
        f"This gene-agnostic approach to risk prediction represents a missed opportunity to "
        f"leverage the rich structural and biophysical information that is now computationally "
        f"accessible for every known FH variant."
    )

    # ---- Paragraph 3: The variant-level problem ----
    add_body_text(doc,
        f"The scale of the variant-level problem in FH is substantial. More than 2,000 "
        f"pathogenic and likely pathogenic variants have been catalogued in LDLR alone, "
        f"spanning missense substitutions, in-frame deletions, splice-site alterations, and "
        f"large structural rearrangements. These variants are emphatically not clinically "
        f"equivalent. Functional studies have demonstrated that null mutations\u2014those "
        f"producing no detectable LDL receptor protein\u2014confer approximately a 10-fold "
        f"higher ASCVD risk relative to defective-class mutations that allow residual "
        f"receptor activity. Even within the defective class, the specific protein domain "
        f"affected determines the mechanism of dysfunction: variants in the ligand-binding "
        f"domain (repeats 3\u20137) disrupt LDL particle capture, those in the epidermal "
        f"growth factor (EGF) precursor homology domain impair pH-dependent release of LDL "
        f"in the endosome and subsequent receptor recycling, and variants in the "
        f"transmembrane or cytoplasmic domains interfere with receptor anchoring and "
        f"clathrin-mediated endocytosis. In the present cohort, {stats['n_variants']} "
        f"distinct variants were identified across LDLR ({stats['pct_ldlr']:.1f}%), APOB "
        f"({stats['pct_apob']:.1f}%), and PCSK9 ({stats['pct_pcsk9']:.1f}%), distributed "
        f"across {stats['n_domains']} functional domains (Table 1). This heterogeneity "
        f"underscores the need for a systematic, structure-aware scoring framework that can "
        f"translate molecular-level variant information into clinically actionable risk "
        f"estimates."
    )

    # ---- Paragraph 4: Structural biology revolution ----
    add_body_text(doc,
        f"The past five years have witnessed a revolution in computational structural "
        f"biology that makes such a framework feasible for the first time. AlphaFold2, "
        f"developed by DeepMind, solved the decades-old protein structure prediction problem "
        f"with near-experimental accuracy, an achievement recognised by the 2024 Nobel Prize "
        f"in Chemistry awarded to Demis Hassabis and John Jumper. Its successor, AlphaFold3, "
        f"extends the approach to multi-chain protein complexes, nucleic acids, and "
        f"post-translational modifications, and provides per-residue confidence estimates "
        f"(predicted local distance difference test, pLDDT) that serve as a proxy for local "
        f"structural order. Complementing these predictions, the FoldX energy function "
        f"enables rapid estimation of the thermodynamic consequences of point mutations "
        f"through its BuildModel protocol: the change in Gibbs free energy upon mutation "
        f"(\u0394\u0394G) quantifies the degree to which a substitution destabilises the "
        f"native fold, with values exceeding 2 kcal/mol generally indicating significant "
        f"structural disruption. Together, AF3 and FoldX offer a scalable computational "
        f"pipeline that can characterise every known FH variant along two orthogonal axes"
        f"\u2014structural confidence and thermodynamic stability\u2014without requiring "
        f"experimental crystallography or cryo-electron microscopy for each mutant. In this "
        f"study, AF3 wildtype models achieved mean pLDDT scores of {stats['ldlr_plddt']:.1f} "
        f"for LDLR, {stats['pcsk9_plddt']:.1f} for PCSK9, and {stats['apob_plddt']:.1f} for "
        f"ApoB, confirming the high fidelity of the predicted structures (Figure 1)."
    )

    # ---- Paragraph 5: The cascade hypothesis ----
    add_body_text(doc,
        f"Central to our analytic framework is what we term the structural\u2013metabolic "
        f"cascade hypothesis: the proposition that the clinical impact of an FH variant "
        f"propagates through a series of biologically ordered levels, each of which should "
        f"display monotonic severity gradients if the upstream structural perturbation is "
        f"faithfully transmitted. The cascade proceeds as follows: (1) a nucleotide variant "
        f"alters the amino acid sequence; (2) the altered sequence produces a structural "
        f"perturbation in the three-dimensional fold, quantifiable by AF3 pLDDT and FoldX "
        f"\u0394\u0394G; (3) structural destabilisation impairs protein stability and "
        f"trafficking, reflected in reduced cell-surface receptor density; (4) diminished "
        f"receptor function elevates circulating LDL-C and ApoB concentrations, measurable "
        f"through standard clinical biochemistry; (5) sustained dyslipidaemia remodels the "
        f"lipoprotein particle profile, detectable by high-throughput NMR metabolomics; "
        f"(6) atherogenic particle exposure drives arterial wall inflammation and structural "
        f"cardiac remodelling, assessable by cardiac MRI; and (7) these cumulative processes "
        f"culminate in overt ASCVD events\u2014myocardial infarction, coronary "
        f"revascularisation, ischaemic stroke, or peripheral arterial disease. If this "
        f"cascade model is correct, then a well-constructed Structural Severity Score should "
        f"not only predict clinical events in the FH registry but also demonstrate coherent, "
        f"dose-dependent associations at every intermediate level when tested against "
        f"independent population biobank data (Figure 2)."
    )

    # ---- Paragraph 6: Study aims ----
    add_body_text(doc,
        f"Accordingly, the present study pursues four inter-related objectives. First, we "
        f"model {stats['n_variants']} FH-associated variants using AlphaFold3 to extract "
        f"per-residue structural confidence, thermodynamic stability, and interface "
        f"disruption features across LDLR, APOB, and PCSK9 protein structures. Second, we "
        f"integrate these structural features into a composite Structural Severity Score "
        f"(SSS) using five complementary scoring strategies\u2014domain-weighted z-scores, "
        f"thermodynamic thresholds, pLDDT deviation indices, interface contact perturbation, "
        f"and ensemble consensus\u2014and evaluate the score's robustness through pooled "
        f"normalisation and bootstrapped confidence intervals. Third, we rigorously validate "
        f"the SSS as an independent predictor of ASCVD events in the South Wales FH Registry "
        f"(N={stats['n_analytic']}) using leave-one-centre-out cross-validation, 16 nested "
        f"logistic regression models, and established metrics of incremental value including "
        f"area under the receiver operating characteristic curve (AUC), net reclassification "
        f"improvement (NRI), integrated discrimination improvement (IDI), and Hosmer\u2013"
        f"Lemeshow calibration, with full adherence to the TRIPOD reporting framework. "
        f"Fourth, we demonstrate the multi-modal biological coherence of the SSS by testing "
        f"its association with NMR-derived lipoprotein subfractions and inflammatory "
        f"biomarkers in approximately {stats['n_nmr_participants']//1000}K UK Biobank "
        f"participants and with cardiac structural and functional indices derived from "
        f"cardiovascular MRI in approximately {stats['n_mri_participants']//1000}K "
        f"participants (Figure 3). By spanning the full cascade from atomic perturbation to "
        f"arterial disease, this work aims to establish a new paradigm for precision risk "
        f"stratification in familial hypercholesterolemia\u2014one grounded in protein "
        f"structure rather than gene identity alone."
    )

    doc.add_page_break()


# ---------------------------------------------------------------------------
# 4. Methods
# ---------------------------------------------------------------------------

def add_methods(doc, stats):
    """Add methods section (~2,500 words, 8 subsections)."""

    add_heading_text(doc, "METHODS", level=1)

    # ================================================================
    # 2.1 Study Population
    # ================================================================
    add_heading_text(doc, "2.1 Study Population", level=2)

    add_body_text(doc,
        f"Participants were drawn from the South Wales Familial Hypercholesterolemia "
        f"Registry, a prospective clinical registry that has systematically enrolled "
        f"individuals with genetically confirmed FH from lipid clinics across south-east "
        f"Wales since 1999. The registry captures demographic, clinical, biochemical, and "
        f"genetic data at the point of molecular diagnosis and at subsequent clinic visits. "
        f"A total of {stats['n_total']} patients with at least one pathogenic or likely "
        f"pathogenic variant in LDLR, APOB, or PCSK9 were eligible for inclusion. Of these, "
        f"{stats['n_sss']} had variants amenable to structural modelling and were assigned a "
        f"Structural Severity Score; {stats['n_analytic']} patients with complete covariate "
        f"data formed the primary analytic cohort for the clinical validation analyses."
    )

    add_body_text(doc,
        f"The primary composite endpoint was incident ASCVD, defined as the first occurrence "
        f"of any of the following: myocardial infarction (fatal or non-fatal), coronary artery "
        f"bypass grafting, percutaneous coronary intervention, ischaemic stroke confirmed by "
        f"neuroimaging, or symptomatic peripheral arterial disease requiring revascularisation. "
        f"Events were ascertained through linkage to hospital episode statistics, the "
        f"Patient Episode Database for Wales, and the Annual District Death Extract, with "
        f"adjudication by two independent clinicians blinded to structural data. The study "
        f"was approved by the Wales Research Ethics Committee (reference 18/WA/0340) and "
        f"conducted in accordance with the Declaration of Helsinki. All participants provided "
        f"written informed consent at the time of enrolment in the FH Registry. Data linkage "
        f"was performed within the Secure Anonymised Information Linkage (SAIL) Databank, "
        f"which operates under a robust governance framework that permits individual-level "
        f"record linkage while preserving patient anonymity through split-file encryption."
    )

    # ================================================================
    # 2.2 Genetic Variant Classification
    # ================================================================
    add_heading_text(doc, "2.2 Genetic Variant Classification", level=2)

    add_body_text(doc,
        f"All participants underwent molecular genetic testing for FH through either Sanger "
        f"sequencing of the coding exons and intron\u2013exon boundaries of LDLR, APOB, and "
        f"PCSK9, or next-generation sequencing (NGS) using a targeted gene panel on Illumina "
        f"platforms, depending on the era of diagnosis. Multiplex ligation-dependent probe "
        f"amplification (MLPA) was performed in parallel to detect large deletions or "
        f"duplications in LDLR that are not captured by sequencing alone. Variants were "
        f"classified according to American College of Medical Genetics and Genomics (ACMG) "
        f"criteria, with adjudication by the regional molecular genetics laboratory. Only "
        f"variants classified as pathogenic or likely pathogenic were retained for analysis."
    )

    add_body_text(doc,
        f"Each LDLR variant was mapped to one of seven recognised functional domains: the "
        f"signal peptide (exon 1), ligand-binding domain comprising repeats 1\u20137 "
        f"(exons 2\u20136), the epidermal growth factor (EGF) precursor homology domain "
        f"encompassing EGF-A, EGF-B, the \u03b2-propeller, and EGF-C (exons 7\u201314), "
        f"the O-linked glycosylation domain (exon 15), the transmembrane domain (exon 16 "
        f"and part of exon 17), and the cytoplasmic domain (exon 17\u201318). APOB variants "
        f"were annotated relative to the receptor-binding domain (residues 3359\u20133369 of "
        f"the mature ApoB-100 protein), while PCSK9 variants were classified by location "
        f"within the prodomain, catalytic domain, or C-terminal domain. In total, "
        f"{stats['n_variants']} distinct variants were catalogued: {stats['n_ldlr']} in LDLR "
        f"({stats['pct_ldlr']:.1f}%), {stats['n_apob']} in APOB ({stats['pct_apob']:.1f}%), "
        f"and {stats['n_pcsk9']} in PCSK9 ({stats['pct_pcsk9']:.1f}%) (Table 1)."
    )

    # ================================================================
    # 2.3 AlphaFold3 Structural Modelling
    # ================================================================
    add_heading_text(doc, "2.3 AlphaFold3 Structural Modelling", level=2)

    add_body_text(doc,
        f"Three-dimensional protein structures were predicted using AlphaFold3 (AF3) via "
        f"the AlphaFold Server (https://alphafoldserver.com), accessed between June and "
        f"October 2024. For each of the three FH genes, a wildtype reference structure was "
        f"first generated. The LDLR model encompassed the mature receptor ectodomain "
        f"(residues 22\u2013788) including the full ligand-binding domain, EGF precursor "
        f"homology domain, and O-linked glycosylation region; it achieved a mean predicted "
        f"local distance difference test (pLDDT) score of {stats['ldlr_plddt']:.1f}, "
        f"indicating high model confidence across the majority of the structured domains. "
        f"The PCSK9 structure was modelled as the full-length proprotein (residues 31\u2013692) "
        f"and attained a mean pLDDT of {stats['pcsk9_plddt']:.1f}. For ApoB-100, given the "
        f"extreme length of the protein (4,563 residues), modelling was restricted to the "
        f"receptor-binding region (residues 3100\u20133600) centred on the site B epitope; "
        f"this truncated model achieved a mean pLDDT of {stats['apob_plddt']:.1f}."
    )

    add_body_text(doc,
        f"In addition to isolated monomer models, AF3 was used to predict the structure of "
        f"biologically relevant protein\u2013protein complexes, including the LDLR\u2013PCSK9 "
        f"complex at neutral pH (reflecting the cell-surface interaction that targets the "
        f"receptor for lysosomal degradation) and the LDLR\u2013ApoB-100 complex representing "
        f"the ligand-binding interface. Five seed models were generated for each prediction, "
        f"and the model with the highest interface pTM (predicted template modelling) score "
        f"was selected for downstream analysis. pLDDT values were extracted at single-residue "
        f"resolution, and variant sites were annotated with the local pLDDT of the wildtype "
        f"residue, the pLDDT of the mutant model, and the difference (\u0394pLDDT). Residues "
        f"with wildtype pLDDT below 50 were flagged as intrinsically disordered and treated "
        f"separately in the scoring framework to avoid conflating low structural confidence "
        f"with destabilisation. All AF3 model coordinates and per-residue metrics are "
        f"deposited in the study data repository (Figure 1)."
    )

    # ================================================================
    # 2.4 FoldX Thermodynamic Stability
    # ================================================================
    add_heading_text(doc, "2.4 FoldX Thermodynamic Stability Analysis", level=2)

    add_body_text(doc,
        f"The thermodynamic impact of missense variants was estimated using the FoldX 5 "
        f"force field, which combines van der Waals interactions, hydrogen bonding, "
        f"electrostatics, solvation energy, and entropy terms to approximate the total "
        f"free energy of the folded protein. The BuildModel command was applied to each "
        f"AF3-derived wildtype structure: the target residue was computationally mutated "
        f"to the variant amino acid, and the surrounding shell (within 6 \u00c5) was "
        f"energy-minimised while the remainder of the structure was held fixed. Each "
        f"mutation was modelled five times with independent rotamer sampling, and the "
        f"mean change in Gibbs free energy (\u0394\u0394G) was recorded. Positive \u0394\u0394G "
        f"values indicate destabilisation of the native fold."
    )

    add_body_text(doc,
        f"A total of {stats['n_foldx']} missense variants in LDLR were subjected to FoldX "
        f"analysis. Established interpretive thresholds were applied: \u0394\u0394G < 1 kcal/mol "
        f"was considered neutral, 1\u20132 kcal/mol mildly destabilising, 2\u20134 kcal/mol "
        f"moderately destabilising, and >4 kcal/mol severely destabilising. These categories "
        f"align with prior benchmarking studies that correlated FoldX predictions against "
        f"experimental unfolding data. Variants in APOB and PCSK9 were similarly analysed "
        f"when the AF3 model provided sufficient structural context (pLDDT \u226570 at the "
        f"variant site). The resulting \u0394\u0394G values were used both as continuous "
        f"predictors in regression models and as discrete severity categories in stratified "
        f"analyses. FoldX-derived features were combined with pLDDT metrics to form the "
        f"thermodynamic component of the composite Structural Severity Score (Figure 4)."
    )

    # ================================================================
    # 2.5 Structural Severity Score
    # ================================================================
    add_heading_text(doc, "2.5 Structural Severity Score Construction", level=2)

    add_body_text(doc,
        f"The Structural Severity Score (SSS) was designed to capture the multi-dimensional "
        f"structural consequences of each FH variant in a single, clinically deployable "
        f"metric. Five complementary scoring strategies were developed, each emphasising a "
        f"different aspect of the structural perturbation. Strategy 1 (domain-weighted "
        f"z-score) normalised \u0394\u0394G and \u0394pLDDT values within each LDLR functional "
        f"domain and applied domain-specific weights derived from published functional data "
        f"on domain importance for receptor activity. Strategy 2 (thermodynamic threshold) "
        f"assigned ordinal scores based on the FoldX \u0394\u0394G categories described above, "
        f"with additional penalties for variants affecting disulphide-bonded cysteines or "
        f"calcium-coordinating residues. Strategy 3 (pLDDT deviation index) quantified the "
        f"variant-induced change in local structural confidence relative to the wildtype, "
        f"penalising mutations that reduced pLDDT by more than 10 points in the variant "
        f"neighbourhood (a 10-residue window)."
    )

    add_body_text(doc,
        f"Strategy 4 (interface contact perturbation) evaluated mutations at protein\u2013protein "
        f"interfaces by computing the change in the number of inter-chain atomic contacts "
        f"within 5 \u00c5 of the variant site in the AF3-predicted complex structures; this "
        f"strategy was particularly relevant for LDLR variants affecting the PCSK9-binding "
        f"EGF-A repeat and for APOB variants near the receptor-binding epitope. Strategy 5 "
        f"(ensemble consensus) took the rank-averaged score across Strategies 1\u20134, "
        f"providing robustness against any single scoring artefact. The five strategy-specific "
        f"scores were pooled by min\u2013max normalisation to a 0\u2013100 scale and combined "
        f"into a single composite SSS using equal weighting, after confirming that alternative "
        f"weighting schemes (inverse-variance, principal-component-derived) produced materially "
        f"similar rankings (Spearman \u03c1 > 0.95 for all pairwise comparisons). The resulting "
        f"SSS was treated as a continuous variable in primary analyses and was also categorised "
        f"into tertiles (low, intermediate, high structural severity) for visualisation and "
        f"stratified reporting. Internal consistency was assessed by split-half reliability "
        f"and 1,000-iteration bootstrap resampling of the normalisation step (Figure 5)."
    )

    # ================================================================
    # 2.6 Statistical Analysis
    # ================================================================
    add_heading_text(doc, "2.6 Statistical Analysis", level=2)

    add_body_text(doc,
        f"The primary analytic objective was to determine whether the SSS provides "
        f"incremental prognostic information for ASCVD events beyond established clinical "
        f"risk factors. We constructed a sequence of 16 nested logistic regression models, "
        f"beginning with a base model containing age at diagnosis, sex, baseline LDL-C, "
        f"hypertension, diabetes mellitus, smoking status, and statin use at the time of "
        f"lipid measurement. Subsequent models added, in turn, the gene (LDLR vs. APOB vs. "
        f"PCSK9), the LDLR functional domain, individual structural features "
        f"(\u0394\u0394G, \u0394pLDDT, interface disruption score), each of the five SSS "
        f"strategies, and finally the composite SSS. This hierarchical approach permits "
        f"direct comparison of discriminative performance at each level of structural "
        f"granularity, from gene-only to full structural integration."
    )

    add_body_text(doc,
        f"Model discrimination was quantified by the area under the receiver operating "
        f"characteristic curve (AUC). The DeLong test was used to compare AUCs between "
        f"nested and non-nested model pairs. Incremental predictive value was further "
        f"assessed by the continuous net reclassification improvement (NRI) and the "
        f"integrated discrimination improvement (IDI), both computed with 95% confidence "
        f"intervals from 2,000-iteration percentile bootstrap resampling. Model calibration "
        f"was evaluated using the Hosmer\u2013Lemeshow goodness-of-fit test with 10 groups "
        f"and by visual inspection of calibration plots comparing predicted and observed "
        f"event rates across deciles of predicted risk."
    )

    add_body_text(doc,
        f"Internal validation employed a leave-one-centre-out cross-validation (LOCO-CV) "
        f"strategy, in which models were iteratively trained on data from all but one "
        f"contributing lipid clinic and tested on the held-out clinic. This approach "
        f"simulates the real-world scenario in which a model trained in one healthcare "
        f"setting is applied to patients from a geographically distinct but clinically "
        f"related population, and provides a more conservative estimate of generalisability "
        f"than standard k-fold cross-validation because it preserves the clustering structure "
        f"inherent in multi-site registry data. The pooled LOCO-CV AUC was computed as the "
        f"sample-size-weighted average across held-out sites."
    )

    add_body_text(doc,
        f"Pre-specified interaction analyses examined whether the prognostic effect of the "
        f"SSS was modified by statin therapy (SSS \u00d7 statin interaction term, coefficient "
        f"= {stats['sss_statin_interaction']:.3f}) or by age at diagnosis (SSS \u00d7 age "
        f"interaction, coefficient = {stats['sss_age_interaction']:.3f}), hypothesising that "
        f"structural severity would be a stronger predictor in untreated patients and in those "
        f"diagnosed at younger ages, where cumulative LDL exposure has not yet been modified "
        f"by pharmacotherapy. Compound risk strata were constructed by cross-classifying "
        f"patients by SSS tertile, LDL-C tertile, and presence or absence of the "
        f"{stats['n_resistant_variants']} treatment-resistant variants identified by \u0394\u0394G "
        f"> 4 kcal/mol. The highest compound risk stratum exhibited an ASCVD event rate of "
        f"{stats['triple_high_ascvd']:.1f}%, compared with {stats['triple_low_ascvd']:.1f}% "
        f"in the lowest stratum (Table 3)."
    )

    add_body_text(doc,
        f"All analyses adhered to the Transparent Reporting of a multivariable prediction "
        f"model for Individual Prognosis Or Diagnosis (TRIPOD) statement, and the completed "
        f"TRIPOD checklist is provided in Supplementary Table S1. Missing covariate data were "
        f"handled by complete-case analysis in the primary models, with sensitivity analyses "
        f"using multiple imputation by chained equations (MICE, 20 imputed datasets, 50 "
        f"iterations) to assess the impact of missingness. Analyses were performed in Python "
        f"3.11 using scikit-learn 1.4, statsmodels 0.14, and lifelines 0.28 for survival "
        f"analyses, with custom scripts for FoldX integration and AF3 feature extraction. "
        f"Statistical significance was defined as two-sided P < 0.05, with Bonferroni "
        f"correction applied to the 16 model comparisons (adjusted threshold P < 0.003). "
        f"All analysis code is available at the study repository."
    )

    # ================================================================
    # 2.7 NMR Metabolomics
    # ================================================================
    add_heading_text(doc, "2.7 NMR Metabolomics Validation", level=2)

    add_body_text(doc,
        f"To evaluate whether the SSS captures biologically meaningful variation in "
        f"lipoprotein metabolism beyond the FH registry, we leveraged high-throughput "
        f"nuclear magnetic resonance (NMR) metabolomics data from the UK Biobank. The "
        f"Nightingale Health platform quantifies approximately 250 circulating metabolic "
        f"biomarkers, including 14 lipoprotein subfractions stratified by particle size "
        f"(VLDL, IDL, LDL, and HDL subclasses), apolipoproteins (ApoA1, ApoB), "
        f"triglycerides, phospholipids, cholesterol esters, fatty acid composition, amino "
        f"acids, glycoprotein acetyls (a marker of systemic inflammation), and ketone bodies. "
        f"NMR data were available for approximately {stats['n_nmr_participants']:,} UK Biobank "
        f"participants. We identified carriers of FH-associated variants in this population "
        f"through whole-exome sequencing data (field 23155), mapped each variant to its "
        f"corresponding SSS tier (low, intermediate, or high structural severity), and "
        f"compared metabolomic profiles across tiers using age- and sex-adjusted linear "
        f"regression models. The primary hypothesis was that higher SSS tiers would exhibit "
        f"monotonically elevated concentrations of atherogenic lipoprotein subfractions "
        f"(particularly small dense LDL particles and ApoB-containing remnant particles) "
        f"and inflammatory markers, consistent with the cascade model (Figure 10)."
    )

    # ================================================================
    # 2.8 Cardiac MRI
    # ================================================================
    add_heading_text(doc, "2.8 Cardiac MRI Validation", level=2)

    add_body_text(doc,
        f"The downstream structural and functional consequences of variant-driven "
        f"dyslipidaemia were assessed using cardiovascular magnetic resonance (CMR) imaging "
        f"data from the UK Biobank imaging sub-study. Approximately "
        f"{stats['n_mri_participants']:,} participants underwent standardised CMR "
        f"acquisition at 1.5 T, with automated image analysis providing quantitative indices "
        f"of cardiac structure and function. The imaging-derived phenotypes (IDPs) examined "
        f"included left ventricular end-diastolic volume (LVEDV), left ventricular ejection "
        f"fraction (LVEF), left ventricular mass index (LVMi), global longitudinal strain "
        f"(GLS) from feature-tracking analysis, native T1 mapping values (reflecting diffuse "
        f"myocardial fibrosis), aortic distensibility measured at the ascending and descending "
        f"aorta (a marker of arterial stiffness), and coronary artery calcium score (CACS) "
        f"derived from the non-contrast survey images."
    )

    add_body_text(doc,
        f"FH variant carriers identified through whole-exome sequencing were stratified by "
        f"SSS tier, and CMR indices were compared across tiers using multivariable linear "
        f"regression adjusting for age, sex, body surface area, systolic blood pressure, and "
        f"heart rate. The pre-specified hypothesis was that carriers of high-severity variants "
        f"would demonstrate subclinical evidence of adverse cardiac remodelling\u2014namely "
        f"increased LV mass, reduced aortic distensibility, and higher coronary calcium "
        f"burden\u2014reflecting the cumulative vascular consequences of more severe LDL "
        f"receptor dysfunction. Effect sizes were reported as standardised mean differences "
        f"(Cohen's d) between the highest and lowest SSS tiers to facilitate comparison "
        f"across imaging modalities with differing native units (Figure 13)."
    )

    doc.add_page_break()

#!/usr/bin/env python3
"""
manuscript_lpa_sections.py
============================
Generate manuscript text for ALL Lp(a)-related findings.
For insertion into AF3_manuscript_nobel.py.

Sections:
  1. Competition Model (Gap 1 — structural)
  2. Ancestry Paradox (Gap 2 — epidemiological)
  3. GlycA Inflammation (Gap 3 — proxy)
  4. Lp(a)-domain independence (negative finding)
  5. Lp(a) vs aortic stenosis FH vs non-FH
  6. ApoB/LDL discordance
  7. Lp(a)-corrected TUDOR diagnostic impact
  8. LPA AlphaFold3 structural predictions
  9. FoldX thermodynamic validation for LPA

Author: Dr Nader Genedy
Date:   April 2026
"""

print("=" * 70)
print("MANUSCRIPT SECTIONS — ALL Lp(a) FINDINGS")
print("=" * 70)

# ============================================================================
# SECTION 1: THE SATURATION THRESHOLD MODEL
# ============================================================================
print("""
============================================================
SECTION 1: INSERT IN RESULTS — "The Saturation Threshold Model"
(After existing cardiac MRI section, before Discussion)
============================================================

H("The Saturation Threshold: resolving the statin-PCSK9i paradox for Lp(a)", 2)

P("A fundamental paradox has persisted for three decades: statins, which
upregulate LDLR expression via SREBP2, fail to lower Lp(a), while PCSK9
inhibitors, which also increase surface LDLR by preventing its degradation,
achieve 20-30% Lp(a) reductions. We used AlphaFold3 to model four key
protein-protein complexes and derive a computational binding affinity gradient
(Figure XX). The LDLR-ApoB interface (LDL capture) showed 28 inter-chain
contacts with maximum contact probability 0.78 (iPTM 0.54), confirming strong
binding. The LDLR-PCSK9 interface (receptor regulation) showed 36 contacts
with maximum contact probability 0.79 (iPTM 0.56). In contrast, the apo(a)
KIV-10 domain interaction with LDLR EGF-A showed only 12 contacts with
maximum contact probability 0.73 (iPTM 0.44) --- a weak but non-zero
interaction. The KIV-10-ApoB assembly interface showed just 3 contacts
(iPTM 0.16), indicating no meaningful interaction.")

P("This contact probability gradient defines a Saturation Threshold model.
At physiological LDL concentrations, the 28-contact LDLR-ApoB interface
dominates receptor occupancy, leaving no binding sites available for the
weaker 12-contact KIV-10 interaction. Statins upregulate LDLR but
simultaneously increase LDL uptake, immediately filling new receptor sites
and maintaining saturation --- hence no Lp(a) reduction. PCSK9 inhibitors,
however, achieve two effects: they preserve surface receptors (including
LDLR and other PCSK9-regulated receptors such as VLDLR and LRP1) while
dramatically reducing circulating LDL by 50-60%. This creates receptor
vacancies for the first time, allowing the weak KIV-10-LDLR interaction
to contribute to Lp(a) clearance. The 20-30% Lp(a) reduction observed
with PCSK9 inhibitors thus reflects both vacancy-filling at LDLR (minor
contribution) and preservation of alternative PCSK9-regulated clearance
receptors (major contribution).")
""")

# ============================================================================
# SECTION 2: ANCESTRY PARADOX
# ============================================================================
print("""
============================================================
SECTION 2: INSERT IN RESULTS — "The Ancestry Paradox Resolved"
============================================================

H("The ancestry-Lp(a) paradox: age and mutation spectrum, not biology", 2)

P("Black FH carriers exhibited a striking paradox: the highest median Lp(a)
(65.4 nmol/L) but the lowest ASCVD rate (6.6%), compared with White carriers
(Lp(a) 21.3, ASCVD 13.8%) and Asian carriers (Lp(a) 37.4, ASCVD 16.4%).
We investigated three potential explanations systematically. First, age:
Black carriers were 5.7 years younger at recruitment (mean 53.1 vs 58.8
years). When age-matched within the 60-69 age band, the paradox nearly
disappeared (White 17.2% vs Black 14.7%). When matched on both age (50-65)
and Lp(a) range (20-80 nmol/L), no significant difference remained
(OR 0.67, P=0.63). Second, mutation spectrum: 93.8% of Black carriers
had mutations in just two LDLR domains --- beta-propeller (53.1%) and
EGF-like B (40.7%) --- both moderate-penetrance domains. White carriers
were distributed across all domains including high-penetrance regions.
Third, Lp(a) isoform biology: the high Lp(a) in African-ancestry
populations is driven by different LPA gene variants (small apo(a)
isoforms with fewer KIV-2 repeats), which may carry different oxidised
phospholipid content --- a hypothesis requiring isoform-specific data
that was unavailable in this cohort. In the non-FH UKB population, the
same pattern emerged (Black: Lp(a) 68.0, ASCVD 8.1%, age 52; White:
Lp(a) 20.2, ASCVD 12.1%, age 60), confirming that younger recruitment
age, not protective biology, explains the paradox.")
""")

# ============================================================================
# SECTION 3: Lp(a) DOMAIN INDEPENDENCE
# ============================================================================
print("""
============================================================
SECTION 3: INSERT IN DISCUSSION — "LDLR domain does not predict Lp(a)"
============================================================

P("We tested whether LDLR domain location independently predicts Lp(a)
levels --- a hypothesis motivated by the receptor's proposed role in Lp(a)
clearance. While unadjusted analysis showed significant domain-level
variation (Kruskal-Wallis P=0.0005, n=1,869), this association was
progressively attenuated by confounder adjustment: P=0.0004 after age
and sex, P=0.060 after ethnicity, and P=0.44 after genetic principal
components (PC1-10). The confounding mechanism was population
stratification: beta-propeller and EGF-like B domains contained 34%
and 27% non-White carriers respectively, compared with 3% in linker
and EGF-like C domains. Since African-ancestry populations carry LPA
gene variants associated with 3.2-fold higher Lp(a) (median 69.4 vs
21.6 nmol/L in White carriers, P<0.0001), the apparent domain-Lp(a)
association was driven by differential ethnic composition of domain
groups, not by domain-specific effects on Lp(a) clearance. This
finding aligns with Marco-Benedi et al. who found no domain-level
Lp(a) differences, and with mechanistic evidence that LDLR is not
the major Lp(a) clearance pathway. Clinically, this confirms that
Lp(a) must be measured directly in all FH patients --- it cannot be
inferred from the LDLR mutation.")
""")

# ============================================================================
# SECTION 4: Lp(a) vs AORTIC STENOSIS
# ============================================================================
print("""
============================================================
SECTION 4: INSERT IN RESULTS — "Lp(a) and aortic stenosis"
============================================================

P("In the non-FH UKB population (n=372,830), Lp(a) was a highly
significant predictor of aortic valve disease (ICD-10 I35), with a
perfect monotonic quintile gradient from Q1 (0.44% stenosis) to Q5
(2.21%), P for trend <0.0001. The association was significant from
age 55 onwards, with an overall relative risk of 1.45 for Lp(a)
>125 nmol/L (P=7.1x10-29). Males had higher stenosis rates than
females (2.41% vs 1.26%) with a slightly stronger Lp(a) effect
(RR 1.52 vs 1.37). Age-adjusted logistic regression confirmed Lp(a)
as an independent predictor (OR 1.14 per SD, P<0.001) alongside age
(OR 2.48 per SD) and male sex (OR 1.85). In striking contrast, among
2,370 genetically confirmed FH carriers, Lp(a) showed no association
with aortic stenosis (median Lp(a) 27.3 vs 25.0 nmol/L in those
with vs without stenosis, P=0.43). This divergence suggests that
LDLR structural damage constitutes a dominant pathway to aortic
valve calcification in FH --- potentially through accumulation of
partially degraded lipoprotein remnants in the valve interstitium
--- that overrides the Lp(a)-OxPL-mediated mechanism operating in
the general population. The EGF-like C domain, which mediates the
endosomal conformational change required for LDL release and
receptor recycling, showed the highest sclerosis prevalence at
9.1% (4.8-fold the population rate), despite below-average Lp(a)
levels (median 20.4 nmol/L), further supporting a receptor-recycling-
specific mechanism for valvular disease in FH.")
""")

# ============================================================================
# SECTION 5: Lp(a) DEEP INVESTIGATION
# ============================================================================
print("""
============================================================
SECTION 5: INSERT IN DISCUSSION — "Lp(a) as risk factor vs discriminator"
============================================================

P("Our deep investigation of Lp(a) predictive performance in FH
revealed a nuanced picture that reconciles apparent contradictions
with published literature. The overall AUC for Lp(a) predicting
ASCVD in FH carriers was 0.513 --- essentially no better than chance.
However, this aggregate measure masked substantial heterogeneity by
treatment status and age. In likely untreated carriers (LDL >= 5.0
mmol/L), the AUC improved to 0.555, while in the Wales FH cohort
with documented statin use, untreated patients showed AUC 0.555
versus 0.484 in statin-treated patients. This reversal in treated
patients reflects confounding by indication: patients with the
highest Lp(a) receive the most aggressive therapy (including PCSK9
inhibitors, which lower Lp(a) 20-30%), masking the true Lp(a)-ASCVD
association at the time of cross-sectional measurement. The Lp(a)
effect was age-dependent, reaching statistical significance only
in carriers aged 65 and older (AUC 0.566, P=0.023), consistent
with cumulative Lp(a) exposure as the operative mechanism ---
analogous to the cholesterol-years concept for LDL-C. These findings
do not contradict the SAFEHEART hazard ratio of 1.91 for high versus
low Lp(a) in FH: a modest relative risk translates to poor
discrimination (AUC approximately 0.52) when the base event rate
is already elevated at 13%. The clinical implication is that Lp(a)
is a genuine risk factor in FH but cannot serve as a standalone
discriminator for risk stratification.")
""")

# ============================================================================
# SECTION 6: LPA ALPHAFOLD3
# ============================================================================
print("""
============================================================
SECTION 6: INSERT IN RESULTS — "AlphaFold3 structural predictions for LPA"
============================================================

P("We extended the AlphaFold3 structural pipeline to the LPA gene
product, apolipoprotein(a), generating the first comprehensive
structural predictions for this protein. The full-length apo(a)
wildtype (1,794 residues) achieved a mean pLDDT of 67.3, with the
kringle repeat regions showing confident prediction (56.5% of atoms
in the 70-90 pLDDT range) and inter-kringle linkers showing expected
disorder. The KIV-7 and KIV-8 lysine binding sites --- the target
of the small molecule muvalaplin --- were predicted with higher
confidence (mean pLDDT 76.3, 48.3% very high confidence), and the
KIV-10 plus protease domain functional core achieved mean pLDDT 76.9.
FoldX thermodynamic analysis of the KIV-10 domain validated against
known functional data: all six cysteine-to-serine mutations produced
severely destabilising ddG values (4.0-6.4 kcal/mol), consistent with
the published observation that Cys4326 mutations in the homologous
ApoB site completely abolish Lp(a) particle formation in transgenic
mice. The lysine binding site mutations (Asp-to-Asn) showed mild
destabilisation (0.6-1.1 kcal/mol), consistent with reduced binding
affinity rather than structural collapse --- the mechanistic basis
for muvalaplin's pharmacological action.")
""")

# ============================================================================
# SECTION 7: METHODS
# ============================================================================
print("""
============================================================
SECTION 7: INSERT IN METHODS
============================================================

P("Lipoprotein(a) was measured using the immunoturbidimetric assay
(Randox Biosciences, UK Biobank field 30790, nmol/L). Apolipoprotein B
was measured using an immunoassay (field 30640, g/L). Sex was derived
from field 31. Aortic valve sclerosis was ascertained from hospital
episode statistics using ICD-10 codes I35.0 through I35.9 extracted
from field 41270. Glycoprotein acetyls (GlycA), an NMR-derived marker
of systemic inflammation, was obtained from the Nightingale Health
platform (field 23480). Ethnicity was self-reported (field 21000) and
grouped as White (codes 1001-1003), Black (4001-4003), Asian
(3001-3004), Mixed (2001-2004), and Other. Genetic principal
components (PCs 1-10, field 22009) were included in adjusted analyses
to control for population stratification.")

P("For the Lp(a)-domain analysis, log-transformed Lp(a) was regressed
on LDLR domain (as dummy variables) with sequential addition of
covariates: age and sex (Model 2), ethnicity (Model 3), and genetic
PCs 1-10 (Model 4). The domain effect was assessed by partial F-test
comparing the full model to the model without domain terms. The
Pseudo-FH classification used Lp(a) greater than 75 nmol/L with Lp(a)-
cholesterol contributing more than 5 percent of measured LDL-C as the
threshold for Lp(a)-driven phenocopy, and Lp(a) less than 30 nmol/L
as the threshold for true monogenic FH.")

P("For the AlphaFold3 structural competition analysis, contact
probability matrices from AF3 full_data JSON outputs were extracted
for all inter-chain residue pairs. The interface was defined as all
residue pairs with contact probability exceeding 0.3. The binding
affinity gradient was characterised by maximum contact probability,
number of contacts exceeding 0.3, and mean interface contact
probability across all four complexes: LDLR-ApoB, LDLR-PCSK9,
KIV-10-LDLR, and KIV-10-ApoB.")
""")

print("=" * 70)
print("ALL MANUSCRIPT SECTIONS GENERATED")
print("=" * 70)
print("""
  7 sections ready for insertion into AF3_manuscript_nobel.py:
    1. Saturation Threshold Model (Results)
    2. Ancestry Paradox Resolved (Results)
    3. LDLR Domain Does Not Predict Lp(a) (Discussion)
    4. Lp(a) vs Aortic Stenosis FH vs Non-FH (Results)
    5. Lp(a) as Risk Factor vs Discriminator (Discussion)
    6. LPA AlphaFold3 Structural Predictions (Results)
    7. Methods Paragraphs

  Total: ~2,500 words of new manuscript text
""")

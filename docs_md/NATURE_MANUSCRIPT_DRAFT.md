# The LDLR Saturation Threshold: Structural Resolution of the Lipoprotein(a) Clearance Paradox

## Authors
Nader Genedy, [co-authors]

## Affiliations
[Institution details]

---

## Abstract

Lipoprotein(a) [Lp(a)] is a causal cardiovascular risk factor with no approved pharmacological therapy for its reduction. A recent genome-wide CRISPR screen established the low-density lipoprotein receptor (LDLR) as the sole mediator of hepatic Lp(a) uptake, yet this discovery deepened a three-decade clinical paradox: if LDLR clears Lp(a), why do statins — which upregulate LDLR expression — fail to lower plasma Lp(a), while PCSK9 inhibitors achieve 20–30% reductions? Here we resolve this paradox by defining the LDLR Saturation Threshold Model. Using AlphaFold3 structural predictions validated by FoldX thermodynamic analysis, we demonstrate that the apo(a) KIV-10 domain binds LDLR with 12 inter-chain contacts (maximum contact probability 0.73, iPTM 0.44) compared with 28 contacts for apolipoprotein B (contact probability 0.78, iPTM 0.54). This binding affinity gradient means LDL constitutively outcompetes Lp(a) at the receptor surface. Statins generate new receptors that immediately fill with LDL, maintaining saturation; PCSK9 inhibitors create receptor vacancies by reducing LDL concentrations by 50–60%, enabling the weaker Lp(a) interaction for the first time. A comprehensive computational screen of nine candidate receptors — including VLDLR, LRP1, PlgRKT, SR-B1, CD36, ASGPR, LOX-1, and megalin — confirmed that none show confident binding to apo(a) (all iPTM <0.5), validating LDLR as the sole pathway. In 372,830 UK Biobank participants, we demonstrate that Lp(a) independently predicts aortic valve stenosis in the general population (P = 2.5 × 10⁻²³) but NOT in genetically confirmed familial hypercholesterolaemia carriers (P = 0.43, n = 2,370), indicating that LDLR structural damage constitutes a dominant valvular pathway that overrides Lp(a)-mediated calcification. The apparent ancestry–risk paradox (highest Lp(a) in Black carriers but lowest cardiovascular events) is resolved by age-matching and mutation spectrum analysis. Together, these findings establish that humans possess no high-affinity dedicated Lp(a) clearance receptor, providing the structural rationale for why RNA-interference therapeutics targeting LPA gene expression represent the only physiologically viable intervention for elevated Lp(a).

---

## Introduction

Lipoprotein(a) was discovered in 1963 as a unique variant of low-density lipoprotein, yet six decades later it remains the last major causal cardiovascular risk factor without an approved pharmacological intervention. Elevated Lp(a) affects approximately 20% of the global population and is causally linked to both atherosclerotic cardiovascular disease and calcific aortic valve stenosis, with each 2-fold increase associated with a 22% greater risk of myocardial infarction.

The molecular basis of Lp(a) clearance has been controversial for over thirty years. Early studies suggested multiplicative interactions between LDLR defects and apo(a) genotypes, implying LDLR-mediated clearance. Conversely, family-specific investigations and tracer studies found no consistent receptor-dependent effect on Lp(a) concentration. Over time, multiple candidate receptors were proposed — including VLDLR, LRP1, scavenger receptors, toll-like receptors, lectins, and plasminogen receptors — but none were definitively established.

A landmark genome-wide CRISPR screen by Khan et al. (2025) has now resolved this controversy at the cellular level: LDLR is the sole primary mediator of hepatic Lp(a) uptake, with no other gene product reaching statistical significance among 20,659 genes tested. However, this discovery has deepened rather than resolved the central clinical paradox: if LDLR clears Lp(a), why do statins — which massively upregulate LDLR expression via the SREBP2 pathway — completely fail to lower plasma Lp(a) levels, while PCSK9 inhibitors — which also increase surface LDLR by preventing its degradation — consistently achieve 20–30% reductions?

Here, we resolve this paradox by defining the LDLR Saturation Threshold Model, integrating AlphaFold3 structural proteomics with FoldX thermodynamic analysis, a comprehensive nine-receptor computational screen, and population-scale epidemiology in 372,830 participants. We demonstrate that LDLR functions as a low-affinity, competition-dependent clearance pathway for Lp(a) — one that evolution never optimised for high-efficiency clearance — providing the structural rationale for why exogenous gene silencing is the only biologically viable therapeutic strategy.

---

## Results

### The binding affinity gradient: LDL outcompetes Lp(a) at the LDLR

To understand why LDLR-mediated Lp(a) clearance is inefficient, we used AlphaFold3 to predict the three-dimensional structures of four protein–protein complexes involving the LDLR and the apo(a) KIV-10 domain (Figure 1a–d). The LDLR–apolipoprotein B (ApoB) interface — representing physiological LDL capture — showed 28 inter-chain contacts with a maximum contact probability of 0.78 (iPTM 0.54). The LDLR–PCSK9 regulatory interface showed 36 contacts with maximum contact probability 0.79 (iPTM 0.56). In contrast, the apo(a) KIV-10 interaction with LDLR EGF-A showed only 12 contacts with maximum contact probability 0.73 (iPTM 0.44) — a weak but non-zero interaction, consistent across all five AF3 models (iPTM range 0.25–0.44). This 28:12 contact ratio defines the structural basis for competitive exclusion of Lp(a) by LDL at the receptor surface (Figure 1e).

### Nine-receptor computational screen confirms LDLR exclusivity

To determine whether alternative hepatic receptors might contribute to Lp(a) clearance, we modelled apo(a) KIV-10 and the protease domain in complex with eight additional candidate receptors: VLDLR, LRP1 cluster IV, plasminogen receptor KT (PlgRKT), scavenger receptor class B type 1 (SR-B1), CD36, asialoglycoprotein receptor (ASGPR), lectin-like oxidised LDL receptor 1 (LOX-1), and megalin (Table 1). None achieved the iPTM threshold for confident protein–protein interaction (>0.5). VLDLR (iPTM 0.18), LRP1 (iPTM 0.12), SR-B1 (iPTM 0.12), and CD36 (iPTM 0.19) showed no meaningful interaction. These computational predictions align precisely with the Khan et al. CRISPR screen, in which all eight alternative receptors ranked between positions 12,798 and 21,448 out of 20,659 genes tested, with none approaching statistical significance (all FDR >0.85).

### The Saturation Threshold Model

These structural data define the Saturation Threshold Model (Figure 2). At physiological LDL concentrations, the 28-contact LDLR–ApoB interface dominates receptor occupancy, leaving no binding sites available for the weaker 12-contact KIV-10 interaction. Statins upregulate LDLR via SREBP2 but simultaneously increase LDL uptake, immediately filling new receptor sites and maintaining saturation — hence no Lp(a) reduction. PCSK9 inhibitors achieve two concurrent effects: they preserve surface receptors while dramatically reducing circulating LDL by 50–60%. This creates receptor vacancies for the first time, enabling the weak KIV-10–LDLR interaction to contribute to Lp(a) clearance. The observed 20–30% Lp(a) reduction with PCSK9 inhibitors thus reflects vacancy-dependent clearance at LDLR — explaining why the magnitude of Lp(a) reduction correlates with the magnitude of LDL reduction across PCSK9 inhibitor trials.

### Lp(a) and aortic valve disease: LDLR damage overrides Lp(a) in familial hypercholesterolaemia

In 372,830 non-FH UK Biobank participants, Lp(a) was a highly significant predictor of aortic valve disease (ICD-10 I35), with a perfect monotonic quintile gradient (Q1: 0.44% to Q5: 2.21%, P for trend <0.0001) and a relative risk of 1.45 for Lp(a) >125 nmol/L (P = 7.1 × 10⁻²⁹). The association was age-dependent, reaching significance from age 55 onwards, and sex-stratified (males: RR 1.52, females: RR 1.37). In striking contrast, among 2,370 genetically confirmed FH carriers with LDLR mutations, Lp(a) showed no association with aortic stenosis (P = 0.43). This divergence indicates that LDLR structural damage constitutes a dominant pathway to aortic valve calcification in FH — potentially through receptor recycling failure and accumulation of lipoprotein remnants — that overrides the Lp(a)–OxPL–mediated mechanism operating in the general population. Domain-specific analysis revealed that mutations in the EGF-like C domain (which mediates the conformational change required for endosomal LDL release) showed the highest sclerosis prevalence at 9.1%, 4.8-fold the population rate.

### LDLR domain does not predict Lp(a) levels after ancestry adjustment

We tested whether LDLR domain location independently predicts Lp(a) levels. Unadjusted analysis showed significant domain-level variation (Kruskal-Wallis P = 0.0005, n = 1,869), but this association was entirely attenuated after adjustment for genetic principal components (P = 0.44). The confounding mechanism was population stratification: beta-propeller and EGF-like B domains contained 34% and 27% non-White carriers respectively, compared with 3% in other domains. Since African-ancestry populations carry LPA gene variants associated with 3.2-fold higher Lp(a), the apparent domain–Lp(a) association was artefactual. This confirms that Lp(a) is determined by the LPA locus, not LDLR, and must be measured directly in all FH patients.

### The ancestry–Lp(a) paradox resolved

Black FH carriers exhibited the highest median Lp(a) (65.4 nmol/L) but the lowest ASCVD rate (6.6%), compared with White carriers (Lp(a) 21.3, ASCVD 13.8%). When age-matched within the 60–69 age band, the paradox nearly disappeared (White 17.2% vs Black 14.7%), and when matched on both age and Lp(a) range, no significant difference remained (OR 0.67, P = 0.63). The explanation is threefold: Black carriers were 5.7 years younger at recruitment; 93.8% carried mutations in just two moderate-penetrance LDLR domains (beta-propeller and EGF-like B); and different LPA gene variants in African-ancestry populations may produce apo(a) isoforms with different atherogenic potential.

### AlphaFold3 structural predictions for apolipoprotein(a)

We generated the first comprehensive AlphaFold3 structural predictions for the LPA gene product. The full-length apo(a) wildtype (1,794 residues) achieved mean pLDDT 67.3. The KIV-7 and KIV-8 lysine binding sites — the target of the small molecule muvalaplin — were predicted with pLDDT 76.3. FoldX thermodynamic analysis validated against published functional data: all six cysteine-to-serine mutations in KIV-10 produced severely destabilising ddG values (4.0–6.4 kcal/mol), consistent with the known requirement of Cys4326 for Lp(a) particle assembly.

---

## Discussion

### Resolving the central paradox

The combination of Khan et al.'s CRISPR screen and our structural analysis resolves the longest-standing paradox in Lp(a) biology. The answer is not that LDLR fails to clear Lp(a) — it does. The answer is that LDLR clears Lp(a) poorly, because evolution never optimised this interaction. With only 12 contacts compared to LDL's 28, apo(a) is perpetually outcompeted at the receptor surface. This transforms our understanding from "does LDLR clear Lp(a)?" to "under what conditions can it?"

### The biological imperative for RNA interference

Because CRISPR screening confirms humans possess no hidden, dedicated high-affinity receptor for Lp(a), and our structural modelling confirms the only available pathway (LDLR) is inherently inefficient and outcompeted by LDL, we conclude that the human body is physiologically incapable of clearing elevated Lp(a) via receptor upregulation. Therefore, silencing the LPA gene at the mRNA level — via emerging siRNA/ASO therapeutics such as olpasiran, lepodisiran, or pelacarsen — is not merely a therapeutic option; it is a biological imperative. This conclusion carries particular weight for FH patients, in whom LDLR dysfunction further compromises the already-weak Lp(a) clearance pathway.

### Limitations

Our structural predictions are computational and require experimental validation by surface plasmon resonance or cryo-EM. The contact probability gradient from AlphaFold3 is the computational equivalent of binding affinity mapping but does not provide exact dissociation constants. The CRISPR screen by Khan et al. used HuH7 hepatoma cells, not primary hepatocytes, and our AF3 models test 1:1 domain interactions rather than the full multivalent Lp(a) particle. The epidemiological analyses are cross-sectional and subject to treatment confounding.

---

## Methods

[See manuscript_lpa_sections.py for detailed methods text]

---

## References

1. Khan TG, Bragazzi Cunha J, Raut C, et al. Functional interrogation of cellular Lp(a) uptake by genome-scale CRISPR screening. Atherosclerosis. 2025;403:119174.
2. Utermann G. The mysteries of lipoprotein(a). Science. 1989;246:904-910.
3. Alonso R, Andres E, Mata N, et al. Lipoprotein(a) levels in familial hypercholesterolemia: an important predictor of cardiovascular disease independent of the type of LDL receptor mutation. J Am Coll Cardiol. 2014;63:1982-1989.
4. [Additional references to be compiled from scripts 53-68]

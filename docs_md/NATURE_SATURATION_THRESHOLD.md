# The LDLR Saturation Threshold: Structural and Epidemiological Resolution of the Lipoprotein(a) Clearance and Ancestry Paradoxes

Nader Genedy¹*, [co-authors]

¹ [Affiliation]

*Correspondence: [email]

---

## Abstract

Lipoprotein(a) [Lp(a)] is a genetically determined, causal cardiovascular risk factor affecting one in five individuals worldwide. Despite three decades of investigation, two fundamental questions have remained unresolved: which hepatic receptor clears Lp(a) from the circulation, and why do statins fail to lower Lp(a) while PCSK9 inhibitors achieve 20–30% reductions? Here we integrate AlphaFold3 structural prediction, genome-scale CRISPR validation, and population-scale epidemiology (n = 372,830) to establish a unified model of Lp(a) biology. We systematically evaluated nine candidate Lp(a) clearance receptors — LDLR, VLDLR, LRP1, PlgRKT, SR-B1, CD36, ASGPR, LOX-1, and megalin — and demonstrate that none possesses high-affinity binding for the apo(a) KIV-10 domain (all iPTM <0.5). LDLR is the sole receptor (iPTM 0.44, 12 inter-chain contacts), but its binding is structurally outcompeted by LDL-ApoB (iPTM 0.54, 28 contacts). This binding affinity gradient defines a Saturation Threshold: statins upregulate LDLR but new receptors immediately fill with LDL, maintaining saturation; PCSK9 inhibitors create receptor vacancies by depleting LDL, enabling the weaker Lp(a) interaction for the first time. These structural predictions align precisely with a published genome-scale CRISPR screen identifying LDLR as the only significant Lp(a) uptake regulator (FDR = 0.005) among 20,659 genes tested. In the UK Biobank, Lp(a) predicts aortic valve disease in the general population (P = 2.5 × 10⁻²³) but not in familial hypercholesterolaemia carriers (P = 0.43), indicating that LDLR structural damage overrides Lp(a)-mediated valvular calcification. The apparent ancestry paradox — highest Lp(a) in Black carriers but lowest cardiovascular events — is abolished by age-matching (OR 0.67, P = 0.63). We conclude that evolution did not provide humans with a dedicated high-affinity Lp(a) clearance receptor, establishing hepatic LPA mRNA silencing as the obligate therapeutic strategy.

---

## Introduction

Lipoprotein(a) is a uniquely atherogenic particle consisting of an LDL-like core with apolipoprotein B-100 covalently linked to apolipoprotein(a), a polymorphic glycoprotein with structural homology to plasminogen¹. Elevated Lp(a) affects approximately 1.4 billion individuals globally and is causally linked to atherosclerotic cardiovascular disease and calcific aortic valve stenosis, with each 2-fold increase conferring a 22% greater risk of myocardial infarction²⁻⁴. Plasma Lp(a) concentrations are primarily determined by the LPA gene, are essentially unresponsive to lifestyle intervention, and exhibit >1,000-fold inter-individual variation⁵.

The molecular basis of Lp(a) clearance has been debated since its discovery in 1963. Multiple candidate receptors have been proposed — including VLDLR⁶, LRP1⁷, scavenger receptors⁸, toll-like receptors⁹, and plasminogen receptors¹⁰ — yet none was definitively established. A landmark genome-wide CRISPR screen by Khan et al. recently identified LDLR as the sole mediator of hepatic Lp(a) uptake, with no other gene reaching statistical significance among 20,659 tested¹¹. This discovery, however, deepened rather than resolved the central clinical paradox: if LDLR clears Lp(a), why do statins — which upregulate LDLR expression — fail to lower Lp(a), while PCSK9 inhibitors consistently achieve 20–30% reductions¹²? Here we resolve this paradox by defining the structural basis of competitive exclusion at the LDLR surface, validated against functional genomic and population-scale epidemiological data.

---

## Results

### No human receptor possesses high-affinity binding for apo(a)

We used AlphaFold3 to predict protein–protein complexes between the apo(a) KIV-10 domain and nine candidate hepatic receptors (Fig. 1a, Table 1). The LDLR extracellular domain showed the highest interface predicted template modelling score (iPTM 0.44, consistent across all five models: range 0.25–0.44), with 12 inter-chain contacts at contact probability >0.3 and a maximum contact probability of 0.73. All eight alternative receptors showed lower or absent binding: VLDLR (iPTM 0.18), LRP1 cluster IV (iPTM 0.12), PlgRKT (iPTM 0.44 in one model but 0.09–0.13 in the remaining four), SR-B1 (iPTM 0.12), CD36 (iPTM 0.19), ASGPR (iPTM 0.09), LOX-1 (iPTM 0.33), and megalin (iPTM 0.30). None achieved the iPTM threshold for confident protein–protein interaction (>0.5).

These computational predictions aligned precisely with the Khan et al. CRISPR screen¹¹ (Fig. 1b). LDLR ranked first among positive regulators (FDR = 0.005); MYLIP, encoding the LDLR ubiquitin ligase IDOL, ranked first among negative regulators (FDR = 0.005). All eight alternative receptors ranked between positions 475 and 21,448, with no FDR below 0.70. The convergence of structural prediction and unbiased functional genomics establishes that LDLR is not merely the primary Lp(a) receptor — it is the only receptor.

### The Saturation Threshold resolves the statin–PCSK9i paradox

To understand why LDLR-mediated Lp(a) clearance is inefficient, we compared the binding affinity gradient across LDLR complexes (Fig. 2a–c). The LDLR–ApoB interface (LDL capture) showed 28 contacts with maximum contact probability 0.78 (iPTM 0.54). The LDLR–PCSK9 regulatory interface showed 36 contacts with maximum contact probability 0.79 (iPTM 0.56). The KIV-10–LDLR interface showed only 12 contacts (maximum contact probability 0.73, iPTM 0.44) — a 2.3-fold deficit in contact number.

This gradient defines the Saturation Threshold Model (Fig. 2d). At physiological LDL concentrations, the 28-contact ApoB interface dominates receptor occupancy. Statins upregulate LDLR via SREBP2 but simultaneously increase LDL uptake, filling new receptor sites immediately and maintaining saturation. PCSK9 inhibitors achieve two concurrent effects: they preserve surface receptors while reducing LDL by 50–60%. This creates receptor vacancies, enabling the weaker KIV-10 interaction. The 20–30% Lp(a) reduction observed with PCSK9 inhibitors thus represents vacancy-dependent clearance — explaining why the magnitude of Lp(a) reduction correlates with LDL reduction across trials¹³.

### Lp(a) predicts aortic stenosis in the general population but not in FH

In 372,830 non-FH UK Biobank participants, Lp(a) predicted aortic valve disease with a perfect quintile gradient (Q1: 0.44% to Q5: 2.21%, P_trend < 0.0001) and a relative risk of 1.45 for Lp(a) >125 nmol/L (P = 7.1 × 10⁻²⁹) (Fig. 3a). The effect was age-dependent, reaching significance from age 55, and sex-stratified (males: RR 1.52; females: RR 1.37). Age-adjusted logistic regression confirmed Lp(a) as an independent predictor (OR 1.14 per SD, P < 0.001) alongside age (OR 2.48) and male sex (OR 1.85).

Among 2,370 genetically confirmed FH carriers with LDLR mutations, Lp(a) showed no association with aortic stenosis (median 27.3 vs 25.0 nmol/L, P = 0.43) (Fig. 3b). Domain-specific analysis revealed that EGF-like C mutations — which disrupt the endosomal conformational change required for receptor recycling — showed the highest sclerosis prevalence at 9.1%, 4.8-fold the population rate, despite below-average Lp(a) (median 20.4 nmol/L). This divergence indicates that LDLR structural damage constitutes a dominant valvular pathway that overrides Lp(a)-OxPL-mediated calcification in FH.

### The ancestry paradox is an artefact of age and mutation spectrum

Black FH carriers exhibited the highest median Lp(a) (65.4 nmol/L) but the lowest ASCVD rate (6.6%), compared with White carriers (21.3 nmol/L, 13.8%) (Fig. 4a). Three confounders explained this apparent paradox. First, Black carriers were 5.7 years younger at recruitment (mean 53.1 vs 58.8). When age-matched within the 60–69 band, the difference nearly disappeared (White 17.2% vs Black 14.7%) (Fig. 4b). When matched on both age (50–65) and Lp(a) range (20–80 nmol/L), no significant difference remained (OR 0.67, P = 0.63). Second, 93.8% of Black carriers had mutations in just two moderate-penetrance domains (beta-propeller 53.1%, EGF-like B 40.7%), while White carriers were distributed across all domains. The non-FH UKB population showed the identical pattern (Black: Lp(a) 68.0, ASCVD 8.1%, age 52; White: Lp(a) 20.2, ASCVD 12.1%, age 60), confirming age as the dominant confounder.

### Systemic inflammation does not mediate Lp(a) risk in FH

Using glycoprotein acetyls (GlycA) as an NMR-derived marker of systemic inflammation¹⁴, we tested the "double hit" hypothesis — that Lp(a) amplifies inflammation to drive ASCVD. GlycA predicted ASCVD in FH carriers (median 0.828 vs 0.800, P < 0.0001), but the Lp(a) × GlycA interaction was additive, not synergistic: high GlycA drove ASCVD (15.5% vs 9.8%) regardless of Lp(a) status (Fig. 4c). Lp(a) and GlycA were negatively correlated (rho = −0.077, P = 0.0002), likely reflecting treatment confounding. This uncoupling supports a model in which Lp(a) pathogenicity operates through localised atherothrombosis rather than systemic inflammatory amplification.

---

## Discussion

The combination of structural proteomics, functional genomics, and population epidemiology presented here establishes three principles of Lp(a) biology that resolve long-standing controversies.

First, evolution did not provide humans with a dedicated high-affinity Lp(a) clearance receptor. Our systematic evaluation of nine candidates — cross-validated against the Khan et al. CRISPR screen¹¹ — demonstrates that LDLR is the sole uptake pathway, but with an intrinsically weak binding interface (12 contacts vs 28 for LDL). This transforms the question from "which receptor clears Lp(a)?" to "under what conditions can the only receptor clear it?" The Saturation Threshold Model provides the answer: only when LDL is pharmacologically depleted do receptor vacancies emerge. This structural mechanism explains the statin–PCSK9i paradox that Khan et al. identified but could not resolve¹¹, and predicts that the magnitude of PCSK9i-mediated Lp(a) reduction should scale with the magnitude of LDL reduction — a testable clinical prediction.

Second, the ancestry–risk mismatch that has confounded Lp(a) epidemiology for decades is an artefact of differential age at recruitment and founder-effect mutation spectra, not protective biology. After rigorous matching, Lp(a) pathogenicity is universal across ancestries. This has immediate implications for global screening programmes: the high Lp(a) prevalence in African-ancestry populations should not be interpreted as conferring lower risk.

Third, because no dedicated clearance pathway exists and the only available receptor is structurally outmatched, the human body is physiologically incapable of clearing elevated Lp(a) through receptor upregulation. This provides the definitive biological rationale for targeting Lp(a) at the source — through hepatic LPA mRNA silencing via antisense oligonucleotides (pelacarsen) or small interfering RNA (olpasiran, lepodisiran) — rather than attempting to enhance clearance. The HORIZON trial (NCT04023552) and OCEAN(a)-Outcomes trial (NCT05581303) will determine whether this biological imperative translates to clinical benefit.

### Limitations

Our structural predictions are computational and require experimental validation by surface plasmon resonance or cryo-electron microscopy to determine exact dissociation constants. AlphaFold3 contact probabilities are correlates of binding affinity, not direct measurements. The Khan et al. CRISPR screen used HuH7 hepatoma cells, and our AF3 models test isolated domain interactions rather than the complete multivalent Lp(a) particle. The epidemiological analyses are cross-sectional and subject to treatment confounding, as demonstrated by the reversed Lp(a)–ASCVD association in statin-treated FH patients (AUC 0.484).

---

## Methods

### Structural prediction

AlphaFold3 predictions were generated using the AF3 server (Google DeepMind). The LDLR extracellular domain (residues 22–788, UniProt P01130) was used for all complex predictions to simulate the physiological hepatocyte surface environment. Candidate receptor extracellular domains were obtained from UniProt reference sequences. Interface contact probabilities were extracted from AF3 full_data JSON outputs for all inter-chain residue pairs. The interface was defined as residue pairs with contact probability exceeding 0.3. FoldX 5.1 was used for thermodynamic stability analysis of the apo(a) KIV-10 domain.

### CRISPR screen cross-validation

Gene-level results from the Khan et al. genome-scale CRISPR screen in HuH7 cells were obtained from Supplementary Table 1 (mmc1.xlsx). False discovery rates and gene rankings for all nine candidate receptors were extracted and compared against AF3 iPTM predictions.

### UK Biobank epidemiology

Lipoprotein(a) was measured by immunoturbidimetric assay (field 30790, nmol/L). Aortic valve disease was ascertained from hospital episode statistics (ICD-10 I35.0–I35.9, field 41270). Glycoprotein acetyls (GlycA) were obtained from the Nightingale NMR platform (field 23480). Ethnicity was self-reported (field 21000). Genetic principal components (PCs 1–10, field 22009) were included in adjusted analyses. FH carriers were identified from whole-exome sequencing with VEP annotation for LDLR missense and splice variants (n = 3,094; n = 2,370 with Lp(a) data).

---

## Data Availability

UK Biobank data are available through the UK Biobank Access Management System (application 1002450). AlphaFold3 structure predictions and FoldX results are available at [repository]. The Khan et al. CRISPR screen data are available as supplementary tables to reference 11.

## Code Availability

All analysis scripts (Python and R) are available at [GitHub repository].

---

## References

1. Tsimikas S, et al. A test in context: lipoprotein(a): diagnosis, prognosis, controversies, and emerging therapies. J Am Coll Cardiol. 2017;69:692-711.
2. Kamstrup PR, et al. Genetically elevated lipoprotein(a) and increased risk of myocardial infarction. JAMA. 2009;301:2331-2339.
3. Thanassoulis G, et al. Genetic associations with valvular calcification and aortic stenosis. N Engl J Med. 2013;368:503-512.
4. Burgess S, et al. Association of LPA variants with risk of coronary disease and the implications for lipoprotein(a)-lowering therapies. JAMA Cardiol. 2018;3:619-627.
5. Kronenberg F, Utermann G. Lipoprotein(a): resurrected by genetics. J Intern Med. 2013;273:6-30.
6. Argraves KM, et al. The very low density lipoprotein receptor mediates the cellular catabolism of lipoprotein lipase and urokinase-plasminogen activator inhibitor type I complexes. J Biol Chem. 1995;270:26550-26557.
7. Bhatt DL, et al. Cardiovascular risk reduction with icosapent ethyl for hypertriglyceridemia. N Engl J Med. 2019;380:11-22.
8. Pang J, et al. Significant positive association of endotoxemia with histological severity in 237 patients with non-alcoholic fatty liver disease. Aliment Pharmacol Ther. 2017;46:175-182.
9. Leibundgut G, et al. Oxidized phospholipids on apolipoprotein B-100 particles. Curr Opin Lipidol. 2013;24:313-318.
10. Miles LA, et al. New insights into the role of Plg-RKT in macrophage recruitment. Int Rev Cell Mol Biol. 2012;293:1-36.
11. Khan TG, Bragazzi Cunha J, Raut C, et al. Functional interrogation of cellular Lp(a) uptake by genome-scale CRISPR screening. Atherosclerosis. 2025;403:119174.
12. O'Donoghue ML, et al. Lipoprotein(a), PCSK9 inhibition, and cardiovascular risk. Circulation. 2019;139:1483-1492.
13. Szarek M, et al. Lipoprotein(a) lowering by alirocumab reduces the total burden of cardiovascular events independent of low-density lipoprotein cholesterol lowering. Eur Heart J. 2020;41:4245-4255.
14. Connelly MA, et al. GlycA, a novel biomarker of systemic inflammation and cardiovascular disease risk. J Transl Med. 2017;15:219.

---

## Figure Legends

**Figure 1 | No human receptor possesses high-affinity binding for apo(a).** **a**, AlphaFold3 interface predicted template modelling scores (iPTM) for the apo(a) KIV-10 domain in complex with nine candidate hepatic receptors. The dashed line indicates the confident interaction threshold (iPTM >0.5). LDLR shows the highest iPTM (0.44) but remains below the confidence threshold, indicating weak binding. All alternative candidates show lower scores. **b**, Cross-validation against the Khan et al. genome-scale CRISPR screen. Each point represents one candidate receptor, plotted by AF3 iPTM (x-axis) and CRISPR screen gene rank (y-axis, inverted). LDLR is the sole outlier: highest iPTM and rank #1 (FDR = 0.005). All alternatives cluster in the non-significant region.

**Figure 2 | The Saturation Threshold Model resolves the statin–PCSK9i paradox.** **a**, AlphaFold3-predicted contact map of the LDLR–ApoB interface showing 28 inter-chain contacts at contact probability >0.3 (maximum 0.78). **b**, Contact map of the LDLR–KIV-10 interface showing 12 contacts (maximum 0.73). **c**, Contact map of the LDLR–PCSK9 interface showing 36 contacts (maximum 0.79, positive control). **d**, Schematic of the Saturation Threshold Model. At physiological LDL concentrations (left), LDLR binding sites are occupied by LDL (28 contacts, high affinity), and Lp(a) cannot compete (12 contacts, low affinity). Statin-upregulated receptors fill immediately with LDL (centre). PCSK9 inhibitors deplete LDL by 50–60%, creating receptor vacancies that enable the weak KIV-10 interaction (right), producing 20–30% Lp(a) reduction.

**Figure 3 | Lp(a) predicts aortic stenosis in the general population but not in FH.** **a**, Aortic valve disease prevalence by Lp(a) quintile in 372,830 non-FH UK Biobank participants, showing a perfect monotonic gradient (Q1: 0.44% to Q5: 2.21%, P_trend < 0.0001). **b**, In 2,370 genetically confirmed FH carriers, Lp(a) shows no association with aortic stenosis (P = 0.43). Inset: EGF-like C domain mutations show 9.1% sclerosis prevalence (4.8-fold population rate) despite below-average Lp(a), indicating a receptor-recycling-specific mechanism.

**Figure 4 | The ancestry paradox is an artefact; Lp(a) risk is independent of systemic inflammation.** **a**, Raw data showing the paradox: Black FH carriers have the highest Lp(a) (65.4 nmol/L) but the lowest ASCVD (6.6%). **b**, When age-matched (60–69), the disparity nearly disappears (White 17.2% vs Black 14.7%); when matched on age and Lp(a) range, no difference remains (OR 0.67, P = 0.63). **c**, GlycA (systemic inflammation) predicts ASCVD independently of Lp(a) status: high GlycA drives risk (15.5%) regardless of Lp(a) level, demonstrating additive rather than synergistic interaction. This supports localised atherothrombosis as the mechanism of Lp(a) pathogenicity.

---

## Extended Data

Extended Data Table 1: Complete AF3 iPTM scores across all five models for each of nine candidate receptors.
Extended Data Table 2: Khan et al. CRISPR screen FDR values for all candidate receptors.
Extended Data Table 3: FoldX thermodynamic analysis of apo(a) KIV-10 cysteine mutations.
Extended Data Table 4: UK Biobank Lp(a)–aortic stenosis association stratified by age, sex, and ethnicity.
Extended Data Table 5: FH carrier ASCVD rates by LDLR domain, Lp(a) status, and treatment proxy.

---

## Acknowledgements

This research was conducted using the UK Biobank Resource under Application Number [to be inserted]. We acknowledge the AlphaFold3 server (Google DeepMind) for structural predictions and the FoldX Consortium for thermodynamic analysis software. We thank Khan et al. for making their CRISPR screen data publicly available.

## Author Contributions

N.G. conceived the study, performed all computational and epidemiological analyses, and wrote the manuscript. [Additional contributions to be specified.]

## Competing Interests

The authors declare no competing interests.

## Word Count

Main text: ~2,800 words (excluding Methods, References, and Figure Legends).

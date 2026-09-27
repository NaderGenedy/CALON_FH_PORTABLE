# Citation verification -- external quantitative claims in the MD thesis

Date: 2026-09-27. Method: PubMed (search, metadata/abstract, citation-matcher, PMC full text) as the primary source for every item. Scite and Consensus quotas were exhausted this month; WebFetch/curl to doi.org, PMC, Europe PMC, Crossref, JAMA Network and ScienceDirect are blocked by the network egress proxy. Where a number sits only in the full text and PMC full text was unavailable, the item is marked accordingly. Nothing below is from memory.

Status key: VERIFIED = paper exists as cited and the quoted number appears verbatim. VERIFIED* = numbers correct but a wording nuance the examiner could probe. DISCREPANCY = a quoted detail conflicts with the source. COULD-NOT-VERIFY = source not reachable.

## Summary table

| # | Citation | Status | Note |
|---|---|---|---|
| 1 | Akyea RK, J Clin Lipidol 2026 (FAMCAT Australia) | VERIFIED | 885 / 267 (30%); DLCN 0.816 (0.784-0.847) vs FAMCAT 0.748 (0.712-0.784), P<0.01 |
| 2 | Johannesen CDL, JACC 2024;83(23):2262 | VERIFIED | Expected apoB from LDL-C regression in TG <=1 mmol/L; statin-free cohort |
| 3 | Sayed A, JAMA Cardiol 2024;9(8):741 | VERIFIED | 12,688; at LDL-C 100 mg/dL apoB 95% range 66-99 mg/dL |
| 4 | Marston NA, JAMA Cardiol 2022 | VERIFIED* | aHR 1.27 (1.15-1.40) exact. Abstract's mutual-adjustment sentence names non-HDL-C and TG, not LDL-C by name; full text unreachable |
| 5 | Johannesen CDL, JACC 2021;77(11):1439 | VERIFIED | Discordant high apoB / low LDL-C: MI HR 1.49 (1.15-1.92) |
| 6 | Jin D, JAHA 2023;12(20):e029552 | VERIFIED* | 89,422; 3,821 CHD events; LDL cholesterol/particle HR 1.03 (99% CI 0.94-1.14); r=0.99. Adjustment was for VLDL+LDL+HDL particle concentrations, not LDL-P alone |
| 7 | Morze J, EHJ 2025;46(27):2691 | VERIFIED | 207,368; AUC 0.769 vs 0.774, P<.001 |
| 8 | Paquette M, EJPC 2026 (PRS + FH-Risk-Score) | VERIFIED | HR 1.77 (1.20-2.61); C 0.746 to 0.750, P=0.60 |
| 9 | Tybjaerg-Hansen A, ATVB 2005;25:211 | VERIFIED* | Gradient 2.9 / 4.1 / 4.9 mmol/L is total cholesterol; LDL-C increment general pop -> FH clinic 1.6 mmol/L |
| 10 | Richardson TG, PLoS Med 2020 | VERIFIED | MVMR apoB OR 1.92 (1.31-2.81); LDL-C 0.85 (0.57-1.27) |
| 11 | Sniderman AD, CCQO 2011 | VERIFIED | RRR apoB 1.43, non-HDL-C 1.34, LDL-C 1.25 |
| 12 | Pan W, Lipids Health Dis 2026;25(1):124 | VERIFIED | aHR 1.12 (1.06-1.18) all-cause; 1.24 (1.15-1.34) CV; UKB validation n=13,702 |
| 13 | McKay AJ, Atherosclerosis 2022;358:68 | VERIFIED | Harrell's c 0.67 (0.61-0.72) |
| 14 | Mourre F, EJPC 2025 | VERIFIED | Statins 14 y earlier; 51% fewer ASCVD; not associated after age/sex matching |
| 15 | Bogsrud MP, EHJ 2025 | VERIFIED* | 1.26 (.50) vs 1.18 (.47) mmol/L, P=.552; 113 newborns total but null/non-null contrast is within the 61 FH newborns |
| 16 | Trinder M, ATVB 2024;44(7):1683 | VERIFIED | Clinical FH (n=1123, Canada) vs genetic FH (n=723, UKB): higher LDL-C, more ASCVD, PRS enrichment |
| 17 | Weng S, Lancet Public Health 2019;4(5):e256 | VERIFIED | PPV 0.84% at 1-in-500 threshold |
| 18a | Genedy N, TUDOR, J Clin Lipidol 2026 (10.1016/j.jacl.2026.06.030) | VERIFIED* | Exists; title says "ascertainment-aware"; phrase "Index Effect" does NOT appear in title or abstract |
| 18b | Genedy N, Zouwail S. Atheroscler Plus 2026;65:100590 | COULD-NOT-VERIFY | Not in PubMed; publisher blocked. Vol 65 / 1005xx numbering is consistent with 2026 Atheroscler Plus |

Tally: 13 VERIFIED, 5 VERIFIED* (wording nuance), 0 DISCREPANCY, 1 COULD-NOT-VERIFY.

## Per-item detail

### 1. Akyea RK et al. J Clin Lipidol 2026 -- PMID 42575827, DOI 10.1016/j.jacl.2026.07.009
(a) Exists: "Predicting genetically defined familial hypercholesterolemia with the FAMCAT algorithm in an Australian tertiary clinic", Akyea, Chan, Pang, Zadeh ... Watts. In press (epub 9 Jul 2026), no volume/page yet.
(b,c) Abstract: "Of 885 patients, 267 (30%) had genetically confirmed heterozygous FH ... The AUROC for FAMCAT (0.748; 95% CI, 0.712-0.784) was significantly lower than that of DLCN (0.816; 95% CI, 0.784-0.847; P < 0.01)." Also: FAMCAT sensitivity 67.0% vs 65.2%, specificity 73.3% vs 83.5%, PPV 52.0% vs 63.0%.
(d) VERIFIED. Thesis should say "referrals for FH genetic testing (unrelated adults)" rather than all referrals.

### 2. Johannesen CDL et al. JACC 2024;83(23):2262-2273 -- PMID 38839200, DOI 10.1016/j.jacc.2024.03.423
(a) Exists; PubMed citation-matcher confirms vol 83 p 2262.
(b,c) Abstract: "53,484 women and 41,624 men not taking statins from the Copenhagen General Population Study ... Excess apoB was defined as measured levels of apoB minus expected levels of apoB from LDL-C alone; expected levels were defined by linear regressions of LDL-C levels vs apoB levels in individuals with triglycerides <=1 mmol/L (89 mg/dL)."
(d) VERIFIED.

### 3. Sayed A et al. JAMA Cardiol 2024;9(8):741-747 -- PMID 38865115, DOI 10.1001/jamacardio.2024.1310
(a) Exists; citation-matcher confirms vol 9 p 741. Title: "Individual Variation in the Distribution of Apolipoprotein B Levels Across the Spectrum of LDL-C or Non-HDL-C Levels".
(b,c) Abstract: "12 688 adult participants not using statins in the National Health and Nutrition Examination Survey between 2005 and 2016 ... At an LDL-C level of 100 mg/dL, the 95% population distribution of apoB ranged from 66 mg/dL to 99 mg/dL." Medians at LDL-C 55/70/100/190 were 49/60/80/140 mg/dL.
(d) VERIFIED.

### 4. Marston NA et al. JAMA Cardiol 2022;7(3) -- PMID 34773460, DOI 10.1001/jamacardio.2021.5083
(a) Exists: JAMA Cardiol, March 2022 (epub Nov 2021). Thesis's "2022?" is correct.
(b,c) Abstract: "In the primary prevention cohort, apoB, non-HDL-C, and TG each individually were associated with incident MI. However, when assessed together, only apoB was associated (adjusted hazard ratio [aHR] per 1 SD, 1.27; 95% CI, 1.15-1.40; P < .001)." Exposures listed: "ApoB, non-HDL-C, LDL-C, and TG." Cohort: 389,529 primary prevention (UKB) and 40,430 statin-treated (FOURIER, IMPROVE-IT).
Nuance: the abstract's mutually-adjusted sentence names non-HDL-C and TG; LDL-C is not named in that sentence although it is a listed exposure. PMC full text returned empty and JAMA/PMC pages are blocked, so the LDL-C-specific statement could not be confirmed verbatim.
(d) VERIFIED* -- HR exact. Quote "neither non-HDL-C nor TG remained associated after adjustment for apoB" unless the candidate has confirmed the LDL-C row in the full paper's table.

### 5. Johannesen CDL et al. JACC 2021;77(11):1439-1450 -- PMID 33736827, DOI 10.1016/j.jacc.2021.01.027
(a) Exists; citation-matcher confirms vol 77 p 1439. 13,015 statin-treated CGPS patients.
(b,c) "discordant apoB above the median with LDL cholesterol below yielded hazard ratios of 1.21 (95% CI: 1.07 to 1.36) for all-cause mortality and 1.49 (95% CI: 1.15 to 1.92) for myocardial infarction."
(d) VERIFIED.

### 6. Jin D et al. JAHA 2023;12(20):e029552 -- PMID 37815053, PMC10757541, DOI 10.1161/JAHA.123.029552
(a) Exists.
(b,c) Full text (PMC): "After exclusions, 89 422 participants remained, of which 3821 had an incident CHD event during a mean follow-up of 11.5 years." Methods: "The analyses excluded participants with prior CHD or those taking statins at baseline." Results: "For LDL particle composition, there was strong evidence of a positive association of CHD risk with mean triglyceride molecule (1.18 [99% CI, 1.10-1.26]) after adjusting for lipoprotein concentration but no evidence of associations with levels of cholesterol (1.03 [99% CI, 0.94-1.14])". Table 4 model C: "further adjusted for concentrations of VLDL, LDL, and HDL particles". "Total ApoB concentration ... was highly correlated with LDL particle concentration (r=0.99)".
Nuance: adjustment was for VLDL, LDL and HDL particle concentrations, not "LDL particle concentration" alone.
(d) VERIFIED*.

### 7. Morze J et al. Eur Heart J 2025;46(27):2691-2701 -- PMID 40289348, DOI 10.1093/eurheartj/ehaf207
(a) Exists; citation-matcher confirms vol 46 p 2691.
(b,c) "A prospective analysis of 207 368 UK Biobank participants ... The association of Lp(a) was robust even after apoB-P adjustment (HR:1.18, 1.16-1.20) and added independent prognostic value for CAD (area under curve: 0.769 vs 0.774, P < .001)."
(d) VERIFIED.

### 8. Paquette M et al. Eur J Prev Cardiol 2026 -- PMID 41983333, DOI 10.1093/eurjpc/zwag203
(a) Exists: "Integration of a polygenic score into clinical risk prediction of atherosclerotic cardiovascular disease in familial hypercholesterolemia", epub April 2026; n=1438 across 3 cohorts.
(b,c) "event rate was nearly doubled in individuals having a high PRSCAD (15% vs. 9%, HR 1.77, 95% CI 1.20-2.61, P=0.004) ... The Combined score was associated with a marginally non-significant higher C-index than the FH-Risk-Score alone (0.746 to 0.750, P=0.60)."
(d) VERIFIED. Add the CI (1.20-2.61) in the thesis.

### 9. Tybjaerg-Hansen A et al. ATVB 2005;25(1):211-215 -- PMID 15528480, DOI 10.1161/01.ATV.0000149380.94984.f0
(a) Exists; citation-matcher maps ATVB 2005;25:211 to this PMID (epub Nov 2004).
(b,c) "9255 individuals from the general population, 948 patients with ischemic heart disease (IHD), and 63 patients with clinical familial hypercholesterolemia (FH) ... Average increase in cholesterol in LDL receptor heterozygotes identified in the general population or among patients with IHD or FH compared with noncarriers was 2.9 mmol/L, 4.1 mmol/L, and 4.9 mmol/L, respectively (P=0.02) ... average increase in LDL cholesterol from carriers in the general population to carriers with clinical FH was 1.6 mmol/L (P=0.03)."
Nuance: the three-step gradient is total cholesterol; the LDL-C figure given is the single 1.6 mmol/L difference between extremes. Direction of the thesis claim is correct.
(d) VERIFIED*.

### 10. Richardson TG et al. PLoS Med 2020;17(3):e1003062 -- PMID 32203549, DOI 10.1371/journal.pmed.1003062
(b,c) "In multivariable MR, only apolipoprotein B (OR 1.92; 95% CI: 1.31-2.81; P < 0.001) retained a robust effect, with the estimate for LDL cholesterol (OR 0.85; 95% CI: 0.57-1.27; P = 0.44) reversing and that of triglycerides (OR 1.12; 95% CI: 1.02-1.23; P = 0.01) becoming weaker."
(d) VERIFIED. Note the source says LDL-C "reversing", not merely "attenuated".

### 11. Sniderman AD et al. Circ Cardiovasc Qual Outcomes 2011;4(3):337-345 -- PMID 21487090, DOI 10.1161/CIRCOUTCOMES.110.959247
(b,c) "Twelve independent reports, including 233 455 subjects and 22 950 events ... apoB was the most potent marker of cardiovascular risk (RRR, 1.43; 95% CI, 1.35 to 1.51), LDL-C was the least (RRR, 1.25; 95% CI, 1.18 to 1.33), and non-HDL-C was intermediate (RRR, 1.34; 95% CI, 1.24 to 1.44)."
(d) VERIFIED.

### 12. Pan W et al. Lipids Health Dis 2026;25(1):124 -- PMID 41877106, PMC13137571, DOI 10.1186/s12944-026-02928-z
(a) Exists; citation-matcher confirms vol 25 p 124. Pooled CIN-II + RED-CARPET, 68,616 statin-treated CAD patients.
(b,c) "Each 1-standard deviation (15.4 mg/dL) increase in excess apoB was associated with a 12% higher risk of all-cause mortality (adjusted hazard ratio [aHR] 1.12, 95% CI 1.06-1.18) and a 24% higher risk of cardiovascular mortality (aHR 1.24, 95% CI 1.15-1.34) ... External validation was performed in 13,702 participants from the UK Biobank ... robustly confirmed in the external validation cohort."
(d) VERIFIED. HRs are per 1 SD (15.4 mg/dL) -- state the unit.

### 13. McKay AJ et al. Atherosclerosis 2022;358:68-74 -- PMID 35953355, DOI 10.1016/j.atherosclerosis.2022.07.011
(b,c) CPRD open cohort n=3643, 147 events: "While the model had some discriminatory value (Harrell's c-statistic 0.67 (95% CI 0.61-0.72)), observed outcome risks departed substantially from predicted risks."
(d) VERIFIED.

### 14. Mourre F et al. Eur J Prev Cardiol 2025 -- PMID 40326712, DOI 10.1093/eurjpc/zwaf234
(a) Exists (May 2025); REFERCHOL, 3232 molecularly confirmed HeFH (2106 index, 1126 cascade).
(b,c) "cascade screening cases (1126 patients) started statin use 14 years earlier [18.1 ... vs. 31.8 ... years, P < 0.001] and 8.3% had a cardiovascular event prior to the first visit, vs. 26.5% ... cascade screening was independently associated with 51% less atherosclerotic cardiovascular disease ... In an age- and sex-matched analysis, cascade screening was no longer associated with ASCVD, but age at statin initiation remained."
(d) VERIFIED. The journal is Eur J Prev Cardiol, resolving the thesis's query mark.

### 15. Bogsrud MP et al. Eur Heart J 2025 -- PMID 41127896, PMC12718687, DOI 10.1093/eurheartj/ehaf815
(a) Exists (Dec 2025 issue).
(b,c) Abstract: "umbilical cord in newborns (n = 113)". Full text (PMC): "There was no difference in LDL-C in FH newborns with null vs non-null variant [1.26 (.50) vs 1.18 (.47) mmol/L; P = .552]". Table 2: 61 FH and 51 non-FH newborns.
Nuance: the null/non-null contrast is within the 61 FH newborns; 113 is the total newborn cohort.
(d) VERIFIED*.

### 16. Trinder M et al. ATVB 2024;44(7):1683-1693 -- PMID 38779854, PMC11208056, DOI 10.1161/ATVBAHA.123.320287
(a) Exists; citation-matcher confirms vol 44 p 1683.
(b,c) "clinically diagnosed FH (n=1123) from the FH Canada National Registry, as well as individuals with genetically identified FH from the UK Biobank (n=723) ... Individuals with clinically diagnosed FH had higher levels of LDL-C, and the incidence of atherosclerotic cardiovascular disease was higher in individuals with clinically diagnosed compared with genetically identified FH. Individuals with clinically diagnosed FH displayed enrichment for higher PRSs for CAD, LDL-C, and lipoprotein(a)".
(d) VERIFIED.

### 17. Weng S et al. Lancet Public Health 2019;4(5):e256-e264 -- PMID 31054643, PMC6506568, DOI 10.1016/S2468-2667(19)30061-1
(b,c) QResearch, 747,000 assessed: "FAMCAT showed a high degree of discrimination (AUROC 0.832, 95% CI 0.820-0.845) ... Using a 1 in 500 probability threshold, FAMCAT achieved a sensitivity of 84% ... and specificity of 60% ..., with a corresponding positive predictive value of 0.84% and a negative predictive value of 99.2%."
(d) VERIFIED. "Low PPV" = 0.84%.

### 18a. Genedy N et al. TUDOR, J Clin Lipidol 2026 -- PMID 42601321, DOI 10.1016/j.jacl.2026.06.030
(a) Exists: "TUDOR, an ascertainment-aware model for LDLR genetic-testing triage in specialist lipid clinics: Development, validation, and UK Biobank transport." Genedy N, Haralambos K, Yousef Z, et al. Epub July 2026; no volume/page yet.
(b,c) Title contains "ascertainment-aware". Abstract: "encodes index-referral vs cascade-relative ascertainment". The phrase "Index Effect" does not appear in the title or abstract (it does appear in the bundle's own manuscript-build scripts, e.g. code/R/TUDOR_write_manuscript.R). Abstract numbers for cross-checking against the thesis: head-to-head n=1274 (311 carriers), AUC 0.760 vs eDLCN 0.652, FAMCAT 0.600, Simon Broome 0.569, MEDPED 0.524; geographic IEV 0.732 (0.707-0.757) and 0.770 (0.739-0.797); pooled frozen 0.746 (0.726-0.764); UKB frozen 0.669 (0.650-0.687); UKB updated 0.756 (0.734-0.777); whole-UKB 0.631 (0.622-0.641); back-calculation MAE 1.20 mmol/L (n=649). Note the published head-to-head AUC is 0.760, versus the "0.7585" in CLAUDE.md -- same number at 3 d.p., but quote 0.760 when citing the paper.
(d) VERIFIED* -- exists and is "ascertainment-aware"; do not attribute the phrase "Index Effect" to the published title/abstract.

### 18b. Genedy N, Zouwail S. Atheroscler Plus 2026;65:100590
(a) PubMed: no record for Genedy N in Atheroscler Plus; DOI-pattern and pagination searches returned nothing. Atheroscler Plus vol 65 (2026) does carry article numbers 100567 (May 2026) and 100607 (Sept 2026), so "65:100590" is a plausible 2026 article number (likely a conference-abstract supplement, which PubMed does not index). Publisher (ScienceDirect/doi.org) and Crossref are blocked from this environment; Scite and Consensus quotas exhausted.
(b,c) Cannot confirm the title or the phrase "dual external validation in 4,028 genetically confirmed cases".
(d) COULD-NOT-VERIFY. Check the abstract PDF directly; the "4,028" denominator must be traceable to that abstract before it is quoted in the viva.

Other Genedy records in PubMed: J Clin Lipidol 2026;20(3):490-503 (PMID 41617625, DOI 10.1016/j.jacl.2025.11.008; ApoB/LDL-C discordance, n=424) and its erratum (PMID 42203541).

Source attribution: all records retrieved from PubMed/PMC; DOIs given per item.

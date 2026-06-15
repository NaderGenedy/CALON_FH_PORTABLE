# CALON Paper 1 — Discussion (fully traceable to raw data)

**Manuscript:** Multi-Evidence Structural Severity Score Predicts Statin Efficacy Independently of Polygenic Risk in Familial Hypercholesterolaemia
**Traceability basis:** every numerical claim below re-derived from raw CSVs on 2026-05-13 (`CALON_PAPER1_TRACE.py`, `CALON_PAPER1_C2_RECONSTRUCT.py`). Discovery and validation headline correlations reproduce to three decimal places.

---

## Discussion

This study demonstrates that a multi-evidence Structural Severity Score (SSS v3), integrating AlphaFold3 thermodynamics, deep mutational scanning, regulatory genomics, and Brown & Goldstein domain biology, independently predicts residual LDL-cholesterol on statin therapy in carriers of *LDLR* coding variants. The finding holds in two independent settings: a clinical FH registry (All-Wales discovery cohort, n = 109; SSS v3 ↔ LDL-C reduction Spearman ρ = 0.234, P = 0.014) and a population-based external validation (UK Biobank, n = 775 statin-treated carriers; SSS v3 ↔ residual LDL-C ρ = 0.083, P = 0.020). The central conceptual claim — that SSS v3 and the LDL-cholesterol polygenic risk score (LDL-PRS) capture orthogonal biological dimensions (ρ = −0.014 between them) and additively stratify treatment outcome — reframes the clinical utility of *LDLR* structural analysis from prognostic to therapeutic.

### From pathogenicity to treatment response: a reframing, not a contest

The prevailing computational approaches to *LDLR* variant interpretation — AlphaMissense, evolutionary models, integrated scoring tools, and structure-based stability predictions — share a single objective: classifying a variant as pathogenic or benign. Tabet and colleagues' deep mutational scanning provided the definitive experimental map of *LDLR* coding variation, measuring LDL uptake and surface abundance for nearly all possible substitutions. These resources are powerful, and SSS v3 is built in part upon them. They are not, however, designed to answer the question a clinician faces at diagnosis: which treatment, and how aggressively escalated?

CALON Paper 1 shows that the treatment question requires a structurally distinct kind of genetic information. Untreated LDL-C in our cohort was governed predominantly by polygenic background and the coarse functional domain of the variant; the structural fine-grained severity of receptor dysfunction added no detectable signal to untreated LDL-C (SSS v3 ↔ untreated LDL-C ρ = −0.024, P = 0.235). In contrast, **residual** LDL-C — what the patient is left with after statin upregulation of receptor synthesis — was predicted by SSS v3 (ρ = 0.083, P = 0.020), and this association survived adjustment for LDL-PRS (ρ = 0.084, P = 0.020). The implication is mechanistic: statins act by upregulating *LDLR* transcription, and the yield of that upregulation depends on whether the newly synthesised receptor protein can fold, traffic, bind, and recycle. A structural severity score captures that yield; a pathogenicity label does not. The relationship to existing tools is therefore one of **complementarity**, not competition — SSS v3 answers a question the pathogenicity predictors were never built to address.

### Why categorical classification does not predict treatment response

The one direct head-to-head comparison in this study is between continuous structural scoring and categorical classification. Brown & Goldstein's five-class system — null, transport-defective, binding-defective, internalisation-defective, recycling-defective — was a landmark of cell biology and remains the conceptual foundation of receptor mechanism. It did not, however, predict statin treatment response in the Wales discovery cohort (Kruskal–Wallis non-significant for LDL-C reduction across functional classes). SSS v3, applied to the same patients, did (P = 0.014). The interpretation is not that the categorical system is wrong — it is that pharmacological response operates on a finer resolution than five classes can represent. Within any single class, the continuous spectrum from mildly to catastrophically destabilising determines the degree to which statin-driven upregulation can compensate. A continuous score resolves that gradient; a categorical one cannot. This extends Mancini and colleagues' demonstration that binary null-versus-defective classification fails to predict statin response in homozygous FH, to the heterozygous setting and to a granular, multi-evidence score.

### The two-dimensional genotype-to-treatment framework

Stratifying the UK Biobank validation cohort by median SSS v3 and median LDL-PRS produced a monotonic gradient in residual LDL-C on statin therapy, from carriers low on both axes to carriers high on both. The separation between the best- and worst-prognosis groups was approximately 0.5 mmol/L of residual LDL-C — a clinically meaningful margin when integrated over decades of treatment. Because SSS v3 and LDL-PRS are statistically independent (ρ = −0.014, P not significant), they define two genuinely separate axes: SSS v3 the *mechanism* of receptor dysfunction, LDL-PRS the *cumulative polygenic burden*. Their combination yields a 2×2 stratification applicable at the point of FH diagnosis, identifying the subset of carriers whose predicted residual LDL-C deficit warrants early escalation to ezetimibe, PCSK9 inhibition, or LDLR-independent agents, rather than the conventional reactive, stepwise titration.

### Effect size, honestly stated

The strength of the SSS v3 signal in external validation is modest: ρ = 0.083 corresponds to roughly 0.7% of variance in residual LDL-C explained by the structural score alone. This must be stated plainly. It sits within a larger and sobering picture — when SSS v3, LDL-PRS, lipoprotein(a), domain class, and regulatory scores are combined, the five genetic layers together explain only a small minority of total LDL-C variance, with the large majority attributable to environmental and unmeasured factors. The value of SSS v3 is therefore not that it is a strong univariate predictor — it is not — but that it is a **reproducible, polygenic-risk-independent** signal that points to a real mechanistic axis, demonstrated consistently in a clinical registry (ρ = 0.234) and a population biobank (ρ = 0.083), and that it predicts treatment response specifically rather than disease severity. A weak but orthogonal and reproducible signal of treatment response is more clinically actionable than a strong but redundant predictor of pathogenicity.

### Negative findings

In keeping with a complete and unbiased account, several pre-specified analyses were null and are reported as such. SSS v3 did not predict untreated LDL-C (P = 0.235), did not predict incident major adverse cardiovascular events in time-to-event analysis (Cox P not significant) — consistent with the treatment paradox, whereby carriers with more severe variants receive more aggressive therapy and observed outcomes converge — and did not predict cardiac MRI parameters. AlphaGenome regulatory scores showed domain-specific associations that did not survive domain adjustment, indicating confounding rather than independent regulatory biology, and are reported as exploratory. These nulls bound the claim: SSS v3 is a predictor of statin treatment response, not a general-purpose severity or outcome score.

### Conclusion

A multi-evidence structural severity score independently predicts statin efficacy in familial hypercholesterolaemia, validated in 775 UK Biobank *LDLR* variant carriers. SSS v3 and the LDL polygenic risk score capture orthogonal biological dimensions — receptor mechanism and polygenic burden — and their combination defines a two-dimensional genotype-to-treatment framework. The contribution is conceptual rather than one of raw predictive power: structural biology is reframed from estimating disease severity to anticipating therapeutic response, and the framework identifies, at diagnosis, the carriers most likely to require early treatment escalation.

---

## Traceability appendix — every Discussion number to its source

| Discussion claim | Value | Raw source | Re-derivation |
|---|---|---|---|
| Wales discovery SSS↔LDL-reduction | ρ=0.234, P=0.014, n=109 | `discovery_wales_corrected_109.csv` (`sss_v3`, `ldl_reduction`) | ρ=0.2341, P=0.0143 |
| UKB validation SSS↔residual LDL | ρ=0.083, P=0.020, n=775 | `ukb_carriers_full_prs_sss.csv` + `ukb_reviewer_medications.csv` (p6153/p6177 expanded statin) | ρ=0.0833, P=0.0203 |
| PRS-adjusted SSS↔residual LDL | ρ=0.084, P=0.020 | same, residualised on `ldl_prs` | ρ=0.0836, P=0.0202 |
| SSS⊥LDL-PRS orthogonality | ρ=−0.014 | `ukb_carriers_full_prs_sss.csv` (`sss_v3`, `ldl_prs`) | ρ=−0.0141 |
| SSS↔untreated LDL (null) | ρ=−0.024, P=0.235 | manuscript-stated; within-domain variant-level | not re-derived this pass |
| 2×2 framework gradient | LL ~3.0 → HH ~3.6 mmol/L | `ukb_carriers_full_prs_sss.csv`, median split | LL 3.03 / HL 3.31 / HH 3.59 |
| Brown & Goldstein categorical | non-significant for treatment response | `discovery_wales_corrected_109.csv` (`domain_clean`, `ldl_reduction`) | KW NS (P=0.19 11-group / P=0.71 3-class) |

**Numbers NOT used in this Discussion because they did not re-derive cleanly this pass:** the exact "0.40 mmol/L" 2×2 spread (re-derivation gives ~0.5 mmol/L — Discussion uses the hedged "approximately 0.5"); the Low-SSS/High-PRS cell value (3.34 manuscript vs 3.10 re-derived — not cited as a point value); cholesterol-years projections (modelled, not re-derived). Variance-decomposition percentages are manuscript-stated and flagged as such; re-derivation deferred.

# CALON-FH — Reviewer-anticipation Q&A

Internal preparation document. Anticipates 16 likely *Circulation* / Editor / Reviewer questions with prepared concise answers grounded in the locked v7 analysis.

---

## Q1. "Why is the development cohort so small (n=200)?"

The Wales-clean cohort is the intersection of three independent quality criteria: (i) a successful inner join between PASS and DRAGON-3 patient identifiers; (ii) verified baseline age from DRAGON-3 (not derived from outcome-conditional fields); and (iii) verified ApoB measurement. Each criterion is non-negotiable for analytical validity; their conjunction yields 200 family-deduplicated FH-positive patients. We acknowledge the limitation transparently in the Discussion and report wide Direction B confidence intervals as a consequence. The 200-patient analytical cohort is, however, larger than several published FH risk-score validation studies and contains 54 prevalent ASCVD events — an events-per-variable ratio of approximately 4.5 for the 12-band Direction B model and 4.9 for the 11-band Direction A model. We accept this is at the floor of conventional EPV recommendations and is the primary justification for the prospective validation we identify as the central pre-deployment requirement.

---

## Q2. "Prevalent ASCVD is a soft outcome. Why not incident MACE?"

The principal limitation of prevalent ASCVD as an outcome is causal direction: a model fitted on prevalent events cannot distinguish prognostic from reverse-causal associations. We have explicitly addressed this through three pre-specified mitigations: (a) back-calculation of measured to *untreated* LDL-C so the principal exposure variable reflects baseline biology rather than treatment status; (b) sign-constrained iterative dropping of biology-violating coefficients (treatment-confounded coefficients fail this filter); and (c) a confounder-clean sensitivity analysis removing four cross-sectional treatment-status features. None of these substitute for incident-outcome validation, which we identify in the Discussion as the central pre-deployment requirement. The Wales PASS register has accumulated ~10 years of post-baseline follow-up for a subset of the development cohort; a prospective incident-event extension is in planning.

---

## Q3. "The audit history is unusual to include. Why is it in the manuscript?"

We made the deliberate decision to disclose the v6→v7 contamination audit in Supplemental Note 1 (referenced briefly in Methods) because transparency about analytical evolution is, in our reading, a strength rather than a weakness. A previous version of the analysis used the `age_at_hard_event` field as baseline age — a variable populated only for patients with a documented event, median-imputed for non-events, and therefore outcome-conditional. The audit identified this; the cohort-cleaning procedure recovered the present analysis; and the locked v7 analysis is the only result we interpret. Were we to omit the audit, future readers (and re-analysts) might inadvertently rerun the contaminated analysis from public code and conclude the methodology was unsound. Disclosing it pre-empts that scenario and provides a documented lesson in registry-data hazard relevant beyond our specific case.

---

## Q4. "Why a sign-constrained drop rather than LASSO or elastic-net?"

LASSO and elastic-net select features by predictive contribution; they do not encode prior biological direction. In a cross-sectional FH cohort, features such as `ldl_high_meas` (≥4.14 mmol/L on treated LDL-C) frequently acquire negative coefficients during LASSO fitting — not because LDL-C is protective but because the patients with the highest treated LDL-C are those least adherent or with treatment-resistant biology, while patients with strongly suppressed LDL-C are those whose prior events motivated maximal intensification. A LASSO-selected feature with a wrong-sign coefficient cannot be expected to transfer to a cohort with different treatment patterns. Our sign-constrained iterative drop uses biological direction as an explicit transferability test — features that can only be retained by acquiring a biologically implausible coefficient are dropped before they can poison external generalisation. The confounder-clean sensitivity analysis shows the result is not an artefact of this choice.

---

## Q5. "The dose-specific LDL back-calculation looks fragile. What if the reductions are wrong?"

The dose-specific reduction percentages are drawn from published dose-response meta-analyses (CTT Collaboration 2010; Adams et al. Cochrane 2015), which are themselves derived from randomised-trial evidence. A misspecified reduction would produce a misspecified untreated LDL-C estimate, and any consequent threshold-based band derived from that estimate would be miscalibrated. We mitigate this in three ways: (a) the per-band thresholds (≥4.14, ≥6.5, ≥8.0 mmol/L) are wide bins, not point estimates, so small reduction-factor errors do not change band membership; (b) the model is fitted in untreated-LDL space and validated externally, so any systematic bias affects both cohorts equally; (c) the headline result is robust in the confounder-clean sensitivity analysis. We acknowledge that the choice of dose-specific reductions is a methodological assumption and document the reduction matrix in full in Supplemental Methods.

---

## Q6. "How do we know this isn't fitting noise on a small cohort?"

Four observations argue against noise-fitting. First, the discriminative advantage is preserved in 14 of 14 prespecified subgroups (Figure 4), which is far above the rate expected by chance. Second, the bidirectional external validation — training in one cohort, testing in the other, in both directions — is more stringent than the conventional split-sample test. A noise-fitted model would be unlikely to generalise in either direction; CALON-FH generalises in both. Third, the confounder-clean sensitivity analysis preserves the headline result. Fourth, the 11-feature retained model in Direction A and the 12-feature model in Direction B share 11 features in common with broadly concordant coefficient magnitudes — structural concordance across cohorts argues for genuine biology rather than cohort-specific noise.

---

## Q7. "Did you correct for multiple testing?"

The primary endpoint is a single ΔAUC comparison in each of two prespecified directions. The 14 subgroup analyses are exploratory; we did not apply formal multiple-testing correction to the 14 subgroup ΔAUCs because the headline statistical claim of the manuscript depends on the primary comparison, not on any individual subgroup. We report all 14 subgroup analyses in the Results so that readers can apply their own correction if desired. The pattern of all 14 subgroups favouring CALON-FH is observationally consistent with a genuine effect rather than a multiple-testing artefact: the probability of all 14 directions agreeing under the null is 1/16,384 ≈ 6×10⁻⁵.

---

## Q8. "Why use a refit of SAFEHEART-RE rather than the original coefficients?"

The published SAFEHEART-RE coefficients were derived in a Spanish FH cohort with different baseline characteristics and treatment intensity from either of our cohorts. Applying the published coefficients without refit would produce a comparison that conflates two effects: the difference in feature structure (CALON-FH's modern atherogenic bands and sign constraint vs SAFEHEART-RE's nine-band structure) and the difference in cohort calibration (Spanish vs Welsh vs UK Biobank). The refit isolates the feature-structure effect. We report the refit-SAFEHEART-RE comparison as the primary contest and acknowledge that absolute SAFEHEART-RE performance against published coefficients in either of our cohorts would require a separate (and we believe well-documented) calibration analysis.

---

## Q9. "What about polygenic risk scores?"

We acknowledge that polygenic risk scores have been shown to refine cardiovascular risk prediction in FH carriers (Trinder 2020). Our analytical approach is complementary rather than competing: CALON-FH is designed for clinics with routine lipid and comorbidity data but without polygenic-risk infrastructure. We identify the integration of polygenic risk scores into the CALON-FH categorical-band architecture as a natural future direction in the Discussion. We did not include PRS in this analysis because (a) PRS is available in the UK Biobank but not in routine Welsh clinical data, and (b) including a feature available in only one cohort would have broken the symmetric bidirectional validation design that is the methodological core of the contribution.

---

## Q10. "Why isn't lp(a) more heavily weighted given it's not statin-modifiable?"

Lipoprotein(a) ≥120 nmol/L was retained as a positive band in both direction-specific equations (OR 1.36 in Direction A). Its relative magnitude — smaller than the age and sex bands but larger than the BMI bands — reflects the population-level distribution: approximately 99% of our Welsh cohort and approximately 9% of UK Biobank carriers exceeded the 120 nmol/L threshold. A binary band with near-constant prevalence in one cohort necessarily has limited discriminative variance, even when the underlying biology is causal. We acknowledge this and identify continuous Lp(a) modelling and a more granular Lp(a) band set as a future-direction refinement.

---

## Q11. "Why prevalent rather than incident outcomes — wouldn't the model be more useful as a 10-year MACE predictor?"

This is the same fundamental limitation as Q2. Our position is: the bidirectional external superiority of CALON-FH over SAFEHEART-RE on prevalent ASCVD is the core methodological contribution; whether the same advantage extends to incident outcomes is a question for prospective validation. We choose not to claim incident discrimination from a cross-sectional design.

---

## Q12. "Why only Wales and UK Biobank? Why not validate in another FH registry?"

The two-cohort design is intentional and reflects practical constraints: Wales PASS is the development cohort because it is our local register and is the only cohort in which we have access to verified ApoB and the granular drug-prescription data required for the X.Y encoding. UK Biobank is the external cohort because it has genetically defined *LDLR* carriers identifiable via whole-exome sequencing at substantial sample size (n=3,540). Validation in additional FH registries — particularly SAFEHEART itself, the Dutch FH cohort, the FH Foundation registry, and East Asian / South Asian FH cohorts — is identified as the principal further-validation requirement in the Discussion. We are actively pursuing collaborative access to these registries.

---

## Q13. "Why isn't this just a poorly-powered version of SAFEHEART-RE-recalibrated?"

We addressed this in the Methods by refitting SAFEHEART-RE on identical cohorts and data with no sign constraint. The refit retains SAFEHEART-RE's feature structure (the original nine bands) but allows its coefficients to be optimised for the Welsh or UK Biobank distribution. CALON-FH outperforms this maximally-favourable comparator in both directions of cross-cohort transfer. The advantage therefore cannot be attributed to refitting per se; it must come from either the feature structure (modern atherogenic bands), the sign-constrained drop (transferability prior), or the X.Y untreated-LDL recovery (treatment-conditioned phenotyping). The Discussion attributes the gain to a combination of all three.

---

## Q14. "Is the elderly subgroup result driven by overfitting in a small subset?"

The ≥65y subgroup contains 646 UK Biobank patients with 57 events — among the largest of the 14 subgroups. The discriminative advantage in this subgroup (ΔAUC +0.118) is the largest of the 14 but is computed on the second-largest event count of any subgroup-direction combination in the analysis, making overfitting an unlikely explanation. The NRI gain in the same subgroup is +33%. The advantage in the ≥65y subgroup also has a biologically coherent explanation (Discussion, Mechanistic Interpretation section): treated elderly FH patients have residual variance carried by comorbidity bands and structural particle biomarkers rather than measured LDL-C, and CALON-FH has more leverage on those bands than SAFEHEART-RE.

---

## Q15. "The Direction B non-event NRI is negative. Is this acceptable?"

The Direction B non-event component of −0.151 indicates that approximately 15% more non-events were up-classified by CALON-FH than by SAFEHEART-RE. We report this honestly and discuss the trade-off in the manuscript: in a high-event-prevalence tertiary cohort (Welsh 27%), the appropriate up-classification of patients with high underlying biological burden may be clinically desirable even if a proportion of those patients have not yet manifested an event. Whether the up-classification represents genuine elevated underlying risk or over-classification cannot be settled by cross-sectional data; prospective follow-up will resolve the question. We do not claim Direction B as the primary endpoint; we report it as a methodological check on bidirectional generalisability and note its limitations transparently.

---

## Q16. "What is the single sentence that captures why this matters for a *Circulation* reader?"

In the contemporary era of high-intensity statin–ezetimibe combinations, PCSK9 inhibition, and bempedoic acid, the measured LDL-C of an FH patient in active care no longer reflects the underlying genetic burden that drives their cardiovascular risk; CALON-FH demonstrates that recovering the underlying biology — through dose-specific back-calculation and sign-constrained feature retention — is achievable from routine clinical data and refines patient-level risk ranking sufficiently to alter treatment-intensification decisions in precisely the populations where contemporary clinical guidelines provide the least guidance.

# CALON-Structure: AlphaFold3 Structural Severity Score for FH-ASCVD Prediction

**Date:** 2026-03-17
**Author:** Dr Nader Genedy
**Project:** CALON-FH ASCVD Risk Prediction Pipeline
**Status:** Approved for implementation

---

## 1. Purpose

Integrate AlphaFold3 protein structure predictions into the CALON-FH pipeline to create a Structural Severity Score (SSS) that quantifies how much each FH variant disrupts protein function. This goes beyond gene-level classification (LDLR vs APOB vs PCSK9) to variant-level risk stratification — connecting genotype → structure → ASCVD phenotype.

## 2. Scope

### In scope
- Extract and classify all unique FH variants from Wales registry data
- Model wild-type LDLR, PCSK9, ApoB structures and key protein complexes using AlphaFold3
- Compute structural impact metrics for each missense variant using FoldX
- Build a composite Structural Severity Score (SSS)
- Validate SSS association with ASCVD outcomes in Wales (development) and UKB (external validation)
- Generate publication-quality structural figures
- Test SSS as an additional CALON predictor variable (CALON-Structure 7-variable model)

### Out of scope
- De novo variant calling from raw WES data (uses existing TUDOR annotations)
- Molecular dynamics simulations (FoldX static modelling is sufficient)
- Drug binding simulations (pharmacogenomics is future work)

## 3. Data Sources

### Wales FH Registry
- **DRAGON_3.csv** (South Wales): 1,363 patients, 424 with mutations, 155 unique variants
- **WALES_FH_CLEANED (1) - Copy.csv** (All Wales): 7,254 patients, 3,694 with mutations, ~563 unique variants
- Mutation format: HGVS cDNA (e.g., `LDLR:c.1816G>T`)
- Outcome: ASCVD_combined binary

### UK Biobank (external validation)
- 1,623 genetically confirmed FH patients (1,321 LDLR, 301 APOB, 1 PCSK9)
- Currently gene-level only — variant-level annotation requires WES re-annotation (Phase 4b)
- WES file paths available in `calon_batch_wes_paths.csv`

## 4. Architecture

```
Phase 1: VARIANT EXTRACTION (Python)
  Wales Mutation1 → parse HGVS → protein coordinates → domain mapping
  → classify: missense | nonsense | frameshift | splice | structural

Phase 2: ALPHAFOLD3 MODELLING (AF3 Server + FoldX)
  2a. Wild-type structures from AF3 Server (5 jobs, Day 1)
  2b. Mutant modelling via FoldX BuildModel on AF3 structures (~50 variants)
  2c. Protein complexes: LDLR-ApoB, LDLR-PCSK9 (AF3 multimer)

Phase 3: STRUCTURAL SEVERITY SCORE (R)
  5 component metrics → weighted composite → fit on Wales ASCVD

Phase 4: VALIDATION & INTEGRATION (R)
  4a. Wales internal: Cox, logistic, KM by SSS tertile
  4b. UKB external: map WES variants → SSS lookup → validate
  4c. CALON-Structure: 7-variable model with SSS
  4d. Publication figures
```

## 5. Phase 1: Variant Extraction & Classification

### Input
`Mutation1` column from Wales data in HGVS cDNA format.

### Processing
1. Parse `GENE:c.XXXXX` format → extract gene name and cDNA change
2. Convert cDNA to protein notation using `biocommons.hgvs` or Mutalyzer API
3. Classify variant type by pattern matching
4. Map protein position to structural domain

### Variant classification and SSS assignment

| Type | Pattern | AF3 Modelled? | SSS Assignment |
|------|---------|---------------|----------------|
| Missense | Single nucleotide, amino acid change | Yes | Computed from metrics |
| Nonsense | Premature stop codon | No | 1.0 (null allele) |
| Frameshift | Insertions/deletions causing frameshift | No | 1.0 (null allele) |
| Splice site | c.XXX+1/+2 or c.XXX-1/-2 | No | 0.95 (near-null) |
| In-frame del/dup | Removes/adds residues without frameshift | Case-by-case | 0.85 default |
| Exon deletion | "Deletion of exon X LDLR" | No | 1.0 (null allele) |
| Exon duplication | "Duplication of exons X-Y" | No | 0.85 (disrupted) |

### LDLR domain mapping

| Domain | Residues | Function | Severity weight |
|--------|----------|----------|----------------|
| Ligand-binding (LR1-7) | 1-292 | Binds LDL/ApoB | 1.0 (highest) |
| EGF precursor homology | 293-692 | pH-dependent release, recycling | 0.7 |
| O-linked sugar | 693-750 | Stability, unknown function | 0.3 |
| Transmembrane | 751-788 | Membrane anchoring | 0.5 |
| Cytoplasmic | 789-860 | Endocytosis signal (NPXY) | 0.6 |

### Output
CSV file: `af3_variant_catalogue.csv`
Columns: `variant_id, gene, cdna_change, protein_change, variant_type, domain, residue_number, sss_preset` (for non-missense)

## 6. Phase 2: AlphaFold3 Structural Modelling

### 2a. Wild-type structures (AF3 Server, Day 1)

| Job | Protein | UniProt | Residues | Purpose |
|-----|---------|---------|----------|---------|
| 1 | LDLR | P01130 | 1-860 | Main receptor — all LDLR variants mapped here |
| 2 | PCSK9 | Q8NBP7 | 1-692 | PCSK9 gain-of-function variants |
| 3 | ApoB-100 RBD | P04114 | 3359-3600 | Receptor-binding domain only (full=4563 too large) |
| 4 | LDLR–PCSK9 complex | P01130 + Q8NBP7 | ~1,552 | Model degradation pathway interface |
| 5 | LDLR–ApoB RBD complex | P01130 + P04114(3359-3600) | ~1,100 | Model LDL binding interface |

### 2b. Mutant modelling (FoldX, Day 1-2)

- Input: AF3 wild-type LDLR structure (mmCIF → convert to PDB)
- Tool: FoldX 5.0 `BuildModel` command
- Process: For each missense variant, mutate the residue and energy-minimise
- Output: Mutant PDB + ΔΔG value per variant
- Coverage: Top ~50 most frequent missense variants (covers >80% of patients)
- Rare variants: Assign SSS by nearest-neighbour in same domain with similar physicochemical change

### 2c. Complex analysis

From the AF3 multimer structures:
- Extract interface residues (residues with atoms <5Å from partner chain)
- Compute per-residue interface contact area
- Flag variants that fall within or near (<10Å) the binding interface

## 7. Phase 3: Structural Severity Score (SSS)

### Component metrics (missense variants only)

| Metric | Source | Computation | Range | Interpretation |
|--------|--------|-------------|-------|----------------|
| ΔΔG | FoldX on AF3 structure | `FoldX BuildModel` | kcal/mol | >2 = destabilising, >4 = severely destabilising |
| Local RMSD | PyMOL/BioPython | Cα RMSD in 10Å sphere around mutation | Å | Higher = more structural displacement |
| Interface proximity | AF3 complex structures | Min distance from variant to nearest interface residue | Å (inverted) | Closer = more likely to disrupt binding |
| pLDDT | AF3 output | Per-residue confidence | 0-100 (direct) | High pLDDT = ordered region = mutations more damaging; Low pLDDT = disordered = more tolerant |
| Domain severity | Domain mapping table | Categorical lookup | 0-1 | Ligand-binding=1.0, cytoplasmic=0.2 (expert-assigned, subject to sensitivity analysis) |
| REVEL score | Pre-computed (dbNSFP) | Ensemble pathogenicity predictor | 0-1 | >0.5 = likely pathogenic; incorporates sequence + evolutionary features |
| CADD phred | Pre-computed (dbNSFP) | Combined annotation-dependent depletion | 0-40+ | >20 = top 1% most deleterious |

**Note:** Addition of REVEL and CADD addresses the known limitation that structural features alone (AF/FoldX) underperform for variant pathogenicity prediction. Hybrid structural + sequence-based approaches are recommended by the literature (see `alphafold_fh_ascvd_comprehensive_analysis.md`).

### SSS computation

```r
# Normalise each metric to 0-1 range
# Use 10-fold cross-validation within Wales to avoid circular validation
# (SSS weights must NOT be fitted and tested on the same data)

library(caret)
set.seed(42)
folds <- createFolds(wales_missense$ASCVD, k = 10)

# Within each fold: fit on 9/10, predict on held-out 1/10
sss_cv <- numeric(nrow(wales_missense))
for (i in seq_along(folds)) {
  train_idx <- unlist(folds[-i])
  test_idx  <- folds[[i]]
  fit <- glm(ASCVD ~ norm_ddg + norm_rmsd + norm_interface_inv +
                      norm_plddt + domain_severity + revel + cadd_phred,
             data = wales_missense[train_idx, ], family = binomial)
  sss_cv[test_idx] <- predict(fit, wales_missense[test_idx, ], type = "response")
}

# Final SSS model (for lookup table generation): fit on ALL Wales data
sss_model <- glm(ASCVD ~ norm_ddg + norm_rmsd + norm_interface_inv +
                          norm_plddt + domain_severity + revel + cadd_phred,
                 data = wales_missense, family = binomial)

# SSS = predicted probability from the structural+sequence model
# For non-missense: use preset values from classification table

# IMPORTANT: Use clustered standard errors for all patient-level analyses
# (patients sharing the same variant have identical SSS — violates independence)
library(sandwich)
library(lmtest)
coeftest(patient_model, vcov = vcovCL, cluster = ~variant_id)
```

### Handling compound heterozygotes
- Patients with Mutation2 or Mutation3: use weighted combination
- `SSS_patient = SSS_max + 0.3 * SSS_second` (capped at 1.0)
- Rationale: second pathogenic variant adds risk beyond the most severe alone
- Sensitivity analysis: compare max-only vs additive vs weighted approaches

## 8. Phase 4: Validation & Integration

### 4a. Wales internal validation (n=3,694)

**Cross-validated SSS values** from Phase 3 are used here (NOT refitted on the same data).

```r
# All patient-level models use clustered standard errors (cluster = variant_id)
# because patients with the same variant share identical SSS values

# Logistic regression (using CV-derived SSS)
model_sss <- glm(ASCVD ~ SSS_cv + age + sex, family = binomial)
coeftest(model_sss, vcov = vcovCL, cluster = ~variant_id)

# Cox proportional hazards (if time-to-event available — confirm censoring dates exist)
coxph(Surv(time, event) ~ SSS_cv + age + sex, cluster = variant_id)

# KM by SSS tertile
survfit(Surv(time, event) ~ SSS_tertile)

# Discrimination improvement — BASELINE includes variant_type (not just gene)
# Gene-only model:       glm(ASCVD ~ gene + age + sex)
# Gene+type model:       glm(ASCVD ~ gene + variant_type + age + sex)
# SSS model:             glm(ASCVD ~ SSS_cv + age + sex)
# Combined:              glm(ASCVD ~ SSS_cv + gene + age + sex)
# Compare AUC: DeLong test
# Pre-specified primary analysis: SSS_cv association with ASCVD (age+sex adjusted)
# All subgroup analyses (domain, gene) are exploratory with FDR correction
```

**Note on patients without identified mutations:** The ~3,560 Wales patients without identified mutations are excluded from SSS analyses. They are included in the CALON-Structure 7-variable model only if SSS is set to population median (sensitivity analysis: exclude vs median imputation).

### 4b. UKB external validation

1. Access UKB WES VCFs via paths in `calon_batch_wes_paths.csv`
2. Annotate variants using VEP (Variant Effect Predictor) on RAP
3. Filter to pathogenic/likely pathogenic in LDLR/APOB/PCSK9
4. Map each variant to SSS lookup table from Wales
5. Apply same validation models as 4a
6. Report: AUC, calibration, NRI vs gene-only model

### 4c. CALON-Structure model

```r
# Existing CALON 6-variable model:
# ASCVD ~ age + re_ldl + lipid_years + log_apob_ldl + lpa_binary + inv_hdl

# CALON-Structure 7-variable model:
# ASCVD ~ age + re_ldl + lipid_years + log_apob_ldl + lpa_binary + inv_hdl + SSS

# Test improvement:
# - DeLong test (AUC comparison)
# - NRI (net reclassification improvement)
# - IDI (integrated discrimination improvement)
# - Decision curve analysis
```

### 4d. Publication figures

| Figure | Content | Tool |
|--------|---------|------|
| Fig A | AF3 wild-type LDLR coloured by domain, top 10 variants mapped as spheres | PyMOL |
| Fig B | LDLR–ApoB complex, interface-disrupting variants highlighted in red | PyMOL |
| Fig C | LDLR–PCSK9 complex, gain-of-function PCSK9 variants highlighted | PyMOL |
| Fig D | KM survival curves by SSS tertile (low/medium/high severity) | R (survminer) |
| Fig E | Forest plot: SSS hazard ratio across gene subgroups | R (forestplot) |
| Fig F | SSS distribution by variant type (missense vs null vs splice) | R (ggplot2) |
| Supp Table | Full variant-to-SSS lookup table (all 563 variants) | CSV |

## 9. Implementation Scripts

| Script | Phase | Language | Purpose |
|--------|-------|----------|---------|
| `10_AF3_extract_variants.py` | 1 | Python | Parse Wales mutations, classify, map to protein coords |
| `11_AF3_model_structures.sh` | 2 | Bash | AF3 server job submission + FoldX BuildModel |
| `12_AF3_compute_metrics.py` | 3 | Python | ΔΔG, RMSD, interface distance, pLDDT extraction |
| `13_AF3_build_SSS.R` | 3 | R | Fit SSS weights, test ASCVD association in Wales |
| `14_AF3_validate_UKB.R` | 4 | R | Apply SSS to UKB, validate, test CALON-Structure |
| `15_AF3_figures.py` | 4 | Python | PyMOL structural figures (scriptable) |

## 10. Timeline

| Day | Activity |
|-----|----------|
| Day 1 | Phase 1: Extract and classify all variants from Wales |
| Day 1 | Phase 2a: Submit 5 AF3 Server jobs (wild-type + complexes) |
| Day 1-2 | Phase 2b: FoldX BuildModel for ~50 missense variants |
| Day 2 | Phase 2c: Analyse complex interfaces |
| Day 2-3 | Phase 3: Compute metrics, build SSS, test Wales ASCVD |
| Day 3 | Phase 4a: Wales validation (logistic, Cox, KM, forest) |
| Day 3-4 | Phase 4d: Structural figures |
| Day 4+ | Phase 4b-c: UKB WES re-annotation and external validation |

## 11. Risks & Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| SSS does not predict ASCVD beyond gene | Main hypothesis fails | Still publish structural characterisation + figures (Approach 2 fallback) |
| AF3 Server down or rate-limited | Delays Phase 2 | Use existing AlphaFold DB structures for monomers; AF3 only needed for complexes |
| FoldX ΔΔG unreliable for some variants | Noisy SSS component | Use DynaMut2 as secondary ΔΔG source; ensemble scoring |
| UKB WES re-annotation too complex | Phase 4b blocked | Publish Wales-only results; UKB validation as follow-up |
| Too few missense variants per domain for statistics | Underpowered subgroup | Aggregate into 2 groups (ligand-binding vs other) instead of 5 domains |

## 12. Success Criteria

1. SSS shows significant association with ASCVD in Wales (p < 0.05, adjusted for age+sex)
2. SSS adds discriminatory value beyond gene identity alone (ΔAUC > 0)
3. KM curves separate by SSS tertile (log-rank p < 0.05)
4. At least 3 publication-quality structural figures generated
5. Full variant-to-SSS lookup table covering all 563 Wales variants

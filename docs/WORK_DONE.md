# Work Done — the CALON-FH / TUDOR programme, from scratch

A complete narrative of what has been built, why, and how the pieces fit. Read this to
understand the project end-to-end before touching any code. Numbers here are the validated
headline values; trace any specific claim to a CSV in `data/` before re-quoting it.

---

## 0. The problem

Familial hypercholaesterolaemia (FH) is the commonest monogenic cardiovascular disorder
(~1 in 250). Two clinical problems motivate this work:

1. **Diagnosis** — most FH carriers are never identified. Existing clinical scores (DLCN,
   Simon Broome, MEDPED, FAMCAT) rely on phenotype that is blurred by treatment. → **TUDOR**.
2. **Severity & prognosis** — not all *LDLR* variants are equal, and not all confirmed-FH
   patients carry the same ASCVD risk. Variant pathogenicity is usually binary (path/benign),
   discarding quantitative loss-of-function information. → **CALON-FH / SSS v3**.

The unifying idea across the programme is **intervention-conditioned phenotyping**: once a
patient is treated, the biomarker (LDL) no longer reads baseline biology — it reads
*response*. Genetic and structural information must therefore be recovered from the genome and
the prescription record, not the treated lipid panel.

---

## 1. CALON-FH — Structural Severity Score (SSS v3)

### 1.1 What it is
A **9-layer, a-priori-weighted composite** scoring the functional severity of *LDLR* missense
variants, computed for **16,340 variants**. Layers and weights:

| Layer | Weight | Source |
|---|---|---|
| DMS Uptake (Tabet 2025) | 0.25 | deep mutational scanning, LDL uptake |
| Boyd & Goldstein domain class | 0.18 | domain-level prior |
| FoldX ddG | 0.15 | AlphaFold3 structure → folding stability |
| AlphaMissense | 0.12 | pathogenicity predictor (AUROC 0.94) |
| DMS Abundance (Tabet) | 0.10 | surface protein abundance |
| Islam cell activity | 0.10 | functional cell assay |
| AlphaGenome L2 | 0.05 | regulatory effect |
| AlphaGenome Splice | 0.05 | splicing effect |

The map lives in `data/alphafold/analysis/ldlr_therapeutic_map_16340.csv` (16,340 × 17) and
`data/alphafold/analysis/sss_v2_catalogue.csv`. Per-variant UKB carrier scores:
`data/alphafold/analysis/ukb_carriers_with_sss.csv` (36,217 × 17) and the annotated carrier
file `ukb_carriers_annotated.csv` (3,094 × 34).

### 1.2 Structural pipeline
AlphaFold3 models of LDLR (and mutants) → FoldX / Rosetta saturation mutagenesis → ddG per
substitution → folded into the SSS composite. The structural code is in `code/python/`
(`37_foldx_saturation_mutagenesis.py`, `38_parse_af3_mutants.py`, `40_alphafold3_mutant_sequences.py`,
`54_LPA_alphafold3_foldx.py`, `63_LPA_foldx_chimerax_pipeline.py`, `64_run_lpa_saturation_parallel.py`)
with structural outputs and figures under `assets/alphafold/`.

### 1.3 Clinical validation (Paper 1 — Orthogonal Decomposition)
Validated in UK Biobank (3,544 *LDLR* coding-variant carriers; 3,085 missense) and Wales
(Dragon-3 / All-Wales PASS). Headline findings:

| Finding | Value | P | n |
|---|---|---|---|
| Production R² (untreated LDL ~ SSS) | 0.54 | <10⁻³⁰⁰ | 2,398 |
| SSS vs treated LDL (UKB) | rho=0.083 | 0.020 | 775 |
| SSS vs LDL reduction (Wales) | rho=0.234 | 0.014 | 109 |
| Statin response is SSS-dependent | rho=0.244 | 0.027 | 82 |
| Ezetimibe is SSS-independent | rho=0.051 | 0.80 | 27 |
| SSS ⊥ polygenic risk (orthogonality) | r=−0.017 | NS | 2,228 |
| CHIP × SSS ASCVD interaction | 38.5% | 0.0004 | 52 |
| HDL sub-3 (CETP axis) | rho=+0.077 | 0.0006 | 1,988 |
| HDL sub-4 (CETP axis) | rho=−0.087 | 0.0001 | 1,988 |

**Interpretation (data-supported):** SSS explains most variance in *untreated* LDL (R²=0.54),
is orthogonal to polygenic risk, and predicts statin (but not ezetimibe) response — consistent
with the orthogonal-decomposition thesis that structural severity is a distinct, treatable axis.

Status: **submission-ready** (6,860 words, 7 figures, 6 tables; QC'd). Tables in
`data/Table1*.csv … Table5_biexternal_validation.csv`; model outputs in
`data/CALON_*_odds_ratios.csv`, `CALON_X3_PRS_coefficients.csv`.

### 1.4 Generalisability & extensions
- **BRCA1** generalisability test (`51_brca1_generalisability.py`; data in
  `data/alphafold/analysis/brca1/`) — applies the structural-score framework beyond LDLR.
- **Lp(a)** structural and epidemiological sub-studies (`54–73_*.py`) — AlphaFold3/FoldX of
  apo(a), aortic-stenosis stratification, ancestry analyses.

---

## 2. TUDOR — UK Biobank FH diagnostic algorithm

### 2.1 What it is
A **diagnostic** risk algorithm that flags probable FH from routinely available lipids
(threshold-based, no family history), built and validated across three cohorts (Wales,
South Wales, UK Biobank). Model: **Elastic Net**, with **leave-one-cohort-out (LOCO)** CV.
Distinct from CALON: TUDOR *diagnoses FH*; CALON *predicts ASCVD in confirmed FH*.

### 2.2 Data & code
- Predictions: `data/tudor_loco_output/loco_predictions_complete.csv` (113,538 × 18; FH+=3,136).
- Full-cohort supplements: `data/tudor_loco_output/tudor_flaw_results/` —
  `ukb_supplement_clean.csv` (501,936 × 26), `predictions_with_supplement.csv` (501,936 × 23),
  `ukb_supplement_lc_cohort*.csv` (lipid-clinic-eligible secondary cohort).
- R modelling pipeline: `code/R/` (`01_*.R`–`16_*.R`, TUDOR scorers, figures).
- Python: `code/python/CALON_FINAL*.py`, `TUDOR_*.py`, head-to-head + traceability scripts.
- Reproduced coefficients/AUC: `data/TUDOR_reproduced_coefficients.csv`,
  `data/TUDOR_reproduced_AUC.csv`; live ledger `data/TUDOR_LIVE_LEDGER.csv`.

### 2.3 Validated metrics (traceability ledger)
| Metric | Value | Status |
|---|---|---|
| Wales index AUC | 0.7585 | PASS |
| Wales cascade AUC | 0.791 | PASS |
| Wales full-cohort AUC (LOCO Fold 2) | live 0.7816 (ms 0.7725) | DRIFT (tracked) |

Full ledger: `docs_md/TUDOR_TRACEABILITY_FINAL.md` (11 PASS / 16 DRIFT / 2 NO_DATA / 29 total),
`data/TUDOR_LIVE_LEDGER.csv`. DRIFT items are documented, not hidden.

### 2.4 Submission history
Submitted to *Journal of Clinical Lipidology* (**JCLINLIPID-D-25-01142**), through R1→R2→R3.
The R3 bundle added a **frozen-coefficient FAMCAT external validation**, a **cascade-leakage
discovery**, and a full **TRIPOD-AI head-to-head** vs eDLCN / FAMCAT-approx / MEDPED /
Simon Broome (comparators as *published rules*, not refit on our data — `comparator-model-validation`
discipline). Submission package index: `docs_md/TUDOR_SUBMISSION_PACKAGE_README.md`. Manuscripts,
cover letters, response-to-reviewers, tables, highlights: `manuscripts/TUDOR_*`.

### 2.5 Sub-studies (recent, per git history)
- **LDL-C equation** sub-study — pre-treatment LDL equation agreement (Bland-Altman +
  Passing-Bablok), Wales v1 baseline n=5,364.
- **Temporal-trends** sub-study (design spec committed).
- MEDPED / Simon Broome scorers, Little MCAR gate + MICE imputation.

---

## 3. CALON-FH Atlas — interactive website

Deployed at https://nadergenedy.github.io/calon-fh-atlas/. Source in `website/`
(`index.html`, `variant_pen60.html`, `website_data.json`, `website_data_v3b.json`). Features:
traffic-light variant classification, penetrance curves, LLT (lipid-lowering therapy) response,
VUS reclassification, and a 3D Mol* structural viewer. Build scripts:
`code/python/41_build_website_v3b.py`, `39_prepare_website_data.py`, `40_atlas_v3b_literature_update.py`.

---

## 4. How the code is organised

Scripts are numbered roughly in execution order:

- **`code/R/01_*.R` … `16_*.R`** — the core R pipeline: cohort assembly, TUDOR Elastic Net,
  LOCO CV, calibration, figures (ggplot2 + patchwork).
- **`code/python/17_*.py` … `36_*.py`** — clinical analysis, cohort processing, validation.
- **`code/python/37_*.py` … `53_*.py`** — AlphaFold3/FoldX structural pipeline, SSS construction,
  website build, reviewer analyses, ClinVar validation, head-to-head pathogenicity, BRCA1.
- **`code/python/54_*.py` … `73_*.py`** — Lp(a) structural + epidemiological deep-dives.
- **`code/python/CALON_*.py`** — the CALON clinical models (full / lite / perfect / external).
- **`code/python/TUDOR_*.py`** — TUDOR reproducers, traceability, R1→R2 pipeline.
- **`code/shell/*.sh`** — UKB-RAP extraction (`00_CALON_extract_ukbrap.sh`, `07d_extract_lpa_correct.sh`).

Profile of every data file: `docs/DATA_DICTIONARY.md`. UKB field meanings:
`docs_md/CALON_UKB_Field_Dictionary.md`.

---

## 5. Methodological backbone (why results are trustworthy)

Every analysis honours the non-negotiables in `docs/INSTRUCTIONS.md`:
family-level deduplication (Dragon-3 ⊂ All-Wales PASS), NoAgeLDL sensitivity variant,
ICD-10/OPCS code-list completeness (incl. K611), TRIPOD reporting, BH-FDR for >10 tests,
E-values for headline effects, leak-free splits, and a **traceability ledger** asserting every
manuscript number against a live CSV. Bugs trigger an end-to-end locked rerun, not a hot-patch.

The historic motivating bug: UKB first-occurrence date fields p131286–p131294 were initially
mis-mapped (coding I10–I20 hypertension rather than the intended ASCVD I-codes). This is now a
standing audit check (`tudor-qc`, `ukb-data-audit`).

---

## 6. Current status (2026-05-30)

| Project | Status |
|---|---|
| CALON-FH Paper 1 (Orthogonal Decomposition) | Submission-ready; QC'd |
| TUDOR diagnostic algorithm | Submitted JCL, R2/R3 in progress |
| CALON-FH Atlas | Deployed |
| Lp(a) structural sub-studies | In progress (exploratory) |
| BRCA1 generalisability | Exploratory |

---

## 7. Reproducing a result

1. Read `CLAUDE.md` (root) for context and the headline numbers.
2. Find the data file in `docs/DATA_DICTIONARY.md`.
3. Open the relevant script in `code/`; update its input path to this bundle's `data/`
   (see the Path note in `README.md`).
4. Run with Python 3.12 / R 4.5.x.
5. Cross-check output against the traceability ledgers (`data/*_TRACEABILITY*.csv`,
   `data/TUDOR_LIVE_LEDGER.csv`, `data/CALON_PAPER1_TRACE_LEDGER.csv`).

# CLAUDE.md — CALON-FH / TUDOR project context (paste-ready)

> Paste this whole file into a new Claude conversation to bootstrap full project context.
> It is a self-contained extract of Dr Genedy's working context, scoped to this bundle.

## Who

**Dr Nader Genedy** — Cardiologist and clinical researcher, University Hospital of Wales,
Cardiff. Specialist in familial hypercholesterolaemia (FH), structural biology, and precision
cardiovascular medicine. British English throughout. Communication: direct and clinical — no
filler, show data and numbers before interpretation, every claim carries a number + CI + P-value.

## The work in this bundle

Three linked projects:

1. **CALON-FH** — a 9-layer **structural severity score (SSS v3)** for LDLR missense variants
   (16,340 variants scored), validated clinically in UK Biobank (3,544 coding-variant carriers;
   3,085 missense) and Wales cohorts (Dragon-3 / All-Wales PASS). Core thesis: SSS captures
   LDLR loss-of-function severity orthogonally to polygenic risk and predicts statin response.
2. **TUDOR** — a UK Biobank lipid-threshold FH **diagnostic** algorithm (Elastic Net, leave-one-
   cohort-out CV across Wales / South Wales / UKB). Submitted to *J Clin Lipidol* (R2/R3).
   Distinct from CALON: TUDOR *diagnoses* FH; CALON *predicts ASCVD* in confirmed FH.
3. **CALON-FH Atlas** — deployed interactive website (`website/`).

## Validated headline findings (cite with these exact numbers)

**CALON / SSS v3 (Paper 1 — Orthogonal Decomposition):**
| Finding | Value | P | n |
|---|---|---|---|
| Production R² (untreated LDL) | 0.54 | <10⁻³⁰⁰ | 2,398 |
| SSS vs treated LDL (UKB) | rho=0.083 | 0.020 | 775 |
| SSS vs LDL reduction (Wales) | rho=0.234 | 0.014 | 109 |
| Statin response SSS-dependent | rho=0.244 | 0.027 | 82 |
| Ezetimibe SSS-independent | rho=0.051 | 0.80 | 27 |
| PRS ⊥ SSS (orthogonality) | r=−0.017 | NS | 2,228 |
| CHIP × SSS ASCVD | 38.5% | 0.0004 | 52 |

**TUDOR (diagnostic algorithm):**
- Wales index AUC 0.7585; Wales cascade AUC 0.791 (both PASS in traceability ledger).
- LOCO validation across Wales / South Wales / UKB; Elastic Net.
- Head-to-head vs eDLCN / FAMCAT-approx / MEDPED / Simon Broome (published rules, not refit).
- Submitted JCLINLIPID-D-25-01142 (R2/R3); frozen-coefficient FAMCAT external validation.

> When unsure whether a number is current, trace it to a CSV in `data/` — do not quote prose.
> The TUDOR traceability ledger (`data/TUDOR_LIVE_LEDGER.csv`, `docs_md/TUDOR_TRACEABILITY_FINAL.md`)
> lists PASS/DRIFT status per claim.

## SSS v3 architecture (a-priori weights)

| Layer | Weight | Coverage |
|---|---|---|
| DMS Uptake (Tabet) | 0.25 | 94.9% |
| B&G Domain Class | 0.18 | 100% |
| FoldX ddG (AF3) | 0.15 | 62.5% |
| AlphaMissense | 0.12 | 100% |
| DMS Abundance (Tabet) | 0.10 | 91.9% |
| Islam Cell Activity | 0.10 | 23.3% |
| AlphaGenome L2 | 0.05 | 26.0% |
| AlphaGenome Splice | 0.05 | 26.0% |

## Data in this bundle

All CSVs are in `data/` (266 files, paths preserved). Key files:
- `data/CALON_clinical_analysis_dataset.csv` — 428,406 × 21, clinical analysis cohort.
- `data/output/v2/CALON_FH_perfect_dataset.csv` — 8,110 × 66, modelling dataset.
- `data/tudor_loco_output/loco_predictions_complete.csv` — 113,538 × 18, TUDOR LOCO predictions (FH+=3,136).
- `data/alphafold/analysis/ldlr_therapeutic_map_16340.csv` — 16,340 × 17, SSS map.
- `data/DRAGON_3.csv` — Wales discovery cohort.
- Table1*/Table2*/Table5* — manuscript tables. *_odds_ratios / *_coefficients / *_results — model outputs.

Full profile: `docs/DATA_DICTIONARY.md`. Field meanings: `docs_md/CALON_UKB_Field_Dictionary.md`.

## Statistical non-negotiables (apply by default — see docs/INSTRUCTIONS.md)

1. **Family-level deduplication** (South Wales/Dragon-3 ⊂ All-Wales PASS by FamilyNumber).
2. **NoAgeLDL sensitivity variant** for every FH risk model.
3. **ICD-10/OPCS code-list completeness** — ASCVD composite = I20–I25, I63, G45, I70 + first-occurrence
   fields p131286/8/90/2/4/6. **K611 (balloon valvuloplasty)** mandatory for severe-AS. Never use I35 for ASCVD.
4. **TRIPOD** reporting: C-stat + 95% CI, calibration intercept/slope, Brier, NRI, IDI, DCA.
5. **BH-FDR** for >10 simultaneous tests; **E-value** for every headline OR/HR (robust if >2.0).
6. **Leak-free splits** — train/validation share zero family IDs.
7. **Validate every numerical claim against source CSV** before "submission-ready".

## Style

British English (colour, analyse, organisation). Nature/Lancet/NEJM register. One claim per
sentence, one hedge word maximum. Distinguish "supported by data" / "consistent with data" /
"speculative" explicitly. Report numbers first, mechanisms after. Realistic journal targets:
JCL, Circulation, JAMA Cardiology (not Nature Med/NEJM without multi-cohort/mechanistic data).

## Environment

Python **3.12** only. R 4.5.x. Windows console is cp1252 — never print Unicode arrows/em-dashes;
use ASCII (`->`, `--`). Always `open(..., encoding='utf-8')`. Dual-format deliverables (md + docx).

## Skills

- Project skills (tudor-qc, cox-analysis, manuscript-qc, ...): `skills/`, documented in `docs/SKILLS.md`.
- Google DeepMind Science Skills (40, pinned v1.2.1): `.claude/skills/` -- ClinVar, gnomAD, Ensembl,
  UniProt, AlphaFold, AlphaGenome, PubMed/EuropePMC/OpenAlex, ClinicalTrials.gov, ChEMBL, etc.
  Run via `uv run scripts/<x>.py`; keys in `~/.env` (see `.claude/skills/credentials/SKILL.md`).
  Rebuild / verify with `python build_science_skills.py`; index in `.claude/skills/SCIENCE_SKILLS_INDEX.md`.

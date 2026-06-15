# Projects Overview — everything in this bundle

This bundle is **not just TUDOR / CALON pipeline**. It now covers the whole research programme.
The original pipeline repo (`calon_ukb_pipeline`) sits at the bundle **root** (`data/`, `code/`,
`manuscripts/`, `docs_md/`, `assets/`, `website/`). Every **other** project is mirrored under
**`projects/<ProjectName>/`** with its internal structure preserved (code + docs + CSV ≤ 300 MB).

Giant raw extracts (> 300 MB CSV) are **indexed, not copied** — see `docs/LARGE_DATA_INDEX.csv`
(project, relative path, size, original D: path). They remain on D: at their source location.

---

## Root project — CALON-FH / TUDOR pipeline
`github.com/NaderGenedy/calon-ukb-pipeline`. SSS v3 structural scoring, TUDOR diagnostic
algorithm, CALON-FH Atlas website. Full detail in `docs/WORK_DONE.md`.

---

## Projects under `projects/`

### Structural / SSS core
- **`CALON_AlphaFold_Rebuild/`** — the central data + analysis hub. SSS v3 construction,
  16,340-variant therapeutic map, UKB carrier files, Paper1/2/3 scaffolds, figure data.
  (Largest CSV footprint — raw UKB masters and the Paper3 GP-prescription extract are indexed,
  not copied.)
- **`AlphaFold_SSS/`** — AlphaFold3 structural models and SSS scoring artefacts.
- **`MD_Thesis_AlphaFold_Cardio/`** — the MD-by-Research thesis (AlphaFold + cardiology).

### NOBEL series (NB01–NB10) — Intervention-Conditioned Phenotyping programme
- **`CALON_NOBEL_NB01_Discordance/`** — Discordance Paradox (submission-ready). Statin-treatment
  paradox OR 2.05 (1.61–2.60), p=3.8×10⁻⁹; duration–dose gradient; Non-Identifiability Theorem.
- **`CALON_NOBEL_NB02_Penetrance/`** — LDLR penetrance 3D (in progress). SSS-high OR 1.59 (p=0.028)
  at age 45, attenuating by 55 (Biobank Event-Horizon law).
- **`CALON_NOBEL_NB03_TissueAtlas/`** … **`NB10_SSSv31_Reweighting/`** — design scaffolds for the
  forward programme: tissue atlas, ICD-10 PheWAS, universal residual risk, cholesterol-years
  extended, pre-event NMR, VUS reclassification, sex-specific effects, SSS v3.1 reweighting.
- **`calon-fh-nobel-series-repo/`** — the consolidated NOBEL-series repository / docs.

### Lp(a) sub-programme
- **`Lpa/`** — Lp(a) structural + epidemiological analyses (large CSV footprint — indexed).
- **`Lpa_Manuscript/`** — Lp(a) manuscript drafts, figures, tables.
- **`Lpa_Multilevel/`** — multilevel Lp(a) modelling (large CSV footprint — partly indexed).

### Risk modelling & shared data
- **`CALON_Risk_Models/`** — CALON risk-model variants (full / lite / external).
- **`TUDOR/`** — TUDOR project workspace (code + manuscript material; complements the root pipeline).
- **`SHARED_MASTER_DATA/`** / **`Shared_Data/`** — shared master datasets (Dragon-3, Wales PASS,
  UKB masters — large files indexed).

### Other
- **`perimeopause/`** — separate cardiometabolic sub-study (menopause / reproductive lipids).
- **`Teaching/`** — Cardiff teaching materials and student feedback.
- **`CLAUDE_CONTEXT_BUNDLE/`** — prior context-bundle documents.

---

## How to navigate
1. Start with `docs/WORK_DONE.md` (root pipeline) and `CLAUDE.md`.
2. For a specific project, open `projects/<name>/` — its own README / docs are inside.
3. For a CSV: check `docs/DATA_DICTIONARY.md` (root data) or browse `projects/<name>/`.
4. For a missing giant file: look it up in `docs/LARGE_DATA_INDEX.csv` for its original D: path.

> The NOBEL-series numbers above are preliminary/working values carried from the project context;
> trace any specific claim to the relevant project's results CSV before re-quoting.

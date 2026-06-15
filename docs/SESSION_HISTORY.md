# Session History — the work, discussion by discussion

A chronological record of the Claude Code discussions that built this programme, derived by
going through the session log one by one. Each row is one discussion: its title, the working
directory it ran in, its git branch (where it used an isolated worktree), and last activity.
This is the "work done from scratch" provenance — it shows *when* and *in what order* each
piece of the bundle was created.

> Source: local Claude Code session index (26 sessions), read 2026-05-30. Dates are last-activity.
> Sessions ran in isolated git worktrees under `.claude/worktrees/` off the pipeline repo, or
> directly in the relevant project directory.

---

## A. Core CALON-FH / TUDOR / structural programme

| Date | Discussion | Worked in | Maps to bundle |
|---|---|---|---|
| 2026-05-30 | **TUDOR** (current) | `calon_ukb_pipeline` worktree | root pipeline; this bundle |
| 2026-05-25 | **Develop ASCVD prediction model for FH patients** | `calon_ukb_pipeline` worktree | `projects/CALON_Risk_Models/`, root CALON |
| 2026-05-24 | **INVENTORY_ukb_CALONF** | `calon_ukb_pipeline` worktree | UKB data inventory; `docs/DATA_DICTIONARY.md` |
| 2026-05-13 | **LPA_CLINICAL** | `calon_ukb_pipeline` worktree | `projects/Lpa/`, `projects/Lpa_Manuscript/` |
| 2026-05-13 | **Lpa** | `calon_ukb_pipeline` worktree | `projects/Lpa/`, `projects/Lpa_Multilevel/` |
| 2026-05-13 | **apob_ldl** | `D:\Projects` | ApoB/LDL discordance (NB01 / CALON) |
| 2026-04-17 | **Analyze ASCVD events across cohorts with lipid panels** | `calon_ukb_pipeline` worktree | root CALON, Table1b ASCVD CSVs |
| 2026-04-14 | **General Updates and Improvements** | `calon_ukb_pipeline` worktree | root pipeline maintenance |
| 2026-03-29 | **Fix Unsloth installation and execution issue** | `calon_ukb_pipeline` | environment / tooling |
| 2026-03-28 | **Review Tudor research and R code analysis** | `calon_ukb_pipeline` | `code/R/`, TUDOR modelling |
| 2026-03-26 | **Find CVS prediction GitHub link** | Desktop | repo / publication links |
| 2026-03-26 | **Check FoldX saturation progress** | `calon_ukb_pipeline` | FoldX ddG layer (SSS v3) |
| 2026-03-21 | **Check FoldX AlphaFold3 progress** | `calon_ukb_pipeline` worktree | AlphaFold3 + FoldX structural pipeline |
| 2026-03-01 | **Build ASCVD risk prediction model for FH** | `calon_ukb_pipeline` worktree | CALON ASCVD model origin |
| 2026-03-01 | **ASCVD prediction model for HeFH** | `calon_ukb_pipeline` worktree | CALON ASCVD model origin |

## B. Adjacent clinical / teaching / other work

| Date | Discussion | Worked in | Notes |
|---|---|---|---|
| 2026-05-25 | Organize and QC ultrasound screenshots by patient | `D:\Recycle_Bin_Saved` | clinical imaging admin |
| 2026-05-13 | PERIMENOPAUSE | `calon_ukb_pipeline` worktree | `projects/perimeopause/` sub-study |
| 2026-05-11 | Review ECG and generate detailed cardiology report | `calon_ukb_pipeline` worktree | clinical reporting |
| 2026-05-10 | Learna | `D:\Projects` | (learning / misc) |
| 2026-04-01 | Create lipid medicine presentation modules | `Downloads\MODULES` | `projects/Teaching/` (lipid teaching) |
| 2026-03-31 | Remove curriculum involvement claim from documents | `D:\HE Fellowship` | HE Fellowship admin |
| 2026-03-26 | Prepare verified CSV with complete scores | `D:\RSA` | RSA-PACE study |
| 2026-03-22 | Grade acute internal medicine master student | `MASTER MARKING` | student marking |
| 2026-03-17 | Analyze RSA-PACE cardiac study data | `Downloads\RSA` | RSA-PACE study |
| 2026-03-16 | Medical education professor role setup | `Downloads\RSA` | teaching role setup |
| 2026-03-02 | Fix validation and API rate limit handling | `New folder (3)` | tooling |

---

## Reading the timeline

- **March 2026** — foundations: ASCVD prediction model first built (1 Mar), TUDOR R-code review
  (28 Mar), and the AlphaFold3 + FoldX structural pipeline run for the SSS v3 ddG layer (21–26 Mar).
- **April 2026** — cohort ASCVD-event analysis across Wales/UKB (17 Apr); general pipeline
  hardening (14 Apr).
- **May 2026** — Lp(a) sub-programme and clinical analyses (13 May), perimenopause sub-study
  (13 May), full UKB inventory (24 May), ASCVD-for-FH model development (25 May), and the
  current TUDOR submission-prep + this portable-bundle work (30 May).

The structural work (SSS v3) precedes the clinical validation, which precedes the manuscript
and submission work — matching the dependency order described in `docs/WORK_DONE.md`.

> To re-read any discussion in full, open it in Claude Code (session index) or search transcripts.
> This catalogue captures the metadata (title/dir/branch/date); the per-session detail lives in
> the project folders the sessions produced (`code/`, `projects/<name>/`, `manuscripts/`).

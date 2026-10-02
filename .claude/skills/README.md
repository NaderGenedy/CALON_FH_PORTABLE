# Google DeepMind Science Skills (vendored)

This directory holds the **Google DeepMind Science Skills** collection, built
into the CALON-FH / TUDOR repo so that any Claude Code session opened here can
call them without an install step.

| | |
|---|---|
| Upstream | https://github.com/google-deepmind/science-skills |
| Pinned release | `v1.2.1`, commit `68832757cbbf941c620b71df5756cf6e5cc287b0` |
| Licence | Code Apache-2.0, docs CC-BY-4.0 -- see `LICENSE-google-deepmind-science-skills.txt`; per-database terms in `SKILL_LICENSES.md` |
| Index | `SCIENCE_SKILLS_INDEX.md` (generated) |
| Build report | `../../docs/SCIENCE_SKILLS_BUILD_REPORT.md` (generated) |
| Build script | `../../build_science_skills.py` |

The 17 project-specific skills (tudor-qc, cox-analysis, manuscript-qc, ...)
live in `../../skills/` and are documented in `../../docs/SKILLS.md`.

## Layout

Each skill is a flat directory `.claude/skills/<skill>/` with a `SKILL.md`
(YAML frontmatter + instructions), usually `scripts/` (Python CLIs) and
`references/` (citation .bib, API notes). Claude Code discovers them
automatically from this path. Nothing was modified from upstream.

## Prerequisites

1. **uv** -- every script runs as `uv run scripts/<x>.py ...`; dependencies are
   declared inline (PEP 723) or in a per-skill `pyproject.toml` and are resolved
   on first use. Install with `curl -LsSf https://astral.sh/uv/install.sh | sh`
   (see `uv/SKILL.md`). The three AlphaGenome skills need Python 3.14, which
   uv downloads on demand (`uv python install 3.14`).
2. **Credentials** go in `~/.env`, never in the chat or the repo. Follow
   `credentials/SKILL.md` (quiet `grep -sq` check, hidden `read -s` prompt).

| Variable | Needed by | Required? |
|---|---|---|
| `ALPHAGENOME_API_KEY` | alphagenome_single_variant_analysis, alphagenome_variant_impact_score | yes |
| `NCBI_API_KEY` | clinvar, dbsnp, pubmed, ncbi_sequence_fetch | optional (3 -> 10 req/s) |
| `USER_EMAIL` | protein_sequence_msa, protein_sequence_similarity_search (EBI), openalex polite pool | yes for EBI tools |
| `OPENALEX_API_KEY` | literature_search_openalex | optional |
| `FDA_API_KEY` | openfda_database | optional |

Some skills also write a notice file to `.licenses/<skill>_LICENSE.txt` in the
workspace root the first time they run; that folder is git-ignored.

## Rebuilding / updating

```bash
python build_science_skills.py              # regenerate index + smoke-test (uv --help on every script)
python build_science_skills.py --update     # re-sync from upstream at the pinned tag
python build_science_skills.py --update --tag vX.Y.Z --commit <sha>   # bump the pin
```

The smoke test classifies each script PASS / CRED (aborted only for a missing
API key) / FAIL and writes `docs/SCIENCE_SKILLS_BUILD_REPORT.md`. Build
artefacts (`.venv/`, `uv.lock`, `__pycache__/`) are scrubbed and git-ignored.

## Where they bite in this programme

| Need | Skill(s) |
|---|---|
| LDLR variant pathogenicity labels / star ratings for SSS v3 validation | `clinvar_database`, `dbsnp_database` |
| Carrier allele frequencies, LDLR / APOB / PCSK9 constraint (pLI, LOEUF) | `gnomad_database` |
| Transcript / protein coordinates, VEP consequences for the 16,340-variant map | `ensembl_database`, `uniprot_database` |
| AlphaFold pLDDT / PAE for the FoldX ddG and B&G domain layers | `alphafold_database_fetch_and_analyze`, `pdb_database`, `pymol`, `foldseek_structural_search` |
| AlphaGenome L2 / splice layers (SSS weights 0.05 each) | `alphagenome_variant_impact_score`, `alphagenome_single_variant_analysis`, `alphagenome_atlas_website_links` |
| Conservation (phyloP) and regulatory context for non-coding LDLR variants | `ucsc_conservation_and_tfbs`, `encode_ccres_database`, `jaspar_database`, `unibind_database` |
| Literature for TUDOR / CALON manuscripts and response-to-reviewers | `pubmed_database`, `literature_search_europepmc`, `literature_search_openalex`, `literature_search_biorxiv` |
| Statin / PCSK9i / ezetimibe trial landscape and drug data | `clinical_trials_database`, `chembl_database`, `pubchem_database`, `openfda_database`, `opentargets_database` |
| Tissue expression of LDLR and candidate modifiers | `gtex_database`, `human_protein_atlas_database` |
| Pathway / interaction context | `reactome_database`, `string_database`, `quickgo_database`, `interpro_database`, `embl_ebi_ols` |
| Turning a finished analysis into a new skill | `workflow_skill_creator` |

`predictingthepast` (Latin / Greek epigraphy) ships with the upstream bundle and
is kept for completeness only.

Project rules still apply on top of these tools: UK Biobank participant data
never passes through any of these APIs, and every number they return must be
traced to its source before it reaches a manuscript.

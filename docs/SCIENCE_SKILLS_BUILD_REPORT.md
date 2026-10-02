# Science Skills build report

Generated 2026-10-02 11:42 by `build_science_skills.py`.  
Upstream `v1.2.1` (`68832757cbbf`); 40 skills vendored in `.claude/skills/`.

Each helper script was executed as `uv run [--project <skill>] <script> --help`.
PASS = resolved dependencies and printed usage. CRED = aborted only because an API key
or e-mail is not in `~/.env` (add it per `.claude/skills/credentials/SKILL.md`).
FAIL = anything else; see the last line captured.

| PASS | CRED | FAIL | total |
|---|---|---|---|
| 64 | 1 | 0 | 65 |

| Script | Mode | Status | s | Last line |
|---|---|---|---|---|
| `alphafold_database_fetch_and_analyze/scripts/analyze_pae.py` | script | PASS | 0 | -h, --help  show this help message and exit |
| `alphafold_database_fetch_and_analyze/scripts/analyze_plddt.py` | script | PASS | 0 | -h, --help     show this help message and exit |
| `alphafold_database_fetch_and_analyze/scripts/fetch_structure.py` | script | PASS | 0 | Output directory to save the files (required) |
| `alphagenome_atlas_website_links/scripts/alphagenome_atlas_links.py` | project | PASS | 6 | -h, --help            show this help message and exit |
| `alphagenome_single_variant_analysis/scripts/analyze_ism.py` | project | PASS | 7 | Minimum absolute score threshold for motif extraction. |
| `alphagenome_single_variant_analysis/scripts/generate_ontology_mapping.py` | project | CRED | 3 | RuntimeError: ALPHAGENOME_API_KEY not set. Ensure the `.env` file contains ALPHAGENOME_API_KEY=<key> and use uv run to run this script. |
| `alphagenome_single_variant_analysis/scripts/interpret_splicing.py` | project | PASS | 3 | --window WINDOW       Window size around variant for analysis. |
| `alphagenome_single_variant_analysis/scripts/lookup_gene_info.py` | project | PASS | 1 | $ALPHAGENOME_GTF_PATH or the public GCS URL). |
| `alphagenome_single_variant_analysis/scripts/resolve_ontology_terms.py` | project | PASS | 0 | Path to tissue_ontology_mapping.json. |
| `alphagenome_single_variant_analysis/scripts/visualize_genome_tracks.py` | project | PASS | 7 | $ALPHAGENOME_GTF_PATH or the public GCS URL). |
| `alphagenome_single_variant_analysis/scripts/visualize_variant_effects.py` | project | PASS | 6 | $ALPHAGENOME_GTF_PATH or the public GCS URL). |
| `alphagenome_variant_impact_score/scripts/alphagenome_atlas_avi.py` | project | PASS | 5 | -h, --help            show this help message and exit |
| `chembl_database/scripts/chembl_api.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `clinical_trials_database/scripts/clinical_trials_api.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `clinvar_database/scripts/clinvar_api.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `dbsnp_database/scripts/dbsnp_cli.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `embl_ebi_ols/scripts/get_individual.py` | script | PASS | 0 | --output OUTPUT      Output file path |
| `embl_ebi_ols/scripts/get_ontology.py` | script | PASS | 0 | --output OUTPUT  Output file path |
| `embl_ebi_ols/scripts/get_property.py` | script | PASS | 0 | --output OUTPUT       Output file path |
| `embl_ebi_ols/scripts/get_stats.py` | script | PASS | 0 | --output OUTPUT  Output file path |
| `embl_ebi_ols/scripts/get_term.py` | script | PASS | 0 | --output OUTPUT       Output file path |
| `embl_ebi_ols/scripts/ols_utils.py` | script | PASS | 0 | Installed 1 package in 1ms |
| `embl_ebi_ols/scripts/search_ols.py` | script | PASS | 0 | --output OUTPUT       Output file path |
| `embl_ebi_ols/scripts/suggest_ols.py` | script | PASS | 0 | --output OUTPUT      Output file path |
| `encode_ccres_database/scripts/encode_portal_api.py` | script | PASS | 0 | -h, --help  show this help message and exit |
| `encode_ccres_database/scripts/screen_api.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `ensembl_database/scripts/ensembl_api.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `foldseek_structural_search/scripts/search.py` | script | PASS | 0 | Comma-separated list of databases to search |
| `gnomad_database/scripts/get_gene_constraint.py` | script | PASS | 0 | --output, -o OUTPUT  Output file path. Prints to stdout if not specified. |
| `gnomad_database/scripts/get_variant_frequency.py` | script | PASS | 0 | --output, -o OUTPUT   Output file path. |
| `gnomad_database/scripts/search_variants.py` | script | PASS | 0 | --output, -o OUTPUT   Output file path. |
| `gtex_database/scripts/gtex_cli.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `human_protein_atlas_database/scripts/hpa_cli.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `interpro_database/scripts/interpro_client.py` | script | PASS | 0 | -h, --help     show this help message and exit |
| `jaspar_database/scripts/jaspar_api.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `literature_search_arxiv/scripts/download_paper.py` | script | PASS | 0 | --output OUTPUT      Output file path |
| `literature_search_arxiv/scripts/download_paper_source.py` | script | PASS | 0 | --output OUTPUT  Output file path for the tar.gz file |
| `literature_search_arxiv/scripts/search_arxiv.py` | script | PASS | 0 | Sort order |
| `literature_search_biorxiv/scripts/search_by_dates.py` | script | PASS | 0 | --include_abstracts   Include full abstracts in the JSON output |
| `literature_search_biorxiv/scripts/search_by_doi.py` | script | PASS | 0 | --include_abstracts   Include full abstracts in the JSON output |
| `literature_search_europepmc/scripts/europepmc_api.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `literature_search_openalex/scripts/openalex_cli.py` | script | PASS | 0 | variable. |
| `ncbi_sequence_fetch/scripts/ncbi_fetch.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `openfda_database/scripts/openfda_query.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `opentargets_database/scripts/query_opentargets.py` | script | PASS | 0 | --output OUTPUT       Path to write the JSON output file |
| `pdb_database/scripts/download_coordinate_files.py` | script | PASS | 0 | Directory to save files to |
| `pdb_database/scripts/fetch_pdb_metadata.py` | script | PASS | 0 | --output OUTPUT  File to write the output to |
| `pdb_database/scripts/fetch_schema.py` | script | PASS | 0 | --output OUTPUT       File to write the schema to. |
| `pdb_database/scripts/search_pdb.py` | script | PASS | 0 | --output OUTPUT       File to write the output to |
| `predictingthepast/scripts/preprocess.py` | script | PASS | 0 | Path to save cleaned text. |
| `predictingthepast/scripts/run_inference.py` | script | PASS | 6 | --embedding           Generate a text embedding vector. |
| `predictingthepast/scripts/visualize_results.py` | script | PASS | 0 | --output OUTPUT  Path to write HTML dashboard. Prints to stdout if omitted. |
| `protein_sequence_msa/scripts/msa_align.py` | script | PASS | 0 | --dry-run            Dry run: print payload and exit without submitting job |
| `protein_sequence_similarity_search/scripts/mmseqs2_search.py` | script | PASS | 0 | mgnify/environmental database |
| `protein_sequence_similarity_search/scripts/uniprot_blast.py` | script | PASS | 0 | Comma-separated list of databases to search |
| `pubchem_database/scripts/pubchem_api.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `pubmed_database/scripts/pubmed_api.py` | script | PASS | 0 | Available: search_pubmed, fetch_article_abstracts, find_linked_biological_data, discover_available_links, get_full_text_pmc, verify_medical_spelling, global_dat |
| `quickgo_database/scripts/quickgo_tool.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `reactome_database/scripts/reactome_analysis.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `string_database/scripts/string_cli.py` | script | PASS | 0 | STRING API version (default: 12) |
| `ucsc_conservation_and_tfbs/scripts/get_conservation.py` | script | PASS | 0 | acceleration. |
| `ucsc_conservation_and_tfbs/scripts/get_tfbs.py` | script | PASS | 0 | --genome GENOME       Genome assembly. Defaults to hg38. |
| `ucsc_conservation_and_tfbs/scripts/list_tracks.py` | script | PASS | 0 | --output OUTPUT  Path where the output matches will be saved in JSON format. |
| `unibind_database/scripts/unibind_api.py` | script | PASS | 0 | -h, --help            show this help message and exit |
| `uniprot_database/scripts/uniprot_tools.py` | script | PASS | 0 | -h, --help            show this help message and exit |

"""
╔══════════════════════════════════════════════════════════════════════════╗
║  TUDOR PIPELINE RECONSTRUCTION GUIDE                                    ║
║  Clinical Risk Prediction in Familial Hypercholesterolaemia             ║
║  Dr Nader Genedy — UHW Cardiff                                         ║
║                                                                        ║
║  STATUS: IN PROGRESS                                                    ║
║  R scripts done, manuscript v4 drafted, needs final submission prep.    ║
╚══════════════════════════════════════════════════════════════════════════╝

CONCEPT:
  TUDOR is a clinical risk prediction model for FH that combines:
    - Genetic data (LDLR/APOB/PCSK9 variant type)
    - Lipid profile (LDL-C, HDL-C, TG, Lp(a))
    - Clinical features (age, sex, BMI, smoking, BP)
    - Imaging (carotid IMT, coronary calcium)
  Validated in Wales FH Registry + UK Biobank.

TARGET JOURNAL: European Heart Journal or Lancet Digital Health

REPO: github.com/NaderGenedy/calon-ukb-pipeline (public)
LOCAL: C:/Users/nader/Downloads/calon_ukb_pipeline/
"""
import os, sys
from pathlib import Path
sys.stdout.reconfigure(encoding='utf-8')

BASE = Path('C:/Users/nader/Downloads/calon_ukb_pipeline')

# ═══════════════════════════════════════════════════════════════════════
# FILE MAP
# ═══════════════════════════════════════════════════════════════════════
FILES = {
    # === R Scripts (core pipeline) ===
    'extract':          BASE / '00_CALON_extract_ukbrap.sh',
    'build_cohort':     BASE / '01_CALON_build_cohort.R',
    'validate':         BASE / '02_CALON_validate.R',
    'cox_model':        BASE / '03_CALON_cox_model.R',
    'figures':          BASE / '04_CALON_figures.R',
    'nmr':              BASE / '05_CALON_nmr.R',
    'mri':              BASE / '06_CALON_mri.R',

    # === TUDOR-specific R scripts ===
    'tudor_main':       BASE / 'TUDOR_LANCET_COMPLETE.R',
    'tudor_loco':       BASE / 'TUDOR_LOCO_CV_FULL.R',
    'tudor_figures':    BASE / 'TUDOR_MANUSCRIPT_FIGURES.R',
    'tudor_tables':     BASE / 'TUDOR_tables_nature.R',
    'tudor_ukb':        BASE / 'TUDOR_UKB_LIPID_CLINIC.R',
    'tudor_iterate':    BASE / 'TUDOR_ITERATE_UKB.R',
    'tudor_h2h':        BASE / 'TUDOR_HEAD_TO_HEAD.R',
    'tudor_fixed':      BASE / 'TUDOR_FIXED_WEIGHTS.R',
    'tudor_resubmit':   BASE / 'TUDOR_resubmission_package.R',
    'tudor_supplement':  BASE / 'tudor_supplement_new_files.R',

    # === TUDOR manuscripts ===
    'manuscript_v4':    BASE / 'TUDOR_Manuscript_v4.docx',
    'manuscript_v3':    BASE / 'TUDOR_Manuscript_v3.docx',
    'cover_letter':     BASE / 'TUDOR_CoverLetter.docx',
    'cover_resubmit':   BASE / 'TUDOR_CoverLetter_Resubmission.docx',
    'response':         BASE / 'TUDOR_ResponseToReviewers.docx',
    'highlights':       BASE / 'TUDOR_Highlights.docx',
    'tables_doc':       BASE / 'TUDOR_Tables.docx',

    # === Data files ===
    'dragon3':          BASE / 'DRAGON_3.csv',
    'pass_drugs':       BASE / 'PASS_drug_codes.xls',
    'loco_results':     BASE / 'LOCO_CV_results.csv',
    'tudor_model':      BASE / 'TUDOR_v2_model.rds',
    'reproduced_auc':   BASE / 'TUDOR_reproduced_AUC.csv',
    'reproduced_coef':  BASE / 'TUDOR_reproduced_coefficients.csv',

    # === Table CSVs ===
    'table1':           BASE / 'Table1_baseline_characteristics.csv',
    'table2a':          BASE / 'Table2a_gene_distribution.csv',
    'table3':           BASE / 'Table3_index_vs_cascade.csv',
    'table4':           BASE / 'Table4_FH_vs_NonFH_matched.csv',
    'table5':           BASE / 'Table5_biexternal_validation.csv',

    # === Figures ===
    'fig1_pdf':         BASE / 'Figure1_cohort_overview.pdf',
    'fig2_pdf':         BASE / 'Figure2_lipid_profiles.pdf',
    'fig3_pdf':         BASE / 'Figure3_lpa_distribution.pdf',
    'fig4_pdf':         BASE / 'Figure4_apob_vs_ldl.pdf',
    'fig5_pdf':         BASE / 'Figure5_forest_plot.pdf',

    # === Python analysis scripts (numbered) ===
    # 17-73 Python scripts for various analyses
    'af3_manuscript':   BASE / 'AF3_manuscript_nobel.py',

    # === Shared raw data ===
    'shared_dragon':    Path('D:/Projects/Shared_Data/DRAGON_3.csv'),
    'shared_wales':     Path('D:/Projects/Shared_Data/combined_soretd_ldl.csv'),

    # === Website ===
    'atlas_html':       BASE / 'website/index.html',
    'atlas_data':       BASE / 'website/website_data.json',
}

# ═══════════════════════════════════════════════════════════════════════
# HOW TO RECONSTRUCT
# ═══════════════════════════════════════════════════════════════════════
print("""
╔══════════════════════════════════════════════════════════════════════════╗
║  TUDOR RECONSTRUCTION STEPS                                             ║
╚══════════════════════════════════════════════════════════════════════════╝

STEP 1: R Environment Setup
  install.packages(c("tidyverse", "survival", "survminer", "pROC",
                      "gtsummary", "patchwork", "ggforestplot"))

STEP 2: Run Core Pipeline (in R)
  source("TUDOR_LANCET_COMPLETE.R")
  # This builds the full model, runs LOCO-CV, generates AUC

STEP 3: Generate Figures (in R)
  source("TUDOR_MANUSCRIPT_FIGURES.R")
  # Generates Figures 1-8 as PDF + PNG

STEP 4: Generate Tables (in R)
  source("TUDOR_tables_nature.R")
  # Generates formatted tables

STEP 5: UKB External Validation (in R)
  source("TUDOR_UKB_LIPID_CLINIC.R")
  source("TUDOR_ITERATE_UKB.R")

STEP 6: Head-to-Head Comparison
  source("TUDOR_HEAD_TO_HEAD.R")

STEP 7: Manuscript
  # Current: TUDOR_Manuscript_v4.docx
  # Cover letter: TUDOR_CoverLetter.docx
  # Response to reviewers: TUDOR_ResponseToReviewers.docx

REMAINING WORK:
  - Final submission formatting
  - Supplementary materials assembly
  - Address any remaining reviewer comments
  - Format for target journal

NOTE: TUDOR is primarily an R project.
      Python scripts (17-73) are for AF3/FoldX/SSS analyses.
      The R scripts are the primary pipeline.
""")

print("FILE STATUS:")
for name, path in FILES.items():
    exists = "EXISTS" if path.exists() else "MISSING"
    size = f"({path.stat().st_size/1024:.0f} KB)" if path.exists() else ""
    print(f"  {exists:8s}  {name:25s}  {size:>12s}  {path}")

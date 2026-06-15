#!/bin/bash
################################################################################
# PAPER 1 (ApoB/LDL discordance) — Olink Explore 3072 full-panel extraction
#
# PURPOSE: Extract NPX values for all 2,923 proteins in the Olink Explore 3072
#          panel for the ~54,000 UKB participants with proteomic data.
#          Used for PWAS validation of the §3.17 PRS-proxy findings (Path A).
#
# UKB STRUCTURE (Category 1839 — Proteomics):
#   /mnt/project/Bulk/Proteomics/
#     olink_data.tsv           ← long-format NPX (eid × protein × instance)
#     coding143.tsv            ← protein_id → UniProt + gene + assay metadata
#     coding3.tsv              ← Olink panel codes
#
# RUN ON: UKB-RAP JupyterLab Bash terminal (NOT local).
# OUTPUT: olink_npx_wide.tsv  (eid × 2923 protein columns) — ~3 GB
#         olink_protein_meta.tsv  (protein_id, uniprot, gene, panel)
#         olink_instance_dates.tsv (eid, instance_0_date, instance_2_date)
#
# DOWNLOAD TARGET: D:/Projects/CALON_AlphaFold_Rebuild/data/olink/
#
# Author: Dr Nader Genedy
# Date:   2026-05-05
# Project: CALON-FH Paper 1 — ApoB/LDL discordance, NEJM revision
################################################################################

set -euo pipefail

PROJECT=$(dx env --bash | grep DX_PROJECT_CONTEXT_ID | cut -d= -f2 | tr -d '"')
OUTPUT_DIR="./paper1_olink"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "============================================================"
echo "  PAPER 1: OLINK EXPLORE 3072 — FULL PANEL (2,923 proteins)"
echo "  ~54,000 UKB participants × 2,923 proteins"
echo "============================================================"

# ──────────────────────────────────────────────────────────────────
# STEP 1: Locate Olink files on RAP
# ──────────────────────────────────────────────────────────────────
PROTEOMICS_DIR="/mnt/project/Bulk/Proteomics"
echo ""
echo "Step 1: Inventory of /mnt/project/Bulk/Proteomics/"
ls -la "$PROTEOMICS_DIR" 2>/dev/null | head -20

# Olink long-format NPX: eid \t plate \t protein \t instance \t value
OLINK_LONG="$PROTEOMICS_DIR/olink_data.dat"
if [ ! -f "$OLINK_LONG" ]; then
    echo "WARNING: $OLINK_LONG not found. Trying alternative paths..."
    OLINK_LONG=$(find /mnt/project/Bulk/Proteomics -name "olink_data*" 2>/dev/null | head -1)
fi
echo "Olink long-format file: $OLINK_LONG"

# Protein metadata (id → UniProt + gene)
META_FILE=$(find /mnt/project/Bulk/Proteomics -name "coding143*" -o -name "olink_*meta*" 2>/dev/null | head -1)
echo "Protein metadata file: $META_FILE"

# ──────────────────────────────────────────────────────────────────
# STEP 2: Pivot long → wide via Swiss Army Knife (handles 2.5M rows)
# ──────────────────────────────────────────────────────────────────
echo ""
echo "Step 2: Pivoting long → wide format with Swiss Army Knife..."

cat > pivot_olink.py <<'PYEOF'
import pandas as pd
import sys
print("Loading Olink long-format...", flush=True)
long = pd.read_csv("/mnt/project/Bulk/Proteomics/olink_data.dat",
                    sep="\t", low_memory=False)
print(f"  {len(long):,} rows × {long.shape[1]} cols", flush=True)
print(f"  Columns: {list(long.columns)}", flush=True)
print(f"  Unique participants: {long['eid'].nunique():,}", flush=True)
print(f"  Unique proteins:     {long['protein_id'].nunique():,}", flush=True)
print(f"  Instances:           {sorted(long['ins_index'].unique())}", flush=True)

# Pivot — instance 0 (baseline) primary
print("\nPivoting baseline (instance 0)...", flush=True)
wide_i0 = long[long['ins_index'] == 0].pivot_table(
    index='eid', columns='protein_id', values='result', aggfunc='first'
)
wide_i0.columns = [f"olink_{c}_i0" for c in wide_i0.columns]
wide_i0 = wide_i0.reset_index()
print(f"  Wide i0: {wide_i0.shape}", flush=True)
wide_i0.to_csv("olink_npx_wide_i0.tsv", sep="\t", index=False)
print("  Saved: olink_npx_wide_i0.tsv", flush=True)

# Imaging visit (instance 2) for sensitivity
print("\nPivoting imaging visit (instance 2)...", flush=True)
wide_i2 = long[long['ins_index'] == 2].pivot_table(
    index='eid', columns='protein_id', values='result', aggfunc='first'
)
wide_i2.columns = [f"olink_{c}_i2" for c in wide_i2.columns]
wide_i2 = wide_i2.reset_index()
print(f"  Wide i2: {wide_i2.shape}", flush=True)
wide_i2.to_csv("olink_npx_wide_i2.tsv", sep="\t", index=False)
print("  Saved: olink_npx_wide_i2.tsv", flush=True)
PYEOF

python3 pivot_olink.py 2>&1 | tee pivot_olink.log

# ──────────────────────────────────────────────────────────────────
# STEP 3: Extract protein metadata (id → gene name → UniProt)
# ──────────────────────────────────────────────────────────────────
echo ""
echo "Step 3: Protein metadata"
if [ -f "$META_FILE" ]; then
    cp "$META_FILE" olink_protein_meta.tsv
    echo "  Copied: olink_protein_meta.tsv ($(wc -l < olink_protein_meta.tsv) rows)"
else
    echo "  WARNING: metadata file not found; will use protein_id only"
fi

# ──────────────────────────────────────────────────────────────────
# STEP 4: Merge with ApoB/LDL/covariates from analysis-ready data
# ──────────────────────────────────────────────────────────────────
echo ""
echo "Step 4: Merging Olink with ApoB/LDL covariates"
cat > merge_olink_covars.py <<'PYEOF'
import pandas as pd
print("Loading Olink baseline...", flush=True)
olink = pd.read_csv("olink_npx_wide_i0.tsv", sep="\t")

# Pull ApoB (p30640), LDL (p30780), age (p21022), sex (p31), BMI (p21001),
# statin status (will come from medication codes)
print("Loading ApoB/LDL covariates...", flush=True)
covars = pd.read_csv("/mnt/project/calon_data/calon_batch2_lipids.csv")
covars = covars[['participant.eid', 'participant.p30640_i0',
                  'participant.p30780_i0']].copy()
covars.columns = ['eid', 'apob', 'ldl']
covars['apob_ldl_ratio'] = covars['apob'] / covars['ldl']
covars['log_apob_ldl'] = covars['apob_ldl_ratio'].apply(
    lambda x: pd.NA if pd.isna(x) or x <= 0 else __import__('math').log(x))

demo = pd.read_csv("/mnt/project/calon_data/calon_batch1_demographics.csv")
demo = demo[['participant.eid', 'participant.p31',
              'participant.p21022', 'participant.p21001_i0']].copy()
demo.columns = ['eid', 'sex', 'age', 'bmi']

merged = olink.merge(covars, on='eid', how='inner').merge(demo, on='eid', how='inner')
merged.to_csv("olink_apob_ldl_merged.tsv", sep="\t", index=False)
print(f"  Merged: {merged.shape[0]:,} participants × {merged.shape[1]:,} cols", flush=True)
print(f"  With ApoB/LDL ratio non-NA: {merged['log_apob_ldl'].notna().sum():,}", flush=True)
PYEOF

python3 merge_olink_covars.py 2>&1 | tee merge_olink.log

# ──────────────────────────────────────────────────────────────────
# STEP 5: Upload results to project root
# ──────────────────────────────────────────────────────────────────
echo ""
echo "Step 5: Uploading to RAP project root"
for f in olink_npx_wide_i0.tsv olink_npx_wide_i2.tsv olink_protein_meta.tsv \
         olink_apob_ldl_merged.tsv pivot_olink.log merge_olink.log; do
  if [ -f "$f" ]; then
    dx upload "$f" --destination /paper1_olink/ --brief
    echo "  ✓ Uploaded: $f"
  fi
done

echo ""
echo "============================================================"
echo "  EXTRACTION COMPLETE"
echo "  Download these files to D:/Projects/CALON_AlphaFold_Rebuild/data/olink/"
echo "  Then run PAPER1_olink_pwas.ipynb on Colab"
echo "============================================================"

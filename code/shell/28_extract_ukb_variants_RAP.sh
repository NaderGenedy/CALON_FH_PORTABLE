#!/bin/bash
################################################################################
#
#  28_extract_ukb_variants_RAP.sh
#  Extract exact FH variants from UK Biobank WES/WGS on RAP
#
#  PURPOSE: Get variant-level data (cDNA change, protein change, consequence,
#           CADD score) for all LDLR/APOB/PCSK9 variant carriers in UKB
#
#  PLATFORM: UK Biobank Research Analysis Platform (RAP / DNAnexus)
#
#  INPUT:  UKB WES population-level VCF (pVCF) or VEP annotations
#  OUTPUT: ukb_fh_variants.csv (eid, gene, variant_cdna, variant_protein,
#          consequence, CADD, SIFT, PolyPhen, gnomAD_AF)
#
#  Author: Dr Nader Genedy
#  Date:   March 2026
#
#  HOW TO RUN:
#    1. Log into RAP (https://ukbiobank.dnanexus.com)
#    2. Open a JupyterLab or Cloud Workstation session
#    3. Upload this script
#    4. chmod +x 28_extract_ukb_variants_RAP.sh && bash 28_extract_ukb_variants_RAP.sh
#
#  ALTERNATIVELY: Use the Python version below (28_extract_ukb_variants_RAP.py)
#                 which can be run in a RAP JupyterLab notebook
#
################################################################################

echo "============================================================"
echo "  UKB FH VARIANT EXTRACTION"
echo "  Extracting LDLR/APOB/PCSK9 variants from WES data"
echo "============================================================"

# ──────────────────────────────────────────────────────────────
# OPTION 1: Using the UKB RAP Table Exporter (Recommended)
# ──────────────────────────────────────────────────────────────

# The UKB Exome OQFE dataset has pre-annotated variant calls
# Field 23148 = Exome sequencing variant calls
# These are stored as pVCF files on RAP

# Step 1: Find the pVCF files for chromosomes containing our genes
# LDLR = chromosome 19 (19p13.2)
# APOB = chromosome 2 (2p24.1)
# PCSK9 = chromosome 1 (1p32.3)

echo ""
echo "Step 1: Locating pVCF files for target chromosomes..."

# RAP file paths for UKB WES pVCF (200K exomes, OQFE pipeline)
# These paths may vary — check your RAP project structure
PVCF_DIR="/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release"

# Check if files exist
for CHR in 1 2 19; do
    PVCF_FILE="${PVCF_DIR}/ukb23148_c${CHR}_b*.vcf.gz"
    echo "  Chromosome ${CHR}: checking ${PVCF_FILE}"
    ls ${PVCF_FILE} 2>/dev/null | head -1
done

# ──────────────────────────────────────────────────────────────
# OPTION 2: Using bcftools to extract variants (if pVCF available)
# ──────────────────────────────────────────────────────────────

echo ""
echo "Step 2: Extracting variants with bcftools..."

# Define gene regions (GRCh38 coordinates)
# LDLR: chr19:11,089,362-11,133,820 (GRCh38)
# APOB: chr2:21,001,429-21,044,073 (GRCh38)
# PCSK9: chr1:55,039,477-55,064,852 (GRCh38)

LDLR_REGION="19:11089362-11133820"
APOB_REGION="2:21001429-21044073"
PCSK9_REGION="1:55039477-55064852"

OUTPUT_DIR="$(pwd)/ukb_fh_variants"
mkdir -p ${OUTPUT_DIR}

# Extract LDLR variants
echo "  Extracting LDLR variants (chr19)..."
if command -v bcftools &> /dev/null; then
    for PVCF in ${PVCF_DIR}/ukb23148_c19_b*.vcf.gz; do
        if [ -f "$PVCF" ]; then
            bcftools view -r ${LDLR_REGION} \
                --min-ac 1 \
                -Oz -o ${OUTPUT_DIR}/ldlr_variants.vcf.gz \
                ${PVCF}
            echo "    LDLR variants extracted"
            break
        fi
    done

    # Extract APOB variants
    echo "  Extracting APOB variants (chr2)..."
    for PVCF in ${PVCF_DIR}/ukb23148_c2_b*.vcf.gz; do
        if [ -f "$PVCF" ]; then
            bcftools view -r ${APOB_REGION} \
                --min-ac 1 \
                -Oz -o ${OUTPUT_DIR}/apob_variants.vcf.gz \
                ${PVCF}
            echo "    APOB variants extracted"
            break
        fi
    done

    # Extract PCSK9 variants
    echo "  Extracting PCSK9 variants (chr1)..."
    for PVCF in ${PVCF_DIR}/ukb23148_c1_b*.vcf.gz; do
        if [ -f "$PVCF" ]; then
            bcftools view -r ${PCSK9_REGION} \
                --min-ac 1 \
                -Oz -o ${OUTPUT_DIR}/pcsk9_variants.vcf.gz \
                ${PVCF}
            echo "    PCSK9 variants extracted"
            break
        fi
    done
else
    echo "  bcftools not found — use OPTION 3 (Python/pysam) or OPTION 4 (dx extract)"
fi

echo ""
echo "Step 3: Run VEP annotation or use pre-annotated data"
echo "  See Python script below for detailed annotation extraction"
echo ""
echo "============================================================"
echo "  DONE — See 28_extract_ukb_variants_RAP.py for full pipeline"
echo "============================================================"

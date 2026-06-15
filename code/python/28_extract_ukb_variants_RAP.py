#!/usr/bin/env python3
"""
28_extract_ukb_variants_RAP.py
Extract exact FH variants from UK Biobank WES on RAP

PURPOSE: Get variant-level data for all LDLR/APOB/PCSK9 carriers in UKB
         so we can map individual SSS scores for external validation

PLATFORM: UK Biobank Research Analysis Platform (RAP)
          Run in JupyterLab or Cloud Workstation

OUTPUT: ukb_fh_variants.csv with columns:
        eid, gene, chrom, pos, ref, alt, variant_cdna, variant_protein,
        consequence, CADD_phred, SIFT, PolyPhen, gnomAD_AF, is_lof,
        is_missense, is_splice, variant_id_formatted

HOW TO RUN:
  1. Log into RAP (https://ukbiobank.dnanexus.com)
  2. Open JupyterLab (Spark cluster or single-node)
  3. Upload this script to your RAP workspace
  4. python3 28_extract_ukb_variants_RAP.py

  OR copy the cells below into a Jupyter notebook

Author: Dr Nader Genedy
Date:   March 2026
"""

import os
import sys

# ═══════════════════════════════════════════════════════════════
# METHOD 1: Using dxpy + pysam (pVCF approach)
# Best for: RAP Cloud Workstation with bcftools/pysam installed
# ═══════════════════════════════════════════════════════════════

def method1_pvcf_extraction():
    """Extract variants from UKB pVCF files using pysam"""
    try:
        import pysam
    except ImportError:
        print("pysam not installed. Run: pip install pysam")
        print("Falling back to Method 2...")
        return None

    # Gene regions (GRCh38)
    GENES = {
        'LDLR':  {'chrom': 'chr19', 'start': 11089362, 'end': 11133820},
        'APOB':  {'chrom': 'chr2',  'start': 21001429, 'end': 21044073},
        'PCSK9': {'chrom': 'chr1',  'start': 55039477, 'end': 55064852},
    }

    # UKB pVCF paths on RAP (adjust based on your project)
    PVCF_PATHS = {
        'chr1':  '/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release/ukb23148_c1_b0_v1.vcf.gz',
        'chr2':  '/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release/ukb23148_c2_b0_v1.vcf.gz',
        'chr19': '/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release/ukb23148_c19_b0_v1.vcf.gz',
    }

    results = []

    for gene, region in GENES.items():
        chrom_key = region['chrom'].replace('chr', 'chr')
        pvcf_path = PVCF_PATHS.get(chrom_key)

        if not pvcf_path or not os.path.exists(pvcf_path):
            print(f"  pVCF not found for {gene} ({chrom_key}). Trying alternative paths...")
            # Try without 'chr' prefix
            alt_paths = [
                f"/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release/ukb23148_c{region['chrom'].replace('chr','')}_b0_v1.vcf.gz",
                f"/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/ukb23148_c{region['chrom'].replace('chr','')}_b0_v1.bed",
            ]
            for ap in alt_paths:
                if os.path.exists(ap):
                    pvcf_path = ap
                    break

        if not pvcf_path or not os.path.exists(pvcf_path):
            print(f"  WARNING: Cannot find pVCF for {gene}. Skipping.")
            continue

        print(f"  Extracting {gene} variants from {pvcf_path}...")
        vcf = pysam.VariantFile(pvcf_path)

        # Fetch variants in the gene region
        chrom_name = region['chrom'].replace('chr', '')  # pVCF may use '19' not 'chr19'
        try:
            records = vcf.fetch(chrom_name, region['start'], region['end'])
        except ValueError:
            # Try with 'chr' prefix
            try:
                records = vcf.fetch(region['chrom'], region['start'], region['end'])
            except ValueError:
                print(f"  Cannot fetch region for {gene}. Check chromosome naming.")
                continue

        for record in records:
            chrom = str(record.chrom)
            pos = record.pos
            ref = record.ref
            for alt in record.alts:
                # Get VEP annotation from INFO field if available
                vep = record.info.get('vep', record.info.get('CSQ', None))

                # Find carriers (samples with non-ref genotype)
                for sample_name in record.samples:
                    sample = record.samples[sample_name]
                    gt = sample.get('GT', (None, None))
                    if gt and any(a is not None and a > 0 for a in gt):
                        results.append({
                            'eid': sample_name,
                            'gene': gene,
                            'chrom': chrom,
                            'pos': pos,
                            'ref': ref,
                            'alt': alt,
                            'vep_annotation': str(vep)[:200] if vep else '',
                        })

        vcf.close()
        print(f"    Found {sum(1 for r in results if r['gene']==gene)} carrier-variant pairs for {gene}")

    return results


# ═══════════════════════════════════════════════════════════════
# METHOD 2: Using dx extract_dataset (Table Exporter)
# Best for: When pVCF access is difficult
# ═══════════════════════════════════════════════════════════════

def method2_table_exporter():
    """Use UKB RAP Table Exporter to extract variant data"""

    # This generates a dx command to extract variant-level data
    # Run these commands in the RAP terminal

    commands = """
# ═══════════════════════════════════════════════════════════════
# RAP TERMINAL COMMANDS — Copy and run these on RAP
# ═══════════════════════════════════════════════════════════════

# Step 1: Find your dataset ID
dx ls -l "Bulk/Exome sequences/"

# Step 2: Extract LDLR region variants (chr19:11089362-11133820)
dx extract_dataset \\
    "app-Swiss Army Knife" \\
    --brief \\
    -iimage="ghcr.io/dnanexus/swiss-army-knife:v0.0.1" \\
    -icmd="bcftools view -r 19:11089362-11133820 --min-ac 1 \\
           /mnt/project/Bulk/Exome\\ sequences/Population\\ level\\ exome\\ OQFE\\ variants,\\ pVCF\\ format\\ -\\ final\\ release/ukb23148_c19_b0_v1.vcf.gz \\
           -Oz -o ldlr_variants.vcf.gz && \\
           bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CSQ\\n' ldlr_variants.vcf.gz > ldlr_annotations.tsv && \\
           bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\n]' -i 'GT=\"alt\"' ldlr_variants.vcf.gz > ldlr_carriers.tsv"

# Step 3: Extract APOB region variants (chr2:21001429-21044073)
dx extract_dataset \\
    "app-Swiss Army Knife" \\
    --brief \\
    -iimage="ghcr.io/dnanexus/swiss-army-knife:v0.0.1" \\
    -icmd="bcftools view -r 2:21001429-21044073 --min-ac 1 \\
           /mnt/project/Bulk/Exome\\ sequences/Population\\ level\\ exome\\ OQFE\\ variants,\\ pVCF\\ format\\ -\\ final\\ release/ukb23148_c2_b0_v1.vcf.gz \\
           -Oz -o apob_variants.vcf.gz && \\
           bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\n]' -i 'GT=\"alt\"' apob_variants.vcf.gz > apob_carriers.tsv"

# Step 4: Extract PCSK9 region variants (chr1:55039477-55064852)
dx extract_dataset \\
    "app-Swiss Army Knife" \\
    --brief \\
    -iimage="ghcr.io/dnanexus/swiss-army-knife:v0.0.1" \\
    -icmd="bcftools view -r 1:55039477-55064852 --min-ac 1 \\
           /mnt/project/Bulk/Exome\\ sequences/Population\\ level\\ exome\\ OQFE\\ variants,\\ pVCF\\ format\\ -\\ final\\ release/ukb23148_c1_b0_v1.vcf.gz \\
           -Oz -o pcsk9_variants.vcf.gz && \\
           bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\n]' -i 'GT=\"alt\"' pcsk9_variants.vcf.gz > pcsk9_carriers.tsv"

# Step 5: Download the carrier files
dx download ldlr_carriers.tsv
dx download apob_carriers.tsv
dx download pcsk9_carriers.tsv
"""
    print(commands)
    return commands


# ═══════════════════════════════════════════════════════════════
# METHOD 3: Using Hail on Spark (for large-scale extraction)
# Best for: RAP Spark JupyterLab environment
# ═══════════════════════════════════════════════════════════════

def method3_hail_spark():
    """Extract variants using Hail on RAP Spark cluster"""

    hail_code = '''
# ═══════════════════════════════════════════════════════════════
# HAIL EXTRACTION — Run in RAP Spark JupyterLab
# ═══════════════════════════════════════════════════════════════

import hail as hl
hl.init(spark_conf={
    'spark.hadoop.io.compression.codecs': 'org.apache.hadoop.io.compress.DefaultCodec,is.hail.io.compress.BGzipCodec,org.apache.hadoop.io.compress.GzipCodec'
})

# Load the UKB WES data (OQFE pipeline)
# Path may vary — check your RAP project
mt = hl.read_matrix_table(
    'dnax://database-XXXX/ukb23148_c{1,2,19}_b0_v1.mt'
)

# Alternative: load from pVCF
# mt = hl.import_vcf(
#     '/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release/ukb23148_c{1,2,19}_b0_v1.vcf.gz',
#     reference_genome='GRCh38',
#     force_bgz=True
# )

# Define gene intervals (GRCh38)
intervals = [
    hl.parse_locus_interval('chr19:11089362-11133820', reference_genome='GRCh38'),  # LDLR
    hl.parse_locus_interval('chr2:21001429-21044073', reference_genome='GRCh38'),   # APOB
    hl.parse_locus_interval('chr1:55039477-55064852', reference_genome='GRCh38'),   # PCSK9
]

# Filter to gene regions
mt = hl.filter_intervals(mt, intervals)
print(f"Variants in gene regions: {mt.count_rows()}")

# Filter to carriers only (at least one alt allele)
mt = mt.filter_entries(mt.GT.is_non_ref())

# Annotate with VEP (if not already annotated)
# mt = hl.vep(mt, config='vep_config.json')

# Extract carrier information
entries = mt.entries()
entries = entries.select(
    eid=entries.s,
    chrom=entries.locus.contig,
    pos=entries.locus.position,
    ref=entries.alleles[0],
    alt=entries.alleles[1],
    gt=entries.GT,
    gq=entries.GQ,
    dp=entries.DP,
)

# Add gene annotation based on position
entries = entries.annotate(
    gene=hl.case()
        .when((entries.chrom == 'chr19') & (entries.pos >= 11089362) & (entries.pos <= 11133820), 'LDLR')
        .when((entries.chrom == 'chr2') & (entries.pos >= 21001429) & (entries.pos <= 21044073), 'APOB')
        .when((entries.chrom == 'chr1') & (entries.pos >= 55039477) & (entries.pos <= 55064852), 'PCSK9')
        .default('Unknown')
)

# Export to TSV
entries.export('ukb_fh_variant_carriers.tsv')
print("Exported carrier data to ukb_fh_variant_carriers.tsv")

# Count summary
for gene in ['LDLR', 'APOB', 'PCSK9']:
    n = entries.filter(entries.gene == gene).count()
    print(f"  {gene}: {n} carrier-variant pairs")
'''
    print(hail_code)
    return hail_code


# ═══════════════════════════════════════════════════════════════
# METHOD 4: Using the UKB Variant Browser / VEP pre-annotations
# Best for: When you already have VEP-annotated variant lists
# ═══════════════════════════════════════════════════════════════

def method4_vep_preannot():
    """Use UKB's pre-computed VEP annotations"""

    code = '''
# ═══════════════════════════════════════════════════════════════
# VEP PRE-ANNOTATIONS — Using UKB's provided VEP data
# ═══════════════════════════════════════════════════════════════

# UKB provides pre-computed VEP annotations for the WES data
# These are typically in the "Variant annotations" bulk data folder

import pandas as pd
import os

# Path to VEP annotations on RAP
VEP_DIR = "/mnt/project/Bulk/Exome sequences/Exome OQFE variant annotations"

# Check available files
if os.path.exists(VEP_DIR):
    print("VEP annotation files:")
    for f in os.listdir(VEP_DIR):
        print(f"  {f}")

# Alternative: Use the helper files that list consequences per variant
# Field 23149 = Exome OQFE variant consequences

# Load your existing FH patient list
fh_eids = set()
with open('calon_ukb_analysis_ready.csv') as f:
    import csv
    reader = csv.DictReader(f)
    for row in reader:
        fh_eids.add(row.get('eid', '').strip())

print(f"FH patients to match: {len(fh_eids)}")

# If you have the carriers TSV from Method 2:
results = []
for gene_file, gene_name in [
    ('ldlr_carriers.tsv', 'LDLR'),
    ('apob_carriers.tsv', 'APOB'),
    ('pcsk9_carriers.tsv', 'PCSK9')
]:
    if os.path.exists(gene_file):
        with open(gene_file) as f:
            for line in f:
                parts = line.strip().split("\\t")
                if len(parts) >= 6:
                    chrom, pos, ref, alt, sample_eid, gt = parts[:6]
                    if sample_eid in fh_eids:
                        results.append({
                            'eid': sample_eid,
                            'gene': gene_name,
                            'chrom': chrom,
                            'pos': int(pos),
                            'ref': ref,
                            'alt': alt,
                            'gt': gt
                        })

print(f"Matched carriers: {len(results)}")

# Save results
df = pd.DataFrame(results)
df.to_csv('ukb_fh_variants_matched.csv', index=False)
print("Saved to ukb_fh_variants_matched.csv")
'''
    print(code)
    return code


# ═══════════════════════════════════════════════════════════════
# METHOD 5: SIMPLEST — Using dx extract_dataset with SQL-like query
# Best for: Quick extraction without complex bioinformatics tools
# ═══════════════════════════════════════════════════════════════

def method5_dx_extract():
    """Simplest method using RAP's built-in data extraction"""

    code = """
# ═══════════════════════════════════════════════════════════════
# SIMPLEST METHOD — Using RAP built-in tools
# Run these in RAP JupyterLab terminal
# ═══════════════════════════════════════════════════════════════

# Step 1: Find your project's exome dataset
dx describe "Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release"

# Step 2: Use Swiss Army Knife to extract and process
# This creates a single job that:
#   a) Extracts variants in LDLR/APOB/PCSK9 regions
#   b) Filters to carriers only
#   c) Annotates with basic consequence types
#   d) Outputs a clean CSV

cat << 'SCRIPT' > extract_fh_variants.sh
#!/bin/bash
set -e

echo "=== Extracting FH variants from UKB WES ==="

# Install bcftools if not present
apt-get update -qq && apt-get install -y -qq bcftools tabix 2>/dev/null || true

# Define regions
LDLR="19:11089362-11133820"
APOB="2:21001429-21044073"
PCSK9="1:55039477-55064852"

PVCF_BASE="/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release"

OUTPUT="ukb_fh_all_carriers.tsv"
echo -e "CHROM\\tPOS\\tREF\\tALT\\tEID\\tGT\\tGENE" > $OUTPUT

# Extract LDLR
echo "Extracting LDLR..."
for VCF in "${PVCF_BASE}"/ukb23148_c19_b*_v1.vcf.gz; do
    if [ -f "$VCF" ]; then
        bcftools view -r $LDLR "$VCF" 2>/dev/null | \\
        bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\tLDLR\\n]' -i 'GT="alt"' >> $OUTPUT 2>/dev/null || true
    fi
done
LDLR_N=$(grep LDLR $OUTPUT | wc -l)
echo "  LDLR carriers: $LDLR_N"

# Extract APOB
echo "Extracting APOB..."
for VCF in "${PVCF_BASE}"/ukb23148_c2_b*_v1.vcf.gz; do
    if [ -f "$VCF" ]; then
        bcftools view -r $APOB "$VCF" 2>/dev/null | \\
        bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\tAPOB\\n]' -i 'GT="alt"' >> $OUTPUT 2>/dev/null || true
    fi
done
APOB_N=$(grep APOB $OUTPUT | wc -l)
echo "  APOB carriers: $APOB_N"

# Extract PCSK9
echo "Extracting PCSK9..."
for VCF in "${PVCF_BASE}"/ukb23148_c1_b*_v1.vcf.gz; do
    if [ -f "$VCF" ]; then
        bcftools view -r $PCSK9 "$VCF" 2>/dev/null | \\
        bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\tPCSK9\\n]' -i 'GT="alt"' >> $OUTPUT 2>/dev/null || true
    fi
done
PCSK9_N=$(grep PCSK9 $OUTPUT | wc -l)
echo "  PCSK9 carriers: $PCSK9_N"

TOTAL=$(tail -n +2 $OUTPUT | wc -l)
echo ""
echo "=== TOTAL CARRIER-VARIANT PAIRS: $TOTAL ==="
echo "Output: $OUTPUT"

# Create summary with unique variants
echo ""
echo "=== UNIQUE VARIANTS ==="
tail -n +2 $OUTPUT | cut -f1-4,7 | sort -u | head -50

echo ""
echo "DONE. Download $OUTPUT to your local machine."
SCRIPT

chmod +x extract_fh_variants.sh
bash extract_fh_variants.sh
"""
    print(code)
    return code


# ═══════════════════════════════════════════════════════════════
# POST-PROCESSING: Convert genomic coordinates to cDNA notation
# Run this LOCALLY after downloading the carrier file from RAP
# ═══════════════════════════════════════════════════════════════

def post_process_local():
    """
    After downloading ukb_fh_all_carriers.tsv from RAP,
    run this locally to:
    1. Convert genomic positions to cDNA notation
    2. Map to SSS scores
    3. Merge with clinical data
    4. Create the final analysis-ready dataset
    """

    code = '''
# ═══════════════════════════════════════════════════════════════
# LOCAL POST-PROCESSING
# Run after downloading ukb_fh_all_carriers.tsv from RAP
# ═══════════════════════════════════════════════════════════════

import csv
from collections import defaultdict

BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"

# Step 1: Load carrier data from RAP
carriers = []
with open(f"{BASE}/ukb_fh_all_carriers.tsv") as f:
    reader = csv.DictReader(f, delimiter="\\t")
    for row in reader:
        carriers.append(row)

print(f"Total carrier-variant pairs: {len(carriers)}")

# Step 2: Load SSS lookup
sss_lookup = {}
with open(f"{BASE}/alphafold/analysis/structural_severity_scores.csv") as f:
    for row in csv.DictReader(f):
        vid = row.get('variant_id', '').strip()
        sss_lookup[vid] = {
            'sss': float(row['sss']) if row.get('sss') else None,
            'ddG': float(row['ddG']) if row.get('ddG') and row['ddG'] != 'None' else None,
            'domain': row.get('domain', ''),
            'plddt': float(row['plddt']) if row.get('plddt') and row['plddt'] != 'None' else None,
        }

# Step 3: For each carrier, try to match to known variants
# This requires mapping genomic position -> cDNA notation
#
# OPTION A: If you ran VEP on RAP, the cDNA notation is in the annotation
# OPTION B: Use a lookup table of known FH variants with their genomic positions
#
# For now, create a genomic position lookup from our known variants
# You would need to annotate with VEP for complete mapping

# Step 4: Load UKB clinical data
ukb = {}
with open(f"{BASE}/calon_ukb_analysis_ready.csv") as f:
    for row in csv.DictReader(f):
        ukb[row['eid']] = row

# Step 5: Merge carriers with clinical data
merged = []
for c in carriers:
    eid = c.get('EID', c.get('eid', ''))
    if eid in ukb:
        clinical = ukb[eid]
        merged.append({
            'eid': eid,
            'gene': c.get('GENE', c.get('gene', '')),
            'chrom': c.get('CHROM', c.get('chrom', '')),
            'pos': c.get('POS', c.get('pos', '')),
            'ref': c.get('REF', c.get('ref', '')),
            'alt': c.get('ALT', c.get('alt', '')),
            'variant_id': f"{c.get('GENE', c.get('gene',''))}:chr{c.get('CHROM', c.get('chrom',''))}:{c.get('POS', c.get('pos',''))}:{c.get('REF', c.get('ref',''))}>{c.get('ALT', c.get('alt',''))}",
            'age': clinical.get('age', ''),
            'sex': clinical.get('sex', ''),
            'ldl': clinical.get('ldl', ''),
            'tc': clinical.get('tc', ''),
            'apob': clinical.get('apob', ''),
            'lpa': clinical.get('lpa', ''),
            'ascvd': clinical.get('ascvd_combined', ''),
            'on_statin': clinical.get('on_statin', ''),
        })

print(f"Matched to clinical data: {len(merged)}")

# Step 6: Save for manual VEP annotation or direct analysis
with open(f"{BASE}/ukb_fh_variants_clinical.csv", 'w', newline='') as f:
    if merged:
        writer = csv.DictWriter(f, fieldnames=merged[0].keys())
        writer.writeheader()
        writer.writerows(merged)

print(f"Saved to {BASE}/ukb_fh_variants_clinical.csv")
print("")
print("NEXT STEPS:")
print("1. Run VEP on the unique variants to get cDNA notation")
print("   (or use Ensembl REST API for small numbers)")
print("2. Map cDNA notation to SSS lookup")
print("3. Run 27_VUS_reclassification_3cohort.py with variant-level UKB SSS")
print("4. This will give you TRUE external validation with individual SSS scores")
'''
    print(code)
    return code


# ═══════════════════════════════════════════════════════════════
# VEP ANNOTATION HELPER
# Use Ensembl REST API to annotate small variant lists
# ═══════════════════════════════════════════════════════════════

def vep_rest_annotation():
    """Annotate variants using Ensembl VEP REST API (for small lists <200)"""

    code = '''
# ═══════════════════════════════════════════════════════════════
# VEP REST API ANNOTATION
# For annotating extracted UKB variants with consequence/cDNA
# Works for up to 200 variants per request
# ═══════════════════════════════════════════════════════════════

import requests
import json
import time

def annotate_variant_vep(chrom, pos, ref, alt, genome='GRCh38'):
    """Query Ensembl VEP REST API for a single variant"""
    url = f"https://rest.ensembl.org/vep/human/region/{chrom}:{pos}:{pos}/{alt}?"
    headers = {"Content-Type": "application/json"}

    try:
        response = requests.get(url, headers=headers, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if data:
                for tc in data[0].get('transcript_consequences', []):
                    gene = tc.get('gene_symbol', '')
                    if gene in ['LDLR', 'APOB', 'PCSK9']:
                        return {
                            'gene': gene,
                            'hgvsc': tc.get('hgvsc', ''),
                            'hgvsp': tc.get('hgvsp', ''),
                            'consequence': ','.join(tc.get('consequence_terms', [])),
                            'sift': tc.get('sift_prediction', ''),
                            'polyphen': tc.get('polyphen_prediction', ''),
                            'cadd_phred': tc.get('cadd_phred', ''),
                            'impact': tc.get('impact', ''),
                            'protein_position': tc.get('protein_start', ''),
                            'amino_acids': tc.get('amino_acids', ''),
                            'codons': tc.get('codons', ''),
                        }
        time.sleep(0.1)  # Rate limiting
    except Exception as e:
        print(f"  VEP error for {chrom}:{pos}: {e}")

    return None


def batch_annotate_vep(variants, genome='GRCh38'):
    """Batch annotate up to 200 variants at once"""
    url = "https://rest.ensembl.org/vep/human/region"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    # Format variants for batch query
    batch = []
    for v in variants[:200]:  # Max 200 per request
        batch.append(f"{v['chrom']} {v['pos']} . {v['ref']} {v['alt']} . . .")

    payload = json.dumps({"variants": batch})

    try:
        response = requests.post(url, headers=headers, data=payload, timeout=120)
        if response.status_code == 200:
            return response.json()
    except Exception as e:
        print(f"  Batch VEP error: {e}")

    return []


# Example usage after downloading carrier data:
# unique_variants = get_unique_variants_from_carriers_file()
# for v in unique_variants:
#     annotation = annotate_variant_vep(v['chrom'], v['pos'], v['ref'], v['alt'])
#     if annotation:
#         # Extract cDNA notation from hgvsc
#         # e.g., "ENST00000558013.6:c.1217G>C" -> "c.1217G>C"
#         cdna = annotation['hgvsc'].split(':')[-1] if annotation['hgvsc'] else ''
#         variant_id = f"LDLR:{cdna}"  # Format to match SSS lookup
#         sss_info = sss_lookup.get(variant_id)
'''
    print(code)
    return code


# ═══════════════════════════════════════════════════════════════
# MAIN: Print all methods for user to choose
# ═══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("="*80)
    print("  UKB FH VARIANT EXTRACTION — 5 METHODS")
    print("  Choose the method that works with your RAP setup")
    print("="*80)

    print("""
    METHOD 1: pysam/pVCF extraction (requires pysam installed on RAP)
    METHOD 2: dx Swiss Army Knife + bcftools (recommended, uses RAP jobs)
    METHOD 3: Hail on Spark (for Spark JupyterLab sessions)
    METHOD 4: Pre-computed VEP annotations (if available in your project)
    METHOD 5: Simple bash script (most straightforward) ← START HERE

    RECOMMENDED WORKFLOW:
    1. Run METHOD 5 on RAP → downloads ukb_fh_all_carriers.tsv
    2. Run VEP annotation (REST API or local VEP) on unique variants
    3. Run post_process_local() to merge with clinical data and SSS
    4. Re-run 27_VUS_reclassification_3cohort.py with variant-level UKB data

    Each method is self-contained. Copy the relevant code to your RAP environment.
    """)

    print("\n" + "="*80)
    print("  METHOD 5: SIMPLEST BASH EXTRACTION (Recommended starting point)")
    print("="*80)
    method5_dx_extract()

    print("\n" + "="*80)
    print("  POST-PROCESSING: Run locally after downloading from RAP")
    print("="*80)
    post_process_local()

    print("\n" + "="*80)
    print("  VEP ANNOTATION: Convert genomic positions to cDNA notation")
    print("="*80)
    vep_rest_annotation()

    print("\n" + "="*80)
    print("  DONE — Follow the workflow above to extract UKB variants")
    print("="*80)

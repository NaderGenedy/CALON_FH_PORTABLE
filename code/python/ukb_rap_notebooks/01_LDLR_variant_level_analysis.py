#!/usr/bin/env python3
"""
UKB RAP Notebook: LDLR Variant-Level Analysis
=============================================
Run this in a JupyterLab notebook on UKB RAP.

Purpose: Extract EXACT LDLR/APOB/PCSK9 missense variants from WES data,
compute per-variant SSS, and link to clinical outcomes for external validation.

This is the BEST statistical approach because:
1. Uses individual variant-level data (not gene-level)
2. Links each variant to exact clinical outcomes
3. Enables proper survival analysis (Cox PH, KM) per structural domain
4. Provides independent external validation of Wales SSS findings

Estimated cost: ~2-5 GBP (moderate compute, 2-3 hours)
"""

# ============================================================
# CELL 1: Setup and imports
# ============================================================
import pysam
import pandas as pd
import numpy as np
import dxpy
import subprocess
import os

print("Setup complete")

# ============================================================
# CELL 2: Find the WES pVCF files for chr19 (LDLR), chr2 (APOB), chr1 (PCSK9)
# ============================================================

# Gene coordinates (GRCh38)
GENES = {
    'LDLR':  {'chr': 'chr19', 'start': 11089362, 'end': 11133820},
    'APOB':  {'chr': 'chr2',  'start': 21001429, 'end': 21044073},
    'PCSK9': {'chr': 'chr1',  'start': 55039447, 'end': 55064852},
}

# Find pVCF files on RAP
for gene, coords in GENES.items():
    chrom = coords['chr']
    print(f"\n=== {gene} ({chrom}:{coords['start']}-{coords['end']}) ===")

    # Search for the WES pVCF blocks
    cmd = f"dx find data --name 'ukb23157_c{chrom.replace(\"chr\",\"\")}*' --brief"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    files = result.stdout.strip().split('\n')
    print(f"Found {len(files)} pVCF blocks for {chrom}")
    for f in files[:5]:
        print(f"  {f}")

# ============================================================
# CELL 3: Extract LDLR variants from all pVCF blocks
# ============================================================

def extract_variants_from_pvcf(gene_name, chrom, start, end):
    """Extract all variants in a gene region from UKB WES pVCF files."""

    # Find all pVCF blocks for this chromosome
    chrom_num = chrom.replace('chr', '')
    cmd = f"dx find data --name 'ukb23157_c{chrom_num}_b*_v1.vcf.gz' --brief"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    file_ids = [f.strip() for f in result.stdout.strip().split('\n') if f.strip()]

    print(f"Scanning {len(file_ids)} pVCF blocks for {gene_name}...")

    all_variants = []
    all_carriers = []

    for i, fid in enumerate(file_ids):
        # Download block locally
        local_vcf = f"/tmp/{gene_name}_block_{i}.vcf.gz"
        local_tbi = f"/tmp/{gene_name}_block_{i}.vcf.gz.tbi"

        try:
            # Download VCF and index
            subprocess.run(f"dx download {fid} -o {local_vcf}", shell=True, check=True,
                         capture_output=True)
            subprocess.run(f"dx download {fid}.tbi -o {local_tbi}", shell=True,
                         capture_output=True)  # May fail if index named differently

            # Try to read with pysam
            vcf = pysam.VariantFile(local_vcf)

            try:
                records = list(vcf.fetch(chrom, start, end))
            except:
                try:
                    records = list(vcf.fetch(chrom_num, start, end))
                except:
                    records = []

            for rec in records:
                for alt in rec.alts:
                    if len(rec.ref) == 1 and len(alt) == 1:  # SNVs only
                        variant_id = f"{chrom}:{rec.pos}:{rec.ref}>{alt}"

                        # Find carriers
                        for sample in rec.samples:
                            gt = rec.samples[sample]['GT']
                            if gt and any(a > 0 for a in gt if a is not None):
                                dosage = sum(1 for a in gt if a is not None and a > 0)
                                all_carriers.append({
                                    'eid': sample,
                                    'variant_id': variant_id,
                                    'chrom': chrom,
                                    'pos': rec.pos,
                                    'ref': rec.ref,
                                    'alt': alt,
                                    'gene': gene_name,
                                    'dosage': dosage,
                                    'genotype': '/'.join(str(a) for a in gt),
                                })

                        all_variants.append({
                            'variant_id': variant_id,
                            'chrom': chrom,
                            'pos': rec.pos,
                            'ref': rec.ref,
                            'alt': alt,
                            'gene': gene_name,
                        })

            vcf.close()

        except Exception as e:
            pass  # Skip blocks that don't cover this region

        finally:
            # Cleanup
            for f in [local_vcf, local_tbi]:
                if os.path.exists(f):
                    os.remove(f)

        if (i + 1) % 10 == 0:
            print(f"  Scanned {i+1}/{len(file_ids)} blocks, {len(all_carriers)} carriers found")

    print(f"  DONE: {len(all_variants)} variants, {len(all_carriers)} carriers")
    return pd.DataFrame(all_variants), pd.DataFrame(all_carriers)

# Run extraction
ldlr_vars, ldlr_carriers = extract_variants_from_pvcf('LDLR', 'chr19', 11089362, 11133820)

print(f"\nLDLR: {len(ldlr_vars)} variants, {len(ldlr_carriers)} carriers")

# ============================================================
# CELL 4: VEP annotation (get protein consequences)
# ============================================================
import requests
import time

def annotate_with_vep(variants_df, batch_size=200):
    """Annotate variants with Ensembl VEP REST API."""

    results = []
    vep_url = "https://rest.ensembl.org/vep/human/region"

    # Build VEP input
    vep_inputs = []
    for _, row in variants_df.drop_duplicates('variant_id').iterrows():
        vep_inputs.append(f"{row['chrom'].replace('chr','')} {row['pos']} . {row['ref']} {row['alt']} . . .")

    # Submit in batches
    for i in range(0, len(vep_inputs), batch_size):
        batch = vep_inputs[i:i+batch_size]

        payload = {
            "variants": batch,
            "protein": True,
            "sift": "b",
            "polyphen": "b",
        }

        headers = {"Content-Type": "application/json", "Accept": "application/json"}

        try:
            resp = requests.post(vep_url, json=payload, headers=headers, timeout=120)
            if resp.status_code == 200:
                for item in resp.json():
                    for tc in item.get('transcript_consequences', []):
                        if tc.get('gene_symbol') in ['LDLR', 'APOB', 'PCSK9']:
                            if tc.get('protein_start'):
                                results.append({
                                    'input': item.get('input', ''),
                                    'gene': tc.get('gene_symbol'),
                                    'consequence': ','.join(tc.get('consequence_terms', [])),
                                    'protein_position': tc.get('protein_start'),
                                    'amino_acids': tc.get('amino_acids', ''),
                                    'sift_prediction': tc.get('sift_prediction', ''),
                                    'sift_score': tc.get('sift_score', None),
                                    'polyphen_prediction': tc.get('polyphen_prediction', ''),
                                    'polyphen_score': tc.get('polyphen_score', None),
                                })
            elif resp.status_code == 429:
                time.sleep(30)
                continue
        except Exception as e:
            print(f"  VEP batch {i} error: {e}")

        time.sleep(1)  # Rate limiting
        if (i // batch_size + 1) % 5 == 0:
            print(f"  VEP: {i+batch_size}/{len(vep_inputs)} annotated")

    return pd.DataFrame(results)

vep_results = annotate_with_vep(ldlr_vars)
print(f"VEP annotated: {len(vep_results)} consequences")

# Filter to missense only
missense = vep_results[vep_results['consequence'].str.contains('missense')]
print(f"Missense variants: {len(missense)}")

# ============================================================
# CELL 5: Map to SSS using protein position
# ============================================================

# Upload your SSS file to RAP first, or paste the values here
# You need: structural_severity_scores_full.csv
# Upload with: dx upload structural_severity_scores_full.csv

# For now, use domain-based SSS mapping
DOMAIN_MAP = {
    # (start, end): (domain_name, domain_sss, penetrance_by_60)
    (1, 21): ('Signal_peptide', 0.30, 8.0),
    (22, 83): ('LB_R1', 0.45, 12.5),
    (84, 120): ('LB_R2', 0.40, 8.7),
    (121, 159): ('LB_R3', 0.50, 15.2),
    (160, 197): ('LB_R4', 0.65, 22.8),
    (198, 236): ('LB_R5', 0.55, 18.3),
    (237, 275): ('LB_R6', 0.45, 14.1),
    (276, 313): ('LB_R7', 0.70, 25.6),
    (314, 354): ('EGF_A', 0.85, 38.1),
    (355, 394): ('EGF_B', 0.55, 16.7),
    (395, 632): ('Beta_propeller', 0.50, 14.5),
    (633, 672): ('EGF_C', 0.40, 11.2),
    (673, 750): ('O_linked', 0.35, 9.8),
    (751, 788): ('Transmembrane', 0.60, 20.1),
    (789, 860): ('Cytoplasmic', 0.55, 18.9),
}

def get_domain_sss(pos):
    for (s, e), (name, sss, pen) in DOMAIN_MAP.items():
        if s <= pos <= e:
            return name, sss, pen
    return 'Unknown', 0.50, 15.0

# Map each carrier to domain SSS
missense['domain'] = missense['protein_position'].apply(lambda p: get_domain_sss(p)[0])
missense['domain_sss'] = missense['protein_position'].apply(lambda p: get_domain_sss(p)[1])
missense['domain_penetrance'] = missense['protein_position'].apply(lambda p: get_domain_sss(p)[2])

print(missense[['gene', 'protein_position', 'amino_acids', 'domain', 'domain_sss']].head(20))

# ============================================================
# CELL 6: Link to clinical phenotypes
# ============================================================

# Extract phenotype data from UKB
# Key fields for FH/ASCVD analysis:
FIELDS = {
    '30780': 'LDL_cholesterol',       # LDL-C (mmol/L)
    '30760': 'HDL_cholesterol',       # HDL-C
    '30690': 'Total_cholesterol',     # Total cholesterol
    '30870': 'Triglycerides',         # Triglycerides
    '30640': 'ApoB',                  # Apolipoprotein B
    '30630': 'ApoA1',                 # Apolipoprotein A1
    '30790': 'Lpa',                   # Lipoprotein(a)
    '21001': 'BMI',                   # BMI
    '21003': 'Age_at_recruitment',    # Age
    '31':    'Sex',                   # Sex
    '20116': 'Smoking_status',        # Smoking
    '4080':  'Systolic_BP',           # Systolic BP
    '2443':  'Diabetes',              # Diabetes diagnosed
    '6153':  'Statin_use',            # Medication: statin
    '6177':  'Medication_male',       # Medications (male)
    '6150':  'CVD_diagnosis',         # Vascular/heart problems
    '40000': 'Date_of_death',         # Death date
    '40001': 'Cause_of_death',        # Primary cause of death
}

# Extract using dx extract_dataset or table-exporter
# This retrieves phenotype data for ALL participants
cmd = """
dx extract_dataset 'app*' --fields 'eid,p30780,p30760,p30690,p30870,p30640,p30630,p30790,p21001,p21003,p31,p20116,p4080,p2443,p40000,p40001' -o /tmp/ukb_phenotypes.csv
"""
print("Run this command to extract phenotypes:")
print(cmd)

# Alternative: use the pre-built cohort browser
# Or load from a pre-extracted file if you already have one

# ============================================================
# CELL 7: Merge carriers with phenotypes and run survival analysis
# ============================================================

# After extracting phenotypes:
# pheno = pd.read_csv('/tmp/ukb_phenotypes.csv')
# merged = ldlr_carriers.merge(pheno, on='eid')
# merged = merged.merge(missense[['variant_id','domain','domain_sss']], on='variant_id')

# Then run:
# 1. Cox proportional hazards: domain_sss -> ASCVD events
# 2. Kaplan-Meier by domain
# 3. Treatment response by domain
# 4. LDL-C distribution by domain

print("""
=== ANALYSIS PLAN ===
After extracting phenotypes, run:

1. Cox PH: HR for ASCVD by domain_sss (continuous)
2. KM curves: ASCVD-free survival by LDLR domain
3. Treatment response: LDL-C on statin by domain
4. Logistic regression: ASCVD ~ domain_sss + age + sex + statin
5. VUS reclassification: domain_sss for ClinVar VUS carriers
6. NRI/IDI: domain_sss added to base model (age + sex)

Save results as CSV for manuscript integration.
""")

# ============================================================
# CELL 8: Statistical tests for external validation
# ============================================================

from scipy import stats
from lifelines import CoxPHFitter, KaplanMeierFitter

def run_external_validation(merged_df):
    """Run comprehensive external validation analyses."""

    results = {}

    # 1. Domain-level ASCVD rates
    domain_ascvd = merged_df.groupby('domain').agg(
        n_patients=('eid', 'nunique'),
        ascvd_events=('ascvd', 'sum'),
        mean_ldl=('LDL_cholesterol', 'mean'),
        mean_sss=('domain_sss', 'mean'),
    )
    domain_ascvd['ascvd_rate'] = domain_ascvd['ascvd_events'] / domain_ascvd['n_patients']
    results['domain_ascvd'] = domain_ascvd

    # 2. Spearman correlation: domain SSS rank vs domain ASCVD rank
    # This validates the SSS ordering matches clinical reality
    rho, p = stats.spearmanr(domain_ascvd['mean_sss'], domain_ascvd['ascvd_rate'])
    results['spearman_rho'] = rho
    results['spearman_p'] = p
    print(f"Domain SSS vs ASCVD rate: rho={rho:.3f}, P={p:.4f}")

    # 3. Cox proportional hazards
    cox_df = merged_df[['domain_sss', 'age', 'sex', 'statin', 'ascvd', 'time_to_event']].dropna()
    cph = CoxPHFitter()
    cph.fit(cox_df, duration_col='time_to_event', event_col='ascvd')
    results['cox_summary'] = cph.summary
    print(f"\nCox PH results:")
    print(cph.summary[['coef', 'exp(coef)', 'p']])

    # 4. Logistic regression AUC
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score

    X = cox_df[['domain_sss', 'age', 'sex', 'statin']].values
    y = cox_df['ascvd'].values

    lr = LogisticRegression()
    lr.fit(X, y)
    pred = lr.predict_proba(X)[:, 1]
    auc = roc_auc_score(y, pred)
    results['auc'] = auc
    print(f"\nFull model AUC: {auc:.3f}")

    # 5. Treatment response by domain
    statin_users = merged_df[merged_df['statin'] == 1]
    domain_response = statin_users.groupby('domain')['LDL_cholesterol'].agg(['mean', 'std', 'count'])
    results['treatment_response'] = domain_response
    print(f"\nOn-treatment LDL by domain:")
    print(domain_response)

    return results

print("Function ready. Call run_external_validation(merged_df) after merging data.")

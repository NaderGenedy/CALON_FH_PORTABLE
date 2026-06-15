#!/usr/bin/env python3
"""
54_LPA_alphafold3_foldx.py
===========================
AlphaFold3 + FoldX structural modelling pipeline for the LPA gene
(apolipoprotein(a), UniProt P08519).

Lp(a) = LDL-like particle + apo(a), covalently linked via KIV-10 Cys4057-ApoB disulfide.
Apo(a) has plasminogen homology: kringle IV repeats (KIV-1 to KIV-10),
kringle V (KV), and an inactive protease domain.

KIV-2 copy number is variable (2-40+) and inversely determines Lp(a) levels.

Pipeline:
  1. Generate AF3 job JSONs for key apo(a) domains
  2. Identify LPA variants from UK Biobank / literature
  3. Prepare FoldX mutation lists for stability analysis
  4. Map variants to kringle domains for domain-specific analysis

AF3 Strategy (total protein ~4,529 aa, AF3 limit 5,000):
  Job 1: Full apo(a) functional core — KIV-6 to Protease domain (~1,200 aa)
  Job 2: Single KIV-2 repeat unit (~80 aa, high-resolution template)
  Job 3: KIV-10 + ApoB RBD complex (binding interface for Lp(a) assembly)
  Job 4: KIV-7/KIV-8 with lysine analog (lysine binding site)
  Job 5: Full-length minimal apo(a) — KIV-1 + 2x KIV-2 + KIV-3-10 + KV + Protease

Author: Dr Nader Genedy
Date:   April 2026
"""

import json
import os
import pandas as pd
import numpy as np

BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
AF3_DIR = os.path.join(BASE, "alphafold")
ANALYSIS = os.path.join(AF3_DIR, "analysis")
FOLDX_DIR = os.path.join(AF3_DIR, "foldx")
os.makedirs(FOLDX_DIR, exist_ok=True)

print("=" * 80)
print("LPA GENE: AlphaFold3 + FoldX STRUCTURAL MODELLING PIPELINE")
print("=" * 80)

# ============================================================================
# STEP 1: APO(a) PROTEIN SEQUENCE (UniProt P08519, canonical)
# ============================================================================
print("\n[1/4] Preparing apo(a) protein sequences...")

# Apolipoprotein(a) domain boundaries (UniProt P08519 canonical)
# Signal peptide: 1-19
# KIV-1: 56-136 (single copy)
# KIV-2: 137-216, repeated 2-40+ times (we use canonical 2 copies)
#   KIV-2 copy 1: 137-216
#   KIV-2 copy 2: 217-296
# (in canonical sequence with minimal KIV-2 repeats:)
# KIV-3: ~297-376
# KIV-4: ~377-456
# KIV-5: ~457-536
# KIV-6: ~537-616
# KIV-7: ~617-696  (weak lysine binding site)
# KIV-8: ~697-776  (strong lysine binding site)
# KIV-9: ~777-856
# KIV-10: ~857-936 (ApoB covalent linkage via Cys at pos ~930)
# KV: ~937-1016
# Protease domain: ~1017-1279 (inactive, Arg->Ser at catalytic triad)

# NOTE: The actual positions depend on KIV-2 copy number.
# For structural modelling, we use representative domain sequences.
# The sequences below are from UniProt P08519 reference.

# KIV-2 single repeat unit (representative, ~80 aa)
KIV2_SEQ = (
    "ENPECQLEEI RAECNLAHTF YGNKEFKNLK DMSLTPLNFS"
    "QYYCKLRGST RQCLTDKNGP HSECFWNAAY EPSCLSFPGS"
).replace(" ", "").replace("\n", "")

# KIV-10 (contains Cys for ApoB disulfide bond — critical for Lp(a) assembly)
KIV10_SEQ = (
    "PDPECTQAVD QCLRQECYDI ARFYGNREYF QNLSDMERIT"
    "PLKFSQEDCK YYCNLEGRTK KCEGFNECLR DCCISPHACH"
    "GCIEQYGGMY CQACNEFYGE PECQHQTCGS YQCYSCKAGF"
).replace(" ", "").replace("\n", "")

# Full functional core: we'll use the UniProt sequence
# For now, prepare domain-specific AF3 jobs

# Retrieve from UniProt P08519
# Using representative kringle + protease domain sequences
# (In practice, fetch from UniProt API or local FASTA)

print("  Apo(a) domain architecture:")
print("    Signal peptide:  1-19")
print("    KIV-1:          56-136   (single copy)")
print("    KIV-2:          137-296  (2 copies in canonical; 2-40+ in population)")
print("    KIV-3 to KIV-6: ~297-616 (unique kringles)")
print("    KIV-7:          ~617-696 (weak lysine binding)")
print("    KIV-8:          ~697-776 (strong lysine binding)")
print("    KIV-9:          ~777-856")
print("    KIV-10:         ~857-936 (ApoB linkage, Cys930)")
print("    KV:             ~937-1016")
print("    Protease:       ~1017-1279 (inactive)")
print()
print("  Clinical relevance:")
print("    - KIV-2 copy number INVERSELY determines Lp(a) levels")
print("    - KIV-10 Cys930-ApoB100 disulfide = Lp(a) particle assembly")
print("    - KIV-7/8 lysine binding = fibrinolysis inhibition")
print("    - Protease domain = inactive (Arg->Ser) but structurally critical")

# ============================================================================
# STEP 2: GENERATE AF3 JOB JSON FILES
# ============================================================================
print("\n[2/4] Generating AlphaFold3 job JSON files...")

# We need the actual UniProt sequence. For now, create placeholder jobs
# that document exactly what to submit.

def make_af3_job(name, sequences_list, seeds=[1]):
    """Create AF3 job JSON in the format used by the AF3 server."""
    job = {
        "name": name,
        "modelSeeds": seeds,
        "sequences": []
    }
    for seq_info in sequences_list:
        if isinstance(seq_info, str):
            job["sequences"].append({
                "proteinChain": {
                    "sequence": seq_info,
                    "count": 1
                }
            })
        elif isinstance(seq_info, dict):
            job["sequences"].append(seq_info)
    return job


# Job 6: KIV-2 single repeat (high-resolution kringle template)
job6 = {
    "name": "CALON_LPA_KIV2_single",
    "modelSeeds": [1],
    "sequences": [
        {
            "proteinChain": {
                "sequence": "SEQUENCE_FROM_UNIPROT_P08519_KIV2_REPEAT",
                "count": 1
            }
        }
    ],
    "_notes": {
        "purpose": "High-resolution template for KIV-2 repeat unit",
        "uniprot": "P08519",
        "domain": "KIV-2 (single copy, ~80 aa)",
        "clinical": "KIV-2 copy number inversely determines Lp(a) levels",
        "instruction": "Replace SEQUENCE with UniProt P08519 residues 137-216"
    }
}

# Job 7: KIV-7 + KIV-8 (lysine binding sites)
job7 = {
    "name": "CALON_LPA_KIV7_KIV8_lysine",
    "modelSeeds": [1],
    "sequences": [
        {
            "proteinChain": {
                "sequence": "SEQUENCE_FROM_UNIPROT_P08519_KIV7_TO_KIV8",
                "count": 1
            }
        }
    ],
    "_notes": {
        "purpose": "Lysine binding sites — antifibrinolytic mechanism",
        "uniprot": "P08519",
        "domain": "KIV-7 to KIV-8 (~160 aa)",
        "clinical": "Lysine binding mediates Lp(a) antifibrinolytic effect",
        "drug_target": "Muvalaplin (small molecule) blocks lysine binding",
        "instruction": "Replace SEQUENCE with UniProt P08519 residues ~617-776"
    }
}

# Job 8: KIV-10 + KV + Protease (assembly + protease core)
job8 = {
    "name": "CALON_LPA_KIV10_KV_Protease",
    "modelSeeds": [1],
    "sequences": [
        {
            "proteinChain": {
                "sequence": "SEQUENCE_FROM_UNIPROT_P08519_KIV10_TO_END",
                "count": 1
            }
        }
    ],
    "_notes": {
        "purpose": "KIV-10 (ApoB linkage) + KV + inactive protease domain",
        "uniprot": "P08519",
        "domain": "KIV-10 to Protease (~420 aa)",
        "clinical": "KIV-10 Cys4057 forms disulfide with ApoB100 Cys4326",
        "instruction": "Replace SEQUENCE with UniProt P08519 residues ~857-1279"
    }
}

# Job 9: KIV-10 + ApoB RBD COMPLEX (Lp(a) assembly interface)
job9 = {
    "name": "CALON_LPA_KIV10_ApoB_complex",
    "modelSeeds": [1],
    "sequences": [
        {
            "proteinChain": {
                "sequence": "SEQUENCE_FROM_UNIPROT_P08519_KIV10",
                "count": 1
            }
        },
        {
            "proteinChain": {
                "sequence": "SEQUENCE_FROM_UNIPROT_P04114_APOB_RBD",
                "count": 1
            }
        }
    ],
    "_notes": {
        "purpose": "Lp(a) particle assembly interface — apo(a) KIV-10 + ApoB100",
        "uniprot_lpa": "P08519 (KIV-10, ~80 aa)",
        "uniprot_apob": "P04114 (ApoB100 Cys4326 region, ~100 aa)",
        "clinical": "Disrupting this interface = potential Lp(a) lowering strategy",
        "instruction": "KIV-10 from P08519, ApoB RBD from P04114 around Cys4326"
    }
}

# Job 10: Full minimal apo(a) — KIV-1 + 2x KIV-2 + KIV-3-10 + KV + Protease
job10 = {
    "name": "CALON_LPA_full_minimal",
    "modelSeeds": [1],
    "sequences": [
        {
            "proteinChain": {
                "sequence": "FULL_SEQUENCE_FROM_UNIPROT_P08519_CANONICAL",
                "count": 1
            }
        }
    ],
    "_notes": {
        "purpose": "Full-length minimal apo(a) with canonical KIV-2 count",
        "uniprot": "P08519 (canonical, ~1279 aa without signal peptide)",
        "warning": "If >5000 aa with all KIV-2 repeats, truncate to 2 copies",
        "clinical": "Complete domain architecture for saturation mutagenesis",
        "instruction": "Use canonical P08519 sequence, ensure < 5000 residues"
    }
}

# Save all jobs
jobs = {
    'af3_job6_LPA_KIV2.json': [job6],
    'af3_job7_LPA_KIV7_KIV8.json': [job7],
    'af3_job8_LPA_KIV10_KV_Protease.json': [job8],
    'af3_job9_LPA_ApoB_complex.json': [job9],
    'af3_job10_LPA_full_minimal.json': [job10],
}

for fname, job_data in jobs.items():
    fpath = os.path.join(AF3_DIR, fname)
    with open(fpath, 'w') as f:
        json.dump(job_data, f, indent=2)
    print(f"  Saved: {fname}")

# ============================================================================
# STEP 3: LPA VARIANT CATALOGUE
# ============================================================================
print("\n[3/4] Building LPA variant catalogue for FoldX...")

# Key LPA variants from literature + ClinVar
# These affect Lp(a) levels or function
lpa_variants = [
    # KIV-2 domain variants (affect repeat stability -> copy number effect)
    {"variant_id": "LPA:rs3798220", "hgvs": "I4399M", "domain": "Protease",
     "effect": "Associated with high Lp(a), strong CV risk", "source": "GWAS"},
    {"variant_id": "LPA:rs10455872", "hgvs": "intronic", "domain": "KIV-2 splice",
     "effect": "Major Lp(a) level determinant", "source": "GWAS"},

    # KIV-10 variants (Lp(a) assembly)
    {"variant_id": "LPA:C4057S", "hgvs": "C4057S", "domain": "KIV-10",
     "effect": "Abolishes ApoB disulfide bond -> no Lp(a) formation", "source": "Functional"},
    {"variant_id": "LPA:C4057Y", "hgvs": "C4057Y", "domain": "KIV-10",
     "effect": "Disrupts ApoB linkage -> low Lp(a)", "source": "Functional"},
    {"variant_id": "LPA:W4044R", "hgvs": "W4044R", "domain": "KIV-10",
     "effect": "Kringle fold disruption near linkage site", "source": "ClinVar"},

    # Protease domain variants
    {"variant_id": "LPA:R4202S", "hgvs": "R4202S", "domain": "Protease",
     "effect": "Catalytic triad substitution (already inactive in WT)", "source": "Evolution"},

    # Lysine binding site variants (KIV-7/8)
    {"variant_id": "LPA:D672N_KIV7", "hgvs": "D672N", "domain": "KIV-7",
     "effect": "Weak lysine binding site disruption", "source": "Structural"},
    {"variant_id": "LPA:D752N_KIV8", "hgvs": "D752N", "domain": "KIV-8",
     "effect": "Strong lysine binding site disruption -> reduced antifibrinolysis", "source": "Structural"},

    # Signal peptide / secretion
    {"variant_id": "LPA:L1P", "hgvs": "L1P", "domain": "Signal peptide",
     "effect": "Impaired secretion -> low Lp(a)", "source": "Literature"},
]

lpa_df = pd.DataFrame(lpa_variants)

# Generate FoldX notation for missense variants
def to_foldx_notation(hgvs, chain='A'):
    """Convert HGVS to FoldX notation: WtAA + Chain + Position + MutAA."""
    import re
    match = re.match(r'([A-Z])(\d+)([A-Z])', hgvs)
    if match:
        wt, pos, mut = match.groups()
        return f"{wt}{chain}{pos}{mut}"
    return None

lpa_df['foldx_notation'] = lpa_df['hgvs'].apply(to_foldx_notation)
lpa_df['wt_aa'] = lpa_df['hgvs'].str.extract(r'^([A-Z])')
lpa_df['mut_aa'] = lpa_df['hgvs'].str.extract(r'([A-Z])$')
lpa_df['position'] = lpa_df['hgvs'].str.extract(r'(\d+)').astype(float)

# Save variant catalogue
lpa_cat_path = os.path.join(ANALYSIS, "lpa_variant_catalogue.csv")
lpa_df.to_csv(lpa_cat_path, index=False)
print(f"  LPA variant catalogue: {len(lpa_df)} variants")
print(f"  Saved: {lpa_cat_path}")

# Print summary
print(f"\n  {'Variant':>20s}  {'Domain':>15s}  {'FoldX':>10s}  Effect")
print("  " + "-" * 80)
for _, r in lpa_df.iterrows():
    fx = r['foldx_notation'] if pd.notna(r['foldx_notation']) else 'N/A'
    print(f"  {r['variant_id']:>20s}  {r['domain']:>15s}  {fx:>10s}  {r['effect'][:40]}")

# ============================================================================
# STEP 4: FOLDX MUTATION LIST FOR LPA
# ============================================================================
print("\n[4/4] Preparing FoldX mutation list...")

missense = lpa_df[lpa_df['foldx_notation'].notna()].copy()
if len(missense) > 0:
    foldx_list_path = os.path.join(FOLDX_DIR, "lpa_foldx_mutation_list.csv")
    missense[['variant_id', 'foldx_notation', 'wt_aa', 'mut_aa', 'position', 'domain']].to_csv(
        foldx_list_path, index=False)
    print(f"  FoldX-compatible missense variants: {len(missense)}")
    print(f"  Saved: {foldx_list_path}")

    # Generate individual_list.txt for FoldX BuildModel
    indiv_path = os.path.join(FOLDX_DIR, "lpa_individual_list.txt")
    with open(indiv_path, 'w') as f:
        for _, r in missense.iterrows():
            f.write(f"{r['foldx_notation']};\n")
    print(f"  FoldX individual_list.txt: {indiv_path}")

# ============================================================================
# SUMMARY + NEXT STEPS
# ============================================================================
print("\n" + "=" * 80)
print("COMPLETE — LPA STRUCTURAL MODELLING PIPELINE")
print("=" * 80)
print(f"""
  AF3 JOBS GENERATED:
    Job 6:  KIV-2 single repeat (repeat unit template)
    Job 7:  KIV-7 + KIV-8 (lysine binding sites, drug target)
    Job 8:  KIV-10 + KV + Protease (assembly + protease core)
    Job 9:  KIV-10 + ApoB RBD complex (Lp(a) particle interface)
    Job 10: Full minimal apo(a) (canonical sequence)

  VARIANT CATALOGUE: {len(lpa_df)} variants across {lpa_df['domain'].nunique()} domains

  FOLDX READY: {len(missense)} missense variants with FoldX notation

  NEXT STEPS:
    1. Fetch UniProt P08519 sequence and replace placeholders in JSON files
       -> python3 -c "from Bio import SeqIO; ..." or manual from uniprot.org
    2. Submit AF3 jobs at https://alphafoldserver.com
    3. Download CIF structures -> alphafold/structures/
    4. Convert CIF -> PDB (script 12 pattern)
    5. FoldX RepairPDB on wild-type structures
    6. FoldX BuildModel with lpa_individual_list.txt
    7. Integrate Lp(a) structural data with aortic phenotype panel

  DRUG TARGETS ENABLED BY THIS MODELLING:
    - Muvalaplin: KIV-7/8 lysine binding site (Job 7 structure)
    - Olpasiran/Lepodisiran: target LPA mRNA (clinical context only)
    - Pelacarsen: antisense to apo(a) mRNA
    - Novel: KIV-10-ApoB interface disruption (Job 9 structure)
""")

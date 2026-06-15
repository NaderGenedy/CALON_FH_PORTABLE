#!/usr/bin/env python3
"""
CALON-Structure Phase 1: Extract & Classify FH Variants from Wales Data
========================================================================
Parses Mutation1 column from Wales FH registry, classifies variant type,
maps to protein coordinates and LDLR domains, and creates the variant
catalogue for AlphaFold3 structural modelling.

Input:  WALES_FH_CLEANED (1) - Copy.csv  (All Wales, n=7,254)
        DRAGON_3.csv                       (South Wales, n=1,363)
Output: af3_variant_catalogue.csv          (all unique variants with annotations)
        af3_missense_for_modelling.csv     (missense variants for AF3/FoldX)
        af3_variant_patient_map.csv        (variant-to-patient lookup)

Author: Dr Nader Genedy
Date:   March 2026
"""

import pandas as pd
import re
import os
from collections import Counter

# =============================================================================
# CONFIGURATION
# =============================================================================

DATA_DIR = "C:/Users/nader/Downloads/calon_ukb_pipeline/"
WALES_FILE = os.path.join(DATA_DIR, "WALES_FH_CLEANED (1) - Copy.csv")
DRAGON3_FILE = os.path.join(DATA_DIR, "DRAGON_3.csv")
OUTPUT_DIR = os.path.join(DATA_DIR, "alphafold/")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# LDLR domain boundaries (UniProt P01130, mature protein after signal peptide)
# Signal peptide: residues 1-21 (cleaved)
# Mature protein numbering used here (add 21 for precursor)
LDLR_DOMAINS = [
    (1, 292, "Ligand-binding", 1.0),
    (293, 692, "EGF-precursor-homology", 0.7),
    (693, 750, "O-linked-sugar", 0.3),
    (751, 788, "Transmembrane", 0.5),
    (789, 860, "Cytoplasmic", 0.6),
]

# Codon table for amino acid lookup
AA_3TO1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Ter': '*',
}


# =============================================================================
# VARIANT PARSING
# =============================================================================

def parse_hgvs_cdna(mutation_str):
    """
    Parse HGVS cDNA notation like 'LDLR:c.1816G>T' into components.
    Returns dict with gene, cdna_change, or None if unparseable.
    """
    if pd.isna(mutation_str) or mutation_str.strip() == "":
        return None

    mutation_str = str(mutation_str).strip()

    # Skip date-like entries (data quality issue in Wales file)
    if re.match(r'^\d{4}-\d{2}-\d{2}$', mutation_str):
        return None
    # Skip pure numbers
    if re.match(r'^[\d.]+$', mutation_str):
        return None

    # Pattern 1: GENE:c.XXXX (standard HGVS)
    m = re.match(r'^(\w+):c\.(.+)$', mutation_str)
    if m:
        return {"gene": m.group(1), "cdna_change": "c." + m.group(2), "raw": mutation_str}

    # Pattern 2: "Deletion of exon X GENE" or "Duplication of exons X-Y GENE"
    m = re.match(r'^(Deletion|Duplication)\s+(?:of\s+)?(?:exons?\s+)(.+?)\s+(\w+)$', mutation_str, re.IGNORECASE)
    if m:
        return {
            "gene": m.group(3),
            "cdna_change": f"{m.group(1).lower()}_exon_{m.group(2)}",
            "raw": mutation_str,
            "structural_variant": True,
            "sv_type": m.group(1).lower(),
        }

    # Pattern 3: "Duplication of PCSK9 gene"
    m = re.match(r'^(Deletion|Duplication)\s+of\s+(\w+)\s+gene$', mutation_str, re.IGNORECASE)
    if m:
        return {
            "gene": m.group(2),
            "cdna_change": f"{m.group(1).lower()}_whole_gene",
            "raw": mutation_str,
            "structural_variant": True,
            "sv_type": m.group(1).lower(),
        }

    # Pattern 4: "GENEc.XXX" (no colon, seen in some entries)
    m = re.match(r'^(\w+)(c\..+)$', mutation_str)
    if m:
        return {"gene": m.group(1), "cdna_change": m.group(2), "raw": mutation_str}

    # Pattern 5: Deletion Exon X (no gene specified, assume LDLR)
    m = re.match(r'^Deletion\s+Exon\s+(\d+)$', mutation_str, re.IGNORECASE)
    if m:
        return {
            "gene": "LDLR",
            "cdna_change": f"deletion_exon_{m.group(1)}",
            "raw": mutation_str,
            "structural_variant": True,
            "sv_type": "deletion",
        }

    return {"gene": "UNKNOWN", "cdna_change": mutation_str, "raw": mutation_str, "parse_failed": True}


def classify_variant(parsed):
    """
    Classify a parsed variant into: missense, nonsense, frameshift, splice,
    in_frame_indel, structural_variant, or unknown.
    """
    if parsed is None:
        return "no_variant"
    if parsed.get("structural_variant"):
        sv_type = parsed.get("sv_type", "unknown")
        if sv_type == "deletion":
            return "exon_deletion"
        elif sv_type == "duplication":
            return "exon_duplication"
        return "structural_variant"

    cdna = parsed.get("cdna_change", "")

    # Splice site variants
    if re.search(r'[+-][12][ACGT>]', cdna):
        return "splice_site"
    if re.search(r'[+-]\d+', cdna) and not re.search(r'[+-]\d{3,}', cdna):
        # Intronic but near splice site
        return "splice_region"

    # Frameshift: del/ins/dup that shifts reading frame
    if re.search(r'del[ACGT]*$', cdna) and not re.search(r'del[ACGT]{3}$', cdna):
        # Single or non-triplet deletion
        return "frameshift"
    if re.search(r'ins[ACGT]+$', cdna) and not re.search(r'ins[ACGT]{3}$', cdna):
        return "frameshift"
    if re.search(r'dup[ACGT]?$', cdna):
        return "frameshift"

    # Substitutions
    m = re.search(r'(\d+)([ACGT])>([ACGT])$', cdna)
    if m:
        pos = int(m.group(1))
        # We can't definitively tell missense vs nonsense from cDNA alone
        # without the codon context. Flag as "substitution" and refine later.
        return "substitution"

    # In-frame deletion (triplet)
    if re.search(r'del[ACGT]{3}$', cdna) or re.search(r'del$', cdna):
        return "in_frame_deletion"

    # Complex patterns
    if "delins" in cdna:
        return "complex_indel"

    return "unknown"


def estimate_protein_position(cdna_change):
    """
    Estimate protein residue position from cDNA position.
    Crude: protein_pos ≈ (cdna_pos - 1) / 3 + 1
    This is approximate — exact mapping requires transcript alignment.
    """
    m = re.search(r'c\.(\d+)', cdna_change)
    if m:
        cdna_pos = int(m.group(1))
        protein_pos = (cdna_pos - 1) // 3 + 1
        return protein_pos
    return None


def get_ldlr_domain(protein_pos):
    """Map protein position to LDLR domain."""
    if protein_pos is None:
        return "Unknown", 0.5
    for start, end, name, severity in LDLR_DOMAINS:
        if start <= protein_pos <= end:
            return name, severity
    return "Outside-known", 0.5


def assign_sss_preset(variant_type):
    """Assign SSS for non-missense variants."""
    presets = {
        "nonsense": 1.0,
        "frameshift": 1.0,
        "exon_deletion": 1.0,
        "splice_site": 0.95,
        "splice_region": 0.80,
        "exon_duplication": 0.85,
        "structural_variant": 0.90,
        "in_frame_deletion": 0.75,
        "complex_indel": 0.85,
        "unknown": 0.50,
        "no_variant": None,
    }
    return presets.get(variant_type)


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():
    print("=" * 70)
    print("  CALON-Structure Phase 1: Variant Extraction & Classification")
    print("=" * 70)
    print()

    # ── Load Wales data ──────────────────────────────────────────────────
    print("Loading Wales data...")

    wales = pd.read_csv(WALES_FILE, low_memory=False)
    dragon3 = pd.read_csv(DRAGON3_FILE, low_memory=False)

    print(f"  All Wales: {len(wales)} patients")
    print(f"  DRAGON_3 (South Wales): {len(dragon3)} patients")

    # ── Extract all unique mutations ────────────────────────────────────
    print("\nExtracting mutations...")

    all_mutations = []

    # Wales: Mutation1, Mutation2, Mutation3 (columns 127, 142, 157 — 0-indexed: 126, 141, 156)
    wales_cols = [c for c in wales.columns if c.startswith("Mutation")]
    dragon3_cols = [c for c in dragon3.columns if c.startswith("Mutation")]

    print(f"  Wales mutation columns: {wales_cols}")
    print(f"  DRAGON_3 mutation columns: {dragon3_cols}")

    # Collect from Wales
    for col in wales_cols:
        vals = wales[col].dropna().unique()
        for v in vals:
            v_str = str(v).strip()
            if v_str and v_str != "nan":
                all_mutations.append(v_str)

    # Collect from DRAGON_3
    for col in dragon3_cols:
        vals = dragon3[col].dropna().unique()
        for v in vals:
            v_str = str(v).strip()
            if v_str and v_str != "nan":
                all_mutations.append(v_str)

    # Deduplicate
    unique_mutations = sorted(set(all_mutations))
    print(f"  Total unique mutation strings: {len(unique_mutations)}")

    # ── Parse and classify each variant ────────────────────────────────
    print("\nParsing and classifying variants...")

    results = []
    for mut_str in unique_mutations:
        parsed = parse_hgvs_cdna(mut_str)
        if parsed is None:
            continue

        vtype = classify_variant(parsed)
        gene = parsed.get("gene", "UNKNOWN")
        cdna = parsed.get("cdna_change", "")

        # Estimate protein position (for LDLR domain mapping)
        prot_pos = estimate_protein_position(cdna)

        # Domain mapping (LDLR only)
        if gene == "LDLR":
            domain, domain_severity = get_ldlr_domain(prot_pos)
        elif gene == "APOB":
            domain = "ApoB-RBD" if prot_pos and 3300 <= prot_pos <= 3600 else "ApoB-other"
            domain_severity = 0.8 if domain == "ApoB-RBD" else 0.4
        elif gene == "PCSK9":
            domain = "PCSK9-catalytic"
            domain_severity = 0.7
        else:
            domain = "Unknown"
            domain_severity = 0.5

        # SSS preset for non-missense
        sss_preset = assign_sss_preset(vtype)
        needs_af3 = vtype == "substitution"  # Only substitutions need AF3 modelling

        # Refine substitution → check if it's truly missense or nonsense
        # (requires codon context — for now flag as "substitution, needs_refinement")

        results.append({
            "variant_id": mut_str,
            "gene": gene,
            "cdna_change": cdna,
            "protein_position": prot_pos,
            "variant_type": vtype,
            "domain": domain,
            "domain_severity": domain_severity,
            "sss_preset": sss_preset,
            "needs_af3_modelling": needs_af3,
            "parse_failed": parsed.get("parse_failed", False),
        })

    catalogue = pd.DataFrame(results)
    print(f"  Parsed variants: {len(catalogue)}")

    # ── Summary statistics ─────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("  VARIANT CATALOGUE SUMMARY")
    print("=" * 70)

    print(f"\n  By gene:")
    for gene, count in catalogue["gene"].value_counts().items():
        print(f"    {gene}: {count}")

    print(f"\n  By variant type:")
    for vtype, count in catalogue["variant_type"].value_counts().items():
        print(f"    {vtype}: {count}")

    print(f"\n  By LDLR domain (LDLR variants only):")
    ldlr = catalogue[catalogue["gene"] == "LDLR"]
    for domain, count in ldlr["domain"].value_counts().items():
        print(f"    {domain}: {count}")

    print(f"\n  Needs AF3 modelling (substitutions): {catalogue['needs_af3_modelling'].sum()}")
    print(f"  Pre-assigned SSS (non-missense): {catalogue['sss_preset'].notna().sum()}")
    print(f"  Parse failures: {catalogue['parse_failed'].sum()}")

    # ── Create patient-variant mapping ──────────────────────────────────
    print("\nCreating patient-variant mapping...")

    patient_variants = []

    # From Wales
    for _, row in wales.iterrows():
        for col in wales_cols:
            val = row.get(col)
            if pd.notna(val) and str(val).strip():
                patient_variants.append({
                    "source": "Wales",
                    "patient_id": row.get(wales.columns[0], None),  # First col = ID
                    "mutation_column": col,
                    "variant_id": str(val).strip(),
                })

    # From DRAGON_3
    for _, row in dragon3.iterrows():
        for col in dragon3_cols:
            val = row.get(col)
            if pd.notna(val) and str(val).strip():
                patient_variants.append({
                    "source": "DRAGON_3",
                    "patient_id": row.get(dragon3.columns[0], None),
                    "mutation_column": col,
                    "variant_id": str(val).strip(),
                })

    patient_map = pd.DataFrame(patient_variants)
    print(f"  Patient-variant pairs: {len(patient_map)}")

    # Count patients per variant (for prioritising AF3 modelling)
    variant_counts = patient_map["variant_id"].value_counts().reset_index()
    variant_counts.columns = ["variant_id", "patient_count"]

    catalogue = catalogue.merge(variant_counts, on="variant_id", how="left")
    catalogue["patient_count"] = catalogue["patient_count"].fillna(0).astype(int)

    # ── Select top missense variants for AF3 modelling ──────────────────
    missense = catalogue[catalogue["needs_af3_modelling"]].copy()
    missense = missense.sort_values("patient_count", ascending=False)

    # Top 50 cover most patients
    top50 = missense.head(50)
    patients_covered = top50["patient_count"].sum()
    total_missense_patients = missense["patient_count"].sum()

    print(f"\n  Top 50 missense variants cover {patients_covered}/{total_missense_patients} "
          f"patients ({100*patients_covered/max(total_missense_patients,1):.1f}%)")

    # ── Save outputs ──────────────────────────────────────────────────
    print("\nSaving outputs...")

    # Full catalogue
    cat_path = os.path.join(OUTPUT_DIR, "af3_variant_catalogue.csv")
    catalogue.to_csv(cat_path, index=False)
    print(f"  Saved: {cat_path} ({len(catalogue)} variants)")

    # Missense for modelling
    mis_path = os.path.join(OUTPUT_DIR, "af3_missense_for_modelling.csv")
    missense.to_csv(mis_path, index=False)
    print(f"  Saved: {mis_path} ({len(missense)} missense variants)")

    # Top 50 priority list
    top_path = os.path.join(OUTPUT_DIR, "af3_top50_missense.csv")
    top50.to_csv(top_path, index=False)
    print(f"  Saved: {top_path} (top 50 for AF3/FoldX)")

    # Patient-variant map
    map_path = os.path.join(OUTPUT_DIR, "af3_variant_patient_map.csv")
    patient_map.to_csv(map_path, index=False)
    print(f"  Saved: {map_path} ({len(patient_map)} patient-variant pairs)")

    # ── Print AF3 Server submission list ─────────────────────────────
    print("\n" + "=" * 70)
    print("  AF3 SERVER SUBMISSION LIST")
    print("=" * 70)
    print("\n  Wild-type structures to submit to alphafoldserver.com:")
    print("  1. LDLR (UniProt P01130, residues 1-860)")
    print("  2. PCSK9 (UniProt Q8NBP7, residues 1-692)")
    print("  3. ApoB-100 RBD (UniProt P04114, residues 3359-3600)")
    print("  4. LDLR + PCSK9 complex (multimer)")
    print("  5. LDLR + ApoB-100 RBD complex (multimer)")

    print("\n  Top 10 missense variants for FoldX modelling:")
    for i, row in top50.head(10).iterrows():
        print(f"    {row['variant_id']} (n={row['patient_count']} patients, "
              f"domain={row['domain']}, pos={row['protein_position']})")

    print("\n" + "=" * 70)
    print("  PHASE 1 COMPLETE")
    print("  Next: Submit AF3 jobs + apply for FoldX license")
    print("=" * 70)


if __name__ == "__main__":
    main()

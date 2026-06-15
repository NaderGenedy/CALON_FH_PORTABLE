"""
Download Findlay et al. (2018) BRCA1 Saturation Genome Editing functional scores
from MaveDB (urn:mavedb:00000097-0-2) and process into analysis-ready CSV.

Reference: Findlay GM et al. "Accurate classification of BRCA1 variants with
saturation genome editing" Nature 2018. doi:10.1038/s41586-018-0461-z

MaveDB accession: urn:mavedb:00000097-0-2 (normalised meta-analysis scores)
3,893 SNVs across 13 exons (2-5, 15-22) of BRCA1 RING and BRCT domains.
"""

import requests
import csv
import re
import sys
import os
import time

# --- Configuration ---
SCORESET_URN = "urn:mavedb:00000097-0-2"
API_BASE = "https://api.mavedb.org/api/v1"
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "brca1_functional_scores.csv")

# BRCA1 NM_007294.3 CDS: 5592 nt coding for 1863 aa
# Codon table
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Complementary bases
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def fetch_brca1_cds():
    """Fetch BRCA1 CDS from NCBI for NM_007294.3 via Entrez."""
    print("Fetching BRCA1 CDS sequence from NCBI (NM_007294.3)...")
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nucleotide",
        "id": "NM_007294.3",
        "rettype": "fasta_cds_na",
        "retmode": "text"
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()

    # Parse FASTA
    lines = resp.text.strip().split('\n')
    seq_lines = [l.strip() for l in lines if not l.startswith('>')]
    cds = ''.join(seq_lines).upper()
    print(f"  CDS length: {len(cds)} nt ({len(cds)//3} codons)")
    return cds


def fetch_scores_from_mavedb():
    """Download all scores from MaveDB API (returns CSV format)."""
    scores_url = f"{API_BASE}/score-sets/{SCORESET_URN}/scores/"
    print(f"Fetching scores from MaveDB ({SCORESET_URN})...")

    resp = requests.get(scores_url, timeout=120)
    resp.raise_for_status()

    # API returns CSV directly
    import io
    reader = csv.DictReader(io.StringIO(resp.text))
    all_scores = list(reader)

    print(f"  Total variants downloaded: {len(all_scores)}")
    return all_scores


def parse_hgvs_nt(hgvs):
    """Parse HGVS nucleotide notation like NM_007294.3:c.5565A>T.
    Returns (cds_position, ref_base, alt_base) or None."""
    if not hgvs or hgvs == 'NA':
        return None

    m = re.match(r'NM_007294\.\d+:c\.(\d+)([ACGT])>([ACGT])', hgvs)
    if m:
        pos = int(m.group(1))
        ref = m.group(2)
        alt = m.group(3)
        return (pos, ref, alt)
    return None


def cds_to_protein_change(cds_pos, ref_base, alt_base, cds_seq):
    """Convert CDS position + substitution to protein change.
    CDS positions are 1-based. Returns (aa_pos, wt_aa, mut_aa) or None."""
    if cds_pos < 1 or cds_pos > len(cds_seq):
        return None

    # Which codon? (1-based)
    codon_num = (cds_pos - 1) // 3 + 1  # 1-based protein position
    codon_start = (codon_num - 1) * 3    # 0-based index into CDS

    if codon_start + 3 > len(cds_seq):
        return None

    wt_codon = cds_seq[codon_start:codon_start + 3]

    # Position within codon (0-based)
    pos_in_codon = (cds_pos - 1) % 3

    # Verify reference base matches
    if wt_codon[pos_in_codon] != ref_base:
        # Try complement (shouldn't happen for CDS notation but just in case)
        return None

    # Create mutant codon
    mut_codon = list(wt_codon)
    mut_codon[pos_in_codon] = alt_base
    mut_codon = ''.join(mut_codon)

    wt_aa = CODON_TABLE.get(wt_codon, '?')
    mut_aa = CODON_TABLE.get(mut_codon, '?')

    return (codon_num, wt_aa, mut_aa)


def classify_function(score):
    """Classify variant function based on Findlay et al. thresholds.
    From the paper:
    - Functional: score >= -1.328
    - LOF: score <= -2.114
    - Intermediate: between these thresholds
    """
    if score is None:
        return "NA"
    if score >= -1.328:
        return "FUNC"
    elif score <= -2.114:
        return "LOF"
    else:
        return "INT"


def determine_variant_type(wt_aa, mut_aa):
    """Determine variant consequence type."""
    if wt_aa == mut_aa:
        return "synonymous"
    elif mut_aa == '*':
        return "nonsense"
    elif wt_aa == '*':
        return "stop_loss"
    else:
        return "missense"


def main():
    # Step 1: Get CDS sequence
    try:
        cds_seq = fetch_brca1_cds()
    except Exception as e:
        print(f"Failed to fetch CDS from NCBI: {e}")
        print("Using hardcoded CDS length for position mapping...")
        cds_seq = None

    # Step 2: Fetch all scores from MaveDB
    raw_scores = fetch_scores_from_mavedb()

    if not raw_scores:
        print("ERROR: No scores downloaded from MaveDB!")
        sys.exit(1)

    # Step 3: Process each variant
    processed = []
    skipped = 0

    for variant in raw_scores:
        hgvs_nt = variant.get('hgvs_nt', '')
        hgvs_pro = variant.get('hgvs_pro', '')
        score_str = variant.get('score', '')
        score_rep1 = variant.get('score_rep1', '')
        score_rep2 = variant.get('score_rep2', '')
        score_rna = variant.get('score_rna', '')

        if not score_str or score_str == 'NA' or score_str == 'None':
            skipped += 1
            continue

        try:
            score = float(score_str)
        except ValueError:
            skipped += 1
            continue

        # Parse nucleotide change
        nt_parsed = parse_hgvs_nt(hgvs_nt)
        if nt_parsed is None:
            # Check if it's an intronic/splice variant (c.NNN+/-NNN)
            m_splice = re.match(r'NM_007294\.\d+:c\.(\d+)[+\-]', hgvs_nt or '')
            if m_splice:
                nearest_cds = int(m_splice.group(1))
                aa_pos = (nearest_cds - 1) // 3 + 1
                row = {
                    'hgvs_nt': hgvs_nt,
                    'hgvs_pro': 'NA',
                    'cds_position': nearest_cds,
                    'protein_position': aa_pos,
                    'wt_aa': 'NA',
                    'mut_aa': 'NA',
                    'variant_type': 'splice_region',
                    'function_score': round(score, 4),
                    'function_score_rep1': round(float(score_rep1), 4) if score_rep1 and score_rep1 not in ('NA', 'None', '') else 'NA',
                    'function_score_rep2': round(float(score_rep2), 4) if score_rep2 and score_rep2 not in ('NA', 'None', '') else 'NA',
                    'function_score_rna': round(float(score_rna), 4) if score_rna and score_rna not in ('NA', 'None', '') else 'NA',
                    'function_class': classify_function(score),
                }
                processed.append(row)
                continue
            skipped += 1
            continue

        cds_pos, ref_base, alt_base = nt_parsed

        # Convert to protein
        if cds_seq:
            prot = cds_to_protein_change(cds_pos, ref_base, alt_base, cds_seq)
        else:
            prot = None

        if prot:
            aa_pos, wt_aa, mut_aa = prot
            var_type = determine_variant_type(wt_aa, mut_aa)
        else:
            # Fallback: calculate position from CDS pos
            aa_pos = (cds_pos - 1) // 3 + 1
            wt_aa = "?"
            mut_aa = "?"
            var_type = "unknown"

        func_class = classify_function(score)

        row = {
            'hgvs_nt': hgvs_nt if hgvs_nt else 'NA',
            'hgvs_pro': hgvs_pro if hgvs_pro and hgvs_pro != 'NA' else f"p.{wt_aa}{aa_pos}{mut_aa}",
            'cds_position': cds_pos,
            'protein_position': aa_pos,
            'wt_aa': wt_aa,
            'mut_aa': mut_aa,
            'variant_type': var_type,
            'function_score': round(score, 4),
            'function_score_rep1': round(float(score_rep1), 4) if score_rep1 and score_rep1 not in ('NA', 'None', '') else 'NA',
            'function_score_rep2': round(float(score_rep2), 4) if score_rep2 and score_rep2 not in ('NA', 'None', '') else 'NA',
            'function_score_rna': round(float(score_rna), 4) if score_rna and score_rna not in ('NA', 'None', '') else 'NA',
            'function_class': func_class,
        }
        processed.append(row)

    # Sort by protein position, then by mutant AA
    processed.sort(key=lambda x: (x['protein_position'], x['mut_aa']))

    # Step 4: Write CSV
    fieldnames = [
        'hgvs_nt', 'hgvs_pro', 'cds_position', 'protein_position',
        'wt_aa', 'mut_aa', 'variant_type', 'function_score',
        'function_score_rep1', 'function_score_rep2', 'function_score_rna',
        'function_class'
    ]

    with open(OUTPUT_FILE, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(processed)

    # Summary stats
    n_total = len(processed)
    n_func = sum(1 for r in processed if r['function_class'] == 'FUNC')
    n_lof = sum(1 for r in processed if r['function_class'] == 'LOF')
    n_int = sum(1 for r in processed if r['function_class'] == 'INT')
    n_missense = sum(1 for r in processed if r['variant_type'] == 'missense')
    n_syn = sum(1 for r in processed if r['variant_type'] == 'synonymous')
    n_nonsense = sum(1 for r in processed if r['variant_type'] == 'nonsense')

    positions = set(r['protein_position'] for r in processed)

    print(f"\n{'='*60}")
    print(f"BRCA1 SGE Functional Scores - Findlay et al. (2018)")
    print(f"{'='*60}")
    print(f"Output: {OUTPUT_FILE}")
    print(f"Total variants:    {n_total}")
    print(f"Skipped:           {skipped}")
    print(f"Protein positions: {len(positions)} (range {min(positions)}-{max(positions)})")
    print(f"")
    print(f"By function class:")
    print(f"  Functional:      {n_func} ({100*n_func/n_total:.1f}%)")
    print(f"  Intermediate:    {n_int} ({100*n_int/n_total:.1f}%)")
    print(f"  LOF:             {n_lof} ({100*n_lof/n_total:.1f}%)")
    print(f"")
    print(f"By variant type:")
    print(f"  Missense:        {n_missense}")
    print(f"  Synonymous:      {n_syn}")
    print(f"  Nonsense:        {n_nonsense}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

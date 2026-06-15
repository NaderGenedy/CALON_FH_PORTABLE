#!/usr/bin/env python3
"""
Download AlphaMissense predictions for LDLR (P01130) from Zenodo.
Streams the gzipped TSV to avoid downloading the full 1.2GB file to disk.
Filters for LDLR UniProt ID (P01130) and saves to CSV.
"""

import gzip
import io
import csv
import sys
import time
import urllib.request

# AlphaMissense aa_substitutions file from Google Cloud Storage
URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz"
OUTPUT = "C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis/alphamissense_ldlr.csv"

LDLR_UNIPROT = "P01130"

def download_and_filter():
    print(f"Downloading AlphaMissense aa_substitutions from Zenodo...")
    print(f"URL: {URL}")
    print(f"This file is ~1.2 GB compressed. Streaming and filtering for LDLR (P01130)...")
    print()

    start_time = time.time()

    # Download the file
    req = urllib.request.Request(URL, headers={'User-Agent': 'Mozilla/5.0'})

    ldlr_rows = []
    total_lines = 0
    header = None

    try:
        with urllib.request.urlopen(req, timeout=600) as response:
            # Read entire response into memory then decompress
            # (streaming gzip over HTTP is unreliable)
            print("Downloading file... (this may take a few minutes)")
            data = response.read()
            elapsed = time.time() - start_time
            size_mb = len(data) / (1024 * 1024)
            print(f"Downloaded {size_mb:.1f} MB in {elapsed:.0f} seconds")
            print()

            print("Decompressing and filtering for LDLR...")
            with gzip.open(io.BytesIO(data), 'rt') as f:
                for line in f:
                    # Skip comment lines
                    if line.startswith('#'):
                        continue

                    total_lines += 1

                    # First non-comment line is header
                    if header is None:
                        header = line.strip().split('\t')
                        print(f"Header columns: {header}")
                        continue

                    # Check if this row is for LDLR
                    fields = line.strip().split('\t')

                    # The file format has columns:
                    # uniprot_id, protein_variant, am_pathogenicity, am_class
                    if len(fields) >= 1 and fields[0] == LDLR_UNIPROT:
                        ldlr_rows.append(fields)

                    if total_lines % 5000000 == 0:
                        print(f"  Processed {total_lines:,} lines, found {len(ldlr_rows)} LDLR variants so far...")

    except Exception as e:
        print(f"Error during download: {e}")
        print("Trying alternative approach...")
        return False

    elapsed = time.time() - start_time
    print(f"\nProcessed {total_lines:,} total lines in {elapsed:.0f} seconds")
    print(f"Found {len(ldlr_rows)} LDLR variants")

    if len(ldlr_rows) == 0:
        print("ERROR: No LDLR rows found. Checking file format...")
        return False

    # Parse protein_variant into position, wt_aa, mut_aa
    # Format is typically: A1B (single letter AA, position, single letter AA)
    print("\nParsing variant annotations...")

    output_rows = []
    for fields in ldlr_rows:
        uniprot_id = fields[0]
        protein_variant = fields[1]  # e.g., "M1A"
        am_score = float(fields[2])
        am_class = fields[3]

        # Parse variant: first char = wt_aa, last char = mut_aa, middle = position
        wt_aa = protein_variant[0]
        mut_aa = protein_variant[-1]
        position = int(protein_variant[1:-1])

        output_rows.append({
            'position': position,
            'wt_aa': wt_aa,
            'mut_aa': mut_aa,
            'am_score': am_score,
            'am_class': am_class,
            'protein_variant': protein_variant,
            'uniprot_id': uniprot_id
        })

    # Sort by position and mutant amino acid
    output_rows.sort(key=lambda x: (x['position'], x['mut_aa']))

    # Write CSV
    print(f"Writing {len(output_rows)} rows to {OUTPUT}")
    with open(OUTPUT, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['position', 'wt_aa', 'mut_aa', 'am_score', 'am_class', 'protein_variant', 'uniprot_id'])
        writer.writeheader()
        writer.writerows(output_rows)

    # Print summary statistics
    pathogenic = sum(1 for r in output_rows if r['am_class'] == 'likely_pathogenic')
    ambiguous = sum(1 for r in output_rows if r['am_class'] == 'ambiguous')
    benign = sum(1 for r in output_rows if r['am_class'] == 'likely_benign')

    print(f"\n=== Summary ===")
    print(f"Total LDLR variants: {len(output_rows)}")
    print(f"Likely pathogenic: {pathogenic} ({100*pathogenic/len(output_rows):.1f}%)")
    print(f"Ambiguous: {ambiguous} ({100*ambiguous/len(output_rows):.1f}%)")
    print(f"Likely benign: {benign} ({100*benign/len(output_rows):.1f}%)")

    positions = set(r['position'] for r in output_rows)
    print(f"Positions covered: {min(positions)}-{max(positions)} ({len(positions)} unique)")

    scores = [r['am_score'] for r in output_rows]
    print(f"Score range: {min(scores):.4f} - {max(scores):.4f}")
    print(f"Mean score: {sum(scores)/len(scores):.4f}")

    print(f"\nFile saved to: {OUTPUT}")
    return True


if __name__ == '__main__':
    success = download_and_filter()
    if not success:
        sys.exit(1)

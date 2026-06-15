"""
Fetch all LDLR variants from ClinVar via NCBI E-utilities API.
Parses missense variants and saves to CSV for FoldX integration.
"""

import urllib.request
import urllib.parse
import json
import csv
import re
import time
import sys

# Amino acid 3-letter to 1-letter mapping
AA_MAP = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Ter': '*', 'Xaa': 'X'
}

# Regex for missense: p.Ala123Val or p.A123V
MISSENSE_3LETTER = re.compile(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})')
MISSENSE_1LETTER = re.compile(r'p\.([A-Z])(\d+)([A-Z])')

OUTPUT_CSV = r"C:\Users\nader\Downloads\calon_ukb_pipeline\alphafold\analysis\clinvar_ldlr_full.csv"


def fetch_url(url, retries=3, delay=2):
    """Fetch URL with retries."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'ClinVar-LDLR-Fetcher/1.0'})
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read()
        except Exception as e:
            print(f"  Attempt {attempt+1} failed: {e}")
            if attempt < retries - 1:
                time.sleep(delay * (attempt + 1))
    return None


def search_clinvar_ids():
    """Search ClinVar for all LDLR variants, return list of variant IDs."""
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        'db': 'clinvar',
        'term': 'LDLR[gene]',
        'retmax': 0,
        'rettype': 'uilist',
        'retmode': 'json'
    }
    url = f"{base}?{urllib.parse.urlencode(params)}"
    print("Fetching variant count...")
    data = fetch_url(url)
    if not data:
        print("ERROR: Could not fetch count")
        return []

    result = json.loads(data)
    count = int(result['esearchresult']['count'])
    print(f"Found {count} LDLR variants in ClinVar")

    all_ids = []
    batch_size = 5000
    for start in range(0, count, batch_size):
        params['retstart'] = start
        params['retmax'] = batch_size
        url = f"{base}?{urllib.parse.urlencode(params)}"
        print(f"  Fetching IDs {start} to {start + batch_size}...")
        data = fetch_url(url)
        if data:
            result = json.loads(data)
            ids = result['esearchresult']['idlist']
            all_ids.extend(ids)
            print(f"  Got {len(ids)} IDs (total: {len(all_ids)})")
        time.sleep(0.5)

    return all_ids


def fetch_variant_summaries(ids):
    """Fetch variant summaries using esummary in batches."""
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    variants = []
    batch_size = 80  # Small enough to avoid URI too long

    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        id_str = ','.join(batch)
        url = f"{base}?db=clinvar&id={id_str}&retmode=json"
        progress = f"[{i}/{len(ids)}]"
        if i % 400 == 0:
            print(f"  {progress} Fetching summaries...")
        data = fetch_url(url)
        if not data:
            print(f"  WARNING: Failed batch at {i}")
            continue

        try:
            result = json.loads(data)
            doc_sums = result.get('result', {})
            uid_list = doc_sums.get('uids', [])

            for uid in uid_list:
                entry = doc_sums.get(uid, {})
                if not entry or 'error' in entry:
                    continue

                title = entry.get('title', '')

                # Clinical significance from germline_classification
                germ = entry.get('germline_classification', {})
                clin_sig = germ.get('description', '') if isinstance(germ, dict) else ''
                review_status = germ.get('review_status', '') if isinstance(germ, dict) else ''

                # Conditions from germline_classification.trait_set
                conditions = []
                trait_set = germ.get('trait_set', []) if isinstance(germ, dict) else []
                if isinstance(trait_set, list):
                    for trait in trait_set:
                        if isinstance(trait, dict):
                            tname = trait.get('trait_name', '')
                            if tname:
                                conditions.append(tname)
                condition_str = '; '.join(conditions)

                # Protein change field (comma-separated isoforms)
                protein_change = entry.get('protein_change', '')

                # Molecular consequence
                mol_consequence = entry.get('molecular_consequence_list', [])
                mol_type = ', '.join(mol_consequence) if isinstance(mol_consequence, list) else str(mol_consequence)

                variants.append({
                    'clinvar_id': uid,
                    'accession': entry.get('accession', ''),
                    'title': title,
                    'protein_change': protein_change,
                    'mol_consequence': mol_type,
                    'clinical_significance': clin_sig,
                    'review_status': review_status,
                    'condition': condition_str,
                })
        except Exception as e:
            print(f"  ERROR parsing batch at {i}: {e}")

        time.sleep(0.35)

    return variants


def parse_protein_change(protein_change_field, title):
    """Extract position, wt_aa, mut_aa from protein change or title.

    The protein_change field can have multiple isoforms like:
    'D289E, D330E, D416E, D457E'
    We want the canonical (usually longest/last) but let's take all unique.

    The title has format: NM_000527.5(LDLR):c.1371C>G (p.Asp457Glu)
    """
    results = []

    # Try title first for 3-letter code
    if title:
        m = MISSENSE_3LETTER.search(title)
        if m:
            wt = AA_MAP.get(m.group(1))
            pos = int(m.group(2))
            mut = AA_MAP.get(m.group(3))
            if wt and mut and wt != '*' and mut != '*':
                results.append((pos, wt, mut))
                return results

    # Parse protein_change field (1-letter codes like D457E)
    if protein_change_field:
        # Pattern: single letter + digits + single letter
        for match in re.finditer(r'([A-Z])(\d+)([A-Z])', protein_change_field):
            wt = match.group(1)
            pos = int(match.group(2))
            mut = match.group(3)
            if wt != '*' and mut != '*' and wt != mut:
                results.append((pos, wt, mut))

    return results


def main():
    print("=" * 60)
    print("ClinVar LDLR Variant Fetcher")
    print("=" * 60)

    # Step 1: Get all variant IDs
    ids = search_clinvar_ids()
    if not ids:
        print("No variant IDs found. Exiting.")
        sys.exit(1)
    print(f"\nTotal variant IDs: {len(ids)}")

    # Step 2: Fetch summaries
    print("\nFetching variant summaries...")
    variants = fetch_variant_summaries(ids)
    print(f"\nTotal variants fetched: {len(variants)}")

    # Step 3: Parse and write
    fieldnames = ['variant_name', 'position', 'wt_aa', 'mut_aa',
                  'clinical_significance', 'review_status', 'condition',
                  'clinvar_id', 'accession', 'molecular_consequence']

    all_rows = []
    missense_rows = []

    for v in variants:
        parsed = parse_protein_change(v['protein_change'], v['title'])

        if parsed:
            # Use the last (canonical) protein change
            # For LDLR NM_000527.5, the mature protein starts at signal peptide pos 1
            # Take the largest position number as canonical
            pos, wt, mut = max(parsed, key=lambda x: x[0])
            row = {
                'variant_name': v['title'],
                'position': pos,
                'wt_aa': wt,
                'mut_aa': mut,
                'clinical_significance': v['clinical_significance'],
                'review_status': v['review_status'],
                'condition': v['condition'],
                'clinvar_id': v['clinvar_id'],
                'accession': v['accession'],
                'molecular_consequence': v['mol_consequence']
            }
            all_rows.append(row)
            missense_rows.append(row)
        else:
            row = {
                'variant_name': v['title'],
                'position': '',
                'wt_aa': '',
                'mut_aa': '',
                'clinical_significance': v['clinical_significance'],
                'review_status': v['review_status'],
                'condition': v['condition'],
                'clinvar_id': v['clinvar_id'],
                'accession': v['accession'],
                'molecular_consequence': v['mol_consequence']
            }
            all_rows.append(row)

    # Write full CSV (all variants)
    with open(OUTPUT_CSV, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nSaved {len(all_rows)} total variants to: {OUTPUT_CSV}")

    # Write missense-only CSV
    missense_csv = OUTPUT_CSV.replace('.csv', '_missense.csv')
    with open(missense_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(missense_rows)
    print(f"Saved {len(missense_rows)} missense variants to: {missense_csv}")

    # Summary stats
    if missense_rows:
        sig_counts = {}
        for r in missense_rows:
            sig = r['clinical_significance'] or '(not provided)'
            sig_counts[sig] = sig_counts.get(sig, 0) + 1

        print("\n--- Missense Clinical Significance Distribution ---")
        for sig, count in sorted(sig_counts.items(), key=lambda x: -x[1]):
            print(f"  {sig}: {count}")

        positions = [r['position'] for r in missense_rows]
        print(f"\n  Position range: {min(positions)} - {max(positions)}")
        print(f"  Unique positions: {len(set(positions))}")
        print(f"  Unique mutations: {len(set((r['position'], r['wt_aa'], r['mut_aa']) for r in missense_rows))}")

    # Non-missense type breakdown
    non_missense = [r for r in all_rows if not r['position']]
    if non_missense:
        type_counts = {}
        for r in non_missense:
            t = r['molecular_consequence'] or '(unknown)'
            type_counts[t] = type_counts.get(t, 0) + 1
        print(f"\n--- Non-missense variant types ({len(non_missense)} total) ---")
        for t, c in sorted(type_counts.items(), key=lambda x: -x[1])[:10]:
            print(f"  {t}: {c}")


if __name__ == '__main__':
    main()

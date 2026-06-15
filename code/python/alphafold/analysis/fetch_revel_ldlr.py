"""
Fetch REVEL scores for LDLR missense variants from myvariant.info API.
Uses dbNSFP data bundled in myvariant.info.
LDLR is on chromosome 19.

Output: revel_ldlr.csv with columns: position, wt_aa, mut_aa, revel_score
"""

import requests
import pandas as pd
import time
import sys
import json

OUTPUT = "C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis/revel_ldlr.csv"

# Three-letter to one-letter amino acid mapping
AA3TO1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
}

def fetch_myvariant_revel():
    """Fetch REVEL scores from myvariant.info for all LDLR variants."""
    base_url = "https://myvariant.info/v1/query"
    all_results = []
    batch_from = 0
    batch_size = 1000

    print("Fetching LDLR variants from myvariant.info (dbNSFP REVEL scores)...")

    while True:
        params = {
            "q": 'dbnsfp.genename:LDLR',
            "fields": "dbnsfp.revel.score,dbnsfp.aa.alt,dbnsfp.aa.ref,dbnsfp.aa.pos,dbnsfp.aapos,dbnsfp.aaref,dbnsfp.aaalt,dbnsfp.hgvsp,dbnsfp.uniprot.acc,_id",
            "size": batch_size,
            "from": batch_from,
        }

        try:
            r = requests.get(base_url, params=params, timeout=60)
            r.raise_for_status()
            data = r.json()
        except Exception as e:
            print(f"  Error at offset {batch_from}: {e}")
            break

        hits = data.get("hits", [])
        total = data.get("total", 0)
        print(f"  Fetched {batch_from} - {batch_from + len(hits)} of {total}")

        if not hits:
            break

        all_results.extend(hits)
        batch_from += batch_size

        if batch_from >= total:
            break

        time.sleep(0.5)  # Be polite to the API

    print(f"Total raw hits: {len(all_results)}")
    return all_results


def parse_results(hits):
    """Parse myvariant.info results into position/wt_aa/mut_aa/revel_score rows."""
    rows = []

    for hit in hits:
        dbnsfp = hit.get("dbnsfp", {})
        if not dbnsfp:
            continue

        # Get REVEL score
        revel_data = dbnsfp.get("revel", {})
        if isinstance(revel_data, list):
            # Take first entry if list
            revel_data = revel_data[0] if revel_data else {}

        revel_score = revel_data.get("score") if isinstance(revel_data, dict) else revel_data
        if revel_score is None or revel_score == "." or revel_score == "":
            continue

        # Handle case where score is a list
        if isinstance(revel_score, list):
            # Filter out None/empty and take first valid
            valid = [s for s in revel_score if s is not None and s != "." and s != ""]
            if not valid:
                continue
            revel_score = valid[0]

        try:
            revel_score = float(revel_score)
        except (ValueError, TypeError):
            continue

        # Get amino acid info - try multiple field paths
        aapos = dbnsfp.get("aapos")
        aaref = dbnsfp.get("aaref")
        aaalt = dbnsfp.get("aaalt")

        # Try aa sub-dict
        if aapos is None:
            aa_data = dbnsfp.get("aa", {})
            if isinstance(aa_data, dict):
                aapos = aa_data.get("pos")
                aaref = aaref or aa_data.get("ref")
                aaalt = aaalt or aa_data.get("alt")

        # Try parsing from hgvsp if still missing
        if aapos is None:
            hgvsp = dbnsfp.get("hgvsp")
            if hgvsp and isinstance(hgvsp, str):
                # Format like p.Cys352Tyr or NP_000518.1:p.C352Y
                import re
                # Try three-letter
                m = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvsp)
                if m:
                    aaref = AA3TO1.get(m.group(1))
                    aapos = int(m.group(2))
                    aaalt = AA3TO1.get(m.group(3))
                else:
                    # Try one-letter
                    m = re.search(r'p\.([A-Z])(\d+)([A-Z])', hgvsp)
                    if m:
                        aaref = m.group(1)
                        aapos = int(m.group(2))
                        aaalt = m.group(3)

        if aapos is None or aaref is None or aaalt is None:
            continue

        # Handle lists (multiple transcripts)
        if isinstance(aapos, list):
            aapos = aapos[0]
        if isinstance(aaref, list):
            aaref = aaref[0]
        if isinstance(aaalt, list):
            aaalt = aaalt[0]

        try:
            aapos = int(aapos)
        except (ValueError, TypeError):
            continue

        # Convert 3-letter to 1-letter if needed
        if len(str(aaref)) == 3:
            aaref = AA3TO1.get(aaref, aaref)
        if len(str(aaalt)) == 3:
            aaalt = AA3TO1.get(aaalt, aaalt)

        # Filter to mature LDLR protein range (signal peptide 1-21, mature 22-860)
        if aapos < 1 or aapos > 900:
            continue

        # Skip synonymous
        if aaref == aaalt:
            continue

        rows.append({
            "position": aapos,
            "wt_aa": str(aaref),
            "mut_aa": str(aaalt),
            "revel_score": revel_score,
        })

    return rows


def main():
    # Strategy 1: myvariant.info bulk query
    print("=" * 60)
    print("Strategy 1: myvariant.info bulk query for LDLR REVEL scores")
    print("=" * 60)

    hits = fetch_myvariant_revel()

    if hits:
        rows = parse_results(hits)
        if rows:
            df = pd.DataFrame(rows)
            # Remove duplicates (same position/wt/mut), keep highest REVEL score
            df = df.sort_values("revel_score", ascending=False).drop_duplicates(
                subset=["position", "wt_aa", "mut_aa"], keep="first"
            )
            df = df.sort_values(["position", "wt_aa", "mut_aa"]).reset_index(drop=True)

            print(f"\nParsed {len(df)} unique LDLR missense variants with REVEL scores")
            print(f"Position range: {df['position'].min()} - {df['position'].max()}")
            print(f"REVEL score range: {df['revel_score'].min():.4f} - {df['revel_score'].max():.4f}")
            print(f"REVEL score mean: {df['revel_score'].mean():.4f}")

            df.to_csv(OUTPUT, index=False)
            print(f"\nSaved to: {OUTPUT}")
            return True

    print("\nStrategy 1 yielded no results. Trying Strategy 2...")

    # Strategy 2: Ensembl VEP API per-variant
    # Read our saturation mutagenesis file to get all positions
    print("=" * 60)
    print("Strategy 2: Ensembl VEP REST API (batch queries)")
    print("=" * 60)

    try:
        sat_df = pd.read_csv(
            "C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis/foldx_saturation_mutagenesis.csv"
        )
        # Get unique position/wt_aa pairs
        positions = sat_df[["position", "wt_aa"]].drop_duplicates().sort_values("position")
        print(f"Found {len(positions)} unique positions from saturation mutagenesis file")
    except Exception as e:
        print(f"Could not read saturation mutagenesis file: {e}")
        return False

    # Use Ensembl VEP POST endpoint for batch queries
    # Format: LDLR p.Xxx123Yyy
    vep_url = "https://rest.ensembl.org/vep/human/hgvs"
    AA1TO3 = {v: k for k, v in AA3TO1.items()}

    all_rows = []
    all_aas = list("ACDEFGHIKLMNPQRSTVWY")

    # Build all HGVS notations
    hgvs_list = []
    for _, row in positions.iterrows():
        pos = int(row["position"])
        wt = row["wt_aa"]
        wt3 = AA1TO3.get(wt, wt)
        for mut in all_aas:
            if mut == wt:
                continue
            mut3 = AA1TO3.get(mut, mut)
            hgvs = f"NP_000518.1:p.{wt3}{pos}{mut3}"
            hgvs_list.append((pos, wt, mut, hgvs))

    print(f"Total variants to query: {len(hgvs_list)}")
    print("This will take a while with batch VEP queries...")

    # Batch in groups of 200 (Ensembl limit)
    batch_size = 200
    for i in range(0, len(hgvs_list), batch_size):
        batch = hgvs_list[i:i+batch_size]
        hgvs_ids = [x[3] for x in batch]

        try:
            r = requests.post(
                vep_url,
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                json={"hgvs_notations": hgvs_ids},
                timeout=120,
            )
            if r.status_code == 429:
                print("  Rate limited, waiting 10s...")
                time.sleep(10)
                r = requests.post(
                    vep_url,
                    headers={"Content-Type": "application/json", "Accept": "application/json"},
                    json={"hgvs_notations": hgvs_ids},
                    timeout=120,
                )
            r.raise_for_status()
            results = r.json()
        except Exception as e:
            print(f"  Error at batch {i}: {e}")
            time.sleep(5)
            continue

        # Parse results
        for result in results:
            input_hgvs = result.get("input", "")
            # Find matching entry
            match = None
            for pos, wt, mut, hgvs in batch:
                if hgvs == input_hgvs:
                    match = (pos, wt, mut)
                    break

            if not match:
                continue

            # Look for REVEL in transcript_consequences
            for tc in result.get("transcript_consequences", []):
                revel = tc.get("revel_score")
                if revel is not None:
                    all_rows.append({
                        "position": match[0],
                        "wt_aa": match[1],
                        "mut_aa": match[2],
                        "revel_score": float(revel),
                    })
                    break

        if (i // batch_size) % 10 == 0:
            print(f"  Processed {i + len(batch)}/{len(hgvs_list)} variants, found {len(all_rows)} REVEL scores")

        time.sleep(0.5)

    if all_rows:
        df = pd.DataFrame(all_rows)
        df = df.drop_duplicates(subset=["position", "wt_aa", "mut_aa"], keep="first")
        df = df.sort_values(["position", "wt_aa", "mut_aa"]).reset_index(drop=True)

        print(f"\nParsed {len(df)} unique LDLR missense variants with REVEL scores")
        print(f"Position range: {df['position'].min()} - {df['position'].max()}")
        print(f"REVEL score range: {df['revel_score'].min():.4f} - {df['revel_score'].max():.4f}")

        df.to_csv(OUTPUT, index=False)
        print(f"\nSaved to: {OUTPUT}")
        return True

    print("No REVEL scores found from either strategy.")
    return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

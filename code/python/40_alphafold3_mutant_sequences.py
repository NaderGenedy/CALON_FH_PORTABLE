#!/usr/bin/env python3
"""
===================================================================================
ALPHAFOLD3 MUTANT SEQUENCE GENERATOR
Generates mutant LDLR sequences + AF3 submission JSON for 52 priority positions

For each position, generates the most clinically relevant mutation.
For interface positions without known mutations, generates the most common
pathogenic substitution (based on ClinVar/literature patterns).

Output:
  - alphafold/analysis/af3_mutant_sequences/ (individual FASTA files)
  - alphafold/analysis/af3_mutant_jobs.csv (master job list)
  - alphafold/analysis/af3_submission_guide.txt (step-by-step guide)

Author: Dr Nader Genedy
===================================================================================
"""

import csv, json, os

BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"
AF = f"{BASE}/alphafold/analysis"
OUT = f"{AF}/af3_mutant_sequences"
os.makedirs(OUT, exist_ok=True)

# ─── LDLR wildtype sequence (860 residues, NP_000518.1) ───
# Load directly from AF3 job JSON
with open(f"{BASE}/alphafold/af3_job1_LDLR_wildtype.json", 'r') as f:
    job_data = json.load(f)
# AF3 JSON structure: list with one dict containing 'sequences' list
LDLR_WT = job_data[0]['sequences'][0]['proteinChain']['sequence']
print(f"LDLR wildtype sequence: {len(LDLR_WT)} residues")
assert len(LDLR_WT) == 860, f"Expected 860 residues, got {len(LDLR_WT)}"
# Verify key positions match atlas
assert LDLR_WT[45] == 'C', f"Position 46 should be C, got {LDLR_WT[45]}"
assert LDLR_WT[100] == 'E', f"Position 101 should be E, got {LDLR_WT[100]}"
assert LDLR_WT[405] == 'R', f"Position 406 should be R, got {LDLR_WT[405]}"
assert LDLR_WT[680] == 'C', f"Position 681 should be C, got {LDLR_WT[680]}"
print("  Sequence verified against atlas positions")

# ─── Standard amino acid substitution patterns ───
# For positions without known mutations, use the most pathogenic common substitution
COMMON_PATHOGENIC = {
    'C': 'Y',  # Cysteine -> Tyrosine (breaks disulphide bonds)
    'D': 'N',  # Aspartate -> Asparagine (charge loss)
    'E': 'K',  # Glutamate -> Lysine (charge reversal)
    'G': 'D',  # Glycine -> Aspartate (steric + charge)
    'R': 'W',  # Arginine -> Tryptophan (size + charge loss)
    'Q': 'P',  # Glutamine -> Proline (helix breaker)
    'L': 'P',  # Leucine -> Proline (helix breaker)
    'A': 'V',  # Alanine -> Valine (mild increase)
    'V': 'M',  # Valine -> Methionine (larger)
    'S': 'P',  # Serine -> Proline (helix breaker)
    'T': 'I',  # Threonine -> Isoleucine (polarity loss)
    'W': 'R',  # Tryptophan -> Arginine (size + polarity change)
    'P': 'L',  # Proline -> Leucine (flexibility change)
    'H': 'R',  # Histidine -> Arginine (charge at neutral pH)
    'M': 'I',  # Methionine -> Isoleucine (start codon loss)
    'N': 'D',  # Asparagine -> Aspartate (charge gain)
    'K': 'E',  # Lysine -> Glutamate (charge reversal)
    'F': 'S',  # Phenylalanine -> Serine (size + polarity)
    'I': 'T',  # Isoleucine -> Threonine (polarity)
    'Y': 'C',  # Tyrosine -> Cysteine (size change)
}

# ─── Load target positions ───
targets = []
with open(f"{AF}/alphafold3_mutant_targets.csv", encoding='utf-8-sig') as f:
    reader = csv.DictReader(f)
    for r in reader:
        pos = int(r['position'])
        wt = r['wildtype_aa']
        suggested = r.get('suggested_mutation', '').strip()
        domain = r['domain']
        reason = r.get('priority_reason', '')
        plddt = r.get('plddt', '')
        ddg = r.get('foldx_max_ddG', '')
        n_patients = int(r.get('n_patients', 0))

        # Determine the mutation to model
        if suggested and len(suggested) >= 3:
            # Parse suggested mutation format: e.g., "C46G" or "C46Y"
            # Extract the mutant amino acid (last character)
            mut_aa = suggested[-1]
            if not mut_aa.isalpha() or mut_aa == wt:
                mut_aa = COMMON_PATHOGENIC.get(wt, 'A')
        elif suggested and '=' in suggested:
            # Format: "M1L=0.01"
            parts = suggested.split('=')
            mut_part = parts[0].strip()
            mut_aa = mut_part[-1] if mut_part else COMMON_PATHOGENIC.get(wt, 'A')
        else:
            mut_aa = COMMON_PATHOGENIC.get(wt, 'A')

        targets.append({
            'position': pos,
            'wt_aa': wt,
            'mut_aa': mut_aa,
            'domain': domain,
            'reason': reason,
            'plddt': plddt,
            'ddg': ddg,
            'n_patients': n_patients,
            'mutation_name': f"{wt}{pos}{mut_aa}",
        })

print(f"\nGenerating {len(targets)} mutant sequences...")

# ─── Generate mutant sequences and FASTA files ───
jobs = []
batch_size = 5  # AF3 server allows ~5 concurrent jobs
batch_num = 1

for i, t in enumerate(targets):
    pos = t['position']
    wt = t['wt_aa']
    mut = t['mut_aa']
    name = t['mutation_name']

    # Verify wildtype matches
    if pos <= len(LDLR_WT):
        actual_wt = LDLR_WT[pos - 1]  # 0-indexed
        if actual_wt != wt:
            print(f"  WARNING: Position {pos} expected {wt} but found {actual_wt} in sequence")
            wt = actual_wt
            name = f"{wt}{pos}{mut}"
            t['mutation_name'] = name
            t['wt_aa'] = wt

    # Create mutant sequence
    mutant_seq = LDLR_WT[:pos-1] + mut + LDLR_WT[pos:]

    # Verify mutation was applied
    assert mutant_seq[pos-1] == mut, f"Mutation at {pos} failed"
    assert len(mutant_seq) == len(LDLR_WT), f"Length mismatch at {pos}"

    # Write FASTA
    fasta_name = f"LDLR_{name}.fasta"
    fasta_path = f"{OUT}/{fasta_name}"
    with open(fasta_path, 'w') as f:
        f.write(f">LDLR_mutant_{name} | Position {pos} | Domain: {t['domain']} | {t['reason']}\n")
        # Write in 80-char lines
        for j in range(0, len(mutant_seq), 80):
            f.write(mutant_seq[j:j+80] + '\n')

    # Write AF3 submission JSON
    json_name = f"af3_job_{name}.json"
    json_path = f"{OUT}/{json_name}"

    # AlphaFold3 server input format
    af3_job = [
        {
            "proteinChain": {
                "sequence": mutant_seq,
                "count": 1
            }
        }
    ]

    with open(json_path, 'w') as f:
        json.dump(af3_job, f, indent=2)

    # Track batch
    batch = (i // batch_size) + 1

    jobs.append({
        'batch': batch,
        'job_name': f"LDLR_{name}",
        'mutation': name,
        'position': pos,
        'wt_aa': wt,
        'mut_aa': mut,
        'domain': t['domain'],
        'priority_reason': t['reason'],
        'plddt': t['plddt'],
        'foldx_ddG': t['ddg'],
        'n_patients': t['n_patients'],
        'fasta_file': fasta_name,
        'json_file': json_name,
        'status': 'pending',
    })

# ─── Write master job list ───
job_csv = f"{AF}/af3_mutant_jobs.csv"
with open(job_csv, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=jobs[0].keys())
    writer.writeheader()
    writer.writerows(jobs)

# ─── Also generate complex prediction jobs for key mutations ───
# These model the mutant LDLR with PCSK9 to see if binding interface changes
print("\nGenerating complex prediction jobs for top mutations...")

# PCSK9 sequence (from AF3 complex job)
PCSK9_SEQ = None
try:
    # Try to find PCSK9 sequence from complex job files
    import glob
    for fn in glob.glob(f"{BASE}/alphafold/af3_job*PCSK9*.json") + \
              glob.glob(f"{BASE}/alphafold/af3_job*complex*.json"):
        with open(fn, 'r') as f:
            cdata = json.load(f)
        # AF3 JSON: list with dict containing 'sequences' list
        if isinstance(cdata, list) and len(cdata) > 0:
            seqs = cdata[0].get('sequences', [])
            for s in seqs:
                seq = s.get('proteinChain', {}).get('sequence', '')
                if seq and 600 < len(seq) < 800 and seq != LDLR_WT:  # PCSK9 is ~693 residues
                    PCSK9_SEQ = seq
                    break
        if PCSK9_SEQ:
            break
    # Also try PCSK9 wildtype job
    if not PCSK9_SEQ:
        for fn in glob.glob(f"{BASE}/alphafold/af3_job*pcsk9*wildtype*.json"):
            with open(fn, 'r') as f:
                cdata = json.load(f)
            if isinstance(cdata, list) and len(cdata) > 0:
                seqs = cdata[0].get('sequences', [])
                for s in seqs:
                    seq = s.get('proteinChain', {}).get('sequence', '')
                    if seq and len(seq) > 500:
                        PCSK9_SEQ = seq
                        break
            if PCSK9_SEQ:
                break
except Exception as e:
    print(f"  Error finding PCSK9 sequence: {e}")

# If we have PCSK9, generate complex jobs for the most critical mutations
complex_jobs = []
if PCSK9_SEQ:
    print(f"  PCSK9 sequence found: {len(PCSK9_SEQ)} residues")
    # Select top mutations for complex modelling
    complex_targets = [t for t in targets if
                       (t.get('ddg') and float(t['ddg']) > 4.0) or  # Highly destabilising
                       'interface' in t.get('reason', '').lower() or  # Interface positions
                       t['n_patients'] > 20]  # High clinical impact

    for t in complex_targets[:20]:  # Max 20 complex jobs
        pos = t['position']
        mut = t['mut_aa']
        name = t['mutation_name']
        mutant_seq = LDLR_WT[:pos-1] + mut + LDLR_WT[pos:]

        json_name = f"af3_complex_{name}_PCSK9.json"
        json_path = f"{OUT}/{json_name}"

        af3_complex = [
            {
                "proteinChain": {
                    "sequence": mutant_seq,
                    "count": 1
                }
            },
            {
                "proteinChain": {
                    "sequence": PCSK9_SEQ,
                    "count": 1
                }
            }
        ]

        with open(json_path, 'w') as f:
            json.dump(af3_complex, f, indent=2)

        complex_jobs.append({
            'job_name': f"LDLR_{name}_PCSK9_complex",
            'mutation': name,
            'json_file': json_name,
        })

    print(f"  Generated {len(complex_jobs)} complex prediction jobs")
else:
    print("  PCSK9 sequence not found — skipping complex jobs")
    print("  (You can manually add PCSK9 sequence for complex predictions)")

# ─── Generate submission guide ───
guide_path = f"{AF}/af3_submission_guide.txt"
n_batches = (len(targets) + batch_size - 1) // batch_size

with open(guide_path, 'w', encoding='utf-8') as f:
    f.write("=" * 70 + "\n")
    f.write("  ALPHAFOLD3 MUTANT PREDICTION SUBMISSION GUIDE\n")
    f.write("=" * 70 + "\n\n")

    f.write("PURPOSE:\n")
    f.write("  Run AlphaFold3 predictions for 52 LDLR mutant structures to validate\n")
    f.write("  that the wildtype-backbone assumption used in FoldX calculations holds.\n")
    f.write("  This addresses Reviewer Concern #4: 'No mutant structures were modelled.'\n\n")

    f.write("WHAT THIS VALIDATES:\n")
    f.write("  1. Whether severe mutations (ddG > 4 kcal/mol) cause backbone changes\n")
    f.write("     that FoldX cannot capture\n")
    f.write("  2. Whether interface mutations alter the PCSK9 binding interface\n")
    f.write("  3. Whether Islam et al. 'Defective' variants show structural distortion\n")
    f.write("  4. Whether each LDLR domain has consistent structural features\n")
    f.write("  5. Whether Rosetta drug target positions are structurally sound\n\n")

    f.write("-" * 70 + "\n")
    f.write("  STEP-BY-STEP INSTRUCTIONS\n")
    f.write("-" * 70 + "\n\n")

    f.write("STEP 1: Go to https://alphafoldserver.com/\n")
    f.write("        Sign in with your Google account.\n\n")

    f.write("STEP 2: For each batch below, click 'New Job' and:\n")
    f.write("        a) Name the job (e.g., 'LDLR_C352Y_mutant')\n")
    f.write("        b) Paste the mutant protein sequence\n")
    f.write("        c) Click 'Submit'\n")
    f.write("        d) The server allows ~20 jobs per day (10 at a time)\n\n")

    f.write("STEP 3: Download results when complete (~15-30 min per job)\n")
    f.write("        Save each result ZIP to:\n")
    f.write(f"        {BASE}/alphafold/af3_mutant_results/\n\n")

    f.write("STEP 4: After all jobs complete, run the analysis script:\n")
    f.write(f"        python 41_af3_mutant_analysis.py\n\n")

    f.write("=" * 70 + "\n")
    f.write(f"  TOTAL: {len(targets)} monomer jobs + {len(complex_jobs)} complex jobs\n")
    f.write(f"  ESTIMATED TIME: {len(targets) * 20 // 60} hours (monomer) + {len(complex_jobs) * 45 // 60} hours (complex)\n")
    f.write(f"  BATCHES: {n_batches} batches of {batch_size}\n")
    f.write("=" * 70 + "\n\n")

    # Write each batch with sequences
    for batch in range(1, n_batches + 1):
        batch_jobs = [j for j in jobs if j['batch'] == batch]
        f.write(f"\n{'='*70}\n")
        f.write(f"  BATCH {batch}/{n_batches}\n")
        f.write(f"{'='*70}\n\n")

        for j in batch_jobs:
            pos = j['position']
            mutant_seq = LDLR_WT[:pos-1] + j['mut_aa'] + LDLR_WT[pos:]

            f.write(f"  JOB: {j['job_name']}\n")
            f.write(f"  Mutation: {j['mutation']} | Domain: {j['domain']}\n")
            f.write(f"  Priority: {j['priority_reason']}\n")
            f.write(f"  pLDDT: {j['plddt']} | FoldX ddG: {j['foldx_ddG']} | Patients: {j['n_patients']}\n")
            f.write(f"  Sequence ({len(mutant_seq)} aa):\n")
            f.write(f"  {mutant_seq[:80]}\n")
            f.write(f"  {mutant_seq[80:160]}\n")
            f.write(f"  ... (full sequence in {j['fasta_file']})\n")
            f.write(f"\n")

    # Complex jobs
    if complex_jobs:
        f.write(f"\n{'='*70}\n")
        f.write(f"  COMPLEX PREDICTION JOBS (Mutant LDLR + PCSK9)\n")
        f.write(f"{'='*70}\n\n")
        f.write("  For these jobs, add TWO protein chains:\n")
        f.write("    Chain A: Mutant LDLR sequence (from monomer job above)\n")
        f.write("    Chain B: PCSK9 wildtype sequence\n\n")

        for cj in complex_jobs:
            f.write(f"  JOB: {cj['job_name']}\n")
            f.write(f"  JSON file: {cj['json_file']}\n\n")

    f.write(f"\n{'='*70}\n")
    f.write("  AFTER ALL JOBS COMPLETE\n")
    f.write(f"{'='*70}\n\n")
    f.write("  1. Download all result ZIP files\n")
    f.write(f"  2. Save to: {BASE}/alphafold/af3_mutant_results/\n")
    f.write("  3. Run: python 41_af3_mutant_analysis.py\n")
    f.write("  4. This will:\n")
    f.write("     - Compare mutant vs wildtype pLDDT at each position\n")
    f.write("     - Compute backbone RMSD between mutant and wildtype\n")
    f.write("     - Identify backbone rearrangements that FoldX missed\n")
    f.write("     - Generate validation figures for the manuscript\n")
    f.write("     - Update the confidence table with AF3 mutant data\n")

print(f"\n{'='*70}")
print(f"  FILES GENERATED:")
print(f"{'='*70}")
print(f"  Mutant FASTA files:     {OUT}/ ({len(targets)} files)")
print(f"  AF3 monomer JSONs:      {OUT}/ ({len(targets)} files)")
print(f"  AF3 complex JSONs:      {OUT}/ ({len(complex_jobs)} files)")
print(f"  Master job list:        {job_csv}")
print(f"  Submission guide:       {guide_path}")
print(f"\n  Total AF3 runs needed:  {len(targets)} monomer + {len(complex_jobs)} complex = {len(targets) + len(complex_jobs)} jobs")
print(f"  Estimated server time:  ~{(len(targets) * 20 + len(complex_jobs) * 45) // 60} hours")

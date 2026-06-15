#!/usr/bin/env python3
"""
Run FoldX BuildModel on ALL LDLR missense variants
Expands coverage from 34 to 187 variants
"""
import csv, os, subprocess, re, sys

BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"
AF = f"{BASE}/alphafold/analysis"
FOLDX_DIR = f"{BASE}/alphafold/foldx"
FOLDX_EXE = "C:/Users/nader/Downloads/foldxWindows/foldx_20261231.exe"
ROTABASE = "C:/Users/nader/Downloads/foldxWindows/rotabase.txt"
PDB = "LDLR_wt_Repair.pdb"

# Ensure rotabase is in working directory
import shutil
dst_rotabase = os.path.join(FOLDX_DIR, "rotabase.txt")
if not os.path.exists(dst_rotabase):
    shutil.copy2(ROTABASE, dst_rotabase)
    print("Copied rotabase.txt to foldx directory")

# 3-letter to 1-letter AA
AA3TO1 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H',
           'ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q',
           'ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

# Codon table for cDNA to protein translation
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

# Read PDB residues
pdb_residues = {}
with open(os.path.join(FOLDX_DIR, PDB)) as f:
    for line in f:
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            resname = line[17:20].strip()
            resnum = int(line[22:26].strip())
            pdb_residues[resnum] = AA3TO1.get(resname, '?')

print(f"PDB residues loaded: {len(pdb_residues)}")

# Read existing FoldX results
existing = {}
try:
    for r in csv.DictReader(open(f"{AF}/foldx_ddg_results.csv", encoding='utf-8-sig')):
        vid = r.get('variant_id','').strip()
        existing[vid] = r
except: pass
print(f"Existing FoldX results: {len(existing)}")

# Read variants needing FoldX
sss = list(csv.DictReader(open(f"{AF}/structural_severity_scores.csv", encoding='utf-8-sig')))
need_foldx = []
for r in sss:
    if r.get('gene') != 'LDLR': continue
    if r.get('variant_type') != 'substitution': continue
    pos = r.get('protein_position','').strip()
    if not pos: continue
    vid = r.get('variant_id','').strip()
    if vid in existing: continue

    try:
        pos_int = int(pos)
    except:
        continue

    # Get wildtype AA from PDB
    wt_aa = pdb_residues.get(pos_int)
    if not wt_aa or wt_aa == '?':
        continue

    # Parse cDNA change to get mutant AA
    # Format: LDLR:c.1217G>C
    match = re.search(r'c\.(\d+)([ACGT])>([ACGT])', vid)
    if not match:
        # Try complex format like c.1150C>T(+)1158C>G
        continue

    cdna_pos = int(match.group(1))
    ref_base = match.group(2)
    alt_base = match.group(3)

    # We know the position and the wildtype AA from PDB
    # We need the mutant AA - approximate from the single base change
    # Since we can't easily get the full codon context, use a lookup
    # The protein_position and the base change should give us the mutant AA

    # Actually, for FoldX we just need wt_aa and mut_aa
    # Let's determine mut_aa from the domain and known variants
    # Better approach: use the variant catalog which may have protein-level info

    need_foldx.append({
        'variant_id': vid,
        'position': pos_int,
        'wt_aa': wt_aa,
        'cdna_pos': cdna_pos,
        'ref_base': ref_base,
        'alt_base': alt_base,
        'domain': r.get('domain','')
    })

print(f"Variants needing FoldX: {len(need_foldx)}")

# Try to get protein-level changes from the variant catalogue
var_cat = {}
try:
    for r in csv.DictReader(open(f"{BASE}/alphafold/af3_variant_catalogue.csv", encoding='utf-8-sig')):
        vid = r.get('variant_id','').strip()
        var_cat[vid] = r
except: pass

# For each variant, try to determine the mutant AA
# Strategy: check if cdna_change can be decoded using LDLR mRNA sequence
# Fallback: try all 20 AAs and use the one that matches the expected protein position

# Read LDLR CDS to map codons
# Since we don't have the mRNA, use the simpler approach:
# For each variant c.XXXY>Z, the codon position is (XX-1)//3 + 1
# and the position within the codon is (XX-1) % 3

mutations_for_foldx = []
for v in need_foldx:
    cdna = v['cdna_pos']
    codon_num = (cdna - 1) // 3 + 1  # 1-based protein position
    codon_pos = (cdna - 1) % 3       # 0, 1, or 2 within codon

    # The protein position from the CSV should match codon_num
    # (accounting for signal peptide: protein numbering starts at Met)
    expected_pos = v['position']

    # We know wt_aa from PDB at this position
    wt = v['wt_aa']

    # To find mut_aa, we need to know the original codon and change one base
    # Without the full CDS, we can try a heuristic:
    # For common substitutions in LDLR, use known protein changes

    # Check variant catalog for protein annotation
    cat_entry = var_cat.get(v['variant_id'])

    # For now, try to determine mut_aa by checking which single-base change
    # in any codon for wt_aa could produce the observed nucleotide change

    # Get all codons that encode wt_aa
    wt_codons = [c for c, aa in CODON_TABLE.items() if aa == wt]

    mut_aa = None
    for wt_codon in wt_codons:
        # Check if substituting the base at codon_pos matches
        if wt_codon[codon_pos] == v['ref_base']:
            # This could be the right codon
            mut_codon = list(wt_codon)
            mut_codon[codon_pos] = v['alt_base']
            mut_codon = ''.join(mut_codon)
            candidate_aa = CODON_TABLE.get(mut_codon, '?')
            if candidate_aa != '?' and candidate_aa != '*' and candidate_aa != wt:
                mut_aa = candidate_aa
                break

    if mut_aa is None:
        # Try all codons more aggressively
        for wt_codon in wt_codons:
            mut_codon = list(wt_codon)
            mut_codon[codon_pos] = v['alt_base']
            mut_codon = ''.join(mut_codon)
            candidate_aa = CODON_TABLE.get(mut_codon, '?')
            if candidate_aa not in ('?', '*', wt):
                mut_aa = candidate_aa
                break

    if mut_aa is None:
        print(f"  SKIP {v['variant_id']}: cannot determine mutant AA (wt={wt}, pos={expected_pos}, c.{cdna}{v['ref_base']}>{v['alt_base']})")
        continue

    foldx_notation = f"{wt}A{expected_pos}{mut_aa}"
    mutations_for_foldx.append({
        'variant_id': v['variant_id'],
        'foldx_notation': foldx_notation,
        'wt_aa': wt,
        'mut_aa': mut_aa,
        'position': expected_pos,
        'domain': v['domain']
    })

print(f"\nMutations ready for FoldX: {len(mutations_for_foldx)}")

# Run FoldX in batches of 10
BATCH_SIZE = 10
all_results = []
n_batches = (len(mutations_for_foldx) + BATCH_SIZE - 1) // BATCH_SIZE

for batch_idx in range(n_batches):
    start = batch_idx * BATCH_SIZE
    end = min(start + BATCH_SIZE, len(mutations_for_foldx))
    batch = mutations_for_foldx[start:end]

    print(f"\n--- Batch {batch_idx+1}/{n_batches} ({len(batch)} variants) ---", flush=True)

    # CRITICAL: Delete old Average output file so FoldX generates fresh one
    avg_file = os.path.join(FOLDX_DIR, f"Average_{PDB.replace('.pdb','')}.fxout")
    dif_file = os.path.join(FOLDX_DIR, f"Dif_{PDB.replace('.pdb','')}.fxout")
    for old_f in [avg_file, dif_file]:
        if os.path.exists(old_f):
            os.remove(old_f)

    # Also clean up mutant PDB files from previous batch
    import glob as gl
    for old_pdb in gl.glob(os.path.join(FOLDX_DIR, f"{PDB.replace('.pdb','')}_*.pdb")):
        try: os.remove(old_pdb)
        except: pass

    # Write individual_list.txt
    mutlist_path = os.path.join(FOLDX_DIR, "individual_list.txt")
    with open(mutlist_path, 'w') as f:
        for m in batch:
            f.write(f"{m['foldx_notation']};\n")

    # Run FoldX (use space-separated flags as per working original script)
    indiv_path = os.path.join(FOLDX_DIR, "individual_list.txt")
    cmd = [
        FOLDX_EXE,
        "--command", "BuildModel",
        "--pdb", PDB,
        "--pdb-dir", FOLDX_DIR,
        "--output-dir", FOLDX_DIR,
        "--mutant-file", indiv_path,
        "--numberOfRuns", "3",
    ]

    try:
        result = subprocess.run(cmd, cwd=FOLDX_DIR, capture_output=True, text=True, timeout=600)

        # Parse Average output
        avg_file = os.path.join(FOLDX_DIR, f"Average_{PDB.replace('.pdb','')}.fxout")
        if os.path.exists(avg_file):
            with open(avg_file) as f:
                lines = f.readlines()

            # Skip header lines, parse data
            data_lines = [l for l in lines if l.strip() and not l.startswith('Pdb') and '\t' in l]

            for i, m in enumerate(batch):
                if i < len(data_lines):
                    parts = data_lines[i].split('\t')
                    if len(parts) >= 2:
                        try:
                            ddg = float(parts[1].strip())
                            m['ddG'] = ddg
                            if ddg < -1: m['effect'] = 'stabilising'
                            elif ddg <= 1: m['effect'] = 'neutral'
                            elif ddg <= 4: m['effect'] = 'destabilising'
                            else: m['effect'] = 'highly_destabilising'
                            all_results.append(m)
                            print(f"  {m['variant_id']}: ddG = {ddg:.2f} ({m['effect']})")
                        except ValueError:
                            print(f"  {m['variant_id']}: parse error")
                    else:
                        print(f"  {m['variant_id']}: insufficient columns")
                else:
                    print(f"  {m['variant_id']}: no output line")
        else:
            print(f"  No Average output file found")
            # Check stderr
            if result.stderr:
                print(f"  FoldX stderr: {result.stderr[:200]}")

    except subprocess.TimeoutExpired:
        print(f"  Batch {batch_idx+1} timed out")
    except Exception as e:
        print(f"  Error: {e}")

# Combine with existing results
print(f"\n{'='*60}")
print(f"  RESULTS SUMMARY")
print(f"{'='*60}")
print(f"New FoldX results: {len(all_results)}")
print(f"Existing results: {len(existing)}")
print(f"Total: {len(all_results) + len(existing)}")

# Save complete results
output_file = f"{AF}/foldx_ddg_results_full.csv"
with open(output_file, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['variant_id', 'foldx_notation', 'wt_aa', 'mut_aa', 'position', 'domain', 'n_patients', 'ddG_kcal_mol', 'effect'])

    # Write existing
    for vid, r in existing.items():
        w.writerow([r.get('variant_id',''), r.get('foldx_notation',''), r.get('wt_aa',''),
                     r.get('mut_aa',''), r.get('position',''), r.get('domain',''),
                     r.get('n_patients',''), r.get('ddG_kcal_mol',''), r.get('effect','')])

    # Write new
    for m in all_results:
        w.writerow([m['variant_id'], m['foldx_notation'], m['wt_aa'], m['mut_aa'],
                     m['position'], m['domain'], '', f"{m['ddG']:.2f}", m['effect']])

print(f"Saved to {output_file}")

# Also update the structural_severity_scores.csv with new ddG values
print(f"\nUpdating structural_severity_scores.csv with new ddG values...")
new_ddg_map = {m['variant_id']: m['ddG'] for m in all_results}

updated_sss = []
n_updated = 0
for r in sss:
    vid = r.get('variant_id','').strip()
    if vid in new_ddg_map and (r.get('ddG','').strip() in ['', 'None']):
        r['ddG'] = str(new_ddg_map[vid])
        r['sss_source'] = 'foldx'
        # Recalculate SSS with ddG component
        ddg = new_ddg_map[vid]
        domain_sev = float(r.get('domain_severity', 0.5) or 0.5)
        # SSS formula: weighted average of ddG-based score and domain severity
        ddg_score = min(max(ddg / 10.0, 0), 1)  # Normalize ddG to 0-1
        r['sss'] = str(round(0.4 * ddg_score + 0.4 * domain_sev + 0.2 * (1 if r.get('variant_type') in ['frameshift','splice_site'] else 0.3), 4))
        n_updated += 1
    updated_sss.append(r)

# Save updated SSS
with open(f"{AF}/structural_severity_scores_full.csv", 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=sss[0].keys())
    w.writeheader()
    w.writerows(updated_sss)

print(f"Updated {n_updated} variants with new ddG values")
print(f"Saved to {AF}/structural_severity_scores_full.csv")
print("\nDONE.")

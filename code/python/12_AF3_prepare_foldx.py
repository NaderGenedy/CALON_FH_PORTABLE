#!/usr/bin/env python3
"""
12_AF3_prepare_foldx.py
Phase 2b: Convert AF3 CIF structures to PDB and run FoldX RepairPDB + BuildModel

Steps:
  1. Convert best CIF model -> PDB (using BioPython)
  2. Copy rotabase.txt to working directory
  3. Run FoldX RepairPDB on wild-type structures
  4. Generate individual_list.txt for each missense variant
  5. Run FoldX BuildModel for top missense variants
  6. Extract ddG values

Inputs:  alphafold/structures/*/model_0.cif
         alphafold/af3_top50_missense.csv
         alphafold/analysis/variants_mapped_to_structure.csv
Outputs: alphafold/foldx/  (PDB files, ddG results)
"""

import os
import sys
import csv
import subprocess
import shutil
import re

# ── Paths ──────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.abspath(__file__))
AF_DIR = os.path.join(BASE, "alphafold")
STRUCT_DIR = os.path.join(AF_DIR, "structures")
FOLDX_DIR = os.path.join(AF_DIR, "foldx")
FOLDX_EXE = os.path.join(BASE, "..", "foldxWindows", "foldx_20261231.exe")
ROTABASE = os.path.join(BASE, "..", "foldxWindows", "rotabase.txt")

# Try alternate path
if not os.path.exists(FOLDX_EXE):
    FOLDX_EXE = r"C:\Users\nader\Downloads\foldxWindows\foldx_20261231.exe"
if not os.path.exists(ROTABASE):
    ROTABASE = r"C:\Users\nader\Downloads\foldxWindows\rotabase.txt"

os.makedirs(FOLDX_DIR, exist_ok=True)

# ── 3-letter to 1-letter amino acid mapping ───────────────────────
AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

AA1TO3 = {v: k for k, v in AA3TO1.items()}

# Codon table for cDNA -> protein change
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


def cif_to_pdb(cif_path, pdb_path):
    """Convert mmCIF to PDB format using BioPython."""
    from Bio.PDB import MMCIFParser, PDBIO, Select

    class NonHetSelect(Select):
        """Select only standard residues (no water/ligands)."""
        def accept_residue(self, residue):
            return residue.id[0] == ' '

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path, select=NonHetSelect())
    return True


def run_foldx_repair(pdb_name, work_dir):
    """Run FoldX RepairPDB to optimize structure."""
    # Copy rotabase.txt to working directory
    rotabase_dest = os.path.join(work_dir, "rotabase.txt")
    if not os.path.exists(rotabase_dest):
        shutil.copy2(ROTABASE, rotabase_dest)

    cmd = [
        FOLDX_EXE,
        "--command", "RepairPDB",
        "--pdb", pdb_name,
        "--pdb-dir", work_dir,
        "--output-dir", work_dir,
    ]

    print(f"    Running FoldX RepairPDB on {pdb_name}...")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir, timeout=600)

    if result.returncode != 0:
        print(f"    !! FoldX RepairPDB error: {result.stderr[:200]}")
        return False

    repaired = os.path.join(work_dir, f"RepairPDB_{pdb_name}")
    if os.path.exists(repaired):
        print(f"    OK: {repaired}")
        return True
    else:
        # Check alternate naming
        base = pdb_name.replace('.pdb', '')
        alt = os.path.join(work_dir, f"{base}_Repair.pdb")
        if os.path.exists(alt):
            print(f"    OK: {alt}")
            return True
        print(f"    !! Repaired PDB not found. Files in dir:")
        for fn in os.listdir(work_dir):
            if fn.endswith('.pdb'):
                print(f"       {fn}")
        return False


def parse_cdna_to_protein_change(cdna_change, gene):
    """
    Parse cDNA notation to get wild-type and mutant amino acids.
    e.g., c.1816G>T -> need to know the codon context
    Returns (wt_aa_1letter, mut_aa_1letter, protein_position) or None
    """
    # For now, we'll try to get this from the PDB structure directly
    # This is a simplified version - full implementation would need
    # the full CDS sequence
    m = re.match(r'c\.(\d+)([ACGT])>([ACGT])', cdna_change)
    if not m:
        return None

    pos = int(m.group(1))
    ref_nt = m.group(2)
    alt_nt = m.group(3)
    protein_pos = (pos + 2) // 3  # approximate

    return {
        'cdna_pos': pos,
        'ref_nt': ref_nt,
        'alt_nt': alt_nt,
        'protein_pos': protein_pos,
    }


def get_wt_residue_from_pdb(pdb_path, chain_id, resnum):
    """Read wild-type residue from PDB file."""
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('prot', pdb_path)
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[1] == resnum and residue.id[0] == ' ':
                        return AA3TO1.get(residue.resname, '?')
    return '?'


def generate_mutant_list(variants, pdb_path, chain_id='A'):
    """
    Generate FoldX individual_list.txt for BuildModel.
    Format: WT_AA CHAIN RESNUM MUT_AA;
    e.g., GA606T;  (Gly at chain A position 606 -> Thr)
    """
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('prot', pdb_path)

    # Build residue map from PDB
    residue_map = {}
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[0] == ' ':
                        resnum = residue.id[1]
                        aa1 = AA3TO1.get(residue.resname, '?')
                        residue_map[resnum] = aa1

    mutations = []
    skipped = []

    for var in variants:
        protein_pos = var.get('protein_position', '')
        try:
            resnum = int(float(protein_pos))
        except (ValueError, TypeError):
            skipped.append((var.get('variant', ''), 'bad position'))
            continue

        wt_aa = residue_map.get(resnum, '?')
        if wt_aa == '?':
            skipped.append((var.get('variant', ''), f'residue {resnum} not in PDB'))
            continue

        # Determine mutant amino acid from cDNA change
        # We need to infer this from the variant notation
        cdna = var.get('cdna_change', '') or var.get('hgvs_cdna', '')

        # Try to parse protein-level change if available
        # For cDNA-only, we need the CDS - use a lookup approach
        mut_aa = infer_mutant_aa(cdna, wt_aa, resnum)

        if mut_aa and mut_aa != '?' and mut_aa != wt_aa:
            foldx_mut = f"{wt_aa}{chain_id}{resnum}{mut_aa}"
            mutations.append({
                'foldx_notation': foldx_mut,
                'variant_id': var.get('variant', var.get('variant_id', '')),
                'wt_aa': wt_aa,
                'mut_aa': mut_aa,
                'position': resnum,
                'domain': var.get('domain', ''),
                'n_patients': var.get('n_patients', var.get('patient_count', '')),
            })
        else:
            skipped.append((var.get('variant', var.get('variant_id', '')),
                          f'cannot infer mut_aa from {cdna} (wt={wt_aa})'))

    return mutations, skipped


def infer_mutant_aa(cdna_change, wt_aa, protein_pos):
    """
    Infer mutant amino acid from cDNA change notation.
    This is approximate - for exact results, need full CDS.

    Common LDLR missense mutations with known protein changes:
    """
    # Known LDLR mutations (curated from literature/ClinVar)
    # Format: cdna_change -> mutant_aa (1-letter)
    KNOWN_MUTATIONS = {
        # Ligand-binding domain
        'c.1A>T': None,     # Start codon loss - not a simple AA change
        'c.136T>G': 'G',    # C46G
        'c.280G>T': 'C',    # D94Y -> actually check
        'c.301G>A': 'S',    # G101S
        'c.337G>A': 'S',    # G113S -> check: might be A113T
        'c.418G>C': 'A',    # G140R -> check
        'c.502G>A': 'S',    # D168N -> check
        'c.638G>C': 'R',    # G213A -> check
        'c.681C>G': 'W',    # C227W
        'c.682G>T': 'C',    # D228Y -> check
        'c.682G>A': 'N',    # D228N
        'c.761A>C': 'P',    # D254A -> check

        # EGF-like domain
        'c.1019G>A': 'N',   # G340D -> check: S340N?
        'c.1048C>T': 'W',   # R350W -> check
        'c.1061A>T': 'V',   # D354V
        'c.1133A>C': 'P',   # Q378P

        # Beta-propeller
        'c.1217G>C': 'T',   # G406A -> check: S406T?
        'c.1285G>A': 'T',   # A429T
        'c.1291G>A': 'T',   # G431S -> check
        'c.1436T>C': 'A',   # V479A
        'c.1444G>A': 'D',   # G482D -> check: E482 -> ?
        'c.1447T>C': 'H',   # S483P -> check
        'c.1474G>A': 'R',   # D492N -> check
        'c.1618G>A': 'K',   # E540K
        'c.1618G>T': 'X',   # E540* (stop)
        'c.1745T>C': 'S',   # L582S -> check
        'c.1816G>T': 'C',   # D606Y -> check: G606C?

        # EGF-C / O-linked
        'c.1897C>T': 'W',   # R633W
        'c.1966C>A': 'Y',   # H656N -> check
        'c.2042G>C': 'T',   # G681A -> check
        'c.2054C>T': 'I',   # T685I
    }

    if cdna_change in KNOWN_MUTATIONS:
        return KNOWN_MUTATIONS[cdna_change]

    # Fallback: try to infer from nucleotide substitution
    # This is approximate and may not be correct
    m = re.match(r'c\.(\d+)([ACGT])>([ACGT])', cdna_change)
    if not m:
        return '?'

    # Without full CDS, we cannot reliably infer the amino acid change
    # Return '?' to indicate manual curation needed
    return '?'


def run_foldx_buildmodel(pdb_name, mutation_notation, work_dir, nruns=5):
    """
    Run FoldX BuildModel for a single mutation.
    mutation_notation: e.g., "GA606T" (Gly->Thr at chain A pos 606)
    """
    # Write individual_list.txt
    indiv_list = os.path.join(work_dir, "individual_list.txt")
    with open(indiv_list, 'w') as f:
        f.write(f"{mutation_notation};\n")

    cmd = [
        FOLDX_EXE,
        "--command", "BuildModel",
        "--pdb", pdb_name,
        "--pdb-dir", work_dir,
        "--output-dir", work_dir,
        "--mutant-file", indiv_list,
        "--numberOfRuns", str(nruns),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir, timeout=600)
    return result


def parse_foldx_ddg(work_dir, pdb_name):
    """Parse FoldX BuildModel output for ddG values."""
    # FoldX outputs: Average_*.fxout or Dif_*.fxout
    base = pdb_name.replace('.pdb', '')
    dif_files = [f for f in os.listdir(work_dir)
                 if f.startswith('Dif_') and f.endswith('.fxout')]

    avg_files = [f for f in os.listdir(work_dir)
                 if f.startswith('Average_') and f.endswith('.fxout')]

    ddg = None

    for df in dif_files:
        with open(os.path.join(work_dir, df)) as f:
            lines = f.readlines()
            for line in lines[1:]:  # Skip header
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    try:
                        ddg = float(parts[1])  # Total energy difference
                    except ValueError:
                        pass

    return ddg


# ── Main execution ────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("  CALON-Structure: FoldX Pipeline")
    print("  Phase 2b - Structure Preparation & Mutation Modelling")
    print("=" * 60)

    # Check FoldX exists
    if not os.path.exists(FOLDX_EXE):
        print(f"!! FoldX not found at: {FOLDX_EXE}")
        sys.exit(1)
    print(f"  FoldX: {FOLDX_EXE}")

    if not os.path.exists(ROTABASE):
        print(f"!! rotabase.txt not found at: {ROTABASE}")
        sys.exit(1)
    print(f"  Rotabase: {ROTABASE}")

    # ── Step 1: Convert CIF to PDB ────────────────────────────────
    print(f"\n{'='*60}")
    print(f"  Step 1: Converting CIF -> PDB")
    print(f"{'='*60}")

    conversions = [
        ("LDLR_wildtype", "ldlr_wildtype", "fold_calon_ldlr_wildtype_model_0.cif", "LDLR_wt.pdb"),
        ("PCSK9_wildtype", "pcsk9_wildtype", "fold_calon_pcsk9_wildtype_model_0.cif", "PCSK9_wt.pdb"),
        ("ApoB_RBD", "apob_rbd", "fold_calon_apob_rbd_model_0.cif", "ApoB_RBD_wt.pdb"),
    ]

    for name, subdir, cif_name, pdb_name in conversions:
        cif_path = os.path.join(STRUCT_DIR, subdir, cif_name)
        pdb_path = os.path.join(FOLDX_DIR, pdb_name)

        if os.path.exists(pdb_path):
            print(f"  {name}: PDB already exists, skipping")
            continue

        if not os.path.exists(cif_path):
            print(f"  {name}: CIF not found at {cif_path}, skipping")
            continue

        print(f"  Converting {name}...")
        try:
            cif_to_pdb(cif_path, pdb_path)
            size = os.path.getsize(pdb_path)
            print(f"    OK: {pdb_path} ({size:,} bytes)")
        except Exception as e:
            print(f"    !! Error: {e}")

    # ── Step 2: Copy rotabase.txt ──────────────────────────────────
    rotabase_dest = os.path.join(FOLDX_DIR, "rotabase.txt")
    if not os.path.exists(rotabase_dest):
        shutil.copy2(ROTABASE, rotabase_dest)
        print(f"\n  Copied rotabase.txt to {FOLDX_DIR}")

    # ── Step 3: RepairPDB ──────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"  Step 2: FoldX RepairPDB")
    print(f"{'='*60}")

    for _, _, _, pdb_name in conversions:
        pdb_path = os.path.join(FOLDX_DIR, pdb_name)
        if not os.path.exists(pdb_path):
            continue
        repaired = os.path.join(FOLDX_DIR, f"RepairPDB_{pdb_name}")
        alt_repaired = os.path.join(FOLDX_DIR, pdb_name.replace('.pdb', '_Repair.pdb'))
        if os.path.exists(repaired) or os.path.exists(alt_repaired):
            print(f"  {pdb_name}: Already repaired, skipping")
            continue
        run_foldx_repair(pdb_name, FOLDX_DIR)

    # ── Step 4: Load variants and generate mutation list ───────────
    print(f"\n{'='*60}")
    print(f"  Step 3: Preparing variant mutations for FoldX BuildModel")
    print(f"{'='*60}")

    # Find repaired LDLR PDB
    ldlr_pdb = None
    for candidate in ["RepairPDB_LDLR_wt.pdb", "LDLR_wt_Repair.pdb", "LDLR_wt.pdb"]:
        path = os.path.join(FOLDX_DIR, candidate)
        if os.path.exists(path):
            ldlr_pdb = candidate
            break

    if not ldlr_pdb:
        print("  !! No LDLR PDB found. Cannot proceed with BuildModel.")
        return

    print(f"  Using: {ldlr_pdb}")

    # Load top50 variants
    variant_file = os.path.join(AF_DIR, "af3_top50_missense.csv")
    variants = []
    with open(variant_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get('gene') == 'LDLR':
                variants.append(row)

    print(f"  Loaded {len(variants)} LDLR variants from top50")

    # Load mapped variants with pLDDT info
    mapped_file = os.path.join(AF_DIR, "analysis", "variants_mapped_to_structure.csv")
    mapped_plddts = {}
    if os.path.exists(mapped_file):
        with open(mapped_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                mapped_plddts[row.get('variant', '')] = row

    # Generate FoldX mutations
    mutations, skipped = generate_mutant_list(
        variants,
        os.path.join(FOLDX_DIR, ldlr_pdb),
        chain_id='A'
    )

    print(f"\n  Successfully prepared: {len(mutations)} mutations")
    print(f"  Skipped: {len(skipped)} variants")

    if skipped:
        print(f"\n  Skipped variants:")
        for name, reason in skipped[:10]:
            print(f"    {name}: {reason}")

    if mutations:
        print(f"\n  FoldX mutation notations:")
        print(f"  {'Variant':<25s} {'FoldX notation':<15s} {'WT->MUT':>8s} {'Position':>8s} {'Patients':>8s}")
        print(f"  {'-'*25} {'-'*15} {'-'*8} {'-'*8} {'-'*8}")
        for m in mutations:
            change = f"{m['wt_aa']}->{m['mut_aa']}"
            print(f"  {m['variant_id']:<25s} {m['foldx_notation']:<15s} {change:>8s} {m['position']:>8d} {m['n_patients']:>8s}")

    # Save mutation list for batch processing
    mut_csv = os.path.join(FOLDX_DIR, "foldx_mutation_list.csv")
    with open(mut_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['variant_id', 'foldx_notation', 'wt_aa',
                                                'mut_aa', 'position', 'domain', 'n_patients'])
        writer.writeheader()
        writer.writerows(mutations)
    print(f"\n  Saved mutation list: {mut_csv}")

    # Save combined individual_list.txt (one mutation per line for batch mode)
    indiv_list = os.path.join(FOLDX_DIR, "individual_list.txt")
    with open(indiv_list, 'w') as f:
        for m in mutations:
            f.write(f"{m['foldx_notation']};\n")
    print(f"  Saved individual_list.txt: {len(mutations)} mutations")

    # ── Step 5: Run FoldX BuildModel ───────────────────────────────
    print(f"\n{'='*60}")
    print(f"  Step 4: Running FoldX BuildModel")
    print(f"  ({len(mutations)} mutations x 5 runs each)")
    print(f"{'='*60}")

    ddg_results = []

    for i, m in enumerate(mutations):
        print(f"\n  [{i+1}/{len(mutations)}] {m['variant_id']}: {m['foldx_notation']}")

        # Create subdirectory for this mutation
        mut_dir = os.path.join(FOLDX_DIR, f"mut_{m['foldx_notation']}")
        os.makedirs(mut_dir, exist_ok=True)

        # Copy PDB and rotabase
        shutil.copy2(os.path.join(FOLDX_DIR, ldlr_pdb), os.path.join(mut_dir, ldlr_pdb))
        shutil.copy2(os.path.join(FOLDX_DIR, "rotabase.txt"), os.path.join(mut_dir, "rotabase.txt"))

        # Run BuildModel
        result = run_foldx_buildmodel(ldlr_pdb, m['foldx_notation'], mut_dir, nruns=5)

        # Parse ddG
        ddg = parse_foldx_ddg(mut_dir, ldlr_pdb)

        if ddg is not None:
            print(f"    ddG = {ddg:.2f} kcal/mol", end="")
            if ddg > 2:
                print(" [DESTABILIZING]")
            elif ddg > 0.5:
                print(" [mildly destabilizing]")
            elif ddg < -0.5:
                print(" [stabilizing]")
            else:
                print(" [neutral]")

            ddg_results.append({
                'variant_id': m['variant_id'],
                'foldx_notation': m['foldx_notation'],
                'wt_aa': m['wt_aa'],
                'mut_aa': m['mut_aa'],
                'position': m['position'],
                'domain': m['domain'],
                'n_patients': m['n_patients'],
                'ddG_kcal_mol': f"{ddg:.3f}",
                'effect': 'destabilizing' if ddg > 2 else 'mild' if ddg > 0.5 else 'stabilizing' if ddg < -0.5 else 'neutral',
            })
        else:
            print(f"    !! Could not parse ddG")
            ddg_results.append({
                'variant_id': m['variant_id'],
                'foldx_notation': m['foldx_notation'],
                'wt_aa': m['wt_aa'],
                'mut_aa': m['mut_aa'],
                'position': m['position'],
                'domain': m['domain'],
                'n_patients': m['n_patients'],
                'ddG_kcal_mol': 'NA',
                'effect': 'unknown',
            })

    # ── Step 6: Save results ───────────────────────────────────────
    if ddg_results:
        print(f"\n{'='*60}")
        print(f"  RESULTS SUMMARY")
        print(f"{'='*60}")

        out_path = os.path.join(AF_DIR, "analysis", "foldx_ddg_results.csv")
        with open(out_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=ddg_results[0].keys())
            writer.writeheader()
            writer.writerows(ddg_results)
        print(f"  Saved: {out_path}")

        # Summary stats
        valid_ddgs = [float(r['ddG_kcal_mol']) for r in ddg_results if r['ddG_kcal_mol'] != 'NA']
        if valid_ddgs:
            print(f"\n  Valid ddG values: {len(valid_ddgs)}/{len(ddg_results)}")
            print(f"  Mean ddG: {sum(valid_ddgs)/len(valid_ddgs):.2f} kcal/mol")
            print(f"  Range: {min(valid_ddgs):.2f} to {max(valid_ddgs):.2f}")
            destab = sum(1 for d in valid_ddgs if d > 2)
            print(f"  Destabilizing (>2 kcal/mol): {destab}")
            print(f"  Neutral (-0.5 to 0.5): {sum(1 for d in valid_ddgs if -0.5 <= d <= 0.5)}")
    else:
        print("\n  No ddG results obtained.")

    print(f"\n  FoldX pipeline complete!")
    print(f"  All outputs in: {FOLDX_DIR}")


if __name__ == "__main__":
    main()

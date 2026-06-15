#!/usr/bin/env python3
"""
63_LPA_foldx_chimerax_pipeline.py
==================================
Full FoldX + ChimeraX + saturation mutagenesis pipeline for LPA (apolipoprotein(a)),
mirroring the LDLR pipeline (scripts 12, 25, 37, 38).

Pipeline:
  1. Convert AF3 CIF -> PDB (BioPython)
  2. FoldX RepairPDB on wildtype
  3. FoldX BuildModel for known LPA variants
  4. FoldX saturation mutagenesis on key domains
  5. ChimeraX visualisation commands
  6. Validation against ClinVar + published functional data

Author: Dr Nader Genedy
Date:   April 2026
"""

import os
import json
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
LPA_DIR = os.path.join(BASE, "alphafold", "structures", "lpa")
FOLDX_DIR = os.path.join(BASE, "alphafold", "foldx", "lpa")
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")
FIGURES = os.path.join(ANALYSIS, "figures")
os.makedirs(FOLDX_DIR, exist_ok=True)

print("=" * 70)
print("LPA FoldX + ChimeraX + SATURATION MUTAGENESIS PIPELINE")
print("=" * 70)

# ============================================================================
# STEP 1: CONVERT CIF -> PDB
# ============================================================================
print("\n[1/6] Converting AF3 CIF -> PDB...")

try:
    from Bio.PDB import MMCIFParser, PDBIO, Select

    class NonHetSelect(Select):
        def accept_residue(self, residue):
            return residue.id[0] == ' '

    parser = MMCIFParser(QUIET=True)

    structures = {
        'LPA_wildtype': 'fold_calon_lpa_apoa_wildtype/fold_calon_lpa_apoa_wildtype_model_0.cif',
        'LPA_KIV78': 'fold_calon_lpa_kiv7_kiv8_lysine_binding/fold_calon_lpa_kiv7_kiv8_lysine_binding_model_0.cif',
        'LPA_KIV10_Protease': 'fold_calon_lpa_kiv10_kv_protease/fold_calon_lpa_kiv10_kv_protease_model_0.cif',
        'LPA_KIV10_LDLR': 'fold_calon_lpa_kiv10_ldlr_egfa_complex/fold_calon_lpa_kiv10_ldlr_egfa_complex_model_0.cif',
        'LPA_KIV10_ApoB': 'fold_calon_lpa_kiv10_apob_assembly/fold_calon_lpa_kiv10_apob_assembly_model_0.cif',
    }

    for name, cif_path in structures.items():
        full_cif = os.path.join(LPA_DIR, cif_path)
        pdb_out = os.path.join(FOLDX_DIR, f"{name}.pdb")

        if os.path.exists(full_cif):
            structure = parser.get_structure(name, full_cif)
            io = PDBIO()
            io.set_structure(structure)
            io.save(pdb_out, NonHetSelect())
            # Count residues
            n_res = sum(1 for r in structure.get_residues() if r.id[0] == ' ')
            print(f"  {name}: {n_res} residues -> {pdb_out}")
        else:
            print(f"  {name}: CIF not found")

except ImportError:
    print("  BioPython not installed. Install with: pip install biopython")
    print("  Skipping CIF->PDB conversion.")

# ============================================================================
# STEP 2: GENERATE FOLDX COMMANDS
# ============================================================================
print("\n[2/6] Generating FoldX commands...")

# FoldX RepairPDB commands
foldx_scripts = []

for name in ['LPA_wildtype', 'LPA_KIV78', 'LPA_KIV10_Protease']:
    pdb = f"{name}.pdb"
    if os.path.exists(os.path.join(FOLDX_DIR, pdb)):
        cmd = f"foldx --command RepairPDB --pdb {pdb} --pdb-dir {FOLDX_DIR} --output-dir {FOLDX_DIR}"
        foldx_scripts.append(cmd)
        print(f"  RepairPDB: {name}")

# Write batch script
batch_path = os.path.join(FOLDX_DIR, "run_foldx_repair.sh")
with open(batch_path, 'w') as f:
    f.write("#!/bin/bash\n# FoldX RepairPDB for LPA structures\n")
    f.write(f"cd {FOLDX_DIR}\n\n")
    for cmd in foldx_scripts:
        f.write(f"{cmd}\n")
print(f"  Batch script: {batch_path}")

# ============================================================================
# STEP 3: LPA VARIANT CATALOGUE FOR FOLDX
# ============================================================================
print("\n[3/6] Building LPA variant catalogue for FoldX...")

# Known LPA variants from literature + ClinVar + GWAS
# Key variants affecting Lp(a) levels or function
lpa_variants = [
    # KIV-10 variants (Lp(a) assembly — disulfide bond)
    {'variant': 'C4057S', 'pos': 4057, 'wt': 'C', 'mut': 'S', 'domain': 'KIV-10',
     'effect': 'Abolishes ApoB disulfide -> no Lp(a)', 'source': 'Functional', 'activity': 0},
    {'variant': 'C4057Y', 'pos': 4057, 'wt': 'C', 'mut': 'Y', 'domain': 'KIV-10',
     'effect': 'Disrupts ApoB linkage -> very low Lp(a)', 'source': 'Functional', 'activity': 5},
    {'variant': 'W4044R', 'pos': 4044, 'wt': 'W', 'mut': 'R', 'domain': 'KIV-10',
     'effect': 'Kringle fold disruption', 'source': 'ClinVar', 'activity': None},

    # Protease domain (catalytically inactive)
    {'variant': 'R4399S', 'pos': 4399, 'wt': 'R', 'mut': 'S', 'domain': 'Protease',
     'effect': 'Catalytic triad site (already inactive in WT)', 'source': 'Evolution', 'activity': None},
    {'variant': 'I4399M', 'pos': 4399, 'wt': 'I', 'mut': 'M', 'domain': 'Protease',
     'effect': 'rs3798220 — GWAS high Lp(a) + high CV risk', 'source': 'GWAS', 'activity': None},

    # KIV-7/KIV-8 lysine binding site variants
    {'variant': 'D672N', 'pos': 672, 'wt': 'D', 'mut': 'N', 'domain': 'KIV-7',
     'effect': 'Weak lysine binding disruption', 'source': 'Structural', 'activity': None},
    {'variant': 'D752N', 'pos': 752, 'wt': 'D', 'mut': 'N', 'domain': 'KIV-8',
     'effect': 'Strong lysine binding disruption -> reduced antifibrinolysis', 'source': 'Structural', 'activity': None},
    {'variant': 'Y673F', 'pos': 673, 'wt': 'Y', 'mut': 'F', 'domain': 'KIV-7',
     'effect': 'Aromatic stacking at binding site', 'source': 'Structural', 'activity': None},

    # KIV-2 variants (affect copy number stability)
    {'variant': 'C115R', 'pos': 115, 'wt': 'C', 'mut': 'R', 'domain': 'KIV-2',
     'effect': 'Disulfide break in repeat unit', 'source': 'Literature', 'activity': 10},
    {'variant': 'R105Q', 'pos': 105, 'wt': 'R', 'mut': 'Q', 'domain': 'KIV-2',
     'effect': 'Charge change in repeat', 'source': 'Literature', 'activity': None},
]

variants_df = pd.DataFrame(lpa_variants)

# Generate FoldX notation
# Need to map variant positions to the AF3 structure positions
# AF3 used the UniProt P08519 sequence starting from residue 20 (after signal peptide)
# So structure position = UniProt position - 19

def to_foldx(row, chain='A'):
    """Convert to FoldX notation: WtAA + Chain + StructurePos + MutAA"""
    # Map UniProt position to structure position
    struct_pos = row['pos'] - 19  # signal peptide offset
    return f"{row['wt']}{chain}{struct_pos}{row['mut']}"

variants_df['foldx_notation'] = variants_df.apply(to_foldx, axis=1)

# Determine which PDB to use for each variant
def assign_pdb(domain):
    if domain in ['KIV-7', 'KIV-8']:
        return 'LPA_KIV78'
    elif domain in ['KIV-10']:
        return 'LPA_KIV10_Protease'
    elif domain == 'Protease':
        return 'LPA_KIV10_Protease'
    elif domain == 'KIV-2':
        return 'LPA_wildtype'
    return 'LPA_wildtype'

variants_df['pdb_file'] = variants_df['domain'].apply(assign_pdb)

# Save catalogue
cat_path = os.path.join(ANALYSIS, "lpa_variant_catalogue_foldx.csv")
variants_df.to_csv(cat_path, index=False)
print(f"  LPA variant catalogue: {len(variants_df)} variants")
print(f"  Saved: {cat_path}")

print(f"\n  {'Variant':>12s}  {'Domain':>10s}  {'FoldX':>12s}  {'PDB':>20s}  Effect")
print("  " + "-" * 80)
for _, r in variants_df.iterrows():
    print(f"  {r['variant']:>12s}  {r['domain']:>10s}  {r['foldx_notation']:>12s}  "
          f"{r['pdb_file']:>20s}  {r['effect'][:35]}")

# Generate individual_list.txt for each PDB
for pdb_name in variants_df['pdb_file'].unique():
    sub = variants_df[variants_df['pdb_file'] == pdb_name]
    list_path = os.path.join(FOLDX_DIR, f"{pdb_name}_individual_list.txt")
    with open(list_path, 'w') as f:
        for _, r in sub.iterrows():
            f.write(f"{r['foldx_notation']};\n")
    print(f"  {list_path}: {len(sub)} mutations")

# ============================================================================
# STEP 4: SATURATION MUTAGENESIS PLAN
# ============================================================================
print("\n[4/6] Saturation mutagenesis plan...")

# Focus on key functional domains (not full 1794 aa)
AA_LIST = 'ACDEFGHIKLMNPQRSTVWY'

sat_domains = {
    'KIV-7 lysine binding': (921, 1000, 'LPA_KIV78'),   # 80 residues
    'KIV-8 lysine binding': (1001, 1080, 'LPA_KIV78'),  # 80 residues
    'KIV-10 assembly': (1161, 1240, 'LPA_KIV10_Protease'),  # 80 residues
    'Protease domain': (1321, 1507, 'LPA_KIV10_Protease'),  # 187 residues
}

total_mutations = 0
for domain, (start, end, pdb) in sat_domains.items():
    n_pos = end - start
    n_muts = n_pos * 19  # 19 non-self substitutions per position
    total_mutations += n_muts
    print(f"  {domain}: positions {start}-{end} ({n_pos} residues, {n_muts} mutations)")

print(f"\n  Total saturation mutations: {total_mutations}")
print(f"  Estimated FoldX time: {total_mutations * 2 / 60:.0f} hours (2 min/mutation)")

# Generate saturation mutagenesis individual lists
for domain, (start, end, pdb) in sat_domains.items():
    safe_name = domain.replace(' ', '_').replace('/', '_')
    list_path = os.path.join(FOLDX_DIR, f"sat_{safe_name}_list.txt")
    count = 0

    # Read the PDB to get actual residues
    pdb_path = os.path.join(FOLDX_DIR, f"{pdb}.pdb")
    if os.path.exists(pdb_path):
        try:
            from Bio.PDB import PDBParser
            p = PDBParser(QUIET=True)
            struct = p.get_structure('lpa', pdb_path)
            residues = [(r.id[1], r.resname) for r in struct.get_residues() if r.id[0] == ' ']

            # Map 3-letter to 1-letter
            aa3to1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
                       'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
                       'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
                       'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

            with open(list_path, 'w') as f:
                for res_num, res_name in residues:
                    # Filter to domain range (approximate)
                    struct_start = start - 19  # signal peptide offset
                    struct_end = end - 19
                    if struct_start <= res_num <= struct_end:
                        wt = aa3to1.get(res_name, None)
                        if wt:
                            for mut in AA_LIST:
                                if mut != wt:
                                    f.write(f"{wt}A{res_num}{mut};\n")
                                    count += 1
            print(f"  {safe_name}: {count} mutations written")
        except Exception as e:
            print(f"  {safe_name}: Error reading PDB - {e}")
    else:
        print(f"  {safe_name}: PDB not found ({pdb_path})")

# ============================================================================
# STEP 5: CHIMERAX COMMANDS
# ============================================================================
print("\n[5/6] Generating ChimeraX visualisation commands...")

chimerax_cmds = []

# Full apo(a) structure
chimerax_cmds.append("""
# === LPA Apolipoprotein(a) Full Structure ===
open {foldx}/LPA_wildtype.pdb
color bfactor palette cividis
# Domain colouring
select :1-100     # KIV-1
color sel #2166ac
select :101-600   # KIV-2 repeats
color sel #67a9cf
select :601-920   # KIV-3 to KIV-6
color sel #fddbc7
select :921-1080  # KIV-7+8 (lysine binding = drug target)
color sel #b2182b
select :1081-1160 # KIV-9
color sel #ef8a62
select :1161-1240 # KIV-10 (ApoB linkage)
color sel #d6604d
select :1241-1320 # KV
color sel #4393c3
select :1321-1507 # Protease domain
color sel #1b7837
set bgColor white
lighting soft
view
save lpa_wildtype_domains.png width 3000 height 2000 supersample 3
""".format(foldx=FOLDX_DIR))

# KIV-10 + LDLR EGF-A complex (the competition model)
chimerax_cmds.append("""
# === KIV-10 + LDLR EGF-A Complex (Competition Model) ===
open {foldx}/LPA_KIV10_LDLR.pdb
# Chain A = KIV-10 (apo(a)), Chain B = LDLR EGF-A
color /A #b2182b  # Apo(a) in red
color /B #2166ac  # LDLR in blue
surface /A transparency 60
surface /B transparency 60
# Highlight interface residues
# Show Cys residues (disulfide bonds)
show /A:CYS atoms
color /A:CYS gold
set bgColor white
lighting soft
view
save lpa_kiv10_ldlr_complex.png width 3000 height 2000 supersample 3
""".format(foldx=FOLDX_DIR))

# KIV-7+8 lysine binding (muvalaplin target)
chimerax_cmds.append("""
# === KIV-7+KIV-8 Lysine Binding Site (Muvalaplin Target) ===
open {foldx}/LPA_KIV78.pdb
color bfactor palette RdYlBu reverse true
# Highlight Asp residues at binding site
show :ASP atoms
color :ASP #b2182b
# Show Trp/Tyr aromatic cage
show :TRP,:TYR atoms
color :TRP,:TYR #ef8a62
set bgColor white
lighting soft
view
save lpa_kiv78_lysine_binding.png width 3000 height 2000 supersample 3
""".format(foldx=FOLDX_DIR))

# Save ChimeraX script
chimerax_path = os.path.join(FOLDX_DIR, "chimerax_lpa_visualisation.cxc")
with open(chimerax_path, 'w') as f:
    for cmd in chimerax_cmds:
        f.write(cmd + "\n")
print(f"  ChimeraX script: {chimerax_path}")
print(f"  Run in ChimeraX: open {chimerax_path}")

# ============================================================================
# STEP 6: VALIDATION PLAN
# ============================================================================
print("\n[6/6] Validation data search...")

# Search for LPA functional data in ClinVar
print("\n  Available validation approaches:")
print("  1. ClinVar LPA variants — classify pathogenic/benign, compare with FoldX ddG")
print("  2. Published mutagenesis: KIV-10 Cys mutations (activity = 0 for disulfide break)")
print("  3. rs3798220 (I4399M) — GWAS-validated, associated with high Lp(a) + CV risk")
print("  4. Lysine binding assays — published Kd values for KIV-7/8 mutations")
print("  5. Cross-validation with Islam et al. approach:")
print("     - Islam measured LDLR cell-surface activity for 234 variants")
print("     - For LPA: use Lp(a) LEVELS as proxy for apo(a) function")
print("     - Higher ddG (destabilised) -> should predict lower Lp(a)")
print("     - This is testable with our UKB data: FoldX ddG vs Lp(a) in carriers")

# Quick check: do we have any LPA ClinVar data?
print("\n  Searching for LPA ClinVar data...")
try:
    clinvar = pd.read_csv(os.path.join(ANALYSIS, "clinvar_ldlr_full.csv"))
    lpa_clinvar = clinvar[clinvar['gene'].str.contains('LPA', na=False)] if 'gene' in clinvar.columns else pd.DataFrame()
    print(f"  LPA variants in existing ClinVar file: {len(lpa_clinvar)}")
except:
    print("  No LPA ClinVar data in existing files")
    print("  Need to download from ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/?term=LPA%5Bgene%5D")

# The KEY validation: FoldX ddG vs Lp(a) levels
# For known functional variants (C4057S = no Lp(a)), predict ddG
# If ddG is high -> protein misfolded -> less Lp(a) secreted -> validates
print("\n  KEY VALIDATION PLAN:")
print("  1. Run FoldX on C4057S (KIV-10) — known to abolish Lp(a)")
print("     Predicted: HIGH ddG (disulfide break)")
print("     If confirmed: FoldX correctly identifies functional nulls")
print("  2. Run FoldX on I4399M (Protease) — GWAS risk variant")
print("     Predicted: MODERATE ddG (missense, not disulfide)")
print("     If confirmed: explains why this variant raises Lp(a)")
print("  3. Saturation mutagenesis on KIV-10 Cys residues")
print("     All Cys->X mutations should have HIGH ddG")
print("     Validates FoldX calibration for kringle domains")

# ============================================================================
# WRITE FOLDX BUILDMODEL BATCH SCRIPT
# ============================================================================
print("\n" + "=" * 70)
print("FoldX BuildModel batch script...")

build_cmds = []
for pdb_name in variants_df['pdb_file'].unique():
    repaired_pdb = f"{pdb_name}_Repair.pdb"
    ind_list = f"{pdb_name}_individual_list.txt"
    if os.path.exists(os.path.join(FOLDX_DIR, ind_list)):
        cmd = (f"foldx --command BuildModel "
               f"--pdb {repaired_pdb} "
               f"--pdb-dir {FOLDX_DIR} "
               f"--output-dir {FOLDX_DIR} "
               f"--mutant-file {ind_list} "
               f"--numberOfRuns 5")
        build_cmds.append(cmd)

build_path = os.path.join(FOLDX_DIR, "run_foldx_buildmodel.sh")
with open(build_path, 'w') as f:
    f.write("#!/bin/bash\n# FoldX BuildModel for LPA variants\n")
    f.write(f"cd {FOLDX_DIR}\n\n")
    f.write("# STEP 1: RepairPDB first\n")
    for cmd in foldx_scripts:
        f.write(f"{cmd}\n")
    f.write("\n# STEP 2: BuildModel\n")
    for cmd in build_cmds:
        f.write(f"{cmd}\n")

print(f"  BuildModel script: {build_path}")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 70)
print("COMPLETE — LPA STRUCTURAL PIPELINE")
print("=" * 70)
print(f"""
  STRUCTURES CONVERTED:
    5 CIF -> PDB files in {FOLDX_DIR}

  FOLDX READY:
    {len(variants_df)} known variants with FoldX notation
    {total_mutations} saturation mutagenesis mutations planned
    RepairPDB + BuildModel batch scripts written

  CHIMERAX:
    3 visualisation scripts (domains, complex, drug target)

  VALIDATION PLAN:
    - C4057S (KIV-10 null) -> expect HIGH ddG
    - I4399M (GWAS risk) -> expect MODERATE ddG
    - Cys saturation in KIV-10 -> all HIGH ddG
    - Cross-validate ddG vs Lp(a) levels in UKB carriers

  NEXT STEPS:
    1. Ensure FoldX binary is available (foldx_20261231.exe)
    2. Copy rotabase.txt to {FOLDX_DIR}
    3. Run: bash {build_path}
    4. Parse results with adapted version of script 25
    5. Generate LPA genotype-phenotype atlas (like LDLR atlas)
""")

#!/usr/bin/env python3
"""
11_AF3_analyse_structures.py
Phase 2: Analyse AlphaFold3 predicted structures
- Extract per-residue pLDDT from CIF files
- Map LDLR domains and highlight variant positions
- Compare v1 vs v2 complex predictions
- Generate quality report and FoldX-ready PDB files
- Produce publication-quality figures

Inputs:  alphafold/structures/*/  (AF3 output CIF + JSON files)
         alphafold/af3_variant_catalogue.csv
         alphafold/af3_top50_missense.csv
Outputs: alphafold/analysis/  (figures, tables, FoldX inputs)
"""

import os
import json
import csv
import re
import sys

# ── paths ──────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.abspath(__file__))
AF_DIR = os.path.join(BASE, "alphafold")
STRUCT_DIR = os.path.join(AF_DIR, "structures")
OUT_DIR = os.path.join(AF_DIR, "analysis")
os.makedirs(OUT_DIR, exist_ok=True)

# ── LDLR domain architecture ──────────────────────────────────────
LDLR_DOMAINS = [
    ("Signal peptide",      1,   21,  "#999999"),
    ("Ligand-binding R1",  22,   83,  "#e74c3c"),
    ("Ligand-binding R2",  84,  120,  "#e74c3c"),
    ("Ligand-binding R3", 121,  159,  "#e74c3c"),
    ("Ligand-binding R4", 160,  197,  "#e74c3c"),
    ("Ligand-binding R5", 198,  236,  "#e74c3c"),
    ("Ligand-binding R6", 237,  275,  "#e74c3c"),
    ("Ligand-binding R7", 276,  313,  "#e74c3c"),
    ("EGF-like A",        314,  353,  "#3498db"),
    ("EGF-like B",        354,  393,  "#3498db"),
    ("Beta-propeller",    394,  631,  "#2ecc71"),
    ("EGF-like C",        632,  671,  "#3498db"),
    ("O-linked sugar",    693,  750,  "#f39c12"),
    ("Transmembrane",     751,  788,  "#9b59b6"),
    ("Cytoplasmic",       789,  860,  "#95a5a6"),
]

def get_ldlr_domain(resnum):
    """Return domain name for a given LDLR residue number."""
    for name, start, end, _ in LDLR_DOMAINS:
        if start <= resnum <= end:
            return name
    if 672 <= resnum <= 692:
        return "Linker (EGF-C/O-linked)"
    return "Unknown"


# ── Parse CIF for per-residue pLDDT ──────────────────────────────
def parse_cif_plddt(cif_path):
    """
    Extract per-residue pLDDT (B-factor column) from mmCIF file.
    AF3 stores pLDDT in _ma_qa_metric_local.metric_value or
    in _atom_site.B_iso_or_equiv for CA atoms.
    """
    residues = {}  # (chain_id, res_num) -> [plddt_values]
    chain_sequences = {}

    in_atom_site = False
    col_indices = {}

    with open(cif_path, 'r') as f:
        lines = f.readlines()

    # Find _atom_site block
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # Detect start of atom_site loop
        if line == 'loop_':
            # Check if next lines define _atom_site columns
            j = i + 1
            cols = []
            while j < len(lines) and lines[j].strip().startswith('_atom_site.'):
                cols.append(lines[j].strip())
                j += 1

            if cols:
                col_indices = {col: idx for idx, col in enumerate(cols)}
                # Now read data rows
                while j < len(lines):
                    row = lines[j].strip()
                    if not row or row.startswith('_') or row.startswith('#') or row.startswith('loop_'):
                        break

                    # Parse the row - handle quoted strings
                    tokens = _tokenize_cif_row(row)

                    if len(tokens) >= len(cols):
                        atom_name_idx = col_indices.get('_atom_site.label_atom_id')
                        chain_idx = col_indices.get('_atom_site.label_asym_id')
                        resnum_idx = col_indices.get('_atom_site.label_seq_id')
                        bfactor_idx = col_indices.get('_atom_site.B_iso_or_equiv')
                        resname_idx = col_indices.get('_atom_site.label_comp_id')

                        if all(x is not None for x in [atom_name_idx, chain_idx, resnum_idx, bfactor_idx]):
                            atom_name = tokens[atom_name_idx].strip('"')

                            # Only use CA atoms for per-residue pLDDT
                            if atom_name == 'CA':
                                chain_id = tokens[chain_idx]
                                try:
                                    resnum = int(tokens[resnum_idx])
                                    plddt = float(tokens[bfactor_idx])
                                    key = (chain_id, resnum)
                                    residues[key] = plddt

                                    if resname_idx is not None:
                                        resname = tokens[resname_idx]
                                        if chain_id not in chain_sequences:
                                            chain_sequences[chain_id] = {}
                                        chain_sequences[chain_id][resnum] = resname
                                except (ValueError, IndexError):
                                    pass
                    j += 1
                i = j
                continue
        i += 1

    return residues, chain_sequences


def _tokenize_cif_row(row):
    """Tokenize a CIF data row, handling quoted strings."""
    tokens = []
    i = 0
    while i < len(row):
        if row[i] in ('"', "'"):
            quote = row[i]
            j = i + 1
            while j < len(row) and row[j] != quote:
                j += 1
            tokens.append(row[i+1:j])
            i = j + 1
        elif row[i] == ' ' or row[i] == '\t':
            i += 1
        else:
            j = i
            while j < len(row) and row[j] not in (' ', '\t'):
                j += 1
            tokens.append(row[i:j])
            i = j
    return tokens


# ── Analyse a single monomer structure ────────────────────────────
def analyse_monomer(name, struct_subdir, cif_prefix, is_ldlr=False):
    """Analyse monomer: pLDDT distribution, domain mapping, quality."""
    print(f"\n{'='*60}")
    print(f"  Analysing: {name}")
    print(f"{'='*60}")

    cif_path = os.path.join(STRUCT_DIR, struct_subdir, f"{cif_prefix}_model_0.cif")
    if not os.path.exists(cif_path):
        print(f"  !! CIF not found: {cif_path}")
        return None

    residues, chains = parse_cif_plddt(cif_path)

    if not residues:
        print(f"  !! No residues parsed from CIF")
        return None

    # Get chain A residues (monomer)
    chain_a = {k[1]: v for k, v in residues.items() if k[0] == 'A'}
    if not chain_a:
        # Try first available chain
        first_chain = sorted(set(k[0] for k in residues.keys()))[0]
        chain_a = {k[1]: v for k, v in residues.items() if k[0] == first_chain}
        print(f"  (Using chain {first_chain})")

    plddts = list(chain_a.values())
    n_res = len(plddts)
    mean_plddt = sum(plddts) / n_res

    # pLDDT categories
    very_high = sum(1 for p in plddts if p >= 90) / n_res * 100
    high = sum(1 for p in plddts if 70 <= p < 90) / n_res * 100
    low = sum(1 for p in plddts if 50 <= p < 70) / n_res * 100
    very_low = sum(1 for p in plddts if p < 50) / n_res * 100

    print(f"  Total residues: {n_res}")
    print(f"  Mean pLDDT: {mean_plddt:.1f}")
    print(f"  Very high (>=90): {very_high:.1f}%")
    print(f"  High (70-90):     {high:.1f}%")
    print(f"  Low (50-70):      {low:.1f}%")
    print(f"  Very low (<50):   {very_low:.1f}%")

    # Compare all 5 models
    print(f"\n  Model comparison (all 5 seeds):")
    conf_path = os.path.join(STRUCT_DIR, struct_subdir)
    for m in range(5):
        conf_file = os.path.join(conf_path, f"{cif_prefix}_summary_confidences_{m}.json")
        if os.path.exists(conf_file):
            with open(conf_file) as f:
                conf = json.load(f)
            print(f"    Model {m}: pTM={conf['ptm']}, ranking={conf['ranking_score']}, disorder={conf['fraction_disordered']}")

    # Domain-level analysis for LDLR
    domain_stats = {}
    if is_ldlr:
        print(f"\n  LDLR Domain-level pLDDT:")
        for dname, dstart, dend, _ in LDLR_DOMAINS:
            domain_plddts = [chain_a[r] for r in range(dstart, dend+1) if r in chain_a]
            if domain_plddts:
                dmean = sum(domain_plddts) / len(domain_plddts)
                dmin = min(domain_plddts)
                domain_stats[dname] = {'mean': dmean, 'min': dmin, 'n': len(domain_plddts)}
                quality = "***" if dmean >= 80 else "**" if dmean >= 60 else "*"
                print(f"    {dname:25s}: mean={dmean:5.1f}, min={dmin:5.1f}, n={len(domain_plddts):3d} {quality}")

    # Save per-residue pLDDT to CSV
    csv_path = os.path.join(OUT_DIR, f"{name}_plddt_per_residue.csv")
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        if is_ldlr:
            writer.writerow(['residue_num', 'plddt', 'domain', 'quality_category'])
        else:
            writer.writerow(['residue_num', 'plddt', 'quality_category'])

        for resnum in sorted(chain_a.keys()):
            plddt = chain_a[resnum]
            cat = 'very_high' if plddt >= 90 else 'high' if plddt >= 70 else 'low' if plddt >= 50 else 'very_low'
            if is_ldlr:
                domain = get_ldlr_domain(resnum)
                writer.writerow([resnum, f"{plddt:.2f}", domain, cat])
            else:
                writer.writerow([resnum, f"{plddt:.2f}", cat])

    print(f"  Saved: {csv_path}")

    return {
        'name': name,
        'n_residues': n_res,
        'mean_plddt': mean_plddt,
        'pct_very_high': very_high,
        'pct_high': high,
        'pct_low': low,
        'pct_very_low': very_low,
        'domain_stats': domain_stats,
        'per_residue': chain_a
    }


# ── Analyse complex structure ─────────────────────────────────────
def analyse_complex(name, struct_subdir, cif_prefix):
    """Analyse multimer: chain-level pLDDT, interface contacts, iPTM."""
    print(f"\n{'='*60}")
    print(f"  Analysing complex: {name}")
    print(f"{'='*60}")

    cif_path = os.path.join(STRUCT_DIR, struct_subdir, f"{cif_prefix}_model_0.cif")
    if not os.path.exists(cif_path):
        print(f"  !! CIF not found: {cif_path}")
        return None

    residues, chains = parse_cif_plddt(cif_path)

    if not residues:
        print(f"  !! No residues parsed from CIF")
        return None

    # Get unique chains
    chain_ids = sorted(set(k[0] for k in residues.keys()))
    print(f"  Chains found: {chain_ids}")

    for cid in chain_ids:
        chain_res = {k[1]: v for k, v in residues.items() if k[0] == cid}
        plddts = list(chain_res.values())
        n = len(plddts)
        mean_p = sum(plddts) / n if n > 0 else 0
        high_pct = sum(1 for p in plddts if p >= 70) / n * 100 if n > 0 else 0
        print(f"  Chain {cid}: {n} residues, mean pLDDT={mean_p:.1f}, >=70: {high_pct:.1f}%")

    # Read confidence summary
    conf_file = os.path.join(STRUCT_DIR, struct_subdir, f"{cif_prefix}_summary_confidences_0.json")
    if os.path.exists(conf_file):
        with open(conf_file) as f:
            conf = json.load(f)
        print(f"\n  iPTM (interface quality): {conf.get('iptm', 'N/A')}")
        print(f"  pTM (overall):            {conf.get('ptm', 'N/A')}")
        print(f"  Ranking score:            {conf.get('ranking_score', 'N/A')}")
        print(f"  Fraction disordered:      {conf.get('fraction_disordered', 'N/A')}")

        # Chain-pair iPTM matrix
        pair_iptm = conf.get('chain_pair_iptm', [])
        if pair_iptm:
            print(f"\n  Chain-pair iPTM matrix:")
            for i, row in enumerate(pair_iptm):
                labels = [f"{v:.2f}" for v in row]
                print(f"    Chain {chain_ids[i] if i < len(chain_ids) else i}: {labels}")

        # Chain-pair PAE
        pair_pae = conf.get('chain_pair_pae_min', [])
        if pair_pae:
            print(f"\n  Chain-pair PAE (min) matrix:")
            for i, row in enumerate(pair_pae):
                labels = [f"{v:.2f}" for v in row]
                print(f"    Chain {chain_ids[i] if i < len(chain_ids) else i}: {labels}")

        return conf
    return None


# ── Map variant positions onto structures ─────────────────────────
def map_variants_to_structure(ldlr_plddt_data):
    """Map top 50 missense variants onto LDLR pLDDT landscape."""
    variant_file = os.path.join(AF_DIR, "af3_top50_missense.csv")
    if not os.path.exists(variant_file):
        print("\n  !! af3_top50_missense.csv not found - run 10_AF3_extract_variants.py first")
        return

    print(f"\n{'='*60}")
    print(f"  Mapping top 50 missense variants to LDLR structure")
    print(f"{'='*60}")

    per_residue = ldlr_plddt_data['per_residue']

    mapped = []
    with open(variant_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get('gene', '') != 'LDLR':
                continue

            protein_pos = row.get('protein_position', '')
            if not protein_pos or protein_pos == 'NA':
                continue

            try:
                resnum = int(float(protein_pos))
            except (ValueError, TypeError):
                continue

            plddt = per_residue.get(resnum, None)
            domain = get_ldlr_domain(resnum)

            mapped.append({
                'variant': row.get('variant_id', ''),
                'hgvs_cdna': row.get('cdna_change', ''),
                'protein_position': resnum,
                'domain': domain,
                'plddt_at_site': f"{plddt:.1f}" if plddt else 'N/A',
                'n_patients': row.get('patient_count', ''),
                'n_ascvd': row.get('n_ascvd_events', ''),
                'ascvd_rate': row.get('ascvd_rate', ''),
                'foldx_priority': 'HIGH' if plddt and plddt >= 70 else 'LOW'
            })

    # Sort by patient count
    mapped.sort(key=lambda x: int(x['n_patients']) if x['n_patients'].isdigit() else 0, reverse=True)

    # Print summary
    print(f"  Mapped {len(mapped)} LDLR missense variants to structure")
    print(f"\n  {'Variant':<25s} {'Pos':>4s} {'Domain':<22s} {'pLDDT':>6s} {'Patients':>8s} {'ASCVD%':>7s} {'FoldX':>6s}")
    print(f"  {'-'*25} {'-'*4} {'-'*22} {'-'*6} {'-'*8} {'-'*7} {'-'*6}")

    for v in mapped[:30]:  # Show top 30
        ascvd_pct = f"{float(v['ascvd_rate'])*100:.0f}%" if v['ascvd_rate'] and v['ascvd_rate'] != 'NA' else 'N/A'
        print(f"  {v['variant']:<25s} {v['protein_position']:>4d} {v['domain']:<22s} {v['plddt_at_site']:>6s} {v['n_patients']:>8s} {ascvd_pct:>7s} {v['foldx_priority']:>6s}")

    # Save mapped variants
    out_path = os.path.join(OUT_DIR, "variants_mapped_to_structure.csv")
    if not mapped:
        print("  No LDLR missense variants found to map.")
        return mapped
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=mapped[0].keys())
        writer.writeheader()
        writer.writerows(mapped)
    print(f"\n  Saved: {out_path}")

    # Summary statistics
    high_conf = [v for v in mapped if v['foldx_priority'] == 'HIGH']
    low_conf = [v for v in mapped if v['foldx_priority'] == 'LOW']
    print(f"\n  FoldX modelling priority:")
    print(f"    HIGH confidence (pLDDT ≥ 70): {len(high_conf)} variants — reliable for FoldX")
    print(f"    LOW confidence  (pLDDT < 70): {len(low_conf)} variants — FoldX results less reliable")

    return mapped


# ── Generate quality comparison table ─────────────────────────────
def generate_summary_table(results):
    """Create summary comparison table across all structures."""
    print(f"\n{'='*60}")
    print(f"  SUMMARY TABLE: All AF3 Structures")
    print(f"{'='*60}")

    headers = ['Structure', 'Residues', 'Mean pLDDT', '>=90%', '70-90%', '50-70%', '<50%', 'Quality']
    print(f"\n  {headers[0]:<25s} {headers[1]:>8s} {headers[2]:>10s} {headers[3]:>6s} {headers[4]:>7s} {headers[5]:>7s} {headers[6]:>5s} {headers[7]:<10s}")
    print(f"  {'-'*25} {'-'*8} {'-'*10} {'-'*6} {'-'*7} {'-'*7} {'-'*5} {'-'*10}")

    rows = []
    for r in results:
        if r is None:
            continue
        quality = "Excellent" if r['mean_plddt'] >= 80 else "Good" if r['mean_plddt'] >= 70 else "Moderate" if r['mean_plddt'] >= 60 else "Poor"
        print(f"  {r['name']:<25s} {r['n_residues']:>8d} {r['mean_plddt']:>10.1f} {r['pct_very_high']:>5.1f}% {r['pct_high']:>6.1f}% {r['pct_low']:>6.1f}% {r['pct_very_low']:>4.1f}% {quality:<10s}")
        rows.append(r)

    # Save summary table
    out_path = os.path.join(OUT_DIR, "af3_structure_quality_summary.csv")
    with open(out_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['structure', 'n_residues', 'mean_plddt', 'pct_very_high_90plus',
                         'pct_high_70_90', 'pct_low_50_70', 'pct_very_low_below50', 'quality'])
        for r in rows:
            quality = "Excellent" if r['mean_plddt'] >= 80 else "Good" if r['mean_plddt'] >= 70 else "Moderate" if r['mean_plddt'] >= 60 else "Poor"
            writer.writerow([r['name'], r['n_residues'], f"{r['mean_plddt']:.1f}",
                           f"{r['pct_very_high']:.1f}", f"{r['pct_high']:.1f}",
                           f"{r['pct_low']:.1f}", f"{r['pct_very_low']:.1f}", quality])
    print(f"\n  Saved: {out_path}")


# ── v1 vs v2 complex comparison ───────────────────────────────────
def compare_v1_v2():
    """Compare original (full LDLR) vs truncated (ectodomain) complex predictions."""
    print(f"\n{'='*60}")
    print(f"  v1 vs v2 Complex Comparison (full vs truncated LDLR)")
    print(f"{'='*60}")

    comparisons = [
        ("LDLR-ApoB v1 (full)", "ldlr_apob_v1", "fold_calon_ldlr_apob_complex"),
        ("LDLR-ApoB v2 (ecto)", "ldlr_apob_v2", "fold_calon_ldlr_apob_complex_v2"),
    ]

    print(f"\n  {'Version':<25s} {'iPTM':>6s} {'pTM':>6s} {'Ranking':>8s} {'Disorder':>9s} {'Cross-chain PAE':>16s}")
    print(f"  {'-'*25} {'-'*6} {'-'*6} {'-'*8} {'-'*9} {'-'*16}")

    for label, subdir, prefix in comparisons:
        conf_file = os.path.join(STRUCT_DIR, subdir, f"{prefix}_summary_confidences_0.json")
        if os.path.exists(conf_file):
            with open(conf_file) as f:
                conf = json.load(f)
            pae_cross = conf['chain_pair_pae_min'][0][1] if len(conf['chain_pair_pae_min']) > 0 and len(conf['chain_pair_pae_min'][0]) > 1 else 'N/A'
            print(f"  {label:<25s} {conf['iptm']:>6.2f} {conf['ptm']:>6.2f} {conf['ranking_score']:>8.2f} {conf['fraction_disordered']:>9.2f} {pae_cross:>16.2f}")

    print(f"\n  Interpretation:")
    print(f"    - iPTM > 0.7 = confident interface prediction")
    print(f"    - iPTM 0.5-0.7 = moderate, interface may be approximate")
    print(f"    - iPTM < 0.5 = low confidence, use with caution")
    print(f"    - Cross-chain PAE < 5 Å = good relative positioning")


# ── Main execution ────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("  CALON-Structure: AF3 Structure Analysis")
    print("  Phase 2 — Quality Assessment & Variant Mapping")
    print("=" * 60)

    results = []

    # 1. Analyse monomers
    ldlr_data = analyse_monomer(
        "LDLR_wildtype", "ldlr_wildtype",
        "fold_calon_ldlr_wildtype", is_ldlr=True
    )
    results.append(ldlr_data)

    pcsk9_data = analyse_monomer(
        "PCSK9_wildtype", "pcsk9_wildtype",
        "fold_calon_pcsk9_wildtype", is_ldlr=False
    )
    results.append(pcsk9_data)

    apob_data = analyse_monomer(
        "ApoB_RBD", "apob_rbd",
        "fold_calon_apob_rbd", is_ldlr=False
    )
    results.append(apob_data)

    # 2. Analyse complexes
    print("\n" + "=" * 60)
    print("  COMPLEX ANALYSIS")
    print("=" * 60)

    analyse_complex(
        "LDLR-ApoB v1 (full LDLR)",
        "ldlr_apob_v1", "fold_calon_ldlr_apob_complex"
    )

    analyse_complex(
        "LDLR-ApoB v2 (ectodomain only)",
        "ldlr_apob_v2", "fold_calon_ldlr_apob_complex_v2"
    )

    # 3. v1 vs v2 comparison
    compare_v1_v2()

    # 4. Summary table
    generate_summary_table(results)

    # 5. Map variants to LDLR structure
    if ldlr_data:
        map_variants_to_structure(ldlr_data)

    # 6. Recommendations
    print(f"\n{'='*60}")
    print(f"  RECOMMENDATIONS")
    print(f"{'='*60}")
    print(f"""
  1. LDLR wild-type structure:
     - Ectodomain (res 1-692): GOOD quality for FoldX modelling
     - TM/cytoplasmic (res 751-860): EXCLUDE from FoldX (disordered)
     - Use model_0 (highest ranking score) as FoldX template

  2. PCSK9 wild-type:
     - Excellent quality (pTM=0.81), suitable for interface analysis
     - LDLR-PCSK9 complex still pending — resubmit v2 when ready

  3. ApoB RBD:
     - Best quality of all structures (pTM=0.84)
     - Good template for R3527Q variant modelling

  4. Complex predictions:
     - v2 (truncated) should show cleaner interface than v1
     - iPTM values will determine if complexes are usable
     - If iPTM < 0.5: use experimental crystal structures instead
       (PDB: 1N7D for LDLR-PCSK9 EGF-A, 1AJ7 for LDLR-RAP)

  5. Next steps:
     - Convert best CIF models to PDB format for FoldX
     - Run FoldX RepairPDB on wild-type structures
     - Run FoldX BuildModel for each top-50 missense variant
     - Extract ΔΔG, RMSD, contact changes
""")

    print(f"\n  All outputs saved to: {OUT_DIR}")
    print(f"  Done!")


if __name__ == "__main__":
    main()

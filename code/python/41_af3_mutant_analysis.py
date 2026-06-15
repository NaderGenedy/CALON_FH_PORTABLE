#!/usr/bin/env python3
"""
===================================================================================
ALPHAFOLD3 MUTANT vs WILDTYPE STRUCTURAL COMPARISON
Analyses AF3 mutant predictions to validate the wildtype-backbone assumption

After downloading AF3 mutant results, this script:
  1. Extracts pLDDT per-residue from each mutant prediction
  2. Computes backbone RMSD between mutant and wildtype structures
  3. Identifies positions with significant backbone rearrangement
  4. Compares mutant pLDDT vs wildtype pLDDT (local confidence change)
  5. Validates/invalidates the FoldX wildtype-backbone assumption
  6. Generates validation figures and updates the confidence table

Run AFTER downloading all AF3 mutant results to:
  alphafold/af3_mutant_results/

Author: Dr Nader Genedy
===================================================================================
"""

import csv, json, os, math, zipfile
from collections import defaultdict

BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"
AF = f"{BASE}/alphafold/analysis"
RESULTS_DIR = f"{BASE}/alphafold/af3_mutant_results"
WT_PDB = f"{BASE}/alphafold/relaxed.pdb"  # Rosetta-relaxed wildtype

def sf(x):
    if x is None or str(x).strip() == '':
        return None
    try:
        return float(str(x).strip())
    except:
        return None

def parse_pdb_ca(pdb_path):
    """Extract CA atom coordinates and B-factors (pLDDT) from PDB file."""
    atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                try:
                    res_num = int(line[22:26].strip())
                    res_name = line[17:20].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    bfactor = float(line[60:66])  # pLDDT stored as B-factor
                    atoms.append({
                        'res_num': res_num,
                        'res_name': res_name,
                        'x': x, 'y': y, 'z': z,
                        'plddt': bfactor,
                    })
                except (ValueError, IndexError):
                    continue
    return atoms

def compute_rmsd(atoms1, atoms2, window=None):
    """Compute RMSD between two sets of CA atoms. Optionally for a window around a position."""
    if window:
        centre, radius = window
        atoms1 = [a for a in atoms1 if abs(a['res_num'] - centre) <= radius]
        atoms2 = [a for a in atoms2 if abs(a['res_num'] - centre) <= radius]

    # Match by residue number
    pos1 = {a['res_num']: a for a in atoms1}
    pos2 = {a['res_num']: a for a in atoms2}
    common = sorted(set(pos1.keys()) & set(pos2.keys()))

    if len(common) < 3:
        return None, 0

    sq_sum = 0
    for p in common:
        dx = pos1[p]['x'] - pos2[p]['x']
        dy = pos1[p]['y'] - pos2[p]['y']
        dz = pos1[p]['z'] - pos2[p]['z']
        sq_sum += dx*dx + dy*dy + dz*dz

    rmsd = math.sqrt(sq_sum / len(common))
    return rmsd, len(common)

def extract_af3_result(zip_path, output_dir):
    """Extract AF3 result ZIP and find the best model PDB."""
    with zipfile.ZipFile(zip_path, 'r') as z:
        z.extractall(output_dir)

    # Find the ranked model PDBs or CIF files
    pdb_files = []
    for root, dirs, files in os.walk(output_dir):
        for f in files:
            if f.endswith('.pdb') or f.endswith('.cif'):
                pdb_files.append(os.path.join(root, f))

    # Prefer ranked_0 (best model)
    for pf in pdb_files:
        if 'ranked_0' in pf or 'model_0' in pf or 'best' in pf:
            return pf

    # Otherwise return first PDB
    return pdb_files[0] if pdb_files else None

# ─── MAIN ANALYSIS ───
print("=" * 70)
print("  ALPHAFOLD3 MUTANT vs WILDTYPE STRUCTURAL COMPARISON")
print("=" * 70)

# Check if results directory exists
if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    print(f"\n  Results directory created: {RESULTS_DIR}")
    print(f"  Please download AF3 mutant results to this directory.")
    print(f"  Each result should be a ZIP file named like: LDLR_C352Y.zip")
    print(f"\n  After downloading, re-run this script.")

# Check for downloaded results
result_zips = [f for f in os.listdir(RESULTS_DIR) if f.endswith('.zip')] if os.path.exists(RESULTS_DIR) else []
print(f"\n  Found {len(result_zips)} result ZIP files in {RESULTS_DIR}")

if len(result_zips) == 0:
    print("\n  No results found yet. Generating template for expected files...")

    # Load job list
    jobs = []
    job_csv = f"{AF}/af3_mutant_jobs.csv"
    if os.path.exists(job_csv):
        with open(job_csv, encoding='utf-8-sig') as f:
            jobs = list(csv.DictReader(f))

    print(f"\n  Expected {len(jobs)} monomer result files:")
    for j in jobs:
        print(f"    {j['job_name']}.zip")

    print(f"\n  Download each from https://alphafoldserver.com/ after submission")
    print(f"  Save to: {RESULTS_DIR}/")

else:
    print("\n  Processing results...")

    # Load wildtype structure
    wt_atoms = None
    if os.path.exists(WT_PDB):
        wt_atoms = parse_pdb_ca(WT_PDB)
        print(f"  Wildtype structure: {len(wt_atoms)} CA atoms from {WT_PDB}")
    else:
        print(f"  WARNING: Wildtype PDB not found at {WT_PDB}")
        print(f"  Will compare pLDDT only (no RMSD)")

    # Load wildtype pLDDT from CSV
    wt_plddt = {}
    plddt_csv = f"{AF}/LDLR_wildtype_plddt_per_residue.csv"
    with open(plddt_csv, encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            pos = int(r['residue_num'])
            wt_plddt[pos] = float(r['plddt'])

    # Load job list for metadata
    jobs = {}
    job_csv = f"{AF}/af3_mutant_jobs.csv"
    if os.path.exists(job_csv):
        with open(job_csv, encoding='utf-8-sig') as f:
            for r in csv.DictReader(f):
                jobs[r['job_name']] = r

    # Process each result
    results = []
    for zip_file in result_zips:
        job_name = zip_file.replace('.zip', '')
        print(f"\n  Processing: {job_name}")

        # Extract ZIP
        extract_dir = f"{RESULTS_DIR}/{job_name}"
        zip_path = f"{RESULTS_DIR}/{zip_file}"

        try:
            pdb_path = extract_af3_result(zip_path, extract_dir)
            if not pdb_path:
                print(f"    ERROR: No PDB found in {zip_file}")
                continue

            mut_atoms = parse_pdb_ca(pdb_path)
            print(f"    PDB: {os.path.basename(pdb_path)} ({len(mut_atoms)} CA atoms)")

            # Get mutation position from job metadata
            job_meta = jobs.get(job_name, {})
            mut_pos = int(job_meta.get('position', 0))
            mut_name = job_meta.get('mutation', job_name)
            domain = job_meta.get('domain', '')
            foldx_ddg = sf(job_meta.get('foldx_ddG', ''))

            result = {
                'job_name': job_name,
                'mutation': mut_name,
                'position': mut_pos,
                'domain': domain,
                'foldx_ddG': foldx_ddg,
                'n_ca_atoms': len(mut_atoms),
            }

            # 1. pLDDT comparison
            mut_plddt = {a['res_num']: a['plddt'] for a in mut_atoms}
            if mut_pos > 0 and mut_pos in mut_plddt and mut_pos in wt_plddt:
                result['wt_plddt_at_site'] = wt_plddt[mut_pos]
                result['mut_plddt_at_site'] = mut_plddt[mut_pos]
                result['plddt_change'] = mut_plddt[mut_pos] - wt_plddt[mut_pos]
                print(f"    pLDDT at mut site: WT={wt_plddt[mut_pos]:.1f} -> MUT={mut_plddt[mut_pos]:.1f} (delta={result['plddt_change']:+.1f})")

            # Global pLDDT comparison
            common_pos = sorted(set(mut_plddt.keys()) & set(wt_plddt.keys()))
            if common_pos:
                wt_mean = sum(wt_plddt[p] for p in common_pos) / len(common_pos)
                mut_mean = sum(mut_plddt[p] for p in common_pos) / len(common_pos)
                result['wt_mean_plddt'] = wt_mean
                result['mut_mean_plddt'] = mut_mean
                result['global_plddt_change'] = mut_mean - wt_mean
                print(f"    Global pLDDT: WT={wt_mean:.1f} -> MUT={mut_mean:.1f} (delta={result['global_plddt_change']:+.1f})")

                # Count positions with >10 pLDDT drop
                n_affected = sum(1 for p in common_pos if mut_plddt[p] - wt_plddt[p] < -10)
                result['n_positions_affected'] = n_affected
                result['pct_positions_affected'] = n_affected / len(common_pos) * 100
                print(f"    Positions with >10 pLDDT drop: {n_affected} ({result['pct_positions_affected']:.1f}%)")

            # 2. RMSD analysis (if wildtype PDB available)
            if wt_atoms:
                # Global RMSD
                global_rmsd, n_matched = compute_rmsd(wt_atoms, mut_atoms)
                result['global_rmsd'] = global_rmsd
                print(f"    Global backbone RMSD: {global_rmsd:.3f} A ({n_matched} atoms matched)")

                # Local RMSD (10-residue window around mutation)
                if mut_pos > 0:
                    local_rmsd, n_local = compute_rmsd(wt_atoms, mut_atoms, window=(mut_pos, 10))
                    result['local_rmsd'] = local_rmsd
                    result['local_rmsd_n'] = n_local
                    print(f"    Local RMSD (+-10 residues): {local_rmsd:.3f} A ({n_local} atoms)")

                    # Extended local RMSD (25-residue window)
                    ext_rmsd, n_ext = compute_rmsd(wt_atoms, mut_atoms, window=(mut_pos, 25))
                    result['extended_rmsd'] = ext_rmsd
                    print(f"    Extended RMSD (+-25 residues): {ext_rmsd:.3f} A ({n_ext} atoms)")

                # Classify backbone change
                if global_rmsd is not None:
                    if global_rmsd < 1.0:
                        result['backbone_status'] = 'PRESERVED'
                        result['foldx_valid'] = 'Yes'
                    elif global_rmsd < 2.0:
                        result['backbone_status'] = 'MINOR_CHANGE'
                        result['foldx_valid'] = 'Partial'
                    else:
                        result['backbone_status'] = 'MAJOR_REARRANGEMENT'
                        result['foldx_valid'] = 'No'
                    print(f"    Backbone status: {result['backbone_status']} (FoldX valid: {result['foldx_valid']})")

            results.append(result)

        except Exception as e:
            print(f"    ERROR: {e}")
            continue

    # ─── Save results ───
    if results:
        output_csv = f"{AF}/af3_mutant_validation_results.csv"
        fieldnames = list(results[0].keys())
        # Ensure all results have all keys
        all_keys = set()
        for r in results:
            all_keys.update(r.keys())
        fieldnames = sorted(all_keys)

        with open(output_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in results:
                clean = {k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()}
                writer.writerow(clean)

        print(f"\n{'='*70}")
        print(f"  RESULTS SAVED: {output_csv}")
        print(f"  Processed: {len(results)} mutant structures")
        print(f"{'='*70}")

        # Summary statistics
        preserved = sum(1 for r in results if r.get('backbone_status') == 'PRESERVED')
        minor = sum(1 for r in results if r.get('backbone_status') == 'MINOR_CHANGE')
        major = sum(1 for r in results if r.get('backbone_status') == 'MAJOR_REARRANGEMENT')
        foldx_valid = sum(1 for r in results if r.get('foldx_valid') == 'Yes')

        print(f"\n  BACKBONE VALIDATION SUMMARY:")
        print(f"    Preserved (RMSD < 1.0 A):        {preserved:3d} ({preserved/len(results)*100:.1f}%)")
        print(f"    Minor change (1.0-2.0 A):         {minor:3d} ({minor/len(results)*100:.1f}%)")
        print(f"    Major rearrangement (> 2.0 A):    {major:3d} ({major/len(results)*100:.1f}%)")
        print(f"    FoldX assumption validated:        {foldx_valid:3d} ({foldx_valid/len(results)*100:.1f}%)")

        # pLDDT analysis
        plddt_drops = [r.get('plddt_change', 0) for r in results if r.get('plddt_change') is not None]
        if plddt_drops:
            print(f"\n  pLDDT CHANGE AT MUTATION SITE:")
            print(f"    Mean change: {sum(plddt_drops)/len(plddt_drops):+.1f}")
            print(f"    Range: {min(plddt_drops):+.1f} to {max(plddt_drops):+.1f}")
            n_confident = sum(1 for d in plddt_drops if d > -10)
            print(f"    Still confident (drop < 10): {n_confident}/{len(plddt_drops)} ({n_confident/len(plddt_drops)*100:.1f}%)")

        # Correlation: FoldX ddG vs backbone RMSD
        ddg_rmsd_pairs = [(r['foldx_ddG'], r['global_rmsd'])
                          for r in results
                          if r.get('foldx_ddG') is not None and r.get('global_rmsd') is not None]
        if len(ddg_rmsd_pairs) >= 5:
            ddg_vals = [p[0] for p in ddg_rmsd_pairs]
            rmsd_vals = [p[1] for p in ddg_rmsd_pairs]
            n = len(ddg_vals)
            mx = sum(ddg_vals) / n
            my = sum(rmsd_vals) / n
            sxy = sum((x-mx)*(y-my) for x, y in zip(ddg_vals, rmsd_vals))
            sxx = sum((x-mx)**2 for x in ddg_vals)
            syy = sum((y-my)**2 for y in rmsd_vals)
            if sxx > 0 and syy > 0:
                r_val = sxy / (sxx * syy) ** 0.5
                print(f"\n  FoldX ddG vs backbone RMSD: r = {r_val:+.3f} (n={n})")
                if abs(r_val) < 0.3:
                    print(f"    Interpretation: Weak correlation - backbone changes do NOT scale with ddG")
                    print(f"    This VALIDATES the FoldX approach for most mutations")
                elif r_val > 0.3:
                    print(f"    Interpretation: Positive correlation - higher ddG = more backbone change")
                    print(f"    Caution: FoldX may underestimate destabilisation for severe mutations")

        # Update confidence table
        print(f"\n  Updating confidence table with AF3 mutant data...")
        conf_table = f"{AF}/LDLR_860_position_confidence_table.csv"
        if os.path.exists(conf_table):
            conf_rows = []
            with open(conf_table, encoding='utf-8-sig') as f:
                reader = csv.DictReader(f)
                conf_fieldnames = reader.fieldnames + ['af3_mutant_plddt', 'af3_backbone_rmsd', 'af3_backbone_status', 'af3_foldx_valid']
                for r in reader:
                    pos = int(r['position'])
                    # Check if we have AF3 mutant data for this position
                    for res in results:
                        if res.get('position') == pos:
                            r['af3_mutant_plddt'] = f"{res.get('mut_plddt_at_site', '')}"
                            r['af3_backbone_rmsd'] = f"{res.get('local_rmsd', '')}"
                            r['af3_backbone_status'] = res.get('backbone_status', '')
                            r['af3_foldx_valid'] = res.get('foldx_valid', '')
                            # Upgrade confidence tier if AF3 validated
                            if res.get('foldx_valid') == 'Yes':
                                current_tier = r.get('confidence_tier', '')
                                if current_tier in ['MODERATE', 'MODERATE+']:
                                    r['confidence_tier'] = 'HIGH'
                                    r['evidence_basis'] += ' + AF3 mutant validated'
                                elif current_tier in ['HIGH']:
                                    r['confidence_tier'] = 'HIGH+'
                                    r['evidence_basis'] += ' + AF3 mutant validated'
                            break
                    else:
                        r['af3_mutant_plddt'] = ''
                        r['af3_backbone_rmsd'] = ''
                        r['af3_backbone_status'] = ''
                        r['af3_foldx_valid'] = ''
                    conf_rows.append(r)

            # Write updated confidence table
            with open(conf_table, 'w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=conf_fieldnames)
                writer.writeheader()
                writer.writerows(conf_rows)
            print(f"    Updated: {conf_table}")

    else:
        print("\n  No results processed. Check that ZIP files are in the correct directory.")

print(f"\n  Done.")

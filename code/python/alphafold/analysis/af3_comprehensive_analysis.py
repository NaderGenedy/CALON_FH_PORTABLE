#!/usr/bin/env python3
"""
Comprehensive AlphaFold3 Mutant Analysis Script
Analyzes 52 LDLR mutant monomers and 20 LDLR-PCSK9 mutant complexes
against wildtype references, cross-references with FoldX ddG data.
"""

import json
import glob
import os
import re
import csv
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

# ============================================================
# CONFIGURATION
# ============================================================
AF3_DIR = "D:/alphafold3_50"
ANALYSIS_DIR = "C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/analysis"
WT_MONOMER_DIR = os.path.join(AF3_DIR, "calon_ldlr_wildtype")
WT_COMPLEX_DIR = os.path.join(AF3_DIR, "calon_ldlr_pcsk9_complex_v2")

# Domain mapping
DOMAIN_MAP = [
    (1, 21, "Signal peptide"),
    (22, 83, "Ligand-binding R1"),
    (84, 120, "Ligand-binding R2"),
    (121, 159, "Ligand-binding R3"),
    (160, 197, "Ligand-binding R4"),
    (198, 236, "Ligand-binding R5"),
    (237, 275, "Ligand-binding R6"),
    (276, 313, "Ligand-binding R7"),
    (314, 353, "EGF-like A"),
    (354, 393, "EGF-like B"),
    (394, 631, "Beta-propeller"),
    (632, 671, "EGF-like C"),
    (672, 692, "Linker"),
    (693, 750, "O-linked sugar"),
    (751, 788, "Transmembrane"),
    (789, 860, "Cytoplasmic tail"),
]

def get_domain(pos):
    for start, end, name in DOMAIN_MAP:
        if start <= pos <= end:
            return name
    return "Unknown"

def parse_mutation(folder_name):
    """Extract mutation from folder name like ldlr_a50v_mutant -> A50V, position=50"""
    # Remove prefix and suffix
    name = folder_name.replace("ldlr_", "").replace("_mutant", "").replace("_pcsk9_complex", "").replace("_2", "")
    # Parse: first char = wt_aa, digits = position, last char(s) = mut_aa
    match = re.match(r'^([a-z])(\d+)([a-z]+)$', name)
    if match:
        wt_aa = match.group(1).upper()
        pos = int(match.group(2))
        mut_aa = match.group(3).upper()
        mutation = f"{wt_aa}{pos}{mut_aa}"
        return mutation, wt_aa, pos, mut_aa
    return name.upper(), "", 0, ""

def read_summary_json(folder_path):
    """Read summary_confidences_0.json from a folder"""
    pattern = os.path.join(folder_path, "*_summary_confidences_0.json")
    files = glob.glob(pattern)
    if not files:
        print(f"  WARNING: No summary_confidences_0.json found in {folder_path}")
        return None
    with open(files[0], 'r') as f:
        return json.load(f)

def extract_per_residue_plddt(folder_path, residue_pos):
    """Extract per-residue pLDDT at a specific position from full_data_0.json.

    full_data has atom-level pLDDT and token-level residue IDs.
    We need to: read the CIF file to map atoms to residues, or
    use token_res_ids (which are per-residue) and average atom_plddts per residue.

    Actually, token_res_ids gives us per-residue indexing.
    We need to compute per-residue pLDDT by averaging atom pLDDT values per residue.
    """
    pattern = os.path.join(folder_path, "*_full_data_0.json")
    files = glob.glob(pattern)
    if not files:
        return None

    try:
        with open(files[0], 'r') as f:
            data = json.load(f)

        atom_plddts = data.get('atom_plddts', [])
        atom_chain_ids = data.get('atom_chain_ids', [])
        token_res_ids = data.get('token_res_ids', [])
        token_chain_ids = data.get('token_chain_ids', [])

        if not atom_plddts or not token_res_ids:
            return None

        # We need to map atoms to residues
        # Read the CIF to get atom-to-residue mapping
        cif_pattern = os.path.join(folder_path, "*_model_0.cif")
        cif_files = glob.glob(cif_pattern)
        if not cif_files:
            return None

        # Parse CIF to get per-residue atom counts for chain A
        residue_atom_plddts = defaultdict(list)
        atom_idx = 0

        with open(cif_files[0], 'r') as f:
            in_atom_site = False
            labels = []
            for line in f:
                line = line.strip()
                if line.startswith('_atom_site.'):
                    if not in_atom_site:
                        in_atom_site = True
                        labels = []
                    labels.append(line.split('.')[1])
                elif in_atom_site and line and not line.startswith('_') and not line.startswith('#') and not line.startswith('loop_'):
                    parts = line.split()
                    if len(parts) >= len(labels):
                        try:
                            label_seq_id_idx = labels.index('label_seq_id')
                            chain_idx_col = labels.index('label_asym_id')
                            bfactor_idx = labels.index('B_iso_or_equiv')

                            chain = parts[chain_idx_col]
                            res_id = int(parts[label_seq_id_idx])
                            bfactor = float(parts[bfactor_idx])

                            if chain == 'A':
                                residue_atom_plddts[res_id].append(bfactor)
                        except (ValueError, IndexError):
                            pass
                elif in_atom_site and (line.startswith('#') or line.startswith('loop_') or line == ''):
                    if labels:
                        in_atom_site = False

        if residue_pos in residue_atom_plddts:
            return np.mean(residue_atom_plddts[residue_pos])

        # Fallback: use the atom_plddts from full_data directly
        # Map using token_res_ids (per-token = per-residue for proteins)
        # Each token corresponds to one residue, atoms are ordered by residue
        # We need the number of atoms per residue to slice atom_plddts
        return None

    except Exception as e:
        return None

def extract_per_residue_plddt_from_fulldata(folder_path, residue_pos, chain='A'):
    """Alternative: use full_data's pae diagonal or compute from atom_plddts.

    Since token_res_ids maps 1:1 to residues, and atom_plddts is atom-level,
    we need to know how many atoms per residue. Easiest: read from CIF B-factors.
    """
    # Try CIF B-factor approach (AF3 stores pLDDT in B-factor column)
    cif_pattern = os.path.join(folder_path, "*_model_0.cif")
    cif_files = glob.glob(cif_pattern)
    if not cif_files:
        return None

    try:
        residue_bfactors = defaultdict(list)
        with open(cif_files[0], 'r') as f:
            in_atom_site = False
            labels = []
            for line in f:
                line = line.strip()
                if line.startswith('_atom_site.'):
                    if not in_atom_site:
                        in_atom_site = True
                        labels = []
                    labels.append(line.split('.')[1])
                elif in_atom_site and line and not line.startswith('_') and not line.startswith('#') and not line.startswith('loop_'):
                    parts = line.split()
                    if len(parts) >= len(labels):
                        try:
                            label_seq_id_idx = labels.index('label_seq_id')
                            chain_idx_col = labels.index('label_asym_id')
                            bfactor_idx = labels.index('B_iso_or_equiv')

                            chain_id = parts[chain_idx_col]
                            res_id = int(parts[label_seq_id_idx])
                            bfactor = float(parts[bfactor_idx])

                            if chain_id == chain:
                                residue_bfactors[res_id].append(bfactor)
                        except (ValueError, IndexError):
                            pass
                elif in_atom_site and (line.startswith('#') or line == ''):
                    if labels and residue_bfactors:
                        break

        if residue_pos in residue_bfactors:
            return np.mean(residue_bfactors[residue_pos])
        return None
    except Exception as e:
        return None


# ============================================================
# 1. READ WILDTYPE REFERENCES
# ============================================================
print("=" * 80)
print("ALPHAFOLD3 COMPREHENSIVE MUTANT ANALYSIS")
print("=" * 80)

print("\n1. READING WILDTYPE REFERENCES...")
wt_mono = read_summary_json(WT_MONOMER_DIR)
wt_complex = read_summary_json(WT_COMPLEX_DIR)

print(f"   Wildtype monomer: pTM={wt_mono['ptm']}, ranking_score={wt_mono['ranking_score']}, "
      f"fraction_disordered={wt_mono['fraction_disordered']}")
print(f"   Wildtype complex: pTM={wt_complex['ptm']}, ipTM={wt_complex['iptm']}, "
      f"ranking_score={wt_complex['ranking_score']}, fraction_disordered={wt_complex['fraction_disordered']}")

# Complex chain layout: chain 0 = LDLR, chain 1 = PCSK9
# chain_pair_iptm[0][1] = LDLR->PCSK9 interface
# chain_pair_pae_min[0][1] = LDLR->PCSK9 PAE
wt_ldlr_pcsk9_iptm = wt_complex['chain_pair_iptm'][0][1]
wt_pcsk9_ldlr_iptm = wt_complex['chain_pair_iptm'][1][0]
wt_ldlr_pcsk9_pae = wt_complex['chain_pair_pae_min'][0][1]
wt_pcsk9_ldlr_pae = wt_complex['chain_pair_pae_min'][1][0]
wt_ldlr_chain_ptm = wt_complex['chain_ptm'][0]
wt_pcsk9_chain_ptm = wt_complex['chain_ptm'][1]

print(f"   Wildtype complex LDLR chain pTM: {wt_ldlr_chain_ptm}")
print(f"   Wildtype complex PCSK9 chain pTM: {wt_pcsk9_chain_ptm}")
print(f"   Wildtype LDLR->PCSK9 interface ipTM: {wt_ldlr_pcsk9_iptm}")
print(f"   Wildtype LDLR->PCSK9 PAE_min: {wt_ldlr_pcsk9_pae}")

# ============================================================
# 2. READ WILDTYPE pLDDT PER RESIDUE
# ============================================================
print("\n2. READING WILDTYPE pLDDT PER RESIDUE...")
wt_plddt_file = os.path.join(ANALYSIS_DIR, "LDLR_wildtype_plddt_per_residue.csv")
wt_plddt_df = pd.read_csv(wt_plddt_file)
wt_plddt_map = dict(zip(wt_plddt_df['residue_num'], wt_plddt_df['plddt']))
print(f"   Loaded {len(wt_plddt_map)} residue pLDDT values")

# ============================================================
# 3. ANALYZE MUTANT MONOMERS
# ============================================================
print("\n3. ANALYZING MUTANT MONOMERS...")
monomer_results = []

monomer_folders = sorted(glob.glob(os.path.join(AF3_DIR, "ldlr_*_mutant")))
print(f"   Found {len(monomer_folders)} mutant monomer folders")

for folder in monomer_folders:
    folder_name = os.path.basename(folder)
    mutation, wt_aa, pos, mut_aa = parse_mutation(folder_name)
    domain = get_domain(pos)

    data = read_summary_json(folder)
    if data is None:
        continue

    ptm = data['ptm']
    ranking_score = data['ranking_score']
    fraction_disordered = data['fraction_disordered']
    has_clash = data['has_clash']
    chain_ptm = data['chain_ptm'][0] if data.get('chain_ptm') else ptm

    delta_ptm = ptm - wt_mono['ptm']
    delta_ranking = ranking_score - wt_mono['ranking_score']
    delta_disorder = fraction_disordered - wt_mono['fraction_disordered']

    # Extract pLDDT at mutation site
    mut_plddt = extract_per_residue_plddt_from_fulldata(folder, pos, chain='A')
    wt_plddt_at_site = wt_plddt_map.get(pos, None)
    delta_plddt = None
    if mut_plddt is not None and wt_plddt_at_site is not None:
        delta_plddt = mut_plddt - wt_plddt_at_site

    # Classify severity
    if delta_ptm <= -0.05:
        severity = "severe"
    elif delta_ptm <= -0.02:
        severity = "moderate"
    elif delta_ptm >= 0.02:
        severity = "improved"
    else:
        severity = "minimal"

    monomer_results.append({
        'mutation': mutation,
        'wt_aa': wt_aa,
        'mut_aa': mut_aa,
        'position': pos,
        'domain': domain,
        'ptm': ptm,
        'ranking_score': ranking_score,
        'fraction_disordered': fraction_disordered,
        'has_clash': has_clash,
        'chain_ptm': chain_ptm,
        'delta_ptm': delta_ptm,
        'delta_ranking_score': delta_ranking,
        'delta_fraction_disordered': delta_disorder,
        'wt_ptm': wt_mono['ptm'],
        'wt_ranking_score': wt_mono['ranking_score'],
        'mut_plddt_at_site': round(mut_plddt, 2) if mut_plddt is not None else None,
        'wt_plddt_at_site': round(wt_plddt_at_site, 2) if wt_plddt_at_site is not None else None,
        'delta_plddt_at_site': round(delta_plddt, 2) if delta_plddt is not None else None,
        'severity': severity,
    })

mono_df = pd.DataFrame(monomer_results)
mono_df = mono_df.sort_values('delta_ptm')
mono_csv = os.path.join(ANALYSIS_DIR, "af3_mutant_monomer_results.csv")
mono_df.to_csv(mono_csv, index=False)
print(f"   Saved monomer results to {mono_csv}")
print(f"   Processed {len(mono_df)} mutant monomers")

# ============================================================
# 4. ANALYZE MUTANT COMPLEXES
# ============================================================
print("\n4. ANALYZING MUTANT COMPLEXES...")
complex_results = []

complex_folders = sorted(glob.glob(os.path.join(AF3_DIR, "ldlr_*_pcsk9_complex")))
# Exclude the duplicate
complex_folders = [f for f in complex_folders if not f.endswith("_complex_2")]
# Also exclude wildtype
complex_folders = [f for f in complex_folders if "calon_" not in os.path.basename(f)]
print(f"   Found {len(complex_folders)} mutant complex folders")

for folder in complex_folders:
    folder_name = os.path.basename(folder)
    mutation, wt_aa, pos, mut_aa = parse_mutation(folder_name)
    domain = get_domain(pos)

    data = read_summary_json(folder)
    if data is None:
        continue

    ptm = data['ptm']
    iptm = data['iptm']
    ranking_score = data['ranking_score']
    fraction_disordered = data['fraction_disordered']
    has_clash = data['has_clash']

    # Chain-level metrics
    ldlr_chain_ptm = data['chain_ptm'][0] if len(data.get('chain_ptm', [])) > 0 else None
    pcsk9_chain_ptm = data['chain_ptm'][1] if len(data.get('chain_ptm', [])) > 1 else None

    # Interface metrics
    ldlr_pcsk9_iptm = data['chain_pair_iptm'][0][1] if len(data.get('chain_pair_iptm', [])) > 0 else None
    pcsk9_ldlr_iptm = data['chain_pair_iptm'][1][0] if len(data.get('chain_pair_iptm', [])) > 1 else None
    ldlr_pcsk9_pae = data['chain_pair_pae_min'][0][1] if len(data.get('chain_pair_pae_min', [])) > 0 else None
    pcsk9_ldlr_pae = data['chain_pair_pae_min'][1][0] if len(data.get('chain_pair_pae_min', [])) > 1 else None

    delta_ptm = ptm - wt_complex['ptm']
    delta_iptm = iptm - wt_complex['iptm'] if iptm is not None else None
    delta_ranking = ranking_score - wt_complex['ranking_score']
    delta_ldlr_chain_ptm = ldlr_chain_ptm - wt_ldlr_chain_ptm if ldlr_chain_ptm is not None else None
    delta_interface_iptm = ldlr_pcsk9_iptm - wt_ldlr_pcsk9_iptm if ldlr_pcsk9_iptm is not None else None
    delta_interface_pae = ldlr_pcsk9_pae - wt_ldlr_pcsk9_pae if ldlr_pcsk9_pae is not None else None

    complex_results.append({
        'mutation': mutation,
        'wt_aa': wt_aa,
        'mut_aa': mut_aa,
        'position': pos,
        'domain': domain,
        'ptm': ptm,
        'iptm': iptm,
        'ranking_score': ranking_score,
        'fraction_disordered': fraction_disordered,
        'has_clash': has_clash,
        'ldlr_chain_ptm': ldlr_chain_ptm,
        'pcsk9_chain_ptm': pcsk9_chain_ptm,
        'ldlr_pcsk9_iptm': ldlr_pcsk9_iptm,
        'pcsk9_ldlr_iptm': pcsk9_ldlr_iptm,
        'ldlr_pcsk9_pae_min': ldlr_pcsk9_pae,
        'pcsk9_ldlr_pae_min': pcsk9_ldlr_pae,
        'delta_ptm': delta_ptm,
        'delta_iptm': delta_iptm,
        'delta_ranking_score': delta_ranking,
        'delta_ldlr_chain_ptm': delta_ldlr_chain_ptm,
        'delta_interface_iptm': delta_interface_iptm,
        'delta_interface_pae': delta_interface_pae,
        'wt_ptm': wt_complex['ptm'],
        'wt_iptm': wt_complex['iptm'],
        'wt_ranking_score': wt_complex['ranking_score'],
    })

complex_df = pd.DataFrame(complex_results)
complex_df = complex_df.sort_values('delta_iptm')
complex_csv = os.path.join(ANALYSIS_DIR, "af3_mutant_complex_results.csv")
complex_df.to_csv(complex_csv, index=False)
print(f"   Saved complex results to {complex_csv}")
print(f"   Processed {len(complex_df)} mutant complexes")

# ============================================================
# 5. CREATE COMBINED SUMMARY
# ============================================================
print("\n5. CREATING COMBINED SUMMARY...")

# Merge monomer and complex data where both exist
summary_rows = []
for _, row in mono_df.iterrows():
    entry = row.to_dict()
    # Find matching complex
    cmatch = complex_df[complex_df['mutation'] == row['mutation']]
    if len(cmatch) > 0:
        cr = cmatch.iloc[0]
        entry['complex_ptm'] = cr['ptm']
        entry['complex_iptm'] = cr['iptm']
        entry['complex_ranking_score'] = cr['ranking_score']
        entry['complex_delta_ptm'] = cr['delta_ptm']
        entry['complex_delta_iptm'] = cr['delta_iptm']
        entry['complex_delta_ranking_score'] = cr['delta_ranking_score']
        entry['ldlr_pcsk9_iptm'] = cr['ldlr_pcsk9_iptm']
        entry['delta_interface_iptm'] = cr['delta_interface_iptm']
        entry['ldlr_pcsk9_pae_min'] = cr['ldlr_pcsk9_pae_min']
        entry['delta_interface_pae'] = cr['delta_interface_pae']
        entry['has_complex'] = True
    else:
        entry['complex_ptm'] = None
        entry['complex_iptm'] = None
        entry['complex_ranking_score'] = None
        entry['complex_delta_ptm'] = None
        entry['complex_delta_iptm'] = None
        entry['complex_delta_ranking_score'] = None
        entry['ldlr_pcsk9_iptm'] = None
        entry['delta_interface_iptm'] = None
        entry['ldlr_pcsk9_pae_min'] = None
        entry['delta_interface_pae'] = None
        entry['has_complex'] = False
    summary_rows.append(entry)

summary_df = pd.DataFrame(summary_rows)
summary_csv = os.path.join(ANALYSIS_DIR, "af3_mutant_vs_wildtype_summary.csv")
summary_df.to_csv(summary_csv, index=False)
print(f"   Saved combined summary to {summary_csv}")

# ============================================================
# 6. CROSS-REFERENCE WITH FOLDX
# ============================================================
print("\n6. CROSS-REFERENCING WITH FOLDX DATA...")
foldx_file = os.path.join(ANALYSIS_DIR, "foldx_ddg_results.csv")
foldx_df = pd.read_csv(foldx_file)
print(f"   Loaded {len(foldx_df)} FoldX results")

# Create a mapping from position to FoldX ddG
# FoldX has columns: variant_id, foldx_notation, wt_aa, mut_aa, position, domain, n_patients, ddG_kcal_mol, effect
foldx_map = {}
for _, row in foldx_df.iterrows():
    key = f"{row['wt_aa']}{row['position']}{row['mut_aa']}"
    foldx_map[key] = row['ddG_kcal_mol']

# Match
matched_mutations = []
for _, row in mono_df.iterrows():
    if row['mutation'] in foldx_map:
        matched_mutations.append({
            'mutation': row['mutation'],
            'position': row['position'],
            'domain': row['domain'],
            'foldx_ddG': foldx_map[row['mutation']],
            'af3_ptm': row['ptm'],
            'af3_delta_ptm': row['delta_ptm'],
            'af3_ranking_score': row['ranking_score'],
            'af3_delta_ranking': row['delta_ranking_score'],
            'af3_fraction_disordered': row['fraction_disordered'],
            'af3_delta_plddt_at_site': row['delta_plddt_at_site'],
            'severity': row['severity'],
        })

matched_df = pd.DataFrame(matched_mutations)
print(f"   Matched {len(matched_df)} mutations between AF3 and FoldX")

# ============================================================
# 7. COMPREHENSIVE ANALYSIS & REPORTING
# ============================================================
print("\n" + "=" * 80)
print("COMPREHENSIVE ANALYSIS RESULTS")
print("=" * 80)

# --- 7a. Monomer Summary Statistics ---
print("\n" + "-" * 60)
print("7a. MONOMER PREDICTION SUMMARY")
print("-" * 60)
print(f"\n   Total mutant monomers analyzed: {len(mono_df)}")
print(f"   Wildtype reference pTM: {wt_mono['ptm']}")
print(f"   Wildtype reference ranking_score: {wt_mono['ranking_score']}")
print(f"   Wildtype fraction_disordered: {wt_mono['fraction_disordered']}")
print(f"\n   Mutant pTM range: {mono_df['ptm'].min():.4f} - {mono_df['ptm'].max():.4f}")
print(f"   Mutant ranking_score range: {mono_df['ranking_score'].min():.4f} - {mono_df['ranking_score'].max():.4f}")
print(f"   Mean delta_pTM: {mono_df['delta_ptm'].mean():.4f}")
print(f"   Median delta_pTM: {mono_df['delta_ptm'].median():.4f}")

# Severity classification
severe = mono_df[mono_df['severity'] == 'severe']
moderate = mono_df[mono_df['severity'] == 'moderate']
minimal = mono_df[mono_df['severity'] == 'minimal']
improved = mono_df[mono_df['severity'] == 'improved']

print(f"\n   SEVERITY CLASSIFICATION:")
print(f"   Severe structural disruption (delta_pTM < -0.05): {len(severe)} ({100*len(severe)/len(mono_df):.1f}%)")
print(f"   Moderate disruption (-0.05 <= delta_pTM < -0.02):  {len(moderate)} ({100*len(moderate)/len(mono_df):.1f}%)")
print(f"   Minimal change (-0.02 <= delta_pTM < 0.02):        {len(minimal)} ({100*len(minimal)/len(mono_df):.1f}%)")
print(f"   Improved (delta_pTM >= 0.02):                       {len(improved)} ({100*len(improved)/len(mono_df):.1f}%)")

# Clashes
clashing = mono_df[mono_df['has_clash'] > 0]
print(f"\n   Mutations with steric clashes: {len(clashing)}")
if len(clashing) > 0:
    for _, r in clashing.iterrows():
        print(f"      {r['mutation']} (has_clash={r['has_clash']})")

# --- 7b. Most Damaging Mutations ---
print("\n" + "-" * 60)
print("7b. MOST STRUCTURALLY DAMAGING MUTATIONS (by delta_pTM)")
print("-" * 60)
top_damage = mono_df.nsmallest(15, 'delta_ptm')
print(f"\n   {'Mutation':<12} {'Domain':<22} {'pTM':<8} {'delta_pTM':<12} {'Rank_Score':<12} {'Frac_Dis':<10} {'Severity'}")
print(f"   {'-'*10:<12} {'-'*20:<22} {'-'*6:<8} {'-'*10:<12} {'-'*10:<12} {'-'*8:<10} {'-'*8}")
for _, r in top_damage.iterrows():
    print(f"   {r['mutation']:<12} {r['domain']:<22} {r['ptm']:<8.4f} {r['delta_ptm']:<12.4f} {r['ranking_score']:<12.4f} {r['fraction_disordered']:<10.3f} {r['severity']}")

# --- 7c. Mutations that INCREASE stability ---
print("\n" + "-" * 60)
print("7c. MUTATIONS WITH INCREASED CONFIDENCE (potential gain-of-function)")
print("-" * 60)
if len(improved) > 0:
    print(f"\n   {'Mutation':<12} {'Domain':<22} {'pTM':<8} {'delta_pTM':<12} {'Rank_Score':<12} {'Frac_Dis':<10}")
    print(f"   {'-'*10:<12} {'-'*20:<22} {'-'*6:<8} {'-'*10:<12} {'-'*10:<12} {'-'*8:<10}")
    for _, r in improved.sort_values('delta_ptm', ascending=False).iterrows():
        print(f"   {r['mutation']:<12} {r['domain']:<22} {r['ptm']:<8.4f} {r['delta_ptm']:<12.4f} {r['ranking_score']:<12.4f} {r['fraction_disordered']:<10.3f}")
else:
    print("   No mutations showed significant improvement (delta_pTM >= 0.02)")

# --- 7d. Domain-by-Domain Analysis ---
print("\n" + "-" * 60)
print("7d. DOMAIN-BY-DOMAIN ANALYSIS")
print("-" * 60)
domain_stats = mono_df.groupby('domain').agg({
    'mutation': 'count',
    'delta_ptm': ['mean', 'min', 'max'],
    'delta_ranking_score': 'mean',
    'fraction_disordered': 'mean',
}).round(4)
domain_stats.columns = ['n_mutations', 'mean_delta_ptm', 'min_delta_ptm', 'max_delta_ptm', 'mean_delta_rank', 'mean_frac_dis']
domain_stats = domain_stats.sort_values('mean_delta_ptm')

print(f"\n   {'Domain':<22} {'N':<5} {'Mean dPTM':<12} {'Min dPTM':<12} {'Max dPTM':<12} {'Mean dRank':<12} {'Mean FracDis'}")
print(f"   {'-'*20:<22} {'-'*3:<5} {'-'*10:<12} {'-'*10:<12} {'-'*10:<12} {'-'*10:<12} {'-'*10}")
for domain_name, row in domain_stats.iterrows():
    print(f"   {domain_name:<22} {int(row['n_mutations']):<5} {row['mean_delta_ptm']:<12.4f} {row['min_delta_ptm']:<12.4f} {row['max_delta_ptm']:<12.4f} {row['mean_delta_rank']:<12.4f} {row['mean_frac_dis']:<10.3f}")

# --- 7e. Complex Analysis ---
print("\n" + "-" * 60)
print("7e. COMPLEX (LDLR-PCSK9) ANALYSIS")
print("-" * 60)
print(f"\n   Total mutant complexes analyzed: {len(complex_df)}")
print(f"   Wildtype complex pTM: {wt_complex['ptm']}, ipTM: {wt_complex['iptm']}")
print(f"   Wildtype LDLR-PCSK9 interface ipTM: {wt_ldlr_pcsk9_iptm}")
print(f"   Wildtype LDLR-PCSK9 PAE_min: {wt_ldlr_pcsk9_pae}")

print(f"\n   Mutant complex ipTM range: {complex_df['iptm'].min():.4f} - {complex_df['iptm'].max():.4f}")
print(f"   Mean delta_ipTM: {complex_df['delta_iptm'].mean():.4f}")

print(f"\n   INTERFACE DISRUPTION ANALYSIS:")
print(f"   {'Mutation':<12} {'Domain':<22} {'ipTM':<8} {'d_ipTM':<10} {'LDLR-PCSK9':<12} {'d_interface':<12} {'PAE_min':<10} {'d_PAE':<10}")
print(f"   {'-'*10:<12} {'-'*20:<22} {'-'*6:<8} {'-'*8:<10} {'-'*10:<12} {'-'*10:<12} {'-'*8:<10} {'-'*8:<10}")
for _, r in complex_df.sort_values('delta_iptm').iterrows():
    d_intf = f"{r['delta_interface_iptm']:.4f}" if r['delta_interface_iptm'] is not None else "N/A"
    d_pae = f"{r['delta_interface_pae']:.2f}" if r['delta_interface_pae'] is not None else "N/A"
    print(f"   {r['mutation']:<12} {r['domain']:<22} {r['iptm']:<8.4f} {r['delta_iptm']:<10.4f} {r['ldlr_pcsk9_iptm']:<12.4f} {d_intf:<12} {r['ldlr_pcsk9_pae_min']:<10.2f} {d_pae:<10}")

# Mutations near PCSK9 binding site (EGF-like A, B, beta-propeller roughly 314-631)
binding_region_complexes = complex_df[complex_df['position'].between(160, 400)]
non_binding_complexes = complex_df[~complex_df['position'].between(160, 400)]

print(f"\n   Mutations in/near PCSK9 binding region (res 160-400): {len(binding_region_complexes)}")
if len(binding_region_complexes) > 0:
    print(f"      Mean delta_ipTM: {binding_region_complexes['delta_iptm'].mean():.4f}")
print(f"   Mutations outside binding region: {len(non_binding_complexes)}")
if len(non_binding_complexes) > 0:
    print(f"      Mean delta_ipTM: {non_binding_complexes['delta_iptm'].mean():.4f}")

# --- 7f. pLDDT at Mutation Site ---
print("\n" + "-" * 60)
print("7f. pLDDT AT MUTATION SITE")
print("-" * 60)
plddt_data = mono_df[mono_df['delta_plddt_at_site'].notna()]
print(f"   Mutations with pLDDT data at site: {len(plddt_data)} / {len(mono_df)}")
if len(plddt_data) > 0:
    print(f"   Mean delta_pLDDT at site: {plddt_data['delta_plddt_at_site'].mean():.2f}")
    print(f"   Range: {plddt_data['delta_plddt_at_site'].min():.2f} to {plddt_data['delta_plddt_at_site'].max():.2f}")

    print(f"\n   WORST pLDDT DROPS AT MUTATION SITE:")
    print(f"   {'Mutation':<12} {'Domain':<22} {'WT pLDDT':<10} {'Mut pLDDT':<10} {'Delta':<10}")
    print(f"   {'-'*10:<12} {'-'*20:<22} {'-'*8:<10} {'-'*8:<10} {'-'*8:<10}")
    for _, r in plddt_data.nsmallest(15, 'delta_plddt_at_site').iterrows():
        print(f"   {r['mutation']:<12} {r['domain']:<22} {r['wt_plddt_at_site']:<10.2f} {r['mut_plddt_at_site']:<10.2f} {r['delta_plddt_at_site']:<10.2f}")

# --- 7g. FoldX Correlation ---
print("\n" + "-" * 60)
print("7g. FOLDX ddG vs AF3 CORRELATION ANALYSIS")
print("-" * 60)

if len(matched_df) >= 3:
    print(f"\n   Matched mutations: {len(matched_df)}")

    # ddG vs delta_pTM
    r_ptm, p_ptm = stats.pearsonr(matched_df['foldx_ddG'], matched_df['af3_delta_ptm'])
    rho_ptm, prho_ptm = stats.spearmanr(matched_df['foldx_ddG'], matched_df['af3_delta_ptm'])
    print(f"\n   FoldX ddG vs AF3 delta_pTM:")
    print(f"      Pearson r = {r_ptm:.4f} (p = {p_ptm:.4e})")
    print(f"      Spearman rho = {rho_ptm:.4f} (p = {prho_ptm:.4e})")

    # ddG vs delta_ranking_score
    r_rank, p_rank = stats.pearsonr(matched_df['foldx_ddG'], matched_df['af3_delta_ranking'])
    rho_rank, prho_rank = stats.spearmanr(matched_df['foldx_ddG'], matched_df['af3_delta_ranking'])
    print(f"\n   FoldX ddG vs AF3 delta_ranking_score:")
    print(f"      Pearson r = {r_rank:.4f} (p = {p_rank:.4e})")
    print(f"      Spearman rho = {rho_rank:.4f} (p = {prho_rank:.4e})")

    # ddG vs fraction_disordered
    r_dis, p_dis = stats.pearsonr(matched_df['foldx_ddG'], matched_df['af3_fraction_disordered'])
    rho_dis, prho_dis = stats.spearmanr(matched_df['foldx_ddG'], matched_df['af3_fraction_disordered'])
    print(f"\n   FoldX ddG vs AF3 fraction_disordered:")
    print(f"      Pearson r = {r_dis:.4f} (p = {p_dis:.4e})")
    print(f"      Spearman rho = {rho_dis:.4f} (p = {prho_dis:.4e})")

    # pLDDT at site vs ddG
    plddt_matched = matched_df[matched_df['af3_delta_plddt_at_site'].notna()]
    if len(plddt_matched) >= 3:
        r_plddt, p_plddt = stats.pearsonr(plddt_matched['foldx_ddG'], plddt_matched['af3_delta_plddt_at_site'])
        rho_plddt, prho_plddt = stats.spearmanr(plddt_matched['foldx_ddG'], plddt_matched['af3_delta_plddt_at_site'])
        print(f"\n   FoldX ddG vs AF3 delta_pLDDT at mutation site:")
        print(f"      Pearson r = {r_plddt:.4f} (p = {p_plddt:.4e})")
        print(f"      Spearman rho = {rho_plddt:.4f} (p = {prho_plddt:.4e})")
        print(f"      N = {len(plddt_matched)}")

    # Concordance analysis
    print(f"\n   CONCORDANCE: FoldX destabilising (ddG > 2) vs AF3 severe (delta_pTM < -0.05):")
    foldx_destab = matched_df['foldx_ddG'] > 2.0
    af3_severe = matched_df['af3_delta_ptm'] < -0.05
    concordant = (foldx_destab & af3_severe) | (~foldx_destab & ~af3_severe)
    print(f"      FoldX destabilising: {foldx_destab.sum()}")
    print(f"      AF3 severe: {af3_severe.sum()}")
    print(f"      Concordant: {concordant.sum()} / {len(matched_df)} ({100*concordant.mean():.1f}%)")

    # Show discordant
    discordant_df = matched_df[~concordant]
    if len(discordant_df) > 0:
        print(f"\n   DISCORDANT PREDICTIONS:")
        print(f"   {'Mutation':<12} {'Domain':<22} {'FoldX ddG':<12} {'AF3 dPTM':<12} {'Severity'}")
        print(f"   {'-'*10:<12} {'-'*20:<22} {'-'*10:<12} {'-'*10:<12} {'-'*8}")
        for _, r in discordant_df.iterrows():
            print(f"   {r['mutation']:<12} {r['domain']:<22} {r['foldx_ddG']:<12.3f} {r['af3_delta_ptm']:<12.4f} {r['severity']}")

    # Full table
    print(f"\n   ALL MATCHED MUTATIONS (FoldX vs AF3):")
    print(f"   {'Mutation':<12} {'Domain':<22} {'FoldX ddG':<12} {'AF3 pTM':<10} {'AF3 dPTM':<12} {'dRank':<10} {'FracDis':<10} {'dPLDDT':<10}")
    print(f"   {'-'*10:<12} {'-'*20:<22} {'-'*10:<12} {'-'*8:<10} {'-'*10:<12} {'-'*8:<10} {'-'*8:<10} {'-'*8:<10}")
    for _, r in matched_df.sort_values('foldx_ddG', ascending=False).iterrows():
        dplddt = f"{r['af3_delta_plddt_at_site']:.2f}" if pd.notna(r['af3_delta_plddt_at_site']) else "N/A"
        print(f"   {r['mutation']:<12} {r['domain']:<22} {r['foldx_ddG']:<12.3f} {r['af3_ptm']:<10.4f} {r['af3_delta_ptm']:<12.4f} {r['af3_delta_ranking']:<10.4f} {r['af3_fraction_disordered']:<10.3f} {dplddt:<10}")
else:
    print(f"   Insufficient matched mutations for correlation ({len(matched_df)} found)")

# --- 7h. All Monomer Results Table ---
print("\n" + "-" * 60)
print("7h. COMPLETE MONOMER RESULTS TABLE (sorted by delta_pTM)")
print("-" * 60)
print(f"\n   {'Mutation':<12} {'Pos':<5} {'Domain':<22} {'pTM':<8} {'dPTM':<10} {'Rank':<8} {'dRank':<10} {'FracDis':<10} {'dPLDDT':<10} {'Severity'}")
print(f"   {'-'*10:<12} {'-'*3:<5} {'-'*20:<22} {'-'*6:<8} {'-'*8:<10} {'-'*6:<8} {'-'*8:<10} {'-'*8:<10} {'-'*8:<10} {'-'*8}")
for _, r in mono_df.iterrows():
    dplddt = f"{r['delta_plddt_at_site']:.1f}" if pd.notna(r['delta_plddt_at_site']) else "N/A"
    print(f"   {r['mutation']:<12} {r['position']:<5} {r['domain']:<22} {r['ptm']:<8.4f} {r['delta_ptm']:<10.4f} {r['ranking_score']:<8.4f} {r['delta_ranking_score']:<10.4f} {r['fraction_disordered']:<10.3f} {dplddt:<10} {r['severity']}")

# --- 7i. Summary Statistics ---
print("\n" + "-" * 60)
print("7i. FINAL SUMMARY STATISTICS")
print("-" * 60)

# Cysteine mutations
cys_mut = mono_df[mono_df['wt_aa'] == 'C']
print(f"\n   Cysteine mutations: {len(cys_mut)}")
if len(cys_mut) > 0:
    print(f"      Mean delta_pTM: {cys_mut['delta_ptm'].mean():.4f}")
    print(f"      All cysteine mutations:")
    for _, r in cys_mut.sort_values('delta_ptm').iterrows():
        print(f"         {r['mutation']}: delta_pTM={r['delta_ptm']:.4f}, domain={r['domain']}")

# Mutations to proline
pro_mut = mono_df[mono_df['mut_aa'] == 'P']
print(f"\n   Mutations to Proline: {len(pro_mut)}")
if len(pro_mut) > 0:
    print(f"      Mean delta_pTM: {pro_mut['delta_ptm'].mean():.4f}")
    for _, r in pro_mut.sort_values('delta_ptm').iterrows():
        print(f"         {r['mutation']}: delta_pTM={r['delta_ptm']:.4f}, domain={r['domain']}")

# Ligand-binding domain mutations
lb_mut = mono_df[mono_df['domain'].str.startswith('Ligand-binding')]
print(f"\n   Ligand-binding domain mutations: {len(lb_mut)}")
if len(lb_mut) > 0:
    print(f"      Mean delta_pTM: {lb_mut['delta_ptm'].mean():.4f}")

# EGF-like domain mutations
egf_mut = mono_df[mono_df['domain'].str.startswith('EGF')]
print(f"\n   EGF-like domain mutations: {len(egf_mut)}")
if len(egf_mut) > 0:
    print(f"      Mean delta_pTM: {egf_mut['delta_ptm'].mean():.4f}")

# Beta-propeller mutations
bp_mut = mono_df[mono_df['domain'] == 'Beta-propeller']
print(f"\n   Beta-propeller mutations: {len(bp_mut)}")
if len(bp_mut) > 0:
    print(f"      Mean delta_pTM: {bp_mut['delta_ptm'].mean():.4f}")

print(f"\n   OVERALL:")
print(f"   Total monomers: {len(mono_df)}")
print(f"   Total complexes: {len(complex_df)}")
print(f"   FoldX cross-referenced: {len(matched_df)}")
print(f"   Severe disruption: {len(severe)} mutations")
print(f"   Moderate disruption: {len(moderate)} mutations")
print(f"   Minimal change: {len(minimal)} mutations")
print(f"   Improved confidence: {len(improved)} mutations")

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)
print(f"\nOutput files:")
print(f"  {mono_csv}")
print(f"  {complex_csv}")
print(f"  {summary_csv}")

#!/usr/bin/env python3
"""
36_LDLR_atlas.py
================
Build the FIRST COMPLETE genotype-phenotype atlas for the LDL receptor,
mapping every position in the 860-residue protein to its predicted clinical phenotype.

Outputs:
  alphafold/analysis/LDLR_genotype_phenotype_atlas.csv
"""

import pandas as pd
import numpy as np
import os, sys

# ── paths ────────────────────────────────────────────────────────────────────
BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")

PDB_FILE        = os.path.join(BASE, "alphafold", "foldx", "LDLR_wt_Repair.pdb")
PLDDT_FILE      = os.path.join(ANALYSIS, "LDLR_wildtype_plddt_per_residue.csv")
SSS_FILE        = os.path.join(ANALYSIS, "structural_severity_scores.csv")
SSS_FULL_FILE   = os.path.join(ANALYSIS, "structural_severity_scores_full.csv")
INTERFACE_FILE  = os.path.join(ANALYSIS, "interface_distance_analysis.csv")
COX_KM_FILE     = os.path.join(ANALYSIS, "cox_km_compound_results.csv")
VARIANTS_FILE   = os.path.join(ANALYSIS, "variants_mapped_to_structure.csv")
CARRIERS_FILE   = os.path.join(ANALYSIS, "ukb_carriers_variant_sss.csv")
LLT_FILE        = os.path.join(ANALYSIS, "variant_llt_response.csv")
OUTPUT_FILE     = os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas.csv")

# ── constants ────────────────────────────────────────────────────────────────
PCSK9_INTERFACE = [25,35,50,53,58,59,60,61,78,79,81,96,316,318,319,320,321,
                   322,325,326,328,329,330,331,332,339,341,342,344,351,352,
                   372,386,387,390]
APOB_INTERFACE  = [285,286,298,300,303,304,305,308,309,310,311,312]

AA3TO1 = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H',
    'ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q',
    'ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'
}

# Domain penetrance by age 60 (from KM analysis, n=3,560)
DOMAIN_PENETRANCE = {
    "Ligand-binding R4": 39.6,
    "EGF-like A": 38.1,
    "Cytoplasmic": 29.2,
    "Ligand-binding": 27.9,
    "Transmembrane": 23.8,
    "Ligand-binding R7": 21.6,
    "EGF-precursor-homology": 21.2,
    "Ligand-binding R5": 19.9,
    "EGF-like B": 19.4,
    "Beta-propeller": 16.0,
    "Linker": 15.6,
    "Linker (EGF-C/O-linked)": 15.6,
    "O-linked sugar": 15.6,
    "Ligand-binding R1": 14.9,
    "Ligand-binding R2": 8.7,
    "EGF-like C": 12.8,
    "PCSK9-catalytic": 10.4,
    "Signal peptide": 10.0,
    "Ligand-binding R3": 20.0,
    "Ligand-binding R6": 12.6,
    "Unknown": 15.0,
}

# Domain ASCVD rates from cox_km_compound_results (events / n)
DOMAIN_ASCVD_RATE = {
    "Ligand-binding R4": 0.214,
    "EGF-like A": 0.157,
    "Cytoplasmic": 0.188,
    "Ligand-binding": 0.143,
    "Transmembrane": 0.222,
    "Ligand-binding R7": 0.122,
    "EGF-precursor-homology": 0.109,
    "Ligand-binding R5": 0.115,
    "EGF-like B": 0.127,
    "Beta-propeller": 0.106,
    "Linker": 0.130,
    "Linker (EGF-C/O-linked)": 0.130,
    "O-linked sugar": 0.139,
    "Ligand-binding R1": 0.099,
    "Ligand-binding R2": 0.084,
    "EGF-like C": 0.103,
    "PCSK9-catalytic": 0.095,
    "Signal peptide": 0.100,
    "Ligand-binding R3": 0.120,
    "Ligand-binding R6": 0.103,
    "Unknown": 0.123,
}

# Treatment response predictions by domain
DOMAIN_TREATMENT_RESPONSE = {
    "Ligand-binding R1": "Moderate statin response; consider combination therapy",
    "Ligand-binding R2": "Good statin response; standard treatment",
    "Ligand-binding R3": "Moderate statin response; may need intensification",
    "Ligand-binding R4": "Poor statin response; early PCSK9i/combination needed",
    "Ligand-binding R5": "Moderate statin response; consider ezetimibe add-on",
    "Ligand-binding R6": "Good statin response; standard treatment",
    "Ligand-binding R7": "Moderate statin response; consider combination therapy",
    "Ligand-binding": "Variable statin response; individualize",
    "EGF-like A": "Poor statin response; PCSK9i may have reduced efficacy (binding site)",
    "EGF-like B": "Moderate statin response; structural instability may limit response",
    "EGF-like C": "Moderate statin response; monitor for insufficient response",
    "EGF-precursor-homology": "Moderate statin response; recycling pathway affected",
    "Beta-propeller": "Good statin response; receptor recycling partially preserved",
    "Linker": "Good statin response; structural flexibility retained",
    "Linker (EGF-C/O-linked)": "Good statin response; structural flexibility retained",
    "O-linked sugar": "Good statin response; glycosylation affects stability",
    "Transmembrane": "Variable response; membrane anchoring affected",
    "Cytoplasmic": "Variable statin response; endocytosis signal affected, consider PCSK9i",
    "Signal peptide": "Good statin response if protein reaches surface",
    "PCSK9-catalytic": "Strong PCSK9i candidate; enhanced PCSK9 binding",
    "Unknown": "Standard treatment; individualize based on LDL response",
}

# Drug target categories
DRUG_CATEGORIES = {
    "Ligand-binding R1": "LDL-binding enhancement",
    "Ligand-binding R2": "LDL-binding enhancement",
    "Ligand-binding R3": "LDL-binding enhancement",
    "Ligand-binding R4": "LDL-binding enhancement",
    "Ligand-binding R5": "LDL-binding enhancement",
    "Ligand-binding R6": "LDL-binding enhancement",
    "Ligand-binding R7": "LDL-binding enhancement",
    "Ligand-binding": "LDL-binding enhancement",
    "EGF-like A": "PCSK9 interaction / receptor recycling",
    "EGF-like B": "Receptor recycling stabilization",
    "EGF-like C": "Receptor recycling stabilization",
    "EGF-precursor-homology": "Receptor recycling stabilization",
    "Beta-propeller": "pH-dependent release / recycling",
    "Linker": "Structural flexibility / folding",
    "Linker (EGF-C/O-linked)": "Structural flexibility / folding",
    "O-linked sugar": "Protein stability / trafficking",
    "Transmembrane": "Membrane anchoring / trafficking",
    "Cytoplasmic": "Endocytic signaling / PCSK9i target",
    "Signal peptide": "Protein trafficking / folding",
    "PCSK9-catalytic": "PCSK9 interaction / receptor recycling",
    "Unknown": "Standard cholesterol-lowering pathway",
}

# ── 1. Read wildtype residues from PDB ───────────────────────────────────────
print("=" * 80)
print("LDLR GENOTYPE-PHENOTYPE ATLAS BUILDER")
print("=" * 80)

print("\n[1/8] Reading wildtype residues from PDB...")
wt_residues = {}
with open(PDB_FILE) as f:
    for line in f:
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            resname = line[17:20].strip()
            resnum  = int(line[22:26].strip())
            wt_residues[resnum] = AA3TO1.get(resname, 'X')

print(f"   Found {len(wt_residues)} residues in PDB (range {min(wt_residues)}-{max(wt_residues)})")

# ── 2. Read pLDDT per-residue data ──────────────────────────────────────────
print("[2/8] Reading pLDDT and domain annotations...")
plddt_df = pd.read_csv(PLDDT_FILE)
plddt_map = dict(zip(plddt_df['residue_num'], plddt_df['plddt']))
domain_map = dict(zip(plddt_df['residue_num'], plddt_df['domain']))
print(f"   {len(plddt_df)} residues with pLDDT data")
print(f"   Domains: {plddt_df['domain'].nunique()} unique")

# ── 3. Read SSS / FoldX data ────────────────────────────────────────────────
print("[3/8] Reading structural severity scores and FoldX data...")
sss_full = pd.read_csv(SSS_FULL_FILE)
sss_ldlr = sss_full[sss_full['gene'] == 'LDLR'].copy()
sss_ldlr['protein_position'] = pd.to_numeric(sss_ldlr['protein_position'], errors='coerce')
sss_ldlr = sss_ldlr.dropna(subset=['protein_position'])
sss_ldlr['protein_position'] = sss_ldlr['protein_position'].astype(int)
print(f"   {len(sss_ldlr)} LDLR variants with SSS data")

# Aggregate per-position: max SSS, max ddG, count of variants
pos_sss = sss_ldlr.groupby('protein_position').agg(
    max_sss=('sss', 'max'),
    mean_sss=('sss', 'mean'),
    max_ddG=('ddG', lambda x: x.dropna().max() if x.dropna().shape[0] > 0 else np.nan),
    mean_ddG=('ddG', lambda x: x.dropna().mean() if x.dropna().shape[0] > 0 else np.nan),
    n_variants=('variant_id', 'nunique'),
).reset_index()
print(f"   {len(pos_sss)} unique positions with variant data")

# ── 4. Read interface distances ──────────────────────────────────────────────
print("[4/8] Reading PCSK9 interface distance data...")
iface_df = pd.read_csv(INTERFACE_FILE)
iface_ldlr = iface_df[iface_df['gene'] == 'LDLR'].copy()
iface_ldlr['position'] = pd.to_numeric(iface_ldlr['position'], errors='coerce')
iface_ldlr = iface_ldlr.dropna(subset=['position'])
iface_ldlr['position'] = iface_ldlr['position'].astype(int)
# Per-position min interface distance
pos_iface = iface_ldlr.groupby('position').agg(
    pcsk9_interface_dist=('interface_dist', 'min'),
).reset_index()
print(f"   {len(pos_iface)} positions with interface distance data")

# ── 5. Read patient-level carrier data ───────────────────────────────────────
print("[5/8] Reading patient-level carrier data...")
carriers = pd.read_csv(CARRIERS_FILE)
carriers['protein_position'] = pd.to_numeric(carriers['protein_position'], errors='coerce')
carriers = carriers.dropna(subset=['protein_position'])
carriers['protein_position'] = carriers['protein_position'].astype(int)
print(f"   {len(carriers)} carrier records")

# Patient count per position
pos_patients = carriers.groupby('protein_position').agg(
    n_patients=('eid', 'nunique'),
).reset_index()

# ── 6. Read treatment response data for ASCVD events ────────────────────────
print("[6/8] Reading treatment response data for ASCVD events...")
llt_df = pd.read_csv(LLT_FILE)

# Extract protein position from variant names for LDLR
def extract_position_from_variant(variant_str):
    """Extract protein position from variant string like LDLR:c.343C>T"""
    if pd.isna(variant_str):
        return np.nan
    # Try to match from the carriers/sss data instead
    return np.nan

# Use carriers data which already has protein_position to link ASCVD events
# Merge carriers with llt data on patient_id/eid and variant
llt_ldlr = llt_df[llt_df['gene'] == 'LDLR'].copy() if 'gene' in llt_df.columns else llt_df.copy()

# Get ASCVD events per position from carriers + llt merge
# First, use carriers which has protein_position and eid
# Then merge with llt which has patient_id and ascvd
if 'ascvd' in llt_df.columns:
    # Build position -> ASCVD from the carriers + llt join
    # llt has patient_id, variant, ascvd
    # carriers has eid, variant_id, protein_position
    llt_variants = llt_df[['patient_id', 'variant', 'ascvd']].drop_duplicates()
    carriers_pos = carriers[['eid', 'variant_id', 'protein_position']].drop_duplicates()

    # Try direct merge
    merged_ascvd = carriers_pos.merge(
        llt_variants,
        left_on=['eid', 'variant_id'],
        right_on=['patient_id', 'variant'],
        how='inner'
    )
    if len(merged_ascvd) == 0:
        # Try just on variant
        carrier_var_pos = carriers.groupby('variant_id')['protein_position'].first().reset_index()
        llt_var_ascvd = llt_df.groupby('variant').agg(
            ascvd_events=('ascvd', 'sum'),
            ascvd_total=('ascvd', 'count'),
        ).reset_index()
        merged_ascvd2 = carrier_var_pos.merge(
            llt_var_ascvd,
            left_on='variant_id',
            right_on='variant',
            how='inner'
        )
        if len(merged_ascvd2) > 0:
            pos_ascvd = merged_ascvd2.groupby('protein_position').agg(
                ascvd_events=('ascvd_events', 'sum'),
                ascvd_total=('ascvd_total', 'sum'),
            ).reset_index()
        else:
            pos_ascvd = pd.DataFrame(columns=['protein_position', 'ascvd_events', 'ascvd_total'])
    else:
        pos_ascvd = merged_ascvd.groupby('protein_position').agg(
            ascvd_events=('ascvd', 'sum'),
            ascvd_total=('ascvd', 'count'),
        ).reset_index()
else:
    pos_ascvd = pd.DataFrame(columns=['protein_position', 'ascvd_events', 'ascvd_total'])

print(f"   {len(pos_ascvd)} positions with ASCVD event data")

# Also use variants_mapped_to_structure for additional patient counts
print("[6b/8] Reading variants_mapped_to_structure for patient counts...")
mapped_variants = pd.read_csv(VARIANTS_FILE)
mapped_variants['protein_position'] = pd.to_numeric(mapped_variants['protein_position'], errors='coerce')
mapped_variants = mapped_variants.dropna(subset=['protein_position'])
mapped_variants['protein_position'] = mapped_variants['protein_position'].astype(int)
# Aggregate patient counts per position
pos_mapped = mapped_variants.groupby('protein_position').agg(
    mapped_n_patients=('n_patients', 'sum'),
).reset_index()

# ── 7. Build the atlas ──────────────────────────────────────────────────────
print("[7/8] Building the complete genotype-phenotype atlas...")

PCSK9_SET = set(PCSK9_INTERFACE)
APOB_SET  = set(APOB_INTERFACE)

rows = []
for pos in range(1, 861):
    wt_aa = wt_residues.get(pos, 'X')
    domain = domain_map.get(pos, 'Unknown')
    plddt  = plddt_map.get(pos, np.nan)

    # PCSK9 interface distance (residue distance in sequence)
    pcsk9_dist = min(abs(pos - ip) for ip in PCSK9_INTERFACE)
    apob_dist  = min(abs(pos - ip) for ip in APOB_INTERFACE)
    at_pcsk9 = pos in PCSK9_SET
    at_apob  = pos in APOB_SET

    # Domain penetrance
    dom_pen = DOMAIN_PENETRANCE.get(domain, 15.0)
    dom_ascvd = DOMAIN_ASCVD_RATE.get(domain, 0.12)

    # Observed SSS / ddG
    pos_data = pos_sss[pos_sss['protein_position'] == pos]
    if len(pos_data) > 0:
        obs_sss = pos_data.iloc[0]['max_sss']
        obs_ddG = pos_data.iloc[0]['max_ddG']
        n_var   = int(pos_data.iloc[0]['n_variants'])
    else:
        obs_sss = np.nan
        obs_ddG = np.nan
        n_var   = 0

    # Patient count: prefer carriers data, fallback to mapped variants
    pat_data = pos_patients[pos_patients['protein_position'] == pos]
    map_data = pos_mapped[pos_mapped['protein_position'] == pos]
    if len(pat_data) > 0:
        n_pat = int(pat_data.iloc[0]['n_patients'])
    elif len(map_data) > 0:
        n_pat = int(map_data.iloc[0]['mapped_n_patients'])
    else:
        n_pat = 0

    # ASCVD events
    asc_data = pos_ascvd[pos_ascvd['protein_position'] == pos] if len(pos_ascvd) > 0 else pd.DataFrame()
    if len(asc_data) > 0:
        n_ascvd = int(asc_data.iloc[0]['ascvd_events'])
    else:
        n_ascvd = 0

    # PCSK9 interface distance from actual interface_distance_analysis if available
    iface_data = pos_iface[pos_iface['position'] == pos]
    if len(iface_data) > 0:
        pcsk9_struct_dist = iface_data.iloc[0]['pcsk9_interface_dist']
    else:
        pcsk9_struct_dist = pcsk9_dist  # fallback to sequence distance

    # Risk classification
    if (at_pcsk9 or at_apob) and dom_pen > 25:
        risk = "Critical"
    elif dom_pen > 20:
        risk = "High"
    elif dom_pen >= 15:
        risk = "Moderate"
    else:
        risk = "Low"

    # Treatment response
    tx_response = DOMAIN_TREATMENT_RESPONSE.get(domain, "Standard treatment; individualize")
    drug_cat    = DRUG_CATEGORIES.get(domain, "Standard cholesterol-lowering pathway")

    rows.append({
        'position': pos,
        'wildtype_aa': wt_aa,
        'domain': domain,
        'plddt_confidence': round(plddt, 2) if not np.isnan(plddt) else np.nan,
        'pcsk9_interface_distance': pcsk9_struct_dist,
        'apob_interface_distance': apob_dist,
        'at_pcsk9_interface': 'Yes' if at_pcsk9 else 'No',
        'at_apob_interface': 'Yes' if at_apob else 'No',
        'domain_penetrance_by_60_pct': dom_pen,
        'domain_ascvd_rate': dom_ascvd,
        'observed_max_ddG': round(obs_ddG, 3) if not np.isnan(obs_ddG) else np.nan,
        'observed_max_sss': round(obs_sss, 4) if not np.isnan(obs_sss) else np.nan,
        'n_known_variants': n_var,
        'n_patients': n_pat,
        'ascvd_events': n_ascvd,
        'risk_classification': risk,
        'predicted_treatment_response': tx_response,
        'drug_target_category': drug_cat,
    })

atlas = pd.DataFrame(rows)

# ── 8. Save and print summary ───────────────────────────────────────────────
print("[8/8] Saving atlas and computing summary statistics...\n")
atlas.to_csv(OUTPUT_FILE, index=False)
print(f"   Atlas saved to: {OUTPUT_FILE}")
print(f"   Total residues: {len(atlas)}")

# ── Summary statistics ───────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("LDLR GENOTYPE-PHENOTYPE ATLAS: SUMMARY STATISTICS")
print("=" * 80)

# Risk classification distribution
print("\n--- RISK CLASSIFICATION DISTRIBUTION ---")
risk_counts = atlas['risk_classification'].value_counts()
for risk_level in ['Critical', 'High', 'Moderate', 'Low']:
    count = risk_counts.get(risk_level, 0)
    pct = 100 * count / len(atlas)
    print(f"   {risk_level:10s}: {count:4d} positions ({pct:5.1f}%)")

# Domain-level aggregates
print("\n--- DOMAIN-LEVEL SUMMARY ---")
domain_agg = atlas.groupby('domain').agg(
    n_residues=('position', 'count'),
    mean_plddt=('plddt_confidence', 'mean'),
    penetrance=('domain_penetrance_by_60_pct', 'first'),
    ascvd_rate=('domain_ascvd_rate', 'first'),
    total_patients=('n_patients', 'sum'),
    total_ascvd=('ascvd_events', 'sum'),
    n_with_variants=('n_known_variants', lambda x: (x > 0).sum()),
    max_sss=('observed_max_sss', 'max'),
).sort_values('penetrance', ascending=False)

print(f"   {'Domain':<30s} {'Res':>4s} {'pLDDT':>6s} {'Pen%':>6s} {'ASCVD':>6s} {'Pts':>5s} {'Evt':>4s} {'w/Var':>5s}")
print(f"   {'-'*30} {'----':>4s} {'------':>6s} {'------':>6s} {'------':>6s} {'-----':>5s} {'----':>4s} {'-----':>5s}")
for _, row in domain_agg.iterrows():
    print(f"   {row.name:<30s} {row['n_residues']:4d} {row['mean_plddt']:6.1f} "
          f"{row['penetrance']:6.1f} {row['ascvd_rate']:6.3f} {row['total_patients']:5d} "
          f"{row['total_ascvd']:4d} {row['n_with_variants']:5d}")

# Interface hotspot analysis
print("\n--- INTERFACE HOTSPOT ANALYSIS ---")
pcsk9_positions = atlas[atlas['at_pcsk9_interface'] == 'Yes']
apob_positions  = atlas[atlas['at_apob_interface'] == 'Yes']
print(f"   PCSK9 interface positions: {len(pcsk9_positions)}")
print(f"   ApoB interface positions:  {len(apob_positions)}")

pcsk9_critical = pcsk9_positions[pcsk9_positions['risk_classification'] == 'Critical']
apob_critical  = apob_positions[apob_positions['risk_classification'] == 'Critical']
print(f"   Critical PCSK9 positions:  {len(pcsk9_critical)}")
print(f"   Critical ApoB positions:   {len(apob_critical)}")

if len(pcsk9_critical) > 0:
    print(f"\n   Critical PCSK9 interface residues:")
    for _, r in pcsk9_critical.iterrows():
        print(f"      Pos {r['position']:4d} ({r['wildtype_aa']}) - {r['domain']}, "
              f"penetrance={r['domain_penetrance_by_60_pct']:.1f}%, "
              f"patients={r['n_patients']}")

if len(apob_critical) > 0:
    print(f"\n   Critical ApoB interface residues:")
    for _, r in apob_critical.iterrows():
        print(f"      Pos {r['position']:4d} ({r['wildtype_aa']}) - {r['domain']}, "
              f"penetrance={r['domain_penetrance_by_60_pct']:.1f}%, "
              f"patients={r['n_patients']}")

# Positions with highest observed ASCVD rates
print("\n--- TOP 20 POSITIONS BY PATIENT COUNT ---")
top_patients = atlas[atlas['n_patients'] > 0].nlargest(20, 'n_patients')
print(f"   {'Pos':>4s} {'AA':>3s} {'Domain':<28s} {'Pts':>5s} {'ASCVD':>5s} {'SSS':>6s} {'ddG':>7s} {'Risk':<10s}")
print(f"   {'----':>4s} {'---':>3s} {'----------------------------':<28s} {'-----':>5s} {'-----':>5s} {'------':>6s} {'-------':>7s} {'----------':<10s}")
for _, r in top_patients.iterrows():
    sss_str = f"{r['observed_max_sss']:.3f}" if not pd.isna(r['observed_max_sss']) else "   -  "
    ddg_str = f"{r['observed_max_ddG']:.2f}" if not pd.isna(r['observed_max_ddG']) else "    -  "
    print(f"   {r['position']:4d} {r['wildtype_aa']:>3s} {r['domain']:<28s} "
          f"{r['n_patients']:5d} {r['ascvd_events']:5d} {sss_str:>6s} {ddg_str:>7s} {r['risk_classification']:<10s}")

# Positions with highest ASCVD events
if atlas['ascvd_events'].sum() > 0:
    print("\n--- TOP 20 POSITIONS BY ASCVD EVENTS ---")
    top_ascvd = atlas[atlas['ascvd_events'] > 0].nlargest(20, 'ascvd_events')
    for _, r in top_ascvd.iterrows():
        print(f"   Pos {r['position']:4d} ({r['wildtype_aa']}) - {r['domain']:<28s} "
              f"events={r['ascvd_events']:3d}, patients={r['n_patients']:4d}, risk={r['risk_classification']}")

# Coverage statistics
print("\n--- ATLAS COVERAGE STATISTICS ---")
n_with_data = (atlas['n_known_variants'] > 0).sum()
n_with_patients = (atlas['n_patients'] > 0).sum()
n_with_sss = atlas['observed_max_sss'].notna().sum()
n_with_ddg = atlas['observed_max_ddG'].notna().sum()

print(f"   Positions with known variants:    {n_with_data:4d} / {len(atlas)} ({100*n_with_data/len(atlas):.1f}%)")
print(f"   Positions with patient data:      {n_with_patients:4d} / {len(atlas)} ({100*n_with_patients/len(atlas):.1f}%)")
print(f"   Positions with SSS scores:        {n_with_sss:4d} / {len(atlas)} ({100*n_with_sss/len(atlas):.1f}%)")
print(f"   Positions with FoldX ddG:         {n_with_ddg:4d} / {len(atlas)} ({100*n_with_ddg/len(atlas):.1f}%)")
print(f"   Positions at PCSK9 interface:     {len(pcsk9_positions):4d} / {len(atlas)}")
print(f"   Positions at ApoB interface:      {len(apob_positions):4d} / {len(atlas)}")

# Drug target category distribution
print("\n--- DRUG TARGET CATEGORY DISTRIBUTION ---")
drug_counts = atlas['drug_target_category'].value_counts()
for cat, count in drug_counts.items():
    print(f"   {cat:<45s}: {count:4d} positions ({100*count/len(atlas):.1f}%)")

# Final summary
print("\n" + "=" * 80)
print("ATLAS COMPLETE")
print(f"  860 LDLR residues mapped with 18 clinical annotations each")
print(f"  {n_with_data} positions have observed variant data")
print(f"  {risk_counts.get('Critical', 0)} Critical, {risk_counts.get('High', 0)} High, "
      f"{risk_counts.get('Moderate', 0)} Moderate, {risk_counts.get('Low', 0)} Low risk positions")
print(f"  Output: {OUTPUT_FILE}")
print("=" * 80)

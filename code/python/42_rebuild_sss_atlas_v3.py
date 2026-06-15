#!/usr/bin/env python3
"""
42_rebuild_sss_atlas_v3.py
==========================
Rebuild the LDLR genotype-phenotype atlas v3 using COMPLETE FoldX saturation
mutagenesis data (16,340 mutations). Updates SSS scores, risk classifications,
and generates publication-ready outputs.

Key improvements over v2:
  - Domain-level imputation for 278 FoldX-failed positions
  - SSS v3 recalculated with position-specific mutational sensitivity
  - Complete 860-position risk classification
  - Publication tables and statistics

Outputs:
  alphafold/analysis/LDLR_genotype_phenotype_atlas_v3.csv  (860 positions)
  alphafold/analysis/sss_v3_catalogue.csv                  (all variants rescored)
  alphafold/analysis/atlas_v3_summary_statistics.csv        (publication table)
  alphafold/analysis/saturation_heatmap_data.csv            (for figure)
"""

import pandas as pd
import numpy as np
import os
import sys
from scipy import stats

# ── paths ────────────────────────────────────────────────────────────────────
BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")

# Input files
SATURATION_CSV = os.path.join(ANALYSIS, "foldx_saturation_mutagenesis.csv")
ATLAS_V2 = os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas_v2.csv")
SSS_V2_CAT = os.path.join(ANALYSIS, "sss_v2_catalogue.csv")
PLDDT_FILE = os.path.join(ANALYSIS, "LDLR_wildtype_plddt_per_residue.csv")
INTERFACE_FILE = os.path.join(ANALYSIS, "interface_distance_analysis.csv")

# Output files
ATLAS_V3 = os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas_v3.csv")
SSS_V3_CAT = os.path.join(ANALYSIS, "sss_v3_catalogue.csv")
SUMMARY_STATS = os.path.join(ANALYSIS, "atlas_v3_summary_statistics.csv")
HEATMAP_DATA = os.path.join(ANALYSIS, "saturation_heatmap_data.csv")

# ── Domain definitions (residue ranges) ──────────────────────────────────────
DOMAIN_RANGES = {
    "Signal peptide": (1, 21),
    "Ligand-binding R1": (22, 63),
    "Ligand-binding R2": (64, 105),
    "Ligand-binding R3": (106, 147),
    "Ligand-binding R4": (148, 189),
    "Ligand-binding R5": (190, 231),
    "Ligand-binding R6": (232, 273),
    "Ligand-binding R7": (274, 315),
    "EGF-like A": (316, 353),
    "EGF-like B": (354, 393),
    "Beta-propeller": (394, 632),
    "EGF-like C": (633, 672),
    "Linker (EGF-C/O-linked)": (673, 692),
    "O-linked sugar": (693, 749),
    "Transmembrane": (750, 789),
    "Cytoplasmic": (790, 860),
}

def get_domain(pos):
    for domain, (start, end) in DOMAIN_RANGES.items():
        if start <= pos <= end:
            return domain
    return "Unknown"

# Domain penetrance by age 60 (from KM analysis)
DOMAIN_PENETRANCE = {
    "Signal peptide": 10.0,
    "Ligand-binding R1": 14.9,
    "Ligand-binding R2": 8.7,
    "Ligand-binding R3": 20.0,
    "Ligand-binding R4": 39.6,
    "Ligand-binding R5": 19.9,
    "Ligand-binding R6": 12.6,
    "Ligand-binding R7": 21.6,
    "EGF-like A": 38.1,
    "EGF-like B": 19.4,
    "Beta-propeller": 16.0,
    "EGF-like C": 12.8,
    "Linker (EGF-C/O-linked)": 15.6,
    "O-linked sugar": 15.6,
    "Transmembrane": 23.8,
    "Cytoplasmic": 29.2,
}

DOMAIN_ASCVD_RATE = {
    "Signal peptide": 0.100,
    "Ligand-binding R1": 0.099,
    "Ligand-binding R2": 0.084,
    "Ligand-binding R3": 0.120,
    "Ligand-binding R4": 0.214,
    "Ligand-binding R5": 0.115,
    "Ligand-binding R6": 0.103,
    "Ligand-binding R7": 0.122,
    "EGF-like A": 0.157,
    "EGF-like B": 0.127,
    "Beta-propeller": 0.106,
    "EGF-like C": 0.103,
    "Linker (EGF-C/O-linked)": 0.130,
    "O-linked sugar": 0.139,
    "Transmembrane": 0.222,
    "Cytoplasmic": 0.188,
}

DOMAIN_TREATMENT_RESPONSE = {
    "Signal peptide": "Good statin response if protein reaches surface",
    "Ligand-binding R1": "Moderate statin response; consider combination therapy",
    "Ligand-binding R2": "Good statin response; standard treatment",
    "Ligand-binding R3": "Moderate statin response; may need intensification",
    "Ligand-binding R4": "Poor statin response; early PCSK9i/combination needed",
    "Ligand-binding R5": "Moderate statin response; consider ezetimibe add-on",
    "Ligand-binding R6": "Good statin response; standard treatment",
    "Ligand-binding R7": "Moderate statin response; consider combination therapy",
    "EGF-like A": "Poor statin response; PCSK9i may have reduced efficacy",
    "EGF-like B": "Moderate statin response; structural instability may limit response",
    "Beta-propeller": "Good statin response; receptor recycling partially preserved",
    "EGF-like C": "Moderate statin response; monitor for insufficient response",
    "Linker (EGF-C/O-linked)": "Good statin response; structural flexibility retained",
    "O-linked sugar": "Good statin response; glycosylation affects stability",
    "Transmembrane": "Variable response; membrane anchoring affected",
    "Cytoplasmic": "Variable statin response; endocytosis signal affected, consider PCSK9i",
}

DRUG_CATEGORIES = {
    "Signal peptide": "Protein folding / trafficking",
    "Ligand-binding R1": "LDL-binding enhancement",
    "Ligand-binding R2": "LDL-binding enhancement",
    "Ligand-binding R3": "LDL-binding enhancement",
    "Ligand-binding R4": "LDL-binding enhancement",
    "Ligand-binding R5": "LDL-binding enhancement",
    "Ligand-binding R6": "LDL-binding enhancement",
    "Ligand-binding R7": "LDL-binding enhancement",
    "EGF-like A": "PCSK9 interaction / receptor recycling",
    "EGF-like B": "Receptor recycling stabilization",
    "Beta-propeller": "Receptor recycling / pH-dependent release",
    "EGF-like C": "Receptor recycling stabilization",
    "Linker (EGF-C/O-linked)": "Structural flexibility",
    "O-linked sugar": "Glycosylation / protein stability",
    "Transmembrane": "Membrane anchoring",
    "Cytoplasmic": "Endocytosis / intracellular trafficking",
}

# PCSK9 and ApoB interface residues
PCSK9_INTERFACE = [25,35,50,53,58,59,60,61,78,79,81,96,316,318,319,320,321,
                   322,325,326,328,329,330,331,332,339,341,342,344,351,352,
                   372,386,387,390]
APOB_INTERFACE = [285,286,298,300,303,304,305,308,309,310,311,312]


def main():
    print("=" * 80)
    print("REBUILDING LDLR GENOTYPE-PHENOTYPE ATLAS v3")
    print("Complete saturation mutagenesis integration")
    print("=" * 80)

    # ── 1. Load saturation mutagenesis data ──────────────────────────────────
    print("\n[1/7] Loading saturation mutagenesis data...")
    sat = pd.read_csv(SATURATION_CSV)
    print(f"   Total mutations: {len(sat)}")
    print(f"   Valid ddG values: {sat['ddG'].notna().sum()}")
    print(f"   Failed: {sat['ddG'].isna().sum()}")

    # Per-position statistics
    pos_stats = sat.groupby('position').agg(
        n_computed=('ddG', lambda x: x.notna().sum()),
        max_ddG=('ddG', 'max'),
        mean_ddG=('ddG', 'mean'),
        median_ddG=('ddG', 'median'),
        min_ddG=('ddG', 'min'),
        std_ddG=('ddG', 'std'),
        n_destabilizing=('ddG', lambda x: (x > 2).sum()),
        n_mild=('ddG', lambda x: ((x > 0.5) & (x <= 2)).sum()),
        n_neutral=('ddG', lambda x: ((x >= -0.5) & (x <= 0.5)).sum()),
        n_stabilizing=('ddG', lambda x: (x < -0.5).sum()),
    ).reset_index()

    # Worst mutation per position
    valid_sat = sat.dropna(subset=['ddG'])
    if len(valid_sat) > 0:
        worst_idx = valid_sat.groupby('position')['ddG'].idxmax()
        worst_muts = valid_sat.loc[worst_idx, ['position', 'foldx_notation', 'mut_aa', 'ddG']].rename(
            columns={'foldx_notation': 'worst_mutation', 'mut_aa': 'worst_mut_aa', 'ddG': 'worst_ddG'}
        )
        pos_stats = pos_stats.merge(worst_muts, on='position', how='left')

    valid_positions = pos_stats[pos_stats['n_computed'] > 0]
    failed_positions = pos_stats[pos_stats['n_computed'] == 0]
    print(f"   Positions with valid data: {len(valid_positions)}")
    print(f"   Positions with all failures: {len(failed_positions)}")

    # ── 2. Domain-level imputation for failed positions ──────────────────────
    print("\n[2/7] Imputing failed positions using domain-level statistics...")

    # Calculate domain-level averages from successful positions
    pos_stats['domain'] = pos_stats['position'].apply(get_domain)
    domain_avg = valid_positions.groupby(
        valid_positions['position'].apply(get_domain)
    ).agg(
        domain_max_ddG=('max_ddG', 'median'),
        domain_mean_ddG=('mean_ddG', 'median'),
        domain_median_ddG=('median_ddG', 'median'),
        domain_pct_destab=('n_destabilizing', lambda x: x.median()),
        domain_n_positions=('position', 'count'),
    ).reset_index().rename(columns={'position': 'domain'})

    # Rename the groupby key properly
    domain_avg.columns = ['domain', 'domain_max_ddG', 'domain_mean_ddG',
                          'domain_median_ddG', 'domain_pct_destab', 'domain_n_positions']

    print("   Domain-level saturation statistics:")
    for _, row in domain_avg.iterrows():
        print(f"     {row['domain']:30s}  max_ddG={row['domain_max_ddG']:6.2f}  "
              f"mean_ddG={row['domain_mean_ddG']:6.2f}  n_pos={int(row['domain_n_positions'])}")

    # Impute failed positions
    for idx in pos_stats.index:
        if pos_stats.loc[idx, 'n_computed'] == 0:
            domain = pos_stats.loc[idx, 'domain']
            match = domain_avg[domain_avg['domain'] == domain]
            if len(match) > 0:
                pos_stats.loc[idx, 'max_ddG'] = match.iloc[0]['domain_max_ddG']
                pos_stats.loc[idx, 'mean_ddG'] = match.iloc[0]['domain_mean_ddG']
                pos_stats.loc[idx, 'median_ddG'] = match.iloc[0]['domain_median_ddG']
                pos_stats.loc[idx, 'imputed'] = True
            else:
                pos_stats.loc[idx, 'imputed'] = True
        else:
            pos_stats.loc[idx, 'imputed'] = False

    imputed_count = pos_stats['imputed'].sum()
    print(f"   Imputed {int(imputed_count)} positions from domain averages")

    # ── 3. Load structural data ──────────────────────────────────────────────
    print("\n[3/7] Loading structural data (pLDDT, interfaces)...")

    # pLDDT
    try:
        plddt = pd.read_csv(PLDDT_FILE)
        plddt_dict = dict(zip(plddt['residue_number'], plddt['plddt']))
        print(f"   pLDDT loaded for {len(plddt_dict)} residues")
    except Exception:
        plddt_dict = {}
        print("   WARNING: pLDDT file not found, using defaults")

    # Interface distances
    try:
        iface = pd.read_csv(INTERFACE_FILE)
        pcsk9_dist = dict(zip(iface['residue_number'], iface['pcsk9_min_distance']))
        apob_dist = dict(zip(iface['residue_number'], iface['apob_min_distance']))
        print(f"   Interface distances loaded for {len(pcsk9_dist)} residues")
    except Exception:
        pcsk9_dist = {}
        apob_dist = {}
        print("   WARNING: Interface file not found, using defaults")

    # ── 4. Load existing clinical data from atlas v2 ─────────────────────────
    print("\n[4/7] Loading clinical data from atlas v2...")
    atlas_v2 = pd.read_csv(ATLAS_V2)
    clinical_cols = ['position', 'observed_max_sss', 'n_known_variants',
                     'n_patients', 'ascvd_events']
    clinical = atlas_v2[clinical_cols].copy()
    print(f"   Clinical data for {len(clinical)} positions")

    # ── 5. Build atlas v3 ────────────────────────────────────────────────────
    print("\n[5/7] Building atlas v3 (860 positions)...")

    atlas = pd.DataFrame({'position': range(1, 861)})
    atlas['wildtype_aa'] = atlas['position'].map(
        dict(zip(sat['position'], sat['wt_aa']))
    )
    atlas['domain'] = atlas['position'].apply(get_domain)

    # Structural features
    atlas['plddt'] = atlas['position'].map(plddt_dict)
    atlas['plddt_confidence'] = atlas['plddt'].apply(
        lambda x: 'Very High' if pd.notna(x) and x >= 90
        else 'High' if pd.notna(x) and x >= 70
        else 'Low' if pd.notna(x) and x >= 50
        else 'Very Low' if pd.notna(x)
        else 'Unknown'
    )

    # Interface features
    atlas['pcsk9_interface_distance'] = atlas['position'].map(pcsk9_dist)
    atlas['apob_interface_distance'] = atlas['position'].map(apob_dist)
    atlas['at_pcsk9_interface'] = atlas['position'].isin(PCSK9_INTERFACE)
    atlas['at_apob_interface'] = atlas['position'].isin(APOB_INTERFACE)

    # Domain-level clinical data
    atlas['domain_penetrance_by_60_pct'] = atlas['domain'].map(DOMAIN_PENETRANCE)
    atlas['domain_ascvd_rate'] = atlas['domain'].map(DOMAIN_ASCVD_RATE)
    atlas['predicted_treatment_response'] = atlas['domain'].map(DOMAIN_TREATMENT_RESPONSE)
    atlas['drug_target_category'] = atlas['domain'].map(DRUG_CATEGORIES)

    # Saturation mutagenesis data
    sat_merge = pos_stats[['position', 'n_computed', 'max_ddG', 'mean_ddG',
                           'median_ddG', 'min_ddG', 'std_ddG',
                           'n_destabilizing', 'n_mild', 'n_neutral',
                           'n_stabilizing', 'imputed']].copy()
    sat_merge.columns = ['position', 'sat_n_computed', 'sat_max_ddG', 'sat_mean_ddG',
                         'sat_median_ddG', 'sat_min_ddG', 'sat_std_ddG',
                         'sat_n_destabilizing', 'sat_n_mild', 'sat_n_neutral',
                         'sat_n_stabilizing', 'sat_imputed']

    if 'worst_mutation' in pos_stats.columns:
        sat_merge = sat_merge.merge(
            pos_stats[['position', 'worst_mutation', 'worst_mut_aa', 'worst_ddG']],
            on='position', how='left'
        )
        sat_merge.rename(columns={
            'worst_mutation': 'sat_worst_mutation',
            'worst_mut_aa': 'sat_worst_mut_aa',
            'worst_ddG': 'sat_worst_ddG'
        }, inplace=True)

    atlas = atlas.merge(sat_merge, on='position', how='left')

    # Clinical data from v2
    atlas = atlas.merge(clinical, on='position', how='left')

    # ── Mutational sensitivity score ─────────────────────────────────────────
    # Normalised: how sensitive is this position to any mutation?
    max_possible_ddG = atlas['sat_max_ddG'].quantile(0.99)  # 99th percentile cap
    atlas['mutational_sensitivity'] = np.clip(
        atlas['sat_max_ddG'] / max_possible_ddG, 0, 1
    )

    # ── 6. Risk classification v3 ────────────────────────────────────────────
    print("\n[6/7] Calculating risk classification v3...")

    def classify_risk_v3(row):
        """
        Multi-factor risk classification:
        - sat_max_ddG: positional mutational sensitivity
        - plddt: structural confidence
        - domain penetrance: clinical outcome data
        - interface proximity: functional importance
        """
        score = 0

        # Factor 1: Mutational sensitivity (0-40 points)
        max_ddg = row.get('sat_max_ddG', 0)
        if pd.isna(max_ddg):
            max_ddg = 0
        if max_ddg > 10:
            score += 40
        elif max_ddg > 5:
            score += 30
        elif max_ddg > 2:
            score += 20
        elif max_ddg > 1:
            score += 10

        # Factor 2: pLDDT confidence (0-20 points)
        plddt_val = row.get('plddt', 50)
        if pd.isna(plddt_val):
            plddt_val = 50
        if plddt_val >= 90:
            score += 20  # High confidence = we trust the ddG
        elif plddt_val >= 70:
            score += 15
        elif plddt_val >= 50:
            score += 5

        # Factor 3: Domain penetrance (0-25 points)
        pen = row.get('domain_penetrance_by_60_pct', 15)
        if pd.isna(pen):
            pen = 15
        if pen > 30:
            score += 25
        elif pen > 20:
            score += 18
        elif pen > 15:
            score += 10
        else:
            score += 5

        # Factor 4: Interface proximity (0-15 points)
        if row.get('at_pcsk9_interface', False):
            score += 15
        elif row.get('at_apob_interface', False):
            score += 15
        else:
            pcsk9_d = row.get('pcsk9_interface_distance', 50)
            apob_d = row.get('apob_interface_distance', 50)
            if pd.isna(pcsk9_d):
                pcsk9_d = 50
            if pd.isna(apob_d):
                apob_d = 50
            min_d = min(pcsk9_d, apob_d)
            if min_d < 5:
                score += 12
            elif min_d < 10:
                score += 8
            elif min_d < 15:
                score += 4

        # Classify
        if score >= 70:
            return 'Critical'
        elif score >= 50:
            return 'High'
        elif score >= 35:
            return 'Moderate'
        elif score >= 20:
            return 'Low-Moderate'
        else:
            return 'Low'

    atlas['risk_classification'] = atlas.apply(classify_risk_v3, axis=1)
    atlas['risk_score'] = atlas.apply(
        lambda row: classify_risk_v3.__code__.co_consts[0] if False else 0, axis=1
    )

    # Recalculate risk_score properly
    def calc_risk_score(row):
        score = 0
        max_ddg = row.get('sat_max_ddG', 0) or 0
        if max_ddg > 10: score += 40
        elif max_ddg > 5: score += 30
        elif max_ddg > 2: score += 20
        elif max_ddg > 1: score += 10

        plddt_val = row.get('plddt', 50) or 50
        if plddt_val >= 90: score += 20
        elif plddt_val >= 70: score += 15
        elif plddt_val >= 50: score += 5

        pen = row.get('domain_penetrance_by_60_pct', 15) or 15
        if pen > 30: score += 25
        elif pen > 20: score += 18
        elif pen > 15: score += 10
        else: score += 5

        if row.get('at_pcsk9_interface', False): score += 15
        elif row.get('at_apob_interface', False): score += 15
        return score

    atlas['risk_score'] = atlas.apply(calc_risk_score, axis=1)

    print("   Risk classification distribution:")
    for cat in ['Critical', 'High', 'Moderate', 'Low-Moderate', 'Low']:
        n = (atlas['risk_classification'] == cat).sum()
        pct = n / len(atlas) * 100
        print(f"     {cat:15s}: {n:4d} ({pct:5.1f}%)")

    # ── 7. Rebuild SSS v3 catalogue ──────────────────────────────────────────
    print("\n[7/7] Rebuilding SSS v3 catalogue...")

    sss_v2 = pd.read_csv(SSS_V2_CAT)
    sss_v3 = sss_v2.copy()

    # For LDLR missense variants with a position, update ddG_score using
    # position-specific saturation data
    ldlr_mask = (sss_v3['gene'] == 'LDLR') & sss_v3['position'].notna()
    updated = 0

    for idx in sss_v3[ldlr_mask].index:
        pos = int(sss_v3.loc[idx, 'position'])
        pos_data = pos_stats[pos_stats['position'] == pos]

        if len(pos_data) > 0 and pos_data.iloc[0]['n_computed'] > 0:
            max_ddg = pos_data.iloc[0]['max_ddG']
            mean_ddg = pos_data.iloc[0]['mean_ddG']

            # Update ddG_score: normalised to 0-1 scale
            # Using position-specific max_ddG relative to global distribution
            ddg_score = np.clip(max_ddg / 20, 0, 1)  # 20 kcal/mol = max severity

            # Recalculate SSS v3
            domain_risk = sss_v3.loc[idx, 'domain_risk']
            type_severity = sss_v3.loc[idx, 'type_severity']
            interface_prox = sss_v3.loc[idx, 'interface_proximity']
            plddt_disrupt = sss_v3.loc[idx, 'plddt_disruption']

            # SSS v3 = weighted combination
            # Domain risk: 25%, Type severity: 20%, ddG: 30%, pLDDT: 15%, Interface: 10%
            sss_v3_score = (0.25 * domain_risk +
                           0.20 * type_severity +
                           0.30 * ddg_score +
                           0.15 * plddt_disrupt +
                           0.10 * interface_prox)

            sss_v3.loc[idx, 'ddG_score'] = round(ddg_score, 4)
            sss_v3.loc[idx, 'sss_v3'] = round(sss_v3_score, 4)
            sss_v3.loc[idx, 'sat_max_ddG'] = round(max_ddg, 2)
            sss_v3.loc[idx, 'sat_mean_ddG'] = round(mean_ddg, 2)
            sss_v3.loc[idx, 'mutational_sensitivity'] = round(
                np.clip(max_ddg / atlas['sat_max_ddG'].quantile(0.99), 0, 1), 4
            )
            updated += 1

    # For variants without position-specific data, keep sss_v2 as sss_v3
    sss_v3['sss_v3'] = sss_v3['sss_v3'] if 'sss_v3' in sss_v3.columns else sss_v3['sss_v2']
    sss_v3.loc[sss_v3['sss_v3'].isna(), 'sss_v3'] = sss_v3.loc[sss_v3['sss_v3'].isna(), 'sss_v2']

    print(f"   Updated {updated} LDLR variants with saturation-derived ddG scores")
    print(f"   SSS v3 range: {sss_v3['sss_v3'].min():.4f} - {sss_v3['sss_v3'].max():.4f}")
    print(f"   SSS v3 mean:  {sss_v3['sss_v3'].mean():.4f}")
    print(f"   SSS v2->v3 correlation: {sss_v3[['sss_v2','sss_v3']].corr().iloc[0,1]:.4f}")

    # ── Save outputs ─────────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("SAVING OUTPUTS")
    print("=" * 80)

    # Atlas v3
    atlas.to_csv(ATLAS_V3, index=False)
    print(f"\n   Atlas v3: {ATLAS_V3}")
    print(f"   Shape: {atlas.shape}")

    # SSS v3 catalogue
    sss_v3.to_csv(SSS_V3_CAT, index=False)
    print(f"\n   SSS v3 catalogue: {SSS_V3_CAT}")
    print(f"   Shape: {sss_v3.shape}")

    # Summary statistics table (for publication)
    summary_rows = []
    for domain, (start, end) in DOMAIN_RANGES.items():
        d = atlas[(atlas['position'] >= start) & (atlas['position'] <= end)]
        d_sat = d[d['sat_n_computed'] > 0]
        summary_rows.append({
            'Domain': domain,
            'Residue_range': f"{start}-{end}",
            'N_residues': len(d),
            'Mean_pLDDT': round(d['plddt'].mean(), 1) if d['plddt'].notna().any() else np.nan,
            'Mean_max_ddG': round(d['sat_max_ddG'].mean(), 2) if d['sat_max_ddG'].notna().any() else np.nan,
            'Pct_critical': round((d['risk_classification'] == 'Critical').mean() * 100, 1),
            'Pct_high': round((d['risk_classification'] == 'High').mean() * 100, 1),
            'Pct_destabilizing': round(
                d['sat_n_destabilizing'].sum() / max(d['sat_n_computed'].sum(), 1) * 100, 1
            ) if d['sat_n_computed'].sum() > 0 else np.nan,
            'N_interface_residues': int(d['at_pcsk9_interface'].sum() + d['at_apob_interface'].sum()),
            'Penetrance_by_60': DOMAIN_PENETRANCE.get(domain, np.nan),
            'ASCVD_rate': DOMAIN_ASCVD_RATE.get(domain, np.nan),
            'N_known_variants': int(d['n_known_variants'].sum()) if d['n_known_variants'].notna().any() else 0,
            'N_patients': int(d['n_patients'].sum()) if d['n_patients'].notna().any() else 0,
            'FoldX_coverage_pct': round(len(d_sat) / len(d) * 100, 1),
        })

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(SUMMARY_STATS, index=False)
    print(f"\n   Summary statistics: {SUMMARY_STATS}")

    # Heatmap data (position × amino acid ddG matrix)
    print("\n   Generating heatmap data...")
    heatmap = sat.pivot_table(index='position', columns='mut_aa', values='ddG', aggfunc='first')
    heatmap.to_csv(HEATMAP_DATA)
    print(f"   Heatmap: {HEATMAP_DATA} ({heatmap.shape})")

    # ── Final summary ────────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("ATLAS v3 COMPLETE")
    print("=" * 80)
    corr_val = sss_v3[['sss_v2','sss_v3']].corr().iloc[0,1]
    print(f"   Positions:                860 (complete LDLR)")
    print(f"   Saturation mutations:     {sat['ddG'].notna().sum()} valid / {len(sat)} total")
    print(f"   FoldX-computed positions: {len(valid_positions)} ({len(valid_positions)/860*100:.1f}%)")
    print(f"   Domain-imputed positions: {int(imputed_count)} ({imputed_count/860*100:.1f}%)")
    print(f"   Critical:      {(atlas['risk_classification']=='Critical').sum()}")
    print(f"   High:          {(atlas['risk_classification']=='High').sum()}")
    print(f"   Moderate:      {(atlas['risk_classification']=='Moderate').sum()}")
    print(f"   Low-Moderate:  {(atlas['risk_classification']=='Low-Moderate').sum()}")
    print(f"   Low:           {(atlas['risk_classification']=='Low').sum()}")
    print(f"   SSS v3 variants: {len(sss_v3)} ({updated} updated with saturation data)")
    print(f"   SSS v2->v3 correlation: {corr_val:.4f}")


if __name__ == "__main__":
    main()

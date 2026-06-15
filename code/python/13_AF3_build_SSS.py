#!/usr/bin/env python3
"""
13_AF3_build_SSS.py
Phase 3: Build the Structural Severity Score (SSS)

Combines:
  1. FoldX ddG (structural stability)
  2. pLDDT at variant site (prediction confidence)
  3. LDLR domain severity (biological importance)
  4. Interface proximity (distance to PCSK9 binding site)
  5. Variant type preset (nonsense/frameshift = 1.0)

Outputs a per-variant SSS and maps to patients in Wales + UKB cohorts.

Inputs:
  alphafold/analysis/foldx_ddg_results.csv
  alphafold/analysis/LDLR_wildtype_plddt_per_residue.csv
  alphafold/analysis/variants_mapped_to_structure.csv
  alphafold/af3_variant_catalogue.csv
  alphafold/af3_variant_patient_map.csv
  DRAGON_3.csv, WALES_FH_CLEANED (1) - Copy.csv

Outputs:
  alphafold/analysis/structural_severity_scores.csv
  alphafold/analysis/sss_patient_level_wales.csv
"""

import os
import csv
import math
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
AF_DIR = os.path.join(BASE, "alphafold")
ANALYSIS_DIR = os.path.join(AF_DIR, "analysis")

# ── LDLR Domain severity weights (expert-derived) ────────────────
# Based on functional importance for LDL uptake
DOMAIN_SEVERITY = {
    'Signal peptide': 0.3,
    'Ligand-binding R1': 1.0,
    'Ligand-binding R2': 1.0,
    'Ligand-binding R3': 1.0,
    'Ligand-binding R4': 1.0,
    'Ligand-binding R5': 1.0,
    'Ligand-binding R6': 1.0,
    'Ligand-binding R7': 1.0,
    'Ligand-binding': 1.0,  # generic
    'EGF-like A': 0.9,      # PCSK9 binding site
    'EGF-like B': 0.7,
    'EGF-like C': 0.7,
    'Beta-propeller': 0.6,
    'Beta-prop': 0.6,
    'Linker (EGF-C/O-linked)': 0.4,
    'Linker': 0.4,
    'O-linked sugar': 0.3,
    'O-linked-sugar': 0.3,
    'Transmembrane': 0.5,
    'Cytoplasmic': 0.4,
    'Unknown': 0.3,
    '?': 0.3,
    # Non-LDLR genes
    'EGF-precursor-homology': 0.7,
    'ApoB-RBD': 0.8,
    'PCSK9-catalytic': 0.8,
    'PCSK9-prodomain': 0.5,
    'PCSK9-V-domain': 0.4,
}


def get_domain_from_position(pos):
    """Map residue position to LDLR domain."""
    domains = [
        ('Signal peptide', 1, 21),
        ('Ligand-binding R1', 22, 83),
        ('Ligand-binding R2', 84, 120),
        ('Ligand-binding R3', 121, 159),
        ('Ligand-binding R4', 160, 197),
        ('Ligand-binding R5', 198, 236),
        ('Ligand-binding R6', 237, 275),
        ('Ligand-binding R7', 276, 313),
        ('EGF-like A', 314, 353),
        ('EGF-like B', 354, 393),
        ('Beta-propeller', 394, 631),
        ('EGF-like C', 632, 671),
        ('O-linked sugar', 693, 750),
        ('Transmembrane', 751, 788),
        ('Cytoplasmic', 789, 860),
    ]
    for name, s, e in domains:
        if s <= pos <= e:
            return name
    if 672 <= pos <= 692:
        return 'Linker (EGF-C/O-linked)'
    return 'Unknown'


def sigmoid(x, midpoint=3.0, steepness=1.0):
    """Sigmoid transform to map ddG to 0-1 range."""
    return 1.0 / (1.0 + math.exp(-steepness * (x - midpoint)))


def compute_sss(ddg, plddt, domain_severity, interface_dist=None, variant_type='missense'):
    """
    Compute Structural Severity Score (SSS).

    Components (pre-specified weights, NOT data-driven):
      - ddG component (40%): sigmoid-transformed FoldX ddG
      - Domain component (25%): expert-assigned domain severity
      - Interface component (15%): proximity to PCSK9 binding interface
      - pLDDT reliability (10%): confidence weighting
      - Variant type (10%): missense vs cysteine-disrupting

    Returns SSS in [0, 1] range where 1 = most severe.
    """
    # 1. ddG component: sigmoid transform, midpoint at 3 kcal/mol
    if ddg is not None:
        ddg_score = sigmoid(ddg, midpoint=3.0, steepness=0.8)
    else:
        ddg_score = 0.5  # unknown = neutral prior

    # 2. Domain severity (already 0-1)
    domain_score = domain_severity

    # 3. Interface proximity
    if interface_dist is not None:
        if interface_dist < 5:
            interface_score = 1.0
        elif interface_dist < 10:
            interface_score = 0.7
        elif interface_dist < 15:
            interface_score = 0.4
        else:
            interface_score = 0.1
    else:
        interface_score = 0.1  # unknown = assume far

    # 4. pLDDT reliability weight
    # High pLDDT = we trust the ddG more
    # Low pLDDT = ddG less reliable, rely more on domain/type
    if plddt is not None and plddt > 0:
        plddt_weight = min(plddt / 100.0, 1.0)
    else:
        plddt_weight = 0.5

    # 5. Variant type modifier
    type_score = 0.5  # baseline for missense
    if variant_type == 'cysteine_disrupting':
        type_score = 0.8
    elif variant_type == 'charged_to_neutral':
        type_score = 0.6

    # Combine with pre-specified weights
    # Adjust ddG weight by pLDDT confidence
    w_ddg = 0.40 * plddt_weight + 0.10 * (1 - plddt_weight)  # if low pLDDT, downweight ddG
    w_domain = 0.25
    w_interface = 0.15
    w_type = 0.10
    w_plddt_bonus = 0.10  # residual weight redistributed

    # Normalize weights
    w_total = w_ddg + w_domain + w_interface + w_type + w_plddt_bonus

    sss = (w_ddg * ddg_score +
           w_domain * domain_score +
           w_interface * interface_score +
           w_type * type_score +
           w_plddt_bonus * (1 - plddt_weight)) / w_total  # low pLDDT adds uncertainty penalty

    return round(sss, 4)


def main():
    print("=" * 70)
    print("  CALON-Structure: Structural Severity Score (SSS) Construction")
    print("  Phase 3 - Pre-specified weights (no data-driven fitting)")
    print("=" * 70)

    # ── Load FoldX ddG results ────────────────────────────────────
    ddg_data = {}
    ddg_file = os.path.join(ANALYSIS_DIR, "foldx_ddg_results.csv")
    if os.path.exists(ddg_file):
        with open(ddg_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                pos = int(row['position'])
                ddg_data[pos] = {
                    'ddG': float(row['ddG_mean']),
                    'mutation': row['foldx_mutation'],
                    'effect': row['effect'],
                }
        print(f"  Loaded {len(ddg_data)} ddG values")

    # ── Load pLDDT per residue ────────────────────────────────────
    plddt_data = {}
    plddt_file = os.path.join(ANALYSIS_DIR, "LDLR_wildtype_plddt_per_residue.csv")
    if os.path.exists(plddt_file):
        with open(plddt_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                plddt_data[int(row['residue_num'])] = float(row['plddt'])
        print(f"  Loaded {len(plddt_data)} pLDDT values")

    # ── Load variant catalogue ────────────────────────────────────
    catalogue_file = os.path.join(AF_DIR, "af3_variant_catalogue.csv")
    variants = []
    if os.path.exists(catalogue_file):
        with open(catalogue_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                variants.append(row)
        print(f"  Loaded {len(variants)} variants from catalogue")

    # ── PCSK9 interface distances (from our analysis) ─────────────
    # Pre-computed from LDLR-PCSK9 complex v2
    # Only positions near the interface (< 15 Angstrom)
    interface_distances = {
        46: 7.5, 48: 5.7, 49: 7.5, 50: 7.6, 54: 7.5, 55: 7.0,
        58: 5.0, 59: 5.0, 60: 5.9, 78: 6.6, 79: 5.7, 94: 12.8,
        96: 7.3, 318: 6.4, 319: 6.8, 320: 6.5, 321: 4.6, 322: 6.7,
        324: 7.6, 326: 5.5, 327: 6.1, 328: 5.0, 329: 5.0, 330: 4.0,
        331: 5.5, 332: 7.2, 340: 9.3, 343: 13.2, 350: 11.2, 354: 11.1,
    }

    # ── Compute SSS for each variant ──────────────────────────────
    print(f"\n{'='*70}")
    print(f"  Computing SSS for all variants")
    print(f"{'='*70}")

    sss_results = []

    for var in variants:
        variant_id = var.get('variant_id', '')
        gene = var.get('gene', '')
        vtype = var.get('variant_type', '')

        # Get protein position
        pos_str = var.get('protein_position', '')
        try:
            pos = int(float(pos_str))
        except (ValueError, TypeError):
            pos = None

        # Preset SSS for loss-of-function variants
        sss_preset = var.get('sss_preset', '')
        if sss_preset and sss_preset not in ('', 'NA', 'nan'):
            try:
                sss = float(sss_preset)
                sss_results.append({
                    'variant_id': variant_id,
                    'gene': gene,
                    'variant_type': vtype,
                    'protein_position': pos if pos else '',
                    'domain': var.get('domain', ''),
                    'sss': sss,
                    'sss_source': 'preset_lof',
                    'ddG': '',
                    'plddt': '',
                    'interface_dist': '',
                    'domain_severity': '',
                })
                continue
            except ValueError:
                pass

        # For missense variants in LDLR
        if gene == 'LDLR' and pos is not None:
            # Get ddG
            ddg_info = ddg_data.get(pos, None)
            ddg = ddg_info['ddG'] if ddg_info else None

            # Get pLDDT
            plddt = plddt_data.get(pos, None)

            # Get domain
            domain = get_domain_from_position(pos)
            dom_sev = DOMAIN_SEVERITY.get(domain, 0.5)

            # Get interface distance
            iface_dist = interface_distances.get(pos, None)

            # Determine variant subtype
            subtype = 'missense'
            # Check if cysteine is disrupted
            if ddg_info and ddg_info['mutation'][0] == 'C':
                subtype = 'cysteine_disrupting'

            sss = compute_sss(ddg, plddt, dom_sev, iface_dist, subtype)

            sss_results.append({
                'variant_id': variant_id,
                'gene': gene,
                'variant_type': vtype,
                'protein_position': pos,
                'domain': domain,
                'sss': sss,
                'sss_source': 'foldx' if ddg is not None else 'domain_only',
                'ddG': f"{ddg:.2f}" if ddg is not None else '',
                'plddt': f"{plddt:.1f}" if plddt is not None else '',
                'interface_dist': f"{iface_dist:.1f}" if iface_dist is not None else '',
                'domain_severity': f"{dom_sev:.2f}",
            })

        elif gene == 'APOB':
            # ApoB: use domain severity preset
            dom_sev = 0.8  # RBD
            sss = dom_sev * 0.5  # no structural data for most
            if 'R3527' in variant_id or 'c.10580' in variant_id:
                sss = 0.65  # well-known pathogenic

            sss_results.append({
                'variant_id': variant_id,
                'gene': gene,
                'variant_type': vtype,
                'protein_position': pos if pos else '',
                'domain': var.get('domain', 'ApoB-RBD'),
                'sss': sss,
                'sss_source': 'gene_preset',
                'ddG': '',
                'plddt': '',
                'interface_dist': '',
                'domain_severity': f"{dom_sev:.2f}",
            })

        elif gene == 'PCSK9':
            dom_sev = 0.4  # gain-of-function, different mechanism
            sss = 0.45  # moderate - PCSK9 GoF variants have variable penetrance

            sss_results.append({
                'variant_id': variant_id,
                'gene': gene,
                'variant_type': vtype,
                'protein_position': pos if pos else '',
                'domain': var.get('domain', 'PCSK9'),
                'sss': sss,
                'sss_source': 'gene_preset',
                'ddG': '',
                'plddt': '',
                'interface_dist': '',
                'domain_severity': f"{dom_sev:.2f}",
            })

        else:
            # Unknown gene or no position
            sss_results.append({
                'variant_id': variant_id,
                'gene': gene,
                'variant_type': vtype,
                'protein_position': pos if pos else '',
                'domain': var.get('domain', ''),
                'sss': 0.5,  # uninformative prior
                'sss_source': 'default',
                'ddG': '',
                'plddt': '',
                'interface_dist': '',
                'domain_severity': '',
            })

    # ── Print summary ─────────────────────────────────────────────
    print(f"\n  Total variants with SSS: {len(sss_results)}")

    # Distribution
    sss_vals = [r['sss'] for r in sss_results]
    bins = [(0, 0.2, 'Very low'), (0.2, 0.4, 'Low'), (0.4, 0.6, 'Moderate'),
            (0.6, 0.8, 'High'), (0.8, 1.0, 'Very high')]
    print(f"\n  SSS Distribution:")
    for lo, hi, label in bins:
        count = sum(1 for s in sss_vals if lo <= s < hi or (hi == 1.0 and s == 1.0))
        print(f"    {label:<12s} ({lo:.1f}-{hi:.1f}): {count:>4d} variants")

    # By source
    sources = {}
    for r in sss_results:
        src = r['sss_source']
        sources[src] = sources.get(src, 0) + 1
    print(f"\n  SSS source breakdown:")
    for src, count in sorted(sources.items()):
        print(f"    {src:<15s}: {count}")

    # Top destabilizing variants (with FoldX data)
    foldx_results = [r for r in sss_results if r['sss_source'] == 'foldx']
    foldx_results.sort(key=lambda x: -x['sss'])

    print(f"\n  Top 15 most severe LDLR missense variants (FoldX-modelled):")
    print(f"  {'Variant':<25s} {'Pos':>5s} {'Domain':<20s} {'ddG':>7s} {'pLDDT':>6s} {'SSS':>6s}")
    print(f"  {'-'*25} {'-'*5} {'-'*20} {'-'*7} {'-'*6} {'-'*6}")
    for r in foldx_results[:15]:
        print(f"  {r['variant_id']:<25s} {r['protein_position']:>5d} {r['domain']:<20s} {r['ddG']:>7s} {r['plddt']:>6s} {r['sss']:>6.3f}")

    # ── Save SSS results ──────────────────────────────────────────
    out_file = os.path.join(ANALYSIS_DIR, "structural_severity_scores.csv")
    with open(out_file, 'w', newline='') as f:
        if sss_results:
            writer = csv.DictWriter(f, fieldnames=sss_results[0].keys())
            writer.writeheader()
            writer.writerows(sss_results)
    print(f"\n  Saved: {out_file}")

    # ── Map SSS to patients ───────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"  Mapping SSS to patient level")
    print(f"{'='*70}")

    # Build variant -> SSS lookup
    sss_lookup = {r['variant_id']: r['sss'] for r in sss_results}

    # Load patient-variant map
    patient_map_file = os.path.join(AF_DIR, "af3_variant_patient_map.csv")
    if os.path.exists(patient_map_file):
        patient_sss = []
        with open(patient_map_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                vid = row.get('variant_id', '')
                sss = sss_lookup.get(vid, None)
                if sss is not None:
                    patient_sss.append({
                        'patient_id': row.get('patient_id', ''),
                        'variant_id': vid,
                        'gene': row.get('gene', ''),
                        'sss': sss,
                        'source': row.get('source', ''),
                    })

        # For patients with multiple variants, use max SSS
        patient_max_sss = {}
        for ps in patient_sss:
            pid = ps['patient_id']
            if pid not in patient_max_sss or ps['sss'] > patient_max_sss[pid]['sss']:
                patient_max_sss[pid] = ps

        print(f"  Patient-variant pairs: {len(patient_sss)}")
        print(f"  Unique patients with SSS: {len(patient_max_sss)}")

        # Save patient-level SSS
        patient_file = os.path.join(ANALYSIS_DIR, "sss_patient_level_wales.csv")
        with open(patient_file, 'w', newline='') as f:
            if patient_sss:
                writer = csv.DictWriter(f, fieldnames=['patient_id', 'variant_id', 'gene', 'sss', 'source'])
                writer.writeheader()
                for pid in sorted(patient_max_sss.keys()):
                    writer.writerow(patient_max_sss[pid])
        print(f"  Saved: {patient_file}")

        # SSS distribution for patients
        patient_sss_vals = [v['sss'] for v in patient_max_sss.values()]
        print(f"\n  Patient-level SSS distribution:")
        for lo, hi, label in bins:
            count = sum(1 for s in patient_sss_vals if lo <= s < hi or (hi == 1.0 and s == 1.0))
            pct = count / len(patient_sss_vals) * 100 if patient_sss_vals else 0
            print(f"    {label:<12s}: {count:>5d} ({pct:>5.1f}%)")
    else:
        print(f"  !! Patient map not found. Run 10_AF3_extract_variants.py first.")

    print(f"\n  SSS construction complete!")
    print(f"  Key design decisions:")
    print(f"    - Pre-specified weights (NOT data-driven) to avoid circularity")
    print(f"    - ddG sigmoid midpoint at 3 kcal/mol (standard FoldX threshold)")
    print(f"    - Domain severity from LDLR functional literature")
    print(f"    - Interface distances from AF3 LDLR-PCSK9 complex")
    print(f"    - Loss-of-function variants preset to SSS=1.0 (worst)")


if __name__ == "__main__":
    main()

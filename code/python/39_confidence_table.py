#!/usr/bin/env python3
"""
===================================================================================
LDLR 860-POSITION CONFIDENCE ASSESSMENT TABLE
Supplementary Table S1 — Transparent reporting of data availability per position

Merges: pLDDT, FoldX ddG, Rosetta energy, Rosetta stabilize-PM, Islam et al.
        functional data, atlas clinical data, interface distances

Each position receives a confidence tier:
  HIGH     = position-specific FoldX ddG + pLDDT > 70 + clinical data
  MODERATE = position-specific Rosetta energy + pLDDT > 50 + domain inference
  LOW      = domain-level inference only (no position-specific structural data)
  MINIMAL  = pLDDT < 50, no structural data, domain inference only

Author: Dr Nader Genedy
===================================================================================
"""

import csv, os
from collections import defaultdict

BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"
AF = f"{BASE}/alphafold/analysis"

def load(path):
    try:
        with open(path, encoding='utf-8-sig') as f:
            return list(csv.DictReader(f))
    except Exception as e:
        print(f"  WARNING: Cannot load {path}: {e}")
        return []

def sf(x):
    if x is None or str(x).strip() == '':
        return None
    try:
        return float(str(x).strip())
    except:
        return None

print("="*70)
print("  LDLR 860-POSITION CONFIDENCE ASSESSMENT TABLE")
print("="*70)

# ─── Load all data sources ───
print("\nLoading data sources...")
plddt = load(f"{AF}/LDLR_wildtype_plddt_per_residue.csv")
atlas = load(f"{AF}/LDLR_genotype_phenotype_atlas.csv")
foldx_28 = load(f"{AF}/foldx_ddg_results.csv")
foldx_162 = load(f"{AF}/foldx_ddg_results_full.csv")
rosetta = load(f"{AF}/rosetta_energy_breakdown.csv")
stab_pm = load(f"{AF}/rosetta_stabilize_pm_results.csv")
islam = load(f"{AF}/islam_et_al_315_variants.csv")
sss_all = load(f"{AF}/structural_severity_scores.csv")
interface = load(f"{AF}/interface_distance_analysis.csv")
nobel = load(f"{AF}/nobel_multimodal_integration.csv")
comp = load(f"{AF}/comprehensive_sss_analysis.csv")

print(f"  pLDDT: {len(plddt)} residues")
print(f"  Atlas: {len(atlas)} positions")
print(f"  FoldX (original): {len(foldx_28)} variants")
print(f"  FoldX (expanded): {len(foldx_162)} variants")
print(f"  Rosetta energy: {len(rosetta)} residues")
print(f"  Rosetta stabilize-PM: {len(stab_pm)} positions")
print(f"  Islam et al.: {len(islam)} variants")
print(f"  SSS scores: {len(sss_all)} variants")
print(f"  Interface: {len(interface)} variants")
print(f"  Nobel integration: {len(nobel)} variants")
print(f"  Comprehensive clinical: {len(comp)} patients")

# ─── Build lookup dictionaries ───
print("\nBuilding position-level lookups...")

# pLDDT by position
plddt_by_pos = {}
for r in plddt:
    pos = sf(r.get('residue_num'))
    if pos is not None:
        plddt_by_pos[int(pos)] = {
            'plddt': sf(r.get('plddt')),
            'domain': r.get('domain', ''),
            'quality': r.get('quality_category', '')
        }

# Atlas by position
atlas_by_pos = {}
for r in atlas:
    pos = sf(r.get('position'))
    if pos is not None:
        atlas_by_pos[int(pos)] = r

# FoldX ddG by position (collect all variants at each position)
foldx_by_pos = defaultdict(list)
for r in foldx_162:
    pos = sf(r.get('position'))
    if pos is not None:
        foldx_by_pos[int(pos)].append({
            'variant': r.get('variant_id', ''),
            'wt': r.get('wt_aa', ''),
            'mut': r.get('mut_aa', ''),
            'ddG': sf(r.get('ddG_kcal_mol')),
            'effect': r.get('effect', ''),
            'n_patients': sf(r.get('n_patients')) or 0,
            'source': 'original' if any(
                r.get('variant_id') == o.get('variant_id') for o in foldx_28
            ) else 'expanded'
        })

# Rosetta energy by position
rosetta_by_pos = {}
for i, r in enumerate(rosetta):
    pos = i + 1  # 1-indexed
    rosetta_by_pos[pos] = {
        'total': sf(r.get('total')),
        'fa_atr': sf(r.get('fa_atr')),
        'fa_rep': sf(r.get('fa_rep')),
        'fa_sol': sf(r.get('fa_sol')),
        'hbond_sr_bb': sf(r.get('hbond_sr_bb')),
        'hbond_lr_bb': sf(r.get('hbond_lr_bb')),
        'hbond_bb_sc': sf(r.get('hbond_bb_sc')),
        'hbond_sc': sf(r.get('hbond_sc')),
    }

# Rosetta stabilize-PM by position
stab_by_pos = {}
for r in stab_pm:
    pos = sf(r.get('position'))
    if pos is not None:
        stab_by_pos[int(pos)] = {
            'best_mutation': r.get('best_mutation', ''),
            'best_score': sf(r.get('best_score')),
            'worst_mutation': r.get('worst_mutation', ''),
            'worst_score': sf(r.get('worst_score')),
            'n_stabilising': sf(r.get('n_stabilising')),
            'n_destabilising': sf(r.get('n_destabilising')),
        }

# Islam et al. by position (may have multiple variants per position)
islam_by_pos = defaultdict(list)
for r in islam:
    pos = sf(r.get('position'))
    if pos is not None:
        islam_by_pos[int(pos)].append({
            'variant': r.get('variant', ''),
            'activity_pct': sf(r.get('activity_pct')),
            'functional_group': r.get('functional_group', ''),
            'mean_gfp': sf(r.get('mean_gfp')),
            'clinvar': r.get('clinvar', ''),
        })

# Clinical data by position (from comp — aggregate patients per position)
clinical_by_pos = defaultdict(lambda: {
    'n_patients': 0, 'ascvd_events': 0, 'ldl_values': [], 'apob_values': [],
    'lpa_values': [], 'treatment_responses': [], 'variants': set()
})
for r in comp:
    if r.get('gene') != 'LDLR':
        continue
    # Extract position from mutation string
    mut = r.get('mutation', '')
    # Try to find position from SSS data
    pos = None
    for s in sss_all:
        if s.get('variant_id') == f"LDLR:{mut}" or s.get('variant_id') == mut:
            pos = sf(s.get('protein_position'))
            break
    if pos is None:
        # Try to extract from atlas matching
        continue
    pos = int(pos)
    clinical_by_pos[pos]['n_patients'] += 1
    clinical_by_pos[pos]['variants'].add(mut)
    if sf(r.get('ascvd')) == 1:
        clinical_by_pos[pos]['ascvd_events'] += 1
    ldl = sf(r.get('ldl1'))
    if ldl is not None:
        clinical_by_pos[pos]['ldl_values'].append(ldl)
    apob = sf(r.get('apob'))
    if apob is not None:
        clinical_by_pos[pos]['apob_values'].append(apob)
    lpa = sf(r.get('lpa'))
    if lpa is not None:
        clinical_by_pos[pos]['lpa_values'].append(lpa)
    resp = sf(r.get('ldl_pct_change'))
    if resp is not None:
        clinical_by_pos[pos]['treatment_responses'].append(resp)

# Interface distances from atlas
interface_by_pos = {}
for r in atlas:
    pos = sf(r.get('position'))
    if pos is not None:
        interface_by_pos[int(pos)] = {
            'pcsk9_dist': sf(r.get('pcsk9_interface_distance')),
            'apob_dist': sf(r.get('apob_interface_distance')),
            'at_pcsk9': r.get('at_pcsk9_interface', '') == 'Yes',
            'at_apob': r.get('at_apob_interface', '') == 'Yes',
        }

# ─── Check for saturation mutagenesis progress ───
sat_file = f"{AF}/foldx_saturation_mutagenesis.csv"
sat_by_pos = defaultdict(list)
if os.path.exists(sat_file):
    sat_data = load(sat_file)
    print(f"  Saturation mutagenesis: {len(sat_data)} mutations computed")
    for r in sat_data:
        pos = sf(r.get('position'))
        if pos is not None:
            sat_by_pos[int(pos)].append(sf(r.get('ddG')))

# ─── Build the 860-row confidence table ───
print("\nBuilding confidence table...")

rows = []
confidence_counts = defaultdict(int)

for pos in range(1, 861):
    row = {'position': pos}

    # 1. Wildtype amino acid and domain
    if pos in plddt_by_pos:
        row['wildtype_aa'] = atlas_by_pos.get(pos, {}).get('wildtype_aa', '')
        row['domain'] = plddt_by_pos[pos]['domain']
        row['plddt'] = plddt_by_pos[pos]['plddt']
        row['plddt_quality'] = plddt_by_pos[pos]['quality']
    elif pos in atlas_by_pos:
        row['wildtype_aa'] = atlas_by_pos[pos].get('wildtype_aa', '')
        row['domain'] = atlas_by_pos[pos].get('domain', '')
        row['plddt'] = sf(atlas_by_pos[pos].get('plddt_confidence'))
        row['plddt_quality'] = ''
    else:
        row['wildtype_aa'] = ''
        row['domain'] = ''
        row['plddt'] = None
        row['plddt_quality'] = ''

    # 2. FoldX data
    if pos in foldx_by_pos:
        variants = foldx_by_pos[pos]
        row['has_foldx'] = True
        row['n_foldx_variants'] = len(variants)
        ddgs = [v['ddG'] for v in variants if v['ddG'] is not None]
        row['foldx_max_ddG'] = max(ddgs) if ddgs else None
        row['foldx_min_ddG'] = min(ddgs) if ddgs else None
        row['foldx_mean_ddG'] = sum(ddgs)/len(ddgs) if ddgs else None
        row['foldx_variants'] = '; '.join(f"{v['wt']}{pos}{v['mut']}={v['ddG']:.2f}"
                                           for v in variants if v['ddG'] is not None)
        row['foldx_source'] = 'original' if any(v['source']=='original' for v in variants) else 'expanded'
    else:
        row['has_foldx'] = False
        row['n_foldx_variants'] = 0
        row['foldx_max_ddG'] = None
        row['foldx_min_ddG'] = None
        row['foldx_mean_ddG'] = None
        row['foldx_variants'] = ''
        row['foldx_source'] = ''

    # 3. Saturation mutagenesis data
    if pos in sat_by_pos:
        row['has_saturation'] = True
        row['n_saturation_muts'] = len(sat_by_pos[pos])
        svals = [v for v in sat_by_pos[pos] if v is not None]
        row['sat_max_ddG'] = max(svals) if svals else None
        row['sat_min_ddG'] = min(svals) if svals else None
        row['sat_mean_ddG'] = sum(svals)/len(svals) if svals else None
    else:
        row['has_saturation'] = False
        row['n_saturation_muts'] = 0
        row['sat_max_ddG'] = None
        row['sat_min_ddG'] = None
        row['sat_mean_ddG'] = None

    # 4. Rosetta energy
    if pos in rosetta_by_pos:
        row['has_rosetta_energy'] = True
        row['rosetta_total'] = rosetta_by_pos[pos]['total']
        row['rosetta_fa_atr'] = rosetta_by_pos[pos]['fa_atr']
        row['rosetta_fa_rep'] = rosetta_by_pos[pos]['fa_rep']
        row['rosetta_hbond_total'] = sum(filter(None, [
            rosetta_by_pos[pos]['hbond_sr_bb'],
            rosetta_by_pos[pos]['hbond_lr_bb'],
            rosetta_by_pos[pos]['hbond_bb_sc'],
            rosetta_by_pos[pos]['hbond_sc'],
        ]))
    else:
        row['has_rosetta_energy'] = False
        row['rosetta_total'] = None
        row['rosetta_fa_atr'] = None
        row['rosetta_fa_rep'] = None
        row['rosetta_hbond_total'] = None

    # 5. Rosetta stabilize-PM
    if pos in stab_by_pos:
        row['has_stabilize_pm'] = True
        row['stab_best_mutation'] = stab_by_pos[pos]['best_mutation']
        row['stab_best_score'] = stab_by_pos[pos]['best_score']
        row['stab_n_stabilising'] = stab_by_pos[pos]['n_stabilising']
    else:
        row['has_stabilize_pm'] = False
        row['stab_best_mutation'] = ''
        row['stab_best_score'] = None
        row['stab_n_stabilising'] = None

    # 6. Islam et al. functional data
    if pos in islam_by_pos:
        variants = islam_by_pos[pos]
        row['has_islam'] = True
        row['n_islam_variants'] = len(variants)
        acts = [v['activity_pct'] for v in variants if v['activity_pct'] is not None]
        row['islam_mean_activity'] = sum(acts)/len(acts) if acts else None
        row['islam_functional_groups'] = '; '.join(set(v['functional_group'] for v in variants))
        row['islam_variants'] = '; '.join(v['variant'] for v in variants)
    else:
        row['has_islam'] = False
        row['n_islam_variants'] = 0
        row['islam_mean_activity'] = None
        row['islam_functional_groups'] = ''
        row['islam_variants'] = ''

    # 7. Interface data
    if pos in interface_by_pos:
        row['pcsk9_interface_dist'] = interface_by_pos[pos]['pcsk9_dist']
        row['apob_interface_dist'] = interface_by_pos[pos]['apob_dist']
        row['at_pcsk9_interface'] = interface_by_pos[pos]['at_pcsk9']
        row['at_apob_interface'] = interface_by_pos[pos]['at_apob']
    else:
        row['pcsk9_interface_dist'] = None
        row['apob_interface_dist'] = None
        row['at_pcsk9_interface'] = False
        row['at_apob_interface'] = False

    # 8. Clinical data
    clin = clinical_by_pos.get(pos, None)
    if clin and clin['n_patients'] > 0:
        row['has_clinical'] = True
        row['n_patients'] = clin['n_patients']
        row['ascvd_events'] = clin['ascvd_events']
        row['ascvd_rate'] = clin['ascvd_events'] / clin['n_patients'] * 100 if clin['n_patients'] > 0 else None
        row['mean_ldl'] = sum(clin['ldl_values'])/len(clin['ldl_values']) if clin['ldl_values'] else None
        row['mean_apob'] = sum(clin['apob_values'])/len(clin['apob_values']) if clin['apob_values'] else None
        row['mean_lpa'] = sum(clin['lpa_values'])/len(clin['lpa_values']) if clin['lpa_values'] else None
        row['mean_treatment_response'] = sum(clin['treatment_responses'])/len(clin['treatment_responses']) if clin['treatment_responses'] else None
        row['n_unique_variants'] = len(clin['variants'])
    else:
        row['has_clinical'] = False
        row['n_patients'] = 0
        row['ascvd_events'] = 0
        row['ascvd_rate'] = None
        row['mean_ldl'] = None
        row['mean_apob'] = None
        row['mean_lpa'] = None
        row['mean_treatment_response'] = None
        row['n_unique_variants'] = 0

    # Also pull from atlas
    if pos in atlas_by_pos:
        a = atlas_by_pos[pos]
        if not row['has_clinical']:
            row['n_patients'] = int(sf(a.get('n_patients', 0)) or 0)
            row['ascvd_events'] = int(sf(a.get('ascvd_events', 0)) or 0)
            row['n_unique_variants'] = int(sf(a.get('n_known_variants', 0)) or 0)
        row['atlas_risk'] = a.get('risk_classification', '')
        row['atlas_treatment'] = a.get('predicted_treatment_response', '')
        row['atlas_drug_target'] = a.get('drug_target_category', '')
        row['domain_penetrance_60'] = sf(a.get('domain_penetrance_by_60_pct'))
        row['domain_ascvd_rate'] = sf(a.get('domain_ascvd_rate'))
    else:
        row['atlas_risk'] = ''
        row['atlas_treatment'] = ''
        row['atlas_drug_target'] = ''
        row['domain_penetrance_60'] = None
        row['domain_ascvd_rate'] = None

    # ─── CONFIDENCE TIER ASSIGNMENT ───
    plddt_val = row.get('plddt')
    has_foldx = row['has_foldx']
    has_saturation = row['has_saturation']
    has_rosetta = row['has_rosetta_energy']
    has_islam = row['has_islam']
    has_stab = row['has_stabilize_pm']
    has_clinical = row['has_clinical'] or row['n_patients'] > 0

    # Count data sources
    n_sources = sum([
        has_foldx or has_saturation,  # thermodynamic data
        has_rosetta,                   # Rosetta energy
        has_islam,                     # experimental validation
        has_stab,                      # Rosetta design
        has_clinical,                  # clinical outcomes
        bool(row['at_pcsk9_interface'] or row['at_apob_interface']),  # interface
    ])
    row['n_data_sources'] = n_sources

    # Tier assignment
    if (has_foldx or has_saturation) and plddt_val and plddt_val > 70 and has_clinical:
        if has_islam:
            tier = 'HIGH+'  # Position-specific + experimentally validated
        else:
            tier = 'HIGH'
        evidence = 'Position-specific FoldX ddG + confident pLDDT + clinical data'
    elif (has_foldx or has_saturation) and plddt_val and plddt_val > 70:
        tier = 'HIGH'
        evidence = 'Position-specific FoldX ddG + confident pLDDT (no clinical data at this position)'
    elif has_rosetta and plddt_val and plddt_val > 50:
        if has_islam:
            tier = 'MODERATE+'
            evidence = 'Rosetta energy + pLDDT > 50 + Islam et al. functional data'
        elif has_clinical:
            tier = 'MODERATE+'
            evidence = 'Rosetta energy + pLDDT > 50 + clinical data'
        else:
            tier = 'MODERATE'
            evidence = 'Rosetta wildtype energy + pLDDT > 50 + domain inference'
    elif plddt_val and plddt_val > 50:
        tier = 'LOW'
        evidence = 'Domain-level inference only (pLDDT > 50)'
    elif plddt_val and plddt_val <= 50:
        tier = 'MINIMAL'
        evidence = 'Low-confidence structure (pLDDT < 50), domain inference only'
    else:
        tier = 'MINIMAL'
        evidence = 'No structural data available'

    row['confidence_tier'] = tier
    row['evidence_basis'] = evidence
    confidence_counts[tier] += 1

    rows.append(row)

# ─── Write output CSV ───
print("\nWriting confidence table...")

output_path = f"{AF}/LDLR_860_position_confidence_table.csv"
fieldnames = [
    'position', 'wildtype_aa', 'domain',
    'plddt', 'plddt_quality',
    'confidence_tier', 'n_data_sources', 'evidence_basis',
    # FoldX
    'has_foldx', 'n_foldx_variants', 'foldx_max_ddG', 'foldx_min_ddG', 'foldx_mean_ddG',
    'foldx_source', 'foldx_variants',
    # Saturation
    'has_saturation', 'n_saturation_muts', 'sat_max_ddG', 'sat_min_ddG', 'sat_mean_ddG',
    # Rosetta
    'has_rosetta_energy', 'rosetta_total', 'rosetta_fa_atr', 'rosetta_fa_rep', 'rosetta_hbond_total',
    'has_stabilize_pm', 'stab_best_mutation', 'stab_best_score', 'stab_n_stabilising',
    # Islam et al
    'has_islam', 'n_islam_variants', 'islam_mean_activity', 'islam_functional_groups', 'islam_variants',
    # Interface
    'pcsk9_interface_dist', 'apob_interface_dist', 'at_pcsk9_interface', 'at_apob_interface',
    # Clinical
    'has_clinical', 'n_patients', 'ascvd_events', 'ascvd_rate',
    'mean_ldl', 'mean_apob', 'mean_lpa', 'mean_treatment_response', 'n_unique_variants',
    # Atlas
    'atlas_risk', 'atlas_treatment', 'atlas_drug_target',
    'domain_penetrance_60', 'domain_ascvd_rate',
]

with open(output_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        # Clean up None values
        clean = {}
        for k in fieldnames:
            v = row.get(k, '')
            if v is None:
                clean[k] = ''
            elif isinstance(v, bool):
                clean[k] = 'Yes' if v else 'No'
            elif isinstance(v, float):
                clean[k] = f"{v:.4f}" if abs(v) < 100 else f"{v:.2f}"
            else:
                clean[k] = v
        writer.writerow(clean)

# ─── Print summary ───
print(f"\n{'='*70}")
print(f"  CONFIDENCE TABLE SAVED: {output_path}")
print(f"  Total positions: {len(rows)}")
print(f"{'='*70}")

print(f"\n  CONFIDENCE TIER DISTRIBUTION:")
for tier in ['HIGH+', 'HIGH', 'MODERATE+', 'MODERATE', 'LOW', 'MINIMAL']:
    n = confidence_counts.get(tier, 0)
    pct = n / len(rows) * 100
    bar = '#' * int(pct / 2)
    print(f"    {tier:12s}: {n:4d} ({pct:5.1f}%)  {bar}")

# ─── Data source coverage ───
print(f"\n  DATA SOURCE COVERAGE:")
print(f"    pLDDT available:        {sum(1 for r in rows if r.get('plddt') is not None):4d} / 860 ({sum(1 for r in rows if r.get('plddt') is not None)/860*100:.1f}%)")
print(f"    FoldX ddG (any):        {sum(1 for r in rows if r['has_foldx']):4d} / 860 ({sum(1 for r in rows if r['has_foldx'])/860*100:.1f}%)")
print(f"    Saturation mutagenesis: {sum(1 for r in rows if r['has_saturation']):4d} / 860 ({sum(1 for r in rows if r['has_saturation'])/860*100:.1f}%)")
print(f"    Rosetta energy:         {sum(1 for r in rows if r['has_rosetta_energy']):4d} / 860 ({sum(1 for r in rows if r['has_rosetta_energy'])/860*100:.1f}%)")
print(f"    Rosetta stabilize-PM:   {sum(1 for r in rows if r['has_stabilize_pm']):4d} / 860 ({sum(1 for r in rows if r['has_stabilize_pm'])/860*100:.1f}%)")
print(f"    Islam et al. functional:{sum(1 for r in rows if r['has_islam']):4d} / 860 ({sum(1 for r in rows if r['has_islam'])/860*100:.1f}%)")
print(f"    Clinical patient data:  {sum(1 for r in rows if r['n_patients'] > 0):4d} / 860 ({sum(1 for r in rows if r['n_patients'] > 0)/860*100:.1f}%)")
print(f"    At PCSK9 interface:     {sum(1 for r in rows if r['at_pcsk9_interface']):4d} / 860")
print(f"    At ApoB interface:      {sum(1 for r in rows if r['at_apob_interface']):4d} / 860")

# ─── Identify positions that NEED AlphaFold3 mutant runs ───
print(f"\n{'='*70}")
print(f"  POSITIONS REQUIRING AlphaFold3 MUTANT PREDICTIONS")
print(f"{'='*70}")

# Priority 1: HIGH clinical impact + FoldX data (validate backbone assumption)
priority1 = []
for r in rows:
    if r['has_foldx'] and r.get('foldx_max_ddG') and r['foldx_max_ddG'] > 4.0:
        priority1.append(r)
priority1.sort(key=lambda x: x.get('foldx_max_ddG', 0), reverse=True)

print(f"\n  PRIORITY 1 — Highly destabilising (ddG > 4 kcal/mol) — validate backbone assumption:")
for r in priority1[:15]:
    print(f"    Pos {r['position']:4d} ({r['wildtype_aa']}) | {r['domain']:30s} | ddG={r['foldx_max_ddG']:+6.2f} | pLDDT={r.get('plddt', 'N/A')} | patients={r['n_patients']}")

# Priority 2: Islam et al. defective variants (experimental ground truth)
priority2 = []
for r in rows:
    if r['has_islam'] and 'Defective' in r.get('islam_functional_groups', ''):
        priority2.append(r)
priority2.sort(key=lambda x: x.get('islam_mean_activity', 100))

print(f"\n  PRIORITY 2 — Islam et al. 'Defective' positions (experimental ground truth):")
for r in priority2[:10]:
    print(f"    Pos {r['position']:4d} ({r['wildtype_aa']}) | {r['domain']:30s} | activity={r.get('islam_mean_activity', 'N/A'):.1f}% | groups={r['islam_functional_groups']}")

# Priority 3: Interface residues (validate complex predictions)
priority3 = []
for r in rows:
    if r['at_pcsk9_interface'] or r['at_apob_interface']:
        priority3.append(r)

print(f"\n  PRIORITY 3 — Interface residues (validate complex predictions):")
for r in priority3[:10]:
    iface = []
    if r['at_pcsk9_interface']:
        iface.append('PCSK9')
    if r['at_apob_interface']:
        iface.append('ApoB')
    print(f"    Pos {r['position']:4d} ({r['wildtype_aa']}) | {r['domain']:30s} | Interface: {', '.join(iface)} | patients={r['n_patients']}")

# Priority 4: One representative per domain (domain validation)
priority4 = []
domains_covered = set()
for r in sorted(rows, key=lambda x: x['n_patients'], reverse=True):
    d = r.get('domain', '')
    if d and d not in domains_covered and r['n_patients'] > 0:
        domains_covered.add(d)
        priority4.append(r)

print(f"\n  PRIORITY 4 — Domain representatives (one per domain, highest patient count):")
for r in priority4:
    print(f"    Pos {r['position']:4d} ({r['wildtype_aa']}) | {r['domain']:30s} | patients={r['n_patients']} | ddG={r.get('foldx_max_ddG', 'N/A')}")

# Priority 5: Drug target positions (Rosetta stabilize-PM)
priority5 = [r for r in rows if r['has_stabilize_pm']]
priority5.sort(key=lambda x: x.get('stab_best_score') or 0)

print(f"\n  PRIORITY 5 — Rosetta stabilize-PM drug targets:")
for r in priority5:
    print(f"    Pos {r['position']:4d} ({r['wildtype_aa']}) | {r['domain']:30s} | stab={r.get('stab_best_mutation', '')} score={r.get('stab_best_score', 'N/A')}")

# ─── Compile master list of positions for AF3 mutant runs ───
af3_positions = set()
af3_details = {}

# Add all priority categories
for r in priority1[:15]:
    af3_positions.add(r['position'])
    af3_details[r['position']] = f"P1-destabilising(ddG={r['foldx_max_ddG']:+.2f})"

for r in priority2[:10]:
    if r['position'] not in af3_positions:
        af3_positions.add(r['position'])
        af3_details[r['position']] = f"P2-islam-defective(act={r.get('islam_mean_activity', 0):.1f}%)"

for r in priority3[:10]:
    if r['position'] not in af3_positions:
        af3_positions.add(r['position'])
        iface = 'PCSK9' if r['at_pcsk9_interface'] else 'ApoB'
        af3_details[r['position']] = f"P3-interface({iface})"

for r in priority4:
    if r['position'] not in af3_positions:
        af3_positions.add(r['position'])
        af3_details[r['position']] = f"P4-domain-rep({r['domain']})"

for r in priority5:
    if r['position'] not in af3_positions:
        af3_positions.add(r['position'])
        af3_details[r['position']] = f"P5-drug-target({r.get('stab_best_mutation', '')})"

print(f"\n{'='*70}")
print(f"  TOTAL UNIQUE POSITIONS FOR AF3 MUTANT RUNS: {len(af3_positions)}")
print(f"{'='*70}")

# Save AF3 target list
af3_list_path = f"{AF}/alphafold3_mutant_targets.csv"
with open(af3_list_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(['position', 'wildtype_aa', 'domain', 'plddt', 'priority_reason',
                     'foldx_max_ddG', 'rosetta_total', 'islam_activity', 'n_patients',
                     'suggested_mutation'])
    for pos in sorted(af3_positions):
        r = rows[pos - 1]  # 0-indexed
        # Determine the most clinically relevant mutation at this position
        suggested = ''
        if r['has_foldx'] and r.get('foldx_variants'):
            # Use the most destabilising known variant
            suggested = r['foldx_variants'].split(';')[0].strip().split('=')[0] if r['foldx_variants'] else ''
        elif r['has_islam'] and r.get('islam_variants'):
            suggested = r['islam_variants'].split(';')[0].strip()

        writer.writerow([
            pos, r['wildtype_aa'], r['domain'],
            f"{r['plddt']:.1f}" if r.get('plddt') else '',
            af3_details.get(pos, ''),
            f"{r['foldx_max_ddG']:.2f}" if r.get('foldx_max_ddG') else '',
            f"{r['rosetta_total']:.2f}" if r.get('rosetta_total') else '',
            f"{r['islam_mean_activity']:.1f}" if r.get('islam_mean_activity') else '',
            r['n_patients'],
            suggested
        ])

print(f"  AF3 target list saved: {af3_list_path}")

print(f"\n  Done.")

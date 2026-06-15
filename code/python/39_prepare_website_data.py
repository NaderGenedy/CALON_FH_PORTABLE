"""
CALON-FH Website Data Preparation
Converts 12 analysis CSVs into a single JSON for the interactive clinical atlas website.
"""
import csv, json, os

DATA_DIR = r'C:\Users\nader\Downloads\calon_ukb_pipeline\alphafold\analysis'
OUT_DIR = r'C:\Users\nader\Downloads\calon_ukb_pipeline\website'
os.makedirs(OUT_DIR, exist_ok=True)

def read_csv(filename):
    path = os.path.join(DATA_DIR, filename)
    with open(path, 'r', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

def safe_float(val, default=None):
    try:
        return round(float(val), 4) if val and val.strip() else default
    except (ValueError, TypeError):
        return default

def safe_int(val, default=None):
    try:
        return int(float(val)) if val and val.strip() else default
    except (ValueError, TypeError):
        return default

print("Loading CSVs...")

# 1. Atlas (860 positions)
atlas_raw = read_csv('LDLR_genotype_phenotype_atlas.csv')
atlas = {}
for r in atlas_raw:
    pos = safe_int(r.get('position'))
    if pos is None:
        continue
    atlas[pos] = {
        'wt_aa': r.get('wildtype_aa', ''),
        'domain': r.get('domain', ''),
        'plddt': safe_float(r.get('plddt_confidence')),
        'pcsk9_dist': safe_float(r.get('pcsk9_interface_distance')),
        'apob_dist': safe_float(r.get('apob_interface_distance')),
        'at_pcsk9': r.get('at_pcsk9_interface', '') == 'True',
        'at_apob': r.get('at_apob_interface', '') == 'True',
        'penetrance_60': safe_float(r.get('domain_penetrance_by_60_pct')),
        'ascvd_rate': safe_float(r.get('domain_ascvd_rate')),
        'max_ddg': safe_float(r.get('observed_max_ddG')),
        'max_sss': safe_float(r.get('observed_max_sss')),
        'n_variants': safe_int(r.get('n_known_variants'), 0),
        'n_patients': safe_int(r.get('n_patients'), 0),
        'ascvd_events': safe_int(r.get('ascvd_events'), 0),
        'risk': r.get('risk_classification', ''),
        'treatment': r.get('predicted_treatment_response', ''),
        'drug_target': r.get('drug_target_category', ''),
    }
print(f"  Atlas: {len(atlas)} positions")

# 2. Confidence table (860 positions)
conf_raw = read_csv('LDLR_860_position_confidence_table.csv')
confidence = {}
for r in conf_raw:
    pos = safe_int(r.get('position'))
    if pos is None:
        continue
    confidence[pos] = {
        'tier': r.get('confidence_tier', ''),
        'tier_v2': r.get('confidence_tier_v2', ''),
        'n_sources': safe_int(r.get('n_data_sources'), 0),
        'evidence': r.get('evidence_basis', ''),
        'has_foldx': r.get('has_foldx', '') == 'True',
        'foldx_max': safe_float(r.get('foldx_max_ddG')),
        'foldx_mean': safe_float(r.get('foldx_mean_ddG')),
        'has_islam': r.get('has_islam', '') == 'True',
        'islam_activity': safe_float(r.get('islam_mean_activity')),
        'n_islam': safe_int(r.get('n_islam_variants'), 0),
        'has_clinical': r.get('has_clinical', '') == 'True',
        'sat_sensitivity': safe_float(r.get('mutational_sensitivity')),
        'sat_pct_destab': safe_float(r.get('sat_pct_destabilising')),
    }
print(f"  Confidence: {len(confidence)} positions")

# 3. pLDDT per residue
plddt_raw = read_csv('LDLR_wildtype_plddt_per_residue.csv')
plddt_track = []
for r in plddt_raw:
    plddt_track.append({
        'pos': safe_int(r.get('residue_num')),
        'plddt': safe_float(r.get('plddt')),
        'domain': r.get('domain', ''),
        'quality': r.get('quality_category', ''),
    })
print(f"  pLDDT track: {len(plddt_track)} residues")

# 4. AF3 mutant validation
af3_raw = read_csv('af3_mutant_validation.csv')
af3_mutants = {}
for r in af3_raw:
    pos = safe_int(r.get('position'))
    if pos:
        af3_mutants[r.get('mutation', '')] = {
            'pos': pos,
            'wt': r.get('wt_aa', ''),
            'mut': r.get('mut_aa', ''),
            'ptm': safe_float(r.get('ptm')),
            'delta_ptm': safe_float(r.get('delta_ptm')),
        }
print(f"  AF3 monomers: {len(af3_mutants)}")

# 5. AF3 complex validation
af3c_raw = read_csv('af3_complex_validation.csv')
af3_complex = {}
for r in af3c_raw:
    mut = r.get('mutation', '')
    if mut:
        af3_complex[mut] = {
            'pos': safe_int(r.get('position')),
            'ptm': safe_float(r.get('complex_ptm')),
            'iptm': safe_float(r.get('complex_iptm')),
            'delta_ptm': safe_float(r.get('delta_complex_ptm')),
            'delta_iptm': safe_float(r.get('delta_complex_iptm')),
        }
print(f"  AF3 complexes: {len(af3_complex)}")

# 6. FoldX clinical variants
foldx_raw = read_csv('foldx_ddg_results.csv')
foldx_variants = []
for r in foldx_raw:
    foldx_variants.append({
        'variant': r.get('variant_id', ''),
        'pos': safe_int(r.get('position')),
        'wt': r.get('wt_aa', ''),
        'mut': r.get('mut_aa', ''),
        'domain': r.get('domain', ''),
        'n_patients': safe_int(r.get('n_patients'), 0),
        'ddg': safe_float(r.get('ddG_kcal_mol')),
        'effect': r.get('effect', ''),
    })
print(f"  FoldX variants: {len(foldx_variants)}")

# 7. FoldX saturation (summarise per position)
sat_raw = read_csv('foldx_saturation_mutagenesis.csv')
saturation = {}
for r in sat_raw:
    pos = safe_int(r.get('position'))
    ddg = safe_float(r.get('ddG'))
    if pos is None or ddg is None:
        continue
    if pos not in saturation:
        saturation[pos] = {'values': [], 'mutations': []}
    saturation[pos]['values'].append(ddg)
    saturation[pos]['mutations'].append({
        'from': r.get('wt_aa', ''),
        'to': r.get('mut_aa', ''),
        'ddg': ddg,
        'effect': r.get('effect_class', ''),
    })

# Summarise
sat_summary = {}
for pos, data in saturation.items():
    vals = data['values']
    destab = [v for v in vals if v > 2.0]
    sat_summary[pos] = {
        'n': len(vals),
        'mean': round(sum(vals)/len(vals), 2),
        'max': round(max(vals), 2),
        'min': round(min(vals), 2),
        'n_destab': len(destab),
        'pct_destab': round(100 * len(destab) / len(vals), 1),
        'worst': data['mutations'][vals.index(max(vals))],
        'best': data['mutations'][vals.index(min(vals))],
    }
print(f"  Saturation: {len(sat_summary)} positions")

# 8. Domain penetrance (KM results)
km_raw = read_csv('cox_km_compound_results.csv')
domain_penetrance = []
for r in km_raw:
    if r.get('analysis', '').startswith('KM_domain'):
        domain_penetrance.append({
            'domain': r.get('group', ''),
            'n': safe_int(r.get('n'), 0),
            'events': safe_int(r.get('events'), 0),
            'rate': safe_float(r.get('rate')),
            'penetrance_60': safe_float(r.get('value')),
        })
# Sort by penetrance
domain_penetrance.sort(key=lambda x: x['penetrance_60'] or 0, reverse=True)
print(f"  Domain penetrance: {len(domain_penetrance)} domains")

# 9. AUC forest plot
auc_raw = read_csv('sss_v2_forest_plot_data.csv')
forest_plot = []
for r in auc_raw:
    forest_plot.append({
        'cohort': r.get('cohort', ''),
        'model': r.get('subgroup', ''),
        'n': safe_int(r.get('N'), 0),
        'events': safe_int(r.get('events'), 0),
        'auc': safe_float(r.get('AUC')),
        'ci_lo': safe_float(r.get('CI_lower')),
        'ci_hi': safe_float(r.get('CI_upper')),
    })
print(f"  Forest plot: {len(forest_plot)} models")

# 10. Treatment response
tx_raw = read_csv('external_treatment_response_validation.csv')
treatment_response = []
for r in tx_raw:
    n_un = safe_int(r.get('n_untreated'), 0)
    n_tr = safe_int(r.get('n_treated'), 0)
    if n_un + n_tr < 3:
        continue
    treatment_response.append({
        'domain': r.get('domain', ''),
        'n_untreated': n_un,
        'n_treated': n_tr,
        'ldl_untreated': safe_float(r.get('mean_ldl_untreated')),
        'ldl_treated': safe_float(r.get('mean_ldl_treated')),
        'ldl_diff': safe_float(r.get('ldl_diff')),
        'p_value': safe_float(r.get('p_value')),
    })
print(f"  Treatment response: {len(treatment_response)} domains")

# 11. Islam functional data
islam_raw = read_csv('islam_et_al_315_variants.csv')
islam = {}
for r in islam_raw:
    variant = r.get('variant', '')
    if variant:
        islam[variant] = {
            'pos': safe_int(r.get('position')),
            'domain': r.get('domain', ''),
            'activity_pct': safe_float(r.get('activity_pct')),
            'group': r.get('functional_group', ''),
            'clinvar': r.get('clinvar', ''),
        }
print(f"  Islam variants: {len(islam)}")

# 13. c. notation mapping from SSS catalogue and FoldX results
# Build lookup: "pos_wt_mut" -> c. notation string
c_notation_map = {}

# From FoldX results (most reliable - has variant_id, position, wt_aa, mut_aa)
for r in foldx_raw:
    vid = r.get('variant_id', '')
    pos = safe_int(r.get('position'))
    wt = r.get('wt_aa', '')
    mut = r.get('mut_aa', '')
    if vid.startswith('LDLR:c.') and pos and wt and mut:
        c_not = vid.replace('LDLR:', '')  # e.g. "c.1055G>A"
        key = f"{pos}_{wt}_{mut}"
        c_notation_map[key] = c_not

# Also from SSS catalogue (wider coverage)
sss_raw = read_csv('sss_v2_catalogue.csv')
for r in sss_raw:
    vid = r.get('variant_id', '')
    pos = safe_int(r.get('position'))
    gene = r.get('gene', '')
    vtype = r.get('variant_type', '')
    if gene == 'LDLR' and vid.startswith('LDLR:c.') and pos and vtype == 'substitution':
        c_not = vid.replace('LDLR:', '')  # e.g. "c.1055G>A"
        # We don't have wt_aa/mut_aa in catalogue, but we can store by position
        # Multiple variants at same position will be stored as list
        pos_key = str(pos)
        if pos_key not in c_notation_map:
            # Store position-only entries for fallback
            pass
        # Store the c. notation keyed by variant_id for direct lookup
        c_notation_map[vid] = c_not

# Build a position → list of c. notations lookup for position-only searches
c_notation_by_pos = {}
for r in foldx_raw:
    vid = r.get('variant_id', '')
    pos = safe_int(r.get('position'))
    wt = r.get('wt_aa', '')
    mut = r.get('mut_aa', '')
    if vid.startswith('LDLR:c.') and pos:
        c_not = vid.replace('LDLR:', '')
        pos_key = str(pos)
        if pos_key not in c_notation_by_pos:
            c_notation_by_pos[pos_key] = []
        entry = {'c': c_not, 'wt': wt, 'mut': mut}
        if entry not in c_notation_by_pos[pos_key]:
            c_notation_by_pos[pos_key].append(entry)

# Add from SSS catalogue too
for r in sss_raw:
    vid = r.get('variant_id', '')
    pos = safe_int(r.get('position'))
    gene = r.get('gene', '')
    vtype = r.get('variant_type', '')
    if gene == 'LDLR' and vid.startswith('LDLR:c.') and pos and vtype == 'substitution':
        c_not = vid.replace('LDLR:', '')
        pos_key = str(pos)
        if pos_key not in c_notation_by_pos:
            c_notation_by_pos[pos_key] = []
        # Check if this c. notation is already there
        existing = [e['c'] for e in c_notation_by_pos[pos_key]]
        if c_not not in existing:
            c_notation_by_pos[pos_key].append({'c': c_not, 'wt': '', 'mut': ''})

print(f"  c. notation: {len(c_notation_map)} variant mappings, {len(c_notation_by_pos)} positions")

# 14. UKB ClinVar + VEP (SIFT/PolyPhen) — keyed by "pos_wt_mut"
vep_raw = read_csv('ukb_ldlr_vep_annotations.csv')
vep_data = {}
for r in vep_raw:
    ppos = safe_int(r.get('protein_position'))
    aas = r.get('amino_acids', '')
    if ppos and '/' in aas:
        wt, mut = aas.split('/')
        key = f"{ppos}_{wt}_{mut}"
        vep_data[key] = {
            'sift': r.get('sift', ''),
            'polyphen': r.get('polyphen', ''),
            'consequence': r.get('consequence', ''),
            'impact': r.get('impact', ''),
        }
print(f"  VEP annotations: {len(vep_data)} variants")

# UKB ClinVar — need to map genomic pos to protein pos via VEP
ukb_clinvar_raw = read_csv('ukb_clinvar_full.csv')
ukb_clinvar = {}
for r in ukb_clinvar_raw:
    cv = r.get('clinvar', '').strip()
    gpos = r.get('pos', '')
    ref = r.get('ref', '')
    alt = r.get('alt', '')
    if cv:
        # Match to VEP by genomic coordinates
        for vr in vep_raw:
            if vr.get('pos') == gpos and vr.get('ref') == ref and vr.get('alt') == alt:
                ppos = safe_int(vr.get('protein_position'))
                aas = vr.get('amino_acids', '')
                if ppos and '/' in aas:
                    wt, mut = aas.split('/')
                    key = f"{ppos}_{wt}_{mut}"
                    ukb_clinvar[key] = cv
                break
print(f"  UKB ClinVar: {len(ukb_clinvar)} variants mapped")

# 15. Tabet DMS scores — keyed by "pos_wt_mut" (single letter)
tabet_path = os.path.join(r'C:\Users\nader\Downloads\calon_ukb_pipeline\alphafold\tabet_data', 'science.ady7186_data_s1.csv')
tabet_data = {}
AA3TO1 = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G','His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'}
if os.path.exists(tabet_path):
    with open(tabet_path, 'r', encoding='utf-8-sig') as f:
        tabet_raw = list(csv.DictReader(f))
    for r in tabet_raw:
        hgvsp = r.get('hgvsp', '')
        pos = safe_int(r.get('aapos'))
        score = safe_float(r.get('score'))
        if pos and score is not None and len(hgvsp) >= 6:
            # Parse e.g. "Gly2Leu" -> wt=G, mut=L
            import re as _re
            m = _re.match(r'^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', hgvsp)
            if m:
                wt1 = AA3TO1.get(m.group(1), '')
                mut1 = AA3TO1.get(m.group(3), '')
                if wt1 and mut1:
                    key = f"{pos}_{wt1}_{mut1}"
                    tabet_data[key] = score
    print(f"  Tabet DMS: {len(tabet_data)} variant scores")
else:
    print("  Tabet DMS: file not found, skipping")

# Build combined evidence lookup for classification
# Keys: "pos_wt_mut" -> dict of evidence sources
evidence_lookup = {}
# Collect all unique variant keys
all_keys = set()
for key in c_notation_map:
    if '_' in key:
        all_keys.add(key)
for key in vep_data:
    all_keys.add(key)
for key in ukb_clinvar:
    all_keys.add(key)
for key in tabet_data:
    all_keys.add(key)
# Islam variants (keyed as "C352Y" -> need to convert)
for var, info in islam.items():
    if info['pos'] and len(var) >= 3:
        wt = var[0]
        mut = var[-1]
        key = f"{info['pos']}_{wt}_{mut}"
        all_keys.add(key)

for key in all_keys:
    parts = key.split('_')
    if len(parts) != 3:
        continue
    pos, wt, mut = parts
    ev = {}

    # ClinVar (Islam source)
    variant_str = f"{wt}{pos}{mut}"
    if variant_str in islam:
        cv = islam[variant_str].get('clinvar', '')
        if cv:
            ev['clinvar_islam'] = cv

    # ClinVar (UKB source)
    if key in ukb_clinvar:
        ev['clinvar_ukb'] = ukb_clinvar[key]

    # SIFT / PolyPhen
    if key in vep_data:
        ev['sift'] = vep_data[key]['sift']
        ev['polyphen'] = vep_data[key]['polyphen']

    # Tabet DMS
    if key in tabet_data:
        ev['tabet_score'] = tabet_data[key]

    if ev:
        evidence_lookup[key] = ev

print(f"  Evidence lookup: {len(evidence_lookup)} variants with external evidence")

# 12. Rosetta energy (top-level per residue)
rosetta_raw = read_csv('rosetta_energy_breakdown.csv')
rosetta = []
for r in rosetta_raw:
    rosetta.append({
        'pos': safe_int(r.get('label', '').replace('residue_', '')),
        'total': safe_float(r.get('total')),
        'hbond': safe_float(r.get('hbond_sr_bb')),
        'fa_rep': safe_float(r.get('fa_rep')),
    })
print(f"  Rosetta: {len(rosetta)} residues")

# Build domain colour map
DOMAIN_COLOURS = {
    'Signal peptide': '#999999',
    'Ligand-binding': '#E41A1C',
    'Ligand-binding R1': '#E41A1C',
    'Ligand-binding R2': '#FF6B6B',
    'Ligand-binding R3': '#E41A1C',
    'Ligand-binding R4': '#FF6B6B',
    'Ligand-binding R5': '#E41A1C',
    'Ligand-binding R6': '#FF6B6B',
    'Ligand-binding R7': '#E41A1C',
    'EGF-like A': '#377EB8',
    'EGF-like B': '#4DBEEE',
    'EGF-like C': '#377EB8',
    'Beta-propeller': '#4DAF4A',
    'EGF-precursor-homology': '#984EA3',
    'EGF_homology': '#984EA3',
    'O-linked sugar': '#FF7F00',
    'Linker (EGF-C/O-linked)': '#A65628',
    'Transmembrane': '#F781BF',
    'Cytoplasmic': '#BDBDBD',
}

# Domain boundaries for the linear diagram
DOMAIN_REGIONS = [
    {'name': 'Signal', 'start': 1, 'end': 21, 'colour': '#999999'},
    {'name': 'R1', 'start': 22, 'end': 64, 'colour': '#E41A1C'},
    {'name': 'R2', 'start': 65, 'end': 107, 'colour': '#FF6B6B'},
    {'name': 'R3', 'start': 108, 'end': 148, 'colour': '#E41A1C'},
    {'name': 'R4', 'start': 149, 'end': 190, 'colour': '#FF6B6B'},
    {'name': 'R5', 'start': 191, 'end': 232, 'colour': '#E41A1C'},
    {'name': 'R6', 'start': 233, 'end': 274, 'colour': '#FF6B6B'},
    {'name': 'R7', 'start': 275, 'end': 314, 'colour': '#E41A1C'},
    {'name': 'EGF-A', 'start': 315, 'end': 353, 'colour': '#377EB8'},
    {'name': 'EGF-B', 'start': 354, 'end': 393, 'colour': '#4DBEEE'},
    {'name': 'β-prop', 'start': 394, 'end': 582, 'colour': '#4DAF4A'},
    {'name': 'EGF-C', 'start': 583, 'end': 626, 'colour': '#377EB8'},
    {'name': 'O-linked', 'start': 627, 'end': 694, 'colour': '#FF7F00'},
    {'name': 'TM', 'start': 695, 'end': 767, 'colour': '#F781BF'},
    {'name': 'Cyto', 'start': 768, 'end': 860, 'colour': '#BDBDBD'},
]

# Assemble final JSON
website_data = {
    'atlas': atlas,
    'confidence': confidence,
    'plddt_track': plddt_track,
    'af3_mutants': af3_mutants,
    'af3_complex': af3_complex,
    'foldx_variants': foldx_variants,
    'saturation_summary': {str(k): v for k, v in sat_summary.items()},
    'domain_penetrance': domain_penetrance,
    'forest_plot': forest_plot,
    'treatment_response': treatment_response,
    'islam': islam,
    'c_notation': c_notation_map,
    'c_notation_by_pos': c_notation_by_pos,
    'evidence': evidence_lookup,
    'domain_colours': DOMAIN_COLOURS,
    'domain_regions': DOMAIN_REGIONS,
    'stats': {
        'n_positions': len(atlas),
        'n_af3_monomers': len(af3_mutants),
        'n_af3_complexes': len(af3_complex),
        'n_foldx_variants': len(foldx_variants),
        'n_saturation_positions': len(sat_summary),
        'n_islam_variants': len(islam),
        'n_domains_penetrance': len(domain_penetrance),
    },
}

out_path = os.path.join(OUT_DIR, 'website_data.json')
with open(out_path, 'w') as f:
    json.dump(website_data, f, separators=(',', ':'))

file_size = os.path.getsize(out_path) / 1024
print(f"\nSaved: {out_path}")
print(f"Size: {file_size:.0f} KB")
print("DONE")

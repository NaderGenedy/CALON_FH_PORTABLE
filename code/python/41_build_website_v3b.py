#!/usr/bin/env python3
"""
41_build_website_v3b.py — CALON-FH Clinical Atlas Website v3b
Integrates: ClinVar classifications, variant-level LLT response (DRAGON3/Wales),
Islam et al. functional validation, traffic light, penetrance, 3D Mol* viewer.
"""

import pandas as pd
import json
import math
import numpy as np
import os

BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")

# ── Load all data ──────────────────────────────────────────────────────
print("[1/5] Loading data...")
pos = pd.read_csv(os.path.join(ANALYSIS, "LDLR_860_position_confidence_table_v3b.csv"))
var = pd.read_csv(os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas_v3b.csv"))
treatment_df = pd.read_csv(os.path.join(ANALYSIS, "external_treatment_response_validation.csv"))
km_df = pd.read_csv(os.path.join(ANALYSIS, "cox_km_compound_results.csv"))
vus_df = pd.read_csv(os.path.join(ANALYSIS, "vus_reclassification_results.csv"))
ext_vus = pd.read_csv(os.path.join(ANALYSIS, "external_vus_validation.csv"))
llt_df = pd.read_csv(os.path.join(ANALYSIS, "variant_llt_response.csv"))
islam_df = pd.read_csv(os.path.join(ANALYSIS, "islam_et_al_315_variants.csv"))
clinvar_df = pd.read_csv(os.path.join(ANALYSIS, "ukb_clinvar_full.csv"))

AA1TO3 = {'A':'Ala','C':'Cys','D':'Asp','E':'Glu','F':'Phe','G':'Gly','H':'His','I':'Ile',
           'K':'Lys','L':'Leu','M':'Met','N':'Asn','P':'Pro','Q':'Gln','R':'Arg',
           'S':'Ser','T':'Thr','V':'Val','W':'Trp','Y':'Tyr'}

def safe(v):
    if pd.isna(v): return None
    if isinstance(v, (np.floating, float)):
        if math.isnan(v) or math.isinf(v): return None
        return round(float(v), 3)
    if isinstance(v, (np.integer, int)):
        return int(v)
    return str(v)

# ── NEW: AlphaGenome + Tabet DMS data ────────────────────────────────
AG_DATA = r"D:/Projects/CALON_AlphaFold_Rebuild/data"
TABET_DATA = r"D:/Projects/CALON_AlphaFold_Rebuild/data/tabet_2025"

# AlphaGenome scores (224 positions)
ag_path = os.path.join(AG_DATA, "alphagenome_sss_359_scores.csv")
if os.path.exists(ag_path):
    ag_df = pd.read_csv(ag_path)
    ag_df['protein_pos'] = pd.to_numeric(ag_df['protein_pos'], errors='coerce')
    ag_by_pos = {}
    for _, ar in ag_df.iterrows():
        p = ar['protein_pos']
        if pd.notna(p):
            ag_by_pos[int(p)] = {
                'cage_diff': safe(ar.get('cage_diff_mean')),
                'splice_max': safe(ar.get('splice_junction_max')),
                'atac_diff': safe(ar.get('atac_diff_mean')),
                'l2_impact': safe(ar.get('l2_impact')),
            }
    print(f"   AlphaGenome: {len(ag_by_pos)} positions with regulatory scores")
else:
    ag_by_pos = {}
    print("   AlphaGenome: not available")

# AlphaGenome promoter data
promoter_path = os.path.join(AG_DATA, "ldlr_promoter_importance.csv")
if os.path.exists(promoter_path):
    promoter_df = pd.read_csv(promoter_path)
    print(f"   Promoter ISM: {len(promoter_df)} positions")
else:
    promoter_df = pd.DataFrame()

# Tabet DMS scores (17,390 variants)
tabet_path = os.path.join(TABET_DATA, "science.ady7186_data_s1.csv")
if os.path.exists(tabet_path):
    tabet_df = pd.read_csv(tabet_path)
    # Aggregate by position: mean DMS, % damaging, best/worst substitution
    tabet_df['aapos'] = pd.to_numeric(tabet_df['aapos'], errors='coerce')
    tabet_by_pos = {}
    for pos_val, grp in tabet_df.groupby('aapos'):
        if pd.isna(pos_val):
            continue
        scores = grp['score'].dropna()
        n_damaging = (scores < 0.118).sum()
        tabet_by_pos[int(pos_val)] = {
            'dms_mean': safe(scores.mean()),
            'dms_min': safe(scores.min()),
            'dms_max': safe(scores.max()),
            'dms_n': int(len(scores)),
            'dms_pct_damaging': safe(n_damaging / len(scores) * 100) if len(scores) > 0 else None,
        }
    print(f"   Tabet DMS: {len(tabet_by_pos)} positions with functional scores")

    # Also load individual variant DMS for variant-level display
    import re
    from Bio.Data.IUPACData import protein_letters_3to1
    aa3to1 = {k.capitalize(): v for k, v in protein_letters_3to1.items()}
    aa3to1['*'] = '*'

    def parse_hgvsp_tabet(hgvsp):
        if pd.isna(hgvsp) or '=' in str(hgvsp):
            return None, None, None
        m = re.match(r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)', str(hgvsp))
        if m:
            wt3, p, mut3 = m.groups()
            wt = aa3to1.get(wt3, '?')
            mut = aa3to1.get(mut3, '?') if mut3 != '*' else '*'
            return int(p), wt, mut
        return None, None, None

    tabet_df[['t_pos','t_wt','t_mut']] = pd.DataFrame(
        tabet_df['hgvsp'].apply(parse_hgvsp_tabet).tolist(), index=tabet_df.index
    )
    tabet_variant_map = {}
    for _, tr in tabet_df.dropna(subset=['t_pos']).iterrows():
        key = f"{int(tr['t_pos'])}_{tr['t_wt']}_{tr['t_mut']}"
        tabet_variant_map[key] = safe(tr['score'])
    print(f"   Tabet variant-level DMS: {len(tabet_variant_map)} individual scores")
else:
    tabet_by_pos = {}
    tabet_variant_map = {}
    print("   Tabet DMS: not available")

# Tabet VCEP classifications
tabet_vcep_path = os.path.join(TABET_DATA, "science.ady7186_data_s4.csv")
if os.path.exists(tabet_vcep_path):
    tabet_vcep = pd.read_csv(tabet_vcep_path)
    tabet_vcep_map = dict(zip(tabet_vcep['hgvs_pro'], tabet_vcep['VCEP classification']))
    print(f"   Tabet VCEP: {len(tabet_vcep_map)} expert-classified variants")
else:
    tabet_vcep_map = {}

# Enhanced therapeutic map with eSSS
esss_path = os.path.join(AG_DATA, "enhanced_therapeutic_map_16340.csv")
if os.path.exists(esss_path):
    esss_df = pd.read_csv(esss_path, usecols=['position','wt_aa','mut_aa','eSSS','eSSS_class','dms_uptake','dms_abundance'])
    esss_by_pos = {}
    for pos_val, grp in esss_df.groupby('position'):
        esss_vals = grp['eSSS'].dropna()
        esss_by_pos[int(pos_val)] = {
            'esss_mean': safe(esss_vals.mean()),
            'esss_max': safe(esss_vals.max()),
            'n_severe': int((grp['eSSS_class'] == 'Severe').sum()),
            'n_moderate': int((grp['eSSS_class'] == 'Moderate').sum()),
        }
    print(f"   eSSS: {len(esss_by_pos)} positions with enhanced scores")
else:
    esss_by_pos = {}

# ── Process ClinVar classifications ───────────────────────────────────
print("[2/5] Processing ClinVar + Islam + LLT data...")

# 1. ClinVar from external_vus_validation.csv (variant_id → classification)
vus_class_map = ext_vus.groupby('variant_id')['classification'].first().to_dict()
print(f"   ClinVar from ext_vus: {len(vus_class_map)} variants")

# 2. ClinVar from ukb_clinvar_full.csv — simplify multi-label classifications
def simplify_clinvar(cv):
    if pd.isna(cv): return None
    cv = str(cv).lower()
    if 'pathogenic' in cv and 'uncertain' not in cv and 'benign' not in cv:
        return 'Pathogenic'
    if 'likely_pathogenic' in cv and 'benign' not in cv:
        return 'Likely pathogenic'
    if 'benign' in cv and 'pathogenic' not in cv:
        return 'Benign'
    if 'likely_benign' in cv and 'pathogenic' not in cv:
        return 'Likely benign'
    if 'uncertain' in cv:
        return 'VUS'
    return None

clinvar_df['clinvar_simple'] = clinvar_df['clinvar'].apply(simplify_clinvar)
# Map by genomic position (chrom:pos:ref:alt)
clinvar_by_pos = {}
for _, cr in clinvar_df.iterrows():
    if pd.notna(cr.get('clinvar_simple')):
        key = f"{cr['chrom']}:{cr['pos']}:{cr['ref']}:{cr['alt']}"
        clinvar_by_pos[key] = cr['clinvar_simple']
print(f"   ClinVar from UKB full: {len(clinvar_by_pos)} variants with classification")

# 3. Islam et al. ClinVar + functional groups (by variant name → position)
islam_by_pos = {}
for _, ir in islam_df.iterrows():
    p = int(ir['position']) if pd.notna(ir.get('position')) else None
    if p:
        if p not in islam_by_pos:
            islam_by_pos[p] = []
        islam_by_pos[p].append({
            'variant': str(ir['variant']),
            'clinvar': str(ir['clinvar']) if pd.notna(ir.get('clinvar')) else None,
            'functional_group': str(ir['functional_group']) if pd.notna(ir.get('functional_group')) else None,
            'activity_pct': safe(ir.get('activity_pct')),
        })
print(f"   Islam et al: {len(islam_by_pos)} positions with functional data")

# 4. Variant-level LLT response from OUR work (DRAGON3/Wales FH)
llt_treated = llt_df[llt_df['intensity'] != 'none'].copy()
llt_by_variant = {}
for vid, grp in llt_treated.groupby('variant'):
    llt_by_variant[vid] = {
        'n_treated': len(grp),
        'ldl_before': round(float(grp['ldl_1'].mean()), 2),
        'ldl_after': round(float(grp['ldl_2'].mean()), 2) if grp['ldl_2'].notna().any() else None,
        'ldl_pct_reduction': round(float(grp['ldl_pct_reduction'].mean()), 1),
        'intensity': grp['intensity'].mode().iloc[0] if len(grp) > 0 else None,
    }
print(f"   LLT variant-level: {len(llt_by_variant)} variants with treatment data")

# 5. Domain-level treatment and KM data (for positions without variant-level data)
treatment_data = {}
for _, tr in treatment_df.iterrows():
    d = tr['domain']
    if pd.notna(tr.get('mean_ldl_untreated')) and pd.notna(tr.get('mean_ldl_treated')):
        treatment_data[d] = {
            'n_untreated': safe(tr['n_untreated']),
            'n_treated': safe(tr['n_treated']),
            'ldl_untreated': safe(tr['mean_ldl_untreated']),
            'ldl_treated': safe(tr['mean_ldl_treated']),
            'ldl_diff': safe(tr['ldl_diff']),
            'p_value': safe(tr['p_value']),
        }

km_domain = {}
for _, kr in km_df[km_df['analysis'] == 'KM_domain'].iterrows():
    km_domain[kr['group']] = {
        'n': safe(kr['n']), 'events': safe(kr['events']),
        'rate': safe(kr['rate']), 'penetrance_60': safe(kr['value']),
    }
km_tertile = {}
for _, kr in km_df[km_df['analysis'] == 'KM_SSS_tertile'].iterrows():
    km_tertile[kr['group']] = {
        'n': safe(kr['n']), 'events': safe(kr['events']),
        'penetrance_60': safe(kr['value']),
    }

# VUS summary stats
vus_summary = {}
for _, vr in vus_df.iterrows():
    a = str(vr['analysis']) if pd.notna(vr.get('analysis')) else ''
    g = str(vr['group']) if pd.notna(vr.get('group')) else ''
    if a == 'clinical_comparison' and g == 'VUS_HIGH_SSS':
        vus_summary['high_sss'] = {
            'n': safe(vr['n']), 'ascvd_rate': safe(vr['ascvd_rate']),
            'ldl_mean': safe(vr['ldl_mean']),
            'simon_broome_definite_rate': safe(vr['simon_broome_definite_rate']),
        }
    elif a == 'clinical_comparison' and g == 'VUS_LOW_SSS':
        vus_summary['low_sss'] = {
            'n': safe(vr['n']), 'ascvd_rate': safe(vr['ascvd_rate']),
            'ldl_mean': safe(vr['ldl_mean']),
            'simon_broome_definite_rate': safe(vr['simon_broome_definite_rate']),
        }
    elif a == 'reclassification_summary':
        vus_summary['reclass'] = {
            'to_pathogenic': safe(vr.get('reclass_to_pathogenic')),
            'to_benign': safe(vr.get('reclass_to_benign')),
        }
    elif a == 'SSS_predicts_reclass':
        vus_summary['prediction'] = {
            'sensitivity': safe(vr.get('sensitivity')),
            'specificity': safe(vr.get('specificity')),
            'ppv': safe(vr.get('ppv')),
            'npv': safe(vr.get('npv')),
        }

# Domain name mappings
DOMAIN_TO_TREAT = {
    'LA1':'Ligand-binding','LA2':'Ligand-binding','LA3':'Ligand-binding R3',
    'LA4':'Ligand-binding','LA5':'Ligand-binding R5','LA6':'Ligand-binding R6',
    'LA7':'Ligand-binding R7','EGF-A':'EGF-like A','EGF-B':'EGF-like B',
    'Beta-propeller':'Beta-propeller','EGF-C':'Linker (EGF-C/O-linked)',
    'O-linked sugars':'O-linked sugar','Transmembrane':'Transmembrane',
    'Cytoplasmic':'Cytoplasmic','Signal peptide':'Unknown',
}
DOMAIN_TO_KM = {
    'LA1':'Ligand-binding R1','LA2':'Ligand-binding R2','LA3':'Ligand-binding R3',
    'LA4':'Ligand-binding R4','LA5':'Ligand-binding R5','LA6':'Ligand-binding R6',
    'LA7':'Ligand-binding R7','EGF-A':'EGF-like A','EGF-B':'EGF-like B',
    'Beta-propeller':'Beta-propeller','EGF-C':'EGF-like C',
    'O-linked sugars':'O-linked sugar','Transmembrane':'Transmembrane',
    'Cytoplasmic':'Cytoplasmic','Signal peptide':'Unknown',
}

def get_traffic_light(row):
    """Evidence-based traffic light using ALL available data sources.

    Sources considered:
    1. Structural: ddG (FoldX), mechanistic class, interface proximity
    2. Functional: Islam et al. cell-based activity, Tabet DMS score
    3. Clinical: ClinVar, ClinGen VCEP, domain penetrance
    4. Computational: AlphaMissense, AlphaGenome, eSSS

    GREEN override: multiple independent sources agree variant is benign
    RED: converging evidence of severe pathogenicity
    AMBER: insufficient or conflicting evidence
    """
    p = int(row.get('position', 0))
    pen = row.get('domain_penetrance_60')
    ddg = row.get('foldx_max_ddG')
    mc = str(row.get('mechanistic_class', '')) if pd.notna(row.get('mechanistic_class')) else ''
    contact = bool(row.get('ldl_contact', False))
    risk = str(row.get('atlas_risk', '')) if pd.notna(row.get('atlas_risk')) else ''
    aa = str(row.get('wildtype_aa', ''))

    # Gather ALL evidence for this position
    islam_act = row.get('islam_mean_activity')
    tabet = tabet_by_pos.get(p, {})
    dms_mean = tabet.get('dms_mean')
    dms_pct_dam = tabet.get('dms_pct_damaging')
    ag = ag_by_pos.get(p, {})
    esss = esss_by_pos.get(p, {})

    # Count benign evidence sources
    benign_votes = 0
    pathogenic_votes = 0

    # Islam: activity >70% = non-defective → benign
    if pd.notna(islam_act) and islam_act > 70:
        benign_votes += 1
    elif pd.notna(islam_act) and islam_act < 40:
        pathogenic_votes += 1

    # Tabet DMS: mean >0.8 = mostly normal → benign; <0.3 = mostly damaging → pathogenic
    if dms_mean is not None and dms_mean > 0.8:
        benign_votes += 1
    elif dms_mean is not None and dms_mean < 0.3:
        pathogenic_votes += 1

    # FoldX ddG: <1 = stable; >4 = severe destabilisation
    if pd.notna(ddg) and ddg < 1:
        benign_votes += 1
    elif pd.notna(ddg) and ddg > 4:
        pathogenic_votes += 1

    # eSSS: n_severe > 5 → concerning; n_severe == 0 and mean < 0.25 → benign
    if esss.get('n_severe', 0) > 5:
        pathogenic_votes += 1
    elif esss.get('esss_mean') is not None and esss['esss_mean'] < 0.25 and esss.get('n_severe', 0) == 0:
        benign_votes += 1

    # GREEN OVERRIDE: ≥3 independent benign sources AND no pathogenic votes
    if benign_votes >= 3 and pathogenic_votes == 0:
        return 'green'

    # RED: severe structural/functional evidence
    if risk == 'Critical':
        return 'red'
    if pd.notna(pen) and pen > 25:
        return 'red'
    if pd.notna(ddg) and ddg > 4 and pathogenic_votes >= 2:
        return 'red'
    if mc == 'recycling_failure':
        return 'red'
    if mc == 'binding_failure' and contact and aa == 'C':
        return 'red'
    # Tabet DMS confirms pathogenicity at position level
    if dms_pct_dam is not None and dms_pct_dam > 50 and pathogenic_votes >= 2:
        return 'red'

    # GREEN: multiple lines of benign evidence
    if benign_votes >= 2 and pathogenic_votes == 0:
        return 'green'
    if (risk in ['Low', ''] and (pd.isna(pen) or pen < 12) and
        (pd.isna(ddg) or ddg < 1) and not contact and
        mc in ['structural_misfolding', ''] and benign_votes >= 1):
        return 'green'

    return 'amber'

# ── Build JSON atlas ──────────────────────────────────────────────────
print("[3/5] Building JSON data...")
atlas = {}
for _, r in pos.iterrows():
    p = int(r['position'])
    wt = r['wildtype_aa']
    p_not = f"p.{AA1TO3.get(wt,'?')}{p}"
    c_pos = (p - 1) * 3 + 1
    c_not = f"c.{c_pos}"
    dom_v3b = str(r.get('domain_v3b', r.get('domain', '')))

    pos_vars = var[(var['position'] == p) & var['variant_id'].notna()]
    variant_list = []
    for _, vr in pos_vars.iterrows():
        vid = str(vr['variant_id'])
        vtype = str(vr.get('variant_type', '')) if pd.notna(vr.get('variant_type')) else None
        sss_val = safe(vr.get('sss_v2'))

        # Get ClinVar classification
        clinvar_class = vus_class_map.get(vid)  # from external_vus_validation
        if clinvar_class:
            # Map VUS-like → VUS, Pathogenic-like → Pathogenic, Benign-like → Benign
            clinvar_label = clinvar_class.replace('-like', '')
        else:
            clinvar_label = None

        # VUS reclassification logic
        vus_reclass = None
        if clinvar_label and clinvar_label == 'VUS':
            if sss_val is not None and sss_val > 0.5:
                vus_reclass = 'likely_pathogenic'
            elif sss_val is not None and sss_val < 0.3:
                vus_reclass = 'likely_benign'
            else:
                vus_reclass = 'uncertain'

        # Get variant-level LLT response
        llt_data = llt_by_variant.get(vid)

        variant_list.append({
            'id': vid,
            'c_notation': str(vr['c_notation']) if pd.notna(vr.get('c_notation')) else None,
            'sss': sss_val,
            'ddG': safe(vr.get('ddG')),
            'type': vtype,
            'clinvar': clinvar_label,
            'vus_reclass': vus_reclass,
            'llt': llt_data,
        })

    # Islam et al data at this position
    islam_data = islam_by_pos.get(p)

    # Treatment/KM lookups
    treat_key = DOMAIN_TO_TREAT.get(dom_v3b, dom_v3b)
    km_key = DOMAIN_TO_KM.get(dom_v3b, dom_v3b)
    treat = treatment_data.get(treat_key)
    km = km_domain.get(km_key)
    tl = get_traffic_light(r)

    atlas[str(p)] = {
        'wt_aa': wt, 'p_notation': p_not, 'c_notation': c_not,
        'amino_acid_pos': p, 'nucleotide_pos': c_pos, 'domain': dom_v3b,
        'mech_class': str(r.get('mechanistic_class', '')) if pd.notna(r.get('mechanistic_class')) else '',
        'binding_site': str(r.get('binding_site', '')) if pd.notna(r.get('binding_site')) else '',
        'ldl_contact': bool(r.get('ldl_contact', False)),
        'plddt': safe(r.get('plddt')),
        'pcsk9_dist': safe(r.get('pcsk9_interface_dist')),
        'at_pcsk9': str(r.get('at_pcsk9_interface', '')) == 'Yes',
        'penetrance_60': safe(r.get('domain_penetrance_60')),
        'ascvd_rate': safe(r.get('domain_ascvd_rate')),
        'max_ddg': safe(r.get('foldx_max_ddG')),
        'n_patients': safe(r.get('n_patients')),
        'risk': str(r.get('atlas_risk', '')) if pd.notna(r.get('atlas_risk')) else '',
        'treatment': str(r.get('atlas_treatment', '')) if pd.notna(r.get('atlas_treatment')) else '',
        'drug_target': str(r.get('atlas_drug_target', '')) if pd.notna(r.get('atlas_drug_target')) else '',
        'rosetta': safe(r.get('rosetta_total')),
        'islam_activity': safe(r.get('islam_mean_activity')),
        'sat_mean_ddG': safe(r.get('sat_mean_ddG')),
        'sat_worst_ddG': safe(r.get('sat_worst_ddG')),
        'sat_worst_mut': str(r.get('sat_worst_mutation', '')) if pd.notna(r.get('sat_worst_mutation')) else None,
        'sat_pct_destab': safe(r.get('sat_pct_destabilising')),
        'jensson_years': safe(r.get('jensson_lifespan_years')),
        'jensson_sig': bool(r.get('jensson_significant')) if pd.notna(r.get('jensson_significant')) else None,
        'reimund_site': str(r.get('reimund_binding_site', '')) if pd.notna(r.get('reimund_binding_site')) else None,
        'reimund_func': str(r.get('reimund_functional_effect', '')) if pd.notna(r.get('reimund_functional_effect')) else None,
        'traffic_light': tl,
        'treat': treat,
        'km': km,
        'islam': islam_data,
        'variants': variant_list,
        # NEW: AlphaGenome regulatory data
        'alphagenome': ag_by_pos.get(p),
        # NEW: Tabet DMS at position level
        'tabet_dms': tabet_by_pos.get(p),
        # NEW: enhanced SSS
        'esss': esss_by_pos.get(p),
        # NEW: Confidence interval on penetrance
        'penetrance_ci': (lambda n, e: {
            'ci_low': round(max(0, e/n - 1.96 * ((e/n * (1-e/n)) / n)**0.5) * 100, 1),
            'ci_high': round(min(1, e/n + 1.96 * ((e/n * (1-e/n)) / n)**0.5) * 100, 1),
        })(km['n'], km['events']) if km and km.get('n') and km.get('events') and km['n'] > 0 else None,
        # NEW: Evidence summary for traffic light
        'evidence_summary': {
            'clinvar': any(v.get('clinvar') for v in variant_list),
            'islam': islam_data is not None,
            'tabet_dms': p in tabet_by_pos,
            'alphagenome': p in ag_by_pos,
            'foldx': pd.notna(r.get('foldx_max_ddG')),
            'esss': p in esss_by_pos,
            'n_sources': sum([
                any(v.get('clinvar') for v in variant_list),
                islam_data is not None,
                p in tabet_by_pos,
                p in ag_by_pos,
                pd.notna(r.get('foldx_max_ddG')),
            ]),
        },
    }

domain_summary = {}
for d in pos['domain_v3b'].unique():
    sub = pos[pos['domain_v3b'] == d]
    treat_key = DOMAIN_TO_TREAT.get(d, d)
    km_key = DOMAIN_TO_KM.get(d, d)
    domain_summary[d] = {
        'start': int(sub['position'].min()), 'end': int(sub['position'].max()),
        'n_residues': len(sub),
        'ldl_contact': bool(sub['ldl_contact'].iloc[0]),
        'mech_class': str(sub['mechanistic_class'].iloc[0]),
        'binding_site': str(sub['binding_site'].iloc[0]),
        'mean_plddt': round(float(sub['plddt'].mean()), 1),
        'penetrance': safe(sub['domain_penetrance_60'].iloc[0]),
        'treat': treatment_data.get(treat_key),
        'km': km_domain.get(km_key),
    }

data = {
    'atlas': atlas, 'domains': domain_summary, 'vus': vus_summary,
    'km_tertile': km_tertile, 'treatment': treatment_data,
}
data_json = json.dumps(data, separators=(',', ':'))
print(f"   JSON size: {len(data_json):,} bytes")
print(f"   Positions: {len(atlas)}")
n_clinvar = sum(1 for p in atlas.values() for v in p['variants'] if v.get('clinvar'))
n_llt = sum(1 for p in atlas.values() for v in p['variants'] if v.get('llt'))
n_islam = sum(1 for p in atlas.values() if p.get('islam'))
print(f"   Variants with ClinVar: {n_clinvar}")
print(f"   Variants with LLT data: {n_llt}")
print(f"   Positions with Islam validation: {n_islam}")

# ── Build HTML ────────────────────────────────────────────────────────
print("[4/5] Building website HTML...")

html_template = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>CALON-FH: Clinical Structural Atlas of LDLR</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/pdbe-molstar@3.10.1/build/pdbe-molstar.css">
<script src="https://cdn.jsdelivr.net/npm/pdbe-molstar@3.10.1/build/pdbe-molstar-plugin.js"></script>
<style>
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:'Segoe UI',Arial,Helvetica,sans-serif;background:#0a0a1a;color:#e0e0e0;line-height:1.5}
.header{background:linear-gradient(135deg,#1a1a3e 0%,#0d2137 100%);padding:30px;text-align:center;border-bottom:2px solid #2196F3}
.header h1{font-size:24px;color:#64b5f6;margin-bottom:6px;letter-spacing:0.5px}
.header p{color:#90a4ae;font-size:13px;max-width:750px;margin:0 auto}
.badge{display:inline-block;background:#1565c0;color:#fff;padding:3px 10px;border-radius:12px;font-size:11px;margin:3px 2px}
.badge.g{background:#2e7d32}.badge.r{background:#c62828}.badge.y{background:#e65100}
.container{max-width:1400px;margin:0 auto;padding:20px}
.search-box{background:#111128;border:1px solid #333;border-radius:8px;padding:20px;margin-bottom:20px}
.search-box input{width:100%;padding:12px 16px;font-size:16px;background:#1a1a2e;border:1px solid #444;border-radius:6px;color:#fff;outline:none}
.search-box input:focus{border-color:#2196F3}
.search-box .hint{color:#666;font-size:12px;margin-top:6px}
.domain-map{background:#111128;border:1px solid #333;border-radius:8px;padding:20px;margin-bottom:20px}
.domain-map h2{color:#64b5f6;font-size:16px;margin-bottom:12px}
.domain-bar{display:flex;height:48px;border-radius:6px;overflow:hidden;margin-bottom:8px}
.domain-bar .segment{position:relative;display:flex;align-items:center;justify-content:center;font-size:9px;font-weight:bold;color:#fff;transition:all 0.2s;border-right:1px solid rgba(0,0,0,0.3);cursor:pointer;touch-action:pan-y}
.domain-bar .segment:hover{filter:brightness(1.3);transform:scaleY(1.1);z-index:1}
.domain-bar .segment .tooltip{display:none;position:absolute;bottom:55px;left:50%;transform:translateX(-50%);background:#1a1a3e;border:1px solid #2196F3;padding:8px 12px;border-radius:6px;font-size:11px;white-space:nowrap;z-index:100;font-weight:normal}
.domain-bar .segment:hover .tooltip{display:block}
.domain-legend{display:flex;flex-wrap:wrap;gap:6px;margin-top:8px}
.domain-legend .item{font-size:11px;display:flex;align-items:center;gap:4px}
.domain-legend .dot{width:10px;height:10px;border-radius:3px}
.mech-bar{display:flex;height:28px;border-radius:4px;overflow:hidden;margin:8px 0}
.mech-bar .seg{display:flex;align-items:center;justify-content:center;font-size:10px;color:#fff}
.stats-row{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));gap:10px;margin-bottom:20px}
.stat-card{background:#111128;border:1px solid #333;border-radius:8px;padding:12px;text-align:center}
.stat-card .value{font-size:24px;font-weight:bold;color:#64b5f6}
.stat-card .label{font-size:10px;color:#888;margin-top:3px}
.result{display:none}.result.active{display:block}
.card{background:#111128;border:1px solid #333;border-radius:8px;padding:20px;margin-bottom:15px}
.card h3{color:#64b5f6;font-size:16px;margin-bottom:10px}
.card-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:12px}
.field{margin-bottom:8px}
.field .label{font-size:10px;color:#888;text-transform:uppercase;letter-spacing:0.5px}
.field .val{font-size:14px;color:#e0e0e0;margin-top:2px}
.field .val.hi{color:#ff9800;font-weight:bold}
.field .val.good{color:#4caf50}.field .val.bad{color:#f44336}.field .val.warn{color:#ff9800}
.tag{display:inline-block;padding:2px 8px;border-radius:4px;font-size:11px;font-weight:bold;margin:2px}
.tag.bs1{background:#c62828;color:#fff}.tag.bs2{background:#e65100;color:#fff}
.tag.contact{background:#2e7d32;color:#fff}.tag.no-contact{background:#455a64;color:#fff}
.tag.binding{background:#b71c1c;color:#fff}.tag.recycling{background:#4a148c;color:#fff}
.tag.misfolding{background:#1565c0;color:#fff}.tag.pivot{background:#6a1b9a;color:#fff}
.tag.interface{background:#e65100;color:#fff}
.tag.path{background:#c62828;color:#fff}.tag.lpath{background:#e65100;color:#fff}
.tag.vus{background:#7c4dff;color:#fff}.tag.benign{background:#2e7d32;color:#fff}
.tag.lbenign{background:#00897b;color:#fff}
.var-table{width:100%;border-collapse:collapse;margin-top:8px}
.var-table th{background:#1a1a3e;padding:6px 8px;text-align:left;font-size:11px;color:#90a4ae;border-bottom:1px solid #333}
.var-table td{padding:5px 8px;font-size:12px;border-bottom:1px solid #222}
.var-table tr:hover{background:rgba(33,150,243,0.1)}
.notation{font-family:'Courier New',monospace;background:#1a1a3e;padding:2px 6px;border-radius:3px;color:#64b5f6;font-size:13px}
.footer{text-align:center;padding:30px;color:#555;font-size:11px;border-top:1px solid #222;margin-top:40px}
.footer a{color:#2196F3;text-decoration:none}
.pen-bar{height:18px;border-radius:4px;background:#222;position:relative;margin-top:4px}
.pen-bar .fill{height:100%;border-radius:4px}
.pen-bar .text{position:absolute;right:6px;top:1px;font-size:11px;color:#fff}
/* Traffic light */
.tl-box{display:flex;align-items:center;gap:18px;padding:18px 22px;border-radius:10px;margin-bottom:15px}
.tl-box.red{background:linear-gradient(135deg,#2d0a0a 0%,#1a0505 100%);border:2px solid #c62828}
.tl-box.amber{background:linear-gradient(135deg,#2d1f0a 0%,#1a1205 100%);border:2px solid #f9a825}
.tl-box.green{background:linear-gradient(135deg,#0a2d0f 0%,#051a08 100%);border:2px solid #2e7d32}
.tl-light{width:50px;height:50px;border-radius:50%;flex-shrink:0;display:flex;align-items:center;justify-content:center;font-size:22px}
.tl-light.red{background:radial-gradient(circle,#ff1744 30%,#c62828 100%);box-shadow:0 0 18px rgba(255,23,68,0.5)}
.tl-light.amber{background:radial-gradient(circle,#ffc400 30%,#f9a825 100%);box-shadow:0 0 18px rgba(255,196,0,0.5)}
.tl-light.green{background:radial-gradient(circle,#69f0ae 30%,#2e7d32 100%);box-shadow:0 0 18px rgba(105,240,174,0.4)}
.tl-text{flex:1}
.tl-text h2{font-size:18px;margin-bottom:4px}
.tl-text.red h2{color:#ff5252}.tl-text.amber h2{color:#ffc400}.tl-text.green h2{color:#69f0ae}
.tl-text p{font-size:13px;color:#b0bec5;line-height:1.5}
/* Panels */
.panel{background:#0d1b2a;border:1px solid #1e3a5f;border-radius:10px;padding:18px;margin-bottom:14px}
.panel h3{color:#64b5f6;font-size:15px;margin-bottom:10px;display:flex;align-items:center;gap:6px}
.panel-grid{display:grid;grid-template-columns:1fr 1fr;gap:12px}
@media(max-width:700px){.panel-grid{grid-template-columns:1fr}}
.metric{text-align:center;padding:10px;background:#111128;border-radius:8px;border:1px solid #222}
.metric .big{font-size:26px;font-weight:bold}
.metric .sub{font-size:10px;color:#888;margin-top:3px}
.metric.red .big{color:#ff5252}.metric.amber .big{color:#ffc400}.metric.green .big{color:#69f0ae}
/* LLT bars */
.llt-bar-wrap{display:flex;align-items:center;gap:8px;margin:5px 0}
.llt-bar{flex:1;height:22px;background:#1a1a2e;border-radius:4px;position:relative;overflow:hidden}
.llt-bar .fill{height:100%;border-radius:4px;transition:width 0.5s}
.llt-bar .txt{position:absolute;right:8px;top:2px;font-size:11px;color:#fff;font-weight:bold}
.llt-label{width:70px;font-size:11px;color:#90a4ae;text-align:right}
/* VUS banner */
.vus-banner{background:linear-gradient(135deg,#1a0d2e 0%,#0d1a2e 100%);border:2px solid #7c4dff;border-radius:10px;padding:16px;margin-bottom:14px}
.vus-banner h3{color:#b388ff;margin-bottom:8px;font-size:15px}
.vus-banner .reclass{font-size:18px;font-weight:bold;margin:4px 0}
.vus-banner .reclass.pathogenic{color:#ff5252}
.vus-banner .reclass.benign{color:#69f0ae}
.vus-banner .reclass.uncertain{color:#ffc400}
.vus-stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(90px,1fr));gap:8px;margin-top:10px}
.vus-stat{text-align:center;padding:6px;background:rgba(0,0,0,0.3);border-radius:6px}
.vus-stat .num{font-size:16px;font-weight:bold;color:#b388ff}
.vus-stat .lab{font-size:9px;color:#888}
/* Islam validation */
.islam-box{background:linear-gradient(135deg,#0d2a1b 0%,#0d1a2e 100%);border:1px solid #2e7d32;border-radius:10px;padding:16px;margin-bottom:14px}
.islam-box h3{color:#81c784;font-size:15px;margin-bottom:10px}
.islam-row{display:flex;align-items:center;gap:12px;padding:6px 0;border-bottom:1px solid rgba(255,255,255,0.05)}
.islam-row:last-child{border-bottom:none}
.func-badge{padding:3px 10px;border-radius:12px;font-size:11px;font-weight:bold;display:inline-block}
.func-badge.lof{background:#c62828;color:#fff}
.func-badge.def{background:#e65100;color:#fff}
.func-badge.mild{background:#f9a825;color:#000}
.func-badge.nondef{background:#2e7d32;color:#fff}
.func-badge.g7090{background:#ff9800;color:#000}
/* 3D viewer */
.viewer-container{width:100%;height:380px;border-radius:8px;overflow:hidden;background:#000;position:relative;margin-bottom:12px}
.viewer-overlay{position:absolute;top:0;left:0;right:0;bottom:0;display:flex;align-items:center;justify-content:center;background:rgba(0,0,0,0.7);z-index:10;cursor:pointer;transition:opacity 0.3s}
.viewer-overlay span{color:#64b5f6;font-size:15px}
.viewer-overlay:hover span{text-decoration:underline}
/* Explanation */
.expl-box{background:linear-gradient(135deg,#0d1b2a 0%,#1b2838 100%);border:1px solid #1e88e5;border-radius:10px;padding:18px;margin-top:14px}
.expl-pathway{display:flex;flex-direction:column;gap:0}
.expl-step{display:flex;gap:12px;padding:12px 0;border-left:3px solid #1e88e5;margin-left:12px;padding-left:18px;position:relative}
.expl-step:last-child{border-left-color:transparent}
.expl-step::before{content:'';position:absolute;left:-7px;top:12px;width:11px;height:11px;border-radius:50%;background:#1e88e5;border:2px solid #0d1b2a}
.expl-num{font-size:11px;font-weight:bold;color:#1e88e5;flex-shrink:0;width:20px;text-align:right}
.expl-text{font-size:13px;line-height:1.6;color:#cfd8dc}
.expl-text b{color:#e0e0e0}
.src{font-size:11px;color:#78909c;font-style:italic}
</style>
</head>
<body>
<div class="header">
<h1>CALON-FH Clinical Atlas &mdash; LDL Receptor Structural-Clinical Integration</h1>
<p>Precision tool for cardiologists: mutation location &rarr; structural mechanism &rarr; ASCVD penetrance &rarr; treatment prediction. Built from AlphaFold3, FoldX, UK Biobank WES, Wales FH Registry, and DRAGON3 cohort.</p>
<div style="margin-top:10px">
<span class="badge">860 LDLR Positions</span>
<span class="badge">331 Characterised Variants</span>
<span class="badge">~4,000 FH Patients</span>
<span class="badge r">ClinVar Integrated</span>
<span class="badge g">LLT Response Validated</span>
<span class="badge y">78 VUS Reclassified</span>
<span class="badge" style="background:#7b1fa2">Tabet DMS (Science 2025)</span>
<span class="badge" style="background:#ff6f00">AlphaGenome Regulatory</span>
</div>
</div>
<div class="container">
<div class="stats-row">
<div class="stat-card"><div class="value">860</div><div class="label">LDLR Residues</div></div>
<div class="stat-card"><div class="value">331</div><div class="label">Variants Mapped</div></div>
<div class="stat-card"><div class="value">113</div><div class="label">ClinVar Classified</div></div>
<div class="stat-card"><div class="value">56</div><div class="label">Variants with LLT Data</div></div>
<div class="stat-card"><div class="value">234</div><div class="label">Islam et al. Validated</div></div>
<div class="stat-card"><div class="value">78</div><div class="label">VUS Reclassified</div></div>
<div class="stat-card"><div class="value">38.1%</div><div class="label">Peak Penetrance (EGF-A)</div></div>
</div>
<div class="domain-map">
<h2>LDLR Domain Architecture (Reimund et al., Nature 2025)</h2>
<div class="domain-bar" id="domainBar"></div>
<div style="display:flex;justify-content:space-between;font-size:10px;color:#666;margin-top:3px">
<span>N-terminus (1)</span><span>C-terminus (860)</span>
</div>
<div style="margin-top:12px">
<div style="font-size:12px;color:#90a4ae;margin-bottom:5px">Mechanistic Classification</div>
<div class="mech-bar" id="mechBar"></div>
</div>
<div class="domain-legend" id="domainLegend"></div>
</div>
<div class="search-box">
<input type="text" id="searchInput" placeholder="Search: 352, p.Cys352, c.1055, EGF-A, binding_failure, VUS, Pathogenic, red...">
<div class="hint">Position &bull; p./c. notation &bull; Domain &bull; Mechanism &bull; ClinVar class &bull; Traffic light colour</div>
</div>
<div class="result" id="resultCard"></div>
<div id="searchResults"></div>
</div>
<div class="footer">
<p><b>CALON-FH Clinical Atlas v3b</b> &mdash; Dr Nader Genedy</p>
<p>Domain boundaries: Reimund et al., <i>Nature</i> 638:829 (2025) &bull; Lifespan: Jensson et al., <i>NEJM</i> 389:1741 (2023)</p>
<p>Functional validation: Islam et al. (2024) &bull; Tabet et al., <i>Science</i> 2025;391:eady7186 &bull; Structure: AlphaFold3 + FoldX 5.0 + Rosetta</p>
<p>Regulatory: AlphaGenome v0.6.1 (Google DeepMind) &bull; Clinical data: UK Biobank WES + Wales FH Registry + DRAGON3 Cohort</p>
<p><a href="https://github.com/NaderGenedy/calon-fh-atlas">GitHub</a></p>
</div>
<script>
const DATA = __DATA_PLACEHOLDER__;

const DOMAIN_COLORS = {
'Signal peptide':'#546e7a','LA1':'#78909c','LA2':'#90a4ae',
'LA3':'#c62828','LA4':'#d32f2f','LA5':'#e53935','LA6':'#ef5350','LA7':'#b71c1c',
'EGF-A':'#ff6f00','EGF-B':'#7b1fa2','Beta-propeller':'#e65100',
'EGF-C':'#827717','O-linked sugars':'#00695c','Transmembrane':'#1565c0','Cytoplasmic':'#283593'
};
const MECH_COLORS = {
'binding_failure':'#c62828','recycling_failure':'#4a148c',
'interface_BS2':'#e65100','recycling_pivot':'#6a1b9a',
'structural_misfolding':'#1565c0'
};
const MECH_LABELS = {
'binding_failure':'Binding Failure (BS1)',
'recycling_failure':'Recycling Failure (pH switch)',
'interface_BS2':'Interface BS2',
'recycling_pivot':'Recycling Pivot (EGF-B)',
'structural_misfolding':'Structural Misfolding'
};

var _touchMoved=false,_molViewer=null,_molReady=false,_molPos=null;
document.addEventListener('touchmove',function(){_touchMoved=true},{passive:true});
document.addEventListener('touchstart',function(){_touchMoved=false},{passive:true});

function initMolstar(pos){
  var c=document.getElementById('mol-viewer');
  if(!c)return;
  var o=document.getElementById('mol-overlay');
  if(o)o.style.display='none';
  _molPos=pos;
  if(_molViewer&&_molReady){highlightResidue(pos);return;}
  if(_molViewer)return;
  _molViewer=new PDBeMolstarPlugin();
  _molViewer.render(c,{
    customData:{url:'https://alphafold.ebi.ac.uk/files/AF-P01130-F1-model_v6.cif',format:'cif'},
    bgColor:{r:10,g:10,b:26},hideControls:true,
    hideCanvasControls:['expand','selection','animation'],
    landscape:true,alphafoldView:true,subscribeEvents:true
  });
  _molViewer.events.loadComplete.subscribe(function(){_molReady=true;highlightResidue(_molPos);});
}
function highlightResidue(pos){
  if(!_molViewer||!_molReady)return;
  var p=parseInt(pos);
  _molViewer.visual.clearSelection();
  _molViewer.visual.select({data:[{struct_asym_id:'A',start_residue_number:p,end_residue_number:p,color:{r:255,g:23,b:68},focus:true}]});
}

function buildDomainBar(){
  var bar=document.getElementById('domainBar'),legend=document.getElementById('domainLegend'),domains=DATA.domains;
  for(var name in domains){
    var d=domains[name],pct=d.n_residues/860*100;
    var seg=document.createElement('div');
    seg.className='segment';seg.style.width=pct+'%';seg.style.background=DOMAIN_COLORS[name]||'#555';
    if(pct>3)seg.textContent=name;
    var tt=document.createElement('div');tt.className='tooltip';
    var penText=d.km?d.km.penetrance_60.toFixed(1)+'%':'N/A';
    var treatText=d.treat?'\u2212'+((d.treat.ldl_untreated-d.treat.ldl_treated)/d.treat.ldl_untreated*100).toFixed(0)+'% expected reduction':'No data';
    tt.innerHTML='<b>'+name+'</b><br>Residues '+d.start+'\u2013'+d.end+' ('+d.n_residues+')<br>Mechanism: '+d.mech_class+'<br>ASCVD penetrance (60y): '+penText+'<br>Expected LLT response: '+treatText;
    seg.appendChild(tt);seg.setAttribute('data-domain',name);
    seg.onclick=function(e){if(_touchMoved){e.preventDefault();return;}searchDomain(this.getAttribute('data-domain'))};
    bar.appendChild(seg);
    var item=document.createElement('div');item.className='item';
    item.innerHTML='<div class="dot" style="background:'+(DOMAIN_COLORS[name]||'#555')+'"></div>'+name;
    legend.appendChild(item);
  }
  var mechBar=document.getElementById('mechBar'),mechCounts={};
  for(var pos in DATA.atlas){var mc=DATA.atlas[pos].mech_class||'unknown';mechCounts[mc]=(mechCounts[mc]||0)+1;}
  for(var mc in mechCounts){
    if(!mc)continue;var pct2=mechCounts[mc]/860*100;
    var seg2=document.createElement('div');seg2.className='seg';seg2.style.width=pct2+'%';seg2.style.background=MECH_COLORS[mc]||'#555';
    if(pct2>5)seg2.textContent=(MECH_LABELS[mc]||mc)+' ('+mechCounts[mc]+')';
    mechBar.appendChild(seg2);
  }
}
function searchDomain(name){document.getElementById('searchInput').value=name;doSearch(name);}
function getMechTag(mc){
  var cls=mc==='binding_failure'?'binding':mc==='recycling_failure'?'recycling':mc==='recycling_pivot'?'pivot':mc==='interface_BS2'?'interface':'misfolding';
  return '<span class="tag '+cls+'">'+(MECH_LABELS[mc]||mc)+'</span>';
}
function getClinvarTag(cv){
  if(!cv)return '<span class="tag no-contact">No ClinVar</span>';
  if(cv==='Pathogenic')return '<span class="tag path">Pathogenic</span>';
  if(cv==='Likely pathogenic')return '<span class="tag lpath">Likely Pathogenic</span>';
  if(cv==='VUS')return '<span class="tag vus">VUS</span>';
  if(cv==='Benign')return '<span class="tag benign">Benign</span>';
  if(cv==='Likely benign')return '<span class="tag lbenign">Likely Benign</span>';
  return '<span class="tag no-contact">'+cv+'</span>';
}
function getPenBar(pct,label){
  if(!pct)return '';
  var color=pct>30?'#c62828':pct>20?'#e65100':pct>10?'#f9a825':'#2e7d32';
  return '<div class="pen-bar"><div class="fill" style="width:'+Math.min(pct,100)+'%;background:'+color+'"></div><div class="text">'+(label||'')+pct.toFixed(1)+'%</div></div>';
}
function getFuncBadge(fg){
  if(!fg)return '';
  var cls=fg==='Loss-of-function'?'lof':fg==='Defective'?'def':fg==='Mildly-defective'?'mild':fg==='Non-defective'?'nondef':'g7090';
  return '<span class="func-badge '+cls+'">'+fg+'</span>';
}

function getTrafficLightHTML(d){
  var tl=d.traffic_light;
  var title=tl==='red'?'HIGH RISK':tl==='amber'?'MODERATE RISK':'BENIGN / LOWER RISK';
  var desc='',evid='';

  // Build evidence badge bar
  var es=d.evidence_summary||{};
  evid='<div style="display:flex;flex-wrap:wrap;gap:4px;margin-top:8px">';
  if(es.clinvar)evid+='<span class="badge">ClinVar</span>';
  if(es.islam)evid+='<span class="badge g">Islam et al.</span>';
  if(es.tabet_dms)evid+='<span class="badge" style="background:#7b1fa2">Tabet DMS</span>';
  if(es.alphagenome)evid+='<span class="badge" style="background:#ff6f00">AlphaGenome</span>';
  if(es.foldx)evid+='<span class="badge" style="background:#1565c0">FoldX</span>';
  if(es.esss)evid+='<span class="badge" style="background:#00695c">eSSS</span>';
  evid+='<span style="font-size:10px;color:#78909c;margin-left:6px">('+es.n_sources+' independent source'+(es.n_sources!==1?'s':'')+')</span>';
  evid+='</div>';

  if(tl==='red'){
    desc='<b>Converging evidence of clinical severity</b> from multiple independent sources. ';
    if(d.mech_class==='recycling_failure')desc+='Critical pH-switch residue required for receptor recycling. ';
    else if(d.mech_class==='binding_failure'&&d.ldl_contact)desc+='Direct LDL-contact residue at binding interface. ';
    if(d.tabet_dms&&d.tabet_dms.dms_pct_damaging>30)desc+='Tabet DMS: '+d.tabet_dms.dms_pct_damaging.toFixed(0)+'% of substitutions damaging. ';
    if(d.max_ddg&&d.max_ddg>4)desc+='\u0394\u0394G = '+d.max_ddg.toFixed(1)+' kcal/mol (severe destabilisation). ';
    desc+='Consider early aggressive lipid-lowering therapy and cascade screening.';
  }else if(tl==='green'){
    desc='<b>Multiple independent sources support benign/low-risk classification.</b> ';
    if(d.tabet_dms&&d.tabet_dms.dms_mean>0.8)desc+='Tabet DMS: mean score '+d.tabet_dms.dms_mean.toFixed(2)+' (normal function). ';
    if(d.islam){var bestIslam=d.islam.reduce(function(a,b){return(a.activity_pct||0)>(b.activity_pct||0)?a:b},{});if(bestIslam.activity_pct&&bestIslam.activity_pct>70)desc+='Islam et al.: '+bestIslam.activity_pct.toFixed(0)+'% activity (non-defective). ';}
    if(d.max_ddg!=null&&d.max_ddg<1)desc+='\u0394\u0394G < 1 kcal/mol (structurally stable). ';
    desc+='Standard risk-based management. This variant is unlikely to cause FH.';
  }else{
    desc='Moderate or uncertain structural significance. ';
    if(d.penetrance_60){desc+='Domain penetrance: '+d.penetrance_60.toFixed(1)+'%';if(d.penetrance_ci)desc+=' (95% CI: '+d.penetrance_ci.ci_low.toFixed(1)+'\u2013'+d.penetrance_ci.ci_high.toFixed(1)+'%)';desc+='. ';}
    if(d.max_ddg&&d.max_ddg>2)desc+='\u0394\u0394G = '+d.max_ddg.toFixed(1)+' kcal/mol. ';
    desc+='Standard high-intensity statin \u00B1 ezetimibe; escalate based on LDL-C response.';
  }
  return '<div class="tl-box '+tl+'"><div class="tl-light '+tl+'">'+(tl==='red'?'\u26A0':tl==='amber'?'\u26A0':'\u2714')+'</div><div class="tl-text '+tl+'"><h2>'+title+'</h2><p>'+desc+'</p>'+evid+'</div></div>';
}

function getPenetrancePanel(d){
  var km=d.km;
  if(!km)return '';
  var h='<div class="panel"><h3>ASCVD Penetrance &mdash; Wales FH Registry Kaplan-Meier Analysis</h3>';
  h+='<div class="panel-grid"><div>';
  h+='<div class="metric '+(km.penetrance_60>25?'red':km.penetrance_60>15?'amber':'green')+'"><div class="big">'+km.penetrance_60.toFixed(1)+'%</div><div class="sub">ASCVD by age 60</div></div>';
  h+='</div><div>';
  h+='<div class="metric"><div class="big">'+km.n+'</div><div class="sub">carriers ('+km.events+' ASCVD events)</div></div>';
  h+='</div></div>';
  // Show confidence interval
  if(d.penetrance_ci){h+='<div style="margin-top:6px;font-size:12px;color:#b0bec5">95% CI: '+d.penetrance_ci.ci_low.toFixed(1)+'% \u2013 '+d.penetrance_ci.ci_high.toFixed(1)+'%</div>';}
  h+=getPenBar(km.penetrance_60,'');
  var t=DATA.km_tertile;
  if(t){h+='<div style="margin-top:8px;font-size:11px;color:#78909c">SSS tertile reference: Low SSS '+t.T1_Low.penetrance_60.toFixed(1)+'% | Mid '+t.T2_Mid.penetrance_60.toFixed(1)+'% | High '+t.T3_High.penetrance_60.toFixed(1)+'%</div>';}
  h+='<div style="margin-top:6px;font-size:10px;color:#546e7a">\u26A0 Domain-level penetrance from treated Northern European cohort (Wales FH Registry, n=3,562). True untreated risk may be higher. Not validated in non-European populations.</div>';
  h+='</div>';
  return h;
}
function getAlphaGenomePanel(d){
  var ag=d.alphagenome;
  if(!ag)return '';
  var h='<div class="panel" style="border-color:#ff6f00"><h3 style="color:#ffb74d">\uD83E\uDDEC AlphaGenome &mdash; mRNA-Level Regulatory Analysis</h3>';
  h+='<div style="font-size:11px;color:#78909c;margin-bottom:10px">Google DeepMind AlphaGenome v0.6.1 \u2014 predicts variant effects on gene expression, splicing, and chromatin accessibility</div>';
  h+='<div class="panel-grid">';
  if(ag.cage_diff!=null){var ec=Math.abs(ag.cage_diff)>0.1?'warn':'good';h+='<div class="metric"><div class="big '+ec+'">'+(ag.cage_diff>=0?'+':'')+ag.cage_diff.toFixed(3)+'</div><div class="sub">Expression (CAGE)</div></div>';}
  if(ag.splice_max!=null){var sc=ag.splice_max>0.5?'bad':ag.splice_max>0.1?'warn':'good';h+='<div class="metric"><div class="big '+sc+'">'+ag.splice_max.toFixed(3)+'</div><div class="sub">Splice disruption (max)</div></div>';}
  if(ag.atac_diff!=null){var ac=Math.abs(ag.atac_diff)>0.1?'warn':'good';h+='<div class="metric"><div class="big '+ac+'">'+(ag.atac_diff>=0?'+':'')+ag.atac_diff.toFixed(3)+'</div><div class="sub">Accessibility (ATAC)</div></div>';}
  if(ag.l2_impact!=null){var lc=ag.l2_impact>0.3?'bad':ag.l2_impact>0.1?'warn':'good';h+='<div class="metric"><div class="big '+lc+'">'+ag.l2_impact.toFixed(3)+'</div><div class="sub">Overall impact (L2)</div></div>';}
  h+='</div>';
  if(ag.splice_max>0.5)h+='<div style="margin-top:8px;padding:8px;background:rgba(244,67,54,0.1);border-radius:6px;font-size:12px;color:#ef9a9a">\u26A0 <b>Strong splice disruption detected.</b> This variant may trigger nonsense-mediated mRNA decay (NMD) or create aberrant splice products, reducing functional LDLR mRNA in a tissue-specific manner.</div>';
  h+='<div style="margin-top:6px;font-size:10px;color:#546e7a">AlphaGenome predictions are computational (not experimentally validated for LDLR). Treat as hypothesis-generating.</div>';
  h+='</div>';
  return h;
}
function getTabetPanel(d){
  var tb=d.tabet_dms;
  if(!tb)return '';
  var h='<div class="panel" style="border-color:#7b1fa2"><h3 style="color:#ce93d8">\uD83E\uDDEA Tabet et al. 2025 &mdash; Deep Mutational Scanning (Science)</h3>';
  h+='<div style="font-size:11px;color:#78909c;margin-bottom:10px">Cell-based LDL uptake functional scores for all possible substitutions at this position (~17,000 LDLR variants)</div>';
  h+='<div class="panel-grid">';
  if(tb.dms_mean!=null){var mc=tb.dms_mean<0.3?'bad':tb.dms_mean<0.5?'warn':'good';h+='<div class="metric"><div class="big '+mc+'">'+tb.dms_mean.toFixed(2)+'</div><div class="sub">Mean DMS score</div></div>';}
  if(tb.dms_pct_damaging!=null){var dc=tb.dms_pct_damaging>30?'bad':tb.dms_pct_damaging>10?'warn':'good';h+='<div class="metric"><div class="big '+dc+'">'+tb.dms_pct_damaging.toFixed(0)+'%</div><div class="sub">Substitutions damaging</div></div>';}
  if(tb.dms_n!=null){h+='<div class="metric"><div class="big">'+tb.dms_n+'</div><div class="sub">Variants tested</div></div>';}
  if(tb.dms_min!=null){h+='<div class="metric"><div class="big">'+tb.dms_min.toFixed(2)+' \u2013 '+tb.dms_max.toFixed(2)+'</div><div class="sub">Score range</div></div>';}
  h+='</div>';
  if(tb.dms_mean>0.8)h+='<div style="margin-top:8px;padding:8px;background:rgba(76,175,80,0.1);border-radius:6px;font-size:12px;color:#a5d6a7">\u2714 <b>Position is tolerant.</b> Most substitutions retain normal LDL uptake function in cell-based assay.</div>';
  else if(tb.dms_pct_damaging>30)h+='<div style="margin-top:8px;padding:8px;background:rgba(244,67,54,0.1);border-radius:6px;font-size:12px;color:#ef9a9a">\u26A0 <b>Functionally sensitive position.</b> >30% of substitutions impair LDL uptake in cell-based assay.</div>';
  h+='<div style="margin-top:6px;font-size:10px;color:#546e7a">Tabet et al., Science 2025;391:eady7186. DMS threshold: <0.118 = damaging (15th percentile). Gold-standard experimental validation.</div>';
  h+='</div>';
  return h;
}

function getRxRec(d){
  // Mechanism-based medication recommendation from our cohort data
  var mc=d.mech_class||'';var dom=d.domain||'';
  var first='',escalate='',rationale='',caution='';
  if(mc==='binding_failure'&&d.ldl_contact){
    first='Atorvastatin 40\u201380 mg + Ezetimibe 10 mg';
    escalate='Add PCSK9i (evolocumab 140 mg Q2W or alirocumab 150 mg Q2W). If LDL-C remains >2.6 mmol/L: consider evinacumab (LDLR-independent, ANGPTL3 pathway).';
    rationale='Direct LDL-contact residue at BS1. Mutant allele cannot bind LDL regardless of LDLR expression level. Statin benefit comes primarily from upregulating the <b>wild-type allele</b>. Ezetimibe provides additive 10.8% reduction (our cohort: statin alone 3.7% vs statin+ezetimibe 10.8%).';
    caution='Statin monotherapy has limited effect via the mutant allele. Early combination therapy recommended.';
  }else if(mc==='recycling_failure'){
    first='Atorvastatin 40\u201380 mg + Ezetimibe 10 mg';
    escalate='PCSK9i may have <b>reduced incremental benefit</b> (receptor recycling already impaired). If inadequate response: evinacumab or lipoprotein apheresis.';
    rationale='pH-switch residue mutation. Receptor binds and internalises LDL normally but is degraded after a single cycle (vs \u223C100 cycles normally). Statins partially compensate by driving new receptor synthesis, but each receptor molecule is only used once.';
    caution='PCSK9 inhibitors prevent PCSK9-mediated degradation, but these receptors are already lost via the pH-switch recycling failure pathway \u2014 different mechanism.';
  }else if(mc==='recycling_pivot'){
    first='Atorvastatin 40\u201380 mg + Ezetimibe 10 mg';
    escalate='Add PCSK9i. If LDL-C >3.5 mmol/L on triple therapy: lipoprotein apheresis referral.';
    rationale='EGF-B hinge domain mutation impairs the conformational change required for LDL release at endosomal pH. Similar to recycling failure but may retain partial function.';
    caution='Monitor LDL-C response at 6\u20138 weeks. Partial recycling may preserve some statin responsiveness.';
  }else if(mc==='interface_BS2'||dom.indexOf('Beta-propeller')>=0||dom.indexOf('beta_propeller')>=0){
    first='Atorvastatin 40\u201380 mg (or Rosuvastatin 20\u201340 mg)';
    escalate='Add Ezetimibe 10 mg if LDL-C >2.6 mmol/L. Then PCSK9i if LDL-C remains >1.8 mmol/L (very high risk) or >2.6 mmol/L (high risk).';
    rationale='<b>Best documented treatment response in our cohort.</b> Beta-propeller domain: expected \u221233% LDL-C reduction (5.06\u21923.41 mmol/L, P=0.001, n=63). BS1 (primary LDL binding) is intact. PCSK9 binding site (EGF-A) is structurally intact \u2192 PCSK9i expected to be effective.';
    caution='';
  }else if(dom.indexOf('Cytoplasmic')>=0||dom.indexOf('NPXY')>=0){
    first='Evinacumab 15 mg/kg IV Q4W (LDLR-independent)';
    escalate='Lipoprotein apheresis if evinacumab unavailable. Statins + ezetimibe as adjuncts for wild-type allele benefit.';
    rationale='<b>Class 5 internalisation defect.</b> Receptor reaches the cell surface and binds LDL but cannot be recruited to clathrin-coated pits (NPXY motif disrupted). Statins upregulate a non-internalising receptor. PCSK9i prevents degradation of a receptor that is already surface-stable but non-functional for endocytosis.';
    caution='Neither statins nor PCSK9i directly address the endocytic defect. LDLR-independent pathways (ANGPTL3 inhibition) bypass the receptor entirely.';
  }else if(dom.indexOf('EGF-A')>=0||dom.indexOf('EGF-like A')>=0){
    first='Atorvastatin 40\u201380 mg + Ezetimibe 10 mg';
    escalate='PCSK9i: <b>check specific residue</b> \u2014 EGF-A contains the PCSK9 binding epitope. If mutation disrupts PCSK9 binding, PCSK9i may be ineffective (the drug cannot bind its target). Alternative: inclisiran (siRNA, different mechanism) or evinacumab.';
    rationale='EGF-A has a dual role: LDL binding (BS1 component) AND PCSK9 binding. Mutations here may simultaneously impair LDL clearance and alter the efficacy of PCSK9-targeted therapies.';
    caution='AlphaFold3 PCSK9 complex modelling recommended to assess PCSK9i binding impact at this specific position.';
  }else if(d.sss_mean>=0.7){
    first='Atorvastatin 40\u201380 mg + Ezetimibe 10 mg';
    escalate='Add PCSK9i early. If SSS \u22650.9 (likely null): refer for evinacumab or apheresis evaluation.';
    rationale='High structural severity (SSS \u22650.7) suggests significant protein destabilisation or loss-of-function. In our cohort, high-SSS variants in the ligand-binding domain showed only \u22121.9% mean LDL reduction on statin monotherapy.';
    caution='Null/severe mutations: mutant allele contributes minimally to LDL clearance regardless of treatment.';
  }else{
    first='Atorvastatin 20\u201340 mg';
    escalate='Uptitrate to high-intensity (Atorvastatin 80 mg or Rosuvastatin 40 mg). Add Ezetimibe 10 mg if LDL-C >2.6 mmol/L. Add PCSK9i if LDL-C remains above target after 6\u20138 weeks.';
    rationale='Moderate structural risk. Standard FH escalation pathway. Monitor LDL-C at 6\u20138 weeks after each change. Target: <1.8 mmol/L (very high risk) or <2.6 mmol/L (high risk) per ESC/EAS 2019.';
    caution='';
  }
  return {first:first,escalate:escalate,rationale:rationale,caution:caution};
}

function getLLTPanel(d){
  // Check for variant-level LLT first (our actual patient data)
  var hasVariantLLT=false;
  var variantLLT=[];
  if(d.variants){
    for(var i=0;i<d.variants.length;i++){
      if(d.variants[i].llt){hasVariantLLT=true;variantLLT.push(d.variants[i]);}
    }
  }

  var h='<div class="panel"><h3>\uD83D\uDC8A Treatment Strategy &mdash; Mechanism-Based</h3>';

  // Expected response section
  if(hasVariantLLT){
    h+='<div style="margin-bottom:14px;padding:10px;background:rgba(76,175,80,0.06);border-radius:8px;border-left:3px solid #4caf50">';
    h+='<div style="font-size:12px;color:#4caf50;font-weight:600;margin-bottom:8px">OBSERVED RESPONSE (Our Cohort)</div>';
    for(var i=0;i<variantLLT.length;i++){
      var v=variantLLT[i],llt=v.llt;
      var pct=llt.ldl_pct_reduction;
      var resp=pct>35?'Good':'Partial';
      var respCol=pct>35?'#4caf50':'#ff9800';
      h+='<div style="display:flex;align-items:center;gap:12px;flex-wrap:wrap;margin-bottom:6px">';
      h+='<span class="notation">'+v.id+'</span>';
      h+='<span style="font-size:18px;font-weight:700;color:'+respCol+'">\u2212'+pct.toFixed(0)+'%</span>';
      h+='<span style="font-size:13px"><span style="display:inline-block;width:8px;height:8px;border-radius:50%;background:'+respCol+';margin-right:4px"></span>'+resp+' responder (n='+llt.n_treated+', '+llt.intensity+' intensity)</span>';
      h+='</div>';
    }
    h+='<div class="src" style="margin-top:6px">Source: DRAGON3/Wales FH Registry. Validated against Islam et al.</div></div>';
  }else if(d.treat){
    var t=d.treat;
    var pct=((t.ldl_untreated-t.ldl_treated)/t.ldl_untreated*100);
    var resp=pct>30?'Good':'Partial';
    var respCol=pct>30?'#4caf50':'#ff9800';
    h+='<div style="margin-bottom:14px;padding:10px;background:rgba(76,175,80,0.06);border-radius:8px;border-left:3px solid #66bb6a">';
    h+='<div style="font-size:12px;color:#66bb6a;font-weight:600;margin-bottom:8px">EXPECTED RESPONSE (Domain-Level)</div>';
    h+='<div style="display:flex;align-items:center;gap:12px;flex-wrap:wrap">';
    h+='<span style="font-size:18px;font-weight:700;color:'+respCol+'">\u2212'+pct.toFixed(0)+'%</span>';
    h+='<span style="font-size:13px"><span style="display:inline-block;width:8px;height:8px;border-radius:50%;background:'+respCol+';margin-right:4px"></span>'+resp+' responder ('+d.domain+', n='+t.n_treated+', P='+(t.p_value<0.001?'<0.001':t.p_value)+')</span>';
    h+='</div>';
    h+='<div class="src" style="margin-top:6px">Source: UK Biobank + Wales FH. Validated against Islam et al.</div></div>';
  }

  // Medication recommendation
  var rx=getRxRec(d);
  h+='<div style="margin-bottom:10px">';
  h+='<div style="font-size:12px;color:#64b5f6;font-weight:600;margin-bottom:8px">RECOMMENDED FIRST-LINE</div>';
  h+='<div style="padding:10px;background:rgba(100,181,246,0.08);border-radius:8px;border-left:3px solid #64b5f6;font-size:14px;font-weight:500">'+rx.first+'</div>';
  h+='</div>';

  h+='<div style="margin-bottom:10px">';
  h+='<div style="font-size:12px;color:#ffb74d;font-weight:600;margin-bottom:8px">ESCALATION IF LDL-C NOT AT TARGET</div>';
  h+='<div style="padding:10px;background:rgba(255,183,77,0.08);border-radius:8px;border-left:3px solid #ffb74d;font-size:13px">'+rx.escalate+'</div>';
  h+='</div>';

  h+='<div style="margin-bottom:10px">';
  h+='<div style="font-size:12px;color:#ce93d8;font-weight:600;margin-bottom:8px">STRUCTURAL RATIONALE</div>';
  h+='<div style="padding:10px;background:rgba(206,147,216,0.06);border-radius:8px;font-size:12px;line-height:1.5;color:#ccc">'+rx.rationale+'</div>';
  h+='</div>';

  if(rx.caution){
    h+='<div style="padding:8px 10px;background:rgba(239,83,80,0.08);border-radius:8px;border-left:3px solid #ef5350;font-size:12px;color:#ef9a9a">';
    h+='<b>\u26A0 Clinical note:</b> '+rx.caution+'</div>';
  }

  h+='</div>';
  return h;
}

function getVUSPanel(d){
  var hasVUS=false,vusVariants=[];
  if(d.variants){for(var i=0;i<d.variants.length;i++){if(d.variants[i].vus_reclass){hasVUS=true;vusVariants.push(d.variants[i]);}}}
  if(!hasVUS)return '';
  var vus=DATA.vus;
  var h='<div class="vus-banner"><h3>VUS Reclassification &mdash; SSS-Based Structural Assessment</h3>';
  for(var i=0;i<vusVariants.length;i++){
    var v=vusVariants[i];
    var cls=v.vus_reclass==='likely_pathogenic'?'pathogenic':v.vus_reclass==='likely_benign'?'benign':'uncertain';
    var label=v.vus_reclass==='likely_pathogenic'?'SSS suggests LIKELY PATHOGENIC':v.vus_reclass==='likely_benign'?'SSS suggests LIKELY BENIGN':'SSS INDETERMINATE \u2014 remains VUS';
    h+='<div style="margin-bottom:8px"><span class="notation">'+v.id+'</span> '+getClinvarTag(v.clinvar)+' (SSS = '+(v.sss!==null?v.sss.toFixed(3):'N/A')+')<div class="reclass '+cls+'">\u27A1 '+label+'</div></div>';
  }
  if(vus&&vus.high_sss){
    h+='<div style="margin-top:10px;font-size:12px;color:#b0bec5"><b>Population validation (UK Biobank):</b><br>';
    h+='VUS with high SSS (n='+vus.high_sss.n+'): ASCVD '+vus.high_sss.ascvd_rate+'%, LDL-C '+vus.high_sss.ldl_mean+' mmol/L, Simon Broome definite '+vus.high_sss.simon_broome_definite_rate+'%<br>';
    if(vus.low_sss)h+='VUS with low SSS (n='+vus.low_sss.n+'): ASCVD '+vus.low_sss.ascvd_rate+'%, LDL-C '+vus.low_sss.ldl_mean+' mmol/L</div>';
  }
  if(vus&&vus.prediction){
    h+='<div class="vus-stats">';
    h+='<div class="vus-stat"><div class="num">'+vus.prediction.sensitivity+'%</div><div class="lab">Sensitivity</div></div>';
    h+='<div class="vus-stat"><div class="num">'+vus.prediction.specificity+'%</div><div class="lab">Specificity</div></div>';
    h+='<div class="vus-stat"><div class="num">'+vus.prediction.ppv+'%</div><div class="lab">PPV</div></div>';
    h+='<div class="vus-stat"><div class="num">'+vus.prediction.npv+'%</div><div class="lab">NPV</div></div>';
    h+='</div>';
  }
  h+='<div class="src">Reclassification threshold: SSS &gt;0.5 = likely pathogenic, &lt;0.3 = likely benign. Derived from 232 VUS carriers in UK Biobank.</div>';
  h+='</div>';
  return h;
}

function getIslamPanel(d){
  if(!d.islam||d.islam.length===0)return '';
  var h='<div class="islam-box"><h3>Independent Functional Validation &mdash; Islam et al. Cell-Based Assay</h3>';
  for(var i=0;i<d.islam.length;i++){
    var iv=d.islam[i];
    h+='<div class="islam-row"><div style="flex:1"><b>'+iv.variant+'</b>';
    if(iv.clinvar)h+=' &mdash; ClinVar: '+getClinvarTag(iv.clinvar);
    h+='</div><div>'+getFuncBadge(iv.functional_group)+'</div>';
    if(iv.activity_pct!==null)h+='<div style="width:80px;text-align:right;font-size:13px;font-weight:bold;color:'+(iv.activity_pct<30?'#ff5252':iv.activity_pct<70?'#ff9800':'#4caf50')+'">'+iv.activity_pct+'% activity</div>';
    h+='</div>';
  }
  h+='<div class="src" style="margin-top:8px">Islam et al. (2024): Cell-surface LDLR activity assay for 234 variants. Functional groups: Loss-of-function (&lt;15%), Defective (15\u201340%), Mildly-defective (40\u201370%), Non-defective (&gt;70%).</div>';
  h+='</div>';
  return h;
}

function get3DViewer(pos){
  return '<div class="panel"><h3>3D Structure &mdash; AlphaFold Predicted (Position '+pos+')</h3><div class="viewer-container" id="mol-viewer"><div class="viewer-overlay" id="mol-overlay" onclick="initMolstar('+pos+')"><span>Click to load AlphaFold structure and highlight residue '+pos+'</span></div></div><div style="font-size:10px;color:#666;margin-top:3px">Human LDLR (UniProt P01130). AlphaFold structure AF-P01130-F1, coloured by pLDDT confidence. Selected residue highlighted. Powered by PDBe Mol*.</div></div>';
}

function getMolExplanation(d){
  var pos=d.amino_acid_pos,aa=d.wt_aa,dom=d.domain,mc=d.mech_class;
  var bs=d.binding_site,contact=d.ldl_contact,ddg=d.max_ddg,pen=d.penetrance_60;

  // LOCATION
  var loc='';
  if(mc==='binding_failure'){
    loc='Position '+pos+' lies within the <b>'+dom+'</b> domain, a component of <b>Binding Site 1</b> (BS1). BS1 comprises LA3\u2013LA7 and EGF-A (residues 107\u2013353) and constitutes the primary LDL-capture interface, directly engaging the LDL receptor-binding domain of apolipoprotein B-100.';
    if(dom==='EGF-A')loc='Position '+pos+' lies within the <b>EGF-A domain</b> (residues 314\u2013353), which serves a dual function: it contributes to BS1 LDL binding AND contains the <b>PCSK9 binding epitope</b>. Mutations in EGF-A may simultaneously impair LDL clearance and reduce the efficacy of PCSK9-targeted therapies (evolocumab, alirocumab).';
  }else if(mc==='recycling_failure'){
    loc='Position '+pos+' is one of <b>4 pH-switch histidine residues</b> (H211, H285, H583, H607) essential for receptor recycling. At endosomal pH 5.5\u20136.0, protonation of these histidines drives a conformational change in the beta-propeller domain, causing it to fold over the ligand-binding repeats and displace bound LDL. Without this switch, the receptor\u2013LDL complex cannot dissociate and both are targeted for lysosomal degradation.';
  }else if(mc==='recycling_pivot'){
    loc='Position '+pos+' is in the <b>EGF-B domain</b> (residues 354\u2013393), which functions as the mechanical hinge enabling the pH-dependent conformational change required for LDL release and receptor recycling.';
  }else if(mc==='interface_BS2'){
    loc='Position '+pos+' is in the <b>beta-propeller domain</b> (residues 394\u2013663) at <b>Binding Site 2</b>. Cryo-EM (Reimund et al., Nature 2025) confirmed BS2 contacts the N-terminal domain of apoB-100 on LDL. This domain also mediates receptor recycling by folding over the ligand-binding repeats at acidic pH.';
  }else if(dom==='LA1'||dom==='LA2'){
    loc='Position '+pos+' is in <b>'+dom+'</b>. Importantly, cryo-EM (Reimund et al., 2025) demonstrated that LA1\u2013LA2 do <b>not directly contact LDL particles</b>. Pathogenicity here is driven by <b>protein misfolding and ER retention</b> (Class 2 defect), not by impaired LDL binding.';
  }else if(dom==='Cytoplasmic'){
    loc='Position '+pos+' is in the <b>cytoplasmic tail</b> (residues 771\u2013860), containing the NPXY internalisation signal (Class 5 defect in Hobbs classification). The receptor reaches the cell surface and binds LDL but cannot be recruited to clathrin-coated pits for endocytosis.';
  }else if(dom==='Transmembrane'){
    loc='Position '+pos+' is in the <b>transmembrane helix</b> (residues 750\u2013770), anchoring the receptor in the hepatocyte plasma membrane. Mutations may prevent membrane insertion or alter receptor orientation.';
  }else if(dom==='O-linked sugars'){
    loc='Position '+pos+' is in the <b>O-linked sugar domain</b> (residues 713\u2013749), a glycosylated stalk that elevates the extracellular domains above the glycocalyx and provides protection against metalloproteinase cleavage.';
  }else{
    loc='Position '+pos+' is in the <b>'+dom+'</b> domain of the LDL receptor (UniProt P01130).';
  }

  // STRUCTURAL CONSEQUENCE
  var struc='';
  if(ddg!==null&&ddg>4){
    struc='FoldX 5.0 analysis predicts <b>severe thermodynamic destabilisation</b> (\u0394\u0394G = '+ddg.toFixed(1)+' kcal/mol, threshold &gt;2 kcal/mol = pathogenic). This magnitude of destabilisation is consistent with <b>ER-associated degradation</b> (ERAD) of the misfolded receptor \u2014 functionally equivalent to a null allele. The mutant receptor is unlikely to reach the hepatocyte cell surface.';
  }else if(ddg!==null&&ddg>2){
    struc='FoldX predicts <b>significant destabilisation</b> (\u0394\u0394G = '+ddg.toFixed(1)+' kcal/mol). This likely impairs folding efficiency in the endoplasmic reticulum, resulting in reduced surface receptor density. Some functional receptor may still be expressed (partial loss-of-function).';
  }else if(ddg!==null&&ddg>1){
    struc='FoldX predicts <b>mild destabilisation</b> (\u0394\u0394G = '+ddg.toFixed(1)+' kcal/mol). The receptor is likely to fold and traffic to the cell surface, though with potentially altered local stability or dynamics.';
  }else if(contact&&bs==='BS1'){
    struc='This is a <b>direct LDL-contact residue</b> in BS1. Even without large-scale protein destabilisation, mutations here can disrupt the electrostatic and van der Waals interactions with apoB-100 on the LDL particle, directly reducing <b>binding affinity and LDL clearance rate</b>.';
  }else if(aa==='C'){
    struc='Wildtype residue is <b>cysteine</b>, forming a disulfide bond essential for domain structural integrity. In the ligand-binding repeats, each ~40-residue module is stabilised by 3 disulfide bonds coordinating a calcium ion. Loss of any cysteine typically causes severe misfolding (Class 2 defect).';
  }else{
    struc='The structural impact at this position depends on the specific amino acid substitution. The wildtype residue ('+aa+') contributes to local protein stability and packing interactions.';
  }

  // CLINICAL
  var clin='';
  if(pen!==null&&pen>30){
    clin='In the UK Biobank whole-exome sequencing cohort, carriers of variants in this domain show <b>'+pen.toFixed(1)+'% ASCVD penetrance by age 60</b>. This represents one of the highest-risk domains and supports early initiation of intensive lipid-lowering therapy, ideally before age 40.';
  }else if(pen!==null&&pen>20){
    clin='Domain-level penetrance: <b>'+pen.toFixed(1)+'% ASCVD by age 60</b> in UK Biobank carriers. This is substantially elevated above general population risk and supports combination therapy (high-intensity statin + ezetimibe \u00B1 PCSK9 inhibitor).';
  }else if(pen!==null&&pen>10){
    clin='Domain-level penetrance: <b>'+pen.toFixed(1)+'% ASCVD by age 60</b>. While moderate, this exceeds general population risk several-fold. Statin therapy with monitoring and dose optimisation is recommended.';
  }else{
    clin='Limited domain-level penetrance data for this region. Clinical risk should be assessed using family history, lipid phenotype (LDL-C, Lp(a)), coronary artery calcium scoring, and CTCA where indicated.';
  }

  // TREATMENT
  var treat='';
  if(mc==='binding_failure'&&dom!=='EGF-A'){
    treat='<b>Therapeutic consideration:</b> If this represents a null or severe loss-of-function allele, upregulation of LDLR expression via statins has limited benefit from the mutant allele (the wild-type allele still responds). PCSK9 inhibitors protect the remaining functional receptors from degradation. For homozygous or compound heterozygous null: <b>evinacumab</b> (anti-ANGPTL3, LDLR-independent mechanism, ELIPSE HoFH trial) or lipoprotein apheresis.';
  }else if(dom==='EGF-A'){
    treat='<b>EGF-A domain mutations require special therapeutic consideration.</b> PCSK9 binds to EGF-A residues 314\u2013353 to promote receptor degradation. Mutations here may alter PCSK9 binding kinetics, potentially <b>reducing the efficacy of PCSK9 monoclonal antibodies</b> (evolocumab, alirocumab). Alternative: <b>inclisiran</b> (PCSK9 siRNA \u2014 reduces PCSK9 production rather than blocking binding) or <b>evinacumab</b>.';
  }else if(mc==='recycling_failure'){
    treat='<b>Recycling-deficient receptors</b> bind and internalise LDL normally but are degraded after a single cycle (vs ~100 cycles for normal receptors). This dramatically reduces receptor half-life. Statins partially compensate by driving new receptor synthesis. PCSK9 inhibitors may have <b>reduced incremental benefit</b> because the primary problem is pH-switch failure, not PCSK9-mediated degradation.';
  }else if(dom==='Cytoplasmic'){
    treat='<b>Class 5 internalisation defect.</b> Receptor reaches the surface and binds LDL but cannot be endocytosed. Neither statins (upregulate expression of a non-internalising receptor) nor PCSK9 inhibitors (prevent degradation of a receptor that is already surface-stable) fully address the defect. Consider <b>evinacumab</b> or <b>apheresis</b>.';
  }else if(dom==='Beta-propeller'){
    treat='<b>Beta-propeller domain: best documented treatment response.</b> In our cohort, LDL-C decreased from 5.06 to 3.41 mmol/L on lipid-lowering therapy (P=0.001, n=63 treated). PCSK9 inhibitors are expected to be effective because the PCSK9 binding site (EGF-A) is structurally intact.';
  }else{
    treat='<b>Standard approach:</b> High-intensity statin (atorvastatin 40\u201380mg or rosuvastatin 20\u201340mg) + ezetimibe 10mg as first-line. If LDL-C remains &gt;2.6 mmol/L (or &gt;1.8 mmol/L in very high risk), add PCSK9 inhibitor. Monitor LDL-C response at 6\u20138 weeks and adjust.';
  }

  // Jensson
  var jensson='';
  if(d.jensson_years!==null){
    jensson='<div class="expl-step"><div class="expl-num">6</div><div class="expl-text"><b>Population lifespan impact (Jensson et al., NEJM 2023)</b><br>In the Icelandic deCODE study (n=166,281), carriers of a variant at this position lived <b>'+Math.abs(d.jensson_years).toFixed(1)+' years '+(d.jensson_years<0?'shorter':'longer')+'</b> than non-carriers'+(d.jensson_sig?' (P&lt;0.05)':' (not significant at current sample size)')+'.<br><span class="src">Jensson et al. N Engl J Med 2023;389:1741\u20131752.</span></div></div>';
  }

  return '<div class="expl-box"><h3 style="color:#64b5f6;margin-bottom:12px">Clinical Molecular Summary</h3><div class="expl-pathway"><div class="expl-step"><div class="expl-num">1</div><div class="expl-text"><b>Structural location</b><br>'+loc+'<br><span class="src">Reimund et al. Nature 2025; 638:829\u2013836.</span></div></div><div class="expl-step"><div class="expl-num">2</div><div class="expl-text"><b>Predicted structural consequence</b><br>'+struc+'</div></div><div class="expl-step"><div class="expl-num">3</div><div class="expl-text"><b>Pathophysiology</b><br>'+(mc==='binding_failure'?'Impaired LDL capture at the hepatocyte surface \u2192 reduced hepatic LDL clearance \u2192 elevated circulating LDL-C \u2192 accelerated subendothelial retention and oxidation \u2192 foam cell formation \u2192 atherosclerotic plaque progression.':mc==='recycling_failure'?'Receptor\u2013LDL complex fails to dissociate in the endosome \u2192 both receptor and LDL diverted to lysosomal degradation \u2192 receptor half-life reduced from ~20 hours to &lt;1 hour \u2192 progressive loss of surface receptor density \u2192 impaired LDL clearance.':mc==='recycling_pivot'?'Impaired mechanical hinge for pH-dependent conformational change \u2192 reduced recycling efficiency \u2192 accelerated receptor turnover \u2192 fewer functional receptors available for LDL clearance.':'Misfolded receptor retained in endoplasmic reticulum \u2192 ER-associated degradation (ERAD) \u2192 reduced surface receptor expression \u2192 impaired hepatic LDL clearance \u2192 elevated circulating LDL-C.')+'<br><span class="src">Boren et al. Nat Rev Cardiol 2025; Brown &amp; Goldstein, Nobel Lecture 1985.</span></div></div><div class="expl-step"><div class="expl-num">4</div><div class="expl-text"><b>Clinical significance</b><br>'+clin+'<br><span class="src">UK Biobank WES cohort + Wales FH Registry + DRAGON3.</span></div></div><div class="expl-step"><div class="expl-num">5</div><div class="expl-text"><b>Therapeutic implications</b><br>'+treat+'</div></div>'+jensson+'</div></div>';
}

function showPosition(pos){
  var d=DATA.atlas[pos];
  if(!d)return;
  var card=document.getElementById('resultCard');
  card.className='result active';
  document.getElementById('searchResults').innerHTML='';
  _molReady=false;
  var html='';

  // Traffic Light
  html+=getTrafficLightHTML(d);

  // Title + Tags
  html+='<div class="card"><h3>Position '+pos+': '+d.wt_aa+' \u2014 '+d.domain+'</h3>';
  html+='<div style="margin-bottom:10px"><span class="notation" style="font-size:16px">'+d.p_notation+'</span> &nbsp; <span class="notation" style="font-size:16px">'+d.c_notation+'</span></div>';
  html+='<div style="margin-bottom:10px">'+getMechTag(d.mech_class)+' ';
  if(d.binding_site==='BS1')html+='<span class="tag bs1">BS1</span> ';
  else if(d.binding_site==='BS2')html+='<span class="tag bs2">BS2</span> ';
  html+=(d.ldl_contact?'<span class="tag contact">LDL Contact</span>':'<span class="tag no-contact">No LDL Contact</span>');
  if(d.at_pcsk9)html+=' <span class="tag bs1">PCSK9 Interface</span>';
  html+='</div></div>';

  // 3D Viewer
  html+=get3DViewer(pos);

  // Penetrance
  html+=getPenetrancePanel(d);

  // LLT Response
  html+=getLLTPanel(d);

  // VUS Reclassification
  html+=getVUSPanel(d);

  // Islam Validation
  html+=getIslamPanel(d);

  // NEW: Tabet DMS Panel
  html+=getTabetPanel(d);

  // NEW: AlphaGenome Panel
  html+=getAlphaGenomePanel(d);

  // Structural Details
  html+='<div class="panel"><h3>Structural and Functional Data</h3><div class="card-grid"><div>';
  html+='<div class="field"><div class="label">Domain (Reimund 2025)</div><div class="val hi">'+d.domain+'</div></div>';
  html+='<div class="field"><div class="label">Mechanistic Class</div><div class="val">'+(MECH_LABELS[d.mech_class]||d.mech_class||'N/A')+'</div></div>';
  html+='<div class="field"><div class="label">PCSK9 Interface</div><div class="val">'+(d.at_pcsk9?'Yes ('+d.pcsk9_dist+' \u00C5)':'No ('+d.pcsk9_dist+' \u00C5)')+'</div></div>';
  html+='<div class="field"><div class="label">AlphaFold pLDDT</div><div class="val">'+(d.plddt||'-')+'</div></div>';
  html+='</div><div>';
  html+='<div class="field"><div class="label">FoldX Max \u0394\u0394G</div><div class="val '+(d.max_ddg>2?'bad':d.max_ddg>1?'warn':'good')+'">'+(d.max_ddg!==null?d.max_ddg+' kcal/mol':'Not computed')+'</div></div>';
  html+='<div class="field"><div class="label">Rosetta Energy</div><div class="val">'+(d.rosetta!==null?d.rosetta+' REU':'-')+'</div></div>';
  html+='<div class="field"><div class="label">Islam et al. Activity</div><div class="val">'+(d.islam_activity!==null?d.islam_activity+'%':'-')+'</div></div>';
  html+='<div class="field"><div class="label">Risk Classification</div><div class="val '+(d.risk==="Critical"?"bad":d.risk==="High"?"warn":"good")+'">'+d.risk+'</div></div>';
  if(d.sat_mean_ddG!==null)html+='<div class="field"><div class="label">Saturation Mutagenesis (mean)</div><div class="val">'+d.sat_mean_ddG+' kcal/mol</div></div>';
  if(d.sat_worst_mut)html+='<div class="field"><div class="label">Worst Mutation</div><div class="val bad">'+d.sat_worst_mut+' ('+d.sat_worst_ddG+')</div></div>';
  html+='</div><div>';
  if(d.jensson_years!==null)html+='<div class="field"><div class="label">Jensson NEJM 2023 \u2014 Lifespan</div><div class="val '+(d.jensson_sig?'bad':'warn')+'">'+d.jensson_years.toFixed(2)+' years '+(d.jensson_sig?'(P&lt;0.05)':'(NS)')+'</div></div>';
  if(d.reimund_site)html+='<div class="field"><div class="label">Reimund Nature 2025 \u2014 Cryo-EM</div><div class="val hi">'+d.reimund_site+'</div><div class="val">'+(d.reimund_func||'N.D.')+'</div></div>';
  html+='<div class="field"><div class="label">Drug Target</div><div class="val">'+d.drug_target+'</div></div>';
  html+='</div></div></div>';

  // Variant Table
  if(d.variants&&d.variants.length>0){
    html+='<div class="panel"><h3>Known Variants at Position '+pos+'</h3><table class="var-table"><thead><tr><th>Variant</th><th>ClinVar</th><th>SSS</th><th>\u0394\u0394G</th><th>Type</th><th>LLT Data</th><th>VUS Reclass</th></tr></thead><tbody>';
    for(var i=0;i<d.variants.length;i++){
      var v=d.variants[i];
      var reclassCell='-';
      if(v.vus_reclass==='likely_pathogenic')reclassCell='<span style="color:#ff5252;font-weight:bold">\u2191 Likely Path.</span>';
      else if(v.vus_reclass==='likely_benign')reclassCell='<span style="color:#69f0ae;font-weight:bold">\u2193 Likely Benign</span>';
      else if(v.vus_reclass==='uncertain')reclassCell='<span style="color:#ffc400">Uncertain</span>';
      var lltCell='-';
      if(v.llt)lltCell='<span style="color:#4caf50">'+v.llt.ldl_pct_reduction.toFixed(0)+'% reduction (n='+v.llt.n_treated+')</span>';
      html+='<tr><td><span class="notation">'+v.id+'</span></td><td>'+getClinvarTag(v.clinvar)+'</td><td>'+(v.sss!==null?v.sss.toFixed(3):'-')+'</td><td>'+(v.ddG!==null?v.ddG.toFixed(1):'-')+'</td><td>'+(v.type||'-')+'</td><td>'+lltCell+'</td><td>'+reclassCell+'</td></tr>';
    }
    html+='</tbody></table></div>';
  }

  // Molecular Explanation
  html+=getMolExplanation(d);

  card.innerHTML=html;
  if(!_touchMoved)window.scrollTo({top:card.offsetTop-20,behavior:'smooth'});
}

function doSearch(query){
  var q=query.trim().toLowerCase();
  if(!q){document.getElementById('resultCard').className='result';document.getElementById('searchResults').innerHTML='';return;}
  var pMatch=q.match(/^p\.?([a-z]{3})(\d+)/i);
  if(pMatch){var pp=pMatch[2];if(DATA.atlas[pp]){showPosition(pp);return;}}
  var cMatch=q.match(/^c\.?(\d+)/i);
  if(cMatch){var cpos=parseInt(cMatch[1]);var aaPos=Math.floor((cpos-1)/3)+1;if(DATA.atlas[String(aaPos)]){showPosition(String(aaPos));return;}}
  if(DATA.atlas[q]&&q.match(/^\d+$/)){showPosition(q);return;}
  var results=[];
  for(var pos in DATA.atlas){
    var d=DATA.atlas[pos];
    var searchStr=[d.domain,d.mech_class,d.binding_site,d.wt_aa,d.p_notation,d.c_notation,d.risk,d.traffic_light,pos].filter(Boolean).join(' ').toLowerCase();
    if(d.variants){for(var i=0;i<d.variants.length;i++){var cv=d.variants[i].clinvar;if(cv)searchStr+=' '+cv.toLowerCase();var vt=d.variants[i].type;if(vt)searchStr+=' '+vt.toLowerCase();if(d.variants[i].vus_reclass)searchStr+=' vus';}}
    if(searchStr.indexOf(q)!==-1)results.push({pos:pos,d:d});
  }
  if(results.length===1){showPosition(results[0].pos);return;}
  if(results.length>0){
    document.getElementById('resultCard').className='result';
    var html='<div class="card"><h3>'+results.length+' positions matching "'+query+'"</h3><table class="var-table"><thead><tr><th>Pos</th><th>p.</th><th>AA</th><th>Domain</th><th>Risk</th><th>Mechanism</th><th>Penetrance</th><th>\u0394\u0394G</th></tr></thead><tbody>';
    for(var i=0;i<Math.min(results.length,200);i++){
      var r=results[i];
      var tlDot=r.d.traffic_light==='red'?'\u{1F534}':r.d.traffic_light==='amber'?'\u{1F7E0}':'\u{1F7E2}';
      html+='<tr style="cursor:pointer" onclick="if(!_touchMoved)showPosition(\''+r.pos+'\')"><td>'+r.pos+'</td><td><span class="notation">'+r.d.p_notation+'</span></td><td>'+r.d.wt_aa+'</td><td>'+r.d.domain+'</td><td>'+tlDot+' '+r.d.risk+'</td><td>'+getMechTag(r.d.mech_class)+'</td><td>'+(r.d.penetrance_60?r.d.penetrance_60.toFixed(1)+'%':'-')+'</td><td>'+(r.d.max_ddg!==null?r.d.max_ddg:'-')+'</td></tr>';
    }
    html+='</tbody></table></div>';
    document.getElementById('searchResults').innerHTML=html;
  }else{
    document.getElementById('resultCard').className='result';
    document.getElementById('searchResults').innerHTML='<div class="card"><p style="color:#888">No results for "'+query+'"</p></div>';
  }
}

var _searchTimer=null;
document.getElementById('searchInput').addEventListener('input',function(e){
  clearTimeout(_searchTimer);_searchTimer=setTimeout(function(){doSearch(e.target.value)},400);
});
document.getElementById('searchInput').addEventListener('keydown',function(e){
  if(e.key==='Enter'){clearTimeout(_searchTimer);doSearch(e.target.value);}
});
buildDomainBar();
</script>
</body>
</html>'''

# Insert data
html_final = html_template.replace('__DATA_PLACEHOLDER__', data_json)

output_dir = os.path.join(BASE, "website")
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, "index.html")
with open(output_path, 'w', encoding='utf-8') as f:
    f.write(html_final)

print(f"\n[5/5] Done!")
print(f"   Website: {len(html_final):,} bytes")
print(f"   Output: {output_path}")
print(f"   Data sources integrated:")
print(f"     - ClinVar (external_vus_validation): {n_clinvar} variants classified")
print(f"     - LLT response (DRAGON3/Wales): {n_llt} variants with treatment data")
print(f"     - Islam et al. functional: {n_islam} positions validated")
print(f"     - UK Biobank KM penetrance: {len(km_domain)} domains")
print(f"     - FoldX, Rosetta, AF3, Saturation: per-position")

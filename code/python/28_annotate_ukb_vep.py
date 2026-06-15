#!/usr/bin/env python3
"""
Annotate UKB LDLR exonic variants via Ensembl VEP REST API
Then map to SSS catalogue for variant-level external validation
"""
import csv, json, time, requests, os
from collections import Counter, defaultdict

BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"
AF = f"{BASE}/alphafold"

def sf(x):
    try: return float(str(x).strip())
    except: return None

# LDLR exon boundaries (GRCh38)
LDLR_EXONS = [
    (11133306, 11133820, 1), (11131193, 11131388, 2), (11129528, 11129667, 3),
    (11128103, 11128254, 4), (11120108, 11120346, 5), (11116921, 11117028, 6),
    (11116779, 11116899, 7), (11116248, 11116394, 8), (11113744, 11113862, 9),
    (11113462, 11113607, 10), (11110682, 11110879, 11), (11107282, 11107564, 12),
    (11105579, 11105752, 13), (11105398, 11105509, 14), (11100051, 11100340, 15),
    (11099365, 11099476, 16), (11098479, 11098575, 17), (11089362, 11091039, 18),
]

# ═══════════════════════════════════════════════════════════════
# STEP 1: Get unique exonic variants
# ═══════════════════════════════════════════════════════════════
print("Step 1: Loading rare carriers...", flush=True)
rare = list(csv.DictReader(open(f'{AF}/ukb_fh_rare_carriers.tsv', encoding='utf-8-sig'), delimiter='\t'))

unique_variants = {}
for r in rare:
    if r['GENE'] != 'LDLR': continue
    pos = int(r['POS'])
    if not any(s <= pos <= e for s, e, _ in LDLR_EXONS): continue
    key = (r['CHROM'], r['POS'], r['REF'], r['ALT'])
    if key not in unique_variants:
        unique_variants[key] = 0
    unique_variants[key] += 1

print(f"  {len(unique_variants)} unique exonic LDLR variants to annotate", flush=True)

# ═══════════════════════════════════════════════════════════════
# STEP 2: VEP REST API annotation
# ═══════════════════════════════════════════════════════════════
print("\nStep 2: Annotating via Ensembl VEP REST API...", flush=True)

url = 'https://rest.ensembl.org/vep/homo_sapiens/region'
headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}

variant_list = list(unique_variants.keys())
annotations = {}

for i in range(0, len(variant_list), 200):
    batch = variant_list[i:i+200]
    vep_input = []
    for chrom, pos, ref, alt in batch:
        p = int(pos)
        if len(ref) == 1 and len(alt) == 1:
            vep_input.append(f"{chrom} {p} {p} {ref}/{alt} 1")
        elif len(ref) > len(alt):
            # Deletion
            del_start = p + 1
            del_end = p + len(ref) - 1
            vep_input.append(f"{chrom} {del_start} {del_end} {ref}/{alt} 1")
        else:
            # Insertion
            vep_input.append(f"{chrom} {p} {p} {ref}/{alt} 1")

    batch_num = i // 200 + 1
    total_batches = (len(variant_list) + 199) // 200
    print(f"  Batch {batch_num}/{total_batches} ({len(batch)} variants)...", end=' ', flush=True)

    retries = 0
    while retries < 3:
        try:
            resp = requests.post(url, headers=headers,
                               data=json.dumps({"variants": vep_input}),
                               timeout=120)
            if resp.status_code == 200:
                results = resp.json()
                for r in results:
                    for tc in r.get('transcript_consequences', []):
                        if tc.get('gene_symbol') == 'LDLR' and tc.get('biotype') == 'protein_coding':
                            hgvsc = tc.get('hgvsc', '').split(':')[-1] if tc.get('hgvsc') else ''
                            hgvsp = tc.get('hgvsp', '').split(':')[-1] if tc.get('hgvsp') else ''
                            inp = r.get('input', '')
                            parts = inp.split()
                            if len(parts) >= 4:
                                ref_alt = parts[3].split('/')
                                gkey = (parts[0], parts[1], ref_alt[0], ref_alt[1] if len(ref_alt) > 1 else '')
                                annotations[gkey] = {
                                    'cdna': hgvsc,
                                    'protein': hgvsp,
                                    'consequence': ','.join(tc.get('consequence_terms', [])),
                                    'impact': tc.get('impact', ''),
                                    'sift': tc.get('sift_prediction', ''),
                                    'polyphen': tc.get('polyphen_prediction', ''),
                                    'protein_position': str(tc.get('protein_start', '')),
                                    'amino_acids': tc.get('amino_acids', ''),
                                }
                            break
                print(f"{len(results)} annotated", flush=True)
                break
            elif resp.status_code == 429:
                print("rate limited, waiting 10s...", end=' ', flush=True)
                time.sleep(10)
                retries += 1
            else:
                print(f"HTTP {resp.status_code}", flush=True)
                break
        except Exception as e:
            print(f"Error: {e}", flush=True)
            time.sleep(5)
            retries += 1

    time.sleep(1)

print(f"\nTotal annotated: {len(annotations)}", flush=True)

# ═══════════════════════════════════════════════════════════════
# STEP 3: Save annotations
# ═══════════════════════════════════════════════════════════════
print("\nStep 3: Saving annotations...", flush=True)

with open(f'{AF}/analysis/ukb_ldlr_vep_annotations.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['chrom', 'pos', 'ref', 'alt', 'n_carriers', 'cdna', 'protein',
                'consequence', 'impact', 'sift', 'polyphen', 'protein_position', 'amino_acids'])
    for (chrom, pos, ref, alt), count in unique_variants.items():
        ann = annotations.get((chrom, pos, ref, alt), {})
        w.writerow([chrom, pos, ref, alt, count,
                   ann.get('cdna', ''), ann.get('protein', ''),
                   ann.get('consequence', ''), ann.get('impact', ''),
                   ann.get('sift', ''), ann.get('polyphen', ''),
                   ann.get('protein_position', ''), ann.get('amino_acids', '')])

print(f"  Saved: {AF}/analysis/ukb_ldlr_vep_annotations.csv")

# ═══════════════════════════════════════════════════════════════
# STEP 4: Consequence summary
# ═══════════════════════════════════════════════════════════════
cons = Counter()
for ann in annotations.values():
    for c in ann['consequence'].split(','):
        if c: cons[c] += 1

print(f"\nConsequence types:")
for c, n in cons.most_common():
    print(f"  {c}: {n}")

# FH-causing consequence types
pathogenic_types = {'missense_variant', 'stop_gained', 'frameshift_variant',
                   'splice_donor_variant', 'splice_acceptor_variant',
                   'start_lost', 'stop_lost', 'inframe_deletion', 'inframe_insertion'}

fh_annotations = {k: v for k, v in annotations.items()
                  if any(c in pathogenic_types for c in v['consequence'].split(','))}

print(f"\nFH-causing variants (protein-altering): {len(fh_annotations)}")

# ═══════════════════════════════════════════════════════════════
# STEP 5: Map to SSS catalogue
# ═══════════════════════════════════════════════════════════════
print("\nStep 5: Mapping to SSS catalogue...", flush=True)

sss_lookup = {}
sss_scores = list(csv.DictReader(open(f'{AF}/analysis/structural_severity_scores.csv', encoding='utf-8-sig')))
for r in sss_scores:
    vid = r.get('variant_id', '').strip()
    if vid:
        sss_lookup[vid] = {
            'sss': sf(r.get('sss')),
            'ddG': sf(r.get('ddG')),
            'domain': r.get('domain', ''),
            'variant_type': r.get('variant_type', ''),
        }

# Also try full SSS file
try:
    for r in csv.DictReader(open(f'{AF}/analysis/structural_severity_scores_full.csv', encoding='utf-8-sig')):
        vid = r.get('variant_id', '').strip()
        if vid and vid not in sss_lookup:
            sss_lookup[vid] = {
                'sss': sf(r.get('sss')),
                'ddG': sf(r.get('ddG')),
                'domain': r.get('domain', ''),
                'variant_type': r.get('variant_type', ''),
            }
except:
    pass

print(f"  SSS catalogue: {len(sss_lookup)} variants")

# Match VEP annotations to SSS
matched = 0
unmatched = 0
matched_variants = {}

for gkey, ann in fh_annotations.items():
    cdna = ann.get('cdna', '')
    variant_id = f"LDLR:{cdna}" if cdna else ""

    sss_info = sss_lookup.get(variant_id)
    if sss_info:
        matched += 1
        matched_variants[gkey] = {
            'variant_id': variant_id,
            'n_carriers': unique_variants.get(gkey, 0),
            **ann,
            **sss_info,
        }
    else:
        unmatched += 1

print(f"  Matched to SSS: {matched}")
print(f"  Unmatched: {unmatched}")

if matched_variants:
    print(f"\n  Matched FH variants with SSS:")
    for gkey, info in sorted(matched_variants.items(), key=lambda x: -x[1]['n_carriers'])[:30]:
        print(f"    {info['variant_id']}: {info['n_carriers']} UKB carriers, SSS={info['sss']}, "
              f"ddG={info['ddG']}, domain={info['domain']}, {info['consequence']}")

# ═══════════════════════════════════════════════════════════════
# STEP 6: Build carrier-level file with SSS for external validation
# ═══════════════════════════════════════════════════════════════
print("\nStep 6: Building carrier-level SSS file...", flush=True)

# Build genomic-to-SSS mapping
gpos_to_sss = {}
for gkey, info in matched_variants.items():
    gpos_to_sss[gkey] = info

# Also map all VEP-annotated variants (not just matched to SSS)
gpos_to_ann = {}
for gkey, ann in annotations.items():
    gpos_to_ann[gkey] = ann

# Process rare carriers
carrier_sss = []
for r in rare:
    if r['GENE'] != 'LDLR': continue
    gkey = (r['CHROM'], r['POS'], r['REF'], r['ALT'])
    ann = gpos_to_ann.get(gkey, {})
    sss_info = gpos_to_sss.get(gkey, {})

    carrier_sss.append({
        'eid': r['EID'],
        'gene': r['GENE'],
        'chrom': r['CHROM'],
        'pos': r['POS'],
        'ref': r['REF'],
        'alt': r['ALT'],
        'gt': r['GT'],
        'cdna': ann.get('cdna', ''),
        'protein': ann.get('protein', ''),
        'consequence': ann.get('consequence', ''),
        'impact': ann.get('impact', ''),
        'sift': ann.get('sift', ''),
        'polyphen': ann.get('polyphen', ''),
        'variant_id': sss_info.get('variant_id', ''),
        'sss': sss_info.get('sss', ''),
        'ddG': sss_info.get('ddG', ''),
        'domain': sss_info.get('domain', ''),
    })

with open(f'{AF}/analysis/ukb_carriers_with_sss.csv', 'w', newline='') as f:
    fields = ['eid', 'gene', 'chrom', 'pos', 'ref', 'alt', 'gt',
              'cdna', 'protein', 'consequence', 'impact', 'sift', 'polyphen',
              'variant_id', 'sss', 'ddG', 'domain']
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    for c in carrier_sss:
        w.writerow(c)

carriers_with_sss = sum(1 for c in carrier_sss if c['sss'])
print(f"  Total carrier records: {len(carrier_sss)}")
print(f"  Carriers with SSS: {carriers_with_sss}")
print(f"  Saved: {AF}/analysis/ukb_carriers_with_sss.csv")

print("\n" + "=" * 60)
print("  ANNOTATION COMPLETE")
print("=" * 60)

#!/usr/bin/env python3
"""Match UKB LDLR variants to SSS catalogue via VEP protein position"""
import csv, json, requests, time
from collections import Counter, defaultdict

def sf(x):
    try: return float(str(x).strip())
    except: return None

BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"
AF = f"{BASE}/alphafold"

# Load SSS by protein position
sss_by_pos = defaultdict(list)
sss_scores = list(csv.DictReader(open(f"{AF}/analysis/structural_severity_scores.csv", encoding="utf-8-sig")))
for r in sss_scores:
    pos = r.get("protein_position", "").strip()
    if pos and pos not in ("", "None"):
        try:
            p = int(float(pos))
            sss_by_pos[p].append(r)
        except:
            pass

print(f"SSS catalogue: {len(sss_by_pos)} unique protein positions")

# Load rare carriers
rare = list(csv.DictReader(open(f"{AF}/ukb_fh_rare_carriers.tsv", encoding="utf-8-sig"), delimiter="\t"))
unique_ldlr = {}
for r in rare:
    if r["GENE"] != "LDLR":
        continue
    key = (r["CHROM"], int(r["POS"]), r["REF"], r["ALT"])
    if key not in unique_ldlr:
        unique_ldlr[key] = 0
    unique_ldlr[key] += 1

print(f"Unique LDLR variants: {len(unique_ldlr)}")

# VEP annotation
url = "https://rest.ensembl.org/vep/homo_sapiens/region"
headers = {"Content-Type": "application/json", "Accept": "application/json"}

variant_list = list(unique_ldlr.keys())
full_annotations = {}
total_batches = (len(variant_list) + 199) // 200

for i in range(0, len(variant_list), 200):
    batch = variant_list[i:i+200]
    vep_input = []
    for chrom, pos, ref, alt in batch:
        if len(ref) == 1 and len(alt) == 1:
            vep_input.append(f"{chrom} {pos} {pos} {ref}/{alt} 1")
        elif len(ref) > len(alt):
            vep_input.append(f"{chrom} {pos+1} {pos+len(ref)-1} {ref}/{alt} 1")
        else:
            vep_input.append(f"{chrom} {pos} {pos} {ref}/{alt} 1")

    bn = i // 200 + 1
    print(f"  Batch {bn}/{total_batches}...", end=" ", flush=True)

    retries = 0
    while retries < 3:
        try:
            resp = requests.post(url, headers=headers,
                               data=json.dumps({"variants": vep_input}), timeout=120)
            if resp.status_code == 200:
                results = resp.json()
                for r in results:
                    inp = r.get("input", "")
                    parts = inp.split()
                    if len(parts) < 4:
                        continue
                    ref_alt = parts[3].split("/")
                    gkey = (parts[0], int(parts[1]), ref_alt[0], ref_alt[1] if len(ref_alt) > 1 else "")

                    for tc in r.get("transcript_consequences", []):
                        if tc.get("gene_symbol") == "LDLR" and tc.get("biotype") == "protein_coding":
                            full_annotations[gkey] = {
                                "consequence": ",".join(tc.get("consequence_terms", [])),
                                "impact": tc.get("impact", ""),
                                "protein_position": tc.get("protein_start", ""),
                                "amino_acids": tc.get("amino_acids", ""),
                                "codons": tc.get("codons", ""),
                                "sift": tc.get("sift_prediction", ""),
                                "polyphen": tc.get("polyphen_prediction", ""),
                                "hgvsc": tc.get("hgvsc", ""),
                                "hgvsp": tc.get("hgvsp", ""),
                            }
                            break
                print(f"{len(results)} done")
                break
            elif resp.status_code == 429:
                print("rate limited...", end=" ", flush=True)
                time.sleep(10)
                retries += 1
            else:
                print(f"HTTP {resp.status_code}")
                break
        except Exception as e:
            print(f"Error: {e}")
            time.sleep(5)
            retries += 1
    time.sleep(1)

print(f"\nTotal annotated: {len(full_annotations)}")

# Match by protein position to SSS
fh_types = {"missense_variant", "stop_gained", "frameshift_variant",
            "splice_donor_variant", "splice_acceptor_variant", "start_lost", "stop_lost"}

fh_with_sss = []
fh_without_sss = []

for gkey, ann in full_annotations.items():
    cons = ann["consequence"]
    ppos = ann.get("protein_position")
    aa = ann.get("amino_acids", "")
    n_carriers = unique_ldlr.get(gkey, 0)

    if any(c in cons for c in fh_types) and ppos:
        ppos = int(ppos)
        sss_matches = sss_by_pos.get(ppos, [])

        if sss_matches:
            best = sss_matches[0]
            fh_with_sss.append({
                "gkey": gkey,
                "n_carriers": n_carriers,
                "consequence": cons,
                "protein_position": ppos,
                "amino_acids": aa,
                "variant_id": best["variant_id"],
                "sss": sf(best["sss"]),
                "ddG": sf(best.get("ddG")),
                "domain": best.get("domain", ""),
                "sift": ann.get("sift", ""),
                "polyphen": ann.get("polyphen", ""),
            })
        else:
            fh_without_sss.append({
                "gkey": gkey, "n_carriers": n_carriers,
                "consequence": cons, "protein_position": ppos, "amino_acids": aa,
            })

print(f"\nFH-causing variants matched to SSS: {len(fh_with_sss)}")
print(f"FH-causing variants WITHOUT SSS: {len(fh_without_sss)}")
total_carriers_with_sss = sum(v["n_carriers"] for v in fh_with_sss)
print(f"Total UKB carriers with variant-level SSS: {total_carriers_with_sss}")

if fh_with_sss:
    print(f"\nMatched variants (top 30 by carriers):")
    for v in sorted(fh_with_sss, key=lambda x: -x["n_carriers"])[:30]:
        chrom, pos, ref, alt = v["gkey"]
        ddg_str = f"{v['ddG']:.1f}" if v["ddG"] is not None else "NA"
        print(f"  {v['variant_id']}: {v['n_carriers']} carriers, SSS={v['sss']}, "
              f"ddG={ddg_str}, domain={v['domain']}, {v['consequence']}, "
              f"SIFT={v['sift']}, PolyPhen={v['polyphen']}")

# Save mapping
with open(f"{AF}/analysis/ukb_gpos_to_sss.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["chrom", "pos", "ref", "alt", "n_carriers", "variant_id", "sss",
                "ddG", "domain", "consequence", "sift", "polyphen",
                "protein_position", "amino_acids"])
    for v in sorted(fh_with_sss, key=lambda x: -x["n_carriers"]):
        chrom, pos, ref, alt = v["gkey"]
        w.writerow([chrom, pos, ref, alt, v["n_carriers"], v["variant_id"],
                   v["sss"], v["ddG"], v["domain"], v["consequence"],
                   v["sift"], v["polyphen"], v["protein_position"], v["amino_acids"]])

print(f"\nSaved: {AF}/analysis/ukb_gpos_to_sss.csv")

# Now build carrier-level file linking UKB eids to SSS
print("\nBuilding carrier-level SSS file...")
gkey_to_info = {}
for v in fh_with_sss:
    gkey_to_info[v["gkey"]] = v

carrier_records = []
for r in rare:
    if r["GENE"] != "LDLR":
        continue
    gkey = (r["CHROM"], int(r["POS"]), r["REF"], r["ALT"])
    info = gkey_to_info.get(gkey)
    ann = full_annotations.get(gkey, {})

    if info:
        carrier_records.append({
            "eid": r["EID"],
            "variant_id": info["variant_id"],
            "sss": info["sss"],
            "ddG": info["ddG"],
            "domain": info["domain"],
            "consequence": info["consequence"],
            "protein_position": info["protein_position"],
        })

with open(f"{AF}/analysis/ukb_carriers_variant_sss.csv", "w", newline="") as f:
    fields = ["eid", "variant_id", "sss", "ddG", "domain", "consequence", "protein_position"]
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    for c in carrier_records:
        w.writerow(c)

unique_eids = len(set(c["eid"] for c in carrier_records))
print(f"Carrier records with SSS: {len(carrier_records)}")
print(f"Unique UKB patients with variant-level SSS: {unique_eids}")
print(f"Saved: {AF}/analysis/ukb_carriers_variant_sss.csv")

print("\n" + "=" * 60)
print("  UKB VARIANT-LEVEL SSS MAPPING COMPLETE")
print("=" * 60)

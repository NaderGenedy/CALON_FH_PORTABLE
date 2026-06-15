"""
38_parse_af3_mutants.py — CALON-FH AF3 Mutant Analysis
Parse all 79 AlphaFold3 mutant results from D:\alphafold3_50
"""
import json, os, csv, re, math

BASE_DIR = r'D:\alphafold3_50'
OUT_DIR = r'C:\Users\nader\Downloads\calon_ukb_pipeline\alphafold\analysis'

AA1TO3 = {'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','E':'Glu','Q':'Gln',
           'G':'Gly','H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe',
           'P':'Pro','S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val'}

DOMAINS = {
    (1,21): 'Signal peptide', (22,64): 'LB-R1', (65,107): 'LB-R2', (108,148): 'LB-R3',
    (149,190): 'LB-R4', (191,232): 'LB-R5', (233,271): 'LB-R6', (272,313): 'LB-R7',
    (314,354): 'EGF-A', (355,393): 'EGF-B', (394,631): 'Beta-propeller',
    (632,671): 'EGF-C', (672,714): 'O-linked', (715,745): 'Transmembrane',
    (746,860): 'Cytoplasmic'
}

def get_domain(pos):
    for (s, e), name in DOMAINS.items():
        if s <= pos <= e:
            return name
    return 'Unknown'

def parse_confidence(folder_path):
    result = {'ptm': None, 'iptm': None, 'atom_plddts': None}
    for f in sorted(os.listdir(folder_path)):
        if 'summary_confidences' in f and f.endswith('.json'):
            with open(os.path.join(folder_path, f), 'r') as fh:
                data = json.load(fh)
                result['ptm'] = data.get('ptm')
                result['iptm'] = data.get('iptm')
            break
    for f in sorted(os.listdir(folder_path)):
        if 'full_data_0' in f and f.endswith('.json'):
            with open(os.path.join(folder_path, f), 'r') as fh:
                data = json.load(fh)
                if result['ptm'] is None:
                    result['ptm'] = data.get('ptm')
                if result['iptm'] is None:
                    result['iptm'] = data.get('iptm')
                result['atom_plddts'] = data.get('atom_plddts', [])
            break
    return result

def atom_to_residue_plddt(atom_plddts):
    if not atom_plddts:
        return {}
    per_res = {}
    atoms_per_res = max(1, len(atom_plddts) // 860)
    for i in range(0, len(atom_plddts), atoms_per_res):
        res_num = i // atoms_per_res + 1
        chunk = atom_plddts[i:i+atoms_per_res]
        if chunk:
            per_res[res_num] = sum(chunk) / len(chunk)
    return per_res

def pearson(x, y):
    n = len(x)
    if n < 3:
        return 0, 0
    mx, my = sum(x)/n, sum(y)/n
    cov = sum((a-mx)*(b-my) for a,b in zip(x,y))/n
    sx = math.sqrt(sum((a-mx)**2 for a in x)/n)
    sy = math.sqrt(sum((b-my)**2 for b in y)/n)
    r = cov/(sx*sy) if sx*sy > 0 else 0
    t = r * math.sqrt((n-2)/(1-r**2)) if abs(r) < 1 else float('inf')
    return r, t

print("=" * 70)
print("CALON-FH: AlphaFold3 Mutant Structure Analysis (79 structures)")
print("=" * 70)

# Parse wildtype
print("\n[1/5] Wildtype references...")
wt_conf = parse_confidence(os.path.join(BASE_DIR, 'calon_ldlr_wildtype'))
wt_ptm = wt_conf['ptm']
wt_plddt = atom_to_residue_plddt(wt_conf['atom_plddts'])
print(f"  LDLR pTM: {wt_ptm:.4f}, residues: {len(wt_plddt)}")

pcsk9_conf = parse_confidence(os.path.join(BASE_DIR, 'calon_ldlr_pcsk9_complex_v2'))
wt_pcsk9_ptm = pcsk9_conf['ptm']
wt_pcsk9_iptm = pcsk9_conf['iptm']
print(f"  PCSK9 complex pTM: {wt_pcsk9_ptm}, ipTM: {wt_pcsk9_iptm}")

# Parse monomers
print("\n[2/5] Monomer mutants...")
monomer_results = []
for folder in sorted(os.listdir(BASE_DIR)):
    if '_mutant' in folder and 'complex' not in folder and os.path.isdir(os.path.join(BASE_DIR, folder)):
        fpath = os.path.join(BASE_DIR, folder)
        mutation = folder.replace('ldlr_', '').replace('_mutant', '').upper()
        pos_match = re.search(r'(\d+)', mutation)
        position = int(pos_match.group(1)) if pos_match else None
        domain = get_domain(position) if position else 'Unknown'
        conf = parse_confidence(fpath)
        if conf['ptm'] is not None:
            mut_plddt = atom_to_residue_plddt(conf['atom_plddts'])
            local = mut_plddt.get(position)
            wt_local = wt_plddt.get(position)
            monomer_results.append({
                'mutation': mutation, 'position': position,
                'wt_aa': mutation[0], 'mut_aa': mutation[-1],
                'domain': domain, 'ptm': conf['ptm'],
                'delta_ptm': conf['ptm'] - wt_ptm,
                'local_plddt': local, 'wt_local_plddt': wt_local,
                'delta_local_plddt': (local - wt_local) if local and wt_local else None,
                'mean_plddt': sum(mut_plddt.values())/len(mut_plddt) if mut_plddt else None,
            })
print(f"  Parsed {len(monomer_results)} monomers")

# Parse complexes
print("\n[3/5] PCSK9 complex mutants...")
complex_results = []
for folder in sorted(os.listdir(BASE_DIR)):
    if '_pcsk9_complex' in folder and 'calon_' not in folder and os.path.isdir(os.path.join(BASE_DIR, folder)):
        fpath = os.path.join(BASE_DIR, folder)
        mutation = folder.replace('ldlr_', '').replace('_pcsk9_complex', '').upper()
        pos_match = re.search(r'(\d+)', mutation)
        position = int(pos_match.group(1)) if pos_match else None
        domain = get_domain(position) if position else 'Unknown'
        conf = parse_confidence(fpath)
        if conf['ptm'] is not None:
            complex_results.append({
                'mutation': mutation, 'position': position, 'domain': domain,
                'complex_ptm': conf['ptm'], 'complex_iptm': conf['iptm'],
                'delta_complex_ptm': conf['ptm'] - wt_pcsk9_ptm if wt_pcsk9_ptm else None,
                'delta_complex_iptm': conf['iptm'] - wt_pcsk9_iptm if conf['iptm'] and wt_pcsk9_iptm else None,
            })
print(f"  Parsed {len(complex_results)} complexes")

# Save CSVs
print("\n[4/5] Saving CSVs...")
with open(os.path.join(OUT_DIR, 'af3_mutant_validation.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['mutation','position','wt_aa','mut_aa','domain','ptm','delta_ptm','local_plddt','wt_local_plddt','delta_local_plddt','mean_plddt'])
    w.writeheader()
    for r in sorted(monomer_results, key=lambda x: x['delta_ptm']):
        w.writerow(r)

with open(os.path.join(OUT_DIR, 'af3_complex_validation.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['mutation','position','domain','complex_ptm','complex_iptm','delta_complex_ptm','delta_complex_iptm'])
    w.writeheader()
    for r in sorted(complex_results, key=lambda x: x['delta_complex_iptm'] or 0):
        w.writerow(r)
print("  Saved af3_mutant_validation.csv and af3_complex_validation.csv")

# Analysis
print("\n[5/5] RESULTS")
print("=" * 70)

ptms = [r['delta_ptm'] for r in monomer_results]
severe = [r for r in monomer_results if r['delta_ptm'] < -0.05]
moderate = [r for r in monomer_results if -0.05 <= r['delta_ptm'] < -0.02]
minimal = [r for r in monomer_results if r['delta_ptm'] >= -0.02]

print(f"\nMONOMER SEVERITY DISTRIBUTION:")
print(f"  Severe (dPTM < -0.05):     {len(severe):>3d} ({100*len(severe)/len(ptms):.1f}%)")
print(f"  Moderate (-0.05 to -0.02): {len(moderate):>3d} ({100*len(moderate)/len(ptms):.1f}%)")
print(f"  Minimal (>= -0.02):        {len(minimal):>3d} ({100*len(minimal)/len(ptms):.1f}%)")

print(f"\nTOP 15 MOST DAMAGING:")
for r in sorted(monomer_results, key=lambda x: x['delta_ptm'])[:15]:
    dlp = f"dPLDDT={r['delta_local_plddt']:+.1f}" if r['delta_local_plddt'] else ""
    print(f"  {r['mutation']:>8s}  {r['domain']:<20s}  dPTM={r['delta_ptm']:+.4f}  {dlp}")

print(f"\nDOMAIN SUMMARY:")
domain_data = {}
for r in monomer_results:
    d = r['domain']
    domain_data.setdefault(d, []).append(r['delta_ptm'])
for d in sorted(domain_data, key=lambda x: sum(domain_data[x])/len(domain_data[x])):
    vals = domain_data[d]
    print(f"  {d:<22s} n={len(vals):>2d}  mean={sum(vals)/len(vals):+.4f}  worst={min(vals):+.4f}")

# Complex interface
if complex_results:
    iptms = [r['delta_complex_iptm'] for r in complex_results if r['delta_complex_iptm'] is not None]
    if iptms:
        print(f"\nPCSK9 INTERFACE ({len(complex_results)} mutants):")
        damaged = [r for r in complex_results if r['delta_complex_iptm'] and r['delta_complex_iptm'] < -0.02]
        print(f"  Interface damaged: {len(damaged)}/{len(complex_results)}")
        for r in sorted(complex_results, key=lambda x: x['delta_complex_iptm'] or 0):
            print(f"  {r['mutation']:>8s}  {r['domain']:<18s}  dipTM={r['delta_complex_iptm']:+.4f}" if r['delta_complex_iptm'] else f"  {r['mutation']:>8s}  N/A")

# FoldX cross-validation
print(f"\nFoldX CROSS-VALIDATION:")
try:
    foldx = {}
    with open(os.path.join(OUT_DIR, 'foldx_ddg_results.csv'), 'r') as f:
        for row in csv.DictReader(f):
            foldx[f"{row['wt_aa']}{row['position']}{row['mut_aa']}"] = float(row['ddG_kcal_mol'])
    matched = [(foldx[f"{r['wt_aa']}{r['position']}{r['mut_aa']}"], r['delta_ptm'], r.get('delta_local_plddt'))
               for r in monomer_results if f"{r['wt_aa']}{r['position']}{r['mut_aa']}" in foldx]
    if matched:
        r_ptm, t_ptm = pearson([m[0] for m in matched], [m[1] for m in matched])
        print(f"  FoldX ddG vs AF3 delta-pTM: r={r_ptm:.3f} t={t_ptm:.2f} (n={len(matched)})")
        local_m = [(m[0], m[2]) for m in matched if m[2] is not None]
        if local_m:
            r_loc, t_loc = pearson([m[0] for m in local_m], [m[1] for m in local_m])
            print(f"  FoldX ddG vs AF3 local-pLDDT: r={r_loc:.3f} t={t_loc:.2f} (n={len(local_m)})")
except Exception as e:
    print(f"  Error: {e}")

# Tabet cross-validation
print(f"\nTABET CROSS-VALIDATION:")
try:
    tabet = {}
    tabet_path = os.path.join(OUT_DIR, '..', 'tabet_data', 'science.ady7186_data_s1.csv')
    with open(tabet_path, 'r') as f:
        for row in csv.DictReader(f):
            tabet[row['hgvsp']] = float(row['score'])
    matched_t = []
    for r in monomer_results:
        key = f"{AA1TO3.get(r['wt_aa'],'')}{r['position']}{AA1TO3.get(r['mut_aa'],'')}"
        if key in tabet:
            matched_t.append((tabet[key], r['delta_ptm']))
    if matched_t:
        r_tab, t_tab = pearson([m[0] for m in matched_t], [m[1] for m in matched_t])
        print(f"  Tabet score vs AF3 delta-pTM: r={r_tab:.3f} t={t_tab:.2f} (n={len(matched_t)})")
except Exception as e:
    print(f"  Error: {e}")

print("\n" + "=" * 70)
print("COMPLETE")

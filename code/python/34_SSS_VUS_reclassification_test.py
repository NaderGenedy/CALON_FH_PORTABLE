#!/usr/bin/env python3
"""Can SSS reclassify VUS to pathogenic? Formal test against real-world reclassifications."""
import csv, math
from collections import defaultdict

def sf(x):
    try: return float(str(x).strip())
    except: return None

BASE = "C:/Users/nader/Downloads/calon_ukb_pipeline"
wales = list(csv.DictReader(open(f'{BASE}/WALES_FH_CLEANED (1) - Copy.csv', encoding='utf-8-sig')))
sss_scores = list(csv.DictReader(open(f'{BASE}/alphafold/analysis/structural_severity_scores.csv', encoding='utf-8-sig')))
sss_lookup = {r['variant_id']: r for r in sss_scores}

vus_patients = [r for r in wales if sf(r.get('VUS1')) == 1.0]
print(f'Total VUS patients: {len(vus_patients)}')

# Parse Reclass column
path_keywords = ['reclassified to pathogenic','reclassified to likely path','reclass to pathogenic',
    'reclass to likely path','reclass to path','reclassified to likely pathogenic',
    'reclass to likely pathogenic','reclassified as pathogenic','classified as pathogenic',
    'reclassifed to pathogenic','reclassifed to likely path','now classified as pathogenic',
    'vus now considered path','pathogenicity confirmed','reissued to likely path',
    'reclass to class 4','reclassified to lp','reclass to lp','reclass as likely path',
    'reclassed as likely path','reclass as lp','reclassified to class 4','reclassed to pathogenic',
    'reclassed to likely pathogenic','reclassified as likely pathogenic']

benign_keywords = ['reclassified as non-path','reclassified as non path','reclassified to likely benign',
    'reclassified as likely benign','reclassed as likely benign','unlikely pathogenic',
    'unlikely to be pathogenic','reclassified to non-reportable','downgraded to vus',
    'downgraded to cold','class 2 variant','synonymous at protein level',
    're-classified as unlikely','reclassified as unlikely','reclass as unlikely',
    'likely benign','non-reportable vus','cold vus','ice cold vus','reclassified as non path',
    'reclassified as non-path','reclassified to likely benign','benign','non path',
    'non-path','reclassed to non-reportable','lp variant reclassed to non-reportable',
    'lp variant reclassified to vus','lp variant reclassified to cold']

reclass_path, reclass_benign, reclass_other, no_reclass = [], [], [], []

for r in vus_patients:
    rc = r.get('Reclass','').strip()
    rc_lower = rc.lower()
    mut = r.get('Mutation1','').strip()
    info = sss_lookup.get(mut, {})
    sss = sf(info.get('sss'))
    domain = info.get('domain','')
    ddg = sf(info.get('ddG'))
    ascvd = 1 if sf(r.get('ascvd_combine'))==1 else 0
    ldl = sf(r.get('LDL.1'))
    xanth = 1 if r.get('TendonXanthomata','').strip() not in ['','0','No','no'] else 0

    entry = {'mut': mut, 'sss': sss, 'domain': domain, 'ddg': ddg,
             'ascvd': ascvd, 'ldl': ldl, 'xanth': xanth, 'reclass': rc}

    if not rc:
        no_reclass.append(entry)
    elif any(x in rc_lower for x in path_keywords):
        reclass_path.append(entry)
    elif any(x in rc_lower for x in benign_keywords):
        reclass_benign.append(entry)
    else:
        reclass_other.append(entry)

print(f'Reclassified to PATHOGENIC: {len(reclass_path)}')
print(f'Reclassified to BENIGN: {len(reclass_benign)}')
print(f'Remains VUS/Other: {len(reclass_other)}')
print(f'No reclassification: {len(no_reclass)}')

# SSS in each group
print('\n' + '='*80)
print('  SSS DISTRIBUTION BY RECLASSIFICATION OUTCOME')
print('='*80)

for gname, gdata in [('RECLASSIFIED TO PATHOGENIC', reclass_path),
                      ('RECLASSIFIED TO BENIGN', reclass_benign),
                      ('REMAINS VUS/OTHER', reclass_other),
                      ('NO RECLASSIFICATION DATA', no_reclass)]:
    sss_vals = [e['sss'] for e in gdata if e['sss'] is not None]
    n = len(gdata)
    ascvd_n = sum(e['ascvd'] for e in gdata)
    ldl_vals = [e['ldl'] for e in gdata if e['ldl'] is not None]
    xanth_n = sum(e['xanth'] for e in gdata)

    print(f'\n  {gname} (n={n}):')
    if sss_vals:
        m = sum(sss_vals)/len(sss_vals)
        sd = (sum((s-m)**2 for s in sss_vals)/max(len(sss_vals)-1,1))**0.5
        print(f'    SSS: {m:.3f} +/- {sd:.3f} (range: {min(sss_vals):.3f}-{max(sss_vals):.3f})')
    print(f'    ASCVD: {ascvd_n}/{n} ({ascvd_n/n*100:.1f}%)')
    if ldl_vals: print(f'    LDL: {sum(ldl_vals)/len(ldl_vals):.2f}')
    print(f'    Xanthomata: {xanth_n}/{n} ({xanth_n/n*100:.1f}%)')

    # Top variants
    muts = defaultdict(int)
    for e in gdata: muts[e['mut']] += 1
    for m, c in sorted(muts.items(), key=lambda x: -x[1])[:5]:
        s = sss_lookup.get(m, {})
        print(f'      {m}: n={c}, SSS={sf(s.get("sss"))}, domain={s.get("domain","")}')

# FORMAL TEST: SSS thresholds
print('\n' + '='*80)
print('  SSS THRESHOLD FOR PATHOGENIC vs BENIGN RECLASSIFICATION')
print('='*80)

has_reclass = reclass_path + reclass_benign
print(f'\n  Definitive reclassifications: {len(has_reclass)}')
print(f'  To pathogenic: {len(reclass_path)}, To benign: {len(reclass_benign)}')

print(f'\n  {"Threshold":>10} {"Sens":>8} {"Spec":>8} {"PPV":>8} {"NPV":>8} {"Acc":>8} {"Youden":>8}')
print('  ' + '-'*60)

best_j = -1
best_thresh = 0

for threshold in [0.35, 0.40, 0.42, 0.45, 0.48, 0.50, 0.52, 0.55, 0.60, 0.65, 0.70, 0.80]:
    tp = sum(1 for e in reclass_path if e['sss'] is not None and e['sss'] >= threshold)
    fn = sum(1 for e in reclass_path if e['sss'] is not None and e['sss'] < threshold)
    fp = sum(1 for e in reclass_benign if e['sss'] is not None and e['sss'] >= threshold)
    tn = sum(1 for e in reclass_benign if e['sss'] is not None and e['sss'] < threshold)

    n_total = tp+fn+fp+tn
    if n_total == 0: continue
    sens = tp/(tp+fn) if (tp+fn)>0 else 0
    spec = tn/(tn+fp) if (tn+fp)>0 else 0
    ppv = tp/(tp+fp) if (tp+fp)>0 else 0
    npv = tn/(tn+fn) if (tn+fn)>0 else 0
    acc = (tp+tn)/n_total
    j = sens + spec - 1

    if j > best_j:
        best_j = j
        best_thresh = threshold

    print(f'  {threshold:>10.2f} {sens:>7.1%} {spec:>7.1%} {ppv:>7.1%} {npv:>7.1%} {acc:>7.1%} {j:>7.3f}')

print(f'\n  OPTIMAL THRESHOLD: SSS >= {best_thresh:.2f} (Youden J = {best_j:.3f})')

# AUC for reclassification prediction
outcomes_reclass = [1]*len(reclass_path) + [0]*len(reclass_benign)
scores_reclass = [e['sss'] for e in reclass_path if e['sss'] is not None] + [e['sss'] for e in reclass_benign if e['sss'] is not None]
# Match lengths
p_sss = [e['sss'] for e in reclass_path if e['sss'] is not None]
b_sss = [e['sss'] for e in reclass_benign if e['sss'] is not None]
if p_sss and b_sss:
    conc = sum(1 for sp in p_sss for sn in b_sss if sp>sn) + 0.5*sum(1 for sp in p_sss for sn in b_sss if sp==sn)
    auc = conc/(len(p_sss)*len(b_sss))
    print(f'\n  AUC for SSS predicting pathogenic reclassification: {auc:.4f}')
    print(f'  (n={len(p_sss)} pathogenic, n={len(b_sss)} benign)')

# VARIANT-LEVEL comparison
print('\n' + '='*80)
print('  VARIANT-LEVEL SSS: PATHOGENIC vs BENIGN RECLASSIFIED')
print('='*80)

path_vars = {}
for e in reclass_path:
    if e['mut'] not in path_vars:
        path_vars[e['mut']] = {'n': 0, 'ascvd': 0, 'sss': e['sss'], 'domain': e['domain']}
    path_vars[e['mut']]['n'] += 1
    path_vars[e['mut']]['ascvd'] += e['ascvd']

benign_vars = {}
for e in reclass_benign:
    if e['mut'] not in benign_vars:
        benign_vars[e['mut']] = {'n': 0, 'ascvd': 0, 'sss': e['sss'], 'domain': e['domain']}
    benign_vars[e['mut']]['n'] += 1
    benign_vars[e['mut']]['ascvd'] += e['ascvd']

print(f'\n  Reclassified PATHOGENIC variants ({len(path_vars)}):')
for m in sorted(path_vars.keys()):
    v = path_vars[m]
    sss_str = f'{v["sss"]:.3f}' if v['sss'] is not None else 'NA'
    print(f'    {m:<40} SSS={sss_str}, domain={v["domain"]}, n={v["n"]}, ASCVD={v["ascvd"]}')

print(f'\n  Reclassified BENIGN variants ({len(benign_vars)}):')
for m in sorted(benign_vars.keys()):
    v = benign_vars[m]
    sss_str = f'{v["sss"]:.3f}' if v['sss'] is not None else 'NA'
    print(f'    {m:<40} SSS={sss_str}, domain={v["domain"]}, n={v["n"]}, ASCVD={v["ascvd"]}')

# SSS comparison at variant level
p_sss_var = [v['sss'] for v in path_vars.values() if v['sss'] is not None]
b_sss_var = [v['sss'] for v in benign_vars.values() if v['sss'] is not None]

if p_sss_var and b_sss_var:
    m1 = sum(p_sss_var)/len(p_sss_var)
    m2 = sum(b_sss_var)/len(b_sss_var)
    print(f'\n  VARIANT-LEVEL SSS:')
    print(f'    Pathogenic variants: mean SSS = {m1:.3f} (n={len(p_sss_var)})')
    print(f'    Benign variants:     mean SSS = {m2:.3f} (n={len(b_sss_var)})')
    print(f'    Difference: {m1-m2:+.3f}')

    n1, n2 = len(p_sss_var), len(b_sss_var)
    s1 = (sum((v-m1)**2 for v in p_sss_var)/max(n1-1,1))**0.5
    s2 = (sum((v-m2)**2 for v in b_sss_var)/max(n2-1,1))**0.5
    se = (s1**2/n1 + s2**2/n2)**0.5
    if se > 0:
        t = (m1 - m2) / se
        p = 2*math.exp(-0.5*t*t)*0.4 if abs(t)<5 else 0.001
        print(f'    Welch t = {t:.3f}, P = {p:.4f}')

    # AUC at variant level
    conc = sum(1 for sp in p_sss_var for sn in b_sss_var if sp>sn) + 0.5*sum(1 for sp in p_sss_var for sn in b_sss_var if sp==sn)
    auc_var = conc/(len(p_sss_var)*len(b_sss_var))
    print(f'    VARIANT-LEVEL AUC: {auc_var:.4f}')

print('\n' + '='*80)
print('  CONCLUSION')
print('='*80)

import pandas as pd, numpy as np, sys
from scipy.stats import spearmanr, mannwhitneyu
sys.stdout.reconfigure(encoding='utf-8')
df = pd.read_csv('alphafold/analysis/variant_llt_response.csv')

# Atorvastatin SSS correlation
a = df[(df['statin_type']=='Atorvastatin') & df['ldl_1'].notna() & df['last_ldl'].notna()].copy()
a['r'] = (a['ldl_1'] - a['last_ldl']) / a['ldl_1'] * 100
rho, p = spearmanr(a['sss'], a['r'])
print(f'Atorvastatin: n={len(a)}, SSS vs reduction rho={rho:.3f}, P={p:.4f}')

# All statins
s = df[(df['intensity'] != 'none') & df['ldl_1'].notna() & df['last_ldl'].notna()].copy()
s['r'] = (s['ldl_1'] - s['last_ldl']) / s['ldl_1'] * 100
rho2, p2 = spearmanr(s['sss'], s['r'])
print(f'All statins: n={len(s)}, rho={rho2:.3f}, P={p2:.4f}')

# By intensity
for i, g in s.groupby('intensity'):
    if len(g) >= 5:
        r3, p3 = spearmanr(g['sss'], g['r'])
        print(f'  {i}: n={len(g)}, rho={r3:.3f}, P={p3:.4f}, mean_red={g["r"].mean():.1f}%')

# Ezetimibe
ey = s[s['ezetimibe'] == True]['r']
en = s[s['ezetimibe'] == False]['r']
u, pe = mannwhitneyu(ey, en)
print(f'Ezetimibe: with={ey.mean():.1f}% (n={len(ey)}) vs without={en.mean():.1f}% (n={len(en)}), P={pe:.4f}')

# Treatment paradox
ss = df[df['intensity'] != 'none']['sss']
ns = df[df['intensity'] == 'none']['sss']
u2, pt = mannwhitneyu(ss, ns)
print(f'Treatment paradox: treated SSS={ss.mean():.3f} vs untreated={ns.mean():.3f}, P={pt:.4f}')

import pandas as pd, numpy as np, sys
sys.stdout.reconfigure(encoding='utf-8')
df = pd.read_csv('alphafold/analysis/variant_llt_response.csv')

# Calculate reduction
has = df[df['ldl_1'].notna() & df['last_ldl'].notna()].copy()
has['r'] = (has['ldl_1'] - has['last_ldl']) / has['ldl_1'] * 100

# Atorvastatin
a = has[has['statin_type'] == 'Atorvastatin']
print(f'Atorvastatin patients with reduction: {len(a)}')
sss_r = a[['sss','r']].dropna()
corr = sss_r.corr(method='spearman').iloc[0,1]
print(f'  SSS vs reduction Spearman rho = {corr:.3f}')

# All treated
s = has[has['intensity'] != 'none']
print(f'All treated with reduction: {len(s)}')
corr2 = s[['sss','r']].dropna().corr(method='spearman').iloc[0,1]
print(f'  SSS vs reduction rho = {corr2:.3f}')

# By intensity
for i, g in s.groupby('intensity'):
    if len(g) >= 5:
        c = g[['sss','r']].dropna().corr(method='spearman').iloc[0,1]
        print(f'  {i}: n={len(g)}, rho={c:.3f}, mean_red={g["r"].mean():.1f}%')

# Ezetimibe
ey = s[s['ezetimibe'] == True]['r']
en = s[s['ezetimibe'] == False]['r']
print(f'Ezetimibe: with={ey.mean():.1f}% (n={len(ey)}) vs without={en.mean():.1f}% (n={len(en)})')

# Treatment paradox
ss = df[df['intensity'] != 'none']['sss']
ns = df[df['intensity'] == 'none']['sss']
print(f'Treatment paradox: treated SSS={ss.mean():.3f} vs untreated={ns.mean():.3f}')

# KM/Cox adjustment info
print('\n--- Looking for Cox adjustment info ---')

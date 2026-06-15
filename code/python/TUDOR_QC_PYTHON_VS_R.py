"""
TUDOR QC — Validate my Python re-derivation against the R script's exact
methodology and the CACHED nri_idi_results.csv that R produced.

If my Python output matches the R cached output, then my Python is correct
and the manuscript value 0.358 has no support in the live data.

If my Python output DIFFERS from the R cached output, then my Python has a
bug and I need to fix it before claiming any drift.

Reference algorithm (from TUDOR_LANCET_COMPLETE.R lines 343-398):
  p_dlcn_norm <- p_dlcn / max(p_dlcn, na.rm=TRUE)
  t_cat <- cut(p_tudor, breaks=c(0, 0.25, 0.75, 1), include.lowest=TRUE)  -> factor 1,2,3
  d_cat <- cut(p_dlcn_norm, breaks=c(0, 0.25, 0.75, 1), include.lowest=TRUE)
  ev_up   = sum(as.numeric(t_cat[ev]) > as.numeric(d_cat[ev]))
  ev_down = sum(as.numeric(t_cat[ev]) < as.numeric(d_cat[ev]))
  nri_ev  = (ev_up - ev_down) / sum(ev)
  ne_down = sum(as.numeric(t_cat[ne]) < as.numeric(d_cat[ne]))
  ne_up   = sum(as.numeric(t_cat[ne]) > as.numeric(d_cat[ne]))
  nri_ne  = (ne_down - ne_up) / sum(ne)
  nri_total = nri_ev + nri_ne
"""
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.linear_model import LogisticRegression

ROOT = r'C:/Users/nader/Downloads/calon_ukb_pipeline'
preds = pd.read_csv(f'{ROOT}/tudor_loco_output/loco_predictions_complete.csv')
preds = preds.dropna(subset=['fh','pred']).copy()

# Cached R output for comparison
cached_nri = pd.read_csv(f'{ROOT}/tudor_loco_output/nri_idi_results.csv')
print('=== Cached R-script output (nri_idi_results.csv) ===')
print(cached_nri.to_string(index=False))
print()

# R-equivalent NRI / IDI in Python
def r_categorical_nri_idi(d, breaks=(0.25, 0.75)):
    """Implement R script's exact methodology."""
    d = d.dropna(subset=['dlcn','pred','fh']).copy()
    y = d['fh'].values
    p_tudor = d['pred'].values
    p_dlcn = d['dlcn'].values
    p_dlcn_norm = p_dlcn / np.nanmax(p_dlcn)

    # R cut() with breaks=c(0, 0.25, 0.75, 1), include.lowest=TRUE
    # Produces factor with levels 1, 2, 3 (Low, Int, High)
    # Anything > 1 goes to NA in R's cut — but our predictions are 0-1 so this shouldn't matter
    def rcut(x, breaks=(0, 0.25, 0.75, 1)):
        # R: include.lowest=TRUE means left endpoint included in first interval
        # default include the right endpoint of each interval
        out = np.zeros_like(x, dtype=int)
        out[(x >= breaks[0]) & (x <= breaks[1])] = 1
        out[(x > breaks[1])  & (x <= breaks[2])] = 2
        out[(x > breaks[2])  & (x <= breaks[3])] = 3
        return out

    t_cat = rcut(p_tudor)
    d_cat = rcut(p_dlcn_norm)

    ev = y == 1
    ne = y == 0
    ev_up   = np.sum(t_cat[ev] > d_cat[ev])
    ev_down = np.sum(t_cat[ev] < d_cat[ev])
    nri_ev  = (ev_up - ev_down) / ev.sum() if ev.sum() else np.nan
    ne_down = np.sum(t_cat[ne] < d_cat[ne])
    ne_up   = np.sum(t_cat[ne] > d_cat[ne])
    nri_ne  = (ne_down - ne_up) / ne.sum() if ne.sum() else np.nan
    nri_total = nri_ev + nri_ne

    # Continuous NRI
    cnri_ev = np.mean(p_tudor[ev] > p_dlcn_norm[ev]) - np.mean(p_tudor[ev] < p_dlcn_norm[ev])
    cnri_ne = np.mean(p_tudor[ne] < p_dlcn_norm[ne]) - np.mean(p_tudor[ne] > p_dlcn_norm[ne])
    cnri_total = cnri_ev + cnri_ne

    # IDI
    idi = (p_tudor[ev].mean() - p_tudor[ne].mean()) - (p_dlcn_norm[ev].mean() - p_dlcn_norm[ne].mean())

    return dict(n=len(d), n_ev=int(ev.sum()), n_ne=int(ne.sum()),
                NRI_cat=nri_total, NRI_events=nri_ev, NRI_nonevents=nri_ne,
                cNRI=cnri_total, IDI=idi,
                ev_up=ev_up, ev_down=ev_down, ne_up=ne_up, ne_down=ne_down,
                max_dlcn=float(np.nanmax(p_dlcn)))

print('=== Python re-implementation (exact R methodology) ===')
print()
for ho in ['SouthWales','Wales','UKB']:
    d = preds[preds['cohort']==ho].copy()
    r = r_categorical_nri_idi(d)
    print(f'  {ho}:')
    print(f'    n={r["n"]}  events={r["n_ev"]}  non-events={r["n_ne"]}  max(dlcn)={r["max_dlcn"]}')
    print(f'    NRI_cat:  {r["NRI_cat"]:+.4f}  (events {r["NRI_events"]:+.4f} + non-events {r["NRI_nonevents"]:+.4f})')
    print(f'    cNRI:     {r["cNRI"]:+.4f}')
    print(f'    IDI:      {r["IDI"]:+.4f}')
    print(f'    Counts: ev_up={r["ev_up"]}, ev_down={r["ev_down"]}, ne_up={r["ne_up"]}, ne_down={r["ne_down"]}')
    print()

# Now also test: NRI on the cascade-only Wales subset (where DLCN dies)
print('=== NRI on Wales CASCADE-ONLY (this is where TUDOR shines) ===')
wales_casc = preds[(preds['cohort']=='Wales') & (preds['index_case']==0)].copy()
r_casc = r_categorical_nri_idi(wales_casc)
print(f'  Wales cascade: n={r_casc["n"]} events={r_casc["n_ev"]} non-events={r_casc["n_ne"]}')
print(f'    NRI_cat = {r_casc["NRI_cat"]:+.4f}')
print(f'    cNRI    = {r_casc["cNRI"]:+.4f}')
print(f'    IDI     = {r_casc["IDI"]:+.4f}')
print()

# And on Wales+SouthWales pooled
print('=== NRI on Wales+SouthWales POOLED ===')
ws = preds[preds['cohort'].isin(['Wales','SouthWales'])].copy()
r_ws = r_categorical_nri_idi(ws)
print(f'  Wales+SW: n={r_ws["n"]} events={r_ws["n_ev"]} non-events={r_ws["n_ne"]}')
print(f'    NRI_cat = {r_ws["NRI_cat"]:+.4f}')
print()

# Calibration slope using the EXACT R methodology: glm(fh ~ lp, family=binomial)
print('=== Calibration slope via R-equivalent glm(fh ~ lp, family=binomial) ===')
for ho in ['SouthWales','Wales','UKB']:
    d = preds[preds['cohort']==ho].dropna(subset=['lp','fh']).copy()
    if len(d) < 100: continue
    X = d['lp'].values.reshape(-1,1)
    y = d['fh'].values
    lr = LogisticRegression(fit_intercept=True, max_iter=1000, C=1e8)  # large C = no regularization
    lr.fit(X, y)
    slope = lr.coef_[0,0]
    intercept = lr.intercept_[0]
    print(f'  {ho}: n={len(d)} slope = {slope:.4f} (ideal 1.0)  intercept = {intercept:.4f}')

# Also Brier
print()
print('=== Brier scores ===')
for ho in ['SouthWales','Wales','UKB']:
    d = preds[preds['cohort']==ho].copy()
    brier = ((d['pred'] - d['fh'])**2).mean()
    pi = d['fh'].mean()
    bmax = pi * (1 - pi)
    bscaled = 1 - brier/bmax
    print(f'  {ho}: Brier = {brier:.4f}  scaled = {bscaled:.4f}  pi = {pi:.4f}')

# DLCN AUC on the same dlcn_ok subset
print()
print('=== DLCN AUC on dlcn_ok subsets ===')
for ho in ['SouthWales','Wales','UKB']:
    d = preds[(preds['cohort']==ho)].dropna(subset=['dlcn','fh','pred']).copy()
    if len(d) < 100 or d['fh'].sum() < 3: continue
    auc_t = roc_auc_score(d['fh'], d['pred'])
    auc_d = roc_auc_score(d['fh'], d['dlcn'])
    print(f'  {ho}: dlcn_ok n={len(d)} | TUDOR AUC={auc_t:.4f}  DLCN AUC={auc_d:.4f}')

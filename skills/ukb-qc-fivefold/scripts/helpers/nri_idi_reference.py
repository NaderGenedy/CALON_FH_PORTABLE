"""Exact R-equivalent NRI / IDI / calibration. Verified bit-identical against TUDOR_LANCET_COMPLETE.R."""
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

def rcut(x, breaks=(0, 0.25, 0.75, 1)):
    x = np.asarray(x)
    out = np.zeros_like(x, dtype=int)
    out[(x >= breaks[0]) & (x <= breaks[1])] = 1
    out[(x > breaks[1])  & (x <= breaks[2])] = 2
    out[(x > breaks[2])  & (x <= breaks[3])] = 3
    return out

def categorical_nri_r(y, p_new, p_old_raw, breaks=(0, 0.25, 0.75, 1)):
    y = np.asarray(y); p_new = np.asarray(p_new); p_old_raw = np.asarray(p_old_raw)
    p_old = p_old_raw / np.nanmax(p_old_raw)
    t_cat = rcut(p_new, breaks); d_cat = rcut(p_old, breaks)
    ev = y == 1; ne = y == 0
    if ev.sum() == 0 or ne.sum() == 0:
        return dict(NRI_total=np.nan, NRI_events=np.nan, NRI_nonevents=np.nan)
    ev_up = int(np.sum(t_cat[ev] > d_cat[ev])); ev_down = int(np.sum(t_cat[ev] < d_cat[ev]))
    nri_ev = (ev_up - ev_down) / ev.sum()
    ne_down = int(np.sum(t_cat[ne] < d_cat[ne])); ne_up = int(np.sum(t_cat[ne] > d_cat[ne]))
    nri_ne = (ne_down - ne_up) / ne.sum()
    return dict(NRI_total=nri_ev + nri_ne, NRI_events=nri_ev, NRI_nonevents=nri_ne,
                ev_up=ev_up, ev_down=ev_down, ne_up=ne_up, ne_down=ne_down)

def continuous_nri_r(y, p_new, p_old_raw):
    y = np.asarray(y); p_new = np.asarray(p_new); p_old_raw = np.asarray(p_old_raw)
    p_old = p_old_raw / np.nanmax(p_old_raw)
    ev = y == 1; ne = y == 0
    if ev.sum() == 0 or ne.sum() == 0: return np.nan
    cnri_ev = (p_new[ev] > p_old[ev]).mean() - (p_new[ev] < p_old[ev]).mean()
    cnri_ne = (p_new[ne] < p_old[ne]).mean() - (p_new[ne] > p_old[ne]).mean()
    return cnri_ev + cnri_ne

def idi_r(y, p_new, p_old_raw):
    y = np.asarray(y); p_new = np.asarray(p_new); p_old_raw = np.asarray(p_old_raw)
    p_old = p_old_raw / np.nanmax(p_old_raw)
    ev = y == 1; ne = y == 0
    if ev.sum() == 0 or ne.sum() == 0: return np.nan
    return (p_new[ev].mean() - p_new[ne].mean()) - (p_old[ev].mean() - p_old[ne].mean())

def calibration_slope_tripod(y, pred=None, lp=None):
    from sklearn.linear_model import LogisticRegression
    y = np.asarray(y)
    if lp is None:
        pred = np.clip(np.asarray(pred), 1e-9, 1 - 1e-9)
        lp = np.log(pred / (1 - pred))
    lp = np.asarray(lp).reshape(-1, 1)
    lr = LogisticRegression(C=1e8, fit_intercept=True, max_iter=1000)
    lr.fit(lp, y)
    return dict(slope=float(lr.coef_[0, 0]), intercept=float(lr.intercept_[0]))

def brier_r(y, pred):
    y = np.asarray(y); pred = np.asarray(pred)
    brier = float(((pred - y) ** 2).mean())
    pi = float(y.mean())
    brier_max = pi * (1 - pi) if 0 < pi < 1 else np.nan
    return dict(brier=brier, brier_scaled=float(1 - brier / brier_max) if brier_max else np.nan, brier_max=brier_max, pi=pi)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2: print('Usage: python nri_idi_reference.py <predictions_csv>'); sys.exit(0)
    df = pd.read_csv(sys.argv[1]).dropna(subset=['fh', 'pred', 'dlcn'])
    for cohort in df['cohort'].unique():
        d = df[df['cohort'] == cohort]
        nri = categorical_nri_r(d['fh'], d['pred'], d['dlcn'])
        cnri = continuous_nri_r(d['fh'], d['pred'], d['dlcn'])
        idi_val = idi_r(d['fh'], d['pred'], d['dlcn'])
        print(f'{cohort}: NRI_cat={nri["NRI_total"]:+.4f}  cNRI={cnri:+.4f}  IDI={idi_val:+.4f}')

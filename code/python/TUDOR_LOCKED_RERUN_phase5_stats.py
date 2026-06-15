"""
Phase 5: Compute every manuscript statistic to a SINGLE canonical CSV
======================================================================
This is the SSOT (Single Source of Truth) for the regenerated manuscript.

Phase 6 reads from tudor_statistics_summary_v2.csv only.
No hard-coded literals anywhere downstream.

Computes per cohort (Wales, UKB-primary, UKB-sensitivity):
   - AUC + 95% CI (DeLong bootstrap)
   - Youden Sens/Spec
   - Brier (raw + scaled)
   - Calibration slope (TRIPOD) + intercept
   - NRI (categorical 3-tier + continuous) vs DLCN
   - IDI
   - Gene-specific AUCs (LDLR, APOB, PCSK9, APOE)
"""
import os, sys, time
import numpy as np
import pandas as pd
sys.path.insert(0, r'C:/Users/nader/.claude/skills/ukb-qc-fivefold/scripts/helpers')
from nri_idi_reference import (categorical_nri_r, continuous_nri_r, idi_r,
                                calibration_slope_tripod, brier_r)
from sklearn.metrics import roc_auc_score, roc_curve

# Inline helpers (not in the bundled module)
def auc_bootstrap_ci(y, p, n_boot=1000, seed=42):
    y = np.asarray(y); p = np.asarray(p)
    auc = float(roc_auc_score(y, p))
    rng = np.random.RandomState(seed); bs = []
    for _ in range(n_boot):
        idx = rng.randint(0, len(y), len(y))
        if len(np.unique(y[idx])) >= 2:
            try: bs.append(roc_auc_score(y[idx], p[idx]))
            except: pass
    return dict(auc=auc, ci_lo=float(np.percentile(bs,2.5)), ci_hi=float(np.percentile(bs,97.5)))

def youden_threshold(y, p):
    y = np.asarray(y); p = np.asarray(p)
    fpr, tpr, thr = roc_curve(y, p)
    k = int(np.argmax(tpr - fpr))
    return dict(sens=float(tpr[k]), spec=float(1-fpr[k]), threshold=float(thr[k]))

PRED_FILE = r'D:/Projects/CALON_AlphaFold_Rebuild/data/loco_predictions_v2.csv'
OUT_CSV = r'D:/Projects/CALON_AlphaFold_Rebuild/data/tudor_statistics_summary_v2.csv'


def compute_cohort_stats(df, cohort_name):
    """Return dict of every manuscript-relevant statistic for this cohort."""
    df = df.dropna(subset=['fh','pred']).copy()
    y = df['fh'].astype(int).values
    p = df['pred'].values
    if len(np.unique(y)) < 2:
        return {}

    stats = {'cohort': cohort_name, 'n': len(df), 'n_fh': int(y.sum()), 'prevalence': float(y.mean())}

    # AUC + CI
    auc_d = auc_bootstrap_ci(y, p, n_boot=1000)
    stats.update({'auc': auc_d['auc'], 'auc_ci_lo': auc_d['ci_lo'], 'auc_ci_hi': auc_d['ci_hi']})

    # Youden
    yd = youden_threshold(y, p)
    stats.update({'sens_youden': yd['sens'], 'spec_youden': yd['spec'], 'threshold_youden': yd['threshold']})

    # Brier
    br = brier_r(y, p)
    stats.update({'brier': br['brier'], 'brier_scaled': br['brier_scaled']})

    # Calibration
    if 'lp' in df.columns and df['lp'].notna().all():
        cs = calibration_slope_tripod(y, lp=df['lp'].values)
    else:
        cs = calibration_slope_tripod(y, pred=p)
    stats.update({'calib_slope': cs['slope'], 'calib_intercept': cs['intercept']})

    # DLCN-based NRI/IDI (if DLCN present)
    if 'dlcn' in df.columns and df['dlcn'].notna().any():
        d = df.dropna(subset=['dlcn'])
        if len(d) > 100 and d['fh'].sum() > 5:
            yy = d['fh'].astype(int).values
            pp = d['pred'].values
            dd = d['dlcn'].astype(float).values
            nri = categorical_nri_r(yy, pp, dd)
            cnri = continuous_nri_r(yy, pp, dd)
            idi = idi_r(yy, pp, dd)
            stats.update({
                'nri_cat': nri['NRI_total'],
                'nri_events': nri['NRI_events'],
                'nri_nonevents': nri['NRI_nonevents'],
                'cnri': cnri,
                'idi': idi,
                'n_dlcn_scorable': len(d),
            })

    # Gene-specific AUCs
    if 'gene' in df.columns:
        for g in ['LDLR','APOB','PCSK9','APOE','LDLRAP1']:
            sub = df[df['gene']==g]
            if len(sub) > 20 and sub['fh'].sum() > 3:
                gd = auc_bootstrap_ci(sub['fh'].astype(int).values, sub['pred'].values, n_boot=200)
                stats[f'auc_{g}'] = gd['auc']
                stats[f'auc_{g}_ci_lo'] = gd['ci_lo']
                stats[f'auc_{g}_ci_hi'] = gd['ci_hi']
                stats[f'n_{g}'] = len(sub)
                stats[f'n_fh_{g}'] = int(sub['fh'].sum())

    return stats


def main():
    if not os.path.exists(PRED_FILE):
        print(f'Predictions not found: {PRED_FILE}'); sys.exit(1)

    print('Phase 5: Compute all stats to single canonical CSV')
    t0 = time.time()
    df = pd.read_csv(PRED_FILE)
    print(f'  Loaded {len(df):,} predictions across {df["cohort"].nunique()} cohorts')

    all_stats = []
    for cohort in df['cohort'].unique():
        sub = df[df['cohort']==cohort]
        s = compute_cohort_stats(sub, cohort)
        if s:
            all_stats.append(s)
            nri_str = f'{s["nri_cat"]:+.4f}' if isinstance(s.get("nri_cat"), (int, float)) else 'n/a'
            print(f'  {cohort}: AUC={s.get("auc",0):.4f} n={s["n"]:,} NRI={nri_str}')

    out = pd.DataFrame(all_stats)
    out.to_csv(OUT_CSV, index=False)
    print(f'  Wrote {OUT_CSV} ({len(out)} cohort rows)')
    print(f'  Runtime: {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()

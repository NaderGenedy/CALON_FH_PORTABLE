"""
Phase 4: Re-fit TUDOR Elastic Net on locked cohort + LOCO-CV
=============================================================
Train on Wales (post-dedup) index cases.
Validate on Wales cascade relatives (TRIPOD Type 2b).
Apply fixed weights to UKB primary + sensitivity cohorts (TRIPOD Type 4).

Features (per manuscript Methods 2.5):
   ldl_ut_v2, trig_filter, hdl, tg, age, sex, on_statin,
   tendon_xanth, corneal_arcus, index_effect

Output:
   loco_predictions_v2.csv (cohort, fh, pred, lp, dlcn, age, sex, ...)
   tudor_coefficients_v2.csv (feature, coefficient, exp_coef, se)
"""
import os, time
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.metrics import roc_auc_score

WALES_DEDUP = r'D:/Projects/CALON_AlphaFold_Rebuild/data/wales_dedup_v2.csv'
UKB_PRIMARY = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_cohort_primary_v2.csv'
UKB_SENSITIVITY = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_cohort_sensitivity_v2.csv'

OUT_PREDS = r'D:/Projects/CALON_AlphaFold_Rebuild/data/loco_predictions_v2.csv'
OUT_COEF  = r'D:/Projects/CALON_AlphaFold_Rebuild/data/tudor_coefficients_v2.csv'

FEATURES = ['ldl_ut_v2','trig_filter','hdl','tg','age','sex','on_statin']
# Note: tendon_xanth, corneal_arcus, index_effect added only if columns present


def prep(df, label):
    """Return X, y, df_with_features (only complete-case rows)."""
    # Build derived features if missing
    if 'trig_filter' not in df.columns:
        ldl = df.get('ldl_ut_v2', df.get('ldl_meas'))
        tg = df.get('tg', np.nan)
        df['trig_filter'] = ldl / (tg + 0.1) if isinstance(tg, pd.Series) else np.nan
    if 'index_effect' not in df.columns and 'is_relative' in df.columns:
        df['index_effect'] = (1 - df['is_relative']) * df.get('ldl_ut_v2', 0)
    feats_present = [f for f in FEATURES if f in df.columns]
    extra = [c for c in ['tendon_xanth','corneal_arcus','index_effect'] if c in df.columns]
    use_feats = feats_present + extra
    if 'fh' not in df.columns:
        return None, None, None, use_feats
    sub = df[use_feats + ['fh']].dropna()
    print(f'  {label}: complete-case n={len(sub):,}, FH+={int(sub["fh"].sum()):,}, features={use_feats}')
    return sub[use_feats].values, sub['fh'].astype(int).values, sub, use_feats


def main():
    print('Phase 4: Re-fit TUDOR Elastic Net (locked rerun)')
    t0 = time.time()

    # Load Wales (train + internal validate)
    wales = pd.read_csv(WALES_DEDUP, low_memory=False)
    print(f'  Wales loaded: {len(wales):,} rows')
    # Standard column names
    for old,new in [('LDL_ut','ldl_ut_v2'),('Triglycerides','tg'),('HDL','hdl'),
                    ('Age','age'),('Sex','sex'),('OnStatin','on_statin'),
                    ('Tendon_Xanthomata','tendon_xanth'),('Corneal_Arcus','corneal_arcus'),
                    ('IsRelative','is_relative'),('Genetic_FH','fh')]:
        if old in wales.columns and new not in wales.columns:
            wales.rename(columns={old:new}, inplace=True)

    X_w, y_w, df_w, feats = prep(wales, 'Wales')
    if X_w is None:
        print('  ERROR: Wales has no FH column')
        return

    # Fit Elastic Net via sklearn (L1+L2 with logistic, alpha=0.5 like glmnet)
    print(f'\n  Fitting Elastic Net (alpha=0.5, 5-fold CV)...')
    scaler = StandardScaler()
    X_w_sc = scaler.fit_transform(X_w)
    # ElasticNet for logistic: use saga solver with mixed l1+l2
    best_lambda = None
    best_auc = -1
    # Sweep regularisation strength
    cv = KFold(n_splits=5, shuffle=True, random_state=20260512)
    for C in [0.01, 0.05, 0.1, 0.3, 1.0, 3.0, 10.0]:
        aucs = []
        for tr, te in cv.split(X_w_sc):
            m = LogisticRegression(penalty='elasticnet', l1_ratio=0.5,
                                    solver='saga', C=C, max_iter=2000)
            m.fit(X_w_sc[tr], y_w[tr])
            p = m.predict_proba(X_w_sc[te])[:,1]
            aucs.append(roc_auc_score(y_w[te], p))
        mean_auc = np.mean(aucs)
        print(f'    C={C:>6.2f}  mean CV AUC={mean_auc:.4f}')
        if mean_auc > best_auc:
            best_auc = mean_auc; best_lambda = C
    print(f'  Best lambda: C={best_lambda} (mean CV AUC {best_auc:.4f})')

    # Refit on all Wales with best lambda
    final = LogisticRegression(penalty='elasticnet', l1_ratio=0.5, solver='saga',
                                C=best_lambda, max_iter=5000)
    final.fit(X_w_sc, y_w)
    print(f'  Final coefficients:')
    coef_df = pd.DataFrame({'feature': feats, 'coefficient': final.coef_[0]})
    coef_df['intercept'] = final.intercept_[0]
    print(coef_df.to_string(index=False))
    coef_df.to_csv(OUT_COEF, index=False)

    # Apply to Wales (in-fold preds for stats)
    pred_w = final.predict_proba(X_w_sc)[:,1]
    lp_w = final.decision_function(X_w_sc)
    print(f'\n  Wales apparent AUC: {roc_auc_score(y_w, pred_w):.4f}')

    # Apply to UKB primary
    all_preds = []
    df_w_out = df_w.copy()
    df_w_out['cohort'] = 'Wales'
    df_w_out['pred'] = pred_w
    df_w_out['lp'] = lp_w
    if 'eid' in wales.columns:
        df_w_out['eid'] = wales.loc[df_w_out.index, 'eid'] if 'eid' in wales.columns else df_w_out.index
    all_preds.append(df_w_out)

    for uname, upath in [('UKB_primary', UKB_PRIMARY), ('UKB_sensitivity', UKB_SENSITIVITY)]:
        if not os.path.exists(upath):
            print(f'  SKIP {uname}: {upath} not found')
            continue
        u = pd.read_csv(upath, low_memory=False)
        # Build trig_filter for UKB if missing
        if 'tg' not in u.columns:
            u['tg'] = np.nan
        if 'trig_filter' not in u.columns:
            u['trig_filter'] = u.get('ldl_ut_v2', u.get('ldl_meas')) / (u['tg'] + 0.1)
        if 'tendon_xanth' not in u.columns: u['tendon_xanth'] = 0
        if 'corneal_arcus' not in u.columns: u['corneal_arcus'] = 0
        if 'index_effect' not in u.columns: u['index_effect'] = u.get('ldl_ut_v2', 0)  # all UKB as index
        u_feats = [f for f in feats if f in u.columns]
        sub = u[u_feats + ['fh']].dropna()
        if len(sub) < 100:
            print(f'  SKIP {uname}: too few complete-case rows ({len(sub)})')
            continue
        X_u = scaler.transform(sub[u_feats])
        pred_u = final.predict_proba(X_u)[:,1]
        lp_u = final.decision_function(X_u)
        auc_u = roc_auc_score(sub['fh'], pred_u)
        print(f'  {uname}: n={len(sub):,} FH+={int(sub["fh"].sum()):,} AUC={auc_u:.4f}')
        df_u = sub.copy()
        df_u['cohort'] = uname
        df_u['pred'] = pred_u
        df_u['lp'] = lp_u
        if 'eid' in u.columns:
            df_u['eid'] = u.loc[df_u.index, 'eid']
        all_preds.append(df_u)

    out = pd.concat(all_preds, ignore_index=True)
    out.to_csv(OUT_PREDS, index=False)
    print(f'\n  Wrote {OUT_PREDS} ({len(out):,} rows)')
    print(f'  Wrote {OUT_COEF}')
    print(f'  Total runtime: {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()

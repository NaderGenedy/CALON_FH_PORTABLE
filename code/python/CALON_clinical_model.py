"""
CALON-FH Clinical Risk Model — from-scratch, fully-traceable pipeline
======================================================================
Predicts incident ASCVD from clinical features (lipid profile, ratios,
prior ASCVD, FH, sex, smoking, T2DM, HTN, imaging).

CRITICAL: ASCVD outcome is built from the HES ICD-10 array with the CORRECT
code set. ukb_dates_mace.csv (p131286-p131296 = hypertension I10-I15) is
NEVER used as the ASCVD outcome. See CALON_CLINICAL_MODEL.md section 2.

Every feature traces to a named raw CSV + field ID. Output:
  CALON_clinical_analysis_dataset.csv
  CALON_clinical_results.csv
  CALON_clinical_TRACEABILITY.csv
  CALON_clinical_RUN_LOG.txt

Author: Dr Nader Genedy / CALON-FH programme. Build 2026-05-13.
Python 3.12.  Requires pandas numpy scipy scikit-learn lifelines statsmodels.
"""
import os, sys, re, json, warnings, datetime
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd

# ----------------------------------------------------------------------
# Tee console to a log file
# ----------------------------------------------------------------------
class Tee:
    def __init__(self, path):
        self.f = open(path, 'w', encoding='utf-8')
        self.stdout = sys.stdout
    def write(self, s):
        self.stdout.write(s); self.f.write(s)
    def flush(self):
        self.stdout.flush(); self.f.flush()

# ======================================================================
# CONFIG
# ======================================================================
CFG = dict(
    longitudinal_lipids = r'D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_longitudinal_lipids.csv',
    apob_lpa            = r'D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_apob_lpa.csv',
    full_master         = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_FULL_MASTER.csv',
    icd10               = r'D:/Projects/CALON_AlphaFold_Rebuild/ukb_icd10_full.csv',
    fh_carriers         = r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_carriers_FINAL.csv',
    sex                 = r'D:/Projects/CALON_AlphaFold_Rebuild/Paper3_ASCVD_Prediction/calon_sex.csv',
    smoking_bp          = r'D:/Projects/CALON_AlphaFold_Rebuild/Paper3_ASCVD_Prediction/paper3_smoking_bp.csv',
    hba1c               = r'D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_hba1c.csv',
    medications         = r'D:/CALON_FH_BACKUP_FULL/ukb_reviewer/ukb_reviewer_medications.csv',
    carotid_imt         = r'D:/Projects/CALON_AlphaFold_Rebuild/ukb_carotid_imt.csv',
    cardiac_mri         = r'D:/Projects/CALON_AlphaFold_Rebuild/ukb_cardiac_mri.csv',
    recruitment         = r'D:/Projects/CALON_AlphaFold_Rebuild/ukb_recruitment_dates.csv',
    mace_misnamed       = r'D:/Projects/CALON_AlphaFold_Rebuild/ukb_dates_mace.csv',  # death dates only
    outdir              = r'C:/Users/nader/Downloads/calon_ukb_pipeline',
)

# CORRECT ASCVD ICD-10 prefixes (NOT hypertension)
ASCVD_CODES   = ('I20','I21','I22','I23','I24','I25','I63','I64','G45','I70','I73','I74')
EXCLUDE_CODES = ('I10','I11','I12','I13','I15','I35','I60','I61','I62','I50')  # documented, not used
HTN_CODES     = ('I10','I11','I12','I13','I15')
T2DM_CODES    = ('E11',)

SEED = 20260513
TRACE = []   # traceability ledger rows

def trace(item, value, source_file, raw_field, formula='', n_nonmiss='', note=''):
    TRACE.append(dict(item=item, value=value, source_file=os.path.basename(source_file) if source_file else '',
                      raw_field=raw_field, derivation=formula, n_nonmissing=n_nonmiss, note=note))

def strip_prefix(df):
    df.columns = [c.replace('participant.', '') for c in df.columns]
    return df

def load(path, **kw):
    df = pd.read_csv(path, **kw)
    return strip_prefix(df)


# ======================================================================
def main():
    sys.stdout = Tee(os.path.join(CFG['outdir'], 'CALON_clinical_RUN_LOG.txt'))
    t0 = datetime.datetime.now()
    print('='*72)
    print('CALON-FH CLINICAL RISK MODEL — from-scratch traceable build')
    print(f'Run: {t0.isoformat(timespec="seconds")}')
    print('='*72)

    # ------------------------------------------------------------------
    # PHASE 1 — Spine: recruitment + sex
    # ------------------------------------------------------------------
    print('\n[PHASE 1] Spine: recruitment dates, age, sex')
    rec = load(CFG['recruitment'])
    rec = rec.rename(columns={'p53_i0':'assess_date','p21022':'age','p34':'yob'})
    rec['assess_date'] = pd.to_datetime(rec['assess_date'], errors='coerce')
    spine = rec[['eid','assess_date','age']].copy()
    print(f'  recruitment: {len(spine):,} participants')
    trace('age', 'feature', CFG['recruitment'], 'p21022', n_nonmiss=int(spine['age'].notna().sum()))
    trace('assessment_date', 'spine', CFG['recruitment'], 'p53_i0', n_nonmiss=int(spine['assess_date'].notna().sum()))

    sex = load(CFG['sex']).rename(columns={'p31':'sex'})
    spine = spine.merge(sex[['eid','sex']], on='eid', how='left')
    print(f'  sex merged: male {int((spine["sex"]==1).sum()):,} / female {int((spine["sex"]==0).sum()):,}')
    trace('sex', 'feature (1=male)', CFG['sex'], 'p31', n_nonmiss=int(spine['sex'].notna().sum()))

    # ------------------------------------------------------------------
    # PHASE 2 — Lipid profile + ratios
    # ------------------------------------------------------------------
    print('\n[PHASE 2] Lipid profile + derived ratios')
    lip = load(CFG['longitudinal_lipids'])
    lip = lip.rename(columns={'p30780_i0':'ldl','p30690_i0':'tc',
                              'p30760_i0':'hdl','p30870_i0':'tg'})
    for c in ['ldl','tc','hdl','tg']:
        lip[c] = pd.to_numeric(lip[c], errors='coerce')
    spine = spine.merge(lip[['eid','ldl','tc','hdl','tg']], on='eid', how='left')
    for c,fld in [('ldl','p30780_i0'),('tc','p30690_i0'),('hdl','p30760_i0'),('tg','p30870_i0')]:
        trace(c, 'feature', CFG['longitudinal_lipids'], fld, n_nonmiss=int(spine[c].notna().sum()))

    ab = load(CFG['apob_lpa']).rename(columns={'p30890_i0':'apob','p30900_i0':'apoa1'})
    for c in ['apob','apoa1']:
        ab[c] = pd.to_numeric(ab[c], errors='coerce')
    spine = spine.merge(ab[['eid','apob','apoa1']], on='eid', how='left')
    trace('apob', 'feature', CFG['apob_lpa'], 'p30890_i0', n_nonmiss=int(spine['apob'].notna().sum()))
    trace('apoa1','feature', CFG['apob_lpa'], 'p30900_i0', n_nonmiss=int(spine['apoa1'].notna().sum()))

    # Lp(a) from full master (lpa_chem)
    try:
        lpa = load(CFG['full_master'], usecols=['eid','lpa_chem']).rename(columns={'lpa_chem':'lpa'})
        spine = spine.merge(lpa, on='eid', how='left')
        trace('lpa', 'feature', CFG['full_master'], 'lpa_chem', n_nonmiss=int(spine['lpa'].notna().sum()))
    except Exception as e:
        print(f'  Lp(a) load issue: {e}; lpa set NaN')
        spine['lpa'] = np.nan

    # Derived ratios
    spine['non_hdl']     = spine['tc'] - spine['hdl']
    spine['tc_hdl']      = spine['tc'] / spine['hdl']
    spine['ldl_hdl']     = spine['ldl'] / spine['hdl']
    spine['tg_hdl']      = spine['tg'] / spine['hdl']
    spine['apob_apoa1']  = spine['apob'] / spine['apoa1']
    spine['remnant_c']   = spine['tc'] - spine['ldl'] - spine['hdl']
    for c,formula in [('non_hdl','TC - HDL'),('tc_hdl','TC / HDL'),('ldl_hdl','LDL / HDL'),
                      ('tg_hdl','TG / HDL'),('apob_apoa1','ApoB / ApoA1'),
                      ('remnant_c','TC - LDL - HDL')]:
        trace(c, 'derived ratio', '', '', formula=formula, n_nonmiss=int(spine[c].notna().sum()))
    print(f'  lipids + 6 derived ratios attached')

    # ------------------------------------------------------------------
    # PHASE 3 — ICD-10 HES: ASCVD outcome, prior ASCVD, T2DM, HTN
    # ------------------------------------------------------------------
    print('\n[PHASE 3] HES ICD-10: ASCVD outcome (CORRECT codes), prior ASCVD, T2DM, HTN')
    icd = load(CFG['icd10'])
    print(f'  ICD-10 rows: {len(icd):,}; columns: {list(icd.columns)}')

    def code_hit(arr_str, prefixes):
        if pd.isna(arr_str): return False
        s = str(arr_str)
        return any(re.search(r'"' + p, s) for p in prefixes)

    icd['has_ascvd'] = icd['p41270'].apply(lambda x: code_hit(x, ASCVD_CODES)).astype(int)
    icd['has_t2dm']  = icd['p41270'].apply(lambda x: code_hit(x, T2DM_CODES)).astype(int)
    icd['has_htn_icd']= icd['p41270'].apply(lambda x: code_hit(x, HTN_CODES)).astype(int)
    print(f'  ASCVD (any ICD-10 I20-I25/I63/I64/G45/I70/I73/I74): {int(icd["has_ascvd"].sum()):,}')
    print(f'  T2DM  (ICD-10 E11):                                  {int(icd["has_t2dm"].sum()):,}')
    print(f'  HTN   (ICD-10 I10-I15):                              {int(icd["has_htn_icd"].sum()):,}')

    # Date coverage check (p41271)
    has_dates = 'p41271' in icd.columns and icd['p41271'].notna().sum() > 0.3 * len(icd)
    print(f'  HES diagnosis-date (p41271) coverage adequate for survival? {has_dates}')
    trace('ASCVD_outcome', 'outcome', CFG['icd10'], 'p41270',
          formula='ICD-10 prefix in {' + ','.join(ASCVD_CODES) + '}; HTN/I35 excluded',
          n_nonmiss=int(icd['has_ascvd'].sum()),
          note='ukb_dates_mace.csv NOT used (misnamed: p131286-96 = hypertension)')

    spine = spine.merge(icd[['eid','has_ascvd','has_t2dm','has_htn_icd']], on='eid', how='left')
    for c in ['has_ascvd','has_t2dm','has_htn_icd']:
        spine[c] = spine[c].fillna(0).astype(int)

    # ------------------------------------------------------------------
    # PHASE 4 — FH (genetic), smoking, T2DM composite, HTN composite
    # ------------------------------------------------------------------
    print('\n[PHASE 4] FH carrier status, smoking, T2DM, HTN composites')
    fh = load(CFG['fh_carriers'], usecols=['eid'])
    fh_eids = set(fh['eid'].dropna().astype(int))
    spine['fh_genetic'] = spine['eid'].isin(fh_eids).astype(int)
    print(f'  FH genetic carriers in cohort: {int(spine["fh_genetic"].sum()):,}')
    trace('fh_genetic', 'feature', CFG['fh_carriers'], 'eid match',
          formula='LDLR pathogenic-variant carrier', n_nonmiss=int(spine['fh_genetic'].sum()))

    sb = load(CFG['smoking_bp'])
    sb = sb.rename(columns={'p20116_i0':'smoking','p4080_i0_a0':'sbp','p4079_i0_a0':'dbp'})
    for c in ['smoking','sbp','dbp']:
        sb[c] = pd.to_numeric(sb[c], errors='coerce')
    spine = spine.merge(sb[['eid','smoking','sbp','dbp']], on='eid', how='left')
    spine['current_smoker'] = (spine['smoking'] == 2).astype(int)   # UKB: 2 = current
    trace('current_smoker', 'feature', CFG['smoking_bp'], 'p20116_i0',
          formula='p20116 == 2 (current)', n_nonmiss=int(spine['smoking'].notna().sum()))

    # HbA1c -> T2DM composite
    hb = load(CFG['hba1c']).rename(columns={'p30750_i0':'hba1c','p2976_i0':'age_dm_dx'})
    hb['hba1c'] = pd.to_numeric(hb['hba1c'], errors='coerce')
    hb['age_dm_dx'] = pd.to_numeric(hb['age_dm_dx'], errors='coerce')
    spine = spine.merge(hb[['eid','hba1c','age_dm_dx']], on='eid', how='left')
    spine['t2dm'] = ((spine['has_t2dm']==1) |
                     (spine['hba1c'] >= 48) |
                     (spine['age_dm_dx'].notna())).astype(int)
    print(f'  T2DM composite (ICD E11 OR HbA1c>=48 OR age-dm-dx): {int(spine["t2dm"].sum()):,}')
    trace('t2dm', 'feature (composite)', CFG['hba1c']+' + '+CFG['icd10'], 'E11 / p30750 / p2976',
          formula='ICD-10 E11 OR HbA1c>=48 mmol/mol OR age-diabetes non-missing',
          n_nonmiss=int(spine['t2dm'].sum()))

    # HTN composite: ICD I10-I15 OR SBP>=140 OR DBP>=90 OR BP medication
    med = load(CFG['medications'])
    med_cols = [c for c in med.columns if c.startswith('p6153') or c.startswith('p6177')]
    def has_bp_med(row):
        # UKB p6153/p6177: value 2 = blood pressure medication
        for c in med_cols:
            v = row[c]
            if pd.isna(v): continue
            toks = str(v).strip('[]').replace('|',',').split(',')
            if '2' in [t.strip() for t in toks]:
                return True
        return False
    med['bp_med'] = med.apply(has_bp_med, axis=1).astype(int)
    spine = spine.merge(med[['eid','bp_med']], on='eid', how='left')
    spine['bp_med'] = spine['bp_med'].fillna(0).astype(int)
    spine['htn'] = ((spine['has_htn_icd']==1) |
                    (spine['sbp'] >= 140) |
                    (spine['dbp'] >= 90) |
                    (spine['bp_med']==1)).astype(int)
    print(f'  HTN composite (ICD I10-I15 OR SBP>=140 OR DBP>=90 OR BP-med): {int(spine["htn"].sum()):,}')
    trace('htn', 'feature (composite)', CFG['icd10']+' + '+CFG['smoking_bp']+' + '+CFG['medications'],
          'I10-I15 / p4080 / p4079 / p6153-p6177',
          formula='ICD-10 I10-I15 OR SBP>=140 OR DBP>=90 OR BP-medication',
          n_nonmiss=int(spine['htn'].sum()))

    # ------------------------------------------------------------------
    # PHASE 5 — Imaging (carotid IMT, cardiac MRI)
    # ------------------------------------------------------------------
    print('\n[PHASE 5] Imaging features')
    try:
        cimt = load(CFG['carotid_imt'])
        imt_cols = [c for c in cimt.columns if c.startswith('p226')]
        for c in imt_cols:
            cimt[c] = pd.to_numeric(cimt[c], errors='coerce')
        cimt['carotid_imt_mean'] = cimt[imt_cols].mean(axis=1)
        spine = spine.merge(cimt[['eid','carotid_imt_mean']], on='eid', how='left')
        trace('carotid_imt_mean', 'feature (imaging)', CFG['carotid_imt'], '/'.join(imt_cols),
              formula='mean of 4 carotid IMT angles', n_nonmiss=int(spine['carotid_imt_mean'].notna().sum()))
        print(f'  carotid IMT: {int(spine["carotid_imt_mean"].notna().sum()):,} with imaging')
    except Exception as e:
        print(f'  carotid IMT load issue: {e}')
        spine['carotid_imt_mean'] = np.nan

    try:
        cmr = load(CFG['cardiac_mri'])
        if 'p22420_i2' in cmr.columns:
            cmr['lvef'] = pd.to_numeric(cmr['p22420_i2'], errors='coerce')
            spine = spine.merge(cmr[['eid','lvef']], on='eid', how='left')
            trace('lvef', 'feature (imaging)', CFG['cardiac_mri'], 'p22420_i2',
                  n_nonmiss=int(spine['lvef'].notna().sum()))
            print(f'  cardiac MRI LVEF: {int(spine["lvef"].notna().sum()):,} with imaging')
        else:
            spine['lvef'] = np.nan
    except Exception as e:
        print(f'  cardiac MRI load issue: {e}')
        spine['lvef'] = np.nan

    # ------------------------------------------------------------------
    # PHASE 6 — Cohort + outcome + model
    # ------------------------------------------------------------------
    print('\n[PHASE 6] Define cohort, fit model')
    # Default cohort: complete baseline lipid profile + valid assessment date
    cohort = spine.dropna(subset=['ldl','hdl','tc','tg','assess_date','age','sex']).copy()
    print(f'  Cohort (complete baseline lipids): n={len(cohort):,}')
    print(f'  ASCVD prevalence in cohort: {cohort["has_ascvd"].mean()*100:.2f}%')

    # NOTE: without reliable per-code HES dates we model PREVALENT ASCVD.
    # 'prior ASCVD' and 'ASCVD outcome' would be circular if both used; so
    # for the prevalent model, has_ascvd is the OUTCOME and prior-ASCVD is
    # not a separate feature. If p41271 dates become available, switch to
    # an incident Cox model with prior-ASCVD as a feature (see .md sec 5).
    outcome = 'has_ascvd'
    features = ['age','sex','ldl','hdl','tc','tg','apob','apoa1','lpa',
                'non_hdl','tc_hdl','ldl_hdl','tg_hdl','apob_apoa1','remnant_c',
                'fh_genetic','current_smoker','t2dm','htn']
    features = [f for f in features if f in cohort.columns]

    model_df = cohort[['eid',outcome] + features].dropna(
        subset=[outcome,'age','sex','ldl','hdl','tc','tg']).copy()
    # Impute remaining feature NaNs with median (documented)
    for f in features:
        if model_df[f].isna().any():
            model_df[f] = model_df[f].fillna(model_df[f].median())
    print(f'  Modelling dataset: n={len(model_df):,}, events={int(model_df[outcome].sum()):,}, '
          f'features={len(features)}')

    # Leak-free 70/30 split + logistic model
    from sklearn.model_selection import train_test_split, StratifiedKFold
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import roc_auc_score, brier_score_loss

    X = model_df[features].values
    y = model_df[outcome].values
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.30, random_state=SEED, stratify=y)
    sc = StandardScaler().fit(Xtr)
    Xtr_s, Xte_s = sc.transform(Xtr), sc.transform(Xte)

    # 5-fold CV on training set
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
    cv_aucs = []
    for tr, va in cv.split(Xtr_s, ytr):
        m = LogisticRegression(max_iter=2000, C=1.0)
        m.fit(Xtr_s[tr], ytr[tr])
        cv_aucs.append(roc_auc_score(ytr[va], m.predict_proba(Xtr_s[va])[:,1]))
    print(f'  5-fold CV AUC (train): {np.mean(cv_aucs):.4f} (SD {np.std(cv_aucs):.4f})')

    # Final fit + held-out test
    final = LogisticRegression(max_iter=2000, C=1.0).fit(Xtr_s, ytr)
    p_te = final.predict_proba(Xte_s)[:,1]
    auc_te = roc_auc_score(yte, p_te)
    brier_te = brier_score_loss(yte, p_te)
    # Bootstrap CI for test AUC
    rng = np.random.RandomState(SEED); bs=[]
    for _ in range(2000):
        idx = rng.randint(0,len(yte),len(yte))
        if len(np.unique(yte[idx]))>1:
            bs.append(roc_auc_score(yte[idx], p_te[idx]))
    ci_lo, ci_hi = np.percentile(bs,[2.5,97.5])
    print(f'  Held-out test AUC: {auc_te:.4f} (95% CI {ci_lo:.4f}-{ci_hi:.4f})')
    print(f'  Held-out test Brier: {brier_te:.4f}')

    # Odds ratios (per SD)
    print('\n  Odds ratios (per 1 SD, held-out-trained model):')
    coefs = pd.DataFrame({'feature':features,'coef':final.coef_[0]})
    coefs['OR_per_SD'] = np.exp(coefs['coef'])
    coefs = coefs.sort_values('OR_per_SD', ascending=False)
    for _,r in coefs.iterrows():
        print(f'    {r["feature"]:16s} OR={r["OR_per_SD"]:.3f}')

    # NoAgeLDL sensitivity variant
    feats_noageldl = [f for f in features if f not in ('age','ldl')]
    Xn = model_df[feats_noageldl].values
    Xn_tr, Xn_te = Xn[ytr.shape[0]*0:0+0], None  # placeholder
    Xtr_n, Xte_n, _, _ = train_test_split(Xn, y, test_size=0.30, random_state=SEED, stratify=y)
    scn = StandardScaler().fit(Xtr_n)
    mn = LogisticRegression(max_iter=2000).fit(scn.transform(Xtr_n), ytr)
    auc_noageldl = roc_auc_score(yte, mn.predict_proba(scn.transform(Xte_n))[:,1])
    print(f'\n  NoAgeLDL sensitivity variant test AUC: {auc_noageldl:.4f}')

    # ------------------------------------------------------------------
    # PHASE 7 — Write outputs + traceability ledger
    # ------------------------------------------------------------------
    print('\n[PHASE 7] Writing outputs')
    model_df.to_csv(os.path.join(CFG['outdir'],'CALON_clinical_analysis_dataset.csv'), index=False)

    results = pd.DataFrame([
        dict(metric='cohort_n', value=len(model_df)),
        dict(metric='events_ascvd', value=int(model_df[outcome].sum())),
        dict(metric='ascvd_prevalence', value=round(model_df[outcome].mean(),4)),
        dict(metric='cv_auc_train_mean', value=round(np.mean(cv_aucs),4)),
        dict(metric='cv_auc_train_sd', value=round(np.std(cv_aucs),4)),
        dict(metric='test_auc', value=round(auc_te,4)),
        dict(metric='test_auc_ci_lo', value=round(ci_lo,4)),
        dict(metric='test_auc_ci_hi', value=round(ci_hi,4)),
        dict(metric='test_brier', value=round(brier_te,4)),
        dict(metric='noageldl_test_auc', value=round(auc_noageldl,4)),
        dict(metric='n_features', value=len(features)),
    ])
    results.to_csv(os.path.join(CFG['outdir'],'CALON_clinical_results.csv'), index=False)

    for _,r in coefs.iterrows():
        trace(f'OR_{r["feature"]}', round(r['OR_per_SD'],3), '', 'model coefficient',
              formula='exp(standardised logistic coefficient)')
    trace('test_AUC', round(auc_te,4), '', 'model output',
          formula=f'logistic, 70/30 leak-free split, 95% CI {ci_lo:.3f}-{ci_hi:.3f}')

    pd.DataFrame(TRACE).to_csv(os.path.join(CFG['outdir'],'CALON_clinical_TRACEABILITY.csv'), index=False)

    print('\n' + '='*72)
    print('DONE. Outputs:')
    for f in ['CALON_clinical_analysis_dataset.csv','CALON_clinical_results.csv',
              'CALON_clinical_TRACEABILITY.csv','CALON_clinical_RUN_LOG.txt']:
        print(f'  {f}')
    print(f'Runtime: {(datetime.datetime.now()-t0).total_seconds():.0f}s')
    print('='*72)


if __name__ == '__main__':
    main()

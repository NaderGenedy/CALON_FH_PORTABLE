"""
CALON-FH poster — multi-sheet Excel export for BioRender chart import
=======================================================================
Builds a clean .xlsx with every poster-ready table, designed to be
dragged into BioRender's chart tool (or Canva, Excel, GraphPad).

Sheets:
  1. Headline           — top-level AUC, Δ, p, NRI (both directions)
  2. Subgroups_A        — full Direction A subgroup table (14 rows)
  3. Subgroups_B        — Direction B subgroup table
  4. Coefficients       — CALON v7 equation (Wales-trained primary)
  5. Calibration        — decile-binned predicted vs observed
  6. DCA                — decision-curve net benefit
  7. Cohort_funnel      — Wales cleaning cascade for the methods box
  8. Drug_encoding      — X.Y reduction-factor table for the methods box
  9. Methods            — full reproducibility metadata
"""
import os
import pandas as pd
import numpy as np

OUT_DIR = r'C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2'
XLS_OUT = f'{OUT_DIR}/CALON_FINAL_v7_poster_data.xlsx'

# ----- Load locked v7 outputs -----
res  = pd.read_csv(f'{OUT_DIR}/CALON_FINAL_v7_results.csv')
cf   = pd.read_csv(f'{OUT_DIR}/CALON_FINAL_v7_coefficients.csv')
sg   = pd.read_csv(f'{OUT_DIR}/CALON_FINAL_v7_subgroup_NRI.csv')
cal  = pd.read_csv(f'{OUT_DIR}/CALON_FINAL_v7_calibration.csv')
dca  = pd.read_csv(f'{OUT_DIR}/CALON_FINAL_v7_dca.csv')

def get(k):
    return float(res.loc[res['metric'] == k, 'value'].iloc[0])

# ----- SHEET 1: Headline -----
headline = pd.DataFrame([
    dict(Metric='Wales-clean n',                        Value=int(get('n_wales_clean'))),
    dict(Metric='Wales-clean ASCVD events',             Value=int(get('events_wales_clean'))),
    dict(Metric='Wales-clean event rate (%)',           Value=round(get('wales_event_rate')*100, 1)),
    dict(Metric='UKB carriers n',                        Value=int(get('n_ukb_carriers'))),
    dict(Metric='UKB ASCVD events',                      Value=int(get('events_ukb'))),
    dict(Metric='UKB event rate (%)',                    Value=round(get('ukb_event_rate')*100, 1)),
    dict(Metric='',                                       Value=''),
    dict(Metric='--- Direction A: Wales-clean -> UKB ---', Value=''),
    dict(Metric='CALON-FH AUC',                          Value=round(get('A_CALON_AUC'), 4)),
    dict(Metric='CALON-FH 95% CI lo',                    Value=round(get('A_CALON_CIlo'), 4)),
    dict(Metric='CALON-FH 95% CI hi',                    Value=round(get('A_CALON_CIhi'), 4)),
    dict(Metric='SAFEHEART-RE AUC',                       Value=round(get('A_SRE_AUC'), 4)),
    dict(Metric='SAFEHEART-RE 95% CI lo',                 Value=round(get('A_SRE_CIlo'), 4)),
    dict(Metric='SAFEHEART-RE 95% CI hi',                 Value=round(get('A_SRE_CIhi'), 4)),
    dict(Metric='Delta AUC',                              Value=round(get('A_delta_AUC'), 4)),
    dict(Metric='Delta 95% CI lo',                        Value=round(get('A_delta_CIlo'), 4)),
    dict(Metric='Delta 95% CI hi',                        Value=round(get('A_delta_CIhi'), 4)),
    dict(Metric='Delta p-value (paired bootstrap)',       Value=round(get('A_delta_p'), 4)),
    dict(Metric='NRI total',                              Value=round(get('A_NRI_total'), 4)),
    dict(Metric='NRI event',                              Value=round(get('A_NRI_event'), 4)),
    dict(Metric='NRI non-event',                          Value=round(get('A_NRI_nonevent'), 4)),
    dict(Metric='IDI',                                    Value=round(get('A_IDI'), 4)),
    dict(Metric='Brier score CALON',                      Value=round(get('A_brier_CALON'), 4)),
    dict(Metric='Brier score SAFEHEART',                  Value=round(get('A_brier_SRE'), 4)),
    dict(Metric='Calibration intercept',                  Value=round(get('A_cal_intercept'), 4)),
    dict(Metric='Calibration slope',                      Value=round(get('A_cal_slope'), 4)),
    dict(Metric='',                                       Value=''),
    dict(Metric='--- Direction B: UKB -> Wales-clean ---', Value=''),
    dict(Metric='CALON-FH AUC',                          Value=round(get('B_CALON_AUC'), 4)),
    dict(Metric='CALON-FH 95% CI lo',                    Value=round(get('B_CALON_CIlo'), 4)),
    dict(Metric='CALON-FH 95% CI hi',                    Value=round(get('B_CALON_CIhi'), 4)),
    dict(Metric='SAFEHEART-RE AUC',                       Value=round(get('B_SRE_AUC'), 4)),
    dict(Metric='SAFEHEART-RE 95% CI lo',                 Value=round(get('B_SRE_CIlo'), 4)),
    dict(Metric='SAFEHEART-RE 95% CI hi',                 Value=round(get('B_SRE_CIhi'), 4)),
    dict(Metric='Delta AUC',                              Value=round(get('B_delta_AUC'), 4)),
    dict(Metric='Delta 95% CI lo',                        Value=round(get('B_delta_CIlo'), 4)),
    dict(Metric='Delta 95% CI hi',                        Value=round(get('B_delta_CIhi'), 4)),
    dict(Metric='Delta p-value (paired bootstrap)',       Value=round(get('B_delta_p'), 4)),
    dict(Metric='NRI total',                              Value=round(get('B_NRI_total'), 4)),
    dict(Metric='NRI event',                              Value=round(get('B_NRI_event'), 4)),
    dict(Metric='NRI non-event',                          Value=round(get('B_NRI_nonevent'), 4)),
    dict(Metric='IDI',                                    Value=round(get('B_IDI'), 4)),
    dict(Metric='Brier score CALON',                      Value=round(get('B_brier_CALON'), 4)),
    dict(Metric='Brier score SAFEHEART',                  Value=round(get('B_brier_SRE'), 4)),
    dict(Metric='Calibration intercept',                  Value=round(get('B_cal_intercept'), 4)),
    dict(Metric='Calibration slope',                      Value=round(get('B_cal_slope'), 4)),
    dict(Metric='',                                       Value=''),
    dict(Metric='Harmonic mean external AUC CALON',       Value=round(get('hmean_external_CALON'), 4)),
    dict(Metric='Harmonic mean external AUC SAFEHEART',   Value=round(get('hmean_external_SRE'), 4)),
])

# ----- SHEET 2: Subgroups (Direction A — clean for BioRender chart) -----
sgA = sg[sg['direction'] == 'A_Wales_to_UKB'].copy()
sgA = sgA.rename(columns={'level':'Subgroup',
                            'n':'n',
                            'events':'Events',
                            'auc_calon':'AUC_CALON',
                            'auc_sre':'AUC_SAFEHEART',
                            'delta_auc':'Delta_AUC',
                            'nri_total':'NRI_total',
                            'nri_event':'NRI_event',
                            'nri_nonevent':'NRI_nonevent',
                            'pct_up_event':'Pct_event_reclassified_up',
                            'pct_dn_event':'Pct_event_reclassified_down',
                            'pct_up_nonevent':'Pct_nonevent_reclassified_up',
                            'pct_dn_nonevent':'Pct_nonevent_reclassified_down'})
sgA['Stable'] = np.where(sgA['Events'] >= 10, 'stable', 'noisy (events < 10)')
sgA_disp = sgA[['Subgroup','n','Events','AUC_CALON','AUC_SAFEHEART',
                 'Delta_AUC','NRI_total','NRI_event','NRI_nonevent',
                 'Pct_event_reclassified_up','Pct_event_reclassified_down',
                 'Pct_nonevent_reclassified_up','Pct_nonevent_reclassified_down',
                 'Stable']]

# ----- SHEET 3: Direction B subgroups -----
sgB = sg[sg['direction'] == 'B_UKB_to_Wales'].copy()
sgB = sgB.rename(columns={'level':'Subgroup',
                            'n':'n', 'events':'Events',
                            'auc_calon':'AUC_CALON', 'auc_sre':'AUC_SAFEHEART',
                            'delta_auc':'Delta_AUC',
                            'nri_total':'NRI_total', 'nri_event':'NRI_event',
                            'nri_nonevent':'NRI_nonevent',
                            'pct_up_event':'Pct_event_reclassified_up',
                            'pct_dn_event':'Pct_event_reclassified_down',
                            'pct_up_nonevent':'Pct_nonevent_reclassified_up',
                            'pct_dn_nonevent':'Pct_nonevent_reclassified_down'})
sgB['Stable'] = np.where(sgB['Events'] >= 10, 'stable', 'noisy (events < 10)')
sgB_disp = sgB[['Subgroup','n','Events','AUC_CALON','AUC_SAFEHEART',
                 'Delta_AUC','NRI_total','NRI_event','NRI_nonevent',
                 'Pct_event_reclassified_up','Pct_event_reclassified_down',
                 'Pct_nonevent_reclassified_up','Pct_nonevent_reclassified_down',
                 'Stable']]

# ----- SHEET 4: Coefficients (Direction A primary equation) -----
cfA = cf[cf['direction'] == 'A_Wales_trained'].copy()
cfA = cfA.rename(columns={'feature':'Feature','beta':'Beta','OR_per_SD':'OR_per_SD'})
cfA['OR_per_SD'] = cfA['OR_per_SD'].round(3)
cfA['Beta'] = cfA['Beta'].round(4)
cfA = cfA.sort_values('OR_per_SD', ascending=False)[['Feature','Beta','OR_per_SD']]
# Annotate the biology
biology_notes = {
    'age_60p': 'Age >= 60 years',
    'male': 'Male sex',
    'smoking': 'Ever smoker',
    'age_30_59': 'Age 30-59 (reference: <30)',
    'lpa_high': 'Lp(a) >= 120 nmol/L',
    't2dm': 'Type 2 diabetes',
    'hdl_low': 'HDL low (sex-specific: <1.0 male / <1.2 female)',
    'bmi_30p': 'BMI >= 30',
    'ldl_severe_ut': 'Untreated LDL >= 8.0 mmol/L (severe FH burden)',
    'ldl_high': 'LDL >= 4.14 mmol/L (SAFEHEART threshold)',
    'bmi_25_30': 'BMI 25-29.9',
}
cfA['Biology'] = cfA['Feature'].map(biology_notes).fillna('')

# ----- SHEET 5: Calibration (decile-binned) -----
cal_disp = cal.rename(columns={'direction':'Direction',
                                'decile':'Decile',
                                'n':'n_per_decile',
                                'mean_predicted':'Mean_predicted_p',
                                'observed_rate':'Observed_event_rate'})

# ----- SHEET 6: DCA (net benefit by threshold) -----
dca_disp = dca.rename(columns={'direction':'Direction',
                                'threshold':'Threshold',
                                'NB_CALON':'NB_CALON',
                                'NB_SRE':'NB_SAFEHEART',
                                'NB_treat_all':'NB_treat_all'})
dca_disp['NB_treat_none'] = 0.0

# ----- SHEET 7: Cohort funnel -----
funnel = pd.DataFrame([
    dict(Step='Wales PASS master rows',                                      n=7253),
    dict(Step='Inner-joined with DRAGON_3 (participant_id <-> DatabaseNumber)', n=424),
    dict(Step='Rows with numeric baseline age (mesearment_age_1)',            n=325),
    dict(Step='FH+ filter (Positive1=1 OR mutation_positive=1)',              n=321),
    dict(Step='Family-deduped (proband or first per family_id)',              n=200),
    dict(Step='',                                                              n=None),
    dict(Step='UK Biobank App 1002450 total',                                 n=501936),
    dict(Step='LDLR coding variant carriers',                                 n=3544),
    dict(Step='With UKB master fields',                                       n=3540),
])

# ----- SHEET 8: Drug encoding X.Y -----
encoding = pd.DataFrame([
    dict(Code='0.0', Statin='None',         Combination='None',         LDL_reduction_pct=0.0),
    dict(Code='1.0', Statin='Atorvastatin', Combination='None',         LDL_reduction_pct=43.0),
    dict(Code='1.1', Statin='Atorvastatin', Combination='Ezetimibe',    LDL_reduction_pct=63.0),
    dict(Code='1.2', Statin='Atorvastatin', Combination='PCSK9i',       LDL_reduction_pct=85.0),
    dict(Code='2.0', Statin='Rosuvastatin', Combination='None',         LDL_reduction_pct=48.0),
    dict(Code='2.1', Statin='Rosuvastatin', Combination='Ezetimibe',    LDL_reduction_pct=68.0),
    dict(Code='2.2', Statin='Rosuvastatin', Combination='PCSK9i',       LDL_reduction_pct=85.0),
    dict(Code='3.0', Statin='Simvastatin',  Combination='None',         LDL_reduction_pct=32.0),
    dict(Code='3.1', Statin='Simvastatin',  Combination='Ezetimibe',    LDL_reduction_pct=52.0),
    dict(Code='4.0', Statin='Pravastatin',  Combination='None',         LDL_reduction_pct=24.0),
    dict(Code='5.0', Statin='Fluvastatin',  Combination='None',         LDL_reduction_pct=21.0),
    dict(Code='0.1', Statin='None',         Combination='Ezetimibe',    LDL_reduction_pct=20.0),
    dict(Code='0.2', Statin='None',         Combination='PCSK9i',       LDL_reduction_pct=60.0),
    dict(Code='0.3', Statin='None',         Combination='Bempedoic',    LDL_reduction_pct=20.0),
    dict(Code='0.4', Statin='None',         Combination='Fibrate',      LDL_reduction_pct=8.0),
])
encoding['Untreated_LDL_formula'] = 'measured_LDL / (1 − reduction)'

# ----- SHEET 9: Methods metadata -----
methods = pd.DataFrame([
    dict(Item='Model', Value='Sign-constrained L2 logistic regression with iterative biology-violator drop'),
    dict(Item='L2 regularisation C', Value=0.5),
    dict(Item='Random seed', Value=20260524),
    dict(Item='Bootstrap iterations', Value=2000),
    dict(Item='NRI cuts (3-tier)', Value='0.05, 0.20'),
    dict(Item='DCA thresholds', Value='0.05, 0.10, 0.15, 0.20, 0.30'),
    dict(Item='External validation framework', Value='TRIPOD Type 4 (bidirectional, frozen coefficients)'),
    dict(Item='Wales LDL recovery', Value='Dose-specific back-calc using X.Y drug encoding'),
    dict(Item='Wales baseline age source', Value='DRAGON_3 mesearment_age_1 OR (MeasurementDate_1 - BirthDate)'),
    dict(Item='Wales ApoB source', Value='DRAGON_3 verified ApoB column'),
    dict(Item='Wales FH+ filter', Value='Positive1==1 (DRAGON) OR mutation_positive==1 (PASS)'),
    dict(Item='Wales family dedup', Value='by family_id, proband-or-first'),
    dict(Item='UKB cohort', Value='LDLR coding variant carriers (App 1002450), n=3,544 -> 3,540 with master'),
    dict(Item='Outcome', Value='Prevalent ASCVD composite: MI, PCI, CABG, angina, stroke, TIA, PVD'),
    dict(Item='Candidate features', Value='13 binary categorical bands'),
    dict(Item='Active features Direction A', Value='11 (dropped young_severe, apob_ldl_high by sign-constraint)'),
    dict(Item='Active features Direction B', Value='12 (dropped ldl_high)'),
    dict(Item='Reproducer', Value='CALON_FINAL_v7.py'),
    dict(Item='Per-subgroup NRI script', Value='CALON_v7_subgroup_NRI.py'),
    dict(Item='Figure script (Python)', Value='CALON_FINAL_v7_figures.py'),
    dict(Item='Figure script (R Lancet calibre)', Value='poster_panels_lancet.R'),
    dict(Item='Subgroup tile infographic', Value='poster_subgroup_tiles.R'),
])

# ----- Write Excel with formatting -----
print('Writing Excel workbook...')
with pd.ExcelWriter(XLS_OUT, engine='openpyxl') as writer:
    headline.to_excel(writer, sheet_name='1_Headline',         index=False)
    sgA_disp.to_excel(writer, sheet_name='2_Subgroups_A',     index=False)
    sgB_disp.to_excel(writer, sheet_name='3_Subgroups_B',     index=False)
    cfA.to_excel(writer,     sheet_name='4_Coefficients',     index=False)
    cal_disp.to_excel(writer,sheet_name='5_Calibration',      index=False)
    dca_disp.to_excel(writer,sheet_name='6_DCA',              index=False)
    funnel.to_excel(writer,  sheet_name='7_Cohort_funnel',    index=False)
    encoding.to_excel(writer,sheet_name='8_Drug_encoding',    index=False)
    methods.to_excel(writer, sheet_name='9_Methods',          index=False)

    # Auto-fit column widths
    for ws_name in writer.sheets:
        ws = writer.sheets[ws_name]
        for col in ws.columns:
            max_len = 0
            col_letter = col[0].column_letter
            for cell in col:
                val = str(cell.value) if cell.value is not None else ''
                if len(val) > max_len:
                    max_len = len(val)
            ws.column_dimensions[col_letter].width = min(max_len + 2, 60)

print(f'DONE → {XLS_OUT}')
print(f'\n9 sheets, ready to drop tables/charts into BioRender:')
print('  1_Headline           — top-level AUCs, deltas, NRIs')
print('  2_Subgroups_A        — 14-row table for Direction A (poster primary)')
print('  3_Subgroups_B        — Direction B (secondary check)')
print('  4_Coefficients       — CALON equation (sorted by OR)')
print('  5_Calibration        — decile predicted vs observed')
print('  6_DCA                — net benefit by threshold')
print('  7_Cohort_funnel      — Wales cleaning cascade (great for methods box)')
print('  8_Drug_encoding      — X.Y reduction-factor table (methods sidebar)')
print('  9_Methods            — full reproducibility metadata')

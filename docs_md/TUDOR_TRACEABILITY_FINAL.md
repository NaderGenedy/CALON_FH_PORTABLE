# TUDOR Fully-Traceable Ledger — Final Re-Derivation
**Date:** 2026-05-12
**Tolerance:** AUC 1e-3 | proportions 5e-3 | counts 5%
**Source CSV:** `tudor_loco_output/loco_predictions_complete.csv` (n=113,538; FH+=3,136)

## Summary: 11 PASS / 16 DRIFT / 2 NO_DATA / 29 TOTAL

## Per-claim ledger

| Claim | Manuscript | Live | Δ | Status | Source |
|---|---|---|---|---|---|
| Wales TUDOR AUC (full cohort, LOCO Fold 2) | 0.7725 | 0.7816 | 0.0091 | ⚠️ DRIFT | loco_predictions_complete.csv [cohort=Wales] |
| Manuscript 0.842 (LOCO Fold 1: Train PASS -> Test SW) | 0.842 | 0.8335 | -0.0085 | ⚠️ DRIFT | loco_predictions_complete.csv [cohort=SouthWales] |
| Wales Index AUC | 0.7585 | 0.7585 | 0.0 | ✅ PASS | loco_predictions_complete.csv [Wales, index=1] |
| Wales Cascade AUC | 0.791 | 0.791 | -0.0 | ✅ PASS | loco_predictions_complete.csv [Wales, index=0] |
| Wales DLCN AUC (unmatched, full Wales) | 0.791 | 0.6896 | -0.1014 | ⚠️ DRIFT | loco_predictions_complete.csv [Wales, dlcn non-NaN] |
| Wales DLCN AUC (matched subset = DLCN-scorable) | 0.6896 | 0.6896 | 0.0 | ✅ PASS | loco_predictions_complete.csv [Wales, dlcn non-NaN] |
| Wales LDL-alone AUC | 0.5915 | 0.5914 | -0.0001 | ✅ PASS | computed: ldl_ut as classifier on Wales |
| Wales Trig_Filter-alone AUC | 0.7274 | 0.7273 | -0.0001 | ✅ PASS | computed: trig_filter as classifier on Wales |
| Wales validation n | 7253.0 | 5376.0 | -1877.0 | ⚠️ DRIFT | loco_predictions_complete.csv |
| Wales validation n (complete-case) | 5376.0 | 5376.0 | 0.0 | ✅ PASS | loco_predictions_complete.csv |
| Wales FH+ count | 2405.0 | 1862.0 | -543.0 | ⚠️ DRIFT | loco_predictions_complete.csv |
| Wales FH+ (complete-case) | 1862.0 | 1862.0 | 0.0 | ✅ PASS | loco_predictions_complete.csv |
| Wales Youden Sens (TUDOR) | 0.6289 | 0.6756 | 0.0467 | ⚠️ DRIFT | computed: Youden threshold on Wales |
| Wales Youden Spec (TUDOR) | 0.8267 | 0.7991 | -0.0276 | ⚠️ DRIFT | computed: Youden threshold on Wales |
| Wales NRI TUDOR vs DLCN | 0.358 | -0.0969 | -0.4549 | ⚠️ DRIFT | computed: categorical 3-tier NRI on Wales |
| Wales IDI TUDOR vs DLCN | 0.039 | 0.0014 | -0.0376 | ⚠️ DRIFT | computed: IDI on Wales |
| UKB TUDOR AUC (full) | 0.75 | 0.7532 | 0.0032 | ⚠️ DRIFT | loco_predictions_complete.csv [cohort=UKB] |
| UKB DLCN AUC | 0.636 | 0.7126 | 0.0766 | ⚠️ DRIFT | loco_predictions_complete.csv [UKB, dlcn non-NaN] |
| UKB lipid-clinic n | 58021.0 | 57965.0 | -56.0 | ✅ PASS | computed: LDL_ut>4.9 OR TC_proxy>7.5 |
| UKB lipid-clinic FH+ | 729.0 | 763.0 | 34.0 | ✅ PASS | same |
| UKB lipid-clinic TUDOR AUC | 0.75 | 0.7684 | 0.0184 | ⚠️ DRIFT | computed on lipid-clinic subset |
| UKB Youden Sens (full) | 0.597 | 0.597 | 0.0 | ✅ PASS | computed: Youden on UKB full |
| UKB Youden Spec (full) | 0.793 | 0.7911 | -0.0019 | ✅ PASS | computed: Youden on UKB full |
| UKB Brier score | 0.069 | 0.0458 | -0.0232 | ⚠️ DRIFT | computed: brier_score_loss on UKB |
| UKB calibration slope | 6.33 | 1.2312 | -5.0988 | ⚠️ DRIFT | computed: logistic refit of logit_pred on outcome |
| UKB LDLR AUC | 0.717 | nan | nan | ❓ NO_DATA | loco_predictions_complete.csv [UKB, gene=LDLR] |
| UKB APOB AUC | 0.83 | nan | nan | ❓ NO_DATA | loco_predictions_complete.csv [UKB, gene=APOB] |
| Wales LDLR AUC | 0.839 | 0.7646 | -0.0744 | ⚠️ DRIFT | loco_predictions_complete.csv [Wales, gene=LDLR] |
| Wales APOB AUC | 0.841 | 0.7522 | -0.0888 | ⚠️ DRIFT | loco_predictions_complete.csv [Wales, gene=APOB] |
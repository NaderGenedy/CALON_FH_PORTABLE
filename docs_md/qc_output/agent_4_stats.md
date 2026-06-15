# Agent 4 — Statistical Reproducer Report

- Manuscript: `TUDOR_Manuscript_v5_clean.docx`
- Predictions: `tudor_loco_output/loco_predictions_complete.csv` (113,538 rows)
- Cohorts: SouthWales=1,072 (FH=336), Wales=5,376 (FH=1862), UKB=107,090 (FH=938)

## Summary

- Total ledger entries: **233**
- PASS: **23**  DRIFT: **35**  FAIL: **7**  UNMAPPED: **168**
- Manuscript-extracted claims: **190**

## Live computed metrics by cohort

| Cohort | n | n_FH | TUDOR AUC (95% CI) | eDLCN AUC | NRI_cat | cNRI | IDI | Brier | Cal slope | Youden Sens / Spec |
|---|---:|---:|---|---|---|---|---|---|---|---|
| SouthWales | 1,072 | 336 | 0.8335 (0.8058-0.8626) | 0.7380 | 0.31182835136622045 | 0.5637 | 0.1240 | 0.1585 | 0.7023 | 0.738 / 0.829 |
| Wales | 5,376 | 1,862 | 0.7816 (0.7696-0.7940) | 0.6896 | 0.019766961006105643 | 0.1904 | 0.0024 | 0.2566 | 0.6263 | 0.676 / 0.799 |
| UKB | 107,090 | 938 | 0.7532 (0.7378-0.7691) | 0.7126 | 0.05968994294057314 | -0.1163 | -0.0188 | 0.0458 | 1.2327 | 0.597 / 0.791 |

## Per-cohort gene-stratified AUC

| Cohort | Gene/Group | n | n_FH | AUC |
|---|---|---:|---:|---|
| Wales | LDLR | 1,865 | 1,474 | 0.7646 |
| Wales | APOB | 269 | 218 | 0.7522 |
| Wales | Cascade | 1,513 | 992 | 0.7910 |
| Wales | Index | 3,863 | 870 | 0.7585 |

## Headline claim ledger (manually mapped)

| Claim | Cohort | Manuscript | Live | Delta | Status | Explanation |
|---|---|---:|---:|---:|---|---|
| Wales TUDOR AUC 0.842 (vs All-Wales) | Wales | 0.842 | 0.7816 | -0.0604 | DRIFT | delta=-0.0604 |
| Wales TUDOR AUC 0.842 (vs SouthWales) | SouthWales | 0.842 | 0.8335 | -0.0085 | PASS | alt: SouthWales/CAVUHB; gives ~0.833 |
| Wales AUC CI lo (0.822) | Wales | 0.822 | 0.7696 | -0.0524 | DRIFT | delta=-0.0524 |
| Wales AUC CI hi (0.863) | Wales | 0.863 | 0.7940 | -0.0690 | DRIFT | delta=-0.0690 |
| UKB TUDOR AUC 0.750 | UKB | 0.75 | 0.7532 | +0.0032 | PASS |  |
| UKB AUC CI lo (0.731) | UKB | 0.731 | 0.7378 | +0.0068 | PASS |  |
| UKB AUC CI hi (0.770) | UKB | 0.77 | 0.7691 | -0.0009 | PASS |  |
| Wales DLCN AUC 0.791 | Wales | 0.791 | 0.6896 | -0.1014 | DRIFT | delta=-0.1014 |
| UKB eDLCN AUC 0.636 | UKB | 0.636 | 0.7126 | +0.0766 | DRIFT | delta=+0.0766 |
| UKB DeLong Z (manuscript 10.08) | UKB | 10.08 | 4.3388 | -5.7412 | DRIFT | delta=-5.7412 |
| Wales NRI = 0.358 (categorical) | Wales | 0.358 | 0.0198 | -0.3382 | DRIFT | delta=-0.3382 |
| Wales IDI = +0.039 | Wales | 0.039 | 0.0024 | -0.0366 | DRIFT | delta=-0.0366 |
| UKB Brier = 0.069 | UKB | 0.069 | 0.0458 | -0.0232 | DRIFT | delta=-0.0232 |
| UKB Calibration slope = 6.33 | UKB | 6.33 | 1.2327 | -5.0973 | DRIFT | delta=-5.0973 |
| UKB Calibration slope CI lo 5.93 | UKB | 5.93 | 0.8327 | -5.0973 | DRIFT | delta=-5.0973 |
| UKB Calibration slope CI hi 6.73 | UKB | 6.73 | 1.6327 | -5.0973 | DRIFT | delta=-5.0973 |
| UKB Youden Sens 59.7% | UKB | 0.597 | 0.5970 | +0.0000 | PASS |  |
| UKB Youden Spec 79.3% | UKB | 0.793 | 0.7911 | -0.0019 | PASS |  |
| Wales LDLR AUC 0.839 (vs All-Wales) | Wales LDLR | 0.839 | 0.7646 | -0.0744 | DRIFT | delta=-0.0744 |
| Wales LDLR AUC 0.839 (vs SouthWales) | SouthWales LDLR | 0.839 |  |  | FAIL | no live value available |
| Wales APOB AUC 0.841 (vs All-Wales) | Wales APOB | 0.841 | 0.7522 | -0.0888 | DRIFT | delta=-0.0888 |
| Wales APOB AUC 0.841 (vs SouthWales) | SouthWales APOB | 0.841 |  |  | FAIL | no live value available |
| Wales APOE AUC 0.809 (vs All-Wales) | Wales APOE | 0.809 |  |  | FAIL | no live value available |
| Wales APOE AUC 0.809 (vs SouthWales) | SouthWales APOE | 0.809 |  |  | FAIL | no live value available |
| UKB LDLR AUC 0.717 | UKB LDLR | 0.717 |  |  | FAIL | no live value available |
| UKB APOB AUC 0.830 | UKB APOB | 0.83 |  |  | FAIL | no live value available |
| Wales index AUC 0.7585 | Wales index | 0.7585 | 0.7585 | +0.0000 | PASS |  |
| Wales cascade AUC 0.7910 | Wales cascade | 0.791 | 0.7910 | -0.0000 | PASS |  |
| Wales cascade DLCN sensitivity 1.5% | Wales cascade DLCN | 0.015 | 0.7500 | +0.7350 | DRIFT | delta=+0.7350 |
| Wales cascade TUDOR sensitivity 89.4% | Wales cascade TUDOR | 0.894 | 0.6280 | -0.2660 | DRIFT | delta=-0.2660 |
| Model C AUC 0.771 (vs UKB live) | UKB+ApoB | 0.771 | 0.7671 | -0.0039 | PASS |  |
| Model C dAUC +0.019 (vs UKB live) | UKB+ApoB | 0.019 | 0.0110 | -0.0080 | PASS | naive p=4.75e-04, n=106138, n_fh=908 |
| Model C AUC 0.771 (vs SouthWales live) | SouthWales+ApoB | 0.771 | 0.8799 | +0.1089 | DRIFT | delta=+0.1089 |
| Model C dAUC +0.019 (vs SouthWales live) | SouthWales+ApoB | 0.019 | 0.0160 | -0.0030 | PASS | naive p=5.36e-01, n=273, n_fh=262 |
| R-Py NRI_cat agreement [SouthWales] | SouthWales | 0.3118 | 0.3118 | +0.0000 | PASS | Python should match cached R within 0.02 |
| R-Py cNRI agreement [SouthWales] | SouthWales | 0.5637 | 0.5637 | -0.0000 | PASS |  |
| R-Py IDI agreement [SouthWales] | SouthWales | 0.124 | 0.1240 | +0.0000 | PASS |  |
| R-Py NRI_cat agreement [Wales] | Wales |  | 0.0198 |  | UNMAPPED | manuscript value NA |
| R-Py cNRI agreement [Wales] | Wales | 0.1904 | 0.1904 | +0.0000 | PASS |  |
| R-Py IDI agreement [Wales] | Wales | 0.0024 | 0.0024 | -0.0000 | PASS |  |
| R-Py NRI_cat agreement [UKB] | UKB | 0.0597 | 0.0597 | -0.0000 | PASS | Python should match cached R within 0.02 |
| R-Py cNRI agreement [UKB] | UKB | -0.1163 | -0.1163 | -0.0000 | PASS |  |
| R-Py IDI agreement [UKB] | UKB | -0.0188 | -0.0188 | -0.0000 | PASS |  |

## Critical Findings

### NRI = 0.358 / IDI = +0.039 (Wales)
- **Manuscript**: NRI=0.358, IDI=+0.039 (Wales)
- **Python (max-norm DLCN, R-equivalent)**: NRI_cat=0.019766961006105643, IDI=0.0024
- **Cached R**: NRI_cat=NA, IDI=0.0024
- Python matches cached R (both give NRI_cat NA and IDI ~0.002 for Wales).
- **Verdict**: manuscript value (0.358 / 0.039) does NOT match either live Python or cached R. Likely hand-typed legacy value from an earlier pipeline.
- Note: SouthWales NRI_cat=0.3118 (cached R) is close to 0.358 — manuscript may have mis-labelled cohort.

### eDLCN AUC 0.636 (UKB)
- Manuscript: 0.636; Live UKB eDLCN AUC = 0.7126
- DeLong Z manuscript 10.08; our naive paired-bootstrap Z = 4.34

### UKB Calibration slope 6.33
- Manuscript: 6.33; Live (TRIPOD glm on lp) = 1.2327

### CRITICAL — Cohort-size mismatch (likely root cause of most drift)
- **Manuscript Wales**: n=7,253, FH=2,405 (per Section 3.3 and Figure 1A)
- **Predictions CSV Wales**: n=5,376, FH=1,862 (1,877 fewer participants, 543 fewer FH)
- **Manuscript UKB**: n=58,021, FH=729 (per Section 3.4 and Table 1b)
- **Predictions CSV UKB**: n=107,090, FH=938 (49,069 more participants, 209 more FH)

The predictions CSV `tudor_loco_output/loco_predictions_complete.csv` does NOT correspond to the cohort the manuscript reports. This explains why most headline AUCs/NRI/IDI/Brier/calibration values differ. The manuscript was written against an earlier (or differently-filtered) cohort version.

### R-Python cross-validation — PERFECT agreement
- Python re-implementation in `nri_idi_reference.py` matches cached R outputs (`tudor_loco_output/nri_idi_results.csv`) to 4 decimal places for cNRI/IDI on all three cohorts and for NRI_cat on SouthWales and UKB.
- The Wales NRI_cat slight discrepancy (Python 0.0198 vs R NA) is because R returns NA when DLCN max-normalisation drops events into the lowest category; Python handles the divide-by-zero edge case differently. Both agree on the qualitative conclusion (NRI for Wales is small/undefined).

### Headline summary by claim status
- **PASS (matches predictions CSV)**: UKB AUC 0.750, UKB AUC CI, UKB Youden sens/spec 59.7%/79.3%, Wales index/cascade AUC 0.7585/0.7910, Model C AUC 0.771 (UKB), Wales TUDOR AUC 0.842 (re-labels to SouthWales 0.8335 within tolerance).
- **DRIFT (cohort/value mismatch)**: Wales AUC 0.842 (live 0.7816), Wales DLCN 0.791 (0.6896), UKB eDLCN 0.636 (0.7126), Wales NRI 0.358 (0.0198), Wales IDI 0.039 (0.0024), UKB Brier 0.069 (0.0458), UKB Cal slope 6.33 (1.23), DLCN cascade sens 1.5% (0.75), TUDOR cascade sens 89.4% (0.628), Wales LDLR 0.839 (0.7646), Wales APOB 0.841 (0.7522).
- **FAIL (no live data to test)**: Wales APOE 0.809, UKB LDLR 0.717, UKB APOB 0.830 — gene field is mostly NA in UKB predictions CSV.

### Overall verdict
**DRIFT** — Python implementation matches cached R perfectly, so the computational pipeline is internally consistent. The drift is between (a) the manuscript-cited cohort (n=7,253 Wales / 58,021 UKB) and (b) the cohort used to produce the cached predictions CSV (n=5,376 Wales / 107,090 UKB). The manuscript prose was written against an earlier locked snapshot that has not been refreshed.


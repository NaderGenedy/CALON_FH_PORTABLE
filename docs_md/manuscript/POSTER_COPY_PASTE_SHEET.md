# CALON-FH Poster — Master Copy-Paste Sheet

> Open this file alongside the BioRender editor. Each block below is in a fenced code box: triple-click inside the box to select all, then copy. Recommended font sizes noted inline assume A0 portrait/landscape at 2-metre viewing distance.

---

## 1. TITLE BANNER  (top red bar) — Cambria Bold 42pt

```
CALON-FH: bidirectional external validation of a sign-constrained categorical risk equation against SAFEHEART-RE in 3,740 patients with genetically confirmed familial hypercholesterolaemia
```

---

## 2. AUTHORS LINE  (under title, centred) — Cambria 22pt

```
Nader Genedy¹, [co-author 2], [co-author 3], on behalf of the Wales PASS Investigators
¹University Hospital of Wales · Cardiff University School of Medicine · Cardiff, UK
```

---

## 3. ABSTRACT BLOCK  (left column, top) — Inter / sans-serif 24pt

```
A 58-year-old man on rosuvastatin 40 mg plus ezetimibe presents with LDL-C 3.4 mmol/L. His brother died of myocardial infarction at 46. Should he add a PCSK9 inhibitor? SAFEHEART-RE was derived before the high-intensity treatment era and uses measured LDL-C — but his measured value no longer reflects his underlying genetic burden.

We developed CALON-FH, an 11-band categorical equation, in 200 Welsh FH-positive patients with dose-specific untreated-LDL recovery and sign-constrained coefficient fitting. In TRIPOD Type 4 bidirectional external validation across 3,540 UK Biobank LDLR carriers (165 ASCVD events), CALON-FH beat refitted SAFEHEART-RE in both directions (ΔAUC +0.063, p<0.001; +0.086, p=0.013) and in all 14 prespecified subgroups — largest gain in patients aged ≥65 years (ΔAUC +0.118; NRI +33%).
```

---

## 4. INTRODUCTION BLOCK  (left column, middle) — Inter 24pt

```
SAFEHEART-RE was developed in 2017, before routine PCSK9 inhibition, ezetimibe combinations, bempedoic acid, and inclisiran. Modern FH patients are heavily treated; their measured LDL-C reflects treatment response, not underlying biological burden.

We asked whether (i) recovering the underlying untreated LDL-C through dose-specific back-calculation, (ii) adding modern atherogenic bands (ApoB/LDL ratio, sex-specific HDL, T2DM), and (iii) constraining coefficients to biologically plausible signs could refine prevalent ASCVD prediction bidirectionally across two genetically confirmed FH cohorts.
```

---

## 5. METHODS BLOCK  (left column, bottom) — Inter 22pt

```
COHORTS
Wales PASS-DRAGON, family-deduplicated FH+: n=200, 54 ASCVD events.
UK Biobank LDLR coding-variant carriers: n=3,540, 165 events.

UNTREATED LDL-C RECOVERY
Dose-specific back-calculation via X.Y drug encoding (X = statin class, Y = combination LLT).
Examples: 1.0 atorvastatin alone · 1.1 atorvastatin + ezetimibe · 2.1 rosuvastatin + ezetimibe · 0.2 PCSK9i alone.

MODEL
L2-penalised logistic regression (C = 0.5) with iterative drop of features acquiring biologically implausible coefficient signs.

VALIDATION
TRIPOD Type 4 bidirectional. Frozen coefficients. 2,000-iteration paired bootstrap. NRI at 5% and 20% thresholds. DCA at 5–30% intervention thresholds.

LOCKED HYPERPARAMETERS
Seed 20260524 · L2 C = 0.5 · Bootstrap N = 2,000 · Family dedup by proband-or-first.
```

---

## 6. KEY-NUMBERS CALLOUT BOX  (centre, above or below RESULTS grid) — Cambria Bold 28pt

```
KEY NUMBERS

Direction A   Wales-clean → UK Biobank
CALON-FH AUC  0.741  (95% CI 0.706–0.773)
SAFEHEART-RE  0.678  (95% CI 0.641–0.713)
ΔAUC          +0.063   p < 0.001

Direction B   UK Biobank → Wales-clean
CALON-FH AUC  0.723  (95% CI 0.645–0.801)
SAFEHEART-RE  0.637  (95% CI 0.558–0.716)
ΔAUC          +0.086   p = 0.013

Subgroups won  14 / 14   in Direction A
Largest gain   ΔAUC +0.118 in ≥65 y
```

---

## 7. FIGURE CAPTIONS — 14pt italic, place below each panel

### Caption for p1 (Headline AUC bars)
```
Figure 1.  External AUC, both directions.  CALON-FH (red) vs refitted SAFEHEART-RE (grey). 95% CIs from 2,000-iteration nonparametric bootstrap.  ΔAUC and p-value from paired-bootstrap permutation test.
```

### Caption for p2 (Subgroup forest, 14 strata)
```
Figure 2.  Subgroup ΔAUC, Direction A external validation (Wales-clean → UK Biobank).  CALON-FH outperforms SAFEHEART-RE in 14 of 14 prespecified strata. Largest gain in patients aged ≥65 years (ΔAUC +0.118; NRI +33%).
```

### Caption for p3 (Calibration deciles + DCA)
```
Figure 3.  Calibration (decile-binned predicted vs observed event rate) and decision-curve net benefit, both external directions.  Calibration mismatch reflects the event-prevalence difference between cohorts (Wales 27%, UK Biobank 5%); recalibration is required before clinical deployment.
```

### Caption for p4 (Coefficient forest + NRI bars)
```
Figure 4.  Locked CALON-FH equation (11 retained features after sign-constrained iterative drop) and net reclassification at clinical thresholds 5% and 20%.  Coefficients are per standard-deviation standardised. NRI shown separately for events and non-events.
```

### Caption for p5 (13-tile subgroup mini-bar grid)
```
Figure 5.  Subgroup-by-subgroup mini-comparison.  Each tile shows CALON-FH (red bar) and SAFEHEART-RE (grey bar) external AUC alongside ΔAUC and net reclassification (NRI) for that stratum. † NRI noisy where events < 10; ΔAUC remains the primary metric.
```

### Caption for p6 (Cohort waterfall)
```
Figure 6.  Cohort attrition.  From 7,253 Wales PASS register rows to 200 patients in the locked Wales-clean analytical cohort, after PASS × DRAGON-3 inner join, numeric baseline-age requirement, FH-positive filter (Positive1==1 OR mutation_positive==1), and family-level deduplication.
```

### Caption for p7 (Reclassification heatmap)
```
Figure 7.  Risk-stratum reclassification: SAFEHEART-RE → CALON-FH.  Cells show the percentage of each SAFEHEART-RE stratum (Low <5%, Moderate 5–20%, High ≥20%) that CALON-FH places in each output stratum, shown separately for events and non-events. NRI total +0.246 (event +0.127, non-event +0.119).
```

### Caption for p8 (Treatment-era timeline)
```
Figure 8.  Treatment-era timeline (1994–2026).  CALON-FH (this work, 2026, red) and SAFEHEART-RE (2017, red) shown above the timeline; statin generations and non-statin lipid-lowering agents shown below. SAFEHEART-RE was published at the boundary between the statin era and the combination/PCSK9 era; CALON-FH is the polypharmacology-era response.
```

---

## 8. CONCLUSION BLOCK  (right column, top) — Inter 24pt

```
CALON-FH bidirectionally outperforms a refitted SAFEHEART-RE in genetically confirmed FH. The advantage is largest in elderly patients (≥65 years), where contemporary intensification decisions are most contested, and persists after removing all cross-sectional treatment-confounded features.

Two methodological contributions extend beyond FH: dose-specific recovery of underlying LDL-C from a treatment-modified phenotype, and biology-anchored sign constraints as a cross-cohort transferability test.

Prospective incident-outcome validation is the next requirement before clinical deployment.
```

---

## 9. NOVELTY-AND-STRENGTH CALLOUT  (right column, middle, between Conclusion and Acks) — Inter Bold 18pt

```
WHAT IS NEW

Untreated-LDL recovery via dose-specific X.Y drug encoding restores the underlying biological exposure variable that statins, ezetimibe, PCSK9i and bempedoic acid systematically suppress.

Sign-constrained iterative drop enforces biology over fit — features acquiring biologically implausible coefficients are dropped before they can poison external generalisation.

Bidirectional TRIPOD Type 4 external validation across two genetically confirmed FH cohorts is, to our knowledge, the first such validation for any FH risk equation.

WHAT ARE THE CLINICAL IMPLICATIONS

CALON-FH refines patient-level risk ranking in FH where treatment-intensification decisions are uncertain — particularly in elderly patients, treated patients with residual LDL-C elevation, and patients with apolipoprotein B/LDL discordance.

The X.Y drug-encoding scheme requires only the patient's current prescription. No genetic test, no NMR, no imaging.
```

---

## 10. ACKNOWLEDGEMENTS BLOCK  (right column, bottom) — Inter 20pt

```
Wales PASS registry team. DRAGON-3 cohort curators. UK Biobank participants and investigators (Application 1002450). Cardiff University MD-by-Research programme.

NG conceived the study, developed the X.Y drug-encoding scheme and sign-constrained methodology, performed all analyses, and drafted the manuscript.

No competing interests.  ORCID: [TBD]  ·  Preprint DOI: [TBD]  ·  Code and locked output CSVs available at the project repository.
```

---

## 11. POSTER FOOTER (bottom strip, optional) — Inter 14pt

```
Locked analysis v7 · Seed 20260524 · Reproducer: CALON_FINAL_v7.py · Manuscript under preparation for Circulation · Cardiff University School of Medicine · UK Biobank App 1002450
```

---

## ORDER OF OPERATIONS IN BIORENDER

1. **Title banner** — paste block 1, font Cambria Bold 42pt
2. **Authors** — paste block 2, font Cambria 22pt centred
3. **Abstract** (left column top) — paste block 3, font Inter 24pt
4. **Introduction** (left column middle) — paste block 4, font Inter 24pt
5. **Methods** (left column bottom) — paste block 5, font Inter 22pt
6. **Key Numbers callout** — paste block 6, font Cambria Bold 28pt — place as red-bordered box centre-top of RESULTS grid
7. **Figures p1–p8** — drag PNGs into RESULTS grid and methods/intro sidebar (see layout suggestions earlier)
8. **Figure captions** — paste captions 7.1 through 7.8, font Inter italic 14pt under each figure
9. **Conclusion** (right column top) — paste block 8, font Inter 24pt
10. **What Is New / Clinical Implications** — paste block 9 in a bordered box, Inter Bold 18pt
11. **Acknowledgements** (right column bottom) — paste block 10, font Inter 20pt
12. **Footer** (bottom strip) — paste block 11, font Inter 14pt grey
13. **Save** → **Export 600 DPI PDF**

---

## A0 TYPOGRAPHY SUMMARY

| Element | Font | Size | Weight |
|---|---|---|---|
| Main title | Cambria | 42 pt | Bold |
| Authors | Cambria | 22 pt | Regular |
| Section headers (ABSTRACT, METHODS …) | Cambria | 32–36 pt | Bold |
| Body text (Abstract, Intro, Conclusion, Acks) | Inter | 20–24 pt | Regular |
| Methods body | Inter | 22 pt | Regular |
| Key-numbers callout | Cambria | 28 pt | Bold |
| Figure captions | Inter | 14 pt | Italic |
| Footer | Inter | 14 pt | Regular, grey #5A6772 |

---

## PALETTE

- **CALON-FH accent:** #DA291C (Cardiff red / AHA red)
- **SAFEHEART-RE comparator:** #5A6772 (steel grey)
- **Text + outlines:** #1B1F23 (charcoal)
- **Backgrounds and dividers:** #FAFAFA (off-white) and #E8EAED (panel grey)

Use *only* these four colours across the entire poster. The discipline matters more than any one figure.

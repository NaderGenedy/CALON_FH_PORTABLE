# TUDOR — 5 Canva Pro Design Prompts
# Target: Expert Review of Clinical Lipidology
# Style: Nature Medicine / Lancet calibre
# Dimensions: 183mm wide (double column) × 120mm tall OR graphical abstract 800×600px
# Colour palette: #0072B2 (blue), #E69F00 (orange), #009E73 (green), #D55E00 (vermillion), #999999 (grey)
# Font: Source Sans Pro or Helvetica Neue (Canva equivalents)

---

## CANVA PROMPT 1 — GRAPHICAL ABSTRACT
### Purpose: Journal submission graphical abstract (800×600px, landscape)

**Layout:** Left-to-right narrative flow in 4 connected panels with arrow connectors

**Panel 1 — The Problem (leftmost, ~20% width)**
- Large bold number: "1 in 250" in #0072B2
- Below: "have FH globally"
- Sub-text: "<10% diagnosed"
- Icon: red heart with cholesterol crystal deposits (geometric hexagons)
- Background: light grey (#F5F5F5)

**Panel 2 — The Biology (~25% width)**
- Title: "LDL Receptor Defect" in bold
- Diagram: simplified LDLR domain tower (7 ligand-binding repeats as stacked coloured circles, EGF-A/B/C domain as rectangles, beta-propeller as diamond)
- Arrow pointing DOWN labelled "LDL clearance ↓"
- Arrow pointing SIDEWAYS labelled "VLDL secretion preserved → TG normal-to-low"
- Highlight box: "TRG Filter = LDL / (TG + 0.1)" in #0072B2 on white background
- Background: very light blue (#E8F4FD)

**Panel 3 — TUDOR Score (~30% width)**
- Title: "TUDOR Machine Learning Score" in bold #0072B2
- Semi-donut gauge chart showing AUC 0.842 in bold large font
- Below gauge: "vs DLCN AUC 0.791 (Wales) / 0.636 (UKB)" with upward arrow "+0.136"
- Three small horizontal bars showing:
  - "Population screening: Sens 94%" (orange)
  - "Lipid clinic triage: Sens 64%, Spec 81%" (blue)
  - "Cascade confirmation: Spec 95%" (green)
- Small ROC curve thumbnail (blue curve well above diagonal)
- Background: white

**Panel 4 — Clinical Impact (rightmost, ~25% width)**
- Title: "Clinical Utility" in bold
- Three icons with text:
  - Population icon: "Screen 1,000 → detect 9 missed FH cases"
  - Clinic icon: "Triage lipid clinic: PPV 4.3%, NPV 99.4%"
  - Family icon: "Cascade confirmation: PPV 88.5%"
- Bottom: warning box: "T2DM adjustment mandatory — lowers TRG by 0.681 units"
- Background: light green (#E8F8F2)

**Bottom strip:** CALON-FH registry | UK Biobank Wales n=7,253 | UK Biobank n=58,021 | 4,028 FH cases | Wales validation AUC 0.842
**Colour treatment:** Clean white panels with coloured header bars. Arrow connectors in #0072B2.

---

## CANVA PROMPT 2 — LDL RECEPTOR BIOLOGY + TRG SHIELD MECHANISM
### Purpose: Manuscript Figure panel or supplementary educational figure (183×120mm)

**Layout:** Two-column comparison diagram

**LEFT COLUMN — "Normal LDLR function"**
- Hepatocyte cell membrane (grey wavy line at top)
- LDLR protein structure standing upright:
  - 7 red circles (R1–R7 ligand-binding domain) stacked
  - 2 yellow rectangles (EGF-A, EGF-B)
  - Blue diamond (beta-propeller — pH sensor)
  - Green rectangle (EGF-C)
  - Purple bar (transmembrane)
- LDL particle (brown/orange sphere with "LDL" label) approaching R domain
- Arrow: LDL particle entering cell (endocytosis)
- Arrow down from liver to blood: "↓ LDL-C (NORMAL)"
- VLDL secretion from liver upward: "TG secretion normal"
- Label: "TRG Filter HIGH"
- TRG formula box: "LDL / (TG+0.1) = HIGH" in green

**RIGHT COLUMN — "LDLR dysfunction (FH)"**
- Same hepatocyte membrane
- LDLR with red X on R domain (mutation site)
- LDL particles floating in blood NOT being cleared
- Red arrow: "LDL clearance ↓↓"
- VLDL secretion still intact (green arrow): "VLDL/TG: PRESERVED"
- Blood space shows elevated LDL particles
- Label: "TRG Filter HIGH in FH" (paradoxically high)
- TRG box: "LDL_elevated / (TG_normal+0.1) = HIGH" in blue

**BOTTOM SECTION — T2DM Confounder**
- Third panel showing T2DM effect:
  - Adipocyte releasing FFAs → liver increases VLDL → TG elevated
  - Result: TRG filter LOWERED even though LDL still elevated
  - Warning triangle: "T2DM raises TG by ~0.5–1.0 mmol/L → false-low TRG filter"
  - β = −0.48 to −0.66 (p=8.2×10⁻⁷⁰) label

**Colour:** White background, coloured domain blocks as above. Nature Medicine style.

---

## CANVA PROMPT 3 — CLINICAL DECISION ALGORITHM (FLOWCHART)
### Purpose: Clinical guidance figure for the Discussion section (89×140mm portrait)

**TITLE:** "TUDOR Score — Clinical Decision Algorithm"

**FLOWCHART NODES (top to bottom):**

Node 1 (entry, orange): "Patient referred to lipid clinic\nor population LDL screen"
↓ arrow

Node 2 (blue rounded rect): "Calculate TUDOR score\n[LDL_untreated, TG, ApoB, age, sex, statin]"
↓

Diamond decision: "TUDOR score?"
→ LEFT branch (TUDOR < 0.003): "Score <0.003 (Sensitivity 94% threshold)"
   → Box (green): "FH unlikely\nNPV 99.9%\nRoutine follow-up"

→ MIDDLE branch (0.003–0.045): "Score 0.003–0.045"
   → Box (yellow): "Intermediate risk\nRepeat lipid profile\nConsider cascade screening"

→ RIGHT branch (TUDOR ≥ 0.045): "Score ≥0.045 (Youden threshold)"
   → Box (orange): "TUDOR positive\nSensitivity 64.2%, Specificity 81.1%\nPPV 4.3% (lipid clinic)"
   ↓
   Diamond: "Genetic testing available?"
   → YES → Box (blue): "LDLR/APOB/PCSK9 sequencing\nDefinitive diagnosis"
   → NO  → Box (red): "Clinical FH criteria (DLCN/Simon Broome)\nIntensify statin therapy\nFamily cascade screening"

**Additional note boxes on right side:**
- "If cascade screening:\nUse threshold 0.159\n(Specificity 95%, PPV 88.5%)"
- "T2DM adjustment:\nCheck T2DM status\nT2DM raises TG → may lower TUDOR"
- "Treatment paradox:\nStatins do NOT invalidate score\n(uses untreated-equivalent LDL)"

**Footer:** "TUDOR score available at [calculator URL] | Validation: n=6,448, AUC 0.842"
**Colour:** Canva flowchart template in blue/orange palette. Arial font 10pt minimum.

---

## CANVA PROMPT 4 — KEY FINDINGS SUMMARY INFOGRAPHIC
### Purpose: Visual abstract for social media / journal website (1080×1080px square)

**LAYOUT:** 2×2 grid of findings cards on dark blue (#0A2A4A) background

**CARD 1 (top left) — Performance**
- Large white text: "AUC 0.842"
- Sub: "External validation, South Wales"
- Comparison bar (white): TUDOR ████████▌ 0.842
- Comparison bar (grey): DLCN  ██████░░░ 0.791
- Bottom: "Wales AUC 0.842, UKB AUC 0.750 | DLCN 0.791 / 0.636"
- Accent colour: #56B4E9 (sky blue)

**CARD 2 (top right) — The TRG Shield**
- Icon: liver + TG molecule
- Large text: "Cohen's d = 1.21"
- Sub: "LDLR vs Non-FH separation"
- Small molecule diagram: LDL↑ ÷ TG↔ = TRG Filter ↑
- Caption: "Novel pathophysiology-informed feature"
- Accent colour: #009E73 (green)

**CARD 3 (bottom left) — T2DM Alert**
- Warning triangle icon
- Large text: "β = −0.48 to −0.66"
- Sub: "T2DM lowers TRG Filter"
- p-value: "3 definitions: self-report, ICD-10, HbA1c ≥48"
- Quote: "T2DM adjustment is mandatory"
- Bar comparison: TRG filter T2DM− 2.75 ██████ vs T2DM+ 2.07 █████
- Accent colour: #D55E00 (vermillion)

**CARD 4 (bottom right) — Three Clinical Uses**
- Title: "One score. Three thresholds."
- Table:
  | Use | Threshold | Key metric |
  | Population screen | 0.003 | Sens 94% |
  | Lipid clinic | 0.045 | Youden |
  | Cascade confirm | 0.159 | Spec 95% |
- Bottom: "Wales n=7,253 | UK Biobank n=58,021 | 4,028 FH cases | UK Biobank | Wales"
- Accent colour: #E69F00 (orange)

**Central element:** TUDOR logo in white bold with "Machine Learning FH Detection"
**Font:** Montserrat Bold (headers), Open Sans (body)

---

## CANVA PROMPT 5 — FIGURE 1 EDITORIAL COVER ILLUSTRATION
### Purpose: Potential journal cover image or highlighted figure (183×240mm portrait)

**CONCEPT:** "Finding the needle in the haystack" — visualising FH detection in a lipid clinic population

**COMPOSITION:**

**BACKGROUND:** Deep midnight blue gradient (#0A2A4A → #1A3A5C)

**FOREGROUND ELEMENTS:**

Layer 1 — Population scatter (bottom 60% of image):
- 500 small white/grey dots representing non-FH patients
- 15–20 bright orange dots representing FH cases hidden within
- Dots arranged in a realistic scatter cloud pattern
- Most orange dots visually blending into the white dots

Layer 2 — TUDOR score filter (middle 20%):
- Horizontal gradient sweep bar (blue to green) across the image
- Label: "TUDOR score threshold = 0.045"
- Below the bar: grey dots remain (correctly excluded, grey)
- Above the bar: orange dots clearly highlighted (correctly identified)
- A few grey dots above the bar (false positives, shown as outline circles)

Layer 3 — Result panel (top 20%):
- White panel overlay
- Text: "64 correctly identified" (bold orange, large)
- Small text: "of 100 true FH cases in lipid clinic"
- ROC curve miniature thumbnail in white panel
- AUC badge: "0.842" in large blue bold

**Typography:**
- Main title: "TUDOR" in Montserrat ExtraBold, 72pt, white
- Subtitle: "Precision FH Detection" in 24pt white
- All other text: Source Sans Pro

**Bottom credit strip:** "CALON-FH | UK Biobank | TUDOR machine learning diagnostic score"

**Art direction:** Elegant, data-driven aesthetic. No clip art. Clean geometric shapes.
The image should feel like it belongs on the cover of Nature Medicine or Lancet.

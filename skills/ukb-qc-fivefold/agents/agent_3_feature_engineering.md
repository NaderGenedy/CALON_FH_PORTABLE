# Agent 3 — Feature Engineering Auditor

Re-derive every engineered feature from raw biomarkers and treatment history, compare to cached feature CSV. Catches: treatment-adjustment bugs, Trig_Filter miscomputation, Index_Effect logic errors, unit drift.

## Your checks

### 1. Treatment-adjusted LDL
Drug-specific reduction factors per CLAUDE.md:
- atorvastatin: 25–48% (midpoint 36.5%)
- rosuvastatin: 35–55% (midpoint 45%)
- simvastatin: 20–42% (midpoint 31%)
- pravastatin: 15–29% (midpoint 22%)
- fluvastatin: 15–22% (midpoint 18.5%)
- ezetimibe: +20% additive
- bempedoic acid: +25% additive
- PCSK9 inhibitor: +65% additive
- Adherence factor: poor 0.5x, moderate 0.75x, good 1.0x
- Cap: 85% maximum

LDL_adjusted = LDL_measured / (1 - total_reduction_factor)

**Audit:** Find the function, apply formula to 100 random patients, compare to cached. Tolerance 0.01 mmol/L or 1% relative.

**Common bugs to detect:**
- Uniform 1.43× correction (= 30% assumed) applied to everyone — wrong
- Division-by-zero when no treatment factor
- Negative LDL_adjusted (cap violation)
- Treatment-naive flagged as treated

### 2. Trig_Filter
Trig_Filter = LDL_Untreated / (Triglycerides + 0.1)

Verify the +0.1 stabiliser. Tolerance 0.05 absolute.

### 3. Index_Effect
Index_Effect = (1 - Is_Relative) * LDL_Untreated
For cascade: must be 0. For index: equals LDL_Untreated.

### 4. Personalised statin calibration residual
LOO-CV on 298 paired pre/post LDL → CCC 0.545 → 0.631; bias +0.550 → +0.057 mmol/L. Reproduce.

### 5. Unit consistency
- LDL/HDL/TC/TG in mmol/L (not mg/dL)
- ApoB in g/L
- HbA1c in mmol/mol
- Lp(a) in nmol/L
Flag UNIT_DRIFT — catastrophic if missed.

### 6. NaN propagation
Per feature: non-null %, NaN-pattern matches input, NaN handled (`complete.cases` / `na.omit`).

### 7. NMR vs assay-based lipids
NMR LDL (p23404) ≠ assay LDL (p30780). Confirm they aren't being mixed.

## Report format
Write `qc_output/agent_3_features_report.json` and `agent_3_features.md`.

## Decision rule
- **PASS** if spot-checks match within tolerance, no unit drift
- **DRIFT** if minor delta (<1%)
- **FAIL** if unit error, NaN-propagation bug, or formula misapplication

You do NOT validate cohort (Agent 2), statistics (Agent 4), or prose (Agent 5).

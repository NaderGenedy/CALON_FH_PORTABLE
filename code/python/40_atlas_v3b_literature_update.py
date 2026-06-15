#!/usr/bin/env python3
"""
40_atlas_v3b_literature_update.py
==================================
Literature-informed update of the LDLR genotype-phenotype atlas.

Incorporates:
  - Reimund et al. Nature 2025: cryo-EM confirmed domain boundaries (~4 Å) and
    binding-site functional data (Extended Data Table 2)
  - Jensson et al. NEJM 2023 (Iceland deCODE): lifespan effect estimates per variant

New columns added to BOTH output tables:
  domain_v3b              — Updated domain label (Reimund 2025 boundaries)
  mechanistic_class       — binding_failure | recycling_failure | interface_BS2 |
                            structural_misfolding | recycling_pivot
  binding_site            — BS1 | BS2 | none
  ldl_contact             — True / False
  p_notation              — HGVS protein (3-letter AA codes)
  c_notation              — Approximate HGVS cDNA
  jensson_lifespan_years  — Effect on lifespan (Jensson 2023), NaN if not reported
  jensson_significant     — True if p < 0.05 in Jensson 2023
  reimund_binding_site    — BS1/BS2/propeller label from Reimund 2025 Ext. Data
  reimund_functional_effect — Functional consequence text from Reimund 2025

Outputs:
  alphafold/analysis/LDLR_genotype_phenotype_atlas_v3b.csv
  alphafold/analysis/LDLR_860_position_confidence_table_v3b.csv

Author  : Dr Nader Genedy
Version : v3b — 2026-03-25
"""

import os
import sys
import numpy as np
import pandas as pd
from collections import defaultdict

# Force UTF-8 output so Unicode characters in print() work on Windows terminals
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ("utf-8", "utf8"):
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ── paths ─────────────────────────────────────────────────────────────────────
BASE     = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")

ATLAS_IN  = os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas.csv")
POS_IN    = os.path.join(ANALYSIS, "LDLR_860_position_confidence_table.csv")
ATLAS_OUT = os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas_v3b.csv")
POS_OUT   = os.path.join(ANALYSIS, "LDLR_860_position_confidence_table_v3b.csv")

print("=" * 80)
print("LDLR ATLAS v3b — LITERATURE-INFORMED UPDATE")
print("Reimund et al. Nature 2025 + Jensson et al. NEJM 2023")
print("=" * 80)

# ── 1-letter → 3-letter amino acid mapping ────────────────────────────────────
AA1TO3 = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    '*': 'Ter', 'X': 'Xaa',
}

# ── DOMAIN BOUNDARIES (Reimund et al. Nature 2025, ~4 Å cryo-EM) ─────────────
# Format: (start, end, label)
DOMAIN_BOUNDARIES_V3B = [
    (1,   21,  "Signal peptide"),
    (22,  65,  "LA1"),
    (66,  106, "LA2"),
    (107, 145, "LA3"),
    (146, 186, "LA4"),
    (187, 233, "LA5"),
    (234, 272, "LA6"),
    (274, 313, "LA7"),         # Note: position 273 is in the gap — see below
    (314, 353, "EGF-A"),
    (354, 393, "EGF-B"),
    (394, 663, "Beta-propeller"),
    (664, 712, "EGF-C"),
    (713, 749, "O-linked sugars"),
    (750, 770, "Transmembrane"),
    (771, 860, "Cytoplasmic"),
]

# pH-switch residues (Reimund 2025 confirmed)
PH_SWITCH_RESIDUES = {211, 285, 583, 607}

# LDL-contact domains (Reimund 2025 BS1/BS2 designation)
# BS1: LA3-LA7 + EGF-A (direct LDL contact at binding site 1)
BS1_DOMAINS   = {"LA3", "LA4", "LA5", "LA7", "EGF-A"}
# LA6 and EGF-B: minimal/no direct contact per Reimund 2025
# BS2: Beta-propeller contacts LDL N-terminal domain (NTD)
BS2_DOMAINS   = {"Beta-propeller"}
LDL_CONTACT_DOMAINS = BS1_DOMAINS | BS2_DOMAINS  # LA5+LA7 pH switch included

# Recycling pivot domain
RECYCLING_PIVOT_DOMAINS = {"EGF-B"}

# Domains with NO or minimal LDL contact
NO_CONTACT_DOMAINS = {
    "Signal peptide", "LA1", "LA2", "LA6",
    "EGF-B", "EGF-C", "O-linked sugars", "Transmembrane", "Cytoplasmic",
}

# ── Build position → domain_v3b lookup ───────────────────────────────────────
def assign_domain_v3b(position: int) -> str:
    """Return the Reimund 2025 domain label for a given position (1-860)."""
    for start, end, label in DOMAIN_BOUNDARIES_V3B:
        if start <= position <= end:
            return label
    # Position 273 falls in the gap between LA6 (234-272) and LA7 (274-313)
    # Conservatively assigned to the adjacent LA6/LA7 linker
    if position == 273:
        return "LA6"   # gap residue — treat with LA6 characteristics
    # Anything beyond 860 should not exist, but guard defensively
    return "Cytoplasmic"

DOMAIN_V3B_MAP = {pos: assign_domain_v3b(pos) for pos in range(1, 861)}

# ── Mechanistic classification ────────────────────────────────────────────────
def mechanistic_class(position: int, domain: str,
                      at_pcsk9_interface: bool,
                      interface_distance) -> str:
    """
    Priority order:
      1. recycling_failure — pH-switch residues (experimental Reimund 2025)
      2. interface_BS2     — beta-propeller + PCSK9 interface or dist < 8 Å
      3. binding_failure   — BS1 domains (LA3-LA7, EGF-A)
      4. recycling_pivot   — EGF-B
      5. structural_misfolding — everything else
    """
    if position in PH_SWITCH_RESIDUES:
        return "recycling_failure"

    if domain in BS2_DOMAINS:
        try:
            dist = float(interface_distance)
        except (TypeError, ValueError):
            dist = float("inf")
        if at_pcsk9_interface or dist < 8.0:
            return "interface_BS2"
        return "structural_misfolding"   # propeller residue, not at interface

    if domain in BS1_DOMAINS:
        return "binding_failure"

    if domain in RECYCLING_PIVOT_DOMAINS:
        return "recycling_pivot"

    # LA1, LA2, LA6, EGF-C, O-linked, TM, Cytoplasmic, Signal peptide
    return "structural_misfolding"

# ── Binding site annotation ────────────────────────────────────────────────────
def binding_site(domain: str) -> str:
    if domain in BS1_DOMAINS:
        return "BS1"
    if domain in BS2_DOMAINS:
        return "BS2"
    return "none"

# ── LDL contact ───────────────────────────────────────────────────────────────
def ldl_contact(domain: str) -> bool:
    return domain in LDL_CONTACT_DOMAINS

# ── p. notation ───────────────────────────────────────────────────────────────
def make_p_notation(wt_aa: str, position: int, mut_aa: str) -> str:
    """
    Returns e.g. p.Cys352Tyr
    Uses 3-letter codes; handles stop codon as Ter.
    Returns empty string if wt or mut are missing/unknown.
    """
    wt  = str(wt_aa).strip() if pd.notna(wt_aa) else ""
    mut = str(mut_aa).strip() if pd.notna(mut_aa) else ""
    if not wt or not mut or wt in ("X", "") or mut in ("X", ""):
        return ""
    wt3  = AA1TO3.get(wt.upper(),  wt.upper())
    mut3 = AA1TO3.get(mut.upper(), mut.upper())
    return f"p.{wt3}{position}{mut3}"

# ── c. notation (approximate — position only, no exact nucleotide) ────────────
def make_c_notation(position: int) -> str:
    """
    Approximate cDNA position.
    Codon start position in CDS = (position - 1) * 3 + 1
    (position is protein position; no exact nucleotide substitution available)
    Format: c.<codon_start> (e.g. c.1054 for protein position 352)
    """
    codon_start = (position - 1) * 3 + 1
    return f"c.{codon_start}"

# ── Jensson et al. NEJM 2023 lifespan data ────────────────────────────────────
# Key: (position, wt_aa, mut_aa) using 1-letter codes; value: (years, significant)
JENSSON_DATA = {
    # c.694+2T>C — splice variant, protein position approximated at 231
    # (intronic splice donor at end of exon 4, last codon ~p.Cys231)
    (231, "C", "splice"): (-6.47, True),
    # p.Tyr576Ser
    (576, "Y", "S"):       (-8.70, False),
    # p.Asp307Asn
    (307, "D", "N"):       (-3.72, False),
    # p.Cys231Ter — separate entry from the splice variant
    (231, "C", "*"):       (+2.72, False),   # NS, possible survivor bias
    # p.Asp707Val
    (707, "D", "V"):       (-2.56, False),
    # p.Ala540Thr
    (540, "A", "T"):       (-0.17, False),
}

# Lookup by position+mutation only (for variant tables without alt allele col)
JENSSON_BY_POS_MUT = {}
for (pos, wt, mut), (yrs, sig) in JENSSON_DATA.items():
    JENSSON_BY_POS_MUT[(pos, mut)] = (yrs, sig)
    JENSSON_BY_POS_MUT[(pos, wt, mut)] = (yrs, sig)

# ── Reimund et al. Nature 2025 functional data ────────────────────────────────
# Key: (position, wt_aa, mut_aa) 1-letter codes
# Value: (reimund_binding_site_label, reimund_functional_effect)
REIMUND_DATA = {
    (131, "D", "H"): ("BS1/LA3",     "N.D."),
    (149, "G", "R"): ("BS1/LA3",     "N.D."),
    (160, "C", "R"): ("BS1/LA4",     "N.D."),
    (161, "I", "T"): ("BS1/LA4",     "N.D."),
    (165, "W", "G"): ("BS1/LA4",     "decreased LDL binding"),
    (172, "D", "N"): ("BS1/LA4",     "N.D."),
    (211, "H", "Y"): ("BS1/LA5",     "decreased LDL binding"),
    (214, "W", "R"): ("BS1/LA5",     "decreased LDL uptake"),
    (221, "D", "N"): ("BS1/LA5",     "N.D."),
    (303, "R", "W"): ("BS1/LA7",     "N.D."),
    (308, "E", "K"): ("BS1/LA7",     "N.D."),
    (406, "R", "P"): ("BS2/propeller","N.D."),
    (406, "R", "W"): ("BS2/propeller","decreased LDL binding and uptake"),
    (428, "N", "K"): ("BS2/propeller","N.D."),
    (520, "R", "G"): ("BS2/propeller","N.D."),
    (583, "H", "Y"): ("BS2/propeller","decreased LDL uptake"),
    (583, "H", "R"): ("BS2/propeller","N.D."),
    (607, "H", "D"): ("BS2/propeller","decreased LDL binding"),
    (625, "N", "K"): ("BS2/propeller","decreased LDL binding and uptake"),
}

# ─────────────────────────────────────────────────────────────────────────────
# LOAD INPUT FILES
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n[1/4] Loading atlas v3a from: {ATLAS_IN}")
if not os.path.exists(ATLAS_IN):
    sys.exit(f"ERROR: Atlas file not found: {ATLAS_IN}")
atlas = pd.read_csv(ATLAS_IN)
print(f"      Loaded {len(atlas)} rows, {atlas.shape[1]} columns")

print(f"\n[2/4] Loading position confidence table from: {POS_IN}")
if not os.path.exists(POS_IN):
    sys.exit(f"ERROR: Position table not found: {POS_IN}")
pos_table = pd.read_csv(POS_IN)
print(f"      Loaded {len(pos_table)} rows, {pos_table.shape[1]} columns")

# Normalise column name differences between the two tables
# Atlas uses 'wildtype_aa'; position table also uses 'wildtype_aa'
# Atlas uses 'pcsk9_interface_distance'; position table uses 'pcsk9_interface_dist'
def get_pcsk9_dist(row: pd.Series, table_type: str):
    if table_type == "atlas":
        return row.get("pcsk9_interface_distance", np.nan)
    return row.get("pcsk9_interface_dist", np.nan)

def get_at_pcsk9(row: pd.Series, table_type: str):
    val = row.get("at_pcsk9_interface", False)
    if isinstance(val, bool):
        return val
    return str(val).strip().lower() in ("yes", "true", "1")

# ─────────────────────────────────────────────────────────────────────────────
# HELPER: APPLY LITERATURE COLUMNS TO A DATAFRAME
# ─────────────────────────────────────────────────────────────────────────────

def apply_literature_columns(df: pd.DataFrame, table_type: str) -> pd.DataFrame:
    """
    Adds all new columns to either the atlas or the position table.
    table_type: "atlas" | "position"
    """
    df = df.copy()

    positions = df["position"].astype(int).values
    wt_aas    = df["wildtype_aa"].astype(str).str.strip().values

    # ── domain_v3b ────────────────────────────────────────────────────────────
    df["domain_v3b"] = [DOMAIN_V3B_MAP[p] for p in positions]

    # ── Columns that require per-row logic ────────────────────────────────────
    mech_classes   = []
    binding_sites  = []
    ldl_contacts   = []

    for i, pos in enumerate(positions):
        dom  = DOMAIN_V3B_MAP[pos]
        row  = df.iloc[i]
        dist = get_pcsk9_dist(row, table_type)
        at_pcsk9 = get_at_pcsk9(row, table_type)

        mech_classes.append(mechanistic_class(pos, dom, at_pcsk9, dist))
        binding_sites.append(binding_site(dom))
        ldl_contacts.append(ldl_contact(dom))

    df["mechanistic_class"] = mech_classes
    df["binding_site"]      = binding_sites
    df["ldl_contact"]       = ldl_contacts

    # ── p. and c. notation ───────────────────────────────────────────────────
    # These are per-position summaries (wildtype only); variant-specific notation
    # is computed where mut_aa is available in the atlas.
    if table_type == "atlas":
        # Atlas has no mut_aa column — leave p_notation blank at position level
        df["p_notation"] = ""
        df["c_notation"] = [make_c_notation(p) for p in positions]
    else:
        # Position table: same — position-level cDNA annotation
        df["p_notation"] = ""
        df["c_notation"] = [make_c_notation(p) for p in positions]

    # ── Jensson lifespan data ─────────────────────────────────────────────────
    # For position tables we can only populate rows where the position + mutation
    # uniquely matches. Here we annotate by position only (using the first match
    # if multiple mutations at the same position exist in Jensson 2023).
    jensson_years = []
    jensson_sig   = []

    # Build a position-level default (NaN / False)
    for pos in positions:
        # Collect all Jensson entries for this position
        matches = [(yrs, sig) for (k, val) in JENSSON_BY_POS_MUT.items()
                   if isinstance(k, tuple) and len(k) >= 2 and k[0] == pos
                   for yrs, sig in [val]]
        # Deduplicate (tuple keys have lengths 2 and 3, some overlap)
        seen = set()
        unique_matches = []
        for k, val in JENSSON_BY_POS_MUT.items():
            if isinstance(k, tuple) and k[0] == pos and val not in seen:
                seen.add(val)
                unique_matches.append(val)

        if unique_matches:
            # Use most extreme (largest absolute) effect if multiple
            best = min(unique_matches, key=lambda x: x[0])   # most negative
            jensson_years.append(best[0])
            jensson_sig.append(best[1])
        else:
            jensson_years.append(np.nan)
            jensson_sig.append(np.nan)

    df["jensson_lifespan_years"] = jensson_years
    df["jensson_significant"]    = [bool(s) if not (isinstance(s, float) and np.isnan(s)) else np.nan
                                    for s in jensson_sig]

    # ── Reimund functional data ────────────────────────────────────────────────
    reimund_bs    = []
    reimund_func  = []

    for pos in positions:
        # Collect all Reimund entries at this position
        matches = [(v[0], v[1]) for (k, v) in REIMUND_DATA.items() if k[0] == pos]
        if matches:
            # Join if multiple mutations at same position
            bs_vals   = "; ".join(sorted(set(m[0] for m in matches)))
            func_vals = "; ".join(sorted(set(m[1] for m in matches)))
            reimund_bs.append(bs_vals)
            reimund_func.append(func_vals)
        else:
            reimund_bs.append("")
            reimund_func.append("")

    df["reimund_binding_site"]     = reimund_bs
    df["reimund_functional_effect"] = reimund_func

    return df


# ─────────────────────────────────────────────────────────────────────────────
# APPLY TO BOTH TABLES
# ─────────────────────────────────────────────────────────────────────────────
print("\n[3/4] Applying literature updates...")

atlas_v3b     = apply_literature_columns(atlas,     "atlas")
pos_table_v3b = apply_literature_columns(pos_table, "position")

# For the ATLAS table, also build variant-level p_notation where possible.
# The atlas is a POSITION-level table (no mut_aa), so p_notation stays blank
# unless there is a stored variant string. We annotate the c_notation column
# at minimum (useful for lookup).
# No further change needed — c_notation is already populated per position.

# ─────────────────────────────────────────────────────────────────────────────
# SAVE OUTPUTS
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n[4/4] Saving updated tables...")

atlas_v3b.to_csv(ATLAS_OUT, index=False)
print(f"      Atlas v3b saved    : {ATLAS_OUT}")
print(f"      Rows: {len(atlas_v3b)}, Columns: {atlas_v3b.shape[1]}")

pos_table_v3b.to_csv(POS_OUT, index=False)
print(f"      Position table v3b : {POS_OUT}")
print(f"      Rows: {len(pos_table_v3b)}, Columns: {pos_table_v3b.shape[1]}")

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY STATISTICS
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("SUMMARY STATISTICS — LDLR ATLAS v3b")
print("=" * 80)

# ── A. Domain boundary comparison (old vs new) ────────────────────────────────
print("\n--- A. DOMAIN BOUNDARY COMPARISON (v3a to v3b) ---")
old_counts = atlas["domain"].value_counts().sort_index()
new_counts = atlas_v3b["domain_v3b"].value_counts().sort_index()

all_domains = sorted(set(list(old_counts.index) + list(new_counts.index)))
print(f"  {'Domain':<35s}  {'v3a (old)':>9s}  {'v3b (new)':>9s}  {'Delta':>6s}")
print(f"  {'-'*35}  {'-'*9}  {'-'*9}  {'-'*6}")
for d in all_domains:
    old_n = int(old_counts.get(d, 0))
    new_n = int(new_counts.get(d, 0))
    delta = new_n - old_n
    delta_str = f"{delta:+d}" if delta != 0 else "  —"
    print(f"  {d:<35s}  {old_n:9d}  {new_n:9d}  {delta_str:>6s}")

# ── B. Variants per mechanistic class ─────────────────────────────────────────
print("\n--- B. POSITIONS PER MECHANISTIC CLASS ---")
mech_counts = atlas_v3b["mechanistic_class"].value_counts()
for cls in ["binding_failure", "recycling_failure", "interface_BS2",
            "recycling_pivot", "structural_misfolding"]:
    n = int(mech_counts.get(cls, 0))
    pct = 100.0 * n / len(atlas_v3b)
    print(f"  {cls:<25s}: {n:4d} positions ({pct:5.1f}%)")

# ── C. Variants per binding site ──────────────────────────────────────────────
print("\n--- C. POSITIONS PER BINDING SITE ---")
bs_counts = atlas_v3b["binding_site"].value_counts()
for bs in ["BS1", "BS2", "none"]:
    n = int(bs_counts.get(bs, 0))
    pct = 100.0 * n / len(atlas_v3b)
    print(f"  {bs:<6s}: {n:4d} positions ({pct:5.1f}%)")

# ── D. LDL contact distribution ───────────────────────────────────────────────
n_contact     = int(atlas_v3b["ldl_contact"].sum())
n_no_contact  = len(atlas_v3b) - n_contact
print(f"\n--- D. LDL CONTACT DISTRIBUTION ---")
print(f"  LDL contact    : {n_contact:4d} positions ({100*n_contact/len(atlas_v3b):.1f}%)")
print(f"  No LDL contact : {n_no_contact:4d} positions ({100*n_no_contact/len(atlas_v3b):.1f}%)")

# ── E. Mean ddG per mechanistic class ─────────────────────────────────────────
print("\n--- E. MEAN ddG (kcal/mol) PER MECHANISTIC CLASS ---")
ddg_col = "observed_max_ddG" if "observed_max_ddG" in atlas_v3b.columns else None
if ddg_col:
    ddg_by_class = (atlas_v3b
                    .dropna(subset=[ddg_col])
                    .groupby("mechanistic_class")[ddg_col]
                    .agg(["mean", "median", "count"]))
    for cls, row in ddg_by_class.iterrows():
        print(f"  {cls:<25s}: mean={row['mean']:+6.2f}, median={row['median']:+6.2f}, "
              f"n={int(row['count'])} positions with ddG data")
else:
    print("  (observed_max_ddG column not found)")

# ── F. Mean penetrance per mechanistic class ──────────────────────────────────
pen_col = "domain_penetrance_by_60_pct" if "domain_penetrance_by_60_pct" in atlas_v3b.columns else None
print("\n--- F. MEAN PENETRANCE (% ASCVD BY AGE 60) PER MECHANISTIC CLASS ---")
if pen_col:
    pen_by_class = atlas_v3b.groupby("mechanistic_class")[pen_col].agg(["mean", "count"])
    for cls, row in pen_by_class.iterrows():
        print(f"  {cls:<25s}: mean={row['mean']:5.1f}%, n={int(row['count'])} positions")
else:
    print("  (domain_penetrance_by_60_pct column not found)")

# ── G. Correlation: ddG vs mechanistic class severity ─────────────────────────
print("\n--- G. CORRELATION: ddG vs MECHANISTIC CLASS SEVERITY ---")
MECH_SEVERITY = {
    "recycling_failure":   4,
    "interface_BS2":       3,
    "binding_failure":     3,
    "recycling_pivot":     2,
    "structural_misfolding": 1,
}
atlas_v3b["mech_severity_score"] = atlas_v3b["mechanistic_class"].map(MECH_SEVERITY)

if ddg_col:
    corr_data = atlas_v3b[[ddg_col, "mech_severity_score"]].dropna()
    if len(corr_data) >= 10:
        corr = corr_data[ddg_col].corr(corr_data["mech_severity_score"])
        print(f"  Pearson r(ddG, severity_score) = {corr:.3f}  (n={len(corr_data)} positions)")
    else:
        print(f"  Insufficient data for correlation (n={len(corr_data)})")
else:
    print("  (ddG column unavailable)")

# ── H. Jensson data ───────────────────────────────────────────────────────────
print("\n--- H. JENSSON et al. 2023 LIFESPAN DATA ---")
jensson_annotated_atlas = atlas_v3b["jensson_lifespan_years"].notna().sum()
jensson_annotated_pos   = pos_table_v3b["jensson_lifespan_years"].notna().sum()
print(f"  Positions annotated with Jensson lifespan effect (atlas)         : {jensson_annotated_atlas}")
print(f"  Positions annotated with Jensson lifespan effect (position table): {jensson_annotated_pos}")
jensson_subset = atlas_v3b.dropna(subset=["jensson_lifespan_years"])
if len(jensson_subset) > 0:
    print(f"\n  {'Pos':>4s}  {'Domain v3b':<18s}  {'Mech class':<22s}  {'Years':>7s}  {'Sig':>5s}")
    print(f"  {'----':>4s}  {'------------------':<18s}  {'----------------------':<22s}  {'-------':>7s}  {'-----':>5s}")
    for _, r in jensson_subset.sort_values("jensson_lifespan_years").iterrows():
        sig_str = "YES" if r.get("jensson_significant") is True else "NS"
        print(f"  {int(r['position']):4d}  {str(r['domain_v3b']):<18s}  "
              f"{str(r['mechanistic_class']):<22s}  "
              f"{r['jensson_lifespan_years']:+7.2f}  {sig_str:>5s}")

# ── I. Reimund functional data ─────────────────────────────────────────────────
print("\n--- I. REIMUND et al. 2025 FUNCTIONAL DATA ---")
reimund_annotated_atlas = (atlas_v3b["reimund_binding_site"] != "").sum()
reimund_annotated_pos   = (pos_table_v3b["reimund_binding_site"] != "").sum()
print(f"  Positions with Reimund binding-site annotation (atlas)          : {reimund_annotated_atlas}")
print(f"  Positions with Reimund binding-site annotation (position table) : {reimund_annotated_pos}")

# Breakdown by functional effect
func_effects = [e for e in atlas_v3b["reimund_functional_effect"] if e != ""]
from collections import Counter
effect_counts = Counter(func_effects)
print(f"\n  Reimund functional effect distribution (position level):")
for eff, cnt in sorted(effect_counts.items(), key=lambda x: -x[1]):
    print(f"    {eff:<45s}: {cnt:3d} positions")

# ── J. New column overview ────────────────────────────────────────────────────
print("\n--- J. NEW COLUMNS ADDED (v3a to v3b) ---")
new_cols = ["domain_v3b", "mechanistic_class", "binding_site", "ldl_contact",
            "p_notation", "c_notation",
            "jensson_lifespan_years", "jensson_significant",
            "reimund_binding_site", "reimund_functional_effect",
            "mech_severity_score"]
for col in new_cols:
    if col in atlas_v3b.columns:
        n_non_null = atlas_v3b[col].notna().sum() if atlas_v3b[col].dtype != object else (atlas_v3b[col] != "").sum()
        print(f"  {col:<35s}: present, {n_non_null} / {len(atlas_v3b)} non-empty")
    else:
        print(f"  {col:<35s}: MISSING — check script logic")

# ── Final ─────────────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("ATLAS v3b BUILD COMPLETE")
print(f"  Atlas   : {ATLAS_OUT}")
print(f"  Positions: {POS_OUT}")
print(f"  New columns: {len(new_cols)} added to each table")
print(f"  Domain boundaries: Reimund et al. Nature 2025 (~4 Å cryo-EM)")
print(f"  Lifespan data    : Jensson et al. NEJM 2023 (Iceland deCODE)")
print(f"  Functional data  : Reimund et al. Nature 2025 Ext. Data Table 2")
print("=" * 80)

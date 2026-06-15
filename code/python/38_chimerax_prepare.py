#!/usr/bin/env python3
"""
38_chimerax_prepare.py
Generates ChimeraX attribute files (.defattr) and individual figure command scripts
from the LDLR genotype-phenotype atlas for publication-quality structural visualizations.

Usage:
    python 38_chimerax_prepare.py

Outputs:
    - alphafold/analysis/figures/chimerax/*.defattr  (attribute files)
    - alphafold/analysis/figures/chimerax/figure_CX01.cxc through figure_CX10.cxc
"""

import csv
import os
import math

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ATLAS_CSV = os.path.join(BASE_DIR, "alphafold", "analysis", "LDLR_genotype_phenotype_atlas.csv")
PLDDT_CSV = os.path.join(BASE_DIR, "alphafold", "analysis", "LDLR_wildtype_plddt_per_residue.csv")
OUT_DIR = os.path.join(BASE_DIR, "alphafold", "analysis", "figures", "chimerax")
FIG_DIR = OUT_DIR

# PDB paths (relative to BASE_DIR for portability in CXC scripts)
LDLR_PDB = os.path.join(BASE_DIR, "alphafold", "foldx", "LDLR_wt_Repair.pdb")
PCSK9_PDB = os.path.join(BASE_DIR, "alphafold", "foldx", "PCSK9_wt_Repair.pdb")
LDLR_PCSK9_PDB = os.path.join(BASE_DIR, "alphafold", "foldx", "LDLR_PCSK9_complex_v2.pdb")
LDLR_APOB_PDB = os.path.join(BASE_DIR, "alphafold", "foldx", "LDLR_ApoB_complex.pdb")

# Interface positions
PCSK9_INTERFACE = [25, 35, 50, 53, 58, 59, 60, 61, 78, 79, 81, 96,
                   316, 318, 319, 320, 321, 322, 325, 326, 328, 329,
                   330, 331, 332, 339, 341, 342, 344, 351, 352, 372,
                   386, 387, 390]
APOB_INTERFACE = [285, 286, 298, 300, 303, 304, 305, 308, 309, 310, 311, 312]

# Domain ranges
DOMAINS = {
    "Signal peptide": (1, 21),
    "Ligand-binding R1": (22, 83),
    "Ligand-binding R2": (84, 120),
    "Ligand-binding R3": (121, 159),
    "Ligand-binding R4": (160, 197),
    "Ligand-binding R5": (198, 236),
    "Ligand-binding R6": (237, 275),
    "Ligand-binding R7": (276, 313),
    "EGF-like A": (314, 354),
    "EGF-like B": (355, 394),
    "Beta-propeller": (395, 632),
    "EGF-like C": (633, 672),
    "O-linked sugar": (693, 750),
    "Transmembrane": (751, 788),
    "Cytoplasmic": (789, 860),
}

# Domain color palette (Nature-calibre)
DOMAIN_COLORS = {
    "Signal peptide": "#999999",
    "Ligand-binding R1": "#1B4F72",
    "Ligand-binding R2": "#2471A3",
    "Ligand-binding R3": "#2E86C1",
    "Ligand-binding R4": "#3498DB",
    "Ligand-binding R5": "#5DADE2",
    "Ligand-binding R6": "#85C1E9",
    "Ligand-binding R7": "#AED6F1",
    "EGF-like A": "#C0392B",
    "EGF-like B": "#E74C3C",
    "Beta-propeller": "#7D3C98",
    "EGF-like C": "#EC7063",
    "O-linked sugar": "#F4D03F",
    "Transmembrane": "#8B6914",
    "Cytoplasmic": "#27AE60",
}

# Drug target category colors (9 categories)
DRUG_CATEGORY_COLORS = {
    "Protein trafficking / folding": "#E64B35",
    "LDL-binding enhancement": "#4DBBD5",
    "PCSK9 interaction / receptor recycling": "#00A087",
    "Protein stability / trafficking": "#3C5488",
    "Receptor recycling stabilization": "#F39B7F",
    "Structural flexibility / folding": "#8491B4",
    "pH-dependent release / recycling": "#91D1C2",
    "Membrane anchoring / trafficking": "#B09C85",
    "Endocytic signaling / PCSK9i target": "#7E6148",
}

# Risk classification colors
RISK_COLORS = {
    "Critical": "#CC0000",
    "High": "#FF6600",
    "Moderate": "#FFB300",
    "Low": "#00AA44",
}


def read_atlas():
    """Read the atlas CSV and return list of dicts."""
    rows = []
    with open(ATLAS_CSV, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def write_defattr_categorical(filepath, attr_name, data_dict):
    """
    Write a categorical .defattr file.
    data_dict: {position_int: category_string}
    """
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(f"attribute: {attr_name}\n")
        f.write("match mode: 1-to-1\n")
        f.write("recipient: residues\n")
        for pos in sorted(data_dict.keys()):
            val = data_dict[pos]
            f.write(f"\t:{pos}.A\t{val}\n")


def write_defattr_continuous(filepath, attr_name, data_dict):
    """
    Write a continuous .defattr file.
    data_dict: {position_int: float_value}  (None values are skipped)
    """
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(f"attribute: {attr_name}\n")
        f.write("match mode: 1-to-1\n")
        f.write("recipient: residues\n")
        for pos in sorted(data_dict.keys()):
            val = data_dict[pos]
            if val is not None:
                f.write(f"\t:{pos}.A\t{val:.4f}\n")


def generate_defattr_files(atlas_rows):
    """Generate all .defattr files from atlas data."""
    risk_data = {}
    plddt_data = {}
    penetrance_data = {}
    ddg_data = {}
    n_patients_data = {}
    drug_cat_data = {}

    for row in atlas_rows:
        pos = int(row["position"])
        risk_data[pos] = row["risk_classification"]
        drug_cat_data[pos] = row["drug_target_category"]

        try:
            plddt_data[pos] = float(row["plddt_confidence"])
        except (ValueError, KeyError):
            plddt_data[pos] = None

        try:
            penetrance_data[pos] = float(row["domain_penetrance_by_60_pct"])
        except (ValueError, KeyError):
            penetrance_data[pos] = None

        try:
            val = row["observed_max_ddG"]
            ddg_data[pos] = float(val) if val else None
        except (ValueError, KeyError):
            ddg_data[pos] = None

        try:
            val = row["n_patients"]
            n_patients_data[pos] = int(val) if val else 0
        except (ValueError, KeyError):
            n_patients_data[pos] = 0

    # Write risk classification
    write_defattr_categorical(
        os.path.join(OUT_DIR, "risk_classification.defattr"),
        "risk_level", risk_data
    )

    # Write drug target category
    write_defattr_categorical(
        os.path.join(OUT_DIR, "drug_target_category.defattr"),
        "drug_category", drug_cat_data
    )

    # Write pLDDT
    write_defattr_continuous(
        os.path.join(OUT_DIR, "plddt_confidence.defattr"),
        "plddt", plddt_data
    )

    # Write penetrance
    write_defattr_continuous(
        os.path.join(OUT_DIR, "domain_penetrance.defattr"),
        "penetrance", penetrance_data
    )

    # Write ddG
    write_defattr_continuous(
        os.path.join(OUT_DIR, "observed_max_ddG.defattr"),
        "ddG", ddg_data
    )

    # Write n_patients as continuous for heatmap
    patients_float = {k: float(v) if v is not None else None for k, v in n_patients_data.items()}
    write_defattr_continuous(
        os.path.join(OUT_DIR, "n_patients.defattr"),
        "npatients", patients_float
    )

    print(f"[OK] Generated .defattr files in {OUT_DIR}")


def _pdb_path_cxc(path):
    """Convert Windows path to forward-slash for ChimeraX."""
    return path.replace("\\", "/")


def generate_figure_cx01():
    """CX01: LDLR Risk Atlas Overview."""
    lines = []
    lines.append("# Figure CX01: LDLR Risk Atlas Overview")
    lines.append(f"open {_pdb_path_cxc(LDLR_PDB)}")
    lines.append("set bgColor white")
    lines.append("lighting soft")
    lines.append("graphics silhouettes true width 1.5")
    lines.append("cartoon style protein arrows true")
    lines.append("hide atoms")
    lines.append("show cartoons")
    lines.append("")
    lines.append("# Color all residues by risk classification")

    atlas = read_atlas()
    risk_groups = {"Critical": [], "High": [], "Moderate": [], "Low": []}
    for row in atlas:
        pos = int(row["position"])
        risk = row["risk_classification"]
        if risk in risk_groups:
            risk_groups[risk].append(pos)

    for risk, color in RISK_COLORS.items():
        positions = risk_groups.get(risk, [])
        if positions:
            pos_str = ",".join(str(p) for p in positions)
            lines.append(f"color /A:{pos_str} {color}")

    lines.append("")
    lines.append("# Domain labels at midpoints")
    domain_labels = {
        "SP": 11, "R1": 52, "R2": 102, "R3": 140, "R4": 178,
        "R5": 217, "R6": 256, "R7": 294, "EGF-A": 334,
        "EGF-B": 374, "Beta-prop": 513, "EGF-C": 652,
        "O-link": 721, "TM": 769, "Cyto": 824,
    }
    for label, pos in domain_labels.items():
        lines.append(f'label /A:{pos} text "{label}" height 1.2 color black offset 0,2,0')

    lines.append("")
    lines.append("view")
    lines.append(f'save {_pdb_path_cxc(os.path.join(FIG_DIR, "Figure_CX01_risk_atlas.png"))} width 3000 height 2000 supersample 4')

    with open(os.path.join(OUT_DIR, "figure_CX01.cxc"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def generate_figure_cx02():
    """CX02: pLDDT Confidence Map."""
    atlas = read_atlas()
    lines = []
    lines.append("# Figure CX02: pLDDT Confidence Map")
    lines.append(f"open {_pdb_path_cxc(LDLR_PDB)}")
    lines.append("set bgColor white")
    lines.append("lighting soft")
    lines.append("graphics silhouettes true width 1.5")
    lines.append("cartoon style protein arrows true")
    lines.append("hide atoms")
    lines.append("show cartoons")
    lines.append("")
    lines.append("# Color by pLDDT using AlphaFold convention")
    lines.append("# Very high (>90): blue #0053D6, High (70-90): cyan #65CBF3")
    lines.append("# Low (50-70): yellow #FFDB13, Very low (<50): orange #FF7D45")

    very_high = []
    high = []
    low = []
    very_low = []
    for row in atlas:
        pos = int(row["position"])
        try:
            plddt = float(row["plddt_confidence"])
        except (ValueError, KeyError):
            very_low.append(pos)
            continue
        if plddt > 90:
            very_high.append(pos)
        elif plddt > 70:
            high.append(pos)
        elif plddt > 50:
            low.append(pos)
        else:
            very_low.append(pos)

    if very_high:
        lines.append(f"color /A:{','.join(str(p) for p in very_high)} #0053D6")
    if high:
        lines.append(f"color /A:{','.join(str(p) for p in high)} #65CBF3")
    if low:
        lines.append(f"color /A:{','.join(str(p) for p in low)} #FFDB13")
    if very_low:
        lines.append(f"color /A:{','.join(str(p) for p in very_low)} #FF7D45")

    lines.append("")
    lines.append("view")
    lines.append(f'save {_pdb_path_cxc(os.path.join(FIG_DIR, "Figure_CX02_plddt_view1.png"))} width 3000 height 2000 supersample 4')
    lines.append("")
    lines.append("# Second orientation: top view")
    lines.append("turn y 90")
    lines.append(f'save {_pdb_path_cxc(os.path.join(FIG_DIR, "Figure_CX02_plddt_view2.png"))} width 3000 height 2000 supersample 4')

    with open(os.path.join(OUT_DIR, "figure_CX02.cxc"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def generate_figure_cx03():
    """CX03: PCSK9 Interface Highlight."""
    lines = []
    lines.append("# Figure CX03: PCSK9 Interface Highlight")
    lines.append(f"open {_pdb_path_cxc(LDLR_PCSK9_PDB)}")
    lines.append("set bgColor white")
    lines.append("lighting soft")
    lines.append("graphics silhouettes true width 1.5")
    lines.append("cartoon style protein arrows true")
    lines.append("hide atoms")
    lines.append("show cartoons")
    lines.append("")
    lines.append("# LDLR (chain A) in light grey cartoon")
    lines.append("color /A #C8C8C8")
    lines.append("")
    lines.append("# PCSK9 (chain C) as transparent light blue surface")
    lines.append("color /C #ADD8E6")
    lines.append("surface /C")
    lines.append("transparency /C 70 target s")
    lines.append("")
    lines.append("# Highlight PCSK9 interface residues on LDLR in magenta spheres")
    pcsk9_pos = ",".join(str(p) for p in PCSK9_INTERFACE)
    lines.append(f"select /A:{pcsk9_pos}")
    lines.append("show sel atoms")
    lines.append("style sel sphere")
    lines.append("color sel #CC00CC")
    lines.append("~select")
    lines.append("")
    lines.append("# Label key interface residues")
    key_residues = [25, 50, 316, 320, 325, 330, 339, 344, 351, 386]
    for pos in key_residues:
        aa = None
        atlas = read_atlas()
        for row in atlas:
            if int(row["position"]) == pos:
                aa = row["wildtype_aa"]
                break
        if aa:
            lines.append(f'label /A:{pos} text "{aa}{pos}" height 1.0 color black offset 0,1.5,0')

    lines.append("")
    lines.append("view")
    lines.append(f'save {_pdb_path_cxc(os.path.join(FIG_DIR, "Figure_CX03_pcsk9_interface.png"))} width 3000 height 2000 supersample 4')

    with open(os.path.join(OUT_DIR, "figure_CX03.cxc"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def generate_figure_cx04():
    """CX04: ApoB Interface Highlight."""
    lines = []
    lines.append("# Figure CX04: ApoB Interface Highlight")
    lines.append(f"open {_pdb_path_cxc(LDLR_APOB_PDB)}")
    lines.append("set bgColor white")
    lines.append("lighting soft")
    lines.append("graphics silhouettes true width 1.5")
    lines.append("cartoon style protein arrows true")
    lines.append("hide atoms")
    lines.append("show cartoons")
    lines.append("")
    lines.append("# LDLR (chain A) in light grey cartoon")
    lines.append("color /A #C8C8C8")
    lines.append("")
    lines.append("# ApoB (chain C) as transparent light green surface")
    lines.append("color /C #90EE90")
    lines.append("surface /C")
    lines.append("transparency /C 70 target s")
    lines.append("")
    lines.append("# Highlight ApoB interface residues on LDLR in cyan spheres")
    apob_pos = ",".join(str(p) for p in APOB_INTERFACE)
    lines.append(f"select /A:{apob_pos}")
    lines.append("show sel atoms")
    lines.append("style sel sphere")
    lines.append("color sel #00CED1")
    lines.append("~select")
    lines.append("")
    lines.append("# Label key interface residues")
    atlas = read_atlas()
    for pos in APOB_INTERFACE:
        aa = None
        for row in atlas:
            if int(row["position"]) == pos:
                aa = row["wildtype_aa"]
                break
        if aa:
            lines.append(f'label /A:{pos} text "{aa}{pos}" height 1.0 color black offset 0,1.5,0')

    lines.append("")
    lines.append("view")
    lines.append(f'save {_pdb_path_cxc(os.path.join(FIG_DIR, "Figure_CX04_apob_interface.png"))} width 3000 height 2000 supersample 4')

    with open(os.path.join(OUT_DIR, "figure_CX04.cxc"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def generate_figure_cx05():
    """CX05: Zero-Overlap Dual Interface."""
    lines = []
    lines.append("# Figure CX05: Zero-Overlap Dual Interface")
    lines.append(f"open {_pdb_path_cxc(LDLR_PDB)}")
    lines.append("set bgColor white")
    lines.append("lighting soft")
    lines.append("graphics silhouettes true width 1.5")
    lines.append("cartoon style protein arrows true")
    lines.append("hide atoms")
    lines.append("show cartoons")
    lines.append("")
    lines.append("# Base color: light grey")
    lines.append("color /A #E0E0E0")
    lines.append("")
    lines.append("# PCSK9 interface residues in magenta")
    pcsk9_pos = ",".join(str(p) for p in PCSK9_INTERFACE)
    lines.append(f"color /A:{pcsk9_pos} #CC00CC")
    lines.append(f"select /A:{pcsk9_pos}")
    lines.append("show sel atoms")
    lines.append("style sel sphere")
    lines.append("color sel #CC00CC")
    lines.append("~select")
    lines.append("")
    lines.append("# ApoB interface residues in cyan")
    apob_pos = ",".join(str(p) for p in APOB_INTERFACE)
    lines.append(f"color /A:{apob_pos} #00CED1")
    lines.append(f"select /A:{apob_pos}")
    lines.append("show sel atoms")
    lines.append("style sel sphere")
    lines.append("color sel #00CED1")
    lines.append("~select")
    lines.append("")
    lines.append("# Labels for interface clusters")
    lines.append('label /A:330 text "PCSK9 interface" height 1.5 color #CC00CC offset 0,3,0')
    lines.append('label /A:300 text "ApoB interface" height 1.5 color #00CED1 offset 0,3,0')
    lines.append("")
    lines.append("# Distance measurement between interface cluster centroids")
    lines.append("# Measure from a representative PCSK9 interface residue to ApoB interface residue")
    lines.append("distance /A:330@CA /A:300@CA")
    lines.append("")
    lines.append("view")
    lines.append(f'save {_pdb_path_cxc(os.path.join(FIG_DIR, "Figure_CX05_zero_overlap.png"))} width 3000 height 2000 supersample 4')

    with open(os.path.join(OUT_DIR, "figure_CX05.cxc"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def generate_figure_cx06():
    """CX06: Domain Architecture."""
    lines = []
    lines.append("# Figure CX06: Domain Architecture")
    lines.append(f"open {_pdb_path_cxc(LDLR_PDB)}")
    lines.append("set bgColor white")
    lines.append("lighting soft")
    lines.append("graphics silhouettes true width 1.5")
    lines.append("cartoon style protein arrows true")
    lines.append("hide atoms")
    lines.append("show cartoons")
    lines.append("")
    lines.append("# Color each domain")
    for domain, (start, end) in DOMAINS.items():
        color = DOMAIN_COLORS[domain]
        lines.append(f"color /A:{start}-{end} {color}")

    lines.append("")
    lines.append("# Semi-transparent surface with cartoon inside")
    lines.append("surface /A")
    lines.append("transparency /A 75 target s")
    lines.append("")
    lines.append("# Domain labels")
    domain_labels = {
        "SP": 11, "R1": 52, "R2": 102, "R3": 140, "R4": 178,
        "R5": 217, "R6": 256, "R7": 294, "EGF-A": 334,
        "EGF-B": 374, "Beta-prop": 513, "EGF-C": 652,
        "O-link": 721, "TM": 769, "Cyto": 824,
    }
    for label, pos in domain_labels.items():
        lines.append(f'label /A:{pos} text "{label}" height 1.2 color black offset 0,2,0')

    lines.append("")
    lines.append("view")
    lines.append(f'save {_pdb_path_cxc(os.path.join(FIG_DIR, "Figure_CX06_domain_architecture.png"))} width 3000 height 2000 supersample 4')

    with open(os.path.join(OUT_DIR, "figure_CX06.cxc"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def generate_figure_cx07():
    """CX07: Critical Residues Close-up."""
    critical_positions = [316, 318, 319, 320, 321, 322, 325, 326, 328, 329,
                          330, 331, 332, 339, 341, 342, 344, 351, 352]

    atlas = read_atlas()
    pos_to_aa = {}
    for row in atlas:
        pos = int(row["position"])
        if pos in critical_positions:
            pos_to_aa[pos] = row["wildtype_aa"]

    lines = []
    lines.append("# Figure CX07: Critical Residues Close-up")
    lines.append(f"open {_pdb_path_cxc(LDLR_PCSK9_PDB)}")
    lines.append("set bgColor white")
    lines.append("lighting soft")
    lines.append("graphics silhouettes true width 1.5")
    lines.append("cartoon style protein arrows true")
    lines.append("hide atoms")
    lines.append("show cartoons")
    lines.append("")
    lines.append("# LDLR in light grey, PCSK9 in light blue")
    lines.append("color /A #D0D0D0")
    lines.append("color /C #ADD8E6")
    lines.append("")
    lines.append("# Show critical residues as ball-and-stick in red")
    crit_pos = ",".join(str(p) for p in critical_positions)
    lines.append(f"select /A:{crit_pos}")
    lines.append("show sel atoms")
    lines.append("style sel ball")
    lines.append("color sel #CC0000")
    lines.append("~select")
    lines.append("")
    lines.append("# Show PCSK9 contact residues near the interface")
    lines.append("show /C atoms")
    lines.append("style /C ball")
    lines.append("color /C #6699CC")
    lines.append("hide /C cartoons")
    lines.append("show /C cartoons")
    lines.append("hide /C atoms")
    lines.append("# Show PCSK9 residues close to interface")
    lines.append("select zone /A:{} 5 /C".format(crit_pos))
    lines.append("show sel atoms")
    lines.append("style sel ball")
    lines.append("color sel #6699CC target a")
    lines.append("~select")
    lines.append("")
    lines.append("# Label each critical residue")
    for pos in critical_positions:
        aa = pos_to_aa.get(pos, "X")
        lines.append(f'label /A:{pos} text "{aa}{pos}" height 0.8 color black offset 0,1,0')

    lines.append("")
    lines.append("# Zoom to the EGF-A region")
    lines.append(f"view /A:{crit_pos}")
    lines.append(f'save {_pdb_path_cxc(os.path.join(FIG_DIR, "Figure_CX07_critical_closeup.png"))} width 3000 height 2000 supersample 4')

    with open(os.path.join(OUT_DIR, "figure_CX07.cxc"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def generate_figure_cx08():
    """CX08: Pathogenic Hotspot Heatmap (n_patients)."""
    atlas = read_atlas()
    lines = []
    lines.append("# Figure CX08: Pathogenic Hotspot Heatmap")
    lines.append(f"open {_pdb_path_cxc(LDLR_PDB)}")
    lines.append("set bgColor white")
    lines.append("lighting soft")
    lines.append("graphics silhouettes true width 1.5")
    lines.append("cartoon style protein arrows true")
    lines.append("hide atoms")
    lines.append("show cartoons")
    lines.append("")
    lines.append("# Color by number of patients at each position")
    lines.append("# White = 0 patients, yellow = low, orange = medium, red = high")

    # Build bins
    zero_pos = []
    low_pos = []      # 1-5
    med_pos = []      # 6-20
    high_pos = []     # 21-50
    vhigh_pos = []    # >50

    for row in atlas:
        pos = int(row["position"])
        try:
            np_val = int(row["n_patients"])
        except (ValueError, KeyError):
            np_val = 0
        if np_val == 0:
            zero_pos.append(pos)
        elif np_val <= 5:
            low_pos.append(pos)
        elif np_val <= 20:
            med_pos.append(pos)
        elif np_val <= 50:
            high_pos.append(pos)
        else:
            vhigh_pos.append(pos)

    if zero_pos:
        lines.append(f"color /A:{','.join(str(p) for p in zero_pos)} #F0F0F0")
    if low_pos:
        lines.append(f"color /A:{','.join(str(p) for p in low_pos)} #FFEDA0")
    if med_pos:
        lines.append(f"color /A:{','.join(str(p) for p in med_pos)} #FEB24C")
    if high_pos:
        lines.append(f"color /A:{','.join(str(p) for p in high_pos)} #F03B20")
    if vhigh_pos:
        lines.append(f"color /A:{','.join(str(p) for p in vhigh_pos)} #BD0026")

    lines.append("")
    lines.append("view")
    lines.append(f'save {_pdb_path_cxc(os.path.join(FIG_DIR, "Figure_CX08_hotspot_heatmap.png"))} width 3000 height 2000 supersample 4')

    with open(os.path.join(OUT_DIR, "figure_CX08.cxc"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def generate_figure_cx09():
    """CX09: FoldX ddG Structural Map."""
    atlas = read_atlas()
    lines = []
    lines.append("# Figure CX09: FoldX ddG Structural Map")
    lines.append(f"open {_pdb_path_cxc(LDLR_PDB)}")
    lines.append("set bgColor white")
    lines.append("lighting soft")
    lines.append("graphics silhouettes true width 1.5")
    lines.append("cartoon style protein arrows true")
    lines.append("hide atoms")
    lines.append("show cartoons")
    lines.append("")
    lines.append("# Color by observed_max_ddG")
    lines.append("# Blue = stabilizing (ddG < -1), white = neutral, red = destabilizing (ddG > 2)")
    lines.append("# Transparent grey = no data")

    no_data = []
    stabilizing = []     # ddG < -1
    mild_stab = []       # -1 <= ddG < 0
    neutral = []         # 0 <= ddG < 1
    mild_destab = []     # 1 <= ddG < 2
    destabilizing = []   # 2 <= ddG < 5
    severe = []          # ddG >= 5

    for row in atlas:
        pos = int(row["position"])
        val = row["observed_max_ddG"]
        if not val:
            no_data.append(pos)
            continue
        try:
            ddg = float(val)
        except ValueError:
            no_data.append(pos)
            continue
        if ddg < -1:
            stabilizing.append(pos)
        elif ddg < 0:
            mild_stab.append(pos)
        elif ddg < 1:
            neutral.append(pos)
        elif ddg < 2:
            mild_destab.append(pos)
        elif ddg < 5:
            destabilizing.append(pos)
        else:
            severe.append(pos)

    if no_data:
        lines.append(f"color /A:{','.join(str(p) for p in no_data)} #CCCCCC")
    if stabilizing:
        lines.append(f"color /A:{','.join(str(p) for p in stabilizing)} #2166AC")
    if mild_stab:
        lines.append(f"color /A:{','.join(str(p) for p in mild_stab)} #67A9CF")
    if neutral:
        lines.append(f"color /A:{','.join(str(p) for p in neutral)} #F7F7F7")
    if mild_destab:
        lines.append(f"color /A:{','.join(str(p) for p in mild_destab)} #FDDBC7")
    if destabilizing:
        lines.append(f"color /A:{','.join(str(p) for p in destabilizing)} #EF8A62")
    if severe:
        lines.append(f"color /A:{','.join(str(p) for p in severe)} #B2182B")

    # Make no-data positions semi-transparent
    if no_data:
        lines.append("# Make no-data positions semi-transparent")
        lines.append(f"transparency /A:{','.join(str(p) for p in no_data)} 60 target c")

    lines.append("")
    lines.append("view")
    lines.append(f'save {_pdb_path_cxc(os.path.join(FIG_DIR, "Figure_CX09_ddG_map.png"))} width 3000 height 2000 supersample 4')

    with open(os.path.join(OUT_DIR, "figure_CX09.cxc"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def generate_figure_cx10():
    """CX10: Drug Target Categories."""
    atlas = read_atlas()
    lines = []
    lines.append("# Figure CX10: Drug Target Categories")
    lines.append(f"open {_pdb_path_cxc(LDLR_PDB)}")
    lines.append("set bgColor white")
    lines.append("lighting soft")
    lines.append("graphics silhouettes true width 1.5")
    lines.append("cartoon style protein arrows true")
    lines.append("hide atoms")
    lines.append("show cartoons")
    lines.append("")
    lines.append("# Color by drug target category (9 categories)")

    cat_groups = {}
    for row in atlas:
        pos = int(row["position"])
        cat = row["drug_target_category"]
        if cat not in cat_groups:
            cat_groups[cat] = []
        cat_groups[cat].append(pos)

    for cat, positions in cat_groups.items():
        color = DRUG_CATEGORY_COLORS.get(cat, "#888888")
        pos_str = ",".join(str(p) for p in positions)
        lines.append(f"# {cat}")
        lines.append(f"color /A:{pos_str} {color}")

    lines.append("")
    lines.append("# Labels for each drug target region")
    # Place a label at the midpoint of each category's first contiguous block
    cat_label_positions = {
        "Protein trafficking / folding": 11,
        "LDL-binding enhancement": 140,
        "PCSK9 interaction / receptor recycling": 334,
        "Protein stability / trafficking": 374,
        "Receptor recycling stabilization": 513,
        "Structural flexibility / folding": 652,
        "pH-dependent release / recycling": 425,
        "Membrane anchoring / trafficking": 769,
        "Endocytic signaling / PCSK9i target": 824,
    }
    for cat, pos in cat_label_positions.items():
        short_label = cat.split("/")[0].strip()
        if len(short_label) > 20:
            short_label = short_label[:18] + ".."
        lines.append(f'label /A:{pos} text "{short_label}" height 1.0 color black offset 0,2,0')

    lines.append("")
    lines.append("view")
    lines.append(f'save {_pdb_path_cxc(os.path.join(FIG_DIR, "Figure_CX10_drug_targets.png"))} width 3000 height 2000 supersample 4')

    with open(os.path.join(OUT_DIR, "figure_CX10.cxc"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def main():
    # Create output directory
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"Output directory: {OUT_DIR}")

    # Read atlas
    atlas_rows = read_atlas()
    print(f"Read {len(atlas_rows)} residue entries from atlas")

    # Generate .defattr files
    generate_defattr_files(atlas_rows)

    # Generate individual figure scripts
    generate_figure_cx01()
    print("[OK] Generated figure_CX01.cxc")
    generate_figure_cx02()
    print("[OK] Generated figure_CX02.cxc")
    generate_figure_cx03()
    print("[OK] Generated figure_CX03.cxc")
    generate_figure_cx04()
    print("[OK] Generated figure_CX04.cxc")
    generate_figure_cx05()
    print("[OK] Generated figure_CX05.cxc")
    generate_figure_cx06()
    print("[OK] Generated figure_CX06.cxc")
    generate_figure_cx07()
    print("[OK] Generated figure_CX07.cxc")
    generate_figure_cx08()
    print("[OK] Generated figure_CX08.cxc")
    generate_figure_cx09()
    print("[OK] Generated figure_CX09.cxc")
    generate_figure_cx10()
    print("[OK] Generated figure_CX10.cxc")

    print(f"\nAll files generated in {OUT_DIR}")
    print("Run 38_chimerax_visualizations.cxc in ChimeraX to generate all figures,")
    print("or run individual figure_CX##.cxc scripts for specific figures.")


if __name__ == "__main__":
    main()

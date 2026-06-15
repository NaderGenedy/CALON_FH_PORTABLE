#!/usr/bin/env python3
"""
37_foldx_saturation_mutagenesis.py
===================================
Systematic saturation mutagenesis using FoldX BuildModel on the AlphaFold3-predicted
LDLR wildtype structure. Computes ddG for ALL possible single amino acid substitutions
at ALL 860 positions (860 x 19 = 16,340 mutations).

Outputs:
  alphafold/analysis/foldx_saturation_mutagenesis.csv
  alphafold/analysis/LDLR_genotype_phenotype_atlas_v2.csv  (updated with position-specific ddG)
"""

import pandas as pd
import numpy as np
import os
import sys
import subprocess
import shutil
import time
import glob as globmod

# ── paths ─────────────────────────────────────────────────────────────────────
BASE = r"C:/Users/nader/Downloads/calon_ukb_pipeline"
ANALYSIS = os.path.join(BASE, "alphafold", "analysis")

# FoldX paths
FOLDX_BIN = r"C:\Users\nader\Downloads\foldxWindows\foldx_20261231.exe"
FOLDX_ROTABASE = os.path.join(BASE, "alphafold", "foldx", "rotabase.txt")
WT_PDB = os.path.join(BASE, "alphafold", "foldx", "LDLR_wt_Repair.pdb")

# Working directory for FoldX runs — use E: (external SSD) to avoid filling C: drive
WORK_DIR = r"E:/CALON_FH_BACKUP/foldx_saturation_work"
if not os.path.exists("E:/CALON_FH_BACKUP"):
    # Fallback to original location if SSD not attached
    WORK_DIR = os.path.join(BASE, "foldx_saturation_work")

# Output files
OUTPUT_CSV = os.path.join(ANALYSIS, "foldx_saturation_mutagenesis.csv")
ATLAS_V1 = os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas.csv")
ATLAS_V2 = os.path.join(ANALYSIS, "LDLR_genotype_phenotype_atlas_v2.csv")

# FoldX parameters
NUM_RUNS = 3
BATCH_SIZE = 20  # mutations per FoldX call to avoid memory issues

# All 20 standard amino acids (single letter)
ALL_AA = list("ACDEFGHIKLMNPQRSTVWY")

AA3TO1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}


def classify_effect(ddg):
    """Classify ddG into effect categories."""
    if pd.isna(ddg):
        return "unknown"
    if ddg < -1.0:
        return "stabilizing"
    elif ddg <= 1.0:
        return "neutral"
    elif ddg <= 2.0:
        return "mild"
    else:
        return "destabilizing"


def read_wt_residues(pdb_path):
    """Read wildtype residue identities from PDB (CA atoms only, chain A)."""
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                chain = line[21].strip()
                if chain != 'A':
                    continue
                resname = line[17:20].strip()
                resnum = int(line[22:26].strip())
                aa = AA3TO1.get(resname, None)
                if aa is not None:
                    residues[resnum] = aa
    return residues


def generate_all_mutations(wt_residues):
    """Generate all 860 x 19 = 16,340 possible single amino acid substitutions."""
    mutations = []
    for pos in sorted(wt_residues.keys()):
        wt_aa = wt_residues[pos]
        for mut_aa in ALL_AA:
            if mut_aa == wt_aa:
                continue
            # FoldX notation: {WT_AA}{chain}{position}{MUT_AA}
            notation = f"{wt_aa}A{pos}{mut_aa}"
            mutations.append({
                'position': pos,
                'wt_aa': wt_aa,
                'mut_aa': mut_aa,
                'foldx_notation': notation,
            })
    return mutations


def load_completed_mutations(output_csv):
    """Load already-computed mutations for resume capability."""
    if not os.path.exists(output_csv):
        return set()
    try:
        df = pd.read_csv(output_csv)
        completed = set(df['foldx_notation'].tolist())
        return completed
    except Exception:
        return set()


def write_individual_list(mutations_batch, list_file):
    """Write FoldX individual_list.txt file for a batch of mutations.

    Format: one mutation per line, terminated with semicolon.
    Each line: {WT_AA}{chain}{position}{MUT_AA};
    """
    with open(list_file, 'w') as f:
        for mut in mutations_batch:
            f.write(f"{mut['foldx_notation']};\n")


def run_foldx_batch(mutations_batch, work_dir, pdb_name):
    """Run FoldX BuildModel for a batch of mutations.

    Returns list of dicts with ddG results, or None for failures.
    """
    list_file = os.path.join(work_dir, "individual_list.txt")
    write_individual_list(mutations_batch, list_file)

    # Run FoldX BuildModel
    cmd = [
        FOLDX_BIN,
        "--command=BuildModel",
        f"--pdb={pdb_name}",
        f"--mutant-file={list_file}",
        f"--numberOfRuns={NUM_RUNS}",
        "--out-pdb=false",
    ]

    try:
        result = subprocess.run(
            cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=1800,  # 30 min timeout per batch (increased for difficult residues)
        )
    except subprocess.TimeoutExpired:
        print(f"    [TIMEOUT] FoldX timed out for this batch")
        return None
    except Exception as e:
        print(f"    [ERROR] FoldX execution failed: {e}")
        return None

    # Parse Average output file
    # FoldX creates Average_<pdb_name_without_ext>.fxout
    pdb_stem = pdb_name.replace('.pdb', '')
    avg_file = os.path.join(work_dir, f"Average_{pdb_stem}.fxout")

    if not os.path.exists(avg_file):
        # Try alternative naming patterns
        avg_candidates = globmod.glob(os.path.join(work_dir, "Average_*.fxout"))
        if avg_candidates:
            avg_file = avg_candidates[0]
        else:
            print(f"    [ERROR] No Average output file found")
            return None

    return parse_average_fxout(avg_file, mutations_batch)


def parse_average_fxout(avg_file, mutations_batch):
    """Parse the Average_*.fxout file to extract ddG values.

    The file has a header, then one row per mutation in the batch.
    Row naming: {pdb_stem}_{mutation_index} (1-based)
    Column 'total energy' (index 2, after Pdb and SD) is the ddG.
    """
    results = []
    try:
        with open(avg_file) as f:
            lines = f.readlines()

        # Find the header line (starts with "Pdb")
        header_idx = None
        for i, line in enumerate(lines):
            if line.startswith('Pdb\t'):
                header_idx = i
                break

        if header_idx is None:
            print(f"    [ERROR] Could not find header in {avg_file}")
            return None

        data_lines = []
        for line in lines[header_idx + 1:]:
            line = line.strip()
            if line and not line.startswith('#'):
                data_lines.append(line)

        if len(data_lines) != len(mutations_batch):
            print(f"    [WARNING] Expected {len(mutations_batch)} results, got {len(data_lines)}")

        for i, line in enumerate(data_lines):
            if i >= len(mutations_batch):
                break
            fields = line.split('\t')
            # Fields: Pdb, SD, total energy, ...
            # total energy is the ddG value
            try:
                ddg = float(fields[2])  # total energy column
            except (IndexError, ValueError):
                ddg = np.nan

            mut = mutations_batch[i].copy()
            mut['ddG'] = ddg
            mut['effect_class'] = classify_effect(ddg)
            results.append(mut)

        return results

    except Exception as e:
        print(f"    [ERROR] Failed to parse {avg_file}: {e}")
        return None


def cleanup_intermediate_files(work_dir, pdb_stem):
    """Remove intermediate PDB files and temporary FoldX output to save disk space."""
    patterns = [
        f"{pdb_stem}_*.pdb",
        f"WT_{pdb_stem}_*.pdb",
        f"Dif_{pdb_stem}.fxout",
        f"Raw_{pdb_stem}.fxout",
        f"PdbList_{pdb_stem}.fxout",
        f"Average_{pdb_stem}.fxout",
    ]
    for pattern in patterns:
        for fpath in globmod.glob(os.path.join(work_dir, pattern)):
            try:
                os.remove(fpath)
            except OSError:
                pass


def format_time(seconds):
    """Format seconds to human readable string."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        return f"{seconds / 60:.1f}min"
    else:
        h = int(seconds // 3600)
        m = int((seconds % 3600) // 60)
        return f"{h}h {m}min"


def progress_bar(current, total, width=50, prefix=""):
    """Print a simple progress bar."""
    pct = current / total if total > 0 else 1.0
    filled = int(width * pct)
    bar = "#" * filled + "-" * (width - filled)
    sys.stdout.write(f"\r   {prefix}[{bar}] {current}/{total} ({100*pct:.1f}%)")
    sys.stdout.flush()


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN EXECUTION
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 80)
    print("FOLDX SATURATION MUTAGENESIS - LDLR (860 positions x 19 mutations)")
    print("=" * 80)

    # ── 1. Read wildtype residues from PDB ────────────────────────────────────
    print("\n[1/6] Reading wildtype residues from PDB...")
    wt_residues = read_wt_residues(WT_PDB)
    n_residues = len(wt_residues)
    print(f"   Found {n_residues} residues (range {min(wt_residues)}-{max(wt_residues)})")

    if n_residues == 0:
        print("   ERROR: No residues found in PDB. Aborting.")
        sys.exit(1)

    # ── 2. Generate all possible mutations ────────────────────────────────────
    print("\n[2/6] Generating all possible single amino acid substitutions...")
    all_mutations = generate_all_mutations(wt_residues)
    total_mutations = len(all_mutations)
    print(f"   Total mutations to compute: {total_mutations}")
    print(f"   ({n_residues} positions x 19 substitutions each)")

    # ── 3. Check for resume capability ────────────────────────────────────────
    print("\n[3/6] Checking for previously computed mutations (resume)...")
    completed = load_completed_mutations(OUTPUT_CSV)
    n_completed = len(completed)

    if n_completed > 0:
        print(f"   Found {n_completed} previously computed mutations")
        # Filter out already-computed mutations
        remaining_mutations = [m for m in all_mutations if m['foldx_notation'] not in completed]
        print(f"   Remaining mutations to compute: {len(remaining_mutations)}")
    else:
        print(f"   No previous results found - starting from scratch")
        remaining_mutations = all_mutations

    if len(remaining_mutations) == 0:
        print("   All mutations already computed! Skipping to atlas rebuild.")
    else:
        # ── 4. Setup working directory ────────────────────────────────────────
        print("\n[4/6] Setting up FoldX working directory...")
        os.makedirs(WORK_DIR, exist_ok=True)

        # Copy PDB to working directory
        pdb_name = os.path.basename(WT_PDB)
        pdb_dest = os.path.join(WORK_DIR, pdb_name)
        if not os.path.exists(pdb_dest):
            shutil.copy2(WT_PDB, pdb_dest)
            print(f"   Copied PDB to working directory")

        # Copy rotabase.txt to working directory (required by FoldX)
        rotabase_dest = os.path.join(WORK_DIR, "rotabase.txt")
        if not os.path.exists(rotabase_dest):
            shutil.copy2(FOLDX_ROTABASE, rotabase_dest)
            print(f"   Copied rotabase.txt to working directory")

        # Also copy molecules directory if it exists (FoldX5 may need it)
        molecules_src = os.path.join(os.path.dirname(FOLDX_BIN), "molecules")
        molecules_dest = os.path.join(WORK_DIR, "molecules")
        if os.path.isdir(molecules_src) and not os.path.isdir(molecules_dest):
            shutil.copytree(molecules_src, molecules_dest)
            print(f"   Copied molecules directory to working directory")

        print(f"   Working directory: {WORK_DIR}")

        # ── 5. Run FoldX in batches ───────────────────────────────────────────
        print(f"\n[5/6] Running FoldX BuildModel in batches of {BATCH_SIZE}...")
        n_remaining = len(remaining_mutations)
        n_batches = (n_remaining + BATCH_SIZE - 1) // BATCH_SIZE

        # Time estimation: ~10-30 seconds per mutation with 3 runs
        est_seconds_per_mut = 15  # conservative estimate
        est_total = est_seconds_per_mut * n_remaining
        print(f"   Estimated total time: {format_time(est_total)} "
              f"(~{est_seconds_per_mut}s per mutation)")
        print(f"   Number of batches: {n_batches}")
        print()

        pdb_stem = pdb_name.replace('.pdb', '')
        all_results = []
        n_success = 0
        n_fail = 0
        start_time = time.time()

        # Load existing results to append to
        if n_completed > 0 and os.path.exists(OUTPUT_CSV):
            existing_df = pd.read_csv(OUTPUT_CSV)
            all_results = existing_df.to_dict('records')

        for batch_idx in range(n_batches):
            batch_start = batch_idx * BATCH_SIZE
            batch_end = min(batch_start + BATCH_SIZE, n_remaining)
            batch = remaining_mutations[batch_start:batch_end]
            batch_num = batch_idx + 1

            # Progress update
            elapsed = time.time() - start_time
            mutations_done = batch_start + n_completed
            total_done_pct = 100 * mutations_done / total_mutations

            if batch_idx > 0 and elapsed > 0:
                rate = batch_start / elapsed  # mutations per second
                remaining_time = (n_remaining - batch_start) / rate if rate > 0 else 0
                eta_str = format_time(remaining_time)
            else:
                eta_str = "calculating..."

            progress_bar(mutations_done, total_mutations, prefix="Overall: ")
            print()
            print(f"   Batch {batch_num}/{n_batches} "
                  f"(mutations {batch_start+1}-{batch_end} of {n_remaining} remaining) "
                  f"| ETA: {eta_str}")

            # Print batch contents summary
            positions_in_batch = sorted(set(m['position'] for m in batch))
            if len(positions_in_batch) <= 5:
                pos_str = ", ".join(str(p) for p in positions_in_batch)
            else:
                pos_str = (f"{positions_in_batch[0]}-{positions_in_batch[-1]} "
                           f"({len(positions_in_batch)} positions)")
            print(f"   Positions: {pos_str}")

            # Run FoldX for this batch
            try:
                batch_results = run_foldx_batch(batch, WORK_DIR, pdb_name)

                if batch_results is not None:
                    all_results.extend(batch_results)
                    n_success += len(batch_results)
                    # Print summary of this batch
                    ddgs = [r['ddG'] for r in batch_results if not np.isnan(r.get('ddG', np.nan))]
                    if ddgs:
                        print(f"   Results: {len(ddgs)} ddG values, "
                              f"range [{min(ddgs):.2f}, {max(ddgs):.2f}] kcal/mol")
                else:
                    # Record failures
                    for mut in batch:
                        fail_result = mut.copy()
                        fail_result['ddG'] = np.nan
                        fail_result['effect_class'] = 'failed'
                        all_results.append(fail_result)
                    n_fail += len(batch)
                    print(f"   FAILED: {len(batch)} mutations in this batch")

            except Exception as e:
                print(f"   [EXCEPTION] Batch {batch_num} failed: {e}")
                for mut in batch:
                    fail_result = mut.copy()
                    fail_result['ddG'] = np.nan
                    fail_result['effect_class'] = 'failed'
                    all_results.append(fail_result)
                n_fail += len(batch)

            # Clean up intermediate files to save disk space
            cleanup_intermediate_files(WORK_DIR, pdb_stem)

            # Save intermediate results every 10 batches for safety
            if (batch_num % 10 == 0) or (batch_num == n_batches):
                interim_df = pd.DataFrame(all_results)
                interim_df.to_csv(OUTPUT_CSV, index=False)
                # Also save backup copy to E: drive if available
                backup_csv = r"E:/CALON_FH_BACKUP/03_COMPUTATIONAL_EXPENSIVE/analysis_outputs/foldx_saturation_mutagenesis.csv"
                try:
                    if os.path.exists("E:/CALON_FH_BACKUP"):
                        interim_df.to_csv(backup_csv, index=False)
                except Exception:
                    pass
                print(f"   [CHECKPOINT] Saved {len(interim_df)} results to {OUTPUT_CSV}")

        # Final progress
        total_elapsed = time.time() - start_time
        progress_bar(total_mutations, total_mutations, prefix="Overall: ")
        print()
        print(f"\n   FoldX runs complete!")
        print(f"   Total time: {format_time(total_elapsed)}")
        print(f"   Successful: {n_success}")
        print(f"   Failed: {n_fail}")
        if n_success > 0:
            print(f"   Rate: {total_elapsed / n_success:.1f}s per mutation")

    # ── 6. Save final results and rebuild atlas ───────────────────────────────
    print(f"\n[6/6] Saving final results and rebuilding atlas v2...")

    # Load the complete saturation mutagenesis results
    if os.path.exists(OUTPUT_CSV):
        sat_df = pd.read_csv(OUTPUT_CSV)
    else:
        print("   ERROR: No saturation mutagenesis results found!")
        sys.exit(1)

    # ── Summary statistics ────────────────────────────────────────────────────
    print(f"\n{'=' * 80}")
    print("SATURATION MUTAGENESIS RESULTS SUMMARY")
    print(f"{'=' * 80}")

    n_total = len(sat_df)
    n_valid = sat_df['ddG'].notna().sum()
    n_failed = sat_df['ddG'].isna().sum()

    print(f"   Total mutations:     {n_total}")
    print(f"   Valid ddG values:    {n_valid}")
    print(f"   Failed/missing:     {n_failed}")

    if n_valid > 0:
        print(f"\n   ddG statistics (kcal/mol):")
        print(f"     Mean:   {sat_df['ddG'].mean():.3f}")
        print(f"     Median: {sat_df['ddG'].median():.3f}")
        print(f"     Std:    {sat_df['ddG'].std():.3f}")
        print(f"     Min:    {sat_df['ddG'].min():.3f}")
        print(f"     Max:    {sat_df['ddG'].max():.3f}")

        # Effect classification
        print(f"\n   Effect classification:")
        effect_counts = sat_df['effect_class'].value_counts()
        for effect in ['stabilizing', 'neutral', 'mild', 'destabilizing', 'failed', 'unknown']:
            count = effect_counts.get(effect, 0)
            if count > 0:
                pct = 100 * count / n_total
                print(f"     {effect:<15s}: {count:6d} ({pct:5.1f}%)")

        # Per-position aggregates
        pos_agg = sat_df[sat_df['ddG'].notna()].groupby('position').agg(
            max_ddG=('ddG', 'max'),
            mean_ddG=('ddG', 'mean'),
            min_ddG=('ddG', 'min'),
            std_ddG=('ddG', 'std'),
            n_destabilizing=('effect_class', lambda x: (x == 'destabilizing').sum()),
            n_neutral=('effect_class', lambda x: (x == 'neutral').sum()),
            n_stabilizing=('effect_class', lambda x: (x == 'stabilizing').sum()),
        ).reset_index()

        # Find the worst mutation at each position
        worst_per_pos = sat_df[sat_df['ddG'].notna()].loc[
            sat_df[sat_df['ddG'].notna()].groupby('position')['ddG'].idxmax()
        ][['position', 'foldx_notation', 'ddG', 'mut_aa']].rename(
            columns={'foldx_notation': 'worst_mutation', 'ddG': 'worst_ddG', 'mut_aa': 'worst_mut_aa'}
        )
        pos_agg = pos_agg.merge(worst_per_pos[['position', 'worst_mutation', 'worst_mut_aa']], on='position', how='left')

        print(f"\n   Per-position statistics ({len(pos_agg)} positions):")
        print(f"     Positions with max ddG > 2.0 (destabilizing): "
              f"{(pos_agg['max_ddG'] > 2.0).sum()}")
        print(f"     Positions with max ddG > 5.0 (severely destabilizing): "
              f"{(pos_agg['max_ddG'] > 5.0).sum()}")
        print(f"     Positions where all mutations are neutral: "
              f"{(pos_agg['max_ddG'] <= 1.0).sum()}")

        # Top 20 most sensitive positions
        print(f"\n   Top 20 most mutationally sensitive positions:")
        print(f"   {'Pos':>4s} {'WT':>3s} {'Max ddG':>8s} {'Mean ddG':>9s} "
              f"{'#Destab':>8s} {'Worst':>10s}")
        print(f"   {'----':>4s} {'---':>3s} {'--------':>8s} {'---------':>9s} "
              f"{'--------':>8s} {'----------':>10s}")
        top20 = pos_agg.nlargest(20, 'max_ddG')
        for _, r in top20.iterrows():
            pos = int(r['position'])
            wt = sat_df[sat_df['position'] == pos].iloc[0]['wt_aa']
            print(f"   {pos:4d} {wt:>3s} {r['max_ddG']:8.2f} {r['mean_ddG']:9.2f} "
                  f"{int(r['n_destabilizing']):8d} {r['worst_mutation']:>10s}")

    # ── Rebuild Atlas v2 ──────────────────────────────────────────────────────
    print(f"\n{'=' * 80}")
    print("REBUILDING GENOTYPE-PHENOTYPE ATLAS v2")
    print(f"{'=' * 80}")

    if not os.path.exists(ATLAS_V1):
        print(f"   WARNING: Atlas v1 not found at {ATLAS_V1}")
        print(f"   Creating atlas v2 from scratch...")
        atlas = pd.DataFrame({'position': range(1, 861)})
    else:
        atlas = pd.read_csv(ATLAS_V1)
        print(f"   Loaded atlas v1: {len(atlas)} positions")

    # Compute position-level saturation mutagenesis metrics
    sat_valid = sat_df[sat_df['ddG'].notna()].copy()

    if len(sat_valid) > 0:
        pos_metrics = sat_valid.groupby('position').agg(
            sat_max_ddG=('ddG', 'max'),
            sat_mean_ddG=('ddG', 'mean'),
            sat_min_ddG=('ddG', 'min'),
            sat_median_ddG=('ddG', 'median'),
            sat_std_ddG=('ddG', 'std'),
            sat_n_destabilizing=('effect_class', lambda x: (x == 'destabilizing').sum()),
            sat_n_mild=('effect_class', lambda x: (x == 'mild').sum()),
            sat_n_neutral=('effect_class', lambda x: (x == 'neutral').sum()),
            sat_n_stabilizing=('effect_class', lambda x: (x == 'stabilizing').sum()),
            sat_n_computed=('ddG', 'count'),
        ).reset_index()

        # Find worst mutation per position
        worst_idx = sat_valid.groupby('position')['ddG'].idxmax()
        worst_muts = sat_valid.loc[worst_idx][['position', 'foldx_notation', 'mut_aa', 'ddG']].rename(
            columns={
                'foldx_notation': 'sat_worst_mutation',
                'mut_aa': 'sat_worst_mut_aa',
                'ddG': 'sat_worst_ddG',
            }
        )
        pos_metrics = pos_metrics.merge(worst_muts, on='position', how='left')

        # Compute mutational sensitivity score (fraction of destabilizing mutations)
        pos_metrics['mutational_sensitivity'] = (
            pos_metrics['sat_n_destabilizing'] / pos_metrics['sat_n_computed']
        ).round(4)

        # Merge into atlas
        atlas = atlas.merge(pos_metrics, on='position', how='left')

        # Update observed_max_ddG with saturation mutagenesis max ddG
        if 'observed_max_ddG' in atlas.columns:
            # Replace observed_max_ddG: use saturation data where available,
            # keep original where saturation data is missing
            atlas['observed_max_ddG_original'] = atlas['observed_max_ddG']
            atlas['observed_max_ddG'] = atlas['sat_max_ddG'].fillna(atlas['observed_max_ddG'])

        # Update risk_classification based on actual position-level data
        def classify_risk_v2(row):
            """Enhanced risk classification using saturation mutagenesis data."""
            at_pcsk9 = row.get('at_pcsk9_interface', 'No') == 'Yes'
            at_apob = row.get('at_apob_interface', 'No') == 'Yes'
            dom_pen = row.get('domain_penetrance_by_60_pct', 15.0)
            max_ddg = row.get('sat_max_ddG', np.nan)
            sensitivity = row.get('mutational_sensitivity', 0)

            # Position-specific classification using FoldX data
            if pd.notna(max_ddg):
                if max_ddg > 5.0 and (at_pcsk9 or at_apob or sensitivity > 0.5):
                    return "Critical"
                elif max_ddg > 5.0 or (max_ddg > 3.0 and dom_pen > 25):
                    return "High"
                elif max_ddg > 2.0 or (sensitivity > 0.3 and dom_pen > 15):
                    return "Moderate"
                elif max_ddg > 1.0:
                    return "Low-Moderate"
                else:
                    return "Low"
            else:
                # Fallback to domain-based classification
                if (at_pcsk9 or at_apob) and dom_pen > 25:
                    return "Critical"
                elif dom_pen > 20:
                    return "High"
                elif dom_pen >= 15:
                    return "Moderate"
                else:
                    return "Low"

        atlas['risk_classification_v1'] = atlas.get('risk_classification', 'Unknown')
        atlas['risk_classification'] = atlas.apply(classify_risk_v2, axis=1)

        print(f"\n   Position-level metrics added:")
        print(f"     Positions with saturation data: {pos_metrics['position'].nunique()}")
        print(f"     Mean max_ddG across positions: {pos_metrics['sat_max_ddG'].mean():.2f}")
        print(f"     Mean sensitivity score: {pos_metrics['mutational_sensitivity'].mean():.3f}")

        # Risk classification comparison
        print(f"\n   Risk classification comparison (v1 vs v2):")
        if 'risk_classification_v1' in atlas.columns:
            changed = atlas[atlas['risk_classification'] != atlas['risk_classification_v1']]
            print(f"     Positions with changed classification: {len(changed)}")

            # New distribution
            print(f"\n   Atlas v2 risk classification:")
            risk_counts = atlas['risk_classification'].value_counts()
            for risk_level in sorted(risk_counts.index):
                count = risk_counts[risk_level]
                pct = 100 * count / len(atlas)
                print(f"     {risk_level:<15s}: {count:4d} positions ({pct:5.1f}%)")

            # Show upgrades/downgrades
            if len(changed) > 0:
                print(f"\n   Classification changes (first 20):")
                for _, r in changed.head(20).iterrows():
                    print(f"     Pos {int(r['position']):4d}: "
                          f"{r['risk_classification_v1']} -> {r['risk_classification']} "
                          f"(max_ddG={r.get('sat_max_ddG', 'N/A')})")

    # Save atlas v2
    atlas.to_csv(ATLAS_V2, index=False)
    print(f"\n   Atlas v2 saved to: {ATLAS_V2}")
    print(f"   Total columns: {len(atlas.columns)}")
    print(f"   Total rows: {len(atlas)}")

    # ── Final summary ─────────────────────────────────────────────────────────
    print(f"\n{'=' * 80}")
    print("SATURATION MUTAGENESIS COMPLETE")
    print(f"{'=' * 80}")
    print(f"   Saturation mutagenesis CSV: {OUTPUT_CSV}")
    print(f"   Atlas v2 CSV:               {ATLAS_V2}")
    print(f"   Total mutations computed:   {n_valid}")
    if n_valid > 0:
        n_destab = (sat_df['effect_class'] == 'destabilizing').sum()
        n_neutral_all = (sat_df['effect_class'] == 'neutral').sum()
        print(f"   Destabilizing mutations:    {n_destab} ({100*n_destab/n_valid:.1f}%)")
        print(f"   Neutral mutations:          {n_neutral_all} ({100*n_neutral_all/n_valid:.1f}%)")
    print(f"{'=' * 80}")


if __name__ == "__main__":
    main()

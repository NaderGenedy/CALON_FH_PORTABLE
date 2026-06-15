#!/usr/bin/env python3
"""
64_run_lpa_saturation_parallel.py
==================================
Parallel FoldX saturation mutagenesis for LPA.
Splits mutation lists into chunks, runs N parallel FoldX processes.

8,037 mutations across KIV-7+8, KIV-10, and Protease domains.
With 8 parallel jobs on 16 cores: ~3 hours total.

Author: Dr Nader Genedy
Date:   April 2026
"""

import os
import subprocess
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import pandas as pd

FOLDX_BIN = r"D:/foldx_backup/foldx5_Windows_2/foldx_1_20270131.exe"
FOLDX_DIR = r"C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/foldx/lpa"
N_PARALLEL = 8  # Number of parallel FoldX jobs
N_RUNS = 3      # FoldX runs per mutation (for averaging)

def split_mutation_list(list_file, n_chunks):
    """Split a mutation list into N chunks."""
    with open(list_file) as f:
        mutations = [l.strip() for l in f if l.strip()]

    chunk_size = max(1, len(mutations) // n_chunks)
    chunks = []
    for i in range(0, len(mutations), chunk_size):
        chunk = mutations[i:i + chunk_size]
        chunk_file = list_file.replace('.txt', f'_chunk{len(chunks)}.txt')
        with open(chunk_file, 'w') as f:
            for m in chunk:
                f.write(m + '\n')
        chunks.append((chunk_file, len(chunk)))

    return chunks


def run_foldx_chunk(args):
    """Run FoldX BuildModel on a single chunk."""
    pdb_file, chunk_file, work_dir, chunk_id = args
    chunk_name = os.path.basename(chunk_file)

    # Create a subdirectory for this chunk to avoid file conflicts
    chunk_dir = os.path.join(work_dir, f"chunk_{chunk_id}")
    os.makedirs(chunk_dir, exist_ok=True)

    # Copy required files
    import shutil
    shutil.copy(os.path.join(work_dir, pdb_file), chunk_dir)
    shutil.copy(os.path.join(work_dir, 'rotabase.txt'), chunk_dir)
    shutil.copy(chunk_file, os.path.join(chunk_dir, 'individual_list.txt'))

    # Run FoldX
    cmd = [
        FOLDX_BIN,
        '--command', 'BuildModel',
        '--pdb', pdb_file,
        '--mutant-file', 'individual_list.txt',
        '--numberOfRuns', str(N_RUNS),
        '--out-pdb', 'false',  # Don't save mutant PDBs (saves disk space)
    ]

    start = time.time()
    try:
        result = subprocess.run(cmd, cwd=chunk_dir, capture_output=True, text=True, timeout=7200)
        elapsed = time.time() - start

        # Check for output
        dif_files = [f for f in os.listdir(chunk_dir) if f.startswith('Dif_')]
        avg_files = [f for f in os.listdir(chunk_dir) if f.startswith('Average_')]

        n_muts = sum(1 for l in open(os.path.join(chunk_dir, 'individual_list.txt')) if l.strip())
        return {
            'chunk_id': chunk_id,
            'chunk_file': chunk_name,
            'n_mutations': n_muts,
            'elapsed': elapsed,
            'success': len(dif_files) > 0,
            'chunk_dir': chunk_dir,
            'dif_files': dif_files,
            'avg_files': avg_files,
        }
    except subprocess.TimeoutExpired:
        return {'chunk_id': chunk_id, 'success': False, 'elapsed': 7200, 'error': 'timeout'}
    except Exception as e:
        return {'chunk_id': chunk_id, 'success': False, 'elapsed': 0, 'error': str(e)}


def collect_results(chunk_dirs, output_prefix):
    """Collect and merge results from all chunks."""
    all_dif = []
    all_avg = []

    for chunk_dir in chunk_dirs:
        for f in os.listdir(chunk_dir):
            if f.startswith('Dif_'):
                dif_path = os.path.join(chunk_dir, f)
                with open(dif_path) as fh:
                    lines = fh.readlines()
                    # Skip header lines (first 8)
                    data_lines = [l for l in lines if not l.startswith('FoldX') and
                                  not l.startswith('by ') and not l.startswith('---') and
                                  not l.startswith('Jesper') and not l.startswith('Luis') and
                                  not l.startswith('Peter') and not l.startswith('Lies') and
                                  not l.startswith('\n') and not l.startswith('PDB') and
                                  not l.startswith('Output')]
                    all_dif.extend(data_lines)

            if f.startswith('Average_'):
                avg_path = os.path.join(chunk_dir, f)
                with open(avg_path) as fh:
                    lines = fh.readlines()
                    data_lines = [l for l in lines if not l.startswith('FoldX') and
                                  not l.startswith('by ') and not l.startswith('---') and
                                  not l.startswith('Jesper') and not l.startswith('Luis') and
                                  not l.startswith('Peter') and not l.startswith('Lies') and
                                  not l.startswith('\n') and not l.startswith('PDB') and
                                  not l.startswith('Output')]
                    all_avg.extend(data_lines)

    # Write combined results
    with open(f"{output_prefix}_dif_combined.txt", 'w') as f:
        f.writelines(all_dif)

    with open(f"{output_prefix}_avg_combined.txt", 'w') as f:
        f.writelines(all_avg)

    return len(all_dif), len(all_avg)


if __name__ == '__main__':
    print("=" * 70)
    print(f"LPA SATURATION MUTAGENESIS — {N_PARALLEL} PARALLEL FoldX JOBS")
    print("=" * 70)

    jobs_to_run = [
        ('LPA_KIV78.pdb', 'sat_KIV78_list.txt', 'KIV-7+8 lysine binding'),
        ('LPA_KIV10_Protease.pdb', 'sat_KIV10_list.txt', 'KIV-10 assembly'),
        ('LPA_KIV10_Protease.pdb', 'sat_Protease_list.txt', 'Protease domain'),
    ]

    for pdb, mut_list, label in jobs_to_run:
        list_path = os.path.join(FOLDX_DIR, mut_list)
        n_muts = sum(1 for l in open(list_path) if l.strip())

        if n_muts == 0:
            print(f"\n  {label}: 0 mutations, skipping")
            continue

        print(f"\n{'='*70}")
        print(f"  {label}: {n_muts} mutations")
        print(f"  PDB: {pdb}")
        print(f"  Splitting into {N_PARALLEL} chunks...")

        # Split
        chunks = split_mutation_list(list_path, N_PARALLEL)
        print(f"  Created {len(chunks)} chunks: {[c[1] for c in chunks]}")

        # Prepare job arguments
        job_args = [
            (pdb, chunk_file, FOLDX_DIR, i)
            for i, (chunk_file, _) in enumerate(chunks)
        ]

        # Run in parallel
        print(f"  Starting {len(job_args)} parallel FoldX processes...")
        start_time = time.time()

        with ProcessPoolExecutor(max_workers=N_PARALLEL) as executor:
            futures = {executor.submit(run_foldx_chunk, args): args for args in job_args}

            for future in as_completed(futures):
                result = future.result()
                if result.get('success'):
                    print(f"    Chunk {result['chunk_id']}: {result['n_mutations']} mutations, "
                          f"{result['elapsed']:.0f}s ({result['n_mutations']/max(result['elapsed'],1):.1f} mut/s)")
                else:
                    print(f"    Chunk {result['chunk_id']}: FAILED - {result.get('error', 'unknown')}")

        total_time = time.time() - start_time
        print(f"  Total time: {total_time/60:.1f} minutes")

        # Collect results
        chunk_dirs = [os.path.join(FOLDX_DIR, f"chunk_{i}") for i in range(len(chunks))]
        output_prefix = os.path.join(FOLDX_DIR, label.replace(' ', '_').replace('+', '').replace('-', ''))
        n_dif, n_avg = collect_results(chunk_dirs, output_prefix)
        print(f"  Collected: {n_dif} Dif lines, {n_avg} Average lines")

    print("\n" + "=" * 70)
    print("ALL SATURATION MUTAGENESIS COMPLETE")
    print("=" * 70)

"""
TUDOR Locked Rerun Orchestrator
================================
Runs phases 2B, 2C, 2D, 3, 4, 5, 6 in sequence.
Each phase produces a checkpoint CSV; the next phase skips if its checkpoint exists.

Usage:
    python TUDOR_LOCKED_RERUN_orchestrator.py
    python TUDOR_LOCKED_RERUN_orchestrator.py --resume  # skip completed phases
    python TUDOR_LOCKED_RERUN_orchestrator.py --skip 2B 2C  # skip specific phases
"""
import os, sys, time, subprocess, argparse

ROOT = r'C:/Users/nader/Downloads/calon_ukb_pipeline'
PYTHON = r'C:/Users/nader/AppData/Local/Programs/Python/Python312/python.exe'

PHASES = [
    ('2B', 'TUDOR_LOCKED_RERUN_phase2B_mpr.py',         r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_mpr_adherence.csv'),
    ('2C', 'TUDOR_LOCKED_RERUN_phase2C_ldl_ut.py',      r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_ldl_ut_v2.csv'),
    ('2D', 'TUDOR_LOCKED_RERUN_phase2D_cohorts.py',     r'D:/Projects/CALON_AlphaFold_Rebuild/data/ukb_cohort_primary_v2.csv'),
    ('3',  'TUDOR_LOCKED_RERUN_phase3_wales_dedup.py',  r'D:/Projects/CALON_AlphaFold_Rebuild/data/wales_dedup_v2.csv'),
    ('4',  'TUDOR_LOCKED_RERUN_phase4_refit_enet.py',   r'D:/Projects/CALON_AlphaFold_Rebuild/data/loco_predictions_v2.csv'),
    ('5',  'TUDOR_LOCKED_RERUN_phase5_stats.py',        r'D:/Projects/CALON_AlphaFold_Rebuild/data/tudor_statistics_summary_v2.csv'),
    ('6',  'TUDOR_LOCKED_RERUN_phase6_regenerate.py',   r'C:/Users/nader/Downloads/calon_ukb_pipeline/TUDOR_Manuscript_v6_TRACEABLE.docx'),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--skip', nargs='+', default=[])
    args = ap.parse_args()

    env = os.environ.copy()
    env['PYTHONIOENCODING'] = 'utf-8'

    t_start = time.time()
    for name, script, checkpoint in PHASES:
        if name in args.skip:
            print(f'\n=== Phase {name}: SKIPPED (--skip) ===')
            continue
        if args.resume and os.path.exists(checkpoint):
            print(f'\n=== Phase {name}: SKIPPED (checkpoint exists: {checkpoint}) ===')
            continue
        print(f'\n{"="*70}')
        print(f'=== Phase {name}: {script} ===')
        print(f'{"="*70}')
        t0 = time.time()
        result = subprocess.run([PYTHON, os.path.join(ROOT, script)],
                                env=env, cwd=ROOT)
        if result.returncode != 0:
            print(f'\nPhase {name} FAILED with exit code {result.returncode}')
            print('Stopping orchestrator.')
            sys.exit(result.returncode)
        if not os.path.exists(checkpoint):
            print(f'\nWARNING: Phase {name} completed but checkpoint {checkpoint} not present')
        else:
            sz = os.path.getsize(checkpoint) / 1e6
            print(f'\nPhase {name} OK in {time.time()-t0:.0f}s; checkpoint: {checkpoint} ({sz:.1f} MB)')
    print(f'\n{"="*70}')
    print(f'All phases complete in {(time.time()-t_start)/60:.1f} min')
    print(f'{"="*70}')


if __name__ == '__main__':
    main()

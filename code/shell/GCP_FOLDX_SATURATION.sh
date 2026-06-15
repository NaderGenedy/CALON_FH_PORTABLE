#!/bin/bash
# ==============================================================================
# GCP FoldX SATURATION MUTAGENESIS FOR LPA
# ==============================================================================
# INSTRUCTIONS:
#   1. Create a GCP VM: e2-standard-16 (16 vCPU, 64 GB RAM)
#      Region: europe-west2 (London) or us-central1
#      OS: Ubuntu 22.04
#      Disk: 50 GB SSD
#
#   2. Upload this script + required files to the VM:
#      gcloud compute scp GCP_FOLDX_SATURATION.sh <VM_NAME>:~/
#      gcloud compute scp foldx_bundle.tar.gz <VM_NAME>:~/
#
#   3. SSH into VM and run:
#      chmod +x GCP_FOLDX_SATURATION.sh
#      nohup bash GCP_FOLDX_SATURATION.sh > foldx_run.log 2>&1 &
#
#   4. Monitor: tail -f foldx_run.log
#
#   5. Download results when done:
#      gcloud compute scp <VM_NAME>:~/results/*.csv ./
#
# ESTIMATED TIME: ~45 minutes with 16 cores
# ESTIMATED COST: ~$0.50 (e2-standard-16 at $0.67/hr)
# ==============================================================================

set -e

echo "=============================================="
echo "LPA FoldX SATURATION MUTAGENESIS (GCP)"
echo "Date: $(date)"
echo "Cores: $(nproc)"
echo "=============================================="

# --- Setup ---
WORK_DIR=~/foldx_lpa
mkdir -p $WORK_DIR/results
cd $WORK_DIR

# Unpack FoldX bundle (uploaded separately due to licence)
if [ -f ~/foldx_bundle.tar.gz ]; then
    tar xzf ~/foldx_bundle.tar.gz -C $WORK_DIR/
    echo "FoldX bundle unpacked"
else
    echo "ERROR: Upload foldx_bundle.tar.gz first"
    echo "Bundle must contain: foldx binary, rotabase.txt, PDB files, mutation lists"
    exit 1
fi

FOLDX=$WORK_DIR/foldx
chmod +x $FOLDX 2>/dev/null

# --- Verify files ---
echo ""
echo "=== Checking files ==="
for f in rotabase.txt LPA_KIV78.pdb LPA_KIV10_Protease.pdb sat_KIV78_list.txt sat_KIV10_list.txt sat_Protease_list.txt; do
    if [ -f "$f" ]; then
        echo "  OK: $f ($(wc -l < $f) lines)"
    else
        echo "  MISSING: $f"
    fi
done

# --- Function: split and run parallel ---
run_saturation() {
    local PDB=$1
    local MUT_LIST=$2
    local LABEL=$3
    local N_CORES=$(nproc)
    local N_JOBS=$((N_CORES - 2))  # Leave 2 cores for OS

    local N_MUTS=$(wc -l < $MUT_LIST)
    echo ""
    echo "=============================================="
    echo "  $LABEL: $N_MUTS mutations, $N_JOBS parallel jobs"
    echo "=============================================="

    # Split mutation list into chunks
    local CHUNK_SIZE=$(( (N_MUTS + N_JOBS - 1) / N_JOBS ))
    split -l $CHUNK_SIZE -d -a 2 $MUT_LIST chunk_${LABEL}_

    # Run each chunk in parallel
    local PIDS=()
    for chunk in chunk_${LABEL}_*; do
        local CHUNK_DIR="${WORK_DIR}/run_${chunk}"
        mkdir -p $CHUNK_DIR
        cp $PDB $CHUNK_DIR/
        cp rotabase.txt $CHUNK_DIR/
        cp $chunk $CHUNK_DIR/individual_list.txt

        (
            cd $CHUNK_DIR
            $FOLDX --command BuildModel \
                --pdb $PDB \
                --mutant-file individual_list.txt \
                --numberOfRuns 3 \
                --out-pdb false \
                > foldx.log 2>&1
            echo "  Chunk $chunk done ($(wc -l < individual_list.txt) mutations)"
        ) &
        PIDS+=($!)
    done

    # Wait for all
    echo "  Waiting for ${#PIDS[@]} jobs..."
    local START=$(date +%s)
    for pid in "${PIDS[@]}"; do
        wait $pid
    done
    local END=$(date +%s)
    local ELAPSED=$(( END - START ))
    echo "  All chunks done in ${ELAPSED}s ($(( ELAPSED / 60 )) min)"

    # Collect results
    echo "  Collecting results..."
    local OUT_FILE="$WORK_DIR/results/${LABEL}_saturation_ddg.csv"
    echo "mutation,ddG,sd" > $OUT_FILE

    for chunk_dir in ${WORK_DIR}/run_chunk_${LABEL}_*; do
        if [ -f "$chunk_dir/Average_${PDB%.pdb}.fxout" ]; then
            # Parse Average file: extract mutation name and total energy (ddG)
            awk 'NR>8 && !/^$/ && !/^FoldX/ && !/^by/ && !/^---/ && !/^Jesper/ && !/^Luis/ && !/^PDB/ && !/^Output/ {
                split($1, a, "_");
                mut_num = a[length(a)];
                printf "%s,%.4f,%.4f\n", $1, $3, $2
            }' "$chunk_dir/Average_${PDB%.pdb}.fxout" >> $OUT_FILE
        fi
    done

    local N_RESULTS=$(( $(wc -l < $OUT_FILE) - 1 ))
    echo "  Results: $N_RESULTS mutations in $OUT_FILE"

    # Cleanup chunk files
    rm -f chunk_${LABEL}_*
}

# --- Run all three domains ---
START_ALL=$(date +%s)

run_saturation "LPA_KIV78.pdb" "sat_KIV78_list.txt" "KIV78"
run_saturation "LPA_KIV10_Protease.pdb" "sat_KIV10_list.txt" "KIV10"
run_saturation "LPA_KIV10_Protease.pdb" "sat_Protease_list.txt" "Protease"

END_ALL=$(date +%s)
TOTAL=$(( END_ALL - START_ALL ))

echo ""
echo "=============================================="
echo "ALL SATURATION MUTAGENESIS COMPLETE"
echo "Total time: ${TOTAL}s ($(( TOTAL / 60 )) min)"
echo "=============================================="
echo ""
echo "Results in $WORK_DIR/results/"
ls -lh $WORK_DIR/results/
echo ""
echo "Download with:"
echo "  gcloud compute scp $(hostname):~/foldx_lpa/results/*.csv ./"
echo ""
echo "NEXT: Run locally:"
echo "  python 65_parse_lpa_saturation.py"

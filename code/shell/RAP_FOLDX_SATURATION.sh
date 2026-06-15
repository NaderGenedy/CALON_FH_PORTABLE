#!/bin/bash
# ==============================================================================
# UKB-RAP: FoldX SATURATION MUTAGENESIS FOR LPA
# ==============================================================================
# RUN ON: UKB-RAP JupyterLab terminal (already have compute access)
#
# STEPS:
#   1. Upload foldx_bundle.tar.gz to RAP: dx upload foldx_bundle.tar.gz
#   2. Download to workspace: dx download foldx_bundle.tar.gz
#   3. Run this script: bash RAP_FOLDX_SATURATION.sh
#
# The bundle must contain:
#   - foldx binary (Linux version: foldx_20261231)
#   - rotabase.txt
#   - LPA_KIV78.pdb
#   - LPA_KIV10_Protease.pdb
#   - sat_KIV78_list.txt (3,040 mutations)
#   - sat_KIV10_list.txt (1,520 mutations)
#   - sat_Protease_list.txt (4,636 mutations)
#
# ESTIMATED: ~1 hour on RAP (8 cores typical)
# ==============================================================================

set -e

echo "=============================================="
echo "LPA SATURATION MUTAGENESIS (UKB-RAP)"
echo "$(date)"
echo "=============================================="

WORK=$HOME/foldx_lpa
mkdir -p $WORK/results
cd $WORK

# Unpack
if [ -f ~/foldx_bundle.tar.gz ]; then
    tar xzf ~/foldx_bundle.tar.gz
    echo "Bundle unpacked"
elif [ -f foldx_bundle.tar.gz ]; then
    tar xzf foldx_bundle.tar.gz
    echo "Bundle unpacked"
else
    echo "ERROR: foldx_bundle.tar.gz not found"
    echo "Upload with: dx upload foldx_bundle.tar.gz"
    exit 1
fi

# Find FoldX binary
FOLDX=$(find . -name "foldx*" -type f -executable | head -1)
if [ -z "$FOLDX" ]; then
    FOLDX=$(find . -name "foldx*" -type f | head -1)
    chmod +x "$FOLDX"
fi
echo "FoldX: $FOLDX"

N_CORES=$(nproc)
N_JOBS=$((N_CORES > 2 ? N_CORES - 2 : 1))
echo "Cores: $N_CORES, Jobs: $N_JOBS"

# Function: run one domain's saturation
run_domain() {
    PDB=$1; LIST=$2; LABEL=$3
    N=$(wc -l < $LIST)
    echo ""
    echo "--- $LABEL: $N mutations ---"

    CHUNK_SIZE=$(( (N + N_JOBS - 1) / N_JOBS ))
    split -l $CHUNK_SIZE -d $LIST chunk_

    PIDS=()
    for c in chunk_*; do
        DIR=$WORK/run_${LABEL}_${c}
        mkdir -p $DIR
        cp $PDB rotabase.txt $DIR/
        cp $c $DIR/individual_list.txt

        (cd $DIR && "$WORK/$FOLDX" --command BuildModel --pdb $PDB \
            --mutant-file individual_list.txt --numberOfRuns 3 \
            --out-pdb false > log.txt 2>&1 && echo "  $c done") &
        PIDS+=($!)
    done

    for p in "${PIDS[@]}"; do wait $p; done

    # Collect
    cat $WORK/run_${LABEL}_*/Average_*.fxout 2>/dev/null | \
        grep -v "^FoldX\|^by \|^---\|^Jesper\|^Luis\|^PDB\|^Output\|^$" \
        > $WORK/results/${LABEL}_avg.txt

    N_RES=$(wc -l < $WORK/results/${LABEL}_avg.txt)
    echo "  $LABEL: $N_RES results collected"

    rm -f chunk_*
}

START=$(date +%s)
run_domain "LPA_KIV78.pdb" "sat_KIV78_list.txt" "KIV78"
run_domain "LPA_KIV10_Protease.pdb" "sat_KIV10_list.txt" "KIV10"
run_domain "LPA_KIV10_Protease.pdb" "sat_Protease_list.txt" "Protease"
END=$(date +%s)

echo ""
echo "=============================================="
echo "DONE in $(( (END-START)/60 )) minutes"
ls -lh $WORK/results/
echo ""
echo "Upload results: dx upload $WORK/results/*.txt --destination /"
echo "=============================================="

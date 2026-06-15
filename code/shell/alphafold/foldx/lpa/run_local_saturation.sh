#!/bin/bash
FOLDX_BIN="D:/foldx_backup/foldx5_Windows_2/foldx_1_20270131.exe"
cd "C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/foldx/lpa"

> results_KIV78_local.txt
COUNT=0
TOTAL=$(wc -l < sat_KIV78_list.txt)

while IFS= read -r line; do
  echo "$line" > individual_list.txt
  "$FOLDX_BIN" --command BuildModel --pdb LPA_KIV78.pdb \
    --mutant-file individual_list.txt --numberOfRuns 1 \
    --out-pdb false > /dev/null 2>&1
  
  if [ -f Average_LPA_KIV78.fxout ]; then
    tail -1 Average_LPA_KIV78.fxout >> results_KIV78_local.txt
    rm -f Average_* Dif_* Raw_* PdbList_* LPA_KIV78_*.pdb 2>/dev/null
  fi
  
  COUNT=$((COUNT + 1))
  if [ $((COUNT % 50)) -eq 0 ]; then
    echo "Progress: $COUNT / $TOTAL $(date)"
  fi
done < sat_KIV78_list.txt

echo "COMPLETE: $(wc -l < results_KIV78_local.txt) results $(date)"

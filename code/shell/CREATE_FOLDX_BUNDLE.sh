#!/bin/bash
# ==============================================================================
# CREATE FoldX BUNDLE FOR UKB-RAP UPLOAD
# ==============================================================================
# After downloading Linux FoldX binary, run this to create the upload bundle.
#
# USAGE:
#   1. Download Linux FoldX from foldx.crg.eu
#   2. Place the binary in this directory (e.g., foldx_20261231)
#   3. Run: bash CREATE_FOLDX_BUNDLE.sh
#   4. Upload to RAP: dx upload foldx_bundle.tar.gz --destination /
#   5. On RAP terminal: dx download foldx_bundle.tar.gz && bash RAP_FOLDX_SATURATION.sh
# ==============================================================================

FOLDX_DIR="C:/Users/nader/Downloads/calon_ukb_pipeline/alphafold/foldx/lpa"
BUNDLE_DIR="D:/foldx_rap_bundle"

mkdir -p "$BUNDLE_DIR"

# Copy PDB structures
cp "$FOLDX_DIR/LPA_KIV78.pdb" "$BUNDLE_DIR/"
cp "$FOLDX_DIR/LPA_KIV10_Protease.pdb" "$BUNDLE_DIR/"

# Copy rotabase
cp "$FOLDX_DIR/rotabase.txt" "$BUNDLE_DIR/"

# Copy mutation lists
cp "$FOLDX_DIR/sat_KIV78_list.txt" "$BUNDLE_DIR/"
cp "$FOLDX_DIR/sat_KIV10_list.txt" "$BUNDLE_DIR/"
cp "$FOLDX_DIR/sat_Protease_list.txt" "$BUNDLE_DIR/"

# Copy the RAP script
cp "C:/Users/nader/Downloads/calon_ukb_pipeline/RAP_FOLDX_SATURATION.sh" "$BUNDLE_DIR/"

echo "=== Bundle contents ==="
ls -lh "$BUNDLE_DIR/"

echo ""
echo "NEXT: Copy your Linux FoldX binary into $BUNDLE_DIR/"
echo "Then: cd $BUNDLE_DIR && tar czf ../foldx_bundle.tar.gz *"
echo "Then: dx upload D:/foldx_bundle.tar.gz --destination /"
echo ""
echo "Files needed in bundle:"
echo "  foldx (Linux binary)     <- DOWNLOAD FROM foldx.crg.eu"
echo "  rotabase.txt             <- included"
echo "  LPA_KIV78.pdb            <- included"
echo "  LPA_KIV10_Protease.pdb   <- included"
echo "  sat_KIV78_list.txt       <- 3,040 mutations"
echo "  sat_KIV10_list.txt       <- 1,520 mutations"
echo "  sat_Protease_list.txt    <- 4,636 mutations"
echo "  RAP_FOLDX_SATURATION.sh  <- included"

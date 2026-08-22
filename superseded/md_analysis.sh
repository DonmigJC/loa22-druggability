#!/bin/bash
# Complete MD analysis for all three leads
# Run from: bash ~/Research/md_analysis.sh

BASE=~/Research/04_Molecular_Dynamics
LEADS=("Lead_1" "Lead_4" "Lead_5")

for LEAD in "${LEADS[@]}"; do
    echo "========================================"
    echo "Analyzing $LEAD"
    echo "========================================"
    
    MD_DIR=$BASE/$LEAD/md_run
    OUT_DIR=$BASE/$LEAD/analysis
    mkdir -p $OUT_DIR
    cd $MD_DIR

    # --- RMSD: Protein backbone ---
    echo "Computing protein backbone RMSD..."
    echo -e "Backbone\nBackbone" | gmx rms \
        -s md.tpr \
        -f md.xtc \
        -o $OUT_DIR/rmsd_backbone.xvg \
        -tu ns 2>/dev/null
    
    # --- RMSD: Ligand ---
    echo "Computing ligand RMSD..."
    echo -e "UNL\nUNL" | gmx rms \
        -s md.tpr \
        -f md.xtc \
        -o $OUT_DIR/rmsd_ligand.xvg \
        -tu ns 2>/dev/null

    # --- RMSF: Per-residue flexibility ---
    echo "Computing RMSF..."
    echo "Protein" | gmx rmsf \
        -s md.tpr \
        -f md.xtc \
        -o $OUT_DIR/rmsf.xvg \
        -res 2>/dev/null

    # --- Radius of gyration ---
    echo "Computing radius of gyration..."
    echo "Protein" | gmx gyrate \
        -s md.tpr \
        -f md.xtc \
        -o $OUT_DIR/gyrate.xvg 2>/dev/null

    # --- Hydrogen bonds between ligand and protein ---
    echo "Computing hydrogen bonds..."
    echo -e "Protein\nUNL" | gmx hbond \
        -s md.tpr \
        -f md.xtc \
        -num $OUT_DIR/hbnum.xvg 2>/dev/null

    echo "$LEAD analysis complete. Files in $OUT_DIR"
    echo ""
done

echo "All analysis complete."

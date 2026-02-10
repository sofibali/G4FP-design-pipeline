#!/bin/bash
# Regenerate all LigandMPNN outputs with corrected parameters
# This will DELETE existing outputs and regenerate with proper sequence diversity

set -e

echo "=========================================================================="
echo "REGENERATE ALL LIGANDMPNN OUTPUTS"
echo "=========================================================================="
echo ""
echo "This script will:"
echo "  1. Delete all existing LigandMPNN outputs (01_ligandmpnn directories)"
echo "  2. Regenerate sequences with FIXED parameters"
echo "  3. Generate ~10,000 DIVERSE sequences per structure"
echo ""
echo "⚠️  WARNING: This will DELETE existing outputs for all 10 structures!"
echo ""
read -p "Continue? (yes/no): " -r
echo ""

if [[ ! $REPLY =~ ^[Yy][Ee][Ss]$ ]]; then
    echo "Aborted."
    exit 0
fi

# Count output directories
OUTPUT_DIRS=(output_G4FP_*)
TOTAL=${#OUTPUT_DIRS[@]}

echo "Found $TOTAL output directories"
echo ""

# Delete existing LigandMPNN outputs
echo "Step 1: Deleting existing LigandMPNN outputs..."
DELETED=0
for dir in "${OUTPUT_DIRS[@]}"; do
    if [ -d "$dir/01_ligandmpnn" ]; then
        echo "  Removing $dir/01_ligandmpnn/"
        rm -rf "$dir/01_ligandmpnn"
        ((DELETED++))
    fi
done
echo "  ✓ Deleted $DELETED directories"
echo ""

# Run batch script with force rerun
echo "Step 2: Running batch LigandMPNN generation..."
echo "  This will take 3-5 hours for all 10 structures"
echo "  Progress will be logged to: batch_ligandmpnn_progress.log"
echo "  Detailed output in: batch_ligandmpnn_detailed.log"
echo ""

# Run in background with nohup
nohup bash batch_run_ligandmpnn.sh --force-rerun > regenerate_ligandmpnn.log 2>&1 &
PID=$!

echo "✓ Batch job started (PID: $PID)"
echo ""
echo "Monitor progress with:"
echo "  tail -f batch_ligandmpnn_progress.log"
echo ""
echo "Or check this output:"
echo "  tail -f regenerate_ligandmpnn.log"
echo ""
echo "To check if still running:"
echo "  ps -p $PID"
echo ""

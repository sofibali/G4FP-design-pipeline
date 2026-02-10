#!/bin/bash
# Batch run LigandMPNN for multiple PDB files
# Usage: nohup ./batch_run_ligandmpnn.sh [--force-rerun] &
#
# Options:
#   --force-rerun    Delete existing outputs and regenerate all sequences

# DO NOT exit on error - we want to continue processing other PDFs even if one fails
set +e

# Activate LigandMPNN environment
eval "$(conda shell.bash hook)"
conda activate unified_mpnn

# ============================================================================
# Configuration
# ============================================================================
INPUTS_DIR="/home/sbali/LigandMPNN/designG4FP/inputs"
CONFIG_TEMPLATE="pipeline_config_template.yaml"
CONFIG_FILE="pipeline_config.yaml"
DETAILED_LOG="batch_ligandmpnn_detailed.log"
PROGRESS_LOG="batch_ligandmpnn_progress.log"

# Check for force rerun flag
FORCE_RERUN=0
if [[ "$1" == "--force-rerun" ]]; then
    FORCE_RERUN=1
    echo "⚠️  FORCE RERUN MODE: Will delete and regenerate all outputs"
    echo ""
fi

# ============================================================================
# Create template config backup if it doesn't exist
# ============================================================================
if [ ! -f "$CONFIG_TEMPLATE" ]; then
    echo "Creating template config backup..."
    cp "$CONFIG_FILE" "$CONFIG_TEMPLATE"
fi

# ============================================================================
# Initialize logs
# ============================================================================
echo "============================================================================" > "$PROGRESS_LOG"
echo "Batch LigandMPNN Run - Started at $(date)" >> "$PROGRESS_LOG"
echo "============================================================================" >> "$PROGRESS_LOG"
echo "" >> "$PROGRESS_LOG"

echo "============================================================================" > "$DETAILED_LOG"
echo "Batch LigandMPNN Detailed Output - Started at $(date)" >> "$DETAILED_LOG"
echo "============================================================================" >> "$DETAILED_LOG"
echo "" >> "$DETAILED_LOG"

# ============================================================================
# Find all PDB files
# ============================================================================
PDB_FILES=($(ls "$INPUTS_DIR"/*.pdb 2>/dev/null))
TOTAL_PDBS=${#PDB_FILES[@]}

if [ $TOTAL_PDBS -eq 0 ]; then
    echo "ERROR: No PDB files found in $INPUTS_DIR" | tee -a "$PROGRESS_LOG"
    exit 1
fi

echo "Found $TOTAL_PDBS PDB files to process:" | tee -a "$PROGRESS_LOG"
for pdb in "${PDB_FILES[@]}"; do
    echo "  - $(basename $pdb)" | tee -a "$PROGRESS_LOG"
done
echo "" | tee -a "$PROGRESS_LOG"

# ============================================================================
# Function to check if LigandMPNN output is complete
# ============================================================================
check_output_complete() {
    local output_dir=$1
    local seqs_dir="${output_dir}/01_ligandmpnn/seqs"
    
    # Check if directory exists
    if [ ! -d "$seqs_dir" ]; then
        return 1  # Not complete
    fi
    
    # Check if FASTA files exist
    local fasta_files=($(ls "$seqs_dir"/*.fa 2>/dev/null))
    if [ ${#fasta_files[@]} -eq 0 ]; then
        return 1  # No FASTA files
    fi
    
    # Check if FASTA files have content (at least one sequence)
    local has_sequences=0
    for fasta in "${fasta_files[@]}"; do
        if grep -q "^>" "$fasta" 2>/dev/null; then
            has_sequences=1
            break
        fi
    done
    
    if [ $has_sequences -eq 0 ]; then
        return 1  # FASTA files exist but are empty
    fi
    
    # Check sequence count
    local seq_count=$(grep -c "^>" "$seqs_dir"/*.fa 2>/dev/null || echo "0")
    if [ "$seq_count" -lt 100 ]; then
        echo "    ⚠️  Warning: Only $seq_count sequences found (expected ~10000)" | tee -a "$PROGRESS_LOG"
        return 1  # Incomplete run
    fi
    
    return 0  # Complete
}

# ============================================================================
# Process each PDB file
# ============================================================================
COMPLETED=0
FAILED=0
SKIPPED=0

for PDB_PATH in "${PDB_FILES[@]}"; do
    PDB_NAME=$(basename "$PDB_PATH" .pdb)
    OUTPUT_DIR="/home/sbali/LigandMPNN/designG4FP/output_${PDB_NAME}"
    
    echo "------------------------------------------------------------------------" | tee -a "$PROGRESS_LOG"
    echo "[$(date +%H:%M:%S)] Processing $((COMPLETED + FAILED + SKIPPED + 1))/$TOTAL_PDBS: $PDB_NAME" | tee -a "$PROGRESS_LOG"
    echo "------------------------------------------------------------------------" | tee -a "$PROGRESS_LOG"
    
    # Check if output already exists and is complete
    if [ $FORCE_RERUN -eq 1 ] && [ -d "$OUTPUT_DIR/01_ligandmpnn" ]; then
        echo "  🗑️  Removing existing output (force rerun mode)..." | tee -a "$PROGRESS_LOG"
        rm -rf "$OUTPUT_DIR/01_ligandmpnn"
    elif check_output_complete "$OUTPUT_DIR"; then
        SEQ_COUNT=$(grep -c "^>" "$OUTPUT_DIR/01_ligandmpnn/seqs"/*.fa 2>/dev/null || echo "0")
        echo "  ⏭️  SKIPPED: Output already exists ($SEQ_COUNT sequences)" | tee -a "$PROGRESS_LOG"
        echo "  Output directory: $OUTPUT_DIR" | tee -a "$PROGRESS_LOG"
        echo "  💡 Use --force-rerun to regenerate" | tee -a "$PROGRESS_LOG"
        ((SKIPPED++))
        echo "" | tee -a "$PROGRESS_LOG"
        continue
    fi
    
    # Restore template config for this run
    cp "$CONFIG_TEMPLATE" "$CONFIG_FILE"
    
    # Update config file with current PDB
    echo "  Updating config for $PDB_NAME..." | tee -a "$PROGRESS_LOG"
    python3 -c "
import yaml
try:
    with open('$CONFIG_FILE', 'r') as f:
        config = yaml.safe_load(f)
    
    # Update input_pdb path
    config['input_pdb'] = '$PDB_PATH'
    
    # Update output directory to be unique per PDB
    config['output_dir'] = '$OUTPUT_DIR/'
    
    with open('$CONFIG_FILE', 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    
    print('Config updated successfully')
except Exception as e:
    print(f'ERROR: {e}')
    exit(1)
" 2>&1 | tee -a "$DETAILED_LOG"
    
    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo "  ❌ ERROR: Failed to update config for $PDB_NAME" | tee -a "$PROGRESS_LOG"
        ((FAILED++))
        echo "" | tee -a "$PROGRESS_LOG"
        continue
    fi
    
    echo "  ✓ Config updated:" | tee -a "$PROGRESS_LOG"
    echo "    Input:  $PDB_PATH" | tee -a "$PROGRESS_LOG"
    echo "    Output: $OUTPUT_DIR/" | tee -a "$PROGRESS_LOG"
    echo "" | tee -a "$PROGRESS_LOG"
    
    # Run LigandMPNN
    echo "  🚀 Starting LigandMPNN (this may take 10-30 minutes)..." | tee -a "$PROGRESS_LOG"
    START_TIME=$(date +%s)
    
    # Add separator to detailed log
    echo "" >> "$DETAILED_LOG"
    echo "========================================================================" >> "$DETAILED_LOG"
    echo "Processing: $PDB_NAME" >> "$DETAILED_LOG"
    echo "Started at: $(date)" >> "$DETAILED_LOG"
    echo "========================================================================" >> "$DETAILED_LOG"
    
    # Run and capture output
    ./01_run_ligandmpnn.sh >> "$DETAILED_LOG" 2>&1
    EXIT_CODE=$?
    
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    MINUTES=$((DURATION / 60))
    SECONDS=$((DURATION % 60))
    
    if [ $EXIT_CODE -eq 0 ] && check_output_complete "$OUTPUT_DIR"; then
        SEQ_COUNT=$(grep -c "^>" "$OUTPUT_DIR/01_ligandmpnn/seqs"/*.fa 2>/dev/null || echo "0")
        echo "  ✅ SUCCESS: Generated $SEQ_COUNT sequences in ${MINUTES}m ${SECONDS}s" | tee -a "$PROGRESS_LOG"
        ((COMPLETED++))
    else
        echo "  ❌ FAILED after ${MINUTES}m ${SECONDS}s" | tee -a "$PROGRESS_LOG"
        echo "     Exit code: $EXIT_CODE" | tee -a "$PROGRESS_LOG"
        echo "     Check detailed log for errors: $DETAILED_LOG" | tee -a "$PROGRESS_LOG"
        ((FAILED++))
    fi
    
    echo "" | tee -a "$PROGRESS_LOG"
done

# ============================================================================
# Restore template config
# ============================================================================
cp "$CONFIG_TEMPLATE" "$CONFIG_FILE"

# ============================================================================
# Final summary
# ============================================================================
echo "============================================================================" | tee -a "$PROGRESS_LOG"
echo "Batch Processing Complete - Finished at $(date)" | tee -a "$PROGRESS_LOG"
echo "============================================================================" | tee -a "$PROGRESS_LOG"
echo "" | tee -a "$PROGRESS_LOG"
echo "Summary:" | tee -a "$PROGRESS_LOG"
echo "  Total PDBs:     $TOTAL_PDBS" | tee -a "$PROGRESS_LOG"
echo "  ✅ Completed:   $COMPLETED" | tee -a "$PROGRESS_LOG"
echo "  ⏭️  Skipped:    $SKIPPED (already complete)" | tee -a "$PROGRESS_LOG"
echo "  ❌ Failed:      $FAILED" | tee -a "$PROGRESS_LOG"
echo "" | tee -a "$PROGRESS_LOG"

# List all output directories
echo "Output directories created:" | tee -a "$PROGRESS_LOG"
for PDB_PATH in "${PDB_FILES[@]}"; do
    PDB_NAME=$(basename "$PDB_PATH" .pdb)
    OUTPUT_DIR="output_${PDB_NAME}"
    if [ -d "$OUTPUT_DIR" ]; then
        if check_output_complete "$OUTPUT_DIR"; then
            SEQ_COUNT=$(grep -c "^>" "$OUTPUT_DIR/01_ligandmpnn/seqs"/*.fa 2>/dev/null || echo "0")
            echo "  ✓ $OUTPUT_DIR ($SEQ_COUNT sequences)" | tee -a "$PROGRESS_LOG"
        else
            echo "  ⚠ $OUTPUT_DIR (incomplete/empty)" | tee -a "$PROGRESS_LOG"
        fi
    fi
done
echo "" | tee -a "$PROGRESS_LOG"

echo "Log files:" | tee -a "$PROGRESS_LOG"
echo "  Progress: $PROGRESS_LOG" | tee -a "$PROGRESS_LOG"
echo "  Detailed: $DETAILED_LOG" | tee -a "$PROGRESS_LOG"
echo "" | tee -a "$PROGRESS_LOG"

if [ $FAILED -gt 0 ]; then
    echo "⚠️  Some jobs failed. Check $DETAILED_LOG for details." | tee -a "$PROGRESS_LOG"
    exit 1
elif [ $COMPLETED -eq 0 ] && [ $SKIPPED -eq $TOTAL_PDBS ]; then
    echo "ℹ️  All outputs already exist. Nothing to process." | tee -a "$PROGRESS_LOG"
else
    echo "🎉 All jobs completed successfully!" | tee -a "$PROGRESS_LOG"
fi

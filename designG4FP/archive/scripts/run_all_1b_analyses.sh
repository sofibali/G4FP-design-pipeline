#!/bin/bash
# Batch run 01b_analyze_pre_filter.py for all design directories
# Usage: bash run_all_1b_analyses.sh

set -e  # Exit on error

echo "=========================================="
echo "Batch Analysis: Pre-Filter Analysis (Step 1b)"
echo "=========================================="
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Find all output directories
OUTPUT_DIRS=(output_G4FP_*)

# Check if any directories were found
if [ ${#OUTPUT_DIRS[@]} -eq 0 ]; then
    echo "❌ ERROR: No output_G4FP_* directories found!"
    exit 1
fi

echo "Found ${#OUTPUT_DIRS[@]} design directories to process:"
for dir in "${OUTPUT_DIRS[@]}"; do
    echo "  - $dir"
done
echo ""

# Counter for tracking progress
TOTAL=${#OUTPUT_DIRS[@]}
CURRENT=0
SUCCESS=0
FAILED=0

# Log file
LOG_FILE="batch_analyze_1b_$(date +%Y%m%d_%H%M%S).log"
echo "Logging to: $LOG_FILE"
echo ""

# Process each directory
for OUTPUT_DIR in "${OUTPUT_DIRS[@]}"; do
    CURRENT=$((CURRENT + 1))
    
    echo "=========================================="
    echo "Processing [$CURRENT/$TOTAL]: $OUTPUT_DIR"
    echo "=========================================="
    
    # Check if the directory exists and has LigandMPNN output
    if [ ! -d "$OUTPUT_DIR/01_ligandmpnn/seqs" ]; then
        echo "⚠️  WARNING: No LigandMPNN output found in $OUTPUT_DIR/01_ligandmpnn/seqs/"
        echo "   Skipping..."
        FAILED=$((FAILED + 1))
        echo "" | tee -a "$LOG_FILE"
        continue
    fi
    
    # Create a temporary config file for this directory
    TEMP_CONFIG="${OUTPUT_DIR}_temp_config.yaml"
    cat > "$TEMP_CONFIG" <<EOF
# Temporary config for $OUTPUT_DIR analysis
output_dir: "$OUTPUT_DIR"
EOF
    
    # Run the analysis
    echo "Running analysis..."
    if python3 01b_analyze_pre_filter.py "$TEMP_CONFIG" 2>&1 | tee -a "$LOG_FILE"; then
        echo "✓ Successfully analyzed $OUTPUT_DIR" | tee -a "$LOG_FILE"
        SUCCESS=$((SUCCESS + 1))
    else
        echo "❌ Failed to analyze $OUTPUT_DIR" | tee -a "$LOG_FILE"
        FAILED=$((FAILED + 1))
    fi
    
    # Clean up temporary config
    rm -f "$TEMP_CONFIG"
    
    echo ""
done

# Final summary
echo "=========================================="
echo "Batch Analysis Complete!"
echo "=========================================="
echo "Total directories:  $TOTAL"
echo "Successful:         $SUCCESS"
echo "Failed/Skipped:     $FAILED"
echo ""
echo "Results saved in each directory under: 01b_pre_filter_analysis/"
echo "Full log saved to: $LOG_FILE"
echo ""

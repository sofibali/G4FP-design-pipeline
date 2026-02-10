#!/bin/bash
# filepath: /home/sbali/LigandMPNN/designG4FP/batch_analyze_pre_filter.sh
# Batch run pre-filter analysis for all LigandMPNN outputs
# Usage: ./batch_analyze_pre_filter.sh

set -e  # Exit on error

# ============================================================================
# Configuration
# ============================================================================
CONFIG_TEMPLATE="pipeline_config_template.yaml"
CONFIG_FILE="pipeline_config.yaml"
PROGRESS_LOG="batch_prefilter_analysis.log"

# ============================================================================
# Initialize log
# ============================================================================
echo "============================================================================" | tee "$PROGRESS_LOG"
echo "Batch Pre-Filter Analysis - Started at $(date)" | tee -a "$PROGRESS_LOG"
echo "============================================================================" | tee -a "$PROGRESS_LOG"
echo "" | tee -a "$PROGRESS_LOG"

# ============================================================================
# Find all LigandMPNN output directories
# ============================================================================
OUTPUT_DIRS=($(ls -d output_*/ 2>/dev/null | sed 's:/$::'))
TOTAL_DIRS=${#OUTPUT_DIRS[@]}

if [ $TOTAL_DIRS -eq 0 ]; then
    echo "ERROR: No output directories found (output_*)" | tee -a "$PROGRESS_LOG"
    echo "       Please run batch_run_ligandmpnn.sh first" | tee -a "$PROGRESS_LOG"
    exit 1
fi

echo "Found $TOTAL_DIRS output directories to analyze:" | tee -a "$PROGRESS_LOG"
for dir in "${OUTPUT_DIRS[@]}"; do
    echo "  - $dir" | tee -a "$PROGRESS_LOG"
done
echo "" | tee -a "$PROGRESS_LOG"

# ============================================================================
# Function to check if analysis already exists
# ============================================================================
check_analysis_exists() {
    local output_dir=$1
    local analysis_dir="${output_dir}/01b_pre_filter_analysis"
    
    # Check if directory exists and has plots
    if [ ! -d "$analysis_dir" ]; then
        return 1  # Doesn't exist
    fi
    
    # Check if all 4 plots exist
    local plot_count=$(ls "$analysis_dir"/*.png 2>/dev/null | wc -l)
    if [ $plot_count -lt 4 ]; then
        return 1  # Incomplete
    fi
    
    # Check if summary CSV exists
    if [ ! -f "$analysis_dir/pre_filter_summary.csv" ]; then
        return 1  # No summary
    fi
    
    return 0  # Complete
}

# ============================================================================
# Process each output directory
# ============================================================================
COMPLETED=0
FAILED=0
SKIPPED=0

for OUTPUT_DIR in "${OUTPUT_DIRS[@]}"; do
    PDB_NAME=$(echo "$OUTPUT_DIR" | sed 's/output_//')
    
    echo "------------------------------------------------------------------------" | tee -a "$PROGRESS_LOG"
    echo "[$(date +%H:%M:%S)] Processing $((COMPLETED + FAILED + SKIPPED + 1))/$TOTAL_DIRS: $PDB_NAME" | tee -a "$PROGRESS_LOG"
    echo "------------------------------------------------------------------------" | tee -a "$PROGRESS_LOG"
    
    # Check if LigandMPNN output exists
    if [ ! -d "$OUTPUT_DIR/01_ligandmpnn/seqs" ]; then
        echo "  ⚠️  WARNING: No LigandMPNN sequences found in $OUTPUT_DIR" | tee -a "$PROGRESS_LOG"
        echo "              Skipping..." | tee -a "$PROGRESS_LOG"
        ((SKIPPED++))
        echo "" | tee -a "$PROGRESS_LOG"
        continue
    fi
    
    # Check if analysis already exists
    if check_analysis_exists "$OUTPUT_DIR"; then
        echo "  ⏭️  SKIPPED: Analysis already exists" | tee -a "$PROGRESS_LOG"
        echo "     Location: $OUTPUT_DIR/01b_pre_filter_analysis" | tee -a "$PROGRESS_LOG"
        ((SKIPPED++))
        echo "" | tee -a "$PROGRESS_LOG"
        continue
    fi
    
    # Create temporary config for this analysis
    TEMP_CONFIG="${OUTPUT_DIR}/temp_config.yaml"
    
    # Check if template exists, otherwise use current config
    if [ -f "$CONFIG_TEMPLATE" ]; then
        cp "$CONFIG_TEMPLATE" "$TEMP_CONFIG"
    else
        cp "$CONFIG_FILE" "$TEMP_CONFIG"
    fi
    
    # Update config with correct output directory
    python3 -c "
import yaml
try:
    with open('$TEMP_CONFIG', 'r') as f:
        config = yaml.safe_load(f)
    
    # Update output directory
    config['output_dir'] = '$OUTPUT_DIR/'
    
    with open('$TEMP_CONFIG', 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    
    print('Config created successfully')
except Exception as e:
    print(f'ERROR: {e}')
    exit(1)
" 2>&1
    
    if [ $? -ne 0 ]; then
        echo "  ❌ ERROR: Failed to create config for $PDB_NAME" | tee -a "$PROGRESS_LOG"
        rm -f "$TEMP_CONFIG"
        ((FAILED++))
        echo "" | tee -a "$PROGRESS_LOG"
        continue
    fi
    
    echo "  ✓ Config created: $TEMP_CONFIG" | tee -a "$PROGRESS_LOG"
    echo "  🚀 Running pre-filter analysis..." | tee -a "$PROGRESS_LOG"
    
    # Run analysis
    START_TIME=$(date +%s)
    
    if python 01b_analyze_pre_filter.py "$TEMP_CONFIG" 2>&1 | tee -a "$PROGRESS_LOG"; then
        END_TIME=$(date +%s)
        DURATION=$((END_TIME - START_TIME))
        
        echo "  ✅ SUCCESS: Analysis completed in ${DURATION}s" | tee -a "$PROGRESS_LOG"
        
        # Show summary if available
        SUMMARY_FILE="$OUTPUT_DIR/01b_pre_filter_analysis/pre_filter_summary.csv"
        if [ -f "$SUMMARY_FILE" ]; then
            # Extract key metrics
            TOTAL_SEQS=$(python3 -c "import pandas as pd; df=pd.read_csv('$SUMMARY_FILE'); print(int(df['total_sequences'].iloc[0]))")
            
            echo "     Total sequences: $TOTAL_SEQS" | tee -a "$PROGRESS_LOG"
        fi
        
        ((COMPLETED++))
    else
        END_TIME=$(date +%s)
        DURATION=$((END_TIME - START_TIME))
        
        echo "  ❌ FAILED after ${DURATION}s" | tee -a "$PROGRESS_LOG"
        ((FAILED++))
    fi
    
    # Clean up temporary config
    rm -f "$TEMP_CONFIG"
    
    echo "" | tee -a "$PROGRESS_LOG"
done

# ============================================================================
# Final summary
# ============================================================================
echo "============================================================================" | tee -a "$PROGRESS_LOG"
echo "Batch Pre-Filter Analysis Complete - Finished at $(date)" | tee -a "$PROGRESS_LOG"
echo "============================================================================" | tee -a "$PROGRESS_LOG"
echo "" | tee -a "$PROGRESS_LOG"
echo "Summary:" | tee -a "$PROGRESS_LOG"
echo "  Total directories:  $TOTAL_DIRS" | tee -a "$PROGRESS_LOG"
echo "  ✅ Completed:       $COMPLETED" | tee -a "$PROGRESS_LOG"
echo "  ⏭️  Skipped:        $SKIPPED (already complete or no data)" | tee -a "$PROGRESS_LOG"
echo "  ❌ Failed:          $FAILED" | tee -a "$PROGRESS_LOG"
echo "" | tee -a "$PROGRESS_LOG"

# List all analysis directories
echo "Analysis directories:" | tee -a "$PROGRESS_LOG"
for OUTPUT_DIR in "${OUTPUT_DIRS[@]}"; do
    ANALYSIS_DIR="${OUTPUT_DIR}/01b_pre_filter_analysis"
    if [ -d "$ANALYSIS_DIR" ]; then
        PLOT_COUNT=$(ls "$ANALYSIS_DIR"/*.png 2>/dev/null | wc -l)
        if [ -f "$ANALYSIS_DIR/pre_filter_summary.csv" ]; then
            echo "  ✓ $ANALYSIS_DIR ($PLOT_COUNT plots)" | tee -a "$PROGRESS_LOG"
        else
            echo "  ⚠ $ANALYSIS_DIR (incomplete)" | tee -a "$PROGRESS_LOG"
        fi
    fi
done
echo "" | tee -a "$PROGRESS_LOG"

if [ $FAILED -gt 0 ]; then
    echo "⚠️  Some analyses failed. Check log for details." | tee -a "$PROGRESS_LOG"
    exit 1
elif [ $COMPLETED -eq 0 ] && [ $SKIPPED -eq $TOTAL_DIRS ]; then
    echo "ℹ️  All analyses already exist. Nothing to process." | tee -a "$PROGRESS_LOG"
else
    echo "🎉 All analyses completed successfully!" | tee -a "$PROGRESS_LOG"
fi
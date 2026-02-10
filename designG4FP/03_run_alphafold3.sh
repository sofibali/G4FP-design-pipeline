#!/bin/bash
# Step 3: Run AlphaFold3 predictions
# Usage: ./03_run_alphafold3.sh [config_file] [bound|apo|both]

set -e  # Exit on error

# Load configuration
CONFIG_FILE="${1:-pipeline_config.yaml}"
RUN_MODE="${2:-both}"  # Options: bound, apo, or both

echo "============================================================================"
echo "STEP 3: Running AlphaFold3 Structure Predictions"
echo "============================================================================"
echo ""

# Parse YAML config
parse_config() {
    python3 -c "
import yaml
with open('$CONFIG_FILE') as f:
    config = yaml.safe_load(f)
print(f\"OUTPUT_DIR='{config['output_dir']}'\")
print(f\"DB_DIR='{config['alphafold3']['db_dir']}'\")
print(f\"MODEL_DIR='{config['alphafold3']['model_dir']}'\")
print(f\"GPU_ID={config['alphafold3']['gpu_id']}\")
"
}

# Load config into environment
eval $(parse_config)

# Setup directories
BOUND_INPUT_DIR="$OUTPUT_DIR/02_alphafold3_inputs_bound"
APO_INPUT_DIR="$OUTPUT_DIR/02_alphafold3_inputs_apo"
BOUND_OUTPUT_DIR="$OUTPUT_DIR/03_alphafold3_predictions_bound"
APO_OUTPUT_DIR="$OUTPUT_DIR/03_alphafold3_predictions_apo"

echo "Configuration:"
echo "  Run mode:       $RUN_MODE"
echo "  GPU ID:         $GPU_ID"
echo "  Database dir:   $DB_DIR"
echo "  Model dir:      $MODEL_DIR"
echo ""

# Function to run AlphaFold3 on a single JSON
run_af3() {
    local json_path=$1
    local output_dir=$2
    local design_name=$(basename "$json_path" .json)
    
    echo "  Running: $design_name"
    
    # Source SBGrid environment
    source /programs/sbgrid.shrc
    
    # Set CUDA device
    export CUDA_VISIBLE_DEVICES=$GPU_ID
    export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas
    
    # Run AlphaFold3
    /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
        --db_dir "$DB_DIR" \
        --model_dir "$MODEL_DIR" \
        --output_dir "$output_dir/${design_name}" \
        --json_path "$json_path" 2>&1 | tee "$output_dir/${design_name}.log"
    
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "    ✓ Completed"
    else
        echo "    ❌ Failed (see log: $output_dir/${design_name}.log)"
        return 1
    fi
}

# Function to run batch of AF3 predictions
run_af3_batch() {
    local input_dir=$1
    local output_dir=$2
    local state_name=$3
    
    echo "----------------------------------------"
    echo "Running $state_name predictions"
    echo "----------------------------------------"
    
    # Check if input directory exists
    if [ ! -d "$input_dir" ]; then
        echo "⚠ Input directory not found: $input_dir"
        echo "  Skipping $state_name predictions"
        return 0
    fi
    
    # Create output directory
    mkdir -p "$output_dir"
    
    # Count JSON files
    json_count=$(ls -1 "$input_dir"/*.json 2>/dev/null | wc -l)
    
    if [ $json_count -eq 0 ]; then
        echo "⚠ No JSON files found in $input_dir"
        return 0
    fi
    
    echo "Found $json_count input files"
    echo ""
    
    # Loop through JSON files
    current=0
    for json_file in "$input_dir"/*.json; do
        current=$((current + 1))
        echo "[$current/$json_count] Processing $(basename "$json_file")"
        
        # Run AF3 (continue on failure)
        run_af3 "$json_file" "$output_dir" || echo "  Continuing with next design..."
        echo ""
    done
    
    echo "✓ Completed $state_name predictions ($current/$json_count processed)"
    echo ""
}

# Main execution
echo "Starting AlphaFold3 predictions..."
echo "NOTE: Each design takes ~30min-2hrs depending on GPU and sequence length"
echo ""

# Run bound (with ligand) predictions
if [ "$RUN_MODE" = "bound" ] || [ "$RUN_MODE" = "both" ]; then
    run_af3_batch "$BOUND_INPUT_DIR" "$BOUND_OUTPUT_DIR" "BOUND (protein+DNA)"
fi

# Run apo (protein-only) predictions
if [ "$RUN_MODE" = "apo" ] || [ "$RUN_MODE" = "both" ]; then
    run_af3_batch "$APO_INPUT_DIR" "$APO_OUTPUT_DIR" "APO (protein-only)"
fi

echo "============================================================================"
echo "✓ AlphaFold3 predictions complete!"
echo "============================================================================"
echo ""

# Summary statistics
if [ "$RUN_MODE" = "bound" ] || [ "$RUN_MODE" = "both" ]; then
    if [ -d "$BOUND_OUTPUT_DIR" ]; then
        bound_count=$(ls -1d "$BOUND_OUTPUT_DIR"/design_* 2>/dev/null | wc -l)
        echo "Bound predictions:  $bound_count designs in $BOUND_OUTPUT_DIR"
    fi
fi

if [ "$RUN_MODE" = "apo" ] || [ "$RUN_MODE" = "both" ]; then
    if [ -d "$APO_OUTPUT_DIR" ]; then
        apo_count=$(ls -1d "$APO_OUTPUT_DIR"/design_* 2>/dev/null | wc -l)
        echo "Apo predictions:    $apo_count designs in $APO_OUTPUT_DIR"
    fi
fi

echo ""
echo "Next step: Run analysis scripts (04a, 04b, 04c)"
echo ""

# ============================================================================
# ALTERNATIVE: Simple loop for manual execution
# ============================================================================
# If you prefer to run AF3 manually, uncomment and use this loop:
#
# for json in $BOUND_INPUT_DIR/*.json; do
#     design=$(basename "$json" .json)
#     echo "Running $design..."
#     run_alphafold.py \
#         --db_dir "$DB_DIR" \
#         --model_dir "$MODEL_DIR" \
#         --output_dir "$BOUND_OUTPUT_DIR/$design" \
#         --json_path "$json"
# done

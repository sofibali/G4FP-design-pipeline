#!/bin/bash
# Step 3: Run AlphaFold3 predictions in parallel
# Usage: ./03_run_alphafold3_parallel.sh [config_file] [bound|apo|both] [num_parallel]

set -e  # Exit on error

# Load configuration
CONFIG_FILE="${1:-pipeline_config.yaml}"
RUN_MODE="${2:-both}"  # Options: bound, apo, or both
NUM_PARALLEL="${3:-2}"  # Number of parallel jobs (default: 2 - conservative for safety)

echo "============================================================================"
echo "STEP 3: Running AlphaFold3 Structure Predictions (PARALLEL MODE)"
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
echo "  Parallel jobs:  $NUM_PARALLEL"
echo "  Database dir:   $DB_DIR"
echo "  Model dir:      $MODEL_DIR"
echo ""

# Function to run AlphaFold3 on a single JSON
run_af3() {
    local json_path=$1
    local output_dir=$2
    local design_name=$(basename "$json_path" .json)
    local job_output_dir="$output_dir/${design_name}"
    
    # Skip if already completed
    if [ -f "$job_output_dir/model_cif/model_0.cif" ]; then
        echo "  ⏭  Skipping $design_name (already completed)"
        return 0
    fi
    
    echo "[$(date '+%H:%M:%S')] 🔄 Running: $design_name (PID: $$)"
    
    # Source SBGrid environment (suppress verbose output)
    source /programs/sbgrid.shrc 2>/dev/null
    
    # Set CUDA device
    export CUDA_VISIBLE_DEVICES=$GPU_ID
    export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas
    
    # Run AlphaFold3
    /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
        --db_dir "$DB_DIR" \
        --model_dir "$MODEL_DIR" \
        --output_dir "$job_output_dir" \
        --json_path "$json_path" &> "$job_output_dir.log"
    
    local exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        echo "[$(date '+%H:%M:%S')] ✓ Completed: $design_name"
    else
        echo "[$(date '+%H:%M:%S')] ❌ Failed: $design_name (exit code: $exit_code)"
        echo "     Check log: $job_output_dir.log"
        # Check for common errors
        if grep -q "dialect.*version" "$job_output_dir.log" 2>/dev/null; then
            echo "     Error: JSON missing 'dialect' and 'version' fields"
        fi
        return 1
    fi
}

# Export function and ALL required variables for parallel execution
export -f run_af3
export GPU_ID DB_DIR MODEL_DIR

# Function to run batch of AF3 predictions in parallel
run_af3_batch_parallel() {
    local input_dir=$1
    local output_dir=$2
    local state_name=$3
    
    echo "----------------------------------------"
    echo "Running $state_name predictions (PARALLEL)"
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
    echo "Running $NUM_PARALLEL jobs in parallel on GPU $GPU_ID"
    echo ""
    
    # Check if GNU parallel is available
    if command -v parallel &> /dev/null; then
        echo "✓ Using GNU Parallel for job distribution"
        echo ""
        
        # Run with GNU parallel - much more efficient
        ls "$input_dir"/*.json | \
            parallel -j $NUM_PARALLEL \
                     --bar \
                     --joblog "$output_dir/parallel_jobs.log" \
                     --halt soon,fail=10% \
                     run_af3 {} "$output_dir"
        
    else
        echo "⚠ GNU Parallel not found, using bash background jobs"
        echo "  Install with: conda install -c conda-forge parallel"
        echo "  or: sudo apt-get install parallel (much better performance)"
        echo ""
        
        # Fallback: Use bash background jobs with manual semaphore
        local running=0
        local completed=0
        local failed=0
        
        for json_file in "$input_dir"/*.json; do
            # Wait if we've hit the parallel limit
            while [ $(jobs -r | wc -l) -ge $NUM_PARALLEL ]; do
                sleep 2
            done
            
            # Launch job in background
            (
                if run_af3 "$json_file" "$output_dir"; then
                    exit 0
                else
                    exit 1
                fi
            ) &
            
            completed=$((completed + 1))
            echo "  Launched: $(basename "$json_file") [$completed/$json_count]"
            
            # Small delay to avoid overwhelming the system
            sleep 1
        done
        
        # Wait for all remaining jobs
        echo ""
        echo "Waiting for all jobs to complete..."
        
        # Wait and count failures
        for job in $(jobs -p); do
            if ! wait $job; then
                failed=$((failed + 1))
            fi
        done
        
        if [ $failed -gt 0 ]; then
            echo "⚠ $failed jobs failed"
        fi
    fi
    
    echo ""
    echo "✓ Completed $state_name predictions"
    echo ""
}

# Main execution
echo "Starting AlphaFold3 predictions in parallel..."
echo ""
echo "⚠ IMPORTANT NOTES:"
echo "  • Each job loads the model into GPU memory (~5-8GB per job)"
echo "  • During MSA search: jobs are CPU-bound, GPU mostly idle"
echo "  • During inference: GPU utilization spikes to 80-100%"
echo "  • If you get OOM errors, reduce NUM_PARALLEL"
echo ""
echo "Recommended parallel jobs by GPU memory:"
echo "  • 12GB VRAM: 1-2 jobs"
echo "  • 16-24GB:   2-4 jobs (current: $NUM_PARALLEL)"
echo "  • 32-48GB:   4-6 jobs"
echo "  • 80GB A100: 6-10 jobs"
echo ""

# Run bound (with ligand) predictions
if [ "$RUN_MODE" = "bound" ] || [ "$RUN_MODE" = "both" ]; then
    run_af3_batch_parallel "$BOUND_INPUT_DIR" "$BOUND_OUTPUT_DIR" "BOUND (protein+DNA)"
fi

# Run apo (protein-only) predictions
if [ "$RUN_MODE" = "apo" ] || [ "$RUN_MODE" = "both" ]; then
    run_af3_batch_parallel "$APO_INPUT_DIR" "$APO_OUTPUT_DIR" "APO (protein-only)"
fi

echo "============================================================================"
echo "✓ AlphaFold3 predictions complete!"
echo "============================================================================"
echo ""

# Summary statistics
total_success=0
total_failed=0

if [ "$RUN_MODE" = "bound" ] || [ "$RUN_MODE" = "both" ]; then
    if [ -d "$BOUND_OUTPUT_DIR" ]; then
        bound_success=$(find "$BOUND_OUTPUT_DIR" -name "model_0.cif" -type f 2>/dev/null | wc -l)
        bound_total=$(ls -1 "$BOUND_INPUT_DIR"/*.json 2>/dev/null | wc -l)
        bound_failed=$((bound_total - bound_success))
        total_success=$((total_success + bound_success))
        total_failed=$((total_failed + bound_failed))
        
        echo "Bound predictions:  $bound_success/$bound_total successful"
        [ $bound_failed -gt 0 ] && echo "                    $bound_failed failed (check logs in $BOUND_OUTPUT_DIR)"
    fi
fi

if [ "$RUN_MODE" = "apo" ] || [ "$RUN_MODE" = "both" ]; then
    if [ -d "$APO_OUTPUT_DIR" ]; then
        apo_success=$(find "$APO_OUTPUT_DIR" -name "model_0.cif" -type f 2>/dev/null | wc -l)
        apo_total=$(ls -1 "$APO_INPUT_DIR"/*.json 2>/dev/null | wc -l)
        apo_failed=$((apo_total - apo_success))
        total_success=$((total_success + apo_success))
        total_failed=$((total_failed + apo_failed))
        
        echo "Apo predictions:    $apo_success/$apo_total successful"
        [ $apo_failed -gt 0 ] && echo "                    $apo_failed failed (check logs in $APO_OUTPUT_DIR)"
    fi
fi

echo ""
echo "Total: $total_success successful, $total_failed failed"
echo ""

if [ $total_failed -gt 0 ]; then
    echo "⚠ Some predictions failed. Check individual log files for errors."
    echo ""
    echo "  Common issues:"
    echo "    • GPU out of memory → reduce NUM_PARALLEL"
    echo "    • JSON missing dialect/version fields → fix input generation script"
    echo "    • MSA search timeout → increase database search time limits"
    echo "    • Invalid sequence → check FASTA formatting"
    echo ""
    echo "  Quick diagnosis:"
    echo "    # Check first failure:"
    if [ "$RUN_MODE" = "bound" ] || [ "$RUN_MODE" = "both" ]; then
        echo "    tail -50 $BOUND_OUTPUT_DIR/*.log | grep -A 5 'Error\\|Failed\\|Traceback'"
    fi
fi

echo ""
echo "Next step: Run analysis scripts (04a, 04b, 04c)"
echo ""

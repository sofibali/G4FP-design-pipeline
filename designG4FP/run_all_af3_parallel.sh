#!/bin/bash
# Run AlphaFold3 efficiently for all G4FP structures across two GPUs
#
# Two-phase approach for maximum GPU efficiency:
#   Phase 1: Data pipeline (MSA/template search) -- CPU only, highly parallel
#   Phase 2: Structure inference -- GPU, model loads ONCE per GPU via --input_dir
#
# Key optimizations vs naive approach:
#   - Model loads once per GPU (not once per prediction)
#   - JAX compilation cached to disk (~2.5x speedup after first compile)
#   - Flash attention via Triton (fastest for Ampere+ GPUs)
#   - Data pipeline separated from inference (CPU parallelism)
#   - Bucket-based compilation reuse (all inputs ~300 tokens -> 512 bucket)
#
# Usage:
#   ./run_all_af3_parallel.sh                         # defaults: 2 GPUs
#   ./run_all_af3_parallel.sh --gpu0 0 --gpu1 1       # specify GPU IDs
#   ./run_all_af3_parallel.sh --msa-parallel 8         # MSA jobs in parallel
#   ./run_all_af3_parallel.sh --skip-msa               # skip phase 1 (already done)
#   ./run_all_af3_parallel.sh --node1-only              # run only GPU 0's share
#   ./run_all_af3_parallel.sh --node2-only              # run only GPU 1's share
#   ./run_all_af3_parallel.sh --dry-run                 # show plan without running
#   ./run_all_af3_parallel.sh --mode bound              # only bound state
#
# Time estimate (A100 80GB, ~300 residue protein):
#   Phase 1 (MSA): ~30 min/seq, 8 parallel = ~62 hours total -> ~8 hours wall
#   Phase 2 (inference): ~5 min/prediction (5 seeds), model loads once
#     1000 predictions per GPU -> ~83 hours / GPU
#   Total: ~3-4 days with 2 GPUs (dominated by inference)
#
# For monitoring:
#   tail -f run_af3_phase1.log                    # MSA progress
#   tail -f run_af3_gpu0_inference.log             # GPU 0 inference
#   tail -f run_af3_gpu1_inference.log             # GPU 1 inference
#   find output_G4FP_*/03_alphafold3_predictions_* -name '*.cif' | wc -l

set +e  # Continue even if individual predictions fail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================================
# Defaults
# ============================================================================
GPU0=0
GPU1=1
MSA_PARALLEL=20         # Parallel MSA/data pipeline jobs (CPU only, ~80GB RAM each)
RUN_MODE="both"          # bound, apo, both
NODE1_ONLY=0
NODE2_ONLY=0
DRY_RUN=0
SKIP_MSA=0
SKIP_INFERENCE=0
JAX_CACHE_DIR="$SCRIPT_DIR/.jax_cache"
NUM_DIFFUSION_SAMPLES=5
NUM_SEEDS=""             # empty = use seeds from JSON

# ============================================================================
# Parse arguments
# ============================================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        --gpu0)            GPU0="$2"; shift 2 ;;
        --gpu1)            GPU1="$2"; shift 2 ;;
        --msa-parallel)    MSA_PARALLEL="$2"; shift 2 ;;
        --mode)            RUN_MODE="$2"; shift 2 ;;
        --node1-only)      NODE1_ONLY=1; shift ;;
        --node2-only)      NODE2_ONLY=1; shift ;;
        --skip-msa)        SKIP_MSA=1; shift ;;
        --skip-inference)  SKIP_INFERENCE=1; shift ;;
        --msa-only)        SKIP_INFERENCE=1; shift ;;
        --inference-only)  SKIP_MSA=1; shift ;;
        --dry-run)         DRY_RUN=1; shift ;;
        --diffusion-samples) NUM_DIFFUSION_SAMPLES="$2"; shift 2 ;;
        --num-seeds)       NUM_SEEDS="$2"; shift 2 ;;
        --jobs-per-node)   shift 2 ;;  # ignored (legacy compat)
        *)                 echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# ============================================================================
# Environment setup
# ============================================================================
source /programs/sbgrid.shrc 2>/dev/null
export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas

# Create JAX compilation cache directory
mkdir -p "$JAX_CACHE_DIR"

# ============================================================================
# Discover all structures with AF3 JSONs ready
# ============================================================================
ALL_STRUCTURES=()
for d in "$SCRIPT_DIR"/output_G4FP_*/; do
    if [ -d "$d/02_alphafold3_inputs_bound" ]; then
        ALL_STRUCTURES+=("$(basename "$d")")
    fi
done

N_STRUCTS=${#ALL_STRUCTURES[@]}
if [ "$N_STRUCTS" -eq 0 ]; then
    echo "ERROR: No output directories with AF3 inputs found."
    echo "       Run 02_filter_and_prepare_af3.py first."
    exit 1
fi

# Split structures between two nodes (roughly equal)
HALF=$(( (N_STRUCTS + 1) / 2 ))
NODE0_STRUCTS=("${ALL_STRUCTURES[@]:0:$HALF}")
NODE1_STRUCTS=("${ALL_STRUCTURES[@]:$HALF}")

# Count total predictions
count_predictions() {
    local total=0
    for s in "$@"; do
        local d="$SCRIPT_DIR/$s"
        if [ "$RUN_MODE" = "bound" ] || [ "$RUN_MODE" = "both" ]; then
            total=$((total + $(ls "$d/02_alphafold3_inputs_bound/"*.json 2>/dev/null | wc -l)))
        fi
        if [ "$RUN_MODE" = "apo" ] || [ "$RUN_MODE" = "both" ]; then
            total=$((total + $(ls "$d/02_alphafold3_inputs_apo/"*.json 2>/dev/null | wc -l)))
        fi
    done
    echo $total
}

TOTAL_PREDS=$(count_predictions "${ALL_STRUCTURES[@]}")
NODE0_PREDS=$(count_predictions "${NODE0_STRUCTS[@]}")
NODE1_PREDS=$(count_predictions "${NODE1_STRUCTS[@]}")

# ============================================================================
# Time estimate (optimized)
# ============================================================================
# Phase 1: MSA ~30 min per sequence, parallelized
MSA_MIN_PER_SEQ=30
MSA_TOTAL_MIN=$(( TOTAL_PREDS * MSA_MIN_PER_SEQ ))
MSA_WALL_MIN=$(( MSA_TOTAL_MIN / MSA_PARALLEL ))
MSA_WALL_HRS=$(( MSA_WALL_MIN / 60 ))

# Phase 2: Inference ~5 min per prediction (5 seeds @ ~60s each) on A100
# First prediction per bucket has ~5 min compilation overhead
INF_MIN_PER_PRED=5
NODE0_INF_HRS=$(( NODE0_PREDS * INF_MIN_PER_PRED / 60 ))
NODE1_INF_HRS=$(( NODE1_PREDS * INF_MIN_PER_PRED / 60 ))
INF_WALL_HRS=$(( NODE0_INF_HRS > NODE1_INF_HRS ? NODE0_INF_HRS : NODE1_INF_HRS ))

echo "============================================================================"
echo "AlphaFold3 Batch Predictions -- Optimized 2-Phase Pipeline"
echo "============================================================================"
echo ""
echo "Structures found:      $N_STRUCTS"
echo "Run mode:              $RUN_MODE"
echo "Total predictions:     $TOTAL_PREDS"
echo "JAX cache:             $JAX_CACHE_DIR"
echo ""
echo "--- Phase 1: Data Pipeline (CPU, $MSA_PARALLEL parallel) ---"
if [ $SKIP_MSA -eq 1 ]; then
    echo "  SKIPPED (--skip-msa)"
else
    echo "  Est. wall time: ~${MSA_WALL_HRS} hours ($TOTAL_PREDS seqs / $MSA_PARALLEL parallel)"
fi
echo ""
echo "--- Phase 2: Inference (2 GPUs, model loads once per GPU) ---"
echo "  GPU 0 ($GPU0): ${NODE0_PREDS} predictions"
echo "    Structures: ${NODE0_STRUCTS[*]}" | sed 's/output_G4FP_//g'
echo "    Est. time: ~${NODE0_INF_HRS} hours"
echo "  GPU 1 ($GPU1): ${NODE1_PREDS} predictions"
echo "    Structures: ${NODE1_STRUCTS[*]}" | sed 's/output_G4FP_//g'
echo "    Est. time: ~${NODE1_INF_HRS} hours"
echo ""
echo "============================================================================"
echo "ESTIMATED TOTAL WALL TIME:"
if [ $SKIP_MSA -eq 0 ]; then
    echo "  Phase 1 (MSA):       ~${MSA_WALL_HRS} hours"
fi
echo "  Phase 2 (inference): ~${INF_WALL_HRS} hours"
echo "  Combined:            ~$(( MSA_WALL_HRS + INF_WALL_HRS )) hours (~$(( (MSA_WALL_HRS + INF_WALL_HRS + 23) / 24 )) days)"
echo "============================================================================"
echo ""

if [ $DRY_RUN -eq 1 ]; then
    echo "[DRY RUN] Would run the above. Exiting."
    exit 0
fi

START_TIME=$(date +%s)

# ============================================================================
# PHASE 1: Data Pipeline (MSA/template search) -- CPU only
# ============================================================================
# Runs --norun_inference to do only MSA search. Produces _data.json files
# that contain pre-computed MSAs for the inference phase.
# This is CPU-bound so we run many in parallel without needing GPUs.
# ============================================================================

run_data_pipeline_for_json() {
    local json_file=$1
    local output_dir=$2
    local design_name=$(basename "$json_file" .json)

    # Check if data pipeline already completed (look for _data.json in output)
    if ls "$output_dir/$design_name"/*_data.json &>/dev/null 2>&1; then
        return 0
    fi

    mkdir -p "$output_dir/$design_name"

    /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
        --db_dir /mnt/alphafold3 \
        --model_dir /mnt/alphafold3 \
        --output_dir "$output_dir/$design_name" \
        --json_path "$json_file" \
        --norun_inference \
        &> "$output_dir/${design_name}_msa.log"
}

if [ $SKIP_MSA -eq 0 ]; then
    echo ""
    echo "============================================================"
    echo "PHASE 1: Data Pipeline (MSA search, $MSA_PARALLEL parallel)"
    echo "Started: $(date)"
    echo "============================================================"
    echo ""

    export -f run_data_pipeline_for_json

    ACTIVE_JOBS=0
    TOTAL_SUBMITTED=0
    TOTAL_SKIPPED=0

    for struct_name in "${ALL_STRUCTURES[@]}"; do
        for state in bound apo; do
            if [ "$RUN_MODE" != "both" ] && [ "$RUN_MODE" != "$state" ]; then
                continue
            fi

            input_dir="$SCRIPT_DIR/$struct_name/02_alphafold3_inputs_${state}"
            output_dir="$SCRIPT_DIR/$struct_name/03_alphafold3_predictions_${state}"

            if [ ! -d "$input_dir" ]; then
                continue
            fi

            mkdir -p "$output_dir"

            for json_file in "$input_dir"/*.json; do
                design_name=$(basename "$json_file" .json)

                # Skip if data pipeline already done
                if ls "$output_dir/$design_name"/*_data.json &>/dev/null 2>&1; then
                    ((TOTAL_SKIPPED++))
                    continue
                fi

                # Also skip if inference already completed (CIF exists)
                if ls "$output_dir/$design_name"/*_model*.cif &>/dev/null 2>&1 || \
                   ls "$output_dir/$design_name"/model_cif/model_0.cif &>/dev/null 2>&1; then
                    ((TOTAL_SKIPPED++))
                    continue
                fi

                # Launch data pipeline in background
                run_data_pipeline_for_json "$json_file" "$output_dir" &
                ((ACTIVE_JOBS++))
                ((TOTAL_SUBMITTED++))

                # Throttle to MSA_PARALLEL concurrent jobs
                if [ $ACTIVE_JOBS -ge $MSA_PARALLEL ]; then
                    wait -n 2>/dev/null || wait
                    ((ACTIVE_JOBS--))
                fi

                # Progress every 50 submissions
                if [ $((TOTAL_SUBMITTED % 50)) -eq 0 ]; then
                    echo "[$(date '+%H:%M:%S')] Phase 1: $TOTAL_SUBMITTED submitted, $TOTAL_SKIPPED skipped"
                fi
            done
        done
    done

    # Wait for remaining MSA jobs
    wait
    echo ""
    echo "[$(date '+%H:%M:%S')] Phase 1 complete: $TOTAL_SUBMITTED MSA jobs run, $TOTAL_SKIPPED skipped"
    echo ""
else
    echo ""
    echo "PHASE 1: SKIPPED (--skip-msa)"
    echo ""
fi

# ============================================================================
# PHASE 2: Structure Inference -- GPU, model loads ONCE per GPU
# ============================================================================
# For each GPU, we create a temporary directory of symlinks to all _data.json
# files assigned to that GPU. Then run AF3 with --input_dir so the model loads
# once and processes all inputs sequentially.
#
# If no _data.json exists (MSA skipped or failed), we fall back to running
# the full pipeline for those inputs.
# ============================================================================

prepare_inference_batch() {
    # Collects _data.json paths OR original JSON paths for a set of structures
    # Returns: populates a temp dir with symlinks
    local batch_dir=$1
    shift
    local structs=("$@")
    local n_ready=0
    local n_fallback=0

    mkdir -p "$batch_dir"

    for struct_name in "${structs[@]}"; do
        for state in bound apo; do
            if [ "$RUN_MODE" != "both" ] && [ "$RUN_MODE" != "$state" ]; then
                continue
            fi

            local input_dir="$SCRIPT_DIR/$struct_name/02_alphafold3_inputs_${state}"
            local output_dir="$SCRIPT_DIR/$struct_name/03_alphafold3_predictions_${state}"

            if [ ! -d "$input_dir" ]; then
                continue
            fi

            for json_file in "$input_dir"/*.json; do
                local design_name=$(basename "$json_file" .json)
                local pred_dir="$output_dir/$design_name"

                # Skip if inference already completed
                if ls "$pred_dir"/*_model*.cif &>/dev/null 2>&1 || \
                   ls "$pred_dir"/model_cif/model_0.cif &>/dev/null 2>&1; then
                    continue
                fi

                # Prefer _data.json (MSA pre-computed), fall back to original
                local data_json=$(ls "$pred_dir"/*_data.json 2>/dev/null | head -1)
                if [ -n "$data_json" ]; then
                    ln -sf "$data_json" "$batch_dir/${struct_name}_${design_name}_data.json"
                    ((n_ready++))
                else
                    ln -sf "$json_file" "$batch_dir/${struct_name}_${design_name}.json"
                    ((n_fallback++))
                fi
            done
        done
    done

    echo "$n_ready pre-computed, $n_fallback need full pipeline"
}

run_inference_batch() {
    local gpu_id=$1
    local node_label=$2
    local batch_dir=$3

    local n_inputs=$(ls "$batch_dir"/*.json 2>/dev/null | wc -l)
    if [ "$n_inputs" -eq 0 ]; then
        echo "[$(date '+%H:%M:%S')] GPU $node_label: No inputs to process"
        return 0
    fi

    echo "[$(date '+%H:%M:%S')] GPU $node_label ($gpu_id): Starting inference on $n_inputs inputs"
    echo "  Model loads once, JAX cache: $JAX_CACHE_DIR"
    echo "  Flash attention: triton"
    echo ""

    local inf_output_dir="$SCRIPT_DIR/.inference_output_gpu${node_label}"
    mkdir -p "$inf_output_dir"

    # Build extra flags
    local extra_flags=""
    if [ -n "$NUM_SEEDS" ]; then
        extra_flags="$extra_flags --num_seeds $NUM_SEEDS"
    fi

    # Check if all inputs are _data.json (pre-computed MSA) -> skip data pipeline
    local n_data=$(ls "$batch_dir"/*_data.json 2>/dev/null | wc -l)
    local n_raw=$(ls "$batch_dir"/*.json 2>/dev/null | wc -l)
    n_raw=$((n_raw - n_data))

    local pipeline_flag=""
    if [ "$n_raw" -eq 0 ] && [ "$n_data" -gt 0 ]; then
        pipeline_flag="--norun_data_pipeline"
        echo "  All inputs have pre-computed MSAs, skipping data pipeline"
    elif [ "$n_data" -gt 0 ] && [ "$n_raw" -gt 0 ]; then
        echo "  WARNING: Mixed inputs ($n_data pre-computed, $n_raw raw). Running full pipeline."
        echo "  (Pre-computed MSAs will be reused where available)"
        pipeline_flag=""
    fi

    export CUDA_VISIBLE_DEVICES=$gpu_id

    /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
        --db_dir /mnt/alphafold3 \
        --model_dir /mnt/alphafold3 \
        --output_dir "$inf_output_dir" \
        --input_dir "$batch_dir" \
        --jax_compilation_cache_dir "$JAX_CACHE_DIR" \
        --flash_attention_implementation triton \
        --num_diffusion_samples $NUM_DIFFUSION_SAMPLES \
        $pipeline_flag \
        $extra_flags

    local exit_code=$?

    # Move results back to the proper per-structure output directories
    echo "[$(date '+%H:%M:%S')] GPU $node_label: Moving results to output directories..."
    for result_dir in "$inf_output_dir"/*/; do
        local result_name=$(basename "$result_dir")
        # Parse structure name and design name from the symlink name
        # Format: output_G4FP_<name>_design_NNNN_<state>[_data]
        # We need to find the original structure and place results correctly

        for struct_name in "${ALL_STRUCTURES[@]}"; do
            for state in bound apo; do
                local design_name="${result_name#${struct_name}_}"
                # Remove _data suffix if present
                design_name="${design_name%_data}"

                local target_dir="$SCRIPT_DIR/$struct_name/03_alphafold3_predictions_${state}/$design_name"
                if [ -d "$target_dir" ] || echo "$design_name" | grep -q "_${state}"; then
                    mkdir -p "$target_dir"
                    cp -rn "$result_dir"/* "$target_dir/" 2>/dev/null
                    break 2
                fi
            done
        done
    done

    echo "[$(date '+%H:%M:%S')] GPU $node_label: Inference complete (exit code: $exit_code)"
    return $exit_code
}

if [ $SKIP_INFERENCE -eq 1 ]; then
    echo ""
    echo "PHASE 2: SKIPPED (--skip-inference / --msa-only)"
    echo ""
else

echo ""
echo "============================================================"
echo "PHASE 2: Structure Inference (2 GPUs, model loads once each)"
echo "Started: $(date)"
echo "============================================================"
echo ""

# Prepare batch directories for each GPU
BATCH_DIR_0="$SCRIPT_DIR/.af3_batch_gpu0"
BATCH_DIR_1="$SCRIPT_DIR/.af3_batch_gpu1"
rm -rf "$BATCH_DIR_0" "$BATCH_DIR_1"

echo "Preparing GPU 0 batch..."
GPU0_SUMMARY=$(prepare_inference_batch "$BATCH_DIR_0" "${NODE0_STRUCTS[@]}")
echo "  GPU 0: $GPU0_SUMMARY"

echo "Preparing GPU 1 batch..."
GPU1_SUMMARY=$(prepare_inference_batch "$BATCH_DIR_1" "${NODE1_STRUCTS[@]}")
echo "  GPU 1: $GPU1_SUMMARY"
echo ""

# Launch both GPUs in parallel
echo "Launching inference..."
echo "  Monitor: tail -f $SCRIPT_DIR/run_af3_gpu0_inference.log"
echo "           tail -f $SCRIPT_DIR/run_af3_gpu1_inference.log"
echo ""

FAILURES=0
if [ $NODE2_ONLY -eq 0 ]; then
    run_inference_batch "$GPU0" "0" "$BATCH_DIR_0" \
        > "$SCRIPT_DIR/run_af3_gpu0_inference.log" 2>&1 &
    GPU0_PID=$!
    echo "GPU 0: PID $GPU0_PID"
fi

if [ $NODE1_ONLY -eq 0 ]; then
    run_inference_batch "$GPU1" "1" "$BATCH_DIR_1" \
        > "$SCRIPT_DIR/run_af3_gpu1_inference.log" 2>&1 &
    GPU1_PID=$!
    echo "GPU 1: PID $GPU1_PID"
fi

# Wait for both
echo ""
echo "Waiting for both GPUs to finish..."
if [ -n "$GPU0_PID" ]; then
    wait $GPU0_PID || ((FAILURES++))
fi
if [ -n "$GPU1_PID" ]; then
    wait $GPU1_PID || ((FAILURES++))
fi

# Cleanup batch dirs
rm -rf "$BATCH_DIR_0" "$BATCH_DIR_1"

fi  # end SKIP_INFERENCE check

# ============================================================================
# Summary
# ============================================================================
END_TIME=$(date +%s)
ELAPSED_SEC=$(( END_TIME - START_TIME ))
ELAPSED_HRS=$(( ELAPSED_SEC / 3600 ))
ELAPSED_MIN=$(( ELAPSED_SEC / 60 ))

echo ""
echo "============================================================================"
echo "AlphaFold3 Batch Complete"
echo "============================================================================"
echo "  Elapsed: ${ELAPSED_HRS}h ${ELAPSED_MIN}m total"
echo ""

# Count results
TOTAL_CIF=0
for d in "$SCRIPT_DIR"/output_G4FP_*/; do
    CIF=$(find "$d" -path "*/03_alphafold3_predictions_*/*.cif" -type f 2>/dev/null | wc -l)
    TOTAL_CIF=$((TOTAL_CIF + CIF))
    if [ $CIF -gt 0 ]; then
        echo "  $(basename "$d"): $CIF CIF files"
    fi
done
echo ""
echo "  Total CIF structures: $TOTAL_CIF / $TOTAL_PREDS expected"
if [ $FAILURES -gt 0 ]; then
    echo "  WARNING: $FAILURES GPU(s) had errors. Check logs."
fi
echo ""
echo "Next: run analysis (steps 5-8)"
echo "  ./run_all.sh --skip-ligandmpnn --skip-af3"

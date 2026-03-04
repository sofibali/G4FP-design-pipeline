#!/bin/bash
# Run AlphaFold3 on the 10 template wild-type sequences (apo + bound = 20 predictions)
#
# Uses the same 2-phase approach as run_all_af3_parallel.sh but simplified for 20 inputs.
# Output goes to template_af3_outputs/{apo,bound}/
#
# Usage:
#   ./run_template_af3.sh                  # defaults: GPUs 2,3
#   ./run_template_af3.sh --gpu0 0 --gpu1 1
#   ./run_template_af3.sh --skip-msa       # skip phase 1 if MSA already done
#   ./run_template_af3.sh --dry-run

set +e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Defaults
GPU0=2
GPU1=3
MSA_PARALLEL=20
SKIP_MSA=0
DRY_RUN=0
JAX_CACHE_DIR="$SCRIPT_DIR/.jax_cache"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --gpu0)          GPU0="$2"; shift 2 ;;
        --gpu1)          GPU1="$2"; shift 2 ;;
        --msa-parallel)  MSA_PARALLEL="$2"; shift 2 ;;
        --skip-msa)      SKIP_MSA=1; shift ;;
        --dry-run)       DRY_RUN=1; shift ;;
        *)               echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# Environment setup
source /programs/sbgrid.shrc 2>/dev/null
export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas
mkdir -p "$JAX_CACHE_DIR"

INPUT_DIR="$SCRIPT_DIR/template_af3_inputs"
OUTPUT_DIR="$SCRIPT_DIR/template_af3_outputs"

mkdir -p "$OUTPUT_DIR/apo" "$OUTPUT_DIR/bound"

# Count inputs
N_APO=$(ls "$INPUT_DIR/apo/"*.json 2>/dev/null | wc -l)
N_BOUND=$(ls "$INPUT_DIR/bound/"*.json 2>/dev/null | wc -l)
N_TOTAL=$((N_APO + N_BOUND))

echo "============================================================================"
echo "AlphaFold3 Template Predictions (Wild-Type Sequences)"
echo "============================================================================"
echo ""
echo "  Apo inputs:    $N_APO  ($INPUT_DIR/apo/)"
echo "  Bound inputs:  $N_BOUND ($INPUT_DIR/bound/)"
echo "  Total:         $N_TOTAL"
echo "  Output:        $OUTPUT_DIR/"
echo "  GPUs:          $GPU0, $GPU1"
echo "  JAX cache:     $JAX_CACHE_DIR"
echo ""

if [ $DRY_RUN -eq 1 ]; then
    echo "[DRY RUN] Would run the above. Exiting."
    exit 0
fi

START_TIME=$(date +%s)

# ============================================================================
# PHASE 1: Data Pipeline (MSA) -- CPU only, all 20 in parallel
# ============================================================================
run_data_pipeline() {
    local json_file=$1
    local state=$2
    local design_name=$(basename "$json_file" .json)
    local out="$OUTPUT_DIR/$state/$design_name"

    # Skip if data pipeline already completed
    if [ -n "$(find "$out" -name '*_data.json' -print -quit 2>/dev/null)" ]; then
        return 0
    fi
    # Skip if inference already completed
    if [ -n "$(find "$out" -name '*.cif' -print -quit 2>/dev/null)" ]; then
        return 0
    fi

    mkdir -p "$out"

    /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
        --db_dir /mnt/alphafold3 \
        --model_dir /mnt/alphafold3 \
        --output_dir "$out" \
        --json_path "$json_file" \
        --norun_inference \
        &> "$OUTPUT_DIR/${state}/${design_name}_msa.log"
}

if [ $SKIP_MSA -eq 0 ]; then
    echo "============================================================"
    echo "PHASE 1: Data Pipeline (MSA, $MSA_PARALLEL parallel)"
    echo "Started: $(date)"
    echo "============================================================"
    echo ""

    export -f run_data_pipeline
    export OUTPUT_DIR

    SUBMITTED=0
    SKIPPED=0

    for state in apo bound; do
        for json_file in "$INPUT_DIR/$state/"*.json; do
            design_name=$(basename "$json_file" .json)
            out="$OUTPUT_DIR/$state/$design_name"

            if [ -n "$(find "$out" -name '*_data.json' -print -quit 2>/dev/null)" ] || \
               [ -n "$(find "$out" -name '*.cif' -print -quit 2>/dev/null)" ]; then
                ((SKIPPED++))
                continue
            fi

            run_data_pipeline "$json_file" "$state" &
            ((SUBMITTED++))

            # Throttle
            if [ $SUBMITTED -ge $MSA_PARALLEL ]; then
                wait -n 2>/dev/null || wait
            fi
        done
    done

    wait
    echo "[$(date '+%H:%M:%S')] Phase 1 complete: $SUBMITTED submitted, $SKIPPED skipped"
    echo ""
else
    echo "PHASE 1: SKIPPED (--skip-msa)"
    echo ""
fi

# ============================================================================
# PHASE 2: Inference -- 2 GPUs, model loads once per GPU
# ============================================================================
# GPU 0: apo (10 predictions)
# GPU 1: bound (10 predictions)
# ============================================================================

echo "============================================================"
echo "PHASE 2: Structure Inference"
echo "  GPU $GPU0: apo templates (10 predictions)"
echo "  GPU $GPU1: bound templates (10 predictions)"
echo "Started: $(date)"
echo "============================================================"
echo ""

run_inference_state() {
    local gpu_id=$1
    local state=$2
    local batch_dir="$SCRIPT_DIR/.template_batch_${state}"

    rm -rf "$batch_dir"
    mkdir -p "$batch_dir"

    local n_inputs=0
    for json_file in "$INPUT_DIR/$state/"*.json; do
        local design_name=$(basename "$json_file" .json)
        local out="$OUTPUT_DIR/$state/$design_name"

        # Skip if inference already completed
        if [ -n "$(find "$out" -name '*.cif' -print -quit 2>/dev/null)" ]; then
            echo "  Skipping $design_name (already done)"
            continue
        fi

        # Prefer _data.json if available
        local data_json=$(find "$out" -name '*_data.json' -print -quit 2>/dev/null)
        if [ -n "$data_json" ]; then
            ln -sf "$data_json" "$batch_dir/"
        else
            ln -sf "$json_file" "$batch_dir/"
        fi
        ((n_inputs++))
    done

    if [ "$n_inputs" -eq 0 ]; then
        echo "  GPU $gpu_id ($state): All predictions already complete"
        rm -rf "$batch_dir"
        return 0
    fi

    echo "  GPU $gpu_id ($state): $n_inputs predictions"

    # Check if all are _data.json
    local n_data=$(ls "$batch_dir"/*_data.json 2>/dev/null | wc -l)
    local pipeline_flag=""
    if [ "$n_data" -eq "$n_inputs" ]; then
        pipeline_flag="--norun_data_pipeline"
    fi

    export CUDA_VISIBLE_DEVICES=$gpu_id

    /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
        --db_dir /mnt/alphafold3 \
        --model_dir /mnt/alphafold3 \
        --output_dir "$OUTPUT_DIR/$state" \
        --input_dir "$batch_dir" \
        --jax_compilation_cache_dir "$JAX_CACHE_DIR" \
        --flash_attention_implementation triton \
        --num_diffusion_samples 5 \
        $pipeline_flag

    local exit_code=$?
    rm -rf "$batch_dir"

    echo "[$(date '+%H:%M:%S')] GPU $gpu_id ($state): Done (exit code: $exit_code)"
    return $exit_code
}

# Launch both GPUs in parallel
run_inference_state "$GPU0" "apo" \
    > "$SCRIPT_DIR/run_template_af3_apo.log" 2>&1 &
APO_PID=$!
echo "  Apo PID: $APO_PID (GPU $GPU0)"

run_inference_state "$GPU1" "bound" \
    > "$SCRIPT_DIR/run_template_af3_bound.log" 2>&1 &
BOUND_PID=$!
echo "  Bound PID: $BOUND_PID (GPU $GPU1)"

echo ""
echo "  Monitor: tail -f $SCRIPT_DIR/run_template_af3_apo.log"
echo "           tail -f $SCRIPT_DIR/run_template_af3_bound.log"
echo ""

# Wait
FAILURES=0
wait $APO_PID || ((FAILURES++))
wait $BOUND_PID || ((FAILURES++))

# ============================================================================
# Summary
# ============================================================================
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
ELAPSED_MIN=$((ELAPSED / 60))

echo ""
echo "============================================================================"
echo "Template AF3 Predictions Complete"
echo "============================================================================"
echo "  Elapsed: ${ELAPSED_MIN} minutes"
echo ""

# Count results
N_CIF_APO=$(find "$OUTPUT_DIR/apo" -name '*.cif' -type f 2>/dev/null | wc -l)
N_CIF_BOUND=$(find "$OUTPUT_DIR/bound" -name '*.cif' -type f 2>/dev/null | wc -l)
echo "  Apo CIF files:   $N_CIF_APO / $N_APO"
echo "  Bound CIF files: $N_CIF_BOUND / $N_BOUND"
if [ $FAILURES -gt 0 ]; then
    echo "  WARNING: $FAILURES GPU(s) had errors. Check logs."
fi
echo ""
echo "Output: $OUTPUT_DIR/"
echo ""

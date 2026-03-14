#!/bin/bash
# Re-run AF3 inference after chromophore fix.
#
# This script:
#   1. Scans for completed designs (those with *_model.cif in output)
#   2. Creates per-batch input dirs with only incomplete designs
#   3. Re-runs AF3 with --norun_data_pipeline (inference only)
#      on 2 GPUs with CPU thread limits (max 16 cores total)
#
# Prerequisites:
#   - Run inject_msa_into_inputs.py first (MSA must be in input JSONs)
#
# Usage:
#   bash rerun_af3_inference.sh                    # all templates, 2 GPUs
#   bash rerun_af3_inference.sh --template G4FP_des1_cro_mod0  # one template
#   bash rerun_af3_inference.sh --dry-run          # preview only
#   bash rerun_af3_inference.sh --test-one         # one design only

set -euo pipefail

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# ── Configuration ──
MAX_CORES=16
GPU0=0
GPU1=1

# ── Parse arguments ──
DRY_RUN=false
TEST_ONE=false
TEMPLATE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --test-one) TEST_ONE=true; shift ;;
        --template) TEMPLATE="$2"; shift 2 ;;
        --gpu0) GPU0="$2"; shift 2 ;;
        --gpu1) GPU1="$2"; shift 2 ;;
        --max-cores) MAX_CORES="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Limit CPU threads: split evenly between 2 GPU processes
CORES_PER_GPU=$((MAX_CORES / 2))
export OMP_NUM_THREADS=$CORES_PER_GPU
export OPENBLAS_NUM_THREADS=$CORES_PER_GPU
export MKL_NUM_THREADS=$CORES_PER_GPU
export NUMEXPR_MAX_THREADS=$CORES_PER_GPU
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false"

# AF3 configuration
AF3_RUNNER="/programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py"
AF3_DB="/mnt/alphafold3"
JAX_CACHE="$SCRIPT_DIR/.jax_cache"
export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}AF3 Inference (skip completed designs)${NC}"
echo -e "${GREEN}========================================${NC}"
echo "  Max CPU cores: $MAX_CORES ($CORES_PER_GPU per GPU)"
echo "  GPUs: $GPU0, $GPU1"

if [ -n "$TEMPLATE" ]; then
    OUTPUT_DIRS=("output_${TEMPLATE}")
else
    OUTPUT_DIRS=(output_G4FP_*)
fi

# ══════════════════════════════════════════════════════════════════════
# Step 1: Identify incomplete designs
# ══════════════════════════════════════════════════════════════════════
echo ""
echo -e "${BLUE}Step 1: Scanning for incomplete designs...${NC}"

# is_design_complete <pred_dir> <design_name>
# Checks if a design already has a *_model.cif in its output dir
is_design_complete() {
    local pred_dir="$1"
    local name="$2"  # e.g. design_0_bound

    # Check non-timestamped dir first (canonical)
    if [ -d "$pred_dir/$name" ]; then
        if ls "$pred_dir/$name"/*_model.cif &>/dev/null; then
            return 0  # complete
        fi
    fi

    # Check timestamped dirs
    for d in "$pred_dir"/${name}_[0-9]*_[0-9]*/; do
        [ ! -d "$d" ] && continue
        if ls "$d"/*_model.cif &>/dev/null; then
            return 0  # complete
        fi
    done

    return 1  # incomplete
}

TOTAL_COMPLETE=0
TOTAL_INCOMPLETE=0
TOTAL_DESIGNS=0

# Build list of (template, state, incomplete_count) for GPU scheduling
declare -a WORK_ITEMS=()

for odir in "${OUTPUT_DIRS[@]}"; do
    [ ! -d "$odir" ] && continue

    for state in bound apo; do
        json_dir="$odir/02_alphafold3_inputs_${state}"
        [ ! -d "$json_dir" ] && continue
        pred_dir="$odir/03_alphafold3_predictions_${state}"
        mkdir -p "$pred_dir"

        # Create temp dir for incomplete designs
        tmp_input="$SCRIPT_DIR/.incomplete_inputs/${odir}_${state}"
        rm -rf "$tmp_input"
        mkdir -p "$tmp_input"

        n_complete=0
        n_incomplete=0

        # Batch-extract names from all JSONs at once (much faster than per-file python3)
        while IFS=$'\t' read -r jf name; do
            [ -z "$name" ] && continue

            ((TOTAL_DESIGNS++)) || true

            if is_design_complete "$pred_dir" "$name"; then
                ((n_complete++)) || true
            else
                # Copy to temp input dir
                if ! $DRY_RUN; then
                    cp "$jf" "$tmp_input/"
                fi
                ((n_incomplete++)) || true

                if $TEST_ONE && [ "$n_incomplete" -ge 1 ]; then
                    break
                fi
            fi
        done < <(python3 -c "
import json, glob, os
d = '$json_dir'
for f in sorted(glob.glob(os.path.join(d, '*.json'))):
    if '.bak' in f:
        continue
    try:
        name = json.load(open(f)).get('name','')
        if name:
            print(f'{f}\t{name}')
    except: pass
")

        TOTAL_COMPLETE=$((TOTAL_COMPLETE + n_complete))
        TOTAL_INCOMPLETE=$((TOTAL_INCOMPLETE + n_incomplete))

        if [ "$n_incomplete" -gt 0 ]; then
            WORK_ITEMS+=("$odir|$state|$n_incomplete|$tmp_input|$pred_dir")
            echo "  $odir/$state: $n_complete done, $n_incomplete remaining"
        else
            echo -e "  $odir/$state: ${GREEN}all $n_complete complete${NC}"
            rm -rf "$tmp_input"
        fi

        if $TEST_ONE && [ "$TOTAL_INCOMPLETE" -ge 1 ]; then
            break 2
        fi
    done
done

echo ""
echo -e "${GREEN}  Total: $TOTAL_COMPLETE complete, $TOTAL_INCOMPLETE remaining (of $TOTAL_DESIGNS)${NC}"

if [ "$TOTAL_INCOMPLETE" -eq 0 ]; then
    echo -e "${GREEN}All designs already complete!${NC}"
    rm -rf "$SCRIPT_DIR/.incomplete_inputs"
    exit 0
fi

if $DRY_RUN; then
    echo ""
    echo -e "${YELLOW}DRY RUN complete. Re-run without --dry-run to proceed.${NC}"
    rm -rf "$SCRIPT_DIR/.incomplete_inputs"
    exit 0
fi

# ══════════════════════════════════════════════════════════════════════
# Step 2: Run AF3 inference (2 GPUs, inference only)
# ══════════════════════════════════════════════════════════════════════
echo ""
echo -e "${BLUE}Step 2: Running AF3 inference (inference only, 2 GPUs)...${NC}"
echo "  CPU threads per GPU: $CORES_PER_GPU"
echo "  JAX cache: $JAX_CACHE"

source /programs/sbgrid.shrc 2>/dev/null || true
mkdir -p "$JAX_CACHE"

if $TEST_ONE; then
    # Single design test — use first work item
    IFS='|' read -r odir state n_inc tmp_input pred_dir <<< "${WORK_ITEMS[0]}"

    echo -e "${BLUE}  TEST: $(ls "$tmp_input"/*.json | head -1 | xargs basename)${NC}"
    CUDA_VISIBLE_DEVICES=$GPU0 "$AF3_RUNNER" \
        --input_dir="$tmp_input" \
        --output_dir="$pred_dir" \
        --model_dir="$AF3_DB" \
        --db_dir="$AF3_DB" \
        --norun_data_pipeline \
        --jax_compilation_cache_dir="$JAX_CACHE" \
        --flash_attention_implementation=triton \
        --buckets=256,512,768,1024,1280

    rm -rf "$SCRIPT_DIR/.incomplete_inputs"
    echo -e "${GREEN}  Test complete!${NC}"
    exit 0
fi

# ── Full run: split work items between 2 GPUs ──

N_ITEMS=${#WORK_ITEMS[@]}
HALF=$(( (N_ITEMS + 1) / 2 ))

GPU0_ITEMS=("${WORK_ITEMS[@]:0:$HALF}")
GPU1_ITEMS=()
if [ "$N_ITEMS" -gt "$HALF" ]; then
    GPU1_ITEMS=("${WORK_ITEMS[@]:$HALF}")
fi

echo "  GPU $GPU0 batches: ${#GPU0_ITEMS[@]}"
echo "  GPU $GPU1 batches: ${#GPU1_ITEMS[@]}"

# Function to run inference for a set of work items on one GPU
run_gpu_inference() {
    local gpu_id=$1
    local gpu_label=$2
    shift 2
    local items=("$@")

    # Per-process thread limits
    export OMP_NUM_THREADS=$CORES_PER_GPU
    export OPENBLAS_NUM_THREADS=$CORES_PER_GPU
    export MKL_NUM_THREADS=$CORES_PER_GPU
    export NUMEXPR_MAX_THREADS=$CORES_PER_GPU
    export CUDA_VISIBLE_DEVICES=$gpu_id

    for item in "${items[@]}"; do
        IFS='|' read -r odir state n_inc tmp_input pred_dir <<< "$item"
        local tpl=$(echo "$odir" | sed 's/output_//')

        echo "[$(date '+%H:%M:%S')] GPU $gpu_label: $tpl / $state ($n_inc designs)"

        "$AF3_RUNNER" \
            --input_dir="$tmp_input" \
            --output_dir="$pred_dir" \
            --model_dir="$AF3_DB" \
            --db_dir="$AF3_DB" \
            --norun_data_pipeline \
            --jax_compilation_cache_dir="$JAX_CACHE" \
            --flash_attention_implementation=triton \
            --buckets=256,512,768,1024,1280

        echo "[$(date '+%H:%M:%S')] GPU $gpu_label: done $tpl / $state"
    done
}

# Launch both GPUs in parallel
echo ""
echo -e "${BLUE}Launching inference on 2 GPUs...${NC}"
echo "  Monitor: tail -f rerun_gpu0.log"
echo "           tail -f rerun_gpu1.log"
echo ""

run_gpu_inference "$GPU0" "0" "${GPU0_ITEMS[@]}" \
    > "$SCRIPT_DIR/rerun_gpu0.log" 2>&1 &
PID0=$!
echo "  GPU 0: PID $PID0"

if [ ${#GPU1_ITEMS[@]} -gt 0 ]; then
    run_gpu_inference "$GPU1" "1" "${GPU1_ITEMS[@]}" \
        > "$SCRIPT_DIR/rerun_gpu1.log" 2>&1 &
    PID1=$!
    echo "  GPU 1: PID $PID1"
fi

# Wait for both
echo ""
echo "Waiting for inference to complete..."
FAILURES=0

wait $PID0 || ((FAILURES++))
echo -e "${GREEN}  GPU 0 finished${NC}"

if [ ${#GPU1_ITEMS[@]} -gt 0 ]; then
    wait $PID1 || ((FAILURES++))
    echo -e "${GREEN}  GPU 1 finished${NC}"
fi

# Cleanup temp dirs
rm -rf "$SCRIPT_DIR/.incomplete_inputs"

echo ""
echo -e "${GREEN}========================================${NC}"
if [ $FAILURES -eq 0 ]; then
    echo -e "${GREEN}AF3 inference complete!${NC}"
else
    echo -e "${RED}AF3 inference finished with $FAILURES GPU failure(s)${NC}"
    echo -e "${RED}Check rerun_gpu0.log and rerun_gpu1.log${NC}"
fi
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Next steps:"
echo "  cd $SCRIPT_DIR"
echo "  python3 05_analyze_af3_structures.py --force"
echo "  python3 06_compare_ligand_states.py --force"
echo "  python3 07_aggregate_and_rank.py --force"

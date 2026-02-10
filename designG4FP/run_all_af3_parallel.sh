#!/bin/bash
# Run AlphaFold3 for all G4FP structures across two GPU nodes
#
# Splits the 10 structures across 2 nodes (5 per node), running 4 AF3 jobs
# in parallel on each node. Each job predicts one design (bound or apo).
#
# Usage:
#   ./run_all_af3_parallel.sh                       # defaults: 2 nodes, 4 jobs each
#   ./run_all_af3_parallel.sh --jobs-per-node 4     # explicit
#   ./run_all_af3_parallel.sh --gpu0 0 --gpu1 1     # specify GPU IDs
#   ./run_all_af3_parallel.sh --node1-only           # run only node 0's share
#   ./run_all_af3_parallel.sh --node2-only           # run only node 1's share
#   ./run_all_af3_parallel.sh --dry-run              # show plan without running
#
# Time estimate (printed at startup):
#   2,000 predictions / (2 nodes x 4 parallel) = 250 serial batches
#   At ~20-40 min per prediction => ~83-167 hours wall time per node
#   With 4 parallel jobs => ~21-42 hours per node
#
# For monitoring:
#   tail -f run_af3_node*.log
#   find output_G4FP_*/03_alphafold3_predictions_* -name '*.cif' | wc -l

set +e  # Continue even if individual predictions fail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_TEMPLATE="$SCRIPT_DIR/pipeline_config_template.yaml"

# ============================================================================
# Defaults
# ============================================================================
JOBS_PER_NODE=4
GPU0=0
GPU1=1
RUN_MODE="both"      # bound, apo, both
NODE1_ONLY=0
NODE2_ONLY=0
DRY_RUN=0

# ============================================================================
# Parse arguments
# ============================================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        --jobs-per-node) JOBS_PER_NODE="$2"; shift 2 ;;
        --gpu0)          GPU0="$2"; shift 2 ;;
        --gpu1)          GPU1="$2"; shift 2 ;;
        --mode)          RUN_MODE="$2"; shift 2 ;;
        --node1-only)    NODE1_ONLY=1; shift ;;
        --node2-only)    NODE2_ONLY=1; shift ;;
        --dry-run)       DRY_RUN=1; shift ;;
        *)               echo "Unknown arg: $1"; exit 1 ;;
    esac
done

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
# Time estimate
# ============================================================================
# AF3 per prediction: ~20 min (apo/small) to ~40 min (bound/large)
# Average ~30 min per prediction
AVG_MIN=30
NODE0_SERIAL_HRS=$(( NODE0_PREDS * AVG_MIN / 60 ))
NODE1_SERIAL_HRS=$(( NODE1_PREDS * AVG_MIN / 60 ))
NODE0_WALL_HRS=$(( NODE0_SERIAL_HRS / JOBS_PER_NODE ))
NODE1_WALL_HRS=$(( NODE1_SERIAL_HRS / JOBS_PER_NODE ))
TOTAL_WALL_HRS=$(( NODE0_WALL_HRS > NODE1_WALL_HRS ? NODE0_WALL_HRS : NODE1_WALL_HRS ))

echo "============================================================================"
echo "AlphaFold3 Batch Predictions -- 2 GPU Nodes"
echo "============================================================================"
echo ""
echo "Structures found:      $N_STRUCTS"
echo "Run mode:              $RUN_MODE"
echo "Total predictions:     $TOTAL_PREDS"
echo "Jobs per node:         $JOBS_PER_NODE"
echo ""
echo "--- Node 0 (GPU $GPU0) ---"
echo "  Structures: ${NODE0_STRUCTS[*]}" | sed 's/output_G4FP_//g'
echo "  Predictions: $NODE0_PREDS"
echo "  Estimated wall time: ~${NODE0_WALL_HRS} hours (${NODE0_SERIAL_HRS}h serial / $JOBS_PER_NODE parallel)"
echo ""
echo "--- Node 1 (GPU $GPU1) ---"
echo "  Structures: ${NODE1_STRUCTS[*]}" | sed 's/output_G4FP_//g'
echo "  Predictions: $NODE1_PREDS"
echo "  Estimated wall time: ~${NODE1_WALL_HRS} hours (${NODE1_SERIAL_HRS}h serial / $JOBS_PER_NODE parallel)"
echo ""
echo "============================================================================"
echo "ESTIMATED TOTAL WALL TIME: ~${TOTAL_WALL_HRS} hours"
echo "  (both nodes running simultaneously, limited by the slower node)"
echo "  Range: ~$((TOTAL_WALL_HRS * 2 / 3))h (fast, 20 min/pred) to ~$((TOTAL_WALL_HRS * 4 / 3))h (slow, 40 min/pred)"
echo "============================================================================"
echo ""

if [ $DRY_RUN -eq 1 ]; then
    echo "[DRY RUN] Would run the above. Exiting."
    exit 0
fi

# ============================================================================
# AF3 runner function (one structure, one state)
# ============================================================================
run_structure_state() {
    local struct_name=$1
    local state=$2       # bound or apo
    local gpu_id=$3
    local n_parallel=$4
    local struct_dir="$SCRIPT_DIR/$struct_name"
    local input_dir="$struct_dir/02_alphafold3_inputs_${state}"
    local output_dir="$struct_dir/03_alphafold3_predictions_${state}"

    if [ ! -d "$input_dir" ]; then
        return 0
    fi

    mkdir -p "$output_dir"

    local json_files=("$input_dir"/*.json)
    local n_total=${#json_files[@]}
    local n_done=0
    local n_skip=0
    local n_fail=0

    echo "[$(date '+%H:%M:%S')] Starting $struct_name $state ($n_total predictions, GPU $gpu_id, $n_parallel parallel)"

    # Process in batches of n_parallel
    local batch=()
    for json_file in "${json_files[@]}"; do
        local design_name=$(basename "$json_file" .json)
        local job_output="$output_dir/$design_name"

        # Skip completed
        if ls "$job_output"/*_model*.cif &>/dev/null 2>&1 || ls "$job_output"/model_cif/model_0.cif &>/dev/null 2>&1; then
            ((n_skip++))
            continue
        fi

        batch+=("$json_file")

        # When batch is full, run in parallel
        if [ ${#batch[@]} -ge $n_parallel ]; then
            for jf in "${batch[@]}"; do
                local dn=$(basename "$jf" .json)
                local jo="$output_dir/$dn"
                mkdir -p "$jo"
                (
                    source /programs/sbgrid.shrc 2>/dev/null
                    export CUDA_VISIBLE_DEVICES=$gpu_id
                    export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas
                    /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
                        --db_dir /mnt/alphafold3 \
                        --model_dir /mnt/alphafold3 \
                        --output_dir "$jo" \
                        --json_path "$jf" &> "$jo.log"
                ) &
            done
            wait
            n_done=$((n_done + ${#batch[@]}))
            echo "[$(date '+%H:%M:%S')]   $struct_name $state: $n_done/$n_total done ($n_skip skipped)"
            batch=()
        fi
    done

    # Run remaining partial batch
    if [ ${#batch[@]} -gt 0 ]; then
        for jf in "${batch[@]}"; do
            local dn=$(basename "$jf" .json)
            local jo="$output_dir/$dn"
            mkdir -p "$jo"
            (
                source /programs/sbgrid.shrc 2>/dev/null
                export CUDA_VISIBLE_DEVICES=$gpu_id
                export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas
                /programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py \
                    --db_dir /mnt/alphafold3 \
                    --model_dir /mnt/alphafold3 \
                    --output_dir "$jo" \
                    --json_path "$jf" &> "$jo.log"
            ) &
        done
        wait
        n_done=$((n_done + ${#batch[@]}))
    fi

    # Count actual successes
    local n_cif=$(find "$output_dir" -name "*.cif" -type f 2>/dev/null | wc -l)
    echo "[$(date '+%H:%M:%S')] Finished $struct_name $state: $n_cif CIF files ($n_skip skipped)"
}

# ============================================================================
# Node runner: processes a list of structures on one GPU
# ============================================================================
run_node() {
    local gpu_id=$1
    local node_label=$2
    shift 2
    local structs=("$@")

    echo ""
    echo "===== NODE $node_label (GPU $gpu_id) starting at $(date) ====="
    echo ""

    for struct_name in "${structs[@]}"; do
        if [ "$RUN_MODE" = "bound" ] || [ "$RUN_MODE" = "both" ]; then
            run_structure_state "$struct_name" "bound" "$gpu_id" "$JOBS_PER_NODE"
        fi
        if [ "$RUN_MODE" = "apo" ] || [ "$RUN_MODE" = "both" ]; then
            run_structure_state "$struct_name" "apo" "$gpu_id" "$JOBS_PER_NODE"
        fi
    done

    echo ""
    echo "===== NODE $node_label (GPU $gpu_id) finished at $(date) ====="
}

# ============================================================================
# Launch both nodes in parallel (as background processes)
# ============================================================================
echo "Launching AF3 predictions..."
echo "  Monitor with: tail -f $SCRIPT_DIR/run_af3_node0.log"
echo "                tail -f $SCRIPT_DIR/run_af3_node1.log"
echo ""

START_TIME=$(date +%s)

if [ $NODE2_ONLY -eq 0 ]; then
    run_node "$GPU0" "0" "${NODE0_STRUCTS[@]}" > "$SCRIPT_DIR/run_af3_node0.log" 2>&1 &
    NODE0_PID=$!
    echo "Node 0 (GPU $GPU0): PID $NODE0_PID -> run_af3_node0.log"
fi

if [ $NODE1_ONLY -eq 0 ]; then
    run_node "$GPU1" "1" "${NODE1_STRUCTS[@]}" > "$SCRIPT_DIR/run_af3_node1.log" 2>&1 &
    NODE1_PID=$!
    echo "Node 1 (GPU $GPU1): PID $NODE1_PID -> run_af3_node1.log"
fi

# Wait for both
echo ""
echo "Waiting for both nodes to finish..."
FAILURES=0
if [ -n "$NODE0_PID" ]; then
    wait $NODE0_PID || ((FAILURES++))
fi
if [ -n "$NODE1_PID" ]; then
    wait $NODE1_PID || ((FAILURES++))
fi

END_TIME=$(date +%s)
ELAPSED=$(( (END_TIME - START_TIME) / 3600 ))
ELAPSED_MIN=$(( (END_TIME - START_TIME) / 60 ))

echo ""
echo "============================================================================"
echo "AlphaFold3 Batch Complete"
echo "============================================================================"
echo "  Elapsed: ${ELAPSED}h ${ELAPSED_MIN}m"
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
echo ""
echo "Next: run analysis (steps 5-8)"
echo "  ./run_all.sh --skip-ligandmpnn --skip-af3"

#!/bin/bash
# Re-run AF3 inference after chromophore fix.
#
# This script:
#   1. Archives current results (predictions + analysis CSVs)
#   2. Clears old inference outputs (CIFs, confidences, rankings)
#      but preserves the fixed _data.json files
#   3. Re-runs AF3 with --norun_data_pipeline (inference only)
#      on 2 GPUs with CPU thread limits (max 16 cores total)
#
# Usage:
#   bash rerun_af3_inference.sh                    # all templates, 2 GPUs
#   bash rerun_af3_inference.sh --template G4FP_des1_cro_mod0  # one template
#   bash rerun_af3_inference.sh --dry-run          # preview only
#   bash rerun_af3_inference.sh --test-one         # one design only
#   bash rerun_af3_inference.sh --skip-archive     # don't archive old results
#   bash rerun_af3_inference.sh --archive-only     # archive without re-running

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
SKIP_ARCHIVE=false
ARCHIVE_ONLY=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --test-one) TEST_ONE=true; shift ;;
        --template) TEMPLATE="$2"; shift 2 ;;
        --gpu0) GPU0="$2"; shift 2 ;;
        --gpu1) GPU1="$2"; shift 2 ;;
        --max-cores) MAX_CORES="$2"; shift 2 ;;
        --skip-archive) SKIP_ARCHIVE=true; shift ;;
        --archive-only) ARCHIVE_ONLY=true; shift ;;
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
echo -e "${GREEN}AF3 Inference Re-run (chromophore fix)${NC}"
echo -e "${GREEN}========================================${NC}"
echo "  Max CPU cores: $MAX_CORES ($CORES_PER_GPU per GPU)"
echo "  GPUs: $GPU0, $GPU1"

if [ -n "$TEMPLATE" ]; then
    OUTPUT_DIRS=("output_${TEMPLATE}")
else
    OUTPUT_DIRS=(output_G4FP_*)
fi

# ══════════════════════════════════════════════════════════════════════
# Step 0: Archive current results
# ══════════════════════════════════════════════════════════════════════
if ! $SKIP_ARCHIVE; then
    ARCHIVE_DIR="$SCRIPT_DIR/archive_pre_chromophore_fix_$(date +%Y%m%d_%H%M%S)"
    echo ""
    echo -e "${BLUE}Step 0: Archiving current results to ${ARCHIVE_DIR}${NC}"

    if ! $DRY_RUN; then
        mkdir -p "$ARCHIVE_DIR"

        # Archive analysis output
        if [ -d "analysis_output" ]; then
            cp -r analysis_output "$ARCHIVE_DIR/"
            echo "  Archived: analysis_output/"
        fi

        # Archive 07_results
        if [ -d "07_results" ]; then
            cp -r 07_results "$ARCHIVE_DIR/"
            echo "  Archived: 07_results/"
        fi

        # Archive per-template result CSVs and analysis dirs
        for odir in "${OUTPUT_DIRS[@]}"; do
            [ ! -d "$odir" ] && continue
            tpl_archive="$ARCHIVE_DIR/$odir"
            mkdir -p "$tpl_archive"

            for f in "$odir"/05_structure_analysis_results.csv \
                     "$odir"/06_ligand_state_comparison_results.csv; do
                [ -f "$f" ] && cp "$f" "$tpl_archive/"
            done

            for d in "$odir"/05_structure_analysis \
                     "$odir"/06_ligand_state_comparison; do
                [ -d "$d" ] && cp -r "$d" "$tpl_archive/"
            done
        done

        echo -e "${GREEN}  Archive complete: $ARCHIVE_DIR${NC}"
    else
        echo -e "${YELLOW}  (DRY RUN — would archive to $ARCHIVE_DIR)${NC}"
    fi
fi

if $ARCHIVE_ONLY; then
    echo ""
    echo -e "${GREEN}Archive-only mode complete.${NC}"
    exit 0
fi

# ══════════════════════════════════════════════════════════════════════
# Step 1: Clear old inference outputs
# ══════════════════════════════════════════════════════════════════════
echo ""
echo -e "${BLUE}Step 1: Clearing old inference outputs...${NC}"

TOTAL_CLEARED=0
for odir in "${OUTPUT_DIRS[@]}"; do
    [ ! -d "$odir" ] && continue

    for state in bound apo; do
        pred_dir="$odir/03_alphafold3_predictions_${state}"
        [ ! -d "$pred_dir" ] && continue

        n_cleared=0
        for design_dir in "$pred_dir"/design_*/; do
            [ ! -d "$design_dir" ] && continue

            for inner_dir in "$design_dir"/*/; do
                [ ! -d "$inner_dir" ] && continue

                # Check _data.json exists
                data_json=$(ls "$inner_dir"/*_data.json 2>/dev/null | head -1)
                [ -z "$data_json" ] && continue

                if ! $DRY_RUN; then
                    rm -rf "$inner_dir"/seed-*_sample-*/
                    rm -rf "$inner_dir"/seed-*/
                    rm -f "$inner_dir"/*_model.cif
                    rm -f "$inner_dir"/*_confidences.json
                    rm -f "$inner_dir"/*_summary_confidences.json
                    rm -f "$inner_dir"/ranking_scores.csv
                fi
                ((n_cleared++)) || true

                if $TEST_ONE && [ "$n_cleared" -ge 1 ]; then
                    break 2
                fi
            done

            if $TEST_ONE && [ "$n_cleared" -ge 1 ]; then
                break
            fi
        done

        echo "  $odir/$state: $n_cleared designs cleared"
        TOTAL_CLEARED=$((TOTAL_CLEARED + n_cleared))

        if $TEST_ONE && [ "$TOTAL_CLEARED" -ge 1 ]; then
            break
        fi
    done

    if $TEST_ONE && [ "$TOTAL_CLEARED" -ge 1 ]; then
        break
    fi
done

echo -e "${GREEN}  Total cleared: $TOTAL_CLEARED designs${NC}"

# Clear per-template analysis CSVs
for odir in "${OUTPUT_DIRS[@]}"; do
    [ ! -d "$odir" ] && continue
    if ! $DRY_RUN; then
        rm -f "$odir"/05_structure_analysis_results.csv
        rm -f "$odir"/06_ligand_state_comparison_results.csv
    fi
done

if $DRY_RUN; then
    echo ""
    echo -e "${YELLOW}DRY RUN complete. Re-run without --dry-run to proceed.${NC}"
    exit 0
fi

# ══════════════════════════════════════════════════════════════════════
# Step 2: Re-run AF3 inference (2 GPUs, inference only)
# ══════════════════════════════════════════════════════════════════════
echo ""
echo -e "${BLUE}Step 2: Re-running AF3 inference (inference only, 2 GPUs)...${NC}"
echo "  CPU threads per GPU: $CORES_PER_GPU"
echo "  JAX cache: $JAX_CACHE"

source /programs/sbgrid.shrc 2>/dev/null || true
mkdir -p "$JAX_CACHE"

if $TEST_ONE; then
    # Single design test
    FIRST_JSON=""
    FIRST_PRED=""
    for odir in "${OUTPUT_DIRS[@]}"; do
        for state in bound apo; do
            json_dir="$odir/02_alphafold3_inputs_${state}"
            [ ! -d "$json_dir" ] && continue
            for jf in "$json_dir"/*.json; do
                [[ "$jf" == *.bak ]] && continue
                [ -f "$jf" ] || continue
                FIRST_JSON="$jf"
                FIRST_PRED="$odir/03_alphafold3_predictions_${state}"
                break 3
            done
        done
    done

    if [ -z "$FIRST_JSON" ]; then
        echo -e "${RED}  No inputs found!${NC}"
        exit 1
    fi

    echo -e "${BLUE}  TEST: $(basename $FIRST_JSON)${NC}"
    CUDA_VISIBLE_DEVICES=$GPU0 python3 "$AF3_RUNNER" \
        --json_path="$FIRST_JSON" \
        --output_dir="$FIRST_PRED" \
        --model_dir="$AF3_DB" \
        --db_dir="$AF3_DB" \
        --norun_data_pipeline \
        --jax_compilation_cache_dir="$JAX_CACHE" \
        --flash_attention_implementation=triton \
        --buckets=256,512,768,1024,1280

    echo -e "${GREEN}  Test complete!${NC}"
    exit 0
fi

# ── Full run: split templates between 2 GPUs ──

# Collect all template output dirs
ALL_TEMPLATES=()
for odir in "${OUTPUT_DIRS[@]}"; do
    [ -d "$odir" ] && ALL_TEMPLATES+=("$odir")
done

N_TEMPLATES=${#ALL_TEMPLATES[@]}
HALF=$(( (N_TEMPLATES + 1) / 2 ))

GPU0_TEMPLATES=("${ALL_TEMPLATES[@]:0:$HALF}")
GPU1_TEMPLATES=("${ALL_TEMPLATES[@]:$HALF}")

echo "  GPU $GPU0 templates (${#GPU0_TEMPLATES[@]}): ${GPU0_TEMPLATES[*]}"
echo "  GPU $GPU1 templates (${#GPU1_TEMPLATES[@]}): ${GPU1_TEMPLATES[*]}"

# Function to run inference for a set of templates on one GPU
run_gpu_inference() {
    local gpu_id=$1
    local gpu_label=$2
    shift 2
    local templates=("$@")

    # Per-process thread limits
    export OMP_NUM_THREADS=$CORES_PER_GPU
    export OPENBLAS_NUM_THREADS=$CORES_PER_GPU
    export MKL_NUM_THREADS=$CORES_PER_GPU
    export NUMEXPR_MAX_THREADS=$CORES_PER_GPU
    export CUDA_VISIBLE_DEVICES=$gpu_id

    for odir in "${templates[@]}"; do
        local tpl=$(echo "$odir" | sed 's/output_//')

        for state in bound apo; do
            local json_dir="$odir/02_alphafold3_inputs_${state}"
            [ ! -d "$json_dir" ] && continue
            local pred_dir="$odir/03_alphafold3_predictions_${state}"
            mkdir -p "$pred_dir"

            local n_jsons=$(ls "$json_dir"/*.json 2>/dev/null | grep -vc '.bak$' || true)
            [ "$n_jsons" -eq 0 ] && continue

            echo "[$(date '+%H:%M:%S')] GPU $gpu_label: $tpl / $state ($n_jsons designs)"

            python3 "$AF3_RUNNER" \
                --input_dir="$json_dir" \
                --output_dir="$pred_dir" \
                --model_dir="$AF3_DB" \
                --db_dir="$AF3_DB" \
                --norun_data_pipeline \
                --jax_compilation_cache_dir="$JAX_CACHE" \
                --flash_attention_implementation=triton \
                --buckets=256,512,768,1024,1280

            echo "[$(date '+%H:%M:%S')] GPU $gpu_label: ✓ $tpl / $state complete"
        done
    done
}

# Launch both GPUs in parallel
echo ""
echo -e "${BLUE}Launching inference on 2 GPUs...${NC}"
echo "  Monitor: tail -f rerun_gpu0.log"
echo "           tail -f rerun_gpu1.log"
echo ""

run_gpu_inference "$GPU0" "0" "${GPU0_TEMPLATES[@]}" \
    > "$SCRIPT_DIR/rerun_gpu0.log" 2>&1 &
PID0=$!
echo "  GPU 0: PID $PID0"

if [ ${#GPU1_TEMPLATES[@]} -gt 0 ]; then
    run_gpu_inference "$GPU1" "1" "${GPU1_TEMPLATES[@]}" \
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

if [ ${#GPU1_TEMPLATES[@]} -gt 0 ]; then
    wait $PID1 || ((FAILURES++))
    echo -e "${GREEN}  GPU 1 finished${NC}"
fi

echo ""
echo -e "${GREEN}========================================${NC}"
if [ $FAILURES -eq 0 ]; then
    echo -e "${GREEN}AF3 inference re-run complete!${NC}"
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

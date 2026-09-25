#!/bin/bash
# Re-run template AF3 predictions with fixed chromophore (TYG).
#
# Runs two parallel jobs on separate GPUs:
#   GPU 0: Inference-only with patched MSA (fast)
#   GPU 1: Full pipeline with new MSA (slow but independent)
#
# Prerequisites:
#   python3 fix_template_chromophore.py   (already run)
#
# Usage:
#   bash rerun_template_af3.sh                  # both paths, 2 GPUs
#   bash rerun_template_af3.sh --patched-only   # inference-only on GPU 0
#   bash rerun_template_af3.sh --full-only      # full pipeline on GPU 1

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Configuration ──
AF3_RUNNER="/programs/x86_64-linux/system/sbgrid_bin/run_alphafold.py"
AF3_DB="/mnt/alphafold3"
JAX_CACHE="$SCRIPT_DIR/.jax_cache"
MAX_CORES=16

GPU0=0
GPU1=1
RUN_PATCHED=true
RUN_FULL=true

# ── Parse arguments ──
while [[ $# -gt 0 ]]; do
    case "$1" in
        --patched-only) RUN_FULL=false; shift ;;
        --full-only) RUN_PATCHED=false; shift ;;
        --gpu0) GPU0="$2"; shift 2 ;;
        --gpu1) GPU1="$2"; shift 2 ;;
        --max-cores) MAX_CORES="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# ── Environment ──
source /programs/sbgrid.shrc 2>/dev/null || true
export TRITON_PTXAS_PATH=/usr/local/cuda-12.4/bin/ptxas

CORES_PER_GPU=$((MAX_CORES / 2))
export OMP_NUM_THREADS=$CORES_PER_GPU
export OPENBLAS_NUM_THREADS=$CORES_PER_GPU
export MKL_NUM_THREADS=$CORES_PER_GPU
export NUMEXPR_MAX_THREADS=$CORES_PER_GPU
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false"

mkdir -p "$JAX_CACHE"

echo "========================================"
echo "Re-run Template AF3 (chromophore fix)"
echo "========================================"
echo "  Max CPU cores: $MAX_CORES ($CORES_PER_GPU per GPU)"
echo "  GPUs: $GPU0 (patched), $GPU1 (full)"
echo "  Run patched: $RUN_PATCHED"
echo "  Run full:    $RUN_FULL"
echo ""

# ── Verify inputs exist ──
if [ ! -d "template_af3_inputs/bound" ] || [ ! -d "template_af3_inputs/apo" ]; then
    echo "ERROR: template_af3_inputs/{bound,apo}/ not found"
    echo "  Run: python3 fix_template_chromophore.py"
    exit 1
fi

N_BOUND=$(ls template_af3_inputs/bound/*.json 2>/dev/null | grep -cv '.bak$' || true)
N_APO=$(ls template_af3_inputs/apo/*.json 2>/dev/null | grep -cv '.bak$' || true)
echo "  Input JSONs: $N_BOUND bound, $N_APO apo"
echo ""

# ── GPU 0: Inference-only (patched MSA) ──
if $RUN_PATCHED; then
    echo "Starting patched MSA run (GPU $GPU0)..."
    (
        export CUDA_VISIBLE_DEVICES=$GPU0
        export OMP_NUM_THREADS=$CORES_PER_GPU
        export OPENBLAS_NUM_THREADS=$CORES_PER_GPU
        export MKL_NUM_THREADS=$CORES_PER_GPU
        export NUMEXPR_MAX_THREADS=$CORES_PER_GPU

        echo "[$(date '+%H:%M:%S')] Patched: bound templates ($N_BOUND)"
        "$AF3_RUNNER" \
            --input_dir="$SCRIPT_DIR/template_af3_inputs/bound" \
            --output_dir="$SCRIPT_DIR/template_af3_outputs_patched/bound" \
            --model_dir="$AF3_DB" \
            --db_dir="$AF3_DB" \
            --norun_data_pipeline \
            --jax_compilation_cache_dir="$JAX_CACHE" \
            --flash_attention_implementation=triton \
            --buckets=256,512,768,1024,1280

        echo "[$(date '+%H:%M:%S')] Patched: apo templates ($N_APO)"
        "$AF3_RUNNER" \
            --input_dir="$SCRIPT_DIR/template_af3_inputs/apo" \
            --output_dir="$SCRIPT_DIR/template_af3_outputs_patched/apo" \
            --model_dir="$AF3_DB" \
            --db_dir="$AF3_DB" \
            --norun_data_pipeline \
            --jax_compilation_cache_dir="$JAX_CACHE" \
            --flash_attention_implementation=triton \
            --buckets=256,512,768,1024,1280

        echo "[$(date '+%H:%M:%S')] Patched: DONE"
    ) > "$SCRIPT_DIR/rerun_tpl_patched.log" 2>&1 &
    PID_PATCHED=$!
    echo "  PID: $PID_PATCHED  Log: rerun_tpl_patched.log"
fi

# ── GPU 1: Full pipeline (new MSA + inference) ──
if $RUN_FULL; then
    echo "Starting full pipeline run (GPU $GPU1)..."
    (
        export CUDA_VISIBLE_DEVICES=$GPU1
        export OMP_NUM_THREADS=$CORES_PER_GPU
        export OPENBLAS_NUM_THREADS=$CORES_PER_GPU
        export MKL_NUM_THREADS=$CORES_PER_GPU
        export NUMEXPR_MAX_THREADS=$CORES_PER_GPU

        echo "[$(date '+%H:%M:%S')] Full: bound templates ($N_BOUND)"
        "$AF3_RUNNER" \
            --input_dir="$SCRIPT_DIR/template_af3_inputs/bound" \
            --output_dir="$SCRIPT_DIR/template_af3_outputs_full/bound" \
            --model_dir="$AF3_DB" \
            --db_dir="$AF3_DB" \
            --jax_compilation_cache_dir="$JAX_CACHE" \
            --flash_attention_implementation=triton \
            --buckets=256,512,768,1024,1280

        echo "[$(date '+%H:%M:%S')] Full: apo templates ($N_APO)"
        "$AF3_RUNNER" \
            --input_dir="$SCRIPT_DIR/template_af3_inputs/apo" \
            --output_dir="$SCRIPT_DIR/template_af3_outputs_full/apo" \
            --model_dir="$AF3_DB" \
            --db_dir="$AF3_DB" \
            --jax_compilation_cache_dir="$JAX_CACHE" \
            --flash_attention_implementation=triton \
            --buckets=256,512,768,1024,1280

        echo "[$(date '+%H:%M:%S')] Full: DONE"
    ) > "$SCRIPT_DIR/rerun_tpl_full.log" 2>&1 &
    PID_FULL=$!
    echo "  PID: $PID_FULL  Log: rerun_tpl_full.log"
fi

# ── Wait ──
echo ""
echo "Monitor:"
$RUN_PATCHED && echo "  tail -f rerun_tpl_patched.log"
$RUN_FULL && echo "  tail -f rerun_tpl_full.log"
echo ""
echo "Waiting for completion..."

FAILURES=0
$RUN_PATCHED && { wait $PID_PATCHED || ((FAILURES++)); echo "  Patched: finished"; }
$RUN_FULL && { wait $PID_FULL || ((FAILURES++)); echo "  Full: finished"; }

# ── Analyze results ──
echo ""
echo "========================================"
if [ $FAILURES -eq 0 ]; then
    echo "All template AF3 runs complete!"
else
    echo "WARNING: $FAILURES run(s) failed — check logs"
fi
echo "========================================"

echo ""
echo "Next steps:"
echo "  # Analyze patched results"
echo "  python3 analyze_template_af3.py --input-dir template_af3_outputs_patched -o template_af3_metrics_patched.csv"
echo ""
echo "  # Analyze full results"
echo "  python3 analyze_template_af3.py --input-dir template_af3_outputs_full -o template_af3_metrics_full.csv"
echo ""
echo "  # Compare"
echo "  python3 -c \"import pandas as pd; p=pd.read_csv('template_af3_metrics_patched.csv'); f=pd.read_csv('template_af3_metrics_full.csv'); m=p.merge(f,on='template_name',suffixes=('_patched','_full')); print(m[['template_name','fitness_score_patched','fitness_score_full']].to_string())\""

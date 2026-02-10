#!/bin/bash
# run_all.sh -- Non-interactive batch pipeline for all G4FP input structures
#
# Runs the complete pipeline (LigandMPNN -> Filter -> AF3 -> Analysis -> Rank)
# across all 10 template structures without prompts, suitable for nohup/HPC.
#
# Usage:
#   nohup ./run_all.sh 2>&1 | tee run_all.log &
#   ./run_all.sh --parallel 4            # 4 parallel AF3 jobs per GPU node (2 nodes)
#   ./run_all.sh --skip-ligandmpnn       # skip step 1 if already done
#   ./run_all.sh --skip-af3              # skip step 3 (AF3 predictions)
#   ./run_all.sh --pareto                # use Pareto filtering in step 2

set +e  # Continue even if individual steps fail

# Parse arguments
PARALLEL=3
SKIP_LMPNN=0
SKIP_AF3=0
PARETO_FLAG=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --parallel)   PARALLEL="$2"; shift 2 ;;
        --skip-ligandmpnn) SKIP_LMPNN=1; shift ;;
        --skip-af3)   SKIP_AF3=1; shift ;;
        --pareto)     PARETO_FLAG="--pareto"; shift ;;
        *)            echo "Unknown arg: $1"; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUTS_DIR="$SCRIPT_DIR/inputs"
CONFIG_TEMPLATE="$SCRIPT_DIR/pipeline_config_template.yaml"

eval "$(conda shell.bash hook)"
conda activate unified_mpnn

echo "============================================================================"
echo "G4FP Design Pipeline -- Full Batch Run"
echo "Started: $(date)"
echo "Parallel AF3 jobs: $PARALLEL"
echo "============================================================================"
echo ""

# ============================================================================
# STEP 1: LigandMPNN (all structures)
# ============================================================================
if [ $SKIP_LMPNN -eq 0 ]; then
    echo "=== STEP 1: LigandMPNN sequence design ==="
    "$SCRIPT_DIR/batch_run_ligandmpnn.sh"
    echo ""
else
    echo "=== STEP 1: SKIPPED (--skip-ligandmpnn) ==="
    echo ""
fi

# ============================================================================
# STEP 1b: Pre-filter analysis (all structures, parallel)
# ============================================================================
echo "=== STEP 1b: Pre-filter analysis ==="
for PDB_PATH in "$INPUTS_DIR"/*.pdb; do
    PDB_NAME=$(basename "$PDB_PATH" .pdb)
    OUTPUT_DIR="$SCRIPT_DIR/output_${PDB_NAME}"

    if [ ! -d "$OUTPUT_DIR/01_ligandmpnn/seqs" ]; then
        continue
    fi
    if [ -f "$OUTPUT_DIR/01b_pre_filter_analysis/pre_filter_summary.csv" ]; then
        echo "  SKIP $PDB_NAME (already done)"
        continue
    fi

    # Write temp config
    TEMP_CFG="$SCRIPT_DIR/_temp_1b_${PDB_NAME}.yaml"
    python3 -c "
import yaml
with open('$CONFIG_TEMPLATE') as f:
    c = yaml.safe_load(f)
c['input_pdb'] = '$PDB_PATH'
c['output_dir'] = '$OUTPUT_DIR/'
with open('$TEMP_CFG', 'w') as f:
    yaml.dump(c, f, default_flow_style=False, sort_keys=False)
"
    python3 "$SCRIPT_DIR/01b_analyze_pre_filter.py" "$TEMP_CFG" &
    rm -f "$TEMP_CFG"
done
wait
echo ""

# ============================================================================
# STEP 2: Filter and prepare AF3 inputs (all structures)
# ============================================================================
echo "=== STEP 2: Filter and prepare AF3 inputs ==="
for PDB_PATH in "$INPUTS_DIR"/*.pdb; do
    PDB_NAME=$(basename "$PDB_PATH" .pdb)
    OUTPUT_DIR="$SCRIPT_DIR/output_${PDB_NAME}"

    if [ ! -d "$OUTPUT_DIR/01_ligandmpnn/seqs" ]; then
        continue
    fi
    if [ -f "$OUTPUT_DIR/02_filtered_sequences.csv" ] && [ $SKIP_LMPNN -eq 1 ]; then
        echo "  SKIP $PDB_NAME (already filtered)"
        continue
    fi

    python3 "$SCRIPT_DIR/02_filter_and_prepare_af3.py" \
        --output-dir "$OUTPUT_DIR" --input-pdb "$PDB_PATH" $PARETO_FLAG
done
echo ""

# ============================================================================
# STEP 2.5: Validate AF3 JSONs
# ============================================================================
echo "=== STEP 2.5: Validate AF3 JSONs ==="
python3 "$SCRIPT_DIR/utils/validate_af3_jsons.py" --fix
echo ""

# ============================================================================
# STEP 3: AlphaFold3 predictions (parallel)
# ============================================================================
if [ $SKIP_AF3 -eq 0 ]; then
    echo "=== STEP 3: AlphaFold3 predictions ($PARALLEL parallel jobs per node, 2 GPUs) ==="
    "$SCRIPT_DIR/run_all_af3_parallel.sh" --jobs-per-node "$PARALLEL"
    echo ""
else
    echo "=== STEP 3: SKIPPED (--skip-af3) ==="
    echo ""
fi

# ============================================================================
# STEP 4a: LigandMPNN analysis (all structures)
# ============================================================================
echo "=== STEP 4a: LigandMPNN analysis ==="
for PDB_PATH in "$INPUTS_DIR"/*.pdb; do
    PDB_NAME=$(basename "$PDB_PATH" .pdb)
    OUTPUT_DIR="$SCRIPT_DIR/output_${PDB_NAME}"

    if [ ! -f "$OUTPUT_DIR/02_filtered_sequences.csv" ]; then
        continue
    fi

    TEMP_CFG="$SCRIPT_DIR/_temp_4a_${PDB_NAME}.yaml"
    python3 -c "
import yaml
with open('$CONFIG_TEMPLATE') as f:
    c = yaml.safe_load(f)
c['input_pdb'] = '$PDB_PATH'
c['output_dir'] = '$OUTPUT_DIR/'
with open('$TEMP_CFG', 'w') as f:
    yaml.dump(c, f, default_flow_style=False, sort_keys=False)
"
    python3 "$SCRIPT_DIR/04a_analyze_ligandmpnn.py" "$TEMP_CFG"
    rm -f "$TEMP_CFG"
done
echo ""

# ============================================================================
# STEP 5-6: Structure analysis + bound vs apo comparison
# ============================================================================
echo "=== STEPS 5-6: Structure analysis and bound vs apo comparison ==="
"$SCRIPT_DIR/run_structure_analysis.sh"
echo ""

# ============================================================================
# STEP 7: Aggregate and rank
# ============================================================================
echo "=== STEP 7: Aggregate and rank across all templates ==="
python3 "$SCRIPT_DIR/07_aggregate_and_rank.py"
echo ""

# ============================================================================
# STEP 8: Export for synthesis
# ============================================================================
echo "=== STEP 8: Export for synthesis ==="
python3 "$SCRIPT_DIR/08_export_for_synthesis.py"
echo ""

# ============================================================================
# DONE
# ============================================================================
echo "============================================================================"
echo "G4FP Pipeline Complete!"
echo "Finished: $(date)"
echo ""
echo "Final outputs:"
echo "  07_results/07_final_candidates.csv    -- Ranked candidates"
echo "  07_results/07_final_candidates.fa     -- FASTA sequences"
echo "  08_synthesis_order/                   -- Synthesis-ready files"
echo "============================================================================"

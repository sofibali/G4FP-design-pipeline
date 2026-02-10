#!/bin/bash
# Batch run 02_filter_and_prepare_af3.py for all structures missing filtered output
# Usage: ./batch_filter_all.sh [--force]
#
# Options:
#   --force    Re-filter even if output already exists

set +e

eval "$(conda shell.bash hook)"
conda activate unified_mpnn

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUTS_DIR="$SCRIPT_DIR/inputs"
CONFIG_TEMPLATE="$SCRIPT_DIR/pipeline_config_template.yaml"
CONFIG_FILE="$SCRIPT_DIR/pipeline_config.yaml"

FORCE=0
if [[ "$1" == "--force" ]]; then
    FORCE=1
fi

echo "============================================================================"
echo "Batch Filter and Prepare AF3 Inputs"
echo "============================================================================"
echo ""

COMPLETED=0
SKIPPED=0
FAILED=0
TOTAL=0

for PDB_PATH in "$INPUTS_DIR"/*.pdb; do
    PDB_NAME=$(basename "$PDB_PATH" .pdb)
    OUTPUT_DIR="$SCRIPT_DIR/output_${PDB_NAME}"
    ((TOTAL++))

    echo "--- [$TOTAL] $PDB_NAME ---"

    # Check LigandMPNN output exists
    if [ ! -d "$OUTPUT_DIR/01_ligandmpnn/seqs" ]; then
        echo "  SKIP: No LigandMPNN output"
        ((SKIPPED++))
        continue
    fi

    # Check if already filtered
    if [ $FORCE -eq 0 ] && [ -f "$OUTPUT_DIR/02_filtered_sequences.csv" ]; then
        N_BOUND=$(ls "$OUTPUT_DIR/02_alphafold3_inputs_bound/"*.json 2>/dev/null | wc -l)
        echo "  SKIP: Already filtered ($N_BOUND bound JSONs)"
        ((SKIPPED++))
        continue
    fi

    # Write a temporary config for this structure
    TEMP_CONFIG="$SCRIPT_DIR/_temp_config_${PDB_NAME}.yaml"
    python3 -c "
import yaml
with open('$CONFIG_TEMPLATE') as f:
    config = yaml.safe_load(f)
config['input_pdb'] = '$PDB_PATH'
config['output_dir'] = '$OUTPUT_DIR/'
with open('$TEMP_CONFIG', 'w') as f:
    yaml.dump(config, f, default_flow_style=False, sort_keys=False)
"
    if [ $? -ne 0 ]; then
        echo "  FAIL: Could not write temp config"
        ((FAILED++))
        rm -f "$TEMP_CONFIG"
        continue
    fi

    # Run filter
    python3 "$SCRIPT_DIR/02_filter_and_prepare_af3.py" "$TEMP_CONFIG"
    EXIT_CODE=$?
    rm -f "$TEMP_CONFIG"

    if [ $EXIT_CODE -eq 0 ] && [ -f "$OUTPUT_DIR/02_filtered_sequences.csv" ]; then
        N_BOUND=$(ls "$OUTPUT_DIR/02_alphafold3_inputs_bound/"*.json 2>/dev/null | wc -l)
        N_APO=$(ls "$OUTPUT_DIR/02_alphafold3_inputs_apo/"*.json 2>/dev/null | wc -l)
        echo "  OK: $N_BOUND bound + $N_APO apo JSONs"
        ((COMPLETED++))
    else
        echo "  FAIL: exit code $EXIT_CODE"
        ((FAILED++))
    fi
    echo ""
done

echo "============================================================================"
echo "Summary: $TOTAL total | $COMPLETED filtered | $SKIPPED skipped | $FAILED failed"
echo "============================================================================"

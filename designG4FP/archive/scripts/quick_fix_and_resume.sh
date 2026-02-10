
#!/bin/bash
# Quick fix and resume script for the LigandMPNN design_G4FP pipeline
# This script fixes the JSON dialect issue and prepares next steps

set -e

BASE_DIR="/home/sbali/LigandMPNN/designG4FP"
cd "$BASE_DIR"

echo "============================================================================"
echo "LigandMPNN Design G4FP Pipeline - Quick Fix and Resume"
echo "============================================================================"
echo ""

# Step 1: Fix JSON dialect issue for existing AF3 inputs
echo "Step 1: Fixing JSON dialect for existing AlphaFold3 inputs..."
echo "------------------------------------------------------------"

if [ -d "output_G4FP_des1_cro_mod0/02_alphafold3_inputs_bound" ]; then
    echo "Fixing bound predictions JSON files..."
    python fix_af3_json_dialect.py output_G4FP_des1_cro_mod0/02_alphafold3_inputs_bound
    echo ""
fi

if [ -d "output_G4FP_des1_cro_mod0/02_alphafold3_inputs_apo" ]; then
    echo "Fixing apo predictions JSON files..."
    python fix_af3_json_dialect.py output_G4FP_des1_cro_mod0/02_alphafold3_inputs_apo
    echo ""
fi

echo "✓ JSON files fixed!"
echo ""

# Step 2: Filter completed LigandMPNN outputs that haven't been filtered
echo "Step 2: Filtering completed LigandMPNN outputs..."
echo "------------------------------------------------------------"

NEEDS_FILTERING=(
    "output_G4FP_des1_cro_mod4"
    "output_G4FP_r20_l1_cro_mod0"
    "output_G4FP_r20_l2_cro_mod2"
)

for OUTPUT_DIR in "${NEEDS_FILTERING[@]}"; do
    if [ -d "$OUTPUT_DIR/01_ligandmpnn/seqs" ]; then
        if [ ! -f "$OUTPUT_DIR/02_filtered_sequences.csv" ]; then
            echo "Filtering $OUTPUT_DIR..."
            
            # Create a temporary config for this structure
            INPUT_NAME="${OUTPUT_DIR#output_}"
            CONFIG_FILE="$OUTPUT_DIR/pipeline_config_temp.yaml"
            
            # Copy template and update paths
            sed "s|input_pdb:.*|input_pdb: inputs/${INPUT_NAME}.pdb|g" pipeline_config.yaml | \
            sed "s|output_dir:.*|output_dir: ${OUTPUT_DIR}/|g" > "$CONFIG_FILE"
            
            # Run filtering
            python 02_filter_and_prepare_af3.py "$CONFIG_FILE"
            
            echo "✓ $OUTPUT_DIR filtered"
            echo ""
        else
            echo "✓ $OUTPUT_DIR already filtered"
            echo ""
        fi
    else
        echo "⚠ $OUTPUT_DIR: No LigandMPNN output found, skipping"
        echo ""
    fi
done

# Step 3: Generate summary of what's ready
echo ""
echo "============================================================================"
echo "Current Status"
echo "============================================================================"
echo ""

# Count ready for AF3
READY_FOR_AF3=0
for dir in output_G4FP_*/02_alphafold3_inputs_bound; do
    if [ -d "$dir" ]; then
        JSON_COUNT=$(ls -1 "$dir"/*.json 2>/dev/null | wc -l)
        if [ "$JSON_COUNT" -gt 0 ]; then
            READY_FOR_AF3=$((READY_FOR_AF3 + 1))
            OUTPUT_NAME=$(basename $(dirname "$dir"))
            echo "✓ $OUTPUT_NAME: $JSON_COUNT sequences ready for AlphaFold3"
        fi
    fi
done

echo ""
echo "Total structures ready for AlphaFold3: $READY_FOR_AF3"
echo ""

# Step 4: List what still needs to be done
echo "============================================================================"
echo "Next Steps"
echo "============================================================================"
echo ""

echo "1. Run AlphaFold3 predictions for ready structures:"
echo "   cd output_G4FP_des1_cro_mod0"
echo "   bash ../03_run_alphafold3_parallel.sh"
echo ""

echo "2. Complete LigandMPNN for missing structures:"
MISSING_STRUCTURES=(
    "G4FP_r20_l2_cro_mod3"
    "G4FP_r28_l1_cro_mod2"
    "G4FP_r28_l2_cro_mod2"
    "G4FP_r28_l2_cro_mod3"
    "G4FP_r28_l3_cro_mod0"
    "G4FP_r28_l3_cro_mod2"
)

for STRUCT in "${MISSING_STRUCTURES[@]}"; do
    if [ ! -d "output_${STRUCT}/01_ligandmpnn" ]; then
        echo "   - $STRUCT"
    fi
done

echo ""
echo "3. Check progress anytime:"
echo "   python check_pipeline_status.py"
echo ""

echo "============================================================================"
echo "Setup Complete!"
echo "============================================================================"

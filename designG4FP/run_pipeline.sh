#!/bin/bash
# Master Pipeline Script: Simplified Protein Design Pipeline
# Usage: ./run_pipeline.sh [config_file]

set -e  # Exit on error

CONFIG_FILE="${1:-pipeline_config.yaml}"

echo ""
echo "###############################################################################"
echo "#                                                                             #"
echo "#              SIMPLIFIED PROTEIN DESIGN PIPELINE                            #"
echo "#                                                                             #"
echo "#  LigandMPNN → Filter → AlphaFold3 → Analysis                               #"
echo "#                                                                             #"
echo "###############################################################################"
echo ""

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ ERROR: Configuration file not found: $CONFIG_FILE"
    echo ""
    echo "Please create a configuration file or use the default: pipeline_config.yaml"
    exit 1
fi

echo "Using configuration: $CONFIG_FILE"
echo ""

# Make scripts executable
chmod +x 01_run_ligandmpnn.sh
chmod +x 03_run_alphafold3.sh

# ============================================================================
# STEP 1: Run LigandMPNN
# ============================================================================
echo ""
read -p "Run STEP 1 (LigandMPNN sequence design)? [Y/n] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]] || [[ -z $REPLY ]]; then
    ./01_run_ligandmpnn.sh "$CONFIG_FILE"
else
    echo "Skipping Step 1"
fi

# ============================================================================
# STEP 1b: Analyze LigandMPNN output before filtering
# ============================================================================
echo ""
read -p "Run STEP 1b (Analyze LigandMPNN output before filtering)? [Y/n] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]] || [[ -z $REPLY ]]; then
    python3 01b_analyze_pre_filter.py "$CONFIG_FILE"
else
    echo "Skipping Step 1b"
fi

# ============================================================================
# STEP 2: Filter and prepare AF3 inputs
# ============================================================================
echo ""
read -p "Run STEP 2 (Filter sequences and prepare AF3 inputs)? [Y/n] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]] || [[ -z $REPLY ]]; then
    python3 02_filter_and_prepare_af3.py "$CONFIG_FILE"
else
    echo "Skipping Step 2"
fi

# ============================================================================
# STEP 3: Run AlphaFold3
# ============================================================================
echo ""
echo "STEP 3: AlphaFold3 Predictions"
echo "This step can take many hours depending on the number of designs."
echo ""
read -p "Run STEP 3 (AlphaFold3 predictions)? [y/N] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo ""
    echo "Choose prediction mode:"
    echo "  1) Bound only (protein+DNA)"
    echo "  2) Apo only (protein alone)"
    echo "  3) Both (bound and apo)"
    echo ""
    read -p "Select option [1/2/3]: " -n 1 -r MODE
    echo ""
    
    case $MODE in
        1)
            ./03_run_alphafold3.sh "$CONFIG_FILE" bound
            ;;
        2)
            ./03_run_alphafold3.sh "$CONFIG_FILE" apo
            ;;
        3)
            ./03_run_alphafold3.sh "$CONFIG_FILE" both
            ;;
        *)
            echo "Invalid option. Skipping Step 3"
            ;;
    esac
else
    echo "Skipping Step 3"
    echo ""
    echo "NOTE: You can run AlphaFold3 predictions manually later with:"
    echo "  ./03_run_alphafold3.sh $CONFIG_FILE [bound|apo|both]"
fi

# ============================================================================
# STEP 4a: Analyze LigandMPNN
# ============================================================================
echo ""
read -p "Run STEP 4a (Analyze LigandMPNN sequences)? [Y/n] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]] || [[ -z $REPLY ]]; then
    python3 04a_analyze_ligandmpnn.py "$CONFIG_FILE"
else
    echo "Skipping Step 4a"
fi

# ============================================================================
# STEP 4b: Analyze AF3 structures (independent)
# ============================================================================
echo ""
echo "STEP 4b: Analyze AlphaFold3 structures (independent analysis)"
echo ""
read -p "Run STEP 4b? [y/N] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    # Parse output directory from config
    OUTPUT_DIR=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG_FILE')); print(c['output_dir'])")
    TEMPLATE=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG_FILE')); print(c['input_pdb'])")
    CHROM_START=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG_FILE')); print(c['analysis']['chromophore_start'])")
    CHROM_END=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG_FILE')); print(c['analysis']['chromophore_end'])")
    
    # Run analysis for bound predictions
    if [ -d "$OUTPUT_DIR/03_alphafold3_predictions_bound" ]; then
        echo "Analyzing bound (protein+DNA) predictions..."
        python3 analyze_af3_structures.py \
            --pipeline-dir "$OUTPUT_DIR" \
            --af3-dir "$OUTPUT_DIR/03_alphafold3_predictions_bound" \
            --template "$TEMPLATE" \
            --chromophore-start "$CHROM_START" \
            --chromophore-end "$CHROM_END" \
            --output-dir "$OUTPUT_DIR/04b_structure_analysis_bound"
    fi
    
    # Run analysis for apo predictions
    if [ -d "$OUTPUT_DIR/03_alphafold3_predictions_apo" ]; then
        echo ""
        echo "Analyzing apo (protein-only) predictions..."
        python3 analyze_af3_structures.py \
            --pipeline-dir "$OUTPUT_DIR" \
            --af3-dir "$OUTPUT_DIR/03_alphafold3_predictions_apo" \
            --template "$TEMPLATE" \
            --chromophore-start "$CHROM_START" \
            --chromophore-end "$CHROM_END" \
            --output-dir "$OUTPUT_DIR/04b_structure_analysis_apo"
    fi
else
    echo "Skipping Step 4b"
fi

# ============================================================================
# STEP 4c: Compare bound vs apo
# ============================================================================
echo ""
echo "STEP 4c: Compare bound vs apo AlphaFold3 predictions"
echo ""
read -p "Run STEP 4c? [y/N] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    OUTPUT_DIR=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG_FILE')); print(c['output_dir'])")
    CHROM_START=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG_FILE')); print(c['analysis']['chromophore_start'])")
    CHROM_END=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG_FILE')); print(c['analysis']['chromophore_end'])")
    
    python3 compare_ligand_states.py \
        --bound-dir "$OUTPUT_DIR/03_alphafold3_predictions_bound" \
        --apo-dir "$OUTPUT_DIR/03_alphafold3_predictions_apo" \
        --chromophore-start "$CHROM_START" \
        --chromophore-end "$CHROM_END" \
        --output-dir "$OUTPUT_DIR/04c_ligand_state_comparison"
else
    echo "Skipping Step 4c"
fi

# ============================================================================
# SUMMARY
# ============================================================================
echo ""
echo "###############################################################################"
echo "#                                                                             #"
echo "#                          PIPELINE COMPLETE                                  #"
echo "#                                                                             #"
echo "###############################################################################"
echo ""

OUTPUT_DIR=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG_FILE')); print(c['output_dir'])" 2>/dev/null || echo "pipeline_outputs")

echo "Results directory: $OUTPUT_DIR/"
echo ""
echo "Pipeline outputs:"
echo "  01_ligandmpnn/                    - LigandMPNN generated sequences"
echo "  01b_pre_filter_analysis/          - Pre-filter sequence analysis"
echo "  02_filtered_sequences.csv         - Top filtered sequences"
echo "  02_alphafold3_inputs_bound/       - AF3 inputs (protein+DNA)"
echo "  02_alphafold3_inputs_apo/         - AF3 inputs (protein only)"
echo "  03_alphafold3_predictions_bound/  - AF3 structures (protein+DNA)"
echo "  03_alphafold3_predictions_apo/    - AF3 structures (protein only)"
echo "  04a_ligandmpnn_analysis/          - Sequence diversity analysis"
echo "  04b_structure_analysis_bound/     - Structure analysis (bound)"
echo "  04b_structure_analysis_apo/       - Structure analysis (apo)"
echo "  04c_ligand_state_comparison/      - Bound vs apo comparison"
echo ""
echo "For detailed usage instructions, see PIPELINE_README.md"
echo ""

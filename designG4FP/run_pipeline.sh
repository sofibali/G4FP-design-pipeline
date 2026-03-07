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
    # Auto-discovers all output_G4FP_* dirs; chromophore range is auto-detected
    # per template from the PDB unless overridden via --chromophore-range.
    # To override: --chromophore-range "197,199"
    python3 05_analyze_af3_structures.py \
        --state both
else
    echo "Skipping Step 4b"
fi

# ============================================================================
# STEP 4c: Compare bound vs apo (all templates combined)
# ============================================================================
echo ""
echo "STEP 4c: Compare bound vs apo AlphaFold3 predictions (all templates)"
echo ""
read -p "Run STEP 4c? [y/N] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    # Auto-discovers all output_G4FP_* dirs; chromophore range is auto-detected
    # per template from the PDB unless overridden via --chromophore-range.
    python3 06_compare_ligand_states.py
else
    echo "Skipping Step 4c"
fi

# ============================================================================
# STEP 5: Aggregate and rank across ALL output directories
# ============================================================================
echo ""
echo "STEP 5: Aggregate and rank all designs across all templates"
echo "This combines results from all output_G4FP_* directories into a single ranking."
echo ""
read -p "Run STEP 5 (Aggregate all 1000 sequences and rank)? [y/N] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    python3 07_aggregate_and_rank.py
else
    echo "Skipping Step 5"
    echo ""
    echo "NOTE: You can run aggregation manually later with:"
    echo "  python3 07_aggregate_and_rank.py"
    echo "  python3 07_aggregate_and_rank.py --n-select 500"
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
echo "Per-template outputs:"
echo "  01_ligandmpnn/                    - LigandMPNN generated sequences"
echo "  01b_pre_filter_analysis/          - Pre-filter sequence analysis"
echo "  02_filtered_sequences.csv         - Top filtered sequences"
echo "  02_alphafold3_inputs_bound/       - AF3 inputs (protein+DNA)"
echo "  02_alphafold3_inputs_apo/         - AF3 inputs (protein only)"
echo "  03_alphafold3_predictions_bound/  - AF3 structures (protein+DNA)"
echo "  03_alphafold3_predictions_apo/    - AF3 structures (protein only)"
echo "  04a_ligandmpnn_analysis/          - Sequence diversity analysis"
echo "  05_structure_analysis/            - Structure analysis (bound+apo)"
echo "  06_ligand_state_comparison/       - Bound vs apo comparison"
echo ""
echo "Aggregated outputs (all templates combined):"
echo "  07_results/                       - Final ranked candidates across all templates"
echo "    07_final_candidates.csv         - Ranked candidates with all metrics"
echo "    07_final_candidates.fa          - FASTA of selected sequences"
echo "    07_pareto_frontier.csv          - Pareto-optimal designs"
echo "    07_selection_summary.png        - Overview plots"
echo "    07_bound_vs_apo_scatter.png     - Bound vs apo visualization"
echo ""
echo "For detailed usage instructions, see PIPELINE_README.md"
echo ""

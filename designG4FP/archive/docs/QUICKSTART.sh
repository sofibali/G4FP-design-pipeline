#!/bin/bash
# QUICK START: Simplified Protein Design Pipeline
# Copy and paste these commands to run the pipeline

# ============================================================================
# SETUP (One-time)
# ============================================================================

# 1. Edit configuration file
nano pipeline_config.yaml
# - Update input_pdb path
# - Update design_residues
# - Update AlphaFold3 db_dir and model_dir

# 2. Install Python dependencies (if needed)
pip install biopython pandas numpy matplotlib seaborn pyyaml

# ============================================================================
# RUN COMPLETE PIPELINE (Interactive)
# ============================================================================

./run_pipeline.sh

# ============================================================================
# RUN INDIVIDUAL STEPS
# ============================================================================

# STEP 1: Generate sequences with LigandMPNN (~10-30 min for 10K sequences)
./01_run_ligandmpnn.sh

# STEP 1b: Analyze all sequences before filtering (~1 min)
python 01b_analyze_pre_filter.py

# STEP 2: Filter top 10% and prepare AF3 inputs (~1 min)
python 02_filter_and_prepare_af3.py

# STEP 3: Run AlphaFold3 predictions (SLOW - hours to days)
./03_run_alphafold3.sh bound    # Protein+DNA structures
./03_run_alphafold3.sh apo      # Protein-only structures  
./03_run_alphafold3.sh both     # Both (sequential)

# STEP 4a: Analyze LigandMPNN sequences (~1 min)
python 04a_analyze_ligandmpnn.py

# STEP 4b: Analyze AF3 structures - Bound state (~5-10 min)
python analyze_af3_structures.py \
    --pipeline-dir ./pipeline_outputs \
    --af3-dir ./pipeline_outputs/03_alphafold3_predictions_bound \
    --template ./inputs/G4FP_designs/G4FP_r28_l3_cro_mod2.pdb \
    --chromophore-start 197 \
    --chromophore-end 199 \
    --output-dir ./pipeline_outputs/04b_structure_analysis_bound

# STEP 4b: Analyze AF3 structures - Apo state (~5-10 min)
python analyze_af3_structures.py \
    --pipeline-dir ./pipeline_outputs \
    --af3-dir ./pipeline_outputs/03_alphafold3_predictions_apo \
    --template ./inputs/G4FP_designs/G4FP_r28_l3_cro_mod2.pdb \
    --chromophore-start 197 \
    --chromophore-end 199 \
    --output-dir ./pipeline_outputs/04b_structure_analysis_apo

# STEP 4c: Compare bound vs apo (~5-10 min)
python compare_ligand_states.py \
    --bound-dir ./pipeline_outputs/03_alphafold3_predictions_bound \
    --apo-dir ./pipeline_outputs/03_alphafold3_predictions_apo \
    --chromophore-start 197 \
    --chromophore-end 199 \
    --output-dir ./pipeline_outputs/04c_ligand_state_comparison

# ============================================================================
# QUICK TESTS (Reduce compute time for testing)
# ============================================================================

# Edit config for fast test:
# - num_sequences: 100  (instead of 10000)
# - top_percent: 10     (gives ~10 designs)
# - modelseeds: [0]     (1 prediction instead of 5)

# ============================================================================
# OUTPUT LOCATIONS
# ============================================================================

# All outputs in: ./pipeline_outputs/
#
# Key files:
#   02_filtered_sequences.csv           - Filtered sequences table
#   04a_ligandmpnn_analysis/            - Sequence diversity plots
#   04b_structure_analysis_bound/       - Structure metrics (bound)
#   04b_structure_analysis_apo/         - Structure metrics (apo)  
#   04c_ligand_state_comparison/        - Bound vs apo comparison
#   */plot_data_csvs/                   - CSV exports for Prism/GraphPad

# ============================================================================
# HELP
# ============================================================================

# See detailed documentation:
less PIPELINE_README.md

# Check individual script help:
python 02_filter_and_prepare_af3.py --help
python analyze_af3_structures.py --help
python compare_ligand_states.py --help

#!/bin/bash
# Run Comprehensive Structure Analysis for designG4FP Outputs
# 
# This script runs both:
#   1. Individual structure analysis (05_analyze_af3_structures.py)
#   2. Ligand state comparison (06_compare_ligand_states.py)
#
# Similar to previous pipeline's 07_structure_analysis and 08_ligand_state_comparison
#
# Usage:
#   bash run_structure_analysis.sh [OUTPUT_DIR] [TEMPLATE_PDB]
#
# Arguments:
#   OUTPUT_DIR: Specific output directory (default: all output_* directories)
#   TEMPLATE_PDB: Path to template structure (required)
#
# Example:
#   bash run_structure_analysis.sh output_G4FP_des1_cro_mod0 inputs/G4FP_des1.pdb

set -euo pipefail

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Configuration
CHROMOPHORE_RANGE="175,225"  # Adjust if needed

# Parse arguments
if [ $# -lt 1 ]; then
    echo -e "${RED}Error: Template structure required${NC}"
    echo "Usage: $0 [OUTPUT_DIR] TEMPLATE_PDB"
    echo "Example: $0 output_G4FP_des1_cro_mod0 inputs/G4FP_des1.pdb"
    exit 1
fi

if [ $# -eq 1 ]; then
    # Only template provided, analyze all outputs
    TEMPLATE_PDB="$1"
    OUTPUT_DIRS=(output_*)
    
    if [ ${#OUTPUT_DIRS[@]} -eq 0 ]; then
        echo -e "${RED}Error: No output directories found${NC}"
        exit 1
    fi
    
    echo -e "${BLUE}Found ${#OUTPUT_DIRS[@]} output directories to analyze${NC}"
else
    # Specific output directory provided
    OUTPUT_DIR="$1"
    TEMPLATE_PDB="$2"
    OUTPUT_DIRS=("$OUTPUT_DIR")
fi

# Check template exists
if [ ! -f "$TEMPLATE_PDB" ]; then
    echo -e "${RED}Error: Template structure not found: $TEMPLATE_PDB${NC}"
    exit 1
fi

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}G4FP Structure Analysis Pipeline${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Template: $TEMPLATE_PDB"
echo "Chromophore range: $CHROMOPHORE_RANGE"
echo "Output directories: ${#OUTPUT_DIRS[@]}"
echo ""

# Process each output directory
for OUTPUT_DIR in "${OUTPUT_DIRS[@]}"; do
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}Processing: $OUTPUT_DIR${NC}"
    echo -e "${BLUE}========================================${NC}"
    
    # Check if AF3 predictions exist
    BOUND_DIR="$OUTPUT_DIR/03_alphafold3_predictions_bound"
    APO_DIR="$OUTPUT_DIR/03_alphafold3_predictions_apo"
    
    BOUND_EXISTS=false
    APO_EXISTS=false
    
    if [ -d "$BOUND_DIR" ] && [ "$(ls -A $BOUND_DIR 2>/dev/null)" ]; then
        BOUND_EXISTS=true
    fi
    
    if [ -d "$APO_DIR" ] && [ "$(ls -A $APO_DIR 2>/dev/null)" ]; then
        APO_EXISTS=true
    fi
    
    if [ "$BOUND_EXISTS" = false ] && [ "$APO_EXISTS" = false ]; then
        echo -e "${YELLOW}⚠ No AF3 predictions found, skipping${NC}"
        echo ""
        continue
    fi
    
    # Determine which states to analyze
    if [ "$BOUND_EXISTS" = true ] && [ "$APO_EXISTS" = true ]; then
        STATE="both"
        echo -e "${GREEN}✓ Found both bound and apo predictions${NC}"
    elif [ "$BOUND_EXISTS" = true ]; then
        STATE="bound"
        echo -e "${YELLOW}⚠ Only bound predictions found${NC}"
    else
        STATE="apo"
        echo -e "${YELLOW}⚠ Only apo predictions found${NC}"
    fi
    
    # 1. Run individual structure analysis (05)
    echo ""
    echo -e "${BLUE}[1/2] Running structure analysis...${NC}"
    python 05_analyze_af3_structures.py \
        --output-dir "$OUTPUT_DIR" \
        --template "$TEMPLATE_PDB" \
        --state "$STATE" \
        --chromophore-range "$CHROMOPHORE_RANGE"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Structure analysis complete${NC}"
        echo "  Results: $OUTPUT_DIR/05_structure_analysis_results.csv"
        echo "  Plots: $OUTPUT_DIR/05_structure_analysis/"
    else
        echo -e "${RED}✗ Structure analysis failed${NC}"
    fi
    
    # 2. Run ligand state comparison (06) - only if both states exist
    if [ "$BOUND_EXISTS" = true ] && [ "$APO_EXISTS" = true ]; then
        echo ""
        echo -e "${BLUE}[2/2] Running ligand state comparison...${NC}"
        python 06_compare_ligand_states.py \
            --output-dir "$OUTPUT_DIR" \
            --template "$TEMPLATE_PDB" \
            --chromophore-range "$CHROMOPHORE_RANGE"
        
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}✓ Ligand state comparison complete${NC}"
            echo "  Results: $OUTPUT_DIR/06_ligand_state_comparison_results.csv"
            echo "  Plots: $OUTPUT_DIR/06_ligand_state_comparison/"
        else
            echo -e "${RED}✗ Ligand state comparison failed${NC}"
        fi
    else
        echo ""
        echo -e "${YELLOW}⚠ Skipping ligand state comparison (need both bound and apo)${NC}"
    fi
    
    echo ""
    echo -e "${GREEN}✓ Analysis complete for $OUTPUT_DIR${NC}"
    echo ""
done

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}All analyses complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Output structure for each directory:"
echo "  - 05_structure_analysis_results.csv"
echo "  - 05_structure_analysis/*.png"
echo "  - 05_structure_analysis/plot_data_csvs/*.csv"
echo "  - 06_ligand_state_comparison_results.csv"
echo "  - 06_ligand_state_comparison/*.png"
echo "  - 06_ligand_state_comparison/plot_data_csvs/*.csv"
echo ""

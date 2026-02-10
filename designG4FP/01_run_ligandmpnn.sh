#!/bin/bash
# Step 1: Run LigandMPNN for sequence design
# Usage: ./01_run_ligandmpnn.sh [config_file]

set -e  # Exit on error

# Activate LigandMPNN environment
eval "$(conda shell.bash hook)"
conda activate unified_mpnn

# Load configuration
CONFIG_FILE="${1:-pipeline_config.yaml}"

echo "============================================================================"
echo "STEP 1: Running LigandMPNN Sequence Design"
echo "============================================================================"
echo ""

# Parse YAML config and expand residue ranges (requires python and PyYAML)
parse_config() {
    python3 -c "
import yaml
import re

with open('$CONFIG_FILE') as f:
    config = yaml.safe_load(f)

# Function to expand residue ranges like 'A15-60' to 'A15 A16 A17 ... A60'
def expand_residues(residue_spec):
    parts = residue_spec.strip().split()
    expanded = []
    
    for part in parts:
        # Match patterns like A15-60 or A15
        match = re.match(r'([A-Z]+)(\d+)-(\d+)', part)
        if match:
            chain, start, end = match.groups()
            expanded.extend([f'{chain}{i}' for i in range(int(start), int(end) + 1)])
        else:
            # Single residue like A15
            expanded.append(part)
    
    return ' '.join(expanded)

# Extract and expand residues
design_residues = expand_residues(config['ligandmpnn']['design_residues'])

# Extract values (with proper quoting for bash)
print(f\"INPUT_PDB='{config['input_pdb']}'\")
print(f\"OUTPUT_DIR='{config['output_dir']}'\")
print(f\"DESIGN_RESIDUES='{design_residues}'\")
print(f\"NUM_SEQUENCES={config['ligandmpnn']['num_sequences']}\")
print(f\"BATCH_SIZE={config['ligandmpnn']['batch_size']}\")
print(f\"SEED={config['ligandmpnn']['seed']}\")
print(f\"MODEL_CKPT='{config['ligandmpnn']['model_checkpoint']}'\")
print(f\"SC_CKPT='{config['ligandmpnn']['sc_checkpoint']}'\")
print(f\"CUTOFF={config['ligandmpnn']['cutoff_for_score']}\")
"
}

# Load config into environment
eval $(parse_config)

# Verify input PDB exists
if [ ! -f "$INPUT_PDB" ]; then
    echo "❌ ERROR: Input PDB not found: $INPUT_PDB"
    echo "   Please check the 'input_pdb' path in $CONFIG_FILE"
    exit 1
fi

# Create output directory
LIGANDMPNN_OUTPUT="$OUTPUT_DIR/01_ligandmpnn"
mkdir -p "$LIGANDMPNN_OUTPUT"

# Calculate number of batches
NUM_BATCHES=$(python3 -c "import math; print(math.ceil($NUM_SEQUENCES / $BATCH_SIZE))")

echo "Configuration:"
echo "  Input PDB:        $INPUT_PDB"
echo "  Output directory: $LIGANDMPNN_OUTPUT"
echo "  Design residues:  $DESIGN_RESIDUES"
echo "  Num sequences:    $NUM_SEQUENCES"
echo "  Batch size:       $BATCH_SIZE"
echo "  Num batches:      $NUM_BATCHES"
echo ""

# Check if output already exists
if [ -d "$LIGANDMPNN_OUTPUT/seqs" ] && [ -n "$(ls -A $LIGANDMPNN_OUTPUT/seqs/*.fa 2>/dev/null)" ]; then
    echo "✓ LigandMPNN output already exists - SKIPPING"
    echo "  Found existing FASTA files in $LIGANDMPNN_OUTPUT/seqs/"
    echo "  To regenerate, delete the directory and re-run this script."
    echo ""
    exit 0
fi

echo "Starting LigandMPNN (this may take 10-30 minutes)..."
echo ""

# Run LigandMPNN
python /home/sbali/LigandMPNN/run.py \
    --model_type ligand_mpnn \
    --seed $SEED \
    --pdb_path "$INPUT_PDB" \
    --out_folder "$LIGANDMPNN_OUTPUT" \
    --checkpoint_ligand_mpnn "$MODEL_CKPT" \
    --checkpoint_path_sc "$SC_CKPT" \
    --redesigned_residues "$DESIGN_RESIDUES" \
    --batch_size $BATCH_SIZE \
    --number_of_batches $NUM_BATCHES \
    --pack_side_chains 1 \
    --number_of_packs_per_design 1 \
    --ligand_mpnn_use_atom_context 1 \
    --ligand_mpnn_cutoff_for_score $CUTOFF \
    --pack_with_ligand_context 1 \
    --save_stats 1 \
    --verbose 1

echo ""
echo "============================================================================"
echo "✓ LigandMPNN completed successfully!"
echo "============================================================================"
echo "Output saved to: $LIGANDMPNN_OUTPUT"
echo ""

# Count generated sequences
if [ -d "$LIGANDMPNN_OUTPUT/seqs" ]; then
    FASTA_COUNT=$(ls -1 "$LIGANDMPNN_OUTPUT/seqs"/*.fa 2>/dev/null | wc -l)
    echo "Generated FASTA files: $FASTA_COUNT"
    
    # Count total sequences in FASTA files
    TOTAL_SEQS=$(grep -c "^>" "$LIGANDMPNN_OUTPUT/seqs"/*.fa 2>/dev/null || echo "0")
    echo "Total sequences: $TOTAL_SEQS"
fi

echo ""
echo "Next step: Run 02_filter_and_prepare_af3.py to filter and prepare AlphaFold3 inputs"
echo ""

#!/bin/bash
# Test script to verify LigandMPNN generates diverse sequences
# This runs just 10 sequences to quickly verify the fix

set -e

# Activate LigandMPNN environment
eval "$(conda shell.bash hook)"
conda activate unified_mpnn

CONFIG_FILE="pipeline_config.yaml"
TEST_OUTPUT="/home/sbali/LigandMPNN/designG4FP/test_output"

echo "=========================================="
echo "Testing LigandMPNN Fix - Small Test Run"
echo "=========================================="
echo ""

# Parse config with residue expansion
parse_config() {
    python3 -c "
import yaml
import re

with open('$CONFIG_FILE') as f:
    config = yaml.safe_load(f)

def expand_residues(residue_spec):
    parts = residue_spec.strip().split()
    expanded = []
    for part in parts:
        match = re.match(r'([A-Z]+)(\d+)-(\d+)', part)
        if match:
            chain, start, end = match.groups()
            expanded.extend([f'{chain}{i}' for i in range(int(start), int(end) + 1)])
        else:
            expanded.append(part)
    return ' '.join(expanded)

design_residues = expand_residues(config['ligandmpnn']['design_residues'])

print(f\"INPUT_PDB='{config['input_pdb']}'\")
print(f\"DESIGN_RESIDUES='{design_residues}'\")
print(f\"MODEL_CKPT='{config['ligandmpnn']['model_checkpoint']}'\")
print(f\"SC_CKPT='{config['ligandmpnn']['sc_checkpoint']}'\")
print(f\"CUTOFF={config['ligandmpnn']['cutoff_for_score']}\")
"
}

eval $(parse_config)

echo "Input PDB: $INPUT_PDB"
echo "Number of design positions: $(echo $DESIGN_RESIDUES | wc -w)"
echo "Test output: $TEST_OUTPUT"
echo ""

mkdir -p "$TEST_OUTPUT"

echo "Running LigandMPNN with 10 test sequences..."
python /home/sbali/LigandMPNN/run.py \
    --model_type ligand_mpnn \
    --seed 123 \
    --pdb_path "$INPUT_PDB" \
    --out_folder "$TEST_OUTPUT" \
    --checkpoint_ligand_mpnn "$MODEL_CKPT" \
    --checkpoint_path_sc "$SC_CKPT" \
    --redesigned_residues "$DESIGN_RESIDUES" \
    --batch_size 5 \
    --number_of_batches 2 \
    --pack_side_chains 1 \
    --number_of_packs_per_design 1 \
    --ligand_mpnn_use_atom_context 1 \
    --ligand_mpnn_cutoff_for_score $CUTOFF \
    --pack_with_ligand_context 1 \
    --save_stats 1 \
    --verbose 1

echo ""
echo "=========================================="
echo "Analyzing results..."
echo "=========================================="

# Check if sequences are diverse
python3 << 'PYEOF'
from pathlib import Path

fasta_files = list(Path("test_output/seqs").glob("*.fa"))
if not fasta_files:
    print("❌ No FASTA files found!")
    exit(1)

sequences = []
with open(fasta_files[0]) as f:
    for line in f:
        if line.startswith('>'):
            header = line.strip()
            # Extract num_res from header
            if 'num_res=' in header:
                num_res = int(header.split('num_res=')[1].split(',')[0])
                print(f"✓ num_res = {num_res} (should be > 0)")
        else:
            sequences.append(line.strip())

unique_seqs = len(set(sequences))
total_seqs = len(sequences)

print(f"\nSequence diversity:")
print(f"  Total sequences: {total_seqs}")
print(f"  Unique sequences: {unique_seqs}")
print(f"  Diversity: {unique_seqs/total_seqs*100:.1f}%")

if unique_seqs > 1:
    print("\n✓✓✓ SUCCESS! LigandMPNN is generating diverse sequences!")
else:
    print("\n❌ PROBLEM! All sequences are identical!")
    
PYEOF

echo ""
echo "Test output saved to: $TEST_OUTPUT"
echo "If successful, you can now re-run the full pipeline."

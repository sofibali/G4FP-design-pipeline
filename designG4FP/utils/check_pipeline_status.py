#!/usr/bin/env python3
"""
Check the status of the LigandMPNN design_G4FP pipeline.
Reports missing LigandMPNN outputs and AlphaFold3 predictions.
"""

import os
import sys
from pathlib import Path
import glob

# Base directory
BASE_DIR = Path("/home/sbali/LigandMPNN/designG4FP")
INPUTS_DIR = BASE_DIR / "inputs"

# Find all input PDB files
input_pdbs = sorted(INPUTS_DIR.glob("*.pdb"))
input_names = [p.stem for p in input_pdbs]

print("=" * 80)
print("LigandMPNN Design G4FP Pipeline Status Check")
print("=" * 80)
print()

# Expected output directories
expected_outputs = [f"output_{name}" for name in input_names]
existing_outputs = sorted([d.name for d in BASE_DIR.glob("output_*") if d.is_dir()])

print(f"Total input PDBs: {len(input_pdbs)}")
print(f"Expected output directories: {len(expected_outputs)}")
print(f"Existing output directories: {len(existing_outputs)}")
print()

# Find missing LigandMPNN runs
missing_ligandmpnn = []
for expected in expected_outputs:
    output_dir = BASE_DIR / expected
    ligandmpnn_dir = output_dir / "01_ligandmpnn"
    
    if not output_dir.exists():
        missing_ligandmpnn.append(expected)
    elif not ligandmpnn_dir.exists():
        missing_ligandmpnn.append(expected)
    elif not (ligandmpnn_dir / "seqs").exists():
        missing_ligandmpnn.append(expected)
    elif not list((ligandmpnn_dir / "seqs").glob("*.fa")):
        missing_ligandmpnn.append(expected)

print("=" * 80)
print("1. MISSING LIGANDMPNN OUTPUTS")
print("=" * 80)
print()

if missing_ligandmpnn:
    print(f"Missing LigandMPNN runs: {len(missing_ligandmpnn)}")
    print()
    for name in missing_ligandmpnn:
        input_name = name.replace("output_", "")
        input_file = INPUTS_DIR / f"{input_name}.pdb"
        print(f"  ❌ {name}")
        print(f"     Input: {input_file}")
        print()
else:
    print("✓ All LigandMPNN runs completed!")
    print()

# Check existing LigandMPNN outputs in detail
print("=" * 80)
print("2. LIGANDMPNN OUTPUT DETAILS")
print("=" * 80)
print()

ligandmpnn_summary = []
for output_name in existing_outputs:
    output_dir = BASE_DIR / output_name
    ligandmpnn_dir = output_dir / "01_ligandmpnn"
    
    status = {
        'name': output_name,
        'has_seqs': False,
        'fasta_count': 0,
        'total_sequences': 0,
        'has_scores': False,
        'has_backbones': False,
        'has_sidechains': False
    }
    
    if ligandmpnn_dir.exists():
        seqs_dir = ligandmpnn_dir / "seqs"
        if seqs_dir.exists():
            status['has_seqs'] = True
            fasta_files = list(seqs_dir.glob("*.fa"))
            status['fasta_count'] = len(fasta_files)
            
            # Count total sequences
            total_seqs = 0
            for fa in fasta_files:
                with open(fa) as f:
                    total_seqs += sum(1 for line in f if line.startswith('>'))
            status['total_sequences'] = total_seqs
        
        if (ligandmpnn_dir / "scores").exists():
            status['has_scores'] = True
        if (ligandmpnn_dir / "backbones").exists():
            status['has_backbones'] = True
        if (ligandmpnn_dir / "sidechains").exists():
            status['has_sidechains'] = True
    
    ligandmpnn_summary.append(status)

for s in ligandmpnn_summary:
    print(f"Directory: {s['name']}")
    if s['has_seqs']:
        print(f"  ✓ FASTA files: {s['fasta_count']}")
        print(f"  ✓ Total sequences: {s['total_sequences']}")
    else:
        print(f"  ❌ No sequences generated")
    
    if s['has_scores']:
        print(f"  ✓ Scores available")
    if s['has_backbones']:
        print(f"  ✓ Backbones available")
    if s['has_sidechains']:
        print(f"  ✓ Sidechains available")
    print()

# Check filtered sequences and AF3 inputs
print("=" * 80)
print("3. FILTERING AND AF3 INPUT PREPARATION")
print("=" * 80)
print()

filtering_summary = []
for output_name in existing_outputs:
    output_dir = BASE_DIR / output_name
    
    status = {
        'name': output_name,
        'has_filtered_csv': False,
        'has_filtered_fasta': False,
        'filtered_count': 0,
        'has_af3_inputs_bound': False,
        'has_af3_inputs_apo': False,
        'af3_inputs_bound_count': 0,
        'af3_inputs_apo_count': 0
    }
    
    filtered_csv = output_dir / "02_filtered_sequences.csv"
    filtered_fasta = output_dir / "02_filtered_sequences.fa"
    
    if filtered_csv.exists():
        status['has_filtered_csv'] = True
        with open(filtered_csv) as f:
            status['filtered_count'] = sum(1 for line in f) - 1  # Subtract header
    
    if filtered_fasta.exists():
        status['has_filtered_fasta'] = True
    
    af3_inputs_bound = output_dir / "02_alphafold3_inputs_bound"
    af3_inputs_apo = output_dir / "02_alphafold3_inputs_apo"
    
    if af3_inputs_bound.exists():
        status['has_af3_inputs_bound'] = True
        status['af3_inputs_bound_count'] = len(list(af3_inputs_bound.glob("*.json")))
    
    if af3_inputs_apo.exists():
        status['has_af3_inputs_apo'] = True
        status['af3_inputs_apo_count'] = len(list(af3_inputs_apo.glob("*.json")))
    
    filtering_summary.append(status)

for s in filtering_summary:
    print(f"Directory: {s['name']}")
    if s['has_filtered_csv']:
        print(f"  ✓ Filtered CSV: {s['filtered_count']} sequences")
    else:
        print(f"  ❌ No filtered CSV")
    
    if s['has_filtered_fasta']:
        print(f"  ✓ Filtered FASTA available")
    else:
        print(f"  ❌ No filtered FASTA")
    
    if s['has_af3_inputs_bound']:
        print(f"  ✓ AF3 inputs (bound): {s['af3_inputs_bound_count']} files")
    else:
        print(f"  ⚠ No AF3 inputs (bound)")
    
    if s['has_af3_inputs_apo']:
        print(f"  ✓ AF3 inputs (apo): {s['af3_inputs_apo_count']} files")
    else:
        print(f"  ⚠ No AF3 inputs (apo)")
    print()

# Check AlphaFold3 predictions
print("=" * 80)
print("4. ALPHAFOLD3 PREDICTIONS STATUS")
print("=" * 80)
print()

af3_summary = []
missing_af3_predictions = []

for output_name in existing_outputs:
    output_dir = BASE_DIR / output_name
    
    status = {
        'name': output_name,
        'expected_bound': 0,
        'completed_bound': 0,
        'expected_apo': 0,
        'completed_apo': 0,
        'missing_bound': [],
        'missing_apo': []
    }
    
    # Count expected predictions
    af3_inputs_bound = output_dir / "02_alphafold3_inputs_bound"
    af3_inputs_apo = output_dir / "02_alphafold3_inputs_apo"
    
    if af3_inputs_bound.exists():
        input_files = list(af3_inputs_bound.glob("*.json"))
        status['expected_bound'] = len(input_files)
        
        # Check which predictions exist
        af3_pred_bound = output_dir / "03_alphafold3_predictions_bound"
        if af3_pred_bound.exists():
            for input_file in input_files:
                seq_name = input_file.stem
                pred_dir = af3_pred_bound / seq_name
                if pred_dir.exists() and list(pred_dir.glob("*fold_model_*.cif")):
                    status['completed_bound'] += 1
                else:
                    status['missing_bound'].append(seq_name)
        else:
            status['missing_bound'] = [f.stem for f in input_files]
    
    if af3_inputs_apo.exists():
        input_files = list(af3_inputs_apo.glob("*.json"))
        status['expected_apo'] = len(input_files)
        
        # Check which predictions exist
        af3_pred_apo = output_dir / "03_alphafold3_predictions_apo"
        if af3_pred_apo.exists():
            for input_file in input_files:
                seq_name = input_file.stem
                pred_dir = af3_pred_apo / seq_name
                if pred_dir.exists() and list(pred_dir.glob("*fold_model_*.cif")):
                    status['completed_apo'] += 1
                else:
                    status['missing_apo'].append(seq_name)
        else:
            status['missing_apo'] = [f.stem for f in input_files]
    
    af3_summary.append(status)
    
    if status['missing_bound'] or status['missing_apo']:
        missing_af3_predictions.append(status)

for s in af3_summary:
    print(f"Directory: {s['name']}")
    
    if s['expected_bound'] > 0:
        completion_pct = (s['completed_bound'] / s['expected_bound']) * 100
        print(f"  Bound predictions: {s['completed_bound']}/{s['expected_bound']} ({completion_pct:.1f}%)")
        if s['missing_bound']:
            print(f"    ❌ Missing: {len(s['missing_bound'])} predictions")
    else:
        print(f"  ⚠ No bound predictions expected")
    
    if s['expected_apo'] > 0:
        completion_pct = (s['completed_apo'] / s['expected_apo']) * 100
        print(f"  Apo predictions: {s['completed_apo']}/{s['expected_apo']} ({completion_pct:.1f}%)")
        if s['missing_apo']:
            print(f"    ❌ Missing: {len(s['missing_apo'])} predictions")
    else:
        print(f"  ⚠ No apo predictions expected")
    print()

# Save detailed missing predictions
if missing_af3_predictions:
    print("=" * 80)
    print("5. DETAILED MISSING AF3 PREDICTIONS")
    print("=" * 80)
    print()
    
    # Save to file
    missing_file = BASE_DIR / "missing_af3_predictions.txt"
    with open(missing_file, 'w') as f:
        f.write("Missing AlphaFold3 Predictions\n")
        f.write("=" * 80 + "\n\n")
        
        for s in missing_af3_predictions:
            print(f"{s['name']}:")
            f.write(f"{s['name']}:\n")
            
            if s['missing_bound']:
                print(f"  Missing bound predictions ({len(s['missing_bound'])}):")
                f.write(f"  Missing bound predictions ({len(s['missing_bound'])}):\n")
                for seq in s['missing_bound'][:10]:  # Show first 10
                    print(f"    - {seq}")
                    f.write(f"    - {seq}\n")
                if len(s['missing_bound']) > 10:
                    print(f"    ... and {len(s['missing_bound']) - 10} more")
                    f.write(f"    ... and {len(s['missing_bound']) - 10} more\n")
            
            if s['missing_apo']:
                print(f"  Missing apo predictions ({len(s['missing_apo'])}):")
                f.write(f"  Missing apo predictions ({len(s['missing_apo'])}):\n")
                for seq in s['missing_apo'][:10]:  # Show first 10
                    print(f"    - {seq}")
                    f.write(f"    - {seq}\n")
                if len(s['missing_apo']) > 10:
                    print(f"    ... and {len(s['missing_apo']) - 10} more")
                    f.write(f"    ... and {len(s['missing_apo']) - 10} more\n")
            print()
            f.write("\n")
    
    print(f"Detailed list saved to: {missing_file}")
    print()

# Summary
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()

total_missing_ligandmpnn = len(missing_ligandmpnn)
total_missing_af3 = sum(len(s['missing_bound']) + len(s['missing_apo']) for s in missing_af3_predictions)

print(f"Total input PDBs: {len(input_pdbs)}")
print(f"Missing LigandMPNN runs: {total_missing_ligandmpnn}")
print()

if total_missing_ligandmpnn > 0:
    print("⚠ ACTION REQUIRED: Run LigandMPNN for missing structures")
    print()

total_expected_af3 = sum(s['expected_bound'] + s['expected_apo'] for s in af3_summary)
total_completed_af3 = sum(s['completed_bound'] + s['completed_apo'] for s in af3_summary)

print(f"Total AlphaFold3 predictions expected: {total_expected_af3}")
print(f"Total AlphaFold3 predictions completed: {total_completed_af3}")
print(f"Missing AlphaFold3 predictions: {total_missing_af3}")
print()

if total_missing_af3 > 0:
    print("⚠ ACTION REQUIRED: Run AlphaFold3 for missing predictions")
    print()

# Create task lists
if missing_ligandmpnn:
    task_file = BASE_DIR / "run_missing_ligandmpnn.sh"
    with open(task_file, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("# Run missing LigandMPNN predictions\n\n")
        f.write("set -e\n\n")
        
        for output_name in missing_ligandmpnn:
            input_name = output_name.replace("output_", "")
            f.write(f"# {input_name}\n")
            f.write(f"echo 'Running LigandMPNN for {input_name}...'\n")
            f.write(f"# TODO: Update config and run pipeline for {input_name}\n\n")
    
    print(f"Script template created: {task_file}")
    print()

print("=" * 80)

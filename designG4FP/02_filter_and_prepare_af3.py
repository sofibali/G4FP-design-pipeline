#!/usr/bin/env python3
"""
Step 2: Filter top sequences by confidence and prepare AlphaFold3 inputs

Supports two filtering modes:
  - Top-N by combined confidence (default, fast)
  - Pareto-optimal filtering across overall + ligand confidence (--pareto)

Usage:
  python 02_filter_and_prepare_af3.py                           # from config
  python 02_filter_and_prepare_af3.py --output-dir DIR --input-pdb PDB  # CLI
  python 02_filter_and_prepare_af3.py --pareto --top-n 200      # Pareto mode
"""

import sys
import json
import yaml
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

def load_config(config_file: str = "pipeline_config.yaml") -> dict:
    """Load pipeline configuration"""
    with open(config_file) as f:
        return yaml.safe_load(f)

def parse_fasta_headers(fasta_file: Path) -> pd.DataFrame:
    """Parse FASTA headers to extract confidence scores"""
    sequences = []
    
    with open(fasta_file) as f:
        lines = f.readlines()
    
    for i in range(0, len(lines), 2):
        if i + 1 < len(lines):
            header = lines[i].strip('>\n')
            seq = lines[i + 1].strip()
            
            # Parse header: ">name, T=X, seed=Y, overall_confidence=Z, ..."
            header_parts = [x.strip() for x in header.split(',')]
            
            data = {'header': header, 'sequence': seq}
            for part in header_parts:
                if '=' in part:
                    key, val = part.split('=', 1)
                    key = key.strip()
                    try:
                        data[key] = float(val)
                    except:
                        data[key] = val
            
            sequences.append(data)
    
    return pd.DataFrame(sequences)

def filter_top_sequences(ligandmpnn_dir: Path, top_n: int = 100) -> pd.DataFrame:
    """Filter top N sequences by combined ligand and overall confidence"""
    
    print(f"Loading sequences from {ligandmpnn_dir}...")
    
    # Find FASTA files
    fasta_files = list(ligandmpnn_dir.glob("seqs/*.fa"))
    
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in {ligandmpnn_dir}/seqs/")
    
    print(f"  Found {len(fasta_files)} FASTA file(s)")
    
    # Parse all sequences
    all_seqs = []
    for fasta_file in fasta_files:
        df = parse_fasta_headers(fasta_file)
        all_seqs.append(df)
        print(f"    {fasta_file.name}: {len(df)} sequences")
    
    df_all = pd.concat(all_seqs, ignore_index=True)
    print(f"\n  Total sequences: {len(df_all)}")
    
    # Calculate combined confidence score
    if 'overall_confidence' in df_all.columns and 'ligand_confidence' in df_all.columns:
        # Combined score: average of both confidences
        df_all['combined_confidence'] = (df_all['overall_confidence'] + df_all['ligand_confidence']) / 2
        print(f"\n  Calculating combined confidence (average of overall + ligand)")
        
        # Sort by combined confidence
        df_all = df_all.sort_values('combined_confidence', ascending=False)
        
        num_to_keep = min(top_n, len(df_all))
        df_top = df_all.head(num_to_keep).reset_index(drop=True)
        
        print(f"\n  Filtering top {num_to_keep} sequences by combined confidence")
        print(f"  Combined confidence range: {df_top['combined_confidence'].min():.4f} - {df_top['combined_confidence'].max():.4f}")
        print(f"  Overall confidence range: {df_top['overall_confidence'].min():.4f} - {df_top['overall_confidence'].max():.4f}")
        print(f"  Ligand confidence range: {df_top['ligand_confidence'].min():.4f} - {df_top['ligand_confidence'].max():.4f}")
    else:
        print(f"  Warning: confidence columns not found, taking first {top_n} sequences")
        num_to_keep = min(top_n, len(df_all))
        df_top = df_all.head(num_to_keep).reset_index(drop=True)
    
    # Add sequence IDs
    df_top['seq_id'] = range(len(df_top))
    
    return df_top


def pareto_filter_sequences(ligandmpnn_dir: Path, top_n: int = 200) -> pd.DataFrame:
    """
    Select sequences on the Pareto frontier of overall_confidence vs ligand_confidence.

    Pareto-optimal = no other sequence beats it on BOTH metrics.
    Returns up to top_n sequences: all Pareto-optimal first, then fill by combined score.
    """
    print(f"Loading sequences from {ligandmpnn_dir} (Pareto mode)...")

    fasta_files = list(ligandmpnn_dir.glob("seqs/*.fa"))
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in {ligandmpnn_dir}/seqs/")

    all_seqs = []
    for fasta_file in fasta_files:
        df = parse_fasta_headers(fasta_file)
        all_seqs.append(df)
    df_all = pd.concat(all_seqs, ignore_index=True)
    print(f"  Total sequences: {len(df_all)}")

    if 'overall_confidence' not in df_all.columns or 'ligand_confidence' not in df_all.columns:
        print("  WARNING: confidence columns missing, falling back to top-N")
        return filter_top_sequences(ligandmpnn_dir, top_n)

    # Find Pareto frontier
    vals = df_all[['overall_confidence', 'ligand_confidence']].values
    n = len(vals)
    is_pareto = np.ones(n, dtype=bool)
    for i in range(n):
        if not is_pareto[i]:
            continue
        # Vectorized dominance check: does any other point dominate i?
        others = vals[is_pareto]
        other_idx = np.where(is_pareto)[0]
        dominated_by = (np.all(others >= vals[i], axis=1) &
                        np.any(others > vals[i], axis=1))
        # Exclude self
        self_mask = other_idx == i
        dominated_by[self_mask] = False
        if np.any(dominated_by):
            is_pareto[i] = False

    n_pareto = int(np.sum(is_pareto))
    pareto_df = df_all[is_pareto].copy()
    print(f"  Pareto-optimal sequences: {n_pareto}")

    # Combined score for ranking
    df_all['combined_confidence'] = (df_all['overall_confidence'] + df_all['ligand_confidence']) / 2
    pareto_df['combined_confidence'] = (pareto_df['overall_confidence'] + pareto_df['ligand_confidence']) / 2
    pareto_df = pareto_df.sort_values('combined_confidence', ascending=False)

    if len(pareto_df) >= top_n:
        df_top = pareto_df.head(top_n).reset_index(drop=True)
    else:
        non_pareto = df_all[~is_pareto].sort_values(
            (df_all[~is_pareto]['overall_confidence'] + df_all[~is_pareto]['ligand_confidence']) / 2,
            ascending=False
        )
        fill = non_pareto.head(top_n - len(pareto_df))
        df_top = pd.concat([pareto_df, fill], ignore_index=True)

    print(f"  Selected {len(df_top)} sequences")
    df_top['seq_id'] = range(len(df_top))
    return df_top


def extract_ions_from_pdb(pdb_path: str) -> List[Dict]:
    """Extract ion information from PDB file"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('pdb', pdb_path)
    
    k_count = 0
    
    for model in structure:
        for chain in model:
            for residue in chain:
                res_name = residue.get_resname().strip()
                
                # Check for potassium ions
                if res_name == 'K':
                    # Get coordinates from the K atom
                    if 'K' in residue:
                        k_count += 1
    
    # Return single entry with total count
    ions = []
    if k_count > 0:
        ions.append({
            'type': 'K',
            'count': k_count
        })
    
    return ions

def extract_dna_sequence(pdb_path: str) -> str:
    """Extract DNA sequence from PDB file"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('pdb', pdb_path)
    
    dna_seq = ""
    
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            
            # Skip protein chain (assuming A is protein)
            if chain_id == 'A':
                continue
            
            nucleotides = []
            for residue in chain:
                res_name = residue.get_resname().strip()
                
                # Map 3-letter codes to 1-letter
                nuc_map = {
                    'DA': 'A', 'A': 'A',
                    'DG': 'G', 'G': 'G',
                    'DC': 'C', 'C': 'C',
                    'DT': 'T', 'T': 'T', 'U': 'U',
                }
                
                if res_name in nuc_map:
                    nucleotides.append(nuc_map[res_name])
            
            if nucleotides:
                dna_seq = ''.join(nucleotides)
                print(f"    Extracted DNA from chain {chain_id}: {dna_seq}")
                break
    
    return dna_seq

def create_af3_input_json(protein_seq: str, dna_seq: str, seq_id: int, 
                          modelseeds: List[int], ions: List[Dict] = None, 
                          include_dna: bool = True) -> Dict:
    """Create AlphaFold3 input JSON"""
    
    sequences = [{"protein": {"id": "1", "sequence": protein_seq}}]

    current_id = 2
    if include_dna and dna_seq:
        sequences.append({"dna": {"id": str(current_id), "sequence": dna_seq}})
        current_id += 1

    # Add ions if present - must wrap in 'ligand' key for AlphaFold3 format
    if ions:
        for ion in ions:
            sequences.append({
                "ligand": {
                    "id": str(current_id),
                    "type": ion['type'],
                    "count": ion['count']
                }
            })
            current_id += 1
    
    af3_input = {
        "dialect": "alphafold3",
        "version": 1,
        "name": f"design_{seq_id}{'_bound' if include_dna else '_apo'}",
        "modelSeeds": modelseeds,
        "sequences": sequences
    }
    
    return af3_input

def prepare_alphafold3_inputs(df_top: pd.DataFrame, config: dict, pdb_input: str) -> Tuple[Path, Path]:
    """Prepare AlphaFold3 input JSON files"""
    
    output_dir = Path(config['output_dir'])
    af3_bound_dir = output_dir / "02_alphafold3_inputs_bound"
    af3_apo_dir = output_dir / "02_alphafold3_inputs_apo"
    
    af3_bound_dir.mkdir(parents=True, exist_ok=True)
    af3_apo_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract DNA sequence from original PDB
    print("\nExtracting DNA sequence from input PDB...")
    dna_seq = extract_dna_sequence(pdb_input)
    
    if not dna_seq:
        print("  ⚠ No DNA found in PDB. Will create protein-only inputs.")
    
    # Extract ions from original PDB
    print("\nExtracting ions from input PDB...")
    ions = extract_ions_from_pdb(pdb_input)
    if ions:
        print(f"    Found {len(ions)} potassium ion(s)")
    else:
        print("    No ions found in PDB")
    
    print(f"\nPreparing AlphaFold3 input JSONs...")
    
    modelseeds = config['alphafold3']['modelseeds']
    create_ligand_free = config['alphafold3']['create_ligand_free']
    
    for idx, row in df_top.iterrows():
        prot_seq = row['sequence']
        seq_id = row['seq_id']
        
        # Create bound (with DNA + ions) input
        af3_bound = create_af3_input_json(prot_seq, dna_seq, seq_id, modelseeds, ions=ions, include_dna=True)
        bound_json_path = af3_bound_dir / f"design_{seq_id:04d}_bound.json"
        with open(bound_json_path, 'w') as f:
            json.dump(af3_bound, f, indent=2)
        
        # Create apo (protein-only) input if requested
        if create_ligand_free:
            af3_apo = create_af3_input_json(prot_seq, "", seq_id, modelseeds, ions=None, include_dna=False)
            apo_json_path = af3_apo_dir / f"design_{seq_id:04d}_apo.json"
            with open(apo_json_path, 'w') as f:
                json.dump(af3_apo, f, indent=2)
    
    print(f"  ✓ Created {len(df_top)} bound (protein+DNA) input files in {af3_bound_dir.name}/")
    if create_ligand_free:
        print(f"  ✓ Created {len(df_top)} apo (protein-only) input files in {af3_apo_dir.name}/")
    
    return af3_bound_dir, af3_apo_dir

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Filter LigandMPNN sequences and prepare AlphaFold3 inputs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python 02_filter_and_prepare_af3.py                                  # from config
  python 02_filter_and_prepare_af3.py my_config.yaml                   # custom config
  python 02_filter_and_prepare_af3.py --output-dir DIR --input-pdb PDB # direct CLI
  python 02_filter_and_prepare_af3.py --pareto --top-n 200             # Pareto mode
        """
    )
    parser.add_argument('config', nargs='?', default='pipeline_config.yaml',
                       help='Path to configuration file (default: pipeline_config.yaml)')
    parser.add_argument('--output-dir', type=str, default=None,
                       help='Override output directory (bypasses config)')
    parser.add_argument('--input-pdb', type=str, default=None,
                       help='Override input PDB path (bypasses config)')
    parser.add_argument('--top-n', type=int, default=None,
                       help='Override number of sequences to keep')
    parser.add_argument('--pareto', action='store_true',
                       help='Use Pareto filtering (overall + ligand confidence)')
    args = parser.parse_args()

    print("="*80)
    print("STEP 2: Filter Sequences and Prepare AlphaFold3 Inputs")
    print("="*80)
    print()

    config = load_config(args.config)

    # CLI overrides take precedence over config
    if args.output_dir:
        config['output_dir'] = args.output_dir
    if args.input_pdb:
        config['input_pdb'] = args.input_pdb
    if args.top_n:
        config['filtering']['top_n'] = args.top_n

    # Set up paths
    output_dir = Path(config['output_dir'])
    ligandmpnn_dir = output_dir / "01_ligandmpnn"

    # Check if LigandMPNN output exists
    if not ligandmpnn_dir.exists():
        print(f"ERROR: LigandMPNN output not found: {ligandmpnn_dir}")
        print("   Please run 01_run_ligandmpnn.sh first")
        sys.exit(1)

    # Filter sequences
    top_n = config['filtering'].get('top_n', 100)

    if args.pareto:
        print(f"  Mode: Pareto filtering (overall + ligand confidence)")
        df_top = pareto_filter_sequences(ligandmpnn_dir, top_n)
    else:
        print(f"  Mode: Top-N by combined confidence")
        df_top = filter_top_sequences(ligandmpnn_dir, top_n)

    # Save filtered sequences
    filtered_csv = output_dir / "02_filtered_sequences.csv"
    df_top.to_csv(filtered_csv, index=False)
    print(f"\n  Saved filtered sequences to {filtered_csv}")

    # Save FASTA
    filtered_fasta = output_dir / "02_filtered_sequences.fa"
    with open(filtered_fasta, 'w') as f:
        for _, row in df_top.iterrows():
            f.write(f">{row['header']}\n{row['sequence']}\n")
    print(f"  Saved FASTA to {filtered_fasta}")

    # Prepare AlphaFold3 inputs
    pdb_input = config['input_pdb']
    af3_bound_dir, af3_apo_dir = prepare_alphafold3_inputs(df_top, config, pdb_input)

    mode = "Pareto" if args.pareto else "combined confidence"
    print()
    print("="*80)
    print("Filtering and AF3 input preparation complete!")
    print("="*80)
    print(f"\nFiltered sequences: {len(df_top)} (top {top_n} by {mode})")
    print(f"Bound inputs: {af3_bound_dir}/")
    print(f"Apo inputs: {af3_apo_dir}/")
    print()
    print("Next step: Run 03_run_alphafold3.sh or validate with utils/validate_af3_jsons.py")
    print()

if __name__ == "__main__":
    main()

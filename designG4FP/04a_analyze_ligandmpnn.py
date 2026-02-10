#!/usr/bin/env python3
"""
Step 4a: Analyze LigandMPNN sequence diversity and design confidence metrics
Usage: python 04a_analyze_ligandmpnn.py [config_file]
"""

import sys
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

def load_config(config_file: str = "pipeline_config.yaml") -> dict:
    """Load pipeline configuration"""
    with open(config_file) as f:
        return yaml.safe_load(f)

def parse_fasta_file(fasta_file: Path) -> List[Dict]:
    """Parse FASTA file with headers"""
    sequences = []
    with open(fasta_file) as f:
        lines = f.readlines()
    
    for i in range(0, len(lines), 2):
        if i + 1 < len(lines):
            header = lines[i].strip('>\n')
            seq = lines[i + 1].strip()
            
            # Parse header metrics
            data = {'header': header, 'sequence': seq, 'length': len(seq)}
            for part in header.split(','):
                if '=' in part:
                    key, val = part.split('=', 1)
                    key = key.strip()
                    try:
                        data[key] = float(val)
                    except:
                        data[key] = val
            
            sequences.append(data)
    
    return sequences

def calculate_sequence_identity(seq1: str, seq2: str) -> float:
    """Calculate sequence identity between two sequences"""
    if len(seq1) != len(seq2):
        return 0.0
    matches = sum(a == b for a, b in zip(seq1, seq2))
    return 100.0 * matches / len(seq1)

def calculate_positional_entropy(sequences: List[str]) -> np.ndarray:
    """Calculate Shannon entropy at each position"""
    if not sequences:
        return np.array([])
    
    seq_len = len(sequences[0])
    entropy = np.zeros(seq_len)
    
    for pos in range(seq_len):
        # Get amino acids at this position
        amino_acids = [seq[pos] for seq in sequences if pos < len(seq)]
        
        # Calculate frequency
        counts = Counter(amino_acids)
        total = len(amino_acids)
        
        # Calculate Shannon entropy
        pos_entropy = 0.0
        for count in counts.values():
            if count > 0:
                p = count / total
                pos_entropy -= p * np.log2(p)
        
        entropy[pos] = pos_entropy
    
    return entropy

def analyze_sequence_diversity(df: pd.DataFrame, template_seq: str = None) -> Dict:
    """Analyze sequence diversity metrics"""
    
    sequences = df['sequence'].tolist()
    
    # Basic statistics
    unique_seqs = len(set(sequences))
    diversity_ratio = unique_seqs / len(sequences)
    
    # Sequence identity matrix (sample if too many sequences)
    sample_size = min(100, len(sequences))
    sample_indices = np.random.choice(len(sequences), sample_size, replace=False)
    sample_seqs = [sequences[i] for i in sample_indices]
    
    identity_matrix = np.zeros((sample_size, sample_size))
    for i, seq1 in enumerate(sample_seqs):
        for j, seq2 in enumerate(sample_seqs):
            identity_matrix[i, j] = calculate_sequence_identity(seq1, seq2)
    
    mean_identity = np.mean(identity_matrix[np.triu_indices(sample_size, k=1)])
    
    # Positional entropy
    entropy = calculate_positional_entropy(sequences)
    
    # Identity to template
    template_identity = []
    if template_seq:
        for seq in sequences:
            identity = calculate_sequence_identity(seq, template_seq)
            template_identity.append(identity)
    
    results = {
        'total_sequences': len(sequences),
        'unique_sequences': unique_seqs,
        'diversity_ratio': diversity_ratio,
        'mean_pairwise_identity': mean_identity,
        'positional_entropy': entropy,
        'template_identity': template_identity,
        'identity_matrix': identity_matrix,
    }
    
    return results

def create_visualizations(df: pd.DataFrame, diversity_results: Dict, 
                         output_dir: Path, template_seq: str = None):
    """Create visualization plots"""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_dir = output_dir / "plot_data_csvs"
    csv_dir.mkdir(exist_ok=True)
    
    # Set style
    sns.set_style("whitegrid")
    
    # 1. Confidence score distribution
    if 'overall_confidence' in df.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(df['overall_confidence'], bins=50, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Overall Confidence Score', fontsize=12)
        ax.set_ylabel('Number of Designs', fontsize=12)
        ax.set_title('LigandMPNN Confidence Score Distribution', fontsize=14, fontweight='bold')
        ax.axvline(df['overall_confidence'].median(), color='red', linestyle='--', 
                   label=f"Median: {df['overall_confidence'].median():.3f}")
        ax.legend()
        plt.tight_layout()
        plt.savefig(output_dir / "01_confidence_distribution.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export CSV
        df[['overall_confidence']].to_csv(csv_dir / "01_confidence_distribution.csv", index=False)
    
    # 2. Ligand confidence vs overall confidence
    if 'ligand_confidence' in df.columns and 'overall_confidence' in df.columns:
        fig, ax = plt.subplots(figsize=(10, 8))
        scatter = ax.scatter(df['overall_confidence'], df['ligand_confidence'], 
                           alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
        ax.set_xlabel('Overall Confidence', fontsize=12)
        ax.set_ylabel('Ligand Confidence', fontsize=12)
        ax.set_title('Overall vs Ligand Confidence Scores', fontsize=14, fontweight='bold')
        
        # Correlation
        corr = df[['overall_confidence', 'ligand_confidence']].corr().iloc[0, 1]
        ax.text(0.05, 0.95, f'Correlation: {corr:.3f}', transform=ax.transAxes,
                fontsize=11, verticalalignment='top', bbox=dict(boxstyle='round', 
                facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        plt.savefig(output_dir / "02_confidence_correlation.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export CSV
        df[['overall_confidence', 'ligand_confidence']].to_csv(
            csv_dir / "02_confidence_correlation.csv", index=False)
    
    # 3. Positional entropy (sequence diversity)
    if len(diversity_results['positional_entropy']) > 0:
        fig, ax = plt.subplots(figsize=(14, 6))
        positions = np.arange(1, len(diversity_results['positional_entropy']) + 1)
        ax.plot(positions, diversity_results['positional_entropy'], linewidth=1.5)
        ax.fill_between(positions, diversity_results['positional_entropy'], alpha=0.3)
        ax.set_xlabel('Residue Position', fontsize=12)
        ax.set_ylabel('Shannon Entropy (bits)', fontsize=12)
        ax.set_title('Sequence Diversity per Position', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_dir / "03_positional_entropy.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export CSV
        entropy_df = pd.DataFrame({
            'position': positions,
            'entropy': diversity_results['positional_entropy']
        })
        entropy_df.to_csv(csv_dir / "03_positional_entropy.csv", index=False)
    
    # 4. Identity to template distribution
    if diversity_results['template_identity']:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(diversity_results['template_identity'], bins=50, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Sequence Identity to Template (%)', fontsize=12)
        ax.set_ylabel('Number of Designs', fontsize=12)
        ax.set_title('Sequence Identity to Template Structure', fontsize=14, fontweight='bold')
        ax.axvline(np.mean(diversity_results['template_identity']), color='red', 
                   linestyle='--', label=f"Mean: {np.mean(diversity_results['template_identity']):.1f}%")
        ax.legend()
        plt.tight_layout()
        plt.savefig(output_dir / "04_template_identity.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export CSV
        identity_df = pd.DataFrame({'template_identity': diversity_results['template_identity']})
        identity_df.to_csv(csv_dir / "04_template_identity.csv", index=False)
    
    # 5. Sequence identity heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(diversity_results['identity_matrix'], cmap='viridis', aspect='auto')
    ax.set_xlabel('Design Index (sample)', fontsize=12)
    ax.set_ylabel('Design Index (sample)', fontsize=12)
    ax.set_title('Pairwise Sequence Identity Matrix', fontsize=14, fontweight='bold')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Sequence Identity (%)', fontsize=12)
    plt.tight_layout()
    plt.savefig(output_dir / "05_identity_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Export CSV
    identity_matrix_df = pd.DataFrame(diversity_results['identity_matrix'])
    identity_matrix_df.to_csv(csv_dir / "05_identity_matrix.csv", index=False)
    
    print(f"  ✓ Generated 5 visualization plots")
    print(f"  ✓ Exported CSV data files to {csv_dir.name}/")

def load_template_sequence(pdb_path: str) -> str:
    """Load template sequence from PDB"""
    try:
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('template', pdb_path)
        
        # Get chain A sequence
        for model in structure:
            for chain in model:
                if chain.get_id() == 'A':
                    seq = ""
                    for residue in chain:
                        if residue.get_id()[0] == ' ':  # Standard residue
                            from Bio.PDB.Polypeptide import three_to_one
                            try:
                                seq += three_to_one(residue.get_resname())
                            except:
                                seq += 'X'
                    return seq
    except Exception as e:
        print(f"  Warning: Could not load template sequence: {e}")
    
    return None

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Analyze LigandMPNN sequence diversity and design confidence',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python 04a_analyze_ligandmpnn.py
  python 04a_analyze_ligandmpnn.py my_config.yaml
        """
    )
    parser.add_argument('config', nargs='?', default='pipeline_config.yaml',
                       help='Path to configuration file (default: pipeline_config.yaml)')
    args = parser.parse_args()
    
    config_file = args.config
    
    print("="*80)
    print("STEP 4a: Analyze LigandMPNN Sequences")
    print("="*80)
    print()
    
    config = load_config(config_file)
    
    # Set up paths
    output_dir = Path(config['output_dir'])
    ligandmpnn_dir = output_dir / "01_ligandmpnn"
    analysis_output = output_dir / "04a_ligandmpnn_analysis"
    
    # Check if LigandMPNN output exists
    if not ligandmpnn_dir.exists():
        print(f"❌ ERROR: LigandMPNN output not found: {ligandmpnn_dir}")
        sys.exit(1)
    
    # Load all sequences
    print("Loading LigandMPNN sequences...")
    fasta_files = list(ligandmpnn_dir.glob("seqs/*.fa"))
    
    if not fasta_files:
        print(f"❌ ERROR: No FASTA files found in {ligandmpnn_dir}/seqs/")
        sys.exit(1)
    
    all_sequences = []
    for fasta_file in fasta_files:
        seqs = parse_fasta_file(fasta_file)
        all_sequences.extend(seqs)
    
    df = pd.DataFrame(all_sequences)
    print(f"  Loaded {len(df)} sequences")
    
    # Load template sequence
    template_seq = None
    if Path(config['input_pdb']).exists():
        print("\nLoading template sequence...")
        template_seq = load_template_sequence(config['input_pdb'])
        if template_seq:
            print(f"  Template length: {len(template_seq)} residues")
    
    # Analyze diversity
    print("\nAnalyzing sequence diversity...")
    diversity_results = analyze_sequence_diversity(df, template_seq)
    
    print(f"  Total sequences: {diversity_results['total_sequences']}")
    print(f"  Unique sequences: {diversity_results['unique_sequences']}")
    print(f"  Diversity ratio: {diversity_results['diversity_ratio']:.3f}")
    print(f"  Mean pairwise identity: {diversity_results['mean_pairwise_identity']:.1f}%")
    
    if diversity_results['template_identity']:
        print(f"  Mean identity to template: {np.mean(diversity_results['template_identity']):.1f}%")
    
    # Create visualizations
    print("\nGenerating visualizations...")
    create_visualizations(df, diversity_results, analysis_output, template_seq)
    
    # Save summary statistics
    summary = {
        'total_sequences': diversity_results['total_sequences'],
        'unique_sequences': diversity_results['unique_sequences'],
        'diversity_ratio': diversity_results['diversity_ratio'],
        'mean_pairwise_identity': diversity_results['mean_pairwise_identity'],
    }
    
    if 'overall_confidence' in df.columns:
        summary.update({
            'mean_confidence': df['overall_confidence'].mean(),
            'median_confidence': df['overall_confidence'].median(),
            'std_confidence': df['overall_confidence'].std(),
        })
    
    if diversity_results['template_identity']:
        summary['mean_template_identity'] = np.mean(diversity_results['template_identity'])
    
    summary_df = pd.DataFrame([summary])
    summary_csv = analysis_output / "ligandmpnn_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    
    print(f"\n  ✓ Saved summary to {summary_csv}")
    
    print()
    print("="*80)
    print("✓ LigandMPNN analysis complete!")
    print("="*80)
    print(f"\nResults saved to: {analysis_output}")
    print(f"  - Visualization plots (PNG)")
    print(f"  - CSV data files for Prism/GraphPad")
    print(f"  - Summary statistics")
    print()

if __name__ == "__main__":
    main()

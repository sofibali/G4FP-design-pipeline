#!/usr/bin/env python3
"""
Step 1b: Analyze LigandMPNN output before filtering
Usage: python 01b_analyze_pre_filter.py [config_file]

Generates:
- Overall confidence distribution
- Ligand confidence distribution  
- Sequence diversity per residue position
- Ligand vs overall confidence scatter (colored by sequence identity)
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

def calculate_positional_entropy(sequences: List[str]) -> np.ndarray:
    """Calculate Shannon entropy at each position"""
    if not sequences:
        return np.array([])
    
    seq_len = len(sequences[0])
    entropy = np.zeros(seq_len)
    
    for pos in range(seq_len):
        amino_acids = [seq[pos] for seq in sequences if pos < len(seq)]
        counts = Counter(amino_acids)
        total = len(amino_acids)
        
        pos_entropy = 0.0
        for count in counts.values():
            if count > 0:
                p = count / total
                pos_entropy -= p * np.log2(p)
        
        entropy[pos] = pos_entropy
    
    return entropy

def create_analysis_plots(df: pd.DataFrame, output_dir: Path):
    """Create comprehensive pre-filter analysis plots"""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    sns.set_style("whitegrid")
    
    print(f"\nGenerating analysis plots...")
    print(f"  Total sequences: {len(df)}")
    
    # ========================================================================
    # Plot 1: Overall Confidence Distribution
    # ========================================================================
    if 'overall_confidence' in df.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(df['overall_confidence'], bins=50, edgecolor='black', alpha=0.7, color='steelblue')
        ax.set_xlabel('Overall Confidence', fontsize=12)
        ax.set_ylabel('Number of Designs', fontsize=12)
        ax.set_title('LigandMPNN Overall Confidence Distribution\n(All Sequences Before Filtering)', 
                     fontsize=14, fontweight='bold')
        
        # Add statistics
        mean_val = df['overall_confidence'].mean()
        median_val = df['overall_confidence'].median()
        ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, 
                   label=f'Mean: {mean_val:.3f}')
        ax.axvline(median_val, color='orange', linestyle='--', linewidth=2,
                   label=f'Median: {median_val:.3f}')
        ax.legend()
        
        plt.tight_layout()
        plt.savefig(output_dir / "01_overall_confidence_distribution.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Plot 1: Overall confidence distribution")
    
    # ========================================================================
    # Plot 2: Ligand Confidence Distribution
    # ========================================================================
    if 'ligand_confidence' in df.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(df['ligand_confidence'], bins=50, edgecolor='black', alpha=0.7, color='coral')
        ax.set_xlabel('Ligand Confidence', fontsize=12)
        ax.set_ylabel('Number of Designs', fontsize=12)
        ax.set_title('LigandMPNN Ligand Confidence Distribution\n(All Sequences Before Filtering)', 
                     fontsize=14, fontweight='bold')
        
        # Add statistics
        mean_val = df['ligand_confidence'].mean()
        median_val = df['ligand_confidence'].median()
        ax.axvline(mean_val, color='darkred', linestyle='--', linewidth=2,
                   label=f'Mean: {mean_val:.3f}')
        ax.axvline(median_val, color='orange', linestyle='--', linewidth=2,
                   label=f'Median: {median_val:.3f}')
        ax.legend()
        
        plt.tight_layout()
        plt.savefig(output_dir / "02_ligand_confidence_distribution.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Plot 2: Ligand confidence distribution")
    
    # ========================================================================
    # Plot 3: Sequence Diversity per Residue Position
    # ========================================================================
    sequences = df['sequence'].tolist()
    entropy = calculate_positional_entropy(sequences)
    
    if len(entropy) > 0:
        fig, ax = plt.subplots(figsize=(14, 6))
        positions = np.arange(1, len(entropy) + 1)
        
        ax.plot(positions, entropy, linewidth=1.5, color='darkgreen', alpha=0.8)
        ax.fill_between(positions, entropy, alpha=0.3, color='green')
        
        ax.set_xlabel('Residue Position', fontsize=12)
        ax.set_ylabel('Shannon Entropy (bits)', fontsize=12)
        ax.set_title('Sequence Diversity per Residue Position\n(Higher entropy = more diverse)', 
                     fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Add mean entropy line
        mean_entropy = entropy.mean()
        ax.axhline(mean_entropy, color='red', linestyle='--', linewidth=1.5,
                   label=f'Mean entropy: {mean_entropy:.3f}')
        ax.legend()
        
        plt.tight_layout()
        plt.savefig(output_dir / "03_sequence_diversity_per_position.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Plot 3: Sequence diversity per position")
    
    # ========================================================================
    # Plot 4: Ligand Confidence vs Overall Confidence (colored by seq identity)
    # ========================================================================
    if 'ligand_confidence' in df.columns and 'overall_confidence' in df.columns:
        # Filter out nan values for seq_rec
        plot_df = df.dropna(subset=['seq_rec']) if 'seq_rec' in df.columns else df
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        if 'seq_rec' in plot_df.columns and len(plot_df) > 0:
            scatter = ax.scatter(plot_df['overall_confidence'], 
                               plot_df['ligand_confidence'],
                               c=plot_df['seq_rec'], 
                               cmap='viridis', 
                               s=30, 
                               alpha=0.6,
                               edgecolors='black', 
                               linewidth=0.3)
            
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('Sequence Identity to Template', fontsize=12)
        else:
            ax.scatter(df['overall_confidence'], 
                      df['ligand_confidence'],
                      s=30, 
                      alpha=0.6,
                      edgecolors='black', 
                      linewidth=0.3)
        
        ax.set_xlabel('Overall Confidence', fontsize=12)
        ax.set_ylabel('Ligand Confidence', fontsize=12)
        ax.set_title('Ligand vs Overall Confidence\n(Colored by Sequence Identity)', 
                     fontsize=14, fontweight='bold')
        
        # Add correlation
        corr = df[['overall_confidence', 'ligand_confidence']].corr().iloc[0, 1]
        ax.text(0.05, 0.95, f'Correlation: {corr:.3f}', 
                transform=ax.transAxes,
                fontsize=11, 
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_dir / "04_ligand_vs_overall_confidence.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Plot 4: Ligand vs overall confidence scatter")
    
    print(f"\n  All plots saved to: {output_dir}")

def generate_summary_stats(df: pd.DataFrame, output_dir: Path):
    """Generate summary statistics"""
    
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    
    stats = {
        'total_sequences': len(df),
        'unique_sequences': len(df['sequence'].unique()),
    }
    
    if 'overall_confidence' in df.columns:
        stats.update({
            'overall_confidence_mean': df['overall_confidence'].mean(),
            'overall_confidence_median': df['overall_confidence'].median(),
            'overall_confidence_std': df['overall_confidence'].std(),
            'overall_confidence_min': df['overall_confidence'].min(),
            'overall_confidence_max': df['overall_confidence'].max(),
        })
    
    if 'ligand_confidence' in df.columns:
        stats.update({
            'ligand_confidence_mean': df['ligand_confidence'].mean(),
            'ligand_confidence_median': df['ligand_confidence'].median(),
            'ligand_confidence_std': df['ligand_confidence'].std(),
            'ligand_confidence_min': df['ligand_confidence'].min(),
            'ligand_confidence_max': df['ligand_confidence'].max(),
        })
    
    if 'seq_rec' in df.columns:
        valid_seq_rec = df['seq_rec'].dropna()
        if len(valid_seq_rec) > 0:
            stats.update({
                'sequence_identity_mean': valid_seq_rec.mean(),
                'sequence_identity_median': valid_seq_rec.median(),
                'sequence_identity_std': valid_seq_rec.std(),
                'sequence_identity_min': valid_seq_rec.min(),
                'sequence_identity_max': valid_seq_rec.max(),
            })
    
    # Save summary
    summary_df = pd.DataFrame([stats])
    summary_csv = output_dir / "pre_filter_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    
    return stats

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Analyze LigandMPNN output before filtering',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python 01b_analyze_pre_filter.py
  python 01b_analyze_pre_filter.py my_config.yaml
        """
    )
    parser.add_argument('config', nargs='?', default='pipeline_config.yaml',
                       help='Path to configuration file (default: pipeline_config.yaml)')
    args = parser.parse_args()
    
    config_file = args.config
    
    print("="*80)
    print("STEP 1b: Pre-Filter Analysis of LigandMPNN Output")
    print("="*80)
    print()
    
    config = load_config(config_file)
    
    # Set up paths
    output_dir = Path(config['output_dir'])
    ligandmpnn_dir = output_dir / "01_ligandmpnn"
    analysis_output = output_dir / "01b_pre_filter_analysis"
    
    # Check if LigandMPNN output exists
    if not ligandmpnn_dir.exists():
        print(f"❌ ERROR: LigandMPNN output not found: {ligandmpnn_dir}")
        print("   Please run 01_run_ligandmpnn.sh first")
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
    print(f"  Loaded {len(df)} sequences from {len(fasta_files)} file(s)")
    
    # Generate summary statistics
    print("\nCalculating statistics...")
    stats = generate_summary_stats(df, analysis_output)
    
    print(f"\nSummary Statistics:")
    print(f"  Total sequences:        {stats['total_sequences']}")
    print(f"  Unique sequences:       {stats['unique_sequences']}")
    
    if 'overall_confidence_mean' in stats:
        print(f"\n  Overall Confidence:")
        print(f"    Mean:   {stats['overall_confidence_mean']:.4f}")
        print(f"    Median: {stats['overall_confidence_median']:.4f}")
        print(f"    Range:  {stats['overall_confidence_min']:.4f} - {stats['overall_confidence_max']:.4f}")
    
    if 'ligand_confidence_mean' in stats:
        print(f"\n  Ligand Confidence:")
        print(f"    Mean:   {stats['ligand_confidence_mean']:.4f}")
        print(f"    Median: {stats['ligand_confidence_median']:.4f}")
        print(f"    Range:  {stats['ligand_confidence_min']:.4f} - {stats['ligand_confidence_max']:.4f}")
    
    if 'sequence_identity_mean' in stats:
        print(f"\n  Sequence Identity:")
        print(f"    Mean:   {stats['sequence_identity_mean']:.4f} ({stats['sequence_identity_mean']*100:.1f}%)")
        print(f"    Median: {stats['sequence_identity_median']:.4f} ({stats['sequence_identity_median']*100:.1f}%)")
        print(f"    Range:  {stats['sequence_identity_min']:.4f} - {stats['sequence_identity_max']:.4f}")
    
    # Create visualizations
    create_analysis_plots(df, analysis_output)
    
    print()
    print("="*80)
    print("✓ Pre-filter analysis complete!")
    print("="*80)
    print(f"\nResults saved to: {analysis_output}")
    print(f"  - 4 visualization plots (PNG)")
    print(f"  - pre_filter_summary.csv (statistics)")
    print()
    print("Next step: Run 02_filter_and_prepare_af3.py to select top sequences")
    print()

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Analyze AlphaFold3 Structure Predictions for G4FP Designs

Comprehensive analysis of AF3 predictions including:
1. Overall confidence distributions (LigandMPNN scores)
2. AlphaFold3 quality metrics (pTM, ipTM, pLDDT, PAE)
3. RMSD analysis (global and chromophore-specific)
4. Confidence-pLDDT correlations
5. Per-residue analysis
6. Top design identification

Adapted for designG4FP pipeline structure where each design has:
- output_X/03_alphafold3_predictions_bound/design_NNNN_bound/
- output_X/03_alphafold3_predictions_apo/design_NNNN_apo/

Usage:
    python 05_analyze_af3_structures.py [--output-dir OUTPUT] [--state STATE]

Arguments:
    --output-dir: Specific output directory to analyze (e.g., output_G4FP_des1_cro_mod0)
    --state: Which predictions to analyze: bound, apo, or both (default: both)
"""

import os
import sys
import json
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from collections import defaultdict

# BioPython for structure analysis
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import PPBuilder

warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10


class AlphaFoldStructureAnalyzer:
    """Analyze AlphaFold3 predictions for G4FP designs"""
    
    def __init__(self, output_dir: Path, template_structure: str, 
                 chromophore_range: Tuple[int, int] = (175, 225)):
        """
        Initialize analyzer
        
        Args:
            output_dir: Path to specific output directory (e.g., output_G4FP_des1_cro_mod0)
            template_structure: Path to template/reference structure for RMSD
            chromophore_range: Tuple of (start, end) residue numbers for chromophore region
        """
        self.output_dir = Path(output_dir)
        self.template_path = Path(template_structure)
        self.chromophore_range = chromophore_range
        
        # Initialize parsers
        self.pdb_parser = PDBParser(QUIET=True)
        self.cif_parser = MMCIFParser(QUIET=True)
        
        # Load template structure
        self.template_structure = self._load_structure(self.template_path)
        self.template_ca_atoms = self._extract_ca_atoms(self.template_structure)
        
        print(f"Output directory: {self.output_dir}")
        print(f"Template structure: {self.template_path}")
        print(f"Template has {len(self.template_ca_atoms)} CA atoms")
        print(f"Chromophore region: {chromophore_range[0]}-{chromophore_range[1]}")
        
    def _load_structure(self, filepath: Path):
        """Load PDB or CIF structure"""
        if filepath.suffix == '.pdb':
            return self.pdb_parser.get_structure('structure', str(filepath))
        elif filepath.suffix == '.cif':
            return self.cif_parser.get_structure('structure', str(filepath))
        else:
            raise ValueError(f"Unknown file format: {filepath.suffix}")
    
    def _extract_ca_atoms(self, structure, chain_id='A'):
        """Extract CA atoms from structure"""
        ca_atoms = []
        for model in structure:
            for chain in model:
                if chain.get_id() == chain_id:
                    for residue in chain:
                        if residue.has_id('CA'):
                            ca_atoms.append(residue['CA'])
        return ca_atoms
    
    def calculate_rmsd(self, model_file: Path, chain_id='A') -> Dict[str, float]:
        """
        Calculate RMSD between model and template
        
        Returns:
            Dict with 'global_rmsd', 'chromophore_rmsd'
        """
        # Load model structure
        model_structure = self._load_structure(model_file)
        model_ca_atoms = self._extract_ca_atoms(model_structure, chain_id)
        
        # Check length match
        if len(model_ca_atoms) != len(self.template_ca_atoms):
            print(f"    Warning: Length mismatch - Model: {len(model_ca_atoms)}, Template: {len(self.template_ca_atoms)}")
            # Use minimum length
            min_len = min(len(model_ca_atoms), len(self.template_ca_atoms))
            model_ca_atoms = model_ca_atoms[:min_len]
            template_ca_atoms = self.template_ca_atoms[:min_len]
        else:
            template_ca_atoms = self.template_ca_atoms
        
        # Global RMSD
        sup = Superimposer()
        sup.set_atoms(template_ca_atoms, model_ca_atoms)
        global_rmsd = sup.rms
        
        # Chromophore RMSD
        chrom_start, chrom_end = self.chromophore_range
        if len(model_ca_atoms) >= chrom_end:
            chrom_template = template_ca_atoms[chrom_start-1:chrom_end]
            chrom_model = model_ca_atoms[chrom_start-1:chrom_end]
            
            sup_chrom = Superimposer()
            sup_chrom.set_atoms(chrom_template, chrom_model)
            chromophore_rmsd = sup_chrom.rms
        else:
            chromophore_rmsd = None
        
        return {
            'global_rmsd': global_rmsd,
            'chromophore_rmsd': chromophore_rmsd
        }
    
    def calculate_per_residue_rmsd(self, model_file: Path, chain_id='A') -> np.ndarray:
        """Calculate per-residue RMSD"""
        model_structure = self._load_structure(model_file)
        model_ca_atoms = self._extract_ca_atoms(model_structure, chain_id)
        
        # Match lengths
        min_len = min(len(model_ca_atoms), len(self.template_ca_atoms))
        model_ca_atoms = model_ca_atoms[:min_len]
        template_ca_atoms = self.template_ca_atoms[:min_len]
        
        # Superimpose structures
        sup = Superimposer()
        sup.set_atoms(template_ca_atoms, model_ca_atoms)
        sup.apply(model_structure)
        
        # Calculate per-residue RMSD
        per_res_rmsd = []
        for template_ca, model_ca in zip(template_ca_atoms, model_ca_atoms):
            distance = np.linalg.norm(
                np.array(template_ca.get_coord()) - np.array(model_ca.get_coord())
            )
            per_res_rmsd.append(distance)
        
        return np.array(per_res_rmsd)
    
    def load_confidence_scores(self, design_dir: Path) -> Dict:
        """Load confidence scores from AlphaFold3 summary JSON"""
        summary_file = design_dir / "summary_confidences_0.json"
        
        if not summary_file.exists():
            # Try alternative locations
            json_files = list(design_dir.glob("**/summary_confidences*.json"))
            if json_files:
                summary_file = json_files[0]
            else:
                raise FileNotFoundError(f"No confidence file found in {design_dir}")
        
        with open(summary_file, 'r') as f:
            confidence = json.load(f)
        
        return confidence
    
    def analyze_design(self, design_dir: Path, state: str = 'bound') -> Dict:
        """
        Analyze a single AlphaFold3 design
        
        Args:
            design_dir: Path to design directory (e.g., design_0000_bound/)
            state: 'bound' or 'apo'
            
        Returns:
            Dictionary with all metrics
        """
        design_name = design_dir.name
        seq_id = design_name.replace(f'design_', '').replace(f'_{state}', '')
        
        results = {
            'seq_id': int(seq_id),
            'design_name': design_name,
            'state': state,
        }
        
        # Find model file (should be design_NNNN_state_model.cif)
        model_files = list(design_dir.glob("*_model.cif"))
        if not model_files:
            print(f"  ⚠ No model file found for {design_name}")
            return results
        
        model_file = model_files[0]
        
        # Calculate RMSD
        print(f"  Analyzing {design_name}...")
        try:
            rmsd_data = self.calculate_rmsd(model_file)
            results.update(rmsd_data)
            
            # Per-residue RMSD
            per_res_rmsd = self.calculate_per_residue_rmsd(model_file)
            results['mean_per_res_rmsd'] = np.mean(per_res_rmsd)
            results['max_per_res_rmsd'] = np.max(per_res_rmsd)
            
            # Store chromophore-specific per-residue RMSD
            chrom_start, chrom_end = self.chromophore_range
            if len(per_res_rmsd) >= chrom_end:
                chrom_rmsd = per_res_rmsd[chrom_start-1:chrom_end]
                results['chromophore_mean_rmsd'] = np.mean(chrom_rmsd)
                results['chromophore_max_rmsd'] = np.max(chrom_rmsd)
        
        except Exception as e:
            print(f"    Error calculating RMSD: {str(e)[:100]}")
        
        # Load confidence scores
        try:
            confidence = self.load_confidence_scores(design_dir)
            
            # Check for pLDDT (try multiple keys)
            plddt_key = None
            for key in ['plddt', 'atom_plddts', 'per_residue_plddt']:
                if key in confidence:
                    plddt_key = key
                    break
            
            if plddt_key:
                plddt = np.array(confidence[plddt_key])
                results['mean_plddt'] = np.mean(plddt)
                results['min_plddt'] = np.min(plddt)
                results['max_plddt'] = np.max(plddt)
                
                # Calculate fraction disordered (pLDDT < 50)
                results['fraction_disordered'] = np.sum(plddt < 50) / len(plddt)
                
                # Chromophore pLDDT
                chrom_start, chrom_end = self.chromophore_range
                if len(plddt) >= chrom_end:
                    chrom_plddt = plddt[chrom_start-1:chrom_end]
                    results['chromophore_mean_plddt'] = np.mean(chrom_plddt)
                    results['chromophore_min_plddt'] = np.min(chrom_plddt)
            
            # Other confidence metrics
            for key in ['ptm', 'iptm', 'ranking_score', 'mean_pae', 'max_pae']:
                if key in confidence:
                    results[key] = confidence[key]

            # PAE interface analysis: extract protein-DNA interface quality
            pae_key = None
            for key in ['pae', 'predicted_aligned_error']:
                if key in confidence:
                    pae_key = key
                    break
            if pae_key is not None:
                try:
                    pae_matrix = np.array(confidence[pae_key])
                    if pae_matrix.ndim == 2:
                        results['mean_pae'] = float(np.mean(pae_matrix))
                        # For bound state with DNA: interface PAE is the
                        # protein-rows x DNA-cols block. Protein is chain A
                        # (first N residues), DNA follows.
                        n_protein = len(self.template_ca_atoms) if self.template_ca_atoms else pae_matrix.shape[0]
                        if state == 'bound' and pae_matrix.shape[0] > n_protein:
                            # Interface block: protein rows, DNA columns
                            interface_pae = pae_matrix[:n_protein, n_protein:]
                            results['interface_pae'] = float(np.mean(interface_pae))
                            results['interface_pae_median'] = float(np.median(interface_pae))
                except Exception:
                    pass

        except Exception as e:
            print(f"    Error loading confidence: {str(e)[:100]}")

        # Multi-seed consensus: check for multiple model files and compute
        # pLDDT variance across seeds (lower variance = more reliable)
        try:
            all_model_files = sorted(design_dir.glob("*_model*.cif"))
            if len(all_model_files) > 1:
                seed_plddts = []
                for mf in all_model_files:
                    # Find corresponding confidence JSON
                    stem = mf.stem  # e.g., design_0000_bound_model_0
                    conf_file = mf.parent / f"{stem}_summary_confidences.json"
                    if not conf_file.exists():
                        conf_file = mf.parent / f"summary_confidences_{stem.split('_')[-1]}.json"
                    if not conf_file.exists():
                        continue
                    with open(conf_file) as f:
                        seed_conf = json.load(f)
                    for key in ['plddt', 'atom_plddts', 'per_residue_plddt']:
                        if key in seed_conf:
                            seed_plddts.append(np.mean(seed_conf[key]))
                            break

                if len(seed_plddts) >= 2:
                    results['n_seeds'] = len(seed_plddts)
                    results['seed_plddt_mean'] = float(np.mean(seed_plddts))
                    results['seed_plddt_std'] = float(np.std(seed_plddts))
                    # Consensus score: penalize high inter-seed variance
                    # Range 0-1: 1 = perfect consensus, 0 = maximum disagreement
                    max_possible_std = 50  # pLDDT range 0-100, std of uniform ~29
                    results['consensus_score'] = max(0, 1 - results['seed_plddt_std'] / max_possible_std)
        except Exception:
            pass

        return results
    
    def analyze_all_designs(self, state: str = 'both') -> pd.DataFrame:
        """
        Analyze all AlphaFold3 predictions in output directory
        
        Args:
            state: 'bound', 'apo', or 'both'
        """
        states_to_analyze = []
        if state in ['bound', 'both']:
            states_to_analyze.append('bound')
        if state in ['apo', 'both']:
            states_to_analyze.append('apo')
        
        all_results = []
        
        for state_name in states_to_analyze:
            af3_dir = self.output_dir / f"03_alphafold3_predictions_{state_name}"
            
            if not af3_dir.exists():
                print(f"⚠ Predictions directory not found: {af3_dir}")
                continue
            
            # Find all design directories
            design_dirs = [d for d in af3_dir.iterdir() 
                          if d.is_dir() and d.name.startswith('design_')]
            
            print(f"\nFound {len(design_dirs)} {state_name} designs to analyze")
            
            for i, design_dir in enumerate(sorted(design_dirs), 1):
                print(f"[{i}/{len(design_dirs)}]", end=" ")
                result = self.analyze_design(design_dir, state_name)
                if result:
                    all_results.append(result)
        
        df = pd.DataFrame(all_results)
        print(f"\nTotal designs analyzed: {len(df)}")
        return df
    
    def load_ligandmpnn_scores(self) -> pd.DataFrame:
        """Load LigandMPNN confidence scores from FASTA file"""
        fasta_file = self.output_dir / "01_ligandmpnn_outputs" / "parsed_sequences.fasta"
        
        if not fasta_file.exists():
            print(f"⚠ LigandMPNN FASTA not found: {fasta_file}")
            return pd.DataFrame()
        
        records = []
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Parse header: >seq_id, score=X, ligand_confidence=Y, ...
                    header = line.strip()[1:]
                    parts = [p.strip() for p in header.split(',')]
                    
                    record = {}
                    record['seq_id'] = int(parts[0])
                    
                    for part in parts[1:]:
                        if '=' in part:
                            key, value = part.split('=', 1)
                            key = key.strip()
                            value = value.strip()
                            try:
                                record[key] = float(value)
                            except ValueError:
                                record[key] = value
                    
                    records.append(record)
        
        return pd.DataFrame(records)
    
    def create_visualizations(self, df: pd.DataFrame, ligandmpnn_df: pd.DataFrame,
                            output_dir: Path):
        """Create comprehensive visualizations"""
        output_dir.mkdir(exist_ok=True, parents=True)
        csv_dir = output_dir / "plot_data_csvs"
        csv_dir.mkdir(exist_ok=True)
        
        # Merge LigandMPNN scores with AF3 data
        if not ligandmpnn_df.empty:
            df = df.merge(ligandmpnn_df, on='seq_id', how='left')
        
        print("\nGenerating visualizations...")
        
        # 1. Overall confidence distribution (LigandMPNN scores)
        if 'score' in df.columns:
            self._plot_confidence_distribution(df, output_dir, csv_dir)
        
        # 2. AlphaFold3 metrics
        self._plot_af3_metrics(df, output_dir, csv_dir)
        
        # 3. Confidence vs pLDDT correlation
        if 'score' in df.columns and 'mean_plddt' in df.columns:
            self._plot_confidence_plddt_correlation(df, output_dir, csv_dir)
        
        # 4. RMSD distributions
        self._plot_rmsd_distributions(df, output_dir, csv_dir)
        
        # 5. Correlation matrix
        self._plot_correlation_matrix(df, output_dir, csv_dir)
        
        # 6. Global vs chromophore RMSD
        self._plot_global_vs_chromophore_rmsd(df, output_dir, csv_dir)
        
        # 7. Top designs summary
        self._export_top_designs(df, output_dir, csv_dir)
        
        # 8. Per-residue analysis for top design
        self._plot_per_residue_analysis(df, output_dir, csv_dir)
        
        print(f"\nAll plots saved to: {output_dir}")
        print(f"CSV data exported to: {csv_dir}")
    
    def _plot_confidence_distribution(self, df: pd.DataFrame, output_dir: Path, csv_dir: Path):
        """Plot LigandMPNN confidence distributions"""
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Overall confidence
        if 'score' in df.columns:
            for state in df['state'].unique():
                state_df = df[df['state'] == state]
                axes[0].hist(state_df['score'], bins=30, alpha=0.6, label=state)
            axes[0].set_xlabel('LigandMPNN Score')
            axes[0].set_ylabel('Frequency')
            axes[0].set_title('Overall Confidence Distribution')
            axes[0].legend()
            
            # Export data
            df[['seq_id', 'state', 'score']].to_csv(
                csv_dir / '01a_overall_confidence_data.csv', index=False
            )
        
        # Ligand confidence (if available)
        if 'ligand_confidence' in df.columns:
            for state in df['state'].unique():
                state_df = df[df['state'] == state]
                axes[1].hist(state_df['ligand_confidence'], bins=30, alpha=0.6, label=state)
            axes[1].set_xlabel('Ligand Confidence')
            axes[1].set_ylabel('Frequency')
            axes[1].set_title('Ligand-Specific Confidence Distribution')
            axes[1].legend()
            
            # Export data
            df[['seq_id', 'state', 'ligand_confidence']].to_csv(
                csv_dir / '01b_ligand_confidence_data.csv', index=False
            )
        
        plt.tight_layout()
        plt.savefig(output_dir / '01_confidence_distributions.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("  ✓ 01_confidence_distributions.png")
    
    def _plot_af3_metrics(self, df: pd.DataFrame, output_dir: Path, csv_dir: Path):
        """Plot AlphaFold3 quality metrics"""
        metrics = ['mean_plddt', 'ptm', 'iptm', 'ranking_score']
        available_metrics = [m for m in metrics if m in df.columns]
        
        if not available_metrics:
            print("  ⚠ No AF3 metrics available")
            return
        
        n_metrics = len(available_metrics)
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        axes = axes.flatten()
        
        for idx, metric in enumerate(available_metrics):
            for state in df['state'].unique():
                state_df = df[df['state'] == state]
                axes[idx].hist(state_df[metric].dropna(), bins=30, alpha=0.6, label=state)
            axes[idx].set_xlabel(metric.replace('_', ' ').title())
            axes[idx].set_ylabel('Frequency')
            axes[idx].set_title(f'{metric.replace("_", " ").title()} Distribution')
            axes[idx].legend()
        
        # Hide unused subplots
        for idx in range(n_metrics, 4):
            axes[idx].axis('off')
        
        plt.tight_layout()
        plt.savefig(output_dir / '02_alphafold3_metrics.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export data
        export_cols = ['seq_id', 'state'] + available_metrics
        df[export_cols].to_csv(csv_dir / '02_alphafold3_metrics_data.csv', index=False)
        print("  ✓ 02_alphafold3_metrics.png")
    
    def _plot_confidence_plddt_correlation(self, df: pd.DataFrame, output_dir: Path, csv_dir: Path):
        """Plot confidence vs pLDDT correlation"""
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Overall correlation
        for state in df['state'].unique():
            state_df = df[df['state'] == state].dropna(subset=['score', 'mean_plddt'])
            if len(state_df) > 0:
                axes[0].scatter(state_df['score'], state_df['mean_plddt'], 
                              alpha=0.5, label=state, s=20)
                
                # Add correlation coefficient
                corr = state_df[['score', 'mean_plddt']].corr().iloc[0, 1]
                axes[0].text(0.05, 0.95 - 0.05*list(df['state'].unique()).index(state), 
                           f'{state}: r={corr:.3f}',
                           transform=axes[0].transAxes, fontsize=10)
        
        axes[0].set_xlabel('LigandMPNN Score')
        axes[0].set_ylabel('Mean pLDDT')
        axes[0].set_title('Confidence vs pLDDT Correlation')
        axes[0].legend()
        
        # Chromophore-specific
        if 'chromophore_mean_plddt' in df.columns:
            for state in df['state'].unique():
                state_df = df[df['state'] == state].dropna(subset=['score', 'chromophore_mean_plddt'])
                if len(state_df) > 0:
                    axes[1].scatter(state_df['score'], state_df['chromophore_mean_plddt'],
                                  alpha=0.5, label=state, s=20)
                    
                    corr = state_df[['score', 'chromophore_mean_plddt']].corr().iloc[0, 1]
                    axes[1].text(0.05, 0.95 - 0.05*list(df['state'].unique()).index(state),
                               f'{state}: r={corr:.3f}',
                               transform=axes[1].transAxes, fontsize=10)
            
            axes[1].set_xlabel('LigandMPNN Score')
            axes[1].set_ylabel('Chromophore Mean pLDDT')
            axes[1].set_title('Confidence vs Chromophore pLDDT')
            axes[1].legend()
        
        plt.tight_layout()
        plt.savefig(output_dir / '03_confidence_vs_plddt.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export data
        df[['seq_id', 'state', 'score', 'mean_plddt']].to_csv(
            csv_dir / '03_confidence_vs_plddt_data.csv', index=False
        )
        if 'chromophore_mean_plddt' in df.columns:
            df[['seq_id', 'state', 'score', 'chromophore_mean_plddt']].to_csv(
                csv_dir / '03b_chromophore_confidence_vs_plddt_data.csv', index=False
            )
        print("  ✓ 03_confidence_vs_plddt.png")
    
    def _plot_rmsd_distributions(self, df: pd.DataFrame, output_dir: Path, csv_dir: Path):
        """Plot RMSD distributions"""
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Global RMSD
        if 'global_rmsd' in df.columns:
            for state in df['state'].unique():
                state_df = df[df['state'] == state]
                axes[0].hist(state_df['global_rmsd'].dropna(), bins=30, alpha=0.6, label=state)
            axes[0].set_xlabel('Global RMSD (Å)')
            axes[0].set_ylabel('Frequency')
            axes[0].set_title('Global RMSD Distribution')
            axes[0].legend()
        
        # Chromophore RMSD
        if 'chromophore_rmsd' in df.columns:
            for state in df['state'].unique():
                state_df = df[df['state'] == state]
                axes[1].hist(state_df['chromophore_rmsd'].dropna(), bins=30, alpha=0.6, label=state)
            axes[1].set_xlabel('Chromophore RMSD (Å)')
            axes[1].set_ylabel('Frequency')
            axes[1].set_title('Chromophore Region RMSD Distribution')
            axes[1].legend()
        
        plt.tight_layout()
        plt.savefig(output_dir / '04_rmsd_distributions.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export data
        df[['seq_id', 'state', 'global_rmsd', 'chromophore_rmsd']].to_csv(
            csv_dir / '04_rmsd_distributions_data.csv', index=False
        )
        print("  ✓ 04_rmsd_distributions.png")
    
    def _plot_correlation_matrix(self, df: pd.DataFrame, output_dir: Path, csv_dir: Path):
        """Plot correlation matrix of all metrics"""
        # Select numeric columns
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        exclude_cols = ['seq_id']
        numeric_cols = [col for col in numeric_cols if col not in exclude_cols]
        
        if len(numeric_cols) < 2:
            print("  ⚠ Not enough numeric columns for correlation matrix")
            return
        
        corr_matrix = df[numeric_cols].corr()
        
        plt.figure(figsize=(12, 10))
        sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm', 
                   center=0, square=True, linewidths=1)
        plt.title('Correlation Matrix of All Metrics')
        plt.tight_layout()
        plt.savefig(output_dir / '05_correlation_matrix.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export data
        corr_matrix.to_csv(csv_dir / '05_correlation_matrix_data.csv')
        print("  ✓ 05_correlation_matrix.png")
    
    def _plot_global_vs_chromophore_rmsd(self, df: pd.DataFrame, output_dir: Path, csv_dir: Path):
        """Plot global vs chromophore RMSD"""
        if 'global_rmsd' not in df.columns or 'chromophore_rmsd' not in df.columns:
            print("  ⚠ RMSD data not available")
            return
        
        plt.figure(figsize=(10, 8))
        
        for state in df['state'].unique():
            state_df = df[df['state'] == state].dropna(subset=['global_rmsd', 'chromophore_rmsd'])
            if len(state_df) > 0:
                plt.scatter(state_df['global_rmsd'], state_df['chromophore_rmsd'],
                          alpha=0.5, label=state, s=30)
                
                # Add correlation
                corr = state_df[['global_rmsd', 'chromophore_rmsd']].corr().iloc[0, 1]
                plt.text(0.05, 0.95 - 0.05*list(df['state'].unique()).index(state),
                        f'{state}: r={corr:.3f}',
                        transform=plt.gca().transAxes, fontsize=10)
        
        plt.xlabel('Global RMSD (Å)')
        plt.ylabel('Chromophore RMSD (Å)')
        plt.title('Global vs Chromophore RMSD')
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_dir / '06_global_vs_chromophore_rmsd.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export data
        df[['seq_id', 'state', 'global_rmsd', 'chromophore_rmsd']].to_csv(
            csv_dir / '06_global_vs_chromophore_rmsd_data.csv', index=False
        )
        print("  ✓ 06_global_vs_chromophore_rmsd.png")
    
    def _export_top_designs(self, df: pd.DataFrame, output_dir: Path, csv_dir: Path):
        """Export top 20 designs by different metrics"""
        metrics_to_rank = []
        
        if 'score' in df.columns:
            metrics_to_rank.append(('score', False))  # Higher is better
        if 'mean_plddt' in df.columns:
            metrics_to_rank.append(('mean_plddt', False))
        if 'global_rmsd' in df.columns:
            metrics_to_rank.append(('global_rmsd', True))  # Lower is better
        if 'chromophore_rmsd' in df.columns:
            metrics_to_rank.append(('chromophore_rmsd', True))
        
        top_designs = {}
        
        for metric, ascending in metrics_to_rank:
            for state in df['state'].unique():
                state_df = df[df['state'] == state]
                top_20 = state_df.nsmallest(20, metric) if ascending else state_df.nlargest(20, metric)
                top_designs[f'top_20_{state}_{metric}'] = top_20
        
        # Save to single CSV
        with pd.ExcelWriter(csv_dir / '07_top20_designs_data.xlsx', engine='openpyxl') as writer:
            for sheet_name, data in top_designs.items():
                data.to_excel(writer, sheet_name=sheet_name[:31], index=False)  # Excel limit
        
        # Also save as CSV (combined)
        all_top = pd.concat([df.assign(metric=key) for key, df in top_designs.items()], ignore_index=True)
        all_top.to_csv(csv_dir / '07_top20_designs_data.csv', index=False)
        print("  ✓ 07_top20_designs_data exported")
    
    def _plot_per_residue_analysis(self, df: pd.DataFrame, output_dir: Path, csv_dir: Path):
        """Plot per-residue RMSD for top design"""
        # Find top design by mean_plddt
        if 'mean_plddt' not in df.columns:
            print("  ⚠ Cannot plot per-residue analysis without pLDDT data")
            return
        
        top_design = df.loc[df['mean_plddt'].idxmax()]
        design_name = top_design['design_name']
        state = top_design['state']
        
        # Find design directory
        design_dir = self.output_dir / f"03_alphafold3_predictions_{state}" / design_name
        model_file = list(design_dir.glob("*_model.cif"))[0]
        
        # Calculate per-residue RMSD
        per_res_rmsd = self.calculate_per_residue_rmsd(model_file)
        
        # Plot
        plt.figure(figsize=(14, 6))
        plt.plot(range(1, len(per_res_rmsd)+1), per_res_rmsd, linewidth=1)
        plt.xlabel('Residue Number')
        plt.ylabel('RMSD (Å)')
        plt.title(f'Per-Residue RMSD for Top Design: {design_name}')
        
        # Highlight chromophore region
        chrom_start, chrom_end = self.chromophore_range
        plt.axvspan(chrom_start, chrom_end, alpha=0.2, color='red', label='Chromophore')
        plt.legend()
        
        plt.tight_layout()
        plt.savefig(output_dir / f'08_per_residue_rmsd_{design_name}.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export data
        per_res_df = pd.DataFrame({
            'residue': range(1, len(per_res_rmsd)+1),
            'rmsd': per_res_rmsd
        })
        per_res_df.to_csv(csv_dir / f'08_per_residue_rmsd_{design_name}_data.csv', index=False)
        print(f"  ✓ 08_per_residue_rmsd_{design_name}.png")


def main():
    parser = argparse.ArgumentParser(description='Analyze AlphaFold3 structure predictions')
    parser.add_argument('--output-dir', type=str, 
                       help='Specific output directory to analyze (e.g., output_G4FP_des1_cro_mod0)')
    parser.add_argument('--template', type=str, required=True,
                       help='Path to template/reference structure (PDB or CIF)')
    parser.add_argument('--state', type=str, default='both', choices=['bound', 'apo', 'both'],
                       help='Which predictions to analyze: bound, apo, or both')
    parser.add_argument('--chromophore-range', type=str, default='175,225',
                       help='Chromophore residue range as "start,end" (default: 175,225)')
    
    args = parser.parse_args()
    
    # Parse chromophore range
    chrom_start, chrom_end = map(int, args.chromophore_range.split(','))
    chromophore_range = (chrom_start, chrom_end)
    
    # If no specific output directory, find all
    if args.output_dir:
        output_dirs = [Path(args.output_dir)]
    else:
        # Find all output directories in current directory
        output_dirs = [d for d in Path('.').iterdir() 
                      if d.is_dir() and d.name.startswith('output_')]
    
    if not output_dirs:
        print("No output directories found. Please specify --output-dir")
        sys.exit(1)
    
    print(f"Found {len(output_dirs)} output directories to analyze")
    
    # Analyze each output directory
    for output_dir in output_dirs:
        print(f"\n{'='*80}")
        print(f"Analyzing: {output_dir}")
        print(f"{'='*80}")
        
        try:
            analyzer = AlphaFoldStructureAnalyzer(
                output_dir=output_dir,
                template_structure=args.template,
                chromophore_range=chromophore_range
            )
            
            # Analyze all designs
            df = analyzer.analyze_all_designs(state=args.state)
            
            if df.empty:
                print(f"⚠ No designs analyzed for {output_dir}")
                continue
            
            # Load LigandMPNN scores
            ligandmpnn_df = analyzer.load_ligandmpnn_scores()
            
            # Create visualizations
            viz_output_dir = output_dir / "05_structure_analysis"
            analyzer.create_visualizations(df, ligandmpnn_df, viz_output_dir)
            
            # Save combined data
            df.to_csv(output_dir / "05_structure_analysis_results.csv", index=False)
            print(f"\n✓ Analysis complete for {output_dir}")
            print(f"  Results: {output_dir}/05_structure_analysis_results.csv")
            print(f"  Plots: {viz_output_dir}/")
            
        except Exception as e:
            print(f"✗ Error analyzing {output_dir}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue


if __name__ == "__main__":
    main()

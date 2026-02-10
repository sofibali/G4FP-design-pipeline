#!/usr/bin/env python3
"""
Compare AlphaFold3 Predictions: Ligand-Bound vs Apo States

Analyzes structural differences between chromophore-bound and apo states:
1. RMSD between bound and apo structures
2. pLDDT differences between states
3. Per-residue conformational changes
4. Chromophore region analysis
5. Comprehensive visualizations

Usage:
    python 06_compare_ligand_states.py --output-dir OUTPUT --template TEMPLATE

Arguments:
    --output-dir: Specific output directory to analyze (e.g., output_G4FP_des1_cro_mod0)
    --template: Path to template structure for reference
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

# BioPython for structure analysis
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import PPBuilder

warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10


class LigandStateComparator:
    """Compare AlphaFold3 predictions with and without chromophore"""
    
    def __init__(self, output_dir: Path, template_structure: str,
                 chromophore_range: Tuple[int, int] = (175, 225)):
        """
        Initialize comparator
        
        Args:
            output_dir: Path to output directory
            template_structure: Path to template structure
            chromophore_range: Tuple of (start, end) residue numbers
        """
        self.output_dir = Path(output_dir)
        self.template_path = Path(template_structure)
        self.chromophore_range = chromophore_range
        
        # Initialize parsers
        self.pdb_parser = PDBParser(QUIET=True)
        self.cif_parser = MMCIFParser(QUIET=True)
        
        print(f"Output directory: {self.output_dir}")
        print(f"Template structure: {self.template_path}")
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
    
    def calculate_structure_rmsd(self, struct1_path: Path, struct2_path: Path,
                                 chain_id='A') -> Dict[str, float]:
        """
        Calculate RMSD between two structures
        
        Returns:
            Dict with 'global_rmsd', 'chromophore_rmsd', and per-residue RMSD
        """
        try:
            # Load structures
            struct1 = self._load_structure(struct1_path)
            struct2 = self._load_structure(struct2_path)
            
            # Extract CA atoms
            ca1 = self._extract_ca_atoms(struct1, chain_id)
            ca2 = self._extract_ca_atoms(struct2, chain_id)
            
            # Match lengths
            min_len = min(len(ca1), len(ca2))
            ca1 = ca1[:min_len]
            ca2 = ca2[:min_len]
            
            # Global RMSD
            sup = Superimposer()
            sup.set_atoms(ca1, ca2)
            global_rmsd = sup.rms
            
            # Apply superposition
            sup.apply(struct2)
            
            # Per-residue RMSD
            per_res_rmsd = []
            for atom1, atom2 in zip(ca1, ca2):
                distance = np.linalg.norm(
                    np.array(atom1.get_coord()) - np.array(atom2.get_coord())
                )
                per_res_rmsd.append(distance)
            per_res_rmsd = np.array(per_res_rmsd)
            
            # Chromophore RMSD
            chrom_start, chrom_end = self.chromophore_range
            if len(ca1) >= chrom_end:
                chrom_ca1 = ca1[chrom_start-1:chrom_end]
                chrom_ca2 = ca2[chrom_start-1:chrom_end]
                
                sup_chrom = Superimposer()
                sup_chrom.set_atoms(chrom_ca1, chrom_ca2)
                chromophore_rmsd = sup_chrom.rms
            else:
                chromophore_rmsd = None
            
            return {
                'global_rmsd': global_rmsd,
                'chromophore_rmsd': chromophore_rmsd,
                'per_residue_rmsd': per_res_rmsd,
                'mean_per_res_rmsd': np.mean(per_res_rmsd),
                'max_per_res_rmsd': np.max(per_res_rmsd)
            }
            
        except Exception as e:
            print(f"    Error calculating RMSD: {str(e)[:100]}")
            return {}
    
    def load_confidence_scores(self, design_dir: Path) -> Dict:
        """Load confidence scores from AlphaFold3 summary JSON"""
        summary_file = design_dir / "summary_confidences_0.json"
        
        if not summary_file.exists():
            json_files = list(design_dir.glob("**/summary_confidences*.json"))
            if json_files:
                summary_file = json_files[0]
            else:
                raise FileNotFoundError(f"No confidence file found in {design_dir}")
        
        with open(summary_file, 'r') as f:
            confidence = json.load(f)
        
        return confidence
    
    def extract_plddt(self, confidence: Dict) -> np.ndarray:
        """Extract pLDDT array from confidence dict"""
        for key in ['plddt', 'atom_plddts', 'per_residue_plddt']:
            if key in confidence:
                return np.array(confidence[key])
        return None
    
    def compare_design_pair(self, seq_id: int) -> Dict:
        """
        Compare bound and apo predictions for a single design
        
        Args:
            seq_id: Design sequence ID (e.g., 0000)
            
        Returns:
            Dictionary with comparison metrics
        """
        print(f"  Comparing design {seq_id:04d}...")
        
        results = {
            'seq_id': seq_id,
        }
        
        # Find bound and apo directories
        bound_dir = self.output_dir / "03_alphafold3_predictions_bound" / f"design_{seq_id:04d}_bound"
        apo_dir = self.output_dir / "03_alphafold3_predictions_apo" / f"design_{seq_id:04d}_apo"
        
        if not bound_dir.exists():
            print(f"    ⚠ Bound directory not found: {bound_dir}")
            return results
        
        if not apo_dir.exists():
            print(f"    ⚠ Apo directory not found: {apo_dir}")
            return results
        
        # Find model files
        bound_model = list(bound_dir.glob("*_model.cif"))
        apo_model = list(apo_dir.glob("*_model.cif"))
        
        if not bound_model or not apo_model:
            print(f"    ⚠ Model files not found")
            return results
        
        bound_model = bound_model[0]
        apo_model = apo_model[0]
        
        # Calculate RMSD between bound and apo
        try:
            rmsd_data = self.calculate_structure_rmsd(bound_model, apo_model)
            results.update(rmsd_data)
            
            # Store per-residue RMSD for later plotting
            if 'per_residue_rmsd' in rmsd_data:
                results['per_residue_rmsd_array'] = rmsd_data['per_residue_rmsd']
        
        except Exception as e:
            print(f"    Error calculating RMSD: {str(e)[:100]}")
        
        # Load pLDDT scores for both states
        try:
            bound_conf = self.load_confidence_scores(bound_dir)
            apo_conf = self.load_confidence_scores(apo_dir)
            
            bound_plddt = self.extract_plddt(bound_conf)
            apo_plddt = self.extract_plddt(apo_conf)
            
            if bound_plddt is not None and apo_plddt is not None:
                # Match lengths
                min_len = min(len(bound_plddt), len(apo_plddt))
                bound_plddt = bound_plddt[:min_len]
                apo_plddt = apo_plddt[:min_len]
                
                # Calculate differences
                plddt_diff = bound_plddt - apo_plddt
                results['mean_plddt_diff'] = np.mean(plddt_diff)
                results['max_plddt_diff'] = np.max(np.abs(plddt_diff))
                
                results['bound_mean_plddt'] = np.mean(bound_plddt)
                results['apo_mean_plddt'] = np.mean(apo_plddt)
                
                # Chromophore region
                chrom_start, chrom_end = self.chromophore_range
                if len(bound_plddt) >= chrom_end:
                    chrom_bound = bound_plddt[chrom_start-1:chrom_end]
                    chrom_apo = apo_plddt[chrom_start-1:chrom_end]
                    
                    results['chromophore_bound_plddt'] = np.mean(chrom_bound)
                    results['chromophore_apo_plddt'] = np.mean(chrom_apo)
                    results['chromophore_plddt_diff'] = np.mean(chrom_bound - chrom_apo)
                
                # Store arrays for later plotting
                results['plddt_diff_array'] = plddt_diff
                results['bound_plddt_array'] = bound_plddt
                results['apo_plddt_array'] = apo_plddt
            
            # Other metrics
            for key in ['ptm', 'iptm', 'ranking_score']:
                if key in bound_conf:
                    results[f'bound_{key}'] = bound_conf[key]
                if key in apo_conf:
                    results[f'apo_{key}'] = apo_conf[key]

            # PAE interface analysis for bound state
            for pae_key in ['pae', 'predicted_aligned_error']:
                if pae_key in bound_conf:
                    try:
                        pae_matrix = np.array(bound_conf[pae_key])
                        if pae_matrix.ndim == 2:
                            results['bound_mean_pae'] = float(np.mean(pae_matrix))
                            # Interface PAE: protein-DNA block
                            # Estimate protein length from pLDDT if available
                            n_prot = len(bound_plddt) if bound_plddt is not None else pae_matrix.shape[0]
                            if pae_matrix.shape[0] > n_prot:
                                interface_block = pae_matrix[:n_prot, n_prot:]
                                results['bound_interface_pae'] = float(np.mean(interface_block))
                    except Exception:
                        pass
                    break
            for pae_key in ['pae', 'predicted_aligned_error']:
                if pae_key in apo_conf:
                    try:
                        pae_matrix = np.array(apo_conf[pae_key])
                        if pae_matrix.ndim == 2:
                            results['apo_mean_pae'] = float(np.mean(pae_matrix))
                    except Exception:
                        pass
                    break

            # PAE difference
            if 'bound_mean_pae' in results and 'apo_mean_pae' in results:
                results['pae_diff'] = results['apo_mean_pae'] - results['bound_mean_pae']

        except Exception as e:
            print(f"    Error loading confidence: {str(e)[:100]}")

        # Multi-seed consensus scoring
        for state_label, state_dir in [('bound', bound_dir), ('apo', apo_dir)]:
            try:
                model_files = sorted(state_dir.glob("*_model*.cif"))
                if len(model_files) > 1:
                    seed_plddts = []
                    for mf in model_files:
                        stem = mf.stem
                        conf_file = mf.parent / f"{stem}_summary_confidences.json"
                        if not conf_file.exists():
                            conf_file = mf.parent / f"summary_confidences_{stem.split('_')[-1]}.json"
                        if not conf_file.exists():
                            continue
                        with open(conf_file) as f:
                            sc = json.load(f)
                        for k in ['plddt', 'atom_plddts', 'per_residue_plddt']:
                            if k in sc:
                                seed_plddts.append(np.mean(sc[k]))
                                break
                    if len(seed_plddts) >= 2:
                        results[f'{state_label}_n_seeds'] = len(seed_plddts)
                        results[f'{state_label}_seed_plddt_std'] = float(np.std(seed_plddts))
                        results[f'{state_label}_consensus'] = max(0, 1 - np.std(seed_plddts) / 50)
            except Exception:
                pass

        return results
    
    def compare_all_designs(self) -> pd.DataFrame:
        """Compare all design pairs (bound vs apo)"""
        # Find all bound designs
        bound_dir = self.output_dir / "03_alphafold3_predictions_bound"
        
        if not bound_dir.exists():
            raise FileNotFoundError(f"Bound predictions not found: {bound_dir}")
        
        design_dirs = [d for d in bound_dir.iterdir() 
                      if d.is_dir() and d.name.startswith('design_')]
        
        print(f"\nFound {len(design_dirs)} design pairs to compare")
        
        results = []
        for i, design_dir in enumerate(sorted(design_dirs), 1):
            seq_id = int(design_dir.name.replace('design_', '').replace('_bound', ''))
            print(f"[{i}/{len(design_dirs)}]", end=" ")
            result = self.compare_design_pair(seq_id)
            if result:
                results.append(result)
        
        df = pd.DataFrame(results)
        print(f"\nTotal design pairs compared: {len(df)}")
        return df
    
    def create_visualizations(self, df: pd.DataFrame, output_dir: Path):
        """Create comprehensive visualizations"""
        output_dir.mkdir(exist_ok=True, parents=True)
        csv_dir = output_dir / "plot_data_csvs"
        csv_dir.mkdir(exist_ok=True)
        
        print("\nGenerating visualizations...")
        
        # 1. RMSD distributions (bound vs apo)
        self._plot_rmsd_distributions(df, output_dir, csv_dir)
        
        # 2. pLDDT comparison
        self._plot_plddt_comparison(df, output_dir, csv_dir)
        
        # 3. Chromophore pLDDT comparison
        self._plot_chromophore_plddt(df, output_dir, csv_dir)
        
        # 4. RMSD vs pLDDT change
        self._plot_rmsd_vs_plddt_change(df, output_dir, csv_dir)
        
        # 5. Top conformational changes
        self._plot_top_changes(df, output_dir, csv_dir)
        
        # 6. Correlation matrix
        self._plot_correlation_matrix(df, output_dir, csv_dir)
        
        print(f"\nAll plots saved to: {output_dir}")
        print(f"CSV data exported to: {csv_dir}")
    
    def _plot_rmsd_distributions(self, df: pd.DataFrame, output_dir: Path, csv_dir: Path):
        """Plot RMSD distributions between bound and apo"""
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Global RMSD
        if 'global_rmsd' in df.columns:
            axes[0].hist(df['global_rmsd'].dropna(), bins=30, alpha=0.7, color='steelblue')
            axes[0].set_xlabel('Global RMSD (Å)')
            axes[0].set_ylabel('Frequency')
            axes[0].set_title('Bound vs Apo: Global RMSD Distribution')
            axes[0].axvline(df['global_rmsd'].median(), color='red', linestyle='--', 
                          label=f'Median: {df["global_rmsd"].median():.2f}Å')
            axes[0].legend()
        
        # Chromophore RMSD
        if 'chromophore_rmsd' in df.columns:
            axes[1].hist(df['chromophore_rmsd'].dropna(), bins=30, alpha=0.7, color='coral')
            axes[1].set_xlabel('Chromophore RMSD (Å)')
            axes[1].set_ylabel('Frequency')
            axes[1].set_title('Bound vs Apo: Chromophore RMSD Distribution')
            axes[1].axvline(df['chromophore_rmsd'].median(), color='red', linestyle='--',
                          label=f'Median: {df["chromophore_rmsd"].median():.2f}Å')
            axes[1].legend()
        
        plt.tight_layout()
        plt.savefig(output_dir / '01_bound_apo_rmsd.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export data
        df[['seq_id', 'global_rmsd', 'chromophore_rmsd']].to_csv(
            csv_dir / '01_bound_apo_rmsd_data.csv', index=False
        )
        print("  ✓ 01_bound_apo_rmsd.png")
    
    def _plot_plddt_comparison(self, df: pd.DataFrame, output_dir: Path, csv_dir: Path):
        """Plot pLDDT comparison between states"""
        if 'bound_mean_plddt' not in df.columns or 'apo_mean_plddt' not in df.columns:
            print("  ⚠ pLDDT data not available")
            return
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Scatter plot
        axes[0].scatter(df['bound_mean_plddt'], df['apo_mean_plddt'], alpha=0.5, s=30)
        
        # Add diagonal line (y=x)
        min_val = min(df['bound_mean_plddt'].min(), df['apo_mean_plddt'].min())
        max_val = max(df['bound_mean_plddt'].max(), df['apo_mean_plddt'].max())
        axes[0].plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, label='y=x')
        
        axes[0].set_xlabel('Bound Mean pLDDT')
        axes[0].set_ylabel('Apo Mean pLDDT')
        axes[0].set_title('pLDDT: Bound vs Apo')
        axes[0].legend()
        
        # Add correlation
        corr = df[['bound_mean_plddt', 'apo_mean_plddt']].corr().iloc[0, 1]
        axes[0].text(0.05, 0.95, f'r = {corr:.3f}', transform=axes[0].transAxes, fontsize=10)
        
        # Distribution of differences
        if 'mean_plddt_diff' in df.columns:
            axes[1].hist(df['mean_plddt_diff'].dropna(), bins=30, alpha=0.7, color='green')
            axes[1].set_xlabel('pLDDT Difference (Bound - Apo)')
            axes[1].set_ylabel('Frequency')
            axes[1].set_title('pLDDT Change Distribution')
            axes[1].axvline(0, color='red', linestyle='--', label='No change')
            axes[1].axvline(df['mean_plddt_diff'].median(), color='blue', linestyle='--',
                          label=f'Median: {df["mean_plddt_diff"].median():.2f}')
            axes[1].legend()
        
        plt.tight_layout()
        plt.savefig(output_dir / '02_plddt_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export data
        df[['seq_id', 'bound_mean_plddt', 'apo_mean_plddt', 'mean_plddt_diff']].to_csv(
            csv_dir / '02_plddt_comparison_data.csv', index=False
        )
        print("  ✓ 02_plddt_comparison.png")
    
    def _plot_chromophore_plddt(self, df: pd.DataFrame, output_dir: Path, csv_dir: Path):
        """Plot chromophore region pLDDT comparison"""
        if 'chromophore_bound_plddt' not in df.columns or 'chromophore_apo_plddt' not in df.columns:
            print("  ⚠ Chromophore pLDDT data not available")
            return
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Scatter plot
        axes[0].scatter(df['chromophore_bound_plddt'], df['chromophore_apo_plddt'], 
                       alpha=0.5, s=30, color='purple')
        
        # Add diagonal
        min_val = min(df['chromophore_bound_plddt'].min(), df['chromophore_apo_plddt'].min())
        max_val = max(df['chromophore_bound_plddt'].max(), df['chromophore_apo_plddt'].max())
        axes[0].plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, label='y=x')
        
        axes[0].set_xlabel('Bound Chromophore pLDDT')
        axes[0].set_ylabel('Apo Chromophore pLDDT')
        axes[0].set_title('Chromophore pLDDT: Bound vs Apo')
        axes[0].legend()
        
        # Correlation
        corr = df[['chromophore_bound_plddt', 'chromophore_apo_plddt']].corr().iloc[0, 1]
        axes[0].text(0.05, 0.95, f'r = {corr:.3f}', transform=axes[0].transAxes, fontsize=10)
        
        # Distribution of differences
        if 'chromophore_plddt_diff' in df.columns:
            axes[1].hist(df['chromophore_plddt_diff'].dropna(), bins=30, alpha=0.7, color='orange')
            axes[1].set_xlabel('Chromophore pLDDT Difference (Bound - Apo)')
            axes[1].set_ylabel('Frequency')
            axes[1].set_title('Chromophore pLDDT Change Distribution')
            axes[1].axvline(0, color='red', linestyle='--', label='No change')
            axes[1].axvline(df['chromophore_plddt_diff'].median(), color='blue', linestyle='--',
                          label=f'Median: {df["chromophore_plddt_diff"].median():.2f}')
            axes[1].legend()
        
        plt.tight_layout()
        plt.savefig(output_dir / '03_chromophore_plddt_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export data
        df[['seq_id', 'chromophore_bound_plddt', 'chromophore_apo_plddt', 'chromophore_plddt_diff']].to_csv(
            csv_dir / '03_chromophore_plddt_comparison_data.csv', index=False
        )
        print("  ✓ 03_chromophore_plddt_comparison.png")
    
    def _plot_rmsd_vs_plddt_change(self, df: pd.DataFrame, output_dir: Path, csv_dir: Path):
        """Plot RMSD vs pLDDT change"""
        if 'global_rmsd' not in df.columns or 'mean_plddt_diff' not in df.columns:
            print("  ⚠ Required data not available for RMSD vs pLDDT plot")
            return
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Global RMSD vs pLDDT change
        data = df.dropna(subset=['global_rmsd', 'mean_plddt_diff'])
        axes[0].scatter(data['global_rmsd'], data['mean_plddt_diff'], alpha=0.5, s=30)
        axes[0].set_xlabel('Global RMSD (Å)')
        axes[0].set_ylabel('pLDDT Difference (Bound - Apo)')
        axes[0].set_title('RMSD vs pLDDT Change')
        axes[0].axhline(0, color='red', linestyle='--', alpha=0.5)
        
        # Add correlation
        corr = data[['global_rmsd', 'mean_plddt_diff']].corr().iloc[0, 1]
        axes[0].text(0.05, 0.95, f'r = {corr:.3f}', transform=axes[0].transAxes, fontsize=10)
        
        # Chromophore RMSD vs chromophore pLDDT change
        if 'chromophore_rmsd' in df.columns and 'chromophore_plddt_diff' in df.columns:
            data_chrom = df.dropna(subset=['chromophore_rmsd', 'chromophore_plddt_diff'])
            axes[1].scatter(data_chrom['chromophore_rmsd'], data_chrom['chromophore_plddt_diff'],
                          alpha=0.5, s=30, color='coral')
            axes[1].set_xlabel('Chromophore RMSD (Å)')
            axes[1].set_ylabel('Chromophore pLDDT Difference (Bound - Apo)')
            axes[1].set_title('Chromophore: RMSD vs pLDDT Change')
            axes[1].axhline(0, color='red', linestyle='--', alpha=0.5)
            
            corr_chrom = data_chrom[['chromophore_rmsd', 'chromophore_plddt_diff']].corr().iloc[0, 1]
            axes[1].text(0.05, 0.95, f'r = {corr_chrom:.3f}', transform=axes[1].transAxes, fontsize=10)
        
        plt.tight_layout()
        plt.savefig(output_dir / '04_rmsd_vs_plddt_change.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export data
        df[['seq_id', 'global_rmsd', 'mean_plddt_diff', 'chromophore_rmsd', 'chromophore_plddt_diff']].to_csv(
            csv_dir / '04_rmsd_vs_plddt_change_data.csv', index=False
        )
        print("  ✓ 04_rmsd_vs_plddt_change.png")
    
    def _plot_top_changes(self, df: pd.DataFrame, output_dir: Path, csv_dir: Path):
        """Plot top conformational changes"""
        if 'global_rmsd' not in df.columns:
            print("  ⚠ RMSD data not available for top changes plot")
            return
        
        # Get top 20 by RMSD
        top_20 = df.nlargest(20, 'global_rmsd')
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Bar plot
        x = range(len(top_20))
        ax.bar(x, top_20['global_rmsd'], alpha=0.7, color='steelblue')
        ax.set_xlabel('Design Rank')
        ax.set_ylabel('Global RMSD (Å)')
        ax.set_title('Top 20 Designs by Conformational Change (Bound vs Apo)')
        ax.set_xticks(x)
        ax.set_xticklabels([f"{sid:04d}" for sid in top_20['seq_id']], rotation=45)
        
        plt.tight_layout()
        plt.savefig(output_dir / '05_top20_conformational_changes.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export data
        top_20.to_csv(csv_dir / '05_top20_conformational_changes_data.csv', index=False)
        print("  ✓ 05_top20_conformational_changes.png")
    
    def _plot_correlation_matrix(self, df: pd.DataFrame, output_dir: Path, csv_dir: Path):
        """Plot correlation matrix of comparison metrics"""
        # Select relevant numeric columns
        cols_to_include = []
        possible_cols = ['global_rmsd', 'chromophore_rmsd', 'mean_per_res_rmsd',
                        'bound_mean_plddt', 'apo_mean_plddt', 'mean_plddt_diff',
                        'chromophore_bound_plddt', 'chromophore_apo_plddt', 'chromophore_plddt_diff',
                        'bound_ptm', 'apo_ptm', 'bound_iptm', 'apo_iptm']
        
        for col in possible_cols:
            if col in df.columns:
                cols_to_include.append(col)
        
        if len(cols_to_include) < 2:
            print("  ⚠ Not enough metrics for correlation matrix")
            return
        
        corr_matrix = df[cols_to_include].corr()
        
        plt.figure(figsize=(12, 10))
        sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm',
                   center=0, square=True, linewidths=1, cbar_kws={"shrink": 0.8})
        plt.title('Correlation Matrix: Bound vs Apo Comparison')
        plt.tight_layout()
        plt.savefig(output_dir / '06_correlation_matrix.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Export data
        corr_matrix.to_csv(csv_dir / '06_correlation_matrix_data.csv')
        print("  ✓ 06_correlation_matrix.png")


def main():
    parser = argparse.ArgumentParser(description='Compare bound and apo AlphaFold3 predictions')
    parser.add_argument('--output-dir', type=str, required=True,
                       help='Output directory to analyze (e.g., output_G4FP_des1_cro_mod0)')
    parser.add_argument('--template', type=str, required=True,
                       help='Path to template structure (PDB or CIF)')
    parser.add_argument('--chromophore-range', type=str, default='175,225',
                       help='Chromophore residue range as "start,end" (default: 175,225)')
    
    args = parser.parse_args()
    
    # Parse chromophore range
    chrom_start, chrom_end = map(int, args.chromophore_range.split(','))
    chromophore_range = (chrom_start, chrom_end)
    
    output_dir = Path(args.output_dir)
    
    if not output_dir.exists():
        print(f"✗ Output directory not found: {output_dir}")
        sys.exit(1)
    
    print(f"\n{'='*80}")
    print(f"Comparing Ligand States: {output_dir}")
    print(f"{'='*80}")
    
    try:
        comparator = LigandStateComparator(
            output_dir=output_dir,
            template_structure=args.template,
            chromophore_range=chromophore_range
        )
        
        # Compare all designs
        df = comparator.compare_all_designs()
        
        if df.empty:
            print("✗ No designs compared")
            sys.exit(1)
        
        # Create visualizations
        viz_output_dir = output_dir / "06_ligand_state_comparison"
        comparator.create_visualizations(df, viz_output_dir)
        
        # Save results
        # Remove array columns before saving CSV
        df_export = df.copy()
        array_cols = [col for col in df_export.columns if 'array' in col]
        df_export = df_export.drop(columns=array_cols)
        
        df_export.to_csv(output_dir / "06_ligand_state_comparison_results.csv", index=False)
        
        print(f"\n✓ Analysis complete!")
        print(f"  Results: {output_dir}/06_ligand_state_comparison_results.csv")
        print(f"  Plots: {viz_output_dir}/")
        
    except Exception as e:
        print(f"✗ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

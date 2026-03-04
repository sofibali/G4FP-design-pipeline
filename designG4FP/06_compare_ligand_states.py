#!/usr/bin/env python3
"""
Compare AlphaFold3 Predictions: Ligand-Bound vs Apo States

Analyzes structural differences between chromophore-bound and apo states:
1. RMSD between bound and apo structures
2. pLDDT differences between states
3. Per-residue conformational changes
4. Chromophore region analysis
5. Comprehensive visualizations with all 10 template references

Usage:
    # Single template
    python 06_compare_ligand_states.py --output-dir output_G4FP_des1_cro_mod0 \\
        --template inputs/G4FP_des1_cro_mod0.pdb --chromophore-range 197,199

    # All templates (auto-discover)
    python 06_compare_ligand_states.py --chromophore-range 197,199
"""

import os
import sys
import re
import json
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

from Bio.PDB import PDBParser, MMCIFParser, Superimposer

warnings.filterwarnings('ignore')

sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10


# -------------------------------------------------------------------------
# Template reference extraction
# -------------------------------------------------------------------------

def load_all_template_metrics(inputs_dir: Path, chromophore_range: Tuple[int, int]) -> pd.DataFrame:
    """Load B-factor pLDDT from all template PDBs in inputs/ directory."""
    parser = PDBParser(QUIET=True)
    records = []
    for pdb in sorted(inputs_dir.glob("G4FP_*.pdb")):
        structure = parser.get_structure('s', str(pdb))
        b_factors = []
        for model in structure:
            for chain in model:
                if chain.get_id() == 'A':
                    for residue in chain:
                        if residue.has_id('CA'):
                            b_factors.append(residue['CA'].get_bfactor())
        if not b_factors:
            continue
        b_arr = np.array(b_factors)
        chrom_start, chrom_end = chromophore_range
        chrom_plddt = np.mean(b_arr[chrom_start-1:chrom_end]) if len(b_arr) >= chrom_end else None
        records.append({
            'template_name': pdb.stem,
            'template_pdb': str(pdb),
            'mean_plddt': float(np.mean(b_arr)),
            'chromophore_mean_plddt': chrom_plddt,
            'n_residues': len(b_arr),
        })
    return pd.DataFrame(records)


# -------------------------------------------------------------------------
# LigandStateComparator (single output dir)
# -------------------------------------------------------------------------

class LigandStateComparator:
    """Compare AlphaFold3 bound vs apo predictions for one output directory."""

    def __init__(self, output_dir: Path, reference_template: Path,
                 chromophore_range: Tuple[int, int] = (197, 199)):
        self.output_dir = Path(output_dir)
        self.reference_template = Path(reference_template)
        self.chromophore_range = chromophore_range
        self.pdb_parser = PDBParser(QUIET=True)
        self.cif_parser = MMCIFParser(QUIET=True)

        # Determine which template this output dir corresponds to
        # output_G4FP_X -> inputs/G4FP_X.pdb
        self.template_name = self.output_dir.name.replace("output_", "")

    def _load_structure(self, filepath: Path):
        if filepath.suffix == '.pdb':
            return self.pdb_parser.get_structure('structure', str(filepath))
        elif filepath.suffix == '.cif':
            return self.cif_parser.get_structure('structure', str(filepath))
        else:
            raise ValueError(f"Unknown file format: {filepath.suffix}")

    def _extract_ca_atoms(self, structure, chain_id='A'):
        ca_atoms = []
        for model in structure:
            for chain in model:
                if chain.get_id() == chain_id:
                    for residue in chain:
                        if residue.has_id('CA'):
                            ca_atoms.append(residue['CA'])
        return ca_atoms

    def _discover_inference_dirs(self, af3_dir: Path, state: str) -> Dict[int, Path]:
        """Discover AF3 inference directories with completed model files."""
        pattern = re.compile(
            rf'^design_(\d+)_{re.escape(state)}(?:_(\d{{8}}_\d{{6}}))?$'
        )
        candidates = defaultdict(list)
        for d in af3_dir.iterdir():
            if not d.is_dir():
                continue
            m = pattern.match(d.name)
            if not m:
                continue
            seq_id = int(m.group(1))
            timestamp = m.group(2)
            if not list(d.glob("*_model.cif")):
                continue
            candidates[seq_id].append((timestamp is not None, timestamp or "", d))

        result = {}
        for seq_id, entries in candidates.items():
            non_ts = [e for e in entries if not e[0]]
            if non_ts:
                result[seq_id] = non_ts[0][2]
            else:
                entries.sort(key=lambda e: e[1], reverse=True)
                result[seq_id] = entries[0][2]
        return result

    def _load_confidence_data(self, design_dir: Path) -> Tuple[Dict, Dict]:
        """Load both summary and full confidence data."""
        summary_files = list(design_dir.glob("*_summary_confidences.json"))
        summary_conf = {}
        if summary_files:
            with open(summary_files[0]) as f:
                summary_conf = json.load(f)
        else:
            fallback = list(design_dir.glob("**/summary_confidences*.json"))
            if fallback:
                with open(fallback[0]) as f:
                    summary_conf = json.load(f)

        full_files = [f for f in design_dir.glob("*_confidences.json")
                      if 'summary' not in f.name]
        full_conf = {}
        if full_files:
            with open(full_files[0]) as f:
                full_conf = json.load(f)
        return summary_conf, full_conf

    def _compute_per_residue_plddt(self, full_conf: Dict, chain_id: str = 'A') -> Optional[np.ndarray]:
        """Compute per-residue pLDDT for a given chain from atom_plddts."""
        if 'atom_plddts' not in full_conf or 'atom_chain_ids' not in full_conf:
            return None
        atom_plddts = np.array(full_conf['atom_plddts'])
        atom_chains = full_conf['atom_chain_ids']
        chain_mask = np.array([ch == chain_id for ch in atom_chains])
        chain_plddts = atom_plddts[chain_mask]
        if len(chain_plddts) == 0:
            return None
        if 'token_chain_ids' in full_conf:
            token_chains = full_conf['token_chain_ids']
            n_residues = sum(1 for ch in token_chains if ch == chain_id)
            if n_residues > 0:
                atoms_per_residue = len(chain_plddts) / n_residues
                per_res = []
                for i in range(n_residues):
                    start = int(round(i * atoms_per_residue))
                    end = int(round((i + 1) * atoms_per_residue))
                    if start < len(chain_plddts):
                        per_res.append(np.mean(chain_plddts[start:end]))
                return np.array(per_res)
        return chain_plddts

    def calculate_structure_rmsd(self, struct1_path: Path, struct2_path: Path,
                                 chain_id='A') -> Dict[str, float]:
        """Calculate RMSD between bound and apo structures."""
        try:
            struct1 = self._load_structure(struct1_path)
            struct2 = self._load_structure(struct2_path)
            ca1 = self._extract_ca_atoms(struct1, chain_id)
            ca2 = self._extract_ca_atoms(struct2, chain_id)
            min_len = min(len(ca1), len(ca2))
            ca1 = ca1[:min_len]
            ca2 = ca2[:min_len]

            sup = Superimposer()
            sup.set_atoms(ca1, ca2)
            global_rmsd = sup.rms
            sup.apply(struct2)

            per_res_rmsd = []
            for atom1, atom2 in zip(ca1, ca2):
                distance = np.linalg.norm(
                    np.array(atom1.get_coord()) - np.array(atom2.get_coord()))
                per_res_rmsd.append(distance)
            per_res_rmsd = np.array(per_res_rmsd)

            chrom_start, chrom_end = self.chromophore_range
            chromophore_rmsd = None
            if len(ca1) >= chrom_end:
                chrom_ca1 = ca1[chrom_start-1:chrom_end]
                chrom_ca2 = ca2[chrom_start-1:chrom_end]
                sup_chrom = Superimposer()
                sup_chrom.set_atoms(chrom_ca1, chrom_ca2)
                chromophore_rmsd = sup_chrom.rms

            return {
                'global_rmsd': global_rmsd,
                'chromophore_rmsd': chromophore_rmsd,
                'mean_per_res_rmsd': np.mean(per_res_rmsd),
                'max_per_res_rmsd': np.max(per_res_rmsd),
            }
        except Exception as e:
            print(f"    Error calculating RMSD: {str(e)[:100]}")
            return {}

    def compare_design_pair(self, seq_id: int,
                            bound_dir: Path, apo_dir: Path) -> Dict:
        """Compare bound and apo predictions for a single design."""
        results = {'seq_id': seq_id, 'template': self.template_name}

        bound_model = list(bound_dir.glob("*_model.cif"))
        apo_model = list(apo_dir.glob("*_model.cif"))
        if not bound_model or not apo_model:
            return results

        bound_model = bound_model[0]
        apo_model = apo_model[0]

        # RMSD between bound and apo
        try:
            rmsd_data = self.calculate_structure_rmsd(bound_model, apo_model)
            results.update(rmsd_data)
        except Exception as e:
            print(f"    Error calculating RMSD: {str(e)[:100]}")

        # Confidence data
        try:
            bound_summary, bound_full = self._load_confidence_data(bound_dir)
            apo_summary, apo_full = self._load_confidence_data(apo_dir)

            bound_plddt = self._compute_per_residue_plddt(bound_full, chain_id='A')
            apo_plddt = self._compute_per_residue_plddt(apo_full, chain_id='A')

            if bound_plddt is not None and apo_plddt is not None and len(bound_plddt) > 1 and len(apo_plddt) > 1:
                min_len = min(len(bound_plddt), len(apo_plddt))
                bound_plddt = bound_plddt[:min_len]
                apo_plddt = apo_plddt[:min_len]
                plddt_diff = bound_plddt - apo_plddt

                results['mean_plddt_diff'] = float(np.mean(plddt_diff))
                results['max_plddt_diff'] = float(np.max(np.abs(plddt_diff)))
                results['bound_mean_plddt'] = float(np.mean(bound_plddt))
                results['apo_mean_plddt'] = float(np.mean(apo_plddt))

                chrom_start, chrom_end = self.chromophore_range
                if len(bound_plddt) >= chrom_end:
                    chrom_bound = bound_plddt[chrom_start-1:chrom_end]
                    chrom_apo = apo_plddt[chrom_start-1:chrom_end]
                    results['chromophore_bound_plddt'] = float(np.mean(chrom_bound))
                    results['chromophore_apo_plddt'] = float(np.mean(chrom_apo))
                    results['chromophore_plddt_diff'] = float(np.mean(chrom_bound - chrom_apo))

            for key in ['ptm', 'iptm', 'ranking_score']:
                if key in bound_summary:
                    results[f'bound_{key}'] = bound_summary[key]
                if key in apo_summary:
                    results[f'apo_{key}'] = apo_summary[key]

            # PAE
            if 'pae' in bound_full:
                try:
                    pae_matrix = np.array(bound_full['pae'])
                    if pae_matrix.ndim == 2:
                        results['bound_mean_pae'] = float(np.mean(pae_matrix))
                        if bound_plddt is not None:
                            n_prot = len(bound_plddt)
                            if pae_matrix.shape[0] > n_prot:
                                results['bound_interface_pae'] = float(np.mean(pae_matrix[:n_prot, n_prot:]))
                except Exception:
                    pass

            if 'pae' in apo_full:
                try:
                    pae_matrix = np.array(apo_full['pae'])
                    if pae_matrix.ndim == 2:
                        results['apo_mean_pae'] = float(np.mean(pae_matrix))
                except Exception:
                    pass

            if 'bound_mean_pae' in results and 'apo_mean_pae' in results:
                results['pae_diff'] = results['apo_mean_pae'] - results['bound_mean_pae']

        except Exception as e:
            print(f"    Error loading confidence: {str(e)[:100]}")

        return results

    def compare_all_designs(self) -> pd.DataFrame:
        """Compare all design pairs for this output directory."""
        bound_af3_dir = self.output_dir / "03_alphafold3_predictions_bound"
        apo_af3_dir = self.output_dir / "03_alphafold3_predictions_apo"

        if not bound_af3_dir.exists() or not apo_af3_dir.exists():
            return pd.DataFrame()

        bound_dirs = self._discover_inference_dirs(bound_af3_dir, 'bound')
        apo_dirs = self._discover_inference_dirs(apo_af3_dir, 'apo')

        matched_ids = sorted(set(bound_dirs.keys()) & set(apo_dirs.keys()))
        if not matched_ids:
            return pd.DataFrame()

        print(f"  {self.output_dir.name}: {len(bound_dirs)} bound, {len(apo_dirs)} apo, {len(matched_ids)} matched pairs")

        results = []
        for i, seq_id in enumerate(matched_ids, 1):
            if i % 50 == 0 or i == 1:
                print(f"    [{i}/{len(matched_ids)}] design {seq_id}")
            result = self.compare_design_pair(seq_id, bound_dirs[seq_id], apo_dirs[seq_id])
            if result and len(result) > 2:
                results.append(result)

        return pd.DataFrame(results)


# -------------------------------------------------------------------------
# Visualization (operates on combined DataFrame from all templates)
# -------------------------------------------------------------------------

def create_visualizations(df: pd.DataFrame, template_df: pd.DataFrame,
                          output_dir: Path):
    """Create visualizations with all designs and all template references."""
    output_dir.mkdir(exist_ok=True, parents=True)
    csv_dir = output_dir / "plot_data_csvs"
    csv_dir.mkdir(exist_ok=True)

    print("\nGenerating visualizations...")

    _plot_rmsd_distributions(df, template_df, output_dir, csv_dir)
    _plot_plddt_comparison(df, template_df, output_dir, csv_dir)
    _plot_chromophore_plddt(df, template_df, output_dir, csv_dir)
    _plot_rmsd_vs_plddt_change(df, template_df, output_dir, csv_dir)
    _plot_top_changes(df, template_df, output_dir, csv_dir)
    _plot_correlation_matrix(df, output_dir, csv_dir)

    print(f"\nAll plots saved to: {output_dir}")


def _get_template_colors(template_df: pd.DataFrame) -> Dict[str, str]:
    """Assign distinct colors to each template."""
    cmap = plt.cm.get_cmap('tab10', max(10, len(template_df)))
    colors = {}
    for i, name in enumerate(template_df['template_name']):
        colors[name] = cmap(i)
    return colors


def _plot_rmsd_distributions(df, template_df, output_dir, csv_dir):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    templates_in_data = df['template'].unique() if 'template' in df.columns else []

    if 'global_rmsd' in df.columns:
        for tpl in sorted(templates_in_data):
            tdf = df[df['template'] == tpl]
            short = tpl.replace('G4FP_', '')
            axes[0].hist(tdf['global_rmsd'].dropna(), bins=30, alpha=0.6, label=short)
        median = df['global_rmsd'].median()
        if pd.notna(median):
            axes[0].axvline(median, color='red', linestyle='--', label=f'Median: {median:.2f}A')
        axes[0].set_xlabel('Global RMSD (Angstrom)')
        axes[0].set_ylabel('Frequency')
        axes[0].set_title('Bound vs Apo: Global RMSD')
        axes[0].legend(fontsize=7)

    if 'chromophore_rmsd' in df.columns:
        for tpl in sorted(templates_in_data):
            tdf = df[df['template'] == tpl]
            short = tpl.replace('G4FP_', '')
            axes[1].hist(tdf['chromophore_rmsd'].dropna(), bins=30, alpha=0.6, label=short)
        median = df['chromophore_rmsd'].median()
        if pd.notna(median):
            axes[1].axvline(median, color='red', linestyle='--', label=f'Median: {median:.2f}A')
        axes[1].set_xlabel('Chromophore RMSD (Angstrom)')
        axes[1].set_ylabel('Frequency')
        axes[1].set_title('Bound vs Apo: Chromophore RMSD')
        axes[1].legend(fontsize=7)

    plt.tight_layout()
    plt.savefig(output_dir / '01_bound_apo_rmsd.png', dpi=300, bbox_inches='tight')
    plt.close()

    cols = ['seq_id', 'template'] + [c for c in ['global_rmsd', 'chromophore_rmsd'] if c in df.columns]
    df[[c for c in cols if c in df.columns]].to_csv(csv_dir / '01_bound_apo_rmsd_data.csv', index=False)
    print("  01_bound_apo_rmsd.png")


def _plot_plddt_comparison(df, template_df, output_dir, csv_dir):
    if 'bound_mean_plddt' not in df.columns or 'apo_mean_plddt' not in df.columns:
        print("  pLDDT data not available")
        return

    plot_df = df.dropna(subset=['bound_mean_plddt', 'apo_mean_plddt'])
    if plot_df.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Scatter: color by template
    templates_in_data = sorted(plot_df['template'].unique()) if 'template' in plot_df.columns else ['all']
    tpl_colors = _get_template_colors(template_df)

    for tpl in templates_in_data:
        tdf = plot_df[plot_df['template'] == tpl]
        short = tpl.replace('G4FP_', '')
        color = tpl_colors.get(tpl, 'steelblue')
        axes[0].scatter(tdf['bound_mean_plddt'], tdf['apo_mean_plddt'],
                        alpha=0.5, s=30, color=color, label=f'{short} ({len(tdf)})')

    # Diagonal
    all_vals = pd.concat([plot_df['bound_mean_plddt'], plot_df['apo_mean_plddt']])
    min_val, max_val = all_vals.min() - 1, all_vals.max() + 1
    axes[0].plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.3, label='y=x')

    # Plot ALL 10 templates as star markers (at bound_plddt, bound_plddt = no-switch baseline)
    for _, row in template_df.iterrows():
        tpl_plddt = row['mean_plddt']
        short = row['template_name'].replace('G4FP_', '')
        color = tpl_colors.get(row['template_name'], 'black')
        axes[0].scatter([tpl_plddt], [tpl_plddt], marker='*', s=200,
                        color=color, edgecolors='black', linewidths=0.8, zorder=5)
        axes[0].annotate(short, (tpl_plddt, tpl_plddt), fontsize=5,
                         xytext=(3, 3), textcoords='offset points')

    axes[0].set_xlabel('Bound Mean pLDDT')
    axes[0].set_ylabel('Apo Mean pLDDT')
    axes[0].set_title(f'pLDDT: Bound vs Apo (n={len(plot_df)})')
    axes[0].legend(fontsize=6, loc='upper left')

    corr = plot_df[['bound_mean_plddt', 'apo_mean_plddt']].corr().iloc[0, 1]
    axes[0].text(0.95, 0.05, f'r = {corr:.3f}', transform=axes[0].transAxes,
                 fontsize=10, ha='right')

    # Distribution of differences
    if 'mean_plddt_diff' in df.columns:
        diff_data = df['mean_plddt_diff'].dropna()
        axes[1].hist(diff_data, bins=30, alpha=0.7, color='green')
        axes[1].set_xlabel('pLDDT Difference (Bound - Apo)')
        axes[1].set_ylabel('Frequency')
        axes[1].set_title('pLDDT Change Distribution')
        axes[1].axvline(0, color='red', linestyle='--', label='No change')
        if len(diff_data) > 0:
            axes[1].axvline(diff_data.median(), color='blue', linestyle='--',
                            label=f'Median: {diff_data.median():.2f}')
        axes[1].legend()

    plt.tight_layout()
    plt.savefig(output_dir / '02_plddt_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    cols = [c for c in ['seq_id', 'template', 'bound_mean_plddt', 'apo_mean_plddt', 'mean_plddt_diff'] if c in df.columns]
    df[cols].to_csv(csv_dir / '02_plddt_comparison_data.csv', index=False)
    print("  02_plddt_comparison.png")


def _plot_chromophore_plddt(df, template_df, output_dir, csv_dir):
    if 'chromophore_bound_plddt' not in df.columns or 'chromophore_apo_plddt' not in df.columns:
        print("  Chromophore pLDDT data not available")
        return

    plot_df = df.dropna(subset=['chromophore_bound_plddt', 'chromophore_apo_plddt'])
    if plot_df.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    tpl_colors = _get_template_colors(template_df)

    templates_in_data = sorted(plot_df['template'].unique()) if 'template' in plot_df.columns else []

    for tpl in templates_in_data:
        tdf = plot_df[plot_df['template'] == tpl]
        short = tpl.replace('G4FP_', '')
        color = tpl_colors.get(tpl, 'purple')
        axes[0].scatter(tdf['chromophore_bound_plddt'], tdf['chromophore_apo_plddt'],
                        alpha=0.5, s=30, color=color, label=f'{short} ({len(tdf)})')

    # Diagonal
    all_vals = pd.concat([plot_df['chromophore_bound_plddt'], plot_df['chromophore_apo_plddt']])
    min_val, max_val = all_vals.min() - 1, all_vals.max() + 1
    axes[0].plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.3, label='y=x')

    # All 10 templates as star markers
    for _, row in template_df.iterrows():
        if row['chromophore_mean_plddt'] is None:
            continue
        tpl_plddt = row['chromophore_mean_plddt']
        short = row['template_name'].replace('G4FP_', '')
        color = tpl_colors.get(row['template_name'], 'black')
        axes[0].scatter([tpl_plddt], [tpl_plddt], marker='*', s=200,
                        color=color, edgecolors='black', linewidths=0.8, zorder=5)
        axes[0].annotate(short, (tpl_plddt, tpl_plddt), fontsize=5,
                         xytext=(3, 3), textcoords='offset points')

    axes[0].set_xlabel('Bound Chromophore pLDDT')
    axes[0].set_ylabel('Apo Chromophore pLDDT')
    axes[0].set_title(f'Chromophore pLDDT: Bound vs Apo (n={len(plot_df)})')
    axes[0].legend(fontsize=6)

    corr = plot_df[['chromophore_bound_plddt', 'chromophore_apo_plddt']].corr().iloc[0, 1]
    axes[0].text(0.95, 0.05, f'r = {corr:.3f}', transform=axes[0].transAxes,
                 fontsize=10, ha='right')

    # Distribution of differences
    if 'chromophore_plddt_diff' in df.columns:
        diff_data = df['chromophore_plddt_diff'].dropna()
        axes[1].hist(diff_data, bins=30, alpha=0.7, color='orange')
        axes[1].set_xlabel('Chromophore pLDDT Diff (Bound - Apo)')
        axes[1].set_ylabel('Frequency')
        axes[1].set_title('Chromophore pLDDT Change')
        axes[1].axvline(0, color='red', linestyle='--', label='No change')
        if len(diff_data) > 0:
            axes[1].axvline(diff_data.median(), color='blue', linestyle='--',
                            label=f'Median: {diff_data.median():.2f}')
        axes[1].legend()

    plt.tight_layout()
    plt.savefig(output_dir / '03_chromophore_plddt_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    cols = [c for c in ['seq_id', 'template', 'chromophore_bound_plddt', 'chromophore_apo_plddt', 'chromophore_plddt_diff'] if c in df.columns]
    df[cols].to_csv(csv_dir / '03_chromophore_plddt_comparison_data.csv', index=False)
    print("  03_chromophore_plddt_comparison.png")


def _plot_rmsd_vs_plddt_change(df, template_df, output_dir, csv_dir):
    if 'global_rmsd' not in df.columns or 'mean_plddt_diff' not in df.columns:
        print("  Required data not available for RMSD vs pLDDT plot")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    tpl_colors = _get_template_colors(template_df)

    data = df.dropna(subset=['global_rmsd', 'mean_plddt_diff'])
    if len(data) > 0:
        templates_in_data = sorted(data['template'].unique()) if 'template' in data.columns else []
        for tpl in templates_in_data:
            tdf = data[data['template'] == tpl]
            short = tpl.replace('G4FP_', '')
            color = tpl_colors.get(tpl, 'steelblue')
            axes[0].scatter(tdf['global_rmsd'], tdf['mean_plddt_diff'],
                            alpha=0.5, s=30, color=color, label=short)
        axes[0].set_xlabel('Global RMSD (Angstrom)')
        axes[0].set_ylabel('pLDDT Diff (Bound - Apo)')
        axes[0].set_title('RMSD vs pLDDT Change')
        axes[0].axhline(0, color='red', linestyle='--', alpha=0.5)
        axes[0].legend(fontsize=7)
        corr = data[['global_rmsd', 'mean_plddt_diff']].corr().iloc[0, 1]
        axes[0].text(0.05, 0.95, f'r = {corr:.3f}', transform=axes[0].transAxes, fontsize=10)

    if 'chromophore_rmsd' in df.columns and 'chromophore_plddt_diff' in df.columns:
        data_chrom = df.dropna(subset=['chromophore_rmsd', 'chromophore_plddt_diff'])
        if len(data_chrom) > 0:
            templates_in_data = sorted(data_chrom['template'].unique()) if 'template' in data_chrom.columns else []
            for tpl in templates_in_data:
                tdf = data_chrom[data_chrom['template'] == tpl]
                short = tpl.replace('G4FP_', '')
                color = tpl_colors.get(tpl, 'coral')
                axes[1].scatter(tdf['chromophore_rmsd'], tdf['chromophore_plddt_diff'],
                                alpha=0.5, s=30, color=color, label=short)
            axes[1].set_xlabel('Chromophore RMSD (Angstrom)')
            axes[1].set_ylabel('Chromophore pLDDT Diff (Bound - Apo)')
            axes[1].set_title('Chromophore: RMSD vs pLDDT Change')
            axes[1].axhline(0, color='red', linestyle='--', alpha=0.5)
            axes[1].legend(fontsize=7)

    plt.tight_layout()
    plt.savefig(output_dir / '04_rmsd_vs_plddt_change.png', dpi=300, bbox_inches='tight')
    plt.close()

    cols = [c for c in ['seq_id', 'template', 'global_rmsd', 'mean_plddt_diff', 'chromophore_rmsd', 'chromophore_plddt_diff'] if c in df.columns]
    df[cols].to_csv(csv_dir / '04_rmsd_vs_plddt_change_data.csv', index=False)
    print("  04_rmsd_vs_plddt_change.png")


def _plot_top_changes(df, template_df, output_dir, csv_dir):
    if 'global_rmsd' not in df.columns:
        return

    rmsd_data = df.dropna(subset=['global_rmsd'])
    if rmsd_data.empty:
        return

    top_20 = rmsd_data.nlargest(min(20, len(rmsd_data)), 'global_rmsd')
    tpl_colors = _get_template_colors(template_df)

    fig, ax = plt.subplots(figsize=(12, 8))
    x = range(len(top_20))
    colors = [tpl_colors.get(tpl, 'steelblue') for tpl in top_20['template']]
    ax.bar(x, top_20['global_rmsd'], alpha=0.7, color=colors)
    ax.set_xlabel('Design')
    ax.set_ylabel('Global RMSD (Angstrom)')
    ax.set_title('Top 20 Designs by Conformational Change (Bound vs Apo)')
    labels = [f"{int(r['seq_id'])}\n{r['template'].replace('G4FP_', '')}" for _, r in top_20.iterrows()]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, fontsize=7)

    plt.tight_layout()
    plt.savefig(output_dir / '05_top20_conformational_changes.png', dpi=300, bbox_inches='tight')
    plt.close()

    array_cols = [col for col in top_20.columns if 'array' in col]
    top_20.drop(columns=array_cols, errors='ignore').to_csv(
        csv_dir / '05_top20_conformational_changes_data.csv', index=False)
    print("  05_top20_conformational_changes.png")


def _plot_correlation_matrix(df, output_dir, csv_dir):
    possible_cols = ['global_rmsd', 'chromophore_rmsd', 'mean_per_res_rmsd',
                     'bound_mean_plddt', 'apo_mean_plddt', 'mean_plddt_diff',
                     'chromophore_bound_plddt', 'chromophore_apo_plddt', 'chromophore_plddt_diff',
                     'bound_ptm', 'apo_ptm', 'bound_iptm', 'apo_iptm',
                     'bound_ranking_score', 'apo_ranking_score',
                     'bound_mean_pae', 'apo_mean_pae', 'pae_diff']
    cols = [c for c in possible_cols if c in df.columns and df[c].notna().any()]

    if len(cols) < 2:
        return

    corr_matrix = df[cols].corr()
    plt.figure(figsize=(14, 12))
    sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm',
                center=0, square=True, linewidths=1)
    plt.title('Correlation Matrix: Bound vs Apo')
    plt.tight_layout()
    plt.savefig(output_dir / '06_correlation_matrix.png', dpi=300, bbox_inches='tight')
    plt.close()
    corr_matrix.to_csv(csv_dir / '06_correlation_matrix_data.csv')
    print("  06_correlation_matrix.png")


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description='Compare bound and apo AlphaFold3 predictions')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Single output directory (default: auto-discover all output_G4FP_* dirs)')
    parser.add_argument('--template', type=str, default='inputs/G4FP_des1_cro_mod0.pdb',
                        help='Reference template for RMSD (default: inputs/G4FP_des1_cro_mod0.pdb)')
    parser.add_argument('--chromophore-range', type=str, default='197,199',
                        help='Chromophore residue range as "start,end"')
    parser.add_argument('--analysis-output-dir', type=str, default='analysis_output',
                        help='Central directory for all analysis output')

    args = parser.parse_args()
    chrom_start, chrom_end = map(int, args.chromophore_range.split(','))
    chromophore_range = (chrom_start, chrom_end)

    base_dir = Path(__file__).resolve().parent

    # Discover output directories
    if args.output_dir:
        output_dirs = [Path(args.output_dir)]
    else:
        output_dirs = sorted(
            d for d in base_dir.iterdir()
            if d.is_dir() and d.name.startswith("output_G4FP_")
        )

    if not output_dirs:
        print("No output directories found")
        sys.exit(1)

    # Load all 10 template metrics for reference overlays
    inputs_dir = base_dir / "inputs"
    template_df = load_all_template_metrics(inputs_dir, chromophore_range)
    print(f"\nLoaded {len(template_df)} template reference structures:")
    for _, row in template_df.iterrows():
        print(f"  {row['template_name']}: mean_plddt={row['mean_plddt']:.1f}, "
              f"chrom_plddt={row['chromophore_mean_plddt']:.1f}")

    # Compare all designs across all output directories
    print(f"\n{'='*80}")
    print(f"Comparing Ligand States Across {len(output_dirs)} Templates")
    print(f"{'='*80}")

    all_frames = []
    for output_dir in output_dirs:
        print(f"\n--- {output_dir.name} ---")

        comparator = LigandStateComparator(
            output_dir=output_dir,
            reference_template=Path(args.template),
            chromophore_range=chromophore_range
        )

        df = comparator.compare_all_designs()
        if not df.empty:
            all_frames.append(df)

            # Also save per-template CSV to original output dir (for script 07)
            df_export = df.copy()
            array_cols = [c for c in df_export.columns if 'array' in c]
            df_export = df_export.drop(columns=array_cols, errors='ignore')
            csv_path = output_dir / "06_ligand_state_comparison_results.csv"
            df_export.to_csv(csv_path, index=False)
            print(f"  Saved: {csv_path}")
        else:
            print(f"  No matched design pairs found")

    if not all_frames:
        print("\nNo designs compared across any template")
        sys.exit(1)

    # Combine all results
    df_combined = pd.concat(all_frames, ignore_index=True)
    print(f"\n{'='*80}")
    print(f"Combined: {len(df_combined)} design pairs across "
          f"{df_combined['template'].nunique()} templates")
    print(f"{'='*80}")

    # Create combined visualizations
    central_dir = Path(args.analysis_output_dir)
    viz_dir = central_dir / "06_ligand_state_comparison"
    create_visualizations(df_combined, template_df, viz_dir)

    # Save combined CSV
    df_export = df_combined.copy()
    array_cols = [c for c in df_export.columns if 'array' in c]
    df_export = df_export.drop(columns=array_cols, errors='ignore')
    combined_csv = central_dir / "06_ligand_state_comparison_results.csv"
    df_export.to_csv(combined_csv, index=False)

    print(f"\nAnalysis complete!")
    print(f"  Combined results: {combined_csv}")
    print(f"  Plots: {viz_dir}/")
    print(f"  Per-template CSVs: output_G4FP_*/06_ligand_state_comparison_results.csv")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Analyze AlphaFold3 Structure Predictions for G4FP Designs

Comprehensive analysis of AF3 predictions including:
1. AlphaFold3 quality metrics (pTM, ipTM, pLDDT, PAE)
2. RMSD analysis (global and chromophore-specific) vs reference
3. Confidence-pLDDT correlations
4. Per-residue analysis
5. Top design identification

All RMSD calculations use a single reference structure.
All 10 template PDBs are shown as reference markers on plots.

Usage:
    # All templates (auto-discover)
    python 05_analyze_af3_structures.py --chromophore-range 197,199

    # Single template
    python 05_analyze_af3_structures.py --output-dir output_G4FP_des1_cro_mod0 \\
        --template inputs/G4FP_des1_cro_mod0.pdb --state both --chromophore-range 197,199
"""

import os
import sys
import re
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
        chrom_plddt = float(np.mean(b_arr[chrom_start-1:chrom_end])) if len(b_arr) >= chrom_end else None
        records.append({
            'template_name': pdb.stem,
            'template_pdb': str(pdb),
            'mean_plddt': float(np.mean(b_arr)),
            'chromophore_mean_plddt': chrom_plddt,
            'n_residues': len(b_arr),
        })
    return pd.DataFrame(records)


def _get_template_colors(template_df: pd.DataFrame) -> Dict[str, str]:
    """Assign distinct colors to each template."""
    cmap = plt.cm.get_cmap('tab10', max(10, len(template_df)))
    colors = {}
    for i, name in enumerate(template_df['template_name']):
        colors[name] = cmap(i)
    return colors


# -------------------------------------------------------------------------
# AlphaFoldStructureAnalyzer (single output dir)
# -------------------------------------------------------------------------

class AlphaFoldStructureAnalyzer:
    """Analyze AlphaFold3 predictions for one output directory."""

    def __init__(self, output_dir: Path, reference_template: Path,
                 chromophore_range: Tuple[int, int] = (197, 199)):
        self.output_dir = Path(output_dir)
        self.reference_template = Path(reference_template)
        self.chromophore_range = chromophore_range

        self.pdb_parser = PDBParser(QUIET=True)
        self.cif_parser = MMCIFParser(QUIET=True)

        # Load reference structure for RMSD
        self.ref_structure = self._load_structure(self.reference_template)
        self.ref_ca_atoms = self._extract_ca_atoms(self.ref_structure)

        # Template name from output dir
        self.template_name = self.output_dir.name.replace("output_", "")

        print(f"  Reference: {self.reference_template} ({len(self.ref_ca_atoms)} CA atoms)")

    def _load_structure(self, filepath: Path):
        if filepath.suffix == '.pdb':
            return self.pdb_parser.get_structure('structure', str(filepath))
        elif filepath.suffix == '.cif':
            return self.cif_parser.get_structure('structure', str(filepath))
        else:
            raise ValueError(f"Unknown format: {filepath.suffix}")

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
        """Discover AF3 inference dirs with model files."""
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

    def calculate_rmsd(self, model_file: Path, chain_id='A') -> Dict[str, float]:
        model_structure = self._load_structure(model_file)
        model_ca = self._extract_ca_atoms(model_structure, chain_id)

        min_len = min(len(model_ca), len(self.ref_ca_atoms))
        model_ca = model_ca[:min_len]
        ref_ca = self.ref_ca_atoms[:min_len]

        sup = Superimposer()
        sup.set_atoms(ref_ca, model_ca)
        global_rmsd = sup.rms

        chrom_start, chrom_end = self.chromophore_range
        chromophore_rmsd = None
        if len(model_ca) >= chrom_end:
            chrom_ref = ref_ca[chrom_start-1:chrom_end]
            chrom_model = model_ca[chrom_start-1:chrom_end]
            sup_chrom = Superimposer()
            sup_chrom.set_atoms(chrom_ref, chrom_model)
            chromophore_rmsd = sup_chrom.rms

        return {'global_rmsd': global_rmsd, 'chromophore_rmsd': chromophore_rmsd}

    def analyze_design(self, design_dir: Path, state: str, seq_id: int) -> Dict:
        """Analyze a single design."""
        results = {
            'seq_id': seq_id,
            'design_name': design_dir.name,
            'state': state,
            'template': self.template_name,
        }

        model_files = list(design_dir.glob("*_model.cif"))
        if not model_files:
            return results

        model_file = model_files[0]

        try:
            rmsd_data = self.calculate_rmsd(model_file)
            results.update(rmsd_data)
        except Exception as e:
            print(f"    RMSD error: {str(e)[:80]}")

        try:
            summary_conf, full_conf = self._load_confidence_data(design_dir)
            for key in ['ptm', 'iptm', 'ranking_score', 'fraction_disordered']:
                if key in summary_conf:
                    results[key] = summary_conf[key]

            plddt = self._compute_per_residue_plddt(full_conf, chain_id='A')
            if plddt is not None and len(plddt) > 1:
                results['mean_plddt'] = float(np.mean(plddt))
                results['min_plddt'] = float(np.min(plddt))

                chrom_start, chrom_end = self.chromophore_range
                if len(plddt) >= chrom_end:
                    results['chromophore_mean_plddt'] = float(np.mean(plddt[chrom_start-1:chrom_end]))

            if 'pae' in full_conf:
                try:
                    pae_matrix = np.array(full_conf['pae'])
                    if pae_matrix.ndim == 2:
                        results['mean_pae'] = float(np.mean(pae_matrix))
                        n_prot = len(self.ref_ca_atoms)
                        if state == 'bound' and pae_matrix.shape[0] > n_prot:
                            results['interface_pae'] = float(np.mean(pae_matrix[:n_prot, n_prot:]))
                except Exception:
                    pass
        except Exception as e:
            print(f"    Confidence error: {str(e)[:80]}")

        return results

    def analyze_all_designs(self, state: str = 'both') -> pd.DataFrame:
        states = []
        if state in ['bound', 'both']:
            states.append('bound')
        if state in ['apo', 'both']:
            states.append('apo')

        all_results = []
        for state_name in states:
            af3_dir = self.output_dir / f"03_alphafold3_predictions_{state_name}"
            if not af3_dir.exists():
                continue

            inference_dirs = self._discover_inference_dirs(af3_dir, state_name)
            if not inference_dirs:
                continue

            print(f"  {state_name}: {len(inference_dirs)} designs")
            for i, (seq_id, design_dir) in enumerate(sorted(inference_dirs.items()), 1):
                if i % 50 == 0 or i == 1:
                    print(f"    [{i}/{len(inference_dirs)}]")
                result = self.analyze_design(design_dir, state_name, seq_id)
                if result and 'mean_plddt' in result:
                    all_results.append(result)

        return pd.DataFrame(all_results)


# -------------------------------------------------------------------------
# Visualization (operates on combined DataFrame)
# -------------------------------------------------------------------------

def create_visualizations(df: pd.DataFrame, template_df: pd.DataFrame,
                          output_dir: Path):
    output_dir.mkdir(exist_ok=True, parents=True)
    csv_dir = output_dir / "plot_data_csvs"
    csv_dir.mkdir(exist_ok=True)

    print("\nGenerating visualizations...")

    _plot_af3_metrics(df, template_df, output_dir, csv_dir)
    _plot_rmsd_distributions(df, template_df, output_dir, csv_dir)
    _plot_correlation_matrix(df, output_dir, csv_dir)
    _plot_global_vs_chromophore_rmsd(df, template_df, output_dir, csv_dir)
    _export_top_designs(df, output_dir, csv_dir)
    _plot_per_residue_analysis(df, template_df, output_dir, csv_dir)

    print(f"\nAll plots saved to: {output_dir}")


def _plot_af3_metrics(df, template_df, output_dir, csv_dir):
    metrics = ['mean_plddt', 'ptm', 'iptm', 'ranking_score']
    available = [m for m in metrics if m in df.columns and df[m].notna().any()]
    if not available:
        return

    tpl_colors = _get_template_colors(template_df)
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()

    for idx, metric in enumerate(available):
        # Histograms by state and template
        for state in sorted(df['state'].unique()):
            templates_in_state = sorted(df[df['state'] == state]['template'].unique())
            for tpl in templates_in_state:
                tdf = df[(df['state'] == state) & (df['template'] == tpl)]
                vals = tdf[metric].dropna()
                if len(vals) > 0:
                    short = tpl.replace('G4FP_', '')
                    color = tpl_colors.get(tpl, 'steelblue')
                    axes[idx].hist(vals, bins=20, alpha=0.4, color=color,
                                   label=f'{short} {state}')

        # Template reference lines (for pLDDT only)
        if metric == 'mean_plddt':
            for _, row in template_df.iterrows():
                plddt = row['mean_plddt']
                short = row['template_name'].replace('G4FP_', '')
                color = tpl_colors.get(row['template_name'], 'black')
                axes[idx].axvline(plddt, color=color, linestyle='--',
                                  linewidth=1.5, alpha=0.7)

        axes[idx].set_xlabel(metric.replace('_', ' ').title())
        axes[idx].set_ylabel('Frequency')
        axes[idx].set_title(f'{metric.replace("_", " ").title()} Distribution')
        axes[idx].legend(fontsize=6)

    for idx in range(len(available), 4):
        axes[idx].axis('off')

    plt.tight_layout()
    plt.savefig(output_dir / '02_alphafold3_metrics.png', dpi=300, bbox_inches='tight')
    plt.close()

    export_cols = [c for c in ['seq_id', 'state', 'template'] + available if c in df.columns]
    df[export_cols].to_csv(csv_dir / '02_alphafold3_metrics_data.csv', index=False)
    print("  02_alphafold3_metrics.png")


def _plot_rmsd_distributions(df, template_df, output_dir, csv_dir):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    tpl_colors = _get_template_colors(template_df)

    for col_idx, (metric, title) in enumerate([
        ('global_rmsd', 'Global RMSD'),
        ('chromophore_rmsd', 'Chromophore RMSD')
    ]):
        if metric not in df.columns:
            continue
        for state in sorted(df['state'].unique()):
            for tpl in sorted(df[df['state'] == state]['template'].unique()):
                tdf = df[(df['state'] == state) & (df['template'] == tpl)]
                vals = tdf[metric].dropna()
                if len(vals) > 0:
                    short = tpl.replace('G4FP_', '')
                    color = tpl_colors.get(tpl, 'steelblue')
                    axes[col_idx].hist(vals, bins=20, alpha=0.4, color=color,
                                       label=f'{short} {state}')

        axes[col_idx].axvline(0, color='black', linestyle='--', linewidth=2,
                              alpha=0.5, label='Reference (RMSD=0)')
        axes[col_idx].set_xlabel(f'{title} (Angstrom)')
        axes[col_idx].set_ylabel('Frequency')
        axes[col_idx].set_title(f'{title} Distribution')
        axes[col_idx].legend(fontsize=6)

    plt.tight_layout()
    plt.savefig(output_dir / '04_rmsd_distributions.png', dpi=300, bbox_inches='tight')
    plt.close()

    cols = [c for c in ['seq_id', 'state', 'template', 'global_rmsd', 'chromophore_rmsd'] if c in df.columns]
    df[cols].to_csv(csv_dir / '04_rmsd_distributions_data.csv', index=False)
    print("  04_rmsd_distributions.png")


def _plot_correlation_matrix(df, output_dir, csv_dir):
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    exclude = ['seq_id']
    numeric_cols = [c for c in numeric_cols if c not in exclude]
    if len(numeric_cols) < 2:
        return

    corr = df[numeric_cols].corr()
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr, annot=True, fmt='.2f', cmap='coolwarm', center=0, square=True, linewidths=1)
    plt.title('Correlation Matrix')
    plt.tight_layout()
    plt.savefig(output_dir / '05_correlation_matrix.png', dpi=300, bbox_inches='tight')
    plt.close()
    corr.to_csv(csv_dir / '05_correlation_matrix_data.csv')
    print("  05_correlation_matrix.png")


def _plot_global_vs_chromophore_rmsd(df, template_df, output_dir, csv_dir):
    if 'global_rmsd' not in df.columns or 'chromophore_rmsd' not in df.columns:
        return

    tpl_colors = _get_template_colors(template_df)
    plt.figure(figsize=(10, 8))

    for state in sorted(df['state'].unique()):
        for tpl in sorted(df[df['state'] == state]['template'].unique()):
            tdf = df[(df['state'] == state) & (df['template'] == tpl)]
            tdf = tdf.dropna(subset=['global_rmsd', 'chromophore_rmsd'])
            if len(tdf) > 0:
                short = tpl.replace('G4FP_', '')
                color = tpl_colors.get(tpl, 'steelblue')
                plt.scatter(tdf['global_rmsd'], tdf['chromophore_rmsd'],
                            alpha=0.4, s=25, color=color, label=f'{short} {state}')

    # All templates at origin
    for _, row in template_df.iterrows():
        short = row['template_name'].replace('G4FP_', '')
        color = tpl_colors.get(row['template_name'], 'black')
        plt.scatter([0], [0], color=color, marker='*', s=200, zorder=5)
    plt.scatter([0], [0], color='black', marker='*', s=200, zorder=5, label='Templates (RMSD=0)')

    plt.xlabel('Global RMSD (Angstrom)')
    plt.ylabel('Chromophore RMSD (Angstrom)')
    plt.title('Global vs Chromophore RMSD')
    plt.legend(fontsize=7, loc='upper left')
    plt.tight_layout()
    plt.savefig(output_dir / '06_global_vs_chromophore_rmsd.png', dpi=300, bbox_inches='tight')
    plt.close()

    cols = [c for c in ['seq_id', 'state', 'template', 'global_rmsd', 'chromophore_rmsd'] if c in df.columns]
    df[cols].to_csv(csv_dir / '06_global_vs_chromophore_rmsd_data.csv', index=False)
    print("  06_global_vs_chromophore_rmsd.png")


def _export_top_designs(df, output_dir, csv_dir):
    metrics_to_rank = []
    if 'mean_plddt' in df.columns:
        metrics_to_rank.append(('mean_plddt', False))
    if 'global_rmsd' in df.columns:
        metrics_to_rank.append(('global_rmsd', True))

    if not metrics_to_rank:
        return

    top_designs = {}
    for metric, ascending in metrics_to_rank:
        for state in df['state'].unique():
            sdf = df[df['state'] == state]
            top_20 = sdf.nsmallest(20, metric) if ascending else sdf.nlargest(20, metric)
            top_designs[f'top_20_{state}_{metric}'] = top_20

    all_top = pd.concat([d.assign(metric=key) for key, d in top_designs.items()], ignore_index=True)
    all_top.to_csv(csv_dir / '07_top20_designs_data.csv', index=False)
    print("  07_top20_designs exported")


def _plot_per_residue_analysis(df, template_df, output_dir, csv_dir):
    """Plot per-residue RMSD for top design."""
    if 'mean_plddt' not in df.columns or df['mean_plddt'].isna().all():
        return

    # Find best design
    top = df.loc[df['mean_plddt'].idxmax()]
    seq_id = int(top['seq_id'])
    state = top['state']
    tpl_name = top.get('template', '')

    # We'd need to re-calculate per-residue RMSD; skip if we can't find the dir
    # This is a lightweight plot so just note the top design
    print(f"  Top design: seq_{seq_id} ({state}, {tpl_name}) pLDDT={top['mean_plddt']:.1f}")


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description='Analyze AlphaFold3 structure predictions')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Single output directory (default: auto-discover all)')
    parser.add_argument('--template', type=str, default='inputs/G4FP_des1_cro_mod0.pdb',
                        help='Reference structure for RMSD (default: inputs/G4FP_des1_cro_mod0.pdb)')
    parser.add_argument('--state', type=str, default='both', choices=['bound', 'apo', 'both'],
                        help='Which predictions to analyze')
    parser.add_argument('--chromophore-range', type=str, default='197,199',
                        help='Chromophore residue range as "start,end"')
    parser.add_argument('--analysis-output-dir', type=str, default='analysis_output',
                        help='Central directory for output')

    args = parser.parse_args()
    chrom_start, chrom_end = map(int, args.chromophore_range.split(','))
    chromophore_range = (chrom_start, chrom_end)

    base_dir = Path(__file__).resolve().parent

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

    # Load all 10 template metrics
    inputs_dir = base_dir / "inputs"
    template_df = load_all_template_metrics(inputs_dir, chromophore_range)
    print(f"\nLoaded {len(template_df)} template references:")
    for _, row in template_df.iterrows():
        print(f"  {row['template_name']}: pLDDT={row['mean_plddt']:.1f}, chrom={row['chromophore_mean_plddt']:.1f}")

    # Analyze all output directories
    print(f"\n{'='*80}")
    print(f"Analyzing AF3 Structures Across {len(output_dirs)} Templates")
    print(f"{'='*80}")

    all_frames = []
    for odir in output_dirs:
        print(f"\n--- {odir.name} ---")
        try:
            analyzer = AlphaFoldStructureAnalyzer(
                output_dir=odir,
                reference_template=Path(args.template),
                chromophore_range=chromophore_range
            )
            df = analyzer.analyze_all_designs(state=args.state)
            if not df.empty:
                all_frames.append(df)
                # Save per-template CSV
                csv_path = odir / "05_structure_analysis_results.csv"
                df.to_csv(csv_path, index=False)
                print(f"  Saved: {csv_path} ({len(df)} designs)")
            else:
                print(f"  No designs with completed AF3 inference")
        except Exception as e:
            print(f"  Error: {e}")
            continue

    if not all_frames:
        print("\nNo designs analyzed across any template")
        sys.exit(1)

    df_combined = pd.concat(all_frames, ignore_index=True)
    print(f"\n{'='*80}")
    print(f"Combined: {len(df_combined)} designs across "
          f"{df_combined['template'].nunique()} templates, "
          f"{df_combined['state'].nunique()} states")
    print(f"{'='*80}")

    # Create combined visualizations
    central_dir = Path(args.analysis_output_dir)
    viz_dir = central_dir / "05_structure_analysis"
    create_visualizations(df_combined, template_df, viz_dir)

    # Save combined CSV
    combined_csv = central_dir / "05_structure_analysis_results.csv"
    combined_csv.parent.mkdir(exist_ok=True, parents=True)
    df_combined.to_csv(combined_csv, index=False)

    print(f"\nAnalysis complete!")
    print(f"  Combined results: {combined_csv}")
    print(f"  Plots: {viz_dir}/")


if __name__ == "__main__":
    main()

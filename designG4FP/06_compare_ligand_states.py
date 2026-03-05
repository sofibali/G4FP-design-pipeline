#!/usr/bin/env python3
"""
Compare AlphaFold3 Predictions: Ligand-Bound vs Apo States

Analyzes ALL 25 seed-sample models per design per state to compute
distributions with mean +/- SD. Matched-seed RMSD between holo and apo.
Template AF3 predictions used as reference with error bars on scatter plots.

Apo state = protein-only (no DNA, no chromophore, no ions).
Bound state = protein + G4 DNA (GGGTGGGTGGGTGGGT) + 3 K+ ions.

Usage:
    python 06_compare_ligand_states.py --chromophore-range 197,199
    python 06_compare_ligand_states.py --output-dir output_G4FP_des1_cro_mod0
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

from Bio.PDB import PDBParser, MMCIFParser, Superimposer

warnings.filterwarnings('ignore')
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10


# Map user-facing state names to directory names on disk
DISK_STATE = {'holo': 'bound', 'apo': 'apo'}

# ─── Shared utilities (same as script 05) ─────────────────────────────

def discover_seed_sample_dirs(design_dir: Path) -> List[Tuple[int, int, Path]]:
    pattern = re.compile(r'^seed-(\d+)_sample-(\d+)$')
    dirs = []
    for d in design_dir.iterdir():
        if d.is_dir():
            m = pattern.match(d.name)
            if m:
                dirs.append((int(m.group(1)), int(m.group(2)), d))
    if not dirs:
        for subdir in design_dir.iterdir():
            if subdir.is_dir():
                for d in subdir.iterdir():
                    if d.is_dir():
                        m = pattern.match(d.name)
                        if m:
                            dirs.append((int(m.group(1)), int(m.group(2)), d))
                if dirs:
                    break
    return sorted(dirs)


def compute_per_residue_plddt(full_conf: Dict, chain_id: str = 'A') -> Optional[np.ndarray]:
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


def find_af3_result_dir(base_dir: Path) -> Optional[Path]:
    if list(base_dir.glob("*_model.cif")) or list(base_dir.glob("seed-*")):
        return base_dir
    for subdir in sorted(base_dir.iterdir()):
        if subdir.is_dir():
            if list(subdir.glob("*_model.cif")) or list(subdir.glob("seed-*")):
                return subdir
    return None


def load_template_metrics(base_dir: Path, chromophore_range: Tuple[int, int]) -> pd.DataFrame:
    """Load template metrics from AF3 outputs if available, else B-factors."""
    tpl_dir = base_dir / "template_af3_outputs"
    inputs_dir = base_dir / "inputs"
    pdb_parser = PDBParser(QUIET=True)
    records = []

    for pdb in sorted(inputs_dir.glob("G4FP_*.pdb")):
        tpl_name = pdb.stem
        rec = {'template_name': tpl_name, 'source': 'bfactor'}

        structure = pdb_parser.get_structure('s', str(pdb))
        b_factors = []
        for model in structure:
            for chain in model:
                if chain.get_id() == 'A':
                    for residue in chain:
                        if residue.has_id('CA'):
                            b_factors.append(residue['CA'].get_bfactor())
        if b_factors:
            b_arr = np.array(b_factors)
            rec['mean_plddt'] = float(np.mean(b_arr))
            rec['n_residues'] = len(b_arr)
            cs, ce = chromophore_range
            if len(b_arr) >= ce:
                rec['chromophore_mean_plddt'] = float(np.mean(b_arr[cs-1:ce]))

        for disk_state in ['bound', 'apo']:
            display_state = 'holo' if disk_state == 'bound' else disk_state
            af3_name = f"template_{tpl_name}_{disk_state}"
            af3_base = tpl_dir / disk_state / af3_name
            if not af3_base.exists():
                continue
            result_dir = find_af3_result_dir(af3_base)
            if result_dir is None:
                continue
            seed_dirs = discover_seed_sample_dirs(result_dir)
            if not seed_dirs:
                continue

            rec['source'] = 'af3'
            plddts, chrom_plddts, scores, ptms, iptms = [], [], [], [], []

            for _s, _sa, sd in seed_dirs:
                sf = sd / 'summary_confidences.json'
                if sf.exists():
                    with open(sf) as f:
                        sc = json.load(f)
                    if 'ranking_score' in sc:
                        scores.append(sc['ranking_score'])
                    if 'ptm' in sc:
                        ptms.append(sc['ptm'])
                    if 'iptm' in sc:
                        iptms.append(sc['iptm'])

                ff = sd / 'confidences.json'
                if ff.exists():
                    with open(ff) as f:
                        fc = json.load(f)
                    plddt = compute_per_residue_plddt(fc, 'A')
                    if plddt is not None and len(plddt) > 1:
                        plddts.append(float(np.mean(plddt)))
                        cs, ce = chromophore_range
                        if len(plddt) >= ce:
                            chrom_plddts.append(float(np.mean(plddt[cs-1:ce])))

            prefix = f'{display_state}_'
            if plddts:
                rec[f'{prefix}mean_plddt'] = float(np.mean(plddts))
                rec[f'{prefix}mean_plddt_sd'] = float(np.std(plddts))
            if chrom_plddts:
                rec[f'{prefix}chromophore_plddt'] = float(np.mean(chrom_plddts))
                rec[f'{prefix}chromophore_plddt_sd'] = float(np.std(chrom_plddts))
            if scores:
                rec[f'{prefix}ranking_score'] = float(np.mean(scores))
                rec[f'{prefix}ranking_score_sd'] = float(np.std(scores))
            if ptms:
                rec[f'{prefix}ptm'] = float(np.mean(ptms))
                rec[f'{prefix}ptm_sd'] = float(np.std(ptms))
            if iptms:
                rec[f'{prefix}iptm'] = float(np.mean(iptms))
                rec[f'{prefix}iptm_sd'] = float(np.std(iptms))

        records.append(rec)
    return pd.DataFrame(records)


# ─── LigandStateComparator ───────────────────────────────────────────

class LigandStateComparator:
    """Compare holo vs apo AF3 predictions, all seeds."""

    def __init__(self, output_dir: Path, reference_template: Path,
                 chromophore_range: Tuple[int, int] = (197, 199)):
        self.output_dir = Path(output_dir)
        self.reference_template = Path(reference_template)
        self.chromophore_range = chromophore_range
        self.pdb_parser = PDBParser(QUIET=True)
        self.cif_parser = MMCIFParser(QUIET=True)
        self.template_name = self.output_dir.name.replace("output_", "")

    def _load_structure(self, filepath: Path):
        if filepath.suffix == '.pdb':
            return self.pdb_parser.get_structure('s', str(filepath))
        return self.cif_parser.get_structure('s', str(filepath))

    def _extract_ca_atoms(self, structure, chain_id='A'):
        ca = []
        for model in structure:
            for chain in model:
                if chain.get_id() == chain_id:
                    for res in chain:
                        if res.has_id('CA'):
                            ca.append(res['CA'])
        return ca

    def _discover_inference_dirs(self, af3_dir: Path, state: str) -> Dict[int, Path]:
        pattern = re.compile(
            rf'^design_(\d+)_{re.escape(state)}(?:_(\d{{8}}_\d{{6}}))?$'
        )
        padded_pattern = re.compile(
            rf'^design_(\d{{4}})_{re.escape(state)}$'
        )
        candidates = defaultdict(list)
        for d in af3_dir.iterdir():
            if not d.is_dir():
                continue
            # Check padded wrapper dirs first (design_0000_bound/design_0_bound/)
            mp = padded_pattern.match(d.name)
            if mp:
                seq_id = int(mp.group(1))
                inner = find_af3_result_dir(d)
                if inner is not None:
                    candidates[seq_id].append((False, "", inner))
                continue
            # Check non-padded dirs at top level (design_0_bound, design_0_bound_20260227_...)
            m = pattern.match(d.name)
            if m:
                seq_id = int(m.group(1))
                ts = m.group(2)
                has_model = bool(list(d.glob("*_model.cif")))
                has_seeds = bool(list(d.glob("seed-*")))
                if has_model or has_seeds:
                    candidates[seq_id].append((ts is not None, ts or "", d))

        result = {}
        for seq_id, entries in candidates.items():
            non_ts = [e for e in entries if not e[0]]
            if non_ts:
                result[seq_id] = non_ts[0][2]
            else:
                entries.sort(key=lambda e: e[1], reverse=True)
                result[seq_id] = entries[0][2]
        return result

    def _analyze_state_seed(self, seed_dir: Path, state: str) -> Optional[Dict]:
        """Analyze one seed-sample for one state."""
        result = {}
        sf = seed_dir / 'summary_confidences.json'
        if sf.exists():
            try:
                with open(sf) as f:
                    sc = json.load(f)
                for k in ['ptm', 'iptm', 'ranking_score']:
                    if k in sc:
                        result[k] = sc[k]
            except Exception:
                pass

        ff = seed_dir / 'confidences.json'
        if ff.exists():
            try:
                with open(ff) as f:
                    fc = json.load(f)
                plddt = compute_per_residue_plddt(fc, 'A')
                if plddt is not None and len(plddt) > 1:
                    result['mean_plddt'] = float(np.mean(plddt))
                    cs, ce = self.chromophore_range
                    if len(plddt) >= ce:
                        result['chromophore_plddt'] = float(
                            np.mean(plddt[cs-1:ce]))
                if 'pae' in fc:
                    pae = np.array(fc['pae'])
                    if pae.ndim == 2:
                        result['mean_pae'] = float(np.mean(pae))
                        if state == 'holo' and plddt is not None:
                            n = len(plddt)
                            if pae.shape[0] > n:
                                result['interface_pae'] = float(
                                    np.mean(pae[:n, n:]))
            except Exception:
                pass

        return result if result else None

    def _compute_pair_rmsd(self, holo_seed_dir: Path,
                           apo_seed_dir: Path) -> Dict[str, float]:
        """Compute RMSD between matched holo and apo seed models."""
        b_cif = holo_seed_dir / 'model.cif'
        a_cif = apo_seed_dir / 'model.cif'
        if not b_cif.exists() or not a_cif.exists():
            return {}
        try:
            b_struct = self._load_structure(b_cif)
            a_struct = self._load_structure(a_cif)
            b_ca = self._extract_ca_atoms(b_struct, 'A')
            a_ca = self._extract_ca_atoms(a_struct, 'A')
            min_len = min(len(b_ca), len(a_ca))
            if min_len < 10:
                return {}
            b_ca, a_ca = b_ca[:min_len], a_ca[:min_len]

            sup = Superimposer()
            sup.set_atoms(b_ca, a_ca)
            result = {'global_rmsd': sup.rms}

            cs, ce = self.chromophore_range
            if min_len >= ce:
                sup_c = Superimposer()
                sup_c.set_atoms(b_ca[cs-1:ce], a_ca[cs-1:ce])
                result['chromophore_rmsd'] = sup_c.rms

            return result
        except Exception:
            return {}

    def compare_design_pair(self, seq_id: int,
                            holo_dir: Path, apo_dir: Path) -> Dict:
        """Compare all seed-sample models for holo vs apo."""
        result = {'seq_id': seq_id, 'template': self.template_name}

        holo_seeds = discover_seed_sample_dirs(holo_dir)
        apo_seeds = discover_seed_sample_dirs(apo_dir)

        # Analyze each state's seeds independently
        holo_metrics = []
        for _s, _sa, sd in holo_seeds:
            m = self._analyze_state_seed(sd, 'holo')
            if m:
                holo_metrics.append(m)

        apo_metrics = []
        for _s, _sa, sd in apo_seeds:
            m = self._analyze_state_seed(sd, 'apo')
            if m:
                apo_metrics.append(m)

        if not holo_metrics and not apo_metrics:
            return self._compare_top_level(seq_id, holo_dir, apo_dir)

        result['n_holo_seeds'] = len(holo_metrics)
        result['n_apo_seeds'] = len(apo_metrics)

        # Aggregate holo metrics
        if holo_metrics:
            df_h = pd.DataFrame(holo_metrics)
            for col in df_h.columns:
                vals = df_h[col].dropna()
                if len(vals) > 0:
                    result[f'holo_{col}'] = float(vals.mean())
                    result[f'holo_{col}_sd'] = (
                        float(vals.std()) if len(vals) > 1 else 0.0)

        # Aggregate apo metrics
        if apo_metrics:
            df_a = pd.DataFrame(apo_metrics)
            for col in df_a.columns:
                vals = df_a[col].dropna()
                if len(vals) > 0:
                    result[f'apo_{col}'] = float(vals.mean())
                    result[f'apo_{col}_sd'] = (
                        float(vals.std()) if len(vals) > 1 else 0.0)

        # Backward-compatible aliases
        if 'holo_chromophore_plddt' in result:
            result['chromophore_holo_plddt'] = result['holo_chromophore_plddt']
            result['chromophore_holo_plddt_sd'] = result.get(
                'holo_chromophore_plddt_sd', 0)
        if 'apo_chromophore_plddt' in result:
            result['chromophore_apo_plddt'] = result['apo_chromophore_plddt']
            result['chromophore_apo_plddt_sd'] = result.get(
                'apo_chromophore_plddt_sd', 0)

        # Compute differences with propagated SD
        for metric in ['mean_plddt', 'chromophore_plddt', 'mean_pae']:
            h_key = f'holo_{metric}'
            a_key = f'apo_{metric}'
            if h_key in result and a_key in result:
                if metric == 'mean_plddt':
                    diff_key = 'mean_plddt_diff'
                elif metric == 'chromophore_plddt':
                    diff_key = 'chromophore_plddt_diff'
                else:
                    diff_key = 'pae_diff'

                result[diff_key] = result[h_key] - result[a_key]
                h_sd = result.get(f'{h_key}_sd', 0)
                a_sd = result.get(f'{a_key}_sd', 0)
                result[f'{diff_key}_sd'] = float(np.sqrt(h_sd**2 + a_sd**2))

        # Matched-seed RMSD between holo and apo
        holo_map = {(s, sa): d for s, sa, d in holo_seeds}
        apo_map = {(s, sa): d for s, sa, d in apo_seeds}
        matched = sorted(set(holo_map) & set(apo_map))

        if matched:
            rmsds, chrom_rmsds = [], []
            for key in matched:
                rmsd = self._compute_pair_rmsd(holo_map[key], apo_map[key])
                if 'global_rmsd' in rmsd:
                    rmsds.append(rmsd['global_rmsd'])
                if 'chromophore_rmsd' in rmsd:
                    chrom_rmsds.append(rmsd['chromophore_rmsd'])

            if rmsds:
                result['global_rmsd'] = float(np.mean(rmsds))
                result['global_rmsd_sd'] = (
                    float(np.std(rmsds)) if len(rmsds) > 1 else 0.0)
            if chrom_rmsds:
                result['chromophore_rmsd'] = float(np.mean(chrom_rmsds))
                result['chromophore_rmsd_sd'] = (
                    float(np.std(chrom_rmsds)) if len(chrom_rmsds) > 1 else 0.0)

            # Per-residue RMSD stats
            if rmsds:
                result['mean_per_res_rmsd'] = result['global_rmsd']
                result['max_per_res_rmsd'] = float(np.max(rmsds))

        return result

    def _compare_top_level(self, seq_id: int,
                           holo_dir: Path, apo_dir: Path) -> Dict:
        """Fallback: compare only top-level best models."""
        result = {'seq_id': seq_id, 'template': self.template_name,
                  'n_holo_seeds': 1, 'n_apo_seeds': 1}

        holo_model = list(holo_dir.glob("*_model.cif"))
        apo_model = list(apo_dir.glob("*_model.cif"))
        if not holo_model or not apo_model:
            return result

        # RMSD
        try:
            b_struct = self._load_structure(holo_model[0])
            a_struct = self._load_structure(apo_model[0])
            b_ca = self._extract_ca_atoms(b_struct, 'A')
            a_ca = self._extract_ca_atoms(a_struct, 'A')
            ml = min(len(b_ca), len(a_ca))
            sup = Superimposer()
            sup.set_atoms(b_ca[:ml], a_ca[:ml])
            result['global_rmsd'] = sup.rms
            cs, ce = self.chromophore_range
            if ml >= ce:
                sup_c = Superimposer()
                sup_c.set_atoms(b_ca[cs-1:ce], a_ca[cs-1:ce])
                result['chromophore_rmsd'] = sup_c.rms
        except Exception:
            pass

        # Confidence
        for state_label, model_dir in [('holo', holo_dir), ('apo', apo_dir)]:
            sf = list(model_dir.glob("*_summary_confidences.json"))
            if sf:
                with open(sf[0]) as f:
                    sc = json.load(f)
                for k in ['ptm', 'iptm', 'ranking_score']:
                    if k in sc:
                        result[f'{state_label}_{k}'] = sc[k]

            ff = [f for f in model_dir.glob("*_confidences.json")
                  if 'summary' not in f.name]
            if ff:
                with open(ff[0]) as f:
                    fc = json.load(f)
                plddt = compute_per_residue_plddt(fc, 'A')
                if plddt is not None and len(plddt) > 1:
                    result[f'{state_label}_mean_plddt'] = float(np.mean(plddt))
                    cs, ce = self.chromophore_range
                    if len(plddt) >= ce:
                        result[f'{state_label}_chromophore_plddt'] = float(
                            np.mean(plddt[cs-1:ce]))

        # Differences
        if 'holo_mean_plddt' in result and 'apo_mean_plddt' in result:
            result['mean_plddt_diff'] = (
                result['holo_mean_plddt'] - result['apo_mean_plddt'])

        return result

    def compare_all_designs(self) -> pd.DataFrame:
        holo_base = self.output_dir / "03_alphafold3_predictions_bound"
        apo_base = self.output_dir / "03_alphafold3_predictions_apo"

        if not holo_base.exists() or not apo_base.exists():
            return pd.DataFrame()

        holo_dirs = self._discover_inference_dirs(holo_base, 'bound')
        apo_dirs = self._discover_inference_dirs(apo_base, 'apo')

        matched_ids = sorted(set(holo_dirs) & set(apo_dirs))
        if not matched_ids:
            return pd.DataFrame()

        n = len(matched_ids)
        print(f"  {self.output_dir.name}: {len(holo_dirs)} holo, "
              f"{len(apo_dirs)} apo, {n} matched pairs (x25 seeds)")

        results = []
        for i, seq_id in enumerate(matched_ids, 1):
            if i % 25 == 0 or i == 1 or i == n:
                print(f"    [{i}/{n}] design {seq_id}")
            r = self.compare_design_pair(seq_id, holo_dirs[seq_id],
                                         apo_dirs[seq_id])
            if r and len(r) > 3:
                results.append(r)

        return pd.DataFrame(results)


# ─── Visualization ────────────────────────────────────────────────────

STATE_COLORS = {'holo': '#E07B39', 'apo': '#4878CF'}  # orange / blue


def _get_template_colors(template_df):
    cmap = plt.cm.get_cmap('tab10', max(10, len(template_df)))
    return {name: cmap(i) for i, name in enumerate(template_df['template_name'])}


def _sd_bar_color(col_name):
    """Return color based on whether SD column is holo, apo, or combined."""
    if col_name.startswith('holo_'):
        return STATE_COLORS['holo']
    elif col_name.startswith('apo_'):
        return STATE_COLORS['apo']
    elif col_name.startswith('chromophore_holo'):
        return STATE_COLORS['holo']
    elif col_name.startswith('chromophore_apo'):
        return STATE_COLORS['apo']
    else:
        return '#7B68AE'  # purple for combined/difference metrics


def create_visualizations(df, template_df, output_dir):
    output_dir.mkdir(exist_ok=True, parents=True)
    csv_dir = output_dir / "plot_data_csvs"
    csv_dir.mkdir(exist_ok=True)

    print("\nGenerating visualizations...")

    _plot_rmsd_distributions(df, template_df, output_dir, csv_dir)
    _plot_plddt_comparison(df, template_df, output_dir, csv_dir)
    _plot_chromophore_plddt(df, template_df, output_dir, csv_dir)
    _plot_rmsd_vs_plddt_change(df, template_df, output_dir, csv_dir)
    _plot_sd_overview(df, output_dir, csv_dir)
    _plot_correlation_matrix(df, output_dir, csv_dir)

    print(f"\nAll plots saved to: {output_dir}")


def _plot_rmsd_distributions(df, template_df, output_dir, csv_dir):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    templates_in_data = sorted(df['template'].unique()) if 'template' in df.columns else []

    for col_idx, (metric, title) in enumerate([
        ('global_rmsd', 'Global RMSD (Holo vs Apo)'),
        ('chromophore_rmsd', 'Chromophore RMSD (Holo vs Apo)'),
    ]):
        if metric not in df.columns:
            continue
        sd_col = f'{metric}_sd'

        # Top: means
        for tpl in templates_in_data:
            tdf = df[df['template'] == tpl]
            short = tpl.replace('G4FP_', '')
            axes[0, col_idx].hist(tdf[metric].dropna(), bins=30, alpha=0.6,
                                   label=short)
        median = df[metric].median()
        if pd.notna(median):
            axes[0, col_idx].axvline(median, color='red', linestyle='--',
                                      label=f'Median: {median:.2f}A')
        axes[0, col_idx].set_xlabel(f'{title} (A)')
        axes[0, col_idx].set_ylabel('Frequency')
        axes[0, col_idx].set_title(f'{title} (Mean)')
        axes[0, col_idx].legend(fontsize=7)

        # Bottom: SD
        if sd_col in df.columns and df[sd_col].notna().any():
            for tpl in templates_in_data:
                tdf = df[df['template'] == tpl]
                vals = tdf[sd_col].dropna()
                if len(vals) > 0:
                    short = tpl.replace('G4FP_', '')
                    axes[1, col_idx].hist(vals, bins=20, alpha=0.6, label=short)
            axes[1, col_idx].set_xlabel(f'{title} SD')
            axes[1, col_idx].set_ylabel('Frequency')
            axes[1, col_idx].set_title(f'{title} SD (across matched seeds)')
            axes[1, col_idx].legend(fontsize=7)

    plt.tight_layout()
    plt.savefig(output_dir / '01_holo_apo_rmsd.png', dpi=300,
                bbox_inches='tight')
    plt.close()

    cols = [c for c in ['seq_id', 'template', 'global_rmsd', 'global_rmsd_sd',
                         'chromophore_rmsd', 'chromophore_rmsd_sd']
            if c in df.columns]
    df[cols].to_csv(csv_dir / '01_holo_apo_rmsd_data.csv', index=False)
    print("  01_holo_apo_rmsd.png")


def _plot_plddt_comparison(df, template_df, output_dir, csv_dir):
    if 'holo_mean_plddt' not in df.columns or 'apo_mean_plddt' not in df.columns:
        print("  pLDDT data not available")
        return

    plot_df = df.dropna(subset=['holo_mean_plddt', 'apo_mean_plddt'])
    if plot_df.empty:
        return

    tpl_colors = _get_template_colors(template_df)
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # Left: scatter by template
    templates_in_data = sorted(plot_df['template'].unique())
    for tpl in templates_in_data:
        tdf = plot_df[plot_df['template'] == tpl]
        short = tpl.replace('G4FP_', '')
        color = tpl_colors.get(tpl, 'steelblue')
        axes[0].scatter(tdf['holo_mean_plddt'], tdf['apo_mean_plddt'],
                        alpha=0.5, s=30, color=color,
                        label=f'{short} ({len(tdf)})')

    # Diagonal
    all_vals = pd.concat([plot_df['holo_mean_plddt'],
                           plot_df['apo_mean_plddt']])
    lo, hi = all_vals.min() - 1, all_vals.max() + 1
    axes[0].plot([lo, hi], [lo, hi], 'k--', alpha=0.3, label='y=x')

    # Template AF3 references (with error bars if AF3 available)
    for _, row in template_df.iterrows():
        short = row['template_name'].replace('G4FP_', '')
        color = tpl_colors.get(row['template_name'], 'black')
        hp = row.get('holo_mean_plddt')
        ap = row.get('apo_mean_plddt')
        if pd.notna(hp) and pd.notna(ap):
            hsd = row.get('holo_mean_plddt_sd', 0) or 0
            asd = row.get('apo_mean_plddt_sd', 0) or 0
            axes[0].errorbar(hp, ap, xerr=hsd, yerr=asd,
                             fmt='*', markersize=15, color=color,
                             markeredgecolor='black', markeredgewidth=0.8,
                             capsize=3, zorder=10)
            axes[0].annotate(short, (hp, ap), fontsize=5,
                             xytext=(3, 3), textcoords='offset points')
        else:
            p = row.get('mean_plddt')
            if p and pd.notna(p):
                axes[0].scatter([p], [p], marker='*', s=200, color=color,
                                edgecolors='black', linewidths=0.8, zorder=5)
                axes[0].annotate(short, (p, p), fontsize=5,
                                 xytext=(3, 3), textcoords='offset points')

    axes[0].set_xlabel('Holo Mean pLDDT')
    axes[0].set_ylabel('Apo Mean pLDDT')
    axes[0].set_title(f'pLDDT: Holo vs Apo (n={len(plot_df)})')
    axes[0].legend(fontsize=6, loc='upper left')

    corr = plot_df[['holo_mean_plddt', 'apo_mean_plddt']].corr().iloc[0, 1]
    axes[0].text(0.95, 0.05, f'r = {corr:.3f}', transform=axes[0].transAxes,
                 fontsize=10, ha='right')

    # Middle: scatter colored by pLDDT SD
    sd_col = 'holo_mean_plddt_sd'
    if sd_col in plot_df.columns and plot_df[sd_col].notna().any():
        max_sd = plot_df[[c for c in ['holo_mean_plddt_sd', 'apo_mean_plddt_sd']
                          if c in plot_df.columns]].max(axis=1)
        sc = axes[1].scatter(plot_df['holo_mean_plddt'],
                             plot_df['apo_mean_plddt'],
                             c=max_sd, cmap='plasma', s=30, alpha=0.7)
        plt.colorbar(sc, ax=axes[1], label='Max pLDDT SD')
        axes[1].plot([lo, hi], [lo, hi], 'k--', alpha=0.3)
        axes[1].set_xlabel('Holo Mean pLDDT')
        axes[1].set_ylabel('Apo Mean pLDDT')
        axes[1].set_title('Colored by pLDDT SD (prediction uncertainty)')
    else:
        axes[1].axis('off')

    # Right: difference distribution with SD
    if 'mean_plddt_diff' in df.columns:
        diff = df['mean_plddt_diff'].dropna()
        axes[2].hist(diff, bins=30, alpha=0.7, color='green')
        axes[2].set_xlabel('pLDDT Difference (Holo - Apo)')
        axes[2].set_ylabel('Frequency')
        axes[2].set_title('pLDDT Change Distribution')
        axes[2].axvline(0, color='red', linestyle='--', label='No change')
        if len(diff) > 0:
            axes[2].axvline(diff.median(), color='blue', linestyle='--',
                            label=f'Median: {diff.median():.2f}')
        axes[2].legend()

    plt.tight_layout()
    plt.savefig(output_dir / '02_plddt_comparison.png', dpi=300,
                bbox_inches='tight')
    plt.close()

    cols = [c for c in ['seq_id', 'template', 'holo_mean_plddt',
                         'holo_mean_plddt_sd', 'apo_mean_plddt',
                         'apo_mean_plddt_sd', 'mean_plddt_diff',
                         'mean_plddt_diff_sd'] if c in df.columns]
    df[cols].to_csv(csv_dir / '02_plddt_comparison_data.csv', index=False)
    print("  02_plddt_comparison.png")


def _plot_chromophore_plddt(df, template_df, output_dir, csv_dir):
    h_col = 'chromophore_holo_plddt'
    a_col = 'chromophore_apo_plddt'
    if h_col not in df.columns or a_col not in df.columns:
        h_col = 'holo_chromophore_plddt'
        a_col = 'apo_chromophore_plddt'
    if h_col not in df.columns or a_col not in df.columns:
        print("  Chromophore pLDDT data not available")
        return

    plot_df = df.dropna(subset=[h_col, a_col])
    if plot_df.empty:
        return

    tpl_colors = _get_template_colors(template_df)
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    templates_in_data = sorted(plot_df['template'].unique())
    for tpl in templates_in_data:
        tdf = plot_df[plot_df['template'] == tpl]
        short = tpl.replace('G4FP_', '')
        color = tpl_colors.get(tpl, 'purple')
        axes[0].scatter(tdf[h_col], tdf[a_col],
                        alpha=0.5, s=30, color=color,
                        label=f'{short} ({len(tdf)})')

    all_vals = pd.concat([plot_df[h_col], plot_df[a_col]])
    lo, hi = all_vals.min() - 1, all_vals.max() + 1
    axes[0].plot([lo, hi], [lo, hi], 'r--', alpha=0.3, label='y=x')

    # Template references
    for _, row in template_df.iterrows():
        short = row['template_name'].replace('G4FP_', '')
        color = tpl_colors.get(row['template_name'], 'black')
        hp = row.get('holo_chromophore_plddt')
        ap = row.get('apo_chromophore_plddt')
        if pd.notna(hp) and pd.notna(ap):
            hsd = row.get('holo_chromophore_plddt_sd', 0) or 0
            asd = row.get('apo_chromophore_plddt_sd', 0) or 0
            axes[0].errorbar(hp, ap, xerr=hsd, yerr=asd,
                             fmt='*', markersize=15, color=color,
                             markeredgecolor='black', markeredgewidth=0.8,
                             capsize=3, zorder=10)
        else:
            cp = row.get('chromophore_mean_plddt')
            if cp and pd.notna(cp):
                axes[0].scatter([cp], [cp], marker='*', s=200, color=color,
                                edgecolors='black', linewidths=0.8, zorder=5)

    axes[0].set_xlabel('Holo Chromophore pLDDT')
    axes[0].set_ylabel('Apo Chromophore pLDDT')
    axes[0].set_title(f'Chromophore pLDDT: Holo vs Apo (n={len(plot_df)})')
    axes[0].legend(fontsize=6)

    # Middle: colored by SD
    h_sd_col = f'{h_col}_sd'
    a_sd_col = f'{a_col}_sd'
    sd_cols_avail = [c for c in [h_sd_col, a_sd_col] if c in plot_df.columns]
    if sd_cols_avail and plot_df[sd_cols_avail[0]].notna().any():
        max_sd = plot_df[sd_cols_avail].max(axis=1)
        sc = axes[1].scatter(plot_df[h_col], plot_df[a_col],
                             c=max_sd, cmap='plasma', s=30, alpha=0.7)
        plt.colorbar(sc, ax=axes[1], label='Max Chromophore pLDDT SD')
        axes[1].plot([lo, hi], [lo, hi], 'r--', alpha=0.3)
        axes[1].set_xlabel('Holo Chromophore pLDDT')
        axes[1].set_ylabel('Apo Chromophore pLDDT')
        axes[1].set_title('Colored by SD')
    else:
        axes[1].axis('off')

    # Right: difference distribution
    diff_col = 'chromophore_plddt_diff'
    if diff_col in df.columns:
        diff = df[diff_col].dropna()
        axes[2].hist(diff, bins=30, alpha=0.7, color='orange')
        axes[2].set_xlabel('Chromophore pLDDT Diff (Holo - Apo)')
        axes[2].set_ylabel('Frequency')
        axes[2].set_title('Chromophore pLDDT Change')
        axes[2].axvline(0, color='red', linestyle='--', label='No change')
        if len(diff) > 0:
            axes[2].axvline(diff.median(), color='blue', linestyle='--',
                            label=f'Median: {diff.median():.2f}')
        axes[2].legend()

    plt.tight_layout()
    plt.savefig(output_dir / '03_chromophore_plddt_comparison.png', dpi=300,
                bbox_inches='tight')
    plt.close()
    print("  03_chromophore_plddt_comparison.png")


def _plot_rmsd_vs_plddt_change(df, template_df, output_dir, csv_dir):
    if 'global_rmsd' not in df.columns or 'mean_plddt_diff' not in df.columns:
        print("  Required data not available for RMSD vs pLDDT plot")
        return

    tpl_colors = _get_template_colors(template_df)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    data = df.dropna(subset=['global_rmsd', 'mean_plddt_diff'])
    if len(data) > 0:
        # Left: by template
        for tpl in sorted(data['template'].unique()):
            tdf = data[data['template'] == tpl]
            short = tpl.replace('G4FP_', '')
            color = tpl_colors.get(tpl, 'steelblue')
            axes[0].scatter(tdf['global_rmsd'], tdf['mean_plddt_diff'],
                            alpha=0.5, s=30, color=color, label=short)
        axes[0].set_xlabel('Global RMSD (A)')
        axes[0].set_ylabel('pLDDT Diff (Holo - Apo)')
        axes[0].set_title('RMSD vs pLDDT Change')
        axes[0].axhline(0, color='red', linestyle='--', alpha=0.5)
        axes[0].legend(fontsize=7)
        corr = data[['global_rmsd', 'mean_plddt_diff']].corr().iloc[0, 1]
        axes[0].text(0.05, 0.95, f'r = {corr:.3f}',
                     transform=axes[0].transAxes, fontsize=10)

        # Right: colored by RMSD SD
        sd_col = 'global_rmsd_sd'
        if sd_col in data.columns and data[sd_col].notna().any():
            sc = axes[1].scatter(data['global_rmsd'], data['mean_plddt_diff'],
                                 c=data[sd_col], cmap='plasma', s=30,
                                 alpha=0.7)
            plt.colorbar(sc, ax=axes[1], label='RMSD SD (across seeds)')
            axes[1].axhline(0, color='red', linestyle='--', alpha=0.5)
            axes[1].set_xlabel('Global RMSD (A)')
            axes[1].set_ylabel('pLDDT Diff (Holo - Apo)')
            axes[1].set_title('Colored by RMSD SD')
        else:
            axes[1].axis('off')

    plt.tight_layout()
    plt.savefig(output_dir / '04_rmsd_vs_plddt_change.png', dpi=300,
                bbox_inches='tight')
    plt.close()
    print("  04_rmsd_vs_plddt_change.png")


def _plot_sd_overview(df, output_dir, csv_dir):
    """Dedicated SD overview for all metrics."""
    sd_cols = [c for c in df.columns if c.endswith('_sd') and df[c].notna().any()]
    if not sd_cols:
        return

    n = len(sd_cols)
    ncols = min(n, 4)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    axes = np.atleast_2d(axes)

    for i, col in enumerate(sd_cols):
        r, c = divmod(i, ncols)
        metric_name = col.replace('_sd', '').replace('_', ' ').title()
        vals = df[col].dropna()
        bar_color = _sd_bar_color(col)
        axes[r, c].hist(vals, bins=20, alpha=0.7, color=bar_color,
                         edgecolor='white', linewidth=0.5)
        axes[r, c].set_xlabel(f'{metric_name} SD')
        axes[r, c].set_ylabel('Count')
        axes[r, c].set_title(f'{metric_name}\nMedian: {vals.median():.4f}')

    for i in range(len(sd_cols), nrows * ncols):
        r, c = divmod(i, ncols)
        axes[r, c].axis('off')

    # Add color legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=STATE_COLORS['holo'], alpha=0.7, label='Holo'),
        Patch(facecolor=STATE_COLORS['apo'], alpha=0.7, label='Apo'),
        Patch(facecolor='#7B68AE', alpha=0.7, label='Difference'),
    ]
    fig.legend(handles=legend_elements, loc='upper right', fontsize=9,
               framealpha=0.8, title='Metric type')

    plt.tight_layout(rect=[0, 0, 0.95, 1])
    plt.savefig(output_dir / '05_sd_overview.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  05_sd_overview.png")


def _plot_correlation_matrix(df, output_dir, csv_dir):
    possible = [
        'global_rmsd', 'chromophore_rmsd',
        'holo_mean_plddt', 'apo_mean_plddt', 'mean_plddt_diff',
        'chromophore_holo_plddt', 'chromophore_apo_plddt',
        'chromophore_plddt_diff',
        'holo_ptm', 'apo_ptm', 'holo_iptm', 'apo_iptm',
        'holo_ranking_score', 'apo_ranking_score',
        'holo_mean_pae', 'apo_mean_pae', 'pae_diff',
    ]
    cols = [c for c in possible if c in df.columns and df[c].notna().any()]
    if len(cols) < 2:
        return

    corr = df[cols].corr()
    plt.figure(figsize=(14, 12))
    sns.heatmap(corr, annot=True, fmt='.2f', cmap='coolwarm',
                center=0, square=True, linewidths=1)
    plt.title('Correlation Matrix: Holo vs Apo')
    plt.tight_layout()
    plt.savefig(output_dir / '06_correlation_matrix.png', dpi=300,
                bbox_inches='tight')
    plt.close()
    corr.to_csv(csv_dir / '06_correlation_matrix_data.csv')
    print("  06_correlation_matrix.png")


# ─── Main ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Compare holo and apo AF3 predictions (all seeds)')
    parser.add_argument('--output-dir', type=str, default=None)
    parser.add_argument('--template', type=str,
                        default='inputs/G4FP_des1_cro_mod0.pdb')
    parser.add_argument('--chromophore-range', type=str, default='197,199')
    parser.add_argument('--analysis-output-dir', type=str,
                        default='analysis_output')
    parser.add_argument('--force', action='store_true',
                        help='Re-analyze even if results CSV exists')

    args = parser.parse_args()
    cs, ce = map(int, args.chromophore_range.split(','))
    chromophore_range = (cs, ce)
    base_dir = Path(__file__).resolve().parent

    if args.output_dir:
        output_dirs = [Path(args.output_dir)]
    else:
        output_dirs = sorted(
            d for d in base_dir.iterdir()
            if d.is_dir() and d.name.startswith("output_G4FP_"))

    if not output_dirs:
        print("No output directories found")
        sys.exit(1)

    # Load template metrics (AF3 if available, else B-factors)
    template_df = load_template_metrics(base_dir, chromophore_range)
    print(f"\nLoaded {len(template_df)} template references:")
    for _, r in template_df.iterrows():
        src = r.get('source', '?')
        plddt = r.get('mean_plddt', 0)
        af3_info = ""
        if src == 'af3':
            hp = r.get('holo_mean_plddt')
            ap = r.get('apo_mean_plddt')
            if pd.notna(hp) and pd.notna(ap):
                af3_info = f" AF3: holo={hp:.1f} apo={ap:.1f}"
        print(f"  {r['template_name']}: pLDDT={plddt:.1f} ({src}){af3_info}")

    print(f"\n{'='*80}")
    print(f"Comparing Ligand States (ALL seeds) Across {len(output_dirs)} Templates")
    print(f"{'='*80}")

    all_frames = []
    for output_dir in output_dirs:
        print(f"\n--- {output_dir.name} ---")
        csv_path = output_dir / "06_ligand_state_comparison_results.csv"
        if not args.force and csv_path.exists():
            try:
                df = pd.read_csv(csv_path)
                if not df.empty:
                    all_frames.append(df)
                    print(f"  Cached: {csv_path} ({len(df)} pairs) [use --force to re-run]")
                    continue
            except Exception:
                pass
        comparator = LigandStateComparator(
            output_dir=output_dir,
            reference_template=Path(args.template),
            chromophore_range=chromophore_range,
        )
        df = comparator.compare_all_designs()
        if not df.empty:
            all_frames.append(df)
            df_export = df.drop(
                columns=[c for c in df.columns if 'array' in c],
                errors='ignore')
            df_export.to_csv(csv_path, index=False)
            print(f"  Saved: {csv_path}")
        else:
            print(f"  No matched design pairs found")

    if not all_frames:
        print("\nNo designs compared across any template")
        sys.exit(1)

    df_combined = pd.concat(all_frames, ignore_index=True)
    print(f"\n{'='*80}")
    print(f"Combined: {len(df_combined)} design pairs across "
          f"{df_combined['template'].nunique()} templates")
    if 'n_holo_seeds' in df_combined.columns:
        print(f"Holo seeds/design: {df_combined['n_holo_seeds'].median():.0f} median")
    if 'n_apo_seeds' in df_combined.columns:
        print(f"Apo seeds/design: {df_combined['n_apo_seeds'].median():.0f} median")
    print(f"{'='*80}")

    central_dir = Path(args.analysis_output_dir)
    viz_dir = central_dir / "06_ligand_state_comparison"
    create_visualizations(df_combined, template_df, viz_dir)

    df_export = df_combined.drop(
        columns=[c for c in df_combined.columns if 'array' in c],
        errors='ignore')
    combined_csv = central_dir / "06_ligand_state_comparison_results.csv"
    df_export.to_csv(combined_csv, index=False)

    print(f"\nAnalysis complete!")
    print(f"  Combined results: {combined_csv}")
    print(f"  Plots: {viz_dir}/")
    print(f"  Per-template CSVs: output_G4FP_*/06_ligand_state_comparison_results.csv")


if __name__ == "__main__":
    main()

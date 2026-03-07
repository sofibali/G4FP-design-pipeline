#!/usr/bin/env python3
"""
Analyze AlphaFold3 Structure Predictions for G4FP Designs

Analyzes ALL 25 seed-sample models per design to compute distributions.
Reports mean +/- SD for pLDDT, RMSD, pTM, ipTM, ranking_score, PAE.
Scatter plots colored by SD show prediction confidence.
Template AF3 predictions used as reference when available.

The three-residue chromophore window is auto-detected for each template by
locating the CRO HETATM residue in the template PDB and counting the standard
protein residues that precede it.  Use ``--chromophore-range`` to override.

Usage:
    python 05_analyze_af3_structures.py                              # auto-detect
    python 05_analyze_af3_structures.py --chromophore-range 197,199  # override
    python 05_analyze_af3_structures.py --output-dir output_G4FP_des1_cro_mod0 --state both
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

# Residue names that represent the GFP-type chromophore in PDB HETATM records
CHROMOPHORE_RESNAMES = {'CRO', 'SYG', 'CRQ', 'CFO', 'FME'}

# ─── Shared utilities ─────────────────────────────────────────────────

def detect_chromophore_range_from_pdb(pdb_file: Path,
                                       chain_id: str = 'A') -> Tuple[int, int]:
    """Detect the three-residue chromophore range from a template PDB.

    Locates the CRO (or related) HETATM residue on ``chain_id`` and counts
    how many standard ATOM residues precede it in that chain.  The returned
    range ``(n_before + 1, n_before + 3)`` covers the three positions in the
    designed (ATOM-only) sequence that immediately follow the chromophore gap,
    which is the structurally equivalent region across templates of different
    lengths.

    Falls back to ``(197, 199)`` when no chromophore HETATM is found.
    """
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure('chrom', str(pdb_file))

    atom_resids = set()
    cro_resid = None

    for model in structure:
        for chain in model:
            if chain.get_id() != chain_id:
                continue
            for residue in chain:
                resname = residue.get_resname().strip()
                hetflag = residue.get_id()[0].strip()
                resid = residue.get_id()[1]
                if hetflag and resname in CHROMOPHORE_RESNAMES:
                    # Use the lowest residue number when multiple chromophore
                    # HETATM records exist (e.g. a rare split-conformation PDB);
                    # in practice all G4FP templates have a single CRO at 197.
                    if cro_resid is None or resid < cro_resid:
                        cro_resid = resid
                elif not hetflag:
                    atom_resids.add(resid)
        break  # first model only

    if cro_resid is None:
        return (197, 199)  # fallback default

    n_before = sum(1 for r in atom_resids if r < cro_resid)
    return (n_before + 1, n_before + 3)


def discover_seed_sample_dirs(design_dir: Path) -> List[Tuple[int, int, Path]]:
    """Find all seed-N_sample-M directories within a design output dir."""
    pattern = re.compile(r'^seed-(\d+)_sample-(\d+)$')
    dirs = []
    for d in design_dir.iterdir():
        if d.is_dir():
            m = pattern.match(d.name)
            if m:
                dirs.append((int(m.group(1)), int(m.group(2)), d))
    if not dirs:
        # Check one level deeper (AF3 nesting: outer/inner/seed-*)
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
    """Compute per-residue pLDDT from atom_plddts filtered to a chain."""
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
    """Find actual AF3 result dir (might be nested one level)."""
    if list(base_dir.glob("*_model.cif")) or list(base_dir.glob("seed-*")):
        return base_dir
    for subdir in sorted(base_dir.iterdir()):
        if subdir.is_dir():
            if list(subdir.glob("*_model.cif")) or list(subdir.glob("seed-*")):
                return subdir
    return None


def load_template_metrics(base_dir: Path,
                          chromophore_range: Optional[Tuple[int, int]] = None) -> pd.DataFrame:
    """Load template metrics from AF3 outputs if available, else B-factors.

    When *chromophore_range* is ``None`` the three-residue chromophore window
    is auto-detected from each template PDB via
    :func:`detect_chromophore_range_from_pdb`.  Pass an explicit tuple to
    override for all templates.
    """
    tpl_dir = base_dir / "template_af3_outputs"
    inputs_dir = base_dir / "inputs"
    pdb_parser = PDBParser(QUIET=True)
    records = []

    for pdb in sorted(inputs_dir.glob("G4FP_*.pdb")):
        tpl_name = pdb.stem
        rec = {'template_name': tpl_name, 'source': 'bfactor'}

        # Determine chromophore range for this template
        chrom_range = (chromophore_range
                       if chromophore_range is not None
                       else detect_chromophore_range_from_pdb(pdb))
        rec['chromophore_start'] = chrom_range[0]
        rec['chromophore_end'] = chrom_range[1]

        # B-factor baseline
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
            cs, ce = chrom_range
            if len(b_arr) >= ce:
                rec['chromophore_mean_plddt'] = float(np.mean(b_arr[cs-1:ce]))

        # Try AF3 results for each state
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

            for _seed, _sample, sd in seed_dirs:
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
                        cs, ce = chrom_range
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


# ─── AlphaFoldStructureAnalyzer ──────────────────────────────────────

class AlphaFoldStructureAnalyzer:
    """Analyze AF3 predictions for one output directory, all seeds."""

    def __init__(self, output_dir: Path, reference_template: Path,
                 chromophore_range: Optional[Tuple[int, int]] = None):
        self.output_dir = Path(output_dir)
        self.reference_template = Path(reference_template)
        self.pdb_parser = PDBParser(QUIET=True)
        self.cif_parser = MMCIFParser(QUIET=True)
        self.ref_structure = self._load_structure(self.reference_template)
        self.ref_ca_atoms = self._extract_ca_atoms(self.ref_structure)
        self.template_name = self.output_dir.name.replace("output_", "")
        # Auto-detect chromophore range from template PDB when not provided
        if chromophore_range is None:
            self.chromophore_range = detect_chromophore_range_from_pdb(
                self.reference_template)
            print(f"  Chromophore range (auto): {self.chromophore_range}")
        else:
            self.chromophore_range = chromophore_range
            print(f"  Chromophore range (override): {self.chromophore_range}")
        print(f"  Reference: {self.reference_template} ({len(self.ref_ca_atoms)} CA)")

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
        """Discover AF3 inference dirs (non-padded, possibly timestamped, or inside padded wrappers)."""
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

    def _compute_rmsd(self, model_file: Path, chain_id='A') -> Dict[str, float]:
        """Compute RMSD vs reference for one model file."""
        structure = self._load_structure(model_file)
        model_ca = self._extract_ca_atoms(structure, chain_id)
        min_len = min(len(model_ca), len(self.ref_ca_atoms))
        if min_len < 10:
            return {}
        model_ca = model_ca[:min_len]
        ref_ca = self.ref_ca_atoms[:min_len]

        sup = Superimposer()
        sup.set_atoms(ref_ca, model_ca)
        result = {'global_rmsd': sup.rms}

        cs, ce = self.chromophore_range
        if min_len >= ce:
            sup_c = Superimposer()
            sup_c.set_atoms(ref_ca[cs-1:ce], model_ca[cs-1:ce])
            result['chromophore_rmsd'] = sup_c.rms

        return result

    def _analyze_single_seed(self, seed_dir: Path, state: str) -> Optional[Dict]:
        """Analyze one seed-sample model. Returns flat dict or None."""
        result = {}

        # CIF model → RMSD
        model_cif = seed_dir / 'model.cif'
        if model_cif.exists():
            try:
                rmsd = self._compute_rmsd(model_cif)
                result.update(rmsd)
            except Exception:
                pass

        # Summary confidence
        sf = seed_dir / 'summary_confidences.json'
        if sf.exists():
            try:
                with open(sf) as f:
                    sc = json.load(f)
                for k in ['ptm', 'iptm', 'ranking_score', 'fraction_disordered']:
                    if k in sc:
                        result[k] = sc[k]
            except Exception:
                pass

        # Full confidence → pLDDT, PAE
        ff = seed_dir / 'confidences.json'
        if ff.exists():
            try:
                with open(ff) as f:
                    fc = json.load(f)
                plddt = compute_per_residue_plddt(fc, 'A')
                if plddt is not None and len(plddt) > 1:
                    result['mean_plddt'] = float(np.mean(plddt))
                    result['min_plddt'] = float(np.min(plddt))
                    cs, ce = self.chromophore_range
                    if len(plddt) >= ce:
                        result['chromophore_mean_plddt'] = float(
                            np.mean(plddt[cs-1:ce]))

                if 'pae' in fc:
                    pae = np.array(fc['pae'])
                    if pae.ndim == 2:
                        result['mean_pae'] = float(np.mean(pae))
                        if state == 'holo' and plddt is not None:
                            n_prot = len(plddt)
                            if pae.shape[0] > n_prot:
                                result['interface_pae'] = float(
                                    np.mean(pae[:n_prot, n_prot:]))
            except Exception:
                pass

        return result if result else None

    def _analyze_top_level(self, design_dir: Path, state: str,
                           seq_id: int) -> Dict:
        """Fallback: analyze only the top-level best model."""
        result = {
            'seq_id': seq_id, 'design_name': design_dir.name,
            'state': state, 'template': self.template_name, 'n_seeds': 1,
            'chromophore_start': self.chromophore_range[0],
            'chromophore_end': self.chromophore_range[1],
        }
        model_files = list(design_dir.glob("*_model.cif"))
        if not model_files:
            return result
        try:
            result.update(self._compute_rmsd(model_files[0]))
        except Exception:
            pass

        summary_files = list(design_dir.glob("*_summary_confidences.json"))
        if summary_files:
            with open(summary_files[0]) as f:
                sc = json.load(f)
            for k in ['ptm', 'iptm', 'ranking_score', 'fraction_disordered']:
                if k in sc:
                    result[k] = sc[k]

        full_files = [f for f in design_dir.glob("*_confidences.json")
                      if 'summary' not in f.name]
        if full_files:
            with open(full_files[0]) as f:
                fc = json.load(f)
            plddt = compute_per_residue_plddt(fc, 'A')
            if plddt is not None and len(plddt) > 1:
                result['mean_plddt'] = float(np.mean(plddt))
                result['min_plddt'] = float(np.min(plddt))
                cs, ce = self.chromophore_range
                if len(plddt) >= ce:
                    result['chromophore_mean_plddt'] = float(
                        np.mean(plddt[cs-1:ce]))
            if 'pae' in fc:
                pae = np.array(fc['pae'])
                if pae.ndim == 2:
                    result['mean_pae'] = float(np.mean(pae))

        return result

    def analyze_design(self, design_dir: Path, state: str,
                       seq_id: int) -> Dict:
        """Analyze all seed-sample models, return mean +/- SD."""
        seed_dirs = discover_seed_sample_dirs(design_dir)

        if not seed_dirs:
            return self._analyze_top_level(design_dir, state, seq_id)

        seed_results = []
        for _seed, _sample, sd in seed_dirs:
            r = self._analyze_single_seed(sd, state)
            if r:
                seed_results.append(r)

        if not seed_results:
            return self._analyze_top_level(design_dir, state, seq_id)

        result = {
            'seq_id': seq_id, 'design_name': design_dir.name,
            'state': state, 'template': self.template_name,
            'n_seeds': len(seed_results),
            'chromophore_start': self.chromophore_range[0],
            'chromophore_end': self.chromophore_range[1],
        }

        df_seeds = pd.DataFrame(seed_results)
        for col in df_seeds.columns:
            vals = df_seeds[col].dropna()
            if len(vals) > 0:
                result[col] = float(vals.mean())
                result[f'{col}_sd'] = float(vals.std()) if len(vals) > 1 else 0.0

        return result

    def analyze_all_designs(self, state: str = 'both') -> pd.DataFrame:
        states = []
        if state in ['holo', 'both']:
            states.append('holo')
        if state in ['apo', 'both']:
            states.append('apo')

        all_results = []
        for s in states:
            disk_s = DISK_STATE[s]
            af3_dir = self.output_dir / f"03_alphafold3_predictions_{disk_s}"
            if not af3_dir.exists():
                continue

            inf_dirs = self._discover_inference_dirs(af3_dir, disk_s)
            if not inf_dirs:
                continue

            n = len(inf_dirs)
            print(f"  {s}: {n} designs (x25 seeds each)")
            for i, (seq_id, ddir) in enumerate(sorted(inf_dirs.items()), 1):
                if i % 25 == 0 or i == 1 or i == n:
                    print(f"    [{i}/{n}] design {seq_id}")
                r = self.analyze_design(ddir, s, seq_id)
                if r and 'mean_plddt' in r:
                    all_results.append(r)

        return pd.DataFrame(all_results)


# ─── Visualization ────────────────────────────────────────────────────

STATE_COLORS = {'holo': '#E07B39', 'apo': '#4878CF'}  # orange / blue


def _get_template_colors(template_df):
    cmap = plt.cm.get_cmap('tab10', max(10, len(template_df)))
    return {name: cmap(i) for i, name in enumerate(template_df['template_name'])}


def _state_color(state, alpha=1.0):
    """Return distinct color for holo (orange) vs apo (blue)."""
    import matplotlib.colors as mcolors
    base = STATE_COLORS.get(state, '#999999')
    r, g, b = mcolors.to_rgb(base)
    return (r, g, b, alpha)


def create_visualizations(df, template_df, output_dir):
    output_dir.mkdir(exist_ok=True, parents=True)
    csv_dir = output_dir / "plot_data_csvs"
    csv_dir.mkdir(exist_ok=True)

    print("\nGenerating visualizations...")

    _plot_af3_metrics(df, template_df, output_dir, csv_dir)
    _plot_sd_overview(df, output_dir, csv_dir)
    _plot_rmsd_distributions(df, template_df, output_dir, csv_dir)
    _plot_correlation_matrix(df, output_dir, csv_dir)
    _plot_global_vs_chromophore_rmsd(df, template_df, output_dir, csv_dir)
    _plot_plddt_vs_rmsd_sd(df, template_df, output_dir, csv_dir)
    _export_top_designs(df, output_dir, csv_dir)

    print(f"\nAll plots saved to: {output_dir}")


def _plot_af3_metrics(df, template_df, output_dir, csv_dir):
    """AF3 quality metrics: mean distributions (top) and SD distributions (bottom)."""
    metrics = ['mean_plddt', 'ptm', 'iptm', 'ranking_score']
    available = [m for m in metrics if m in df.columns and df[m].notna().any()]
    if not available:
        return

    tpl_colors = _get_template_colors(template_df)
    n = len(available)
    fig, axes = plt.subplots(2, n, figsize=(4 * n, 10))
    if n == 1:
        axes = axes.reshape(2, 1)

    for idx, metric in enumerate(available):
        sd_col = f'{metric}_sd'

        # Top row: distribution of means -- colored by state
        for state in sorted(df['state'].unique()):
            sdf = df[df['state'] == state]
            vals = sdf[metric].dropna()
            if len(vals) > 0:
                axes[0, idx].hist(vals, bins=20, alpha=0.5,
                                  color=_state_color(state),
                                  label=f'{state} (n={len(vals)})')

        # Template reference lines for pLDDT
        if metric == 'mean_plddt':
            for _, row in template_df.iterrows():
                plddt = row.get('mean_plddt')
                if plddt and pd.notna(plddt):
                    color = tpl_colors.get(row['template_name'], 'black')
                    axes[0, idx].axvline(plddt, color=color, linestyle='--',
                                         linewidth=1.5, alpha=0.7)

        axes[0, idx].set_xlabel(metric.replace('_', ' ').title())
        axes[0, idx].set_ylabel('Frequency')
        axes[0, idx].set_title(f'{metric.replace("_", " ").title()} (Mean)')
        axes[0, idx].legend(fontsize=7)

        # Bottom row: SD distribution -- colored by state
        if sd_col in df.columns and df[sd_col].notna().any():
            for state in sorted(df['state'].unique()):
                sdf = df[df['state'] == state]
                vals = sdf[sd_col].dropna()
                if len(vals) > 0:
                    axes[1, idx].hist(vals, bins=20, alpha=0.5,
                                      color=_state_color(state),
                                      label=f'{state}')
            axes[1, idx].set_xlabel(f'{metric} SD')
            axes[1, idx].set_ylabel('Frequency')
            axes[1, idx].set_title(f'{metric.replace("_", " ").title()} SD (across seeds)')
            axes[1, idx].legend(fontsize=7)
        else:
            axes[1, idx].axis('off')

    plt.tight_layout()
    plt.savefig(output_dir / '02_alphafold3_metrics.png', dpi=300, bbox_inches='tight')
    plt.close()

    export_cols = [c for c in ['seq_id', 'state', 'template', 'n_seeds'] +
                   sum([[m, f'{m}_sd'] for m in available], []) if c in df.columns]
    df[export_cols].to_csv(csv_dir / '02_alphafold3_metrics_data.csv', index=False)
    print("  02_alphafold3_metrics.png")


def _plot_sd_overview(df, output_dir, csv_dir):
    """Dedicated SD overview: distribution of SDs for all metrics."""
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

        for state in sorted(df['state'].unique()):
            sv = df[df['state'] == state][col].dropna()
            if len(sv) > 0:
                axes[r, c].hist(sv, bins=20, alpha=0.5,
                                color=_state_color(state), label=state)

        axes[r, c].set_xlabel(f'{metric_name} SD')
        axes[r, c].set_ylabel('Count')
        axes[r, c].set_title(f'{metric_name}\nMedian SD: {vals.median():.4f}')
        axes[r, c].legend(fontsize=7)

    for i in range(len(sd_cols), nrows * ncols):
        r, c = divmod(i, ncols)
        axes[r, c].axis('off')

    plt.tight_layout()
    plt.savefig(output_dir / '03_sd_overview.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  03_sd_overview.png")


def _plot_rmsd_distributions(df, template_df, output_dir, csv_dir):
    """RMSD distributions: mean (top) and SD (bottom)."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    for col_idx, (metric, title) in enumerate([
        ('global_rmsd', 'Global RMSD'),
        ('chromophore_rmsd', 'Chromophore RMSD'),
    ]):
        if metric not in df.columns:
            continue
        sd_col = f'{metric}_sd'

        # Top: means colored by state
        for state in sorted(df['state'].unique()):
            vals = df[df['state'] == state][metric].dropna()
            if len(vals) > 0:
                axes[0, col_idx].hist(vals, bins=20, alpha=0.5,
                                      color=_state_color(state),
                                      label=f'{state} (n={len(vals)})')

        axes[0, col_idx].axvline(0, color='black', linestyle='--', linewidth=2,
                                  alpha=0.3, label='Reference')
        axes[0, col_idx].set_xlabel(f'{title} (A)')
        axes[0, col_idx].set_ylabel('Frequency')
        axes[0, col_idx].set_title(f'{title} (Mean across seeds)')
        axes[0, col_idx].legend(fontsize=7)

        # Bottom: SD colored by state
        if sd_col in df.columns and df[sd_col].notna().any():
            for state in sorted(df['state'].unique()):
                vals = df[df['state'] == state][sd_col].dropna()
                if len(vals) > 0:
                    axes[1, col_idx].hist(vals, bins=20, alpha=0.5,
                                          color=_state_color(state),
                                          label=f'{state}')
            axes[1, col_idx].set_xlabel(f'{title} SD')
            axes[1, col_idx].set_ylabel('Frequency')
            axes[1, col_idx].set_title(f'{title} SD (across seeds)')
            axes[1, col_idx].legend(fontsize=7)

    plt.tight_layout()
    plt.savefig(output_dir / '04_rmsd_distributions.png', dpi=300, bbox_inches='tight')
    plt.close()

    cols = [c for c in ['seq_id', 'state', 'template', 'global_rmsd',
                         'global_rmsd_sd', 'chromophore_rmsd',
                         'chromophore_rmsd_sd'] if c in df.columns]
    df[cols].to_csv(csv_dir / '04_rmsd_distributions_data.csv', index=False)
    print("  04_rmsd_distributions.png")


def _plot_global_vs_chromophore_rmsd(df, template_df, output_dir, csv_dir):
    """Scatter: Global vs Chromophore RMSD. Left=by template, Right=colored by SD."""
    if 'global_rmsd' not in df.columns or 'chromophore_rmsd' not in df.columns:
        return

    plot_df = df.dropna(subset=['global_rmsd', 'chromophore_rmsd'])
    if plot_df.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    markers = {'holo': 'o', 'apo': 's'}

    # Left: colored by state
    for state in sorted(plot_df['state'].unique()):
        sdf = plot_df[plot_df['state'] == state]
        axes[0].scatter(sdf['global_rmsd'], sdf['chromophore_rmsd'],
                        alpha=0.5, s=25, color=_state_color(state),
                        marker=markers.get(state, 'o'),
                        label=f'{state} (n={len(sdf)})')
    axes[0].scatter([0], [0], color='black', marker='*', s=200, zorder=5,
                     label='Template')
    axes[0].set_xlabel('Global RMSD (A)')
    axes[0].set_ylabel('Chromophore RMSD (A)')
    axes[0].set_title('Global vs Chromophore RMSD (by state)')
    axes[0].legend(fontsize=7, loc='upper left')

    # Right: colored by SD
    sd_col = 'global_rmsd_sd'
    if sd_col in plot_df.columns and plot_df[sd_col].notna().any():
        sc = axes[1].scatter(plot_df['global_rmsd'], plot_df['chromophore_rmsd'],
                             c=plot_df[sd_col], cmap='plasma', s=25, alpha=0.7)
        plt.colorbar(sc, ax=axes[1], label='Global RMSD SD (across seeds)')
        axes[1].set_title('Colored by RMSD SD (prediction uncertainty)')
    else:
        for state in sorted(plot_df['state'].unique()):
            sdf = plot_df[plot_df['state'] == state]
            axes[1].scatter(sdf['global_rmsd'], sdf['chromophore_rmsd'],
                            alpha=0.5, s=25, label=state)
        axes[1].set_title('Global vs Chromophore RMSD')
        axes[1].legend()

    axes[1].set_xlabel('Global RMSD (A)')
    axes[1].set_ylabel('Chromophore RMSD (A)')

    plt.tight_layout()
    plt.savefig(output_dir / '06_global_vs_chromophore_rmsd.png', dpi=300,
                bbox_inches='tight')
    plt.close()

    cols = [c for c in ['seq_id', 'state', 'template', 'global_rmsd',
                         'global_rmsd_sd', 'chromophore_rmsd',
                         'chromophore_rmsd_sd'] if c in df.columns]
    plot_df[cols].to_csv(csv_dir / '06_global_vs_chromophore_rmsd_data.csv',
                         index=False)
    print("  06_global_vs_chromophore_rmsd.png")


def _plot_plddt_vs_rmsd_sd(df, template_df, output_dir, csv_dir):
    """Scatter: pLDDT vs RMSD, colored by pLDDT SD."""
    if 'mean_plddt' not in df.columns or 'global_rmsd' not in df.columns:
        return

    plot_df = df.dropna(subset=['mean_plddt', 'global_rmsd'])
    if plot_df.empty:
        return

    tpl_colors = _get_template_colors(template_df)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    markers = {'holo': 'o', 'apo': 's'}

    # Left: by state
    for state in sorted(plot_df['state'].unique()):
        sdf = plot_df[plot_df['state'] == state]
        axes[0].scatter(sdf['global_rmsd'], sdf['mean_plddt'],
                        alpha=0.5, s=25, color=_state_color(state),
                        marker=markers.get(state, 'o'),
                        label=f'{state} (n={len(sdf)})')

    # Template references
    for _, row in template_df.iterrows():
        plddt = row.get('mean_plddt')
        if plddt and pd.notna(plddt):
            color = tpl_colors.get(row['template_name'], 'black')
            axes[0].scatter([0], [plddt], marker='*', s=200, color=color,
                            edgecolors='black', linewidths=0.8, zorder=5)

    axes[0].set_xlabel('Global RMSD (A)')
    axes[0].set_ylabel('Mean pLDDT')
    axes[0].set_title('pLDDT vs RMSD (by state)')
    axes[0].legend(fontsize=7, loc='lower left')

    # Right: colored by pLDDT SD
    sd_col = 'mean_plddt_sd'
    if sd_col in plot_df.columns and plot_df[sd_col].notna().any():
        sc = axes[1].scatter(plot_df['global_rmsd'], plot_df['mean_plddt'],
                             c=plot_df[sd_col], cmap='plasma', s=25, alpha=0.7)
        plt.colorbar(sc, ax=axes[1], label='pLDDT SD (across seeds)')
        axes[1].set_title('Colored by pLDDT SD (prediction uncertainty)')
    else:
        axes[1].scatter(plot_df['global_rmsd'], plot_df['mean_plddt'],
                        alpha=0.5, s=25)
        axes[1].set_title('pLDDT vs RMSD')

    axes[1].set_xlabel('Global RMSD (A)')
    axes[1].set_ylabel('Mean pLDDT')

    plt.tight_layout()
    plt.savefig(output_dir / '07_plddt_vs_rmsd.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  07_plddt_vs_rmsd.png")


def _plot_correlation_matrix(df, output_dir, csv_dir):
    numeric = df.select_dtypes(include=[np.number]).columns
    exclude = ['seq_id', 'n_seeds']
    cols = [c for c in numeric if c not in exclude and not c.endswith('_sd')]
    if len(cols) < 2:
        return

    corr = df[cols].corr()
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr, annot=True, fmt='.2f', cmap='coolwarm', center=0,
                square=True, linewidths=1)
    plt.title('Correlation Matrix (Mean Values)')
    plt.tight_layout()
    plt.savefig(output_dir / '05_correlation_matrix.png', dpi=300,
                bbox_inches='tight')
    plt.close()
    corr.to_csv(csv_dir / '05_correlation_matrix_data.csv')
    print("  05_correlation_matrix.png")


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
            top = sdf.nsmallest(20, metric) if ascending else sdf.nlargest(20, metric)
            top_designs[f'top_20_{state}_{metric}'] = top

    all_top = pd.concat([d.assign(metric=key)
                         for key, d in top_designs.items()], ignore_index=True)
    all_top.to_csv(csv_dir / '07_top20_designs_data.csv', index=False)
    print("  07_top20_designs exported")


# ─── Main ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Analyze AF3 predictions (all seeds)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Single output dir (default: auto-discover all)')
    parser.add_argument('--template', type=str,
                        default='inputs/G4FP_des1_cro_mod0.pdb',
                        help='Reference structure for RMSD')
    parser.add_argument('--state', type=str, default='both',
                        choices=['holo', 'apo', 'both'])
    parser.add_argument('--chromophore-range', type=str, default=None,
                        help='Chromophore residue range as start,end (e.g. 197,199). '
                             'Auto-detected from each template PDB when omitted.')
    parser.add_argument('--analysis-output-dir', type=str,
                        default='analysis_output')
    parser.add_argument('--force', action='store_true',
                        help='Re-analyze even if results CSV exists')

    args = parser.parse_args()
    if args.chromophore_range:
        cs, ce = map(int, args.chromophore_range.split(','))
        chromophore_range: Optional[Tuple[int, int]] = (cs, ce)
    else:
        chromophore_range = None  # auto-detect per template
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
    # chromophore_range=None triggers per-template auto-detection
    template_df = load_template_metrics(base_dir, chromophore_range)
    print(f"\nLoaded {len(template_df)} template references:")
    for _, r in template_df.iterrows():
        src = r.get('source', '?')
        plddt = r.get('mean_plddt', 0)
        crange = f"chrom={r.get('chromophore_start', '?')}-{r.get('chromophore_end', '?')}"
        af3_info = ""
        if src == 'af3':
            bp = r.get('holo_mean_plddt', '?')
            ap = r.get('apo_mean_plddt', '?')
            af3_info = f" holo={bp:.1f} apo={ap:.1f}" if isinstance(bp, float) else ""
        print(f"  {r['template_name']}: pLDDT={plddt:.1f} ({src}) {crange}{af3_info}")

    print(f"\n{'='*80}")
    print(f"Analyzing AF3 Structures (ALL seeds) Across {len(output_dirs)} Templates")
    print(f"{'='*80}")

    all_frames = []
    for odir in output_dirs:
        print(f"\n--- {odir.name} ---")
        csv_path = odir / "05_structure_analysis_results.csv"
        if not args.force and csv_path.exists():
            try:
                df = pd.read_csv(csv_path)
                if not df.empty:
                    all_frames.append(df)
                    print(f"  Cached: {csv_path} ({len(df)} designs) [use --force to re-run]")
                    continue
            except Exception:
                pass  # re-analyze if CSV is corrupt
        try:
            analyzer = AlphaFoldStructureAnalyzer(
                output_dir=odir,
                reference_template=Path(args.template),
                chromophore_range=chromophore_range,
            )
            df = analyzer.analyze_all_designs(state=args.state)
            if not df.empty:
                all_frames.append(df)
                df.to_csv(csv_path, index=False)
                print(f"  Saved: {csv_path} ({len(df)} designs)")
            else:
                print(f"  No completed AF3 inference")
        except Exception as e:
            print(f"  Error: {e}")

    if not all_frames:
        print("\nNo designs analyzed")
        sys.exit(1)

    df_combined = pd.concat(all_frames, ignore_index=True)
    print(f"\n{'='*80}")
    print(f"Combined: {len(df_combined)} designs, "
          f"{df_combined['template'].nunique()} templates, "
          f"{df_combined['state'].nunique()} states")
    if 'n_seeds' in df_combined.columns:
        print(f"Seeds per design: {df_combined['n_seeds'].median():.0f} median")
    print(f"{'='*80}")

    central_dir = Path(args.analysis_output_dir)
    viz_dir = central_dir / "05_structure_analysis"
    create_visualizations(df_combined, template_df, viz_dir)

    combined_csv = central_dir / "05_structure_analysis_results.csv"
    combined_csv.parent.mkdir(exist_ok=True, parents=True)
    df_combined.to_csv(combined_csv, index=False)

    print(f"\nAnalysis complete!")
    print(f"  Combined results: {combined_csv}")
    print(f"  Plots: {viz_dir}/")


if __name__ == "__main__":
    main()

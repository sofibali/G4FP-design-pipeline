#!/usr/bin/env python3
"""
Step 7: Aggregate results across all template structures and rank candidates.

Pools per-structure holo-vs-apo comparison data from 06_compare_ligand_states.py,
computes a G4FP fitness score, and selects the top 100-1000 candidates using:
  1. Hard filters (holo must fold, apo should be disordered)
  2. Pareto-optimal frontier across multiple objectives
  3. Diversity-weighted selection to avoid redundant sequences

Now handles per-seed SD columns from the updated 06 script and loads template
AF3 predictions (holo + apo separately) for proper reference overlays with
error bars.

Outputs:
  07_results/
    07_final_candidates.csv      -- ranked candidates with all metrics + SD
    07_final_candidates.fa       -- FASTA of selected sequences
    07_pareto_frontier.csv       -- Pareto-optimal designs
    07_selection_summary.png     -- overview plots
    07_holo_vs_apo_scatter.png   -- holo pLDDT vs apo pLDDT (colored by fitness)
    07_holo_vs_apo_sd.png        -- holo pLDDT vs apo pLDDT (colored by SD)
    07_sd_overview.png           -- SD distributions for key metrics
    07_pareto_plot.png           -- Pareto frontier visualization
    plot_data_csvs/              -- raw data for Prism

Usage:
    python 07_aggregate_and_rank.py
    python 07_aggregate_and_rank.py --n-select 500
    python 07_aggregate_and_rank.py --output-dir output_G4FP_des1_cro_mod0
"""

import sys
import json
import argparse
import warnings
from pathlib import Path
from typing import List, Tuple, Optional, Dict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster, linkage

warnings.filterwarnings("ignore")
sns.set_style("whitegrid")
plt.rcParams["figure.figsize"] = (12, 8)
plt.rcParams["font.size"] = 10


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_comparison_results(output_dir: Path, base_dir: Path = None) -> Optional[pd.DataFrame]:
    """Load 06_ligand_state_comparison_results.csv for one output directory.

    Checks both the output dir itself and analysis_output/output_dir_name/.
    Only loads per-template CSV (not the combined one).
    """
    csv_path = output_dir / "06_ligand_state_comparison_results.csv"
    if not csv_path.exists() and base_dir is not None:
        csv_path = base_dir / "analysis_output" / output_dir.name / "06_ligand_state_comparison_results.csv"
    if not csv_path.exists():
        return None
    df = pd.read_csv(csv_path)
    if "template" not in df.columns:
        df["template"] = output_dir.name
    return df


def load_filtered_sequences(output_dir: Path) -> Optional[pd.DataFrame]:
    """Load 02_filtered_sequences.csv to get actual protein sequences."""
    csv_path = output_dir / "02_filtered_sequences.csv"
    if not csv_path.exists():
        return None
    return pd.read_csv(csv_path)


def aggregate_all(base_dir: Path, output_dirs: Optional[List[Path]] = None) -> pd.DataFrame:
    """Merge comparison results across all output directories."""
    if output_dirs is None:
        output_dirs = sorted(
            d for d in base_dir.iterdir()
            if d.is_dir() and d.name.startswith("output_G4FP_")
        )

    frames = []
    for odir in output_dirs:
        comparison_df = load_comparison_results(odir, base_dir)
        seq_df = load_filtered_sequences(odir)

        if comparison_df is None:
            print(f"  SKIP {odir.name}: no 06_ligand_state_comparison_results.csv")
            continue

        # Merge in sequences if available
        if seq_df is not None and "sequence" in seq_df.columns:
            comparison_df = comparison_df.merge(
                seq_df[["seq_id", "sequence"]].drop_duplicates(),
                on="seq_id", how="left"
            )

        frames.append(comparison_df)
        n_with_sd = sum(1 for c in comparison_df.columns if c.endswith("_sd"))
        print(f"  Loaded {odir.name}: {len(comparison_df)} designs, {n_with_sd} SD columns")

    if not frames:
        print("ERROR: No comparison results found in any output directory.")
        print("       Run 06_compare_ligand_states.py first.")
        sys.exit(1)

    df = pd.concat(frames, ignore_index=True)
    df["global_id"] = df["template"] + "_design_" + df["seq_id"].astype(str).str.zfill(4)
    print(f"\nTotal designs across all templates: {len(df)}")
    return df


# ---------------------------------------------------------------------------
# Template AF3 reference loading
# ---------------------------------------------------------------------------

def _find_af3_result_dir(base_dir: Path) -> Optional[Path]:
    """Find AF3 result directory, handling nested outputs."""
    if not base_dir.exists():
        return None
    # Check for CIF directly
    cifs = list(base_dir.glob("*_model.cif"))
    if cifs:
        return base_dir
    # Check one level down (AF3 sometimes nests)
    for sub in base_dir.iterdir():
        if sub.is_dir():
            cifs = list(sub.glob("*_model.cif"))
            if cifs:
                return sub
    return None


def _discover_seed_dirs(result_dir: Path) -> List[Path]:
    """Find all seed-N_sample-M directories under a result dir."""
    seeds = []
    for d in sorted(result_dir.iterdir()):
        if d.is_dir() and d.name.startswith("seed-"):
            if list(d.glob("model.cif")) or list(d.glob("*_model.cif")):
                seeds.append(d)
    return seeds


def _analyze_template_seed(seed_dir: Path, chain_id: str = 'A') -> Optional[Dict]:
    """Analyze a single seed-sample model from template AF3 output."""
    cif = list(seed_dir.glob("model.cif")) or list(seed_dir.glob("*_model.cif"))
    if not cif:
        return None

    result = {}

    # Load confidence files
    conf_files = list(seed_dir.glob("confidences.json")) + \
                 list(seed_dir.glob("*_confidences.json"))
    conf_files = [f for f in conf_files if 'summary' not in f.name]

    summ_files = list(seed_dir.glob("summary_confidences.json")) + \
                 list(seed_dir.glob("*_summary_confidences.json"))

    if summ_files:
        with open(summ_files[0]) as f:
            sc = json.load(f)
        for k in ['ptm', 'iptm', 'ranking_score', 'fraction_disordered']:
            if k in sc:
                result[k] = sc[k]

    if conf_files:
        with open(conf_files[0]) as f:
            fc = json.load(f)
        plddts = fc.get('atom_plddts', [])
        chain_ids = fc.get('atom_chain_ids', [])
        if plddts and chain_ids:
            chain_plddts = [p for p, c in zip(plddts, chain_ids)
                            if c == chain_id]
            if chain_plddts:
                result['mean_plddt'] = float(np.mean(chain_plddts))

        # PAE
        pae = fc.get('pae', [])
        if pae:
            result['mean_pae'] = float(np.mean(pae))

    return result if result else None


def load_template_metrics(base_dir: Path) -> pd.DataFrame:
    """Load template AF3 predictions (holo + apo) for reference overlays.

    Checks template_af3_outputs/{bound,apo}/ for AF3 results with per-seed
    analysis. Falls back to B-factor pLDDT from input PDBs.
    """
    template_af3_dir = base_dir / "template_af3_outputs"
    records = []

    if template_af3_dir.exists():
        # Process each template's holo and apo AF3 outputs
        for disk_state in ['bound', 'apo']:
            display_state = 'holo' if disk_state == 'bound' else disk_state
            state_dir = template_af3_dir / disk_state
            if not state_dir.exists():
                continue
            for template_dir in sorted(state_dir.iterdir()):
                if not template_dir.is_dir():
                    continue
                result_dir = _find_af3_result_dir(template_dir)
                if result_dir is None:
                    continue

                seed_dirs = _discover_seed_dirs(result_dir)
                if not seed_dirs:
                    # Try top-level only
                    m = _analyze_template_seed(result_dir)
                    if m:
                        seed_dirs = [result_dir]
                    else:
                        continue

                # Analyze all seeds
                metrics = []
                for sd in seed_dirs:
                    m = _analyze_template_seed(sd)
                    if m:
                        metrics.append(m)

                if not metrics:
                    continue

                # Extract template name from dir name
                # e.g., template_G4FP_des1_cro_mod0_bound -> G4FP_des1_cro_mod0
                tname = template_dir.name
                for suffix in ['_bound', '_apo']:
                    if tname.endswith(suffix):
                        tname = tname[:-len(suffix)]
                if tname.startswith('template_'):
                    tname = tname[len('template_'):]

                df_m = pd.DataFrame(metrics)
                rec = {'template_name': tname}
                for col in df_m.columns:
                    vals = df_m[col].dropna()
                    if len(vals) > 0:
                        rec[f'{display_state}_{col}'] = float(vals.mean())
                        rec[f'{display_state}_{col}_sd'] = (
                            float(vals.std()) if len(vals) > 1 else 0.0)

                # Find existing record for this template and merge
                found = False
                for r in records:
                    if r['template_name'] == tname:
                        r.update(rec)
                        found = True
                        break
                if not found:
                    records.append(rec)

    # If no AF3 data, fall back to B-factor pLDDT from input PDBs
    if not records:
        try:
            from Bio.PDB import PDBParser
            parser = PDBParser(QUIET=True)
            inputs_dir = base_dir / "inputs"
            for pdb in sorted(inputs_dir.glob("G4FP_*.pdb")):
                structure = parser.get_structure('s', str(pdb))
                b_factors = []
                for model in structure:
                    for chain in model:
                        if chain.get_id() == 'A':
                            for residue in chain:
                                if residue.has_id('CA'):
                                    b_factors.append(residue['CA'].get_bfactor())
                if b_factors:
                    mean_plddt = float(np.mean(b_factors))
                    records.append({
                        'template_name': pdb.stem,
                        'holo_mean_plddt': mean_plddt,
                        'holo_mean_plddt_sd': 0.0,
                        'apo_mean_plddt': mean_plddt,  # same: no separate apo AF3
                        'apo_mean_plddt_sd': 0.0,
                        'source': 'bfactor',
                    })
        except ImportError:
            pass

    df = pd.DataFrame(records)
    if not df.empty:
        has_af3 = 'source' not in df.columns or (df.get('source', '') != 'bfactor').any()
        src = "AF3" if has_af3 else "B-factor"
        print(f"  Loaded {len(df)} template references ({src})")
    return df


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def compute_relative_apo_plddt(df: pd.DataFrame,
                                template_df: pd.DataFrame) -> pd.DataFrame:
    """Compute apo pLDDT relative to each design's parent template.

    Positive = design apo is more stable than template (bad for switch).
    Negative = design apo is less stable than template (good).
    """
    df = df.copy()
    df["apo_plddt_vs_template"] = np.nan

    if template_df is None or template_df.empty:
        return df

    # Build template -> apo_mean_plddt map
    tpl_apo = {}
    for _, row in template_df.iterrows():
        tname = row["template_name"]
        ap = row.get("apo_mean_plddt")
        if pd.notna(ap):
            tpl_apo[tname] = ap

    for idx, row in df.iterrows():
        tpl = str(row.get("template", "")).replace("output_", "")
        ap = row.get("apo_mean_plddt")
        if pd.isna(ap):
            continue
        ref = tpl_apo.get(tpl)
        if ref is not None:
            df.at[idx, "apo_plddt_vs_template"] = ap - ref

    return df


def compute_fitness_score(df: pd.DataFrame,
                          w_chrom_holo: float = 0.25,
                          w_iptm: float = 0.25,
                          w_rmsd: float = 0.15,
                          w_holo_ptm: float = 0.15,
                          w_rel_apo: float = 0.10,
                          w_confidence: float = 0.05,
                          w_chrom_diff: float = 0.05) -> pd.DataFrame:
    """
    Compute G4FP fitness score prioritizing chromophore and ligand interface.

    Components (all normalized 0-1 before weighting):
      - chrom_holo:   holo chromophore pLDDT, higher = better positioned
      - iptm:         holo iPTM, higher = better protein-DNA interface
      - rmsd:         global RMSD holo vs apo, higher = more structural change
      - holo_ptm:     holo pTM, higher = better fold
      - rel_apo:      apo pLDDT relative to template (lower = better switch)
      - confidence:   1 - normalized(SD), reproducible predictions
      - chrom_diff:   chromophore pLDDT diff (holo - apo)
    """
    df = df.copy()

    def _norm(series):
        smin, smax = series.min(), series.max()
        if smax == smin:
            return pd.Series(0.5, index=series.index)
        return (series - smin) / (smax - smin)

    # Chromophore holo pLDDT (higher = better)
    chrom_col = "holo_chromophore_plddt"
    if chrom_col not in df.columns:
        chrom_col = "chromophore_holo_plddt"
    if chrom_col in df.columns:
        df["_n_chrom_holo"] = _norm(df[chrom_col].fillna(df[chrom_col].median()))
    else:
        df["_n_chrom_holo"] = 0.5

    # iPTM: protein-DNA interface (higher = better)
    if "holo_iptm" in df.columns:
        df["_n_iptm"] = _norm(df["holo_iptm"].fillna(df["holo_iptm"].median()))
    else:
        df["_n_iptm"] = 0.5

    # RMSD (higher = more structural change)
    if "global_rmsd" in df.columns:
        df["_n_rmsd"] = _norm(df["global_rmsd"].fillna(0))
    else:
        df["_n_rmsd"] = 0

    # Holo pTM (higher = better fold)
    if "holo_ptm" in df.columns:
        df["_n_holo_ptm"] = _norm(df["holo_ptm"].fillna(0))
    else:
        df["_n_holo_ptm"] = 0

    # Relative apo pLDDT (lower = better switch, inverted)
    if "apo_plddt_vs_template" in df.columns:
        rel = df["apo_plddt_vs_template"].fillna(0)
        df["_n_rel_apo"] = 1.0 - _norm(rel)  # invert: lower relative = higher score
    elif "apo_mean_plddt" in df.columns:
        df["_n_rel_apo"] = 1.0 - _norm(df["apo_mean_plddt"].fillna(100))
    else:
        df["_n_rel_apo"] = 0.5

    # Confidence (lower SD = better)
    if "mean_plddt_diff_sd" in df.columns:
        sd_vals = df["mean_plddt_diff_sd"].fillna(df["mean_plddt_diff_sd"].median())
        df["_n_confidence"] = 1.0 - _norm(sd_vals)
    else:
        df["_n_confidence"] = 0.5

    # Chromophore pLDDT diff
    if "chromophore_plddt_diff" in df.columns:
        df["_n_chrom_diff"] = _norm(df["chromophore_plddt_diff"].fillna(0))
    else:
        df["_n_chrom_diff"] = 0

    df["fitness_score"] = (
        w_chrom_holo * df["_n_chrom_holo"]
        + w_iptm * df["_n_iptm"]
        + w_rmsd * df["_n_rmsd"]
        + w_holo_ptm * df["_n_holo_ptm"]
        + w_rel_apo * df["_n_rel_apo"]
        + w_confidence * df["_n_confidence"]
        + w_chrom_diff * df["_n_chrom_diff"]
    )

    # Keep normalized components for diagnostic plots
    component_cols = [c for c in df.columns if c.startswith("_n_")]
    df.rename(columns={c: f"score_component{c[2:]}" for c in component_cols},
              inplace=True)

    return df


# ---------------------------------------------------------------------------
# Hard filters
# ---------------------------------------------------------------------------

def apply_hard_filters(df: pd.DataFrame,
                       min_holo_plddt: float = 70.0,
                       min_holo_ptm: float = 0.5,
                       min_holo_iptm: float = 0.3,
                       max_apo_plddt: float = None) -> pd.DataFrame:
    """Apply hard cutoffs. No absolute apo filter by default (use relative ranking)."""
    n_before = len(df)
    mask = pd.Series(True, index=df.index)

    if "holo_mean_plddt" in df.columns:
        mask &= df["holo_mean_plddt"] >= min_holo_plddt
    if "holo_ptm" in df.columns:
        mask &= df["holo_ptm"] >= min_holo_ptm
    if "holo_iptm" in df.columns and min_holo_iptm > 0:
        mask &= df["holo_iptm"] >= min_holo_iptm
    if max_apo_plddt is not None and "apo_mean_plddt" in df.columns:
        mask &= df["apo_mean_plddt"] <= max_apo_plddt

    df_filtered = df[mask].copy()
    n_after = len(df_filtered)
    filters_str = f"holo_plddt >= {min_holo_plddt}, holo_ptm >= {min_holo_ptm}"
    if min_holo_iptm > 0:
        filters_str += f", holo_iptm >= {min_holo_iptm}"
    if max_apo_plddt is not None:
        filters_str += f", apo_plddt <= {max_apo_plddt}"
    print(f"\nHard filters: {n_before} -> {n_after} designs "
          f"(removed {n_before - n_after})")
    print(f"  {filters_str}")

    if n_after == 0:
        print("\nWARNING: All designs filtered out! Relaxing to top 50%.")
        df_sorted = df.sort_values("fitness_score", ascending=False)
        df_filtered = df_sorted.head(len(df_sorted) // 2).copy()
        print(f"  Kept {len(df_filtered)} designs by fitness score fallback.")

    return df_filtered


# ---------------------------------------------------------------------------
# Pareto frontier
# ---------------------------------------------------------------------------

def pareto_frontier(df: pd.DataFrame,
                    objectives: List[Tuple[str, bool]]) -> pd.DataFrame:
    """Find Pareto-optimal (non-dominated) designs across objectives."""
    available = [(col, mx) for col, mx in objectives if col in df.columns]
    if not available:
        print("  WARNING: no Pareto objectives available, returning all")
        return df

    vals = np.zeros((len(df), len(available)))
    for j, (col, maximize) in enumerate(available):
        v = df[col].fillna(df[col].median()).values
        if not maximize:
            v = -v
        vals[:, j] = v

    is_pareto = np.ones(len(df), dtype=bool)
    for i in range(len(df)):
        if not is_pareto[i]:
            continue
        for j in range(len(df)):
            if i == j or not is_pareto[j]:
                continue
            if np.all(vals[j] >= vals[i]) and np.any(vals[j] > vals[i]):
                is_pareto[i] = False
                break

    pareto_df = df[is_pareto].copy()
    pareto_df["is_pareto"] = True
    print(f"\nPareto frontier: {len(pareto_df)} non-dominated designs "
          f"(from {len(df)} total)")
    obj_strs = [f"{c} ({'max' if m else 'min'})" for c, m in available]
    print(f"  Objectives: {obj_strs}")
    return pareto_df


# ---------------------------------------------------------------------------
# Diversity-weighted selection
# ---------------------------------------------------------------------------

def sequence_hamming_distance(seqs: List[str]) -> np.ndarray:
    """Compute pairwise Hamming distance matrix for aligned sequences."""
    n = len(seqs)
    max_len = max(len(s) for s in seqs)
    padded = [s.ljust(max_len, "X") for s in seqs]
    arr = np.array([list(s) for s in padded])
    dist = np.zeros((n, n))
    for i in range(n):
        dist[i, :] = np.sum(arr[i] != arr, axis=1) / max_len
    return dist


def diversity_weighted_selection(df: pd.DataFrame,
                                  n_select: int,
                                  score_col: str = "fitness_score",
                                  seq_col: str = "sequence",
                                  n_clusters: int = 0) -> pd.DataFrame:
    """Select top designs while maximizing sequence diversity."""
    if seq_col not in df.columns or df[seq_col].isna().all():
        print("  WARNING: no sequence data, using score only")
        return df.nlargest(n_select, score_col).copy()

    df = df.copy()
    seqs = df[seq_col].tolist()
    n = len(seqs)

    if n <= n_select:
        print(f"  Only {n} designs available, returning all")
        df["cluster"] = 0
        return df

    print(f"  Computing pairwise distances for {n} sequences...")
    dist_matrix = sequence_hamming_distance(seqs)

    condensed = squareform(dist_matrix)
    Z = linkage(condensed, method="average")

    if n_clusters <= 0:
        n_clusters = min(n_select, max(10, n // 5))
    n_clusters = min(n_clusters, n)

    labels = fcluster(Z, t=n_clusters, criterion="maxclust")
    df["cluster"] = labels
    print(f"  Clustered into {len(set(labels))} groups")

    selected_indices = []
    remaining_by_cluster = {}

    for cluster_id in sorted(set(labels)):
        cluster_df = df[df["cluster"] == cluster_id].sort_values(
            score_col, ascending=False)
        indices = cluster_df.index.tolist()
        if indices:
            selected_indices.append(indices[0])
            remaining_by_cluster[cluster_id] = indices[1:]

    if len(selected_indices) >= n_select:
        selected_indices = selected_indices[:n_select]
    else:
        while len(selected_indices) < n_select:
            added = False
            for cid in sorted(remaining_by_cluster.keys()):
                if len(selected_indices) >= n_select:
                    break
                if remaining_by_cluster[cid]:
                    selected_indices.append(remaining_by_cluster[cid].pop(0))
                    added = True
            if not added:
                break

    result = df.loc[selected_indices].copy()
    result["diversity_rank"] = range(1, len(result) + 1)
    print(f"  Selected {len(result)} diverse designs from "
          f"{len(set(result['cluster']))} clusters")
    return result


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def _plot_template_references(ax, template_df, x_col, y_col,
                              x_sd_col=None, y_sd_col=None):
    """Overlay template AF3 references as stars with optional error bars."""
    if template_df is None or template_df.empty:
        return
    if x_col not in template_df.columns or y_col not in template_df.columns:
        return

    cmap = plt.cm.get_cmap('tab10', max(10, len(template_df)))
    for i, (_, row) in enumerate(template_df.iterrows()):
        x = row[x_col]
        y = row[y_col]
        xerr = row.get(x_sd_col, 0) if x_sd_col else 0
        yerr = row.get(y_sd_col, 0) if y_sd_col else 0
        short = row['template_name'].replace('G4FP_', '')

        if xerr > 0 or yerr > 0:
            ax.errorbar(x, y, xerr=xerr if xerr > 0 else None,
                        yerr=yerr if yerr > 0 else None,
                        fmt='none', ecolor=cmap(i), elinewidth=1.5,
                        capsize=3, capthick=1, alpha=0.7, zorder=6)

        ax.scatter([x], [y], marker='*', s=200,
                   color=cmap(i), edgecolors='black', linewidths=0.8,
                   zorder=7)
        ax.annotate(short, (x, y), fontsize=5,
                    xytext=(4, 4), textcoords='offset points')

    ax.scatter([], [], marker='*', s=100, color='black', label='Templates')


def _plot_top10_per_template_heatmap(df_all, df_selected, df_pareto,
                                     output_dir, csv_dir):
    """Heatmap of top 10 designs per template, columns sorted by fitness weight.

    Only includes metrics used in the fitness score + fraction_disordered.
    """
    if "fitness_score" not in df_all.columns:
        return

    # Collect top 10 per template by fitness score
    frames = []
    for tpl in sorted(df_all["template"].unique()):
        tpl_df = df_all[df_all["template"] == tpl].nlargest(10, "fitness_score")
        frames.append(tpl_df)
    if not frames:
        return
    top_per_tpl = pd.concat(frames).sort_values("fitness_score", ascending=False)

    # Columns ordered by fitness weight importance, then diagnostics
    # Resolve chromophore column name
    chrom_col = ("holo_chromophore_plddt" if "holo_chromophore_plddt" in top_per_tpl.columns
                 else "chromophore_holo_plddt")
    metrics_ordered = [
        ('fitness_score', 'Fitness\nScore', False),
        (chrom_col, 'Chrom\npLDDT\n(25%)', False),
        ('holo_iptm', 'iPTM\n(25%)', False),
        ('global_rmsd', 'RMSD\n(15%)', False),
        ('holo_ptm', 'pTM\n(15%)', False),
        ('apo_plddt_vs_template', 'Apo pLDDT\nvs Tpl\n(10%)', True),
        ('mean_plddt_diff_sd', 'pLDDT Diff\nSD (5%)', True),
        ('chromophore_plddt_diff', 'Chrom\nDiff (5%)', False),
        ('holo_fraction_disordered', 'Holo\nDisorder', True),
        ('apo_fraction_disordered', 'Apo\nDisorder', True),
    ]
    # Filter to available columns
    available = [(c, l, inv) for c, l, inv in metrics_ordered
                 if c in top_per_tpl.columns and top_per_tpl[c].notna().any()]
    if len(available) < 3:
        return

    all_cols = [c for c, _, _ in available]
    col_labels = [l for _, l, _ in available]
    invert_flags = [inv for _, _, inv in available]

    # Row labels
    labels = []
    for _, row in top_per_tpl.iterrows():
        tpl = str(row.get('template', '')).replace('output_G4FP_', '').replace('G4FP_', '')
        sid = int(row['seq_id'])
        labels.append(f"d{sid:03d} {tpl}")

    mat = top_per_tpl[all_cols].values.astype(float)

    # Normalize 0-1 across ALL 1000 designs (not just top) for fair color
    col_min = np.array([df_all[c].min() for c in all_cols])
    col_max = np.array([df_all[c].max() for c in all_cols])
    col_range = col_max - col_min
    col_range[col_range == 0] = 1
    mat_norm = (mat - col_min) / col_range
    for j, inv in enumerate(invert_flags):
        if inv:
            mat_norm[:, j] = 1.0 - mat_norm[:, j]

    # Draw template separators
    n_rows = len(labels)
    fig_h = max(10, n_rows * 0.35)
    fig, ax = plt.subplots(figsize=(max(12, len(all_cols) * 1.2), fig_h))
    im = ax.imshow(mat_norm, aspect='auto', cmap='RdYlGn', vmin=0, vmax=1)

    ax.set_xticks(range(len(all_cols)))
    ax.set_xticklabels(col_labels, fontsize=8, ha='center')
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(labels, fontsize=7)

    # Annotate cells with values
    for i in range(n_rows):
        for j in range(len(all_cols)):
            val = mat[i, j]
            if np.isnan(val):
                continue
            nv = mat_norm[i, j]
            text_color = 'white' if nv < 0.25 or nv > 0.85 else 'black'
            fmt = '.0f' if abs(val) > 10 else '.2f' if abs(val) > 1 else '.3f'
            ax.text(j, i, f'{val:{fmt}}', ha='center', va='center',
                    fontsize=5.5, color=text_color)

    # Draw horizontal lines between templates
    prev_tpl = None
    for i, (_, row) in enumerate(top_per_tpl.iterrows()):
        tpl = row.get('template', '')
        if prev_tpl is not None and tpl != prev_tpl:
            ax.axhline(i - 0.5, color='black', linewidth=1.5)
        prev_tpl = tpl

    # Mark Pareto and selected designs
    pareto_ids = set(df_pareto.get('global_id', pd.Series())) if not df_pareto.empty else set()
    selected_ids = set(df_selected.get('global_id', pd.Series())) if not df_selected.empty else set()
    for i, (_, row) in enumerate(top_per_tpl.iterrows()):
        gid = row.get('global_id', '')
        markers = []
        if gid in pareto_ids:
            markers.append('P')
        if gid in selected_ids:
            markers.append('*')
        if markers:
            ax.text(-0.7, i, ''.join(markers), ha='center', va='center',
                    fontsize=7, fontweight='bold',
                    color='red' if 'P' in markers else 'blue')

    ax.set_title('Top 10 Per Template (sorted by fitness weight, green=better)\n'
                 'P=Pareto  *=Selected', fontsize=11)
    plt.colorbar(im, ax=ax, label='Normalized (1=best)', shrink=0.5)
    plt.tight_layout()
    plt.savefig(output_dir / "07_top10_per_template_heatmap.png",
                dpi=300, bbox_inches="tight")
    plt.close()

    top_per_tpl[['global_id', 'template', 'seq_id'] + all_cols].to_csv(
        csv_dir / "07_top10_per_template_heatmap_data.csv", index=False)
    print("  07_top10_per_template_heatmap.png")


def _plot_top50_global_heatmap(df_all, df_selected, df_pareto,
                                output_dir, csv_dir):
    """Heatmap of top 50 designs globally by fitness score."""
    if "fitness_score" not in df_all.columns:
        return

    top50 = df_all.nlargest(50, "fitness_score")

    chrom_col = ("holo_chromophore_plddt" if "holo_chromophore_plddt" in top50.columns
                 else "chromophore_holo_plddt")
    metrics_ordered = [
        ('fitness_score', 'Fitness\nScore', False),
        (chrom_col, 'Chrom\npLDDT\n(25%)', False),
        ('holo_iptm', 'iPTM\n(25%)', False),
        ('global_rmsd', 'RMSD\n(15%)', False),
        ('holo_ptm', 'pTM\n(15%)', False),
        ('apo_plddt_vs_template', 'Apo pLDDT\nvs Tpl\n(10%)', True),
        ('mean_plddt_diff_sd', 'pLDDT Diff\nSD (5%)', True),
        ('chromophore_plddt_diff', 'Chrom\nDiff (5%)', False),
        ('holo_fraction_disordered', 'Holo\nDisorder', True),
        ('apo_fraction_disordered', 'Apo\nDisorder', True),
    ]
    available = [(c, l, inv) for c, l, inv in metrics_ordered
                 if c in top50.columns and top50[c].notna().any()]
    if len(available) < 3:
        return

    all_cols = [c for c, _, _ in available]
    col_labels = [l for _, l, _ in available]
    invert_flags = [inv for _, _, inv in available]

    labels = []
    for i, (_, row) in enumerate(top50.iterrows()):
        tpl = str(row.get('template', '')).replace('output_G4FP_', '').replace('G4FP_', '')
        sid = int(row['seq_id'])
        labels.append(f"#{i+1} d{sid:03d} {tpl}")

    mat = top50[all_cols].values.astype(float)

    # Normalize across all designs for fair color
    col_min = np.array([df_all[c].min() for c in all_cols])
    col_max = np.array([df_all[c].max() for c in all_cols])
    col_range = col_max - col_min
    col_range[col_range == 0] = 1
    mat_norm = (mat - col_min) / col_range
    for j, inv in enumerate(invert_flags):
        if inv:
            mat_norm[:, j] = 1.0 - mat_norm[:, j]

    fig, ax = plt.subplots(figsize=(max(12, len(all_cols) * 1.2), 18))
    im = ax.imshow(mat_norm, aspect='auto', cmap='RdYlGn', vmin=0, vmax=1)

    ax.set_xticks(range(len(all_cols)))
    ax.set_xticklabels(col_labels, fontsize=8, ha='center')
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)

    for i in range(len(labels)):
        for j in range(len(all_cols)):
            val = mat[i, j]
            if np.isnan(val):
                continue
            nv = mat_norm[i, j]
            text_color = 'white' if nv < 0.25 or nv > 0.85 else 'black'
            fmt = '.0f' if abs(val) > 10 else '.2f' if abs(val) > 1 else '.3f'
            ax.text(j, i, f'{val:{fmt}}', ha='center', va='center',
                    fontsize=5.5, color=text_color)

    # Mark Pareto designs
    pareto_ids = set(df_pareto.get('global_id', pd.Series())) if not df_pareto.empty else set()
    for i, (_, row) in enumerate(top50.iterrows()):
        gid = row.get('global_id', '')
        if gid in pareto_ids:
            ax.text(-0.7, i, 'P', ha='center', va='center',
                    fontsize=7, fontweight='bold', color='red')

    ax.set_title('Top 50 Global by Fitness Score (sorted by weight, green=better)\n'
                 'P=Pareto', fontsize=11)
    plt.colorbar(im, ax=ax, label='Normalized (1=best)', shrink=0.4)
    plt.tight_layout()
    plt.savefig(output_dir / "07_top50_global_heatmap.png",
                dpi=300, bbox_inches="tight")
    plt.close()

    top50[['global_id', 'template', 'seq_id'] + all_cols].to_csv(
        csv_dir / "07_top50_global_heatmap_data.csv", index=False)
    print("  07_top50_global_heatmap.png")


def create_visualizations(df_all: pd.DataFrame, df_selected: pd.DataFrame,
                          df_pareto: pd.DataFrame, output_dir: Path,
                          template_df: pd.DataFrame = None):
    """Generate summary plots with SD coloring and template references."""
    output_dir.mkdir(exist_ok=True, parents=True)
    csv_dir = output_dir / "plot_data_csvs"
    csv_dir.mkdir(exist_ok=True)

    # ---- 1. Selection summary (6-panel) ----
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))

    # (a) Fitness score distribution: all vs selected vs pareto
    ax = axes[0, 0]
    ax.hist(df_all["fitness_score"].dropna(), bins=40,
            alpha=0.35, label="All", color="gray")
    ax.hist(df_selected["fitness_score"].dropna(), bins=40,
            alpha=0.6, label="Selected", color="steelblue")
    if not df_pareto.empty and "fitness_score" in df_pareto.columns:
        ax.hist(df_pareto["fitness_score"].dropna(), bins=20,
                alpha=0.7, label="Pareto", color="red")
    ax.set_xlabel("Fitness Score")
    ax.set_ylabel("Count")
    ax.set_title("Fitness Score Distribution")
    ax.legend(fontsize=8)

    # (b) Fitness score violin per template
    ax = axes[0, 1]
    if "template" in df_all.columns:
        templates = sorted(df_all["template"].unique())
        tpl_short = [t.replace("output_G4FP_", "").replace("G4FP_", "")
                     for t in templates]
        data_all = [df_all[df_all["template"] == t]["fitness_score"].dropna().values
                    for t in templates]
        data_sel = [df_selected[df_selected["template"] == t]["fitness_score"].dropna().values
                    for t in templates]
        positions = np.arange(len(templates))
        # All designs (gray violins)
        parts = ax.violinplot([d for d in data_all if len(d) > 1],
                              positions=[p for p, d in zip(positions, data_all)
                                         if len(d) > 1],
                              showmedians=True, widths=0.7)
        for pc in parts.get('bodies', []):
            pc.set_facecolor('lightgray')
            pc.set_alpha(0.5)
        # Selected (colored box on top)
        sel_data = [(p, d) for p, d in zip(positions, data_sel) if len(d) > 0]
        if sel_data:
            bp = ax.boxplot([d for _, d in sel_data],
                            positions=[p for p, _ in sel_data],
                            widths=0.3, patch_artist=True,
                            showfliers=False, zorder=3)
            for patch in bp['boxes']:
                patch.set_facecolor('steelblue')
                patch.set_alpha(0.7)
        ax.set_xticks(positions)
        ax.set_xticklabels(tpl_short, rotation=45, ha='right', fontsize=7)
        ax.set_ylabel("Fitness Score")
        ax.set_title("Fitness by Template (gray=all, blue=selected)")

    # (c) Score component radar: selected median vs all median
    ax = axes[0, 2]
    comp_cols = [c for c in df_all.columns if c.startswith("score_component_")]
    if len(comp_cols) >= 3:
        comp_labels = [c.replace("score_component_", "").replace("_", "\n")
                       for c in comp_cols]
        all_medians = [df_all[c].median() for c in comp_cols]
        sel_medians = [df_selected[c].median() for c in comp_cols]
        n_comp = len(comp_cols)
        angles = np.linspace(0, 2 * np.pi, n_comp, endpoint=False).tolist()
        angles += angles[:1]
        all_medians += all_medians[:1]
        sel_medians += sel_medians[:1]

        ax.set_visible(False)
        ax_r = fig.add_subplot(2, 3, 3, polar=True)
        ax_r.plot(angles, all_medians, 'o-', color='gray', alpha=0.6,
                  label='All median', linewidth=1.5)
        ax_r.fill(angles, all_medians, alpha=0.1, color='gray')
        ax_r.plot(angles, sel_medians, 'o-', color='steelblue', alpha=0.8,
                  label='Selected median', linewidth=2)
        ax_r.fill(angles, sel_medians, alpha=0.15, color='steelblue')
        ax_r.set_xticks(angles[:-1])
        ax_r.set_xticklabels(comp_labels, fontsize=7)
        ax_r.set_ylim(0, 1)
        ax_r.set_title("Score Components\n(0-1 normalized)", fontsize=9, y=1.08)
        ax_r.legend(fontsize=7, loc='lower right')
    else:
        ax.text(0.5, 0.5, "Score components\nnot available",
                transform=ax.transAxes, ha='center', va='center')

    # (d) iPTM vs Chromophore pLDDT scatter colored by template
    ax = axes[1, 0]
    chrom_c = ("holo_chromophore_plddt" if "holo_chromophore_plddt" in df_all.columns
               else "chromophore_holo_plddt")
    if "holo_iptm" in df_all.columns and chrom_c in df_all.columns:
        cmap_t = plt.cm.get_cmap('tab10', max(10, len(templates)))
        for i, t in enumerate(templates):
            mask_all = df_all["template"] == t
            mask_sel = df_selected["template"] == t
            short = t.replace("output_G4FP_", "").replace("G4FP_", "")
            ax.scatter(df_all.loc[mask_all, "holo_iptm"],
                       df_all.loc[mask_all, chrom_c],
                       alpha=0.15, s=10, color=cmap_t(i))
            ax.scatter(df_selected.loc[mask_sel, "holo_iptm"],
                       df_selected.loc[mask_sel, chrom_c],
                       alpha=0.7, s=25, color=cmap_t(i),
                       edgecolors='black', linewidths=0.3, label=short)
        ax.set_xlabel("Holo iPTM (Interface)")
        ax.set_ylabel("Holo Chromophore pLDDT")
        ax.set_title("Key Metrics by Template (large=selected)")
        ax.legend(fontsize=6, ncol=2, loc='lower right')
    else:
        ax.text(0.5, 0.5, "iPTM / chromophore data\nnot available",
                transform=ax.transAxes, ha='center', va='center')

    # (e) Rank vs fitness with template coloring (top 50)
    ax = axes[1, 1]
    if "final_rank" in df_selected.columns:
        top50 = df_selected.nsmallest(50, "final_rank")
        cmap_t = plt.cm.get_cmap('tab10', max(10, len(templates)))
        tpl_to_i = {t: i for i, t in enumerate(templates)}
        colors = [cmap_t(tpl_to_i.get(t, 0)) for t in top50["template"]]
        ax.barh(top50["final_rank"], top50["fitness_score"],
                color=colors, edgecolor='white', linewidth=0.3)
        ax.set_ylabel("Rank")
        ax.set_xlabel("Fitness Score")
        ax.set_title("Top 50 Candidates by Rank")
        ax.invert_yaxis()
        # Template legend
        for t in templates:
            short = t.replace("output_G4FP_", "").replace("G4FP_", "")
            if t in top50["template"].values:
                ax.barh([], [], color=cmap_t(tpl_to_i[t]), label=short)
        ax.legend(fontsize=6, loc='lower right', ncol=2)

    # (f) Key metric distributions: selected boxplots
    ax = axes[1, 2]
    box_metrics = [
        ("holo_iptm", "iPTM"),
        (chrom_c, "Chrom\npLDDT"),
        ("global_rmsd", "RMSD"),
        ("holo_ptm", "pTM"),
        ("mean_plddt_diff", "pLDDT\ndiff"),
    ]
    box_data = [(l, df_selected[c].dropna().values)
                for c, l in box_metrics
                if c in df_selected.columns and df_selected[c].notna().any()]
    if box_data:
        # Normalize each to 0-1 for comparison
        norm_data = []
        raw_ranges = []
        for label, vals in box_data:
            vmin, vmax = vals.min(), vals.max()
            raw_ranges.append((vmin, vmax))
            if vmax > vmin:
                norm_data.append((label, (vals - vmin) / (vmax - vmin)))
            else:
                norm_data.append((label, np.full_like(vals, 0.5)))
        bp = ax.boxplot([v for _, v in norm_data], patch_artist=True,
                        showfliers=True, flierprops=dict(markersize=2))
        colors_box = ['#E07B39', '#D62728', '#2CA02C', '#FF7F0E', '#9467BD']
        for i, patch in enumerate(bp['boxes']):
            patch.set_facecolor(colors_box[i % len(colors_box)])
            patch.set_alpha(0.7)
        ax.set_xticklabels([l for l, _ in norm_data], fontsize=8)
        # Add raw range as annotation
        for i, ((vmin, vmax), (label, _)) in enumerate(zip(raw_ranges, norm_data)):
            ax.text(i + 1, -0.08, f"[{vmin:.1f}-{vmax:.1f}]",
                    ha='center', fontsize=6, color='gray',
                    transform=ax.get_xaxis_transform())
        ax.set_ylabel("Normalized Value")
        ax.set_title("Selected: Key Metric Distributions")
        ax.set_ylim(-0.05, 1.15)

    plt.suptitle(f"Selection Summary: {len(df_selected)} designs from "
                 f"{df_selected['template'].nunique()} templates",
                 fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig(output_dir / "07_selection_summary.png", dpi=300,
                bbox_inches="tight")
    plt.close()
    print("  07_selection_summary.png")

    # ---- 2. Holo vs Apo scatter (fitness-colored) ----
    has_holo_apo = ("holo_mean_plddt" in df_all.columns and
                     "apo_mean_plddt" in df_all.columns)
    if has_holo_apo:
        fig, ax = plt.subplots(figsize=(10, 8))

        ax.scatter(df_all["holo_mean_plddt"], df_all["apo_mean_plddt"],
                   alpha=0.15, s=15, color="gray", label="All designs")

        sc = ax.scatter(df_selected["holo_mean_plddt"],
                        df_selected["apo_mean_plddt"],
                        c=df_selected["fitness_score"], cmap="viridis",
                        alpha=0.7, s=30, edgecolors="black", linewidths=0.3,
                        label="Selected")
        plt.colorbar(sc, ax=ax, label="Fitness Score")

        if not df_pareto.empty and "holo_mean_plddt" in df_pareto.columns:
            pareto_sorted = df_pareto.sort_values("holo_mean_plddt")
            ax.plot(pareto_sorted["holo_mean_plddt"],
                    pareto_sorted["apo_mean_plddt"],
                    "r-", alpha=0.5, linewidth=1.5, label="Pareto frontier")
            ax.scatter(df_pareto["holo_mean_plddt"],
                       df_pareto["apo_mean_plddt"],
                       c="red", s=40, marker="D", alpha=0.6, zorder=5)

        lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lims, lims, "k--", alpha=0.3, label="y=x")

        _plot_template_references(
            ax, template_df,
            'holo_mean_plddt', 'apo_mean_plddt',
            'holo_mean_plddt_sd', 'apo_mean_plddt_sd')

        ax.annotate("IDEAL\n(high holo,\nlow apo)",
                     xy=(85, 35), fontsize=10, color="green",
                     ha="center", fontstyle="italic",
                     bbox=dict(boxstyle="round,pad=0.3",
                               fc="lightgreen", alpha=0.3))

        ax.set_xlabel("Holo Mean pLDDT")
        ax.set_ylabel("Apo Mean pLDDT")
        ax.set_title("G4FP Design Space: Holo vs Apo (colored by fitness)")
        ax.legend(loc="upper left", fontsize=8)
        plt.tight_layout()
        plt.savefig(output_dir / "07_holo_vs_apo_scatter.png", dpi=300,
                    bbox_inches="tight")
        plt.close()
        print("  07_holo_vs_apo_scatter.png")

    # ---- 3. Holo vs Apo scatter (SD-colored) ----
    has_sd = "mean_plddt_diff_sd" in df_all.columns
    if has_holo_apo and has_sd:
        fig, ax = plt.subplots(figsize=(10, 8))

        # Compute combined SD for color
        holo_sd = df_all.get("holo_mean_plddt_sd",
                               pd.Series(0, index=df_all.index))
        apo_sd = df_all.get("apo_mean_plddt_sd",
                             pd.Series(0, index=df_all.index))
        combined_sd = np.sqrt(holo_sd.fillna(0)**2 + apo_sd.fillna(0)**2)

        sc = ax.scatter(df_all["holo_mean_plddt"], df_all["apo_mean_plddt"],
                        c=combined_sd, cmap="plasma",
                        alpha=0.5, s=20, edgecolors="none",
                        vmin=0, vmax=max(combined_sd.quantile(0.95), 1))
        plt.colorbar(sc, ax=ax,
                     label="Combined SD (holo + apo pLDDT)")

        # Highlight selected with black edge
        if not df_selected.empty:
            sel_sd_b = df_selected.get("holo_mean_plddt_sd",
                                        pd.Series(0, index=df_selected.index))
            sel_sd_a = df_selected.get("apo_mean_plddt_sd",
                                        pd.Series(0, index=df_selected.index))
            sel_csd = np.sqrt(sel_sd_b.fillna(0)**2 + sel_sd_a.fillna(0)**2)
            ax.scatter(df_selected["holo_mean_plddt"],
                       df_selected["apo_mean_plddt"],
                       c=sel_csd, cmap="plasma",
                       alpha=0.8, s=35, edgecolors="black", linewidths=0.5,
                       vmin=0, vmax=max(combined_sd.quantile(0.95), 1))

        lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lims, lims, "k--", alpha=0.3, label="y=x")

        _plot_template_references(
            ax, template_df,
            'holo_mean_plddt', 'apo_mean_plddt',
            'holo_mean_plddt_sd', 'apo_mean_plddt_sd')

        ax.set_xlabel("Holo Mean pLDDT")
        ax.set_ylabel("Apo Mean pLDDT")
        ax.set_title("G4FP Design Space: Holo vs Apo "
                      "(colored by prediction SD)")
        ax.legend(loc="upper left", fontsize=8)
        plt.tight_layout()
        plt.savefig(output_dir / "07_holo_vs_apo_sd.png", dpi=300,
                    bbox_inches="tight")
        plt.close()
        print("  07_holo_vs_apo_sd.png")

    # ---- 4. SD overview ----
    sd_cols = [c for c in df_all.columns if c.endswith("_sd") and
               df_all[c].notna().sum() > 0]
    if sd_cols:
        n_sd = len(sd_cols)
        n_cols = min(3, n_sd)
        n_rows = (n_sd + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols,
                                  figsize=(5 * n_cols, 4 * n_rows))
        axes = np.atleast_2d(axes)

        for idx, col in enumerate(sorted(sd_cols)):
            r, c = divmod(idx, n_cols)
            ax = axes[r, c]
            vals = df_all[col].dropna()
            if len(vals) == 0:
                ax.set_visible(False)
                continue

            # Color by holo/apo/difference
            if col.startswith("holo_") or col.startswith("chromophore_holo"):
                bar_color = "#E07B39"  # orange = holo
            elif col.startswith("apo_") or col.startswith("chromophore_apo"):
                bar_color = "#4878CF"  # blue = apo
            else:
                bar_color = "#7B68AE"  # purple = difference/combined

            ax.hist(vals, bins=30, alpha=0.6, color=bar_color,
                    edgecolor="white")
            # Selected overlay
            if col in df_selected.columns:
                sel_vals = df_selected[col].dropna()
                if len(sel_vals) > 0:
                    ax.hist(sel_vals, bins=30, alpha=0.5, color="#2CA02C",
                            edgecolor="white", label="Selected")

            med = vals.median()
            ax.axvline(med, color="red", linestyle="--", alpha=0.7,
                        label=f"median={med:.2f}")

            label = col.replace("_sd", " SD").replace("_", " ").title()
            ax.set_xlabel(label)
            ax.set_ylabel("Count")
            ax.set_title(label)
            ax.legend(fontsize=7)

        # Hide unused axes
        for idx in range(n_sd, n_rows * n_cols):
            r, c = divmod(idx, n_cols)
            axes[r, c].set_visible(False)

        plt.suptitle("Standard Deviation Distributions (across seeds)",
                      fontsize=13, y=1.01)
        plt.tight_layout()
        plt.savefig(output_dir / "07_sd_overview.png", dpi=300,
                    bbox_inches="tight")
        plt.close()
        print("  07_sd_overview.png")

    # ---- 5. Pareto frontier plots ----
    if not df_pareto.empty:
        chrom_col = ("holo_chromophore_plddt" if "holo_chromophore_plddt" in df_all.columns
                     else "chromophore_holo_plddt")
        plot_pairs = [
            ("holo_iptm", chrom_col,
             "Holo iPTM (Interface)", "Holo Chromophore pLDDT",
             "Pareto: Interface Quality vs Chromophore Stability"),
            ("holo_iptm", "global_rmsd",
             "Holo iPTM (Interface)", "Global RMSD (A)",
             "Pareto: Interface Quality vs Structural Change"),
        ]

        fig, axes = plt.subplots(1, 2, figsize=(18, 7))
        for ax_idx, (obj1, obj2, xlabel, ylabel, title) in enumerate(plot_pairs):
            ax = axes[ax_idx]
            if obj1 not in df_all.columns or obj2 not in df_all.columns:
                ax.set_visible(False)
                continue

            sc = ax.scatter(df_all[obj1], df_all[obj2],
                           c=df_all["fitness_score"], cmap="viridis",
                           alpha=0.3, s=15, edgecolors="none")
            plt.colorbar(sc, ax=ax, label="Fitness Score")

            # Draw Pareto frontier as a connected line (step envelope)
            if obj1 in df_pareto.columns and obj2 in df_pareto.columns:
                pf = df_pareto[[obj1, obj2]].dropna().copy()
                # Sort by obj1 ascending, then trace the upper-right envelope
                pf = pf.sort_values(obj1)
                # Build step-line: at each x, the frontier y is the max y
                # among Pareto points with x >= current x
                xs = pf[obj1].values
                ys = pf[obj2].values
                # Running max from right gives the envelope
                envelope_y = np.maximum.accumulate(ys[::-1])[::-1]
                ax.step(xs, envelope_y, where='post', color='red',
                        linewidth=2, alpha=0.8, label="Pareto frontier")
                ax.scatter(xs, ys, c="red", s=50, marker="D",
                           alpha=0.7, zorder=5)

            # Highlight overlap (selected) designs
            if not df_selected.empty and obj1 in df_selected.columns:
                ax.scatter(df_selected[obj1], df_selected[obj2],
                           alpha=0.7, s=35, color="steelblue",
                           edgecolors="black", linewidths=0.5,
                           label=f"Selected ({len(df_selected)})", zorder=4)

            _plot_template_references(ax, template_df, obj1, obj2)

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            ax.legend(fontsize=7)

        plt.tight_layout()
        plt.savefig(output_dir / "07_pareto_plot.png", dpi=300,
                    bbox_inches="tight")
        plt.close()
        print("  07_pareto_plot.png")

    # ---- 6. Fitness component histograms ----
    # Match the new fitness score components
    chrom_col = ("holo_chromophore_plddt" if "holo_chromophore_plddt" in df_all.columns
                 else "chromophore_holo_plddt")
    components = [
        (chrom_col, "Holo Chromophore pLDDT (25%)", "#E07B39"),
        ("holo_iptm", "Holo iPTM - Interface (25%)", "#D62728"),
        ("global_rmsd", "Global RMSD (15%)", "#2CA02C"),
        ("holo_ptm", "Holo pTM (15%)", "#FF7F0E"),
        ("apo_plddt_vs_template", "Apo pLDDT vs Template (10%)", "#4878CF"),
        ("chromophore_plddt_diff", "Chromophore pLDDT Diff (5%)", "#9467BD"),
        ("mean_plddt_diff_sd", "pLDDT Diff SD (5%, lower=better)", "#7B68AE"),
        ("holo_fraction_disordered", "Holo Fraction Disordered", "#8C564B"),
        ("apo_fraction_disordered", "Apo Fraction Disordered", "#17BECF"),
    ]
    avail_comp = [(c, l, clr) for c, l, clr in components if c in df_all.columns]
    if avail_comp:
        n_comp = len(avail_comp)
        ncols = min(3, n_comp)
        nrows = (n_comp + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4.5 * nrows))
        axes = np.atleast_2d(axes).flatten()

        for idx, (col, label, color) in enumerate(avail_comp):
            ax = axes[idx]
            vals_all = df_all[col].dropna()
            vals_sel = df_selected[col].dropna() if col in df_selected.columns else pd.Series()
            vals_par = df_pareto[col].dropna() if col in df_pareto.columns else pd.Series()

            ax.hist(vals_all, bins=30, alpha=0.3, color="gray", label="All")
            if len(vals_sel) > 0:
                ax.hist(vals_sel, bins=30, alpha=0.6, color=color, label="Selected")
            if len(vals_par) > 0:
                for v in vals_par:
                    ax.axvline(v, color="red", alpha=0.15, linewidth=0.5)
                ax.axvline(vals_par.median(), color="red", linestyle="--",
                           alpha=0.7, label=f"Pareto median={vals_par.median():.2f}")

            ax.axvline(vals_all.median(), color="black", linestyle=":",
                       alpha=0.5, label=f"All median={vals_all.median():.2f}")
            ax.set_xlabel(label)
            ax.set_ylabel("Count")
            ax.set_title(label)
            ax.legend(fontsize=7)

        for idx in range(len(avail_comp), len(axes)):
            axes[idx].set_visible(False)

        plt.suptitle("Fitness Score Components: All vs Selected vs Pareto",
                     fontsize=13, y=1.01)
        plt.tight_layout()
        plt.savefig(output_dir / "07_fitness_components.png", dpi=300,
                    bbox_inches="tight")
        plt.close()
        print("  07_fitness_components.png")

    # ---- 7. Top 10 per template heatmap ----
    _plot_top10_per_template_heatmap(df_all, df_selected, df_pareto,
                                     output_dir, csv_dir)

    # ---- 8. Top 50 global heatmap ----
    _plot_top50_global_heatmap(df_all, df_selected, df_pareto,
                               output_dir, csv_dir)

    # ---- Export CSVs ----
    df_all.to_csv(csv_dir / "07_all_designs_scored.csv", index=False)
    df_selected.to_csv(csv_dir / "07_selected_designs.csv", index=False)
    if not df_pareto.empty:
        df_pareto.to_csv(csv_dir / "07_pareto_frontier.csv", index=False)
    print(f"  CSVs exported to {csv_dir}/")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Aggregate G4FP designs across templates and rank",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Selection strategy:
  1. Load 06_ligand_state_comparison_results.csv from all output directories
  2. Compute fitness score (weighted components + confidence penalty for high SD)
  3. Apply hard filters (holo must fold)
  4. Find Pareto-optimal designs across iPTM, chromophore pLDDT, RMSD
  5. Take top N by fitness score
  6. Output the OVERLAP of Pareto-optimal AND top-N fitness
        """
    )
    parser.add_argument("--output-dir", type=str,
                        help="Analyze a single output directory instead of all")
    parser.add_argument("--min-holo-plddt", type=float, default=70.0,
                        help="Minimum holo state pLDDT (default: 70)")
    parser.add_argument("--min-holo-ptm", type=float, default=0.5,
                        help="Minimum holo state pTM (default: 0.5)")
    parser.add_argument("--min-holo-iptm", type=float, default=0.3,
                        help="Minimum holo iPTM (default: 0.3)")
    parser.add_argument("--max-apo-plddt", type=float, default=None,
                        help="Maximum apo pLDDT (default: None, use relative)")

    args = parser.parse_args()
    base_dir = Path(__file__).resolve().parent
    results_dir = base_dir / "07_results"

    print("=" * 80)
    print("STEP 7: Aggregate and Rank G4FP Designs")
    print("=" * 80)

    # 1. Load data
    print("\nLoading comparison results...")
    if args.output_dir:
        output_dirs = [Path(args.output_dir)]
    else:
        output_dirs = None
    df = aggregate_all(base_dir, output_dirs)

    # Report SD availability
    sd_cols = [c for c in df.columns if c.endswith("_sd")]
    if sd_cols:
        print(f"\n  SD columns found: {len(sd_cols)}")
        for sc in sorted(sd_cols)[:8]:
            med = df[sc].median()
            print(f"    {sc}: median={med:.2f}")
        if len(sd_cols) > 8:
            print(f"    ... and {len(sd_cols) - 8} more")

    # Load template references early (needed for relative apo pLDDT)
    print("\nLoading template references...")
    template_df = load_template_metrics(base_dir)

    # Compute relative apo pLDDT (vs each design's parent template)
    df = compute_relative_apo_plddt(df, template_df)
    if "apo_plddt_vs_template" in df.columns:
        valid = df["apo_plddt_vs_template"].dropna()
        if len(valid) > 0:
            print(f"  Relative apo pLDDT: {valid.min():.1f} to {valid.max():.1f} "
                  f"(median={valid.median():.1f})")

    # 2. Compute fitness score (includes confidence penalty)
    print("\nComputing fitness scores...")
    df = compute_fitness_score(df)
    df = df.sort_values("fitness_score", ascending=False)
    print(f"  Score range: {df['fitness_score'].min():.4f} - "
          f"{df['fitness_score'].max():.4f}")
    print(f"  Median: {df['fitness_score'].median():.4f}")

    # 3. Pareto frontier
    print("\nFinding Pareto-optimal designs...")
    # 3 objectives for meaningful Pareto frontier (more -> all non-dominated)
    chrom_pareto_col = ("holo_chromophore_plddt" if "holo_chromophore_plddt" in df.columns
                        else "chromophore_holo_plddt")
    pareto_objectives = [
        ("holo_iptm", True),             # interface quality
        (chrom_pareto_col, True),         # chromophore stability
        ("global_rmsd", True),            # structural change
    ]
    df_pareto = pareto_frontier(df, pareto_objectives)

    # 4. Hard filters
    df_filtered = apply_hard_filters(
        df,
        min_holo_plddt=args.min_holo_plddt,
        min_holo_ptm=args.min_holo_ptm,
        min_holo_iptm=args.min_holo_iptm,
        max_apo_plddt=args.max_apo_plddt,
    )

    # 5. Build two selected pools
    pareto_ids = set(df_pareto["global_id"])

    # Pool A: Top 10 per template
    top10_frames = []
    for tpl in sorted(df_filtered["template"].unique()):
        tpl_df = df_filtered[df_filtered["template"] == tpl].nlargest(
            10, "fitness_score")
        top10_frames.append(tpl_df)
    df_top10_per_tpl = pd.concat(top10_frames).sort_values(
        "fitness_score", ascending=False).reset_index(drop=True)
    df_top10_per_tpl["rank_in_pool"] = range(1, len(df_top10_per_tpl) + 1)

    # Pool B: Top 50 global by fitness
    df_top50_global = df_filtered.nlargest(50, "fitness_score").copy()
    df_top50_global = df_top50_global.reset_index(drop=True)
    df_top50_global["rank_in_pool"] = range(1, len(df_top50_global) + 1)

    # Combined selected = union of both pools (for summary plots)
    combined_ids = set(df_top10_per_tpl["global_id"]) | set(df_top50_global["global_id"])
    df_selected = df_filtered[df_filtered["global_id"].isin(combined_ids)].copy()
    df_selected = df_selected.sort_values("fitness_score", ascending=False).reset_index(drop=True)
    df_selected["final_rank"] = range(1, len(df_selected) + 1)

    # Tag designs in full df
    df["in_pareto"] = df["global_id"].isin(pareto_ids)
    df["in_top10_per_tpl"] = df["global_id"].isin(set(df_top10_per_tpl["global_id"]))
    df["in_top50_global"] = df["global_id"].isin(set(df_top50_global["global_id"]))

    n_overlap = len(set(df_top10_per_tpl["global_id"]) & set(df_top50_global["global_id"]))
    print(f"\n  Pool A (top 10/template): {len(df_top10_per_tpl)} designs "
          f"from {df_top10_per_tpl['template'].nunique()} templates")
    print(f"  Pool B (top 50 global):   {len(df_top50_global)} designs "
          f"from {df_top50_global['template'].nunique()} templates")
    print(f"  Overlap A & B:            {n_overlap}")
    print(f"  Union (total selected):   {len(df_selected)}")
    print(f"  Pareto-optimal:           {len(df_pareto)}")

    # 6. Save outputs
    results_dir.mkdir(exist_ok=True, parents=True)

    def _save_pool(pool_df, name, results_dir):
        """Save CSV and FASTA for a selection pool."""
        export_cols = [c for c in pool_df.columns if "array" not in c]
        pool_df[export_cols].to_csv(results_dir / f"{name}.csv", index=False)
        if "sequence" in pool_df.columns:
            with open(results_dir / f"{name}.fa", "w") as f:
                for _, row in pool_df.iterrows():
                    rank = int(row.get('rank_in_pool', row.get('final_rank', 0)))
                    header = (f">{row['global_id']} "
                              f"rank={rank} "
                              f"fitness={row['fitness_score']:.4f} "
                              f"template={row['template']}")
                    if pd.notna(row.get("mean_plddt_diff")):
                        header += f" delta_plddt={row['mean_plddt_diff']:.1f}"
                    f.write(f"{header}\n{row['sequence']}\n")
        print(f"  {name}.csv ({len(pool_df)} designs)")

    _save_pool(df_top10_per_tpl, "07_top10_per_template", results_dir)
    _save_pool(df_top50_global, "07_top50_global", results_dir)
    _save_pool(df_selected, "07_final_candidates", results_dir)

    if not df_pareto.empty:
        pareto_export = [c for c in df_pareto.columns if "array" not in c]
        df_pareto[pareto_export].to_csv(
            results_dir / "07_pareto_frontier.csv", index=False)
        print(f"  07_pareto_frontier.csv ({len(df_pareto)} designs)")

    # Visualizations
    print("\nGenerating plots...")
    create_visualizations(df, df_selected, df_pareto, results_dir, template_df)

    # Summary
    print("\n" + "=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)
    print(f"  Total designs analyzed:  {len(df)}")
    print(f"  Pareto-optimal designs:  {len(df_pareto)}")
    print(f"  Top 10 per template:     {len(df_top10_per_tpl)}")
    print(f"  Top 50 global:           {len(df_top50_global)}")
    print(f"  Combined selected:       {len(df_selected)}")
    print(f"  Templates represented:   {df_selected['template'].nunique()}")

    print(f"\n  Output directory: {results_dir}/")
    print(f"    07_top10_per_template.csv  -- top 10 per template + FASTA")
    print(f"    07_top50_global.csv        -- top 50 by fitness + FASTA")
    print(f"    07_final_candidates.csv    -- union of both pools")
    print(f"    07_pareto_frontier.csv     -- Pareto-optimal set")
    print(f"    07_selection_summary.png   -- overview plots")
    print(f"    07_top10_per_template_heatmap.png")
    print(f"\nNext step: python 08_export_for_synthesis.py")


if __name__ == "__main__":
    main()

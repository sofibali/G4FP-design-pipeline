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

def compute_fitness_score(df: pd.DataFrame,
                          w_plddt_diff: float = 0.30,
                          w_holo_ptm: float = 0.20,
                          w_apo_disorder: float = 0.20,
                          w_chrom_plddt_diff: float = 0.15,
                          w_rmsd: float = 0.10,
                          w_confidence: float = 0.05) -> pd.DataFrame:
    """
    Compute a composite G4FP fitness score.

    Higher score = better G4-dependent folding switch.

    Components (all normalized 0-1 before weighting):
      - delta_plddt:  mean_plddt_diff (holo - apo), higher = better
      - holo_ptm:    holo state pTM, higher = better fold
      - apo_disorder: (100 - apo_mean_plddt)/100, higher = more disordered apo
      - chrom_diff:   chromophore_plddt_diff, higher = better local switch
      - rmsd:         global_rmsd between states, higher = more structural change
      - confidence:   1 - normalized(mean_plddt_diff_sd), penalizes high variability
    """
    df = df.copy()

    def _norm(series):
        """Min-max normalize to 0-1."""
        smin, smax = series.min(), series.max()
        if smax == smin:
            return pd.Series(0.5, index=series.index)
        return (series - smin) / (smax - smin)

    # Delta pLDDT (higher = bigger switch)
    if "mean_plddt_diff" in df.columns:
        df["_n_plddt_diff"] = _norm(df["mean_plddt_diff"].fillna(0))
    else:
        df["_n_plddt_diff"] = 0

    # Holo pTM (higher = better fold when holo)
    if "holo_ptm" in df.columns:
        df["_n_holo_ptm"] = _norm(df["holo_ptm"].fillna(0))
    else:
        df["_n_holo_ptm"] = 0

    # Apo disorder (higher = more disordered without ligand)
    if "apo_mean_plddt" in df.columns:
        df["_n_apo_disorder"] = _norm(100 - df["apo_mean_plddt"].fillna(100))
    else:
        df["_n_apo_disorder"] = 0

    # Chromophore pLDDT difference
    if "chromophore_plddt_diff" in df.columns:
        df["_n_chrom_diff"] = _norm(df["chromophore_plddt_diff"].fillna(0))
    else:
        df["_n_chrom_diff"] = 0

    # RMSD between states
    if "global_rmsd" in df.columns:
        df["_n_rmsd"] = _norm(df["global_rmsd"].fillna(0))
    else:
        df["_n_rmsd"] = 0

    # Confidence penalty: lower SD across seeds = more reproducible prediction
    if "mean_plddt_diff_sd" in df.columns:
        sd_vals = df["mean_plddt_diff_sd"].fillna(df["mean_plddt_diff_sd"].median())
        df["_n_confidence"] = 1.0 - _norm(sd_vals)
    else:
        df["_n_confidence"] = 0.5

    df["fitness_score"] = (
        w_plddt_diff * df["_n_plddt_diff"]
        + w_holo_ptm * df["_n_holo_ptm"]
        + w_apo_disorder * df["_n_apo_disorder"]
        + w_chrom_plddt_diff * df["_n_chrom_diff"]
        + w_rmsd * df["_n_rmsd"]
        + w_confidence * df["_n_confidence"]
    )

    df.drop(columns=[c for c in df.columns if c.startswith("_n_")], inplace=True)
    return df


# ---------------------------------------------------------------------------
# Hard filters
# ---------------------------------------------------------------------------

def apply_hard_filters(df: pd.DataFrame,
                       min_holo_plddt: float = 70.0,
                       min_holo_ptm: float = 0.5,
                       max_apo_plddt: float = 60.0) -> pd.DataFrame:
    """Apply hard cutoffs. Designs that fail any filter are excluded."""
    n_before = len(df)
    mask = pd.Series(True, index=df.index)

    if "holo_mean_plddt" in df.columns:
        mask &= df["holo_mean_plddt"] >= min_holo_plddt
    if "holo_ptm" in df.columns:
        mask &= df["holo_ptm"] >= min_holo_ptm
    if "apo_mean_plddt" in df.columns:
        mask &= df["apo_mean_plddt"] <= max_apo_plddt

    df_filtered = df[mask].copy()
    n_after = len(df_filtered)
    print(f"\nHard filters: {n_before} -> {n_after} designs "
          f"(removed {n_before - n_after})")
    print(f"  holo_plddt >= {min_holo_plddt}, holo_ptm >= {min_holo_ptm}, "
          f"apo_plddt <= {max_apo_plddt}")

    if n_after == 0:
        print("\nWARNING: All designs filtered out! Relaxing filters to keep top 50%.")
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


def create_visualizations(df_all: pd.DataFrame, df_selected: pd.DataFrame,
                          df_pareto: pd.DataFrame, output_dir: Path,
                          template_df: pd.DataFrame = None):
    """Generate summary plots with SD coloring and template references."""
    output_dir.mkdir(exist_ok=True, parents=True)
    csv_dir = output_dir / "plot_data_csvs"
    csv_dir.mkdir(exist_ok=True)

    # ---- 1. Selection summary ----
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    axes[0, 0].hist(df_all["fitness_score"].dropna(), bins=40,
                     alpha=0.5, label="All designs", color="gray")
    axes[0, 0].hist(df_selected["fitness_score"].dropna(), bins=40,
                     alpha=0.7, label="Selected", color="steelblue")
    axes[0, 0].set_xlabel("Fitness Score")
    axes[0, 0].set_ylabel("Count")
    axes[0, 0].set_title("Fitness Score Distribution")
    axes[0, 0].legend()

    template_counts = df_selected["template"].value_counts()
    axes[0, 1].barh(range(len(template_counts)), template_counts.values,
                     color="steelblue")
    axes[0, 1].set_yticks(range(len(template_counts)))
    labels = [t.replace("output_G4FP_", "") for t in template_counts.index]
    axes[0, 1].set_yticklabels(labels, fontsize=8)
    axes[0, 1].set_xlabel("Selected Designs")
    axes[0, 1].set_title("Selected Designs per Template")

    if "mean_plddt_diff" in df_selected.columns:
        axes[1, 0].hist(df_selected["mean_plddt_diff"].dropna(), bins=30,
                         alpha=0.7, color="green")
        axes[1, 0].set_xlabel("pLDDT Difference (Holo - Apo)")
        axes[1, 0].set_ylabel("Count")
        axes[1, 0].set_title("Selected: pLDDT Switch Magnitude")
        axes[1, 0].axvline(0, color="red", linestyle="--", alpha=0.5)
        # SD shading if available
        if "mean_plddt_diff_sd" in df_selected.columns:
            med_sd = df_selected["mean_plddt_diff_sd"].median()
            axes[1, 0].axvline(0 + med_sd, color="orange", linestyle=":",
                                alpha=0.4, label=f"median SD={med_sd:.1f}")
            axes[1, 0].axvline(0 - med_sd, color="orange", linestyle=":",
                                alpha=0.4)
            axes[1, 0].legend(fontsize=8)

    if "cluster" in df_selected.columns:
        cluster_sizes = df_selected["cluster"].value_counts().sort_index()
        axes[1, 1].bar(range(len(cluster_sizes)), cluster_sizes.values,
                        alpha=0.7, color="purple")
        axes[1, 1].set_xlabel("Cluster ID")
        axes[1, 1].set_ylabel("Designs Selected")
        axes[1, 1].set_title("Diversity: Designs per Sequence Cluster")

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

    # ---- 5. Pareto frontier plot ----
    if not df_pareto.empty:
        obj1 = "mean_plddt_diff"
        obj2 = "holo_ptm"

        if obj1 in df_all.columns and obj2 in df_all.columns:
            fig, ax = plt.subplots(figsize=(10, 8))

            # Color all by SD if available
            if "mean_plddt_diff_sd" in df_all.columns:
                sc = ax.scatter(
                    df_all[obj1], df_all[obj2],
                    c=df_all["mean_plddt_diff_sd"].fillna(0),
                    cmap="plasma", alpha=0.3, s=15, edgecolors="none",
                    vmin=0,
                    vmax=max(
                        df_all["mean_plddt_diff_sd"].quantile(0.95), 1))
                plt.colorbar(sc, ax=ax, label="pLDDT Diff SD")
            else:
                ax.scatter(df_all[obj1], df_all[obj2],
                           alpha=0.15, s=15, color="gray", label="All")

            ax.scatter(df_pareto[obj1], df_pareto[obj2],
                       c="red", s=50, marker="D", alpha=0.7,
                       label="Pareto optimal")
            ax.scatter(df_selected[obj1], df_selected[obj2],
                       alpha=0.4, s=20, color="steelblue",
                       edgecolors="black", linewidths=0.2,
                       label="Selected")

            ax.set_xlabel("pLDDT Difference (Holo - Apo)")
            ax.set_ylabel("Holo pTM")
            ax.set_title("Pareto Frontier: pLDDT Switch vs "
                          "Holo Fold Quality")
            ax.legend()
            plt.tight_layout()
            plt.savefig(output_dir / "07_pareto_plot.png", dpi=300,
                        bbox_inches="tight")
            plt.close()
            print("  07_pareto_plot.png")

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
  3. Apply hard filters (holo must fold, apo should be disordered)
  4. Find Pareto-optimal designs across multiple objectives
  5. Diversity-weighted selection to avoid redundant sequences
  6. Export top N candidates for experimental testing
        """
    )
    parser.add_argument("--output-dir", type=str,
                        help="Analyze a single output directory instead of all")
    parser.add_argument("--n-select", type=int, default=500,
                        help="Number of final candidates (default: 500)")
    parser.add_argument("--min-holo-plddt", type=float, default=70.0,
                        help="Minimum holo state pLDDT (default: 70)")
    parser.add_argument("--min-holo-ptm", type=float, default=0.5,
                        help="Minimum holo state pTM (default: 0.5)")
    parser.add_argument("--max-apo-plddt", type=float, default=60.0,
                        help="Maximum apo state pLDDT (default: 60)")
    parser.add_argument("--no-diversity", action="store_true",
                        help="Skip diversity-weighted selection")
    parser.add_argument("--pareto-only", action="store_true",
                        help="Only output Pareto-optimal designs")

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

    # 2. Compute fitness score (includes confidence penalty)
    print("\nComputing fitness scores...")
    df = compute_fitness_score(df)
    df = df.sort_values("fitness_score", ascending=False)
    print(f"  Score range: {df['fitness_score'].min():.4f} - "
          f"{df['fitness_score'].max():.4f}")
    print(f"  Median: {df['fitness_score'].median():.4f}")

    # 3. Pareto frontier
    print("\nFinding Pareto-optimal designs...")
    pareto_objectives = [
        ("mean_plddt_diff", True),
        ("holo_ptm", True),
        ("apo_mean_plddt", False),
        ("chromophore_plddt_diff", True),
        ("global_rmsd", True),
    ]
    df_pareto = pareto_frontier(df, pareto_objectives)

    if args.pareto_only:
        df_selected = df_pareto.sort_values("fitness_score", ascending=False)
    else:
        # 4. Hard filters
        df_filtered = apply_hard_filters(
            df,
            min_holo_plddt=args.min_holo_plddt,
            min_holo_ptm=args.min_holo_ptm,
            max_apo_plddt=args.max_apo_plddt,
        )

        # 5. Diversity-weighted selection
        if args.no_diversity or "sequence" not in df_filtered.columns:
            df_selected = df_filtered.nlargest(
                args.n_select, "fitness_score").copy()
            print(f"\nSelected top {len(df_selected)} by fitness score")
        else:
            print(f"\nDiversity-weighted selection (target: {args.n_select})...")
            df_selected = diversity_weighted_selection(
                df_filtered, n_select=args.n_select)

    # 6. Final ranking
    df_selected = df_selected.sort_values(
        "fitness_score", ascending=False).reset_index(drop=True)
    df_selected["final_rank"] = range(1, len(df_selected) + 1)

    # 7. Save outputs
    results_dir.mkdir(exist_ok=True, parents=True)

    export_cols = [c for c in df_selected.columns if "array" not in c]
    df_selected[export_cols].to_csv(
        results_dir / "07_final_candidates.csv", index=False)

    if "sequence" in df_selected.columns:
        with open(results_dir / "07_final_candidates.fa", "w") as f:
            for _, row in df_selected.iterrows():
                header = (f">{row['global_id']} "
                          f"rank={row['final_rank']} "
                          f"fitness={row['fitness_score']:.4f} "
                          f"template={row['template']}")
                if "mean_plddt_diff" in row:
                    header += f" delta_plddt={row['mean_plddt_diff']:.1f}"
                if "mean_plddt_diff_sd" in row and pd.notna(
                        row.get("mean_plddt_diff_sd")):
                    header += f" sd={row['mean_plddt_diff_sd']:.1f}"
                f.write(f"{header}\n{row['sequence']}\n")

    if not df_pareto.empty:
        pareto_export = [c for c in df_pareto.columns if "array" not in c]
        df_pareto[pareto_export].to_csv(
            results_dir / "07_pareto_frontier.csv", index=False)

    # Load template references
    print("\nLoading template references...")
    template_df = load_template_metrics(base_dir)

    # Visualizations
    print("\nGenerating plots...")
    create_visualizations(df, df_selected, df_pareto, results_dir, template_df)

    # Summary
    print("\n" + "=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)
    print(f"  Total designs analyzed:  {len(df)}")
    print(f"  Pareto-optimal designs:  {len(df_pareto)}")
    print(f"  Final selected designs:  {len(df_selected)}")
    print(f"  Templates represented:   {df_selected['template'].nunique()}")
    print(f"\n  Score range (selected):  "
          f"{df_selected['fitness_score'].min():.4f} - "
          f"{df_selected['fitness_score'].max():.4f}")

    if "mean_plddt_diff" in df_selected.columns:
        print(f"  pLDDT diff range:        "
              f"{df_selected['mean_plddt_diff'].min():.1f} - "
              f"{df_selected['mean_plddt_diff'].max():.1f}")
        if "mean_plddt_diff_sd" in df_selected.columns:
            print(f"  pLDDT diff SD range:     "
                  f"{df_selected['mean_plddt_diff_sd'].min():.1f} - "
                  f"{df_selected['mean_plddt_diff_sd'].max():.1f}")

    print(f"\n  Output directory: {results_dir}/")
    print(f"    07_final_candidates.csv    -- ranked candidates + SD")
    print(f"    07_final_candidates.fa     -- FASTA sequences")
    print(f"    07_pareto_frontier.csv     -- Pareto-optimal set")
    print(f"    07_selection_summary.png   -- overview plots")
    print(f"    07_holo_vs_apo_scatter.png -- fitness-colored")
    print(f"    07_holo_vs_apo_sd.png     -- SD-colored")
    print(f"    07_sd_overview.png         -- SD distributions")
    print(f"\nNext step: python 08_export_for_synthesis.py")


if __name__ == "__main__":
    main()

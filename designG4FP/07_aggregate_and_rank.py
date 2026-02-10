#!/usr/bin/env python3
"""
Step 7: Aggregate results across all template structures and rank candidates.

Pools per-structure bound-vs-apo comparison data from 06_compare_ligand_states.py,
computes a G4FP fitness score, and selects the top 100-1000 candidates using:
  1. Hard filters (bound must fold, apo should be disordered)
  2. Pareto-optimal frontier across multiple objectives
  3. Diversity-weighted selection to avoid redundant sequences

Outputs:
  07_results/
    07_final_candidates.csv      -- ranked candidates with all metrics
    07_final_candidates.fa       -- FASTA of selected sequences
    07_pareto_frontier.csv       -- Pareto-optimal designs
    07_selection_summary.png     -- overview plots
    07_bound_vs_apo_scatter.png  -- bound pLDDT vs apo pLDDT
    07_pareto_plot.png           -- Pareto frontier visualization
    07_diversity_clusters.png    -- sequence clustering
    plot_data_csvs/              -- raw data for Prism

Usage:
    python 07_aggregate_and_rank.py
    python 07_aggregate_and_rank.py --n-select 500
    python 07_aggregate_and_rank.py --output-dir output_G4FP_des1_cro_mod0  # single structure
"""

import sys
import argparse
import warnings
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster, linkage

warnings.filterwarnings("ignore")
sns.set_style("whitegrid")
plt.rcParams["figure.figsize"] = (12, 8)
plt.rcParams["font.size"] = 10


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_comparison_results(output_dir: Path) -> Optional[pd.DataFrame]:
    """Load 06_ligand_state_comparison_results.csv from one output directory."""
    csv_path = output_dir / "06_ligand_state_comparison_results.csv"
    if not csv_path.exists():
        return None
    df = pd.read_csv(csv_path)
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
        comparison_df = load_comparison_results(odir)
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
        print(f"  Loaded {odir.name}: {len(comparison_df)} designs")

    if not frames:
        print("ERROR: No comparison results found in any output directory.")
        print("       Run 06_compare_ligand_states.py first.")
        sys.exit(1)

    df = pd.concat(frames, ignore_index=True)
    # Create globally unique identifier
    df["global_id"] = df["template"] + "_design_" + df["seq_id"].astype(str).str.zfill(4)
    print(f"\nTotal designs across all templates: {len(df)}")
    return df


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def compute_fitness_score(df: pd.DataFrame,
                          w_plddt_diff: float = 0.35,
                          w_bound_ptm: float = 0.20,
                          w_apo_disorder: float = 0.20,
                          w_chrom_plddt_diff: float = 0.15,
                          w_rmsd: float = 0.10) -> pd.DataFrame:
    """
    Compute a composite G4FP fitness score.

    Higher score = better G4-dependent folding switch.

    Components (all normalized 0-1 before weighting):
      - delta_plddt:  mean_plddt_diff (bound - apo), higher = better
      - bound_ptm:    bound state pTM, higher = better fold
      - apo_disorder: (100 - apo_mean_plddt)/100, higher = more disordered apo
      - chrom_diff:   chromophore_plddt_diff, higher = better local switch
      - rmsd:         global_rmsd between states, moderate values preferred
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

    # Bound pTM (higher = better fold when bound)
    if "bound_ptm" in df.columns:
        df["_n_bound_ptm"] = _norm(df["bound_ptm"].fillna(0))
    else:
        df["_n_bound_ptm"] = 0

    # Apo disorder (higher = more disordered when unbound)
    if "apo_mean_plddt" in df.columns:
        df["_n_apo_disorder"] = _norm(100 - df["apo_mean_plddt"].fillna(100))
    else:
        df["_n_apo_disorder"] = 0

    # Chromophore pLDDT difference
    if "chromophore_plddt_diff" in df.columns:
        df["_n_chrom_diff"] = _norm(df["chromophore_plddt_diff"].fillna(0))
    else:
        df["_n_chrom_diff"] = 0

    # RMSD between states (some conformational change is good)
    if "global_rmsd" in df.columns:
        df["_n_rmsd"] = _norm(df["global_rmsd"].fillna(0))
    else:
        df["_n_rmsd"] = 0

    df["fitness_score"] = (
        w_plddt_diff * df["_n_plddt_diff"]
        + w_bound_ptm * df["_n_bound_ptm"]
        + w_apo_disorder * df["_n_apo_disorder"]
        + w_chrom_plddt_diff * df["_n_chrom_diff"]
        + w_rmsd * df["_n_rmsd"]
    )

    # Drop normalization intermediates
    df.drop(columns=[c for c in df.columns if c.startswith("_n_")], inplace=True)

    return df


# ---------------------------------------------------------------------------
# Hard filters
# ---------------------------------------------------------------------------

def apply_hard_filters(df: pd.DataFrame,
                       min_bound_plddt: float = 70.0,
                       min_bound_ptm: float = 0.5,
                       max_apo_plddt: float = 60.0) -> pd.DataFrame:
    """
    Apply hard cutoffs. Designs that fail any filter are excluded.

    Defaults are intentionally moderate -- tighten for stricter selection.
    """
    n_before = len(df)
    mask = pd.Series(True, index=df.index)

    if "bound_mean_plddt" in df.columns:
        mask &= df["bound_mean_plddt"] >= min_bound_plddt
    if "bound_ptm" in df.columns:
        mask &= df["bound_ptm"] >= min_bound_ptm
    if "apo_mean_plddt" in df.columns:
        mask &= df["apo_mean_plddt"] <= max_apo_plddt

    df_filtered = df[mask].copy()
    n_after = len(df_filtered)
    print(f"\nHard filters: {n_before} -> {n_after} designs "
          f"(removed {n_before - n_after})")
    print(f"  bound_plddt >= {min_bound_plddt}, bound_ptm >= {min_bound_ptm}, "
          f"apo_plddt <= {max_apo_plddt}")

    if n_after == 0:
        print("\nWARNING: All designs filtered out! Relaxing filters to keep top 50% by fitness_score.")
        df_sorted = df.sort_values("fitness_score", ascending=False)
        df_filtered = df_sorted.head(len(df_sorted) // 2).copy()
        print(f"  Kept {len(df_filtered)} designs by fitness score fallback.")

    return df_filtered


# ---------------------------------------------------------------------------
# Pareto frontier
# ---------------------------------------------------------------------------

def pareto_frontier(df: pd.DataFrame,
                    objectives: List[Tuple[str, bool]]) -> pd.DataFrame:
    """
    Find the Pareto-optimal designs across multiple objectives.

    Args:
        df: DataFrame with design metrics
        objectives: list of (column_name, maximize). True=maximize, False=minimize.

    Returns:
        DataFrame subset of Pareto-optimal (non-dominated) designs.
    """
    # Build objective matrix (all converted to "higher is better")
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

    # Find non-dominated set
    is_pareto = np.ones(len(df), dtype=bool)
    for i in range(len(df)):
        if not is_pareto[i]:
            continue
        for j in range(len(df)):
            if i == j or not is_pareto[j]:
                continue
            # j dominates i if j >= i on all objectives and j > i on at least one
            if np.all(vals[j] >= vals[i]) and np.any(vals[j] > vals[i]):
                is_pareto[i] = False
                break

    pareto_df = df[is_pareto].copy()
    pareto_df["is_pareto"] = True
    print(f"\nPareto frontier: {len(pareto_df)} non-dominated designs "
          f"(from {len(df)} total)")
    obj_strs = ["{} ({})".format(c, "max" if m else "min") for c, m in available]
    print(f"  Objectives: {obj_strs}")
    return pareto_df


# ---------------------------------------------------------------------------
# Diversity-weighted selection
# ---------------------------------------------------------------------------

def sequence_hamming_distance(seqs: List[str]) -> np.ndarray:
    """Compute pairwise Hamming distance matrix for aligned sequences."""
    n = len(seqs)
    max_len = max(len(s) for s in seqs)
    # Pad shorter sequences
    padded = [s.ljust(max_len, "X") for s in seqs]
    arr = np.array([list(s) for s in padded])
    # Pairwise distance
    dist = np.zeros((n, n))
    for i in range(n):
        dist[i, :] = np.sum(arr[i] != arr, axis=1) / max_len
    return dist


def diversity_weighted_selection(df: pd.DataFrame,
                                  n_select: int,
                                  score_col: str = "fitness_score",
                                  seq_col: str = "sequence",
                                  n_clusters: int = 0) -> pd.DataFrame:
    """
    Select top designs while maximizing sequence diversity.

    Strategy:
      1. Cluster sequences by Hamming distance
      2. Within each cluster, pick the top-scoring design
      3. Fill remaining slots by iterating clusters round-robin, picking next-best

    If n_clusters=0, auto-select as min(n_select, n_unique_sequences // 2).
    """
    if seq_col not in df.columns or df[seq_col].isna().all():
        print("  WARNING: no sequence data for diversity selection, using score only")
        return df.nlargest(n_select, score_col).copy()

    df = df.copy()
    seqs = df[seq_col].tolist()
    n = len(seqs)

    if n <= n_select:
        print(f"  Only {n} designs available, returning all")
        df["cluster"] = 0
        return df

    # Compute distance matrix
    print(f"  Computing pairwise distances for {n} sequences...")
    dist_matrix = sequence_hamming_distance(seqs)

    # Hierarchical clustering
    condensed = squareform(dist_matrix)
    Z = linkage(condensed, method="average")

    if n_clusters <= 0:
        n_clusters = min(n_select, max(10, n // 5))
    n_clusters = min(n_clusters, n)

    labels = fcluster(Z, t=n_clusters, criterion="maxclust")
    df["cluster"] = labels

    print(f"  Clustered into {len(set(labels))} groups")

    # Select: top from each cluster, then round-robin
    selected_indices = []
    remaining_by_cluster = {}

    for cluster_id in sorted(set(labels)):
        cluster_df = df[df["cluster"] == cluster_id].sort_values(score_col, ascending=False)
        indices = cluster_df.index.tolist()
        if indices:
            selected_indices.append(indices[0])
            remaining_by_cluster[cluster_id] = indices[1:]

    # If we already have enough from one-per-cluster
    if len(selected_indices) >= n_select:
        selected_indices = selected_indices[:n_select]
    else:
        # Round-robin fill
        while len(selected_indices) < n_select:
            added_this_round = False
            for cluster_id in sorted(remaining_by_cluster.keys()):
                if len(selected_indices) >= n_select:
                    break
                if remaining_by_cluster[cluster_id]:
                    selected_indices.append(remaining_by_cluster[cluster_id].pop(0))
                    added_this_round = True
            if not added_this_round:
                break

    result = df.loc[selected_indices].copy()
    result["diversity_rank"] = range(1, len(result) + 1)
    print(f"  Selected {len(result)} diverse designs from {len(set(result['cluster']))} clusters")
    return result


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def create_visualizations(df_all: pd.DataFrame, df_selected: pd.DataFrame,
                          df_pareto: pd.DataFrame, output_dir: Path):
    """Generate summary plots."""
    output_dir.mkdir(exist_ok=True, parents=True)
    csv_dir = output_dir / "plot_data_csvs"
    csv_dir.mkdir(exist_ok=True)

    # 1. Selection summary: score distributions
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Fitness score distribution
    axes[0, 0].hist(df_all["fitness_score"].dropna(), bins=40, alpha=0.5, label="All designs", color="gray")
    axes[0, 0].hist(df_selected["fitness_score"].dropna(), bins=40, alpha=0.7, label="Selected", color="steelblue")
    axes[0, 0].set_xlabel("Fitness Score")
    axes[0, 0].set_ylabel("Count")
    axes[0, 0].set_title("Fitness Score Distribution")
    axes[0, 0].legend()

    # Designs per template
    template_counts = df_selected["template"].value_counts()
    axes[0, 1].barh(range(len(template_counts)), template_counts.values, color="steelblue")
    axes[0, 1].set_yticks(range(len(template_counts)))
    labels = [t.replace("output_G4FP_", "") for t in template_counts.index]
    axes[0, 1].set_yticklabels(labels, fontsize=8)
    axes[0, 1].set_xlabel("Selected Designs")
    axes[0, 1].set_title("Selected Designs per Template")

    # Delta pLDDT distribution
    if "mean_plddt_diff" in df_selected.columns:
        axes[1, 0].hist(df_selected["mean_plddt_diff"].dropna(), bins=30, alpha=0.7, color="green")
        axes[1, 0].set_xlabel("pLDDT Difference (Bound - Apo)")
        axes[1, 0].set_ylabel("Count")
        axes[1, 0].set_title("Selected: pLDDT Switch Magnitude")
        axes[1, 0].axvline(0, color="red", linestyle="--", alpha=0.5)

    # Cluster sizes
    if "cluster" in df_selected.columns:
        cluster_sizes = df_selected["cluster"].value_counts().sort_index()
        axes[1, 1].bar(range(len(cluster_sizes)), cluster_sizes.values, alpha=0.7, color="purple")
        axes[1, 1].set_xlabel("Cluster ID")
        axes[1, 1].set_ylabel("Designs Selected")
        axes[1, 1].set_title("Diversity: Designs per Sequence Cluster")

    plt.tight_layout()
    plt.savefig(output_dir / "07_selection_summary.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  07_selection_summary.png")

    # 2. Bound vs Apo scatter
    if "bound_mean_plddt" in df_all.columns and "apo_mean_plddt" in df_all.columns:
        fig, ax = plt.subplots(figsize=(10, 8))

        # All designs (gray)
        ax.scatter(df_all["bound_mean_plddt"], df_all["apo_mean_plddt"],
                   alpha=0.15, s=15, color="gray", label="All designs")

        # Selected (colored by fitness)
        sc = ax.scatter(df_selected["bound_mean_plddt"], df_selected["apo_mean_plddt"],
                        c=df_selected["fitness_score"], cmap="viridis",
                        alpha=0.7, s=30, edgecolors="black", linewidths=0.3,
                        label="Selected")
        plt.colorbar(sc, ax=ax, label="Fitness Score")

        # Pareto frontier
        if not df_pareto.empty and "bound_mean_plddt" in df_pareto.columns:
            pareto_sorted = df_pareto.sort_values("bound_mean_plddt")
            ax.plot(pareto_sorted["bound_mean_plddt"], pareto_sorted["apo_mean_plddt"],
                    "r-", alpha=0.5, linewidth=1.5, label="Pareto frontier")
            ax.scatter(df_pareto["bound_mean_plddt"], df_pareto["apo_mean_plddt"],
                       c="red", s=40, marker="D", alpha=0.6, zorder=5)

        # y=x line
        lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lims, lims, "k--", alpha=0.3, label="y=x")

        # Ideal region annotation
        ax.annotate("IDEAL\n(high bound,\nlow apo)",
                     xy=(85, 35), fontsize=10, color="green",
                     ha="center", style="italic",
                     bbox=dict(boxstyle="round,pad=0.3", fc="lightgreen", alpha=0.3))

        ax.set_xlabel("Bound Mean pLDDT")
        ax.set_ylabel("Apo Mean pLDDT")
        ax.set_title("G4FP Design Space: Bound vs Apo Structural Confidence")
        ax.legend(loc="upper left")
        plt.tight_layout()
        plt.savefig(output_dir / "07_bound_vs_apo_scatter.png", dpi=300, bbox_inches="tight")
        plt.close()
        print("  07_bound_vs_apo_scatter.png")

    # 3. Pareto frontier plot (two main objectives)
    if not df_pareto.empty:
        fig, ax = plt.subplots(figsize=(10, 8))
        obj1 = "mean_plddt_diff"
        obj2 = "bound_ptm"

        if obj1 in df_all.columns and obj2 in df_all.columns:
            ax.scatter(df_all[obj1], df_all[obj2],
                       alpha=0.15, s=15, color="gray", label="All")
            ax.scatter(df_pareto[obj1], df_pareto[obj2],
                       c="red", s=50, marker="D", alpha=0.7, label="Pareto optimal")
            ax.scatter(df_selected[obj1], df_selected[obj2],
                       alpha=0.4, s=20, color="steelblue", label="Selected")
            ax.set_xlabel("pLDDT Difference (Bound - Apo)")
            ax.set_ylabel("Bound pTM")
            ax.set_title("Pareto Frontier: pLDDT Switch vs Bound Fold Quality")
            ax.legend()
            plt.tight_layout()
            plt.savefig(output_dir / "07_pareto_plot.png", dpi=300, bbox_inches="tight")
            plt.close()
            print("  07_pareto_plot.png")

    # Export CSVs
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
        description="Aggregate G4FP designs across templates and rank by bound-vs-apo switch",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Selection strategy:
  1. Load 06_ligand_state_comparison_results.csv from all output directories
  2. Compute fitness score (weighted: delta_pLDDT, bound_pTM, apo_disorder, chromophore, RMSD)
  3. Apply hard filters (bound must fold, apo should be disordered)
  4. Find Pareto-optimal designs across multiple objectives
  5. Diversity-weighted selection to avoid redundant sequences
  6. Export top N candidates for experimental testing

Pareto filtering selects designs that are non-dominated across objectives
(no other design is better on ALL metrics simultaneously).

Diversity-weighted selection clusters sequences by similarity and picks
the best from each cluster, ensuring broad coverage of sequence space.
        """
    )
    parser.add_argument("--output-dir", type=str,
                        help="Analyze a single output directory instead of all")
    parser.add_argument("--n-select", type=int, default=500,
                        help="Number of final candidates to select (default: 500)")
    parser.add_argument("--min-bound-plddt", type=float, default=70.0,
                        help="Minimum bound state pLDDT (default: 70)")
    parser.add_argument("--min-bound-ptm", type=float, default=0.5,
                        help="Minimum bound state pTM (default: 0.5)")
    parser.add_argument("--max-apo-plddt", type=float, default=60.0,
                        help="Maximum apo state pLDDT (default: 60)")
    parser.add_argument("--no-diversity", action="store_true",
                        help="Skip diversity-weighted selection, rank by score only")
    parser.add_argument("--pareto-only", action="store_true",
                        help="Only output Pareto-optimal designs (ignore n-select)")

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

    # 2. Compute fitness score
    print("\nComputing fitness scores...")
    df = compute_fitness_score(df)
    df = df.sort_values("fitness_score", ascending=False)
    print(f"  Score range: {df['fitness_score'].min():.4f} - {df['fitness_score'].max():.4f}")
    print(f"  Median: {df['fitness_score'].median():.4f}")

    # 3. Pareto frontier
    print("\nFinding Pareto-optimal designs...")
    pareto_objectives = [
        ("mean_plddt_diff", True),        # maximize: bigger switch
        ("bound_ptm", True),               # maximize: good fold when bound
        ("apo_mean_plddt", False),         # minimize: disordered when unbound
        ("chromophore_plddt_diff", True),  # maximize: chromophore switches
        ("global_rmsd", True),             # maximize: structural change
    ]
    df_pareto = pareto_frontier(df, pareto_objectives)

    if args.pareto_only:
        df_selected = df_pareto.sort_values("fitness_score", ascending=False)
    else:
        # 4. Hard filters
        df_filtered = apply_hard_filters(
            df,
            min_bound_plddt=args.min_bound_plddt,
            min_bound_ptm=args.min_bound_ptm,
            max_apo_plddt=args.max_apo_plddt,
        )

        # 5. Diversity-weighted selection
        if args.no_diversity or "sequence" not in df_filtered.columns:
            df_selected = df_filtered.nlargest(args.n_select, "fitness_score").copy()
            print(f"\nSelected top {len(df_selected)} by fitness score (no diversity weighting)")
        else:
            print(f"\nApplying diversity-weighted selection (target: {args.n_select})...")
            df_selected = diversity_weighted_selection(
                df_filtered, n_select=args.n_select
            )

    # 6. Final ranking
    df_selected = df_selected.sort_values("fitness_score", ascending=False).reset_index(drop=True)
    df_selected["final_rank"] = range(1, len(df_selected) + 1)

    # 7. Save outputs
    results_dir.mkdir(exist_ok=True, parents=True)

    # CSV with all metrics
    export_cols = [c for c in df_selected.columns if "array" not in c]
    df_selected[export_cols].to_csv(results_dir / "07_final_candidates.csv", index=False)

    # FASTA
    if "sequence" in df_selected.columns:
        with open(results_dir / "07_final_candidates.fa", "w") as f:
            for _, row in df_selected.iterrows():
                header = (f">{row['global_id']} "
                          f"rank={row['final_rank']} "
                          f"fitness={row['fitness_score']:.4f} "
                          f"template={row['template']}")
                if "mean_plddt_diff" in row:
                    header += f" delta_plddt={row['mean_plddt_diff']:.1f}"
                f.write(f"{header}\n{row['sequence']}\n")

    # Pareto frontier
    if not df_pareto.empty:
        pareto_export = [c for c in df_pareto.columns if "array" not in c]
        df_pareto[pareto_export].to_csv(results_dir / "07_pareto_frontier.csv", index=False)

    # Visualizations
    print("\nGenerating plots...")
    create_visualizations(df, df_selected, df_pareto, results_dir)

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

    print(f"\n  Output directory: {results_dir}/")
    print(f"    07_final_candidates.csv   -- ranked candidates")
    print(f"    07_final_candidates.fa    -- FASTA sequences")
    print(f"    07_pareto_frontier.csv    -- Pareto-optimal set")
    print(f"    07_selection_summary.png  -- overview plots")
    print(f"    07_bound_vs_apo_scatter.png")
    print(f"\nNext step: python 08_export_for_synthesis.py")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Generate a publication-quality multi-panel overview figure for the G4FP
design-analysis pipeline.

Panels:
  A  (top, full width)   Pipeline flow chart
  B  (bottom left)       Structure comparison schematic
  C  (bottom right)      Multi-objective selection strategy

Requires only matplotlib (no external images).

Usage:
    python3 make_model_figure.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
from pathlib import Path


# ---------------------------------------------------------------------------
# Colour palette
# ---------------------------------------------------------------------------
TAB10 = plt.cm.tab10.colors  # 10 colours for templates
STAGE_COLORS = {
    "input":    "#4e79a7",  # steel blue
    "design":   "#59a14f",  # green
    "filter":   "#edc948",  # gold
    "predict":  "#e15759",  # coral red
    "analyse":  "#76b7b2",  # teal
    "rank":     "#b07aa1",  # purple
}
HOLO_CLR = "#2a7bbe"
APO_CLR  = "#d95f02"
PARETO_CLR = "#e41a1c"
BG_CLR   = "white"


# ---------------------------------------------------------------------------
# Helper: draw a rounded box with centred text
# ---------------------------------------------------------------------------
def _box(ax, x, y, w, h, text, color, fontsize=9, text_color="white",
         alpha=1.0, style="round,pad=0.08", lw=1.2, zorder=3):
    """Draw a FancyBboxPatch and centred text, return patch centre."""
    box = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle=style,
        facecolor=color, edgecolor="0.25", linewidth=lw, alpha=alpha,
        zorder=zorder,
    )
    ax.add_patch(box)
    ax.text(x, y, text, ha="center", va="center",
            fontsize=fontsize, color=text_color, fontweight="bold",
            zorder=zorder + 1)
    return (x, y)


def _arrow(ax, xy_from, xy_to, color="0.35", lw=1.8, style="->",
           connectionstyle="arc3,rad=0", shrinkA=8, shrinkB=8):
    """Draw a FancyArrowPatch between two (x,y) points."""
    arrow = FancyArrowPatch(
        xy_from, xy_to,
        arrowstyle=style, color=color, linewidth=lw,
        connectionstyle=connectionstyle,
        shrinkA=shrinkA, shrinkB=shrinkB,
        zorder=2,
    )
    ax.add_patch(arrow)
    return arrow


def _multiline(ax, x, y, lines, fontsize=8, va="top", ha="left",
               color="0.15", linespacing=1.35, **kw):
    ax.text(x, y, "\n".join(lines), fontsize=fontsize, va=va, ha=ha,
            color=color, linespacing=linespacing, **kw)


# ===================================================================
# Panel A -- Pipeline flow chart
# ===================================================================
def draw_panel_a(ax):
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 3.6)
    ax.axis("off")

    # -- label --
    ax.text(0.15, 3.4, "A", fontsize=18, fontweight="bold", va="top")
    ax.text(1.0, 3.4, "G4FP Design Pipeline Overview", fontsize=13,
            va="top", fontstyle="italic", color="0.25")

    bw, bh = 2.0, 1.05  # box width / height
    y_mid = 1.55

    # Stage positions (evenly spaced)
    xs = np.linspace(1.3, 14.7, 6)

    stages = [
        ("Template PDBs\n(10 structures)", STAGE_COLORS["input"]),
        ("LigandMPNN\n(10K designs\nper template)",  STAGE_COLORS["design"]),
        ("Filter top 1%\n(100 per\ntemplate)",        STAGE_COLORS["filter"]),
        ("AlphaFold 3\n(bound + apo\n5 seeds each)", STAGE_COLORS["predict"]),
        ("Structure\nAnalysis\n(steps 05-06)",        STAGE_COLORS["analyse"]),
        ("Ranking &\nSelection\n(step 07)",           STAGE_COLORS["rank"]),
    ]

    centres = []
    for i, (label, clr) in enumerate(stages):
        c = _box(ax, xs[i], y_mid, bw, bh, label, clr, fontsize=8.5)
        centres.append(c)

    # Arrows between consecutive boxes
    for i in range(len(centres) - 1):
        _arrow(ax, (centres[i][0] + bw / 2, y_mid),
               (centres[i + 1][0] - bw / 2, y_mid),
               shrinkA=2, shrinkB=2, lw=2.2, style="-|>")

    # Numeric counts below arrow midpoints
    counts = ["100K total", "1,000 total", "2,000 jobs", "metrics + SD", "top 100-150"]
    for i, txt in enumerate(counts):
        mx = (xs[i] + xs[i + 1]) / 2
        ax.text(mx, y_mid - 0.72, txt, ha="center", va="top",
                fontsize=7.5, color="0.4", fontstyle="italic")


# ===================================================================
# Panel B -- Structure comparison schematic
# ===================================================================
def draw_panel_b(ax):
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 7)
    ax.axis("off")

    ax.text(0.15, 6.85, "B", fontsize=18, fontweight="bold", va="top")
    ax.text(0.85, 6.85, "Structural Metrics", fontsize=13, va="top",
            fontstyle="italic", color="0.25")

    # -- Holo / Apo protein schematics --
    # Holo (left)
    holo_cx, holo_cy = 2.2, 4.7
    _box(ax, holo_cx, holo_cy, 3.0, 1.4, "", HOLO_CLR, alpha=0.15,
         lw=1.5, style="round,pad=0.12")
    ax.text(holo_cx, holo_cy + 0.42, "Holo structure", ha="center",
            fontsize=9.5, fontweight="bold", color=HOLO_CLR)
    ax.text(holo_cx, holo_cy - 0.05, "(protein + G4-DNA)", ha="center",
            fontsize=8, color="0.35")

    # Apo (right)
    apo_cx, apo_cy = 6.0, 4.7
    _box(ax, apo_cx, apo_cy, 3.0, 1.4, "", APO_CLR, alpha=0.15,
         lw=1.5, style="round,pad=0.12")
    ax.text(apo_cx, apo_cy + 0.42, "Apo structure", ha="center",
            fontsize=9.5, fontweight="bold", color=APO_CLR)
    ax.text(apo_cx, apo_cy - 0.05, "(protein only)", ha="center",
            fontsize=8, color="0.35")

    # Double-headed arrow between holo and apo
    _arrow(ax, (holo_cx + 1.5, holo_cy), (apo_cx - 1.5, holo_cy),
           color="0.3", lw=1.8, style="<|-|>", shrinkA=4, shrinkB=4)
    ax.text((holo_cx + apo_cx) / 2, holo_cy + 0.22, "compare",
            ha="center", va="bottom", fontsize=7.5, fontstyle="italic",
            color="0.4")

    # -- Metrics table --
    metrics_y = 3.2
    ax.text(0.4, metrics_y, "Metrics computed per design (5 AF3 seeds):",
            fontsize=9, fontweight="bold", color="0.2")

    metric_rows = [
        ("Global RMSD",           "holo vs apo backbone C\u03b1 overlay"),
        ("Chromophore RMSD",      "TYG region (auto-detected position)"),
        ("Mean pLDDT",            "per-residue confidence (holo & apo)"),
        ("Chromophore pLDDT",     "TYG pLDDT in holo vs apo"),
        ("pTM",                   "predicted TM-score (fold quality)"),
        ("iPTM",                  "interface TM-score (protein\u2013DNA)"),
        ("Fraction disordered",   "from AF3 summary_confidences.json"),
        ("Standard deviation",    "across 5 seeds (reproducibility)"),
    ]

    for j, (name, desc) in enumerate(metric_rows):
        yy = metrics_y - 0.34 * (j + 1)
        ax.plot(0.55, yy, marker="s", markersize=4.5,
                color=STAGE_COLORS["analyse"], clip_on=False)
        ax.text(0.80, yy, name, fontsize=8, fontweight="bold", va="center",
                color="0.15")
        ax.text(3.4, yy, desc, fontsize=7.5, va="center", color="0.4")

    # -- Chromophore callout box --
    chrom_y = 0.35
    chrom_box = FancyBboxPatch(
        (0.3, chrom_y - 0.25), 7.3, 0.55,
        boxstyle="round,pad=0.1",
        facecolor="#fef3c7", edgecolor="#d97706", linewidth=1.2,
        zorder=3,
    )
    ax.add_patch(chrom_box)
    ax.text(4.0, chrom_y, "Chromophore (TYG) position auto-detected per template "
            "from LigandMPNN reference FASTA  (residues ~160\u2013172)",
            ha="center", va="center", fontsize=7.5, color="#92400e", zorder=4)


# ===================================================================
# Panel C -- Selection strategy
# ===================================================================
def draw_panel_c(ax):
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 7)
    ax.axis("off")

    ax.text(0.15, 6.85, "C", fontsize=18, fontweight="bold", va="top")
    ax.text(0.85, 6.85, "Selection Strategy", fontsize=13, va="top",
            fontstyle="italic", color="0.25")

    # -- Three selection paths --
    path_y = 5.65
    path_h = 0.80

    # Path 1: Top 10 per template
    _box(ax, 2.0, path_y, 3.3, path_h,
         "Top 10 per template\n(100 designs)", STAGE_COLORS["input"],
         fontsize=8.5)

    # Path 2: Top 50 global
    _box(ax, 6.0, path_y, 3.3, path_h,
         "Top 50 global\n(by fitness score)", STAGE_COLORS["design"],
         fontsize=8.5)

    # Union arrow
    _arrow(ax, (2.0, path_y - path_h / 2), (4.0, path_y - 1.3),
           color="0.4", lw=1.5, style="-|>",
           connectionstyle="arc3,rad=-0.15", shrinkA=6, shrinkB=6)
    _arrow(ax, (6.0, path_y - path_h / 2), (4.0, path_y - 1.3),
           color="0.4", lw=1.5, style="-|>",
           connectionstyle="arc3,rad=0.15", shrinkA=6, shrinkB=6)

    # Union box
    _box(ax, 4.0, path_y - 1.65, 3.0, 0.7,
         "Union (~140 designs)", STAGE_COLORS["rank"],
         fontsize=9)

    # -- Pareto frontier --
    pareto_y = 2.95
    pareto_box = FancyBboxPatch(
        (0.3, pareto_y - 0.55), 7.3, 1.1,
        boxstyle="round,pad=0.12",
        facecolor="#fee2e2", edgecolor=PARETO_CLR, linewidth=1.5,
        zorder=2,
    )
    ax.add_patch(pareto_box)
    ax.text(4.0, pareto_y + 0.35, "Pareto Frontier (3 objectives)",
            ha="center", fontsize=10, fontweight="bold", color=PARETO_CLR,
            zorder=3)
    pareto_objs = [
        ("iPTM",              "protein\u2013DNA interface quality"),
        ("Chromophore pLDDT",  "TYG confidence in holo state"),
        ("Global RMSD",        "holo vs apo structural change"),
    ]
    for j, (name, desc) in enumerate(pareto_objs):
        yy = pareto_y - 0.04 - 0.27 * j
        ax.text(1.0, yy, f"{j+1}. {name}", fontsize=8.5, fontweight="bold",
                va="center", color="0.15", zorder=3)
        ax.text(4.2, yy, desc, fontsize=7.5, va="center",
                color="0.45", zorder=3)

    # -- Fitness function --
    fit_y = 1.55
    ax.text(0.4, fit_y, "Fitness score (weighted sum of normalised components):",
            fontsize=9, fontweight="bold", color="0.2")

    weights = [
        ("Chrom. holo pLDDT", 0.25, STAGE_COLORS["predict"]),
        ("iPTM",              0.25, STAGE_COLORS["predict"]),
        ("Global RMSD",       0.15, STAGE_COLORS["analyse"]),
        ("Holo pTM",          0.15, STAGE_COLORS["analyse"]),
        ("Rel. apo pLDDT",    0.10, APO_CLR),
        ("Reproducibility",   0.05, "0.55"),
        ("Chrom. pLDDT diff", 0.05, "0.55"),
    ]

    bar_x0, bar_w_total = 0.5, 6.8
    bar_y, bar_h = 1.15, 0.30
    cum_x = bar_x0
    for name, w, clr in weights:
        seg_w = bar_w_total * w
        rect = FancyBboxPatch(
            (cum_x, bar_y), seg_w, bar_h,
            boxstyle="square,pad=0",
            facecolor=clr, edgecolor="white", linewidth=0.8,
            alpha=0.85, zorder=3,
        )
        ax.add_patch(rect)
        # Label inside or below
        if w >= 0.12:
            ax.text(cum_x + seg_w / 2, bar_y + bar_h / 2,
                    f"{name}\n{int(w*100)}%",
                    ha="center", va="center", fontsize=6.5,
                    color="white", fontweight="bold", zorder=4)
        else:
            ax.text(cum_x + seg_w / 2, bar_y - 0.08,
                    f"{name} {int(w*100)}%",
                    ha="center", va="top", fontsize=5.8,
                    color="0.3", zorder=4, rotation=0)
        cum_x += seg_w

    # Hard-filter footnote
    ax.text(0.4, 0.30,
            "Hard filters:  holo pLDDT \u2265 70  |  holo pTM \u2265 0.5  "
            "|  holo iPTM \u2265 0.3",
            fontsize=7.5, color="0.45", fontstyle="italic")


# ===================================================================
# Main: compose panels
# ===================================================================
def main():
    fig = plt.figure(figsize=(16, 10), facecolor=BG_CLR)

    # Panel A: top row, full width
    ax_a = fig.add_axes([0.02, 0.68, 0.96, 0.30])   # [left, bottom, w, h]

    # Panel B: bottom-left
    ax_b = fig.add_axes([0.02, 0.02, 0.48, 0.64])

    # Panel C: bottom-right
    ax_c = fig.add_axes([0.52, 0.02, 0.46, 0.64])

    draw_panel_a(ax_a)
    draw_panel_b(ax_b)
    draw_panel_c(ax_c)

    # Light horizontal separator between top and bottom rows
    fig.patches.append(mpatches.FancyBboxPatch(
        (0.02, 0.665), 0.96, 0.003,
        boxstyle="square,pad=0", facecolor="0.80",
        transform=fig.transFigure, zorder=0,
    ))

    out_path = Path(__file__).resolve().parent / "07_results" / "pipeline_overview.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight",
                facecolor=BG_CLR, edgecolor="none")
    plt.close(fig)
    print(f"Saved  {out_path}  ({out_path.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()

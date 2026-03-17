#!/usr/bin/env python3
"""
Generate a multi-panel cartoon schematic showing what the G4FP analysis
measures and how top designs are selected.

Panels:
  A  (top-left)     Holo structure cartoon with annotated metrics
  B  (top-right)    Apo structure cartoon + holo-vs-apo comparison
  C  (bottom-left)  Pareto frontier cartoon (3 objectives)
  D  (bottom-right) Selection funnel: how top designs are picked

All drawings use matplotlib primitives (no external images).

Usage:
    python3 make_model_figure.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Arc
from matplotlib.collections import PatchCollection
import numpy as np
from pathlib import Path

# ── Colours ──
HOLO_CLR = "#2a7bbe"
APO_CLR = "#d95f02"
DNA_CLR = "#7b9e3a"
CHROM_CLR = "#fbbf24"  # yellow
ION_CLR = "#9333ea"    # purple
PARETO_CLR = "#e41a1c"
SELECTED_CLR = "#4e79a7"
BG_CLR = "white"
METRIC_CLR = "#1e3a5f"
LIGHT_BLUE = "#dbeafe"
LIGHT_ORANGE = "#ffedd5"
LIGHT_GREEN = "#dcfce7"
LIGHT_RED = "#fee2e2"


def _box(ax, x, y, w, h, text, color, fontsize=9, text_color="white",
         alpha=1.0, style="round,pad=0.08", lw=1.2, zorder=3):
    box = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h, boxstyle=style,
        facecolor=color, edgecolor="0.25", linewidth=lw, alpha=alpha,
        zorder=zorder)
    ax.add_patch(box)
    ax.text(x, y, text, ha="center", va="center",
            fontsize=fontsize, color=text_color, fontweight="bold",
            zorder=zorder + 1)
    return (x, y)


def _arrow(ax, xy_from, xy_to, color="0.35", lw=1.8, style="-|>",
           connectionstyle="arc3,rad=0", shrinkA=8, shrinkB=8, zorder=2):
    arrow = FancyArrowPatch(
        xy_from, xy_to, arrowstyle=style, color=color, linewidth=lw,
        connectionstyle=connectionstyle, shrinkA=shrinkA, shrinkB=shrinkB,
        zorder=zorder)
    ax.add_patch(arrow)


def _callout(ax, text, xy_point, xy_text, fontsize=7.5, color=METRIC_CLR,
             arrowcolor="0.5"):
    """Draw an annotation with an arrow from xy_text to xy_point."""
    ax.annotate(text, xy=xy_point, xytext=xy_text,
                fontsize=fontsize, color=color, fontweight="bold",
                arrowprops=dict(arrowstyle="-|>", color=arrowcolor,
                                lw=1.0, shrinkB=3),
                bbox=dict(boxstyle="round,pad=0.15", fc="white",
                          ec=arrowcolor, alpha=0.9, lw=0.8),
                zorder=10)


def _draw_protein_blob(ax, cx, cy, rx, ry, color, alpha=0.3, lw=2):
    """Draw an irregular protein-shaped blob using a deformed ellipse."""
    theta = np.linspace(0, 2 * np.pi, 200)
    # Add bumps to make it look organic
    r_mod = (1.0
             + 0.08 * np.sin(5 * theta)
             + 0.05 * np.cos(7 * theta)
             + 0.04 * np.sin(11 * theta))
    x = cx + rx * r_mod * np.cos(theta)
    y = cy + ry * r_mod * np.sin(theta)
    ax.fill(x, y, color=color, alpha=alpha, zorder=2)
    ax.plot(x, y, color=color, lw=lw, alpha=0.7, zorder=3)


def _draw_helix(ax, x0, y0, length, angle=0, color="#4488cc", lw=2):
    """Draw a cartoon alpha helix as a wavy ribbon."""
    t = np.linspace(0, length, 100)
    amp = 0.12
    freq = 8
    # Ribbon edges
    x1 = x0 + t * np.cos(angle) - amp * np.sin(freq * t) * np.sin(angle)
    y1 = y0 + t * np.sin(angle) + amp * np.sin(freq * t) * np.cos(angle)
    x2 = x0 + t * np.cos(angle) + amp * np.sin(freq * t) * np.sin(angle)
    y2 = y0 + t * np.sin(angle) - amp * np.sin(freq * t) * np.cos(angle)
    ax.fill_between(x1, y1, y2, alpha=0.5, color=color, zorder=4)
    ax.plot(x1, y1, color=color, lw=lw, alpha=0.8, zorder=5)
    ax.plot(x2, y2, color=color, lw=lw * 0.6, alpha=0.5, zorder=5)


def _draw_beta_strand(ax, x0, y0, length, angle=0, color="#4488cc", lw=1.5):
    """Draw a cartoon beta strand as a flat arrow."""
    dx = length * np.cos(angle)
    dy = length * np.sin(angle)
    w = 0.08
    # Arrow body
    nx, ny = -np.sin(angle) * w, np.cos(angle) * w
    body = plt.Polygon([
        (x0 + nx, y0 + ny),
        (x0 + dx * 0.7 + nx, y0 + dy * 0.7 + ny),
        (x0 + dx * 0.7 + nx * 2, y0 + dy * 0.7 + ny * 2),
        (x0 + dx, y0 + dy),
        (x0 + dx * 0.7 - nx * 2, y0 + dy * 0.7 - ny * 2),
        (x0 + dx * 0.7 - nx, y0 + dy * 0.7 - ny),
        (x0 - nx, y0 - ny),
    ], closed=True, fc=color, ec=color, alpha=0.6, lw=lw, zorder=4)
    ax.add_patch(body)


def _draw_dna_quadruplex(ax, cx, cy, size=0.5, color=DNA_CLR):
    """Draw a G-quadruplex as four stacked squares with connecting lines."""
    offsets = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
    s = size * 0.3
    layers = 3
    for layer in range(layers):
        ly = cy + (layer - 1) * size * 0.35
        for ox, oy in offsets:
            x = cx + ox * s
            y = ly + oy * s * 0.3
            rect = plt.Rectangle((x - s * 0.3, y - s * 0.15), s * 0.6, s * 0.3,
                                  fc=color, ec="0.3", lw=0.8, alpha=0.6, zorder=5)
            ax.add_patch(rect)
        # Connect corners
        for i in range(4):
            j = (i + 1) % 4
            ax.plot([cx + offsets[i][0] * s, cx + offsets[j][0] * s],
                    [ly + offsets[i][1] * s * 0.3, ly + offsets[j][1] * s * 0.3],
                    color=color, lw=0.8, alpha=0.5, zorder=4)


def _draw_ions(ax, positions, color=ION_CLR, size=0.06):
    """Draw K+ ions as small spheres."""
    for x, y in positions:
        circle = Circle((x, y), size, fc=color, ec="0.3", lw=0.6,
                        alpha=0.8, zorder=6)
        ax.add_patch(circle)
        ax.text(x, y, "K\u207a", ha="center", va="center", fontsize=4,
                color="white", fontweight="bold", zorder=7)


# ═══════════════════════════════════════════════════════════════════
# Panel A: Holo structure with annotated metrics
# ═══════════════════════════════════════════════════════════════════
def draw_panel_a(ax):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis("off")

    ax.text(0.2, 7.7, "A", fontsize=20, fontweight="bold", va="top")
    ax.text(0.9, 7.7, "Holo Structure (Bound State)", fontsize=14,
            va="top", fontstyle="italic", color="0.25")

    # Draw protein blob
    pcx, pcy = 4.0, 4.2
    _draw_protein_blob(ax, pcx, pcy, 2.2, 2.0, HOLO_CLR, alpha=0.15)

    # Draw secondary structure elements inside
    _draw_helix(ax, 2.5, 5.0, 1.5, angle=0.3, color=HOLO_CLR)
    _draw_helix(ax, 3.0, 3.2, 1.2, angle=-0.2, color=HOLO_CLR)
    _draw_beta_strand(ax, 4.5, 5.2, 1.0, angle=-0.5, color=HOLO_CLR)
    _draw_beta_strand(ax, 4.8, 3.8, 0.9, angle=0.4, color=HOLO_CLR)
    _draw_beta_strand(ax, 3.8, 4.5, 0.8, angle=-0.1, color=HOLO_CLR)

    # Chromophore (TYG) highlighted region
    chrom_cx, chrom_cy = 4.2, 4.3
    chrom = FancyBboxPatch((chrom_cx - 0.3, chrom_cy - 0.2), 0.6, 0.4,
                            boxstyle="round,pad=0.05",
                            fc=CHROM_CLR, ec="#b45309", lw=1.5,
                            alpha=0.9, zorder=6)
    ax.add_patch(chrom)
    ax.text(chrom_cx, chrom_cy, "TYG", ha="center", va="center",
            fontsize=6, fontweight="bold", color="#78350f", zorder=7)

    # G-quadruplex DNA
    dna_cx, dna_cy = 6.2, 4.5
    _draw_dna_quadruplex(ax, dna_cx, dna_cy, size=0.7, color=DNA_CLR)
    ax.text(dna_cx, dna_cy - 0.8, "G4-DNA", ha="center", fontsize=7,
            color=DNA_CLR, fontweight="bold")

    # K+ ions
    _draw_ions(ax, [(dna_cx, dna_cy - 0.15), (dna_cx, dna_cy + 0.15),
                    (dna_cx, dna_cy + 0.45)])

    # Interface region
    ax.annotate("", xy=(5.6, 4.5), xytext=(5.0, 4.5),
                arrowprops=dict(arrowstyle="<->", color="0.4", lw=1.2))

    # ── Metric callouts ──
    _callout(ax, "pLDDT\n(per-residue\nconfidence)",
             (2.8, 5.3), (0.3, 6.8), fontsize=7)
    _callout(ax, "pTM\n(fold quality\nTM-score)",
             (3.5, 2.8), (0.3, 1.5), fontsize=7)
    _callout(ax, "iPTM\n(protein\u2013DNA\ninterface)",
             (5.3, 4.5), (7.5, 6.8), fontsize=7)
    _callout(ax, "Chromophore\npLDDT (TYG\nresidues 197-199)",
             (chrom_cx + 0.3, chrom_cy + 0.2), (7.5, 5.3), fontsize=7,
             color="#b45309")
    _callout(ax, "Fraction\ndisordered",
             (2.5, 3.5), (0.3, 2.8), fontsize=7)

    # SD note
    ax.text(5.0, 0.5, "All metrics computed across\n25 AF3 seed replicates\n"
            "\u2192 mean \u00b1 SD per design",
            ha="center", va="center", fontsize=8, color="0.35",
            fontstyle="italic",
            bbox=dict(boxstyle="round,pad=0.3", fc=LIGHT_BLUE,
                      ec="0.7", alpha=0.7))


# ═══════════════════════════════════════════════════════════════════
# Panel B: Holo vs Apo comparison
# ═══════════════════════════════════════════════════════════════════
def draw_panel_b(ax):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis("off")

    ax.text(0.2, 7.7, "B", fontsize=20, fontweight="bold", va="top")
    ax.text(0.9, 7.7, "Holo vs Apo Comparison", fontsize=14,
            va="top", fontstyle="italic", color="0.25")

    # ── Holo (left) ──
    hcx, hcy = 2.5, 5.0
    _draw_protein_blob(ax, hcx, hcy, 1.3, 1.2, HOLO_CLR, alpha=0.15)
    _draw_helix(ax, 1.8, 5.3, 0.8, angle=0.2, color=HOLO_CLR)
    _draw_beta_strand(ax, 2.5, 4.6, 0.6, angle=-0.3, color=HOLO_CLR)

    # G4 DNA next to holo
    _draw_dna_quadruplex(ax, 3.7, 5.2, size=0.4, color=DNA_CLR)

    # Chromophore
    chrom_h = FancyBboxPatch((hcx - 0.15, hcy - 0.12), 0.3, 0.24,
                              boxstyle="round,pad=0.03",
                              fc=CHROM_CLR, ec="#b45309", lw=1, alpha=0.9, zorder=6)
    ax.add_patch(chrom_h)

    ax.text(hcx, 3.5, "HOLO\n(+ G4-DNA)", ha="center", fontsize=9,
            fontweight="bold", color=HOLO_CLR)

    # ── Apo (right) ──
    acx, acy = 7.5, 5.0
    _draw_protein_blob(ax, acx, acy, 1.3, 1.2, APO_CLR, alpha=0.15)
    # Apo is slightly unfolded — offset helices
    _draw_helix(ax, 6.9, 5.4, 0.7, angle=0.4, color=APO_CLR)
    _draw_beta_strand(ax, 7.6, 4.5, 0.5, angle=-0.1, color=APO_CLR)

    # Chromophore (dimmer in apo)
    chrom_a = FancyBboxPatch((acx - 0.15, acy - 0.12), 0.3, 0.24,
                              boxstyle="round,pad=0.03",
                              fc=CHROM_CLR, ec="#b45309", lw=1, alpha=0.4, zorder=6)
    ax.add_patch(chrom_a)

    ax.text(acx, 3.5, "APO\n(no DNA)", ha="center", fontsize=9,
            fontweight="bold", color=APO_CLR)

    # ── Comparison arrows ──
    mid_x = 5.0
    _arrow(ax, (3.9, 5.0), (6.1, 5.0), color="0.3", lw=2.5,
           style="<|-|>", shrinkA=15, shrinkB=15)

    # RMSD label
    ax.text(mid_x, 5.35, "Global RMSD", ha="center", fontsize=8,
            fontweight="bold", color=METRIC_CLR)
    ax.text(mid_x, 5.05, "(C\u03b1 backbone)", ha="center", fontsize=7,
            color="0.4")

    # ── What we're looking for ──
    goal_y = 2.2
    goal_box = FancyBboxPatch((0.5, goal_y - 0.8), 9.0, 1.6,
                               boxstyle="round,pad=0.15",
                               fc=LIGHT_GREEN, ec="#16a34a", lw=1.5,
                               zorder=2)
    ax.add_patch(goal_box)
    ax.text(5.0, goal_y + 0.55, "Design Goal: G4-dependent folding",
            ha="center", fontsize=10, fontweight="bold", color="#15803d", zorder=3)

    goals = [
        "\u2713  High holo pLDDT + pTM  (folds well with G4-DNA)",
        "\u2713  Low apo pLDDT  (disordered without G4-DNA)",
        "\u2713  Large RMSD between holo and apo  (structural switch)",
        "\u2713  High chromophore pLDDT in holo  (fluorescent when bound)",
    ]
    for j, txt in enumerate(goals):
        ax.text(1.0, goal_y + 0.15 - 0.32 * j, txt, fontsize=7.5,
                color="#166534", va="center", zorder=3)


# ═══════════════════════════════════════════════════════════════════
# Panel C: Pareto frontier cartoon
# ═══════════════════════════════════════════════════════════════════
def draw_panel_c(ax):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis("off")

    ax.text(0.2, 7.7, "C", fontsize=20, fontweight="bold", va="top")
    ax.text(0.9, 7.7, "Pareto Frontier (3 Objectives)", fontsize=14,
            va="top", fontstyle="italic", color="0.25")

    # Draw a fake scatter plot with Pareto frontier
    plot_ax = ax.inset_axes([0.1, 0.35, 0.55, 0.55])
    np.random.seed(42)

    # Random designs
    n = 80
    x = np.random.beta(2, 5, n) * 0.8 + 0.1
    y = np.random.beta(2, 5, n) * 0.8 + 0.1

    # Pareto front (upper-right corner)
    px = np.array([0.95, 0.85, 0.72, 0.55, 0.40, 0.30])
    py = np.array([0.30, 0.55, 0.70, 0.82, 0.88, 0.95])

    plot_ax.scatter(x, y, s=15, c="0.65", alpha=0.5, zorder=2,
                    label="All designs")
    plot_ax.scatter(px, py, s=60, c=PARETO_CLR, marker="*", edgecolors="0.2",
                    linewidths=0.5, zorder=5, label="Pareto-optimal")
    plot_ax.plot(px, py, color=PARETO_CLR, lw=1.5, ls="--", alpha=0.7, zorder=4)

    # Fill the dominated region
    fill_x = np.concatenate([[0], px, [0]])
    fill_y = np.concatenate([[py[0]], py, [py[-1]]])
    plot_ax.fill(fill_x, fill_y, color=PARETO_CLR, alpha=0.06)

    plot_ax.set_xlabel("iPTM  \u2192", fontsize=8, color="0.3")
    plot_ax.set_ylabel("Chromophore pLDDT  \u2192", fontsize=8, color="0.3")
    plot_ax.set_xlim(0, 1.05)
    plot_ax.set_ylim(0, 1.05)
    plot_ax.set_xticks([])
    plot_ax.set_yticks([])
    plot_ax.spines['top'].set_visible(False)
    plot_ax.spines['right'].set_visible(False)
    plot_ax.legend(fontsize=6, loc="lower right")

    ax.text(3.2, 7.0, "(cartoon — axes show\n2 of 3 objectives)",
            fontsize=6, color="0.5", fontstyle="italic")

    # Objective labels on the right
    obj_x = 7.0
    ax.text(obj_x, 6.5, "3 Objectives:", fontsize=10, fontweight="bold",
            color=METRIC_CLR)

    objectives = [
        ("1. iPTM", "protein\u2013DNA\ninterface quality", HOLO_CLR),
        ("2. Chromophore pLDDT", "TYG confidence\nin holo state", "#b45309"),
        ("3. Global RMSD", "structural change\nholo \u2194 apo", "0.4"),
    ]
    for j, (name, desc, clr) in enumerate(objectives):
        yy = 5.8 - j * 1.3
        ax.plot(obj_x - 0.3, yy, marker="o", markersize=10, color=clr,
                zorder=5)
        ax.text(obj_x - 0.3, yy, str(j + 1), ha="center", va="center",
                fontsize=7, color="white", fontweight="bold", zorder=6)
        ax.text(obj_x + 0.1, yy + 0.15, name, fontsize=9, fontweight="bold",
                color=clr, va="center")
        ax.text(obj_x + 0.1, yy - 0.2, desc, fontsize=7, color="0.45",
                va="top")

    # Non-dominated explanation
    ax.text(obj_x, 2.2,
            "Pareto-optimal = not\ndominated on ALL 3\nobjectives simultaneously",
            fontsize=7.5, color=PARETO_CLR, fontstyle="italic",
            bbox=dict(boxstyle="round,pad=0.3", fc=LIGHT_RED,
                      ec=PARETO_CLR, alpha=0.7, lw=1))


# ═══════════════════════════════════════════════════════════════════
# Panel D: Selection funnel
# ═══════════════════════════════════════════════════════════════════
def draw_panel_d(ax):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis("off")

    ax.text(0.2, 7.7, "D", fontsize=20, fontweight="bold", va="top")
    ax.text(0.9, 7.7, "Design Selection & Ranking", fontsize=14,
            va="top", fontstyle="italic", color="0.25")

    # ── Funnel stages ──
    stages = [
        ("100K sequences\n(LigandMPNN)", 8.0, "0.75", "white"),
        ("1,000 filtered\n(top 1% confidence)", 6.5, "#4e79a7", "white"),
        ("All 1,000 scored\n(fitness function)", 5.0, "#59a14f", "white"),
    ]

    for text, width, color, tc in stages:
        idx = stages.index((text, width, color, tc))
        y = 6.8 - idx * 1.3
        _box(ax, 5.0, y, width, 0.85, text, color, fontsize=8.5,
             text_color=tc, lw=1.5)
        if idx < len(stages) - 1:
            _arrow(ax, (5.0, y - 0.43), (5.0, y - 0.87),
                   color="0.4", lw=2, style="-|>", shrinkA=2, shrinkB=2)

    # ── Fitness function breakdown ──
    fit_y = 3.0
    ax.text(1.0, fit_y + 0.4, "Fitness Score (weighted sum):", fontsize=9,
            fontweight="bold", color="0.2")

    weights = [
        ("Chrom. holo pLDDT", 0.25, "#e15759"),
        ("iPTM",              0.25, "#e15759"),
        ("Global RMSD",       0.15, "#76b7b2"),
        ("Holo pTM",          0.15, "#76b7b2"),
        ("Rel. apo pLDDT",    0.10, APO_CLR),
        ("Reproducibility",   0.05, "0.55"),
        ("Chrom. diff",       0.05, "0.55"),
    ]

    bar_x0, bar_w = 1.0, 8.0
    bar_y, bar_h = fit_y - 0.3, 0.35
    cum_x = bar_x0
    for name, w, clr in weights:
        seg_w = bar_w * w
        rect = FancyBboxPatch(
            (cum_x, bar_y), seg_w, bar_h,
            boxstyle="square,pad=0",
            facecolor=clr, edgecolor="white", linewidth=0.8,
            alpha=0.85, zorder=3)
        ax.add_patch(rect)
        if w >= 0.12:
            ax.text(cum_x + seg_w / 2, bar_y + bar_h / 2,
                    f"{name}\n{int(w * 100)}%",
                    ha="center", va="center", fontsize=6,
                    color="white", fontweight="bold", zorder=4)
        else:
            ax.text(cum_x + seg_w / 2, bar_y - 0.08,
                    f"{name} {int(w * 100)}%",
                    ha="center", va="top", fontsize=5.5, color="0.3", zorder=4)
        cum_x += seg_w

    # ── Selection outputs ──
    out_y = 1.5
    _arrow(ax, (5.0, fit_y - 0.7), (5.0, out_y + 0.55),
           color="0.4", lw=2, style="-|>", shrinkA=2, shrinkB=2)

    # Three output boxes
    boxes = [
        (2.0, out_y, "Top 10 per\ntemplate\n(100 total)", "#4e79a7"),
        (5.0, out_y, "Top 50\nglobal\n(by fitness)", "#59a14f"),
        (8.0, out_y, "Pareto\nfrontier\n(~30-40)", PARETO_CLR),
    ]
    for bx, by, txt, clr in boxes:
        _box(ax, bx, by, 2.5, 0.9, txt, clr, fontsize=7.5, lw=1.5)

    # Hard filter note
    ax.text(5.0, 0.3,
            "Hard filters:  holo pLDDT \u2265 70  |  holo pTM \u2265 0.5  "
            "|  holo iPTM \u2265 0.3",
            ha="center", fontsize=7, color="0.45", fontstyle="italic")


# ═══════════════════════════════════════════════════════════════════
# Main: compose 2x2 panels
# ═══════════════════════════════════════════════════════════════════
def main():
    fig = plt.figure(figsize=(18, 14), facecolor=BG_CLR)

    # 2x2 grid
    ax_a = fig.add_axes([0.01, 0.52, 0.48, 0.46])
    ax_b = fig.add_axes([0.51, 0.52, 0.48, 0.46])
    ax_c = fig.add_axes([0.01, 0.02, 0.48, 0.48])
    ax_d = fig.add_axes([0.51, 0.02, 0.48, 0.48])

    draw_panel_a(ax_a)
    draw_panel_b(ax_b)
    draw_panel_c(ax_c)
    draw_panel_d(ax_d)

    # Separators
    for y in [0.51]:
        fig.patches.append(mpatches.FancyBboxPatch(
            (0.02, y), 0.96, 0.002,
            boxstyle="square,pad=0", facecolor="0.82",
            transform=fig.transFigure, zorder=0))
    for x in [0.50]:
        fig.patches.append(mpatches.FancyBboxPatch(
            (x, 0.03), 0.002, 0.94,
            boxstyle="square,pad=0", facecolor="0.82",
            transform=fig.transFigure, zorder=0))

    out_dir = Path(__file__).resolve().parent / "07_results"
    out_dir.mkdir(parents=True, exist_ok=True)

    png_path = out_dir / "pipeline_overview.png"
    svg_path = out_dir / "pipeline_overview.svg"
    fig.savefig(png_path, dpi=300, bbox_inches="tight",
                facecolor=BG_CLR, edgecolor="none")
    fig.savefig(svg_path, bbox_inches="tight",
                facecolor=BG_CLR, edgecolor="none")
    plt.close(fig)
    print(f"Saved {png_path} ({png_path.stat().st_size / 1024:.0f} KB)")
    print(f"Saved {svg_path}")


if __name__ == "__main__":
    main()

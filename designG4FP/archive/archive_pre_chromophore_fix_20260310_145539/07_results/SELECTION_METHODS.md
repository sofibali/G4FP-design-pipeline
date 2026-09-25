# G4FP Design Selection: Methods and Interpretation

## Overview

Script `07_aggregate_and_rank.py` selects the best G4FP candidates from 1,000 designs (100 per template x 10 templates) using a multi-stage pipeline:

1. **Fitness scoring** — weighted combination of normalized metrics
2. **Hard filters** — remove designs that clearly won't work
3. **Pareto frontier** — find designs that are optimal across multiple objectives simultaneously
4. **Two selection pools** — top 10 per template + top 50 global by fitness

### Output Files

| File | Contents |
|------|----------|
| `07_top10_per_template.csv` + `.fa` | 100 designs: best 10 from each of the 10 templates |
| `07_top50_global.csv` + `.fa` | 50 designs: highest fitness scores regardless of template |
| `07_final_candidates.csv` + `.fa` | 140 designs: union of both pools (with 10 overlap) |
| `07_pareto_frontier.csv` | 37 Pareto-optimal designs |

---

## 1. Fitness Score

Each design gets a single fitness score (0–1) computed as a weighted sum of normalized components. Every raw metric is min-max normalized to [0, 1] across all 1,000 designs before weighting.

### Components and Weights

| Weight | Component | Raw Column | Direction | Why It Matters |
|--------|-----------|-----------|-----------|----------------|
| **25%** | Chromophore holo pLDDT | `holo_chromophore_plddt` | Higher = better | Chromophore (residues 197-199) must be well-positioned in the holo state for fluorescence |
| **25%** | iPTM | `holo_iptm` | Higher = better | Measures predicted quality of the protein-DNA interface. Higher iPTM = AF3 is more confident the protein binds the G4 DNA correctly |
| **15%** | Global RMSD | `global_rmsd` | Higher = better | RMSD between holo and apo structures. Larger structural change means the protein genuinely restructures upon G4 binding |
| **15%** | Holo pTM | `holo_ptm` | Higher = better | Overall fold quality of the holo state. The protein must fold well when bound to G4 |
| **10%** | Relative apo pLDDT | `apo_plddt_vs_template` | Lower = better | Apo pLDDT *relative to the parent template*. Negative means the design's apo state is less stable than the template's — good for a switch. We use relative (not absolute) because all designs have similar absolute apo pLDDT (~82-85) |
| **5%** | Confidence | `mean_plddt_diff_sd` | Lower SD = better | Reproducibility across the 25 AF3 seed-sample predictions. Low SD means AF3 consistently predicts the same structure |
| **5%** | Chromophore pLDDT diff | `chromophore_plddt_diff` | Higher = better | Difference in chromophore pLDDT between holo and apo. Larger diff = chromophore is specifically stabilized by G4 binding |

### Normalization

For each component:
```
normalized = (value - min) / (max - min)
```
where min/max are computed across all 1,000 designs. A score of 1.0 is the best design for that metric, 0.0 is the worst. For metrics where lower is better (relative apo pLDDT, SD), the normalization is inverted: `1 - normalized`.

### Score Component Columns

The final CSV includes `score_component_*` columns showing each design's normalized (0-1) value for each component. These are what get multiplied by the weights and summed to produce `fitness_score`.

---

## 2. Hard Filters

Before selection, designs that fail basic quality thresholds are removed:

| Filter | Default | Rationale |
|--------|---------|-----------|
| `holo_mean_plddt >= 70` | The holo state must fold | A protein that doesn't fold when bound to G4 is useless |
| `holo_ptm >= 0.5` | Minimum fold quality | AF3 pTM below 0.5 indicates poor predicted structure |
| `holo_iptm >= 0.3` | Minimum interface quality | Must have some predicted interaction with G4 DNA |

**No absolute apo filter** — we don't filter by absolute apo pLDDT because all designs have similar values (82-85). Instead, we rank by *relative* apo pLDDT in the fitness score.

If all designs get filtered out (shouldn't happen with current thresholds), the script falls back to the top 50% by fitness score.

**Current result:** 999 of 1,000 designs pass (1 removed for low iPTM).

---

## 3. Pareto Frontier

### What is Pareto Optimality?

A design is **Pareto-optimal** (non-dominated) if no other design is better in *all* objectives simultaneously. To improve any one metric of a Pareto-optimal design, you'd have to sacrifice at least one other metric.

### Example

Consider two objectives: iPTM and chromophore pLDDT.

- Design A: iPTM=0.8, chromophore=92
- Design B: iPTM=0.6, chromophore=98
- Design C: iPTM=0.5, chromophore=90

A and B are both Pareto-optimal: A has better iPTM, B has better chromophore, neither dominates the other. C is **dominated** by A (A is better in both metrics), so C is not Pareto-optimal.

### Our Pareto Objectives

We use 3 objectives (more objectives = nearly everything becomes non-dominated):

| Objective | Direction | Meaning |
|-----------|-----------|---------|
| `holo_iptm` | Maximize | Best protein-G4 DNA interface |
| `holo_chromophore_plddt` | Maximize | Best chromophore positioning |
| `global_rmsd` | Maximize | Largest structural change upon binding |

### How to Read the Pareto Plots

**`07_pareto_plot.png`** shows two 2D projections of the 3D Pareto frontier:

- **Left panel: iPTM vs Chromophore pLDDT** — The red step-line traces the Pareto frontier (the upper-right envelope). Red diamonds are the individual Pareto-optimal designs. Anything below/left of the red line is dominated. Blue dots are the selected candidates.

- **Right panel: iPTM vs Global RMSD** — Same idea. The red step-line traces the frontier. Top-right = high interface quality AND large structural change.

The background scatter (colored by fitness score) shows all 1,000 designs. The Pareto frontier defines the absolute best trade-off boundary — you can't improve one axis without worsening the other.

**Current result:** 37 Pareto-optimal designs out of 1,000.

### Pareto vs Fitness Score

These are complementary:
- **Fitness score** combines all metrics into one number using predefined weights. Good for ranking.
- **Pareto frontier** makes no assumptions about relative importance. It finds designs that represent the *best possible trade-offs* regardless of weighting.

A design can be Pareto-optimal but have a mediocre fitness score (e.g., extreme in one metric but poor in others). And a high-fitness design might not be Pareto-optimal if another design slightly beats it in every dimension.

---

## 4. Selection Pools

After hard filtering, two independent selection pools are built:

### Pool A: Top 10 Per Template (100 designs)

The 10 highest-fitness designs from each of the 10 templates. This guarantees representation from every template, since different templates may explore different regions of sequence/structure space.

**Why this matters:** If selection were purely by fitness, a single dominant template could take all slots (as happens with Pool B — all 50 come from `r20_l1_cro_mod0`). Pool A ensures we don't miss good designs from less dominant templates that might still perform well experimentally.

### Pool B: Top 50 Global (50 designs)

The 50 highest-fitness designs across all templates, regardless of origin. This is the pure "best predicted" set with no template balancing.

**Current result:** All 50 come from template `r20_l1_cro_mod0`, which dominates the fitness ranking.

### Combined Output

`07_final_candidates.csv` contains the **union** of Pool A and Pool B (140 designs, since 10 overlap). This gives you both the globally best candidates and diverse template coverage for experimental testing.

---

## 5. Plot Guide

### `07_selection_summary.png` (6 panels)

| Panel | What It Shows |
|-------|---------------|
| (a) Fitness Score Distribution | Gray = all 1,000 designs, blue = 140 selected, red = 37 Pareto. Selected should be shifted right |
| (b) Fitness by Template | Gray violins = all designs per template, blue boxes = selected. Shows which templates produce the best candidates |
| (c) Score Components Radar | Gray = median of all designs, blue = median of selected. Blue should be larger (better) on all axes |
| (d) iPTM vs Chromophore by Template | Each color = one template. Large dots = selected. Shows whether certain templates dominate the ideal (top-right) corner |
| (e) Top 50 by Rank | Horizontal bars colored by template. Shows rank-ordered fitness scores and which templates contribute |
| (f) Key Metric Distributions | Box plots of the 5 most important metrics for selected designs, normalized to [0,1] with raw ranges below |

### `07_top10_per_template_heatmap.png`

100 rows (10 per template), grouped by template with black separator lines. Columns sorted left-to-right by fitness weight importance:

1. Fitness Score
2. Chromophore pLDDT (25%)
3. iPTM (25%)
4. RMSD (15%)
5. pTM (15%)
6. Apo pLDDT vs Template (10%)
7. pLDDT Diff SD (5%)
8. Chromophore Diff (5%)
9. Holo Fraction Disordered
10. Apo Fraction Disordered

Color: **green = better, red = worse**. Normalized across all 1,000 designs (not just the displayed rows) so colors are comparable across templates. Raw values printed in each cell.

- **`P`** = Pareto-optimal design
- **`*`** = in the selected pool (union of A + B)

**How to use it:** Compare template sections. A template whose top 10 are mostly green is producing better candidates overall. Look for rows green in the first few columns (highest-weight components).

### `07_top50_global_heatmap.png`

50 rows ranked by fitness score. Same columns and color scheme as above. Shows the absolute best candidates across all templates. `P` marks Pareto-optimal.

**How to use it:** This is your shortlist. Look for consistently green rows. Designs near the top with some red cells may have a specific weakness worth investigating.

### `07_pareto_plot.png`

Two panels showing 2D projections of the 3D Pareto frontier. Red step-line = frontier boundary. Red diamonds = Pareto-optimal designs. Blue dots = selected. Background colored by fitness score.

### `07_fitness_components.png`

Histograms of each raw metric feeding into the fitness score, plus fraction_disordered. Gray = all, colored = selected, red = Pareto medians. Shows how selection enriches for favorable values.

### `07_holo_vs_apo_scatter.png`

X = holo pLDDT, Y = apo pLDDT. Colored by fitness. Stars = template references. Ideal G4FP is bottom-right (high holo, low apo).

### `07_holo_vs_apo_sd.png`

Same axes, colored by prediction uncertainty (SD). Darker = more reproducible.

### `07_sd_overview.png`

Histograms of all SD columns, showing prediction reproducibility distributions.

---

## 6. Key Metrics Glossary

| Metric | Source | Range | Meaning |
|--------|--------|-------|---------|
| `holo_mean_plddt` | AF3 confidences.json | 0-100 | Per-residue confidence of the holo (G4-bound) structure |
| `apo_mean_plddt` | AF3 confidences.json | 0-100 | Per-residue confidence of the apo (no G4) structure |
| `holo_chromophore_plddt` | AF3 confidences.json | 0-100 | pLDDT specifically at the CRO chromophore (residues 197-199) |
| `holo_ptm` | AF3 summary_confidences.json | 0-1 | Predicted TM-score of the holo state (fold quality) |
| `holo_iptm` | AF3 summary_confidences.json | 0-1 | Predicted interface TM-score (protein-DNA interaction quality) |
| `global_rmsd` | Superimposition of holo vs apo CA atoms | 0-20+ A | How much the backbone moves between bound and unbound states |
| `chromophore_rmsd` | Superimposition of chromophore CA atoms | 0-5+ A | Chromophore-specific structural change |
| `fraction_disordered` | AF3 summary_confidences.json | 0-1 | Fraction of residues predicted as disordered |
| `apo_plddt_vs_template` | Computed | -10 to +10 | Design's apo pLDDT minus its parent template's apo pLDDT |
| `*_sd` columns | Std dev across 25 seeds | varies | Prediction reproducibility (lower = more confident) |

---

## 7. Tuning the Selection

### Adjusting Weights

Edit the weight arguments in `compute_fitness_score()`:

```python
df = compute_fitness_score(df,
    w_chrom_holo=0.25,   # chromophore stability
    w_iptm=0.25,         # protein-DNA interface
    w_rmsd=0.15,         # structural change
    w_holo_ptm=0.15,     # fold quality
    w_rel_apo=0.10,      # relative apo destabilization
    w_confidence=0.05,   # prediction reproducibility
    w_chrom_diff=0.05,   # chromophore-specific switch
)
```

### Adjusting Hard Filters

```bash
python 07_aggregate_and_rank.py \
    --min-holo-plddt 75 \
    --min-holo-ptm 0.6 \
    --min-holo-iptm 0.4 \
    --max-apo-plddt 80      # optional absolute apo cutoff
```

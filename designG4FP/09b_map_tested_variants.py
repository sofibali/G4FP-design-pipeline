#!/usr/bin/env python3
"""Map the 21 experimentally tested round-1 designs onto literature-supported metrics.

Reads 09_literature_metrics.csv (all 1000 designs) and round2/tested_variants.csv
(the 21 ordered designs and their outcomes), and asks a single question:

    on metrics the literature actually validates, where do the tested designs
    land -- and is there any separation between the 2 that fluoresced and the
    19 that did not?

Statistics are exact permutation tests over all C(21,2) = 210 ways of choosing
which 2 of the 21 were positive, so the smallest attainable p-value is
1/210 = 0.0048. With 2 positives nothing here can be significant after
multiplicity correction; the p-values are reported to bound the claim, not to
support one.
"""

from __future__ import annotations

import argparse
import itertools
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# CVD-validated palette (OKLab dE x100: normal 31.2, deut 28.0, prot 22.3)
C_DARK = '#0072B2'
C_FLUO = '#D55E00'
C_BG = '#9AA0A6'
INK, INK2, GRID = '#202124', '#5F6368', '#DADCE0'

# metric -> (label, reference line, reference label, better direction)
PANELS = [
    ('holo_pae_interaction_ordered', 'Interface PAE, ordered residues (Å)', 10.0,
     'Bennett 2023 threshold (<10)', 'lower'),
    ('apo_mean_plddt', 'Apo mean pLDDT', 60.0,
     'documented filter, never implemented (<60)', 'lower'),
    ('holo_chromophore_cyclization_distance', 'Chromophore cyclization N→C (Å)', None,
     None, 'lower'),
    ('holo_min_contact_distance', 'Min protein–DNA heavy-atom distance (Å)', 2.2,
     'steric floor', 'higher'),
    ('holo_n_core_contacts', 'Ordered-residue contacts, 3.0–4.5 Å', None, None, 'higher'),
    ('holo_chromophore_plddt', 'Chromophore pLDDT (corrected)', None, None, 'higher'),
]


def exact_p(values_pos, values_all):
    """Two-sided exact permutation p over all C(n,2) choices of 2 positives."""
    n = len(values_all)
    obs = np.mean(values_pos)
    stats = [np.mean([values_all[i], values_all[j]])
             for i, j in itertools.combinations(range(n), 2)]
    stats = np.array(stats)
    centre = np.mean(stats)
    return float(np.mean(np.abs(stats - centre) >= abs(obs - centre)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--metrics', default='09_literature_metrics.csv')
    ap.add_argument('--tested', default='round2/tested_variants.csv')
    ap.add_argument('--out-prefix', default='09_results/09_tested_map')
    args = ap.parse_args()

    df = pd.read_csv(args.metrics)
    tv = pd.read_csv(args.tested)
    tv['outcome'] = tv['outcome'].str.strip()
    fluo = set(tv.loc[tv.outcome.str.startswith('fluorescent'), 'global_id'])
    dark = set(tv.loc[tv.outcome == 'dark', 'global_id'])
    tested = fluo | dark

    df['group'] = np.where(df.global_id.isin(fluo), 'fluorescent',
                  np.where(df.global_id.isin(dark), 'dark', 'untested'))
    t = df[df.global_id.isin(tested)].copy()
    print(f"matched {len(t)} of {len(tested)} tested designs "
          f"({(t.group=='fluorescent').sum()} fluorescent, {(t.group=='dark').sum()} dark)")
    missing = tested - set(df.global_id)
    if missing:
        print("  NOT FOUND in metrics table:", sorted(missing))

    Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)

    # ---------------- statistics table ----------------
    stats_rows = []
    for col, label, *_ in PANELS:
        if col not in df.columns:
            continue
        sub = t.dropna(subset=[col])
        if len(sub) < 10 or (sub.group == 'fluorescent').sum() < 2:
            continue
        vals = sub[col].to_numpy()
        pos = sub.loc[sub.group == 'fluorescent', col].to_numpy()
        neg = sub.loc[sub.group == 'dark', col].to_numpy()
        ranks = pd.Series(vals).rank().to_numpy()
        pos_ranks = sorted(ranks[(sub.group == 'fluorescent').to_numpy()])
        stats_rows.append({
            'metric': col,
            'label': label,
            'fluorescent_mean': pos.mean(),
            'dark_mean': neg.mean(),
            'difference': pos.mean() - neg.mean(),
            'positive_ranks_of_21': '/'.join(str(int(r)) for r in pos_ranks),
            'exact_p': exact_p(pos, vals),
            'all1000_median': df[col].median(),
            'all1000_p5': df[col].quantile(0.05),
            'all1000_p95': df[col].quantile(0.95),
        })
    st = pd.DataFrame(stats_rows).sort_values('exact_p')
    st.to_csv(f'{args.out_prefix}_stats.csv', index=False)
    print(f"\nwrote {args.out_prefix}_stats.csv")
    print(st[['metric', 'fluorescent_mean', 'dark_mean',
              'positive_ranks_of_21', 'exact_p']].to_string(index=False))

    # ---------------- figure ----------------
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.2))
    fig.patch.set_facecolor('white')

    for ax, (col, label, ref, reflabel, better) in zip(axes.ravel(), PANELS):
        if col not in df.columns:
            ax.axis('off')
            continue
        allv = df[col].dropna()
        ax.hist(allv, bins=45, color=C_BG, alpha=0.45,
                label=f'all designs (n={len(allv)})', zorder=1)
        ymax = ax.get_ylim()[1]

        d = t.loc[t.group == 'dark', col].dropna()
        f = t.loc[t.group == 'fluorescent', col].dropna()
        ax.scatter(d, np.full(len(d), ymax * 0.52), marker='o', s=34,
                   facecolor=C_DARK, edgecolor='white', linewidth=0.8,
                   label=f'tested, dark (n={len(d)})', zorder=3)
        ax.scatter(f, np.full(len(f), ymax * 0.66), marker='*', s=260,
                   facecolor=C_FLUO, edgecolor='white', linewidth=0.9,
                   label=f'tested, fluorescent (n={len(f)})', zorder=4)

        if ref is not None:
            ax.axvline(ref, color=INK2, linestyle='--', linewidth=1.3, zorder=2)
            ax.text(ref, ymax * 0.96, f'  {reflabel}', fontsize=7.5,
                    color=INK2, va='top', ha='left', rotation=0)

        ax.set_xlabel(label, fontsize=9.5, color=INK)
        ax.set_ylabel('designs', fontsize=9, color=INK2)
        ax.tick_params(labelsize=8, colors=INK2)
        ax.grid(axis='y', color=GRID, linewidth=0.6, alpha=0.8)
        ax.set_axisbelow(True)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        for s in ('left', 'bottom'):
            ax.spines[s].set_color(GRID)

    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=3, frameon=False,
               fontsize=9.5, bbox_to_anchor=(0.5, 1.005))
    fig.suptitle('Round-1 G4FP designs on literature-supported metrics',
                 fontsize=13, color=INK, y=0.945, x=0.085, ha='left')
    fig.text(0.085, 0.905,
             'Stars are the two designs that fluoresced; neither responded to G4. '
             'Exact permutation p-values over C(21,2)=210 are in the stats CSV.',
             fontsize=8.8, color=INK2, ha='left')
    fig.tight_layout(rect=[0, 0, 1, 0.885])
    for ext in ('png', 'svg'):
        fig.savefig(f'{args.out_prefix}.{ext}', dpi=200,
                    bbox_inches='tight', facecolor='white')
    print(f"wrote {args.out_prefix}.png / .svg")


if __name__ == '__main__':
    main()

# What the corrected re-run changed, and what still stands

Companion to `20_RERUN_DELTA.md`, which holds the generated tables. This file is
written by hand; the numbers in it are quoted from that one.

## The situation

The committed result tables were last written on **2026-09-24** (`79aa6b3ed`,
"updated changes final results"). The chromophore-localisation and per-residue
pLDDT defect was fixed on **2026-11-15** in `91cb028b9` and `022d820e4`, and the
`token_chains` regression in `6bc4b5a18`. The tables were never regenerated, so
for seven weeks the repo carried **corrected code beside uncorrected numbers** —
the worst of the two states, because the code review passes and the tables still
lie.

This commit replaces the tables with the re-run, and publishes the delta so the
size of the correction is visible rather than silently overwritten.

## 1. The corrected columns were not off by a constant — they were a different quantity

| column | mean before | mean after | shift | r(old, new) |
|---|---:|---:|---:|---:|
| `apo_chromophore_plddt` | 91.25 | 70.48 | **−20.78** | 0.43 |
| `holo_chromophore_plddt` | 93.97 | 78.30 | **−15.67** | 0.37 |
| `chromophore_plddt_diff` | 2.72 | 7.83 | +5.11 | 0.17 |

The correlation is the important number, not the shift. If the old column had
been the right quantity with a constant offset, r would be ~1.0. At r = 0.37–0.43
the old column was **measuring something else** — it was reporting pLDDT over a
far wider set of atoms than the chromophore triad, which is why every design
looked like it had a 91–94 pLDDT chromophore. Nothing in a designed barrel is
that well predicted; the implausibility was visible in the old numbers and was
not caught.

`chromophore_plddt_diff` at r = 0.17 is effectively uncorrelated with its own
previous version.

## 2. The ranking it fed is about half new

| | value |
|---|---:|
| Spearman rho, old vs new ordering | **0.779** |
| top-10 overlap | **3 / 10** |
| top-50 overlap | **22 / 50** |
| top-100 overlap | 66 / 100 |
| `in_pareto` | 36 → **27** (20 kept, 16 dropped, 7 added) |

A rho of 0.78 across 1,000 designs sounds reassuring and is not: at the only end
of the distribution anyone uses, **seven of the old top ten are gone**. Any
shortlist drawn off the pre-fix tables is more than half different now.

## 3. Where the two fluorescent designs went — the only external check

| design | rank before | rank after | score before | score after |
|---|---:|---:|---:|---:|
| `G4FP_r28_l2_cro_mod2_design_0016` | 139 | **179** | 0.5684 | 0.4540 |
| `G4FP_r28_l2_cro_mod2_design_0072` | 169 | **325** | 0.5524 | 0.3882 |

Out of 1,000. **The correction moved both measured positives down.** Fixing a
real bug made the score agree with experiment *less*, which is the opposite of
what a repaired metric is supposed to do, and it is consistent with the n=21
result in `19_ROUND1_RESCORE_RESULTS.md` (AUC 0.158, exact permutation p = 0.076).

Two points of honesty about that. The fix is still correct — `chromophore_plddt`
now measures the chromophore, and that is true independently of which direction
it moves the ranking. And two positives cannot establish a direction; this is an
observation that the score has no demonstrated signal, not a finding that it is
reliably anti-correlated.

## 4. The RMSD correction was conditional on the template

`chromophore_rmsd` changed on **exactly six of ten templates, 100/100 designs in
each, and on none of the other four**:

| changed | unchanged |
|---|---|
| `r20_l1_cro_mod0`, `r20_l2_cro_mod2`, `r20_l2_cro_mod3`, `r28_l1_cro_mod2`, `r28_l2_cro_mod2`, `r28_l2_cro_mod3` | `des1_cro_mod0`, `des1_cro_mod4`, `r28_l3_cro_mod0`, `r28_l3_cro_mod2` |

All-or-nothing per template means the defect was conditional on something
template-specific — residue numbering or chain order — not a numerical drift.
The four untouched templates' RMSD numbers in the old tables still stand.

Note that `des1_cro_mod4`, which supplies 47 of the 71 survivors, is one of the
untouched four. Its dominance of the survivor set is not an artifact of this bug.

## 5. What this does not change

- **The gate architecture.** The cascade still reconciles at 999 → 298 → 197 →
  101 → 71 (`19_ROUND1_RESCORE_RESULTS.md`).
- **The diversity ceiling.** 6 of 10 templates yield a survivor; 47 of 71 from one
  template. That is a property of the library, not of the scoring.
- **The round-2 framing.** Nothing here depends on the fitness score being right —
  the argument in `17_ROUND2_DESIGN.md` is that a new campaign is forecastable and
  an iteration is not.

## 6. One thing to fix separately

`07_results/07_final_candidates.csv` contains **all 1,000 designs**, not a
candidate shortlist — the shortlist is the `in_pareto` / `in_top50_global` flag
columns inside it. The name has misled at least one reader of this repo. Renaming
it is a separate change; it is recorded here so it is not rediscovered.

## 7. Reproducing the delta

```
python3 -I round2/20_delta_rerun.py \
    --repo <repo root> \
    --tested round2/tested_variants.csv \
    --out round2/20_RERUN_DELTA.md
```

The script reads the pre-fix tables from git and the corrected ones from the
working tree, so it runs only while the re-run is uncommitted. After this commit,
compare against the parent instead: `git show <commit>^:<path>`.

No AlphaFold3, LigandMPNN or wet-lab step was run to produce this document or the
re-run it describes; the re-run is CPU analysis over structures predicted in
March 2026.

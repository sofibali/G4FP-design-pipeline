# G4FP Round 2

Everything for the second design round lives in this directory. Nothing outside
`designG4FP/round2/` is modified by this work.

Prepared November 2026. Division of labour unchanged: LigandMPNN is run on
Kielen's side and committed here; AlphaFold3 runs on Sofia's lab server.

---

## 1. Round-1 defects found during the review

Three problems in the round-1 scoring were identified by auditing the committed
code against the stored results. All three are reproduced and quantified below;
`09_rescore_and_select.py` corrects them.

### 1.1 Chromophore pLDDT was measured on the wrong residues

`compute_per_residue_plddt()` in `../05_analyze_af3_structures.py` converts
atom-level pLDDT to per-residue values by assuming every residue has the same
number of atoms:

```python
atoms_per_residue = len(chain_plddts) / n_residues
start = int(round(i * atoms_per_residue))
end   = int(round((i + 1) * atoms_per_residue))
```

Real residues have 4-14 atoms, so the mapping drifts cumulatively along the
chain. For a 305-residue G4FP construct (2403 chain-A atoms, 7.88 atoms per
residue assumed), the slice `plddt[196:199]` spans atoms 1544-1568, which belong
to **residues 201-204 (Cys-Phe-Ser-Arg)** — four residues past the chromophore,
in a well-ordered stretch of barrel.

Verification: re-implementing the original arithmetic reproduces the stored
`holo_chromophore_plddt` to within 0.2 for all 14 tested-variant structures.
Against correctly mapped values it overestimates by **12.8 pLDDT points on
average (range 6.6-19.8)**.

`holo_chromophore_plddt` carries 25% of the round-1 fitness weight and
`chromophore_plddt_diff` another 5%, so roughly **30% of the composite was
computed on the wrong residues**.

### 1.2 The chromophore is not at residues 197-199

`chromophore_start: 197` / `chromophore_end: 199` is hardcoded in
`pipeline_config.yaml`. Across the 14 tested-variant structures the Thr-Tyr-Gly
triad actually begins at residue **187, 195, 196 or 197** depending on the
construct architecture. `09_rescore_and_select.py` locates it by the invariant
cpEGFP sequence motif `VTTL(TYG)VQ`, which is construct-length independent and
works on AF3 models that contain no `CRO` heteroresidue.

This supersedes the approach on the unmerged branch
`copilot/update-chromophore-residue-selection`, which auto-detects the `CRO`
HETATM. Round-2 inputs are AF3 models with no CRO at all, so CRO detection
would find nothing in them.

### 1.3 The apo term did no work

`apo_plddt_vs_template` is NaN for all 1000 designs in the committed
`07_all_designs_scored.csv`, so `score_component_rel_apo` is exactly 0.5 for
every row — a constant. The 10% weight meant to penalise a self-stabilising apo
state contributed nothing to the ranking.

The Drive copy of that table (1020 rows) does have the metric populated, and the
repo copy is stale relative to it. **Reconcile before scoring anything off the
committed copy.**

---

## 2. What the corrected metrics show

Running `09_rescore_and_select.py` over the 36 designs that have stored
structures, validated against the 14 experimentally tested variants among them
(2 fluorescent, 12 dark):

| Score component | mean for the 2 hits | mean for the 12 dark |
|---|---|---|
| holo iPTM | 0.289 | 0.406 |
| chromophore pLDDT (corrected) | 0.345 | 0.570 |
| holo pTM | 0.045 | 0.459 |
| apo pLDDT vs template | **0.653** | **0.345** |
| global RMSD | 0.571 | 0.490 |
| seed reproducibility | 0.678 | 0.644 |

**Every structural quality metric runs backwards.** The fluorescent variants
have worse pTM, worse iPTM and lower chromophore pLDDT than the dark ones.
`apo_plddt_vs_template` is the only component that favours them, and it has to
carry about 70% of the total weight before both hits reach the top two ranks:

| weight on apo term | rank of the 2 hits (of 14) |
|---|---|
| 0.00 | 11, 14 |
| 0.15 (default) | 10, 12 |
| 0.30 | 6, 8 |
| 0.50 | 4, 6 |
| 0.70 | **1, 2** |

Read this as a warning, not a recipe. Fitting a six-term score to two positives
is overfitting, and the apo metric is partly a property of the parent template
rather than of the design. The defensible conclusion is that **the round-1 AF3
structure metrics have no demonstrated predictive value for fluorescence**, and
round-2 selection should not lean heavily on them.

### 2.1 Every AF3 holo model has protein-DNA steric clashes

All 36 stored design models — and the 10 template predictions as well — contain
protein-DNA heavy-atom pairs below 2.5 A: **10 to 37 pairs per model, minimum
separations 0.7-1.8 A**. A real hydrogen bond is 2.8-3.0 A. Example from
`r104_r28_l2_cro_mod2_d0072_s0.cif`: Lys27 CE to DT12 C7 at 1.05 A.

This is not a parsing artifact — backbone CA-CA spacing in the same files is a
correct 3.83 A. AF3 reports `has_clash = 0` for these models regardless.

Consequence: the clash gate in `09_rescore_and_select.py` is **off by default**,
because `--max-clash-pairs 0` removes the entire dataset. The count is reported
per design as a diagnostic. More importantly, it means the round-1 holo
predictions are unreliable at precisely the protein-DNA interface the whole
design depends on, which is consistent with the iPTM values clustering near 0.5.

---

## 3. Round-2 design parameters

Agreed with Kielen, November 2026.

| Parameter | Value |
|---|---|
| Design inputs | the experimentally tested round-1 variants, not the original templates |
| Designs per input | 10,000 (so 20,000 across the `r28_l2_cro_mod2` pair) |
| Redesign window | first 40 and last 40 residues of each input, symmetric by construction |
| Chromophore | AF3 models used as-is; no CRO grafting |
| Model | `ligandmpnn_v_32_010_25.pt`, ligand context on, 8.0 A score cutoff |

### Why 40 + 40

Round 1 applied one hardcoded window, `A15-60 A240-285`, to all ten templates
even though they differ in length by up to 21 residues, so the window landed
differently in each. Measured against residues within 5 A of the G4:

| Window | contacts captured (of 283 across 10 templates) |
|---|---|
| round 1 (`A15-60 A240-285`, 92 positions) | 140 (49%) |
| first 40 + last 40 (80 positions) | 251 (89%) |

On the 14 actual round-2 input structures the 40/40 window captures 341 of 357
contacts. Everything it misses is mid-chain barrel surface (residues 48-208),
which should stay fixed anyway.

Consequence worth tracking: for the pKZ169 lineage the windows become residues
1-40 and 266-305, which brings the C-terminal native RHAU half — including the
stacking Trp at position 295 — into the designable set for the first time. Both
stacking tryptophans are now mutable in the same run. That is the intended
symmetry, but it means some round-2 designs will lose both.

---

## 4. Files

| File | Purpose |
|---|---|
| `09_rescore_and_select.py` | corrected chromophore metrics, geometry gate, fitness, Pareto, top-N per template |
| `MISSING_STRUCTURES.md` | the 7 tested variants with no AF3 structure anywhere |
| `09_results/` | output of the partial local run (36 designs only — see below) |

### Running it

```bash
# Sofia: over the full round-1 prediction set on the lab server
python 09_rescore_and_select.py \
    --af3-root  /path/to/output_G4FP_*/03_alphafold3_predictions_bound \
    --apo-root  /path/to/output_G4FP_*/03_alphafold3_predictions_apo \
    --round1-table ../07_results/plot_data_csvs/07_all_designs_scored.csv \
    --validate tested_variants.csv \
    --out 09_results --top-n 5
```

**The committed `09_results/` is a partial run.** Only 36 of the 1000 round-1
designs have stored structures (3 per template, 9 for `r20_l1_cro_mod0`), and
the repo has never contained a single `.cif` — AF3 predictions are excluded by
`designG4FP/.gitignore`. A genuine "top 5 per template" selection needs the full
prediction set and can only be produced on Sofia's server.


---

## 5. Restructured pipeline (after Chen, Wu et al., Nature 2026)

The weighted composite in `09_rescore_and_select.py` is superseded by a staged
gate cascade, `11_cascade.py`, modelled on the metalloprotease pipeline in
doi 10.1038/s41586-026-11167-2. Their Supplementary Methods "Filtering"
thresholds, quoted for reference:

| State | Criteria |
|---|---|
| apo protease | pLDDT > 86, catalytic sidechain RMSD < 1.5 A, Ca RMSD < 1.8 A |
| ES complex | pLDDT > 92, mean interface PAE < 5, min interface PAE < 1.2, catalytic sidechain RMSD < 1.5 A, Ca RMSD < 1.8 A |
| TSA complex | pLDDT > 92, mean interface PAE < 5, catalytic sidechain RMSD < 1.5 A, Ca RMSD < 1.8 A |
| geometry | Zn-chelator distances within 0.5 A of PDB 4qhp; Zn to scissile carbonyl within 0.5 A of the astacin ES model |
| iteration | AF3 TSA models fed back to MPNN for one further round; only designs passing **both before and after** were ordered |

Earlier rounds use looser thresholds that tighten progressively: apo pLDDT
> 82 after the first MPNN round, > 85 after partial diffusion at T=40, > 86 at
the final filter. MPNN was run with `--temperature 0.1 --omit_AA C
--chains_to_design A,B --pack_side_chains 1`, with constrained relax iterated
ten times (their "EnhancedMPNN-FR").

**Their numbers are not transferable.** Those thresholds are calibrated
against a theozyme extracted from a crystal structure; G4FP has no equivalent
ground truth for the intended holo sensor. What transfers is the architecture:
non-compensatory AND-gates per state, cheap criteria first, progressive
tightening, a self-consistency iteration, and the mechanistically essential
species present in the model.

### 5.1 The chromophore is our catalytic water

Modelling Zn(II) alone gave 3/89 active designs and could not tell active from
inactive apart; including the catalytic water gave 31/38. For G4FP the
equivalent omission is the mature chromophore. AF3 predicts the Thr-Tyr-Gly
triad as an ordinary extended tripeptide: across the 14 round-2 input
structures the Thr-Ca to Gly-Ca distance is 4.66-5.35 A, against 3.10 A in the
cyclized CRO of the templates. `10_graft_chromophore.py` restores it.

Measured consequence on the 14 tested variants, ranking by the clearance
available to the grafted chromophore (`cro_min_dist`), the two fluorescent
variants came 1st and 2nd (3.00 and 2.85 A; the twelve dark spanning 2.56 down
to 1.99 A).

**That result is an artifact and is retracted.** `cro_min_dist` measures how
far the *transplanted* CRO ends up from the model's own atoms, which depends
on where the superposition puts it. Within the 14 tested variants it
correlates with the superposition residual `graft_fit_rmsd` at r = 0.849
(p = 0.0001), and the two bright designs are the only ones whose barrels
deviate appreciably from their template around the chromophore (fit 0.98 and
1.26 A, against 0.16-0.32 A for all twelve dark). Their apparent extra
clearance is imprecise placement, not a roomier pocket.

Measured directly in the model with no superposition at all, the pocket is the
same size in both classes:

| pocket measure (triad excluded) | bright (n=2) | dark (n=12) |
|---|---|---|
| protein atoms within 6 A of triad centroid | 23.0 | 23.3 |
| protein atoms within 8 A | 63.0 | 64.3 |
| mean distance to the nearest 20 atoms | 4.93 A | 4.93 A |
| closest environment atom to the triad | 2.76 A | 2.53 A |

On the last measure the bright designs rank 4th and 5th of 14, inside the dark
distribution. **The dark designs are not dark because they compress the
chromophore pocket.** No chromophore-geometry measure examined separates the
two classes.

What survives: restoring CRO is still correct, because AF3 genuinely omits the
mature chromophore and LigandMPNN loses that context without it (5 of the 80
designed positions lie within the 8 A ligand cutoff of the chromophore, since
circular permutation places the cpGFP termini against it). But the graft is an
input-preparation step, not a discriminator.

Consequence for the cascade: the chromophore gate is **off by default**. The
2-of-2 recall reported below comes entirely from stage B -- disabling stage A
leaves 8 survivors of 36, still containing both bright variants and none of
the twelve dark ones. The CRO clash count remains available as an opt-in
criterion, but only where `graft_fit_rmsd` is small enough for the placement
to be trusted: `r20_l2_cro_mod3` designs clash 11-15 times at fit 0.25-0.31 A,
which is interpretable; a clash count at fit > 0.5 A is not.

### 5.2 Cascade stages

| Stage | Question | Needs AF3? |
|---|---|---|
| A chromophore | can the barrel host a mature chromophore? | no |
| B apo foldability | can it fold and mature at all? | apo prediction |
| C holo complex | does it engage the G4 with a credible interface? | holo prediction |
| D self-consistency | does LigandMPNN on the AF3 model recover the sequence? | both |

Run on the 36 designs with stored structures: 36 -> 12 (A) -> 4 (B) -> 4 (C).
Of the 14 tested variants among them, exactly 2 survive and both are the
fluorescent ones. The other two survivors are untested `des1_cro_mod0`
designs -- the template where 100/100 designs clear the apo floor and which
round 1 never sampled.

### 5.3 Selection: stratify, do not maximise

Round 1 only ever tested designs below their template's apo order. The upper
edge of the viable window is uncharacterised, and pKZ168 suggests that far
above the floor you get ligand-independent folding. `11_cascade.py --select
stratify` samples bands of `apo_plddt_vs_template` instead of maximising it,
so round 2 locates that edge rather than assuming it.

### 5.4 What we are not copying

Their hit rates rest on 135 designs tested through a deliberately cheap
primary screen (a cis fusion of protease and substrate in one gene), with a
stricter trans screen on the 22 hits; the two correlate at Pearson r = 0.91.
G4FP tested 21. With two positives from a single template, no computational
metric can be validated on round-1 data -- the binding constraint is the size
of the experimental screen, not the sophistication of the score.

# Rescoring round 1 with confidence-contrast gates

**Status:** gate implemented (`13_delta_ipae_gate.py`), run over all 999
apo/holo pairs on the lab server. Results and what they decide are at the end.

---

## 1. What changes, and why

The round-2 SWITCH gate asks whether the chromophore *moves* between the apo
and holo models. `14_KNOWN_LIMITATIONS.md` records the answer it gets:
**0.19–0.22 Å, at a superposition fit residual of 0.16–0.24 Å.** The signal is
smaller than the error of the operation that produces it. That gate is ranking
noise, and tightening its percentile only ranks noise more confidently.

Chang & Polizzi (bioRxiv 2026, [10.64898/2026.08.03.742366](https://doi.org/10.64898/2026.08.03.742366))
designed a working ligand-driven conformational switch — a lid that closes on
exatecan — and did not gate on a structural difference at all. They gated on a
**confidence contrast**: keep designs the predictor is *confident* about with
ligand and *uncertain* about without it. 24 designs tested, 4 clear switches,
leads at 200 pM and 510 pM. Hits had ΔiPAE > 3.

Translated to this campaign:

| | round 2 as written | this gate |
|---|---|---|
| quantity | chromophore apo–holo RMSD | mean PAE, module ↔ core, apo minus holo |
| needs superposition | **yes** | **no** |
| can encode its own fit error | yes — and does | structurally cannot |
| dynamic range observed | 0.19–0.22 Å | see §5 |
| new predictions required | none | none |

The last row is the reason to run it now: **both PAE matrices are already in
the `*_confidences.json` files AlphaFold3 wrote in March.** Nothing is
re-predicted. The gate is arithmetic over data on disk.

## 2. Definitions

```
MODULE  chain A residues 15–60 ∪ 240–285
        = ligandmpnn.design_residues from pipeline_config.yaml,
          identical across all ten templates (verified, not assumed)

CORE    the rest of chain A, minus ±4 residues around each module boundary
        (PAE between sequence neighbours is low whatever the fold does;
         without the buffer the contrast is diluted by trivial pairs)

iPAE(state) = mean PAE over all MODULE↔CORE token pairs, both directions
              (PAE is asymmetric, so both blocks are averaged)

ΔiPAE       = iPAE(apo) − iPAE(holo)
```

A design that switches is floppy at that interface without G4 and pinned with
it, so ΔiPAE is large and positive. A design rigid in both states scores ~0 —
which is what round 1's composite selected for, since every term in it rewarded
confidence.

**Cross-check, free:** the same quantity computed over a *contact-derived*
module (AF3 `contact_probs` ≥ 0.5 against the G4 chain) instead of the design
window. If the two disagree, the result depends on where the window was drawn
and neither should be trusted.

## 3. The gate order

Non-compensatory, cheapest and most decisive first. Nothing is averaged into
anything else.

| # | Gate | Quantity | Cost | New? |
|---|---|---|---|---|
| 1 | **FOLD** | holo pLDDT, pTM, chromophore pLDDT (corrected localiser) | free | no |
| 2 | **APO window** | apo pLDDT inside a band, not above a floor | free | no |
| 3 | **BIND** | protein↔G4 interface PAE, holo | free | recomputed here |
| 4 | **SWITCH** | **ΔiPAE > threshold** | free | **yes** |
| 5 | **NEGATIVE** | reject designs confident in *both* states | free | **yes** |

Gate 5 is the one round 1 had no analogue for. Chang & Polizzi's lid carries a
designed Glu that would bury an uncompensated charge if it closed without
ligand; nothing in round 1 penalised a design that is bright without G4. In
rescoring it is the cheap half of that idea: a design with low iPAE in *both*
states is pre-closed and is removed, not ranked lower.

## 4. Thresholds — percentile, not absolute

Chang & Polizzi's numbers (iPAE < 5 bound, ΔiPAE > 3) were calibrated on two
protein domains tethered by a GGS linker around a 50 Å² small molecule. This
system is one chain against a 16-nt nucleic acid. **Do not transplant the
cutoffs.** Use percentiles of the measured round-1 distribution for iteration
1, and tighten over iterations as round 2's plan already does. Report the
absolute value each percentile corresponds to, so the two systems can be
compared later.

## 5. What the rescoring can and cannot decide

It answers exactly one question: **does this architecture have a switch to
find?**

If ΔiPAE is flat across all 999 designs, the barrel-plus-grafted-module
architecture does not switch in AF3's view, and no re-gating of that population
will produce one. That is an argument for rebuilding rather than re-ranking,
and it is worth far more than a corrected top-10.

It cannot validate the gate. The gate predicts **switching**. Of 21 tested, 2
were fluorescent and **0 were G4-responsive** — the switch phenotype has zero
positives, so there is no label to test against. Even for the fluorescence
label, an exact permutation test over C(21,2) = 210 arrangements has a floor of
p = 0.0048; no metric can reach significance at that sample size. The gate is
adopted on mechanism and on a prospective published result, and validated in
round 2 by testing designs chosen *because* of it.

## 6. Provenance

`13_delta_ipae.csv` carries its generating script, git commit, root path, the
module and core definitions, the buffer width, the contact threshold, and the
row and failure counts in its header comments. `13_delta_ipae_gate.py` prints a
per-column participation report (non-empty / distinct / standard deviation) on
every run; any column with zero variance is flagged as not participating,
whatever weight a downstream score gives it.

## 7. Results — run over all 999 apo/holo pairs, 17 November 2026

`13_delta_ipae_module.csv`, 999 rows, 0 failures. One design in
`output_G4FP_r28_l1_cro_mod2` has a holo prediction with no matching apo and is
excluded; every other pair scored.

### 7.0 A correction, recorded rather than deleted

The first run used `ligandmpnn.design_residues` (A15–60 ∪ A240–285) as MODULE.
**That is the wrong region.** The collaborator's
`13_fixed_core_definition.csv` records where the grafted RHAU DSM module
actually sits — on most templates residues 1–28 ∪ 279–306, with the barrel at
29–278. The design window overlaps the module only at its edges and is
otherwise barrel, so the first run measured **barrel against barrel**: iPAE
~10 Å in both states, ΔiPAE median −0.55, 3/999 clearing the threshold. That
result is an artifact of the region definition, not a property of the designs,
and it is superseded by everything below. It is left on the record because it
is exactly the failure the claim test exists to catch — a quantity correctly
computed over the wrong atoms, carrying a plausible name.

The corrected run uses `--module-map`, which reads the per-template module and
barrel boundaries from that file. Two templates (`r20_l1_cro_mod0`,
`r20_l2_cro_mod3`) have designs whose module position varies; the modal
definition is used and the script prints a NOTE for each. Treat those two
templates' rows as provisional.

### 7.1 There is a switch signal, and it is large

| | design-window module *(wrong)* | RHAU module vs barrel *(correct)* |
|---|---|---|
| ΔiPAE median | −0.55 Å | **−0.04 Å** |
| interquartile width | 1.26 Å | **3.28 Å** |
| p95 | +1.21 Å | **+5.64 Å** |
| max | +4.11 Å | **+11.68 Å** |
| clearing the published ΔiPAE > 3 | 3 / 999 (0.3%) | **168 / 999 (16.8%)** |

For comparison, the SWITCH gate this replaces has a dynamic range of
0.19–0.22 Å at a fit residual of 0.16–0.24 Å. This one has an interquartile
width of 3.28 Å and needs no fit at all.

The module is genuinely loose relative to the barrel — median iPAE 21.2 Å apo,
19.8 Å holo — which is what a terminal graft on a flexible linker should look
like, and is the precondition for a lid-style switch.

### 7.2 Nearly a third of the library is pre-closed

| apo iPAE threshold | designs below it |
|---|---|
| < 6 Å | 100 (10.0%) |
| < 8 Å | 215 (21.5%) |
| < 12 Å | **300 (30.0%)** |
| < 15 Å | 300 (30.0%) |

The distribution is sharply bimodal with a clean gap between 12 and 15 Å:
**30% of designs have the module already packed against the barrel without any
G4 present**, and 70% have it loose. A pre-closed design cannot switch — there
is nothing for the ligand to close. This is precisely the negative-design
failure Chang & Polizzi designed against with their buried-charge Glu, and
round 1 had no gate for it. Removing those 300 is free.

### 7.3 For most designs the interface does not exist

| | value |
|---|---|
| protein↔G4 interface PAE, median | **20.3 Å** |
| designs with ≤ 2 confident protein–G4 contacts (`contact_probs` ≥ 0.5) | **701 / 999 (70%)** |
| designs with **zero** such contacts | **152 / 999 (15%)** |

Round 1 ran AlphaFold3 on all 1,000 designs and removed none of them. A gate
asking only "does the model think this protein touches the G4 at all" would
have removed 70% of the library before the composite score ever ran.

### 7.4 The cascade on the existing population

| Step | Gate | Survivors |
|---|---|---:|
| 0 | start | 999 |
| 1 | **CONTACT** ≥ 3 predicted protein–G4 contacts | 298 |
| 2 | **NEGATIVE** apo iPAE ≥ 8 Å (module not pre-closed) | 197 |
| 3 | **BIND** interface PAE ≤ P40 (19.7 Å) | 101 |
| 4 | **SWITCH** ΔiPAE > 3 Å | **71** |

Survivors by template: `des1_cro_mod4` 47, `r28_l2_cro_mod2` 9,
`r28_l3_cro_mod2` 9, `r28_l2_cro_mod3` 3, `r28_l3_cro_mod0` 2,
`r28_l1_cro_mod2` 1. Better spread than the artifact run gave, but still
template-dominated — apply percentiles **within template** before treating any
of the 71 as a design-level result.

Best few, all four gates passed:

| design | ΔiPAE | interface PAE | apo iPAE | contacts |
|---|---:|---:|---:|---:|
| `r28_l3_cro_mod2_design_38` | +11.68 | 11.35 | 24.62 | 14 |
| `r28_l3_cro_mod2_design_25` | +11.11 | 8.54 | 23.00 | 19 |
| `r28_l3_cro_mod2_design_98` | +10.75 | 10.27 | 24.31 | 17 |
| `r28_l2_cro_mod2_design_40` | +8.93 | 8.75 | 22.77 | 9 |

### 7.5 The 21 tested designs cannot settle this

All 21 matched. Both fluorescent designs sit mid-distribution —
`r28_l2_cro_mod2_design_0072` (the lead) at ΔiPAE −0.30, 47th percentile, and
`design_0016` at −0.44, 44th. Dark mean +3.28, fluorescent mean −0.37.

That looks bad for the gate until you ask what the 19 dark designs tell you
about switching: **nothing at all.** A protein that never fluoresces cannot be
scored for whether its fluorescence changes. 19 of 21 were dark, so the switch
phenotype was unobservable in 90% of the experiment, and 0 of 21 were
G4-responsive. There are no positives for the phenotype this gate predicts.

`r28_l3_cro_mod2_design_0025` is the sharpest case: it passes every gate here
(ΔiPAE 11.1, 99.8th percentile; interface PAE 8.5, 2.7th) and it was dark. That
is a chromophore-maturation failure, not a switching failure, and the two are
not separable in this dataset.

Even for the fluorescence label, an exact permutation test over C(21,2) = 210
arrangements has a floor of p = 0.0048; no metric can reach significance at
this sample size.

### 7.6 What this decides

1. **Adopt the gate**, with the module map, not the design window. It costs
   nothing, has 15× the dynamic range of the metric it replaces, and cannot
   encode its own fit error.
2. **Add CONTACT and NEGATIVE floors to every cascade.** Together they remove
   ~80% of the population for free, on grounds no reasonable reviewer would
   argue with: no predicted contact with the ligand, or already closed without
   it.
3. **Re-gating is worth running — it yields 71 candidates the composite never
   surfaced — but it is not the round-2 plan.** The binding constraint is not
   the filter. It is that 19 of 21 designs never turned on, which made the
   switch untestable. Until the maturation failure mode is removed, no switch
   gate can be validated, however good it is. That is the argument for the
   non-covalent fluorophore in Option E: it deletes the step that made 90% of
   round 1 uninformative.

# Backbone diversification without touching the fluorescent core

## The constraint, measured

AlphaFold3 predicts a **closed, sealed beta-barrel in every G4FP model we have** - apo
and holo alike.

Chromophore triad solvent-accessible surface area, protein chain only, mean of 5 AF3
models per state:

| construct | apo (A^2) | holo (A^2) | apo - holo |
|---|---|---|---|
| pKZ168 | 9.5 +- 2.6 | 9.2 +- 3.1 | +0.3 |
| pKZ169 | 8.3 +- 1.3 | 11.6 +- 2.1 | -3.3 |
| pKZ170 | 7.7 +- 2.7 | 8.1 +- 3.4 | -0.4 |

That is 1.4-2.2% of the triad's theoretical maximum exposure. Across the 36 grafted
round-2 inputs the mature CRO ranges 1.2-30.3 A^2 (median 8.3) - sealed in all of them.

Two consequences:

1. **There is no modelled OFF state.** The differences above are inside the
   model-to-model spread. AF3 returns a folded barrel whatever we give it, so the
   SWITCH gate in `12_iterate.py` is scoring prediction *confidence*, not an open
   barrel. No scoring change fixes this - the conformation we are selecting for is
   not present in the data being scored.
2. **The chromophore shell is buried, not a gate.** Of the 59 residues within 8 A of
   the chromophore in pKZ169, 46 are buried and 13 are exposed - and those 13 sit on
   the *outer* barrel surface, not lining an opening. Redesigning them buys no
   transduction path and risks the barrel. This retracts the suggestion to widen the
   design window into the barrel.

So the switch has to be built into the geometry, not selected for after the fact. The
place to build it is everywhere except the barrel.

## What is fixed and what is free

Per-input boundaries are in `13_fixed_core_definition.csv`. For the pKZ169 family
(`r28_l2_cro_mod2`, 302 residues):

| segment | residues | treatment |
|---|---|---|
| RHAU-N | 1-28 | rigid body, internal structure preserved |
| linker N | (none in this template) | buildable |
| **cpEGFP barrel + CRO** | **29-277** | **frozen, zero internal change** |
| linker C | - | buildable |
| RHAU-C | 278-305 | rigid body, preserved (stacking Trp295) |

All 46 buried chromophore-shell residues lie inside 29-277.

Note: for the `r20_*` templates the N-terminal boundary in the CSV is approximate - it
is found by aligning the RHAU28 motif, and those templates carry RHAU20, so the match
slides. Set those boundaries by hand before use.

## Three options, in increasing risk

### Tier 0 - linker enumeration. No diffusion, no risk, available now.

Vary the length and composition of the two linkers only. In GCaMP-lineage sensors the
linkers are the historical tuning variable (the GCaMP6 `LKI`/`LP` pair came from
exactly this). Changing linker length changes where the RHAU domains sit relative to
the barrel, which is the degree of freedom that sets dynamic range. Implement with
LigandMPNN `--fixed_residues` covering the barrel, one job per linker length.

### Tier 1 - multi-motif scaffolding. The principled route.

Hold three motifs and rebuild only what joins them:

    contigmap.contigs=[A1-28/2-8/A29-277/2-8/A278-305]

Separate contig blocks are not held in fixed relative geometry, so the two RHAU
domains are repositioned around the barrel while every atom that determines
fluorescence is carried through untouched.

Because the barrel is rigid, the chromophore is protected by construction - you do not
strictly need all-atom diffusion to keep it intact. You do want the G4 in context for
sensible binder placement, which argues for RFdiffusionAA over vanilla RFdiffusion.

### Tier 2 - partial diffusion. Cheapest to set up, highest risk.

`diffuser.partial_T` noises the *whole* input, barrel included. Only acceptable at low
T and behind the guard filters below. Prefer Tier 1 if RFdiffusionAA is available.

## Guard filters - the backbone-filter box

Apply after any diversification, before sequence design:

1. Barrel motif Ca RMSD to parent < 0.5 A.
2. The buried chromophore-shell residues (column `shell_buried`) unchanged in
   sequence - enforce via LigandMPNN `--fixed_residues`, then verify.
3. Chromophore triad relative SASA < 5%. Measured baseline 1.4-2.2%.
4. Pre-cyclization Thr-C to Gly-N distance inside the established range.

Filters 1 and 3 are what make diversification safe; neither exists in the pipeline yet.

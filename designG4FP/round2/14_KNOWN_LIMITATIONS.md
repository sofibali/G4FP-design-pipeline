# Known limitations of the round-2 pipeline

Recorded so they are not rediscovered later. None of these is a blocker; all are
bounds on what the pipeline can be expected to deliver.

## 1. AlphaFold3 never predicts the OFF state

The chromophore is sealed in every model we have, apo and holo alike. Triad
solvent-accessible surface area, protein chain only, mean of 5 AF3 models:

| construct | apo (A^2) | holo (A^2) | apo - holo |
|---|---|---|---|
| pKZ168 | 9.5 +- 2.6 | 9.2 +- 3.1 | +0.3 |
| pKZ169 | 8.3 +- 1.3 | 11.6 +- 2.1 | -3.3 |
| pKZ170 | 7.7 +- 2.7 | 8.1 +- 3.4 | -0.4 |

1.4-2.2% of maximum exposure, and the differences sit inside the model-to-model
spread. Across the 36 grafted round-2 inputs the mature CRO spans 1.2-30.3 A^2
(median 8.3) - buried in all of them.

The SWITCH gate therefore scores a *confidence* change and a small triad
displacement, not an opening barrel. Measured on the three tested constructs the
chromophore-local RMSD between apo and holo is 0.19-0.22 A (local superposition
residual 0.16-0.24 A, so the numbers are interpretable). That is near zero: there
is no chromophore-local conformational switch in any construct tested so far.

## 2. The chromophore shell is not designable, and should not be

Of the 59 residues within 8 A of the chromophore in pKZ169, 46 are buried and 13
are solvent-exposed - and those 13 sit on the OUTER barrel surface, not lining an
opening. Widening the design window into the barrel would risk fluorescence
without buying a transduction path. Per-input boundaries and the buried shell
residue list are in `13_fixed_core_definition.csv`.

Backbone diversification (partial diffusion / motif scaffolding) was considered
as a way around this and deliberately not adopted for round 2.

## 3. Interface confidence is low across the board

Mean protein-DNA interface PAE on the tested constructs is 13.2-16.7 A
(min 2.0-3.3 A). For reference the metalloprotease pipeline gated at mean
interface PAE < 5. Combined with the steric clashes present in every holo model
(10-37 atom pairs below 2.5 A), the BIND gate is ranking interfaces that AF3
itself is not confident about. It is used as a relative filter, not as evidence
that any design binds.

## 4. Two positives from one template

Round 1 returned 2 fluorescent variants out of 21, both from `r28_l2_cro_mod2`.
No gate threshold can be validated on that, which is why the thresholds are
percentiles of the current population rather than absolute values transferred
from the literature.

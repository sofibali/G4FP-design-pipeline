# G4-binding protein modules — what exists, and what it implies for the switch

Literature retrieved from PubMed; structures from the RCSB PDB. Measurements in
§3 were made here from the deposited coordinates.

---

## 1. The candidates

| Module | Size | Target | Selectivity | Structure | Works in cells |
|---|---|---|---|---|---|
| **Nb55** | nanobody, ~123 aa | TBA G4 | **selective** — of 15 sequences tested, binds only AGCGA and TBA | **9GXH** (bound, 1.90 Å) + **9GV4** (free, 1.53 Å) | not reported |
| **SG4** | nanobody, ~15 kDa | raised against human **MYC** G4 | **pan-G4** — low nM across a wide range of folded G4s | AlphaFold + MD model only | **yes** — expressed intracellularly, images endogenous G4s |
| **BG4** | scFv, 27 kDa | pan-G4 | pan-G4 | — | no intracellular expression; poor in vitro yield |
| **G4P** | artificial, 6.7 kDa | pan-G4 | pan-G4 | — | yes (ChIP-seq in live cells) |
| **RHAU DSM** | peptide | parallel G4s | topology-biased | several | yes |

Sources: Pevec *et al.*, *Nucleic Acids Research* 2025,
[10.1093/nar/gkaf453](https://doi.org/10.1093/nar/gkaf453) (Nb55);
Galli *et al.*, *JACS* 2022, [10.1021/jacs.2c10656](https://doi.org/10.1021/jacs.2c10656) (SG4);
Zheng *et al.*, *Nucleic Acids Research* 2020,
[10.1093/nar/gkaa841](https://doi.org/10.1093/nar/gkaa841) (G4P).

**Nb55 sequence** (from 9GV4, `_entity_poly`, C-terminal His6 included):

```
QVQLQESGGGLVQAGGSLRLSCAASGSRFSSNTMTWYRQAPGKQREWVATMRSIGTTRYASSVEGRFTLSRDNAKNTVYL
QMNSLKPEDTAVYYCNLRRGGGIYWGQGTQVTVSSHHHHHH
```

### Availability

- **SG4** — Addgene **#196071** (`pHEN2-SG4`, wild type) and **#196072**
  (`pHEN2-SG4 R105A`), both deposited by the Balasubramanian lab, academic /
  non-profit only. The R105A mutant is the paper's own loss-of-binding control:
  building the same sensor on R105A gives a matched dead-binder negative
  control, which is the cleanest way to show a fluorescence change is actually
  G4-driven rather than an expression or folding artifact. The SG4 sequence is
  in Figure 1C of the paper, which is open access at PMC9782783.
- **Nb55** — sequence above, taken from the 9GV4 deposition. No plasmid
  distribution found; it would be synthesised.

## 2. Which one for which job

- **If the target is c-MYC and pan-G4 binding is acceptable → SG4.** It was
  raised against the MYC G4, it is the only one of these demonstrated to fold
  and function inside a cell, and the authors mapped the key contacts to
  specific arginines and confirmed them by Arg→Ala mutagenesis. The cost is
  that it binds everything G4-shaped, so it cannot drive the discrimination
  gate that Option D is built around.
- **If the sensor must discriminate topologies → Nb55.** It is the only
  protein module here with demonstrated *sequence* selectivity among G4s, and
  the only one with a crystal structure of the complex. Its paratope is an
  unusual one: a concave framework β-sheet with two β-hairpins and very little
  CDR involvement, mimicking thrombin's own interface with TBA by isosteric
  substitution at key positions. A rigid β-sheet face is a far better substrate
  for vdM-style interface design than CDR loops are.
- **If size is the constraint → G4P**, at 6.7 kDa, roughly half a nanobody.

## 3. The measurement that matters: does a G4 module switch?

9GXH (bound) and 9GV4 (free) are an apo/holo crystallographic pair of the same
nanobody, so they answer directly what conformational change G4 binding
actually produces. Measured here, Cα only, with the two independent copies in
the 9GXH asymmetric unit used as a noise floor:

| Quantity | Value |
|---|---|
| global Cα RMSD, apo vs holo copy A | **0.83 Å** |
| global Cα RMSD, apo vs holo copy B | **1.10 Å** |
| global Cα RMSD, holo copy A vs copy B *(same crystal, both bound — noise)* | **0.74 Å** |
| median per-residue deviation, apo vs holo | **0.60 Å** |
| median per-residue deviation, copy A vs copy B *(noise)* | **0.43 Å** |
| residues moving > 1.0 Å **and** > 2× the copy-to-copy noise | **8 of 107** |

**Globally, Nb55 does not switch.** 0.83 Å apo–holo against 0.74 Å between two
copies of the *same* bound state is a ratio of 1.1. A selective, structurally
characterised, high-affinity G4-binding protein is, in bulk, rigid.

**Locally, it switches hard.** The motion is concentrated in one short
stretch:

| Residue | apo vs holo A | apo vs holo B | copy A vs copy B (noise) | signal / noise |
|---:|---:|---:|---:|---:|
| 42 | 3.85 Å | 3.79 Å | 0.50 Å | **7.6** |
| 43 | 2.09 Å | 2.42 Å | 0.48 Å | 4.7 |
| 44 | 1.56 Å | 1.99 Å | 0.55 Å | 3.2 |
| 41 | 2.33 Å | 3.17 Å | 0.97 Å | 2.8 |
| 107 | 1.28 Å | 1.11 Å | 0.29 Å | 4.1 |

Residues 41–44 are `A-P-G-K` in framework 2, immediately before CDR2 —
consistent with the paper's finding that the paratope is scaffold-based rather
than CDR-based.

## 4. What this decides

1. **Grafting any of these modules will not give an induced-fit switch.** The
   module's own internal rearrangement on G4 binding is ~0.6 Å, inside the same
   noise band as the 0.19–0.22 Å the round-2 SWITCH gate measures. This is the
   same conclusion already reached for RHAU DSM, now with a number attached and
   a crystallographic noise floor to compare it against.

2. **So the switch has to be rigid-body domain reorientation, not local
   rearrangement** — the module moving *as a unit* relative to the fluorescent
   body. That is the lid mechanism, and it is invisible to any
   chromophore-local RMSD. It is what `13_delta_ipae_gate.py` measures.

3. **But there is a real 4 Å signal if you couple to the right eight
   residues.** The usable amplitude is not spread over the module; it sits in
   one loop. A design that couples the chromophore or dye pocket to residues
   41–44 of an Nb55-derived module has genuine mechanical range to work with.
   A design that couples to the module's centre of mass has none. For SG4 the
   equivalent loop is unknown, because there is no experimental structure —
   which is an argument for Nb55 as the design template even if SG4 is the
   better reagent.

4. **Neither module is a drop-in.** Nb55 binds TBA and AGCGA, not MYC; SG4
   binds MYC but binds everything else too. Choosing between them is choosing
   between selectivity you have to retarget and targeting you have to make
   selective.

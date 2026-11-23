# Round 2: a Rap1 hinge carrying a disorder-to-order insertion loop

Written 23 November 2026. Supersedes the module choice in `17_ROUND2_DESIGN.md`
§3; the cost model and campaign framing there still stand.

Literature retrieved from PubMed. Measurements marked *measured here* were
computed in this directory from deposited coordinates.

---

## 0. The objection this answers

Two objections, both correct, and they kill the previous plan together:

1. **A single grafted domain leaves nothing to act on the barrel.** Round 1
   grafted a module at residues 1-28 + 279-306 and hoped the barrel would
   notice. Nothing in the construct was responsible for moving a chromophore
   contact.
2. **Asking a generator for two opposing interfaces produces a rigid binder.**
   An optimiser that must satisfy an apo interface and a holo interface at once
   finds a backbone that satisfies both without moving. This is not a risk, it
   is the expected outcome -- `16_G4_MODULE_SURVEY.md` §3 shows every G4 binder
   with an apo/holo pair sits at 0.75-1.02 A, two unrelated scaffolds, two
   independent selection campaigns.

So the switch cannot come from the binder's own backbone, and it cannot be
asked for as a constraint. It has to be borrowed from something that already
moves, and coupled through a mechanism that is known to work.

## 1. What NCaMP7 actually shows

Subach et al., *Int J Mol Sci* 2020, 10.3390/ijms21051644, **PDB 6XW2** at
1.75 A. Read in full via PubMed Central 23 Nov 2026.

Two facts that change the architecture:

**It is an insertion, not a circular permutation.** CaM/M13 sits inside
mNeonGreen **between residues 145 and 146**. Both FP termini stay free. This is
the NTnC-like design, and it is *not* what G4-Flame does (terminal fusion on a
cpYFP) nor what round 1 did.

**Transmission is one side chain.** The chromophore (68-GYG-70) sits on the
central helix, fixed by P65, R98, Y114, E401. What reaches it from the sensing
domain is **Y225 -- a single tyrosine at the end of CaM helix alpha-4 -- which
protrudes into a surface "hole" in the barrel** and makes a water-mediated
H-bond to the deprotonated phenolate. In GCaMP6m (3WLD) the equivalent pair is
R376/Y380; Y380 makes no H-bond at all and simply blocks solvent access with a
bulky phenol.

> The entire allosteric transmission in both sensor families is one aromatic
> side chain entering one hole.

That is a far smaller target than "couple a domain to a barrel", and it is the
reason this architecture is worth attempting at all.

**How they found it: screening, not design.** Two libraries (M13/CaM and
CaM/M13 order) at the 145/146 site with **both 3-residue linkers fully
randomised**, screened on bacterial colonies. The final sensors carry **25-26
mutations**: 11 in the FP, 8-9 in CaM/M13, **6 in the linkers**. Delta-F/F 27,
Kd 125 nM.

*Caution:* the paper never uses the word "allosteric", and hedges the
mechanism -- "we speculate that Y225 ... plays a crucial role in the
translation of the Ca2+-dependent conformational change". Treat the Y225 story
as a well-supported structural interpretation, not a measured mechanism.

## 2. The construct

| part | residues | role |
|---|---|---|
| Block 1 | Rap1 361-440, incl. Myb1 recognition helix 398-412 | rigid, anchored to the barrel |
| **Hinge** | **440-460** | designable; source of the 104 deg rotation |
| Block 2 | Rap1 501-573 | moves as a unit; centroid 14.7 A, residues to 35.7 A |
| Insertion loop | de novo, or the RHAU DSM | spans block 2 -> barrel wall; slack in apo, taut in holo |
| The tyrosine | one position in the loop | enters the chromophore hole; vdM-anchored |

The loop does the job CaM's alpha-4 helix does in NCaMP7. The hinge does the job
Ca2+ binding does in CaM. Neither half is designed from scratch.

**Why Rap1 and nothing else.** It is the only natural G4-binding domain with
real amplitude, and the amplitude is inter-domain rather than internal.
*Measured here* (`16_G4_MODULE_SURVEY.md` §4): sliding 40-residue windows give
5.99 A internal RMSD at 421-460 and under 0.9 A either side, so the motion is
in one place. The recognition helix itself is rigid at **0.46 A** -- it is
*repositioned*, not deformed. Global Ca RMSD 11.33 A over 175 residues.

## 3. What to do with RHAU

Three options were considered. One is rejected on published evidence.

### (a) RHAU DSM replaces a Rap1 helix -- REJECT

Isolated Rap1 recognition-helix peptides bind G4 at 1.2-5.4 uM, **40-200x
weaker than the folded DBD**. The fold is load-bearing, and that coupling is
precisely what makes Rap1 worth using. Substituting a peptide that *does* work
alone (the DSM does) removes the dependence of binding on the fold -- and with
it, any reason for binding to move the hinge. The product would be RHAU binding
with a Rap1-shaped decoration.

### (b) RHAU DSM becomes the insertion loop -- RUN

The DSM is natively disordered and **folds into an alpha-helix on binding**
(`16_G4_MODULE_SURVEY.md` §6). That is exactly the disorder-to-order element
the architecture needs, and it exists in nature instead of being designed.

**Round 1's real error, stated plainly: the DSM was modelled pre-folded.** The
folding transition was spent before the sensor ever saw a G4. See
`g4fp-round2-binding-module` and `13_RESCORING_PLAN.md` §130.

*The geometric hypothesis that makes it work.* Rap1's Myb1 hydrophobic patch
(Ser402/His405/Val409) docks flat on the **3' G-tetrad**. The DSM end-stacks on
a terminal tetrad. A 16-nt parallel G4 has exactly **two** faces, so:

- DSM on the **5' face**, Rap1 on the **3' face** -> compatible, and the two
  modules brace the same G4 from opposite ends.
- DSM on the 3' face -> direct competition, and the construct is worse than
  either module alone.

This is a hard constraint, not a preference, and it is testable in AF3 before
anything is built. **Make it a pass/fail gate, not a score.**

### (c) No RHAU -- RUN

Rap1 plus a de novo loop. Simplest, avoids the G4-Flame prior art entirely
(10.1093/nar/gkag179), and makes the loop a pure design problem. (b) and (c)
share every input except the loop, so run them in parallel.

## 4. How to express this to the pipeline

**The hinge cannot be given to a generator as a flexibility parameter.** It can
be given as two motif-placement problems that use the same two motifs.

1. **Two conditioning states.**
   - Motif A: block 1 + recognition helix, coordinates from the structure.
   - Motif B: block 2, coordinates from the structure.
   - Run once with A and B in the **6LDM** (G4-bound) relative orientation.
   - Run once with A and B in the **1IGN** (duplex-bound) relative orientation.
   - Keep scaffolds compatible with **both** placements.

   The amplitude lives in the difference between the two runs. It is never
   expressed inside a single input, which is why the previous "two interfaces"
   framing could not work.

2. **Hinge left unconstrained.** Residues 440-460 are *not* motif. The
   generator's job there is a linker that tolerates both orientations without
   preferring either.

3. **No nucleic acid at generation.** RFD4-Proteina's DNA conditioning is
   duplex-only -- G-quadruplex, ssDNA and non-duplex folds appear nowhere in its
   README, cookbook, input-prep guide or model card, and its DNA benchmark is 3
   targets on a previous release (see `17_ROUND2_DESIGN.md` §4.1). So this is
   motif scaffolding plus a small-molecule pocket, which are its validated
   modes. **The G4 enters only at AF3 scoring.**

4. **The loop is a second motif problem.** Anchor one end on block 2 and one end
   in the barrel wall at the NCaMP7-equivalent insertion site; let the generator
   build the span. For option (b) the DSM helix is itself a motif, placed with
   its stacking face toward the 5' tetrad.

5. **The tyrosine is placed, not designed.** Use a vdM for the Tyr-OH ->
   chromophore-phenolate contact (`npolizzi/Combs`). Note the published limit:
   vdMs cover protein-small-molecule and protein-protein contacts only, so the
   G4 interface cannot be vdM-anchored -- only this one contact can.

## 5. What to score

Four quantities, and **the cascade must stay non-compensatory**. Round 1's
lesson was that a weighted score lets a strong term mask a failed one; the
resulting ranking scored below a coin flip against experiment (AUC 0.158, n=21,
`19_ROUND1_RESCORE_RESULTS.md`).

| quantity | measures | status |
|---|---|---|
| Delta-iPAE, apo minus holo | rigid-body switch | built, validated |
| **Delta-pLDDT on the loop residues** | **disorder -> order** | **new, free** |
| Tyr-OH to phenolate distance, apo vs holo | transmission | new, geometric |
| DSM face assignment (5' vs 3') | competition | new, pass/fail |

**Delta-pLDDT on the loop costs nothing and should be added regardless of which
option wins.** AF3 already writes per-residue pLDDT in both states, and the
round-1 tree has it for all 999 designs. A loop that is genuinely
disorder-to-order shows low apo pLDDT and high holo pLDDT *on the same
residues*. This is the one gate round 1 never had and could have had for free.

Keep using Delta-iPAE for the rigid-body half: it needs no superposition, so it
cannot encode its own fit error -- the defect that left the round-1 RMSD switch
gate dead at 0.19-0.22 A.

## 6. The first experiment, before any design

**An AF3 apo ensemble of the Rap1 DBD.** There is no apo structure of the yeast
Rap1 DBD anywhere, so the 104 deg rotation is *duplex-vs-G4*, not *apo-vs-G4*.
Whether this is a clean two-state hinge or a disorder-to-order collapse from a
heterogeneous apo state is unknown, and it decides the architecture:

- clean two-state -> gate on Delta-iPAE, select for reorientation.
- heterogeneous apo -> gate on Delta-pLDDT, select for ordering.

Hours, not weeks, and it reuses the machinery in `13_delta_ipae_gate.py`.

## 7. The screening fallback

If the design cannot place the loop, do it the way NCaMP7 was actually built:
randomise the insertion position and both linkers and screen.

- **DIP-seq** (Nadler et al., *Nat Commun* 2016, 10.1038/ncomms12266) maps
  allosteric hotspots by random transposon insertion + NGS.
- **Sort-seq** (Koberstein et al., *ACS Chem Biol* 2021,
  10.1021/acschembio.1c00423) quantifies brightness and dynamic range per
  variant by FACS + NGS + maximum likelihood, and trains a model that predicts
  Delta-F/F from linker sequence.
- **Benchmark to beat:** 50 of ~220 CaM-glucose dehydrogenase chimeras were
  ligand-activated, a **23% hit rate** (Ergun Ayva et al., *JMB* 2022,
  10.1016/j.jmb.2022.167678). That work also found dynamic range and response
  time are anticorrelated -- a thermodynamic trade-off, not a tuning failure.
- **Where to run it:** CF2H, on hardware the deck already has, at EUR 0.10 per
  reaction. Cheaper than the design round it would replace.

Ostermeier's framing (*PEDS* 2005, 10.1093/protein/gzi048): the recombination
space of unrelated protein domains is "rich in switching behavior." Switches are
found, not designed.

## 8. What is genuinely unknown

1. **No apo Rap1 structure exists.** The whole amplitude argument rests on a
   duplex-bound structure standing in for an apo state.
2. **No split-G4-binder complementation sensor exists** in the literature.
3. **No DIP-seq or sort-seq campaign has ever run on a nucleic-acid-binding
   domain.** Every published example uses a small molecule or a peptide ligand.
   A 16-nt G4 is large and highly charged; nothing says the insertion-hotspot
   logic transfers.
4. **Rap1 is only ~2x sequence-selective** (28 nM cognate vs 45 nM for
   `GGGTGGGTGGGTGGGT`) and binds RNA G4 equally. Discrimination has to live in
   the readout, not the affinity.
5. **vdMs do not cover nucleobases**, so only the Tyr-to-chromophore contact can
   be statistically anchored.
6. **The DSM has no deposited structure as an isolated peptide** in this repo's
   records; round 1 grafted it from a predicted model.

---

No AlphaFold3, LigandMPNN or wet-lab step was run to produce this document.

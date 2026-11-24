# What RFD4-Proteina actually accepts, and what we can give it

Read 23 November 2026 from the repository's own documentation
(`RosettaCommons/RFD4-Proteina`, main branch): `docs/concepts/conditioning.md`,
`docs/concepts/annotated-cif.md`, `docs/reference/annotations.md`,
`docs/how-to/input_prep/`, and the two runnable condition-spec examples.

**This corrects two claims in `17_ROUND2_DESIGN.md` §4.1 and
`21_HINGE_LOOP_DESIGN.md` §4.** Both were too pessimistic. See §5.

---

## 1. The three authoring routes

The model consumes an **annotated CIF**: an ordinary mmCIF plus per-atom and
per-residue condition columns. You can author it three ways, all equivalent:

| route | when |
|---|---|
| **condition-spec YAML/JSON** | declarative, version-controllable, what we should use |
| Python API (`atomworks.ml.conditions`, `rfproteina.conditions`) | programmatic sweeps |
| hand-written annotated CIF | debugging |

A spec has three parts: an `input` structure, a `contig` string, and a
`conditions` block.

## 2. The contig string

Declares what chains exist in the design and how long each designed piece is.

```
contig: "A1-150, /0, B386, /0, EX1:15-30"
```

- `A1-150` — keep chain A residues 1-150 from the input
- `/0` — chain break
- `B386` — keep that residue (here a ligand)
- `EX1:15-30` — an **expandable segment**: design a new chain of 15-30
  residues, length resampled per design
- a bare range like `50-120` — a de novo chain of sampled length

Expandable segments also go *inside* a chain, which is how a variable-length
loop or linker is declared. Only polypeptide expandable segments work today.

## 3. Selection syntax

```
sel('<chain>/<res_name>/<res_id>/<atom>')     # trailing fields default to *
sel('B')                       # whole chain
sel('A/*/[1-30, 120-150]')     # two residue ranges at once
sel('A/*/[10, 20, 130]')       # a residue list
sel('A/*/56/[CG,OH]')          # named atoms of one residue
sel('A/*/*/[N,CA,C,O]')        # backbone only, whole chain
is_ligand                      # general query over any annotation
and not has_nan_coord()        # composable
```

Atom-level selection is the important part: conditions attach **per atom**, not
only per residue.

## 4. Every condition the model acts on

From `docs/reference/annotations.md`. All are listed Supported except symmetry.

| condition | what it does | values |
|---|---|---|
| `C_CRD` | fixed absolute xyz | per-atom xyz, NaN = free |
| `C_INTCRD` | internal geometry, centred and randomly rotated — shape without a global dock | per-atom xyz in own frame |
| `C_CRDGRP` | partitions atoms into rigid groups: relative geometry known **within** a group, unknown **across** groups | `-1` unknown, `0` absolute, `≥1` one group per rigid body |
| `C_SEQ` | reveal residue identity | res_name, `MASKED` = design |
| seq-probability | constrain per-identity probabilities over 28 canonical polymer tokens; explicit `0` excludes | float, NaN = unset |
| `C_IDX` | `False` = **unindexed**, model places the residue | bool, default all on |
| `C_IDXGRP` | group residues whose relative offsets are known | int |
| `C_HOT` | interface atoms the design should engage | bool |
| `C_RAS` | burial target | UNSPECIFIED / BURIED / INTERMEDIATE / EXPOSED |
| `C_HB` | request an H-bond role | NONE / DONOR / ACCEPTOR / BOTH |
| `C_PISTACK` | **request aromatic stacking with the generated partner** | bool, select whole rings |
| `C_PICATION` | π-cation role | NONE / AROMATIC / CATION / BOTH |
| `C_SEC` | secondary structure | NONE / HELIX / SHEET / LOOP |
| `C_CATH` | request a CATH fold class / architecture / topology | `1.x.x`, `2.60.40` |
| `C_CENTROID` | target centre of mass for a chain, without fixing any atom | per-atom xyz, broadcast per chain |

### The two rules that bite

1. **Conditioning any sidechain coordinate forces sequence conditioning.** The
   atom37 slot-occupancy pattern fingerprints the residue. Backbone-only
   (`N, CA, C, O`) is the exception and leaks nothing — that is how you fix a
   fold and redesign its sequence.
2. **Atomized tokens are always sequence conditioned and are never generated.**
   Ligands, non-standard residues and PTMs are context.

## 5. Two corrections to our earlier plan

### 5.1 The hinge does NOT need two separate runs

`21_HINGE_LOOP_DESIGN.md` §4 proposed generating once at the 6LDM orientation
and once at 1IGN, keeping scaffolds compatible with both. **That is
unnecessary.** `C_CRDGRP` expresses it natively in one run:

> "it can't express 'these atoms are rigid relative to each other, but I don't
> know or care where that rigid body sits in space' — which is exactly what an
> unindexed motif, or a design with several independent motifs, needs ...
> This is how one design can carry several motifs that are each internally
> rigid but placed independently of each other."

So: Rap1 block 1 in group `1`, block 2 in group `2`, both backbone-only. Each
keeps its internal fold; their **relative placement is free**, which is the
hinge. Add `C_IDX: false` with a shared `C_IDXGRP` so the chain spacing between
them is also free, and make the hinge itself an expandable segment so its length
resamples. The model was trained for exactly this case.

The two-run approach is still worth doing afterwards as a *validation*: generate
freely, then check whether the sampled relative orientations bracket the 6LDM
and 1IGN poses. That is a measurement, not an input.

### 5.2 The G4 can be present during generation, as context

The earlier claim "no nucleic acid at generation" was too strong. What the docs
actually say is that the model **generates** only the 20 standard amino acids —
"nucleic acid and small molecule generation are coming" — and that atomized
tokens are "context the model gets told about and never generates."

Giving it a G4 as fixed, sequence-conditioned context is a different thing from
asking it to design one, and the former is what we want. The sequence-probability
condition covers **28 canonical polymer identities** (20 amino acids plus, by
inference, 4 DNA and 4 RNA), so nucleic acid residues are polymer tokens in the
model's vocabulary rather than atomized blobs.

The paper's abstract says the latent space "encompasses standard and modified
amino acids, sugars, DNA, and RNA."

**What remains genuinely unvalidated is the model's *skill* at placing protein
against a four-stranded Hoogsteen-paired fold.** On the DNA-binding benchmark
RFD4 scores **3.5%** against 0.8% for BoltzGen and 0.5% for fixed-DNA RFD3, and
it succeeds on **all three** DNA targets where fixed-DNA RFD3 manages two. So it
is the best of the three and still measured on **3 duplex targets with no
wet-lab DNA result anywhere**, and with no G-quadruplex or RNA target mentioned.
That is a competence question, not an expressiveness one, and the two were
conflated before. (Correcting an earlier note that put the 3-target benchmark on
a previous release: the 3.5% is RFD4's own; the previous release is the RFD3
comparator.)

## 6. The conditions that matter for a G4 sensor

Three of these are near-exact matches for what this project needs, and none of
them were in our plan before today.

- **`C_PISTACK`** — "requests aromatic stacking contacts with the generated
  partner; select every atom of the intended aromatic ring." **G4 recognition
  by end-stacking on a terminal G-tetrad is literally aromatic stacking.** Every
  binder in `16_G4_MODULE_SURVEY.md` works this way. Select the four guanine
  ring systems of the 5' tetrad and the model is asked to build a stacking
  partner against them.
- **`C_HB`** — the Y225-equivalent contact to the chromophore phenolate is an
  H-bond with a defined role. Set ACCEPTOR on the phenolate oxygen.
- **`C_RAS: BURIED`** — the transmitting tyrosine must sit *in* the barrel hole,
  which is a burial statement, not a distance restraint.

Plus `C_SEC: HELIX` on the DSM segment (option b of §3 in
`21_HINGE_LOOP_DESIGN.md`), and seq-probability `0` on cysteine to keep the
sensor disulfide-free.

### The limitation to design around

`C_HOT` is trained as "a **target-chain** atom within 4.5 Å of any
**binder-chain** atom", count sampled `geometric(p=0.2)` and capped at 10% of
eligible atoms. `C_HB` likewise "uses **interchain** H-bond roles". So:

| contact | chains | on-distribution? |
|---|---|---|
| protein -> G4 tetrad | separate chains | **yes** |
| loop Tyr -> chromophore | same chain (CRO atomized) | **uncertain** |

The G4 half is exactly what these conditions were trained for. The
tyrosine-to-chromophore half is an *intra*chain contact, and whether atomization
of CRO makes it read as interchain is undocumented. **Test this before relying
on it** — generate a small batch with and without `C_HB` on the phenolate and
see whether the condition moves anything.

Also: hotspots must stay **sparse**. The training distribution is a small patch,
not every exposed atom; marking the whole tetrad face would move us off it.

### The chromophore gotcha

> "Do not leave any bond type as `AROMATIC`, because those all get cast to
> `SINGLE`." Kekulize by hand to `AROMATIC_SINGLE` / `AROMATIC_DOUBLE`.

CRO is an aromatic non-canonical residue, so this applies to us directly. A CIF
that leaves its ring bonds as `AROMATIC` silently hands the model a chromophore
with the wrong bond orders. Formal charges matter the same way: "leaving a
charge as 'unknown' is not equivalent to neutral."

## 7. The draft spec

Residue numbers marked `<>` are placeholders to be read off the assembled input
CIF; everything else is literal.

```yaml
g4fp_hinge_loop:
  input: inputs/rap1_barrel_g4_assembled.cif
  # A = Rap1 blocks + hinge; B = FP barrel + CRO; C = 16-nt parallel G4
  contig: "A361-440, EXHINGE:12-24, A501-573, /0, B1-<end>, /0, C1-16"

  conditions:
    # --- the hinge: two internally rigid blocks, placed independently -------
    C_CRDGRP:
      1:
        - select: "sel('A/*/[361-440]/[N,CA,C,O]') and not has_nan_coord()"
      2:
        - select: "sel('A/*/[501-573]/[N,CA,C,O]') and not has_nan_coord()"
    C_IDX:
      false:
        - select: "sel('A/*/[361-440]')"
          group: rap1
        - select: "sel('A/*/[501-573]')"
          group: rap1
    C_SEQ:
      true:
        - select: "sel('A/*/[361-440, 501-573]')"   # Rap1 sequence is known
        - select: "sel('B')"                        # barrel sequence is known
        - select: "sel('C')"                        # G4 sequence is known

    # --- the barrel is a fixed object --------------------------------------
    C_CRD:
      true:
        - select: "sel('B')"

    # --- G4 recognition: end-stacking on the 5' tetrad ----------------------
    C_PISTACK:
      true:
        - select: "sel('C/DG/[<5p tetrad resids>]/[N9,C8,N7,C5,C6,N1,C2,N3,C4,O6,N2]')"
    C_HOT:
      true:
        - select: "sel('C/DG/<one or two tetrad atoms>')"   # keep SPARSE

    # --- transmission: the loop tyrosine into the barrel hole ---------------
    C_HB:
      ACCEPTOR:
        - select: "sel('B/CRO/<cro resid>/[OH]')"
    C_RAS:
      BURIED:
        - select: "sel('A/*/<loop Tyr resid>')"

    # --- option (b): the DSM segment is a binding-induced helix -------------
    C_SEC:
      HELIX:
        - select: "sel('*/EXHINGE')"      # or the DSM segment, per variant

    # --- housekeeping -------------------------------------------------------
    # seq-probability CYS = 0 across designed positions (disulfide-free sensor)
```

Run option (c), the no-RHAU variant, by dropping the `C_SEC: HELIX` block and
letting the loop be designed freely. Everything else is shared.

## 8. What to settle before the first run

1. **Does `C_HB` act on an intrachain / atomized-residue contact?** Small batch,
   with and without. This decides whether transmission can be conditioned at all
   or has to be screened for.
2. **Is CRO kekulized** in our input CIF, and are its formal charges set?
3. **Does free `C_CRDGRP` sampling bracket the 6LDM and 1IGN orientations?**
   If it does, the hinge is being modelled; if every sample lands in one pose,
   it is not.
4. **How sparse should the tetrad hotspots be?** Start at one or two atoms,
   following `geometric(p=0.2)`.

Items 1 and 3 are the real experiments. Both are cheap and both are
falsifiable.

---

Documentation read 23 Nov 2026. No AlphaFold3, LigandMPNN, RFD4-Proteina or
wet-lab step was run to produce this document.

---

## 9. The prose route is real, but it is an interview

From `docs/how-to/input_prep/agent.md` and the project site.

> "Tell your coding agent what you want to build. It interviews you about the
> parts you left open, checks the structure for the usual traps, and hands back
> a model ready annotated CIF."

> **"The conversion is not the product. The interview is."**

So "just tell it what you want" is close but not the whole claim. The mechanism:

- a shared skill **`rfd4-proteina-annotated-cif`**, run inside Codex or Claude
  Code — `Use rfd4-proteina-annotated-cif with <artifact>`;
- Claude Code's **`/design2cif`** command wraps the same skill;
- `rfd4-proteina-generate` hands input prep to it.

It opens the structure rather than reasoning about it, restates each intent as
*"I read this as X, is that right?"*, then edits the annotated CIF **directly
through the AtomWorks Python API**. It will **not** create a YAML/JSON
intermediate unless explicitly asked — so if we want the spec in §7 as a
version-controlled artifact, we have to ask for it.

Evidence it accepts: "A paragraph, lab note, source structure, existing
annotated CIF, or condition-spec JSON/YAML."

What it refuses to do: guess past a real ambiguity, silently choose structural
preparation, silently drop covalent bonds / charges / PTMs, or claim a check it
did not run. Its own warning: "An agent can misread intent. Before you run the
model, inspect the annotated CIF and its validation report."

### Additional gotchas it screens for

Beyond the kekulization issue already recorded:

| gotcha | why it matters to us |
|---|---|
| `"100-155"` is a distribution, not a number; `"180-180"` is one length | our hinge and loop lengths must be deliberate, not accidental ranges |
| length counts only what is **designed** — 180 plus a 4-residue motif is 184 | affects every contig we write |
| **leaving ligand coordinate conditioning unset conditions no pose** — the model invents the ligand's pose *and* conformer, getting only identity, bonds and charge | **applies to the G4 and to CRO.** If we do not condition G4 coordinates we are asking it to invent a quadruplex fold |
| `?` in `pdbx_formal_charge` means unknown, not neutral | CRO and the G4 phosphates |
| PDB-derived ligands often arrive single-bonded | our CRO almost certainly will |
| **covalent polymer-to-ligand bonds must remain present** | CRO is covalently in the chain; losing that bond silently detaches the chromophore |

## 10. The calmodulin result — our exact problem, and the part nobody documents

The site's headline example is a **conformational-state-selective binder**, on
the very protein family our sensor mechanism is borrowed from:

> "Calmodulin switches conformation when it binds Ca2+, and many of its partners
> recognize only one state. We asked RFD4-Proteina for binders to the
> calcium-bound form only."

- **12 of 18 designs** bound with apparent Kd below 1 uM **and showed no
  detectable binding without calcium** -- a 67% hit rate on state selectivity.
- Apparent Kd from 69 pM (with a published caveat that this is a 1:1 fit at the
  off-rate fitting floor, so its precise value and ranking are unresolved) to
  37.0 nM; the site's interactive example binds at 2.48 nM.

**The site gives no mechanism for how the state selectivity was specified.** No
condition, parameter or input format is named for it, and there is nothing in
the annotations table (§4) that expresses "bind state A but not state B". The
nearest adjacent statements are "multi-state AF3 filtering" used on the ligase
example -- a post-generation filter, not a specification -- and a general
mention of "extrinsic reasoning: search steered by any score".

**The reading that fits the evidence: you specify the state you want by giving
that conformer's structure as the target, and selectivity against the other
state is emergent, then filtered.** You do not negative-design against state B
in the generator.

If that is right, it transfers cleanly and it simplifies our plan:

- give the **G4-bound** Rap1 conformer as the conditioned target;
- do **not** attempt to encode "and not the duplex/apo state" as a condition;
- recover the selectivity afterwards, where we already have the machinery --
  Delta-iPAE for occupancy contrast and Delta-pLDDT for the loop.

**Worth confirming before relying on it**, since it rests on absence of evidence
rather than a stated method. The paper (10.64898/2026.11.12.766037) may say;
the site does not.

Encouraging either way: a 67% hit rate on exactly the class of problem we have,
and no other campaign on that page needed a special mechanism for it.

## 11. The target brief

The prose to hand to `/design2cif`. It is written to be interviewed about, so
the open decisions are marked rather than guessed.

> **Target.** Build a genetically encoded G-quadruplex sensor. Three components
> are supplied in one input CIF: (A) the yeast Rap1 DNA-binding domain split
> into two blocks, residues 361-440 and 501-573; (B) a fluorescent-protein
> beta-barrel carrying its mature chromophore CRO, covalently in the chain; and
> (C) a 16-nt parallel G-quadruplex, sequence GGGTGGGTGGGTGGGT.
>
> **What is known and must not be designed.** The barrel is a fixed object --
> condition its backbone coordinates and its full sequence. Both Rap1 blocks
> have known sequence and known internal fold, so condition them backbone-only
> so their sequences stay fixed but their side chains are free. The G4's
> coordinates and sequence are known and must be conditioned, because leaving
> them unset would ask the model to invent a quadruplex fold.
>
> **What must be designed.** First, the hinge between the two Rap1 blocks: a
> variable-length segment, roughly 12-24 residues, resampled per design.
> Second, an insertion loop that runs from Rap1 block 2 into the barrel wall,
> entering between two beta-strands. Third, one aromatic residue near the end of
> that loop which must sit buried in the barrel's surface cavity and hydrogen
> bond to the chromophore's phenolate oxygen.
>
> **The motion the design has to allow.** The two Rap1 blocks are each
> internally rigid but their relative placement must be free -- put them in
> separate coordinate groups and unindex them as one index group. Their relative
> orientation is the switch: in the G4-bound structure 6LDM they sit roughly 104
> degrees apart from where they sit in the duplex-bound structure 1IGN. Do not
> pin either block in absolute space.
>
> **How the protein should recognise the G4.** By aromatic end-stacking on a
> terminal G-tetrad, which is how every characterised G4-binding protein works.
> Request pi-stacking against the complete guanine ring systems of the 5' tetrad.
> Add a small, sparse set of hotspots on that face -- one or two atoms, not the
> whole surface. Keep the 3' tetrad free: Rap1's own recognition helix
> (residues 398-412, hydrophobic patch Ser402/His405/Val409) docks there.
>
> **Variant to run alongside.** Replace the designed insertion loop with the
> RHAU/DHX36 G4-recognition motif, conditioned as a helix, positioned so its
> stacking face points at the 5' tetrad. In that variant the DSM supplies the
> disorder-to-order transition instead of a de novo loop. Everything else is
> identical.
>
> **Housekeeping.** Exclude cysteine from all designed positions. CRO's
> aromatic bonds must be kekulized and its formal charges set explicitly, not
> left unknown. Preserve the covalent bond joining CRO to the chain.
>
> **Open decisions I expect to be asked about.** Whether the hinge length range
> should resample per design or be fixed; whether the barrel insertion point is
> fixed or should be sampled across candidate beta-strand positions; whether the
> loop tyrosine's identity should be conditioned or designed; whether the
> hydrogen-bond condition can act on an intrachain, atomized-residue contact at
> all (see section 6); and how sparse the tetrad hotspots should be.

---

Site and documentation read 23 Nov 2026. No RFD4-Proteina, AlphaFold3,
LigandMPNN or wet-lab step was run to produce this document.

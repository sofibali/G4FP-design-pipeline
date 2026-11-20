# G4-binding modules, surveyed — SUD-M, and what else exists

Literature retrieved from PubMed; structures from the RCSB PDB. **All numbers in
the tables marked *measured here* were computed in this directory from deposited
coordinates**, not taken from the papers — where a paper's number and mine
disagree, both are given.

Extends `13_G4_BINDING_MODULES.md`, which covers Nb55, SG4, G4P, BG4 and the
RHAU DSM. Written 19 November 2026.

---

## 0. The question this answers

Round 1's module is residues 1–28 + 279–306 — the two chain ends, packed into
one compact unit against the barrel. Anything that replaces it has to (a) fit
that volume, (b) bind a G4 selectively, and (c) **change shape when it does**,
because `13_G4_BINDING_MODULES.md` established that Nb55 does not.

### The slot, measured

| | residues | Rg | max span | termini |
|---|---:|---:|---:|---|
| **G4FP module slot** | 56 | **13.5 Å** | 32.3 Å | attaches to barrel at two points **5.2 Å** apart |
| Nb55 folded core (1–105) | 105 | 13.3 Å | 39.3 Å | N-to-C 12.0 Å |
| SUD-M alone (2JZD) | 129 | 13.7 Å | 45.9 Å | N-to-C 26.8 Å |
| SUD-NM, the selective one (8GQC) | 264 | **19.3 Å** | **60.0 Å** | N-to-C 26.2 Å |
| *the barrel itself, for scale* | 250 | 18.2 Å | 60.6 Å | — |

A ~105-residue folded domain is a near-exact volumetric match for the slot
(Rg 13.3 vs 13.5 Å). **So the architecture can host one without redesigning the
barrel.** Beyond ~130 residues it stops being a module and starts being a second
protein.

*Caution recorded:* measured over the full deposited chain, Nb55's N-to-C
distance is 46.1 Å, which would rule the graft out. That is an artefact —
residues 106→119 move outward at a constant 2.9 Å per residue, a fully extended
disordered His-tag. The folded domain ends near residue 105.

---

## 1. SUD-M — rejected, with reasons

The SARS-unique-domain middle macrodomain of coronavirus nsp3.

**It fails on selectivity, and the failure is published.** Lavigne *et al.*,
*Nucleic Acids Research* 2021 ([10.1093/nar/gkab571](https://doi.org/10.1093/nar/gkab571),
PMID 34232992): the isolated SUD-M "interacted with all the different
oligonucleotides which suggests that this domain is too short or too positively
charged (pI > 9), to discriminate between G4 and non-G4 folded structures."
Selectivity appears only in the two-domain SUD-NM — which the table above shows
is the size of the whole barrel, and whose interdomain linker is **unmodelled
for 9 residues even in a 1.35 Å structure** (8GQC, verified here).

Three further disqualifiers:

1. **No complex structure exists, for either virus, by any method.** Verified
   against RCSB: 2W2G, 2WCT, 2JZD, 2RNK, 2KQV, 2KQW, 2KAF, 8GQC, 8HBL, 7THH —
   every one `homomeric protein`, zero nucleic-acid entities. Every published
   "binding site" is docking on an apo structure.
2. **The G4 interpretation is disputed.** Johnson, Chatterjee, Neuman &
   Wüthrich, *J Mol Biol* 2010 ([10.1016/j.jmb.2010.05.027](https://doi.org/10.1016/j.jmb.2010.05.027),
   PMID 20493876) found K⁺ *reduces* binding — the opposite of G4-stabilising —
   and that SUD-M binds poly(A)₁₀ and random purine RNA with no dependence on
   G-runs. Tan *et al.* 2009 ([10.1371/journal.ppat.1000428](https://doi.org/10.1371/journal.ppat.1000428))
   reported the opposite K⁺ dependence. The two labs do not reproduce each other.
3. **No transferable dead-binder control.** SARS-CoV mut4 (K565A/K568A/E571A)
   abolishes binding, but its SARS-CoV-2 homologue has "very minor effect"
   (Lavigne), and Qin *et al.* 2023 ([10.1038/s41467-023-39709-6](https://doi.org/10.1038/s41467-023-39709-6))
   found E595A negligible. Compare SG4, where Addgene #196072 is a ready-made
   matched R105A control.

Plus a sensor-specific liability: **SUD oligomerises on G4** (DLS, Tan 2009;
NMR line-broadening to precipitation, Johnson 2010). No 1:1 stoichiometry.

Useful residue: Addgene **#156470** (`pETM33_Nsp3c_SUD-M`, SARS-CoV-2), with
#156469 and #156471 as the SUD-N and SUD-C siblings, if the domain is ever
wanted in hand.

---

## 2. How many graftable options actually exist

Protein PDB entries with a genuinely resolved, protein-contacted G4: **~20,
across ~14 families.** Made by a folded domain rather than an RGG/peptide:
**8 families.** Natural: **6.** Natural, folded and under ~300 aa: **2** —
yeast Rap1 and nucleolin RBD12. **Of those two, neither has an apo structure.**

Every apo/holo pair in the entire G4 field is on an *engineered* binder.

---

## 3. The rigidity result is now a pattern, not an anecdote

| Binder | apo / holo | Cα RMSD | *measured here* |
|---|---|---|---|
| Nb55 nanobody | 9GV4 / 9GXH | **0.83 Å** (noise floor 0.74 Å) | yes |
| DARPin 2E4 | 8IPF / 8IPP | **0.75 – 1.02 Å** across all 8 apo copies | yes |

Two unrelated scaffold families, two independent selection campaigns, the same
answer. The shared mechanism explains it: both bind by **flat-surface
end-stacking on a terminal G-tetrad** plus backbone electrostatics — the same
mode as the RHAU DSM and as small-molecule G4 ligands. A flat tetrad face
selects for a pre-organised flat complementary surface.

**Treat this as a design constraint, not a search failure. Selecting harder for
G4 binders will keep producing rigid ones.** The switch has to come from
somewhere other than the binder's own backbone.

DARPin 2E4: Ngo, Liew, Heddi & Phan, *JACS* 2024
([10.1021/jacs.4c01971](https://doi.org/10.1021/jacs.4c01971), PMID 38738955);
selection in Scholz, Hansen & Plückthun, *NAR* 2014
([10.1093/nar/gku571](https://doi.org/10.1093/nar/gku571), PMID 25053846).
Note 8IPP resolves 2 clean tetrad layers, not 3, at 2.0 Å — verified here
(12 guanines, 11 Hoogsteen pairs).

---

## 4. Rap1 — the one large conformational change, and what it really is

Traczyk, Liew, Gill & Rhodes, *Nucleic Acids Research* 2020
([10.1093/nar/gkaa171](https://doi.org/10.1093/nar/gkaa171), PMID 32187364).
*S. cerevisiae* Rap1 DNA-binding domain: **6LDM** G4-bound (2.4 Å) vs **1IGN**
duplex-bound (2.25 Å).

**Both structures verified here.** 6LDM: 12 guanines, 14 Hoogsteen G-pairs →
3 tetrad layers, a real G4. 1IGN: 24 guanines, **0** Hoogsteen pairs — duplex,
as expected.

### Where the motion is — measured here, and it is not where the paper's summary puts it

Sliding 40-residue windows, G4-bound vs duplex-bound, internal Cα RMSD:

| window | internal RMSD | |
|---|---:|---|
| 361–400 | 0.69 Å | rigid |
| 381–420 | 0.67 Å | rigid |
| 401–440 | 0.90 Å | rigid |
| **421–460** | **5.99 Å** | **the hinge** |
| 441–480 | 3.60 Å | still spanning it |
| 501–540 | 1.40 Å | rigid |
| 521–560 | 0.84 Å | rigid |

So Rap1 is **two rigid blocks — roughly 361–440 and 501–573 — joined by a hinge
at ~440–460.** Global Cα RMSD over the 175 shared residues is **11.33 Å**
(the paper's figure, reproduced). Superposing on the N-terminal block alone, the
C-terminal block's centroid moves **14.7 Å**, individual residues up to
**35.7 Å**, with a relative rotation of **104°** — consistent with the paper's
"twisted 117°, displaced 36 Å" within the choice of domain window.

**Correction worth recording:** the claim that *both* Myb domains are internally
rigid at 2.66 / 1.26 Å does not reproduce. Myb2 is rigid (1.31 Å measured), but
a 386–460 "Myb1" window gives **6.58 Å**, because that window contains the
hinge. The rigid part of Myb1 is 361–440.

### The design-relevant fact

**The Myb1 recognition helix (398–412) has an internal RMSD of 0.46 Å.** The
helix that reads the G-tetrad does not deform — it gets *repositioned*. Its
mechanism is a single bimodal α-helix: Ser402/His405/Val409 form a hydrophobic
patch that docks flat on the 3′ G-tetrad, while Asn401/Arg404/Arg408/Ser410 on a
face rotated ~120° make the base-specific duplex contacts. One helix, two faces,
two DNA structures.

This is the same lesson as Nb55's 41–44 loop, at ten times the amplitude:
**couple the dye pocket across the 440–460 hinge and you get a ~35 Å throw;
couple it inside either rigid block and you get nothing.**

### Three caveats that decide whether it is usable

1. **This is duplex-bound vs G4-bound, not apo vs holo.** There is no DNA-free
   structure of the yeast Rap1 DBD. So it is a *structure-discrimination*
   switch, not an *occupancy* switch — a different quantity from the one
   `13_delta_ipae_gate.py` measures.
2. **The apo state is probably floppy.** The authors infer the interdomain
   linker is unstructured in solution, so apo→G4 may be disorder-to-order with a
   heterogeneous starting state rather than a clean two-state hinge.
3. **It is not sequence-selective.** EMSA apparent Kd: cognate yeast telomeric
   `4-G3` 28 nM; non-telomeric `GGGTGGGTGGGTGGGT` core 45 nM — only ~2× worse.
   It binds RNA G4 as well as DNA G4. It reads the tetrad face, not the
   sequence, and its day job is high-affinity telomeric *duplex* binding.

**But caveat 3 may not matter for this sensor.** A G4 sensor's real job is
telling a G4 from the duplex and ssDNA around it, and Rap1 is the only natural
folded domain documented to adopt two grossly different conformations on those
two substrates. The discrimination lives in the readout rather than in the
affinity.

### Grafting evidence, which is unusually strong

Peptides of the isolated recognition helices bind G4 at 1.2–5.4 µM — **40–200×
weaker than the folded DBD.** Unlike the RHAU DSM, this helix does *not* work as
a free peptide, so the fold is load-bearing. The Myb2 helix peptide is a further
4–5× weaker, traceable to His405→Asp543, Val409→Lys547, Ser402→Ala540 — a
natural built-in negative-control pair.

Sizes: DBD construct 360–598 = 239 aa (~27 kDa); **Myb1 alone ≈ 90 aa, ~10 kDa**,
which fits the slot. The isolated DBD binds G4 at affinity comparable to
full-length Rap1; the BRCT and RCT domains show no binding at all.

**Gap: no loss-of-binding point mutants exist.** H405A / V409A / S402A on the
folded DBD are obvious and unclaimed.

---

## 5. Nucleolin RBD12 — the selectivity option

Chen *et al.*, *Science* 2025 ([10.1126/science.adr1752](https://doi.org/10.1126/science.adr1752),
PMID 40245140). **PDB 9CB5**, 2.6 Å, RRM1–RRM4 + Myc G4, intact 3-tetrad
parallel G4 with 2 channel K⁺ per copy.

The best *selectivity* in this survey, and it directly settles the
"nucleolin RBDs versus its RGG" question: **the folded RRMs do the work.** RBD1
and RBD2 make all the G4 contacts; RBD3 and RBD4 make none in the crystal.
Recognition is conformation-based, through the G4's 6-nt central loop and 5′
flank — it reads loop and flank identity, not just a tetrad face. It binds Myc
G4 more tightly than its own cognate NRE RNA, and CUT&Tag confirms it in cells.

RBD12 ≈ 180 aa, ~20 kDa — larger than the slot wants. Two problems: high
affinity **requires all four RBDs** (multivalency, not a single rigid module),
and the crystal needed a Fab chaperone, hinting the complex is not rigid alone.
Apo RRM structures exist (2FC8, 2FC9, 2KRR) but with different domain
boundaries, so no clean apo/holo comparison is available.

---

## 6. Also real, but not graftable

- **ToPif1 "G4-recognising surface"** (7OAR, 2.58 Å) — intact 3-layer parallel
  G4, 2 channel K⁺, and **the best mutant set in the field**: R355A 133 nM,
  R135A 119 nM, R150A 152 nM vs WT 29.6 nM, all CD/DLS-verified as folded, all
  with k_cat/K_m → 0 for unwinding. Domain 2B shifts 5.40 Å / rotates 25.7°, and
  loop3 becomes an ordered β-hairpin — genuine local disorder-to-order at the
  interface. But 444 aa, and the epitope is discontinuous across two subdomains.
  Dai *et al.*, *EMBO Rep* 2022 ([10.15252/embr.202153874](https://doi.org/10.15252/embr.202153874)).
  **The loop3→β-hairpin ordering is the design-relevant half, and it is local.**
- **DHX36 helicase core** — yes, the folded core contacts the G4, answering that
  question directly; the DSM folds into a binding-induced α-helix that works
  with the OB-fold subdomain to select parallel G4s. Cleanest structures are
  8VX1/8VX8 (DNA G4) and 8VV2/8VVD (RNA G4), Banco *et al.*, *Cell Rep* 2025
  ([10.1016/j.celrep.2025.116136](https://doi.org/10.1016/j.celrep.2025.116136)).
  Apo 5VHA vs G4-bound 5VHE differ by ~6.9 Å over ~800 residues, but the apo
  constructs **lack the DSM**, so an unknown share of that is construct
  difference. ~800 residues; not a module.
- **RECQ1** (9I22, 2.42 Å) — documented closed→open, but only 8–15 interface
  residues per G4. An encounter complex, not a recognition module.
- **CMG replisome** (9E2Z) and **PRC2** (8FYH) — real G4s, far too large.

**Viral G4-binding domains: the category is empty.** No experimental structure
of any viral protein bound to a G4 exists, for any virus. Not worth more time.

**No de novo / RFdiffusion / BindCraft-designed G4 binder exists.** The nearest
historical attempt is phage-display three-finger zinc fingers (Isalan, Patel,
Balasubramanian & Choo, *Biochemistry* 2001,
[10.1021/bi001728v](https://doi.org/10.1021/bi001728v)) — ~90 aa, folded, binding
reported structure-specific, but no structure, no Kd, no mutants, never followed
up in 25 years. Given that an ankyrin scaffold worked, this is an open
opportunity.

---

## 7. Traps — entries whose titles promise a G4 complex and do not deliver

Recorded because two of them were on our own candidate list.

1. **1JB7 / 1PH1, Oxytricha TEBP α/β.** A genuine 4-tetrad antiparallel G4 *is*
   present — but it is not the protein's ligand. Contact analysis: the ssDNA
   chain makes 42 residue contacts in the α/β cleft; one G4 chain makes 4
   peripheral contacts and the other makes **zero**. Horvath & Schultz say so
   themselves: "the G-quartet:protein contact surfaces are modest and might
   simply represent crystal packing interactions"
   ([10.1006/jmbi.2001.4766](https://doi.org/10.1006/jmbi.2001.4766)).
   **The most misleading entry in the field.**
2. **6CRM**, *C. sakazakii* RecQ "bound to an unfolded G-quadruplex" — the entire
   G4-forming region (nt 1–18) is absent from the coordinates. Zero tetrads.
3. **1VYA** ("first plant G4-binding protein, ZmNDPK1") and **6TAZ** (Timeless /
   DDX11) — both contain **zero nucleic acid entities**. Apo protein.
4. **Thrombin/TBA series** (4DIH, 4LZ1, 5CMX, 6EVV, 6GN7, 7NTU, 8BW5), and
   **HupB / DevR** — direction reversed: these are G4 *aptamers* binding a
   protein, not G4-binding domains. Thrombin has no G4-reading surface.
5. **8QN2** — the "protein" entity is an aromatic oligoamide foldamer.
6. **POT1, TRF1, TRF2, CST, Cdc13, BLM, WRN, RPA, PARP1, CNBP, hnRNP A1, Rif1,
   TERT** — all checked, all nucleic-acid entries are ssDNA or duplex. **No
   WRN–G4 or human Pif1–G4 structure exists anywhere.** CNBP has no PDB entry
   with nucleic acid at all.

---

## 8. What this decides

1. **SUD-M is out**, and so is the whole "find a natural viral G4 domain" line.
2. **Rigidity is a property of G4 binders, not of our choices.** Two verified
   apo/holo pairs, two unrelated scaffolds, 0.75–1.02 Å. Stop looking for an
   induced-fit G4 binder; design the switch somewhere else.
3. **Rap1 is the only natural lead with real amplitude**, and the amplitude is
   inter-domain, not internal — a 104° rotation and up to 35.7 Å across a hinge
   at residues ~440–460, with both flanking blocks and the recognition helix
   itself rigid to under 1 Å. That is a better sensor geometry than induced fit,
   because two rigid bodies changing relative orientation is what a split
   reporter or a circularly-permuted FP actually reads.
4. **But it is duplex-vs-G4, not apo-vs-G4**, and there is no apo structure. The
   cheapest way to find out whether it is a usable two-state switch is an AF3
   apo ensemble of the Rap1 DBD — hours, not weeks, and it reuses the ΔiPAE
   machinery already built.
5. **If selectivity matters more than switching**, nucleolin RBD12 is the only
   folded domain with demonstrated G4-conformation selectivity — at the cost of
   multivalency and size.
6. **PDB 9KCQ** (1.77 Å ankyrin + human telomeric G4, 3 tetrads) is public and
   unpublished, higher resolution than 8IPP and on a more relevant sequence, with
   interface residues at nearly the same ankyrin positions as 2E4. Worth mining.

# Round 2 — campaign plan, costed

**Status: planning. Nothing in this document has been run.** Every number is
either a measurement taken from the round-1 tree on this server (marked
*measured*) or a projection from one (marked *predicted*, with the assumption
that produces it stated). No new prediction, design or wet-lab step was
executed to produce any of it.

Produced 17 November 2026 by `14_regroup.py` and `15_round2_forecast.py`.
Reproduce with:

```bash
python3 -I 14_regroup.py && python3 -I 15_round2_forecast.py
```

---

## 1. Regrouping — where round 1 actually died

`13_RESCORING_PLAN.md` reported the cascade as a funnel: 999 → 298 → 197 → 101
→ 71. That tells you how many survive each stage. It does not tell you how the
928 non-survivors are distributed across causes, because a design that fails
gate 1 is discarded before gates 2–4 ever see it. The cause is what round 2 has
to fix, so `14_regroup.py` assigns every design exactly one failure mode — the
*first* gate it fails — and regroups.

*measured*, from `14_regroup_failuremode.csv`:

| First gate failed | n | share | What it means |
|---|---:|---:|---|
| **CONTACT** | **701** | **70.2%** | AF3 predicts < 3 module–G4 contacts |
| NEGATIVE | 101 | 10.1% | module already shut with no G4 present |
| BIND | 96 | 9.6% | interface PAE in the worst 60% |
| SWITCH | 30 | 3.0% | no apo/holo confidence contrast |
| **PASS** | **71** | **7.1%** | clears all four |

**Seven in ten designs die because the model does not think the protein touches
the DNA at all.** Not because the interface is weak, not because the switch is
small — because there is no interface. That one number reorganises the round-2
priority list: the switch gate is not the bottleneck, contact is.

And the gates are not independent. Taken alone, SWITCH passes 168/999 (16.8%)
and BIND passes 400/999 (40.0%); if they were independent, 67 designs would
pass both. 130 do. Designs that bind well also switch well, which is consistent
with the partial correlation already measured (ρ = −0.321 controlling for model
quality) and means the cascade is less selective than its stage counts suggest.

### 1.1 The library is four templates wide, not ten

*measured*:

| Template | n | CONTACT | pre-closed | weak | no switch | **PASS** | rate |
|---|---:|---:|---:|---:|---:|---:|---:|
| `des1_cro_mod4` | 100 | 30 | 0 | 21 | 2 | **47** | **47.0%** |
| `r28_l2_cro_mod2` | 100 | 21 | 0 | 52 | 18 | **9** | 9.0% |
| `r28_l3_cro_mod2` | 100 | 91 | 0 | 0 | 0 | **9** | 9.0% |
| `r28_l2_cro_mod3` | 100 | 83 | 0 | 14 | 0 | **3** | 3.0% |
| `r28_l3_cro_mod0` | 100 | 88 | 0 | 9 | 1 | **2** | 2.0% |
| `r28_l1_cro_mod2` | 99 | 98 | 0 | 0 | 0 | **1** | 1.0% |
| `des1_cro_mod0` | 100 | 99 | 0 | 0 | 1 | **0** | 0.0% |
| `r20_l1_cro_mod0` | 100 | 12 | **88** | 0 | 0 | **0** | 0.0% |
| `r20_l2_cro_mod2` | 100 | 87 | 13 | 0 | 0 | **0** | 0.0% |
| `r20_l2_cro_mod3` | 100 | 92 | 0 | 0 | 8 | **0** | 0.0% |

Four of the ten templates produce **zero** survivors. One, `des1_cro_mod4`,
produces 47 of the 71. Round 1 spent a tenth of its compute on each template;
on this accounting it should have spent it on two or three.

`r20_l1_cro_mod0` is a different animal from everything else and deserves a
name. Its median apo interface PAE is **3.8 Å** against 20–24 Å for every other
template; its median iPTM is the highest in the set (0.750) and its median
interface PAE the lowest (9.2 Å). It builds a protein that is confidently
packed against itself with or without the G4 — 88 of 100 pre-closed. It looks
best on every binding metric and cannot switch, which is exactly the failure
mode the NEGATIVE gate was added to catch. All three of its tested designs were
dark.

---

## 2. The result that should govern round 2

Regrouping the 21 experimentally tested designs by the same gates
(*measured*):

| First gate failed | tested | fluorescent |
|---|---:|---:|
| CONTACT | 5 | 0 |
| NEGATIVE | 3 | 0 |
| BIND | 2 | **1** (`r28_l2_cro_mod2_design_0072`, the lead) |
| SWITCH | 1 | **1** (`r28_l2_cro_mod2_design_0016`) |
| **PASS** | **10** | **0** |

**Ten of the tested designs clear all four gates. All ten were dark. Both
designs that ever fluoresced are rejected — one at BIND, one at SWITCH.**

On the fluorescence phenotype this cascade is not merely uninformative, it is
inverted. That has to be said plainly before any of the compute below is spent.

Two readings, and the plan has to hold both:

1. **The defensive reading, which is correct as far as it goes.** The cascade
   predicts *switching*, not *brightness*. 19 of 21 were dark and 0 of 21 were
   G4-responsive, so the phenotype it predicts has zero positives in this
   dataset and cannot be scored against it. The two fluorescent designs are
   rejected for coherent mechanistic reasons — `design_0016` is the
   self-stabilised pKZ168 construct, which *should* fail a switch gate.
2. **The reading that actually changes the plan.** Whatever the cascade selects
   for, it demonstrably does not select for chromophore maturation, and
   maturation is where 19 of 21 round-1 constructs failed. Adding gates that
   are silent about the dominant observed failure mode does not raise the hit
   rate. Round 2 must either gate on maturation separately, or delete the
   maturation step.

Deleting it is what Option E does, and this table is the strongest argument for
it that the data contains.

---

## 3. Measured cost anchors

All taken from this tree. None is a vendor figure.

| Quantity | Value | Source |
|---|---|---|
| AF3 holo, 5 seeds × 5 samples, 326 tokens | **141.6 s / prediction** | median of 10 per-template medians of inter-directory mtime intervals, `output_G4FP_*/03_alphafold3_predictions_bound/` |
| AF3 apo, 307 tokens | **138.2 s / prediction** | same, `…_apo/` |
| AF3 concurrency | **2 GPUs**, two interleaved streams ≈ 1 h offset | the mtime pattern across the 20 template × state runs |
| AF3 round 1 total | 2,002 predictions, **40.6 h** calendar | first/last prediction mtime, 10 Mar 16:19 → 12 Mar 08:53 2026 |
| LigandMPNN | **32.3 min** per 10,001 sequences per input | `batch_ligandmpnn_progress.log`, the two runs not skipped (29m43s, 34m56s) |
| Downstream analysis | **11 min** per 100 designs, both states, CPU | step 05 → 06 output mtimes, `r28_l3_cro_mod2` |
| Sequence filter | 10,001 → 100 (top 1%) | `pipeline_config.yaml` |

Consistency check: 2,002 × mean(141.6, 138.2) s ÷ 2 GPUs = **38.9 h** predicted
against **40.6 h** measured calendar. The 4% gap is scheduling overhead between
the 20 runs. The model is sound.

### 3.1 What the AF3 inputs actually contain

Read from `design_0_bound_data.json`, all ten templates:

- chain A — protein, **287–307 residues**, CRO as a bonded modified residue
- chain B — DNA, **`GGGTGGGTGGGTGGGT`**, 16 nt, identical across all ten
- chain C — **3 K⁺**

and the apo state is **chain A alone** — so "apo/holo" in round 1 means
*without / with the G4*, not without/with the chromophore. The ΔiPAE gate is
measuring G4-dependence, which is the right thing, and this confirms it.

**Flag.** `MISSING_STRUCTURES.md` and the round-2 plan both specify a **22-nt
Myc G4**. Round 1 was predicted against a 16-nt GGGT repeat. These are
different targets with different topologies. If round 2 switches to Myc, every
percentile threshold in `13_RESCORING_PLAN.md` is calibrated against a ligand
that is no longer being used and must be recomputed — cheap, but it must
happen. The choice also decides the binding module: Nb55 binds TBA and AGCGA,
SG4 was raised against Myc (`13_G4_BINDING_MODULES.md`).

**Second flag.** `MISSING_STRUCTURES.md` says 7 of the 21 tested variants have
no AF3 structure. That was true of Drive and of git. It is **not** true here:
`14_regroup.py` matched **21 of 21** against the full prediction set on this
server. The document is stale and the seven re-predictions it asks for are not
needed.

---

## 4. Predicted cost per prediction, by round-2 construct

*predicted*. Assumption **A1**: round 1 contains exactly one token-scaling step
— 307 → 326 tokens (+6.2%) cost 138.2 → 141.6 s (+2.5%), an apparent exponent
of 0.41. A single 6% step is too short a lever to extrapolate from, so both the
measured 0.41 and a pessimistic quadratic 2.0 are carried.

| Construct | tokens | holo s/prediction | vs round 1 |
|---|---:|---|---|
| round 1 as run (16-nt G4) | 326 | 141.6 | — |
| same protein, 22-nt Myc G4 | 332 | 142.7 – 146.9 | +1% to +4% |
| + rhodamine (Option E) | 365 | 148.3 – 177.5 | +5% to +25% |
| Nb55 graft + rhodamine | 341 | 144.2 – 154.9 | +2% to +9% |

**The construct barely changes size.** Even the worst case is a quarter more
per prediction. Round-2 compute is set by *how many* predictions, not by what
is in them — which means the scenario choice in §5 is the whole budget
decision.

---

## 5. Predicted campaign cost and predicted funnel

*predicted*, from `15_round2_forecast.csv`. Each scenario keeps the round-1
shape: *n* inputs × 10,000 LigandMPNN sequences → top 100 each → AF3 holo +
apo → the four gates.

| | Scenario | inputs | designs | AF3 preds | MPNN | AF3 GPU-h | **wall, 2 GPUs** | predicted survivors |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| **A** | Kielen plan as agreed — the 21 tested variants | 21 | 2,100 | 4,200 | 11.3 h | 163 | **3.8 d** | 149 – 593 |
| **B** | Iterate all 71 survivors | 71 | 7,100 | 14,200 | 38.2 h | 552 | **12.8 d** | 505 – 4,276 |
| **B-lite** | 12 survivors, stratified across templates | 12 | 1,200 | 2,400 | 6.5 h | 93 | **52 h** | 85 – 723 |
| **E** | Option E rebuild — new backbones + rhodamine | 8 | 800 | 1,600 | 4.3 h | 62 | **35 h** | not predictable |

### 5.1 The gate rates behind those numbers

*measured*, conditional, over the full 999:

| Gate | Question | pass | |
|---|---|---:|---|
| CONTACT | ≥ 3 module–G4 contacts | 29.8% | 999 → 298 |
| NEGATIVE | apo iPAE ≥ 8 Å | 66.1% | 298 → 197 |
| BIND | interface PAE ≤ P40 (19.72 Å) | 51.3% | 197 → 101 |
| SWITCH | ΔiPAE > 3 Å | 70.3% | 101 → 71 |
| | **all four** | **7.1%** | 999 → 71 |

### 5.2 Why the survivor ranges are so wide, and what to do about it

Assumption **A2 (inheritance)**: a LigandMPNN child of a design that passed
gate X may or may not pass it. Round 1 cannot measure this — it contains no
parent/child pairs. The band is bounded by:

- **pessimistic** — redesign randomises completely: the library rate, 7.1%.
- **optimistic** — children inherit their parents' template rate exactly
  (33.6% weighted across the survivor templates), times the window uplift below.

Assumption **A3 (window)**: round 1 used `A15-60 A240-285`, capturing 140 of
283 G4 contacts (49%). Round 2 uses first-40 + last-40, capturing 251 (89%) —
`README.md` §3. CONTACT is the gate that kills 70% of round 1 and is precisely
what that window change targets, so the ratio 251/140 = **1.79** is carried as
the optimistic uplift and 1.0 as the pessimistic one.

Scenario B's optimistic bound — 4,276 survivors from 7,100 designs — is not a
credible forecast. It is the arithmetic consequence of assuming perfect
inheritance, and no design round has ever shown that. Its real value is to show
that **the band is too wide to budget against**: an 8× uncertainty on the
output of a 12.8-day run is not a plan. §7 closes it for 3 hours of compute.

---

## 6. How many to order — the constraint that actually binds

*measured*: round 1 ordered **21**, got **2 fluorescent** (9.5%) and **0
G4-responsive**.

*predicted*, from the binomial:

| If the true switch rate is | n for an 80% chance of ≥1 hit | P(0 hits) at n = 21 |
|---:|---:|---:|
| 3% | 53 | 53% |
| 5% | 32 | 34% |
| 10% | 16 | 11% |
| 17% *(Chang & Polizzi: 4 of 24)* | 9 | 2% |
| 25% | 6 | 0.2% |

Read precisely. The 21-design null **does** rule something out: a true switch
rate above **13.3%** is excluded at 95% (exact, 1 − 0.05^(1/21)). What it
cannot do is distinguish a 3% design set from a dead one — below **3.2%**,
seeing nothing was more likely than not *even if the designs work*.

So the round-1 null constrains the architecture, not the score. **The useful
round-2 question is how many constructs to order, not which metric to rank by.**
At Chang & Polizzi's measured 17% a 9-design order suffices; at a more sober 5%
it takes 32. Round 2 should decide its order size from the rate it is willing
to miss, and that decision deserves more attention than any further gate
engineering.

---

## 7. The pilot — run this before committing to anything above

*predicted cost*: **3.1 h on 2 GPUs** (LigandMPNN 2.7 h + 20 AF3 predictions,
23 min).

Take 10 of the 71 survivors, spread across their parent templates. Redesign
each with the 40 + 40 window. Keep the single best sequence per parent. Predict
holo and apo. Score with `13_delta_ipae_gate.py --module-map` and
`14_regroup.py`.

It measures all three assumptions directly:

- **A1** — the round-2 construct's own wall clock, ending the 0.41-vs-2.0 band
- **A2** — how many of 10 children keep their parent's PASS, which collapses
  every "low – high" in §5 to a number
- **A3** — whether the 40 + 40 window actually lifts the CONTACT rate from 29.8%

**3.1 hours de-risks up to 12.8 days.** There is no version of this campaign
where skipping it is correct.

---

## 8. Open decisions, for the record

1. **Which G4?** 16-nt GGGT repeat (what round 1 measured) or 22-nt Myc (what
   the round-2 docs assume). Decides the binding module and forces a threshold
   recalibration if it changes. §3.1.
2. **Which scenario?** A (agreed with Kielen, 3.8 d), B-lite (52 h), B (12.8 d),
   or E (35 h plus a rebuild). §5.
3. **Does round 2 keep the matured chromophore?** §2 is the argument against.
4. **Order size.** §6 — pick the switch rate you are willing to miss first.

Items 1–3 are prerequisites to spending the compute in §5. Item 4 is
prerequisite to the compute being worth spending.

---

## 9. Provenance

| File | What it is |
|---|---|
| `14_regroup.py` | assigns each design its first-failed gate; regroups by template |
| `14_regroup_failuremode.csv` | 999 rows, every design with its failure mode and per-gate flags |
| `14_regroup_survivors.csv` | the 71, full metrics, sorted by ΔiPAE |
| `15_round2_forecast.py` | the cost and funnel model; all constants sourced in its header |
| `15_round2_forecast.csv` | 4 scenarios, predicted compute and predicted survivors |
| `13_delta_ipae_module.csv` | upstream — the ΔiPAE gate output these are built from |

All PAE columns are in Ångström. `n_module_contact` is a residue count.
`*_h` columns in the forecast CSV are hours. Rows marked *predicted* are
projections and are labelled as such in the CSV header.

No AlphaFold3, LigandMPNN or wet-lab step was run to produce this document.

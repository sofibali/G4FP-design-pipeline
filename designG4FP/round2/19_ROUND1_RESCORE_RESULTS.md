# Round 1, rescored and re-gated over the whole library

**Cheap analysis only. No LigandMPNN, no AlphaFold3, no model building, no
new structures.** Everything below reads predictions that already existed.

Cost: **~1.2 CPU-hours**, ten templates in parallel, 7m16s each. For comparison
the predictions being read cost 40.6 GPU-hours to make.

Produced by `19_select_top.py` over `19_all_designs_rescored_full.csv`
(`09_rescore_and_select.py` run per template and concatenated) joined to
`14_regroup_failuremode.csv`. Full console output in `19_RESULTS.txt`.

---

## 1. What changed: coverage

| | before | now |
|---|---:|---:|
| Designs rescored | 36 | **1,000** |
| Tested variants validated against | 14 | **21** |
| Gate rows joined to a rescored row | 0 | **999 / 999** |

The rescoring had only ever run on 36 designs because the AF3 predictions were
not findable from the repo (see `18_COLLABORATOR_HANDOFF.md`). They were on the
server the whole time. 25 AF3 models are read per design per state — about
50,000 structures.

**A join gotcha worth keeping.** `14_regroup_failuremode.csv` and
`tested_variants.csv` write `G4FP_des1_cro_mod4_design_0096`;
`09_rescore_and_select.py` writes `des1_cro_mod4_design_0096`, and the design
number is zero-padded inconsistently (`design_6` vs `design_0006`). Joining on
the raw string gives **zero** matches between the 999-row and 36-row tables, and
every downstream rate would have read 0 %. `19_select_top.py` normalises all
three spellings and refuses to report rates silently if the join is incomplete.

---

## 2. The cascade reconciles exactly

| gate | reached | passed | conditional |
|---|---:|---:|---:|
| CONTACT | 999 | 298 | 29.8 % |
| NEGATIVE | 298 | 197 | 66.1 % |
| BIND | 197 | 101 | 51.3 % |
| SWITCH | 101 | 71 | 70.3 % |

**71 pass all four (7.11 %).** First gate failed: CONTACT 701 (70.2 %),
pre-closed 101 (10.1 %), weak BIND 96 (9.6 %), no SWITCH 30 (3.0 %).

These are the same numbers `14_regroup.py` produced, now reached from an
independently re-derived rescoring of every structure. That is a real check, not
a restatement: the two paths share the gate thresholds but not the metric
computation.

---

## 3. The fitness score does not rank — and may be inverted

Now testable at n = 21 instead of n = 14.

| | rank of the two fluorescent designs | of | mean rank | expected by chance |
|---|---|---:|---:|---:|
| at n = 14 | 10, 12 | 14 | 11.0 | 7.5 |
| **at n = 21** | **16, 19** | 21 | **17.5** | **11.0** |

**AUC = 0.158** — the probability that a random fluorescent design outscores a
random dark one. A coin flip is 0.5. Exact permutation **P(AUC ≤ 0.158) = 0.076**.

Stated honestly: **with 2 positives this does not reach significance.** What it
does do is remove any basis for using the score to order anything. The point
estimate is far on the wrong side of chance, it got worse rather than better
when the sample grew from 14 to 21, and the score's top two picks at n = 14 were
`r20_l1_cro_mod0` designs — the template that looks best on every binding metric,
cannot switch, and was dark in all three constructs tested.

The separate hard gates inside `09_rescore_and_select.py` removed **0 designs of
100** on every template: pLDDT ≥ 70, pTM ≥ 0.5, ipTM ≥ 0.3 and the chromophore
geometry window are all satisfied by the entire library. They are not filters.

---

## 4. A diversity ceiling that limits any order from round 1

Only **6 of 10 templates** produce a single survivor, and **47 of the 71 come
from `des1_cro_mod4`**.

| template | survivors / designs |
|---|---:|
| des1_cro_mod4 | **47 / 100** |
| r28_l2_cro_mod2 | 9 / 100 |
| r28_l3_cro_mod2 | 9 / 100 |
| r28_l2_cro_mod3 | 3 / 100 |
| r28_l3_cro_mod0 | 2 / 100 |
| r28_l1_cro_mod2 | 1 / 99 |
| des1_cro_mod0, r20_l1_cro_mod0, r20_l2_cro_mod2, r20_l2_cro_mod3 | **0** |

Capping at 4 per template to avoid betting the order on one backbone yields
**18 designs, not the 32 the power analysis asks for**. Getting to 32 from round
1 means taking roughly two-thirds of the order from `des1_cro_mod4`. That is the
diversity constraint and the power constraint colliding, and neither more
compute nor better scoring resolves it — only new backbones do.

---

## 5. The top set

`19_top_set.csv` — 18 designs, 6 templates. Selection rule, stated rather than
implied: **must pass all four gates; at most 4 per template; within a template,
spread evenly across the ΔPAE range.** The fitness score is reported beside each
design and is **not** used to order anything, for the reason in §3.

`19_survivors_full.csv` holds all 71 with their gate values and sequences.

Section 5 of `19_RESULTS.txt` prints the same survivors ordered by fitness score,
clearly labelled as comparison only. The two lists overlap in **8 of 18**, and
the score-ordered list is 20 of 32 `des1_cro_mod4`.

---

## 6. Caveats

- **70 of 1,000 designs have a local-superposition residual > 0.5 Å**, so their
  `chromophore_rmsd` measures barrel deformation rather than triad motion and
  must not be gated on. Concentrated in two templates (30 and 20).
- The gates still measure binding and switching, not maturation. Ten of the 21
  tested designs passed all four and every one was dark.
- Thresholds are percentiles of this library against the 16-nt
  `GGGTGGGTGGGTGGGT` target. They do not transfer to a 22-nt Myc G4.

---

## 7. What was deliberately not run

No LigandMPNN iteration, no AlphaFold3, no grafting, no new backbones.
`09_rescore_and_select.py` emits `09_round2_ligandmpnn_inputs.csv` as a side
effect; it is a candidate list and was not fed to anything.

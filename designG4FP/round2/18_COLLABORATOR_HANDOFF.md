# Handoff: what to push, and what it unlocks

**Status: nothing has been committed or pushed. This is the proposal.**

Remote: `origin git@github.com:sofibali/G4FP-design-pipeline.git`.
`upstream` is `dauparas/LigandMPNN` — the tool we forked. Never push there.

---

## 1. The headline the collaborator needs first

**The seven structures in `MISSING_STRUCTURES.md` were never missing.** All
seven are on `pollux.fraserlab.com`, dated 2026-03-11, alongside the other 993:

```
output_G4FP_des1_cro_mod4/03_alphafold3_predictions_bound/design_6_bound/design_6_bound_model.cif
output_G4FP_des1_cro_mod4/03_alphafold3_predictions_bound/design_30_bound/...
output_G4FP_des1_cro_mod4/03_alphafold3_predictions_bound/design_38_bound/...
output_G4FP_r28_l2_cro_mod3/03_alphafold3_predictions_bound/design_51_bound/...
output_G4FP_r28_l2_cro_mod3/03_alphafold3_predictions_bound/design_56_bound/...
output_G4FP_r28_l2_cro_mod3/03_alphafold3_predictions_bound/design_94_bound/...
output_G4FP_r28_l3_cro_mod0/03_alphafold3_predictions_bound/design_27_bound/...
```

Nobody was careless. They searched all branches, all 4,872 git objects and a
516-file Drive folder and correctly found zero `.cif`. **`MANIFEST.md` did not
list the AF3 predictions** — 162,077 files, 59.56 GB, the other half of the
campaign by volume. The repo described the data it does not contain, and then
left out the half that mattered. That is fixed in this push.

Do not re-predict anything.

---

## 2. What the push changes for them

| Their analysis | Ran at | Can now run at |
|---|---:|---:|
| Round-2 validation against tested variants | n = 14 | **n = 21** |
| Gate calibration population | 36 designs | **999 designs** |
| Chromophore-burial survey (`14_KNOWN_LIMITATIONS.md` §1) | 36 grafted inputs | **999**, both states |
| Threshold percentiles | subset of a subset | the full measured distribution |

One conclusion changes outright. `14_KNOWN_LIMITATIONS.md` §4 says "two
positives from one template… no gate threshold can be validated on that." That
is still true, but at n=21 a stronger statement is available and it points the
other way: **10 of the 21 tested designs pass all four gates and every one of
the 10 was dark, while both designs that did fluoresce are rejected** —
`design_0072` at BIND, `design_0016` at SWITCH. The cascade is not merely
unvalidated; on the only data that exists it is anti-correlated with the
phenotype. That is the result that should drive round 2, and it was not visible
at n=14.

A second thing to flag: `MISSING_STRUCTURES.md` asks for re-prediction against
"22-nt Myc G4". **Round 1 never used Myc.** All ten templates used the 16-nt
`GGGTGGGTGGGTGGGT`. Every percentile threshold in `13_RESCORING_PLAN.md` is
calibrated against the 16-nt target, so switching to Myc invalidates them. This
needs to be settled before round 2 starts (`17_ROUND2_DESIGN.md` §3.1).

---

## 3. Files to commit

All small — the largest is 162 KB. Nothing here is bulk data.

**New analysis, round-1 regrouping**
| File | What it is |
|---|---|
| `14_regroup.py` | assigns every design exactly one failure mode (first gate failed) |
| `14_regroup_failuremode.csv` | **999 rows** — the full library with its failure mode. The table that replaces n=36. |
| `14_regroup_survivors.csv` | the 71 that pass all four gates |

**The ΔiPAE switch gate**
| File | What it is |
|---|---|
| `13_delta_ipae_gate.py` | mean module↔barrel PAE, apo − holo, both directions, no superposition |
| `13_analyze_delta_ipae.py` | distribution and threshold analysis |
| `13_delta_ipae_module.csv` | per-design ΔiPAE on the corrected module definition |
| `13_delta_ipae_designwindow_SUPERSEDED.csv` | the earlier, wrong module definition — kept so the correction is auditable, not deleted |
| `13_RESCORING_PLAN.md` | gate definitions and where each threshold comes from |

**Round-2 planning**
| File | What it is |
|---|---|
| `17_ROUND2_DESIGN.md` | **the operative round-2 brief** — new design campaign, not an iteration |
| `17_forecast_denovo.py` / `.csv` | cost and funnel for the new-design framing |
| `15_ROUND2_CAMPAIGN_PLAN.md` | the costing record; its scenarios A/B/B-lite are withdrawn, see `17_ROUND2_DESIGN.md` §1 |
| `15_round2_forecast.py` / `.csv` | the iteration cost model those withdrawn scenarios came from |
| `13_G4_BINDING_MODULES.md` | first-pass module survey |
| `16_G4_MODULE_SURVEY.md` | full module survey: graft envelope, rigidity as a class property, the Rap1 hinge, SUD-M rejected, the traps |

**Amended**
| File | Change |
|---|---|
| `../MANIFEST.md` | adds the AF3 prediction inventory — host, path, counts, sizes, regeneration command, and an explicit note that checksums for those 162,077 files are not yet computed |
| `MISSING_STRUCTURES.md` | marked resolved in place, with where each file is |

Suggested commits, each a claim rather than a filename list:

1. `Inventory the AF3 predicted structures the manifest omitted (59.56 GB)`
2. `Resolve MISSING_STRUCTURES: all 21 of 21 are on the server, none re-predicted`
3. `Regroup round 1 by first gate failed: 70% die at CONTACT, not at the switch`
4. `Add the ΔiPAE confidence-contrast gate on the corrected module definition`
5. `Survey graftable G4-binding modules; rigidity is a class property`
6. `Reframe round 2 as a new design campaign; withdraw the iteration scenarios`

---

## 4. What stays out of git, and where it is

Unchanged tier rule: code and small tables are committed, bulk is not.

| Class | Files | Size | Where |
|---|---:|---:|---|
| AF3 predictions (`.cif`, confidences) | 162,077 | 59.56 GB | `pollux.fraserlab.com:/home/sbali/LigandMPNN/designG4FP/output_G4FP_*/03_alphafold3_predictions_{bound,apo}/` |
| LigandMPNN backbones / packed / stats | 245,269 | 59.90 GB | same host, `output_G4FP_*/01_ligandmpnn/` |

Verify a copy against `MANIFEST.sha256` with `sha256sum -c` from the repo root.
That list currently covers the LigandMPNN tier only; the AF3 tier is inventoried
by count and size but not yet checksummed, and `MANIFEST.md` says so.

**Absence from git means gitignored. It has never meant absent.**

---

## 5. Before pushing

- [ ] `git remote -v` — pushing to `origin`, not `upstream`
- [ ] `git diff --cached --stat` — nothing over ~1 MB staged
- [ ] 87 modified files under `07_results/` and `analysis_output/` are the
      corrected-analysis reruns and are a **separate** commit from this handoff;
      they need their own before/after delta reported, not a silent overwrite
- [ ] `MANIFEST.md` numbers match what is on disk
- [ ] no number in `README.md` still reflects the pre-regrouping funnel

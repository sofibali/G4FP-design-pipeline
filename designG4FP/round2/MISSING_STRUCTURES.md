# RESOLVED — the structures were never missing

**Resolved 2026-11-19. Do not act on the request below; no re-prediction is
needed.** Kept rather than deleted so the episode is not rediscovered.

All seven are on `pollux.fraserlab.com`, dated 2026-03-11, at
`designG4FP/output_G4FP_<template>/03_alphafold3_predictions_bound/design_<N>_bound/design_<N>_bound_model.cif`
(`design_0006` → `design_6`). All 21 of the 21 tested variants matched, not 14.

The search that produced this file was correct. `MANIFEST.md` simply did not
inventory the AF3 predictions — 162,077 files, 59.56 GB — so the repo gave no
way to learn they existed. That omission is fixed; see `MANIFEST.md` and
`18_COLLABORATOR_HANDOFF.md`.

One further correction: the request below asks for prediction against a **22-nt
Myc G4**. Round 1 never used Myc. All ten templates used the 16-nt
`GGGTGGGTGGGTGGGT`, which is why every threshold in `13_RESCORING_PLAN.md` is
calibrated against that target.

---

## Original request, superseded


Seven of the 21 experimentally tested round-1 variants have **no AlphaFold3
structure** in the Google Drive project folder or anywhere in this repository.

Searched and confirmed absent:
- Drive `Projects/G4FP` — all 516 `.cif` files checked by design id
- This repo — all branches, all 4,872 objects in git history, zero `.cif` ever
  committed (AF3 predictions are excluded by `designG4FP/.gitignore`)

They exist only on the lab server, if at all. **Sofia: please re-predict these
seven in the holo state (protein + 22-nt Myc G4 + 3 K+), 5 seeds, same settings
as round 1.** Round 2 proceeds with the other 14 meanwhile.

| # | design id | parent template |
|---|---|---|
| 1 | `G4FP_des1_cro_mod4_design_0006` | des1_cro_mod4 |
| 2 | `G4FP_des1_cro_mod4_design_0030` | des1_cro_mod4 |
| 3 | `G4FP_des1_cro_mod4_design_0038` | des1_cro_mod4 |
| 4 | `G4FP_r28_l2_cro_mod3_design_0051` | r28_l2_cro_mod3 |
| 5 | `G4FP_r28_l2_cro_mod3_design_0056` | r28_l2_cro_mod3 |
| 6 | `G4FP_r28_l2_cro_mod3_design_0094` | r28_l2_cro_mod3 |
| 7 | `G4FP_r28_l3_cro_mod0_design_0027` | r28_l3_cro_mod0 |

Sequences for all seven are in
`../07_results/plot_data_csvs/07_all_designs_scored.csv` (column `sequence`),
matched on `global_id`.

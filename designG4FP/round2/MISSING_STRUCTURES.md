# Missing AF3 structures — needed for round 2

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

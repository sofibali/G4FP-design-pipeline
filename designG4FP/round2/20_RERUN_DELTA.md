# Re-run delta: what the corrected analysis changed

Generated 2026-11-22T02:00:13Z by `round2/20_delta_rerun.py`.
Baseline is git `2bb74e6be` (the committed, pre-fix tables); "new" is the
working tree after the corrected re-run. 40 changed CSV tables.

**Read the shifts, not just the new tables.** A corrected table that
replaces the old one in place hides the size of the correction from
anyone who already used the old numbers.

## 1. The master scored table

`designG4FP/07_results/plot_data_csvs/07_all_designs_scored.csv`

1000 rows before, 1000 after, 1000 joined on `global_id`.

### 1.1 Ranking stability of `fitness_score`

| quantity | value |
|---|---:|
| designs compared | 1000 |
| Pearson r (old vs new score) | 0.789 |
| Spearman rho (old vs new ordering) | 0.779 |
| top-10 overlap | 3 / 10 |
| top-50 overlap | 22 / 50 |
| top-100 overlap | 66 / 100 |

### 1.2 Where the two fluorescent designs moved

| design | rank before | rank after | score before | score after |
|---|---:|---:|---:|---:|
| `G4FP_r28_l2_cro_mod2_design_0016` | 139 | 179 | 0.5684 | 0.4540 |
| `G4FP_r28_l2_cro_mod2_design_0072` | 169 | 325 | 0.5524 | 0.3882 |

Ranks are out of 1000. These are the only two designs in the
campaign with a measured positive phenotype, so where they sit
is the only part of the ranking that has an external check.

### 1.3 Selection-flag churn

| flag | before | after | kept | dropped | added |
|---|---:|---:|---:|---:|---:|
| `in_pareto` | 36 | 27 | 20 | 16 | 7 |
| `in_top10_per_tpl` | 100 | 100 | 47 | 53 | 53 |
| `in_top50_global` | 50 | 50 | 22 | 28 | 28 |

### 1.5 Which templates the correction touched

| template | `chromophore_rmsd` | `holo_chromophore_plddt` | `fitness_score` |
|---|---|---|---|
| `G4FP_des1_cro_mod0` | 0 / 100 | 100 / 100 | 100 / 100 |
| `G4FP_des1_cro_mod4` | 0 / 100 | 100 / 100 | 100 / 100 |
| `G4FP_r20_l1_cro_mod0` | 100 / 100 | 100 / 100 | 100 / 100 |
| `G4FP_r20_l2_cro_mod2` | 100 / 100 | 100 / 100 | 100 / 100 |
| `G4FP_r20_l2_cro_mod3` | 100 / 100 | 100 / 100 | 100 / 100 |
| `G4FP_r28_l1_cro_mod2` | 100 / 100 | 100 / 100 | 100 / 100 |
| `G4FP_r28_l2_cro_mod2` | 100 / 100 | 100 / 100 | 100 / 100 |
| `G4FP_r28_l2_cro_mod3` | 100 / 100 | 100 / 100 | 100 / 100 |
| `G4FP_r28_l3_cro_mod0` | 0 / 100 | 100 / 100 | 100 / 100 |
| `G4FP_r28_l3_cro_mod2` | 0 / 100 | 100 / 100 | 100 / 100 |

A column that changed on some templates and not others means the
defect was conditional on something template-specific, so the
untouched templates' earlier numbers still stand.

### 1.4 Columns that moved

| column | n changed | mean before | mean after | shift | r(old,new) | max abs change |
|---|---:|---:|---:|---:|---:|---:|
| `apo_chromophore_plddt` | 1000 | 91.251 | 70.475 | -20.776 | 0.431 | 32.918 |
| `chromophore_apo_plddt` | 1000 | 91.251 | 70.475 | -20.776 | 0.431 | 32.918 |
| `holo_chromophore_plddt` | 1000 | 93.969 | 78.300 | -15.669 | 0.372 | 23.788 |
| `chromophore_holo_plddt` | 1000 | 93.969 | 78.300 | -15.669 | 0.372 | 23.788 |
| `chromophore_plddt_diff` | 1000 | 2.718 | 7.825 | 5.107 | 0.166 | 20.298 |
| `holo_chromophore_plddt_sd` | 1000 | 0.293 | 1.432 | 1.139 | -0.081 | 4.470 |
| `chromophore_holo_plddt_sd` | 1000 | 0.293 | 1.432 | 1.139 | -0.081 | 4.470 |
| `chromophore_plddt_diff_sd` | 1000 | 0.487 | 1.559 | 1.072 | -0.078 | 4.449 |
| `score_component_chrom_holo` | 1000 | 0.659 | 0.320 | -0.339 | 0.372 | 0.791 |
| `apo_mean_plddt` | 1000 | 81.997 | 82.223 | 0.226 | 0.994 | 0.616 |
| `apo_chromophore_plddt_sd` | 1000 | 0.364 | 0.561 | 0.197 | 0.372 | 1.305 |
| `chromophore_apo_plddt_sd` | 1000 | 0.364 | 0.561 | 0.197 | 0.372 | 1.305 |
| `mean_plddt_diff` | 1000 | 2.316 | 2.186 | -0.129 | 0.999 | 0.290 |
| `holo_mean_plddt` | 1000 | 84.313 | 84.409 | 0.097 | 0.996 | 0.467 |
| `fitness_score` | 1000 | 0.469 | 0.381 | -0.087 | 0.789 | 0.214 |
| `score_component_chrom_diff` | 1000 | 0.451 | 0.387 | -0.064 | 0.166 | 0.639 |
| `apo_plddt_vs_template` | 1000 | 0.500 | 0.564 | 0.064 | 0.993 | 0.568 |
| `chromophore_rmsd` | 600 | 0.040 | 0.085 | 0.045 | -0.170 | 0.422 |
| `chromophore_rmsd_sd` | 600 | 0.016 | 0.041 | 0.025 | -0.374 | 0.216 |
| `mean_plddt_diff_sd` | 1000 | 0.564 | 0.550 | -0.014 | 0.998 | 0.060 |
| `apo_mean_plddt_sd` | 1000 | 0.401 | 0.390 | -0.011 | 0.999 | 0.057 |
| `holo_mean_plddt_sd` | 1000 | 0.367 | 0.359 | -0.008 | 0.998 | 0.048 |
| `score_component_confidence` | 998 | 0.799 | 0.803 | 0.004 | 0.998 | 0.032 |
| `score_component_rel_apo` | 998 | 0.560 | 0.563 | 0.004 | 0.993 | 0.061 |
| `final_rank` | 998 | 500.500 | 500.500 | 0.000 | 0.779 | 724.000 |

## 2. Every changed table

One line per table; the column with the largest absolute mean shift.

| table | rows | cols changed | largest shift | mean before | mean after | r |
|---|---:|---:|---|---:|---:|---:|
| `designG4FP/07_results/07_final_candidates.csv` | 1000 | 25 | `apo_chromophore_plddt` -20.776 | 91.251 | 70.475 | 0.431 |
| `designG4FP/07_results/07_pareto_frontier.csv` | 27 | 24 | `apo_chromophore_plddt` -19.501 | 90.763 | 71.262 | 0.270 |
| `designG4FP/07_results/07_top10_per_template.csv` | 100 | 25 | `apo_chromophore_plddt` -21.159 | 91.776 | 70.618 | 0.358 |
| `designG4FP/07_results/07_top50_global.csv` | 50 | 25 | `apo_chromophore_plddt` -20.820 | 90.044 | 69.224 | 0.149 |
| `designG4FP/07_results/plot_data_csvs/07_all_designs_scored.csv` | 1000 | 25 | `apo_chromophore_plddt` -20.776 | 91.251 | 70.475 | 0.431 |
| `designG4FP/07_results/plot_data_csvs/07_pareto_frontier.csv` | 27 | 24 | `apo_chromophore_plddt` -19.501 | 90.763 | 71.262 | 0.270 |
| `designG4FP/07_results/plot_data_csvs/07_selected_designs.csv` | 133 | 25 | `apo_chromophore_plddt` -21.156 | 91.459 | 70.303 | 0.377 |
| `designG4FP/07_results/plot_data_csvs/07_top10_per_template_heatmap_data.csv` | 100 | 5 | `holo_chromophore_plddt` -13.302 | 95.619 | 82.316 | 0.402 |
| `designG4FP/07_results/plot_data_csvs/07_top50_global_heatmap_data.csv` | 50 | 5 | `chromophore_plddt_diff` 12.554 | 4.119 | 16.672 | 0.144 |
| `designG4FP/analysis_output/05_structure_analysis/plot_data_csvs/02_alphafold3_metrics_data.csv` | 200 | 2 | `mean_plddt` 0.120 | 82.121 | 82.241 | 0.998 |
| `designG4FP/analysis_output/05_structure_analysis/plot_data_csvs/04_rmsd_distributions_data.csv` | 200 | 0 | (no numeric change) | | | |
| `designG4FP/analysis_output/05_structure_analysis/plot_data_csvs/05_correlation_matrix_data.csv` | 11 | 0 | (no numeric change) | | | |
| `designG4FP/analysis_output/05_structure_analysis/plot_data_csvs/06_global_vs_chromophore_rmsd_data.csv` | 200 | 0 | (no numeric change) | | | |
| `designG4FP/analysis_output/05_structure_analysis/plot_data_csvs/07_top20_designs_data.csv` | 80 | 22 | `chromophore_mean_plddt` -17.538 | 93.291 | 75.753 | 0.594 |
| `designG4FP/analysis_output/05_structure_analysis_results.csv` | 200 | 6 | `chromophore_mean_plddt` -19.365 | 94.274 | 74.909 | 0.793 |
| `designG4FP/analysis_output/06_ligand_state_comparison/plot_data_csvs/01_holo_apo_rmsd_data.csv` | 100 | 0 | (no numeric change) | | | |
| `designG4FP/analysis_output/06_ligand_state_comparison/plot_data_csvs/02_plddt_comparison_data.csv` | 100 | 6 | `apo_mean_plddt` 0.203 | 80.812 | 81.016 | 0.995 |
| `designG4FP/analysis_output/06_ligand_state_comparison/plot_data_csvs/06_correlation_matrix_data.csv` | 16 | 0 | (no numeric change) | | | |
| `designG4FP/analysis_output/06_ligand_state_comparison/plot_data_csvs/08_top20_heatmap_data.csv` | 20 | 19 | `holo_chromophore_plddt` -17.428 | 95.601 | 78.173 | 0.285 |
| `designG4FP/analysis_output/06_ligand_state_comparison_results.csv` | 100 | 16 | `apo_chromophore_plddt` -20.495 | 92.822 | 72.327 | 0.432 |
| `designG4FP/output_G4FP_des1_cro_mod0/05_structure_analysis_results.csv` | 200 | 6 | `chromophore_mean_plddt` -16.866 | 88.310 | 71.444 | 0.802 |
| `designG4FP/output_G4FP_des1_cro_mod0/06_ligand_state_comparison_results.csv` | 100 | 16 | `apo_chromophore_plddt` -17.519 | 85.604 | 68.085 | 0.260 |
| `designG4FP/output_G4FP_des1_cro_mod4/05_structure_analysis_results.csv` | 200 | 6 | `chromophore_mean_plddt` -17.279 | 92.969 | 75.690 | 0.910 |
| `designG4FP/output_G4FP_des1_cro_mod4/06_ligand_state_comparison_results.csv` | 100 | 16 | `apo_chromophore_plddt` -20.106 | 90.021 | 69.915 | 0.571 |
| `designG4FP/output_G4FP_r20_l1_cro_mod0/05_structure_analysis_results.csv` | 200 | 8 | `chromophore_mean_plddt` -15.545 | 91.333 | 75.788 | 0.708 |
| `designG4FP/output_G4FP_r20_l1_cro_mod0/06_ligand_state_comparison_results.csv` | 100 | 18 | `apo_chromophore_plddt` -21.712 | 90.734 | 69.022 | 0.414 |
| `designG4FP/output_G4FP_r20_l2_cro_mod2/05_structure_analysis_results.csv` | 200 | 8 | `chromophore_mean_plddt` -18.680 | 91.634 | 72.954 | 0.284 |
| `designG4FP/output_G4FP_r20_l2_cro_mod2/06_ligand_state_comparison_results.csv` | 100 | 18 | `apo_chromophore_plddt` -22.443 | 91.692 | 69.249 | 0.721 |
| `designG4FP/output_G4FP_r20_l2_cro_mod3/05_structure_analysis_results.csv` | 200 | 8 | `chromophore_mean_plddt` -20.534 | 91.469 | 70.936 | 0.440 |
| `designG4FP/output_G4FP_r20_l2_cro_mod3/06_ligand_state_comparison_results.csv` | 100 | 18 | `apo_chromophore_plddt` -25.251 | 91.248 | 65.998 | 0.648 |
| `designG4FP/output_G4FP_r28_l1_cro_mod2/05_structure_analysis_results.csv` | 200 | 8 | `chromophore_mean_plddt` -17.519 | 94.601 | 77.082 | 0.707 |
| `designG4FP/output_G4FP_r28_l1_cro_mod2/06_ligand_state_comparison_results.csv` | 100 | 18 | `apo_chromophore_plddt` -20.280 | 92.878 | 72.598 | 0.097 |
| `designG4FP/output_G4FP_r28_l2_cro_mod2/05_structure_analysis_results.csv` | 200 | 8 | `chromophore_mean_plddt` -19.256 | 94.271 | 75.016 | 0.454 |
| `designG4FP/output_G4FP_r28_l2_cro_mod2/06_ligand_state_comparison_results.csv` | 100 | 18 | `apo_chromophore_plddt` -20.935 | 93.173 | 72.238 | 0.081 |
| `designG4FP/output_G4FP_r28_l2_cro_mod3/05_structure_analysis_results.csv` | 200 | 8 | `chromophore_mean_plddt` -19.951 | 94.913 | 74.962 | 0.742 |
| `designG4FP/output_G4FP_r28_l2_cro_mod3/06_ligand_state_comparison_results.csv` | 100 | 18 | `apo_chromophore_plddt` -21.618 | 93.903 | 72.285 | -0.171 |
| `designG4FP/output_G4FP_r28_l3_cro_mod0/05_structure_analysis_results.csv` | 200 | 6 | `chromophore_mean_plddt` -17.232 | 92.326 | 75.094 | 0.746 |
| `designG4FP/output_G4FP_r28_l3_cro_mod0/06_ligand_state_comparison_results.csv` | 100 | 16 | `apo_chromophore_plddt` -17.402 | 90.435 | 73.032 | 0.424 |
| `designG4FP/output_G4FP_r28_l3_cro_mod2/05_structure_analysis_results.csv` | 200 | 6 | `chromophore_mean_plddt` -19.365 | 94.274 | 74.909 | 0.793 |
| `designG4FP/output_G4FP_r28_l3_cro_mod2/06_ligand_state_comparison_results.csv` | 100 | 16 | `apo_chromophore_plddt` -20.495 | 92.822 | 72.327 | 0.432 |

## 3. Tables whose row membership changed

- `designG4FP/07_results/07_pareto_frontier.csv`: 16 rows dropped, 7 added
- `designG4FP/07_results/07_top10_per_template.csv`: 53 rows dropped, 53 added
- `designG4FP/07_results/07_top50_global.csv`: 28 rows dropped, 28 added
- `designG4FP/07_results/plot_data_csvs/07_pareto_frontier.csv`: 16 rows dropped, 7 added
- `designG4FP/07_results/plot_data_csvs/07_selected_designs.csv`: 67 rows dropped, 71 added
- `designG4FP/07_results/plot_data_csvs/07_top10_per_template_heatmap_data.csv`: 53 rows dropped, 53 added
- `designG4FP/07_results/plot_data_csvs/07_top50_global_heatmap_data.csv`: 28 rows dropped, 28 added
- `designG4FP/analysis_output/05_structure_analysis/plot_data_csvs/07_top20_designs_data.csv`: 46 rows dropped, 37 added
- `designG4FP/analysis_output/06_ligand_state_comparison/plot_data_csvs/08_top20_heatmap_data.csv`: 16 rows dropped, 16 added

## 4. Biggest individual rank moves

| design | rank before | rank after | delta |
|---|---:|---:|---:|
| `G4FP_r20_l2_cro_mod2_design_0085` | 942 | 218 | -724 |
| `G4FP_r20_l1_cro_mod0_design_0054` | 852 | 156 | -696 |
| `G4FP_r20_l2_cro_mod3_design_0039` | 923 | 283 | -640 |
| `G4FP_r20_l2_cro_mod2_design_0055` | 824 | 228 | -596 |
| `G4FP_r20_l1_cro_mod0_design_0000` | 849 | 306 | -543 |
| `G4FP_r20_l2_cro_mod2_design_0014` | 791 | 256 | -535 |
| `G4FP_r20_l2_cro_mod2_design_0027` | 759 | 226 | -533 |
| `G4FP_r20_l2_cro_mod2_design_0069` | 865 | 334 | -531 |
| `G4FP_r20_l2_cro_mod2_design_0025` | 816 | 291 | -525 |
| `G4FP_r20_l2_cro_mod2_design_0051` | 794 | 281 | -513 |

## 5. How to reproduce this

```
python3 -I round2/20_delta_rerun.py \
    --repo <repo root> \
    --tested round2/tested_variants.csv \
    --out round2/20_RERUN_DELTA.md
```

It reads the pre-fix tables from git and the corrected ones from the
working tree, so it only works while the re-run is still uncommitted.
After the commit, compare against the parent: `git show <commit>^:<path>`.


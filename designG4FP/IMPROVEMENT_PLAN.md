# G4FP Pipeline Improvement Plan

Assessment of the current pipeline and concrete improvements to make it more robust, easier to run, and better suited to the goal of selecting 100-1000 sequences with maximal bound-vs-apo stability difference.

## Current State Assessment

### What works well
- Clean numbered step convention (01, 02, 03, ...) makes the pipeline order obvious
- YAML config centralizes parameters
- All analysis scripts export CSV for external plotting (Prism/GraphPad)
- LigandMPNN step is complete: 100,000 sequences across 10 templates
- Pre-filter analysis catches problems before expensive AF3 runs

### What needs improvement

**1. Pipeline is not end-to-end automated**
- `run_pipeline.sh` is interactive (requires `read -p` prompts) -- cannot be submitted as a batch job
- The config file is per-structure, so running all 10 structures requires manually editing `pipeline_config.yaml` each time or using the separate `batch_run_ligandmpnn.sh`
- Steps 2-6 have no batch wrapper for all 10 structures (only step 1 has `batch_run_ligandmpnn.sh`)

**2. Five of 10 structures are missing filtering/AF3 prep**
- `output_G4FP_r20_l2_cro_mod3`, `r28_l2_cro_mod2`, `r28_l2_cro_mod3`, `r28_l3_cro_mod0`, `r28_l3_cro_mod2` have LigandMPNN output but no filtered sequences or AF3 JSONs

**3. No final aggregation step**
- After running analysis per-structure, there is no script that pools results across all 10 templates and selects the top 100-1000 overall designs
- The bound-vs-apo comparison (`06_compare_ligand_states.py`) works per-structure but doesn't merge across structures

**4. Step numbering is inconsistent**
- Script files: 01, 01b, 02, 03, 04a, 05, 06
- `run_pipeline.sh` calls them: steps 1, 1b, 2, 3, 4a, 4b, 4c
- The 05/06 scripts are called as "step 4b/4c" in the master script

**5. Too many fix/utility scripts cluttering the directory**
- `fix_af3_json_ids.py`, `fix_af3_json_dialect.py`, `fix_af3_json_structure.py`, `fix_af3_json.py` -- four JSON fixers
- `regenerate_all_ligandmpnn.sh`, `run_missing_ligandmpnn.sh`, `quick_fix_and_resume.sh` -- ad hoc recovery scripts
- `test_ligandmpnn_fix.sh`, `test_analysis_setup.py` -- test scripts
- These were useful during development but make the directory hard to navigate

**6. No single-command batch execution**
- Cannot do `./run_all.sh` to process all 10 structures end-to-end without interaction

**7. AF3 JSON generation is fragile**
- Multiple fix scripts suggest the JSON format has been a recurring pain point
- The `id` field, `dialect` field, and structure format have all needed post-hoc fixes

---

## Improvement Plan

### Phase 1: Immediate fixes (complete the current run)

#### 1.1 Batch-filter the remaining 5 structures
Create a simple loop script that runs `02_filter_and_prepare_af3.py` for all structures that are missing filtered output. This unblocks AF3 for the remaining 5.

```bash
# batch_filter_all.sh
for pdb in inputs/*.pdb; do
    name=$(basename "$pdb" .pdb)
    outdir="output_${name}"
    if [ ! -f "$outdir/02_filtered_sequences.csv" ]; then
        # generate temp config and run filter
        ...
    fi
done
```

#### 1.2 Validate all AF3 JSONs before submitting
Add a `validate_af3_jsons.py` script that checks all JSON files have the correct structure (`id` fields, `dialect`, sequence types) before the expensive AF3 step. This replaces the need for multiple fix scripts after the fact.

#### 1.3 Run AF3 predictions
With all 10 structures filtered, submit parallel AF3 jobs via `run_all_af3_parallel.sh`.

### Phase 2: Add the missing final selection step

#### 2.1 Create `07_aggregate_and_rank.py`
This is the most critical missing piece. After per-structure analysis, this script should:

1. Load `06_ligand_state_comparison_results.csv` from all 10 output directories
2. Merge into a single DataFrame with columns:
   - `template` (which of the 10 input structures)
   - `design_id`
   - `sequence`
   - `bound_plddt`, `apo_plddt`, `delta_plddt`
   - `bound_rmsd`, `apo_rmsd`, `delta_rmsd`
   - `bound_ptm`, `apo_ptm`
   - `chromophore_bound_plddt`, `chromophore_apo_plddt`
   - `ligandmpnn_confidence`, `ligandmpnn_ligand_confidence`
3. Compute a **G4FP fitness score**:
   ```
   score = w1 * delta_plddt + w2 * delta_rmsd + w3 * bound_ptm + w4 * (1 - apo_ptm)
   ```
   Where delta = bound - apo. Higher score = bigger ligand-dependent stabilization.
4. Apply hard filters:
   - bound_plddt > 70 (must fold when bound)
   - bound_ptm > 0.5 (plausible fold)
   - apo_plddt < 60 (should be disordered without ligand)
5. Rank by fitness score, select top 100-1000
6. Export:
   - `07_final_candidates.csv` -- ranked list with all metrics
   - `07_final_candidates.fa` -- FASTA of selected sequences
   - `07_selection_summary.png` -- visualization of score distributions
   - `07_bound_vs_apo_scatter.png` -- bound pLDDT vs apo pLDDT colored by rank

#### 2.2 Create `08_export_for_synthesis.py`
Format final sequences for gene synthesis ordering:
- Add codon optimization (or flag for external tools)
- Add cloning overhangs/restriction sites
- Export in plate-map format if needed
- Generate a summary table for lab notebook

### Phase 3: Streamline the pipeline structure

#### 3.1 Consolidate to a single `run_all.sh` non-interactive batch script
```bash
#!/bin/bash
# run_all.sh -- Run the complete G4FP design pipeline for all input structures
# Usage: nohup ./run_all.sh [--parallel N] &

PARALLEL=${1:-3}

# Step 1: LigandMPNN (all structures)
./batch_run_ligandmpnn.sh

# Step 1b: Pre-filter analysis (all structures)
for dir in output_G4FP_*/; do
    python 01b_analyze_pre_filter.py --output-dir "$dir" &
done
wait

# Step 2: Filter and prepare AF3 (all structures)
for dir in output_G4FP_*/; do
    python 02_filter_and_prepare_af3.py --output-dir "$dir"
done

# Step 3: AF3 predictions (parallel)
./run_all_af3_parallel.sh $PARALLEL

# Step 4-6: Analysis (all structures)
./run_structure_analysis.sh

# Step 7: Aggregate and rank
python 07_aggregate_and_rank.py
```

#### 3.2 Make all scripts accept `--output-dir` directly
Currently `02_filter_and_prepare_af3.py` reads from the YAML config. It should also accept `--output-dir` and `--input-pdb` as CLI args so the batch wrapper doesn't need to rewrite the config file for each structure. The config file should be a fallback default, not required.

#### 3.3 Move utility/fix scripts to a `utils/` subdirectory
```
designG4FP/
├── utils/
│   ├── fix_af3_json_ids.py
│   ├── fix_af3_json_dialect.py
│   ├── fix_af3_json_structure.py
│   ├── fix_af3_json.py
│   ├── validate_af3_jsons.py
│   ├── check_pipeline_status.py
│   ├── test_analysis_setup.py
│   └── test_ligandmpnn_fix.sh
```

#### 3.4 Move ad hoc scripts to an `archive/` subdirectory
```
designG4FP/
├── archive/
│   ├── regenerate_all_ligandmpnn.sh
│   ├── run_missing_ligandmpnn.sh
│   ├── quick_fix_and_resume.sh
│   ├── batch_run_ligandmpnn.sh.bak
│   ├── LIGANDMPNN_FIX_SUMMARY.md
│   ├── AF3_FIX_SUMMARY.md
│   ├── REGENERATE_GUIDE.md
│   ├── PIPELINE_STATUS_SUMMARY.md
│   └── *.log
```

#### 3.5 Consolidate documentation
Current state: 8+ markdown files (`PIPELINE_README.md`, `IMPLEMENTATION_SUMMARY.md`, `PIPELINE_STRUCTURE.md`, `STRUCTURE_ANALYSIS_README.md`, `ANALYSIS_SUMMARY.md`, `AF3_FIX_SUMMARY.md`, `LIGANDMPNN_FIX_SUMMARY.md`, `REGENERATE_GUIDE.md`, `PIPELINE_STATUS_SUMMARY.md`).

Consolidate to:
- `README.md` -- main documentation (created)
- `IMPROVEMENT_PLAN.md` -- this file
- Archive the rest to `archive/docs/`

### Phase 4: Robustness improvements

#### 4.1 Add error handling to `02_filter_and_prepare_af3.py`
- Validate that the input PDB actually contains DNA chain and ions before generating AF3 JSONs
- Validate JSON output with a schema check before writing
- Report if any sequences have unusual lengths or compositions

#### 4.2 Add AF3 output validation
After AF3 predictions complete, automatically check:
- CIF files are parseable
- Confidence JSON files exist and contain expected fields
- Flag failed/incomplete predictions before analysis

#### 4.3 Make chromophore residues auto-detected
Currently `chromophore_start: 197` and `chromophore_end: 199` are hardcoded in the config. These could be auto-detected from the input PDB by searching for the CRO residue.

#### 4.4 Add a Snakemake or Nextflow workflow (optional, longer-term)
For full reproducibility and dependency tracking, the pipeline could be expressed as a Snakemake workflow where each step's inputs/outputs are declared explicitly. This would handle:
- Automatic re-runs of failed steps
- Parallelization across structures
- No re-computation of completed steps
- Provenance tracking

### Phase 5: Scientific improvements

#### 5.1 Explore alternative filtering strategies
Current: top 1% by `overall_confidence` alone.
Consider:
- **Pareto filtering**: select sequences on the Pareto frontier of overall_confidence vs ligand_confidence (both matter)
- **Diversity-weighted selection**: cluster sequences by similarity, pick top N from each cluster to avoid redundancy
- **Multi-temperature sampling**: run LigandMPNN at multiple temperatures (T=0.1, 0.2, 0.5) and combine pools

#### 5.2 Add ESMFold as a fast pre-screen
Before expensive AF3 runs, use ESMFold (single-sequence structure prediction, ~seconds per sequence) to do a quick structural screen:
- Predict structure from sequence alone (no DNA context = simulates apo)
- Filter out sequences that fold well without ligand (those are NOT what we want)
- Only send sequences with low ESMFold confidence to AF3 (these are more likely to be ligand-dependent)
- This could reduce AF3 load by 50-80%

#### 5.3 Use AF3 pAE for interface quality
The current analysis uses pTM and pLDDT but doesn't fully exploit the predicted aligned error (PAE) matrix. The protein-DNA interface PAE could be a better discriminator of ligand binding quality than ipTM alone.

#### 5.4 Consider multiple AF3 seeds for consensus
Currently using 5 model seeds per design. Designs where all 5 models agree (low variance in pLDDT/RMSD across seeds) are more reliable than designs where only 1/5 models looks good. Add a **consensus score** that penalizes high inter-seed variance.

---

## Priority Order

| Priority | Task | Effort | Impact |
|----------|------|--------|--------|
| **P0** | Batch-filter remaining 5 structures | 30 min | Unblocks AF3 |
| **P0** | Run AF3 predictions | Days (GPU) | Core pipeline |
| **P1** | Create `07_aggregate_and_rank.py` | 2-3 hrs | **Critical** -- the whole point |
| **P1** | Make scripts accept `--output-dir` CLI args | 1-2 hrs | Enables batch runs |
| **P2** | Create `run_all.sh` non-interactive script | 1 hr | Ease of use |
| **P2** | Move fix/utility scripts to `utils/` | 30 min | Clean directory |
| **P2** | Consolidate docs | 30 min | Clarity |
| **P3** | Add AF3 JSON validation | 1 hr | Prevents failures |
| **P3** | Add ESMFold pre-screen | 3-4 hrs | Reduces AF3 cost |
| **P3** | Diversity-weighted filtering | 2 hrs | Better candidates |
| **P4** | Snakemake workflow | 1-2 days | Reproducibility |
| **P4** | Export for synthesis script | 1-2 hrs | Lab integration |

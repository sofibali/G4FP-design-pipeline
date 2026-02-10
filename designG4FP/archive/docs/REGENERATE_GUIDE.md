# Quick Start: Regenerate LigandMPNN Outputs

## Problem Fixed
✅ LigandMPNN now correctly designs **92 residues** per structure  
✅ Generates **diverse sequences** (not identical)  
✅ Uses correct conda environment (`unified_mpnn`)

## How to Regenerate All 10 Structures

### Option 1: Interactive Script (Recommended)
```bash
cd /home/sbali/LigandMPNN/designG4FP
bash regenerate_all_ligandmpnn.sh
```

This will:
- Ask for confirmation
- Delete all existing LigandMPNN outputs
- Run batch generation for all 10 structures
- Run in background with logging

### Option 2: Manual Batch Run
```bash
cd /home/sbali/LigandMPNN/designG4FP

# Force regenerate all
bash batch_run_ligandmpnn.sh --force-rerun

# Or run in background
nohup bash batch_run_ligandmpnn.sh --force-rerun &
```

### Option 3: Single Structure
```bash
cd /home/sbali/LigandMPNN/designG4FP

# Delete old output
rm -rf output_G4FP_des1_cro_mod0/01_ligandmpnn

# Edit pipeline_config.yaml to point to your PDB
bash 01_run_ligandmpnn.sh pipeline_config.yaml
```

## Monitor Progress

### Check progress log:
```bash
tail -f batch_ligandmpnn_progress.log
```

### Check detailed output:
```bash
tail -f batch_ligandmpnn_detailed.log
```

### Check if job is still running:
```bash
ps aux | grep batch_run_ligandmpnn
```

## Expected Results

For each structure, you should see:
- **~10,000 sequences generated**
- **High sequence diversity** (thousands of unique sequences)
- **Varying confidence scores** (not all 1.0)
- **92 residues designed** per sequence

## Time Estimates
- Per structure: ~20-30 minutes
- All 10 structures: ~3-5 hours

## Verify Success

After completion, check one structure:
```bash
# Count sequences
grep -c "^>" output_G4FP_des1_cro_mod0/01_ligandmpnn/seqs/*.fa

# Check diversity (should show num_res=92)
head -5 output_G4FP_des1_cro_mod0/01_ligandmpnn/seqs/*.fa

# Run analysis
bash run_all_1b_analyses.sh
```

## Files Updated
- ✅ [`01_run_ligandmpnn.sh`](01_run_ligandmpnn.sh) - Fixed residue expansion + conda env
- ✅ [`batch_run_ligandmpnn.sh`](batch_run_ligandmpnn.sh) - Added conda env + force rerun flag
- ✅ [`regenerate_all_ligandmpnn.sh`](regenerate_all_ligandmpnn.sh) - New convenience script

## Next Steps After Regeneration

1. Run pre-filter analysis:
   ```bash
   bash run_all_1b_analyses.sh
   ```

2. Verify diversity in analysis plots

3. Continue with pipeline (filter & AlphaFold3)

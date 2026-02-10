# Comprehensive Structure Analysis for designG4FP

This directory contains comprehensive analysis scripts for AlphaFold3 predictions, matching the analysis paradigms from the previous pipeline (`/home/sbali/LigandMPNN/pipeline_outputs`).

## Analysis Scripts

### 05_analyze_af3_structures.py
Comprehensive individual structure analysis including:
- Overall confidence distributions (LigandMPNN scores)
- AlphaFold3 quality metrics (pTM, ipTM, pLDDT, PAE)
- RMSD analysis (global and chromophore-specific)
- Confidence-pLDDT correlations
- Per-residue analysis
- Top design identification

Similar to previous pipeline's `07_structure_analysis/`.

### 06_compare_ligand_states.py
Ligand-bound vs apo state comparison including:
- RMSD between bound and apo structures
- pLDDT differences between states
- Per-residue conformational changes
- Chromophore region analysis
- Comprehensive visualizations

Similar to previous pipeline's `08_ligand_state_comparison/`.

## Quick Start

### Analyze a Single Output Directory

```bash
# If you know the template structure
bash run_structure_analysis.sh output_G4FP_des1_cro_mod0 ../inputs/G4FP_design1_cro.pdb

# Or run scripts individually:
python 05_analyze_af3_structures.py \
    --output-dir output_G4FP_des1_cro_mod0 \
    --template ../inputs/G4FP_design1_cro.pdb \
    --state both \
    --chromophore-range 175,225

python 06_compare_ligand_states.py \
    --output-dir output_G4FP_des1_cro_mod0 \
    --template ../inputs/G4FP_design1_cro.pdb \
    --chromophore-range 175,225
```

### Analyze All Output Directories

```bash
# Analyze all output_* directories
bash run_structure_analysis.sh ../inputs/G4FP_design1_cro.pdb
```

## Output Structure

For each output directory, the analysis creates:

```
output_X/
├── 05_structure_analysis_results.csv              # Combined metrics CSV
├── 05_structure_analysis/
│   ├── 01_confidence_distributions.png            # LigandMPNN score distributions
│   ├── 02_alphafold3_metrics.png                  # pTM, ipTM, pLDDT, ranking_score
│   ├── 03_confidence_vs_plddt.png                 # Correlation plots
│   ├── 04_rmsd_distributions.png                  # Global and chromophore RMSD
│   ├── 05_correlation_matrix.png                  # All metrics correlation heatmap
│   ├── 06_global_vs_chromophore_rmsd.png         # RMSD scatter
│   ├── 07_top20_designs_data.xlsx                 # Top designs by different metrics
│   ├── 08_per_residue_rmsd_design_XXXX.png       # Per-residue analysis for top design
│   └── plot_data_csvs/                            # All plot data as CSV
│       ├── 01a_overall_confidence_data.csv
│       ├── 01b_ligand_confidence_data.csv
│       ├── 02_alphafold3_metrics_data.csv
│       ├── 03_confidence_vs_plddt_data.csv
│       ├── 03b_chromophore_confidence_vs_plddt_data.csv
│       ├── 04_rmsd_distributions_data.csv
│       ├── 05_correlation_matrix_data.csv
│       ├── 06_global_vs_chromophore_rmsd_data.csv
│       ├── 07_top20_designs_data.csv
│       └── 08_per_residue_rmsd_design_XXXX_data.csv
│
├── 06_ligand_state_comparison_results.csv         # Bound vs apo comparison CSV
└── 06_ligand_state_comparison/
    ├── 01_bound_apo_rmsd.png                      # RMSD distributions
    ├── 02_plddt_comparison.png                    # pLDDT scatter and differences
    ├── 03_chromophore_plddt_comparison.png        # Chromophore-specific pLDDT
    ├── 04_rmsd_vs_plddt_change.png               # Correlation plots
    ├── 05_top20_conformational_changes.png        # Top designs by conformational change
    ├── 06_correlation_matrix.png                  # Comparison metrics heatmap
    └── plot_data_csvs/
        ├── 01_bound_apo_rmsd_data.csv
        ├── 02_plddt_comparison_data.csv
        ├── 03_chromophore_plddt_comparison_data.csv
        ├── 04_rmsd_vs_plddt_change_data.csv
        ├── 05_top20_conformational_changes_data.csv
        └── 06_correlation_matrix_data.csv
```

## Requirements

```bash
# Activate environment
conda activate unified_mpnn  # or your Python environment

# Required packages
pip install biopython pandas numpy matplotlib seaborn scipy openpyxl
```

## Configuration

### Template Structure
The template structure is used as reference for RMSD calculations. For designG4FP, use:
- `G4FP_design1_cro.pdb` - For design1 variants
- Or the appropriate backbone from `output_X/01_ligandmpnn/backbones/`

### Chromophore Range
Default: residues 175-225

Adjust if your chromophore region is different:
```bash
--chromophore-range START,END
```

### Analysis States
- `--state bound`: Only analyze bound (with chromophore) predictions
- `--state apo`: Only analyze apo (without chromophore) predictions
- `--state both`: Analyze both states (default)

## Current Status (as of pipeline check)

### Completed AF3 Predictions
- ✓ `output_G4FP_des1_cro_mod0` - Both bound and apo (200 designs each)

### Pending AF3 Predictions
- ⏳ `output_G4FP_des1_cro_mod4` - 200 bound, 200 apo (ready to run)
- ⏳ `output_G4FP_r20_l1_cro_mod0` - 200 bound, 200 apo (ready to run)
- ⏳ `output_G4FP_r20_l2_cro_mod2` - 200 bound, 200 apo (ready to run)

Run AF3 predictions first:
```bash
# Option 1: Run all in parallel
bash run_all_af3_parallel.sh 3  # 3 parallel jobs

# Option 2: Run specific structure
bash 03_run_alphafold3_parallel.sh output_G4FP_des1_cro_mod4/02_alphafold3_inputs_bound/config.json 3
```

## Comparison to Previous Pipeline

### Previous: 07_structure_analysis
Generated:
- Confidence histograms
- AF3 quality metrics
- RMSD distributions
- Confidence-pLDDT correlations (labeled and unlabeled)
- Chromophore-specific analysis
- Per-residue RMSD
- Top20 designs summary
- Correlation heatmap

### Current: 05_analyze_af3_structures.py
Generates **identical analysis** adapted for designG4FP structure:
- Works with `output_X/03_alphafold3_predictions_{bound,apo}/` directories
- Handles both bound and apo states
- Integrates with LigandMPNN `parsed_sequences.fasta` scores
- Exports all plot data as CSV for Prism/GraphPad

### Previous: 08_ligand_state_comparison
Generated:
- Bound/apo RMSD distributions
- pLDDT comparison
- Chromophore pLDDT comparison
- RMSD vs pLDDT change
- Conformational changes
- Correlation heatmap

### Current: 06_compare_ligand_states.py
Generates **identical analysis** adapted for G4FP chromophore:
- Direct structure-to-structure comparison (bound vs apo)
- Per-design conformational change analysis
- Chromophore region focus
- Top designs by largest conformational change

## Usage Examples

### Example 1: Analyze completed predictions
```bash
cd /home/sbali/LigandMPNN/designG4FP

# Single directory
bash run_structure_analysis.sh output_G4FP_des1_cro_mod0 ../inputs/G4FP_design1_cro.pdb

# Results will be in:
# - output_G4FP_des1_cro_mod0/05_structure_analysis/
# - output_G4FP_des1_cro_mod0/06_ligand_state_comparison/
```

### Example 2: Only analyze bound state
```bash
python 05_analyze_af3_structures.py \
    --output-dir output_G4FP_des1_cro_mod0 \
    --template ../inputs/G4FP_design1_cro.pdb \
    --state bound

# No ligand comparison will run (need both states)
```

### Example 3: Custom chromophore range
```bash
bash run_structure_analysis.sh output_G4FP_des1_cro_mod0 ../inputs/G4FP_design1_cro.pdb

# Edit run_structure_analysis.sh line 22:
CHROMOPHORE_RANGE="197,199"  # GFP chromophore (Ser65-Tyr66-Gly67)
```

### Example 4: Batch analysis after AF3 completion
```bash
# Run all pending AF3 predictions first
bash run_all_af3_parallel.sh 3

# Once complete, analyze all
bash run_structure_analysis.sh ../inputs/G4FP_design1_cro.pdb
```

## Troubleshooting

### Issue: "No model file found"
**Cause**: AF3 predictions haven't completed yet
**Fix**: Check `03_alphafold3_predictions_{bound,apo}/design_XXXX_{bound,apo}/` for `*_model.cif` files

### Issue: "No confidence file found"
**Cause**: AF3 output incomplete or different naming
**Fix**: Check for `summary_confidences_0.json` in design directories

### Issue: "Template structure not found"
**Cause**: Wrong path to template
**Fix**: Use absolute path: `/home/sbali/LigandMPNN/inputs/G4FP_design1_cro.pdb`

### Issue: Length mismatch warnings
**Cause**: Template and model have different lengths (expected for multi-state designs)
**Fix**: Warnings are informational only - analysis uses minimum length

## Data Export for Prism/GraphPad

All plot data is exported as CSV in `plot_data_csvs/` directories:
- Import CSVs into Prism for publication-quality figures
- Use correlation matrices for statistical analysis
- Top designs lists for experimental validation

## Citation

If using these scripts, please cite the original pipeline tools:
- AlphaFold3: Abramson et al. Nature 2024
- LigandMPNN: Dauparas et al. Science 2022
- BioPython: Cock et al. Bioinformatics 2009

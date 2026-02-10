# Simplified Protein Design Pipeline

A modular, command-line pipeline for protein sequence design using LigandMPNN and structure prediction with AlphaFold3.

## Pipeline Overview

```
Input PDB → LigandMPNN → Pre-Filter → Filter Top 10% → AlphaFold3 → Analysis
            (design)      Analysis    (by confidence)   (bound/apo)   (metrics)
```

### Steps

1. **LigandMPNN**: Generate diverse sequences with ligand context
2. **Pre-Filter Analysis**: Analyze all sequences before filtering (NEW!)
3. **Filter**: Select top 10% by confidence, prepare AF3 inputs (with/without ligands)
4. **AlphaFold3**: Predict structures for bound and apo states
5. **Analysis**:
   - 5a: Sequence diversity and confidence metrics
   - 5b: Independent structure analysis (RMSD, pLDDT)
   - 5c: Comparison of bound vs apo states

---

## Quick Start

### 1. Configure Pipeline

Edit `pipeline_config.yaml`:

```yaml
input_pdb: "./inputs/your_protein.pdb"
output_dir: "./pipeline_outputs"

ligandmpnn:
  design_residues: "A15-60 A240-285"  # Residues to redesign
  num_sequences: 10000                # Total sequences

filtering:
  top_percent: 10                     # Keep top 10%

alphafold3:
  db_dir: "/path/to/af3/databases"    # UPDATE THIS
  model_dir: "/path/to/af3/models"    # UPDATE THIS
  gpu_id: "0"
```

### 2. Run Full Pipeline

```bash
# Interactive mode (prompts for each step)
./run_pipeline.sh

# With custom config
./run_pipeline.sh my_config.yaml
```

### 3. Run Individual Steps

```bash
# Step 1: LigandMPNN sequence design
./01_run_ligandmpnn.sh

# Step 1b: Analyze all sequences before filtering (NEW!)
python 01b_analyze_pre_filter.py

# Step 2: Filter and prepare AF3 inputs
python 02_filter_and_prepare_af3.py

# Step 3: Run AlphaFold3
./03_run_alphafold3.sh              # Interactive
./03_run_alphafold3.sh bound        # Bound only
./03_run_alphafold3.sh apo          # Apo only
./03_run_alphafold3.sh both         # Both states

# Step 4a: Analyze LigandMPNN sequences
python 04a_analyze_ligandmpnn.py

# Step 4b: Analyze AF3 structures
python analyze_af3_structures.py \
    --pipeline-dir ./pipeline_outputs \
    --af3-dir ./pipeline_outputs/03_alphafold3_predictions_bound \
    --template ./inputs/your_protein.pdb \
    --output-dir ./pipeline_outputs/04b_structure_analysis_bound

# Step 4c: Compare bound vs apo
python compare_ligand_states.py \
    --bound-dir ./pipeline_outputs/03_alphafold3_predictions_bound \
    --apo-dir ./pipeline_outputs/03_alphafold3_predictions_apo \
    --output-dir ./pipeline_outputs/04c_ligand_state_comparison
```

---

## Output Structure

```
pipeline_outputs/
├── 01_ligandmpnn/
│   └── seqs/                          # Generated FASTA files
├── 01b_pre_filter_analysis/           # Analysis before filtering (NEW!)
│   ├── 01_overall_confidence_distribution.png
│   ├── 02_ligand_confidence_distribution.png
│   ├── 03_sequence_diversity_per_position.png
│   ├── 04_ligand_vs_overall_confidence.png
│   └── pre_filter_summary.csv
├── 02_filtered_sequences.csv          # Top 10% sequences
├── 02_filtered_sequences.fa           # FASTA of filtered sequences
├── 02_alphafold3_inputs_bound/        # JSON inputs (protein+DNA)
├── 02_alphafold3_inputs_apo/          # JSON inputs (protein only)
├── 03_alphafold3_predictions_bound/   # AF3 structures (bound)
│   └── design_XXXX/
│       ├── fold_*.cif                 # Structure files
│       └── fold_*_summary_confidences.json
├── 03_alphafold3_predictions_apo/     # AF3 structures (apo)
├── 04a_ligandmpnn_analysis/           # Sequence analysis
│   ├── 01_confidence_distribution.png
│   ├── 03_positional_entropy.png
│   └── plot_data_csvs/                # CSV exports for Prism
├── 04b_structure_analysis_bound/      # Structure metrics (bound)
│   ├── structure_analysis_results.csv
│   ├── 01_plddt_distribution.png
│   ├── 06_global_vs_chromophore_rmsd.png
│   └── plot_data_csvs/
├── 04b_structure_analysis_apo/        # Structure metrics (apo)
└── 04c_ligand_state_comparison/       # Bound vs apo comparison
    ├── ligand_state_comparison_results.csv
    ├── 02_plddt_comparison.png
    ├── 04_rmsd_vs_plddt_change.png
    └── plot_data_csvs/
```

---

## Dependencies

### Python Packages
```bash
pip install biopython pandas numpy matplotlib seaborn pyyaml
```

### Software
- **LigandMPNN**: Already in this directory
- **AlphaFold3**: Requires SBGrid installation or standalone AF3
- **PyYAML**: For config file parsing

---

## Configuration Details

### Input PDB Requirements
- Must contain protein chain A
- Optional: DNA/RNA chains and ions
- Standard PDB format

### Design Residues Format
- Space-separated ranges: `"A15-60 A240-285"`
- Chain letter + residue numbers
- Can specify multiple regions

### AlphaFold3 Settings
- `modelseeds`: Number of predictions per design (default: 5)
- `create_ligand_free`: Generate apo structures (default: true)

---

## Analysis Outputs

### 1b: Pre-Filter Analysis (NEW!)
- **Overall confidence distribution**: All sequences before filtering
- **Ligand confidence distribution**: Ligand binding confidence scores
- **Sequence diversity per position**: Shannon entropy at each residue
- **Confidence correlation**: Ligand vs overall confidence, colored by sequence identity

### 4a: LigandMPNN Analysis
- **Confidence distribution**: Overall and ligand-specific scores
- **Sequence diversity**: Positional entropy, pairwise identity
- **Template comparison**: Identity to input structure

### 4b: Structure Analysis
- **pLDDT scores**: Global and chromophore-specific confidence
- **RMSD metrics**: Alignment to template structure
- **Per-residue analysis**: Deviation and confidence profiles
- **Composite scores**: Combined ranking metric

### 4c: Bound vs Apo Comparison
- **Structural changes**: RMSD between states
- **Confidence changes**: pLDDT difference (bound - apo)
- **Chromophore analysis**: Local flexibility
- **Top conformational changes**: Designs with largest differences

---

## Tips and Best Practices

### 1. LigandMPNN
- Use `num_sequences: 10000` for good diversity
- Adjust `batch_size` based on GPU memory
- `cutoff_for_score: 8.0` defines ligand interaction distance

### 2. Filtering
- `top_percent: 10` typically gives 1000 sequences for 10K designs
- Can adjust to get desired number of AF3 predictions
- Lower percent = higher quality but less diversity

### 3. AlphaFold3
- Each prediction takes ~30min-2hrs on GPU
- Consider running bound and apo separately
- Monitor GPU memory usage
- Use `modelseeds: [0, 1, 2]` for faster testing (3 models instead of 5)

### 4. Analysis
- Run 4a early to check sequence quality
- 4b and 4c require completed AF3 predictions
- All analyses export CSV files for external plotting

---

## Troubleshooting

### LigandMPNN Fails
- Check input PDB path in config
- Verify model checkpoint paths
- Ensure design residues are valid

### AF3 Out of Memory
- Reduce `modelseeds` to fewer predictions
- Run bound and apo separately
- Check sequence length (very long proteins need more memory)

### Analysis Errors
- Ensure AF3 predictions completed successfully
- Check that .cif files exist in output directories
- Verify template PDB path is correct

### Missing Dependencies
```bash
# Install BioPython
pip install biopython

# Install plotting libraries
pip install matplotlib seaborn

# Install YAML parser
pip install pyyaml
```

---

## Example Workflow

```bash
# 1. Edit configuration
nano pipeline_config.yaml

# 2. Run pipeline interactively
./run_pipeline.sh

# Or run steps individually:

# Generate 10,000 sequences (~10-30 min)
./01_run_ligandmpnn.sh

# Filter to top 1,000 and prepare AF3 inputs (~1 min)
python 02_filter_and_prepare_af3.py

# Analyze LigandMPNN sequences immediately
python 04a_analyze_ligandmpnn.py

# Run AlphaFold3 (bound) - can take many hours
./03_run_alphafold3.sh bound

# Run AlphaFold3 (apo) - while bound is running or after
./03_run_alphafold3.sh apo

# Analyze structures as they complete
python analyze_af3_structures.py \
    --af3-dir ./pipeline_outputs/03_alphafold3_predictions_bound \
    --template ./inputs/your_protein.pdb \
    --output-dir ./pipeline_outputs/04b_structure_analysis_bound

# Compare bound vs apo (once both are complete)
python compare_ligand_states.py \
    --bound-dir ./pipeline_outputs/03_alphafold3_predictions_bound \
    --apo-dir ./pipeline_outputs/03_alphafold3_predictions_apo \
    --output-dir ./pipeline_outputs/04c_ligand_state_comparison
```

---

## Modifying the Pipeline

### Change Number of Sequences
Edit `pipeline_config.yaml`:
```yaml
ligandmpnn:
  num_sequences: 5000  # Fewer sequences for faster testing
```

### Adjust Filtering Threshold
```yaml
filtering:
  top_percent: 5  # More stringent (top 5%)
  top_percent: 20 # Less stringent (top 20%)
```

### Skip Apo Predictions
```yaml
alphafold3:
  create_ligand_free: false  # Only create bound inputs
```

### Change Chromophore Region
```yaml
analysis:
  chromophore_start: 65
  chromophore_end: 67
```

---

## Citations

If you use this pipeline, please cite:

- **LigandMPNN**: Dauparas et al., Science (2022)
- **AlphaFold3**: Abramson et al., Nature (2024)
- **Analysis tools**: See individual script headers

---

## Support

For issues or questions:
1. Check this README
2. Review individual script help: `python script.py --help`
3. Check output logs in pipeline_outputs/

---

## License

See LICENSE file in this directory.

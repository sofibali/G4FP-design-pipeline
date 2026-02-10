# Simplified Protein Design Pipeline - File Summary

## Created Pipeline Files

### Configuration
- **pipeline_config.yaml** - Central configuration file for all parameters

### Execution Scripts
1. **01_run_ligandmpnn.sh** - Run LigandMPNN sequence design
2. **02_filter_and_prepare_af3.py** - Filter sequences and create AF3 inputs
3. **03_run_alphafold3.sh** - Execute AlphaFold3 predictions
4. **04a_analyze_ligandmpnn.py** - Analyze sequence diversity

### Orchestration
- **run_pipeline.sh** - Master script to run full pipeline interactively

### Documentation
- **PIPELINE_README.md** - Comprehensive usage guide
- **QUICKSTART.sh** - Quick reference commands

### Existing Analysis Scripts (Reused)
- **analyze_af3_structures.py** - Structure analysis (Step 4b)
- **compare_ligand_states.py** - Bound vs apo comparison (Step 4c)

---

## Pipeline Flow Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                         INPUT                                    │
│                                                                   │
│  • PDB file with protein + ligands                               │
│  • Design residues specification                                 │
│  • Configuration (pipeline_config.yaml)                          │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                    STEP 1: LigandMPNN                            │
│                   01_run_ligandmpnn.sh                           │
│                                                                   │
│  Input:  PDB file                                                │
│  Output: ~10,000 sequences with confidence scores                │
│  Time:   10-30 minutes                                           │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│              STEP 1b: Pre-Filter Analysis (NEW!)                 │
│              01b_analyze_pre_filter.py                           │
│                                                                   │
│  Input:  All LigandMPNN sequences                                │
│  Output: 4 plots + summary statistics                            │
│          - Confidence distributions                              │
│          - Sequence diversity per position                       │
│          - Ligand vs overall confidence scatter                  │
│  Time:   <1 minute                                               │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│              STEP 2: Filter & Prepare AF3 Inputs                 │
│              02_filter_and_prepare_af3.py                        │
│                                                                   │
│  Input:  LigandMPNN sequences                                    │
│  Filter: Top 10% by confidence (~1,000 sequences)                │
│  Output: JSON files for AlphaFold3                               │
│          - Bound state (protein+DNA)                             │
│          - Apo state (protein only)                              │
│  Time:   <1 minute                                               │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                  STEP 3: AlphaFold3                              │
│                03_run_alphafold3.sh                              │
│                                                                   │
│  Input:  JSON files (bound and/or apo)                          │
│  Output: Predicted structures (.cif)                             │
│          - 5 models per design (modelseeds)                      │
│          - Confidence scores (pLDDT, pTM, ipTM)                  │
│  Time:   30 min - 2 hrs per design                               │
│          (HOURS to DAYS total)                                   │
└─────────────┬───────────────────────┬──────────────────────────┘
              │                       │
              ▼                       ▼
    Bound Structures         Apo Structures
    (protein+DNA)           (protein only)
              │                       │
              └───────────┬───────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────────────┐
│                   STEP 4: Analysis                               │
└─────────────────────────────────────────────────────────────────┘
              │
              ├──────────────────────────┐
              │                          │
              ▼                          ▼
┌──────────────────────────┐  ┌──────────────────────────┐
│    4a: LigandMPNN        │  │  4b: Structure Analysis  │
│    Sequence Analysis     │  │  analyze_af3_structures  │
│                          │  │                          │
│  • Confidence scores     │  │  • RMSD to template      │
│  • Sequence diversity    │  │  • pLDDT scores          │
│  • Positional entropy    │  │  • Per-residue metrics   │
│  • Template identity     │  │  • Chromophore analysis  │
│                          │  │                          │
│  Time: ~1 minute         │  │  Time: 5-10 minutes      │
└──────────────────────────┘  └────────┬─────────────────┘
                                       │
                                       ▼
                          ┌──────────────────────────┐
                          │  4c: Bound vs Apo        │
                          │  compare_ligand_states   │
                          │                          │
                          │  • Structural changes    │
                          │  • RMSD (bound-apo)      │
                          │  • pLDDT differences     │
                          │  • Conformational shifts │
                          │                          │
                          │  Time: 5-10 minutes      │
                          └──────────┬───────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────┐
│                        OUTPUTS                                   │
│                                                                   │
│  • CSV files (for Prism/GraphPad)                                │
│  • Visualization plots (PNG)                                     │
│  • Summary statistics                                            │
│  • Ranked design candidates                                      │
└─────────────────────────────────────────────────────────────────┘
```

---

## File Dependencies

```
pipeline_config.yaml (required by all scripts)
    │
    ├── 01_run_ligandmpnn.sh
    │       │
    │       └── Outputs: 01_ligandmpnn/seqs/*.fa
    │                    └── Used by ↓
    │
    ├── 02_filter_and_prepare_af3.py
    │       │
    │       ├── Reads: 01_ligandmpnn/seqs/*.fa
    │       │        input_pdb (for DNA sequence)
    │       │
    │       └── Outputs: 02_filtered_sequences.csv
    │                    02_alphafold3_inputs_bound/*.json
    │                    02_alphafold3_inputs_apo/*.json
    │                    └── Used by ↓
    │
    ├── 03_run_alphafold3.sh
    │       │
    │       ├── Reads: 02_alphafold3_inputs_bound/*.json
    │       │         02_alphafold3_inputs_apo/*.json
    │       │
    │       └── Outputs: 03_alphafold3_predictions_bound/design_*/
    │                    03_alphafold3_predictions_apo/design_*/
    │                    └── Used by ↓
    │
    ├── 04a_analyze_ligandmpnn.py
    │       │
    │       ├── Reads: 01_ligandmpnn/seqs/*.fa
    │       │
    │       └── Outputs: 04a_ligandmpnn_analysis/
    │                    ├── plots (PNG)
    │                    └── plot_data_csvs/
    │
    ├── analyze_af3_structures.py (Step 4b)
    │       │
    │       ├── Reads: 03_alphafold3_predictions_bound/ (or _apo)
    │       │         template PDB
    │       │
    │       └── Outputs: 04b_structure_analysis_bound/ (or _apo)
    │                    ├── structure_analysis_results.csv
    │                    ├── plots (PNG)
    │                    └── plot_data_csvs/
    │
    └── compare_ligand_states.py (Step 4c)
            │
            ├── Reads: 03_alphafold3_predictions_bound/
            │         03_alphafold3_predictions_apo/
            │
            └── Outputs: 04c_ligand_state_comparison/
                         ├── ligand_state_comparison_results.csv
                         ├── plots (PNG)
                         └── plot_data_csvs/
```

---

## Key Features

### Modularity
- Each step is an independent script
- Can run steps individually or together
- Easy to resume from any point

### Flexibility
- Configure everything via YAML file
- Run bound-only, apo-only, or both
- Adjust filtering thresholds easily

### Comprehensive Analysis
- Sequence diversity metrics
- Structure quality assessment  
- Conformational change analysis
- CSV exports for external tools

### User-Friendly
- Interactive master script
- Clear progress indicators
- Detailed documentation
- Quick start guide

---

## Usage Examples

### Example 1: Full Pipeline
```bash
# Edit config, then run everything
nano pipeline_config.yaml
./run_pipeline.sh
```

### Example 2: Test Run (Fast)
```bash
# Edit config for small test:
# num_sequences: 100
# modelseeds: [0]

./run_pipeline.sh
# Complete in ~1-2 hours instead of days
```

### Example 3: Stepwise Execution
```bash
# Run design
./01_run_ligandmpnn.sh

# Analyze sequences before continuing
python 04a_analyze_ligandmpnn.py

# If satisfied, continue to AF3
python 02_filter_and_prepare_af3.py
./03_run_alphafold3.sh both
```

### Example 4: Bound Only Workflow
```bash
# Edit config: create_ligand_free: false
./01_run_ligandmpnn.sh
python 02_filter_and_prepare_af3.py
./03_run_alphafold3.sh bound
python analyze_af3_structures.py \
    --af3-dir ./pipeline_outputs/03_alphafold3_predictions_bound \
    --template ./inputs/your_protein.pdb \
    --output-dir ./pipeline_outputs/04b_structure_analysis
```

---

## Time Estimates

| Step | Typical Time | Adjustable Parameters |
|------|--------------|----------------------|
| 1. LigandMPNN | 10-30 min | num_sequences, batch_size |
| 1b. Pre-Filter Analysis | <1 min | - |
| 2. Filter | <1 min | top_percent |
| 3. AlphaFold3 | 30min-2hrs per design | modelseeds (predictions per design) |
| 4a. Sequence Analysis | 1 min | - |
| 4b. Structure Analysis | 5-10 min | Number of designs |
| 4c. Comparison | 5-10 min | Number of designs |

**Total Pipeline Time**: 
- Test (100 seqs, 10 designs, 1 model): 2-3 hours
- Small (1K seqs, 100 designs, 5 models): 1-2 days
- Full (10K seqs, 1000 designs, 5 models): 1-2 weeks

---

## Next Steps

1. **Edit pipeline_config.yaml** with your parameters
2. **Review PIPELINE_README.md** for detailed instructions
3. **Run test**: Start with small num_sequences to validate
4. **Scale up**: Increase to production parameters
5. **Analyze results**: Use CSV exports in external tools

---

## Support Files

- See **PIPELINE_README.md** for comprehensive documentation
- See **QUICKSTART.sh** for copy-paste commands
- See **pipeline_config.yaml** for all configurable parameters
- Each Python script has `--help` option for detailed usage

# G4FP Design Pipeline

Computational pipeline for designing G-quadruplex-dependent fluorescent proteins (G4FPs) -- proteins that adopt a GFP fold and fluoresce **only** when bound to the DNA ligand G-quadruplex (G4), and are unstructured/non-fluorescent when unbound.

## Scientific Rationale

G4FP is based on the circularly permuted GFP (cpGFP) scaffold used in calcium sensors like GCaMP. The key idea is to couple protein folding/fluorescence to G4 DNA binding:

- **Bound state (holo):** Protein + G4 DNA + K+ ions. The G4 stabilizes the barrel fold, the chromophore (CRO, residues 197-199) matures, and the protein fluoresces.
- **Unbound state (apo):** Protein alone. Without G4 stabilization, the barrel is partially or fully unfolded and the chromophore cannot mature -- no fluorescence.

The design goal is to maximize the **stability difference** between bound and unbound states: high pLDDT/pTM when bound, low pLDDT/disordered when apo. The final target is **100-1000 sequences** with the largest bound-vs-apo structural difference for experimental testing.

> **Note:** ThermoMPNN was tested in an earlier iteration for stability prediction but showed no useful discriminating features for this system and was removed from the current pipeline.

## Pipeline Overview

```
10 template PDBs
       |
       v
[01] LigandMPNN ──> 10,000 sequences each (100K total)
       |
[01b] Pre-filter analysis (confidence distributions, entropy)
       |
[02] Filter (top-N or Pareto) + generate AF3 JSONs ──> 100 bound + 100 apo per template
       |
[02.5] Validate AF3 JSONs
       |
[03] AlphaFold3 predictions (bound + apo, 5 seeds each)
       |
[04a] LigandMPNN sequence analysis
[05]  AF3 structure analysis (pLDDT, RMSD, PAE, consensus)
[06]  Bound vs apo comparison (delta pLDDT, interface PAE, conformational change)
       |
[07] Aggregate across all 10 templates + Pareto + diversity-weighted selection
       |
[08] Export for gene synthesis (plate maps, codon optimization)
```

## Pipeline Status

| Structure | LigandMPNN | Filtered | AF3 JSONs | AF3 Predictions |
|-----------|:----------:|:--------:|:---------:|:---------------:|
| G4FP_des1_cro_mod0 | 10,000 | 100 | 200 | pending |
| G4FP_des1_cro_mod4 | 10,000 | 100 | 200 | pending |
| G4FP_r20_l1_cro_mod0 | 10,000 | 100 | 200 | pending |
| G4FP_r20_l2_cro_mod2 | 10,000 | 100 | 200 | pending |
| G4FP_r20_l2_cro_mod3 | 10,000 | 100 | 200 | pending |
| G4FP_r28_l1_cro_mod2 | 10,000 | 100 | 200 | pending |
| G4FP_r28_l2_cro_mod2 | 10,000 | 100 | 200 | pending |
| G4FP_r28_l2_cro_mod3 | 10,000 | 100 | 200 | pending |
| G4FP_r28_l3_cro_mod0 | 10,000 | 100 | 200 | pending |
| G4FP_r28_l3_cro_mod2 | 10,000 | 100 | 200 | pending |

**Total:** 100,000 designed sequences. All 10 filtered (100 each = 1,000 designs). 2,000 AF3 JSONs ready (1,000 bound + 1,000 apo). AF3 predictions pending.

## Directory Structure

```
designG4FP/
├── inputs/                          # 10 template PDB structures
│
├── pipeline_config.yaml             # Central configuration
├── pipeline_config_template.yaml    # Clean config backup
│
│   ── Core pipeline scripts ──
├── 01_run_ligandmpnn.sh             # Step 1: single-structure LigandMPNN
├── 01b_analyze_pre_filter.py        # Step 1b: pre-filter sequence analysis
├── 02_filter_and_prepare_af3.py     # Step 2: filter + AF3 JSON generation
├── 03_run_alphafold3.sh             # Step 3: AF3 (sequential)
├── 03_run_alphafold3_parallel.sh    # Step 3: AF3 (parallel)
├── 04a_analyze_ligandmpnn.py        # Step 4a: sequence diversity analysis
├── 05_analyze_af3_structures.py     # Step 5: structure quality + PAE + consensus
├── 06_compare_ligand_states.py      # Step 6: bound vs apo comparison
├── 07_aggregate_and_rank.py         # Step 7: cross-template ranking (Pareto + diversity)
├── 08_export_for_synthesis.py       # Step 8: synthesis order export
│
│   ── Batch runners ──
├── run_all.sh                       # Non-interactive full pipeline (all structures)
├── run_pipeline.sh                  # Interactive single-structure pipeline
├── batch_run_ligandmpnn.sh          # Batch LigandMPNN (all structures)
├── batch_filter_all.sh              # Batch filtering (all structures)
├── run_all_af3_parallel.sh          # Batch AF3 (all structures)
├── run_structure_analysis.sh        # Batch analysis (all structures)
│
│   ── Utilities ──
├── utils/
│   ├── validate_af3_jsons.py        # Validate AF3 JSONs before submission
│   ├── fix_af3_json_ids.py          # Fix missing id fields
│   ├── fix_af3_json_dialect.py      # Fix dialect field
│   ├── fix_af3_json_structure.py    # Fix JSON structure
│   ├── check_pipeline_status.py     # Pipeline completion checker
│   └── test_analysis_setup.py       # Dependency checker
│
│   ── Archive (historical fixes/docs) ──
├── archive/
│   ├── docs/                        # Old documentation files
│   └── scripts/                     # One-off fix scripts
│
│   ── Outputs (per template, x10) ──
├── output_G4FP_<name>/
│   ├── 01_ligandmpnn/seqs/*.fa      # 10,000 designed sequences
│   ├── 01b_pre_filter_analysis/     # Confidence plots, entropy
│   ├── 02_filtered_sequences.csv    # Top 100 sequences
│   ├── 02_alphafold3_inputs_{bound,apo}/  # AF3 input JSONs
│   ├── 03_alphafold3_predictions_{bound,apo}/  # AF3 CIF structures
│   ├── 04a_ligandmpnn_analysis/     # Sequence analysis plots
│   ├── 05_structure_analysis/       # Structure quality plots
│   └── 06_ligand_state_comparison/  # Bound vs apo comparison
│
│   ── Final outputs ──
├── 07_results/                      # Aggregated ranking across all templates
│   ├── 07_final_candidates.csv
│   ├── 07_final_candidates.fa
│   ├── 07_pareto_frontier.csv
│   └── *.png
└── 08_synthesis_order/              # Gene synthesis ordering files
    ├── 08_synthesis_sequences.csv
    ├── 08_plate_map.csv
    └── 08_sequences.fa
```

## Quick Start

### Prerequisites

```bash
conda activate unified_mpnn
pip install biopython pandas numpy matplotlib seaborn pyyaml scipy openpyxl
```

### Run everything (non-interactive, all 10 structures)

```bash
cd /home/sbali/LigandMPNN/designG4FP

# Full pipeline, 3 parallel AF3 jobs, with Pareto filtering
nohup ./run_all.sh --parallel 3 --pareto 2>&1 | tee run_all.log &

# Skip LigandMPNN (already done) and AF3 (not ready yet)
./run_all.sh --skip-ligandmpnn --skip-af3
```

### Run individual steps

```bash
# Step 1: Design sequences (~10-30 min per structure)
./01_run_ligandmpnn.sh pipeline_config.yaml

# Step 2: Filter (top-N or Pareto) and generate AF3 inputs
python 02_filter_and_prepare_af3.py                              # from config
python 02_filter_and_prepare_af3.py --output-dir DIR --input-pdb PDB  # direct
python 02_filter_and_prepare_af3.py --pareto --top-n 200         # Pareto mode

# Step 2.5: Validate AF3 JSONs before submitting
python utils/validate_af3_jsons.py          # check all
python utils/validate_af3_jsons.py --fix    # auto-fix missing ids

# Step 3: AlphaFold3 (parallel)
./03_run_alphafold3_parallel.sh pipeline_config.yaml both 3

# Steps 5-6: Analysis (after AF3 completes)
python 05_analyze_af3_structures.py --output-dir output_G4FP_des1_cro_mod0 --template inputs/G4FP_des1_cro_mod0.pdb --state both
python 06_compare_ligand_states.py --output-dir output_G4FP_des1_cro_mod0 --template inputs/G4FP_des1_cro_mod0.pdb

# Step 7: Aggregate and rank across all templates
python 07_aggregate_and_rank.py --n-select 500

# Step 8: Export for synthesis
python 08_export_for_synthesis.py --organism ecoli --prefix G4FP_v1
```

## Filtering Modes

### Top-N (default)
Ranks all sequences by `(overall_confidence + ligand_confidence) / 2` and takes the top N. Simple and fast.

### Pareto filtering (`--pareto`)
Finds sequences on the Pareto frontier of `overall_confidence` vs `ligand_confidence` -- designs where you cannot improve one metric without worsening the other. This captures sequences with exceptional ligand confidence that might be missed by simple averaging, and vice versa.

### Diversity-weighted selection (Step 7)
After scoring, sequences are clustered by Hamming distance. The top-scoring sequence from each cluster is selected first, then remaining slots are filled round-robin across clusters. This prevents the final set from being dominated by near-identical sequences.

## Key Metrics

### LigandMPNN (per sequence)
- `overall_confidence` -- model confidence over redesigned positions (0-1)
- `ligand_confidence` -- confidence for residues near the DNA (0-1)

### AlphaFold3 (per structure)
- `pLDDT` -- per-residue confidence (0-100, >70 = confident fold)
- `pTM` -- predicted TM-score (0-1, >0.5 = plausible fold)
- `ipTM` -- interface pTM for protein-DNA quality
- `PAE` -- predicted aligned error matrix (lower = better)
- `interface_pae` -- PAE specifically at the protein-DNA interface
- `consensus_score` -- agreement across AF3 model seeds (1 = all agree)

### Final selection (Step 7)
- `fitness_score` -- weighted composite: delta_pLDDT (35%) + bound_pTM (20%) + apo_disorder (20%) + chromophore_diff (15%) + RMSD (10%)
- Hard filters: bound_pLDDT > 70, bound_pTM > 0.5, apo_pLDDT < 60
- Pareto frontier across 5 objectives
- Diversity-weighted to cover sequence space

## Configuration

All parameters in `pipeline_config.yaml`:

```yaml
input_pdb: ./inputs/G4FP_des1_cro_mod0.pdb
output_dir: ./output_G4FP_des1_cro_mod0/

ligandmpnn:
  design_residues: "A15-60 A240-285"   # 92 positions
  num_sequences: 10000
  cutoff_for_score: 8.0                # ligand context distance

filtering:
  top_n: 100
  top_percent: 1
  sort_by: overall_confidence

alphafold3:
  create_ligand_free: true
  modelseeds: [0, 1, 2, 3, 4]
  db_dir: /mnt/alphafold3
  model_dir: /mnt/alphafold3

analysis:
  chromophore_start: 197
  chromophore_end: 199
```

## Dependencies

- **LigandMPNN** -- sequence design with ligand context (this repo)
- **AlphaFold3** -- structure prediction (via SBGrid at `/programs/sbgrid.shrc`)
- **BioPython** -- PDB/CIF parsing, RMSD calculations
- **pandas/numpy/scipy** -- data analysis, clustering, Pareto optimization
- **matplotlib/seaborn** -- visualization
- **PyYAML** -- configuration

## References

- Dauparas et al. "Robust deep learning-based protein sequence design using ProteinMPNN" *Science* (2022)
- Abramson et al. "Accurate structure prediction of biomolecular interactions with AlphaFold 3" *Nature* (2024)

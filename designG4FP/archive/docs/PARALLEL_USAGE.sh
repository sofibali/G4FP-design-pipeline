#!/bin/bash
# Quick reference for running AlphaFold3 predictions in parallel

cat << 'EOF'
============================================================================
AlphaFold3 Parallel Prediction - Quick Reference
============================================================================

BASIC USAGE:
  ./03_run_alphafold3_parallel.sh [config] [mode] [num_parallel]

PARAMETERS:
  config       - YAML config file (default: pipeline_config.yaml)
  mode         - bound|apo|both (default: both)
  num_parallel - Number of parallel jobs (default: 2)

EXAMPLES:

  # Run with default settings (2 parallel jobs, both modes)
  ./03_run_alphafold3_parallel.sh

  # Run only bound predictions with 4 parallel jobs
  ./03_run_alphafold3_parallel.sh pipeline_config.yaml bound 4

  # Run apo predictions with conservative parallelism
  ./03_run_alphafold3_parallel.sh pipeline_config.yaml apo 1

GPU MEMORY GUIDANCE:
  • 12GB VRAM:  Use 1-2 parallel jobs
  • 24GB VRAM:  Use 2-4 parallel jobs  ✓ Safe default
  • 48GB VRAM:  Use 4-6 parallel jobs
  • 80GB A100:  Use 6-10 parallel jobs

PERFORMANCE:
  Sequential:  100 structures × 30min = ~50 hours
  2 parallel:  ~25 hours
  4 parallel:  ~12-15 hours (with 24GB+ GPU)

MONITORING:
  # Watch GPU usage in real-time
  watch -n 1 nvidia-smi

  # Check progress
  ls output_G4FP*/03_alphafold3_predictions_bound/*/model_cif/model_0.cif | wc -l

  # View recent failures
  tail -20 output_G4FP*/03_alphafold3_predictions_bound/*.log | grep -B 5 'Error'

TROUBLESHOOTING:
  • All jobs fail with JSON errors:
    → Run: python fix_af3_json_dialect.py <input_json_dir>

  • GPU out of memory:
    → Reduce num_parallel to 1 or 2

  • Jobs not running in parallel:
    → Install GNU Parallel: conda install -c conda-forge parallel

  • Want to resume after interrupt:
    → Just re-run the script - it skips completed structures

INSTALL GNU PARALLEL (recommended for better performance):
  conda install -c conda-forge parallel
  # or
  sudo apt-get install parallel

============================================================================
EOF

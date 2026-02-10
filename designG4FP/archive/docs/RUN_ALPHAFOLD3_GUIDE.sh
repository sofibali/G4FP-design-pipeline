#!/bin/bash
# LigandMPNN G4FP - AlphaFold3 Batch Runner
# Quick reference for running all ready AF3 predictions

cd /home/sbali/LigandMPNN/designG4FP

cat << 'EOF'
================================================================================
  LigandMPNN G4FP - AlphaFold3 Predictions Ready Status
================================================================================

📊 READY FOR ALPHAFOLD3:

  Structure                      | Bound | Apo  | Total | Status
  -------------------------------|-------|------|-------|--------
  output_G4FP_des1_cro_mod0      | 100   | 100  | 200   | ✅ Ready
  output_G4FP_des1_cro_mod4      | 100   | 100  | 200   | ✅ Ready  
  output_G4FP_r20_l1_cro_mod0    | 100   | 100  | 200   | ✅ Ready
  output_G4FP_r20_l2_cro_mod2    | 100   | 100  | 200   | ✅ Ready
  -------------------------------|-------|------|-------|--------
  TOTAL                          | 400   | 400  | 800   | ✅ Ready

✅ JSON FORMAT: All files have correct 'dialect' and 'version' fields
✅ PARALLEL SCRIPT: 03_run_alphafold3_parallel.sh is ready and executable

================================================================================
  QUICK START OPTIONS
================================================================================

OPTION 1: Run ALL structures in sequence (recommended for monitoring)
────────────────────────────────────────────────────────────────────────
for output_dir in output_G4FP_des1_cro_mod0 output_G4FP_des1_cro_mod4 \
                  output_G4FP_r20_l1_cro_mod0 output_G4FP_r20_l2_cro_mod2; do
    echo "Processing $output_dir..."
    cd "$output_dir"
    
    # Run with 3 parallel jobs on GPU 0
    ../03_run_alphafold3_parallel.sh pipeline_config_temp.yaml both 3
    
    cd ..
done

OPTION 2: Run ONE structure at a time (easier to monitor)
────────────────────────────────────────────────────────────────────────
# Start with the first one
cd output_G4FP_des1_cro_mod0

# Run both bound and apo with 3 parallel jobs
../03_run_alphafold3_parallel.sh pipeline_config_temp.yaml both 3

# Or run bound first, then apo
../03_run_alphafold3_parallel.sh pipeline_config_temp.yaml bound 3
../03_run_alphafold3_parallel.sh pipeline_config_temp.yaml apo 3

OPTION 3: Background execution with logging
────────────────────────────────────────────────────────────────────────
cd output_G4FP_des1_cro_mod0

nohup bash ../03_run_alphafold3_parallel.sh pipeline_config_temp.yaml both 3 \
    > af3_parallel_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# Monitor progress
tail -f af3_parallel_*.log

================================================================================
  PARALLEL JOB TUNING
================================================================================

Adjust the number based on your GPU memory:

  GPU Memory | Recommended Parallel Jobs | Command
  -----------|---------------------------|----------------------------------
  12 GB      | 1-2 jobs                  | ...parallel.sh config both 2
  16-24 GB   | 2-4 jobs (recommended)    | ...parallel.sh config both 3
  32-48 GB   | 4-6 jobs                  | ...parallel.sh config both 5
  80 GB A100 | 6-10 jobs                 | ...parallel.sh config both 8

Current GPU setup: Check with 'nvidia-smi'

⚠️  If you get OOM errors, reduce the number of parallel jobs

================================================================================
  MONITORING PROGRESS
================================================================================

# Watch GPU utilization in real-time
watch -n 1 nvidia-smi

# Count completed predictions
find output_G4FP_*/03_alphafold3_predictions_* -name "model_0.cif" | wc -l

# Check for failures
grep -r "Error\|Failed" output_G4FP_*/03_alphafold3_predictions_*/*.log | head -20

# Monitor a specific structure
watch -n 10 'ls output_G4FP_des1_cro_mod0/03_alphafold3_predictions_bound/ | wc -l'

# Check overall status anytime
python check_pipeline_status.py

================================================================================
  ESTIMATED COMPLETION TIME
================================================================================

Per structure (100 sequences × 2 states = 200 predictions):
  - Sequential: ~200-300 minutes (~4-5 hours)
  - 2 parallel:  ~100-150 minutes (~2-2.5 hours)
  - 3 parallel:  ~70-100 minutes (~1.5-2 hours)
  - 4 parallel:  ~50-75 minutes (~1-1.5 hours)

For all 4 structures (800 total predictions):
  - 2 parallel:  ~8-10 hours
  - 3 parallel:  ~6-8 hours
  - 4 parallel:  ~4-6 hours

💡 Recommendation: Start with 2-3 parallel jobs, monitor GPU memory with 
   nvidia-smi, and increase if you have headroom.

================================================================================
  TROUBLESHOOTING
================================================================================

Problem: GPU Out of Memory
Solution: Reduce NUM_PARALLEL parameter

Problem: Jobs not starting
Solution: Check GPU is available with 'nvidia-smi'
          Ensure CUDA is set up: echo $CUDA_VISIBLE_DEVICES

Problem: Very slow progress
Solution: Check if GNU parallel is installed for better job management
          conda install -c conda-forge parallel

Problem: Predictions failing
Solution: Check logs in output_*/03_alphafold3_predictions_*/*.log
          Common issues: invalid sequence, JSON format, database access

================================================================================
  NEXT STEPS AFTER COMPLETION
================================================================================

1. Verify all predictions completed:
   python check_pipeline_status.py

2. Run analysis on completed structures:
   bash 04a_analyze_ligandmpnn.py

3. Compare bound vs apo predictions:
   bash 04b_compare_states.py

4. Generate summary reports:
   bash 04c_summarize_results.py

================================================================================

EOF

echo ""
echo "💡 To execute, copy commands from above or run:"
echo ""
echo "   cd output_G4FP_des1_cro_mod0"
echo "   ../03_run_alphafold3_parallel.sh pipeline_config_temp.yaml both 3"
echo ""

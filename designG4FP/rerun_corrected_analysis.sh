#!/bin/bash
# Re-run 05/06/07 with the corrected chromophore localisation and per-residue
# pLDDT (see chromophore_utils.py).  Round-1 outputs are preserved in git at
# d5e056e49 and snapshotted to round1_snapshot/ before anything is overwritten.
set -u
cd /home/sbali/LigandMPNN/designG4FP || exit 1
PY=/home/sbali/LigandMPNN/.venv/bin/python
LOG_DIR=analysis_logs/rerun_$(date +%Y%m%d_%H%M%S)
mkdir -p "$LOG_DIR"
echo "logs -> $LOG_DIR"

# ---- 0. snapshot round-1 tables for side-by-side comparison --------------
SNAP=round1_snapshot
if [ ! -d "$SNAP" ]; then
    mkdir -p "$SNAP"
    cp -p 07_results/07_final_candidates.csv  "$SNAP/" 2>/dev/null
    cp -p 07_results/07_pareto_frontier.csv   "$SNAP/" 2>/dev/null
    cp -p 07_results/07_top50_global.csv      "$SNAP/" 2>/dev/null
    cp -p 07_results/07_top10_per_template.csv "$SNAP/" 2>/dev/null
    cp -p template_af3_metrics.csv            "$SNAP/" 2>/dev/null
    cp -p analysis_output/06_ligand_state_comparison_results.csv "$SNAP/" 2>/dev/null
    echo "snapshotted round-1 tables -> $SNAP/"
else
    echo "$SNAP/ already exists, leaving it alone"
fi

# ---- 1. template reference metrics (feeds the 10% rel_apo term) ----------
echo "=== [1/4] analyze_template_af3.py ==="
$PY analyze_template_af3.py > "$LOG_DIR/01_template.log" 2>&1
echo "    exit $?  ($(wc -l < "$LOG_DIR/01_template.log") lines)"

# ---- 2+3. per-template structure analysis and holo/apo comparison --------
n=0
for d in output_G4FP_*; do
    n=$((n+1))
    echo "=== [2/4] 05 on $d  ($n/10) ==="
    $PY 05_analyze_af3_structures.py --output-dir "$d" --force \
        > "$LOG_DIR/05_${d}.log" 2>&1
    echo "    exit $?"
    echo "=== [3/4] 06 on $d  ($n/10) ==="
    $PY 06_compare_ligand_states.py --output-dir "$d" --force \
        > "$LOG_DIR/06_${d}.log" 2>&1
    echo "    exit $?"
    grep -h "Chromophore:" "$LOG_DIR/05_${d}.log" | head -1
done

# ---- 4. aggregate and rank ----------------------------------------------
echo "=== [4/4] 07_aggregate_and_rank.py ==="
$PY 07_aggregate_and_rank.py > "$LOG_DIR/07_aggregate.log" 2>&1
echo "    exit $?"

echo "DONE -- logs in $LOG_DIR"

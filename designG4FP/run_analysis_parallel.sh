#!/bin/bash
# Run 05, 06, 07 analysis in parallel across templates.
#
# Steps 05 and 06 are per-template and can run fully in parallel.
# Step 07 aggregates all results and runs after 05+06 complete.
#
# Usage:
#   bash run_analysis_parallel.sh           # all templates, max parallelism
#   bash run_analysis_parallel.sh -j 5      # limit to 5 parallel jobs
#   bash run_analysis_parallel.sh --dry-run # preview only

set -euo pipefail

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

MAX_JOBS=10
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j) MAX_JOBS="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

TEMPLATES=($(ls -d output_G4FP_* 2>/dev/null | sed 's/output_//'))

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Parallel Analysis Pipeline${NC}"
echo -e "${GREEN}========================================${NC}"
echo "  Templates: ${#TEMPLATES[@]}"
echo "  Max parallel jobs: $MAX_JOBS"
echo ""

if $DRY_RUN; then
    echo -e "${YELLOW}DRY RUN — would run:${NC}"
    for tpl in "${TEMPLATES[@]}"; do
        echo "  python3 05_analyze_af3_structures.py --template $tpl --force"
        echo "  python3 06_compare_ligand_states.py --template $tpl --force"
    done
    echo "  python3 07_aggregate_and_rank.py --force"
    exit 0
fi

LOG_DIR="$SCRIPT_DIR/analysis_logs"
mkdir -p "$LOG_DIR"

# ── Step 05: Structure analysis (parallel) ──
echo -e "${BLUE}Step 05: Running structure analysis (${#TEMPLATES[@]} templates parallel)...${NC}"

PIDS_05=()
for tpl in "${TEMPLATES[@]}"; do
    # Wait if at max jobs
    while [ ${#PIDS_05[@]} -ge $MAX_JOBS ]; do
        new_pids=()
        for pid in "${PIDS_05[@]}"; do
            if kill -0 "$pid" 2>/dev/null; then
                new_pids+=("$pid")
            fi
        done
        PIDS_05=("${new_pids[@]}")
        if [ ${#PIDS_05[@]} -ge $MAX_JOBS ]; then
            sleep 2
        fi
    done

    echo "  Starting 05: $tpl"
    python3 05_analyze_af3_structures.py --template "$tpl" --force \
        > "$LOG_DIR/05_${tpl}.log" 2>&1 &
    PIDS_05+=($!)
done

# Wait for all 05 jobs
FAIL_05=0
for pid in "${PIDS_05[@]}"; do
    wait "$pid" || ((FAIL_05++))
done

if [ $FAIL_05 -gt 0 ]; then
    echo -e "${RED}  Step 05: $FAIL_05 template(s) failed${NC}"
    echo "  Check logs in $LOG_DIR/05_*.log"
else
    echo -e "${GREEN}  Step 05: All ${#TEMPLATES[@]} templates complete${NC}"
fi

# ── Step 05b: Re-run 05 combined (merge all per-template results + plots) ──
echo ""
echo -e "${BLUE}Step 05b: Generating combined 05 plots (all designs)...${NC}"

python3 05_analyze_af3_structures.py --force \
    > "$LOG_DIR/05_combined.log" 2>&1

if [ $? -eq 0 ]; then
    echo -e "${GREEN}  Step 05b: Complete${NC}"
else
    echo -e "${RED}  Step 05b: Failed — check $LOG_DIR/05_combined.log${NC}"
fi

# ── Step 06: Ligand state comparison (parallel) ──
echo ""
echo -e "${BLUE}Step 06: Running ligand state comparison (${#TEMPLATES[@]} templates parallel)...${NC}"

PIDS_06=()
for tpl in "${TEMPLATES[@]}"; do
    while [ ${#PIDS_06[@]} -ge $MAX_JOBS ]; do
        new_pids=()
        for pid in "${PIDS_06[@]}"; do
            if kill -0 "$pid" 2>/dev/null; then
                new_pids+=("$pid")
            fi
        done
        PIDS_06=("${new_pids[@]}")
        if [ ${#PIDS_06[@]} -ge $MAX_JOBS ]; then
            sleep 2
        fi
    done

    echo "  Starting 06: $tpl"
    python3 06_compare_ligand_states.py --template "$tpl" --force \
        > "$LOG_DIR/06_${tpl}.log" 2>&1 &
    PIDS_06+=($!)
done

FAIL_06=0
for pid in "${PIDS_06[@]}"; do
    wait "$pid" || ((FAIL_06++))
done

if [ $FAIL_06 -gt 0 ]; then
    echo -e "${RED}  Step 06: $FAIL_06 template(s) failed${NC}"
    echo "  Check logs in $LOG_DIR/06_*.log"
else
    echo -e "${GREEN}  Step 06: All ${#TEMPLATES[@]} templates complete${NC}"
fi

# ── Step 06b: Re-run 06 combined (merge all per-template results + plots) ──
echo ""
echo -e "${BLUE}Step 06b: Generating combined 06 plots (all 1000 designs)...${NC}"

python3 06_compare_ligand_states.py --force \
    > "$LOG_DIR/06_combined.log" 2>&1

if [ $? -eq 0 ]; then
    echo -e "${GREEN}  Step 06b: Complete${NC}"
else
    echo -e "${RED}  Step 06b: Failed — check $LOG_DIR/06_combined.log${NC}"
fi

# ── Step 07: Aggregate and rank (single job, needs all results) ──
echo ""
echo -e "${BLUE}Step 07: Aggregating and ranking...${NC}"

python3 07_aggregate_and_rank.py \
    > "$LOG_DIR/07_aggregate.log" 2>&1

if [ $? -eq 0 ]; then
    echo -e "${GREEN}  Step 07: Complete${NC}"
else
    echo -e "${RED}  Step 07: Failed — check $LOG_DIR/07_aggregate.log${NC}"
fi

# ── Summary ──
echo ""
echo -e "${GREEN}========================================${NC}"
TOTAL_FAIL=$((FAIL_05 + FAIL_06))
if [ $TOTAL_FAIL -eq 0 ]; then
    echo -e "${GREEN}Analysis pipeline complete!${NC}"
else
    echo -e "${RED}Analysis finished with $TOTAL_FAIL failure(s)${NC}"
    echo "  Check logs: ls $LOG_DIR/"
fi
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Logs: $LOG_DIR/"
echo "Results: 07_results/"

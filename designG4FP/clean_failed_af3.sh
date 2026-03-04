#!/bin/bash
# Clean up failed AF3 prediction outputs (empty dirs, failed logs)
# Run this before re-running AF3 to ensure a clean start.
#
# Usage:
#   ./clean_failed_af3.sh              # show what would be deleted
#   ./clean_failed_af3.sh --delete     # actually delete

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DELETE=0
if [ "$1" = "--delete" ]; then
    DELETE=1
fi

echo "Scanning for failed AF3 outputs..."
echo ""

TOTAL_EMPTY_DIRS=0
TOTAL_LOGS=0
TOTAL_BATCH=0

# Clean empty prediction subdirectories (MSA ran but produced nothing)
for pred_dir in "$SCRIPT_DIR"/output_G4FP_*/03_alphafold3_predictions_*/; do
    [ -d "$pred_dir" ] || continue

    for design_dir in "$pred_dir"/design_*/; do
        [ -d "$design_dir" ] || continue

        # Check if dir is empty or has no useful output
        n_files=$(find "$design_dir" -type f 2>/dev/null | wc -l)
        has_cif=$(find "$design_dir" -name '*.cif' -type f 2>/dev/null | wc -l)

        if [ "$n_files" -eq 0 ] || [ "$has_cif" -eq 0 ]; then
            if [ $DELETE -eq 1 ]; then
                rm -rf "$design_dir"
            else
                echo "  EMPTY: $design_dir"
            fi
            ((TOTAL_EMPTY_DIRS++))
        fi
    done

    # Clean failed log files
    for log_file in "$pred_dir"/*.log; do
        [ -f "$log_file" ] || continue
        if [ $DELETE -eq 1 ]; then
            rm -f "$log_file"
        else
            echo "  LOG:   $log_file"
        fi
        ((TOTAL_LOGS++))
    done
done

# Clean batch temp directories
for batch_dir in "$SCRIPT_DIR"/.af3_batch_gpu* "$SCRIPT_DIR"/.inference_output_gpu*; do
    if [ -d "$batch_dir" ]; then
        if [ $DELETE -eq 1 ]; then
            rm -rf "$batch_dir"
        else
            echo "  BATCH: $batch_dir"
        fi
        ((TOTAL_BATCH++))
    fi
done

echo ""
echo "Summary: $TOTAL_EMPTY_DIRS empty/failed dirs, $TOTAL_LOGS log files, $TOTAL_BATCH batch dirs"

if [ $DELETE -eq 0 ]; then
    echo ""
    echo "Run with --delete to remove these files."
else
    echo "Cleaned up."
fi

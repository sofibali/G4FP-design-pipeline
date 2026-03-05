#!/usr/bin/env python3
"""Sort AF3 inference outputs from hidden GPU dirs back to per-template directories.

AF3 used the 'name' field inside each JSON (e.g. 'design_0_bound') for output dirs,
so results from all 10 templates collided into timestamped dirs. This script matches
each result to its template by comparing sequences, then copies the results to the
correct 03_alphafold3_predictions_{bound,apo}/ directory.
"""

import json
import glob
import os
import shutil
import sys
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)

GPU_DIRS = [".inference_output_gpu0", ".inference_output_gpu1"]

# Discover templates
TEMPLATES = sorted(
    os.path.basename(d.rstrip("/"))
    for d in glob.glob(os.path.join(SCRIPT_DIR, "output_G4FP_*/"))
    if os.path.isdir(os.path.join(d, "02_alphafold3_inputs_bound"))
)
print(f"Found {len(TEMPLATES)} templates: {', '.join(t.replace('output_G4FP_','') for t in TEMPLATES)}")


def get_sequence_from_json(path):
    """Extract protein sequence from an AF3 JSON or _data.json."""
    try:
        with open(path) as f:
            j = json.load(f)
        for s in j.get("sequences", []):
            if "protein" in s:
                return s["protein"]["sequence"]
    except Exception:
        return None
    return None


def build_template_seq_map():
    """Build mapping: (design_num, state) -> {sequence: template_name}."""
    seq_map = {}  # (design_num, state, sequence) -> template
    for template in TEMPLATES:
        for state in ["bound", "apo"]:
            input_dir = f"{template}/02_alphafold3_inputs_{state}"
            if not os.path.isdir(input_dir):
                continue
            for json_file in sorted(glob.glob(f"{input_dir}/design_*_{state}.json")):
                basename = os.path.basename(json_file).replace(".json", "")
                # Extract design number: design_0000_bound -> 0
                parts = basename.split("_")
                design_num = int(parts[1])
                seq = get_sequence_from_json(json_file)
                if seq:
                    seq_map[(design_num, state, seq)] = template
    return seq_map


def identify_result_dir(result_path, seq_map):
    """Identify which template a result directory belongs to."""
    dirname = os.path.basename(result_path)

    # Parse design number and state from dirname like "design_0_bound" or "design_0_bound_20260227_105952"
    # Also handle "design_0_apo", "design_0_apo_20260227_..."
    parts = dirname.split("_")
    if len(parts) < 3 or parts[0] != "design":
        return None, None, None

    design_num = int(parts[1])
    state = parts[2]  # "bound" or "apo"
    if state not in ("bound", "apo"):
        return None, None, None

    # Get sequence from _data.json in this result dir
    data_files = glob.glob(os.path.join(result_path, "*_data.json"))
    if not data_files:
        return None, None, None

    seq = get_sequence_from_json(data_files[0])
    if not seq:
        return None, None, None

    template = seq_map.get((design_num, state, seq))
    return template, design_num, state


def main():
    dry_run = "--dry-run" in sys.argv
    if dry_run:
        print("DRY RUN - no files will be moved\n")

    print("Building template sequence map...")
    seq_map = build_template_seq_map()
    print(f"  Indexed {len(seq_map)} (design, state, sequence) entries\n")

    # Track results
    stats = defaultdict(lambda: {"bound": 0, "apo": 0})
    errors = []
    skipped = 0
    moved = 0

    for gpu_dir in GPU_DIRS:
        if not os.path.isdir(gpu_dir):
            print(f"  {gpu_dir}: not found, skipping")
            continue

        result_dirs = sorted(glob.glob(f"{gpu_dir}/design_*/"))
        # Also catch timestamped dirs
        result_dirs += sorted(glob.glob(f"{gpu_dir}/design_*_20*/"))
        # Deduplicate
        result_dirs = sorted(set(d.rstrip("/") for d in result_dirs))

        print(f"Processing {gpu_dir}: {len(result_dirs)} result directories")

        for result_path in result_dirs:
            template, design_num, state = identify_result_dir(result_path, seq_map)

            if template is None:
                errors.append(result_path)
                continue

            # Target: template/03_alphafold3_predictions_{state}/design_{NNNN}_{state}/design_{N}_{state}/
            padded_name = f"design_{design_num:04d}_{state}"
            inner_name = f"design_{design_num}_{state}"
            target_dir = os.path.join(
                template,
                f"03_alphafold3_predictions_{state}",
                padded_name,
                inner_name,
            )

            # Check if target already has model files
            existing_cifs = glob.glob(os.path.join(target_dir, "*_model.cif"))
            existing_seeds = glob.glob(os.path.join(target_dir, "seed-*"))
            if existing_cifs or existing_seeds:
                skipped += 1
                continue

            if dry_run:
                print(f"  WOULD COPY: {os.path.basename(result_path)} -> {target_dir}")
            else:
                os.makedirs(target_dir, exist_ok=True)
                # Copy all contents from result dir to target
                for item in os.listdir(result_path):
                    src = os.path.join(result_path, item)
                    dst = os.path.join(target_dir, item)
                    if os.path.isdir(src):
                        if os.path.exists(dst):
                            # Merge into existing dir
                            for sub in os.listdir(src):
                                s2 = os.path.join(src, sub)
                                d2 = os.path.join(dst, sub)
                                if not os.path.exists(d2):
                                    if os.path.isdir(s2):
                                        shutil.copytree(s2, d2)
                                    else:
                                        shutil.copy2(s2, d2)
                        else:
                            shutil.copytree(src, dst)
                    else:
                        if not os.path.exists(dst):
                            shutil.copy2(src, dst)

            stats[template][state] += 1
            moved += 1

    # Summary
    print(f"\n{'=' * 70}")
    print(f"SORT COMPLETE {'(DRY RUN)' if dry_run else ''}")
    print(f"{'=' * 70}")
    print(f"  Moved:   {moved}")
    print(f"  Skipped: {skipped} (already had results)")
    print(f"  Errors:  {len(errors)} (could not match)")
    print()

    for template in TEMPLATES:
        s = stats[template]
        total = s["bound"] + s["apo"]
        if total > 0:
            print(f"  {template}: {s['bound']} bound + {s['apo']} apo = {total}")

    if errors:
        print(f"\nUnmatched directories ({len(errors)}):")
        for e in errors[:10]:
            print(f"  {e}")
        if len(errors) > 10:
            print(f"  ... and {len(errors) - 10} more")

    # Verify final counts
    print(f"\n{'=' * 70}")
    print("VERIFICATION: model.cif counts per template")
    print(f"{'=' * 70}")
    for template in TEMPLATES:
        for state in ["bound", "apo"]:
            pred_dir = f"{template}/03_alphafold3_predictions_{state}"
            cifs = glob.glob(f"{pred_dir}/**/*_model.cif", recursive=True)
            seeds = glob.glob(f"{pred_dir}/**/seed-*/model.cif", recursive=True)
            print(f"  {template} {state:5s}: {len(cifs)} top + {len(seeds)} seeds")


if __name__ == "__main__":
    main()

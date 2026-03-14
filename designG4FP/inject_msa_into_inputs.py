#!/usr/bin/env python3
"""
Inject MSA data from existing _data.json files into AF3 input JSONs.

AF3's --norun_data_pipeline doesn't read from existing _data.json files.
Instead, it expects the input JSON to already contain unpairedMsa/pairedMsa
fields. This script copies MSA data from the _data.json (produced by a
previous AF3 data pipeline run + chromophore fix) into the input JSONs.

After running this, AF3 can be called with --norun_data_pipeline and will
see the MSA data directly in the input, avoiding the "missing unpaired MSA"
error.

Usage:
    python3 inject_msa_into_inputs.py                    # all templates
    python3 inject_msa_into_inputs.py --template G4FP_des1_cro_mod0
    python3 inject_msa_into_inputs.py --dry-run           # preview only
    python3 inject_msa_into_inputs.py --test-one          # one design only
"""

import argparse
import json
import sys
from pathlib import Path


def find_data_json(pred_dir: Path, name: str) -> Path | None:
    """Find _data.json for a given design name in the predictions dir.

    Checks multiple directory layouts:
      - padded:      design_0000_bound/design_0_bound/<name>_data.json
      - non-padded:  design_0_bound/<name>_data.json
      - timestamped: design_0_bound_YYYYMMDD_HHMMSS/<name>_data.json
    """
    # Extract seq_id and state from name like "design_0_bound"
    # Try padded first (most reliable - fix_chromophore fixed these)
    for d in sorted(pred_dir.iterdir()):
        if not d.is_dir():
            continue
        # Check for _data.json directly in this dir
        candidate = d / f"{name}_data.json"
        if candidate.exists():
            return candidate
        # Check one level deeper (padded layout)
        for sub in d.iterdir():
            if sub.is_dir():
                candidate = sub / f"{name}_data.json"
                if candidate.exists():
                    return candidate
    return None


def inject_msa(input_json_path: Path, data_json_path: Path, dry_run: bool) -> bool:
    """Copy MSA fields from _data.json into input JSON."""
    with open(data_json_path) as f:
        data = json.load(f)

    with open(input_json_path) as f:
        inp = json.load(f)

    # Find protein sequence in both
    data_protein = None
    for seq in data['sequences']:
        if 'protein' in seq:
            data_protein = seq['protein']
            break

    if data_protein is None:
        print(f"  WARNING: No protein in {data_json_path}")
        return False

    # Check MSA exists
    if not data_protein.get('unpairedMsa'):
        print(f"  WARNING: No unpairedMsa in {data_json_path}")
        return False

    # Find protein in input JSON
    inp_protein = None
    for seq in inp['sequences']:
        if 'protein' in seq:
            inp_protein = seq['protein']
            break

    if inp_protein is None:
        print(f"  WARNING: No protein in {input_json_path}")
        return False

    # Verify sequences match
    if inp_protein['sequence'] != data_protein['sequence']:
        print(f"  WARNING: Sequence mismatch for {input_json_path.name}")
        print(f"    Input:  {inp_protein['sequence'][:40]}... ({len(inp_protein['sequence'])})")
        print(f"    Data:   {data_protein['sequence'][:40]}... ({len(data_protein['sequence'])})")
        return False

    # Inject MSA fields
    fields_injected = []
    for field in ['unpairedMsa', 'pairedMsa', 'templates']:
        if field in data_protein and data_protein[field]:
            inp_protein[field] = data_protein[field]
            fields_injected.append(field)

    if not fields_injected:
        print(f"  WARNING: No MSA fields found in {data_json_path}")
        return False

    if not dry_run:
        # Backup original if not already backed up
        bak = input_json_path.with_suffix('.json.pre_msa_bak')
        if not bak.exists():
            input_json_path.rename(bak)
        with open(input_json_path, 'w') as f:
            json.dump(inp, f, indent=2)

    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--template', type=str, default='',
                        help='Process only this template')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--test-one', action='store_true')
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent

    if args.template:
        output_dirs = [script_dir / f"output_{args.template}"]
    else:
        output_dirs = sorted(script_dir.glob("output_G4FP_*"))

    total_ok = 0
    total_fail = 0
    total_skip = 0

    for odir in output_dirs:
        if not odir.is_dir():
            continue
        tpl = odir.name.replace('output_', '')

        for state in ['bound', 'apo']:
            input_dir = odir / f"02_alphafold3_inputs_{state}"
            pred_dir = odir / f"03_alphafold3_predictions_{state}"

            if not input_dir.exists() or not pred_dir.exists():
                continue

            n_ok = 0
            n_fail = 0
            n_skip = 0

            for jf in sorted(input_dir.glob("*.json")):
                if jf.name.endswith('.bak') or jf.name.endswith('.pre_msa_bak'):
                    continue

                # Check if already has MSA
                with open(jf) as f:
                    inp = json.load(f)
                has_msa = False
                for seq in inp['sequences']:
                    if 'protein' in seq:
                        if seq['protein'].get('unpairedMsa'):
                            has_msa = True
                        break
                if has_msa:
                    n_skip += 1
                    continue

                name = inp.get('name', jf.stem)

                # Find corresponding _data.json
                data_json = find_data_json(pred_dir, name)
                if data_json is None:
                    print(f"  MISSING _data.json for {name}")
                    n_fail += 1
                    continue

                if inject_msa(jf, data_json, args.dry_run):
                    n_ok += 1
                else:
                    n_fail += 1

                if args.test_one and n_ok >= 1:
                    break

            print(f"  {tpl}/{state}: {n_ok} injected, {n_skip} already had MSA, {n_fail} failed")
            total_ok += n_ok
            total_fail += n_fail
            total_skip += n_skip

            if args.test_one and total_ok >= 1:
                break

        if args.test_one and total_ok >= 1:
            break

    print(f"\nTotal: {total_ok} injected, {total_skip} skipped, {total_fail} failed")
    if total_fail > 0:
        print("WARNING: Some designs failed — check output above")
        sys.exit(1)


if __name__ == '__main__':
    main()

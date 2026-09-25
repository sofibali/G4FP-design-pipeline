#!/usr/bin/env python3
"""
Fix template AF3 files: insert missing TYG chromophore into sequences.

Same issue as the designs: LigandMPNN omitted CRO, so template sequences
have VTTLVQCF instead of VTTLTYGVQCF (or VTTLQCF → VTTLTYGVQCF for
r28_l2_cro_mod2 which needs TYGV).

Two outputs:
  1. Fixed input JSONs in template_af3_inputs/{bound,apo}/
     → for full AF3 re-run (MSA + inference)
  2. Fixed _data.json + inject MSA into input JSONs
     → for inference-only re-run (--norun_data_pipeline)

Usage:
    python3 fix_template_chromophore.py              # fix everything
    python3 fix_template_chromophore.py --dry-run    # preview only
"""

import json
import re
import sys
import argparse
import numpy as np
from pathlib import Path
from copy import deepcopy


def insert_chromophore(seq):
    """Insert TYG into the sequence at the chromophore site.

    Returns (fixed_seq, insertion_index, inserted_residues) or None.
    """
    idx = seq.find('VTTLVQCF')
    if idx >= 0:
        return seq[:idx+4] + 'TYG' + seq[idx+4:], idx+4, 'TYG'
    idx = seq.find('VTTLQCF')
    if idx >= 0:
        return seq[:idx+4] + 'TYGV' + seq[idx+4:], idx+4, 'TYGV'
    return None


def fix_msa(msa_text, insert_pos, insert_len):
    """Fix an A3M-format MSA by inserting aligned columns at insert_pos.

    For the query sequence: insert uppercase residues (already done in seq).
    For hit sequences: promote any lowercase insertions at that position
    to uppercase, or insert gaps.
    """
    if not msa_text or not msa_text.strip():
        return msa_text

    lines = msa_text.strip().split('\n')
    fixed_lines = []

    for i, line in enumerate(lines):
        if line.startswith('>'):
            fixed_lines.append(line)
            continue

        if i == 1:
            # Query sequence — already has the insertion from seq fix
            # Just pass through (will be replaced by the fixed seq)
            fixed_lines.append(line)
            continue

        # Hit sequence — need to handle insertions
        # Walk through alignment columns to find insertion point
        aln_col = 0
        char_idx = 0
        insert_char_idx = None

        while char_idx < len(line) and aln_col < insert_pos:
            c = line[char_idx]
            if c.isupper() or c == '-':
                aln_col += 1
            # lowercase = insertion, doesn't count as alignment column
            char_idx += 1

        insert_char_idx = char_idx

        # Collect any lowercase insertions at this position
        insertions = []
        while insert_char_idx < len(line) and line[insert_char_idx].islower():
            insertions.append(line[insert_char_idx])
            insert_char_idx += 1

        # Promote insertions to uppercase (or pad with gaps)
        promoted = [c.upper() for c in insertions[:insert_len]]
        while len(promoted) < insert_len:
            promoted.append('-')
        promoted = promoted[:insert_len]

        # Remove the original lowercase insertions and insert promoted
        before = line[:char_idx]
        # Skip past original lowercase insertions
        skip_idx = char_idx
        while skip_idx < len(line) and line[skip_idx].islower():
            skip_idx += 1
        after = line[skip_idx:]

        fixed_line = before + ''.join(promoted) + after
        fixed_lines.append(fixed_line)

    return '\n'.join(fixed_lines)


def fix_template_indices(templates, insert_pos, insert_len):
    """Shift template queryIndices by insert_len for indices >= insert_pos."""
    if not templates:
        return templates
    fixed = []
    for tpl in templates:
        tpl = deepcopy(tpl)
        if 'queryIndices' in tpl:
            tpl['queryIndices'] = [
                idx + insert_len if idx >= insert_pos else idx
                for idx in tpl['queryIndices']
            ]
        fixed.append(tpl)
    return fixed


def fix_json_file(json_path, dry_run=False):
    """Fix a single input JSON or _data.json file."""
    with open(json_path) as f:
        data = json.load(f)

    protein = None
    for seq_entry in data.get('sequences', []):
        if 'protein' in seq_entry:
            protein = seq_entry['protein']
            break

    if protein is None:
        return None, "no protein chain"

    seq = protein['sequence']
    result = insert_chromophore(seq)
    if result is None:
        if 'VTTLTYGVQCF' in seq or 'VTTLTYGVQCF' in seq:
            return None, "already fixed"
        return None, "no VTTLVQCF or VTTLQCF motif found"

    fixed_seq, insert_pos, inserted = result
    protein['sequence'] = fixed_seq

    # Fix MSA fields if present
    for msa_field in ['unpairedMsa', 'pairedMsa']:
        if msa_field in protein and protein[msa_field]:
            msa = protein[msa_field]
            fixed_msa = fix_msa(msa, insert_pos, len(inserted))
            # Replace query line with fixed sequence
            msa_lines = fixed_msa.split('\n')
            if len(msa_lines) >= 2:
                msa_lines[1] = fixed_seq
            protein[msa_field] = '\n'.join(msa_lines)

    # Fix template queryIndices
    if 'templates' in protein:
        protein['templates'] = fix_template_indices(
            protein['templates'], insert_pos, len(inserted))

    if not dry_run:
        # Backup
        bak = json_path.with_suffix('.json.bak')
        if not bak.exists():
            json_path.rename(bak)
            with open(json_path, 'w') as f:
                json.dump(data, f, indent=2)
        else:
            with open(json_path, 'w') as f:
                json.dump(data, f, indent=2)

    return {
        'old_len': len(seq),
        'new_len': len(fixed_seq),
        'insert_pos': insert_pos,
        'inserted': inserted,
    }, None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent

    print("=" * 70)
    print("Fix Template Chromophore (TYG insertion)")
    print("=" * 70)

    # ── Phase 1: Fix input JSONs ──
    print("\nPhase 1: Fixing template input JSONs...")
    inputs_dir = base_dir / "template_af3_inputs"
    n_fixed_inputs = 0
    n_skip_inputs = 0

    for state_dir in sorted(inputs_dir.iterdir()):
        if not state_dir.is_dir():
            continue
        for jf in sorted(state_dir.glob("*.json")):
            if '.bak' in jf.name:
                continue
            result, err = fix_json_file(jf, args.dry_run)
            if result:
                n_fixed_inputs += 1
                print(f"  FIXED: {jf.name} ({result['old_len']}→{result['new_len']}, "
                      f"+{result['inserted']} at pos {result['insert_pos']})")
            else:
                n_skip_inputs += 1
                if err != "already fixed":
                    print(f"  SKIP: {jf.name} ({err})")

    print(f"  Input JSONs: {n_fixed_inputs} fixed, {n_skip_inputs} skipped")

    # ── Phase 2: Fix _data.json files (MSA + sequences) ──
    print("\nPhase 2: Fixing template _data.json files (MSA repair)...")
    outputs_dir = base_dir / "template_af3_outputs"
    n_fixed_data = 0
    n_skip_data = 0

    for dj in sorted(outputs_dir.rglob("*_data.json")):
        if '.bak' in dj.name:
            continue
        result, err = fix_json_file(dj, args.dry_run)
        if result:
            n_fixed_data += 1
            print(f"  FIXED: {dj.relative_to(outputs_dir)} "
                  f"({result['old_len']}→{result['new_len']}, +{result['inserted']})")
        else:
            n_skip_data += 1

    print(f"  _data.json: {n_fixed_data} fixed, {n_skip_data} skipped")

    # ── Phase 3: Inject MSA into input JSONs for inference-only ──
    print("\nPhase 3: Injecting MSA into input JSONs for inference-only re-run...")
    n_injected = 0

    for state in ['bound', 'apo']:
        input_state_dir = inputs_dir / state
        output_state_dir = outputs_dir / state
        if not input_state_dir.exists() or not output_state_dir.exists():
            continue

        for jf in sorted(input_state_dir.glob("*.json")):
            if '.bak' in jf.name:
                continue

            with open(jf) as f:
                inp = json.load(f)

            protein = None
            for seq_entry in inp['sequences']:
                if 'protein' in seq_entry:
                    protein = seq_entry['protein']
                    break
            if protein is None or protein.get('unpairedMsa'):
                continue

            # Find matching _data.json
            name = inp.get('name', jf.stem)
            data_json = None
            for candidate in output_state_dir.rglob("*_data.json"):
                if '.bak' in candidate.name:
                    continue
                try:
                    with open(candidate) as f:
                        ddata = json.load(f)
                    if ddata['sequences'][0]['protein'].get('unpairedMsa'):
                        # Check sequence matches
                        if ddata['sequences'][0]['protein']['sequence'] == protein['sequence']:
                            data_json = candidate
                            break
                except Exception:
                    continue

            if data_json is None:
                print(f"  WARNING: No MSA source for {jf.name}")
                continue

            with open(data_json) as f:
                ddata = json.load(f)
            data_protein = ddata['sequences'][0]['protein']

            for field in ['unpairedMsa', 'pairedMsa', 'templates']:
                if field in data_protein and data_protein[field]:
                    protein[field] = data_protein[field]

            if not args.dry_run:
                bak = jf.with_suffix('.json.pre_msa_bak')
                if not bak.exists():
                    # Already backed up from phase 1, just overwrite
                    pass
                with open(jf, 'w') as f:
                    json.dump(inp, f, indent=2)

            n_injected += 1

    print(f"  Injected MSA into {n_injected} input JSONs")

    # ── Summary ──
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}")
    print(f"  Input JSONs fixed:  {n_fixed_inputs}")
    print(f"  _data.json fixed:   {n_fixed_data}")
    print(f"  MSA injected:       {n_injected}")

    if args.dry_run:
        print("\nDRY RUN — no files modified. Re-run without --dry-run.")
    else:
        print(f"\nDone! Next steps:")
        print(f"  Full re-run (MSA + inference):")
        print(f"    Use input JSONs from template_af3_inputs/{{bound,apo}}/")
        print(f"    WITHOUT --norun_data_pipeline")
        print(f"")
        print(f"  Inference-only re-run (uses patched MSA):")
        print(f"    Use input JSONs from template_af3_inputs/{{bound,apo}}/")
        print(f"    WITH --norun_data_pipeline")
        print(f"    Output to a separate dir (e.g. template_af3_outputs_inference_only/)")


if __name__ == '__main__':
    main()

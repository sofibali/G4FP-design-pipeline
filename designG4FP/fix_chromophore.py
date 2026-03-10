#!/usr/bin/env python3
"""
Fix missing chromophore TYG in all G4FP design sequences.

LigandMPNN omitted the CRO residue instead of expanding it back to
Thr-Tyr-Gly.  This script inserts TYG (or TYGV for r28_l2_cro_mod2)
so all sequences read ...VTTLTYGVQCF... instead of ...VTTLVQCF...

Three modes:
  --step sequences   Fix CSVs, FASTAs, and AF3 input JSONs
  --step msa         Fix _data.json MSA files (for inference-only rerun)
  --step all         Both (default)

Usage:
    python fix_chromophore.py                    # fix everything
    python fix_chromophore.py --step sequences   # just input files
    python fix_chromophore.py --step msa         # just MSA data
    python fix_chromophore.py --dry-run          # preview without writing
    python fix_chromophore.py --test-one         # fix 1 design, then stop
"""

import argparse
import json
import re
import shutil
import sys
from pathlib import Path
from typing import List, Optional, Tuple


# ── Chromophore insertion logic ──────────────────────────────────────

def insert_chromophore(seq: str) -> Tuple[str, int, str]:
    """Insert TYG into sequence. Returns (fixed_seq, insert_col, insert_str).

    9 templates:  VTTLVQCF  → VTTLTYGVQCF  (+3 at col after VTTL)
    1 template:   VTTLQCF   → VTTLTYGVQCF  (+4 at col after VTTL)
    """
    idx = seq.find('VTTLVQCF')
    if idx >= 0:
        insert_col = idx + 4  # after VTTL
        fixed = seq[:insert_col] + 'TYG' + seq[insert_col:]
        return fixed, insert_col, 'TYG'

    idx = seq.find('VTTLQCF')
    if idx >= 0:
        insert_col = idx + 4
        fixed = seq[:insert_col] + 'TYGV' + seq[insert_col:]
        return fixed, insert_col, 'TYGV'

    raise ValueError(f"Cannot find VTTLVQCF or VTTLQCF in sequence "
                     f"(len={len(seq)}, context: ...{seq[180:210]}...)")


def verify_fixed(seq: str) -> bool:
    """Check that a sequence has the chromophore."""
    return 'VTTLTYGVQCF' in seq


# ── A3M MSA repair ──────────────────────────────────────────────────

def fix_a3m_msa(msa_text: str, insert_col: int, n_insert: int) -> str:
    """Fix an A3M-format MSA by inserting aligned columns at insert_col.

    For each hit sequence:
      - Walk uppercase/gap chars to find position insert_col
      - Collect any lowercase insertions at that position
      - If n lowercase == n_insert: promote to uppercase (aligned)
      - If n lowercase < n_insert: promote what exists, pad with gaps
      - If n lowercase > n_insert: promote first n_insert, keep rest lowercase
      - If n lowercase == 0: insert n_insert gap chars
    """
    if not msa_text.strip():
        return msa_text

    lines = msa_text.rstrip('\n').split('\n')
    fixed_lines = []

    for i, line in enumerate(lines):
        if line.startswith('>') or line.startswith('#'):
            fixed_lines.append(line)
            continue

        # Process sequence line
        result = []
        col = 0  # counts uppercase + gap positions (alignment columns)
        j = 0
        inserted = False

        while j < len(line):
            ch = line[j]

            if ch == '-' or ch.isupper():
                if col == insert_col and not inserted:
                    inserted = True
                    # Collect any lowercase insertions right before this position
                    # (they were already added to result as we encountered them)
                    # Now we need to handle lowercase chars AFTER this col
                    # Actually, in A3M, lowercase chars come AFTER the preceding
                    # uppercase/gap char and BEFORE the next one.
                    # So at col == insert_col, any lowercase chars between
                    # the previous uppercase and this one were already processed.

                    # Check if there are lowercase chars following from a
                    # previous iteration — they're already in result.
                    # We need to look at what lowercase chars were collected
                    # between col insert_col-1 and col insert_col.

                    # Let me re-approach: scan backwards in result to find
                    # lowercase chars that belong to position insert_col
                    lc_chars = []
                    while result and result[-1].islower():
                        lc_chars.append(result.pop())
                    lc_chars.reverse()
                    lc_str = ''.join(lc_chars)

                    if len(lc_str) == n_insert:
                        # Perfect match: promote to uppercase
                        result.append(lc_str.upper())
                    elif len(lc_str) < n_insert:
                        # Fewer insertions: promote + pad with gaps
                        result.append(lc_str.upper() + '-' * (n_insert - len(lc_str)))
                    else:
                        # More insertions: promote first n_insert, keep rest lowercase
                        result.append(lc_str[:n_insert].upper())
                        result.append(lc_str[n_insert:])  # stays lowercase

                    # Now add the current character
                    result.append(ch)
                    col += 1
                else:
                    result.append(ch)
                    col += 1
            else:
                # Lowercase insertion character
                result.append(ch)

            j += 1

        # If insert_col is at/past the end (sequence shorter than query)
        if not inserted:
            # Sequence is shorter; add gaps at the end
            result.append('-' * n_insert)

        fixed_lines.append(''.join(result))

    return '\n'.join(fixed_lines) + '\n'


def fix_query_sequence_in_msa(msa_text: str, insert_col: int,
                               insert_str: str) -> str:
    """Fix the query (first) sequence in an A3M MSA.

    The query has no lowercase chars, so just insert uppercase at the position.
    """
    if not msa_text.strip():
        return msa_text

    lines = msa_text.rstrip('\n').split('\n')
    # Find the first non-header line (query sequence)
    for i, line in enumerate(lines):
        if not line.startswith('>') and not line.startswith('#'):
            # This is the query sequence - insert at column position
            lines[i] = line[:insert_col] + insert_str + line[insert_col:]
            break

    return '\n'.join(lines) + '\n'


def fix_msa_complete(msa_text: str, insert_col: int,
                     insert_str: str) -> str:
    """Fix both query and hit sequences in an A3M MSA."""
    n_insert = len(insert_str)

    # Step 1: Fix query sequence (simple insertion)
    fixed = fix_query_sequence_in_msa(msa_text, insert_col, insert_str)

    # Step 2: Fix all hit sequences (handle lowercase insertions)
    # We need to re-parse after query fix since line indices shifted
    lines = fixed.rstrip('\n').split('\n')
    result_lines = []
    query_done = False

    for line in lines:
        if line.startswith('>') or line.startswith('#'):
            result_lines.append(line)
            continue

        if not query_done:
            # Query sequence already fixed
            result_lines.append(line)
            query_done = True
            continue

        # Hit sequence: handle lowercase insertions at insert_col
        result = []
        col = 0
        inserted = False

        for ch in line:
            if ch == '-' or ch.isupper():
                if col == insert_col and not inserted:
                    inserted = True
                    # Collect preceding lowercase chars from result
                    lc_chars = []
                    while result and result[-1].islower():
                        lc_chars.append(result.pop())
                    lc_chars.reverse()
                    lc_str = ''.join(lc_chars)

                    if len(lc_str) >= n_insert:
                        result.append(lc_str[:n_insert].upper())
                        if len(lc_str) > n_insert:
                            result.append(lc_str[n_insert:])
                    else:
                        result.append(lc_str.upper() + '-' * (n_insert - len(lc_str)))

                    result.append(ch)
                    col += 1
                else:
                    result.append(ch)
                    col += 1
            else:
                result.append(ch)

        if not inserted:
            result.append('-' * n_insert)

        result_lines.append(''.join(result))

    return '\n'.join(result_lines) + '\n'


def fix_template_indices(templates: List[dict], insert_col: int,
                         n_insert: int) -> List[dict]:
    """Shift queryIndices >= insert_col by n_insert."""
    fixed = []
    for t in templates:
        t_new = dict(t)
        qi = t.get('queryIndices', [])
        t_new['queryIndices'] = [q + n_insert if q >= insert_col else q
                                  for q in qi]
        fixed.append(t_new)
    return fixed


# ── Phase 1: Fix sequences in CSVs, FASTAs, AF3 input JSONs ─────────

def fix_sequences(base_dir: Path, dry_run: bool = False,
                  test_one: bool = False) -> int:
    """Fix all sequence files. Returns count of files modified."""
    count = 0

    for output_dir in sorted(base_dir.glob("output_G4FP_*")):
        if not output_dir.is_dir():
            continue
        tpl = output_dir.name.replace("output_", "")
        print(f"\n{'='*60}")
        print(f"Template: {tpl}")
        print(f"{'='*60}")

        # 1a. Fix filtered sequences CSV
        csv_path = output_dir / "02_filtered_sequences.csv"
        if csv_path.exists():
            import csv
            with open(csv_path) as f:
                reader = csv.DictReader(f)
                rows = list(reader)
                fieldnames = reader.fieldnames

            n_fixed = 0
            for row in rows:
                if 'sequence' in row and not verify_fixed(row['sequence']):
                    try:
                        row['sequence'], _, _ = insert_chromophore(row['sequence'])
                        n_fixed += 1
                    except ValueError as e:
                        print(f"  WARNING: {e}")

            if n_fixed > 0:
                print(f"  CSV: {n_fixed}/{len(rows)} sequences fixed")
                if not dry_run:
                    shutil.copy2(csv_path, csv_path.with_suffix('.csv.bak'))
                    with open(csv_path, 'w', newline='') as f:
                        writer = csv.DictWriter(f, fieldnames=fieldnames)
                        writer.writeheader()
                        writer.writerows(rows)
                count += 1
            else:
                print(f"  CSV: already correct or not found")

        # 1b. Fix filtered sequences FASTA
        fa_path = output_dir / "02_filtered_sequences.fa"
        if fa_path.exists():
            lines = fa_path.read_text().splitlines()
            n_fixed = 0
            for i, line in enumerate(lines):
                if not line.startswith('>') and 'VTTL' in line:
                    if not verify_fixed(line):
                        try:
                            lines[i], _, _ = insert_chromophore(line)
                            n_fixed += 1
                        except ValueError:
                            pass
            if n_fixed > 0:
                print(f"  FASTA: {n_fixed} sequences fixed")
                if not dry_run:
                    shutil.copy2(fa_path, fa_path.with_suffix('.fa.bak'))
                    fa_path.write_text('\n'.join(lines) + '\n')
                count += 1

        # 1c. Fix AF3 input JSONs
        for state in ['bound', 'apo']:
            json_dir = output_dir / f"02_alphafold3_inputs_{state}"
            if not json_dir.exists():
                continue
            n_fixed = 0
            for jf in sorted(json_dir.glob("*.json")):
                if jf.name.endswith('.bak'):
                    continue
                data = json.loads(jf.read_text())
                seq = data['sequences'][0]['protein']['sequence']
                if verify_fixed(seq):
                    continue
                try:
                    fixed_seq, _, _ = insert_chromophore(seq)
                    data['sequences'][0]['protein']['sequence'] = fixed_seq
                    n_fixed += 1
                    if not dry_run:
                        shutil.copy2(jf, jf.with_suffix('.json.bak'))
                        jf.write_text(json.dumps(data, indent=2) + '\n')
                    if test_one and n_fixed >= 1:
                        break
                except ValueError as e:
                    print(f"  WARNING: {jf.name}: {e}")

            print(f"  AF3 {state} JSONs: {n_fixed} fixed")
            count += n_fixed

        if test_one:
            print("\n  --test-one: stopping after first template")
            break

    return count


# ── Phase 2: Fix _data.json MSA files ────────────────────────────────

def fix_msa_files(base_dir: Path, dry_run: bool = False,
                  test_one: bool = False) -> int:
    """Fix all _data.json MSA files. Returns count of files modified."""
    count = 0

    for output_dir in sorted(base_dir.glob("output_G4FP_*")):
        if not output_dir.is_dir():
            continue
        tpl = output_dir.name.replace("output_", "")
        print(f"\n{'='*60}")
        print(f"Template: {tpl}")
        print(f"{'='*60}")

        for state in ['bound', 'apo']:
            pred_dir = output_dir / f"03_alphafold3_predictions_{state}"
            if not pred_dir.exists():
                continue

            n_fixed = 0
            n_skip = 0
            n_err = 0

            # Find all _data.json files (nested: design_NNNN_state/design_N_state/)
            data_jsons = sorted(pred_dir.glob("*/*/design_*_data.json"))
            if not data_jsons:
                data_jsons = sorted(pred_dir.glob("*/design_*_data.json"))

            for djf in data_jsons:
                try:
                    data = json.loads(djf.read_text())
                    prot = data['sequences'][0]['protein']
                    seq = prot['sequence']

                    if verify_fixed(seq):
                        n_skip += 1
                        continue

                    fixed_seq, insert_col, insert_str = insert_chromophore(seq)
                    n_insert = len(insert_str)

                    # Fix protein sequence
                    prot['sequence'] = fixed_seq

                    # Fix unpaired MSA
                    if prot.get('unpairedMsa', '').strip():
                        prot['unpairedMsa'] = fix_msa_complete(
                            prot['unpairedMsa'], insert_col, insert_str)

                    # Fix paired MSA
                    if prot.get('pairedMsa', '').strip():
                        prot['pairedMsa'] = fix_msa_complete(
                            prot['pairedMsa'], insert_col, insert_str)

                    # Fix template queryIndices
                    if prot.get('templates'):
                        prot['templates'] = fix_template_indices(
                            prot['templates'], insert_col, n_insert)

                    n_fixed += 1

                    if not dry_run:
                        bak = djf.with_suffix('.json.bak')
                        if not bak.exists():
                            shutil.copy2(djf, bak)
                        djf.write_text(json.dumps(data, indent=2) + '\n')

                    if test_one and n_fixed >= 1:
                        # Validate the fix
                        print(f"\n  TEST: {djf.name}")
                        print(f"    Old seq len: {len(seq)}")
                        print(f"    New seq len: {len(fixed_seq)}")
                        print(f"    Insert: '{insert_str}' at col {insert_col}")
                        print(f"    Verify: {'VTTLTYGVQCF' in fixed_seq}")

                        # Check MSA query
                        msa_lines = prot['unpairedMsa'].strip().split('\n')
                        for ml in msa_lines:
                            if not ml.startswith('>'):
                                print(f"    MSA query len: {len(ml)}")
                                idx = ml.find('VTTLTYGVQCF')
                                print(f"    MSA query VTTLTYGVQCF at: {idx}")
                                break

                        # Check a hit
                        hit_count = 0
                        for ml in msa_lines[2:6]:
                            if not ml.startswith('>'):
                                hit_count += 1
                                print(f"    Hit {hit_count} len: {len(ml)} "
                                      f"(should be {len(fixed_seq)}+ with lowercase)")
                        break

                except Exception as e:
                    n_err += 1
                    print(f"  ERROR: {djf.name}: {e}")

            print(f"  {state}: {n_fixed} fixed, {n_skip} already OK, "
                  f"{n_err} errors (of {len(data_jsons)} files)")
            count += n_fixed

            if test_one and n_fixed > 0:
                break

        if test_one:
            print("\n  --test-one: stopping after first template")
            break

    return count


# ── Main ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Fix missing chromophore TYG in G4FP design sequences')
    parser.add_argument('--step', choices=['sequences', 'msa', 'all'],
                        default='all',
                        help='Which step to run (default: all)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Preview changes without writing files')
    parser.add_argument('--test-one', action='store_true',
                        help='Fix only 1 design per template for testing')
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent
    print(f"Base directory: {base_dir}")
    if args.dry_run:
        print("DRY RUN — no files will be modified")

    total = 0

    if args.step in ('sequences', 'all'):
        print(f"\n{'#'*60}")
        print(f"Phase 1: Fixing sequence files")
        print(f"{'#'*60}")
        n = fix_sequences(base_dir, dry_run=args.dry_run,
                          test_one=args.test_one)
        total += n
        print(f"\nPhase 1 complete: {n} files modified")

    if args.step in ('msa', 'all'):
        print(f"\n{'#'*60}")
        print(f"Phase 2: Fixing _data.json MSA files")
        print(f"{'#'*60}")
        n = fix_msa_files(base_dir, dry_run=args.dry_run,
                          test_one=args.test_one)
        total += n
        print(f"\nPhase 2 complete: {n} files modified")

    print(f"\n{'='*60}")
    print(f"Total files modified: {total}")
    if args.dry_run:
        print("(DRY RUN — no files were actually modified)")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Validate AlphaFold3 input JSON files before submitting expensive predictions.

Checks:
  - Required fields: dialect, version, name, modelSeeds, sequences
  - Each sequence entry has an 'id' field (required by AF3 v2+)
  - Protein sequences are valid amino acids
  - DNA sequences are valid nucleotides
  - Ion/ligand entries have type and count
  - No duplicate IDs within a file
  - Bound files contain protein + DNA; apo files contain protein only

Usage:
    python utils/validate_af3_jsons.py                     # validate all output dirs
    python utils/validate_af3_jsons.py --output-dir DIR    # validate specific dir
    python utils/validate_af3_jsons.py --fix               # auto-fix missing id fields
"""

import json
import sys
import argparse
from pathlib import Path
from typing import List, Tuple

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
VALID_DNA = set("ACGT")

def validate_json(filepath: Path) -> List[str]:
    """Validate a single AF3 input JSON. Returns list of error strings."""
    errors = []

    try:
        with open(filepath) as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        return [f"Invalid JSON: {e}"]

    # Top-level fields
    for field in ("dialect", "version", "name", "modelSeeds", "sequences"):
        if field not in data:
            errors.append(f"Missing top-level field: {field}")

    if data.get("dialect") != "alphafold3":
        errors.append(f"dialect should be 'alphafold3', got '{data.get('dialect')}'")

    if not isinstance(data.get("sequences"), list) or len(data.get("sequences", [])) == 0:
        errors.append("sequences must be a non-empty list")
        return errors

    seen_ids = set()
    has_protein = False
    has_dna = False

    for i, seq in enumerate(data["sequences"]):
        # Determine the sequence type dict (protein, dna, ligand, rna)
        type_key = None
        type_dict = None
        for k in ("protein", "dna", "rna", "ligand"):
            if k in seq:
                type_key = k
                type_dict = seq[k]
                break

        # Check id field -- AF3 3.0.1 expects id INSIDE the type dict, not at top level
        if "id" in seq and type_key:
            errors.append(f"sequences[{i}]: 'id' is at top level, must be inside '{type_key}' dict")
        if type_dict is not None and isinstance(type_dict, dict):
            sid = type_dict.get("id")
            if sid is None:
                errors.append(f"sequences[{i}].{type_key}: missing 'id' field")
            else:
                if sid in seen_ids:
                    errors.append(f"sequences[{i}].{type_key}: duplicate id '{sid}'")
                seen_ids.add(sid)
        elif type_key is None:
            if "id" in seq:
                errors.append(f"sequences[{i}]: has 'id' but no sequence type (protein/dna/ligand)")
            else:
                errors.append(f"sequences[{i}]: unknown entry type (keys: {list(seq.keys())})")

        # Check sequence type content
        if "protein" in seq:
            has_protein = True
            aa_seq = seq["protein"].get("sequence", "")
            bad_chars = set(aa_seq.upper()) - VALID_AA
            if bad_chars:
                errors.append(f"sequences[{i}] protein: invalid chars {bad_chars}")
            if len(aa_seq) == 0:
                errors.append(f"sequences[{i}] protein: empty sequence")

        elif "dna" in seq:
            has_dna = True
            dna_seq = seq["dna"].get("sequence", "")
            bad_chars = set(dna_seq.upper()) - VALID_DNA
            if bad_chars:
                errors.append(f"sequences[{i}] DNA: invalid chars {bad_chars}")

        elif "ligand" in seq:
            lig = seq["ligand"]
            if "type" not in lig and "ccd" not in lig and "smiles" not in lig:
                errors.append(f"sequences[{i}] ligand: needs type, ccd, or smiles")

    if not has_protein:
        errors.append("No protein sequence found")

    # Bound/apo consistency check based on filename
    name = filepath.stem
    if "_bound" in name and not has_dna:
        errors.append("Bound file missing DNA sequence")
    if "_apo" in name and has_dna:
        errors.append("Apo file should not contain DNA sequence")

    return errors


def fix_json(filepath: Path) -> bool:
    """Fix common AF3 JSON issues: dialect, version, and id placement.

    AF3 3.0.1 requires id INSIDE each sequence type dict (protein/dna/ligand),
    not as a sibling key at the top level of the sequence entry.
    """
    with open(filepath) as f:
        data = json.load(f)

    changed = False

    # Fix missing dialect
    if "dialect" not in data:
        data["dialect"] = "alphafold3"
        changed = True

    # Fix missing version
    if "version" not in data:
        data["version"] = 1
        changed = True

    # Fix id placement: move top-level id into the type dict
    for i, seq in enumerate(data.get("sequences", [])):
        # Find the type key (protein, dna, rna, ligand)
        type_key = None
        for k in ("protein", "dna", "rna", "ligand"):
            if k in seq:
                type_key = k
                break

        if type_key is None:
            continue

        type_dict = seq[type_key]
        if not isinstance(type_dict, dict):
            continue

        # Case 1: id at top level but not inside type dict -> move it
        if "id" in seq and "id" not in type_dict:
            type_dict["id"] = seq.pop("id")
            # Reorder type dict to put id first
            reordered = {"id": type_dict["id"]}
            for k, v in type_dict.items():
                if k != "id":
                    reordered[k] = v
            seq[type_key] = reordered
            changed = True

        # Case 2: id at top level AND inside type dict -> remove top-level
        elif "id" in seq and "id" in type_dict:
            del seq["id"]
            changed = True

        # Case 3: no id anywhere -> assign one
        elif "id" not in seq and "id" not in type_dict:
            type_dict["id"] = str(i + 1)
            reordered = {"id": type_dict["id"]}
            for k, v in type_dict.items():
                if k != "id":
                    reordered[k] = v
            seq[type_key] = reordered
            changed = True

    if changed:
        # Reorder top-level keys
        ordered = {}
        for key in ("dialect", "version", "name", "modelSeeds", "sequences"):
            if key in data:
                ordered[key] = data[key]
        for key in data:
            if key not in ordered:
                ordered[key] = data[key]

        with open(filepath, "w") as f:
            json.dump(ordered, f, indent=2)

    return changed


def main():
    parser = argparse.ArgumentParser(description="Validate AlphaFold3 input JSONs")
    parser.add_argument("--output-dir", type=str,
                        help="Specific output directory (default: all output_G4FP_* dirs)")
    parser.add_argument("--fix", action="store_true",
                        help="Auto-fix missing id fields")
    args = parser.parse_args()

    base = Path(__file__).resolve().parent.parent  # designG4FP/

    if args.output_dir:
        output_dirs = [Path(args.output_dir)]
    else:
        output_dirs = sorted(d for d in base.iterdir()
                             if d.is_dir() and d.name.startswith("output_G4FP_"))

    if not output_dirs:
        print("No output directories found.")
        sys.exit(1)

    total_files = 0
    total_errors = 0
    total_fixed = 0

    for odir in output_dirs:
        for subdir_name in ("02_alphafold3_inputs_bound", "02_alphafold3_inputs_apo"):
            json_dir = odir / subdir_name
            if not json_dir.exists():
                continue

            json_files = sorted(json_dir.glob("*.json"))
            dir_errors = 0

            for jf in json_files:
                total_files += 1
                errs = validate_json(jf)

                if errs and args.fix:
                    if fix_json(jf):
                        total_fixed += 1
                        errs = validate_json(jf)  # re-validate

                if errs:
                    dir_errors += len(errs)
                    total_errors += len(errs)
                    for err in errs:
                        print(f"  ERROR {jf.name}: {err}")

            status = "OK" if dir_errors == 0 else f"{dir_errors} errors"
            print(f"{odir.name}/{subdir_name}: {len(json_files)} files, {status}")

    print(f"\nSummary: {total_files} files checked, {total_errors} errors"
          + (f", {total_fixed} fixed" if total_fixed else ""))

    sys.exit(1 if total_errors > 0 else 0)


if __name__ == "__main__":
    main()

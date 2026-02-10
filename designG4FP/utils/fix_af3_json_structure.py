#!/usr/bin/env python3
"""
Fix AlphaFold3 JSON inputs by moving 'id' fields to correct location.

The id field must be INSIDE the protein/dna/ligand object, not at the sequence level.

Correct structure:
{
  "sequences": [
    {
      "protein": {
        "id": "A",
        "sequence": "..."
      }
    }
  ]
}

NOT:
{
  "sequences": [
    {
      "protein": {"sequence": "..."},
      "id": "1"
    }
  ]
}
"""

import json
from pathlib import Path
import sys

def fix_json_structure(json_path: Path, dry_run: bool = False) -> bool:
    """Fix JSON by moving id fields to correct location"""
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        needs_fix = False
        
        if 'sequences' in data:
            for seq in data['sequences']:
                # Check if id is at wrong level (sequence level instead of inside protein/dna/etc)
                if 'id' in seq:
                    needs_fix = True
                    break
                
                # Also check if id is missing from protein/dna/ligand level
                for key in ['protein', 'dna', 'ligand', 'rna']:
                    if key in seq and isinstance(seq[key], dict):
                        if 'id' not in seq[key]:
                            needs_fix = True
                            break
        
        if not needs_fix:
            return False
        
        # Fix the structure
        for i, seq in enumerate(data['sequences']):
            # If there's an id at sequence level, move it
            if 'id' in seq:
                seq_id = seq.pop('id')
            else:
                seq_id = str(i + 1)
            
            # Move id to the correct location
            if 'protein' in seq:
                if 'id' not in seq['protein']:
                    seq['protein']['id'] = seq_id if i == 0 else chr(65 + i)  # A, B, C...
            elif 'dna' in seq:
                if 'id' not in seq['dna']:
                    seq['dna']['id'] = seq_id if seq_id.isalpha() else chr(65 + i)
            elif 'rna' in seq:
                if 'id' not in seq['rna']:
                    seq['rna']['id'] = seq_id if seq_id.isalpha() else chr(65 + i)
            elif 'ligand' in seq:
                if 'id' not in seq['ligand']:
                    # Ligands can use numbers
                    seq['ligand']['id'] = seq_id
            elif 'type' in seq:  # Ion
                # Ions don't get an id in the container, just the type
                pass
        
        if dry_run:
            return True
        
        # Write back
        with open(json_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        return True
        
    except Exception as e:
        print(f"  ✗ Error: {str(e)}")
        return False


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--output-dir', type=str)
    args = parser.parse_args()
    
    # Find output directories
    if args.output_dir:
        output_dirs = [Path(args.output_dir)]
    else:
        output_dirs = [d for d in Path('.').iterdir() 
                      if d.is_dir() and d.name.startswith('output_')]
    
    print(f"{'DRY RUN - ' if args.dry_run else ''}Fixing AlphaFold3 JSON structure...")
    print(f"Moving 'id' fields to correct location\n")
    
    total_fixed = 0
    total_checked = 0
    
    for output_dir in sorted(output_dirs):
        print(f"Processing: {output_dir.name}")
        
        for state in ['bound', 'apo']:
            input_dir = output_dir / f"02_alphafold3_inputs_{state}"
            
            if not input_dir.exists():
                continue
            
            json_files = list(input_dir.glob("*.json"))
            if not json_files:
                continue
            
            print(f"  {state}: {len(json_files)} files")
            
            fixed_count = 0
            for json_file in json_files:
                if fix_json_structure(json_file, args.dry_run):
                    fixed_count += 1
            
            if fixed_count > 0:
                print(f"    {'Would fix' if args.dry_run else 'Fixed'}: {fixed_count}")
            else:
                print(f"    ✓ Already correct")
            
            total_fixed += fixed_count
            total_checked += len(json_files)
        
        print()
    
    print("=" * 60)
    print(f"Total: {total_checked} checked, {total_fixed} {'need' if args.dry_run else ''} fixed")
    
    if args.dry_run and total_fixed > 0:
        print("\nRun without --dry-run to apply fixes")
    elif total_fixed > 0:
        print("\n✓ All JSON files fixed! Ready to run AlphaFold3.")
    else:
        print("\n✓ All files already have correct structure.")


if __name__ == "__main__":
    main()

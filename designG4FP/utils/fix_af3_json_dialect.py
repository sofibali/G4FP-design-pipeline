#!/usr/bin/env python3
"""
Fix AlphaFold3 JSON input files by adding required 'dialect' and 'version' fields.

AlphaFold3 requires these fields:
  "dialect": "alphafold3",
  "version": 1

Usage: python fix_af3_json_dialect.py <input_dir>
"""

import json
import sys
from pathlib import Path


def fix_json_file(json_path: Path, backup: bool = True) -> bool:
    """Fix a single JSON file by adding dialect and version fields."""
    try:
        # Read existing JSON
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        # Check if already has required fields
        if 'dialect' in data and 'version' in data:
            print(f"  ✓ {json_path.name} - already has required fields")
            return True
        
        # Backup original file
        if backup:
            backup_path = json_path.with_suffix('.json.bak')
            if not backup_path.exists():
                with open(backup_path, 'w') as f:
                    json.dump(data, f, indent=2)
        
        # Add required fields at the beginning
        fixed_data = {
            "dialect": "alphafold3",
            "version": 1
        }
        fixed_data.update(data)
        
        # Write fixed JSON
        with open(json_path, 'w') as f:
            json.dump(fixed_data, f, indent=2)
        
        print(f"  ✓ {json_path.name} - fixed")
        return True
    
    except Exception as e:
        print(f"  ❌ {json_path.name} - error: {e}")
        return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python fix_af3_json_dialect.py <input_dir>")
        print("\nExample:")
        print("  python fix_af3_json_dialect.py output_G4FP_des1_cro_mod0/02_alphafold3_inputs_bound/")
        sys.exit(1)
    
    input_dir = Path(sys.argv[1])
    
    if not input_dir.exists():
        print(f"Error: Directory not found: {input_dir}")
        sys.exit(1)
    
    # Find all JSON files
    json_files = sorted(input_dir.glob("*.json"))
    
    if not json_files:
        print(f"No JSON files found in {input_dir}")
        sys.exit(1)
    
    print(f"Found {len(json_files)} JSON files in {input_dir}")
    print(f"Adding 'dialect' and 'version' fields...")
    print()
    
    fixed_count = 0
    for json_path in json_files:
        if fix_json_file(json_path):
            fixed_count += 1
    
    print()
    print(f"✓ Successfully fixed {fixed_count}/{len(json_files)} files")
    print()
    print("Backups saved with .json.bak extension")
    print("You can now re-run AlphaFold3 predictions")


if __name__ == "__main__":
    main()

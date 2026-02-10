#!/usr/bin/env python3
"""
Fix AlphaFold3 JSON input files - wrap ligand type in proper structure.

The issue: ligand entries like {"type": "K", "count": 3, "id": "3"}
Should be: {"ligand": {"type": "K", "count": 3}, "id": "3"}
"""

import json
import sys
from pathlib import Path

def fix_json_file(filepath):
    """Fix a single JSON file."""
    with open(filepath, 'r') as f:
        data = json.load(f)
    
    fixed = False
    if 'sequences' in data:
        for i, seq_entry in enumerate(data['sequences']):
            # Check if this is a malformed ligand entry
            if 'type' in seq_entry and 'ligand' not in seq_entry:
                # Extract the ligand info
                ligand_type = seq_entry.pop('type')
                count = seq_entry.pop('count', 1)
                
                # Restructure properly
                seq_entry['ligand'] = {
                    'type': ligand_type,
                    'count': count
                }
                fixed = True
    
    if fixed:
        # Write back with proper formatting
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)
        return True
    return False

def main():
    if len(sys.argv) < 2:
        print("Usage: python fix_af3_json.py <directory_or_file>")
        print("  Fixes AlphaFold3 JSON ligand format issues")
        sys.exit(1)
    
    path = Path(sys.argv[1])
    
    if path.is_file():
        files = [path]
    elif path.is_dir():
        files = list(path.glob("*.json"))
    else:
        print(f"Error: {path} is not a valid file or directory")
        sys.exit(1)
    
    fixed_count = 0
    for filepath in files:
        if fix_json_file(filepath):
            fixed_count += 1
            print(f"✓ Fixed: {filepath.name}")
    
    print(f"\nFixed {fixed_count}/{len(files)} files")

if __name__ == '__main__':
    main()
